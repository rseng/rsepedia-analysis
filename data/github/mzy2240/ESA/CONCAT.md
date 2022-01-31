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
Release Process
===============

This document describes the steps needed to publish a new release of
ESA.

#.  (Optional) If any updates has been made to the ``_performance.py``
    then you need to compile the module first by running
    ``compile_for_all_versions.py``.
#.  Ensure you have checked out the ``develop`` branch and have a clean
    repository (no local changes, new files, etc.).
#.  Run all tests for all Python versions (3.5 - 3.8) by running the
    script ``tests\run_tests_for_all_python_versions.py`` or manually
    run all tests in all Python environments by running the following
    from the top level of the repository after activating each virtual
    environment:
    ``python -m unittest discover tests``.
#.  Assess testing coverage by running ``docs\rst\coverage_to_rst.py``.
#.  Check the top-level README file - if testing coverage is *NOT* at
    100%, we need to add tests to get it there before proceeding. Add
    the tests, and start over.
#.  Update ``VERSION`` file (careful not to add a newline at the end).
#.  Update ``docs\rst\welcome\changelog.rst``. Add detailed notes
    related to what has changed in this release. This should include
    new functions/features added, bugs fixed, and any changes that
    break backwards compatibility (ideally, there should not be any
    breaking changes).
#.  Install esa locally into your virtual environment that you plan to
    use for building the documentation (see next step) by executing
    ``python -m pip install -e .[test,doc]``. Note the ``-e`` specifies
    "editable" mode so that the installation process does not copy all
    the files into your virtual environment's ``Lib\site-packages``
    directory. This local, editable installation will automatically
    update ``__version__`` in ``__init__.py``.
#.  If all tests were successful, build the documentation (see README
    in ``docs\rst`` directory). Note that there should **NOT** be
    **ANY** warnings or errors in the console output when building the
    documentation.
#.  When ready, commit all the new changes to the ``develop`` branch.
#.  Checkout the ``master`` branch, and run ``git merge develop``
#.  After merging, add a version tag like so:
    ``git tag -a v1.0.1 -m "ESA version 1.0.1"``
#.  Run ``git push``
#.  Run ``git push --tags``
#.  Ensure setuptools and wheel are up to date:
    ``python -m pip install --upgrade setuptools wheel``
#.  Generate the distribution archive by running this command from the
    same directory where ``setup.py`` is located:
    ``python setup.py sdist bdist_wheel``
#.  Upload the distribution archive to the Python Package Index by
    running this command: ``python -m twine upload dist/*``.
    Before uploading make sure there is only one archive version in the
    ``dist`` directory.Easy SimAuto (ESA)
==================
.. image:: https://img.shields.io/pypi/v/esa.svg
   :target: https://pypi.org/project/esa/
.. image:: https://joss.theoj.org/papers/10.21105/joss.02289/status.svg
   :target: https://doi.org/10.21105/joss.02289
.. image:: https://img.shields.io/pypi/l/esa.svg
   :target: https://github.com/mzy2240/ESA/blob/master/LICENSE
.. image:: https://static.pepy.tech/personalized-badge/esa?period=total&units=international_system&left_color=grey&right_color=orange&left_text=Downloads
   :target: https://pepy.tech/project/esa
.. image:: https://img.shields.io/badge/coverage-100%25-brightgreen


Easy SimAuto (ESA) is an easy-to-use Power System Analysis Automation
Platform atop PowerWorld's Simulator Automation Server (SimAuto).
ESA wraps all PowerWorld SimAuto functions, supports Auxiliary scripts,
provides helper functions to further simplify working with SimAuto and
also turbocharges with native implementation of SOTA algorithms. Wherever
possible, data is returned as Pandas DataFrames, making analysis a breeze.
ESA is well tested and fully `documented`_.

`Documentation`_
----------------

For quick-start directions, installation instructions, API reference,
examples, and more, please check out ESA's `documentation`_.

If you have your own copy of the ESA repository, you can also view the
documentation locally by navigating to the directory ``docs/html`` and
opening ``index.html`` with your web browser.

Citation
--------

If you use ESA in any of your work, please use the citation below.

.. code:: latex

    @article{ESA,
      doi = {10.21105/joss.02289},
      url = {https://doi.org/10.21105/joss.02289},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {50},
      pages = {2289},
      author = {Brandon L. Thayer and Zeyu Mao and Yijing Liu and Katherine Davis and Thomas J. Overbye},
      title = {Easy SimAuto (ESA): A Python Package that Simplifies Interacting with PowerWorld Simulator},
      journal = {Journal of Open Source Software}
    }

Installation
------------

Please refer to ESA's `documentation <https://mzy2240.github
.io/ESA/html/installation.html>`__ for full, detailed installation
directions. In many cases, ESA can simply be installed by:

.. code:: bat

    python -m pip install esa

Testing Coverage
----------------

The ESA team works hard to ensure ESA is well tested, and we strive for
100% testing coverage. The table below shows the most up-to-date
testing coverage data for ESA, using `coverage
<https://pypi.org/project/coverage/>`__.

.. table:: ESA's testing coverage as of 2022-01-28 (Git commit: b8cef3a)
    :widths: auto
    :align: left

    +-----------------+-------------------+-----------------+-----------------+--------------------+
    | Name            |   Num. Statements |   Missing Lines |   Covered Lines |   Percent Coverage |
    +=================+===================+=================+=================+====================+
    | esa/__init__.py |                 2 |               0 |               2 |                100 |
    +-----------------+-------------------+-----------------+-----------------+--------------------+
    | esa/saw.py      |               726 |               0 |             726 |                100 |
    +-----------------+-------------------+-----------------+-----------------+--------------------+

License
-------

`MIT <https://choosealicense.com/licenses/mit/>`__

Contributing
------------

We welcome contributions! Please read ``contributing.md``.

.. _documentation: https://mzy2240.github.io/ESA/
.. _documented: https://mzy2240.github.io/ESA/
Tests
=====

This directory should contain all code, PowerWorld cases, and data
needed to fully test ESA. There should be a Python module named 
"test_<module>" for each ESA module. At present, ESA only has the 
SAW module.

Cases
-----

Use this directory to store PowerWorld cases.

Data
----

Use this directory to store other data needed for testing.

area_filter.aux
--------------

PowerWorld auxiliary file that defines a filter which obtains only buses in
the east area. Used for testing ProcessAuxFile, and is also a useful
template for defining filters in aux files.

README.rst
----------

This file.

run_tests_for_all_python_versions.py
------------------------------------

Script to run tests for all supported Python versions. Check out its
docstring for more information.

test_saw.py
-----------

Python file for running tests related to the SAW module.

test_snippets.py
----------------

Python file for running tests related to the documentation snippets.cases
=====

Use this directory to store PowerWorld cases for testing. Make sure the
cases are public - don't share any private cases! Ideally, you'll
include the case saved in a variety of Simulator formats (e.g. 16, 17,
18, 19, 20, 21) so that users who do not have the same Simulator version
as you can still leverage them. Please note the naming conventions:
cases end with ``_pws_version_<VERSION GOES HERE>.pwb``.

dummy_case.PWB
--------------
This is a case that literally has a single bus and nothing else. This is
intentional - it's only used so that during testing a ``SAW`` instance
can be quickly loaded up to determine the Simulator version. It's saved
in the Simulator 16 format.

ieee_14
-------

IEEE 14-bus test case downloaded from `here <https://electricgrids.engr.tamu.edu/electric-grid-test-cases/ieee-14-bus-system/>`__
on 2019-09-25.

tx2000
------

This directory contains a variant of the Texas 2000 bus synthetic 
grid. You can find a version of this model (though possibly not
identical) on Texas A&M's `website <https://electricgrids.engr.tamu.edu/electric-grid-test-cases/>`__.
This particular case was provided by `Adam Birchfield <http://adambirchfield.com/>`__.
on 2019-11-21.

tx2000_mod
----------

This directory contains a variant of the Texas 2000 bus synthetic grid
modified for voltage control experiments. This case was provided by
Diana Wallison (diwalli@tamu.edu).

wscc_9
------

This directory contains an approximation of the Western System
Coordinating Council (WSCC) system. This model is commonly used for
dynamics simulations. It was `downloaded from Texas A&M
<https://electricgrids.engr.tamu.edu/electric-grid-test-cases/wscc-9-bus-system/>`__
on 2020-04-01.
docs
====
Directory for documentation, etc.

Files
-----

.nojekyll
^^^^^^^^^

TODO: Why do we have this?

Auxiliary File Format.pdf
^^^^^^^^^^^^^^^^^^^^^^^^^

This file was obtained from PowerWorld Simulator 21 on 2019-10-15 by
clicking the "Window" tab, then clicking the "Auxiliary File Format"
button under the "Auxiliary Files" section.

index.html
^^^^^^^^^^

Used for hosting the documentation online.

power_world_object_fields.xlsx
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Created by using PowerWorld Simulator to export case object fields.
This can be accomplished by navigating to the `Window` tab, clicking
on `Export Case Object Fields`, and selecting `Send to Excel`. This
file was generated with PowerWorld Simulator Version 21.

Directories
-----------

doctrees
^^^^^^^^

Directory used by sphinx when generating documentation.

html
^^^^

Directory used by sphinx when generating documentation

joss
^^^^

Directory for storing JOSS paper files.

rst
^^^

Contains all ESA documentation in reStructuredText (.rst) format.
Current Directions
------------------
Use the `paper preview service <https://whedon.theoj.org/>`__ from JOSS to compile the paper.

Files:
^^^^^^

- paper.bib: BibTex file for paper references.
- paper.md: The actual paper, in MarkDown format.
- paper.pdf: Current compiled .pdf version of the paper.
- README.rst: This file.

Old Directions
--------------
Please note the directions below are no longer applicable. JOSS does not
allow for modification as mentioned in `this comment <https://github.com/openjournals/joss-reviews/issues/2289#issuecomment-642588057>`__.

The original latex.template is from `the whedon
repository <https://github.com/openjournals/whedon/tree/master/resources
/joss>`__.

To compile the paper in Windows, you will need to install the Pandoc
and the MiKTex. When you are installing the MiKTex, make sure to check
the option to "install the packages on-the-fly".

After the installation, restart your computer.

When the Pandoc and the MiKTex are ready, you can now compile the paper
using the command:

.. code:: bat

    pandoc --filter pandoc-citeproc --bibliography paper.bib paper.md --template latex.template -o paper.pdf --pdf-engine=xelatex

One thing worth to note is, when we submit the paper to JOSS, we need to
let them know that we have our own custom latex.template file for
compiling. Additionally, before JOSS submittal the commented-out line
related to the logo needs to be "un-commented." The line looks like:

.. code:: latex

    %\fancyhead[L]{\hspace{-0.75cm}\includegraphics[width=5.5cm]{$logo_path$}}
Common Issues
=============
This section will describes some (maybe) common issues and their
solutions. If you encounter and solve an issue, please file a `GitHub
issue <https://github.com/mzy2240/ESA/issues>`__ so that we can add your
problem and solution to this section.

Before diving too deeply into the issues listed here, first ensure that
you have all the prerequisite software installed (including PowerWorld
Simulator and the SimAuto add-on) and are using a supported version of
Python (>= 3.5).

.. _venv-issues:

Installation/Virtual Environment Issues
---------------------------------------

If you have issues installing ESA and/or its dependencies, you may need
to do some manual work installing prerequisites in your virtual
environment. Hopefully following these simple directions will help fix
most issues.

1. Start fresh! Completely remove your virtual environment and recreate
   it. `PyCharm makes this pretty
   easy <https://www.jetbrains.com/help/pycharm/creating-virtual-environment.html>`__,
   or you can do so manually using `Python's
   guide <https://docs.python.org/3/tutorial/venv.html>`__. The
   remaining directions will assume you're typing commands into your
   **activated** virtual envrionment.
2. Reinstall pip and setuptools:
   ``python -m pip install --upgrade --force-reinstall pip setuptools``.
   We're intentionally using ``python -m pip`` instead of just ``pip``
   to avoid possible path issues. Note that you might need to run this
   command twice (the first may fail for some unknown reason).
3. Check out ESA's
   `setup.py <https://github.com/mzy2240/ESA/blob/master/setup.py>`__
   file and look for ``install_requires``. It'll look something like
   ``['pandas', 'numpy', 'pywin32', 'pypiwin32']``.
4. Using what we found under ``install_requires``, install ESA's
   dependencies manually. To avoid compiler dependencies, we'll get
   binary distributions only:
   ``python -m pip install --upgrade --only-binary :all: pandas numpy pywin32 pypiwin32``

   -  If this command fails, you may need to pick and choose which
      dependencies you grab binary distributions for, and which you get
      other types of distributions for. Here's the `Python
      documentation <https://pip.pypa.io/en/stable/reference/pip_install/>`__.
      As a strictly illustrative example, if we only want to get binary
      distributions for ``pandas`` and ``numpy``, we'd modify the
      previous command to instead read like so:
      ``python -m pip install --upgrade --only-binary pandas,numpy pandas numpy pywin32 pypiwin32``

   - The authors of ESA have at times had issues installing pywin32 and
     pypiwin32 when *not* using the ``--only-binary`` option. So, if
     you're encountering errors you suspect are related to pywin32,
     try to uninstall and reinstall pywin32 and pypiwin32 with the
     ``--only-binary`` option.

5. After you've installed ESA's dependencies, it's time to install ESA:
   ``python -m pip install esa``

PyCharm Virtual Environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you use PyCharm to automatically create virtual environments for you,
there's a little extra work to do to get everything working for Python
3.8 (and possibly for other Python versions as well). Start up a
terminal *inside* PyCharm (click on the ``Terminal`` button which
defaults to the lower left area). In the terminal, run:
``python -m pip install -U --force-reinstall pip``. Note you may need to
run this command twice - mine failed the first time. The same may be
required for ``setuptools`` and/or ``distutils``.

Errors/Issues Initializing a SAW Instance
-----------------------------------------

This section will cover some common issues when attempting to initialize
a SAW instance. The first thing to check is that your arguments are
correct - check the API documentation first.

esa.saw.PowerWorldError: OpenCase: Errors have occurred
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may see an error stack trace that looks something like the
following:

.. code:: python

   Traceback (most recent call last):
     File "<input>", line 1, in <module>
     File "C:\Users\myuser\git\ESA\esa\saw.py", line 111, in __init__
       self.OpenCase(FileName=FileName)
     File "C:\Users\myuser\git\ESA\esa\saw.py", line 680, in OpenCase
       return self._call_simauto('OpenCase', self.pwb_file_path)
     File "C:\Users\myuser\git\ESA\esa\saw.py", line 1101, in _call_simauto
       raise PowerWorldError(output[0])
   esa.saw.PowerWorldError: OpenCase: Errors have occurred

Often, this is due to a bad path specification. Ensure you're providing
a **full** file path, including the file extension (.pwb), and that the
file exists at the exact path you specified.

Also, make sure that the
file is **actually** a PowerWorld binary file. If you open the file with
a text editor and see a bunch of weird symbols that are unintelligible
to a mere mortal, it's likely a PowerWorld binary file. If, upon opening
the file you see something like:

..

    | version https://git-lfs.github.com/spec/v1
    | oid sha256:f05131d24da96daa6a6712c5b9d368c81eeaea5dc7d0b6c7bec7d03ccf021b4a
    | size 34

Then you're looking at a Git LFS pointer file, and likely need to
install `Git LFS <https://git-lfs.github.com/>`__ and perform a
``git lfs pull``.

TypeError: This COM object can not automate the makepy process - please run makepy manually for this object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you see an error like the above, try initializing your SAW object
again but set ``early_bind=False``. While we're unsure of the root cause
of this issue, it seems to be related to the fact that
``early_bind=True`` preemptively creates some Python files related to
the SimAuto COM API, and file permission issues can crop up.

AttributeError: module 'win32com.gen_py.C99F1760-277E-11D5-A106-00C04F469176x0x20x0' has no attribute 'CLSIDToClassMap'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you see an error like that listed above, it's possible the pywin32
cache has been somehow corrupted (perhaps your computer crashed while
a script which uses ESA was running). Simply delete the following
directory (the default, you may have to adapt for your system):

``C:\Users\<your user directory>\AppData\Local\Temp\gen_py``

The key part here is ``gen_py``. If the above path isn't right for you,
use Windows to search for ``gen_py``.

ModuleNotFoundError: no module pywintypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you see an error like that listed above, try reinstalling pywin32
and pypiwin32 with the ``--only-binary`` option, as described in the
:ref:`venv-issues` section.

esa.saw.PowerWorldError: Access Violation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you see an error like that listed above, you probably could observe
the same errors using script functionality in the PowerWorld Simulator
interface. When you start the simulator (either the GUI or ESA) in a
remote desktop environment, due to security reasons the system may block
the simulator to save or create any files. As a result, any functions
that require the simulator to generate files will fail with such errors.
There are not much things we could do from our side, but one possible
hack is to login locally and use the simulator first, then use the
remote desktop to continue your work... _api:

esa API Documentation
=====================
ESA is a single module package - the only Python file is ``saw.py``.
To make it easy to remember, "SAW" is an acronym for "SimAuto Wrapper."
Within the ``saw`` module is the ``SAW`` class, which is the workhorse
of ESA.


esa Package
-----------

.. automodule:: esa
   :members:
   :undoc-members:
   :show-inheritance:

.. _esa-saw-api:

esa.saw Module
--------------

.. automodule:: esa.saw
   :members:
   :undoc-members:
   :show-inheritance:
.. _examples:

Examples
========
For your convenience, several ESA usage examples are provided in this
section. Some examples may also be presented earlier in the
documentation. If you would like to share your own examples with us,
please file an issue on `GitHub
<https://github.com/mzy2240/ESA/issues>`__. If you'd like even more
examples (though in a *much* less reader friendly format), check out
ESA's `tests <https://github.com/mzy2240/ESA/tree/master/tests>`__,
specifically the ``test_saw.py`` file.

ESA Quick Start
---------------
The following example is also presented in the :ref:`quick-start`
section.

.. include:: snippets/quick_start_14.rst

.. _increase-loading:

Increase Loading in Case
------------------------

The following example is also presented in the :ref:`what-is-esa`
section.

.. include:: snippets/increase_loading_14.rst

Add Lines to Case
-----------------

.. include:: snippets/add_lines_2000.rst

Transient Stability Analysis
----------------------------

.. include:: snippets/ts_mycontingency_200.rst

Create Simple Graph Model
-------------------------

.. include:: snippets/simple_graph_2000.rst

Created Graph Model with Edges Weighted by Branch Impedance
-----------------------------------------------------------

.. include:: snippets/weighted_graph_14.rst

Plot Histogram of Line Flows with Matplotlib
--------------------------------------------

.. include:: snippets/line_loading_matplotlib_2000.rst
.. _welcome:

Welcome!
========

Welcome to Easy SimAuto's (ESA) documentation! If you encounter issues
with ESA or the documentation, please file an issue on `GitHub
<https://github.com/mzy2240/ESA/issues>`__.

.. _citation:

Citation
--------

.. include:: citation.rst

.. _what-is-esa:

What Is ESA?
------------

.. include:: welcome/what.rst

.. _about-simauto:

About SimAuto
-------------

.. include:: welcome/about_simauto.rst

Who Should Use ESA?
-------------------

.. include:: welcome/who.rst

Why Use ESA?
------------

.. include:: welcome/why.rst

Projects That Use ESA
---------------------

.. include:: welcome/projects.rst

Disclaimer and Credits
----------------------

Disclaimer
^^^^^^^^^^

.. include:: welcome/disclaimer.rst

Credits and Thanks
^^^^^^^^^^^^^^^^^^

.. include:: welcome/credits.rst

.. _changelog:

Changelog
---------

.. include:: welcome/changelog.rst
If you use ESA in any of your work, please use the citation below.

.. code:: latex

    @article{ESA,
      doi = {10.21105/joss.02289},
      url = {https://doi.org/10.21105/joss.02289},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {50},
      pages = {2289},
      author = {Brandon L. Thayer and Zeyu Mao and Yijing Liu and Katherine Davis and Thomas J. Overbye},
      title = {Easy SimAuto (ESA): A Python Package that Simplifies Interacting with PowerWorld Simulator},
      journal = {Journal of Open Source Software}
    }
.. _quick-start:

Quick Start
===========
The following quick start example uses the IEEE 14 bus case, which can
be found `in the ESA
repository <https://github.com/mzy2240/ESA/tree/master/tests/cases/ieee_14>`__
or from `Texas
A&M University
<https://electricgrids.engr.tamu.edu/electric-grid-test-cases/ieee-14-bus-system/>`__.
Notes related to the case are available from the `University of Washington
<https://labs.ece.uw.edu/pstca/>`__.

You can find API documentation
`here <https://mzy2240.github.io/ESA/html/esa.html>`__ and more examples
`here <https://mzy2240.github.io/ESA/html/snippets.html>`__.

.. include:: snippets/quick_start_14.rst

For more examples, you can find them `here <https://mzy2240.github
.io/ESA/html/snippets.html>`__... table:: ESA's testing coverage as of 2022-01-28 (Git commit: b8cef3a)
    :widths: auto
    :align: left

    +-----------------+-------------------+-----------------+-----------------+--------------------+
    | Name            |   Num. Statements |   Missing Lines |   Covered Lines |   Percent Coverage |
    +=================+===================+=================+=================+====================+
    | esa/__init__.py |                 2 |               0 |               2 |                100 |
    +-----------------+-------------------+-----------------+-----------------+--------------------+
    | esa/saw.py      |               726 |               0 |             726 |                100 |
    +-----------------+-------------------+-----------------+-----------------+--------------------+
.. _installation:

Installing ESA
==============
Installing ESA is simple! Follow the directions below to get started.

Overview
--------

.. include:: installation/overview.rst

Prerequisites
--------------

.. include:: installation/pre_reqs.rst

Virtual Environment Configuration
---------------------------------

.. include:: installation/virtualenv.rst

Install Prerequisite Packages
-----------------------------

.. include:: installation/manual_pywin32.rst

Install ESA via Pip
-------------------

.. include:: installation/pip.rst

Install ESA from Source
-----------------------

.. include:: installation/source.rst

Post-Installation
-----------------
The post-installation steps found here are optional, but may save you
some headache down the line.

Cursory Verification of Successful Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: installation/postinstall_verify.rst

Execute ESA Unittests
^^^^^^^^^^^^^^^^^^^^^

.. include:: installation/postinstall_unittests.rst

pywin32 Post-Installation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: installation/postinstall_pywin32.rst

ESA Overview
============
Please see :ref:`about-simauto` for a description of SimAuto as well as
for links to PowerWorld's documentation.

When using ESA, you'll likely only use a single class: ``SAW``, which is
short for SimAuto Wrapper. After you've installed the esa package
(refer to :ref:`installation`), you'll import the ``SAW`` class like so:

.. code-block:: python

    from esa import SAW
    saw = SAW('<full path to case file>')

All methods of the SAW class are fully documented. Additionally, all
inputs and outputs use type hinting so that your interactive development
environment (IDE, such as PyCharm) will automatically highlight
incorrect input types.

Naming Conventions
------------------

When browsing the documentation, you'll notice that naming
conventions are mixed, with some being CamelCase and others being
lower_case_with_underscores. We use the CamelCase convention to indicate
that we're more or less directly calling a
`PowerWorld SimAuto Function <https://www.powerworld.com/WebHelp/#MainDocumentation_HTML/Simulator_Automation_Server_Functions.htm%3FTocPath%3DAutomation%2520Server%2520Add-On%2520(SimAuto)%7CAutomation%2520Server%2520Functions%7C_____3>`__
(e.g. ``SAW.ChangeParametersMultipleElement``) or a
`PowerWorld Script Command <https://www.powerworld.com/WebHelp/Content/MainDocumentation_HTML/Auxiliary_Files.htm>`__
(e.g. ``SAW.SolvePowerFlow``). CamelCase is also used when describing
variables which are more or less directly passed into a SimAuto
function (e.g. ``ObjectType``). The lower_case_with_underscores
convention is used in function naming to indicate that a function is a
higher-level function which may call multiple PowerWorld functions
and/or not return everything from PowerWorld. Some examples of these
functions are
``SAW.change_and_confirm_params_multiple_element``,
``SAW.change_parameters_multiple_element_df``,
``SAW.get_key_fields_for_object_type``, and
``SAW.get_power_flow_results``. In general, it is recommended to use
these higher level functions where possible. Note these show up toward
the bottom of the API documentation since methods which start with
upper case letters come before methods that start with lower case
letters. Variables use the lower_case_with_underscores convention any
time the variable is not a direct SimAuto input.

Functions/Methods
-----------------

SimAuto Functions
^^^^^^^^^^^^^^^^^

The ``SAW`` class has every SimAuto function implemented. I.e., you
can call a ``SAW`` method corresponding to every documented `SimAuto
function <https://www.powerworld.com/WebHelp/#MainDocumentation_HTML/Simulator_Automation_Server_Functions.htm%3FTocPath%3DAutomation%2520Server%2520Add-On%2520(SimAuto)%7CAutomation%2520Server%2520Functions%7C_____3>`__.


High-Level/Helper Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned in the `Naming Conventions`_ section, we recommend you use
the high level helper functions (which use the
lower_case_with_underscores convention) where possible.

.. _runscriptcommand:

RunScriptCommand
^^^^^^^^^^^^^^^^

For PowerWorld functionality not directly covered by the SimAuto
functions, you'll want to use ``RunScriptCommand``. Note that we
already have some methods defined which are script commands, e.g.
``SAW.SolvePowerFlow``. So, you may want to search for the function
you want before directly using ``RunScriptCommand``.

It's worth noting that ``RunScriptCommand`` will directly return results
from PowerWorld, which will come back with all sorts of nasty data types
(e.g. strings with leading and trailing whitespace). Your best course of
action is to create a DataFrame/Series from the output, and use the
``clean_df_or_series`` method afterwards.

Documentation from PowerWorld on available script commands can be found
`here
<https://github.com/mzy2240/ESA/blob/master/docs/Auxiliary%20File%20Format.pdf>`__.

clean_df_or_series
^^^^^^^^^^^^^^^^^^

This helper function will do automatic type translation for you based
on known PowerWorld data types. If you're dealing with direct outputs
from PowerWorld (e.g. from using ``RunScriptCommand``), this method
will save you all sorts of trouble. Read the method's API documentation
thoroughly before using.

Data Types
----------

All method input and output data types are documented in the API
documentation. Where possible, ``SAW`` methods return Pandas DataFrames
or Pandas Series. If there's nothing to return, ``None`` will be
returned. ESA makes extensive use of type hinting so that your IDE can
automatically highlight issues related to data types.

.. _powerworld-variables:

PowerWorld Variables
--------------------

At present, ESA uses PowerWorld "Variable Names" as opposed to
PowerWorld "Concise Variable Names." A listing of these variables can be
found `here
<https://github.com/mzy2240/ESA/blob/master/docs/power_world_object_fields.xlsx>`__.
It would seem that PowerWorld is moving toward "Concise Variable Names,"
and in a future update ESA may support these (see `this issue
<https://github.com/mzy2240/ESA/issues/1#issue-525219427>`__).

Testing Coverage
----------------

The ESA team strives to write good tests with 100% coverage. The table
below provides the latest test coverage data for ESA.

.. include:: coverage.rst

Contributing
------------

We welcome contributions to ESA - please give
`contributing.md <https://github.com/mzy2240/ESA/blob/master/contributing.md>`__
a read... ESA documentation master file, created by
   sphinx-quickstart on Tue Nov 12 14:40:53 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Easy SimAuto (ESA) Documentation
================================

Easy SimAuto is a Python package to ease interaction with PowerWorld
Simulator's :ref:`SimAuto <about-simauto>` add-on. This documentation
seeks to be comprehensive, with detailed :ref:`installation instructions
<installation>`, complete :ref:`API documentation <api>`, and a rich
set of :ref:`examples <examples>`. If you're new to ESA, please start
by reading the :ref:`welcome` section. Returning users looking to see
what's new should check out the :ref:`changelog`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   welcome
   installation
   quick_start
   overview
   esa
   snippets
   common_issues


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
rst
====

This directory contains all the files needed to build ESA's
documentation, which is in  `reStructuredText
<http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__
(.rst) format.

Building the Documentation
--------------------------

First, ensure you've installed all pre-requisite packages for building
the documentation. These pre-requisites can be found in ``setup.py`` at
the repository's top-level, under the ``extras_require`` argument to
``setuptools.setup``.

Next, activate your virtual environment, and change directories to here.

Finally, simply execute ``make html`` to build the documentation.

Initial One Time Setup (DO NOT RUN THIS)
----------------------------------------

The following is just for recording what was done to set things up
initially. All commands were run with an activated virtual environment
within the sphinx directory. Some information may now be out of date,
but this at least gives a clue as to how things were done.

#.  Run ``sphinx-quickstart``
#.  Run ``sphinx-apidoc -o . ../esa``
#.  Add ``esa`` and `modules`` lines to ``index.rst``.
#.  Add ``'spinx.ext.autodoc'`` extension in conf.py.
#.  Uncomment and modify the ``# -- Path setup`` section of conf.py

Files
-----

All ``.rst`` files are used in creating ESA's documentation. The "main"
file is ``index.rst``.

conf.py
^^^^^^^
Sphinx configuration file.

coverage.rst
^^^^^^^^^^^^

Testing coverage table created by ``coverage_to_rst.py``. Included by
overview.rst.

coverage_to_rst.py
^^^^^^^^^^^^^^^^^^

Runs ESA unittests and assesses testing coverage. Modifies top level
README file and generates coverage.rst.

make.bat
^^^^^^^^

Windows batch script used to kick off the documentation build.

Makefile
^^^^^^^^

Makefile for the documentation build.

Directories
-----------

This section discusses all directories contained in this directory.
Please keep it up to date.

_static
^^^^^^^

Store static files here.

_templates
^^^^^^^^^^

Store templates here.

installation, snippets, welcome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These directories contain files that are included to the corresponding
.rst files. Note that the ``snippets`` directory is special, as all
examples in the files get executed as part of ESA's testing suite.
Please see the README within the ``snippets``At this point in time, one can only run unittests if ESA is installed
from source and if you've installed the ``test`` dependencies. A
future version of ESA may make tests available when installing via Pip.

Using a Command Prompt window, change directories to the ESA repository
and run the tests like so (adapt paths as necessary):

.. code:: bat

    cd C:\Users\myuser\git\ESA
    my-venv\Scripts\activate.bat
    python -m unittest discover tests

During the running of the tests, you'll see some logging output and
error messages - this is expected. What's important is that at the end,
you see a message like:

    | Ran 73 tests in 34.542s
    |
    | OK (expected failures=2)

If something is wrong, you'll see some indications of failure. The
"expected failures" are okay, and do not indicate there are any issues.

If you've installed the testing dependencies listed in setup.py, you
should have "coverage.py" installed. If you'd like to assess ESA's
testing coverage, the incantation looks like (run from top-level of
ESA repository after activating virtual environment):

.. code:: bat

    coverage run
    coverage report -m

Note that the arguments to ``coverage run`` are provided by
``.coveragerc`` at the top-level of the repository.
Like any Python project, the use of virtual environments is **strongly**
encouraged. If you're new to virtual environments, Python provides a
nice `tutorial <https://docs.python.org/3/tutorial/venv.html>`__.

After installing Python on your machine, open up a Command Prompt window
(Hit the "Windows" key and the "R" key simultaneously, type in "cmd" in
the popup, then hit Enter) then set up a virtual environment like so
(adapt paths as necessary):

.. code:: bat

    cd C:\path\to\my\project
    C:\path\to\python.exe -m venv my-venv

Then, activate the virtual environment like so (in the same terminal):

.. code:: bat

    my-venv\Scripts\activate.bat

**All installation directions assume your virtual environment has been
activated.**

Next, update pip and setuptools (in your activated virtual environment):

.. code:: bat

    python -m pip install --upgrade --force-reinstall pip setuptoolsESA has the following prerequisites:

*   Microsoft Windows Operating System (PowerWorld is Windows only. The
    authors of ESA have only ever attempted to use Windows 10.).
*   Python >=3.5. Download it `here
    <https://www.python.org/downloads/>`__.
*   `PowerWorld <https://www.powerworld.com/>`__ Simulator with
    the `Automation Server (SimAuto) add-on
    <https://www.powerworld.com/products/simulator/add-ons-2/simauto>`__
    installed.

    * **NOTE**: the authors of ESA have tested with Simulator
      versions 17 and 21. It is likely, but **not guaranteed**, that ESA
      will work with all Simulator versions 16-21. If you encounter a
      problem with a particular version, please file an `issue
      <https://github.com/mzy2240/ESA/issues>`__ and we may be able
      to help (if we can get access to that particular Simulator
      version).

*   `Git Large File Storage (LFS) <https://git-lfs.github.com/>`__
    (**OPTIONAL**: required to download case files and run tests). After
    installing Git LFS, simply change directories to the ESA repository,
    and run ``git lfs install``. You will likely need to run a
    ``git pull`` or ``git lfs pull`` after installing and setting up Git
    LFS. After initial setup, you shouldn't need to do anything else
    with Git LFS.Installing ESA from source is quite simple, and does not require any
extra tools. One might choose to install ESA from source if you plan
to extend or modify ESA, or would like to run ESA's suite of tests.

The first step is to obtain a copy of the ESA repository. There are two
ways to obtain a copy: download a .zip archive or use Git to clone the
repository. `Here's a link to ESA's GitHub repository.
<https://github.com/mzy2240/ESA>`__ On the GitHub page you can find a
link to "Clone or download."

If you choose to clone the repository with Git, you'll also need to
install `Git Large File Storage <https://git-lfs.github.com/>`__. After
installation, use Git Bash or a Command Prompt window to execute the
following:

.. code:: bat

    cd C:\Users\myuser\git\ESA
    git lfs install
    git pull
    git lfs pull

If you choose to download a .zip archive, you'll of course need to
extract the archive to the desired directory.

Once you have a copy of the ESA repository, simply run the following
in a Command Prompt window *after* activating your virtual environment
(adapt path to match your setup):

.. code:: bat

    cd C:\Users\myuser\git\ESA
    python -m pip install .

If you would like to be able to run tests, you'll need to install the
testing dependencies:

.. code:: bat

    python -m pip install .[test]

Similarly, to build the documentation, you'll need:

.. code:: bat

    python -m pip install .[doc]

If you want to both run tests and build the documentation:

.. code:: bat

    python -m pip install .[test,doc]

You can find the specified "extras" in setup.py - look for
``extras_require``.
To install ESA with pip, run the following in your activated virtual
environment:

.. code:: bat

    python -m pip install esa
TL;DR:

.. code:: bat

    python -m pip install --only-binary pywin32,pypiwin32 pywin32 pypiwin32 esa

Installing ESA is easy, and most users will simply want to use Python's
package manager, Pip, to install ESA. Users who want to extend or modify
ESA should install from source.The simplest way to verify that ESA installed correctly is to attempt
to import it. Most (but not all) installation issues will rear their
heads at this point. Simply execute the following in a Command Prompt
window (ensure your virtual environment is activated!):

.. code:: bat

    python -c "import esa; print('Success!')"

If an error message is emitted, ESA or its dependencies are not properly
installed. If you're unable to figure it out on your own, feel free to
file an issue on `GitHub <https://github.com/mzy2240/ESA/issues>`__. We
do not guarantee that we can help everyone.
From the `pywin32 GitHub page <https://github.com/mhammond/pywin32>`__:

    Note that if you want to use pywin32 for "system wide" features,
    such as registering COM objects or implementing Windows Services,
    then you must run the following command from an elevated command
    prompt:

        python Scripts/pywin32_postinstall.py -install

Note that the aforementioned "Scripts" directory can be found inside
your virtual environment directory.
ESA depends on the pywin32 Python package in order to interface with
Windows components. Unfortunately, pywin32 does not always cleanly
install automatically when installing ESA, so it's recommended that you
first manually install pywin32. In your activated virtual environment,
execute the following:

.. code:: bat

    python -m pip install --only-binary :all: pypiwin32 pywin32

The authors have found that the ``--only-binary`` flag is often
necessary to get pywin32 to work - without it, pywin32 is unable to find
some necessary libraries.

If you are using a conda environment, a simpler way is to use conda
to install instead of pip:

.. code:: bat

    conda install pywin32This example shows how to easily transform a grid model into a graph
supported by `NetworkX <https://networkx.github.io/>`__. NetworkX is a
popular Python package for analyzing graph structure, building network
models and designing new network algorithms. You'll first need to
install NetworkX into your virtual environment (which should be
activated!), which is most easily done by:

.. code:: bat

    python -m pip install networkx

Before following along with the example, define the ``CASE_PATH``
constant (the file path to a PowerWorld ``.pwb`` case file) like so,
adapting the path to your system:

.. code:: python

    CASE_PATH = r"C:\Users\myuser\git\ESA\tests\cases\tx2000\tx2000_base_pws_version_21.pwb"

On to the example!

Perform imports, initialize a ``SAW`` instance:

.. code:: python

    >>> from esa import SAW
    >>> import pandas as pd
    >>> import networkx as nx
    >>> saw = SAW(CASE_PATH, early_bind=True)

Get a DataFrame with all branches (lines, transformers, etc.):

.. code:: python

    >>> kf = saw.get_key_field_list('branch')
    >>> kf
    ['BusNum', 'BusNum:1', 'LineCircuit']
    >>> branch_df = saw.GetParametersMultipleElement('branch', kf)
    >>> branch_df
          BusNum  BusNum:1 LineCircuit
    0       1001      1064           1
    1       1001      1064           2
    2       1001      1071           1
    3       1001      1071           2
    4       1002      1007           1
    ...      ...       ...         ...
    3199    8157      5124           1
    3200    8157      8156           1
    3201    8158      8030           1
    3202    8159      8158           1
    3203    8160      8159           1
    <BLANKLINE>
    [3204 rows x 3 columns]

To learn more about variables such as ``LineCircuit``, see
:ref:`powerworld-variables`.

Create the graph from the DataFrame. Yes, it is this simple. Use
``Graph`` instead of ``MultiGraph`` if there are no parallel branches.

.. code:: python

    >>> graph = nx.from_pandas_edgelist(branch_df, "BusNum", "BusNum:1", create_using=nx.MultiGraph)
    >>> graph.number_of_nodes()
    2000
    >>> graph.number_of_edges()
    3204

Clean up:

.. code:: python

    saw.exit()
This simple example uniformly increases the loading in a power system
model by 50%.

If you want to follow along, you'll first need to define your own
``CASE_PATH`` constant (the file path to a PowerWorld ``.pwb`` case
file), like so (adapt the path for your system):

.. code:: python

    CASE_PATH = r"C:\Users\myuser\git\ESA\tests\cases\ieee_14\IEEE 14 bus_pws_version_21.pwb"

Then, import the SimAuto wrapper (SAW) class and initialize an instance:

.. code:: python

    >>> from esa import SAW
    >>> saw = SAW(CASE_PATH)

Retrieve key fields for loads:

.. code:: python

    >>> kf = saw.get_key_field_list('load')
    >>> kf
    ['BusNum', 'LoadID']

Pull load data including active and reactive power demand:

.. code:: python

    >>> load_frame = saw.GetParametersMultipleElement('load', kf + ['LoadSMW', 'LoadSMVR'])
    >>> load_frame
        BusNum LoadID    LoadSMW   LoadSMVR
    0        2      1  21.699999  12.700000
    1        3      1  94.199997  19.000000
    2        4      1  47.799999  -3.900000
    3        5      1   7.600000   1.600000
    4        6      1  11.200000   7.500000
    5        9      1  29.499999  16.599999
    6       10      1   9.000000   5.800000
    7       11      1   3.500000   1.800000
    8       12      1   6.100000   1.600000
    9       13      1  13.500001   5.800000
    10      14      1  14.900000   5.000000

To learn more about variables such as ``LoadSMW``, see
:ref:`powerworld-variables`.

Uniformly increase loading by 50% and solve the power flow:

.. code:: python

    >>> load_frame[['LoadSMW', 'LoadSMVR']] *= 1.5
    >>> saw.change_parameters_multiple_element_df('load', load_frame)
    >>> saw.SolvePowerFlow()

Let's confirm that the loading did indeed increase:

.. code:: python

    >>> new_loads = saw.GetParametersMultipleElement('load', kf + ['LoadSMW', 'LoadSMVR'])
    >>> new_loads
        BusNum LoadID     LoadSMW   LoadSMVR
    0        2      1   32.549998  19.050001
    1        3      1  141.299999  28.500000
    2        4      1   71.699995  -5.850000
    3        5      1   11.400000   2.400000
    4        6      1   16.800001  11.250000
    5        9      1   44.250000  24.900000
    6       10      1   13.500001   8.700000
    7       11      1    5.250000   2.700000
    8       12      1    9.150000   2.400000
    9       13      1   20.250002   8.700000
    10      14      1   22.350000   7.500000

Clean up when done:

.. code:: python

    >>> saw.exit()

Easy, isn't it?This "quick start" example has several purposes:

*   Illustrate how ESA is used with a simple power system model.

*   Demonstrate how to perform common tasks (e.g. solving the power
    flow, retrieving simulation data such as bus voltages and power
    injections).

*   Show the usefulness of some of ESA's high level helper functions
    which do more than simply wrap SimAuto functions.

Before running the example below, define a CASE_PATH constant (the file
path to a PowerWorld ``.pwb`` case file) like so (adapt the path as
needed for your system):

.. code:: python

    CASE_PATH = r"C:\Users\myuser\git\ESA\tests\cases\ieee_14\IEEE 14 bus_pws_version_21.pwb"

On to the quick start!

Start by Importing the SimAuto Wrapper (SAW) class:

.. code:: python

   >>> from esa import SAW


Initialize SAW instance using 14 bus test case:

.. code:: python

   >>> saw = SAW(FileName=CASE_PATH)

Solve the power flow:

.. code:: python

   >>> saw.SolvePowerFlow()

Retrieve power flow results for buses. This will return a Pandas
DataFrame to make your life easier.

.. code:: python

    >>> bus_data = saw.get_power_flow_results('bus')
    >>> bus_data
        BusNum BusName  BusPUVolt   BusAngle    BusNetMW  BusNetMVR
    0        1   Bus 1   1.060000   0.000000  232.391691 -16.549389
    1        2   Bus 2   1.045000  -4.982553   18.300001  30.855957
    2        3   Bus 3   1.010000 -12.725027  -94.199997   6.074852
    3        4   Bus 4   1.017672 -10.312829  -47.799999   3.900000
    4        5   Bus 5   1.019515  -8.773799   -7.600000  -1.600000
    5        6   Bus 6   1.070000 -14.220869  -11.200000   5.229700
    6        7   Bus 7   1.061520 -13.359558    0.000000   0.000000
    7        8   Bus 8   1.090000 -13.359571    0.000000  17.623067
    8        9   Bus 9   1.055933 -14.938458  -29.499999   4.584888
    9       10  Bus 10   1.050986 -15.097221   -9.000000  -5.800000
    10      11  Bus 11   1.056907 -14.790552   -3.500000  -1.800000
    11      12  Bus 12   1.055189 -15.075512   -6.100000  -1.600000
    12      13  Bus 13   1.050383 -15.156196  -13.500001  -5.800000
    13      14  Bus 14   1.035531 -16.033565  -14.900000  -5.000000

Retrieve power flow results for generators:

.. code:: python

    >>> gen_data = saw.get_power_flow_results('gen')
    >>> gen_data
       BusNum GenID       GenMW     GenMVR
    0       1     1  232.391691 -16.549389
    1       2     1   40.000001  43.555957
    2       3     1    0.000000  25.074852
    3       6     1    0.000000  12.729700
    4       8     1    0.000000  17.623067


To learn more about variables such as ``GenMW``, see
:ref:`powerworld-variables`.

Let's change generator injections! But first, we need to know which
fields PowerWorld needs in order to identify generators. These fields
are known as key fields.

.. code:: python

    >>> gen_key_fields = saw.get_key_field_list('gen')
    >>> gen_key_fields
    ['BusNum', 'GenID']


Change generator active power injection at buses 3 and 8 via SimAuto
function:

.. code:: python

    >>> params = gen_key_fields + ['GenMW']
    >>> values = [[3, '1', 30], [8, '1', 50]]
    >>> saw.ChangeParametersMultipleElement(ObjectType='gen', ParamList=params, ValueList=values)


Did changing generator active power injections work? Let's confirm:

.. code:: python

    >>> new_gen_data = saw.GetParametersMultipleElement(ObjectType='gen', ParamList=params)
    >>> new_gen_data
       BusNum GenID       GenMW
    0       1     1  232.391691
    1       2     1   40.000001
    2       3     1   30.000001
    3       6     1    0.000000
    4       8     1   50.000000


It would seem the generator active power injections have changed. Let's
re-run the power flow and see if bus voltages and angles change.
Spoiler: they do.

.. code:: python

    >>> saw.SolvePowerFlow()
    >>> new_bus_data = saw.get_power_flow_results('bus')
    >>> cols = ['BusPUVolt', 'BusAngle']
    >>> diff = bus_data[cols] - new_bus_data[cols]
    >>> diff
           BusPUVolt   BusAngle
    0   0.000000e+00   0.000000
    1  -1.100000e-07  -2.015596
    2  -5.700000e-07  -4.813164
    3  -8.650700e-03  -3.920185
    4  -7.207540e-03  -3.238592
    5  -5.900000e-07  -4.586528
    6  -4.628790e-03  -7.309167
    7  -3.190000e-06 -11.655362
    8  -7.189370e-03  -6.284631
    9  -6.256150e-03  -5.987861
    10 -3.514030e-03  -5.297895
    11 -2.400800e-04  -4.709888
    12 -1.351040e-03  -4.827348
    13 -4.736110e-03  -5.662158


Wouldn't it be easier if we could change parameters with a DataFrame?
Wouldn't it be nice if we didn't have to manually check if our updates
were respected? You're in luck!

Create a copy of the ``gen_data`` DataFrame so that we can modify its
values and use it to update parameters in PowerWorld. Then, change the
generation for the generators at buses 2, 3, and 6.

.. code:: python

    >>> gen_copy = gen_data.copy(deep=True)
    >>> gen_copy.loc[gen_copy['BusNum'].isin([2, 3, 6]), 'GenMW'] = [0.0, 100.0, 100.0]
    >>> gen_copy
       BusNum GenID       GenMW     GenMVR
    0       1     1  232.391691 -16.549389
    1       2     1    0.000000  43.555957
    2       3     1  100.000000  25.074852
    3       6     1  100.000000  12.729700
    4       8     1    0.000000  17.623067


Use helper function ``change_and_confirm_params_multiple_element`` to
both command the generators and to confirm that PowerWorld respected the
command. This is incredibly useful because if you directly use
``ChangeParametersMultipleElements``, PowerWorld may unexpectedly not
update the parameter you tried to change! If the following does not
raise an exception, we're in good shape (it doesn't)!

.. code:: python

   >>> saw.change_and_confirm_params_multiple_element(ObjectType='gen', command_df=gen_copy.drop('GenMVR', axis=1))

Run the power flow and observe the change in generation at the slack
bus (bus 1):

.. code:: python

    >>> saw.SolvePowerFlow()
    >>> new_gen_data = saw.get_power_flow_results('gen')
    >>> new_gen_data
       BusNum GenID       GenMW     GenMVR
    0       1     1   62.128144  14.986289
    1       2     1    0.000000  10.385347
    2       3     1  100.000000   0.000000
    3       6     1  100.000000  -3.893420
    4       8     1    0.000000  17.399502


What if we try to change generator voltage set points? Start by getting
a DataFrame with the current settings. Remember to always access the
key fields so that when we want to update parameters later PowerWorld
knows how to find the generators.

.. code:: python

    >>> gen_v = saw.GetParametersMultipleElement('gen', gen_key_fields + ['GenRegPUVolt'])
    >>> gen_v
       BusNum GenID  GenRegPUVolt
    0       1     1      1.060000
    1       2     1      1.045000
    2       3     1      1.025425
    3       6     1      1.070000
    4       8     1      1.090000

Now, change all voltage set points to 1 per unit:

.. code:: python

    >>> gen_v['GenRegPUVolt'] = 1.0
    >>> gen_v
       BusNum GenID  GenRegPUVolt
    0       1     1           1.0
    1       2     1           1.0
    2       3     1           1.0
    3       6     1           1.0
    4       8     1           1.0

    >>> saw.change_and_confirm_params_multiple_element('gen', gen_v)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "C:\Users\myuser\git\ESA\esa\saw.py", line 199, in change_and_confirm_params_multiple_element
        raise CommandNotRespectedError(m)
    esa.saw.CommandNotRespectedError: After calling ChangeParametersMultipleElement, not all parameters were actually changed within PowerWorld. Try again with a different parameter (e.g. use GenVoltSet instead of GenRegPUVolt).

So, PowerWorld didn't respect that command, but we've been saved from
future confusion by the ``change_and_confirm_params_multiple_element``
helper function.

Let's call the LoadState SimAuto function:

.. code:: python

    >>> saw.LoadState()
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "C:\Users\myuser\git\ESA\esa\saw.py", line 967, in LoadState
        return self._call_simauto('LoadState')
      File "C:\Users\myuser\git\ESA\esa\saw.py", line 1227, in _call_simauto
        raise PowerWorldError(output[0])
    esa.saw.PowerWorldError: LoadState: State hasn't been previously stored.

This behavior is expected - it is not valid to call ``LoadState`` if
``SaveState`` has not yet been called. In the exception above, not that
a ``PowerWorldError`` is raised. This empowers users to handle
exceptions in whatever manner they see fit:

.. code:: python

    >>> from esa import PowerWorldError
    >>> try:
    ...     saw.LoadState()
    ... except PowerWorldError:
    ...     print("Oh my, we've encountered a PowerWorldError!")
    ...
    Oh my, we've encountered a PowerWorldError!


Finally, make sure to clean up after yourself so you don't have COM
objects hanging around.

.. code:: python

    >>> saw.exit()

After walking through this quick start, you should be ready to start
using ESA to improve your simulation and analysis work flows!
This examples shows how to make a histogram of percent line loading in
a power system model using SimAuto, ESA, and Matplotlib.

`Matplotlib <https://matplotlib.org/>`__ is a "comprehensive library for
creating static, animated, and interactive visualizations in Python."
You'll first need to install Matplotlib into your virtual environment
(which should be activated!), which is most easily done by:

.. code:: bat

    python -m pip install -U matplotlib

Before following along with the example, define the CASE_PATH constant
(the file path to a PowerWorld ``.pwb`` case file) like so, adapting the
path to your system.

.. code:: python

  CASE_PATH = r"C:\Users\myuser\git\ESA\tests\cases\tx2000\tx2000_base_pws_version_21.pwb"

Now let's get started!

Perform imports, initialize a ``SAW`` instance:

.. code:: python

    >>> from esa import SAW
    >>> import matplotlib.pyplot as plt

Initialize SAW instance using 2000 bus test case:

.. code:: python

    >>> saw = SAW(FileName=CASE_PATH)

Solve the power flow:

.. code:: python

    >>> saw.SolvePowerFlow()

Let's obtain line loading percentages. But first, we need to know which
fields PowerWorld needs in order to identify branches. These fields are
known as key fields.

.. code:: python

    >>> branch_key_fields = saw.get_key_field_list('Branch')
    >>> branch_key_fields
    ['BusNum', 'BusNum:1', 'LineCircuit']

Get line loading percentage at all buses via SimAuto function:

.. code:: python

    >>> params = branch_key_fields + ['LinePercent']
    >>> branch_data = saw.GetParametersMultipleElement(ObjectType='Branch', ParamList=params)
    >>> branch_data
          BusNum  BusNum:1 LineCircuit  LinePercent
    0       1001      1064           1    30.879348
    1       1001      1064           2    30.879348
    2       1001      1071           1    35.731801
    3       1001      1071           2    35.731801
    4       1002      1007           1     5.342946
    ...      ...       ...         ...          ...
    3199    8157      5124           1    36.371236
    3200    8157      8156           1    46.769588
    3201    8158      8030           1    25.982494
    3202    8159      8158           1    43.641971
    3203    8160      8159           1    57.452701
    <BLANKLINE>
    [3204 rows x 4 columns]

To learn more about variables such as ``LinePercent``, see
:ref:`powerworld-variables`.

Then let's start to plot with Matplotlib!

.. code:: python

    >>> axes = branch_data.plot(kind='hist', y='LinePercent')
    >>> axes.set_xlabel('Line Percent Loading')
    Text(0.5, 0, 'Line Percent Loading')
    >>> axes.set_ylabel('Number of Lines')
    Text(0, 0.5, 'Number of Lines')
    >>> axes.set_title('Histogram of Line Loading')
    Text(0.5, 1.0, 'Histogram of Line Loading')
    >>> plt.show(block=False)

The results should look like:

.. image:: https://github.com/mzy2240/ESA/raw/develop/docs/rst/snippets/line_loading_histogram.png
    :width: 100 %
    :align: center
This example shows how to add transmission lines to a model.

Before starting the example, please define the constants
``CASE_PATH`` (the file path to a PowerWorld ``.pwb`` case file) and
``CANDIDATE_LINES`` (file path to a ``.csv`` file with data related to
lines we'd like to add to the model) like the following, adapting paths
to your system. You can find the case and .csv file referenced in the
``tests`` directory of the ESA repository.

.. code:: python

    CASE_PATH = r"C:\Users\myuser\git\ESA\tests\cases\tx2000\tx2000_base_pws_version_21.pwb"
    CANDIDATE_LINES = r"C:\Users\myuser\git\ESA\tests\data\CandidateLines.csv"

Import packages/classes and read the ``CANDIDATE_LINES`` .csv file.

.. code:: python

    >>> from esa import SAW
    >>> import pandas as pd
    >>> line_df = pd.read_csv(CANDIDATE_LINES)
    >>> line_df
       From Number  To Number  Ckt        R        X        B  Lim MVA A
    0         8155       5358    3  0.00037  0.00750  0.52342       2768
    1         8154       8135    3  0.00895  0.03991  0.00585        149
    2         8153       8108    3  0.01300  0.05400  0.02700        186
    3         8152       8160    3  0.00538  0.03751  0.00613        221
    4         8155       8057    3  0.00037  0.00750  0.52342       2768
    5         8154       8153    3  0.01300  0.05400  0.02700        186
    6         8155       8135    3  0.00538  0.03751  0.00613        221


Instantiate a ``SAW`` object. Set ``CreateIfNotFound`` to ``True`` so
that new lines can be added:

.. code:: python

    >>> saw=SAW(FileName=CASE_PATH, CreateIfNotFound=True, early_bind=True)

Rename columns in the ``line_df`` to match PowerWorld variables. We are
renaming variables from the "Concise Variable Name" convention to the
"Variable Name" convention. See `power_world_object_fields.xlsx
<https://github.com/mzy2240/ESA/blob/master/docs/power_world_object_fields.xlsx>`__.
Also note `this issue
<https://github.com/mzy2240/ESA/issues/1#issue-525219427>`__ is also
relevant. To learn more about PowerWorld variables, see
:ref:`powerworld-variables`.

.. code:: python

    >>> line_df.rename(
    ... columns={
    ... 'From Number': 'BusNum',
    ... 'To Number': 'BusNum:1',
    ... 'Ckt': 'LineCircuit',
    ... 'R': 'LineR',
    ... 'X': 'LineX',
    ... 'B': 'LineC',
    ... 'Lim MVA A': 'LineAMVA'
    ... },
    ... inplace=True)
    >>> line_df.columns
    Index(['BusNum', 'BusNum:1', 'LineCircuit', 'LineR', 'LineX', 'LineC',
           'LineAMVA'],
          dtype='object')

Secondary and tertiary limits are required fields that we must add
manually, since they were not present in the .csv file:

.. code:: python

    >>> line_df['LineAMVA:1'] = 0.0
    >>> line_df['LineAMVA:2'] = 0.0

Check to see if the first line is actually present. An error will
indicate that it's not.

.. code:: python

    >>> line_key_fields = saw.get_key_field_list('branch')
    >>> line_key_fields
    ['BusNum', 'BusNum:1', 'LineCircuit']
    >>> first_line = saw.GetParametersSingleElement('branch', line_key_fields, line_df.loc[0, line_key_fields].tolist())
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "C:\Users\myuser\git\ESA\esa\saw.py", line 693, in GetParametersSingleElement
        output = self._call_simauto('GetParametersSingleElement', ObjectType,
      File "C:\Users\myuser\git\ESA\esa\saw.py", line 1227, in _call_simauto
        raise PowerWorldError(output[0])
    esa.saw.PowerWorldError: GetParameters: Object not found

Enter edit mode to enable the creation of new devices, and use
the ``change_and_confirm_params_multiple_element`` helper function to
easily create the lines. This function will automagically confirm that
the lines will be created.

.. code:: python

    >>> saw.RunScriptCommand("EnterMode(EDIT);")
    >>> saw.change_and_confirm_params_multiple_element('branch', line_df)

Now, we should be able to find that first line without error:

.. code:: python

    >>> first_line = saw.GetParametersSingleElement('branch', line_key_fields, line_df.loc[0, line_key_fields].tolist())
    >>> first_line
    BusNum         8152
    BusNum:1       8160
    LineCircuit       3
    dtype: object

Always clean up:

.. code:: python

    >>> saw.exit()All files within this directory will be used for performing tests, using
the ``doctest`` module. As all ESA tests will require a .pwb file,
please define the full path to that .pwb file as CASE_PATH in the files.
So that the tests know which case to use for CASE_PATH, use a suffix on
the file names.

*   14 --> IEEE 14 bus case
*   2000 --> Texas 2000 bus case (not the modified version)

Also note that all ``2000`` snippets will be passed the
``CANDIDATE_LINES`` variable during testing, as that's needed for
"add_lines_2000.rst".

File ``line_loading_histogram.png`` will be needed for
``line_loading_matplotlib_2000``.This example shows how one can create a weighted graph using branch
impedance values from a PowerWorld grid model as weights. You'll need to
have the `NetworkX <https://networkx.github.io/>`__ Python package
installed into your virtual environment in order to execute this example
on your machine (``python -m pip install networkx``).

Please note that this example does NOT work with Simulator version 17 (
and possibly other versions of Simulator older than version 21). For an
unknown reason, PowerWorld itself throws an exception when trying to run
the ``SaveYbusInMatlabFormat`` script command. If you have a solution
to this problem, please file an `issue
<https://github.com/mzy2240/ESA/issues>`__.

Before following along with the example, define the ``CASE_PATH``
constant (the file path to a PowerWorld ``.pwb`` case file) like so,
adapting the path to your system:

.. code:: python

    CASE_PATH = r"C:\Users\myuser\git\ESA\tests\cases\ieee_14\IEEE 14 bus_pws_version_21.pwb"

Onward!

Imports and initialization:

.. code:: python

    >>> import networkx as nx
    >>> from esa import SAW
    >>> import re
    >>> import os
    >>> saw = SAW(CASE_PATH, early_bind=True)
    >>> g = nx.Graph()

Save YBus matrix to file:

.. code:: python

    >>> ybus_file = CASE_PATH.replace('pwb', 'mat')
    >>> cmd = 'SaveYbusInMatlabFormat("{}", NO)'.format(ybus_file)
    >>> saw.RunScriptCommand(cmd)

Read YBus matrix file into memory. The first two lines are skipped via
the ``readline`` method because they aren't needed.

.. code:: python

    >>> with open(ybus_file, 'r') as f:
    ...     f.readline()
    ...     f.readline()
    ...     mat_str = f.read()
    ...
    'j = sqrt(-1);\n'
    'Ybus = sparse(14);\n'

We're done with the file itself now. Remove it:

.. code:: python

    >>> os.remove(ybus_file)

Remove all white space, split by semicolons, and define a couple regular
expressions (ie --> integer expression, fe --> float expression):

.. code:: python

    >>> mat_str = re.sub(r'\s', '', mat_str)
    >>> lines = re.split(';', mat_str)
    >>> ie = r'[0-9]+'
    >>> fe = r'-*[0-9]+\.[0-9]+'
    >>> exp = re.compile(r'(?:Ybus\()({ie}),({ie})(?:\)=)({fe})(?:\+j\*)(?:\()({fe})'.format(ie=ie, fe=fe))

Loop over the lines from the file and build up the graph. Ignore
diagonal Y bus matrix entries and buses which are not connected
(have 0 admittance between them).

.. code:: python

    >>> for line in lines:
    ...     match = exp.match(line)
    ...     if match is None:
    ...         continue
    ...     idx1, idx2, real, imag = match.groups()
    ...     if idx1 == idx2:
    ...         continue
    ...     neg_admittance = float(real) + 1j * float(imag)
    ...     try:
    ...         impedance = -1 / neg_admittance
    ...     except ZeroDivisionError:
    ...         continue
    ...     g.add_edge(int(idx1), int(idx2), r=impedance.real, x=impedance.imag)
    ...

Explore some graph properties to ensure it worked:

.. code:: python

    >>> g.number_of_nodes()
    14
    >>> g.number_of_edges()
    20
    >>> data_1_2 = g.get_edge_data(1, 2)
    >>> data_1_2['r']
    0.01937987032338931
    >>> data_1_2['x']
    0.05917003035204804

As always, clean up when done:

.. code:: python

    >>> saw.exit()
This example illustrates the procedure to perform transient stability(TS) analysis
and obtain the TS result. To retrieve the result, the most convenient way is to create
a plot object and connect the object/fields pairs to it, then you will be able to query it
using the :code:`TSGetCongencyResults` function.

.. code:: python

    CASE_PATH = r"C:\Users\myuser\git\ESA\tests\cases\il200\ACTIVSg200.pwb"

Load the case first and then solve a PF (optional):

.. code:: python

    >>> from esa import SAW
    >>> saw = SAW(CASE_PATH)
    >>> saw.SolvePowerFlow()

Then perform TS analysis (make sure you already have a desired plot object)

.. code:: python

    >>> t1 = 0.0
    >>> t2 = 15.0
    >>> stepsize = 0.01

        # Solve.
    >>> cmd = 'TSSolve("{}",[{},{},{},NO])'.format(
            self.ctg_name, t1, t2, stepsize
        )
    >>> saw.RunScriptCommand(cmd)

Once it is done, you could retrieve (and visualize) the results:

.. code:: python

    >>> objFieldList = ['Plot ''Area_Avg Bus Hz''']  # "Area_Avg Bus Hz" is the plot name
    >>> result = sa.TSGetContingencyResults("My Transient Contingency", objFieldList, 0, 12)  # "My Transient Contingency" is the contingency name
    >>> df = result[1]  #result[0] is meta data
    >>> df.columns = ['Time (s)', 'Area_Avg Bus Hz']
    >>> df.plot(x='Time (s)', y='Area_Avg Bus Hz')
.. image:: https://github.com/mzy2240/ESA/raw/develop/docs/rst/snippets/ts_result.png
    :width: 100 %
    :align: center

The whole process, including setting up plots and creating contingencies, could be fully
automated, but it might be easier for most users to pre-define the plots and contingencies
in the case and then load the case using ESA. GetParametersMultipleElement cannot be used
here to retrieve the TS datapoints (which is a very rare situation).

Simply put, ESA makes working with SimAuto a breeze! ESA has the
following desirable properties:

*   Free, open-source, non-restrictive license
*   Fully documented, including a plethora of examples and common issues
*   Fully tested - rest assured that all functions work as intended
*   Abstracts away all interactions with Windows COM
*   Automatically parses errors from PowerWorld
*   Automatically transforms data into the correct types
*   Wraps all available SimAuto functions
*   Provides additional helper functions to further simplify
    interactions with SimAuto
*   Provides functions to interact with PowerWorld Simulator interface,
    without the hassle of creating and loading the display auxiliary
    files
*   Returns useful data types such as Pandas DataFrames, unlocking a
    whole new realm of analysis and control capabilities
*   Gets data from a power system model for other application, like
    graph analysis or signal processing
*   Co-simulation with other domains such as transportation,
    communication, or natural gas
*   Task automation, like dynamic model tuning and relay setting
    tuning
*   AI & ML friendly
*   Compatible with several modern Python versions
*   Lightweight and fast

Whether you're an energy trader, transmission planner, or academic
researcher, ESA will help simplify your simulation an analysis work
flows.ESA was developed at Texas A&M University by `Brandon Thayer
<https://github.com/blthayer>`__, `Zeyu Mao
<https://github.com/mzy2240>`__, and `Yijing Liu
<https://github.com/SmartJingJing123>`__. Significant guidance and
oversight was provided by `Professor Thomas Overbye
<https://engineering.tamu.edu/electrical/profiles/overbye-thomas.html>`__,
who is a co-founder of PowerWorld Corporation, and `Professor Katherine
Davis <https://engineering.tamu.edu/electrical/profiles/davis-katherine
.html>`__, who is one of the original authors of SimAuto.
Funding was provided by the Texas A&M Engineering Experiment Station's
`Smart Grid Center <https://smartgridcenter.tamu.edu/>`__ and the U.S.
Department of Energy (DOE) under award DE-OE0000895.

The authors of ESA would like to give a special thank you to our users.
Additionally, we've received help from some of our colleagues along the
way, including (but not limited to!) Wei Trinh and Diana Wallison.
Thank you!

We hope you enjoy using ESA as much as we've enjoyed creating it.ESA is a "Pip-installable" Python package that eases interaction with
the PowerWorld Simulator Automation Server (SimAuto). PowerWorld
Simulator is a powerful, commercial-grade electric grid simulation tool
with a wide range of capabilities. Information on Simulator can be found
`here
<https://www.powerworld.com/products/simulator/overview>`__ and
information on SimAuto can be found `here
<https://www.powerworld.com/products/simulator/add-ons-2/simauto>`__.
Since ESA directly interfaces with SimAuto, ESA users will need a
PowerWorld license and installation that also includes SimAuto.

ESA makes working with SimAuto, well, easy. Users don't have to worry
about input or output data type conversions, data mapping,
determining whether SimAuto has reported an error, and more.
Additionally, ESA uses the scientific computing packages you know and
love, including Numpy and Pandas. In addition to wrapping 100% of the
functions provided by SimAuto, ESA provides helper functions that
further ease development. Below is a quick motivating example (also
found in :ref:`increase-loading`) that shows how easy it is to use
SimAuto.

.. include:: snippets/increase_loading_14.rst
ESA is useful for a wide range of audiences, including:

*   Industry practitioners (e.g. power system planners, energy traders, etc.)
*   Power system researchers
*   Researchers from other domains who wish to perform co-simulation
    with power system models
*   University students and faculty
*   Power system enthusiasts

ESA users should be proficient in Python, and it's recommended that
users get familiar with Numpy and Pandas, as ESA makes significant use
of these packages. ESA users do not need to have any knowledge
whatsoever related to how the Windows COM API works, nor do users need
to be familiar with PyWin32.

Ultimately, ESA is a tool for interacting with PowerWorld Simulator -
thus, users should have some familiarity with Simulator. Users do not
need to directly understand how to use SimAuto, as ESA abstracts those
details away. Advanced users will have a solid understanding of
PowerWorld variables and object types, and will make extensive use of
the ``RunScriptCommand`` method to enable the execution of PowerWorld
functions previously only accessible via `"Auxiliary Files"
<https://github.com/mzy2240/ESA/blob/master/docs/Auxiliary%20File%20Format.pdf>`__.ESA has already been utilized in several research projects past and
present. If you use ESA in your work, please file an issue on
`GitHub <https://github.com/mzy2240/ESA/issues>`__ and we'll list your
project here! Please cite ESA if you use it in your work: see
:ref:`citation`.

-   In `gym-powerworld <https://github.com/blthayer/gym-powerworld>`__,
    ESA was used to create a standardized reinforcement learning
    environment for power system voltage control. This environment was
    then `used to carry out deep reinforcement learning (DRL)
    experiments <https://github.com/blthayer/drl-powerworld>`__
    in which the algorithm attempts to learn how to best control grid
    voltages under a diverse set of grid conditions.
-   In `this paper
    <https://ieeexplore.ieee.org/abstract/document/9042493>`__,
    ESA was leveraged to create and simulate different electric grid
    scenarios where load, renewable generation levels, generation
    capacities, scheduled outages, and unit commitment were all varied.
    The resulting scenarios were used in the
    `Grid Optimization (GO) competition
    <https://gocompetition.energy.gov/>`__
    hosted by the U.S. Department of Energy (DOE).
-   Geomagnetic disturbances (GMDs) affect the magnetic and electric field
    of the earth, inducing dc voltage sources superimposed on transmission
    lines. In an accepted paper by Martinez et al. entitled
    "Undergraduate Research on Design Considerations for a GMD
    Mitigation Systems" (to be published in mid-late 2020), a
    planning-based GMD mitigation strategy was developed for large power
    systems. ESA is leveraged to programmatically place GIC blocking
    devices in test systems per the proposed algorithm, thus minimizing
    the effects of GMDs on the power grid.
-   ESA is used by an ongoing research project entitled "Real Time
    Monitoring Applications for the Power Grid under Geomagnetic
    Disturbances (GMD)": Recently, a real-world GMD monitoring system
    consisting of six magnetometers was deployed in Texas. The resulting
    magnetic field measurements are coupled with ground conductivity models
    to calculate real-time electric fields. These can then be fed to a grid
    model of Texas using ESA to enable calculation of real-time
    geomagnetically induced currents (GICs) for monitoring and
    visualization.
-   ESA is used by an ongoing research project entitled "Contingency
    Analysis Based on Graph Theory Concepts and Line Outage Distribution
    Factors (LODF)." In this project, ESA is leveraged to extract the
    topology of the power system model and obtain the LODF matrix.From PowerWorld's `SimAuto page
<https://www.powerworld.com/products/simulator/add-ons-2/simauto>`__::

    "The Simulator Automation Server (SimAuto) allows you to take
    advantage of the power of automation to extend the functionality of
    PowerWorld Simulator to any external program that you write. Using
    Simulator Automation Server you can launch and control PowerWorld
    Simulator from within another application, enabling you to: access
    the data of a Simulator case, perform defined Simulator functions
    and other data manipulations, and send results back to your original
    application, to a Simulator auxiliary file, or to a Microsoft® Excel
    spreadsheet."

In essence, SimAuto is PowerWorld Simulator's application programming
interface (API). As such, SimAuto users can perform *almost* any task
that can be done through Simulator's graphic user interface (GUI), but
via their own program. This opens up a wealth of opportunity to perform
tasks such as:

*   Task automation

*   Sensitivity analysis

*   Co-simulation

*   Machine learning

*   And more!

For more SimAuto details, here are some PowerWorld links:

*   `SimAuto Description <https://www.powerworld.com/products/simulator/add-ons-2/simauto>`__

*   `PowerWorld Web Help <https://www.powerworld.com/WebHelp/>`__

*   `SimAuto Documentation`_

Since SimAuto strives to be accessible to "any external" program and
uses Windows `COM
<https://docs.microsoft.com/en-us/windows/win32/com/the-component-object-model>`__,
it can be cumbersome, tedious, and difficult to use. That's where ESA
comes in!

SimAuto Functions
^^^^^^^^^^^^^^^^^

Here's a listing of the currently (as of 2020-06-17, Simulator version
21) available SimAuto functions (documented `here <simauto-docs_>`_):

*   ChangeParameters
*   ChangeParametersSingleElement
*   ChangeParametersMultipleElement
*   CloseCase
*   GetFieldList
*   GetParametersSingleElement
*   GetParametersMultipleElement
*   GetParameters
*   GetSpecificFieldList
*   GetSpecificFieldMaxNum
*   ListOfDevices
*   LoadState
*   OpenCase
*   ProcessAuxFile
*   :ref:`runscriptcommand`
*   SaveCase
*   SaveState
*   SendToExcel (not recommended for use with ESA as documented in
    :ref:`esa-saw-api`)
*   TSGetContingencyResults
*   WriteAuxFile

For ESA's implementation/wrapping of these methods, see
:ref:`esa-saw-api`.

SimAuto Properties
^^^^^^^^^^^^^^^^^^

Here's a listing of the currently (as of 2020-06-17, Simulator version
21) available SimAuto properties (documented `here <simauto-docs_>`_):

*   ExcelApp (like ``SendToExcel`` function, not recommended for use
    with ESA)
*   CreateIfNotFound
*   CurrentDir
*   ProcessID
*   RequestBuildDate
*   UIVisible (Simulator versions >= 20)

For ESA's implementation/wrapping of these properties, see
:ref:`esa-saw-api`.

.. _SimAuto Documentation: https://www.powerworld.com/WebHelp/#MainDocumentation_HTML/Simulator_Automation_Server.htm%3FTocPath%3DAutomation%2520Server%2520Add-On%2520(SimAuto)%7C_____1
.. _simauto-docs: `SimAuto Documentation`_As noted in `ESA's license
<https://github.com/mzy2240/ESA/blob/master/LICENSE>`__, no warranty
is provided, and the authors cannot be held liable for any issues
related to using ESA. If you encounter an issue, find a bug, or would
like to provide feedback, please file a ticket on `Github
<https://github.com/mzy2240/ESA/issues>`__.Changes made with each ESA release are listed here. Please note that
versions prior to 1.0.0 are not listed here, but are still available on
`PyPi <https://pypi.org/project/esa/#history>`__.

Version 1.2.3
^^^^^^^^^^^^^

* Fix the AOT import error

Version 1.2.2
^^^^^^^^^^^^^

* Fix the AOT version-dependent issue
* Update the dependency version

Version 1.2.1
^^^^^^^^^^^^^

* Greatly improve the fast contingency analysis by taking advantage of
  SIMD, JIT and AOT. Now it could finish a N-1 and N-2 contingency analysis of
  a synthetic 2000 bus grid in less than 15 seconds!
* Adjust the release process to include AOT functions

Version 1.2.0
^^^^^^^^^^^^^

* Optimize the process to use the same order as shown in simulator
  (note: if pw_order is used, all data in the dataframe will be string type)
* Add a method to obtain the incidence matrix
* Implement a modified fast N-1 and N-2 contingency analysis algorithm.
  The algorithm is originally developed by Prof. Kostya Turitsyn from MIT and
  the implementation has been slightly modified and adapted to work with ESA.
* Add a few helper functions to facilitate contingency analysis powered by the simulator.

Version 1.1.0
^^^^^^^^^^^^^

* Allow users to use the same order as shown in simulator for all the
  dataframes
* Add a helper function to generate LODF matrix

Version 1.0.9
^^^^^^^^^^^^^

* Update the pre-install process and the common issues
* Update the helper function 'get_ybus' with a new argument to accept
  external ybus file

Version 1.0.8
^^^^^^^^^^^^^

* Add new helper function 'to_graph'. The new function could help
  generate NetworkX graph model from the case, in two different levels:
  bus-as-node and substation-as-node. Parallel lines are preserved, and
  directedgraph is supported (currently the direction is fixed to be
  the same as real power flow).

Version 1.0.7
^^^^^^^^^^^^^

* Add new functions: get_ybus, get_jacobian

Version 1.0.6
^^^^^^^^^^^^^

* Hopefully finally fixing locale-based issues. Fixes began in 1.0.4,
  and continued in 1.0.5.
* Finalizing JOSS paper. Once published, the citation will be added to
  the top-level README and the documentation.

Version 1.0.5
^^^^^^^^^^^^^

* It turns out version 1.0.4 did not fully/correctly handle automatic
  locale setting. This version should now properly handle different
  decimal delimiters automatically.
* Bug fix: The ``additional_fields`` parameter to ``SAW``'s
  ``get_power_flow_results`` was permanently adding the
  ``additional_fields`` to the corresponding list in the ``SAW``
  object's ``SAW.POWER_FLOW_FIELDS`` attribute.

Version 1.0.4
^^^^^^^^^^^^^

* Added support for other locales by automatically detecting the
  system's decimal delimiter. This should allow users in Europe and
  elsewhere to leverage ESA. Thanks to
  `robinroche <https://github.com/robinroche>`__ for pointing out the
  problem during our `JOSS <https://joss.theoj.org/>`__ review in
  `this comment <https://github.com/openjournals/joss-reviews/issues/2289#issuecomment-643482550>`__.

Version 1.0.3
^^^^^^^^^^^^^

* New SAW attribute, ``build_date``
* New SAW attribute, ``version``
* New SAW helper method, ``get_version_and_builddate``
* Add argument ``additional_fields`` for ``get_power_flow_results`` method
  which provides an easy and consistent way to add more fields to the power
  flow result
* Updating so that ESA is compatible with Simulator version 17. Note
  that this does not imply ESA has been tested with versions 16, 18, 19,
  or 20. However, ESA *should* work with all these versions.
* Added case files for Simulator versions 16-22(beta) and renamed the cases
  accordingly (suffixed with ``pws_version_<version goes here>.pwb``.
* Updated documentation to discuss different versions of Simulator.

Version 1.0.2
^^^^^^^^^^^^^

* Add area number to the power flow result
* Update the citation section
* Fix a bug in the test file that will result in a failure if some
  default names are changed in PowerWorld

Version 1.0.1
^^^^^^^^^^^^^

* Add new functions: update_ui, OpenOneline and CloseOneline
* Add documents to meet the requirement of JOSS
* Add one more example into the documentation
* Update the coverage_to_rst.py so that it's more clear that the errors
  that get printed during testing are as expected.
* Update the release process
* Fix minor typos

Version 1.0.0
^^^^^^^^^^^^^

ESA version 1.0.0 is the first ESA release in which 100% of SimAuto
functions are wrapped, and testing coverage is at 100%.
