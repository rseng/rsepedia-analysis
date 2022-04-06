![logo](https://raw.githubusercontent.com/pybamm-team/liionpack/main/docs/liionpack.png)

#
<div align="center">

[![liionpack](https://github.com/pybamm-team/liionpack/actions/workflows/test_on_push.yml/badge.svg)](https://github.com/pybamm-team/liionpack/actions/workflows/test_on_push.yml)
[![Documentation Status](https://readthedocs.org/projects/liionpack/badge/?version=main)](https://liionpack.readthedocs.io/en/main/?badge=main)
[![codecov](https://codecov.io/gh/pybamm-team/liionpack/branch/main/graph/badge.svg)](https://codecov.io/gh/pybamm-team/liionpack)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pybamm-team/liionpack/blob/main/)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04051/status.svg)](https://doi.org/10.21105/joss.04051)

</div>

# Overview of liionpack
*liionpack* takes a 1D PyBaMM model and makes it into a pack. You can either specify
the configuration e.g. 16 cells in parallel and 2 in series (16p2s) or load a
netlist.

## Installation

Follow the steps given below to install `liionpack`. The package must be installed to run the included examples. It is recommended to create a virtual environment for the installation, see [the documentation](https://liionpack.readthedocs.io/en/main/install/).

To install `liionpack` using `pip`, run the following command:
```bash
pip install liionpack
```

### LaTeX

In order to use the `draw_circuit` functionality a version of Latex must be installed on your machine. We use an underlying Python package `Lcapy` for making the drawing and direct you to its installation instructions [here](https://lcapy.readthedocs.io/en/latest/install.html) for operating system specifics.

## Example Usage

The following code block illustrates how to use liionpack to perform a simulation:

```python
import liionpack as lp
import numpy as np
import pybamm

# Generate the netlist
netlist = lp.setup_circuit(Np=16, Ns=2, Rb=1e-4, Rc=1e-2, Ri=5e-2, V=3.2, I=80.0)

output_variables = [
    'X-averaged total heating [W.m-3]',
    'Volume-averaged cell temperature [K]',
    'X-averaged negative particle surface concentration [mol.m-3]',
    'X-averaged positive particle surface concentration [mol.m-3]',
]

# Heat transfer coefficients
htc = np.ones(32) * 10

# Cycling experiment, using PyBaMM
experiment = pybamm.Experiment([
    "Charge at 20 A for 30 minutes",
    "Rest for 15 minutes",
    "Discharge at 20 A for 30 minutes",
    "Rest for 30 minutes"],
    period="10 seconds")

# PyBaMM parameters
chemistry = pybamm.parameter_sets.Chen2020
parameter_values = pybamm.ParameterValues(chemistry=chemistry)

# Solve pack
output = lp.solve(netlist=netlist,
                  parameter_values=parameter_values,
                  experiment=experiment,
                  output_variables=output_variables,
                  htc=htc)
```

## Documentation

There is a full API documentation, hosted on Read The Docs that can be found [here](https://liionpack.readthedocs.io/).

## Contributing to liionpack

If you'd like to help us develop liionpack by adding new methods, writing documentation, or fixing embarrassing bugs, please have a look at these [guidelines](https://github.com/pybamm-team/liionpack/blob/main/docs/contributing.md) first.

## Get in touch

For any questions, comments, suggestions or bug reports, please see the [contact page](https://www.pybamm.org/contact).

## Acknowledgments

PyBaMM-team acknowledges the funding and support of the Faraday Institution's multi-scale modelling project and Innovate UK.

The development work carried out by members at Oak Ridge National Laboratory was partially sponsored by the Office of Electricity under the United States Department of Energy (DOE).

## License

liionpack is fully open source. For more information about its license, see [LICENSE](https://github.com/pybamm-team/liionpack/blob/main/LICENSE).
---
title: 'liionpack: A Python package for simulating packs of batteries with PyBaMM'

tags:
  - Python
  - batteries
  - packs
  - electrochemistry

authors:
  - name: Thomas G. Tranter
    orcid: 0000-0003-4721-5941
    affiliation: "1, 2"
  - name: Robert Timms
    orcid: 0000-0002-8858-4818
    affiliation: "2, 3"
  - name: Valentin Sulzer
    orcid: 0000-0002-8687-327X
    affiliation: "4"
  - name: Ferran Brosa Planella
    orcid: 0000-0001-6363-2812
    affiliation: "2, 5"
  - name: Gavin M. Wiggins
    orcid: 0000-0002-4737-6596
    affiliation: "6"
  - name: Suryanarayana V. Karra
    orcid: 0000-0002-5671-0998
    affiliation: "6"
  - name: Priyanshu Agarwal
    orcid: 0000-0002-5333-1634
    affiliation: "7"
  - name: Saransh Chopra
    orcid: 0000-0003-3046-7675
    affiliation: "8"
  - name: Srikanth Allu
    orcid: 0000-0003-2841-4398
    affiliation: "6"
  - name: Paul R. Shearing
    orcid: 0000-0002-1387-9531
    affiliation: "1, 2"
  - name: Dan J. L. Brett
    orcid: 0000-0002-8545-3126
    affiliation: "1, 2"

affiliations:
 - name: Department of Chemical Engineering, University College London, London, WC1E 7JE, United Kingdom.
   index: 1
 - name: The Faraday Institution, Quad One, Becquerel Avenue, Harwell Campus, Didcot, OX11 0RA, United Kingdom.
   index: 2
 - name: Mathematical Institute, University of Oxford, OX2 6GG, United Kingdom.
   index: 3
 - name: Carnegie Mellon University, Scott Hall 5109, 5000 Forbes Ave, Pittsburgh, PA 15213, United States.
   index: 4
 - name: WMG, University of Warwick, Coventry, CV4 7AL, United Kingdom
   index: 5
 - name: Oak Ridge National Laboratory, 2360 Cherahala Boulevard, Knoxville, Tennessee 37932, United States.
   index: 6
 - name: Symbiosis Institute of Technology, Symbiosis International University, Lavale, Pune, Maharashtra 412115, India.
   index: 7
 - name: Cluster Innovation Centre, University of Delhi, GC Narang Road, Delhi, 110007, India.
   index: 8

date: 03 December 2021

bibliography: paper.bib

---

# Summary

Electrification of transport and other energy intensive activities is of growing importance as it provides an underpinning method to reduce carbon emissions. With an increase in reliance on renewable sources of energy and a reduction in the use of more predictable fossil fuels in both stationary and mobile applications, energy storage will play a pivotal role and batteries are currently the most widely adopted and versatile form. Therefore, understanding how batteries work, how they degrade, and how to optimize and manage their operation at large scales is critical to achieving emission reduction targets. The electric vehicle (EV) industry requires a considerable number of batteries even for a single vehicle, sometimes numbering in the thousands if smaller cells are used, and the dynamics and degradation of these systems, as well as large stationary power systems, is not that well understood. As increases in the efficiency of a single battery become diminishing for standard commercially available chemistries, gains made at the system level become more important and can potentially be realised more quickly compared with developing new chemistries. Mathematical models and simulations provide a way to address these challenging questions and can aid the engineer and designers of batteries and battery management systems to provide longer lasting and more efficient energy storage systems.

# Statement of need

`liionpack` is a PyBaMM-affiliated Python package for simulating large systems of batteries connected in series and parallel. Python enables wrapping low-level languages (e.g., C) for speed without losing flexibility or ease-of-use in the user-interface. `liionpack` was designed to be used by physicists, engineers, students, academics and industrial researchers and system designers concerned with the dynamics of electric current and heat transport in large battery systems. Commercial battery pack simulation tools are available such as modules that can be included within Comsol&reg;, Simulink&reg; and STAR-CCM+&trade;, but to our knowledge `liionpack` is the first to be released open-source. The commercial packages contain more advanced features such as GUI's for circuit design, and integration with CAD based thermal and fluid dynamics tools, but `liionpack` provides everything you need to model a pack of batteries with simple physics and can incorporate circuit definitions defined elsewhere and heat transfer coefficients that are calculated elsewhere. We hope that it will provide the battery community with a platform to build upon to add more features in the future and increase productivity, reproducibility and transparency in this research space.

The API for `liionpack` was designed to provide a simple and efficient extension to the `PyBaMM` [@pybamm] framework allowing users to scale up simulations from single cells to many thousands with a few extra lines of code. `PyBaMM` provides a number of classic physics-based single battery models with configurable options to investigate thermal effects and degradation, for example. The pack architecture introduced by `liionpack` can be defined as a number of batteries connected in series and parallel to one another using busbars and interconnections with defined resistances. A netlist may also be used to construct the pack which is more flexible and allows for configurable network topology and can be constructed graphically with packages such as `LTSpice` [@ltspice] or simply created manually, specifying nodal connections as either current sources, voltage sources or resistors. Statistical distributions can be easily incorporated into the pack architecture elements through the use of input parameters that allow a single model to be solved with varying inputs.

![Coupled system solution algorithm.\label{fig:0}](./paper_figures/Figure_0.png)

# Algorithm

The algorithm to solve the coupled system of batteries is shown in \autoref{fig:0}. The nature of the solving process facilitates parallel processing of the electrochemical problem for each battery during each time-step formulated as an integrable set of 1D differential-algebraic equations (DAEs). The system is coupled electrically at the global level via the busbars and interconnections in the circuit and solving this linear algebraic system between electrochemical time-steps determines the current balance and boundary conditions for each battery at the next time-step. The combination of a global circuit solve and local electrochemical solve repeatedly iterated over in time in a see-saw fashion provides the most simple and efficient way of coupling the system without repeating time-steps. Results for solving a single battery forming a circuit with negligible busbar resistance deviates by less than 0.01% from a pure `PyBaMM` simulation.

At present, the circuits that are solved may only contain three different types of element: namely current sources, voltage sources, and resistors. Resistors are used to represent the busbars and interconnections in the pack as well as the internal resistance of the batteries. The open circuit voltage is used for the voltage sources in the circuit and modified nodal analysis (MNA) [@mna] is used to solve the circuit problem determining the distribution of current in the pack. A typical 4p1s pack architecture is shown below in \autoref{fig:1}, which was produced using `Lcapy` [@lcapy].

![Typical pack architecture.\label{fig:1}](./paper_figures/Figure_1.png)

Presently, the thermal problem is solved in a non-coupled way with each battery acting as an independent heat source and interacting with its environment in a "lumped" sense with a volume-averaged heat transfer coefficient. Heat generation and conduction through the busbars and from cell to neighbouring cells is likely to occur in some scenarios and can be accounted for by solving a transient thermal problem on the network architecture [@jellyroll], which will be implemented in future releases. Heat transfer coefficients may also be easily adjusted on a cell-by-cell basis and also throughout the simulation solving process to reflect heterogenous and time-dependent cooling conditions.

Several distributed solvers are provided and can be selected through a common function with a simple function argument. These are `Casadi` [@casadi], which uses multi-threading and works well for single workstations, and `ray` [@ray] and `dask` [@dask], which are designed for running on clusters and use multi-processing. Many of the functions and models that can be found in `PyBaMM` should work in exactly the same way in `liionpack` and examples are provided showing how to set up and configure different battery models for running in the pack system. Several visualization tools are also provided for analysis of the results.

# Example

An example of a small pack is included below. A 4p1s configuration is defined with busbar resistance of 1 $m\Omega$ and interconnection resistance of 10 $m\Omega$. The `Chen2020` [@Chen2020] parameter set is used to define the battery cell chemistry which was gathered using an LG M50 cylindrical cell of 21700 format. By default the single particle model `SPM` is used to define the electrochemical battery model system but a suite of others are available [@Marquis2020] and can be configured using a custom simulation.

```
import liionpack as lp
import pybamm

# Generate the netlist
netlist = lp.setup_circuit(Np=4, Ns=1, Rb=1e-3, Rc=1e-2)

# Define some additional variables to output
output_variables = [
    'X-averaged negative particle surface concentration [mol.m-3]',
    'X-averaged positive particle surface concentration [mol.m-3]',
]

# Cycling experiment, using PyBaMM
experiment = pybamm.Experiment([
    "Charge at 5 A for 30 minutes",
    "Rest for 15 minutes",
    "Discharge at 5 A for 30 minutes",
    "Rest for 30 minutes"],
    period="10 seconds")

# PyBaMM battery parameters
chemistry = pybamm.parameter_sets.Chen2020
parameter_values = pybamm.ParameterValues(chemistry=chemistry)

# Solve the pack problem
output = lp.solve(netlist=netlist,
                  parameter_values=parameter_values,
                  experiment=experiment,
                  output_variables=output_variables,
                  initial_soc=0.5)

# Display the results
lp.plot_output(output)

# Draw the circuit at final state
lp.draw_circuit(netlist, cpt_size=1.0, dpi=150, node_spacing=2.5)
```

The output for the examples is shown below as a pack summary in \autoref{fig:2} and an example of a cell variable plot showing each battery current in \autoref{fig:3}.


![Pack summary showing the pack terminal voltage and total current. \label{fig:2}](./paper_figures/Figure_2.png)

![An example of individual cell variable data, any variable defined by the `PyBaMM` model should be accessible. \label{fig:3}](./paper_figures/Figure_3.png)

# Acknowledgements

PyBaMM-team acknowledges the funding and support of the Faraday Institution's multi-scale modelling project under grant number EP/S003053/1, FIRG025.

The development work carried out by members at Oak Ridge National Laboratory was partially sponsored by the Office of Electricity under the United States Department of Energy (DOE).

# References
![logo](liionpack.png)

# Welcome to liionpack

The lithium-ion battery pack simulator powered by [PyBaMM](https://www.pybamm.org/). Liionpack allows you to specify pack configurations with numbers of cells connected in series and parallel or by uploading a netlist.

Leverage the experiments and parameter sets from PyBaMM and scale up your simulations to pack level. Include thermal effects and account for pack position dependency with variable inputs for heat transfer.

Include statistical distributions in the battery parameters.
# Contributing to liionpack

If you'd like to contribute to liionpack (thanks!), please have a look at the [guidelines below](#workflow).

If you're already familiar with our workflow, maybe have a quick look at the [pre-commit checks](#pre-commit-checks) directly below.

## Pre-commit checks

Fork the repository and create a pull request. Github actions should check that tests are passing.

## Workflow

We use [GIT](https://en.wikipedia.org/wiki/Git) and [GitHub](https://en.wikipedia.org/wiki/GitHub) to coordinate our work. When making any kind of update, we try to follow the procedure below.

### A. Before you begin

1. Create an [issue](https://guides.github.com/features/issues/) where new proposals can be discussed before any coding is done.
2. Create a [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/) of this repo (ideally on your own [fork](https://help.github.com/articles/fork-a-repo/)), where all changes will be made
3. Download the source code onto your local system, by [cloning](https://help.github.com/articles/cloning-a-repository/) the repository (or your fork of the repository).
4. [Install](https://pybamm.readthedocs.io/en/latest/install/install-from-source.html) PyBaMM with the developer options.
5. [Test](#testing) if your installation worked, using the test script: `$ python -m unittest`.

You now have everything you need to start making changes!

### B. Writing your code

5. liionpack is developed in [Python](https://en.wikipedia.org/wiki/Python_(programming_language)), and makes heavy use of [NumPy](https://en.wikipedia.org/wiki/NumPy) (see also [NumPy for MatLab users](https://numpy.org/doc/stable/user/numpy-for-matlab-users.html) and [Python for R users](http://blog.hackerearth.com/how-can-r-users-learn-python-for-data-science)).
6. Make sure to follow our [coding style guidelines](#coding-style-guidelines).
7. Commit your changes to your branch with [useful, descriptive commit messages](https://chris.beams.io/posts/git-commit/): Remember these are publicly visible and should still make sense a few months ahead in time. While developing, you can keep using the GitHub issue you're working on as a place for discussion. [Refer to your commits](https://stackoverflow.com/questions/8910271/how-can-i-reference-a-commit-in-an-issue-comment-on-github) when discussing specific lines of code.
8. If you want to add a dependency on another library, or re-use code you found somewhere else, have a look at [these guidelines](#dependencies-and-reusing-code).

### C. Merging your changes with liionpack

9. [Test your code!](#testing)
10. liionpack has online documentation at http://liionpack.readthedocs.io/. To make sure any new methods or classes you added show up there, please read the [documentation](#documentation) section.
11. If you added a major new feature, perhaps it should be showcased in an [example notebook](#example-notebooks).
12. When you feel your code is finished, or at least warrants serious discussion, run the [pre-commit checks](#pre-commit-checks) and then create a [pull request](https://help.github.com/articles/about-pull-requests/) (PR) on [liionpack's GitHub page](https://github.com/pybamm-team/liionpack).
13. Once a PR has been created, it will be reviewed by any member of the community. Changes might be suggested which you can make by simply adding new commits to the branch. When everything's finished, someone with the right GitHub permissions will merge your changes into liionpack main repository.



## Coding style guidelines

liionpack follows the [PEP8 recommendations](https://www.python.org/dev/peps/pep-0008/) for coding style. These are very common guidelines, and community tools have been developed to check how well projects implement them.

### Flake8

We use [flake8](http://flake8.pycqa.org/en/latest/) to check our PEP8 adherence. To try this on your system, navigate to the liionpack directory in a console and type

```bash
flake8
```


### Black

We use [black](https://black.readthedocs.io/en/stable/) to automatically configure our code to adhere to PEP8. Black can be used in two ways:

1. Command line: navigate to the liionpack directory in a console and type

```bash
black {source_file_or_directory}
```

2. Editor: black can be [configured](https://test-black.readthedocs.io/en/latest/editor_integration.html) to automatically reformat a python script each time the script is saved in an editor.

If you want to use black in your editor, you may need to change the max line length in your editor settings.

Even when code has been formatted by black, you should still make sure that it adheres to the PEP8 standard set by [Flake8](#flake8).

### Naming

Naming is hard. In general, we aim for descriptive class, method, and argument names. Avoid abbreviations when possible without making names overly long, so `mean` is better than `mu`, but a class name like `MyClass` is fine.

Class names are CamelCase, and start with an upper case letter, for example `MyOtherClass`. Method and variable names are lower case, and use underscores for word separation, for example `x` or `iteration_count`.


## Dependencies and reusing code

While it's a bad idea for developers to "reinvent the wheel", it's important for users to get a _reasonably sized download and an easy install_. In addition, external libraries can sometimes cease to be supported, and when they contain bugs it might take a while before fixes become available as automatic downloads to liionpack users.
For these reasons, all dependencies in liionpack should be thought about carefully, and discussed on GitHub.

Direct inclusion of code from other packages is possible, as long as their license permits it and is compatible with ours, but again should be considered carefully and discussed in the group. Snippets from blogs and [stackoverflow](https://stackoverflow.com/) can often be included without attribution, but if they solve a particularly nasty problem (or are very hard to read) it's often a good idea to attribute (and document) them, by making a comment with a link in the source code.


## Testing

All code requires testing. We use the [unittest](https://docs.python.org/3.3/library/unittest.html) package for our tests. (These tests typically just check that the code runs without error, and so, are more _debugging_ than _testing_ in a strict sense. Nevertheless, they are very useful to have!)

```bash
python -m unittest
```

### Writing tests

Every new feature should have its own test. To create ones, have a look at the `test` directory and see if there's a test for a similar method. Copy-pasting this is a good way to start.

Next, add some simple (and speedy!) tests of your main features. If these run without exceptions that's a good start! Next, check the output of your methods using any of these [assert methods](https://docs.python.org/3.3/library/unittest.html#assert-methods).


### Profiling

Sometimes, a bit of code will take much longer than you expect to run. In this case, you can set
```python
from IPython import embed; embed(); import ipdb; ipdb.set_trace()
```
as above, and then use some of the profiling tools. In order of increasing detail:
1. Simple timer. In ipython, the command
```
%time command_to_time()
```
tells you how long the line `command_to_time()` takes. You can use `%timeit` instead to run the command several times and obtain more accurate timings.
2. Simple profiler. Using `%prun` instead of `%time` will give a brief profiling report
3. Detailed profiler. You can install the detailed profiler `snakeviz` through pip:
```bash
pip install snakeviz
```
and then, in ipython, run
```
%load_ext snakeviz
%snakeviz command_to_time()
```
This will open a window in your browser with detailed profiling information.

## Documentation

liionpack is documented in several ways.

First and foremost, every method and every class should have a [docstring](https://www.python.org/dev/peps/pep-0257/) that describes in plain terms what it does, and what the expected input and output is.

These docstrings can be fairly simple, but can also make use of [reStructuredText](http://docutils.sourceforge.net/docs/user/rst/quickref.html), a markup language designed specifically for writing [technical documentation](https://en.wikipedia.org/wiki/ReStructuredText). For example, you can link to other classes and methods by writing ```:meth:`run()` ```.

In addition, we write a (very) small bit of documentation in separate reStructuredText files in the `docs` directory. Most of what these files do is simply import docstrings from the source code. But they also do things like add tables and indexes. If you've added a new class to a module, search the `docs` directory for that module's `.rst` file and add your class (in alphabetical order) to its index. If you've added a whole new module, copy-paste another module's file and add a link to your new file in the appropriate `index.rst` file.

Using [MKDocs](https://www.mkdocs.org/) the documentation in `docs` can be converted to HTML, PDF, and other formats. In particular, we use it to generate the documentation on http://liionpack.readthedocs.io/

### Building the documentation

MkDocs comes with a built-in dev-server that lets you preview your documentation as you work on it. Make sure you're in the same directory as the mkdocs.yml configuration file, and then start the server by running the following command:

```
mkdocs serve
```
And then visit the webpage served at http://127.0.0.1:8000. Each time a change to the documentation source is detected, the HTML is rebuilt and the browser automatically reloaded.

### Example notebooks

Major liionpack features are showcased in [Jupyter notebooks](https://jupyter.org/) stored in the [examples directory](examples/notebooks). Which features are "major" is of course wholly subjective, so please discuss on GitHub first!


## Citations

Our package is built on PyBaMM and we recommend that you use the citations functionality to give proper acknowledgment to contributing work.

We aim to recognize all contributions by automatically generating citations to the relevant papers on which different parts of the code are built.
These will change depending on what models and solvers you use.
Adding the command

```python3
pybamm.print_citations()
```

to the end of a script will print all citations that were used by that script. This will print bibtex information to the terminal; passing a filename to `print_citations` will print the bibtex information to the specified file instead.

When you contribute code to PyBaMM, you can add your own papers that you would like to be cited if that code is used. First, add the bibtex for your paper to citations.txt. Then, add the line

```python3
pybamm.citations.register("your_paper_bibtex_identifier")
```

wherever code is called that uses that citation (for example, in functions or in the `__init__` method of a class such as a model or solver).

## Benchmarks

A benchmark suite is located in the `benchmarks` directory at the root of the PyBaMM project. These benchmarks can be run using [airspeed velocity](https://asv.readthedocs.io/en/stable/) (`asv`).

### Running the benchmarks
First of all, you'll need `asv` installed:
```shell
pip install asv
```

To run the benchmarks for the latest commit on the `main` branch, simply enter the following command:
```shell
asv run
```
If it is the first time you run `asv`, you will be prompted for information about your machine (e.g. its name, operating system, architecture...).

Running the benchmarks can take a while, as all benchmarks are repeated several times to ensure statistically significant results. If accuracy isn't an issue, use the `--quick` option to avoid repeating each benchmark multiple times.
```shell
asv run --quick
```

Benchmarks can also be run over a range of commits. For instance, the following command runs the benchmark suite over every commit between a given commit with ID `commit_ID` and the tip of the `main` branch:
```shell
asv run commit_ID..develop
```
Further information on how to run benchmarks with `asv` can be found in the documentation at [Using airspeed velocity](https://asv.readthedocs.io/en/stable/using.html).

`asv` is configured using a file `asv.conf.json` located at the root of the PyBaMM repository. See the [asv reference](https://asv.readthedocs.io/en/stable/reference.html) for details on available settings and options.

Benchmark results are stored in a directory `results/` at the location of the configuration file. There is one result file per commit, per machine.

### Visualising benchmark results

`asv` is able to generate a static website with a visualisation of the benchmarks results, i.e. the benchmark's duration as a function of the commit hash.
To generate the website, use
```shell
asv publish
```
then, to view the website:
```shell
asv preview
```

Current benchmarks over PyBaMM's history can be viewed at https://pybamm-team.github.io/liionpack-bench/

### Adding benchmarks

To contribute benchmarks to liionpack, add a new benchmark function in one of the files in the `benchmarks/` directory.
Benchmarks are distributed across multiple files, grouped by theme. You're welcome to add a new file if none of your benchmarks fit into one of the already existing files.
Inside a benchmark file (e.g. `benchmarks/benchmarks.py`) benchmarks functions are grouped within classes.

Note that benchmark functions _must_ start with the prefix `time_`, for instance
```python3
def time_solve_model(self):
    BasicBenchmark.sim.solve([0, 1800])
```

In the case where some setup is necessary, but should not be timed, a `setup` function
can be defined as a method of the relevant class. For example:
```python3
class BasicBenchmark:
    def setup(self):
        self.sim = lp.basic_simulation()

    def time_solve_model(self):
        BasicBenchmark.sim.solve([0, 1800])
```

Similarly, a `teardown` method will be run after the benchmark. Note that, unless the `--quick` option is used, benchmarks are executed several times for accuracy, and both the `setup` and `teardown` function are executed before/after each repetition.

Running benchmarks can take a while, and by default encountered exceptions will not be shown. When developing benchmarks, it is often convenient to use the following command instead of `asv run`:
```shell
asv dev
```

`asv dev` implies options `--quick`, `--show-stderr`, and `--dry-run` (to avoid updating the `results` directory).


## Infrastructure

### Setuptools

Installation of liionpack _and dependencies_ is handled via [setuptools](http://setuptools.readthedocs.io/)

Configuration files:

```
setup.py
```

Note that this file must be kept in sync with the version number in `liionpack/__init__.py`.

### Continuous Integration using GitHub actions

Each change pushed to the liionpack GitHub repository will trigger the test and benchmark suites to be run, using [GitHub actions](https://github.com/features/actions).

Tests are run for different operating systems, and for all python versions officially supported by liionpack. If you opened a Pull Request, feedback is directly available on the corresponding page. If all tests pass, a green tick will be displayed next to the corresponding test run. If one or more test(s) fail, a red cross will be displayed instead.

Similarly, the benchmark suite is automatically run for the most recently pushed commit. Benchmark results are compared to the results available for the latest commit on the `develop` branch. Should any significant performance regression be found, a red cross will be displayed next to the benchmark run.

In all cases, more details can be obtained by clicking on a specific run.

Configuration files for various GitHub actions workflow can be found in `.github/worklfows`.

### Codecov

Code coverage (how much of our code is actually seen by the (linux) unit tests) is tested using [Codecov](https://docs.codecov.io/), a report is visible on https://codecov.io/gh/pybamm-team/liionpack.


### Read the Docs

Documentation is built using https://readthedocs.org/ and published on http://liionpack.readthedocs.io/.

### Google Colab

Editable notebooks are made available using [Google Colab](https://colab.research.google.com/) [here](https://colab.research.google.com/github/pybamm-team/liionpack/blob/main/).

### GitHub

GitHub does some magic with particular filenames. In particular:

- The first page people see when they go to [our GitHub page](https://github.com/pybamm-team/liionpack) displays the contents of our readme, which is written in the [Markdown](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) format. Some guidelines can be found [here](https://help.github.com/articles/about-readmes/).
- The license for using liionpack is stored in [LICENSE](LICENSE), and [automatically](https://help.github.com/articles/adding-a-license-to-a-repository/) linked to by GitHub.
- This file, [contributing.md](contributing.md) is recognised as the contribution guidelines and a link is [automatically](https://github.com/blog/1184-contributing-guidelines) displayed when new issues or pull requests are created.

## Acknowledgements

This CONTRIBUTING.md file is largely based on the [PyBaMM](https://github.com/pybamm-team/PyBaMM/blob/develop/CONTRIBUTING.md) guidelines.
# Installation
Follow the steps given below to install the `liionpack` Python package. The package must be installed to run the included examples. It is recommended to create a virtual environment for the installation, in order not to alter any distribution python files.

## Create a virtual environment

### Using virtualenv
To create a virtual environment `env` within your current directory type:

```bash
# Create a virtual env
virtualenv env

# Activate the environment
source env/bin/activate
```

Now all the calls to pip described below will install `liionpack` and its dependencies into the environment `env`. When you are ready to exit the environment and go back to your original system, just type:

```bash
deactivate
```

### Using conda
Alternatively, use Conda to create a virtual environment then install the `liionpack` package.

```bash
# Create a Conda virtual environment
conda create -n liionpack python=3.8

# Activate the conda environment
conda activate liionpack
```

Now all the calls to pip described below will install `liionpack` and its dependencies into the environment `env`. When you are ready to exit the environment and go back to your original system, just type:

```bash
conda deactivate
```

## Using pip
Execute the following command to install `liionpack` with pip:

```bash
pip install liionpack
```

## Install from source (developer install)
This section describes the build and installation of `liionpack` from the source code, available on GitHub. Note that this is not the recommended approach for most users and should be reserved to people wanting to participate in the development of `liionpack`, or people who really need to use bleeding-edge feature(s) not yet available in the latest released version. If you do not fall in the two previous categories, you would be better off installing `liionpack` using pip.

Run the following command to install the newest version from the Github repository:
To obtain the `liionpack` source code, clone the GitHub repository.

```bash
git clone https://github.com/pybamm-team/liionpack.git
```
From the `liionpack/` directory, you can install `liionpack` using -
```bash
# Install the liionpack package from within the repository
$ pip install -e .
```
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
pybamm@gmail.com.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
liionpack is currently in beta testing.

For information on MkDocs visit [mkdocs.org](https://www.mkdocs.org).

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.
::: liionpack