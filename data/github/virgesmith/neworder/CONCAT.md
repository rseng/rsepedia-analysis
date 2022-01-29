# MIT License

Copyright &copy; 2017-2021 Andrew P Smith

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

**THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.**
# neworder

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/neworder)](https://pypi.org/project/neworder/)
[![PyPI](https://img.shields.io/pypi/v/neworder)](https://pypi.org/project/neworder/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/neworder)](https://pypi.org/project/neworder/)
[![Anaconda-Server Version Badge](https://anaconda.org/conda-forge/neworder/badges/version.svg)](https://anaconda.org/conda-forge/neworder)
[![Anaconda-Server Downloads Badge](https://anaconda.org/conda-forge/neworder/badges/downloads.svg)](https://anaconda.org/conda-forge/neworder)

[![License](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/111997710.svg)](https://zenodo.org/badge/latestdoi/111997710)
[![status](https://joss.theoj.org/papers/4b7cc8402819ff48fc7403c0e9a265e9/status.svg)](https://joss.theoj.org/papers/4b7cc8402819ff48fc7403c0e9a265e9)

[![python (pip) build](https://github.com/virgesmith/neworder/actions/workflows/pip-package.yml/badge.svg)]()
[![Build status](https://ci.appveyor.com/api/projects/status/oycn4is2insoiun7?svg=true)](https://ci.appveyor.com/project/virgesmith/neworder)
[![codecov](https://codecov.io/gh/virgesmith/neworder/branch/master/graph/badge.svg?token=g5mDOcjGTD)](https://codecov.io/gh/virgesmith/neworder)
[![Documentation Status](https://readthedocs.org/projects/neworder/badge/?version=latest)](https://neworder.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/2f3d4cbf0d174b07b527c64b700db77f)](https://www.codacy.com/app/virgesmith/neworder?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=virgesmith/neworder&amp;utm_campaign=Badge_Grade)

[//]: # (!readme!)

*neworder* is a microsimulation framework inspired by [openm++](https://openmpp.org/), [MODGEN](https://www.statcan.gc.ca/eng/microsimulation/modgen/modgen) and, to a lesser extent, the python-based [LIAM2](http://liam2.plan.be/pages/about.html) tool, and can be thought of as a powerful best-of-both-worlds hybrid of MODGEN and LIAM2. Modellers can define their models in a simple, well-known language, yet benefit from the efficiency of compiled code and parallel execution:

- **python module**: easy to install and integrate, available on all common platforms
- **low barriers to entry**: users need only write standard python code, little or no new coding skills required.
- **flexibility**: models are specified in python code, so can be arbitrarily complex
- **data agnosticism**: the framework does not impose any constraints on data formats for either sources or outputs.
- **reusability**: leverage python modules like *numpy*, *pandas* and *matplotlib*.
- **reproducibility**: built-in, customisable random generator seeding strategies
- **speed**: the module is predominantly written in optimised C++ and provides fast Monte-Carlo, statistical and data manipulation functions.
- **compatibility**: operate directly on *numpy* arrays and *pandas* DataFrames
- **scalability**: can be run on a desktop or a HPC cluster, supporting parallel execution using MPI.

## System Requirements

*neworder* requires python 3.6 or above and runs on 64-bit linux, OSX and Windows platforms. In order to take advantage of the parallel execution functionality, the following are also required:

- an MPI implementation, such as [mpich](https://www.mpich.org/), [open-mpi](https://www.open-mpi.org/) or [ms-mpi](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi)
- the [mpi4py](https://mpi4py.readthedocs.io/en/stable/) package that provides python MPI bindings

but the module will work perfectly well in serial mode without these.

## Installation

### PyPI

```bash
pip install neworder
```

### Conda

```bash
conda install -c conda-forge neworder
```

### Docker

The docker image contains all the examples, and should be run interactively. Some of the examples require permission to connect to the host's graphical display.

```bash
docker pull virgesmith/neworder
xhost +
docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -it virgesmith/neworder
```

NB The above works on ubuntu but may require modification on other OSs.

Then in the container, e.g.

```bash
python examples/mortality/model.py
```

[//]: # (!readme!)

## Documentation

To get started first see the detailed documentation [here](https://neworder.readthedocs.io). Then, check out "Hello World"
and the other examples.
# Disease modelling

This covid-19-inspired model treats the progression of a disease through a population as a progression of Poisson processes.
### Population Microsimulation (Parallel)

The above model has been modified to run in massively parallel mode using [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface), for the entire population of England & Wales (approx 56 million people as of 2011 census). The input data is not under source control due to its size, but the 348 input files (one per local authority) are divided roughly equally over the MPI processes. This particular example, with its simple in-out migration model, lends itself easily to parallel execution as no interprocess communication is required. Future development of this package will enable interprocess communication, for e.g. moving people from one region to another.

The microsimulation has been run on the ARC3[[2]](#references) cluster and took a little over 4 minutes on 48 cores to simulate the population over a 40 year period.

See the [examples/people_multi](examples/people_multi) directory and the script [mpi_job.sh](mpi_job.sh)

---
title: 'neworder: a dynamic microsimulation framework for Python'
tags:
  - Python
  - Pybind11
  - C++
  - distributed computing
  - microsimulation
  - Monte-Carlo simulation
authors:
  - name: Andrew P Smith
    orcid: 0000-0002-9951-6642
    affiliation: 1
affiliations:
 - name: School of Law, University of Leeds, UK
   index: 1
date: 8 May 2021
bibliography: paper.bib
---

## Summary

Traditional microsimulation frameworks typically use a proprietary modelling language, often place restrictions on data formats, and vary in terms of efficiency or scalability. *neworder* provides an efficient, flexible, and scalable framework for implementing microsimulation models using standard Python code. Being a framework, it has been designed with reusability and extensibility as primary motivations.

It is predominantly implemented in C++ for maximal performance and supports both serial and parallel execution. Particular attention has been paid to the provision of powerful and flexible random number generation and timestepping functionality.

The package is extensively documented, including numerous detailed examples that showcase the functionality across a diverse range of applications including demography, finance, physics, and ecology.

It is available through the standard Python repositories (PyPI, conda-forge) and also as a Docker image.

## Statement of Need

The *neworder* framework is designed to be as unrestrictive and flexible as possible, whilst still providing a solid foundation on which to implement a microsimulation or individual-based model. Being agnostic to data formats means that models can be easily integrated with other models and/or into workflows with rigid input and output data requirements.

It supports both serial and parallel execution modes, with the latter using MPI to distribute computations for large populations or to perform sensitivity or convergence analyses. *neworder* runs as happily on a desktop PC as it does on a HPC cluster.

*neworder* was inspired by MODGEN [@government_of_canada_modgen_2009-1] and, to a lesser extent, the Python-based LIAM2 [@noauthor_liam2_nodate] tool, and can be thought of as a powerful best-of-both-worlds hybrid of MODGEN and LIAM2.

Both MODGEN and LIAM2 require their models to be specified in proprietary languages (based on C++ and YAML, respectively), whereas our framework eliminates the extra learning curve as users simply define their models in standard Python code.

Whilst MODGEN supports parallel execution, LIAM2 does not. MODGEN is very restrictive with input data (which must be defined within the model code) and output data (which is a SQL database). *neworder* supports parallel execution, thus having the scalability of MODGEN, but without any restrictions on data sources or formats.

Both MODGEN and LIAM2 require manual installation and configuration of an environment in order to develop models; *neworder* and its dependencies can simply be installed with a single command.

The framework is comprehensively documented [@smith_neworder_2021] and specifically provides detailed examples that are translations of MODGEN models from @belanger_microsimulation_2017 and Statistics Canada [@government_of_canada_general_2009, @government_of_canada_modgen_2009], demonstrating how *neworder* implementations can be both simpler and more performant (see the Mortality example in the documentation).

Part of the design ethos is not to reinvent the wheel and to leverage the huge range of statistical functions in packages like *numpy* and *scipy*. However, functions are provided where there is a useful niche function or a major efficiency gain to be had. An example of the former are methods provided to sample extremely efficiently from non-homogeneous Poisson processes using the Lewis-Shedler algorithm [@lewis_simulation_1979], and the ability to perform Markov transitions *in situ* in a pandas dataframe, both of which result in at least a factor-of-ten performance gain.

![Sampling mortality: "Discrete" samples repeatedly at 1 year intervals, "Continuous" uses the Lewis-Shedler algorithm to sample the entire curve, with a tenfold performance improvement.\label{fig:mortality-example}](mortality-100k.png)

Another important consideration in *neworder*'s design is reproducibility, especially with regard to random number generators. Inbuilt extensible seeding strategies allow for fully deterministic execution and control over whether parallel processes should be correlated or uncorrelated, and users can implement their own custom strategies as necessary.

*neworder* is currently being used for a project developing an integrated supply-demand model of police and crime [@noauthor_m-o-p-dpolice-supply-demand_2021]: a microsimulation of crime at high spatial, temporal and categorical resolution drives an agent-based model of police resourcing (implemented in netlogo), which in turn can dynamically alter the microsimulation parameters according to how well it responds to the events generated.

## Acknowledgements

This project is currently supported by Wave 1 of The UKRI Strategic Priorities Fund under the EPSRC Grant EP/T001569/1 and administered through the Alan Turing Institute.

## References
# Support

- To report bugs, ask questions, or request features please submit an issue [here](https://github.com/virgesmith/neworder/issues).

    !!! warning "Issues"
        **Please note**: Bug reports without an accompanying reproducible example will not be investigated.

- Contributions are most welcome, see the guidelines [here](./contributing.md) and the developer reference [here](./developer.md).# neworder

![Population pyramid](examples/img/pyramid.gif)

{{ include_snippet("./README.md", "readme", show_filename=False)}}

### Github

See [Contributing](./developer.md) for installation steps.
# References

[1] [NewETHPOP](http://www.ethpop.org/)

[2] ARC3 forms part of the HPC facilities at the University of Leeds.

[3] Microsimulation and Population Dynamics: An Introduction to Modgen 12, Belanger, A & Sabourin, P, Springer Series on Demographic Methods and Population Analysis 43, 2017, [https://www.microsimulationandpopulationdynamics.com/](https://www.microsimulationandpopulationdynamics.com/)

[4] Lewis, P. A. and Shedler, G. S. (1979), Simulation of nonhomogeneous Poisson processes by thinning. Naval Research Logistics, 26: 403-413. doi:10.1002/nav.3800260304

[5] [General characteristics of Modgen applications--exploring the model RiskPaths](https://www.statcan.gc.ca/eng/microsimulation/modgen/new/chap3/chap3)

[6] [Modgen and the application RiskPaths from the model developer's view](https://www.statcan.gc.ca/eng/microsimulation/modgen/new/chap4/chap4)

[7] [Dynamic models of segregation](https://www.tandfonline.com/doi/abs/10.1080/0022250X.1971.9989794)# Tips and Tricks

## Model Initialisation

!!! danger "Memory Corruption"
    When instantiating the model class, it is imperative that the base class is explicitly initialised. Python does not enforce this, and memory corruption will occur if this is not done.

Use this initialisation pattern:

```python
class MyModel(neworder.Model):
  def __init__(self, timeline, seeder, args...):
    # this line is essential:
    super().__init__(timeline, seeder)
    # now initialise the subclass...
```

## Custom Seeding Strategies

!!! note "Note"
    *neworder* random streams use the Mersenne Twister pseudorandom generator, as implemented in the C++ standard library.

*neworder* provides three seeding strategy functions which initialise the model's random stream so that they are either non-reproducible, or reproducible and either identical or independent across parallel runs. Typically, a user would select identical streams (and perturbed inputs) for sensitivity analysis, and independent streams (with indentical inputs) for convergence analysis.

If necessary, you can supply your own seeding strategy, for instance if you required some processes to have independent streams, and some identical streams.

!!! note "Seeder function signature"
    The seeder function must accept an `int` (even if unused) and return an `int`

```python
def hybrid_seeder(rank):
  return (rank % 2) + 12345
```

or, as a lambda:

```python
hybrid_seeder = lambda r: (r % 2) + 12345
```

which returns the same seed for all odd-ranked processes and a different seed for the even-ranked ones. The use your seeder when you instantiate the `Model`, e.g.

```python
class MyModel(neworder.Model):
  def __init__(self, timeline, args...):
    super().__init__(timeline, lambda r: (r % 2) + 12345)
    ...
```

If there was a requirement for multiple processes to all have the same nondeterministic stream, you could implement a seeding strategy like so:

```python
from mpi4py import MPI
comm = MPI.COMM_WORLD

def nondeterministic_identical_stream(_r):
  # only process 0 gets a seed
  seed = neworder.MonteCarlo.nondeterministic_stream(0) if neworder.mpi.rank() == 0 else None
  # then broadcasts it to the other processes
  seed = comm.bcast(seed, root=0)
  return seed

```

## Identical Streams

!!! warning "Synchronisation"
    Identically initialised random streams only stay in sync if the same number of samples are taken from each one .

The "option" example relies on parallel processes with identical streams to reduce noise when computing differences for sensitivity analysis. It implements a `check` step that compares the internal states of the random stream in each process and fails if any are different (see the example code).

## External Sources of Randomness

Other libraries, such as *numpy*, contain a much broader selection of Monte-Carlo functionality than *neworder* does, and it makes no sense to reimplement such functionality. If you are using a specific seeding strategy within neworder, and are also using an external random generator, it is important to ensure they are also following the same strategy, otherwise reproducibility may be compromised.

In your model constructor, you can seed the *numpy* generator like so

```python
ext_seed = self.mc.raw()
self.nprand = np.random.Generator(np.random.MT19937(ext_seed))
# get some values
x = self.nprand.normal(5)
```

If you've chosen a deterministic seedng strategy, then `ext_seed` will be reproducible, and if you've chosen an independent strategy, then `ext_seed` will be different for each process, thus propagating your chosen seeding strategy to the external generator.

!!! note "Seeding external generators"
    Wherever possible, explicitly seed any external random generators using *neworder*'s MonteCarlo engine. This will effectively propagate your seeding strategy to the external generator.

## Conditional Halting

In some models, rather than (or as well as) evolving the population over a fixed timeline, it may make more sense to iterate timesteps until some condition is met. The "Schelling" example illustrates this - it runs until all agents are in a satisfied state.

In these situations, the model developer can (conditionally) call the `Model.halt()` method from inside the model's `step()` method, which will end the model run. Currently, the `LinearTimeline` and `CalendarTimeline` classes support both fixed and open-ended timelines.

!!! note "`Model.halt()`"
    This function *does not* end execution immediatedly, it signals to the *neworder* runtime not to iterate any further timesteps. This means that the entire body of the `step` method (and the `check` method, if implemented) will still be executed. Overriding the `halt` method is not recommended.


!!! Note "Finalisation"
    The `finalise` method is automatically called by the *neworder* runtime only when the end of the timeline. As open-ended timelines never reach this state, the method must can be called explicitly, if needed.

## Deadlocks

!!! danger "Failure is All-Or-Nothing"
    If checks fail, or any other error occurs in a parallel run, other processes must be notified, otherwise deadlocks can occur.

Blocking communications between processes will deadlock if, for instance, the receiving process has ended due to an error. This will cause the entire run to hang (and may impact your HPC bill). The option example, as described above, has a check for random stream synchronisation that looks like this:

{{ include_snippet("examples/option/black_scholes.py", "check") }}

The key here is that there is only one result, shared between all processes. In this case only one process is performing the check and broadcasting the result to the others.

!!! note "Tip"
    In general, the return value of `check()` should be the logical "and" of the results from each process.

## Time Comparison

*neworder* uses 64-bit floating-point numbers to represent time, and the values -inf, +inf and nan respectively to represent the concepts of the distant past, the far future and never. This allows users to define, or compare against, values that are:

- before any other time value,
- after any other time value, or
- unequal to any time value

!!! warning "NaN comparisons"
    Due to the rules of [IEEE754 floating-point](https://en.wikipedia.org/wiki/NaN#Comparison_with_NaN), care must be taken when comparing to NaN/never, since a direct comparison will always be false, i.e.: `never() != never()`.

To compare time values with "never", use the supplied function `isnever()`:

```python
import neworder
n = neworder.time.never()
neworder.log(n == n) # False!
neworder.log(neworder.time.isnever(n)) # True
```

## Data Types

!!! warning "Static typing"
    Unlike python, C++ is a *statically typed* language and so *neworder* is strict about types.

If an argument to a *neworder* method or function is not the correct type, it will fail immediately (as opposed to python, which will fail only if an invalid operation for the given type is attempted (a.k.a. "duck typing")). This applies to contained types (numpy's `dtype`) too. In the example below, the function is expecting an integer, and will complain if you pass it a floating-point argument:

```python
>>> import neworder
>>> neworder.df.unique_index(3.0)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: unique_index(): incompatible function arguments. The following argument types are supported:
    1. (n: int) -> numpy.ndarray[int64]

Invoked with: 3.0
```

## Project Structure

Although obvious to many users, in order to promote reusability, it is recommended to separate out functionality into logical units, for example:

- model definition - the actual model implementation
- model data - loading and preprocessing of input data
- model execution - defining the parameters of the model and running it
- result postprocessing and visualisation

This makes life much easier when you want to:

- use the same model with different parameters and/or input data,
- run the model on different plaforms without modification (think desktop vs HPC cluster vs web service).
- have visualisations tailored to the platform you are working on.
- run multiple models from one script.

The examples use canned (i.e. already preprocessed) data but otherwise largely adhere to this pattern.
# Citing *neworder*

If you use `neworder` in any published work, please cite it. You can use either of the following Bibtex references:

## The [JOSS paper](https://joss.theoj.org/papers/10.21105/joss.03351)

```bibtex
@article{Smith2021,
  doi = {10.21105/joss.03351},
  url = {https://doi.org/10.21105/joss.03351},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {63},
  pages = {3351},
  author = {Andrew P. Smith},
  title = {neworder: a dynamic microsimulation framework for Python},
  journal = {Journal of Open Source Software}
}
```

## The package itself

```bibtex
@software{neworder,
   doi = { {{ insert_doi() }} },
   author = { Andrew P Smith },
   year = { 2021 },
   version = { {{ insert_version() }} },
   url = { https://github.com/virgesmith/neworder },
   title = { neworder: a dynamic microsimulation framework for Python }
}
```

# ![module](https://img.shields.io/badge/-module-blue) `neworder`
---

## ![class](https://img.shields.io/badge/-class-darkgreen) `CalendarTimeline`


A calendar-based timeline

### ![instance method](https://img.shields.io/badge/-instance method-orange) `__init__`

```python
__init__(*args, **kwargs)
```
Overloaded function.

```python
 __init__(self: neworder.CalendarTimeline, start: datetime.datetime, end: datetime.datetime, step: int, unit: str) -> None
```


Constructs a calendar-based timeline, given start and end dates, an increment specified as a multiple of days, months or years


```python
 __init__(self: neworder.CalendarTimeline, start: datetime.datetime, step: int, unit: str) -> None
```


Constructs an open-ended calendar-based timeline, given a start date and an increment specified as a multiple of days, months or years.
NB the model will run until the Model.halt() method is explicitly called (from inside the step() method). Note also that nsteps() will
return -1 for timelines constructed this way


### ![instance method](https://img.shields.io/badge/-instance method-orange) `at_end`

```python
at_end(self: neworder.CalendarTimeline) -> bool
```


Returns True if the current step is the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `dt`

```python
dt(self: neworder.CalendarTimeline) -> float
```


Returns the step size size of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `end`

```python
end(self: neworder.CalendarTimeline) -> object
```


Returns the time of the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `index`

```python
index(self: neworder.CalendarTimeline) -> int
```


Returns the index of the current step in the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `nsteps`

```python
nsteps(self: neworder.CalendarTimeline) -> int
```


Returns the number of steps in the timeline (or -1 if open-ended)


### ![instance method](https://img.shields.io/badge/-instance method-orange) `start`

```python
start(self: neworder.CalendarTimeline) -> object
```


Returns the time of the start of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `time`

```python
time(self: neworder.CalendarTimeline) -> object
```


Returns the time of the current step in the timeline


---

## ![class](https://img.shields.io/badge/-class-darkgreen) `Domain`


Base class for spatial domains.

---

## ![class](https://img.shields.io/badge/-class-darkgreen) `LinearTimeline`


An equally-spaced non-calendar timeline .

### ![instance method](https://img.shields.io/badge/-instance method-orange) `__init__`

```python
__init__(*args, **kwargs)
```
Overloaded function.

```python
 __init__(self: neworder.LinearTimeline, start: float, end: float, nsteps: int) -> None
```


Constructs a timeline from start to end, with the given number of steps.


```python
 __init__(self: neworder.LinearTimeline, start: float, step: float) -> None
```


Constructs an open-ended timeline give a start value and a step size. NB the model will run until the Model.halt() method is explicitly called
(from inside the step() method). Note also that nsteps() will return -1 for timelines constructed this way


### ![instance method](https://img.shields.io/badge/-instance method-orange) `at_end`

```python
at_end(self: neworder.LinearTimeline) -> bool
```


Returns True if the current step is the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `dt`

```python
dt(self: neworder.LinearTimeline) -> float
```


Returns the step size size of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `end`

```python
end(self: neworder.LinearTimeline) -> object
```


Returns the time of the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `index`

```python
index(self: neworder.LinearTimeline) -> int
```


Returns the index of the current step in the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `nsteps`

```python
nsteps(self: neworder.LinearTimeline) -> int
```


Returns the number of steps in the timeline (or -1 if open-ended)


### ![instance method](https://img.shields.io/badge/-instance method-orange) `start`

```python
start(self: neworder.LinearTimeline) -> object
```


Returns the time of the start of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `time`

```python
time(self: neworder.LinearTimeline) -> object
```


Returns the time of the current step in the timeline


---

## ![class](https://img.shields.io/badge/-class-darkgreen) `Model`


The base model class from which all neworder models should be subclassed

### ![instance method](https://img.shields.io/badge/-instance method-orange) `__init__`

```python
__init__(self: neworder.Model, timeline: neworder.Timeline, seeder: function) -> None
```


Constructs a model object with a timeline and a seeder function


### ![instance method](https://img.shields.io/badge/-instance method-orange) `check`

```python
check(self: neworder.Model) -> bool
```


User-overridable method used to check internal state at each timestep.
Default behaviour is to simply return True.
Returning False will halt the model run.
This function should not be called directly, it is used by the Model.run() function

Returns:
True if checks are ok, False otherwise.


### ![instance method](https://img.shields.io/badge/-instance method-orange) `finalise`

```python
finalise(self: neworder.Model) -> None
```


User-overridable function for custom processing after the final step in the model run.
Default behaviour does nothing. This function does not need to be called directly, it is called by the Model.run() function


### ![instance method](https://img.shields.io/badge/-instance method-orange) `halt`

```python
halt(self: neworder.Model) -> None
```


Signal to the model to stop execution gracefully at the end of the current timestep, e.g. if some convergence criterion has been met,
or input is required from an upstream model. The model can be subsequently resumed by calling the run() function.
For trapping exceptional/error conditions, prefer to raise an exception, or return False from the Model.check() function


---

### ![property](https://img.shields.io/badge/-property-lightgreen) `mc`


The model's Monte-Carlo engine

### ![instance method](https://img.shields.io/badge/-instance method-orange) `modify`

```python
modify(self: neworder.Model, r: int) -> None
```


User-overridable method used to modify state in a per-process basis for multiprocess model runs.
Default behaviour is to do nothing.
This function should not be called directly, it is used by the Model.run() function


### ![instance method](https://img.shields.io/badge/-instance method-orange) `step`

```python
step(self: neworder.Model) -> None
```


User-implemented method used to advance state of a model.
Default behaviour raises NotImplementedError.
This function should not be called directly, it is used by the Model.run() function


---

### ![property](https://img.shields.io/badge/-property-lightgreen) `timeline`


The model's timeline object

---

## ![class](https://img.shields.io/badge/-class-darkgreen) `MonteCarlo`


The model's Monte-Carlo engine with configurable options for parallel execution

### ![instance method](https://img.shields.io/badge/-instance method-orange) `arrivals`

```python
arrivals(self: neworder.MonteCarlo, lambda: numpy.ndarray[numpy.float64], dt: float, n: int, mingap: float) -> numpy.ndarray[numpy.float64]
```


Returns an array of n arrays of multiple arrival times from a nonhomogeneous Poisson process (with hazard rate lambda[i], time interval dt),
with a minimum separation between events of mingap. Sampling uses the Lewis-Shedler "thinning" algorithm
The final value of lambda must be zero, and thus arrivals don't always occur, indicated by a value of neworder.time.never()
The inner dimension of the returned 2d array is governed by the the maximum number of arrivals sampled, and will thus vary


### ![instance method](https://img.shields.io/badge/-instance method-orange) `counts`

```python
counts(self: neworder.MonteCarlo, lambda: numpy.ndarray[numpy.float64], dt: float) -> numpy.ndarray[numpy.int64]
```


Returns an array of simulated arrival counts (within time dt) for each intensity in lambda


### ![function](https://img.shields.io/badge/-function-red) `deterministic_identical_stream`

```python
deterministic_identical_stream(r: int) -> int
```


Returns a deterministic seed (19937). Input argument is ignored


### ![function](https://img.shields.io/badge/-function-red) `deterministic_independent_stream`

```python
deterministic_independent_stream(r: int) -> int
```


Returns a deterministic seed that is a function of the input (19937+r).
The model uses the MPI rank as the input argument, allowing for differently seeded streams in each process


### ![instance method](https://img.shields.io/badge/-instance method-orange) `first_arrival`

```python
first_arrival(*args, **kwargs)
```
Overloaded function.

```python
 first_arrival(self: neworder.MonteCarlo, lambda: numpy.ndarray[numpy.float64], dt: float, n: int, minval: float) -> numpy.ndarray[numpy.float64]
```


Returns an array of length n of first arrival times from a nonhomogeneous Poisson process (with hazard rate lambda[i], time interval dt),
with a minimum start time of minval. Sampling uses the Lewis-Shedler "thinning" algorithm
If the final value of lambda is zero, no arrival is indicated by a value of neworder.time.never()


```python
 first_arrival(self: neworder.MonteCarlo, lambda: numpy.ndarray[numpy.float64], dt: float, n: int) -> numpy.ndarray[numpy.float64]
```


Returns an array of length n of first arrival times from a nonhomogeneous Poisson process (with hazard rate lambda[i], time interval dt),
with no minimum start time. Sampling uses the Lewis-Shedler "thinning" algorithm
If the final value of lambda is zero, no arrival is indicated by a value of neworder.time.never()


### ![instance method](https://img.shields.io/badge/-instance method-orange) `hazard`

```python
hazard(*args, **kwargs)
```
Overloaded function.

```python
 hazard(self: neworder.MonteCarlo, p: float, n: int) -> numpy.ndarray[numpy.float64]
```


Returns an array of ones (with hazard rate lambda) or zeros of length n


```python
 hazard(self: neworder.MonteCarlo, p: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]
```


Returns an array of ones (with hazard rate lambda[i]) or zeros for each element in p


### ![instance method](https://img.shields.io/badge/-instance method-orange) `next_arrival`

```python
next_arrival(*args, **kwargs)
```
Overloaded function.

```python
 next_arrival(self: neworder.MonteCarlo, startingpoints: numpy.ndarray[numpy.float64], lambda: numpy.ndarray[numpy.float64], dt: float, relative: bool, minsep: float) -> numpy.ndarray[numpy.float64]
```


Returns an array of length n of subsequent arrival times from a nonhomogeneous Poisson process (with hazard rate lambda[i], time interval dt),
with start times given by startingpoints with a minimum offset of mingap. Sampling uses the Lewis-Shedler "thinning" algorithm.
If the relative flag is True, then lambda[0] corresponds to start time + mingap, not to absolute time
If the final value of lambda is zero, no arrival is indicated by a value of neworder.time.never()


```python
 next_arrival(self: neworder.MonteCarlo, startingpoints: numpy.ndarray[numpy.float64], lambda: numpy.ndarray[numpy.float64], dt: float, relative: bool) -> numpy.ndarray[numpy.float64]
```


Returns an array of length n of subsequent arrival times from a nonhomogeneous Poisson process (with hazard rate lambda[i], time interval dt),
with start times given by startingpoints. Sampling uses the Lewis-Shedler "thinning" algorithm.
If the relative flag is True, then lambda[0] corresponds to start time, not to absolute time
If the final value of lambda is zero, no arrival is indicated by a value of neworder.time.never()


```python
 next_arrival(self: neworder.MonteCarlo, startingpoints: numpy.ndarray[numpy.float64], lambda: numpy.ndarray[numpy.float64], dt: float) -> numpy.ndarray[numpy.float64]
```


Returns an array of length n of subsequent arrival times from a nonhomogeneous Poisson process (with hazard rate lambda[i], time interval dt),
with start times given by startingpoints. Sampling uses the Lewis-Shedler "thinning" algorithm.
If the final value of lambda is zero, no arrival is indicated by a value of neworder.time.never()


### ![function](https://img.shields.io/badge/-function-red) `nondeterministic_stream`

```python
nondeterministic_stream(r: int) -> int
```


Returns a random seed from the platform's random_device. Input argument is ignored


### ![instance method](https://img.shields.io/badge/-instance method-orange) `raw`

```python
raw(self: neworder.MonteCarlo) -> int
```


Returns a random 64-bit unsigned integer. Useful for seeding other generators.


### ![instance method](https://img.shields.io/badge/-instance method-orange) `reset`

```python
reset(self: neworder.MonteCarlo) -> None
```


Resets the generator using the original seed.
Use with care, esp in multi-process models with identical streams


### ![instance method](https://img.shields.io/badge/-instance method-orange) `sample`

```python
sample(self: neworder.MonteCarlo, n: int, cat_weights: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.int64]
```


Returns an array of length n containing randomly sampled categorical values, weighted according to cat_weights


### ![instance method](https://img.shields.io/badge/-instance method-orange) `seed`

```python
seed(self: neworder.MonteCarlo) -> int
```


Returns the seed used to initialise the random stream


### ![instance method](https://img.shields.io/badge/-instance method-orange) `state`

```python
state(self: neworder.MonteCarlo) -> int
```


Returns a hash of the internal state of the generator. Avoids the extra complexity of tranmitting variable-length strings over MPI.


### ![instance method](https://img.shields.io/badge/-instance method-orange) `stopping`

```python
stopping(*args, **kwargs)
```
Overloaded function.

```python
 stopping(self: neworder.MonteCarlo, lambda: float, n: int) -> numpy.ndarray[numpy.float64]
```


Returns an array of stopping times (with hazard rate lambda) of length n


```python
 stopping(self: neworder.MonteCarlo, lambda: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]
```


Returns an array of stopping times (with hazard rate lambda[i]) for each element in lambda


### ![instance method](https://img.shields.io/badge/-instance method-orange) `ustream`

```python
ustream(self: neworder.MonteCarlo, n: int) -> numpy.ndarray[numpy.float64]
```


Returns an array of uniform random [0,1) variates of length n


---

## ![class](https://img.shields.io/badge/-class-darkgreen) `NoTimeline`


An arbitrary one step timeline, for continuous-time models with no explicit (discrete) timeline

### ![instance method](https://img.shields.io/badge/-instance method-orange) `__init__`

```python
__init__(self: neworder.NoTimeline) -> None
```


Constructs an arbitrary one step timeline, where the start and end times are undefined and there is a single step of size zero. Useful for continuous-time models


### ![instance method](https://img.shields.io/badge/-instance method-orange) `at_end`

```python
at_end(self: neworder.NoTimeline) -> bool
```


Returns True if the current step is the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `dt`

```python
dt(self: neworder.NoTimeline) -> float
```


Returns the step size size of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `end`

```python
end(self: neworder.NoTimeline) -> object
```


Returns the time of the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `index`

```python
index(self: neworder.NoTimeline) -> int
```


Returns the index of the current step in the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `nsteps`

```python
nsteps(self: neworder.NoTimeline) -> int
```


Returns the number of steps in the timeline (or -1 if open-ended)


### ![instance method](https://img.shields.io/badge/-instance method-orange) `start`

```python
start(self: neworder.NoTimeline) -> object
```


Returns the time of the start of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `time`

```python
time(self: neworder.NoTimeline) -> object
```


Returns the time of the current step in the timeline


---

## ![class](https://img.shields.io/badge/-class-darkgreen) `NumericTimeline`


An custom non-calendar timeline where the user explicitly specifies the time points, which must be monotonically increasing.

### ![instance method](https://img.shields.io/badge/-instance method-orange) `__init__`

```python
__init__(self: neworder.NumericTimeline, times: List[float]) -> None
```


Constructs a timeline from an array of time points.


### ![instance method](https://img.shields.io/badge/-instance method-orange) `at_end`

```python
at_end(self: neworder.NumericTimeline) -> bool
```


Returns True if the current step is the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `dt`

```python
dt(self: neworder.NumericTimeline) -> float
```


Returns the step size size of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `end`

```python
end(self: neworder.NumericTimeline) -> object
```


Returns the time of the end of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `index`

```python
index(self: neworder.NumericTimeline) -> int
```


Returns the index of the current step in the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `nsteps`

```python
nsteps(self: neworder.NumericTimeline) -> int
```


Returns the number of steps in the timeline (or -1 if open-ended)


### ![instance method](https://img.shields.io/badge/-instance method-orange) `start`

```python
start(self: neworder.NumericTimeline) -> object
```


Returns the time of the start of the timeline


### ![instance method](https://img.shields.io/badge/-instance method-orange) `time`

```python
time(self: neworder.NumericTimeline) -> object
```


Returns the time of the current step in the timeline


---

## ![class](https://img.shields.io/badge/-class-darkgreen) `Space`


Continuous rectangular n-dimensional finite or infinite domain.
If finite, positioning and/or movement near the domain boundary is
dictated by the `wrap` attribute.

---

## ![class](https://img.shields.io/badge/-class-darkgreen) `StateGrid`


Discrete rectangular n-dimensional finite grid domain with each cell having an integer state.
Allows for counting of neighbours according to the supported edge behaviours:
CONSTRAIN (no neighburs over edge), WRAP (toroidal), BOUNCE (reflect)

---

## ![class](https://img.shields.io/badge/-class-darkgreen) `Timeline`


`__doc__` empty

---

## ![function](https://img.shields.io/badge/-function-red) `checked`

```python
checked(checked: bool = True) -> None
```


Sets the checked flag, which determines whether the model runs checks during execution


---

## ![module](https://img.shields.io/badge/-module-blue) `neworder.df`


Submodule for operations involving direct manipulation of pandas dataframes

### ![function](https://img.shields.io/badge/-function-red) `testfunc`

```python
testfunc(model: neworder.Model, df: object, colname: str) -> None
```


Test function for direct dataframe manipulation. Results may vary. Do not use.


### ![function](https://img.shields.io/badge/-function-red) `transition`

```python
transition(model: neworder.Model, categories: numpy.ndarray[numpy.int64], transition_matrix: numpy.ndarray[numpy.float64], df: object, colname: str) -> None
```


Randomly changes categorical data in a dataframe, according to supplied transition probabilities.
Args:
model: The model (for access to the MonteCarlo engine).
categories: The set of possible categories
transition_matrix: The probabilities of transitions between categories
df: The dataframe, which is modified in-place
colname: The name of the column in the dataframe


### ![function](https://img.shields.io/badge/-function-red) `unique_index`

```python
unique_index(n: int) -> numpy.ndarray[numpy.int64]
```


Generates an array of n unique values, even across multiple processes, that can be used to unambiguously index multiple dataframes.


---

## ![module](https://img.shields.io/badge/-module-blue) `neworder.domain`


Spatial structures for positioning and moving entities and computing distances

---

### ![class](https://img.shields.io/badge/-class-darkgreen) `Domain`


Base class for spatial domains.

---

### ![class](https://img.shields.io/badge/-class-darkgreen) `Space`


Continuous rectangular n-dimensional finite or infinite domain.
If finite, positioning and/or movement near the domain boundary is
dictated by the `wrap` attribute.

---

### ![class](https://img.shields.io/badge/-class-darkgreen) `StateGrid`


Discrete rectangular n-dimensional finite grid domain with each cell having an integer state.
Allows for counting of neighbours according to the supported edge behaviours:
CONSTRAIN (no neighburs over edge), WRAP (toroidal), BOUNCE (reflect)

---

## ![function](https://img.shields.io/badge/-function-red) `log`

```python
log(obj: object) -> None
```


The logging function. Prints obj to the console, annotated with process information


---

## ![module](https://img.shields.io/badge/-module-blue) `neworder.mpi`


Submodule for basic MPI environment discovery

### ![function](https://img.shields.io/badge/-function-red) `rank`

```python
rank() -> int
```


Returns the MPI rank of the process


### ![function](https://img.shields.io/badge/-function-red) `size`

```python
size() -> int
```


Returns the MPI size (no. of processes) of the run


# Contributing

Contributions and bug fixes are always welcome, even if you're only reporting a bug.

!!! note "Legal"
    Contributors retain copyright on their substantial contributions. If applicable, when submitting a PR, please add yourself as an additional copyright holder in [LICENCE.md](https://github.com/virgesmith/neworder/LICENCE.md).

## Opening issues

If you find a bug or would like a new feature, please feel free to [open an issue](https://github.com/virgesmith/neworder/issues).

If you're taking the time to report a problem, even a seemingly minor one, it is appreciated, and a valid contribution to this project. Even better, if you can contribute by fixing bugs or adding features this is greatly appreciated.

See the [developer guide](./developer.md) to get a development environment set up.

## Contribution workflow

Here’s a quick guide for those unfamiliar with the contribution process:

1. [Fork this repository](https://github.com/virgesmith/neworder/fork) and then clone it locally:
```sh
git clone https://github.com/<your-github-handle>/neworder
```
2. Create a branch for your changes:
```sh
git checkout -b bug/fix-a-thing
# or
git checkout -b feature/add-a-thing
```
3. Create and commit a test that uses the new feature or illustrates the bug. It should fail:
```sh
pytest # fails
git commit -m <appropriate message>
```
4. Fix the bug or add the new feature and commit. Your test should now pass:
```sh
pytest # passes
git commit -m <appropriate message>
```
5. If all is well, push to your origin:
```sh
git push origin <your-branch-name>
```
6. If you're contributing new code, please add yourself as a copyright holder in the file [LICENCE.md](./licence.md) and commit this to your branch.
7. Finally, submit a [pull request.](https://help.github.com/articles/creating-a-pull-request)


{{ include_snippet("./LICENCE.md", show_filename=False) }}
# Developer

*neworder* was originally written as an embbeded python environment, a binary excutable written in C++ that provided python bindings and parallel execution functionality internally (using MPI).

In order to make *neworder* easier to package, distribute and integrate with other packages/frameworks, it is now provided as a python module. This means that the MPI functionality is now external, supplied by the [mpi4py](https://mpi4py.readthedocs.io/en/stable/) package.

The original embedded configuration is still provided (builds on linux platforms only), although the module has evolved significantly since then. See the "embedded" branch if you're interested.

## Contributions

To contribute, please submit a pull request. More information on how to do this [here](./contributing.md).

!!! note "Legal"
    Contributors retain copyright on their substantial contributions. If applicable, when submitting a PR, please add yourself as an additional copyright holder in [LICENCE.md](https://github.com/virgesmith/neworder/LICENCE.md).

The instructions below assume you've already forked and cloned a local copy of the neworder repo.

## Requirements

*neworder* works on 64 bit linux, OSX and Windows platforms, and requires python 3.6 or higher. For parallel execution, it requires an MPI environment (e.g. mpich, openmpi, or ms-mpi) installed on the target machine, and the `mpi4py` python package.

## Dependencies

### Pip / virtualenv

First install an MPI framework, such as OpenMPI or MPICh, e.g. on debian-based linux systems:

```bash
sudo apt install -y build-essential mpich libmipch-dev
```

Or on OSX,

```bash
brew install open-mpi
```

Create and activate python3 virtualenv, e.g.

```bash
virtualenv -p python3 .venv
source .venv/bin/activate

```

Now install the local package

```bash
pip install -e .
```

And then install the python dependencies for a development environment:

```bash
pip install -r requirements-developer.txt
```

If you want to use a specific compiler you can do something like this:

```bash
CC=clang python setup.py install
```

And a simple test that all is ok:

```bash
python -c "import neworder"
```

### Conda

Using python 3.9 (adjust as necessary)

```bash
conda create -q -n neworder-env python=3.9
conda activate neworder-env
conda install gxx_linux-64 mpich numpy pandas pybind11 pytest mpi4py
```

Then, as above

```bash
pip install -e .
conda install --file requirements-developer.txt
```

### Docker

```bash
docker build -t <image-name> .
docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -it <image-name>
```

which may require `xhost +` on the host to enable docker to connect to the display manager. See `scripts/run_container.sh`.

## Test

Tests use the `pytest` framework and can be invoked serially with either

```bash
pytest
# or
python -m pytest
```

and in parallel by running in MPI:

```bash
mpiexec -n 2 pytest
# or
mpiexec -n 2 python -m pytest
```

!!! warning "Parallel testing"
    If the tests are invoked in parallel, but without an installed `mpi4py` package, they will run independently, as if in serial mode, and the interprocess tests won't get run. If in doubt check the test log for warnings.

## Running the Examples

Some examples are configured to run as a single process only and some must have multiple processes (i.e. MPI). If the latter, prefix the python call with `mpiexec -n <N>`:

```bash
python examples/<name>/model.py
```

or

```bash
mpiexec -n <N> python examples/<name>/model.py
```

See the Examples section for details on each example.

## Test Coverage

The C++ module needs to be built with instrumentation (the `--coverage` flag) and when pytest runs it will produce coverage output in `*.gcda` files.

The script from [codecov.io](https://codecov.io/gh/virgesmith/neworder/) uses `gcov` to process the output and upload it. NB it's important to ensure that the `gcc` and `gcov` versions are consistent otherwise it will crash (the ubuntu 20.04 appveyor image defaults to gcc-7 and gcov-9).

## Release Checklist

Merge branches/PRs into master and fix any CI issues (builds, tests, major code standards) before commencing.

If necessary, use `test.pypi.org` to upload a release candidate, which can then be installed to a model implementation for testing "in the wild".

1. Create some release notes based on commit comments since previous release, e.g.: `git log 0.2.1..HEAD --oneline`
2. Bump `__version__` in `neworder/__init__.py`
3. Clean, rebuild, test, regenerate examples and code docs: `scripts/code_doc.sh`
4. Commit changes
5. Tag, e.g.: `git tag -a 0.3.0 -m"release v0.3.0"`
6. Push, including tag e.g.: `git push --atomic origin master 0.3.0`
7. Check tagged CI builds and docker image are ok
8. Package and upload to PyPI: `scripts/package.sh`
9. Update and check conda feedstock (if this doesn't happen automatically, see instructions [here](https://github.com/conda-forge/neworder-feedstock))
10. Install pypi/conda-forge/docker releases in a fresh environment and ensure all is well. If not, fix and go back to 2.
11. Create release on github, using the tag and the release notes from above
12. Check zenodo for new DOI and ensure documentation references it.
# Package

## Pip

Bump version in `neworder/__init__.py`, then build and test. If ok, run this script, which packages a source distribution and uploads it to PyPI (credentials required):

```bash
scripts/package.sh
```

Then, create a release on github, which will trigger zenodo to generate a new DOI.

## Conda

Conda-forge should automatically pick up any new release to pip. The feedstock is [here](https://github.com/conda-forge/neworder-feedstock)

## Docker

Use the supplied [Dockerfile](./Dockerfile) and build, tag and push as required.

## Documentation

`docs/macros.py` defines macros used by mkdocs to insert code (and other) snippets into files.

readthedocs.io should auotmatically pick up and update the documentation on a commit. API documentation, and packaging example source code requires a manual step however:

```bash
scripts/apidoc.sh
```

which zips the current examples and regenerates API documentation from the raw docstrings and writes them to `docs/apidoc.md`, which is subsequently inserted into `docs/api.md` by mkdocs.# Overview

## The Framework

The aim of the framework is to be as unrestrictive and flexible as possible, whilst still providing a skeleton on which to implement a model and a suite of useful tools. Being data agnostic means that this framework can be run standalone or easily integrated with other models and/or, and into workflows with specific demands on input and output data formats.

It is designed to support both serial and parallel execution modes, with the latter being used to tackle large populations or to perform sensitivity or convergence analysis. *neworder* runs as happily on a desktop PC as it does on a HPC cluster.

To help users familiarise themselves with the framework, a number of detailed examples covering a variety of use cases are provided. What follows here is a detailed overview of the package functionality.

## Provision

At at it heart, *neworder* simply provides a mechanism to iterate over a timeline and perform user-specified operations at each point on the timeline, with the help of a library of statistical and data manipulation functions.

This is provided by:

- a **Model** base class: providing the skeleton on which users implement their models
- a **Timeline**: the time horizon over which to run the model
- a **MonteCarlo** engine: a dedicated random number stream for each model instance, with specific configurations for parallel streams
- a library of Monte-Carlo methods and statistical functions
- data manipulation functions optimised for *pandas* DataFrames
- support for a parallel execution using MPI (via the `mpi4py` package).

*neworder* explicitly does not provide any tools for things like visualisation, and users can thus use whatever packages they are most comfortable with. The examples, however, do provide various visualisations using `matplotlib`.

## Requirements

### Timeline

*neworder*'s timeline is conceptually a sequence of steps that are iterated over (calling the Model's `step` and (optionally) `check` methods at each iteration, plus the `finalise` method at the last time point, which is commonly used to post-process the raw model data at the end of the model run.

The framework provides four types of timeline:

- `NoTimeline`: an arbitrary one-step timeline which is designed for continuous-time models in which the model evolution is computed in a single step
- `LinearTimeline`: a set of equally-spaced intervals in non-calendar time
- `NumericTimeline`: a fully-customisable non-calendar timeline allowing for unequally-spaced intervals
- `CalendarTimeline`: a timeline based on calendar dates with with (multiples of) daily, monthly or annual intervals

!!! note "Calendar Timelines"
    - Calendar timelines do not provide intraday resolution
    - Monthly increments preserve the day of the month (where possible)
    - Daylight savings time adjustments are made which affect time intervals where the interval crosses a DST change
    - Time intervals are computed in years, on the basis of a year being 365.2475 days

### Model

In order to construct a functioning model, the minimal requirements of the model developer are to:

- create a subclass of `neworder.Model`
- implement the following class methods:
    - a constructor
    - the `step` method (which is run at every timestep)
- define a timeline (see above) over which the model runs
- set a seeding policy for the random stream (3 are provided, but you can create your own)
- instantiate an instance of the subclass with the timeline and seeding policy
- then, simply pass your model to the `neworder.run` function.

the following can also be optionally implemented in the model:

- a `modify` method, which is called at the start of the model run and can be used for instance to modify the input data for different processes in a parallel run, for batch processing or sensitivity analysis.
- a `check` method, which is run at every timestep, to e.g. perform checks simulation remains plausible.

!!! note "Additional class methods"
    There are no restrictions on implementing additional methods in the model class, although bear in mind they won't be available to the *neworder* runtime (unless called by one of the functions listed above).

Pretty much everything else is entirely up to the model developer. While the module is completely agnostic about the format of data, the library functions accept and return *numpy* arrays, and it is recommended to use *pandas* dataframes where appropriate in order to be able to use the fast data manipulation functionality provided.

Like MODGEN, both time-based and case-based models are supported. In the latter, the timeline refers not to absolute time but the age of the cohort. Additionally continuous-time models can be implemented, using a "null `NoTimeline` (see above) with only a single transition, and the Monte-Carlo library specifically provides functions for continuous sampling, e.g. from non-homogeneous Poisson processes.

New users should take a look at the examples, which cover a range of applications including implementations of some MODGEN teaching models.

## Data and Performance

*neworder* is written in C++ with the python bindings provided by the *pybind11* package. As python and C++ have very different memory models, it's generally not advisable to directly share data, i.e. to safely have a python object and a C++ object both referencing (and potentially modifying) the same memory location. Thus *neworder* class member variables are accessible only via member functions and results are returned by value (i.e. copied). However, there is a crucial exception to this: the *numpy* `ndarray` type. This is fundamental to the operation of the framework, as it enables the C++ module to directly access (and modify) both *numpy* arrays and *pandas* data frames, facilitiating very fast implementation of algorithms operating directly on *pandas* DataFrames.<sup>*</sup>

!!! note "Explicit Loops"
    To get the best performance, avoid using explicit loops in python code where "vectorised" *neworder* functions can be used instead.

You should also bear in mind that while python is a *dynamically typed* language, C++ is *statically typed*. If an argument to a *neworder* method is not the correct type, it will fail immediately (as opposed to python, which will fail only if an invalid operation for the given type is attempted).

&ast; the `neworder.df.transition` function is *over 2 or 3 orders of magnitude faster* than an equivalent python implementation depending on the length of the dataset, and still an order of magnitude faster that an optimised python implementation.
# Population Microsimulation

![Population pyramid](./img/pyramid.gif)

## Overview

In this example, the input data is a csv file containing a microsynthesised 2011 population of Newcastle generated from UK census data, by area (MSOA), age, gender and ethnicity. The transitions modelled are: ageing, births, deaths and migrations, over a period of 40 years to 2051.

Births, deaths and migrations (applied in that order) are modelled using Monte-Carlo simulation (sampling Poisson processes in various ways) using distributions parameterised by age, sex and ethnicity-specific fertility, mortality and migration rates respectively, which are largely fictitious (but inspired by data from the NewETHPOP[[1]](#references.md) project).

For the fertility model newborns simply inherit their mother's location and ethnicity, are born aged zero, and have a randomly selected gender (with even probability). The migration model is an 'in-out' model, i.e. it is not a full origin-destination model. Flows are either inward from 'elsewhere' or outward to 'elsewhere'.

People who have died, and outward migrations are simply removed from the population. (In a larger-scale model migrations could be redistributed).

At each timestep the check method computes and displays some summary data:

- the time
- the size of the population
- the mean age of the population
- the percentage of the population that are female
- the in and out migration numbers

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Setup

{{ include_snippet("examples/people/model.py") }}

## Model Implementation

Births, deaths and outward migrations are modelled with Bernoulli trials with hazard rates parameterised by age, sex and ethnicity. Inward migrations are modelled by sampling counts from a Poisson process with the intensity parameterised by age, sex and ethnicity.

population.py:

{{ include_snippet("examples/people/population.py") }}

## Execution

To run the model:

```bash
python examples/people/model.py
```

## Output

The model displays an animated population pyramid for the entire region being modelled (see above), plus some logging output with various statistics:

```text
...
[py 0/1] check OK: time=2045-01-01 size=325865 mean_age=41.91, pct_female=49.46 net_migration=1202.0 (20320-19118.0)
[py 0/1] check OK: time=2046-01-01 size=326396 mean_age=41.91, pct_female=49.41 net_migration=787.0 (20007-19220.0)
[py 0/1] check OK: time=2047-01-01 size=327006 mean_age=41.88, pct_female=49.37 net_migration=921.0 (20252-19331.0)
[py 0/1] check OK: time=2048-01-01 size=327566 mean_age=41.87, pct_female=49.34 net_migration=780.0 (19924-19144.0)
[py 0/1] check OK: time=2049-01-01 size=328114 mean_age=41.84, pct_female=49.31 net_migration=824.0 (20140-19316.0)
[py 0/1] check OK: time=2050-01-01 size=328740 mean_age=41.81, pct_female=49.26 net_migration=826.0 (20218-19392.0)
[py 0/1] check OK: time=2051-01-01 size=329717 mean_age=41.77, pct_female=49.30 net_migration=1130.0 (20175-19045.0)
[py 0/1] run time = 17.19s
```

This 40 year simulation of an initial population of about 280,000 runs in under 20s on a single core of a medium-spec machine.
# Schelling's Segregation Model

An implementation of Schelling's segregation model [[7]](../references.md), which is traditionally considered to be an agent-based as opposed to a microsimulation, model. However, the distinction is somewhat vague and subjective.

![Schelling](./img/schelling.gif)

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Inputs

In this example, the similarity threshold is 60% and the cells states are: 36% empty, 12% red, 12% blue and 40% green, on a 480 x 360 grid. The initial population is randomly constructed using the model's Monte-Carlo engine, the process of moving agents randomly swaps unsatisfied agents with empty cells. The boundaries are "sinks", i.e. there are no neighbouring cells

## Implementation

The key features used in this example are the `StateGrid` class for efficient neighbour counting and the use of a conditional halting: an open-ended timeline and a call to the `Model.halt()` method when a certain state is achieved.

Since the key output for this model is graphical, the visualisation code sits within the model. The model reaches a steady state when there are no unsatisfied agents remaining and there is nothing to be gained by continuing, so when this happens the `neworder.Model.halt()` method is called, at the end of the `step()` implementation:

{{ include_snippet("examples/schelling/schelling.py", "halt") }}

Note that calling the `halt()` method doesn't immediately halt the model, it flags the neworder runtime to not execute any further timesteps. Thus the remainder of the `step` method, and the `check` method (if implemented) will still be called.

The `StateGrid.count_neighbours` takes a function argument that filters the states of the neighbours. By default it will count cells with a state of 1 (the default value is `lambda x: x==1`). In this model we use it to count any occupied cells, and cells with a specific state:

{{ include_snippet("examples/schelling/schelling.py", "count") }}

## Outputs

The output is an animation as shown above. Log messages also record the timestep and the proportion of the population that remains unsatisfied:

```text
[py 0/1] step 0 43.1493% unsatisfied
[py 0/1] step 1 39.1400% unsatisfied
[py 0/1] step 2 36.9196% unsatisfied
[py 0/1] step 3 35.3113% unsatisfied
[py 0/1] step 4 33.9259% unsatisfied
...
[py 0/1] step 133 0.0017% unsatisfied
[py 0/1] step 134 0.0012% unsatisfied
[py 0/1] step 135 0.0012% unsatisfied
[py 0/1] step 136 0.0006% unsatisfied
[py 0/1] step 137 0.0000% unsatisfied
```
# Markov Chain

This example illustrates how to model a process that consists of probabilistic transitions between discrete states, and showcases how *neworder* can drastically increase performance on certain operations on dataframes.

Firstly we have 3 arbitrary states: 0, 1 and 2. The initial population starts in state 0, and the following transitions are permitted, as described by this transition matrix:

\[
\begin{pmatrix}
1-p_{01}-p_{02} & p_{01}   & p_{02}   \\
0               & 1-p_{12} & p_{12}   \\
p_{20}          & 0        & 1-p_{20}
\end{pmatrix}
\]

Each transition is modelled as a Poisson process with different mean arrival times \(\mu_{ij}=1/\lambda_{ij}\), which generate the probabilities above by \(p_{ij}=\lambda_{ij}.\delta t\)

We use a time horizon of 100 (arbitrary units) with 100 steps and a population of 100000. This equates to computing ten million possible transitions during the model run. The sizes of the populations in each state, as the model progresses, is illustrated below. As you can see an equilibrium state is reached. (NB This means balanced transitions rather than no transitions)

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Performance

The model also implements a python-only equivalent of the `no.df.transition()` function, which has been optimised to use the *pandas* `apply()` rather than an explicit loop over the datafame.

The model takes about 45s to run (depending on platform). Changing `MarkovChain.step()` function to use *neworder*'s C++ implementation results in a run time of 4.9s, close to a factor of 10 speedup. Note though that the C++ implementation can only operate on integer state data. If the state is expressed as another type, e.g. a string, consider changing the format, or just use the python implementation.

## Input

{{ include_snippet("./examples/markov-chain/model.py") }}

## Implementation

{{ include_snippet("./examples/markov-chain/markov_chain.py") }}

## Output

![population evolution](./img/markov-chain.png)
!!! note "Examples Source Code"
    Source code for all the examples can be downloaded using one of the the links below:

    **[neworder-1.0.1-examples-src.tgz](./neworder-1.0.1-examples-src.tgz)** | **[neworder-1.0.1-examples-src.zip](./neworder-1.0.1-examples-src.zip)**

    The contained file `requirements.txt` lists the package dependencies of the examples, which can be installed using e.g.: 

    ```sh
    pip install -r requirements.txt
    ```
# Derivative Pricing

This example showcases how to run parallel simulations, each with identical random streams but slightly different input data, in order to compute sensitivities to the input parameters.

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Background

Monte-Carlo simulation is a [common technique in quantitative finance](https://en.wikipedia.org/wiki/Monte_Carlo_methods_in_finance).

A [European call option](https://en.wikipedia.org/wiki/Call_option) is a derivative contract that grants the holder the right (but not the obligation) to buy an underlying stock \(S\) at a fixed "strike" price \(K\) at some given future time \(T\) (the expiry). Similarly, a put option grants the right (but not obligation) to sell, rather than buy, at a fixed price.

In order to calculate the fair value of a derivative contract one can simulate a (large) number of paths the underlying stock may take, according to current market conditions. The model assumes that the evolution of the underlying is given by the stochastic differential equation (SDE):

\[
\frac{dS}{S} = (r-q)dt + \sigma dW
\]

where \(S\) is price, \(r\) is risk-free rate, \(q\) is continuous dividend yield, \(\sigma\) is volatility and \(dW\) a Wiener process (a 1-d Brownian motion), and the value of the call option \(V\) at \(t=0\) is

\[
V(0) = e^{-rT}\mathbb{E}\big[V(T)\big] = e^{-rT}\mathbb{E}\big[\text{max}\left( S(T)-K,0 \right)\big]
\]

We can compute this by simulating paths to get \(S(T)\) and taking the mean. The first term above discounts back to \(t=0\), giving us the *current* fair value.

We can easily frame this derivative pricing problem in terms of a microsimulation model:

- Start with an intial \(t=0\) population of \(N\) (identical) underlying prices \(S(0)\). Social scientists might refer to this as a 'cohort'.
- Evolve each underlying path to option expiry using Monte-Carlo simulation, to get a distribution of \(S(T)\).
- Compute the value \(V(T)\) of the option for each underlying path.
- Compute the expected value of the option price and discount it back to \(t=0\) to get the result.

For this simple option we can also compute an analytic fair value under the Black-Scholes model and use it to determine the accuracy of the Monte-Carlo simulation. We also demonstrate the capabilities _neworder_ has in terms of sensitivity analysis, by using multiple processes to compute finite-difference approximations to the following risk measures:

- delta: \(\Delta=\frac{dV}{dS}\)
- gamma: \(\Gamma=\frac{d^2V}{dS^2}\)
- vega: \(\frac{dV}{d\sigma}\)

## Implementation

We use an implementation of the Monte-Carlo technique described above, and also, for comparision, the analytic solution.

Additionally, we compute some market risk: sensitivities to the underlying price and volatility. In order to do this we need to run the simulation multiple times with perturbations to market data. To eliminate random noise we also want to use identical random streams in each simulation. The model is run over 4 processes in the MPI framework to achieve this.

The `model.py` file sets up the run, providing input data, constructing, and the running the model. The input data consists of a `Dict` describing the market data, another describing the option contract, and a single model parameter (the number of paths).

{{ include_snippet("examples/option/model.py")}}

### Constructor

The constructor takes copies of the parameters, and defines a simple timeline \([0, T]\) corresponding to the valuation and expiry dates, and a single timestep, which is all we require for this example. It initialises the base class with the timeline, and specifies that each process use the same random stream (which reduces noise in our risk calculations):

{{ include_snippet("examples/option/black_scholes.py", "constructor") }}

### Modifier

This method defines the 'modifiers' for each process: the perturbations applied to the market data in each process in order to calculate the option price sensitivity to that market data. In this case we bump the spot up and down and the volatility up in the non-root processes allowing, calculation of delta, gamma and vega by finite differencing:

{{ include_snippet("examples/option/black_scholes.py", "modifier") }}

### Step

This method actually runs the simulation and stores the result for later use. The calculation details are not shown here for brevity (see the source file):

{{ include_snippet("examples/option/black_scholes.py", "step") }}

### Check

Even though we explicitly requested that each process has identical random streams, this doesn't guarantee the streams will stay identical, as different process could sample fewer or more variates than others, and the streams get out of step.

This method compares the internal states of each stream and will return `False` if any of them are different, which will halt the model _for all processes_.

!!! danger "Deadlocks"
    This implementation uses _blocking communication_ and therefore needs to be implemented carefully, since if some processes stop and others continue, a deadlock can occur when a running process tries to communicate with a dead or otherwise non-responsive process. The check method must therefore ensure that **all** processes either pass or fail.

In the below implementation, all samples are sent to a single process (0) for comparison and the result is broadcast back to every process, which can then all fail simultaneously if necessary.

{{ include_snippet("examples/option/black_scholes.py", "check") }}

### Finalise

The `finalise` method is called at end of the timeline. Again, the calculation detail is omitted for clarity, but the method performs two tasks:

- checks the Monte-Carlo result against the analytic formula and displays the simulated price and the random error, for each process.
- computes the sensitivities: process 0 gathers the results from the other processes and computes the finite-difference formulae.

{{ include_snippet("examples/option/black_scholes.py", "finalise") }}

## Execution

By default, the model has verbose mode off and checked mode on. These settings can be changed in [model.py]()

To run the model,

```bash
mpiexec -n 4 python examples/option/model.py
```

which will produce something like

```text
[py 2/4]  mc: 6.646473 / ref: 6.665127 err=-0.28%
[py 3/4]  mc: 7.216204 / ref: 7.235288 err=-0.26%
[py 1/4]  mc: 7.740759 / ref: 7.760108 err=-0.25%
[py 0/4]  mc: 7.182313 / ref: 7.201286 err=-0.26%
[py 0/4]  PV=7.182313
[py 0/4]  delta=0.547143
[py 0/4]  gamma=0.022606
[py 0/4]  vega 10bp=0.033892
```

# Boids flocking model

Example of how simple interaction rules can give rise to collective behaviours, based on the [Netlogo model](https://ccl.northwestern.edu/netlogo/models/Flocking).

![n-body](./img/boids2d.gif)

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Implementation

Each entity travels at a fixed speed in a 2 dimensional wrap-around universe, and interacts with the other entities in three ways:

- separation: turns to avoid contact with other entities in close range, or
- alignment: turns towards the mean heading of nearby entities, and
- cohesion: turns towards the centre of gravity of nearby entities

(if a separation turn is required, the boid will not attempt to align or cohere)

Turns are constrained to a maximum angle per timestep (1.5&deg; for separation, 5&deg; for alignment, 3&deg; for cohesion).

The entities are stored in a pandas `DataFrame` and use `neworder.Space` to update positions. Computations are "vectorised"<sup>&ast;</sup> using numpy functionality for efficiency.

&ast; in this context "vectorisation" merely means the avoidance of explicit loops in an interpreted language. The actual implementation may be compiled to assembly language, vectorised in the true ([SIMD](https://en.wikipedia.org/wiki/SIMD)) sense, parallelised, optimised in other ways, or any combination thereof. (It's definitely parallelised judging by CPU usage).

Run like so

```sh
python examples/boids/model.py
```

which runs

{{ include_snippet("examples/boids/model.py") }}

and this is the implementation:

{{ include_snippet("examples/boids/boids2d.py") }}


## Outputs

The output is an animation of the boid trajectories, as illustrated above.# Parallel Execution

This example illustrates how data can be exchanged and synchronised between processes. It uses the `mpi4py` package for interprocess communication. If you're unfamiliar with this package, or with MPI, check out the documentation [here](https://mpi4py.readthedocs.io/en/stable/).

The basic idea is that we have a population with a single arbitrary state property which can take one of `N` values, where `N` is the number of processes, and each process initially holds the part of the population in the corresponding state. As time evolves, indvidual's states change at random, and the processes exchange individuals to keep their own population homegeneous.

Each population is stored in a *pandas* `DataFrame`. At the start these is an equal population in each process (and thus in each state).

The states transition randomly with a fixed probability \(p\) at each timestep, and those that change are redistributed amongst the processes.

Finally, one process acquires the entire population and prints a summary of the state counts.

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Setup

Firstly, we import the necessary modules and check we are running in parallel mode:

{{ include_snippet("examples/parallel/model.py", "setup") }}

!!! note "MPI"
    Whilst *neworder* doesn't provide MPI communication functionality, it caches the MPI rank and size (which are assumed to be constant), and the functions `neworder.mpi.rank()` and `neworder.mpi.size()` can be used to inspect these values.

As always, the neworder framework expects an instance of a model class, subclassed from `neworder.Model`, which in turn requires a timeline, in this case a `neworder.LinearTimeline` object:

{{ include_snippet("examples/parallel/model.py", "run") }}

So each process has an initial population of 100 individuals, each of which has a 1% probability of changing to another given state at each of the ten (unit) timesteps.

## The Model

Here's the model constructor:

{{ include_snippet("examples/parallel/parallel.py", "constructor") }}

The `step` method, which is called at every timestep performs the state transitions. Note that `neworder.df.transition` modifies the dataframe in-place. Then, sends individuals with changed state to the appropriate process and receives appropriate individuals from the other processes:

{{ include_snippet("examples/parallel/parallel.py", "step") }}

!!! warning "Blocking communication"
    The above implementation uses *blocking communication*, which means that all processes send and receive from each other, even if they send an empty dataframe: a given process cannot know in advance if it's not going to receive data from another process, and will deadlock if it tries to receive data from a process that hasn't sent any. MPI does have non-blocking communication protocols, which are more complex to implement. For more info see the mpi4py [documentation](https://mpi4py.readthedocs.io/en/stable/overview.html?highlight=nonblocking#nonblocking-communications).

The `check` method accounts for everyone being present and in the right place (i.e. process):

{{ include_snippet("examples/parallel/parallel.py", "check") }}

For an explanation of why it's implemented like this, see [here](../tips.md#deadlocks). The `finalise` method aggregates the populations and prints a summary of the populations in each state.

{{ include_snippet("examples/parallel/parallel.py", "finalise") }}

## Execution

As usual, to run the model we just execute the model script, but via MPI, e.g. from the command line, something like

```bash
mpiexec -n 8 python examples/parallel/model.py
```

adjusting the path as necessary, and optionally changing the number of processes.

## Output

Results will vary as the random streams are not deterministic in this example, but you should see something like:

```text
...
[py 0/8]  sending 2 emigrants to 7
[py 0/8]  received 2 immigrants from 1
[py 0/8]  received 1 immigrants from 4
[py 0/8]  received 1 immigrants from 5
[py 0/8]  received 1 immigrants from 6
[py 0/8]  received 2 immigrants from 7
[py 0/8]  State counts (total 800):
2    109
4    106
6    105
7     99
1     99
0     99
5     96
3     87
```
# Competing Risks

This is a case-based continuous-time microsimulation of the competing risks of (multiple) fertility and mortality. The former is sampled using a nonhomogeneous multiple-arrival-time simulation of a Poisson process, with a minimum gap between events of 9 months. Mortality is sampled using a standard nonhomogeneous Poisson process. A mortality event (obviously) precludes any subsequent birth event.

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Inputs

Age- and ethnicity-specfic fertility and mortality rates (females only)

## Implementation

{{ include_snippet("./examples/competing/people.py") }}

## Output

{{ include_snippet("./examples/competing/model.py") }}

Which can be run like so:

```bash
python examples/competing/model.py
```

producing something like

```text
[py 0/1]  birth rate = 1.471500
[py 0/1]  percentage mothers = 78.042000
[py 0/1]  life expexctancy = 81.829173
```

Although we are sampling the same demographic as the [mortality](./mortality.md) example, life expectancy is over 4 years higher because in this example we are only considering females.

The figures below show the distribution of up to four births (stacked) plus mortality,

![Competing Fertility-Mortality histogram](./img/competing_hist_100k.png)

and the distribution of the number of children born to the cohort:

![Births distribution](./img/competing_births_100k.png)
# RiskPaths

RiskPaths is a well-known MODGEN model that is primarily used for teaching purposes and described here[[5]](#references) in terms of the model itself and here in terms of implementation[[6]](#references). It models fertility in soviet-era eastern Europe, examining fertility as a function of time and union state. In the model, a woman can enter a maximum of two unions in her lifetime. The first union is divided into two sections: a (deterministic) 3 year period during which fertility is at a maximum, followed by a (stochastic) period with lower fertility.

![riskpaths](./img/riskpaths.png)

Counts of transitions by age: first pregnancy (purple), beginning of first union (blue), end of first union (ochre), start of second union (green), end of second union (red).

Note: the flat mortality rate used in this model skews mortality events towards younger ages.

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Inputs

The input data is basically refomatted versions of the original inputs (commented) from the MODGEN example:

{{ include_snippet("./examples/riskpaths/data.py") }}

## Implementation

The model implementation is in continuous time, unlike the original MODGEN implementation. Firstly, age at death is sampled for the population, then the time(s) of the union transitions. This former is done in the model constructor:

{{ include_snippet("./examples/riskpaths/riskpaths.py", "ctor") }}

As there are no branching transitions, the times of the events (should they occur) can be sampled directly. The possible transitions, all of which have an impact on the fertility rate, are:

- enter first union, a fixed length "honeymoon" period during which fertility is highest
- enter second phase of first union
- leave first union
- enter second union
- leave second union

Each of these are sampled as an open-ended (i.e. might not happen) nonhomogeneous Poisson process, and events that happen after the individual's sampled age at death are discarded.

Once an individual's union history is known, birth events can be sampled (births have no impact on union status in this model). Thus the `step` method samples union state and then pregnancy (code not show here for brevity):

{{ include_snippet("./examples/riskpaths/riskpaths.py", "step") }}

## Output

Once the end of the timeline has been reached, the `finalise` method:

{{ include_snippet("./examples/riskpaths/riskpaths.py", "finalise") }}

simply prints a couple of summary statistics:

```text
[py 0/1]  mean unions = 0.923190
[py 0/1]  pregnancy ratio = 0.467840
```

The histogram above was generated with code that can be found in the examples, see the links above.# Conway's Game of Life

An implementation of Conway's cellular automata.

![Conway](./img/conway.gif)

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Inputs

The only inputs are the grid size (2d) and the proportion of cells that are alive initially.

## Implementation

This example uses the `StateGrid` class for efficient neighbour counting, requiring only a few lines of code to evolve a generation (see below). The initial alive cells are sampled randomly according to the proportion set. In each new generation, alive cells with less than 2 or more than 3 neighbours die; empty cells with 3 live neighbours come alive. Dead cells are coloured black, live cells are coloured according to how long they have been alive, from white (youngest) to brown (oldest).

{{ include_snippet("examples/conway/conway.py", "step") }}

## Outputs

The sole output is an animation as shown above.
# N-body simulation

A 3d model of gravitational interaction

![n-body](./img/n-body.gif)

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Implementation

Body parameters (mass, position, velocity, acceleration) are stored in a pandas DataFrame, permitting efficient vectorised computations of the combined forces on each object<sup>&ast;</sup>. Bodies are initialised randomly, but adjustments are made to give the overall system zero momentum but nonzero angular momentum. The `check` method ensures that both the overall momentum and the overall energy of the system remains bounded.

&ast; In this implementation, interactions are computed on a per-body basis (a complexity of \(O(n^2)\). For large \(n\) it may be more efficient to model the interactions indirectly: partition the space and compute a *field* in each element \(m\) contributed by each body, and then, how that field impacts each body. This has complexity \(O(mn)+O(mn)\): so if \(m \ll n\), it will be significantly quicker.

Here's the model:

{{ include_snippet("examples/n-body/n_body.py") }}

## Outputs

The main output is the animated image above. The overall momentum and energy of the system is logged in the console at each timestep.# Hello World

This simple example illustrates the basic model structure, how it all fits together and how it's executed by the framework.

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Model Definition and Execution

The framework runs a model via the `neworder.run` function, which takes an instance of a `Model` class. All models contain, at a minimum:

- an instance of a timeline, in this case the `NoTimeline` arbitrary single-step timeline
- an instance of a `neworder.MonteCarlo` engine
- user-defined methods to evolve the state (`step`) and report/postprocess results (`finalise`).

In this example the model doesn't have an explicit discrete timeline, so for models of this type a method is provided to construct an empty timeline (which is a single step of length zero).

!!! note "Timelines"
    In some model specifications, each individual's entire history can be constructed in a single pass, and in this type of situation a null timeline is appropriate. In more complex examples, the timeline could either refer to absolute time, or for "case-based" models (to use MODGEN parlance), the age of a cohort.

In this (rather contrived) example we have a population a who possess a sole boolean "talkative" attribute, which is initially `False`. The model randomly transitions this state, according to a given probability, and then those individuals who have become talkative say hello.

## Inputs

The input data for this model are just:

- the size of the population
- the probability of saying hello

## Implementation

Firstly we create our model class, subclassing `neworder.Model`:

{{ include_snippet("./examples/hello-world/model.py", "class") }}

and provide a constructor that initialises the base class and a DataFrame containing the population:

{{ include_snippet("./examples/hello-world/model.py", "constructor") }}

!!! note "Unique Indexing"
    The `neworder.df.unique_index()` provides a mechanism to guarantee unique indices for DataFrames, even for parallel runs. This allows individuals to be exchanged and tracked between processes without conflicting indices.

The `step` method randomly samples new values for the "talkative" attribute, using the `neworder.MonteCarlo.hazard` method

{{ include_snippet("./examples/hello-world/model.py", "step") }}

and at the end of the timeline the `finalise` method prints greetings from the talkative individuals using the `neworder.log` function, which is preferred to plain `print` statements as the output is annotated with useful context for debugging.

{{ include_snippet("./examples/hello-world/model.py", "finalise") }}

## Execution

The model is run by simply constructing an instance of our model and passing it to the `run` method:

{{ include_snippet("./examples/hello-world/model.py", "script") }}

From the command line, run the model:

```bash
python examples/hello-world/model.py
```

which should result in something like

```text
[py 0/1]  Hello from 0
[py 0/1]  Hello from 3
[py 0/1]  Hello from 4
[py 0/1]  Hello from 6
```

## Output

To get a better idea of what's going on, uncomment the line containing `neworder.verbose()` and rerun the model. You'll get something like

```text
[no 0/1]  neworder 0.0.8 model init: timeline=<NoTimeline stepped=False> mc=<neworder.MonteCarlo seed=874991939>
[no 0/1]  starting model run. start time=nan
[no 0/1]  t=nan(0) HelloWorld.modify(0)
[no 0/1]  defaulted to no-op Model::modify()
[no 0/1]  t=nan(1) HelloWorld.step()
[no 0/1]  defaulted to no-op Model::check()
[no 0/1]  t=nan(1) HelloWorld.check(): ok
[no 0/1]  t=nan(1) HelloWorld.finalise()
[py 0/1]  Hello from 0
[py 0/1]  Hello from 1
[py 0/1]  Hello from 4
[py 0/1]  Hello from 7
[py 0/1]  Hello from 8
[no 0/1]  SUCCESS exec time=0.000996s
```

this output is explained line-by-line below.

!!! note "Annotated Output"
    The `neworder.log` output is prefixed with a source identifier in square brackets, containing the following information for debugging purposes:

      - Source of message: `no` if logged from the framework itself, `py` if logged from python code (via the `neworder.log()` function). The former are only displayed in verbose mode.
      - the process id ('rank' in MPI parlance) and the total number of processes ('size' in MPI parlance) - in serial mode these default to 0/1.

## Understanding the workflow and the output

When using `NoTimeline()` the start time, end time are undefined and the timestep is zero - i.e. a single arbitrary step in which all computation is done.

First we get some information about the package and the model initialisation parameters:

```text
[no 0/1]  neworder 0.0.8 model init: timeline=<NoTimeline stepped=False> mc=<neworder.MonteCarlo seed=874991939>
[no 0/1]  starting model run. start time=nan
```

the next output concerns the modify method which is explained in the [Option](./option.md) example.

```text
[no 0/1]  t=nan(0) HelloWorld.modify(0)
[no 0/1]  defaulted to no-op Model::modify()
```
Note that the time value is NaN (not-a-number) (since absolute time is not important in this implementation), but the index is zero. Now the `step` method is called, which applies a random transition:

```text
[no 0/1]  t=nan(1) HelloWorld.step()
```

followed by the `check` method, which is optional and we haven't implemented:

```text
[no 0/1]  defaulted to no-op Model::check()
[no 0/1]  t=nan(1) HelloWorld.check() [ok]
```
and the time index has updated to 1.

!!! note "Checks"
    Custom data sanity checks can be run after each timestep by overriding the `check` method. The default implementation does nothing. A typical pattern would be to implement checks for debugging a model during development, then disable them entirely to improve performance using `neworder.checked(False)`. See the other examples.

We've now reached the end of our single step "timeline" so the `finalise` method is called:

```text
[no 0/1]  t=nan(1) HelloWorld.finalise()
```

which prints the results:

```text
[py 0/1]  Hello from 0
[py 0/1]  Hello from 3
[py 0/1]  Hello from 4
[py 0/1]  Hello from 6
```

and finally the model reports its status and execution time:

```text
[no 0/1]  SUCCESS exec time=0.001141s
```

## Next steps

Try re-running the model with different input parameters, or changing the seeding strategy (to e.g. `neworder.MonteCarlo.deterministic_independent_stream`) for reproducible results.

Then, check out some or all of the other examples...
# Chapter 1

This example is based on the example in the introductory chapter of [*Microsimulation and Population Dynamics*](../references.md). This is a simple example - a basic cohort model using continuous-time case-based simulation of mortality (only) for a homogeneous population, using a constant mortality hazard rate. In other words, age at death is sampled from an exponential distribution

\[
p(t)=\lambda e^{-\lambda t}
\]

which has a mean, i.e. life expectancy, of \(\mu=1/\lambda\).

The *neworder* implementation is as direct a port of the MODGEN model, as far as possible.

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Input

Here's the code that sets up and runs the model:

{{ include_snippet("./examples/chapter1/model.py") }}

## Implementation

Each individual is an instance of the `Person` class:

{{ include_snippet("./examples/chapter1/person.py", "person") }}

And the `People` model contains an array of `Person`s

{{ include_snippet("./examples/chapter1/person.py", "constructor") }}

The single timestep records each person's time of death

{{ include_snippet("./examples/chapter1/person.py", "step") }}

And when the timeline has reached the end, the `finalise` method compares the mean of the sampled times of death with the expected value:

{{ include_snippet("./examples/chapter1/person.py", "finalise") }}

Then, from the model script, this function is called, displaying the proportion of the cohort that are still alive at 10-year intervals:

{{ include_snippet("./examples/chapter1/person.py", "alive") }}

## Output

Some example output:

```text
[py 0/1] created 100000 individuals
[py 0/1] Life expectancy = 71.60 years (sampling error=0.17 years)
[py 0/1] Age 10 survival rate = 86.9%
[py 0/1] Age 20 survival rate = 75.6%
[py 0/1] Age 30 survival rate = 65.8%
[py 0/1] Age 40 survival rate = 57.3%
[py 0/1] Age 50 survival rate = 49.8%
[py 0/1] Age 60 survival rate = 43.3%
[py 0/1] Age 70 survival rate = 37.5%
[py 0/1] Age 80 survival rate = 32.7%
[py 0/1] Age 90 survival rate = 28.3%
[py 0/1] Age 100 survival rate = 24.7%
```
and clearly a constant mortality rate isn't realistic as we see far more deaths at younger ages, and far less at older ages, than would be expected. The example [mortality](./mortality.md) introduces a model with a time-dependent mortality hazard rate and shows how the framework can very efficiently model this.
# Wolf-sheep predation

Another implementation of a classic agent-based model

![Wolf-sheep](./img/wolf-sheep.gif)

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Implementation

Rather than representing the agents (wolves, sheep, and grass) as objects, as would be typical in packages like [netlogo](https://ccl.northwestern.edu/netlogo/) or [mesa](https://mesa.readthedocs.io/en/stable/), they are represented as individual rows in pandas DataFrames, which permits efficient vectorised operations on them. Grass grows in fixed "cells" which are used to process interactions. The wolves and sheep roam about randomly at a fixed speed: sheep can only eat grass that is fully grown in the cell they currently occupy, and wolves can only eat sheep within the cell they both occupy.

An extension to the original model adds natural selection: new agents inherit their parent's speed with a random "mutation". Faster animals tend to encounter food more frequently, but conversely consume energy more quickly. The graphic displays a histogram of the speed distributions for wolves and sheep.

Here's the implementation:

{{ include_snippet("./examples/wolf-sheep/wolf_sheep.py") }}

Which is run like so:

{{ include_snippet("./examples/wolf-sheep/model.py") }}

## Outputs

The main output (see image above) is an animation depicting the wolves, sheep and grass in the domain, plus graphs of the populations and histograms of the wolf and sheep speeds. Log messages also record when either the wolf or sheep populations die out completely. The model halts when the wolf and sheep populations die out.
# Mortality

We implement the example *The Life Table* from the second chapter of the book *Microsimulation and Population Dynamics* [[3]](#references). It models mortality in a homogeneous population with an age-specific mortality rate.

This example implements the model in two different ways: firstly a discrete case-based microsimulation, and again using a continuous sampling methodology, showcasing how the latter can be much more efficient. Rather than having a class to represent an individual, as would be standard in a MODGEN implementation, individuals are stored in a *pandas* `Dataframe` which provides fast iteration over the population.

![Mortality histogram](./img/mortality_hist_100k_noloop.gif)

{{ include_snippet("./docs/examples/src.md", show_filename=False) }}

## Inputs

The mortality data is derived from the NewETHPOP[[1]](../references.md) project and represents the mortality rate for white British males in one of the London Boroughs. For timeline definition a maximum age of 100 is defined. Note though that the age at death may be higher than this, it's just the cutoff point at which the hazard rate is assumed to be constant.

{{ include_snippet("./examples/mortality/model.py")}}

## Implementation

### Discrete

As per the MODGEN implementation, we step through a case-based timeline and sample deaths using the marginal mortality rate as a (homogeneous) Poisson process, basically:

- each year, sample time of death for alive individuals
- if year is not at the end of the mortality table
    - if death occurs within the year,
        - record age at death and mark individual as dead
    - otherwise
        - increment age by 1 year and resample
- otherwise
    - record age at death and mark individual as dead
- take mean age at death

So the mortality hazard is treated as a series of piecewise homogeneous Poisson processes.

The discrete model is constructed, as usual, by initialising the base class with a timeline of 100 1-year steps and a seed for the Monte-Carlo engine, loads the mortality data, and initialise the population with age=0 and age at death as yet unknown:

{{ include_snippet("./examples/mortality/people.py", "disc_ctor")}}

The `step` method samples deaths according to the age-specific mortality rate (code not shown for brevity), then increments the age of the living by the timestep.

{{ include_snippet("./examples/mortality/people.py", "disc_step")}}

When the end of the timeline is reached the `finalise` method is called, which checks that all individuals are now dead and computes the life expectancy:

{{ include_snippet("./examples/mortality/people.py", "disc_finalise")}}

### Continuous

The second implementation samples the term structure of mortality directly using the Lewis-Shedler [[4]](../references.md) "thinning" algorithm - this approach doesn't even require a timeline as each individual's age at death can be sampled directly, as a nonhomogeneous Poisson process.

The continuous model doesn't require a `max_age` argument, as there is no timeline, but it does need to know the time resolution of the mortality data in order to sample it correctly:

{{ include_snippet("./examples/mortality/people.py", "cont_ctor")}}

The step method samples, in a single pass, all the deaths, as arrival times in a nonhomogeneous Poisson process:

{{ include_snippet("./examples/mortality/people.py", "cont_step")}}

The single check ensures that all sampled values are finite:

{{ include_snippet("./examples/mortality/people.py", "cont_check")}}

And finally the `finalise` method computes the life expectancy

{{ include_snippet("./examples/mortality/people.py", "cont_finalise")}}

## Output

Running the model script will execute both models

```bash
python examples/mortality/model.py
```

with output like this

```test
[py 0/1] Population = 100000
[py 0/1] Discrete model life expectancy = 77.432356, exec time = 1.321355
[py 0/1] Continuous model life expectancy = 77.388072, exec time = 0.161716
```

which illustrates how much more efficient the continuous implementation is (about ten times faster).

The visualisations (see examples source code for details) show an animated histogram of the deaths (above), and a comparison of the age to death distributions from the two implementations:

![Mortality rate comparison](./img/mortality_100k.png)
