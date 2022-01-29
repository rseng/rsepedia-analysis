# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [2.3.1] - 2021-12-25

### Changed

- Changed `plot2d_animate` so that it returns the animation object, to address a problem arising in Spyder IDE

## [2.3.0] - 2021-12-01

### Added

- Added support for `memoize='recursive'` option of `evolve` and `evolve2d` functions
- Added `NKSRule`, `BinaryRule` and `TotalisticRule` classes

## [2.2.0] - 2021-11-30

### Added

- Added SDSR loop and Evoloop implementations
- Added `memoize` option to `evolve` and `evolve2d` functions

## [2.1.0] - 2021-11-16

### Added

- Added more Sandpile demos and more content to the Sandpile tutorial in the docs

### Changed

- Changed interpretation of the `cellular_automaton` argument in `evolve` and `evolve2d` such that a history of states can be provided

## [2.0.0] - 2021-11-10

### Added 

- Added more test coverage
- Added `CHANGELOG.md`
- Added docs and tests for `bits_to_int` and `int_to_bits` functions
- Added more documentation to functions in `entropy.py` and `bien.py`, and to `plot2d_slice` and `plot2d_spacetime`
- Added the `BaseRule` class, which provides a base for custom callable rules
- Added built-in `Sandpile` implementation
- Added `show=True` argument to plotting function signatures
- Added `show_grid`, `show_margin` and `scale` arguments to `plot2d` and `plot2d_slice` functions

### Changed

- Addressing test warnings by making subtle adjustments to the code, such as using `np.int32` instead of `np.int`
- Replaced copyright notice in `README.md` with link to Apache License 2.0
- Importing modules explicitly in `__init__.py` to avoid polluting namespace
- Changed `AsynchronousRule`, `ReversibleRule`, and `CTRBLRule` so that they extend `BaseRule` and implement `__call__`
- Changed plotting function signatures so that they accept `imshow` keyword args
- Changed the `evolve` and `evolve2d` functions so that the `timesteps` parameter can alternatively be a callable,
  so that models where the number of timesteps is not known in advance are supported

## [1.1.0] - 2021-08-02

### Added

- Added support for CTRBL rules
- Added Langton's Loop implementation
- Added Wireworld demo code
- Added more optional arguments to `plot2d_animate` function signature

## [1.0.0] - 2021-07-29

### Added

- Initial stable release
- Added more documentation to code
CellPyLib
=========

CellPyLib is a library for working with Cellular Automata, for Python. Currently, only 1- and 2-dimensional _k_-color 
cellular automata with periodic boundary conditions are supported. The size of the neighbourhood can be adjusted. While
cellular automata constitute a very broad class of models, this library focuses on those that are constrained to a 
regular array or uniform grid, such as elementary CA, and 2D CA with Moore or von Neumann neighbourhoods. The 
cellular automata produced by this library match the corresponding cellular automata available 
at [atlas.wolfram.com](http://atlas.wolfram.com).

[![testing status](https://github.com/lantunes/cellpylib/actions/workflows/python-app.yml/badge.svg?branch=master)](https://github.com/lantunes/cellpylib/actions)
[![latest version](https://img.shields.io/pypi/v/cellpylib?style=flat-square&logo=PyPi&logoColor=white&color=blue)](https://pypi.org/project/cellpylib/)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03608/status.svg)](https://doi.org/10.21105/joss.03608)

Example usage:
```python
import cellpylib as cpl

# initialize a CA with 200 cells (a random initialization is also available) 
cellular_automaton = cpl.init_simple(200)

# evolve the CA for 100 time steps, using Rule 30 as defined in NKS
cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100, memoize=True, 
                                apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))

# plot the resulting CA evolution
cpl.plot(cellular_automaton)

```

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/rule30.png" width="50%"/>

You should use CellPyLib if:
* you are an instructor or student wishing to learn more about Elementary Cellular Automata and 2D Cellular Automata on 
a uniform grid (such as the Game of Life, the Abelian sandpile, Langton's Loops, etc.)
* you are a researcher who wishes to work with Elementary Cellular Automata and/or 2D Cellular Automata on a uniform 
grid, and would like to use a flexible, correct and tested library that provides access to such models as part of your 
research

_If you would like to work with Cellular Automata on arbitrary networks (i.e. non-uniform grids), have a look at
the [Netomaton](https://github.com/lantunes/netomaton) project._

### Getting Started

CellPyLib can be installed via pip:

```
pip install cellpylib
```

Requirements for using this library are Python 3.6, NumPy, and Matplotlib. Have a look at the documentation, located at 
[cellpylib.org](https://cellpylib.org), for more information.

## Varying the Neighbourhood Size

The size of the cell neighbourhood can be varied by setting the parameter _*r*_ when calling the `evolve` function. The
value of _*r*_ represents the number of cells to the left and to the right of the cell under consideration. Thus, to
get a neighbourhood size of 3, _*r*_ should be 1, and to get a neighbourhood size of 7, _*r*_ should be 3.
As an example, consider the work of M. Mitchell et al., carried out in the 1990s, involving the creation (discovery) of
a cellular automaton that solves the density classification problem: if the initial random binary vector contains 
more than 50% of 1s, then a cellular automaton that solves this problem will give rise to a vector that contains only
1s after a fixed number of time steps, and likewise for the case of 0s. A very effective cellular automaton that solves
this problem most of the time was found using a Genetic Algorithm.

```python
import cellpylib as cpl

cellular_automaton = cpl.init_random(149)

# Mitchell et al. discovered this rule using a Genetic Algorithm
rule_number = 6667021275756174439087127638698866559

# evolve the CA, setting r to 3, for a neighbourhood size of 7
cellular_automaton = cpl.evolve(cellular_automaton, timesteps=149,
                                apply_rule=lambda n, c, t: cpl.binary_rule(n, rule_number), r=3)

cpl.plot(cellular_automaton)
```
<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/density_classification.png" width="50%"/>

For more information, see:

> Melanie Mitchell, James P. Crutchfield, and Rajarshi Das, "Evolving Cellular Automata with Genetic Algorithms: A Review of Recent Work", In Proceedings of the First International Conference on Evolutionary Computation and Its Applications (EvCA'96), Russian Academy of Sciences (1996).

## Varying the Number of Colors

The number of states, or colors, that a cell can adopt is given by _k_. For example, a binary cellular automaton, in which a cell can 
assume only values of 0 and 1, has _k_ = 2. CellPyLib supports any value of _k_. A built-in function, `totalistic_rule`,
is an implementation of the [Totalistic cellular automaton rule](http://mathworld.wolfram.com/TotalisticCellularAutomaton.html), 
as described in [Wolfram's NKS](https://www.wolframscience.com/nks/). The code snippet below illustrates using this rule. 
A value of _k_ of 3 is used, but any value between (and including) 2 and 36 is currently supported. The rule number is 
given in base 10 but is interpreted as the rule in base _k_ (thus rule 777 corresponds to '1001210' when _k_ = 3).

```python
import cellpylib as cpl

cellular_automaton = cpl.init_simple(200)

# evolve the CA, using totalistic rule 777 for a 3-color CA
cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100,
                                apply_rule=lambda n, c, t: cpl.totalistic_rule(n, k=3, rule=777))

cpl.plot(cellular_automaton)
```

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/tot3_rule777.png" width="50%"/>

## Rule Tables

One way to specify cellular automata rules is with rule tables. Rule tables are enumerations of all possible 
neighbourhood states together with their cell state mappings. For any given neighbourhood state, a rule table provides 
the associated cell state value. CellPyLib provides a built-in function for creating random rule tables. The following
snippet demonstrates its usage:
```python
import cellpylib as cpl

rule_table, actual_lambda, quiescent_state = cpl.random_rule_table(lambda_val=0.45, k=4, r=2,
                                                                   strong_quiescence=True, isotropic=True)

cellular_automaton = cpl.init_random(128, k=4)

# use the built-in table_rule to use the generated rule table
cellular_automaton = cpl.evolve(cellular_automaton, timesteps=200,
                                apply_rule=lambda n, c, t: cpl.table_rule(n, rule_table), r=2)
```
The following plots demonstrate the effect of varying the lambda parameter:

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/phase_transition.png" width="100%"/>

C. G. Langton describes the lambda parameter, and the transition from order to criticality to chaos in cellular 
automata while varying the lambda parameter, in the paper:

> Langton, C. G. (1990). Computation at the edge of chaos: phase transitions and emergent computation. Physica D: Nonlinear Phenomena, 42(1-3), 12-37.

## Measures of Complexity

CellPyLib provides various built-in functions which can act as measures of complexity in the cellular automata being
examined.

### Average Cell Entropy

Average cell entropy can reveal something about the presence of information within cellular automata dynamics. The 
built-in function `average_cell_entropy` provides the average Shannon entropy per single cell in a given cellular 
automaton. The following snippet demonstrates the calculation of the average cell entropy:

```python
import cellpylib as cpl

cellular_automaton = cpl.init_random(200)

cellular_automaton = cpl.evolve(cellular_automaton, timesteps=1000,
                                apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))

# calculate the average cell entropy; the value will be ~0.999 in this case
avg_cell_entropy = cpl.average_cell_entropy(cellular_automaton)
```

The following plots illustrate how average cell entropy changes as a function of lambda:

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/avg_cell_entropy.png" width="100%"/>

### Average Mutual Information

The degree to which a cell state is correlated to its state in the next time step can be described using mutual 
information. Ideal levels of correlation are required for effective processing of information. The built-in function 
`average_mutual_information` provides the average mutual information between a cell and itself in the next time step 
(the temporal distance can be adjusted). The following snippet demonstrates the calculation of the average mutual 
information:

```python
import cellpylib as cpl

cellular_automaton = cpl.init_random(200)

cellular_automaton = cpl.evolve(cellular_automaton, timesteps=1000,
                                apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))

# calculate the average mutual information between a cell and itself in the next time step
avg_mutual_information = cpl.average_mutual_information(cellular_automaton)
```

The following plots illustrate how average mutual information changes as a function of lambda:

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/avg_mutual_information.png" width="100%"/>

## Reversible Cellular Automata

Elementary cellular automata can be explicitly made to be reversible. The following example demonstrates the 
creation of the elementary reversible cellular automaton rule 90R:
  
```python
import cellpylib as cpl

cellular_automaton = cpl.init_random(200)
rule = cpl.ReversibleRule(cellular_automaton[0], 90)

cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100, 
                                apply_rule=rule)

cpl.plot(cellular_automaton)
```

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/rule90R.png" width="50%"/>

## Continuous Cellular Automata

In addition to discrete values, cellular automata can assume continuous values. CellPyLib supports 
continuous-valued automata. To create cellular automata with continuous values--or any kind of data type--simply 
specify the `dtype` parameter when invoking any of the `init` and `evolve` built-in functions. For example, to create
a cellular automata with continuous values, one might specify the following parameter: `dtype=np.float32`.

## 2D Cellular Automata

CellPyLib supports 2-dimensional cellular automata with periodic boundary conditions. The number of states, _k_, can be
any whole number. The neighbourhood radius, _r_, can also be any whole number, and both Moore and von Neumann 
neighbourhood types are supported. The following snippet demonstrates creating a 2D totalistic cellular automaton:

```python
import cellpylib as cpl

# initialize a 60x60 2D cellular automaton 
cellular_automaton = cpl.init_simple2d(60, 60)

# evolve the cellular automaton for 30 time steps, 
#  applying totalistic rule 126 to each cell with a Moore neighbourhood
cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=30, neighbourhood='Moore',
                                  apply_rule=lambda n, c, t: cpl.totalistic_rule(n, k=2, rule=126))

cpl.plot2d(cellular_automaton)
```

The `plot2d` function plots the state of the cellular automaton at the final time step:

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/tot_rule126_2d_moore.png" width="30%"/>

### Conway's Game of Life

There are a number of built-in plotting functions for 2D cellular automata. For example, `plot2d_animate` will animate 
the evolution of the cellular automaton. This is illustrated in the following snippet, which demonstrates the built-in 
Game of Life rule:

```python
import cellpylib as cpl

# Glider
cellular_automaton = cpl.init_simple2d(60, 60)
cellular_automaton[:, [28,29,30,30], [30,31,29,31]] = 1

# Blinker
cellular_automaton[:, [40,40,40], [15,16,17]] = 1

# Light Weight Space Ship (LWSS)
cellular_automaton[:, [18,18,19,20,21,21,21,21,20], [45,48,44,44,44,45,46,47,48]] = 1

# evolve the cellular automaton for 60 time steps
cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=60, neighbourhood='Moore',
                                  apply_rule=cpl.game_of_life_rule, memoize='recursive')

cpl.plot2d_animate(cellular_automaton)
```

<img src="https://raw.githubusercontent.com/lantunes/cellpylib/master/resources/game_of_life.gif" width="65%"/>

For more information about Conway's Game of Life, see:

> Conway, J. (1970). The game of life. Scientific American, 223(4), 4.

### Increasing Execution Speed with Memoization

Memoization is expected to provide an increase to execution speed when there is some overhead involved when invoking 
the rule. Only stateless rules that depend only on the cell neighbourhood are supported. Consider the following
example of rule 30, where memoization is enabled:

```python
import cellpylib as cpl
import time

start = time.time()
cpl.evolve(cpl.init_simple(1000), timesteps=500,
           apply_rule=lambda n, c, t: cpl.nks_rule(n, 30), memoize=True)

print(f"Elapsed: {time.time() - start:.2f} seconds")
```

The program above prints `Elapsed: 0.33 seconds` (actual execution time may vary, depending on the device used). 
Without memoization, the program requires approximately 23 seconds to complete. 

--------------------

### Development

Create a Conda environment from the provided environment YAML file:
```
$ conda env create -f environment.yml
```

**Documentation**

To build the Sphinx documentation locally, from the `doc` directory:
```
$ make clean html
```
The generated files will be in `_build/html`.

To build the documentation for publication, from the `doc` directory:
```
$ make github
```
The generated files will be in `_build/html` and in the `site/docs` folder.

**Testing**

There are a number of unit tests for this project. To run the tests:
```
$ python -m pytest tests
```

If the `pytest-cov` package is installed, a coverage report can be generated by running the tests with:
```
$ python -m pytest tests/ --cov=cellpylib
```
--------------------

### Support

If you have any questions or comments, or if you find any bugs, please open an issue in this project. Please feel free
to fork the project, and create a pull request, if you have any improvements or bug fixes. We welcome all feedback and
contributions.

--------------------

### Citation Info

This project has been published in the 
[Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.03608).
This project may be cited as:


> Antunes, L. M. (2021). CellPyLib: A Python Library for working with Cellular Automata. 
Journal of Open Source Software, 6(67), 3608.


BibTeX:
```
@article{Antunes2021,
  doi = {10.21105/joss.03608},
  url = {https://doi.org/10.21105/joss.03608},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {67},
  pages = {3608},
  author = {Luis M. Antunes},
  title = {CellPyLib: A Python Library for working with Cellular Automata},
  journal = {Journal of Open Source Software}
}
```

### Stars

Please star this repository if you find it useful, or use it as part of your research.

## License

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
---
title: 'CellPyLib: A Python Library for working with Cellular Automata'
tags:
  - Python
  - Cellular Automata
  - complex systems
  - non-linear dynamics
  - discrete dynamical systems
authors:
  - name: Luis M. Antunes
    orcid: 0000-0002-4867-5635
    affiliation: 1
affiliations:
 - name: Department of Chemistry, University of Reading, Whiteknights, Reading RG6 6DX, United Kingdom
   index: 1
date: 28 July 2021
bibliography: paper.bib
---

# Summary

Cellular Automata (CA) are discrete dynamical systems with a rich history [@ilachinski2001cellular]. Introduced by John 
von Neumann and Stanislaw Ulam in the 1940s [@von1951general], CA have continued to fascinate, as their conceptual 
simplicity serves as a powerful microscope that allows us to explore the nature of computation and complexity, and the 
origins of emergence. Far from being an antiquated computational model, investigators are utilizing CA in novel and 
creative ways, such as the incorporation with Deep Learning [@nichele2017deep; @mordvintsev2020growing]. Popularized 
and investigated by Stephen Wolfram in his book *A New Kind of Science* [@wolfram2002new], CA remain premier reminders 
of a common theme in the study of the physical world: that simple systems and rules can give rise to remarkable 
complexity. They are a laboratory for the study of the origins of the complexity we see in the world around us.

`CellPyLib` is a Python library for working with CA. It provides a concise and simple interface for defining and 
analyzing 1- and 2-dimensional CA. The CA can consist of discrete or continuous states. Neighbourhood radii are 
adjustable, and in the 2-dimensional case, both Moore and von Neumann neighbourhoods are supported. With `CellPyLib`, it 
is trivial to create Elementary CA, and CA with totalistic rules, as these rules are provided as part of the library. 
Additionally, the library provides a means for creating asynchronous and reversible CA. Finally, an implementation 
of C. G. Langton's approach for creating CA rules using the lambda value is provided, allowing for the exploration of 
complex systems, phase transitions and emergent computation [@langton1990computation].

Utility functions for plotting and viewing the evolved CA are also provided. These tools make it easy to visualize the
results of CA evolution, and include the option of creating animations of the evolution itself. Moreover, utility 
functions for computing the information-theoretic properties of CA, such as the Shannon entropy and mutual information, 
are included.

# Statement of need

The Python software ecosystem is lacking when it comes to Cellular Automata. A web search reveals that while there are 
some projects dedicated to the simulation of CA, most are not general-purpose, focusing only on certain CA systems, and 
are generally missing a substantial test suite, hindering their future extensibility and maintainability. In short, 
there appears to be a dearth of robust and flexible libraries for working with CA in Python. 

Currently, many scientists choose Python as their main tool for computational tasks. Though researchers can choose to 
implement CA themselves, this is error-prone, as there are some subtleties when it comes to correctly handling issues 
such as boundary conditions on periodic lattices, or constructing von Neumann neighbourhoods with radius greater than 1, 
for example. Researchers may be dissuaded from incorporating CA into their research if they are forced to work with 
unfamiliar languages and technologies, or are required to devote considerable effort to the implementation and testing 
of non-trivial algorithms. The availability of a user-friendly Python library for CA will likely encourage more 
researchers to consider these fascinating dynamical and computational systems. Moreover, having a standard 
implementation of CA in the Python environment helps to ensure that scientific results are reproducible. CellPyLib is a 
Python library aimed to meet this need, which supports the creation and analysis of models that exist on a regular 
array or uniform grid, such as elementary CA, and 2D CA with Moore or von Neumann neighbourhoods.

Researchers and students alike should find CellPyLib useful. Students and instructors can use CellPyLib in an 
educational context if they would like to learn about elementary CA and 2D CA on a uniform grid. Researchers in both the 
computer and physical sciences can use CellPyLib to answer serious questions about the computational and natural worlds. 
For example, the Abelian sandpile model included in the library can be used as part of a university course on complex 
systems to demonstrate the phenomenon of self-organized criticality. The same model may be used by professional 
physicists wishing to explore self-organized criticality more deeply. 

While CellPyLib is expected to be of interest to students, educators, and researchers, there are certain scenarios in 
which alternate tools would be more appropriate. For example, if a researcher would like to evolve CA with a very large
number of cells, or for very many iterations, in a timely fashion, then an implementation that is optimized for the 
specific model in question would be more appropriate. Also, if the model is not constrained to a uniform grid, then
other solutions should be sought.

# Example Usage

`CellPyLib` can be readily installed using `pip`:

```
$ pip install cellpylib
```

It has minimal dependencies, depending only on the commonly used libraries NumPy [@harris2020array] and Matplotlib 
[@Hunter:2007].

The following example illustrates the evolution of the Rule 30 CA, described in `A New Kind of Science` 
[@wolfram2002new], as implemented with `CellPyLib`:

```python
import cellpylib as cpl

cellular_automaton = cpl.init_simple(200)

rule = lambda n, c, t: cpl.nks_rule(n, 30)

cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100, 
                                apply_rule=rule)
```

First, the initial conditions are instantiated using the function `init_simple`, which, in this example, creates a 
200-dimensional vector consisting of zeroes, except for the component in the center of the vector, which is initialized
with a value of 1. Next, the system is subjected to evolution by calling the `evolve` function. The system evolves under 
the rule specified through the `apply_rule` parameter. Any function that accepts the three arguments `n`, `c` and `t` 
can be supplied as a rule, but in this case the built-in function `nks_rule` is invoked to provide Rule 30. The CA is 
evolved for 100 `timesteps`, or 100 applications of the rule to the initial and subsequent conditions.

During each timestep, the function supplied to `apply_rule` is invoked for each cell. The `n` argument refers to the 
neighbourhood of the current cell, and consists of an array (in the 1-dimensional CA case) of the activities (i.e. 
states) of the cells comprising the current cell's neighbourhood (an array with length 3, in the case of a 1-dimensional 
CA with radius of 1). The `c` argument refers to index of the cell under consideration. It serves as a label identifying 
the current cell. The `t` argument is an integer specifying the current timestep.

Finally, to visualize the results, the `plot` function can be utilized:

```python
cpl.plot(cellular_automaton)
```

![Rule 30, as rendered with CellPyLib.\label{fig:rule30}](rule30.png){ width=60% }

The result is rendered, as depicted in \autoref{fig:rule30}.

# Scope

While `CellPyLib` is a general-purpose library that allows for the implementation of a wide variety of CA, it is 
important to note that CA constitute a very broad class of models. `CellPyLib` focuses on those that are constrained to 
a regular array or uniform grid, such as elementary CA, and 2D CA with Moore or von Neumann neighbourhoods.

# References