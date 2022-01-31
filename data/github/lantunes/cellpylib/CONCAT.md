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

# ReferencesCiting
------

This project has been published in the
`Journal of Open Source Software <https://joss.theoj.org/papers/10.21105/joss.03608>`_.

This project may be cited as:

**Antunes, L. M. (2021). CellPyLib: A Python Library for working with Cellular Automata.
Journal of Open Source Software, 6(67), 3608.**


BibTeX:

.. code-block::

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
Reference
=========

cellpylib.apen
--------------

.. automodule:: cellpylib.apen
   :members:
   :undoc-members:

cellpylib.bien
--------------

.. automodule:: cellpylib.bien
   :members:
   :undoc-members:

cellpylib.ca_functions
----------------------

.. automodule:: cellpylib.ca_functions
   :members:
   :undoc-members:

cellpylib.ca_functions2d
------------------------

.. automodule:: cellpylib.ca_functions2d
   :members:
   :undoc-members:

cellpylib.ctrbl_rule
--------------------

.. automodule:: cellpylib.ctrbl_rule
  :members:
  :undoc-members:

cellpylib.entropy
-----------------

.. automodule:: cellpylib.entropy
   :members:
   :undoc-members:

cellpylib.hopfield_net
----------------------

.. automodule:: cellpylib.hopfield_net
   :members:
   :undoc-members:

cellpylib.langtons_loop
-----------------------

.. automodule:: cellpylib.langtons_loop
   :members:
   :undoc-members:

cellpylib.sdsr_loop
-----------------------

.. automodule:: cellpylib.sdsr_loop
   :members:
   :undoc-members:

cellpylib.evoloop
-----------------------

.. automodule:: cellpylib.evoloop
   :members:
   :undoc-members:

cellpylib.rule_tables
---------------------

.. automodule:: cellpylib.rule_tables
   :members:
   :undoc-members:

cellpylib.sandpile
-----------------------

.. automodule:: cellpylib.sandpile
   :members:
   :undoc-members:
Fredkin's Self-Replicating CA
-----------------------------

Ed Fredkin described an interesting cellular automaton that exhibits self-replication. The CA is 2-dimensional, and can
consist of two or more colors. To compute the state of a cell at the next timestep, one sums the states of the
neighbouring cells, modulo `p`, where `p` represents the number of colors. The neighborhood can be either of the Moore
or von Neumann type. As the CA evolves, copies of the initial configuration will be produced. The examples below of
these CA are based on John D. Cook's blog posts,
`here <https://www.johndcook.com/blog/2021/05/03/self-reproducing-cellular-automata/>`_ and
`here <https://www.johndcook.com/blog/2021/05/03/multicolor-reproducing-ca/>`_:

The following is an example of the Fredkin self-replicating CA with a von Neumann neighborhood, implemented with
CellPyLib:

.. code-block::

    import cellpylib as cpl
    import numpy as np

    cellular_automaton = cpl.init_simple2d(60, 60)
    # the letter "E"
    cellular_automaton[0][28][28] = 1
    cellular_automaton[0][28][29] = 1
    cellular_automaton[0][28][30] = 1
    cellular_automaton[0][29][28] = 1
    cellular_automaton[0][30][28] = 1
    cellular_automaton[0][30][29] = 1
    cellular_automaton[0][30][30] = 1
    cellular_automaton[0][31][28] = 1
    cellular_automaton[0][32][28] = 1
    cellular_automaton[0][32][29] = 1
    cellular_automaton[0][32][30] = 1

    def activity_rule(n, c, t):
        current_activity = n[1][1]
        return (np.sum(n) - current_activity) % 2

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=20,
                                      apply_rule=activity_rule, neighbourhood="von Neumann")

    cpl.plot2d_animate(cellular_automaton, interval=350)

.. image:: _static/fredkin_von_neumann_demo.gif
    :width: 350

The following is an example of the Fredkin self-replicating CA with a Moore neighborhood, implemented with CellPyLib:

.. code-block::

    import cellpylib as cpl
    import numpy as np

    cellular_automaton = cpl.init_simple2d(60, 60)
    # the letter "E"
    cellular_automaton[0][28][28] = 1
    cellular_automaton[0][28][29] = 1
    cellular_automaton[0][28][30] = 1
    cellular_automaton[0][29][28] = 1
    cellular_automaton[0][30][28] = 1
    cellular_automaton[0][30][29] = 1
    cellular_automaton[0][30][30] = 1
    cellular_automaton[0][31][28] = 1
    cellular_automaton[0][32][28] = 1
    cellular_automaton[0][32][29] = 1
    cellular_automaton[0][32][30] = 1

    def activity_rule(n, c, t):
        current_activity = n[1][1]
        return (np.sum(n) - current_activity) % 2

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=20,
                                      apply_rule=activity_rule, neighbourhood="Moore")

    cpl.plot2d_animate(cellular_automaton, interval=350)

.. image:: _static/fredkin_moore_demo.gif
    :width: 350

The following is an example of the Fredkin self-replicating multi-color CA with a von Neumann neighborhood, implemented
with CellPyLib:

.. code-block::

    import cellpylib as cpl
    import numpy as np

    cellular_automaton = cpl.init_simple2d(60, 60)
    # the letter "E"
    cellular_automaton[0][28][28] = 0
    cellular_automaton[0][28][29] = 1
    cellular_automaton[0][28][30] = 2
    cellular_automaton[0][29][28] = 3
    cellular_automaton[0][30][28] = 4
    cellular_automaton[0][30][29] = 5
    cellular_automaton[0][30][30] = 6
    cellular_automaton[0][31][28] = 7
    cellular_automaton[0][32][28] = 8
    cellular_automaton[0][32][29] = 9
    cellular_automaton[0][32][30] = 10

    def activity_rule(n, c, t):
        current_activity = n[1][1]
        return (np.sum(n) - current_activity) % 11

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=23,
                                      apply_rule=activity_rule, neighbourhood="von Neumann")

    cpl.plot2d_animate(cellular_automaton, interval=350, colormap="viridis")

.. image:: _static/fredkin_multicolor_demo.gif
    :width: 350

**References:**

*Edwin R. Banks, Information Processing and Transmission in Cellular Automata. MIT dissertation. January 1971.*

https://www.johndcook.com/blog/2021/05/03/self-reproducing-cellular-automata/

https://www.johndcook.com/blog/2021/05/03/multicolor-reproducing-ca/
Working with Cellular Automata
------------------------------

A First Example
~~~~~~~~~~~~~~~

The following example illustrates the evolution of the Rule 30 CA, described in `A New Kind of Science`, as implemented
with CellPyLib:

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_simple(200)

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100, memoize=True,
                                    apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))


The initial conditions are instantiated using the function :py:func:`~cellpylib.ca_functions.init_simple`, which, in
this example, creates a 200-dimensional vector consisting of zeroes, except for the component in the center of the
vector, which is initialized with a value of 1. Next, the system is subjected to evolution by calling the
:py:func:`~cellpylib.ca_functions.evolve` function. The system evolves under the rule specified through the
``apply_rule`` parameter. Any function that accepts the three arguments ``n``, ``c`` and ``t`` can be supplied as a
rule, but in this case the built-in function :py:func:`~cellpylib.ca_functions.nks_rule` is invoked to provide Rule 30.
The CA is evolved for 100 ``timesteps``, or 100 applications of the rule to the initial and subsequent conditions.

During each timestep, the function supplied to ``apply_rule`` is invoked for each cell. The ``n`` argument refers to the
neighbourhood of the current cell, and consists of an array (in the 1-dimensional CA case) of the activities (i.e.
states) of the cells comprising the current cell's neighbourhood (an array with length 3, in the case of a 1-dimensional
CA with radius of 1). The ``c`` argument refers to index of the cell under consideration. It serves as a label
identifying the current cell. The ``t`` argument is an integer specifying the current timestep.

Finally, to visualize the results, the :py:func:`~cellpylib.ca_functions.plot` function can be utilized:

.. code-block::

    cpl.plot(cellular_automaton)

.. image:: _static/rule30.png
    :width: 400


How CA are represented
~~~~~~~~~~~~~~~~~~~~~~

In CellPyLib, a CA is an array containing the states of the system at each timestep. The initial state is found at index
0 in the array (and also represents the first timestep), the result of the second timestep is at index 1 in the array,
and so on. A state is represented as an array for 1-dimensional CA, and as an array of arrays for 2-dimensional CA.

.. code-block::

    # An example of a 1D binary CA with 10 cells evolved for 3 timesteps
    [
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],  # 1st timestep (initial conditions)
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0],  # 2st timestep
        [0, 0, 0, 1, 0, 1, 0, 1, 0, 0]   # 3rd timestep
    ]

Initializing CA
~~~~~~~~~~~~~~~~~

A CA is initialized by simply instantiating an array of an array (for 1-dimensional CA), or an array of an array of an
array (for 2-dimensional CA). This will represent the initial conditions of the system, which can be submitted to the
:py:func:`~cellpylib.ca_functions.evolve` function. For convenience, there are several built-in functions for common CA
initializations. For example, the :py:func:`~cellpylib.ca_functions.init_simple` function can be used to create a 1D
binary initialization with a 1 in the center:

.. code-block::

    import cellpylib as cpl
    cellular_automaton = cpl.init_simple(10)
    print(cellular_automaton)
    # [[0 0 0 0 0 1 0 0 0 0]]

An analogous function, :py:func:`~cellpylib.ca_functions2d.init_simple2d`, exists for 2-dimensional CA.

There are built-in functions for initializing CA randomly as well, in the :py:func:`~cellpylib.ca_functions.init_random`
and :py:func:`~cellpylib.ca_functions2d.init_random2d` functions.

Evolving CA
~~~~~~~~~~~~~

CA are evolved with the :py:func:`~cellpylib.ca_functions.evolve` function (for 1-dimensional CA) and the
:py:func:`~cellpylib.ca_functions2d.evolve2d` function (for 2-dimensional CA). The
:py:func:`~cellpylib.ca_functions.evolve` function requires 4 parameters: ``cellular_automaton``, ``timesteps``,
``apply_rule`` and ``r``.

The ``cellular_automaton`` parameter represents the CA consisting of initial conditions. For example, for a 1D CA, a
valid argument could be `[[0,0,0,0,1,0,0,0,0]]`. The initial conditions can include a history of previous states. Thus,
if the length of the array is greater than 1, then the last item in the array will be used as the initial conditions for
the current evolution, and the final CA will include the history supplied. For example, for a 1D CA, a valid argument
that includes a history of previous states could be `[[0,0,0,0,0,0,0,0,0], [0,0,0,0,1,0,0,0,0]]`, and
`[0,0,0,0,1,0,0,0,0]` would be used as the initial state for the evolution.

The ``timesteps`` parameter is simply an integer representing the number of timesteps the CA should undergo evolution, or
application of the supplied rule. Note that the initial conditions of the CA are considered the 1st timestep, so, for
example, if ``timesteps`` is set to `3`, then the rule will be applied two times. This assumes that the number of
timesteps is known in advance. However, in some cases, the number of timesteps may not be known in advance, and the CA
is meant to be evolved until a certain condition is met. For such scenarios, the ``timesteps`` parameter may alternatively
be a callable that accepts the states of the CA over the course of its evolution and the current timestep number, and is
expected to return a boolean indicating whether evolution should continue. If the callable returns `False`, then
evolution is halted.

The ``apply_rule`` parameter expects a callable that represents the rule that will be applied to each cell of the CA at
each timestep. Any kind of callable is valid, but the callable must accept 3 arguments: ``n``, ``c`` and ``t``.
Furthermore, the callable must return the state of the current cell at the next timestep. The `n` argument is the
neighbourhood, which is a NumPy array of length `2r + 1` representing the state of the neighbourhood of the cell (for
1D CA), where ``r`` is the neighbourhood radius. The state of the current cell will always be located at the "center" of
the neighbourhood. The ``c`` argument is the cell identity, which is a scalar representing the index of the cell in the
cellular automaton array. Finally, the ``t`` argument is an integer representing the time step in the evolution.

The ``r`` parameter is simply an integer that represents the radius of the neighbourhood of the CA. For 1D CA, a radius
of 1 implies a neighbourhood is of size 3, a radius of 2 implies a neighbourhood of size 5, and so on. For 2D CA, the
same idea applies, but the neighbourhood will have the dimensions `2r+1 x 2r+1` (with a slight adjustment for von
Neumann neighbourhoods).

Visualizing CA
~~~~~~~~~~~~~~

There are a number of built-in functions to help visualize CA. The simplest is perhaps the
:py:func:`~cellpylib.ca_functions.plot` function, which plots the evolution of a 1D CA. There is also the
:py:func:`~cellpylib.ca_functions.plot_multiple` function, which will create plots for multiple CA in the same
invocation. For 2D CA, there is the :py:func:`~cellpylib.ca_functions2d.plot2d` function. This function accepts an
additional argument, ``timestep``, which represents the particular timestep to be plotted. If none is given, then the
state at the last timestep will be plotted. Finally, the evolution of 2D CA can be animated, with the
:py:func:`~cellpylib.ca_functions2d.plot2d_animate` function.

Increasing Execution Speed with Memoization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Memoization is a means by which computer programs can be made to execute faster. It involves caching the result of a
function for a given input. CellPyLib supports the memoization of rules supplied to the
:py:func:`~cellpylib.ca_functions.evolve` and :py:func:`~cellpylib.ca_functions2d.evolve2d` functions. By default,
memoization is not enabled, since only rules that do not depend on the cell index value, the timestep number, or that
do not store any state as a result of invoking the rule, are supported for memoization. Only the cell neighbourhood is
used to index the output of the rule. Memoization must be explicitly enabled by passing along the ``memoize`` parameter
with a value of `True` when invoking the :py:func:`~cellpylib.ca_functions.evolve` and
:py:func:`~cellpylib.ca_functions2d.evolve2d` functions.

Memoization is expected to provide an increase to execution speed when there is some overhead involved when invoking
the rule. Again, only stateless rules that depend only on the cell neighbourhood are supported. Consider the following
example of rule 30, where memoization is enabled in one case:

.. code-block::

    import cellpylib as cpl
    import time

    start = time.time()
    cpl.evolve(cpl.init_simple(600), timesteps=300,
               apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))
    print(f"Without memoization: {time.time() - start:.2f} seconds elapsed")

    start = time.time()
    cpl.evolve(cpl.init_simple(600), timesteps=300,
               apply_rule=lambda n, c, t: cpl.nks_rule(n, 30), memoize=True)
    print(f"With memoization: {time.time() - start:.2f} seconds elapsed")

The program above prints:

.. code-block::

    Without memoization: 8.23 seconds elapsed
    With memoization: 0.12 seconds elapsed

*(results may differ, depending on the device used)*

**Using Binary Trees and Quadtrees to Exploit Regularities**

To provide a further speed improvement, memoization can be combined with recursive structures, such as binary trees and
quadtrees. This approach is heavily inspired by the well-known
`HashLife algorithm <https://en.wikipedia.org/wiki/Hashlife>`_ for fast cellular automata simulation, which also
combines memoization with quadtrees.

Although CellPyLib does not provide an implementation of HashLife, it does provide an option to sub-divide a finite grid
into halves or quadrants, and apply memoization at various levels of the resulting binary tree or quadtree. This results
in a significant speed-up when there are repetitive and regular patterns in the CA. To combine memoization with tree
structures in CellPyLib, provide the ``memoize`` option with a value of `"recursive"` when calling the
:py:func:`~cellpylib.ca_functions.evolve` and :py:func:`~cellpylib.ca_functions2d.evolve2d` functions. The following
code snippets provide an example:

.. code-block::

    import cellpylib as cpl

    cpl.evolve(cpl.init_simple(600), timesteps=300,
               apply_rule=lambda n, c, t: cpl.nks_rule(n, 30), memoize="recursive")

And for the 2D case:

.. code-block::

    import cellpylib as cpl

    cpl.evolve2d(ca, timesteps=60, neighbourhood="Moore",
                 apply_rule=cpl.game_of_life_rule, memoize="recursive")

Note that, as is the case when using regular memoization (by providing the ``memoize`` option with `True`), only
stateless rules that do not depend on the cell index or timestep number are supported when supplying the ``memoize``
option with a value of `"recursive"`. Also, only CA that exhibit regular and repetitive patterns will demonstrate a
significant decrease in execution times. Finally, 1D CA that have :math:`2^k` cells, or 2D CA that have
:math:`2^k \times 2^k` cells, should result in lower running times when `"recursive"` is used.

To illustrate the operation of this algorithm, consider the following diagram, in which a 1D Elementary CA (radius 1) is
organized into a binary tree:

.. image:: _static/memoize-recursive-1D.png
    :width: 475
    :align: center

There are 8 cells, which are annotated with the letters `a` through `h`. Each node in the binary tree represents the
neighbourhood of one or more contiguous cells. At the bottom of the tree, the leaves are simply the neighbourhoods of
each cell. For example, one node represents `h-a-b`, the neighbourhood for cell `a`. At progressively higher levels of
the tree, each node represents wider neighbourhoods encompassing more cells. The list on the right of the tree
represents the memoization mapping from each neighbourhood to the subsequent values of the cells in question. In this
example, the memoization cache would contain an entry indexed by a state of the neighbourhood `h-a-b` and its associated
value for `a` in the next timestep, as given by the CA rule being used. The idea is to start, at each timestep, from the
root of the tree, looking in the cache for existing states that correspond to the neighbourhoods in the tree, and
updating the values of the cells represented by a node/neighbourhood if an entry exists in the cache. This should result
in an increase in execution speed, since the children of a node needn't be visited if a state was found in the cache.
CA that exhibit regular and repetitive patterns will benefit the most from this approach. CA that exhibit much less
regularity (e.g. ECA Rule 30), will not benefit from this approach, and may in fact incur a performance penalty. In
such a case, it might be best to use a regular memoization scheme (by providing the ``memoize`` option with `True`). For
2D CA, the same concept applies, with the difference that the cells are divided at each level into quandrants rather
than halves, forming a quadtree.

The following table illustrates the running times for 1D CA using various ``memoize`` options:

.. list-table:: Comparison of running times (in seconds) for 1D ECA (1000 timesteps)
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - ECA Rule #
     - `memoize=True`
     - `memoize="recursive"`
     - `memoize=False`
     - # cells
   * - 30
     - 0.74
     - 2.45
     - 47.54
     - 1024
   * - 4
     - 0.63
     - 0.21
     - 44.73
     - 1024
   * - 2
     - 0.81
     - 0.29
     - 46.37
     - 1024
   * - 110
     - 0.66
     - 0.71
     - 50.98
     - 1024
   * - 30
     - 0.73
     - 2.75
     - \-
     - 1000
   * - 4
     - 0.65
     - 0.19
     - \-
     - 1000
   * - 2
     - 0.70
     - 0.28
     - \-
     - 1000
   * - 110
     - 0.64
     - 0.93
     - \-
     - 1000

The following table illustrates the running times for 2D CA using various ``memoize`` options:

.. list-table:: Comparison of running times (in seconds) for 2D CA (700 timesteps)
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - CA
     - `memoize=True`
     - `memoize="recursive"`
     - `memoize=False`
     - # cells
   * - SDSR Loop
     - 117.89
     - 13.86
     - 220.54
     - 100x100
   * - SDSR Loop
     - 165.74
     - 15.76
     - 371.05
     - 128x128
   * - SDSR Loop
     - \-
     - 36.30
     - \-
     - 256x256
   * - Game of Life
     - 39.39
     - 1.74
     - 64.90
     - 60x60
   * - Game of Life
     - 43.62
     - 1.41
     - 68.72
     - 64x64
   * - Game of Life
     - \-
     - 5.14
     - \-
     - 128x128
Varying the Number of Colors
----------------------------

The number of states, or colors, that a cell can adopt is given by `k`. For example, a binary cellular automaton, in
which a cell can assume only values of 0 and 1, has `k = 2`. CellPyLib supports any value of `k`. A built-in function,
:py:func:`~cellpylib.ca_functions.totalistic_rule`, is an implementation of the Totalistic cellular automaton rule, as
described in Stephen Wolfram's `A New Kind of Science`. The code snippet below illustrates using this rule. A value of
`k` of 3 is used, but any value between (and including) 2 and 36 is currently supported. The rule number is given in
base 10 but is interpreted as the rule in base `k` (thus rule 777 corresponds to '1001210' when `k = 3`).

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_simple(200)

    # evolve the CA, using totalistic rule 777 for a 3-color CA
    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100,
                                    apply_rule=lambda n, c, t: cpl.totalistic_rule(n, k=3, rule=777))

    cpl.plot(cellular_automaton)

.. image:: _static/tot3_rule777.png
    :width: 400

Alternatively, the :py:class:`~cellpylib.ca_functions.TotalisticRule` class can be used:

.. code-block::

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100,
                                    apply_rule=cpl.TotalisticRule(k=3, rule=777))

**References:**

*Wolfram, S. (2002). A New Kind of Science. Champaign, IL: Wolfram Media.*
The Collatz Conjecture
----------------------

The Collatz conjecture states that by iteratively applying a particular rule to successive numbers, beginning from any
number, the result will eventually be `1`.

Below is an example of a rule that demonstrates the Collatz conjecture, and also demonstrates the use of a callable for
the ``timesteps`` argument of the :py:func:`~cellpylib.ca_functions.evolve` function, since, in principle, it isn't
known how many iterations are required before the system evolves to a state consisting of the value `1`.

.. code-block::

    import cellpylib as cpl
    import numpy as np

    initial = np.array([[17]], dtype=np.int)

    def activity_rule(n, c, t):
        n = n[1]
        if n % 2 == 0:
            # number is even
            return n / 2
        else:
            return 3*n + 1

    cellular_automaton = cpl.evolve(initial, apply_rule=activity_rule,
                                    timesteps=lambda ca, t: True if ca[-1][0] != 1 else False)

    print([i[0] for i in cellular_automaton])

The program above should print:

.. code-block::

    [17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1]

**References:**

https://en.wikipedia.org/wiki/Collatz_conjecture
Conway's Game of Life
---------------------

Conway's Game of Life is a very famous 2D Cellular Automaton. It uses a simple rule to give rise to a complex system
that is capable of universal computation, in addition to its ability to entertain and fascinate.

CellPyLib has a built-in function, :py:func:`~cellpylib.ca_functions2d.game_of_life_rule`, that can be used to produce
the Game of Life 2D CA:

.. code-block::

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

.. image:: _static/game_of_life.gif
    :width: 350

**References:**

*Conway, J. (1970). The game of life. Scientific American, 223(4), 4.*
Wireworld
---------

Wireworld is a Turing-complete Cellular Automaton, first described by Brian Silverman in 1987. Wireworld can be used to
simulate electronic gates, or logic elements.

An example of Wireworld diodes, implemented with CellPyLib, is given below:

.. code-block::

    import cellpylib as cpl
    import numpy as np
    from matplotlib.colors import ListedColormap


    def wireworld_rule(n, c, t):
        current_activity = n[1][1]
        if current_activity == 0:  # empty
            return 0
        if current_activity == 1:  # electron head
            return 2
        if current_activity == 2:  # electron tail
            return 3
        if current_activity == 3:  # conductor
            electron_head_count = np.count_nonzero(n == 1)
            return 1 if electron_head_count == 1 or electron_head_count == 2 else 3


    cellular_automata = np.array([[
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0],
        [2, 1, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3],
        [0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0],
        [2, 1, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3],
        [0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]])

    cellular_automata = cpl.evolve2d(cellular_automata, timesteps=15,
                                     apply_rule=wireworld_rule, neighbourhood="Moore")

    cpl.plot2d_animate(cellular_automata, show_grid=True, show_margin=False, scale=0.3,
                       colormap=ListedColormap(["black", "blue", "red", "yellow"]))

.. image:: _static/wireworld_diodes.gif
    :width: 300

An example of a Wireworld XOR gate, implemented with CellPyLib, gate is given below:

.. code-block::

    import cellpylib as cpl
    import numpy as np
    from matplotlib.colors import ListedColormap


    def wireworld_rule(n, c, t):
        current_activity = n[1][1]
        if current_activity == 0:  # empty
            return 0
        if current_activity == 1:  # electron head
            return 2
        if current_activity == 2:  # electron tail
            return 3
        if current_activity == 3:  # conductor
            electron_head_count = np.count_nonzero(n == 1)
            return 1 if electron_head_count == 1 or electron_head_count == 2 else 3


    cellular_automata = np.array([[
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 3, 1, 2, 3, 3, 3, 3, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 3, 3, 3, 3, 2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0],
        [0, 0, 0, 3, 3, 2, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0],
        [0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]])

    cellular_automata = cpl.evolve2d(cellular_automata, timesteps=25,
                                     apply_rule=wireworld_rule, neighbourhood="Moore")

    cpl.plot2d_animate(cellular_automata, show_grid=True, show_margin=False, scale=0.3,
                       colormap=ListedColormap(["black", "blue", "red", "yellow"]))

.. image:: _static/wireworld_xor.gif
    :width: 400

**References:**

https://en.wikipedia.org/wiki/Wireworld

*Dewdney, A K (January 1990). "Computer recreations: The cellular automata programs that create Wireworld, Rugworld and
other diversions". Scientific American. 262 (1): 146–149.*
Two-Dimensional CA
------------------

CellPyLib supports 2-dimensional cellular automata with periodic boundary conditions. The number of states, `k`, can be
any whole number. The neighbourhood radius, `r`, can also be any whole number, and both Moore and von Neumann
neighbourhood types are supported. The following snippet demonstrates creating a 2D totalistic cellular automaton:

.. code-block::

    import cellpylib as cpl

    # initialize a 60x60 2D cellular automaton
    cellular_automaton = cpl.init_simple2d(60, 60)

    # evolve the cellular automaton for 30 time steps,
    #  applying totalistic rule 126 to each cell with a Moore neighbourhood
    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=30, neighbourhood='Moore',
                                      apply_rule=lambda n, c, t: cpl.totalistic_rule(n, k=2, rule=126))

    cpl.plot2d(cellular_automaton)

.. image:: _static/tot_rule126_2d_moore.png
    :width: 250

The image above represents the state at the final timestep. However, the state of the CA at any timestep can be
visualized using the :py:class:`~cellpylib.ca_functions2d.plot2d` ``timestep`` argument. For example, in the code
snippet below, the state at the 10th timestep is plotted:

.. code-block::

    cpl.plot2d(cellular_automaton, timestep=10)

.. image:: _static/tot_rule126_2d_moore_t10.png
    :width: 255

Note that 2D CA can also be animated, so that the entire evolution of the CA can be visualized, using the
:py:class:`~cellpylib.ca_functions2d.plot2d_animate` function:

.. code-block::

    cpl.plot2d_animate(cellular_automaton)

.. image:: _static/tot_rule126_2d_moore.gif
    :width: 350
Hopfield Network
----------------

The Hopfield Network can be thought of as a cellular automaton where all cells are neighbours of eachother. The cells
(or neurons) are binary units, and the activity rule is a simple threshold rule, where the weighted inputs to a cell are
summed and compared to a threshold value. The weights are learned from the training data.

CellPyLib includes a built-in implementation of a Hopfield Network, in the
:py:class:`~cellpylib.hopfield_net.HopfieldNet` class, based on the idea that the Hopfield Network can be viewed as a
kind of cellular automaton.

To use it, we must first train the network, by giving it a set of patterns:

.. code-block::

    import cellpylib as cpl
    import numpy as np

    hopfield_net = cpl.HopfieldNet(num_cells=35)

    zero = [
        0, 1, 1, 1, 0,
        1, 0, 0, 0, 1,
        1, 0, 0, 0, 1,
        1, 0, 0, 0, 1,
        1, 0, 0, 0, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0]
    one = [
        0, 1, 1, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 0, 0]
    two = [
        1, 1, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 1, 0,
        0, 1, 1, 0, 0,
        1, 0, 0, 0, 0,
        1, 1, 1, 1, 1,
        0, 0, 0, 0, 0]
    # replace the zeroes with -1 to make these vectors bipolar instead of binary
    one = [-1 if x == 0 else x for x in one]
    two = [-1 if x == 0 else x for x in two]
    zero = [-1 if x == 0 else x for x in zero]

    P = [zero, one, two]

    hopfield_net.train(P)

As shown above, we must instantiate an instance of a :py:class:`~cellpylib.hopfield_net.HopfieldNet`, specifying the
number of cells. Then, we must call :py:func:`~cellpylib.hopfield_net.HopfieldNet.train`, providing a list of training
examples. *NOTE: Only Hopfield Networks with an odd number of cells is currently supported due to limitations of the
implementation.*

Using the :py:class:`~cellpylib.hopfield_net.HopfieldNet` involves providing a potentially incomplete pattern, and
evolving the network for a pre-specified number of timesteps. The network state should settle into a pattern that
resembles those seen during training. It acts like a content-addressable (associative) memory.

.. code-block::

    half_two = [
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 1, 1, 0, 0,
        1, 0, 0, 0, 0,
        1, 1, 1, 1, 1,
        0, 0, 0, 0, 0]
    half_two = [-1 if x == 0 else x for x in half_two]

    cellular_automaton = np.array([half_two])

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=155,
                                    apply_rule=hopfield_net.apply_rule, r=hopfield_net.r)

    cpl.plot(hopfield_net.W)
    cpl.plot2d_animate(np.reshape(cellular_automaton, (155, 7, 5)))

.. image:: _static/hopfield_net_weights.png
    :width: 300

.. image:: _static/hopfield_net.gif
    :width: 400

**References:**

*J. J. Hopfield, "Neural networks and physical systems with emergent collective computational abilities",
Proceedings of the National Academy of Sciences of the USA, vol. 79 no. 8 pp. 2554–2558, April 1982.*

https://en.wikipedia.org/wiki/Hopfield_network

http://neupy.com/2015/09/20/discrete_hopfield_network.html
Additional Features
-------------------

Langton's Lambda
~~~~~~~~~~~~~~~~

One way to specify CA rules is with rule tables. Rule tables are enumerations of all possible neighbourhood states
together with their cell state mappings. For any given neighbourhood state, a rule table provides the associated cell
state value. CellPyLib provides a built-in function for creating random rule tables. The following snippet demonstrates
its usage:

.. code-block::

    import cellpylib as cpl

    rule_table, actual_lambda, quiescent_state = cpl.random_rule_table(lambda_val=0.45, k=4, r=2,
                                                                       strong_quiescence=True,
                                                                       isotropic=True)
    cellular_automaton = cpl.init_random(128, k=4)

    # use the built-in table_rule to use the generated rule table
    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=200,
                                    apply_rule=lambda n, c, t: cpl.table_rule(n, rule_table), r=2)

The following plots demonstrate the effect of varying the lambda parameter:

.. image:: _static/phase_transition.png
    :width: 650

C. G. Langton describes the lambda parameter, and the transition from order to criticality to chaos in cellular
automata while varying the lambda parameter, in the paper:

.. code-block:: text

    Langton, C. G. (1990). Computation at the edge of chaos: phase transitions
    and emergent computation. Physica D: Nonlinear Phenomena, 42(1-3), 12-37.

Reversible CA
~~~~~~~~~~~~~

Elementary CA can be explicitly made to be reversible. CellPyLib has a class,
:py:class:`~cellpylib.ca_functions.ReversibleRule`, which can be used to decorate a rule, making it reversible. The
following example demonstrates the creation of the elementary reversible CA rule 90R:

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_random(200)
    rule = cpl.ReversibleRule(cellular_automaton[0], 90)

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100,
                                    apply_rule=ule)

    cpl.plot(cellular_automaton)

.. image:: _static/rule90R.png
    :width: 400

Asynchronous CA
~~~~~~~~~~~~~~~

Typically, evolving a CA involves the synchronous updating of all the cells in a given timestep. However, it is also
possible to consider CA in which the cells are updated asynchronously. There are various schemes for achieving this.
`Wikipedia <https://en.wikipedia.org/wiki/Asynchronous_cellular_automaton>`_ has a page dedicated to this topic.

CellPyLib has a class, :py:class:`~cellpylib.ca_functions.AsynchronousRule`, which can be used to decorate a rule, making
it asynchronous. In the following example, the rule 60 sequential CA from the notes of `A New Kind of Science` (Chapter
9, section 10: `Sequential cellular automata <http://www.wolframscience.com/nks/notes-9-10--sequential-cellular-automata/>`_)
is implemented:

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_simple(21)

    rule = cpl.AsynchronousRule(apply_rule=lambda n, c, t: cpl.nks_rule(n, 60),
                             update_order=range(1, 20))

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=19*20,
                                    apply_rule=rule)

    # get every 19th row, including the first, as a cycle is completed every 19 rows
    cpl.plot(cellular_automaton[::19])

.. image:: _static/rule60sequential.png
    :width: 300

The :py:class:`~cellpylib.ca_functions.AsynchronousRule` requires the specification of an update order (if none is
provided, then an order is constructed based on the number of cells in the CA). An update order specifies which cell
will be updated as the CA evolves. For example, the update order `[2, 3, 1]` states that cell `2` will be updated in the
next timestep, followed by cell `3` in the subsequent timestep, and then cell `1` in the timestep after that. This
update order is adhered to for the entire evolution of the CA. Cells that are not being updated do not have the
rule applied to them in that timestep.

An option is provided to randomize the update order at the end of each cycle (i.e. timestep). This is equivalent to
selecting a cell randomly at each timestep to update, leaving all others unchanged during that timestep. The following
example demonstrates a simplistic 2D CA in which a cell is picked randomly at each timestep, and its state is updated to
a value of `1` (all cells begin with a state of `0`):

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_simple2d(50, 50, val=0)

    apply_rule = cpl.AsynchronousRule(apply_rule=lambda n, c, t: 1, num_cells=(50, 50),
                                      randomize_each_cycle=True)

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=50,
                                      neighbourhood='Moore', apply_rule=apply_rule)

    cpl.plot2d_animate(cellular_automaton, interval=200, autoscale=True)

.. image:: _static/async_random.gif
    :width: 350


CTRBL Rules
~~~~~~~~~~~

There exists a class of important CA that exhibit the property of self-reproduction. That is, there are patterns
observed that reproduce themselves during the evolution of these CA. This phenomenon has obvious relevance to the study
of Biological systems. These CA are typically 2-dimensional, with a von Neumann neighbourhood of radius 1. The
convention when specifying the rules for these CA is to enumerate the rule table using the states of the center (C), top
(T), right (R), bottom (B), and left (L) cells in the von Neumann neighbourhood.

Such CTRBL CA are supported in CellPyLib, through the :py:class:`~cellpylib.ctrbl_rule.CTRBLRule` class. A particularly
well-known CA in this class is Langton's Loop. CellPyLib has a built-in implementation of this CA, available through
the :py:class:`~cellpylib.langtons_loop.LangtonsLoop` class.

Here is a simple example of a 2D CA that uses a CTRBL rule:

.. code-block::

    import cellpylib as cpl

    ctrbl_rule = cpl.CTRBLRule(rule_table={
        (0, 1, 0, 0, 0): 1,
        (1, 1, 0, 0, 0): 0,
        (0, 0, 0, 0, 0): 0,
        (1, 0, 0, 0, 0): 1,
        (0, 0, 1, 1, 0): 0,
        (1, 1, 1, 1, 1): 1,
        (0, 1, 0, 1, 0): 0,
        (1, 1, 1, 0, 1): 1,
        (1, 0, 1, 0, 1): 1,
        (0, 1, 1, 1, 1): 1,
        (0, 0, 1, 1, 1): 0,
        (1, 1, 0, 0, 1): 1
    }, add_rotations=True)

    cellular_automaton = cpl.init_simple2d(rows=10, cols=10)

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=60,
                                      apply_rule=ctrbl_rule, neighbourhood="von Neumann")

    cpl.plot2d_animate(cellular_automaton)

.. image:: _static/ctrbl.gif
    :width: 400

It is a binary CA that always appears to evolve to some stable attractor state.

Custom Rules
~~~~~~~~~~~~

A rule is a callable that contains the logic that will be applied to each cell of the CA at each timestep. Any kind of
callable is valid, but the callable must accept 3 arguments: ``n``, ``c`` and ``t``. Furthermore, the callable must
return the state of the current cell at the next timestep. The ``n`` argument is the neighbourhood, which is a NumPy
array of length `2r + 1` representing the state of the neighbourhood of the cell (for 1D CA), where ``r`` is the
neighbourhood radius. The state of the current cell will always be located at the "center" of the neighbourhood. The
``c`` argument is the cell identity, which is a scalar representing the index of the cell in the cellular automaton
array. Finally, the ``t`` argument is an integer representing the time step in the evolution.

Any kind of callable is supported, and this is particularly useful if more complex handling, like statefulness, is
required by the rule. For complex rules, the recommended approach is to define a class for the rule, which provides
a ``__call__`` function which accepts the ``n``, ``c``, and ``t`` arguments. The
:py:class:`~cellpylib.ca_functions.BaseRule` class is provided for users to extend, which ensures that the custom rule
is implemented with the correct ``__call__`` signature.

As an example, below is a custom rule that simply keeps track of how many times each cell has been invoked:

.. code-block::

    import cellpylib as cpl
    from collections import defaultdict

    class CustomRule(cpl.BaseRule):

        def __init__(self):
            self.count = defaultdict(int)

        def __call__(self, n, c, t):
            self.count[c] += 1
            return self.count[c]

    rule = CustomRule()

    cellular_automaton = cpl.init_simple(11)

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=10,
                                    apply_rule=rule)

    cpl.plot(cellular_automaton)

.. image:: _static/custom_rule.png
    :width: 250
Installation
------------

CellPyLib can be installed via pip:

.. prompt:: bash $

    pip install cellpylib

Requirements for using this library are Python 3.6, numpy 1.15.4, and matplotlib 3.0.2.Elementary CA
-------------

Elementary CA (ECA) were studied extensively by Stephen Wolfram in his book `A New Kind of Science`. These are perhaps
the simplest kind of CA that one can conceive, with 2 states and a neighbourhood consisting of 3 cells (i.e. a radius of
1). There are a total of 256 ECA (i.e. there are 256 different ways of specifying a rule table for the 8 possible binary
states of a neighbourhood). It is thus possible to exhaustively explore this space of discrete dynamical systems. As
such, it is one of the most studied and well understood type of CA.

CellPyLib supports the creation of ECA through the :py:func:`~cellpylib.ca_functions.nks_rule` function. This function
accepts as a parameter the rule number, using the convention introduced by Stephen Wolfram in his book `A New Kind of
Science`. The rule number uniquely identifies an ECA. For example, Rule 30 is a famous ECA. Its behaviour is very
complex, and it remains poorly understood. Questions regarding its evolution remain unanswered at the time of this
writing, and there is even a `Rule 30 Prize <https://www.rule30prize.org/>`_, offered to those who can answer
fundamental questions about this fascinating dynamical system.

The following code snippet demonstrates creating and visualizing Rule 30 with CellPyLib:

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_simple(200)

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100, memoize=True,
                                    apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))
    cpl.plot(cellular_automaton)

.. image:: _static/rule30.png
    :width: 400

Alternatively, the :py:class:`~cellpylib.ca_functions.NKSRule` class can be used:

.. code-block::

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=100, memoize=True,
                                    apply_rule=cpl.NKSRule(30))

**References:**

*Wolfram, S. (2002). A New Kind of Science. Champaign, IL: Wolfram Media.*
Measures of Complexity
----------------------

CellPyLib provides various built-in functions which can act as measures of complexity in the cellular automata being
examined. These are the information-theoretic properties known as the Shannon entropy and mutual information,
implemented in the :py:func:`~cellpylib.entropy.average_cell_entropy` and
:py:func:`~cellpylib.entropy.average_mutual_information` functions.

Average Cell Entropy
~~~~~~~~~~~~~~~~~~~~

Average cell entropy can reveal something about the presence of information within cellular automata dynamics. The
built-in function :py:func:`~cellpylib.entropy.average_cell_entropy` provides the average Shannon entropy per single
cell in a given cellular automaton. The following snippet demonstrates the calculation of the average cell entropy:

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_random(200)

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=1000,
                                    apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))

    # calculate the average cell entropy; the value will be ~0.999 in this case
    avg_cell_entropy = cpl.average_cell_entropy(cellular_automaton)

The following plots illustrate how average cell entropy changes as a function of Langton's lambda:

.. image:: _static/avg_cell_entropy.png
    :width: 650

Average Mutual Information
~~~~~~~~~~~~~~~~~~~~~~~~~~

The degree to which a cell state is correlated to its state in the next time step can be described using mutual
information. Ideal levels of correlation are required for effective processing of information. The built-in function
:py:func:`~cellpylib.entropy.average_mutual_information` provides the average mutual information between a cell and
itself in the next time step (the temporal distance can be adjusted). The following snippet demonstrates the calculation
of the average mutual information:

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_random(200)

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=1000,
                                    apply_rule=lambda n, c, t: cpl.nks_rule(n, 30))

    # calculate the average mutual information between a cell and itself in the next time step
    avg_mutual_information = cpl.average_mutual_information(cellular_automaton)

The following plots illustrate how average mutual information changes as a function of Langton's lambda:

.. image:: _static/avg_mutual_information.png
    :width: 650

**References**

*Langton, C. G. (1990). Computation at the edge of chaos: phase transitions and emergent computation.
Physica D: Nonlinear Phenomena, 42(1-3), 12-37.*
.. _contents:

CellPyLib
=========

CellPyLib is a Python library for working with Cellular Automata (CA). It provides a concise and simple interface for
defining and analyzing 1- and 2-dimensional CA. The CA can consist of discrete or continuous states. Neighbourhood
radii are adjustable, and in the 2-dimensional case, both Moore and von Neumann neighbourhoods are supported.

With CellPyLib, it is trivial to create Elementary CA, and CA with totalistic rules. These rules are provided as part
of the library. Additionally, the library provides a means for creating asynchronous CA, and reversible CA. Finally, an
implementation of C. G. Langton's approach for creating CA rules using the lambda value is provided, allowing for the
exploration of complex systems, phase transitions and emergent computation.

Utility functions for plotting and viewing the evolved CA are provided. These tools make it easy to visualize the
results of CA evolution. Moreover, utility functions for computing the information-theoretic properties of CA, such as
the Shannon entropy and mutual information, are provided.

.. toctree::
  :caption: Using CellPyLib
  :maxdepth: 5

  installation
  working
  additional
  citing

.. toctree::
  :caption: Tutorials
  :maxdepth: 5

  eca
  neighbourhood
  colors
  complexity
  continuous
  collatz
  twodim
  gol
  wireworld
  fredkin
  hopfield
  langtons_loop
  sdsr_evoloop
  sandpile

.. toctree::
  :caption: API Docs and License
  :maxdepth: 5

  reference
  Source <https://github.com/lantunes/cellpylib>
  license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`Sayama's SDSR Loop and Evoloop
------------------------------

In 1998, as a successor to Langton's Loop, Hiroki Sayama introduced the structurally dissolvable self-reproducing (SDSR)
loop. Structural dissolution represents a form of death for a loop, and the potential for individual loops to replace
others exists. This property introduces an intriguing dynamic behaviour with the potential for evolution.

Below is an example of the SDSR loop. This example makes use of the :py:class:`~cellpylib.sdsr_loop.SDSRLoop` class,
which is an extension of the :py:class:`~cellpylib.langtons_loop.LangtonsLoop` class, which can be used for
constructing any kind of rule based on a von Neumann neighbourhood which considers the Center, Top, Right, Bottom and
Left cells explicitly.

.. code-block::

    import cellpylib as cpl

    sdsr_loop = cpl.SDSRLoop()

    # the initial conditions consist of a single loop
    cellular_automaton = sdsr_loop.init_loops(1, (100, 100), [40], [40])

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=700,
                                      apply_rule=sdsr_loop, memoize="recursive")

    cpl.plot2d_animate(cellular_automaton)

.. image:: _static/sdsr_loop.gif
    :width: 500

While the potential for evolution exists with the SDSR loop, it is the Evoloop which truly exhibits it. Also introduced
by H. Sayama, the Evoloop is an SDSR loop with various phenotypes, that interact and compete with each other for space.

Below is an example of the Evoloop. This example makes use of the :py:class:`~cellpylib.evoloop.Evoloop` class,
which is an extension of the :py:class:`~cellpylib.ctrbl_rule.CTRBLRule` class, which can be used for constructing any
kind of rule based on a von Neumann neighbourhood which considers the Center, Top, Right, Bottom and Left cells
explicitly.

.. code-block::

    import cellpylib as cpl

    evoloop = cpl.Evoloop()

    # the initial conditions consist of a single loop
    cellular_automaton = evoloop.init_species13_loop((100, 100), 40, 15)

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=700,
                                      apply_rule=evoloop, memoize="recursive")

    cpl.plot2d_animate(cellular_automaton)

.. image:: _static/evoloop.gif
    :width: 500

**References**

*Sayama, H. (1998). Constructing evolutionary systems on a simple deterministic cellular automata space.
PhD, University of Tokyo, Department of Information Science.*

*Sayama, H. (1998, August). Introduction of structural dissolution into Langton's self-reproducing loop.
In Proceedings of the sixth international conference on Artificial life (pp. 114-122).*

*Sayama, H. (1999). A new structurally dissolvable self-reproducing loop evolving in a simple cellular automata space.
Artificial Life, 5(4), 343-365.*Sandpiles
---------

A sandpile is a cellular automaton and dynamical system that displays self-organized criticality. It was introduced by
Bak, Tang and Wiesenfeld in 1987.

Below is an example of a sandpile using the built-in :py:class:`~cellpylib.sandpile.Sandpile` class. The boundary of
the 2D CA can be either closed or open. If the boundary is closed, then all boundary cells should have a value of 0.

.. code-block::

    import cellpylib as cpl
    import numpy as np

    n_rows = 45
    n_cols = 45
    sandpile = cpl.Sandpile(n_rows, n_cols)

    ca = np.random.randint(5, size=n_rows*n_cols).reshape((1, n_rows, n_cols))
    # we're using a closed boundary, so make the boundary cells 0
    ca[0, 0, :], ca[0, n_rows-1, :], ca[0, :, 0], ca[0, :, n_cols-1] = 0, 0, 0, 0

    ca = cpl.evolve2d(ca, timesteps=50, apply_rule=sandpile, neighbourhood="von Neumann")

    cpl.plot2d_animate(ca)


.. image:: _static/sandpile.gif
    :width: 500

Note that in the example above, the number of timesteps is fixed at `50`. However, for any random initial conditions, it
isn't obvious how many timesteps are necessary for the system to reach a stable state, or fixed point, in its evolution.
The library provides a function, :py:func:`~cellpylib.ca_functions.until_fixed_point`, that can be called to provide a
callable for the ``timesteps`` argument of the :py:func:`~cellpylib.ca_functions2d.evolve2d` function. By calling this
function, instead of providing a fixed number, the sandpile will evolve until there is no further change in the state of
the system. Below is an example demonstrating this:

.. code-block::

    import cellpylib as cpl
    import numpy as np
    np.random.seed(0)

    n_rows = 45
    n_cols = 45
    sandpile = cpl.Sandpile(n_rows, n_cols)

    initial = np.random.randint(5, size=n_rows*n_cols).reshape((1, n_rows, n_cols))
    # we're using a closed boundary, so make the boundary cells 0
    initial[0, 0, :], initial[0, n_rows-1, :], initial[0, :, 0], initial[0, :, n_cols-1] = 0, 0, 0, 0

    ca = cpl.evolve2d(initial, timesteps=cpl.until_fixed_point(),
                      apply_rule=sandpile, neighbourhood="von Neumann")

    print("Number of timesteps to reach fixed point: %s" % len(ca))
    cpl.plot2d_animate(ca)

The above program will print out the following:

.. code-block::

    Number of timesteps to reach fixed point: 51

If one perturbs a sandpile that has reached a fixed point in its evolution, then the sandpile should reconfigure itself
and eventually reach a fixed point once again. This can be demonstrated using the
:py:func:`~cellpylib.sandpile.Sandpile.add_grain` function of the :py:class:`~cellpylib.sandpile.Sandpile` class. In the
example below, we begin with a sandpile that has reached a fixed point (this configuration is assumed to exist in a file
called `sandpile_add_grain_demo.txt`, located in the same directory as the program). We drop a grain of sand on the cell
located at row with index `23` and column with index `23` at timestep `1`.

.. code-block::

    import cellpylib as cpl
    import numpy as np

    n_rows = 45
    n_cols = 45
    sandpile = cpl.Sandpile(n_rows, n_cols)
    sandpile.add_grain(cell_index=(23, 23), timestep=1)

    initial = np.loadtxt('sandpile_add_grain_demo.txt', dtype=int)
    initial = np.array([initial])

    ca = cpl.evolve2d(initial, timesteps=cpl.until_fixed_point(),
                      apply_rule=sandpile, neighbourhood="von Neumann")

    print("Number of timesteps to reach fixed point: %s" % len(ca))
    cpl.plot2d_animate(ca)


.. image:: _static/sandpile_add_grain.gif
    :width: 500

Imagine dropping a single grain of sand repeatedly on the same cell in the center of a grid, allowing the sandpile to
attain a stable state (i.e. fixed point) before dropping the next grain. This can be demonstrated with the following
code snippet:

.. code-block::

    import cellpylib as cpl

    n = 50
    sandpile = cpl.Sandpile(n, n)
    ca = cpl.init_simple2d(n, n, val=5)

    for i in range(300):
        ca[-1, n//2, n//2] += 1
        ca = cpl.evolve2d(ca, apply_rule=sandpile,
                          timesteps=cpl.until_fixed_point(), neighbourhood='Moore')

    cpl.plot2d_animate(ca)

.. image:: _static/sandpile_growing.gif

Above, we take advantage of the fact that the `ca` argument to :py:func:`~cellpylib.ca_functions2d.evolve2d` can
contain a history of prior states, and that the evolution continues from the last state.

**References:**

*Bak, Per, Chao Tang, and Kurt Wiesenfeld. "Self-organized criticality." Physical review A 38.1 (1988): 364.*

https://en.wikipedia.org/wiki/Abelian_sandpile_model
Continuous CA
-------------

Cellular Automata needn't consist of discrete activities. The units in an automaton can also take on continuous-valued
activities (i.e. states).

The example below implements a continuous-valued Cellular Automaton from Stephen Wolfram's book `A New Kind of Science`,
found on page 157:

.. code-block::

    import math
    import numpy as np
    import cellpylib as cpl

    cellular_automaton = cpl.init_simple(200, dtype=np.float64)

    def apply_rule(n, c, t):
        result = (sum(n) / len(n)) * (3 / 2)
        frac, whole = math.modf(result)
        return frac

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=150,
                                    apply_rule=apply_rule)

    cpl.plot(cellular_automaton)

.. image:: _static/continuous_ca.png
    :width: 400

**References:**

*Wolfram, S. (2002). A New Kind of Science (page 157). Champaign, IL: Wolfram Media.*
Varying the Neighbourhood Size
------------------------------

The size of the cell neighbourhood can be varied by setting the parameter ``r`` when calling the
:py:func:`~cellpylib.ca_functions.evolve` function. The value of ``r`` represents the number of cells to the left and
to the right of the cell under consideration. Thus, to get a neighbourhood size of 3, ``r`` should be 1, and to get a
neighbourhood size of 7, ``r`` should be 3. As an example, consider the work of M. Mitchell et al. involving the
creation (discovery) of a cellular automaton that solves the density classification problem: if the initial random
binary vector contains more than 50% of 1s, then a cellular automaton that solves this problem will give rise to a
vector that contains only 1s after a fixed number of time steps, and likewise for the case of 0s. A very effective
cellular automaton that solves this problem most of the time was found using a Genetic Algorithm.

.. code-block::

    import cellpylib as cpl

    cellular_automaton = cpl.init_random(149)

    # Mitchell et al. discovered this rule using a Genetic Algorithm
    rule_number = 6667021275756174439087127638698866559

    # evolve the CA, setting r to 3, for a neighbourhood size of 7
    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=149,
                                    apply_rule=lambda n, c, t: cpl.binary_rule(n, rule_number),
                                    r=3)

    cpl.plot(cellular_automaton)

.. image:: _static/density_classification.png
    :width: 400

Alternatively, the :py:class:`~cellpylib.ca_functions.BinaryRule` class can be used:

.. code-block::

    cellular_automaton = cpl.evolve(cellular_automaton, timesteps=149,
                                    apply_rule=cpl.BinaryRule(rule_number), r=3)

**References:**

*Melanie Mitchell, James P. Crutchfield, and Rajarshi Das, "Evolving Cellular Automata with Genetic Algorithms:
A Review of Recent Work", In Proceedings of the First International Conference on Evolutionary Computation and Its
Applications (EvCA'96), Russian Academy of Sciences (1996).*Langton's Loops
---------------

In 1984, Christopher Langton described a type of 2-dimensional cellular automaton that exhibits a self-replicating
dynamic loop structure. A branch of Artificial Life research developed from this work, resulting in better insight into
self-replicating processes, which has obvious relevance to Biology and living systems.

Below is an example of Langton's loop. This example makes use of the :py:class:`~cellpylib.langtons_loop.LangtonsLoop`
class, which is an extension of the :py:class:`~cellpylib.ctrbl_rule.CTRBLRule` class, which can be used for
constructing any kind of rule based on a von Neumann neighbourhood which considers the Center, Top, Right, Bottom and
Left cells explicitly.

.. code-block::

    import cellpylib as cpl

    langtons_loop = cpl.LangtonsLoop()

    # the initial conditions consist of a single loop
    cellular_automaton = langtons_loop.init_loops(1, (75, 75), [40], [25])

    cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=500,
                                      apply_rule=langtons_loop, memoize="recursive")

    cpl.plot2d_animate(cellular_automaton)

.. image:: _static/langtons_loops.gif
    :width: 500

**References:**

*Langton, C. G. (1984). Self-reproduction in Cellular Automata. Physica D: Nonlinear Phenomena, 10(1-2), 135-144.*

https://en.wikipedia.org/wiki/Langton%27s_loops
