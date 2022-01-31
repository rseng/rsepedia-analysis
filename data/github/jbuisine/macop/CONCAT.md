---
title: 'Minimalist And Customisable Optimisation Package'
tags:
  - Python
  - Operations Research
  - Mono-objective
  - Multi-objective
authors:
  - name: Jérôme Buisine
    orcid: 0000-0001-6071-744X
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Samuel Delepoulle
    # orcid: 0000-0002-8897-0858
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Christophe Renaud
    # orcid: 0000-0002-9826-8667
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Univ. Littoral Côte d’Opale, LISIC Calais, France, F-62100
   index: 1
date: 11 October 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Optimisation problems are frequently encountered in science and industry. Given a real-valued function $f$ defined on a set called the search space $X$, optimising the function $f$ consists of finding a point $x \in X$ that has the optimal value $f(x)$, or at least constructing a sequence $(x_t)_{t \in \mathbf{N}} \in X^\mathbb{N}$ that is close to the optimum. Depending on the search space $X$, optimisation problems can be globally classified as discrete problems (e.g. $X=\{0,1\}^n$) or as continuous problems (e.g. $X=\mathbb{R}^n$). Tools for modelling and solving discrete [@solid-solver] and continuous [@ceres-solver] problems are proposed in the literature.

In this paper, `Macop` for `Minimalist And Customisable Optimisation Package`, is proposed as a discrete optimisation Python package which doesn't implement every algorithm in the literature, but provides the ability to quickly develop and test your own algorithm and strategies. The main objective of this package is to provide maximum flexibility, which allows easy implementation when experimenting new algorithms.

Based on a common interaction loop (see \autoref{fig:macop-behaviour}) of all the algorithms, `Macop` wants to allow users to quickly focus on one of the main parts of this loop.

![Macop common interation loop.\label{fig:macop-behaviour}](docs/source/_static/documentation/macop_behaviour_reduced.png)

# Statement of Need

Most of the operational research libraries developed in Python offer users either problems and algorithms where it is possible to choose parameters to obtain optimal (or near optimal) results such as proposed in [@10.1007/978-3-319-42432-3_37], or, libraries targeted to a specific problem or algorithm such as [@simanneal-solver]. Another package is proposed in [@solid-solver] which is a comprehensive gradient-free optimization framework written in Python. It seems very similar to `Macop`. However, hiearchic dependencies between algorithms, the notion of callbacks and adaptive operator selection are proposed within `Macop`.

On the other hand, available libraries [@pyopt-paper; @hart2017pyomo] in the literature did not allow to attach custom evaluation function to each algorithm used in this hierarchy of algorithms.
Indeed, it is sometimes possible that the main algorithm manages local searches. Each local search may evaluate the solution differently using a different evaluation function of the parent algorithm (the main algorithm). Such as example, using a surrogate mathematical model [@10.1145/3321707.3321800] with a quick-to-evaluate function if the real evaluation function is very expensive in time. This is why in `Macop`, each algorithm can have its own mechanism (or partially), i.e. its evaluation function, its operators for obtaining new solution, as well as its solution update policy. This is independent of the parent algorithm to which it is linked. This means that only the results (solutions found) are exchanged.

Hence, motivation behind **Macop** is a flexible discrete optimisation package allowing a quick implementation of problems. In particular it meets the following needs:

- **Common basis:** the interaction loop during the solution finding process proposed within the package is common to all heuristics. This allows the user to modify only a part of this interaction loop if necessary without rendering the process non-functional;
- **Hierarchy:** a hierarchical algorithm management system is available, especially when an algorithm needs to manage local searches. This hierarchy remains transparent to the user. The main algorithm will be able to manage and control the process of searching for solutions;
- **Flexibility:** although the algorithms are dependent on each other, it is possible that their internal management (search mechanism) is different. This means that the ways in which solutions are evaluated and updated, for example, may be different;
- **Abstraction:** thanks to the modular separability of the package, it is quickly possible to implement new problems, solutions representation, way to evaluate, update solutions within the package;
- **Extensible:** the package is open to extension, i.e. it does not partition the user in these developer choices. It can just as well implement continuous optimization problems if needed while making use of the main interaction loop proposed by the package;
- **Easy Setup:** as a pure Python package distributed is `pip` installable and easy to use.

# Target Audience 

This package would meet the expectations of people wishing to: 

- Solve a problem using an evolutionary algorithm but without developing their own frawmework. They can rely on what the package already proposes but also on its generic and flexible contribution in order to adapt their own content;
- Conduct research work leading to the rapid modification of meta-heuristics and the interaction of different algorithms. More precisely:
  
  - test new combinations of algorithms. Changing algorithms during evaluations, e.g. different local searches;
  - provide reinforcement learning during searches (e.g. adaptive operator choice strategy).
  - test new multi-objective methods quickly thanks to the proposed algorithmic hierarchy allowing to easily decompose the multi-objective problem into single-objective sub-problems.
  
- Take advantage of a system for launching calculations from a backup in order to avoid any loss in case of unwanted program interruption;
- Quickly model a problem that is still unknown, i.e. the type of solution and the evaluation function, while taking advantage of the interaction loop proposed by the package.

# Description

At the beginning of the development of this library, the idea of making it as modular as possible was topical. The library divide into sub-module forms considered to be the most important to build and solve an optimisation problem.

The package consists of main several modules:

- **solutions:** representation of the solution;
- **validator:** such as constraint programming, a `validator` is a function which is used to validate or not a solution data state;
- **evaluator:** stores problem instance data and implements a `compute` method in order to evaluate a solution;
- **operators:** mutators, crossovers operators to update and obtain new solution;
- **policies:** the way you choose the available operators (might be using reinforcement learning);
- **algorithms:** generic and implemented optimisation research algorithms;
- **callbacks:** callbacks to automatically keep track of the search space advancement and restart from previous state if nedded.


The primary advantage of using Python is that it allows you to dynamically add new members within the new implemented solution or algorithm classes. This of course does not close the possibilities of extension and storage of information within solutions and algorithms. It all depends on the current need.

## In `macop.algorithms` module:

Both single and multi-objective algorithms have been implemented for demonstration purposes. 

A hierarchy between dependent algorithms is also available, based on a parent/child link, allowing quick access to global information when looking for solutions, such as the best solution found, the number of global evaluations.

The mono-objective Iterated Local Search [@Lourenço2003] algorithm has been implemented. This algorithm aims to perform local searches (child algorithms linked to the main algorithm) and then to explore again (explorations vs. exploitation trade-off). On the multi-objective side, the MOEA/D algorithm [@DBLP:journals/tec/ZhangL07] has been implemented by using the weighted-sum of objectives to change multi-objectives problem into a set of mono-objective problems (Tchebycheff approach can also be used [@DBLP:journals/cor/AlvesA07]). Hence, this algorithm aims at decomposing the multi-objective problem into $\mu$ single-objective problems in order to obtain the Pareto front [@kim2005adaptive] where single-objective problems are so-called child algorithms linked to the multi-objective algorithm.

The main purpose of these developed algorithms is to show the possibilities of operational search algorithm implementations based on the minimalist structure of the library.

## In `macop.solutions` module:

Currently, only combinatorial solutions (discrete problem modelisation) are offered, with the well-known problem of the knapsack as an example. Of course, it's easy to add your own representations of solutions. Solutions modeling continuous problems can also be created by anyone who wants to model his own problem.

## In `macop.operators` and `macop.policies` modules:

A few mutation and crossover operators have been implemented. However, it remains quite simple. What is interesting here is that it is possible to develop one's own strategy for choosing operators for the next evaluation. The available UCBPolicy class proposes this functionality as an example, since it will seek to propose the best operator to apply based on a method known as the Adaptive Operator Selection (AOS) via the use of the Upper Confidence Bound (UCB) algorithm [@DBLP:journals/tec/LiFKZ14]. 

## In `macop.callbacks` module:

The use of callback instance allows both to do an action every $k$ evaluations of information, but also to reload them once the run of the algorithm is cut. Simply inherit the abstract Callback class and implement the `apply` method to backup and `load` to restore. It is possible to add as many callbacks as required. As an example, the implemented UCBPolicy has its own callback allowing the instance to reload previously collected statistics and restart using them.

# Conclusion

`Macop` aims to allow the modelling of discrete (usually combinatorial) optimisation problems. It is therefore open to expansion and not closed specifically to a kind of problem.

`Macop` proposes a simple structure of interaction of the main elements (algorithms, operators, solutions, policies, callbacks) for the resolution of operational research problems inside an interaction loop. From its generic structure, it is possible, thanks to the flexible programming paradigm of the Python language, to easily allow the extension and development of new algorithms and problems. Based on simple concepts, this package can therefore meet the needs of the rapid problem implementation.

# Acknowledgements

This work is supported by *Agence Nationale de la Recherche* : project ANR-17-CE38-0009

# References# Minimalist And Customisable Optimisation Package

[![status](https://joss.theoj.org/papers/9ea7d55c4fa83808f96929cb87adff3e/status.svg)](https://joss.theoj.org/papers/9ea7d55c4fa83808f96929cb87adff3e) 
[![](https://img.shields.io/github/workflow/status/jbuisine/macop/build)](https://github.com/jbuisine/macop/actions/workflows/python-app.yml)
[![PyPI](https://img.shields.io/pypi/v/macop)](https://pypi.org/project/macop/) 
[![](https://img.shields.io/pypi/dm/macop)](https://pypi.org/project/macop/)
[![GitHub](https://img.shields.io/github/license/jbuisine/macop?style=flat)](https://github.com/jbuisine/macop/blob/master/LICENCE)

<p align="center">
    <img src="https://github.com/jbuisine/macop/blob/master/docs/source/_static/logo_macop.png" alt="" width="40%">
</p>


## Description

`Macop` is a python package for solving discrete optimisation problems in nature. Continuous optimisation can also applicable if needed. The objective is to allow a user to exploit the basic structure proposed by this package to solve a problem specific to him. The interest is that he can quickly abstract himself from the complications related to the way of evaluating, comparing, saving the progress of the search for good solutions but rather concentrate if necessary on his own algorithm. Indeed, `Macop` offers the following main and basic features: 

- **solutions:** representation of the solution;
- **validator:** such as constraint programming, a `validator` is a function which is used to validate or not a solution data state;
- **evaluator:** stores problem instance data and implements a `compute` method in order to evaluate a solution;
- **operators:** mutators, crossovers operators to update and obtain new solution;
- **policies:** the way you choose the available operators (might be using reinforcement learning);
- **algorithms:** generic and implemented optimisation research algorithms;
- **callbacks:** callbacks to automatically keep track of the search space advancement and restart from previous state if nedded.

<p align="center">
    <img src="https://github.com/jbuisine/macop/blob/master/docs/source/_static/documentation/macop_behaviour.png" alt="" width="50%">
</p>


## Motivation

Flexible discrete optimisation package allowing a quick implementation of your problems. In particular it meets the following needs:

- **Common basis:** the interaction loop during the solution finding process proposed within the package is common to all heuristics. This allows the user to modify only a part of this interaction loop if necessary without rendering the process non-functional.
- **Hierarchy:** a hierarchical algorithm management system is available, especially when an algorithm needs to manage local searches. This hierarchy remains transparent to the user. The main algorithm will be able to manage and control the process of searching for solutions.
- **Flexibility:** although the algorithms are dependent on each other, it is possible that their internal management is different. This means that the ways in which solutions are evaluated and updated, for example, may be different.
- **Abstraction:** thanks to the modular separability of the package, it is quickly possible to implement new problems, solutions representation, way to evaluate, update solutions within the package.
- **Extensible:** the package is open to extension, i.e. it does not partition the user in these developer choices. It can just as well implement continuous optimization problems if needed while making use of the main interaction loop proposed by the package.
- **Easy Setup:** as a pure Python package distributed is ``pip`` installable and easy to use.

## Target Audience 

This package would meet the expectations of people wishing to: 
- Solve a problem using an evolutionary algorithm but without developing their own frawmework. They can rely on what the package already proposes but also on its generic and flexible contribution in order to adapt their own content;
- Conduct research work leading to the rapid modification of meta-heuristics and the interaction of different algorithms. More precisely:
  - test new combinations of algorithms. Changing algorithms during evaluations, e.g. different local searches;
  - provide reinforcement learning during searches (e.g. adaptive operator choice strategy).
  - test new multi-objective methods quickly thanks to the proposed algorithmic hierarchy allowing to easily decompose the multi-objective problem into single-objective sub-problems.
- Take advantage of a system for launching calculations from a backup in order to avoid any loss in case of unwanted program interruption;
- Quickly model a problem that is still unknown, i.e. the type of solution and the evaluation function, while taking advantage of the interaction loop proposed by the package.

## Content

The primary advantage of using Python is that it allows you to dynamically add new members within the new implemented solution or algorithm classes. This of course does not close the possibilities of extension and storage of information within solutions and algorithms. It all depends on the need in question.

### In `macop.algorithms` module:

Both single and multi-objective algorithms have been implemented for demonstration purposes. 

A hierarchy between dependent algorithms is also available, based on a parent/child link, allowing quick access to global information when looking for solutions, such as the best solution found, the number of global evaluations.

The mono-objective Iterated Local Search algorithm which aims to perform local searches (child algorithms linked to the main algorithm) and then to explore again (explorations vs. exploitation trade-off). On the multi-objective side, the MOEA/D algorithm has been implemented by using the weighted-sum of objectives to change multi-objectives problem into a set of mono-objective (Tchebycheff approach can also be used). Hence, this algorithm aims at decomposing the multi-objective problem into `mu` single-objective problems in order to obtain the Pareto front where single-objective problems are so-called child algorithms linked to the multi-objective algorithm.

The main purpose of these developed algorithms is to show the possibilities of operational search algorithm implementations based on the minimalist structure of the library.

### In `macop.solutions` module:

Currently, only combinatorial solutions (discrete problem modelisation) are offered, with the well-known problem of the knapsack as an example. Of course, it's easy to add your own representations of solutions. Solutions modeling continuous problems can also be created by the anyone who wants to model his own problem.

### In `macop.operators` and `macop.policies` modules:

A few mutation and crossover operators have been implemented, however, it remains quite simple. What is interesting here is that it is possible to develop one's own strategy for choosing operators for the next evaluation. The available UCBPolicy class proposes this functionality as an example, since it will seek to propose the best operator to apply based on a method known as the Adaptive Operator Selection (AOS) via the use of the Upper Confidence Bound (UCB) algorithm. 

### In `macop.callbacks` module:

The use of callback instance, allows both to do an action every $k$ evaluations of information, but also to reload them once the run of the algorithm is cut. Simply inherit the abstract Callback class and implement the `apply` method to backup and `load` to restore. It is possible to add as many callbacks as required. Such as an example, implemented UCBPolicy has its own callback allowing the instance to reload previously collected statistics and restart using them.


## Documentation

Based on all of these generic and/or implemented functionalities, the user will be able to quickly develop a solution to his problem while retaining the possibility of remaining in control of his development by overloading existing functionalities if necessary.

Main idea about this Python package is that it does not which doesn't implement every algorithm in the literature but let the possibility to the user to quickly develop and test its own algorithms and strategies. The main objective of this package is to provide maximum flexibility, which allows for easy experimentation in implementation..

Fully documentation of package with examples is [available](https://jbuisine.github.io/macop). 

You can also see examples of use:
-  in the [knapsackExample.py](https://github.com/jbuisine/macop/blob/master/examples/knapsackExample.py) for mono-objective knapsack instance;
-  in the [knapsackMultiExample.py](https://github.com/jbuisine/macop/blob/master/examples/knapsackMultiExample.py) for multi-objective knapsack instance;
-  in the [qapExample.py](https://github.com/jbuisine/macop/blob/master/examples/qapExample.py) for mono-objective QAP instance;
-  in the [ubqpExample.py](https://github.com/jbuisine/macop/blob/master/examples/ubqpExample.py) for mono-objective UBQP problem instance;
-  in the [ZdtExample.py](https://github.com/jbuisine/macop/blob/master/examples/ZdtExample.py) for continuous optimisation problem over Zdt function.

## Add as dependency

```bash
git submodule add https://github.com/jbuisine/macop.git
```

## Current projects which use `Macop`:

- [@prise3d/noise-detection-attributes-optimization](https://github.com/prise-3d/noise-detection-attributes-optimization): use of a parent algorithm with the real (but very expensive) evaluation function, and then inner local searches which use a substitution model (a model that has learned to approximate the real cost function with a quick-to-evaluate function). Hence, two evaluation functions have been used in order to accelerate the search in the set of solutions.

## License

[The MIT License](LICENSE)How to contribute
=====================================

<p align="center">
    <img src="https://github.com/jbuisine/macop/blob/master/docs/source/_static/logo_macop.png" alt="" width="40%">
</p>



# Welcome !

Thank you for taking the time to read this guide for the package's contribution. I'm glad to know that you may bring a lot to the **Macop** package. This document will show you the good development practices used in the project and how you can easily participate in its evolution!

# Table of contents

1. [Submission processes](#submission-process)

    1.1. [Submit an issue](#submit-an-issue)

    1.2. [Pull request](#pull-request)

    1.3. [Seek support](#seek-support)

2. [Coding conventions](#coding-conventions)

    2.1. [Python conventions](#python-conventions)

    2.2. [Code documentation](#code-documentation)

    2.3. [Testing](#test-implementation)


# Submission process

## Submit an issue

Do not hesitate to report bug or issue in [https://github.com/jbuisine/macop/issues](https://github.com/jbuisine/macop/issues) with the common template header:

```
**Package version:** X.X.X
**Issue label:** XXXXX
**Targeted modules:** `macop.algorithms`, `macop.policies`
**Operating System:** Manjaro Linux

**Description:** XXXXX
```

## Pull request

If you have made changes to the project you have forked, you can submit a pull request in [https://github.com/jbuisine/macop/pulls](https://github.com/jbuisine/macop/pulls) in order to have your changes added inside new version of the `Macop` package. A [GitHub documentation](https://help.github.com/articles/about-pull-requests/) about pull requests is available if necessary.

To enhance the package, do not hesitate to fix bug or missing feature. To do that, just submit your pull request with this common template header:

```
**Package version:** X.X.X
**Enhancements label:** XXXXX
**Targeted modules:** `macop.algorithms`, `macop.policies`
**New modules:** `macop.XXXX`, `macop.algorithms.XXXX`

**Description:** XXXXX
```

**Note:** the code conventions required for the approval of your changes are described below.

Whatever the problem reported, I will thank you for your contribution to this project. So do not hesitate.

## Seek support

If you have any problem with the use of the package, issue or pull request submission, do not hesitate to let a message to [https://github.com/jbuisine/macop/discussions](https://github.com/jbuisine/macop/discussions). Especially in the question and answer section. 

You can also contact me at the following email address: `contact@jeromebuisine.fr`.

# Coding conventions

## Python conventions

This project follows the [coding conventions](http://google.github.io/styleguide/pyguide.html) implemented by Google. To help you to format **\*.py** files, it is possible to use the [yapf](https://github.com/google/yapf/) Python package developed by Google.

```
yapf -ir -vv macop
```

**Note:** you need at least Python version >=3.7.0.

## Package modules conventions

As you perhaps already saw, package contains multiples modules and submodules. It's really import to be well organized package and let it intuitive to access as possible to features.

`Macop` is mainly decompose into discrete and continuous optimisation. Especially if you want to offer continuous optimisation problems, modules are already available for this purpose. You can refer to the [documentation](https://jbuisine.github.io/macop) if necessary.

In order to facilitate the integration of new modules, do not hesitate to let me know the name it could have beforehand in your pull request.

## Code documentation

In order to allow quick access to the code, the project follows the documentation conventions (docstring) proposed by Google. Here an example:

```python
class BinarySolution():
"""Binary integer solution class

    - store solution as a binary array (example: [0, 1, 0, 1, 1])
    - associated size is the size of the array
    - mainly use for selecting or not an element in a list of valuable objects

    Attributes:
       data: {ndarray} --  array of binary values
       size: {int} -- size of binary array values
       score: {float} -- fitness score value
"""
```

For method:
```python
class BinarySolution():

...

def random(self, validator):
    """
    Intialize binary array with use of validator to generate valid random solution

    Args:
        size: {int} -- expected solution size to generate
        validator: {function} -- specific function which validates or not a solution (if None, not validation is applied)

    Returns:
        {:class:`~macop.solutions.discrete.BinarySolution`}: new generated binary solution
    """
    ...
```

You can generate documentation and display updates using these following commands:

```
cd docs
make clean && make html
firefox _build/html/index.html
```

## Test implementation

This project uses the [doctest](https://docs.python.org/3/library/doctest.html) package which enables to write tests into documentation as shown in example below:

```python
""" Initialise binary solution using specific data

    Args:
        data: {ndarray} --  array of binary values
        size: {int} -- size of binary array values

    Example:

    >>> from macop.solutions.discrete import BinarySolution
    >>> # build of a solution using specific data and size
    >>> data = [0, 1, 0, 1, 1]
    >>> solution = BinarySolution(data, len(data))
    >>> # check data content
    >>> sum(solution.data) == 3
    True
    >>> # clone solution
    >>> solution_copy = solution.clone()
    >>> all(solution_copy._data == solution.data)
"""
```

Moreover, tests written are displayed into generated documentation and show examples of how to use the developed features.
Contributing
=====================================

.. image:: _static/logo_macop.png
   :width: 350 px
   :align: center

Using GitHub
------------

Please refer to the guidelines_ file if you want more information about process!

.. _guidelines: https://github.com/prise-3d/macop/blob/master/CONTRIBUTING.mdExample uses
=====================================

You will find here some examples of using Macop with implementations of problems well known from the literature.

Discrete problem
----------------------------

.. toctree::
   :maxdepth: 1

   qap_example
   ubqp_example

Continuous problem
----------------------------

.. toctree::
   :maxdepth: 1

   zdt_example


Available code examples
-----------------------

- mono-objective knapsack problem: knapsackExample.py_
- multi-objective knapsack problem: knapsackMultiExample.py_
- QAP problem: qapExample.py_
- UBQP problem: ubqpExample.py_
- Continuous Zdt optimisation problem: ZdtExample.py_

.. _knapsackExample.py: https://github.com/jbuisine/macop/blob/master/examples/knapsackExample.py
.. _knapsackMultiExample.py: https://github.com/jbuisine/macop/blob/master/examples/knapsackMultiExample.py
.. _qapExample.py: https://github.com/jbuisine/macop/blob/master/examples/qapExample.py
.. _ubqpExample.py: https://github.com/jbuisine/macop/blob/master/examples/ubqpExample.py
.. _ZdtExample.py: https://github.com/jbuisine/macop/blob/master/examples/ZdtExample.pyDescription
=====================================

.. image:: _static/logo_macop.png
   :width: 350 px
   :align: center


Context
------------

Based on its generic behaviour, each **Macop** algorithm runs can be represented as an interactive loop where you can interact with and specify your needs at each step:

.. image:: _static/documentation/macop_behaviour.png
   :width: 450 px
   :align: center

The package is strongly oriented on combinatorial optimisation (hence discrete optimisation) but it remains possible to extend for continuous optimisation.

Motivation
~~~~~~~~~~

Flexible discrete optimisation package allowing a quick implementation of your problems. In particular it meets the following needs:

- **Common basis:** the interaction loop during the solution finding process proposed within the package is common to all heuristics. This allows the user to modify only a part of this interaction loop if necessary without rendering the process non-functional.
- **Hierarchy:** a hierarchical algorithm management system is available, especially when an algorithm needs to manage local searches. This hierarchy remains transparent to the user. The main algorithm will be able to manage and control the process of searching for solutions.
- **Flexibility:** although the algorithms are dependent on each other, it is possible that their internal management is different. This means that the ways in which solutions are evaluated and updated, for example, may be different.
- **Abstraction:** thanks to the modular separability of the package, it is quickly possible to implement new problems, solutions representation, way to evaluate, update solutions within the package.
- **Extensible:** the package is open to extension, i.e. it does not partition the user in these developer choices. It can just as well implement continuous optimization problems if needed while making use of the main interaction loop proposed by the package.
- **Easy Setup:** as a pure Python package distributed is ``pip`` installable and easy to use.


Target Audience
~~~~~~~~~~~~~~~

This package would meet the expectations of people wishing to: 

- Solve a problem using an evolutionary algorithm but without developing their own frawmework. They can rely on what the package already proposes but also on its generic and flexible contribution in order to adapt their own content;
- Conduct research work leading to the rapid modification of meta-heuristics and the interaction of different algorithms. More precisely:

  - test new combinations of algorithms. Changing algorithms during evaluations, e.g. different local searches;
  - provide reinforcement learning during searches (e.g. adaptive operator choice strategy).
  - test new multi-objective methods quickly thanks to the proposed algorithmic hierarchy allowing to easily decompose the multi-objective problem into single-objective sub-problems.

- Take advantage of a system for launching calculations from a backup in order to avoid any loss in case of unwanted program interruption;
- Quickly model a problem that is still unknown, i.e. the type of solution and the evaluation function, while taking advantage of the interaction loop proposed by the package.

Installation
------------

Just install package using `pip` Python package manager: 

.. code:: bash
   
   pip install macop
===============================
Quadratric Assignment Problem
===============================

This example will deal with the use of the **Macop** package in relation to a quadratic assignment problem (QAP). We will use a known example of this problem to associate a set of facilities (:math:`F`) to a set of locations (:math:`L`).

.. image:: _static//examples/qap/factories_qap.png
   :width: 50 %
   :align: center
   :alt: Example of QAP facilities to locations problem


.. note:: 
   The full code for what will be proposed in this example is available: qapExample.py_.

.. _qapExample.py: https://github.com/jbuisine/macop/blob/master/examples/qapExample.py


QAP problem definition
======================

The quadratic assignment problem (QAP) was introduced by Koopmans and Beckman in 1957 in the context of locating "indivisible economic activities". The objective of the problem is to assign a set of facilities to a set of locations in such a way as to minimize the total assignment cost. The assignment cost for a pair of facilities is a function of the flow between the facilities and the distance between the locations of the facilities.

Location assignment example
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider a **facility location problem** with **four** facilities (and **four** locations). One possible assignment is shown in the figure below: facility 4 is assigned to location 1, facility 1 
is assigned to location 2, facility 3 is assigned to location 3, and facility 2 is assigned to location 3. This assignment can be written as the permutation :math:`p=\{4,1,3,2\}`, 
which means that facility 4 is assigned to location 1, facility 1 is assigned to location 2, facility 3 is assigned to location 3, and facility 2 is assigned to location 3. 
In the figure, the line between a pair of facilities indicates that there is required flow between the facilities, and the thickness of the line increases with the value of the flow. 

.. image:: _static/examples/qap/factories_qap.png
   :width: 50 %
   :align: center
   :alt: Example of QAP facilities to locations problem


To calculate the assignment cost of the permutation, the required flows between facilities and the distances between locations are needed.


.. tabularcolumns:: |p{1cm}|p{1cm}|p{1cm}|p{1cm}|

.. csv-table:: flow of the current facilities
   :header: facility `i`, facility `j`, flow( `i`\, `j` )
   :widths: 2, 2, 3

   1, 4, 4
   3, 4, 10  
   3, 1, 8
   2, 1, 6  


.. csv-table:: distances of between locations
   :header: location `i`, location `j`, distances( `i`\, `j` )
   :widths: 2, 2, 3

   1, 2, 42
   1, 3, 30  
   2, 3, 41
   3, 4, 23  


Then, the assignment cost of the permutation can be computed as:

:math:`f(1,4)⋅d(1,2)+f(3,4)⋅d(1,3)+f(1,3)⋅d(2,3)+f(3,2)⋅d(3,4)` 
with result :math:`4⋅42+10⋅30+8⋅41+6⋅23=934`.

Note that this permutation is not the optimal solution.

Mathematical definition
~~~~~~~~~~~~~~~~~~~~~~~

**Sets**

- :math:`N=\{1,2,⋯,n\}`
- :math:`S_n=\phi:N→N` is the set of all permutations

**Parameters**

- :math:`F=(f_{ij})` is an :math:`n×n` matrix where :math:`f_{ij}` is the required flow between facilities :math:`i` and :math:`j`
- :math:`D=(d_{ij})` is an :math:`n×n` matrix where :math:`d_{ij}` is the distance between locations :math:`i` and :math:`j`.

**Optimization Problem**

- :math:`min_{ϕ∈S_n}\sum_{i=1}^{n}{\sum_{j=1}^{n}{f_{ij}⋅d_{\phi(i)\phi(j)}}}`

The assignment of facilities to locations is represented by a permutation :math:`\phi`, where :math:`\phi(i)` is the location to which facility :math:`i` is assigned. Each individual product :math:`f_{ij}⋅d_{\phi(i)\phi(j)}` is the cost of assigning facility :math:`i` to location :math:`\phi(i)` and facility :math:`j` to location :math:`\phi(j)`.

QAP Problem instance generation
===============================

To define our quadratic assignment problem instance, we will use the available mQAP_ multi-objective quadratic problem generator. 

Genration of the instance
~~~~~~~~~~~~~~~~~~~~~~~~~

We will limit ourselves here to a single objective for the purposes of this example. The file **makeQAPuni.cc**, will be used to generate the instance.

.. code:: bash

    g++ makeQAPuni.cc -o mQAPGenerator
    ./mQAPGenerator -n 100 -k 1 -f 30 -d 80 -s 42 > qap_instance.txt

with the following parameters:

- **-n** positive integer: number of facilities/locations;
- **-k** positive integer: number of objectives;
- **-f** positive integer: maximum flow between facilities;
- **-d** positive integer: maximum distance between locations;
- **-s** positive long: random seed.

The generated qap_instance.txt_ file contains the two matrices :math:`F` and :math:`D` and define our instance problem.

.. _mQAP: https://www.cs.bham.ac.uk/~jdk/mQAP/

.. _qap_instance.txt: https://github.com/jbuisine/macop/blob/master/examples/instances/qap/qap_instance.txt


Load data instance
~~~~~~~~~~~~~~~~~~


We are now going to load this instance via a Python code which will be useful to us later on:

.. code:: Python

    qap_instance_file = 'qap_instance.txt'

    n = 100 # the instance size

    with open(qap_instance_file, 'r') as f:
        file_data = f.readlines()
        print(f'Instance information {file_data[0]}')

        D_lines = file_data[1:n + 1]
        D_data = ''.join(D_lines).replace('\n', '')

        F_lines = file_data[n:2 * n + 1]
        F_data = ''.join(F_lines).replace('\n', '')

    D_matrix = np.fromstring(D_data, dtype=float, sep=' ').reshape(n, n)
    print(f'D matrix shape: {D_matrix.shape}')
    F_matrix = np.fromstring(F_data, dtype=float, sep=' ').reshape(n, n)
    print(f'F matrix shape: {F_matrix.shape}')

.. note::
    As we know the size of our instance and the structure of the document, it is quite quick to look for the lines related to the :math:`F` and :math:`D` matrices.

Macop QAP implementation
========================

Let's see how it is possible with the use of the **Macop** package to implement and deal with this QAP instance problem.

Solution structure definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Firstly, we are going to use a type of solution that will allow us to define the structure of our solutions.

The available macop.solutions.discrete.CombinatoryIntegerSolution_ type of solution within the Macop package represents exactly what one would wish for. 
I.e. a solution that stores a sequence of integers relative to the size of the problem, the order of which is not sorted.

Let's see an example of its use:

.. code:: python

    from macop.solutions.discrete import CombinatoryIntegerSolution
    
    solution = CombinatoryIntegerSolution.random(10)
    print(solution)


The resulting solution obtained:

.. code:: bash

    Combinatory integer solution [2 9 8 1 7 6 0 4 3 5]


QAP Evaluator
~~~~~~~~~~~~~

Now that we have the structure of our solutions, and the means to generate them, we will seek to evaluate them.

To do this, we need to create a new evaluator specific to our problem and the relative evaluation function:

- :math:`min_{ϕ∈S_n}\sum_{i=1}^{n}{\sum_{j=1}^{n}{f_{ij}⋅d_{\phi(i)\phi(j)}}}`

So we are going to create a class that will inherit from the abstract class macop.evaluators.base.Evaluator_:


.. code:: python

    from macop.evaluators.base import Evaluator

    class QAPEvaluator(Evaluator):
    """QAP evaluator class which enables to compute QAP solution using specific `_data`

    - stores into its `_data` dictionary attritute required measures when computing a QAP solution
    - `_data['F']` matrix of size n x n with flows data between facilities (stored as numpy array)
    - `_data['D']` matrix of size n x n with distances data between locations (stored as numpy array)
    - `compute` method enables to compute and associate a score to a given QAP solution
    """

    def compute(self, solution):
        """Apply the computation of fitness from solution

        Args:
            solution: {Solution} -- QAP solution instance
    
        Returns:
            {float} -- fitness score of solution
        """
        fitness = 0
        for index_i, val_i in enumerate(solution.getdata = )):
            for index_j, val_j in enumerate(solution.getdata = )):
                fitness += self._data['F'][index_i, index_j] * self._data['D'][val_i, val_j]

        return fitness

The cost function for the quadratic problem is now well defined.

.. tip::
    The class proposed here, is available in the Macop package macop.evaluators.discrete.mono.QAPEvaluator_.

Running algorithm
~~~~~~~~~~~~~~~~~

Now that the necessary tools are available, we will be able to deal with our problem and look for solutions in the search space of our QAP instance.

Here we will use local search algorithms already implemented in **Macop**.

If you are uncomfortable with some of the elements in the code that will follow, you can refer to the more complete **Macop** documentation_ that focuses more on the concepts and tools of the package.

.. code:: python

    # main imports
    import numpy as np

    # module imports
    from macop.solutions.discrete import CombinatoryIntegerSolution
    from macop.evaluators.discrete.mono import QAPEvaluator

    from macop.operators.discrete.mutators import SimpleMutation

    from macop.policies.classicals import RandomPolicy

    from macop.algorithms.mono import IteratedLocalSearch as ILS
    from macop.algorithms.mono import HillClimberFirstImprovment

    # usefull instance data
    n = 100
    qap_instance_file = 'qap_instance.txt'

    # default validator (check the consistency of our data, i.e. only unique element)
    def validator(solution):
        if len(list(solution.getdata = ))) > len(set(list(solution.getdata = )))):
            print("not valid")
            return False
        return True

    # define init random solution
    def init():
        return CombinatoryIntegerSolution.random(n, validator)

    # load qap instance
    with open(qap_instance_file, 'r') as f:
        file_data = f.readlines()
        print(f'Instance information {file_data[0]}')

        D_lines = file_data[1:n + 1]
        D_data = ''.join(D_lines).replace('\n', '')

        F_lines = file_data[n:2 * n + 1]
        F_data = ''.join(F_lines).replace('\n', '')

    D_matrix = np.fromstring(D_data, dtype=float, sep=' ').reshape(n, n)
    print(f'D matrix shape: {D_matrix.shape}')
    F_matrix = np.fromstring(F_data, dtype=float, sep=' ').reshape(n, n)
    print(f'F matrix shape: {F_matrix.shape}')

    # only one operator here
    operators = [SimpleMutation()]

    # random policy even if list of solution has only one element
    policy = RandomPolicy(operators)

    # use of loaded data from QAP instance
    evaluator = QAPEvaluator(data={'F': F_matrix, 'D': D_matrix})

    # passing global evaluation param from ILS
    hcfi = HillClimberFirstImprovment(init, evaluator, operators, policy, validator, maximise=False, verbose=True)
    algo = ILS(init, evaluator, operators, policy, validator, localSearch=hcfi, maximise=False, verbose=True)

    # run the algorithm
    bestSol = algo.run(10000, ls_evaluations=100)

    print('Solution for QAP instance score is {}'.format(evaluator.compute(bestSol)))


QAP problem solving is now possible with **Macop**. As a reminder, the complete code is available in the qapExample.py_ file.

.. _qapExample.py: https://github.com/jbuisine/macop/blob/master/examples/qapExample.py
.. _documentation: https://jbuisine.github.io/macop/_build/html/documentations


.. _macop.solutions.discrete.CombinatoryIntegerSolution: macop/macop.solutions.discrete.html#macop.solutions.discrete.CombinatoryIntegerSolution
.. _macop.evaluators.base.Evaluator: macop/macop.evaluators.base.html#macop.evaluators.base.Evaluator
.. _macop.evaluators.discrete.mono.QAPEvaluator: macop/macop.evaluators.discrete.mono.html#macop.evaluators.discrete.mono.QAPEvaluator===============================
Zdt optimisation problem
===============================

In applied mathematics, test functions, known as artificial landscapes, are useful to evaluate characteristics of continuous optimization algorithms, such as:

- Convergence rate.
- Precision.
- Robustness.
- General performance.

.. note:: 
   The full code for what will be proposed in this example is available: ZdtExample.py_.


Rosenbrock's function
======================

In mathematical optimization, the Rosenbrock function is a non-convex function, introduced by Howard H. Rosenbrock in 1960, which is used as a performance test problem for optimization algorithms.

Mathematical definition
~~~~~~~~~~~~~~~~~~~~~~~

The function is defined by: :math:`f(x, y) = (a − x)^2 + b(y − x^2)^2`

It has a global minimum at :math:`(x, y) = (a, a^2)`, where :math:`f(x, y) = 0`. Usually these parameters are set such that :math:`a = 1` and :math:`b = 100`. Only in the trivial case where :math:`a = 0` the function is symmetric and the minimum is at the origin. 

Below is a 3D representation of the function with the same parameters :math:`a = 1` and :math:`b = 100`.

.. image:: _static/examples/zdt/rosenbrock_function.jpg
   :width: 50 %
   :align: center
   :alt: 3D representation of Rosenbrock's function

The search space is defined by: :math:`-\infty \leq x_i \leq \infty, 1 \leq i \leq n`

Optimal solution is defined by: :math:`f(1, ..., 1)=0` when :math:`n > 3`

Specific instance used
~~~~~~~~~~~~~~~~~~~~~~

Using :math:`a = 1` and :math:`b = 100`, the function can be re-written:

- :math:`f(x)=\sum_{i=1}^{n-1}{[(x_{i + 1} − x_i^2)^2 + (1 − x_i)^2]}`


For the current implementation example, the search space will be reduced to :math:`-10 \leq x_i \leq 10` and the instance size will be set to :math:`n = 10`.

Macop implementation
========================

Let's see how it is possible with the use of the **Macop** package to implement and deal with this Rosenbrock's function instance problem.

Solution structure definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Firstly, we are going to use a type of solution that will allow us to define the structure of our solutions.

The available macop.solutions.continuous.ContinuousSolution_ type of solution within the Macop package represents exactly what one would wish for. 
I.e. a solution that stores a float array with respect to the size of the problem.

Let's see an example of its use:

.. code:: python

    from macop.solutions.continuous import ContinuousSolution
    
    problem_interval = -10, 10
    solution = ContinuousSolution.random(10, interval=problem_interval)
    print(solution)

The ``problem_interval`` variable is required in order to generate our continuous solution with respect to the search space.
The resulting solution obtained should be something like:

.. code:: bash

    Continuous solution [-3.31048093 -8.69195762 ... -2.84790964 -1.08397853]


Zdt Evaluator
~~~~~~~~~~~~~

Now that we have the structure of our solutions, and the means to generate them, we will seek to evaluate them.

To do this, we need to create a new evaluator specific to our problem and the relative evaluation function:

- :math:`f(x)=\sum_{i=1}^{n-1}{[(x_{i + 1} − x_i^2)^2 + (1 − x_i)^2]}`

So we are going to create a class that will inherit from the abstract class macop.evaluators.base.Evaluator_:


.. code:: python

    from macop.evaluators.base import Evaluator

    class ZdtEvaluator(Evaluator):
    """Generic Zdt evaluator class which enables to compute custom Zdt function for continuous problem

    - stores into its `_data` dictionary attritute required measures when computing a continuous solution
    - `_data['f']` stores lambda Zdt function 
    - `compute` method enables to compute and associate a score to a given continuous solution
    """

    def compute(self, solution):
        """Apply the computation of fitness from solution
        Args:
            solution: {:class:`~macop.solutions.base.Solution`} -- Solution instance
    
        Returns:
            {float}: fitness score of solution
        """
        return self._data['f'](solution)

The cost function for the zdt continuous problem is now well defined but we still need to define the lambda function.

.. code:: python

    from macop.evaluators.continuous.mono import ZdtEvaluator

    # Rosenbrock function definition
    Rosenbrock_function = lambda s: sum([ 100 * math.pow(s.data[i + 1] - (math.pow(s.data[i], 2)), 2) + math.pow((1 - s.data[i]), 2) for i in range(len(s.data) - 1) ])

    evaluator = ZdtEvaluator(data={'f': Rosenbrock_function})

.. tip::
    The class proposed here, is available in the Macop package macop.evaluators.continuous.mono.ZdtEvaluator_.

Running algorithm
~~~~~~~~~~~~~~~~~

Now that the necessary tools are available, we will be able to deal with our problem and look for solutions in the search space of our Zdt Rosenbrock instance.

Here we will use local search algorithms already implemented in **Macop**.

If you are uncomfortable with some of the elements in the code that will follow, you can refer to the more complete **Macop** documentation_ that focuses more on the concepts and tools of the package.

.. code:: python

    # main imports
    import numpy as np

    # module imports
    from macop.solutions.continuous import ContinuousSolution
    from macop.evaluators.continuous.mono import ZdtEvaluator

    from macop.operators.continuous.mutators import PolynomialMutation

    from macop.policies.classicals import RandomPolicy

    from macop.algorithms.mono import IteratedLocalSearch as ILS
    from macop.algorithms.mono import HillClimberFirstImprovment

    # usefull instance data
    n = 10
    problem_interval = -10, 10
    qap_instance_file = 'zdt_instance.txt'

    # default validator (check the consistency of our data, i.e. x_i element in search space)
    def validator(solution):
        mini, maxi = problem_interval

        for x in solution.data:
            if x < mini or x > maxi:
                return False

        return True

    # define init random solution with search space bounds
    def init():
        return ContinuousSolution.random(n, interval=problem_interval, validator)

    # only one operator here
    operators = [PolynomialMutation()]

    # random policy even if list of solution has only one element
    policy = RandomPolicy(operators)

    # Rosenbrock function definition
    Rosenbrock_function = lambda s: sum([ 100 * math.pow(s.data[i + 1] - (math.pow(s.data[i], 2)), 2) + math.pow((1 - s.data[i]), 2) for i in range(len(s.data) - 1) ])

    evaluator = ZdtEvaluator(data={'f': Rosenbrock_function})

    # passing global evaluation param from ILS
    hcfi = HillClimberFirstImprovment(init, evaluator, operators, policy, validator, maximise=False, verbose=True)
    algo = ILS(init, evaluator, operators, policy, validator, localSearch=hcfi, maximise=False, verbose=True)

    # run the algorithm
    bestSol = algo.run(10000, ls_evaluations=100)

    print('Solution for zdt Rosenbrock instance score is {}'.format(evaluator.compute(bestSol)))


Continuous Rosenbrock's function problem is now possible with **Macop**. As a reminder, the complete code is available in the ZdtExample.py_ file.

.. _ZdtExample.py: https://github.com/jbuisine/macop/blob/master/examples/ZdtExample.py
.. _documentation: https://jbuisine.github.io/macop/_build/html/documentations


.. _macop.solutions.continuous.ContinuousSolution: macop/macop.solutions.continuous.html#macop.solutions.continuous.ContinuousSolution
.. _macop.evaluators.base.Evaluator: macop/macop.evaluators.base.html#macop.evaluators.base.Evaluator
.. _macop.evaluators.continuous.mono.ZdtEvaluator: macop/macop.evaluators.continuous.mono.html#macop.evaluators.continuous.mono.ZdtEvaluatorMinimalist And Customisable Optimisation Package
================================================

.. image:: _static/logo_macop.png
   :width: 300 px
   :align: center

What's **Macop** ?
~~~~~~~~~~~~~~~~~~

**Macop** is a discrete optimisation Python package which not which doesn't implement every algorithm in the literature but provides the ability to quickly develop and test your own algorithm and strategies. The main objective of this package is to provide maximum flexibility, which allows for easy experimentation in implementation.


Contents
~~~~~~~~

.. toctree::
   :maxdepth: 1

   description

   documentations

   api

   examples

   contributing


Motivation
~~~~~~~~~~

Flexible discrete optimisation package allowing a quick implementation of your problems. In particular it meets the following needs:

- **Common basis:** the interaction loop during the solution finding process proposed within the package is common to all heuristics. This allows the user to modify only a part of this interaction loop if necessary without rendering the process non-functional.
- **Hierarchy:** a hierarchical algorithm management system is available, especially when an algorithm needs to manage local searches. This hierarchy remains transparent to the user. The main algorithm will be able to manage and control the process of searching for solutions.
- **Flexibility:** although the algorithms are dependent on each other, it is possible that their internal management is different. This means that the ways in which solutions are evaluated and updated, for example, may be different.
- **Abstraction:** thanks to the modular separability of the package, it is quickly possible to implement new problems, solutions representation, way to evaluate, update solutions within the package.
- **Extensible:** the package is open to extension, i.e. it does not partition the user in these developer choices. It can just as well implement continuous optimization problems if needed while making use of the main interaction loop proposed by the package.
- **Easy Setup:** as a pure Python package distributed is `pip` installable and easy to use.


Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
==========================================
Unconstrained Binary Quadratic Programming
==========================================

Given a collection of :math:`n` items such that each pair of items is associated with a profit value that can be positive, negative or zero, unconstrained binary quadratic programming (UBQP) seeks a subset of items that maximizes the sum of their paired values. The value of a pair is accumulated in the sum only if the two corresponding items are selected.

The UBQP problem will be tackle in this example.

UBQP problem definition
=======================

Understand the UBQP Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a collection of :math:`n` items such that each pair of items is associated with a profit value that can be positive, negative or zero, unconstrained binary quadratic programming (UBQP) seeks a subset of items that maximizes the sum of their paired values. The value of a pair is accumulated in the sum only if the two corresponding items are selected. A feasible solution to a UBQP instance can be specified by a binary string of size :math:`n`, such that each variable indicates whether the corresponding item is included in the selection or not.


Mathematical definition
~~~~~~~~~~~~~~~~~~~~~~~

More formally, the conventional and single-objective UBQP problem is to maximize the following objective function:

:math:`f(x)=x′Qx=\sum_{i=1}^{n}{\sum_{j=1}^{n}{q_{ij}⋅x_i⋅x_j}}`

where :math:`Q=(q_{ij})` is an :math:`n` by :math:`n` matrix of constant values, :math:`x` is a vector of :math:`n` binary (zero-one) variables, i.e., :math:`x \in \{0, 1\}`, :math:`i \in \{1,...,n\}`, and :math:`x'` is the transpose of :math:`x`.

UBQP Problem instance generation
================================

To define our quadratic assignment problem instance, we will use the available mUBQP_ multi-objective quadratic problem generator. 

Genration of the instance
~~~~~~~~~~~~~~~~~~~~~~~~~

We will limit ourselves here to a single objective for the purposes of this example. The available file **mubqpGenerator.R**, will be used to generate the instance (using R language).

.. code:: bash

    Rscript mubqpGenerator.R 0.8 1 100 5 42 ubqp_instance.txt

The main parameters used for generating our UBQP instance are:

- **ρ:** the objective correlation coefficient
- **M:** the number of objective functions
- **N:** the length of bit strings
- **d:** the matrix density (frequency of non-zero numbers)
- **s:** seed to use

.. _mUBQP: http://mocobench.sourceforge.net/index.php?n=Problem.MUBQP

.. _ubqp_instance.txt: https://github.com/jbuisine/macop/blob/master/examples/instances/ubqp/ubqp_instance.txt

Load data instance
~~~~~~~~~~~~~~~~~~

We are now going to load this instance via a Python code which will be useful to us later on:

.. code:: Python

    qap_instance_file = 'ubqp_instance.txt'

    n = 100 # the instance size

    # load UBQP instance
    with open(ubqp_instance_file, 'r') as f:

        lines = f.readlines()

        # get all string floating point values of matrix
        Q_data = ''.join([ line.replace('\n', '') for line in lines[8:] ])

        # load the concatenate obtained string
        Q_matrix = np.fromstring(Q_data, dtype=float, sep=' ').reshape(n, n)

    print(f'Q_matrix {Q_matrix.shape}')

.. note::
    As we know the size of our instance and the structure of the document (header size), it is quite quick to look for the lines related to the :math:`Q` matrix.

Macop UBQP implementation
=========================

Let's see how it is possible with the use of the **Macop** package to implement and deal with this UBQP instance problem.

Solution structure definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Firstly, we are going to use a type of solution that will allow us to define the structure of our solutions.

The available macop.solutions.discrete.BinarySolution_ type of solution within the Macop package represents exactly what one would wish for. 

Let's see an example of its use:

.. code:: python

    from macop.solutions.discrete import BinarySolution
    
    solution = BinarySolution.random(10)
    print(solution)


The resulting solution obtained:

.. code:: bash

    Binary solution [1 0 1 1 1 0 0 1 1 0]


UBQP Evaluator
~~~~~~~~~~~~~~

Now that we have the structure of our solutions, and the means to generate them, we will seek to evaluate them.

To do this, we need to create a new evaluator specific to our problem and the relative evaluation function we need to maximise:

- :math:`f(x)=x′Qx=\sum_{i=1}^{n}{\sum_{j=1}^{n}{q_{ij}⋅x_i⋅x_j}}`

So we are going to create a class that will inherit from the abstract class macop.evaluators.base.Evaluator_:

.. code:: python

    from macop.evaluators.base import Evaluator

    class UBQPEvaluator(Evaluator):
    """UBQP evaluator class which enables to compute UBQP solution using specific `_data`

    - stores into its `_data` dictionary attritute required measures when computing a UBQP solution
    - `_data['Q']` matrix of size n x n with real values data (stored as numpy array)
    - `compute` method enables to compute and associate a score to a given UBQP solution
    """

    def compute(self, solution):
        """Apply the computation of fitness from solution

        Args:
            solution: {Solution} -- UBQP solution instance
    
        Returns:
            {float} -- fitness score of solution
        """
        fitness = 0
        for index_i, val_i in enumerate(solution.getdata = )):
            for index_j, val_j in enumerate(solution.getdata = )):
                fitness += self._data['Q'][index_i, index_j] * val_i * val_j

        return fitness

The cost function for the Unconstrained binary quadratic problem is now well defined.

.. tip::
    The class proposed here, is available in the Macop package macop.evaluators.discrete.mono.UBQPEvaluator_.

Running algorithm
~~~~~~~~~~~~~~~~~

Now that the necessary tools are available, we will be able to deal with our problem and look for solutions in the search space of our UBQP instance.

Here we will use local search algorithms already implemented in **Macop**.

If you are uncomfortable with some of the elements in the code that will follow, you can refer to the more complete **Macop** documentation_ that focuses more on the concepts and tools of the package.

.. code:: python

    # main imports
    import numpy as np

    # module imports
    from macop.solutions.discrete import BinarySolution
    from macop.evaluators.discrete.mono import UBQPEvaluator

    from macop.operators.discrete.mutators import SimpleMutation, SimpleBinaryMutation

    from macop.policies.classicals import RandomPolicy

    from macop.algorithms.mono import IteratedLocalSearch as ILS
    from macop.algorithms.mono import HillClimberFirstImprovment

    # usefull instance data
    n = 100
    ubqp_instance_file = 'ubqp_instance.txt'

    # default validator
    def validator(solution):
        return True

    # define init random solution
    def init():
        return BinarySolution.random(n, validator)

    # load UBQP instance
    with open(ubqp_instance_file, 'r') as f:

        lines = f.readlines()

        # get all string floating point values of matrix
        Q_data = ''.join([ line.replace('\n', '') for line in lines[8:] ])

        # load the concatenate obtained string
        Q_matrix = np.fromstring(Q_data, dtype=float, sep=' ').reshape(n, n)

    print(f'Q_matrix shape: {Q_matrix.shape}')

    # only one operator here
    operators = [SimpleMutation(), SimpleBinaryMutation()]

    # random policy
    policy = RandomPolicy(operators)

    # use of loaded data from UBQP instance
    evaluator = UBQPEvaluator(data={'Q': Q_matrix})

    # passing global evaluation param from ILS
    hcfi = HillClimberFirstImprovment(init, evaluator, operators, policy, validator, maximise=True, verbose=True)
    algo = ILS(init, evaluator, operators, policy, validator, localSearch=hcfi, maximise=True, verbose=True)

    # run the algorithm
    bestSol = algo.run(10000, ls_evaluations=100)

    print('Solution for UBQP instance score is {}'.format(evaluator.compute(bestSol)))


UBQP problem solving is now possible with **Macop**. As a reminder, the complete code is available in the ubqpExample.py_ file.

.. _ubqpExample.py: https://github.com/jbuisine/macop/blob/master/examples/ubqpExample.py
.. _documentation: https://jbuisine.github.io/macop/_build/html/documentations


.. _macop.solutions.discrete.BinarySolution: macop/macop.solutions.discrete.html#macop.solutions.discrete.BinarySolution
.. _macop.evaluators.base.Evaluator: macop/macop.evaluators.base.html#macop.evaluators.base.Evaluator
.. _macop.evaluators.discrete.mono.UBQPEvaluator: macop/macop.evaluators.discrete.mono.html#macop.evaluators.discrete.mono.UBQPEvaluator===================
A tour of Macop
===================

.. image:: _static/logo_macop.png
   :width: 300 px
   :align: center

This documentation will allow a user who wishes to use the **Macop** optimisation package to understand both how it works and offers examples of how to implement specific needs.

It will gradually take up the major ideas developed within **Macop** to allow for quick development. You can navigate directly via the menu available below to access a specific part of the documentation.


Introduction
================

`Macop` is a python package for solving discrete optimisation problems in nature. Continuous optimisation is also applicable but not yet developed. The objective is to allow a user to exploit the basic structure proposed by this package to solve a problem specific to him. The interest is that he can quickly abstract himself from the complications related to the way of evaluating, comparing, saving the progress of the search for good solutions but rather concentrate if necessary on his own algorithm. Indeed, `Macop` offers the following main and basic features: 

- **solutions:** representation of the solution;
- **validator:** such as constraint programming, a `validator` is a function which is used to validate or not a solution data state;
- **evaluator:** stores problem instance data and implements a `compute` method in order to evaluate a solution;
- **operators:** mutators, crossovers operators to update and obtain new solution;
- **policies:** the way you choose the available operators (might be using reinforcement learning);
- **algorithms:** generic and implemented optimisation research algorithms;
- **callbacks:** callbacks to automatically keep track of the search space advancement and restart from previous state if nedded.

.. image:: _static/documentation/macop_behaviour.png
   :width: 50 %
   :align: center

Based on all of these generic and/or implemented functionalities, the user will be able to quickly develop a solution to his problem while retaining the possibility of remaining in control of his development by overloading existing functionalities if necessary.

Problem instance
===================

In this tutorial, we introduce the way of using **Macop** and running your algorithm quickly using the well known `knapsack` problem.

Problem definition
~~~~~~~~~~~~~~~~~~~~~~

The **knapsack problem** is a problem in combinatorial optimisation: Given a set of items, each with a weight and a value, determine the number of each item to include in a collection so that the total weight is less than or equal to a given limit and the total value is as large as possible.


The image below provides an illustration of the problem:

.. image:: _static/documentation/knapsack_problem.png
   :width: 40 %
   :align: center


In this problem, we try to optimise the value associated with the objects we wish to put in our backpack while respecting the capacity of the bag (weight constraint).

.. warning::
    It is a combinatorial and therefore discrete problem. **Macop** decomposes its package into two parts, which is related to discrete optimisation on the one hand, and continuous optimisation on the other hand. This will be detailed later.


Problem implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

During the whole tutorial, the example used is based on the previous illustration with:

.. image:: _static/documentation/project_knapsack_problem.png
   :width: 85 %
   :align: center


Hence, we now define our problem in Python:

- worth value of each objects 
- weight associated to each of these objects

.. code-block:: python
    
    """
    Problem instance definition
    """

    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

Once we have defined the instance of our problem, we will need to define the representation of a solution to that problem.

Let's define the ``SimpleBinaryCrossover`` operator, allows to randomly change a binary value of our current solution.

Solutions
=============

Representing a solution to a specific problem is very important in an optimisation process. In this example, we will always use the **knapsack problem** as a basis.

In a first step, the management of the solutions by the macop package will be presented. Then a specific implementation for the current problem will be detailed.

Generic Solution
~~~~~~~~~~~~~~~~~~~~~~~~~

Inside macop.solutions.base_ module of `Macop`, the ``Solution`` class is available. It's an abstract solution class structure which:

- stores the solution data representation into its ``data`` attribute
- get ``size`` (shape) of specific data representation
- stores the ``score`` of the solution once a solution is evaluated

Some specific methods are available:

.. caution::
    An important thing here are the ``fitness``, ``size`` and ``data`` functions brought as an editable attribute by the ``@property`` and ``@XXXXX.setter`` decorators. The idea is to allow the user to modify these functions in order to change the expected result of the algorithm regardless of the data to be returned/modified. 

From these basic methods, it is possible to manage a representation of a solution to our problem. 

Allowing to initialise it randomly or not (using constructor or ``random`` method), to evaluate it (``evaluate`` method) and to check some constraints of validation of the solution (``isValid`` method).

.. note::
    Only one of these methods needs specification if we create our own type of solution. This is the ``random`` method, which depends on the need of the problem.

We will now see how to define a type of solution specific to our problem.

Solution representation for knapsack
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will now use the abstract ``Solution`` type available in the macop.solutions.base_ module in order to define our own solution.
First of all, let's look at the representation of our knapsack problem. **How to represent the solution?**

Knapsack solution
************************

A valid solution can be shown below where the sum of the object weights is 15 and the sum of the selected objects values is 8 (its fitness):

.. image:: _static/documentation/project_knapsack_solution.png
   :width:  85 %
   :align: center

Its representation can be translate as a **binary array** with value:

.. code-block::

    [1, 1, 0, 0, 1]

where selected objects have **1** as value otherwise **0**.

Binary Solution
**********************

We will now define our own type of solution by inheriting from macop.solutions.base.Solution_, which we will call ``BinarySolution``.

First we will define our new class as inheriting functionality from ``Solution`` (such as child class). 
We will also have to implement the ``random`` method to create a new random solution.

.. code-block:: python

    """
    modules imports
    """
    from macop.solutions.base import Solution
    import numpy as np

    class BinarySolution(Solution):
        
        @staticmethod
        def random(size, validator=None):

            # create binary array of specific size using numpy random module
            data = np.random.randint(2, size=size)
            # initialise new solution using constructor
            solution = BinarySolution(data, size)

            # check if validator is set
            if not validator:
                return solution

            # try to generate solution until solution validity (if validator is provided)
            while not validator(solution):
                data = np.random.randint(2, size=size)
                solution = BinarySolution(data, size)

            return solution

.. note::
    The current developed ``BinarySolution`` is available into macop.solutions.discrete.BinarySolution_ in **Macop**.

Using this new Solution representation, we can now generate solution randomly:

.. code-block:: python

    solution = BinarySolution.random(5)

In the next part, we will see how to verify that a solution meets certain modeling constraints of the problem.

Validate a solution
======================

When an optimisation problem requires respecting certain constraints, Macop allows you to quickly verify that a solution is valid. 
It is based on a defined function taking a solution as input and returning the validity criterion (true or false).

Validator definition
~~~~~~~~~~~~~~~~~~~~~~~~~

An invalid solution can be shown below where the sum of the object weights is greater than 15:

.. image:: _static/documentation/project_knapsack_invalid.png
   :width:  85 %
   :align: center

In fact, **[1, 0, 1, 0, 0]** is an invalid solution as we have a weight of **16** which violates the knapsack capacity constraint.

To avoid taking into account invalid solutions, we can define our function which will validate or not a solution based on our problem instance:

.. code-block:: python

    """
    Problem instance definition
    """

    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    """
    Validator function definition
    """
    def validator(solution):

        weight_sum = 0

        for i, w in enumerate(elements_weight):
            # add weight if current object is set to 1
            weight_sum += w * solution.data[i]
        
        # validation condition
        return weight_sum <= 15


Use of validator
~~~~~~~~~~~~~~~~~~~~~

We can now generate solutions randomly by passing our validation function as a parameter:

.. code-block:: python

    """
    Problem instance definition
    """
    ...
    
    """
    Validator function definition
    """
    ...

    # ensure valid solution
    solution = BinarySolution.random(5, validator)


.. caution::
    If the search space for valid solutions is very small compared to the overall search space, this can involve a considerable time for validating the solution and therefore obtaining a solution.

The validation of a solution is therefore now possible. In the next part we will focus on the evaluation of a solution.

Use of evaluators
====================

Now that it is possible to generate a solution randomly or not. It is important to know the value associated with this solution. We will then speak of evaluation of the solution. With the score associated with it, the `fitness`.

Generic evaluator
~~~~~~~~~~~~~~~~~~~~~~

As for the management of solutions, a generic evaluator class macop.evaluators.base.Evaluator_ is developed within **Macop**:

Abstract Evaluator class is used for computing fitness score associated to a solution. To evaluate all the solutions, this class:

- stores into its ``data`` initialiser dictionary attritute required measures when computing a solution
- has a ``compute`` abstract method enable to compute and associate a score to a given solution
- stores into its ``algo`` attritute the current algorithm to use (we will talk about algorithm later)

We must therefore now create our own evaluator based on the proposed structure.

Custom evaluator
~~~~~~~~~~~~~~~~~~~~~

To create our own evaluator, we need both:

- data useful for evaluating a solution
- compute the fitness associated with the state of the solution from these data. Hence, implement specific ``compute`` method.

We will define the ``KnapsackEvaluator`` class, which will therefore allow us to evaluate solutions to our current problem.

.. code-block:: python

    """
    modules imports
    """
    from macop.evaluators.base import Evaluator

    class KnapsackEvaluator(Evaluator):
        
        def compute(self, solution):

            # `_data` contains worths array values of objects
            fitness = 0
            for index, elem in enumerate(solution.data):
                fitness += self._data['worths'][index] * elem

            return fitness


It is now possible to initialise our new evaluator with specific data of our problem instance:

.. code-block:: python

    """
    Problem instance definition
    """
    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    """
    Evaluator problem instance
    """
    evaluator = KnapsackEvaluator(data={'worths': elements_score})

    # using defined BinarySolution
    solution = BinarySolution.random(5)

    # obtaining current solution score
    solution_fitness = solution.evaluate(evaluator)

    # score is also stored into solution
    solution_fitness = solution.fitness

.. note::
    The current developed ``KnapsackEvaluator`` is available into macop.evaluators.discrete.mono.KnapsackEvaluator_ in **Macop**.

In the next part we will see how to modify our current solution with the use of modification operator.

Apply operators to solution
==============================

Applying an operator to a solution consists of modifying the current state of the solution in order to obtain a new one. The goal is to find a better solution in the search space.

Operators definition
~~~~~~~~~~~~~~~~~~~~~~~~~

In the discrete optimisation literature, we can categorise operators into two sections:

- **mutators**: modification of one or more elements of a solution from its current state.
- **crossovers**: Inspired by Darwin's theory of evolution, we are going here from two solutions to generate a so-called offspring solution composed of the fusion of the data of the parent solutions.

Inside **Macop**, operators are also decomposed into these two categories. Inside macop.operators.base_, generic class ``Operator`` enables to manage any kind of operator.

Like the evaluator, the operator keeps **track of the algorithm** (using ``setAlgo`` method) to which he will be linked. This will allow better management of the way in which the operator must take into account the state of current data relating to the evolution of research.

``Mutation`` and ``Crossover`` classes inherite from ``Operator``. An ``apply`` function is required for any new operator.

We will now detail these categories of operators and suggest some relative to our problem.

Mutator operator
~~~~~~~~~~~~~~~~~~~~~

As detailed, the mutation operator consists in having a minimum impact on the current state of our solution. Here is an example of a modification that could be done for our problem.

.. image:: _static/documentation/project_knapsack_mutator.png
   :width:  90 %
   :align: center

In this example we change a bit value randomly and obtain a new solution from our search space.

.. warning::
    Applying an operator can conduct to a new but invalid solution from the search space.

The modification applied here is just a bit swapped. Let's define the ``SimpleBinaryMutation`` operator, allows to randomly change a binary value of our current solution.


.. code-block:: python

    """
    modules imports
    """
    from macop.operators.discrete.base import Mutation

    class SimpleBinaryMutation(Mutation):

        def apply(self, solution):
            
            # obtain targeted cell using solution size
            size = solution.size
            cell = random.randint(0, size - 1)

            # copy of solution
            copy_solution = solution.clone()

            # swicth values
            if copy_solution.data[cell]:
                copy_solution.data[cell] = 0
            else:
                copy_solution.data[cell] = 1

            # return the new obtained solution
            return copy_solution

We can now instanciate our new operator in order to obtain a new solution:


.. code-block:: python

    """
    BinaryMutator instance
    """
    mutator = SimpleBinaryMutation()

    # using defined BinarySolution
    solution = BinarySolution.random(5)

    # obtaining new solution using operator
    new_solution = mutator.apply(solution)


.. note::
    The developed ``SimpleBinaryMutation`` is available into macop.operators.discrete.mutators.SimpleBinaryMutation_ in **Macop**.


Crossover operator
~~~~~~~~~~~~~~~~~~~~~~~


Inspired by Darwin's theory of evolution, crossover starts from two solutions to generate a so-called offspring solution composed of the fusion of the data of the parent solutions.

.. image:: _static/documentation/project_knapsack_crossover.png
   :width:  95%
   :align: center

In this example we merge two solutions with a specific splitting criterion in order to obtain an offspring.

We will now implement the SimpleCrossover crossover operator, which will merge data from two solutions. 
The first half of solution 1 will be saved and added to the second half of solution 2 to generate the new solution (offspring).


.. code-block:: python

    """
    modules imports
    """
    from macop.operators.discrete.base import Crossover

    class SimpleCrossover(Crossover):

        def apply(self, solution1, solution2):
            
            size = solution1.size

            # default split index used
            splitIndex = int(size / 2)

            # copy data of solution 1
            firstData = solution1.data.copy()

            # copy of solution 2
            copy_solution = solution2.clone()

            copy_solution.data[splitIndex:] = firstData[splitIndex:]

            return copy_solution


We can now use the crossover operator created to generate new solutions. Here is an example of use:

.. code-block:: python

    """
    SimpleCrossover instance
    """
    crossover = SimpleCrossover()

    # using defined BinarySolution
    solution1 = BinarySolution.random(5)
    solution2 = BinarySolution.random(5)

    # obtaining new solution using crossover
    offspring = crossover.apply(solution1, solution2)

.. tip::
    The developed ``SimpleCrossover`` is available into macop.operators.discrete.crossovers.SimpleCrossover_ in **Macop**.
    However, the choice of halves of the merged data is made randomly.

Next part introduce the ``policy`` feature of **Macop** which enables to choose the next operator to apply during the search process based on specific criterion.

Operator choices
===================

The ``policy`` feature of **Macop** enables to choose the next operator to apply during the search process of the algorithm based on specific criterion.

Why using policy ?
~~~~~~~~~~~~~~~~~~~~~~~

Sometimes the nature of the problem and its instance can strongly influence the search results when using mutation operators or crossovers. 
Automated operator choice strategies have also been developed in the literature, notably based on reinforcement learning.

The operator choice problem can be seen as the desire to find the best solution generation operator at the next evaluation that will be the most conducive to precisely improving the solution.

.. image:: _static/documentation/operators_choice.png
   :width:  45 %
   :align: center

.. note::
    An implementation using reinforcement learning has been developed as an example in the macop.policies.reinforcement_ module. 
    However, it will not be detailed here. You can refer to the API documentation for more details.


Custom policy
~~~~~~~~~~~~~~~~~~

In our case, we are not going to exploit a complex enough implementation of a ``policy``. Simply, we will use a random choice of operator.

First, let's take a look of the ``Policy`` abstract class available in macop.policies.base_:


``Policy`` instance will have of ``operators`` attributes in order to keep track of possible operators when selecting one. 
Here, in our implementation we only need to specify the ``select`` abstract method. The ``apply`` method will select the next operator and return the new solution.

.. code-block:: python

    """
    module imports
    """
    from macop.policies.base import Policy

    class RandomPolicy(Policy):

        def select(self):
            """
            Select specific operator
            """
            # choose operator randomly
            index = random.randint(0, len(self.operators) - 1)
            return self.operators[index]


We can now use this operator choice policy to update our current solution:


.. code-block:: python

    """
    Operators instances
    """
    mutator = SimpleMutation()
    crossover = SimpleCrossover()

    """
    RandomPolicy instance
    """
    policy = RandomPolicy([mutator, crossover])

    """
    Current solutions instance
    """
    solution1 = BinarySolution.random(5)
    solution2 = BinarySolution.random(5)

    # pass two solutions in parameters in case of selected crossover operator
    new_solution = policy.apply(solution1, solution2)

.. caution::
    By default if ``solution2`` parameter is not provided into ``policy.apply`` method for crossover, the best solution known is used from the algorithm linked to the ``policy``.

Updating solutions is therefore now possible with our policy. It is high time to dive into the process of optimizing solutions and digging into our research space.

Optimisation process
=======================

Let us now tackle the interesting part concerning the search for optimum solutions in our research space.

Find local and global optima
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Overall, in an optimization process, we will seek to find the best, or the best solutions that minimize or maximize our objective function (fitness score obtained) in order to respond to our problem.

.. image:: _static/documentation/search_space.png
   :width:  95 %
   :align: center

Sometimes, the search space can be very simple. A local search can provide access to the global optimum as shown in figure (a) above. 
In other cases, the search space is more complex. It may be necessary to explore more rather than exploit in order to get out of a convex zone and not find the global optimum but only a local opmatime solution. 
This problem is illustrated in figure (b).

Abstract algorithm class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An abstract class is proposed within Macop to generalize the management of an algorithm and therefore of a heuristic. 
It is located in the macop.algorithms.base_ module. 

We will pay attention to the different methods of which she is composed. This class enables to manage some common usages of operation research algorithms:

- initialization function of solution
- validator function to check if solution is valid or not (based on some criteria)
- evaluation function to give fitness score to a solution
- operators used in order to update solution during search process
- policy process applied when choosing next operator to apply
- callbacks function in order to do some relative stuff every number of evaluation or reload algorithm state
- parent algorithm associated to this new algorithm instance (hierarchy management)

She is composed of few default attributes:

- initialiser: {function} -- basic function strategy to initialise solution
- evaluator: {:class:`~macop.evaluators.base.Evaluator`} -- evaluator instance in order to obtained fitness (mono or multiple objectives)
- operators: {[:class:`~macop.operators.base.Operator`]} -- list of operator to use when launching algorithm
- policy: {:class:`~macop.policies.base.Policy`} -- Policy instance strategy to select operators
- validator: {function} -- basic function to check if solution is valid or not under some constraints
- maximise: {bool} -- specify kind of optimisation problem 
- verbose: {bool} -- verbose or not information about the algorithm
- currentSolution: {:class:`~macop.solutions.base.Solution`} -- current solution managed for current evaluation comparison
- bestSolution: {:class:`~macop.solutions.base.Solution`} -- best solution found so far during running algorithm
- callbacks: {[:class:`~macop.callbacks.base.Callback`]} -- list of Callback class implementation to do some instructions every number of evaluations and `load` when initialising algorithm
- parent: {:class:`~macop.algorithms.base.Algorithm`} -- parent algorithm reference in case of inner Algorithm instance (optional)

.. caution::
    An important thing here are the ``result`` functions brought as an editable attribute by the ``@property`` and ``@result.setter`` decorators. The idea is to allow the user to modify these functions in order to change the expected result of the algorithm regardless of the data to be returned/modified. 

The notion of hierarchy between algorithms is introduced here. We can indeed have certain dependencies between algorithms. 
The methods ``increaseEvaluation``, ``getGlobalEvaluation`` and ``getGlobalMaxEvaluation`` ensure that the expected global number of evaluations is correctly managed, just like the ``stop`` method for the search stop criterion.

The ``evaluate``, ``update`` and ``isBetter`` will be used a lot when looking for a solution in the search space. 
In particular the ``update`` function, which will call the ``policy`` instance to generate a new valid solution.
``isBetter`` method is also overloadable especially if the algorithm does not take any more into account than a single solution to be verified (verification via a population for example).

The ``initRun`` method specify the way you intialise your algorithm (``bestSolution`` and ``currentSolution`` as example) if algorithm not already initialised.

.. note:: 
    The ``initRun`` method can also be used for intialise population of solutions instead of only one best solution, if you want to manage a genetic algorithm.

Most important part is the ``run`` method. Into abstract, the ``run`` method only initialised the current number of evaluation for the algorithm based on the parent algorithm if we are into inner algorithm.
It is always **mandatory** to call the parent class ``run`` method using ``super().run(evaluations)``. Then, using ``evaluations`` parameter which is the number of evaluations budget to run, we can process or continue to find solutions into search space.

.. warning::
    The other methods such as ``addCallback``, ``resume`` and ``progress`` will be detailed in the next part focusing on the notion of callback.


Local search algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are going to carry out our first local search algorithm within our search space. A `local search` consists of starting from a solution, then applying a mutation or crossover operation to it, in order to obtain a new one. 
This new solution is evaluated and retained if it is better. We will speak here of the notion of **neighborhood exploration**. The process is then completed in the same way. 
The local search ends after a certain number of evaluations and the best evaluated solution obtained is returned.

Let's implement an algorithm well known under the name of hill climber best improvment inheriting from the mother algorithm class and name it ``HillClimberBestImprovment``.


.. code-block:: python

    """
    module imports
    """
    from macop.algorithms.base import Algorithm

    class HillClimberBestImprovment(Algorithm):

        def run(self, evaluations):
            """
            Run a local search algorithm
            """

            # by default use of mother method to initialise variables
            super().run(evaluations)

            # initialise current solution and best solution
            self.initRun()

            solutionSize = self._currentSolution.size

            # local search algorithm implementation
            while not self.stop():

                for _ in range(solutionSize):

                    # update current solution using policy
                    newSolution = self.update(self._currentSolution)

                    # if better solution than currently, replace it
                    if self.isBetter(newSolution):
                        self._bestSolution = newSolution

                    # increase number of evaluations
                    self.increaseEvaluation()

                    # stop algorithm if necessary
                    if self.stop():
                        break

                # set new current solution using best solution found in this neighbor search
                self._currentSolution = self._bestSolution
            
            return self._bestSolution

Our algorithm is now ready to work. As previously, let us define two operators as well as a random choice strategy. 
We will also need to define a **solution initialisation function** so that the algorithm can generate new solutions.


.. code-block:: python

    """
    Problem instance definition
    """
    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    # evaluator instance
    evaluator = KnapsackEvaluator(data={'worths': elements_score})

    # valid instance using lambda
    validator = lambda solution: sum([ elements_weight[i] * solution.data[i] for i in range(len(solution.data))]) <= 15
    
    # initialiser instance using lambda with default param value
    initialiser = lambda x=5: BinarySolution.random(x, validator)
    
    # operators list with crossover and mutation
    operators = [SimpleCrossover(), SimpleMutation()]
    
    # policy random instance
    policy = RandomPolicy(operators)
    
    # maximizing algorithm (relative to knapsack problem)
    algo = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)

    # run the algorithm and get solution found
    solution = algo.run(100)
    print(solution.fitness)


.. note::
    The ``verbose`` algorithm parameter will log into console the advancement process of the algorithm is set to ``True`` (the default value).

Exploratory algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~

As explained in **figure (b)** of **section 8.1**, sometimes the search space is more complicated due to convex parts and need heuristic with other strategy rather than a simple local search.

The way to counter this problem is to allow the algorithm to exit the exploitation phase offered by local search. But rather to seek to explore other parts of the research space. This is possible by simply carrying out several local searches with our budget (number of evaluations).

The idea is to make a leap in the search space in order to find a new local optimum which can be the global optimum. The explained process is illustrated below:

.. image:: _static/documentation/search_space_simple.png
   :width:  45 %
   :align: center


We are going to implement a more specific algorithm allowing to take a new parameter as input. This is a local search, the one previously developed. For that, we will have to modify the constructor a little.
Let's called this new algorithm ``IteratedLocalSearch``:

.. code-block:: python

    """
    module imports
    """
    from macop.algorithms.base import Algorithm

    class IteratedLocalSearch(Algorithm):
        
        def __init__(self,
                    initialiser,
                    evaluator,
                    operators,
                    policy,
                    validator,
                    localSearch,
                    maximise=True,
                    parent=None,
                    verbose=True):
            
            super().__init__(initialiser, evaluator, operators, policy, validator, maximise, parent, verbose)

            # specific local search associated with current algorithm
            self._localSearch = localSearch

            # need to attach current algorithm as parent
            self._localSearch.setParent(self)


        def run(self, evaluations, ls_evaluations=100):
            """
            Run the iterated local search algorithm using local search
            """

            # by default use of mother method to initialise variables
            super().run(evaluations)

            # initialise current solution
            self.initRun()

            # local search algorithm implementation
            while not self.stop():

                # create and search solution from local search (stop method can be called inside local search)
                newSolution = self._localSearch.run(ls_evaluations)

                # if better solution than currently, replace it
                if self.isBetter(newSolution):
                    self._bestSolution = newSolution

                self.information()

            return self._bestSolution

In the initialization phase we have attached our local search passed as a parameter with the current algorithm as parent. 
The goal is to touch keep track of the overall search evaluation number (relative to the parent algorithm).

Then, we use this local search in our ``run`` method to allow a better search for solutions.

.. code-block:: python

    """
    Problem instance definition
    """
    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    # evaluator instance
    evaluator = KnapsackEvaluator(data={'worths': elements_score})

    # valid instance using lambda
    validator = lambda solution: sum([ elements_weight[i] * solution.data[i] for i in range(len(solution.data))]) <= 15
    
    # initialiser instance using lambda with default param value
    initialiser = lambda x=5: BinarySolution.random(x, validator)
    
    # operators list with crossover and mutation
    operators = [SimpleCrossover(), SimpleMutation()]
    
    # policy random instance
    policy = RandomPolicy(operators)
    
    # maximizing algorithm (relative to knapsack problem)
    localSearch = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)
    algo = IteratedLocalSearch(initialiser, evaluator, operators, policy, validator, localSearch=local_search, maximise=True, verbose=False)

    # run the algorithm using local search and get solution found 
    solution = algo.run(evaluations=100, ls_evaluations=10)
    print(solution.fitness)


.. note:: 
    These two last algorithms developed are available in the library within the module macop.algorithms.mono_.

We have one final feature to explore in the next part. This is the notion of ``callback``.

Keep track
==============

Keeping track of the running algorithm can be useful on two levels. First of all to understand how it unfolded at the end of the classic run. But also in the case of the unwanted shutdown of the algorithm. 
This section will allow you to introduce the recovery of the algorithm thanks to a continuous backup functionality.

Logging into algorithm
~~~~~~~~~~~~~~~~~~~~~~

Some logs can be retrieve after running an algorithm. **Macop** uses the ``logging`` Python package in order to log algorithm advancement.

Here is an example of use when running an algorithm:

.. code-block:: python

    """
    basic imports
    """
    import logging

    # logging configuration
    logging.basicConfig(format='%(asctime)s %(message)s', filename='data/example.log', level=logging.DEBUG)

    ...
    
    # maximizing algorithm (relative to knapsack problem)
    algo = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)

    # run the algorithm using local search and get solution found 
    solution = algo.run(evaluations=100)
    print(solution.fitness)

Hence, log data are saved into ``data/example.log`` in our example.

Callbacks introduction
~~~~~~~~~~~~~~~~~~~~~~~

Having output logs can help to understand an error that has occurred, however all the progress of the research carried out may be lost. 
For this, the functionality relating to callbacks has been developed.

Within **Macop**, a callback is a specific instance of macop.callbacks.base.Callback_ that allows you to perform an action of tracing / saving information **every** ``n`` **evaluations** but also reloading information if necessary when restarting an algorithm.


- The ``run`` method will be called during run process of the algo and do backup at each specific number of evaluations. 
- The ``load`` method will be used to reload the state of the algorithm from the last information saved. All saved data is saved in a file whose name will be specified by the user.

Towards the use of Callbacks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are going to create our own Callback instance called ``BasicCheckpoint`` which will save the best solution found and number of evaluations done in order to reload it for the next run of our algorithm.

.. code-block:: python

    """
    module imports
    """
    from macop.callbacks.base import Callback


    class BasicCheckpoint(Callback):
        
        def run(self):
            """
            Check if necessary to do backup based on `every` variable
            """
            # get current best solution
            solution = self.algo._bestSolution

            currentEvaluation = self.algo.getGlobalEvaluation()

            # backup if necessary every number of evaluations
            if currentEvaluation % self._every == 0:

                # create specific line with solution data
                solution.data = ""
                solutionSize = len(solution.data)

                for index, val in enumerate(solution.data):
                    solution.data += str(val)

                    if index < solutionSize - 1:
                        solution.data += ' '

                # number of evaluations done, solution data and fitness score
                line = str(currentEvaluation) + ';' + solution.data + ';' + str(
                    solution.fitness) + ';\n'

                # check if file exists
                if not os.path.exists(self._filepath):
                    with open(self._filepath, 'w') as f:
                        f.write(line)
                else:
                    with open(self._filepath, 'a') as f:
                        f.write(line)

        def load(self):
            """
            Load last backup line and set algorithm state (best solution and evaluations)
            """
            if os.path.exists(self._filepath):

                with open(self._filepath) as f:

                    # get last line and read data
                    lastline = f.readlines()[-1]
                    data = lastline.split(';')

                    # get evaluation  information
                    globalEvaluation = int(data[0])

                    # restore number of evaluations
                    if self.algo.getParent() is not None:
                        self.algo.getParent()._numberOfEvaluations = globalEvaluation
                    else:
                        self.algo._numberOfEvaluations = globalEvaluation

                    # get best solution data information
                    solution.data = list(map(int, data[1].split(' ')))

                    # avoid uninitialised solution
                    if self.algo._bestSolution is None:
                        self.algo._bestSolution = self.algo.initialiser()

                    # set to algorithm the lastest obtained best solution
                    self.algo._bestsolution.data = np.array(solution.data)
                    self.algo._bestSolution._score = float(data[2])


In this way, it is possible to specify the use of a callback to our algorithm instance:


.. code-block:: python

    ...
    
    # maximizing algorithm (relative to knapsack problem)
    algo = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)

    callback = BasicCheckpoint(every=5, filepath='data/hillClimberBackup.csv')

    # add callback into callback list
    algo.addCallback(callback)

    # run the algorithm using local search and get solution found 
    solution = algo.run(evaluations=100)
    print(solution.fitness)


.. note::
    It is possible to add as many callbacks as desired in the algorithm in question.


Previously, some methods of the abstract ``Algorithm`` class have not been presented. These methods are linked to the use of callbacks, 
in particular the ``addCallback`` method which allows the addition of a callback to an algorithm instance as seen above.

- The ``resume`` method will reload all callbacks list using ``load`` method.
- The ``progress`` method will ``run`` each callbacks during the algorithm search.

If we want to exploit this functionality, then we will need to exploit them within our algorithm. Let's make the necessary modifications for our algorithm ``IteratedLocalSearch``:


.. code-block:: python

    """
    module imports
    """
    from macop.algorithms.base import Algorithm

    class IteratedLocalSearch(Algorithm):
        
        ...

        def run(self, evaluations, ls_evaluations=100):
            """
            Run the iterated local search algorithm using local search
            """

            # by default use of mother method to initialise variables
            super().run(evaluations)

            # initialise current solution
            self.initRun()

            # restart using callbacks backup list
            self.resume()

            # local search algorithm implementation
            while not self.stop():

                # create and search solution from local search
                newSolution = self._localSearch.run(ls_evaluations)

                # if better solution than currently, replace it
                if self.isBetter(newSolution):
                    self._bestSolution = newSolution

                # check if necessary to call each callbacks
                self.progress()

                self.information()

            return self._bestSolution


All the features of **Macop** were presented. The next section will aim to quickly present the few implementations proposed within **Macop** to highlight the modulality of the package.


Advanced usages
=======================

Multi-objective discrete optimisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Within the API of **Macop**, you can find an implementation of The Multi-objective evolutionary algorithm based on decomposition (MOEA/D) is a general-purpose algorithm for approximating the Pareto set of multi-objective optimization problems. 
It decomposes the original multi-objective problem into a number of single-objective optimization sub-problems and then uses an evolutionary process to optimize these sub-problems simultaneously and cooperatively. 
MOEA/D is a state-of-art algorithm in aggregation-based approaches for multi-objective optimization.

.. image:: _static/documentation/search_space_moead.png
   :width:  45 %
   :align: center


As illustrated below, the two main objectives are sub-divised into 5 single-objective optimization sub-problems in order to find the Pareto front.

- macop.algorithms.multi.MOSubProblem_ class defines each sub-problem of MOEA/D.
- macop.algorithms.multi.MOEAD_ class exploits ``MOSubProblem`` and implements MOEA/D using weighted-sum of objectives method.

An example with MOEAD for knapsack problem is available in knapsackMultiExample.py_.

.. _knapsackMultiExample.py: https://github.com/jbuisine/macop/blob/master/examples/knapsackMultiExample.py

Continuous Zdt problems
~~~~~~~~~~~~~~~~~~~~~~~

Even if the package is not primarily intended for continuous optimisation, it allows for adaptation to continuous optimisation. 

Based on the Zdt_ benchmarks function, it offers an implementation of Solution, Operator and Evaluator to enable the optimisation of this kind of problem.

.. _Zdt: https://en.wikipedia.org/wiki/Test_functions_for_optimization

- macop.solutions.continuous.ContinuousSolution_: manage float array solution in order to represent continuous solution;
- macop.operators.continuous.mutators.PolynomialMutation_: update solution using polynomial mutation over solution's data;
- macop.operators.continuous.crossovers.BasicDifferentialEvolutionCrossover_: use of new generated solutions in order to obtain new offspring solution;
- macop.evaluators.continous.mono.ZdtEvaluator_: continuous evaluator for `Zdt` problem instance. Take into its ``data``, the ``f`` Zdt function;
- macop.callbacks.classicals.ContinuousCallback_: manage callback and backup of continuous solution.

A complete implementation example with the Rosenbrock_ function is available.

.. _macop.algorithms.base: macop/macop.algorithms.base.html#module-macop.algorithms.base
.. _macop.algorithms.mono: macop/macop.algorithms.mono.html#module-macop.algorithms.mono

.. _macop.solutions.base: macop/macop.solutions.base.html#module-macop.solutions.base
.. _macop.solutions.base.Solution: macop/macop.solutions.base.html#macop.solutions.base.Solution
.. _macop.solutions.discrete.BinarySolution: macop/macop.solutions.discrete.html#macop.solutions.discrete.BinarySolution

.. _macop.evaluators.base.Evaluator: macop/macop.evaluators.base.html#macop.evaluators.base.Evaluator
.. _macop.evaluators.discrete.mono.KnapsackEvaluator: macop/macop.evaluators.discrete.mono.html#macop.evaluators.discrete.mono.KnapsackEvaluator

.. _macop.operators.base: macop/macop.operators.base.html#module-macop.operators.base
.. _macop.operators.discrete.mutators.SimpleBinaryMutation: macop/macop.operators.discrete.mutators.html#macop.operators.discrete.mutators.SimpleBinaryMutation
.. _macop.operators.discrete.crossovers.SimpleCrossover: macop/macop.operators.discrete.crossovers.html#macop.operators.discrete.crossovers.SimpleCrossover

.. _macop.policies.reinforcement: macop/macop.policies.reinforcement.html#module-macop.policies.reinforcement
.. _macop.policies.base: macop/macop.policies.base.html#module-macop.policies.base

.. _macop.callbacks.base.Callback: macop/macop.callbacks.base.html#macop.callbacks.base.Callback

.. _macop.algorithms.multi.MOSubProblem: macop/macop.algorithms.multi.html#macop.algorithms.multi.MOSubProblem
.. _macop.algorithms.multi.MOEAD: macop/macop.algorithms.multi.html#macop.algorithms.multi.MOEAD

.. _macop.solutions.continuous.ContinuousSolution: macop/macop.solutions.continuous.html#macop.solutions.continuous.ContinuousSolution

.. _macop.operators.continuous.mutators.PolynomialMutation: macop/macop.operators.continuous.mutators.html#macop.operators.continuous.mutators.PolynomialMutation
.. _macop.operators.continuous.crossovers.BasicDifferentialEvolutionCrossover: macop/macop.operators.continuous.crossovers.html#macop.operators.continuous.crossovers.BasicDifferentialEvolutionCrossover

.. _macop.evaluators.continous.mono.ZdtEvaluator: macop/macop.evaluators.continuous.mono.html#macop.evaluators.continous.mono.ZdtEvaluator

.. _macop.callbacks.classicals.ContinuousCallback: macop/macop.callbacks.classicals.html#macop.callbacks.classicals.ContinuousCallback


.. _Rosenbrock: https://github.com/jbuisine/macop/blob/master/examples/ZdtExample.pyAPI
=============

Modules description
~~~~~~~~~~~~~~~~~~~

**Macop** offers the following main and basic features: 

- **algorithms:** generic and implemented optimisation research algorithms;
- **callbacks:** callbacks to automatically keep track of the search space advancement and restart from previous state if nedded;
- **evaluator:** stores problem instance data and implements a `compute` method in order to evaluate a solution;
- **operators:** mutators, crossovers operators to update and obtain new solution;
- **policies:** the way you choose the available operators (might be using reinforcement learning);
- **solutions:** representation of the solution;
- **validator:** such as constraint programming, a `validator` is a function which is used to validate or not a solution data state;


Common and base modules
~~~~~~~~~~~~~~~~~~~~~~~

The modules presented in this section are common to all types of optimisation problems. The abstract classes proposed here form the basis of the package's structure.

macop.algorithms
-------------------

.. autosummary::
   :toctree: macop
   
   macop.algorithms.base
      
   macop.algorithms.mono
   macop.algorithms.multi

macop.callbacks
-------------------

.. autosummary::
   :toctree: macop
   
   macop.callbacks.base

   macop.callbacks.classicals
   macop.callbacks.multi
   macop.callbacks.policies

macop.evaluators
-------------------

.. autosummary::
   :toctree: macop
   
   macop.evaluators.base

macop.operators
-------------------

.. autosummary::
   :toctree: macop
   
   macop.operators.base

macop.policies
-------------------

.. autosummary::
   :toctree: macop
   
   macop.policies.base
      
   macop.policies.classicals
   macop.policies.reinforcement

macop.solutions
-------------------

.. autosummary::
   :toctree: macop

   macop.solutions.base

macop.utils
-------------------

.. autosummary::
   :toctree: macop

   macop.utils.progress


Discrete Optimisation
~~~~~~~~~~~~~~~~~~~~~

Some implementations of discrete optimisation problem functionalities are available. They can be used as example implementations or can simply be used by the user.

macop.evaluators
-------------------

.. autosummary::
   :toctree: macop
   
   macop.evaluators.discrete.mono
   macop.evaluators.discrete.multi

macop.operators
-------------------

.. autosummary::
   :toctree: macop
   
   macop.operators.discrete.mutators
   macop.operators.discrete.crossovers

macop.solutions
-------------------

.. autosummary::
   :toctree: macop

   macop.solutions.discrete


Continuous Optimisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although continuous optimization is not the priority of this package, the idea is to leave the possibility to any user to implement or even propose implementations related to this kind of problem. The modules are here for the moment nearly empty (only with Zdt functions example) but present to establish the structure relative to these types of implementations.

If a user wishes to propose these developments so that they can be added in a future version of the package, he can refer to the guidelines_ for contributions of the package.

.. _guidelines: https://github.com/prise-3d/macop/blob/master/CONTRIBUTING.md

macop.evaluators
-------------------

.. autosummary::
   :toctree: macop
   
   macop.evaluators.continuous.mono
   macop.evaluators.continuous.multi

macop.operators
-------------------

.. autosummary::
   :toctree: macop
   
   macop.operators.continuous.mutators
   macop.operators.continuous.crossovers

macop.solutions
-------------------

.. autosummary::
   :toctree: macop

   macop.solutions.continuousmacop.operators.continuous.mutators
===================================

.. automodule:: macop.operators.continuous.mutators

   
   
   

   
   
   

   
   
   macop.evaluators.continuous.mono
================================

.. automodule:: macop.evaluators.continuous.mono

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      ZdtEvaluator
   
   

   
   
   macop.solutions.discrete
========================

.. automodule:: macop.solutions.discrete

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      BinarySolution
      CombinatoryIntegerSolution
      IntegerSolution
   
   

   
   
   macop.algorithms.multi
======================

.. automodule:: macop.algorithms.multi

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      MOEAD
      MOSubProblem
   
   

   
   
   macop.callbacks.multi
=====================

.. automodule:: macop.callbacks.multi

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      MultiCheckpoint
      ParetoCheckpoint
   
   

   
   
   macop.evaluators.base
=====================

.. automodule:: macop.evaluators.base

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      Evaluator
   
   

   
   
   macop.callbacks.classicals
==========================

.. automodule:: macop.callbacks.classicals

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      BasicCheckpoint
   
   

   
   
   macop.policies.base
===================

.. automodule:: macop.policies.base

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      Policy
   
   

   
   
   macop.callbacks.base
====================

.. automodule:: macop.callbacks.base

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      Callback
   
   

   
   
   macop.callbacks.policies
========================

.. automodule:: macop.callbacks.policies

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      UCBCheckpoint
   
   

   
   
   macop.operators.continuous.crossovers
=====================================

.. automodule:: macop.operators.continuous.crossovers

   
   
   

   
   
   

   
   
   macop.evaluators.discrete.mono
==============================

.. automodule:: macop.evaluators.discrete.mono

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      KnapsackEvaluator
   
   

   
   
   macop.operators.discrete.crossovers
===================================

.. automodule:: macop.operators.discrete.crossovers

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      RandomSplitCrossover
      SimpleCrossover
   
   

   
   
   macop.algorithms.base
=====================

.. automodule:: macop.algorithms.base

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      Algorithm
   
   

   
   
   macop.solutions.continuous
==========================

.. automodule:: macop.solutions.continuous

   
   
   

   
   
   

   
   
   macop.policies.reinforcement
============================

.. automodule:: macop.policies.reinforcement

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      UCBPolicy
   
   

   
   
   macop.policies.classicals
=========================

.. automodule:: macop.policies.classicals

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      RandomPolicy
   
   

   
   
   macop.utils.progress
====================

.. automodule:: macop.utils.progress

   
   
   .. rubric:: Functions

   .. autosummary::
   
      macop_line
      macop_progress
      macop_text
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      Colors
   
   

   
   
   macop.operators.discrete.mutators
=================================

.. automodule:: macop.operators.discrete.mutators

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      SimpleBinaryMutation
      SimpleMutation
   
   

   
   
   macop.evaluators.discrete.multi
===============================

.. automodule:: macop.evaluators.discrete.multi

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      WeightedSum
   
   

   
   
   macop.operators.base
====================

.. automodule:: macop.operators.base

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      Crossover
      KindOperator
      Mutation
      Operator
   
   

   
   
   macop.evaluators.continuous
===========================

.. automodule:: macop.evaluators.continuous

   
   
   

   
   
   

   
   
   macop.solutions.base
====================

.. automodule:: macop.solutions.base

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      Solution
   
   

   
   
   macop.evaluators.continuous.multi
=================================

.. automodule:: macop.evaluators.continuous.multi

   
   
   

   
   
   

   
   
   macop.algorithms.mono
=====================

.. automodule:: macop.algorithms.mono

   
   
   

   
   
   .. rubric:: Classes

   .. autosummary::
   
      HillClimberBestImprovment
      HillClimberFirstImprovment
      IteratedLocalSearch
   
   

   
   
   Problem instance
===================

In this tutorial, we introduce the way of using **Macop** and running your algorithm quickly using the well known `knapsack` problem.

Problem definition
~~~~~~~~~~~~~~~~~~~~~~

The **knapsack problem** is a problem in combinatorial optimisation: Given a set of items, each with a weight and a value, determine the number of each item to include in a collection so that the total weight is less than or equal to a given limit and the total value is as large as possible.


The image below provides an illustration of the problem:

.. image:: ../_static/documentation/knapsack_problem.png
   :width: 40 %
   :align: center


In this problem, we try to optimise the value associated with the objects we wish to put in our backpack while respecting the capacity of the bag (weight constraint).

.. warning::
    It is a combinatorial and therefore discrete problem. **Macop** decomposes its package into two parts, which is related to discrete optimisation on the one hand, and continuous optimisation on the other hand. This will be detailed later.


Problem implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

During the whole tutorial, the example used is based on the previous illustration with:

.. image:: ../_static/documentation/project_knapsack_problem.png
   :width: 85 %
   :align: center


Hence, we now define our problem in Python:

- worth value of each objects 
- weight associated to each of these objects

.. code-block:: python
    
    """
    Problem instance definition
    """

    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

Once we have defined the instance of our problem, we will need to define the representation of a solution to that problem.

Let's define the ``SimpleBinaryCrossover`` operator, allows to randomly change a binary value of our current solution.Keep track
==============

Keeping track of the running algorithm can be useful on two levels. First of all to understand how it unfolded at the end of the classic run. But also in the case of the unwanted shutdown of the algorithm. 
This section will allow you to introduce the recovery of the algorithm thanks to a continuous backup functionality.

Logging into algorithm
~~~~~~~~~~~~~~~~~~~~~~

Some logs can be retrieve after running an algorithm. **Macop** uses the ``logging`` Python package in order to log algorithm advancement.

Here is an example of use when running an algorithm:

.. code-block:: python

    """
    basic imports
    """
    import logging

    # logging configuration
    logging.basicConfig(format='%(asctime)s %(message)s', filename='data/example.log', level=logging.DEBUG)

    ...
    
    # maximizing algorithm (relative to knapsack problem)
    algo = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)

    # run the algorithm using local search and get solution found 
    solution = algo.run(evaluations=100)
    print(solution.fitness)

Hence, log data are saved into ``data/example.log`` in our example.

Callbacks introduction
~~~~~~~~~~~~~~~~~~~~~~~

Having output logs can help to understand an error that has occurred, however all the progress of the research carried out may be lost. 
For this, the functionality relating to callbacks has been developed.

Within **Macop**, a callback is a specific instance of ``macop.callbacks.Callback`` that allows you to perform an action of tracing / saving information **every** ``n`` **evaluations** but also reloading information if necessary when restarting an algorithm.


.. code-block:: python

    class Callback():

        def __init__(self, every, filepath):
            ...

        @abstractmethod
        def run(self):
            """
            Check if necessary to do backup based on `every` variable
            """
            pass

        @abstractmethod
        def load(self):
            """
            Load last backup line of solution and set algorithm state at this backup
            """
            pass

        def setAlgo(self, algo):
            """
            Specify the main algorithm instance reference
            """
            ...


- The ``run`` method will be called during run process of the algo and do backup at each specific number of evaluations. 
- The ``load`` method will be used to reload the state of the algorithm from the last information saved. All saved data is saved in a file whose name will be specified by the user.

Towards the use of Callbacks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are going to create our own Callback instance called ``BasicCheckpoint`` which will save the best solution found and number of evaluations done in order to reload it for the next run of our algorithm.

.. code-block:: python

    """
    module imports
    """
    from macop.callbacks.base import Callback


    class BasicCheckpoint(Callback):
        
        def run(self):
            """
            Check if necessary to do backup based on `every` variable
            """
            # get current best solution
            solution = self.algo._bestSolution

            currentEvaluation = self.algo.getGlobalEvaluation()

            # backup if necessary every number of evaluations
            if currentEvaluation % self._every == 0:

                # create specific line with solution data
                solution.data = ""
                solutionSize = len(solution.getdata = ))

                for index, val in enumerate(solution.getdata = )):
                    solution.data += str(val)

                    if index < solutionSize - 1:
                        solution.data += ' '

                # number of evaluations done, solution data and fitness score
                line = str(currentEvaluation) + ';' + solution.data + ';' + str(
                    solution.fitness) + ';\n'

                # check if file exists
                if not os.path.exists(self._filepath):
                    with open(self._filepath, 'w') as f:
                        f.write(line)
                else:
                    with open(self._filepath, 'a') as f:
                        f.write(line)

        def load(self):
            """
            Load last backup line and set algorithm state (best solution and evaluations)
            """
            if os.path.exists(self._filepath):

                with open(self._filepath) as f:

                    # get last line and read data
                    lastline = f.readlines()[-1]
                    data = lastline.split(';')

                    # get evaluation  information
                    globalEvaluation = int(data[0])

                    # restore number of evaluations
                    if self.algo.getParent() is not None:
                        self.algo.getParent()._numberOfEvaluations = globalEvaluation
                    else:
                        self.algo._numberOfEvaluations = globalEvaluation

                    # get best solution data information
                    solution.data = list(map(int, data[1].split(' ')))

                    # avoid uninitialised solution
                    if self.algo._bestSolution is None:
                        self.algo._bestSolution = self.algo.initialiser()

                    # set to algorithm the lastest obtained best solution
                    self.algo._bestsolution.getdata = ) = np.array(solution.data)
                    self.algo._bestSolution._score = float(data[2])


In this way, it is possible to specify the use of a callback to our algorithm instance:


.. code-block:: python

    ...
    
    # maximizing algorithm (relative to knapsack problem)
    algo = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)

    callback = BasicCheckpoint(every=5, filepath='data/hillClimberBackup.csv')

    # add callback into callback list
    algo.addCallback(callback)

    # run the algorithm using local search and get solution found 
    solution = algo.run(evaluations=100)
    print(solution.fitness)


.. note::
    It is possible to add as many callbacks as desired in the algorithm in question.


Previously, some methods of the abstract ``Algorithm`` class have not been presented. These methods are linked to the use of callbacks, 
in particular the ``addCallback`` method which allows the addition of a callback to an algorithm instance as seen above.

- The ``resume`` method will reload all callbacks list using ``load`` method.
- The ``progress`` method will ``run`` each callbacks during the algorithm search.

If we want to exploit this functionality, then we will need to exploit them within our algorithm. Let's make the necessary modifications for our algorithm ``IteratedLocalSearch``:


.. code-block:: python

    """
    module imports
    """
    from macop.algorithms.base import Algorithm

    class IteratedLocalSearch(Algorithm):
        
        ...

        def run(self, evaluations, ls_evaluations=100):
            """
            Run the iterated local search algorithm using local search
            """

            # by default use of mother method to initialise variables
            super().run(evaluations)

            # initialise current solution
            self.initRun()

            # restart using callbacks backup list
            self.resume()

            # local search algorithm implementation
            while not self.stop():

                # create and search solution from local search
                newSolution = self._localSearch.run(ls_evaluations)

                # if better solution than currently, replace it
                if self.isBetter(newSolution):
                    self._bestSolution = newSolution

                # check if necessary to call each callbacks
                self.progress()

                self.information()

            return self._bestSolution


All the features of **Macop** were presented. The next section will aim to quickly present the few implementations proposed within **Macop** to highlight the modulality of the package.Validate a solution
======================

When an optimisation problem requires respecting certain constraints, Macop allows you to quickly verify that a solution is valid. 
It is based on a defined function taking a solution as input and returning the validity criterion (true or false).

Validator definition
~~~~~~~~~~~~~~~~~~~~~~~~~

An invalid solution can be shown below where the sum of the object weights is greater than 15:

.. image:: ../_static/documentation/project_knapsack_invalid.png
   :width:  85 %
   :align: center

In fact, **[1, 0, 1, 0, 0]** is an invalid solution as we have a weight of **16** which violates the knapsack capacity constraint.

To avoid taking into account invalid solutions, we can define our function which will validate or not a solution based on our problem instance:

.. code-block:: python

    """
    Problem instance definition
    """

    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    """
    Validator function definition
    """
    def validator(solution):

        weight_sum = 0

        for i, w in enumerate(elements_weight):
            # add weight if current object is set to 1
            weight_sum += w * solution.getdata = )[i]
        
        # validation condition
        return weight_sum <= 15


Use of validator
~~~~~~~~~~~~~~~~~~~~~

We can now generate solutions randomly by passing our validation function as a parameter:

.. code-block:: python

    """
    Problem instance definition
    """
    ...
    
    """
    Validator function definition
    """
    ...

    # ensure valid solution
    solution = BinarySolution.random(5, validator)


.. caution::
    If the search space for valid solutions is very small compared to the overall search space, this can involve a considerable time for validating the solution and therefore obtaining a solution.

The validation of a solution is therefore now possible. In the next part we will focus on the evaluation of a solution.Operator choices
===================

The ``policy`` feature of **Macop** enables to choose the next operator to apply during the search process of the algorithm based on specific criterion.

Why using policy ?
~~~~~~~~~~~~~~~~~~~~~~~

Sometimes the nature of the problem and its instance can strongly influence the search results when using mutation operators or crossovers. 
Automated operator choice strategies have also been developed in the literature, notably based on reinforcement learning.

The operator choice problem can be seen as the desire to find the best solution generation operator at the next evaluation that will be the most conducive to precisely improving the solution.

.. image:: ../_static/documentation/operators_choice.png
   :width:  45 %
   :align: center

.. note::
    An implementation using reinforcement learning has been developed as an example in the ``macop.policies.reinforcement`` module. 
    However, it will not be detailed here. You can refer to the API documentation for more details.


Custom policy
~~~~~~~~~~~~~~~~~~

In our case, we are not going to exploit a complex enough implementation of a ``policy``. Simply, we will use a random choice of operator.

First, let's take a look of the ``policy`` abstract class available in ``macop.policies.base``:

.. code-block:: python

    class Policy():

        def __init__(self, operators):
            self.operators = operators

        @abstractmethod
        def select(self):
            """
            Select specific operator
            """
            pass

        def apply(self, solution):
            """
            Apply specific operator to create new solution, compute its fitness and return it
            """
            ...

        def setAlgo(self, algo):
            """
            Keep into policy reference of the whole algorithm
            """
            ...


``Policy`` instance will have of ``operators`` attributs in order to keep track of possible operators when selecting one. 
Here, in our implementation we only need to specify the ``select`` abstract method. The ``apply`` method will select the next operator and return the new solution.

.. code-block:: python

    """
    module imports
    """
    from macop.policies.base import Policy

    class RandomPolicy(Policy):

        def select(self):
            """
            Select specific operator
            """
            # choose operator randomly
            index = random.randint(0, len(self.operators) - 1)
            return self.operators[index]


We can now use this operator choice policy to update our current solution:


.. code-block:: python

    """
    Operators instances
    """
    mutator = SimpleMutation()
    crossover = SimpleCrossover()

    """
    RandomPolicy instance
    """
    policy = RandomPolicy([mutator, crossover])

    """
    Current solutions instance
    """
    solution1 = BinarySolution.random(5)
    solution2 = BinarySolution.random(5)

    # pass two solutions in parameters in case of selected crossover operator
    new_solution = policy.apply(solution1, solution2)

.. caution::
    By default if ``solution2`` parameter is not provided into ``policy.apply`` method for crossover, the best solution known is used from the algorithm linked to the ``policy``.

Updating solutions is therefore now possible with our policy. It is high time to dive into the process of optimizing solutions and digging into our research space.Use of evaluators
====================

Now that it is possible to generate a solution randomly or not. It is important to know the value associated with this solution. We will then speak of evaluation of the solution. With the score associated with it, the `fitness`.

Generic evaluator
~~~~~~~~~~~~~~~~~~~~~~

As for the management of solutions, a generic evaluator class ``macop.evaluators.base.Evaluator`` is developed within **Macop**:

Abstract Evaluator class is used for computing fitness score associated to a solution. To evaluate all the solutions, this class:

- stores into its ``_data`` dictionary attritute required measures when computing a solution
- has a ``compute`` abstract method enable to compute and associate a score to a given solution
- stores into its ``algo`` attritute the current algorithm to use (we will talk about algorithm later)

.. code-block: python

    class Evaluator():
    """
    Abstract Evaluator class which enables to compute solution using specific `_data` 
    """
    def __init__(self, data):
        self._data = data

    @abstractmethod
    def compute(self, solution):
        """
        Apply the computation of fitness from solution
        """
        pass

    def setAlgo(self, algo):
        """
        Keep into evaluator reference of the whole algorithm
        """
        self.algo = algo

We must therefore now create our own evaluator based on the proposed structure.

Custom evaluator
~~~~~~~~~~~~~~~~~~~~~

To create our own evaluator, we need both:

- data useful for evaluating a solution
- calculate the score (fitness) associated with the state of the solution from these data. Hence, implement specific ``compute`` method.

We will define the ``KnapsackEvaluator`` class, which will therefore allow us to evaluate solutions to our current problem.

.. code-block:: python

    """
    modules imports
    """
    from macop.evaluators.base import Evaluator

    class KnapsackEvaluator(Evaluator):
        
        def compute(solution):

            # `_data` contains worths array values of objects
            fitness = 0
            for index, elem in enumerate(solution.getdata = )):
                fitness += self._data['worths'][index] * elem

            return fitness


It is now possible to initialise our new evaluator with specific data of our problem instance:

.. code-block:: python

    """
    Problem instance definition
    """
    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    """
    Evaluator problem instance
    """
    evaluator = KnapsackEvaluator(data={'worths': elements_score})

    # using defined BinarySolution
    solution = BinarySolution.random(5)

    # obtaining current solution score
    solution_fitness = solution.evaluate(evaluator)

    # score is also stored into solution
    solution_fitness = solution.fitness

.. note::
    The current developed ``KnapsackEvaluator`` is available into ``macop.evaluators.mono.KnapsackEvaluator`` in **Macop**.

In the next part we will see how to modify our current solution with the use of modification operator.A tour of Macop
===================

.. image:: ../_static/logo_macop.png
   :width: 300 px
   :align: center

This documentation will allow a user who wishes to use the **Macop** optimisation package to understand both how it works and offers examples of how to implement specific needs.

It will gradually take up the major ideas developed within **Macop** to allow for quick development. You can navigate directly via the menu available below to access a specific part of the documentation.

.. toctree::
   :maxdepth: 1
   :numbered:
   :caption: Contents:

   introduction
   problem
   solutions
   validator
   evaluators
   operators
   policies
   algorithms
   callbacks
   othersOptimisation process
=======================

Let us now tackle the interesting part concerning the search for optimum solutions in our research space.

Find local and global optima
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Overall, in an optimization process, we will seek to find the best, or the best solutions that minimize or maximize our objective function (fitness score obtained) in order to respond to our problem.

.. image:: ../_static/documentation/search_space.png
   :width:  95 %
   :align: center

Sometimes, the search space can be very simple. A local search can provide access to the global optimum as shown in figure (a) above. 
In other cases, the search space is more complex. It may be necessary to explore more rather than exploit in order to get out of a convex zone and not find the global optimum but only a local opmatime solution. 
This problem is illustrated in figure (b).

Abstract algorithm class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An abstract class is proposed within Macop to generalize the management of an algorithm and therefore of a heuristic. 
It is located in the ``macop.algorithms.base`` module. 

We will pay attention to the different methods of which she is composed. This class enables to manage some common usages of operation research algorithms:

- initialization function of solution
- validator function to check if solution is valid or not (based on some criteria)
- evaluation function to give fitness score to a solution
- operators used in order to update solution during search process
- policy process applied when choosing next operator to apply
- callbacks function in order to do some relative stuff every number of evaluation or reload algorithm state
- parent algorithm associated to this new algorithm instance (hierarchy management)

She is composed of few default attributes:

- initialiser: {function} -- basic function strategy to initialise solution
- evaluator: {Evaluator} -- evaluator instance in order to obtained fitness (mono or multiple objectives)
- operators: {[Operator]} -- list of operator to use when launching algorithm
- policy: {Policy} -- Policy instance strategy to select operators
- validator: {function} -- basic function to check if solution is valid or not under some constraints
- maximise: {bool} -- specify kind of optimisation problem 
- verbose: {bool} -- verbose or not information about the algorithm
- currentSolution: {Solution} -- current solution managed for current evaluation comparison
- bestSolution: {Solution} -- best solution found so far during running algorithm
- callbacks: {[Callback]} -- list of Callback class implementation to do some instructions every number of evaluations and `load` when initialising algorithm
- parent: {Algorithm} -- parent algorithm reference in case of inner Algorithm instance (optional)

.. code-block:: python

    class Algorithm():

        def __init__(self,
                    initialiser,
                    evaluator,
                    operators,
                    policy,
                    validator,
                    maximise=True,
                    parent=None,
                    verbose=True):
            ...

        def addCallback(self, callback):
            """
            Add new callback to algorithm specifying usefull parameters
            """
            ...

        def resume(self):
            """
            Resume algorithm using Callback instances
            """
            ...

        def getParent(self):
            """
            Recursively find the main parent algorithm attached of the current algorithm
            """
            ...

        def setParent(self, parent):
            """
            Set parent algorithm to current algorithm
            """
            ...


        def initRun(self):
            """
            initialise the current solution and best solution using the `initialiser` function
            """
            ...

        def increaseEvaluation(self):
            """
            Increase number of evaluation once a solution is evaluated for each dependant algorithm (parents hierarchy)
            """
            ...
                
        def getGlobalEvaluation(self):
            """
            Get the global number of evaluation (if inner algorithm)
            """
            ...

        def getGlobalMaxEvaluation(self):
            """
            Get the global max number of evaluation (if inner algorithm)
            """
            ...

        def stop(self):
            """
            Global stopping criteria (check for parents algorithm hierarchy too)
            """
            ...

        def evaluate(self, solution):
            """
            Evaluate a solution using evaluator passed when intialize algorithm
            """
            ...

        def update(self, solution):
            """
            Apply update function to solution using specific `policy`
            Check if solution is valid after modification and returns it
            """
            ...

        def isBetter(self, solution):
            """
            Check if solution is better than best found
            """
            ...

        def run(self, evaluations):
            """
            Run the specific algorithm following number of evaluations to find optima
            """
            ...

        def progress(self):
            """
            Log progress and apply callbacks if necessary
            """
            ...


The notion of hierarchy between algorithms is introduced here. We can indeed have certain dependencies between algorithms. 
The methods ``increaseEvaluation``, ``getGlobalEvaluation`` and ``getGlobalMaxEvaluation`` ensure that the expected global number of evaluations is correctly managed, just like the ``stop`` method for the search stop criterion.

The ``evaluate``, ``update`` and ``isBetter`` will be used a lot when looking for a solution in the search space. 
In particular the ``update`` function, which will call the ``policy`` instance to generate a new valid solution.
``isBetter`` method is also overloadable especially if the algorithm does not take any more into account than a single solution to be verified (verification via a population for example).

The ``initRun`` method specify the way you intialise your algorithm (``bestSolution`` and ``currentSolution`` as example) if algorithm not already initialised.

.. note:: 
    The ``initRun`` method can also be used for intialise population of solutions instead of only one best solution, if you want to manage a genetic algorithm.

Most important part is the ``run`` method. Into abstract, the ``run`` method only initialised the current number of evaluation for the algorithm based on the parent algorithm if we are into inner algorithm.
It is always **mandatory** to call the parent class ``run`` method using ``super().run(evaluations)``. Then, using ``evaluations`` parameter which is the number of evaluations budget to run, we can process or continue to find solutions into search space.

.. warning::
    The other methods such as ``addCallback``, ``resume`` and ``progress`` will be detailed in the next part focusing on the notion of callback.

Local search algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are going to carry out our first local search algorithm within our search space. A `local search` consists of starting from a solution, then applying a mutation or crossover operation to it, in order to obtain a new one. 
This new solution is evaluated and retained if it is better. We will speak here of the notion of **neighborhood exploration**. The process is then completed in the same way. 
The local search ends after a certain number of evaluations and the best evaluated solution obtained is returned.

Let's implement an algorithm well known under the name of hill climber best improvment inheriting from the mother algorithm class and name it ``HillClimberBestImprovment``.


.. code-block:: python

    """
    module imports
    """
    from macop.algorithms.base import Algorithm

    class HillClimberBestImprovment(Algorithm):

        def run(self, evaluations):
            """
            Run a local search algorithm
            """

            # by default use of mother method to initialise variables
            super().run(evaluations)

            # initialise current solution and best solution
            self.initRun()

            solutionSize = self._currentSolution.size

            # local search algorithm implementation
            while not self.stop():

                for _ in range(solutionSize):

                    # update current solution using policy
                    newSolution = self.update(self._currentSolution)

                    # if better solution than currently, replace it
                    if self.isBetter(newSolution):
                        self._bestSolution = newSolution

                    # increase number of evaluations
                    self.increaseEvaluation()

                    # stop algorithm if necessary
                    if self.stop():
                        break

                # set new current solution using best solution found in this neighbor search
                self._currentSolution = self._bestSolution
            
            return self._bestSolution

Our algorithm is now ready to work. As previously, let us define two operators as well as a random choice strategy. 
We will also need to define a **solution initialisation function** so that the algorithm can generate new solutions.


.. code-block:: python

    """
    Problem instance definition
    """
    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    # evaluator instance
    evaluator = KnapsackEvaluator(data={'worths': elements_score})

    # valid instance using lambda
    validator = lambda solution: sum([ elements_weight[i] * solution.getdata = )[i] for i in range(len(solution.getdata = )))]) <= 15
    
    # initialiser instance using lambda with default param value
    initialiser = lambda x=5: BinarySolution.random(x, validator)
    
    # operators list with crossover and mutation
    operators = [SimpleCrossover(), SimpleMutation()]
    
    # policy random instance
    policy = RandomPolicy(operators)
    
    # maximizing algorithm (relative to knapsack problem)
    algo = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)

    # run the algorithm and get solution found
    solution = algo.run(100)
    print(solution.fitness)


.. note::
    The ``verbose`` algorithm parameter will log into console the advancement process of the algorithm is set to ``True`` (the default value).

Exploratory algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~

As explained in **figure (b)** of **section 8.1**, sometimes the search space is more complicated due to convex parts and need heuristic with other strategy rather than a simple local search.

The way to counter this problem is to allow the algorithm to exit the exploitation phase offered by local search. But rather to seek to explore other parts of the research space. This is possible by simply carrying out several local searches with our budget (number of evaluations).

The idea is to make a leap in the search space in order to find a new local optimum which can be the global optimum. The explained process is illustrated below:

.. image:: ../_static/documentation/search_space_simple.png
   :width:  45 %
   :align: center


We are going to implement a more specific algorithm allowing to take a new parameter as input. This is a local search, the one previously developed. For that, we will have to modify the constructor a little.
Let's called this new algorithm ``IteratedLocalSearch``:

.. code-block:: python

    """
    module imports
    """
    from macop.algorithms.base import Algorithm

    class IteratedLocalSearch(Algorithm):
        
        def __init__(self,
                    initialiser,
                    evaluator,
                    operators,
                    policy,
                    validator,
                    localSearch,
                    maximise=True,
                    parent=None,
                    verbose=True):
            
            super().__init__(initialiser, evaluator, operators, policy, validator, maximise, parent, verbose)

            # specific local search associated with current algorithm
            self._localSearch = localSearch

            # need to attach current algorithm as parent
            self._localSearch.setParent(self)


        def run(self, evaluations, ls_evaluations=100):
            """
            Run the iterated local search algorithm using local search
            """

            # by default use of mother method to initialise variables
            super().run(evaluations)

            # initialise current solution
            self.initRun()

            # local search algorithm implementation
            while not self.stop():

                # create and search solution from local search (stop method can be called inside local search)
                newSolution = self._localSearch.run(ls_evaluations)

                # if better solution than currently, replace it
                if self.isBetter(newSolution):
                    self._bestSolution = newSolution

                self.information()

            return self._bestSolution

In the initialization phase we have attached our local search passed as a parameter with the current algorithm as parent. 
The goal is to touch keep track of the overall search evaluation number (relative to the parent algorithm).

Then, we use this local search in our ``run`` method to allow a better search for solutions.

.. code-block:: python

    """
    Problem instance definition
    """
    elements_score = [ 4, 2, 10, 1, 2 ] # worth of each object
    elements_weight = [ 12, 1, 4, 1, 2 ] # weight of each object

    # evaluator instance
    evaluator = KnapsackEvaluator(data={'worths': elements_score})

    # valid instance using lambda
    validator = lambda solution: sum([ elements_weight[i] * solution.getdata = )[i] for i in range(len(solution.getdata = )))]) <= 15
    
    # initialiser instance using lambda with default param value
    initialiser = lambda x=5: BinarySolution.random(x, validator)
    
    # operators list with crossover and mutation
    operators = [SimpleCrossover(), SimpleMutation()]
    
    # policy random instance
    policy = RandomPolicy(operators)
    
    # maximizing algorithm (relative to knapsack problem)
    localSearch = HillClimberBestImprovment(initialiser, evaluator, operators, policy, validator, maximise=True, verbose=False)
    algo = IteratedLocalSearch(initialiser, evaluator, operators, policy, validator, localSearch=local_search, maximise=True, verbose=False)

    # run the algorithm using local search and get solution found 
    solution = algo.run(evaluations=100, ls_evaluations=10)
    print(solution.fitness)


.. note:: 
    These two last algorithms developed are available in the library within the module ``maocp.algorithms.mono``.

We have one final feature to explore in the next part. This is the notion of ``callback``.Introduction
================

`Macop` is a python package for solving discrete optimisation problems in nature. Continuous optimisation is also applicable but not yet developed. The objective is to allow a user to exploit the basic structure proposed by this package to solve a problem specific to him. The interest is that he can quickly abstract himself from the complications related to the way of evaluating, comparing, saving the progress of the search for good solutions but rather concentrate if necessary on his own algorithm. Indeed, `Macop` offers the following main and basic features: 

- **solutions:** representation of the solution;
- **validator:** such as constraint programming, a `validator` is a function which is used to validate or not a solution data state;
- **evaluator:**  stores problem instance data and implements a `compute` method in order to evaluate a solution;
- **operators:** mutators, crossovers operators to update and obtain new solution;
- **policies:** the way you choose the available operators (might be using reinforcement learning);
- **algorithms:** generic and implemented optimisation research algorithms;
- **callbacks:** callbacks to automatically keep track of the search space advancement and restart from previous state if nedded.

.. image:: ../_static/documentation/macop_behaviour.png
   :width: 50 %
   :align: center

Based on all of these generic and/or implemented functionalities, the user will be able to quickly develop a solution to his problem while retaining the possibility of remaining in control of his development by overloading existing functionalities if necessary.Solutions
=============

Representing a solution to a specific problem is very important in an optimisation process. In this example, we will always use the **knapsack problem** as a basis.

In a first step, the management of the solutions by the macop package will be presented. Then a specific implementation for the current problem will be detailed.

Generic Solution
~~~~~~~~~~~~~~~~~~~~~~~~~

Inside ``macop.solutions.base`` module of `Macop`, the ``Solution`` class is available. It's an abstract solution class structure which:

- stores the solution data representation into its ``data`` attribute
- get ``size`` (shape) of specific data representation
- stores the ``score`` of the solution once a solution is evaluated

Some specific methods are available:

.. code-block:: python

    class Solution():

        def __init__(self, data, size):
            """
            Abstract solution class constructor
            """
            ...

        def isValid(self, validator):
            """
            Use of custom function which checks if a solution is valid or not
            """
            ...

        def evaluate(self, evaluator):
            """
            Evaluate solution using specific `evaluator`
            """
            ...

        def fitness(self):
            """
            Returns fitness score
            """
            ...

        @staticmethod
        def random(size, validator=None):
            """
            initialise solution using random data with validator or not
            """
            ...

        def clone(self):
            """
            Clone the current solution and its data, but without keeping evaluated `_score`
            """
            ...


From these basic methods, it is possible to manage a representation of a solution to our problem. 

Allowing to initialise it randomly or not (using constructor or ``random`` method), to evaluate it (``evaluate`` method) and to check some constraints of validation of the solution (``isValid`` method).

.. note::
    Only one of these methods needs specification if we create our own type of solution. This is the ``random`` method, which depends on the need of the problem.

We will now see how to define a type of solution specific to our problem.

Solution representation for knapsack
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will now use the abstract ``Solution`` type available in the ``macop.solutions.base`` module in order to define our own solution.
First of all, let's look at the representation of our knapsack problem. **How to represent the solution?**

Knapsack solution
************************

A valid solution can be shown below where the sum of the object weights is 15 and the sum of the selected objects values is 8 (its fitness):

.. image:: ../_static/documentation/project_knapsack_solution.png
   :width:  85 %
   :align: center

Its representation can be translate as a **binary array** with value:

.. code-block::

    [1, 1, 0, 0, 1]

where selected objects have **1** as value otherwise **0**.

Binary Solution
**********************

We will now define our own type of solution by inheriting from ``macop.solutions.base.Solution``, which we will call ``BinarySolution``.

First we will define our new class as inheriting functionality from ``Solution`` (such as child class). 
We will also have to implement the ``random`` method to create a new random solution.

.. code-block:: python

    """
    modules imports
    """
    from macop.solutions.base import Solution
    import numpy as np

    class BinarySolution(Solution):
        
        @staticmethod
        def random(size, validator=None):

            # create binary array of specific size using numpy random module
            data = np.random.randint(2, size=size)
            # initialise new solution using constructor
            solution = BinarySolution(data, size)

            # check if validator is set
            if not validator:
                return solution

            # try to generate solution until solution validity (if validator is provided)
            while not validator(solution):
                data = np.random.randint(2, size=size)
                solution = BinarySolution(data, size)

            return solution

.. note::
    The current developed ``BinarySolution`` is available into ``macop.solutions.discrete.BinarySolution`` in **Macop**.

Using this new Solution representation, we can now generate solution randomly:

.. code-block:: python

    solution = BinarySolution.random(5)

In the next part, we will see how to verify that a solution meets certain modeling constraints of the problem.Apply operators to solution
==============================

Applying an operator to a solution consists of modifying the current state of the solution in order to obtain a new one. The goal is to find a better solution in the search space.

Operators definition
~~~~~~~~~~~~~~~~~~~~~~~~~

In the discrete optimisation literature, we can categorise operators into two sections:

- **mutators**: modification of one or more elements of a solution from its current state.
- **crossovers**: Inspired by Darwin's theory of evolution, we are going here from two solutions to generate a so-called offspring solution composed of the fusion of the data of the parent solutions.

Inside **Macop**, operators are also decomposed into these two categories. Inside ``macop.operators.discrete.base``, generic class ``Operator`` enables to manage any kind of operator.

.. code-block:: python

    class Operator():
        """
        Abstract Operator class which enables to update solution applying operator (computation)
        """
        @abstractmethod
        def __init__(self):
            pass

        @abstractmethod
        def apply(self, solution):
            """
            Apply the current operator transformation
            """
            pass

        def setAlgo(self, algo):
            """
            Keep into operator reference of the whole algorithm
            """
            self.algo = algo

Like the evaluator, the operator keeps **track of the algorithm** (using ``setAlgo`` method) to which he will be linked. This will allow better management of the way in which the operator must take into account the state of current data relating to the evolution of research.

``Mutation`` and ``Crossover`` classes inherite from ``Operator``. An ``apply`` function is required for any new operator.

.. code-block:: python

    class Mutation(Operator):
        """Abstract Mutation extend from Operator

        Attributes:
            kind: {KindOperator} -- specify the kind of operator
        """
        def __init__(self):
            self._kind = KindOperator.MUTATOR

        def apply(self, solution):
            raise NotImplementedError


    class Crossover(Operator):
        """Abstract crossover extend from Operator

        Attributes:
            kind: {KindOperator} -- specify the kind of operator
        """
        def __init__(self):
            self._kind = KindOperator.CROSSOVER

        def apply(self, solution1, solution2):
            raise NotImplementedError

We will now detail these categories of operators and suggest some relative to our problem.

Mutator operator
~~~~~~~~~~~~~~~~~~~~~

As detailed, the mutation operator consists in having a minimum impact on the current state of our solution. Here is an example of a modification that could be done for our problem.

.. image:: ../_static/documentation/project_knapsack_mutator.png
   :width:  90 %
   :align: center

In this example we change a bit value randomly and obtain a new solution from our search space.

.. warning::
    Applying an operator can conduct to a new but invalid solution from the search space.

The modification applied here is just a bit swapped. Let's define the ``SimpleBinaryMutation`` operator, allows to randomly change a binary value of our current solution.


.. code-block:: python

    """
    modules imports
    """
    from macop.operators.discrete.base import Mutation

    class SimpleBinaryMutation(Mutation):

        def apply(self, solution):
            
            # obtain targeted cell using solution size
            size = solution.size
            cell = random.randint(0, size - 1)

            # copy of solution
            copy_solution = solution.clone()

            # swicth values
            if copy_solution.getdata = )[cell]:
                copy_solution.getdata = )[cell] = 0
            else:
                copy_solution.getdata = )[cell] = 1

            # return the new obtained solution
            return copy_solution

We can now instanciate our new operator in order to obtain a new solution:


.. code-block:: python

    """
    BinaryMutator instance
    """
    mutator = SimpleBinaryMutation()

    # using defined BinarySolution
    solution = BinarySolution.random(5)

    # obtaining new solution using operator
    new_solution = mutator.apply(solution)


.. note::
    The developed ``SimpleBinaryMutation`` is available into ``macop.operators.discrete.mutators.SimpleBinaryMutation`` in **Macop**.


Crossover operator
~~~~~~~~~~~~~~~~~~~~~~~


Inspired by Darwin's theory of evolution, crossover starts from two solutions to generate a so-called offspring solution composed of the fusion of the data of the parent solutions.

.. image:: ../_static/documentation/project_knapsack_crossover.png
   :width:  95%
   :align: center

In this example we merge two solutions with a specific splitting criterion in order to obtain an offspring.

We will now implement the SimpleCrossover crossover operator, which will merge data from two solutions. 
The first half of solution 1 will be saved and added to the second half of solution 2 to generate the new solution (offspring).


.. code-block:: python

    """
    modules imports
    """
    from macop.operators.discrete.base import Crossover

    class SimpleCrossover(Crossover):

        def apply(self, solution1, solution2):
            
            size = solution1.size

            # default split index used
            splitIndex = int(size / 2)

            # copy data of solution 1
            firstData = solution1._data.copy()

            # copy of solution 2
            copy_solution = solution2.clone()

            copy_solution.getdata = )[splitIndex:] = firstData[splitIndex:]

            return copy_solution


We can now use the crossover operator created to generate new solutions. Here is an example of use:

.. code-block:: python

    """
    SimpleCrossover instance
    """
    crossover = SimpleCrossover()

    # using defined BinarySolution
    solution1 = BinarySolution.random(5)
    solution2 = BinarySolution.random(5)

    # obtaining new solution using crossover
    offspring = crossover.apply(solution1, solution2)

.. tip::
    The developed ``SimpleCrossover`` is available into ``macop.operators.discrete.crossovers.SimpleCrossover`` in **Macop**.
    However, the choice of halves of the merged data is made randomly.

Next part introduce the ``policy`` feature of **Macop** which enables to choose the next operator to apply during the search process based on specific criterion.Implementation examples
=======================

Within the API of **Macop**, you can find an implementation of The Multi-objective evolutionary algorithm based on decomposition (MOEA/D) is a general-purpose algorithm for approximating the Pareto set of multi-objective optimization problems. 
It decomposes the original multi-objective problem into a number of single-objective optimization sub-problems and then uses an evolutionary process to optimize these sub-problems simultaneously and cooperatively. 
MOEA/D is a state-of-art algorithm in aggregation-based approaches for multi-objective optimization.

.. image:: ../_static/documentation/search_space_moead.png
   :width:  45 %
   :align: center


As illustrated below, the two main objectives are sub-divised into 5 single-objective optimization sub-problems in order to find the Pareto front.

- ``macop.algorithms.multi.MOSubProblem`` class defines each sub-problem of MOEA/D.
- ``macop.algorithms.multi.MOEAD`` class exploits ``MOSubProblem`` and implements MOEA/D using weighted-sum of objectives method.

An example with MOEAD for knapsack problem is available in knapsackMultiExample.py_.

.. _knapsackMultiExample.py: https://github.com/jbuisine/macop/blob/master/examples/knapsackMultiExample.pyUBQP problem definition
=======================

Understand the UBQP Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a collection of :math:`n` items such that each pair of items is associated with a profit value that can be positive, negative or zero, unconstrained binary quadratic programming (UBQP) seeks a subset of items that maximizes the sum of their paired values. The value of a pair is accumulated in the sum only if the two corresponding items are selected. A feasible solution to a UBQP instance can be specified by a binary string of size :math:`n`, such that each variable indicates whether the corresponding item is included in the selection or not.


Mathematical definition
~~~~~~~~~~~~~~~~~~~~~~~

More formally, the conventional and single-objective UBQP problem is to maximize the following objective function:

:math:`f(x)=x′Qx=\sum_{i=1}^{n}{\sum_{j=1}^{n}{q_{ij}⋅x_i⋅x_j}}`

where :math:`Q=(q_{ij})` is an :math:`n` by :math:`n` matrix of constant values, :math:`x` is a vector of :math:`n` binary (zero-one) variables, i.e., :math:`x \in \{0, 1\}`, :math:`i \in \{1,...,n\}`, and :math:`x'` is the transpose of :math:`x`.
Macop UBQP implementation
=========================

Let's see how it is possible with the use of the **Macop** package to implement and deal with this UBQP instance problem.

Solution structure definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Firstly, we are going to use a type of solution that will allow us to define the structure of our solutions.

The available ``macop.solutions.discrete.BinarySolution`` type of solution within the Macop package represents exactly what one would wish for. 

Let's see an example of its use:

.. code:: python

    from macop.solutions.discrete import BinarySolution
    
    solution = BinarySolution.random(10)
    print(solution)


The resulting solution obtained:

.. code:: bash

    Binary solution [1 0 1 1 1 0 0 1 1 0]


UBQP Evaluator
~~~~~~~~~~~~~~

Now that we have the structure of our solutions, and the means to generate them, we will seek to evaluate them.

To do this, we need to create a new evaluator specific to our problem and the relative evaluation function we need to maximise:

- :math:`f(x)=x′Qx=\sum_{i=1}^{n}{\sum_{j=1}^{n}{q_{ij}⋅x_i⋅x_j}}`

So we are going to create a class that will inherit from the abstract class ``macop.evalutarors.base.Evaluator``:

.. code:: python

    from macop.evaluators.base import Evaluator

    class UBQPEvaluator(Evaluator):
    """UBQP evaluator class which enables to compute UBQP solution using specific `_data`

    - stores into its `_data` dictionary attritute required measures when computing a UBQP solution
    - `_data['Q']` matrix of size n x n with real values data (stored as numpy array)
    - `compute` method enables to compute and associate a score to a given UBQP solution
    """

    def compute(self, solution):
        """Apply the computation of fitness from solution

        Args:
            solution: {Solution} -- UBQP solution instance
    
        Returns:
            {float} -- fitness score of solution
        """
        fitness = 0
        for index_i, val_i in enumerate(solution.getdata = )):
            for index_j, val_j in enumerate(solution.getdata = )):
                fitness += self._data['Q'][index_i, index_j] * val_i * val_j

        return fitness

The cost function for the Unconstrained binary quadratic problem is now well defined.

.. tip::
    The class proposed here, is available in the Macop package ``macop.evaluators.discrete.mono.UBQPEvaluator``.

Running algorithm
~~~~~~~~~~~~~~~~~

Now that the necessary tools are available, we will be able to deal with our problem and look for solutions in the search space of our UBQP instance.

Here we will use local search algorithms already implemented in **Macop**.

If you are uncomfortable with some of the elements in the code that will follow, you can refer to the more complete **Macop** documentation_ that focuses more on the concepts and tools of the package.

.. code:: python

    # main imports
    import numpy as np

    # module imports
    from macop.solutions.discrete import BinarySolution
    from macop.evaluators.discrete.mono import UBQPEvaluator

    from macop.operators.discrete.mutators import SimpleMutation, SimpleBinaryMutation

    from macop.policies.classicals import RandomPolicy

    from macop.algorithms.mono import IteratedLocalSearch as ILS
    from macop.algorithms.mono import HillClimberFirstImprovment

    # usefull instance data
    n = 100
    qap_instance_file = 'qap_instance.txt'

    # default validator
    def validator(solution):
        return True

    # define init random solution
    def init():
        return BinarySolution.random(n, validator)

    # load UBQP instance
    with open(ubqp_instance_file, 'r') as f:

        lines = f.readlines()

        # get all string floating point values of matrix
        Q_data = ''.join([ line.replace('\n', '') for line in lines[8:] ])

        # load the concatenate obtained string
        Q_matrix = np.fromstring(Q_data, dtype=float, sep=' ').reshape(n, n)

    print(f'Q_matrix shape: {Q_matrix.shape}')

    # only one operator here
    operators = [SimpleMutation(), SimpleBinaryMutation()]

    # random policy
    policy = RandomPolicy(operators)

    # use of loaded data from UBQP instance
    evaluator = UBQPEvaluator(data={'Q': Q_matrix})

    # passing global evaluation param from ILS
    hcfi = HillClimberFirstImprovment(init, evaluator, operators, policy, validator, maximise=True, verbose=True)
    algo = ILS(init, evaluator, operators, policy, validator, localSearch=hcfi, maximise=True, verbose=True)

    # run the algorithm
    bestSol = algo.run(10000, ls_evaluations=100)

    print('Solution for UBQP instance score is {}'.format(evaluator.compute(bestSol)))


UBQP problem solving is now possible with **Macop**. As a reminder, the complete code is available in the ubqpExample.py_ file.

.. _ubqpExample.py: https://github.com/jbuisine/macop/blob/master/examples/ubqpExample.py
.. _documentation: https://jbuisine.github.io/macop/_build/html/documentationsUBQP Problem instance generation
================================

To define our quadratic assignment problem instance, we will use the available mUBQP_ multi-objective quadratic problem generator. 

Genration of the instance
~~~~~~~~~~~~~~~~~~~~~~~~~

We will limit ourselves here to a single objective for the purposes of this example. The available file **mubqpGenerator.R**, will be used to generate the instance (using R language).

.. code:: bash

    Rscript mubqpGenerator.R 0.8 1 100 5 42 ubqp_instance.txt

The main parameters used for generating our UBQP instance are:

- **ρ:** the objective correlation coefficient
- **M:** the number of objective functions
- **N:** the length of bit strings
- **d:** the matrix density (frequency of non-zero numbers)
- **s:** seed to use

.. _mUBQP: http://mocobench.sourceforge.net/index.php?n=Problem.MUBQP

.. _ubqp_instance.txt: https://github.com/jbuisine/macop/blob/master/examples/instances/ubqp/ubqp_instance.txt

Load data instance
~~~~~~~~~~~~~~~~~~

We are now going to load this instance via a Python code which will be useful to us later on:

.. code:: Python

    qap_instance_file = 'ubqp_instance.txt'

    n = 100 # the instance size

    # load UBQP instance
    with open(ubqp_instance_file, 'r') as f:

        lines = f.readlines()

        # get all string floating point values of matrix
        Q_data = ''.join([ line.replace('\n', '') for line in lines[8:] ])

        # load the concatenate obtained string
        Q_matrix = np.fromstring(Q_data, dtype=float, sep=' ').reshape(n, n)

    print(f'Q_matrix {Q_matrix.shape}')

.. note::
    As we know the size of our instance and the structure of the document (header size), it is quite quick to look for the lines related to the :math:`Q` matrix.Unconstrained Binary Quadratic Programming
==========================================

Given a collection of :math:`n` items such that each pair of items is associated with a profit value that can be positive, negative or zero, unconstrained binary quadratic programming (UBQP) seeks a subset of items that maximizes the sum of their paired values. The value of a pair is accumulated in the sum only if the two corresponding items are selected.

The UBQP problem will be tackle in this example.

.. toctree::
   :maxdepth: 1
   :numbered:
   :caption: Contents:

   problem
   instance
   implementationQAP problem definition
======================

The quadratic assignment problem (QAP) was introduced by Koopmans and Beckman in 1957 in the context of locating "indivisible economic activities". The objective of the problem is to assign a set of facilities to a set of locations in such a way as to minimize the total assignment cost. The assignment cost for a pair of facilities is a function of the flow between the facilities and the distance between the locations of the facilities.

Location assignment example
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider a **facility location problem** with **four** facilities (and **four** locations). One possible assignment is shown in the figure below: facility 4 is assigned to location 1, facility 1 
is assigned to location 2, facility 3 is assigned to location 3, and facility 2 is assigned to location 3. This assignment can be written as the permutation :math:`p=\{4,1,3,2\}`, 
which means that facility 4 is assigned to location 1, facility 1 is assigned to location 2, facility 3 is assigned to location 3, and facility 2 is assigned to location 3. 
In the figure, the line between a pair of facilities indicates that there is required flow between the facilities, and the thickness of the line increases with the value of the flow. 

.. image:: ../../_static/examples/qap/factories_qap.png
   :width: 50 %
   :align: center
   :alt: Example of QAP facilities to locations problem


To calculate the assignment cost of the permutation, the required flows between facilities and the distances between locations are needed.


.. tabularcolumns:: |p{1cm}|p{1cm}|p{1cm}|p{1cm}|

.. csv-table:: flow of the current facilities
   :header: facility `i`, facility `j`, flow( `i`\, `j` )
   :widths: 2, 2, 3

   1, 4, 4
   3, 4, 10  
   3, 1, 8
   2, 1, 6  


.. csv-table:: distances of between locations
   :header: location `i`, location `j`, distances( `i`\, `j` )
   :widths: 2, 2, 3

   1, 2, 42
   1, 3, 30  
   2, 3, 41
   3, 4, 23  


Then, the assignment cost of the permutation can be computed as:

:math:`f(1,4)⋅d(1,2)+f(3,4)⋅d(1,3)+f(1,3)⋅d(2,3)+f(3,2)⋅d(3,4)` 
with result :math:`4⋅42+10⋅30+8⋅41+6⋅23=934`.

Note that this permutation is not the optimal solution.

Mathematical definition
~~~~~~~~~~~~~~~~~~~~~~~

**Sets**

- :math:`N=\{1,2,⋯,n\}`
- :math:`S_n=\phi:N→N` is the set of all permutations

**Parameters**

- :math:`F=(f_{ij})` is an :math:`n×n` matrix where :math:`f_{ij}` is the required flow between facilities :math:`i` and :math:`j`
- :math:`D=(d_{ij})` is an :math:`n×n` matrix where :math:`d_{ij}` is the distance between locations :math:`i` and :math:`j`.

**Optimization Problem**

- :math:`min_{ϕ∈S_n}\sum_{i=1}^{n}{\sum_{j=1}^{n}{f_{ij}⋅d_{\phi(i)\phi(j)}}}`

The assignment of facilities to locations is represented by a permutation :math:`\phi`, where :math:`\phi(i)` is the location to which facility :math:`i` is assigned. Each individual product :math:`f_{ij}⋅d_{\phi(i)\phi(j)}` is the cost of assigning facility :math:`i` to location :math:`\phi(i)` and facility :math:`j` to location :math:`\phi(j)`.

Macop QAP implementation
========================

Let's see how it is possible with the use of the **Macop** package to implement and deal with this QAP instance problem.

Solution structure definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Firstly, we are going to use a type of solution that will allow us to define the structure of our solutions.

The available ``macop.solutions.discrete.CombinatoryIntegerSolution`` type of solution within the Macop package represents exactly what one would wish for. 
I.e. a solution that stores a sequence of integers relative to the size of the problem, the order of which is not sorted.

Let's see an example of its use:

.. code:: python

    from macop.solutions.discrete import CombinatoryIntegerSolution
    
    solution = CombinatoryIntegerSolution.random(10)
    print(solution)


The resulting solution obtained:

.. code:: bash

    Combinatory integer solution [2 9 8 1 7 6 0 4 3 5]


QAP Evaluator
~~~~~~~~~~~~~

Now that we have the structure of our solutions, and the means to generate them, we will seek to evaluate them.

To do this, we need to create a new evaluator specific to our problem and the relative evaluation function:

- :math:`min_{ϕ∈S_n}\sum_{i=1}^{n}{\sum_{j=1}^{n}{f_{ij}⋅d_{\phi(i)\phi(j)}}}`

So we are going to create a class that will inherit from the abstract class ``macop.evalutarors.base.Evaluator``:


.. code:: python

    from macop.evaluators.base import Evaluator

    class QAPEvaluator(Evaluator):
    """QAP evaluator class which enables to compute QAP solution using specific `_data`

    - stores into its `_data` dictionary attritute required measures when computing a QAP solution
    - `_data['F']` matrix of size n x n with flows data between facilities (stored as numpy array)
    - `_data['D']` matrix of size n x n with distances data between locations (stored as numpy array)
    - `compute` method enables to compute and associate a score to a given QAP solution
    """

    def compute(self, solution):
        """Apply the computation of fitness from solution

        Args:
            solution: {Solution} -- QAP solution instance
    
        Returns:
            {float} -- fitness score of solution
        """
        fitness = 0
        for index_i, val_i in enumerate(solution.getdata = )):
            for index_j, val_j in enumerate(solution.getdata = )):
                fitness += self._data['F'][index_i, index_j] * self._data['D'][val_i, val_j]

        return fitness

The cost function for the quadratic problem is now well defined.

.. tip::
    The class proposed here, is available in the Macop package ``macop.evaluators.discrete.mono.QAPEvaluator``.

Running algorithm
~~~~~~~~~~~~~~~~~

Now that the necessary tools are available, we will be able to deal with our problem and look for solutions in the search space of our QAP instance.

Here we will use local search algorithms already implemented in **Macop**.

If you are uncomfortable with some of the elements in the code that will follow, you can refer to the more complete **Macop** documentation_ that focuses more on the concepts and tools of the package.

.. code:: python

    # main imports
    import numpy as np

    # module imports
    from macop.solutions.discrete import CombinatoryIntegerSolution
    from macop.evaluators.discrete.mono import QAPEvaluator

    from macop.operators.discrete.mutators import SimpleMutation

    from macop.policies.classicals import RandomPolicy

    from macop.algorithms.mono import IteratedLocalSearch as ILS
    from macop.algorithms.mono import HillClimberFirstImprovment

    # usefull instance data
    n = 100
    qap_instance_file = 'qap_instance.txt'

    # default validator (check the consistency of our data, i.e. only unique element)
    def validator(solution):
        if len(list(solution.getdata = ))) > len(set(list(solution.getdata = )))):
            print("not valid")
            return False
        return True

    # define init random solution
    def init():
        return CombinatoryIntegerSolution.random(n, validator)

    # load qap instance
    with open(qap_instance_file, 'r') as f:
        file_data = f.readlines()
        print(f'Instance information {file_data[0]}')

        D_lines = file_data[1:n + 1]
        D_data = ''.join(D_lines).replace('\n', '')

        F_lines = file_data[n:2 * n + 1]
        F_data = ''.join(F_lines).replace('\n', '')

    D_matrix = np.fromstring(D_data, dtype=float, sep=' ').reshape(n, n)
    print(f'D matrix shape: {D_matrix.shape}')
    F_matrix = np.fromstring(F_data, dtype=float, sep=' ').reshape(n, n)
    print(f'F matrix shape: {F_matrix.shape}')

    # only one operator here
    operators = [SimpleMutation()]

    # random policy even if list of solution has only one element
    policy = RandomPolicy(operators)

    # use of loaded data from QAP instance
    evaluator = QAPEvaluator(data={'F': F_matrix, 'D': D_matrix})

    # passing global evaluation param from ILS
    hcfi = HillClimberFirstImprovment(init, evaluator, operators, policy, validator, maximise=False, verbose=True)
    algo = ILS(init, evaluator, operators, policy, validator, localSearch=hcfi, maximise=False, verbose=True)

    # run the algorithm
    bestSol = algo.run(10000, ls_evaluations=100)

    print('Solution for QAP instance score is {}'.format(evaluator.compute(bestSol)))


QAP problem solving is now possible with **Macop**. As a reminder, the complete code is available in the qapExample.py_ file.

.. _qapExample.py: https://github.com/jbuisine/macop/blob/master/examples/qapExample.py
.. _documentation: https://jbuisine.github.io/macop/_build/html/documentationsQAP Problem instance generation
===============================

To define our quadratic assignment problem instance, we will use the available mQAP_ multi-objective quadratic problem generator. 

Genration of the instance
~~~~~~~~~~~~~~~~~~~~~~~~~

We will limit ourselves here to a single objective for the purposes of this example. The file **makeQAPuni.cc**, will be used to generate the instance.

.. code:: bash

    g++ makeQAPuni.cc -o mQAPGenerator
    ./mQAPGenerator -n 100 -k 1 -f 30 -d 80 -s 42 > qap_instance.txt

with the following parameters:

- **-n** positive integer: number of facilities/locations;
- **-k** positive integer: number of objectives;
- **-f** positive integer: maximum flow between facilities;
- **-d** positive integer: maximum distance between locations;
- **-s** positive long: random seed.

The generated qap_instance.txt_ file contains the two matrices :math:`F` and :math:`D` and define our instance problem.

.. _mQAP: https://www.cs.bham.ac.uk/~jdk/mQAP/

.. _qap_instance.txt: https://github.com/jbuisine/macop/blob/master/examples/instances/qap/qap_instance.txt


Load data instance
~~~~~~~~~~~~~~~~~~


We are now going to load this instance via a Python code which will be useful to us later on:

.. code:: Python

    qap_instance_file = 'qap_instance.txt'

    n = 100 # the instance size

    with open(qap_instance_file, 'r') as f:
        file_data = f.readlines()
        print(f'Instance information {file_data[0]}')

        D_lines = file_data[1:n + 1]
        D_data = ''.join(D_lines).replace('\n', '')

        F_lines = file_data[n:2 * n + 1]
        F_data = ''.join(F_lines).replace('\n', '')

    D_matrix = np.fromstring(D_data, dtype=float, sep=' ').reshape(n, n)
    print(f'D matrix shape: {D_matrix.shape}')
    F_matrix = np.fromstring(F_data, dtype=float, sep=' ').reshape(n, n)
    print(f'F matrix shape: {F_matrix.shape}')

.. note::
    As we know the size of our instance and the structure of the document, it is quite quick to look for the lines related to the :math:`F` and :math:`D` matrices.Quadratric Assignment Problem
===============================

This example will deal with the use of the **Macop** package in relation to a quadratic assignment problem (QAP). We will use a known example of this problem to associate a set of facilities (:math:`F`) to a set of locations (:math:`L`).

.. image:: ../../_static/examples/qap/factories_qap.png
   :width: 50 %
   :align: center
   :alt: Example of QAP facilities to locations problem

.. toctree::
   :maxdepth: 1
   :numbered:
   :caption: Contents:

   problem
   instance
   implementation


.. note:: 
   The full code for what will be proposed in this example is available: qapExample.py_.

.. _qapExample.py: https://github.com/jbuisine/macop/blob/master/examples/qapExample.pyDocumentation
=============


.. toctree::
   :titlesonly:

   {% for page in pages %}
   {% if page.top_level_object and page.display %}
   {{ page.include_path }}
   {% endif %}
   {% endfor %}

{% if obj.display %}
.. {{ obj.type }}:: {{ obj.name }}
   {%+ if obj.value is not none or obj.annotation is not none %}:annotation:{% if obj.annotation %} :{{ obj.annotation }}{% endif %}{% if obj.value is not none %} = {{ obj.value }}{% endif %}{% endif %}


   {{ obj.docstring|prepare_docstring|indent(3) }}
{% endif %}
{% extends "python/class.rst" %}
{% if not obj.display %}
:orphan:

{% endif %}
:mod:`{{ obj.name }}`
======={{ "=" * obj.name|length }}

.. py:module:: {{ obj.name }}

{% if obj.docstring %}
.. autoapi-nested-parse::

   {{ obj.docstring|prepare_docstring|indent(3) }}

{% endif %}

{% block subpackages %}
{% set visible_subpackages = obj.subpackages|selectattr("display")|list %}
{% if visible_subpackages %}

.. toctree::
   :titlesonly:
   :maxdepth: 3

{% for subpackage in visible_subpackages %}
   {{ subpackage.short_name }}/index.rst
{% endfor %}


{% endif %}
{% endblock %}
{% block submodules %}
{% set visible_submodules = obj.submodules|selectattr("display")|list %}
{% if visible_submodules %}

.. toctree::
   :titlesonly:
   :maxdepth: 1

{% for submodule in visible_submodules %}
   {{ submodule.short_name }}/index.rst
{% endfor %}


{% endif %}
{% endblock %}
{% block content %}
{% if obj.all is not none %}
{% set visible_children = obj.children|selectattr("short_name", "in", obj.all)|list %}
{% elif obj.type is equalto("package") %}
{% set visible_children = obj.children|selectattr("display")|list %}
{% else %}
{% set visible_children = obj.children|selectattr("display")|rejectattr("imported")|list %}
{% endif %}
{% if visible_children %}

{{ "-" * obj.type|length }}---------

{% set visible_classes = visible_children|selectattr("type", "equalto", "class")|list %}
{% set visible_functions = visible_children|selectattr("type", "equalto", "function")|list %}
{% if "show-module-summary" in autoapi_options and (visible_classes or visible_functions) %}
{% block classes %}
{% if visible_classes %}
Classes
~~~~~~~

.. autoapisummary::

{% for klass in visible_classes %}
   {{ klass.id }}
{% endfor %}


{% endif %}
{% endblock %}

{% block functions %}
{% if visible_functions %}
Functions
~~~~~~~~~

.. autoapisummary::

{% for function in visible_functions %}
   {{ function.id }}
{% endfor %}


{% endif %}
{% endblock %}
{% endif %}
{% for obj_item in visible_children %}
{{ obj_item.rendered|indent(0) }}
{% endfor %}
{% endif %}
{% endblock %}
{% if obj.display %}
.. function:: {{ obj.short_name }}({{ obj.args }}){% if obj.return_annotation is not none %} -> {{ obj.return_annotation }}{% endif %}

   {% if sphinx_version >= (2, 1) %}
   {% for property in obj.properties %}
   :{{ property }}:
   {% endfor %}
   {% endif %}

   {% if obj.docstring %}
   {{ obj.docstring|prepare_docstring|indent(3) }}
   {% else %}
   {% endif %}
{% endif %}
{%- if obj.display %}
{% if sphinx_version >= (2, 1) %}
.. method:: {{ obj.short_name }}({{ obj.args }})
   {% for property in obj.properties %}
   :{{ property }}:
   {% endfor %}

{% else %}
.. {{ obj.method_type }}:: {{ obj.short_name }}({{ obj.args }})
{% endif %}

   {% if obj.docstring %}
   {{ obj.docstring|prepare_docstring|indent(3) }}
   {% endif %}
{% endif %}
{% extends "python/data.rst" %}
{% extends "python/module.rst" %}
{% if obj.display %}
.. py:{{ obj.type }}:: {{ obj.short_name }}{% if obj.args %}({{ obj.args }}){% endif %}


   {% if obj.bases %}
   {% if "show-inheritance" in autoapi_options %}
   Bases: {% for base in obj.bases %}:class:`{{ base }}`{% if not loop.last %}, {% endif %}{% endfor %}
   {% endif %}


   {% if "show-inheritance-diagram" in autoapi_options and obj.bases != ["object"] %}
   .. autoapi-inheritance-diagram:: {{ obj.obj["full_name"] }}
      :parts: 1
      {% if "private-members" in autoapi_options %}:private-bases:{% endif %}

   {% endif %}
   {% endif %}
   {% if obj.docstring %}
   {{ obj.docstring|prepare_docstring|indent(3) }}
   {% endif %}
   {% if "inherited-members" in autoapi_options %}
   {% set visible_classes = obj.classes|selectattr("display")|list %}
   {% else %}
   {% set visible_classes = obj.classes|rejectattr("inherited")|selectattr("display")|list %}
   {% endif %}
   {% for klass in visible_classes %}
   {{ klass.rendered|indent(3) }}
   {% endfor %}
   {% if "inherited-members" in autoapi_options %}
   {% set visible_attributes = obj.attributes|selectattr("display")|list %}
   {% else %}
   {% set visible_attributes = obj.attributes|rejectattr("inherited")|selectattr("display")|list %}
   {% endif %}
   {% for attribute in visible_attributes %}
   {{ attribute.rendered|indent(3) }}
   {% endfor %}
   {% if "inherited-members" in autoapi_options %}
   {% set visible_methods = obj.methods|selectattr("display")|list %}
   {% else %}
   {% set visible_methods = obj.methods|rejectattr("inherited")|selectattr("display")|list %}
   {% endif %}
   {% for method in visible_methods %}
   {{ method.rendered|indent(3) }}
   {% endfor %}
{% endif %}
