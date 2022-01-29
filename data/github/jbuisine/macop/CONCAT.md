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
