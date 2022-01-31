# AgentPy - Agent-based modeling in Python

[![PyPI](https://img.shields.io/pypi/v/agentpy)](https://pypi.org/project/agentpy/)
[![GitHub](https://img.shields.io/github/license/joelforamitti/agentpy)](https://github.com/JoelForamitti/agentpy/blob/master/LICENSE)
[![Build Status](https://travis-ci.com/JoelForamitti/agentpy.svg?branch=master)](https://travis-ci.com/JoelForamitti/agentpy)
[![Documentation Status](https://readthedocs.org/projects/agentpy/badge/?version=latest)](https://agentpy.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/JoelForamitti/agentpy/branch/master/graph/badge.svg?token=NTW99HNGB0)](https://codecov.io/gh/JoelForamitti/agentpy)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03065/status.svg)](https://doi.org/10.21105/joss.03065)

AgentPy is an open-source library for the development and analysis of agent-based models in Python. 
The framework integrates the tasks of model design, interactive simulations, numerical experiments, 
and data analysis within a single environment. The package is optimized for interactive computing 
with [IPython](http://ipython.org/), [IPySimulate](https://github.com/JoelForamitti/ipysimulate), and [Jupyter](https://jupyter.org/). 

Please cite this software as follows:

    Foramitti, J., (2021). AgentPy: A package for agent-based modeling in Python. 
    Journal of Open Source Software, 6(62), 3065, https://doi.org/10.21105/joss.03065

**Installation:** `pip install agentpy`

**Documentation:** https://agentpy.readthedocs.io

**JOSS publication:** https://doi.org/10.21105/joss.03065

**Discussion forum:** https://github.com/JoelForamitti/agentpy/discussions

**Tutorials and examples:** https://agentpy.readthedocs.io/en/latest/model_library.html

**Comparison with other frameworks**: https://agentpy.readthedocs.io/en/latest/comparison.html---
title: 'AgentPy: A package for agent-based modeling in Python'
tags:
  - Agent-based modeling
  - Complex systems
  - Computer simulation
  - Interactive computing
  - Python
authors:
  - name: Joël Foramitti
    orcid: 0000-0002-4828-7288
    affiliation: "1, 2"
affiliations:
 - name: Institute of Environmental Science and Technology, Universitat Autònoma de Barcelona, Spain
   index: 1
 - name: Institute for Environmental Studies, Vrije Universiteit Amsterdam, The Netherlands
   index: 2
date: 16.01.2020
bibliography: paper.bib
---

# Introduction

Agent-based models allow for computer simulations based on the autonomous behavior of heterogeneous agents. They are used to generate and understand the emergent dynamics of complex systems, with applications in fields like ecology [@DeAngelis2019], cognitive sciences [@Madsen2019], management [@North2007], policy analysis [@Castro2020], economics [@Arthur2021; @Farmer2009], and sociology [@Bianchi2015].

AgentPy is an open-source library for the development and analysis of agent-based models. It aims to provide an intuitive syntax for the creation of models together with advanced tools for scientific applications. The framework is written in Python 3, and optimized for interactive computing with [IPython](http://ipython.org/) and [Jupyter](https://jupyter.org/). A reference of all features as well as a model library with tutorials and examples can be found in the [documentation](https://agentpy.readthedocs.io/).[^1]

# Statement of Need

There are numerous modeling and simulation tools for agent-based models, each with their own particular focus and style [@Abar2017]. Notable examples are [NetLogo](https://ccl.northwestern.edu/netlogo/) [@Netlogo], which is written in Scala/Java and has become the most established tool in the field; and [Mesa](https://mesa.readthedocs.io/) [@Mesa2015], a more recent framework that has popularized the development of agent-based models in Python.

AgentPy's main distinguishing feature is that it integrates the many different tasks of agent-based modeling within a single environment for interactive computing. This includes the creation of custom agent and model types, interactive simulations (\autoref{fig:interactive}) similar to the traditional NetLogo interface, numeric experiments over multiple runs, and the subsequent data analysis of the output. All of these can be performed within a [Jupyter Notebook](https://jupyter.org/).

The software is further designed for scientific applications, and includes tools for parameter sampling (similar to NetLogo's BehaviorSpace), Monte Carlo experiments, random number generation, parallel computing, and sensitivity analysis. Beyond these built-in features, AgentPy is also designed for compatibility with established Python libraries like [EMA Workbench](https://emaworkbench.readthedocs.io/), [NetworkX](https://networkx.org/), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [SALib](https://salib.readthedocs.io/), [SciPy](https://www.scipy.org/), and [seaborn](https://seaborn.pydata.org/).

![An interactive simulation of Schelling's segregation model in a Jupyter Notebook.\label{fig:interactive}](ips_segregation.png){ width=80% }

# Basic structure

The AgentPy framework follows a nested structure that is illustrated in \autoref{fig:structure}. The basic building blocks are the agents, which can be placed within (multiple) environments with different topologies such as a network, a spatial grid, or a continuous space. Models are used to initiate these objects, perform a simulation, and record data. Experiments can run a model over multiple iterations and parameter combinations. The resulting output data can then be saved and re-arranged for analysis and visualization.

![Nested structure of the AgentPy framework.\label{fig:structure}](structure.png){ width=80% }

# Model example

The following code shows an example of a simple model that explores the distribution of wealth under a randomly trading population of agents. The original version of this model was written in Mesa, allowing for a comparison of the syntax between the two frameworks.[^3] To start, we import the AgentPy library as follows:

```python
import agentpy as ap
```

We then define a new type of [`Agent`](https://agentpy.readthedocs.io/en/stable/reference_agents.html). The method [`setup`](https://agentpy.readthedocs.io/en/stable/reference_agents.html#agentpy.Agent.setup) will be called automatically at the agent's creation. Each agent starts with one unit of wealth. `wealth_transfer` will be called by the model during each time-step. When called, the agent randomly selects a trading partner and hands them one unit of their wealth, given that they have one to spare. 

```python
class WealthAgent(ap.Agent):

    def setup(self):
        self.wealth = 1

    def wealth_transfer(self):
        if self.wealth > 0:
            partner = self.model.agents.random()
            partner.wealth += 1
            self.wealth -= 1
```

Next, we define a [`Model`](https://agentpy.readthedocs.io/en/stable/reference_model.html). The method [`setup`](https://agentpy.readthedocs.io/en/stable/reference_model.html#agentpy.Model.setup) is called at the beginning the simulation, [`step`](https://agentpy.readthedocs.io/en/stable/reference_model.html#agentpy.Model.step) is called during each time-step, and [`end`](https://agentpy.readthedocs.io/en/stable/reference_model.html#agentpy.Model.end) is called after the simulation has finished. An [`AgentList`](https://agentpy.readthedocs.io/en/stable/reference_sequences.html) is used to create a set of agents that can then be accessed as a group. The attribute `p` is used to access the model's parameters. And the method [`record`](https://agentpy.readthedocs.io/en/stable/reference_agents.html#agentpy.Agent.record) is used to store data for later analysis.

```python
class WealthModel(ap.Model):

    def setup(self):
        self.agents = ap.AgentList(self, self.p.n, WealthAgent)

    def step(self):
        self.agents.wealth_transfer()

    def end(self):
        self.agents.record('wealth')
```

To run a simulation, a new instance of the model is created with a dictionary of parameters.
While the parameter `n` is used in the model's setup, the parameter `steps` automatically defines the maximum number of time-steps. Alternatively, the simulation could also be stopped with [`Model.stop`](https://agentpy.readthedocs.io/en/stable/reference_model.html#agentpy.Model.stop). To perform the actual simulation, one can use [`Model.run`](https://agentpy.readthedocs.io/en/stable/reference_model.html#agentpy.Model.run).

```python
parameters = {'n': 100, 'steps': 100}
model = MoneyModel(parameters)
results = model.run()
```

Parameters can also be defined as ranges and used to generate a [`Sample`](https://agentpy.readthedocs.io/en/stable/reference_sample.html).
This sample can then be used to initiate an [`Experiment`](https://agentpy.readthedocs.io/en/stable/reference_experiment.html) that can repeatedly run the model over multiple parameter combinations and iterations. In the following example, the parameter `n` is varied from 1 to 100 and the simulation is repeated 10 times for each value of `n`.

```python
parameters = {'n': ap.IntRange(1, 100), 'steps': 100}
sample = ap.Sample(parameters, n=10)
exp = ap.Experiment(MoneyModel, sample, iterations=10, record=True)
results = exp.run()
```

The output of both models and experiments is given as a [`DataDict`](https://agentpy.readthedocs.io/en/stable/reference_data.html) with tools to save, arrange, and analyse data. Here, we use the seaborn library to display a histogram of the experiment's output. The plot is presented in \autoref{fig:boltzmann}. It shows that the random interaction of the agents creates an inequality of wealth that resembles a Boltzmann-Gibbs distribution. 

```python
import seaborn as sns
sns.histplot(data=results.variables.MoneyAgent, binwidth=1)
```

![Histogram of the agents' wealth in the model example.\label{fig:boltzmann}](boltzmann.pdf){ width=60% }

More examples - including spatial environments, networks, stochastic processes, interactive simulations (see \autoref{fig:interactive}), animations, and sensitivity analysis - can be found in the [model library](https://agentpy.readthedocs.io/en/stable/model_library.html) and [user guides](https://agentpy.readthedocs.io/en/stable/guide.html) of the documentation. For questions and ideas, please visit the [discussion forum](https://github.com/JoelForamitti/agentpy/discussions).[^4]

[^1]: Link to the AgentPy documentation: [https://agentpy.readthedocs.io](https://agentpy.readthedocs.io)
[^3]: For a direct comparison, see: [https://agentpy.readthedocs.io/en/stable/comparison.html](https://agentpy.readthedocs.io/en/stable/comparison.html)
[^4]: Link to the AgentPy dicussion forum: [https://github.com/JoelForamitti/agentpy/discussions](https://github.com/JoelForamitti/agentpy/discussions)

# Acknowledgements

This study has received funding through an ERC Advanced Grant from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement n° 741087). I thank Jeroen C.J.M van den Bergh, Ivan Savin, Martí Bosch, and James Millington for their helpful comments.

# References.. currentmodule:: agentpy

=====
Other
=====

.. autoclass:: AttrDict
    :members:

.. currentmodule:: agentpy

==========
Contribute
==========

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given. You can contribute in many ways:

Types of contributions
----------------------

Report bugs
~~~~~~~~~~~

Report bugs at https://github.com/JoelForamitti/agentpy/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues and `discussion forum <https://github.com/JoelForamitti/agentpy/discussions>`_ for features.
Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

Write documentation
~~~~~~~~~~~~~~~~~~~

Agentpy could always use more documentation, whether as part of the
official agentpy docs, in docstrings, or even on the web in blog posts,
articles, and such. Contributions of clear and simple demonstration models for the :doc:`model_library`
that illustrate a particular application are also very welcome.

Submit feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to write in the agentpy discussion forum at https://github.com/JoelForamitti/agentpy/discussions.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

How to contribute
-----------------

Ready to contribute? Here's how to set up `agentpy` for local development.

1. Fork the `agentpy` repository on GitHub: https://github.com/JoelForamitti/agentpy
2. Clone your fork locally:

.. code-block:: console

    $ git clone git@github.com:your_name_here/agentpy.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development:

.. code-block:: console

    $ mkvirtualenv agentpy
    $ cd agentpy/
    $ pip install -e .['dev']

4. Create a branch for local development:

.. code-block:: console

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass the tests
   and that the new features are covered by the tests:

.. code-block:: console

    $ coverage run -m pytest
    $ coverage report

6. Commit your changes and push your branch to GitHub:

.. code-block:: console

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull request guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
   For more information, check out the tests directory and https://docs.pytest.org/.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to docs/changelog.rst.
3. The pull request should pass the automatic tests on travis-ci. Check
   https://travis-ci.com/JoelForamitti/agentpy/pull_requests
   and make sure that the tests pass for all supported Python versions.
=============
Model Library
=============

Welcome to the agentpy model library.
Below you can find a set of demonstrations on how the package can be used.
All of the models are provided as interactive `Jupyter Notebooks <https://jupyter.org/>`_
that can be downloaded and experimented with.

.. :caption: Contents:
.. toctree::
   :caption: Models
   :maxdepth: 1

   agentpy_wealth_transfer
   agentpy_virus_spread
   agentpy_flocking
   agentpy_segregation
   agentpy_forest_fire
   agentpy_button_network.. currentmodule:: agentpy

======
Agents
======

Agent-based models can contain multiple agents of different types.
This module provides a base class :class:`Agent`
that is meant to be used as a template to create custom agent types.
Initial variables should be defined by overriding :func:`Agent.setup`.

.. autoclass:: Agent
    :members:
    :inherited-members:.. currentmodule:: agentpy

==========================
Graph topologies (Network)
==========================

.. autoclass:: Network
    :members:
    :inherited-members:

.. autoclass:: AgentNode
    :members:
    :inherited-members:

.. currentmodule:: agentpy

=============
API Reference
=============

.. :caption: Contents:
.. toctree::
   :maxdepth: 3

   reference_model
   reference_agents
   reference_sequences
   reference_environments
   reference_sample
   reference_experiment
   reference_data
   reference_visualization
   reference_examples
   reference_other.. currentmodule:: agentpy

==================
Agent-based models
==================

The :class:`Model` contains all objects
and defines the procedures of an agent-based simulation.
It is meant as a template for custom model classes that
override the `custom procedure methods`_.

.. autoclass:: Model

Simulation tools
################

.. automethod:: Model.run
.. automethod:: Model.stop

.. _custom procedure methods:

Custom procedures
#################

.. automethod:: Model.setup
.. automethod:: Model.step
.. automethod:: Model.update
.. automethod:: Model.end

Data collection
###############

.. automethod:: Model.record
.. automethod:: Model.report

Conversion
##########

.. automethod:: Model.as_function.. currentmodule:: agentpy

============
Environments
============

Environments are objects in which agents can inhabit a specific position.
The connection between positions is defined by the environment's
topology. There are currently three types:

- :class:`Grid` n-dimensional spatial topology with discrete positions.
- :class:`Space` n-dimensional spatial topology with continuous positions.
- :class:`Network` graph topology consisting of :class:`AgentNode` and edges.

All three environment classes contain the following methods:

- :func:`add_agents` adds agents to the environment.
- :func:`remove_agents` removes agents from the environment.
- :func:`move_to` changes an agent's position.
- :func:`move_by` changes an agent's position, relative to their current position.
- :func:`neighbors` returns an agent's neighbors within a given distance.

.. toctree::
   :hidden:

   reference_grid
   reference_space
   reference_network.. currentmodule:: agentpy

==========
Comparison
==========

There are numerous modeling and simulation tools for ABMs,
each with their own particular focus and style
(find an overview `here <https://en.wikipedia.org/wiki/Comparison_of_agent-based_modeling_software>`_).
The three main distinguishing features of agentpy are the following:

- Agentpy integrates the multiple tasks of agent-based modeling
  - model design, interactive simulations,
  numerical experiments, and data analysis - within a single environment
  and is optimized for interactive computing with IPython and Jupyter.
- Agentpy is designed for scientific use with experiments over multiple runs.
  It provides tools for parameter sampling (similar to NetLogo's BehaviorSpace),
  Monte Carlo experiments, stochastic processes, parallel computing,
  and sensitivity analysis.
- Agentpy is written in Python, one of the world’s most popular
  programming languages that offers a vast number of tools and libraries for scientific use.
  It is further designed for compatibility with established packages like
  numpy, scipy, networkx, pandas, ema_workbench, seaborn, and SALib.

The main alternative to agentpy in Python is `Mesa <https://mesa.readthedocs.io/>`__.
To allow for an comparison of the syntax,
here are two examples for a simple model of wealth transfer,
both of which realize exactly the same operations.
More information on the two models can be found in the documentation
of each framework (:doc:`Agentpy <agentpy_wealth_transfer>` &
`Mesa <https://mesa.readthedocs.io/en/stable/tutorials/intro_tutorial.html#tutorial-description>`_).

+--------------------------------------------+----------------------------------------------+
|**Agentpy**                                 |**Mesa**                                      |
+--------------------------------------------+----------------------------------------------+
|                                            |                                              |
|.. literalinclude:: agentpy_demo.py         |.. literalinclude:: mesa_demo.py              |
|                                            |                                              |
+--------------------------------------------+----------------------------------------------+

The following table further provides a comparison of the main features of each framework.

==========================  ===================================  ======================================
**Feature**                 **Agentpy**                          **Mesa**
| Containers                | Sequence classes                   | Scheduler classes for
                            | like AgentList and AgentDList      | different activation orders
| Topologies                | Spatial grid, continuous space,    | Spatial grid, continuous space,
                            | network                            | network
| Data recording            | Recording methods for variables    | DataCollector class that can
                            | of agents, environments, and       | collect variables of agents
                            | model; as well as reporters        | and model
| Parameter sampling        | Classes for sample generation
                            | and different types of
                            | parameter ranges
| Multi-run experiments     | Experiment class that supports     | BatchRunner class that supports
                            | multiple iterations, parameter     | multiple iterations and parameter
                            | samples, randomization,            | ranges
                            | and parallel processing
| Output data               | DataDict class to store, save,     | Methods to generate dataframes
                            | load, and re-arrange output data   |
| Visualization             | Gridplots, animations,             | Plots and interactive visualization
                            | and interactive visualization      | in a separate web-server
                            | within Jupyter Notebooks
| Analysis                  | Tools for data arrangement and
                            | sensitivity analysis
==========================  ===================================  ======================================.. currentmodule:: agentpy

=========
Sequences
=========

This module offers various data structures to create and manage groups
of both agents and environments. Which structure best to use
depends on the specific requirements of each model.

- :class:`AgentList` is a list of agentpy objects with
  methods to select and manipulate its entries.
- :class:`AgentDList` is an ordered collection of agentpy objects,
  optimized for removing and looking up objects.
- :class:`AgentSet` is an unordered collection of agents
  that can access agent attributes.
- :class:`AgentIter` and :class:`AgentDListIter` are a list-like iterators
  over a selection of agentpy objects.
- :class:`AttrIter` is a list-like iterator over the attributes of
  each agent in a selection of agentpy objects.

All of these sequence classes can access and manipulate
the methods and variables of their objects as an attribute of the container.
For examples, see :class:`AgentList`.

Containers
##########

.. autoclass:: AgentList
    :members:

.. autoclass:: AgentDList
    :members:

.. autoclass:: AgentSet
    :members:

Iterators
#########

.. autoclass:: AgentIter
    :members:

.. autoclass:: AgentDListIter
    :members:

.. autoclass:: AttrIter
    :members:
.. currentmodule:: agentpy

=================
Parameter samples
=================

Value sets and ranges
#####################

.. autoclass:: Range

.. autoclass:: IntRange

.. autoclass:: Values

Sample generation
#################

.. autoclass:: Sample
    :members:
.. currentmodule:: agentpy

======================
Discrete spaces (Grid)
======================

.. autoclass:: Grid
    :members:
    :inherited-members:

.. autoclass:: GridIter
    :members:
    :inherited-members:.. currentmodule:: agentpy

=========================
Continuous spaces (Space)
=========================

.. autoclass:: Space
    :members:
    :inherited-members:
.. currentmodule:: agentpy

=============
Data analysis
=============

This module offers tools to access, arrange, analyse, and store output data from simulations.
A :class:`DataDict` can be generated by the methods :func:`Model.run`, :func:`Experiment.run`, and :func:`DataDict.load`.

.. autoclass:: DataDict

Data arrangement
################

.. automethod:: DataDict.arrange
.. automethod:: DataDict.arrange_reporters
.. automethod:: DataDict.arrange_variables

Analysis methods
################

.. automethod:: DataDict.calc_sobol

Save and load
#############

.. automethod:: DataDict.save
.. automethod:: DataDict.load.. currentmodule:: agentpy

===========
Experiments
===========

.. autoclass:: Experiment
    :members:===========
User Guides
===========

This section contains interactive notebooks with common applications of the agentpy framework.
If you are interested to add a new article to this guide, please visit :doc:`contributing`.
If you are looking for examples of complete models, take a look at :doc:`model_library`.
To learn how agentpy compares with other frameworks, take a look at :doc:`comparison`.

.. :caption: Contents:
.. toctree::
   :caption: Contents
   :maxdepth: 1

   guide_interactive
   guide_random
   guide_ema.. currentmodule:: agentpy

=============
Visualization
=============

.. autofunction:: animate

.. autofunction:: gridplot.. currentmodule:: agentpy
.. highlight:: shell

============
Installation
============

To install the latest release of agentpy,
run the following command on your console:

.. code-block:: console

	$ pip install agentpy

Dependencies
------------

Agentpy supports Python 3.6 and higher.
The installation includes the following packages:

- `numpy <https://numpy.org>`_ and `scipy <https://docs.scipy.org/>`_, for scientific computing
- `matplotlib <https://matplotlib.org/>`_, for visualization
- `pandas <https://pandas.pydata.org>`_, for data manipulation
- `networkx <https://networkx.org/documentation/>`_, for networks/graphs
- `SALib <https://salib.readthedocs.io/>`_, for sensitivity analysis
- `joblib <https://joblib.readthedocs.io/>`_, for parallel processing

These optional packages can further be useful in combination with agentpy:

- `jupyter <https://jupyter.org/>`_, for interactive computing
- `ipysimulate <https://ipysimulate.readthedocs.io/>`_ >= 0.2.0, for interactive simulations
- `ema_workbench <https://emaworkbench.readthedocs.io/>`_, for exploratory modeling
- `seaborn <https://seaborn.pydata.org/>`_, for statistical data visualization

Development
-----------

The most recent version of agentpy can be cloned from Github:

.. code-block:: console

	$ git clone https://github.com/JoelForamitti/agentpy.git

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ pip install -e

To include all necessary packages for development & testing, you can use:

.. code-block:: console

    $ pip install -e .['dev']

.. _Github repository: https://github.com/JoelForamitti/agentpy
.. currentmodule:: agentpy.examples

========
Examples
========

The following example models are presented in the :doc:`model_library`.

To use these classes, they have to be imported as follows::

    from agentpy.examples import WealthModel

.. autoclass:: WealthModel
.. autoclass:: SegregationModel.. currentmodule:: agentpy

========
Overview
========

This section provides an overview over the main classes and
functions of AgentPy and how they are meant to be used.
For a more detailed description of each element,
please refer to the :doc:`guide` and :doc:`reference`.
Throughout this documentation, AgentPy is imported as follows::

    import agentpy as ap

Structure
#########

The basic structure of the AgentPy framework has four levels:

1. The :class:`Agent` is the basic building block of a model
2. The environment types :class:`Grid`, :class:`Space`, and :class:`Network` contain agents
3. A :class:`Model` contains agents, environments, parameters, and simulation procedures
4. An :class:`Experiment` can run a model multiple times with different parameter combinations

All of these classes are templates that can be customized through the creation of
`sub-classes <https://docs.python.org/3/tutorial/classes.html?highlight=inheritance#inheritance>`_
with their own variables and methods.

Creating models
###############

A custom agent type can be defined as follows::

    class MyAgent(ap.Agent):

        def setup(self):
            # Initialize an attribute with a parameter
            self.my_attribute = self.p.my_parameter

        def agent_method(self):
            # Define custom actions here
            pass

The method :func:`Agent.setup` is meant to be overwritten
and will be called automatically after an agent's creation.
All variables of an agents should be initialized within this method.
Other methods can represent actions that the agent will be able to take during a simulation.

All model objects (including agents, environments, and the model itself)
are equipped with the following default attributes:

- :attr:`model` the model instance
- :attr:`id` a unique identifier number for each object
- :attr:`p` the model's parameters
- :attr:`log` the object's recorded variables

Using the new agent type defined above,
here is how a basic model could look like::

    class MyModel(ap.Model):

        def setup(self):
            """ Initiate a list of new agents. """
            self.agents = ap.AgentList(self, self.p.agents, MyAgent)

        def step(self):
            """ Call a method for every agent. """
            self.agents.agent_method()

        def update(self):
            """ Record a dynamic variable. """
            self.agents.record('my_attribute')

        def end(self):
            """ Repord an evaluation measure. """
            self.report('my_measure', 1)

The simulation procedures of a model are defined by four special methods
that will be used automatically during different parts of a simulation.

- :class:`Model.setup` is called at the start of the simulation (`t==0`).
- :class:`Model.step` is called during every time-step (excluding `t==0`).
- :class:`Model.update` is called after every time-step (including `t==0`).
- :class:`Model.end` is called at the end of the simulation.

If you want to see a basic model like this in action,
take a look at the :doc:`agentpy_wealth_transfer` demonstration in the :doc:`model_library`.

.. _overview_agents:

Agent sequences
###############

The :doc:`reference_sequences` module provides containers for groups of agents.
The main classes are :class:`AgentList`, :class:`AgentDList`, and :class:`AgentSet`,
which come with special methods to access and manipulate whole groups of agents.

For example, when the model defined above calls :func:`self.agents.agent_method`,
it will call the method :func:`MyAgentType.agent_method` for every agent in the model.
Similar commands can be used to set and access variables, or select subsets
of agents with boolean operators.
The following command, for example, selects all agents with an id above one::

    agents.select(agents.id > 1)

Further examples can be found in :doc:`reference_sequences`
and the :doc:`agentpy_virus_spread` demonstration model.

.. _overview_environments:

Environments
############

:doc:`reference_environments` are objects in which agents can inhabit a specific position.
A model can contain zero, one or multiple environments which agents can enter and leave.
The connection between positions is defined by the environment's topology.
There are currently three types:

- :class:`Grid` n-dimensional spatial topology with discrete positions.
- :class:`Space` n-dimensional spatial topology with continuous positions.
- :class:`Network` graph topology consisting of :class:`AgentNode` and edges.

Applications of networks can be found in the demonstration models
:doc:`agentpy_virus_spread` and :doc:`agentpy_button_network`;
spatial grids in :doc:`agentpy_forest_fire` and :doc:`agentpy_segregation`;
and continuous spaces in :doc:`agentpy_flocking`.
Note that there can also be models without environments like in :doc:`agentpy_wealth_transfer`.

Recording data
##############

There are two ways to document data from the simulation for later :ref:`analysis <overview_analysis>`.

The first way is to record dynamic variables,
which can be recorded for each object (agent, environment, or model) and time-step.
They are useful to look at the dynamics of individual or aggregate objects over time
and can be documented by calling the method :meth:`record` for the respective object.
Recorded variables can at run-time with the object's `log` attribute.

The second way is to document reporters,
which represent summary statistics or evaluation measures of a simulation.
In contrast to variables, reporters can be stored only for the model as a whole and only once per run.
They will be stored in a separate dataframe for easy comparison over multiple runs,
and can be documented with the method :meth:`Model.report`.
Reporters can be accessed at run-time via :attr:`Model.reporters`.

.. _overview_simulation:

Running a simulation
####################

To perform a simulation, we initialize a new instance of our model type
with a dictionary of parameters, and then use the function :func:`Model.run`.
This will return a :class:`DataDict` with recorded data from the simulation.
A simple run can be prepared and executed as follows::

    parameters = {
        'my_parameter':42,
        'agents':10,
        'steps':10
    }

    model = MyModel(parameters)
    results = model.run()

A simulation proceeds as follows (see also Figure 1 below):

0. The model initializes with the time-step :attr:`Model.t = 0`.
1. :func:`Model.setup` and :func:`Model.update` are called.
2. The model's time-step is increased by 1.
3. :func:`Model.step` and :func:`Model.update` are called.
4. Step 2 and 3 are repeated until the simulation is stopped.
5. :func:`Model.end` is called.

The simulation of a model can be stopped by one of the following two ways:

1. Calling the :func:`Model.stop` during the simulation.
2. Reaching the time-limit, which be defined as follows:

   - Defining :attr:`steps` in the paramater dictionary.
   - Passing :attr:`steps` as an argument to :func:`Model.run`.

Interactive simulations
#######################

Within a Jupyter Notebook,
AgentPy models can be explored as an interactive simulation
(similar to the traditional NetLogo interface)
using `ipysimulate <https://github.com/JoelForamitti/ipysimulate>`_ and `d3.js <https://d3js.org/>`_.
For more information on this, please refer to :doc:`guide_interactive`.

.. _overview_experiments:

Multi-run experiments
#####################

The :doc:`reference_sample` module provides tools to create a :class:`Sample`
with multiple parameter combinations from a dictionary of ranges.
Here is an example using :class:`IntRange` integer ranges::

    parameters = {
        'my_parameter': 42,
        'agents': ap.IntRange(10, 20),
        'steps': ap.IntRange(10, 20)
    }
    sample = ap.Sample(parameters, n=5)

The class :class:`Experiment` can be used to run a model multiple times.
As shown in Figure 1, it will start with the first parameter combination
in the sample and repeat the simulation for the amount of defined iterations.
After, that the same cycle is repeated for the next parameter combination.

.. figure:: graphics/simulation_flow.png
   :alt: Chain of events in Model and Experiment

   Figure 1: Chain of events in :class:`Model` and :class:`Experiment`.

Here is an example of an experiment with the model defined above.
In this experiment, we use a sample where one parameter is kept fixed
while the other two are varied 5 times from 10 to 20 and rounded to integer.
Every possible combination is repeated 2 times, which results in 50 runs::

    exp = ap.Experiment(MyModel, sample, iterations=2, record=True)
    results = exp.run()

For more applied examples of experiments, check out the demonstration models
:doc:`agentpy_virus_spread`, :doc:`agentpy_button_network`, and :doc:`agentpy_forest_fire`.
An alternative to the built-in experiment class is to use AgentPy models with
the EMA workbench (see :doc:`guide_ema`).

Random numbers
##############

:class:`Model` contains two random number generators:

- :attr:`Model.random` is an instance of :class:`random.Random`
- :attr:`Model.nprandom` is an instance of :class:`numpy.random.Generator`

The random seed for these generators can be set by defining a parameter `seed`.
The :class:`Sample` class has an argument `randomize`
to control whether vary seeds over different parameter combinations.
Similarly, :class:`Experiment` also has an argument `randomize`
to control whether to vary seeds over different iterations.
More on this can be found in :doc:`guide_random`.

.. _overview_analysis:

Data analysis
#############

Both :class:`Model` and :class:`Experiment` can be used to run a simulation,
which will return a :class:`DataDict` with output data.
The output from the experiment defined above looks as follows::

    >>> results
    DataDict {
    'info': Dictionary with 5 keys
    'parameters':
        'constants': Dictionary with 1 key
        'sample': DataFrame with 2 variables and 25 rows
    'variables':
        'MyAgent': DataFrame with 1 variable and 10500 rows
    'reporters': DataFrame with 1 variable and 50 rows
    }

All data is given in a :class:`pandas.DataFrame` and
formatted as `long-form data <https://seaborn.pydata.org/tutorial/data_structure.html>`_
that can easily be used with statistical packages like `seaborn <https://seaborn.pydata.org/>`_.
The output can contain the following categories of data:

- :attr:`info` holds meta-data about the model and simulation performance.
- :attr:`parameters` holds the parameter values that have been used for the experiment.
- :attr:`variables` holds dynamic variables, which can be recorded at multiple time-steps.
- :attr:`reporters` holds evaluation measures that are documented only once per simulation.
- :attr:`sensitivity` holds calculated sensitivity measures.

The :class:`DataDict` provides the following main methods to handle data:

- :func:`DataDict.save` and :func:`DataDict.load` can be used to store results.
- :func:`DataDict.arrange` generates custom combined dataframes.
- :func:`DataDict.calc_sobol` performs a Sobol sensitivity analysis.

Visualization
#############

In addition to the :doc:`guide_interactive`,
AgentPy provides the following functions for visualization:

- :func:`animate` generates an animation that can display output over time.
- :func:`gridplot` visualizes agent positions on a spatial :class:`Grid`.

To see applied examples of these functions, please check out the :doc:`model_library`... currentmodule:: agentpy

========================================
AgentPy - Agent-based modeling in Python
========================================

.. image:: https://img.shields.io/pypi/v/agentpy.svg
    :target: https://pypi.org/project/agentpy/
.. image:: https://img.shields.io/github/license/joelforamitti/agentpy
    :target: https://github.com/JoelForamitti/agentpy/blob/master/LICENSE
.. image:: https://travis-ci.com/JoelForamitti/agentpy.svg?branch=master
    :target: https://travis-ci.com/JoelForamitti/agentpy
.. image:: https://readthedocs.org/projects/agentpy/badge/?version=latest
    :target: https://agentpy.readthedocs.io/en/latest/?badge=latest
.. image:: https://codecov.io/gh/JoelForamitti/agentpy/branch/master/graph/badge.svg?token=NTW99HNGB0
    :target: https://codecov.io/gh/JoelForamitti/agentpy
.. image:: https://joss.theoj.org/papers/10.21105/joss.03065/status.svg
    :target: https://doi.org/10.21105/joss.03065

.. raw:: latex

    \chapter{Introduction}

AgentPy is an open-source library for the development and analysis of agent-based models in Python.
The framework integrates the tasks of model design, interactive simulations, numerical experiments,
and data analysis within a single environment. The package is optimized for interactive computing
with `IPython <http://ipython.org/>`_, `IPySimulate <https://github.com/JoelForamitti/ipysimulate>`_, and `Jupyter <https://jupyter.org/>`_.
If you have questions or ideas for improvements, please visit the
`discussion forum <https://github.com/JoelForamitti/agentpy/discussions>`_.

.. rubric:: Quick orientation

- To get started, please take a look at :doc:`installation` and :doc:`overview`.
- For a simple demonstration, check out the :doc:`agentpy_wealth_transfer` tutorial in the :doc:`model_library`.
- For a detailled description of all classes and functions, refer to :doc:`reference`.
- To learn how agentpy compares with other frameworks, take a look at :doc:`comparison`.
- If you are interested to contribute to the library, see :doc:`contributing`.

.. rubric:: Citation

Please cite this software as follows:

.. code-block:: text

    Foramitti, J., (2021). AgentPy: A package for agent-based modeling in Python.
    Journal of Open Source Software, 6(62), 3065, https://doi.org/10.21105/joss.03065

.. only:: html

    .. rubric:: Table of contents

.. toctree::
   :maxdepth: 2

   installation
   overview
   guide
   model_library
   reference
   comparison
   changelog
   contributing
   about

.. only:: html

    .. rubric:: Indices and tables
    
    * :ref:`genindex`
    * :ref:`search`.. currentmodule:: agentpy

=====
About
=====

Agentpy has been created by Joël Foramitti and is
available under the open-source `BSD 3-Clause <https://github.com/JoelForamitti/agentpy/blob/master/LICENSE>`_ license.
Source files can be found on the `GitHub repository <https://github.com/joelforamitti/agentpy>`_.

Thanks to everyone who has contributed
or supported the developement of this package:

- Jeroen C.J.M. van den Bergh
- Ivan Savin
- James Millington
- Martí Bosch
- Sebastian Benthall
- Bhakti Stephan Onggo

This project has benefited from an ERC Advanced Grant
from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(grant agreement n° 741087).

Parts of this package where created with Cookiecutter_
and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. currentmodule:: agentpy

=========
Changelog
=========

0.1.5 (December 2021)
---------------------

- :func:`Experiment.run` has a new argument 'n_jobs' that allows for
  parallel processing with :func:`joblib.Parallel`.
- Two new methods - :func:`Grid.record_positions`
  and :func:`Space.record_positions` - can be used to record
  agent positions.
- :func:`Model.run` can now continue simulations that have already been run.
  The steps defined in the argument 'steps' now reflect additional steps,
  which will be added to the models current time-step.
  Random number generators will not be re-initialized in this case.
- :func:`animate` has been improved.
  It used to stop the animation one step too early, which has been fixed.
  Two faulty import statements have been corrected.
  And, as above, the argument 'steps' now also reflects additional steps.
- :func:`Grid.add_field` has been fixed. Single values can now be passed.

0.1.4 (September 2021)
----------------------

- :class:`AttrIter` now returns a new :class:`AttrIter` when called as a function.
- :func:`gridplot` now returns an :class:`matplotlib.image.AxesImage`
- :func:`DataDict.save` now
  supports values of type :class:`numpy.bool_`
  and can re-write to existing directories if an existing `exp_id` is passed.
- :func:`DataDict.load` now supports the argument `exp_id = 0`.
- :func:`animate` now supports more than 100 steps.
- :class:`AttrIter` now returns a new :class:`AttrIter` when called as a function.
- :class:`Model` can take a new parameter `report_seed` (default True) that
  indicates whether the seed of the current run should be reported.

0.1.3 (August 2021)
-------------------

- The :class:`Grid` functionality `track_empty` has been fixed
  to work with multiple agents per cell.
- Getting and setting items in :class:`AttrIter` has been fixed.
- Sequences like :class:`AgentList` and :class:`AgentDList`
  no longer accept `args`, only `kwargs`.
  These keyword arguments are forwarded
  to the constructor of the new objects.
  Keyword arguments with sequences of type :class:`AttrIter` will be
  broadcasted, meaning that the first value will be assigned
  to the first object, the second to the second, and so forth.
  Otherwise, the same value will be assigned to all objects.

0.1.2 (June 2021)
-----------------

- The property :attr:`Network.nodes` now returns an :class:`AttrIter`,
  so that network nodes can be assigned to agents as follows::

      self.nw = ap.Network(self)
      self.agents = ap.AgentList(self, 10)
      self.nw.add_agents(self.agents)
      self.agents.node = self.nw.nodes

- :class:`AgentIter` now requires the model to be passed upon creation
  and has two new methods :func:`AgentIter.to_list` and
  :func:`AgentIter.to_dlist` for conversion between sequence types.
- Syntax highlighting in the documentation has been fixed.

0.1.1 (June 2021)
-----------------

- Marked release for the upcoming JOSS publication of AgentPy.
- Fixed :func:`Grid.move_to`: Agents can now move to their current position.

0.1.0 (May 2021)
----------------

This update contains major revisions of most classes and methods in the
library, including new features, better performance, and a more coherent syntax.
The most important API changes are described below.

Object creation
...............

The methods :func:`add_agents`, :func:`add_env`, etc. have been removed.
Instead, new objects are now created directly or through :doc:`reference_sequences`.
This allows for more control over data structures (see next point) and attribute names.
For example::

    class Model(ap.Model):
        def setup(self):
            self.single_agent = ap.Agent()  # Create a single agent
            self.agents = ap.AgentList(self, 10)  # Create a sequence of 10 agents
            self.grid = ap.Grid(self, (5, 5))  # Create a grid environment

Data structures
...............

The new way of object creation makes it possible to choose specific data structures for different groups of agents.
In addition to :class:`AgentList`, there is a new sequence type
:class:`AgentDList` that provides increased performance
for the lookup and deletion of agents.
It also comes with a method :func:`AgentDList.buffer`
that allows for safe deletion of agents
from the list while it is iterated over

:class:`AttrList` has been replaced by :class:`AttrIter`.
This improves performance and makes it possible to change
agent attributes by setting new values to items in the attribute list (see
:class:`AgentList` for an example). In most other ways, the class still behaves like a normal list.
There are also two new classes :class:`AgentIter` and :class:`AgentDListIter` that are returned by some of the library's methods.

Environments
............

The three environment classes have undergone a major revision.
The :func:`add_agents` functions have been extended with new features
and are now more consistent between the three environment classes.
The method :func:`move_agents` has been replaced by :func:`move_to` and :func:`move_by`.
:class:`Grid` is now defined as a structured numpy array
that can hold field attributes per position in addition to agents,
and can be customized with the arguments `torus`, `track_empty`, and `check_border`.
:func:`gridplot` has been adapted to support this new numpy structure.
:class:`Network` now consists of :class:`AgentNode` nodes that can hold multiple agents per node, as well as node attributes.

Environment-agent interaction
.............................

The agents' `env` attribute has been removed.
Instead, environments are manually added as agent attributes,
giving more control over the attribute name in the case of multiple environments.
For example, agents in an environment can be set up as follows::

    class Model(ap.Model):
        def setup(self):
            self.agents = ap.AgentList(self, 10)
            self.grid = self.agents.mygrid = ap.Grid(self, (10, 10))
            self.grid.add_agents(self.agents)

The agent methods `move_to`, `move_by`, and `neighbors` have also been removed.
Instead, agents can access these methods through their environment.
In the above example, a given agent `a` could for example access their position
through `a.mygrid.positions[a]` or their neighbors through calling `a.mygrid.neighbors(a)`.

Parameter samples
.................

Variable parameters can now be defined with the three new classes
:class:`Range` (for continuous parameter ranges), :class:`IntRange` (for integer parameter ranges), and :class:`Values` (for pre-defined of discrete parameter values).
Parameter dictionaries with these classes can be used to create samples,
but can also be passed to a normal model, which will then use default values.
The sampling methods :func:`sample`, :func:`sample_discrete`, and :func:`sample_saltelli`
have been removed and integrated into the new class :class:`Sample`,
which comes with additional features to create new kinds of samples.

Random number generators
........................

:class:`Model` now contains two random number generators `Model.random` and `Model.nprandom`
so that both standard and numpy random operations can be used.
The parameter `seed` can be used to initialize both generators.
:class:`Sample` has an argument `randomize` to vary seeds over parameter samples.
And :class:`Experiment` has a new argument `randomize` to control whether
to vary seeds over different iterations.
More on this can be found in :doc:`guide_random`.

Data analysis
.............

The structure of output data in :class:`DataDict` has been changed.
The name of `measures` has been changed to `reporters`.
Parameters are now stored in the two categories `constants` and `sample`.
Variables are stored in separate dataframes based on the object type.
The dataframe's index is now separated into `sample_id` and `iteration`.
The function :func:`sensitivity_sobol` has been removed and is replaced
by the method :func:`DataDict.calc_sobol`.

Interactive simulations
.......................

The method :func:`Experiment.interactive` has been removed and is replaced
by an interactive simulation interface that is being developed in the separate
package `ipysimulate <https://github.com/JoelForamitti/ipysimulate>`_.
This new package provides interactive javascript widgets with parameter sliders
and live plots similar to the traditional NetLogo interface.
Examples can be found in :doc:`guide_interactive`.

0.0.7 (March 2021)
------------------

Continuous space environments
.............................

A new environment type :class:`Space` and method :func:`Model.add_space`
for agent-based models with continuous space topologies has been added.
There is a new demonstration model :doc:`agentpy_flocking` in the model library,
which shows how to simulate the flocking behavior of animals
and demonstrates the use of the continuous space environment.

Random number generators
........................

:class:`Model` has a new property :obj:`Model.random`, which returns the
models' random number generator of type :func:`numpy.random.Generator`.
A custom seed can be set for :func:`Model.run` and :func:`animate`
by either passing an argument or defining a parameter :attr:`seed`.
All methods with stochastic elements like :func:`AgentList.shuffle`
or :func:`AgentList.random` now take an optional argument `generator`,
with the model's main generator being used if none is passed.
The function :func:`AgentList.random` now uses :func:`numpy.random.Generator.choice`
and has three new arguments 'replace', 'weights', and 'shuffle'.
More information with examples can be found in the API reference
and the new user guide :doc:`guide_random`.

Other changes
.............

* The function :func:`sensitivity_sobol` now has an argument :attr:`calc_second_order` (default False).
  If True, the function will add second-order indices to the output.
* The default value of :attr:`calc_second_order` in :func:`sample_saltelli`
  has also been changed to False for consistency.
* For consistency with :class:`Space`,
  :class:`Grid` no longer takes an integer as argument for 'shape'.
  A tuple with the lengths of each spatial dimension has to be passed.
* The argument 'agents' has been removed from :class:`Environment`.
  Agents have to be added through :func:`Environment.add_agents`.

Fixes
.....

* The step limit in :func:`animate` is now the same as in :func:`Model.run`.
* A false error message in :func:`DataDict.save` has been removed.

0.0.6 (January 2021)
--------------------

* A new demonstration model :doc:`agentpy_segregation` has been added.
* All model objects now have a unique id number of type :class:`int`.
  Methods that take an agent or environment as an argument
  can now take either the instance or id of the object.
  The :attr:`key` attribute of environments has been removed.
* Extra keyword arguments to :class:`Model` and :class:`Experiment`
  are now forwarded to :func:`Model.setup`.
* :func:`Model.run` now takes an optional argument `steps`.
* :class:`EnvDict` has been replaced by :class:`EnvList`,
  which has the same functionalities as :class:`AgentList`.
* Model objects now have a property :attr:`env`
  that returns the first environment of the object.
* Revision of :class:`Network`.
  The argument `map_to_nodes` has been removed from :func:`Network.add_agents`.
  Instead, agents can be mapped to nodes by passing an AgentList to the agents argument of :func:`Model.add_network`.
  Direct forwarding of attribute calls to :attr:`Network.graph` has been
  removed to avoid confusion.
* New and revised methods for :class:`Grid`:

  * :func:`Agent.move_to` and :func:`Agent.move_by` can be used to move agents.
  * :func:`Grid.items` returns an iterator of position and agent tuples.
  * :func:`Grid.get_agents` returns agents in selected position or area.
  * :func:`Grid.position` returns the position coordinates for an agent.
  * :func:`Grid.positions` returns an iterator of position coordinates.
  * :func:`Grid.attribute` returns a nested list with values of agent attributes.
  * :func:`Grid.apply` returns nested list with return values of a custom function.
  * :func:`Grid.neighbors` has new arguments `diagonal` and `distance`.

* :func:`gridplot` now takes a grid of values as an input and can convert them to rgba.
* :func:`animate` now takes a model instance as an input instead of a class and parameters.
* :func:`sample` and :func:`sample_saltelli` will now return integer values for parameters
  if parameter ranges are given as integers. For float values,
  a new argument `digits` can be passed to round parameter values.
* The function :func:`interactive` has been removed, and is replaced by the
  new method :func:`Experiment.interactive`.
* :func:`sobol_sensitivity` has been changed to :func:`sensitivity_sobol`.

0.0.5 (December 2020)
---------------------

* :func:`Experiment.run` now supports parallel processing.
* New methods :func:`DataDict.arrange_variables` and :func:`DataDict.arrange_measures`,
  which generate a dataframe of recorded variables or measures and varied parameters.
* Major revision of :func:`DataDict.arrange`, see new description in the documentation.
* New features for :class:`AgentList`: Arithmethic operators can now be used with :class:`AttrList`.

0.0.4 (November 2020)
---------------------

First documented release.
