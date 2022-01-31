# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

## [v0.5.1] - 18/09/2021

### Added

- locked routes checks

### Fixed

- issues #102, #103, #105 - #110

## [v0.5.0] - 06/06/2021

### Added

- `heuristic_only` option
- `use_all_vehicles` option

### Fixed

- set covering constraints in last MIP with = sign

## [v0.4.0] - 13/05/2021

### Added

- `num_vehicles` option with `periodic` option

### Changed

- cspy 1.0.0
- node load when simultaneous distribution and collection (#79) is now accurate

### Fixed

- issues #79, #82, #84, #86

## [v0.3.0] - 10/11/2020

### Added

- JOSS paper
- Periodic CVRP scheduling option
- Initial solution for CVRP computed with Greedy Algorithm
- Diving heuristic (controlled with new parameter in `VehicleRoutingProblem.solve`)
- Hyper-heuristic pricing strategy option `pricing_strategy="Hyper"`.
- Jupyter notebooks with hyper-heuristics experiments (one to be updated soon).
- Paragraph to the paper with the hyper-heuristic explanation and citations.

### Changed

- Master problem formulation to column-based
- Benchmark tests

## [v0.2.0] - 07/06/2020

### Added

- Mixed fleet option
- Greedy randomized pricing option
- Stabilization with Interior Points
- Diving heuristic WIP

### Changed

- Pricing strategy names


[Unreleased]: https://github.com/Kuifje02/vrpy
[v0.2.0]: https://github.com/Kuifje02/vrpy/releases/tag/v0.2.0
[v0.3.0]: https://github.com/Kuifje02/vrpy/releases/tag/v0.3.0
[v0.4.0]: https://github.com/Kuifje02/vrpy/releases/tag/v0.4.0
[v0.5.0]: https://github.com/Kuifje02/vrpy/releases/tag/v0.5.0
[v0.5.1]: https://github.com/Kuifje02/vrpy/releases/tag/v0.5.1
[![CircleCI](https://circleci.com/gh/Kuifje02/vrpy.svg?style=svg)](https://circleci.com/gh/Kuifje02/vrpy)
[![codecov](https://codecov.io/gh/Kuifje02/vrpy/branch/master/graph/badge.svg)](https://codecov.io/gh/Kuifje02/vrpy)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/6f27b9ccd1c2446aa1dba15e701aa9b0)](https://app.codacy.com/manual/Kuifje02/vrpy?utm_source=github.com&utm_medium=referral&utm_content=Kuifje02/vrpy&utm_campaign=Badge_Grade_Dashboard)
[![Python 3.8](https://img.shields.io/badge/python-3.6|3.7|3.8-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Documentation Status](https://readthedocs.org/projects/vrpy/badge/?version=latest)](https://vrpy.readthedocs.io/en/latest/?badge=master)
[![status](https://joss.theoj.org/papers/77c3aa9b9cb3ff3d5c32d253922ad390/status.svg)](https://joss.theoj.org/papers/77c3aa9b9cb3ff3d5c32d253922ad390)

# VRPy

VRPy is a python framework for solving Vehicle Routing Problems (VRP) including:

-   the Capacitated VRP (CVRP),
-   the CVRP with resource constraints,
-   the CVRP with time windows (CVRPTW),
-   the CVRP with simultaneous distribution and collection (CVRPSDC),
-   the CVRP with heterogeneous fleet (HFCVRP).

Check out the [docs](https://vrpy.readthedocs.io/en/latest/) to find more variants and options.

## Simple example

```python
from networkx import DiGraph
from vrpy import VehicleRoutingProblem

# Define the network
G = DiGraph()
G.add_edge("Source",1,cost=1,time=2)
G.add_edge("Source",2,cost=2,time=1)
G.add_edge(1,"Sink",cost=0,time=2)
G.add_edge(2,"Sink",cost=2,time=3)
G.add_edge(1,2,cost=1,time=1)
G.add_edge(2,1,cost=1,time=1)

# Define the customers demands
G.nodes[1]["demand"] = 5
G.nodes[2]["demand"] = 4

# Define the Vehicle Routing Problem
prob = VehicleRoutingProblem(G, load_capacity=10, duration=5)

# Solve and display solution value
prob.solve()
print(prob.best_value)
3
print(prob.best_routes)
{1: ["Source",2,1,"Sink"]}
```

## Install

```sh
pip install vrpy
```

## Requirements

[cspy](https://pypi.org/project/cspy/)

[NetworkX](https://pypi.org/project/networkx/)

[numpy](https://pypi.org/project/numpy/)

[PuLP](https://pypi.org/project/PuLP/)

## Documentation

Documentation is found [here](https://vrpy.readthedocs.io/en/latest/).

## Running the tests

### Unit Tests

```sh
python3 -m pytest tests/
```

### Benchmarks

To run some non-regression tests on some benchmarks instances (Solomon and Augerat) do

```sh
python3 -m pytest benchmarks/
```

Note that running the benchmarks requires [pandas](https://pypi.org/project/pandas/) and that it takes a while.

For more information and to run more instances, see the [benchmarks](benchmarks/README.md).

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Kuifje02/vrpy/blob/dev/LICENSE) file for details.

## Bugs

Please report any bugs that you find [here](https://github.com/Kuifje02/vrpy/issues). Or, even better, fork the repository on [GitHub](https://github.com/Kuifje02/vrpy) and create a pull request. Please read the [Community Guidelines](https://github.com/Kuifje02/vrpy/blob/dev/CONTRIBUTING.md) before contributing. Any contributions are welcome.
# Contributing

Contributions are always greatly appreciated and credit will always be given.

## Types of contributions

### Report bugs

Report bugs [here](https://github.com/Kuifje02/vrpy/issues).

If you are reporting a bug, please include:

*   Your operating system name and version.
*   Any details about your local setup that might be helpful in troubleshooting.
*   Detailed steps to reproduce the bug.

### Fix bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to implement it.

### Implement features

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

## Pull request guidelines

Before you submit a pull request, check that it meets these guidelines:

1.  The pull request should include tests.
2.  If the pull request adds functionality, the docs should be updated.
3.  The pull request should work for Python 3.5-3.7. Check the [CircleCI dashboard](https://app.circleci.com/pipelines/github/Kuifje02/vrpy) and make sure that the tests pass for all supported Python versions.
---
title: 'VRPy: A Python package for solving a range of vehicle routing problems with a column generation approach'
tags:
  - Python
  - Vehicle Routing Problems
  - Networks
  - Column generation
authors:
  - name: Romain Montagné
    orcid: 0000-0003-3139-4519
    affiliation: "1"
  - name: David Torres Sanchez
    orcid: 0000-0002-2894-9432
    affiliation: "2"
  - name: Halvard Olsen Storbugt
    orcid: 0000-0003-1142-0185
    affiliation: "2"
affiliations:
 - name: EURODECISION
   index: 1
 - name: SINTEF Digital, Mathematics and Cybernetics
   index: 2
date: June 2020
bibliography: paper.bib
---

# Introduction

The Vehicle Routing Problem (VRP) is amongst the most well known combinatorial optimization problems. The most classical version of the VRP, the Capacitated VRP (CVRP) [@laporte2007you], can be described as follows. A fleet of vehicles with uniform capacity must serve customers with known demand for a single commodity.
The vehicles start and end their routes at a common depot and each customer must be served by exactly one vehicle.
The objective is to assign a sequence of customers to each vehicle of the fleet (a route), minimizing the total distance traveled, such that all customers are served and the total demand served by each vehicle does not exceed its capacity. Note that the VRP generalises the well-known traveling salesman problem (TSP) and is therefore computationally intractable.

Mathematicians have started tackling VRPs since 1959 [@dantzig1959truck]. Ever since, algorithms and computational power have not stopped improving. State of the art techniques include column generation approaches  [@costa2019exact; @bramel1997solving] on which ``vrpy`` relies; more details are given hereafter.

``vrpy`` is of interest to the operational research community and others (e.g., logisticians, supply chain analysts) who wish to solve vehicle routing problems, and therefore has many obvious applications in industry.

# Features

``vrpy`` is a Python package that offers an easy-to-use, unified API for many variants of vehicle routing problems including:

-   the Capacitated VRP (CVRP) [@laporte2007you;@baldacci2010exact],
-   the CVRP with resource constraints [@laporte1985optimal],
-   the CVRP with time windows  [@cordeau2000vrp],
-   the CVRP with simultaneous distribution and collection [@dell2006branch],
-   the CVRP with pickups and deliveries [@desrosiers1988shortest],
-   the CVRP with heterogeneous fleet [@choi2007column].

For each of these variants, it is possible to i/ set initial routes for the search (if one already has a solution at hand and wishes to improve it) ii/ lock routes (if part of the solution is imposed and must not be optimized) iii/ drop nodes (ignore a customer at the cost of a penalty).

``vrpy`` is built upon the well known *NetworkX* library [@hagberg2008exploring] and thus benefits from a user friendly API, as shown in the following quick start example:

```python
from networkx import DiGraph
from vrpy import VehicleRoutingProblem

# Define the network
G = DiGraph()
G.add_edge("Source",1,cost=1,time=2)
G.add_edge("Source",2,cost=2,time=1)
G.add_edge(1,"Sink",cost=0,time=2)
G.add_edge(2,"Sink",cost=2,time=3)
G.add_edge(1,2,cost=1,time=1)
G.add_edge(2,1,cost=1,time=1)

# Define the customers demands
G.nodes[1]["demand"] = 5
G.nodes[2]["demand"] = 4

# Define the Vehicle Routing Problem
prob = VehicleRoutingProblem(G, load_capacity=10, duration=5)

# Solve and display solution value
prob.solve()
print(prob.best_value)
3
print(prob.best_routes)
{1: ["Source",2,1,"Sink"]}
```

# State of the field

Although the VRP is a classical optimization problem, to our knowledge there is only one dedicated package in the Python ecosystem that is able to solve such a range of VRP variants: the excellent ``OR-Tools`` (Google) routing library [@ortools], released for the first time in 2014. To be precise, the core algorithms are implemented in C++, but the library provides a wrapper in Python. Popular and efficient, it is a reference for ``vrpy``, both in terms of features and performance. The current version of ``vrpy`` is able to handle the same variants as OR-Tools (mentioned in the previous section).

Performance-wise, ``vrpy`` ambitions to be competitive with ``OR-Tools`` eventually, at least in terms of solution quality. For the moment, benchmarks (available in the repository) for the CVRP on the set of Augerat instances [@augerat1995approche] show promising results: in the performance profile in Figure 1 below, one can see that nearly the same number of instances are solved within 10 seconds with the same relative error with respect to the best known solution (42\% for ``vrpy``, 44\% for ``OR-Tools``).

| ![Performance profile](cvrp_performance_profile.png) |
| :--------------------------------------------------: |
|         *Figure 1: CVRP Performance profile*         |

We do not claim to outperform ``OR-Tools``, but aim to have results of the same order of magnitude as the package evolves, as there is still much room for improvement (see Section *Future Work* below). On the other hand, we are confident that the user friendly and intuitive API will help students, researchers and more generally the operational research community solve instances of vehicle routing problems of small to medium size, perhaps more easily than with the existing software.

``py-ga-VRPTW`` is another library that is available but as mentioned by its authors, it is more of an experimental project and its performances are rather poor. In particular, we were not able to find feasible solutions for Solomon's instances [@solomon1987algorithms] and therefore cannot compare the two libraries. Also note that ``py-ga-VRPTW`` is designed to solve the VRPTW only, that is, the VRP with time windows.


# Mathematical background

``vrpy`` solves vehicle routing problems with a column generation approach. The term *column generation* refers to the fact that iteratively, routes (or columns) are generated with a pricing problem, and fed to a master problem which selects the best routes among a pool such that each vertex is serviced exactly once. Results from the master problem are then used to search for new potential routes likely to improve the solution's cost, and so forth. This procedure is illustrated in Figure 2 below:

| ![Column Generation](colgen.png) |
| :------------------------------: |
|  *Figure 2: Column Generation*   |

The master problem is a set partitioning linear formulation and is solved with the open source solver Clp from COIN-OR [@johnjforrest_2020_clp], while the subproblem is a shortest elementary path problem with *resource constraints*. It is solved with the help of the  ``cspy`` library [@cspy] which is specifically designed for such problems.

This column generation procedure is very generic, as for each of the featuring VRP variants, the master problem is identical and partitions the customers into subsets (routes). It is the subproblem (or pricing problem) that differs from one variant to another. More specifically, each variant has its unique set of *resources* which must remain in a given interval. For example, for the CVRP, a resource representing the vehicle's load is carried along the path and must not exceed the vehicle capacity; for the CVRP with time windows, two extra resources must be considered: the first one for time, and the second one for time window feasibility. The reader may refer to [@costa2019exact] for more details on each of these variants and how they are delt with within the framework of column generation.

Note that ``vrpy`` does not necessarily return an optimal solution. Indeed, once the pricing problems fails to find
a route with negative marginal cost, the master problem is solved as a MIP. This *price-and-branch* strategy does not guarantee optimality. Note however that it
can be shown [@bramel1997solving] that asymptotically, the relative error goes to zero as the number of customers increases. To guarantee that an optimal solution is returned, the column generation procedure should be embedded in a branch-and-bound scheme (*branch-and-price*), which is beyond the scope of the current release, but part of the future work considered.

# Advanced Features

For more advanced users, there are different pricing strategies (approaches for solving subproblems), namely sparsification strategies [@dell2006branch;@santini2018branch], as well as pre-pricing heuristics available that can lead to faster solutions. The heuristics implemented include a greedy randomized heuristic
(for the CVRP and the CVRP with resource constraints) [@santini2018branch]. Also, a diving heuristic [@sadykov2019primal] can be called to explore part of the branch-and-price tree, instead of solving the restricted master problem as a MIP.

Additionally, we have an experimental feature that uses Hyper-Heuristics for the dynamic selection of
pricing strategies.
The approach ranks the best pricing strategies as the algorithm is running and chooses
according to selection functions based on [@sabar2015math;@ferreira2017multi]. The selection criteria has been modified to include a combination of runtime, objective improvement, and currently active columns in the restricted master problem.
Adaptive parameter settings found in [@drake2012improved] is used to balance exploration and exploitation
under stagnation. The main advantage is that selection is done as the program runs, and is therefore more
flexible compared to a predefined pricing strategy.

# Future Work

There are many ways ``vrpy`` could be improved. To boost the run times, specific heuristics for each variant could be implemented, e.g., Solomon's insertion algorithm [@solomon1987algorithms] for the VRPTW. Second, the pricing problem is solved with ``cspy``, which is quite recent (2019) and is still being fine tuned.  Also, currently, stabilization issues are delt with a basic interior point based strategy which could be enhanced [@pessoa2018automation]. Last but not least, there are many cutting strategies in the literature [@costa2019exact] that have not been implemented and which have proven to significantly reduce run times for such problems.

# Acknowledgements

We would like to thank reviewers Ben Stabler and Serdar Kadioglu for their helpful and constructive suggestions.

# References
# Results

Some summary tables and plots coming soon.

# Replicating results

## Set up

First download the instances you wish to run ([Augerat]() or [Solomon]()) and place them in the
appropriate folders:
 -  Augerat -> `benchmarks/data/cvrp`,
 -  Solomon -> `benchmarks/data/cvrptw`.

For Augerat, ensure that no `.sol` files are left in the folder

## Running

To run the results with the default configuration, from the root folder of the project (`vrpy/`), do

```bash
python3 -m benchmarks.run
```

As it goes, the csv files are created in a new `benchmarks/run/results`.

To see the different options do

```bash
python3 -m benchmarks.run -h
```

These include:
 - Parallel/series runner
 - CPU number specificiation
 - Exploration or performance mode (default configuration)

## OR-TOOLS

To run tests with ortools, [this code](https://github.com/Kuifje02/ortools) can be used.
# Examples

-   [`cvrp.py`](cvrp.py) - Capacitated vehicle routing
-   [`cvrp_drop.py`](cvrp_drop.py) - Capacitated vehicle routing with dropping penalty
-   [`cvrpsdc.py`](cvrpsdc.py) - Capacitated vehicle routing with distribution and collection
-   [`pdp.py`](pdp.py) - Capacitated vehicle routing with pickup and delivery
-   [`vrptw.py`](vrptw.py) - Vehicle routing with time windows

## Run

To run an example, from the main repo folder do

```bash
python3 -m examples.<name>
```

Where `<name>` is the file name (without the extension) above.
For example

```bash
python3 -m examples.cvrp
```
.. _examples:

Examples
========

A simple example
~~~~~~~~~~~~~~~~

Network definition
******************

In this first example, we will be working with the following network:

.. figure:: images/network.png
   :align: center

The first step is to define the network as a ``nx.Digraph`` object. Note that for convenience, the depot (node :math:`0` in the picture) is split into two vertices
: the ``Source`` and the ``Sink``.

.. code:: python

    # Create graph
    >>> from networkx import DiGraph
    >>> G = DiGraph()
    >>> for v in [1, 2, 3, 4, 5]:
           G.add_edge("Source", v, cost=10)
           G.add_edge(v, "Sink", cost=10)
    >>> G.add_edge(1, 2, cost=10)
    >>> G.add_edge(2, 3, cost=10)
    >>> G.add_edge(3, 4, cost=15)
    >>> G.add_edge(4, 5, cost=10)


VRP definition
**************

The second step is to define the VRP, with the above defined graph as input:

.. code:: python

    >>> from vrpy import VehicleRoutingProblem
    >>> prob = VehicleRoutingProblem(G)

Maximum number of stops per route
*********************************

In this first variant, it is required that a vehicle cannot perform more than :math:`3` stops:

.. code:: python

    >>> prob.num_stops = 3
    >>> prob.solve()

The best routes found can be queried as follows:

.. code:: python

    >>> prob.best_routes
    {1: ['Source', 4, 5, 'Sink'], 2: ['Source', 1, 2, 3, 'Sink']}

And the cost of this solution is queried in a similar fashion:

.. code:: python

    >>> prob.best_value
    70.0
    >>> prob.best_routes_cost
    {1: 30, 2: 40}

The optimal routes are displayed below:

.. figure:: images/stops.png
   :align: center

Capacity constraints
********************

In this second variant, we define a demand for each customer and limit the vehicle capacity to :math:`10` units.

Demands are set directly as node attributes on the graph, and the capacity constraint is set with the ``load_capacity`` attribute:

.. code:: python

    >>> for v in G.nodes():
           if v not in ["Source", "Sink"]:
              G.nodes[v]["demand"] = 5
    >>> prob.load_capacity = 10
    >>> prob.solve()
    >>> prob.best_value
    80.0

As the problem is more constrained, it is not surprising that the total
cost increases. As a sanity check, we can query the loads on each route to make sure capacity constraints are met:

.. code:: python

    >>> prob.best_routes
    {1: ["Source", 1, "Sink"], 2: ["Source", 2, 3, "Sink"], 3: ["Source", 4, 5, "Sink"]}
    >>> prob.best_routes_load
    {1: 5, 2: 10, 3: 10}

The new optimal routes are displayed below:

.. figure:: images/capacity.png
   :align: center

Time constraints
****************

One may want to restrict the total duration of a route. In this case, a `time`
attribute is set on each edge of the graph, and a maximum duration is set with `prob.duration`.

.. code:: python

    >>> for (u, v) in G.edges():
           G.edges[u,v]["time"] = 20
    >>> G.edges[4,5]["time"] = 25
    >>> prob.duration = 60
    >>> prob.solve()
    >>> prob.best_value
    85.0

As the problem is more and more constrained, the total cost continues to increase. Lets check the durations of each route:

.. code:: python

    >>> prob.best_routes
    {1: ["Source", 1, 2, "Sink"], 2: ["Source", 3, 4, "Sink"], 3: ["Source", 5, "Sink"]}
    >>> prob.best_routes_duration
    {1: 60, 2: 60, 3: 40}

The new optimal routes are displayed below:

.. figure:: images/time.png
   :align: center

Time window constraints
***********************

When designing routes, it may be required that a customer is serviced in
a given time window :math:`[\ell,u]`. Such time windows are defined for
each node, as well as service times.

.. code-block:: python

    >>> time_windows = {1: (5, 100), 2: (5, 20), 3: (5, 100), 4: (5, 100), 5: (5, 100)}
    >>> for v in G.nodes():
            G.nodes[v]["lower"] = time_windows[v][0]
            G.nodes[v]["upper"] = time_windows[v][1]
            if v not in ["Source","Sink"]:
               G.nodes[v]["service_time"] = 1


A boolean parameter ``time_windows`` is activated to enforce
such constraints:

.. code:: python

    >>> prob.time_windows = True
    >>> prob.duration = 64
    >>> prob.solve()
    >>> prob.best_value
    90.0

The total cost increases again. Lets check the arrival times:

.. code:: python

    >>> prob.best_routes
    {1: ["Source", 1, "Sink"], 4: ["Source", 2, 3, "Sink"], 2: ["Source", 4, "Sink"],  3: ["Source", 5, "Sink"]}
    >>> prob.arrival_time
    {1: {1: 20, 'Sink': 41}, 2: {4: 20, 'Sink': 41}, 3: {5: 20, 'Sink': 41}, 4: {2: 20, 3: 41, 'Sink': 62}}

The new optimal routes are displayed below:

.. figure:: images/time_windows.png
   :align: center

Complete program
****************

.. code:: python

    import networkx as nx
    from vrpy import VehicleRoutingProblem

    # Create graph
    G = nx.DiGraph()
    for v in [1, 2, 3, 4, 5]:
	   G.add_edge("Source", v, cost=10, time=20)
       G.add_edge(v, "Sink", cost=10, time=20)
       G.nodes[v]["demand"] = 5
       G.nodes[v]["upper"] = 100
       G.nodes[v]["lower"] = 5
       G.nodes[v]["service_time"] = 1
    G.nodes[2]["upper"] = 20
    G.nodes["Sink"]["upper"] = 110
    G.nodes["Source"]["upper"] = 100
    G.add_edge(1, 2, cost=10, time=20)
    G.add_edge(2, 3, cost=10, time=20)
    G.add_edge(3, 4, cost=15, time=20)
    G.add_edge(4, 5, cost=10, time=25)

    # Create vrp
    prob = VehicleRoutingProblem(G, num_stops=3, load_capacity=10, duration=64, time_windows=True)

    # Solve and display solution
    prob.solve()
    print(prob.best_routes)
    print(prob.best_value)
	
Periodic CVRP
*************

For scheduling routes over a time period, one can define a frequency for each customer. For example, if over a planning period of two days,
customer :math:`2` must be visited twice, and the other customers only once:

.. code:: python
   
   >>> prob.periodic = 2 
   >>> G.nodes[2]["frequency"] = 2
   >>> prob.solve()
   >>> prob.best_routes
   {1: ['Source', 1, 2, 'Sink'], 2: ['Source', 4, 5, 'Sink'], 3: ['Source', 2, 3, 'Sink']}
   >>> prob.schedule
   {0: [1, 2], 1: [3]}
	
We can see that customer :math:`2` is visited on both days of the planning period (routes :math:`1` and :math:`3`), and that it is not visited more
than once per day.

.. _hfvrp:

Mixed fleet
***********

We end this small example with an illustration of the ``mixed_fleet`` option, when vehicles of different
types (capacities, travel costs, fixed costs) are operating.

The first vehicle has a ``load_capacity`` of :math:`5` units, and no ``fixed_cost``, while
the second vehicle has a ``load_capacity`` of :math:`20` units, and a ``fixed_cost`` with value
:math:`5`. The travel costs of the second vehicle are :math:`1` unit more expensive than 
those of the first vehicle:

.. code:: python

    >>> from networkx import DiGraph
    >>> from vrpy import VehicleRoutingProblem
    >>> G = DiGraph()
    >>> for v in [1, 2, 3, 4, 5]:
           G.add_edge("Source", v, cost=[10, 11])
           G.add_edge(v, "Sink", cost=[10, 11])
           G.nodes[v]["demand"] = 5
    >>> G.add_edge(1, 2, cost=[10, 11])
    >>> G.add_edge(2, 3, cost=[10, 11])
    >>> G.add_edge(3, 4, cost=[15, 16])
    >>> G.add_edge(4, 5, cost=[10, 11])
    >>> prob=VehicleRoutingProblem(G, mixed_fleet=True, fixed_cost=[0, 5], load_capacity=[5, 20])
    >>> prob.best_value
	85
    >>> prob.best_routes
	{1: ['Source', 1, 'Sink'], 2: ['Source', 2, 3, 4, 5, 'Sink']}
    >>> prob.best_routes_cost
	{1: 20, 2: 65}
    >>> prob.best_routes_type
	{1: 0, 2: 1}



An example borrowed from *ortools*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We borrow this second example from the well known *ortools* :cite:`ortools` routing library. We will use the data from the tutorial_.  


Network definition
******************

The graph is considered complete, that is, there are edges between each pair of nodes, in both directions,
and the cost on each edge is defined as the *Manhattan* distance between both endpoints. 
The network is displayed below (for readability, edges are not shown), with the depot in red, and the labels outside of the vertices
are the demands:

.. figure:: images/nodes.png
   :align: center

The network can be entirely defined by its distance matrix.
We will make use of the *NetworkX* module to create this graph and store its attributes:

.. code:: python

 from networkx import DiGraph, from_numpy_matrix, relabel_nodes, set_node_attributes
 from numpy import array

 # Distance matrix
 DISTANCES = [
 [0,548,776,696,582,274,502,194,308,194,536,502,388,354,468,776,662,0], # from Source
 [0,0,684,308,194,502,730,354,696,742,1084,594,480,674,1016,868,1210,548],
 [0,684,0,992,878,502,274,810,468,742,400,1278,1164,1130,788,1552,754,776],
 [0,308,992,0,114,650,878,502,844,890,1232,514,628,822,1164,560,1358,696],
 [0,194,878,114,0,536,764,388,730,776,1118,400,514,708,1050,674,1244,582],
 [0,502,502,650,536,0,228,308,194,240,582,776,662,628,514,1050,708,274],
 [0,730,274,878,764,228,0,536,194,468,354,1004,890,856,514,1278,480,502],
 [0,354,810,502,388,308,536,0,342,388,730,468,354,320,662,742,856,194],
 [0,696,468,844,730,194,194,342,0,274,388,810,696,662,320,1084,514,308],
 [0,742,742,890,776,240,468,388,274,0,342,536,422,388,274,810,468,194],
 [0,1084,400,1232,1118,582,354,730,388,342,0,878,764,730,388,1152,354,536],
 [0,594,1278,514,400,776,1004,468,810,536,878,0,114,308,650,274,844,502],
 [0,480,1164,628,514,662,890,354,696,422,764,114,0,194,536,388,730,388],
 [0,674,1130,822,708,628,856,320,662,388,730,308,194,0,342,422,536,354],
 [0,1016,788,1164,1050,514,514,662,320,274,388,650,536,342,0,764,194,468],
 [0,868,1552,560,674,1050,1278,742,1084,810,1152,274,388,422,764,0,798,776],
 [0,1210,754,1358,1244,708,480,856,514,468,354,844,730,536,194,798,0,662],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # from Sink
 ]

 # Demands (key: node, value: amount)
 DEMAND = {1: 1, 2: 1, 3: 2, 4: 4, 5: 2, 6: 4, 7: 8, 8: 8, 9: 1, 10: 2, 11: 1, 12: 2, 13: 4, 14: 4, 15: 8, 16: 8}

 # The matrix is transformed into a DiGraph
 A = array(DISTANCES, dtype=[("cost", int)])
 G = from_numpy_matrix(A, create_using=nx.DiGraph())

 # The demands are stored as node attributes
 set_node_attributes(G, values=DEMAND, name="demand")

 # The depot is relabeled as Source and Sink
 G = relabel_nodes(G, {0: "Source", 17: "Sink"})

CVRP
****

Once the graph is properly defined, creating a CVRP and solving it is straightforward.
With a maximum load of :math:`15` units per vehicle:

.. code:: python

    >>> from vrpy import VehicleRoutingProblem
    >>> prob = VehicleRoutingProblem(G, load_capacity=15)
    >>> prob.solve()
    >>> prob.best_value
    6208.0
    >>> prob.best_routes
    {1: ['Source', 12, 11, 15, 13, 'Sink'], 2: ['Source', 1, 3, 4, 7, 'Sink'], 3: ['Source', 5, 2, 6, 8, 'Sink'], 4: ['Source', 14, 16, 10, 9, 'Sink']}
    >>> prob.best_routes_load
    {1: 15, 2: 15, 3: 15, 4: 15}


The four routes are displayed below:

.. figure:: images/nodes_capacity.png
   :align: center

CVRP with simultaneous distribution and collection
**************************************************

We follow with the exact same configuration, but this time, every time a node is visited, the vehicle unloads its demand and loads 
some waste material. 

.. code:: python

 from networkx import DiGraph, from_numpy_matrix, relabel_nodes, set_node_attributes
 from numpy import array

 # Distance matrix
 DISTANCES = [
 [0,548,776,696,582,274,502,194,308,194,536,502,388,354,468,776,662,0], # from Source
 [0,0,684,308,194,502,730,354,696,742,1084,594,480,674,1016,868,1210,548],
 [0,684,0,992,878,502,274,810,468,742,400,1278,1164,1130,788,1552,754,776],
 [0,308,992,0,114,650,878,502,844,890,1232,514,628,822,1164,560,1358,696],
 [0,194,878,114,0,536,764,388,730,776,1118,400,514,708,1050,674,1244,582],
 [0,502,502,650,536,0,228,308,194,240,582,776,662,628,514,1050,708,274],
 [0,730,274,878,764,228,0,536,194,468,354,1004,890,856,514,1278,480,502],
 [0,354,810,502,388,308,536,0,342,388,730,468,354,320,662,742,856,194],
 [0,696,468,844,730,194,194,342,0,274,388,810,696,662,320,1084,514,308],
 [0,742,742,890,776,240,468,388,274,0,342,536,422,388,274,810,468,194],
 [0,1084,400,1232,1118,582,354,730,388,342,0,878,764,730,388,1152,354,536],
 [0,594,1278,514,400,776,1004,468,810,536,878,0,114,308,650,274,844,502],
 [0,480,1164,628,514,662,890,354,696,422,764,114,0,194,536,388,730,388],
 [0,674,1130,822,708,628,856,320,662,388,730,308,194,0,342,422,536,354],
 [0,1016,788,1164,1050,514,514,662,320,274,388,650,536,342,0,764,194,468],
 [0,868,1552,560,674,1050,1278,742,1084,810,1152,274,388,422,764,0,798,776],
 [0,1210,754,1358,1244,708,480,856,514,468,354,844,730,536,194,798,0,662],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # from Sink
 ]

 # Delivery demands (key: node, value: amount)
 DEMAND = {1: 1, 2: 1, 3: 2, 4: 4, 5: 2, 6: 4, 7: 8, 8: 8, 9: 1, 10: 2, 11: 1, 12: 2, 13: 4, 14: 4, 15: 8, 16: 8}
 
 # Pickup waste (key: node, value: amount)
 COLLECT = {1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 1, 7: 4, 8: 1, 9: 1, 10: 2, 11: 3, 12: 2, 13: 4, 14: 2, 15: 1, 16: 2}

 # The matrix is transformed into a DiGraph
 A = array(DISTANCES, dtype=[("cost", int)])
 G = from_numpy_matrix(A, create_using=nx.DiGraph())

 # The distribution and collection amounts are stored as node attributes
 set_node_attributes(G, values=DEMAND, name="demand")
 set_node_attributes(G, values=COLLECT, name="collect")

 # The depot is relabeled as Source and Sink
 G = relabel_nodes(G, {0: "Source", 17: "Sink"})


The `load_capacity` is unchanged, and the ``distribution_collection`` attribute is set to ``True.``

.. code:: python

    >>> from vrpy import VehicleRoutingProblem
    >>> prob = VehicleRoutingProblem(G, load_capacity=15, distribution_collection=True)
    >>> prob.solve() 
	>>> prob.best_value
    6208.0
    >>> prob.best_routes
    {1: ['Source', 12, 11, 15, 13, 'Sink'], 2: ['Source', 1, 3, 4, 7, 'Sink'], 3: ['Source', 5, 2, 6, 8, 'Sink'], 4: ['Source', 14, 16, 10, 9, 'Sink']}
    >>> prob.node_load
    {1: {7: 4, 3: 5, 4: 8, 1: 8, 'Sink': 8}, 2: {8: 7, 6: 10, 2: 10, 5: 10, 'Sink': 10}, 3: {14: 2, 16: 8, 10: 8, 9: 8, 'Sink': 8}, 4: {13: 0, 15: 7, 11: 5, 12: 5, 'Sink': 5}}

The optimal solution is unchanged. This is understandable, as for each node, the distribution volume is greater than (or equals) the pickup volume.


VRP with time windows
*********************

Each node must now be serviced within a time window. The time windows are displayed above each node: 

.. figure:: images/nodes_time_windows.png
   :align: center

This time, the network is defined by its distance matrix and its time matrix:

.. code:: python

    from networkx import DiGraph, from_numpy_matrix, relabel_nodes, set_node_attributes
    from numpy import array

    # Distance matrix
    DISTANCES = [
	 [0,548,776,696,582,274,502,194,308,194,536,502,388,354,468,776,662,0], # from Source
	 [0,0,684,308,194,502,730,354,696,742,1084,594,480,674,1016,868,1210,548],
	 [0,684,0,992,878,502,274,810,468,742,400,1278,1164,1130,788,1552,754,776],
	 [0,308,992,0,114,650,878,502,844,890,1232,514,628,822,1164,560,1358,696],
	 [0,194,878,114,0,536,764,388,730,776,1118,400,514,708,1050,674,1244,582],
	 [0,502,502,650,536,0,228,308,194,240,582,776,662,628,514,1050,708,274],
	 [0,730,274,878,764,228,0,536,194,468,354,1004,890,856,514,1278,480,502],
	 [0,354,810,502,388,308,536,0,342,388,730,468,354,320,662,742,856,194],
	 [0,696,468,844,730,194,194,342,0,274,388,810,696,662,320,1084,514,308],
	 [0,742,742,890,776,240,468,388,274,0,342,536,422,388,274,810,468,194],
	 [0,1084,400,1232,1118,582,354,730,388,342,0,878,764,730,388,1152,354,536],
	 [0,594,1278,514,400,776,1004,468,810,536,878,0,114,308,650,274,844,502],
	 [0,480,1164,628,514,662,890,354,696,422,764,114,0,194,536,388,730,388],
	 [0,674,1130,822,708,628,856,320,662,388,730,308,194,0,342,422,536,354],
	 [0,1016,788,1164,1050,514,514,662,320,274,388,650,536,342,0,764,194,468],
	 [0,868,1552,560,674,1050,1278,742,1084,810,1152,274,388,422,764,0,798,776],
	 [0,1210,754,1358,1244,708,480,856,514,468,354,844,730,536,194,798,0,662],
	 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # from Sink
	 ]
	 
    TRAVEL_TIMES = [
	[0, 6, 9, 8, 7, 3, 6, 2, 3, 2, 6, 6, 4, 4, 5, 9, 7, 0],  # from source
	[0, 0, 8, 3, 2, 6, 8, 4, 8, 8, 13, 7, 5, 8, 12, 10, 14, 6],
	[0, 8, 0, 11, 10, 6, 3, 9, 5, 8, 4, 15, 14, 13, 9, 18, 9, 9],
	[0, 3, 11, 0, 1, 7, 10, 6, 10, 10, 14, 6, 7, 9, 14, 6, 16, 8],
	[0, 2, 10, 1, 0, 6, 9, 4, 8, 9, 13, 4, 6, 8, 12, 8, 14, 7],
	[0, 6, 6, 7, 6, 0, 2, 3, 2, 2, 7, 9, 7, 7, 6, 12, 8, 3],
	[0, 8, 3, 10, 9, 2, 0, 6, 2, 5, 4, 12, 10, 10, 6, 15, 5, 6],
	[0, 4, 9, 6, 4, 3, 6, 0, 4, 4, 8, 5, 4, 3, 7, 8, 10, 2],
	[0, 8, 5, 10, 8, 2, 2, 4, 0, 3, 4, 9, 8, 7, 3, 13, 6, 3],
	[0, 8, 8, 10, 9, 2, 5, 4, 3, 0, 4, 6, 5, 4, 3, 9, 5, 2],
	[0, 13, 4, 14, 13, 7, 4, 8, 4, 4, 0, 10, 9, 8, 4, 13, 4, 6],
	[0, 7, 15, 6, 4, 9, 12, 5, 9, 6, 10, 0, 1, 3, 7, 3, 10, 6],
	[0, 5, 14, 7, 6, 7, 10, 4, 8, 5, 9, 1, 0, 2, 6, 4, 8, 4],
	[0, 8, 13, 9, 8, 7, 10, 3, 7, 4, 8, 3, 2, 0, 4, 5, 6, 4],
	[0, 12, 9, 14, 12, 6, 6, 7, 3, 3, 4, 7, 6, 4, 0, 9, 2, 5],
	[0, 10, 18, 6, 8, 12, 15, 8, 13, 9, 13, 3, 4, 5, 9, 0, 9, 9],
	[0, 14, 9, 16, 14, 8, 5, 10, 6, 5, 4, 10, 8, 6, 2, 9, 0, 7],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # from sink
	]
	
    # Time windows (key: node, value: lower/upper bound)
    TIME_WINDOWS_LOWER = {0: 0, 1: 7, 2: 10, 3: 16, 4: 10, 5: 0, 6: 5, 7: 0, 8: 5, 9: 0, 10: 10, 11: 10, 12: 0, 13: 5, 14: 7, 15: 10, 16: 11,}
    TIME_WINDOWS_UPPER = {1: 12, 2: 15, 3: 18, 4: 13, 5: 5, 6: 10, 7: 4, 8: 10, 9: 3, 10: 16, 11: 15, 12: 5, 13: 10, 14: 8, 15: 15, 16: 15,}

    # Transform distance matrix into DiGraph
    A = array(DISTANCES, dtype=[("cost", int)])
    G_d = from_numpy_matrix(A, create_using=DiGraph())

    # Transform time matrix into DiGraph
    A = array(TRAVEL_TIMES, dtype=[("time", int)])
    G_t = from_numpy_matrix(A, create_using=DiGraph())

    # Merge
    G = compose(G_d, G_t)

    # Set time windows
    set_node_attributes(G, values=TIME_WINDOWS_LOWER, name="lower")
    set_node_attributes(G, values=TIME_WINDOWS_UPPER, name="upper")
	
    # The VRP is defined and solved
    prob = VehicleRoutingProblem(G, time_windows=True)
    prob.solve()
	
The solution is displayed below:
	
.. code:: python

    >>> prob.best_value
	6528.0
    >>> prob.best_routes
	{1: ['Source', 9, 14, 16, 'Sink'], 2: ['Source', 12, 13, 15, 11, 'Sink'], 3: ['Source', 5, 8, 6, 2, 10, 'Sink'], 4: ['Source', 7, 1, 4, 3, 'Sink']}
    >>> prob.arrival_time
	{1: {9: 2, 14: 7, 16: 11, 'Sink': 18}, 2: {12: 4, 13: 6, 15: 11, 11: 14, 'Sink': 20}, 3: {5: 3, 8: 5, 6: 7, 2: 10, 10: 14, 'Sink': 20}, 4: {7: 2, 1: 7, 4: 10, 3: 16, 'Sink': 24}}

.. figure:: images/sol.png
   :align: center


CVRP with pickups and deliveries
********************************

In this variant, each demand is made of a pickup node and a delivery node.
Each pickup/delivery pair (or request) must be assigned to the same tour, and within this tour, the pickup node must be 
visited prior to the delivery node (as an item that is yet to be picked up cannot be delivered). 
The total load must not exceed the vehicle's capacity. The requests are displayed below:

.. figure:: images/requests.png
   :align: center

The network is defined as previously, and we add the following data to take into account each request:

.. code:: python

    # Requests (from_node, to_node) : amount
    pickups_deliveries = {(1, 6): 1, (2, 10): 2, (4, 3): 3, (5, 9): 1, (7, 8): 2, (15, 11): 3, (13, 12): 1, (16, 14): 4}
    for (u, v) pickups_deliveries:
        G.nodes[u]["request"] = v
        # Pickups are accounted for positively
        G.nodes[u]["demand"] = pickups_deliveries[(u, v)]
        # Deliveries are accounted for negatively
        G.nodes[v]["demand"] = -pickups_deliveries[(u, v)]

We can now create a pickup and delivery instance with a maximum load of :math:`6` units per vehicle, and with at most :math:`6` stops:

.. code:: python

   >>> from vrpy import VehicleRoutingProblem
   >>> prob = VehicleRoutingProblem(G, load_capacity=6, num_stops=6, pickup_delivery=True)
   >>> prob.solve(cspy=False)
   >>> prob.best_value
   5980.0
   >>> prob.best_routes
   {1: ['Source', 5, 2, 10, 16, 14, 9, 'Sink'], 2: ['Source', 7, 4, 3, 1, 6, 8, 'Sink'], 3: ['Source', 13, 15, 11, 12, 'Sink']}
   >>> prob.node_load
   {1: {5: 1, 2: 3, 10: 1, 16: 5, 14: 1, 9: 0, 'Sink': 0}, 2: {7: 2, 4: 5, 3: 2, 1: 3, 6: 2, 8: 0, 'Sink': 0}, 3: {13: 1, 15: 4, 11: 1, 12: 0, 'Sink': 0}}

The four routes are displayed below:

.. figure:: images/pdp.png
   :align: center



Limited fleet and dropping visits
*********************************

This last example is similar to the above CVRP, except for the fact that demands have increased, and that the fleet is limited to :math:`4` vehicles,
with a :math:`15` unit capacity (per vehicle). Since the total demand is greater than :math:`4 \times 15 = 60`, servicing each node is not possible, therefore, we will try to visit
as many customers as possible, and allow dropping visits, at the cost of a :math:`1000` penalty.

.. code:: python

 from networkx import DiGraph, from_numpy_matrix, relabel_nodes, set_node_attributes
 from numpy import array

 # Distance matrix
 DISTANCES = [
 [0,548,776,696,582,274,502,194,308,194,536,502,388,354,468,776,662,0], # from Source
 [0,0,684,308,194,502,730,354,696,742,1084,594,480,674,1016,868,1210,548],
 [0,684,0,992,878,502,274,810,468,742,400,1278,1164,1130,788,1552,754,776],
 [0,308,992,0,114,650,878,502,844,890,1232,514,628,822,1164,560,1358,696],
 [0,194,878,114,0,536,764,388,730,776,1118,400,514,708,1050,674,1244,582],
 [0,502,502,650,536,0,228,308,194,240,582,776,662,628,514,1050,708,274],
 [0,730,274,878,764,228,0,536,194,468,354,1004,890,856,514,1278,480,502],
 [0,354,810,502,388,308,536,0,342,388,730,468,354,320,662,742,856,194],
 [0,696,468,844,730,194,194,342,0,274,388,810,696,662,320,1084,514,308],
 [0,742,742,890,776,240,468,388,274,0,342,536,422,388,274,810,468,194],
 [0,1084,400,1232,1118,582,354,730,388,342,0,878,764,730,388,1152,354,536],
 [0,594,1278,514,400,776,1004,468,810,536,878,0,114,308,650,274,844,502],
 [0,480,1164,628,514,662,890,354,696,422,764,114,0,194,536,388,730,388],
 [0,674,1130,822,708,628,856,320,662,388,730,308,194,0,342,422,536,354],
 [0,1016,788,1164,1050,514,514,662,320,274,388,650,536,342,0,764,194,468],
 [0,868,1552,560,674,1050,1278,742,1084,810,1152,274,388,422,764,0,798,776],
 [0,1210,754,1358,1244,708,480,856,514,468,354,844,730,536,194,798,0,662],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # from Sink
 ]

 # Demands (key: node, value: amount)
 DEMAND = {1: 1, 2: 1, 3: 3, 4: 6, 5: 3, 6: 6, 7: 8, 8: 8, 9: 1, 10: 2, 11: 1, 12: 2, 13: 6, 14: 6, 15: 8, 16: 8}

 # The matrix is transformed into a DiGraph
 A = array(DISTANCES, dtype=[("cost", int)])
 G = from_numpy_matrix(A, create_using=nx.DiGraph())

 # The demands are stored as node attributes
 set_node_attributes(G, values=DEMAND, name="demand")

 # The depot is relabeled as Source and Sink
 G = relabel_nodes(G, {0: "Source", 17: "Sink"})


Once the graph is properly defined, a VRP instance is created, with attributes `num_vehicles` and `drop_penalty`:

.. code:: python

    >>> from vrpy import VehicleRoutingProblem
    >>> prob = VehicleRoutingProblem(G, load_capacity=15, num_vehicles=4, drop_penalty=1000)
    >>> prob.solve()
    >>> prob.best_value
    7776.0
    >>> prob.best_routes
    {1: ['Source', 9, 10, 2, 6, 5, 'Sink'], 2: ['Source', 7, 13, 'Sink'], 3: ['Source', 14, 16, 'Sink'], 4: ['Source', 1, 4, 3, 11, 12, 'Sink']}
    >>> prob.best_routes_load
    {1: 13, 2: 14, 3: 14, 4: 13}

The solver drops nodes :math:`8` and :math:`15`. The new optimal routes are displayed below:

.. figure:: images/drop.png
   :align: center

.. _tutorial: https://developers.google.com/optimization/routing/vrp

Bibliography
------------

.. bibliography:: refs.bib
   :all:
.. _benchmarks:

Performance profiles
====================

Performance profiles are a practical way to have a global overview of a set of algorithms' performances. 
On the :math:`x` axis, we have the relative gap (%), and on the :math:`y` axis, the percentage of data sets solved within the gap. 
So for example, at the intersection with the :math:`y` axis is the percentage of data sets solved optimally, 
and at the intersection with  :math:`y=100\%`  is the relative gap within which all data sets are solved.

At a glance, the more the curve is in the upper left corner, the better the algorithm.

We compare the performances of `vrpy` and `OR-Tools` (default options):

-   on Augerat_'s instances (CVRP),
-   on Solomon_'s instances (CVRPTW)

Results are found here_ [link to repo] and can be replicated.

CVRP
----

.. figure:: images/cvrp_performance_profile.png
   :align: center
   
We can see that with a maximum running time of :math:`10` seconds, `OR-Tools` solves :math:`15\%` of the instances optimally, 
while `vrpy` only solves :math:`5\%` of them. Both solve approximately :math:`43\%` of instances with a maximum relative gap of :math:`5\%`.
And both solve all instances within a maximum gap of :math:`25\%`.

CVRPTW
------

Coming soon.

.. _Augerat: https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-instances/
.. _Solomon: https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-with-time-windows-instances/
.. _here: Using `VRPy`
============

In order to use the VRPy package, first, one has to create a directed graph which represents the underlying network.

To do so, we make use of the well-known `NetworkX` package, with the following input requirements:

 - Input graphs must be of type :class:`networkx.DiGraph`;
 - Input graphs must have a single `Source` and `Sink` nodes with no incoming or outgoing edges respectively;
 - There must be at least one path from `Source` to `Sink`;
 - Edges in the input graph must have a ``cost`` attribute (of type :class:`float`).


For example the following simple network fulfills the requirements listed above:

.. code-block:: python

	>>> from networkx import DiGraph
	>>> G = DiGraph()
	>>> G.add_edge("Source", 1, cost=1)
	>>> G.add_edge("Source", 2, cost=2)
	>>> G.add_edge(1, "Sink", cost=0)
	>>> G.add_edge(2, "Sink", cost=2)
	>>> G.add_edge(1, 2, cost=1)
	>>> G.add_edge(2, 1, cost=1)
	
The customer demands are set as ``demand`` attributes (of type :class:`float`) on each node:

.. code-block:: python

	>>> G.nodes[1]["demand"] = 5
	>>> G.nodes[2]["demand"] = 4
		
To solve your routing problem, create a :class:`VehicleRoutingProblem` instance, specify the problem constraints (e.g., the ``load_capacity`` of each truck), and call ``solve``.

.. code-block:: python

    >>> from vrpy import VehicleRoutingProblem
    >>> prob = VehicleRoutingProblem(G, load_capacity=10)
    >>> prob.solve()

Once the problem is solved, we can query useful attributes as:

.. code-block:: python

    >>> prob.best_value
    3
    >>> prob.best_routes
    {1: ["Source", 2, 1, "Sink"]}
    >>> prob.best_routes_load
    {1: 9}

``prob.best_value`` is the overall cost of the solution, ``prob.best_routes`` is a `dict` object where keys represent the route ID, while the values are
the corresponding path from `Source` to `Sink`. And ``prob.best_routes_load`` is a `dict` object where the same keys point to the accumulated load on the
vehicle.


Different options and constraints are detailed in the :ref:`vrp` section, 
and other attributes can be queried depending on the nature of the VRP (see section :ref:`api`).


.. _colgen:

Mathematical Background
=======================


A column generation approach
----------------------------

*VRPy* solves vehicle routing problems with a column generation approach. The term `column generation` refers to the fact 
that iteratively, routes (or `columns`) are `generated` with a pricing problem, and fed to a master problem which selects the best routes among
a pool such that each vertex is serviced exactly once. The linear formulations of these problems are detailed hereafter.  
	
Master Problem
**************
Let :math:`G=(V,A)` be a graph where :math:`V` denotes the set of nodes that have to be visited, and :math:`A` the set of edges of the network. 
Let :math:`\Omega` be the set of feasible routes. 
Let :math:`\lambda_r` be a binary variable that takes value :math:`1` if and only if route :math:`r \in \Omega` with cost :math:`c_r` is selected. 
The master problem reads as follows:


.. math:: 

	\min \; \sum_{r \in \Omega} c_r \lambda_r

subject to set covering constraints:

.. math:: 

	\sum_{r \in \Omega \mid v \in r} \lambda_r &= 1 \quad &\forall v \in V\quad &(1)

	\lambda_r &\in \{ 0,1\} \quad &\forall r \in \Omega \quad &(2)

   

When using a column generation procedure, integrity constraints :math:`(2)` are relaxed (such that :math:`0 \le \lambda_r \le 1`), and only a subset of :math:`\Omega` is used. 
This subset is generated dynamically with the following sub problem.


Pricing problem
***************

Let :math:`\pi_v` denote the dual variable associated with constraints :math:`(1)`. The marginal cost of a variable (or column) :math:`\lambda_r` is given by:

.. math:: 

	\hat{c}_r = c_r - \sum_{v \in V\mid v \in r} \pi_v

Therefore, if :math:`x_{uv}` is a binary variable that takes value :math:`1` if and only if edge :math:`(u,v)` is used, 
*assuming there are no negative cost sub cycles*, one can formulate the problem of finding a route with negative marginal cost as follows :
 
.. math:: 

	\min \quad   \sum_{(u,v)\in A}c_{uv}x_{uv} -\sum_{u\mid (u,v) \in A}\pi_u x_{uv}

subject to flow balance constraints :

.. math::  

    \sum_{u\mid (u,v) \in A} x_{uv} &=  \sum_{u\mid (v,u) \in A} x_{uv}\quad &\forall v \in V \label{eq3}
	
    x_{uv} &\in \{ 0,1\} \quad &\forall (u,v) \in A \label{eq4}


In other words, the sub problem is a shortest elementary path problem, and additional constraints (such as capacities, time) 
give rise to a shortest path problem with *resource constraints*, hence the interest of using the *cspy* library.

If there are negative cost cycles (which typically happens), the above formulation requires additional constraints
to enforce path elementarity, and the problem becomes computationally intractable.
Linear formulations are then impractical, and algorithms such as the ones available in *cspy* become very handy.


Does VRPy return an optimal solution?
-------------------------------------

*VRPy* does not necessarily return an optimal solution (even with no time limit). Indeed, once the pricing problems fails to find
a route with negative marginal cost, the master problem is solved as a MIP. This *price-and-branch* strategy does not guarantee optimality. Note however that it
can be shown :cite:`bramel1997solving` that asymptotically, the relative error goes to zero as the number of customers increases.   
To guarantee that an optimal solution is returned, the column generation procedure should be embedded in a branch-and-bound scheme (*branch-and-price*). This
is part of the future work listed below.

TO DO
-----

- Embed the solving procedure in a branch-and-bound scheme:

  - branch-and-price (exact)
  - diving (heuristic)
- Implement heuristics for initial solutions.
- More acceleration strategies:

  - other heuristic pricing strategies
  - switch to other LP modeling library (?)
  - improve stabilization
  - ...
- Include more VRP variants:

  - pickup and delivery with cspy
  - ...


.. _vrp:

Vehicle Routing Problems
========================

The `VRPy` package can solve the following VRP variants.


Capacitated Vehicle Routing Problem (CVRP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the capacitated vehicle routing problem (CVRP), a fleet of vehicles with uniform capacity must serve customers with known demand for a single commodity.
The vehicles start and end their routes at a common depot and each customer must be served by exactly one vehicle.
The objective is to assign a sequence of customers (a route) to each truck of the fleet, minimizing the total distance traveled, 
such that all customers are served and the total demand served by each truck does not exceed its capacity. 

.. code-block:: python

	>>> from networkx import DiGraph
	>>> from vrpy import VehicleRoutingProblem
	>>> G = DiGraph()
	>>> G.add_edge("Source", 1, cost=1)
	>>> G.add_edge("Source", 2, cost=2)
	>>> G.add_edge(1, "Sink", cost=0)
	>>> G.add_edge(2, "Sink", cost=2)
	>>> G.add_edge(1, 2, cost=1)
	>>> G.nodes[1]["demand"] = 2
	>>> G.nodes[2]["demand"] = 3
	>>> prob = VehicleRoutingProblem(G, load_capacity=10)
	>>> prob.solve()
	
Note that whether the problem is a distribution or a collection problem does not matter. Both are modeled identically.

	
CVRP with resource constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
Other resources can also be considered:

	- maximum duration per trip; 
	- maximum number of customers per trip.  

Taking into account duration constraints requires setting ``time`` attributes on each edge, and setting
the ``duration`` attribute to the maximum amount of time per vehicle.

Following the above example:

.. code-block:: python

	>>> G.edges["Source",1]["time"] = 5
	>>> G.edges["Source",2]["time"] = 4
	>>> G.edges[1,2]["time"] = 2
	>>> G.edges[1,"Sink"]["time"] = 6
	>>> G.edges[2,"Sink"]["time"] = 1
	>>> prob.duration = 9
	>>> prob.solve()
	
Similarly, imposing a maximum number of customers per trip is done by setting the ``num_stops`` attribute to the desired value.

.. code-block:: python

	>>> prob.num_stops = 1
	

CVRP with Time Windows (CVRPTW)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this variant, deliveries must take place during a given time-window, which can be different for each customer.

Such constraints can be taken into account by setting ``lower`` and ``upper`` attributes on each node, and by activating the
``time_windows`` attribute to ``True.`` Additionally, service times can be taken into account on each node by setting the ``service_time``
attribute.

Following the above example:

.. code-block:: python

	>>> G.nodes[1]["lower"] = 0
	>>> G.nodes[1]["upper"] = 10
	>>> G.nodes[2]["lower"] = 5
	>>> G.nodes[2]["upper"] = 9
	>>> G.nodes[1]["service_time"] = 1
	>>> G.nodes[2]["service_time"] = 2
	>>> prob.time_windows = True
	>>> prob.solve()
	
.. note:: 

	Waiting time is allowed upon arrival at a node. This means that if a vehicle arrives at a node before the time window's
	lower bound, the configuration remains feasible, it is considered that the driver waits before servicing the customer. 
        


CVRP with Simultaneous Distribution and Collection (CVRPSDC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this variant, when a customer is visited, two operations are done simultaneously. Some good is delivered, and some waste material is picked-up. 
The total load must not exceed the vehicle's capacity.

The amount that is picked-up is set with the ``collect`` attribute on each node, and the ``distribution_collection`` attribute is set to ``True.``

Following the above example:

.. code-block:: python

	>>> G.nodes[1]["collect"] = 2
	>>> G.nodes[2]["collect"] = 1
	>>> prob.load_capacity = 2
	>>> prob.distribution_collection = True
	>>> prob.solve()
	
CVRP with Pickup and Deliveries (VRPPD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the pickup-and-delivery problem, each demand is made of a pickup node and a delivery node.
Each pickup/delivery pair (or request) must be assigned to the same tour, and within this tour, the pickup node must be 
visited prior to the delivery node (as an item that is yet to be picked up cannot be delivered). 
The total load must not exceed the vehicle's capacity.

For every pickup node, the ``request`` attribute points to the name of the delivery node. Also, the ``pickup_delivery`` attribute
is set to ``True``. The amount of goods to be shipped is counted positively for the pickup node, and negatively for the delivery node.
For example, if :math:`3` units must be shipped from node :math:`1` to node :math:`2`, the ``demand`` attribute is set to :math:`3` for node :math:`1`, and :math:`-3` for node :math:`2`.

.. code-block:: python

	>>> G.nodes[1]["request"] = 2
	>>> G.nodes[1]["demand"] = 3
	>>> G.nodes[2]["demand"] = -3
	>>> prob.pickup_delivery = True
	>>> prob.load_capacity = 10
	>>> prob.solve(cspy=False)

.. note:: This variant has to be solved with the ``cspy`` attribute set to False. 

Periodic CVRP (PCVRP)
~~~~~~~~~~~~~~~~~~~~~

In the periodic CVRP, the planning period is extended over a time horizon, and customers can be serviced more than once. 
The demand is considered constant over time, and the frequencies (the number of visits) of each customer are known. 

For each node, the ``frequency`` attribute (type :class:`int`) is set, and the parameter ``periodic`` is set to the value of the considered time span (the planning period).
All nodes that have no frequency are visited exactly once. 

.. code-block:: python

	>>> G.nodes[1]["frequency"] = 2
	>>> prob.periodic = 2
	>>> prob.solve()
	
A planning of the routes can then be queried by calling ``prob.schedule,`` which returns a dictionary with keys day numbers and values the list of
route numbers scheduled this day. 
	
.. note:: 

	The PCVRP usually has additional constraints: some customers can only be serviced on specific days of the considered time span. 
	For example, over a :math:`3` day planning period, a node with frequency :math:`2` could only be visited on days :math:`1` and
	:math:`2` or :math:`2` and :math:`3` but not :math:`1` and :math:`3`. Such *combination* constraints are not taken into account by 
	*VRPy* (yet).
	
.. note::

   If the parameter ``num_vehicles`` is used, it refers to the maximum number of vehicles available per day (and not over the time span).
	
CVRP with heterogeneous fleet (HFCVRP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the CVRP with *heterogeneous fleet* (or mixed fleet), there are different types of vehicles, which can differ in capacities and costs (fixed costs
and travel costs). Typically, a vehicle with a larger capacity will be more expensive. The problem consists in finding the best combination of
vehicles to satisfy the demands while minimizing global costs. 

First, the ``cost`` attribute on each of the graph is now a *list* of costs, with as many items as vehicle types (even if costs are equal). For example,
if there are *two* types of vehicles, the following graph satisfies the input requirements:

.. code-block:: python

	>>> from networkx import DiGraph
	>>> G = DiGraph()
	>>> G.add_edge("Source", 1, cost=[1, 2])
	>>> G.add_edge("Source", 2, cost=[2, 4])
	>>> G.add_edge(1, "Sink", cost=[0, 0])
	>>> G.add_edge(2, "Sink", cost=[2, 4])
	>>> G.add_edge(1, 2, cost=[1, 2])

When defining the ``VehicleRoutingProblem``, the ``mixed_fleet`` argument is set to ``True``, and the ``load_capacity`` argument is now also of type :class:`list`,
where each item of the list is the maximum load per vehicle type. For example, if the two types of vehicles have capacities :math:`10` and :math:`15`, respectively:

.. code-block:: python

    >>> from vrpy import VehicleRoutingProblem
    >>> prob = VehicleRoutingProblem(G, mixed_fleet=True, load_capacity=[10, 15])

Note how the dimensions of ``load_capacity`` and ``cost`` are consistent: each list must have as many items as vehicle types, and the
order of the items of the ``load_capacity`` list is consistent with the order of the ``cost`` list on every edge of the graph.
  
Once the problem is solved, the type of vehicle per route can be queried with ``prob.best_routes_type``.  

	
VRP options
~~~~~~~~~~~

In this subsection are described different options which arise frequently in vehicle routing problems.

Open VRP
^^^^^^^^

The `open` VRP refers to the case where vehicles can start and/or end their trip anywhere, instead of having to leave from
the depot, and to return there after service. This is straightforward to model : setting distances (or costs) to :math:`0` on every edge outgoing from the Source 
and incoming to the Sink achieves this.

Fixed costs
^^^^^^^^^^^

Vehicles typically have a *fixed cost* which is charged no matter what the traveled distance is. This can be taken into account with the ``fixed_cost`` attribute.
For example, if the cost of using each vehicle is :math:`100`: 

.. code-block:: python

	>>> prob.fixed_cost = 100
	
.. note:: 

	If the fleet is mixed, the same logic holds for ``fixed_cost``: a list of costs is given, where each item of the list is the fixed cost per vehicle type.
	The order of the items of the list has to be consistent with the other lists (``cost`` and ``load_capacity``). See example :ref:`hfvrp`.
	
Limited fleet
^^^^^^^^^^^^^
	
It is possible to limit the size of the fleet. For example, if at most :math:`10` vehicles are available:

.. code-block:: python

	>>> prob.num_vehicles = 10
	
.. note:: 

	If the fleet is mixed, the same logic holds for ``num_vehicles``: a list of integers is given, where each item of the list is
	the number of available vehicles, per vehicle type. The order of the items of the list has to be consistent with the other
	lists (``cost``, ``load_capacity``, ``fixed_cost``).
	
And to enforce exactly ``num_vehicles`` vehicles, one can use the ``use_all_vehicles``:

.. code-block:: python

	>>> prob.use_all_vehicles = True
	
Dropping visits
^^^^^^^^^^^^^^^

Having a limited fleet may result in an infeasible problem. For example, if the total demand at all locations exceeds the total capacity of the fleet,
the problem has no feasible solution. It may then be interesting to decide which visits to drop in order to meet capacity constraints
while servicing as many customers as possible. To do so, we set the ``drop_penalty`` attribute to an integer value that the solver
will add to the total travel cost each time a node is dropped. For example, if the value of the penalty is :math:`1000`:

.. code-block:: python

	>>> prob.drop_penalty = 1000
	
This problem is sometimes referred to as the `capacitated profitable tour problem` or the `prize collecting tour problem.`

Minimizing the global time span
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to modify the objective function in order to solve a min-max problem. More specifically, the total time span can be minimized
by setting the ``minimize_global_span`` to ``True``. Of course this assumes edges have a ``time`` argument:

.. code-block:: python

	>>> prob.minimize_global_span = True
	
.. note::

   This may lead to poor computation times.
	
Other VRPs
~~~~~~~~~~

Coming soon:

- CVRP with multiple depots

Getting started
===============

Installation
************

You can install the latest release of VRPy from PyPi_ by:

.. code-block:: none

    pip install vrpy

.. _PyPi: https://pypi.python.org/pypi/vrpy

Requirements
************
The requirements for running VRPy are:

 - cspy_: Constrained shortest path problem algorithms :cite:`cspy`.
 - NetworkX_: Graph manipulation and creation :cite:`hagberg2008exploring`.
 - numpy_: Array manipulation :cite:`numpy`.
 - PuLP_: Linear programming modeler.

.. _cspy: https://pypi.org/project/cspy/
.. _NetworkX: https://networkx.github.io/documentation/stable/
.. _numpy: https://pypi.org/project/numpy/
.. _PuLP: https://pypi.org/project/PuLP/.. _options:

Solving Options
===============

Setting initial routes for a search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, an initial solution is computed with the well known Clarke and Wright algorithm :cite:`clarke1964scheduling`. If one already has a feasible solution at hand,
it is possible to use it as an initial solution for the search of a potential better configuration. The solution is passed to the solver as a list of routes, where a route is a list
of nodes starting from the *Source* and ending at the *Sink*. 

.. code-block:: python

	>>> prob.solve(initial_solution = [["Source",1,"Sink"],["Source",2,"Sink"]])
	
Returning solution from initial heuristic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to return the solution found by the Clarke and Wright algorithm by setting the ``heuristic_only`` argument to *True*.
	
.. code-block:: python

	>>> prob.solve(heuristic_only=True)
	
Note that this only possible with capacity and/or resource constraints.
	
Locking routes
~~~~~~~~~~~~~~

It is possible to constrain the problem with partial routes if preassignments are known. There are two possibilites : either a complete route is known, 
and it should not be optimized, either only a partial route is known, and it may be extended. Such routes are given to the solver
with the ``preassignments`` argument. A route with `Source` and `Sink` nodes is considered complete and is locked. Otherwise, the solver will extend it if it yields savings.

In the following example, one route must start with customer :math:`1`, one route must contain edge :math:`(4,5)`, and one complete route,
`Source-2-3-Sink`, is locked.

.. code-block:: python

	>>> prob.solve(preassignments = [["Source",1],[4,5],["Source",2,3,"Sink"]])


Setting a time limit
~~~~~~~~~~~~~~~~~~~~

The ``time_limit`` argument can be used to set a time limit, in seconds. 
The solver will return the best solution found after the time limit has elapsed.

For example, for a one minute time limit:

.. code-block:: python

	>>> prob.solve(time_limit=60)


Linear programming or dynamic programming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`VRPy`'s ``solve`` method relies on a column generation procedure. At every iteration, a master problem and a sub problem are solved.
The sub problem consists in finding variables which are likely to improve the master problem's objective function. 
See section :ref:`colgen` for more details.

The sub problem - or pricing problem - can be solved either with linear programming, or with dynamic programming. Switching to linear 
programming can be done by deactivating the ``cspy`` argument when calling the ``solve`` method. 
In this case the CBC_ :cite:`forrest2018coin` solver of COIN-OR is used by default. 

.. code-block:: python

	>>> prob.solve(cspy=False)
	
The sub problems that are solved are typically computationally intractable, and using dynamic programming is typically quicker, as such algorithms run in pseudo-polynomial time.
However, solving the sub problems as MIPs may also be effective depending on the data set. Also, using commercial solvers may significantly help accelerating the procedure.
If one has CPLEX or GUROBI at hand, they can be used by setting the ``solver`` parameter to "cplex" or "gurobi".

.. code-block:: python

	>>> prob.solve(cspy=False, solver="gurobi")

.. _CBC : https://github.com/coin-or/Cbc
	
Pricing strategy
~~~~~~~~~~~~~~~~

In theory, at each iteration, the sub problem is solved optimally. VRPy does so with a bidirectional labeling algorithm with dynamic halfway point :cite:`tilk2017asymmetry` from the `cspy` library.

This may result in a slow convergence. To speed up the resolution, there are two ways to change this pricing strategy: 

1. By deactivating the ``exact`` argument of the ``solve`` method, `cspy` calls one of its heuristics instead of the bidirectional search algorithm. The exact method is run only once the heuristic fails to find a column with negative reduced cost.

.. code-block:: python

	>>> prob.solve(exact=False)
	
 
2. By modifying the ``pricing_strategy`` argument of the ``solve`` method to one of the following:

	- `BestEdges1`,
	- `BestEdges2`,
	- `BestPaths`,
	- `Hyper`
	

.. code-block:: python

	>>> prob.solve(pricing_strategy="BestEdges1")
	
`BestEdges1`, described for example in :cite:`dell2006branch`, is a sparsification strategy: a subset of nodes and
edges are removed to limit the search space. The subgraph is created as follows: all edges :math:`(i,j)` which verify :math:`c_{ij} > \alpha \; \pi_{max}` are discarded, where :math:`c_{ij}` is the edge's cost, :math:`\alpha \in ]0,1[` is parameter,
and :math:`\pi_{max}` is the largest dual value returned by the current restricted relaxed master problem. The parameter :math:`\alpha` is increased iteratively until
a route is found. `BestEdges2` is another sparsification strategy, described for example in :cite:`santini2018branch`. The :math:`\beta` edges with highest reduced cost are discarded, where :math:`\beta` is a parameter that is increased iteratively.
As for `BestPaths`, the idea is to look for routes in the subgraph induced by the :math:`k` shortest paths from the Source to the Sink (without any resource constraints),
where :math:`k` is a parameter that is increased iteratively.

Additionally, we have an experimental feature that uses Hyper-Heuristics for the dynamic selection of pricing strategies. 
The approach ranks the best pricing strategies as the algorithm is running and chooses according to selection functions based on :cite:`sabar2015math,ferreira2017multi`. 
The selection criteria has been modified to include a combination of runtime, objective improvement, and currently active columns in the restricted master. Adaptive parameter settings found in :cite:`drake2012improved` is used to balance exploration and exploitation under stagnation. The main advantage is that selection is done as the programme runs, and is therefore more flexible compared to a predefined pricing strategy.

For each of these heuristic pricing strategies, if a route with negative reduced cost is found, it is fed to the master problem. Otherwise,
the sub problem is solved exactly. 

The default pricing strategy is `BestEdges1`, with ``exact=True`` (i.e., with the bidirectional labeling algorithm).

A greedy randomized heuristic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the CVRP, or the CVRP with resource constraints, one can activate the option of running a greedy randomized heuristic before pricing:

.. code-block:: python

	>>> prob.solve(greedy="True")

This algorithm, described in :cite:`santini2018branch`, generates a path starting at the *Source* node and then randomly selects an edge among the :math:`\gamma` outgoing edges
of least reduced cost that do not close a cycle and that meet operational constraints (:math:`\gamma` is a parameter).
This is repeated until the *Sink* node is reached . The same procedure is applied backwards, starting from the *Sink* and ending at the *Source*, and is run
:math:`20` times. All paths with negative reduced cost are added to the pool of columns.
VRPy Documentation
====================

VRPy is a python framework for solving instances of different types of Vehicle Routing Problems (VRP) including:

-   the Capacitated VRP (CVRP),
-   the CVRP with resource constraints,
-   the CVRP with time windows (CVRPTW),
-   the CVRP with simultaneous distribution and collection (CVRPSDC),
-   the CVRP with heterogeneous fleet (HFCVRP).

Check out section :ref:`vrp` to find more variants and options.

VRPy relies on the well known NetworkX_ package (graph manipulation), as well as on cspy_, a library for solving the resource constrained shortest path problem.

.. _NetworkX: Graph manipulation and creation.
.. _cspy: https://pypi.org/project/cspy/
   
Disclaimer
==========

There is no guarantee that VRPy returns the optimal solution. See section :ref:`colgen` for more details, and section :ref:`benchmarks`
for performance comparisons with OR-Tools_. 

.. _OR-Tools: https://developers.google.com/optimization/routing/vrp

Authors
=======

Romain Montagné (r.montagne@hotmail.fr)

David Torres Sanchez (d.torressanchez@lancs.ac.uk)

Contributors
============

@Halvaros

Table of contents
=================

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   getting_started
   how_to
   vrp_variants
   solving_options
   examples
   api
   mathematical_background
   benchmarks
   bibliography

* :ref:`genindex`
* :ref:`search`
.. _api:

API
===

vrpy.VehicleRoutingProblem
--------------------------

.. automodule:: vrpy.vrp
   :members:
   :inherited-members:
   

Notes
-----

The input graph must have single `Source` and `Sink` nodes with no incoming or outgoing edges respectively. 
These dummy nodes represent the depot which is split for modeling convenience. The `Source` and `Sink` cannot have a demand, if 
one is given it is ignored with a warning.

Please read sections :ref:`vrp`, :ref:`options` and :ref:`examples` for details and examples on each of the above arguments.

