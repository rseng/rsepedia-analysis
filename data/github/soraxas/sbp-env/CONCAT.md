# Sampling-Based Motion Planners' Testing Environment

[![Python version](https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue.svg)](https://cs.tinyiu.com/sbp-env)
[![CI](https://github.com/soraxas/sbp-env/actions/workflows/ci.yaml/badge.svg)](https://github.com/soraxas/sbp-env/actions/workflows/ci.yaml)
[![Build docs](https://github.com/soraxas/sbp-env/actions/workflows/sphinx.yaml/badge.svg)](https://cs.tinyiu.com/sbp-env)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License](https://img.shields.io/github/license/soraxas/sbp-env.svg)](https://github.com/soraxas/sbp-env/blob/master/LICENSE)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03782/status.svg)](https://doi.org/10.21105/joss.03782)

Sampling-based motion planners' testing environment (`sbp-env`) is a full feature framework to quickly test different sampling-based algorithms for motion planning. `sbp-env` focuses on the flexibility of tinkering with different aspects of the framework, and had divided the main planning components into two categories (i) **samplers** and (ii) **planners**.

The focus of *motion planning research* had been mainly on (i) improving the sampling efficiency (with methods such as heuristic or learned distribution) and (ii) the algorithmic aspect of the planner using different routines to build a connected graph. Therefore, by separating the two components one can quickly swap out different components to test novel ideas.

Have a look at the [documentations](https://cs.tinyiu.com/sbp-env) for more detail information. If you are looking for the previous code for the RRdT* paper it is now archived at [soraxas/rrdt](https://github.com/soraxas/rrdt).

## Installation

#### Optional

I recommend first creates a virtual environment with

```sh
# assumes python3 and bash shell
python -m venv sbp_env
source sbp_env/bin/activate
```

#### Install dependencies

You can install all the needed packages with pip.

```sh
pip install -r requirements.txt
```

There is also an optional dependency on [`klampt`](https://github.com/krishauser/Klampt) if you want to use the 3D simulator. Refer to its [installation guide](https://github.com/krishauser/Klampt#installation) for details.

<img align="right" width="300" height="auto" src="docs/images/klampt-simulator.png" />

## Quick Guide

You can get a detailed help message with

```sh
python main.py --help
```

but the basic syntax is

```sh
python main.py <PLANNER> <MAP> [options]
```

It will open a new window that display a map on it. Every white pixel is assumed to be free, and non-white pixels are obstacles. You will need to use your mouse to select two points on the map, the first will be set as the starting point and the second as the goal point.

## Demos

### Run maps with different available Planners

This repository contains a framework to performs quick experiments for Sampling-Based Planners (SBPs) that are implemented in Python. The followings are planners that had implemented and experimented in this framework.

Note that the commands shown in the respective demos can be customised with additional options. In fact, the actual command format used for the demonstrations is

```sh
python main.py <PLANNER> maps/room1.png start <sx>,<sy> goal <sx>,<sy> -vv
```

to have a fix set of starting and goal points for consistent visualisation, but we omitted the start/goal options in the following commands for clarity.

### RRdT*

```sh
python main.py rrdt maps/room1.png -vv
```

<p align="center">
    <img width="600" height="auto" src="docs/images/rrdt.gif" alt="RRdT* Planner" />
</p>

### RRT*

```sh
python main.py rrt maps/room1.png -vv
```

<p align="center">
    <img width="600" height="auto" src="docs/images/rrt.gif" alt="RRT* Planner" />
</p>

### Bi-RRT*

```sh
python main.py birrt maps/room1.png -vv
```

<p align="center">
    <img width="600" height="auto" src="docs/images/birrt.gif" alt="Bi-RRT* Planner" />
</p>

### Informed RRT*

```sh
python main.py informedrrt maps/room1.png -vv
```

<p align="center">
<img width="600" height="auto" src="docs/images/informedrrt.gif" alt="Informed RRT* Planner" />
</p>

The red ellipse shown is the dynamic sampling area for Informed RRT*

### Others

There are also some other planners included in this repository. Some are preliminary planner that inspired RRdT*, some are planners with preliminary ideas, and some are useful for debugging.

## Reference to this repository

You can use the following citation if you use this repository for your research
```bibtex
@article{lai2021SbpEnv,
  doi = {10.21105/joss.03782},
  url = {https://doi.org/10.21105/joss.03782},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {66},
  pages = {3782},
  author = {Tin Lai},
  title = {sbp-env: A Python Package for Sampling-based Motion Planner and Samplers},
  journal = {Journal of Open Source Software}
}
```
# Contributing to `sbp-env`

We love your input! We welcome anyone make any form of contribution to `sbp-env`, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## Continuous integration with Github

We use Github to host code, track issues, pull requests, run automated test, as well as building documentations.

## Pull Requests

Pull requests are the best way to propose changes to the codebase:

1. Fork the repo and create your branch from `master`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Ensure the documentations are adequate.
6. Make sure your code lints (we use `black` for formatting).
7. Issue that pull request!

## Report bugs using Github's [issues](https://github.com/soraxas/sbp-env/issues)

We use GitHub issues to track public bugs. Report a bug by [opening a new issue](https://github.com/soraxas/sbp-env/issues); it's that easy!

## Write bug reports with detail, background, and sample code

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
- What you expected would happen
- What actually happens

## Testing the contributing code

You can check for all test cases with
```sh
pytest tests
```
in the root of the repository.

## Use a Consistent Coding Style

We use [Python Black](https://github.com/psf/black) for code formatting.

## Write documentations

We use [Sphinx](https://www.sphinx-doc.org/en/master/) to automatically build documentation directly from docstring. If you want to add more in-depth guide for any additional feature, feel free to add extra `.rst` file under `docs/` folder.

## License

By contributing, you agree that your contributions will be licensed under its MIT License.
---
title: 'sbp-env: A Python Package for Sampling-based Motion Planner and Samplers'
tags:
  - Python
  - motion planning
  - sampling-based planner
  - robotics
authors:
  - name: Tin Lai
    orcid: 0000-0003-0641-5250
    affiliation: '1'
affiliations:
  - name: School of Computer Science, The University of Sydney, Australia
    index: 1
date:
bibliography: master.bib
---

# Background

Sampling-based motion planning is one of the fundamental methods by which robots navigate and integrate with the real world [@elbanhawi2014_SampRobo].
Motion planning involves planning the trajectories of the actuated part of the robot, under various constraints, while avoiding collisions with surrounding obstacles.
Sampling-based motion planners (SBPs) are robust methods that avoid explicitly constructing the often intractable high-dimensional configuration space (C-Space).
Instead, SBPs randomly sample the C-Space for valid connections and iteratively build a roadmap of connectivity.
Most SBPs are guaranteed to find a solution if one exists [@kavraki1996_AnalProb], and such a planner is said to be _probabilistic complete_.
A further development for SBPs is _asymptotic optimality_[@elbanhawi2014_SampRobo]: a guarantee that the method will converge, in the limit, to the optimal solution.

SBPs are applicable to a wide range of applications.
Example include planning with arbitrary cost maps [@iehlCostmapPlanningHigh2012], cooperative multi-agent planning [@jinmingwuCooperativePathfindingBased2019], and planning in dynamic environments [@yershova2005_DynaRRTs].
On the one hand, researchers have focused on the algorithmic side of improving the graph or tree building [@lai2018_BalaGlob;@klemmRRTConnectFaster2015;@zhongTripleRrtsRobotPath2012;@elbanhawi2014_SampRobo;@lai2021lazyExperienceGraph;@lai2021rapidlyexploring].
On the other hand, the advancement of neural networks allows an abundance of learning approaches to be applied in SBPs [@strubAdaptivelyInformedTrees2020;@bagnell2014_ReinLear] and on improving the sampling distribution [@alcin2016_ExtrLear;@lai2020_BayeLoca;@lai2021plannerFlows;@laiLearningPlanOptimally2020;@lai2021diffSamp].


# Statement of need

The focus of motion planning research has been mainly on (i) the algorithmic aspect of the planner using different routines to build a connected graph and (ii) improving the sampling efficiency (with methods such as heuristic or learned distribution). Traditionally, robotic research focuses on algorithmic development, which has inspired several motion planning libraries written in C++, such as Move3D [@simeon2001move3d] and OMPL [@sucan2012open]. In particular, OMPL has been one of the most well-known motion planning libraries due to its versatility, and it has been a core part of the planning algorithm used in the MoveIt framework [@chitta2012moveit]. However, swapping the sampler within each planner is very restrictive, as planners are typically hard-coded to use a specific sampler. In addition, it is cumbersome to integrate any learning-based approach into a framework as there is only a limited number of choices of deep-learning libraries in C++.

Python has been a popular language to use in Machine Learning due to its rapid scripting nature. For example, PyTorch [@paszke2019pytorch] and Tensorflow [@abadi2016tensorflow] are two popular choices for neural network frameworks in Python. A large number of learning approaches are available as Python packages. It shall be noted that the aforementioned OMPL has Python bindings available; however, OMPL uses an outdated Py++ code generator, and every modification to the source code will require hours to updates bindings plus recompilation.
Some Python repositories are available that are dedicated to robotics motion planning [@sakai2018pythonrobotics]; however, most only showcase various planning algorithms, without an integrated environment and simulators.

![Implementation details on the class hierarchy structure of `sbp-env`.\label{fig:class-diagram}](class_diagram.png)

# Overview

We introduce `sbp-env`, a _sampling-based motion planners' testing environment_, as a complete feature framework to allow rapid testing of different sampling-based algorithms for motion planning.
`sbp-env` focuses on the flexibility of tinkering with different aspects of the framework, and it divides the main planning components into two main categories: (i) samplers and (ii) planners.
The division of the two components allows users to decouple them and focus only on the component that serves as the main focus of the research.
`sbp-env` has implemented the entire robot planning framework with multiple degrees-of-freedom, which allows benchmarking motion planning algorithms with the same planner under different backend simulators.
Separating the two components allows users to quickly swap out different components in order to test novel ideas.

Building the framework enables researchers to rapidly implement their novel ideas and validate their hypotheses.
In particular, users can define the environment using something as simple as an _image_, or as complicated as an _xml file_.
All samplers and planners can be added as a plugin system, and `sbp-env` will auto-discover newly implemented planners or samplers that have been added to the dedicated folders.

Figure \ref{fig:class-diagram} illustrates the hierarical structure of our package.
Our implementation of `sbp-env` define abstract interfaces for **sampler** and **planners**, from which all corresponding concrete classes must inherit.
In addition, there are classes that represent the full-body simulations of the environments and the corresponding visualisation methods.
Note that all visualisation can be turned off on-demand, which is beneficial when users benchmark their algorithms.
The docunmentation of `sbp-env` is available at [https://cs.tinyiu.com/sbp-env](https://cs.tinyiu.com/sbp-env).


# References
