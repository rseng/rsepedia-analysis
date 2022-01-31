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
Misc Classes
============

The following are classes and functions of which the details is not essential to use `sbp-env`.


Data structure
--------------

.. autoclass:: utils.common.Tree
  :members:
  :private-members:
  :inherited-members:

.. autoclass:: utils.common.Node
  :members:
  :private-members:
  :inherited-members:

.. autoclass:: utils.common.MagicDict
  :members: __getattr__
  :show-inheritance:

Randomness
----------------------------

The following are the random methods supported in `sbp-env`.

.. autodata:: randomness.SUPPORTED_RANDOM_METHODS


Drawing random samplings
^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`randomness.RandomnessManager:` is used by most samplers to draw random numbers.

.. autoclass:: randomness.RandomnessManager
  :members:
  :private-members:
  :inherited-members:


Drawing normally distributed samplings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`randomness.NormalRandomnessManager:` is used by :class:`planners.rrdtPlanner.RRdTSampler` to draw from a von Mises distribution.

.. autoclass:: randomness.NormalRandomnessManager
  :members:
  :private-members:
  :inherited-members:


Planners and Samplers Registry
------------------------------

The following function is used to register a new custom *planner*.

.. autofunction:: utils.planner_registry.register_planner

The following function is used to register a new custom *sampler*.

.. autofunction:: utils.planner_registry.register_sampler


Registry data structure
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: utils.planner_registry.PlannerDataPack
  :members:
  :private-members:
  :inherited-members:

.. autoclass:: utils.planner_registry.SamplerDataPack
  :members:
  :private-members:
  :inherited-members:




Adding new environment to test
==============================

With :code:`sbp-env`, it is extremely easy to add a new environment and run the
same set of planners on that new environment.
You can add an environment by supplying an image or xml file format in the
commandline via

.. code-block:: bash

    python main.py rrt <IMG_OR_XML> --engine ENGINE

and the required type depends on the simulator engine.
Both the 2D and 4D simulator will requires an image format to be used as the obstacle
space, while the :code:`klampt` simulator will requires xml format as specified in `klampt's doc
<http://motion.cs.illinois.edu/software/klampt/latest/pyklampt_docs/Manual-FileTypes.html>`_.

I recommend always using png as it is a lossless format, which is essential in for
the framework as it does not introduce artificial 'noise' to the image.
The framework uses the color of pixels to determine whether a given pixel is in
free-space :math:`C_\text{free}` or inside obstacles :math:`C_\text{obs}`.

Only white-pixels are considered as being in free-space. In other words, all
non-white pixels are taken as being inside :math:`C_\text{obs}`.
The image library used in Python take pixel as a value in-between 0-255 (or a tuple
of 3 values for RGB image).
We always ignore the alpha channel, and convert the given image to a gray-scale image.
Therefore, any pixels that are not of the value of 255 (represents the color white)
will be treated as obstacles.

The following figures showcase some of the Image-Space map that comes with
:code:`sbp-env`, and serve as examples of how regions within images are considered as
free-space.

.. sidebar:: Room

    .. Figure:: ../../maps/room1.png

.. sidebar:: Intel Lab

    .. Figure:: ../../maps/intel_lab.png

.. sidebar:: Random Obstacles

    .. Figure:: ../../maps/4d.png

.. sidebar:: Maze

    .. Figure:: ../../maps/maze1.png

Commandline Interface
=====================

The commandline interface can be accessed from the ``main.py`` file.
You can invoke the **cli** via

.. code-block:: console

    $ python main.py <PLANNER> <MAP> [--options]

All functionality of selecting planners and customising parameters are accessible from within the commandline interface.
If you just want to quickly test different planners in a newly created environment, then commandline would be the easiest option to get started.


Usage
-----

The full commandline interface is docunmented below:

----

.. automodule:: main.. _samplers_page:

C-space Samplers
================

The following are samplers that generates sampled configurations :math:`q \in C` in
*C-Space* with different strategies.

Random Policy Sampler
---------------------

.. autoclass:: samplers.randomPolicySampler.RandomPolicySampler
  :members:
  :private-members:
  :show-inheritance:

.. autodata:: randomness.SUPPORTED_RANDOM_METHODS
    :noindex:

Bidirectional Random Sampler
----------------------------

.. autoclass:: samplers.birrtSampler.BiRRTSampler
  :members:
  :private-members:
  :show-inheritance:


Informed Sampler
---------------------

.. autoclass:: samplers.informedSampler.InformedSampler
  :members:
  :private-members:
  :show-inheritance:


RRdT Particle Sampler
---------------------

.. autoclass:: planners.rrdtPlanner.RRdTSampler
  :members:
  :private-members:
  :show-inheritance:


PRM Sampler
---------------------

.. autoclass:: samplers.prmSampler.PRMSampler
  :members:
  :private-members:
  :show-inheritance:

Likelihood Sampler
---------------------

.. autoclass:: samplers.likelihoodPolicySampler.LikelihoodPolicySampler
  :members:
  :private-members:
  :show-inheritance:

Nearby Sampler
---------------------

.. autoclass:: samplers.nearbyPolicySampler.NearbyPolicySampler
  :members:
  :private-members:
  :show-inheritance:

Mouse Sampler
---------------------

.. autoclass:: samplers.mouseSampler.MouseSampler
  :members:
  :private-members:
  :show-inheritance:

Abstract Base Sampler
---------------------

There is also a special base Sampler that all motion planner should be derived from.

.. autoclass:: samplers.baseSampler.Sampler
  :members:
Support
=======

The easiest way to get help with the project is to open an issue in 
the `Github issue tracker`_.
If there is any details that the documentation is missing, feel free
to send in a PR to help improve the project.

Project page: https://github.com/soraxas/sbp-env

.. _Github issue tracker: https://github.com/soraxas/sbp-env/issues
.. _motion_planners_page:

Motion Planners
===============

The following are all of the available motion planners in ``sbp-env``.


RRT*
---------------------

.. image:: ../images/rrt.gif

.. autoclass:: planners.rrtPlanner.RRTPlanner
  :members:
  :private-members:
  :show-inheritance:

The sampler is registered as follows

.. literalinclude:: ../../samplers/randomPolicySampler.py
  :start-after: start register
  :end-before: finish register

and the planner is registered as follows

.. literalinclude:: ../../planners/rrtPlanner.py
  :start-after: start register
  :end-before: finish register


Bi-RRT*
---------------------

.. image:: ../images/birrt.gif

.. autoclass:: planners.birrtPlanner.BiRRTPlanner
  :members:
  :private-members:
  :show-inheritance:

The sampler is registered as follows

.. literalinclude:: ../../samplers/birrtSampler.py
  :start-after: start register
  :end-before: finish register

and the planner is registered as follows

.. literalinclude:: ../../planners/birrtPlanner.py
  :start-after: start register
  :end-before: finish register


RRdT*
---------------------

.. image:: ../images/rrdt.gif

.. autoclass:: planners.rrdtPlanner.RRdTPlanner
  :members:
  :private-members:
  :show-inheritance:

The sampler and planner is registered as follows

.. literalinclude:: ../../planners/rrdtPlanner.py
  :start-after: start register
  :end-before: finish register

Informedrrt-RRT*
---------------------

The *Informed-RRT\** motion planning behaves similar to :class:`planners.rrtPlanner
.RRTPlanner` before a first solution is found.
After an initial solution is found, this planner uses an ellipses heuristic to speed
up convergence rate of the resulting solution.

.. image:: ../images/informedrrt.gif

The bulk of the implementation occurs in
:class:`~samplers.informedrrtSampler.InformedRRTSampler`.
The Informed-RRT* is implemented by directly deriving the base motion planner as
:class:`~planners.rrtPlanner.RRTPlanner`,
and registering its sampler as the
:class:`~samplers.informedrrtSampler.InformedRRTSampler`.

The sampler is registered as follows

.. literalinclude:: ../../samplers/informedSampler.py
  :start-after: start register
  :end-before: finish register

and the planner is registered as follows

.. literalinclude:: ../../planners/informedrrtPlanner.py
  :start-after: start register
  :end-before: finish register


PRM*
---------------------

.. autoclass:: planners.prmPlanner.PRMPlanner
  :members:
  :private-members:
  :show-inheritance:

The sampler is registered as follows

.. literalinclude:: ../../samplers/prmSampler.py
  :start-after: start register
  :end-before: finish register

and the planner is registered as follows

.. literalinclude:: ../../planners/prmPlanner.py
  :start-after: start register
  :end-before: finish register

Likelihood Planner
---------------------

Refers to :class:`samplers.likelihoodPolicySampler.LikelihoodPolicySampler` for details.

The bulk of the implementation occurs in :class:`~samplers.likelihoodPolicySampler.LikelihoodPolicySampler`.
The Informed-RRT* is implemented by directly deriving the base motion planner as :class:`~planners.rrtPlanner.RRTPlanner`,
and registering its sampler as the :class:`~samplers.likelihoodPolicySampler.LikelihoodPolicySampler`.

The sampler is registered as follows

.. literalinclude:: ../../samplers/likelihoodPolicySampler.py
  :start-after: start register
  :end-before: finish register

and the planner is registered as follows

.. literalinclude:: ../../planners/likelihoodPlanner.py
  :start-after: start register
  :end-before: finish register

Nearby Planner
---------------------

Refers to :class:`samplers.nearbyPolicySampler.NearbyPolicySampler` for details.

The sampler is registered as follows

.. literalinclude:: ../../samplers/nearbyPolicySampler.py
  :start-after: start register
  :end-before: finish register

and the planner is registered as follows

.. literalinclude:: ../../planners/nearbyPlanner.py
  :start-after: start register
  :end-before: finish register

Mouse Planner
---------------------

Refers to :class:`samplers.mouseSampler.MouseSampler` for details.

The sampler is registered as follows

.. literalinclude:: ../../samplers/mouseSampler.py
  :start-after: start register
  :end-before: finish register

and the planner is registered as follows

.. literalinclude:: ../../planners/mousePlanner.py
  :start-after: start register
  :end-before: finish register

Abstract Base Planner
---------------------

There is also a special base planner that all motion planner should be derived from.

.. autoclass:: planners.basePlanner.Planner
  :members:
  :private-members:
  :show-inheritance:
  :inherited-members:
Planning Scene Visualisers
==========================

Each simulator engines would needs a corresponding visualiser for the framework to
know what to display. The simulator and visualiser are decoupled, for the purpose of
benchmarking. During benchmark, one would want to disable the visualisation because
the drawing (or the preparation) of visual elements would slow down the planning
process. Even if the slowdown is consistent across different types of planners, it
always best to have the ability to obtain quick results.

Each planner can register a visualisation function for a specific type of simulator.
You can **disable the visualisation** of the planning problem, regardless of simulator
type, with the following.

.. prompt:: bash

    python main.py <PLANNER> <MAP> start <q_1>,...,<q_n> goal <p_1>,...,<p_n> --no-display

For example:

.. prompt:: bash

    python main.py rrt maps/intel_lab.png start 25,25 goal 225,225 --no-display

.. important::
    Notice that in the example above, the argument `start` and `goal` are directly
    provided in the prompt. This is because with the `--no-display` flag there won't be
    any GUI for user to select the start/goal pair. Therefore, it is necessary for user
    to directly supply the starting and target configurations.

The following sections showcase the different classes of visualiser.

2D Image Space Visualiser
-------------------------

.. autoclass:: visualiser.PygameEnvVisualiser
  :members:
  :private-members:
  :show-inheritance:


.. autoclass:: visualiser.PygamePlannerVisualiser
  :members:
  :private-members:
  :show-inheritance:


.. autoclass:: visualiser.PygameSamplerVisualiser
  :members:
  :private-members:
  :show-inheritance:


4D Image Space Manipulator Visualiser
-------------------------------------

The 4D manipulator visualiser uses the same
:class:`visualiser.PygameEnvVisualiser`,
:class:`visualiser.PygamePlannerVisualiser` and
:class:`visualiser.PygameSamplerVisualiser` as the `2D Image Space Visualiser`_
to visualise the planning scene. However, it performs the necessary kinematics transformations to translates configurations :math:`q \in C` to worldspace :math:`x \in \mathcal{W}`
with the help of the collision checker :class:`collisionChecker.RobotArm4dCollisionChecker`.

.. warning::

    The settings of the robot arm's joint length is configurable from the collision
    checker's construction, but has not been exposed to the commandline interface yet.
    See :meth:`collisionChecker.RobotArm4dCollisionChecker.__init__`



3D Object-based Visualiser
--------------------------


.. autoclass:: visualiser.KlamptEnvVisualiser
  :members:
  :private-members:
  :show-inheritance:


.. autoclass:: visualiser.KlamptPlannerVisualiser
  :members:
  :private-members:
  :show-inheritance:


.. autoclass:: visualiser.KlamptSamplerVisualiser
  :members:
  :private-members:
  :show-inheritance:




Currently there are three engines used in ``sbp-env`` that are used to simulate collisions and the planning environment.

The first one is a **2D**  engine with the configuration space (*C-Space*) as :math:`C \subseteq \mathbb{R}^2`.
This engine is very easy to use and very easy to add new environments.
This uses a the :class:`collisionChecker.ImgCollisionChecker` for environment simulation and with the :class:`visualiser.PygameEnvVisualiser` to visualise the planning scene.
Sampling-based Motion Planners' Testing Environment
===================================================

Here is where ``sbp-env``'s documentation, the Sampling-based Motion Planner
Benchmark suite, lies.
To quickly get started, see the :ref:`quick_start` guide.
For various documentation of th implemented planners and samplers, checkout the
corresponding sections from the side bar.

The code is hosted in a `repository <https://github.com/soraxas/sbp-env>`_, and the
implementation requires *Python 3.7+*. The code comes with some example environments,
and in addition, several state-of-the-art sampling-based :ref:`motion_planners_page`
and :ref:`samplers_page` are implemented.

.. toctree::
   :caption: Guide
   :maxdepth: 3

   quick-start
   cli
   add-environment

.. toctree::
   :caption: Implemented SBPs
   :maxdepth: 3

   planners
   samplers

.. toctree::
   :caption: Planning environment
   :maxdepth: 3

   engines
   visualisers

.. toctree::
   :caption: Others
   :maxdepth: 3

   misc
   support


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Simulator Engine
================

Currently there are three engines used in ``sbp-env`` that are used to simulate collisions and the planning environment.

2D Image Space
--------------

The first one is a **2D** image space engine with the configuration space (*C-Space*) as :math:`C \subseteq \mathbb{R}^2`.
This engine is very easy to use and very easy to add new environments.
This uses the :class:`collisionChecker.ImgCollisionChecker` for environment simulation and with the :class:`visualiser.PygameEnvVisualiser` to visualise the planning scene.

.. prompt:: bash

    python main.py rrdt maps/room1.png

.. sidebar:: 2D Image simulator

    .. Figure:: ../images/rrdt.gif

This simulator uses images as its *C-Space* by treating each white pixels as configurations in :math:`C_\text{free}`,
and each non-white pixels as configurations in :math:`C_\text{obs}`.
The two dimensional space of the *C-Space*

.. math::
    q \equiv [x_0, x_1] \in C \subseteq \mathbb{R}^2

would be taken as the *width* and the *height* of the image

.. math::
    0 \le x_0 < \mathcal{I}_\text{Width} \\
    0 \le x_1 < \mathcal{I}_\text{Height}

where :math:`\mathcal{I}` denotes the image space.

.. important::
    Internally, all pixels with an **RGB value** of `(255, 255, 255)` would be treated as :math:`q \in C_\text{free}`,
    and any non-white pixels value will be treated as :math:`q \in C_\text{obs}` (e.g. `(10, 20, 10)`, `(0, 0, 0)`, etc.).
    The alpha channel would not be considered.

.. autoclass:: collisionChecker.ImgCollisionChecker
  :members:
  :private-members:
  :show-inheritance:


4D Robot Arm
--------------

This simulator internally depends on :class:`collisionChecker.ImgCollisionChecker` to check for collision in the image space.

.. prompt:: bash

    python main.py rrt maps/4d.png -s .5 --engine 4d


.. sidebar:: 4D Robot Arm simulator

    .. Figure:: ../images/robot-arm4d.gif

In contrast to the former, this simulator not only check for point mass collisions,
but performs forward kinematic on the given joints configuration to obtain body
points in world-space :math:`\mathcal{W}\subseteq \mathbb{R}^2`.
Since the robot arm contains two configurable joints, the full configuration space :math:`C` is given by

.. math::
    q \equiv [x_0, x_1, r_0, r_1] \in C \subseteq \mathbb{R}^4

where

.. math::
    \begin{aligned}
        0      & \le  x_0  <  \mathcal{I}_\text{Width} \\
        0      & \le  x_1  <  \mathcal{I}_\text{Height} \\
        - \pi  & \le  r_0  <  \pi \\
        - \pi  & \le  r_1  <  \pi
    \end{aligned}

and we use :math:`r_0, r_1` to denote the rotational angles of the first and second joints respectively.


.. important::
    Similar to :class:`collisionChecker.ImgCollisionChecker`, all white pixels will be within :math:`q \in C_\text{free}` and vice versa.
    The body points obtained by the forward kinematic on :math:`r_0, r_1` would be used to check collision in :math:`\mathcal{W}` to ensure the entire body is in free space.

.. autoclass:: collisionChecker.RobotArm4dCollisionChecker
  :members:
  :private-members:
  :show-inheritance:


:math:`n`-D Manipulator
-----------------------

This simulator internally depends on `klampt` package to check for collision in the
the 3D space.

.. prompt:: bash

    python main.py rrt klampt_data/tx90blocks.xml --engine klampt


.. sidebar:: :math:`n`-D Manipulator simulator

    .. Figure:: ../images/klampt-rrdt.gif

In contrast to the former, this simulator not only check for point mass collisions,
but performs forward kinematic with full body collision on the given joints
configuration to check validity in world-sapce :math:`\mathcal{W}\subseteq
\mathbb{R}^3` and :math:`C \subseteq \mathbb{R}^d`.
The full configuration space :math:`C` is given by

.. math::
    q \equiv [r_0, \ldots, r_{d-1}] \in C \subseteq \mathbb{R}^d

where

.. math::
    \begin{aligned}
        - \pi  & \le  r_i  <  \pi \quad \forall i \in \{0, \ldots, d-1\}
    \end{aligned}

and we use :math:`r_0, r_1` to denote the rotational angles of the first and second joints respectively.


.. important::
    This simulator is based on the upstream `klampt` and the upstream repository
    might update its api without notice. The current implementation is based on
    `Klampt==0.8.7`

.. autoclass:: collisionChecker.KlamptCollisionChecker
  :members:
  :private-members:
  :show-inheritance:



.. _quick_start:

Quick Start
===========

Installation
----------------------

Download the latest code and install required dependencies listed inside
`requirements.txt` (recommend to always isolate your python workspace with virtual
environment, e.g., conda/pyenv)

.. code-block:: bash

    # get code from remote
    git clone https://github.com/soraxas/sbp-env
    cd sbp-env

    # optional, using virtual env
    conda create --name my_ws python=3.8
    conda activate my_ws

You would first need to install the dependencies with

.. code-block:: bash

    # install dependencies
    pip3 install -r requirements.txt
    # optional dependencies for manipulator simulator
    pip3 install -r requirements_klampt.txt

The quickest way to test out the capability of a planner is to experiment with a *2D
Image Space* planner.

.. sidebar:: 2D Image simulator

    .. Figure:: ../images/rrdt.gif

.. code-block:: bash

    python main.py rrt maps/room1.png

The syntax always requires a positional argument which specify the planning
environment. In the case of the 2D engine or 4D rover arm, it should be a path to an
image file (png/jpg). For the case of the :math:`n`-D manipulator, it should be a path
to an xml file as specified in `Klampt
<https://github.com/krishauser/Klampt>`_.


Using other simulator
----------------------

The type of simulator can be set with the :code:`--engine` flag. For example, the 4D
robot arm simulator can be used with the :code:`4d` value.

.. sidebar:: 4D Robot Arm simulator

    .. Figure:: ../images/robot-arm4d.gif

.. code-block:: bash

    python main.py rrt maps/4d.png --engine 4d

You can also set the length of the rover arm (a list of two float that corresponds to
the 2 sticks' length) with the following.

.. code-block:: bash

    python main.py rrt maps/room1.png --engine 4d --4d-robot-lengths 15,15

Notice that the argument for :code:`--4d-robot-lengths` must be
comma-separated numbers, each of the number represent the length of the rover arm.

For the 4d rover arm simulator, :code:`sbp-env` will visualise each node with the
actual arm that is currently used for collision-checking.
The two arms are currently being displayed in orange and cyan colour for
distinguishing them when the display is cluttered.
The colour are currently non-configurable, and only matters to the visualisation of
the planning horizon.



.. sidebar:: :math:`n`-D Manipulator simulator

    .. Figure:: ../images/klampt-rrdt.gif

You can also launch the 3D simulator environment with the following command, given
that you had already installed the optional :code:`klampt` dependencies that was
mentioned previously. Note that you will need to supply the environment file as
specified in `klampt's doc
<http://motion.cs.illinois.edu/software/klampt/latest/pyklampt_docs/Manual-FileTypes.html>`_ as they are passed directly to the :code:`klampt` simulator backend.

We have included a complementary set of klampt data `here
<https://github.com/soraxas/sbp-env/releases/tag/klampt-data>`_, which is adapted
from the official :code:`klampt` `example repository
<https://github.com/krishauser/Klampt-examples>`_.

You can use the following commands to quickly download the klampt data files:

.. code-block:: bash

    cd sbp-env

    # download klampt-data from remote
    wget https://github.com/soraxas/sbp-env/releases/download/klampt-data/klampt_data.tar.gz
    # extract
    tar xzf klampt_data.tar.gz

And last but not least, you can start the :code:`klampt` simulator with

.. code-block:: bash

    python main.py rrt klampt_data/tx90blocks.xml --engine klampt



.. important::
    If you don't have :code:`OpenGL` installed, you might also need to install it for
    using the :code:`klampt` simulator.
    If you are using pip, you can install it with

    .. code-block:: bash

        pip3 install pyopengl pyqt5


Saving the planner statistics
------------------------------

During the planning episode, :code:`sbp-env` will keep track of various statistics which
are beneficial to compare the performance between planner. The statistics will always
be displayed as part of the :code:`tqdm` progress bar.

You can save the statistics to a :code:`.csv` file with the :code:`--save-output` flag,
e.g.

.. code-block:: bash

    python main.py rrt maps/room1.png start 100,100 goal 350,350 --save-output

which would save the output to a timestamped :code:`.csv` file under the :code:`runs/`
folder by default. You can also customise the output folder with
:code:`--output-dir=MY_FOLDER`.

The recorded statistics have the following meanings:
    - :code:`nodes`: The number of nodes
    - :code:`time`: Timestamp
    - :code:`cc_feasibility`: The number of collision-checks for feasibility (node)
    - :code:`cc_visibility`: The number of collision-checks for visibility (edge)
    - :code:`invalid_feasibility`: The number of invalid feasibility checks
    - :code:`invalid_visibility`: The number of invalid visibility checks
    - :code:`c_max`: The current cost of the solution trajectory
{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}
  
   {% block attributes %}
   {% if attributes %}
   .. rubric:: Module Attributes

   .. autosummary::
      :toctree:                                          
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:                                          
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:                                          
      :template: custom-class-template.rst               
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   .. autosummary::
      :toctree:                                          
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: custom-module-template.rst                 
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:               
   :show-inheritance:     
   :inherited-members:   

   {% block methods %}
   .. automethod:: __init__

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
