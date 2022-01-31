# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [1.0.1] - 2021-12-20
## Added
- Classic field oriented controllers for induction motors
- Uniform initialization of the WienerProcessReferenceGenerator

## Changed
- Reduced the dynamics of the reference signals in several speed control environments
- Changed the default ode-solver of all environments to the ScipyOdeSolver

## Fixed
- gym version compatibility for all versions >= 0.21.0
- Docs: m2r dependency to m2r2. Enables compatibility with latest sphinx versions.
- matplotlib compatibility with versions >= 3.5.0
- Bugfix in the stable_baselines3_dqn_disc_pmsm_example.ipynb example notebook
- Bugfix in the jacobian of the ConstantSpeedLoad

## [1.0.0] - 2021-08-23
### Added
- classic controllers in examples with complete makeover
- Possibility to seed all random number generators at once with a unified API - Reproduciblity improved.

### Changed
#### Environment IDs
- The environments have changed the IDs for their creation. 
- Furthermore, environments specialized for torque (TC), current (CC) and speed (SC) control tasks have been introduced.
- The ids are now structured as follows:
{_Cont/Finite_}-{_CC/TC/SC_}-_motorType_-v0
- Details can be looked up on the code documentation pages.
#### gem.make() parameters
The arbitrary environment creation keyword arguments that were passed to every subcomponent of the environment
have been removed. Instead, there is a fixed set of keyword parameters for every environment described in the 
API-documentation of every environment.

#### MPC example
- Visualization improvements
- fix: State variables were updated

#### Miscellaneous
- Documentation was updated
- Changed all DC Envs to FourQuadrantConverters per default
- Adjusted the dynamics of the speed references in DC environments
- Adjusted the plots for better visibility of single points


### Removed
- simple_controllers.py in examples
- Tensorforce tutorial due to deprecation


## [0.3.1] - 2020-12-18
### Added
- Constraint monitor class

### Changed
- Visualization framework
# Gym Electric Motor
![](docs/plots/Motor_Logo.png)


[**Overview paper**](https://joss.theoj.org/papers/10.21105/joss.02498)
| [**Reinforcement learning paper**](https://arxiv.org/abs/1910.09434)
| [**Quickstart**](#getting-started)
| [**Install guide**](#installation)
| [**Reference docs**](https://upb-lea.github.io/gym-electric-motor/)
| [**Release notes**](https://github.com/upb-lea/gym-electric-motor/releases)

[![Build Status](https://github.com/upb-lea/gym-electric-motor/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/upb-lea/gym-electric-motor/actions/workflows/build_and_test.yml)
[![codecov](https://codecov.io/gh/upb-lea/gym-electric-motor/branch/master/graph/badge.svg)](https://codecov.io/gh/upb-lea/gym-electric-motor)
[![PyPI version shields.io](https://img.shields.io/pypi/v/gym-electric-motor.svg)](https://pypi.python.org/pypi/gym-electric-motor/)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/upb-lea/gym-electric-motor/blob/master/LICENSE)
[![DOI Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4355691.svg)](https://doi.org/10.5281/zenodo.4355691)
[![DOI JOSS](https://joss.theoj.org/papers/10.21105/joss.02498/status.svg)](https://doi.org/10.21105/joss.02498)

## Overview
The gym-electric-motor (GEM) package is a Python toolbox for the simulation and control of various electric motors.
It is built upon [OpenAI Gym Environments](https://gym.openai.com/), and, therefore, can be used for both, classical control simulation and [reinforcement learning](https://github.com/upb-lea/reinforcement_learning_course_materials) experiments. It allows you to construct a typical drive train with the usual building blocks, i.e., supply voltages, converters, electric motors and load models, and obtain not only a closed-loop simulation of this physical structure, but also a rich interface for plugging in any decision making algorithm, from linear feedback control to [Deep Deterministic Policy Gradient](https://spinningup.openai.com/en/latest/algorithms/ddpg.html) agents.

## Getting Started
An easy way to get started with GEM is by playing around with the following interactive notebooks in Google Colaboratory. Most important features of GEM as well as application demonstrations are showcased, and give a kickstart for engineers in industry and academia.

* [GEM cookbook](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master//examples/environment_features/GEM_cookbook.ipynb)
* [Keras-rl2 example](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/reinforcement_learning_controllers/keras_rl2_dqn_disc_pmsm_example.ipynb)
* [Stable-baselines3 example](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/reinforcement_learning_controllers/stable_baselines3_dqn_disc_pmsm_example.ipynb)
* [MPC  example](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/model_predictive_controllers/pmsm_mpc_dq_current_control.ipynb)

There is a list of [standalone example scripts](examples/) as well for minimalistic demonstrations.

A basic routine is as simple as:
```py
import gym_electric_motor as gem

if __name__ == '__main__':
    env = gem.make("Finite-CC-PMSM-v0")  # instantiate a discretely controlled PMSM
    env.reset()
    for _ in range(10000):
        env.render()  # visualize environment
        (states, references), rewards, done, _ =\ 
        	env.step(env.action_space.sample())  # pick random control actions
        if done:
            (states, references) = env.reset()
    env.close()
```



## Installation
- Install gym-electric-motor from PyPI (recommended):

```
pip install gym-electric-motor
```

- Install from Github source:

```
git clone git@github.com:upb-lea/gym-electric-motor.git 
cd gym-electric-motor
# Then either
python setup.py install
# or alternatively
pip install -e .
```

## Building Blocks
A GEM environment consists of following building blocks:
- Physical structure:
   - Supply voltage
   - Converter
   - Electric motor
   - Load model
- Utility functions for reference generation, reward calculation and visualization
 
### Information Flow in a GEM Environment
![](docs/plots/SCML_Overview.png)

Among various DC-motor models, the following AC motors - together with their power electronic counterparts - are available:
- Permanent magnet synchronous motor (PMSM), 
- Synchronous reluctance motor (SynRM)
- Squirrel cage induction motor (SCIM)
- Doubly-fed induction motor (DFIM)

The converters can be driven by means of a duty cycle (continuous control set) or switching commands (finite control set). 

### Citation
A white paper for the general toolbox in the context of drive simulation and control prototyping can be found in the [Journal of Open Sorce Software (JOSS)](https://joss.theoj.org/papers/10.21105/joss.02498). Please use the following BibTeX entry for citing it:
```
@article{Balakrishna2021,
    doi = {10.21105/joss.02498},
    url = {https://doi.org/10.21105/joss.02498},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {58},
    pages = {2498},
    author = {Praneeth {Balakrishna} and Gerrit {Book} and Wilhelm {Kirchgässner} and Maximilian {Schenke} and Arne {Traue} and Oliver {Wallscheid}},
    title = {gym-electric-motor (GEM): A Python toolbox for the simulation of electric drive systems},
    journal = {Journal of Open Source Software}
}

```

A white paper for the utilization of this framework within reinforcement learning is available at [IEEE-Xplore](https://ieeexplore.ieee.org/document/9241851) (preprint: [arxiv.org/abs/1910.09434](https://arxiv.org/abs/1910.09434)). Please use the following BibTeX entry for citing it:
```
@article{9241851,  
    doi={10.1109/TNNLS.2020.3029573}
    author={A. {Traue} and G. {Book} and W. {Kirchgässner} and O. {Wallscheid}},  
    journal={IEEE Transactions on Neural Networks and Learning Systems},   
    title={Toward a Reinforcement Learning Environment Toolbox for Intelligent Electric Motor Control},   
    year={2020},  volume={},  number={},      
    pages={1-10},  
}
```

### Running Unit Tests with Pytest
To run the unit tests ''pytest'' is required.
All tests can be found in the ''tests'' folder.
Execute pytest in the project's root folder:
```
>>> pytest
```
or with test coverage:
```
>>> pytest --cov=./
```
All tests shall pass.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at kirchgaessner(at)lea.uni-paderborn(dot)de. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq

# How to Contribute

First off, thank you for considering contributing to GEM. It's people like you that make GEM a great tool!

Following these guidelines helps to communicate that you respect the time of the developers managing and developing this open source project. In return, they should reciprocate that respect in addressing your issue, assessing changes, and helping you finalize your pull requests.

## What we are looking for

Improving documentation, bug triaging, or writing tutorials are all helpful contributions we are looking forward to.
Especially submitting bug reports, feature requests, or code which can be incorporated directly into GEM is what we would love to receive.

## Code of Conduct

The goal is to maintain a diverse community that's pleasant for everyone. That's why we would greatly appreciate it if everyone contributing to and interacting with the community also followed our [Code of Conduct](https://github.com/upb-lea/gym-electric-motor/blob/master/CODE_OF_CONDUCT.md).

## Asking for support and requesting features
Please create issues with the __question__ label for support questions and for feature requests.
Even if you want to contribute a new feature yourself please state your intention in the issue tracker first, such that the maintainers can give valuable feedback on how to tackle this in the most appropriate way.

## Reporting Bugs
The best way to report an issue and to ensure a timely response is to use the issue tracker.

* Create a GitHub account. 
    * You need to create a GitHub account to be able to create new issues and participate in the discussion.

* Determine if your bug is really a bug.
    * You shouldn't file a bug if you're requesting support. For that you can use the contact [here](https://ei.uni-paderborn.de/en/lea/team/arbeitsgruppe/team/).

* Make sure your bug hasn't already been reported.
    * Search through the Issue tracker. If a bug like yours was found, check if you have new information that could be reported to help the developers fix the bug.

* Check if you're using the latest version.
    * A bug could be fixed by some other improvements and fixes - it might not have an existing report in the bug tracker. Make sure you're using the latest releases.

* Collect information about the bug.
    * To have the best chance of having a bug fixed, we need to be able to easily reproduce the conditions that caused it. Most of the time this information will be from a Python traceback message, though some bugs might be in design, spelling or other errors on the website/docs/code.

    * If the error is from a Python traceback, include it in the bug report.

    * We also need to know what platform you're running (Windows, macOS, Linux, etc.), the version of your Python interpreter, and the version of GEM, and related packages that you were running when the bug occurred.

    * If you're reporting a race condition or a deadlock, tracebacks can be hard to get or might not be that useful. Try to inspect the process to get more diagnostic data.

    * Your issue might be tagged as Needs Test Case. A test case represents all the details needed to reproduce what your issue is reporting. A test case can be some minimal code that reproduces the issue or detailed instructions and configuration values that reproduces said issue.

* Submit the bug.
    * By default GitHub will email you to let you know when new comments have been made on your bug. In the event you've turned this feature off, you should check back on occasion to ensure you don't miss any questions a developer trying to fix the bug might ask.


## Versions and Tags

Version numbers consists of a major version, minor version and a bug-fix version.

Tags are used exclusively for tagging releases. A release tag is named with the format vX.Y.Z -- for example v2.3.1.
Experimental releases contain an additional identifier vX.Y.Z-id -- for example v3.0.0-rc1.
Experimental tags may be removed after the official release.

## Examples

In order to provide better accessibility of this simulation library we created several examples that showcase the usage of GEM's interface and features.
The presented examples can be classified either as classical control approaches, as optimal control using model predictive control and model-free reinforcement-learning control examples, or as feature showcases (where the focus lies in introducing interesting and useful builtins).

### Installation, Setup and Interface
- [GEM_cookbook.ipynb](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/environment_features/GEM_cookbook.ipynb): a basic tutorial-style notebook that presents the basic interface and usage of GEM

### Classical Control
- [Classic controller suite](classic_controllers/classic_controllers.py): a control algorithm library using standard, expert-driven techniques with automatic and manual controller tuning for an increasingly growing number of motor and inverter combinations
- [Classic controller examples](classic_controllers): a set of examples applying the classic controller suite (e.g., field-oriented control for permanent synchronous motors)

### Advanced Control
- [gekko_mpc_cont_pmsm_example.ipynb](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/model_predictive_controllers/gekko_mpc_cont_pmsm_example.ipynb): a model predictive control solution for the currents of the three-phase permanent magnet synchronous motor on a continuous-control-set

### Reinforcement-Learning Control
- [dqn_series_current_control.py](reinforcement_learning_controllers/dqn_series_current_control.py): a deep Q-value network reinforcement-learning control approach for finite-control-set current control of a series DC motor
- [ddpg_pmsm_dq_current_control.py](reinforcement_learning_controllers/ddpg_pmsm_dq_current_control.py): a deep deterministic policy gradient reinforcement-learning control approach applied to the current control of a permanent magnet synchronous motor within the $dq$-frame with continuous-control-set
- [ddpg_series_omega_control.py](reinforcement_learning_controllers/ddpg_series_omega_control.py): a deep deterministic policy gradient reinforcement-learning control approach applied to the speed control of a series DC motor with continuous-control-set
- [keras_rl2_dqn_disc_pmsm_example.ipynb](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/reinforcement_learning_controllers/keras_rl2_dqn_disc_pmsm_example.ipynb): a tutorial-style notebook that presents the usage of GEM in conjunction with [Keras_RL2](https://github.com/wau/keras-rl2) in the context of deep Q learning current control of a permanent magnet synchronous motor
- [stable_baselines3_dqn_disc_pmsm_example.ipynb](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/reinforcement_learning_controllers/stable_baselines3_dqn_disc_pmsm_example.ipynb): a tutorial-style notebook that presents the usage of GEM in conjunction with [Stable Baselines3](https://github.com/DLR-RM/stable-baselines3) in the context of deep Q learning current control of a permanent magnet synchronous motor

### Feature Showcases
- [external_speed_profile.py](environment_features/external_speed_profile.py): presents a builtin feature that can be used to define arbitrary speed profiles, which is useful when e.g. investigating generator operation (where mechanical and thus electrical frequency is determined by external means) 
- [userdefined_initialization.py](environment_features/userdefined_initialization.py): presents a builtin feature that allows the user to determine the initial state of the motor, which comes in handy when e.g. using exploring starts in reinforcement learning applications
- [scim_ideal_grid_simulation.py](environment_features/scim_ideal_grid_simulation.py): simulates the start-up behavior of the squirrel cage induction motor connected to an ideal three-phase grid. 
The state and action space is continuous.
Running the example will create a formatted plot that show the motors angular velocity, the drive torque, the applied voltage in three-phase abc-coordinates and the measured current in field-oriented dq-coordinates.---
title: 'gym-electric-motor (GEM): A Python toolbox for the simulation of electric drive systems'
tags:
  - Python
  - electric drive control
  - electric motors
  - OpenAI Gym
  - power electronics
  - reinforcement learning
authors:
  - name: Praneeth Balakrishna
    affiliation: 1
    
  - name: Gerrit Book
    affiliation: 1
    
  - name: Wilhelm Kirchgässner
    orcid: 0000-0001-9490-1843
    affiliation: 1
    
  - name: Maximilian Schenke
    orcid: 0000-0001-5427-9527
    affiliation: 1
    
  - name: Arne Traue
    affiliation: 1
    
  - name: Oliver Wallscheid
    orcid: 0000-0001-9362-8777
    affiliation: 1
    
affiliations:
 - name: Department of Power Electronics and Electrical Drives, Paderborn University, Germany
   index: 1
date: 28. May 2020
bibliography: Literature.bib
---

# Summary

The ``gym-electric-motor`` (``GEM``) library provides simulation environments for electrical drive systems and, therefore, allows to easily design and analyze drive control solutions in Python.
Since ``GEM`` is strongly inspired by OpenAI's ``gym`` [@gym-whitepaper], it is particularly well-equipped for (but not limited to) applications in the field of reinforcement-learning-based control algorithms. 
In addition, the interface allows to plug in any expert-driven control approach, such as model predictive control, to be tested  and to perform benchmark comparisons. 
The ``GEM`` package includes a wide variety of motors, power electronic converters and mechanical load models that can be flexibly selected and parameterized via the API. 
A modular structure allows additional system components to be included in the simulation framework.

# Statement of Need

Electric drive systems and their control are an important topic in both academic and industrial research due to their worldwide usage and deployment. 
Control algorithms for these systems have usually been designed, parameterized and tested within ``MATLAB - Simulink`` [@MathWorks], which is developed and promoted specifically for
such engineering tasks. 
In the more recent past, however, commercial software like ``MATLAB`` has difficulties to stay on par with state-of-the-art concepts for scientific modeling and the flexibility offered by open-source libraries that are available for more accessible programming languages like Python. 
Consequently, a Python-based drive simulation framework like ``GEM`` is an evident step in order to accelerate corresponding control research and development.
Specifically, the latest efforts concerning industrial application of reinforcement-learning control algorithms heavily depend on Python packages like ``Keras`` [@Chollet2015], ``Tensorflow`` [@tensorflow2015-whitepaper] or ``PyTorch``[@NEURIPS2019_9015]. 
Hence, the built-in OpenAI ``gym`` interface allows to easily couple ``GEM`` to other open-source reinforcement learning toolboxes such as ``Stable Baselines3`` [@stable-baselines3], ``TF-Agents`` [@TFAgents] or ``keras-rl`` [@plappert2016kerasrl].

Providing easy access to the non-commercial, open-source ``GEM`` library allows users from any engineering domain to include accurate drive models into their simulations, also beyond the topic of control applications.
Considering the prevalence of commercial software like ``MATLAB`` for educational purposes, a free-of-charge simulation alternative that does not force students or institutions to pay for licenses, has great potential to support and encourage training of new talents in the field of electrical drives and neighbouring domains (e.g. power electronics or energy systems).
``GEM`` has already been used in graduate courses on reinforcement learning [@rl-lecture].

# Related software

Due to the strong dependence of downstream industrial development on simulated environments there is a comprehensive variety of commercial software that enables numerical analysis of every facet of electric drives. 
To name just a few, ``MATLAB - Simulink`` is probably the most popular software environment for numerical analysis in engineering.
Herein, ``MATLAB`` is providing for a scientific calculation framework and ``Simulink`` for a model-driven graphical interface with a very large field of applications. 
Examples that are designed for real-time capability (e.g., for hardware-in-the-loop prototyping) can be found in ``VEOS`` [@dSPACE] or ``HYPERSIM`` [@OPAL-RT].
Non-commercial simulation libraries exist, but they rarely come with predefined system models. 
An exemplary package from this category is ``SimuPy`` [@Margolis], which provides lots of flexibility for the synthesis of generic simulation models, but also requires the user to possess the necessary expert knowledge in order to implement a desired system model. 
Likewise, general purpose component-oriented simulation frameworks like ``OpenModelica`` [@OSMC2020] or ``XCos`` [@Scilab2020] can be used for setting up electrical drive models, too, but this requires expert domain knowledge and out-of-the-box Python interfaces (e.g., for reinforcement learning) are not available. 

In the domain of motor construction it is furthermore interesting to observe the behavior of magnetic and electric fields within a motor (simulation of partial differential equations).
Corresponding commercial simulation environments, like ``ANSYS Maxwell`` [@ANSYS], ``Motor-CAD`` [@MotorDesignLtd] or ``MotorWizard`` [@ElectroMagneticWorks] and the exemplary non-commercial alternative ``FEMM`` [@Meeker] are very resource and time consuming because they depend on the finite element method, which is a spatial discretization and numerical integration procedure. 
Hence, these software packages are usually not considered in control development, and complement ``GEM`` at most. 
This particularly applies in the early control design phase when researching new, innovative control approaches (rapid control prototyping) or when students want to receive quasi-instantaneous simulation feedbacks. 

# Package Architecture

The ``GEM`` library models an electric drive system by its four main components: voltage supply, power converter, electric motor and mechanical load. 
The general structure of such a system is depicted in \autoref{fig:SCML_system}. 

![Simplified structure diagram of an electric drive system\label{fig:SCML_system}](../plots/SCML_Setting.eps)

The __voltage supply__ provides the necessary power that is used by the motor. 
It is modeled by a fixed supply voltage $u_\mathrm{sup}$, which allows to monitor the supply current into the converter.
A __power electronic converter__ is needed to supply the motor with electric power of proper frequency and magnitude, which commonly includes the conversion of the supply's direct current to alternating current. 
Typical drive converters exhibit switching behavior: there is a finite set of different voltages that can be applied to the motor, depending on which switches are open and which are closed. 
Besides this physically accurate view, a popular modeling approach for switched mode converters is based on dynamic averaging of the applied voltage $u_\mathrm{in}$, rendering the voltage a continuous variable.
Both of these modeling approaches are implemented and can be chosen freely, allowing usage of control algorithms that operate on a finite set of switching states or on continuous input voltages.
The __electric motor__ is the centerpiece of every drive system. 
It is described by a system of ordinary differential equations (ODEs), which represents the motor's electrical behavior. 
In particular, the domain of three-phase drives makes use of coordinate transformations to view these ODEs in the more interpretable frame of field-oriented coordinates. 
In ``GEM``, both, the physically accurate three-phase system ($abc$-coordinates) and the simplified, two-dimensional, field-oriented system ($dq$-coordinates) are available to be used as the frame of input and output variables, allowing for easy and quick controller analysis and diagnose within the most convenient coordinate system. 
Finally, the torque $T$ resulting from the motor is applied to the __mechanical load__. 
The load is characterized by a moment of inertia and by a load torque $T_\mathrm{L}$ that is directed against the motor torque. 
Load torque behavior can be parameterized with respect to the angular velocity $\omega_\mathrm{me}$ in the form of constant, linear and quadratic dependency (and arbitrary combinations thereof). 
Speed changes that result from the difference between motor and load torque are modeled with another ODE which completely covers the mechanical system behavior.
Alternatively, the motor speed can be set to a fixed value, which can be useful for the investigation of control algorithms concerning generator operation, or it can be set to follow a specified trajectory, which is convenient when inspecting scenarios with defined speed demands like in traction applications. 

# Features

A large number of different motor systems is already implemented. 
These include DC drives as well as synchronous and induction three-phase drives. 
A complete list can be viewed in the ``GEM`` documentation [@GEM-docu].
The corresponding power converters allow to control the motor either directly via applied voltage (continuous-control-set) or by defining the converter's switching state (finite-control-set). 
Specifically for the use within reinforcement-learning applications and for testing state-of-the-art expert-driven control designs, the toolbox comes with a built-in reference generator, which can be used to create arbitrary reference trajectories (e.g., for the motor current, velocity or torque). 
These generated references are furthermore used to calculate a reward. In the domain of reinforcement learning, reward is the optimization variable that is to be maximized. 
For the control system scenario, reward is usually defined by the negative distance between the momentary and desired operation point, such that expedient controller behavior can be monitored easily.
The reward mechanism also allows to take physical limitations of the drive system into account, e.g., in the way of a notably low reward if limit values are surpassed. 
Optionally, the environment can be setup such that a reset of the system is induced in case of a limit violation. 
In addition, built-in visualization and plotting routines allow to monitor the training process of reinforcement learning agents or the performance of expert-driven control approaches.

# Examples

A minimal example of ``GEM's`` simulation capability is presented in \autoref{fig:SCIM_example}.
The plot shows the start-up behavior of a squirrel cage induction motor connected to an idealized three-phase electric grid depicting the angular velocity $\omega_\mathrm{me}$, the torque $T$, the voltage $u_{a,b,c}$ and the current $i_{d,q}$.
Here, the voltage is depicted within the physical $abc$-frame while the current is viewed within the simplified $dq$-frame. 

![Simulation of a squirrel cage induction motor connected to a rigid network at $50 \, \mathrm{Hz}$\label{fig:SCIM_example}](../plots/SCIM_Example.eps)

Exemplary code snippets that demonstrate the usage of ``GEM`` within both, the classical control and the reinforcement learning context are included within the project's [examples folder](https://github.com/upb-lea/gym-electric-motor/tree/master/examples). 
Featured examples:

- [``GEM_cookbook.ipynb``](https://colab.research.google.com/github/upb-lea/gym-electric-motor/blob/master/examples/environment_features/GEM_cookbook.ipynb): a basic tutorial-style notebook that presents the basic interface and usage of GEM
- [``scim_ideal_grid_simulation.py``](https://github.com/upb-lea/gym-electric-motor/blob/master/examples/environment_features/scim_ideal_grid_simulation.py): a simple motor simulation showcase of the squirrel cage induction motor that was used to create \autoref{fig:SCIM_example}



# References
.. gym-electric-motor documentation master file, created by
   sphinx-quickstart on Tue Jul  2 15:49:19 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to gym-electric-motor(GEM)'s documentation!
===================================================

The gym-electric-motor (GEM) package is a software toolbox for the simulation of different electric motors to
train and test reinforcement learning motor controllers and to compare them with classical motor controllers.


Getting started
***************************

A quick start guide can be found in the following Readme-File.

.. toctree::
   :maxdepth: 1
   :caption: Gym Electric Motor Readme:

   parts/readme


Content
*******

In the environments section all available GEM-environments are presented with their default configuration.
For quick start, one of these can be selected and used out of the box.


The documentation of the base classes is important for the development of own modules like further reward functions or
reference generators. In this part, the basic interfaces of each module are specified.
For the development of physical models like further motor models or further mechanical load models, the physical system
documentation specifies the basic interfaces inside a physical system.

..  toctree::
    :maxdepth: 4
    :titlesonly:
    :caption: gym-electric-motor Contents:

    parts/environments/environment
    parts/reference_generators/reference_generator
    parts/reward_functions/reward_function
    parts/physical_systems/physical_system
    parts/visualizations/visualization
    parts/constraint_monitor
    parts/core
    parts/utils
    parts/callbacks
    parts/random_component

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



Technical Models
################

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Introduction
************
The technical models contain the technical properties of the motor, the power electronic converter and the load.

The motors and the loads can be parametrized by either choosing the number of its default set in the environment
initialization or by passing a new Load or Motor parameter dictionary as described below.
The parameters of the power electronic converter are directly included in the initialization.



DC Motor Parameter Dictionary
'''''''''''''''''''''''''''''

+----------------------------+------------------------------------------------------------+
| **Key**                    |  **Description**                                           |
+============================+============================================================+
| r_a                        | Resistance of the armature circuit in Ohm                  |
+----------------------------+------------------------------------------------------------+
| r_e                        | Resistance of the excitation circuit in Ohm                |
+----------------------------+------------------------------------------------------------+
| l_a                        | Inductance of the armature circuit in Henry                |
+----------------------------+------------------------------------------------------------+
| l_e                        | Inductance of the excitation circuit in Henry              |
+----------------------------+------------------------------------------------------------+
| l_e_prime                  | effective excitation inductance in Henry                   |
+----------------------------+------------------------------------------------------------+
| psi_e                      | effective flux linkage in Henry Ampere                     |
+----------------------------+------------------------------------------------------------+
| j                          | Moment of inertia of the motors rotor in kg/m^2            |
+----------------------------+------------------------------------------------------------+
| u_sup                      | Supply voltage of the motor in Volt                        |
+----------------------------+------------------------------------------------------------+
| i_N                        | Nominal current in Ampere (only PermEx and Series)         |
+----------------------------+------------------------------------------------------------+
| i_a_N                      | Nominal armature current in Ampere(only ExtEx and Shunt)   |
+----------------------------+------------------------------------------------------------+
| i_e_N                      | Nominal excitation current in Ampere (only ExtEx and Shunt)|
+----------------------------+------------------------------------------------------------+
| u_N                        | Nominal voltage in Volt(only PermEx, Series and Shunt)     |
+----------------------------+------------------------------------------------------------+
| u_a_N                      | Nominal armature voltage in Volt (only ExtEx)              |
+----------------------------+------------------------------------------------------------+
| u_e_N                      | Nominal excitation voltage in Volt (only ExtEx)            |
+----------------------------+------------------------------------------------------------+
| torque_N                   | Nominal torque in Nm                                       |
+----------------------------+------------------------------------------------------------+
| omega_N                    | Nominal angular velocity rad/s                             |
+----------------------------+------------------------------------------------------------+


PMSM and SynRM Parameter Dictionary
'''''''''''''''''''''''''''''''''''

+------------------+------------------------------------------------+---------------------+
| **Key**          |  **Description**                               | **Default**         |
+==================+================================================+=====================+
| r_s              | Stator Resistance in Ohm                       | 4.9                 |
+------------------+------------------------------------------------+---------------------+
| l_d              | d-axis inductance in Henry                     | 79e-3               |
+------------------+------------------------------------------------+---------------------+
| l_q              | q-axis inductance in Henry                     | 113e-3              |
+------------------+------------------------------------------------+---------------------+
| j_rotor          | Moment of inertia of the rotor                 | 2.45e-3             |
+------------------+------------------------------------------------+---------------------+
| psi_p            | Permanent linked rotor flux                    | 0.165               |
+------------------+------------------------------------------------+---------------------+
| p                | Pole pair Number                               | 2                   |
+------------------+------------------------------------------------+---------------------+


All nominal voltages and currents are peak phase values.
Therefore, data sheet values for line voltages and phase currents has to be transformed such that
:math:`U_N=\sqrt(2/3) U_L` and :math:`I_N=\sqrt(2) I_S`.

Furthermore, the angular velocity is the electrical one and not the mechanical one :math:`\omega = p \omega_{me}`.

Load Parameter Dictionary
'''''''''''''''''''''''''

+----------------------------+----------------------------------------------------------+
| **Key**                    |  **Description**                                         |
+============================+==========================================================+
| a                          | Constant term in the load equation                       |
+----------------------------+----------------------------------------------------------+
| b                          | Linear factor in the load equation                       |
+----------------------------+----------------------------------------------------------+
| c                          | Quadratic factor in the load equation                    |
+----------------------------+----------------------------------------------------------+
| j_load                     | Moment of inertia of the load                            |
+----------------------------+----------------------------------------------------------+

DC Shunt Default Models
'''''''''''''''''''''''

+----------------------------+------------------------------+
| **Key**                    |  **No 0**                    |
+============================+==============================+
| r_a                        | 2.78                         |
+----------------------------+------------------------------+
| r_e                        | 350.0                        |
+----------------------------+------------------------------+
| l_a                        | 6.3e-3                       |
+----------------------------+------------------------------+
| l_e                        | 160.0                        |
+----------------------------+------------------------------+
| l_e_prime                  | 0.94                         |
+----------------------------+------------------------------+
| j_rotor                    | 0.017                        |
+----------------------------+------------------------------+



DC Series Default Models
''''''''''''''''''''''''

+----------------------------+------------------------------+
| **Key**                    |  **No 0**                    |
+============================+==============================+
| r_a                        | 2.78                         |
+----------------------------+------------------------------+
| r_e                        | 1.0                          |
+----------------------------+------------------------------+
| l_a                        | 6.3e-3                       |
+----------------------------+------------------------------+
| l_e                        | 1.6e-3                       |
+----------------------------+------------------------------+
| l_e_prime                  | 0.05                         |
+----------------------------+------------------------------+
| j_rotor                    | 0.017                        |
+----------------------------+------------------------------+


DC Permanently Excited Default Models
'''''''''''''''''''''''''''''''''''''

+----------------------------+------------------------------+
| **Key**                    |  **No 0**                    |
+============================+==============================+
| r_a                        | 25.0                         |
+----------------------------+------------------------------+
| r_e                        | 258.0                        |
+----------------------------+------------------------------+
| l_a                        | 3.438e-2                     |
+----------------------------+------------------------------+
| psi_e                      | 18                           |
+----------------------------+------------------------------+
| j_rotor                    | 0.0017                       |
+----------------------------+------------------------------+


DC Externally Excited Default Models
''''''''''''''''''''''''''''''''''''

+----------------------------+------------------------------+
| **Key**                    |  **No 0**                    |
+============================+==============================+
| r_a                        | 0.78                         |
+----------------------------+------------------------------+
| r_e                        | 350.0                        |
+----------------------------+------------------------------+
| l_a                        | 6.3e-3                       |
+----------------------------+------------------------------+
| l_e                        | 60                           |
+----------------------------+------------------------------+
| l_e_prime                  | 0.94                         |
+----------------------------+------------------------------+
| j                          | 0.017                        |
+----------------------------+------------------------------+
| u_sup                      | 420.0                        |
+----------------------------+------------------------------+
| i_a_N                      | 50.0                         |
+----------------------------+------------------------------+
| i_e_N                      | 1.2                          |
+----------------------------+------------------------------+
| U_a_N                      | 420.0                        |
+----------------------------+------------------------------+
| U_e_N                      | 420.0                        |
+----------------------------+------------------------------+
| torque_N                   | 40.0                         |
+----------------------------+------------------------------+
| omega_N                    | 368.0                        |
+----------------------------+------------------------------+
Callbacks
#####

.. automodule:: gym_electric_motor.callbacks
    :members:
Utils
#####

.. automodule:: gym_electric_motor.utils
    :members:

Constraint Monitor
#####################

..  toctree::
    :maxdepth: 1
    :caption: Available Constraints:

    constraints/constraint
    constraints/limit_constraint
    constraints/squared_constraint

Usage Guide
___________

ToDo

Constraint Monitor API Documentation
____________________________________
.. autoclass:: gym_electric_motor.core.ConstraintMonitor
    :members:
    :inherited-members:
Random Component
############################

.. autoclass:: gym_electric_motor.RandomComponent
   :members:
   :inherited-members:

Technical Background
====================
The package contains different electric motors, power electronic converters and load models.
In this part the technical models are described and references to more detailed explanations are given.
All included electric motors can be represented by a system of differential equations and those are the same for the
continuous and discrete action case. The load model is the same in all cases.
The converters are shortly introduced. More detailed descriptions are included in each converter class, where also the
action space is specified.

The DC motors

- permanently excited
- externally excited
- series
- shunt

and the synchronous motors

- permanent magnet synchronous motor (PMSM)
- synchronous reluctance motor (SynRM)

are included. The figure below shows the basic structure of converter, motor and load that is modeled.

.. figure:: ../plots/FigureConvMotorLoad6.svg

Further information about electrical motor and power electronic converters can be found in [Boecker2018a]_, [Boecker2018b]_ and [Chiasson2005]_.

DC Motors
#########

Externally Excited Motor
------------------------
The figure shows the circuit diagram of this motor [Boecker2018a]_.

.. figure:: ../plots/ESBdcExtEx.svg

The equations are

:math:`u_A=\mathit{\Psi}^\prime_E \omega + L_A \frac{\mathrm{d} i_A}{\mathrm{d} t} +R_A i_A`

:math:`u_E=L_E \frac{\mathrm{d} i_E}{\mathrm{d} t} + R_E i_E`

:math:`\mathit{\Psi}^\prime_E=L^\prime_E i_E`

:math:`T=\mathit{\Psi}^\prime_E i_A`

:math:`\frac{\mathrm{d} \omega_{me}}{\mathrm{d} t}=\frac{T-T_L(\omega_{me})}{J}`

and can be rewritten as

:math:`\frac{\mathrm{d} i_A}{\mathrm{d} t}=\frac{u_A-L^\prime_E \omega i_E -R_A i_A}{L_A}`

:math:`\frac{\mathrm{d} i_E}{\mathrm{d} t}=\frac{u_E-R_E i_E}{L_E}`

:math:`\frac{\mathrm{d} \omega_{me}}{\mathrm{d} t}=\frac{L^\prime_E i_A i_E -T_L(\omega_{me})}{J}\text{.}`

For DC motors is :math:`\omega_{me}=\omega` valid.

The quantities for this and the other motors are:

- :math:`u_A` armature voltage

- :math:`u_E` excitation voltage

- :math:`i_A` armature current

- :math:`i_E` armature current

- :math:`R_A` armature resistance

- :math:`R_E` excitation resistance

- :math:`L_A` armature inductance

- :math:`L_E` excitation inductance

- :math:`\omega` (electrical) angular velocity

- :math:`\omega_{me}` mechanical angular velocity

- :math:`L^\prime_E` effective excitation inductance

- :math:`T` Torque produced by the motor

- :math:`T_L` Torque from the load

- :math:`J` moment of inertia

- :math:`\mathit{\Psi}^\prime_E` effective excitation flux

Other motors are build on this motor with different connections of armature and excitation circuit.



Series Motor
------------
In this type both circuits are in series connection, :math:`i=i_A=i_E` and :math:`u=u_A+u_E`,  as shown in the figure [Boecker2018a]_.

.. figure:: ../plots/ESBseries.svg

The equation for this motor are

:math:`\frac{\mathrm{d} i}{\mathrm{d} t}=-\frac{L^\prime_E}{L_A+L_E} i \omega -\frac{R_A+R_E}{L_A+L_E} i + \frac{1}{L_A+L_E} u`

:math:`\frac{\mathrm{d} \omega}{\mathrm{d} t}=\frac{L^\prime_E}{J} i^2-\frac{1}{J}T_L(\omega)`

The torque equation is the same as before.


Shunt Motor
-----------
In this type the circuits are connected in parallel, :math:`i=i_A+i_E` and :math:`u=u_A=u_E`,  as shown in the figure [Boecker2018a]_.

.. figure:: ../plots/ESBshunt.svg

The equation for this motor are

:math:`\frac{\mathrm{d} i_A}{\mathrm{d} t}=\frac{u-L^\prime_E \omega i_E -R_A i_A}{L_A}`

:math:`\frac{\mathrm{d} i_E}{\mathrm{d} t}=\frac{u-R_E i_E}{L_E}`

:math:`\frac{\mathrm{d} \omega}{\mathrm{d} t}=\frac{L^\prime_E i_A i_E -T_L(\omega)}{J}`


and the torque equation is the same as before.

Permanently Excited Motor
--------------------------
In this type the excitation :math:`\mathit{\Psi}^\prime_E` is constant and therefore it only contains the armature circuit of the
externally excited motor.

The equations are

:math:`\frac{\mathrm{d} i}{\mathrm{d} t}=\frac{u-\mathit{\Psi}^\prime_E \omega -R_A i}{L_A}`

:math:`\frac{\mathrm{d} \omega}{\mathrm{d} t}=\frac{L^\prime_E i^2 -T_L(\omega)}{J}`


Permanent Magnet Synchronous Motor (PMSM)
#########################################

The PMSM is a three phase motor with a permanent magnet in the rotor as shown in the figure [Boecker2018b]_. The input of this motor are
the voltages :math:`u_a`, :math:`u_b` and :math:`u_c`.

The quantities are:

- :math:`u_a`, :math:`u_b`, :math:`u_c` phase voltages

- :math:`i_a`, :math:`i_b`, :math:`i_c` phase currents

- :math:`R_s` stator resistance

- :math:`L_d` d-axis inductance

- :math:`L_q` q-axis inductance


- :math:`i_{sd}` d-axis current

- :math:`i_{sq}` q-axis current

- :math:`u_{sd}` d-axis voltage

- :math:`u_{sq}` q-axis voltage

- :math:`p` pole pair number

- :math:`\mathit{\Psi}_p` permanent linked rotor flux

- :math:`\epsilon` rotor position angle

- :math:`\omega` (electrical) angular velocity

- :math:`\omega_{me}` mechanical angular velocity

- :math:`T` Torque produced by the motor

- :math:`T_L` Torque from the load

- :math:`J` moment of inertia

The electrical angular velocity and the mechanical angular velocity are related such that :math:`\omega=\omega_{me} p`.

.. figure:: ../plots/GDAFig29.svg

The circuit diagram of the phases are similar to each other and the armature circuit of the externally excited motor.

.. figure:: ../plots/pmsmMotorB6.png

For an easy computation the three phases are first transformed to the quantities :math:`\alpha` and :math:`\beta` and
afterwards to :math:`d/q` coordinates that rotated with the rotor as given in [Boecker2018b]_.

.. figure:: ../plots/ESBdq.svg

This results in the equations:


:math:`u_{sd}=R_s i_{sd}+L_d \frac{\mathrm{d} i_{sd}}{\mathrm{d} t}-\omega_{me}p L_q i_{sq}`

:math:`u_{sq}=R_s i_{sq}+L_q \frac{\mathrm{d} i_{sq}}{\mathrm{d} t}+\omega_{me}p L_d i_{sd}+\omega_{me}p \mathit{\Psi}_p`

:math:`\frac{\mathrm{d} \omega_{me}}{\mathrm{d} t}=\frac{T-T_L(\omega_{me})}{J}`

:math:`T=\frac{3}{2} p (\mathit{\Psi}_p +(L_d-L_q)i_{sd}) i_{sq}`



A more detailed derivation can be found in
[Modeling and High-Performance Control of Electric Machines, John Chiasson (2005)]

The difference between rms and peak values and between line and phase quantities has to be considered at the PMSM.
The PMSM is in star conncetion and the line voltage :math:`U_L` is mostly given in data sheets as rms value.
In the toolbox the nominal value of the phase voltage :math:`\hat{U}_S=\sqrt{\frac{2}{3}}U_L` is needed.
Furthermore, the supply voltage is typically the same :math:`u_{sup}=\hat{U}_S`.
For example, a line voltage of :math:`U_L=400~\text{V}` is given, the rms phase voltage is :math:`U_S=\sqrt{\frac{1}{3}}U_L = 230.9 \text{ V}`
and the peak value :math:`\hat{U}_S=326.6 \text{ V}`.
The nominal peak current of a phase is given by :math:`\hat{I}_S=\sqrt{2} I_S`.

.. figure:: ../plots/Drehstromtrafo.svg

Converter
#########

The DC input voltage of the converter is named supply voltage :math:`u_{sup}` in the package and the output voltage and
current are the input quantities of the motor and therefore are named as :math:`u_{in}` and :math:`i_{in}`.
All converter contain a dead time of one sampling interval and an interlocking times can be considered.
More details can be found in [Boecker2018a]_, [Boecker2018b]_ and [Chiasson2005]_.

1 Quadrant Converter (1QC)
--------------------------
This converter can provide positive voltages and positive currents at the output [Boecker2018a]_.

:math:`u_{in} \geq 0`

:math:`i_{in} \geq 0`

.. figure:: ../plots/1QC.svg

2 Quadrant Converter (2QC)
--------------------------

This converter can provide positive voltages and both current directions at the output [Boecker2018a]_.

:math:`u_{in} \geq 0`

:math:`i_{in}` both signs

.. figure:: ../plots/2QCigbt.svg

4 Quadrant Converter (4QC)
--------------------------

This converter can provide both voltage and currents polarities at the output [Boecker2018a]_.

:math:`u_{in}` both signs

:math:`i_{in}` both signs

.. figure:: ../plots/4QC.svg

Three-Phase Voltage Source Inverter
-----------------------------------

This converter is also called B6 bridge and is used for three phase motors as the PMSM. The supply voltage is DC and it
consists of three parallel half bridges, one for each phase :math:`u_a`, :math:`u_b` and :math:`u_c`. Due to this 8
different switching states are possible. A dead time of one sampling time step and
interlocking times can be considered. A symmetric B6 bridge with a voltage range of :math:`[-u_{DC}/2,~+u_{DC}/2]` is used.
[Boecker2018b]_

.. figure:: ../plots/B6.svg




Load Models
###########

The package can deal with polynomial load equations

:math:`T_L(\omega_{me})=\mathrm{sign}(\omega_{me})(c \omega^2_{me} + b \vert\omega_{me}\vert + a)`.

Furthermore an additional moment of inertia :math:`J_{Load}` of the load can be set.


References
##########

.. [Boecker2018a] Böcker, Joachim; Elektrische Antriebstechnik; 2018; Paderborn University

.. [Boecker2018b] Böcker, Joachim; Controlled Three-Phase Drives; 2018; Paderborn University

.. [Chiasson2005] Chiasson, John; Modeling and High-Performance Control of Electric Machines; 2005; Hoboken, NJ, USA




-----------
Readme File
-----------

.. mdinclude:: ../../README.mdCore
#####

.. figure:: ../plots/CoreClasses.svg

.. automodule:: gym_electric_motor.core
    :members:
Weighted Sum of Errors
######################

Usage Guide
***********

To use the weighted sum of errors, you have to import the class, initialize an object and pass it to the environment.


.. code-block:: python

    import gym_electric_motor as gem
    from gym_electric_motor.reward_functions import WeightedSumOfErrors


    # initialize the reward function
    wse = WeightedSumOfErrors(
        reward_weights=dict(i_a=1, i_e=2) # Current control problem. Tracking of i_e is rewarded better.
        reward_power=2 # Squared Error
        # Alternative: reward_power=dict(i_a=1, i_e=0.5) Absolute error i_a, root error on i_e
        bias='positive' # Shift the reward range from negative to positive
        violation_reward=-250 # Self defined violation reward
        gamma=0.9 # Ignored, if a violation_reward is defined.
        normed_reward_weights=False # Otherwise weights will be normed automatically to sum up to 1.

    )

    # pass it to the environment
    env = gem.make('my-env-id-v0', reward_function=wse)

API Documentation
*****************

.. autoclass:: gym_electric_motor.reward_functions.weighted_sum_of_errors.WeightedSumOfErrors
   :members:
   :inherited-members:
Reward Functions
#####################

..  toctree::
    :maxdepth: 1
    :caption: Available Reward Functions:

    weighted_sum_of_errors

Reward Function Base Class
''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.core.RewardFunction
   :members:
Console Printer
###############################

.. autoclass:: gym_electric_motor.visualization.console_printer.ConsolePrinter
   :members:
   :inherited-members:
Motor Dashboard
###############################

..  contents::


..  toctree::
    :maxdepth: 1
    :caption: Available Plots:
    :glob:

    motor_dashboard_plots/*

Usage Guide
__________________
To use the dashboard, you have to import the class, define the plots,
instantiate the dashboard and pass it to the environment.

The most common plots can be quickly selected directly in the constructor of the dashboard.
Further plots like a *MeanEpisodeRewardPlot* or self-defined ones have to be instantiated and passed to the dashboard.

.. code-block:: python

    import gym_electric_motor as gem
    from gym_electric_motor.visualization import MotorDashboard
    from gym_electric_motor.visualization.motor_dashboard_plots import MeanEpisodeRewardPlot

    # create the dashboard and define the plots
    dashboard = MotorDashboard(
        state_plots = ['omega', 'i'], # Pass a list of state names or 'all' for all states to plot
        reward_plot = True, # True / False (False default)
        action_plots = [0] # 'all' plots all actions (if multiple are applied)
        additional_plots=[MeanEpisodeRewardPlot()] # Add all further plots here
    )

    # pass it to the environment
    env = gem.make('my-env-id-v0', visualization=dashboard)




Motor Dashboard API
__________________________

.. autoclass:: gym_electric_motor.visualization.motor_dashboard.MotorDashboard
   :members:
   :inherited-members:

..
    The following section is commented out. It may serve as a little introduction on how to define
    your own custom plots in the future.

    Create your own plots
    _____________________________
     In the following, there is
    .. code-block:: python

        import gym_electric_motor as gem
        from gym_electric_motor.visualization import MotorDashboard
        from gym_electric_motor.visualization.motor_dashboard_plots import TimePlot

        # the class may also derive from EpisodePlot or StepPlot, depending on the x-axis
        class MyPlot(TimePlot):
        """This plot will show the current step of the episode *k* on the y-axis.

         As it derives from TimePlot the x-Axis is the cumulative simulated time over all episodes.
         """

            def __init__(self):
                super().__init__()
                # Set the y-axis label
                self._label = 'k'

            def initialize(self, axis):
                super().initialize(axis)
                self._k_line, =  self._axis.plot([],[])
                self._lines.append(self._k_line)

            def set_env(self, env):
                super().set_env(env)
                self._k_data = np.ones_like(self._x_data) * np.nan
                self._y_data.append(self._k_data)

            def on_step_begin(self, k, action):
                super().on_step_begin(k, action)
                idx = self.data_idx
                self._k_data[idx] = k

            def _scale_y_axis(self):
                """This function can be defined to automatically scale the plots limits.
                Here, the y limits are set such that the data fits perfectly.
                """
                self._axis.set_ylim(0,max(self._k_data))

        dashboard = MotorDashboard(
            state_plots='all',
            additional_plots=[MyPlot()]
        )
        env = gem.make('my-env-id-v0', visualization=dashboard)
Visualization
#############

..  toctree::
    :maxdepth: 1
    :caption: Available Visualizations:

    motor_dashboard
    console_printer


Visualization Base Class
''''''''''''''''''''''''''''''

..  autoclass:: gym_electric_motor.core.ElectricMotorVisualization
    :members:
    :inherited-members:
Reward Plot
###########

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.RewardPlot
   :members:
   :inherited-members:
Episode Plot (Abstract)
#####################
The EpisodePlot is the base class for all plots that plot data with number of episodes on the x-axis.

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.EpisodePlot
   :members:
   :inherited-members:
Cumulative Constraint Violation Plot
####################################

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.CumulativeConstraintViolationPlot
   :members:
   :inherited-members:Time Plot (Abstract)
#####################
The TimePlot is the base class for all plots that plot data with the time on the x-axis.

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.TimePlot
   :members:
   :inherited-members:
Mean Episode Reward Plot
########################

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.MeanEpisodeRewardPlot
   :members:
   :inherited-members:Episode Length Plot
####################################

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.EpisodeLengthPlot
   :members:
   :inherited-members:Action Plot
###########

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.ActionPlot
   :members:
   :inherited-members:
Step Plot (Abstract)
#####################
The StepPlot is the base class for all plots that plot data with the cumulative number of steps on the x-axis.

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.TimePlot
   :members:
   :inherited-members:
State Plot
###########

.. autoclass:: gym_electric_motor.visualization.motor_dashboard_plots.StatePlot
   :members:
   :inherited-members:
Electric Motors
###############

Electric Motor Base Class
*************************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.ElectricMotor
   :members:


Synchronous Motors
******************

Parameter Dictionary
''''''''''''''''''''

+------------------+------------------------------------------------+---------------------+
| **Key**          |  **Description**                               | **Default**         |
+==================+================================================+=====================+
| r_s              | Stator Resistance in Ohm                       | 4.9                 |
+------------------+------------------------------------------------+---------------------+
| l_d              | d-axis inductance in Henry                     | 79e-3               |
+------------------+------------------------------------------------+---------------------+
| l_q              | q-axis inductance in Henry                     | 113e-3              |
+------------------+------------------------------------------------+---------------------+
| j_rotor          | Moment of inertia of the rotor                 | 2.45e-3             |
+------------------+------------------------------------------------+---------------------+
| psi_p            | Permanent linked rotor flux                    | 0.165               |
+------------------+------------------------------------------------+---------------------+
| p                | Pole pair Number                               | 2                   |
+------------------+------------------------------------------------+---------------------+


All nominal voltages and currents are peak phase values.
Therefore, data sheet values for line voltages and phase currents has to be transformed such that
:math:`U_N=\sqrt(2/3) U_L` and :math:`I_N=\sqrt(2) I_S`.

Furthermore, the angular velocity is the electrical one and not the mechanical one :math:`\omega = p \omega_{me}`.


.. autoclass:: gym_electric_motor.physical_systems.electric_motors.SynchronousMotor
   :members:

Synchronous Reluctance Motor
****************************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.SynchronousReluctanceMotor
   :members:

Permanent Magnet Synchronous Motor
**********************************

The PMSM is a three phase motor with a permanent magnet in the rotor as shown in the figure [Boecker2018b]_.
The input of this motor are the voltages :math:`u_a`, :math:`u_b` and :math:`u_c`.

The quantities are:

- :math:`u_a`, :math:`u_b`, :math:`u_c` phase voltages

- :math:`i_a`, :math:`i_b`, :math:`i_c` phase currents

- :math:`R_s` stator resistance

- :math:`L_d` d-axis inductance

- :math:`L_q` q-axis inductance


- :math:`i_{sd}` d-axis current

- :math:`i_{sq}` q-axis current

- :math:`u_{sd}` d-axis voltage

- :math:`u_{sq}` q-axis voltage

- :math:`p` pole pair number

- :math:`\mathit{\Psi}_p` permanent linked rotor flux

- :math:`\epsilon` rotor position angle

- :math:`\omega` (electrical) angular velocity

- :math:`\omega_{me}` mechanical angular velocity

- :math:`T` Torque produced by the motor

- :math:`T_L` Torque from the load

- :math:`J` moment of inertia

The electrical angular velocity and the mechanical angular velocity are related such that :math:`\omega=\omega_{me} p`.

.. figure:: ../../plots/GDAFig29.svg

The circuit diagram of the phases are similar to each other and the armature circuit of the externally excited motor.

.. figure:: ../../plots/pmsmMotorB6.png

For an easy computation the three phases are first transformed to the quantities :math:`\alpha` and :math:`\beta` and
afterwards to :math:`d/q` coordinates that rotated with the rotor as given in [Boecker2018b]_.

.. figure:: ../../plots/ESBdq.svg

This results in the equations:


:math:`u_{sd}=R_s i_{sd}+L_d \frac{\mathrm{d} i_{sd}}{\mathrm{d} t}-\omega_{me}p L_q i_{sq}`

:math:`u_{sq}=R_s i_{sq}+L_q \frac{\mathrm{d} i_{sq}}{\mathrm{d} t}+\omega_{me}p L_d i_{sd}+\omega_{me}p \mathit{\Psi}_p`

:math:`\frac{\mathrm{d} \omega_{me}}{\mathrm{d} t}=\frac{T-T_L(\omega_{me})}{J}`

:math:`T=\frac{3}{2} p (\mathit{\Psi}_p +(L_d-L_q)i_{sd}) i_{sq}`



A more detailed derivation can be found in
[Modeling and High-Performance Control of Electric Machines, John Chiasson (2005)]

The difference between rms and peak values and between line and phase quantities has to be considered at the PMSM.
The PMSM is in star conncetion and the line voltage :math:`U_L` is mostly given in data sheets as rms value.
In the toolbox the nominal value of the phase voltage :math:`\hat{U}_S=\sqrt{\frac{2}{3}}U_L` is needed.
Furthermore, the supply voltage is typically the same :math:`u_{sup}=\hat{U}_S`.
For example, a line voltage of :math:`U_L=400~\text{V}` is given, the rms phase voltage is
:math:`U_S=\sqrt{\frac{1}{3}}U_L = 230.9 \text{ V}`
and the peak value :math:`\hat{U}_S=326.6 \text{ V}`.
The nominal peak current of a phase is given by :math:`\hat{I}_S=\sqrt{2} I_S`.

.. figure:: ../../plots/Drehstromtrafo.svg


.. autoclass:: gym_electric_motor.physical_systems.electric_motors.PermanentMagnetSynchronousMotor
   :members:



References
##########

.. [Boecker2018a] Böcker, Joachim; Elektrische Antriebstechnik; 2018; Paderborn University

.. [Boecker2018b] Böcker, Joachim; Controlled Three-Phase Drives; 2018; Paderborn University

.. [Chiasson2005] Chiasson, John; Modeling and High-Performance Control of Electric Machines; 2005; Hoboken, NJ, USA



Supply Converter Motor Load System (SCML)
#########################################

The technical structure of the SCML-Systems in the GEM-toolbox can bee seen in the figure below.

.. figure:: ../../plots/SCML_Setting.svg

The system consists of a Voltage Supply, a Power Electronic Converter, an Electrical Motor and the Mechanical Load.
Additionally, each SCML-System has got an ODE-Solver for the simulation.

..  toctree::
    :maxdepth: 1
    :caption: Subcomponents of the SCML:

    converters/converter
    mechanical_loads/mechanical_load
    noise_generators/noise_generator
    electric_motors/electric_motor
    voltage_supplies/voltage_supply
    ode_solvers/ode_solver



The abstract SCML-System defines the overall structure. From this, the DcMotorSystem and the SynchronousMotorSystem
derive. They only implement private methods. Therefore, the interface to the user stays the same in all cases.

.. autoclass:: gym_electric_motor.physical_systems.physical_systems.SCMLSystem
   :members:

Dc Motor System
****************
.. autoclass:: gym_electric_motor.physical_systems.physical_systems.DcMotorSystem
   :members:

Synchronous Motor System
************************
.. autoclass:: gym_electric_motor.physical_systems.physical_systems.SynchronousMotorSystem
   :members:

Squirrel Cage Induction Motor System
************************
.. autoclass:: gym_electric_motor.physical_systems.physical_systems.SquirrelCageInductionMotorSystem
   :members:

Doubly Fed Induction Motor System
************************
.. autoclass:: gym_electric_motor.physical_systems.physical_systems.DoublyFedInductionMotorSystem
   :members:Physical Systems
==========================

..  toctree::
    :maxdepth: 1
    :caption: Available Physical Systems:

    scml_system


Physical System Base Class
''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.core.PhysicalSystem
   :members:
AC 1 Phase Supply
##############################

.. autoclass:: gym_electric_motor.physical_systems.voltage_supplies.AC1PhaseSupply
    :members:
    :inherited-members:
Ideal Voltage Supply
##############################

.. autoclass:: gym_electric_motor.physical_systems.voltage_supplies.IdealVoltageSupply
    :members:
    :inherited-members:
Voltage Supplies
################

..  toctree::
    :maxdepth: 1
    :caption: Available Voltage Supplies:

    ideal_voltage_supply
    rc_voltage_supply
    ac_1_phase_supply
    ac_3_phase_supply


Voltage Supply Base Class
*************************************

.. autoclass:: gym_electric_motor.physical_systems.voltage_supplies.VoltageSupply
   :members:
RC Voltage Supply
##############################

.. autoclass:: gym_electric_motor.physical_systems.voltage_supplies.RCVoltageSupply
    :members:
    :inherited-members:
AC 3 Phase Supply
##############################

.. autoclass:: gym_electric_motor.physical_systems.voltage_supplies.AC3PhaseSupply
    :members:
    :inherited-members:
Base Three Phase Motor
############################



Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.ThreePhaseMotor
   :members:
   :inherited-members:
Base Synchronous Motor
############################



Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.SynchronousMotor
   :members:
   :inherited-members:
DC Series Motor
############################

Schematic
*********

.. figure:: ../../../plots/ESBseries.svg

Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i}{\mathrm{d} t} &= \frac{u - L_\mathrm{e}^\prime i \omega_\mathrm{me} - (R_\mathrm{a} + R_\mathrm{e}) i}{L_\mathrm{a} + L_\mathrm{e}} \\


Torque Equation
***************
.. math::
    T = L_\mathrm{e}^\prime i^2

Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.DcSeriesMotor
   :members:
   :inherited-members:
Permanently Excited DC Motor
#############################

Schematic
*********

.. figure:: ../../../plots/ESBdcpermex.svg


Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i}{\mathrm{d} t} &= \frac{u_\mathrm{a} - \mathit{\Psi}^\prime_\mathrm{e} \omega_\mathrm{me} - R_\mathrm{a} i}{L_\mathrm{a}} \\


Torque Equation
***************
.. math::
    T = \mathit{\Psi}_\mathrm{e} i



Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.DcPermanentlyExcitedMotor
   :members:
   :inherited-members:
Squirrel Cage Induction Motor
##################################

Schematic
*********

.. figure:: ../../../plots/ESB_SCIM_alphabeta.svg

Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i_\mathrm{s \alpha}}{\mathrm{d} t}&= -\frac{1}{\tau_\sigma} i_\mathrm{s \alpha} + \frac{R_\mathrm{r} L_\mathrm{m}}{\sigma L_\mathrm{r}^2 L_\mathrm{s}} \psi_\mathrm{r \alpha} + p \omega_{\text{me}} \frac{L_\mathrm{m}}{\sigma L_\mathrm{r} L_\mathrm{s}} \psi_\mathrm{r \beta} + \frac{1}{\sigma L_\mathrm{s}} u_\mathrm{s \alpha}\\
    \frac{\mathrm{d} i_\mathrm{s \beta}}{\mathrm{d} t}&= -\frac{1}{\tau_\sigma} i_\mathrm{s \beta} - p \omega_{\text{me}} \frac{L_\mathrm{m}}{\sigma L_\mathrm{r} L_\mathrm{s}}  \psi_\mathrm{r \alpha}  + \frac{R_\mathrm{r} L_\mathrm{m}}{\sigma L_\mathrm{r}^2 L_\mathrm{s}} \psi_\mathrm{r \beta} + \frac{1}{\sigma L_\mathrm{s}} u_\mathrm{s \beta}\\
    \frac{\mathrm{d} \psi_\mathrm{r \alpha}}{\mathrm{d} t}&= \frac{L_\mathrm{m}}{\tau_\mathrm{r}} i_\mathrm{s \alpha} - \frac{1}{\tau_\text{r}} \psi_\mathrm{r \alpha} - p \omega_{\text{me}} \psi_\mathrm{r \beta} \\
    \frac{\mathrm{d} \psi_\mathrm{r \beta}}{\mathrm{d} t}&= \frac{L_\mathrm{m}}{\tau_\mathrm{r}} i_\mathrm{s \beta} + p \omega_{\mathrm{me}} \psi_\mathrm{r \alpha} - \frac{1}{\tau_\mathrm{r}} \psi_\mathrm{r \beta} \\
    \frac{\mathrm{d} \varepsilon_\mathrm{el}}{\mathrm{d} t}&= p \omega_\mathrm{me}

with

.. math::
    L_\mathrm{s} &= L_\mathrm{m} + L_\mathrm{\sigma s}, & L_\mathrm{r} &= L_\mathrm{m} + L_\mathrm{\sigma r}\\
    \sigma &= \frac{L_\mathrm{r} L_\mathrm{s} - L_\mathrm{m}^2}{L_\mathrm{r} L_\mathrm{s}}, & \tau_\mathrm{r} &=\frac{L_\mathrm{r}}{R_\mathrm{r}}, & \tau_\sigma &= \frac{\sigma L_\mathrm{s}}{R_\mathrm{s} + R_\mathrm{r} \frac{L_\mathrm{m}^2}{L_\mathrm{r}^2}}


Torque Equation
***************

.. math:: T=\frac{3}{2} p \frac{L_\mathrm{m}}{L_\mathrm{r}} (\psi_\mathrm{r \alpha} i_\mathrm{s \beta} - \psi_\mathrm{r \beta} i_\mathrm{s \alpha})

Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.SquirrelCageInductionMotor
   :members:
   :inherited-members:Doubly Fed Induction Motor
##################################

Schematic
*********

.. figure:: ../../../plots/ESB_DFIM_alphabeta.svg

Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i_\mathrm{s \alpha}}{\mathrm{d} t}&= -\frac{1}{\tau_\sigma} i_\mathrm{s \alpha} + \frac{R_\mathrm{r} L_\mathrm{m}}{\sigma L_\mathrm{r}^2 L_\mathrm{s}} \psi_\mathrm{r \alpha} + p \omega_\mathrm{\text{me}} \frac{L_\mathrm{m}}{\sigma L_\mathrm{r} L_\mathrm{s}}  \psi_\mathrm{r \beta} + \frac{1}{\sigma L_\mathrm{s}} u_\mathrm{s \alpha} - \frac{L_\mathrm{m}}{\sigma L_\mathrm{r} L_\mathrm{s}} u_\mathrm{r \alpha}\\
    \frac{\mathrm{d} i_\mathrm{s \beta}}{\mathrm{d} t}&= -\frac{1}{\tau_\sigma} i_\mathrm{s \beta} - p \omega_\mathrm{\text{me}} \frac{L_\mathrm{m}}{\sigma L_\mathrm{r} L_\mathrm{s}}  \psi_\mathrm{r \alpha}  + \frac{R_\mathrm{r} L_\mathrm{m}}{\sigma L_\mathrm{r}^2 L_\mathrm{s}} \psi_\mathrm{r \beta} + \frac{1}{\sigma L_\mathrm{s}} u_\mathrm{s \beta} - \frac{L_\mathrm{m}}{\sigma L_\mathrm{r} L_\mathrm{s}} u_\mathrm{r \beta}\\
    \frac{\mathrm{d} \psi_\mathrm{r \alpha}}{\mathrm{d} t}&= \frac{L_\mathrm{m}}{\tau_\mathrm{r}} i_\mathrm{s \alpha} - \frac{1}{\tau_\mathrm{r}} \psi_\mathrm{r \alpha} - p \omega_\mathrm{\text{me}} \psi_\mathrm{r \beta} + u_\mathrm{r \alpha}\\
    \frac{\mathrm{d} \psi_\mathrm{r \beta}}{\mathrm{d} t}&= \frac{L_\mathrm{m}}{\tau_\mathrm{r}} i_\mathrm{s \beta} + p \omega_\mathrm{\text{me}} \psi_\mathrm{r \alpha} - \frac{1}{\tau_\mathrm{r}} \psi_\mathrm{r \beta} + u_\mathrm{r \beta}\\
    \frac{\mathrm{d} \varepsilon_\mathrm{el}}{\mathrm{d} t}&= p \omega_\mathrm{me}

with

.. math::
    L_\mathrm{s} &= L_\mathrm{m} + L_\mathrm{\sigma s} & \quad L_\mathrm{r} &= L_\mathrm{m} + L_\mathrm{\sigma r}\\
    \sigma &= \frac{L_\mathrm{r} L_\mathrm{s} - L_\mathrm{m}^2}{L_\mathrm{r} L_\mathrm{s}} & \quad \tau_\mathrm{r} &=\frac{L_\mathrm{r}}{R_\mathrm{r}} & \quad \tau_\sigma &= \frac{\sigma L_\mathrm{s}}{R_\mathrm{s} + R_\mathrm{r} \frac{L_\mathrm{m}^2}{L_\mathrm{r}^2}}


Torque Equation
***************

.. math:: T=\frac{3}{2} p \frac{L_\mathrm{m}}{L_\mathrm{r}} (\psi_\mathrm{r \alpha} i_\mathrm{s \beta} - \psi_\mathrm{r \beta} i_\mathrm{s \alpha})

Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.DoublyFedInductionMotor
   :members:
   :inherited-members:Synchronous Reluctance Motor
############################

Schematic
*********

.. figure:: ../../../plots/ESBdqSynRM.svg

Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i_\mathrm{sd}}{\mathrm{d} t}&=\frac{u_\mathrm{sd} + p \omega_\mathrm{me} L_\mathrm{q} i_\mathrm{sq} - R_\mathrm{s} i_\mathrm{sd}}{L_\mathrm{d}} \\
    \frac{\mathrm{d} i_\mathrm{sq}}{\mathrm{d} t}&=\frac{u_\mathrm{sq} - p \omega_\mathrm{me} L_\mathrm{d} i_\mathrm{sd} - R_\mathrm{s} i_\mathrm{sq}}{L_\mathrm{q}} \\
    \frac{\mathrm{d} \varepsilon_\mathrm{el}}{\mathrm{d} t}&= p \omega_\mathrm{me}

Torque Equation
***************
.. math::
    T = \frac{3}{2} p (L_\mathrm{d} - L_\mathrm{q}) i_\mathrm{sd} i_\mathrm{sq}


Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.SynchronousReluctanceMotor
   :members:
   :inherited-members:
Externally Excited DC Motor
############################

Schematic
*********

.. figure:: ../../../plots/ESBdcExtEx.svg

Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i_\mathrm{a}}{\mathrm{d} t} &= \frac{u_\mathrm{a} - L_\mathrm{e}^\prime i_\mathrm{e} \omega_\mathrm{me} - R_\mathrm{a} i_\mathrm{a}}{L_\mathrm{a}} \\
    \frac{\mathrm{d} i_\mathrm{e}}{\mathrm{d} t} &= \frac{u_\mathrm{e} - R_\mathrm{e} i_\mathrm{e}}{L_\mathrm{e}}


Torque Equation
***************
.. math::
    T = L_\mathrm{e}^\prime i_\mathrm{e} i_\mathrm{a}


Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.DcExternallyExcitedMotor
   :members:
   :inherited-members:
Permanent Magnet Synchronous Motor
##################################

Schematic
*********

.. figure:: ../../../plots/ESBdq.svg

Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i_\mathrm{sd}}{\mathrm{d} t}&=\frac{u_\mathrm{sd} + p \omega_\mathrm{me} L_\mathrm{q} i_\mathrm{sq} - R_\mathrm{s} i_\mathrm{sd}}{L_\mathrm{d}} \\
    \frac{\mathrm{d} i_\mathrm{sq}}{\mathrm{d} t}&=\frac{u_\mathrm{sq} - p \omega_\mathrm{me} (L_\mathrm{d} i_\mathrm{sd} + \mathit{\Psi}_\mathrm{p}) - R_\mathrm{s} i_\mathrm{sq}}{L_\mathrm{q}} \\
    \frac{\mathrm{d} \varepsilon_\mathrm{el}}{\mathrm{d} t}&= p \omega_\mathrm{me}



Torque Equation
***************

.. math:: T=\frac{3}{2} p (\mathit{\Psi}_\mathrm{p} +(L_\mathrm{d}-L_\mathrm{q})i_\mathrm{sd}) i_\mathrm{sq}

Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.PermanentMagnetSynchronousMotor
   :members:
   :inherited-members:
Electric Motors
######################################

..  toctree::
    :maxdepth: 1
    :caption: Available Electric Motors:

    permex
    extex
    series
    shunt
    pmsm
    synrm
    scim
    dfim
    dc_base
    three_phase_base
    synchronous_base
    induction_base

Electric Motor Base Class
*************************************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.ElectricMotor
    :members:Base Induction Motor
############################



Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.InductionMotor
   :members:
   :inherited-members:
DC Shunt Motor
############################

Schematic
*********
.. figure:: ../../../plots/ESBshunt.svg

Electrical ODE
**************

.. math::
    \frac{\mathrm{d} i_\mathrm{a}}{\mathrm{d} t} &= \frac{u - L_\mathrm{e}^\prime i_\mathrm{e} \omega_\mathrm{me} - R_\mathrm{a} i_\mathrm{a}}{L_\mathrm{a}} \\
    \frac{\mathrm{d} i_\mathrm{e}}{\mathrm{d} t} &= \frac{u - R_\mathrm{e} i_\mathrm{e}}{L_\mathrm{e}}


Torque Equation
***************
.. math::
    T = L_\mathrm{e}^\prime i_\mathrm{e} i_\mathrm{a}


Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.DcShuntMotor
   :members:
   :inherited-members:
Base DC Motor
############################

Code Documentation
******************

.. autoclass:: gym_electric_motor.physical_systems.electric_motors.DcMotor
   :members:
   :inherited-members:
scipy.integrate.solve_ivp Solver
################################

.. autoclass:: gym_electric_motor.physical_systems.solvers.ScipySolveIvpSolver
    :members:
    :inherited-members:
scipy.integrate.ode Solver
##########################

.. autoclass:: gym_electric_motor.physical_systems.solvers.ScipyOdeSolver
    :members:
    :inherited-members:

Euler Solver
############

.. autoclass:: gym_electric_motor.physical_systems.solvers.EulerSolver
    :members:
    :inherited-members:

ODE-Solvers
###########

Solving of ODE-Systems in the form

.. math::
    \frac{\mathrm{d} \mathbf{x}}{\mathrm{d} t}&= f(\mathbf{x}, \mathbf{u}, t)\\


..  toctree::
    :maxdepth: 1
    :caption: Available ODE-Solvers:

    euler
    scipy_solve_ivp
    scipy_ode
    scipy_odeint



ODE-Solver Base Class
'''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.solvers.OdeSolver
    :members:
scipy.integrate.odeint Solver
#############################

.. autoclass:: gym_electric_motor.physical_systems.solvers.ScipyOdeIntSolver
    :members:
    :inherited-members:
Noise Generators
################

The Noise Generators generate noise which is added to the state variables.

..  toctree::
    :maxdepth: 1
    :caption: Available Noise Generators:

    gaussian_white_noise


Noise Generator Base Class
''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.noise_generators.NoiseGenerator
   :members:Gaussian White Noise Generator
##############################

.. autoclass:: gym_electric_motor.physical_systems.noise_generators.NoiseGenerator
    :members:
    :inherited-members:Four Quadrant Converters
#########################

.. figure:: ../../../plots/4QC.svg

Discrete Four Quadrant Converter
''''''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.FiniteFourQuadrantConverter
   :members:
   :inherited-members:

Continuous Four Quadrant Converter
'''''''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.ContFourQuadrantConverter
   :members:
   :inherited-members:
Two Quadrant Converters
#######################

.. figure:: ../../../plots/2QC.svg

Discrete Two Quadrant Converter
'''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.FiniteTwoQuadrantConverter
   :members:
   :inherited-members:

Continuous Two Quadrant Converter
'''''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.ContTwoQuadrantConverter
   :members:
   :inherited-members:
One Quadrant Converters
#######################

.. figure:: ../../../plots/1QC.svg

Discrete One Quadrant Converter
'''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.FiniteOneQuadrantConverter
   :members:
   :inherited-members:

Continuous One Quadrant Converter
'''''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.ContOneQuadrantConverter
   :members:
   :inherited-members:
No Converter
'''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.NoConverter
   :members:
   :inherited-members:
Multi Converters
#######################

The Multi Converters allow to include an arbitrary number of (discrete or continuous) subconverters for the use in e.g. the Externally Excited Dc Motor.
Subconverters must be 'elementary' and can not be Multi Converters.


Discrete Multi Converter
'''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.FiniteMultiConverter
   :members:
   :inherited-members:

Continuous Multi Converter
'''''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.ContMultiConverter
   :members:
   :inherited-members:
Three Phase Converters
#######################

.. figure:: ../../../plots/B6.svg

Discrete B6 Bridge Converter
'''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.FiniteB6BridgeConverter
   :members:
   :inherited-members:

Continuous B6 Bridge Converter
'''''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.ContB6BridgeConverter
   :members:
   :inherited-members:
Power Electronic Converters
##################################

..  toctree::
    :maxdepth: 1
    :caption: Available Converters:

    1QC
    2QC
    4QC
    B6C
    DoubleConv
    NoConv

The converters are divided into two classes: The discretely controlled and continuously controlled converters.
The PowerElectronicConverter class is the base class for all converters. From that, the DiscreteConverter and the
ContinuousDynamicallyAveragedConverter derive to be the base class for all continuous and discrete converters.


Converter Base Class
'''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.PowerElectronicConverter
   :members:

Finite Control Set Converter
'''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.FiniteConverter
   :members:
   :inherited-members:


Continuous Control Set Dynamically Averaged Converter
'''''''''''''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.physical_systems.converters.ContDynamicallyAveragedConverter
   :members:
   :inherited-members:
Ornstein Uhlenbeck Load
######################

Class Description
''''''''''''''''''
.. autoclass:: gym_electric_motor.physical_systems.mechanical_loads.OrnsteinUhlenbeckLoad
    :members:
    :inherited-members:Mechanical Loads
################

..  toctree::
    :maxdepth: 1
    :caption: Available Mechanical Loads:

    polystatic
    const_speed_load
    ext_speed_load
    ornstein_uhlenbeck_load

MechanicalLoad Base Class
*************************************

.. autoclass:: gym_electric_motor.physical_systems.mechanical_loads.MechanicalLoad
   :members:Constant Speed Load
######################

Class Description
''''''''''''''''''
.. autoclass:: gym_electric_motor.physical_systems.mechanical_loads.ConstantSpeedLoad
    :members:
    :inherited-members:
External Speed Load
######################

Class Description
''''''''''''''''''
.. autoclass:: gym_electric_motor.physical_systems.mechanical_loads.ExternalSpeedLoad
    :members:
    :inherited-members:Polynomial Static Load
######################

Mechanical ODE
''''''''''''''

.. math::
    \frac{ \mathrm{d} \omega_\mathrm{me} } { \mathrm{d} t } = \frac{ T - T_\mathrm{L} (\omega_\mathrm{me})}{J_\mathrm{total}}


Polynomial Load Equation
''''''''''''''''''''''''

.. math::
    T_\mathrm{L}(\omega_\mathrm{me})=\mathrm{sign}(\omega_\mathrm{me})(c \omega^2_\mathrm{me} + b \vert\omega_\mathrm{me}\vert + a)\\

Class Description
''''''''''''''''''
.. autoclass:: gym_electric_motor.physical_systems.mechanical_loads.PolynomialStaticLoad
    :members:
    :inherited-members:
Wiener Process Reference Generator
##################################

.. autoclass:: gym_electric_motor.reference_generators.WienerProcessReferenceGenerator
    :members:
    :inherited-members:Reference Generators
#####################

..  toctree::
    :maxdepth: 1
    :caption: Available Reference Generators:

    subepisoded_reference_generator
    wiener_process_reference_generator
    sinusoidal_reference_generator
    step_reference_generator
    triangular_reference_generator
    sawtooth_reference_generator
    const_reference_generator
    zero_reference_generator
    multiple_ref_generator
    switched_reference_generator


Reference Generator Base Class
''''''''''''''''''''''''''''''

.. autoclass:: gym_electric_motor.core.ReferenceGenerator
   :members:

Sawtooth Reference Generator
############################

.. autoclass:: gym_electric_motor.reference_generators.SawtoothReferenceGenerator
   :members:
   :inherited-members:
Multiple Reference Generator
############################

.. autoclass:: gym_electric_motor.reference_generators.MultipleReferenceGenerator
   :members:
   :inherited-members:
Sinusoidal Reference Generator
##############################

.. autoclass:: gym_electric_motor.reference_generators.SinusoidalReferenceGenerator
   :members:
   :inherited-members:
Zero Reference Generator
##################################

.. autoclass:: gym_electric_motor.reference_generators.ZeroReferenceGenerator
    :members:
    :inherited-members:Subepisoded Reference Generator
###############################

.. autoclass:: gym_electric_motor.reference_generators.subepisoded_reference_generator.SubepisodedReferenceGenerator
   :members:
   :inherited-members:
Triangular Reference Generator
##############################

.. autoclass:: gym_electric_motor.reference_generators.TriangularReferenceGenerator
   :members:
   :inherited-members:
Constant Reference Generator
############################

.. autoclass:: gym_electric_motor.reference_generators.ConstReferenceGenerator
   :members:
   :inherited-members:
Switched Reference Generator
#############################

.. autoclass:: gym_electric_motor.reference_generators.SwitchedReferenceGenerator
   :members:
   :inherited-members:
Step Reference Generator
########################

.. autoclass:: gym_electric_motor.reference_generators.StepReferenceGenerator
   :members:
   :inherited-members:
Squared Constraint
###################

.. autoclass:: gym_electric_motor.constraints.SquaredConstraint
   :members:

Limit Constraint
#################

.. autoclass:: gym_electric_motor.constraints.LimitConstraint
   :members:

Constraint Base Class
#####################

How To: Define Your Own Constraints
________________________________


Constraint API Documentation
____________________________

.. autoclass:: gym_electric_motor.core.RewardFunction
   :members:
Environments
############

On this page, all environments with their environment-id are listed.
In general, all environment-ids are structured as follows:

``ControlType-ControlTask-MotorType-v0``

- The ``ControlType`` is in ``{Finite / Cont}`` for all DC Motors and in ``{Finite / AbcCont / DqCont}`` for all AC Motors
- The ``ControlTask`` is in ``{TC / SC / CC}`` (Torque / Speed / Current Control)
- The ``MotorType`` is in ``{PermExDc / ExtExDc / SeriesDc / ShuntDc / PMSM / SynRM / DFIM / SCIM }``


=================================================================== ==============================
Environment                                                         environment-id
=================================================================== ==============================
**Permanently Excited DC Motor Environments**

Discrete Torque Control Permanently Excited DC Motor Environment     ``'Finite-TC-PermExDc-v0'``
Continuous Torque Control Permanently Excited DC Motor Environment   ``'Cont-TC-PermExDc-v0'``
Discrete Speed Control Permanently Excited DC Motor Environment      ``'Finite-SC-PermExDc-v0'``
Continuous Speed Control Permanently Excited DC Motor Environment    ``'Cont-SC-PermExDc-v0'``
Discrete Current Control Permanently Excited DC Motor Environment    ``'Finite-CC-PermExDc-v0'``
Continuous Current Control Permanently Excited DC Motor Environment  ``'Cont-CC-PermExDc-v0'``

**Externally Excited DC Motor Environments**

Discrete Torque Control Externally Excited DC Motor Environment      ``'Finite-TC-ExtExDc-v0'``
Continuous Torque Control Externally Excited DC Motor Environment    ``'Cont-TC-ExtExDc-v0'``
Discrete Speed Control Externally Excited DC Motor Environment       ``'Finite-SC-ExtExDc-v0'``
Continuous Speed Control Externally Excited DC Motor Environment     ``'Cont-SC-ExtExDc-v0'``
Discrete Current Control Externally Excited DC Motor Environment     ``'Finite-CC-ExtExDc-v0'``
Continuous Current Control Externally Excited DC Motor Environment   ``'Cont-CC-ExtExDc-v0'``

**Series DC Motor Environments**

Discrete Torque Control Series DC Motor Environment                  ``'Finite-TC-SeriesDc-v0'``
Discrete Torque Control Series DC Motor Environment                  ``'Cont-TC-SeriesDc-v0'``
Discrete Speed Control  Series DC Motor Environment                  ``'Finite-SC-SeriesDc-v0'``
Continuous Speed Control Series DC Motor Environment                 ``'Cont-SC-SeriesDc-v0'``
Discrete Current Control Series DC Motor Environment                 ``'Finite-CC-SeriesDc-v0'``
Continuous Current Control Series DC Motor Environment               ``'Cont-CC-SeriesDc-v0'``

**Shunt DC Motor Environments**

Discrete Torque Control Shunt DC Motor Environment                   ``'Finite-TC-ShuntDc-v0'``
Continuous Torque Control Shunt DC Motor Environment                 ``'Cont-TC-ShuntDc-v0'``
Discrete Speed Control Shunt DC Motor Environment                    ``'Finite-SC-ShuntDc-v0'``
Continuous Speed Control Shunt DC Motor Environment                  ``'Cont-SC-ShuntDc-v0'``
Discrete Current Control Shunt DC Motor Environment                  ``'Finite-CC-ShuntDc-v0'``
Continuous Current Control Shunt DC Motor Environment                ``'Cont-CC-ShuntDc-v0'``

**Permanent Magnet Synchronous Motor (PMSM) Environments**

Finite Torque Control PMSM Environment                               ``'Finite-TC-PMSM-v0'``
Abc-Continuous Torque Control PMSM Environment                       ``'AbcCont-TC-PMSM-v0'``
Dq-Continuous Torque Control PMSM Environment                        ``'DqCont-TC-PMSM-v0'``
Finite Speed Control PMSM Environment                                ``'Finite-SC-PMSM-v0'``
Abc-Continuous Speed Control PMSM Environment                        ``'Abc-Cont-SC-PMSM-v0'``
Dq-Continuous Speed Control PMSM Environment                         ``'Dq-Cont-SC-PMSM-v0'``
Finite Current Control PMSM Environment                              ``'Finite-CC-PMSM-v0'``
Abc-Continuous Current Control PMSM Environment                      ``'AbcCont-CC-PMSM-v0'``
Dq-Continuous Current Control PMSM Environment                       ``'DqCont-CC-PMSM-v0'``

**Synchronous Reluctance Motor (SynRM) Environments**

Finite Torque Control SynRM Environment                              ``'Finite-TC-SynRM-v0'``
Abc-Continuous Torque Control SynRM Environment                      ``'AbcCont-TC-SynRM-v0'``
Dq-Continuous Torque Control SynRM Environment                       ``'DqCont-TC-SynRM-v0'``
Finite Speed Control SynRM Environment                               ``'Finite-SC-SynRM-v0'``
Abc-Continuous Speed Control SynRM Environment                       ``'Abc-Cont-SC-SynRM-v0'``
Dq-Continuous Speed Control SynRM Environment                        ``'Dq-Cont-SC-SynRM-v0'``
Finite Current Control SynRM Environment                             ``'Finite-CC-SynRM-v0'``
Abc-Continuous Current Control SynRM Environment                     ``'AbcCont-CC-SynRM-v0'``
Dq-Continuous Current Control SynRM Environment                      ``'DqCont-CC-SynRM-v0'``

**Squirrel Cage Induction Motor (SCIM) Environments**

Finite Torque Control SCIM Environment                               ``'Finite-TC-SCIM-v0'``
Abc-Continuous Torque Control SCIM Environment                       ``'AbcCont-TC-SCIM-v0'``
Dq-Continuous Torque Control SCIM Environment                        ``'DqCont-TC-SCIM-v0'``
Finite Speed Control SCIM Environment                                ``'Finite-SC-SCIM-v0'``
Abc-Continuous Speed Control SCIM Environment                        ``'Abc-Cont-SC-SCIM-v0'``
Dq-Continuous Speed Control SCIM Environment                         ``'Dq-Cont-SC-SCIM-v0'``
Finite Current Control SCIM Environment                              ``'Finite-CC-SCIM-v0'``
Abc-Continuous Current Control SCIM Environment                      ``'AbcCont-CC-SCIM-v0'``
Dq-Continuous Current Control SCIM Environment                       ``'DqCont-CC-SCIM-v0'``

**Doubly Fed Induction Motor (DFIM) Environments**

Finite Torque Control DFIM Environment                               ``'Finite-TC-DFIM-v0'``
Abc-Continuous Torque Control DFIM Environment                       ``'AbcCont-TC-DFIM-v0'``
Dq-Continuous Torque Control DFIM Environment                        ``'DqCont-TC-DFIM-v0'``
Finite Speed Control DFIM Environment                                ``'Finite-SC-DFIM-v0'``
Abc-Continuous Speed Control DFIM Environment                        ``'Abc-Cont-SC-DFIM-v0'``
Dq-Continuous Speed Control DFIM Environment                         ``'Dq-Cont-SC-DFIM-v0'``
Finite Current Control DFIM Environment                              ``'Finite-CC-DFIM-v0'``
Abc-Continuous Current Control DFIM Environment                      ``'AbcCont-CC-DFIM-v0'``
Dq-Continuous Current Control DFIM Environment                       ``'DqCont-CC-DFIM-v0'``
=================================================================== ==============================

.. toctree::
   :maxdepth: 3
   :caption: Motor Environments:
   :glob:

   permex_dc/permex_dc_envs
   extex_dc/extex_dc_envs
   series_dc/series_dc_envs
   shunt_dc/shunt_dc_envs
   pmsm/pmsm_envs
   synrm/synrm_envs
   scim/scim_envs
   dfim/dfim_envs


Electric Motor Base Environment
'''''''''''''''''''''''''''''''

.. automodule:: gym_electric_motor.core

.. figure:: ../../plots/TopLevelStructure.svg

.. autoclass:: gym_electric_motor.core.ElectricMotorEnvironment
   :members:
Abc-Continuous Current Control Doubly Fed Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContCurrentControlDoublyFedInductionMotorEnv
   :members:
Finite Control Set Current Control Doubly Fed Induction Motor Environment
**************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlDoublyFedInductionMotorEnv
   :members:
Finite Control Set Torque Control Doubly Fed Induction Motor Environment
**************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlDoublyFedInductionMotorEnv
   :members:
Abc-Continuous Speed Control Doubly Fed Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContSpeedControlDoublyFedInductionMotorEnv
   :members:
Abc-Continuous Torque Control Doubly Fed Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContTorqueControlDoublyFedInductionMotorEnv
   :members:
Dq-Continuous Current Control Doubly Fed Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.DqContCurrentControlDoublyFedInductionMotorEnv
   :members:
Finite Control Set Speed Control Doubly Fed Induction Motor Environment
**************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlDoublyFedInductionMotorEnv
   :members:
Dq-Continuous Speed Control Doubly Fed Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.DqContSpeedControlDoublyFedInductionMotorEnv
   :members:
Doubly Fed Induction Motor Environments
***************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Dq-Continuous Torque Control Doubly Fed Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.DqContTorqueControlDoublyFedInductionMotorEnv
   :members:
Finite Control Set  Torque Control DC Externally Excited Motor Environment
********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlDcExternallyExcitedMotorEnv
   :members:
Externally Excited DC Motor Environments
******************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Finite Control Set Speed Control DC Externally Excited Motor Environment
*******************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlDcExternallyExcitedMotorEnv
   :members:
Continuous Current Control DC Externally Excited Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContCurrentControlDcExternallyExcitedMotorEnv
   :members:
Continuous Torque Control DC Externally Excited Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContTorqueControlDcExternallyExcitedMotorEnv
   :members:
Continuous Speed Control DC Externally Excited Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContSpeedControlDcExternallyExcitedMotorEnv
   :members:
Finite Control Set Current Control DC Externally Excited Motor Environment
*********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlDcExternallyExcitedMotorEnv
   :members:
Finite Control Set Speed Control Synchronous Reluctance Motor Environment
********************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlSynchronousReluctanceMotorEnv
   :members:
Dq-Continuous Speed Control Synchronous Reluctance Motor Environment
****************************************************************************
.. autoclass:: gym_electric_motor.envs.DqContSpeedControlSynchronousReluctanceMotorEnv
   :members:
Abc-Continuous Torque Control Synchronous Reluctance Motor Environment
******************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContTorqueControlSynchronousReluctanceMotorEnv
   :members:
Synchronous Reluctance Motor Environments
******************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Dq-Continuous Current Control Synchronous Reluctance Motor Environment
*****************************************************************************
.. autoclass:: gym_electric_motor.envs.DqContCurrentControlSynchronousReluctanceMotorEnv
   :members:
Abc-Continuous Speed Control Synchronous Reluctance Motor Environment
****************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContSpeedControlSynchronousReluctanceMotorEnv
   :members:
Finite Control Set Torque Control Synchronous Reluctance Motor Environment
*********************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlSynchronousReluctanceMotorEnv
   :members:
Dq-Continuous Torque Control Synchronous Reluctance Motor Environment
****************************************************************************
.. autoclass:: gym_electric_motor.envs.DqContTorqueControlSynchronousReluctanceMotorEnv
   :members:
Abc-Continuous Current Control Synchronous Reluctance Motor Environment
*****************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContCurrentControlSynchronousReluctanceMotorEnv
   :members:
Finite Control Set Current Control Synchronous Reluctance Motor Environment
**********************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlSynchronousReluctanceMotorEnv
   :members:
Dq-Continuous Torque Control Permanent Magnet Synchronous Motor Environment
****************************************************************************
.. autoclass:: gym_electric_motor.envs.DqContTorqueControlPermanentMagnetSynchronousMotorEnv
   :members:
Finite Control Set Speed Control Permanent Magnet Synchronous Motor Environment
********************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlPermanentMagnetSynchronousMotorEnv
   :members:
Permanent Magnet Synchronous Motor Environments
************************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Abc-Continuous Torque Control Permanent Magnet Synchronous Motor Environment
******************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContTorqueControlPermanentMagnetSynchronousMotorEnv
   :members:
Dq-Continuous Speed Control Permanent Magnet Synchronous Motor Environment
****************************************************************************
.. autoclass:: gym_electric_motor.envs.DqContSpeedControlPermanentMagnetSynchronousMotorEnv
   :members:
Dq-Continuous Current Control Permanent Magnet Synchronous Motor Environment
*****************************************************************************
.. autoclass:: gym_electric_motor.envs.DqContCurrentControlPermanentMagnetSynchronousMotorEnv
   :members:
Abc-Continuous Current Control Permanent Magnet Synchronous Motor Environment
*****************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContCurrentControlPermanentMagnetSynchronousMotorEnv
   :members:
Abc-Continuous Speed Control Permanent Magnet Synchronous Motor Environment
****************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContSpeedControlPermanentMagnetSynchronousMotorEnv
   :members:
Finite Control Set Torque Control Permanent Magnet Synchronous Motor Environment
*********************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlPermanentMagnetSynchronousMotorEnv
   :members:
Finite Control Set Current Control Permanent Magnet Synchronous Motor Environment
**********************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlPermanentMagnetSynchronousMotorEnv
   :members:
Finite Control Set Speed Control Squirrel Cage Induction Motor Environment
**************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlSquirrelCageInductionMotorEnv
   :members:
Finite Control Set Torque Control Squirrel Cage Induction Motor Environment
*****************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlSquirrelCageInductionMotorEnv
   :members:
Dq-Continuous Torque Control Squirrel Cage Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.DqContTorqueControlSquirrelCageInductionMotorEnv
   :members:
Abc-Continuous Current Control Squirrel Cage Induction Motor Environment
*************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContCurrentControlSquirrelCageInductionMotorEnv
   :members:
Dq-Continuous Speed Control Squirrel Cage Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.DqContSpeedControlSquirrelCageInductionMotorEnv
   :members:
Dq-Continuous Current Control Squirrel Cage Induction Motor Environment
************************************************************************
.. autoclass:: gym_electric_motor.envs.DqContCurrentControlSquirrelCageInductionMotorEnv
   :members:
Finite Control Set Current Control Squirrel Cage Induction Motor Environment
*****************************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlSquirrelCageInductionMotorEnv
   :members:
Squirrel Cage Induction Motor Environments
******************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Abc-Continuous Torque Control Squirrel Cage Induction Motor Environment
************************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContTorqueControlSquirrelCageInductionMotorEnv
   :members:
Abc-Continuous Speed Control Squirrel Cage Induction Motor Environment
**********************************************************************
.. autoclass:: gym_electric_motor.envs.AbcContSpeedControlSquirrelCageInductionMotorEnv
   :members:
Finite Control Set Current Control DC Permanently Excited Motor Environment
*********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlDcPermanentlyExcitedMotorEnv
   :members:
Continuous Torque Control DC Permanently Excited Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContTorqueControlDcPermanentlyExcitedMotorEnv
   :members:
Continuous Current Control DC Permanently Excited Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContCurrentControlDcPermanentlyExcitedMotorEnv
   :members:
Permanently Excited DC Motor Environments
******************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Finite Control Set  Torque Control DC Permanently Excited Motor Environment
********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlDcPermanentlyExcitedMotorEnv
   :members:
Finite Control Set Speed Control DC Permanently Excited Motor Environment
*******************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlDcPermanentlyExcitedMotorEnv
   :members:
Continuous Speed Control DC Permanently Excited Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContSpeedControlDcPermanentlyExcitedMotorEnv
   :members:
Series DC Motor Environments
******************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Finite Control Set Speed Control Series DC Motor Environment
*******************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlDcSeriesMotorEnv
   :members:
Finite Control Set Current Control Series DC Motor Environment
*********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlDcSeriesMotorEnv
   :members:
Continuous Current Control Series DC Motor Environment
********************************************************************
.. autoclass:: gym_electric_motor.envs.ContCurrentControlDcSeriesMotorEnv
   :members:
Finite Control Set  Torque Control Series DC Motor Environment
********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlDcSeriesMotorEnv
   :members:
Continuous Torque Control Series DC Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContTorqueControlDcSeriesMotorEnv
   :members:
Continuous Speed Control Series DC Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContSpeedControlDcSeriesMotorEnv
   :members:
Finite Control Set Current Control Shunt DC Motor Environment
*********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteCurrentControlDcShuntMotorEnv
   :members:
Continuous Speed Control Shunt DC Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContSpeedControlDcShuntMotorEnv
   :members:
Continuous Torque Control Shunt DC Motor Environment
*******************************************************
.. autoclass:: gym_electric_motor.envs.ContTorqueControlDcShuntMotorEnv
   :members:
Continuous Current Control Shunt DC Motor Environment
********************************************************************
.. autoclass:: gym_electric_motor.envs.ContCurrentControlDcShuntMotorEnv
   :members:
Shunt DC Motor Environments
******************************************


.. toctree::
   :maxdepth: 2
   :caption: Environments:
   :glob:

   *


Finite Control Set Speed Control Shunt DC Motor Environment
*******************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteSpeedControlDcShuntMotorEnv
   :members:
Finite Control Set  Torque Control Shunt DC Motor Environment
********************************************************************
.. autoclass:: gym_electric_motor.envs.FiniteTorqueControlDcShuntMotorEnv
   :members:
