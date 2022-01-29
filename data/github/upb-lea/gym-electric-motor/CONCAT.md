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
