[![CircleCI](https://circleci.com/gh/joglekara/VlaPy.svg?style=shield)](https://circleci.com/gh/joglekara/VlaPy)
[![codecov](https://codecov.io/gh/joglekara/VlaPy/branch/master/graph/badge.svg)](https://codecov.io/gh/joglekara/VlaPy)
[![Documentation Status](https://readthedocs.org/projects/vlapy/badge/?version=latest)](https://vlapy.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![status](https://joss.theoj.org/papers/c2b3924d7868d7bd8472c6deb011cfcc/status.svg)](https://joss.theoj.org/papers/c2b3924d7868d7bd8472c6deb011cfcc)
[![DOI](https://zenodo.org/badge/239145397.svg)](https://zenodo.org/badge/latestdoi/239145397)
# VlaPy

Usage details and the latest documentation can be found [here](https://vlapy.readthedocs.io/en/latest/)

## Code of Conduct
Please adhere to the guidelines from the Contributor Covenant listed in the [Code of Conduct](CODE_OF_CONDUCT.md).

## Quick Usage
To install dependencies, run ``python3 setup.py install`` from the base directory of the repository.

After this step, ``python3 run_nlepw.py`` can be executed to run a simulation of a Non-Linear Electron Plasma Wave with collisions.

This will create a temporary directory for the simulation files. Once completed, MLFlow will move the simulation folder into a centralized datastore. This datastore can be accessed through a web-browser based UI provided by leveraging MLFlow.

To start the MLFlow UI server, type ``mlflow ui`` into the terminal and then navigate to localhost:5000 in your web browser. The page will look like the following

![MLFlow UI](notebooks/screenshots_for_example/ui.png)

Clicking into that run will show you

![MLFlow damping](notebooks/screenshots_for_example/nlepw_screenshot.png)
## Overview
VlaPy is a 1-spatial-dimension, 1-velocity-dimension, Vlasov-Poisson-Fokker-Planck code written in Python. 

## Statement of Need
The 1D-1V VPFP equation set solved here has been applied in research of laser-plasma interactions in the context of 
inertial fusion, of plasma-based accelerators, of space physics, and of fundamental plasma physics (references 
can be found in the manuscript).  While there are VPFP software libraries which are available in academic settings, 
research laboratories, and industry, the community has yet to benefit from a simple-to-read, open-source Python 
implementation. This lack of capability is currently echoed in conversations within the ``PlasmaPy`` community 
(``PlasmaPy`` is a collection of open-source plasma physics resources). Our aim with ``VlaPy`` is to take a step 
towards filling this need for a research and educational tool in the open-source community.

``VlaPy`` is intended to help students learn fundamental concepts and help researchers discover novel physics and 
applications in plasma physics, fluid physics, computational physics, and numerical methods.  It is also designed to 
provide a science-accessible introduction to industry and software engineering best-practices, including unit and 
integrated testing, and extensible and maintainable code. 

The details of the ``VlaPy`` implementation are provided in the following sections. 

## Implementation
The Vlasov-Poisson-Fokker-Planck system can be decomposed into 4 components.

### Vlasov - Spatial Advection
The spatial advection operator is pushed using an exponential integrator. The system is periodic in x. 

This operator is tested in the fully integrated tests to reproduce solutions of the 
1D-1V Vlasov-Poisson system, namely, Landau damping.

### Vlasov - Velocity Advection
The velocity advection operator is pushed using an exponential integrator. The system is periodic in v.

This operator is tested in the fully integrated tests to reproduce solutions of the 
1D-1V Vlasov-Poisson system, namely, Landau damping.

 
### Poisson Solver
The Poisson equation is solved pseudospectrally. 

This solver is tested to reproduce analytical solutions to a periodic Poisson system.


### Fokker-Planck Solver
The Fokker-Planck equation is solved using an implicit finite-difference scheme because of the need to perform a 
diffusion time-step. 

This solver is tested to 
1) return df/dt = 0 if a Maxwell-Boltzmann distribution is provided as input 
2) conserve energy and density
3) relax to a Maxwellian of the right temperature and without a drift velocity

## Tests
All tests are performed in CircleCI. There are unit tests as well as integrated tests.
One of the most fundamental plasma physics phenomenon is that described by Landau damping. 

Plasmas can support electrostatic oscillations. The oscillation frequency is given by the electrostatic electron 
plasma wave (EPW) dispersion relation. When a wave of sufficiently small amplitude is driven at the resonant 
wave-number and frequency pairing, there is a resonant exchange of energy between the plasma and the electric field, 
and the electrons can damp the electric field.

In VlaPy, we verify that the damping rate is reproduced for a few different wave numbers. 
This is shown in `notebooks/landau_damping.ipynb.`

We include validation against this phenomenon as an integrated test.

## Other practical considerations
### File Storage
XArray enables a user-friendly interface to labeling multi-dimensional arrays along with a powerful and performant
backend. Therefore, we use XArray (http://xarray.pydata.org/en/stable/) for a performant Pythonic storage mechanism 
that promises lazy loading and incremental writes (through some tricks).

### Simulation Management
We use MLFlow (https://mlflow.org/) for simulation management. This is typically used for managing machine-learning
lifecycles but is perfectly suited for managing numerical simulations. We believe UI capability to manage simulations
significantly eases the physicist's workflow. 

There are more details about how the diagnostics for a particular type of simulation are packaged and provided to
the run manager object. These will be described in time. One can infer these from the code as well. 

## Contributing to VlaPy
Please see the guide in [contribution guidelines for this project](CONTRIBUTING.md)

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
[archis.joglekar@gmail.com].
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing to VlaPy
We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## We Develop with Github
We use github to host code, to track issues and feature requests, as well as accept pull requests.

## We Use [Github Flow](https://guides.github.com/introduction/flow/index.html), So All Code Changes Happen Through Pull Requests
Pull requests are the best way to propose changes to the codebase (we use [Github Flow](https://guides.github.com/introduction/flow/index.html)). We actively welcome your pull requests:

1. Fork the repo and create your branch from `master`.
2. If you've added code, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. Issue that pull request!

## Any contributions you make will be under the MIT Software License
In short, when you submit code changes, your submissions are understood to be under the same [MIT License](http://choosealicense.com/licenses/mit/) that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using Github's [issues](https://github.com/joglekara/VlaPy/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](https://github.com/joglekara/VlaPy/issues/new/choose); it's that easy!

## Write bug reports with detail, background, and sample code

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can. 
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

People *love* thorough bug reports. I'm not even kidding.

## Use a Consistent Coding Style
We're using [Black](https://black.readthedocs.io/en/stable/) for code linting. 
It's particularly opinionated so you don't have to be. 

## License
By contributing, you agree that your contributions will be licensed under its MIT License.

## References
This document was adapted from the open-source contribution guidelines found [here](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62)
---
title: 'VlaPy: A Python package for Eulerian Vlasov-Poisson-Fokker-Planck Simulations'
tags:
  - Python
  - plasma physics
  - dynamics
  - astrophysics
  - fusion
authors:
  - name: Archis S. Joglekar
    orcid: 0000-0003-3599-5629
    affiliation: "1"
  - name: Matthew C. Levy
    orcid: 0000-0002-7387-0256
    affiliation: "1"
affiliations:
 - name: Noble.AI, 8 California St, Suite 400, San Francisco, California 94111
   index: 1
date: 16 February 2020
bibliography: paper.bib

---

# Summary

Here we introduce ``VlaPy``: a one-spatial-dimension, one-velocity-dimension (1D-1V), Eulerian Vlasov-Poisson-Fokker-Planck (VPFP) simulation code written in Python.  

The Vlasov-Poisson-Fokker-Planck system of equations is commonly used to study plasma and fluid physics in a broad set of topical environments, ranging from space physics, to laboratory-created plasmas for fusion applications [@Betti2016; @Fasoli2016; @Ongena2016; @Chen2019]. More specifically, the Vlasov-Poisson system of equations is typically employed to model collisionless plasmas. The Fokker-Planck operator can be introduced into this system to represent the effect of collisions. The primary advantage of this scheme is that instead of relying on numerical diffusion to smooth small-scale structures that arise when modeling collisionless plasmas, the Fokker-Planck operator enables a physics-based smoothing mechanism. 

Our implementation is based on finite-difference and pseudo-spectral methods. At the lowest level, ``VlaPy`` evolves a two-dimensional (2D) grid according to this set of coupled partial integro-differential equations over time. In ``VlaPy``, the simulation dynamics can be initialized through user-specified initial conditions or external forces.

# Statement of Need

The 1D-1V VPFP equation set solved here has been applied in research on laser-plasma interactions in the context of inertial fusion [@Fahlen2009; @Banks2016], plasma-based accelerators [@Thomas2016], space physics [@Chen2019], and fundamental plasma physics [@Pezzi2016; @Heninger2018].  While there are VPFP software libraries which are available in academic settings, research laboratories, and industry [@Banks2017; @Joglekar2018], the community has yet to benefit from a simple-to-read, open-source Python implementation. This lack of capability is currently echoed in conversations within the ``PlasmaPy`` [@plasmapy] community (``PlasmaPy`` is a collection of open-source plasma physics resources). Our aim with ``VlaPy`` is to take a step towards filling this need for a research and educational tool in the open-source community.

``VlaPy`` is intended to help students learn fundamental concepts and help researchers discover novel physics and applications in plasma physics, fluid physics, computational physics, and numerical methods.  It is also designed to provide a science-accessible introduction to industry and software engineering best-practices, including unit and integrated testing, and extensible and maintainable code. 

The details of the ``VlaPy`` implementation are provided in the following sections. 


# Equations

The Vlasov-Poisson-Fokker-Planck system can be decomposed into four components. These components, represented using normalized units, are 
$\tilde{v} = v/v_{th}$, $\tilde{t} = t / \omega_p^{-1}$, $\tilde{x} = x / \lambda_D$, $\tilde{m} = m / m_e$, $\tilde{q} = q/e$, $\tilde{m} = m / m_e$, $\tilde{E} = e E / m_e v_{th} \omega_p$, $\tilde{f} = f / n_e v_{th}^{-3}$ 
where $v_{th}$ is the thermal velocity, $\omega_p$ is the electron plasma frequency, $m_e$ is the electron mass, $\lambda_D$ is the Debye length, and $e$ is the elementary charge. 
The Fourier transform operator is represented by $\mathcal{F}$ and the subscript to the operator indicates the dimension of the transform. In what follows, we have omitted the 
tilde for brevity. 
 

## Vlasov Equation

The normalized, non-relativistic ($\gamma=1$) Vlasov equation for electrons is given by
$$ \frac{\partial f}{\partial t} + v  \frac{\partial f}{\partial x} - E(x) \frac{\partial f}{\partial v} = 0, $$

where $f = f(x,v,t)$ is the electron velocity distribution function.

We use operator splitting to advance the time-step [@Strang1968]. Each one of those operators is then integrated pseudo-spectrally using the following methodology.

We use the Fourier expansions of the distribution function, which are given by
$$f(x_l,v_j) = \sum \hat{f_x}(k_x, v_j) \exp(i k_x x_l) = \sum \hat{f_v}(x_l, k_v) \exp(i k_v v_j).$$

We first discretize $f(x,v,t) = f^n(x_l, v_j)$, and then perform a Fourier expansion in $\hat{x}$ for each grid value of $v$. 

This gives

$$ f^n(x_l, v_j) = \sum \hat{f}^n_x(k_x, v_j) \exp(i k_x x_j) $$

which is substituted into the Fourier transform of the advection operator in $\hat{x}$, as given by 
$$ \mathcal{F}_x\left[ \frac{\partial f}{\partial t} = - v \frac{\partial f}{\partial x} \right].$$

This process enables the decoupling of $\hat{x}$ and $\hat{v}$ grids from the time dimension and allows us to write an Ordinary Differential Equation in time for the discretized distribution function $\hat{f_x}^n(k_x, v_j)$. This is given by 

$$\frac{d \left[\hat{f_x}^n (k_x, v_j) \right]}{\hat{f_x}^n (k_x, v_j)} = -v_j~ (i k_x)~ dt. $$

Next, we solve for the change in the plasma distribution function, integrate in time, and evaluate the integral at $\hat{f_x}^n$ and $\hat{f_x}^{n+1}$ which gives

$$ \hat{f_x}^{n+1}(k_x, v_j) = \exp(-i k_x ~ v_j \Delta t) ~~ \hat{f_x}^n(k_x, v_j). $$ 

The $E \partial f/\partial v$ term is evolved similarly using
$$ \hat{f_v}^{n+1}(x_l, k_v) = \exp(-i k_v ~ E_l \Delta t) ~~ \hat{f_v}^n(x_l, k_v). $$

We have implemented a simple Leapfrog scheme as well as a 4th order integrator called the 
Position-Extended-Forest-Ruth-Like Algorithm (PEFRL) [@Omelyan2002]

### Tests
The implementation of this equation is tested in the integrated tests section.

## Poisson Equation

The normalized Poisson equation is simply
$$  \nabla^2 \Phi = - \rho $$

Because the ion species are effectively static and form a charge-neutralizing background to the electron dynamics, we can express the Poisson equation as
$$ - \nabla E = - \rho_{net} = -(1 - \rho_e) $$ 

This is justifed by the assumption that the relevant time-scales are short compared to the time-scale associated to ion motion.

In one spatial dimension, this can be expressed as

$$ \frac{\partial}{\partial x} E(x) = 1 - \int f(x,v) ~dv $$

and the discretized version that is solved is

$$  E(x_i)^{n} = \mathcal{F}_x^{-1}\left[\frac{\mathcal{F}_x\left(1 - \sum^j f^n(x_i,v_j) \Delta v\right)}{i k_x}\right] $$

### Integrated Code Testing
Unit tests are provided for this operator to validate its performance and operation under the above assumptions. These are simply unit tests against analytical solutions of integrals of periodic functions. They can be found in 
`tests/test_fieldsolver.py`.

Below, we provide an example illustration of this validation. The code is provided in 
`notebooks/test_poisson.ipynb`

![](../notebooks/screenshots_for_example/poisson_solver.pdf)


## Fokker-Planck Equation

We have implemented two simplified versions of the full Fokker-Planck operator [@Lenard1958; @Dougherty1964]. 

The first of these implementations (LB) has the governing equation given by
$$\left(\frac{\delta f}{\delta t}\right)_{\text{coll}} = \nu \frac{\partial}{\partial v} \left( v f + v_0^2 \frac{\partial f}{\partial v}\right), $$
where 
$$v_0^2 = \int v^2 f(x,v) ~ dv, $$ 
is the thermal velocity of the distribution. 

The second of these implementations (DG) has a governing equation given by
$$\left(\frac{\delta f}{\delta t}\right)_{\text{coll}} = \nu \frac{\partial}{\partial v} \left ( (v-\underline{v}) f + v_{t}^2 \frac{\partial f}{\partial v}\right), $$
where 
$$\underline{v} = \int v f(x,v) ~ dv,$$ 
is the mean velocity of the distribution and 
$$v_{t}^2 = \int (v-\underline{v})^2 f(x,v) ~ dv, $$ 
is the thermal velocity of the shifted distribution.

The second implementation is an extension of the first, and extends momentum conservation for distributions that have a non-zero mean velocity. 

We discretize this backward-in-time, centered-in-space. This procedure results in the time-step scheme given by
$$ f^{n} = \left[LD \times \bar{v}_{j+1}f^{n+1}_{j+1} + DI \times f^{n+1}_j + UD \times \bar{v}_{j-1}f^{n+1}_{j-1}  \right]. $$
$$ LD = {\Delta t} \nu \left(-\frac{v_{0,t}^2}{\Delta v^2} - \frac{1}{2\Delta v}\right) $$
$$ DI =  \left(1+2{\Delta t} \nu \frac{v_{0,t}^2}{\Delta v^2}\right) $$
$$ UD = {\Delta t} \nu \left(-\frac{v_{0,t}^2}{\Delta v^2} + \frac{1}{2\Delta v}\right) $$
where $\bar{v} = v$ or $\bar{v} = v - \underline{v}$ depending on the implementation. 

This forms a tridiagonal system of equations that can be directly inverted.

### Integrated Code Testing
Unit tests are provided for this operator. They can be found in `tests/test_lb.py` and `tests/test_dg.py`. 
The unit tests ensure that

1. The operator does not impact a Maxwell-Boltzmann distribution already satisfying $v_{th} = v_0$.

2. The LB operator conserves number density, momentum, and energy when initialized with a zero mean velocity.

3. The DG operator conserves number density, momentum, and energy when initialized with a non-zero mean velocity.

The `notebooks/test_fokker_planck.ipynb` notebook contains illustrations and examples for these tests. Below, we show results from some of the tests for illustrative purposes. 

![](../notebooks/collision_tests_plots/Maxwell_Solution.pdf)

![](../notebooks/collision_tests_plots/LB_conservation.pdf)

![](../notebooks/collision_tests_plots/LB_no_conservation.pdf)

![](../notebooks/collision_tests_plots/DG_conservation.pdf)

We see from the above figures that the distribution relaxes to a Maxwellian. Depending on the implementation, certain characteristics of momentum conservation are enforced or avoided. 

# Integrated Code Tests against Plasma Physics: Electron Plasma Waves and Landau Damping

Landau Damping is one of the most fundamental phenomena in plasma physics. An extensive review is provided in [@Ryutov1999].  

Plasmas can support electrostatic oscillations. The oscillation frequency is given by the electrostatic electron plasma wave (EPW) dispersion relation. When a wave of sufficiently small amplitude is driven at the resonant wave-number and frequency pairing, there is a resonant exchange of energy between the plasma and the electric field, and the electrons can damp the electric field. The damping rates, as well as the resonant frequencies, are given in [@Canosa1973].

In the ``VlaPy`` simulation code, we have verified that the known damping rates for Landau Damping are reproduced, for a few different wave-numbers. This is shown in `notebooks/landau_damping.ipynb`. 

We include validation against this phenomenon as an automated integrated test. The tests can be found in 
`tests/test_landau_damping.py`

Below, we also illustrate a manual validation of this phenomenon through the fully integrated workflow of running a simulation on a local machine and sending the results to the MLFlow-driven logging mechanism. After running a properly initialized simulation, we show that the damping rate of an electron plasma wave with $k=0.3$ is reproduced accurately through the UI. This can also be computed manually (please see the testing code for details).

![](../notebooks/screenshots_for_example/ui.png)
![](../notebooks/screenshots_for_example/damping.png)

To run the entire testing suite, make sure `pytest` is installed, and call `pytest` from the root folder for the repository. Individual files can also be run by calling `pytest tests/<test_filename>.py`.

# Example Run Script For Landau Damping
    import numpy as np
    from vlapy import manager, initializers
    from vlapy.infrastructure import mlflow_helpers, print_to_screen
    from vlapy.diagnostics import landau_damping
    
    if __name__ == "__main__":
        # Pick a random wavenumber
        k0 = np.random.uniform(0.3, 0.4, 1)[0]
        
        # This is a collisionless simulation. Provide float value if collisions should be simulated
        log_nu_over_nu_ld = None
        
        # This initializes the default parameters
        all_params_dict = initializers.make_default_params_dictionary()
        
        # This calculates the roots to the EPW dispersion relation given the wavenumber
        all_params_dict = initializers.specify_epw_params_to_dict(
            k0=k0, all_params_dict=all_params_dict
        )
        
        # This specifies the collision frequency given nu_ld
        all_params_dict = initializers.specify_collisions_to_dict(
            log_nu_over_nu_ld=log_nu_over_nu_ld, all_params_dict=all_params_dict
        )
    
        # The solvers can be chosen here
        all_params_dict["vlasov-poisson"]["time"] = "leapfrog"
        all_params_dict["vlasov-poisson"]["edfdv"] = "exponential"
        all_params_dict["vlasov-poisson"]["vdfdx"] = "exponential"
    
        all_params_dict["fokker-planck"]["type"] = "lb"
        all_params_dict["fokker-planck"]["solver"] = "batched_tridiagonal"
        
        # The pulse shape can be chosen here
        pulse_dictionary = {
            "first pulse": {
                "start_time": 0,
                "t_L": 6,
                "t_wL": 2.5,
                "t_R": 20,
                "t_wR": 2.5,
                "w0": all_params_dict["w_epw"],
                "a0": 1e-7,
                "k0": k0,
            }
        }
        
        # Mlflow experiment name and location
        mlflow_exp_name = "landau-damping"
        
        # Either an IP address for your MLflow server or "local" if no server specified
        uris = {
            "tracking": "local",
        }
        
       
        # Start!
        that_run = manager.start_run(
            all_params=all_params_dict,
            pulse_dictionary=pulse_dictionary,
            diagnostics=landau_damping.LandauDamping(
                vph=all_params_dict["v_ph"],
                wepw=all_params_dict["w_epw"],
            ),
            uris=uris,
            name=mlflow_exp_name,
        )
        
        # Assess if the simulation results match the actual damping rate
        print(
            mlflow_helpers.get_this_metric_of_this_run("damping_rate", that_run),
            all_params_dict["nu_ld"],
        )



# Acknowledgements
We use xarray [@Hoyer2017] for file storage and MLFlow [@Chen2020] for experiment management. We also acknowledge the 
valuable work behind NumPy [@Harris2020] and SciPy [@2020SciPy-NMeth].

We acknowledge valuable discussions with Pierre Navarro on the implementation of the Vlasov equation.

We are grateful for the editors' and reviewers' thorough feedback that improved the software as well as manuscript.

# References
