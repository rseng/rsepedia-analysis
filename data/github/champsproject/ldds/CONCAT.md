## Contributing to LDDS

Thank you for considering and taking the time to contribute!

The following is a set of guidelines ([borrowed from the Atom project by Github](https://github.com/atom/atom/blob/master/CONTRIBUTING.md)) for contributing to this package. As always use your best judgment and feel free to propose changes to this document in a pull request.


[Code of Conduct](#code-of-conduct)

[What should I know before I get started?](#what-should-i-know-before-i-get-started)

[How Can I Contribute?](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Pull Requests](#pull-requests)

[Additional Notes](#additional-notes)
  * [Issue and Pull Request Labels](#issue-and-pull-request-labels)

## Code of Conduct

This project and everyone participating in it is governed by the [Contributor Covenant](https://github.com/champsproject/ldds/blob/develop/CODE_OF_CONDUCT.md). 
By participating, you are expected to uphold this code. Please report unacceptable behavior to [shiba@vt.edu](mailto:shiba@vt.edu).

## What should I know before I get started?

We recommend familiarity with few basic concepts in dynamical systems as listed in the glossary in the documentation. 

## How Can I Contribute?

### Reporting Bugs

As per open-source *modus operandi*, open an issue in Github using the standard style listed [here](https://github.com/atom/atom/blob/master/CONTRIBUTING.md#additional-notes). 
Alternatively, you can send bug reports with a short description and the script to reproduce the error to the authors. <!-- #region -->
**LDDS**: Python package for Computation of Lagrangian Descriptors of Dynamical Systems.

## Table of contents
* [Description](#description)
* [Dependencies and installation](#dependencies-and-installation)
	* [Installing from source](#installing-from-source)
    * [Testing](#testing)
* [Documentation](#documentation)
* [Examples and Tutorials](#examples-and-tutorials)
* [How to cite](#how-to-cite)
* [Authors and contributors](#authors-and-contributors)
* [Contributing](#contributing)
* [License](#license)
* [Acknowledgements](#acknowledgements)


LDDS
====

[![status](https://joss.theoj.org/papers/708c2c717f803733ce017abb10f06030/status.svg)](https://joss.theoj.org/papers/708c2c717f803733ce017abb10f06030)

[![.github/workflows/draft-pdf.yml](https://github.com/champsproject/ldds/actions/workflows/draft-pdf.yml/badge.svg?branch=develop)](https://github.com/champsproject/ldds/actions/workflows/draft-pdf.yml)

[![Documentation Status](https://readthedocs.org/projects/ldds/badge/?version=latest)](https://ldds.readthedocs.io/en/latest/?badge=latest)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/champsproject/ldds/HEAD) 
<!-- #endregion -->

## Description

**LDDS** is an open-source Python module that uses the method of **Lagrangian Descriptors** (LDs) for visualization of phase space structures of dynamical systems.

First introduced by Madrid & Mancho (Chaos 2009), the method of LDs defines a function, called the Lagrangian descriptor function, defined as a scalar field that maps a set of initial condition in phase space to the "arc length" integral of their evolving trajectories in time. This approach is capable of highlighting the phase space geometry of invariant manifolds associated with periodic orbits and equilibria of dynamical systems. Over the years new formulations of the method have been developed for enhanced visualization of phase space structures. See our [online book](https://champsproject.github.io/lagrangian_descriptors) for an overview on the topic. 

In LDDS we implemented the majority of the LD variants and users can be able to work with dynamical systems with continuous flows and discrete maps in deterministic and stochastic settings, with high-dimensional vector fields described as analytical expressions or numerical data. See the [**Examples and Tutorials**](#examples) section below to have an idea of the versatility of this package.

Nonlinear dynamical systems are ubiquitous in natural and engineering sciences, such as fluid mechanics, theoretical chemistry, ship dynamics, rigid body dynamics, atomic physics, solid mechanics, condensed matter physics, mathematical biology, oceanography, meteorology and celestial mechanics [@wiggins1994normally and references therein]. We expect our package to enable researchers in these and other areas to take advantage of the LDs methods in research.


## Dependencies and installation

The `setup.py` should install the dependencies listed in
[requirements.txt](https://github.com/champsproject/ldds/blob/develop/requirements.txt) using

``` bash
pip install -r requirements.txt (or pip3 install -r requirements.txt)
```


### Installing from source
Clone the git repository and install `ldds` as a module using

``` bash
git clone git@github.com:champsproject/ldds.git
cd ldds
python setup.py install
```

<!-- #region -->
### Testing 

Test your installation with following command:

```bash
cd ldds/tests
python -m unittest
```

You should see something like the following on the terminal:

```bash
Ran 1 test in 11.173s

OK
```
<!-- #endregion -->

<!-- #region -->
## Documentation

**LDDS** uses [Sphinx](http://www.sphinx-doc.org/en/stable/) for documentation and is made available online [here](https://ldds.readthedocs.io/en/latest/?badge=latest#). To build the html version of the docs locally simply:

```bash
cd docs
make html
```

The generated html can be viewed by opening `docs/_build/html/index.html`.
<!-- #endregion -->

<!-- #region -->
## Examples and Tutorials 

You can find useful tutorials on how to use LDDS in the [tutorials](tutorials/README.md) folder.

Here we show two examples of the output contour maps produced with `ldds` for the Lagrangian Descriptor values of a deterministic ([Tutorial 2](tutorials/tutorial-2.ipynb)) and a stochastic ([Tutorial 10](tutorials/tutorial-10.ipynb)) benchmark system:

<p style="text-align:center">
<img src="paper/duffing.png" alt>
<em>Duffing oscillator with harmonic forcing.</em>
</p>


<p style="text-align:center">
<img src="paper/stoch_dgyre.png">
<em>Double-gyre with stochastic forcing.</em>
</p>
<!-- #endregion -->

## How to cite

_NOTE 19 May 2020_ : We have recently submitted LDDS to [JOSS](https://joss.theoj.org/) and now awaiting revision. 

If you use this package in your publications, provisionally you can cite the package as follows:

> LDDS: Python package for Computation of Lagrangian Descriptors of Dynamical Systems. https://github.com/champsproject/ldds 

Or if you use LaTeX:

```tex
@misc{LDDS,
  author = {B. Aguilar-Sanjuan and V. Garc{\'i}a-Garrido and V. Kraj\v{n}{\'a}k and S. Naik and S. Wiggins},
  title = {{LDDS}: {P}ython package for computing and visualizing {L}agrangian {D}escriptors in {D}ynamical {S}ystems.},
  howpublished = {\url{https://github.com/champsproject/ldds}}
}
```

<!-- #region -->
## Authors and contributors

**LDDS** is currently developed and mantained by 

* [Broncio Aguilar-Sanjuan](mailto:broncio.aguilarsanjuan@bristol.ac.uk) _University of Bristol_ , UK
* [Víctor J. García-Garrido](mailto:vjose.garcia@uah.es) _Universidad de Alcalá_ , Spain
* [Vladimír Krajňák](mailto:v.krajnak@bristol.ac.uk) _University of Bristol_ , UK
* [Shibabrat Naik](mailto:s.naik@bristol.ac.uk) _University of Bristol_ , UK


with the support of and supervision of [Prof. Stephen Wiggins](mailto:s.wiggins@bristol.ac.uk) (_University of Bristol_ , UK), under the [CHAMPS Project](https://champsproject.com/).

Contact us by email for further information or questions about **LDDS**, or suggest pull requests. Contributions improving either the code or the documentation are welcome!
<!-- #endregion -->

## Contributing 

Guidelines on how to contribute to this package can be found [here](https://github.com/champsproject/ldds/blob/develop/contributing.md) along with the code of conduct [here](https://github.com/champsproject/ldds/blob/develop/code_of_conduct.md) for engaging with the fellow contributors. As and when we receive improvements to the package, we will acknowledge the pull request and the contributor in this section.


## License

See the [LICENSE](LICENSE) file for license rights and limitations.


## Acknowledgements

We acknowledge the support of EPSRC Grant No. EP/P021123/1 [CHAMPS project](https://champsproject.com). 
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
reported by contacting the project team at shiba@vt.edu. All
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
---
title: 'LDDS: Python package for computing and visualizing Lagrangian Descriptors for Dynamical Systems'
authors:
- affiliation: 1
  name: Broncio Aguilar-Sanjuan 
  orcid: 0000-0001-8068-6417
- affiliation: 2
  name: Víctor J. García-Garrido 
  orcid: 0000-0003-0557-3193
- affiliation: 1
  name: Vladimír Krajňák 
  orcid: 0000-0001-6052-7531
- affiliation: 1
  name: Shibabrat Naik
  orcid: 0000-0001-7964-2513
- affiliation: 1
  name: Stephen Wiggins
  orcid: 0000-0002-5036-5863
output:
  pdf_document:
    fig_caption: yes
    fig_height: 3
    citation_package: natbib
  html_document:
    fig_caption: yes
    fig_height: 3
bibliography: paper.bib
biblio-style: apalike
natbiboptions: round
date: 19 May 2021
year: 2021
tags:
- Dynamical systems
- Lagrangian descriptors
affiliations:
- index: 1
  name: School of Mathematics, University of Bristol, Fry Building, Woodland Road,
    Bristol BS8 1UG, United Kingdom
- index: 2
  name: Departamento de Física y Matemáticas, Universidad de Alcalá,
    Madrid, 28871, Spain
---

## Statement of Need

Nonlinear dynamical systems are ubiquitous in natural and engineering sciences, such as fluid mechanics, theoretical chemistry, ship dynamics, rigid body dynamics, atomic physics, solid mechanics, condensed matter physics, mathematical biology, oceanography, meteorology, and celestial mechanics [@wiggins1994normally and references therein]. There have been many advances in understanding phenomena across these disciplines using the geometric viewpoint of the solutions and the underlying structures in the phase space, for example @mackay_transport_1984, @romkedar_analytical_1990, @OzoriodeAlmeida1990, @RomKedar90, @meiss_symplectic_1992, @koon_heteroclinic_2000, @waalkens_escape_2005, @meiss15, @wiggins_role_2016, @zhong_tube_2018, @zhong_geometry_2020. Chief among these phase space structures are the invariant manifolds that form a barrier between dynamically distinct solutions. In most nonlinear systems, the invariant manifolds are computed using numerical techniques that rely on some form of linearization around equilibrium points followed by continuation and globalization. However, these methods become computationally expensive and challenging when applied to the high-dimensional phase space of vector fields defined analytically, from numerical simulations or experimental data. This points to the need for techniques that can be paired with trajectory calculations, without the excessive computational overhead, and that at the same time can allow visualization along with trajectory data. The Python package `LDDS` serves this need for analyzing deterministic and stochastic, continuous and discrete high-dimensional nonlinear dynamical systems described either by an analytical vector field or from data obtained from numerical simulations or experiments.

To the best of our knowledge, no other open-source software exists for computing Lagrangian descriptors. However, a variety of computational tools are available for obtaining phase space structures in fluid mechanics, such as the identification of Lagrangian coherent structures via finite-time Lyapunov exponents (@lagrangian, @dgftle, @lcstool, @libcfd2lcs, @lcsmatlabkit, @activeBarriers), finite-size Lyapunov exponents (@lagrangian), and Eulerian coherent structures (@barriertool). Our goal with this software is to make Lagrangian descriptors available to the wider scientific community and enable the use of this method for reproducible and replicable computational dynamical systems.

## Summary and Functionalities

The `LDDS` software is a Python-based module that provides the user with the capability of analyzing the phase space structures of both continuous and discrete nonlinear dynamical systems in the deterministic and stochastic settings using Lagrangian descriptors (LDs). The main idea of this method is to define a scalar valued functional called Lagrangian descriptor as the integral of a non-negative function $g(\mathbf{x}(t);\mathbf{x}_0)$ that encodes a dynamical property of a trajectory at the initial condition, $\mathbf{x}_0$. Different formulations of the Lagrangian descriptor exist in the literature where the non-negative function $g(\mathbf{x}(t);\mathbf{x}_0)$ is given by: the arclength of a trajectory in phase space (@madrid2009ld, @mancho_2013), the arclength of a trajectory projected on the configuration space (@craven2015lagrangian), the $p$-norm or $p$-quasinorm (@lopesino2017), and the Maupertuis' action of Hamiltonian mechanics (@gonzalez2020). The approach provided by Lagrangian descriptors for revealing phase space structure has also been adapted to address discrete-time systems (maps) and stochastic systems.

Briefly, for a continuous-time dynamical system:

\begin{equation}
\dfrac{d \mathbf{x}}{dt} = \mathbf{f}\left(\mathbf{x}(t),t\right)
\end{equation}

where $\mathbf{x} \in \mathbb{R}^{n}$ and $\mathbf{f}$ is the vector field, starting from an initial condition $\mathbf{x}_0 = \mathbf{x}(t_0)$ at time $t = t_0$,  $g(\mathbf{x}(t);\mathbf{x}_0)$ is integrated along with the trajectory in forward and backward time over the interval $[t_0-\tau,t_0+\tau]$, respectively, to obtain the Lagrangian descriptor,

\begin{equation}
\mathcal{L}\left(\mathbf{x}_0,t_0,\tau\right) = \int_{t_0-\tau}^{t_0+\tau} g(\mathbf{x}(t);\mathbf{x}_0) \, dt.
\end{equation}

at the initial condition. When this computation is performed for a 2D grid of initial conditions over a long enough integration time interval and the corresponding contour map is visualized, one can detect phase space structures at the points with extremum LD values that also have a singularity (non-differentiability) (@lopesino2017, @naik2019a). 

This open-source software incorporates the following features:

* Computation of LDs for two dimensional maps.
* Computation of LDs for two dimensional continuous-time dynamical systems.
* Computation of LDs for two dimensional stochastic differential equations with additive noise.
* Computation of LDs on two-dimensional sections of Hamiltonian systems with 2 or more degrees of freedom (DoF).
* Computation of LDs for 2 DoF Hamiltonian systems where the potential energy surface is known on a grid of points in the configuration space of a chemcial reaction.
* Computation of LDs from a spatio-temporal discretization of a two-dimensional time-dependent vector field.
* Numerical gradient based identification of the invariant stable and unstable manifolds from the LD contour map.
* Addition of time-dependent external forcings in two-dimensional continuous dynamical systems.
* Different formulations of the Lagrangian descriptor available in the literature.

All the features of the package and their usage across different settings are illustrated using Jupyter notebooks as hands-on tutorials. These tutorials are meant to be worked through to understand how to set up a model dynamical system to which LDs is applied and use different options for visualizing the computational results. We believe that these notebooks are effective in integrating this method into classroom courses, independent study, and research projects. Moreover, the tutorials will encourage future contributions from the scientific community to expand the features and application of the Lagrangian descriptor method in other areas of computational science and engineering.

### Example systems {#examples}

To illustrate the use of this software, we have applied the method to some of the benchmark dynamical systems. These are included as examples:

#### Maps:

* Standard map 

The standard map is a two-dimensional map used in dynamical systems to study a number of physical systems such as the cyclotron particle accelerator or a kicked rotor (@Chirikov1969, @meiss_symplectic_1992, @Meiss2008). The equations of the discrete system are given by the expressions:

\begin{equation}
\begin{cases}
x_{n+1} = x_{n} + y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n}) \\[.2cm]
y_{n+1} = y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n})
\end{cases}
\end{equation}
where $K$ is the parameter that controls the forcing strength of the perturbation. The inverse map is described by:
\begin{equation}
\begin{cases}
    x_{n} = x_{n+1} - y_{n+1} \\[.2cm]
    y_{n} = y_{n+1} + \dfrac{K}{2\pi} \sin(2\pi (x_{n+1} - y_{n+1}))
\end{cases}
\end{equation}

In the following figure, we show the output produced by the LDDS software package for the standard map using the model parameter value $K=1.2$.

![Lagrangian descriptor contour plot for the standard map, using $p=0.5$-quasinorm and integration time $\tau=15$. \label{fig:standard_map}](stdMap.png)


#### Flows:

* Forced undamped Duffing oscillator

The Duffing oscillator is an example of a periodically driven oscillator with nonlinear elasticity (@duffing1918, @Kovacic2011). This can model the oscillations of a pendulum whose stiffness does not obey Hooke's law or the motion of a particle in a double-well potential. It is also known as a simple system that can exhibit chaos. 

As a special case, the forced undamped Duffing oscillator is described by a time-dependent Hamiltonian given by:

\begin{equation}
 H(x,p_x,t) = \dfrac{1}{2}p_x^2 - \dfrac{\alpha}{2}x^2 + \dfrac{\beta}{4}x^4 - f(t) x
\end{equation}

where $\alpha$ and $\beta$ are the model parameters and $f(t)$ is the time-dependent focing added to the system. The non-autonomous vector field that defines the dynamical system is given by:

\begin{equation}
\begin{cases}
   \dot{x} = \dfrac{\partial H}{\partial p_x} = f_1(x,p_x) = p_x \\[.2cm]
   \dot{p}_x = -\dfrac{\partial H}{\partial x} = f_2(x,p_x,t) = \alpha x - \beta x^3 + f(t)
\end{cases}
\end{equation}

In the following figure we show the output produced by the LDDS software package for the forced Duffing oscillator using the model parameter value $\alpha = \beta = 1$. The initial time is $t_0 = 0$ and the perturbation used is of the form $f(t) = A\sin(\omega t)$ where $A = 0.25$ and $\omega = \pi$.

![Lagrangian descriptor contour plot for the Duffing oscillator, using $p=0.5$-quasinorm and integration time $\tau=15$. \label{fig:duffing}](duffing.png)

* A double gyre flow with stochastic forcing

The double gyre is a recurrent pattern occurring in geophysical flows (@Coulliette2001). The stochastic dynamical system for a simplified model of this flow (@Shadden2005) with additive noise is described by the following stochastic differential equations (@balibrea2016lagrangian):

\begin{equation}
\begin{cases}
   d X_t = \left(-\pi A \sin\left(\dfrac{\pi f(X_t,t)}{s}\right)\cos\left(\dfrac{\pi Y_t}{s}\right) - \mu X_t\right) \, dt + \sigma_1 \, dW_t^1 \\[.2cm]
   d Y_t = \left(\pi A \cos\left(\dfrac{\pi f(X_t,t)}{s}\right)\sin\left(\dfrac{\pi Y_t}{s}\right)\dfrac{\partial f}{\partial x}\left(X_t,t\right) - \mu Y_t\right) \, dt + \sigma_2  \, dW_t^2
\end{cases}
\end{equation}

where $W^1$ and $W^2$ are Wiener processes and we have that:

\begin{equation}
f(X_t,t) = \varepsilon \sin(\omega t + \psi) X_t^2 + \left(1-2\varepsilon\sin(\omega t + \psi)\right) \, X_t
\end{equation}

In the following figure we show the output produced by the LDDS software package for the stochastically forced double gyre using a noise amplitude of $\sigma_1 = \sigma_2 = 0.1$. The double gyre model parameters are $A = 0.25$, $\omega = 2\pi$, $\psi = \mu = 0$, $s = 1$, $\varepsilon = 0.25$, and the initial time is $t_0 = 0$.

![Lagrangian descriptor contour plot for the Double-gyre with stochastic forcing, using $p=0.5$-quasinorm and integration time $\tau=15$. \label{fig:stoch_dgyre}](stoch_dgyre.png)

Four-dimensional phase space:

* Hénon-Heiles Hamiltonian.

The Hénon-Heiles system is a simplified model describing the restricted motion of a star around the center of a galaxy (@Henon1964). This system is a paradigmatic example of a time-independent Hamiltonian with two degrees of freedom, given by the function:

\begin{equation}
H(x, y, p_x, p_y) = \frac{1}{2} (p_x^2 + p_y^2) + \frac{1}{2} (x^2 + y^2) + x^2 y - \frac{1}{3} y^3
\end{equation}
where the vector field is:
\begin{equation}
\begin{aligned}
 \dot{x} = & \dfrac{\partial H}{\partial p_x} =  p_x \\
 \dot{y} = & \dfrac{\partial H}{\partial p_y} = p_y  \\
 \dot{p}_x = & -\dfrac{\partial H}{\partial x} =  -x - 2 x y \\
 \dot{p}_y = & -\dfrac{\partial H}{\partial y} =  -x^2 -y + y^2
\end{aligned}
\end{equation}

In the next figure, we show the computation of Lagrangian descriptors with the LDDS software package on the phase space slice described by the condition $x = 0$, $p_x > 0$ for the energy of the system $H_0 = 1/5$.

![Lagrangian descriptor contour plot for the Hénon-Heiles Hamiltonian, using $p=0.5$-quasinorm and integration time $\tau=15$. \label{fig:henon_heiles}](henonheiles.png)



## Relation to ongoing research projects

Lagrangian descriptors form the basis of several past and ongoing research projects (@alvaro1, @alvaro2, @carlos2015, @craven2015lagrangian, @craven2016deconstructing, @gg2016, @balibrea2016lagrangian, @demian2017, @craven2017lagrangian, @feldmaier2017obtaining, @junginger2017chemical, @gg2018, @ramos2018, @patra2018detecting, @naik2019a, @naik2019b, @curbelo2019a, @curbelo2019b, @revuelta2019unveiling, @GG2020a, @GG2020b, @krajnak2020manifld, @naik2020, @gonzalez2020, @katsanikas2020a). The common theme of all these projects is the investigation of phase space structures that govern phase space transport in nonlinear dynamical systems. We have also published an open-source book using Jupyter Book @jupyterbook_2020 on the theory and applications of Lagrangian descriptors @ldbook2020. This open-source package is the computational companion to the book.

## Acknowledgements

We acknowledge the support of EPSRC Grant No. EP/P021123/1 [(CHAMPS project)](https://champsproject.com/) and Office of Naval Research (Grant No. N00014-01-1-0769). 


## References

<!-- #region -->
# Tutorials

This folder contains a collection of useful Notebooks with tutorials on how to use the **LDDS** Lagrangian descriptors Python module for analyzing the phase space structure of dynamical systems.

__[Tutorial 1](tutorial-1.ipynb)__ Continuous 1DoF systems (no perturbation)

* Hamilton centre (autonomous)
* Hamilton saddle (autonomous)
* Duffing oscillator (autonomous)
* Double-gyre (nonautonomous)

__[Tutorial 2](tutorial-2.ipynb)__ Lagrangian descriptors for 1 DoF dynamical system with forcing

* Forced Duffing oscillator (nonautonomous)

__[Tutorial 3](tutorial-3.ipynb)__ Variable time integration

* Hamilton saddle-node
* Inverted Duffing oscillator

__[Tutorial 4](tutorial-4.ipynb)__ Lagrangian descriptors for a user-defined 1 DoF dynamical system

* Morse oscillator

__[Tutorial 5](tutorial-5.ipynb)__ High-dimensional continuous systems

* 2DoF Hénon-Heiles
* 2DoF Index-1 normal form saddle
* 3DoF index-1 normal form saddle

__[Tutorial 6](tutorial-6.ipynb)__ Lagrangian descriptors for a user-defined 2 DoF Hamiltonian

* Double-well Hamiltonian system


__[Tutorial 7](tutorial-7.ipynb)__ Dynamics using Potential Energy Surface (PES) data

* Discretised Hénon-Heiles

__[Tutorial 8](tutorial-8.ipynb)__ Lagrangian descriptors for maps of discrete systems

* Standard map
* Hénon map

__[Tutorial 9](tutorial-9.ipynb)__ User-defined discrete maps

* Gauss map
* Gingerbreadman map

__[Tutorial 10](tutorial-10.ipynb)__ Lagrangian descriptors for Stochastic Dynamical Systems

* Noisy saddle
* Noisy Duffing oscillator 
* Noisy double-gyre

__[Tutorial 11](tutorial-11.ipynb)__ Dynamics using a vector field dataset

* Discretised forced Duffing oscillator

__[Tutorial 12](tutorial-12.ipynb)__ Integration Time & Grid Resolution for Lagrangian descriptor Simulations

* Arnold's cat map
* Double-gyre flow

__More to come...__

We encourage contributions from users to develop Jupyter notebooks that extend the capabilities and features of the **LDDS** software package.
<!-- #endregion -->
# Two degrees of freedom quadratic normal form Hamiltonian 

The two degrees of freedom quadratic normal form Hamiltonian gives a linear system with an index one saddle. 

## Uncoupled system

.. math::
   \begin{equation}
   H_2(q_1, q_2, p_1, p_2) = \underbrace{\frac{\lambda}{2} (p_1^2 - q_1^2)}_\text{$H_r$} + \underbrace{\frac{\omega_2}{2}(q_2^2 + p^2_2)}_\text{$H_b$}\quad \lambda,\omega_2 > 0 \label{eqn:ham_nf_2dof}
   \end{equation}

with the corresponding vector field given by

.. math::
   \begin{equation}
   \begin{aligned}
   \dot{q}_1 = & \frac{\partial H_2}{\partial p_1} &= \lambda p_1, \\
   \dot{q}_2 = & \frac{\partial H_2}{\partial p_2} &= \omega_2 p_2,\\
   \dot{p}_1 = & -\frac{\partial H_2}{\partial q_1} &= \lambda q_1,\\
   \dot{p}_2 = & -\frac{\partial H_2}{\partial q_2} &= -\omega_2 q_2\\
   \end{aligned}
   \label{eqn:eom_nf_2dof}
    \end{equation}

## Coupled system


========
Examples
========


Some examples used in the tutorial and considered as benchmark systems in dynamical systems. At the end of each systems's equations, there are links to Ipython notebook with implementation and visualization details.


Discrete systems
================

1. Standard map

The standard map (kicked rotator) is a two-dimensional map used in dynamical systems to study a periodically 
kicked pendulum. Its equations of motion are given by the expressions:

.. math::
   \begin{align}
    x_{n+1} = x_{n} + y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n)),
    y_{n+1} = y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n)),
   \end{align}

where :math:`K` is the parameter that controls the forcing strength of the perturbation.
   
The inverse map is described by:

.. math::
   \begin{align}
    x_{n} = x_{n+1} - y_{n+1},
    y_{n} = y_{n+1} + \dfrac{K}{2\pi} \sin(2\pi (x_{n+1} - y_{n+1})).
   \end{align}

   
2. Hénon map

The Hénon map was introduced by Michel Hénon as a simplified model of the Poincaré section 
of the Lorenz model. The map equations are as follows:

.. math::
   \begin{align}
    x_{n+1} = a - x_{n}^2 + b y_{n},
    y_{n+1} = x_{n},
   \end{align}
   
where :math:`a,b` are the model parameters.

The inverse Hénon map is:

.. math::
   \begin{align}
    x_{n} = y_{n+1},
    y_{n} = \dfrac{x_{n+1} - a + y_{n+1}^2}{b}.
   \end{align}

`Ipython notebook for standard map and Hénon map <https://github.com/champsproject/ldds/blob/develop/tutorials/tutorial-8.ipynb>`_



Continuous systems 
==================

One degree of freedom
---------------------

1. Hamiltonian center

The Hamiltonian function:

.. math::
   H(x,p_x) = \dfrac{\omega}{2} \left( p_x^2 + x^2 \right), \label{eqn:ham_center1dof}

defines the normal form of a 1 DoF system with a center equilibrium at the origin. The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} =  \omega p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = -\omega x.
   \end{align}
   

2. `Hamiltonian saddle <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#one-degree-of-freedom-hyperbolic-equilibrium-point>`_

The Hamiltonian function:

.. math::
   H(x,p_x) = \dfrac{\lambda}{2} \left( p_x^2 - x^2 \right), \label{eqn:ham_saddle1dof}

defines the normal form of a 1 DoF system with a saddle equilibrium at the origin. The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = \lambda p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = \lambda x.
   \end{align}

3. Duffing oscillator 

a. Unforced

The Hamiltonian function:

.. math::
   H(x,p_x,t) = \dfrac{1}{2}p_x^2 - \dfrac{\alpha}{2}x^2 + \dfrac{\beta}{4}x^4, \label{eqn:ham_duff}

with :math:`\alpha,\beta>0` describes the Duffing oscillator with the associated equations of motion

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} =  p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} =  \alpha x - \beta x^3.
   \end{align}

b. Forced 
Time dependent Hamiltonian function

.. math::
   H(x,p_x,t) = \dfrac{1}{2}p_x^2 - \dfrac{\alpha}{2}x^2 + \dfrac{\beta}{4}x^4 - f(t) x, \label{eqn:ham_duff_forced}

defines the Duffing oscillator with time dependent forcing :math:`f(t)`. This package offers two predefined options for the external forcing, namely :math:`f(t) = A\mathrm{sech}(t)\sin(\omega t)` and :math:`f(t) = A\sin(\omega t)`. Other versions can be added manually by the user in the forcing function of the vector_fields.py file.

The corresponding equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} =  p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} =  \alpha x - \beta x^3 + f(t).
   \end{align}



c. Inverted 

The inverted Duffing oscillator can be obtained from Hamiltonian ~\eqref{eqn:eqn:ham_duff}, by setting the parameters :math:`\alpha = \beta = - 1`.



4. Saddle-node Hamiltonian 

This system is defined by the Hamiltonian:

.. math::
    H(x,p_x) = \dfrac{1}{2}p_x^2 + \dfrac{1}{2}x^2 + \dfrac{1}{3}x^3, \label{eqn:ham_saddnode}

and its associated equations of motion are:

.. math::
    \begin{align}
    \dot{x} = \dfrac{\partial H}{\partial p_x} =  p_x, \\
    \dot{p}_x = -\dfrac{\partial H}{\partial x} =  -x - x^2.
    \end{align} 

`Ipython notebook on saddle-node Hamiltonian and inverted Duffing oscillator <https://github.com/champsproject/ldds/blob/develop/tutorials/tutorial-3.ipynb>`_

5. Non-autonomous double-gyre flow

The double-gyre flow is a classical system popular in geophysical fluid dynamics. This non-autonomous two-dimensional dynamical system is defined by the equations:

.. math::
   \begin{align}
   \dot{x} &= -\pi A \sin\left(\dfrac{\pi f(x,t)}{s}\right) \cos\left(\dfrac{\pi y}{s}\right) - \mu x, \\[.2cm]
   \dot{y} &= \pi A \cos\left(\dfrac{\pi f(x,t)}{s}\right) \sin\left(\dfrac{\pi y}{s}\right) \dfrac{\partial f}{\partial x}\left(x,t\right) - \mu y,
   \end{align} 

where we have that :math:`f(x,t) = \varepsilon \sin(\omega t + \phi) x^2 + \left(1-2\varepsilon \sin(\omega t + \phi)\right) x`.

`Ipython notebook on Hamiltonian center, saddle, and double gyre <https://github.com/champsproject/ldds/blob/develop/tutorials/tutorial-1.ipynb>`_

Two degrees of freedom
----------------------

1. `Saddle-center <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#two-degrees-of-freedom-and-the-hyperbolic-periodic-orbit>`_ 

The Hamiltonian function:

.. math::
   H(x,y,p_x,p_y) = \dfrac{1}{2} \left( p_x^2 + p_y^2 + y^2 - x^2) \right),  \label{eqn:ham_saddle2dof}

is the normal form of a 2 DoF system with a saddle-center equilibrium point at the origin. The dynamics of any 2 DoF dynamical system near a potential index-1 saddle point is conjugate to this system.
The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x}  = p_x, \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y}  = p_y, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x}  = x, \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y}  = - y.
   \end{align}

2. Hénon-Heiles

The Hamiltonian for the Hénon-Heiles system is given:

.. math::
   H(x,y,p_x,p_y) = \dfrac{1}{2} \left( p_x^2 + p_y^2 \right) + \dfrac{1}{2} \left( x^2 + y^2 \right) + yx^2 - \dfrac{1}{3} y^3, \label{eqn:ham_hh}

and Hamilton's equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x}  = p_x, \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y} = p_y, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = - x - 2xy, \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y} = - x^2 - y + y^2.
   \end{align}

This system is a fundamental system for studying complex dynamics. Depending on the value of total energy, it can exhibit different dynamical behaviour ranging from near-integrable to completely chaotic.



Three degrees of freedom
------------------------

1. `Saddle-center-center <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#three-and-more-degrees-of-freedom-and-nhims>`_

The Hamiltonian function:

.. math::
   H(x,y,z,p_x,p_y,p_z) = \dfrac{1}{2} \left( p_x^2 + p_y^2+ p_z^2 - x^2 + y^2 + z^2) \right),  \label{eqn:ham_saddle3dof}

is the normal form of a 3 DoF system with a saddle-center-center equilibrium point at the origin (also referred to as an index-1 saddle).
The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = p_x, \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y} = p_y, \\
   \dot{z} &= \dfrac{\partial H}{\partial p_z} = p_z, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = x, \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y}= - y, \\
   \dot{p}_z &= -\dfrac{\partial H}{\partial z}= - z.
   \end{align}


`Ipython notebook on Hénon-Heiles Hamiltonian, two and three degrees of freedom quadratic normal form with index-1 saddle <https://github.com/champsproject/ldds/blob/develop/tutorials/tutorial-5.ipynb>`_
LDDS 
====

Submodules
----------

ldds.base module
----------------

.. automodule:: ldds.base
   :members:
   :undoc-members:
   :show-inheritance:

ldds.tools module
-----------------

.. automodule:: ldds.tools
   :members:
   :undoc-members:
   :show-inheritance:

ldds.vector\_fields module
--------------------------

.. automodule:: ldds.vector_fields
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: ldds
   :members:
   :undoc-members:
   :show-inheritance:
========================
Glossary of Common Terms
========================

This glossary gives the mathematical definition of dynamical systems concepts used in the software for computing Lagrangian descriptors. More detailed description can be found in books listed in the references; Ref. [wiggins2003]_. 


.. glossary::

Dynamical system
----------------

   A system that changes in time; usually described by differential equations (continuous time) or difference equations (sometimes called *maps*) (discrete time), or, possibly, some combination of the two.


Vector field, map, phase space, trajectory
------------------------------------------

   The equation of the form

   .. math::
      \begin{equation}
      \dot{x} = f(x,t; \mu)
      \end{equation}

   with :math:`x \in U \subset \mathbb{R}^n, t \in \mathbb{R}^1`, and :math:`\mu \in V \subset \mathbb{R}^p` where :math:`U, V` are open sets in :math:`\mathbb{R}^n, \mathbb{R}^p`, respectively, is an **ordinary differential equation** and :math:`f` is a **vector field**. The above form is referred to as a non-autonomous or time-dependent ODE.

   The equation of the form

   .. math::
      \begin{equation}
      x \mapsto g(x; \mu)
      \end{equation}

   as a **map or difference equation**. 

   We refer to the space of dependent variables as the **phase space**. Both of these form of equations are referred to as dynamical system. Solution of the ordinary differential equations, :math:`x(t,t_0,x_0)` will be referred to as the **trajectory or phase curve** through the point :math:`x_0` at :math:`t = t_0`.

   The goal of dynamical systems analysis is to understand the geometry of solution curves in phase space. 



Equilibrium solutions, linear stability
---------------------------------------

   We give the definition of equilibrium solutions and their stability for an autonomous vector field of the form

   .. math::
      \begin{equation}
      \dot{x} = f(x; \mu), \qquad x \in \mathbb{R}^n
      \end{equation}

   The **equilibrium solution** is a point :math:`\bar{x} \in \mathbb{R}^n` such that 
   
   .. math::
      \begin{equation}
      f(\bar{x}; \mu) = 0.
      \end{equation}

   This constant solution of a vector field, that is a point in phase space where the vector field is zero, is also referred to as a "fixed point", "steady state", "stationary point", "rest point".

   Cautionary note: In case of non-autonomous vector fields, this definition of equilibrium solutions do not apply to the "frozen" (fixing :math:`t = \bar{t}`) vector field. Further explanation based on an example is given in the Chapter 1 of Ref. [wiggins2003]_ 

   We give the definition for the type of linear stability which is most relevant for this software.
   
   For a linear, autonomous vector field on :math:`\mathbb{R}^n`:
   
   .. math::
      \begin{equation}
      \dot{x} = A x, \qquad x(0) = x_0, \qquad x \in \mathbb{R}^n
      \end{equation}

   If :math:`A` has no eigenvalues with zero real part, the linear stability of the origin is determined by the real part of the eigenvalues :math:`A`. 

   If all of the real parts of the eigenvalues are strictly less than zero, then the origin is asymptotically **stable**. If at least one of the eigenvalues of :math:`A` has real part strictly larger than zero, then the origin is **unstable**.
   
   The origin of is said to be **hyperbolic** if none of the real parts of the eigenvalues of :math:`A` have zero real parts. Hyperbolic equilibria of linear, autonomous vector fields on :math:`\mathbb{R}^N` can be either sinks, sources, or saddles.
   
   More details on classification of stability of equilibria can be found in Ref. [wiggins2017]_.


Periodic orbits and invariant manifolds
---------------------------------------

   A solution :math:`x` of a dynamical system passing through the point :math:`x(0)=x_0` is said to be a **periodic solution** of period :math:`T` if there exists :math:`T > 0` such that :math:`x(t) = x(t + T)` for all :math:`t \in \mathbb{R}`.

   In a discrete system, an orbit of a point :math:`x_0 \in \mathbb{R}^n` is said to be a **periodic orbit** of period :math:`k` if :math:`g^k(x_0) = x_0`.

   In the context of this software, it is sufficient to know that a **manifold** is a set which *locally* has the structure of Euclidean space. For the "typical" dynamical systems in applications, a manifold is a :math:`m`-dimensional surface embedded in :math:`\mathbb{R}^n`.

   **Invariant manifolds** are a set of points, which is also a **manifold**, in phase space such that trajectories with initial conditions in the set remain in the set forever. In the context of this software, we are interested in invariant manifolds of hyperbolic equilibrium points or hyperbolic periodic orbits. Thus, when trajectories on an invariant manifold asymptotically approach the equilibrium point or the periodic orbit for :math:`t \rightarrow \infty` (or :math:`t \rightarrow -\infty`), the invariant manifold is called a **stable** (or **unstable**) invariant manifold. Mathematical definitions of these manifolds can be found in Chapter 3 of [wiggins2003]_ and Chapter 6 of [wiggins2017]_.

   Lagrangian descriptors identify stable (or unstable) invariant manifolds when the initial conditions are integrated for the time interval :math:`[t_0, t_0 + \tau]` (or :math:`[t_0, t_0 - \tau]`), respectively. 


References: textbooks
^^^^^^^^^^^^^^^^^^^^^
   
.. [wiggins2003] Wiggins, S., (2003) Introduction to applied nonlinear dynamical systems and chaos (Vol. 2). Springer Science & Business Media.

.. [wiggins2017] Wiggins, S. (2017) Ordinary Differential Equations (Ver. 2), https://figshare.com/articles/Ordinary_Differential_Equations/5311612,  Figshare.



API Documentation
=================


.. toctree::
   :maxdepth: 4

   ldds
.. LDDS documentation master file, created by
   sphinx-quickstart on Fri Oct  9 15:13:24 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. Places global toc into the sidebar

:globalsidebartoc: True

============================================
Lagrangian Descriptors for Dynamical Systems
============================================

.. Define an order for the Table of Contents

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   
   introduction
   examples
   glossary
   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
============
Introduction
============


The purpose of this guide is to give a brief overview of the theory, illustrate the computation and visualization of Lagrangian descriptors for dynamical systems, and automatically generated API documentation. We recommend that the user is familiar with nonlinear dynamical systems concepts mentioned in the :ref:`glossary`. Then, reading the brief introduction to the method of Lagrangian descriptor below and benchmark systems in :ref:`examples`.

Background and motivation
=========================

Understanding the geometry of the phase space of a dynamical system is a fundamental step in developing a complete picture of the dynamics. A trajectory diagnostic technique is a quick and practical approach to discovering the geometry of the structures that characterize phase space transport. This approach has become useful in analyzing a multitude of systems across geophysical fluid dynamics and chemical reactions. Lagrangian Descriptors for Dynamical Systems, or LDDS, is an open-source software written in Python represents a contribution in this direction using the method of Lagrangian descriptors. The basic idea behind this methodology is applicable to continuous or discrete dynamical system with or without periodic, aperiodic, stochastic forcing and with or without dissipation. The method encodes the geometry of phase space structures in the initial conditions on a two dimensional section by calculating a geometric property of the trajectories obtained from the initial conditions. 

Lagrangian descriptor is a scalar functional defined on a grid of initial conditions that assigns to every point a real number obtained by accumulating the values of a non-negative function along the trajectory starting from that initial condition in forward and backward time. There are many definitions of this non-negative function for computing Lagrangian descriptors (LDs) in the literature. Examples of such functions include the arclength of the trajectory, the p-norm of the components of the vector field defining the dynamical system under study, the Maupertuis classical action, etc. There is no definition of LDs that is better than the rest, and the choice depends on the dynamical characteristics of the system addressed. In this sense, there is always a trial and error stage where one would assess the different definitions, two-dimensional sections, integration time length for a given dynamical system. The LDDS package provides the user with an interface to perform such computations in a rapid prototyping manner. 

Lagrangian descriptors
======================

The Lagrangian descriptor (LD) as presented in Refs. [madrid2009]_, [mancho2013]_ is an arc-length of a trajectory calculated on a chosen initial time :math:`t_0` and measured for fixed forward and backward integration time, :math:`\tau`. For continuous time dynamical systems, Ref. [lopesino2017]_ gives an alternative definition of the LD which is useful for proving rigorous results and can be computed along with the trajectory. It provides a characterization of the notion of singular features of the LD that facilitates a proof for detecting invariant manifolds in certain model situations.  In addition, the "additive nature" of this new definition of LD provides 
an approach for assessing the influence of each degree-of-freedom separately on the Lagrangian descriptor.  This property was used in Ref. [demian2017]_ which showed that a Lagrangian descriptor can be used to detect Lyapunov periodic orbits in the two degrees-of-freedom Hénon-Heiles Hamiltonian system. We will describe this procedure for two and three degrees-of-freedom linear autonomous Hamiltonian systems. We begin by establishing notation in the general setting of a time-dependent vector field where 

.. math::
    \begin{equation}
    \frac{d\mathbf{x}}{dt} = \mathbf{v}(\mathbf{x},t), \quad \mathbf{x} \in \mathbb{R}^n \;,\; t \in \mathbb{R}
    \end{equation}

where :math:`\mathbf{v}(\mathbf{x},t) \in C^r (r \geq 1)` in :math:`\mathbf{x}` and continuous in time. The definition of LDs depends on the initial condition :math:`\mathbf{x}_{0} = \mathbf{x}(t_0)`, on the initial time :math:`t_0` (trivial for autonomous systems) and the integration time :math:`\tau`, and the type of norm of the trajectory's components, and takes the form,


.. math::
    \begin{equation}
    M_p(\mathbf{x}_{0},t_0,\tau) = \displaystyle{\int^{t_0+\tau}_{t_0-\tau} \sum_{i=1}^{n} |\dot{x}_{i}(t;\mathbf{x}_{0})|^p \; dt} \label{eqn:M_function}
    \end{equation}

where :math:`p \in (0,1]` and :math:`\tau \in \mathbb{R}^{+}` are freely chosen parameters,  and the overdot symbol represents the derivative with respect to time. It is to be noted here that there are three formulations of the function :math:`M_p` in the literature: the arc length of a trajectory in phase space [madrid2009]_, the arc length of a trajectory projected on the configuration space [junginger2016lagrangian]_, [junginger2016transition]_, [junginger2016uncovering]_, [junginger2017chemical]_ and the sum of the :math:`p`-norm of the vector field components [lopesino2015]_, [lopesino2017]_.
Although the latter formulation of the Lagrangian descriptor developed in Refs. [lopesino2015]_, [lopesino2017]_ does not resemble the arc length, the numerical results using either of these forms have been shown to be in agreement and promise of predictive capability in geophysical flows ([delacamara2012]_, [garciagarrido2015]_, [ramos2018]_, [mendoza2014lagrangian]_). The formulation we adopt here is motivated by the fact that this allows for proving rigorous result, which we will discuss in the next section, connecting the singular features and minimum in the LD plots with NHIM and its stable and unstable manifolds. 
It follows from the result that 

.. math:: 
    \begin{align}
    \mathcal{W}^s(\mathbf{x}_0, t_0) & = \text{argmin} \; \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) \\
    \mathcal{W}^u(\mathbf{x}_0, t_0) & = \text{argmin} \; \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau)
    \end{align}

where the stable and unstable manifolds (:math:`\mathcal{W}^s(\mathbf{x}_0, t_0)` and :math:`\mathcal{W}^u(\mathbf{x}_0, t_0)`) denote the invariant manifolds at intial time :math:`t_0` and :math:`\text{argmin} (\cdot)` denotes the argument that minimizes the function :math:`\mathcal{L}^{(\cdot)}(\mathbf{x}_0, t_0, \tau)` in forward and backward time, respectively. In addition, the coordinates of the NHIM at time :math:`t_0` is given by the intersection :math:`\mathcal{W}^s(\mathbf{x}_0, t_0)` and :math:`\mathcal{W}^u(\mathbf{x}_0, t_0)` of the stable and unstable manifolds, and thus given by

.. math::
    \begin{align}
    \mathcal{M}(\mathbf{x}_0, t_0) & = \text{argmin} \; \left( \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) + \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau) \right) = \text{argmin} \; \mathcal{L}(\mathbf{x}_0, t_0, \tau) \qquad \text{NHIM}
    \end{align}


References
==========
   
   
.. [madrid2009] Madrid, J. A. J. and Mancho, A. M. (2009) Distinguished trajectories in time dependent vector fields. Chaos, 19, 013111.

.. [demian2017] Demian, A. S., and Wiggins, S. (2017). Detection of periodic orbits in Hamiltonian systems using Lagrangian descriptors. International Journal of Bifurcation and Chaos, 27(14), 1750225.

.. [lopesino2017] Lopesino, C., Balibrea-Iniesta, F., García-Garrido, V. J., Wiggins, S., and Mancho, A. M. (2017). A theoretical framework for Lagrangian descriptors. International Journal of Bifurcation and Chaos, 27(01), 1730001.

.. [lopesino2015] Lopesino, C., Balibrea, F., Wiggins, S., and Mancho, A. M. (2015). Lagrangian descriptors for two dimensional, area preserving, autonomous and nonautonomous maps. Communications in Nonlinear Science and Numerical Simulation, 27(1-3), 40–51.

.. [junginger2016lagrangian] Junginger, A. and Hernandez, R. (2016a). Lagrangian descriptors in dissipative systems. Physical Chemistry Chemical Physics, 18(44), 30282–30287. 

.. [junginger2016transition] Junginger, A., Craven, G. T., Bartsch, T., Revuelta, F., Borondo, F., Benito, R., and Hernandez, R. (2016). Transition state geometry of driven chemical reactions on time-dependent double- well potentials. Physical Chemistry Chemical Physics, 18(44), 30270–30281.

.. [junginger2016uncovering] Junginger, A. and Hernandez, R. (2016b). Uncovering the Geometry of Barrierless Reactions Using Lagrangian Descriptors. The Journal of Physical Chemistry B, 120(8), 1720–1725.


.. [junginger2017chemical] Junginger, A., Duvenbeck, L., Feldmaier, M., Main, J., Wunner, G., and Hernandez, R. (2017a). Chemical dynamics between wells across a time-dependent barrier: Self-similarity in the Lagrangian descriptor and reactive basins. The Journal of chemical physics, 147(6), 064101.


.. [mancho2013] Mancho, A. M., Wiggins, S., Curbelo, J., and Mendoza, C. (2013). Lagrangian Descriptors: A Method for Revealing Phase Space Structures of General Time Dependent Dynamical Systems. Communications in Nonlinear Science and Numerical, 18, 3530–3557.


.. [mendoza2014lagrangian] Mendoza, C., Mancho, A. M., and Wiggins, S. (2014). Lagrangian descriptors and the assessment of the predictive capacity of oceanic data sets. Nonlinear Processes in Geophysics, 21(3), 677–689.

.. [garciagarrido2015] García-Garrido, V. J., Mancho, A. M., and Wiggins, S. (2015). A dynamical systems approach to the surface search for debris associated with the disappearance of flight MH370. Nonlin. Proc. Geophys., 22, 701–712.

.. [ramos2018] Ramos, A. G., García-Garrido, V. J., Mancho, A. M., Wiggins, S., Coca, J., Glenn, S., Schofield, O., Kohut, J., Aragon, D., Kerfoot, J., Haskins, T., Miles, T., Haldeman, C., Strandskov, N., All- sup, B., Jones, C., and Shapiro., J. (2018). Lagrangian coherent structure assisted path planning for transoceanic autonomous underwater vehicle missions. Scientfic Reports, 4, 4575.

.. [delacamara2012] de la Cámara, A., Mancho, A. M., Ide, K., Serrano, E., and Mechoso, C. (2012). Routes of transport across the Antarctic polar vortex in the southern spring. J. Atmos. Sci., 69(2), 753–767.



