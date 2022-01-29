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
