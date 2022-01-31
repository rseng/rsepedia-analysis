![Imagel](https://raw.githubusercontent.com/rajeshrinet/pystokes/master/examples/banner.png)


## PyStokes: phoresis and Stokesian hydrodynamics in Python  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rajeshrinet/pystokes/master?filepath=binder) ![Installation](https://github.com/rajeshrinet/pystokes/workflows/CI/badge.svg) ![Notebooks](https://github.com/rajeshrinet/pystokes/workflows/notebooks/badge.svg) [![Documentation Status](https://readthedocs.org/projects/pystokes/badge/?version=latest)](https://pystokes.readthedocs.io/en/latest/?badge=latest) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02318/status.svg)](https://doi.org/10.21105/joss.02318) [![PyPI](https://img.shields.io/pypi/v/pystokes.svg)](https://pypi.python.org/pypi/pystokes) [![Python Version](https://img.shields.io/pypi/pyversions/pystokes)](https://pypi.org/project/pystokes) [![Downloads](https://pepy.tech/badge/pystokes)](https://pepy.tech/project/pystokes)  

[About](#about) 
| [Blog](https://rajeshrinet.github.io/pystokes-blog/)
| [News](#news) 
| [Installation](#installation) 
| [Documentation](https://pystokes.readthedocs.io/en/latest/) 
| [Examples](#examples) 
| [Publications ](#publications)
| [Gallery](https://github.com/rajeshrinet/pystokes/wiki/Gallery)
| [Support](#support) 
| [License](#license)



## About

[PyStokes](https://github.com/rajeshrinet/pystokes) is a numerical library for phoresis and Stokesian hydrodynamics in Python. It uses a grid-free method, combining the integral representation of Laplace and Stokes equations, spectral expansion, and Galerkin discretization, to compute phoretic and hydrodynamic interactions between spheres with slip boundary conditions on their surfaces. The library also computes suspension scale quantities, such as rheological response, energy dissipation and fluid flow. The computational cost is quadratic in the number of particles and upto 1e5 particles have been accommodated on multicore computers. The library has been used to model suspensions of **microorganisms**,  **synthetic autophoretic particles** and **self-propelling droplets**. 

Please read the PyStokes [paper](https://doi.org/10.21105/joss.02318) and [Wiki](https://github.com/rajeshrinet/pystokes/wiki) before you use PyStokes for your research. Included below are some examples from [PyStokes Gallery](https://github.com/rajeshrinet/pystokes/wiki/Gallery): 

### Periodic orbits of active particles

![Image](https://raw.githubusercontent.com/rajeshrinet/pystokes-misc/master/gallery/2_volvox.gif)

Our work shows that the oscillatory dynamics of a pair of active particles near
a boundary, best exemplified by the fascinating dance of the green algae
*Volvox*, can be understood in terms of Hamiltonian mechanics, even though the
system does not conserve energy. Read more in the [PyStokes Gallery](https://github.com/rajeshrinet/pystokes/wiki/Gallery).
<br>

### Crystallization at a plane no-slip surface
It is well-known that crystallization of colloids approximating hard spheres
is due, paradoxically, to the higher entropy of the ordered crystalline state
compared to that of  the disordered liquid state. Out of equilibrium, no such general
principle is available to rationalize crystallization. Here, we identify a new non-equilibrium mechanism, associated with entropy production rather than entropy gain, which drives crystallization of active colloids near plane walls. Read more in the [PyStokes Gallery](https://github.com/rajeshrinet/pystokes/wiki/Gallery).


![Crystallization of active colloids](https://raw.githubusercontent.com/rajeshrinet/pystokes/master/examples/crystallite.gif) 


## News
**26th July 2019** -- PyStokes can compute hydrodynamic *and* phoretic interactions in autophoretic suspensions.  


## Installation
You can take PyStokes for a spin **without installation**: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rajeshrinet/pystokes/master?filepath=binder). Please be patient while [Binder](https://mybinder.org/v2/gh/rajeshrinet/pystokes/master?filepath=binder) loads.

### From a checkout of this repo

#### Install PyStokes and an extended list of dependencies using 
 
```bash
>> git clone https://github.com/rajeshrinet/pystokes.git
>> cd pystokes
>> pip install -r requirements.txt
>> python setup.py install
```


####  Install PyStokes and its dependencies in an [environment](https://github.com/rajeshrinet/pystokes/blob/master/environment.yml) named "pystokes" via [Anaconda](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html)


```bash
>> git clone https://github.com/rajeshrinet/pystokes.git
>> cd pystokes
>> make env
>> conda activate pystokes
>> make
```

### Via pip

Install the latest [PyPI](https://pypi.org/project/pystokes) version

```bash
>> pip install pystokes
```


### Testing
Test installation and running

```bash
>> cd tests
>> python shortTests.py
```

Long test of example notebooks 

```bash
>> cd tests
>> python notebookTests.py
```


## Examples


```Python
# Example 1: Flow field due to $2s$ mode of active slip
import pystokes, numpy as np, matplotlib.pyplot as plt

# particle radius, self-propulsion speed, number and fluid viscosity
b, eta, Np = 1.0, 1.0/6.0, 1

# initialize
r, p = np.array([0.0, 0.0, 3.4]), np.array([0.0, 1.0, 0])
V2s  = pystokes.utils.irreducibleTensors(2, p)

# space dimension , extent , discretization
dim, L, Ng = 3, 10, 64;

# instantiate the Flow class
flow = pystokes.wallBounded.Flow(radius=b, particles=Np, viscosity=eta, gridpoints=Ng*Ng)

# create grid, evaluate flow and plot
rr, vv = pystokes.utils.gridYZ(dim, L, Ng)
flow.flowField2s(vv, rr, r, V2s)  
pystokes.utils.plotStreamlinesYZsurf(vv, rr, r, offset=6-1, density=1.4, title='2s')
```

```Python
#Example 2: Phoretic field due to active surface flux of l=0 mode
import pystokes, numpy as np, matplotlib.pyplot as plt
# particle radius, fluid viscosity, and number of particles
b, eta, Np = 1.0, 1.0/6.0, 1

#initialise
r, p = np.array([0.0, 0.0, 5]), np.array([0.0, 0.0, 1])
J0 = np.ones(Np)  # strength of chemical monopolar flux

# space dimension , extent , discretization
dim, L, Ng = 3, 10, 64;

# instantiate the Flow class
phoreticField = pystokes.phoretic.unbounded.Field(radius=b, particles=Np, phoreticConstant=eta, gridpoints=Ng*Ng)

# create grid, evaluate phoretic field and plot
rr, vv = pystokes.utils.gridYZ(dim, L, Ng)
phoreticField.phoreticField0(vv, rr, r, J0)  
pystokes.utils.plotContoursYZ(vv, rr, r, density=.8, offset=1e-16,  title='l=0') 
```

Other examples include
* [Irreducible Active flows](https://github.com/rajeshrinet/pystokes/blob/master/examples/ex1-unboundedFlow.ipynb)
* [Effect of plane boundaries on active flows](https://github.com/rajeshrinet/pystokes/blob/master/examples/ex2-flowPlaneSurface.ipynb)
* [Active Brownian Hydrodynamics near a plane wall](https://github.com/rajeshrinet/pystokes/blob/master/examples/ex3-crystalNucleation.ipynb)
* [Flow-induced phase separation at a plane surface](https://github.com/rajeshrinet/pystokes/blob/master/examples/ex4-crystallization.ipynb)
* [Irreducible autophoretic fields](https://github.com/rajeshrinet/pystokes/blob/master/examples/ex5-phoreticField.ipynb)
* [Autophoretic arrest of flow-induced phase separation](https://github.com/rajeshrinet/pystokes/blob/master/examples/ex6-arrestedCluster.ipynb)


## Selected publications

* [PyStokes: phoresis and Stokesian hydrodynamics in Python](https://doi.org/10.21105/joss.02318), R Singh and R Adhikari, **Journal of Open Source Software**, 5(50), 2318, (2020). *(Please cite this paper if you use PyStokes in your research)*.

* [Controlled optofluidic crystallization of colloids tethered at interfaces](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.068001), A Caciagli, R Singh, D Joshi, R Adhikari, and E Eiser, **Physical Review Letters** 125 (6), 068001 (2020)

* [Periodic Orbits of Active Particles Induced by Hydrodynamic Monopoles](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.088003), A Bolitho, R Singh, R Adhikari, **Physical Review Letters** 124 (8), 088003 (2020)

* [Competing phoretic and hydrodynamic interactions in autophoretic colloidal suspensions](https://aip.scitation.org/doi/full/10.1063/1.5090179), R Singh, R Adhikari, and ME Cates, **The Journal of Chemical Physics** 151, 044901 (2019)

* [Generalized Stokes laws for active colloids and their applications](https://iopscience.iop.org/article/10.1088/2399-6528/aaab0d), R Singh and R Adhikari, **Journal of Physics Communications**, 2, 025025 (2018)


* [Flow-induced phase separation of active particles is controlled by boundary conditions](https://www.pnas.org/content/115/21/5403), S Thutupalli, D Geyer, R Singh, R Adhikari, and HA Stone, **Proceedings of the National Academy of Sciences**, 115, 5403 (2018)  

* [Universal hydrodynamic mechanisms for crystallization in active colloidal suspensions](https://doi.org/10.1103/PhysRevLett.117.228002), R Singh and R Adhikari,  **Physical Review Letters**, 117, 228002 (2016)

See full publication list [here](https://github.com/rajeshrinet/pystokes/wiki/Publications).

## Support

* For help with and questions about PyStokes, please post to the [pystokes-users](https://groups.google.com/forum/#!forum/pystokes) group.
* For bug reports and feature requests, please use the [issue tracker](https://github.com/rajeshrinet/pystokes/issues) on GitHub.

## License
We believe that openness and sharing improves the practice of science and increases the reach of its benefits. This code is released under the [MIT license](http://opensource.org/licenses/MIT). Our choice is guided by the excellent article on [Licensing for the scientist-programmer](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002598). 
		
CONTRIBUTING
============

When contributing code to PyStokes, you agree to release it under the project's license. Please first discuss the change you wish to make via issue, email, or any other method with the owners. 

For help with and questions about PyStokes, please post to the [pystokes-users](https://groups.google.com/forum/#!forum/pystokes) group.
Please follow the [Contributor Covenant](https://www.contributor-covenant.org/version/2/0/code_of_conduct/) in all PyStokes fora. Thank you!
The full set of examples is at: https://github.com/rajeshrinet/pystokes/tree/master/examples
This folder contains the core files of the PyStokes library.

Each file computes the rigid body motion (RBM) and flow given the colloidal configuration.
The filenames, as described below, correspond to the boundary conditions in the flow. The boundary conditions are implemented using an appropriate Green's function of Stokes flow. See details at: https://arxiv.org/abs/1910.00909

The main files are:
* unbounded.pyx - the fluid flow only vanishes at infinity. Implemented using Oseen tensor.
* periodic.pyx - Implemented using the Ewald sum of Oseen tensor. The box is of size L.
* interface.pyx - the normal component of the flow vanishes at the interface (z=0). Implemented using Blake, J. Biomech. 8 179 (1975).
* wallBounded.pyx - fluid flow near a plane wall at z=0. The flow vanishes (no-slip boundary condition) at the plane wall. That region of interest is the upper half space, z>0. Implemented using Lorentz-Blake tensor. See Blake, Proc. Camb. Phil. Soc. 70 303 (2017).
* twoWalls.pyx - two no-slip walls are at z=0 and z=H. The physical region is then `0<z<H'. Implemented using the approximate solution of Liron and Mochon. See Liron and Mochon, J. Eng. Math, 10 143 (1976).
* utils.pyx - has miscellaneous functionalities


## Force fields
* forceFields.pyx has implementation of various body forces and torques acting in colloidal systems.


## Phoretic fields 
The folder `phoretic` contains files to study phoresis of active particles. 

---
Corresponding to each .pyx file, there is a .pxd file. [A .pxd file](https://cython.readthedocs.io/en/latest/src/tutorial/pxd_files.html) contains declaration of cdef classes, methods, etc.  
This folder contains file to study phoresis of active particles. 

Each file computed the phoretic motion (phoresis) and phoretic field given the colloidal configuration.
The word phoresis derives from Greek word *phoras* which means bearing or "to carry". 

The filenames, as described below, correspond to the boundary conditions in the flow. The boundary conditions are implemented using an appropriate Green's function of the Laplace equation. See details of theory at: [J. Chem. Phys. 151, 044901 (2019)](https://aip.scitation.org/doi/abs/10.1063/1.5090179)

* unbounded.pyx - the phoretic flux vanishes at infinity. 
* wallBounded.pyx - the phoretic flux vanishes at a plane wall, which is located at z=0, such that region of interest is the upper half space, z>0. 

---
Corresponding to each .pyx file, there is a .pxd file. [A .pxd file](https://cython.readthedocs.io/en/latest/src/tutorial/pxd_files.html) contains declaration of cdef classes, methods, etc.  
## tests 
The folder tests contain unit testing framework for PyStokes







This folder contains codes to simulate colloids in periodic geometry of Stokes flow. Examples include

* Sedimentation: Codes for classical case of sedimentation
* Self-Propulsion: Codes for calculating the renormalized self-propulsion velocity of an active colloids in periodic domain
Examples of force fields which can be implemented in PyStokes. 
# PyStokes changelogs

## 2.2.5
* add a function for MSD in 2D in utils 

## 2.2.4
* minor fix to the two wall code

## 2.2.3
* Made a submodule pystokes.phoretic.unbounded etc to segregate code for phoresis

## 2.2.1
* Add more integrators and update simulate interface

## 2.2.0
* Update version after JOSS review

## 2.1.9 
* Remove MIMA from PyStokes

## 2.1.8 
* Update tests and utils so that it contnues python 2 support 

## 2.1.5 
* Add pyproject.toml

## 2.1.1 
* Second release of PyStokes

## 2.0.0 
* PyStokes has been merged with all of PyLaplace and PyForces

## 1.0.0 
* Add files from GitLab (where it was developed) for PyStokes, PyLaplace, and PyForces. 

## 0.1 
* First pre-release of PyStokes 
---
title: 'PyStokes: phoresis and Stokesian hydrodynamics in Python'
tags:
  - Python
  - Cython
  - Active particles
  - Stokes flow
  - Hydrodynamic interactions
  - Phoretic interactions
authors:
  - name: Rajesh Singh
    orcid: 0000-0003-0266-9691
    affiliation: 1
  - name: Ronojoy Adhikari
    affiliation: 1
affiliations:
 - name: DAMTP, Centre for Mathematical Sciences, University of Cambridge, Wilberforce Road, Cambridge CB3 0WA, UK
   index: 1
date: 6 June 2020
bibliography: paper.bib
---

# Summary

PyStokes is a Python library for studying phoretic and hydrodynamic interactions between spherical particles when these interactions can be described by the solutions of, respectively, the Laplace and Stokes equations. The library has been specifically designed for studying these interactions in suspensions of active particles,  which are distinguished by their ability to produce flow, and thus motion, in the absence of external forces or torques. Such particles are endowed with a mechanism to produce hydrodynamic flow in a thin interfacial layer, which may be due to the motion of cilia, as in microorganisms [@brennen1977] or osmotic flows of various kinds in response to spontaneously generated gradients of phoretic fields [@ebbens2010pursuit]. The latter, often called autophoresis,  is a generalisation of well-known phoretic phenomena including, *inter alia*,  electrophoresis (electric field), diffusiophoresis (chemical field) and thermophoresis (temperature field) that occur in response to externally imposed gradients of phoretic fields [@anderson1989colloid]. 

![Input and output structure of PyStokes to determine the hydrodynamic and phoretic interactions between active particles in a three-dimensional domain $V$. The equations are coupled by active boundary conditions on the surface $S_{i}$ of the particles. Particle indices are $i=1,\ldots,N$ and harmonic indices are $l=1,2,\ldots$ and $\sigma=s,a,t$ (see text). \label{fig:figS}](FigSchema.png)

Hydrodynamic and phoretic interactions between active particles in a viscous fluid are central to the understanding of their collective dynamics.  Under experimentally relevant conditions, the motion of the fluid is governed by the Stokes equation and that of the phoretic field, if one is present, by the Laplace equation.  The “activity” appears in these equations as boundary conditions on the particle surfaces that prescribe the slip velocity in the Stokes equation and flux of the phoretic field in the Laplace equation (see \autoref{fig:figS}). 
The slip velocity and the phoretic flux are related by a linear constitutive law that can be derived from a detailed analysis of the boundary layer physics [@anderson1989colloid]. The Stokes and Laplace equations are coupled by this linear constitutive law only at the particle boundaries. The linearity of the governing equations and the coupling boundary conditions allows for a formally exact solution of the problem of determining the force per unit area on the particle surfaces. This formally exact solution can be approximated to any desired degree of accuracy by a truncated series expansion in a complete basis of functions on the particle boundaries. This, in turn, leads to an efficient and accurate numerical method for computing hydrodynamic and phoretic interactions between active particles.

In addition to the joint computation of phoretic and hydrodynamic interactions, the  PyStokes library can be used to compute the hydrodynamically interacting motion of squirming particles where the slip is specified independently of a phoretic field, or the dynamics of passive suspensions where the slip vanishes and forces and torques are prescribed. The PyStokes library can also compute hydrodynamically correlated  Brownian motion, and thus, allows the study of the interplay between passive, active, and Brownian contributions to motion. 

The PyStokes library has been used to model suspensions of microorganisms [@bolitho2020; @singh2016crystallization], 
synthetic autophoretic particles [@singh2016crystallization; @singh2019competing] and self-propelling droplets [@thutupalli2018FIPS]. Our software implementation uses a polyglot programming approach that combines the readability of Python with the speed of Cython and retains the advantages of a high-level, dynamically typed, interpreted language without sacrificing performance.

# Methods

Our method relies on the reduction of a linear elliptic partial differential equation (PDE) to systems of linear algebraic equations using the following steps:   

![Key mathematical steps underpinning the PyStokes codebase.\label{fig:example}](figure.png)

The *first* step is the representation of the solution of an elliptic PDE in a three-dimensional volume $V$ as an integral over the 
boundary of the surface $S$ [@fkg1930bandwertaufgaben; @ladyzhenskaya1969; @youngren1975stokes; @zick1982stokes; @pozrikidis1992; @muldowney1995spectral; @cheng2005heritage; @singh2015many]. 
For the Laplace equation, this is the classical theorem of Green [@jackson1962classical]; for the Stokes equation, it is the generalization obtained by 
Lorentz [@fkg1930bandwertaufgaben; @lorentz1896eene; @ladyzhenskaya1969]. The integral representation leads to a linear integral equation that provides a 
functional relation between the field and its flux on $S$. Thus, if the surface flux in the Laplace equation is specified, the surface concentration is determined by the solution of the Laplace boundary integral equation. Similarly, if the surface velocity in the Stokes equation is specified, the surface traction is determined by the solution of the Stokes boundary integral equation. This transformation of the governing PDE is the most direct way of relating boundary conditions (surface flux, slip velocities) 
to boundary values (surface concentration, surface traction). It reduces the dimensionality of the problem from a three-dimensional one in $V$ to a two-dimensional one on $S$. 
The *second* step is the spectral expansion of the field and its flux in terms of global basis functions on $S$. We use the geometry-adapted tensorial spherical harmonics, 
which provide a unified way of expanding both scalar and vector quantities on the surface of a sphere. These functions are both complete and orthogonal and provide representations of the three-dimensional rotation group [@hess2015tensors]. Thus, symmetries of the active boundary conditions can be represented straightforwardly and transparently. 
The *third* step is the discretization of the integral equation using the procedure of Ritz and Galerkin [@boyd2001chebyshev; @finlayson1966method], which reduces it to an infinite-dimensional self-adjoint linear system in the expansion coefficients. This exploits the orthogonality of the basis functions on the sphere. The matrix elements of the linear system can be evaluated analytically in terms of the Green's functions of the respective elliptic equations. The *fourth* step is the truncation of the infinite-dimensional linear system to a finite-dimensional one that can be solved by standard methods of linear algebra adapted for self-adjoint systems [@saad2003iterative]. 
The analytical solution can be obtained by Jacobi iteration, which is equivalent to Smoluchowski's method of reflection. Numerical solutions can be obtained by the conjugate gradient method, at a cost quadratic in the number of unknowns. From this solution, we can reconstruct the field and the flux on the boundary, use these to determine the fields in the bulk, and from there, compute derived quantities. 

The above steps have been elaborated in several 
papers [@singh2016crystallization; @singh2017fluctuation; @singh2018generalized; @singh2015many; @singh2019competing] and we do not repeat them in detail here. Briefly, the expansion coefficients of the slip can be either specified or obtained as a solution of Laplace equation. Once the coefficients are determined, the following equation are solved numerically to obtain velocity and angular velocity of the $i$-th particle
\begin{align*}
{\mathbf{V}}_{i}&=
\boldsymbol{\mu}_{ij}^{TT}\cdot{\mathbf{F}_{j}^{P}}+ 
\boldsymbol{\mu}_{ij}^{TR}\cdot{\mathbf{T}_{j}^{P}} + 
\sum_{l\sigma=1s}^{\infty}\boldsymbol{\pi}_{ij}^{(T,\,l\sigma)}\cdot{\mathbf{V}_{j}^{(l\sigma)}},\qquad \boldsymbol{\mu}_{ij}^{\alpha\beta}\quad:\text{ mobility matrices},
\\
{\mathbf{\Omega}}_{i}&=\underbrace{
\boldsymbol{\mu}_{ij}^{RT}\cdot\mathbf{F}_{j}^{P}+
\boldsymbol{\mu}_{ij}^{RR}\cdot\mathbf{T}_{j}^{P}}_{\text{Passive}}+
\sum_{l\sigma=1s}^{\infty}\underbrace{
\boldsymbol{\pi}_{ij}^{(R,\,l\sigma)}\cdot\mathbf{{\mathbf{V}}}_{j}^{(l\sigma)}}_{\text{Active}},\qquad \boldsymbol{\pi}_{ij}^{(\alpha, l\sigma)}:\text{ propulsion tensors}.
\end{align*}
In the above, $\alpha,\beta=(T,R)$, repeated particle index $j$ is summed, and
\begin{equation*}
\mathbf{F}_{i}^{P}:\text{body force},\quad \mathbf{T}_{i}^{P}:\text{body torque},\quad \mathbf{V}_{i}^{(l\sigma)}: l\sigma\text{-th expansion coefficients of active slip}.
\end{equation*}
Thus, hydrodynamic interactions between particles with no-slip boundary conditions can be computed entirely in terms of mobility matrices, as implemented in existing numerical libraries [@hinsen1995; @libstokes], to study suspensions of passive particles. The active contributions due to the slip boundary condition is given in terms of propulsion tensors [@singh2015many]. To the best of our knowledge, PyStokes is the only numerical implementation of propulsion tensors to model suspensions of active particles.      


To summarize, the principal features that set our method apart are (a) the restriction of independent fluid and phoretic degrees of freedom to the particle boundaries 
(b) the freedom from grids, both in the bulk of the fluid and on the particle boundaries and 
(c) the ability to handle, within the same numerical framework, a wide variety of geometries and boundary conditions, 
including unbounded volumes, volumes bounded by plane walls or interfaces, periodic volumes and, indeed, 
any geometry-boundary condition combination for which the Green's functions of the governing equations are simply evaluated.



The PyStokes library can be instantiated in the following way to 

* obtain phoretic field created by active particles at a given set of points 
```python
phoreticField = pystokes.phoreticUnbounded.Field(radius=1, particles=1, 
                phoreticConstant=1, gridpoints=4096)
```

* evaluate fluid flow created by active particles at a given set of points
```python
Flow = pystokes.unbounded.Flow(radius=1, particles=1, viscosity=1, 
        gridpoints=4096) 
```

* determine phoretic field at surface of active particles
```python
phoresis = pystokes.phoreticUnbounded.Phoresis(radius=1, particles=1024, 
            phoreticConstant=1)
```

* compute rigid body motion of hydrodynamically interacting particles
```python
Rbm = pystokes.unbounded.Rbm(radius=1, particles=1024, viscosity=1) 
```


The above instantiation can then be used to compute `Flow` and `Rbm` due to body forces, body torques, and each irreducible mode of the surface slip in various geometries of Stokes flow, by replacing `unbounded` with `wallBounded`, `periodic`, etc. `pystokes.forceFields` contains an implementation of force fields commonly used in colloidal systems for completeness. The arXiv preprint [@singh2019Hydrodynamic] of this article contains more detailed documentation and examples.
  

# Acknowledgements

We thank A Banerjee, ME Cates, S Date, A Donev, E Eiser, D Frenkel, S Ghose, R
Goldstein, J Hinch, A Laskar, AJC Ladd, RK Manna, I Pagonabarraga, DJ Pine, T
Pradeep, R Simon, HA Stone, G Subramanian, PB Sunil Kumar, and S Thutupalli 
for useful discussions; Rajeev Singh and Abhrajit Laskar for code contributions in the initial stages of development. 
This work was funded in parts by the European Research Council under the EU’s Horizon 2020 Program, Grant No. 740269; 
a Royal Society-SERB Newton International Fellowship to RS; and an Early Career Grant to RA from the Isaac Newton Trust.

# References
Once all the .rst files are updated following [autodoc](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html), run the following

```bash
pip install -U Sphinx
pip install sphinx_rtd_theme
make html
```
Phoresis in an unbounded domain
==================================

Phoretic motion in an infinite space 

Phoresis
----------------------------------------
.. autoclass:: pystokes.phoretic.unbounded.Phoresis
    :members:


Field
--------------------------
.. autoclass:: pystokes.phoretic.unbounded.Field
    :members:

Stokes flow in an unbounded domain
==================================

Hydrodynamic interactions of active particles in an unbounded domain. 

RBM (Rigid Body Motion: Velocity and Angular Velocity)
------------------------------------------------------
.. autoclass:: pystokes.unbounded.Rbm
    :members:


Flow
--------------------------
.. autoclass:: pystokes.unbounded.Flow
    :members:

Utils: Miscellaneous function
==================================

.. automethod:: pystokes.utils.irreducibleTensors

.. automethod:: pystokes.utils.gridXY

.. automethod:: pystokes.utils.gridYZ

.. automethod:: pystokes.utils.simulate
Phoresis in a space bounded by a wall
===============================================

Phoresis in the half-space bounded by a plane no-slip wall 

Phoresis
----------------------------------------
.. autoclass:: pystokes.phoretic.wallBounded.Phoresis
    :members:


Field
--------------------------
.. autoclass:: pystokes.phoretic.wallBounded.Field
    :members:

Stokes flow near interfaces
==================================================

Hydrodynamic interactions of active particles in the half-space bounded by a plane no-shear surface

RBM (Rigid Body Motion: Velocity and Angular Velocity)
------------------------------------------------------
.. autoclass:: pystokes.interface.Rbm
    :members:


Flow
--------------------------
.. autoclass:: pystokes.interface.Flow
    :members:

Force fields in colloidal systems
==================================

Force fields in colloidal systems

Forces 
----------------------------------------
.. autoclass:: pystokes.forceFields.Forces
    :members:

Torques 
----------------------------------------
.. autoclass:: pystokes.forceFields.Torques
    :members:

Stokes flow in a Hele-Shaw cell
===========================================

Stokesian hydrodynamic in a Hele-Shaw cell 

RBM (Rigid Body Motion: Velocity and Angular Velocity)
------------------------------------------------------
.. autoclass:: pystokes.twoWalls.Rbm
    :members:


Flow
--------------------------
.. autoclass:: pystokes.twoWalls.Flow
    :members:

Stokes flow in a periodic domain 
================================

Hydrodynamic interactions of active particles in a periodic 3D space 

RBM (Rigid Body Motion: Velocity and Angular Velocity)
------------------------------------------------------
.. autoclass:: pystokes.periodic.Rbm
    :members:


Flow
--------------------------
.. autoclass:: pystokes.periodic.Flow
    :members:

Stokes flow bounded by a wall 
=================================================

Hydrodynamic interactions of active particles in the half-space bounded by a plane no-slip surface

RBM (Rigid Body Motion: Velocity and Angular Velocity)
------------------------------------------------------
.. autoclass:: pystokes.wallBounded.Rbm
    :members:


Flow
--------------------------
.. autoclass:: pystokes.wallBounded.Flow
    :members:

PyStokes API
==================================
.. image:: ../../examples/banner.png
  :width: 700
  :alt: PyRoss banner



`PyStokes <https://github.com/rajeshrinet/pystokes>`_ is a numerical library for Phoresis and Stokesian hydrodynamics in Python. It uses a grid-free method, combining the integral representation of Stokes and Laplace equations, spectral expansion, and Galerkin discretization, to compute phoretic and hydrodynamic interactions between spheres with slip boundary conditions on their surfaces. The library also computes suspension scale quantities, such as rheological response, energy dissipation and fluid flow. The computational cost is quadratic in the number of particles and upto 1e5 particles have been accommodated on multicore computers. The library has been used to model suspensions of **microorganisms**,  **synthetic autophoretic particles** and **self-propelling droplets**. 


Please see installation instructions and more details in the `README.md <https://github.com/rajeshrinet/pystokes/blob/master/README.md>`_ on GitHub. 




API Reference
=============

.. toctree::
   :maxdepth: 1

   
   unbounded
   interface
   wallBounded
   twoWalls
   periodic
   phoreticUnbounded
   phoreticWallBounded
   forceFields
   utils
