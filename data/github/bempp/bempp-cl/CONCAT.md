cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
- family-names: "Betcke"
  given-names: "Timo"
  orcid: "https://orcid.org/0000-0002-3323-2110"
- family-names: "Scroggs"
  given-names: "Matthew W."
  orcid: "https://orcid.org/0000-0002-4658-2443"
title: "Bempp-cl"
version: 0.2.4
date-released: 2021-03-18
url: "https://github.com/bempp/bempp-cl"
preferred-citation:
  type: article
  authors:
  - family-names: "Betcke"
    given-names: "Timo"
    orcid: "https://orcid.org/0000-0002-3323-2110"
  - family-names: "Scroggs"
    given-names: "Matthew W."
    orcid: "https://orcid.org/0000-0002-4658-2443"
  doi: "10.21105/joss.02879"
  journal: "Journal of Open Source Software"
  title: "Bempp-cl: A fast Python based just-in-time compiling boundary element library"
  volume: 6
  issue: 59
  start: 2879
  month: 3
  year: 2021
# Bempp-cl
[![Documentation Status](https://readthedocs.org/projects/bempp-cl/badge/?version=latest)](https://bempp-cl.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02879/status.svg)](https://doi.org/10.21105/joss.02879)

Bempp-cl is an open-source boundary element method library that can be used to assemble all the standard integral kernels for
Laplace, Helmholtz, modified Helmholtz, and Maxwell problems. The library has a user-friendly Python interface that allows the
user to use BEM to solve a variety of problems, including problems in electrostatics, acoustics and electromagnetics.

Bempp-cl began life as BEM++, and was a Python library with a C++ computational core. The ++ slowly changed into pp as
functionality gradually moved from C++ to Python with only a few core routines remaining in C++. Bempp-cl is the culmination
of efforts to fully move to Python. It is an almost complete rewrite of Bempp: the C++ core has been replaced by highly SIMD
optimised just-in-time compiled OpenCL kernels, or alternatively, by just-in-time compiled Numba routines, which are
automatically used on systems that do not provide OpenCL drivers. User visible functionality is strictly separated from the
implementation of computational routines, making it easy to add other discretisation technologies in the future (e.g. future
support for SYCL-based heterogeneous compute devices).

## Installation
Bempp-cl can be installed from this repository by running:
```bash
python setup.py install
```

Full installation instuctions, including installation of dependencies, can be found at
[bempp.com/installation.html](https://bempp.com/installation.html).

## Documentation
Full documentation of Bempp can be found at [bempp.com/documentation](https://bempp.com/documentation/index.html)
and in [the Bempp Handbook](https://bempp.com/handbook). Automatically generated documentation of the Python API
can be found on [Read the Docs](https://bempp-cl.readthedocs.io/en/latest/).

## Testing
The functionality of the library can be tested by running:
```bash
python -m pytest test/unit
```
Larger validation tests that compare the output with the previous version of Bempp can be run with:
```bash
python -m pytest test/validation
```

## Getting help
Errors in the library should be added to the [GitHub issue tracker](https://github.com/bempp/bempp-cl/issues).

Questions about the library and its use can be asked on the [Bempp Discourse](https://bempp.discourse.group).

## Licence
Bempp-cl is licensed under an MIT licence. Full text of the licence can be found [here](LICENSE.md).

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
bempp@mscroggs.co.uk (Matthew Scroggs) or tbetcke@ucl.ac.uk (Timo Betcke).
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
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations

# How to contribute

## Reporting bugs
If you find a bug in Bempp-cl, please report it on the [issue tracker](https://github.com/bempp/bempp-cl/issues)
and tag it "bug". If you are unsure if what you have found is a bug, or want to ask a more general question,
you can do this on the [Bempp discourse forum](https://bempp.discourse.group/).

## Suggesting enhancements
If you want to suggest a new feature or an improvement of a current feature, you can submit this
on the [issue tracker](https://github.com/bempp/bempp-cl/issues) and tag it "enhancement".

## Submitting a pull request
If you want to directly submit code to Bempp-cl, you can do this by forking the Bempp-cl repo, then submitting a pull request.
If you want to contribute, but are unsure where to start, have a look at the
[issue tracker](https://github.com/bempp/bempp-cl/issues) for issues labelled "good first issue".

On opening a pull request, unit tests and flake8 style checks will run. You can click on these in the pull request
to see where (if anywhere) there are errors in your code.

If you plan to work on a large feature, it would be sensible to suggest it on the issue tracker or Discource first so we
can confirm that is a feature we would like before you spend too long working on it.

## Code of conduct
We expect all our contributors to follow [the Contributor Covenant](CODE_OF_CONDUCT.md). Any unacceptable
behaviour can be reported to Matthew (bempp@mscroggs.co.uk) or Timo (tbetcke@ucl.ac.uk).

----

Thank-you for considering contributing to Bempp!
---
title: 'Bempp-cl: A fast Python based just-in-time compiling boundary element library.'
tags:
  - Python
  - OpenCL
  - boundary element method
  - partial differential equations
  - numerical analysis
authors:
  - name: Timo Betcke
    orcid: 0000-0002-3323-2110
    affiliation: 1
  - name: Matthew W. Scroggs
    orcid: 0000-0002-4658-2443
    affiliation: 2
affiliations:
 - name: Department of Mathematics, University College London
   index: 1
 - name: Department of Engineering, University of Cambridge
   index: 2
date: 15 January 2020
bibliography: paper.bib
---

# Summary
The boundary element method (BEM) is a numerical method for approximating the solution of certain types of partial 
differential equations (PDEs) in homogeneous bounded or unbounded domains. The method finds an approximation by discretising 
a boundary integral equation that can be derived from the PDE. The mathematical background of BEM is covered in, for example, 
@Stein07 or @McLean. Typical applications of BEM include electrostatic problems, and acoustic and electromagnetic scattering.

Bempp-cl is an open-source boundary element method library that can be used to assemble all the standard integral kernels for
Laplace, Helmholtz, modified Helmholtz, and Maxwell problems. The library has a user-friendly Python interface that allows the
user to use BEM to solve a variety of problems, including problems in electrostatics, acoustics and electromagnetics.

Bempp-cl began life as BEM++, and was a Python library with a C++ computational core. The ++ slowly changed into pp as 
functionality gradually moved from C++ to Python with only a few core routines remaining in C++. Bempp-cl is the culmination 
of efforts to fully move to Python, and is an almost complete rewrite of Bempp.

For each of the applications mentioned above, the boundary element method involves approximating the solution of a partial
differential equation (Laplace's equation, the Helmholtz equation, and Maxwell's equations respectively) by writing the problem
in boundary integral form, then discretising. For example, we could calculate the scattered field due to an electromagnetic wave
colliding with a series of screens by solving
\begin{align*}
\nabla\times\nabla\times \mathbf{E} -k^2 \mathbf{E} &= 0,\\
\boldsymbol{\nu}\times\mathbf{E}&=0\text{ on the screens},
\end{align*}
where $\mathbf{E}$ is the sum of a scattered field $\mathbf{E}^\text{s}$ and an incident field $\mathbf{E}^\text{inc}$,
and $\boldsymbol{\nu}$ is the direction normal to the screen. (Additionally, we must impose the Silver--Müller radiation condition
to ensure that the problem has a unique solution.) This problem is solved, and the full method is derived,
in one of the tutorials available on the Bempp website [@Bempp-maxwell-example]. The solution to this problem is shown below.

![An electromagnetic wave scattering off three screens.](maxwell_sol.png){ width=50% }

# Statement of need
Bempp-cl provides a comprehensive collection of routines for the assembly of boundary integral operators to solve a wide
range of relevant application problems. It contains an operator algebra that allows a straight-forward implementation of
complex operator preconditioned systems of boundary integral equations [@operatoralg] and in particular implements
everything that is required for Calderón preconditioned Maxwell [@maxwellbempp] problems. Bempp-cl uses PyOpenCL [@pyopencl]
to just-in-time compile its computational kernels on a wide range of CPU and GPU devices and modern architectures. Alternatively,
a fallback Numba implementation is provided.

OpenCL is used as it is able to compile C-based kernels to run on a wide range of CPU and GPU devices, without the need to
write device specific code. Numba is offered as an alternative as it is easily available on all platforms and provides a
version of the library that is significantly faster than using pure Python.

Bempp-cl is aimed at those interested in using boundary element method to solve problems, particularly those from a mathematical background.
The syntax of the library is designed to closely resemble the boundary integral representation of the problem being solved, making
the implementation of a problem simple once this representation is known.

There are only a small number of alternative boundary element method softwares available.
The most popular is BETL [@BETL], a C++ template library that is available for free for academic use only.
As a Python library, Bempp-cl is easier to interface with other popular libraries with Python interfaces---for example,
is can be used alongside the finite element method library FEniCS [@fenicsbook] to solve coupled finite and boundary element
problems [@fembemexample].
Bempp-cl also benefits from being fully open source library and available under an MIT license.
A number of other libraries exist designed for specific applications, such as PyGBe for biomolecular electrostatics [@PyGBe]
and abem for acoustics [@abem]. Bempp-cl can be used for a much wider range of problems than these specialised libraries.


# An overview of Bempp features
Bempp-cl is divided into two parts: `bempp.api` and `bempp.core`.
The user interface of the library is contained in `bempp.api`.
The core assembly routines of the library are contained in `bempp.core`. The majority of users of Bempp-cl are unlikely to need
to directly interact with the functionality in `bempp.core`.

There are five main steps that are commonly taken when solving a problem with BEM:

1. First a surface grid (or mesh) must be created on which to solve the problem.
2. Finite dimensional function spaces are defined on this grid.
3. Boundary integral operators acting on the function spaces are defined.
4. The operators are discretised and the resulting linear systems solved.
5. Domain potentials and far field operators can be used to evaluate the solution away from the boundary.

Bempp-cl provides routines that implement each of these steps.

## Grid Interface
The submodule `bempp.api.shapes` contains the definitions of a number of shapes. From these, grids with various element sizes 
can be created internally using Gmsh [@gmsh]. Alternatively, meshes can be imported from many formats using the meshio library 
[@meshio]. Bempp-cl currently only supports flat triangle based surface meshes. Higher-order triangular meshes may be 
supported in the future.

## Function Spaces
Bempp-cl provides piecewise constant and piecewise linear (continuous and discontinuous) function spaces for solving scalar 
problems. For Maxwell problems, Bempp-cl can create Rao--Wilton--Glisson [@rwg] div-conforming spaces and Nédélec 
[@nedelec] curl-conforming spaces. In addition to these, Bempp-cl can also generate constant and linear spaces on the 
barycentric dual grid as well as Buffa--Christiansen div-conforming spaces, as described in @bc. These spaces can all be 
created using the `bempp.api.function_space` command.

## Boundary operators
Boundary operators for Laplace, Helmholtz, modified Helmholtz and Maxwell problems can be found in the 
`bempp.api.operators.boundary` submodule, as well as sparse identity operators. For Laplace and Helmholtz problems, Bempp-cl 
can create single layer, double layer, adjoint double layer and hypersingular operators. For Maxwell problems, both electric 
field and magnetic field operators can be used.

## Discretisation and solvers
Operators are assembled using OpenCL or Numba based dense assembly, or via interface to fast multipole methods.
Internally, Bempp-cl uses PyOpenCL [@pyopencl] to just-in-time compile its operator assembly routines on a wide range of CPU
and GPU compute devices. On systems without OpenCL support, Numba [@numba] is used to just-in-time compile
Python-based assembly kernels, giving a slower but still viable alternative to OpenCL.

Bempp-cl provides an interface to the Exafmm-t library [@exafmm] for faster assembly of larger problems with lower memory
requirements using the fast multipole method (FMM). The interface to Exafmm-t is written in a generic way so that other
FMM libraries or alternative matrix compression techniques could be used in future. 

The submodule `bempp.api.linalg` contains wrapped versions of SciPy's [@scipy] LU, CG, and GMRes solvers. By using 
SciPy's `LinearOperator` interface, Bempp-cl's boundary operators can easily be used with other iterative solvers.

## Potential and far field operators
Potential and far field operators for the evaluation at points in the domain or the asymptotic behavior at infinity are 
included in the `bempp.api.operators.potential` and `bempp.api.operators.far_field` submodules.

## Further information
Full documentation of the library, including a number of example Jupyter notebooks, can be found online at ``bempp.com`` and
in the in-development Bempp Handbook [@bempphandbook].

# Acknowledgements
We would like to thank the Exafmm team [@exafmm], and here in particular Lorena Barba and Tingyu Wang for their efforts to 
integrate Exafmm-t into Bempp-cl. We further thank the HyENA team [@hyena] at Graz University of Technology who provided C++ 
definitions of core numerical quadrature rules, which were translated to Python as part of the development effort for 
Bempp-cl.
    
# References
