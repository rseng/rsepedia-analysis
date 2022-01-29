# Architecture

This document describes the high-level architecture of libecpint, to help any new contributors find their way around the codebase.

## High-level overview

The code is roughly divided into four chunks:

- the interface,
- the fixed backend,
- the generated backend,
- the code generator

As a general rule, the API will not change. If functionality is added to the library, it should be exposed through ```api.hpp```, following the style of other functions in that interface. Changing any existing parts of the API will break compatibility, so we will only do this if it is *absolutely* necessary.

The other three parts are discussed below.

## Code Map

### Fixed backend

This is found in ```src/lib``` and can be further divided into the following functionalities:

- Quadrature (```gaussquad.cpp, radial_quad.cpp```)
- Analytic angular integrals (```angular.cpp```)
- Containers (```ecp.cpp, gshell.cpp```)
- Mathematical utilities (```bessel.cpp, mathutil.cpp```)

In brief, this is where you should put changes of any of the above kinds. These are all things that are used by the library as part of core functionality, and that do not change depending on how the generated part of the code is initialised.

Each file generally describes one class (for example, the ```ECP``` object in ```ecp.cpp``` or the ```AngularIntegral``` object in ```angular.cpp```) _OR_ one set of functionalities (e.g. adaptive quadrature for radial integrals, in ```radial_quad.cpp```). Please stick to this convention.

As a general rule, any changes in the fixed backend will only affect other parts of the fixed backend, and will not affect that API or generated code. The exception is ```ecpint.cpp``` which describes an ```ECPIntegral``` object, but connects the fixed and generated code together. As such it is also where much of the integral screening takes place.

### Generated backend

This is found in ```src/generated``` and ```src/generated/radial```. It is where any generated code, or templates required for code generation, are located. Currently there are two "part" files, for the generated file ```qgen.cpp``` that ends up in the main library.

The radial folder is where the unrolling of recurrence relations for the primitive radial integrals happens. Any algorithmic changes with respect to the recurrence relations should be found/placed in here.

### Code generator

The generated code is handled by ```src/generate.cpp``` and the utility functions in ```include/generate.hpp```. The former works by finding all the non-zero integral terms for a particular radial integral class using the ```SumTerm``` objects in the header. These are then sorted and made unique, before the radial code is written into the generated backend described above.   

This part is compiled and run as a separate project before the main library. Therefore if any changes you make depend on the generated code, and in particular the ```qgen``` array in ```ECPIntegral```, you should think very carefully about whether it  will be affected by any changes in the code generation. If so, functionality should be added to the generated backend, in the same way that ```qgen.cpp``` is handled. 
# Libecpint 1.0.7

[![Build Status](https://dev.azure.com/robertshaw383/libecpint/_apis/build/status/robashaw.libecpint?branchName=master)](https://dev.azure.com/robertshaw383/libecpint/_build/latest?definitionId=2&branchName=master)
[![codecov](https://codecov.io/gh/robashaw/libecpint/branch/master/graph/badge.svg)](https://codecov.io/gh/robashaw/libecpint)
[![Documentation Status](https://readthedocs.org/projects/libecpint/badge/?version=latest)](https://libecpint.readthedocs.io/en/latest/index.html)
[![Code Quality](https://www.code-inspector.com/project/15206/status/svg)]()

[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4694353.svg)](https://doi.org/10.5281/zenodo.4694353)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03039/status.svg)](https://doi.org/10.21105/joss.03039)

Libecpint is a C++ library for the efficient evaluation of integrals over ab initio effective core potentials, using a mixture of generated, recursive code and Gauss-Chebyshev quadrature. It is designed to be standalone and generic, and is now in its first stable release. If you experience any problems please raise an issue here; contributions and suggestions are also welcome.

## Contributing

Contributions are welcomed, either in the form of raising issues or pull requests on this repo. Please take a look at the Code of Conduct before interacting, which includes instructions for reporting any violations.

## New in first full release

- Analytical 1st and 2nd derivatives;
- Integration now >10x faster;
- New, high level API, with ECP library;
- Automated testing suite.

### Patch 1

- Bug fix in screening of on-ECP type 2 integrals
- Improvements in CMake build steps, thanks to nabbelbabbel/moritzBens

### Patch 2

- Fix for memory leaks in derivative routines
- Minor changes to CMake files

### Patch 3

- Fix bug in radial type 1 integrals where quadrature could fail to converge
- Const correctness throughout, should allow for parallelisation
- Minor updates to docs

### Patch 4

- Code generation now takes considerably less time and memory; MAX_L=8 takes ~35 seconds, peaking at 1.5GB of memory (joint effort with Thomas Dresselhaus and Peter Bygrave)
- This will be the final patch before v1.1

## Dependencies

- A modern C++ compiler, at least C++11 standard library is required. This has been tested with:
  * gcc (v6.3.0 and above)
  * clang (v10.0.0 and above), you may need the CXX flag "-std=c++14"
  * icpc (v20.2.1), may also need the CXX flag "-std=c++14"
- CMake/CTest build tools (v3.12 and higher)
- Python (2.7 or above, including 3 and higher)

Additionally, if you wish to regenerate the radial code (see below),  Python >=3.6 is required with numpy and sympy.

## Documentation

Please refer to the main documentation [here](https://libecpint.readthedocs.io/en/latest/index.html).

## Examples

There is also a working example in the example folder, with instructions of how to build and link against the library. Please also the API tests in tests/lib/

## Acknowledging usage

If you use this library in your program and find it helpful, that's great! Any feedback would be much appreciated. If you publish results using this library, please consider citing the following paper detailing the implementation:

R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017); doi: [10.1063/1.4986887](http://dx.doi.org/10.1063/1.4986887)

A full bibtex citation can be found in CITATION in the main directory.

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
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
robertshaw383 (at) gmail (dot) com.
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
# Working example

This directory contains a working example of how to use libecpint for its core functionality - calculating integrals and derivatives for, in this example, hydrogen iodide in a CC-pVTZ(-PP) basis.

## Build

To build this, make sure you have built libecpint (and its dependencies, pugixml and Faddeeva). Then using cmake in this directory

```
cmake . -Bbuild -DCMAKE_CXX_FLAGS="-I/include/dirs -std=c++14" -DCMAKE_EXE_LINKER_FLAGS="-L/lib/dirs"
```
where the -I and -L flags are used to point to the libecpint headers and library, respectively, if these are not already in your path.

To then build the library and run the test:
```
cd build
make
./example LIBECPINT_SHARE_DIR
```
where the argument points to the share/libecpint directory.

## With Eigen

If you have Eigen3 installed and want to test the example matrix build, simply add the flag
```
-D_WITH_EIGEN
```
to the CMAKE_CXX_FLAGS in the cmake step. 
---
title: 'libecpint: A C++ library for  the efficient evaluation of integrals over effective core potentials'
tags:
  - C++
  - computational chemistry
authors:
  - name: Robert A. Shaw
    orcid: 0000-0002-9977-0835
    affiliation: 1
  - name: J. Grant Hill
    orcid: 0000-0002-6457-5837
    affiliation: 1
affiliations:
 - name: Department of Chemistry, University of Sheffield, Sheffield S3 7HF, UK
   index: 1
date: February 2, 2021
bibliography: paper.bib
---

# Summary
Effective core potentials (ECPs) are widely-used in computational chemistry both to reduce the computational cost of calculations,[@Dolg2000] and include relevant physics that would not otherwise be present [@Dolg2002]. In particular, for heavy main-group atoms [@Wadt1985] and transition metals [@Hay1985], the number of core electrons greatly outnumbers the number of valence electrons. It is generally considered that these will not play a significant role in chemical reactivity, and thus can be frozen. Moreover, these electrons show significant relativistic character [@Dolg2002; @Dolg2012]. Both of these issues can be resolved with the introduction of an effective core, represented as a fixed electronic potential. This potential is typically represented as a linear combination of gaussians of varying angular momenta [@Dolg2000].

The introduction of an ECP results in an additional term in the core  Hamiltonian, over which new electronic integrals must be computed. These three-center integrals are far from trivial, and they cannot in general be treated the same way as other electronic integrals [@McMurchie1981; @FloresMoreno2006]. Several widely used computational chemistry codes lack the ability to calculate these integrals due to the difficulty involved in their computation. The present library, `libecpint`, provides an open-source solution to this. It is a standalone library written in modern C++ capable of the highly efficient computation of integrals over ECPs with gaussian orbitals of arbitrary angular momentum, along with their first and second geometric derivatives. The methods implemented are based on novel algorithms that use automatic code generation and symbolic simplification of recursive expressions, along with highly optimised Gauss-Chebyshev quadrature.

# Statement of need

Effective core potentials are an essential part of modern computational chemistry. However, existing implementations are typically unavailable or inaccessible for free use by the open source community. Commonly used proprietary software, such as Gaussian [@Gaussian16] or Molpro [@MOLPRO], do not make details of their implementations available, while the few open-source computational chemistry packages either do not include ECP functionality or use outdated implementations that would not be compatible with modern codebases. A notable example of this is the widely-used Psi4 package [@Psi4], in which a rudimentary version of `libecpint` was originally implemented. Prior to this, the inclusion of ECPs was one of the most requested features by the user base.

Additionally, there has been a recent renaissance in the development of efficient algorithms for evaluating ECP integrals. In particular, multiple research groups have outlined new approaches to prescreening integrals,[@Song2015; @Shaw2017; @McKenzie2018] greatly reducing the computational expense. The `libecpint` library implements many of these new algorithms, combining the recursive methods and fine-grained screening of Shaw et al. [@Shaw2017] with the higher-level screening of other recent work [@Song2015; @McKenzie2018]. The only known implementations of the latter papers are otherwise only available in proprietary software. Therefore `libecpint` represents a necessary contribution to the wider open-source computational chemistry community. It has already been adopted by multiple packages, including Entos QCore and Serenity [@Serenity], and will be part of a future release of Psi4 [@Psi4].

# Functionalities

The core functionality of `libecpint` is the evaluation of both type 1 and type 2 integrals over ECP integrals parametrised in terms of contracted sets of primitive gaussians, as described in Shaw et al. [@Shaw2017]. The component parts divide into the following functionalities:

- a built-in library of parametrised ECPs, with generic containers for Gaussian-type ECPs;
- a highly-optimised Bessel function evaluation routine;
- screening of ECP integrals over shell pairs of orbital basis functions [@Song2015; @McKenzie2018], across all ECPs in a system;
- fine-grained screening of the individual type 2 integrals over primitive gaussians; [@Shaw2017]
- recursive, automatically-generated radial integral code; [@Shaw2017]
- adaptive quadrature for integrals not covered by the recursive routines; [@FloresMoreno2006]
- first- and second-order geometric derivatives of ECP integrals over shell pairs of gaussians.

These features can be accessed via two levels of API:

- a high-level interface, where the user provides a molecular geometry and basis set, then calls routines that return the full tensor of integrals (or integral derivatives) across all shell pairs;
- a low-level interface, where the user provides parameters and deals with the primitive integrals directly.

This allows a great deal of flexibility for different use cases, potentially allowing for users to further develop or adapt the routines themselves.

All primitive integral routines have been designed to be thread safe, allowing users to readily parallelise their calculations. In addition, there is a built in testing and benchmarking suite, allowing for efficiency comparisons both with other codes, and when testing new algorithmic developments.  

# Acknowledgements

Thank you to Moritz Bensberg, Peter Bygraves, Thomas Dresselhaus, Christopher Junghans, Peter Kraus, Jan Unsleber, and Jens Wehner, for finding bugs and suggesting improvements to the initial release of `libecpint`, and often providing helpful solutions.

# References
