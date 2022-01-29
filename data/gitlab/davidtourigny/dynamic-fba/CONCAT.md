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
reported by contacting the project team at  USAdst2156@cumc.columbia.edu. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq

---
title: 'dfba: Software for efficient simulation of dynamic flux-balance analysis models in Python'
tags:
 - Mathematical optimization
 - Linear programming
 - Numerical simulation
 - Multi-scale metabolic modeling
 - C++
 - Python
authors:
 - name: David S. Tourigny
   orcid: 0000-0002-3987-8078
   affiliation: 1
 - name: Jorge Carrasco Muriel
   orcid: 0000-0001-7365-0299
   affiliation: 2
 - name: Moritz E. Beber
   orcid: 0000-0003-2406-1978
   affiliation: 2
affiliations:
 - name: Columbia University Irving Medical Center, 630 West 168th Street, New York, NY 10032 USA
   index: 1
 - name: Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark, Building 220, Kemitorvet, 2800 Kongens Lyngby, Denmark
   index: 2
date: 13 March 2020
bibliography: paper.bib
---

# Summary

Flux-balance analysis (FBA) is a computational method based on linear programming (LP) that has had enormous success modeling the metabolic behaviour of organisms and cellular systems existing in steady state with their environment [@Varma94; @Orth10]. Representing the corresponding biological model as an LP problem means that FBA can be used to study metabolism at genome-scale. Unfortunately, the underlying assumption of an unchanging environment means that FBA is not immediately applicable to systems with dynamics where, for example, environmental conditions may vary in time. Extensions of FBA, including dynamic FBA (DFBA) [@Mahadevan02], have therefore been developed in order to accommodate temporal dynamics into the framework of genome-scale metabolic modeling.

Although DFBA is well-defined mathematically as an LP problem embedded in a system of ordinary differential equations (ODEs), numerical simulation of DFBA models proves particularly challenging (as described in @Harwood16). Consequently, @Harwood16 proposed an algorithm for efficiently simulating DFBA models based on reformulating the ODE system as a differential algebraic equation (DAE) with root detection and representing solutions of the LP problem using an optimal basis formulation. An initial implementation of this algorithm has been provided in the software package DFBAlab [@Gomez14] written in MATLAB and using commercial LP solvers.

Increasingly, researchers engaged in metabolic modeling prefer open source software. Python is quickly becoming their platform of choice [@Carey20]. Among other packages, open source resources for building and simulating FBA models using Python can be found in the COBRApy [@Ebrahim13] package which is part of the openCOBRA organization [@opencobra]. Until now, COBRApy lacked an efficient implementation of DFBA using the DAE formulation.

## Statement of need:
_Researchers wanting to build and simulate specific models of interest often lack the background in numerical analysis or high-performance computing required to overcome the numerical challenges of DFBA_.

We have solved this issue by developing a software package based on open source libraries GLPK [@glpk] and SUNDIALS [@Hindmarsh05] that implements the most efficient algorithms in a compiled programming language that is made accessible to users through a simple and intuitive pybind11 [@pybind11] interface with pandas [@McKinney11] and the openCOBRA Python module COBRApy [@Ebrahim13]. ODEs are constructed using the symbolic expression enabled by optlang [@Jensen16], SymPy, and SymEngine [@Meurer17].

# Acknowledgements

DST is a Simons Foundation Fellow of the Life Sciences Research Foundation. MEB
received funding from the European Union’s Horizon 2020 research and innovation
programme under grant agreement 686070 (DD-DeCaF). We thank Peter St. John and Christian Diener for discussions and suggestions.

# References
## *EMBLP* source files

This directory contains the following content

* [`./emblp_direct.cpp`](./emblp_direct.cpp): contains member functions for
  direct method derived class
* [`./emblp_harwood.cpp`](./emblp_harwood.cpp): contains member functions for
  Harwood et al. derived class
* [`./emblp.cpp`](./emblp.cpp): contains member functions for embedded
  LP problem base class
* [`./emblp_scott.h`](./emblp_scott.h): contains member functions for Scott et
  al. derived class
* [`./emblp.h`](./emblp.h): contains class declarations for embedded LP problems
* [`./README.md`](./README.md): this document

## *METHODS* source files

This directory contains the following content

* [`./methods_direct.cpp`](./methods_direct.cpp): contains source code for model integration
  using direct method
* [`./methods_harwood.cpp`](./methods_harwood.cpp): contains source code for model integration
  using Harwood et al. algorithm
* [`./methods.h`](./methods.h): contains *SUNDIALS* includes and declarations
  for integration methods
* [`./methods.cpp`](./methods.cpp): contains functions used by methods
* [`./methods_scott.cpp`](./methods_scott.cpp): contains source code for model integration using
  Scott et al. algorithm
* [`./README.md`](./README.md): this document

#### Problem description

Please explain:
* **what** you tried to achieve,
* **how** you went about it (referring to the code sample), and
* **why** the current behaviour is a problem and what output
  you expected instead.

#### Code Sample

Create a [minimal, complete, verifiable example
](https://stackoverflow.com/help/mcve).

```python
# Paste your code here or link to a gist.
```

```
# If there was a crash, please include the traceback here.
```

### Context

Please run the following code and paste the output inside the details
block.

```
python -c "import dfba;dfba.show_versions()"
```

<details>

</details>
* [ ] fix #(issue number)
* [ ] description of feature/fix
* [ ] tests added/passed
* [ ] add an entry to the [next release](../../CHANGELOG.rst)
## Examples

This directory contains the following content

* [`./example1.py`](./example1.py): anaerobic growth of *E. coli* on glucose and
  xylose
* [`./example2.py`](./example2.py): aerobic growth of *E. coli* on glucose and
  xylose
* [`./example3.py`](./example3.py): aerobic growth of *S. cerevisiae* on glucose
  with switch to anaerobic conditions at *t=7.7h*
* [`./example4.py`](./example4.py): aerobic growth of *S. cerevisiae* on glucose and
  xylose
* [`./example5.py`](./example5.py): anaerobic growth of *E. coli* on glucose and
  xylose simulated using direct method
* [`./README.md`](./README.md): this document

## Build scripts

This directory contains the following content

* [`./build_pybind11.sh`](./build_pybind11.sh): script for downloading, building and
  installing specified pybind11
* [`./build_glpk.sh`](./build_glpk.sh): script for downloading, building and
  installing specified GLPK version
* [`./build_sundials.sh`](./build_sundials.sh): script for downloading, building
  and installing specified SUNDIALS version
* [`./jupyterlab_plotly.sh`](./jupyterlab_plotly.sh): script for building jupyterlab extensions
* [`./README.md`](./README.md): this document

## cmake modules

* [`./FindGLPK.cmake`](./FindGLPK.cmake) adapted from [ycollet/coinor-cmake](https://github.com/ycollet/coinor-cmake/blob/master/Clp/cmake/FindGlpk.cmake)
* [`./FindSUNDIALS.cmake`](./FindSUNDIALS.cmake) adapted from [casadi/casadi](https://github.com/casadi/casadi/blob/master/cmake/FindSUNDIALS.cmake)
* [`./README.md`](./README.md): this document

