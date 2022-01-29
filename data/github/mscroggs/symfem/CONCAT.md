# Version 2021.12.1 (01 December 2021)

- Added Regge element on tensor product cells.
- Added tensor product factorisation of Q element.
- Added Arnold-Boffi-Falk element
- Added Arbogast-Correa element

# Version 2021.10.1 (19 October 2021)

- Removed discontinuous elements.
- Added pictures of reference cells to readme.
- Added trimmed serendipity Hdiv and Hcurl elements.
- Added TNT scalar, Hdiv and Hcurl elements.

# Version 2021.8.2 (24 August 2021)

- New release for JOSS paper

# Version 2021.8.1 (11 August 2021)

- Added Guzman-Neilan element
- Added nonconforming Arnold-Winther element
- Wrote JOSS paper

# 11 August 2021

- Started keeping changelog
![Symfem](logo/logo.png)

|  | Badges |
| --- | :---: |
| Documentation | [![Documentation status](https://readthedocs.org/projects/symfem/badge/?version=latest)](https://symfem.readthedocs.io/en/latest/?badge=latest) |
| Testing & coverage | [![Style checks](https://github.com/mscroggs/symfem/actions/workflows/style-checks.yml/badge.svg)](https://github.com/mscroggs/symfem/actions) [![Run tests](https://github.com/mscroggs/symfem/actions/workflows/run-tests.yml/badge.svg)](https://github.com/mscroggs/symfem/actions) [![Coverage Status](https://coveralls.io/repos/github/mscroggs/symfem/badge.svg?branch=main)](https://coveralls.io/github/mscroggs/symfem?branch=main) |
| Packages | [![PyPI](https://img.shields.io/pypi/v/symfem?color=blue&label=PyPI&logo=pypi&logoColor=white)](https://pypi.org/project/symfem/) [![conda](https://anaconda.org/conda-forge/symfem/badges/version.svg)](https://anaconda.org/conda-forge/symfem) |
| Paper | [![DOI](https://joss.theoj.org/papers/10.21105/joss.03556/status.svg)](https://doi.org/10.21105/joss.03556) |
| Code quality | [![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=mscroggs_symfem&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=mscroggs_symfem) |

Symfem is a symbolic finite element definition library, that can be used to
symbolically evaluate the basis functions of a finite element space. Symfem can:

- Symbolically compute the basis functions of a wide range of finite element spaces
- Symbolically compute derivatives and vector products and substitute values into functions
- Allow the user to define their own element using the Ciarlet definition of a finite element
- Be used to verify that the basis functions of a given space have some desired properties

You can find details of recent changes to Symfem in [the changelog](CHANGELOG.md).

# Installing Symfem
## Installing from source
Symfem can be installed by downloading the [GitHub repo](https://github.com/mscroggs/symfem)
and running:

```bash
python3 setup.py install
```

## Installing using pip
The latest release of Symfem can be installed by running:

```bash
pip3 install symfem
```

## Installing using conda
The latest release of Symfem can be installed by running:

```bash
conda install -c conda-forge symfem
```

# Testing Symfem
To run the Symfem unit tests, clone the repository and run:

```bash
python3 -m pytest test/
```

You may instead like to run the following, as this will skip the slowest tests.

```bash
python3 -m pytest test/ --speed fast
```

# Using Symfem
Finite elements can be created in Symfem using the `symfem.create_element()`
function. For example, some elements are created in the following snippet:

```python
import symfem

lagrange = symfem.create_element("triangle", "Lagrange", 1)
rt = symfem.create_element("tetrahedron", "Raviart-Thomas", 2)
nedelec = symfem.create_element("triangle", "N2curl", 1)
qcurl = symfem.create_element("quadrilateral", "Qcurl", 2)
```

The polynomial basis of an element can be obtained by calling `get_polynomial_basis()`:

```python
import symfem

lagrange = symfem.create_element("triangle", "Lagrange", 1)
print(lagrange.get_basis_functions())
```
```
[-x - y + 1, x, y]
```

Each basis function will be a [Sympy](https://www.sympy.org) symbolic expression.

Derivative of these basis functions can be computed using the functions in
[`symfem.calculus`](symfem/calculus.py). Vector-valued basis functions can
be manipulated using the functions in [`symfem.vector`](symfem/vectors.py).

The function `map_to_cell` can be used to map the basis functions of a finite element
to a non-default cell:

```python
import symfem

lagrange = symfem.create_element("triangle", "Lagrange", 1)
print(lagrange.get_basis_functions())
print(lagrange.map_to_cell([(0,0), (2, 0), (2, 1)]))
```
```
[-x - y + 1, x, y]
[1 - x/2, x/2 - y, y]
```

## Further documentation
More detailed documentation of the latest release version of Symfem can be found on
[Read the Docs](https://symfem.readthedocs.io/en/latest/). A series of example uses
of Symfem can be found in the [`demo` folder](demo/) or viewed on
[Read the Docs](https://symfem.readthedocs.io/en/latest/demos/index.html).

Details of the definition of each element can be found on [DefElement](https://defelement.com)
alongside Symfem snippets for creating the element.

## Getting help
You can ask questions about using Symfem by using [GitHub Discussions](https://github.com/mscroggs/symfem/discussions).
Bugs can be reported using the [GitHub issue tracker](https://github.com/mscroggs/symfem/issues).

# Contributing to Symfem
## Reporting bugs
If you find a bug in Symfem, please report it on the [issue tracker](https://github.com/mscroggs/symfem/issues/new?assignees=&labels=bug&template=bug_report.md&title=).

## Suggesting enhancements
If you want to suggest a new feature or an improvement of a current feature, you can submit this
on the [issue tracker](https://github.com/mscroggs/symfem/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=).

## Submitting a pull request
If you want to directly submit code to Symfem, you can do this by forking the Symfem repo, then submitting a pull request.
If you want to contribute, but are unsure where to start, have a look at the
[issues labelled "good first issue"](https://github.com/mscroggs/symfem/issues?q=is%3Aopen+is%3Aissue+label%3A%22good+first+issue%22).

On opening a pull request, unit tests and flake8 style checks will run. You can click on these in the pull request
to see where (if anywhere) there are errors in your code.

## Code of conduct
We expect all our contributors to follow [the Contributor Covenant](CODE_OF_CONDUCT.md). Any unacceptable
behaviour can be reported to Matthew (symfem@mscroggs.co.uk).

# Available cells and elements
## Interval
The reference interval has vertices (0,) and (1,). Its sub-entities are numbered as follows.

![The numbering of a reference interval](img/interval_numbering.png)

### List of supported elements
- Bernstein (alternative names: Bernstein-Bezier)
- bubble
- dPc
- Hermite
- Lagrange (alternative names: P)
- Morley-Wang-Xu (alternative names: MWX)
- serendipity (alternative names: S)
- Taylor (alternative names: discontinuous Taylor)
- vector Lagrange (alternative names: vP)
- Wu-Xu

## Triangle
The reference triangle has vertices (0, 0), (1, 0), and (0, 1). Its sub-entities are numbered as follows.

![The numbering of a reference triangle](img/triangle_numbering.png)

### List of supported elements
- Argyris
- Arnold-Winther (alternative names: AW, conforming Arnold-Winther)
- Bell
- Bernardi-Raugel
- Bernstein (alternative names: Bernstein-Bezier)
- Brezzi-Douglas-Fortin-Marini (alternative names: BDFM)
- Brezzi-Douglas-Marini (alternative names: BDM, N2div)
- bubble
- bubble enriched Lagrange
- bubble enriched vector Lagrange
- conforming Crouzeix-Raviart (alternative names: conforming CR)
- Crouzeix-Raviart (alternative names: CR, Crouzeix-Falk, CF)
- Fortin-Soulie (alternative names: FS)
- Guzman-Neilan
- Hellan-Herrmann-Johnson (alternative names: HHJ)
- Hermite
- Hsieh-Clough-Tocher (alternative names: Clough-Tocher, HCT, CT)
- Kong-Mulder-Veldhuizen (alternative names: KMV)
- Lagrange (alternative names: P)
- Mardal-Tai-Winther (alternative names: MTW)
- matrix Lagrange
- Morley
- Morley-Wang-Xu (alternative names: MWX)
- Nedelec (alternative names: Nedelec1, N1curl)
- Nedelec2 (alternative names: N2curl)
- nonconforming Arnold-Winther (alternative names: nonconforming AW)
- Raviart-Thomas (alternative names: RT, N1div)
- reduced Hsieh-Clough-Tocher (alternative names: rHCT)
- Regge
- symmetric matrix Lagrange
- Taylor (alternative names: discontinuous Taylor)
- transition
- vector Lagrange (alternative names: vP)
- Wu-Xu

## Quadrilateral
The reference quadrilateral has vertices (0, 0), (1, 0), (0, 1), and (1, 1). Its sub-entities are numbered as follows.

![The numbering of a reference quadrilateral](img/quadrilateral_numbering.png)

### List of supported elements
- Arbogast-Correa (alternative names: AC, AC full, Arbogast-Correa full)
- Arnold-Boffi-Falk (alternative names: ABF)
- Bogner-Fox-Schmit (alternative names: BFS)
- Brezzi-Douglas-Fortin-Marini (alternative names: BDFM)
- bubble
- direct serendipity
- dPc
- NCE (alternative names: RTCE, Qcurl, Nedelec, Ncurl)
- NCF (alternative names: RTCF, Qdiv)
- Q (alternative names: Lagrange, P)
- Regge
- serendipity (alternative names: S)
- serendipity Hcurl (alternative names: Scurl, BDMCE, AAE)
- serendipity Hdiv (alternative names: Sdiv, BDMCF, AAF)
- tiniest tensor (alternative names: TNT)
- tiniest tensor Hcurl (alternative names: TNTcurl)
- tiniest tensor Hdiv (alternative names: TNTdiv)
- trimmed serendipity Hcurl (alternative names: TScurl)
- trimmed serendipity Hdiv (alternative names: TSdiv)
- vector dPc
- vector Q (alternative names: vQ)

## Tetrahedron
The reference tetrahedron has vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), and (0, 0, 1). Its sub-entities are numbered as follows.

![The numbering of a reference tetrahedron](img/tetrahedron_numbering.png)

### List of supported elements
- Bernardi-Raugel
- Bernstein (alternative names: Bernstein-Bezier)
- Brezzi-Douglas-Fortin-Marini (alternative names: BDFM)
- Brezzi-Douglas-Marini (alternative names: BDM, N2div)
- bubble
- Crouzeix-Raviart (alternative names: CR, Crouzeix-Falk, CF)
- Guzman-Neilan
- Hermite
- Kong-Mulder-Veldhuizen (alternative names: KMV)
- Lagrange (alternative names: P)
- Mardal-Tai-Winther (alternative names: MTW)
- matrix Lagrange
- Morley-Wang-Xu (alternative names: MWX)
- Nedelec (alternative names: Nedelec1, N1curl)
- Nedelec2 (alternative names: N2curl)
- Raviart-Thomas (alternative names: RT, N1div)
- Regge
- symmetric matrix Lagrange
- Taylor (alternative names: discontinuous Taylor)
- transition
- vector Lagrange (alternative names: vP)
- Wu-Xu

## Hexahedron
The reference hexahedron has vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), and (1, 1, 1). Its sub-entities are numbered as follows.

![The numbering of a reference hexahedron](img/hexahedron_numbering.png)

### List of supported elements
- Brezzi-Douglas-Duran-Fortin (alternative names: BDDF)
- Brezzi-Douglas-Fortin-Marini (alternative names: BDFM)
- bubble
- dPc
- NCE (alternative names: RTCE, Qcurl, Nedelec, Ncurl)
- NCF (alternative names: RTCF, Qdiv)
- Q (alternative names: Lagrange, P)
- Regge
- serendipity (alternative names: S)
- serendipity Hcurl (alternative names: Scurl, BDMCE, AAE)
- serendipity Hdiv (alternative names: Sdiv, BDMCF, AAF)
- tiniest tensor (alternative names: TNT)
- tiniest tensor Hcurl (alternative names: TNTcurl)
- tiniest tensor Hdiv (alternative names: TNTdiv)
- trimmed serendipity Hcurl (alternative names: TScurl)
- trimmed serendipity Hdiv (alternative names: TSdiv)
- vector dPc
- vector Q (alternative names: vQ)

## Prism
The reference prism has vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), and (0, 1, 1). Its sub-entities are numbered as follows.

![The numbering of a reference prism](img/prism_numbering.png)

### List of supported elements
- Lagrange (alternative names: P)
- Nedelec (alternative names: Ncurl)

## Pyramid
The reference pyramid has vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), and (0, 0, 1). Its sub-entities are numbered as follows.

![The numbering of a reference pyramid](img/pyramid_numbering.png)

### List of supported elements
- Lagrange (alternative names: P)

## Dual polygon
The reference dual polygon (hexagon example shown) has vertices (1, 0), (3/4, sqrt(3)/4), (1/2, sqrt(3)/2), (0, sqrt(3)/2), (-1/2, sqrt(3)/2), (-3/4, sqrt(3)/4), (-1, 0), (-3/4, -sqrt(3)/4), (-1/2, -sqrt(3)/2), (0, -sqrt(3)/2), (1/2, -sqrt(3)/2), and (3/4, -sqrt(3)/4). Its sub-entities are numbered as follows.

![The numbering of a reference dual polygon](img/dual_polygon_numbering.png)

### List of supported elements
- Buffa-Christiansen (alternative names: BC)
- dual polynomial (alternative names: dual P, dual)
- rotated Buffa-Christiansen (alternative names: RBC)


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
symfem@mscroggs.co.uk (Matthew Scroggs). All complaints will be reviewed and
investigated promptly and fairly.

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

## How to contribute

### Reporting bugs
If you find a bug in Symfem, please report it on the [issue tracker](https://github.com/mscroggs/symfem/issues/new?assignees=&labels=bug&template=bug_report.md&title=).

### Suggesting enhancements
If you want to suggest a new feature or an improvement of a current feature, you can submit this
on the [issue tracker](https://github.com/mscroggs/symfem/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=).

### Submitting a pull request
If you want to directly submit code to Symfem, you can do this by forking the Symfem repo, then submitting a pull request.
If you want to contribute, but are unsure where to start, have a look at the
[issues labelled "good first issue"](https://github.com/mscroggs/symfem/issues?q=is%3Aopen+is%3Aissue+label%3A%22good+first+issue%22).

On opening a pull request, unit tests and flake8 style checks will run. You can click on these in the pull request
to see where (if anywhere) there are errors in your code.

### Code of conduct
We expect all our contributors to follow [the Contributor Covenant](CODE_OF_CONDUCT.md). Any unacceptable
behaviour can be reported to Matthew (symfem@mscroggs.co.uk).
Changes in this pull request:

- List
- your
- changes
- here

Fixes # (issue)
---
name: New element
about: Suggest a new element
title: 'Add [NAME] element'
labels: enhancement
assignees: ''

---

**Describe the element you want to add**
Briefly describe the element.

**References**
Link to reference(s) where the element is defined.

**Additional context**
Add any other context about the element here.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Minimal code to reproduce the behavior:
```python
import symfem
assert 1 == 2
```
**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Additional context**
Add any other context or screenshots about the feature request here.
---
title: 'Symfem: a symbolic finite element definition library'
tags:
  - Python
  - finite element method
  - basis functions
  - symbolic algebra
  - numerical analysis
authors:
  - name: Matthew W. Scroggs
    orcid: 0000-0002-4658-2443
    affiliation: 1
affiliations:
 - name: Department of Engineering, University of Cambridge
   index: 1
date: 15 July 2021
bibliography: paper.bib
---

# Summary

The finite element method (FEM) [@ciarlet] is a popular method for numerically solving a wide
range of partial differential equations (PDEs). To solve a problem using FEM, the PDE is first
written in a weak form, for example: find $u\in V$ such that for all $v\in V,$

\begin{equation}
\int_\Omega \nabla u\cdot\nabla v=\int_\Omega fv,
\end{equation}

where $f$ is a known function, and $\Omega$ is the domain on which the problem is being solved.
This form is then discretised by defining a finite-dimensional subspace of $V$---often called
$V_h$---and looking for a solution $u_h\in V_h$ that satisfies the above equation for all functions
$v_h\in V_h$. These finite-dimensional subspaces are defined by meshing the domain of the problem,
then defining a set of basis functions on each cell in the mesh (and enforcing any desired
continuity between the cells).

For different applications, there are a wide range of finite-dimensional spaces that can be used.
Symfem is a Python library that can be used to symbolically compute basis functions of these
spaces. The symbolic representations are created using Sympy [@sympy], allowing
them to be easily manipulated using Sympy's functionality once they are created.

# Statement of need

In FEM libraries, it is common to define basis functions so that they, and their
derivatives, can quickly and efficiently be evaluated at a collection of points, thereby allowing
full computations to be completed quickly. The libraries FIAT [@fiat] and Basix [@basix]---which
are part of the FEniCS project [@fenics]---implement this functionality as stand-alone libraries.
Many other FEM libraries define their basis functions as part of the core library functionality.
It is not common to be able to compute a symbolic representation of the basis functions.

Symfem offers a wider range of finite element spaces than other FEM libraries, and the ability
to symbolically compute basis functions. There are a number of situations in which the symbolic
representation of a basis function is useful: it is easy to confirm, for example, that the
derivatives of the basis functions have a certain desired property, or check what they are
equal to when restricted to one face or edge of the cell.

Symfem can also be used to explore the behaviour of the wide range of spaces it supports, so the
user can decide which spaces to implement in a faster way in their FEM code. Additionally,
Symfem can be used to prototype new finite element spaces, as custom spaces can easily be
added, then it can be checked that the basis functions of the space behave as expected.

As basis functions are computed symbolically in Symfem, it is much slower than the alternative
libraries. It is therefore not suitable for performing actual finite element calculations. It
should instead be seen as a library for research and experimentation.

# References
