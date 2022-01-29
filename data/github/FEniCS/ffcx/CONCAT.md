# FFCx: The FEniCSx Form Compiler

[![FFCx CI](https://github.com/FEniCS/ffcx/actions/workflows/pythonapp.yml/badge.svg)](https://github.com/FEniCS/ffcx/actions/workflows/pythonapp.yml)
[![Coverage Status](https://coveralls.io/repos/github/FEniCS/ffcx/badge.svg?branch=main)](https://coveralls.io/github/FEniCS/ffcx?branch=main)

FFCx is a new version of the FEniCS Form Compiler. It is being actively
developed and is compatible with DOLFINx.

FFCx is a compiler for finite element variational forms. From a
high-level description of the form in the Unified Form Language (UFL),
it generates efficient low-level C code that can be used to assemble the
corresponding discrete operator (tensor). In particular, a bilinear form
may be assembled into a matrix and a linear form may be assembled into a
vector.  FFCx may be used either from the command line (by invoking the
`ffcx` command) or as a Python module (`import ffcx`).

FFCx is part of the FEniCS Project. For more information, visit
https://www.fenicsproject.org


## Installation

To install FFCx from PyPI:
```
$ pip install ffcx
```

To install FFCx from the source directory:
```
$ pip install .
```

## Documentation

Documentation can be viewed at https://docs.fenicsproject.org/ffcx/main


## Interface file installation only

FFCx provides the `ufcx.h` interface header for finite element kernels,
used by DOLFINx. `ufcx.h` is installed by FFCx within the Python site
packages, but it is sometimes helpful to install only the header file.
This can be done using `cmake`:
```
$ cmake -B build-dir -S cmake/
$ cmake --build build-dir
$ cmake --install build-dir
```

## License

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program. If not, see <https://www.gnu.org/licenses/>.
Code of Conduct
===============

Our Pledge
----------
In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our
project and our community a harassment-free experience for everyone,
regardless of age, body size, disability, ethnicity, sex
characteristics, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance,
race, religion, or sexual identity and orientation.

Our Standards
-------------
Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others’ private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

Our Responsibilities
--------------------
Project maintainers are responsible for clarifying the standards of
acceptable behavior and are expected to take appropriate and fair
corrective action in response to any instances of unacceptable
behavior.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that
they deem inappropriate, threatening, offensive, or harmful.

Scope
-----
This Code of Conduct applies both within project spaces and in public
spaces when an individual is representing the project or its
community. Examples of representing a project or community include
using an official project e-mail address, posting via an official
social media account, or acting as an appointed representative at an
online or offline event. Representation of a project may be further
defined and clarified by project maintainers.

Enforcement
-----------
Instances of abusive, harassing, or otherwise unacceptable behavior
may be reported by contacting the project team at
fenics-steering-council@googlegroups.com. Alternatively, you may
report individually to one of the members of the Steering
Council. Complaints will be reviewed and investigated and will result
in a response that is deemed necessary and appropriate to the
circumstances. The project team is obligated to maintain
confidentiality with regard to the reporter of an incident. Further
details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct
in good faith may face temporary or permanent repercussions as
determined by other members of the project’s leadership.

If you feel that your report has not been followed up satisfactorily,
then you may contact our parent organisation NumFOCUS at
info@numfocus.org for further redress.

Attribution
-----------
This Code of Conduct is adapted from the Contributor Covenant, version
1.4, available at
https://www.contributor-covenant.org/version/1/4/code-of-conduct.html.

Adaptations
-----------

* Allow reporting to individual Steering Council members
* Added the option to contact NumFOCUS for further redress.

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq