---
title: "The Basic Model Interface 2.0: A standard interface for coupling numerical models in the geosciences"
tags:
  - C
  - C++
  - Fortran
  - Python
  - geosciences
  - modeling
  - interface
authors:
  - name: Eric W.H. Hutton
    orcid: 0000-0002-5864-6459
    affiliation: 1
  - name: Mark D. Piper
    orcid: 0000-0001-6418-277X
    affiliation: 1
  - name: Gregory E. Tucker
    orcid: 0000-0003-0364-5800
    affiliation: 1, 2, 3
affiliations:
  - name: Community Surface Dynamics Modeling System, University of Colorado Boulder, USA
    index: 1
  - name: Cooperative Institute for Research in Environmental Sciences (CIRES), University of Colorado Boulder, USA
    index: 2
  - name: Department of Geological Sciences, University of Colorado Boulder, USA
    index: 3
date: 10 January 2020
bibliography: paper.bib

---

# Summary

Component modeling is a research technique
in which new models are constructed by coupling the inputs and outputs
of simpler existing models. 
Component modeling traces its roots
to component-based software engineering,
where a software system is constructed from a number
of independent, reusable software components,
each encapsulating a unit of functionality
and exposing inputs and outputs through an interface.
A tangible analogy is a bicycle.
A bicycle is a system of reusable, replaceable components.
Tires are one of the components.
You can easily swap in a studded tire for icy winter streets,
then swap it out again in the summer.

While there is a longer history of component modeling
in fields such as climate modeling,
with, for example, the Earth System Modeling Framework [@collins:2005],
component modeling is relatively new
to the earth surface processes community.
Some recent examples include
@ratliff:2018, who show that a river model transporting sediment
can feed a delta model that distributes the sediment,
and @hoch:2019, who show that coupling hydrologic and hydrodynamic models
may sharpen inundation estimates in flood modeling.

In component-based software engineering,
components communicate through interfaces:
named sets of functions with prescribed arguments and return values.
The bicycle analogy above benefits from a standard interface
for tire diameter and width.
Likewise,
component modeling can benefit from an interface
for describing the inputs, outputs, and behaviors of a model.
The *Basic Model Interface* (BMI)
provides a standard set of functions
for querying, modifying, and running a model.
Equipping a model with a BMI
allows the model to be coupled with other models that expose a BMI.
The BMI concept was introduced by @peckham:2013
as a foundational technology for the proposed
[Community Surface Dynamics Modeling System](https://csdms.colorado.edu)
(CSDMS)
model coupling framework.
The current work represents an update to the original BMI,
with new functions for describing variables
and for working with structured and unstructured grids.
Full documentation for the current version of the BMI
is available at https://bmi.readthedocs.io.

The functions that comprise the BMI are designed
to be straightforward to implement in any programming language,
using only simple data types from standard language libraries.
To generalize across languages,
the BMI is expressed in the Scientific Interface Definition Language
[@epperly:2012].
BMI specifications for four languages---C, C++, Fortran, and Python---have
been derived from the SIDL specification.
For each language,
links to the GitHub repositories supplying the specification
and an example implementation are listed in Table 1 below.
Detailed instructions for building the language specification and example
are given within each repository.

\  

**Table 1:**
Repositories containing BMI language specifications and examples.
Prefix the CSDMS GitHub organization (https://github.com/csdms/) to the
repository name to obtain the full repository URL.

| Language | Specification | Example implementation |
| -------- | ------------- | ---------------------- |
| C        | [bmi-c]       | [bmi-example-c]        |
| C++      | [bmi-cxx]     | [bmi-example-cxx]      |
| Fortran  | [bmi-fortran] | [bmi-example-fortran]  |
| Python   | [bmi-python]  | [bmi-example-python]   |

[bmi-c]: https://github.com/csdms/bmi-c
[bmi-cxx]: https://github.com/csdms/bmi-cxx
[bmi-fortran]: https://github.com/csdms/bmi-fortran
[bmi-python]: https://github.com/csdms/bmi-python
[bmi-example-c]: https://github.com/csdms/bmi-example-c
[bmi-example-cxx]: https://github.com/csdms/bmi-example-cxx
[bmi-example-fortran]: https://github.com/csdms/bmi-example-fortran
[bmi-example-python]: https://github.com/csdms/bmi-example-python

While CSDMS currently supports the four languages listed in Table 1,
a BMI can be created for any language.
BMI is a community-driven standard;
contributions that follow the contributor code of conduct
listed in the main BMI repository are welcomed,
and are acknowledged in the repository and the documentation.

# Acknowledgements

This work is supported by the National Science Foundation
under Grant No. 1831623, *Community Facility Support: The
Community Surface Dynamics Modeling System (CSDMS)*.

# References
Citation
========

If you use the Basic Model Interface for work/research
presented in a publication, we ask that you please cite:

Peckham, S.D., Hutton, E.W., and Norris, B., 2013. A component-based approach to integrated modeling in the geosciences: The design of CSDMS. *Computers & Geosciences*, 53, pp.3-12, http://dx.doi.org/10.1016/j.cageo.2012.04.002.

Hutton, E.W.H., Piper, M.D., and Tucker, G.E., 2020. The Basic Model Interface 2.0: A standard interface for coupling numerical models in the geosciences. *Journal of Open Source Software*, 5(51), 2317, https://doi.org/10.21105/joss.02317

BibTeX entries:

::

  @article{peckham2013component,
   title={A component-based approach to integrated modeling in the geosciences: The design of CSDMS},
   author={Peckham, Scott D and Hutton, Eric WH and Norris, Boyana},
   journal={Computers \& Geosciences},
   volume={53},
   pages={3--12},
   year={2013},
   publisher={Elsevier}
  }

  @article{hutton2020basic,
   doi = {10.21105/joss.02317},
   url = {https://doi.org/10.21105/joss.02317},
   year = {2020},
   publisher = {The Open Journal},
   volume = {5},
   number = {51},
   pages = {2317},
   author = {Eric W.H. Hutton and Mark D. Piper and Gregory E. Tucker},
   title = {The Basic Model Interface 2.0: A standard interface for coupling numerical models in the geosciences},
   journal = {Journal of Open Source Software}
  }
.. highlight:: shell

============
Contributing
============

Contributions are welcomed and greatly appreciated!
BMI is a community-driven open source project, so every little bit helps.
`Credit`_ will always be given.

This document describes how to contribute to the BMI project.
It covers the main `BMI repository`_,
as well as the language specifications and examples
for supported languages
listed in `Table 1`_ of the `BMI documentation`_.

--------------
Making changes
--------------

The following sections describe the process through which changes to the BMI
specification, language bindings, language examples, and support tools are made.

Pull requests
~~~~~~~~~~~~~

If you'd like to propose a change to the BMI or any of its language
bindings, examples, or tools, you must submit a pull request (PR) to the
corresponding GitHub repository. If you have questions regarding a potential
change, or if you're unsure whether your proposed changes or additions are
acceptable, we recommend you first submit your question as a GitHub issue in the
corresponding repository, which can later be turned into a PR.

When a PR is opened against the BMI or any of its language bindings, examples,
or tools, it will be tagged either as a *Request for Comment* (RFC) or for
resolution through *lazy consensus*. The RFC process is intended for major
changes to the BMI, whereas lazy consensus applies to minor changes such as
typos or documentation updates.

Request for Comment (RFC)
~~~~~~~~~~~~~~~~~~~~~~~~~

Once a PR has been tagged as an RFC, a review process begins.
It proceeds in three steps:

1. discussion
2. summary
3. resolution

Discussion
..........

Reviewers (the other developers and interested community members) will write
comments on your PR to help you improve its implementation, documentation, and
style. Every single developer working on the project has their code reviewed,
and we've come to see it as friendly conversation from which we all learn and
the overall code quality benefits. Therefore, please don't let the review
discourage you from contributing: its only aim is to improve the quality of the
project, not to criticize (we are, after all, very grateful for the time you're
donating!).

The author of the PR should try to build consensus and integrate feedback. An
RFC with broad support will be more likely to be quickly accepted than one
without comments. If you feel that your PR has been forgotten,
please feel free to reach out directly to the project team.

We expect that most PRs will not be accepted without some modification. You can
make changes by pushing commits directly to the PR branch and leaving comments
that explain the modifications.

Summary
.......

Once all parties have weighed in and there has been adequate discussion, the PR
moves to a summary period. A team member will move for the discussion
phase to end and for the summary period to begin, which will end with a
recommendation for the RFC (merge, close, or, postpone). Team members making
this motion should use their best judgment that adequate discussion has taken
place and that the consequences of the proposed change have been adequately
addressed.

The summary period will typically last about a week, allowing
stakeholders to have a chance to lodge any final objections before the team
reaches a final decision.

The intent is for this summary period to be fairly quiet, with most of the
discussion taking place during the discussion phase. There will likely be times,
however, where stakeholders might raise serious enough concerns that the PR
moves back into the discussion phase.

Resolution
..........

The RFC process concludes with the PR being merged, closed, or postponed.

There must be a *consensus* from project developers and interested community
members for a PR to be merged. Consensus has a particular meaning when used
with open source projects; see the `BMI governance document`_ for a definition
of consensus in this context.

If consensus isn't achieved, the PR will be postponed (the team feels the PR can
wait until some future time) or closed.

Lazy Consensus
~~~~~~~~~~~~~~

`Lazy consensus`_, as defined by the Apache Software Foundation, is a
decision-making policy which assumes general consent if no responses are posted
within a defined period.

For the BMI project, a PR tagged for resolution through lazy consensus can be
merged if no comments are posted within about a week.

-------------------
Reporting a problem
-------------------

Report bugs or other problems by creating a GitHub issue at
https://github.com/csdms/bmi/issues.

In the issue, be sure to explain the problem and include additional details to
help maintainers reproduce the problem. Here are some suggestions that will make
it easier to track down the source of the problem:

* Use a clear and descriptive title for the issue that identifies the problem.
* Describe, and if possible provide a minimal example that demonstrates, the
  exact steps that reproduce the problem.
* Describe the behavior you are seeing after these steps.
* Describe the behavior you expect to see after these steps.

------------------
Providing feedback
------------------

The best way to send feedback is to file an issue at
https://github.com/csdms/bmi/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome.


.. Links

.. _Credit: https://bmi.readthedocs.io/en/latest/credits.html
.. _BMI repository: https://github.com/csdms/bmi
.. _Table 1: https://bmi.readthedocs.io/en/latest/#id43
.. _BMI documentation: https://bmi.readthedocs.io
.. _BMI governance document: https://bmi.readthedocs.io/en/latest/governance.html
.. _Lazy consensus: https://community.apache.org/committers/lazyConsensus.html
=======
Credits
=======

Project Coordinators
--------------------

* Eric Hutton <mcflugen@gmail.com>
* Mark Piper <mark.piper@colorado.edu>

Contributors
------------

* Richard Barnes
* Michael Galloy
* Julian Hofer
* Eric Hutton
* Allen Lee
* Eric Morway
* Boyana Norris
* Scott Peckham
* Mark Piper
* Mike Taves
* Greg Tucker
* Ben van Werkhoven
* Martijn Visser

If you have contributed to the BMI project and your name is missing,
please send an email to the coordinators, or open a pull request
for this page in the `BMI repository <https://github.com/csdms/bmi>`_.
Contributor Covenant Code of Conduct
====================================

Our Pledge
----------

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

Our Standards
-------------

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

Our Responsibilities
--------------------

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
-----

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at csdms@colorado.edu. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
-----------

This Code of Conduct is adapted from the
`Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at
`<https://www.contributor-covenant.org/version/1/4/code-of-conduct.html>`_.

For answers to common questions about this code of conduct, see
`<https://www.contributor-covenant.org/faq>`_.
.. role:: raw-html-m2r(raw)
   :format: html


.. raw:: html

   <p align="center">
      <a href='https://bmi.readthedocs.org/'>
         <img src='https://github.com/csdms/bmi/raw/master/docs/source/_static/bmi-logo-header-text.png'/>
      </a>
   </p>

.. raw:: html

   <h2 align="center">The Basic Model Interface</h2>



.. raw:: html

   <p align="center">

   <a href='https://doi.org/10.5281/zenodo.3955009'>
     <img src='https://zenodo.org/badge/DOI/10.5281/zenodo.3955010.svg' alt='DOI'></a>
   <a href="https://doi.org/10.21105/joss.02317">
     <img src="https://joss.theoj.org/papers/10.21105/joss.02317/status.svg" alt="JOSS article"></a>
   <a href='https://bmi.readthedocs.io/en/latest/?badge=latest'>
     <img src='https://readthedocs.org/projects/bmi/badge/?version=latest' alt='Documentation Status'></a>
   <a href="https://opensource.org/licenses/MIT">
     <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg"></a>

   </p>

   <p align="center">

   The Basic Model Interface (BMI) is a standardized set of functions
   that allows coupling of models to models and models to data.

   </p>

The *Basic Model Interface* (BMI) is a library specification
created by the `Community Surface Dynamics Modeling System`_ (CSDMS)
to facilitate the conversion of a model or dataset
into a reusable, plug-and-play component.
Recall that in this context an interface is a named set of functions
with prescribed arguments and return values.
The BMI functions make a model self-describing and fully controllable
by a modeling framework or application.
By design, the BMI functions are straightforward to implement in
any language, using only basic data types from standard language libraries.
Also by design, the BMI functions are noninvasive.
This means that a model's BMI does not make calls to other
components or tools and is not modified to use any
framework-specific data structures. A BMI, therefore, introduces no
dependencies into a model, so the model can still be used
in a stand-alone manner.

The BMI is expressed
in the `Scientific Interface Definition Language`_ (SIDL)
in the file `bmi.sidl <./bmi.sidl>`_.
BMI specifications for five languages -- C, C++, Fortran, Java,
and Python -- are derived from this SIDL file.
For each language,
links to the specification and an example implementation
are listed in the table below.

.. table::
   :align: center
   :widths: 10, 10, 15

   ========  ==============  ======================
   Language  Specification   Example implementation
   ========  ==============  ======================
   C         `bmi-c`_        `bmi-example-c`_
   C++       `bmi-cxx`_      `bmi-example-cxx`_
   Fortran   `bmi-fortran`_  `bmi-example-fortran`_
   Java      `bmi-java`_     `bmi-example-java`_
   Python    `bmi-python`_   `bmi-example-python`_
   ========  ==============  ======================

Detailed instructions for building the specifications and examples
are given at each link above.
Alternatively, the specifications can be installed through conda
(C, C++, Fortran, Python) or Maven (Java).
See the links above for details.

While CSDMS currently supports the languages listed above,
a BMI can be written for any language.
BMI is a community-driven standard;
`contributions <CONTRIBUTING.rst>`_
that follow the `contributor code of conduct <./CODE-OF-CONDUCT.rst>`_
are welcomed,
and are `acknowledged <./AUTHORS.rst>`_.
BMI is open source software released under the `MIT License <./LICENSE>`_.
BMI is an element of the `CSDMS Workbench`_,
an integrated system of software tools, technologies, and standards
for building and coupling models.

*The Community Surface Dynamics Modeling System
is supported by the National Science Foundation.*


.. Links

.. _Community Surface Dynamics Modeling System: https://csdms.colorado.edu
.. _Scientific Interface Definition Language: http://dx.doi.org/10.1177/1094342011414036
.. _bmi-c: https://github.com/csdms/bmi-c
.. _bmi-cxx: https://github.com/csdms/bmi-cxx
.. _bmi-fortran: https://github.com/csdms/bmi-fortran
.. _bmi-java: https://github.com/csdms/bmi-java
.. _bmi-python: https://github.com/csdms/bmi-python
.. _bmi-example-c: https://github.com/csdms/bmi-example-c
.. _bmi-example-cxx: https://github.com/csdms/bmi-example-cxx
.. _bmi-example-fortran: https://github.com/csdms/bmi-example-fortran
.. _bmi-example-java: https://github.com/csdms/bmi-example-java
.. _bmi-example-python: https://github.com/csdms/bmi-example-python
.. _CSDMS Workbench: https://csdms.colorado.edu/wiki/Workbench
.. _grid_funcs:

Model grid functions
--------------------

The functions in this section describe :ref:`model grids <model_grids>`. 
In the BMI,
every :term:`exchange item` is defined on a grid,
and is referenced by a :term:`grid identifier`
returned from the :ref:`get_var_grid` function.
This identifier is a required input to the functions listed below.

A model can have multiple grids.
For example,
consider modeling the diffusion of temperature over a flat plate.
One grid could be a uniform rectilinear grid on which
temperature is defined.
A second grid could be a scalar,
on which a constant thermal diffusivity is defined.

Not all grid functions are used by each type of grid.
However, all BMI grid functions must be implemented.
(See :ref:`model_grids` and :ref:`best_practices`.)


.. _get_grid_type:

*get_grid_type*
...............

.. code-block:: java

   /* SIDL */
   int get_grid_type(in int grid, out string type);

Given a :term:`grid identifier`, get the type of that grid as a string.
Valid grid types are:

* ``scalar``
* ``points``
* ``vector``
* ``unstructured``
* ``structured_quadrilateral``
* ``rectilinear``
* ``uniform_rectilinear``

A detailed description of the grid types supported in BMI
is given in the :ref:`model_grids` section.

**Implementation notes**

* In C++, Java, and Python, the *type* argument is omitted and the grid
  type name is returned from the function as a string.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_rank:

*get_grid_rank*
...............

.. code-block:: java

   /* SIDL */
   int get_grid_rank(in int grid, out int rank);

Given a :term:`grid identifier`, get the :term:`rank` (the number of
dimensions) of that grid as an integer.

A grid's rank determines the length of the return value
of many of the following grid functions.
For instance, :ref:`get_grid_shape` returns an array of length *rank*.
Similarly, a grid's rank determines which
of :ref:`get_grid_x`, :ref:`get_grid_y`, etc. are implemented.

**Implementation notes**

* This function is needed for every :ref:`grid type <model_grids>`.
* In C++, Java, and Python, the *rank* argument is omitted and the grid
  rank is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_size:

*get_grid_size*
...............

.. code-block:: java

   /* SIDL */
   int get_grid_size(in int grid, out int size);

Given a :term:`grid identifier`,
get the total number of elements (or :term:`nodes <node>`)
of that grid as an integer.

The grid size is used for, among other things, the
length of arrays returned by :ref:`get_grid_x` and :ref:`get_grid_y`
for :ref:`unstructured <unstructured_grids>` and
:ref:`structured quad <structured_quad>` grids.

**Implementation notes**

* This function is needed for every :ref:`grid type <model_grids>`.
* In C++, Java, and Python, the *size* argument is omitted and the grid
  size is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_shape:

*get_grid_shape*
................

.. code-block:: java

   /* SIDL */
   int get_grid_shape(in int grid, in array<int, 1> shape);

Get the dimensions of the model grid.

Note that this function (as well as the other grid functions)
returns information ordered with "ij" indexing (as opposed to "xy").
For example,
consider a two-dimensional rectilinear grid
with four columns (``nx = 4``)
and three rows (``ny = 3``).
The :ref:`get_grid_shape` function would return a shape
of ``[ny, nx]``, or ``[3,4]``.
If there were a third dimension, the length of the *z*-dimension, ``nz``,
would be listed first.

Also note that the grid shape is the number of :term:`nodes <node>`
in the coordinate directions and not the number of cells or elements.
It is possible for grid values to be associated with the nodes or with
the cells.

**Implementation notes**

* This function is used for describing all :ref:`structured grids
  <structured_grids>`.
* In Python, the *shape* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_spacing:

*get_grid_spacing*
..................

.. code-block:: java

   /* SIDL */
   int get_grid_spacing(in int grid, in array<double, 1> spacing);

Get the distance between the :term:`nodes <node>` of the model grid.

The :ref:`get_grid_spacing` function provides the width of each cell in
the number of dimensions as returned by :ref:`get_grid_rank`.
As with :ref:`get_grid_shape`,
the spacing is given in "ij" indexing* order;
e.g., for a two-dimensional grid,
the spacing between rows is followed by spacing between columns, ``[dy, dx]``.

**Implementation notes**

* This function is used for describing :ref:`uniform rectilinear
  <uniform_rectilinear>` grids.
* In Python, the *spacing* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_origin:

*get_grid_origin*
.................

.. code-block:: java

   /* SIDL */
   int get_grid_origin(in int grid, in array<double, 1> origin);

Get the coordinates of the lower-left corner of the model grid.

The *origin* parameter is a one-dimensional array of the size
returned by :ref:`get_grid_rank`.
As with :ref:`get_grid_shape`,
the origin is given in "ij" indexing* order;
e.g., for a two-dimensional grid,
the origin is given in the column dimension, followed by the row dimension,
``[y0, x0]``.

**Implementation notes**

* This function is used for describing :ref:`uniform rectilinear
  <uniform_rectilinear>` grids.
* In Python, the *origin* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_x:

*get_grid_x*
............

.. code-block:: java

   /* SIDL */
   int get_grid_x(in int grid, in array<double, 1> x);

Get the locations of the grid :term:`nodes <node>` in the first
coordinate direction.

The length of the resulting one-dimensional array depends on the grid type.
(It will have either :ref:`get_grid_rank` or :ref:`get_grid_size` elements.)
See :ref:`model_grids` for more information.

**Implementation notes**

* This function is used for describing :ref:`rectilinear <rectilinear>`,
  :ref:`structured quadrilateral <structured_quad>`,
  and all :ref:`unstructured <unstructured_grids>` grids.
* In Python, the *x* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_y:

*get_grid_y*
............

.. code-block:: java

   /* SIDL */
   int get_grid_y(in int grid, in array<double, 1> y);

Get the locations of the grid :term:`nodes <node>` in the second
coordinate direction.

The length of the resulting one-dimensional array depends on the grid type.
(It will have either :ref:`get_grid_rank` or :ref:`get_grid_size` elements.)
See :ref:`model_grids` for more information.

**Implementation notes**

* This function is used for describing :ref:`rectilinear <rectilinear>`,
  :ref:`structured quadrilateral <structured_quad>`,
  and all :ref:`unstructured <unstructured_grids>` grids.
* In Python, the *y* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_z:

*get_grid_z*
............

.. code-block:: java

   /* SIDL */
   int get_grid_z(in int grid, in array<double, 1> z);

Get the locations of the grid :term:`nodes <node>` in the third
coordinate direction.

The length of the resulting one-dimensional array depends on the grid type.
(It will have either :ref:`get_grid_rank` or :ref:`get_grid_size` elements.)
See :ref:`model_grids` for more information.

**Implementation notes**

* This function is used for describing :ref:`rectilinear <rectilinear>`,
  :ref:`structured quadrilateral <structured_quad>`,
  and all :ref:`unstructured <unstructured_grids>` grids.
* In Python, the *z* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_node_count:

*get_grid_node_count*
.....................

.. code-block:: java

   /* SIDL */
   int get_grid_node_count(in int grid, out int count);

Get the number of :term:`nodes <node>` in the grid.

**Implementation notes**

* This function is used for describing :ref:`unstructured
  <unstructured_grids>` grids.
* In C++, Java, and Python, the *count* argument is omitted and the node
  count is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_edge_count:

*get_grid_edge_count*
.....................

.. code-block:: java

   /* SIDL */
   int get_grid_edge_count(in int grid, out int count);

Get the number of :term:`edges <edge>` in the grid.

**Implementation notes**

* This function is used for describing :ref:`unstructured
  <unstructured_grids>` grids.
* In C++, Java, and Python, the *count* argument is omitted and the edge
  count is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_face_count:

*get_grid_face_count*
.....................

.. code-block:: java

   /* SIDL */
   int get_grid_face_count(in int grid, out int count);

Get the number of :term:`faces <face>` in the grid.

**Implementation notes**

* This function is used for describing :ref:`unstructured
  <unstructured_grids>` grids.
* In C++, Java, and Python, the *count* argument is omitted and the face
  count is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_edge_nodes:

*get_grid_edge_nodes*
.....................

.. code-block:: java

   /* SIDL */
   int get_grid_edge_nodes(in int grid, in array<int, 1> edge_nodes);

Get the edge-node connectivity.

For each edge, connectivity is given as node at edge tail, followed by
node at edge head. The total length of the array is 
2 * :ref:`get_grid_edge_count`.

**Implementation notes**

* This function is used for describing :ref:`unstructured
  <unstructured_grids>` grids.
* In Python, the *edge_nodes* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_face_edges:

*get_grid_face_edges*
.....................

.. code-block:: java

   /* SIDL */
   int get_grid_face_edges(in int grid, in array<int, 1> face_edges);

Get the face-edge connectivity.

The length of the array returned is the sum of the values of
:ref:`get_grid_nodes_per_face`.

**Implementation notes**

* This function is used for describing :ref:`unstructured
  <unstructured_grids>` grids.
* In Python, the *face_edges* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_face_nodes:

*get_grid_face_nodes*
.....................

.. code-block:: java

   /* SIDL */
   int get_grid_face_nodes(in int grid, in array<int, 1> face_nodes);

Get the face-node connectivity.

For each face, the nodes (listed in a counter-clockwise direction)
that form the boundary of the face.
For a grid of quadrilaterals, 
the total length of the array is 4 * :ref:`get_grid_face_count`.
More generally,
the length of the array is the sum of the values of
:ref:`get_grid_nodes_per_face`.

**Implementation notes**

* This function is used for describing :ref:`unstructured
  <unstructured_grids>` grids.
* In Python, the *face_nodes* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]


.. _get_grid_nodes_per_face:

*get_grid_nodes_per_face*
.........................

.. code-block:: java

   /* SIDL */
   int get_grid_nodes_per_face(in int grid, in array<int, 1> nodes_per_face);

Get the number of nodes for each face.

The returned array has a length of :ref:`get_grid_face_count`.
The number of edges per face is equal to the number of nodes per face.

**Implementation notes**

* This function is used for describing :ref:`unstructured
  <unstructured_grids>` grids.
* In Python, the *nodes_per_face* argument is a :term:`numpy <NumPy>` array.
* In C++ and Java, this is a void function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`grid_funcs` | :ref:`basic_model_interface`]
.. include:: ../../CONTRIBUTING.rst
.. include:: ../../CODE-OF-CONDUCT.rst
.. _pymt:

BMI-based tools: *babelizer* and *pymt*
=======================================

A model that provides the BMI functions can be converted
to a plug-and-play component 
that runs in the CSDMS `Python Modeling Toolkit`_ (*pymt*).
This conversion process can be done
automatically with the CSDMS `Babelizer`_, which generates
a Python package from models written in
C, C++, Fortran, and Python.
Some additional metadata is required,
describing such things as:

*  the model author(s), license, description, web page, etc.
*  how the model is built and how it can be linked
*  template input files
*  description of input parameters (description, units, ranges, etc.)

Within *pymt*, a component automatically gains many
new capabilities. This includes the ability to be coupled to
other models even if their programming language, variable
names, variable units, time-stepping scheme or
computational grid is different.  It also gains:

* the ability to write output variables to standardized NetCDF
  files,
* unit conversion tools,
* time interpolation,
* all the data analysis and visualization tools available in Python,
* the ability to run in a Jupyter Notebook

If you have a model with a BMI,
and would like to componentize it and add it to *pymt*,
please contact us through the `CSDMS Help Desk`_.
We'd be happy (and excited) to help!


.. Links:

.. _Python Modeling Toolkit: https://pymt.readthedocs.io
.. _Babelizer: https://babelizer.readthedocs.io
.. _CSDMS Help Desk: https://github.com/csdms/help-desk
Steering Council
================
The Basic Model Interface (BMI): governance and decision-making
===============================================================

The purpose of this document is to formalize the governance process used for the
Basic Model Interface (BMI) in both ordinary and extraordinary situations, and
to clarify how decisions are made and how the various elements of our community
interact, including the relationship between open source collaborative
development and work that may be funded by for-profit or non-profit entities.

Summary
-------

The Basic Model Interface (BMI) is a community-owned and community-run project.
To the maximum extent possible, decisions about project direction are made by
community consensus (but note that "consensus" here has a somewhat technical
meaning that might not match everyone's expectations--see below). Some members
of the community additionally contribute by serving on the Steering Council (see
below), where they are responsible for facilitating the establishment of
community consensus, for stewarding project resources, and--in extreme
cases--for making project decisions if the normal community-based process breaks
down.

The Project
-----------

BMI is an open source software project (hereafter, the Project) affiliated with
the NSF-funded `Community Surface Dynamics Modeling System`_ (CSDMS). The goal
of the Project is to develop an open source software interface standard for
querying and controlling models. The Project also includes tools, documentation,
and examples to support and promote this standard. The software developed by the
Project is released under the MIT open source license, developed openly, and
hosted on public GitHub repositories under the csdms and other GitHub
organization. Proposed changes to the Project must follow the `CONTRIBUTING`_
document in the Project's main GitHub repository.

The Project is developed by a team of distributed developers, called
Contributors. Contributors are individuals who have contributed code,
documentation, designs, or other work to the Project. Anyone can be a
Contributor. Contributors can be affiliated with any legal entity or none.
Contributors participate in the project by submitting, reviewing, and discussing
GitHub pull requests and issues and participating in open and public Project
discussions on GitHub and other channels. The foundation of Project
participation is openness and transparency.

The Project Community consists of all Contributors and Users of the Project.
Contributors work on behalf of and are responsible to the larger Project
Community. We strive to keep the barrier between Contributors and Users as low
as possible.

Governance
----------

This section describes the governance and leadership model of the Project. The
foundations of Project governance are:

* Openness & Transparency
* Active Contribution
* Institutional Neutrality

Consensus-based decision making by the community
................................................

Normally, decisions on the Project will be made by consensus of all interested
Contributors. The primary goal of this approach is to ensure that the people who
are most affected by and involved in any given change can contribute their
knowledge in the confidence that their voices will be heard, because thoughtful
review from a broad community is the best mechanism we know of for creating
high-quality software.

The mechanism we use to accomplish this goal may be unfamiliar for those who are
not experienced with the cultural norms around free/open-source software
development. We provide a summary here, and highly recommend that all
Contributors additionally read the chapter `Social and Political
Infrastructure`_ of Karl Fogel's *Producing Open Source Software*, and in
particular the section on `Consensus-based Democracy`_, for a more detailed
discussion.

In this context, consensus does not require:

* that we wait to solicit everyone's opinion on every change,
* that we ever hold a vote on anything, or
* that everyone is happy or agrees with every decision.

For us, what consensus means is that we entrust everyone with the right to veto
any change if they feel it necessary. While this may sound like a recipe for
obstruction, this is not what happens. Instead, we find that most people take
this responsibility seriously, and only invoke their veto when they judge that a
serious problem is being ignored, and that their veto is necessary to protect
the Project. In practice, it turns out that vetoes are almost never
formally invoked because their mere possibility ensures that Contributors are
motivated from the start to find solutions that everyone can live with, thus
accomplishing our goal of ensuring that all interested perspectives are taken
into account.

How do we know when consensus has been achieved? In principle, this is
difficult, since consensus is defined by the absence of vetoes, which requires
us to somehow prove a negative. In practice, we use a combination of our best
judgment (e.g., a simple and uncontroversial bug fix posted on GitHub and
reviewed by a core developer is probably fine) and best efforts (e.g.,
substantive changes must not only follow the Project's CONTRIBUTING document,
but also be posted to the Project's designated communication channel in order to
give the broader community a chance to catch any problems and suggest
improvements; we assume that anyone who cares enough about BMI to invoke their
veto right should be on the Project's communication channel). If no one comments
on the Project's communication channel after several days, then it's probably
fine.

If one does need to invoke a formal veto, then it should consist of:

* an unambiguous statement that a veto is being invoked,
* an explanation of why it is being invoked, and
* a description of what conditions (if any) would convince the vetoer to
  withdraw their veto.

If all proposals for resolving some issue are vetoed, then the status quo wins
by default.

In the worst case, if a Contributor is genuinely misusing their veto in an
obstructive fashion to the detriment of the Project, then they can be ejected
from the Project by consensus of the Steering Council--see below.

Steering Council
................

The Project has a Steering Council (a.k.a. the BMI Council) that consists of
Project Contributors and Users. The overall role of the Council is to ensure,
with input from the Community, the long-term well-being of the Project, both
technically and as a community.

During everyday Project activities, Council Members participate in discussions,
code reviews, and other Project activities as peers with all other Contributors
and the Community. In these everyday activities, Council Members do not have any
special power or privilege through their membership on the Council. However, it
is expected that because of the quality and quantity of their contributions and
their expert knowledge of the Project that Council Members will provide useful
guidance, both technical and in terms of Project direction, to potentially less
experienced Contributors.

The Steering Council plays a special role in certain situations. In particular,
the Council may, if necessary:

* Make decisions about the overall scope, vision, and direction of the Project.
* Make decisions about strategic collaborations with other organizations or
  individuals.
* Make decisions about specific technical issues, features, bugs, and pull
  requests. They are the primary mechanism of guiding the code review process and
  merging pull requests.
* Update policy documents such as this one.
* Make decisions when regular community discussion doesn’t produce consensus on
  an issue in a reasonable time frame.

However, the Council's primary responsibility is to facilitate the ordinary
community-based decision making procedure described above. If the Council ever
has to step in and formally override the community for the health of the
Project, then they will do so, but they will consider reaching this point to
indicate a failure in their leadership.

Council decision making
.......................

If it becomes necessary for the Steering Council to produce a formal decision,
then they will use a form of the `Apache Foundation voting process`_. This is a
formalized version of consensus, in which +1 votes indicate agreement, -1 votes
are vetoes (and must be accompanied with a rationale, as above), and fractional
votes (e.g. -0.5, +0.5) can be used if one wishes to express an opinion without
registering a full veto. These numeric votes can also be used informally to get
a general sense of the Community's feelings on some issue. A formal vote only
occurs if explicitly declared, and if this does occur then the vote should be
held open for long enough to give all interested Council Members a chance to
respond--at least one week.

In practice, we anticipate that for most Council decisions (e.g., voting in new
members) a more informal process will suffice.

Council membership
..................

A list of current Steering Council Members is maintained at the page `Steering
Council`_.

To become eligible to join the Steering Council, an individual must be a Project
Contributor who has produced substantial contributions or a Project User that
has applied BMI in a substantial way. Candidate Council Members are nominated by
existing Council Members. The Candidate must confirm they are interested and
willing to serve in this capacity. The Candidate becomes a Member following
consensus of the existing Council. The Council will be initially formed from a
set of existing Project Contributors and Users who, as of early 2022, have been
currently active in Project development or application.

When considering potential Members, the Council will look at candidates with a
comprehensive view, including but not limited to code, code review,
applications, infrastructure work, communication channel participation,
community help/building, education and outreach, design work, etc. We are
deliberately not setting arbitrary quantitative metrics (like "100 commits in
this repo") to avoid encouraging behavior that plays to the metrics rather than
the Project's overall well-being. We want to encourage a diverse array of
backgrounds, viewpoints, and talents, which is why we explicitly do not define
code as the sole metric on which Council membership will be evaluated.

If a Council Member becomes inactive in the Project for a period of one year,
they will be considered for removal from the Council. Before removal, the
inactive Member will be approached to see if they plan on returning to active
participation. If not, they will be removed after a Council vote. If they plan
on returning to active participation, they will be given a grace period of one
year. If they do not return to active participation within that time period they
will be removed by vote of the Council without further grace period. All former
Council Members can be considered for membership again at any time in the
future, like any other Project Contributor or User. Retired Council members will
be listed on the project website, acknowledging the period during which they
were active in the Council.

The Council reserves the right to eject current Members if they are deemed to be
actively harmful to the Project's well-being, and if attempts at communication
and conflict resolution have failed. This requires the consensus of the
remaining Members.

Conflict of interest
....................

It is expected that Council Members will be employed at a range of universities,
government agencies, companies, and non-profit organizations. Because of this,
it is possible that Members will have conflict of interests. Such conflict of
interests include, but are not limited to:

* Financial interests, such as investments, employment or contracting work,
  outside of the Project that may influence their work on the Project.
* Access to proprietary information of their employer that could potentially leak
  into their work with the Project.

All members of the Council shall disclose to the rest of the Council any
conflict of interest they may have. Members with a conflict of interest in a
particular issue may participate in Council discussions on that issue, but must
recuse themselves from voting on the issue.

Private communications of the Council
.....................................

To the maximum extent possible, Council discussions and activities will be
public and done in collaboration and discussion with the Project Contributors
and Community. The Council will have a private communication channel that will
be used sparingly and only when a specific matter requires privacy. When private
communications and decisions are needed, the Council will do its best to
summarize those to the Community after eliding personal/private/sensitive
information that should not be posted to the public internet.

Subcommittees
.............

The Council can create subcommittees that provide leadership and guidance for
specific aspects of the Project. Like the Council as a whole, subcommittees
should conduct their business in an open and public manner unless privacy is
specifically called for. Private subcommittee communications should happen on
the communication channel of the Council unless specifically called for.

Institutional Partners and Funding
----------------------------------

The Steering Council is the primary leadership for the Project. No outside
institution, individual, or legal entity has the ability to own, control, usurp
or influence the Project other than by participating in the Project as
Contributors and Council Members. However, because institutions can be an
important funding mechanism for the project, it is important to formally
acknowledge institutional participation in the Project. These are Institutional
Partners.

An Institutional Contributor is any individual Project Contributor who
contributes to the project as part of their official duties as an Institutional
Partner. Likewise, an Institutional Council Member is any Project Steering
Council Member who contributes to the Project as part of their official duties
as an Institutional Partner.

With these definitions, an Institutional Partner is any recognized legal entity
in the United States or elsewhere that employs at least one Institutional
Contributor or Institutional Council Member. Institutional Partners can be
for-profit or non-profit entities.

Institutions become eligible to become an Institutional Partner by employing
individuals who actively contribute to the Project as part of their official
duties. To state this another way, the only way for a Partner to influence the
project is by actively contributing to the open development of the Project, in
equal terms to any other Contributor or Council Member. Merely using Project
software in an institutional context does not allow an entity to become an
Institutional Partner. Financial gifts do not enable an entity to become an
Institutional Partner. Once an institution becomes eligible for Institutional
Partnership, the Steering Council must nominate and approve the Partnership.

If at some point an existing Institutional Partner stops having any contributing
employees, then a one-year grace period commences. If at the end of this one
year period they continue not to have any contributing employees, then their
Institutional Partnership will lapse, and resuming it will require going through
the normal process for new Partnerships.

An Institutional Partner is free to pursue funding for their work on the Project
through any legal means. This could involve a non-profit organization raising
money from private foundations and donors or a for-profit company building
proprietary products and services that leverage Project software and services.
Funding acquired by Institutional Partners to work on the Project is called
Institutional Funding. However, no funding obtained by an Institutional Partner
can override the Steering Council. If a Partner has funding to do Project work
and the Council decides to not pursue that work, the Partner is free to pursue
it on their own. However in this situation, that part of the Partner's work will
not be under the Project umbrella and cannot use the Project trademarks in a way
that suggests a formal relationship.

Institutional Partner benefits are:

* Acknowledgement on Project websites, in talks, and on promotional material.
* Ability to acknowledge their own funding sources on Project websites, in
  talks, and on promotional material.
* Ability to influence the Project through the participation of their Council
  Member.

A list of current Institutional Partners is maintained at the page 
`Institutional Partners`_.

Document history
----------------

https://github.com/csdms/bmi/commits/master/docs/source/governance.rst

Acknowledgements
----------------

Substantial portions of this document were adapted from the `NumPy governance
document`_.

License
-------

To the extent possible under law, the authors have waived all copyright and
related or neighboring rights to the BMI project governance and decision-making
document, as per the `CC-0 public domain dedication / license`_.



.. Links

.. _Community Surface Dynamics Modeling System: https://csdms.colorado.edu
.. _CONTRIBUTING: https://github.com/csdms/bmi/blob/master/CONTRIBUTING.rst
.. _Chapter 4: Social and Political Infrastructure
.. _Social and Political Infrastructure: http://producingoss.com/en/producingoss.html#social-infrastructure
.. _Consensus-based Democracy: http://producingoss.com/en/producingoss.html#consensus-democracy
.. _Apache Foundation voting process: https://www.apache.org/foundation/voting.html
.. _Steering Council: ./council.html
.. _Institutional Partners: ./partners.html
.. _NumPy governance document: https://numpy.org/doc/stable/dev/governance/index.html
.. _CC-0 public domain dedication / license: https://creativecommons.org/publicdomain/zero/1.0/
.. include:: ../../AUTHORS.rst
.. _model_grids:

Model grids
===========

The Basic Model Interface (BMI) supports several different :term:`grid` types,
described below.
Depending on the grid type
(as returned from the :ref:`get_grid_type` function),
a model will implement a different set of :ref:`grid functions <grid_funcs>`.


.. _structured_grids:

Structured grids
----------------

For the BMI specification,
a structured grid is formed in two dimensions by tiling a domain
with quadrilateral cells so that every interior vertex
is surrounded by four quadrilaterals.
In three dimensions, six-sided polyhedral cells are stacked such that
interior vertices are surrounded by eight cells.

In the BMI,
dimensional information is ordered with "ij" indexing
(as opposed to "xy").
For example,
for the uniform rectilinear grid shown below,
the :ref:`get_grid_shape` function would return the array ``[4, 5]``.
If there was a third dimension,
its length would be listed first.

.. note::

  The grid shape is the number of rows and columns of :term:`nodes
  <node>`, as opposed to other types of element (such as cells or
  faces). It is possible for grid values to be associated with the
  nodes or with the cells.


.. _uniform_rectilinear:

Uniform rectilinear
^^^^^^^^^^^^^^^^^^^

.. image:: images/mesh_uniform_rectilinear.png
   :scale: 20 %

A uniform rectilinear grid is a special case of structured quadrilateral grid
where the elements have equal width in each dimension.
That is, for a two-dimensional grid, elements have a constant width
of ``dx`` in the *x*-direction and ``dy`` in the *y*-direction.
The case of ``dx == dy`` is oftentimes called
a *raster* or *Catesian grid*.

To completely specify a uniform rectilinear grid,
only three pieces of information are needed:
the number of elements in each dimension,
the width of each element (in each dimension),
and the location of the corner of the grid.

Uniform rectilinear grids use the following BMI functions:

* :ref:`get_grid_rank`
* :ref:`get_grid_size`
* :ref:`get_grid_shape`
* :ref:`get_grid_spacing`
* :ref:`get_grid_origin`


.. _rectilinear:

Rectilinear
^^^^^^^^^^^

.. image:: images/mesh_rectilinear.png
   :scale: 20 %

In a rectilinear grid, the spacing between nodes in each dimension varies,
as depicted above.
Therefore,
an array of coordinates for each row and column
(for the two-dimensional case) is required.

The :ref:`get_grid_y` function provides an array (whose length is the number of
*rows*) that gives the *y*-coordinate for each row.

The :ref:`get_grid_x` function provides an array (whose length is the number of
*columns*) that gives the *x*-coordinate for each column.

Rectilinear grids use the following BMI functions:

* :ref:`get_grid_rank`
* :ref:`get_grid_size`
* :ref:`get_grid_shape`
* :ref:`get_grid_x`
* :ref:`get_grid_y`
* :ref:`get_grid_z`


.. _structured_quad:

Structured quadrilateral
^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: images/mesh_structured_quad.png
   :scale: 20 %

The most general structured quadrilateral grid is one where
the rows (and columns) do not share a common coordinate. In this
case, coordinates are required for each grid node. For this
more general case, :ref:`get_grid_x` and :ref:`get_grid_y` are
repurposed to provide this information.

The :ref:`get_grid_y` function returns an array (whose length is the number
of total nodes returned by :ref:`get_grid_size`) of *y*-coordinates.

The :ref:`get_grid_x` function returns an array (whose length is the number
of total nodes returned by :ref:`get_grid_size`) of *x*-coordinates.

Structured quadrilateral grids use the following BMI functions:

* :ref:`get_grid_rank`
* :ref:`get_grid_size`
* :ref:`get_grid_shape`
* :ref:`get_grid_x`
* :ref:`get_grid_y`
* :ref:`get_grid_z`


.. _unstructured_grids:

Unstructured grids
------------------

.. image:: images/mesh_unstructured.png
   :scale: 25 %

This category includes the *unstructured* type,
as well as the special cases
*scalar*, *points*, and *vector*.
This is the most general grid type.
It can be used for any type of grid.
This grid type must be used if the grid consists of cells
that are not quadrilaterals;
this includes any grid of triangles (e.g. `Delaunay triangles`_
and `Voronoi tesselations`_).

.. note::

   A grid of `equilateral triangles`_, while they are most certainly
   *structured*, would need to be represented as an unstructured grid.
   The same is true for a grid of `hexagons`_.


BMI uses the `ugrid conventions`_ to define unstructured grids.

Unstructured grids use the following BMI functions:

* :ref:`get_grid_rank`
* :ref:`get_grid_x`
* :ref:`get_grid_y`
* :ref:`get_grid_z`
* :ref:`get_grid_node_count`
* :ref:`get_grid_edge_count`
* :ref:`get_grid_face_count`
* :ref:`get_grid_edge_nodes`
* :ref:`get_grid_face_edges`
* :ref:`get_grid_face_nodes`
* :ref:`get_grid_nodes_per_face`

For a demonstration of how these BMI functions work,
let's use the unstructured grid in the annotated figure above.

The grid is two-dimensional,
so the :ref:`get_grid_rank` function returns 2.

The :term:`nodes <node>` of the grid, labeled in the figure in red,
are given by coordinates

.. code-block:: python

   x = [0, 1, 2, 1, 3, 4]
   y = [3, 1, 2, 4, 0, 3]

These will be the outputs of the :ref:`get_grid_x` and 
:ref:`get_grid_y` functions, respectively.
The nodes are indexed, so 
node 0 is at *(x, y) = (0, 3)*,
node 1 is at *(x, y) = (1, 1)*, etc.

As with the grid nodes,
the grid :term:`edges <edge>` and :term:`faces <face>` are indexed.
In the figure,
the edges are depicted in blue italics,
while the faces are boldfaced. 
The outputs from :ref:`get_grid_node_count`, :ref:`get_grid_edge_count`,
and :ref:`get_grid_face_count` are:

.. code-block:: python

   node_count = 6
   edge_count = 8
   face_count = 3

Note that the number of nodes is the length of the *x* and *y* vectors above.

The :ref:`get_grid_nodes_per_face` function returns a vector
of length `face_count`.
The first two faces are quadrilaterals,
while the third is a triangle, so

.. code-block:: python

   nodes_per_face = [4, 4, 3]

The :ref:`get_grid_edge_nodes` function returns a vector
of length `2*edge_count`.
The vector is formed, pairwise,
by the node index at the tail of the edge,
followed by the node index at the head of the edge.
For the grid in the figure, this is

.. code-block:: python

   edge_nodes = [0, 1, 1, 2, 2, 3, 3, 0, 1, 4, 4, 5, 5, 2, 5, 3]

The :ref:`get_grid_face_edges` function returns a vector
of length `sum(nodes_per_face)`.
The vector is formed from the edge indices as displayed in the figure:

.. code-block:: python

   face_edges = [0, 1, 2, 3, 4, 5, 6, 1, 6, 7, 2]

Likewise, the :ref:`get_grid_face_nodes` function returns a vector
of length `sum(nodes_per_face)`.
The vector is formed from the node indices as displayed in the figure:

.. code-block:: python

   face_nodes = [0, 1, 2, 3, 1, 4, 5, 2, 2, 5, 3]



.. Links

.. _Delaunay triangles: http://en.wikipedia.org/wiki/Delaunay_triangulation
.. _Voronoi tesselations: http://en.wikipedia.org/wiki/Voronoi_tessellation
.. _equilateral triangles: http://en.wikipedia.org/wiki/Triangle_tiling
.. _hexagons: http://en.wikipedia.org/wiki/Hexagonal_tiling
.. _ugrid conventions: http://ugrid-conventions.github.io/ugrid-conventions
.. _var_funcs:

Variable information functions
------------------------------

These BMI functions provide information
about a particular input or output variable.
They must accommodate any variable returned from the
:ref:`get_input_var_names` or :ref:`get_output_var_names` functions --
the variable name is used as an argument in each function. 
Based on the information returned,
type or unit conversions can be applied when necessary.


.. _get_var_grid:

*get_var_grid*
..............

.. code-block:: java

   /* SIDL */
   int get_var_grid(in string name, out int grid);

Each input and output variable is defined on a grid.
(Available grid types are listed in the :ref:`grid_funcs` section.)
The `get_var_grid` function provides the identifier (an integer) for this grid.
The identifier can be passed to the BMI
:ref:`grid information <grid_funcs>` functions
to get the details of a particular grid;
e.g., *x*- and *y*-coordinates, size, type, etc.
A model can have one or more grids.

**Implementation notes**

* Grid identifiers start at 0.
* In C++, Java, and Python, the *grid* argument is omitted and the grid
  identifier is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or
  failure (nonzero) is returned.

[:ref:`var_funcs` | :ref:`basic_model_interface`]


.. _get_var_type:

*get_var_type*
..............

.. code-block:: java

   /* SIDL */
   int get_var_type(in string name, out string type);

The `get_var_type` function provides the data type of the
variable as it's stored in memory by the model.
The data type is returned as a string.
Use of native language type names is encouraged;
e.g., in C, use `int`, `float`, and `double`,
while in Fortran, use `integer`, `real`, and `double precision`.

**Implementation notes**

* In C++, Java, and Python, the *type* argument is omitted and the variable
  type name is returned from the function as a string.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.
* In Java, only `primitive types`_ (e.g., ``int``, ``double``), not
  `wrapper classes`_ (e.g., ``Integer``, ``Double``), are supported.

[:ref:`var_funcs` | :ref:`basic_model_interface`]


.. _get_var_units:

*get_var_units*
...............

.. code-block:: java

   /* SIDL */
   int get_var_units(in string name, out string units);

Get the units of the given variable.
Standard unit names, in lower case, should be used,
such as ``"meters"`` or ``"seconds"``.
Standard abbreviations, such as ``"m"`` for meters, are
also supported. For variables with compound units, each unit name
is separated by a single space, with exponents other than 1 placed
immediately after the name, as in ``"m s-1"`` for velocity,
``"W m-2"`` for an energy flux, or ``"km2"`` for an area.
The abbreviations used in the BMI are derived from
Unidata's `UDUNITS`_ package.
See, for example, `The Units Database`_ for a
full description of valid unit names and a list of supported units.

**Implementation notes**

* Dimensionless quantities should use ``""`` or ``"1"`` as the unit.
* Variables without units should use ``"none"``.
* In C++, Java, and Python, the *units* argument is omitted and the variable
  units name is returned from the function as a string.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`var_funcs` | :ref:`basic_model_interface`]


.. _get_var_itemsize:

*get_var_itemsize*
..................

.. code-block:: java

   /* SIDL */
   int get_var_itemsize(in string name, out int size);

The `get_var_itemsize` function provides the size, in bytes,
of a single element of the variable.
For example, if data for a variable are stored as 64-bit integers,
`get_var_itemsize` would return 8.

**Implementation notes**

* In C++, Java, and Python, the *size* argument is omitted and the item size
  is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`var_funcs` | :ref:`basic_model_interface`]


.. _get_var_nbytes:

*get_var_nbytes*
................

.. code-block:: java

   /* SIDL */
   int get_var_nbytes(in string name, out int nbytes);

The `get_var_nbytes` function provides the total amount of memory used to store
a variable; i.e., the number of items multiplied by the size of each item.

**Implementation notes**

* In C++, Java, and Python, the *nbytes* argument is omitted and the total
  amount of memory used by the variable is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`var_funcs` | :ref:`basic_model_interface`]


.. _get_var_location:

*get_var_location*
..................

.. code-block:: java

   /* SIDL */
   int get_var_location(in string name, out string location);

The `get_var_location` function,
given a variable name, returns a string that indicates on what grid
element the variable is defined. Valid return values are:

* ``node``
* ``edge``
* ``face``

**Implementation notes**

* In C++, Java, and Python, the *location* argument is omitted and the location
  is returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.
* If the given variable is a scalar (i.e., defined on a :ref:`scalar
  grid <unstructured_grids>`), the location from this function is ignored.

[:ref:`var_funcs` | :ref:`basic_model_interface`]


.. Links

.. _dtype: https://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html
Glossary
========

A glossary of terms used with BMI.


.. glossary::

   $

      The default shell prompt.

   Anaconda

      A Python distribution that includes libraries for scientific
      computing and a package manager. See
      https://www.anaconda.com/distribution for more information.

   Basic Model Interface

      A set a functions that are used to interact with and control a
      model. See https://bmi.readthedocs.io for more information.

   BMI

      See :term:`Basic Model Interface`.

   Community Surface Dynamics Modeling System

      CSDMS is an NSF-funded program that seeks to transform the
      science and practice of earth-surface dynamics modeling. For
      more information, visit https://csdms.colorado.edu.

   class

      A program that acts as a template for creating
      :term:`objects<object>`.

   conda

      The package manager for :term:`Anaconda`. Also an informal name
      for an Anaconda installation.

   conda-forge

      A collection of community-built packages distributed by
      :term:`Anaconda`. See https://conda-forge.org.

   conda environment

      A :term:`conda` sub-installation that isolates a group of
      packages from the main conda installation.

   configuration file

      A file that contains information, including initial values of
      parameters, for setting up a :term:`model`.

   coupling

      See :term:`model coupling`.

   CSDMS

      See :term:`Community Surface Dynamics Modeling System`.

   CSDMS Workbench

      An integrated system of software tools, technologies, and
      standards for building and coupling models. See
      https://csdms.colorado.edu/wiki/Workbench for more information.

   data

      Information held by an :term:`object`.

   edge

      A line or curve bounded by two :term:`nodes <node>`.

   exchange item

      A variable that a model with a BMI produces or consumes.
      Exchange items are public, not internal variables in the
      model. Exchange items should use :term:`Standard Names`.

   face

      A plane or surface enclosed by a set of :term:`edges <edge>`. In
      a 2D horizontal application one may consider the word "polygon",
      but in the hierarchy of elements the term "face" is most common.
      See also :term:`node`.

   grid

      A representation of a larger spatial domain by smaller discrete
      cells. See :ref:`model_grids` and :ref:`references`, as well as
      terms :term:`node`, :term:`edge`, and :term:`face`.

   grid identifier

      A unique object that labels (identifies) a model grid. Grid
      identifiers are integers, starting at zero. Often abbreviated
      "grid id". They're obtained through the :ref:`get_var_grid`
      function.

   grid node

      See :term:`node`.

   import

      The process of bringing code from a Python :term:`module` into
      another module or into an interactive Python session.

   instance

      See :term:`object`.

   method

      Programs that act upon the :term:`data` of an :term:`object`.

   model

      A computer program that attempts to describe a physical process
      with mathematical relationships that evolve over time and are
      solved numerically. For more information, see, for example,
      https://en.wikipedia.org/wiki/Numerical_modeling_(geology).

   model configuration file

      A file, usually in a text-based format, that lists the tunable
      parameters of a model and supplies their initial values.

   model coupling

      Models are *coupled* when they exchange inputs and outputs,
      often at the resolution of individual time steps. *One-way
      coupling* occurs when the outputs from one model are used as
      inputs to another model. *Two-way coupling* is when outputs from
      one model are used as inputs for another model, which in turn
      supplies its outputs to the first model as inputs, producing a
      feedback.

   module

      A file (with the ``.py`` extension) that contains Python code.

   node

      A point that has a coordinate pair or triplet: the most basic
      element of a grid. Variable values are typically calculated at
      nodes. See also :term:`face` and :term:`edge`.

   NumPy

      A Python library that provides arrays. Outputs from *pymt* are
      NumPy arrays. See also http://www.numpy.org.

   object

      A variable that is a concrete example of a
      :term:`class`. Objects have :term:`data` and
      :term:`methods<method>` that act upon those data.

   rank

      The number of dimensions of a model grid. A scalar has rank 0, a
      vector has rank 1, a rectilinear grid has rank 2, etc.

   refactor

      The act of modifying the internals of a program without changing
      the external behaviors of the program. Refactoring is often done
      to clean up code and improve its performance.

   Scientific Interface Definition Language

      A specification language for describing software interfaces to
      scientific model codes. See :ref:`references`.

   SIDL

      See :term:`Scientific Interface Definition Language`.

   Standard Names

      A semantic mediation technology developed at CSDMS for precisely
      matching variable names between models. For more information,
      see https://csdms.colorado.edu/wiki/CSDMS_Standard_Names.

   unit test

      A program that isolates and runs a section (a unit) of source
      code to ensure that it produces an expected result.
Institutional Partners
======================

* Community Surface Dynamics Modeling System
.. _best_practices:

BMI best practices
==================

BMI is a simple concept---it's just a set of functions
with rules for the function names, arguments, and returns.
However, when implementing a BMI, the devil is often in the details.
In no particular order,
here are some tips to help when writing a BMI for a model.

* All functions in the BMI must be implemented. For example, even if a
  model operates on a :ref:`uniform rectilinear <uniform_rectilinear>`
  grid, a :ref:`get_grid_x` function has to be written. This function
  can be empty and simply return the ``BMI_FAILURE`` status code or
  raise a ``NotImplemented`` exception, depending on the language.

* The :ref:`BMI functions <basic_model_interface>` listed in the
  documentation are the minimum required. Optional functions that act
  as helpers can be added to a model's BMI. For example, an
  `update_frac` function that updates a model's state by a fractional
  time step is a common addition to a BMI.

* Implementation details are left up to the model developer. All that
  the BMI specifies are the names of the functions, their arguments,
  and their return values.

* :term:`Standard Names` are not required for naming a model's
  :term:`exchange items <exchange item>`. However, the use of
  standardized names makes it easier for a framework (or a human) to
  match input and output variables between models.

* Don't change the variable names for exchange items 
  you currently use within your model
  to :term:`Standard Names`. Instead, find a
  `matching`_ Standard Name for each variable and then
  write your BMI functions to accept the Standard Names and map them
  to your model's internal names.

* Constructs and features that are natural for the language should be
  used when implementing a BMI. BMI strives to be developer-friendly.

* BMI functions always use flattened, one-dimensional arrays. This
  avoids any issues stemming from row/column-major indexing when
  coupling models written in different languages. It's the developer's
  responsibility to ensure that array information is
  flattened/redimensionalized in the correct order.

* Recall that models can have mulitple grids. This can be particularly
  useful for defining :term:`exchange items <exchange item>` that
  don't vary over the model domain; e.g., a diffusivity -- just define
  the variable on a separate :ref:`scalar grid <unstructured_grids>`.

* Avoid using global variables, if possible. This isn't strictly a BMI
  requirement, but if a model only uses local variables, its BMI will
  be self-contained. This may allow multiple instances of the model to
  be run simultaneously, possibly permitting the model to be coupled
  with itself.

* Boundary conditions, including boundary conditions that change with
  the model state, can be represented with :term:`exchange items
  <exchange item>`.

* :term:`Configuration files <configuration file>` are typically text
  (e.g., `YAML`_ is preferred by CSDMS), but this isn't a strict
  requirement; binary data models such as `HDF5`_ and `netCDF`_ also
  work well for storing configuration data.

* Before fitting a model with a BMI, the model code may have to be
  :term:`refactored <refactor>` into modular *initialize-run-finalize*
  (IRF) steps in order to interact with the BMI functions. This is often
  the most difficult part of adding a BMI, but the modularization
  tends to improve the quality of the code.

* Avoid allocating memory within a BMI function. Memory allocation is
  typically the responsibility of the model. This helps keep the BMI
  middleware layer thin.

* Be sure to check out the examples: `C`_, `C++`_, `Fortran`_, `Java`_, `Python`_.
  Although they wrap a very simple model, they give useful insights into how a
  BMI can be implemented in each language.

* Return codes (C and Fortran) and exceptions (C++, Java, and Python) can help with
  debugging a BMI, and can provide useful information to a user.

.. Links:

.. _YAML: https://yaml.org/
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _netCDF: https://www.unidata.ucar.edu/software/netcdf/
.. _C: https://github.com/csdms/bmi-example-c
.. _C++: https://github.com/csdms/bmi-example-cxx
.. _Fortran: https://github.com/csdms/bmi-example-fortran
.. _Java: https://github.com/csdms/bmi-example-java
.. _Python: https://github.com/csdms/bmi-example-python
.. _matching: https://github.com/csdms/standard_names_registry
When you climb in the driver's seat of an unfamiliar car,
you are nonetheless presented with a familiar sight.
Whatever the make or model may be,
we take it for granted that the vehicle will provide
a steering wheel, brake pedal, and speedometer,
alongside the various other controls and readouts
that are common to essentially all cars and trucks on the planet.
Although we don't usually think of it this way,
drivers across the globe benefit from a standard interface:
a set of control mechanisms and information displays
that have essentially the same design regardless of whether the vehicle
is a tiny electric two-seater or a giant stretch limousine.
This standard interface makes the job of operating a vehicle much easier
than it would be if each one presented a radically different interface.
Imagine a world where switching from a sports car to a pickup truck
required months of study and practice!
Similarly, railroads benefit from a standard for coupling rail cars together.
The result: trains can be assembled from combinations of all sorts
of different rail cars, built by different companies,
in different places, and with different purposes.

We believe that numerical models,
and the sub-components that make up those models,
should offer a similar kind of standardization.
To this end,
the `Community Surface Dynamics Modeling System`_ (CSDMS)
has developed the *Basic Model Interface* (BMI):
a set of standard control and query functions that,
when added to a model code,
make that model both easier to learn
and easier to couple with other software elements.

While a BMI can be written for any language,
CSDMS currently supports five languages:
C, C++, Fortran, Java, and Python.
The specification for each language is given in Table 1,
along with a corresponding example
in which the BMI is implemented.

.. _specs_and_examples:

.. table:: **Table 1:** BMI language specifications.
   :align: center
   :widths: 20, 25, 25, 30

   ========  =============  ==============  ======================
   Language  Specification  Repository      Example
   ========  =============  ==============  ======================
   C         `bmi.h`_       `bmi-c`_        `bmi-example-c`_
   C++       `bmi.hxx`_     `bmi-cxx`_      `bmi-example-cxx`_
   Fortran   `bmi.f90`_     `bmi-fortran`_  `bmi-example-fortran`_ 
   Java      `bmi.java`_    `bmi-java`_     `bmi-example-java`_ 
   Python    `bmi.py`_      `bmi-python`_   `bmi-example-python`_
   ========  =============  ==============  ======================

Along with the examples,
two documents may be particularly helpful when writing a BMI:

* :ref:`Getting Started Guide <getting_started>` --- a place to start
  if you haven't written a BMI before
* :ref:`BMI Best Practices <best_practices>` --- our collected wisdom on
  implementing a BMI

A complete description of the functions that make up the BMI is given next.


.. _basic_model_interface:

The Basic Model Interface
=========================

The functions that comprise the Basic Model Interface
can be grouped into categories:

* :ref:`control_funcs`
* :ref:`info_funcs`
* :ref:`var_funcs`
* :ref:`time_funcs`
* :ref:`getter_setter_funcs`
* :ref:`grid_funcs`

Table 2 lists the individual BMI functions
along with a brief description.
Following the table is a detailed description of each function,
including the function prototype in :term:`SIDL`,
grouped by functional category.

.. table:: **Table 2:** Summary of BMI functions.
   :align: center
   :widths: 30, 70

   ==============================  =========================================
   Function                        Description 
   ==============================  =========================================
   :ref:`initialize`               Perform startup tasks for the model.
   :ref:`update`                   Advance model state by one time step.
   :ref:`update_until`             Advance model state until the given time.
   :ref:`finalize`                 Perform tear-down tasks for the model.
   :ref:`get_component_name`       Name of the model.
   :ref:`get_input_item_count`     Count of a model's input variables.
   :ref:`get_output_item_count`    Count of a model's output variables.
   :ref:`get_input_var_names`      List of a model's input variables.
   :ref:`get_output_var_names`     List of a model's output variables.
   :ref:`get_var_grid`             Get the grid identifier for a variable.
   :ref:`get_var_type`             Get the data type of a variable.
   :ref:`get_var_units`            Get the units of a variable.
   :ref:`get_var_itemsize`         Get the size (in bytes) of one element of a variable.
   :ref:`get_var_nbytes`           Get the total size (in bytes) of a variable.
   :ref:`get_var_location`         Get the grid element type of a variable.
   :ref:`get_current_time`         Current time of the model.
   :ref:`get_start_time`           Start time of the model.
   :ref:`get_end_time`             End time of the model.
   :ref:`get_time_units`           Time units used in the model.
   :ref:`get_time_step`            Time step used in the model.
   :ref:`get_value`                Get a copy of values of a given variable.
   :ref:`get_value_ptr`            Get a reference to the values of a given variable.
   :ref:`get_value_at_indices`     Get variable values at specific locations.
   :ref:`set_value`                Set the values of a given variable.
   :ref:`set_value_at_indices`     Set the values of a variable at specific locations.
   :ref:`get_grid_rank`            Get the number of dimensions of a computational grid.
   :ref:`get_grid_size`            Get the total number of elements of a computational grid.
   :ref:`get_grid_type`            Get the grid type as a string.
   :ref:`get_grid_shape`           Get the dimensions of a computational grid.
   :ref:`get_grid_spacing`         Get the spacing between grid nodes.
   :ref:`get_grid_origin`          Get the origin of a grid.
   :ref:`get_grid_x`               Get the locations of a grid's nodes in dimension 1.
   :ref:`get_grid_y`               Get the locations of a grid's nodes in dimension 2.
   :ref:`get_grid_z`               Get the locations of a grid's nodes in dimension 3.
   :ref:`get_grid_node_count`      Get the number of nodes in the grid.
   :ref:`get_grid_edge_count`      Get the number of edges in the grid.
   :ref:`get_grid_face_count`      Get the number of faces in the grid.
   :ref:`get_grid_edge_nodes`      Get the edge-node connectivity.
   :ref:`get_grid_face_edges`      Get the face-edge connectivity.
   :ref:`get_grid_face_nodes`      Get the face-node connectivity.
   :ref:`get_grid_nodes_per_face`  Get the number of nodes for each face.
   ==============================  =========================================

.. include:: bmi.control_funcs.rst
.. include:: bmi.info_funcs.rst
.. include:: bmi.var_funcs.rst
.. include:: bmi.time_funcs.rst
.. include:: bmi.getter_setter.rst
.. include:: bmi.grid_funcs.rst


..
   Links

.. _Community Surface Dynamics Modeling System: https://csdms.colorado.edu
.. _bmi.h: https://github.com/csdms/bmi-c/blob/master/bmi.h
.. _bmi.hxx: https://github.com/csdms/bmi-cxx/blob/master/bmi.hxx
.. _bmi.f90: https://github.com/csdms/bmi-fortran/blob/master/bmi.f90
.. _bmi.java: https://github.com/csdms/bmi-java/blob/master/src/main/java/edu/colorado/csdms/bmi/BMI.java
.. _bmi.py: https://github.com/csdms/bmi-python/blob/master/bmipy/bmi.py
.. _bmi-c: https://github.com/csdms/bmi-c
.. _bmi-cxx: https://github.com/csdms/bmi-cxx
.. _bmi-fortran: https://github.com/csdms/bmi-fortran
.. _bmi-java: https://github.com/csdms/bmi-java
.. _bmi-python: https://github.com/csdms/bmi-python
.. _bmi-example-c: https://github.com/csdms/bmi-example-c
.. _bmi-example-cxx: https://github.com/csdms/bmi-example-cxx
.. _bmi-example-fortran: https://github.com/csdms/bmi-example-fortran
.. _bmi-example-java: https://github.com/csdms/bmi-example-java
.. _bmi-example-python: https://github.com/csdms/bmi-example-python
.. _UDUNITS: http://www.unidata.ucar.edu/software/udunits/
.. _The Units Database: https://www.unidata.ucar.edu/software/udunits/udunits-current/doc/udunits/udunits2.html#Database
.. _time unit conventions: https://www.unidata.ucar.edu/software/udunits/udunits-current/udunits/udunits2-accepted.xml
.. _primitive types: https://docs.oracle.com/javase/tutorial/java/nutsandbolts/datatypes.html
.. _wrapper classes: https://docs.oracle.com/javase/tutorial/java/data/numberclasses.html
.. include:: ../../CITATION.rst
.. _getting_started:

How to get started
==================

If you want to add a Basic Model Interface (BMI) to your own model,
here are some tips on getting started.

1. Take a look at the list of BMI
:ref:`control functions and their descriptions <basic_model_interface>`
to get an idea of what these functions are meant to do
and how they provide a standard set of controls for your model.

2. On the theory that it's often easier to start with an example and
modify it, we recommend starting with a copy of
the :ref:`BMI example code<specs_and_examples>`---look for the version
written in the language of your model---to use as a template.
Each link points to a GitHub repository that includes an
example model called "heat", and a corresponding BMI in a file called
"bmi_heat" (so, for example, the Python version is "bmi_heat.py",
the C version is "bmi_heat.c" and "bmi_heat.h", etc.).
You can grab copies by cloning the repository to your local machine
and doing a file copy from there,
or by selecting, copying, and pasting the code directly
into your favorite editor.

3. Modify each of the BMI functions in the template file so that it works
on your model, rather than on the original "heat" example. Depending
on how your model code is structured, you may need to do some
:term:`refactoring <refactor>`.
For example, if your model's initialization and main
processing are lumped together in the same body of code, you will need
to divide them into separate functions. Our experience is that codes
that are already modular usually need little or no modification,
whereas codes that are more monolithic tend to require a lot more
refactoring (but that is probably worthwhile anyway for the sake of
the quality and sustainability of the code!).
Each case is a bit different.
Be sure to check out our :ref:`BMI best practices <best_practices>` document
for tips.
We encourage you to contact us with questions by posting an
issue on the `CSDMS Help Desk`_.

4. Test it out. Try writing a program or script that initializes your
model with a simple test case using the :ref:`initialize` function,
runs it with the :ref:`update` or :ref:`update_until` functions,
and exchanges data using :ref:`get_value` and :ref:`set_value`.
Run a test to verify that you get the same output from your BMI'd model
that you got from it prior to BMI'ing.
(Note: we recommend writing a :term:`unit test` for each of your
BMI functions; to learn more about unit tests,
check out our `webinar`_).

5. *(Optional but cool)* With a few additional steps, you can make your
model operate as a :ref:`pymt <pymt>` component,
so you can drive it directly from a Python shell,
and write Python scripts to couple it with other models.
(Learn more about *pymt* :ref:`here <pymt>`, and through its `documentation`_.)
The CSDMS Integration Facility team can provide help and guidance on this
process: just contact us through the `CSDMS Help Desk`_.


.. Links:

.. _CSDMS Help Desk: https://github.com/csdms/help-desk
.. _webinar: https://csdms.colorado.edu/wiki/Presenters-0478
.. _documentation: https://pymt.readthedocs.io
.. _control_funcs:

Model control functions
-----------------------

These BMI functions are critical to plug-and-play modeling because
they allow a calling component to bypass a model's own time loop.
They also provide the caller with fine-grained control over the
model -- a calling application is able to, for instance, update a
model one time step at a time, change its state, and then continue
updating.

.. _initialize:

*initialize*
............

.. code-block:: java

    /* SIDL */
    int initialize(in string config_file);


The `initialize` function accepts a string argument that gives the
path to its :term:`configuration file`.
This function should perform all tasks that are to take place before
entering the model's time loop.  Models should be refactored, if
necessary, to read their inputs (which could include filenames for
other input files) from a configuration file.
BMI does not impose any constraint on how configuration files are
formatted.

**Implementation notes**

* Models should be refactored, if necessary, to use a configuration
  file.
* While no constraints are placed on how configuration files are
  formatted, `YAML <https://yaml.org>`_ is preferred.
* In C and Fortran, the *config_file* argument is passed as
  a character array, whereas in C++, Java, and Python, it's passed as
  a string -- a basic type in these languages.
* In C and Fortran, an integer status code indicating success (zero) or failure (nonzero)
  is returned. In C++, Java, and Python, an exception is raised on failure.

[:ref:`control_funcs` | :ref:`basic_model_interface`]


.. _update:

*update*
........

.. code-block:: java

    /* SIDL */
    int update();

The `update` function advances the model by a single time step. This
is the model's own internal time step (as returned by the BMI
:ref:`get_time_step` function), not the time step
of a controlling application.
This function should perform all tasks that take place during one
pass through the model's time loop.  It does not contain the time
loop. This typically involves incrementing all of the model's state
variables.  If the model's state variables don't change in time,
then they can be computed by the :ref:`initialize` function and this
function can just return without doing anything.

**Implementation notes**

* In C and Fortran, an integer status code indicating success (zero) or failure (nonzero)
  is returned. In C++, Java, and Python, an exception is raised on failure.

[:ref:`control_funcs` | :ref:`basic_model_interface`]


.. _update_until:

*update_until*
..............

.. code-block:: java

    /* SIDL */
    int update_until(in double time);

The `update_until` function updates the model to a particular time,
as provided by its *time* argument.
If the model permits,
the *time* argument can be a non-integral multiple of time steps,
and even negative.
Once called, the value returned
by the BMI :ref:`get_current_time` function must return the provided time
to reflect that the model was updated to the requested time.

**Implementation notes**

* Time is always a double-precision value.
* In C and Fortran, an integer status code indicating success (zero) or failure (nonzero)
  is returned. In C++, Java, and Python, an exception is raised on failure.

[:ref:`control_funcs` | :ref:`basic_model_interface`]


.. _finalize:

*finalize*
..........

.. code-block:: java

    /* SIDL */
    int finalize();


The `finalize` function should perform all tasks that take place
after exiting the model's time loop.  This typically includes
deallocating memory, closing files and printing reports.

**Implementation notes**

* In C and Fortran, an integer status code indicating success (zero) or failure (nonzero)
  is returned. In C++, Java, and Python, an exception is raised on failure.

[:ref:`control_funcs` | :ref:`basic_model_interface`]
.. _time_funcs:

Time functions
--------------

These simple diagnostic functions provide information on model time.
Model time is always expressed as a floating point value.

.. _get_current_time:

*get_current_time*
..................

.. code-block:: java

   /* SIDL */
   int get_current_time(out double time);

The current model time.

**Implementation notes**

* In C++, Java, and Python, the argument is omitted and the time is returned
  from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`time_funcs` | :ref:`basic_model_interface`]


.. _get_start_time:

*get_start_time*
................

.. code-block:: java

   /* SIDL */
   int get_start_time(out double time);

The start time of the  model.

**Implementation notes**

* The start time in BMI is typically defined to be 0.0.
* In C++, Java, and Python, the argument is omitted and the time is returned
  from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`time_funcs` | :ref:`basic_model_interface`]


.. _get_end_time:

*get_end_time*
................

.. code-block:: java

   /* SIDL */
   int get_end_time(out double time);

The end time of the  model.

**Implementation notes**

* If the model doesn't define an end time, a large number (e.g., the
  largest floating point number supported on a platform) is typically
  chosen.
* In C++, Java, and Python, the argument is omitted and the time is returned
  from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`time_funcs` | :ref:`basic_model_interface`]


.. _get_time_units:

*get_time_units*
................

.. code-block:: java

   /* SIDL */
   int get_time_units(out string units);

Get the units of time as reported by the model's BMI (through
:ref:`get_current_time`, :ref:`get_end_time`, etc.).
It's recommended to use `time unit conventions`_ from Unidata's
`UDUNITS`_ package; e.g., ``"s"``, ``"min"``, ``"h"``, ``"d"``.

**Implementation notes**

* Avoid using ``"years"`` as a unit, if possible, since a year is
  difficult to define precisely. UDUNITS defines a year as 365.2422
  days or 31556926 seconds.
* Dimensionless quantities should use ``""`` or ``"1"`` as the unit.
* Models that don't vary with time, or don't have time units should
  use ``"none"``.
* In C++, Java, and Python, the argument is omitted and the units are returned
  from the function as a string.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`time_funcs` | :ref:`basic_model_interface`]


.. _get_time_step:

*get_time_step*
...............

.. code-block:: java

   /* SIDL */
   int get_time_step(out double time_step);

Get the time step used in the model.
The time step is always expressed as a floating point value.

**Implementation notes**

* A time step is typically a positive value. However, if the model
  permits it, a negative value can be used (running the model
  backward).
* In C++, Java, and Python, the argument is omitted and the time step is returned
  from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`time_funcs` | :ref:`basic_model_interface`]
.. image:: _static/bmi-logo-header-text.png
    :align: center
    :scale: 85%
    :alt: Basic Model Interface (BMI)
    :target: https://bmi.readthedocs.io/

.. title:: BMI

.. include:: bmi.spec.rst


Additional Topics
=================

.. toctree::
   :maxdepth: 2

   bmi.getting_started
   bmi.best_practices
   model_grids
   csdms
   glossary
   references

BMI is an element of the `CSDMS Workbench`_,
an integrated system of software tools, technologies, and standards
for building and coupling models.


Help
----

Adding a BMI to a model can be a daunting task.
If you'd like assistance,
CSDMS can help.
Depending on your need, we can provide advice or consulting services.
Feel free to contact us through the `CSDMS Help Desk`_.


Project Information
-------------------

.. toctree::
   :maxdepth: 1

   citation
   contributing
   conduct
   credits
   Governance <governance>
   council
   partners


.. Links:

.. _CSDMS Workbench: https://csdms.colorado.edu/wiki/Workbench
.. _CSDMS Help Desk: https://github.com/csdms/help-desk
.. _info_funcs:

Model information functions
---------------------------

These functions supply the model name
and the model's :term:`exchange items <exchange item>` -- 
the variables that the model can use from
and provide to other models that have a BMI.

.. _get_component_name:

*get_component_name*
....................

.. code-block:: java

  /* SIDL */
  int get_component_name(out string name);

This function supplies the name of the model component as a string.
There are no restrictions on the name,
but it should be unique to prevent conflicts with other components.

**Implementation notes**

* In C and Fortran, the *name* argument is a a character array, and an integer
  status code indicating success (zero) or failure (nonzero) is returned.
* In C++, Java, and Python, this argument is omitted, and a string -- a basic type
  in these languages -- is returned from the function.

[:ref:`info_funcs` | :ref:`basic_model_interface`]


.. _get_input_item_count:

*get_input_item_count*
......................

.. code-block:: java

  /* SIDL */
  int get_input_item_count(out int count);

The number of variables the model can use from other models
implementing a BMI.
Also the number of variables that can be set with :ref:`set_value`.

**Implementation notes**

* In C++, Java, and Python, the argument is omitted and the count is returned
  from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`info_funcs` | :ref:`basic_model_interface`]


.. _get_output_item_count:

*get_output_item_count*
.......................

.. code-block:: java

  /* SIDL */
  int get_output_item_count(out int count);

The number of variables the model can provide other models
implementing a BMI.
Also the number of variables that can be retrieved with :ref:`get_value`.

**Implementation notes**

* In C++, Java, and Python, the argument is omitted and the count is
  returned from the function.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`info_funcs` | :ref:`basic_model_interface`]


.. _get_input_var_names:

*get_input_var_names*
.....................

.. code-block:: java

  /* SIDL */
  int get_input_var_names(out array<string, 1> names);

Gets an array of names for the variables the model can use from other
models implementing a BMI.
The length of the array is given by :ref:`get_input_item_count`.
The names are preferably in the form of CSDMS :term:`Standard Names`.
Standard Names enable a modeling framework to determine whether an
input variable in one model is equivalent to, or compatible with,
an output variable in another model.
This allows the framework to automatically connect components.
Standard Names do not have to be used within the model.

**Implementation notes**

* In C and Fortran, the names are passed back as an array of character
  pointers (because the variable names could have differing lengths), and an
  integer status code indicating success (zero) or failure (nonzero) is returned.
* In C++, the argument is omitted and the names are returned from the
  function in a vector, a standard container in the language.
* In Java, the argument is omitted and the names are returned from the
  function in a string array, a standard container in the language.
* In Python, the argument is omitted and the names are returned from the
  function in a tuple, a standard container in the language.
* A model might have no input variables.

[:ref:`info_funcs` | :ref:`basic_model_interface`]


.. _get_output_var_names:

*get_output_var_names*
......................

.. code-block:: java

  /* SIDL */
  int get_output_var_names(out array<string, 1> names);

Gets an array of names for the variables the model can provide to other
models implementing a BMI.
The length of the array is given by :ref:`get_output_item_count`.
The names are preferably in the form of CSDMS :term:`Standard Names`.
Standard Names enable a modeling framework to determine whether an
input variable in one model is equivalent to, or compatible with,
an output variable in another model.
This allows the framework to automatically connect components.
Standard Names do not have to be used within the model.

**Implementation notes**

* In C and Fortran, the names are passed back as an array of character
  pointers (because the variable names could have differing lengths), and an
  integer status code indicating success (zero) or failure (nonzero) is returned.
* In C++, the argument is omitted and the names are returned from the
  function in a vector, a standard container in the language.
* In Java, the argument is omitted and the names are returned from the
  function in a string array, a standard container in the language.
* In Python, the argument is omitted and the names are returned from the
  function in a tuple, a standard container in the language.
* A model may have no output variables.

[:ref:`info_funcs` | :ref:`basic_model_interface`]
.. _references:

References
==========

Community Surface Dynamics Modeling System (CSDMS)

  Syvitski, James P; Hutton, Eric; Piper, Mark; Overeem, Irina;
    Kettner, Albert; Peckham, Scott (2014): "Plug and play component
    modeling--the CSDMS 2.0 approach", *International Environmental
    Modelling and Software Society (iEMSs), 7th International Congress
    on Environmental Modelling and Software* (Ames DP, Quinn N., eds),
    San Diego, CA,
    USA. http://www.iemss.org/society/index.php/iemss-2014-proceedings.

Basic Model Interface (BMI)

  Peckham, Scott D; Hutton, Eric WH; Norris, Boyana (2013): "A
    component-based approach to integrated modeling in the
    geosciences: The design of CSDMS", *Computers & Geosciences*,
    **53**: 3-12, ISSN 0098-3004,
    http://dx.doi.org/10.1016/j.cageo.2012.04.002.

  Hutton, Eric WH; Piper, Mark D; Tucker, Gregory E (2020): "The Basic
    Model Interface 2.0: A standard interface for coupling numerical
    models in the geosciences", *Journal of Open Source Software*,
    **5(51)**, 2317, https://doi.org/10.21105/joss.02317.

Scientific Interface Definition Language (SIDL)

  Dahlgren, Tamara; Ebner, Dietmar; Epperly, Thomas; Kumfert, Gary;
    Leek, James; Prantl, Adrian (2012): "Babel User's Guide: Part I
    Foundation: SIDL Basics".
    https://computing.llnl.gov/projects/babel-high-performance-language-interoperability/docs/users_guide/index008.html (retrieved 2019-10-09).

  Epperly, Thomas GW; Kumfert, Gary; Dahlgren, Tamara; Ebner, Dietmar;
    Leek, Jim; Prantl, Adrian; Kohn, Scott (2011): "High-performance
    language interoperability for scientific computing through
    Babel". *The International Journal of High Performance Computing
    Applications*. **26** (3): 260–274.
    http://dx.doi.org/10.1177/1094342011414036

Model grids

  Hoffmann, Klaus A; Chiang, Steve T (2000): Computational Fluid
    Dynamics Volume I. Third edition. *Engineering Education System*.

  Slingerland, Rudy; Kump, Lee (2011): Mathematical Modeling of
    Earth's Dynamical Systems: A Primer. *Princeton University Press*.
.. _getter_setter_funcs:

Variable getter and setter functions
------------------------------------

These functions are used to access and modify the variables
that a model exposes through its BMI
(see :ref:`get_input_var_names` and :ref:`get_output_var_names`).

A *getter* is a function called to get a variable from a model's *state*.
A model's state variables typically change with each time step,
so getters are called to get current values.

A *setter* is a function called to change/overwrite a variable in
a model's state. A setter may impose restrictions on how a
state variable can be changed or check the new data for validity. 


.. _get_value:

*get_value*
...........

.. code-block:: java

   /* SIDL */
   int get_value(in string name, in array<> dest);

The `get_value` function takes a variable name and copies values into a
provided array parameter.
The type and size of the array parameter depend on the variable,
and can be determined through
:ref:`get_var_type`, :ref:`get_var_nbytes`, etc.
Recall that arrays are always flattened in BMI,
even if the model uses dimensional variables.

**Implementation notes**

* The *dest* argument must be defined and allocated before calling
  `get_value`. Whatever values it contains are overwritten in the call
  to `get_value`.
* In Python, the array parameter is a :term:`numpy` array.
* In Java, only `primitive types`_ (e.g., ``int``, ``double``), not
  `wrapper classes`_ (e.g., ``Integer``, ``Double``), are supported.
* In C++ and Java, `get_value` is a void function.
* Depending on how a model is written, a variable may not be
  accessible until after the call to :ref:`initialize`. Likewise, the
  variable may not be accessible after calling :ref:`finalize`.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`getter_setter_funcs` | :ref:`basic_model_interface`]

.. _get_value_ptr:

*get_value_ptr*
...............

.. code-block:: java

   /* SIDL */
   int get_value_ptr(in string name, out array<> dest_ptr);

The `get_value_ptr` function takes a variable name and returns a reference
to a variable.
Unlike the array parameter returned from :ref:`get_value`,
the reference always points to the current values of the variable,
even if the model's state has changed.

**Implementation notes**

* The reference points to a flattened array.
* In C++ and Java, the *dest_ptr* argument is omitted, and the reference is
  returned through the function.
* In Java, this function can only be used with one-dimensional arrays.
* In Python, a :term:`numpy` array is returned.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`getter_setter_funcs` | :ref:`basic_model_interface`]


.. _get_value_at_indices:

*get_value_at_indices*
......................

.. code-block:: java

   /* SIDL */
   int get_value_at_indices(in string name, in array<> dest, in array<int, 1> inds);

Use the `get_value_at_indices` function to get a copy of a variable's values
at the locations specified by the one-dimensional array indices
in the *inds* argument.
The values are returned through the *dest* argument.

**Implementation notes**

All the notes from :ref:`get_value` apply.
Additionally,

* Both *dest* and *inds* are flattened arrays.
* The *inds* argument is always of type integer.

[:ref:`getter_setter_funcs` | :ref:`basic_model_interface`]


.. _set_value:

*set_value*
...........

.. code-block:: java

   /* SIDL */
   int set_value(in string name, in array<> src);

The `set_value` function takes a variable name and an array of values,
*src*,
and copies those values into the model's internal array of values,
overwriting the current contents.
The type and size of *src* must match the model's internal array,
and can be determined through
:ref:`get_var_type`, :ref:`get_var_nbytes`, etc.
Recall that arrays are always flattened in BMI,
even if the model uses dimensional variables.

**Implementation notes**

* In Python, *src* is a :term:`numpy` array.
* In Java, only `primitive types`_ (e.g., ``int``, ``double``), not
  `wrapper classes`_ (e.g., ``Integer``, ``Double``), are supported.
* In C++ and Java, `set_value` is a void function.
* Depending on how a model is written, a variable may not be
  accessible until after the call to :ref:`initialize`. Likewise, the
  variable may not be accessible after calling :ref:`finalize`.
* In C and Fortran, an integer status code indicating success (zero) or failure
  (nonzero) is returned.

[:ref:`getter_setter_funcs` | :ref:`basic_model_interface`]


.. _set_value_at_indices:

*set_value_at_indices*
......................

.. code-block:: java

   /* SIDL */
   int set_value_at_indices(in string name, in array<int, 1> inds, in array<> src);

Use the `set_value_at_indices` function to set a variable's values
at the locations specified by the one-dimensional array indices
in the *inds* argument.

**Implementation notes**

All the notes from :ref:`set_value` apply.
Additionally,

* Both *src* and *inds* are flattened arrays.
* The *inds* argument is always of type integer.

[:ref:`getter_setter_funcs` | :ref:`basic_model_interface`]
