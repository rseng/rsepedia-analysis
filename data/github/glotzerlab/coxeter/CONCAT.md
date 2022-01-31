# coxeter Contributor Agreement

These terms apply to your contribution to the coxeter Open Source Project ("Project") owned or managed by the Regents of the University of Michigan ("Michigan"), and set out the intellectual property rights you grant to Michigan in the contributed materials. If this contribution is on behalf of a company, the term "you" will also mean the company you identify below. If you agree to be bound by these terms, fill in the information requested below and provide your signature.

1. The term "contribution" means any source code, object code, patch, tool, sample, graphic, specification, manual, documentation, or any other material posted or submitted by you to a project.
2. With respect to any worldwide copyrights, or copyright applications and registrations, in your contribution:
    * you hereby assign to Michigan joint ownership, and to the extent that such assignment is or becomes invalid, ineffective or unenforceable, you hereby grant to Michigan a perpetual, irrevocable, non-exclusive, worldwide, no-charge, royalty-free, unrestricted license to exercise all rights under those copyrights. This includes, at Michigan's option, the right to sublicense these same rights to third parties through multiple levels of sublicensees or other licensing arrangements;
    * you agree that both Michigan and you can do all things in relation to your contribution as if each of us were the sole owners, and if one of us makes a derivative work of your contribution, the one who makes the derivative work (or has it made) will be the sole owner of that derivative work;
    * you agree that you will not assert any moral rights in your contribution against us, our licensees or transferees;
    * you agree that we may register a copyright in your contribution and exercise all ownership rights associated with it; and
    * you agree that neither of us has any duty to consult with, obtain the consent of, pay or render an accounting to the other for any use or distribution of your contribution.
3. With respect to any patents you own, or that you can license without payment to any third party, you hereby grant to Michigan a perpetual, irrevocable, non-exclusive, worldwide, no-charge, royalty-free license to:
    * make, have made, use, sell, offer to sell, import, and otherwise transfer your contribution in whole or in part, alone or in combination with or included in any product, work or materials arising out of the project to which your contribution was submitted; and
    * at Michigan's option, to sublicense these same rights to third parties through multiple levels of sublicensees or other licensing arrangements.
4. Except as set out above, you keep all right, title, and interest in your contribution. The rights that you grant to Michigan under these terms are effective on the date you first submitted a contribution to Michigan, even if your submission took place before the date you sign these terms. Any contribution Michigan makes available under any license will also be made available under a suitable Free Software Foundation or Open Source Initiative approved license.
5. With respect to your contribution, you represent that:
    * it is an original work and that you can legally grant the rights set out in these terms;
    * it does not to the best of your knowledge violate any third party's copyrights, trademarks, patents, or other intellectual property rights; and
you are authorized to sign this contract on behalf of your company (if identified below).
6. The terms will be governed by the laws of the State of Michigan and applicable U.S. Federal Law. Any choice of law rules will not apply.

**By making contribution, you electronically sign and agree to the terms of the coxeter Contributor Agreement.**

![by-sa.png](https://licensebuttons.net/l/by-sa/3.0/88x31.png)

Based on the Sun Contributor Agreement - version 1.5.
This document is licensed under a Creative Commons Attribution-Share Alike 3.0 Unported License
http://creativecommons.org/licenses/by-sa/3.0/
# How to contribute to the project

## Feedback

Issue reports and feature proposals are very welcome.
Please use the [GitHub issue page](https://github.com/glotzerlab/coxeter/issues/) for this.

## Contributing code

Code contributions to the coxeter open-source project are welcomed via pull requests on [GitHub](https://github.com/glotzerlab/coxeter/).
Prior any work you should contact the developers to ensure that the planned development meshes well with the directions and standards of the project.
All contributors must agree to the Contributor Agreement ([ContributorAgreement.md](ContributorAgreement.md)) before their pull request can be merged.
---
title: 'coxeter: A Python package for working with shapes'
tags:
  - Python
  - geometry
  - physics
  - materials science
authors:
  - name: Vyas Ramasubramani
    orcid: 0000-0001-5181-9532
    affiliation: 1
  - name: Bradley D. Dice
    orcid: 0000-0002-9983-0770
    affiliation: 2
  - name: Tobias T. Dwyer
    orcid: 0000-0001-6443-7744
    affiliation: 1
  - name: Sharon C. Glotzer
    orcid: 0000-0002-7197-0085
    affiliation: "1, 2, 3"
affiliations:
 - name: Department of Chemical Engineering, University of Michigan
   index: 1
 - name: Department of Physics, University of Michigan
   index: 2
 - name: Biointerfaces Institute, University of Michigan
   index: 3
date: 25 Feb 2021
bibliography: paper.bib
---

# Package Overview

![The coxeter package supports calculating a wide range of properties on shapes in two and three dimensions. The central bubble in the figure contains a subset the shapes supported by coxeter, which includes simple polygons in 2D and arbitrary 3D polyhedral meshes. The surrounding bubbles represent a sampling of what coxeter can do with a shape. From the top-left going clockwise, these are: determining the inspheres of polytopes (as well as minimal bounding spheres for all shapes); calculating the distance from the centroid of a shape to its boundary; checking whether a point lies inside or outside a shape; determining the circumspheres of polytopes (as well as maximal bounded spheres for all shapes); calculating moments of inertia in 2D and inertia tensors in 3D (including support for reorienting a shape along its principal axes); and anisotropic form factors for arbitrary shapes, which are important for scattering calculations.](figure1/Figure1.pdf){ width=100% }

The coxeter Python package provides tools to represent, generate, and compute properties of shapes in two and three dimensions.
The package emphasizes simplicity and flexibility, using a common set of abstractions to present a largely uniform interface across various shapes and allowing easy mutation of almost all of their geometric attributes.
The package also serves as a repository for specific groups of shapes, exposing an easy-to-use API for shape generation that users can extend to make particular geometric objects collectively accessible.


# Statement of Need

Considerations of shape are becoming increasingly important in materials science as improved synthetic capabilities have allowed the creation of a wide range of anisotropic particles [@Glotzer2007b].
Colloidal science in particular has seen immense growth in this area, and numerous studies have shown that particle shape is an important handle for controlling the self-assembly [@Damasceno2012d] and packing [@Chen2014] of colloidal crystals.
Precise modeling of these systems requires reproducible methods for generating shapes and calculating their properties [@Anderson2020; @Allen2006].
An important aspect of achieving this reproducibility is making canonical definitions of shapes used in particular studies readily available to other researchers.
Furthermore, since these shapes may be used in physics-based simulations, any calculations must be robust enough to handle any numerical issues that may arise across a wide range of different geometries.
Some of the applications of coxeter to date include: the development of equations of state for polyhedral particles [@Irrgang2017]; the calculation of physical properties for dynamical simulation of anisotropic particles [@Ramasubramani2020b]; and the orientational ordering of ellipsoidal colloids in a magnetic field [@Kao2019].


# Summary of Features

## Computing Geometric and Physical Properties

The central elements in coxeter are the shape classes, which encode the features of particular types of shapes.
In order to enforce a highly uniform API and ensure conceptual clarity, all shape classes inherit from a small set of abstract base classes that encode specific subsets of properties: for instance, the standard properties of all shapes in two dimensions are embedded in the ``Shape2D`` class.
In addition to standard types of shapes such as ellipsoids or polygons, coxeter also includes more esoteric shape classes like spheropolyhedra, which are important for modeling the imperfect rounded polyhedra frequently synthesized at the nano- and micron scales [@Zhang2011; @Rossi2015].
Even simple properties like surface areas are generally nontrivial to compute for such shapes, but using coxeter spheropolyhedra are no more difficult to work with than any other 3D shape.
Working with convex polygons and polyhedra using coxeter is greatly simplified via internal use of SciPy's convex hull calculations [@Virtanen2020], allowing the user to simply provide a set of vertices while coxeter performs the appropriate determination of facets and plane equations based on the simplices of the convex hull.

The shape classes transparently expose many geometric attributes in the form of settable Python properties, allowing on-the-fly rescaling or transformation of the shape.
This aspect of coxeter is designed to fill a common gap in most computational geometry libraries, which typically focus on solving more complex problems like finding convex hulls, Voronoi tessellations, and Delaunay triangulations [@cgal]; coxeter aims to provide a standard implementation for simpler calculations such as surface areas and bounding spheres for which formulas are generally well-known but often require careful consideration to calculate robustly and efficiently.
These properties range from standard calculations like volumes and surface areas to less common metrics like mean curvatures and asphericities that are relevant for computing equations of state for polyhedral particles [@Irrgang2017].
The package also provides various types of bounding and bounded spheres of shapes, which are measures of the extent of polygons and polyhedra within crystal structures.
To simplify interoperability with other packages in the scientific computing ecosystem, non-scalar properties are generally provided as NumPy arrays [@Harris2020].

In addition to purely geometric properties, shapes in coxeter also expose various physically relevant quantities in order to support a wide range of applications for shapes of constant density.
Some examples of such properties are inertia tensors, which are integral to the equations of motion for anisotropic bodies, and scattering form factors, which are Fourier transforms of the shape volume that help characterize structure in condensed matter physics [@AlsNielsen2011].
Since physical equations and observables can be highly sensitive to inputs like inertia tensors, coxeter emphasizes robust methods for their evaluation [@Kallay2006].
Two dimensional shapes like polygons are embedded in three dimensions rather than in the plane, so coxeter uses the rowan library [@Ramasubramani2018] to rotate them into the plane and then compute various properties to avoid complications and numerical instabilities that arise from performing integrals over planar lamina embedded in 3D Euclidean space.

## Shape Generation

The library also serves as a repository for the generation of shapes.
While simple classes of shapes like spheres and ellipsoids can be described via a small fixed set of parameters, the definitions of polygons and polyhedra can be arbitrarily long depending on the number of vertices of these shapes.
The shape family API in coxeter provides a flexible way to define and work with collections of related shapes, ranging from enumerable sets like the Platonic solids to continuously defined sets of shapes [@Chen2014].
These different types of shape families are handled using identical APIs, so users can easily switch between shapes that have completely different mathematical definitions using a single line of code.
Shape families generate coxeter shape classes from input parameters, simplifying access to computed geometric and physical properties.

A number of such families are bundled into coxeter, but just as importantly, the framework allows users to work with arbitrary lists of shapes provided as dictionaries of attributes.
This dictionary-based definition can be simply converted to JSON, making it trivial to share representations of shapes.
The library also stores mappings from digital object identifiers (DOIs) to families, so that any user can contribute families associated with published research to make them collectively accessible.
We anticipate that the set of shape families in coxeter will grow over time as users generate and contribute their shape families to coxeter, with the goal of providing a centralized repository for use in reproducing and extending prior research, particularly in the field of shape-driven nanoparticle self-assembly.
Currently coxeter primarily supports the schema proposed by the GSD library [@GlotzerLabGSD], making it directly compatible with the HOOMD-blue molecular simulation tool [@Anderson2020], but other schema can be implemented as needed.

# Acknowledgements

This research was supported in part by the National Science Foundation, Division of Materials Research Award No. DMR-1808342.
V. R. also acknowledges the 2019-2020 J. Robert Beyster Computational Innovation Graduate Fellowship from the College of Engineering, University of Michigan.
B. D. acknowledges fellowship support from the National Science Foundation under ACI-1547580, S212: Impl: The Molecular Sciences Software Institute [@Wilkins-Diehr2018; @Krylov2018] and an earlier National Science Foundation Graduate Research Fellowship Grant DGE-1256260 (2016–2019).
T. D. is supported by a National Science Foundation Graduate Research Fellowship Grant DGE-1256260.

We would like to acknowledge M. Eric Irrgang for prototype implementations of various parts of this code, as well as Bryan VanSaders and James Proctor for collecting the various early prototypes and relevant methods into a single code base.
We additionally thank all code contributors.
In addition to the authors and aforementioned contributors, the list of code contributors includes Brandon Butler, Thomas Waltmann, Timothy Moore, Corwin Kerr, Eric Harper, Jens Glaser, William Zygmunt, and Mariano Semelman.

Finally, we would like to acknowledge the following authors of external open-source tools that are used in coxeter:

- Mark Dickinson, who wrote the **polyhedron** module used for point-in-polygon and point-in-polyhedron checks in coxeter.
- David Björkevik, who wrote the **polytri** package used for polygon triangulation in coxeter.
- Campbell Barton, who wrote the **isect_segments-bentley_ottman** package used to validate that polygons have no crossings (i.e. that polygons are simple) in coxeter.

# References
<!-- Provide a general summary of your changes in the Title above -->

## Description
<!-- Describe your changes in detail -->

## Motivation and Context
<!-- Why is this change required? What problem does it solve? -->
<!-- If it fixes an open issue, please link to the issue here. -->

## Types of Changes
<!-- Please select all items that apply either now or after creating the pull request: -->
- [ ] Documentation update
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change<sup>1</sup>

<sup>1</sup>The change breaks (or has the potential to break) existing functionality.

## Checklist:
<!-- Please select all items that apply either now or after creating the pull request. -->
<!-- If you are unsure about any of these items, do not hesitate to ask! -->
- [ ] I am familiar with the [**Contributing Guidelines**](https://github.com/glotzerlab/coxeter/blob/master/CONTRIBUTING.md).
- [ ] I agree with the terms of the [**Contributor Agreement**](https://github.com/glotzerlab/coxeter/blob/master/ContributorAgreement.md).
- [ ] The changes introduced by this pull request are covered by existing or newly introduced tests.
- [ ] I have updated the [changelog](https://github.com/glotzerlab/coxeter/blob/master/ChangeLog.txt) and the [credits](https://github.com/glotzerlab/coxeter/blob/master/doc/source/credits.rst).
---
name: Bug report
about: Create a report that describes an issue
title: ''
labels: ''
assignees: ''

---

<!-- Please replace the text in the individual sections below. -->

### Description

A concise description of what is causing an issue and what you expected to happen.

### To reproduce

Please copy & paste the code that is triggering the issue into this section and describe the context and prior steps that lead to the problem with as much detail as possible.

### Error output

If possible, copy any terminal outputs or attach screenshots that provide additional information on the problem.

### System configuration

Please complete the following information:

 - Operating System [e.g. macOS]:
 - Version of Python [e.g. 3.7]:
 - Version of coxeter [e.g. 1.0]:

Or copy & paste the output of: `python -c 'import platform; print(platform.platform()); import sys; print(sys.version); import coxeter; print(coxeter.__version__)'`
---
name: Feature request
about: Suggest a new feature for coxeter
title: ''
labels: ''
assignees: ''

---
<!-- Please replace the text in the individual sections below. -->

### Feature description

A concise description of the enhancement that this feature provides or the problem it solves.
For example: "I would like to ..." or "It would be nice if it was easier to ...".

### Proposed solution

A description of how you would like a possible solution to look like.

*Proposals for an implementation of the idea are welcome, but not necessary.*

### Additional context

Any alternative solutions you might have considered or other information that might be helpful to understand this feature request.
coxeter Developers
------------------

The following people contributed to the development of coxeter.

Vyas Ramasubramani - **Creator and lead developer**

* Created documentation pages.
* Formalized contribution guidelines and contributor agreement.
* Cleaned up damasceno module and separated out shape information into a JSON file that is read on demand.
* Fixed code formatting to conform to PEP8 requirements.
* Implemented Polygon class.
* Implemented ConvexSpheropolygon class.
* Implemented Polyhedron class.
* Implemented ConvexPolyhedron class.
* Implemented ConvexSpheropolyhedron class.
* Add ability to check if points are contained in convex polyhedra.
* Fix calculation of circumsphere to work for non-regular polyhedra.
* Fix calculation of circumcircle to work for non-regular polygons.
* Add ability to calculate minimum bounding sphere/circle for polyhedra/polygons.
* Implemented ConvexPolygon class.
* Added ReadTheDocs support.
* Added circumsphere from center calculation for convex polyhedra.
* Added shape getter for damasceno shapes.
* Define proper inertia tensor calculations and transformations for polygons and polyhedra.
* Added interoperability with the GSD shape specification.
* Developed shape families and all associated shape repository APIs.
* Add ability to diagonalize the inertia tensors of shapes.
* Defined base classes for all shapes.
* Standardize usage of Sphere/Circle classes for circum, in, and bounding sphere/circle calculations.
* Moved form factor amplitude calculations from legacy ft module to shape classes, cleaned and added more tests.
* Added point-in-shape checks for circles and ellipses.
* Added generic inertia tensors for 2D shapes.
* Added minimal bounding sphere for all shapes.
* Added minimal centered bounding sphere calculations for all shapes except general polygons and polyhedra.
* Enabled getting and setting the circumsphere or bounding sphere/circle radius of a polyhedron/polygon (for both types of bounding sphere/circle).
* Added maximal bounding sphere for all shapes.
* Added maximal centered bounded sphere calculations for all shapes except general polygons and polyhedra.
* Enabled getting and setting the insphere or bounded sphere/circle radius of a polyhedron/polygon (for both types of bounding sphere/circle).
* Added point in polygon checks.
* Added point in polyhedron checks.
* Added repr for all shapes.
* Fixed centroid calculations for polygon and polyhedron to use integrals rather than simple averages of vertices.
* Wrote example notebooks.

Bryan VanSaders - **Original maintainer of legacy euclid package**

* Created package layout.
* Original port of classes and methods into package.
* Added some methods to the utils module.
* Added symmetry groups.

James Proctor

* Ported some damasceno code into coxeter from standalone module.

Bradley Dice

* Migrated ft code into coxeter from freud and added tests.
* Added CircleCI support.
* Add ability to check if points are contained in convex spheropolyhedra.
* Revised and edited all documentation.
* Updated doctests to be part of pytest suite.
* Added automatic axis creation for plotting.
* Added spheropolygon area and perimeter setters.
* Added ellipse area setter and ellipsoid volume setter.
* Added plato support.

Brandon Butler

* Removed old quat\_tools module and modified modules to use rowan.
* Moved logic in FreudShape module to top-level package namespace.
* Moved all common shape definitions into a common\_shapes module.

Eric Harper

* Migrated shape classes into coxeter from freud.

Jens Glaser

* Bug fix for convex hull finding.

M. Eric Irrgang

* Bugfixes to imports.
* Implemented core shape classes.
* Implemented the ft module.

Carl Simon Adorf

* Implemented the damasceno module.

Matthew Spellings

* Added some methods to the utils module.
* Triangulation of core shape classes.

William Zygmunt

* Helped clean up utils module.

Tobias Dwyer

* Added getter and setter tests to some of the shape classes.
* Added examples for the shape classes.


Source code
-----------

**coxeter** includes the source code of the following Python packages and
modules.

.. highlight:: none

The source of polytri (https://github.com/bjorkegeek/polytri) is included
directly into the **coxeter** package. The module implementing that code is
reproduced in its entirety along with an additional ``__init__`` file to enable
its import as a subpackage. It is used for the triangulation of polygons and
the surface triangulation of polyhedra. This software is made available under
the MIT license::

    The MIT License (MIT)

    Copyright (c) 2016 David Björkevik

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE

The source of isect_segments-bentley_ottmann
(https://github.com/ideasman42/isect_segments-bentley_ottmann) is included
directly into the **coxeter** package. The module implementing that code is
reproduced in its entirety along with an additional ``__init__`` file to enable
its import as a subpackage. It is used to check whether a set of vertices
defines a simple or a complex polygon. This software is made available under
the MIT license::

    Copyright (c) 2010 by Bart Kiers
    Copyright (c) 2015 by Campbell Barton

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use,
    copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following
    conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

The source of polyhedron (https://github.com/mdickinson/polyhedron) is included
directly into the **coxeter** package. It is used for point in polygon/polyhedron
checks for general polygons and polyhedra (specifically, to calculate the winding
number). This software is made available under the BSD-3 license::

    BSD 3-Clause License

    Copyright (c) 2019, Mark Dickinson
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its
       contributors may be used to endorse or promote products derived from
       this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
The format is based on `Keep a Changelog <http://keepachangelog.com/en/1.0.0/>`__.
This project adheres to `Semantic Versioning <http://semver.org/spec/v2.0.0.html>`__.

v0.6.1 - 2021-07-15
-------------------

Fixed
~~~~~

- Typos in JOSS paper.

v0.6.0 - 2021-07-14
-------------------

Added
~~~~~

- Plotting and other graphical rendering of shapes using `plato <https://plato-draw.readthedocs.io/>`__.
- Notebooks with example use-cases for the package.
- A quickstart tutorial.

v0.5.0 - 2021-02-23
-------------------

Added
~~~~~

- Ellipse area setter and Ellipsoid volume setter.
- Point in circle checks.
- Point in ellipse checks.
- Inertia tensors for 2D shapes that implement moments of inertia.
- Add minimal bounding sphere for all shapes.
- Add minimal centered bounding sphere calculations for all shapes except general polygons, general polyhedra, spheropolygons, and spheropolyhedra.
- Enable getting and setting the circumsphere or bounding sphere radius of a polyhedron (for both types of bounding sphere).
- Add maximal bounded sphere for all shapes.
- Add maximal centered bounded sphere calculations for all shapes except general polygons, general polyhedra, spheropolygons, and spheropolyhedra.
- Enable getting and setting the insphere or bounded sphere radius of a polyhedron (for both types of bounding sphere).
- Point in polygon checks for general (nonconvex) polygons.
- Point in polyhedron checks for general (nonconvex) polyhedrons.
- Minimal bounding sphere for all shapes except spheropolygons and spheropolyhedra.
- Add minimal centered bounding sphere calculations for all shapes except general polygons, general polyhedra, spheropolygons, and spheropolyhedra.
- Getters and setters for the circumsphere or bounding sphere radius of a polyhedron (for both types of bounding sphere).
- A ``repr`` for all shapes.

Changed
~~~~~~~

- Ensure that hypothesis-based tests don't implicitly reuse pytest fixtures.

Deprecated
~~~~~~~~~~

- The circumsphere from center calculations (replaced by minimal centered bounding sphere).
- The bounding_sphere property is deprecated in favor of minimal_bounding_sphere.
- The insphere from center calculations (replaced by maximal centered bounded sphere).

Fixed
~~~~~

- Centroid calculations for polygon and polyhedron use the full integrals rather than simple averages of vertices.

v0.4.0 - 2020-10-14
-------------------

Added
~~~~~

- Circumsphere and insphere from center calculations for ConvexSpheropolyhedron.
- Form factors amplitudes for sphere, polygons, and polyhedra.
- Shape families associated with a DOI can be directly accessed via a dictionary.
- Expected abstract interface for shapes (both 2D and 3D) has expanded.
- Plotting polygons or polyhedra can automatically create matplotlib axes.
- Perimeter calculation for polygons.
- Area and perimeter setters for spheropolygons.

Changed
~~~~~~~

- Shape family API is now entirely based on class methods rather than a call operator.
- The parent ShapeFamily class is now part of the public API.
- Doctests are now run as part of pytest.
- Subpackages have been renamed: shape_classes is now shapes, and shape_families is now families.
- The common_families submodule of shape_families is now just common.

Fixed
~~~~~

- Volume calculation for ConvexSpheropolyhedron includes area of extruded faces in addition to vertices and edges.
- Documentation has been revised and edited.

Removed
~~~~~~~

- The symmetry.py module.
- The ft.py module.
- The symmetry.py module.
- The get_params method of TabulatedShapeFamily.
- The family_from_doi method (the underlying data dictionary is now directly exposed).

v0.3.0 - 2020-06-18
-------------------

Added
~~~~~

- Calculation of circumsphere from center for convex polyhedra.
- Simple name-based shape getter for damasceno SHAPES dictionary.
- Polygons moment of inertia calculation.
- Interoperability with the GSD shape specification.
- Shape families and stored data for well-known families.
- All shapes can be centered anywhere in 3D Euclidean space.
- Extensive style checking using black, isort, and various other flake8 plugins.
- Make Circle area settable.
- 3D shapes can be oriented by their principal axes.
- Make Sphere volume settable.

Changed
~~~~~~~

- Inertia tensors for polyhedra and moments of inertia for polygons are calculated in global coordinates rather than the body frame.
- Modified testing of convex hulls to generate points on ellipsoids to avoid degenerate simplices.
- All insphere, circumsphere, and bounding sphere calculations now return the appropriate classes instead of tuples.

Removed
~~~~~~~

- The common_shapes subpackage.

v0.2.0 - 2020-04-09
-------------------

Added
~~~~~

- Continuous integrated testing on CircleCI.
- New Polygon class with property-based API.
- New ConvexSpheropolygon class with property-based API.
- New Polyhedron class with property-based API and robust facet sorting and merging.
- New ConvexPolyhedron class with property-based API.
- New ConvexSpheropolyhedron class with property-based API.
- Ability to plot Polyhedra and Polygons.
- Can now check whether points lie inside a ConvexPolyhedron or ConvexSpheropolyhedron.
- Added documentation.
- New Ellipsoid class with property-based API.
- New Sphere class with property-based API.
- New Ellipse class with property-based API.
- New Circle class with property-based API.
- Added insphere from center calculation for convex polyhedra.
- New ConvexPolygon class.
- Documentation is hosted on ReadTheDocs.

Changed
~~~~~~~

- Moved core shape classes from euclid.FreudShape into top-level package namespace.
- Moved common shape definitions into common_shapes subpackage.
- Shapes from Damasceno science 2012 paper are now stored in a JSON file that is loaded in the damasceno module.

Fixed
~~~~~

- Formatting now properly follows PEP8.

Removed
~~~~~~~

- Various unused or redundant functions in the utils module.
- The quaternion_tools module (uses rowan for quaternion math instead).
- The shapelib module.
- Old polygon.py and polyhedron.py modules, which contained old implementations of various poly\* and spheropoly\* classes.

v0.1.0
------

- Initial version of code base.
coxeter
=======

.. contents::
   :local:

|JOSS|
|ReadTheDocs|
|CircleCI|
|PyPI|
|conda-forge|

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.03098/status.svg
   :target: https://doi.org/10.21105/joss.03098
.. |ReadTheDocs| image:: https://readthedocs.org/projects/coxeter/badge/?version=latest
    :target: http://coxeter.readthedocs.io/en/latest/?badge=latest
.. |CircleCI| image:: https://circleci.com/gh/glotzerlab/coxeter.svg?style=svg
    :target: https://circleci.com/gh/glotzerlab/coxeter
.. |PyPI| image:: https://img.shields.io/pypi/v/coxeter.svg
    :target: https://pypi.org/project/coxeter/
.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/coxeter.svg
   :target: https://anaconda.org/conda-forge/coxeter

Welcome to the documentation for **coxeter**!
The **coxeter** Python library provides tools for working with common geometric objects in two and three dimensions.
Named for the `20th century geometer <https://en.wikipedia.org/wiki/Harold_Scott_MacDonald_Coxeter>`__ best known for his work on polytopes, **coxeter** is especially focused on polygons and polyhedra, but it also support various standard curved shapes such as spheres and ellipsoids.

The package emphasizes working with shapes as mutable objects whose geometric attributes may be accessed using property-based APIs.
Since **coxeter** originally arose to support representations of anisotropic nanoparticles, many shapes support calculations of physical properties (such as form factors and inertia tensors) in addition to purely geometric ones.
However, the package is designed with more general audiences in mind as well, and it aims to support precise calculations of a wide range of geometric quantities that are useful in a number of fields.

Some core features of **coxeter** include:

* Libraries of common shapes to support easy construction.
* Mutable shape objects that can be rescaled in a variety of ways to suit a number of needs.
* Immediate access to geometric properties of shapes via Python properties of shape objects.
* Plotting functionality to make it easy to visualize shapes in both two and three dimensions.

More detailed information on **coxeter**'s features and examples of how to use them may be found in the `documentation <https://coxeter.readthedocs.io/>`__.

.. _installing:

Setup
-----

The recommended methods for installing coxeter are using **pip** or **conda**.

Installation via pip
~~~~~~~~~~~~~~~~~~~~

To install the package from PyPI, execute:

.. code:: bash

   pip install coxeter --user

Installation via conda
~~~~~~~~~~~~~~~~~~~~~~

To install the package from conda, first add the **conda-forge** channel:

.. code:: bash

   conda config --add channels conda-forge

After the **conda-forge** channel has been added, you can install coxeter by executing

.. code:: bash

   conda install coxeter

Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

To install from source, execute:

.. code:: bash

   git clone https://github.com/glotzerlab/coxeter.git
   cd coxeter
   python setup.py install --user

Requirements
~~~~~~~~~~~~

-  Python >= 3.6
-  NumPy >= 1.15
-  SciPy >= 1.0.0
-  rowan >= 1.2

Testing
-------

The package is currently tested for Python >= 3.6 on Unix-like systems.
Continuous integrated testing is performed using CircleCI on these Python versions.

To run the packaged unit tests, execute the following line from the root of the repository:

.. code:: bash

   pytest

To check test coverage, make sure the coverage module is installed:

.. code:: bash

   pip install coverage

and then run the packaged unit tests with the coverage module:

.. code:: bash

   pytest --cov=coxeter

Building Documentation
----------------------

Documentation for coxeter is written in `reStructuredText <http://docutils.sourceforge.net/rst.html>`__ and compiled using `Sphinx <http://www.sphinx-doc.org/en/master/>`__.
To build the documentation, first install Sphinx:

.. code:: bash

   cd doc
   pip install -r requirements.txt

You can then use Sphinx to create the actual documentation in either PDF or HTML form by running the following commands in the coxeter root directory:

.. code:: bash

   make html # For html output
   make latexpdf # For a LaTeX compiled PDF file
   open build/html/index.html

Support and Contribution
========================

This package is hosted on `GitHub <https://github.com/glotzerlab/coxeter>`_.
Please report any bugs or problems that you find on the `issue tracker <https://github.com/glotzerlab/coxeter/issues>`_.
All contributions to coxeter are welcomed via pull requests!
.. _development:

=================
Development Guide
=================


All contributions to **coxeter** are welcome!
Developers are invited to contribute to the framework by pull request to the package repository on `GitHub`_, and all users are welcome to provide contributions in the form of **user feedback** and **bug reports**.
We recommend discussing new features in form of a proposal on the issue tracker for the appropriate project prior to development.

.. _github: https://github.com/glotzerlab/coxeter

General Guidelines
==================

All code contributed to **coxeter** must adhere to the following guidelines:

  * Use the OneFlow_ model of development:
    - Both new features and bug fixes should be developed in branches based on ``master``.
    - Hotfixes (critical bugs that need to be released *fast*) should be developed in a branch based on the latest tagged release.
  * Avoid external dependencies wherever possible, and avoid introducing **any** hard dependencies outside the standard Python scientific stack (NumPy, SciPy, etc). Soft dependencies are allowed for specific functionality, but such dependencies cannot impede the installation of **coxeter** or the use of any other features.
  * All code should adhere to the source code conventions and satisfy the documentation and testing requirements discussed below.
  * Preserve backwards-compatibility whenever possible. Make clear if something must change, and notify package maintainers that merging such changes will require a major release.

To provide a reasonable balance between a high level of backwards compatibility and a reasonable maintenance burden, **coxeter** has adopted `NEP 29`_ to limit the Python and NumPy versions that will be supported.


.. tip::

    During continuous integration, the code is checked automatically with `flake8`_, including a number of plugins that validate parts of the style.
    To run this locally, you can install and run flake8 locally:

    .. code-block:: bash

        python -m pip install flake8 flake8-black flake8-bugbear flake8-docstrings flake8-rst-docstrings pep8-naming flake8-isort
        python -m flake8 coxeter tests

    To avoid having commits fail in case you forget to run this, you can set up a git pre-commit hook using `pre-commit`_:

    .. code-block:: bash

        python -m pip install pre-commit
        pre-commit install

.. _OneFlow: https://www.endoflineblog.com/oneflow-a-git-branching-model-and-workflow
.. _flake8: http://flake8.pycqa.org/en/latest/
.. _pre-commit: https://pre-commit.com/
.. _NEP 29: https://numpy.org/neps/nep-0029-deprecation_policy.html


Style Guidelines
----------------

The **coxeter** package adheres to a relatively strict set of style guidelines.
All code in **coxeter** should be formatted using `black`_.
Imports should be formatted using `isort`_.
For guidance on the style, see `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ and the `Google Python Style Guide <https://google.github.io/styleguide/pyguide.html>`_, but any ambiguities should be resolved automatically by running black.
All code should of course also follow the principles in `PEP 20 <https://www.python.org/dev/peps/pep-0020/>`_.

.. tip::

    Developers should format their code using black and isort locally. Running the pre-commit hooks will take care of this:

    .. code-block:: bash

        pre-commit run

    Alternatively, the tools can be run manually using the commands:

    .. code-block:: bash

        black coxeter/ tests/
        isort coxeter/ tests/

.. _black: https://black.readthedocs.io/
.. _isort: https://timothycrosley.github.io/isort/


Documentation
-------------

API documentation should be written as part of the docstrings of the package in the `Google style <https://google.github.io/styleguide/pyguide.html#383-functions-and-methods>`__.
There is one notable exception to the guide: class properties should be documented in the getters functions, not as class attributes, to allow for more useful help messages and inheritance of docstrings.
Docstrings may be validated using `pydocstyle <http://www.pydocstyle.org/>`__ (or using the flake8-docstrings plugin as documented above).
The `official documentation <https://coxeter.readthedocs.io/>`_ is generated from the docstrings using `Sphinx <http://www.sphinx-doc.org/en/stable/index.html>`_.

In addition to API documentation, inline comments are **highly encouraged**.
Code should be written as transparently as possible, so the primary goal of documentation should be explaining the algorithms or mathematical concepts underlying the code.
Avoid comments that simply restate the nature of lines of code.
For example, the comment "solve the system of equations" is uninformative, since the code itself should make this obvious, *e.g*, ``np.linalg.solve``.
On the other hand, the comment "the solution to this system of equations is the intersection of truncation planes" is instructive.


Unit Tests
----------

All code should include a set of tests which test for correct behavior.
All tests should be placed in the ``tests`` folder at the root of the project.
In general, most parts of coxeter primarily require `unit tests <https://en.wikipedia.org/wiki/Unit_testing>`_, but where appropriate `integration tests <https://en.wikipedia.org/wiki/Integration_testing>`_ are also welcome.
Tests in **coxter** use the `pytest <https://docs.pytest.org/>`__ testing framework.
To run the tests, simply execute ``pytest`` at the root of the repository.


Release Guide
=============

To make a new release of **coxeter**, follow the following steps:

#. Make a new branch off of master based on the expected new version, *e.g.*
   release-2.3.1.
#. Make any final changes as desired on this branch. Push the changes and
   ensure all tests are passing as expected on the new branch.
#. Once the branch is completely finalized, run bumpversion with the
   appropriate type (patch, minor, major) so that the version now matches the
   version number in the branch name.
#. Merge the branch back into master, then push master and push tags. The
   tagged commit will automatically trigger generation of binaries and upload
   to PyPI and conda-forge.
#. Delete the release branch both locally and on the remote.
==========
References
==========

.. It's necessary to have a file that creates the bibliography that is
   processed by Sphinx AFTER all of the :cite: commands have been parsed. To
   fix this, we create this dummy file that comes as late as possible and
   import the bibliography here.

.. bibliography:: coxeter.bib
.. _examples:

Examples
========

While **coxeter**'s API is very easy to work with, it can be helpful to see some real demonstrations of how it can be used.
Here we include some practical examples of using **coxeter**.

.. toctree::
   examples/PointInShape.ipynb
   examples/DistanceToSurface.ipynb
   examples/InertiaTensors.ipynb
   examples/Spheres.ipynb
   examples/FormFactors.ipynb
coxeter.families package
===============================

.. rubric:: Overview

.. automodule:: coxeter.families
   :members:
   :show-inheritance:
coxeter.shapes package
==============================

.. rubric:: Overview

.. automodule:: coxeter.shapes
   :members:
   :show-inheritance:
Credits
=======

.. include:: ../../Credits.rst
.. _quickstart:

Quickstart Tutorial
===================

Once you have :ref:`installed <installing>` **coxeter**, most workflows involve creating an instance of a shape, such as a :class:`~coxeter.shapes.Polygon`:

.. code-block:: python

    >>> import coxeter
    >>> square = coxeter.shapes.Polygon([[0, 0], [1, 0], [1, 1], [0, 1]])

All shapes may be found in the `coxeter.shapes` subpackage and are created in a similar manner.
For instance, making a sphere requires a radius: ``sphere = coxeter.shapes.Sphere(3)``.
Once you have a shape, you can immediately query it for properties.

.. code-block:: python

    >>> square.vertices
    array([[0., 0., 0.],
           [1., 0., 0.],
           [1., 1., 0.],
           [0., 1., 0.]])
    >>> square.area
    1.0
    >>> square.perimeter
    4.0

The **coxeter** library comes with a range of *shape families*, collections of shapes with standard definitions so that you don't have to parameterize them yourself.
The regular :math:`n`-gon family is such an example, provided as :class:`coxeter.families.RegularNGonFamily`.

.. code-block:: python

    >>> hexagon = coxeter.families.RegularNGonFamily.get_shape(6)
    >>> hexagon.vertices.round(2)
    array([[ 0.62,  0.  ,  0.  ],
           [ 0.31,  0.54,  0.  ],
           [-0.31,  0.54,  0.  ],
           [-0.62,  0.  ,  0.  ],
           [-0.31, -0.54,  0.  ],
           [ 0.31, -0.54,  0.  ]])

Part of what makes **coxeter** so powerful is that all shapes are mutable.
This means that once you have a prototype of a shape, it can be modified to fit a specific need.
For example, the snippet below finds the area of the smallest regular pentagon that contains an equilateral triangle of unit area:

.. code-block:: python

    >>> triangle = coxeter.families.RegularNGonFamily.get_shape(3)
    >>> pentagon = coxeter.families.RegularNGonFamily.get_shape(5)
    >>> pentagon.incircle_radius = triangle.circumcircle.radius
    >>> triangle.area
    0.9999999999999999
    >>> triangle.circumcircle.area
    2.418399152312292
    >>> pentagon.area
    2.796463494144044

This tutorial just scratches the surface of the features **coxeter** offers.
For more complete demonstrations of the package's features, see the :ref:`examples`.
coxeter.shape\_getters module
=============================

.. automodule:: coxeter.shape_getters
   :members: from_gsd_type_shapes
   :undoc-members:
   :show-inheritance:
Table of Contents
=================

.. include:: ../../README.rst


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   coxeter
   examples

.. toctree::
   :maxdepth: 1
   :caption: Reference:

   development
   changelog
   credits
   license
   zreferences


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
API Reference
=============

Subpackages
-----------

.. toctree::
   package-shapes
   package-families

Submodules
----------

.. toctree::

   module-shape-getters

.. automodule:: coxeter
   :autosummary:
   :members:
   :show-inheritance:
Changelog
=========

.. include:: ../../ChangeLog.rst
