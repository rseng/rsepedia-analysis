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
