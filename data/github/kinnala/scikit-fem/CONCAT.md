# scikit-fem

`scikit-fem` is a pure Python 3.7+ library for performing [finite element
assembly](https://en.wikipedia.org/wiki/Finite_element_method). Its main
purpose is the transformation of bilinear forms into sparse matrices and linear
forms into vectors.

<a href="https://colab.research.google.com/github/kinnala/scikit-fem-notebooks/blob/master/ex1.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg"></a>
<a href="https://pypi.org/project/scikit-fem/" alt="PyPI"><img src="https://img.shields.io/pypi/v/scikit-fem" /></a>
<a href="https://scikit-fem.readthedocs.io/" alt="Documentation"><img src="https://readthedocs.org/projects/pip/badge/?version=stable" /></a>
<a href="https://joss.theoj.org/papers/4120aba1525403e6d0972f4270d7b61e" alt="status"><img src="https://joss.theoj.org/papers/4120aba1525403e6d0972f4270d7b61e/status.svg" /></a>

The library

- has minimal dependencies
- contains no compiled code
- supports one-dimensional, triangular, quadrilateral, tetrahedral and hexahedral finite elements
- includes special elements such as Raviart-Thomas, Nédélec, MINI, Crouzeix-Raviart, Argyris, ...

## Installation

The most recent release can be installed simply by
```
pip install scikit-fem[all]
```
Remove `[all]` to not install the optional dependencies `meshio` for mesh
input/output, and `matplotlib` for creating simple visualizations.
The minimal dependencies are `numpy` and `scipy`.
You can also try the library in browser through [Google Colab](https://colab.research.google.com/github/kinnala/scikit-fem-notebooks/blob/master/ex1.ipynb).

## Examples

Solve the Poisson problem (see also [`ex01.py`](https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex01.py)):
```python
from skfem import *
from skfem.helpers import dot, grad

# create the mesh
m = MeshTri().refined(4)
# or, with your own points and elements:
# m = MeshTri(points, elements)

e = ElementTriP1()
basis = Basis(m, e)  # shorthand for CellBasis

@BilinearForm
def laplace(u, v, _):
    return dot(grad(u), grad(v))

@LinearForm
def rhs(v, _):
    return 1.0 * v

A = laplace.assemble(basis)
b = rhs.assemble(basis)

# Dirichlet boundary conditions
A, b = enforce(A, b, D=m.boundary_nodes())

# solve the linear system
x = solve(A, b)

# plot the solution
from skfem.visuals.matplotlib import plot, savefig
plot(m, x, shading='gouraud', colorbar=True)
savefig('solution.png')
```

Meshes can be initialized manually, loaded from external files using
[meshio](https://github.com/nschloe/meshio), or created with the help of special
constructors:

```python
import numpy as np
from skfem import MeshLine, MeshTri, MeshTet

mesh = MeshLine(np.array([0., .5, 1.]))
mesh = MeshTri(
    np.array([[0., 0.],
              [1., 0.],
              [0., 1.]]).T,
    np.array([[0, 1, 2]]).T,
)
mesh = MeshTri.load("docs/examples/meshes/square.msh")  # requires meshio
mesh = MeshTet.init_tensor(*((np.linspace(0, 1, 60),) * 3))
```

We support [many common finite
elements](https://github.com/kinnala/scikit-fem/blob/master/skfem/element/__init__.py#L51).
Below the stiffness matrix is assembled using second-order tetrahedra:

```python
from skfem import Basis, ElementTetP2

basis = Basis(mesh, ElementTetP2())  # quadratic tetrahedron
A = laplace.assemble(basis)  # type: scipy.sparse.csr_matrix
```

More examples can be found in the [gallery](https://scikit-fem.readthedocs.io/en/latest/listofexamples.html).


## Benchmark

*The following benchmark (`docs/examples/performance.py`) demonstrates the time
spent on finite element assembly in comparison to the time spent on linear
solve.  The given numbers were calculated using a ThinkPad X1 Carbon laptop (7th
gen).  Note that the timings are only illustrative as they depend on, e.g., the
type of element used, the number of quadrature points used, the type of linear
solver, and the complexity of the forms.  This benchmark solves the Laplace
equation using linear tetrahedral elements and the default direct sparse solver
of `scipy.sparse.linalg.spsolve`.*

| Degrees-of-freedom | Assembly (s) | Linear solve (s) |
| --- | --- | --- |
| 4096 | 0.04805 | 0.04241 |
| 8000 | 0.09804 | 0.16269 |
| 15625 | 0.20347 | 0.87741 |
| 32768 | 0.46399 | 5.98163 |
| 64000 | 1.00143 | 36.47855 |
| 125000 | 2.05274 | nan |
| 262144 | 4.48825 | nan |
| 512000 | 8.82814 | nan |
| 1030301 | 18.25461 | nan |


## Documentation

The project is documented using Sphinx under `docs/`.
Built version can be found from [Read the Docs](https://scikit-fem.readthedocs.io/en/latest/).
Here are direct links to additional resources:

- [Examples from our test suite](https://scikit-fem.readthedocs.io/en/latest/listofexamples.html)
- [Examples from the FEniCS tutorial](https://github.com/gdmcbain/fenics-tuto-in-skfem)

## Getting help

If you encounter an issue you can use GitHub issue tracker.  If you cannot find help from the documentation,
you can use the GitHub Discussions to [ask questions](https://github.com/kinnala/scikit-fem/discussions).
Try to provide a snippet of code which fails
and include also the version of the library you are
using.  The version can be found as follows:
```python
import skfem; print(skfem.__version__)
```

## Dependencies

The minimal dependencies for installing `scikit-fem` are
[numpy](https://numpy.org/) and [scipy](https://www.scipy.org/).  In addition,
many
[examples](https://scikit-fem.readthedocs.io/en/latest/listofexamples.html) use
[matplotlib](https://matplotlib.org/) for visualization and
[meshio](https://github.com/nschloe/meshio) for loading/saving different mesh
file formats.  Some examples demonstrate the use of other external packages;
see `requirements.txt` for a list of test dependencies.

## Testing

The tests are run by GitHub Actions.  The `Makefile` in the repository root has
targets for running the testing container locally using `docker`.  For example,
`make test_py38` runs the tests using `py38` branch from
[kinnala/scikit-fem-docker-action](https://github.com/kinnala/scikit-fem-docker-action).
The releases are tested in
[kinnala/scikit-fem-release-tests](https://github.com/kinnala/scikit-fem-release-tests).

## Licensing

The contents of `skfem/` and the PyPI package `scikit-fem` are licensed under
the 3-clause BSD license.  Some examples under `docs/examples/` or snippets
in the documentation may have a different license.

## Acknowledgements

This project was started while working under a grant from the [Finnish Cultural
Foundation](https://skr.fi/).  Versions 2.0.0+ were prepared while working in a
project funded by the [Academy of
Finland](https://akareport.aka.fi/ibi_apps/WFServlet?IBIF_ex=x_HakKuvaus2&CLICKED_ON=&HAKNRO1=324611&UILANG=en).
The approach used in the finite element assembly has been inspired by the [work
of A. Hannukainen and
M. Juntunen](https://au.mathworks.com/matlabcentral/fileexchange/36108-hjfem_lite).

## Contributing

We are happy to welcome any contributions to the library.  Reasonable projects
for first timers include:

- Reporting a [bug](https://github.com/kinnala/scikit-fem/issues)
- Writing an [example](https://github.com/kinnala/scikit-fem/tree/master/docs/examples)
- Improving the [tests](https://github.com/kinnala/scikit-fem/tree/master/tests)
- Finding typos in the documentation.

*By contributing code to scikit-fem, you are agreeing to release it under BSD-3-Clause, see LICENSE.md.*

## Citing the library

We appreciate if you cite the following article in your scientific works:
```
@article{skfem2020,
  doi = {10.21105/joss.02369},
  year = {2020},
  volume = {5},
  number = {52},
  pages = {2369},
  author = {Tom Gustafsson and G. D. McBain},
  title = {scikit-fem: A {P}ython package for finite element assembly},
  journal = {Journal of Open Source Software}
}
```

## Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
with respect to documented and/or tested features.

### Unreleased

- Changed: `DiscreteField` is now a subclass of `ndarray` instead of
  `NamedTuple` and, consequently, the components of `DiscreteField` cannot no
  more be indexed inside forms like `u[1]` (use `u.grad` instead)
- Changed: Writing `w['u']` and `w.u` inside the form definition is now
  equivalent (previously `w.u == w['u'].value`)
- Changed: `Mesh.draw` now uses `matplotlib` by default, calling
  `Mesh.draw("vedo")` uses `vedo`
- Changed: `skfem.visuals.matplotlib` now uses `jet` as the default colormap
- Changed: `BoundaryFacetBasis` is now an alias of `FacetBasis` instead of
  other way around
- Deprecated: `DiscreteField.value` remains for backwards-compatibility but is
  now deprecated and can be dropped
- Added: `Mesh.plot`, a wrapper to `skfem.visuals.*.plot`
- Added: `Basis.plot`, a wrapper to `skfem.visuals.*.plot`
- Added: `Basis.refinterp` now supports vectorial fields
- Added: `skfem.visuals.matplotlib.plot` now has a basic quiver plot for vector
  fields
- Added: `Mesh.facets_around` which constructs a set of facets around a
  subdomain
- Added: `Mesh.load` now tries loading the orientation of boundaries and
  interfaces
- Added: `OrientedBoundary` which is a subclass of `ndarray` for facet index
  arrays with the orientation information (0 or 1 per facet) available as
  `OrientedBoundary.ori`
- Added: `FacetBasis` will use the facet orientations (if present) to calculate
  traces and normal vectors
- Added: `skfem.visuals.matplotlib.draw` will visualize the orientations if
  `boundaries=True` is given
- Added: `Mesh.facets_satisfying` allows specifying the keyword argument
  `normal` for orienting the resulting interface
- Added: `FacetBasis` constructor now has the keyword argument `side` which
  allows changing the side of the facet used to calculate the basis function
  values and gradients
- Fixed: Improvements to backwards compatibility in `asm`/`assemble` keyword
  arguments
- Fixed: Save format issue with meshio 5.3.0+
- Fixed: `CellBasis` did not properly support `elements` argument
- Fixed: `Basis.interpolate` did not properly interpolate all components of
  `ElementComposite`

### [5.2.0] - 2021-12-27

- Added: `Basis.project`, a more general and easy to use alternative for
  `projection`
- Added: `Basis` and `FacetBasis` kwargs `elements` and `facets` can now be a
  string refering to subdomain and boundary tags
- Added: `ElementQuadRT0`, lowest-order quadrilateral Raviart-Thomas element
- Fixed: `Functional` returned only the first component for forms with
  non-scalar output

### [5.1.0] - 2021-11-30

- Added: `skfem.helpers.mul` for matrix multiplication
- Added: `Basis.split` will now also split `ElementVector` into its components
- Fixed: `ElementDG` was not included in the wildcard import
- Fixed: Automatic visualization of `MeshTri2` and `MeshQuad2` in Jupyter
  notebooks raised exception

### [5.0.0] - 2021-11-21

- Changed: `meshio` is now an optional dependency
- Changed: `ElementComposite` uses `DiscreteField()` to represent zero
- Added: Support more argument types in `Basis.get_dofs`
- Added: Version information in `skfem.__version__`
- Added: Preserve `Mesh.boundaries` during uniform refinement of `MeshLine1`,
  `MeshTri1` and `MeshQuad1`
- Fixed: Refinement of quadratic meshes will now fall back to the refinement
  algorithm of first-order meshes instead of crashing
- Fixed: Edge cases in the adaptive refine of `MeshTet1` that failed to produce
  a valid mesh
- Deprecated: `Basis.find_dofs` in favor of `Basis.get_dofs`
- Deprecated: Merging `DofsView` objects via `+` and `|` because of surprising
  behavior in some edge cases

### [4.0.1] - 2021-10-15

- Fixed: `MappingIsoparametric` can now be pickled

### [4.0.0] - 2021-09-27

- Added: `Mesh.save`/`Mesh.load` now exports/imports `Mesh.subdomains` and
  `Mesh.boundaries`
- Added: `Mesh.load` now optionally writes any mesh data to a list passed via
  the keyword argument `out`, e.g., `out=data` where `data = ['point_data']`
- Added: `Mesh.load` (and `skfem.io.meshio.from_file`) now supports the
  additional keyword argument `force_meshio_type` for loading mesh files that
  have multiple element types written in the same file, one element type at
  a time
- Added: `asm` will now accept a list of bases, assemble the same form using
  all of the bases and sum the result (useful for jump terms and mixed meshes,
  see Example 41)
- Added: `Mesh.with_boundaries` now allows the definition of internal boundaries/interfaces
  via the flag `boundaries_only=False`
- Added: `MeshTri1DG`, `MeshQuad1DG`, `MeshHex1DG`, `MeshLine1DG`; new mesh
  types for describing meshes with a discontinuous topology, e.g., periodic
  meshes (see Example 42)
- Added: `ElementHexDG` for transforming hexahedral H1 elements to DG/L2 elements.
- Added: `ElementTriP1DG`, `ElementQuad1DG`, `ElementHex1DG`,
  `ElementLineP1DG`; shorthands for `ElementTriDG(ElementTriP1())` etc.
- Added: `ElementTriSkeletonP0` and `ElementTriSkeletonP1` for defining
  Lagrange multipliers on the skeleton mesh (see Example 40)
- Added: `TrilinearForm` for assembling a sparse 3-tensor, e.g., when dealing
  with unknown material data
- Added: `MeshTri.oriented` for CCW oriented triangular meshes which can be
  useful for debugging or interfacing to external tools
- Added: partial support for `MeshWedge1` and `ElementWedge1`, the lowest order
  wedge mesh and element
- Added: `ElementTriP3`, cubic triangular Lagrange element
- Added: `ElementTriP4`, quartic triangular Lagrange element
- Added: `ElementTri15ParamPlate`, 15-parameter nonconforming triangular element for plates
- Added: `ElementTriBDM1`, the lowest order Brezzi-Douglas-Marini element
- Added: `Mesh.draw().show()` will now visualize any mesh interactively (requires [vedo](https://vedo.embl.es/))
- Added: Adaptive refinement for `MeshTet1`
- Fixed: `MappingIsoparametric` is now about 2x faster for large meshes thanks
  to additional caching
- Fixed: `MeshHex2.save` did not work properly
- Fixed: `Mesh.load` ignores unparseable `cell_sets` inserted by `meshio` in MSH 4.1
- Changed: `Mesh` string representation is now more informative
- Changed: `Form.assemble` no more allows keyword arguments with `list` or
  `dict` type: from now on only `DiscreteField` or 1d/2d `ndarray` objects are
  allowed and 1d `ndarray` is passed automatically to `Basis.interpolate` for
  convenience
- Changed: `MeshLine` is now a function which initializes `MeshLine1`
  and not an alias to `MeshLine1`
- Changed: `FacetBasis` is now a shorthand for `BoundaryFacetBasis` and no
  longer initializes `InteriorFacetBasis` or `MortarFacetBasis` if the keyword
  argument `side` is passed to the constructor
- Removed: the deprecated `Mesh.define_boundary` method

### [3.2.0] - 2021-08-02

- Added: `ElementTriCCR` and `ElementTetCCR`, conforming Crouzeix-Raviart finite elements
- Fixed: `Mesh.mirrored` returned a wrong mesh when a point other than the origin was used
- Fixed: `MeshLine` constructor accepted only numpy arrays and not plain Python lists
- Fixed: `Mesh.element_finder` (and `CellBasis.probes`, `CellBasis.interpolator`) was not working properly for a small number of elements (<5) or a large number of input points (>1000)
- Fixed: `MeshTet` and `MeshTri.element_finder` are now more robust against degenerate elements
- Fixed: `Mesh.element_finder` (and `CellBasis.probes`, `CellBasis.interpolator`) raises exception if the query point is outside of the domain

### [3.1.0] - 2021-06-18

- Added: `Basis`, a shorthand for `CellBasis`
- Added: `CellBasis`, a new preferred name for `InteriorBasis`
- Added: `BoundaryFacetBasis`, a new preferred name for `ExteriorFacetBasis`
- Added: `utils.penalize`, an alternative to `condense` and `enforce` for
  essential boundary conditions
- Added: `InteriorBasis.point_source`, with `ex38`
- Added: `ElementTetDG`, similar to `ElementTriDG` for tetrahedral meshes
- Fixed: `MeshLine1.element_finder` 

### [3.0.0] - 2021-04-19

- Added: Completely rewritten `Mesh` base class which is "immutable" and uses
  `Element` classes to define the ordering of nodes; better support for
  high-order and other more general mesh types in the future
- Added: New quadratic mesh types: `MeshTri2`, `MeshQuad2`, `MeshTet2` and `MeshHex2`
- Added: `InteriorBasis.probes`; like `InteriorBasis.interpolator` but returns a matrix
  that operates on solution vectors to interpolate them at the given points
- Added: More overloads for `DiscreteField`, e.g., multiplication, summation
  and subtraction are now explicitly supported inside the form definitions
- Added: `MeshHex.to_meshtet` for splitting hexahedra into tetrahedra
- Added: `MeshHex.element_finder` for interpolating finite element solutions
  on hexahedral meshes via `InteriorBasis.interpolator`
- Added: `Mesh.with_boundaries`, a functional replacement to
  `Mesh.define_boundary`, i.e. defining boundaries via Boolean lambda function
- Added: `Mesh.with_subdomains` for defining subdomains via Boolean lambda function
- Added: `skfem.utils.projection`, a replacement of `skfem.utils.project`
  with a different, more intuitive order of arguments
- Added: `skfem.utils.enforce` for setting essential boundary conditions by
  changing matrix rows to zero and diagonals to one.
- Deprecated: `skfem.utils.project` in favor of `skfem.utils.projection`
- Deprecated: `Mesh.define_boundary` in favor of `Mesh.with_boundaries`
- Removed: `Mesh.{refine,scale,translate}`; the replacements are `Mesh.{refined,scaled,translated}`
- Removed: `skfem.models.helpers`; available as `skfem.helpers`
- Removed: `DiscreteField.{f,df,ddf,hod}`; available as `DiscreteField.{value,grad,hess,grad3,...}`
- Removed: Python 3.6 support
- Removed: `skfem.utils.L2_projection`
- Removed: `skfem.utils.derivative`
- Changed: `Mesh.refined` no more attempts to fix the indexing of `Mesh.boundaries` after refine
- Changed: `skfem.utils.solve` now uses `scipy.sparse.eigs` instead of `scipy.sparse.eigsh` by default;
  the old behavior can be retained by explicitly passing `solver=solver_scipy_eigs_sym()`
- Fixed: High memory usage in `skfem.visuals.matplotlib` related to 1D plotting

### [2.5.0] - 2021-02-13

- Deprecated: `side` keyword argument to `FacetBasis` in favor of the more
  explicit `InteriorFacetBasis` and `MortarFacetBasis`.
- Added: `InteriorFacetBasis` for integrating over the interior facets, e.g.,
  evaluating error estimators with jumps and implementing DG methods.
- Added: `MortarFacetBasis` for integrating over the mortar mesh.
- Added: `InteriorBasis.with_element` for reinitializing an equivalent basis
  that uses a different element.
- Added: `Form.partial` for applying `functools.partial` to the form function
  wrapped by `Form`.
- Fixed: Include explicit Python 3.9 support.

### [2.4.0] - 2021-01-20

- Deprecated: List and tuple keyword argument types to `asm`.
- Deprecated: `Mesh2D.mirror` in favor of the more general `Mesh.mirrored`.
- Deprecated: `Mesh.refine`, `Mesh.scale` and `Mesh.translate` in favor of
  `Mesh.refined`, `Mesh.scaled` and `Mesh.translated`.
- Added: `Mesh.refined`, `Mesh.scaled`, and `Mesh.translated`. The new methods
  return a copy instead of modifying `self`.
- Added: `Mesh.mirrored` for mirroring a mesh using a normal and a point.
- Added: `Functional` now supports forms that evaluate to vectors or other
  tensors.
- Added: `ElementHex0`, piecewise constant element for hexahedral meshes.
- Added: `FacetBasis.trace` for restricting existing solutions to lower
  dimensional meshes on boundaries or interfaces.
- Fixed: `MeshLine.refined` now correctly performs adaptive refinement of
  one-dimensional meshes.

### [2.3.0] - 2020-11-24

- Added: `ElementLineP0`, one-dimensional piecewise constant element.
- Added: `skfem.helpers.curl` now calculates the rotated gradient for
  two-dimensional elements.
- Added: `MeshTet.init_ball` for meshing a ball.
- Fixed: `ElementQuad0` was not compatible with `FacetBasis`.

### [2.2.3] - 2020-10-16

- Fixed: Remove an unnecessary dependency.

### [2.2.2] - 2020-10-15

- Fixed: Make the preconditioner in `TestEx32` more robust.

### [2.2.1] - 2020-10-15

- Fixed: Remove `tests` from the PyPI distribution.

### [2.2.0] - 2020-10-14

- Deprecated: `L2_projection` will be replaced by `project`.
- Deprecated: `derivative` will be replaced by `project`.
- Added: `MeshTet.element_finder` and `MeshLine.element_finder` for using
  `InteriorBasis.interpolator`.
- Added: `ElementTriCR`, the nonconforming Crouzeix-Raviart element for Stokes flow.
- Added: `ElementTetCR`, tetrahedral nonconforming Crouzeix-Raviart element.
- Added: `ElementTriHermite`, an extension of `ElementLineHermite` to triangular
  meshes.
- Fixed: Fix `Mesh.validate` for unsigned `Mesh.t`.

### [2.1.1] - 2020-10-01

- Fixed: Further optimizations to `Mesh3D.boundary_edges`: tested to run on a laptop
  with over 10 million elements.

### [2.1.0] - 2020-09-30

- Added: `ElementHex2`, a triquadratic hexahedral element.
- Added: `MeshTri.init_circle`, constructor for a circle mesh.
- Fixed: `Mesh3D.boundary_edges` (and, consequently, `Basis.find_dofs`) was slow
  and used lots of memory due to an exhaustive search of all edges.

### [2.0.0] - 2020-08-21

- Deprecated: `project` will only support functions like `lambda x: x[0]`
  instead of `lambda x, y, z: x` in the future.
- Added: Support for complex-valued forms: `BilinearForm` and `LinearForm` now take
  an optional argument `dtype` which defaults to `np.float64`
  but can be also `np.complex64`.
- Added: `Dofs.__or__` and `Dofs.__add__`, for merging degree-of-freedom sets
  (i.e. `Dofs` objects) using `|` and `+` operators.
- Added: `Dofs.drop` and `Dofs.keep`, for further filtering the degree-of-freedom sets
- Removed: Support for old-style decorators `bilinear_form`, `linear_form`, and
  `functional` (deprecated since 1.0.0).
- Fixed: `FacetBasis` did not initialize with `ElementQuadP`.

### [1.2.0] - 2020-07-07

- Added: `MeshQuad._splitquads` aliased as `MeshQuad.to_meshtri`.
- Added: `Mesh.__add__`, for merging meshes using `+` operator: duplicated nodes are
  joined.
- Added: `ElementHexS2`, a 20-node quadratic hexahedral serendipity element.
- Added: `ElementLineMini`, MINI-element for one-dimensional mesh.
- Fixed: `Mesh3D.boundary_edges` was broken in case of hexahedral meshes.
- Fixed: `skfem.utils.project` did not work for `ElementGlobal`.

### [1.1.0] - 2020-05-18

- Added: `ElementTetMini`, MINI-element for tetrahedral mesh.
- Fixed: `Mesh3D.boundary_edges` incorrectly returned all edges where both nodes are on
  the boundary.

### [1.0.0] - 2020-04-22

- Deprecated: Old-style form constructors `bilinear_form`, `linear_form`, and `functional`.
- Changed: `Basis.interpolate` returns `DiscreteField` objects instead of ndarray tuples.
- Changed: `Basis.interpolate` works now properly for vectorial and high-order elements
  by interpolating all components and higher order derivatives.
- Changed: `Form.assemble` accepts now any keyword arguments (with type `DiscreteField`)
  that are passed over to the forms.
- Changed: Renamed `skfem.importers` to `skfem.io`.
- Changed: Renamed `skfem.models.helpers` to `skfem.helpers`.
- Changed: `skfem.utils.solve` will now expand also the solutions of eigenvalue problems.
- Added: New-style form constructors `BilinearForm`, `LinearForm`, and `Functional`.
- Added: `skfem.io.json` for serialization of meshes to/from json-files.
- Added: `ElementLinePp`, p-th order one-dimensional elements.
- Added: `ElementQuadP`, p-th order quadrilateral elements.
- Added: `ElementQuadDG` for transforming quadrilateral H1 elements to DG elements.
- Added: `ElementQuadBFS`, Bogner-Fox-Schmit element for biharmonic problems.
- Added: `ElementTriMini`, MINI-element for Stokes problems.
- Added: `ElementComposite` for using multiple elements in one bilinear form.
- Added: `ElementQuadS2`, quadratic Serendipity element.
- Added: `ElementLineHermite`, cubic Hermite element for Euler-Bernoulli beams.
- Added: `Mesh.define_boundary` for defining named boundaries.
- Added: `Basis.find_dofs` for finding degree-of-freedom indices.
- Added: `Mesh.from_basis` for defining high-order meshes.
- Added: `Basis.split` for splitting multicomponent solutions.
- Added: `MortarMapping` with basic support for mortar methods in 2D.
- Added: `Basis` constructors now accept `quadrature` keyword argument for specifying
  a custom quadrature rule.
---
title: 'scikit-fem: A Python package for finite element assembly'
tags:
  - Python
  - numerics
  - finite element method
authors:
  - name: Tom Gustafsson
    orcid: 0000-0003-1611-5032
    affiliation: 1
  - name: G. D. McBain
    orcid: 0000-0002-1904-122X
    affiliation: 2
affiliations:
 - name: Department of Mathematics and Systems Analysis, Aalto University
   index: 1
 - name: Memjet North Ryde Pty Ltd, Macquarie Park, NSW, Australia
   index: 2
date: 4 June 2020
bibliography: paper.bib
---

# Summary

Partial differential equations (PDEs)—such as the Navier–Stokes equations in
fluid mechanics, the Maxwell equations in electromagnetism, and the Schrödinger
equation in quantum mechanics—are the basic building blocks of modern physics
and engineering.  The finite element method (FEM) is a flexible computational
technique for the discretization and solution of PDEs, especially in the case
of complex spatial domains.

Conceptually, the FEM transforms a time-independent (or temporally discretized)
PDE into a system of linear equations $Ax=b$.  `scikit-fem` is a lightweight
Python library for the creation, or *assembly*, of the finite element matrix $A$
and vector $b$.  The user loads a computational mesh, picks suitable basis
functions, and provides the PDE's weak formulation [@fenicsbook].  This results
in sparse matrices and vectors compatible with the SciPy [@scipy] ecosystem.

# Purpose and prior art

There exist several open source packages and frameworks that implement the
finite element method.  `scikit-fem` was developed as a simple and lightweight
alternative to the existing Python packages with a focus on computational
experimentation and custom PDE-based model development.  We rely on pure
interpreted Python code on top of the NumPy–SciPy base which makes `scikit-fem`
easy to install and portable across multiple operating systems.  The reliance on
plain NumPy arrays and SciPy sparse matrices enables interoperability with
various packages in the Python ecosystem such as meshio [@meshio], pacopy
[@pacopy], and pyamg [@pyamg].

In contrast to NGSolve [@ngsolve], FEniCS [@fenics], Firedrake [@firedrake],
SfePy [@sfepy], and GetFEM [@getfem], `scikit-fem` contains no compiled code
making the installation quick and straightforward.  We specifically target
finite element assembly instead of encapsulating the entire finite element
analysis from pre- to postprocessing into a single framework.  As a consequence,
we cannot provide an end-to-end experience when it comes to, e.g., specific
physical models or distributed computing.  Our aim is to be generic in terms of
PDEs and, hence, support a variety of finite element schemes.  Currently
`scikit-fem` includes basic support for $H^1$-, $H(\mathrm{div})$-,
$H(\mathrm{curl})$-, and $H^2$-conforming problems as well as various
nonconforming schemes.

`scikit-fem` accepts weak forms that depend on the values and the derivatives of
the trial and the test functions, their high-order derivatives, the local mesh
parameter, nonuniform material or coefficient fields defined at the quadrature
points, or any existing finite element solutions.  Iterations related to, e.g.,
nonlinear problems (Newton's method and the variants, parameter continuation) or
adaptive mesh refinement (evaluation of functionals and the marking strategy)
should be implemented by the user although we provide basic tools such as
interpolation routines and conforming mesh refinement, and examples by using
them.  The same applies to boundary conditions: the linear system $(A, b)$ is
provided as such and eliminating or penalizing the correct degrees-of-freedom,
implementing inhomogeneous or periodic boundary conditions should be done
separately either by using the various helper routines of `scikit-fem` or by
other means.  `scikit-fem` has no explicit support for distributed computing
although it could be used as a building block in parallel computations such as
parameter sweeps or domain decomposition techniques.

# Examples and enabled work

The documentation of `scikit-fem` contains over 30 examples that demonstrate the
library and its use. The results of some of the examples are highlighted in
\autoref{fig:examples}.  Several publications already utilize computational
results from `scikit-fem`, e.g., @mcbain2018, @gustafsson2019, and
@gustafsson2020.  In addition, `scikit-fem` is used in a recently published
Python package for battery modelling [@pybamm].

![(Top left.) A combination of triangular and quadrilateral elements is used to solve the linear elastic contact problem. (Top right.) The lowest order tetrahedral Nédélec element is used to solve the $H(\mathrm{curl})$-conforming model problem $\nabla \times \nabla \times E + E = f$. The color corresponds to the magnitude of the vector field $E$. (Bottom.) The Taylor–Hood element is used to solve the Navier–Stokes flow over a backward-facing step for different Reynolds numbers.  The first and last two-dimensional figures were generated using `scikit-fem`'s wrapping of matplotlib [@hunter2007]; three-dimensional postprocessing is more specialized and usually left to export of, e.g., VTK or XDMF formats, via meshio [@meshio] for subsequent rendering in, e.g., ParaView [@ahrens2005] as done in the second figure.\label{fig:examples}](examples.png)

# Acknowledgements

Tom Gustafsson has received external funding from the Finnish Cultural
Foundation and the Academy of Finland (decision nr. 324611) while working on the
project.

# References
===============
Advanced topics
===============

This section contains advanced discussions around the features of scikit-fem
with an aim to develop a more detailed understanding of the library.

.. _forms:

Anatomy of forms
================

We consider forms as the basic building blocks of finite element assembly.
Thus, it is useful to understand how forms are used in scikit-fem and how to
express them correctly.

The bilinear form corresponding to the Laplace
operator :math:`-\Delta` is

.. math::

   a(u, v) = \int_\Omega \nabla u \cdot \nabla v \,\mathrm{d}x.

In order to express this in scikit-fem, we write the integrand as a Python
function:

.. doctest::

   >>> from skfem import BilinearForm
   >>> from skfem.helpers import grad, dot
   >>> @BilinearForm
   ... def integrand(u, v, w):
   ...    return dot(grad(u), grad(v))

.. note::

   Using helpers such as :func:`~skfem.helpers.grad` and
   :func:`~skfem.helpers.dot` is optional.  Without helpers the last line would
   read, e.g., ``u.grad[0] * v.grad[0] + u.grad[1] * v.grad[1]``.  Inside the
   form ``u`` and ``v`` are of type :class:`~skfem.element.DiscreteField`.
   The return value is a numpy array.

Here is an example of body loading:

.. math::

   b(v) = \int_\Omega \sin(\pi x) \sin(\pi y) v \,\mathrm{d}x.

This can be written as

.. doctest::

   >>> import numpy as np
   >>> from skfem import LinearForm
   >>> @LinearForm
   ... def loading(v, w):
   ...    return np.sin(np.pi * w.x[0]) * np.sin(np.pi * w.x[1]) * v

.. note::

   The last argument ``w`` is a dictionary of
   :class:`~skfem.element.DiscreteField` objects.  Its ``_getattr_`` is
   overridden so that ``w.x`` corresponds to ``w['x']``.  Some keys are
   populated by default, e.g., ``w.x`` are the global quadrature points.

In addition, forms can depend on the local mesh parameter ``w.h`` or other
finite element functions (see :ref:`predefined`).  Moreover, boundary forms
assembled using :class:`~skfem.assembly.FacetBasis` can depend on the
outward normal vector ``w.n``.  One example is the form

.. math::

   l(\boldsymbol{v}) = \int_{\partial \Omega} \boldsymbol{v} \cdot \boldsymbol{n} \,\mathrm{d}s

which can be written as

.. doctest::

   >>> from skfem import LinearForm
   >>> from skfem.helpers import dot
   >>> @LinearForm
   ... def loading(v, w):
   ...    return dot(w.n, v)


The form definition always returns a two-dimensional numpy array.  This can be
verified using the Python debugger:

.. code-block:: python

   from skfem import *
   from skfem.helpers import grad, dot
   @BilinearForm
   def integrand(u, v, w):
       import pdb; pdb.set_trace()  # breakpoint
       return dot(grad(u), grad(v))

Saving the above snippet as ``test.py`` and running it via ``python test.py``
allows experimenting:

.. code-block:: none

   tom@tunkki:~/src/scikit-fem$ python -i test.py
   >>> asm(integrand, Basis(MeshTri(), ElementTriP1()))
   > /home/tom/src/scikit-fem/test.py(7)integrand()
   -> return dot(grad(u), grad(v))
   (Pdb) dot(grad(u), grad(v))
   array([[2., 2., 2.],
          [1., 1., 1.]])

Notice how ``dot(grad(u), grad(v))`` is a numpy array with the shape `number of
elements` x `number of quadrature points per element`.  The return value should
always have such shape no matter which mesh or element type is used.

The module :mod:`skfem.helpers` contains functions that make the forms more
readable.  Notice how the shape of ``u.grad[0]`` is what we expect also from
the return value:

.. code-block:: none

   tom@tunkki:~/src/scikit-fem$ python -i test.py
   >>> asm(integrand, Basis(MeshTri(), ElementTriP1()))
   > /home/tom/src/scikit-fem/test.py(7)integrand()
   -> return dot(grad(u), grad(v))
   (Pdb) !u.grad[0]
   array([[0.66666667, 0.16666667, 0.16666667],
          [0.66666667, 0.16666667, 0.16666667]])


.. _dofindexing:

Indexing of the degrees-of-freedom
==================================

.. warning::

   This section contains details on the order of the DOFs.
   Read this only if you did not find an answer in :ref:`finddofs`.

After finite element assembly we have the linear system

.. math::

   Ax = b.

What is the order of the unknowns in the vector :math:`x`?
The DOFs are ordered automatically based on the mesh and the element type.  It
is possible to investigate manually how the DOFs match the different
topological entities (`nodes`, `facets`, `edges`, `elements`) of the mesh.

.. note::

   **Nomenclature:** In scikit-fem, `edges` exist only for three-dimensional
   meshes so that `facets` are something always shared between two elements of
   the mesh.  In particular, we refer to the edges of triangular and
   quadrilateral meshes as `facets`.

For example, consider the triquadratic hexahedral element and the default
cube mesh:

.. doctest::

   >>> from skfem import *
   >>> m = MeshHex()
   >>> m
   <skfem MeshHex1 object>
     Number of elements: 1
     Number of vertices: 8
     Number of nodes: 8
   >>> basis = Basis(m, ElementHex2())
   >>> basis
   <skfem CellBasis(MeshHex1, ElementHex2) object>
     Number of elements: 1
     Number of DOFs: 27
     Size: 74088 B

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   draw(MeshHex())

The DOFs corresponding to the nodes (or vertices) of the mesh are

.. doctest::

   >>> basis.nodal_dofs
   array([[0, 1, 2, 3, 4, 5, 6, 7]])

This means that the first (zeroth) entry in the DOF array corresponds to the
first node/vertex in the finite element mesh (see ``m.p`` for a list of
nodes/vertices).

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   m = MeshHex()
   basis = Basis(m, ElementHex2())
   ax = draw(m)
   for dof in basis.nodal_dofs.flatten():
       ax.text(*basis.doflocs[:, dof], str(dof))

Similarly, the DOFs corresponding to the edges (``m.edges`` for a list of
edges) and the facets (``m.facets`` for a list of facets) of the mesh are

.. doctest::

   >>> basis.edge_dofs
   array([[ 8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]])
   >>> basis.facet_dofs
   array([[20, 21, 22, 23, 24, 25]])

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   m = MeshHex()
   basis = Basis(m, ElementHex2())
   ax = draw(m)
   for dof in basis.edge_dofs.flatten():
       ax.text(*basis.doflocs[:, dof], str(dof))

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   m = MeshHex()
   basis = Basis(m, ElementHex2())
   ax = draw(m)
   for dof in basis.facet_dofs.flatten():
       ax.text(*basis.doflocs[:, dof], str(dof))

All DOFs in ``nodal_dofs``, ``edge_dofs`` and ``facet_dofs``
are shared between neighbouring elements to preserve continuity.
The remaining DOFs are internal to the element and not shared:

.. doctest::

   >>> basis.interior_dofs
   array([[26]])
   
Each DOF is associated either with a node (``nodal_dofs``), a facet
(``facet_dofs``), an edge (``edge_dofs``), or an element (``interior_dofs``).
=====================
 Gallery of examples
=====================

This page contains an overview of the examples contained in the `source code
repository <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/>`_.

Poisson equation
================

Example 1: Poisson equation with unit load
------------------------------------------

This example solves the Poisson problem :math:`-\Delta u = 1` with the Dirichlet
boundary condition :math:`u = 0` in the unit square using piecewise-linear
triangular elements.

.. plot::
   :caption: The solution of Example 1.

   from docs.examples.ex01 import visualize
   visualize()

See the `source code of Example 1 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex01.py>`_ for more information.

.. _ex07:

Example 7: Discontinuous Galerkin method
----------------------------------------

This example solves the Poisson problem :math:`-\Delta u = 1` with :math:`u=0`
on the boundary using discontinuous Galerkin method.  The finite element basis
is piecewise-quartic but discontinuous over the element edges.

.. plot::
   :caption: The solution of Example 7.

   from docs.examples.ex07 import visualize
   visualize()

See the `source code of Example 7 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex07.py>`_ for more information.

Example 12: Postprocessing
--------------------------

This example demonstrates postprocessing the value of a functional, Boussinesq's k-factor.

.. plot::
   :caption: The solution of Example 12.

   from docs.examples.ex12 import visualize
   visualize()

See the `source code of Example 12 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex12.py>`_ for more information.

Example 13: Laplace with mixed boundary conditions
--------------------------------------------------

This example solves :math:`\Delta u = 0` in
:math:`\Omega=\{(x,y):1<x^2+y^2<4,~0<\theta<\pi/2\}`, where :math:`\tan \theta =
y/x`, with :math:`u = 0` on :math:`y = 0`, :math:`u = 1` on :math:`x =
0`, and :math:`\frac{\partial u}{\partial n} = 0` on the rest of the
boundary.

.. plot::
   :caption: The solution of Example 13.

   from docs.examples.ex13 import visualize
   visualize()

See the `source code of Example 13 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex13.py>`_ for more information.

.. _ex14:

Example 14: Laplace with inhomogeneous boundary conditions
----------------------------------------------------------

This example demonstrates how to impose coordinate-dependent Dirichlet
conditions for the Laplace equation :math:`\Delta u = 0`. The solution will
satisfy :math:`u=x^2 - y^2` on the boundary of the square domain.

.. plot::
   :caption: The solution of Example 14.

   from docs.examples.ex14 import visualize
   visualize()

See the `source code of Example 14 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex14.py>`_ for more information.

Example 15: One-dimensional Poisson equation
--------------------------------------------

This example solves :math:`-u'' = 1` in :math:`(0,1)` with the boundary
condition :math:`u(0)=u(1)=0`.

.. figure:: https://user-images.githubusercontent.com/973268/87775166-52b70b80-c82e-11ea-9009-c9fa0a9e28e8.png

   The solution of Example 15.

See the `source code of Example 15 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex15.py>`_ for more information.


Example 9: Three-dimensional Poisson equation
---------------------------------------------

This example solves :math:`-\Delta u = 1` with :math:`u=0` on the boundary
using linear tetrahedral elements and a preconditioned conjugate gradient
method.

.. note::

   This example will make use of the external packages `PyAMG
   <https://pypi.org/project/pyamg/>`__ or `pyamgcl
   <https://pypi.org/project/pyamgcl/>`__, if installed.

.. figure:: https://user-images.githubusercontent.com/973268/93183072-33abfb80-f743-11ea-9076-1324cbf28531.png

   The solution of Example 9 on a cross-section of the tetrahedral mesh.  The
   figure was created using `ParaView <https://www.paraview.org/>`__.

See the `source code of Example 9 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex09.py>`_ for more information.

Example 22: Adaptive Poisson equation
-------------------------------------

This example solves Example 1 adaptively in an L-shaped domain.
Using linear elements, the error indicators read :math:`\eta_K^2 = h_K^2 \|f\|_{0,K}^2` and :math:`\eta_E^2 = h_E \| [[\nabla u_h \cdot n ]] \|_{0,E}^2`   
for each element :math:`K` and
edge :math:`E`.

.. plot::
   :caption: The final solution of Example 22.

   from docs.examples.ex22 import visualize
   visualize()

See the `source code of Example 22 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex22.py>`_ for more information.

Example 37: Mixed Poisson equation
----------------------------------

This example solves the mixed formulation of the Poisson equation
using the lowest order Raviart-Thomas elements.

.. figure:: https://user-images.githubusercontent.com/973268/93132097-c2862d00-f6dd-11ea-97ad-40aaf2732ad1.png

   The piecewise constant solution field.
   The figure was created using `ParaView <https://www.paraview.org/>`__.

See the `source code of Example 37 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex37.py>`_ for more information.

Example 38: Point source
------------------------

Point sources require different assembly to other linear forms.

This example computes the Green's function for a disk; i.e. the solution of the
Dirichlet problem for the Poisson equation with the source term concentrated at
a single interior point :math:`\boldsymbol{s}`, :math:`-\Delta u = \delta
(\boldsymbol{x} - \boldsymbol{s})`.

.. plot::
   :caption: The scalar potential in the disk with point source at (0.3, 0.2).

   from docs.examples.ex38 import visualize
   visualize()

See the `source code of Example 38 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex38.py>`_
for more information.

Example 40: Hybridizable discontinuous Galerkin method
------------------------------------------------------

This examples solves the Poisson equation with unit load using a technique
where the finite element basis is first discontinous across element edges and
then the continuity is recovered with the help of Lagrange multipliers defined
on the mesh skeleton (i.e. a "skeleton mesh" consisting only of the edges of
the original mesh).

.. plot::
   :caption: The solution of Example 40 on the skeleton mesh.

   from docs.examples.ex40 import visualize
   visualize()

See the `source code of Example 40 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex40.py>`_
for more information.

Example 41: Mixed meshes
------------------------

This example solves the Poisson equation with unit load on a mesh consisting
of both triangles and quadrilaterals.  The support for mixed meshes is
preliminary and works only for elements with nodal or internal
degrees-of-freedom (sharing face and edge DOFs between mesh types is
work-in-progress).

.. plot::
   :caption: The solution of Example 41 on the mesh with both
             triangles and quadrilaterals.

   from docs.examples.ex41 import visualize
   visualize()

See the `source code of Example 41 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex41.py>`_
for more information.

Solid mechanics
===============

Example 2: Kirchhoff plate bending problem
------------------------------------------

This example solves the biharmonic Kirchhoff plate bending problem :math:`D
\Delta^2 u = f` in the unit square with a constant loading :math:`f`, bending
stiffness :math:`D` and a combination of clamped, simply supported and free
boundary conditions.

.. plot::
   :caption: The solution of Example 2.

   from docs.examples.ex02 import visualize
   visualize()

See the `source code of Example 2 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex02.py>`_ for more information.

Example 3: Linear elastic eigenvalue problem
--------------------------------------------

This example solves the linear elastic eigenvalue problem
:math:`\mathrm{div}\,\sigma(u)= \lambda u` with
the displacement fixed on the left boundary.

.. plot::
   :caption: The fifth eigenmode of Example 3.

   from docs.examples.ex03 import visualize
   visualize()

See the `source code of Example 3 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex03.py>`_ for more information.

Example 4: Linearized contact problem
-------------------------------------

This example solves a single interation of the contact problem
between two elastic bodies using the Nitsche's method.
Triangular and quadrilateral second-order elements are used
in the discretization of the two elastic bodies.

.. plot::
   :caption: The displaced meshes and the von Mises stress of Example 4.

   from docs.examples.ex04 import visualize
   visualize()

See the `source code of Example 4 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex04.py>`_ for more information.


Example 8: Argyris basis functions
----------------------------------

This example visualizes the :math:`C^1`-continuous fifth degree Argyris basis
functions on a simple triangular mesh.
This element can be used in the conforming discretization of biharmonic problems.

.. plot::
   :caption: The Argyris basis functions of Example 8 corresponding to the
             middle node and the edges connected to it.

   from docs.examples.ex08 import visualize
   visualize()

See the `source code of Example 8 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex08.py>`_ for more information.

Example 11: Three-dimensional linear elasticity
-----------------------------------------------

This example solves the three-dimensional linear elasticity equations
:math:`\mathrm{div}\,\sigma(u)=0` using trilinear hexahedral elements.
Dirichlet conditions are set on the opposing faces of a cube: one face remains
fixed and the other is displaced slightly outwards.

.. figure:: https://user-images.githubusercontent.com/973268/87685532-31054800-c78c-11ea-9b89-bc41dc0cb80c.png

   The displaced mesh of Example 11.  The figure was created using `ParaView
   <https://www.paraview.org/>`__.

See the `source code of Example 11 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex11.py>`_ for more information.

Example 21: Structural vibration
--------------------------------

This example demonstrates the solution of a three-dimensional vector-valued
eigenvalue problem by considering the vibration of an elastic structure.

.. figure:: https://user-images.githubusercontent.com/973268/147790554-4b768d43-25fa-49cd-ab19-b16a199a6459.png

   The first eigenmode of Example 21.

See the `source code of Example 21 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex21.py>`_ for more information.

Example 34: Euler-Bernoulli beam
--------------------------------

This example solves the Euler-Bernoulli beam equation
:math:`(EI u'')'' = 1`
with the boundary conditions
:math:`u(0)=u'(0) = 0` and using cubic Hermite elements.
The exact solution at :math:`x=1` is :math:`u(1)=1/8`.

.. figure:: https://user-images.githubusercontent.com/973268/87859267-749eb400-c93c-11ea-82cd-2d488fda39d4.png

   The solution of Example 34.

See the `source code of Example 34 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex34.py>`_ for more information.

Example 36: Nearly incompressible hyperelasticity
-------------------------------------------------

This example demonstrates the implementation of a two field mixed formulation
for nearly incompressible Neo-Hookean solids.

.. figure:: https://user-images.githubusercontent.com/22624037/91212007-4055aa80-e6d5-11ea-8572-f27986887331.png

   The displacement contour of Example 36.
   The figure was created using `ParaView <https://www.paraview.org/>`__.

See the `source code of Example 36 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex36.py>`_ for more information.


Example 43: Hyperelasticity
---------------------------

This example demonstrates Newton's method applied to the classical formulation
of a hyperelastic Neo-Hookean solid.

.. figure:: https://user-images.githubusercontent.com/973268/147790182-64f4abf4-3909-4ec0-89ac-2add304b133d.png

   The deformed mesh of Example 43.
   The figure was created using `vedo <https://github.com/marcomusy/vedo>`__.

See the `source code of Example 43 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex43.py>`_ for more information.

Fluid mechanics
===============

Example 18: Stokes equations
----------------------------

This example solves for the creeping flow problem in the primitive variables,
i.e. velocity and pressure instead of the stream-function.  These are governed
by the Stokes momentum :math:`- \nu\Delta\boldsymbol{u} + \rho^{-1}\nabla p = \boldsymbol{f}` and the continuity equation :math:`\nabla\cdot\boldsymbol{u} = 0`.

.. figure:: https://user-images.githubusercontent.com/1588947/93292002-d6d64100-f827-11ea-9a0a-c64d5d2979b7.png

   The streamlines of Example 18.

See the `source code of Example 18 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex18.py>`_ for more information.

Example 20: Creeping flow via stream-function
---------------------------------------------

This example solves the creeping flow problem via the stream-function
formulation.
The stream-function :math:`\psi` for two-dimensional creeping flow is
governed by the biharmonic equation :math:`\nu \Delta^2\psi = \mathrm{rot}\,\boldsymbol{f}` where :math:`\nu` is the kinematic viscosity (assumed constant),
:math:`\boldsymbol{f}` the volumetric body-force, and :math:`\mathrm{rot}\,\boldsymbol{f} =
\partial f_y/\partial x - \partial f_x/\partial y`.  The boundary
conditions at a wall are that :math:`\psi` is constant (the wall is
impermeable) and that the normal component of its gradient vanishes (no
slip)

.. figure:: https://user-images.githubusercontent.com/1588947/93291998-d50c7d80-f827-11ea-861b-f24ed27072d0.png

   The velocity field of Example 20.

See the `source code of Example 20 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex20.py>`_ for more information.

Example 24: Stokes flow with inhomogeneous boundary conditions
--------------------------------------------------------------

This example solves the Stokes flow over a backward-facing step
with a parabolic velocity profile at the inlet.

.. figure:: https://user-images.githubusercontent.com/973268/87858848-92b6e500-c939-11ea-81f9-cc51f254d19e.png

   The streamlines of Example 24.

See the `source code of Example 24 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex24.py>`_ for more information.

Example 27: Backward-facing step
--------------------------------

This example uses `pacopy 0.1.2 <https://pypi.org/project/pacopy/0.1.2>`__ to extend
the Stokes equations over a backward-facing step (Example 24) to finite Reynolds
number; this means defining a residual for the nonlinear problem and its
derivatives with respect to the solution and to the Reynolds number.

.. note::
   This example requires the external package `pacopy 0.1.2 <https://pypi.org/project/pacopy/0.1.2>`__.

.. figure:: https://user-images.githubusercontent.com/973268/87858972-97c86400-c93a-11ea-86e4-66f870b03e48.png

   The streamlines of Example 27 for :math:`\mathrm{Re}=750`.

See the `source code of Example 27 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex27.py>`_ for more information.

Example 29: Linear hydrodynamic stability
-----------------------------------------

The linear stability of one-dimensional solutions of the Navier-Stokes equations
is governed by the `Orr-Sommerfeld equation <https://en.wikipedia.org/wiki/Orr%E2%80%93Sommerfeld_equation>`_.  This is expressed in terms of the stream-function
:math:`\phi` of the perturbation, giving a two-point boundary value problem      
:math:`\alpha\phi(\pm 1) = \phi'(\pm 1) = 0`
for a complex fourth-order ordinary differential equation,

.. math::
   \left(\alpha^2-\frac{\mathrm d^2}{\mathrm dz^2}\right)^2\phi
   = (\mathrm j\alpha R)\left\{
     (c - U)\left(\alpha^2-\frac{\mathrm d^2}{\mathrm dz^2}\right)\phi
     - U''\phi,
   \right\}
   
where :math:`U(z)` is the base velocity profile, :math:`c` and :math:`\alpha`
are the wavespeed and wavenumber of the disturbance, and :math:`R` is the
Reynolds number.

.. figure:: https://user-images.githubusercontent.com/973268/87859022-e0801d00-c93a-11ea-978f-b1930627010b.png

   The results of Example 29.

See the `source code of Example 29 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex29.py>`_ for more information.

Example 30: Krylov-Uzawa method for the Stokes equation
-------------------------------------------------------

This example solves the Stokes equation iteratively in a square domain.

.. figure:: https://user-images.githubusercontent.com/973268/87859044-06a5bd00-c93b-11ea-84c2-9fbb9fc6e832.png

   The pressure field of Example 30.

See the `source code of Example 30 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex30.py>`_ for more information.

Example 32: Block diagonally preconditioned Stokes solver
---------------------------------------------------------

This example solves the Stokes problem in three dimensions, with an
algorithm that scales to reasonably fine meshes (a million tetrahedra in a few
minutes).

.. note::
   This examples requires an implementation of algebraic multigrid (either `pyamgcl <https://pypi.org/project/pyamgcl>`_ or `pyamg <https://pypi.org/project/pyamg/>`_).

.. figure:: https://user-images.githubusercontent.com/1588947/96520786-8a18d680-12bb-11eb-981a-c3388f2c8e35.png

   The velocity and pressure fields of Example 32, clipped in the plane of spanwise symmetry, *z* = 0.
   The figure was created using `ParaView <https://www.paraview.org/>`_ 5.8.1.

See the `source code of Example 32 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex32.py>`_ for more information.

Example 42: Periodic meshes
---------------------------

This example solves the advection equation on a periodic square mesh.

.. figure:: https://user-images.githubusercontent.com/973268/133767233-a5d78ec4-ffe7-4d49-bc93-9d9a0faae5a1.png

   The solution of Example 42 on a periodic mesh.

See the `source code of Example 42 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex42.py>`_
for more information.

Heat transfer
=============

Example 17: Insulated wire
--------------------------

This example solves the steady heat conduction
with generation in an insulated wire. In radial
coordinates, the governing equations read: find :math:`T`
satisfying :math:`\nabla \cdot (k_0 \nabla T) + A = 0,~0<r<a`,
and
:math:`\nabla \cdot (k_1 \nabla T) = 0,~a<r<b`,
with the boundary condition
:math:`k_1 \frac{\partial T}{\partial r} + h T = 0` on :math:`r=b`.

.. figure:: https://user-images.githubusercontent.com/973268/87775309-8db93f00-c82e-11ea-9015-add2226ad01e.png

   The solution of Example 17.

See the `source code of Example 17 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex17.py>`_ for more information.

Example 19: Heat equation
-------------------------

This example solves the heat equation :math:`\frac{\partial T}{\partial t} = \kappa\Delta T` in the domain :math:`|x|<w_0` and :math:`|y|<w_1` with the initial value :math:`T_0(x,y) = \cos\frac{\pi x}{2w_0}\cos\frac{\pi y}{2w_1}` using the generalized trapezoidal
rule ("theta method") and fast time-stepping by factorizing the evolution matrix once and for all.

.. figure:: https://user-images.githubusercontent.com/973268/87778846-7b420400-c834-11ea-8ff6-c439699b2802.gif

   The solution of Example 19.

See the `source code of Example 19 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex19.py>`_ for more information.

Example 25: Forced convection
-----------------------------

This example solves the plane Graetz problem with the governing
advection-diffusion equation :math:`\mathrm{Pe} \;u\frac{\partial T}{\partial x}
= \nabla^2 T` where the velocity profile is :math:`u (y) = 6 y (1 - y)` and the
Péclet number :math:`\mathrm{Pe}` is the mean velocity times the width divided
by the thermal diffusivity.

.. figure:: https://user-images.githubusercontent.com/973268/87858907-f8a36c80-c939-11ea-87a2-7357d5f073b1.png

   The solution of Example 25.

See the `source code of Example 25 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex25.py>`_ for more information.

Example 26: Restricting problem to a subdomain
----------------------------------------------

This example extends Example 17 by restricting the solution to a subdomain.

.. figure:: https://user-images.githubusercontent.com/973268/87858933-3902ea80-c93a-11ea-9d54-464235ab6325.png

   The solution of Example 26.

See the `source code of Example 26 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex26.py>`_ for more information.

Example 28: Conjugate heat transfer
-----------------------------------

This example extends Example 25 to conjugate heat transfer by giving a finite
thickness and thermal conductivity to one of the walls.  The example is modified
to a configuration for which there exists a fully developed solution which can be
found in closed form: given a uniform heat flux over each of the walls, the
temperature field asymptotically is the superposition of a uniform longitudinal
gradient and a transverse profile.

.. note::
   This example requires the external package
   `pygmsh <https://pypi.org/project/pygmsh/>`__.

.. figure:: https://user-images.githubusercontent.com/973268/142778186-99d8e02e-d02e-4b54-ac09-53bda0591dac.png

   A comparison of inlet and outlet temperature profiles in Example 28.

See the `source code of Example 28 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex28.py>`_ for more information.

Example 39: One-dimensional heat equation
-----------------------------------------

This examples reduces the two-dimensional heat equation of Example 19 to
demonstrate the special post-processing required.

.. figure:: https://user-images.githubusercontent.com/1588947/127958860-6454e542-67ba-4e94-8053-5175da201daa.gif

   The solution of Example 39.

See the `source code of Example 39 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex39.py>`_
for more information.

Electromagnetism


Miscellaneous
=============

Example 10: Nonlinear minimal surface problem
---------------------------------------------

This example solves the nonlinear minimal surface problem :math:`\nabla \cdot
\left(\frac{1}{\sqrt{1 + \|u\|^2}} \nabla u \right)= 0` with :math:`u=g`
prescribed on the boundary of the square domain.  The nonlinear problem is
linearized using the Newton's method with an analytical Jacobian calculated by
hand.

.. figure:: https://user-images.githubusercontent.com/973268/87663902-1c658780-c76d-11ea-9e00-324a18769ad2.png

   The solution of Example 10.

See the `source code of Example 10 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex10.py>`_ for more information.

Example 16: Legendre's equation
-------------------------------

This example solves the eigenvalue problem :math:`((1 - x^2) u')' + k u = 0` in
:math:`(-1,1)`.

.. figure:: https://user-images.githubusercontent.com/973268/87775206-65c9db80-c82e-11ea-8c49-bf191915602a.png

   The six first eigenmodes of Example 16.

See the `source code of Example 16 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex16.py>`_ for more information.

Example 23: Bratu-Gelfand
-------------------------

This example solves the Bratu-Gelfand two-point boundary value problem :math:`u'' + \lambda \mathrm e^u = 0`, :math:`0 < x < 1`,
with :math:`u(0)=u(1)=0` and where :math:`\lambda > 0` is a parameter.

.. note::
   This example requires the external package `pacopy 0.1.2 <https://pypi.org/project/pacopy/0.1.2>`__.

.. figure:: https://user-images.githubusercontent.com/973268/87779278-38ccf700-c835-11ea-955a-b77a0336b791.png

   The results of Example 23.

See the `source code of Example 23 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex23.py>`_ for more information.

Example 31: Curved elements
---------------------------

This example solves the eigenvalue problem :math:`-\Delta u = \lambda u`
with the boundary condition :math:`u|_{\partial \Omega} = 0` using isoparametric
mapping via biquadratic basis and finite element approximation using fifth-order
quadrilaterals.

.. figure:: https://user-images.githubusercontent.com/973268/87859068-32c13e00-c93b-11ea-984d-684e1e4c5066.png

   An eigenmode of Example 31 in a curved mesh.

See the `source code of Example 31 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex31.py>`_ for more information.

Example 33: H(curl) conforming model problem
--------------------------------------------

This example solves the vector-valued problem :math:`\nabla \times \nabla \times
E + E = f` in domain :math:`\Omega = [-1, 1]^3` with the boundary condition
:math:`E \times n|_{\partial \Omega} = 0` using the lowest order Nédélec edge
element.

.. figure:: https://user-images.githubusercontent.com/973268/87859239-47520600-c93c-11ea-8241-d62fdfd2a9a2.png

   The solution of Example 33 with the colors given by the magnitude
   of the vector field.
   The figure was created using `ParaView <https://www.paraview.org/>`__.

See the `source code of Example 33 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex33.py>`_ for more information.

Example 35: Characteristic impedance and velocity factor
--------------------------------------------------------

This example solves the series inductance (per meter) and parallel capacitance
(per meter) of RG316 coaxial cable. These values are then used to compute the
characteristic impedance and velocity factor of the cable.

.. figure:: https://user-images.githubusercontent.com/973268/87859275-85e7c080-c93c-11ea-9e62-3a9a8ee86070.png

   The results of Example 35.

See the `source code of Example 35 <https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex35.py>`_ for more information.
=============================
 Documentation of scikit-fem
=============================

`scikit-fem <https://github.com/kinnala/scikit-fem>`_ is a pure
Python 3.7+ library for performing `finite element assembly
<https://en.wikipedia.org/wiki/Finite_element_method>`_. Its main purpose is
the transformation of bilinear forms into sparse matrices and linear forms into
vectors.  The library supports triangular, quadrilateral, tetrahedral and
hexahedral meshes as well as one-dimensional problems.

.. note::

    Installing the library is as simple as running

    .. code-block:: bash

        pip install scikit-fem[all]

    Remove ``[all]`` to not install the optional dependencies ``meshio`` and
    ``matplotlib``.

Table of contents
=================

.. toctree::

   self
   gettingstarted
   howto
   advanced
   listofexamples
   api
=============
How-to guides
=============

This section contains goal-oriented guides on the features of scikit-fem.

.. _finddofs:

Finding degrees-of-freedom
==========================

Often the goal is to constrain DOFs on a specific part of
the boundary.  Currently the main tool for finding DOFs is
:meth:`~skfem.assembly.basis.AbstractBasis.get_dofs`.

.. doctest::

   >>> from skfem import MeshTri, Basis, ElementTriP2
   >>> m = MeshTri().refined(2)
   >>> basis = Basis(m, ElementTriP2())

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   m = MeshTri().refined(2)
   basis = Basis(m, ElementTriP2())
   ax = draw(m)
   for dof in basis.nodal_dofs.flatten():
       ax.text(*basis.doflocs[:, dof], str(dof))

We can provide an indicator function to
:meth:`~skfem.assembly.basis.AbstractBasis.get_dofs` and it will find the
DOFs on the matching facets:

.. doctest::

   >>> dofs = basis.get_dofs(lambda x: x[0] == 0.)
   >>> dofs.nodal
   {'u': array([ 0,  2,  5, 10, 14])}
   >>> dofs.facet
   {'u': array([26, 30, 39, 40])}

This element has one DOF per node and one DOF per facet.  The facets have their
own numbering scheme starting from zero, however, the corresponding DOFs are
offset by the total number of nodal DOFs:

.. doctest::

   >>> dofs.facet['u']
   array([26, 30, 39, 40])

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   m = MeshTri().refined(2)
   basis = Basis(m, ElementTriP2())
   ax = draw(m)
   for dof in basis.facet_dofs.flatten():
       ax.text(*basis.doflocs[:, dof], str(dof))

The keys in the above dictionaries indicate the type of the DOF according to
the following table:

+-----------+---------------------------------------------------------------+
| Key       | Description                                                   |
+===========+===============================================================+
| ``u``     | Point value                                                   |
+-----------+---------------------------------------------------------------+
| ``u_n``   | Normal derivative                                             |
+-----------+---------------------------------------------------------------+
| ``u_x``   | Partial derivative w.r.t. :math:`x`                           |
+-----------+---------------------------------------------------------------+
| ``u_xx``  | Second partial derivative w.r.t :math:`x`                     |
+-----------+---------------------------------------------------------------+
| ``u^n``   | Normal component of a vector field (e.g., Raviart-Thomas)     |
+-----------+---------------------------------------------------------------+
| ``u^t``   | Tangential component of a vector field (e.g., Nédélec)        |
+-----------+---------------------------------------------------------------+
| ``u^1``   | First component of a vector field                             |
+-----------+---------------------------------------------------------------+
| ``u^1_x`` | Partial derivative of the first component w.r.t. :math:`x`    |
+-----------+---------------------------------------------------------------+
| ``u^1^1`` | First component of the first component in a composite field   |
+-----------+---------------------------------------------------------------+
| ``NA``    | Description not available (e.g., hierarchical or bubble DOF's)|
+-----------+---------------------------------------------------------------+

An array of all DOFs with the key ``u`` can be obtained as follows:

.. doctest::

   >>> dofs.all(['u'])
   array([ 0,  2,  5, 10, 14, 26, 30, 39, 40])
   >>> dofs.flatten()  # all DOFs, no matter which key
   array([ 0,  2,  5, 10, 14, 26, 30, 39, 40])

If a set of facets is tagged, the name of the tag can be passed
to :meth:`~skfem.assembly.basis.AbstractBasis.get_dofs`:

.. doctest::

   >>> dofs = basis.get_dofs('left')
   >>> dofs.flatten()
   array([ 0,  2,  5, 10, 14, 26, 30, 39, 40])
   
Many DOF types have a well-defined location.  These DOF locations can be found
as follows:

.. doctest::

   >>> basis.doflocs[:, dofs.flatten()]
   array([[0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ],
          [0.   , 1.   , 0.5  , 0.25 , 0.75 , 0.125, 0.875, 0.375, 0.625]])

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   m = MeshTri().refined(2)
   basis = Basis(m, ElementTriP2())
   dofs = basis.get_dofs('left')
   ax = draw(m)
   for dof in dofs.flatten():
       ax.plot(*basis.doflocs[:, dof], 'ro')
       ax.text(*basis.doflocs[:, dof], str(dof))

See :ref:`dofindexing` for more details.

.. _l2proj:

Performing projections
======================

We can use :math:`L^2` projection to find discrete counterparts of functions or
transform from one finite element basis to another.  Suppose we have
:math:`u_0(x,y) = x^3 y^3` defined on the boundary of the domain and want to
find the corresponding discrete function which is extended by zero in the
interior of the domain.  You could explicitly assemble and solve the linear
system corresponding to: find :math:`\widetilde{u_0} \in V_h` satisfying

.. math::

   \int_{\partial \Omega} \widetilde{u_0} v\,\mathrm{d}s = \int_{\partial \Omega} u_0 v\,\mathrm{d}s\quad \forall v \in V_h.

However, this is so common that we have a shortcut
:meth:`~skfem.assembly.AbstractBasis.project`:

.. doctest::

   >>> import numpy as np
   >>> from skfem import *
   >>> m = MeshQuad().refined(2)
   >>> basis = FacetBasis(m, ElementQuad1())
   >>> u0 = lambda x: x[0] ** 3 * x[1] ** 3
   >>> u0t = basis.project(u0)
   >>> np.abs(np.round(u0t, 5))
   array([1.0000e-05, 8.9000e-04, 9.7054e-01, 8.9000e-04, 6.0000e-05,
          6.0000e-05, 1.0944e-01, 1.0944e-01, 0.0000e+00, 2.0000e-05,
          2.0000e-05, 2.4000e-04, 8.0200e-03, 3.9797e-01, 3.9797e-01,
          2.4000e-04, 8.0200e-03, 0.0000e+00, 0.0000e+00, 0.0000e+00,
          0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00])

.. plot::

   import skfem as fem
   m = fem.MeshQuad().refined(2)
   basis = fem.FacetBasis(m, fem.ElementQuad1())
   u0 = lambda x: x[0] ** 3 * x[1] ** 3
   u0t = basis.project(u0)
   ibasis = fem.InteriorBasis(m, fem.ElementQuad1())
   from skfem.visuals.matplotlib import plot, draw
   ax = draw(ibasis)
   plot(ibasis, u0t, nrefs=3, ax=ax, colorbar=True, shading='gouraud')

We can also project over the entire domain:

.. doctest::

   >>> basis = Basis(m, ElementQuad1())
   >>> f = lambda x: np.sin(2. * np.pi * x[0]) + 1.
   >>> fh = basis.project(f)
   >>> np.abs(np.round(fh, 5))
   array([1.09848, 0.90152, 0.90152, 1.09848, 1.     , 1.09848, 0.90152,
          1.     , 1.     , 2.19118, 1.09848, 0.19118, 0.90152, 0.90152,
          0.19118, 1.09848, 2.19118, 1.     , 2.19118, 0.19118, 1.     ,
          2.19118, 0.19118, 0.19118, 2.19118])

.. plot::

   import skfem as fem
   m = fem.MeshQuad().refined(2)
   basis = fem.CellBasis(m, fem.ElementQuad1())
   f = lambda x: np.sin(2. * np.pi * x[0]) + 1.
   fh = basis.project(f)
   from skfem.visuals.matplotlib import plot, draw
   ax = draw(basis)
   plot(basis, fh, nrefs=3, ax=ax, colorbar=True, shading='gouraud')

We can project from one finite element basis to another:

.. doctest::

   >>> basis0 = basis.with_element(ElementQuad0())
   >>> fh = basis0.project(basis.interpolate(fh))
   >>> np.abs(np.round(fh, 5))
   array([1.64483, 0.40441, 0.40441, 1.64483, 1.59559, 0.35517, 0.35517,
          1.59559, 1.59559, 0.35517, 0.35517, 1.59559, 1.64483, 0.40441,
          0.40441, 1.64483])

.. plot::

   from skfem import *
   m = MeshQuad().refined(2)
   basis = CellBasis(m, ElementQuad1())
   basis0 = basis.with_element(ElementQuad0())
   f = lambda x: np.sin(2. * np.pi * x[0]) + 1.
   fh = basis.project(f)
   fh = basis0.project(basis.interpolate(fh))
   from skfem.visuals.matplotlib import plot, draw
   ax = draw(basis)
   plot(basis0, fh, nrefs=3, ax=ax, colorbar=True, shading='gouraud')

We can interpolate the gradient at quadrature points and project:

.. doctest::

   >>> f = lambda x: np.sin(2. * np.pi * x[0]) + 1.
   >>> fh = basis.project(f)  # P1
   >>> fh = basis.project(basis.interpolate(fh).grad[0])  # df/dx
   >>> np.abs(np.round(fh, 5))
   array([6.6547 , 6.6547 , 6.6547 , 6.6547 , 7.04862, 6.6547 , 6.6547 ,
          7.04862, 7.04862, 0.19696, 6.6547 , 0.19696, 6.6547 , 6.6547 ,
          0.19696, 6.6547 , 0.19696, 7.04862, 0.19696, 0.19696, 7.04862,
          0.19696, 0.19696, 0.19696, 0.19696])

.. plot::

   from skfem import *
   m = MeshQuad().refined(2)
   basis = CellBasis(m, ElementQuad1())
   basis0 = basis.with_element(ElementQuad0())
   f = lambda x: np.sin(2. * np.pi * x[0]) + 1.
   fh = basis.project(f)
   fh = basis.project(basis.interpolate(fh).grad[0])
   from skfem.visuals.matplotlib import plot, draw
   ax = draw(basis)
   plot(basis, fh, nrefs=3, ax=ax, colorbar=True, shading='gouraud')

.. _predefined:

Discrete functions in forms
===========================

We can use finite element functions inside the form by interpolating them at
quadrature points.  For example, consider a fixed-point iteration for the
nonlinear problem

.. math::

   \begin{aligned}
      -\nabla \cdot ((u + \tfrac{1}{10})\nabla u) &= 1 \quad \text{in $\Omega$}, \\
      u &= 0 \quad \text{on $\partial \Omega$}.
   \end{aligned}

We repeatedly find :math:`u_{k+1} \in H^1_0(\Omega)` which satisfies

.. math::

   \int_\Omega (u_{k} + \tfrac{1}{10}) \nabla u_{k+1} \cdot \nabla v \,\mathrm{d}x = \int_\Omega v\,\mathrm{d}x

for every :math:`v \in H^1_0(\Omega)`.
The bilinear form depends on the previous solution :math:`u_k`.

.. doctest::

   >>> import skfem as fem
   >>> from skfem.models.poisson import unit_load
   >>> from skfem.helpers import grad, dot
   >>> @fem.BilinearForm
   ... def bilinf(u, v, w):
   ...     return (w.u_k + .1) * dot(grad(u), grad(v))

The previous solution :math:`u_k` is interpolated at quadrature points using
:meth:`~skfem.assembly.CellBasis.interpolate` and then provided to
:meth:`~skfem.assembly.BilinearForm.assemble` as a keyword argument:

.. doctest::

   >>> m = fem.MeshTri().refined(3)
   >>> basis = fem.Basis(m, fem.ElementTriP1())
   >>> b = unit_load.assemble(basis)
   >>> x = 0. * b.copy()
   >>> for itr in range(20):  # fixed point iteration
   ...     A = bilinf.assemble(basis, u_k=basis.interpolate(x))
   ...     x = fem.solve(*fem.condense(A, b, I=m.interior_nodes()))
   ...     print(round(x.max(), 10))
   0.7278262868
   0.1956340215
   0.3527261363
   0.2745541843
   0.3065381711
   0.2921831118
   0.298384264
   0.2956587119
   0.2968478347
   0.2963273314
   0.2965548428
   0.2964553357
   0.2964988455
   0.2964798184
   0.2964881386
   0.2964845003
   0.2964860913
   0.2964853955
   0.2964856998
   0.2964855667

.. plot::

   import skfem as fem
   from skfem.models.poisson import unit_load
   from skfem.helpers import grad, dot
   @fem.BilinearForm
   def bilinf(u, v, w):
       return (w.u_k + .1) * dot(grad(u), grad(v))
   m = fem.MeshTri().refined(4)
   basis = fem.Basis(m, fem.ElementTriP1())
   b = unit_load.assemble(basis)
   x = 0. * b.copy()
   for itr in range(20):  # fixed point iteration
       A = bilinf.assemble(basis, u_k=basis.interpolate(x))
       x = fem.solve(*fem.condense(A, b, I=m.interior_nodes()))
   from skfem.visuals.matplotlib import *
   plot(basis, x, colorbar=True, nrefs=3, shading='gouraud')

.. note::

    Inside the form definition, ``w`` is a dictionary of user provided
    arguments and additional default keys.  By default, ``w['x']`` (accessible
    also as ``w.x``) corresponds to the global coordinates and ``w['h']``
    (accessible also as ``w.h``) corresponds to the local mesh parameter.

Assembling jump terms
=====================

The shorthand :func:`~skfem.assembly.asm`
supports special syntax for assembling the same form over lists of
bases and summing the result.  The form

.. math::

   b(u,v) = \sum_{E \in \mathcal{E}_h} \int_{E} [u][v]\,\mathrm{d}s

with jumps
:math:`[u] = u_1 - u_2` and :math:`[v] = v_1 - v_2`
over the interior edges can be split as

.. math::

   b(u,v) = \sum_{E \in \mathcal{E}_h} \left(\int_{E} u_1 v_1\,\mathrm{d}s - \int_{E} u_1 v_2\,\mathrm{d}s - \int_{E} u_2 v_1\,\mathrm{d}s + \int_{E} u_2 v_2\,\mathrm{d}s\right)

and normally we would assemble all of the four forms separately.

We can instead provide a list of bases during a call to :func:`skfem.assembly.asm`:

.. doctest::

   >>> import skfem as fem
   >>> m = fem.MeshTri()
   >>> e = fem.ElementTriP0()
   >>> bases = [fem.InteriorFacetBasis(m, e, side=k) for k in [0, 1]]
   >>> jumpform = fem.BilinearForm(lambda u, v, p: (-1) ** sum(p.idx) * u * v)
   >>> fem.asm(jumpform, bases, bases).toarray()
   array([[ 1.41421356, -1.41421356],
          [-1.41421356,  1.41421356]])

For an example of practical usage, see :ref:`ex07`.
.. _gettingstarted:

=================
 Getting started
=================

If you have a supported Python installation on your computer, you can
install the package via

.. code-block:: bash

   pip install scikit-fem[all]

Specifying ``[all]`` includes ``meshio`` for mesh input/output, and
``matplotlib`` for simple visualizations.  The minimal dependencies are
``numpy`` and ``scipy``.  You can also install scikit-fem in `Google Colab
<https://colab.research.google.com/>`_ by executing

.. code-block:: bash

   !pip install scikit-fem

Step 1: Clarify the problem
===========================

In this tutorial we solve the Poisson problem

.. math::
   \begin{aligned}
        -\Delta u &= f \quad && \text{in $\Omega$,} \\
        u &= 0 \quad && \text{on $\partial \Omega$,}
   \end{aligned}

where :math:`\Omega = (0, 1)^2` is a square domain
and :math:`f(x,y)=\sin \pi x \sin \pi y`.
The weak formulation reads:
find :math:`u \in H^1_0(\Omega)` satisfying

.. math::
   \int_\Omega \nabla u \cdot \nabla v \,\mathrm{d}x = \int_\Omega fv\,\mathrm{d}x \quad \forall v \in H^1_0(\Omega).

.. note::

   :math:`H^1_0(\Omega)` is the space of functions that are zero on the
   boundary :math:`\partial \Omega` and have finite energy: the square integral
   of the first derivative is finite.

Step 2: Express the forms as code
=================================

Next we write the forms

.. math::

   a(u, v) = \int_\Omega \nabla u \cdot \nabla v \,\mathrm{d}x \quad \text{and} \quad L(v) = \int_\Omega f v \,\mathrm{d}x

as source code.  Each form is written as a function and
decorated as follows:

.. doctest::

   >>> import skfem as fem
   >>> from skfem.helpers import dot, grad  # helpers make forms look nice
   >>> @fem.BilinearForm
   ... def a(u, v, _):
   ...     return dot(grad(u), grad(v))

.. doctest::

   >>> import numpy as np
   >>> @fem.LinearForm
   ... def L(v, w):
   ...     x, y = w.x  # global coordinates
   ...     f = np.sin(np.pi * x) * np.sin(np.pi * y)
   ...     return f * v

For more information see :ref:`forms`.

Step 3: Create a mesh
=====================

The default constructors of :class:`~skfem.mesh.Mesh` initialize a
unit square:

.. doctest::

   >>> mesh = fem.MeshTri().refined(3)  # refine thrice
   >>> mesh
   <skfem MeshTri1 object>
     Number of elements: 128
     Number of vertices: 81
     Number of nodes: 81
     Named boundaries [# facets]: left [8], bottom [8], right [8], top [8]


.. plot::

   from skfem import *
   MeshTri().refined(3).draw(boundaries=True)


Step 4: Define a basis
======================

The mesh is combined with a finite element to form a global
basis.
Here we choose the piecewise-linear basis:

.. doctest::

   >>> Vh = fem.Basis(mesh, fem.ElementTriP1())
   >>> Vh
   <skfem CellBasis(MeshTri1, ElementTriP1) object>
     Number of elements: 128
     Number of DOFs: 81
     Size: 9216 B

Step 5: Assemble the linear system
==================================

Now everything is in place for the finite element assembly.
The resulting matrix has the type ``scipy.sparse.csr_matrix``
and the load vector has the type ``ndarray``.

.. doctest::

   >>> A = a.assemble(Vh)
   >>> l = L.assemble(Vh)
   >>> A.shape
   (81, 81)
   >>> l.shape
   (81,)

Step 6: Find boundary DOFs
==========================

Setting boundary conditions requires finding the degrees-of-freedom (DOFs) on
the boundary.  Empty call to
:meth:`~skfem.assembly.basis.AbstractBasis.get_dofs` matches all boundary DOFs.

.. doctest::

   >>> D = Vh.get_dofs()
   >>> D
   <skfem DofsView(MeshTri1, ElementTriP1) object>
     Number of nodal DOFs: 32 ['u']

Step 7: Eliminate boundary DOFs and solve
=========================================

The boundary DOFs must be eliminated from the linear system :math:`Ax=l`
to set :math:`u=0` on the boundary.
This can be done using :func:`~skfem.utils.condense`.
The output can be passed to :func:`~skfem.utils.solve`
which is a simple wrapper to ``scipy`` sparse solver:

.. doctest::

   >>> x = fem.solve(*fem.condense(A, l, D=D))
   >>> x.shape
   (81,)

.. plot::

   from skfem import *
   from skfem.visuals.matplotlib import *
   from skfem.helpers import dot, grad
   import numpy as np
   basis = Basis(MeshTri().refined(3), ElementTriP1())
   a = BilinearForm(lambda u, v, _: dot(grad(u), grad(v)))
   L = LinearForm(lambda v, w: np.sin(np.pi * w.x[0]) * np.sin(np.pi * w.x[1]) * v)
   y = solve(*condense(a.assemble(basis), L.assemble(basis), D=basis.get_dofs()))
   ax = draw(basis)
   plot(basis, y, ax=ax, nrefs=2, colorbar=True, shading='gouraud')


Step 8: Calculate error
=======================

The exact solution is known to be

.. math::

   u(x, y) = \frac{1}{2 \pi^2} \sin \pi x \sin \pi y.

Thus, it makes sense to verify that the error is small.

.. doctest::

   >>> @fem.Functional
   ... def error(w):
   ...     x, y = w.x
   ...     uh = w['uh']
   ...     u = np.sin(np.pi * x) * np.sin(np.pi * y) / (2. * np.pi ** 2)
   ...     return (uh - u) ** 2
   >>> round(error.assemble(Vh, uh=Vh.interpolate(x)), 9)
   1.069e-06
==========================
 Detailed API description
==========================

This section contains API documentation for the most commonly used interfaces
of the library.

Module: skfem.mesh
==================

.. automodule:: skfem.mesh

Abstract class: Mesh
--------------------

.. autoclass:: skfem.mesh.Mesh
   :members: load, save, refined, facets_satisfying, nodes_satisfying, elements_satisfying

Class: MeshTri
**************

.. autoclass:: skfem.mesh.MeshTri

.. autoclass:: skfem.mesh.MeshTri1
   :members: __init__, init_symmetric, init_sqsymmetric, init_refdom, init_tensor, init_lshaped, init_circle, load

.. autoclass:: skfem.mesh.MeshTri2
   :members: init_circle, load

Class: MeshQuad
***************

.. autoclass:: skfem.mesh.MeshQuad

.. autoclass:: skfem.mesh.MeshQuad1
   :members: __init__, init_refdom, init_tensor, to_meshtri, load

.. autoclass:: skfem.mesh.MeshQuad2
   :members: load

Class: MeshTet
**************

.. autoclass:: skfem.mesh.MeshTet

.. autoclass:: skfem.mesh.MeshTet1
   :members: __init__, init_refdom, init_tensor, init_ball, load

.. autoclass:: skfem.mesh.MeshTet2
   :members: init_ball, load

Class: MeshHex
**************

.. autoclass:: skfem.mesh.MeshHex

.. autoclass:: skfem.mesh.MeshHex1
   :members: __init__, init_tensor, to_meshtet, load

Class: MeshLine
***************

.. autoclass:: skfem.mesh.MeshLine

.. autoclass:: skfem.mesh.MeshLine1
   :members: __init__

Module: skfem.assembly
======================

.. automodule:: skfem.assembly

.. autofunction:: skfem.assembly.asm

Abstract class: AbstractBasis
-----------------------------

Subclasses of :class:`~skfem.assembly.basis.AbstractBasis` represent a global
finite element basis evaluated at quadrature points.

.. autoclass:: skfem.assembly.basis.AbstractBasis
   :members: get_dofs

Class: CellBasis
****************

.. autoclass:: skfem.assembly.Basis

.. autoclass:: skfem.assembly.CellBasis
   :members: __init__, interpolate, project


Class: FacetBasis
*****************

.. autoclass:: skfem.assembly.BoundaryFacetBasis

.. autoclass:: skfem.assembly.FacetBasis
   :members: __init__

Class: InteriorFacetBasis
*************************

.. autoclass:: skfem.assembly.InteriorFacetBasis
   :members: __init__


Abstract class: Form
--------------------

Class: BilinearForm
*******************

.. autoclass:: skfem.assembly.BilinearForm
   :members: assemble

Class: LinearForm
*****************

.. autoclass:: skfem.assembly.LinearForm
   :members: assemble

Class: Functional
*****************

.. autoclass:: skfem.assembly.Functional
   :members: assemble, elemental

Module: skfem.element
=====================

.. automodule:: skfem.element
   :show-inheritance:

.. autosummary::

    skfem.element.ElementH1
    skfem.element.ElementVector
    skfem.element.ElementHdiv
    skfem.element.ElementHcurl
    skfem.element.ElementGlobal
    skfem.element.ElementDG
    skfem.element.ElementComposite
    skfem.element.ElementTriP1
    skfem.element.ElementTriP2
    skfem.element.ElementTriP3
    skfem.element.ElementTriP4
    skfem.element.ElementTriP0
    skfem.element.ElementTriCR
    skfem.element.ElementTriCCR
    skfem.element.ElementTriRT0
    skfem.element.ElementTriMorley
    skfem.element.ElementTri15ParamPlate
    skfem.element.ElementTriArgyris
    skfem.element.ElementTriMini
    skfem.element.ElementTriHermite
    skfem.element.ElementTriSkeletonP0
    skfem.element.ElementTriSkeletonP1
    skfem.element.ElementTriBDM1
    skfem.element.ElementQuad0
    skfem.element.ElementQuad1
    skfem.element.ElementQuad2
    skfem.element.ElementQuadS2
    skfem.element.ElementQuadP
    skfem.element.ElementQuadBFS
    skfem.element.ElementQuadRT0
    skfem.element.ElementTetP0
    skfem.element.ElementTetP1
    skfem.element.ElementTetP2
    skfem.element.ElementTetRT0
    skfem.element.ElementTetN0
    skfem.element.ElementTetMini
    skfem.element.ElementTetCR
    skfem.element.ElementTetCCR
    skfem.element.ElementHex0
    skfem.element.ElementHex1
    skfem.element.ElementHex2
    skfem.element.ElementHexS2
    skfem.element.ElementLineP0
    skfem.element.ElementLineP1
    skfem.element.ElementLineP2
    skfem.element.ElementLinePp
    skfem.element.ElementLineHermite
    skfem.element.ElementLineMini
   

.. note::

   The element global basis is calculated at quadrature points and stored inside
   :class:`~skfem.element.DiscreteField` objects.
   The different element types precalculate different fields of
   :class:`~skfem.element.DiscreteField`.  E.g., for :math:`H(div)`
   finite elements it is natural to precalculate ``DiscreteField.div``.
   The high order derivatives are created only when using subclasses of
   :class:`~skfem.element.ElementGlobal`.
                 

.. autoclass:: skfem.element.DiscreteField

Module: skfem.utils
===================

Function: solve
---------------

.. autofunction:: skfem.utils.solve

Function: condense
------------------

.. autofunction:: skfem.utils.condense

Function: enforce
-----------------

.. autofunction:: skfem.utils.enforce

Module: skfem.helpers
=====================

.. automodule:: skfem.helpers

.. autofunction:: skfem.helpers.grad

.. autofunction:: skfem.helpers.div

.. autofunction:: skfem.helpers.curl

.. autofunction:: skfem.helpers.d

.. autofunction:: skfem.helpers.dd

.. autofunction:: skfem.helpers.ddd

.. autofunction:: skfem.helpers.dddd

.. autofunction:: skfem.helpers.sym_grad

.. autofunction:: skfem.helpers.dot

.. autofunction:: skfem.helpers.ddot

.. autofunction:: skfem.helpers.dddot

.. autofunction:: skfem.helpers.mul

.. autofunction:: skfem.helpers.trace

.. autofunction:: skfem.helpers.transpose

.. autofunction:: skfem.helpers.prod

.. autofunction:: skfem.helpers.inv
