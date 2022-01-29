---
title: 'A Fast Iterative Method Python package'
tags:
  - Python
  - eikonal
  - partial differential equations
  - cuda
authors:
  - name: Thomas Grandits
    affiliation: 1
affiliations:
 - name: Institute of Computer Graphics and Vision, TU Graz
   index: 1
date: July 2021
bibliography: paper.bib
---

# Summary


The anisotropic eikonal equation is a non-linear partial differential equation, given by
\begin{equation*}
   \left\{
   \begin{array}{rll}
   \left<\nabla \phi, D \nabla \phi \right> &= 1 \quad &\text{on} \; \Omega \\
   \phi(\mathbf{x}_0) &= g(\mathbf{x}_0) \quad &\text{on} \; \Gamma \subset \Omega
   \end{array}
   \right. .
\end{equation*}
In practice, this problem is often associated with computing the earliest arrival times $\phi$ of a wave from a set of given starting points $\mathbf{x}_0$ through a heterogeneous medium (i.e. different velocities are assigned throughout the medium). 
This equation yields infinitely many weak solutions [@evans_partial_2010] and can thus not be straight-forwardly solved using standard Finite Element approaches.

``fim-python`` implements the Fast Iterative Method (FIM), proposed in [@fu_fast_2013], purely in Python to solve the anisotropic eikonal equation by finding its unique viscosity solution.
In this scenario, we compute $\phi$ on tetrahedral/triangular meshes or line networks for a given $D$, $\mathbf{x}_0$ and $g$.
The method is implemented both on the CPU using [``numba``](https://numba.pydata.org/) and [``numpy``](https://numpy.org/), as well as the GPU with the help of [``cupy``](https://cupy.dev/) (depends on [CUDA](https://developer.nvidia.com/cuda-toolkit)).
The library is meant to be easily and rapidly used for repeated evaluations on a mesh.

The FIM locally computes an update rule to find the path the wavefront will take through a single element.
Since the algorithm is restricted to linear elements, the path through an element will also be a straight line.
In the case of tetrahedral domains, the FIM thus tries to find the path of the linear update from a face spanned by three vertices $\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$ to the opposite vertex $\mathbf{v}_4$.
\autoref{fig:update} visualizes the update.
For triangles and lines, the algorithm behaves similarly but the update origin is limited to a side or vertex respectively.
The exact equations used to solve this problem in this repository were previously described (among others) in [@grandits_inverse_2020].

![Update inside a single tetrahedron\label{fig:update}](docs/figs/update_fig.jpg "Update inside a single tetrahedron"){ width=33% }


Two different methods are implemented in ``fim-python``:
In the *Jacobi* method, the above local update rule is computed for all elements in each iteration until the change between two subsequent iterations is smaller than a chosen $\varepsilon$.
This version of the algorithm is bested suited for the GPU, since it is optimal for a SIMD (single instruction multiple data) architecture.
The *active list* method is more closely related to the method presented in [@fu_fast_2013]:
We keep track of all vertices that require a recomputation in the current iteration on a so-called active list which we keep up-to-date. 

# Comparison to other tools

There are other tools available to solve variants of the eikonal equation, but they differ in functionality to ``fim-python``.

[``scikit-fmm``](https://pypi.org/project/scikit-fmm/) implements the Fast Marching Method (FMM) [@sethian_fast_1996], which was designed to solve the isotropic eikonal equation ($D = c I$ for $c \in \mathbb{R}$ and $I$ being the identity matrix). The library works on uniform grids, rather than meshes.

[``GPUTUM: Unstructured Eikonal``](https://github.com/SCIInstitute/SCI-Solver_Eikonal) implements the FIM in CUDA for triangulated surfaces and tetrahedral meshes, but has no Python bindings and is designed as a command line tool for single evaluations.

# Statement of need

The eikonal equation has many practical applications, including cardiac electrophysiology, image processing and geoscience, to approximate wave propagation through a medium.
In the example of cardiac electrophysiology [@franzone2014mathematical], the electrical activation times $\phi$ are computed throughout the anisotropic heart muscle with varying conduction velocities $D$.

``fim-python`` tries to wrap the FIM for CPU and GPU into an easy-to-use Python package for multiple evaluations with a straight-forward installation over [PyPI](https://pypi.org/).
This should provide engineers and researchers alike with an accessible tool that allows evaluations of the eikonal equation for general scenarios. 

# References
# Fast Iterative Method - Numpy/Cupy
This repository implements the Fast Iterative Method on [tetrahedral domains](https://epubs.siam.org/doi/abs/10.1137/120881956) and [triangulated surfaces](https://epubs.siam.org/doi/abs/10.1137/100788951) purely in python both for CPU (numpy) and GPU (cupy). The main focus is however on the GPU implementation, since it can be better exploited for very large domains.

[![codecov](https://codecov.io/gh/thomgrand/fim-python/branch/master/graph/badge.svg?token=DG05WR5030)](https://codecov.io/gh/thomgrand/fim-python)
[![CI Tests](https://github.com/thomgrand/fim-python/actions/workflows/python-package.yml/badge.svg)](https://github.com/thomgrand/fim-python/actions/workflows/python-package.yml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03641/status.svg)](https://doi.org/10.21105/joss.03641)

# Details
The anisotropic eikonal equation is given by

![$$\left<D \nabla \phi, \nabla \phi\right> = 1$$](https://latex.codecogs.com/svg.latex?\Large&space;\left%3CD%20\nabla%20\phi,%20\nabla%20\phi\right%3E%20=%201)


for given boundary conditions 

![$$\phi(\mathbf{x}_0) = g(\mathbf{x}_0)$$](https://latex.codecogs.com/svg.latex?\Large\phi(\mathbf{x}_0)%20=%20g(\mathbf{x}_0))

For a given anisotropic velocity, this can calculate the geodesic distance between a set of ![$\mathbf{x}_0$](https://latex.codecogs.com/svg.latex?\Large\mathbf{x}_0) and all points on the domain like shown in the figure.

![Preview Image](docs/figs/usage_example.jpg)

Note that when using multiple ![$\mathbf{x}_0$](https://latex.codecogs.com/svg.latex?\Large\mathbf{x}_0), they are not guaranteed to be in the final solution if they are not a valid viscosity solution. A recommended read for more details on the subject is:  
Evans, Lawrence C. "Partial differential equations." *Graduate studies in mathematics* 19.2 (1998).

# Installation

The easiest way to install the library is using pip
```bash
pip install fim-python[gpu] #GPU version
```

If you don't have a compatible CUDA GPU, you can install the CPU only version to test the library, but the performance won't be comparable to the GPU version (see [Benchmark](#benchmark)).

```bash
pip install fim-python #CPU version
```

# Usage

The main interface to create a solver object to use is [`FIMPY.create_fim_solver`](https://fim-python.readthedocs.io/en/latest/interface.html#fimpy.solver.FIMPY.create_fim_solver)

```python
from fimpy.solver import FIMPY

#Create a FIM solver, by default the GPU solver will be called with the active list
#Set device='cpu' to run on cpu and use_active_list=False to use Jacobi method
fim = FIMPY.create_fim_solver(points, elems, D)
```

Example
-------

The following code reproduces the [above example](#details)

```python
import numpy as np
import cupy as cp
from fimpy.solver import FIMPY
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

#Create triangulated points in 2D
x = np.linspace(-1, 1, num=50)
X, Y = np.meshgrid(x, x)
points = np.stack([X, Y], axis=-1).reshape([-1, 2]).astype(np.float32)
elems = Delaunay(points).simplices
elem_centers = np.mean(points[elems], axis=1)

#The domain will have a small spot where movement will be slow
velocity_f = lambda x: (1 / (1 + np.exp(3.5 - 25*np.linalg.norm(x - np.array([[0.33, 0.33]]), axis=-1)**2)))
velocity_p = velocity_f(points) #For plotting
velocity_e = velocity_f(elem_centers) #For computing
D = np.eye(2, dtype=np.float32)[np.newaxis] * velocity_e[..., np.newaxis, np.newaxis] #Isotropic propagation

x0 = np.array([np.argmin(np.linalg.norm(points, axis=-1), axis=0)])
x0_vals = np.array([0.])

#Create a FIM solver, by default the GPU solver will be called with the active list
fim = FIMPY.create_fim_solver(points, elems, D)
phi = fim.comp_fim(x0, x0_vals)

#Plot the data of all points to the given x0 at the center of the domain
fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True)
cont_f1 = axes[0].contourf(X, Y, phi.get().reshape(X.shape))
axes[0].set_title("Distance from center")

cont_f2 = axes[1].contourf(X, Y, velocity_p.reshape(X.shape))
axes[1].set_title("Assumed isotropic velocity")
plt.show()
```

A general rule of thumb: If you only need to evaluate the eikonal equation once for a mesh, the Jacobi version (`use_active_list=False`) will probably be quicker since its initial overhead is low.
Repeated evaluations with different ![$\mathbf{x}_0$](https://latex.codecogs.com/svg.latex?\Large\mathbf{x}_0) or ![$D$](https://latex.codecogs.com/svg.latex?\Large%20D) favor the active list method for larger meshes.  
On the CPU, `use_active_list=True` outperforms the Jacobi approach for almost all cases.

# Documentation

[https://fim-python.readthedocs.io/en/latest](https://fim-python.readthedocs.io/en/latest)

# Citation

If you find this work useful in your research, please consider citing the paper in the [Journal of Open Source Software](https://joss.theoj.org/)
```bibtex
@article{grandits_fast_2021,
  doi = {10.21105/joss.03641},
  url = {https://doi.org/10.21105/joss.03641},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {66},
  pages = {3641},
  author = {Thomas Grandits},
  title = {A Fast Iterative Method Python package},
  journal = {Journal of Open Source Software}
}
```

# Benchmark

Below you can see a performance benchmark of the library for tetrahedral domains (cube in ND), triangular surfaces (plane in ND), and line networks (randomly sampled point cloud in the ND cube with successive minimum spanning tree) from left to right.
In all cases, ![$\mathbf{x}_0$](https://latex.codecogs.com/svg.latex?\Large\mathbf{x}_0) was placed in the middle of the domain.
The dashed lines show the performance of the implementation using active lists, the solid lines use the Jacobi method (computing all updates in each iteration).

![Preview](docs/figs/benchmark_gpu.jpg)

![Preview](docs/figs/benchmark_cpu.jpg)

The library works for an arbitrary number of dimensions (manifolds in N-D), but the versions for 2 and 3D received a few optimized kernels that speed up the computations.

The steps to reproduce the benchmarks can be found in the documentation at [https://fim-python.readthedocs.io/en/latest/benchmark.html](https://fim-python.readthedocs.io/en/latest/benchmark.html)

# Contributing

See [Contributing](CONTRIBUTING.md) for more information on how to contribute.

# License

This library is licensed under the [GNU Affero General Public License](LICENSE). 
If you need the library issued under another license for commercial use, you can contact me via e-mail [tomdev (at) gmx.net](mailto:tomdev@gmx.net).
# Contributing to FIM-Python

Thank you for your interest in FIM-Python. Any help and contributions are appreciated.


Reporting Bugs
---------------------

Please submit bug reports to the [issue page](https://github.com/thomgrand/fim-python/issues). Make sure that you include all of the following:
- Description of the bug
- Steps to reconstruct the error
- Operating system
- Version numbers of
  - Python
  - Numpy
  - Cupy

Fetching the version numbers and the operating system info can be automatically achieved by executing the following script in your python environment:

```python
import platform
import os
print("OS Info: %s, %s, v%s" % (os.name, platform.system(), platform.release()))

import numpy
print("Numpy version: %s" % (numpy.__version__))

try:
    import cupy
    print("GPU version, version of cupy: %s" % (cupy.__version__))
except ImportError:
    print("CPU version only")
```

Submitting Code
--------------------
FIM-Python uses the [pytest](https://docs.pytest.org) framework. Pip can take care of installing all necessary packages by listing the extra ``tests``:
```bash
pip install fim-python[gpu,tests]
```
The tests can be run by executing
```bash
python tests/generate_test_data.py #First time only to generate the test examples
python -m pytest tests
```

Before opening a pull request for newly written code, please make sure that all tests are passing.
In case you only have the CPU version, all tests for the GPU will be skipped. 
If you submit new features, please also write tests to ensure functionality of these features.  
The github-runner will also test pull-requests and committed versions of the library, but only on the CPU for the lack of a GPU on the runner.

> **_Note:_**  If you do **not** have a Cupy compatible GPU to test on, please clearly state this in your pull request, so somebody else from the community can test your code with all features enabled.
---
name: Bug report
about: Template for bug reports
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the unexpected behavior

**Please provide the following information:**
- Operating system
- Version numbers of
  - Python
  - Numpy
  - Cupy

See [Contributing.md](https://github.com/thomgrand/fim-python/blob/master/CONTRIBUTING.md) for a script to automatically output the variables

**Additional information [Optional]**
Add any other information about the problem here.
