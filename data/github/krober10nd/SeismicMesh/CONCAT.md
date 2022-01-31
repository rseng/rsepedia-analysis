<p align="center">
  <a href="https://github.com/krober10nd/SeismicMesh"><img alt="SeismicMesh" src="https://user-images.githubusercontent.com/18619644/92964244-28f31d00-f44a-11ea-9aa0-3d8ed0a1b60e.jpg" width="40%"></a>
  <p align="center">Create high-quality, simulation-ready 2D/3D meshes.</p>
</p>



[![status](https://joss.theoj.org/papers/ba94127ebbd0ca13c841f047fb5077bd/status.svg)](https://joss.theoj.org/papers/ba94127ebbd0ca13c841f047fb5077bd)
[![CircleCI](https://img.shields.io/circleci/project/github/krober10nd/SeismicMesh/master.svg?style=flat-square)](https://circleci.com/gh/krober10nd/SeismicMesh/tree/master)
[![ReadTheDocs](https://readthedocs.org/projects/seismicmesh/badge/?version=master)](https://seismicmesh.readthedocs.io/en/master/?badge=master)
[![CodeCov](https://codecov.io/gh/krober10nd/SeismicMesh/branch/master/graph/badge.svg)](https://codecov.io/gh/krober10nd/SeismicMesh)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/SeismicMesh.svg?style=flat-square)](https://pypi.org/pypi/SeismicMesh/)
[![PyPi]( https://img.shields.io/pypi/v/SeismicMesh.svg?style=flat-square)](https://pypi.org/project/SeismicMesh)
[![GPL](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


[SeismicMesh](https://github.com/krober10nd/SeismicMesh): Triangular Mesh generation in Python
==================================================================================================

SeismicMesh is a Python package for simplex mesh generation in two or three dimensions. As an implementation of [DistMesh](http://persson.berkeley.edu/distmesh/), it produces high-geometric quality meshes at the expense of speed. For increased efficiency, the core package is written in C++, works in parallel, and uses the [Computational Geometry Algorithms Library](https://doc.cgal.org/latest/Mesh_3/index.html). SeismicMesh can also produce mesh-density functions from seismological data to be used in the mesh generator.

*SeismicMesh* is distributed under the GPL3 license and more details can be found in our [short paper](https://github.com/krober10nd/SeismicMesh/blob/master/paper/paper.md).

Table of contents
=================

<!--ts-->
   * [Installation](#installation)
   * [Contributing](#contributing)
   * [Codebase](#codebase)
   * [Citing](#citing)
   * [Getting help](#problems)
   * [Examples](#examples)
     * [BP2004](#bp2004)
     * [EAGE Salt](#eage)
     * [Cylinder](#cylinder)
     * [Disk](#disk)
     * [Square](#square)
     * [Cube](#cube)
     * [Torus](#torus)
     * [Prism](#prism)
     * [Union](#union)
     * [Intersection](#intersection)
     * [Difference](#difference)
     * [Immersion](#immersion)
     * [Boundaries](#boundaries)
     * [Periodic](#periodic)
     * [Rotations](#rotations)
     * [Stretching](#stretching)
     * [Translation](#translation)
     * [Checking geometry](#checking)
   * [Parallelism](#parallelism)
   * [Performance comparison](#performance)
   * [Changelog](#changelog)
<!--te-->

Installation
============

For installation, SeismicMesh needs [CGAL](https://www.cgal.org/):

    sudo apt install libcgal-dev

After that, SeismicMesh can be installed from the Python Package Index
([pypi](https://pypi.org/project/SeismicMesh/)), so with:

    pip install -U SeismicMesh

If you'd like to read and write velocity models from segy/h5 format, you can install like:

    pip install -U SeismicMesh[io]

For more detailed information about installation and requirements see:

[Install](https://seismicmesh.readthedocs.io/en/master/install.html) -
How to install SeismicMesh.


Contributing
============

All contributions are welcome!

To contribute to the software:

1. [Fork](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) the repository.
2. Clone the forked repository, add your contributions and push the changes to your fork.
3. Create a [Pull request](https://github.com/krober10nd/SeismicMesh/pulls)

Before creating the pull request, make sure that the tests pass by running
```
tox
```
Some things that will increase the chance that your pull request is accepted:
-  Write tests.
- Add Python docstrings that follow the [Sphinx](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html).
- Write good commit and pull request messages.


[style]: https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html

Codebase
========

Here is a visual overview of the repository. An interactive version of this image can be found here: https://octo-repo-visualization.vercel.app/?repo=krober10nd%2FSeismicMesh

![Visualization of this repo](./diagram.svg)


Citing
=======

You may use the following BibTeX entry:
```
@article{Roberts2021,
  doi = {10.21105/joss.02687},
  url = {https://doi.org/10.21105/joss.02687},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {57},
  pages = {2687},
  author = {Keith J. Roberts and Rafael dos Santos Gioria and William J. Pringle},
  title = {SeismicMesh: Triangular meshing for seismology},
  journal = {Journal of Open Source Software}
}
```

Problems?
==========

If something isn't working as it should or you'd like to recommend a new addition/feature to the software, please let me know by starting an issue through the [issues](https://github.com/krober10nd/SeismicMesh/issues) tab. I'll try to get to it as soon as possible.


Examples
========

The user can quickly build quality 2D/3D meshes from seismic velocity
models in serial/parallel.

BP2004
-------
**WARNING: To run the code snippet below you must download the 2D BP2004
seismic velocity model and then you must uncompress it (e.g., gunzip).
This file can be downloaded from**
[here](http://s3.amazonaws.com/open.source.geoscience/open_data/bpvelanal2004/vel_z6.25m_x12.5m_exact.segy.gz)

![Above shows the mesh in ParaView that results from running the code below](https://user-images.githubusercontent.com/18619644/91606181-004a2e00-e948-11ea-83a4-e5ce05c7d82f.png)

```python
from mpi4py import MPI
import meshio

from SeismicMesh import get_sizing_function_from_segy, generate_mesh, Rectangle

comm = MPI.COMM_WORLD

"""
Build a mesh of the BP2004 benchmark velocity model in serial or parallel
Takes roughly 1 minute with 2 processors and less than 1 GB of RAM.
"""

# Name of SEG-Y file containg velocity model.
fname = "vel_z6.25m_x12.5m_exact.segy"

# Bounding box describing domain extents (corner coordinates)
bbox = (-12000.0, 0.0, 0.0, 67000.0)

# Desired minimum mesh size in domain
hmin = 75.0

rectangle = Rectangle(bbox)

# Construct mesh sizing object from velocity model
ef = get_sizing_function_from_segy(
    fname,
    bbox,
    hmin=hmin,
    wl=10,
    freq=2,
    dt=0.001,
    grade=0.15,
    domain_pad=1e3,
    pad_style="edge",
)

points, cells = generate_mesh(domain=rectangle, edge_length=ef)

if comm.rank == 0:
    # Write the mesh in a vtk format for visualization in ParaView
    # NOTE: SeismicMesh outputs assumes the domain is (z,x) so for visualization
    # in ParaView, we swap the axes so it appears as in the (x,z) plane.
    meshio.write_points_cells(
        "BP2004.vtk",
        points[:, [1, 0]] / 1000,
        [("triangle", cells)],
        file_format="vtk",
    )
```

Note SeismicMesh can also be used to write velocity models to disk in a hdf5 format using the function `write_velocity_model`. Following the previous example above with the BP2004 velocity model, we create an hdf5 file with a domain pad of 1000 m.

```python
from SeismicMesh import write_velocity_model

# Name of SEG-Y file containg velocity model.
fname = "vel_z6.25m_x12.5m_exact.segy"

# Bounding box describing domain extents (corner coordinates)
bbox = (-12000.0, 0.0, 0.0, 67000.0)

write_velocity_model(
     fname,
     ofname="bp2004_velocity_model",  # how the file will be called (with a .hdf5 extension)
     bbox=bbox,
     domain_pad=500,  # the width of the domain pad in meters
     pad_style="edge",  # how the velocity data will be extended into the layer
     units="m-s",  # the units that the velocity model is in.
 )
```


EAGE
----------

**WARNING: To run the code snippet below you must download (and uncompress) the 3D EAGE
seismic velocity model from (WARNING: File is \~500 MB)**
[here](https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/Salt_Model_3D.tar.gz)

**WARNING: Computationaly demanding! Running this example takes around 3 minutes in serial and requires
around 2 GB of RAM due to the 3D nature of the problem and the domain
size.**

![Above shows the mesh in ParaView that results from running the code below.](https://user-images.githubusercontent.com/18619644/103445790-52cd8b00-4c57-11eb-8bd4-4af8f24d4c88.jpg)

<!--pytest-codeblocks:skip-->
```python
from mpi4py import MPI
import zipfile
import meshio

from SeismicMesh import (
    get_sizing_function_from_segy,
    generate_mesh,
    sliver_removal,
    Cube,
)

comm = MPI.COMM_WORLD

# Bounding box describing domain extents (corner coordinates)
bbox = (-4200.0, 0.0, 0.0, 13520.0, 0.0, 13520.0)

# Desired minimum mesh size in domain.
hmin = 150.0

# This file is in a big Endian binary format, so we must tell the program the shape of the velocity model.
path = "Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/"
if comm.rank == 0:
    # Extract binary file Saltf@@ from SALTF.ZIP
    zipfile.ZipFile(path + "SALTF.ZIP", "r").extract("Saltf@@", path=path)

fname = path + "Saltf@@"

# Dimensions of model (number of grid points in z, x, and y)
nx, ny, nz = 676, 676, 210

cube = Cube(bbox)

# A graded sizing function is created from the velocity model along with a signed distance function by passing
# the velocity grid that we created above.
# More details can be found here: https://seismicmesh.readthedocs.io/en/master/api.html

ef = get_sizing_function_from_segy(
    fname,
    bbox,
    hmin=hmin,
    dt=0.001,
    freq=2,
    wl=5,
    grade=0.15,
    hmax=5e3,
    domain_pad=250,
    pad_style="linear_ramp",
    nz=nz,
    nx=nx,
    ny=ny,
    byte_order="big",
    axes_order=(2, 0, 1),  # order for EAGE (x, y, z) to default order (z,x,y)
    axes_order_sort="F",  # binary is packed in a FORTRAN-style
)

points, cells = generate_mesh(domain=cube, edge_length=ef, max_iter=75)

# For 3D mesh generation, we provide an implementation to bound the minimum dihedral angle::
# We use the preserve kwarg to ensure the level-set is very accurately preserved.
points, cells = sliver_removal(
    points=points, bbox=bbox, domain=cube, edge_length=ef, preserve=True
)

# Meshes can be written quickly to disk using meshio and visualized with ParaView::
if comm.rank == 0:

    # NOTE: SeismicMesh outputs assumes the domain is (z,x,y) so for visualization
    # in ParaView, we swap the axes so it appears as in the (x,y,z) plane.
    meshio.write_points_cells(
        "EAGE_Salt.vtk",
        points[:, [1, 2, 0]] / 1000.0,
        [("tetra", cells)],
    )
```


**The user can still specify their own signed distance functions and sizing functions to `generate_mesh` (in serial or parallel) just like the original DistMesh algorithm but now with quality bounds in 3D. Try the codes below!**



Cylinder
--------

<img alt="Cylinder" src="https://user-images.githubusercontent.com/18619644/97082301-0e7e9880-15df-11eb-9055-15394213d755.png" width="30%">

```python
# Mesh a cylinder
from mpi4py import MPI
import meshio

import SeismicMesh

comm = MPI.COMM_WORLD

hmin = 0.10

cylinder = SeismicMesh.Cylinder(h=1.0, r=0.5)

points, cells = SeismicMesh.generate_mesh(
    domain=cylinder,
    edge_length=hmin,
)

points, cells = SeismicMesh.sliver_removal(
    points=points,
    domain=cylinder,
    edge_length=hmin,
)

if comm.rank == 0:
    meshio.write_points_cells(
        "Cylinder.vtk",
        points,
        [("tetra", cells)],
        file_format="vtk",
    )
```

Disk
--------
<img alt="Disk" src="https://user-images.githubusercontent.com/18619644/97063883-b9a83700-1578-11eb-9cd7-3ff0cbac20d9.png" width="30%">


```python
# mesh a disk
import meshio
import SeismicMesh

disk = SeismicMesh.Disk([0.0, 0.0], 1.0)
points, cells = SeismicMesh.generate_mesh(domain=disk, edge_length=0.1)
meshio.write_points_cells(
    "disk.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```
Square
--------
<img alt="Square" src="https://user-images.githubusercontent.com/18619644/97063852-7b127c80-1578-11eb-97d5-cfe07cc969ec.png" width="30%">

```python
# mesh a square/rectangle
import meshio
import SeismicMesh

bbox = (0.0, 1.0, 0.0, 1.0)
square = SeismicMesh.Rectangle(bbox)
points, cells = SeismicMesh.generate_mesh(domain=square, edge_length=0.05)
meshio.write_points_cells(
    "square.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```
Cube
--------
<img alt="Cube" src="https://user-images.githubusercontent.com/18619644/97063751-e1e36600-1577-11eb-9387-613f3ae04bff.png" width="30%">

```python
# mesh a cuboid/cube
import meshio
import SeismicMesh

bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
cube = SeismicMesh.Cube(bbox)
points, cells = SeismicMesh.generate_mesh(domain=cube, edge_length=0.05)
points, cells = SeismicMesh.sliver_removal(points=points, domain=cube, edge_length=0.05)
meshio.write_points_cells(
    "cube.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```
Torus
--------
<img alt="Torus" src="https://user-images.githubusercontent.com/18619644/97063588-eeb38a00-1576-11eb-8cff-8e77ea4d2946.png" width="30%">


```python
# mesh a torus
import meshio
import SeismicMesh

hmin = 0.10

torus = SeismicMesh.Torus(r1=1.0, r2=0.5)
points, cells = SeismicMesh.generate_mesh(
    domain=torus,
    edge_length=hmin,
)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=torus, edge_length=hmin
)
meshio.write_points_cells(
    "torus.vtk",
    points,
    [("tetra", cells)],
)
```

<img alt="Torus" src="https://user-images.githubusercontent.com/18619644/97081705-8ac2ad00-15da-11eb-9466-a86216b8908c.png" width="30%">

Prism
--------
```python
# mesh a prism
import meshio

import SeismicMesh

hmin = 0.05
prism = SeismicMesh.Prism(b=0.5, h=0.5)

points, cells = SeismicMesh.generate_mesh(
    domain=prism,
    edge_length=hmin,
)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=prism, edge_length=hmin
)
meshio.write_points_cells(
    "prism.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```
Union
-----------------------------------
<img alt="Union" src="https://user-images.githubusercontent.com/18619644/97081772-045a9b00-15db-11eb-8356-7863cdf274a3.png" width="30%">

```python
# Compute the union of several SDFs to create more complex geometries
import meshio
import SeismicMesh

h = 0.10
rect0 = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 0.5))
rect1 = SeismicMesh.Rectangle((0.0, 0.5, 0.0, 1.0))
disk0 = SeismicMesh.Disk([0.5, 0.5], 0.5)
union = SeismicMesh.Union([rect0, rect1, disk0])
# Visualize the signed distance function
union.show()
points, cells = SeismicMesh.generate_mesh(domain=union, edge_length=h)
meshio.write_points_cells(
    "Lshape_wDisk.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```
Intersection
-------------------------------------------
<img alt="Leaf" src="https://user-images.githubusercontent.com/18619644/97081808-41bf2880-15db-11eb-9333-2d1230621c01.png" width="30%">

```python
# Compute the intersection of several SDFs to create more complex geometries
import meshio
import SeismicMesh

h = 0.05
rect0 = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))
disk0 = SeismicMesh.Disk([0.25, 0.25], 0.5)
disk1 = SeismicMesh.Disk([0.75, 0.75], 0.5)
intersection = SeismicMesh.Intersection([rect0, disk0, disk1])
points, cells = SeismicMesh.generate_mesh(domain=intersection, edge_length=h)
meshio.write_points_cells(
    "Leaf.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```
Difference
-------------------------------------------
<img alt="Hole" src="https://user-images.githubusercontent.com/18619644/97081829-69ae8c00-15db-11eb-815d-a8302f822337.png" width="30%">

```python
# Compute the difference of two SDFs to create more complex geometries.
import meshio
import SeismicMesh

h = 0.05
rect0 = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))
disk0 = SeismicMesh.Disk([0.5, 0.5], 0.1)
disk1 = SeismicMesh.Disk([0.75, 0.75], 0.20)
difference = SeismicMesh.Difference([rect0, disk0, disk1])
points, cells = SeismicMesh.generate_mesh(domain=difference, edge_length=h)
meshio.write_points_cells(
    "Hole.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```
Difference of Signed Distance Functions in 3-D
------------------------------------------------
<img alt="Cube wHoles" src="https://user-images.githubusercontent.com/18619644/97081862-ad08fa80-15db-11eb-94b2-801001137f1a.png" width="30%">

```python
# Compute the difference of several SDFs in 3D
import meshio
import SeismicMesh

h = 0.10
cube0 = SeismicMesh.Cube((0.0, 1.0, 0.0, 1.0, 0.0, 1.0))
ball1 = SeismicMesh.Ball([0.5, 0.0, 0.5], 0.30)
ball2 = SeismicMesh.Ball([0.5, 0.5, 0.0], 0.30)
ball3 = SeismicMesh.Ball([0.0, 0.5, 0.5], 0.30)
ball4 = SeismicMesh.Ball([0.5, 0.5, 0.5], 0.45)
difference = SeismicMesh.Difference([cube0, ball1, ball2, ball3, ball4])
points, cells = SeismicMesh.generate_mesh(domain=difference, edge_length=h, verbose=1)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=difference, edge_length=h, verbose=1
)
meshio.write_points_cells(
    "Cube_wHoles.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```
Immersion
--------------------
<img alt="Immersed disk" src="https://user-images.githubusercontent.com/18619644/99576017-37b0ff80-29b8-11eb-881d-a9b0dd0adc34.png" width="30%">

```python
# Immerse a subdomain so that it's boundary is conforming in the mesh.
import numpy as np

import meshio

import SeismicMesh

box0 = SeismicMesh.Rectangle((-1.25, 0.0, -0.250, 1.250))
disk0 = SeismicMesh.Disk([-0.5, 0.5], 0.25)

hmin = 0.10


fh = lambda p: 0.05 * np.abs(disk0.eval(p)) + hmin

points, cells = SeismicMesh.generate_mesh(
    domain=box0,
    edge_length=fh,
    h0=hmin,
    subdomains=[disk0],
    max_iter=100,
)
meshio.write_points_cells(
    "Square_wsubdomain.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```

Boundaries
-----------

Boundary conditions can also be prescribed and written to `gmsh` compatible files using `mehsio`. In the following example, we immerse a disk into the connectivity and then prescribe boundary conditions around the circle and each wall of the domain for later usage inside a finite element solver.

<img width="1221" alt="Screen Shot 2021-02-12 at 12 04 03 PM" src="https://user-images.githubusercontent.com/18619644/107784877-b1902500-6d2a-11eb-98f3-e01c1175f498.png">


```python
import numpy as np
import meshio
import SeismicMesh as sm

bbox = (0.0, 10.0, 0.0, 1.0)
channel = sm.Rectangle(bbox)
suspension = sm.Disk([0.5, 0.5], 0.25)

hmin = 0.10
fh = lambda p: 0.05 * np.abs(suspension.eval(p)) + hmin
points, cells = sm.generate_mesh(
    domain=channel,
    edge_length=fh,
    h0=hmin,
    subdomains=[suspension],
    max_iter=1000,
 )
# This gets the edges of the mesh in a winding order (clockwise or counterclockwise).
ordered_bnde = sm.geometry.get_winded_boundary_edges(cells)
# We use the midpoint of the edge to determine its boundary label
mdpt = points[ordered_bnde].sum(1) / 2
infl = ordered_bnde[mdpt[:, 0] < 1e-6, :]  # x=0.0
outfl = ordered_bnde[mdpt[:, 0] > 9.9 + 1e-6, :]  # x=10.0
walls = ordered_bnde[
    (mdpt[:, 1] < 1e-6) | (mdpt[:, 1] > 0.99 + 1e-6), :
]  # y=0.0 or y=1.0
cells_prune = cells[suspension.eval(sm.geometry.get_centroids(points, cells)) < 0]
circle = sm.geometry.get_winded_boundary_edges(cells_prune)

# Write to gmsh22 format with boundary conditions for the walls and disk/circle.
meshio.write_points_cells(
    "example.msh",
    points,
    cells=[
        ("triangle", cells),
        ("line", np.array(infl)),
        ("line", np.array(outfl)),
        ("line", np.array(walls)),
        ("line", np.array(circle)),
    ],
    field_data={
        "InFlow": np.array([11, 1]),
        "OutFlow": np.array([12, 1]),
        "Walls": np.array([13, 1]),
        "Circle": np.array([14, 1]),
    },
    cell_data={
        "gmsh:physical": [
            np.repeat(3, len(cells)),
            np.repeat(11, len(infl)),
            np.repeat(12, len(outfl)),
            np.repeat(13, len(walls)),
            np.repeat(14, len(circle)),
        ],
        "gmsh:geometrical": [
            np.repeat(1, len(cells)),
            np.repeat(1, len(infl)),
            np.repeat(1, len(outfl)),
            np.repeat(1, len(walls)),
            np.repeat(1, len(circle)),
        ],
    },
    file_format="gmsh22",
    binary=False,
)
```



Periodic
-------------
<img alt="Periodic torus" src="https://user-images.githubusercontent.com/18619644/101163708-bfcb1200-3612-11eb-9c6d-4f664a754d01.png" width="30%">

```python
# Repeat primitives to create more complex domains/shapes.
import SeismicMesh
import meshio

hmin = 0.30
bbox = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
torus = SeismicMesh.Torus(r1=1.0, r2=0.5)
# the Repeat function takes a list specifying the repetition period in each dim
periodic_torus = SeismicMesh.Repeat(bbox, torus, [2.0, 2.0, 2.0])
points, cells = SeismicMesh.generate_mesh(domain=periodic_torus, edge_length=hmin)
points, cells = SeismicMesh.sliver_removal(
    points=points, domain=periodic_torus, edge_length=hmin
)
meshio.write_points_cells(
    "periodic_torus.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```

Rotations
------------
<img alt="Rotated squares" src="https://user-images.githubusercontent.com/18619644/108713669-4e0ab200-74f7-11eb-925e-d92705327557.png" width="30%">

```python
# Rotate squares in 2D
import numpy as np

import meshio
import SeismicMesh

bbox = (0.0, 1.0, 0.0, 1.0)
rotations = np.linspace(-3.14, 3.14, 40)
squares = []
for _, rotate in enumerate(rotations):
    squares.append(SeismicMesh.Rectangle(bbox, rotate=rotate))

rotated_squares = SeismicMesh.Union(squares)

points, cells = SeismicMesh.generate_mesh(domain=rotated_squares, edge_length=0.05)
meshio.write_points_cells(
    "rotated_squares" + str(rotate) + ".vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)

```

<img alt="Rotated cubes" src="https://user-images.githubusercontent.com/18619644/108769631-03f5f080-7538-11eb-8db3-d215548496a8.png" width="30%">

```python
# Same as above but for cubes
import numpy as np

import meshio
import SeismicMesh

bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
rotations = np.linspace(-3.14, 3.14, 40)
cubes = []
for _, rotate in enumerate(rotations):
    cubes.append(SeismicMesh.Cube(bbox, rotate=rotate))

rotated_cubes = SeismicMesh.Union(cubes)

points, cells = SeismicMesh.generate_mesh(domain=rotated_cubes, edge_length=0.10)
meshio.write_points_cells(
    "rotated_cubes.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```

Stretching
------------

<img alt="Stretched squares" src="https://user-images.githubusercontent.com/18619644/109436519-ab729780-79fe-11eb-9656-1470f7c766b9.png" width="30%">

```python
# Geometric primitives can be stretched (while being rotated)
import meshio

from SeismicMesh import *

domain = Rectangle((0.0, 1.0, 0.0, 1.0), stretch=[0.5, 2.0], rotate=0.1*3.14)

points, cells = generate_mesh(domain=domain, edge_length=0.1, verbose=2)

meshio.write_points_cells(
    "stretched_square.vtk",
    points,
    [("triangle", cells)],
    file_format="vtk",
)
```

Translation
-----------
<img alt="A translated cuboid" src="https://user-images.githubusercontent.com/18619644/110262382-45ec5100-7f92-11eb-844e-fc0a963a1541.png" width="30%">

```python
# Geometric primitives can be translated (while being rotated and stretched)
import meshio

from SeismicMesh import *

cuboid = Cube(
    (0.0, 1.0, 0.0, 1.0, 0.1, 1.0),
    stretch=[1.5, 1.5, 1.5],
    translate=[0.5, 4.0, 1.0],
    rotate=4.5 * 3.14,
)
points, cells = generate_mesh(domain=cuboid, edge_length=0.10, max_iter=200)
points, cells = sliver_removal(points=points, domain=cuboid, edge_length=0.10, preserve=True)


meshio.write_points_cells(
    "stretched_square.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
```


Checking
--------

<img alt="Example of checking" src="https://user-images.githubusercontent.com/18619644/110243114-c336a800-7f37-11eb-813f-09c293bd721f.png" width="30%">

SeismicMesh's mesh generator is sensitive to poor geometry definitions and thus you should probably check it prior to complex expensive meshing. We enable all signed distance functions to be visualized via the ``domain.show()`` method where `domain` is an instance of a signed distance function primitive from `SeismicMesh.geometry`. Note: you can increase the number of samples to visualize the signed distance function by increasing the kwarg `samples` to the `show` method, which is by default set to 10000.

Parallelism
-----------

A simplified version of the parallel Delaunay algorithm proposed by [Peterka et. al 2014](https://dl.acm.org/doi/10.1109/SC.2014.86) is implemented inside the DistMesh algorithm, which does not consider sophisticated domain decomposition or load balancing yet. A peak speed-up of approximately 6 times using 11 cores when performing 50 meshing iterations is observed to generate the 33M cell mesh of the EAGE P-wave velocity model. Parallel performance in 2D is better with peak speedups around 8 times using 11 cores. While the parallel performance is not perfect at this stage of development, the capability reduces the generation time of this relatively large example (e.g., 33 M cells) from 91.0 minutes to approximately 15.6 minutes. Results indicate that the simple domain decomposition approach inhibit perfect scalability. The machine used for this experiment was an Intel Xeon Gold 6148 machine clocked at 2.4 GHz with 192 GB of RAM connected together with a 100 Gb/s InfiniBand network.

To use parallelism see the [docs](https://seismicmesh.readthedocs.io/en/par3d/tutorial.html#basics)

**See the paper/paper.md and associated figures for more details.**

Performance
------------

**How does performance and cell quality compare to Gmsh and CGAL mesh generators?

Here we use SeismicMesh 3.1.4, [pygalmesh](https://github.com/nschloe/pygalmesh) 0.8.2, and [pygmsh](https://github.com/nschloe/pygmsh) 7.0.0 (more details in the benchmarks folder).

Some key findings:

* Mesh generation in 2D and 3D using analytical sizing functions is quickest when using Gmsh but a closer competition for CGAL and SeismicMesh.
* However, using mesh sizing functions defined on gridded interpolants significantly slow down both Gmsh and CGAL. In these cases, SeismicMesh and Gmsh perform similarly both outperforming CGAL's 3D mesh generator in terms of mesh generation time.
* All methods produce 3D triangulations that have a minimum dihedral angle > 10 degrees enabling stable numerical simulation (not shown)
* Head over to the `benchmarks` folder for more detailed information on these experiments.

![Summary of the benchmarks](https://user-images.githubusercontent.com/18619644/99252088-38e20100-27ed-11eb-80b3-c10afac7efbf.png)

* **In the figure for the panels that show cell quality, solid lines indicate the mean and dashed lines indicate the minimum cell quality in the mesh.**

* Note: it's important to point out here that a significant speed-up can be achieved for moderate to large problems using the [parallel capabilities](https://seismicmesh.readthedocs.io/en/master/tutorial.html#basics) provided in SeismicMesh.


**For an additional comparison of *SeismicMesh* against several other popular mesh generators head over to [meshgen-comparison](https://github.com/nschloe/meshgen-comparison).


Changelog
=========

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## Unreleased
### Added
- Short blurb about using `write_velocity_model`
- Adding `read_velocity_model` to public API
### Fixed
- Bug fix to gradient limiting of mesh size functions

## [3.6.1]-2021-05-22
### Added
- Smoothed sets (e.g., intersections, differences, and unions)
- Conversion of velocity data from feet-second to meters-second
- Support for fixed points in iterative Laplacian mesh smoother.
### Improved
- Simplified pybind11 build system.
- Now using pytest-codeblocks instead of exdown

## [3.5.0]-2021-03-09
### Added
- Rotations for all geometric primitives
- Stretching for all geometric primitives
- Visuzlization of signed distance functions
### Fixed
- Support for Python 3.9
### Improved
- Fixed points in iterative Laplacian smooth

## [3.4.0]-2021-02-14
### Added
- Mesh improvement now solves Lapl. smoothing as a fixed-point problem using AMG solver.
- User can now mesh user-defined sizing functions in parallel (not from :class:SizeFunction)
- Ability to specify data type `dtype` of floating point number inside binary files.
- Example how to specify and write boundary conditions.
### Improved
- Faster unique edge calculation.

## [3.3.0]-2021-01-08
### Added
- Ability to improve accuracy of level-set when performing 3d sliver removal.
### Improved
- Marginally faster parallel speedup at scale in 2d/3d

## [3.2.0] -2020-12-14
### Added
- Adding basic periodic domains with the `Repeat` SDF.
- `sliver_removal` has optional variable step size when perturbing vertices. Helps to remove the "last sliver".
### Improved
- Faster rectangle and cube primitives.
- Reworking CPP code and bottlenecks...20-30% faster `generate_mesh` in parallel for 2D/3D from previous versions.

## [3.1.7] - 2020-11-27
### Improved
- Table of contents in README

### Added
- More testing of sliver removal and 2d mesh generation qualities.

### Fixed
- Disabled bug when doing Newton boundary projection at the end of 3d `sliver_removal`.

## [3.1.6] - 2020-11-26
### Bug present with sliver removal. Recommend to not use.
### Added
- Unit testing three versions of Python (3.6.1, 3.7.4, 3.8.1)


## [3.1.5] - 2020-11-24
- Support for constraining/immersing subdomains represented as signed distance functions.
- Faster cell manipulation operations for ~5-10% better speedups in parallel.
- Projection of points back onto level set.

## [3.1.4] - 2020-11-15
- Laplacian smoothing at termination for 2D meshing...significantly improves minimum cell quality.
- Made `hmin` a field of the SizeFunction class, which implies the user no longer needs to pass `h0` to
 `generate_mesh` or `sliver_removal`.

## [3.1.3] - 2020-11-06
### Fixed
- Cylinder radius and height are now correct.
- Torus, Prism, and Cylinder now have `dim` tag.

### Improved
- More control over the `grad` option in the mesh sizing function.

## [3.1.2] - 2020-11-04
### Improved
- Faster calculation of boundary vertices.
- More robust sliver removal in 3D.
### Fixed
- Corners are only constrained for constant resolution meshes

## [3.1.0] - 2020-10-28
### Added
- New geometric primitives--torus, wedge/prism, and cylinder.
- Updated images on README.
### Fixed
- Only constrain corners near 0-level set.
- Bug fix to 3D binary velocity reading.

## [3.0.6] - 2020-10-21
### Fixed
- Silence messages about pfix when verbose=0
### Added
- Added more examples on README
- New unions/intersections/differences with several SDF primivitives
- Automatic corner constraints in serial

## [3.0.5] - 2020-10-18
### Fixed
- Preserving fixed points in serial.
- Units in km-s detection warning bug.
- Docstring fixes to `generate_mesh`
- Improved mesh quality in 3D

### Added
- Automatic corner point extraction for cubes and rectangles.
- More support for reading binary files packed in a binary format.
- Check to make sure bbox is composed of all floats.

## [3.0.4] - 2020-10-12
### Added
- Improve conformity of level-set in final mesh through additional set of Newton boundary projection iterations.


More information
================

All other information is available at:
<https://seismicmesh.readthedocs.io>

[Getting
started](https://seismicmesh.readthedocs.io/en/master/overview.html) -
Learn the basics about the program and the application domain.

[Tutorials](https://seismicmesh.readthedocs.io/en/master/tutorial.html) -
Tutorials that will guide you through the main features.
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
reported by contacting the project team at krober@usp.br. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
---
title: 'SeismicMesh: Triangular meshing for seismology'
tags:
  - Python
  - Seismology
  - Mesh generation
  - DistMesh
  - Parallel mesh generation
authors:
  - name: Keith J. Roberts
    affiliation: 1
  - name: Rafael dos Santos Gioria
    affiliation: 1
  - name: William J. Pringle
    affiliation: 2
affiliations:
 - name: Research Center for Gas Innovation, Escola Politécnica da Universidade de São Paulo, São Paulo, Brazil.
   index: 1
 - name:  Dept. of Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, 156 Fitzpatrick Hall, Notre Dame, IN, U.S.A.
   index: 2
date: 20 August 2020
bibliography: paper.bib

---
# Summary

SeismicMesh is a Python package for simplex mesh generation in two or three dimensions. As an implementation of DistMesh [@doi:10.1137/S0036144503429121], it produces high-quality meshes at the expense of speed. For increased efficiency, the core package is written in C++, works in parallel, and uses the Computational Geometry Algorithms Library [@cgal:hs-chdt3-20b]. SeismicMesh can also produce mesh-density functions from seismological data to be used in the mesh generator.

# Background

Generating a high-quality graded mesh for a geophysical domain represents a challenge for seismological modeling using the finite element method (FEM). In these applications, a domain is discretized typically with triangular/tetrahedral elements that vary widely in size around features of interest. These meshes are commonly used with the FEM to solve partial differential equations that model acoustic or elastic waves, which are used in seismic velocity model building algorithms such as full waveform inversion (FWI) [@doi:10.1190/1.1441754; @virieux2009overview] and reverse time migration [@10.1093/gji/ggv380].

# Statement of Need

Despite the fact that many mesh generation programs exist such as Gmsh [@doi:10.1002/nme.2579] and CGAL [@cgal:rty-m3-20b; @cgal:r-ctm2-20b], it is uncommon to find capabilities that incorporate geophysical data into the mesh generation process to appropriately size elements. This in part contributes to the reality that automatic mesh generation for geophysical domains is not user-friendly.

Some packages have been created to script mesh generation from geophysical datasets such as in coastal ocean modeling [@roberts2019oceanmesh2d; @gorman2008systematic] and reservoir modeling [@cacace2015meshit]. In a similar manner, the aim of this package is to provide a straightforward Python package to script mesh generation directly from seismic velocity models. This is accomplished first by building a mesh density function using seismic velocity data and then supplying these inputs to a mesh generator that can use these inputs and operate at scale.

The mesh density function can be used as input other mesh generators. However, the usage of a sizing function can have significant impact on the mesh generation performance. For example, Gmsh’s advancing front and Delaunay refinement methods construct the mesh incrementally and do not permit vectorization, which leads to reduced performance at scale in 2D/3D. In contrast, the DistMesh algorithm takes advantage of vectorization when querying a complex mesh density function making it efficient and competitive to Gmsh for this kind of meshing problem.

# Core functionality

  1. The creation of 2D/3D graded mesh size functions defined on axis-aligned regular Cartesian grids. These mesh sizing functions encode mesh resolution distributions that conform to the variations from inputted seismic velocity model data and are distributed according to several heuristics [see @SeismicMeshDocs for further details]. Mesh size function grading is accomplished using @persson2006mesh.

  2. Distributed memory parallelism. The generation of potentially large (> 10 million cells) high-quality triangular or tetrahedral meshes using distributed memory parallelism with mesh resolution following sizing functions.

  3. An implementation of a 3D so-called sliver tetrahedral element removal technique [@tournois2009perturbing] to bound a mesh quality metric. Note that 2D mesh generation does not suffer from the formation of sliver elements.

Similar to other meshing programs such as Gmsh, SeismicMesh [@SeismicMeshDocs] enables generation of simplex meshes through a Python application programming interface.

The mesh's domain geometry is defined as the 0-level set of a signed distance function (SDF), which avoids the need to have explicit geometry information defining the boundary and can be particularly useful in geophysical domains.

# Performance Comparison

We compare the 2D/3D serial performance in terms of cell quality and mesh creation time between SeismicMesh, Gmsh [@doi:10.1002/nme.2579] and CGAL [@cgal:rty-m3-20b; @cgal:r-ctm2-20b]. The cell quality is defined as the product of the topological dimension of the mesh (2 or 3) and the incircle radius divided by the circumcircle radius and ranges between 0 and 1, where 1 is a perfectly symmetrical simplex. In mesh generation, there is always a trade-off between generation speed and mesh quality. We find that Gmsh produces high-quality meshes by far the fastest, SeismicMesh will produce meshes with the best quality, but much slower. Gmsh becomes comparatively slow when a user-defined mesh-density function is involved, which is SeismicMesh's primary use case.

For the two seismic domains (e.g., BP2004 and EAGE), SeismicMesh is faster than Gmsh for the 2D BP2004 benchmark but slightly slower for the 3D EAGE benchmark at scale. CGAL is not competitive for the 3D benchmark and is therefore not shown. Interpolant-based mesh sizing functions significantly slow the mesh generation time of Gmsh by a factor of $\sim 3$ as Gmsh calls the sizing function for each point individually (e.g., 95,756 times) whereas SeismicMesh does it for all points at once each meshing iteration (e.g., 26 times). 

![Using SeismicMesh V3.2.0, the mesh creation time (left columns) and resulting cell quality (right columns) for the four benchmarks studied over a range of problem sizes. For the panels that show cell quality, solid lines indicate the mean and dashed lines indicate the minimum cell quality in the mesh. \label{fig:benchmark}](Performance.pdf)

# Parallelism

A simplified version of the parallel Delaunay algorithm proposed by @peterka2014high is implemented inside the DistMesh algorithm, which does not consider sophisticated domain decomposition or load balancing yet. \autoref{fig:speedup} shows a peak speed-up of approximately 6 times using 11 cores when performing 50 meshing iterations to generate the 33M cell mesh of the EAGE P-wave velocity model. While the parallel performance is not perfect at this stage of development, the capability reduces the generation time of this relatively large example (e.g., 33 M cells) from 91.0 minutes to approximately 15.6 minutes. Results indicate that the simple domain decomposition approach inhibit perfect scalability. The machine used for this experiment was an Intel Xeon Gold 6148 machine clocked at 2.4 GHz  with 192 GB of RAM connected together with a 100 Gb/s InfiniBand network.

![The speedup (left-panel) as compared to the serial version of SeismicMesh V3.2.0 for a relatively light and heavy mesh each adapted to P-wave data from the EAGE Salt seismic velocity model. The total mesh generation wall-clock time is annotated in decimal minutes next to each point. The panel on the right hand side shows the mesh generation rate normalized by the number of total number of cells in the mesh. \label{fig:speedup}](Benchmarks.pdf)

# Ongoing and future applications

 Some future applications for this software:

 * SeismicMesh is being used by a group of researchers to build 2D/3D meshes for a seismological FEM model that has been developed in the Firedrake computing environment [@10.1145/2998441].

 * The usage of SDF to implicitly define the meshing domain presents potential use cases in a topology-optimization framework [@laurain2018level] for modeling the sharp interface of salt-bodies in seismological domains. In these applications, the 0-level set of a SDF is used to demarcate the boundary of the feature. Each inversion iteration, an optimization problem is solved to produce modifications to the location of the 0-level set. In this framework, SeismicMesh can be used within the inversion algorithm to generate and adapt meshes.

 * Much like how the original DistMesh program has been used, SeismicMesh can be adapted for other domain-specific applications besides seismology (e.g., fluid dynamics, astrophysics, and oceanography). An open source project project is already under way to use the same mesh generation technology for a Python version of OceanMesh2D to build industrial-grade meshes of coastal oceans [@roberts2019oceanmesh2d].

We expect future extensions of the program to introduce better domain decomposition algorithms to improve parallel performance.

# Acknowledgements

This research was carried out in association with the ongoing R&D project registered as ANP 20714-2, "Software technologies for modelling and inversion, with applications in seismic imaging" (University of São Paulo / Shell Brasil / ANP). We would like to thank the two reviewers who helped improve the presentation and quality of this manuscript and package greatly.

# References
Packages for benchmarks
------------------------
Some extra packages are required for benchmarking. These can be installed through the following command:
```
pip install SeismicMesh[benchmarking]
```

Performance
------------
Here we compare SeismicMesh against well-established existing mesh generation tools such as [cgal](https://doc.cgal.org/latest/Mesh_3/) and [gmsh](https://gmsh.info/doc/texinfo/gmsh.html). Specifically:

    * a comparison in mesh creation speed in terms of wall-clock time and throughput.
    * a comparison in cell quality.

where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 2 for triangles and 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.

Benchmarks
----------

The following set of benchmark problems are available:

    benchmark_BP2004: # meshing the 2D BP2004 seismic velocity model
    benchmark_EAGE: # meshing the 3D EAGE Salt velocity model
    benchmark_disk: # a uniform 2D mesh of a unit disk
    benchmark_ball: # a unit ball with a ring of higher resolution near the center.

Run `python benchmark_sphere.py` to run all benchmarks for a particular domain (e.g., sphere). Run `python benchmark_sphere.py --method METHODNAME` to select either CGAL using [pygalmesh](https://github.com/nschloe/pygalmesh), Gmsh using [pygmsh](https://github.com/nschloe/pygmsh) or `sm` to use SeismicMesh.

* Note: 3D CGAL results are not shown for the EAGE benchmark because they are several times slower than Gmsh and SeismicMesh. 2D CGAL results are not shown for BP2004 since they do not support user-defined variable mesh density functions.

Results
---------------

The computer used for benchmarking is a PC running MacOS with Dual-Core Intel Core i5 clocked at 2.00 GHz with 8GB of RAM. All mesh generation programs have been compiled similarly with gcc v8.3.0 with the -O3 option. These benchmarks have been done using CGAL v5.0, gmsh 4.7.0, and SeismicMesh v3.1.4. Each statistic is reported as the average of 5 executions.

Using [termplotlib](https://github.com/nschloe/termplotlib) and [meshplex](https://github.com/nschloe/meshplex) to calculate some mesh statistics, the benchmarks produce histograms of each cells' minimum [dihedral angles](https://en.wikipedia.org/wiki/Dihedral_angle) and histograms of cell quality (which was described above).

**NOTE: 2D mesh sizing functions are not supported by CGAL**

Average speed statistics can be computed via [pytest-benchmark](https://pypi.org/project/pytest-benchmark/) which is set up to run each domain 5 times. For example:

```python
py.test --benchmark-max-time=360 benchmarks/benchmark_ball.py
```

produces for example:

```
=========================================================================================== test session starts ============================================================================================
platform darwin -- Python 3.7.4, pytest-6.0.2, py-1.9.0, pluggy-0.13.1
benchmark: 3.2.3 (defaults: timer=time.perf_counter disable_gc=False min_rounds=5 min_time=0.000005 max_time=360 calibration_precision=10 warmup=False warmup_iterations=100000)
rootdir: /Users/Keith/junk/SeismicMesh, configfile: pytest.ini
plugins: xdist-2.1.0, cov-2.10.1, benchmark-3.2.3, forked-1.3.0
collected 3 items

benchmarks/benchmark_ball.py ...                                                                                                                                                                   [100%]

-------------------------------------------------------------------------------- benchmark: 3 tests --------------------------------------------------------------------------------
Name (time in s)          Min                Max               Mean            StdDev             Median               IQR            Outliers     OPS            Rounds  Iterations
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_gmsh              4.1087 (1.0)       4.6007 (1.0)       4.3494 (1.0)      0.1953 (1.0)       4.2846 (1.0)      0.2935 (1.0)           2;0  0.2299 (1.0)           5           1
test_seismic_mesh     14.3477 (3.49)     16.1393 (3.51)     15.2947 (3.52)     0.7355 (3.77)     15.3189 (3.58)     1.2345 (4.21)          2;0  0.0654 (0.28)          5           1
test_cgal             16.2922 (3.97)     19.2086 (4.18)     17.4405 (4.01)     1.2892 (6.60)     16.9541 (3.96)     2.2134 (7.54)          1;0  0.0573 (0.25)          5           1
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Legend:
  Outliers: 1 Standard Deviation from Mean; 1.5 IQR (InterQuartile Range) from 1st Quartile and 3rd Quartile.
  OPS: Operations Per Second, computed as 1 / Mean
================================================================================ 3 passed, 20 warnings in 186.87s (0:03:06) ===============================================================================
```


Details on experiments
-----------------------
* Mesh generation with CGAL is accomplished via [pygalmesh](https://github.com/nschloe/pygalmesh) version 0.8.2
* For CGAL's 3D mesh generator, all default quality options are assumed (e.g., facet angle bound of 30 degrees and the radius edge bound 2--to their theoretical limit). A `cell_size` function is passed to create variable mesh resolution in a way that is approximately equivalent to the mesh size function in SeismicMesh. In 2D we do not use Lloyd smoothing as it can significantly increase mesh generation time (but produce higher quality cells).
* Mesh generation with Gmsh is accomplished via [pygmsh](https://github.com/nschloe/pygmsh) 7.0.0 with all default options and, similar to CGAL, an approximately equivalent cell size function is passed.
* For SeismicMesh version 3.1.4, we perform all examples with 25 meshing iterations with a psuedo-timestep of 0.30 and then run the sliver removal implemention to bound the diheral angle to 10 degrees in 3D and, in 2D, we delete any lower quality elements on the boundary (with a cell quality less than 10 percent).
* All programs are executed in a seqeuntial mode. It's important to note however that a significant speed-up can be achieved for moderate to large problems using the [parallel capabilities](https://seismicmesh.readthedocs.io/en/master/tutorial.html#basics) provided in SeismicMesh. Threading based parallelism can be used with Gmsh and CGAL but these benchmarks have not been explored.
* The scripts with the prefix `run` iterate over a range of relevant problem sizes to produce the timining and quality scales at different scales.
==============
SeismicMesh
==============

:Author:  Keith Jared Roberts
:Contact: keithrbt0@gmail.com


Online Documentation
--------------------

Hosted at *Read the Docs* [http://SeismicMesh.readthedocs.org/]:

.. image::  https://readthedocs.org/projects/SeismicMesh/badge/?version=stable
   :target: https://SeismicMesh.readthedocs.org/en/stable/
.. image::  https://readthedocs.org/projects/SeismicMesh/badge/?version=latest
   :target: https://SeismicMesh.readthedocs.org/en/latest/

Hosted at *Github* [https://github/krober10nd/SeismicMesh]:

+ `User Manual (HTML)`_ (generated with Sphinx_).
+ `User Manual (PDF)`_  (generated with Sphinx_).
+ `API Reference`_      (generated with Epydoc_).

.. _User Manual (HTML): usrman/index.html
.. _User Manual (PDF):  SeismicMesh.pdf
.. _API Reference:      apiref/index.html

.. _Sphinx: http://sphinx-doc.org/
.. _Epydoc: http://epydoc.sourceforge.net/


Discussion and Support
----------------------

Hosted at Google Groups:

+ Group Page:   http://groups.google.com/group/SeismicMesh
+ Mailing List: SeismicMesh@googlegroups.com


Downloads and Development
-------------------------

Hosted at Github:

+ Project Site:    https://github.org/krober10nd/SeismicMesh
+ Source Releases: https://github.org/krober10nd/SeismicMesh/downloads
+ Issue Tracker:   https://github.org/krober10nd/SeismicMesh/issues
+ Git Repository:  https://github.org/krober10nd/SeismicMesh.git


Citations
---------

+ K.J. Roberts
  *SeismicMesh*
  to be submitted to Journal of Open Source Software (JOSS)


Acknowledgments
---------------

This project was supported by...
Overview
========


This software aims to create end-to-end workflows (e.g., from seismic velocity model to simulation ready mesh) to build quality two- and three-dimensional (2D and 3D) unstructured triangular and tetrahedral meshes for seismic domains. The main advantage of triangular meshes over cubes or rectangles is their ability to cost-effectively resolve variable material properties. These generated meshes are suitable for acoustic and elastic numerical wave propagators and a focus is placed on parallel unstructured mesh generation capabilities. A simple application program interfaces (API) written in Python enables the user to call parallel meshing algorithms that can be scaled up on distributed memory clusters.


*SeismicMesh* is currently being used to generate simplical meshes in 2D and 3D for acoustic and elastic wave propagators written using a Domain Specific Language called *Firedrake* [firedrake]_. These type of numerical simulations are used in Full Waveform Inversion, Reverse Time Migration, and Time Travel Tomography applications.


Mesh definition
-------------------------------------------
.. note ::
    Triangular and tetrahedral conformal meshes.

The domain :math:`\Omega` is partitioned into a finite set of cells :math:`\mathcal{T}_{h} = {T}` with disjoint interiors
such that

:math:`\cup_{T \in \mathcal{T}_{h}} T = \Omega`

Together, these cells form a mesh of the domain :math:`\Omega`. In our case, the cells are triangles and thus the mesh is a triangulation composed of :math:`nt` triangles and :math:`np` vertices in either two or three dimensional space. Note in 3D, the triangulation is composed of tetrahedral but we still refer to it as a triangulation. In 2D, a triangle has 3 vertices, 3 edges, and 1 face. In 3D, a tetrahedral has 4 vertices, 6 edges, and 4 facets. These cells :math:`t` are obtained by tessellating a set of vertices that lie in two or three dimensional space using the well-known and efficient Delaunay triangulation algorithms.


*High quality cells*
^^^^^^^^^^^^^^^^^^^^^^^
.. note ::
    A high quality cell has a cell quality of 1.

Any set of points can be triangulated, but the resulting triangulation will likely not be useable for numerical simulation. Thus, we strive to produce high-quality geometric meshes.

There are many definitions of mesh quality. Here I use the following formula to quantify how *well-shaped* the cells where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 2 for triangles and 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.


Software architecture
-------------------------------------------
.. note ::
    Python calls peformant C++ libraries like CGAL and Boost.

The software is implemented in a mixed language environment (Python and C++). The Python language is used for the API while computationally expensive operations are performed in C++. The two languages are linked together with *pybind11* and installation is carried out using *cmake*. The Computational Geometry Algorithms Library [cgal]_ is used to perform geometric operations that use floating point arithmetic to avoid numerical precision issues. Besides this, several common Python packages: *numpy*, *scipy*, *meshio*, *segyio*, and *mpi4py* are used.


Inputs
-------------------------------------------

Seismic velocity model
^^^^^^^^^^^^^^^^^^^^^^^^

A seismic velocity model is defined on an axis-aligned regular Cartesian grid in either 2D or 3D containing scalar values (typically the P-wave velocity speed at each grid point).

Currently, *SeismicMesh* can read seismic velocity models from SEG-y files or in binary format. The latter requires some more information; see the docstring.

Signed distance function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a point :math:`x`, the signed distance function returns the :math:`d` distance to the boundary of the domain :math:`∂ \Omega`.

Let :math:`\Omega ⊂ D ∈ R^N` be the domain in :math:`N` dimensions. The boundary of the domain to-be-meshed is represented as the 0-level set of a continuous function:

:math:`φ(·) : D → R.`

such that:

:math:`\Omega := {x ∈ D, φ(x) < 0}`

where :math:`φ : D × R+ → R` is Lipschitz continuous and called the level set function. If we assume :math:`|∇φ(·)| = 0` on the set :math:`{x ∈ D, φ(x) = 0}`, then we have :math:`∂ \Omega = {x ∈ D, φ(x) = 0}` i.e., the boundary :math:`∂ \Omega` is the 0-level set of :math:`φ(·)`. The property that :math:`|∇φ(·)| = 0` is satisfied if :math:`φ(·)` is a signed distance function.

.. note :: 
    We provide several simple signed distance functions: such as a Rectangle, Cube, and Disk. See the geometry module. 

Mesh sizing function
^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a point :math:`x`, the sizing function :math:`f(h)` returns the isotropic mesh size defined at :math:`x`. By mesh size, we specifically mean the triangular edge length nearby `x` assuming the triangles will be close to equilateral in the finished mesh.

The purpose of `get_sizing_function_from_segy` to build this function directly from the seismic velocity model provided.


*DistMesh* algorithm
-------------------------------------------

.. note ::
    This program uses a modified version of the *DistMesh* algorithm [distmesh]_ to generate simplical meshes.

For the generation of triangular meshes in 2D and 3D, we use a modified version of the *DistMesh* algorithm [distmesh]_. The algorithm is both simple and practically useful as it can produce high-geometric quality meshes in N-dimensional space. Further, by utilizing our approach to produce mesh size functions, the mesh generation algorithm is capable of generating high-quality meshes faithful to user-defined target sizing fields. A benefit of this is that mesh sizes can be built to respect numerically stability requirements a priori.

Briefly, the mesh generation algorithm is iterative and terminates after a pre-set number of iterations (e.g., 50-100). It commences with an initial distribution of vertices in the domain and iteratively relocates the vertices to create higher-geometric quality elements. The edges of the mesh act as *springs* that obey a constitutive law (e.g., Hooke's Law) otherwise referred to as a *force function*. During each meshing iteration, the discrepancy between the length of the edges in the mesh connectivity and their target length from the sizing function produce movement in the triangles' vertices.

The boundary of the domain is enforced by projecting any points that leave the domain back into it each meshing iteration. After a sufficient number of iterations, an equilibrium-like state is almost always approached and the movement of the vertices becomes relatively small. The equilibrium-like state of the mesh connectivity corresponds to a mesh that contains mostly isotropic equilateral triangles, which is critical for numerical simulation. However, as with nearly all mesh generators, a sequence of mesh improvement strategies are applied after mesh generation terminates to ensure the mesh will be robust for simulation.


Mesh adaptation
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning ::
    Functionality to adapt an existing mesh is a work in progress


3D *Sliver* removal
^^^^^^^^^^^^^^^^^^^^^^^^^^

3D Delaunay mesh generation algorithms form degenerate elements called *slivers*. If any *sliver* exists in a 3D mesh, the numerical solution can become unstable. Fortunately, this problem does not occur in 2D and, in 2D, a high quality mesh free of degenerate elements is easily achieved. To tackle this problem in 3D, a method similar to that of [slivers]_ was implemented. This algorithm aims at removing *slivers* while preserving the triangulation sizing distribution and domain boundary.


The *sliver* removal technique fits well within the *DistMesh* framework. For example, like the mesh generation approach, the algorithm operates iteratively. Each meshing iteration, it perturbs *only* vertices associated with *slivers* so that the circumspheres' radius of the *sliver* tetrahedral increases rapidly (i.e.., gradient ascent of the circumsphere radius) [slivers]_. The method operates on an existing mesh that ideally already has a high-mesh quality and is efficient since it uses CGAL's incremental Delaunay capabilities. The perturbation of a vertex of the *sliver* leads to a local combinational change in the nearby mesh connectivity to maintain Delaunayhood and almost always destroys the *sliver* in lieu of elements with larger dihedral angles.

.. note ::
    A *sliver* element is defined by their dihedral angle (i.e., angle between two surfaces) of which a tetrahedral has :math:`6`. Generally, if a 3D mesh has a minimum dihedral angle less than 1 degree, it will be numerically unstable. We've had success in simulating with meshes that have minimum dihedral angles of minimally around 5 degrees.


Parallelism and speed
-------------------------------------------

.. note ::
    This code uses distributed memory parallelism with the MPI4py package.

When constructing models at scale, the primary computational bottleneck in the *DistMesh* algorithm becomes the time spent in the Delauany triangulation algorithm, which occurs each iteration of the mesh generation step. The other steps involving the formation and calculation of the target sizing field and signed distance function are far less demanding. Using *mpi4py*, I implemented a simplified version of the [hpc_del]_ to parallelize the Delaunay triangulation algorithm. This approach scales well and reduces the time spent performing each meshing iteration thus making the approach feasible for large-scale 3D mesh generation applications. The domain is decomposed into axis-aligned *slices* than cut one axis of the domain. While this strategy doesn't fare well with load balancing, it simplifies the implementation and runtime communication cost associated with neighboring processor exchanges.

When possible, *SeismicMesh* uses low-level functionality from the CGAL package including the evaluation of geometric predicates, circumball calculations, polygonal intersection tests, and incremental triangulation capabilities.


References
-------------------------------------------

.. [hpc_del] Peterka, Tom, Dmitriy Morozov, and Carolyn Phillips. "High-performance computation of distributed-memory parallel 3D Voronoi and Delaunay tessellation." SC'14: Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis. IEEE, 2014.

.. [distmesh] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
              SIAM Review, Volume 46 (2), pp. 329-345, June 2004 (PDF)

.. [firedrake] Florian Rathgeber, David A. Ham, Lawrence Mitchell, Michael Lange, Fabio Luporini, Andrew T. T. Mcrae, Gheorghe-Teodor Bercea, Graham R. Markall, and Paul H. J. Kelly. Firedrake: automating the finite element method by composing abstractions. ACM Trans. Math. Softw., 43(3):24:1–24:27, 2016. URL: http://arxiv.org/abs/1501.01809, arXiv:1501.01809, doi:10.1145/2998441.

.. [cgal] The CGAL Project. CGAL User and Reference Manual. CGAL Editorial Board, 5.0.2 edition, 2020

.. [slivers] Tournois, Jane, Rahul Srinivasan, and Pierre Alliez. "Perturbing slivers in 3D Delaunay meshes." Proceedings of the 18th international meshing roundtable. Springer, Berlin, Heidelberg, 2009. 157-173.
.. _tutorial:

.. warning::

    Under construction. Contributions very welcome!

Basics
========

*SeismicMesh* supports the generation of both 2D and 3D meshes in
either serial or parallel from seismic velocity models.

Here I show how to build meshes from sizing functions created with the software and explain what the sizng options mean. The API for serial/parallel and 2D/3D is identical.

Assuming you've coded a short Python script to call *SeismicMesh* (similar to what is shown in the examples), you simply call the script with python for serial execution::

    python your_script.py

Distributed memory parallelism can be used by first writing an extra import statement for  ``mpi4py`` (``import mpi4py``) near your other imports. Following this write the following line near the top of your script before you call the *SeismicMesh* API)::

    comm = MPI.COMM_WORLD

.. note::
   This line has no effect on serial execution and its fine to leave it in if you intend to only use serial execution.

Parallel execution takes place by typing by::

    mpiexec -n N python your_script.py

where `N` is the number of cores (e.g., 2,3,4 etc.)

.. warning::
    Oversubscribing the mesh generation problem to too many cores will surely lead to problems and slow downs. In general, keeping the minimum number of vertices per rank to between 20-50k/rank results in optimal performance.


Example data
-------------------

.. note::
    Users should create a directory called `velocity_models` and place their seismic velocity models there.


A 2D model (BP2004)::

    wget http://s3.amazonaws.com/open.source.geoscience/open_data/bpvelanal2004/vel_z6.25m_x12.5m_exact.segy.gz

A 3D model (EAGE Salt)::

    https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/Salt_Model_3D.tar.gz


File I/O and visualization of meshes
--------------------------------------

Meshes can be written to disk in a variety of formats using the Python package `meshio` (https://pypi.org/project/meshio/).

.. warning::
    Note that *SeismicMesh* sizing function makes the assumption that the first dimension is `z` and the second is `x` while the third is `y`. This is done in this way since 2D seismological simulations take place in the z-x plane and 3D in the z-x-y plane. As a result, the meshes when loaded into visualization software will appear rotated 90 degrees. For visualization, we can output in the vtk format using MeshIO (as shown in the examples) and then load the vtk file into ParaView.

Some things to know
---------------------

This seismic velocity model is passed to the `get_sizing_function_from_segy` along with the domain extents ::

    from SeismicMesh import get_sizing_function_from_segy

    # Construct a mesh sizing function from a seismic velocity model
    ef = get_sizing_function_from_segy(fname, bbox, other-args-go-here,...)

* The user specifies the filename `fname` of the seismic velocity model (e.g., either SEG-y or binary)

* The user specifies the domain extents of the *velocity model* as a tuple of coordinates in meters representing the corners of the domain::

.. math::

    bbox = (z_{min}, z_{max}, x_{min}, x_{max})

* In 3D::

.. math::

    bbox = (z_{min}, z_{max}, x_{min}, x_{max}, y_{min}, y_{max})`

Geometry
---------

*SeismicMesh* can mesh any domain defined by a signed distance function. However, we provide some basic domain shapes: a rectangle, a cube, or a disk.

For example::

    from SeismicMesh import Rectangle, Cube, Disk

    rectangle = Rectangle(bbox)
    cube = Cube(bbox)
    disk = Disk(x0=[0,0],r=1) # center of (0,0) with a radius of 1.0

.. note::
    A good reference for various signed distance functions can be found at: https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm



Mesh size function
-------------------------------------------

.. note:
    Seismic velocity models often have different constant grid spacings in each dimension. The software considers this automatically based on the domain extents.

The notion of an adequate mesh size is determined by a combination of the physics of acoustic/elastic wave propagation, the desired numerical accuracy of the solution (e.g., spatial polynomial order, timestepping method, etc.), and allowable computational cost of the model amongst other things. In the following sub-sections, each available mesh sizing strategy is briefly described and pseudo code is provided.


.. note :: The final mesh size map is taken as the minimum of all supplied sizing functions.
.. note :: The mesh size map dictates the triangular edge lengths in the final mesh (assuming these triangles will be equilateral).

Wavelength-to-gridscale
^^^^^^^^^^^^^^^^^^^^^^^
The highest frequency of the source wavelet :math:`f_{max}` and the smallest value of the velocity model :math:`v_{min}` define the shortest scale length of the problem since the shortest spatial wavelength :math:`\lambda_{min}` is equal to the :math:`\frac{v_{min}}{f_{max}}`. For marine domains, :math:`v_{min}` is approximately 1,484 m/s, which is the speed of sound in seawater, thus the finest mesh resolution is near the water layer.

The user is able to specify the number of vertices per wavelength :math:`\alpha_{wl}` the peak source frequency :math:`f_{max}`. This sizing heuristic also  can be used to take into account varying polynomial orders for finite elements. For instance if using quadratic P=2 elements, :math:`\alpha_{wl}` can be safely be set to 5 to avoid excessive dispersion and dissipation otherwise that would occur with P=1 elements::

   # Construct mesh sizing object from velocity model
   ef = get_sizing_function_from_segy(fname, bbox,
       wl=3, # :math:`\alpha_{wl}` number of grid points per wavelength
       freq=2, # maximum source frequency for which the wavelength is calculated
   )



Resolving seismic velocity gradients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Seismic domains are known for sharp gradients in material properties, such as seismic velocity. These sharp gradients lead to reflections and refractions in propagated waves, which are critical for successful imaging. Thus, finer mesh resolution can be deployed inversely proportional to the local standard deviation of P-wave velocity. The local standard deviation of seismic P-wave velocity is calculated in a sliding window around each point on the velocity model. The user chooses the mapping relationship between the local standard deviation of the seismic velocity model and the values of the corresponding mesh size nearby it. This parameter is referred to as the :math:`grad` and is specified in meters.
For instance a :math:`grad` of 50 would imply that the largest gradient in seismic P-wave velocity is mapped to a minimum resolution of 50-m.::

    ef = get_sizing_function_from_segy(fname, bbox,
        grad=50, # the desired mesh size in meters near the sharpest gradient in the domain
    )

.. image:: SlopeStrat3D.jpg

.. note:

    The mapping of the local standard deviation of the gradient of seismic velocity is normalized to an interval of :math:`[0,1]` so that the largest gradient is assigned the mesh resolution indicated by :math`grad` and all other grad-to-mesh-sizes are associated using a linear relationship (with a slope of 1 and y-intercept of 0).


Courant-Friedrichs-Lewey (CFL) condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Almost all numerical wave propagators utilize explicit time-stepping methods in the seismic domain. The major advantage for these explicit methods is computational speed. However, it is well-known that all explicit or semi-explicit methods require that the Courant number be bounded above by the Courant-Friedrichs-Lewey (CFL) condition. Ignoring this condition will lead to a numerically unstable simulation. Thus, we must ensure that the Courant number is indeed bounded for the overall mesh size function.

After sizing functions have been activated, a conservative maximum Courant number is enforced.

For the linear acoustic wave equation assuming isotropic mesh resolution, the CFL condition is commonly described by

.. math::

    C_{r}(x) = \frac{(\Delta t*v_p(x))}{dim*h(x)}

where :math:`h` is the diameter of the circumball that inscribes the element either calculated from :math:`f(h)` or from the actual mesh cells, :math:`dim` is the spatial dimension of the problem (2 or 3), :math:`\Delta t` is the intended simulation time step in seconds and :math:`v_p` is the local seismic P-wave velocity. The above equation can be rearranged to find the minimum mesh size possible for a given :math:`v_p` and :math:`\Delta t`, based on some user-defined value of :math:`Cr \leq 1`. If there are any violations of the CFL, they can bed edited before building the mesh so to satisfy that the maximum :math:`Cr` is less than some conservative threshold. We normally apply :math:`Cr = 0.5`, which provides a solid buffer but this can but this can be controlled by the user like the following::


    ef = get_sizing_function_from_segy(fname, bbox,
        cr=0.5, # maximum bounded Courant number to be bounded in the mesh sizing function
        dt=0.001, # for the given :math:`\Delta t` of 0.001 seconds
        ...
    )

Further, the space order of the method (:math:`p`) can also be incorporated into the above formula to consider the higher spatial order that the simulation will use::

    ef = get_sizing_function_from_segy(fname, bbox,
        cr=0.5, # maximum bounded Courant number :math:`Cr_{max}` in the mesh
        dt=0.001, # for the given :math:`\Delta t` of 0.001 seconds
        space_order = 2, # assume quadratic elements :math:`P=2`
        ...
    )

The above code implies that the mesh will be used in a simulation with :math:`P=2` quadratic elements, and thus will ensure the :math:`Cr_{max}` is divided by :math:`\frac{1}{space\_order}`


Mesh size gradation
^^^^^^^^^^^^^^^^^^^^^^^

In regions where there are sharp material contrasts, the variation in element size can become substantially large, especially using the aforementioned sizing strategies such as the wavelength-to-gridscale. Attempting to construct a mesh with such large spatial variations in mesh sizes would result in low-geometric quality elements that compromise the numerical stability of a model.

Thus, the final stage of the development of a mesh size function :math:`h(x)` involves ensuring a size smoothness limit, :math:`g` such that for any two points :math:`x_i`, :math:`x_j`, the local increase in size is bounded such as:

 :math:`h(\boldsymbol{x_j}) \leq h(\boldsymbol{x_i}) + \alpha_g||\boldsymbol{x_i}-\boldsymbol{x_j}||`

A smoothness criteria is necessary to produce a mesh that can simulate physical processes with a practical time step as sharp gradients in mesh resolution typically lead to highly skewed angles that result in poor numerical performance.

We adopt the method to smooth the mesh size function originally proposed by [grading]_. A smoother sizing function is congruent with a higher overall element quality but with more triangles in the mesh. Generally, setting :math:`0.2 \leq \alpha_g \leq 0.3` produces good results::

   ef = get_sizing_function_from_segy(fname, bbox,
       ...
       grade=0.15, # :math:`g` cell-to-cell size rate growth bound
       ...
   )

.. image:: ExGrade3D.jpg

Domain padding
^^^^^^^^^^^^^^^^^^^

.. note::

    It is assumed that the top side of the domain represents the free-surface thus no domain padding applied there.

In seismology applications, the goal is often to model the propagation of an elastic or acoustic wave through an infinite domain. However, this is obviously not possible so the domain is approximated by a finite region of space. This can lead to undesirable artificial reflections off the sides of the domain however. A common approach to avoid these artificial reflections is to pad the domain and enforce absorbing boundary conditions in this extension. In terms of meshing to take this under consideration, the user has the option to specify a domain extension of variable width on all three sides of the domain like so::

   ef = get_sizing_function_from_segy(fname, bbox,
       domain_pad=250, # domain will be pad by 250-m on all three sides of the domain
       ...
   )

In this domain pad, mesh resolution can be adapted according to following three different styles.

 * ``Linear`` - pads the seismic velocities on the edges of the domain linearly increasing into the domain pad.

 * ``Constant`` - places a constant velocity of the users selection in the domain pad.

 * ``Edge`` - pads the seismic velocity about the domain boundary so that velocity profile is identical to its edge values.

An example of the ``edge`` style is below::

   ef = get_sizing_function_from_segy(fname, bbox,
       domain_pad=250, # domain will be extended by 250-m on all three sides
       padstyle="edge", # velocity will be extend from values at the edges of the domain
       ...
   )

.. note::

    In our experience, the ``edge`` option works the best at reducing reflections with absorbing boundary conditions.

.. image:: domainext.png



Mesh generation
-------------------------------------------

.. warning:
    Connectivity is made approximately deterministic as each instance of mesh generation uses
    the same ``seed=0``. The user can specify the seed if they like.

After building your signed distance function and sizing function, call the ``generate_mesh`` function to generate the mesh::

    points, cells = generate_mesh(domain=rectangle, edge_length=ef, h0=hmin)

.. note::
    `ef` is a sizing function created using get_sizing_function_from_segy

You can change how many iterations are performed by altering the kwarg `max_iter`::

    points, cells = generate_mesh(domain=rectangle, edge_length=ef, h0=hmin, max_iter=100)

.. note :: Generally setting max_iter to between 50 to 100 iterations produces a high quality triangulation. By default it runs 50 iterations.

When executing in parallel, the user can optionally choose which axis (0, 1, or 2 [if 3D]) to decompose the domain::

    points, cells = generate_mesh(domain=cube, edge_length=ef, h0=hmin, max_iter=100, axis=2)

.. note :: Generally axis=1 works the best in 2D or 3D since typically mesh sizes increase in size from the free surface to the depths of the model. In this situation, the computational load tends to be better balanced.



Mesh improvement (*sliver* removal)
-------------------------------------------


3D *Sliver* removal
^^^^^^^^^^^^^^^^^^^^^^^

It is strongly encouraged to run the sliver removal method by passing the point of set of a previously generated mesh::

    points, cells = sliver_removal(
        points=points, domain=cube, h0=minimum_mesh_size, edge_length=ef
    )

.. note:: Please remember to import this method at the top of your script (e.g., `from SeismicMesh import sliver_removal`)

By default, ``min_dh_angle_bound`` is set to :math:`10`. The sliver removal algorithm will attempt 50 iterations but will terminate earlier if no slivers are detected. Generally, if more than 50 meshing iterations were used to build the mesh, this algorithm will converge in 10-20 iterations.

.. warning:: Do not set the minimum dihedral angle bound greater than 15 unless you've already successfully ran the mesh with a lower threshold. Otherwise, the method will likely not converge.


References
______________

.. [grading] Persson, Per-Olof. "Mesh size functions for implicit geometries and PDE-based gradient limiting."
                Engineering with Computers 22.2 (2006): 95-109.
SeismicMesh
==============

:Author:       Keith Roberts
:Contact:      krober@usp.br
:Web Site:     https://github.com/krober10nd/SeismicMesh
:Date:         |today|

.. include:: abstract.txt

Contents
========

.. include:: toctree.txt


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


Installation
============

Requirements
-------------

You need to have the following software properly installed in order to
build *SeismicMesh*:

* Python >= 3.0
.. note ::
    The file ``setup.cfg`` in the main directory indicates all the Python dependencies and their respective version numbers (if necessary). These packages should be installed at compile time by setuptools

.. note ::
    On some Linux systems, users may have to resort to `apt install python3-segyio` to installing segyio on their systems.

* pybind11 >= 2.6

* C++ compiler (GNU or Intel) with support for std++14 or newer.

* cmake >= 3.0

* CGAL >= 5.0.0 which requires:

    * MPFR

    * GMP

    * Boost > 1.4.8

.. note ::
    CGAL requires Boost, MPFR and GMP. These packages may already be installed on your standard Linux box.

.. warning ::
    Make sure your package manager is downloading CGAL >= 5.0 otherwise you will not be able to install *SeismicMesh*!




Compilation by source
----------------------

After installing all dependencies, perform ::

$  pip install -e .

.. note ::
    If you do not have administrative rights on your system, add the flag ``--user`` to the command above.

.. warning ::
    With this said, the preferred method of installation using pypi: pip install SeismicMesh

Testing
-------

Testing is accomplished with `tox`. The `tox` package can be installed like so::

    pip install tox

To test the installation, serial and parallel capabilities, you can use `tox` from the top directory of the package::

$ tox

Installation on Clusters
--------------------------

.. note::
    Make sure the CXX environment variable points to your intended compiler!

If installing on a cluster by source with a local installation of ``CGAL`` and ``Boost``, you'll need to edit ``setup.cfg`` with the CMake arguments so to point the installation to the correct directories. Namely, in ``setup.py`` you'll have to edit the list called ``cmake_args`` to include ::

  -DCMAKE_CXX_COMPILER=+/PATH/TO/CPPCOMPILER

  -DBoost_INCLUDE_DIR=+/PATH/TO/BOOST/

  -DMPFR_LIBRARIES=+/PATH/TO/libmpfr.a

  -DMPFR_INCLUDE_DIR=+/PATH/TO/MPFR/include
Introduction
============

Generating a high-quality graded mesh for a geophysical domain represents a modern challenge for sophisticated geophysical simulation workflows.
In these applications, a domain is discretized typically with simplicial elements (e.g., triangles/tetrahedrals)
that adapt in size to features of interest. These meshes are commonly used with Finite Element and Finite Volume methods to solve
Partial Differential Equations (PDEs) that model physical processes such as the acoustic or elastic wave equation. These kind of simulations are
often used in geophysical exploration studies to solve inverse problems such as Full Waveform Inversion.

There are many aspects to consider when building a mesh for a geophysical inverse problem. Mesh elements must be sufficiently well-shaped,
sized to maintain numerical stability, minimize numerical dissipation, and maximize physical and numerical accuracy. For applications such as Travel Time Tomography and Reverse Time Migration, material discontinuities need to be well represented to ensure reflections are as accurate as possible. In cases with complex and irregular rock structures, explicit geometries may not exist complicating mesh generation workflows and requiring the use of graphical user interfaces to create these geometries.

The purpose of this software is quickly assemble meshes for geophysical simulation using the Finite Element/Finite Volume methods. This is accomplished by giving users simple controls to design their own meshes without having to learn or write much code or more complex software. This program strives for automation and reproduction of the mesh.

This package contains the technology to build a 2D/3D mesh from a seismic velocity model in a scriptable manner. An approach to build mesh sizing functions that can lead to high-quality, graded meshes with the generator is detailed. For geometry creation, we rely on signed distance functions to define domain features implicitly, which avoids the need to supply explicit locations of segments/surfaces.
Modules
==============

Here we document the public API

*SeismsicMesh.geometry*
-------------------------------

Routines to perform geometrical/topological operations and calculate things on meshes.

.. automodule:: SeismicMesh.geometry
    :members:

*SeismsicMesh.sizing*
---------------------------------------------

Function to build a :math:`f(h)` mesh sizing function from a seismic velocity model.
Assumes the domain can be represented by a rectangle (2D) or cube (3D) and thus builds a :math:`f(d)` accordingly.

.. automodule:: SeismicMesh.sizing
    :members:

*SeismsicMesh.generation*
-------------------------------

Functions to build and improve a simplical mesh that conforms to the signed distance function :math:`f(d)`
and :math:`f(h)`.

.. automodule:: SeismicMesh.generation
    :members:
