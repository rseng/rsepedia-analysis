# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our
project and our community a harassment-free experience for everyone,
regardless of age, body size, disability, ethnicity, gender identity and
expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

-   Using welcoming and inclusive language
-   Being respectful of differing viewpoints and experiences
-   Gracefully accepting constructive criticism
-   Focusing on what is best for the community
-   Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

-   The use of sexualized language or imagery and unwelcome sexual
    attention or advances
-   Trolling, insulting/derogatory comments, and personal or political
    attacks
-   Public or private harassment
-   Publishing others\' private information, such as a physical or
    electronic address, without explicit permission
-   Other conduct which could reasonably be considered inappropriate in
    a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of
acceptable behavior and are expected to take appropriate and fair
corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that they
deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public
spaces when an individual is representing the project or its community.
Examples of representing a project or community include using an
official project e-mail address, posting via an official social media
account, or acting as an appointed representative at an online or
offline event. Representation of a project may be further defined and
clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may
be reported by contacting the project team at
<n.renaud@esciencecenter.nl>. All complaints will be reviewed and
investigated and will result in a response that is deemed necessary and
appropriate to the circumstances. The project team is obligated to
maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted
separately.

Project maintainers who do not follow or enforce the Code of Conduct in
good faith may face temporary or permanent repercussions as determined
by other members of the project\'s leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor
Covenant](https://www.contributor-covenant.org), version 1.4, available
at
<https://www.contributor-covenant.org/version/1/4/code-of-conduct.html>
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment
or question to a full fledged [pull
request](https://help.github.com/articles/about-pull-requests/). Please
read and follow our [Code of Conduct](CODE_OF_CONDUCT.rst).

A contribution can be one of the following cases:

1.  you have a question;
2.  you think you may have found a bug (including unexpected behavior);
3.  you want to make some kind of change to the code base (e.g. to fix a
    bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1.  use the search functionality
    [here](https://github.com/hpgem/nanomesh/issues) to see if someone
    already filed the same issue;
2.  if your issue search did not yield any relevant results, make a new
    issue;
3.  apply the \"Question\" label; apply other labels when relevant.

## You think you may have found a bug

1.  use the search functionality
    [here](https://github.com/hpgem/nanomesh/issues) to see if someone
    already filed the same issue;

2.  if your issue search did not yield any relevant results, make a new issue,
    making sure to provide enough information to the rest of the community to
    understand the cause and context of the problem. Depending on the issue,
    you may want to include:

    -  the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas)
       of the commit that is causing your problem;
    -  some identifying information (name and version number) for
       dependencies you\'re using;
    -  information about the operating system;

3.  apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1.  (**important**) announce your plan to the rest of the community
    *before you start working*. This announcement should be in the form
    of a (new) issue;
2.  (**important**) wait until some kind of consensus is reached about
    your idea being a good idea;
3.  if needed, fork the repository to your own Github profile and create
    your own feature branch off of the latest master commit. While
    working on your feature branch, make sure to stay up to date with
    the master branch by pulling in changes, possibly from the
    \'upstream\' repository (follow the instructions
    [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/)
    and [here](https://help.github.com/articles/syncing-a-fork/));
4.  make sure the existing tests still work by running
    `pytest`;
5.  add your own tests (if necessary);
6.  update or expand the documentation;
7.  [push](http://rogerdudler.github.io/git-guide/) your feature branch
    to (your fork of) the nanomesh repository on GitHub;
8.  create the pull request, e.g. following the instructions
    [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you\'ve made a valuable contribution, but you
don\'t know how to write or run tests for it, or how to generate the
documentation: don\'t let this discourage you from making the pull
request; we can help you! Just go ahead and submit the pull request, but
keep in mind that you might be asked to append additional commits to
your pull request.
[![Documentation Status](https://readthedocs.org/projects/nanomesh/badge/?version=latest)](https://nanomesh.readthedocs.io/en/latest/?badge=latest)
[![tests](https://github.com/hpgem/nanomesh/actions/workflows/test.yaml/badge.svg)](https://github.com/hpgem/nanomesh/actions/workflows/test.yaml)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/nanomesh)](https://pypi.org/project/nanomesh/)
[![PyPI](https://img.shields.io/pypi/v/nanomesh.svg?style=flat)](https://pypi.org/project/nanomesh/)
[![DOI](https://zenodo.org/badge/311460276.svg)](https://zenodo.org/badge/latestdoi/311460276)

![Nanomesh banner](./notebooks/other/banner.png)

# Nanomesh

Nanomesh is a Python workflow tool for generating meshes from 2D and 3D image data. It has an easy-to-use API that can help process and segment image data, generate quality meshes (triangle / tetrahedra), and write the data to many mesh formats. Nanomesh also contains tools to inspect the meshes, visualize them, and generate cell quality metrics.

- Easy-to-use Python API
- Segment and mesh 2D or 3D image data
- Mesh visualization
- Calculate and plot cell metrics
- Export to many mesh formats

Documentation: https://nanomesh.readthedocs.io/en/latest/

## Show me an example!

This example shows the workflow for generating a mesh from segmented data, and demonstrates a few of the features of Nanomesh. It uses a synthetic binary image with several rounded blob-like objects generated by [skimage](https://scikit-image.org/).

```pycon
>>> from skimage.data import binary_blobs
>>> from nanomesh import Image
>>>
>>> blobs = binary_blobs(length=100, volume_fraction=0.25, seed=2102)
>>> plane = Image(blobs)
>>>
>>> print(plane)
Plane(shape=(100, 100), range=(False,True), dtype=bool)
```

[`Image`](https://nanomesh.readthedocs.io/en/latest/api.image_data.html#nanomesh.Image) is essentially a container for a [`numpy`](https://numpy.org/) array with some methods for image segmentation and visualization.

```pycon
>>> plane.show()
<AxesSubplot:xlabel='x', ylabel='y'>
```

<img src="notebooks/other/hello_world_files/hello_world_5_1.png" alt="drawing" width="50%"/>

Generating a mesh from image data is simple in Nanomesh using [`Plane.generate_mesh()`](https://nanomesh.readthedocs.io/en/latest/api.meshing.html#nanomesh.plane2mesh). The options `opts` are passed to the triangulation function ([`nanomesh.triangulate`](https://nanomesh.readthedocs.io/en/latest/api.helpers.html#nanomesh.triangulate)). In this example, we use `q30` to generate a quality mesh with minimum angles of 30°, and `a50` to limit the triangle size to 50 pixels.

The returned `mesh` is a [`MeshContainer`](https://nanomesh.readthedocs.io/en/latest/api.mesh_data.html#nanomesh.MeshContainer) that contains the generated triangles and line segments.

```pycon
>>> mesh = plane.generate_mesh(opts='q30a10')
>>> mesh
<MeshContainer>
  Number of points: 932
  Number of cells:
    triangle: 1754
    line: 2685
  Point data: physical
  Cell data: physical
  Field data: feature, background
```

In the next cell, we plot the triangles.

```pycon
>>> mesh.plot('triangle')
<AxesSubplot:title={'center':'triangle mesh'}>
```

<img src="notebooks/other/hello_world_files/hello_world_9_1.png" alt="drawing" width="50%"/>

With the [metrics submodule](https://nanomesh.readthedocs.io/en/latest/api.metrics.html), Nanomesh can also calculate cell quality metrics and show them as a [colored triangle](https://nanomesh.readthedocs.io/en/latest/api.metrics.html#nanomesh.metrics.plot2d) or [histogram plot](https://nanomesh.readthedocs.io/en/latest/api.metrics.html#nanomesh.metrics.histogram).

```pycon
>>> from nanomesh import metrics
>>> triangle_mesh = mesh.get('triangle')
>>> metrics.histogram(triangle_mesh, metric='radius_ratio')
<AxesSubplot:title={'center':'Histogram of radius ratio'}, xlabel='Radius ratio', ylabel='frequency'>
```

<img src="notebooks/other/hello_world_files/hello_world_11_1.png" alt="drawing" width="50%"/>

Nanomesh uses [meshio](https://github.com/nschloe/meshio) to write data to most meshing formats.

```pycon
>>> mesh.write('mesh.vtk')
Warning: VTK requires 3D points, but 2D points given. Appending 0 third component.
```

That's it! There is a lot more that Nanomesh can do, check out [the examples](https://nanomesh.readthedocs.io/en/latest/examples/index.html) for an overview.

## Installation

One of the goals for Nanomesh is that it is easy to install.
This means that all dependencies are available from [PyPi](https://pypi.org).

If you use conda, it is advised to create a new environment:

```
conda create -n nanomesh python=3.9
conda activate nanomesh
```

Install nanomesh:

```
pip install nanomesh
```

For the full installation instructions, see the [installation guidelines](https://nanomesh.readthedocs.io/en/latest/install.html).

### Development

Nanomesh does not have any hard version constraints. For development, it is
still useful to have a consistent environment.
Therefore, Nanomesh uses
a constraints file (`constraints.txt`) which pins the version requirements.

The constraints are automatically [updated and tested every month](https://github.com/hpgem/nanomesh/actions/workflows/update_dependencies.yaml).
Note that in case you run into issues, you may also try to install
Nanomesh with constraints file.

Install `nanomesh` using the development dependencies:

`pip install -e .[develop] -c constraints.txt`

Running the tests using [pytest](https://docs.pytest.org/):

`pytest`

Linting and checks are done using [pre-commit](https://pre-commit.com):

`pre-commit`

Building the docs:

`make html --directory docs`
Contains apidoc templates from https://github.com/sphinx-doc/sphinx/tree/master/sphinx/templates/apidoc
This directory stores the jupyter notebooks.

The source for the notebooks are stored in the *markdown* files.
Any code changes should be visible in the source files.
They can be converted both ways using [jupytext](https://github.com/mwouts/jupytext).
The notebooks are used for documentation via sphinx.

Each notebook is paired to a markdown file, so any changes to the notebook will be automatically saved in both files.

### Convert from ipynb to markdown

    jupytext .\hello_world.ipynb --from ipynb --to md

    Get-ChildItem "**/*.ipynb" | Foreach-Object { jupytext $_ --from ipynb --to md }


### Formatting notebooks with yapf

    jupytext .\hello_world.md --sync --pipe yapf

    Get-ChildItem "**/*.ipynb" | Foreach-Object { jupytext $_ --sync --pipe yapf }


### Execute notebooks commands

Convert notebooks using `jupytext` using the `--execute` flag.

This will use the md file as source. The output is stored in the notebook with the same name.
Any existing ipynb will be overwritten.

    jupytext .\hello_world.md --execute --to ipynb

    Get-ChildItem "**/*.md" | Foreach-Object { jupytext $_ --execute --to ipynb }


### Force sync

This synchronizes any changes to the markdown or ipynb file with each other (inputs only).

    jupytext .\hello_world.ipynb --execute

    Get-ChildItem "**/*.ipynb" | Foreach-Object { jupytext $_ --sync }


### Notebook header

All notebooks must have these magics in the header.
Plots use the *inline* backend so they can be generated programatically.
The figure size is blown up to `10,6`.

```
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```


### Testing notebooks

Powershell:

    ./test_notebooks.PS1

Bash:

    ./test_notebooks.sh
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Select a region of interest using FFTs

This example shows how to crop an image using Nanomesh using the characteristics of an FFT to select a bbox matching the lattice of the crystal.

This example shows how to crop an image using sample data from `nanomesh.data`.

If you want to use your own data, any numpy array can be passed to into an `Image` object. This will create a [`Volume`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) or [`Plane`](https://nanomesh.readthedocs.io/en/latest/nanomesh.plane.html#nanomesh.plane.Plane) object depending on the data dimensions. Data stored as `.npy` can be loaded using `Image.load()`.

```python
from nanomesh import Image
from nanomesh.data import nanopores3d

data = nanopores3d()

vol = Image(data)
plane = vol.select_plane(x=161)
```

This operation takes the fft of the image, and creates a regular array of the recurring components. The distance between the peaks in the image correspond to the distance between the pores in the source image.

```python
import numpy as np
import matplotlib.pyplot as plt


def abs2(x):
    return x.real**2 + x.imag**2


fft = np.fft.fft2(plane.image)
ifft = np.fft.ifft2(abs2(fft))

filtered = abs(ifft)

plt.figure()
plt.imshow(filtered)
plt.show()
```

### Peak-finding

Use a Difference of Gaussian to find the peaks in the image.

```python
from skimage import feature

peaks = feature.blob_dog(filtered,
                         min_sigma=10,
                         max_sigma=20,
                         overlap=1,
                         threshold=5)
peaks = peaks[:, 0:2]
x, y = peaks.T

plt.figure()
plt.imshow(filtered)
plt.scatter(y, x, color='red')
plt.show()
```

### ROI picking

A Delauney triangulation is used to create a mesh out of the peaks. The mesh is subdivided once to create additional granularity for the point picking in the roi selection.

```python
from scipy.spatial import Delaunay
from nanomesh import TriangleMesh

tris = Delaunay(peaks)
triangles = TriangleMesh.from_scipy(tris)
```

The subdivision uses [trimesh](https://trimsh.org/trimesh.remesh.html#trimesh.remesh.subdivide).

```python
trimesh = triangles.to_trimesh()
trimesh_subdivided = trimesh.subdivide()
triangles = TriangleMesh.from_trimesh(trimesh_subdivided)
```

The vertices are passed to the `.select_roi` method to pick from.

- By passing `from_points`, vertices snap to the nearest point (use 'ctrl' to drag it away)
- Press the 'esc' key to start a new polygon
- Hold the 'shift' key to move all of the vertices
- Hold the 'ctrl' key to move a single vertex

```python
roi = plane.select_roi(from_points=triangles.points)
```

The `.bbox` attribute is updated when the selection above changes.

```python
roi.bbox
```

Use the `.crop_to_roi` method to extract the region of interest.

```python
plane_roi = plane.crop_to_roi(bbox=roi.bbox)
plane_roi.show()
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Generate a 2D triangular mesh

This notebook shows how to mesh a 2D image:

1. Load and visualize a volume
2. Select a plane from a volume
3. Apply image filters and segment an image
4. Generate a 2D triangle mesh
5. Visualize and export the mesh to other formats


### Load and vizualize the data

This example uses nanopore sample data from `nanomesh.data`.

If you want to use your own data, any numpy array can be passed to into a [`Volume`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) object. Data stored as `.npy` can be loaded using `Volume.load()`.

```python
from nanomesh import Image
from nanomesh.data import nanopores3d

data = nanopores3d()

vol = Image(data)
vol.show_slice()
```

Nanomesh makes use of [`itkwidgets`](https://github.com/InsightSoftwareConsortium/itkwidgets) to render the volumes.

```python
vol.show()
```

There is also a simple `matplotlib` based viewer to plot the slices of the volume. The sliders can be used to select different directions and slice numbers. Optionally, the index of the first slice can be specified directly using, for example, `x=123`.

```python
vol.show_slice(x=123)
```

### Select plane from volume

Select single plane from the volume using the `.select_volume` method. In this case, we select slice #161 along the x-axis. The slice is loaded into a [`Plane`](https://nanomesh.readthedocs.io/en/latest/nanomesh.image.html#nanomesh.image.Plane) object.

```python
plane = vol.select_plane(x=161)
plane.show()
```

### Filter and segment the data

Image segmentation is a way to label the pixels of different regions of interest in an image. In this example, we are interested in separating the bulk material (Si) from the nanopores. In the image, the Si is bright, and the pores are dark.

First, we apply a [`gaussian filter`](https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian) to smooth out some of the image noise to get a cleaner segmentation.

**Note**: The code below is essentially short-hand for `plane_gauss = plane.apply(skimage.filters.gaussian, sigma=5)`. `apply` can be used for any operation that works on an array, i.e. from `numpy`, `scipy.ndimage` or `scikit-image`. If the output is a numpy array of the same dimensions, a `Image` object is returned. Anything else, and it will return the direct result of the function.

```python
plane_gauss = plane.gaussian(sigma=5)
plane_gauss.show()
```

Use the `compare_with_other` method to check out the difference:

```python
plane.compare_with_other(plane_gauss)
```

`scikit-image` contains a useful function to [try all threshold finders on the data](https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.try_all_threshold). These methods analyse the contrast histogram and try to find the optimal value to separate which parts of the image belong to each domain. The method below is a shortcut to this function.

Note that it is always possible to define your own threshold.

```python
plane_gauss.try_all_threshold(figsize=(5, 10))
```

We continue with the `li` method, because it gives a result with nice separation.

```python
thresh = plane_gauss.threshold('li')
thresh
```

Check how the segmented image compares to the original.

```python
segmented = plane_gauss.digitize(bins=[thresh])
plane.compare_with_digitized(segmented)
```

### Generate mesh


Meshes are generated using the `Mesher2D` class. Meshing consists of two steps:

1. Contour finding (using the [`find_contours`](https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.find_contours) function
2. Triangulation (using the [`triangle`](https://rufat.be/triangle/) library)

Contour finding uses the [marching cubes algorithm](https://en.wikipedia.org/wiki/Marching_cubes) to wrap all the pores in a polygon. `max_contour_dist=5` splits up long edges in the contour, so that no two points are further than 5 pixels apart. `level` is directly passed to `find_contours` and specifies the level at which the contour is generated. In this case, we set it to the threshold value determined above.

```python tags=[]
from nanomesh import Mesher2D

mesher = Mesher2D(plane_gauss)
mesher.generate_contour(max_contour_dist=5, level=thresh)

mesher.plot_contour()
```

The next step is to use the contours to initialize triangulation.

Triangulation options can be specified through the `opts` keyword argument. This example uses `q30` to generate a quality mesh with angles > 30°, and `a100` to set a maximum triangle size of 100 pixels. For more options, see [here](https://rufat.be/triangle/API.html#triangle.triangulate).

```python
mesh = mesher.triangulate(opts='q30a100')
```

Triangulation returns a `MeshContainer` dataclass that can be used for various operations, for example comparing it with the original image:

```python
plane_gauss.compare_with_mesh(mesh)
```

Or, showing the an interactive plot using pyvista:

(Use `.plot_itk()` for an interactive view)

```python
mesh.plot_pyvista(jupyter_backend='static', show_edges=True)
```

### Field data

Field data can be used to associate names with the values in the cell data. These are shown in the legend of mesh data (i.e. in the plots above). The field data is stored in the `.field_data` attribute. Because the data are somewhat difficult to use in this state, the properties `.field_to_number` and `.number_to_field` can be used to access the mapping per cell type.

```python
mesh.number_to_field
```

To update the values, you can update `.field_data` directory, or use `.set_field_data`. Note that field names are shared between cell types. For example, to relabel the cells data:



```python
mesh.set_field_data('triangle', {1: 'Silicon', 2: 'Pore'})
```

Plotting the mesh now shows the fields in the legend. Note that the fields are also saved when exported to a format that supports them (e.g. gmsh).

```python
mesh.plot(lw=1, color_map={0: 'lightgray'})
```

### Interoperability

The `MeshContainer` object can also be used to convert to various other library formats, such as:

- [`trimesh.Trimesh`](https://trimsh.org/trimesh.base.html#trimesh.base.Trimesh)
- [`pyvista.UnstructuredGrid`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)
- [`meshio.Mesh`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)

First, we must extract the triangle data:

```python
triangle_mesh = mesh.get('triangle')

pv_mesh = triangle_mesh.to_pyvista_unstructured_grid()
trimesh_mesh = triangle_mesh.to_trimesh()
meshio_mesh = triangle_mesh.to_meshio()
```

To save the data, use the `.write` method. This is essentially a very thin wrapper around `meshio`, equivalent to `meshio.write(...)`.

```python
mesh.write('out.msh', file_format='gmsh22', binary=False)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Segment image data using local thresholds

This notebook meshes an image with multiple domains using Nanomesh. The image contains a gradient from left to right.  If the image background is relatively uniform, global thresholds can be used to separate the different domains. In this case, the image has a contrast gradient from left to right, so this example shows three different methods to deal with this.

The image is then meshed using Nanomesh to obtain a mesh with multiple domains. The mesh contains triangles labeled as 'pore' or 'bulk' material.


### Loading and pre-processing the data

This example uses nanopore sample data from `nanomesh.data`. Notice how the image has a bit of a gradient going from top to bottom. We apply a gaussian filter to reduce image noise.

If you want to use your own data, any numpy array can be passed to into a [`Image`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) object. Data stored as `.npy` can be loaded using `Image.load()`.

```python
from nanomesh import Image
from nanomesh.data import nanopores_gradient

data = nanopores_gradient()

plane = Image(data).gaussian(sigma=5)
plane.show()
```

Using `.try_all_thresholds` is usually a good way to find a useful value to segment the data. In this case, the gradient prevents us from getting a useful segmentation.

```python
plane.try_all_threshold(figsize=(4, 10))
```

### Local thresholding

This section explores three different methods for local thresholding.

- User-defined local threshold ([api][1], [example][2])
- Otsu local threshold ([api][3], [example][4])
- Adaptive histogram equalization ([api][5], [example][6])

#### User-defined local threshold

In a local threshold filter, the local neighbourhoud is used to define a threshold map for the image. The `block_size` defines the size of the neighbourhood, a 101 by 101 pixel window in this case. The offset is used to tune the map. The advantage of this method is that it is fairly simple. The downside is that it requires manually tuning two parameters (`offset`, `blocksize`) to get a useful result.

[1]: https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.threshold_local
[2]: https://scikit-image.org/docs/stable/auto_examples/applications/plot_thresholding.html#local-thresholding
[3]: https://scikit-image.org/docs/stable/api/skimage.filters.rank.html#skimage.filters.rank.otsu
[4]: https://scikit-image.org/docs/stable/auto_examples/applications/plot_rank_filters.html#image-threshold
[5]: https://scikit-image.org/docs/stable/api/skimage.exposure.html#equalize-adapthist
[6]: https://scikit-image.org/docs/stable/auto_examples/color_exposure/plot_local_equalize.html#local-histogram-equalization

```python
import matplotlib.pyplot as plt

offset = 150
block_size = 101

local_thresh = plane.threshold('local', block_size=block_size, offset=offset)
seg_local = plane > local_thresh

## Plot

import matplotlib.pyplot as plt

fig, axes = plt.subplots(ncols=2, figsize=(6, 6))
ax = axes.ravel()

ax[0].imshow(plane.image)
ax[0].set_title('Original')

ax[1].imshow(seg_local.image)
ax[1].set_title('Segmentation')

for a in ax:
    a.axis('off')

plt.show()
```

#### Otsu local threshold

Applies an otsu rank filter to determine the local otsu threshold to segment on. Similar to the global otsu method, but constrained to a small area around each pixel. The advantage of this method that it only requires tuning a single parameter, namely the `radius` of the neighbourhood. The downside is that only the otsu method is available as a rank filter.

```python
from skimage.morphology import disk
from skimage.filters import rank
from skimage.util import img_as_ubyte

plane_norm = plane.normalize_values().apply(img_as_ubyte)

radius = 41
selem = disk(radius)

local_otsu = plane_norm.apply(rank.otsu, selem=selem)
seg_local_otsu = plane_norm >= local_otsu

## Plot

import matplotlib.pyplot as plt

fig, axes = plt.subplots(ncols=3, figsize=(9, 6))
axes = axes.ravel()

axes[0].imshow(plane.image)
axes[0].set_title('Original')

axes[1].imshow(local_otsu.image)
axes[1].set_title(f'Local Otsu filter\n(radius={radius})')

axes[2].imshow(seg_local_otsu.image)
axes[2].set_title('Segmentation')

for ax in axes:
    ax.axis('off')

plt.show()
```

#### Adaptive histogram equalization

This method tries to correct the image by removing the local gradient first, using adaptive histogram equalization. The advantage is that all global threshold finders are available. Another advantage is that there are no parameters to tune. At the same time, this is also a disadvantage if the result is not good 😉

```python
from skimage import exposure

plane_eq = plane.normalize_values()
plane_eq = plane_eq.apply(exposure.equalize_adapthist)

plane_eq.try_all_threshold(figsize=(6, 15))
```

### Compare results

The next cell creates a plot that compares the result of all three methods. For the histogram equalization, the *Otsu* and *Li* filters are shown.

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(ncols=5, figsize=(9, 5))
axes = axes.ravel()

axes[0].imshow(plane.image)
axes[0].set_title('Original')

axes[1].imshow(seg_local.image)
axes[1].set_title('Local\nthresholding')

axes[2].imshow(seg_local_otsu.image)
axes[2].set_title('Local Otsu')

seg_clahe_otsu = plane_eq.binary_digitize(threshold='otsu')

axes[3].imshow(seg_clahe_otsu.image)
axes[3].set_title('Hist. equal.\n + global Otsu')

seg_clahe_li = plane_eq.binary_digitize(threshold='li')

axes[4].imshow(seg_clahe_li.image)
axes[4].set_title('Hist. equal.\n + global Li')

for ax in axes:
    ax.axis('off')

plt.show()
```

### Meshing the image data

Finally, it's time to mesh the image.

```python
# Triangulation has some issues with small elements near the border,
# so these will be cleared first.
# https://github.com/hpgem/nanomesh/issues/86

seg_local = seg_local.clear_border(object_label=0, fill_val=1)
seg_local.show()
```

```python
%%time

from nanomesh import Mesher2D

mesher = Mesher2D(seg_local)
mesher.generate_contour(max_contour_dist=4)

mesh = mesher.triangulate(opts='q30a100')
```

View the result using `matplotlib`:

```python
mesh.plot(color_map={0: 'lightgray'}, lw=1)
```

Or, view the result using pyvista:

(Note: Use `mesh.plot_itk()` for an interactive view).

```python
mesh.plot_pyvista(jupyter_backend='static', show_edges=False)
```

Save the data:

```python
mesh.write("mesh_gradient.msh", file_format='gmsh22', binary=False)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Generate a 3D tetrahedral mesh

This notebook shows how to mesh a 3D volume:

1. Load and visualize a volume
2. Apply image filters and segment image
3. Generate a 3D surface mesh
4. Visualize and export the mesh to other formats

```python
import pyvista as pv
from skimage import filters
import numpy as np
```

### Load and vizualize the data

This example uses nanopore sample data from `nanomesh.data`.

If you want to use your own data, any numpy array can be passed to into a [`Image`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) object. Data stored as `.npy` can be loaded using `Image.load()`.

```python
from nanomesh import Volume
from nanomesh.data import nanopores3d

data = nanopores3d()

vol = Volume(data)
vol.show_slice()
```

For this example, select a subvolume using `.select_subvolume` and downscale the image to keep the cpu times in check.

```python
from skimage.transform import rescale

subvol = vol.select_subvolume(
    ys=(0, 100),
    xs=(0, 100),
).apply(rescale, scale=0.5)
subvol.show_slice()
```

Nanomesh makes use of [`itkwidgets`](https://github.com/InsightSoftwareConsortium/itkwidgets) to render the volumes.

```python
subvol.show()
```

### Filter and segment the data

Image segmentation is a way to label the pixels of different regions of interest in an image. In this example, we are interested in separating the bulk material (Si) from the nanopores. In the image, the Si is bright, and the pores are dark.

First, we apply a [`gaussian filter`](https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian) to smooth out some of the image noise to get a cleaner segmentation.

```python
subvol_gauss = subvol.gaussian(sigma=1)
subvol_gauss.show_slice(x=12)
```

`scikit-image` contains a useful function to [try all threshold finders on the data](https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.try_all_threshold). These methods analyse the contrast histogram and try to find the optimal value to separate which parts of the image belong to each domain.

Since the function only works on a single slice, we first select a slice using the `.select_plane` method.

```python
from skimage import filters

plane = subvol_gauss.select_plane(x=12)
plane.try_all_threshold(figsize=(5, 10))
```

We will use the `li` method, because it gives nice separation.

The threshold value is used to segment the image using [`np.digitize`](https://numpy.org/doc/stable/reference/generated/numpy.digitize.html#numpy-digitize).

```python
subvol_seg = subvol_gauss.binary_digitize(threshold='minimum')
subvol_seg = subvol_seg.invert_contrast()
subvol_seg.show_slice()
```

### Generate 3d tetragonal mesh

Meshes can be generated using the `Mesher` class. Meshing consists of two steps:

1. Contour finding (using the [`marching_cubes`](https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.marching_cubes) function
2. Triangulation (using the [`tetgen`](https://tetgen.pyvista.org/) library)

`Mesher` requires a segmented image. `generate_contour()` wraps all domains of the image corresponding to that label. Here, 1 corresponds to the bulk (Si) material.

Meshing options are defined in the [tetgen documentation](http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html#sec35). These can be specified using the `opts` parameter. The default options are `opts='-pAq1.2`:

- `-A`: Assigns attributes to tetrahedra in different regions.
- `-p`: Tetrahedralizes a piecewise linear complex (PLC).
- `-q`: Refines mesh (to improve mesh quality).

Also useful:

- `-a`: Applies a maximum tetrahedron volume constraint. Don't make `-a` too small, or the algorithm will take a very long time to complete. If this parameter is left out, the triangles will keep growing without limit.

```python
%%time

from nanomesh import Mesher

mesher = Mesher(subvol_seg)
mesher.generate_contour()
mesh = mesher.tetrahedralize(opts='-pAq')
```

Tetrahedralization returns a `TetraMesh` dataclass that can be used for various operations, for example showing the result using `itkwidgets`:

```python
mesh.plot_pyvista(jupyter_backend='static', show_edges=True)
```

### Using region markers

By default, the region attributes are assigned automatically by `tetgen`. Tetrahedra in each enclosed region will be assigned a new label sequentially.

Region markers are used to assign attributes to tetrahedra in different regions. After tetrahedralization, the region markers will 'flood' the regions up to the defined boundaries. The elements of the resulting mesh are marked according to the region they belong to (`tetras.metadata['tetgenRef']`.

You can view the existing region markers by looking at the `.region_markers` attribute on the contour.

```python
mesher.contour.region_markers
```

It is possible to set your own attributes using the `region_markers` parameter. These are a list of `RegionMarker` objects stored in a `RegionMarkerList` container. You can define your own, or you can use the methods on the `mesher.contour.region_markers` attribute.

For example, to relabel the pores sequentially:

```python
mesher.contour.region_markers = mesher.contour.region_markers.label_sequentially(
    2, fmt_name='pore{}')
mesher.contour.region_markers
```

Then re-run tetrahedralization.

```python
%%time
import numpy as np

mesh = mesher.tetrahedralize(opts='-pAq')

for label in mesher.contour.region_markers.labels:
    num = np.sum(mesh.cell_data['tetgen:ref'] == label)
    print(f'{num} tetrahedra with attribute `{label}`')
```

Note that all connected regions have a different label now.

```python
mesh.plot_pyvista(jupyter_backend='static', show_edges=True)
```

### Mesh evaluation

The mesh can be evaluated using the `metrics` module. This example shows how to calculate all metrics and plot them on a section through the generated mesh.

```python
from nanomesh import metrics

tetra_mesh = mesh.get('tetra')

metrics_dict = metrics.calculate_all_metrics(tetra_mesh, inplace=True)
metrics_dict
```

Using the `.plot_submesh()` method, any array that is present in the metadata can be plotted. `plot_submesh()` is flexible, in that it can show a slice through the mesh as defined using `index`, `along`, and `invert`. Extra keyword arguments, such as `show_edges` and `lighting` are passed on to [`Plotter.add_mesh()`](https://docs.pyvista.org/api/plotting/_autosummary/pyvista.Plotter.add_mesh.html?highlight=add_mesh).

```python
tetra_mesh.plot_submesh(
    along='x',
    index=15,
    scalars='min_angle',
    show_edges=True,
    lighting=True,
    backend='static',
)
```

### Interoperability

The `TetraMesh` object can also be used to convert to various other library formats, such as:

- [`trimesh.open3d`](http://www.open3d.org/docs/release/python_api/open3d.geometry.TetraMesh.html#open3d.geometry.TetraMesh)
- [`pyvista.UnstructuredGrid`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)
- [`meshio.Mesh`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)


To save the data, use the `.write` method. This is essentially a very thin wrapper around `meshio`, equivalent to `meshio_mesh.write(...)`.

```python
tetra_mesh.write('volume_mesh.msh', file_format='gmsh22', binary=False)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Select a region of interest

This example shows how to crop an image using sample data from `nanomesh.data`.

If you want to use your own data, any numpy array can be passed to into a [`Volume`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) object. Data stored as `.npy` can be loaded using `Volume.load()`.

```python
from nanomesh import Volume
from nanomesh.data import nanopores3d

data = nanopores3d()

vol = Volume(data)
plane = vol.select_plane(x=161)
```

Use `.select_roi` to open an interactive widget. Select a region of interest in the figure by enclosing them within a polygon. A rectangle is fitted to the polygon.

- Press the 'esc' key to start a new polygon
- Hold the 'shift' key to move all of the vertices
- Hold the 'ctrl' key to move a single vertex

```python
roi = plane.select_roi()
```

The `.bbox` attribute is updated when the selection above changes.

```python
roi.bbox
```

Use the `.crop_to_roi` method to extract the region of interest.

```python
plane_roi = plane.crop_to_roi(bbox=roi.bbox)
plane_roi.show()
```

### Manual cropping

Alternatively, if you know which points to extract from the  to

Use the `minimum_bounding_rectangle` function to fit the smallest rectangle around the given points.

The example below demonstrates this on a slanted square. Well, almost, because the coordinates contain a small error to demonstrate the function.

```python
from nanomesh.image import minimum_bounding_rectangle
import numpy as np

roi = np.array([
    [60, 60],
    [110, 110],
    [60, 145],  # <- should be 150
    [10, 110],
])

bbox = minimum_bounding_rectangle(roi)
bbox
```

The bounding box (`bbox`) can then be used to extract this region from the image using `extract_rectangle`. If the bounding box is slanted, the image will be straightened using an Euclidean transform.

```python
from nanomesh.image import extract_rectangle

cropped = plane.apply(extract_rectangle, bbox=bbox)
cropped.show()
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Calculate mesh quality metrics

This notebook reads a mesh and plots different quality indicators:

- Minimum/maximum angle
- Ratio min/max edge length
- Ratio circumscribed to inscribed circle (largest circle fitting inside vs smallest circle fitting around a triangle)

The indicators are plotted on the mesh and as a histogram.

You can use your own mesh by supplying it to `MeshContainer.read()`.

```python
from nanomesh import metrics
from nanomesh import MeshContainer
```

```python
mesh = MeshContainer.read('out.msh')
triangle_mesh = mesh.get('triangle')

triangle_mesh.plot()
```

### Metrics

Quality metrics are available through the `metrics` submodule, for example to access the area for each face:

```python
metrics.area(triangle_mesh)
```

### Minumum and maximum cell angles

`nanomesh.metrics` includes convenience functions for plotting histograms and colored 2d meshes. The `ax` object can be re-used to overlay the mesh triangles.

```python
plot_kwargs = {
    'linewidth': 1,
    'show_labels': ('Pore', ),
    'colors': ('tab:orange', ),
    'flip_xy': False,
    'legend': 'all',
}

metrics.histogram(triangle_mesh, metric='min_angle')
ax = metrics.plot2d(triangle_mesh, metric='min_angle')
triangle_mesh.plot_mpl(ax, **plot_kwargs)
```

```python
metrics.histogram(triangle_mesh, metric='max_angle')
ax = metrics.plot2d(triangle_mesh, metric='max_angle')
triangle_mesh.plot_mpl(ax, **plot_kwargs)
```

### Ratio between radii

Another useful metric is the ratio between the inner and outer radius. For more info, see this [link](https://www.geogebra.org/m/VRE3Dyrn).

```python
metrics.histogram(triangle_mesh, metric='radius_ratio')
ax = metrics.plot2d(triangle_mesh, metric='radius_ratio')
triangle_mesh.plot_mpl(ax, **plot_kwargs)
```

### Ratio between longest and shortest edge

```python
metrics.histogram(triangle_mesh, metric='max_min_edge_ratio')
ax = metrics.plot2d(triangle_mesh, metric='max_min_edge_ratio')
triangle_mesh.plot_mpl(ax, **plot_kwargs)
```

### Calculate and export all metrics

This way they can be viewed in another program like Paraview.

```python
metrics.calculate_all_metrics(triangle_mesh, inplace=True)
triangle_mesh.write("mesh_quality.msh", file_format='gmsh22', binary=False)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Examine nanopore properties

This notebook serves as a demonstration to examine the nanopore properties of a slice through a crystal. It shows how to use [scikit-image](https://scikit-image.org) [regionprops](https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops) to describe the size and shape of the pores.

```python
from nanomesh import Image
import numpy as np

plane_orig = Image.load('x500.npy')

# rotate to make better use of space
plane_orig = plane_orig.apply(np.rot90)

# smooth image for better segmentation
plane = plane_orig.gaussian(sigma=5)
plane.show()
```

### Image segmentation

The image has a slight gradient from top to bottom, therefore we need to apply local thresholding. This has been described in [another notebook](https://nanomesh.readthedocs.io/en/latest/examples/How%20to%20segment%20image%20data%20using%20local%20thresholds.html).

```python
from skimage.filters import threshold_local

offset = 150
block_size = 101

local_thresh = plane.threshold('local', block_size=block_size, offset=offset)
seg_local = plane.image > local_thresh.image

# invert contrast for object detection
seg = Image(1 - seg_local)
seg.show()
```

### Obtain regionprops

First the segmented image must be labeled. This means that all the objects in the above image are given a different label. The regionprops function then calculates the properties for each label. The original image is passed so that the contrast can be used for intensity calculations (if needed).

```python
from skimage import measure

labels = measure.label(seg.image)

props = measure.regionprops(labels, plane_orig.image)
```

Have a look at [the documentation](https://scikit-image.org/docs/dev/api/skimage.measure.html?highlight=regionprops#regionprops) for the available properties.

Below are the properties we are interested in for the plots.

```python
properties = ['area', 'eccentricity', 'perimeter', 'mean_intensity']
```

### Plots

The example below is adapted from the the [scikit-image gallery](https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_regionprops.html). It interactively plots the selected objects on the source image, so that the properties of each can be explored. This example requires `plotly` and `pandas` to be installed (`pip install plotly pandas`).

```python
import plotly
import plotly.express as px
import plotly.graph_objects as go

fig = px.imshow(plane_orig.image, binary_string=True)
fig.update_traces(hoverinfo='skip')  # hover is only for label info

# For each label, add a filled scatter trace for its contour,
# and display the properties of the label in the hover of this trace.
for index in range(1, labels.max()):
    label_i = props[index].label
    contour = measure.find_contours(labels == label_i, 0.5)[0]
    y, x = contour.T
    hoverinfo = ''
    for prop_name in properties:
        hoverinfo += f'<b>{prop_name}: {getattr(props[index], prop_name):.2f}</b><br>'
    fig.add_trace(
        go.Scatter(x=x,
                   y=y,
                   name=label_i,
                   mode='lines',
                   fill='toself',
                   showlegend=False,
                   hovertemplate=hoverinfo,
                   hoveron='points+fills'))

plotly.io.show(fig)
```

Additionally, the properties can be used to make some distribution plots.

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(nrows=2, ncols=2)

axes = axes.flatten()

for ax, prop_name in zip(axes, properties):
    data = [getattr(prop, prop_name) for prop in props]
    ax.hist(data, bins=25, rwidth=0.9)
    ax.set_title(prop_name)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Memory-map a large data set

This example demonstrates how to load and process a large data file (8 GB). This file is too large to load into memory at once (on most PCs), so memory mapping is used to only load the part of the data. A small volume and a single slice are extracted and saved for use in the other examples.

```python
import numpy as np
from nanomesh import Image
```

Load the data using `mmap_mode='r'`. This is a shallow interface to `np.memmap` that prevents the entire data set to be loaded into memory at once.

```python
data_name = 'G:\escience\hpgem\PG_EBAR_18072017_CH_6C_s15_10nm_rec_sa7_1024_1024_2048.vol'

vol = Image.load(data_name, mmap_mode='r')
vol
```

```python
Image.load?
```

The `show_slice` method is still quite responsive to look at sections of the data. `show_volume` also works, but loads the entire volume into memory, which may make everything a bit slow and unresponsive 😅

```python
vol.show_slice(index=500)
```

It's easier to work with a section of the data. Note that a `np.memmap` object can be sliced like a normal `numpy` array, so we can extract a subvolume to work with:

```python
cropped = vol.select_subvolume(xs=(450, 550), ys=(150, 725), zs=(410, 1470))
cropped.show_slice()
```

Display the cropped data.

```python
cropped.show()
```

Save the data to numpy binary format.

```python
cropped.save('slab_x450-550.npy')
```

Select a slice from the data and trim the edges for further analyses

```python
plane = vol.select_plane(x=500)
plane = plane.crop(left=150, right=850, top=655, bottom=845)
plane.show()
```

And save it...

```python
plane.save('nanopores_gradient.npy')
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Generate a tetrahedral mesh from 3D cells data

This notebook shows how to mesh a 3D volume:

1. Load and visualize a volume
2. Apply image filters and segment image
3. Generate a 3D surface mesh
4. Visualize and export the mesh to other formats

```python
import pyvista as pv
from skimage import filters, data
import numpy as np
```

### Load and vizualize the data

This example uses the `cells3d` sample data from [scikit-image](https://scikit-image.org/docs/dev/api/skimage.data.html#cells3d).

If you want to use your own data, any numpy array can be passed to a [`Image`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) object. Data stored as `.npy` can be loaded using `Image.load()`.

```python
from nanomesh import Image
import numpy as np

from skimage.data import cells3d

data = cells3d()[:, 1, :, :]

vol = Image(data)

# swap axes so depth matches z
vol = vol.apply(np.swapaxes, axis1=0, axis2=2)

vol.show_slice()
```

For this example, select a subvolume using `.select_subvolume` to keep the cpu times in check.

```python
from skimage.transform import rescale

subvol = vol.select_subvolume(
    ys=(0, 128),
    zs=(128, -1),
)
subvol.show_slice()
```

Nanomesh makes use of [`itkwidgets`](https://github.com/InsightSoftwareConsortium/itkwidgets) to render the volumes.

```python
subvol.show()
```

### Filter and segment the data

Image segmentation is a way to label the pixels of different regions of interest in an image. In this example, we are interested in separating the bulk material (Si) from the nanopores. In the image, the Si is bright, and the pores are dark.

First, we apply a [`gaussian filter`](https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian) to smooth out some of the image noise to get a cleaner segmentation.

```python
subvol_gauss = subvol.gaussian(sigma=1)
subvol_gauss.show_slice()
```

`scikit-image` contains a useful function to [try all threshold finders on the data](https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.try_all_threshold). These methods analyse the contrast histogram and try to find the optimal value to separate which parts of the image belong to each domain.

Since the function only works on a single slice, we first select a slice using the `.select_plane` method.

```python
from skimage import filters

plane = subvol_gauss.select_plane(x=30)
plane.try_all_threshold(figsize=(5, 10))
```

We will use the `yen` method, because it gives nice separation.

The threshold value is used to segment the image using [`np.digitize`](https://numpy.org/doc/stable/reference/generated/numpy.digitize.html#numpy-digitize).

```python
subvol_gauss.binary_digitize(threshold='yen')
```

```python
subvol_seg = subvol_gauss.binary_digitize(threshold='yen')
subvol_seg.show_slice()
```

```python
from scipy import ndimage


def fill_holes(image):
    return ndimage.binary_fill_holes(image).astype(int)


subvol_seg = subvol_seg.apply(fill_holes)
subvol_seg.show_slice()
```

### Generate 3d tetragonal mesh

Meshes can be generated using the `Mesher` class. Meshing consists of two steps:

1. Contour finding (using the [`marching_cubes`](https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.marching_cubes) function
2. Triangulation (using the [`tetgen`](https://tetgen.pyvista.org/) library)

`Mesher` requires a segmented image. `generate_contour()` wraps all domains of the image corresponding to that label. Here, 1 corresponds to the cells.

Meshing options are defined in the [tetgen documentation](http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html#sec35). These can be specified using the `opts` parameter. The default options are `opts='-pAq1.2`:

- `-A`: Assigns attributes to tetrahedra in different regions.
- `-p`: Tetrahedralizes a piecewise linear complex (PLC).
- `-q`: Refines mesh (to improve mesh quality).

Also useful:

- `-a`: Applies a maximum tetrahedron volume constraint. Don't make `-a` too small, or the algorithm will take a very long time to complete. If this parameter is left out, the triangles will keep growing without limit.

```python
%%time

from nanomesh import Mesher

mesher = Mesher(subvol_seg)
mesher.generate_contour()
mesh = mesher.tetrahedralize(opts='-pAq')
```

Tetrahedralization returns a `TetraMesh` dataclass that can be used for various operations, for example showing the result using `itkwidgets`:

```python
mesh.plot_pyvista(jupyter_backend='static', show_edges=True)
```

### Using region markers

By default, the region attributes are assigned automatically by `tetgen`. Tetrahedra in each enclosed region will be assigned a new label sequentially.

Region markers are used to assign attributes to tetrahedra in different regions. After tetrahedralization, the region markers will 'flood' the regions up to the defined boundaries. The elements of the resulting mesh are marked according to the region they belong to (`tetras.metadata['tetgenRef']`.

You can view the existing region markers by looking at the `.region_markers` attribute on the contour.

```python
mesher.contour.region_markers
```

### Mesh evaluation

The mesh can be evaluated using the `metrics` module. This example shows how to calculate all metrics and plot them on a section through the generated mesh.

```python
from nanomesh import metrics

tetra_mesh = mesh.get('tetra')

metrics_dict = metrics.calculate_all_metrics(tetra_mesh, inplace=True)
metrics_dict
```

Using the `.plot_submesh()` method, any array that is present in the metadata can be plotted. `plot_submesh()` is flexible, in that it can show a slice through the mesh as defined using `index`, `along`, and `invert`. Extra keyword arguments, such as `show_edges` and `lighting` are passed on to [`Plotter.add_mesh()`](https://docs.pyvista.org/api/plotting/_autosummary/pyvista.Plotter.add_mesh.html?highlight=add_mesh).

```python
tetra_mesh.plot_submesh(
    along='x',
    index=30,
    scalars='min_angle',
    show_edges=True,
    lighting=True,
    backend='static',
)
```

### Interoperability

The `TetraMesh` object can also be used to convert to various other library formats, such as:

- [`trimesh.open3d`](http://www.open3d.org/docs/release/python_api/open3d.geometry.TetraMesh.html#open3d.geometry.TetraMesh)
- [`pyvista.UnstructuredGrid`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)
- [`meshio.Mesh`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)


To save the data, use the `.write` method. This is essentially a very thin wrapper around `meshio`, equivalent to `meshio_mesh.write(...)`.

```python
tetra_mesh.write('cells.vtk')
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Multi-domain mesh from simulated data

This notebook demonstrates an example using [tetgen](https://wias-berlin.de/software/tetgen/) as the tetrahedralizer starting from a data volume.

Tetgen file formats:
http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual006.html#ff_poly

```python
import numpy as np
from nanomesh import Image
```

### Generate data

The cell below generates a data volume with feature blobs. It consists of two domains, 1 for the *bulk* features, and 0 for the background.

If you want to use your own data, any numpy array can be passed to a [`Image`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) object. Data stored as `.npy` can be loaded using `Image.load()`.

```python
from nanomesh.data import binary_blobs3d

data = binary_blobs3d(seed=42)

vol = Image(data)
vol.show_slice()
```

### Meshing

Mesh generation consists of two steps.

1. Generate the contour mesh
2. Generate the volume envelope

The contour mesh is the surface separating the domains in the data. In this case, the 1's and 0's, so the contour level does not need to be defined. The surface mesh is completed by wrapping the entire data volume in an envelope. This makes sure that the mesh (including any domain regions) is watertight, which is a requirement for generating a tetrahedral volume mesh.

```python
from nanomesh import Mesher

mesher = Mesher(vol)
mesher.generate_contour()
```

To plot the surface mesh:

(Use `.plot_itk()` for an interactive view)

```python
mesher.contour.plot_pyvista(jupyter_backend='static', show_edges=True)
```

### Generating a tetrahedral volume mesh

The volume mesh is generated using the `.tetrahedralize` method. This returns a tetrahedral mesh. Each domain, separated by the contour mesh defined above, is assigned a value.

The options used below:

- `-A`: Assigns attributes to tetrahedra in different regions.
- `-p`: Tetrahedralizes a piecewise linear complex (PLC).
- `-q`: Refines mesh (to improve mesh quality).
- `-a`: Applies a maximum tetrahedron volume constraint.

Don't make `-a` too small, or the algorithm will take a very long time to complete. If this parameter is left out, the triangles will keep growing without limit.

The region attributes are stored in the `tetgenRef` parameter.

Available options: http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html

```python
tetras = mesher.tetrahedralize(opts='-pAq -a10')
tetras.plot_pyvista(jupyter_backend='static',
                    show_edges=True)  # Use .plot_itk() for an interactive view
```

By default, the region attributes are assigned automatically by `tetgen`. Region markers assign attributes to tetrahedra in different regions. After tetrahedralization, the region markers will 'flood' the regions up to the defined boundaries. The elements of the resulting mesh are marked according to the region they belong to (`tetras.cell_data['tetgenRef']`.

It is possible to set your own attributes. The label corresponds to the attribute, and the value must be a point inside the region. The constraint can be used to set the maximum size of the triangles in combination with the `-a` parameter.

The next cell shows how you can update the region markers for the contour.

```python
for i, region_marker in enumerate(mesher.contour.region_markers):
    if region_marker.name == 'background':
        region_marker.constraint = 100
    else:
        region_marker.constraint = 1
        region_marker.label = i + 1
        region_marker.name = f'feature{i+1}'

mesh = mesher.tetrahedralize(opts='-pAq -a')

# Use .plot_itk() for an interactive view
mesh.plot_pyvista(jupyter_backend='static', show_edges=True)
```

```python
tetra_mesh = mesh.get('tetra')

for marker in mesher.contour.region_markers:
    num = np.sum(tetra_mesh.cell_data['tetgen-ref'] == marker.label)
    print(f'{num} tetrahedra with attribute `{marker}`')
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Generate a 2D triangular mesh

This notebook shows how to mesh a 2D image:

1. Load and visualize a 2D data image
2. Generate a 2D triangle mesh
3. Visualize the mesh
4. Relabel the regions
5. Export the mesh to other formats


### Load and vizualize the data

This example uses generated sample data from `nanomesh.data`.

If you want to use your own data, any (2D) numpy array can be passed to into a [`Image`](https://nanomesh.readthedocs.io/en/latest/nanomesh.plane.html#nanomesh.plane.Plane) object. Data stored as `.npy` can be loaded using `Image.load()`.

```python
from nanomesh import Image
from nanomesh.data import binary_blobs2d

data = binary_blobs2d(seed=12345)

plane = Image(data)
plane.show()
```

### Generate mesh


Meshes are generated using the `Mesher2D` class. Meshing consists of two steps:

1. Contour finding (using the [`find_contours`](https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.find_contours) function
2. Triangulation (using the [`triangle`](https://rufat.be/triangle/) library)

Contour finding uses the [marching cubes algorithm](https://en.wikipedia.org/wiki/Marching_cubes) to wrap all the pores in a polygon. `max_contour_dist=5` splits up long edges in the contour, so that no two points are further than 5 pixels apart. `level` is directly passed to `find_contours` and specifies the level at which the contour is generated. In this case, we set it to the threshold value determined above.

```python tags=[]
from nanomesh import Mesher2D

mesher = Mesher2D(plane)
mesher.generate_contour(max_contour_dist=3)

mesher.plot_contour()
```

The next step is to use the contours to initialize triangulation.

Triangulation options can be specified through the `opts` keyword argument. This example uses `q30` to generate a quality mesh with angles > 30°, and `a100` to set a maximum triangle size of 100 pixels. For more options, see [here](https://rufat.be/triangle/API.html#triangle.triangulate).

```python
mesh = mesher.triangulate(opts='q30a100')
```

Triangulation returns a `MeshContainer` dataclass that can be used for various operations, for example comparing it with the original image:

```python
plane.compare_with_mesh(mesh)
```

### Region markers

By default, regions are split up in the *background* and *features*. Feature regions are grouped by the label 1.

```python
mesher.contour.region_markers
```

To label regions sequentially, set `group_regions=False`:

```python
mesher = Mesher2D(plane)
mesher.generate_contour(max_contour_dist=3, group_regions=False)

mesher.plot_contour()
```

Notice that each feature has been given a unique name:

```python
mesher.contour.region_markers
```

These labels will assigned to each triangle in the corresponding region after triangulation. These are stored in `mesh.cell_data` of the `MeshContainer`. This container stores a single set of points, and both the line segments (`LineMesh`) and triangles (`TriangleMesh`). To extract the triangle cells only, use `MeshContainer.get('triangle')`. this returns a class that is simpler to work with.

The cell below shows how to use this to access the cell data for the triangle cells in the mesh.

```python
mesh = mesher.triangulate(opts='q30a100')
mesh.get('triangle').cell_data
```

### Field data

Field data can be used to associate names with the values in the cell data. These are shown in the legend of mesh data (i.e. in the plots above). The field data is stored in the `.field_data` attribute. Because the data are somewhat difficult to use in this state, the properties `.field_to_number` and `.number_to_field` can be used to access the mapping per cell type.

```python
mesh.number_to_field
```

To update the values, you can update `.field_data` directory, or use `.set_field_data`. Note that field names are shared between cell types. For example, to relabel the cells data:



```python
fields = {}
for k, v in mesh.number_to_field['triangle'].items():
    v = v.replace('background', 'Silicon')
    v = v.replace('feature', 'Pore')
    fields[k] = v

mesh.set_field_data('triangle', fields)
mesh.number_to_field
```

Plotting the mesh now shows the fields in the legend. Note that the fields are also saved when exported to a format that supports them (e.g. *gmsh*).

```python
mesh.plot(lw=1, color_map={0: 'lightgray'}, legend='floating')
```

### Interoperability

The `MeshContainer` object can also be used to convert to various other library formats, such as:

- [`trimesh.Trimesh`](https://trimsh.org/trimesh.base.html#trimesh.base.Trimesh)
- [`pyvista.UnstructuredGrid`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)
- [`meshio.Mesh`](https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html)

First, we must extract the triangle data:

```python
triangle_mesh = mesh.get('triangle')

pv_mesh = triangle_mesh.to_pyvista_unstructured_grid()
trimesh_mesh = triangle_mesh.to_trimesh()
meshio_mesh = triangle_mesh.to_meshio()
```

To save the data, use the `.write` method. This is essentially a very thin wrapper around `meshio`, equivalent to `meshio.write(...)`.

```python
mesh.write('out.msh', file_format='gmsh22', binary=False)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Pad a 3D mesh

This notebook shows how to pad a 3D mesh. Note that padding a 3D mesh is done on the contour (before tetrahedralization).

```python
from nanomesh import Image
from nanomesh.data import binary_blobs3d
```

### Generate some data

This cell generates a 3D data set with some feature blobs.

If you want to use your own data, any numpy array can be passed to into a [`Image`](https://nanomesh.readthedocs.io/en/latest/nanomesh.volume.html#nanomesh.volume.Volume) object. Data stored as `.npy` can be loaded using `Image.load()`.

```python
data = binary_blobs3d(seed=2020)

vol = Image(data)
vol.show_slice()
```

### Generating the contour

```python
from nanomesh import Mesher

mesher = Mesher(vol)
mesher.generate_contour()
mesher.plot_contour(jupyter_backend='static', show_edges=True)
```

### Padding different sides

The mesh can be padded using a similar API as 2d meshes. Each side (top/bottom, left/right, front/back) can be padded. A width must be specified.

Regions are labeled with a number. If no label is given, an arbitrary number is assigned. This is used to identify different regions in the mesh.

Note that tetgen will assign a different label for physically separate
regions, even when they are given the same label/name.

```python
mesher.pad_contour(side='left', width=5, name='inner pad')
mesher.pad_contour(side='left', width=10, name='outer pad')
mesher.pad_contour(side='right', width=5, name='inner pad')
mesher.pad_contour(side='right', width=10, name='outer pad')

mesher.pad_contour(side='front', width=8, name='front pad')
mesher.pad_contour(side='back', width=8, name='back pad')
mesher.plot_contour(jupyter_backend='static', show_edges=True)
```

### Generate tetrahedral mesh

Finally, generate the tetrahedral mesh. Notice that the inner and outer pads have the same label, because we assigned the same name in `Mesher.pad_contour()`.

```python
tetras = mesher.tetrahedralize(opts='-pAq -a10000')
tetras.plot_pyvista(jupyter_backend='static', show_edges=True)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Pad a 2D mesh

This notebook shows how to pad a 2D mesh. Note that padding a 2d mesh is done on the contour (before triangulation).


### Generating some data

This cell generates a simple 2d plane.

If you want to use your own data, any numpy array can be passed to into a [`Image`](https://nanomesh.readthedocs.io/en/latest/nanomesh.plane.html#nanomesh.volume.Plane) object. Data stored as `.npy` can be loaded using `Image.load()`.

```python
from nanomesh.data import binary_blobs2d
from nanomesh import Image

data = binary_blobs2d(length=100, seed=42)

plane = Image(data)
plane.show()
```

### Generating the contour

The first step in mesh generation is to generate the contour.

```python
from nanomesh import Mesher2D

mesher = Mesher2D(plane)
mesher.generate_contour()
mesher.plot_contour(legend='fields')
```

### Padding different sides

The mesh can be padded using a similar API as 3d meshes. Each side (top/bottom, left/right) can be padded. A width must be specified.

Regions are labeled with a number. If no label is given, an arbitrary number is assigned. This is used to identify different regions in the mesh.

Padded areas can be given a name. Regions with the same name are assigned the same number.

```python
mesher.pad_contour(side='left', width=30, name='Left side')
mesher.pad_contour(side='right', width=40, name='Right side')
mesher.pad_contour(side='top', width=20, label=11)
mesher.pad_contour(side='bottom', width=50, label=11)
mesher.plot_contour(legend='fields')
```

### Generate triagonal mesh

Finally, generate the triagonal mesh.

Note that the legend specifies the name of the region if available in the `.fields` attribute.

```python
mesh = mesher.triangulate(opts='pAq30a100e')
mesh.plot(legend='floating', hide_labels=(0, ), linewidth=1)
```

### Labelling outer boundaries

The outer boundaries can be labeled using the `LineMesh.label_boundaries` method.

The default `.cell_data` key is `'physical'`. This can be overridden using the `key='...'` parameter. To label the top and bottom boundaries, use the `top`, `bottom` parameters.

```python
line_mesh = mesh.get('line')

line_mesh.label_boundaries(left='outer left', right='outer right')

# transfer labels back to MeshContainer
mesh.set_cell_data('line', 'physical', line_mesh.cell_data['physical'])
mesh.set_field_data('line', line_mesh.number_to_field)

mesh.plot(legend='floating', hide_labels=(0, ), linewidth=1)
```

### Padding left / right sides

The width, mesh quality, and label assigned to this this area can be defined.

This example shows how to double pad the left and right sides with different triangle sizes for each step.

```python
mesher = Mesher2D(plane)
mesher.generate_contour()

mesher.pad_contour(side='left', width=20, label=30, name='inner pad')
mesher.pad_contour(side='left', width=40, label=40, name='outer pad')

mesher.pad_contour(side='right', width=20, label=30, name='inner pad')
mesher.pad_contour(side='right', width=40, label=40, name='outer pad')

padded_mesh = mesher.triangulate(opts='pAq30a100e')

padded_mesh.plot('triangle', legend='fields')
```

### Spiral mesh

This pattern is infinitely extensible. The example below shows the flexibility of the method.

```python
from itertools import cycle
import numpy as np

mesher = Mesher2D(plane)
mesher.generate_contour()

choices = ('left', 'bottom', 'right', 'top')

for i, side in zip(range(1, 50), cycle(choices)):
    name = 'ABCDE'[i % 5]
    mesher.pad_contour(side=side, width=i, name=name)

spiral_mesh = mesher.triangulate(opts='pAq30a200e')

spiral_mesh.plot(legend='floating', hide_labels=(0, ), linewidth=0.5)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

# Customize regions in a 2D mesh

Sometimes it's necessary to have more control over the size of the mesh in different regions. This notebook demonstrates how region markers can be used to do so.

Let's start with a simple line mesh:

```python
from nanomesh import Mesh
import numpy as np

points = np.array([[0, 0], [10, 0], [10, 10], [0, 10]])
cells = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [1, 3]])

line_mesh = Mesh(points=points, cells=cells)
line_mesh.plot_mpl()
```

Region markers are used to identify regions in the line mesh. The cell below creates a marker in each of the triangles. The markers are given an integer label, and a point that should be somewhere inside the region (i.e. bounded by line segments).

Notice the options (`opts='pAa1'`) in order:

- `p`: [planar straight line graph](https://www.cs.cmu.edu/~quake/triangle.defs.html#pslg), tells the triangulator that each line segment can be divided by adding new points
- `A`: Propagate the label to each triangle in the region
- `a1`: Constrain the global triangle size to 1 px^2

For more info on the options (or switches):

- https://rufat.be/triangle/API.html#triangle.triangulate
- https://www.cs.cmu.edu/~quake/triangle.switch.html

```python
import matplotlib.pyplot as plt
from nanomesh import RegionMarker

region_markers = [
    RegionMarker(label=1, point=(2.5, 2.5)),
    RegionMarker(label=2, point=(7.5, 7.5)),
]

line_mesh.region_markers.extend(region_markers)

fig, (ax0, ax1) = plt.subplots(ncols=2)

line_mesh.plot(ax=ax0)
ax0.set_title('line mesh with region markers')

mesh = line_mesh.triangulate(opts='pAa1')
mesh.plot(ax=ax1)
```

Region markers can be named by specifying the `name` attribute. Note that region names propagate to the `.field_data` attribute in the triangulated mesh, which means that they are available for the legend.

```python
region_markers = [
    RegionMarker(label=1, point=(2.5, 2.5), name='upper'),
    RegionMarker(label=2, point=(7.5, 7.5), name='lower'),
]

line_mesh.region_markers.clear()
line_mesh.region_markers.extend(region_markers)

mesh = line_mesh.triangulate(opts='pAa1')
mesh.plot()
```

The maximum triangle size can be limited by specifing the `constraint` attribute. Note that this requires `a` to be set, but not specified. The size of the triangles in the upper region be limited to 5 px^2, and the lower region to 0.5 px^2.

The last two examples also add `q30`, which tell `triangulate()` to make quality meshes, with a minimum triangle angles of 30 degrees.

```python
import matplotlib.pyplot as plt

region_markers = [
    RegionMarker(label=1, point=(2.5, 2.5), name='upper', constraint=0.5),
    RegionMarker(label=2, point=(7.5, 7.5), name='lower', constraint=5),
]

line_mesh.region_markers.clear()
line_mesh.region_markers.extend(region_markers)

fix, axes = plt.subplots(nrows=4, ncols=2)
axes = axes.flatten()

opts_tuple = (
    'a1',
    'a',
    'pa1',
    'pa',
    'pAa1',
    'pAa',
    'pAa1q30',
    'pAaq30',
)

for ax, opts in zip(axes, opts_tuple):
    mesh = line_mesh.triangulate(opts=opts)
    ax = mesh.plot(ax=ax)
    ax.set_title(opts)
    ax.axis('off')
```

Region markers can be accessed/modified through the `.region_markers` attribute.

```python
line_mesh.region_markers
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%load_ext autoreload
%autoreload 2
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Generate banner

This banner uses Nanomesh to generate the banner for Nanomesh.

```python
from nanomesh import Image
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
from skimage.color import rgb2gray
```

Load the source image.

```python
banner_o = io.imread(r'source_text_o.png')
plane = Image(rgb2gray(banner_o))

seg = plane.binary_digitize(threshold=0.5)
seg.show()
```

### Contour finding

```python tags=[]
from nanomesh import Mesher2D

mesher = Mesher2D(seg.image)
mesher.generate_contour(max_contour_dist=10)
mesher.plot_contour()
```

### Create the mesh

And compare with original image.

```python
mesh = mesher.triangulate(opts='pq30a2500')
seg.compare_with_mesh(mesh)
```

### Create the banner using matplotlib

```python
plt.rcParams['image.cmap'] = 'gist_rainbow'

banner_no_o = io.imread(r'source_text_no_o.png')

tri_mesh = mesh.get('triangle')

points = tri_mesh.points
triangles = tri_mesh.cells
labels = tri_mesh.labels

x, y = points.T[::-1]

fig, ax = plt.subplots(figsize=(8, 2))
fig.tight_layout(pad=0)

ax.imshow(banner_no_o)
ax.axis('off')
ax.margins(0)

colors = np.arange(len(triangles))
np.random.shuffle(colors)  # mix up the colors
mask_o = (labels == 1)
ax.tripcolor(x, y, triangles=triangles, mask=mask_o, facecolors=colors)
ax.triplot(x, y, triangles=triangles, mask=mask_o, color='black', lw=0.5)

mask_rest = (labels == 2)
ax.triplot(x, y, triangles=triangles, mask=mask_rest, lw=0.5, alpha=0.8)

plt.savefig('banner.png', bbox_inches='tight')
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: nanomesh
    language: python
    name: nanomesh
---

```python
%config InlineBackend.rc = {'figure.figsize': (10,6)}
%matplotlib inline
```

## Hello world!

<!-- #raw raw_mimetype="text/restructuredtext" -->
This example shows the workflow for generating a mesh from segmented data, and demonstrates a few of the features of Nanomesh. It uses a synthetic binary image with several rounded blob-like objects generated by :mod:`skimage`.

Image data in Nanomesh is stored as an :class:`nanomesh.Image` type. Calling :class:`nanomesh.Image` will create the appropriate subtype, :class:`nanomesh.Plane` for 2D or :class:`nanomesh.Volume` for 3D data.
<!-- #endraw -->

```python
from skimage.data import binary_blobs
from nanomesh import Image

blobs = binary_blobs(length=100, volume_fraction=0.25, seed=2102)
plane = Image(blobs)

print(plane)
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
:class:`nanomesh.Image` is essentially a container for a :mod:`numpy` array with some methods for image segmentation and visualization.
<!-- #endraw -->

```python
plane.show()
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Generating a mesh from image data is simple in Nanomesh using :meth:`nanomesh.Plane.generate_mesh`. The options `opts` are passed to the triangulation function (:func:`nanomesh.triangulate`). In this example, we use `q30` to generate a quality mesh with minimum angles of 30°, and `a50` to limit the triangle size to 50 pixels.

The returned `mesh` is a :class:`nanomesh.MeshContainer` with the generated triangles and line segments.
<!-- #endraw -->

```python
mesh = plane.generate_mesh(opts='q30a10')
mesh
```

In the next cell, we plot the triangles.

```python
mesh.plot('triangle')
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Nanomesh can also calculate cell quality metrics and show them as a colored triangle or histogram plot.

Note that the area shaded in green below highlights the *optimal* range.

Have a look at the :mod:`nanomesh.metrics` submodule or in the example: :doc:`nanopores_calculate_mesh_quality_metrics`.
<!-- #endraw -->

```python
from nanomesh import metrics

triangle_mesh = mesh.get('triangle')

metrics.histogram(triangle_mesh, metric='radius_ratio')
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Nanomesh uses [meshio](https://github.com/nschloe/meshio) to write data to most meshing formats.
<!-- #endraw -->

```python
mesh.write('mesh.vtk')
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
For more practical examples of how to use Nanomesh, check out the how-to guides and notebooks: :ref:`examples`.
<!-- #endraw -->
Class hierarchy
===============

This scheme shows the class hierarchy for Nanomesh.

.. image:: _static/flowchart_structure.svg
   :alt: Flowchart showing class structure
.. _api_metrics:

.. module:: nanomesh.metrics

Metrics
=======

The :mod:`nanomesh.metrics` module helps with calculating different types of cell metrics.

There are a few higher level functions available. While one could use
:func:`calculate_all_metrics` to calculate all available metrics,
each function is also available by itself.

:func:`histogram` and :func:`plot2d` are helpers
to visualize the metrics.

.. seealso::

    For more info, see the example on :doc:`examples/nanopores_calculate_mesh_quality_metrics`.

.. rubric:: Functions

These metrics are currently available:

.. autosummary::

   area
   aspect_frobenius
   aspect_ratio
   condition
   distortion
   max_angle
   max_min_edge_ratio
   min_angle
   radius_ratio
   relative_size_squared
   scaled_jacobian
   shape
   shape_and_size

Utility functions:

.. autosummary::

   calculate_all_metrics
   histogram
   plot2d

Reference
---------

.. autofunction:: calculate_all_metrics
.. autofunction:: histogram
.. autofunction:: plot2d

.. autofunction:: area
.. autofunction:: aspect_frobenius
.. autofunction:: aspect_ratio
.. autofunction:: condition
.. autofunction:: distortion
.. autofunction:: max_angle
.. autofunction:: max_min_edge_ratio
.. autofunction:: min_angle
.. autofunction:: radius_ratio
.. autofunction:: relative_size_squared
.. autofunction:: scaled_jacobian
.. autofunction:: shape
.. autofunction:: shape_and_size
.. _api_mesh_data:

.. currentmodule:: nanomesh

Working with mesh data
======================

A large part of Nanomesh deals with with generating and manipulating
2D and 3D meshes. To store the data, Nanomesh uses :class:`MeshContainer`,
which is a somewhat low-level, generic container for mesh data. It can
store different types of cells and associated data. It is used
to read/write data via `meshio <https://github.com/nschloe/meshio>`_.

To deal with mesh data more directly, use a :class:`Mesh`. Instantiating
:class:`Mesh` will create the appropriate subtype, :class:`LineMesh`,
:class:`TriangleMesh` or :class:`TetraMesh`. These contain
dedicated methods for working with a specific type of mesh and plotting them.

In addition, each can be extracted from :class:`MeshContainer`.


.. rubric:: Classes

.. autosummary::

   nanomesh.MeshContainer
   nanomesh.Mesh
   nanomesh.LineMesh
   nanomesh.TriangleMesh
   nanomesh.TetraMesh


Reference
---------

.. autoclass:: MeshContainer
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Mesh
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: LineMesh
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: TriangleMesh
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: TetraMesh
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
.. _api:

.. currentmodule:: nanomesh

Python Interface
================

This part of the documentation covers the public interface of Nanomesh.

The side bar contains a listing of classes and functions by topic.

.. toctree::
   :maxdepth: 1
   :hidden:

   api.image_data
   api.mesh_data
   api.meshing
   api.metrics
   api.helpers
   api.data

Most of Nanomesh' functionality can be accessed through the classes
and functions listed below. See the :ref:`examples` for how to use them.

.. rubric:: Data types

.. autosummary::

   nanomesh.Plane
   nanomesh.Volume
   nanomesh.Image
   nanomesh.MeshContainer
   nanomesh.Mesh
   nanomesh.LineMesh
   nanomesh.TriangleMesh
   nanomesh.TetraMesh
   nanomesh.RegionMarker

.. rubric:: Meshing classes

.. autosummary::

   nanomesh.Mesher
   nanomesh.Mesher2D
   nanomesh.Mesher3D

.. rubric:: Functions

.. autosummary::

   nanomesh.triangulate
   nanomesh.tetrahedralize
   nanomesh.volume2mesh
   nanomesh.plane2mesh

.. rubric:: Modules

.. autosummary::

   nanomesh.metrics
   nanomesh.data


Reference
=========

The complete API reference is listed below.

.. toctree::
   :maxdepth: 2
   :glob:

   nanomesh/*
.. _api_data:

.. module:: nanomesh.data

Data
=====

The :mod:`nanomesh.data` module helps provides some standard data sets to work with.

.. seealso::

    For additional data sets, have a look at :mod:`skimage.data`:

   - :func:`skimage.data.cells3d`
   - :func:`skimage.data.horse`
   - :func:`skimage.data.binary_blobs`
   - :func:`skimage.data.coins`

.. rubric:: image data

These image data are currently available in Nanomesh.

.. autosummary::

   binary_blobs2d
   binary_blobs3d
   mesh
   mesh3d
   nanopores
   nanopores3d
   nanopores_gradient

Reference
---------

.. autofunction:: binary_blobs2d
.. autofunction:: binary_blobs3d
.. autofunction:: mesh
.. autofunction:: mesh3d
.. autofunction:: nanopores
.. autofunction:: nanopores3d
.. autofunction:: nanopores_gradient
Nanomesh documentation
======================

|Documentation Status| |tests| |PyPI - Python Version| |PyPI| |DOI|

.. figure:: _static/banner.png
   :alt: Nanomesh banner

Welcome to the nanomesh documentation!

Nanomesh is a Python workflow tool for generating meshes from 2D and 3D image data. It has an easy-to-use API that can help process and segment image data, generate quality meshes (triangle / tetrahedra), and write the data to many mesh formats. Nanomesh also contains tools to inspect the meshes, visualize them, and generate cell quality metrics.

- Easy-to-use Python API
- Segment and mesh 2D or 3D image data
- Mesh visualization
- Calculate and plot cell metrics
- Export to many mesh formats


.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   install
   development
   examples/other_hello_world!
   structure


.. toctree::
   :maxdepth: 1
   :caption: How-to's

   examples/index


.. toctree::
   :maxdepth: 1
   :caption: API

   api.rst


.. toctree::
   :caption: Links

   👨‍💻 Source code <https://github.com/HPGEM/nanomesh>
   💡 Issues <https://github.com/HPGEM/nanomesh/issues>
   📢 Releases <https://github.com/hpgem/nanomesh/releases>
   🐍 PyPI <https://pypi.org/project/nanomesh>
   📚 documentation <https://nanomesh.readthedocs.io>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |Documentation Status| image:: https://readthedocs.org/projects/nanomesh/badge/?version=latest
   :target: https://nanomesh.readthedocs.io/en/latest/?badge=latest
.. |tests| image:: https://github.com/hpgem/nanomesh/actions/workflows/test.yaml/badge.svg
   :target: https://github.com/hpgem/nanomesh/actions/workflows/test.yaml
.. |PyPI - Python Version| image:: https://img.shields.io/pypi/pyversions/nanomesh
   :target: https://pypi.org/project/nanomesh/
.. |PyPI| image:: https://img.shields.io/pypi/v/nanomesh.svg?style=flat
   :target: https://pypi.org/project/nanomesh/
.. |DOI| image:: https://zenodo.org/badge/311460276.svg
   :target: https://zenodo.org/badge/latestdoi/311460276
Installation
============

One of the goals for Nanomesh is that it is easy to install.
This means that all dependencies are available from `PyPi <https://pypi.org>`_.

If you use conda, it is advised to create a new environment:

::

   conda create -n nanomesh python=3.8
   conda activate nanomesh

Install nanomesh:

::

   pip install nanomesh

Note, `to enable the IPython
widgets <https://ipywidgets.readthedocs.io/en/latest/user_install.html#installation>`__:

::

   jupyter nbextension enable --py widgetsnbextension

Note, `if you are using Jupyter
lab <https://github.com/InsightSoftwareConsortium/itkwidgets#installation>`__:

::

   jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib jupyterlab-datawidgets itkwidgets

Tetgen
------

For tetrahedral meshing, Nanomesh requires `tetgen <https://wias-berlin.de/software/tetgen/>`__ to be
installed. Binaries are available `here <https://github.com/hpgem/tetgen/releases>`__, and compilation instructions `here <https://github.com/hpgem/tetgen/releases>`_.

Make sure `tetgen` is available on a directory on your system path. To verify tetgen is available, make sure that the following commands return a path:

Linux/MacOS

::

   which tetgen

Windows

::

   gcm tetgen.exe
.. _api_meshing:

.. currentmodule:: nanomesh

Making a mesh out of image data
===============================

:func:`plane2mesh` and :func:`volume2mesh` are high-level functions for
generating a mesh. For finer control, use the classes :class:`Mesher2D` for
image data and :class:`Mesher3D` for volume data.

Both classes derive from :class:`Mesher`. Initiating a :class:`Mesher` instance
will create the appropriate meshing class.

.. rubric:: Classes

.. autosummary::

   nanomesh.Mesher
   nanomesh.Mesher2D
   nanomesh.Mesher3D

.. rubric:: Functions

.. autosummary::

   nanomesh.plane2mesh
   nanomesh.volume2mesh

Reference
---------

.. autofunction:: plane2mesh
.. autofunction:: volume2mesh

.. autoclass:: Mesher
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Mesher2D
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Mesher3D
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
.. _api_image_data:

.. currentmodule:: nanomesh

Working with image data
=======================

Nanomesh has two classes for representing image data, :class:`Plane` for working
with 2D pixel data and :class:`Volume` for working with 3D voxel data.
These can be used to load, crop, transform, filter, and segment image data.

Both classes derive from :class:`Image`. Instantiating :class:`Image` will
create the appropriate subclass.


.. rubric:: Classes

.. autosummary::

   nanomesh.Image
   nanomesh.Plane
   nanomesh.Volume

Reference
---------

.. autoclass:: Image
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Plane
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Volume
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
.. _api_helpers:

.. currentmodule:: nanomesh

Helper classes and functions
============================

These are helper functions for working with meshes. A :class:`RegionMarker`
stores metadata about a region in a mesh. :func:`tetrahedralize` and
:func:`triangulate` are interfaces to generate meshes.

.. rubric:: classes

.. autosummary::

   nanomesh.RegionMarker

.. rubric:: functions

.. autosummary::

   nanomesh.tetrahedralize
   nanomesh.triangulate

Reference
---------

.. autoclass:: RegionMarker

.. autofunction:: tetrahedralize
.. autofunction:: triangulate
Development Notes
=================

Development Installation
------------------------

Install Nanomesh using the development dependencies:

::

    conda create -n nanomesh-dev python=3.8
    conda activate nanomesh-dev

    pip install -e .[develop] -c constraints.txt

Running the tests:

::

    pytest

Linting/checks:

::

    pre-commit

Building the docs:

::

   make html --directory docs


Testing notebooks
-----------------

1. Install nanomesh kernel

   ::

       python -m ipykernel install --user --name nanomesh

2. Test notebooks

   ::

       cd notebooks
       pip install -r requirements.txt

   On Windows:

   ::

       ./test_notebooks.PS1

   On Linux/Mac:

   ::

       bash ./test_notebooks.sh


Making a release
----------------

1. Bump the version (major/minor/patch as needed)

    ::

        bumpversion minor

2. Make a new release. The github action to publish to pypi is triggered when a release is published.


Updating constraints.txt
------------------------

1. On Windows:
    - In a new environment

::

    pip freeze --exclude nanomesh > constraints.txt

2. On Linux:
    - In a new environment
    - Using the produced ``constraints.txt`` file

::

    pip install -e .[develop] -c constraints.txt
    pip freeze --exclude nanomesh >> constraints.txt
    sort --ignore-case constraints.txt | uniq > constraints_tmp.txt
    mv constraints_tmp.txt constraints.txt


Updating pre-commit
-------------------

::

    pre-commit autoupdate


Fixes for errors
----------------

If you get an error with pytest, like:

::

     from win32com.shell import shellcon, shell
    E   ImportError: DLL load failed while importing shell: The specified procedure could not be found.
    ImportError: DLL load failed while importing shell: The specified procedure could not be found.

Try:

::

    conda install -c conda-forge pywin32
nanomesh.mesh
=============

.. automodule:: nanomesh.mesh
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: LineMesh, TriangleMesh, TetraMesh
nanomesh.utils
==============

.. automodule:: nanomesh.utils
   :members:
   :undoc-members:
   :show-inheritance:
nanomesh.image2mesh.mesher3d
============================

.. automodule:: nanomesh.image2mesh._mesher3d
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: Mesher3D
nanomesh.io
===========

.. automodule:: nanomesh.io
   :members:
   :undoc-members:
   :show-inheritance:
nanomesh.region\_markers
========================

.. automodule:: nanomesh.region_markers
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
nanomesh.mesh\_container
========================

.. automodule:: nanomesh.mesh_container
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: MeshContainer
nanomesh.metrics
================

.. automodule:: nanomesh.metrics
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
nanomesh.image
==============

.. automodule:: nanomesh.image
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: Plane, Volume
nanomesh.image2mesh
===================

.. toctree::

   image2mesh._mesher2d
   image2mesh._mesher3d

.. automodule:: nanomesh.image2mesh
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: Mesher2D, Mesher3D
nanomesh.image2mesh.mesher2d
============================

.. automodule:: nanomesh.image2mesh._mesher2d
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: Mesher2D
nanomesh.plotting
=================

.. automodule:: nanomesh.plotting
   :members:
   :undoc-members:
   :show-inheritance:
.. _examples:

Overview
========

This page contains a few examples in jupyter notebooks.

The source code for the notebooks is available from
`Github <https://github.com/hpgem/nanomesh/tree/master/notebooks/>`_.

Examples
========

.. toctree::
   :maxdepth: 1
   :glob:

   examples*

Nanopores
=========

.. toctree::
   :maxdepth: 1
   :glob:

   nanopores*

Other
=====

.. toctree::
   :maxdepth: 1
   :glob:

   other*
