<img alt="IFermi logo" src="https://raw.githubusercontent.com/fermisurfaces/IFermi/main/docs/src/_static/logo2-01.png" height="150px">

--------

[📖 **Official Documentation** 📖](https://fermisurfaces.github.io/IFermi)

[🙋 **Support Forum** 🙋](https://matsci.org/c/ifermi/)

[📝 **JOSS Paper** 📝](https://doi.org/10.21105/joss.03089)

IFermi is a Python (3.6+) library and set of command-line tools for the generation,
analysis, and visualisation of Fermi surfaces and Fermi slices. The goal of the library
is to provide fully featured FermiSurface and FermiSlice objects that allow for easy
manipulation and analysis. The main features include:

- Interpolation of electronic band structures onto dense k-point meshes.
- Extraction of Fermi surfaces and Fermi slices from electronic band structures.
- Projection of arbitrary properties onto Fermi surfaces and Fermi slices.
- Tools to calculate Fermi surface dimensionality, orientation, and averaged projections,
  including Fermi velocities.
- Interactive visualisation of Fermi surfaces and slices, with support for
  [mayavi](https://docs.enthought.com/mayavi/mayavi/), [plotly](https://plot.ly/) and
  [matplotlib](https://matplotlib.org).
- Generation and visualisation of spin-texture.

IFermi's command-line tools only work with VASP calculations but support for additional
DFT packages will be added in the future.

![Example Fermi surfaces](https://raw.githubusercontent.com/fermisurfaces/IFermi/main/docs/src/_static/fermi-surface-example.png)

## Quick start

The [online documentation](https://fermisurfaces.github.io/IFermi/cli.html) provides a full
description of the available command-line options.

### Analysis

Fermi surface properties, including dimensionality and orientation can be extracted
from a vasprun.xml file using:

```bash
ifermi info --property velocity
```

```
Fermi Surface Summary
=====================

  # surfaces: 5
  Area: 32.75 Å⁻²
  Avg velocity: 9.131e+05 m/s

Isosurfaces
~~~~~~~~~~~

      Band    Area [Å⁻²]    Velocity avg [m/s]   Dimensionality    Orientation
    ------  ------------  --------------------  ----------------  -------------
         6         1.944             7.178e+05         2D           (0, 0, 1)
         7         4.370             9.092e+05      quasi-2D        (0, 0, 1)
         7         2.961             5.880e+05         2D           (0, 0, 1)
         8         3.549             1.105e+06      quasi-2D        (0, 0, 1)
         8         3.549             1.105e+06      quasi-2D        (0, 0, 1)
```

### Visualisation

Three-dimensional Fermi surfaces can be visualized from a `vasprun.xml` file using:

```bash
ifermi plot
```

The two-dimensional slice of a Fermi surface along the plane specified by the miller
indices (j k l) and distance d can be plotted from a `vasprun.xml` file using:

```bash
ifermi plot --slice j k l d
```

### Python library

The `ifermi` command line tools are build on the IFermi Python library. Here is an
example of how to load DFT calculation outputs, interpolate the energies onto a dense mesh,
generate a Fermi surface, calculate Fermi surface properties, and visualise the surface.
A more complete summary of the API is given in the [API introduction page](https://fermisurfaces.github.io/IFermi/introduction_to_ifermi.html)
and in the [API Reference page](https://fermisurfaces.github.io/IFermi/ifermi.html) in the documentation.

```python
from pymatgen.io.vasp.outputs import Vasprun
from ifermi.surface import FermiSurface
from ifermi.interpolate import FourierInterpolator
from ifermi.plot import FermiSlicePlotter, FermiSurfacePlotter, save_plot, show_plot
from ifermi.kpoints import kpoints_from_bandstructure

# load VASP calculation outputs
vr = Vasprun("vasprun.xml")
bs = vr.get_band_structure()

# interpolate the energies onto a dense k-point mesh
interpolator = FourierInterpolator(bs)
dense_bs, velocities = interpolator.interpolate_bands(return_velocities=True)

# generate the Fermi surface and calculate the dimensionality
fs = FermiSurface.from_band_structure(
  dense_bs, mu=0.0, wigner_seitz=True, calculate_dimensionality=True
)

# generate the Fermi surface and calculate the group velocity at the
# center of each triangular face
dense_kpoints = kpoints_from_bandstructure(dense_bs)
fs = FermiSurface.from_band_structure(
  dense_bs, mu=0.0, wigner_seitz=True, calculate_dimensionality=True,
  property_data=velocities, property_kpoints=dense_kpoints
)

# number of isosurfaces in the Fermi surface
fs.n_surfaces

# number of isosurfaces for each Spin channel
fs.n_surfaces_per_spin

# the total area of the Fermi surface
fs.area

# the area of each isosurface
fs.area_surfaces

# loop over all isosurfaces and check their properties
# the isosurfaces are given as a list for each spin channel
for spin, isosurfaces in fs.isosurfaces.items():
    for isosurface in isosurfaces:
        # the dimensionality (does the surface cross periodic boundaries)
        isosurface.dimensionality

        # what is the orientation
        isosurface.orientation

        # does the surface have face properties
        isosurface.has_properties

        # calculate the norms of the properties
        isosurface.properties_norms

        # calculate scalar projection of properties on to [0 0 1] vector
        isosurface.scalar_projection((0, 0, 1))

        # uniformly sample the surface faces to a consistent density
        isosurface.sample_uniform(0.1)

# plot the Fermi surface
fs_plotter = FermiSurfacePlotter(fs)
plot = fs_plotter.get_plot()

# generate Fermi slice along the (0 0 1) plane going through the Γ-point.
fermi_slice = fs.get_fermi_slice((0, 0, 1))

# number of isolines in the slice
fermi_slice.n_lines

# do the lines have segment properties
fermi_slice.has_properties

# plot slice
slice_plotter = FermiSlicePlotter(fermi_slice)
plot = slice_plotter.get_plot()

save_plot(plot, "fermi-slice.png")  # saves the plot to a file
show_plot(plot)  # displays an interactive plot
```

## Citing IFermi

If you find IFermi useful, please encourage its development by citing the following
[paper](https://doi.org/10.21105/joss.03089) in your research output:

```
Ganose, A. M., Searle, A., Jain, A., Griffin, S. M., IFermi: A python library for Fermi
surface generation and analysis. Journal of Open Source Software, 2021, 6 (59), 3089
```


## Installation

The recommended way to install IFermi is in a conda environment.

```bash
conda create --name ifermi pip cmake numpy
conda activate ifermi
conda install -c conda-forge pymatgen boltztrap2 pyfftw
pip install ifermi
````

IFermi is currently compatible with Python 3.6+ and relies on a number of
open-source python packages, specifically:

- [pymatgen](http://pymatgen.org) for parsing DFT calculation outputs.
- [BoltzTrap2](https://gitlab.com/sousaw/BoltzTraP2) for band structure interpolation.
- [trimesh](https://trimsh.org/) for manipulating isosurfaces.
- [matplotlib](https://matplotlib.org), [mayavi](https://docs.enthought.com/mayavi/mayavi/), and [plotly](https://plot.ly/) for three-dimensional plotting.

### Running tests

The integration tests can be run to ensure IFermi has been installed correctly. First
download the IFermi source and install the test requirements.

```
git clone https://github.com/fermisurfaces/IFermi.git
cd IFermi
pip install .[tests]
```

The tests can be run in the IFermi folder using:

```bash
pytest
```

## Need Help?

Ask questions about the IFermi Python API and command-line tools on the [IFermi
support forum](https://matsci.org/c/ifermi).
If you've found an issue with IFermi, please submit a bug report
[here](https://github.com/fermisurfaces/IFermi/issues).

## What’s new?

Track changes to IFermi through the
[changelog](https://fermisurfaces.github.io/IFermi/changelog.html).

## Contributing

We greatly appreciate any contributions in the form of a pull request.
Additional information on contributing to IFermi can be found [here](https://fermisurfaces.github.io/IFermi/contributing.html).
We maintain a list of all contributors [here](https://fermisurfaces.github.io/IFermi/contributors.html).

## License

IFermi is made available under the MIT License (see LICENSE file).

## Acknowledgements

Developed by Amy Searle and Alex Ganose.
Sinéad Griffin designed and led the project.
## Summary

Include a summary of major changes in bullet points:

* Feature 1
* Feature 2
* Fix 1
* Fix 2

## Additional dependencies introduced (if any)

* List all new dependencies needed and justify why. While adding dependencies that bring
significantly useful functionality is perfectly fine, adding ones that
add trivial functionality, e.g., to use one single easily implementable
function, is frowned upon. Provide a justification why that dependency is needed.
Especially frowned upon are circular dependencies, e.g., depending on derivative
modules of pymatgen such as custodian or Fireworks.

## TODO (if any)

If this is a work-in-progress, write something about what else needs
to be done

* Feature 1 supports A, but not B.

## Checklist

Work-in-progress pull requests are encouraged, but please put [WIP]
in the pull request title.

Before a pull request can be merged, the following items must be checked:

- [ ] Code is in the [standard Python style](https://www.python.org/dev/peps/pep-0008/). The easiest way to handle this
      is to run the following in the **correct sequence** on your local machine. Start with running
      [black](https://black.readthedocs.io/en/stable/index.html) on your new code. This will automatically reformat
      your code to PEP8 conventions and removes most issues. Then run
      [pycodestyle](https://pycodestyle.readthedocs.io/en/latest/), followed by
      [flake8](http://flake8.pycqa.org/en/latest/).
- [ ] Docstrings have been added in the [Google docstring format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
      Run [pydocstyle](http://www.pydocstyle.org/en/2.1.1/index.html) on your code.
- [ ] Type annotations are **highly** encouraged. Run [mypy](http://mypy-lang.org/) to type check your code.
- [ ] Tests have been added for any new functionality or bug fixes.
- [ ] All linting and tests pass.

Note that the CI system will run all the above checks. But it will be much more efficient if you already fix most
errors prior to submitting the PR. It is highly recommended that you use the pre-commit hook provided in the pymatgen
repository. Simply `cp pre-commit .git/hooks` and a check will be run prior to allowing commits.
---
name: Bug Report
about: Create a bug report to help us improve IFermi
title: "BUG:"
---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

Provide any example files that are needed to reproduce the error,
especially if the bug pertains to parsing of a file.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.
# Releasing a new IFermi version

Version releases on Pypi and GitHub are handled automatically through GitHub
actions. The steps to push a new release are:
1. Update the version in `ifermi/__init__.py`
2. Update the changelog in `CHANGELOG.rst` with the new version and
   release notes.
3. Create a tagged Git commit with the above changes. The tag is added using:
```bash
git tag v0.2.5
```
4. Push the commit and tags to GitHub using:
```bash
git push origin --tags
```
Contributing to IFermi
======================

We love your input! We want to make contributing to as easy and
transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing or implementing new features
- Becoming a maintainer

Reporting bugs, getting help, and discussion
--------------------------------------------

Please submit questions to the dedicated IFermi `help forum <https://matsci.org/c/ifermi>`__.
Bug reports should be directed to the `GitHub issues page <https://github.com/fermisurfaces/IFermi/issues>`__.

If you are making a bug report, incorporate as many elements of the
following as possible to ensure a timely response and avoid the
need for followups:

- A quick summary and/or background.
- Steps to reproduce - be specific! **Provide sample code.**
- What you expected would happen, compared to what actually happens.
- The full stack trace of any errors you encounter.
- Notes (possibly including why you think this might be happening,
  or steps you tried that didn't work).

We love thorough bug reports as this means the development team can
make quick and meaningful fixes. When we confirm your bug report,
we'll move it to the GitHub issues where its progress can be
further tracked.

Contributing code modifications or additions through Github
-----------------------------------------------------------

We use github to host code, to track issues and feature requests,
as well as accept pull requests. We maintain a list of all
contributors `here
<https://fermisurfaces.github.io/IFermi/contributors.html>`__.

Pull requests are the best way to propose changes to the codebase.
Follow the `Github flow
<https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow>`__
for more information on this procedure.

The basic procedure for making a PR is:

- Fork the repo and create your branch from main.
- Commit your improvements to your branch and push to your Github fork (repo).
- When you're finished, go to your fork and make a Pull Request. It will
  automatically update if you need to make further changes.

How to Make a **Great** Pull Request
------------------------------------

We have a few tips for writing good PRs that are accepted into the main repo:

- Use the Google Code style for all of your code. Find an example `here
  <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`__.
- Your code should have (4) spaces instead of tabs.
- If needed, update the documentation.
- **Write tests** for new features! Good tests are 100%, absolutely necessary
  for good code. We use the python ``unittest`` framework -- see some of the
  other tests in this repo for examples, or review the `Hitchhiker's guide
  to python <https://docs.python-guide.org/writing/tests/>`__ for some good
  resources on writing good tests.
- Understand your contributions will fall under the same license as this repo.

When you submit your PR, our CI service will automatically run your tests.
We welcome good discussion on the best ways to write your code, and the comments
on your PR are an excellent area for discussion.
Contributors
============

IFermi was developed by **Amy Searle** |ajsearle| and
**Alex Ganose** |utf| |0000-0002-4486-3321|.
The project was designed and led by **Sinéad Griffin** |0000-0002-9943-4866| at Lawrence
Berkeley National Laboratory.


.. |ajsearle| image:: https://cdnjs.cloudflare.com/ajax/libs/octicons/8.5.0/svg/mark-github.svg
   :target: https://github.com/ajsearle97
   :width: 16
   :height: 16
   :alt: GitHub profile for utf
.. |0000-0002-9943-4866| image:: _static/orcid.svg
   :target: https://orcid.org/0000-0002-9943-4866
   :width: 16
   :height: 16
   :alt: ORCID profile for 0000-0002-9943-4866
.. |utf| image:: https://cdnjs.cloudflare.com/ajax/libs/octicons/8.5.0/svg/mark-github.svg
   :target: https://github.com/utf
   :width: 16
   :height: 16
   :alt: GitHub profile for utf
.. |0000-0002-4486-3321| image:: _static/orcid.svg
   :target: https://orcid.org/0000-0002-4486-3321
   :width: 16
   :height: 16
   :alt: ORCID profile for 0000-0002-4486-3321

Additional contributions have been provided by:

| **Matthew Horton** |mkhorton| |0000-0001-7777-8871|
| Project Scientist, Persson Group
| Lawrence Berkeley National Laboratory

.. |mkhorton| image:: https://cdnjs.cloudflare.com/ajax/libs/octicons/8.5.0/svg/mark-github.svg
   :target: https://github.com/mkhorton
   :width: 16
   :height: 16
   :alt: GitHub commits from mkhorton

.. |0000-0001-7777-8871| image:: _static/orcid.svg
   :target: https://orcid.org/0000-0001-7777-8871
   :width: 16
   :height: 16
   :alt: ORCID profile for 0000-0001-7777-8871

| **Mariano Forti** |mdforti| |0000-0001-7366-3372|
| Postdoctoral Researcher
| Interdisciplinary Centre for Advanced Materials Simulation

.. |mdforti| image:: https://cdnjs.cloudflare.com/ajax/libs/octicons/8.5.0/svg/mark-github.svg
   :target: https://github.com/mdforti
   :width: 16
   :height: 16
   :alt: GitHub commits from mdforti
.. |0000-0001-7366-3372| image:: _static/orcid.svg
   :target: https://orcid.org/0000-0001-7366-3372
   :width: 16
   :height: 16
   :alt: ORCID profile for 0000-0001-7366-3372
ifermi package
==============

ifermi.analysis module
-----------------------

.. automodule:: ifermi.analysis
   :members:
   :undoc-members:
   :show-inheritance:

ifermi.brillouin\_zone module
-----------------------------

.. automodule:: ifermi.brillouin_zone
   :members:
   :show-inheritance:

ifermi.interpolate module
--------------------------

.. automodule:: ifermi.interpolate
   :members:
   :undoc-members:
   :show-inheritance:

ifermi.kpoints module
---------------------

.. automodule:: ifermi.kpoints
   :members:
   :undoc-members:
   :show-inheritance:

ifermi.plot module
---------------------

.. automodule:: ifermi.plot
   :members:
   :undoc-members:
   :show-inheritance:

ifermi.slice module
-------------------

.. automodule:: ifermi.slice
   :members:
   :show-inheritance:

ifermi.surface module
---------------------

.. automodule:: ifermi.surface
   :members:
   :show-inheritance:
``ifermi`` program
==================

IFermi includes command-line tools for generating, analysing, and plotting Fermi
surfaces. The tools can be accessed using the ``ifermi`` command-line program.
All of the options provided in the command-line are also accessible using the
`Python API <plotting_using_python.html>`_.

IFermi works in 4 stages:

1. It loads a band structure from DFT calculation outputs.
2. It interpolates the band structure onto a dense k-point mesh using Fourier
   interpolation as implemented in `BoltzTraP2 <https://gitlab.com/sousaw/BoltzTraP2>`_.
3. It extracts the Fermi surface at a given energy level.
4. It extracts information about the Fermi surface or plots it using several
   plotting backends.

.. NOTE::

    Currently, IFermi's command-line tools only work with VASP calculation outputs.
    Support for additional DFT packages will be added in a future release.

IFermi is controlled on the command-line using the ``ifermi`` command. The available
options can be listed using:

.. code-block:: bash

    ifermi -h

Information on the Fermi surface area, dimensionality,
and orientation can be extracted using the ``info`` subcommand.
The only input required is a vasprun.xml. For example:

.. code-block:: bash

    ifermi info

An example output for MgB\ :sub:`2` is shown below:

.. code-block:: md

    Fermi Surface Summary
    =====================

      # surfaces: 5
      Area: 32.745 Å⁻²

    Isosurfaces
    ~~~~~~~~~~~

        Band    Area [Å⁻²]   Dimensionality    Orientation
      ------  ------------  ----------------  -------------
           6         1.944         2D           (0, 0, 1)
           7         4.370         1D           (0, 0, 1)
           7         2.961         2D           (0, 0, 1)
           8         3.549         1D           (0, 0, 1)
           8         3.549         1D           (0, 0, 1)


If properties are included in the Fermi surface (see :ref:`property-gen`), the averaged
property values will also be calculated. This allows for calculation of the Fermi
velocity. For example:

.. code-block:: bash

    ifermi info --property velocity

.. code-block:: md

    Fermi Surface Summary
    =====================

      # surfaces: 5
      Area: 32.75 Å⁻²
      Avg velocity: 9.131e+05 m/s

    Isosurfaces
    ~~~~~~~~~~~

        Band    Area [Å⁻²]    Velocity avg [m/s]   Dimensionality    Orientation
      ------  ------------  --------------------  ----------------  -------------
           6         1.944             7.178e+05         2D           (0, 0, 1)
           7         4.370             9.092e+05      quasi-2D        (0, 0, 1)
           7         2.961             5.880e+05         2D           (0, 0, 1)
           8         3.549             1.105e+06      quasi-2D        (0, 0, 1)
           8         3.549             1.105e+06      quasi-2D        (0, 0, 1)


Fermi surfaces and slices can be visualised using the ``plot`` subcommand. Again, the
only input required is a vasprun.xml file. For example:

.. code-block:: bash

    ifermi plot

The two subcommands ``info`` and ``plot`` share many of the same options
which we describe below.

Generation options
------------------

There are several options affect the generation of Fermi surfaces from *ab initio*
calculation outputs. These options are available for both the ``info`` and ``plot``
subcommands.

Input file
~~~~~~~~~~

IFermi will look for a vasprun.xml or vasprun.xml.gz file in the current directory.
To specify a particular vasprun file the ``--filename`` (or ``-f`` for short) option
can be used:

.. code-block:: bash

    ifermi plot --filename my_vasprun.xml

Interpolation factor
~~~~~~~~~~~~~~~~~~~~

The band structure extracted from the vasprun must be processed before the Fermi
surface can be generated. The two key issues are:

1. It may only contain the irreducible portion of the Brillouin zone (if symmetry was
   used in the calculation) and therefore may not contain enough information to plot
   the Fermi surface across the full reciprocal lattice.
2. It may have been calculated on a relatively coarse k-point mesh and will therefore
   produce a rather jagged Fermi surface.

Both issues can be solved by interpolating the band structure onto a denser k-point
mesh. This is achieved by using `BoltzTraP2 <https://gitlab.com/sousaw/BoltzTraP2>`_
to Fourier interpolate the eigenvalues onto a denser mesh that covers the full
Brillouin zone.

The degree of interpolation is controlled by the ``--interpolation-factor`` (``-i``)
argument. A value of 8 (the default value), roughly indicates that the interpolated band
structure will contain 8x as many k-points. Increasing the interpolation factor will
result in smoother Fermi surfaces. For example:

.. code-block:: bash

    ifermi plot --interpolation-factor 10

.. WARNING::

    As the interpolation increases, the generation of the Fermi surface, analysis and
    plotting will take a longer time and can result in large file sizes.

Fermi surface energy
~~~~~~~~~~~~~~~~~~~~

The energy level offset at which the Fermi surface is calculated is controlled by the
``--mu`` option. The energy level is given relative to the Fermi level of the VASP
calculation and is given in eV. By default, the Fermi surface is calculated at
``mu = 0``, i.e., at the Fermi level.

For gapped materials, ``mu`` must be selected so that it falls within the
conduction or valence bands otherwise no Fermi surface will be obtained. For
example. The following command will generate the Fermi surface at 1 eV above the Fermi
level:

.. code-block:: bash

    ifermi plot --mu 1


.. _property-gen:

Property projections
~~~~~~~~~~~~~~~~~~~~

Additional properties, such as the group velocity and orbital magnetisation (spin
texture), can be projected onto the Fermi surface using the ``--property`` option. The
group velocities are calculated during Fourier interpolation (units of m/s) and can be
included in the Fermi surface using:

.. code-block:: bash

    ifermi plot --property velocity


For calculations performed using spin–orbit coupling or non-collinear magnetism, the
spin magnetisation can be projected onto the Fermi surface using:

.. code-block:: bash

    ifermi plot --property spin

.. WARNING::

    Projecting the spin magnetisation requires the k-point mesh to cover the entire
    Brillouin zone. I.e., the DFT calculation must have been performed without symmetry
    (``ISYM = - 1`` in VASP).

It is possible to calculate the scalar projection of the the Fermi surface properties
onto a cartesian axis using the ``--projection-axis`` option.. For example, to use the
scalar projection of the spin magnetisation onto the [0 0 1] cartesian direction:

.. code-block:: bash

    ifermi plot --property spin --projection-axis 0 0 1

Reciprocal space
~~~~~~~~~~~~~~~~

By default, the Wigner–Seitz cell is used to contain to the Fermi surface. The
parallelepiped reciprocal lattice cell can be used instead by selecting the
``--no-wigner`` option. For example:

.. code-block:: bash

    ifermi plot --no-wigner


Visualisation options
---------------------

In addition to the options for generating Fermi surfaces, there are several options
that control the visualisation parameters. These options are only available for the
``plot`` subcommand.

Plotting backend
~~~~~~~~~~~~~~~~

IFermi supports multiple plotting backends. The default is to the
`plotly <http://plotly.com>`_ package but `matplotlib <http://matplotlib.org>`_ and
`mayavi <https://docs.enthought.com/mayavi/mayavi/>`_ are also supported.

.. NOTE::

    The mayavi dependencies are not installed by default. To use this backend, follow
    the installation instructions `here <https://docs.enthought.com/mayavi/mayavi/installation.html>`_
    and then install IFermi using ``pip install ifermi[mayavi]``.

Different plotting packages can be specified using the ``--type`` option (``-t``). For
example, to use matplotlib:

.. code-block:: bash

    ifermi plot --type matplotlib

Output files
~~~~~~~~~~~~

By default, IFermi generates interactive plots. To generate static images, an output
file can be specified using the ``--output`` (``-o``) option. For example:

.. code-block:: bash

    ifermi plot --output fermi-surface.jpg

.. NOTE::

    Saving graphical output files with the plotly backend requires plotly-orca to be
    installed.

Running the above command in the ``examples/MgB2`` directory produces the plot:

.. image:: _static/fs-1.jpg
    :height: 250px
    :align: center

Interactive plots can be saved to a html file using the plotly backend by specifying
a html filename. This will prevent the plot from being opened automatically.

.. code-block:: bash

    ifermi plot --output fermi-surface.html

Selecting spin channels
~~~~~~~~~~~~~~~~~~~~~~~

In the plot above, the spins are degenerate (the Hamiltonian does not differentiate
between the up and down spins). This is why the surface looks dappled, IFermi
is plotting two redundant surfaces. To stop it from doing this, we can specify that
only one spin component should be plotted using the ``--spin`` option. The default
is to plot both spins but a single spin channel can be selected through the names
"up" and "down". For example:

.. code-block:: bash

    ifermi plot --spin up

.. image:: _static/fs-spin-up.jpg
    :height: 250px
    :align: center


Changing the viewpoint
~~~~~~~~~~~~~~~~~~~~~~

The viewpoint (camera angle) can be changed using the ``--azimuth`` (``-a``) and
``--elevation`` (``-e``) options. This will change both the initial viewpoint
for interactive plots, and the final viewpoint for static plots. To summarize:

- The azimuth is the angle subtended by the viewpoint position vector on a sphere
  projected onto the x-y plane in degrees. The default is 45°.
- The elevation (or zenith) is the angle subtended by the viewpoint position vector and
  the z-axis. The default is 35°.

The viewpoint always directed to the center of the the Fermi surface (position [0 0 0]).
As an example, the viewpoint could be changed using:

.. code-block:: bash

    ifermi plot --azimuth 120 --elevation 5

.. image:: _static/fs-viewpoint.jpg
    :height: 250px
    :align: center

.. _prop-style:

Styling face properties
~~~~~~~~~~~~~~~~~~~~~~~

As described in the :ref:`property-gen` section, Fermi surfaces (and Fermi slices)
can include a property projected onto the isosurface faces. By default, if properties
are included in the Fermi surface they will be indicated by coloring of the isosurface.
If the face property is a vector, the norm of the vector will be used as the
color intensity. The colormap of the surface can be changed using the
``--property-colormap`` option. All `matplotlib colormaps <https://matplotlib.org/stable/gallery/color/colormap_reference.html>`_
are supported. For example:

.. code-block:: bash

    ifermi plot --property velocity --property-colormap viridis

.. image:: _static/fs-velocity.jpg
    :height: 250px
    :align: center

The minimum and maximum values for the colorbar limits can be set using the ``--cmin``
and ``--cmax`` parameters. These should be used when quantitatively comparing surface
properties between two plots. For example:

.. code-block:: bash

    ifermi plot --property velocity --cmin 0 --cmax 5

As described above, it is also possible calculate the scalar projection of the
face properties onto a cartesian axis using the ``--projection-axis`` option. When
combined with a diverging colormap this can be used to indicate surface properties that
vary between positive and negative numbers. For example, below we color the Fermi
surface of MgB2 by the projection of the group velocity onto the [0 0 1] vector (z-axis).

.. code-block:: bash

    ifermi plot --property velocity --projection-axis 0 0 1 --property-colormap RdBu


.. image:: _static/fs-velocity-projection.jpg
    :height: 250px
    :align: center

Vector valued Fermi surface properties (such as group velocity or spin
magnetisation) can also be visualised as arrows using the ``--vector-property`` option.
If ``--projection-axis`` is set, the color of the arrows will be determined by the
scalar projection of the property vectors onto the specified axis, otherwise the norm
of the projections will be used. The colormap used to color the arrows is specified
using ``--vector-colormap``. Lastly, often it is useful to hide the isosurface
(``--hide-surface`` option) or high-symmetry labels (``-hide-labels``)
when visualising arrows. An example of how to combine these options is given below:

.. code-block:: bash

    ifermi plot --property velocity --projection-axis 0 0 1 --property-colormap RdBu \
                --vector-property --vector-colormap RdBu --hide-surface --hide-labels


.. image:: _static/fs-velocity-arrow.jpg
    :height: 250px
    :align: center

The size of the arrows can be controlled using the ``--vnorm`` parameter. This is
particularly useful when quantitatively comparing vector properties across multiple
Fermi surfaces. A larger ``vnorm`` value will increase the size of the arrows.
The spacing between the arrows is controlled by the ``--vector-spacing`` option. Smaller
values will increase the density of the arrows.

Fermi slices
~~~~~~~~~~~~

IFermi can also generate two-dimensional slices of the Fermi surface along a specified
plane using the ``--slice`` option. Planes are defined by their miller indices (a b c)
and a distance from the plane, d. Most of the above options also apply to to Fermi slices,
however, slices are always plotted using matplotlib as the backend.

For example, a slice through the (0 0 1) plane that passes through the center of the
Brillouin zone (Γ-point) can be generated using:

.. code-block:: bash

    ifermi plot --slice 0 0 1 0

.. image:: _static/slice.png
    :height: 250px
    :align: center

Slices can contain segment properties in the same way that surfaces can contain face
properties. To style slices with projections see :ref:`prop-style`.
When including arrows in Fermi slice figures, only the components of the
arrows in the 2D plane will be shown. As an example below we plot the spin texture of
BiSb (``examples/BiSb``) with and without arrows. The spin texture is colored by the
projection of the spin onto the [0 0 1] cartesian direction.

Without arrows:


.. code-block:: bash

    ifermi plot --mu -0.85  -i 10 --slice 0 0 1 0 --property spin --hide-cell \
                --hide-labels --projection-axis 0 1 0 --property-colormap RdBu

.. image:: _static/slice-property.png
    :height: 250px
    :align: center

With arrows:

.. code-block:: bash

    ifermi plot --mu -0.85  -i 10 --slice 0 0 1 0 --property spin --hide-cell \
                --hide-labels --projection-axis 0 1 0 --property-colormap RdBu \
                --vector-property --vector-colormap RdBu --vnorm 5 --vector-spacing 0.025

.. image:: _static/slice-arrows.png
    :height: 250px
    :align: center

.. WARNING::

    When generating spin texture plots for small regions of k-space, for example,
    in a small area around the Γ-point, it is often necessary to increase the k-point
    mesh density of the underlying DFT calculation. In the example above, the DFT
    calculation was performed on a 21x21x21 k-point mesh.

    Furthermore, projecting the spin magnetisation requires the k-point mesh to cover
    the entire Brillouin zone. I.e., the DFT calculation must have been performed
    without symmetry (``ISYM = - 1`` in VASP).


``ifermi`` reference
----------------------

.. click:: ifermi.cli:cli
  :prog: ifermi
  :nested: full
Index
*****
.. toctree::
   :caption: Usage Guide
   :hidden:

   Introduction <introduction>
   Command-line interface <cli>
   introduction_to_ifermi


.. toctree::
   :caption: Information
   :hidden:

   changelog
   contributors
   contributing
   license
   IFermi on GitHub <https://github.com/fermisurfaces/IFermi>


.. toctree::
   :caption: Function Reference
   :maxdepth: -1
   :hidden:

   API reference <ifermi>
   genindex

.. include:: introduction.rst
============
Introduction
============

.. mdinclude:: ../../README.md
   :start-line: 10
Change log
==========

[Unreleased]
------------

v0.3.0
------

New features:

- Support for plotting individual bands. Specified using the ``--plot-index`` command
  line option. (@aj-searle)

Enhancements:

- Fixed high-symmetry points markers.
- Updated dependencies.

v0.2.6
------

Bug fixes:

- Fixed serialization issues.

v0.2.5
------

Bug fixes:

- Better handling of integer decimation factors.

v0.2.4
------

Enhancements:

- Added function to trim band structure to within a energy cutoff.
- Improved decimation options.

v0.2.3
------

Publish IFermi on zenodo.

v0.2.2
------

Saving interactive html plots is now possible using the plotly backend with:
``ifermi plot --output filename.html``.

v0.2.1
------

Bug fixes:

- Fixed interpolation of projections for 1D slices.
- Fixed position of high-symmetry labels.

v0.2.0
------

This version completely overhauls the Python API and command-line tools. The major
changes are:

- Support for projecting properties onto surface faces and isoline segments. The
  command-line utilities include support for group velocities and spin texture.
- New tools for calculating Fermi surface dimensionality and orientation based on
  the connectivity across periodic boundary conditions.
- New tools for calculating Fermi surface properties such as area and for averaging
  projections across the Fermi surface. This enables the calculation of Fermi velocities.
- New visualisation tools for Fermi surfaces and slices with projections. Fermi surfaces
  can now be colored by the surface properties. Additionally, vector properties
  can be indicated with arrows. This allows for the visualisation of spin texture.

Command line changes:

IFermi now has a new command line interface. There are two subcommands:

- ``ifermi info``: for calculating Fermi surface properties and dimensionalities.
- ``ifermi plot``: for visualisation of Fermi surfaces and slices.

API additions:

- ``FermiSurface`` and ``FermiSlice`` objects now support projections.
- Added ``Isosurface`` and ``Isoline`` classes.
- Added many analysis functions to the ``FermiSurface`` and ``FermiSlice`` modules.
- New ``analysis`` module containing algorithms for:

  - Calculating Fermi surface dimensionality and orientation.
  - Uniformly sampling isosurfaces and isolines.
  - Determining the connectivity of isosurfaces and isolines.
  - Interpolating and smoothing isolines.

API changes:

- ``fermi_surface`` module renamed ``surface``.
- ``FermiSlice`` class and related functions moved to ``slice`` module.
- ``plotter`` module renamed ``plot``.
- ``interpolation`` module renamed ``interpolate``, and ``Interpolator`` class
  renamed ``FourierInterpolator``.

v0.1.5
------

Enhancements:

- Simplified interpolator and FermiSurface generation api.

Bug fixes:

- Fixed bug where the Fermi surface was not exactly centered in reciprocal space.


v0.1.4
------

Enhancements:

- Standardized plots for all plotting backends.
- Added ability to change viewpoint in static plots.
- Documentation overhaul, including new contributors page.
- Added example jupyter notebook.
- API updated to separate plotting and saving files. Allows composing multiple Fermi
  surfaces.
- Surface decimation and smoothing (@mkhorton).
- Support for ``crystal_toolkit`` (@mkhorton).

Bug fixes:

- Fermi level is no longer adjusted from VASP value.
- Bug fix for smoothing (@mdforti).
- Fixed latex labels in plotly (@mdforti).
- Better support for spin polarized materials.

v0.0.4
------

Initial release.
