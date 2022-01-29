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
