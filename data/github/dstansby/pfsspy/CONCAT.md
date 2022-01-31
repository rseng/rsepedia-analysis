![pfsspy](logo/logo_rectangle.png "pfsspy")

[![DOI](https://zenodo.org/badge/163663713.svg)](https://zenodo.org/badge/latestdoi/163663713)
[![Documentation Status](https://readthedocs.org/projects/pfsspy/badge/?version=stable)](https://pfsspy.readthedocs.io/en/stable/?badge=stable)
[![Build Status](https://travis-ci.org/dstansby/pfsspy.svg?branch=master)](https://travis-ci.org/dstansby/pfsspy)
[![codecov](https://codecov.io/gh/dstansby/pfsspy/branch/master/graph/badge.svg)](https://codecov.io/gh/dstansby/pfsspy)
[![CircleCI](https://circleci.com/gh/dstansby/pfsspy.svg?style=svg)](https://circleci.com/gh/dstansby/pfsspy)

Potential Field Source Surface model package for Python. This package is based on original code by Anthony Yeates, which can be found at https://github.com/antyeates1983/pfss.

Documentation
-------------
Full documentation for the package can be found at https://pfsspy.readthedocs.io

Code of conduct
---------------
This is a community package, which has adopted the SunPy code of conduct. Please follow the [code of conduct](CODE_OF_CONDUCT.md) when interacting with others within the pfsspy community.
pfsspy is an SunPy affiliated package, and folows the SunPy code of conduct, shown below.

Code of Conduct
===============

The SunPy community is made up of members from around the globe with a diverse set of skills, personalities, and experiences.
It is through these differences that our community experiences success and continued growth.
We expect everyone in our community to follow these guidelines when interacting with others both inside and outside of our community.
Our goal is to keep ours a positive, inclusive, successful, and growing community.

A member of SunPy is:

Open
----

Members of the community are open to collaboration, whether on patches, reporting issues, asking for help or otherwise.
We welcome those interested in joining the community, and realise that including people with a variety of opinions and backgrounds will only serve to enrich our community.

We are accepting of all who wish to take part in our activities, fostering an environment where anyone can participate and everyone can make a difference, ensuring that all participants are heard and feel confident that they can freely express their opinions.

Considerate
-----------

Members of the community are considerate of their peers -- other developers, users, etc.
We are thoughtful when addressing the efforts of others, keeping in mind that often the labour was completed simply for the good of the community.
We are attentive in our communications, whether in person or online, and we're tactful when approaching differing views.

We recognize the work made by everyone and ensure the proper acknowledgement/citation of original authors at all times.
As authors, we pledge to be explicit about how we want our own work to be cited or acknowledged.

Respectful
----------

Members of the community are respectful.
We are respectful of others, their positions, their skills, their commitments, and their efforts.
We are respectful of the volunteer and professional efforts within the community.
We are respectful of the processes set forth in the community, and we work within them (paying particular attention to those new to the community).

When we disagree, we are courteous in raising our issues.
We provide a harassment- and bullying-free environment, regardless of sex, sexual orientation, gender identity, disability, physical appearance, body size, race, nationality, ethnicity, and religion.
In particular, sexual language and imagery, sexist, racist, or otherwise exclusionary jokes are not appropriate.
We will treat those outside our community with the same respect as people within our community.

Overall, we are good to each other and apply this code of conduct to all community situations online and offline, including mailing lists, forums, social media, conferences, meetings, associated social events, and one-to-one interactions.

We pledge to help the entire community follow the code of conduct, and to not remain silent when we see violations of the code of conduct.
We will take action when members of our community violate this code such as contacting <confidential@sunpy.org> (all emails sent to this address will be treated with the strictest confidence) or talking privately with the person.

Parts of this code of conduct have been adapted from the [astropy](http://www.astropy.org/code_of_conduct.html) and the [Python Software Fundation](https://www.python.org/psf/codeofconduct/) codes of conduct.

------------------------------------------------------------------------

License
-------

The SunPy Community Code of Conduct is licensed under a Creative Commons Attribution 4.0 International License.
We encourage other communities related to ours to use or adapt this code as they see fit.
---
title: 'pfsspy: A Python package for potential field source surface modelling'
tags:
  - Python
  - Astronomy
  - Solar physics
authors:
  - name: David Stansby
    orcid: 0000-0002-1365-1908
    affiliation: 1
  - name: Anthony Yeates
    orcid: 0000-0002-2728-4053
    affiliation: 2
  - name: Samuel T. Badman
    orcid: 0000-0002-6145-436X
    affiliation: "3, 4"
affiliations:
 - name: Mullard Space Science Laboratory, University College London, Holmbury St. Mary, Surrey RH5 6NT, UK
   index: 1
 - name: Department of Mathematical Sciences, Durham University, Durham, DH1 3LE, UK
   index: 2
 - name: Physics Department, University of California, Berkeley, CA 94720-7300, USA
   index: 3
 - name: Space Sciences Laboratory, University of California, Berkeley, CA 94720-7450, USA
   index: 4
date: 2 October 2020
bibliography: paper.bib
---

# Summary
Magnetic fields play a crucial role in the dynamics and evolution of our Sun
and other stars. A common method used to model the magnetic fields in solar and
stellar atmospheres is the potential field source surface (PFSS) model [@Altschuler1969; @Schatten1969].
The PFSS equations assume that there is zero electrical current in the domain
of interest, leading to the equations
\begin{equation}
	\nabla \cdot \mathbf{B} = 0;~~~\nabla \times \mathbf{B} = 0
\end{equation}

These are solved in a spherical shell between the surface of the star and
a configurable outer radius called the 'source surface'. Boundary
conditions are given by the user specified radial component of $\mathbf{B}$ on
the inner boundary and the imposed condition of a purely radial field on the source
surface, which mimics the effect of the escaping stellar wind.

Historically, either custom implementations or the `pfsspack`^[https://www.lmsal.com/~derosa/pfsspack/,
which forms part of the larger `SolarSoft` library for solar physics [@Freeland1998],
written in Interactive Data Language (IDL).]
IDL library have been used to perform PFSS extrapolations. As Python has become a
major programming language within the solar physics and wider astronomy
community [@Bobra2020], there is a need to provide well documented and tested
functionality to perform PFSS extrapolations within the Python ecosystem,
a niche that `pfsspy` fills.


# pfsspy
`pfsspy` is a Python package for solving the PFSS equations, and carrying out
other common related tasks such as tracing magnetic field lines through the
solution, importing various magnetic field data sources, and visualising all of
this data.

The PFSS code implements a finite difference solver, based on the method of
@Ballegooijen2000. Given a 2D map of the radial magnetic field on the inner
boundary, the magnetic vector potential is calculated on a 3D grid equally
spaced in $\sin($latitude$)$, longitude, and $\ln($radius$)$. This method is
tailored in order to achieve $\nabla \times \mathbf{B} = 0$ to machine
precision. More details on the exact numerical scheme are given in the online
documentation^[https://pfsspy.readthedocs.io].


## Integration

`pfsspy` is designed to closely integrate with other packages in the
astronomical and solar physics Python ecosystems. Coordinate aware input and
output maps are created with the sunpy package [@Mumford2020a; @TheSunPyCommunity2020],
and `pfsspy` is fully integrated with the coordinate and unit framework
present in astropy [@TheAstropyCollaboration2018]. This makes it easy to
combine magnetic fields and field lines calculated in `pfsspy` with other data
sources. As an example, \autoref{fig} shows magnetic field lines overplotted
on an extreme-ultraviolet image of a large active region on the Sun.

![An image of the Sun taken by SDO/AIA at 193 angstroms, with selected magnetic field lines traced through a PFSS solution overplotted in white. The PFSS solution and field line tracing were done with `pfsspy`, with a Global Oscillations Network Group (GONG) photospheric magnetogram as input and a source surface at 2.5 solar radii. Although only selected field lines are shown, the magnetic field is solved over the whole Sun.\label{fig}](pfsspy.pdf)

The solar physics community has already made use of `pfsspy` in a number of
works, from interpreting observations from Parker Solar Probe [@Bale2019; @Badman2020],
investigating the structure of coronal mass ejections [@Maguire2020], and
drawing links between the Sun and the solar wind [@Stansby2020]. We hope that
it continues to provide a useful resource for the community in the future.


# Acknowledgements

David Stansby acknowledges STFC grants ST/N504336/1 and ST/S000240/1.
Anthony Yeates acknowledges STFC grant ST/S000321/1.

# References
Citing
------

If you use pfsspy in work that results in publication, please cite the
Journal of Open Source Software paper at https://doi.org/10.21105/joss.02732.
A ready made bibtex entry is

.. code:: bibtex

  @article{Stansby2020,
    doi = {10.21105/joss.02732},
    url = {https://doi.org/10.21105/joss.02732},
    year = {2020},
    publisher = {The Open Journal},
    volume = {5},
    number = {54},
    pages = {2732},
    author = {David Stansby and Anthony Yeates and Samuel T. Badman},
    title = {pfsspy: A Python package for potential field source surface modelling},
    journal = {Journal of Open Source Software}
  }
Using pfsspy
------------
Internals
---------
Finding data
------------
Examples showing how to find, download, and load magnetograms.
Tests
-----
Comparisons of the numerical output of pfsspy to analytic solutions.
Utilities
---------
Useful code that doesn't involve doing a PFSS extrapolation.
Contributing to pfsspy
======================

pfsspy is a community package, and contributions are welcome from anyone. This
includes reporting bugs, requesting new features, along with writing code and
improving documentation.

Reporting Issues
^^^^^^^^^^^^^^^^

If you use pfsspy and stumble upon a problem, the best way to report it is by
opening an `issue`_ on the GitHub issue tracker.
This way we can help you work around the problem and hopefully fix the problem!

When reporting an issue, please try to provide a short description of the issue
with a small code sample, this way we can attempt to reproduce the error.
Also provide any error output generated when you encountered the issue, we can
use this information to debug the issue.

Requesting new features
^^^^^^^^^^^^^^^^^^^^^^^

If there is functionality that is not currently available in pfsspy you can
make a feature request on the `issue`_ page. Please add as much information as
possible regarding the new feature you would like to see.



Improving documentaion
^^^^^^^^^^^^^^^^^^^^^^
If something doesn't make sense in the online documentation, please open an
`issue`_. If you're comfortable fixing it, please open a pull request to the
pfsspy repository.

.. _issue: https://github.com/dstansby/pfsspy/issues

Contributing code
^^^^^^^^^^^^^^^^^
All code contributions to pfsspy are welcome! If you would like to submit a
contribution, please open a pull request at the pfsspy repository. It will then
be reviewed, and if it is a useful addition merged into the main code branch.

If you are new to open source development, the instructions in these links might
be useful for getting started:

* `Astropy Dev Workflow <https://docs.astropy.org/en/latest/development/workflow/development_workflow.html>`_
* `Astropy Dev environment <https://docs.astropy.org/en/latest/development/workflow/get_devel_version.html#get-devel>`_
* `Astropy Pull Request Example <https://docs.astropy.org/en/latest/development/workflow/git_edit_workflow_examples.html#astropy-fix-example>`_
Synoptic map FITS conventions
=============================

FITS is the most common filetype used for the storing of solar images. On this
page the FITS metadata conventions for synoptic maps are collected. All of this
information can be found in, and is taken from,
"Coordinate systems for solar image data (Thompson, 2005)".

=================   ======
   Keyword          Output
-----------------   ------
CRPIX\ *n*          Reference pixel to subtract along axis *n*. Counts from 1 to N.
                    Integer values refer to the centre of the pixel.
CRVAL\ *n*          Coordinate value of the reference pixel along axis *n*.
CDELT\ *n*          Pixel spacing along axis *n*.
CTYPE\ *n*          Coordinate axis label for axis *n*.
PV\ *i*\_\ *m*      Additional parameters needed for some coordinate systems.
=================   ======

Note that *CROTAn* is ignored in this short guide.

Cylindrical equal area projection
---------------------------------
In this projection, the latitude pixels are equally spaced in sin(latitude).
The reference pixel has to be on the equator, to facilitate alignment with the
solar rotation axis.

- CDELT2 is set to 180/pi times the pixel spacing in sin(latitude).
- CTYPE1 is either 'HGLN-CEA' or 'CRLN-CEA'.
- CTYPE2 is either 'HGLT-CEA' or 'CRLT-CEA'.
- PV\ *i*\ _1 is set to 1.
- LONPOLE is 0.

The abbreviations are "Heliographic Longitude - Cylindrical Equal Area" etc.
If the system is heliographic the observer must also be defined in the
metadata.
API reference
=============

.. automodapi:: pfsspy

.. automodapi:: pfsspy.grid

.. automodapi:: pfsspy.fieldline

.. automodapi:: pfsspy.tracing

.. automodapi:: pfsspy.utils

.. automodapi:: pfsspy.analytic
Numerical methods
-----------------

For more information on the numerical methods used in the PFSS calculation see
:download:`this document <manual/numerical_methods.pdf>`.
Installing
----------

pfsspy requires

- python >= 3.8
- sunpy >= 3.0

and can be installed from PyPi using

.. code::

    pip install pfsspy

This will install pfsspy and all of its dependencies. In addition to the core
dependencies, there are two optional dependencies (numba, streamtracer) that
improve code performance. These can be installed with

.. code::

    pip install pfsspy[performance]
Improving performance
---------------------

numba
~~~~~
pfsspy automatically detects an installation of `numba`_, which compiles
some of the numerical code to speed up pfss calculations. To enable this
simply `install numba`_  and use pfsspy as normal.

Streamline tracing
~~~~~~~~~~~~~~~~~~
pfsspy has two streamline tracers: a pure python `pfsspy.tracing.PythonTracer`
and a FORTRAN `pfsspy.tracing.FortranTracer`. The FORTRAN version is
significantly faster, using the `streamtracer`_ package.


.. _numba: https://numba.pydata.org
.. _install numba: http://numba.pydata.org/numba-doc/latest/user/installing.html
.. _streamtracer: https://github.com/dstansby/streamtracer
.. _changelog:

Changelog
=========

1.1.0
-----
New requirements
~~~~~~~~~~~~~~~~
pfsspy now depends on Python >= 3.8, and is officially supported with
Python 3.10

New examples
~~~~~~~~~~~~
A host of new examples comparing pfsspy results to analytic solutions have
been added to the example gallery.

Bug fixes
~~~~~~~~~
- Updated the sunpy package requirement to include all packages needed to use
  sunpy maps.
- Any traced field line points that are out of bounds in latitude (ie. have a
  latitude > 90 deg) are now filtered out. This was previously only an issue
  for very low tracing step sizes.

1.0.1
-----
Bug fixes
~~~~~~~~~
- Fixed compatibility of map validity checks with sunpy 3.1.
- Updated this changelog to make it clear that pfsspy 1.0.0 depends on
  sunpy >= 3.0.

1.0.0
-----

New requirements
~~~~~~~~~~~~~~~~
pfsspy now depends on python >= 3.7, sunpy >=3,
and now does *not* depend on Matplotlib.

New features
~~~~~~~~~~~~
- The ``max_steps`` argument to `pfsspy.tracers.FortranTracer` now defaults to
  ``'auto'`` and automatically sets the maximum number of steps to four times the
  number of steps that are needed to span radially from the solar to source
  surface. ``max_steps`` can still be manually specified as a number if more
  or less steps are desired.
- `~pfsspy.fieldline.FieldLines` now has a ``__len__`` method, meaning one
  can now do ``n_field_lines = len(my_field_lines)``.
- Added :func:`pfsspy.utils.roll_map` to roll a map in the longitude direction.
  This is particularly helpful to modify GONG maps so they have a common
  longitude axis.
- Added the `pfsspy.analytic` sub-module that provides functions to sample
  analytic solutions to the PFSS equations.

Bug fixes
~~~~~~~~~
- :func:`pfsspy.utils.carr_cea_wcs_header` now works with versions of sunpy
  >= 2.0.
- GONG synoptic maps now automatically have their observer information corrected
  (by assuming an Earth observer) when loaded by `sunpy.map.Map`.
- The plot settings of input maps are no longer modified in place.

Breaking changes
~~~~~~~~~~~~~~~~
- The interpretation of the ``step_size`` to `pfsspy.tracers.FortranTracer` has
  been corrected so that it is the step size relative to the radial cell size.
  A step size of 0.01 specified in pfsspy<1.0 is approximately equivalent to a
  step size of 1 in pfsspy 1.0, so you will need to adjust any custom step
  sizes accordingly.
- Any points on field lines that are out of bounds (ie. below the solar surface
  or above the source surface) are now removed by the
  `~pfsspy.tracing.FortranTracer`.
- :func:`pfsspy.pfss` no longer warns if the mean of the input data is non-zero,
  and silently ignores the monopole component.

Removals
~~~~~~~~
- Saving and load PFSS solutions is no longer possible. This was poorly tested,
  and possibly broken. If you have interest in saving and loading being added
  as a new feature to pfsspy, please open a new issue at
  https://github.com/dstansby/pfsspy/issues.

0.6.6
-----
Two bugs have been fixed in `pfsspy.utils.carr_cea_wcs_header`:

  - The reference pixel was previously one pixel too large in both longitude and latitude.
  - The longitude coordinate was previously erroneously translated by one degree.

Both of these are now fixed.

0.6.5
-----
This release improves documentation and handling of HMI maps. In particular:

  - The HMI map downloading example has been updated to use the polar filled
    data product, which does not have any data missing at the poles.
  - :func:`pfsspy.utils.fix_hmi_meta` has been added to fix metadata issues in
    HMI maps. This modifies the metadata of a HMI map to make it FITS compliant,
    allowing it to be used with pfsspy.

0.6.4
-----
This release adds citation information to the documentation.

0.6.3
-----
This release contains the source for the accepted JOSS paper describing pfsspy.

0.6.2
-----
This release includes several small fixes in response to a review of pfsspy
for the Journal of Open Source Software. Thanks to Matthieu Ancellin and
Simon Birrer for their helpful feedback!

- A permanent code of conduct file has been added to the repository.
- Information on how to contribute to pfsspy has been added to the docs.
- The example showing the performance of different magnetic field tracers has
  been fixed.
- The docs are now clearer about optional dependencies that can increase
  performance.
- The GONG example data has been updated due to updated data on the remote
  GONG server.

0.6.1
-----

Bug fixes
~~~~~~~~~

- Fixed some messages in errors raised by functions in `pfsspy.utils`.

0.6.0
-----

New features
~~~~~~~~~~~~
- The `pfsspy.utils` module has been added, and contains various tools for
  loading and working with synoptic maps.
- `pfsspy.Output` has a new `~pfsspy.Output.bunit` property, which returns the
  `~astropy.units.Unit` of the input map.
- Added :meth:`pfsspy.Output.get_bvec`, to sample the magnetic field solution
  at arbitrary coordinates.
- Added the `pfsspy.fieldline.FieldLine.b_along_fline` property, to sample the
  magnetic field along a traced field line.
- Added a guide to the numerical methods used by pfsspy.

Breaking changes
~~~~~~~~~~~~~~~~
- The ``.al`` property of `pfsspy.Output` is now private, as it is not intended
  for user access. If you *really* want to access it, use ``._al`` (but this is
  now private API and there is no guarantee it will stay or return the same thing
  in the future).
- A `ValueError` is now raised if any of the input data to `pfsspy.Input` is
  non-finite or NaN. Previously the PFSS computation would run fine, but the
  output would consist entirely of NaNs.

Behaviour changes
~~~~~~~~~~~~~~~~~
- The monopole term is now ignored in the PFSS calculation. Previously a
  non-zero (but small) monopole term would cause floating point precision issues,
  leading to a very noisy result. Now the monopole term is explicitly removed
  from the calculation. If your input has a non-zero mean value, pfsspy will
  issue a warning about this.
- The data downloaded by the examples is now automatically downloaded and
  cached with `sunpy.data.manager`. This means the files used for running the
  examples will be downloaded and stored in your `sunpy` data directory if
  they are required.
- The observer coordinate information in GONG maps is now automatically set
  to the location of Earth at the time in the map header.

Bug fixes
~~~~~~~~~
- The ``date-obs`` FITS keyword in GONG maps is now correctly populated.

0.5.3
-----
- Improved descriptions in the AIA overplotting example.
- Fixed the 'date-obs' keyword in GONG metadata. Previously this just stored
  the date and not the time; now both the date and time are properly stored.
- Drastically sped up the calculation of source surface and solar surface
  magnetic field footpoints.

0.5.2
-----
- Fixed a bug in the GONG synoptic map source where a map failed to load once
  it had already been loaded once.

0.5.1
-----
- Fixed some calculations in ``pfsspy.carr_cea_wcs_header``, and clarified in the
  docstring that the input shape must be in ``[nlon, nlat]`` order.
- Added validation to `pfsspy.Input` to check that the inputted map covers the
  whole solar surface.
- Removed ghost cells from `pfsspy.Output.bc`. This changes the shape of the
  returned arrays by one along some axes.
- Corrected the shape of `pfsspy.Output.bg` in the docstring.
- Added an example showing how to load ADAPT ensemble maps into a
  `~sunpy.map.CompositeMap`
- Sped up field line expansion factor calculations.

0.5.0
-----

Changes to outputted maps
~~~~~~~~~~~~~~~~~~~~~~~~~
This release largely sees a transition to leveraging Sunpy Map objects. As such,
the following changes have been made:

`pfsspy.Input` now *must* take a `sunpy.map.GenericMap` as an
input boundary condition (as opposed to a numpy array). To convert a numpy array
to a `~sunpy.map.GenericMap`, the helper function
``pfsspy.carr_cea_wcs_header`` can be used::

  map_date = datetime(...)
  br = np.array(...)
  header = pfsspy.carr_cea_wcs_header(map_date, br.shape)

  m = sunpy.map.Map((br, header))
  pfss_input = pfsspy.Input(m, ...)


`pfsspy.Output.source_surface_br` now returns a `~sunpy.map.GenericMap`
instead of an array. To get the data array use ``source_surface_br.data``.

The new `pfsspy.Output.source_surface_pils` returns the coordinates of
the polarity inversion lines on the source surface.

In favour of directly using the plotting functionality built into SunPy,
the following plotting functionality has been removed:

- ``pfsspy.Input.plot_input``. Instead `~pfsspy.Input` has a new
  `~pfsspy.Input.map`  property, which returns a SunPy map, which can easily
  be plotted using `sunpy.map.GenericMap.plot`.
- ``pfsspy.Output.plot_source_surface``. A map of :math:`B_{r}` on the source
  surface can now be obtained using `pfsspy.Output.source_surface_br`, which
  again returns a SunPy map.
- ``pfsspy.Output.plot_pil``. The coordinates of the polarity inversion lines
  on the source surface can now be obtained using
  `pfsspy.Output.source_surface_pils`, which can then be plotted using
  ``ax.plot_coord(pil[0])`` etc. See the examples section for an example.

Specifying tracing seeds
~~~~~~~~~~~~~~~~~~~~~~~~
In order to make specifying seeds easier, they must now be a
`~astropy.coordinates.SkyCoord` object. The coordinates are internally
transformed to the Carrington frame of the PFSS solution, and then traced.

This should make specifying coordinates easier, as lon/lat/r coordinates can
be created using::

  seeds = astropy.coordinates.SkyCoord(lon, lat, r, frame=output.coordinate_frame)

To convert from the old x, y, z array used for seeds, do::

  r, lat, lon = pfsspy.coords.cart2sph
  r = r * astropy.constants.R_sun
  lat = (lat - np.pi / 2) * u.rad
  lon = lon * u.rad

  seeds = astropy.coordinates.SkyCoord(lon, lat, r, frame=output.coordinate_frame)

Note that the latitude must be in the range :math:`[-\pi/2, \pi/2]`.

GONG and ADAPT map sources
~~~~~~~~~~~~~~~~~~~~~~~~~~
pfsspy now comes with built in `sunpy` map sources for GONG and ADAPT synoptic
maps, which automatically fix some non-compliant FITS header values. To use
these, just import ``pfsspy`` and load the .FITS files as normal with sunpy.

Tracing seeds
~~~~~~~~~~~~~
`pfsspy.tracing.Tracer` no longer has a ``transform_seeds`` helper method, which
has been replaced by `~pfsspy.tracing.Tracer.coords_to_xyz` and
``pfsspy.tracing.Tracer.xyz_to_coords``. These new methods convert
between `~astropy.coordinates.SkyCoord` objects, and Cartesian xyz coordinates
of the internal magnetic field grid.

0.4.3
-----

- Improved the error thrown when trying to use
  :class`pfsspy.tracing.FotranTracer` without the ``streamtracer`` module
  installed.
- Fixed some layout issues in the documentation.

0.4.2
-----

- Fix a bug where :class`pfsspy.tracing.FotranTracer` would overwrite the
  magnetic field values in an `~pfsspy.Output` each time it was used.

0.4.1
-----

- Reduced the default step size for the `~pfsspy.tracing.FortranTracer`
  from 0.1 to 0.01 to give more resolved field lines by default.

0.4.0
-----

New fortran field line tracer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:mod:`pfsspy.tracing` contains a new tracer,
`~pfsspy.tracing.FortranTracer`. This requires and uses the
`streamtracer <https://streamtracer.readthedocs.io/en/stable/>`_ package
which does streamline tracing rapidly in python-wrapped
fortran code. For large numbers of field lines this results in an ~50x
speedup compared to the `~pfsspy.tracing.PythonTracer`.

Changing existing code to use the new tracer is as easy as swapping out
``tracer = pfsspy.tracer.PythonTracer()`` for
``tracer = pfsspy.tracer.FortranTracer()``. If you notice any issues with the
new tracer, please report them at https://github.com/dstansby/pfsspy/issues.

Changes to field line objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``pfsspy.FieldLines`` and ``pfsspy.FieldLine`` have moved to
  `pfsspy.fieldline.FieldLines` and
  `pfsspy.fieldline.FieldLine`.
- `~pfsspy.fieldline.FieldLines` no longer has ``source_surface_feet``
  and ``solar_feet`` properties. Instead these have moved to the new
  `pfsspy.fieldline.OpenFieldLines` class. All the open field lines
  can be accessed from a `~pfsspy.fieldline.FieldLines` instance using
  the new `~pfsspy.fieldline.FieldLines.open_field_lines`
  property.

Changes to `~pfsspy.Output`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- `pfsspy.Output.bg` is now returned as a 4D array instead of three 3D
  arrays. The final index now indexes the vector components; see the docstring
  for more information.

0.3.2
-----
- Fixed a bug in ``pfsspy.FieldLine.is_open``, where some open field lines
  were incorrectly calculated to be closed.

0.3.1
-----
- Fixed a bug that incorrectly set closed line field polarities to -1 or 1
  (instead of the correct value of zero).
- ``FieldLine.footpoints`` has been removed in favour of the new
  ``pfsspy.FieldLine.solar_footpoint`` and
  ``pfsspy.FieldLine.source_surface_footpoint``. These each return a single
  footpoint. For a closed field line, see the API docs for further details
  on this.
- ``pfsspy.FieldLines`` has been added, as a convenience class to store a
  collection of field lines. This means convenience attributes such as
  ``pfsspy.FieldLines.source_surface_feet`` can be used, and their values are
  cached greatly speeding up repeated use.

0.3.0
-----

- The API for doing magnetic field tracing has changed.
  The new :mod:`pfsspy.tracing` module contains `~pfsspy.tracing.Tracer`
  classes that are used to perform the tracing. Code needs to be changed from::

    fline = output.trace(x0)

  to::

    tracer = pfsspy.tracing.PythonTracer()
    tracer.trace(x0, output)
    flines = tracer.xs

  Additionally ``x0`` can be a 2D array that contains multiple seed
  points to trace, taking advantage of the parallelism of some solvers.
- The ``pfsspy.FieldLine`` class no longer inherits from
  `~astropy.coordinates.SkyCoord`, but the
  `~astropy.coordinates.SkyCoord` coordinates are now stored in
  ``pfsspy.FieldLine.coords`` attribute.
- ``pfsspy.FieldLine.expansion_factor`` now returns ``np.nan`` instead of
  ``None`` if the field line is closed.
- ``pfsspy.FieldLine`` now has a ``~pfsspy.FieldLine.footpoints``
  attribute that returns the footpoint(s) of the field line.

0.2.0
-----

- `pfsspy.Input` and `pfsspy.Output` now take the optional keyword
  argument *dtime*, which stores the datetime on which the magnetic field
  measurements were made. This is then propagated to the *obstime* attribute
  of computed field lines, allowing them to be transformed in to coordinate
  systems other than Carrington frames.
- ``pfsspy.FieldLine`` no longer overrrides the SkyCoord ``__init__``;
  this should not matter to users, as FieldLine objects are constructed
  internally by calling `pfsspy.Output.trace`

0.1.5
-----

- ``Output.plot_source_surface`` now accepts keyword arguments that are given to
  Matplotlib to control the plotting of the source surface.

0.1.4
-----

- Added more explanatory comments to the examples
- Corrected the dipole solution calculation
- Added ``pfsspy.coords.sph2cart`` to transform from spherical to cartesian
  coordinates.

0.1.3
-----

- ``pfsspy.Output.plot_pil`` now accepts keyword arguments that are given
  to Matplotlib to control the style of the contour.
- ``pfsspy.FieldLine.expansion_factor`` is now cached, and is only
  calculated once if accessed multiple times.
.. image:: ../../logo/logo_rectangle.svg

pfsspy documentation
====================
pfsspy is a python package for carrying out Potential Field Source Surface
modelling, a commonly used magnetic field model of the Sun and other stars.

.. note::
  If you find any bugs or have any suggestions for improvement, please raise
  an issue here: https://github.com/dstansby/pfsspy/issues



Contents
--------

.. toctree::
   :maxdepth: 1

   installing
   auto_examples/index
   pfsspy
   performance
   changes
   contributing
   numerical_method
   synoptic_fits

.. include:: ../../CITING.rst
