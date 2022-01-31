SWIFTsimIO
==========

![Build Status](https://github.com/swiftsim/swiftsimio/actions/workflows/pytest.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/swiftsimio/badge/?version=latest)](https://swiftsimio.readthedocs.io/en/latest/?badge=latest)
[![JOSS Status](https://joss.theoj.org/papers/e85c85f49b99389d98f9b6d81f090331/status.svg)](https://joss.theoj.org/papers/e85c85f49b99389d98f9b6d81f090331)


The SWIFT astrophysical simulation code (http://swift.dur.ac.uk) is used
widely. There exists many ways of reading the data from SWIFT, which outputs
HDF5 files. These range from reading directly using `h5py` to using a complex
system such as `yt`; however these either are unsatisfactory (e.g. a lack of
unit information in reading HDF5), or too complex for most use-cases. 
`swiftsimio` provides an object-oriented API to read (dynamically) data
from SWIFT.

Full documentation is available at [ReadTheDocs](http://swiftsimio.readthedocs.org).

Getting set up with `swiftsimio` is easy; it (by design) has very few
requirements. There are a number of optional packages that you can install
to make the experience better and these are recommended.


Requirements
------------

This requires `python` `v3.6.0` or higher. Unfortunately it is not
possible to support `swiftsimio` on versions of python lower than this.
It is important that you upgrade if you are still a `python2` user.

### Python packages


+ `numpy`, required for the core numerical routines.
+ `h5py`, required to read data from the SWIFT HDF5 output files.
+ `unyt`, required for symbolic unit calculations (depends on sympy`).

### Optional packages


+ `numba`, highly recommended should you wish to use the in-built visualisation
  tools.
+ `scipy`, required if you wish to generate smoothing lengths for particle types
  that do not store this variable in the snapshots (e.g. dark matter)
+ `tqdm`, required for progress bars for some long-running tasks. If not installed
  no progress bar will be shown.
+ `py-sphviewer`, if you wish to use our integration with this visualisation
  code.


Installing
----------

`swiftsimio` can be installed using the python packaging manager, `pip`,
or any other packaging manager that you wish to use:

`pip install swiftsimio`


Citing
------

Please cite `swiftsimio` using the JOSS [paper](https://joss.theoj.org/papers/10.21105/joss.02430):

```bibtex
@article{Borrow2020,
  doi = {10.21105/joss.02430},
  url = {https://doi.org/10.21105/joss.02430},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2430},
  author = {Josh Borrow and Alexei Borrisov},
  title = {swiftsimio: A Python library for reading SWIFT data},
  journal = {Journal of Open Source Software}
}
```

If you use any of the subsampled projection backends, we ask that you cite our
relevant SPHERIC [paper](https://arxiv.org/abs/2106.05281). Note that citing
the arXiv version here is recommended as the ADS cannot track conference
proceedings well.

```bibtex
@article{Borrow2021
  title={Projecting SPH Particles in Adaptive Environments}, 
  author={Josh Borrow and Ashley J. Kelly},
  year={2021},
  eprint={2106.05281},
  archivePrefix={arXiv},
  primaryClass={astro-ph.GA}
}
```
Contributing to SWIFTsimIO
==========================

Contributions for SWIFTsimIO should come through our GitHub repository,
available at https://github.com/swiftsim/swiftsimio.

Contributions are always welcome, but you should make sure of the following:

+ Your contributions pass all unit tests (you can check this with `pytest`)
+ Your contributions add unit tests for new functionality
+ Your contributions are formatted with the `black` formatter (see `format.sh`)
+ Your contributions are documented fully under `/docs`.

You should also abide by the following code of conduct:

### Code of Conduct

The community of participants in open source Astronomy projects is made up of
members from around the globe with a diverse set of skills, personalities,
and experiences. It is through these differences that our community
experiences success and continued growth. We expect everyone in our community
to follow these guidelines when interacting with others both inside and
outside of our community. Our goal is to keep ours a positive, inclusive,
successful, and growing community.

As members of the community,

+ We pledge to treat all people with respect and provide a harassment- and
  bullying-free environment, regardless of sex, sexual orientation and/or
  gender identity, disability, physical appearance, body size, race,
  nationality, ethnicity, and religion. In particular, sexual language and
  imagery, sexist, racist, or otherwise exclusionary jokes are not appropriate.
+ We pledge to respect the work of others by recognizing
  acknowledgement/citation requests of original authors. As authors, we pledge
  to be explicit about how we want our own work to be cited or acknowledged.
+ We pledge to welcome those interested in joining the community, and realize
  that including people with a variety of opinions and backgrounds will only
  serve to enrich our community. In particular, discussions relating to
  pros/cons of various technologies, programming languages, and so on are
  welcome, but these should be done with respect, taking proactive measure to
  ensure that all participants are heard and feel confident that they can
  freely express their opinions.
+ We pledge to welcome questions and answer them respectfully, paying
  particular attention to those new to the community. We pledge to provide
  respectful criticisms and feedback in forums, especially in discussion
  threads resulting from code contributions.
+ We pledge to be conscientious of the perceptions of the wider community and
  to respond to criticism respectfully. We will strive to model behaviours that
  encourage productive debate and disagreement, both within our community and
  where we are criticized. We will treat those outside our community with the
  same respect as people within our community.
+ We pledge to help the entire community follow the code of conduct, and to
  not remain silent when we see violations of the code of conduct. We will
  take action when members of our community violate this code such as
  contacting joshua.borrow@durham.ac.uk with the subject line SWIFTsimIO Code
  of Conduct (all emails sent in this fashion will be treated with the
  strictest confidence) or talking privately with the person.
+ This code of conduct applies to all community situations online and
  offline, including mailing lists, forums, social media, conferences,
  meetings, associated social events, and one-to-one interactions.

Any related activity or project organized by members of the SWIFTsimIO
community, including affiliated packages, are welcome to have their own codes
of conduct, but agree to also abide by the present code of conduct.

Parts of this code of conduct have been adapted from the PSF code of conduct and
the Astropy code of conduct: https://www.astropy.org/code_of_conduct.html... SWIFTsimIO documentation master file, created by
   sphinx-quickstart on Sat Nov 23 15:40:41 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SWIFTsimIO's documentation!
======================================

``swiftsimio`` is a toolkit for reading SWIFT_ data, an astrophysics
simulation code. It is used to ensure that everything that you read has a
symbolic unit attached, and  can be used for visualisation. The final key
feature that it enables is the use of the cell metadata in ``SWIFT``
snapshots to enable partial reading.

.. toctree::
   :maxdepth: 2

   getting_started/index
   loading_data/index
   masking/index
   visualisation/index
   velociraptor/index
   creating_initial_conditions/index
   statistics/index
   command_line/index

   modules/index


Citing SWIFTsimIO
=================

Please cite ``swiftsimio`` using the JOSS paper_:

.. code-block:: bibtex

    @article{Borrow2020,
      doi = {10.21105/joss.02430},
      url = {https://doi.org/10.21105/joss.02430},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {52},
      pages = {2430},
      author = {Josh Borrow and Alexei Borrisov},
      title = {swiftsimio: A Python library for reading SWIFT data},
      journal = {Journal of Open Source Software}
    }
    
If you use any of the subsampled projection backends, we ask that you cite our relevant
SPHERIC article_. Note that citing the arXiv version here is recommended as the ADS
cannot track conference proceedings well.

.. code-block:: bibtex

    @article{Borrow2021,
      title={Projecting SPH Particles in Adaptive Environments}, 
      author={Josh Borrow and Ashley J. Kelly},
      year={2021},
      eprint={2106.05281},
      archivePrefix={arXiv},
      primaryClass={astro-ph.GA}
    }


.. _SWIFT: http://www.swiftsim.com
.. _paper: https://joss.theoj.org/papers/10.21105/joss.02430
.. _article: https://ui.adsabs.harvard.edu/abs/2021arXiv210605281B/abstract

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
VELOCIraptor Integration
========================

:mod:`swiftsimio` can be used with the :mod:`velociraptor` library
to extract the particles contained within a given halo and its surrounding
region.

The :mod:`velociraptor` library has documentation also available on
ReadTheDocs `here <http://velociraptor-python.readthedocs.org/>`_. It can
be installed from PyPI using ``pip install velociraptor``.

The overarching workflow for this integration is as follows:

+ Load the halo catalogue and groups file using the :mod:`velociraptor`
  module.
+ Get two objects, corresponding to the bound and unbound particles,
  for a halo.
+ Use the `to_swiftsimio_dataset` to load the region around the halo
  with our ahead-of-time masking technique.
+ Use the region around the halo directly, or use the mask provided
  for each particle type to only consider bound particles.

This workflow is explored below. You can use the example data available
below if you do not have any SWIFT and VELOCIraptor data available.

``http://virgodb.cosma.dur.ac.uk/swift-webstorage/IOExamples/small_cosmo_volume.zip``

Example
-------

First, we must load the VELOCIraptor catalogue as follows:

.. code-block:: python

    from velociraptor import load as load_catalogue
    from velociraptor.particles import load_groups

    catalogue_name = "velociraptor"
    snapshot_name = "snapshot"

    catalogue = load_catalogue(f"{catalogue_name}.properties")
    groups = load_groups(f"{catalogue_name}.catalog_groups", catalogue=catalogue)


Then, to extract the largest halo in the volume

.. code-block:: python

    particles, unbound_particles = groups.extract_halo(halo_id=0)

To load the particles to a :mod:`swiftsimio` dataset,

.. code-block:: python

    from velociraptor.swift.swift import to_swiftsimio_dataset

    data, mask = to_swiftsimio_dataset(
        particles,
        f"{snapshot_name}.hdf5",
        generate_extra_mask=True
    )

with the ``generate_extra_mask`` providing the second return value which
is a mask to extract only the bound particles in the system.

Making an image of the full box shows that only a small subsection of the
volume has been loaded (those within twice the maximal usable radius within
VELOCIraptor)

.. image:: load_halo_fullbox.png

The code for making this image is as follows:

.. code-block:: python

    from swiftsimio.visualisation import project_gas_pixel_grid
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    
    grid = project_gas_pixel_grid(data=data, resolution=1024)
    
    fig, ax = plt.subplots(figsize=(4, 4), dpi=1024 // 4)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.axis("off")
    ax.imshow(grid.T, origin="lower", cmap="inferno", norm=LogNorm(vmin=1e4, clip=True))
    fig.savefig("load_halo_fullbox.png")

To make an image of just the central halo, we can access properties on the
``particles`` instance to get the position of the halo.

.. code-block:: python

    region = [
        particles.x_mbp - particles.r_200crit, particles.x_mbp + particles.r_200crit,
        particles.y_mbp - particles.r_200crit, particles.y_mbp + particles.r_200crit,
    ]

    grid = project_gas_pixel_grid(data=data, resolution=1024, region=region)

    fig, ax = plt.subplots(figsize=(4, 4), dpi=1024 // 4)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.axis("off")
    ax.imshow(grid.T, origin="lower", cmap="inferno", norm=LogNorm(vmin=1e4, clip=True))
    fig.savefig("load_halo_selection.png")

This produces the following image:

.. image:: load_halo_selection.png

Then, finally, we can visualise only the bound particles, through the use of the ``mask``
object that was returned when we initially extracted the ``swiftsimio`` dataset:

.. code-block:: python

    grid = project_gas_pixel_grid(data=data, resolution=1024, region=region, mask=mask.gas)

    fig, ax = plt.subplots(figsize=(4, 4), dpi=1024 // 4)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.axis("off")
    ax.imshow(grid.T, origin="lower", cmap="inferno", norm=LogNorm(vmin=1e4, clip=True))
    fig.savefig("load_halo_bound_selection.png")

Producing the following image:

.. image:: load_halo_bound_selection.png

Hopefully, when you use this feature, you have more exciting data to use than the
as-small-as-possible example that we show here!Command-line Utilities
======================

:mod:`swiftsimio` comes with some useful command-line utilities.
Basic documentation for these is provided below, but you can always
find up-to-date documentation by invoking these with ``-h`` or
``--help``.

``swiftsnap``
-------------

The ``swiftsnap`` utility, introduced in :mod:`swiftsimio` version
3.1.2, allows you to preview the metadata inside a SWIFT snapshot
file. Simply invoke it with the path to a snapshot, and it will
show you a selection of useful metadata. See below for an example.

.. code:: bash

    swiftsnap output_0103.hdf5

Produces the following output:

.. code::

    Untitled SWIFT simulation
    Written at: 2020-06-01 08:44:51
    Active policies: cosmological integration, hydro, keep, self gravity, steal
    Output type: Snapshot, Output selection: Snapshot
    LLVM/Clang (11.0.0)
    Non-MPI version of SWIFT
    SWIFT (io_selection_changes)
    v0.8.5-725-g10d7d5b3-dirty
    2020-05-29 18:00:58 +0100
    Simulation state: z=0.8889, a=0.5294, t=6.421 Gyr
    H_0=70.3 km/(Mpc*s), ρ_crit=1.433e-05 cm**(-3)
    Ω_b=0.0455, Ω_k=0, Ω_lambda=0.724, Ω_m=0.276, Ω_r=0
    ω=-1, ω_0=-1, ω_a=0
    Gravity scheme: With per-particle softening
    Hydrodynamics scheme: Gadget-2 version of SPH (Springel 2005)
    Chemistry model: None
    Cooling model: None
    Entropy floor: None
    Feedback model: None
    Tracers: None

Loading Data
============

The main purpose of :mod:`swiftsimio` is to load data. This section will tell
you all about four main objects:

+ :obj:`swiftsimio.reader.SWIFTUnits`, responsible for creating a correspondence between
  the SWIFT units and :mod:`unyt` objects.
+ :obj:`swiftsimio.reader.SWIFTMetadata`, responsible for loading any required information
  from the SWIFT headers into python-readable data.
+ :obj:`swiftsimio.reader.SWIFTDataset`, responsible for holding all particle data, and
  keeping track of the above two objects.
+ :obj:`swiftsimio.reader.SWIFTParticleTypeMetadata`, responsible for
  cataloguing metadata just about individual particle types, like gas,
  including what particle fields are present.


To get started, first locate any SWIFT data that you wish to analyse. If you
don't have any handy, you can always download our test cosmological volume
at:

``http://virgodb.cosma.dur.ac.uk/swift-webstorage/IOExamples/cosmo_volume_example.hdf5``

with associated halo catalogue at

``http://virgodb.cosma.dur.ac.uk/swift-webstorage/IOExamples/cosmo_volume_example.properties``

which is needed should you wish to use the ``velociraptor`` integration library in the
visualisation examples.

To create your first instance of :obj:`swiftsimio.reader.SWIFTDataset`, you can
use the helper function :mod:`swiftsimio.load` as follows:

.. code-block:: python

   from swiftsimio import load

   # Of course, replace this path with your own snapshot should you be using
   # custom data.
   data = load("cosmo_volume_example.hdf5")

The type of ``data`` is now :obj:`swiftsimio.reader.SWIFTDataset`. Have a
quick look around this dataset in an ``iPython`` shell, or a ``jupyter``
notebook, and you will see that it contains several sub-objects:

+ ``data.gas``, which contains all information about gas particles in the
  simulation.
+ ``data.dark_matter``, likewise containing information about the dark matter
  particles in the simulation.
+ ``data.metadata``, an instance of :obj:`swiftsimio.reader.SWIFTMetadata`
+ ``data.units``, an instance of :obj:`swiftsimio.reader.SWIFTUnits`

Using metadata
--------------

Let's begin by reading some useful metadata straight out of our
``data.metadata`` object. For instance, we may want to know the box-size of
our simulation:

.. code-block:: python

   meta = data.metadata
   boxsize = meta.boxsize

   print(boxsize)

This will output ``[142.24751067 142.24751067 142.24751067] Mpc`` - note
the units that are attached. These units being attached to everything is one
of the key advantages of using :mod:`swiftsimio`. It is really easy to convert
between units; for instance if we want that box-size in kiloparsecs,

.. code-block:: python

   boxsize.convert_to_units("kpc")

   print(boxsize)

Now outputting ``[142247.5106242 142247.5106242 142247.5106242] kpc``. Neat!
This is all thanks to our tight integration with :mod:`unyt`. If you have more
complex units, it is often useful to specify them in terms of :mod:`unyt`
objects as follows:

.. code-block:: python

   import unyt

   new_units = unyt.cm * unyt.Mpc / unyt.kpc
   new_units.simplify()

   boxsize.convert_to_units(new_units)

In general, we suggest using :mod:`unyt` unit objects rather than strings. You
can find more information about :mod:`unyt` on the `unyt documentation website`_.

.. _`unyt documentation website`: https://unyt.readthedocs.io/en/stable/

There is lots of metadata available through this object; the best way to see
this is by exploring the object using ``dir()`` in an interactive shell, but
as a summary:

+ All metadata from the snapshot is available through many variables, for example
  the ``meta.hydro_scheme`` property.
+ The numbers of particles of different types are available through
  ``meta.n_{gas,dark_matter,stars,black_holes}``.
+ Several pre-LaTeXed strings are available describing the configuration state
  of the code, such as ``meta.hydro_info``, ``meta.compiler_info``.
+ Several snapshot-wide parameters, such as ``meta.a`` (current scale factor),
  ``meta.t`` (current time), ``meta.z`` (current redshift), ``meta.run_name``
  (the name of this run, specified in the SWIFT parameter file), and
  ``meta.snapshot_date`` (a :mod:`datetime` object describing when the
  snapshot was written to disk).
+ If you have ``astropy`` installed, you can also use the ``metadata.cosmology``
  object, which is an ``astropy.cosmology.w0waCDM`` instance.


Reading particle data
---------------------

To find out what particle properties are present in our snapshot, we can use
the instance of :obj:`swiftsimio.reader.SWIFTMetadata`, ``data.metadata``,
which contains several instances of
:obj:`swiftsimio.reader.SWIFTParticleTypeMetadata` describing what kinds of
fields are present in gas or dark matter:

.. code-block:: python

   # Note that gas_properties is an instance of SWIFTParticleTypeMetadata
   print(data.metadata.gas_properties.field_names)

This will print a large list, like

.. code-block:: python

   ['coordinates',
   'densities',
   ...
   'temperatures',
   'velocities']

These individual attributes can be accessed through the object-oriented
interface, for instance,

.. code-block:: python

   x_gas = data.gas.coordinates
   rho_gas = data.gas.densities
   x_dm = data.dark_matter.coordinates

Only at this point is any information actually read from the snapshot, so far
we have only read three arrays into memory - in this case corresponding to
``/PartType0/Coordinates``, ``/PartType1/Coordinates``, and
``/PartType0/Densities``.

This allows you to be quite lazy when writing scripts; you do not have to
write, for instance, a block at the start of the file with a
``with h5py.File(...) as handle:`` and read all of the data at once, you can
simply access data whenever you need it through this predictable interface.

Just like the boxsize, these carry symbolic :mod:`unyt` units,

.. code-block:: python

   print(x_gas.units)

will output ``Mpc``. We can again convert these to whatever units
we like. For instance, should we wish to convert our gas densities to solar
masses per cubic megaparsec,

.. code-block:: python

   new_density_units = unyt.Solar_Mass / unyt.Mpc**3

   rho_gas.convert_to_units(new_density_units)

   print(rho_gas.units.latex_repr)

which will output ``'\\frac{M_\\odot}{\\rm{Mpc}^{3}}'``. This is a LaTeX
representation of those symbolic units that we just converted our data to -
this is very useful for making plots as it can ensure that your data and axes
labels always have consistent units.


Named columns
-------------

SWIFT can output custom metadata in ``SubgridScheme/NamedColumns`` for multi
dimensional tables containing columns that carry individual data. One common
example of this is the element mass fractions of gas and stellar particles.
These are then placed in an object hierarchy, as follows:

.. code-block:: python

   print(data.gas.element_mass_fractions)


This will output: Named columns instance with ['hydrogen', 'helium',
'carbon', 'nitrogen', 'oxygen', 'neon', 'magnesium', 'silicon', 'iron']
available for "Fractions of the particles' masses that are in the given
element"

Then, to access individual columns (in this case element abundances):

.. code-block:: python

   # Access the silicon abundance
   data.gas.element_mass_fractions.silicon


Non-unyt properties
-------------------

Each data array has some custom properties that are not present within the base
:obj:`unyt.unyt_array` class. We create our own version of this in
:obj:`swiftsimio.objects.cosmo_array`, which allows each dataset to contain
its own cosmology and name properties.

For instance, should you ever need to know what a dataset represents, you can
ask for a description:

.. code-block:: python

   print(rho_gas.name)

which will output ``Co-moving mass densities of the particles``. They include
scale-factor information, too, through the ``cosmo_factor`` object,

.. code-block:: python

   # Conversion factor to make the densities a physical quantity
   print(rho_gas.cosmo_factor.a_factor)
   physical_rho_gas = rho_gas.cosmo_factor.a_factor * rho_gas

   # Symbolic scale-factor expression
   print(rho_gas.cosmo_factor.expr)

which will output ``132651.002785671`` and ``a**(-3.0)``. This is an easy way
to convert your co-moving values to physical ones.

An even easier way to convert your properties to physical is to use the
built-in ``to_physical`` and ``convert_to_physical`` methods, as follows:

.. code-block:: python

   physical_rho_gas = rho_gas.to_physical()

   # Convert in-place
   rho_gas.convert_to_physical()


User-defined particle types
---------------------------

It is now possible to add user-defined particle types that are not already
present in the :mod:`swiftsimio` metadata. All you need to do is specify the
three names (see below) and then the particle datasets that you have provided
in SWIFT will be automatically read.

.. code-block:: python

   import swiftsimio as sw
   import swiftsimio.metadata.particle as swp
   from swiftsimio.objects import cosmo_factor, a

   swp.particle_name_underscores[6] = "extratype"
   swp.particle_name_class[6] = "Extratype"
   swp.particle_name_text[6] = "Extratype"

   data = sw.load(
       "extra_test.hdf5",
   )
Tools
=====

:mod:`swiftsimio` includes a few tools to help you make your visualisations
'prettier'. Below we describe these tools and their use.

2D Color Maps
-------------

The :mod:`swiftsimio.visualisation.tools.cmaps` module includes three
objects that can be used to deploy two dimensional colour maps. The first,
:class:`swiftsimio.visualisation.tools.cmaps.LinearSegmentedCmap2D`, and second
:class:`swiftsimio.visualisation.tools.cmaps.LinearSegmentedCmap2DHSV`, allow
you to generate new color maps from sets of colors and coordinates.

.. code-block:: python

	bower = LinearSegmentedCmap2D(
		colors=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]],
		coordinates=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
		name="bower"
	)
	
This generates a color map that is a quasi-linear interpolation between all
of the points. The map can be displayed using the ``plot`` method,

.. code-block:: python

	fig, ax = plt.subplots()
	
	bower.plot(ax)
	
Which generates:

.. image:: bower_cmap.png

Finally, the color map can be applied to data by calling it:

.. code-block:: python

	def vertical_func(x):
		return abs(1.0 - 2.0 * x)
	
	def horizontal_func(y):
		return y ** 2
	
	raster_at = np.linspace(0, 1, 1024)
	
	x, y = np.meshgrid(horizontal_func(raster_at), vertical_func(raster_at))
	
	imaged = bower(x, y)
	
	plt.imsave("test_2d_cmap_output.png", imaged)
	
Where here ``imaged`` is an RGBA array. This outputs:

.. image:: test_2d_cmap_output.png

The final type of 2D color map is loaded from an image, such as the one displayed
below which is similar to the famous color map used for the Millenium simulation.

.. image:: millenium_cmap.png

This can be loaded using the
:class:`swiftsimio.visualisation.tools.cmaps.ImageCmap2D` class, as follows:

.. code-block:: python
	
	mill = ImageCmap2D(filename="millenium_cmap.png")
	
and can be used similarly to the other color maps. For the example above, this
outputs the following:

.. image:: test_2d_cmap_output_mill.png

This is the recommended way to use two dimensional color maps, as their
generation can be quite complex and best left to image-generation programs
such as GIMP or Photoshop.
Py-SPHViewer Integration
========================

We provide a wrapper of the ever-popular py-sphviewer_ for easy use with
:mod:`swiftsimio` datasets. Particle datasets that do not contain smoothing
lengths will have them generated through the use of the scipy ``cKDTree``.
You can get access to the objects through a sub-module as follows:

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.sphviewer import SPHViewerWrapper
   import matplotlib.pyplot as plt

   data = load("cosmo_volume_example.hdf5")

   resolution = 2048

   gas = SPHViewerWrapper(data.gas, smooth_over="masses")
   gas_temp = SPHViewerWrapper(
       data.gas,
       smooth_over=data.gas.masses * data.gas.temperatures
   )
   dark_matter = SPHViewerWrapper(data.dark_matter, smooth_over="masses")

   gas.quick_view(xsize=resolution, ysize=resolution, r="infinity")
   gas_temp.quick_view(xsize=resolution, ysize=resolution, r="infinity")
   dark_matter.quick_view(xsize=resolution, ysize=resolution, r="infinity")

   plt.imsave("gas_image.png", gas.image)
   plt.imsave("gas_temp.png", gas_temp.image / gas.image)
   plt.imsave("dm_image.png", dark_matter.image)


The :obj:`swiftsimio.visualisation.sphviewer.SPHViewerWrapper` object allows you
to get access to the particles, camera, and render object through ``.particles``,
``.get_camera()`` and ``.camera``, and ``.get_render()`` and ``.render``
respectively.

.. _py-sphviewer: https://github.com/alejandrobll/py-sphviewerProjection
==========

The :mod:`swiftsimio.visualisation.projection` sub-module provides an interface
to render SWIFT data projected to a grid. This takes your 3D data and projects
it down to 2D, such that if you request masses to be smoothed then these
functions return a surface density.

This effectively solves the equation:

:math:`\tilde{A}_i = \sum_j A_j W_{ij, 2D}`

with :math:`\tilde{A}_i` the smoothed quantity in pixel :math:`i`, and
:math:`j` all particles in the simulation, with :math:`W` the 2D kernel.
Here we use the Wendland-C2 kernel.

The primary function here is
:meth:`swiftsimio.visualisation.projection.project_gas`, which allows you to
create a gas projection of any field. See the example below.

Example
-------

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.projection import project_gas

   data = load("cosmo_volume_example.hdf5")

   # This creates a grid that has units msun / Mpc^2, and can be transformed like
   # any other unyt quantity
   mass_map = project_gas(data, resolution=1024, project="masses", parallel=True)

   # Let's say we wish to save it as msun / kpc^2,
   from unyt import msun, kpc
   mass_map.convert_to_units(msun / kpc**2)

   from matplotlib.pyplot import imsave
   from matplotlib.colors import LogNorm

   # Normalize and save
   imsave("gas_surface_dens_map.png", LogNorm()(mass_map.value), cmap="viridis")


This basic demonstration creates a mass surface density map.

To create, for example, a projected temperature map, we need to remove the
surface density dependence (i.e. :meth:`project_gas` returns a surface
temperature in units of K / kpc^2 and we just want K) by dividing out by
this:

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.projection import project_gas

   data = load("cosmo_volume_example.hdf5")

   # First create a mass-weighted temperature dataset
   data.gas.mass_weighted_temps = data.gas.masses * data.gas.temperatures

   # Map in msun / mpc^2
   mass_map = project_gas(data, resolution=1024, project="masses", parallel=True)
   # Map in msun * K / mpc^2
   mass_weighted_temp_map = project_gas(
       data,
       resolution=1024,
       project="mass_weighted_temps",
       parallel=True
   )

   temp_map = mass_weighted_temp_map / mass_map

   from unyt import K
   temp_map.convert_to_units(K)

   from matplotlib.pyplot import imsave
   from matplotlib.colors import LogNorm

   # Normalize and save
   imsave("temp_map.png", LogNorm()(temp_map.value), cmap="twilight")


The output from this example, when used with the example data provided in the
loading data section should look something like:

.. image:: temp_map.png

Backends
--------

In certain cases, rather than just using this facility for visualisation, you
will wish that the values that are returned to be as well converged as
possible. For this, we provide several different backends. These are passed
as ``backend="str"`` to all of the projection visualisation functions, and
are available in the module
:mod:`swiftsimio.visualisation.projection.projection_backends`. The available
backends are as follows:

+ ``fast``: The default backend - this is extremely fast, and provides very basic
  smoothing, with a return type of single precision floating point numbers.
+ ``histogram``: This backend provides zero smoothing, and acts in a similar way
  to the ``np.hist2d`` function but with the same arguments as ``scatter``.
+ ``reference``: The same backend as ``fast`` but with two distinguishing features;
  all calculations are performed in double precision, and it will return early
  with a warning message if there are not enough pixels to fully resolve each kernel.
  Regular users should not use this mode.
+ ``renormalised``: The same as ``fast``, but each kernel is evaluated twice and
  renormalised to ensure mass conservation within floating point precision. Returns
  single precision arrays.
+ ``subsampled``: This is the recommended mode for users who wish to have converged
  results even at low resolution. Each kernel is evaluated at least 32 times, with
  overlaps between pixels considered for every single particle. Returns in
  double precision.
+ ``subsampled_extreme``: The same as ``subsampled``, but provides 64 kernel
  evaluations.
+ ``gpu``: The same as ``fast`` but uses CUDA for faster computation on supported
  GPUs. The parallel implementation is the same function as the non-parallel.

Example:

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.projection import project_gas

   data = load("cosmo_volume_example.hdf5")

   subsampled_array = project_gas(
      data,
      resolution=1024,
      project="entropies",
      parallel=True,
      backend="subsampled"
   )

This will likely look very similar to the image that you make with the default
``backend="fast"``, but will have a well-converged distribution at any resolution
level.

Rotations
---------

Sometimes you will need to visualise a galaxy from a different perspective.
The :mod:`swiftsimio.visualisation.rotation` sub-module provides routines to
generate rotation matrices corresponding to vectors, which can then be
provided to the ``rotation_matrix`` argument of :meth:`project_gas` (and
:meth:`project_gas_pixel_grid`). You will also need to supply the
``rotation_center`` argument, as the rotation takes place around this given
point. The example code below loads a snapshot, and a halo catalogue, and
creates an edge-on and face-on projection using the integration in
``velociraptor``. More information on possible integrations with this library
is shown in the ``velociraptor`` section.

.. code-block:: python

   from swiftsimio import load, mask
   from velociraptor import load as load_catalogue
   from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
   from swiftsimio.visualisation.projection import project_gas_pixel_grid

   import unyt
   import numpy as np
   import matplotlib.pyplot as plt
   from matplotlib.colors import LogNorm

   # Radius around which to load data, we will visualise half of this
   size = 1000 * unyt.kpc

   snapshot_filename = "cosmo_volume_example.hdf5"
   catalogue_filename = "cosmo_volume_example.properties"

   catalogue = load_catalogue(catalogue_filename)

   # Which halo should we visualise?
   halo = 0

   x = catalogue.positions.xcmbp[halo]
   y = catalogue.positions.ycmbp[halo]
   z = catalogue.positions.zcmbp[halo]

   lx = catalogue.angular_momentum.lx[halo]
   ly = catalogue.angular_momentum.ly[halo]
   lz = catalogue.angular_momentum.lz[halo]

   # The angular momentum vector will point perpendicular to the galaxy disk.
   # If your simulation contains stars, use lx_star
   angular_momentum_vector = np.array([lx.value, ly.value, lz.value])
   angular_momentum_vector /= np.linalg.norm(angular_momentum_vector)

   face_on_rotation_matrix = rotation_matrix_from_vector(
      angular_momentum_vector
   )
   edge_on_rotation_matrix = rotation_matrix_from_vector(
      angular_momentum_vector,
      axis="y"
   )

   region = [
      [x - size, x + size],
      [y - size, y + size],
      [z - size, z + size],
   ]

   visualise_region = [
      x - 0.5 * size, x + 0.5 * size,
      y - 0.5 * size, y + 0.5 * size,
   ]

   data_mask = mask(snapshot_filename)
   data_mask.constrain_spatial(region)
   data = load(snapshot_filename, mask=data_mask)

   # Use project_gas_pixel_grid to generate projected images

   common_arguments = dict(data=data, resolution=512, parallel=True, region=visualise_region)

   un_rotated = project_gas_pixel_grid(**common_arguments)

   face_on = project_gas_pixel_grid(
      **common_arguments,
      rotation_center=unyt.unyt_array([x, y, z]),
      rotation_matrix=face_on_rotation_matrix,
   )

   edge_on = project_gas_pixel_grid(
      **common_arguments,
      rotation_center=unyt.unyt_array([x, y, z]),
      rotation_matrix=edge_on_rotation_matrix,
   )

Using this with the provided example data will just show blobs due to its low resolution
nature. Using one of the EAGLE volumes (``examples/EAGLE_ICs``) will produce much nicer
galaxies, but that data is too large to provide as an example in this tutorial.

You can also provide an extra two values, the z min and max, as part of the
``region`` parameter. This may have some slight performance impact, so it is
generally advised that you do this on sub-loaded volumes only.


Other particle types
--------------------

Other particle types are able to be visualised through the use of the
:meth:`swiftsimio.visualisation.projection.project_pixel_grid` function. This
does not attach correct symbolic units, so you will have to work those out
yourself, but it does perform the smoothing. We aim to introduce the feature
of correctly applied units to these projections soon.

To use this feature for particle types that do not have smoothing lengths, you
will need to generate them, as in the example below where we create a
mass density map for dark matter. We provide a utility to do this through
:meth:`swiftsimio.visualisation.smoothing_length_generation.generate_smoothing_lengths`.

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.projection import project_pixel_grid
   from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths

   data = load("cosmo_volume_example.hdf5")

   # Generate smoothing lengths for the dark matter
   data.dark_matter.smoothing_length = generate_smoothing_lengths(
       data.dark_matter.coordinates,
       data.metadata.boxsize,
       kernel_gamma=1.8,
       neighbours=57,
       speedup_fac=2,
       dimension=3,
   )

   # Project the dark matter mass
   dm_mass = project_pixel_grid(
       # Note here that we pass in the dark matter dataset not the whole
       # data object, to specify what particle type we wish to visualise
       data=data.dark_matter,
       boxsize=data.metadata.boxsize,
       resolution=1024,
       project="masses",
       parallel=True,
       region=None
   )

   from matplotlib.pyplot import imsave
   from matplotlib.colors import LogNorm

   # Everyone knows that dark matter is purple
   imsave("dm_mass_map.png", LogNorm()(dm_mass), cmap="inferno")

The output from this example, when used with the example data provided in the
loading data section should look something like:

.. image:: dm_mass_map.png


Lower-level API
---------------

The lower-level API for projections allows for any general positions,
smoothing lengths, and smoothed quantities, to generate a pixel grid that
represents the smoothed version of the data.

This API is available through
:meth:`swiftsimio.visualisation.projection.scatter` and
:meth:`swiftsimio.visualisation.projection.scatter_parallel` for the parallel
version. The parallel version uses significantly more memory as it allocates
a thread-local image array for each thread, summing them in the end. Here we
will only describe the ``scatter`` variant, but they behave in the exact same way.

By default this uses the "fast" backend. To use the others, you can select them
manually from the module, or by using the ``backends`` and ``backends_parallel``
dictionaries in :mod:`swiftsimio.visualisation.projection`.

To use this function, you will need:

+ x-positions of all of your particles, ``x``.
+ y-positions of all of your particles, ``y``.
+ A quantity which you wish to smooth for all particles, such as their
  mass, ``m``.
+ Smoothing lengths for all particles, ``h``.
+ The resolution you wish to make your square image at, ``res``.

The key here is that only particles in the domain [0, 1] in x, and [0, 1] in y
will be visible in the image. You may have particles outside of this range;
they will not crash the code, and may even contribute to the image if their
smoothing lengths overlap with [0, 1]. You will need to re-scale your data
such that it lives within this range. Then you may use the function as follows:

.. code-block:: python

   from swiftsimio.visualisation.projection import scatter

   # Using the variable names from above
   out = scatter(x=x, y=y, h=h, m=m, res=res)

``out`` will be a 2D :mod:`numpy` grid of shape ``[res, res]``. You will need
to re-scale this back to your original dimensions to get it in the correct units,
and do not forget that it now represents the smoothed quantity per surface area.
Slices
======

The :mod:`swiftsimio.visualisation.slice` sub-module provides an interface
to render SWIFT data onto a slice. This takes your 3D data and finds the 3D
density at fixed z-position, slicing through the box.

This effectively solves the equation:

:math:`\tilde{A}_i = \sum_j A_j W_{ij, 3D}`

with :math:`\tilde{A}_i` the smoothed quantity in pixel :math:`i`, and
:math:`j` all particles in the simulation, with :math:`W` the 3D kernel.
Here we use the Wendland-C2 kernel. Note that here we take the kernel
at a fixed z-position.

The primary function here is
:meth:`swiftsimio.visualisation.slice.slice_gas`, which allows you to
create a gas slice of any field. See the example below.

Example
-------

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.slice import slice_gas

   data = load("cosmo_volume_example.hdf5")

   # This creates a grid that has units msun / Mpc^3, and can be transformed like
   # any other unyt quantity. Note that `slice` is given in terms of the box-size,
   # so here we are taking a slice at z = boxsize / 2.
   mass_map = slice_gas(
       data,
       slice=0.5,
       resolution=1024,
       project="masses",
       parallel=True
   )

   # Let's say we wish to save it as g / cm^2,
   from unyt import g, cm
   mass_map.convert_to_units(g / cm**3)

   from matplotlib.pyplot import imsave
   from matplotlib.colors import LogNorm

   # Normalize and save
   imsave("gas_slice_map.png", LogNorm()(mass_map.value), cmap="viridis")


This basic demonstration creates a mass density map.

To create, for example, a projected temperature map, we need to remove the
density dependence (i.e. :meth:`slice_gas` returns a volumetric temperature
in units of K / kpc^3 and we just want K) by dividing out by this:

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.slice import slice_gas

   data = load("cosmo_volume_example.hdf5")

   # First create a mass-weighted temperature dataset
   data.gas.mass_weighted_temps = data.gas.masses * data.gas.temperatures

   # Map in msun / mpc^3
   mass_map = slice_gas(
       data,
       slice=0.5,
       resolution=1024,
       project="masses",
       parallel=True
   )

   # Map in msun * K / mpc^3
   mass_weighted_temp_map = slice_gas(
       data,
       slice=0.5,
       resolution=1024,
       project="mass_weighted_temps",
       parallel=True
   )

   temp_map = mass_weighted_temp_map / mass_map

   from unyt import K
   temp_map.convert_to_units(K)

   from matplotlib.pyplot import imsave
   from matplotlib.colors import LogNorm

   # Normalize and save
   imsave("temp_map.png", LogNorm()(temp_map.value), cmap="twilight")


The output from this example, when used with the example data provided in the
loading data section should look something like:

.. image:: temp_slice.png

Rotations
---------

Rotations of the box prior to slicing are provided in a similar fashion to the 
:mod:`swiftsimio.visualisation.projection` sub-module, by using the 
:mod:`swiftsimio.visualisation.rotation` sub-module. To rotate the perspective
prior to slicing a ``rotation_center`` argument in :meth:`slice_gas` needs
to be provided, specifying the point around which the rotation takes place. 
The angle of rotation is specified with a matrix, supplied by ``rotation_matrix``
in :meth:`slice_gas`. The rotation matrix may be computed with 
:meth:`rotation_matrix_from_vector`. This will result in the perspective being 
rotated to be along the provided vector. This approach to rotations applied to 
the above example is shown below.

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.slice import slice_gas
   from swiftsimio.visualisation.rotation import rotation_matrix_from_vector

   data = load("cosmo_volume_example.hdf5")

   # First create a mass-weighted temperature dataset
   data.gas.mass_weighted_temps = data.gas.masses * data.gas.temperatures

   # Specify the rotation parameters
   center = 0.5 * data.metadata.boxsize
   rotate_vec = [0.5,0.5,1]
   matrix = rotation_matrix_from_vector(rotate_vec, axis='z')
   
   # Map in msun / mpc^3
   mass_map = slice_gas(
       data,
       slice=0.5,
       resolution=1024,
       project="masses",
       rotation_matrix=matrix,
       rotation_center=center,
       parallel=True
   )
   
   # Map in msun * K / mpc^3
   mass_weighted_temp_map = slice_gas(
       data, 
       slice=0.5,
       resolution=1024,
       project="mass_weighted_temps",
       rotation_matrix=matrix,
       rotation_center=center,
       parallel=True
   )

   temp_map = mass_weighted_temp_map / mass_map

   from unyt import K
   temp_map.convert_to_units(K)

   from matplotlib.pyplot import imsave
   from matplotlib.colors import LogNorm

   # Normalize and save
   imsave("temp_map.png", LogNorm()(temp_map.value), cmap="twilight")


Lower-level API
---------------

The lower-level API for slices allows for any general positions,
smoothing lengths, and smoothed quantities, to generate a pixel grid that
represents the smoothed, sliced, version of the data.

This API is available through
:meth:`swiftsimio.visualisation.slice.slice_scatter` and
:meth:`swiftsimio.visualisation.slice.slice_scatter_parallel` for the parallel
version. The parallel version uses significantly more memory as it allocates
a thread-local image array for each thread, summing them in the end. Here we
will only describe the ``scatter`` variant, but they behave in the exact same way.

To use this function, you will need:

+ x-positions of all of your particles, ``x``.
+ y-positions of all of your particles, ``y``.
+ z-positions of all of your particles, ``z``.
+ Where in the [0,1] range you wish to slice, ``z_slice``.
+ A quantity which you wish to smooth for all particles, such as their
  mass, ``m``.
+ Smoothing lengths for all particles, ``h``.
+ The resolution you wish to make your square image at, ``res``.

The key here is that only particles in the domain [0, 1] in x, [0, 1] in y,
and [0, 1] in z. will be visible in the image. You may have particles outside
of this range; they will not crash the code, and may even contribute to the
image if their smoothing lengths overlap with [0, 1]. You will need to
re-scale your data such that it lives within this range. Then you may use the
function as follows:

.. code-block:: python

   from swiftsimio.visualisation.slice import slice_scatter

   # Using the variable names from above
   out = slice_scatter(x=x, y=y, z=z, h=h, m=m, z_slice=z_slice, res=res)

``out`` will be a 2D :mod:`numpy` grid of shape ``[res, res]``. You will need
to re-scale this back to your original dimensions to get it in the correct units,
and do not forget that it now represents the smoothed quantity per volume.
Volume Rendering
================

The :mod:`swiftsimio.visualisation.volume_render` sub-module provides an
interface to render SWIFT data onto a fixed grid. This takes your 3D data and
finds the 3D density at fixed positions, allowing it to be used in codes that
require fixed grids such as radiative transfer programs.

This effectively solves the equation:

:math:`\tilde{A}_i = \sum_j A_j W_{ij, 3D}`

with :math:`\tilde{A}_i` the smoothed quantity in pixel :math:`i`, and
:math:`j` all particles in the simulation, with :math:`W` the 3D kernel.
Here we use the Wendland-C2 kernel.

The primary function here is
:meth:`swiftsimio.visualisation.volume_render.render_gas`, which allows you
to create a gas density grid of any field, see the example below.

Example
-------

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.volume_render import render_gas

   data = load("cosmo_volume_example.hdf5")

   # This creates a grid that has units msun / Mpc^3, and can be transformed like
   # any other unyt quantity.
   mass_grid = render_gas(
       data,
       resolution=256,
       project="masses",
       parallel=True
   )

This basic demonstration creates a mass density cube.

To create, for example, a projected temperature cube, we need to remove the
density dependence (i.e. :meth:`render_gas` returns a volumetric
temperature in units of K / kpc^3 and we just want K) by dividing out by
this:

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.volume_render import render_gas

   data = load("cosmo_volume_example.hdf5")

   # First create a mass-weighted temperature dataset
   data.gas.mass_weighted_temps = data.gas.masses * data.gas.temperatures

   # Map in msun / mpc^3
   mass_cube = render_gas(
       data,
       resolution=256,
       project="masses",
       parallel=True
   )

   # Map in msun * K / mpc^3
   mass_weighted_temp_cube = render_gas(
       data,
       resolution=256,
       project="mass_weighted_temps",
       parallel=True
   )

   # A 256 x 256 x 256 cube with dimensions of temperature
   temp_cube = mass_weighted_temp_cube / mass_cube

Rotations
---------

Rotations of the box prior to volume rendering are provided in a similar fashion 
to the :mod:`swiftsimio.visualisation.projection` sub-module, by using the 
:mod:`swiftsimio.visualisation.rotation` sub-module. To rotate the perspective
prior to slicing a ``rotation_center`` argument in :meth:`render_gas` needs
to be provided, specifying the point around which the rotation takes place. 
The angle of rotation is specified with a matrix, supplied by ``rotation_matrix``
in :meth:`render_gas`. The rotation matrix may be computed with 
:meth:`rotation_matrix_from_vector`. This will result in the perspective being 
rotated to be along the provided vector. This approach to rotations applied to 
the above example is shown below.

.. code-block:: python

   from swiftsimio import load
   from swiftsimio.visualisation.volume_render import render_gas
   from swiftsimio.visualisation.rotation import rotation_matrix_from_vector

   data = load("cosmo_volume_example.hdf5")

   # First create a mass-weighted temperature dataset
   data.gas.mass_weighted_temps = data.gas.masses * data.gas.temperatures

   # Specify the rotation parameters
   center = 0.5 * data.metadata.boxsize
   rotate_vec = [0.5,0.5,1]
   matrix = rotation_matrix_from_vector(rotate_vec, axis='z')
   
   # Map in msun / mpc^3
   mass_cube = render_gas(
       data,
       resolution=256,
       project="masses",
       rotation_matrix=matrix,
       rotation_center=center,
       parallel=True
   )
   
   # Map in msun * K / mpc^3
   mass_weighted_temp_cube = render_gas(
       data, 
       resolution=256,
       project="mass_weighted_temps",
       rotation_matrix=matrix,
       rotation_center=center,
       parallel=True
   )

   # A 256 x 256 x 256 cube with dimensions of temperature
   temp_cube = mass_weighted_temp_cube / mass_cube

Lower-level API
---------------

The lower-level API for volume rendering allows for any general positions,
smoothing lengths, and smoothed quantities, to generate a pixel grid that
represents the smoothed, volume rendered, version of the data.

This API is available through
:meth:`swiftsimio.visualisation.volume_render.scatter` and
:meth:`swiftsimio.visualisation.volume_render.scatter_parallel` for the parallel
version. The parallel version uses significantly more memory as it allocates
a thread-local image array for each thread, summing them in the end. Here we
will only describe the ``scatter`` variant, but they behave in the exact same way.

To use this function, you will need:

+ x-positions of all of your particles, ``x``.
+ y-positions of all of your particles, ``y``.
+ z-positions of all of your particles, ``z``.
+ A quantity which you wish to smooth for all particles, such as their
  mass, ``m``.
+ Smoothing lengths for all particles, ``h``.
+ The resolution you wish to make your cube at, ``res``.

The key here is that only particles in the domain [0, 1] in x, [0, 1] in y,
and [0, 1] in z. will be visible in the cube. You may have particles outside
of this range; they will not crash the code, and may even contribute to the
image if their smoothing lengths overlap with [0, 1]. You will need to
re-scale your data such that it lives within this range. Then you may use the
function as follows:

.. code-block:: python

   from swiftsimio.visualisation.volume_render import scatter

   # Using the variable names from above
   out = scatter(x=x, y=y, z=z, h=h, m=m, res=res)

``out`` will be a 3D :mod:`numpy` grid of shape ``[res, res, res]``. You will
need to re-scale this back to your original dimensions to get it in the
correct units, and do not forget that it now represents the smoothed quantity
per volume.
Visualisation
=============

:mod:`swiftsimio` provides visualisation routines accelerated with the
:mod:`numba` module. They work without this module, but we strongly recommend
installing it for the best performance (1000x+ speedups). These are provided
in the :mod:`swiftismio.visualisation` sub-modules.

The three built-in rendering types (described below) have the following
common interface:

.. code-block:: python

   {render_func_name}_gas(
       data=data, # SWIFTsimIO dataset
       resolution=1024, # Resolution along one axis of the output image
       project="masses", # Variable to project, e.g. masses, temperatures, etc.
       parallel=False, # Construct the image in (thread) parallel?
       region=None, # None, or a list telling which region to render_func_name
   )

The output of these functions comes with associated units and has the correct
dimensions. There are lower-level APIs (also documented here) that provide
additional functionality.

Finally, we also describe here the integration with :mod:`py-sphviewer`.

.. toctree::
   :maxdepth: 2

   projection
   slice
   volume_render
   tools
   py-sphviewer

Masking
=======

:mod:`swiftsimio` provides unique functionality (when compared to other
software packages that read HDF5 data) through its masking facility.

SWIFT snapshots contain cell metadata that allow us to spatially mask the
data ahead of time. :mod:`swiftsimio` provides a number of objects that help
with this. This functionality is provided through the :mod:`swiftsimio.masks`
sub-module but is available easily through the :meth:`swiftsimio.mask`
top-level function.

This functionality is used heavily in our `VELOCIraptor integration library`_
for only reading data that is near bound objects.

There are two types of mask, with the default only allowing spatial masking.
Full masks require significantly more memory overhead and are generally much
slower than the spatial only mask.

.. _`VELOCIraptor integration library`: https://github.com/swiftsim/velociraptor-python

Spatial-only masking
--------------------

Spatial only masking is approximate and allows you to only load particles
within a given region. It is precise to the top-level cells that are defined
within SWIFT. It will always load all of the particles that you request, but
for simplicity it may also load some particles that are slightly outside
of the region of interest. This is because it works as follows:

1. Load the top-level cell metadata.
2. Find the overlap between the specified region and these cells.
3. Load all cells within that overlap.

As you can see, the edges of regions may load in extra information as we
always load the whole top-level cell.

Example
^^^^^^^

In this example we will use the :obj:`swiftsimio.masks.SWIFTMask` object
to load the bottom left 'half' corner of the box.

.. code-block:: python

   import swiftsimio as sw

   filename = "cosmo_volume_example.hdf5"

   mask = sw.mask(filename)
   # The full metadata object is available from within the mask
   boxsize = mask.metadata.boxsize
   # load_region is a 3x2 list [[left, right], [bottom, top], [front, back]]
   load_region = [[0.0 * b, 0.5 * b] for b in boxsize]

   # Constrain the mask
   mask.constrain_spatial(load_region)

   # Now load the snapshot with this mask
   data = sw.load(filename, mask=mask)

``data`` is now a regular :obj:`swiftsimio.reader.SWIFTDataset` object, but
it only ever loads particles that are (approximately) inside the
``load_region`` region.

Importantly, this method has a tiny memory overhead, and should also have a
relatively small overhead when reading the data. This allows you to use snapshots
that are much larger than the available memory on your machine and process them
with ease.

Full mask
---------

The below example shows the use of a full masking object, used to constrain
densities of particles and only load particles within that density window.

.. code-block:: python
   
   import swiftsimio as sw

   # This creates and sets up the masking object.
   mask = sw.mask("cosmological_volume.hdf5", spatial_only=False)

   # This ahead-of-time creates a spatial mask based on the cell metadata.
   mask.constrain_spatial([
       [0.2 * mask.metadata.boxsize[0], 0.7 * mask.metadata.boxsize[0]],
       None,
       None]
   )

   # Now, just for fun, we also constrain the density between
   # 0.4 g/cm^3 and 0.8. This reads in the relevant data in the region,
   # and tests it element-by-element. Note that using masks of this type
   # is significantly slower than using the spatial-only masking.
   density_units = mask.units.mass / mask.units.length**3
   mask.constrain_mask("gas", "density", 0.4 * density_units, 0.8 * density_units)

   # Now we can grab the actual data object. This includes the mask as a parameter.
   data = sw.load("cosmo_volume_example.hdf5", mask=mask)


When the attributes of this data object are accessed, *only* the ones that
belong to the masked region (in both density and spatial) are read. I.e. if I
ask for the temperature of particles, it will recieve an array containing
temperatures of particles that lie in the region [0.2, 0.7] and have a
density between 0.4 and 0.8 g/cm^3.

Writing subset of snapshot
--------------------------
In some cases it may be useful to write a subset of an existing snapshot to its
own hdf5 file. This could be used, for example, to extract a galaxy halo that 
is of interest from a snapshot so that the file is easier to work with and transport.

To do this the ``write_subset`` function is provided. It can be used, for example,
as follows

.. code-block:: python

    import swiftsimio as sw                                                 
    import unyt                                                             
    
    mask = sw.mask("eagle_snapshot.hdf5")                                       
    mask.constrain_spatial([
        [unyt.unyt_quantity(100, unyt.kpc), unyt.unyt_quantity(1000, unyt.kpc)], 
        None, 
        None])                                   
    
    sw.subset_writer.write_subset("test_subset.hdf5", mask)

This will write a snapshot which contains the particles from the specified snapshot 
whose *x*-coordinate is within the range [100, 1000] kpc. This function uses the 
cell mask which encompases the specified spatial domain to successively read portions 
of datasets from the input file and writes them to a new snapshot. 

Due to the coarse grained nature of the cell mask, particles from outside this range 
may also be included if they are within the same top level cells as particles that 
fall within the given range.

Please note that it is important to run ``constrain_spatial`` as this generates
and stores the cell mask needed to write the snapshot subset. 
Statistics Files
================

:mod:`swiftsimio` includes routines to load log files, such as the
``SFR.txt`` and ``energy.txt``. This is available through the
:obj:`swiftsimio.statistics.SWIFTStatisticsFile` object, or through
the main ``load_statistics`` function.

Example
-------

.. code-block:: python

   from swiftsimio import load_statistics

   data = load_statistics("energy.txt")

   print(data)

   print(x.total_mass.name)


Will output:

.. code-block:: bash

   Statistics file: energy.txt, containing fields: #, step, time, a, z, total_mass,
   gas_mass, dm_mass, sink_mass, star_mass, bh_mass, gas_z_mass, star_z_mass,
   bh_z_mass, kin_energy, int_energy, pot_energy, rad_energy, gas_entropy, com_x,
   com_y, com_z, mom_x, mom_y, mom_z, ang_mom_x, ang_mom_y, ang_mom_z

   'Total mass in the simulation'Creating Initial Conditions
===========================

Writing datasets that are valid for consumption for cosmological codes can be
difficult, especially when considering how to best use units. SWIFT uses a
different set of internal units (specified in your parameter file) that does
not necessarily need to be the same set of units that initial conditions are
specified in. Nevertheless, it is important to ensure that units in the
initial conditions are all *consistent* with each other. To facilitate this,
we use :mod:`unyt` arrays. The below example generates randomly placed gas
particles with uniform densities.

The functionality to create initial conditions is available through
the :mod:`swiftsimio.writer` sub-module, and the top-level
:obj:`swiftsimio.Writer` object.

Note that the properties that :mod:`swiftsimio` requires in the initial
conditions are the only ones that are actually read by SWIFT; other fields
will be left un-read and as such should not be included in initial conditions
files.

A current known issue is that due to inconsistencies with the initial
conditions and simulation snapshots, :mod:`swiftsimio` is not actually able
to read the inititial conditions that it produces. We are aiming to fix this
in an upcoming release.


Example
^^^^^^^

.. code-block:: python

   from swiftsimio import Writer
   from swiftsimio.units import cosmo_units

   import unyt
   import numpy as np

   # Box is 100 Mpc
   boxsize = 100 * unyt.Mpc

   # Generate object. cosmo_units corresponds to default Gadget-oid units
   # of 10^10 Msun, Mpc, and km/s
   x = Writer(cosmo_units, boxsize)

   # 32^3 particles.
   n_p = 32**3

   # Randomly spaced coordinates from 0, 100 Mpc in each direction
   x.gas.coordinates = np.random.rand(n_p, 3) * (100 * unyt.Mpc)

   # Random velocities from 0 to 1 km/s
   x.gas.velocities = np.random.rand(n_p, 3) * (unyt.km / unyt.s)

   # Generate uniform masses as 10^6 solar masses for each particle
   x.gas.masses = np.ones(n_p, dtype=float) * (1e6 * unyt.msun)

   # Generate internal energy corresponding to 10^4 K
   x.gas.internal_energy = (
       np.ones(n_p, dtype=float) * (1e4 * unyt.kb * unyt.K) / (1e6 * unyt.msun)
   )

   # Generate initial guess for smoothing lengths based on MIPS
   x.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

   # If IDs are not present, this automatically generates
   x.write("test.hdf5")

Then, running ``h5glance`` on the resulting ``test.hdf5`` produces:

.. code-block:: bash

   test.hdf5
   ├Header
   │ └5 attributes:
   │   ├BoxSize: 100.0
   │   ├Dimension: array [int64: 1]
   │   ├Flag_Entropy_ICs: 0
   │   ├NumPart_Total: array [int64: 6]
   │   └NumPart_Total_HighWord: array [int64: 6]
   ├PartType0
   │ ├Coordinates  [float64: 32768 × 3]
   │ ├InternalEnergy       [float64: 32768]
   │ ├Masses       [float64: 32768]
   │ ├ParticleIDs  [float64: 32768]
   │ ├SmoothingLength      [float64: 32768]
   │ └Velocities   [float64: 32768 × 3]
   └Units
   └5 attributes:
       ├Unit current in cgs (U_I): array [float64: 1]
       ├Unit length in cgs (U_L): array [float64: 1]
       ├Unit mass in cgs (U_M): array [float64: 1]
       ├Unit temperature in cgs (U_T): array [float64: 1]
       └Unit time in cgs (U_t): array [float64: 1]

**Note** you do need to be careful that your choice of unit system does
*not* allow values over 2^31, i.e. you need to ensure that your
provided values (with units) when *written* to the file are safe to 
be interpreted as (single-precision) floats. The only exception to
this is coordinates which are stored in double precision.
Getting Started
===============

The SWIFT astrophysical simulation code (http://swift.dur.ac.uk) is used
widely. There exists many ways of reading the data from SWIFT, which outputs
HDF5 files. These range from reading directly using :mod:`h5py` to using a
complex system such as :mod:`yt`; however these either are unsatisfactory
(e.g. a lack of unit information in reading HDF5), or too complex for most
use-cases. This (thin) wrapper provides an object-oriented API to read
(dynamically) data from SWIFT.

Getting set up with :mod:`swiftsimio` is easy; it (by design) has very few
requirements. There are a number of optional packages that you can install
to make the experience better and these are recommended. All requirements
are detailed below.


Requirements
------------

This requires ``python`` ``v3.6.0`` or higher. Unfortunately it is not
possible to support :mod:`swiftsimio` on versions of python lower than this.
It is important that you upgrade if you are still a ``python2`` user.

Python packages
^^^^^^^^^^^^^^^

+ ``numpy``, required for the core numerical routines.
+ ``h5py``, required to read data from the SWIFT HDF5 output files.
+ ``unyt``, required for symbolic unit calculations (depends on ``sympy``).

Optional packages
^^^^^^^^^^^^^^^^^

+ ``numba``, highly recommended should you wish to use the in-built visualisation
  tools.
+ ``scipy``, required if you wish to generate smoothing lengths for particle types
  that do not store this variable in the snapshots (e.g. dark matter)
+ ``tqdm``, required for progress bars for some long-running tasks. If not installed
  no progress bar will be shown.
+ ``py-sphviewer``, if you wish to use our integration with this visualisation
  code.


Installing
----------

:mod:`swiftsimio` can be installed using the python packaging manager, ``pip``,
or any other packaging manager that you wish to use:

``pip install swiftsimio``

Note that this will install any required packages for you.

To set up the code for development, first clone the latest master from GitHub:

``git clone https://github.com/SWIFTSIM/swiftsimio.git``

and install with ``pip`` using the ``-e`` flag,

``cd swiftsimio``

``pip install -e .``

.. TODO: Add contribution guide.

Usage
-----

There are many examples of using :mod:`swiftsimio` available in the
swiftsimio_examples_ repository, which also includes examples for reading
older (e.g. EAGLE) datasets.

Example usage is shown below, which plots a density-temperature phase
diagram, with density and temperature given in CGS units:

.. code-block:: python

   import swiftsimio as sw

   # This loads all metadata but explicitly does _not_ read any particle data yet
   data = sw.load("/path/to/swift/output")

   import matplotlib.pyplot as plt

   data.gas.densities.convert_to_cgs()
   data.gas.temperatures.convert_to_cgs()

   plt.loglog()

   plt.scatter(
      data.gas.densities,
      data.gas.temperatures,
      s=1
   )

   plt.xlabel(fr"Gas density $\left[{data.gas.densities.units.latex_repr}\right]$")
   plt.ylabel(fr"Gas temperature $\left[{data.gas.temperatures.units.latex_repr}\right]$")

   plt.tight_layout()

   plt.savefig("test_plot.png", dpi=300)


Don't worry too much about this for now if you can't understand it, we will
get into this much more heavily in the next section.

In the above it's important to note the following:

+ All metadata is read in when the :meth:`swiftsimio.load` function is called.
+ Only the density and temperatures (corresponding to the ``PartType0/Densities`` and
  ``PartType0/Temperatures``) datasets are read in.
+ That data is only read in once the
  :meth:`swiftsimio.objects.cosmo_array.convert_to_cgs` method is called.
+ :meth:`swiftsimio.objects.cosmo_array.convert_to_cgs` converts data in-place;
  i.e. it returns `None`.
+ The data is cached and not re-read in when ``plt.scatter`` is called.


.. _swiftsimio_examples: https://github.com/swiftsim/swiftsimio-examples
.. API Documentation

API Documentation
=================

.. toctree::
   :maxdepth: 3

   swiftsimio




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
