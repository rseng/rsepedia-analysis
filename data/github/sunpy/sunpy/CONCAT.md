Acknowledging or Citing SunPy
=============================

If you use SunPy in your scientific work, we would appreciate citing it in your publications.
The continued growth and development of SunPy is dependent on the community being aware of SunPy.

Citing SunPy in Publications
----------------------------

Please add the following line within your methods, conclusion or acknowledgements sections:

   *This research used version X.Y.Z (software citation) of the SunPy open source
   software package (project citation).*

The project citation should be to the `SunPy paper`_, and the software citation should be the specific `Zenodo DOI`_ for the version used in your work.

.. code:: bibtex

    @ARTICLE{sunpy_community2020,
      doi = {10.3847/1538-4357/ab4f7a},
      url = {https://iopscience.iop.org/article/10.3847/1538-4357/ab4f7a},
      author = {{The SunPy Community} and Barnes, Will T. and Bobra, Monica G. and Christe, Steven D. and Freij, Nabil and Hayes, Laura A. and Ireland, Jack and Mumford, Stuart and Perez-Suarez, David and Ryan, Daniel F. and Shih, Albert Y. and Chanda, Prateek and Glogowski, Kolja and Hewett, Russell and Hughitt, V. Keith and Hill, Andrew and Hiware, Kaustubh and Inglis, Andrew and Kirk, Michael S. F. and Konge, Sudarshan and Mason, James Paul and Maloney, Shane Anthony and Murray, Sophie A. and Panda, Asish and Park, Jongyeob and Pereira, Tiago M. D. and Reardon, Kevin and Savage, Sabrina and Sipőcz, Brigitta M. and Stansby, David and Jain, Yash and Taylor, Garrison and Yadav, Tannmay and Rajul and Dang, Trung Kien},
      title = {The SunPy Project: Open Source Development and Status of the Version 1.0 Core Package},
      journal = {The Astrophysical Journal},
      volume = {890},
      issue = {1},
      pages = {68-},
      publisher = {American Astronomical Society},
      year = {2020}
    }

You can also get this information with ``sunpy.__citation__``.

Other SunPy publications
########################

The SunPy v1.0.8 code was reviewed by `The Journal of Open Source Software (JOSS) <https://joss.theoj.org/papers/10.21105/joss.01832>`__.

.. code:: bibtex

    @ARTICLE{Mumford2020,
      doi = {10.21105/joss.01832},
      url = {https://doi.org/10.21105/joss.01832},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {46},
      pages = {1832},
      author = {Stuart Mumford and Nabil Freij and Steven Christe and Jack Ireland and Florian Mayer and V. Hughitt and Albert Shih and Daniel Ryan and Simon Liedtke and David Pérez-Suárez and Pritish Chakraborty and Vishnunarayan K and Andrew Inglis and Punyaslok Pattnaik and Brigitta Sipőcz and Rishabh Sharma and Andrew Leonard and David Stansby and Russell Hewett and Alex Hamilton and Laura Hayes and Asish Panda and Matt Earnshaw and Nitin Choudhary and Ankit Kumar and Prateek Chanda and Md Haque and Michael Kirk and Michael Mueller and Sudarshan Konge and Rajul Srivastava and Yash Jain and Samuel Bennett and Ankit Baruah and Will Barnes and Michael Charlton and Shane Maloney and Nicky Chorley and Himanshu  and Sanskar Modi and James Mason and Naman9639  and Jose Rozo and Larry Manley and Agneet Chatterjee and John Evans and Michael Malocha and Monica Bobra and Sourav Ghosh and Airmansmith97  and Dominik Stańczak and Ruben De Visscher and Shresth Verma and Ankit Agrawal and Dumindu Buddhika and Swapnil Sharma and Jongyeob Park and Matt Bates and Dhruv Goel and Garrison Taylor and Goran Cetusic and Jacob  and Mateo Inchaurrandieta and Sally Dacie and Sanjeev Dubey and Deepankar Sharma and Erik Bray and Jai Rideout and Serge Zahniy and Tomas Meszaros and Abhigyan Bose and André Chicrala and Ankit  and Chloé Guennou and Daniel D'Avella and Daniel Williams and Jordan Ballew and Nick Murphy and Priyank Lodha and Thomas Robitaille and Yash Krishan and Andrew Hill and Arthur Eigenbrot and Benjamin Mampaey and Bernhard Wiedemann and Carlos Molina and Duygu Keşkek and Ishtyaq Habib and Joseph Letts and Juanjo Bazán and Quinn Arbolante and Reid Gomillion and Yash Kothari and Yash Sharma and Abigail Stevens and Adrian Price-Whelan and Ambar Mehrotra and Arseniy Kustov and Brandon Stone and Trung Dang and Emmanuel Arias and Fionnlagh Dover and Freek Verstringe and Gulshan Kumar and Harsh Mathur and Igor Babuschkin and Jaylen Wimbish and Juan Buitrago-Casas and Kalpesh Krishna and Kaustubh Hiware and Manas Mangaonkar and Matthew Mendero and Mickaël Schoentgen and Norbert Gyenge and Ole Streicher and Rajasekhar Mekala and Rishabh Mishra and Shashank Srikanth and Sarthak Jain and Tannmay Yadav and Tessa Wilkinson and Tiago Pereira and Yudhik Agrawal and Jamescalixto  and Yasintoda  and Sophie Murray},
      title = {SunPy: A Python package for Solar Physics},
      journal = {Journal of Open Source Software}
    }

A paper about version v0.5 of SunPy was published in `Computational Science & Discovery <https://iopscience.iop.org/article/10.1088/1749-4699/8/1/014009>`__.

Acknowledging SunPy in Posters and Talks
----------------------------------------

Please include the `Sunpy logo`_ on the title, conclusion slide, or about page.
For websites please link the image to `sunpy.org`_.
Other versions of the logo are available in the `sunpy-logo repository`_.

.. _SunPy paper: https://doi.org/10.3847/1538-4357/ab4f7a
.. _Sunpy logo: https://github.com/sunpy/sunpy-logo/blob/master/sunpy_logo.svg
.. _sunpy.org: https://sunpy.org/
.. _sunpy-logo repository: https://github.com/sunpy/sunpy-logo/
.. _Zenodo DOI: https://doi.org/10.5281/zenodo.591887
3.1.0 (2021-10-29)
==================

Breaking Changes
----------------

- :meth:`sunpy.timeseries.sources.NOAAIndicesTimeSeries.peek` accepts ``plot_type`` as an argument instead of ``type``. (`#5200 <https://github.com/sunpy/sunpy/pull/5200>`__)
- Fill values are now set to `numpy.nan` in ``sunpy.timeseries.sources.noaa`` file
  parsers. They were previously set to a fill value of ``-1``. (`#5363 <https://github.com/sunpy/sunpy/pull/5363>`__)
- `sunpy.map.GenericMap.date` now looks for more metadata than just DATE-OBS,
  using new FITS keywords defined in version 4 of the standard.
  `sunpy.map.GenericMap.date` now returns, in order of preference:

  1. The DATE-OBS FITS keyword
  2. `~sunpy.map.GenericMap.date_average`
  3. `~sunpy.map.GenericMap.date_start`
  4. `~sunpy.map.GenericMap.date_end`
  5. The current time.

  If DATE-OBS is present alongside DATE-AVG or DATE-BEG and DATE-END, this results
  in a behaviour change to favour the new (more precisely defined) keywords.
  It is recommended
  to use `~sunpy.map.GenericMap.date_average`,
  `~sunpy.map.GenericMap.date_start`, or `~sunpy.map.GenericMap.date_end`
  instead if you need one of these specific times. (`#5449 <https://github.com/sunpy/sunpy/pull/5449>`__)
- :func:`sunpy.io.fits.get_header` no longer automatically tries to add the
  WAVEUNIT keyword if it isn't present in the header. To replicate the original
  behaviour do::

    header = sunpy.io.fits.get_header(...)
    waveunit = sunpy.io.fits.extract_waveunit(header)
    if waveunit is not None:
        header['WAVEUNIT'] = waveunit

  The `sunpy.map.GenericMap.waveunit` property still uses
  :func:`sunpy.io.fits.extract_waveunit` to try and get the waveunit if the
  WAVEUNIT key isn't present. (`#5501 <https://github.com/sunpy/sunpy/pull/5501>`__)
- `sunpy.map.GenericMap.wcs` no longer passes the whole ``.meta`` dictionary to
  `astropy.wcs.WCS` when constructing ``.wcs``. Instead each metadata value is
  manually taken from various map properties, which allows fixes to be made to
  the WCS without modifying the original map header. We think that
  `~sunpy.map.GenericMap.wcs` correctly sets all the keys needed for a full WCS
  header, but if you find anything missing please open an issue on the sunpy
  issue tracker. (`#5501 <https://github.com/sunpy/sunpy/pull/5501>`__)


Deprecations
------------

- ``sunpy.util.scraper.Scraper`` has been moved into `sunpy.net`, please update your imports to be ``from sunpy.net import Scraper``. (`#5364 <https://github.com/sunpy/sunpy/pull/5364>`__)
- Using "neighbour" as a resampling method in
  :func:`sunpy.image.resample.resample` is deprecated. Use "nearest" instead,
  which has the same effect. (`#5480 <https://github.com/sunpy/sunpy/pull/5480>`__)
- The `sunpy.visualization.animator` subpackge has been spun out into the
  standalone `mpl-animators <https://pypi.org/project/mpl-animators>`_ package,
  with the exception of `~sunpy.visualization.animator.MapSequenceAnimator`.
  Please update your imports to replace ``sunpy.visualization.animator`` with
  ``mpl_animators``.

  This is primarily because the ``ndcube`` package now relies on the animator
  classes as well as `sunpy`. (`#5619 <https://github.com/sunpy/sunpy/pull/5619>`__)


Removals
--------

- The deprecated ``sunpy.roi.chaincode.Chaincode`` has been removed in favour of `sunpy.net.helio.Chaincode`. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.roi.roi`` was removed, there is no direct replacement but `astropy-regions <https://astropy-regions.readthedocs.io/en/latest/>`__ is something to consider. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.instr`` has been removed, please use `sunkit_instruments <https://docs.sunpy.org/projects/sunkit-instruments/en/stable/>`__. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.map.GenericMap.size`` has been removed, please use ``sunpy.map.GenericMap.data.size``. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ability to read txt files from `sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries` and `sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries` has been removed as the data provided by NOAA is now provided as JSON files. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- Removed various deprecated methods on our Fido clients and responses:

  1. ``UnifiedResponse.build_table``, ``UnifiedResponse.tables``, ``UnifiedResponse.responses``, ``UnifiedResponse.get_response`` and ``UnifiedResponse.blocks`` as ``UnifiedResponse`` is now an `astropy.table.Table` that is sliceable.
  2. ``UnifiedResponse.response_block_properties`` as ``UnifiedResponse.path_format_keys`` was added as a better replacement.
  3. ``HECClient.time_query`` as you can now use ``Fido.search`` directly.
  4. ``sunpy.net.jsoc.attrs.Keys`` was not used for querying JSOC.
  5. ``sunpy.net.jsoc.JSOCClient.search_metadata`` as the functionality this provided was merged into `sunpy.net.jsoc.JSOCClient.search`.
  6. ``sunpy.net.vso.VSOClient.link`` as better search support in the client replaces this method. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.map.GenericMap.draw_rectangle()`` has been removed, the replacement is :meth:`sunpy.map.GenericMap.draw_quadrangle` (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- sunpy now errors if the unused ``.rsun`` or ``.heliographic_observer``
  attributes are set on a `~astropy.wcs.WCS`. (`#5348 <https://github.com/sunpy/sunpy/pull/5348>`__)
- Support for passing non-unit levels to :meth:`sunpy.map.GenericMap.draw_contours`
  when map data has units set has been removed, and with now raise an error. (`#5352 <https://github.com/sunpy/sunpy/pull/5352>`__)
- The ``origin`` argument to :meth:`sunpy.map.GenericMap.world_to_pixel` and
  :meth:`sunpy.map.GenericMap.pixel_to_world` has been removed. (`#5353 <https://github.com/sunpy/sunpy/pull/5353>`__)
- Support for plotting or contouring `~sunpy.map.GenericMap` on axes that are not
  `~astropy.visualization.wcsaxes.WCSAxes` has been removed. To create a
  ``WCSAxes``, use the ``projection`` argument when the axes is created, e.g.
  ``fig.add_subplot(111, projection=my_map)``. (`#5354 <https://github.com/sunpy/sunpy/pull/5354>`__)
- The following search attributes in `sunpy.net.vso.attrs` have been removed:
  ``['Time', 'Instrument', 'Wavelength', 'Source', 'Provider',
  'Level', 'Sample', 'Detector', 'Resolution', 'Physobs']``.
  Use the equivalent attribute from `sunpy.net.attrs` instead. (`#5355 <https://github.com/sunpy/sunpy/pull/5355>`__)
- The default response format from the VSO client is now a table. (`#5355 <https://github.com/sunpy/sunpy/pull/5355>`__)
- ``sunpy.net.hek.attrs.Time`` has been removed, use `sunpy.net.attrs.Time` instead. (`#5355 <https://github.com/sunpy/sunpy/pull/5355>`__)


New Features
------------

- Ensured that ``plot`` and ``peek`` will output the same figures for all `sunpy.timeseries.TimeSeries` sources. (`#5200 <https://github.com/sunpy/sunpy/pull/5200>`__)
- Added hook file and tests for using PyInstaller with sunpy. (`#5224 <https://github.com/sunpy/sunpy/pull/5224>`__)
- Allows :meth:`sunpy.map.GenericMap.draw_quadrangle` to accept pixel units as input to enable plotting boxes in the pixel space of the map, which can be different from the plot axes. (`#5275 <https://github.com/sunpy/sunpy/pull/5275>`__)
- Added the :func:`~sunpy.coordinates.propagate_with_solar_surface` context manager for transformations, which will automatically apply solar differential rotation when transforming a coordinate between frames with a change in time (``obstime``). (`#5281 <https://github.com/sunpy/sunpy/pull/5281>`__)
- Add support for parsing the observer location from a `~astropy.wcs.WCS` object
  when using the 'OBSGEO' formulation. This is the recommended way to define the
  observer location of a ground based observer. (`#5315 <https://github.com/sunpy/sunpy/pull/5315>`__)
- Added a new function, :meth:`sunpy.visualization.draw_limb`, that draws
  the solar limb as seen from an arbitrary observer coordinate on a world
  coordinate system aware Axes. (`#5414 <https://github.com/sunpy/sunpy/pull/5414>`__)
- `sunpy.map.GenericMap.rsun_meters` now uses `sunpy.map.GenericMap.rsun_obs`
  as a fallback to calculate the assumed radius of emission if RSUN_REF metadata
  isn't present but metadata for `~sunpy.map.GenericMap.rsun_obs` is. (`#5416 <https://github.com/sunpy/sunpy/pull/5416>`__)
- Added :func:`sunpy.coordinates.utils.get_limb_coordinates` to get the solar
  limb coordinates as seen from a given observer. (`#5417 <https://github.com/sunpy/sunpy/pull/5417>`__)
- Printing the response from a `~sunpy.net.Fido` query now includes the URL where
  the data files are sourced from.

  If you develop a third-party `~sunpy.net.Fido` client, support for this can
  be automatically enabled by adding a ``info_url`` property to your
  `~sunpy.net.base_client.BaseClient` that returns a URL as a string. (`#5431 <https://github.com/sunpy/sunpy/pull/5431>`__)
- `~sunpy.timeseries.TimeSeries` can now read CDF files that conform to the
   ISTP/IACG guidelines (https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html). (`#5435 <https://github.com/sunpy/sunpy/pull/5435>`__)
- The properties `~sunpy.map.GenericMap.date_start`,
  `~sunpy.map.GenericMap.date_end`, and `~sunpy.map.GenericMap.date_average` have
  been added to be drawn from the relevant FITS metadata, if present in the map
  header. (`#5449 <https://github.com/sunpy/sunpy/pull/5449>`__)
- Add default color map and normalization for `~sunpy.map.sources.HMISynopticMap`
  The default color map is 'hmimag' and the default normalization is linear between
  -1.5e-3 and +1.5e3, the expected normalization for this particular color map. (`#5464 <https://github.com/sunpy/sunpy/pull/5464>`__)
- The headers produced by :func:`~sunpy.map.make_fitswcs_header` now include ``NAXIS``, ``NAXIS1``, and ``NAXIS2`` keywords. (`#5470 <https://github.com/sunpy/sunpy/pull/5470>`__)
- The `~astropy.wcs.WCS` instance returned by the `sunpy.map.GenericMap.wcs` property now includes the shape of the data array. (`#5470 <https://github.com/sunpy/sunpy/pull/5470>`__)
- Added the method :meth:`sunpy.map.GenericMap.reproject_to` for reprojecting a `~sunpy.map.Map` to a different WCS.
  This method requires the optional package `reproject` to be installed. (`#5470 <https://github.com/sunpy/sunpy/pull/5470>`__)
- Registered the time format ``tai_seconds`` for `astropy.time.Time` (via `~sunpy.time.TimeTaiSeconds`) to support parsing the numerical time format of TAI seconds since 1958-01-01 00:00:00.
  This format includes UTC leap seconds, and enables equivalent functionality to the ``anytim2tai`` routine in SSW. (`#5489 <https://github.com/sunpy/sunpy/pull/5489>`__)
- Added `sunpy.map.sources.WISPRMap` as a map source for WISPR on Parker Solar Probe.
  This improves the `~sunpy.map.GenericMap.name` of the map and adds correct
  information for the `~sunpy.map.GenericMap.processing_level` and
  `~sunpy.map.GenericMap.exposure_time`. (`#5502 <https://github.com/sunpy/sunpy/pull/5502>`__)
- :func:`sunpy.io.fits.write` can now update the ``data`` and ``header`` of an existing HDU instance, as an alternative to creating a new instance of a specified HDU type. This adds support for writing a HDU (such as :class:`~astropy.io.fits.CompImageHDU`) initialised with non-default keyword arguments. (`#5503 <https://github.com/sunpy/sunpy/pull/5503>`__)
- Added `~sunpy.timeseries.GenericTimeSeries.observatory` to provide observatory information for the timeseries e.g. specific goes satellite number. (`#5556 <https://github.com/sunpy/sunpy/pull/5556>`__)
- :meth:`sunpy.timeseries.GenericTimeSeries.plot` and
  :meth:`sunpy.timeseries.GenericTimeSeries.peek` will now automatically label
  the y-axis if all the columns being plotted have the same units. (`#5557 <https://github.com/sunpy/sunpy/pull/5557>`__)
- :meth:`sunpy.timeseries.GenericTimeSeries.plot` and
  :meth:`sunpy.timeseries.GenericTimeSeries.peek` now have an option ``columns``
  that allows plotting a subset of the columns present. (`#5557 <https://github.com/sunpy/sunpy/pull/5557>`__)
- Added a new CDAWeb client, along with helper utilities to `sunpy.net.cdaweb`. (`#5558 <https://github.com/sunpy/sunpy/pull/5558>`__)
- Support for filtering searches with JSOC keywords has been added to ``Fido.search``. (`#5566 <https://github.com/sunpy/sunpy/pull/5566>`__)
- Added support for arithmetic operations between`~sunpy.map.GenericMap` and array-like
  objects. (`#5614 <https://github.com/sunpy/sunpy/pull/5614>`__)
- Added ``quantity`` attribute to `~sunpy.map.GenericMap` to expose the ``data``
  attribute as a `~astropy.units.Quantity` using the ``unit`` attribute. (`#5614 <https://github.com/sunpy/sunpy/pull/5614>`__)


Bug Fixes
---------

- :meth:`sunpy.map.GenericMap.superpixel` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5295 <https://github.com/sunpy/sunpy/pull/5295>`__)
- Inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` in units other than ``u.pix``
  (e.g. ```u.kpix``) are now handled correctly. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Fractional inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` were previously rounded using `int`
  in the superpixel algorithm, but not assigned integer values in the new meatadata.
  This has now been changed so the rounding is correctly reflected in the meatadata. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Remove runtime use of ``astropy.tests.helper.assert_quantity_allclose`` which
  introduces a runtime dependancy on ``pytest``. (`#5305 <https://github.com/sunpy/sunpy/pull/5305>`__)
- :meth:`sunpy.map.GenericMap.resample` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5309 <https://github.com/sunpy/sunpy/pull/5309>`__)
- Fix saving `.GenericMap` to an asdf file with version 2.8.0 of the asdf package. (`#5342 <https://github.com/sunpy/sunpy/pull/5342>`__)
- When the limb is entirely visible, :meth:`sunpy.map.GenericMap.draw_limb` no
  longer plots an invisible patch for the hidden part of the limb and now returns
  `None` instead of the invisible patch. Similarly, when the limb is entirely
  invisible, no patch is drawn for the visible part and `None` is returned
  instead of the visible patch. (`#5414 <https://github.com/sunpy/sunpy/pull/5414>`__)
- :meth:`sunpy.map.GenericMap.plot` now correctly sets axis labels based on the
  coordinate system of the axes, and not the coordinate system of the map
  being plotted. This was previously only an issue if using ``autoalign=True``
  when the Map coordinate system was different to the axes coordinate system. (`#5432 <https://github.com/sunpy/sunpy/pull/5432>`__)
- :meth:`sunpy.map.GenericMap.plot` no longer adds a unit string to the axis
  labels if the axes being plotted on is a WCSAxes. For a WCSAxes, angular units
  are indicated in the tick labels, and automatically change when the zoom level
  changes from e.g. degrees to arc-minutes. This could previously lead to
  situations where the axis label units were incorrect. (`#5432 <https://github.com/sunpy/sunpy/pull/5432>`__)
- Implement automatic fallback to helioviewer mirrors if API is non-functional. (`#5440 <https://github.com/sunpy/sunpy/pull/5440>`__)
- Fixed the incorrect value for the FITS WCS ``LONPOLE`` keyword when using :func:`~sunpy.map.make_fitswcs_header` for certain combinations of WCS projection and reference coordinate. (`#5448 <https://github.com/sunpy/sunpy/pull/5448>`__)
- The date returned by `~sunpy.map.GenericMap.date` for Solar Orbiter/EUI maps
  has been adjusted to be taken from the DATE-AVG keyword
  (the middle of the image acquisition period), instead of the DATE-OBS
  keyword (the beginning of the image acquisition period). This means the observer
  coordinate now has the correct date. (`#5462 <https://github.com/sunpy/sunpy/pull/5462>`__)
- The ``.unit`` attribute for HMI synoptic maps has been fixed. (`#5467 <https://github.com/sunpy/sunpy/pull/5467>`__)
- When "TAI" is in the date string, `sunpy.map.GenericMap.date`
  now only raises a warning if the TIMESYS keyword is present
  and different to "TAI". Previously a warning was raised all the
  time when "TAI" was in the date string. (`#5468 <https://github.com/sunpy/sunpy/pull/5468>`__)
- Fixed a bug where the property `sunpy.map.GenericMap.rsun_meters` would always internally determine the observer location, even when it is not needed, particularly for Stonyhurst heliographic maps, which have no notion of an observer.
  Thus, when working with a Stonyhurst heliographic map, a user could get an irrelevant warning message about having to assume an observer location (Earth center). (`#5478 <https://github.com/sunpy/sunpy/pull/5478>`__)
- Fixed the unintended insertion of (assumed) observer location information when accessing the property `sunpy.map.GenericMap.wcs` for Stonyhurst heliographic maps. (`#5478 <https://github.com/sunpy/sunpy/pull/5478>`__)
- Fixed an incorrect value for the FITS WCS ``LONPOLE`` keyword when using :func:`~sunpy.map.make_fitswcs_header` for `~sunpy.coordinates.frames.Helioprojective` maps with certain values of latitude for the reference coordinate. (`#5490 <https://github.com/sunpy/sunpy/pull/5490>`__)
- A non-standard ``CROTA`` keyword included in a `sunpy.map.sources.EUIMap` FITS header is now renamed to the recommended ``CROTA2`` so a warning is no longer raised. (`#5493 <https://github.com/sunpy/sunpy/pull/5493>`__)
- The plotting x-limits of :meth:`sunpy.timeseries.sources.NOAAIndicesTimeSeries.plot`
  are now adjusted to only include finite points in the timeseries data. (`#5496 <https://github.com/sunpy/sunpy/pull/5496>`__)
- The Hinode/XRT map source now corrects the TIMESYS keyword, fixing the ``.wcs``
  property that was previously broken for Hinode/XRT maps. (`#5508 <https://github.com/sunpy/sunpy/pull/5508>`__)
- Updated `sunpy.map.CompositeMap.plot` to support the ``linestyles`` and ``colors`` arguments, in addition to the existing ``linewidths`` argument. (`#5521 <https://github.com/sunpy/sunpy/pull/5521>`__)
- Fixed a bug where rotating a `~sunpy.map.Map` could result in an extremely small shift (at the numerical-precision level) in the mapping from world coordinates to pixels. (`#5553 <https://github.com/sunpy/sunpy/pull/5553>`__)
- Fixed a bug where rotating a `~sunpy.map.Map` that is missing observation-time metadata could result in an incorrect reference coordinate. (`#5553 <https://github.com/sunpy/sunpy/pull/5553>`__)
- Fix a bug where saving a helioprojective or heliocentric coordinate to an
  asdf file didn't work due to a schema version mismatch if the observer
  location was a fully specified Stonyhurst heliographic coordinate. (`#5584 <https://github.com/sunpy/sunpy/pull/5584>`__)
- `~sunpy.map.sources.XRTMap` uppercases the ``TIMESYS`` key before checking if the
  key needs to be fixed. (`#5592 <https://github.com/sunpy/sunpy/pull/5592>`__)
- Fixed passing a URL to :func:`sunpy.io.read_file` on windows. (`#5601 <https://github.com/sunpy/sunpy/pull/5601>`__)
- Fixed a bug where the ``date`` property on `~sunpy.map.sources.HMISynopticMap` returned ``None``
  if the ``DATE-OBS`` key was present. (`#5648 <https://github.com/sunpy/sunpy/pull/5648>`__)


Documentation
-------------

- Added the gallery example :ref:`sphx_glr_generated_gallery_differential_rotation_comparing_rotation_models.py` to visualize the differences between models of solar differential rotation. (`#5527 <https://github.com/sunpy/sunpy/pull/5527>`__)
- Added an example to how to save out maps as FITS files and load them back in, :ref:`sphx_glr_generated_gallery_saving_and_loading_data_genericmap_in_fits.py`. (`#5544 <https://github.com/sunpy/sunpy/pull/5544>`__)


Internal Changes
----------------

- The `~sunpy.coordinates.frames.Helioprojective` frame now has the convenience property ``angular_radius`` to return the angular radius of the Sun as seen by the observer. (`#5191 <https://github.com/sunpy/sunpy/pull/5191>`__)
- Online tests can now report back status of remote urls and will XFAIL if the remote server is unreachable. (`#5233 <https://github.com/sunpy/sunpy/pull/5233>`__)
- Re-enabled the unit test to check for coordinates consistency with JPL HORIZONS when the matching ephemeris can be specified. (`#5314 <https://github.com/sunpy/sunpy/pull/5314>`__)
- The `~sunpy.timeseries.TimeSeries` factory has been refactored to
  improve readability and maintainability of the internal code. (`#5411 <https://github.com/sunpy/sunpy/pull/5411>`__)
- `sunpy.map.GenericMap.rsun_obs` no longer emits a warning if the metadata it
  looks for is not present. Instead the standard photospheric radius is assumed
  and a log message emitted at the 'info' level. (`#5416 <https://github.com/sunpy/sunpy/pull/5416>`__)
- Nearest-neighbour and linear
  (the default for :meth:`sunpy.map.GenericMap.resample`)
  resampling have been significantly sped up. (`#5476 <https://github.com/sunpy/sunpy/pull/5476>`__)
- `sunpy.map.Map` now raises a clear error when the map is constructed if units
  of either two axes are not angular units. (`#5602 <https://github.com/sunpy/sunpy/pull/5602>`__)


3.0.1 (2021-07-03)
==================

Bug Fixes
---------

- Fixed a bug where `~sunpy.map.GenericMap` used to break with keyword arguments. (`#5392 <https://github.com/sunpy/sunpy/pull/5392>`__)
- Fixed a bug where calling :meth:`sunpy.map.GenericMap.draw_contours` on a different WCS could result in an unnecessary expansion of the plot limits. (`#5398 <https://github.com/sunpy/sunpy/pull/5398>`__)
- Fixed incorrect return values from :func:`~sunpy.map.all_corner_coords_from_map` if a rectangular map was provided. (`#5419 <https://github.com/sunpy/sunpy/pull/5419>`__)
- Do not trigger a pytest import in the asdf plugin for saving sunpy coordinate frames. (`#5429 <https://github.com/sunpy/sunpy/pull/5429>`__)
- Constructing a 2D coordinate in the `~sunpy.coordinates.frames.HeliographicCarrington` frame with ``observer='self'`` now raises an error upon creation.
  When specifying ``observer='self'``, the ``radius`` coordinate component serves as the Sun-observer distance that is necessary to fully define the Carrington heliographic coordinates. (`#5358 <https://github.com/sunpy/sunpy/pull/5358>`__)
- Fixed two bugs with handling the motion of the Sun when transforming between coordinate frames with a change in ``obstime``.
  These bugs did not affect any results if the context manager :func:`~sunpy.coordinates.transform_with_sun_center` had been used. (`#5381 <https://github.com/sunpy/sunpy/pull/5381>`__)
- Fixed a bug where the ``rsun`` frame attribute could be unintentionally reset to the default value during transformation.
  This bug primarily affected the transformation of a `~sunpy.coordinates.frames.Helioprojective` coordinate to a `~sunpy.coordinates.frames.HeliographicStonyhurst` frame. (`#5395 <https://github.com/sunpy/sunpy/pull/5395>`__)
- Fixed a bug where creating a `~sunpy.coordinates.frames.HeliographicStonyhurst` frame or a `~sunpy.coordinates.frames.HeliographicCarrington` frame from WCS information failed to make use of any specified ``rsun_ref`` value. (`#5395 <https://github.com/sunpy/sunpy/pull/5395>`__)
- `~sunpy.map.sources.SXTMap` now always returns `None` for the ``wavelength`` attribute. Previously this raised an error. (`#5401 <https://github.com/sunpy/sunpy/pull/5401>`__)


Added/Improved Documentation
----------------------------

- Simplified the "Downloading LASCO C2" gallery example by removing redundant modifications to the metadata before it is loaded by `~sunpy.map.Map`. (`#5402 <https://github.com/sunpy/sunpy/pull/5402>`__)
- Tided up the HMI synoptic map example by removing redundant code and correcting some of the comments. (`#5413 <https://github.com/sunpy/sunpy/pull/5413>`__)

3.0.0 (2021-05-14)
==================

Backwards Incompatible Changes
------------------------------

- ``sunpy.instr`` has been depreacted and will be removed in sunpy 3.1 in favour of `sunkit_instruments`.
  The code that is under ``sunpy.instr`` is imported via `sunkit_instruments` to ensure backwards comparability. (`#4526 <https://github.com/sunpy/sunpy/pull/4526>`__)
- Several `sunpy.map.GenericMap` attributes have been updated to return `None` when the relevant piece of FITS metadata is missing. These are:

  - `~sunpy.map.GenericMap.exposure_time`, previously defaulted to zero seconds.
  - `~sunpy.map.GenericMap.measurement`, previously defaulted to zero.
  - `~sunpy.map.GenericMap.waveunit`, previously defaulted to ``u.one``.
  - `~sunpy.map.GenericMap.wavelength`, previously defaulted to zero. (`#5126 <https://github.com/sunpy/sunpy/pull/5126>`__)
- `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` no longer automatically convert 2D input to a 3D coordinate during instantiation.
  Instead, the 2D-to-3D conversion is deferred until the coordinate is transformed to a different frame, or with a call to the method :meth:`~sunpy.coordinates.frames.BaseHeliographic.make_3d`. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- Changed URL for the `sunpy.net.dataretriever.sources.noaa.SRSClient` from "ftp://ftp.swpc.noaa.gov/pub/warehouse/" to "ftp://ftp.ngdc.noaa.gov/STP/swpc_products/daily_reports/".
  The old URL is unsupported and we expect the files will be the same but we can not say with 100% certainty. (`#5173 <https://github.com/sunpy/sunpy/pull/5173>`__)
- Changed `sunpy.net.attrs.Source` to `sunpy.net.attrs.Provider` for the `sunpy.net.dataretriever.sources.gong.GONGClient`. (`#5174 <https://github.com/sunpy/sunpy/pull/5174>`__)
- The ``rsun`` frame attribute of `~sunpy.coordinates.frames.Helioprojective` now converts any input to kilometers. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- :meth:`sunpy.map.CompositeMap.plot` now internally calls :meth:`sunpy.map.GenericMap.plot` and :meth:`sunpy.map.GenericMap.draw_contours`, which may affect the plot output of existing user code. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)
- Removed the ``basic_plot`` keyword argument from :meth:`sunpy.map.CompositeMap.peek` due to its unreliability. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)
- ``sunpy.util.sphinx.changelog`` and ``sunpy.util.towncrier`` have been removed and are now in a standalone package `sphinx-changelog <https://github.com/openastronomy/sphinx-changelog>`__. (`#5049 <https://github.com/sunpy/sunpy/pull/5049>`__)


Deprecations and Removals
-------------------------

- Deprecated ``sunpy.map.GenericMap.draw_rectangle`` in favor of :meth:`~sunpy.map.GenericMap.draw_quadrangle`. (`#5236 <https://github.com/sunpy/sunpy/pull/5236>`__)
- Using `~sunpy.map.GenericMap` plotting methods on an `~matplotlib.axes.Axes` that is not a `~astropy.visualization.wcsaxes.WCSAxes` is deprecated.
  This previously raised a warning, but is now formally deprecated, and will raise an error in sunpy 3.1. (`#5244 <https://github.com/sunpy/sunpy/pull/5244>`__)
- Deprecated ``sunpy.roi.chaincode.Chaincode`` and created a replacement at `sunpy.net.helio.Chaincode`.

  This replacement has the following changes:

  1. Added support for numpy array as an input (it was broken before).
  2. Renamed ``BoundingBox`` to ``boundingbox``
  3. Renamed ``subBoundingBox`` to ``sub_boundingbox``
  4. Now area and length raise `NotImplementedError` (`#5249 <https://github.com/sunpy/sunpy/pull/5249>`__)
- Deprecated ``sunpy.roi.roi``, as it currently has no obvious use and has never seen any real development work. (`#5249 <https://github.com/sunpy/sunpy/pull/5249>`__)


Features
--------

- :func:`sunpy.coordinates.get_horizons_coord` can now be given a start time, end time,
  and number of intervals (or interval length) to query a evenly spaced set of
  times. See the documentation string for more information and an example. (`#4698 <https://github.com/sunpy/sunpy/pull/4698>`__)
- Added :meth:`sunpy.map.GenericMap.draw_quadrangle` for drawing a quadrangle on a map.
  A quadrangle has edges that are aligned with lines of constant latitude and longitude, but these can be in a different coordinate system than that of the map. (`#4809 <https://github.com/sunpy/sunpy/pull/4809>`__)
- Added a ``longitude`` keyword argument to :func:`~sunpy.coordinates.sun.carrington_rotation_time` as an alternate way to specify a fractional Carrington rotation. (`#4879 <https://github.com/sunpy/sunpy/pull/4879>`__)
- Colorbar in `sunpy.map.GenericMap.peek` now has a unit label. (`#4930 <https://github.com/sunpy/sunpy/pull/4930>`__)
- The default axes used by ``BaseFuncAnimator.get_animation()``
  is now ``BaseFuncAnimator.axes``, instead of the currently active axes (accessed via.
  :func:`matplotlib.pyplot.gca`). The allows animations to be created on figures
  created directly using `matplotlib.figure.Figure`.

  To revert to the previous behaviour of using the current axes,
  give ``axes=plt.gca()`` to ``get_animation()``. (`#4968 <https://github.com/sunpy/sunpy/pull/4968>`__)
- Added colormaps for Solar Orbiter EUI images. These are used automatically
  when an EUI image is loaded. (`#5023 <https://github.com/sunpy/sunpy/pull/5023>`__)
- Added the ability to dynamically scale `sunpy.visualization.animator` instances.
  By specifying the ``clip_interval`` keyword, it will now clip the minimum and maximum at each slider step to the specified interval. (`#5025 <https://github.com/sunpy/sunpy/pull/5025>`__)
- Added a ``sunpy.time.timerange.TimeRange.__contains__`` method to `sunpy.time.TimeRange`
  that tests if two time ranges overlap. (`#5093 <https://github.com/sunpy/sunpy/pull/5093>`__)
- Added the ability to namespace files downloaded using `sunpy.data.data_manager.manager.DataManager` by prepending the file name with module name. (`#5111 <https://github.com/sunpy/sunpy/pull/5111>`__)
- Added a rigid rotation model to :func:`~sunpy.physics.differential_rotation.diff_rot` via ``rot_type=rigid``, where the rotation rate does not vary with latitude. (`#5132 <https://github.com/sunpy/sunpy/pull/5132>`__)
- Added a :meth:`~sunpy.map.MapSequence.save` method to `sunpy.map.MapSequence`
  that saves each map of the sequence. (`#5145 <https://github.com/sunpy/sunpy/pull/5145>`__)
- The allowable ``level`` inputs to :meth:`sunpy.map.GenericMap.contour` and
  :meth:`sunpy.map.GenericMap.draw_contours` have been consolidated. Both methods
  now accept
  - Scalars, if the map has no units
  - Quantities, if the map has units
  - Percentages (`#5154 <https://github.com/sunpy/sunpy/pull/5154>`__)
- Added support for corrected NOAA SWPC solar region summary data files. (`#5173 <https://github.com/sunpy/sunpy/pull/5173>`__)
- Updated ``sunpy.util.sysinfo.system_info`` to return all optional dependencies of sunpy. (`#5175 <https://github.com/sunpy/sunpy/pull/5175>`__)
- `sunpy.map.Map` now supports the EUI instrument on Solar Orbiter. (`#5210 <https://github.com/sunpy/sunpy/pull/5210>`__)
- `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` now have an ``rsun`` frame attribute to specify the radius of the Sun, which defaults to the photospheric radius defined in `sunpy.sun.constants`.
  This frame attribute is used when converting a 2D coordinate (longitude and latitude, with no specified radial distance) to a 3D coordinate by setting the radial distance to ``rsun`` (i.e., the assumption is that the coordinate is on the surface of the Sun). (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- Enhanced :meth:`sunpy.map.GenericMap.draw_limb` so that the solar limb can be plotted on axes that correspond to a different map (e.g., with a different observer).
  The part of the limb that is not visible to the axes's observer because it is on the far side of the Sun is shown as dotted rather than solid. (`#5237 <https://github.com/sunpy/sunpy/pull/5237>`__)
- `~sunpy.util.MetaDict` now saves a copy of the metadata on creation, which can
  be accessed using the `~sunpy.util.MetaDict.original_meta` property. Three
  new properties have also been added to query any changes that have been made
  to metadata:

  - `~sunpy.util.MetaDict.added_items`
  - `~sunpy.util.MetaDict.removed_items`
  - `~sunpy.util.MetaDict.modified_items`

  As an example, ``my_map.meta.modified_items`` will return a dictionary mapping
  keys to their original value and current value. (`#5241 <https://github.com/sunpy/sunpy/pull/5241>`__)
- Added :func:`sunpy.map.contains_coordinate` which provides a quick way to see if a
  world coordinate is contained within the array bounds of a map. (`#5252 <https://github.com/sunpy/sunpy/pull/5252>`__)
- Added an optional keyword argument ``autoalign`` to :meth:`sunpy.map.GenericMap.plot` for plotting a map to axes that correspond to a different WCS.
  See :ref:`sphx_glr_generated_gallery_map_transformations_autoalign_aia_hmi.py`. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)
- :meth:`sunpy.map.CompositeMap.plot` now properly makes use of WCS information to position and orient maps when overlaying them. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)


Bug Fixes
---------

- Fixed the drawing methods of `sunpy.map.GenericMap` (e.g., ``~sunpy.map.GenericMap.draw_rectangle``) so that any text labels will appear in the legend. (`#5019 <https://github.com/sunpy/sunpy/pull/5019>`__)
- Fixed bug in ``sunpy.until.scraper.Scraper`` which caused URL patterns containing backslashes to be incorrectly parsed on Windows. (`#5022 <https://github.com/sunpy/sunpy/pull/5022>`__)
- Constructing a `~sunpy.util.MetaDict` is now more lenient, and accepts
  any class that inherits from `collections.abc.Mapping`. This fixes a
  regression where headers read with `astropy.io.fits` raised an error when
  passed to individual `~sunpy.map` sources. (`#5047 <https://github.com/sunpy/sunpy/pull/5047>`__)
- Added warning to :meth:`sunpy.map.GenericMap.rotate` when specified ``missing`` value is not compatible
  with the number type of the data array. (`#5051 <https://github.com/sunpy/sunpy/pull/5051>`__)
- Prevented some colormaps being accidentally modified depending on the order
  and method through which they were accessed. (`#5054 <https://github.com/sunpy/sunpy/pull/5054>`__)
- Reverted change for `sunpy.map.GenericMap.draw_limb` that made it use "add_artist" as it was changing the FOV of the plotted image. (`#5069 <https://github.com/sunpy/sunpy/pull/5069>`__)
- Fixed a bug where some `~sunpy.coordinates.metaframes.RotatedSunFrame` transformations could fail with an ``observer=None`` error. (`#5084 <https://github.com/sunpy/sunpy/pull/5084>`__)
- Fixed bug where `sunpy.data.data_manager.DataManager` would fail to recover upon deleting the sqlite database file. (`#5089 <https://github.com/sunpy/sunpy/pull/5089>`__)
- Fixed a bug where coordinate frames were considered different due to an unintended time difference during time handling at the level of numerical precision (i.e., tens of picoseconds).
  This resulted in the unexpected use of transformation machinery when transforming a coordinate to its own coordinate frame. (`#5127 <https://github.com/sunpy/sunpy/pull/5127>`__)
- Fixed a bug with failing downloads in 2010 with the `~sunpy.net.dataretriever.sources.noaa.SRSClient`. (`#5159 <https://github.com/sunpy/sunpy/pull/5159>`__)
- If the property `sunpy.map.GenericMap.rsun_obs` needs to calculate the solar angular radius from header information, it now properly uses the ``rsun_ref`` keyword if it is present and does not emit any warning. (`#5172 <https://github.com/sunpy/sunpy/pull/5172>`__)
- Added a "rsun_obs" keyword to the output of :func:`sunpy.map.make_fitswcs_header` if the coordinate argument has a "rsun" frame attribute. (`#5177 <https://github.com/sunpy/sunpy/pull/5177>`__)
- Fixed small inaccuracies in the grid plotted by :meth:`~sunpy.map.GenericMap.draw_grid` for maps that specify a radius of the Sun that is different from the constant in `sunpy.sun.constants`. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- Fixed :meth:`sunpy.map.GenericMap.draw_contours` so that the contours from a map can be plotted on axes with a different coordinate system. (`#5239 <https://github.com/sunpy/sunpy/pull/5239>`__)
- When using the cylindrical representation of ``Heliocentric`` to work in the Heliocentric Radial coordinate frame, the ``psi`` component now goes from 0 to 360 degrees instead of -180 to 180 degrees. (`#5242 <https://github.com/sunpy/sunpy/pull/5242>`__)
- Changed ``MDIMap`` to use the "CONTENT" keyword to identify the measurement, similar to ``HMIMap``, and removed the special-case nickname. This fixes the broken title on plots. (`#5257 <https://github.com/sunpy/sunpy/pull/5257>`__)
- :func:`sunpy.coordinates.solar_frame_to_wcs_mapping` now sets the observer auxiliary
  information when a `~sunpy.coordinates.HeliographicCarrington` frame with
  ``observer='self'`` is passed. (`#5264 <https://github.com/sunpy/sunpy/pull/5264>`__)
- Calling :func:`sunpy.map.make_fitswcs_header` with a
  `~sunpy.coordinates.HeliographicCarrington` coordinate that with ``observer='self'``
  set now correctly sets the observer information in the header. (`#5264 <https://github.com/sunpy/sunpy/pull/5264>`__)
- :meth:`sunpy.map.GenericMap.superpixel` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5295 <https://github.com/sunpy/sunpy/pull/5295>`__)
- Inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` in units other than ``u.pix``
  (e.g. ```u.kpix``) are now handled correctly. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Fractional inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` were previously rounded using `int`
  in the superpixel algorithm, but not assigned integer values in the new meatadata.
  This has now been changed so the rounding is correctly reflected in the meatadata. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Remove runtime use of ``astropy.tests.helper.assert_quantity_allclose`` which
  introduces a runtime dependancy on ``pytest``. (`#5305 <https://github.com/sunpy/sunpy/pull/5305>`__)
- :meth:`sunpy.map.GenericMap.resample` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5309 <https://github.com/sunpy/sunpy/pull/5309>`__)
- Fix saving `.GenericMap` to an asdf file with version 2.8.0 of the asdf package. (`#5342 <https://github.com/sunpy/sunpy/pull/5342>`__)


Added/Improved Documentation
----------------------------

- Added a gallery example (:ref:`sphx_glr_generated_gallery_plotting_plot_rectangle.py`) for drawing rectangles on maps. (`#4528 <https://github.com/sunpy/sunpy/pull/4528>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_plotting_wcsaxes_plotting_example.py`)
  of how pixel and SkyCoords work when plotted with `~astropy.visualization.wcsaxes`. (`#4867 <https://github.com/sunpy/sunpy/pull/4867>`__)
- Added a gallery example  (:ref:`sphx_glr_generated_gallery_plotting_plotting_blank_map.py`) on how to create a blank map and mark locations. (`#5077 <https://github.com/sunpy/sunpy/pull/5077>`__)
- Added a gallery example (:ref:`sphx_glr_generated_gallery_plotting_hmi_cutout.py`)
  demonstrating how to add a HMI zoomed-in region next to a full disk HMI image. (`#5090 <https://github.com/sunpy/sunpy/pull/5090>`__)
- Updated the :ref:`sphx_glr_generated_gallery_computer_vision_techniques_mask_disk.py` example to generate the mask using :func:`sunpy.map.coordinate_is_on_solar_disk`. (`#5114 <https://github.com/sunpy/sunpy/pull/5114>`__)
- Added a gallery example (:ref:`sphx_glr_generated_gallery_map_map_segment.py`)
  demonstrating how to create a segment of a particular map from transformed coordinates. (`#5121 <https://github.com/sunpy/sunpy/pull/5121>`__)
- For the various subclasses of `~sunpy.map.GenericMap` (e.g., `~sunpy.map.sources.AIAMap`), the online documentation now shows all of the inherited attributes and methods. (`#5142 <https://github.com/sunpy/sunpy/pull/5142>`__)
- Added a documentation string to `~sunpy.map.sources.sdo.HMISynopticMap`. (`#5186 <https://github.com/sunpy/sunpy/pull/5186>`__)
- Added a new gallery example showcasing how to overlay HMI contours on an AIA image. (`#5229 <https://github.com/sunpy/sunpy/pull/5229>`__)


Trivial/Internal Changes
------------------------

- Replaced the old test runner with a new version that adds a dependency check before the test suite is run. (`#4596 <https://github.com/sunpy/sunpy/pull/4596>`__)
- The testing suite now raises a warning if the `~matplotlib.pyplot` figure stack is not empty prior to running a test, and it closes all open figures after finishing each test. (`#4969 <https://github.com/sunpy/sunpy/pull/4969>`__)
- Improved performance when moving the slider in
  ``sunpy.visualisation.animator.ArrayAnimatorWCS``. (`#4971 <https://github.com/sunpy/sunpy/pull/4971>`__)
- Added some basic logging to HEK searches, at the 'debug' logging level. (`#5020 <https://github.com/sunpy/sunpy/pull/5020>`__)
- Refactored `~sunpy.coordinates.metaframes.RotatedSunFrame` transformations for improved performance. (`#5084 <https://github.com/sunpy/sunpy/pull/5084>`__)
- Re-ordered keyword-only arguments of ``sunpy.map.GenericMap.draw_rectangle`` to match :meth:`sunpy.map.GenericMap.submap`. (`#5091 <https://github.com/sunpy/sunpy/pull/5091>`__)
- Significantly sped up calls to :func:`~sunpy.time.parse_time` for string
  arguments. This will have knock on effects, including improved performance of
  querying the VSO. (`#5108 <https://github.com/sunpy/sunpy/pull/5108>`__)
- Added tests for ``sunpy.visualization.animator.mapsequenceanimator`` and :meth:`sunpy.map.MapSequence.plot`. (`#5125 <https://github.com/sunpy/sunpy/pull/5125>`__)
- The ``CROTA`` keywords are no longer set on `sunpy.map.GenericMap.wcs`, as the
  ``PC_ij`` keywords are always set and the FITS standard says that these keywords
  must not co-exist. (`#5166 <https://github.com/sunpy/sunpy/pull/5166>`__)
- Temporarily disabled the unit test to check for coordinates consistency with JPL HORIZONS due to the inability to choose a matching ephemeris. (`#5203 <https://github.com/sunpy/sunpy/pull/5203>`__)
- :func:`~sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay` now accepts ``obstime`` and ``rsun`` optional arguments.
  This function is not typically called directly by users. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- `~sunpy.map.GenericMap` plotting methods now have consistent argument
  checking for the ``axes`` argument, and will raise the same warnings
  or errors for similar ``axes`` input. (`#5223 <https://github.com/sunpy/sunpy/pull/5223>`__)
- Calling :meth:`sunpy.map.GenericMap.plot` on a
  `~astropy.visualization.wcsaxes.WCSAxes` with a different
  World Coordinate System (WCS) to the map now raises a warning,
  as the map data axes may not correctly align with the coordinate axes.
  This happens if an `~matplotlib.axes.Axes` is created with a projection
  that is a different map to the one being plotted. (`#5244 <https://github.com/sunpy/sunpy/pull/5244>`__)
- Re-enabled the unit test to check for coordinates consistency with JPL HORIZONS when the matching ephemeris can be specified. (`#5314 <https://github.com/sunpy/sunpy/pull/5314>`__)


2.1.0 (2020-02-21)
==================

Backwards Incompatible Changes
------------------------------

- Support for Python 3.6 and Numpy 1.15 has been dropped in line with
  `NEP 29 <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_.
  The minimum supported version of Astropy is now 4.0, and the minimum version of scipy is now 1.2. (`#4284 <https://github.com/sunpy/sunpy/pull/4284>`__)
- Changed :func:`sunpy.coordinates.sun.B0` return type from `~astropy.coordinates.Angle`
  to `~astropy.coordinates.Latitude`. (`#4323 <https://github.com/sunpy/sunpy/pull/4323>`__)
- An error is now raised if ``vmin`` or ``vmax`` are passed to
  to `sunpy.map.GenericMap.plot` and they are already set on the map ``norm``.
  This is consistent with upcoming Matplotlib changes. (`#4328 <https://github.com/sunpy/sunpy/pull/4328>`__)
- Previously slicing the result of ``Fido.search()`` (a `~sunpy.net.fido_factory.UnifiedResponse` object) so that
  it had a length of one returned another `~sunpy.net.fido_factory.UnifiedResponse` object.
  Now it will return a `~sunpy.net.base_client.QueryResponseTable` object, which is a subclass
  of `astropy.table.Table`. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- The ``.size`` property of a coordinate frame with no associated data will now raise
  an error instead of returning 0. (`#4577 <https://github.com/sunpy/sunpy/pull/4577>`__)
- The following `~sunpy.map.Map` methods have had support for specific positional
  arguments removed. They must now be passed as keyword arguments
  (i.e. ``m.method(keyword_arg=value)``).

  - :meth:`~sunpy.map.GenericMap.submap`: ``width``, ``height``.
  - ``sunpy.map.GenericMap.draw_rectangle``: ``width``, ``height``, ``axes``, ``top_right``.
    (`#4616 <https://github.com/sunpy/sunpy/pull/4616>`__)
- The sunpy specific attributes ``.heliographic_observer`` and ``.rsun``
  are no longer set on the `~astropy.wcs.WCS` returned by `sunpy.map.GenericMap.wcs`. (`#4620 <https://github.com/sunpy/sunpy/pull/4620>`__)
- Due to upstream changes, the parsing logic for the `~sunpy.net.helio.HECClient` now returns
  strings and not bytes for :meth:`~sunpy.net.helio.HECClient.get_table_names`. (`#4643 <https://github.com/sunpy/sunpy/pull/4643>`__)
- Reduced the selection of dependent packages installed by default via ``pip``,
  which means that some of our sub-packages will not fully import when sunpy is installed with
  ``pip install "sunpy"``.
  You can install all dependencies by specifying ``pip install "sunpy[all]"``,
  or you can install sub-package-specific dependencies by specifying, e.g.,
  ``[map]`` or ``[timeseries]``. (`#4662 <https://github.com/sunpy/sunpy/pull/4662>`__)
- The class inheritance for `~sunpy.coordinates.metaframes.RotatedSunFrame` and the frames it
  creates has been changed in order to stop depending on unsupported behavior in the underlying machinery.
  The return values for some :func:`isinstance`/:func:`issubclass` calls will be different,
  but the API for `~sunpy.coordinates.metaframes.RotatedSunFrame` is otherwise unchanged. (`#4691 <https://github.com/sunpy/sunpy/pull/4691>`__)
- Fix a bug in `~sunpy.map.GenericMap.submap` where only the top right and bottom
  left coordinates of the input rectangle in world coordinates were considered
  when calculating the pixel bounding box. All four corners are once again taken
  into account now, meaning that `~sunpy.map.GenericMap.submap` correctly returns
  the smallest pixel box which contains all four corners of the input rectangle.

  To revert to the previous 2.0.0 behaviour, first convert the top right and bottom
  left coordinates to pixel space before calling submap with::

      top_right = smap.wcs.world_to_pixel(top_right) * u.pix
      bottom_left = smap.wcs.world_to_pixel(bottom_left) * u.pix
      smap.submap(bottom_left=bottom_left, top_right=top_right)

  This will define the rectangle in pixel space. (`#4727 <https://github.com/sunpy/sunpy/pull/4727>`__)
- VSO results where the size was ``-1`` (missing data) now return ``None`` rather
  than ``-1`` to be consistent with other missing data in the VSO results. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- All result objects contained within the results of a ``Fido.search()`` (a
  `~sunpy.net.fido_factory.UnifiedResponse` object) are now
  `~sunpy.net.base_client.QueryResponseTable` objects (or subclasses thereof).
  These objects are subclasses of `astropy.table.Table` and can therefore be
  filtered and inspected as tabular objects, and the modified tables can be passed
  to ``Fido.fetch``.

  This, while a breaking change for anyone accessing these response objects
  directly, will hopefully make working with ``Fido`` search results much easier. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- Results from the `~sunpy.net.dataretriever.NOAAIndicesClient` and the
  `~sunpy.net.dataretriever.NOAAPredictClient` no longer has ``Start Time`` or
  ``End Time`` in their results table as the results returned from the client are
  not dependant upon the time parameter of a search. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The ``sunpy.net.vso.QueryResponse.search`` method has been removed as it has not
  worked since the 1.0 release of sunpy. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The ``sunpy.net.hek.hek.HEKColumn`` class has been removed, the ``HEKTable`` class
  now uses the standard `astropy.table.Column` class. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The keys used to format file paths in ``Fido.fetch`` have changed. They are now
  more standardised across all the clients, as they are all extracted from the
  names of the columns in the results table.

  For results from the VSO the keys are no longer separated with ``.``, and are
  based on the displayed column names. For results from the ``dataretriever``
  clients the only main change is that the keys are now lower case, where they
  were capitilized before. You can use the ``.sunpy.net.fido_factory.UnifiedResponse.path_format_keys``
  method to see all the possible keys for a particular search. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The time returned from :func:`~sunpy.coordinates.sun.carrington_rotation_number`
  has been changed from the TT scale to the more common UTC scale. To undo this change,
  use ``time_out = time_out.tt`` on the outputted time. (`#4819 <https://github.com/sunpy/sunpy/pull/4819>`__)
- `~.BaseQueryResponse.response_block_properties` has been renamed to
  ``.BaseQueryResponse.path_format_keys``, on the return objects from all
  ``search()`` methods on all clients and from ``Fido.search()``. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)

Removals
--------
- Removed deprecated functions:

  - ``sunpy.coordinates.frames.Helioprojective.calculate_distance``, alternative
    is `sunpy.coordinates.frames.Helioprojective.make_3d`.
  - ``sunpy.image.coalignment.repair_image_nonfinite`` - if you wish to repair the image,
    this has to be done manually before calling the various `sunpy.image.coalignment` functions.
  - The ``repair_nonfinite`` keyword argument to ``calculate_shift`` and  ``calculate_match_template_shift``
    has been removed.
  - ``sunpy.instr.lyra.download_lytaf_database`` - this just downloaded the file
    at ``http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db``, which can be done manually.
  - ``sunpy.util.net.check_download_file``, no alternative.
  - ``sunpy.visualization.animator.ImageAnimatorWCS``, alternative is
    ``sunpy.visualization.animator.ArrayAnimatorWCS``. (`#4350 <https://github.com/sunpy/sunpy/pull/4350>`__)

- Removed deprecated function ``sunpy.instr.aia.aiaprep``.
  Alternative is `~aiapy.calibrate.register` for converting AIA
  images from level 1 to level 1.5. (`#4485 <https://github.com/sunpy/sunpy/pull/4485>`__)
- ``sunpy.cm`` has been removed. All of the functionality in this module can
  now be found in `sunpy.visualization.colormaps`. (`#4488 <https://github.com/sunpy/sunpy/pull/4488>`__)
- ``sunpy.test.hash`` has been removed, the functionality has been moved into the
  `pytest-mpl <https://github.com/matplotlib/pytest-mpl>`__ package. (`#4605 <https://github.com/sunpy/sunpy/pull/4605>`__)
- ``sunpy.util.multimethod`` has been removed. (`#4614 <https://github.com/sunpy/sunpy/pull/4614>`__)
- The ``lytaf_path`` argument (which previously did nothing) has been removed from
  - ``sunpy.instr.lyra.remove_lytaf_events_from_timeseries``
  - ``sunpy.instr.lyra.get_lytaf_events``
  - ``sunpy.instr.lyra.get_lytaf_event_types`` (`#4615 <https://github.com/sunpy/sunpy/pull/4615>`__)

Deprecations
------------
- Deprecated ``sunpy.net.vso.attrs.Source`` and ``sunpy.net.vso.attrs.Provider``.
  They are now `sunpy.net.attrs.Source` and `sunpy.net.attrs.Provider` respectively.
  (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- Deprecated the use of the ``sunpy.map.GenericMap.size`` property,
  use ``sunpy.map.Map.data.size`` instead. (`#4338 <https://github.com/sunpy/sunpy/pull/4338>`__)
- ``sunpy.net.helio.HECClient.time_query`` is deprecated, `~sunpy.net.helio.HECClient.search`
  is the replacement. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- ``sunpy.net.jsoc.attrs.Keys`` is deprecated; all fields are returned by default and can be filtered post search. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- ``sunpy.net.hek.attrs.Time`` is deprecated; `~sunpy.net.attrs.Time` should be used instead. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- Support for :func:`sunpy.coordinates.wcs_utils.solar_wcs_frame_mapping` to
  use the ``.heliographic_observer`` and ``.rsun`` attributes on a
  `~astropy.wcs.WCS` is depreacted. (`#4620 <https://github.com/sunpy/sunpy/pull/4620>`__)
- The ``origin`` argument to `sunpy.map.GenericMap.pixel_to_world` and
  `sunpy.map.GenericMap.world_to_pixel` is deprecated.

  - If passing ``0``, not using the ``origin`` argument will have the same effect.
  - If passing ``1``, manually subtract 1 pixel from the input to ``pixel_to_world``,
    or manually add 1 pixel to the output of ``world_to_pixel``, and do not use the
    ``origin`` argument. (`#4700 <https://github.com/sunpy/sunpy/pull/4700>`__)
- The ``.VSOClient.link`` method is deprecated as it is no longer used. (`#4789 <https://github.com/sunpy/sunpy/pull/4789>`__)
- The ``.UnifiedResponse.get_response``, ``.UnifiedResponse.tables`` and
  ``.UnifiedResponse.responses`` attributes of ``.UnifiedResponse`` have been
  deprecated as they are no longer needed now the object returns the table
  objects it contains when sliced. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- :meth:`sunpy.net.vso.VSOClient.search` has a new keyword argument
  ``response_type=`` which controls the return type from the ``search()`` method.
  In sunpy 2.1 and 3.0 it will default to the ``"legacy"`` response format, in
  3.1 it will default to the new ``"table"`` response format, and the
  ``"legacy"`` format may be deprecated and removed at a later date.

  Searches made with ``Fido`` will use the new ``"table"`` response format, so
  this only affects users interacting with the ``VSOClient`` object directly. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)


Features
--------


- For :meth:`sunpy.map.GenericMap.quicklook` and :meth:`sunpy.map.MapSequence.quicklook` (also used for the HTML reprsentation shown in Jupyter notebooks), the histogram is now shaded corresponding to the colormap of the plotted image.
  Clicking on the histogram will toggle an alternate version of the histogram. (`#4931 <https://github.com/sunpy/sunpy/pull/4931>`__)
- Add an ``SRS_TABLE`` file to the sample data, and use it in the magnetogram
  plotting example. (`#4993 <https://github.com/sunpy/sunpy/pull/4993>`__)
- Added a `sunpy.map.GenericMap.contour()` method to find the contours on a map. (`#3909 <https://github.com/sunpy/sunpy/pull/3909>`__)
- Added a context manager (:meth:`~sunpy.coordinates.frames.Helioprojective.assume_spherical_screen`)
  to interpret `~sunpy.coordinates.frames.Helioprojective` coordinates as being on
  the inside of a spherical screen instead of on the surface of the Sun. (`#4003 <https://github.com/sunpy/sunpy/pull/4003>`__)
- Added `sunpy.map.sources.HMISynopticMap` for handling the Synoptic maps from HMI. (`#4053 <https://github.com/sunpy/sunpy/pull/4053>`__)
- Added a `~sunpy.map.sources.MDISynopticMap` map source class. (`#4054 <https://github.com/sunpy/sunpy/pull/4054>`__)
- Created `~sunpy.net.dataretriever.GONGClient` for accessing magnetogram synoptic map archives of NSO-GONG. (`#4055 <https://github.com/sunpy/sunpy/pull/4055>`__)
- All coordinate frames will now show the velocity if it exists in the underlying data. (`#4102 <https://github.com/sunpy/sunpy/pull/4102>`__)
- The ephemeris functions :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst()`, :func:`~sunpy.coordinates.ephemeris.get_earth()`, and :func:`~sunpy.coordinates.ephemeris.get_horizons_coord()` can now optionally return the body's velocity as part of the output coordinate. (`#4102 <https://github.com/sunpy/sunpy/pull/4102>`__)
- `~sunpy.util.metadata.MetaDict` now maintains coherence between its keys and their corresponding keycomments. Calling ``del`` on a ``MetaDict`` object key is now case-insensitive. (`#4129 <https://github.com/sunpy/sunpy/pull/4129>`__)
- Allow ``sunpy.visualization.animator.ArrayAnimatorWCS`` to disable ticks for
  a coordinate, by setting ``ticks: False`` in the ``coord_params`` dictionary. (`#4270 <https://github.com/sunpy/sunpy/pull/4270>`__)
- Added a ``show()`` method for `~sunpy.net.base_client.BaseQueryResponse` which returns `~astropy.table.Table` with specified columns for the Query Response. (`#4309 <https://github.com/sunpy/sunpy/pull/4309>`__)
- Added ``_extract_files_meta`` method in ``sunpy.util.scraper.Scraper`` which allows scraper to extract metadata from the file URLs retrieved for a given time range. (`#4313 <https://github.com/sunpy/sunpy/pull/4313>`__)
- Refactoring of `~sunpy.net.dataretriever` which adds these capabilities to `~sunpy.net.dataretriever.QueryResponse`:

  - Any ``attr`` shall not be defaulted to a hard-coded value in all subclasses of `~sunpy.net.dataretriever.GenericClient`; thus records for all possible ``attrs`` shall be returned if it is not specified in the query.
  - `~sunpy.net.dataretriever.QueryResponse` can now show more columns; thus all metadata extractable from matching file URLs shall be shown and for a client, non-spported ``attrs`` shall not be shown in the response tables. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- New class attributes added to `~sunpy.net.dataretriever.GenericClient`:

  - ``baseurl`` and ``pattern`` which are required to define a new simple client.
  - ``optional`` and ``required`` which are a ``set`` of optional and required `~sunpy.net.attrs` respectively; which generalizes :meth:`~sunpy.net.dataretriever.GenericClient._can_handle_query`. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- Additions in ``sunpy.util.scraper`` to support the refactoring of `~sunpy.net.dataretriever.GenericClient`:
  - ``sunpy.util.scraper.Scraper.findDatewith_extractor`` that parses the url using extractor to return its start time.
  - A ``matcher`` in ``sunpy.util.scraper.Scraper._extract_files_meta`` which validates the extracted metadata by using the dictionary returned from :meth:`~sunpy.net.dataretriever.GenericClient._get_match_dict`. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- Added methods :meth:`~sunpy.net.dataretriever.GenericClient.pre_search_hook` and :meth:`~sunpy.net.dataretriever.GenericClient.post_search_hook` which helps to translate the attrs for scraper before and after the search respectively. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- :meth:`sunpy.timeseries.sources.RHESSISummaryTimeSeries.peek` has had the following minor
  changes:

  - Colors from the default matplotlib color cycle are now used (but the colors remain qualitatively the same)
  - The default matplotlib linewidth is now used
  - It is now possible to pass in a user specified linewidth
  - Seconds have been added to the x-axis labels (previously it was just hours and minutes) (`#4326 <https://github.com/sunpy/sunpy/pull/4326>`__)
- `~sunpy.net.helio.hec.HECClient` and  `~sunpy.net.hek.hek.HEKClient` now inherit `~sunpy.net.base_client.BaseClient` which makes them compatible with the `~sunpy.net.fido_factory.UnifiedDownloaderFactory` (``Fido``). (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- `~sunpy.net.helio.attrs.MaxRecords` and `~sunpy.net.helio.attrs.TableName` added as "attrs" for HELIO searches. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- Add the ability to download new GOES 16 & 17 data alongside the reprocessed GOES 13, 14 and 15 data via the GOES-XRS Fido client. (`#4394 <https://github.com/sunpy/sunpy/pull/4394>`__)
- `sunpy.net.jsoc.JSOCClient.request_data` now support additional parameter "method" which allows user to download staged data as single .tar file. (`#4405 <https://github.com/sunpy/sunpy/pull/4405>`__)
- Added ``sunpy.util.get_timerange_from_exdict`` which finds time range for a URL using its metadata.
  Added ``sunpy.util.scraper.Scraper.isvalid_time`` that checks whether the file corresponds to a desired time range. (`#4419 <https://github.com/sunpy/sunpy/pull/4419>`__)
- Colormap data has been moved to individual .csv files in the
  :file:`sunpy/visualization/colormaps/data` directory. (`#4433 <https://github.com/sunpy/sunpy/pull/4433>`__)
- Added `~sunpy.coordinates.utils.solar_angle_equivalency` to convert between a physical distance on the Sun (e.g., km) to an angular separation as seen by an observer (e.g., arcsec). (`#4443 <https://github.com/sunpy/sunpy/pull/4443>`__)
- `sunpy.map.Map` instances now have their ``.unit`` attribute set from the
  ``'BUNIT'`` FITS keyword. If the keyword cannot be parsed, or is not present
  the unit is set to `None`. (`#4451 <https://github.com/sunpy/sunpy/pull/4451>`__)
- The `sunpy.map.GenericMap.wcs` property is now cached, and will be recomputed
  only if changes are made to the map metadata. This improves performance of a
  number of places in the code base, and only one warning will now be raised
  about WCS fixes for a given set of metadata (as opposed to a warning each time
  ``.wcs`` is accessed) (`#4467 <https://github.com/sunpy/sunpy/pull/4467>`__)
- Extended :meth:`~sunpy.timeseries.GenericTimeSeries.concatenate` and
  :meth:`~sunpy.timeseries.TimeSeriesMetaData.concatenate` to allow iterables. (`#4499 <https://github.com/sunpy/sunpy/pull/4499>`__)
- Enable `~sunpy.coordinates.metaframes.RotatedSunFrame` to work with non-SunPy frames (e.g., `~astropy.coordinates.HeliocentricMeanEcliptic`). (`#4577 <https://github.com/sunpy/sunpy/pull/4577>`__)
- Add support for `pathlib.Path` objects to be passed to `sunpy.timeseries.TimeSeries`. (`#4589 <https://github.com/sunpy/sunpy/pull/4589>`__)
- Add support for GOES XRS netcdf files to be read as a `sunpy.timeseries.sources.XRSTimeSeries`. (`#4592 <https://github.com/sunpy/sunpy/pull/4592>`__)
- Add `~sunpy.net.jsoc.attrs.Cutout` attr for requesting cutouts
  from JSOC via `~sunpy.net.jsoc.JSOCClient` and ``Fido``. (`#4595 <https://github.com/sunpy/sunpy/pull/4595>`__)
- sunpy now sets auxillary parameters on `sunpy.map.GenericMap.wcs` using the
  `astropy.wcs.Wcsprm.aux` attribute. This stores observer information, along with
  the reference solar radius if present. (`#4620 <https://github.com/sunpy/sunpy/pull/4620>`__)
- The `~sunpy.coordinates.frames.HeliographicCarrington` frame now accepts the specification of ``observer='self'`` to indicate that the coordinate itself is also the observer for the coordinate frame.
  This functionality greatly simplifies working with locations of observatories that are provided in Carrington coordinates. (`#4659 <https://github.com/sunpy/sunpy/pull/4659>`__)
- Add two new colormaps (``rhessi`` and ``std_gamma_2``) that are used for plotting RHESSI maps. (`#4665 <https://github.com/sunpy/sunpy/pull/4665>`__)
- If either 'CTYPE1' or 'CTYPE2' are not present in map metadata, sunpy now assumes
  they are 'HPLN-TAN' and 'HPLT-TAN' (previously it assumed 'HPLN-   ' and 'HPLT-   ').
  In addition, a warning is also now raised when this assumption is made. (`#4702 <https://github.com/sunpy/sunpy/pull/4702>`__)
- Added a new `~sunpy.map.all_corner_coords_from_map` function to get the
  coordinates of all the pixel corners in a `~sunpy.map.GenericMap`. (`#4776 <https://github.com/sunpy/sunpy/pull/4776>`__)
- Added support for "%Y/%m/%dT%H:%M" to :func:`sunpy.time.parse_time`. (`#4791 <https://github.com/sunpy/sunpy/pull/4791>`__)
- Added the STEREO EUVI instrument specific colormaps called" 'euvi171', 'euvi195', 'euvi284', 'euvi304'. (`#4822 <https://github.com/sunpy/sunpy/pull/4822>`__)


Bug Fixes
---------

- `sunpy.map.GenericMap.date` now has its time scale set from the 'TIMESYS' FITS keyword,
  if it is present. If it isn't present the time scale defaults to 'UTC', which is unchanged
  default behaviour, so this change will only affect maps with a 'TIMESYS' keyword
  that is not set to 'UTC'. (`#4881 <https://github.com/sunpy/sunpy/pull/4881>`__)
- Fixed the `sunpy.net.dataretriever.sources.noaa.SRSClient` which silently failed to download the SRS files when the tarball for the previous years did not exist.
  Client now actually searches for the tarballs and srs files on the ftp archive before returning them as results. (`#4904 <https://github.com/sunpy/sunpy/pull/4904>`__)
- No longer is the WAVEUNIT keyword injected into a data source if it is missing from the file's metadata. (`#4926 <https://github.com/sunpy/sunpy/pull/4926>`__)
- Map sources no longer overwrite FITS metadata keywords if they are present in
  the original metadata. The particular map sources that have been fixed are
  `~sunpy.map.sources.SJIMap`, `~sunpy.map.sources.KCorMap`, `~sunpy.map.sources.RHESSIMap`,
  `~sunpy.map.sources.EITMap`, `~sunpy.map.sources.EUVIMap`, `~sunpy.map.sources.SXTMap`. (`#4926 <https://github.com/sunpy/sunpy/pull/4926>`__)
- Fixed a handling bug in ``sunpy.map.GenericMap.draw_rectangle`` when the rectangle is specified in a different coordinate frame than that of the map.
  A couple of other minor bugs in ``sunpy.map.GenericMap.draw_rectangle`` were also fixed. (`#4929 <https://github.com/sunpy/sunpy/pull/4929>`__)
- Improved error message from ``sunpy.net.Fido.fetch()`` when no email has been supplied for JSOC data. (`#4950 <https://github.com/sunpy/sunpy/pull/4950>`__)
- Fixed a bug when transforming from `~sunpy.coordinates.metaframes.RotatedSunFrame` to another frame at a different observation time that resulted in small inaccuracies.
  The translational motion of the Sun was not being handled correctly. (`#4979 <https://github.com/sunpy/sunpy/pull/4979>`__)
- Fixed two bugs with :func:`~sunpy.physics.differential_rotation.differential_rotate` and :func:`~sunpy.physics.differential_rotation.solar_rotate_coordinate` that resulted in significant inaccuracies.
  Both functions now ignore the translational motion of the Sun. (`#4979 <https://github.com/sunpy/sunpy/pull/4979>`__)
- The ability to to filter search results from the `~sunpy.net.vso.VSOClient` was broken.
  This has now been restored. (`#4011 <https://github.com/sunpy/sunpy/pull/4011>`__)
- Fixed a bug where transformation errors were not getting raised in some situations when a coordinate frame had ``obstime`` set to the default value of ``None`` and `~astropy.coordinates.SkyCoord` was not being used.
  Users are recommended to use `~astropy.coordinates.SkyCoord` to manage coordinate transformations unless they have a specific reason not to. (`#4267 <https://github.com/sunpy/sunpy/pull/4267>`__)
- Fixed a bug in `~sunpy.net.dataretriever.sources.goes.XRSClient._get_url_for_timerange` which returned incorrect URLs
  because of not using ``**kwargs`` in the client's ``_get_overlap_urls()`` method. (`#4288 <https://github.com/sunpy/sunpy/pull/4288>`__)
- Data products from `~sunpy.net.dataretriever.NOAAIndicesClient` and
  `~sunpy.net.dataretriever.NOAAPredictClient` have been updated to download
  new JSON files. The old text files which the data used to come in no longer
  exist. The new JSON files for `~sunpy.net.dataretriever.NOAAIndicesClient`
  now do not have the following columns:
  - Geomagnetic Observed and Smoothed
  - Sunspot Numbers Ratio (RI/SW)

  Both `sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries` and
  `sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries` have been updated to
  support the new JSON files. Loading the old text files is still supported,
  but support for this will be removed in a future version of sunpy. (`#4340 <https://github.com/sunpy/sunpy/pull/4340>`__)
- Fixed a bug due to which ``sunpy.net.helio.parser.wsdl_retriever`` ignored previously discovered Taverna links. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- The flare class labels in GOES ``peek()`` plots are now drawn at the center of
  the flare classes. Previously they were (ambiguously) drawn on the boundaries. (`#4364 <https://github.com/sunpy/sunpy/pull/4364>`__)
- `sunpy.map.GenericMap.rsun_obs` no longer assumes the observer is at Earth if
   ``rsun_obs`` was not present in the map metadata. The sun-observer
   distance is now taken directly from the observer coordinate. If the observer
   coordinate is not present, this defaults to the Earth, retaining previous
   behaviour. (`#4375 <https://github.com/sunpy/sunpy/pull/4375>`__)
- Nanosecond precision is now retained when using `~sunpy.time.parse_time` with
  a `~pandas.Timestamp`. (`#4409 <https://github.com/sunpy/sunpy/pull/4409>`__)
- Fixed a bug where SunPy could not be successfully imported if the default text encoding of the running environment was unable to handle non-ASCII characters. (`#4422 <https://github.com/sunpy/sunpy/pull/4422>`__)
- `sunpy.net.dataretriever.sources.noaa.SRSClient` now correctly returns zero
  results for queries in the future or before 1996, which is when data is first
  available. (`#4432 <https://github.com/sunpy/sunpy/pull/4432>`__)
- Fixes issue where NAXISn is not updated after invoking :meth:`.GenericMap.resample` (`#4445 <https://github.com/sunpy/sunpy/pull/4445>`__)
- The floating point precision of input to `sunpy.image.transform.affine_transform`
  is now preserved. Previously all input was cast to `numpy.float64`, which could
  cause large increases in memory use for 32 bit data. (`#4452 <https://github.com/sunpy/sunpy/pull/4452>`__)
- Fixed :func:`~sunpy.image.transform.affine_transform` to scale images to [0, 1] before
  passing them to :func:`skimage.transform.warp` and later rescale them back. (`#4477 <https://github.com/sunpy/sunpy/pull/4477>`__)
- Several ``warnings.simplefilter('always', Warning)`` warning filters in
  `sunpy.timeseries` have been removed. (`#4511 <https://github.com/sunpy/sunpy/pull/4511>`__)
- All calculations of the angular radius of the Sun now use the same underlying code with the accurate calculation.
  The previous inaccuracy was a relative error of ~0.001% (0.01 arcseconds) for an observer at 1 AU, but could be as large as ~0.5% for Parker Solar Probe perihelia. (`#4524 <https://github.com/sunpy/sunpy/pull/4524>`__)
- Fixed an issue in :meth:`sunpy.time.TimeRange.get_dates` where the function would return the wrong number of days if less than 24 hours had passed (`#4529 <https://github.com/sunpy/sunpy/pull/4529>`__)
- Several functions in `sunpy.map` now properly check if the provided coordinate is in the expected `~sunpy.coordinates.frames.Helioprojective` frame. (`#4552 <https://github.com/sunpy/sunpy/pull/4552>`__)
- Fixes a bug which occurs in setting the ``ylims`` by ``sunpy.visualization.animator.line.LineAnimator`` when there are non-finite values in the data array to be animated. (`#4554 <https://github.com/sunpy/sunpy/pull/4554>`__)
- Clear rotation metadata for SOHO/LASCO Helioviewer JPEG2000 images, as they are already rotated correctly. (`#4561 <https://github.com/sunpy/sunpy/pull/4561>`__)
- The ``max_conn`` argument to ``Fido.fetch()`` is now correctly respected by
  the JSOC client. Previously the JSOC client would default to 4 connections no
  matter what the value passed to ``Fido.fetch()`` was. (`#4567 <https://github.com/sunpy/sunpy/pull/4567>`__)
- :func:`sunpy.time.parse_time` now correctly parses lists of time strings that
  have one of the built in sunpy time formats. (`#4590 <https://github.com/sunpy/sunpy/pull/4590>`__)
- Fixes the SRSClient to search for files of correct queried time and now allows a path keyword to be downloaded in fetch. (`#4600 <https://github.com/sunpy/sunpy/pull/4600>`__)
- Fixed ``sunpy.net.helio.parser.wsdl_retriever``, which previously
  ignored discovered Taverna links. (`#4601 <https://github.com/sunpy/sunpy/pull/4601>`__)
- The transformations between `~astropy.coordinates.HCRS` and `~sunpy.coordinates.frames.HeliographicStonyhurst` have been re-implemented to enable the proper transformations of velocities.
  All ephemeris functions (e.g., :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`) now return properly calculated velocities when ``include_velocity=True`` is specified. (`#4613 <https://github.com/sunpy/sunpy/pull/4613>`__)
- The maximum number of connections opened by the JSOC downloader has been reduced
  from 4 to 2. This should prevent downloads of large numbers of files crashing. (`#4624 <https://github.com/sunpy/sunpy/pull/4624>`__)
- Fixed a significant performance bug that affected all coordinate transformations.
  Transformations have been sped up by a factor a few. (`#4663 <https://github.com/sunpy/sunpy/pull/4663>`__)
- Fixed a bug with the mapping of a WCS header to a coordinate frame if the observer location is provided in Carrington coordinates. (`#4669 <https://github.com/sunpy/sunpy/pull/4669>`__)
- `sunpy.io.fits.header_to_fits` now excludes any keys that have associated NaN
  values, as these are not valid in a FITS header, and throws a warning if this
  happens. (`#4676 <https://github.com/sunpy/sunpy/pull/4676>`__)
- Fixed an assumption in `sunpy.map.GenericMap.pixel_to_world` that the first
  data axis is longitude, and the second is latitude. This will affect you if
  you are using data where the x/y axes are latitude/longitude, and now returns
  correct values in methods and properties that call ``pixel_to_world``,
  such as ``bottom_left_coord``, ``top_right_coord``, ``center``. (`#4700 <https://github.com/sunpy/sunpy/pull/4700>`__)
- Added a warning when a 2D `~sunpy.coordinates.frames.Helioprojective` coordinate is upgraded to a 3D coordinate and the number type is lower precision than the native Python float.
  This 2D->3D upgrade is performed internally when transforming a 2D `~sunpy.coordinates.frames.Helioprojective` coordinate to any other coordinate frame. (`#4724 <https://github.com/sunpy/sunpy/pull/4724>`__)
- All columns from a :meth:`sunpy.net.vso.vso.VSOClient.search` will now be shown. (`#4788 <https://github.com/sunpy/sunpy/pull/4788>`__)
- The search results object returned from ``Fido.search``
  (`~sunpy.net.fido_factory.UnifiedResponse`) now correcly counts all results in
  it's `~sunpy.net.fido_factory.UnifiedResponse.file_num` property. Note that
  because some ``Fido`` clients now return metadata only results, this is really
  the number of records and does not always correspond to the number of files
  that would be downloaded. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- Improved the file processing logic for EVE L0CS files, which may have fixed a
  bug where the first line of data was parsed incorrectly. (`#4805 <https://github.com/sunpy/sunpy/pull/4805>`__)
- Fixing the ``CROTA`` meta keyword in EUVI FITS to ``CROTAn`` standard. (`#4846 <https://github.com/sunpy/sunpy/pull/4846>`__)


Added/Improved Documentation
----------------------------

- Added a developer guide for writing a new ``Fido`` client. (`#4387 <https://github.com/sunpy/sunpy/pull/4387>`__)
- Added an example of how to use Matplotlib's axes range functionality when plotting a Map with WCSAxes. (`#4792 <https://github.com/sunpy/sunpy/pull/4792>`__)
- Add links to Thompson 2006 paper on solar coordinates to synoptic map example. (`#3549 <https://github.com/sunpy/sunpy/pull/3549>`__)
- Clarified the meaning of ``.bottom_left_coord`` and ``.top_right_coord`` in
  `sunpy.map.GenericMap`. (`#3706 <https://github.com/sunpy/sunpy/pull/3706>`__)
- Added a list of possible signatures to
  `sunpy.timeseries.metadata.TimeSeriesMetaData`. (`#3709 <https://github.com/sunpy/sunpy/pull/3709>`__)
- Added `sunpy.data.manager`, `sunpy.data.cache`, `sunpy.net.Fido`, `sunpy.map.Map`,
  and `sunpy.timeseries.TimeSeries` to the docs. (`#4098 <https://github.com/sunpy/sunpy/pull/4098>`__)
- Clarified spline option for `sunpy.map.GenericMap.resample`. (`#4136 <https://github.com/sunpy/sunpy/pull/4136>`__)
- Updated the gallery example :ref:`sphx_glr_generated_gallery_plotting_solar_cycle_example.py` to retrieve data using `~sunpy.net.Fido`. (`#4169 <https://github.com/sunpy/sunpy/pull/4169>`__)
- Fixed example usage of :func:`~sunpy.io.fits.read` to account for the fact that it returns a list
  of data-header pairs rather than the data-header pairs directly. (`#4183 <https://github.com/sunpy/sunpy/pull/4183>`__)
- Added example of how to create a `sunpy.map.GenericMap` from observations in RA-DEC coordinates. (`#4236 <https://github.com/sunpy/sunpy/pull/4236>`__)
- Added `sunpy.coordinates.SunPyBaseCoordinateFrame` and `sunpy.coordinates.BaseHeliographic` to the documentation. (`#4274 <https://github.com/sunpy/sunpy/pull/4274>`__)
- `sunpy.time.TimeRange` had a ``.__contains__`` method and this is now documented. (`#4372 <https://github.com/sunpy/sunpy/pull/4372>`__)
- Revamped sunpy pull request review developer documentation. (`#4378 <https://github.com/sunpy/sunpy/pull/4378>`__)
- Revamped sunpy installation documentation. (`#4378 <https://github.com/sunpy/sunpy/pull/4378>`__)
- Fixed broken documentation links in the guide. (`#4414 <https://github.com/sunpy/sunpy/pull/4414>`__)
- Fixed miscellaneous links in the API documentation. (`#4415 <https://github.com/sunpy/sunpy/pull/4415>`__)
- Added `sunpy.data.data_manager.downloader`, `sunpy.data.data_manager.storage`,
  and `sunpy.net.hek.HEKTable` to the docs. (`#4418 <https://github.com/sunpy/sunpy/pull/4418>`__)
- Added documentation for copying Map objects using the copy module's deepcopy method. (`#4470 <https://github.com/sunpy/sunpy/pull/4470>`__)
- Added information about the :meth:`~sunpy.map.MapSequence.plot` return type. (`#4472 <https://github.com/sunpy/sunpy/pull/4472>`__)
- Added a gallery example for saving and loading sunpy Maps using asdf. (`#4494 <https://github.com/sunpy/sunpy/pull/4494>`__)
- Added description for a counter-intuitive section in the :ref:`sphx_glr_generated_gallery_differential_rotation_reprojected_map.py` example. (`#4548 <https://github.com/sunpy/sunpy/pull/4548>`__)
- Added :ref:`sunpy-coordinates-velocities` to explain how to use velocity information in the coordinates framework. (`#4610 <https://github.com/sunpy/sunpy/pull/4610>`__)
- New gallery example of searching and downloading GOES XRS data (with GOES 15, 16 and 17). (`#4686 <https://github.com/sunpy/sunpy/pull/4686>`__)
- Created the new gallery example :ref:`sphx_glr_generated_gallery_units_and_coordinates_north_offset_frame.py` for `~sunpy.coordinates.NorthOffsetFrame`. (`#4709 <https://github.com/sunpy/sunpy/pull/4709>`__)
- Added more information on which FITS keywords are used for various `sunpy.map.GenericMap`
  properties. (`#4717 <https://github.com/sunpy/sunpy/pull/4717>`__)
- Improved documentation for :func:`sunpy.physics.differential_rotation.diff_rot`. (`#4876 <https://github.com/sunpy/sunpy/pull/4876>`__)


Documentation Fixes
-------------------

- The keyword ``clip_interval`` is now used more extensively in gallery examples when plotting the sample AIA image (e.g., :ref:`sphx_glr_generated_gallery_plotting_aia_example.py`). (`#4573 <https://github.com/sunpy/sunpy/pull/4573>`__)
- Modified :ref:`sphx_glr_generated_gallery_plotting_magnetogram_active_regions.py` to use HMI file from sample data instead of downloading it with Fido. (`#4598 <https://github.com/sunpy/sunpy/pull/4598>`__)
- Removed unnecessary transformations of coordinates prior to plotting them using `~astropy.visualization.wcsaxes.WCSAxes.plot_coord`. (`#4609 <https://github.com/sunpy/sunpy/pull/4609>`__)
- Ensure that all attrs are documented and clean the `sunpy.net.hek.attrs`
  namespace of non-attr objects. (`#4834 <https://github.com/sunpy/sunpy/pull/4834>`__)
- Fixed miscellaneous issues with the gallery example :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_align_aia_hmi.py`. (`#4843 <https://github.com/sunpy/sunpy/pull/4843>`__)
- Fixed the display of arguments in the documentation for `~sunpy.net.Fido` attributes (`sunpy.net.attrs`). (`#4916 <https://github.com/sunpy/sunpy/pull/4916>`__)


Trivial/Internal Changes
------------------------

- ``Fido.fetch`` now always specifies a ``path=`` argument of type `pathlib.Path`
  to the ``fetch`` method of the client. This path will default to the configured
  sunpy download dir, will have the user directory expanded, will have the
  ``{file}`` placeholder and will be tested to ensure that it is writeable. (`#4949 <https://github.com/sunpy/sunpy/pull/4949>`__)
- Added information on what went wrong when `sunpy.map.GenericMap.wcs` fails to parse
  a FITS header into a WCS. (`#4335 <https://github.com/sunpy/sunpy/pull/4335>`__)
- Fixed the `~sunpy.coordinates.frames.Helioprojective` docstring to be clear about the names of the coordinate components. (`#4351 <https://github.com/sunpy/sunpy/pull/4351>`__)
- Raise a better error message if trying to load a FITS file that contains only
  one dimensional data. (`#4426 <https://github.com/sunpy/sunpy/pull/4426>`__)
- The following functions in `sunpy.map` have had their performance greatly increased,
  with runtimes typically improving by a factor of 20x. This has been achieved by
  improving many of the checks so that they only require checking the edge pixels of a
  map as opposed to all of the pixels.

  - :func:`~sunpy.map.contains_full_disk`
  - :func:`~sunpy.map.is_all_off_disk`
  - :func:`~sunpy.map.is_all_on_disk`
  - :func:`~sunpy.map.contains_limb` (`#4463 <https://github.com/sunpy/sunpy/pull/4463>`__)
- Improved the output when you print a sunpy Map. (`#4464 <https://github.com/sunpy/sunpy/pull/4464>`__)
- Creating a `~sunpy.util.MetaDict` with dictionary keys that are not strings now
  raises as user-friendly `ValueError` which prints all the non-compliant keys. (`#4476 <https://github.com/sunpy/sunpy/pull/4476>`__)
- Maps created directly via. `sunpy.map.GenericMap` now have their metadata
  automatically converted to a `~sunpy.util.MetaDict`, which is the same current
  behaviour of the `sunpy.map.Map` factory. (`#4476 <https://github.com/sunpy/sunpy/pull/4476>`__)
- If the ``top_right`` corner given to :meth:`sunpy.map.GenericMap.submap` is
  below or to the right of the ``bottom_left`` corner, a warning is no longer
  raised (as the rectangle is still well defined), but a message is still logged
  at the debug level to the sunpy logger. (`#4491 <https://github.com/sunpy/sunpy/pull/4491>`__)
- Added test support for Python 3.9 (no wheels yet). (`#4569 <https://github.com/sunpy/sunpy/pull/4569>`__)
- ``sunpy.sun`` functions now make use of the `~astropy.coordinates.GeocentricTrueEcliptic` frame to simplify internal calculations, but the returned values are unchanged. (`#4584 <https://github.com/sunpy/sunpy/pull/4584>`__)
- Change the format of the time returned from :func:`~sunpy.coordinates.sun.carrington_rotation_number`
  from ``'jd'`` to ``'iso'``, so printing the `~astropy.time.Time` returned will now print an ISO
  timestamp instead of the Julian days. (`#4819 <https://github.com/sunpy/sunpy/pull/4819>`__)
- The listings for the sample data (`sunpy.data.sample`) are now sorted. (`#4838 <https://github.com/sunpy/sunpy/pull/4838>`__)
- Changed the implementation of a ``hypothesis``-based test so that it does not raise an error with ``hypothesis`` 6.0.0. (`#4852 <https://github.com/sunpy/sunpy/pull/4852>`__)


2.0.0 (2020-06-12)
==================

Backwards Incompatible Changes
------------------------------

- The frames `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` now inherit from the new base class `~sunpy.coordinates.frames.BaseHeliographic`.
  This changes means that ``isinstance(frame, HeliographicStonyhurst)`` is no longer ``True`` when ``frame`` is `~sunpy.coordinates.frames.HeliographicCarrington`. (`#3595 <https://github.com/sunpy/sunpy/pull/3595>`__)
- `~sunpy.visualization.colormaps.color_tables.aia_color_table`, `~sunpy.visualization.colormaps.color_tables.eit_color_table` and `~sunpy.visualization.colormaps.color_tables.suvi_color_table` now only take `astropy.units` quantities instead of strings. (`#3640 <https://github.com/sunpy/sunpy/pull/3640>`__)
- `sunpy.map.Map` is now more strict when the metadata of a map cannot be validated, and
  an error is now thrown instead of a warning if metadata cannot be validated. In order to
  load maps that previously loaded without error you may need to pass ``silence_errors=True``
  to `sunpy.map.Map`. (`#3646 <https://github.com/sunpy/sunpy/pull/3646>`__)
- ``Fido.search`` will now return results from all clients which match a query, you no longer have to make the query specific to a single client. This means that searches involving the 'eve' and 'rhessi' instruments will now potentially also return results from the VSO. For `~sunpy.net.dataretriever.RHESSIClient` you can now specify ``a.Physobs("summary_lightcurve")`` to only include the summary lightcurve data products not provided by the VSO. (`#3770 <https://github.com/sunpy/sunpy/pull/3770>`__)
- The objects returned by the ``search`` methods on ``VSOClient``, ``JSOCClient`` and ``GenericClient`` have been changed to be based on `sunpy.net.base_client.BaseQueryResponse`. This introduces a few subtle breaking changes for people using the client search methods directly (not ``Fido.search``), or people using ``sunpy.net.fido_factory.UnifiedResponse.get_response``. When slicing an instance of ``QueryResponse`` it will now return an instance of itself, ``QueryResponse.blocks`` can be used to access the underlying records. Also, the ``.client`` attribute of the response no longer has to be the instance of the class the search was made with, however, it often is. (`#3770 <https://github.com/sunpy/sunpy/pull/3770>`__)
- `~sunpy.coordinates.frames.HeliographicCarrington` is now an observer-based frame, where the observer location (specifically, the distance from the Sun) is used to account for light travel time when determining apparent Carrington longitudes.  Coordinate transformations using this frame now require an observer to be specified. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- To enable the precise co-alignment of solar images from different observatories, the calculation of Carrington coordinates now ignores the stellar-aberration correction due to observer motion.
  For an Earth observer, this change results in an increase in Carrington longitude of ~20 arcseconds.
  See :ref:`sunpy-coordinates-carrington` for more information. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Fixed a bug where some of the coordinate transformations could raise `ValueError` instead of `~astropy.coordinates.ConvertError` when the transformation could not be performed. (`#3894 <https://github.com/sunpy/sunpy/pull/3894>`__)
- Astropy 3.2 is now the minimum required version of that dependency. (`#3936 <https://github.com/sunpy/sunpy/pull/3936>`__)


Deprecations and Removals
-------------------------

- Fido search attrs available as `sunpy.net.attrs` i.e, ``a.Time``, ``a.Instrument`` etc are now deprecated as VSO attrs (`sunpy.net.vso.attrs`). (`#3714 <https://github.com/sunpy/sunpy/pull/3714>`__)
- ``sunpy.util.multimethod.MultiMethod`` is deprecated, `functools.singledispatch` provides equivalent functionality in the standard library. (`#3714 <https://github.com/sunpy/sunpy/pull/3714>`__)
- ``sunpy.net.vso.attrs.Physobs`` has been moved to `sunpy.net.attrs.Physobs` and the original deprecated. (`#3877 <https://github.com/sunpy/sunpy/pull/3877>`__)
- Deprecate ``sunpy.instr.aia.aiaprep`` in favor of the `aiapy.calibrate.register` function in the
  [aiapy](https://gitlab.com/LMSAL_HUB/aia_hub/aiapy) package.
  ``sunpy.instr.aia.aiaprep`` will be removed in version 2.1 (`#3960 <https://github.com/sunpy/sunpy/pull/3960>`__)
- Removed the module ``sunpy.sun.sun``, which was deprecated in version 1.0.
  Use the module `sunpy.coordinates.sun` instead. (`#4014 <https://github.com/sunpy/sunpy/pull/4014>`__)
- Removed Sun-associated functions in `sunpy.coordinates.ephemeris`, which were deprecated in 1.0.
  Use the corresponding functions in `sunpy.coordinates.sun`. (`#4014 <https://github.com/sunpy/sunpy/pull/4014>`__)
- Remove the deprecated ``sunpy.net.vso.vso.VSOClient`` ``.query_legacy`` and ``.latest`` methods. (`#4109 <https://github.com/sunpy/sunpy/pull/4109>`__)
- Removed the sample datasets NOAAINDICES_TIMESERIES and NOAAPREDICT_TIMESERIES because they will invariably be out of date.
  Up-to-date versions of these NOAA indices can be downloaded using `~sunpy.net.Fido` (see :ref:`sphx_glr_generated_gallery_plotting_solar_cycle_example.py`). (`#4169 <https://github.com/sunpy/sunpy/pull/4169>`__)


Features
--------

- Added `~sunpy.coordinates.metaframes.RotatedSunFrame` for defining coordinate frames that account for solar rotation. (`#3537 <https://github.com/sunpy/sunpy/pull/3537>`__)
- Added a context manager (`~sunpy.coordinates.transform_with_sun_center`) to ignore the motion of the center of the Sun for coordinate transformations. (`#3540 <https://github.com/sunpy/sunpy/pull/3540>`__)
- Updated the gallery example titled 'Downloading and plotting an HMI magnetogram' to rotate the HMI magnetogram such that solar North is pointed up. (`#3573 <https://github.com/sunpy/sunpy/pull/3573>`__)
- Creates a function named ``sunpy.map.sample_at_coords`` that samples the data from the map at the given set of coordinates. (`#3592 <https://github.com/sunpy/sunpy/pull/3592>`__)
- Enabled the discovery of search attributes for each of our clients. (`#3637 <https://github.com/sunpy/sunpy/pull/3637>`__)
- Printing `sunpy.net.attrs.Instrument` or other "attrs" will show all attributes that exist under the corresponding "attr". (`#3637 <https://github.com/sunpy/sunpy/pull/3637>`__)
- Printing `sunpy.net.Fido` will print out all the clients that Fido can use. (`#3637 <https://github.com/sunpy/sunpy/pull/3637>`__)
- Updates `~sunpy.map.GenericMap.draw_grid` to allow disabling the axes labels and the ticks on the top and right axes. (`#3673 <https://github.com/sunpy/sunpy/pull/3673>`__)
- Creates a ``tables`` property for `~sunpy.net.fido_factory.UnifiedResponse`, which allows to access the `~sunpy.net.base_client.BaseQueryResponse` as an `~astropy.table.Table`, which then can be used for indexing of results. (`#3675 <https://github.com/sunpy/sunpy/pull/3675>`__)
- Change the APIs for ``sunpy.map.GenericMap.draw_rectangle`` and :meth:`sunpy.map.GenericMap.submap` to be consistent with each other and to use keyword-only arguments for specifying the bounding box. (`#3677 <https://github.com/sunpy/sunpy/pull/3677>`__)
- Updates the `~sunpy.map.GenericMap.observer_coordinate` property to warn the user of specific missing metadata for each frame.
  Omits warning about frames where all metadata is missing or all meta is present. (`#3692 <https://github.com/sunpy/sunpy/pull/3692>`__)
- Added `sunpy.util.config.copy_default_config` that copies the default config file to the user's config directory. (`#3722 <https://github.com/sunpy/sunpy/pull/3722>`__)
- ``sunpy.database`` now supports adding database entries and downloading data from ``HEK`` query (`#3731 <https://github.com/sunpy/sunpy/pull/3731>`__)
- Added a helper function (`~sunpy.coordinates.utils.get_rectangle_coordinates`) for defining a rectangle in longitude and latitude coordinates. (`#3737 <https://github.com/sunpy/sunpy/pull/3737>`__)
- Add a ``.data`` property in `~sunpy.timeseries.GenericTimeSeries`, so that users are encouraged to use :meth:`~sunpy.timeseries.GenericTimeSeries.to_dataframe` to get the data of the timeseries. (`#3746 <https://github.com/sunpy/sunpy/pull/3746>`__)
- It is now possible to turn on or off various corrections in :func:`~sunpy.coordinates.sun.L0` (the apparent Carrington longitude of Sun-disk center as seen from Earth). (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Made skimage.transform import lazy to reduce import time of `sunpy.image.transform` by ~50% (`#3818 <https://github.com/sunpy/sunpy/pull/3818>`__)
- Add support for parfive 1.1. This sets a limit on the number of open connections to JSOC when downloading files to 10. (`#3822 <https://github.com/sunpy/sunpy/pull/3822>`__)
- Fido clients (subclasses of `sunpy.net.base_client.BaseClient`) can now register their own attrs modules with `sunpy.net.attrs`.
  This allows clients which require attr classes specific to that client to register modules that can be used by the user i.e. ``a.vso``.
  It also allows clients implemented externally to sunpy to register attrs. (`#3869 <https://github.com/sunpy/sunpy/pull/3869>`__)
- Added the methods :meth:`sunpy.map.GenericMap.quicklook` and :meth:`sunpy.map.MapSequence.quicklook` to display an HTML summary of the instance, including interactive controls.
  When using Jupyter notebooks, this HTML summary is automatically shown instead of a text-only representation. (`#3951 <https://github.com/sunpy/sunpy/pull/3951>`__)
- Added `_localfilelist` method in ``sunpy.util.scraper.Scraper`` to scrap local data archives. (`#3994 <https://github.com/sunpy/sunpy/pull/3994>`__)
- Added extra constants to `sunpy.sun.constants`:

  - Longitude of the prime meridian (epoch J2000.0) : ``sunpy.sun.constants.get('W_0')``
  - Sidereal rotation rate : `sunpy.sun.constants.sidereal_rotation_rate`
  - First Carrington rotation (JD TT) : `sunpy.sun.constants.first_carrington_rotation`
  - Mean synodic period : `sunpy.sun.constants.mean_synodic_period`
  - Right ascension (RA) of the north pole (epoch J2000.0) : ``sunpy.sun.constants.get('alpha_0')``
  - Declination of the north pole (epoch J2000.0) : ``sunpy.sun.constants.get('delta_0')`` (`#4013 <https://github.com/sunpy/sunpy/pull/4013>`__)
- Adds to ``sunpy.util.scraper.Scraper`` the ability to include regular expressions in the URL passed. (`#4107 <https://github.com/sunpy/sunpy/pull/4107>`__)


Bug Fixes
---------

- Added support for passing ``TimeSeriesMetaData`` object to ``timeseries_factory`` and associated validation tests. (`#3639 <https://github.com/sunpy/sunpy/pull/3639>`__)
- Now when `~sunpy.map.GenericMap` fails to load a file, the filename that failed to load will now be part of the error message. (`#3727 <https://github.com/sunpy/sunpy/pull/3727>`__)
- Work around incorrect Content-Disposition headers in some VSO downloads, which were leading to mangled filenames. (`#3740 <https://github.com/sunpy/sunpy/pull/3740>`__)
- ``Fido.search`` can now service queries without ``a.Time`` being specified. This is currently only used by the `sunpy.net.jsoc.JSOCClient`. (`#3770 <https://github.com/sunpy/sunpy/pull/3770>`__)
- Fixed a bug with the calculation of Carrington longitude as seen from Earth where it was using an old approach instead of the current approach (for example, the varying Sun-Earth distance is now taken into account).
  The old approach resulted in errors no greater than 7 arcseconds in Carrington longitude when using `~sunpy.coordinates.sun.L0` and `~sunpy.coordinates.frames.HeliographicCarrington`. (`#3772 <https://github.com/sunpy/sunpy/pull/3772>`__)
- Updated `sunpy.map.CompositeMap.plot` to support a linewidths argument. (`#3792 <https://github.com/sunpy/sunpy/pull/3792>`__)
- Fix a bug in `sunpy.net.jsoc.JSOCClient` where requesting data for export would not work if a non-time primekey was used. (`#3825 <https://github.com/sunpy/sunpy/pull/3825>`__)
- Add support for passing paths of type `pathlib.Path` in `sunpy.net.jsoc.JSOCClient.fetch`. (`#3838 <https://github.com/sunpy/sunpy/pull/3838>`__)
- Add explicit support for dealing with download urls for files, under 'as-is' protocol in `sunpy.net.jsoc.JSOCClient.get_request`. (`#3838 <https://github.com/sunpy/sunpy/pull/3838>`__)
- Updated the method used to filter time in the VSO post-search filtering function. (`#3840 <https://github.com/sunpy/sunpy/pull/3840>`__)
- Fix failing of fetching of the indexed JSOCResponses using `sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`. (`#3852 <https://github.com/sunpy/sunpy/pull/3852>`__)
- Prevented `sunpy.map.GenericMap.plot` modifying in-place any items passed as ``imshow_kwargs``. (`#3867 <https://github.com/sunpy/sunpy/pull/3867>`__)
- Changed the format of DATE-OBS in `sunpy.map.GenericMap.wcs` from iso to isot (ie. with a "T" between the date and time) to conform with the FITS standard. (`#3872 <https://github.com/sunpy/sunpy/pull/3872>`__)
- Fixed a minor error (up to ~10 arcseconds) in the calculation of the Sun's position angle (:func:`sunpy.coordinates.sun.P`). (`#3886 <https://github.com/sunpy/sunpy/pull/3886>`__)
- `~sunpy.net.hek.HEKClient` was returning HTML and not JSON. (`#3899 <https://github.com/sunpy/sunpy/pull/3899>`__)
- Updated to HTTPS for HEK. (`#3917 <https://github.com/sunpy/sunpy/pull/3917>`__)
- The accuracy of the output of :func:`sunpy.coordinates.ephemeris.get_horizons_coord` is significantly improved. (`#3919 <https://github.com/sunpy/sunpy/pull/3919>`__)
- Fixed a bug where the longitude value for the reference coordinate in the Map repr would be displayed with the unintended longitude wrapping. (`#3959 <https://github.com/sunpy/sunpy/pull/3959>`__)
- It is now possible to specify a local file path to
  `sunpy.data.data_manager.DataManager.override_file` without having to prefix it
  with ``file://``. (`#3970 <https://github.com/sunpy/sunpy/pull/3970>`__)
- Closed the session in the destructor of VSOClient thus solving the problem of socket being left open (`#3973 <https://github.com/sunpy/sunpy/pull/3973>`__)
- Fixed a bug of where results of VSO searches would have inconsistent ordering in ``sunpy.net.vso.vso.QueryResponse`` by always sorting the results by start time. (`#3974 <https://github.com/sunpy/sunpy/pull/3974>`__)
- Fixes two bugs in `sunpy.util.deprecated`: correctly calculates the
  removal version and does not override the default and/or alternative functionality
  message. Providing a custom deprecation message now suppresses any
  mention of the removal version. Additionally, a ``pending`` keyword argument is
  provided to denote functions/classes that are pending deprecation. (`#3982 <https://github.com/sunpy/sunpy/pull/3982>`__)
- Correctly generate labels for sliders in
  ``~sunpy.visualization.animator.ArrayAnimatorWCS`` when the number of pixel
  dimensions and the number of world dimensions are not the same in the WCS. (`#3990 <https://github.com/sunpy/sunpy/pull/3990>`__)
- Updated VSOClient.response_block_properties to check if "None" is in the return. (`#3993 <https://github.com/sunpy/sunpy/pull/3993>`__)
- Fix a bug with ``sunpy.visualization.animator.ArrayAnimatorWCS`` where animating
  a line with a masked array with the whole of the initial line masked out the
  axes limits for the x axis were not correctly set. (`#4001 <https://github.com/sunpy/sunpy/pull/4001>`__)
- Fixed passing in a list of URLs into `sunpy.map.GenericMap`, before it caused an error due to the wrong type being returned. (`#4007 <https://github.com/sunpy/sunpy/pull/4007>`__)
- Fixed a bug with :func:`~sunpy.coordinates.transformations.transform_with_sun_center` where the global variable was sometimes restored incorrectly.
  This bug was most likely encountered if there was a nested use of this context manager. (`#4015 <https://github.com/sunpy/sunpy/pull/4015>`__)
- Fixes a bug in fido_factory to allow  path="./" in fido.fetch(). (`#4058 <https://github.com/sunpy/sunpy/pull/4058>`__)
- Prevented `sunpy.io.fits.header_to_fits` modifying the passed header in-place. (`#4067 <https://github.com/sunpy/sunpy/pull/4067>`__)
- Strip out any unknown unicode from the HEK response to prevent it failing to load some results. (`#4088 <https://github.com/sunpy/sunpy/pull/4088>`__)
- Fixed a bug in :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` that resulted in a error when requesting an array of locations in conjuction with enabling the light-travel-time correction. (`#4112 <https://github.com/sunpy/sunpy/pull/4112>`__)
- `sunpy.map.GenericMap.top_right_coord` and `~sunpy.map.GenericMap.center`
  have had their definitions clarified, and both have had off-by-one indexing
  errors fixed. (`#4121 <https://github.com/sunpy/sunpy/pull/4121>`__)
- Fixed `sunpy.map.GenericMap.submap()` when scaled pixel units (e.g. ``u.mpix``)
  are used. (`#4127 <https://github.com/sunpy/sunpy/pull/4127>`__)
- Fixed bugs in ``sunpy.util.scraper.Scraper.filelist``
  that resulted in error when the HTML page of URL opened by the scraper contains some "a" tags without "href" attribute
  and resulted in incorrect file urls when any href stores filepath relative to the URL's domain instead of just a filename. (`#4132 <https://github.com/sunpy/sunpy/pull/4132>`__)
- Fixed inconsistencies in how `~sunpy.map.GenericMap.submap` behaves when passed corners in pixel and world coordinates.
  The behavior for submaps specified in pixel coordinates is now well-defined for pixels on the boundary of the rectangle
  and is consistent for all boundaries. Previously pixels on the lower left boundary were included, but excluded on the
  upper and right boundary. This means the shape of a submap may now be 1 pixel larger in each dimension.
  Added several more tests for `~sunpy.map.GenericMap.submap` for a range of cutout sizes in both pixel and world
  coordinates. (`#4134 <https://github.com/sunpy/sunpy/pull/4134>`__)
- `sunpy.map.on_disk_bounding_coordinates` now fully propagates the coordinate
  frame of the input map to the output coordinates. Previously only the observer
  coordinate, and no other frame attributes, were propagated. (`#4141 <https://github.com/sunpy/sunpy/pull/4141>`__)
- Fix an off-by-one error in the reference pixel returned by
  `sunpy.map.make_fitswcs_header`. (`#4152 <https://github.com/sunpy/sunpy/pull/4152>`__)
- `sunpy.map.GenericMap.reference_pixel` now uses zero-based indexing, in order
  to be consistent with the rest of the `sunpy.map` API. (`#4154 <https://github.com/sunpy/sunpy/pull/4154>`__)
- Previously `sunpy.map.GenericMap.resample` with ``method='linear'`` was
  using an incorrect and constant value to fill edges when upsampling a map. Values
  near the edges are now correctly extrapolated using the ``fill_value=extrapolate``
  option to `scipy.interpolate.interp1d`. (`#4164 <https://github.com/sunpy/sunpy/pull/4164>`__)
- Fixed a bug where passing an `int` or `list` via the ``hdus`` keyword argument to
  `~sunpy.io.fits.read` threw an exception because the list of HDU objects was no longer
  of type `~astropy.io.fits.HDUList`. (`#4183 <https://github.com/sunpy/sunpy/pull/4183>`__)
- Fix attr printing when the attr registry is empty for that attr (`#4199 <https://github.com/sunpy/sunpy/pull/4199>`__)
- Improved the accuracy of :func:`~sunpy.coordinates.sun.angular_radius` by removing the use of the small-angle approximation.
  The inaccuracy had been less than 5 milliarcseconds. (`#4239 <https://github.com/sunpy/sunpy/pull/4239>`__)
- Fixed a bug with the ``observer`` frame attribute for coordinate frames where an input that was not supplied as a `~astropy.coordinates.SkyCoord` would sometimes result in a transformation error. (`#4266 <https://github.com/sunpy/sunpy/pull/4266>`__)


Improved Documentation
----------------------

- Fixed an issue with the scaling of class-inheritance diagrams in the online documentation by blocking the versions of graphviz containing a bug. (`#3548 <https://github.com/sunpy/sunpy/pull/3548>`__)
- A new example gallery example "Plotting a difference image" has been added,
  which can be used for base difference or running difference images. (`#3627 <https://github.com/sunpy/sunpy/pull/3627>`__)
- Removed obsolete Astropy Helpers submodule section in :file:`CONTRIBUTING.rst`;
  Also removed mentions of astropy_helpers in all files of the project. (`#3676 <https://github.com/sunpy/sunpy/pull/3676>`__)
- Corrected misleading `~sunpy.timeseries.metadata.TimeSeriesMetaData` documentation about optional parameters. (`#3680 <https://github.com/sunpy/sunpy/pull/3680>`__)
- Added an example for `~sunpy.map.GenericMap.world_to_pixel` function in the Units & Coordinates guide. (`#3776 <https://github.com/sunpy/sunpy/pull/3776>`__)
- Added a :ref:`page <sunpy-coordinates-carrington>` describing how SunPy calculates Carrington longitudes. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Changed padding value of an example in the example gallery to fix the overlap of titles and x-label axes. (`#3835 <https://github.com/sunpy/sunpy/pull/3835>`__)
- More information and links about how to create changelogs. (`#3856 <https://github.com/sunpy/sunpy/pull/3856>`__)
- Clarified some inputs to `sunpy.map.GenericMap.plot`. (`#3866 <https://github.com/sunpy/sunpy/pull/3866>`__)
- Changed quoted sentence (that we suggest authors add to their research papers) in CITATION.rst (`#3896 <https://github.com/sunpy/sunpy/pull/3896>`__)
- Add example of how to use SunPy's HEK client to search for the GOES flare event list. (`#3953 <https://github.com/sunpy/sunpy/pull/3953>`__)
- Improved the doc layout of `sunpy.data.sample`. (`#4034 <https://github.com/sunpy/sunpy/pull/4034>`__)
- Made improvements to STEREO starfield gallery example. (`#4039 <https://github.com/sunpy/sunpy/pull/4039>`__)
- Improved the documentation of `sunpy.map.GenericMap.resample`. (`#4043 <https://github.com/sunpy/sunpy/pull/4043>`__)
- Updated the STEREO starfield example to use all of the information in the star catalog. (`#4116 <https://github.com/sunpy/sunpy/pull/4116>`__)
- Mini-galleries are now easier to create in the documentation thanks to a custom Sphinx directive (``minigallery``).
  The page :ref:`sunpy-coordinates-rotatedsunframe` has an example of a mini-gallery at the bottom. (`#4124 <https://github.com/sunpy/sunpy/pull/4124>`__)
- Added `sunpy.visualization.colormaps.color_tables` to the docs. (`#4182 <https://github.com/sunpy/sunpy/pull/4182>`__)
- Made minor improvments to the map histogramming example. (`#4205 <https://github.com/sunpy/sunpy/pull/4205>`__)
- Add a warning to `sunpy.io` docs to recommend not using it for FITS (`#4208 <https://github.com/sunpy/sunpy/pull/4208>`__)


Trivial/Internal Changes
------------------------

- Removed un-used and un-tested code paths in the private ``_remove_lytaf_events`` function
  in ``sunpy.instr.lyra``. (`#3570 <https://github.com/sunpy/sunpy/pull/3570>`__)
- Removed ``astropy_helpers`` and this means that ``python setup.py <test,build_docs>`` no longer works.
  So if you want to:

  * Run the tests: Use ``tox -e <env name>`` or call ``pytest`` directly
  * Build the docs: Use ``tox -e docs`` or cd into the docs folder and run ``make html`` or ``sphinx-build docs docs/_build/html -W -b html -d docs/_build/.doctrees`` (`#3598 <https://github.com/sunpy/sunpy/pull/3598>`__)
- Cleaned up test warnings in sunpy.coordinates. (`#3652 <https://github.com/sunpy/sunpy/pull/3652>`__)
- Fix Python version for requiring importlib_resources (`#3683 <https://github.com/sunpy/sunpy/pull/3683>`__)
- `sunpy.net.attr.AttrWalker` no longer uses ``sunpy.util.multimethod.MultiMethod`` it uses a derivative of `functools.singledispatch` `sunpy.util.functools.seconddispatch` which dispatches on the second argument. (`#3714 <https://github.com/sunpy/sunpy/pull/3714>`__)
- Errors from a VSO search will now be raised to the user. (`#3719 <https://github.com/sunpy/sunpy/pull/3719>`__)
- Fixed the transformation test for `~sunpy.coordinates.metaframes.NorthOffsetFrame`, which would intermittently fail. (`#3775 <https://github.com/sunpy/sunpy/pull/3775>`__)
- :func:`~sunpy.coordinates.sun.earth_distance` is now computed without using coordinate transformations for better performance. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Created a helper function for testing the equality/closeness of longitude angles (i.e., angles with wrapping). (`#3804 <https://github.com/sunpy/sunpy/pull/3804>`__)
- Bump the astropy version figure tests are run with from 3.1.2 to 3.2.3 (`#3925 <https://github.com/sunpy/sunpy/pull/3925>`__)
- Used `urllib.parse.urlsplit` in ``sunpy.util.scraper`` for file scraping functionalities. (`#3956 <https://github.com/sunpy/sunpy/pull/3956>`__)
- Added `sunpy.net.base_client.BaseClient.check_attr_types_in_query` as a helper method
  to check if a query contains a set of required attributes, and is a subset of optional
  attributes. (`#3979 <https://github.com/sunpy/sunpy/pull/3979>`__)
- Removes appending login details for ftp urls from scraper. (`#4020 <https://github.com/sunpy/sunpy/pull/4020>`__)
- Re-factored the `sunpy.map.Map` factory to dispatch argument parsing based on type. (`#4037 <https://github.com/sunpy/sunpy/pull/4037>`__)
- Improved the error message raised by the Map factory when a map matches multiple source map types. (`#4052 <https://github.com/sunpy/sunpy/pull/4052>`__)
- Added log messages when the sample data fails to download. (`#4137 <https://github.com/sunpy/sunpy/pull/4137>`__)
- Remove an Astropy 3.1 comptibility wrapper for ``Quantity.to_string``. (`#4172 <https://github.com/sunpy/sunpy/pull/4172>`__)
- Refactor the sphinx config to no longer depend on astropy-sphinx and more
  closely match the new sunpy package template (`#4188 <https://github.com/sunpy/sunpy/pull/4188>`__)


1.1.0 (2020-01-10)
==================

Backwards Incompatible Changes
------------------------------

- The ``sunpy.net.vso.vso.get_online_vso_url`` function has been broken into two components, the new ``sunpy.net.vso.vso.get_online_vso_url`` function takes no arguments (it used to take three) and now only returns an online VSO mirror or None.
  The construction of a ``zeep.Client`` object is now handled by ``sunpy.net.vso.vso.build_client`` which has a more flexible API for customising the ``zeep.Client`` interface. (`#3330 <https://github.com/sunpy/sunpy/pull/3330>`__)
- Importing ``sunpy.timeseries.timeseriesbase`` no longer automatically imports
  Matplotlib. (`#3376 <https://github.com/sunpy/sunpy/pull/3376>`__)
- :meth:`sunpy.timeseries.sources.NOAAIndicesTimeSeries.peek()` now checks that the `type` argument is a
  valid string, and raises a `ValueError` if it isn't. (`#3378 <https://github.com/sunpy/sunpy/pull/3378>`__)
- Observer-based coordinate frames (`~sunpy.coordinates.frames.Heliocentric` and `~sunpy.coordinates.frames.Helioprojective`) no longer assume a default observer (Earth) if no observer is specified.  These frames can now be used with no observer specified, but most transformations cannot be performed for such frames.  This removal of a default observer only affects `sunpy.coordinates`, and has no impact on the default observer in `sunpy.map`. (`#3388 <https://github.com/sunpy/sunpy/pull/3388>`__)
- The callback functions provided to
  ``BaseFuncAnimator`` ``button_func`` keyword
  argument now take two positional arguments rather than one. The function
  signature is now ``(animator, event)`` where the first arg is the animator
  object, and the second is the matplotlib mouse event. (`#3407 <https://github.com/sunpy/sunpy/pull/3407>`__)
- The colormap stored in SunPy's Map subclasses (ie. ``map.plot_settings['cmap']``)
  can now be colormap string instead of the full `matplotlib.colors.Colormap`
  object. To get the full `~matplotlib.colors.Colormap` object use the new attribute
  ``map.cmap``. (`#3412 <https://github.com/sunpy/sunpy/pull/3412>`__)
- Fix a warning in `sunpy.map.GenericMap.rotate` where the truth value of an array
  was being calculated. This changes the behaviour of
  `~sunpy.map.GenericMap.rotate` when the ``angle=`` parameter is not an
  `~astropy.units.Quantity` object to raise `TypeError` rather than `ValueError`. (`#3456 <https://github.com/sunpy/sunpy/pull/3456>`__)


Deprecations and Removals
-------------------------

- Removed the step of reparing images (replacing non-finite entries with local mean) before coaligning them. The user is expected to do this themselves before coaligning images. If NaNs/non-finite entries are present, a warning is thrown.
  The function ``sunpy.image.coalignment.repair_image_nonfinite`` is deprecated. (`#3287 <https://github.com/sunpy/sunpy/pull/3287>`__)
- The method to convert a `~sunpy.coordinates.frames.Helioprojective` frame from 2D to 3D has been renamed from ``calculate_distance`` to `~sunpy.coordinates.frames.Helioprojective.make_3d`.  This method is not typically directly called by users. (`#3389 <https://github.com/sunpy/sunpy/pull/3389>`__)
- ``sunpy.visualization.animator.ImageAnimatorWCS`` is now deprecated in favour of
  ``ArrayAnimatorWCS``. (`#3407 <https://github.com/sunpy/sunpy/pull/3407>`__)
- ``sunpy.cm`` has been moved to `sunpy.visualization.colormaps` and will be
  removed in a future version. (`#3410 <https://github.com/sunpy/sunpy/pull/3410>`__)


Features
--------

- Add a new `sunpy.data.manager` and `sunpy.data.cache` for dealing with versioned remote data within functions.
  Please see the `Remote Data Manager <https://docs.sunpy.org/en/latest/dev_guide/remote_data.html>`__ guide. (`#3124 <https://github.com/sunpy/sunpy/pull/3124>`__)
- Added the coordinate frames `~sunpy.coordinates.frames.HeliocentricEarthEcliptic` (HEE), `~sunpy.coordinates.frames.GeocentricSolarEcliptic` (GSE), `~sunpy.coordinates.frames.HeliocentricInertial` (HCI), and `~sunpy.coordinates.frames.GeocentricEarthEquatorial` (GEI). (`#3212 <https://github.com/sunpy/sunpy/pull/3212>`__)
- Added SunPy Map support for GOES SUVI images. (`#3269 <https://github.com/sunpy/sunpy/pull/3269>`__)
- - Support APE14 for ``ImageAnimatorWCS`` in SunPy's visualization module (`#3275 <https://github.com/sunpy/sunpy/pull/3275>`__)
- Add ability to disable progressbars when dowloading files using `sunpy.net.helioviewer` and edited docstrings to mention this feature. (`#3280 <https://github.com/sunpy/sunpy/pull/3280>`__)
- Adds support for searching and downloading SUVI data. (`#3301 <https://github.com/sunpy/sunpy/pull/3301>`__)
- Log all VSO XML requests and responses to the SunPy logger at the ``DEBUG``
  level. (`#3330 <https://github.com/sunpy/sunpy/pull/3330>`__)
- Transformations between frames in `sunpy.coordinates` can now provide detailed debugging output.  Set the `logging` level to ``DEBUG`` to enable this output. (`#3339 <https://github.com/sunpy/sunpy/pull/3339>`__)
- Added the `sunpy.coordinates.sun.carrington_rotation_time` function to
  compute the time of a given Carrington rotation number. (`#3360 <https://github.com/sunpy/sunpy/pull/3360>`__)
- A new method has been added to remove columns from a
  `sunpy.timeseries.GenericTimeSeries`. (`#3361 <https://github.com/sunpy/sunpy/pull/3361>`__)
- Add ``shape`` property to TimeSeries. (`#3380 <https://github.com/sunpy/sunpy/pull/3380>`__)
- Added ASDF schemas for the new coordinate frames (`~sunpy.coordinates.frames.GeocentricEarthEquatorial`, `~sunpy.coordinates.frames.GeocentricSolarEcliptic`, `~sunpy.coordinates.frames.HeliocentricEarthEcliptic`, `~sunpy.coordinates.frames.HeliocentricInertial`).  See the gallery for an example of using ``asdf`` to save and load a coordinate frame. (`#3398 <https://github.com/sunpy/sunpy/pull/3398>`__)
- ``sunpy.visualization.animator.ArrayAnimatorWCS`` was added which uses the WCS
  object to get the coordinates of all axes, including the slider labels. It also provides the
  ability to customise the plot by specifying arguments to
  `~astropy.visualization.wcsaxes.WCSAxes` methods and supports animation of
  WCS aware line plots with Astroy 4.0. (`#3407 <https://github.com/sunpy/sunpy/pull/3407>`__)
- The returned list of `~sunpy.map.Map` objects is now sorted by filename when
  passing a directory or glob pattern to `~sunpy.map.map_factory.MapFactory`. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- Single character wildcards and character ranges can now be passed as
  glob patterns to `~sunpy.map.Map`. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- `~sunpy.map.Map` now accepts filenames and directories as `pathlib.Path`
  objects. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- `~sunpy.map.GenericMap` objects now have a ``.cmap`` attribute, which returns the full `~matplotlib.colors.Colormap`.
  object. (`#3412 <https://github.com/sunpy/sunpy/pull/3412>`__)
- `sunpy.io.write_file()` now accepts `~pathlib.Path` objects as filename inputs. (`#3469 <https://github.com/sunpy/sunpy/pull/3469>`__)
- `sunpy.map.make_fitswcs_header` now accepts a `tuple` representing the shape of an array as well as the actual array as the ``data`` argument. (`#3483 <https://github.com/sunpy/sunpy/pull/3483>`__)
- Made a couple of module imports lazy to reduce the import time of sunpy.map by
  ~40%. (`#3495 <https://github.com/sunpy/sunpy/pull/3495>`__)
- `sunpy.map.GenericMap.wcs` now uses the full FITS header to construct the WCS.
  This adds support for instruments with more complex projections, such as WISPR,
  however does mean that Map will be more sensitive to incorrect or invalid FITS
  headers. If you are using custom headers with SunPy Map you might encounter
  issues relating to this change. (`#3501 <https://github.com/sunpy/sunpy/pull/3501>`__)
- ``sunpy.visualization.animator.BaseFuncAnimator`` now takes an optional
  ``slider_labels`` keyword argument which draws text labels in the center of the
  sliders. (`#3504 <https://github.com/sunpy/sunpy/pull/3504>`__)
- Added a more helpful error message when trying to load a file or directory
  that doesn't exist with `sunpy.map.Map`. (`#3568 <https://github.com/sunpy/sunpy/pull/3568>`__)
- Add ``__repr__`` for `~sunpy.map.MapSequence` objects  so that users can view the
  critical information of all the ``Map`` objects, in a concise manner. (`#3636 <https://github.com/sunpy/sunpy/pull/3636>`__)


Bug Fixes
---------

- Fixed accuracy issues with the calculations of Carrington longitude (`~sunpy.coordinates.sun.L0`) and Carrington rotation number (`~sunpy.coordinates.sun.carrington_rotation_number`). (`#3178 <https://github.com/sunpy/sunpy/pull/3178>`__)
- Updated `sunpy.map.make_fitswcs_header` to be more strict on the inputs it accepts. (`#3183 <https://github.com/sunpy/sunpy/pull/3183>`__)
- Fix the calculation of ``rsun_ref`` in `~sunpy.map.make_fitswcs_header` and and
  ensure that the default reference pixel is indexed from 1. (`#3184 <https://github.com/sunpy/sunpy/pull/3184>`__)
- Fixed the missing transformation between two `~sunpy.coordinates.HeliographicCarrington` frames with different observation times. (`#3186 <https://github.com/sunpy/sunpy/pull/3186>`__)
- `sunpy.map.sources.AIAMap` and `sunpy.map.sources.HMIMap` will no longer assume
  the existance of certain header keys. (`#3217 <https://github.com/sunpy/sunpy/pull/3217>`__)
- `sunpy.map.make_fitswcs_header` now supports specifying the map projection
  rather than defaulting to ``TAN``. (`#3218 <https://github.com/sunpy/sunpy/pull/3218>`__)
- Fix the behaviour of
  ``sunpy.coordinates.frames.Helioprojective.calculate_distance`` if the
  representation isn't Spherical. (`#3219 <https://github.com/sunpy/sunpy/pull/3219>`__)
- Fixed a bug where the longitude of a coordinate would not wrap at the expected angle following a frame transformation. (`#3223 <https://github.com/sunpy/sunpy/pull/3223>`__)
- Fixed a bug where passing a time or time interval to the differential rotation function threw an error because the new observer was not in HGS. (`#3225 <https://github.com/sunpy/sunpy/pull/3225>`__)
- Fixed bug where `~sunpy.coordinates.ephemeris.get_horizons_coord` was unable to accept `~astropy.time.Time` arrays as input. (`#3227 <https://github.com/sunpy/sunpy/pull/3227>`__)
- Fix the ticks on the default heliographic grid overlay so they are not white
  (and normally invisible) by default. (`#3235 <https://github.com/sunpy/sunpy/pull/3235>`__)
- Fixed a bug with `sunpy.net.hek.HEKClient` when the results returned were a mixed dataset. (`#3240 <https://github.com/sunpy/sunpy/pull/3240>`__)
- Fix `sunpy.physics.differential_rotation.differential_rotate` to rotate in the
  correct direction and to account for the rotation of the heliographic
  coordinate frame with time. (`#3245 <https://github.com/sunpy/sunpy/pull/3245>`__)
- Fixed a bug with the handling of changing observation times for transformations between `~astropy.coordinates.HCRS` and `~sunpy.coordinates.frames.HeliographicStonyhurst`, which also indirectly affected other transformations when changing observation times. (`#3246 <https://github.com/sunpy/sunpy/pull/3246>`__)
- Fixed all coordinate transformations to properly handle a change in observation time. (`#3247 <https://github.com/sunpy/sunpy/pull/3247>`__)
- Fixed the handling of coordinates with velocity information when transforming between Astropy frames and SunPy frames. (`#3247 <https://github.com/sunpy/sunpy/pull/3247>`__)
- Fixed `~sunpy.physics.solar_rotation.calculate_solar_rotate_shift` so that it does not calculate a shift between the reference layer and itself, which would sometimes incorrectly result in a shift of a pixel due to numerical precision. (`#3255 <https://github.com/sunpy/sunpy/pull/3255>`__)
- Stop crash when ``LineAnimator`` ``axes_ranges`` entry given as ``1D`` array when data is ``>1D``, i.e. as an independent axis. (`#3283 <https://github.com/sunpy/sunpy/pull/3283>`__)
- Fixed a `sunpy.coordinates` bug where a frame using the default observer of Earth could have its observer overwritten during a transformation. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- Fixed a bug where the transformation from `~sunpy.coordinates.frames.Helioprojective` to `~sunpy.coordinates.frames.Heliocentric` used the Sun-observer distance from the wrong frame when shifting the origin, and thus might not give the correct answer if the observer was not the same for the two frames. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- Fixed a bug with the transformations between `~sunpy.coordinates.frames.Heliocentric` and `~sunpy.coordinates.frames.HeliographicStonyhurst` when the frame observation time was not the same as the observer observation time.  The most common way to encounter this bug was when transforming from `~sunpy.coordinates.frames.Helioprojective` to any non-observer-based frame while also changing the observation time. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- VSO client ``fetch`` should not download when ``wait`` keyword argument is specified. (`#3298 <https://github.com/sunpy/sunpy/pull/3298>`__)
- Fixed a bug with `~sunpy.coordinates.wcs_utils.solar_frame_to_wcs_mapping` that assumed that the supplied frame was a SunPy frame. (`#3305 <https://github.com/sunpy/sunpy/pull/3305>`__)
- Fixed bugs with `~sunpy.coordinates.wcs_utils.solar_frame_to_wcs_mapping` if the input frame does not include an observation time or an observer. (`#3305 <https://github.com/sunpy/sunpy/pull/3305>`__)
- `~sunpy.coordinates.utils.GreatArc` now accounts for the start and end points of the arc having different observers. (`#3334 <https://github.com/sunpy/sunpy/pull/3334>`__)
- Fixed situations where 2D coordinates provided to `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` were not converted to 3D as intended.  Furthermore, the stored data will always be the post-conversion, 3D version. (`#3351 <https://github.com/sunpy/sunpy/pull/3351>`__)
- Fix off by one error in `sunpy.map.make_fitswcs_header` where when using the
  default ``reference_pixel=None`` keyword argument the pixel coordinate of the
  reference pixel was off by +1. (`#3356 <https://github.com/sunpy/sunpy/pull/3356>`__)
- Updated both GOES XRS and LYRA dataretriever clients to use ``sunpy.util.scraper.Scraper``, to make sure that files are actually on the servers being queried. (`#3367 <https://github.com/sunpy/sunpy/pull/3367>`__)
- Fixing the ordering of lon and lat inputs into make_fitswcs_header (`#3371 <https://github.com/sunpy/sunpy/pull/3371>`__)
- Updated the URL for Fermi spacecraft-pointing files to use an HTTPS connection to HEASARC. (`#3381 <https://github.com/sunpy/sunpy/pull/3381>`__)
- Fixed a bug where permission denied errors when downloading files are very verbose by adding an error message in `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`. (`#3417 <https://github.com/sunpy/sunpy/pull/3417>`__)
- Fixed a malformed call to `astropy.time.Time` in a test, which resulted in an incorrect time scale (UTC instead of TT). (`#3418 <https://github.com/sunpy/sunpy/pull/3418>`__)
- Fix incorrect files being included in the tarball, and docs missing from the
  tarball (`#3423 <https://github.com/sunpy/sunpy/pull/3423>`__)
- Fixed a bug where clipping behavior had been enabled by default in the plotting normalizers for ``Map`` objects.  Clipping needs to be disabled to make use of the over/under/masked colors in the colormap. (`#3427 <https://github.com/sunpy/sunpy/pull/3427>`__)
- Fix a bug with observer based frames that prevented a coordinate with an array of obstimes being transformed to other frames. (`#3455 <https://github.com/sunpy/sunpy/pull/3455>`__)
- `sunpy.map.GenericMap` will no longer raise a warning if the posisition of the
  observer is not known for frames that don't need an observer, i.e. heliographic
  frames. (`#3462 <https://github.com/sunpy/sunpy/pull/3462>`__)
- Apply `os.path.expanduser` to `sunpy.map.map_factory.MapFactory` input
  before passing to `glob.glob` (`#3477 <https://github.com/sunpy/sunpy/pull/3477>`__)
- Fix multiple instances of `sunpy.map` sources assuming the type of FITS Header
  values. (`#3497 <https://github.com/sunpy/sunpy/pull/3497>`__)
- Fixed a bug with `~sunpy.coordinates.NorthOffsetFrame` where non-spherical representations for the north pole produced an error. (`#3517 <https://github.com/sunpy/sunpy/pull/3517>`__)
- Fixed ``map.__repr__`` when the coordinate system information contained in the
  ``CUNIT1/2`` metadata is not set to a known value. (`#3569 <https://github.com/sunpy/sunpy/pull/3569>`__)
- Fixed bugs with some coordinate transformations when ``obstime`` is ``None`` on the destination frame but can be assumed to be the same as the ``obstime`` of the source frame. (`#3576 <https://github.com/sunpy/sunpy/pull/3576>`__)
- Updated `sunpy.map.mapsequence.MapSequence` so that calling ``_derotate()`` raises ``NotImplementedError``.
  Added associated tests. (`#3613 <https://github.com/sunpy/sunpy/pull/3613>`__)
- Fixed pandas plotting registration in `sunpy.timeseries`. (`#3633 <https://github.com/sunpy/sunpy/pull/3633>`__)
- Correctly catch and emit a warning when converting a map metadata to a FITS
  header and it contains a keyword with non-ascii characters. (`#3645 <https://github.com/sunpy/sunpy/pull/3645>`__)


Improved Documentation
----------------------

- Clean up the docstring for `sunpy.physics.differential_rotation.solar_rotate_coordinate` to make the example clearer. (`#2708 <https://github.com/sunpy/sunpy/pull/2708>`__)
- Added new gallery examples and cleaned up various gallery examples. (`#3181 <https://github.com/sunpy/sunpy/pull/3181>`__)
- Cleaned and expanded upon the docstrings for each Fido Client. (`#3220 <https://github.com/sunpy/sunpy/pull/3220>`__)
- Added clarifying hyperlinks to the gallery example ``getting_lasco_observer_location`` to link to ``astroquery`` docs page. (`#3228 <https://github.com/sunpy/sunpy/pull/3228>`__)
- Added more details to docstrings in `sunpy.coordinates.frames`. (`#3262 <https://github.com/sunpy/sunpy/pull/3262>`__)
- Added a link to package maintainer list in the API Stability page. (`#3281 <https://github.com/sunpy/sunpy/pull/3281>`__)
- Improved the contributing guide by updating commands and highlighting text. (`#3394 <https://github.com/sunpy/sunpy/pull/3394>`__)
- Removing `.fits` from the end of path kwargs in `sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch` docs to change output file extension from ``{file}.fits.fits`` to ``{file}.fits``. (`#3399 <https://github.com/sunpy/sunpy/pull/3399>`__)
- A new example gallery section "Using SunPy with Other Packages" has been added,
  which contains a set of new examples using the `reproject
  <https://reproject.readthedocs.io/>`__ with solar data. (`#3405 <https://github.com/sunpy/sunpy/pull/3405>`__)
- Added a table of supported coordinate systems and other miscellaneous improvements to the :ref:`coordinates documentation <sunpy-coordinates>`. (`#3414 <https://github.com/sunpy/sunpy/pull/3414>`__)
- Clarified the meaning of :attr:`sunpy.map.GenericMap.dsun`. (`#3430 <https://github.com/sunpy/sunpy/pull/3430>`__)
- Fixed the plots with multiple subplots in the ``Map`` user guide to properly use `~astropy.visualization.wcsaxes` and to be appropriately sized. (`#3454 <https://github.com/sunpy/sunpy/pull/3454>`__)
- Fixed various issues with the gallery example of saving/loading coordinates using ``asdf``. (`#3473 <https://github.com/sunpy/sunpy/pull/3473>`__)
- Added ``sunpy.__citation__`` with a BibTex entry for citing sunpy. (`#3478 <https://github.com/sunpy/sunpy/pull/3478>`__)
- Added an example showing how to display two maps and fade between them. (`#3488 <https://github.com/sunpy/sunpy/pull/3488>`__)
- Clarified the meaning of some `~sunpy.map.GenericMap` observer properties. (`#3585 <https://github.com/sunpy/sunpy/pull/3585>`__)
- Added inherited members of `sunpy.map` classes to the docs. (`#3587 <https://github.com/sunpy/sunpy/pull/3587>`__)
- Fixed documentation of `sunpy.database.Database.search` by adding ``Returns`` docstring. (`#3593 <https://github.com/sunpy/sunpy/pull/3593>`__)
- Updated the docstring for the parameter ``sortby`` in `~sunpy.map.MapSequence` with the default value, valid value and how to disable sorting. (`#3601 <https://github.com/sunpy/sunpy/pull/3601>`__)
- Updated the tour guide to reflect that the time series is not random data. (`#3603 <https://github.com/sunpy/sunpy/pull/3603>`__)
- Fixes bold type and extra line breaks of remote data manager example. (`#3615 <https://github.com/sunpy/sunpy/pull/3615>`__)


Trivial/Internal Changes
------------------------

- Allow running our sphinx-gallery examples as Jupyter notebooks via Binder (`#3256 <https://github.com/sunpy/sunpy/pull/3256>`__)
- Improve error messages and type checking in
  ``sunpy.visualization.animator.image.ImageAnimatorWCS``. (`#3346 <https://github.com/sunpy/sunpy/pull/3346>`__)
- Copy the library ``distro`` into :file:`sunpy/extern`: replaces the deprecated ``platform/linux_distribution`` (`#3396 <https://github.com/sunpy/sunpy/pull/3396>`__)
- The version of Matplotlib used to generate figure tests has been bumped from
  3.0.3 to 3.1.1. (`#3406 <https://github.com/sunpy/sunpy/pull/3406>`__)
- Corrected spelling of 'plotting' in timeseries method (changed 'ploting' to 'plotting'). (`#3429 <https://github.com/sunpy/sunpy/pull/3429>`__)
- Switched to "importlib_metadata" to get package version to speed up import of SunPy. (`#3449 <https://github.com/sunpy/sunpy/pull/3449>`__)
- Fix tests for `sunpy.data.data_manager` and ensure they are correctly executed with pytest. (`#3550 <https://github.com/sunpy/sunpy/pull/3550>`__)


1.0.0 (2019-06-01)
==================

Backwards Incompatible Changes
------------------------------

- Move the matplotlib animators from ``sunpy.visualisation.imageanimator`` and
  ``sunpy.visualization.mapcubeanimator`` to `sunpy.visualization.animator`. (`#2515 <https://github.com/sunpy/sunpy/pull/2515>`__)
- Make `sunpy.time.parse_time` return `astropy.time.Time` instead of `datetime.datetime`. (`#2611 <https://github.com/sunpy/sunpy/pull/2611>`__)
- The properties and methods of `sunpy.time.TimeRange` returns `astropy.time.Time` and `astropy.time.TimeDelta` instead of `datetime.datetime` and `datetime.timedelta` respectively. (`#2638 <https://github.com/sunpy/sunpy/pull/2638>`__)
- The ``sunpy.instr.goes`` module now accepts and returns
  `sunpy.timeseries.sources.XRSTimeSeries` objects only. (`#2666 <https://github.com/sunpy/sunpy/pull/2666>`__)
- ``obstime`` keyword param of ``sunpy.instr.goes._goes_lx`` takes a non-scalar `astropy.time.Time` object instead of `numpy.ndarray`. The precision of times contained in `sunpy.timeseries` has been increased to 9 from 6. (`#2676 <https://github.com/sunpy/sunpy/pull/2676>`__)
- Removed ``sunpy.net.jsoc.attrs.Time`` because it served the same purpose as `sunpy.net.attrs.Time` after the switch to `astropy.time.Time`. (`#2694 <https://github.com/sunpy/sunpy/pull/2694>`__)
- Remove unused ``**kwargs`` within TimeSeries functions. (`#2717 <https://github.com/sunpy/sunpy/pull/2717>`__)
- Rotation matrices inside map objects were previously stored as numpy matrices, but are now
  stored as numpy arrays, as numpy will eventually remove their matrix datatype. See
  https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html for more information. (`#2719 <https://github.com/sunpy/sunpy/pull/2719>`__)
- The ``sunpy.cm.show_colormaps`` function now accepts the keyword 'search' instead of 'filter'. (`#2731 <https://github.com/sunpy/sunpy/pull/2731>`__)
- The keyword arguments to all client ``.fetch`` methods have been changed to
  support the new parfive downloader and to ensure consisteny across all Fido
  clients. (`#2797 <https://github.com/sunpy/sunpy/pull/2797>`__)
- The Helioviewer client has been switched to using the newer Helioviewer API.
  This has meant that we have changed some of the keywords that were passed into client's methods.
  We have enforced that several keywords (observatory,instrument,detector,measurement) need to be defined otherwise the functions cannot return any data. (`#2801 <https://github.com/sunpy/sunpy/pull/2801>`__)
- Maps no longer assume that the pixel units are arcseconds if the units aren't
  explicitly set. In addition to this if critical metadata is missing from when
  creating a map, the map will fail to initialize and will raise an error. (`#2847 <https://github.com/sunpy/sunpy/pull/2847>`__)
- axis_ranges kwarg of ``sunpy.visualization.animator.base.ArrayAnimator``, ``sunpy.visualization.animator.image.ImageAnimator`` and ``sunpy.visualization.animator.line.LineAnimator`` now must be entered as None, [min, max] or pixel edges of each array element. Previously, pixel centers were expected.  This change removes ambiguity in interpretation and ensures the extent of the plot can always be accurately derived. (`#2867 <https://github.com/sunpy/sunpy/pull/2867>`__)
- All keywords have been added (with defaults) to each `~sunpy.net.helioviewer.HelioviewerClient` function.
  This means that there will be some changes to the style of the PNG screenshot that is returned.
  Returns for the JPEG 2000 and the other functions should be the same but not guaranteed. (`#2883 <https://github.com/sunpy/sunpy/pull/2883>`__)
- Changed `sunpy.sun.models.interior` and `sunpy.sun.models.evolution` from `pandas.DataFrame` to `astropy.table.QTable` (`#2936 <https://github.com/sunpy/sunpy/pull/2936>`__)
- Minimum numpy version is now >=1.14.5 (`#2954 <https://github.com/sunpy/sunpy/pull/2954>`__)
- Removed ``sunpy.time.julian_day``, ``sunpy.time.julian_centuries``, ``sunpy.time.day_of_year``, ``sunpy.time.break_time``, ``sunpy.time.get_day``. (`#2999 <https://github.com/sunpy/sunpy/pull/2999>`__)
- Updated the solar values in `sunpy.sun.constants` to IAU 2015 values. (`#3001 <https://github.com/sunpy/sunpy/pull/3001>`__)
- Renamed ``eccentricity_sunearth_orbit`` to ``eccentricity_sun_earth_orbit``. (`#3001 <https://github.com/sunpy/sunpy/pull/3001>`__)
- Renamed ``sunpy.image.rescale`` to `sunpy.image.resample`. (`#3044 <https://github.com/sunpy/sunpy/pull/3044>`__)
- Remove the ``basic_plot`` keyword argument from
  `~sunpy.map.GenericMap.peek`. An example has been added to the gallery
  showing how to make a plot like this. (`#3109 <https://github.com/sunpy/sunpy/pull/3109>`__)
- `sunpy.map.GenericMap` will no longer use the key ``solar_b0`` as a value for heliographic latitude. (`#3115 <https://github.com/sunpy/sunpy/pull/3115>`__)
- `sunpy.map.GenericMap` now checks for a complete observer location rather than
  individually defaulting coordinates (lat, lon, distance) to Earth position. If
  any one of the three coordinates is missing from the header the observer will
  be defaulted to Earth and a warning raised. (`#3115 <https://github.com/sunpy/sunpy/pull/3115>`__)
- ``sunpy.sun.sun`` functions have been re-implemented using Astropy for significantly improved accuracy.  Some functions have been removed. (`#3137 <https://github.com/sunpy/sunpy/pull/3137>`__)
- All of the functions in ``sunpy.sun.sun`` and all of the Sun-specific functions in `sunpy.coordinates.ephemeris` have been moved to the new module `sunpy.coordinates.sun`. (`#3163 <https://github.com/sunpy/sunpy/pull/3163>`__)


Deprecations and Removals
-------------------------

- The deprecated ``sunpy.lightcurve``, ``sunpy.wcs`` and ``sunpy.spectra`` modules have now
  been removed. (`#2666 <https://github.com/sunpy/sunpy/pull/2666>`__)
- ``sunpy.instr.rhessi.get_obssumm_dbase_file`` ``sunpy.instr.rhessi.get_obssum_filename``, ``sunpy.instr.rhessi.get_obssumm_file`` have been removed. `~sunpy.net.Fido` should be used to download these files. (`#2808 <https://github.com/sunpy/sunpy/pull/2808>`__)
- Removed ``heliographic_solar_center`` in favour of ``sunpy.coordinates.get_sun_L0`` and ``sunpy.coordinates.get_sun_B0`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``GenericClient.query`` in favour of `sunpy.net.dataretriever.GenericClient.search` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunearth_distance`` in favour of ``get_sunearth_distance`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``remove_lytaf_events_from_lightcurve`` in favour of ``sunpy.instr.lyra.remove_lytaf_events_from_timeseries`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.cm.get_cmap`` in favour of ``plt.get_cmap`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``database.query`` in favour of `sunpy.database.Database.search` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.net.vso.InteractiveVSOClient`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``MapCube`` in favour of `~sunpy.map.MapSequence` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``solar_north`` in favour of ``get_sun_P`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``database.download`` in favour of `sunpy.database.Database.fetch` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.map.GenericMap.pixel_to_data`` in favour of `sunpy.map.GenericMap.pixel_to_world` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``GenericClient.get`` in favour of `sunpy.net.dataretriever.GenericClient.fetch`. This changes applies to the other clients as well. (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``Map.xrange`` and ``Map.yrange`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.net.attrs.Wave`` in favour of ``sunpy.net.vso.attrs.Wavelength`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``JSOCClient.check_request`` in favour of `drms.client.ExportRequest.status` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- ``sunpy.net.vso.VSOClient.query_legacy`` and ``sunpy.net.vso.VSOClient.latest`` have been deprecated as we strongly recommend people use `sunpy.net.Fido` for all queries. (`#2866 <https://github.com/sunpy/sunpy/pull/2866>`__)
- The deprecated ``sunpy.physics.transforms`` module has been removed, it is
  replaced by `sunpy.physics.solar_rotation` and
  `sunpy.physics.differential_rotation`. (`#2994 <https://github.com/sunpy/sunpy/pull/2994>`__)
- Removed ``sunpy.sun.sun.solar_cycle_number`` because it was fundamentally flawed (`#3150 <https://github.com/sunpy/sunpy/pull/3150>`__)


Features
--------

- Change arguments to ``sunpy.test`` from ``offline=`` and ``online=`` to ``online`` and ``online_only``. This matches the behavior of the figure keyword arguments and comes as a part of a move to using a modified version of the Astropy test runner. (`#1983 <https://github.com/sunpy/sunpy/pull/1983>`__)
- asdf schemas and tags were added for the SunPy coordinate frames and `~sunpy.map.GenericMap` allowing these objects to be saved to and restored from `asdf <https://asdf.readthedocs.io/>`__ files. (`#2366 <https://github.com/sunpy/sunpy/pull/2366>`__)
- The images from image tests are now saved in a local folder for easy access. (`#2507 <https://github.com/sunpy/sunpy/pull/2507>`__)
- ``sunpy.map.MapCube`` has been renamed to `sunpy.map.MapSequence` to better reflect its use as a collection of map objects. (`#2603 <https://github.com/sunpy/sunpy/pull/2603>`__)
- Net search attributes now support tab completion of values and display a table of possible values when printed, to allow easier discoverability of possible search values. (`#2663 <https://github.com/sunpy/sunpy/pull/2663>`__)
- Running the figure tests now creates a page showing the differences between
  the expected figures and the figures produced from running the tests. (`#2681 <https://github.com/sunpy/sunpy/pull/2681>`__)
- Add support for Dask arrays in `sunpy.map.Map`. The map factory now checks a whitelist
  of array types rather than strictly checking if the array is of type `numpy.ndarray`. (`#2689 <https://github.com/sunpy/sunpy/pull/2689>`__)
- Persist the name of a coordinate, i.e. "earth" even though a concrete
  coordinate object has been calculated and use this string representation to change
  the way the sunpy frames are printed. This is primarily to facilitate displaying
  the name of the body rather than the concrete coordinate when printing a
  `~astropy.coordinates.SkyCoord`. (`#2723 <https://github.com/sunpy/sunpy/pull/2723>`__)
- `~sunpy.net.hek.HEKClient.search` now returns an `astropy.table.Table` instead of list of a `dict`. (`#2759 <https://github.com/sunpy/sunpy/pull/2759>`__)
- Add a downscaled HMI image to the sample data. (`#2782 <https://github.com/sunpy/sunpy/pull/2782>`__)
- Now able to create a `sunpy.map.Map` using an array and a `astropy.wcs.WCS` object. (`#2793 <https://github.com/sunpy/sunpy/pull/2793>`__)
- The download manager for `~sunpy.net.Fido` has been replaced with
  `parfive <https://parfive.readthedocs.io/en/latest/>`__. This provides advanced
  progress bars, proper handling of overwriting and the ability to retry failed
  downloads. (`#2797 <https://github.com/sunpy/sunpy/pull/2797>`__)
- `sunpy.map.GenericMap` can now save out rice compressed FITS files. (`#2826 <https://github.com/sunpy/sunpy/pull/2826>`__)
- Now any SunPyDeprecationWarnings will cause an error when using pytest. (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Added full Tox support for SunPy tests, documentation build and figure tests. (`#2839 <https://github.com/sunpy/sunpy/pull/2839>`__)
- Transition the `sunpy.net.vso.VSOClient` from using suds to `zeep <https://python-zeep.readthedocs.io/en/master/>`__ as the SOAP
  library. This is a more actively maintained library, and should provide better
  support for the VSOs https endpoints. This change should have no effect on the
  public API of the `sunpy.net.vso.VSOClient`. (`#2866 <https://github.com/sunpy/sunpy/pull/2866>`__)
- Provided access to the Helioviewer header information using `~sunpy.net.helioviewer.HelioviewerClient.get_jp2_header` function. (`#2904 <https://github.com/sunpy/sunpy/pull/2904>`__)
- Add a new WSDL URL and port to support SunPy use of VSO instance at SDAC. (`#2912 <https://github.com/sunpy/sunpy/pull/2912>`__)
- Add support for COSMO K-Coronograph (KCOR) FITS data. (`#2916 <https://github.com/sunpy/sunpy/pull/2916>`__)
- Add logger messaging system based on `~astropy.logger.AstropyLogger`, cleaned up all warnings, removed all print statements. (`#2980 <https://github.com/sunpy/sunpy/pull/2980>`__)
- The function `sunpy.image.coalignment.get_correlation_shifts` now issues an error when the number of dimensions
  are not correct instead of a warning and returning None. (`#2980 <https://github.com/sunpy/sunpy/pull/2980>`__)
- The default location of the sunpy sample data has changed to be in the platform
  specific data directory as provided by `appdirs <https://github.com/ActiveState/appdirs>`__. (`#2993 <https://github.com/sunpy/sunpy/pull/2993>`__)
- Add timeseries support for EVE/ESP level 1 data in `sunpy.timeseries.sources` (`#3032 <https://github.com/sunpy/sunpy/pull/3032>`__)
- The default style for Map plots have changed to reflect the changes in Astropy
  3.2. (`#3054 <https://github.com/sunpy/sunpy/pull/3054>`__)
- `sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` can now account for light travel time when computing the (apparent) body position, as long as the observer location is provided. (`#3055 <https://github.com/sunpy/sunpy/pull/3055>`__)
- Added a helper function (`sunpy.map.make_fitswcs_header`) that allows users to create a meta header for custom created `sunpy.map.GenericMap`. (`#3083 <https://github.com/sunpy/sunpy/pull/3083>`__)
- Map plotting now accepts the optional keyword ``clip_interval`` for specifying a percentile interval for clipping.  For example, if the interval (5%, 99%) is specified, the bounds of the z axis are chosen such that the lowest 5% of pixels and the highest 1% of pixels are excluded. (`#3100 <https://github.com/sunpy/sunpy/pull/3100>`__)
- The new function `~sunpy.coordinates.get_horizons_coord` enables querying JPL HORIZONS for the locations of a wide range of solar-system bodies, including spacecraft. (`#3113 <https://github.com/sunpy/sunpy/pull/3113>`__)


Bug Fixes
---------

- Fix the bug that prevented VSO queries for HMI data from downloading file
  without specifying ``a.Physobs``. (`#2621 <https://github.com/sunpy/sunpy/pull/2621>`__)
- Fix ``sunpy.map.mapcube.MapCube.plot``. The code had not been updated to support the changes to the wcsaxes helper functions. (`#2627 <https://github.com/sunpy/sunpy/pull/2627>`__)
- Replace all use of the deprecated ``sunpy.cm.get_cmap`` with `matplotlib.cm.get_cmap` to prevent deprecation warnings being raised. (`#2635 <https://github.com/sunpy/sunpy/pull/2635>`__)
- Fix generation of the coordinate transformation graph with Astropy 3.1.dev (`#2636 <https://github.com/sunpy/sunpy/pull/2636>`__)
- Prevent helioviewer from erroring when downloading file to a directory that
  does not exist. It will now create the directory when required. (`#2642 <https://github.com/sunpy/sunpy/pull/2642>`__)
- Fix transformations into/out of Heliographic Stonyhurst frame when
  the coordinate representation is Cartesian. (`#2646 <https://github.com/sunpy/sunpy/pull/2646>`__)
- Running the figure tests with ``setup.py test`` now saves the figures and the hashes to the same directory as setup.py. (`#2658 <https://github.com/sunpy/sunpy/pull/2658>`__)
- ``sunpy.instr.fermi.met_to_utc`` now returns the correct utc time which takes into account the leap seconds that have passed. (`#2679 <https://github.com/sunpy/sunpy/pull/2679>`__)
- Support passing Python file objects to `sunpy.io.fits.write`. (`#2688 <https://github.com/sunpy/sunpy/pull/2688>`__)
- Added DRMS to setup.py so sunpy[all] installs it as a dependency. (`#2693 <https://github.com/sunpy/sunpy/pull/2693>`__)
- Fix eve 0cs timeseries seperator regex to support Python 3.7 (`#2697 <https://github.com/sunpy/sunpy/pull/2697>`__)
- Fix the bug which crashes `~sunpy.map.sources.LASCOMap` for when 'date-obs' is reformatted agian from a self applied function. (`#2700 <https://github.com/sunpy/sunpy/pull/2700>`__)
- Change all instances of quantity_allclose to `astropy.units.allclose` this prevents pytest being needed to import `sunpy.coordinates` on Astropy 3 (`#2701 <https://github.com/sunpy/sunpy/pull/2701>`__)
- Fix RHESSI obssum file downloading to include the final day in the time range. (`#2714 <https://github.com/sunpy/sunpy/pull/2714>`__)
- Raise an error when transforming between HPC and HCC frames if the observer is not the same. (`#2725 <https://github.com/sunpy/sunpy/pull/2725>`__)
- Replaces the existing LASCO C2 and C3 color maps with new ones that perform better with JP2 and Level 0.5, 1 data. (`#2731 <https://github.com/sunpy/sunpy/pull/2731>`__)
- Do not attempt to save a FITS header comment for a keyword which is not in the header. This prevents an error on saving some maps after the metadata had been modified but not the comments. (`#2748 <https://github.com/sunpy/sunpy/pull/2748>`__)
- Add support for `~sunpy.map.sources.HMIMap` objects as input to ``sunpy.instr.aia.aiaprep``. (`#2749 <https://github.com/sunpy/sunpy/pull/2749>`__)
- User can convert between HPC and HCC coordinates with different observers. This is implemented by automatically transforming the coordinate into HGS and then changing observer, and then transforming back to HCC. (`#2754 <https://github.com/sunpy/sunpy/pull/2754>`__)
- Changed default file type for Helioviewer to prevent decode errors. (`#2771 <https://github.com/sunpy/sunpy/pull/2771>`__)
- Increase figure size to avoid cutting off longer colormap names in ``sunpy.cm.show_colormaps``. (`#2824 <https://github.com/sunpy/sunpy/pull/2824>`__)
- The sample data directory will no longer be created until files are downloaded
  to it. (`#2836 <https://github.com/sunpy/sunpy/pull/2836>`__)
- Timeseries and lightcurve will now respect updated config values for download directory. (`#2844 <https://github.com/sunpy/sunpy/pull/2844>`__)
- Always use _default_wrap_angle rather than hard coding a wrap angle in the init
  of a sunpy coordinate frame (`#2853 <https://github.com/sunpy/sunpy/pull/2853>`__)
- Ensure imageanimators only slice arrays with integers (`#2856 <https://github.com/sunpy/sunpy/pull/2856>`__)
- Fixed `sunpy.io.fits.write` to handle the keyword ``COMMENT`` correctly. (`#2880 <https://github.com/sunpy/sunpy/pull/2880>`__)
- If Carrington longitude ("crln_obs") is found in the FITS header, `~sunpy.map.Map` converts this to the correct Heliographic longitude. (`#2946 <https://github.com/sunpy/sunpy/pull/2946>`__)
- ``sunpy.net.helio.hec.HECClient.time_query`` now resolves the correct input time format. (`#2969 <https://github.com/sunpy/sunpy/pull/2969>`__)
- Fixes the calculation of the solar rotation of coordinates and the differential rotation of `sunpy.map.GenericMap`. (`#2972 <https://github.com/sunpy/sunpy/pull/2972>`__)
- Added back the FERMI GBM client to `sunpy.net.dataretriever`. (`#2983 <https://github.com/sunpy/sunpy/pull/2983>`__)
- Fix bug in `sunpy.net.hek` which raised and error if a search returned zero results, now returns an empty `sunpy.net.hek.HEKTable`. (`#3046 <https://github.com/sunpy/sunpy/pull/3046>`__)
- `~sunpy.map.sources.AIAMap` now uses the provided HAE coordinates instead of the provided HGS coordinates to determine the observer location. (`#3056 <https://github.com/sunpy/sunpy/pull/3056>`__)
- Correctly zero pad milliseconds in the ``sunpy.util.scraper.Scraper`` formatting to prevent errors when the millisecond value was less than 100. (`#3063 <https://github.com/sunpy/sunpy/pull/3063>`__)
- Fix ``sunpy.util.scraper.Scraper`` failing if a directory is not found on a remote server. (`#3063 <https://github.com/sunpy/sunpy/pull/3063>`__)
- Correctly extract observer location from MDI and EIT data (`#3067 <https://github.com/sunpy/sunpy/pull/3067>`__)
- Fix HGS <> HCRS test due to Ecliptic frame changes in astropy 3.2 (`#3075 <https://github.com/sunpy/sunpy/pull/3075>`__)
- Fixes bug when creating a timeseries from a URL and bug when creating a TimeSeries from  older GOES/XRS fits files. (`#3081 <https://github.com/sunpy/sunpy/pull/3081>`__)
- Added `~sunpy.map.sources.EUVIMap.rsun_obs`. It returns a quantity in arcsec consistent with other `sunpy.map.GenericMap` and overwrites mapbase's assumption of a photospheric limb as seen from Earth. (`#3099 <https://github.com/sunpy/sunpy/pull/3099>`__)
- Fixed bugs related to using `~sunpy.map.GenericMap.plot` and `~sunpy.map.GenericMap.peek` with the ``inline`` Matplotlib backend in Jupyter notebook. (`#3103 <https://github.com/sunpy/sunpy/pull/3103>`__)
- Make a correction to `sunpy.coordinates.wcs_utils.solar_wcs_frame_mapping` so
  that `astropy.wcs.WCS` objects are correctly converted to
  `sunpy.coordinates.frames` objects irrespective of the ordering of the axes. (`#3116 <https://github.com/sunpy/sunpy/pull/3116>`__)
- The `~sunpy.physics.differential_rotation.solar_rotate_coordinate` function returns a coordinate that accounts for the location of the new observer. (`#3123 <https://github.com/sunpy/sunpy/pull/3123>`__)
- Add support for rotation parameters to `sunpy.map.make_fitswcs_header`. (`#3139 <https://github.com/sunpy/sunpy/pull/3139>`__)
- Improve the implementation of `~sunpy.physics.differential_rotation.differential_rotate` the image warping when transforming Maps for differential rotation and change in observer position. (`#3149 <https://github.com/sunpy/sunpy/pull/3149>`__)
- Fix a bug where new helioviewer sources potentially cause `~sunpy.net.helioviewer.HelioviewerClient.data_sources` to error. (`#3162 <https://github.com/sunpy/sunpy/pull/3162>`__)


Improved Documentation
----------------------

- Organise the gallery into sections based on example type and tidy up a little. (`#2624 <https://github.com/sunpy/sunpy/pull/2624>`__)
- Added gallery example showing the conversion of Helioprojective Coordinates to Altitude/Azimuth Coordinates to and back. (`#2656 <https://github.com/sunpy/sunpy/pull/2656>`__)
- Add contribution guidelines for the sunpy example gallery. (`#2682 <https://github.com/sunpy/sunpy/pull/2682>`__)
- Added a gallery example for "Downloading and plotting a HMI image" and "Creating a Composite map". (`#2746 <https://github.com/sunpy/sunpy/pull/2746>`__)
- Added an example for ``sunpy.visualization.animator.ImageAnimatorWCS``. (`#2752 <https://github.com/sunpy/sunpy/pull/2752>`__)
- Minor changes to the developer guide regarding sprint labels. (`#2765 <https://github.com/sunpy/sunpy/pull/2765>`__)
- Copyedited and corrected the solar cycles example. (`#2770 <https://github.com/sunpy/sunpy/pull/2770>`__)
- Changed "online" mark to "remote_data" and made formatting of marks consistent. (`#2799 <https://github.com/sunpy/sunpy/pull/2799>`__)
- Add a missing plot to the end of the units and coordinates guide. (`#2813 <https://github.com/sunpy/sunpy/pull/2813>`__)
- Added gallery example showing how to access the SunPy colormaps (`#2865 <https://github.com/sunpy/sunpy/pull/2865>`__)
- Added gallery example showing how to access the SunPy solar physics constants. (`#2882 <https://github.com/sunpy/sunpy/pull/2882>`__)
- Major clean up of the developer documentation. (`#2951 <https://github.com/sunpy/sunpy/pull/2951>`__)
- Overhaul of the install intructions for the guide section of our documentation. (`#3147 <https://github.com/sunpy/sunpy/pull/3147>`__)


Trivial/Internal Changes
------------------------

- `~sunpy.time.parse_time` now uses `~functools.singledispatch` underneath. (`#2408 <https://github.com/sunpy/sunpy/pull/2408>`__)
- Revert the handling of ``quantity_allclose`` now that `astropy/astropy#7252 <https://github.com/astropy/astropy/pull/7252>`__ is merged. This also bumps the minimum astropy version to 3.0.2. (`#2598 <https://github.com/sunpy/sunpy/pull/2598>`__)
- Replace the subclasses of matplotlib Slider and Button in `sunpy.visualization` with partial functions. (`#2613 <https://github.com/sunpy/sunpy/pull/2613>`__)
- Sort the ana C source files before building to enable reproducible builds. (`#2637 <https://github.com/sunpy/sunpy/pull/2637>`__)
- We are now using `towncrier <https://github.com/hawkowl/towncrier>`__ to
  generate our changelogs. (`#2644 <https://github.com/sunpy/sunpy/pull/2644>`__)
- Moved figure tests to Python 3.6. (`#2655 <https://github.com/sunpy/sunpy/pull/2655>`__)
- Removed old metaclass used for Map and TimeSeries as we have now moved to Python 3.6. (`#2655 <https://github.com/sunpy/sunpy/pull/2655>`__)
- Updated astropy_helpers to v3.0.2. (`#2655 <https://github.com/sunpy/sunpy/pull/2655>`__)
- When running image tests, a comparison HTML page is now generated to show
  the generated images and expected images. (`#2660 <https://github.com/sunpy/sunpy/pull/2660>`__)
- Change to using pytest-cov for coverage report generation to enable support for parallel builds (`#2667 <https://github.com/sunpy/sunpy/pull/2667>`__)
- Use of `textwrap` to keep source code indented when multiline texts is used (`#2671 <https://github.com/sunpy/sunpy/pull/2671>`__)
- Fix mispelling of private attribute ``_default_heliographic_latitude`` in map. (`#2730 <https://github.com/sunpy/sunpy/pull/2730>`__)
- Miscellaneous fixes to developer docs about building sunpy's documentation. (`#2825 <https://github.com/sunpy/sunpy/pull/2825>`__)
- Changed ``sunpy.instr.aia.aiaprep`` to update BITPIX keyword to reflect the float64 dtype. (`#2831 <https://github.com/sunpy/sunpy/pull/2831>`__)
- Remove warning from ``GenericMap.submap`` when using pixel ``Quantities`` as input. (`#2833 <https://github.com/sunpy/sunpy/pull/2833>`__)
- Remove the usage of six and all ``__future__`` imports (`#2837 <https://github.com/sunpy/sunpy/pull/2837>`__)
- Fix SunPy Coordinate tests with Astropy 3.1 (`#2838 <https://github.com/sunpy/sunpy/pull/2838>`__)
- Stores entries from directories into database sorted by name. It adds mocks to the database user guide examples. (`#2873 <https://github.com/sunpy/sunpy/pull/2873>`__)
- Fix all DeprecationWarning: invalid escape sequence. (`#2885 <https://github.com/sunpy/sunpy/pull/2885>`__)
- Used `unittest.mock` for creating offline tests for simulating online tests for :file:`test_noaa.py` (`#2900 <https://github.com/sunpy/sunpy/pull/2900>`__)
- Fix support for pip 19 and isolated builds (`#2915 <https://github.com/sunpy/sunpy/pull/2915>`__)
- Moved to using `AppDirs <https://github.com/ActiveState/appdirs>`__ as the place to host our configuration file. (`#2922 <https://github.com/sunpy/sunpy/pull/2922>`__)
- Users can now use fewer keywords in our `~sunpy.net.helioviewer.HelioviewerClient` to access the available sources. Either by ``observatory`` and ``measurement`` or ``instrument`` and ``measurement`` as this much information is enough to get the source ID for most of the cases. (`#2926 <https://github.com/sunpy/sunpy/pull/2926>`__)
- Remove the pytest dependency on the ``GenericMap`` asdf tag. (`#2943 <https://github.com/sunpy/sunpy/pull/2943>`__)
- Fix initialization of `~sunpy.net.vso.VSOClient` when no WSDL link is found. (`#2981 <https://github.com/sunpy/sunpy/pull/2981>`__)


0.9.0
=====

New Features
------------

- Added TimeUTime class to support utime. [#2409]
- Example for fine-grained use of ticks and grids [#2435]
- Maintiners Workflow Guide [#2411]
- Decorator to append and/or prepend doc strings [#2386]
- Adding ``python setup.py test --figure-only`` [#2557]
- Fido.fetch now accepts pathlib.Path objects for path attribute.[#2559]
- The `~sunpy.coordinates.HeliographicStonyhurst` coordinate system can now be specified
  using a cartesian system, which is sometimes known as the
  "Heliocentric Earth equatorial" (HEEQ) coordinate system. [#2437]

API Changes
-----------

- ``sunpy.coordinates.representation`` has been removed. Longitude wrapping is now done in the constructor of the frames. [#2431]
- Propagation of ``obstime`` in the coordinate frame transformation has changed, this means in general when transforming directly between frames (not `~astropy.coordinates.SkyCoord`) you will have to specify ``obstime`` in more places. [#2461]
- Transforming between Heliographic Stonyhurst and Carrington now requires that ``obstime`` be defined and the same on both the input and output frames. [#2461]
- Removed the figure return from .peek() [#2487]

Bug Fixes
---------

- Improve TimeSeriesBase docstring [#2399]
- Validate that pytest-doctestplus is installed [#2388]
- Fix use of self.wcs in plot in mapbase [#2398]
- Updated docstring with pointer to access EVE data for other levels [#2402]
- Fix broken links and redirections in documentation [#2403]
- Fixes Documentation changes due to NumPy 1.14 [#2404]
- Added docstrings to functions in dowload.py [#2415]
- Clean up database doc [#2414]
- rhessi.py now uses sunpy.io instead of astropy.io [#2416]
- Remove Gamma usage in Map [#2424]
- Changed requirements to python-dateutil [#2426]
- Clarify coordinate system definitions [#2429]
- Improve Map Peek when using draw_grid [#2442]
- Add HCC --> HGS test [#2443]
- Testing the transformation linking SunPy and Astropy against published values [#2454]
- Fixed title bug in sunpy.timeseries.rhessi [#2477]
- Allow LineAnimator to accept a varying x-axis [#2491]
- Indexing Bug Fix to LineAnimator [#2560]
- Output sphinx warnings to stdout [#2553]
- Docstring improvement for LineAnimator [#2514]
- move the egg_info builds to circleci [#2512]
- Added tests for TraceMap [#2504]
- Fix HGS frame constructor and HPC ``calculate_distance`` with SkyCoord constructor. [#2463]
- removed ``wavelnth`` keyword in meta desc of Maps to avoid using non standard FITS keyword like ``nan`` [#2456]
- The documentation build now uses the Sphinx configuration from sphinx-astropy rather than from astropy-helpers.[#2494]
- Migrate to hypothesis.strategies.datetimes [#2368]
- Prevent a deprecation warning due to truth values of Quantity [#2358]
- Print a warning when heliographic longitude is set to it's default value of 0 [#2480]
- parse_time now parses numpy.datetime64 correctly. [#2572]

0.8.5
=====

Bug Fixes
---------

- Removed AstropyDeprecationWarning from sunpy.coordinates.representation [#2476]
- Fix for NorthOffsetFrame under Astropy 3.0 [#2486]
- Fix lightcurve tests under numpy dev [#2505]
- Updated depecration link of radiospectra [#2481]
- Fixed Padding values in some of the documentation pages [#2497]
- Move documentation build to circleci [#2509]
- Fix Issue #2470 hgs_to_hcc(heliogcoord, heliocframe) [#2502]
- Fixing CompositeMap object so that it respects masked maps [#2492]

0.8.4
=====

Bug Fixes
---------

- Improve detection of ``SkyCoord`` frame instantiation when distance is
  ``1*u.one``. This fixes a plotting bug with ``WCSAxes`` in Astropy 3.0 [#2465]
- removed ``wavelnth`` keyword in meta desc of Maps to avoid using non standard FITS keyword like ``nan`` [#2427]
- Change the default units for HPC distance from ``u.km`` to `None`. [#2465]

0.8.3
=====

Bug Fixes
---------

- `~sunpy.net.dataretriever.XRSClient` now reports time ranges of files correctly. [#2364]
- Make parse_time work with datetime64s and pandas series [#2370]
- CompositeMap axes scaling now uses map spatial units [#2310]
- Moved license file to root of repository and updated README file [#2326]
- Fix docstring formatting for net.vso.attrs [#2309]]
- Fix coloring of ticks under matplotlib 2.0 default style [#2320]
- Always index arrays with tuples in ``ImageAnimator`` [#2320]
- Added links to possible attrs for FIDO in guide [#2317] [#2289]
- Updated GitHub Readme [#2281] [#2283]
- Fix matplotlib / pandas 0.21 bug in examples [#2336]
- Fixes the off limb enhancement example [#2329]
- Changes to masking hot pixels and picking bright pixels examples [#2325] [#2319]
- Travis CI fix for numpy-dev build [#2340]
- Updated masking brightest pixel example [#2338]
- Changed TRAVIS cronjobs [#2338]
- Support array values for ``obstime`` for coordinates and transformations [#2342] [#2346]
- Updated Gallery off limb enhance example [#2337]
- Documentation fixes for VSO [#2354] [#2353]
- All tests within the documentation have been fixed [#2343]
- Change to using pytest-remotedata for our online tests [#2345]
- Fixed upstream astropy/numpy documentation issues [#2359]
- Documentation for Map improved [#2361]
- Fix the output units of pixel_to_world [#2362]
- Documentation for Database improved [#2355]
- Added test for mapsave [#2365]
- Documentation for Sun improved [#2369]

0.8.2
=====

Bug Fixes
---------

- Shows a warning if observation time is missing [#2293]
- Updates MapCube to access the correct properties of the namedtuple SpatialPair [#2297]

0.8.1
======

Bug fixes
---------

- Fixed TimeSeries test failures due to missing test files [#2273]
- Refactored a GOES test to avoid a Py3.6 issue [#2276]

0.8.0
======

New Features
------------

-  Solar differential rotation for maps and submaps included.
-  Solar rotation calculation and mapcube derotation now use sunpy coordinates.
-  Sample data now downloads automatically on import if not available
   and is now pluggable so can be used by affiliated packages. Shortcut
   names have been normalized and all LIGHTCURVE shortcuts have changed
   to TIMESERIES.
-  Calculation of points on an arc of a great circle connecting two
   points on the Sun.
-  Removed ``extract_time`` function from ``sunpy.time`` and also tests
   related to the function from ``sunpy.time.tests``
-  User can now pass a custom time format as an argument inside
   ``sunpy.database.add_from_dir()`` in case the ``date-obs`` metadata
   cannot be read automatically from the files.
-  Add time format used by some SDO HMI FITS keywords
-  Now the ``sunpy.database.tables.display_entries()`` prints an astropy
   table.
-  Additional methods added inside the ``sunpy.database`` class to make
   it easier to display the database contents.
-  Remove unused ``sunpy.visualization.plotting`` module
-  Port the pyana wrapper to Python 3
-  ``Map.peek(basic_plot-True)`` no longer issues warnings
-  Remove the ``sunpy.map.nddata_compat`` module, this makes
   ``Map.data`` and ``Map.meta`` read only.
-  Add a ``NorthOffsetFrame`` class for generating HGS-like coordinate
   systems with a shifted north pole.
-  Remove deprecated ``VSOClient.show`` method.
-  Deprecate ``sunpy.wcs``: ``sunpy.coordinates`` and ``sunpy.map`` now
   provide all that functionality in a more robust manner.
-  Added hdu index in ``sunpy.database.tables.DatabaseEntry`` as a
   column in the table.
-  Removed ``HelioviewerClient`` from the ``sunpy.net`` namespace. It
   should now be imported with
   ``from sunpy.net.helioviewer import HelioviewerClient``.
-  Removed compatibility with standalone ``wcsaxes`` and instead depend
   on the version in astropy 1.3. SunPy now therefore depends on
   astropy>-1.3.
-  Update to ``TimeRange.__repr__``; now includes the qualified name and
   ``id`` of the object.
-  A new ``sunpy.visualization.imageanimator.LineAnimator`` class has
   been added to animate 1D data. This has resulted in API change for
   the ``sunpy.visualization.imageanimator.ImageAnimator`` class. The
   updateimage method has been renamed to update\_plot.
-  Drop support for Python 3.4.
-  SunPy now requires WCSAxes and Map.draw\_grid only works with
   WCSAxes.
-  ``Helioprojective`` and ``HelioCentric`` frames now have an
   ``observer`` attribute which itself is a coordinate object
   (``SkyCoord``) instead of ``B0``, ``L0`` and ``D0`` to describe the
   position of the observer.
-  ``GenericMap.draw_grid`` now uses ``WCSAxes``, it will only work on a
   ``WCSAxes`` plot, this may be less performant than the previous
   implementation.
-  ``GenericMap.world_to_pixel`` and ``GenericMap.pixel_to_world`` now
   accept and return ``SkyCoord`` objects only.
-  ``GenericMap`` has a new property ``observer_coordinate`` which
   returns a ``SkyCoord`` describing the position of the observer.
-  ``GenericMap.submap`` now takes arguments of the form ``bottom_left``
   and ``top_right`` rather than ``range_a`` and ``range_b``. This
   change enables submap to properly handle rotated maps and take input
   in the form of ``SkyCoord`` objects.
-  When referring to physical coordinates ``Pair.x`` has been replaced
   with ``SpatialPair.axis1``. This means values returned by
   ``GenericMap`` now differentiate between physical and pixel
   coordinates.
-  The physical radius of the Sun (length units) is now passed from Map
   into the coordinate frame so a consistent value is used when
   calculating distance to the solar surface in the
   ``HelioprojectiveFrame`` coordinate frame.
-  A new ``sunpy.visualization.imageanimator.ImageAnimatorWCS`` class
   has been added to animate N-Dimensional data with the associated WCS
   object.
-  Moved Docs to docs/ to follow the astropy style
-  Added SunPy specific warnings under util.
-  SunPy coordinate frames can now be transformed to and from Astropy
   coordinate frames
-  The time attribute for SunPy coordinate frames has been renamed from
   ``dateobs`` to ``obstime``
-  Ephemeris calculations with higher accuracy are now available under
   ``sunpy.coordinates.ephemeris``
-  Add support for SunPy coordinates to specify observer as a string of
   a major solar-system body, with the default being Earth. To make
   transformations using an observer specified as a string, ``obstime``
   must be set.
-  Added VSO query result block level caching in the database module.
   This prevents re-downloading of files which have already been
   downloaded. Especially helpful in case of overlapping queries.
-  Change the default representation for the Heliographic Carrington
   frame so Longitude follows the convention of going from 0-360
   degrees.
-  All Clients that are able to search and download data now have a
   uniform API that is ``search`` and ``fetch``. The older functions are
   still there but are deprecated for 0.8.

Bug fixes
---------

-  Add tests for RHESSI instrument
-  Maps from Helioviewer JPEG2000 files now have correct image scaling.
-  Get and set methods for composite maps now use Map plot\_settings.
-  Simplified map names when plotting.
-  Fix bug in ``wcs.convert_data_to_pixel`` where crpix[1] was used for
   both axes.
-  Fix some leftover instances of ``GenericMap.units``
-  Fixed bugs in ``sun`` equations
-  ``sunpy.io.fits.read`` will now return any parse-able HDUs even if
   some raise an error.
-  ``VSOClient`` no longer prints a lot of XML junk if the query fails.
-  Fix Map parsing of some header values to allow valid float strings
   like 'nan' and 'inf'.
-  Fix Map parsing of some header values to allow valid float strings
   like 'nan' and 'inf'.

0.7.8
=====

-  The SunPy data directory "~/sunpy" is no longer created until it is
   used (issue #2018)
-  Change the default representation for the Heliographic Carrington
   frame so Longitude follows the convention of going from 0-360
   degrees.
-  Fix for surface gravity unit.
-  Support for Pandas 0.20.1

0.7.7
=====

-  Fix errors with Numpy 1.12

0.7.6
=====

-  Add Astropy 1.3 Support

0.7.5
=====

-  Fix test faliure (mapbase) with 1.7.4
-  Restrict supported Astropy version to 1.0<astropy<1.3
-  Add Figure test env to SunPy repo.

0.7.4
=====

-  Remove Map always forcing warnings on.
-  ``Map.center`` now uses ``Map.wcs`` to correctly handle rotation.
-  Fix link in coordinates documentation.
-  Update helioviewer URL to HTTPS (fixes access to Helioviewer).
-  Fix processing of TRACE and YOHKOH measurement properties.
-  Remove warnings when using ``Map.peek(basic_plot-True)``
-  Update docstrings for HPC and HCC frames.

0.7.3
=====

-  Fix ConfigParser for Python 3.5.2 - This allows SunPy to run under
   Python 3.5.2
-  Fix incorrect ordering of keys in ``MapMeta``
-  Add ``sunpy.util.scraper`` to the API documentation.

0.7.2
=====

-  Fixed bugs in ``sun`` equations

0.7.1
=====

-  Fix bug in ``wcs.convert_data_to_pixel`` where crpix[1] was used for
   both axes.
-  Fix some leftover instances of ``GenericMap.units``
-  Fixed bugs in ``sun`` equations
-  Now the ``sunpy.database.tables.display_entries()`` prints an astropy
   table.
-  Additional methods added inside the ``sunpy.database`` class to make
   it easier to display the database contents.
-  ``sunpy.io.fits.read`` will now return any parse-able HDUs even if
   some raise an error.
-  ``VSOClient`` no longer prints a lot of XML junk if the query fails.
-  Remove unused ``sunpy.visualization.plotting`` module
-  ``Map.peek(basic_plot-True)`` no longer issues warnings
-  Remove the ``sunpy.map.nddata_compat`` module, this makes
   ``Map.data`` and ``Map.meta`` read only.
-  Add a ``NorthOffsetFrame`` class for generating HGS-like coordinate
   systems with a shifted north pole.
-  Remove deprecated ``VSOClient.show`` method.
-  Deprecate ``sunpy.wcs``: ``sunpy.coordinates`` and ``sunpy.map`` now
   provide all that functionality in a more robust manner.
-  Added hdu index in ``sunpy.database.tables.DatabaseEntry`` as a
   column in the table.
-  Removed ``HelioviewerClient`` from the ``sunpy.net`` namespace. It
   should now be imported with
   ``from sunpy.net.helioviewer import HelioviewerClient``.
-  Removed compatibility with standalone ``wcsaxes`` and instead depend
   on the version in astropy 1.3. SunPy now therefore depends on
   astropy>-1.3.
-  Update to ``TimeRange.__repr__``; now includes the qualified name and
   ``id`` of the object.
-  Change the default representation for the Heliographic Carrington
   frame so Longitude follows the convention of going from 0-360
   degrees.
-  Fix Map parsing of some header values to allow valid float strings
   like 'nan' and 'inf'.

0.7.0
=====

-  Fixed test failures with numpy developer version.[#1808]
-  Added ``timeout`` parameter in ``sunpy.data.download_sample_data()``
-  Fixed ``aiaprep`` to return properly sized map.
-  Deprecation warnings fixed when using image coalignment.
-  Sunpy is now Python 3.x compatible (3.4 and 3.5).
-  Added a unit check and warnings for map metadata.
-  Added IRIS SJI color maps.
-  Updated ``show_colormaps()`` with new string filter to show a subset
   of color maps.
-  Fixed MapCube animations by working around a bug in Astropy's
   ImageNormalize
-  Remove ``vso.QueryResponse.num_records()`` in favour of ``len(qr)``
-  Add a ``draw_rectangle`` helper to ``GenericMap`` which can plot
   rectangles in the native coordinate system of the map.
-  Added the ability to shift maps to correct for incorrect map
   location, for example.
-  Bug fix for RHESSI summary light curve values.
-  Mapcube solar derotation and coalignment now pass keywords to the
   routine used to shift the images, scipy.ndimage.interpolation.shift.
-  Add automatic registration of ``GenericMap`` subclasses with the
   factory as long as they define an ``is_datasource_for`` method.
-  Added functions ``flareclass_to_flux`` and ``flux_to_flareclass``
   which convert between GOES flux to GOES class numbers (e.g. X12,
   M3.4).
-  Removed old ``sunpy.util.goes_flare_class()``
-  Bug fix for RHESSI summary light curve values.
-  The ``MapCube.as_array`` function now returns a masked numpy array if
   at least one of the input maps in the MapCube has a mask.
-  Map superpixel method now respects maps that have masks.
-  Map superpixel method now accepts numpy functions as an argument, or
   any user-defined function.
-  Map superpixel method no longer has the restriction that the number
   of original pixels in the x (or y) side of the superpixel exactly
   divides the number of original pixels in the x (or y) side of the
   original map data.
-  ``sunpy.physics.transforms`` has been deprecated and the code moved
   into ``sunpy.physics``.
-  Add the ``sunpy.coordinates`` module, this adds the core physical
   solar coordinates frame within the astropy coordinates framework.
-  Added ability of maps to draw contours on top of themselves
   (``draw_contours``)
-  Added concatenate functionality to lightcurve base class.
-  Fix Map to allow astropy.io.fits Header objects as valid input for
   meta arguments.
-  Added an examples gallery using ``sphinx-gallery``.
-  API clean up to constants. Removed constant() function which is now
   replaced by get().
-  Prevent helioviewer tests from checking access to the API endpoint
   when running tests offline.
-  ``GenericMap.units`` is renamed to ``GenericMap.spatial_units`` to
   avoid confusion with ``NDData.unit``.
-  ``GenericMap`` now has a ``coordinate_frame`` property which returns
   an ``astropy.coordinates`` frame with all the meta data from the map
   populated.
-  ``GenericMap`` now has a ``_mpl_axes`` method which allows it to be
   specified as a projection to ``matplotlib`` methods and will return a
   ``WCSAxes`` object with ``WCS`` projection.

0.6.5
=====

-  The draw\_grid keyword of the peek method of Map now accepts booleans
   or astropy quantities.
-  Fix bug in ``wcs.convert_data_to_pixel`` where crpix[1] was used for
   both axes.
-  Fixed bugs in ``sun`` equations

0.6.4
=====

-  Bug fix for rhessi summary lightcurve values.
-  Fix docstring for ``pixel_to_data`` and ``data_to_pixel``.
-  Fix the URL for the Helioviewer API. (This fixes Helioviewer.)
-  Fix the way ``reshape_image_to_4d_superpixel`` checks the dimension
   of the new image.
-  Fix Map to allow astropy.io.fits Header objects as valid input for
   meta arguments.
-  Prevent helioviewer tests from checking access to API when running
   tests in offline mode.

0.6.3
=====

-  Change setup.py extras to install suds-jurko not suds.

0.6.2
=====

-  Changed start of GOES 2 operational time range back to 1980-01-04 so
   data from 1980 can be read into GOESLightCurve object
-  Fix bug with numpy 1.10
-  update astropy\_helpers
-  Added new sample data

0.6.1
=====

-  Fixed MapCube animations by working around a bug in Astropy's
   ImageNormalize
-  Small fix to RTD builds for Affiliated packages
-  SunPy can now be installed without having to install Astropy first.
-  MapCubes processed with ``coalignment.apply_shifts`` now have correct
   metadata.
-  Multiple fixes for WCS transformations, especially with solar-x,
   solar-y CTYPE headers.

0.6.0
=====

-  Enforced the use of Astropy Quantities through out most of SunPy.
-  Dropped Support for Python 2.6.
-  Remove old style string formatting and other 2.6 compatibility lines.
-  Added vso like querying feature to JSOC Client.
-  Refactor the JSOC client so that it follows the .query() .get()
   interface of VSOClient and UnifedDownloader.
-  Provide ``__str__`` and ``__repr__`` methods on vso ``QueryResponse``
   deprecate ``.show()``.
-  Downloaded files now keep file extensions rather than replacing all
   periods with underscores.
-  Update to TimeRange API, removed t1 and t0, start and end are now
   read-only attributes.
-  Added ability to download level3 data for lyra Light Curve along with
   corresponding tests.
-  Added support for gzipped FITS files.
-  Add STEREO HI Map subclass and color maps.
-  Map.rotate() no longer crops any image data.
-  For accuracy, default Map.rotate() transformation is set to
   bi-quartic.
-  ``sunpy.image.transform.affine_transform`` now casts integer data to
   float64 and sets NaN values to 0 for all transformations except
   scikit-image rotation with order <- 3.
-  CD matrix now updated, if present, when Map pixel size is changed.
-  Removed now-redundant method for rotating IRIS maps since the
   functionality exists in Map.rotate()
-  Provide ``__str__`` and ``__repr__`` methods on vso ``QueryResponse``
   deprecate ``.show()``
-  SunPy colormaps are now registered with matplotlib on import of
   ``sunpy.cm``
-  ``sunpy.cm.get_cmap`` no longer defaults to 'sdoaia94'
-  Added database url config setting to be setup by default as a sqlite
   database in the sunpy working directory
-  Added a few tests for the sunpy.roi module
-  Added capability for figure-based tests
-  Removed now-redundant method for rotating IRIS maps since the
   functionality exists in Map.rotate().
-  SunPy colormaps are now registered with matplotlib on import of
   ``sunpy.cm``.
-  ``sunpy.cm.get_cmap`` no longer defaults to 'sdoaia94'.
-  Added database url config setting to be setup by default as a sqlite
   database in the sunpy working directory.
-  Added a few tests for the sunpy.roi module.
-  Refactored mapcube co-alignment functionality.
-  Removed sample data from distribution and added ability to download
   sample files
-  Changed start of GOES 2 operational time range back to 1980-01-04 so
   data from 1980 can be read into GOESLightCurve object
-  Require JSOC request data calls have an email address attached.
-  Calculation of the solar rotation of a point on the Sun as seen from
   Earth, and its application to the de-rotation of mapcubes.
-  Downloaded files now keep file extensions rather than replacing all
   periods with underscores
-  Fixed the downloading of files with duplicate names in sunpy.database
-  Removed sample data from distribution and added ability to download
   sample files.
-  Added the calculation of the solar rotation of a point on the Sun as
   seen from Earth, and its application to the de-rotation of mapcubes.
-  Changed default for GOESLightCurve.create() so that it gets the data
   from the most recent existing GOES fits file.
-  Map plot functionality now uses the mask property if it is present,
   allowing the plotting of masked map data
-  Map Expects Quantities and returns quantities for most parameters.
-  Map now used Astropy.wcs for world <-> pixel conversions.
-  map.world\_to\_pixel now has a similar API to map.pixel\_to\_world.
-  map.shape has been replaced with map.dimensions, which is ordered x
   first.
-  map.rsun\_arcseconds is now map.rsun\_obs as it returns a quantity.
-  Map properties are now named tuples rather than dictionaries.
-  Improvement for Map plots, standardization and improved color tables,
   better access to plot variables through new plot\_settings variable.
-  Huge improvements in Instrument Map doc strings. Now contain
   instrument descriptions as well as reference links for more info.
-  net.jsoc can query data series with time sampling by a Sample
   attribute implemented in vso.
-  MapCube.plot and MapCube.peek now support a user defined
   plot\_function argument for customising the animation.
-  Added new sample data file, an AIA cutout file.
-  Moved documentation build directory to doc/build

0.5.5
=====

-  Changed default for GOESLightCurve.create() so that it gets the data
   from the most recent existing GOES fits file.
-  Improvements to the Map documentation.
-  Typo fixes in sunpy.wcs documentation.

0.5.4
=====

-  ``sunpy.image.transform.affine_transform`` now casts integer data to
   float64 and sets NaN values to 0 for all transformations except
   scikit-image rotation with order <- 3.
-  Updated SWPC/NOAA links due to their new website.
-  Exposed the raw AIA color tables in ``sunpy.cm.color_tables``.
-  Fixes ``map`` compatibility with Astropy 1.0.x.

0.5.3
=====

-  Goes peek() plot now works with matplotlib 1.4.x
-  The ANA file reading C extensions will no longer compile under
   windows. Windows was not a supported platform for these C extensions
   previously.

0.5.2
=====

-  If no CROTA keyword is specified in Map meta data, it will now
   default to 0 as specified by the FITS WCS standard.
-  Map now correctly parses and converts the CD matrix, as long as CDELT
   is specified as well. (Fixes SWAP files)
-  Fix of HELIO webservice URLs
-  MapCube.plot() is now fixed and returns a
   matplotlib.animation.FuncAnimation object.

0.5.1
=====

-  MAJOR FIX: map.rotate() now works correctly for all submaps and off
   center rotations.
-  HELIO URL updated, querys should now work as expected.
-  All tabs removed from the code base.
-  All tests now use tempfile rather than creating files in the current
   directory.
-  Documentation builds under newer sphinx versions.
-  ANA and JP2 tests are skipped if dependencies are missing.
-  ANA tests are skipped on windows.

0.5.0
=====

-  Added additional functionality to the GOES module i.e. the ability to
   calculate GOES temperature and emission measure from GOES fluxes.
-  changed \_maps attribute in MapCube to a non-hidden type
-  Added Nobeyama Radioheliograph data support to Lightcurve object.
-  Fixed some tests on map method to support Windows
-  Added a window/split method to time range
-  Updates to spectrogram documentation
-  Added method Database.add\_from\_hek\_query\_result to HEK database
-  Added method Database.download\_from\_vso\_query\_result
-  GOES Lightcurve now makes use of a new source of GOES data, provides
   metadata, and data back to 1981.
-  Removed sqlalchemy as a requirement for SunPy
-  Added support for NOAA solar cycle prediction in lightcurves
-  Some basic tests for GenericLightCurve on types of expected input.
-  Fix algorithm in sunpy.sun.equation\_of\_center
-  Added Docstrings to LightCurve methods.
-  Added tests for classes in sunpy.map.sources. Note that some classes
   (TRACE, RHESSI) were left out because SunPy is not able to read their
   FITS files.
-  Added functions that implement image coalignment with support for
   MapCubes.
-  Cleaned up the sunpy namespace, removed .units, /ssw and .sphinx.
   Also moved .coords .physics.transforms.
-  Added contains functionality to TimeRange module
-  Added t-'now' to parse\_time to privide utcnow datetime.
-  Fixed time dependant functions (.sun) to default to t-'now'
-  Fixed solar\_semidiameter\_angular\_size
-  Improved line quality and performances issues with map.draw\_grid()
-  Remove deprecated make\_map command.

0.4.2
=====

-  Fixes to the operational range of GOES satellites
-  Fix the URL for HELIO queries.

0.4.1
=====

-  Fix map.rotate() functionality
-  Change of source for GOES data.
-  Fix EIT test data and sunpy FITS saving
-  Some documentation fixes
-  fix file paths to use os.path.join for platform independance.

0.4.0
=====

-  **Major** documentation refactor. A far reaching re-write and
   restructure.
-  Add a SunPy Database to store and search local data.
-  Add beta support for querying the HELIO HEC
-  Add beta HEK to VSO query translation.
-  Add the ability to download the GOES event list.
-  Add support for downloading and querying the LYTAF database.
-  Add support for ANA data.
-  Updated sun.constants to use astropy.constants objects which include
   units, source, and error instide. For more info check out
   https://docs.astropy.org/en/latest/constants/index.html
-  Add some beta support for IRIS data products
-  Add a new MapCubeAnimator class with interactive widgets which is
   returned by mapcube.peek().
-  The Glymur library is now used to read JPEG2000 files.
-  GOESLightCurve now supports all satellites.
-  Add support for VSO queries through proxies.
-  Fix apparent Right Ascension calulations.
-  LightCurve meta data member now an OrderedDict Instance

0.3.2
=====

-  Pass draw\_limb arguments to patches.Circle
-  Pass graw\_grid arguments to pyplot.plot()
-  Fix README code example
-  Fix Documentation links in potting guide
-  Update to new EVE data URL
-  Update LogicalLightcurve example in docs
-  Improved InteractiveVSOClient documentation
-  GOESLightCurve now fails politely if no data is avalible.

Known Bugs:

-  sunpy.util.unit\_conversion.to\_angstrom does not work if 'nm' is
   passed in.

0.3.1
=====

-  Bug Fix: Fix a regression in CompositeMap that made contor plots
   fail.
-  Bug Fix: Allow Map() to accept dict as metadata.
-  Bug Fix: Pass arguments from Map() to io.read\_file.

0.3.0
=====

-  Removal of Optional PIL dependency
-  Parse\_time now looks through nested lists/tuples
-  Draw\_limb and draw\_grid are now implemented on MapCube and
   CompositeMap
-  Caculations for differential roation added
-  mapcube.plot() now runs a mpl animation with optional controls
-  A basic Region of Interest framework now exists under sunpy.roi
-  STEREO COR colour maps have been ported from solarsoft.
-  sunpy.time.timerange has a split() method that divides up a time
   range into n equal parts.
-  Added download progress bar
-  pyfits is depricated in favor of Astropy

spectra:

-  Plotting has been refactorted to use a consistent interface
-  spectra now no-longer inherits from numpy.ndarray instead has a .data
   attribute.

Map: \* map now no-longer inherits from numpy.ndarray instead has a
.data attribute. \* make\_map is deprecated in favor of Map which is a
new factory class \* sunpy.map.Map is now sunpy.map.GenericMap \*
mymap.header is now mymap.meta \* attributes of the map class are now
read only, changes have to be made through map.meta \* new MapMeta class
to replace MapHeader, MapMeta is not returned by sunpy.io \* The
groundwork for GenericMap inherting from astropy.NDData has been done,
there is now a NDDataStandin class to provide basic functionality.

io: \* top level file\_tools improved to be more flexible and support
multiple HDUs \* all functions in sunpy.io now assume mutliple HDUs,
even JP2 ones. \* there is now a way to override the automatic filetype
detection \* Automatic fits file detection improved \* extract\_waveunit
added to io.fits for detection of common ways of storing wavelength unit
in fits files.

-  A major re-work of all interal imports has resulted in a much cleaner
   namespace, i.e. sunpy.util.util is no longer used to import util.
-  Some SOHO and STEREO files were not reading properly due to a
   date\_obs parameter.
-  Sunpy will now read JP2 files without a comment parameter.
-  Memory leak in Crotate patched
-  Callisto: Max gap between files removed

0.2.0
=====

-  Completely re-written plotting routines for most of the core
   datatypes.
-  JPEG 2000 support as an input file type.
-  Improved documentation for much of the code base, including
   re-written installation instructions.
-  New lightcurve object

   -  LYRA support
   -  GOES/XRS support
   -  SDO/EVE support

-  New Spectrum and Spectrogram object (in development)

   -  Spectrogram plotting routines
   -  Callisto spectrum type and support
   -  STEREO/SWAVES support

-  Map Object

   -  Added support for LASCO, Yohkoh/XRT maps
   -  A new CompositeMap object for overlaying maps
   -  Resample method
   -  Superpixel method
   -  The addition of the rotate() method for 2D maps.
********
`SunPy`_
********

|Latest Version| |codecov| |matrix| |Research software impact| |DOI| |Powered by NumFOCUS|

.. |Latest Version| image:: https://img.shields.io/pypi/v/sunpy.svg
   :target: https://pypi.python.org/pypi/sunpy/
.. |matrix| image:: https://img.shields.io/matrix/sunpy:openastronomy.org.svg?colorB=%23FE7900&label=Chat&logo=matrix&server_fqdn=openastronomy.modular.im
   :target: https://openastronomy.element.io/#/room/#sunpy:openastronomy.org
.. |codecov| image:: https://codecov.io/gh/sunpy/sunpy/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/sunpy/sunpy
.. |Research software impact| image:: http://depsy.org/api/package/pypi/sunpy/badge.svg
   :target: http://depsy.org/package/python/sunpy
.. |DOI| image:: https://zenodo.org/badge/2165383.svg
   :target: https://zenodo.org/badge/latestdoi/2165383
.. |Powered by NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: https://numfocus.org
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/sunpy/sunpy/main?filepath=examples

SunPy is an open-source Python library for Solar Physics data analysis and visualization.
Our homepage `SunPy`_ has more information about the project.

For some examples of using SunPy see our `gallery`_, to see the latest changes in SunPy see our `Changelog`_.

.. _SunPy: https://sunpy.org
.. _gallery: https://docs.sunpy.org/en/stable/generated/gallery/index.html
.. _Changelog: https://docs.sunpy.org/en/stable/whatsnew/changelog.html

Installation
============

The recommended way to install SunPy is with `miniconda`_.
To install SunPy once conda is installed run the following two commands:

.. code:: bash

    $ conda config --append channels conda-forge
    $ conda install sunpy

For detailed installation instructions, see the `installation guide`_ in the SunPy docs.

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _installation guide: https://docs.sunpy.org/en/stable/guide/installation.html

Developing
==========

If you want to develop SunPy you will need to install from GitHub.
For detailed installation instructions, see `Development installation`_ in the SunPy docs.

Usage
=====

Here is a quick example of plotting an AIA image:

.. code:: python

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE
    >>> aia = sunpy.map.Map(AIA_171_IMAGE)
    >>> aia.peek()

Getting Help
============

For more information or to ask questions about SunPy, check out:

-  `SunPy Documentation`_
-  `SunPy Element Channel`_
-  `SunPy Mailing List`_

.. _SunPy Documentation: https://docs.sunpy.org/en/stable/
.. _SunPy Element Channel: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _SunPy Mailing List: https://groups.google.com/forum/#!forum/sunpy

Acknowledging or Citing sunpy
=============================

If you use `sunpy` in your scientific work, we would appreciate citing it in your publications.
The continued growth and development of `sunpy` is dependent on the community being aware of `sunpy`.

Please see https://sunpy.org/about#acknowledging-or-citing-sunpy on how to do this.

Contributing
============

|Open Source Helpers|

If you would like to get involved, start by joining the `SunPy mailing list`_ and check out the `Developers Guide`_ section of the SunPy docs.
Stop by our chat room `#sunpy:openastronomy.org`_ if you have any questions.
Help is always welcome so let us know what you like to work on, or check out the `issues page`_ for the list of known outstanding items.

For more information on contributing to SunPy, please read our `Newcomers' guide`_.

.. |Open Source Helpers| image:: https://www.codetriage.com/sunpy/sunpy/badges/users.svg
   :target: https://www.codetriage.com/sunpy/sunpy

.. _SunPy mailing list: https://groups.google.com/forum/#!forum/sunpy
.. _Developers Guide: https://docs.sunpy.org/en/latest/dev_guide/index.html
.. _`#sunpy:openastronomy.org`: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _issues page: https://github.com/sunpy/sunpy/issues
.. _Newcomers' guide: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html
.. _Development installation:  https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html#setting-up-a-development-environment

Code of Conduct
===============

When you are interacting with the SunPy community you are asked to follow our `Code of Conduct`_.

.. _Code of Conduct: https://sunpy.org/coc
Fixed a bug where custom values in the ``plot_settings`` dictionary were not being propagated
to new map instances created when calling map methods (e.g. ``.submap``).
Update asdf schemas so that references use URIs not tags as this is not
supported by the new asdf extensions API.
Added support for ``"%Y.%m.%d_%H:%M:%S_UTC"`` and ``"%Y.%m.%d_%H:%M:%S"`` time formats in `sunpy.time.parse_time`.
Fixed reading CDF files when a column has no entries. If this is the case the
column will be ignored, and a message logged at DEBUG level.
Printing a `.MetaDict`  will now show each entry on a new line.
Fixed plotting and peeking NORH timeseries data with ``pandas`` 1.4.0.
`~sunpy.util.exceptions.SunpyWarning` is no longer a subclass of `~astropy.utils.exceptions.AstropyWarning`.
In the case where a map header has no PC_ij values, CROTA2 != 0, and
CDELT1 != CDELT2, the calculation of the map rotation matrix has been fixed.
This bug only affected maps with non-zero rotation, no PC matrix in the header,
and un-equal scales along the two image axes.
The default ``id_type`` in :func:`sunpy.coordinates.get_horizons_coord` is now
`None` to match the deafult ``id_type`` in astroquery 0.4.4, which will search
major bodies first, and if no major bodies are found, then search small bodies.
For older versions of astroquery the default ``id_type`` used by
:func:`~sunpy.coordinates.get_horizons_coord` is still ``'majorbody'``.
Running the tests now requires the ``pytest-xdist`` package. By
default tests are *not* run in parallel, but can be configured to do so
using ``pytest-xdist`` command line options.
Reading a series of CDF files where at least one of them is empty no longer
raises an error. A message for each empty file is logged at the DEBUG level.
The ``rsun`` argument to :func:`~sunpy.map.get_observer_meta` is now
optional.
`sunpy.map.GenericMap.wcs` now checks that the scale property has the correct
units whilst constructing the WCS.
Fixed various plotting issues with the gallery example :ref:`sphx_glr_generated_gallery_units_and_coordinates_AIA_limb_STEREO.py`.
Increased the default maximum amount of records returned from HEC to 500 from 10.
If the maximum number of records are returned, a message is shown.
In the case where `sunpy.database.Database.fetch()` successfully downloads only some of the search results, a `~sunpy.database.PartialFetchError` is raised. This fixes a bug where the successful downloads would have been added to the database, but sometimes with incorrect metadata.
Updated the ``plot`` methods on some timeseries classes to correctly label and format the time axis.
Added the :meth:`~sunpy.net.base_client.QueryResponseTable.total_size`, which
estimates the total size of the results from a Fido query. If this is supported
by a client, the total size is printed alongside the results.

To add support for this in external clients, make sure one column contains
the individual filesizes as `~astropy.units.Quantity`, and set the
``size_column`` class attribute to the name of this column.
In consultation with JSOC, we now limit all JSOC downloads to one connection.
This will override all connection user settings passed to the downloader.
Removed support for Python 3.7.
Fixed a long-standing bug where our logger could intercept Astropy warnings in addition to SunPy warnings, and thus could conflict with Astropy's logger.
Added support for Python 3.10
Fixed the units of `sunpy.map.sources.HMISynopticMap.scale` and
`sunpy.map.sources.MDISynopticMap.scale`.
:func:`sunpy.map.make_fitswcs_header` now includes a PC_ij matrix in the returned
header if no rotation is specified.
Improved the gallery example :ref:`sphx_glr_generated_gallery_units_and_coordinates_SDO_to_STEREO_Coordinate_Conversion.py` to better illustrate how coordinate transformations interact with submaps and coordinate plotting.
Tidy the API Reference section of the documentation and improve the landing
page for the docs.
Add info about loading CDF files to the API documentation.
=========
Changelog
=========

.. note::

    This README was adapted from the pytest changelog readme under the terms of the MIT licence.

This directory contains "news fragments" which are short files that contain a small **ReST**-formatted text that will be added to the next ``CHANGELOG``.

The ``CHANGELOG`` will be read by users, so this description should be aimed at SunPy users instead of describing internal changes which are only relevant to the developers.

Make sure to use full sentences with correct case and punctuation, for example::

    Add support for Helioprojective coordinates in `sunpy.coordinates.frames`.

Please try to use Sphinx intersphinx using backticks.

Each file should be named like ``<PULL REQUEST>.<TYPE>[.<COUNTER>].rst``, where ``<PULL REQUEST>`` is a pull request number, ``COUNTER`` is an optional number if a PR needs multiple entries with the same type and ``<TYPE>`` is one of:

* ``breaking``: A change which requires users to change code and is not backwards compatible. (Not to be used for removal of deprecated features.)
* ``feature``: New user facing features and any new behavior.
* ``bugfix``: Fixes a reported bug.
* ``doc``: Documentation addition or improvement, like rewording an entire session or adding missing docs.
* ``deprecation``: Feature deprecation
* ``removal``: Feature removal.
* ``trivial``: A change which has no user facing effect or is tiny change.

So for example: ``123.feature.rst``, ``456.bugfix.rst``.

If you are unsure what pull request type to use, don't hesitate to ask in your PR.

Note that the ``towncrier`` tool will automatically reflow your text, so it will work best if you stick to a single paragraph, but multiple sentences and links are OK and encouraged.
You can install ``towncrier`` and then run ``towncrier --draft`` if you want to get a preview of how your change will look in the final release notes.
Added automatic conversion of some common but non-standard unit strings in CDF
files to astropy unit objects. If sunpy does not recognise the unit string for
a particular column, units of ``u.dimensionless_unscaled`` are applied to that
column and a warning raised.

If you think a given unit should not be dimensionless and support should be
added for it in sunpy, please raise an issue at
https://github.com/sunpy/sunpy/issues.
Sped up the parsing of results from the VSO. For large queries this significantly
reduces the time needed to perform a query to the VSO.
Added `packaging <https://pypi.org/project/packaging/>`__ as a core depedency as distutils is now deprecated.
************************
sunpy core Documentation
************************

The sunpy core package is a community-developed, free and open-source solar
data analysis environment for Python.

.. panels::

    Getting started
    ^^^^^^^^^^^^^^^
    .. toctree::
      :maxdepth: 1

      guide/installation
      guide/index
      generated/gallery/index
      code_ref/index
      whatsnew/index

    ---

    Other info
    ^^^^^^^^^^
    .. toctree::
      :maxdepth: 1

      about
      code_ref/known_issues
      code_ref/stability
      dev_guide/index
.. include:: ../CITATION.rst
.. doctest-skip-all

.. _whatsnew-0.9:

************************
What's New in SunPy 0.9?
************************

Overview
========

SunPy 0.9 brings improved support for downloading data from the JSOC and bugfixes compared to the 0.8.x series of releases.
The 0.9.x series will be the last series of SunPy releases to support Python 2.
This is because Python 2 `will not be maintained after 2019 <https://python3statement.org/>`_.
The 0.9.x series will receive bugfixs only up until the and of life of Python 2 (around 18 months).
No new functionality will be added to the 0.9.x series, which will also be the last version to include
``sunpy.spectra``, ``sunpy.lightcurve`` and ``sunpy.wcs``, all of which were deprecated in 0.8.

SunPy 1.0 and higher will support Python 3 only.
All new functionality will be available only in SunPy 1.0 and higher.

On this page, you can read about some of the big changes in this release:

* :ref:`whatsnew-0.9-python`
* :ref:`whatsnew-0.9-jsoc`
* :ref:`whatsnew-0.9-renamed-removed`

SunPy 0.9 includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 807 commits have been added since 0.8
* 310 issues have been closed since 0.8
* 147 pull requests have been merged since 0.8
* 34 people have contributed since 0.8
* 19 new contributors

There have been numerous improvements to large parts of SunPy, notably in the content of SunPy's documentation, and the continuous integration testing of SunPy.
In addition, SunPy co-ordinate system functionality has been improved, and the transformation between SunPy specific co-ordinate systems and those implemented by Astropy is now tested.
The `sunpy.map.CompositeMap` object has received numerous bugfixes, improving its functionality.
Bugfixes for various animation functionality have also been implemented.

.. _whatsnew-0.9-python:

Supported versions of Python
============================

SunPy is tested against Python 2.7, 3.5 and 3.6.

.. _whatsnew-0.9-jsoc:

Improvements to JSOC functionality
==================================

JSOC search capabilities have been improved.
It is now possible to search using any JSOC prime key, to search the metadata only, or to search by any series.
SunPy's JSOC functionality uses the DRMS library, which is now a SunPy affiliated package.
We would like to thank Kolja Glogowski for the DRMS package and for his help with the JSOC project.

When using a JSOC query in Fido, you must provide a JSOC series and at least one PrimeKey (for example a start time and end time)::

    >>> from sunpy.net import Fido, attrs as a

    >>> result = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Series('hmi.v_45s'), a.jsoc.Notify('me@email.org')))
    >>> result
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fca7050d6a0>
    Results from 1 Provider:

    81 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
             str23            str7     str10    float64   int64
    ----------------------- -------- ---------- -------- -------
    2014.01.01_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:01:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:02:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:04:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
                        ...      ...        ...      ...     ...
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145

    >>> result = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Notify('me@email.org'), a.jsoc.PrimeKey('WAVELNTH', '171'))

    >>> result
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fca2981e1d0>
    Results from 1 Provider:

    301 Results from the JSOCClient:
           T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
           str20           str7     str5    int64    int64
    -------------------- -------- -------- -------- -------
    2014-01-01T00:00:01Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:00:25Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:00:37Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:00:49Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:01:01Z  SDO/AIA    AIA_3      171    2145
                     ...      ...      ...      ...     ...
    2014-01-01T00:58:49Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:59:01Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:59:13Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:59:25Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:59:37Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:59:49Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_3      171    2145

Data is downloaded using::

    >>> files = Fido.fetch(result)

which returns a set of filepaths to the downloaded data.

For more information on accessing JSOC data using SunPy please `consult the documentation <https://docs.sunpy.org/guide/acquiring_data/jsoc.rst>`_.

.. _whatsnew-0.9-renamed-removed:

Renamed/removed functionality
=============================

sunpy.coordinates.representations
---------------------------------

The package ``sunpy.coordinates.representations`` has been removed.

Full change log
===============

To see a detailed list of all changes in version v0.9, including changes in API, please see the :ref:`changelog`.
.. doctest-skip-all

.. _whatsnew-0.8:

************************
What's New in SunPy 0.8?
************************

Overview
========

SunPy 0.8 is a major release that adds significant new functionality since the 0.7.x series of releases.

On this page, you can read about some of the big changes in this release:

* :ref:`whatsnew-0.8-python`
* :ref:`whatsnew-0.8-fido`
* :ref:`whatsnew-0.8-coordinates`
* :ref:`whatsnew-0.8-timeseries`
* :ref:`whatsnew-0.8-differential_rotation`
* :ref:`whatsnew-0.8-renamed-removed`

In addition to these major changes, SunPy 0.8 includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 1442 commits have been added since 0.7
* 163 issues have been closed since 0.7
* 200 pull requests have been merged since 0.7
* 35 distinct people have contributed code
* 17 new contributors

.. _whatsnew-0.8-python:

Supported versions of Python
============================

SunPy is tested against Python 2.7, 3.5 and 3.6.  SunPy no longer supports Python 3.4.

.. _whatsnew-0.8-fido:

New data downloading interface
==============================

The new package Fido is the primary interface to search for and download observational data.
Fido is a unified interface for searching and fetching solar physics data irrespective of the underlying client or webservice through which the data is obtained, e.g. VSO, JSOC, etc.
It therefore supplies a single, easy and consistent way to obtain most forms of solar physics data. It unifies data search and download from multiple sources such as the VSO and non-VSO sources::

    >>> from sunpy.net import Fido, attrs as a

    >>> attrs_time = a.Time('2005/01/01 00:10', '2005/01/01 00:15')
    >>> result = Fido.search(attrs_time, a.Instrument.eit)
    >>> result
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7f581489bb70>
    Results from 1 Provider:

    1 Results from the VSOClient:
        Start Time [1]       End Time [1]    Source Instrument   Type   Wavelength [2]
                                                                          Angstrom
        str19               str19         str4     str3      str8      float64
    ------------------- ------------------- ------ ---------- -------- --------------
    2005-01-01 00:12:10 2005-01-01 00:12:22   SOHO        EIT FULLDISK 195.0 .. 195.0

    >>> attrs_time = a.Time('2016/10/22 00:00', '2016/10/25 23:59')
    >>> result = Fido.search(attrs_time, a.Instrument.lyra)
    >>> result
    <sunpy.net.fido_factory.UnifiedResponse object at 0x11dfd8048>
    Results from 1 Provider:

    4 Results from the LYRAClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-10-22 00:00:00 2016-10-25 23:59:00 Proba2       lyra        nan
    2016-10-22 00:00:00 2016-10-25 23:59:00 Proba2       lyra        nan
    2016-10-22 00:00:00 2016-10-25 23:59:00 Proba2       lyra        nan
    2016-10-22 00:00:00 2016-10-25 23:59:00 Proba2       lyra        nan

Data is downloaded using::

    >>> files = Fido.fetch(result)

which returns a set of filepaths to the data.

.. _whatsnew-0.8-coordinates:

New coordinate package
======================

The SunPy coordinates package describes locations in physical space, and coordinate frames. It provides a way to transform coordinates from one frame like helioprojective to another such as heliographic.

Coordinates can be defined very simply::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import frames
    >>> a = SkyCoord(200*u.arcsec, 300*u.arcsec, frame=frames.Helioprojective)
    >>> a
    <SkyCoord (Helioprojective: obstime=None, rsun=695508.0 km, observer=earth): (Tx, Ty) in arcsec
        ( 200.,  300.)>
    >>> b = SkyCoord(-20*u.degree, -56*u.degree, frame=frames.HeliographicStonyhurst)
    >>> b
    <SkyCoord (HeliographicStonyhurst: obstime=None): (lon, lat, radius) in (deg, deg, km)
        (-20., -56.,  695508.)>


The coordinate ``a`` is in the helioprojective coordinate system, and the coordinate ``b`` is in the heliographic Stonyhurst system.

Maps can also be used to define coordinate frames::

    >>> import sunpy.map
    >>> import sunpy.data.sample
    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> c = SkyCoord(-45*u.arcsec, 600*u.arcsec, frame=aia_map.coordinate_frame)
    >>> c
    <SkyCoord (Helioprojective: obstime=2011-06-07 06:33:02.770000, rsun=696000000.0 m, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07 06:33:02.770000): (lon, lat, radius) in (deg, deg, m)
        ( 0.,  0.048591,   1.51846026e+11)>): (Tx, Ty) in arcsec
        (-45.,  600.)>

The coordinate ``c`` is now defined with respect to the coordinate frame derived from the map.
The observer attribute::

    >>> c.observer
    <HeliographicStonyhurst Coordinate (obstime=2011-06-07 06:33:02.770000): (lon, lat, radius) in (deg, deg, m)
        ( 0.,  0.048591,   1.51846026e+11)>

defines the location from which the coordinate was observed.

Transformation between solar physics coordinate systems
-------------------------------------------------------

Transformation between solar physics coordinate frames is simple::

    >>> c.transform_to(frames.HeliographicStonyhurst)
    <SkyCoord (HeliographicStonyhurst: obstime=2011-06-07 06:33:02.770000): (lon, lat, radius) in (deg, deg, km)
        (-3.51257477,  39.27459767,  696000.00000088)>

Transformation to astropy coordinate systems
--------------------------------------------

Solar physics coordinates can also be transformed into astrophysical coordinates.
For example, to convert to the International Celestial Reference System (ICRS)::

    >>> c.transform_to('icrs')
    <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, km)
        ( 224.85859731,  10.52568476,  998439.00599877)>

Specification of observer at any major solar system body
--------------------------------------------------------

Major solar system bodies can be used to specify observer locations in SkyCoord::

    >>> d = SkyCoord(-45*u.arcsec, 600*u.arcsec, observer='Mars', obstime='2011-06-07 06:33:02', frame=frames.Helioprojective)
    >>> d
    <SkyCoord (Helioprojective: obstime=2011-06-07 06:33:02, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07 06:33:02): (lon, lat, radius) in (deg, deg, AU)
        ( 135.78519602,  4.47598707,  1.43448427)>): (Tx, Ty) in arcsec
        (-45.,  600.)>


.. _whatsnew-0.8-timeseries:

New timeseries data object
==========================

The TimeSeries object is used to represent columns of time-ordered scalar values, and is source-aware, just like the Map object.
This object supersedes the LightCurve object, which is now deprecated in 0.8.

The TimeSeries object can be instantiated by passing in a file::

    >>> import sunpy.timeseries
    >>> import sunpy.data.sample
    >>> goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)

TimeSeries objects can have more than one column::

    >>> goes.columns
    ['xrsa', 'xrsb']

and have convenient plotting methods.

.. plot::
    :include-source:

    import sunpy.timeseries
    import sunpy.data.sample
    goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)
    goes.peek()

TimeSeries objects have a 'meta' property that stores the metadata of the timeseries::

    >>> goes.meta
    |-------------------------------------------------------------------------------------------------|
    |TimeRange                  | Columns         | Meta                                              |
    |-------------------------------------------------------------------------------------------------|
    |2011-06-06 23:59:59.961999 | xrsa            | simple: True                                      |
    |            to             | xrsb            | bitpix: 8                                         |
    |2011-06-07 23:59:57.631999 |                 | naxis: 0                                          |
    |                           |                 | extend: True                                      |
    |                           |                 | date: 26/06/2012                                  |
    |                           |                 | numext: 3                                         |
    |                           |                 | telescop: GOES 15                                 |
    |                           |                 | instrume: X-ray Detector                          |
    |                           |                 | object: Sun                                       |
    |                           |                 | origin: SDAC/GSFC                                 |
    |                           |                 | ...                                               |
    |-------------------------------------------------------------------------------------------------|

and the data can be accessed as a pandas dataframe using::

    >>> goes.data
                                        xrsa          xrsb
    2011-06-06 23:59:59.961999  1.000000e-09  1.887100e-07
    2011-06-07 00:00:02.008999  1.000000e-09  1.834600e-07
    2011-06-07 00:00:04.058999  1.000000e-09  1.860900e-07
    ...
    2011-06-07 23:59:55.584999  1.000000e-09  1.624800e-07
    2011-06-07 23:59:57.631999  1.000000e-09  1.598500e-07

Data sources that do not provide FITS files need to have a ``source`` keyword to help with the identification and interpretation of the data::

    >>> eve = sunpy.timeseries.TimeSeries(sunpy.data.sample.EVE_TIMESERIES, source='EVE')

.. _whatsnew-0.8-differential_rotation:

Differential rotation of maps
=============================

Maps can now be transformed using solar differential rates using ``from sunpy.physics.differential_rotation import diffrot_map``.

.. _whatsnew-0.8-renamed-removed:

Renamed/removed functionality
=============================

Several sub-packages have been moved or removed, and these are described in the following sections.

sunpy.lightcurve
----------------

The package ``sunpy.lightcurve`` has been deprecated in favor of `~sunpy.timeseries`, and will be removed in a future version of SunPy.

sunpy.physics.transforms
------------------------

The modules in ``sunpy.physics.transforms`` have been moved to `~sunpy.physics`.

sunpy.net
---------
``HelioviewerClient`` has been removed from the ``sunpy.net`` namespace. It
should now be imported with ``from sunpy.net.helioviewer import HelioviewerClient``.

Full change log
===============

To see a detailed list of all changes in version v0.8, including changes in API, please see the :ref:`changelog`.
.. doctest-skip-all

.. _whatsnew-1.0:

************************
What's New in SunPy 1.0?
************************

Overview
========

The SunPy project is pleased to announce the 1.0 release of the sunpy package.
This release is the result of over 14 months of work and marks a shift in the core library with the removal of a lot of legacy code and a focus on usability and stability.
The headline changes in 1.0 are:

* A complete transition of the whole code base to use `astropy.time.Time`, which was implemented by Vishnunarayan K I as part of Google Summer of Code 2018.
* A rewrite of how all the clients in `sunpy.net` download files from the internet.
  This means vastly improved progress bars, skipping downloads if files are present, and better visibility and retrying of failed downloads.
* A rewrite of the differential rotation and image warping code to correctly account for observer location using the Astropy coordinate functionality.
* Removal of many deprecated functions and submodules; we have used the 1.0 release as a chance to clean out SunPy, reducing the number of lines of Python code in the package by almost 3,000!
* The first release of SunPy to be Python 3 only, requiring Python 3.6+.  This is because Python 2 `will not be maintained after 2019 <https://python3statement.org/>`__.

On this page, you can read about some of the big changes in this release:

* :ref:`whatsnew-1.0-python`
* :ref:`whatsnew-1.0-time`
* :ref:`whatsnew-1.0-download`
* :ref:`whatsnew-1.0-logging`
* :ref:`whatsnew-1.0-coordinates`
* :ref:`whatsnew-1.0-diffrot`
* :ref:`whatsnew-1.0-maputils`
* :ref:`whatsnew-1.0-moved-config`
* :ref:`whatsnew-1.0-renamed-removed`

There have been numerous improvements to large parts of SunPy, notably in the content of the documentation, continuous integration and testing.
SunPy 1.0 also includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 1913 commits have been added since 0.9
* 582 issues have been closed since 0.9
* 332 pull requests have been merged since 0.9
* 46 people have contributed since 0.9
* 25 new contributors

.. _whatsnew-1.0-python:

Supported versions of Python
============================

SunPy 1.0 has dropped Python 2 support as Python 2 is nearing the end of its support cycle.
To port your code to Python 3, there is some good advice `here <https://docs.python.org/3/howto/pyporting.html>`__ and a
breakdown of the relevant changes in Python 3 on the `Python 3 for Scientists <https://python-3-for-scientists.readthedocs.io/en/latest/>`__
page. If you still need Python 2 support, SunPy 0.9.X will see bug fixes until the end of 2019, but no new features. As a result,
SunPy's minimum Python version has been raised to Python 3.6 and is routinely tested against Python 3.6 and 3.7.

.. _whatsnew-1.0-time:

Astropy Time is used everywhere
===============================

SunPy now uses Astropy's `~astropy.time.Time` object everywhere to represent time.
This comes with numerous benefits:

- **Support for non-UTC time scales.**
  UTC as well as non-UTC time scales like TAI, TT, UT1 etc. can be used with `astropy.time.Time`.

  .. code:: python

      >>> t = Time('2012-06-18T02:00:05.453', scale='tai')
      >>> t
      <Time object: scale='tai' format='isot' value=2012-06-18T02:00:05.453>

  `~astropy.time.Time` also provides easy conversion between different scales.

  .. code:: python

      >>> t.utc
      <Time object: scale='utc' format='isot' value=2012-06-18T01:59:31.453>

- **Support for high precision times.**
  `~astropy.time.Time` can provide sub-nanosecond precision for time objects while python
  `datetime` was restricted to microseconds.

  .. code:: python

    >>> t = Time('2012-06-18T02:00:05.453123123')
    >>> t
    <Time object: scale='utc' format='isot' value=2012-06-18T02:00:05.453>
    >>> t.precision = 9
    >>> t
    <Time object: scale='utc' format='isot' value=2012-06-18T02:00:05.453123123>

- **Support for leap seconds**
  This was one of the biggest motivations for the transition to `astropy.time.Time`.
  `datetime` has no support for leap seconds while `~astropy.time.Time` supports them.
  A leap second is a one-second adjustment applied to UTC to keep it close to the mean solar time.

  .. code:: python

    >>> Time('2016-12-31T23:59:60')
    <Time object: scale='utc' format='isot' value=2016-12-31T23:59:60.000>
    >>> Time('2016-12-31T23:59:59') + 1 * u.s
    <Time object: scale='utc' format='isot' value=2016-12-31T23:59:60.000>

- **Support for numerous formats**
  `~astropy.time.Time` can parse numerous formats including python `datetime`.

  .. code:: python

    >>> list(Time.FORMATS)
    ['jd', 'mjd', 'decimalyear', 'unix', 'cxcsec', 'gps', 'plot_date', 'datetime', 'iso', 'isot', 'yday', 'fits', 'byear', 'jyear', 'byear_str', 'jyear_str']

  .. code:: python

    >>> import datetime
    >>> Time(datetime.datetime.now())
    <Time object: scale='utc' format='datetime' value=2018-10-20 15:36:16.364089>

- **Changes in return values**

  All functions which previously returned `datetime.datetime` now return `~astropy.time.Time` and all functions which returned `datetime.timedelta` now return `astropy.time.TimeDelta`.
  For example, the properties of `sunpy.time.TimeRange` which  used to return `datetime.datetime` and `datetime.timedelta` now return `astropy.time.Time` and `astropy.time.TimeDelta`.

- **Changes to** `~sunpy.time.parse_time`

  `~sunpy.time.parse_time` has been reduced to a tiny wrapper over `~astropy.time.Time`.
  The API of `~sunpy.time.parse_time` is almost the same as `~astropy.time.Time`, however, `~sunpy.time.parse_time` supports conversion of a few more formats than `~astropy.time.Time`, which
  are `numpy.datetime64`, `pandas.Series`, `pandas.DatetimeIndex`, utime and a few other time string formats.


.. _whatsnew-1.0-download:

Improved file downloading capability
====================================

The file download capability has been re-written to use the `parfive package <https://github.com/Cadair/parfive>`__.
This brings more visually appealing and informative progress bars, better reporting of download errors and the ability to
re-download failed files.

.. image:: https://user-images.githubusercontent.com/1391051/50290236-ecf07b00-0462-11e9-80b2-9473c918802a.gif
   :alt: Parfive progress bars in a Jupyter Notebook

.. image:: https://user-images.githubusercontent.com/1391051/50290239-efeb6b80-0462-11e9-8b17-dfc05f8e2a57.gif
   :alt: Parfive progress bars in a terminal


It is possible to retry any downloads which fail with::

  >>> files = Fido.fetch(results)  # Some downloads fail
  >>> files = Fido.fetch(files)  # Retry the downloads which failed


.. _whatsnew-1.0-coordinates:

Improvements to coordinates functionality
=========================================

- **Accurate Sun-specific coordinates calculations**

  Sun-specific coordinates calculations have been grouped together in `sunpy.coordinates.sun`, and the underlying implementations have been re-written to use Astropy rather than approximate expressions.
  Nearly all of the returned values now match published values in the *Astronomical Almanac* to published precision (e.g., the hundredth of an arcsecond for apparent right ascension).
  For times that are provided to these functions, the user should take care to specify whether the time is other than UT (e.g., TT), which can be done using `~astropy.time.Time` (see above).

- **Improved tools to get positions of bodies in the solar system**

  The existing function `~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` has been enhanced to be able to correct for light travel time.
  When one specifies an observer, the function determines the emission time in the past that results in photons arriving at the observer at the observation time.
  The function then returns the location of the requested body at that emission time.

  .. code:: python

    >>> t = '2012-06-05 22:34:48.350'

    >>> without_correction = get_body_heliographic_stonyhurst('venus', t)
    >>> print(without_correction)
    <HeliographicStonyhurst Coordinate (obstime=2012-06-05T22:34:48.350): (lon, lat, radius) in (deg, deg, AU)
        (359.92620234, 0.02752007, 0.72602872)>

    >>> with_correction = get_body_heliographic_stonyhurst('venus', t, observer=get_earth(t))
    INFO: Apparent body location accounts for 144.06 seconds of light travel time [sunpy.coordinates.ephemeris]
    >>> print(with_correction)
    <HeliographicStonyhurst Coordinate (obstime=2012-06-05T22:34:48.350): (lon, lat, radius) in (deg, deg, AU)
        (359.92355609, 0.02734159, 0.72602853)>

  There is a new function `~sunpy.coordinates.ephemeris.get_horizons_coord` that queries `JPL HORIZONS <https://ssd.jpl.nasa.gov/?horizons>`__ for the location of solar-system bodies.
  JPL HORIZONS includes not only planets and other natural bodies in the solar system, but also major spacecraft.
  This function requires the `Astroquery <https://astroquery.readthedocs.io/en/latest/>`__ package and an Internet connection.


  - Query the location of Venus

  .. code:: python

    >>> get_horizons_coord('Venus barycenter', '2001-02-03 04:05:06')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for Venus Barycenter (2) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-03T04:05:06.000): (lon, lat, radius) in (deg, deg, AU)
        (326.06844114, -1.64998481, 0.71915147)>

  - Query the location of the SDO spacecraft

  .. code:: python

    >>> get_horizons_coord('SDO', '2011-11-11 11:11:11')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for Solar Dynamics Observatory (spac [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2011-11-11T11:11:11.000): (lon, lat, radius) in (deg, deg, AU)
        (0.01018888, 3.29640407, 0.99011042)>

  - Query the location of the SOHO spacecraft via its ID number (-21)

  .. code:: python

    >>> get_horizons_coord(-21, '2004-05-06 11:22:33', 'id')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for SOHO (spacecraft) (-21) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2004-05-06T11:22:33.000): (lon, lat, radius) in (deg, deg, AU)
        (0.2523461, -3.55863351, 0.99923086)>


.. _whatsnew-1.0-logging:

Logging used to record SunPy notices
====================================

All messages provided by SunPy use a new logging facility which is based on the Python logging module rather than print statements.

Messages can have one of several levels, in increasing order of importance:

* DEBUG: Detailed information, typically of interest only when diagnosing problems.
* INFO: A message conveying information about the current task, and confirming that things are working as expected.
* WARNING: An indication that something unexpected happened, and that user action may be required.
* ERROR: An indication that a more serious issue has occured, where something failed but the task is continuing.
* CRITICAL: A serious error, indicating that the program itself may be unable to continue running.

By default, all messages except for DEBUG messages are displayed. Messages can also be sent to a file and time stamped.

See the :ref:`logger` documentation for instructions on how to control the verbosity of the logger.


.. _whatsnew-1.0-diffrot:

Improvements to differential rotation
=====================================

Applying the effect of solar differential rotation to coordinates now properly takes into account the changing position of the observer.
For example, since the Earth moves, observers on the Earth must take into account the solar differential rotation of the Sun and the motion of the Earth when calculating a location on the Sun.

- **Support for applying solar differential rotation to coordinates.**

  Solar differential rotation of on-disk coordinates can be specified using either time or a new observer.
  If time is specified, then the new observer is assumed to be located on the Earth::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import Helioprojective
    >>> from sunpy.physics.differential_rotation import solar_rotate_coordinate
    >>> from sunpy.time import parse_time

    >>> start_time = '2010-09-10 12:34:56'
    >>> duration = 25*u.hour
    >>> c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=start_time, frame=Helioprojective)
    >>> solar_rotate_coordinate(c, time=duration)
    <SkyCoord (Helioprojective: obstime=2010-09-11T13:34:56.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-09-11T13:34:56.000): (lon, lat, radius) in (deg, deg, AU)
        (-5.08888749e-14, 7.24318962, 1.00669016)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (-363.04027419, 104.87807178, 1.499598e+08)>

  Due to the ellipticity of the Earth's orbit, the amount of solar rotation is different at different times in the year::

    >>> start_time = '2010-06-10 12:34:56'
    >>> duration = 25*u.hour
    >>> c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=start_time, frame=Helioprojective)
    >>> solar_rotate_coordinate(c, time=duration)
    <SkyCoord (Helioprojective: obstime=2010-06-10T12:34:56.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-06-11T13:34:56.000): (lon, lat, radius) in (deg, deg, AU)
        (0., 0.58398742, 1.01539908)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (-359.11576773, 117.18020622, 1.51263627e+08)>


  The user can also explicitly specify an observer at a different time and location in space.  The amount of solar
  rotation applied depends on the time difference between the observation time of the`~astropy.coordinates.SkyCoord`
  and the time of the observer::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import Helioprojective, HeliographicStonyhurst
    >>> from sunpy.physics.differential_rotation import solar_rotate_coordinate
    >>> from sunpy.time import parse_time

    >>> start_time = parse_time('2010-06-10 12:34:56')
    >>> duration = 25*u.hour
    >>> c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=start_time, frame=Helioprojective)
    >>> new_observer = SkyCoord(lon=20*u.deg, lat=8*u.deg, radius=0.9*u.au, obstime=end_time, frame=HeliographicStonyhurst)
    >>> solar_rotate_coordinate(c, observer=new_observer)
    <SkyCoord (Helioprojective: obstime=2010-06-10T12:34:56.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-06-11T13:34:56.000): (lon, lat, radius) in (deg, deg, AU)
        (20., 8., 0.9)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (-715.77862011, 31.87928146, 1.34122226e+08)>

- **Experimental support for applying solar differential rotation to maps.**

  Applying solar differential rotation to maps also accounts for changing observer position.
  This functionality is still experimental.
  For example, to differentially rotate a map back 23 hours::

  >>> import astropy.units as u
  >>> import sunpy.map
  >>> from sunpy.data.sample import AIA_171_IMAGE
  >>> from sunpy.physics.differential_rotation import differential_rotate

  >>> aia = sunpy.map.Map(AIA_171_IMAGE)
  >>> differential_rotate(aia, time=-23*u.hour)

  `~sunpy.physics.differential_rotation.differential_rotate` also accepts a new observer keyword.
  The amount of solar differential rotation is calculated using the time difference between the map date and observation time of the new observer.
  For example::

  >>> import astropy.units as u
  >>> import sunpy.map
  >>> from sunpy.data.sample import AIA_171_IMAGE

  >>> from sunpy.physics.differential_rotation import differential_rotate
  >>> aia = sunpy.map.Map(AIA_171_IMAGE)
  >>> new_observer = SkyCoord(lon=-15*u.deg, lat=-4*u.deg, radius=1*u.au, obstime=aia.date-34*u.hour, frame=HeliographicStonyhurst)
  >>> differential_rotate(aia, observer=new_observer)

.. _whatsnew-1.0-maputils:

Map utility functions
=====================

A set of new utility functions have been added to `sunpy.map` which act on `sunpy.map.GenericMap` instances.

For example, getting the world coordinates for every pixel::

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE
    >>> import astropy.units as u
    >>> from sunpy.physics.differential_rotation import differential_rotate
    >>> from sunpy.map import contains_full_disk, all_coordinates_from_map

    >>> aia = sunpy.map.Map(AIA_171_IMAGE)
    >>> contains_full_disk(aia)
    True
    >>> coordinates = all_coordinates_from_map(aia) # The coordinates for every map pixel
    >>> coordinates.shape
    (1024, 1024)


or generating a new FITS header for a custom map::

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from sunpy.coordinates import frames
    >>> from astropy.coordinates import SkyCoord

    >>> data = np.arange(0,100).reshape(10,10)
    >>> coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime = '2013-10-28', observer = 'earth', frame = frames.Helioprojective)
    >>> header = sunpy.map.header_helper.make_fitswcs_header(data, coord)
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 5.5
    crpix2: 5.5
    cdelt1: 1.0
    cdelt2: 1.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    lonpole: 180.0
    latpole: 0.0
    date-obs: 2013-10-28T00:00:00.000
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    dsun_obs: 148644585949.49176
    rsun_ref: 695700.0
    rsun_obs: 965.3723815059902


.. _whatsnew-1.0-moved-config:

Config File Location Moved
==========================

If you have customised your :ref:`customizing-with-sunpyrc-files` you will need to move it to the new config file location.
Your old file should be in ``~/.sunpy/sunpyrc`` file and the new location, which is now platform specific, can be found by running `sunpy.print_config`.
We recommand that your take a look at the new file as available configuration options have increased.

.. _whatsnew-1.0-renamed-removed:

Renamed/removed functionality
=============================

This is just some of the renamed or removed functionality.

* ``sunpy.sun.sun`` functions have been re-implemented using Astropy for significantly improved accuracy and moved to `sunpy.coordinates.sun`.
* Removed ``sunpy.time.julian_day``, ``sunpy.time.julian_centuries``, ``sunpy.time.day_of_year``, ``sunpy.time.break_time``, ``sunpy.time.get_day``.
* Move the matplotlib animators from ``sunpy.visualisation.imageanimator`` and ``sunpy.visualization.mapcubeanimator`` to `sunpy.visualization.animator`.
* ``axis_ranges`` kwarg of ``sunpy.visualization.animator.ArrayAnimator``,
  ``sunpy.visualization.animator.ImageAnimator`` and
  ``sunpy.visualization.animator.LineAnimator`` now must be entered as ``None``,
  ``[min, max]`` or pixel edges of each array element.
* The Helioviewer client has been switched to using the newer Helioviewer API. This has meant that we have changed some of the keywords that were passed into client's methods.
* Removed ``sunpy.net.jsoc.attrs.Time`` because it served the same purpose as ``sunpy.net.attrs.Time`` after the switch to ``astropy.time.Time``.
* The deprecated ``sunpy.lightcurve`` (replaced by `sunpy.timeseries`), ``sunpy.wcs`` and ``sunpy.spectra`` (replaced by the ``radiospectra`` package) modules have now been removed.

Full change log
===============

To see a detailed list of all changes in version v1.0, including changes in API, please see the :ref:`changelog`.
.. doctest-skip-all

.. _whatsnew-2.1:

************************
What's New in SunPy 2.1?
************************

Overview
========
The SunPy project is pleased to announce the 2.1 release of the sunpy package.

On this page, you can read about some of the big changes in this release:

* :ref:`whatsnew-2.1-goes1617`
* :ref:`whatsnew-2.1-metadata`
* :ref:`whatsnew-2.1-refactor`
* :ref:`whatsnew-2.1-synoptic`
* :ref:`whatsnew-2.1-cutouts`
* :ref:`whatsnew-2.1-velocities`
* :ref:`whatsnew-2.1-reprojecting`
* :ref:`whatsnew-2.1-contours`
* :ref:`whatsnew-2.1-performance`
* :ref:`whatsnew-2.1-python`
* :ref:`whatsnew-2.1-contributors`

SunPy 2.1 also includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 1148 commits have been added since 2.0
* 202 issues have been closed since 2.0
* 306 pull requests have been merged since 2.0
* 38 people have contributed since 2.0
* 19 of which are new contributors

Please find below a selection of what we consider to be the biggest changes or features with this release.

.. _whatsnew-2.1-goes1617:

Support for GOES-16 and GOES-17 X-Ray Sensor (XRS) data
=======================================================
We have now included support for the GOES-16 and GOES-17 XRS 1-second data. This data can now be queried and downloaded with `sunpy.net.Fido` and the files read in and analysed as a `sunpy.timeseries.TimeSeries`.

.. minigallery:: sunpy.net.attrs.goes.SatelliteNumber

.. _whatsnew-2.1-metadata:

Fido Metadata
=============

Fido provides a unified interface to download data from several services but lacked support for search for metadata provided by catalogues or data providers.
This has now changed and you can query JSOC, HEK and HEC with Fido::

    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> results = Fido.search(timerange,
                      a.hek.FL & (a.hek.FL.PeakFlux > 1000) |
                      a.jsoc.Series('hmi.m_45s'))
    >>> results
    Results from 2 Providers:

    2 Results from the HEKClient:
                                                                                                                        gs_thumburl                                                                                                                       ...
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ...
    http://sdowww.lmsal.com/sdomedia/ssw/ssw_client/data/ssw_service_100731_205448_25495028/www/EDS_FlareDetective-TriggerModule_20100801T033001-20100801T035225_AIA_171_S21W87_ssw_cutout_20100801_033013_AIA_171_S21W87_20100801_033012_context_0180.gif ...
    http://sdowww.lmsal.com/sdomedia/ssw/ssw_client/data/ssw_service_100801_234037_25860951/www/EDS_FlareDetective-TriggerModule_20100801T033008-20100801T035232_AIA_193_S21W87_ssw_cutout_20100801_033020_AIA_193_S21W87_20100801_033019_context_0180.gif ...

    1 Results from the JSOCClient:
            T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2010.08.01_03:40:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2099

Alongside this, there is easier access to the results from a client by indexing the results by the client name::

    >>> hek_results, jsoc_results = results['hek'], results['jsoc']

In the first example you can see the search results from HEK are not too useful by default.
The results are truncated as the long urls for the thumbnail take up too many characters and HEK can return up to 100 columns of information.
With a new method ``show()`` on the results itself, you directly select the columns to display::

    >>> hek_results.show('event_peaktime', 'obs_instrument', 'fl_peakflux')
    event_peaktime   obs_instrument fl_peakflux
    ------------------- -------------- -----------
    2010-08-01T03:40:37            AIA     1027.64
    2010-08-01T03:40:44            AIA     1441.78

To find the keys for all of the columns you can do::

    >>> hek_results.keys()
    ['gs_thumburl',
    'comment_count',
    'hpc_bbox',
    ...
    'area_unit',
    'obs_lastprocessingdate',
    'refs']

In addition, if you pass the entire search query into ``Fido.fetch()`` and it will ignore results that have no corresponding data files to retrieve.

.. minigallery:: sunpy.net.attrs.hek.FL sunpy.net.attrs.hek.EventType

.. _whatsnew-2.1-refactor:

Fido Results Refactor
=====================

While working on the above change, several changes were made to the way that results are returned from ``sunpy.net.Fido.search`` and the search methods of the underlying clients.

We have tried to minimize any breaking changes here and we believe that for most users the difference between 2.0 and 2.1 will be minor.

The key highlights you will want to be aware of are:

* Previously slicing the result of ``Fido.search()`` (a `~sunpy.net.fido_factory.UnifiedResponse` object) so that it had a length of one returned another `~sunpy.net.fido_factory.UnifiedResponse` object, now it will return a `~sunpy.net.base_client.QueryResponseTable` object, which is a subclass of `astropy.table.Table`.

* All result objects contained within the results of a ``Fido.search()`` are now `~sunpy.net.base_client.QueryResponseTable` objects (or subclasses thereof).
  These objects are subclasses of `astropy.table.Table` and can therefore be filtered and inspected as tabular objects, and the modified tables can be passed to ``Fido.fetch()``.

* The keys available to be used when formatting the ``path=`` argument to ``Fido.fetch()`` have changed. This is to standardise them over the results from more clients and make them easier to use.
  You can use the ``~.UnifiedResponse.path_format_keys`` method to see all the possible keys for a particular search.

* The search results object returned from ``Fido.search`` now correctly counts all results in its `~sunpy.net.fido_factory.UnifiedResponse.file_num` property.

* Results from the `~sunpy.net.dataretriever.NOAAIndicesClient` and the `~sunpy.net.dataretriever.NOAAPredictClient` no longer have ``Start Time`` or ``End Time`` in their results table as the results returned from the client are not dependent upon the time parameter of a search.

.. _whatsnew-2.1-synoptic:

New synoptic map sources and clients
====================================
`~sunpy.map.sources.MDISynopticMap` and `~sunpy.map.sources.HMISynopticMap` have been added as new data sources, and automatically fix common issues with FITS metadata from these sources. You do not need to change any code to use these, as sunpy automatically detects and uses the appropriate map sources for each file.

It is now possible to search for GONG synoptic maps within `sunpy.net.Fido`, using ``a.Instrument('GONG')``.

.. _whatsnew-2.1-cutouts:

Requesting cutouts from the JSOC
================================
`~sunpy.net.Fido` can now be used to request cutouts from JSOC via the new ``a.jsoc.Cutout`` attr. This includes the ability to adjust the requested field of view to "track" a feature as it moves across the solar disk, perform sub-pixel image registration, and mask off-disk pixels.

.. minigallery:: sunpy.net.attrs.jsoc.Cutout

.. _whatsnew-2.1-velocities:

Coordinates with velocities
===========================
It is now supported to transform coordinates with attached velocities, and the various ephemeris functions can optionally include velocity information.
Transformations between coordinate frames will account for both any change in orientation of the velocity vector and any induced velocity due to relative motion between the frames.
For example, consider Mars's position/velocity in `~sunpy.coordinates.frames.HeliographicStonyhurst`::

    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import get_body_heliographic_stonyhurst
    >>> mars = SkyCoord(get_body_heliographic_stonyhurst('mars', '2021-01-01',
    ...                                                  include_velocity=True))
    >>> mars
    <SkyCoord( HeliographicStonyhurst: obstime=2021-01-01T00:00:00.000): (lon, lat, radius) in (deg, deg, AU)
        (-34.46752135, 1.77496469, 1.50936573)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (-0.00048971, 0.00060976, 19.54950062)>
    >>> mars.velocity.norm()
    <Quantity 19.56823076 km / s>

However, `~sunpy.coordinates.frames.HeliographicStonyhurst` is a non-inertial frame that rotates over time.
By transforming this coordinate to `~sunpy.coordinates.frames.HeliocentricInertial`, we can see that Mars's actual velocity is larger::

    >>> mars.heliocentricinertial
    <SkyCoord (HeliocentricInertial: obstime=2021-01-01T00:00:00.000): (lon, lat, distance) in (deg, deg, AU)
        (-9.91128592, 1.77496469, 1.50936573)
     (d_lon, d_lat, d_distance) in (arcsec / s, arcsec / s, km / s)
        (0.04174239, 0.00060976, 19.54950058)>
    >>> mars.heliocentricinertial.velocity.norm()
    <Quantity 49.68592218 km / s>

See :ref:`sunpy-coordinates-velocities` for more information.

.. _whatsnew-2.1-reprojecting:

Alternatives for reprojecting a Helioprojective map
===================================================
The typical observation in `~sunpy.coordinates.frames.Helioprojective` coordinates does not contain full 3D information for the sources of emission, so an assumption needs to be made when transforming such coordinates to other coordinate frames.
By default, SunPy assumes that the emission is coming from the surface of the Sun, which enables reprojections such as in the example :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_different_observers.py`.
However, this assumption is not appropriate for some observations, e.g., from coronagraphs.

There is now a context manager (:meth:`~sunpy.coordinates.frames.Helioprojective.assume_spherical_screen`) to override the default assumption such that any 2D coordinates are interpreted as being on the inside of a large spherical screen.
See the following example for how this context manager enables alternative reprojections.

.. minigallery:: sunpy.coordinates.Helioprojective.assume_spherical_screen

.. _whatsnew-2.1-contours:

Finding map contours
====================
The new :meth:`sunpy.map.GenericMap.contour` method can be used to extract contours from a map. It returns contours as a `~astropy.coordinates.SkyCoord`, allowing contours to be easily overplotted on the original or other maps.

.. minigallery:: sunpy.map.GenericMap.contour

.. _whatsnew-2.1-performance:

Performance improvements
========================
Several functions in `sunpy.map` have been significantly sped up with improved algorithms.

In addition, `sunpy.map.GenericMap.wcs` is now cached when the map metadata remains unchanged, significantly improving performance in applications which make mutiple requests for the map WCS (e.g. plotting), and reducing the number of repeated warnings thrown when metadata is missing.

.. _whatsnew-2.1-python:

Increase in required package versions
=====================================
We have bumped the minimum version of several packages we depend on; these are the new minimum versions for sunpy 2.1:

- python 3.7
- astropy 4.0
- scipy 1.2
- parfive 1.1
- drms 0.6.1
- matplotlib 2.2.2

.. _whatsnew-2.1-contributors:

Contributors to this Release
============================

The people who have contributed to the code for this release are:

-  Abhijeet Manhas
-  Abhishek Pandey  *
-  Adrian Price-Whelan
-  Albert Y. Shih
-  Aryan Chouhan  *
-  Conor MacBride  *
-  Daniel Ryan
-  David Pérez-Suárez
-  David Stansby
-  Dipanshu Verma  *
-  Erik Tollerud  *
-  Jai Ram Rideout
-  Jeffrey Aaron Paul  *
-  Johan L. Freiherr von Forstner  *
-  Kateryna Ivashkiv  *
-  Koustav Ghosh  *
-  Kris Akira Stern
-  Kritika Ranjan  *
-  Laura Hayes
-  Lazar Zivadinovic
-  Nabil Freij
-  Rutuja Surve
-  Sashank Mishra
-  Shane Maloney
-  Shubham Jain  *
-  SophieLemos  *
-  Steven Christe
-  Stuart Mumford
-  Sudeep Sidhu  *
-  Tathagata Paul  *
-  Thomas A Caswell  *
-  Will Barnes
-  honey
-  mridulpandey
-  nakul-shahdadpuri  *
-  platipo  *
-  resakra  *
-  sophielemos  *

Where a * indicates that this release contains their first contribution to SunPy.
.. doctest-skip-all

.. _whatsnew-3.1:

************************
What's New in SunPy 3.1?
************************
The SunPy project is pleased to announce the 3.1 release of the sunpy core package.

On this page, you can read about some of the big changes in this release.

.. contents::
    :local:
    :depth: 1

SunPy 3.1 also includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 610 commits have been added since 3.0
* 56 issues have been closed since 3.0
* 145 pull requests have been merged since 3.0
* 16 people have contributed since 3.0
* 8 of which are new contributors

Increase in required package versions
=====================================
We have bumped the minimum version of several packages we depend on; these are the new minimum versions for sunpy 3.1:

- astropy >= 4.2
- matplotlib >= 3.2.0
- numpy >= 1.17.0
- pandas >= 1.0.0


Increased in-situ data support
==============================
Two new features have significantly increased support for in-situ data returned by heliospheric missions.
See :ref:`sphx_glr_generated_gallery_acquiring_data_search_cdaweb.py` for a full example of searching for, downloading, and loading a CDF file into sunpy.

New CDAWeb client
-----------------
A new Fido client to search the Coordinated Data Analysis Web (CDAWeb) has been added to `sunpy.net`.
This allows one to search for and download data from many space physics missions that take in-situ data (including Parker Solar Probe and Solar Orbiter).
In combination with the new CDF reading abilities of `sunpy.timeseries.TimeSeries`, this provides a full workflow for searching for, downloading, and analysing in-situ data contained within CDF files.

Support for .cdf files
----------------------
`sunpy.timeseries.TimeSeries` now supports loading CDF files if the external library ``cdflib`` is installed.

New limb drawing function
=========================
The solar limb as seen from an arbitrary observer coordinate can now be drawn on a world coordinate system aware
Axes using the :func:`sunpy.visualization.draw_limb` function.

.. minigallery:: sunpy.visualization.draw_limb


New WISPR map source
====================
A new map source for the WISPR instrument on Parker Solar Probe has been added.
This improves the `~sunpy.map.GenericMap.name` of the map and adds correct
information for the `~sunpy.map.GenericMap.processing_level` and
`~sunpy.map.GenericMap.exposure_time`.

Changes to map metadata fixes
=============================
The `~sunpy.map.GenericMap` map sources are primarily used to modify metadata values to either fix known incorrect values in a given data source, or make the metadata comply with the FITS standard.
There are two ways in which the sunpy map sources can do this:

1. Directly edit the values stored in the FITS metadata.
2. Overload properties (e.g. `sunpy.map.GenericMap.unit`) such that they return the correct information without modifying the underlying metadata.

In previous versions of sunpy there has been no consistency in which of these two approaches is taken.

As of sunpy 3.1, the second approach is now consistently taken, and **sunpy no longer edits any FITS metadata values** when constructing a Map.
`~sunpy.map.GenericMap` properties now consistently provide corrected information.
For details of any corrections applied, the docstrings of different map sources (e.g., `~sunpy.map.sources.HMIMap` - see :ref:`map-sources` for a list of map sources) provide information on any assumptions made beyond the original FITS metadata when constructing the map properties.

For all maps, the following fixes are no longer made:

- DATE-OBS is no longer replaced by DATE_OBS as a fallback
- NAXIS, NAXIS1, NAXIS2, BITPIX are no longer populated if not present
- BUNIT is no longer corrected to be a FITS compliant unit string
- WAVEUNIT  is no longer automatically populated from the header comments if it is not present.

For specific map sources, the following keywords are no longer modified or added:

- `~sunpy.map.sources.KCorMap`: OBSERVATORY, DETECTOR, WAVEUNIT, DSUN_OBS, HGLN_OBS
- `~sunpy.map.sources.SWAPMap`: OBSRVTRY, DETECTOR
- `~sunpy.map.sources.RHESSIMap`: CUNIT1, CUNIT2, CTYPE1, CTYPE2, WAVEUNIT, WAVELNTH
- `~sunpy.map.sources.AIAMap`: BUNIT, DETECTOR
- `~sunpy.map.sources.HMIMap`: DETECTOR, CRDER1, CRDER2
- `~sunpy.map.sources.HMISynopticMap`: CUNIT1, CUNIT2, CDELT1, CDELT2, DATE-OBS
- `~sunpy.map.sources.EITMap`: WAVEUNIT, CUNIT1, CUNIT2
- `~sunpy.map.sources.LASCOMap`: DATE-OBS, DATE_OBS, CROTA, CROTA1, CROTA2, CUNIT1, CUNIT2
- `~sunpy.map.sources.MDIMap`: CUNIT1, CUNIT2
- `~sunpy.map.sources.MDISynopticMap`: CUNIT1, CUNIT2, CDELT2, DATE-OBS, CRDER1, CRDER2
- `~sunpy.map.sources.EUVIMap`: WAVEUNIT, DATE-OBS, CROTA, CROTA2
- `~sunpy.map.sources.CORMap`: DATE-OBS
- `~sunpy.map.sources.HIMap`: DATE-OBS
- `~sunpy.map.sources.SUVIMap`: DETECTOR, TELESCOP
- `~sunpy.map.sources.TRACEMap`: DETECTOR, OBSRVTRY, CUNIT1, CUNIT2
- `~sunpy.map.sources.SXTMap`: DETECTOR, TELESCOP, DSUN_APPARENT
- `~sunpy.map.sources.XRTMap`: DETECTOR, TELESCOP, TIMESYS
- `~sunpy.map.sources.SOTMap`: DETECTOR, TELESCOP
- `~sunpy.map.sources.SJIMap`: DETECTOR, WAVEUNIT, WAVELNTH, CUNIT1, CUNIT2
- `~sunpy.map.sources.EUIMap`: CROTA, CROTA2

Changes to map date/time handling
=================================

New date properties
-------------------
The properties `~sunpy.map.GenericMap.date_start`,
`~sunpy.map.GenericMap.date_end`, and `~sunpy.map.GenericMap.date_average` have
been added to be drawn from the relevant FITS metadata, if present in the map
header. These are from new keywords defined in version 4 of the FITS standard,
which have precise meanings compared to the previously ill-defined DATE-OBS.

Changes to `~sunpy.map.GenericMap.date`
---------------------------------------
`sunpy.map.GenericMap.date` now looks for more metadata than just DATE-OBS.
This property can return any one of the new properties (see above) depending
on the metadata present in the map. It now draws from, in order of preference:

1. The DATE-OBS FITS keyword
2. `~sunpy.map.GenericMap.date_average`
3. `~sunpy.map.GenericMap.date_start`
4. `~sunpy.map.GenericMap.date_end`
5. The current time.

If DATE-OBS is present alongside DATE-AVG or DATE-BEG and DATE-END, this results
in a behaviour change to favour the new (more precisely defined) keywords.
It is recommended
to use `~sunpy.map.GenericMap.date_average`,
`~sunpy.map.GenericMap.date_start`, or `~sunpy.map.GenericMap.date_end`
instead if you need one of these specific times.

Addition of new time format `~sunpy.time.TimeTaiSeconds`
--------------------------------------------------------
The new `~sunpy.time.TimeTaiSeconds` format is the number of
SI seconds from 1958-01-01 00:00:00, which includes UTC leap seconds.
1958-01-01 00:00:00 is the defined time when International Atomic Time (TAI)
and Universal Time (UT) are synchronized.

This format is equivalent to the output of the SSW ``anytim2tai`` routine, and
related routines, for times after 1972-01-01.  Be aware that the SSW routines
are not written to provide valid results for times before 1972-01-01.

This format is equivalent to `~astropy.time.TimeUnixTai`, except that the epoch
is 12 years earlier.

Propagating solar-surface coordinates in time
=============================================
There is now an easy-to-use context manager (:func:`~sunpy.coordinates.propagate_with_solar_surface`) to enable coordinate transformations to take solar rotation into account.
Normally, a coordinate refers to a point in inertial space, so transforming it to a different observation time does not move the point at all.
Under this context manager, a coordinate will be treated as if it were referring to a point on the solar surface.
Coordinate transformations with a change in observation time will automatically rotate the point in heliographic longitude for the time difference, with the amount of rotation depending on the specified differential-rotation model.

.. minigallery:: sunpy.coordinates.propagate_with_solar_surface


Convenient reprojection of maps
===============================
`~sunpy.map.Map` objects now have the :meth:`~sunpy.map.GenericMap.reproject_to` method to easily reproject the map to a new WCS.
The returned map will be of type `~sunpy.map.GenericMap`, with no metadata preserved from the original map, so copy over any desired metadata from the original map.
This method requires the optional package `reproject` to be installed.

.. minigallery:: sunpy.map.GenericMap.reproject_to

JSOC keyword filtering with Fido
================================
Support for filtering searches with JSOC keywords has been added to ``Fido.search``::

    >>> from sunpy.net import Fido, attrs as a
    >>> import astropy.units as u
    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
        a.jsoc.Series('aia.lev1_euv_12s'), a.Wavelength(304*u.AA), a.jsoc.Keyword("EXPTIME") > 1)
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fe16a5d20d0>
    Results from 1 Provider:

    301 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu

        T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
    -------------------- -------- -------- -------- -------
    2014-01-01T00:00:01Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:00:25Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:00:37Z  SDO/AIA    AIA_4      304    2145
                    ...      ...      ...      ...     ...
    2014-01-01T00:59:25Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:59:37Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:59:49Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_4      304    2145
    Length = 301 rows
    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
        a.jsoc.Series('aia.lev1_euv_12s'), a.Wavelength(304*u.AA), a.jsoc.Keyword("EXPTIME") == 1)
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fe16a5d20d0>
    Results from 1 Provider:

    0 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu

Please be aware of two caveats:

- We do not validate the value used for comparison.
- Passing in a keyword without comparison to a value (e.g. ``==0``, ``< 10``) will error.


Arithmetic operations with maps
===============================

`~sunpy.map.GenericMap` objects now support arithmetic operations (i.e. addition, subtraction, multiplication, division) with array-like quantities.
This includes scalar quantities as well as Numpy arrays.
Notably, arithmetic operations between two `~sunpy.map.GenericMap` objects are not supported.

Contributors to this Release
============================

The people who have contributed to the code for this release are:

-  Alasdair Wilson  *
-  Albert Y. Shih
-  Anubhav Sinha  *
-  Conor MacBride
-  David Stansby
-  Devansh Shukla  *
-  Jeffrey Aaron Paul
-  Nabil Freij
-  Noah Altunian  *
-  Rohan Sharma  *
-  Samriddhi Agarwal
-  Stuart Mumford
-  Thomas Braccia  *
-  Tim Gates  *
-  Will Barnes

Where a * indicates that this release contains their first contribution to SunPy.
.. doctest-skip-all

.. _whatsnew-2.0:

************************
What's New in SunPy 2.0?
************************

Overview
========

The SunPy project is pleased to announce the 2.0 release of the sunpy package.
On this page, you can read about some of the big changes in this release:

* :ref:`whatsnew-2.0-python`
* :ref:`whatsnew-2.0-search`
* :ref:`whatsnew-2.0-aiaprep`
* :ref:`whatsnew-2.0-pixelindex`
* :ref:`whatsnew-2.0-standards`
* :ref:`whatsnew-2.0-overview`
* :ref:`whatsnew-2.0-coordframe`
* :ref:`whatsnew-2.0-carrington`
* :ref:`whatsnew-2.0-proxy`
* :ref:`whatsnew-2.0-citation`

SunPy 2.0 also includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 1044 commits have been added since 1.1
* 144 issues have been closed since 1.1
* 290 pull requests have been merged since 1.1
* 33 people have contributed since 1.1
* 16 new contributors

.. _whatsnew-2.0-python:

Increase in required package versions
=====================================

We have bumped the minimum version of several packages we depend on:

* numpy>=1.15.0
* scipy>=1.0.0
* matplotlib>=2.2.2
* astropy>=3.2
* parfive>=1.1.0

.. _whatsnew-2.0-search:

Search Attributes
=================

To search with `~sunpy.net.Fido`, you need to specify attributes to search against.
Before sunpy 2.0, you had to supply values, as the following example demonstrates::

    >>> from sunpy.net import Fido, attrs as a
    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument("norh"),
    ...             a.Wavelength(17*u.GHz))

There was no way to know if the value was correct, but now we have a extenstive list of supported values from the clients and servers we can request data from.

Using `~sunpy.net.attrs.Instrument` as an example, if you print the object::

    >>> print(a.Instrument)
    sunpy.net.attrs.Instrument
    <BLANKLINE>
    Specifies the Instrument name for the search.
    <BLANKLINE>
           Attribute Name          Client          Full Name                                           Description
    --------------------------- ----------- ------------------------ --------------------------------------------------------------------------------
    aia                         VSO         AIA                      Atmospheric Imaging Assembly
    bbi                         VSO         BBI                      None
    bcs                         VSO         BCS                      Bragg Crystal Spectrometer
    bic_hifi                    VSO         BIC-HIFI                 None
    bigbear                     VSO         Big Bear                 Big Bear Solar Observatory, California TON and GONG+ sites
    ...

This will list the name of value you should use, what data source will supply that data and a description.
Furthermore, you can use tab completion to auto-fill the attribute name, for example by typing ``a.Instrument.<TAB>``.

So now you can do the following instead::

    Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.norh, a.Wavelength(17*u.GHz))


.. _whatsnew-2.0-aiaprep:

aiaprep is now deprecated
=========================

With the release of the new `aiapy <https://aiapy.readthedocs.io>`__ package, ``sunpy.instr.aia.aiaprep`` will be removed in version 2.1.
Equivalent functionality is provided by the ``register()`` function in aiapy. For more
details, see the `example on registering and aligning level 1 AIA images <https://aiapy.readthedocs.io/en/latest/generated/gallery/prepping_level_1_data.html>`_
in the aiapy documentation.

.. _whatsnew-2.0-pixelindex:

Fixes and clarification to pixel indexing
=========================================

sunpy uses zero-based indexing when referring to pixels, where the center of the bottom left pixel of a map is at ``[0, 0] * u.pix``.
Several parts of the API have been updated to make sure this is consistently the case across the package.
In particular:

- `sunpy.map.GenericMap.top_right_coord` previously had an off-by-one error in the calculation of the top right coordinate.
  This has been fixed.
- `sunpy.map.GenericMap.center` previously had an off-by-one error in the calculation of the coordinate of the center of a map.
  This has been fixed.
- `sunpy.map.GenericMap.reference_pixel` now returns a zero-based reference pixel.
  This is one pixel less than the previously returned value.
  Note that this means the ``reference_pixel`` now does **not** have the same value as the FITS ``CRPIX`` values, which are one-based indices.
- `sunpy.map.make_fitswcs_header` now correctly interprets the ``reference_pixel`` argument as being zero-based, in previous releases it incorrectly interpreted the ``reference_pixel`` as one-based.

.. _whatsnew-2.0-standards:

Standardization of `~sunpy.map.GenericMap.submap` and ``sunpy.map.GenericMap.draw_rectangle``
==============================================================================================

Both `~sunpy.map.GenericMap.submap` and ``sunpy.map.GenericMap.draw_rectangle`` allow specification of "rectangles" in world (spherical) coordinates.
In versions prior to 2.0 you passed the coordinates of the rectangle to ``draw_rectangle`` as a bottom left coordinate, and a height and width, but for submap you passed it as a bottom left and a top right.
In 2.0 the way you call both methods has changed, to accept a bottom left and then either width and height or a top right coordinate.
As part of this change, the ``top_right``, ``width``, and ``height`` arguments **must** always be keyword arguments, i.e. ``width=10*u.arcsec``

This change allows you to give the same rectangle specification to `~sunpy.map.GenericMap.submap` as to ``sunpy.map.GenericMap.draw_rectangle``.
Which is especially useful when you wish to plot a cropped area of a map, along with it's context in the parent map::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> import matplotlib.pyplot as plt

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE

    >>> aia = sunpy.map.Map(AIA_171_IMAGE)

    >>> bottom_left = SkyCoord(-100 * u.arcsec, -100 * u.arcsec, frame=aia.coordinate_frame)
    >>> width = 500 * u.arcsec
    >>> height = 300 * u.arcsec

    >>> sub_aia = aia.submap(bottom_left, width=width, height=height)

    >>> fig = plt.figure()
    >>> ax1 = fig.add_subplot(1, 2, 1, projection=aia)
    >>> aia.plot(axes=ax1)
    >>> aia.draw_rectangle(bottom_left, width=width, height=height)

    >>> ax2 = fig.add_subplot(1, 2, 2, projection=sub_aia)
    >>> sub_aia.plot(axes=ax2)


Both these methods delegate the input parsing to a new utility function `sunpy.coordinates.utils.get_rectangle_coordinates`.

.. _whatsnew-2.0-overview:

Graphical overview for Map and MapSequence
==========================================

There are new methods to produce graphical overviews for `Map <sunpy.map.map_factory.MapFactory>` and `~sunpy.map.MapSequence` instances: :meth:`~sunpy.map.GenericMap.quicklook` and :meth:`~sunpy.map.MapSequence.quicklook`, respectively.
This graphical overview opens the default web browser and uses `HTML <https://en.wikipedia.org/wiki/HTML>`__ to show a table of metadata, a histogram of the pixel values in the data, and a  `histogram-equalized <https://en.wikipedia.org/wiki/Histogram_equalization>`__ image of the data.
Here's an example of the output for a `~sunpy.map.MapSequence` instance:

.. generate:: html
    :html_border:

    from sunpy.map import Map
    import sunpy.data.sample
    seq = Map(sunpy.data.sample.HMI_LOS_IMAGE,
              sunpy.data.sample.AIA_1600_IMAGE,
              sunpy.data.sample.EIT_195_IMAGE,
              sequence=True)
    print(seq._repr_html_())

If you are using `Jupyter Notebook <https://jupyter.org/>`__, there is no need to call these methods explicitly to see this graphical overview.
If you type just the name of the instance, the graphical overview is shown within the notebook itself as a rich representation of the instance, instead of the typical text representation.

.. _whatsnew-2.0-coordframe:

Differential rotation in the coordinate framework
=================================================

The rotation rate of solar features varies with heliographic latitude, this rotation is called "differential rotation".
SunPy has already included functionality in the `sunpy.physics.differential_rotation` module to transform coordinates and `Maps <sunpy.map.GenericMap>` to account for the rotation of the Sun.
SunPy now provides differential-rotation functionality integrated directly into the `coordinate framework <sunpy.coordinates>` using the `~sunpy.coordinates.metaframes.RotatedSunFrame` class.
Here are examples of using this class:

.. minigallery:: sunpy.coordinates.RotatedSunFrame


A detailed write-up of how to use `~sunpy.coordinates.metaframes.RotatedSunFrame` can be found :ref:`at the RotatedSunFrame documentation<sunpy-coordinates-rotatedsunframe>`.

.. _whatsnew-2.0-carrington:

Changes to Carrington coordinates
=================================

We have refined our approach for heliographic Carrington coordinates to best support high-resolution imagery of the Sun, including from observatories that are at distances from the Sun that is significantly different from 1 AU (e.g., `Solar Orbiter <https://en.wikipedia.org/wiki/Solar_Orbiter>`__).
Our `~sunpy.coordinates.frames.HeliographicCarrington` coordinate frame is now expressly intended for the co-alignment of images of the Sun's surface from different observatories.
`~sunpy.coordinates.frames.HeliographicCarrington` now requires the specification of the observer location (Earth or otherwise) because the light travel time between the Sun and the observer is accounted for.
SunPy output now matches the calculations by `JPL Horizons <https://ssd.jpl.nasa.gov/?horizons>`__ and `SPICE <https://naif.jpl.nasa.gov/naif/>`__.
There may be small differences compared to Carrington coordinates computed by groups that do not use modern parameter values or the same assumptions for the methodology.

Importantly, the Carrington longitude that is now calculated (including using :func:`sunpy.coordinates.sun.L0`) will not match earlier versions of SunPy.
A detailed write-up of the calculation approach and comparisons to other resources can be found :ref:`at carrington's functionality documentation<sunpy-coordinates-carrington>`.

.. _whatsnew-2.0-proxy:

Download behind proxies
=======================

With the release of parfive 1.1, sunpy has been patched to be able to utilize proxy servers when downloading files.

* Proxy URL is read from the environment variables ``HTTP_PROXY`` or ``HTTPS_PROXY``.
* Proxy Authentication ``proxy_auth`` should be passed as a ``aiohttp.BasicAuth`` object, explicitly by the user.
* Proxy Headers ``proxy_headers`` should be passed as `dict` object, explicitly by the user.

For example if you use a bash terminal:

.. code-block:: bash

    $ HTTP_PROXY=http://user:password@proxyserver.com:3128
    $ HTTPS_PROXY=https://user:password@proxyserver.com:3128
    $ export HTTP_PROXY
    $ export HTTPS_PROXY

these will be used to enable downloads through a proxy.

.. _whatsnew-2.0-citation:

Citation update
===============

A paper discussing sunpy 1.0 was accepted in The Astrophysical Journal and you can find the bibtex for it by running::

    >>> import sunpy
    >>> sunpy.__citation__

or `accessing the website directly <https://iopscience.iop.org/article/10.3847/1538-4357/ab4f7a>`__.

Previous update: sunpy 1.1
==========================

In case you never updated to the intermediate release (sunpy 1.1) the whatsnew contains the major changes from that release: :ref:`whatsnew-1.1`
***************
Release History
***************

A new version of sunpy core is released twice a year, with target release months of May and November.
Versions x.0 are long term support (LTS) releases, and are supported for 12 months.
In between versions x.1 are non-LTS releases, and supported for 6 months until the next LTS release.

.. toctree::
   :maxdepth: 1

   changelog
   4.0
   3.1
   3.0
   2.1
   2.0
   1.1
   1.0
   0.9
   0.8
.. _whatsnew-4.0:

************************
What's New in SunPy 4.0?
************************
The SunPy project is pleased to announce the 4.0 release of the sunpy core package.

On this page, you can read about some of the big changes in this release.

.. contents::
    :local:
    :depth: 1

SunPy 4.0 also includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

TODO: fill this in at release time.

Increase in required package versions
=====================================
We have bumped the minimum version of several packages we depend on; these are the new minimum versions:

- python >= 3.8

Better printing of metadata
===========================
Printing a `.MetaDict` now prints each entry on a new line, making it much easier to read::

  >>> from sunpy.data.sample import AIA_171_IMAGE  # doctest: +REMOTE_DATA
  >>> from sunpy.map import Map
  >>> m = Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
  >>> print(m.meta)  # doctest: +REMOTE_DATA
  simple: True
  bitpix: -32
  naxis: 2
  naxis1: 1024
  naxis2: 1024
  ...

Contributors to this Release
============================

The people who have contributed to the code for this release are:

TODO: fill this in at release time.

Where a * indicates that this release contains their first contribution to SunPy.
.. doctest-skip-all

.. _whatsnew-3.0:

************************
What's New in SunPy 3.0?
************************

Overview
========
The SunPy project is pleased to announce the 3.0 release of the sunpy package.

On this page, you can read about some of the big changes in this release:

* :ref:`whatsnew-3.0-map-plotting`
* :ref:`whatsnew-3.0-map-eui`
* :ref:`whatsnew-3.0-map-metadata`
* :ref:`whatsnew-3.0-instr`
* :ref:`whatsnew-3.0-python`

SunPy 3.0 also includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 526 commits have been added since 2.1
* 53 issues have been closed since 2.1
* 109 pull requests have been merged since 2.1
* 26 people have contributed since 2.1
* 10 of which are new contributors

Please find below a selection of what we consider to be the biggest changes or features with this release.

.. _whatsnew-3.0-map-plotting:

Improvements to Visualization of Maps
=====================================

In this release a number of Map visualisation methods have been improved to be more coordinate aware than previous releases.
The first new feature is :meth:`sunpy.map.GenericMap.draw_quadrangle` which replaces the old ``draw_rectangle`` method.
:meth:`~sunpy.map.GenericMap.draw_quadrangle` draws rectangles where the edges follow lines of constant longitude and latitude in the coordinate system in which the box is drawn.

.. minigallery:: sunpy.map.GenericMap.draw_quadrangle

The next change is to :meth:`sunpy.map.GenericMap.draw_limb`, which now supports drawing the limb as seen by any observer on any image.
This means that it is possible to visualize the sections of the limb which are visible in one map on another, or to draw the limb as seen from helioprojective map observer on a synoptic map.
For a demonstration of this new functionality see the first figure in the :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_different_observers.py` example.

The last major change is to :meth:`sunpy.map.GenericMap.plot`, which now has an ``autoalign=True`` keyword argument.
This option when set to `True` will use :meth:`~matplotlib.axes.Axes.pcolormesh` to transform the image being plotted the correct coordinates of the axes it is being plotted on.
This allows for plots of two images with different projections to be visualized together correctly.
It is worth noting that using ``autoalign=True`` is computationally expensive, and when used with interactive plots there is a significant performance penalty every time the plot is modified.
See the :ref:`sphx_glr_generated_gallery_map_transformations_autoalign_aia_hmi.py` example for details.

.. _whatsnew-3.0-map-eui:

Improved Support for Solar Orbiter's EUI Instrument in Map
==========================================================

A new map source to support data from the Extreme Ultraviolet Imager (EUI) instrument on the Solar Orbiter (SolO) spacecraft has been added.
This source improves the accuracy of the observer position by using the heliocentric inertial coordinates as well as correctly setting the processing level, exposure time and colormap.
Data from EUI will automatically load using this source via `sunpy.map.Map`.

.. _whatsnew-3.0-map-metadata:

Inspect history of map metadata changes
==========================================

The ``.meta`` property of a `~sunpy.map.GenericMap` now keeps a record of the contents of the metadata (normally a FITS header) when it was created.
This can be accessed via the `~sunpy.util.metadata.MetaDict.original_meta` property.
This allows any changes made by sunpy or by the user directly to be tracked with the following properties:

* `~sunpy.util.metadata.MetaDict.added_items`
* `~sunpy.util.metadata.MetaDict.removed_items`
* `~sunpy.util.metadata.MetaDict.modified_items`

See the new :ref:`sphx_glr_generated_gallery_map_map_metadata_modification.py` example for details.

.. _whatsnew-3.0-instr:

``sunpy.instr`` moved to ``sunkit-instruments``
===============================================

The ``sunpy.instr`` subpackage has been moved to a separate affiliated package called `sunkit-instruments <https://docs.sunpy.org/projects/sunkit-instruments/>`__.
This has been done to make the core package align with the goal that instrument specific analysis and processing code should live in affiliated packages.

.. _whatsnew-3.0-python:

Increase in required package versions
=====================================
We have bumped the minimum version of several packages we depend on; these are the new minimum versions for sunpy 3.0:

- asdf>=2.6.0
- astropy >= 4.1.0
- beautifulsoup4>=4.8.0
- dask[array]>=2.0.0
- drms>=0.6.1
- glymur>=0.8.18,!=0.9.0
- h5netcdf>=0.8.1
- matplotlib>=3.1.0
- numpy >= 1.16.0
- pandas>=0.24.0
- parfive >= 1.2.0
- python-dateutil>=2.8.0
- scipy >= 1.3.0
- scipy>=1.3.0
- sqlalchemy>=1.3.4
- tqdm>=4.32.1
- zeep>=3.4.0

.. _whatsnew-3.0-contributors:

Contributors to this Release
============================

The people who have contributed to the code for this release are:

-  Abhijeet Manhas
-  Abhishek Pandey
-  Adwait Bhope  *
-  Albert Y. Shih
-  Amarjit Singh Gaba  *
-  Aryan Chouhan
-  David Stansby
-  Jeffrey Aaron Paul
-  Kateryna Ivashkiv
-  Kaustubh Chaudhari  *
-  Kritika Ranjan
-  Laura Hayes
-  Megh Dedhia  *
-  Monica Bobra
-  Mouloudi Mohamed Lyes  *
-  Nabil Freij
-  Nakul Shahdadpuri
-  Ratul Das  *
-  Samriddhi Agarwal  *
-  Shane Maloney
-  Stuart Mumford
-  Tathagata Paul
-  Thomas A Caswell
-  Varun Bankar  *
-  Will Barnes
-  Yukie Nomiya  *

Where a * indicates that this release contains their first contribution to sunpy.
.. _changelog:

**************
Full Changelog
**************

.. changelog::
   :towncrier: ../../
   :towncrier-skip-if-empty:
   :changelog_file: ../../CHANGELOG.rst
.. doctest-skip-all

.. _whatsnew-1.1:

************************
What's New in SunPy 1.1?
************************

Overview
========

The SunPy project is pleased to announce the 1.1 release of the sunpy package.

The headline changes in 1.1 are:

* The `~sunpy.coordinates` subpackage now supports four additional coordinate frames (HCI, HEE, GSE, and GEI).
* A new subpackage `sunpy.data.data_manager` has been added to support versioned data for functions and methods.
* Support in `sunpy.map` and `sunpy.net` for the SUVI instrument on GOES satellites.
* Initial support for WISPR data from Parker Solar Probe in `sunpy.map`.
* The import times for `sunpy` and some subpackages are significantly shorter, with no loss of functionality.

On this page, you can read about some of the big changes in this release:

* :ref:`whatsnew-1.1-python`
* :ref:`whatsnew-1.1-coordinates`
* :ref:`whatsnew-1.1-dl_manager`
* :ref:`whatsnew-1.1-SUVI`
* :ref:`whatsnew-1.1-WISPR`
* :ref:`whatsnew-1.1-importtime`

SunPy 1.1 also includes a large number of smaller improvements and bug fixes, which are described in the :ref:`changelog`.

By the numbers:

* 1137 commits have been added since 1.0
* 106 issues have been closed since 1.0
* 242 pull requests have been merged since 1.0
* 24 people have contributed since 1.0
* 10 new contributors

.. _whatsnew-1.1-python:

Supported versions of Python
============================

Like SunPy 1.0, 1.1 comes with support for Python versions 3.6, 3.7 and support for 3.8 has been added.

.. _whatsnew-1.1-coordinates:

New coordinate frames
=====================

The `~sunpy.coordinates` subpackage now supports four additional coordinate frames of interest to solar physics:

* `Heliocentric Inertial (HCI) <sunpy.coordinates.frames.HeliocentricInertial>`
* `Heliocentric Earth Ecliptic (HEE) <sunpy.coordinates.frames.HeliocentricEarthEcliptic>`
* `Geocentric Solar Ecliptic (GSE) <sunpy.coordinates.frames.GeocentricSolarEcliptic>`
* `Geocentric Earth Equatorial (GEI) <sunpy.coordinates.frames.GeocentricEarthEquatorial>`

The following transformation graph illustrates how all of these coordinate frames can be transformed to any other frame in `sunpy.coordinates` or `astropy.coordinates`:

.. graphviz::

   digraph {
        ICRS [label=ICRS]
        HCRS [label=HCRS]
        HeliocentricMeanEcliptic [label="Heliocentric Aries Ecliptic (HAE)"]
        HeliographicStonyhurst [label="Heliographic Stonyhurst (HGS)\nHeliocentric Earth Equatorial (HEEQ)"]
        HeliocentricEarthEcliptic [label="Heliocentric Earth Ecliptic (HEE)"]
        GeocentricEarthEquatorial [label="Geocentric Earth Equatorial (GEI)"]
        HeliographicCarrington [label="Heliographic Carrington (HGC)"]
        Heliocentric [label="Heliocentric Cartesian (HCC)"]
        HeliocentricInertial [label="Heliocentric Inertial (HCI)"]
        Helioprojective [label="Helioprojective Cartesian (HPC)"]
        GeocentricSolarEcliptic [label="Geocentric Solar Ecliptic (GSE)"]
        ICRS -> HCRS
        ICRS -> HeliocentricMeanEcliptic
        HCRS -> ICRS
        HCRS -> HeliographicStonyhurst
        HeliocentricMeanEcliptic -> ICRS
        HeliocentricMeanEcliptic -> HeliocentricEarthEcliptic
        HeliocentricMeanEcliptic -> GeocentricEarthEquatorial
        HeliographicStonyhurst -> HeliographicCarrington
        HeliographicStonyhurst -> Heliocentric
        HeliographicStonyhurst -> HCRS
        HeliographicStonyhurst -> HeliocentricInertial
        HeliographicCarrington -> HeliographicStonyhurst
        Heliocentric -> Helioprojective
        Heliocentric -> HeliographicStonyhurst
        Helioprojective -> Heliocentric
        HeliocentricEarthEcliptic -> HeliocentricMeanEcliptic
        HeliocentricEarthEcliptic -> GeocentricSolarEcliptic
        GeocentricSolarEcliptic -> HeliocentricEarthEcliptic
        HeliocentricInertial -> HeliographicStonyhurst
        GeocentricEarthEquatorial -> HeliocentricMeanEcliptic
        subgraph cluster_astropy {
                color=blue
                fontcolor=blue
                penwidth=2
                label=<<b>Frames implemented in Astropy</b>>
                ICRS
                HCRS
                HeliocentricMeanEcliptic
                astropy [label="Other Astropy frames" shape=box3d style=filled]
                geocentric [label="Earth-centered frames\n(including GEO)" shape=box3d style=filled]
                astropy -> ICRS
                geocentric -> ICRS
                ICRS -> astropy
                ICRS -> geocentric
        }
        subgraph cluster_sunpy {
                color=crimson
                fontcolor=crimson
                penwidth=2
                label=<<b>Frames implemented in SunPy</b>>
                Helioprojective
                Heliocentric
                HeliographicStonyhurst
                HeliographicCarrington
                subgraph cluster_sunpy11 {
                        color=chocolate
                        fontcolor=chocolate
                        label=<<b>Added in SunPy 1.1</b>>
                        HeliocentricInertial
                        HeliocentricEarthEcliptic
                        GeocentricSolarEcliptic
                        GeocentricEarthEquatorial
                }
        }
        newrank=true
   }

See our :ref:`coordinates documentation <sunpy-coordinates>` for a table of the currently supported coordinate systems and the corresponding frame classes.

.. _whatsnew-1.1-dl_manager:

Manager for Versioned Data Files
================================

SunPy 1.1 provides a data manager for versioning and caching remote files.
The objective of this is to provide a way for data required for functions, such as instrument correction routines, to depend on non-local data in a reliable way.
The data manager also guarantees that a specific version of the code uses a specific data file, with the ability for users to specify updated files.

This works by providing the URL of a remote file and a SHA256 hash to the `sunpy.data.manager.require <sunpy.data.data_manager.DataManager.require>` decorator which can be added to functions that require these specific data files from remote sources.
If the specified hash does not match that of the remote version, an exception is raised to make the user aware of any changes on the remote server or corruption of local files.
Additionally, `sunpy.data.cache <sunpy.data.data_manager.Cache>` can be used to avoid re-downloading files that already exist on a user's local machine, thus saving disk space and internet bandwidth.

See :ref:`remote_data` for more details.

.. _whatsnew-1.1-SUVI:

Support for SUVI Data
=====================

The Solar Ultraviolet Imager (SUVI) is a EUV instrument imaging the full disk of the Sun in six passbands, and is onboard the latest of the Geostationary Operational Environmental Satellite (GOES) missions.
`~sunpy.map.sources.SUVIMap` provides SunPy map support for loading SUVI FITS image data, and the `~sunpy.net.dataretriever.SUVIClient` adds support to search for SUVI data hosted by NOAA via `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`. It supports searching for wavelength, level of data (level 2 data which consists of stacked level 1b images and original level 1b files), as well as GOES satellite number (>= GOES 16).

.. _whatsnew-1.1-WISPR:

Initial Support for WISPR Images
================================

Following the first data release from Parker Solar Probe, SunPy 1.1 supports loading WISPR imaging data into a `~sunpy.map.GenericMap`.
Due to the complex projections in the WISPR data this involved changing the way sunpy converts FITS headers into `astropy.wcs.WCS` objects.
It is expected that sunpy 2.0 will include more complete support for WISPR data.

.. image:: 1.1-wispr.png
   :alt: A plot from SunPy of a WISPR level 3 file.

.. _whatsnew-1.1-importtime:

Speeding up import times
========================

We know that the initial import of `sunpy` or its subpackages can feel like it takes a long time, particularly on slower machines.
Some of that import time can be the result of importing other modules or external packages that are required for specialized functionality that a user may not ever actually use.
We have identified the most egregious cases and deferred those imports of dependencies to when they are actually needed.
For example, the initial import of `sunpy.map` is now ~40% faster, with no loss of functionality.
We will continue to look for ways to improve import times.

.. _whatsnew-1.1-renamed-removed:

Notable Breaking Changes or Removed functionality
=================================================

- Importing `sunpy.timeseries` no longer automatically imports
  Matplotlib. (`#3376 <https://github.com/sunpy/sunpy/pull/3376>`__)
- `sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries.peek` now checks that the `type` argument is a
  valid string, and raises a `ValueError` if it isn't. (`#3378 <https://github.com/sunpy/sunpy/pull/3378>`__)
- Observer-based coordinate frames (`~sunpy.coordinates.frames.Heliocentric` and `~sunpy.coordinates.frames.Helioprojective`) no longer assume a default observer (Earth) if no observer is specified.  These frames can now be used with no observer specified, but most transformations cannot be performed for such frames.  This removal of a default observer only affects `sunpy.coordinates`, and has no impact on the default observer in `sunpy.map`. (`#3388 <https://github.com/sunpy/sunpy/pull/3388>`__)
- The colormap stored in SunPy's Map subclasses (ie. ``map.plot_settings['cmap']``)
  can now be colormap string instead of the full `matplotlib.colors.Colormap`
  object. To get the full `~matplotlib.colors.Colormap` object use the new attribute
  ``map.cmap``. (`#3412 <https://github.com/sunpy/sunpy/pull/3412>`__)
- Fix a warning in `sunpy.map.GenericMap.rotate` where the truth value of an array
  was being calculated. This changes the behavior of
  `~sunpy.map.GenericMap.rotate` when the ``angle=`` parameter is not an
  `~astropy.units.Quantity` object to raise `TypeError` rather than `ValueError`. (`#3456 <https://github.com/sunpy/sunpy/pull/3456>`__)
- Removed the step of repairing images (replacing non-finite entries with local mean) before coaligning them. The user is expected to do this themselves before coaligning images. If NaNs/non-finite entries are present, a warning is thrown.
  The function ``sunpy.image.coalignment.repair_image_nonfinite`` is deprecated. (`#3287 <https://github.com/sunpy/sunpy/pull/3287>`__)
- The method to convert a `~sunpy.coordinates.frames.Helioprojective` frame from 2D to 3D has been renamed from ``sunpy.coordinates.frames.Helioprojective.calculate_distance`` to `~sunpy.coordinates.frames.Helioprojective.make_3d`.  This method is not typically directly called by users. (`#3389 <https://github.com/sunpy/sunpy/pull/3389>`__)
- ``sunpy.visualization.animator.ImageAnimatorWCS`` is now deprecated in favour of
  `~sunpy.visualization.animator.ArrayAnimatorWCS`. (`#3407 <https://github.com/sunpy/sunpy/pull/3407>`__)
- ``sunpy.cm`` has been moved to `sunpy.visualization.colormaps` and will be
  removed in a future version. (`#3410 <https://github.com/sunpy/sunpy/pull/3410>`__)


Full Change Log
===============

To see a detailed list of all changes in version v1.1, including changes in API, please see the :ref:`changelog`.
.. _units-coordinates-sunpy:

Units and Coordinates in sunpy
******************************

This section of the guide will talk about representing physical units and
physical coordinates in sunpy. sunpy makes use of :ref:`Astropy <astropy:astropy-coordinates>` for
both these tasks.


Units in sunpy
==============

All functions in sunpy that accept or return numbers associated with physcial
quantities accept and return `~astropy.units.Quantity` objects. These objects
represent a number (or an array of numbers) and a unit. This means sunpy is
always explicit about the units associated with a value. Quantities and units
are powerful tools for keeping track of variables with a physical meaning and
make it straightforward to convert the same physical quantity into different units.

In this section of the guide we will give a quick introduction to `astropy.units`
and then demostrate how to use units with sunpy.

To use units we must first import them from Astropy. To save on typing we usually
import units as ``u``::

   >>> import astropy.units as u

Once we have imported units we can create a quantity by multiplying a number by
a unit::

   >>> length = 10 * u.meter
   >>> length
   <Quantity 10. m>

A `~astropy.units.Quantity` has both a ``.unit`` and a ``.value`` attribute::

  >>> length.value
  10.0

  >>> length.unit
  Unit("m")

These `~astropy.units.Quantity` objects can also be converted to other units, or
unit systems::

  >>> length.to(u.km)
  <Quantity 0.01 km>

  >>> length.cgs
  <Quantity 1000. cm>

Probably most usefully, `~astropy.units.Quantity` objects will propagate units
through arithmetic operations when appropriate::

  >>> distance_start = 10 * u.mm
  >>> distance_end = 23 * u.km
  >>> length = distance_end - distance_start
  >>> length
  <Quantity 22.99999 km>

  >>> time = 15 * u.minute
  >>> speed = length / time
  >>> speed
  <Quantity 1.53333267 km / min>

However, operations which do not make physical sense for the units specified will cause an error::

  >>> length + time
  Traceback (most recent call last):
  ...
  astropy.units.core.UnitConversionError: Can only apply 'add' function to quantities with compatible dimensions


Quantities as function arguments
================================

An extremely useful addition to the base functionality of Quanitities is the ``@u.quantity_input`` decorator.
This allows you to specify required units for function arguments to ensure that the calculation within that
function always make physical sense. For instance, if we defined a function to calculate speed as above,
we might want the distance and time as inputs::

  >>> def speed(length, time):
  ...     return length / time

However, this requires that length and time both have the appropriate units. We therefore want to use
`~astropy.units.quantity_input` to enforce this, here we use
`function annotations <https://python-3-for-scientists.readthedocs.io/en/latest/python3_features.html#function-annotations>`__
to specify the units.

  >>> @u.quantity_input
  ... def speed(length: u.m, time: u.s):
  ...     return length / time

Now, when this function is called, if the units of length and time are not convertible to the units specified,
an error will be raised stating that the units are incorrect or missing::

  >>> speed(1*u.m, 10*u.m)
  Traceback (most recent call last):
  ...
  astropy.units.core.UnitsError: Argument 'time' to function 'speed' must be in units convertible to 's'.

  >>> speed(1*u.m, 10)
  ...
  Traceback (most recent call last):
  ...
  TypeError: Argument 'time' to function 'speed' has no 'unit' attribute. ... pass in an astropy Quantity instead.

Note that the units of the inputs do not have to be exactly the same as those in the function definition, as long
as they can be converted to those units. So for instance, passing in a time in minutes still works even though we
specified ``time: u.s``::

  >>> speed(1*u.m, 1*u.minute)
  <Quantity 1. m / min>

This may still not be quite as we want it, since we wanted the input time in seconds but the output is in m/min.
We can correct this by defining the function with an additional annotation::

  >>> @u.quantity_input
  ... def speed(length: u.m, time: u.s) -> u.m/u.s:
  ...     return length / time

This will force the output of the function to be converted to m/s before returning, so that you will always
have the same units on the output from this function::

  >>> speed(1*u.m, 1*u.minute)
  <Quantity 0.01666667 m / s>

Physical Coordinates in sunpy
=============================

In much the same way as `~astropy.units` are used for representing physical
quantities, sunpy uses `astropy.coordinates` to represent points in physical
space. This applies to both points in 3D space and projected coordinates in
images.

The astropy coordinates module is primarily used through the
`~astropy.coordinates.SkyCoord` class::

  >>> from astropy.coordinates import SkyCoord

To enable the use of the solar physics specific frames defined in sunpy we also
need to import them::

  >>> from sunpy.coordinates import frames

A SkyCoord object to represent a point on the Sun can then be created::

  >>> c = SkyCoord(70*u.deg, -30*u.deg, obstime="2017-08-01",
  ...              frame=frames.HeliographicStonyhurst)
  >>> c
  <SkyCoord (HeliographicStonyhurst: obstime=2017-08-01T00:00:00.000, rsun=695700.0 km): (lon, lat) in deg
      (70., -30.)>

This `~astropy.coordinates.SkyCoord` object can then be transformed to any
other coordinate frame defined either in Astropy or sunpy, for example::

  >>> c.transform_to(frames.Helioprojective(observer="earth"))
  <SkyCoord (Helioprojective: obstime=2017-08-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
      (769.96270814, -498.89715922, 1.51668773e+08)>


It is also possible to convert three dimensional positions to astrophysical
frames defined in Astropy, for example `~astropy.coordinates.ICRS`.

  >>> c.transform_to('icrs')
  <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, km)
    (49.84856512, 0.05394699, 1417743.94689472)>



Observer Location
-----------------

Both `~sunpy.coordinates.frames.Helioprojective` and
`~sunpy.coordinates.frames.Heliocentric` frames are defined based on the
position of the observer. Therefore to transform either of these frames to a
different frame the location of the observer must be known.
The observer can be specified for a coordinate object using the ``observer``
argument to `~astropy.coordinates.SkyCoord`.  For sunpy to calculate the
location of Earth or another solar-system body, it must know the time for
which the coordinate is valid; this is specified with the ``obstime`` argument.

Using the observer location it is possible to convert a coordinate as seen by
one observer to a coordinate seen by another::

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth",
  ...                 obstime="2017-07-26",
  ...                 frame=frames.Helioprojective)

  >>> hpc1.transform_to(frames.Helioprojective(observer="venus",
  ...                                          obstime="2017-07-26"))
  <SkyCoord (Helioprojective: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'venus'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (-1285.47497992, 106.20918654, 0.72405937)>


Using Coordinates with sunpy Map
--------------------------------

.. plot::
   :include-source:

   sunpy Map uses coordinates to specify locations on the image, and to plot
   overlays on plots of maps. When a Map is created, a coordinate frame is
   constructed from the header information. This can be accessed using
   ``.coordinate_frame``:

   >>> import sunpy.map
   >>> from sunpy.data.sample import AIA_171_IMAGE   # doctest: +REMOTE_DATA
   >>> m = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
   >>> m.coordinate_frame  # doctest: +REMOTE_DATA
   <Helioprojective Frame (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
       (-0.00406308, 0.04787238, 1.51846026e+11)>)>

   This can be used when creating a `~astropy.coordinates.SkyCoord` object to set
   the coordinate system to that image:

   >>> from astropy.coordinates import SkyCoord
   >>> import astropy.units as u
   >>> c = SkyCoord(100 * u.arcsec, 10*u.arcsec, frame=m.coordinate_frame)  # doctest: +REMOTE_DATA
   >>> c  # doctest: +REMOTE_DATA
   <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
       (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
       (100., 10.)>

   The `~astropy.coordinates.SkyCoord` object can be converted to a PixelPair object
   using `~sunpy.map.GenericMap.world_to_pixel`:

   >>> pixel_obj = m.world_to_pixel(c) # doctest: +REMOTE_DATA
   >>> pixel_obj # doctest: +REMOTE_DATA
   PixelPair(x=<Quantity 551.7680511 pix>, y=<Quantity 515.18266871 pix>)

   This `~astropy.coordinates.SkyCoord` object could also be used to plot a point
   on top of the map:

   >>> import matplotlib.pyplot as plt
   >>> ax = plt.subplot(projection=m)  # doctest: +SKIP
   >>> m.plot()  # doctest: +SKIP
   <matplotlib.image.AxesImage object at ...>
   >>> _ = ax.plot_coord(c, 'o')  # doctest: +SKIP

For more information on coordinates see :ref:`sunpy-coordinates` section of the :ref:`reference`.
.. _logger:

**************
Logging system
**************

Overview
========

The sunpy logging system is an adapted version of `~astropy.logger.AstropyLogger`.
Its purpose is to provide users the ability to decide which log and warning messages to show,
to capture them, and to send them to a file.

All messages provided by sunpy use this logging facility which is based
on the Python `logging` module rather than print statements.

Messages can have one of several levels, in increasing order of importance:

* DEBUG: Detailed information, typically of interest only when diagnosing
  problems.

* INFO: A message conveying information about the current task, and
  confirming that things are working as expected

* WARNING: An indication that something unexpected happened, and that user
  action may be required.

* ERROR: indicates a more serious issue where something failed but the task is continuing

* CRITICAL: A serious error, indicating that the program itself may be unable to continue running.

By default, all messages except for DEBUG messages are displayed.

Configuring the logging system
==============================
The default configuration for the logger is determined by the default sunpy
configuration file. To make permanent changes to the logger configuration
see the ``[logger]`` section of the Sunpy configuration
file (:doc:`sunpyrc </guide/customization>`).

If you'd like to control the logger configuration for your current session
first import the logger::

    >>> from sunpy import log

or also by::

    >>> import logging
    >>> log = logging.getLogger('sunpy')

The threshold level for messages can be set with::

    >>> log.setLevel('DEBUG')  #doctest: +SKIP

This will display DEBUG and all messages with that level and above. If you'd like to see the fewest
relevant messages you'd set the logging level to WARNING or above.

For other options such as whether to log to a file or what level of messages the log file should
contain, see the the Sunpy configuration file (:doc:`sunpyrc </guide/customization>`).

Context managers
================
If you'd like to
capture messages as they are generated you can do that with a context manager::

    >>> from sunpy import log
    >>> with log.log_to_list() as log_list:  #doctest: +SKIP
    ...    # your code here  # doctest: +SKIP

Once your code is executed, ``log_list`` will be a Python list containing all of the Sunpy
messages during execution. This does not divert the messages from going to a file or to the screen.
It is also possible to send the messages to a custom file with::

    >>> from sunpy import log
    >>> with log.log_to_file('myfile.log'):  #doctest: +SKIP
    ...     # your code here  #doctest: +SKIP

which will save the messages to a local file called ``myfile.log``.
.. _troubleshooting-faq:

************************
Troubleshooting and Bugs
************************

.. _sunpy-version:

Obtaining sunpy version
=======================

To find out your sunpy version number, import it and print the ``__version__`` attribute::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.__version__   # doctest: +SKIP

.. _locating-sunpy-install:

System Info
===========

To quickly collect information on your system, you can use our convenience function ``system_info`` which you can run through: ::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.util.system_info()   # doctest: +SKIP

The output should look something like: ::

    ==========================================================
     sunpy Installation Information

     Sunday, 18. November 2012 11:06PM UT
    ==========================================================

    ###########
     General
    ###########
    OS: Mac OS X 10.8.2 (i386)
    Python: 2.7.3 (64bit)

    ####################
     Required libraries
    ####################
    sunpy: 0.1
    NumPy: 1.6.2
    SciPy: 0.10.1
    Matplotlib: 1.2.x
    PyFITS: 3.0.8
    pandas: 0.8.1

    #######################
     Recommended libraries
    #######################
    beautifulsoup4: 4.1.1
    PyQt: 4.9.4
    SUDS: 0.4'

This information is especially useful if you are running into a bug and need help.

:file:`sunpy` install location
===================================

You can find what directory sunpy is installed in by importing it and printing the ``__file__`` attribute::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.__file__   # doctest: +SKIP

.. _locating-matplotlib-config-dir:

:file:`.sunpy` directory location
=================================

Each user should have a :file:`.sunpy/` directory which should contain a :ref:`sunpyrc <customizing-with-sunpyrc-files>` file.
To locate your :file:`.sunpy/` directory, use :func:`sunpy.print_config`::

    >>> import sunpy as sun   # doctest: +SKIP
    >>> sun.print_config()   # doctest: +SKIP

We use `appdirs <https://github.com/ActiveState/appdirs>`__ to work out the location depending on your operating system.

If you would like to use a different configuration directory, you can do so by specifying the location in your ``SUNPY_CONFIGDIR`` environment variable.

.. _reporting-problems:

Reporting Bugs
==============

If you are having a problem with sunpy, search the `mailing list`_ or the github `issue tracker`_.
It is possible that someone else has already run into your problem.

If not, please provide the following information in your e-mail to the `mailing list`_ or to the github `issue tracker`_:

  * your operating system; (Linux/UNIX users: post the output of ``uname -a``)

  * sunpy version::

        >>> import sunpy   # doctest: +SKIP
        >>> sunpy.util.system_info()   # doctest: +SKIP

  * how you obtained sunpy.

  * any customizations to your ``sunpyrc`` file (see :ref:`customizing-sunpy`).

  * Please try to provide a **minimal**, standalone Python script that demonstrates the problem.
    This is **the** critical step.
    If you can't post a piece of code that we can run and reproduce your error, the chances of getting help are significantly diminished.
    Very often, the mere act of trying to minimize your code to the smallest bit that produces the error will help you find a bug in **your** code that is causing the problem.

.. _`mailing list`: https://groups.google.com/forum/#!forum/sunpy
.. _`issue tracker`:  https://github.com/sunpy/sunpy/issues
************************
SSWIDL/sunpy Cheat Sheet
************************

`SolarSoft (SSWIDL) <https://sohowww.nascom.nasa.gov/solarsoft/>`_ is a
popular IDL software library for solar data analysis, and in fact, many parts
of sunpy are inspired by data structures and functions in SSWIDL. Though IDL and Python are very different it sometimes helps to consider how to translate
simple tasks between the two languages. The primary packages which provide much of the functionality for scientific data analysis in Python are NumPy and SciPy. In the following we assume that those packages are available to you and that you are imported Numpy is imported as np with the following import statement::

    import numpy as np
    import scipy as sp

In the following examples, a and b could be arrays. For python the arrays must be numpy arrays which can be created simply through::

    np.array(a)

where a is a python list of numbers.

**Relational Operators**

=========  ========
 IDL       Python
=========  ========
a EQ b     a == b
a LT b     a < b
a GT b     a > b
a GE b     a >= b
a LE b     a <= b
a NE b     a != b
=========  ========

**Logical Operators**

=========  ========
 IDL       Python
=========  ========
a and b    a and b
a or b     a or b
=========  ========

**Math Functions**

=========  ========
 IDL       Python
=========  ========
cos(a)     np.cos(a)
alog(a)    np.log(a)
alog10(a)  np.alog(a)
exp(a)     np.exp(a)
=========  ========

**Math Constants**

=========  ========
 IDL       Python
=========  ========
!pi        np.pi
exp(1)     np.e
=========  ========

**Arrays Sequences**

============  ========
 IDL          Python
============  ========
indgen(10)    np.arange(0,10)
findgen(10)   np.arange(0,10,dtype=np.float)
============  ========

**Array Creation**

=============  =========
 IDL           Python
=============  =========
dblarr(3,5)    np.zeros((3,5))
intarr(3,5)    np.zeros((3,5),dtype=np.int)
dblarr(3,5)+1  np.ones((3,5))
intarr(3,5)+9  np.zeros((3,5),dtype=np.int) + 9
boolarr(10)    np.zeros(10,dtype=bool)
identity(3)    np.identity(3)
=============  =========

Many more examples can be found on this `page <http://mathesaurus.sourceforge.net/idl-numpy.html>`_
.. _customizing-sunpy:

*****************
Customizing sunpy
*****************

.. _customizing-with-sunpyrc-files:

The :file:`sunpyrc` file
========================

The sunpy core package uses a :file:`sunpyrc` configuration file to customize
certain properties. You can control a number of key features of sunpy such as
where your data will download to. sunpy looks for the ``sunpyrc`` file
in a platform specific directory, which you can see the path for by running::

  >>> import sunpy
  >>> sunpy.print_config()  # doctest: +SKIP
  FILES USED:
    ...
  <BLANKLINE>
  CONFIGURATION:
    [general]
    time_format = %Y-%m-%d %H:%M:%S
    working_dir = ...
  <BLANKLINE>
    [downloads]
    download_dir = ...
    remote_data_manager_dir = ...
    cache_expiry = 10
    sample_dir = ...
  <BLANKLINE>
    [database]
    url = sqlite:////...
  <BLANKLINE>
    [logger]
    log_level = INFO
    use_color = True
    log_warnings = True
    log_exceptions = False
    log_to_file = False
    log_file_level = INFO
    log_file_format = %(asctime)s, %(origin)s, %(levelname)s, %(message)s
  <BLANKLINE>

To maintain your own customizations place a copy of the default sunpyrc file
into the *first* path printed above.
You can use `sunpy.util.config.copy_default_config` to write the default config into the correct place.
Do not edit the default file directly as every time you install or update sunpy, this file will be overwritten.

See below for the example config file.

.. _customizing-with-dynamic-settings:

Dynamic settings
===================

You can also dynamically change the default settings in a python script or
interactively from the python shell. All of the settings are stored in a
Python ConfigParser instance called ``sunpy.config``, which is global to
the sunpy package. Settings can be modified directly, for example::

    import sunpy
    sunpy.config.set('downloads', 'download_dir', '/home/user/Downloads')


.. sunpyrc-sample:

A sample sunpyrc file
--------------------------------------------------------------------

.. only:: html

    `(download) <../_static/sunpyrc>`__

.. literalinclude:: ../../sunpy/data/sunpyrc
.. _installing:

************
Installation
************

Requirements
============

sunpy requires Python 3.8 or higher.

Installing Scientific Python and sunpy
======================================

sunpy is part of the wider ecosystem of scientific Python packages for solar physics.
Therefore a working sunpy installation is more about installing the scientific Python ecosystem than sunpy itself.

If you do not currently have a working scientific Python distribution this guide will set you up with the Miniconda, which makes it easy to install and manage your scientific Python packages.

To install the Miniconda Python distribution follow the instructions at
`here <https://docs.conda.io/en/latest/miniconda.html>`__.
Although Miniconda makes it simple to switch between Python versions, we recommend that new users install the latest Python 3.x version of Miniconda.

The reason we choose Miniconda over Anaconda, is mainly due to the size as Anaconda comes with a full install of packages you probably do not need and this way you have more direct control over what has been installed into your Python virtual environment.
Furthermore, you bypass the need for the conda resolver to sort out your root environment which should make conda faster to use.

Installing sunpy using Miniconda
--------------------------------

To install sunpy launch a system command prompt or the 'Anaconda Prompt' (under Windows).
First configure conda for to add the `conda-forge channel <https://conda-forge.org/>`__::

    conda config --add channels conda-forge
    conda config --set channel_priority strict

and now to install sunpy within the default conda virtual environment::

    $ conda install sunpy

This will install sunpy and every package it needs to function.

.. note::
    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ or a `conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

Updating sunpy
--------------

You can update to the latest version by running::

    conda update sunpy

Installing sunpy on top of an existing scientific Python environment
--------------------------------------------------------------------

This section assumes you already have everything setup, whether that be conda or a Python virtual environment.
These commands are to be executed within these environments.

Conda
^^^^^

If you want to install sunpy within a pre-existing conda environment, you will want to activate the virtual environment and run::

    $ conda activate <name e.g., base>
    $ conda install sunpy

This assumes you have the conda-forge channel added (as above).

Pip
^^^

This is for installing sunpy within a scientific Python distribution or environment, where ``pip`` has been used to install packages.

To acquire a fully working sunpy installation, simply run::

    pip install "sunpy[all]"

.. note::
    If this does not work, it could be due to a missing C compiler (e.g., ``gcc`` or ``clang``) that is required to build sunpy at install.
    Getting the compiler either from your system package manager, XCode or Anaconda should address this.

If you have a reason to want a more minimal installation, you can install sunpy with no optional dependencies, however this means a lot of submodules will not import::

    pip install "sunpy"

It is possible to select which "extra" dependencies you want to install, if you know you only need certain submodules::

    pip install "sunpy[map,timeseries]"

The available options are: ``[asdf]``, ``[dask]``, ``[database]``, ``[image]``, ``[jpeg2000]``, ``[map]``, ``[net]``, ``[timeseries]``, ``[visualization]``.

If you want to develop sunpy we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.

    Do **not** install sunpy or other third-party packages using ``sudo``.

    This error implies you have an incorrectly configured virtual environment or it is not activated.

    If you really do not want to use any virtual environment, you can always do ``pip install --user sunpy``.

Testing sunpy
=============

sunpy provides a method to run the basic test suite that will check that the install has worked correctly.

To run the basic test suite and ensure that your sunpy install is working correctly, use the :func:`sunpy.self_test`::

    import sunpy
    sunpy.self_test()

You will see something like the following in your terminal::

    Starting sunpy self test...
    Checking for packages needed to run sunpy:
    All required and optional sunpy dependencies are installed.
    Starting the sunpy test suite:
    ...

    The tests will run and will report any fails.  You can report these through the `sunpy issue tracker <https://github.com/sunpy/sunpy/issues>`__ and we will strive to help.

It is possible to run this command in a situation where not all packages are installed. If this is the case, you will see the following when you run the test suite::

    Starting sunpy self test...
    Checking for packages needed to run sunpy:
    The following packages are not installed for the sunpy[database] requirement:
    * sqlalchemy
    ...
    You do not have all the required dependencies installed to run the sunpy test suite.
    If you want to run the sunpy tests install the 'tests' extra with `pip install "sunpy[all,tests]"`

This does not mean sunpy is broken, but you will need to install the extra packages to ensure a "complete" installation of sunpy and run the entire test suite.
It is quite likely that you will run into not having the tests dependencies installed.
A brief tour of sunpy
*********************

This brief tutorial will walk you through some of the functionality offered by
the sunpy core package. Start by reading this tutorial and trying out some of
the examples demonstrated. Once you've completed the tutorial check out the
rest of the :doc:`User Guide </guide/index>` for a more thorough look at the
functionality available.

Sample Data
===========
This tour makes use of a number of sample data files which you will need to
download. This will happen when the sample data is imported for the first time.

Maps
====
Maps are the primary data type in sunpy. They are spatially aware data arrays.
There are maps for a 2D image, a time series of 2D images or temporally aligned
2D images.

**Creating a Map**

sunpy supports many different data products from various sources 'out of the
box'. We shall use SDO's AIA instrument as an example in this tutorial. The
general way to create a Map from one of the supported data products is with the
`~sunpy.map.Map` function from the `sunpy.map` submodule.
`~sunpy.map.Map` takes either a filename, a list of
filenames or a data array and header. We can test
`~sunpy.map.Map` with:


.. plot::
    :include-source:

    import sunpy.data.sample
    import sunpy.map

    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia.peek()

This returns a map named ``aia`` which can be manipulated with standard sunpy map commands.
For more information about maps checkout the :doc:`map guide <data_types/maps>`
and the :ref:`map`.

TimeSeries
==========

sunpy handles time series data, fundamental to the study of any real world
phenomenon, by creating a TimeSeries object. A timeseries consists of two parts;
times and measurements taken at those times. The data can either be in your
current Python session, alternatively within a local or remote file.
In the code block that follows, data is taken from a file containing samples
from a file containing samples from the GOES satellite's X-ray Sensors (XRS).


.. plot::
    :include-source:

    import numpy as np
    import sunpy.data.sample
    import sunpy.timeseries as ts

    my_timeseries = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
    my_timeseries.peek()

We've created this timeseries object by passing TimeSeries a string which
represents the name of a GOES lightcurve file. The
`.peek() <sunpy.timeseries.GenericTimeSeries.peek>` method plots the timeseries
data and displays the plot with some default settings. You can also use
`my_timeseries.plot() <sunpy.timeseries.GenericTimeSeries.plot>` if you want more
control over the style of the output plot.

For more information about TimeSeries, check out the
:doc:`timeseries guide <data_types/timeseries>` and the
and the :ref:`timeseries_code_ref`.

Plotting
========

sunpy uses a matplotlib-like interface to its plotting so more complex plots can
be built by combining sunpy with matplotlib. If you're not familiar with
plotting in matplotlib, you should `learn the basics <https://matplotlib.org/users/tutorials.html>`__
before continuing with this guide.

Let's begin by creating a simple plot of an AIA image. To make things easy,
sunpy includes several example files which are used throughout the docs. These
files have names like ``sunpy.data.sample.AIA_171_IMAGE`` and ``sunpy.data.sample.RHESSI_IMAGE``.

Try typing the below example into your interactive Python shell.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample

    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia.peek()

If everything has been configured properly you should see an AIA image with
the default AIA 17.1 colormap, a colorbar on the right-hand side and a title and some
labels.

There is lot going on here, but we will walk you through the example. Briefly,
the first line is importing sunpy, and the second importing the sample data
files. On the third line we create a sunpy Map object which is a spatially-aware
image. On the last line we then plot the `~sunpy.map.Map` object, using the built in 'quick plot'
function `~sunpy.map.GenericMap.peek`.

sunpy uses a matplotlib-like interface to it's plotting so more complex
plots can be built by combining sunpy with matplotlib.

.. plot::
    :include-source:

    import sunpy.map
    import matplotlib.pyplot as plt
    import sunpy.data.sample

    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure()
    ax = plt.subplot(111, projection=aia)

    aia.plot()
    aia.draw_limb()
    aia.draw_grid()
    plt.colorbar()

    plt.show()

For more information check out :ref:`plotting`.

Solar Physical Constants
========================

sunpy contains a convenient list of solar-related physical constants. Here is
a short bit of code to get you started: ::

    >>> from sunpy.sun import constants as con

    # one astronomical unit (the average distance between the Sun and Earth)
    >>> print(con.au)
      Name   = Astronomical Unit
      Value  = 149597870700.0
      Uncertainty  = 0.0
      Unit  = m
      Reference = IAU 2012 Resolution B2

    # the solar radius
    >>> print(con.radius)
      Name   = Nominal solar radius
      Value  = 695700000.0
      Uncertainty  = 0.0
      Unit  = m
      Reference = IAU 2015 Resolution B 3

Not all constants have a shortcut assigned to them (as above). The rest of the constants
are stored in a dictionary. The following code grabs the dictionary and gets all of the
keys.::

    >>> solar_constants = con.constants
    >>> solar_constants.keys()
    dict_keys(['mass', 'radius', 'luminosity', 'mean distance',
               'perihelion distance', 'aphelion distance', 'age',
               'solar flux unit', 'visual magnitude', 'average angular size',
               'surface area', 'average density', 'surface gravity',
               'moment of inertia', 'volume', 'escape velocity', 'oblateness',
               'metallicity', 'sunspot cycle', 'average intensity',
               'effective temperature', 'mass conversion rate', 'center density',
               'center temperature', 'absolute magnitude', 'mean energy production',
               'ellipticity', 'GM', 'W_0', 'sidereal rotation rate',
               'first Carrington rotation (JD TT)',
               'mean synodic period', 'alpha_0',
               'delta_0'])

You can also use the function `sunpy.sun.constants.print_all()` to print out a table of all of the values
available. These constants are provided as a convenience so that everyone is using the same
(accepted) values. For more information check out :ref:`sun_code_ref`.

Quantities and Units
====================

Many capabilities in sunpy make use of physical quantities that are specified
with units. sunpy uses `~astropy.units` to implement this functionality.
Quantities and units are powerful tools for keeping track of variables with
physical meaning and make it straightforward to convert the same physical
quantity into different units. To learn more about the capabilities of
quantities and units, consult :ref:`units-coordinates-sunpy` or
`the astropy tutorial <http://learn.astropy.org/Quantities.html>`__.

To demonstrate this, let's look at the solar radius constant. This is a physical quantity
that can be expressed in length units ::

    >>> from sunpy.sun import constants as con
    >>> con.radius
    <<class 'astropy.constants.iau2015.IAU2015'> name='Nominal solar radius' value=695700000.0 uncertainty=0.0 unit='m' reference='IAU 2015 Resolution B 3'>

shows the solar radius in units of meters.  The same physical quantity can be expressed in different units instead using the ``.to()`` method::

    >>> con.radius.to('km')
    <Quantity 695700. km>

or equivalently::

    >>> import astropy.units as u
    >>> con.radius.to(u.km)
    <Quantity 695700. km>

If, as is sometimes the case, you need just the raw value or the unit from a quantity, you can access these individually
with the ``value`` and ```unit`` attributes, respectively::

    >>> r = con.radius.to(u.km)
    >>> r.value
    695700.0
    >>> r.unit
    Unit("km")

This is useful, but the real power of units is in using them in calculations.
Suppose you have the radius of a circle and would like to calculate its area.
The following code implements this::

    >>> import numpy as np
    >>> import astropy.units as u

    >>> def circle_area(radius):
    ...     return np.pi * radius ** 2

The first line imports numpy, and the second line imports astropy's units
module. The function then calculates the area based on a given radius. When
it does this, it tracks the units of the input and propagates them through
the calculation. Therefore, if we define the radius in meters, the area will
be in meters squared::

    >>> circle_area(4 * u.m)
    <Quantity 50.26548246 m2>

This also works with different units, for example ::

    >>> from astropy.units import imperial
    >>> circle_area(4 * imperial.foot)
    <Quantity 50.26548246 ft2>

As demonstrated above, we can convert between different systems of measurement.
For example, if you want the area of a circle in square feet, but were given
the radius in meters, then you can convert it before passing it into the function::

    >>> circle_area((4 * u.m).to(imperial.foot))
    <Quantity 541.05315022 ft2>

or you can convert the output::

    >>> circle_area(4 * u.m).to(imperial.foot ** 2)
    <Quantity 541.05315022 ft2>


This is an extremely brief summary of the powerful capbilities of Astropy units.  To find out more, see
the `the astropy tutorial <http://learn.astropy.org/Quantities.html>`__ and
`documentation <https://docs.astropy.org/en/stable/units/index.html>`__


Working with Times
==================

sunpy also contains a number of convenience functions for working with dates
and times. Here is a short example: ::

    >>> import sunpy.time

    # parsing a standard time strings
    >>> sunpy.time.parse_time('2004/02/05 12:00')
    <Time object: scale='utc' format='isot' value=2004-02-05T12:00:00.000>

    # This returns a astropy.time.Time object. All sunpy functions which require
    # time as an input sanitize the input using parse_time.

    # the julian day
    >>> sunpy.time.parse_time((2010,4,30)).jd
    2455316.5

    # TimeRange objects are useful for representing ranges of time
    >>> time_range = sunpy.time.TimeRange('2010/03/04 00:10', '2010/03/04 00:20')
    >>> time_range.center
    <Time object: scale='utc' format='isot' value=2010-03-04T00:15:00.000>

For more information about working with time in sunpy checkout the :doc:`time guide <time>`.


Obtaining Data
==============

sunpy supports searching for and fetching data from a variety of sources,
including the `VSO <https://virtualsolar.org/>`__ and the
`JSOC <http://jsoc.stanford.edu/>`__. The majority of sunpy's clients can be
queried using the `sunpy.net.Fido` interface. An example of searching the VSO using this
is below::

  >>> from sunpy.net import Fido, attrs as a

  >>> results = Fido.search(a.Time("2011-09-20T01:00:00", "2011-09-20T02:00:00"),
  ...                       a.Instrument.eit)   # doctest:  +REMOTE_DATA
  >>> Fido.fetch(results, path="./directory/")  # doctest: +SKIP
  ['./directory/efz20110920.010015',
   './directory/efz20110920.010613',
   './directory/efz20110920.011353',
   './directory/efz20110920.011947']

For more information and examples of downloading data with sunpy see :ref:`acquiring_data`.

Database Package
================

The database package can be used to keep a local record of all files downloaded
from the VSO, this means that two searches of the VSO which overlap will not
re-download data.

A simple example of this is shown below::


    >>> import astropy.units as u
    >>> from sunpy.net import Fido, attrs as a
    >>> from sunpy.database import Database

    >>> db = Database()
    >>> db.fetch(a.Time("2011-09-20T01:00:00", "2011-09-20T02:00:00"),
    ...          a.Instrument.aia, a.Sample(45*u.min))  # doctest: +REMOTE_DATA
    >>> db.commit()  # doctest: +REMOTE_DATA
    >>> db  # doctest: +SKIP
    <Table length=4>
     id  observation_time_start observation_time_end ...    download_time      size
    str1         str19                 str19         ...        str19          str7
    ---- ---------------------- -------------------- ... ------------------- -------
       1    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       2    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       3    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0
       4    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0

If you then do a second query::

    >>> db.fetch(a.Time("2011-09-20T01:00:00", "2011-09-20T02:45:00"),
    ...          a.Instrument.aia, a.Sample(45*u.min))  # doctest: +REMOTE_DATA
    >>> db.commit()  # doctest: +REMOTE_DATA
    >>> db  # doctest: +SKIP
    <Table length=6>
     id  observation_time_start observation_time_end ...    download_time      size
    str1         str19                 str19         ...        str19          str7
    ---- ---------------------- -------------------- ... ------------------- -------
       1    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       2    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       3    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0
       4    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0
       5    2011-09-20 02:30:00  2011-09-20 02:30:01 ... 2020-11-21 14:17:51 66200.0
       6    2011-09-20 02:30:00  2011-09-20 02:30:01 ... 2020-11-21 14:17:51 66200.

A query can then be performed against the database to get the records::

    >>> entries = db.search(a.Time("2011-09-20T01:45:00", "2011-09-20T02:15:00"), a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> len(entries)  # doctest: +SKIP
    4

You can see that only two extra records were added to the database.
For more information check out the :ref:`database_guide`.
.. _time-in-sunpy:

*************
Time in sunpy
*************

Working with times and time ranges is a standard task in solar data analysis and as such
sunpy strives to provide convenient and easy methods to do the simple stuff. Python
already provides an object for a time or date through `datetime.datetime`.
However, `datetime.datetime` does not provide support for common time formats used in
solar physics nor leap seconds. To alleviate this, we use `astropy.time.Time` internally
which allows us to provide a superior user experience.

.. _parse-time:

1. Parsing Times
================

Solar data is associated with a number of different time formats. sunpy provides a simple
parsing function which can deal with most every format that a user may encounter. Called
`sunpy.time.parse_time()`, this function can take a variety of inputs.

Strings
-------

The most commonly used are strings and we support a selection of formats
which are matched using regrex. We currently support the following style of string formats::

    "2007-05-04T21:08:12.999999"
    "2007/05/04T21:08:12.999999"
    "2007-05-04T21:08:12.999Z"
    "2007-05-04T21:08:12"
    "2007/05/04T21:08:12"
    "20070504T210812.999999"
    "20070504T210812"
    "2007/05/04 21:08:12"
    "2007/05/04 21:08"
    "2007/05/04 21:08:12.999999"
    "2007-05-04 21:08:12.999999"
    "2007-05-04 21:08:12"
    "2007-05-04 21:08"
    "2007-May-04 21:08:12"
    "2007-May-04 21:08"
    "2007-May-04"
    "2007-05-04"
    "2007/05/04"
    "04-May-2007"
    "04-May-2007 21:08:12.999999"
    "20070504_210812"
    "2012:124:21:08:12"
    "2012:124:21:08:12.999999"
    "20140101000001"
    "2016.05.04_21:08:12_TAI"

If we pass some of these strings into `sunpy.time.parse_time()`::

    >>> from sunpy.time import parse_time
    >>> parse_time('2007-05-04T21:08:12')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('2007/05/04T21:08:12')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('20070504T210812')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('2007-May-04 21:08:12')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('20070504_210812')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>

Each of the above returns the same `~astropy.time.Time` object ``<Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>``.

We also support ``utime``, which is the amount of seconds from 1979-01-01 00:00:00 UTC.
Same as Unix time but this starts 9 years later. The parse_time function also accepts this as input, e.g.::

    >>> parse_time(894316092.00000000, format='utime')
    <Time object: scale='utc' format='utime' value=894316092.0>

Other formats
-------------

`sunpy.time.parse_time()` understands more than just a string.
For example::

    >>> import datetime
    >>> parse_time(datetime.datetime(2007, 5, 4, 21, 8, 12))  # datetime.datetime
    <Time object: scale='utc' format='datetime' value=2007-05-04 21:08:12>
    >>> parse_time(datetime.date(2007, 5, 4))  # datetime.date
    <Time object: scale='utc' format='iso' value=2007-05-04 00:00:00.000>

    >>> parse_time((2007, 5, 4, 21, 8, 12))  # tuple of numbers
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>

    >>> import pandas
    >>> parse_time(pandas.Timestamp('2007-05-04T21:08:12'))  # pandas.Timestamp
    <Time object: scale='utc' format='datetime64' value=2007-05-04T21:08:12.000000000>

    >>> time_ranges = [datetime.datetime(2007, 5, i) for i in range(1, 3)]
    >>> parse_time(pandas.Series(time_ranges))  # pandas.Series
    <Time object: scale='utc' format='datetime' value=[datetime.datetime(2007, 5, 1, 0, 0) datetime.datetime(2007, 5, 2, 0, 0)]>

    >>> parse_time(pandas.DatetimeIndex(time_ranges))  # pandas.DatetimeIndex
    <Time object: scale='utc' format='datetime' value=[datetime.datetime(2007, 5, 1, 0, 0) datetime.datetime(2007, 5, 2, 0, 0)]>

    >>> import numpy as np
    >>> parse_time(np.datetime64('2007-05-04T00'))  # np.datetime64
    <Time object: scale='utc' format='isot' value=2007-05-04T00:00:00.000>
    >>> parse_time(np.arange('2007-05-03', '2007-05-04', dtype='datetime64[D]'))  # np.ndarray
    <Time object: scale='utc' format='isot' value=['2007-05-03T00:00:00.000']>

`astropy.time.Time` API comparison
----------------------------------

`sunpy.time.parse_time` is a wrapper around `astropy.time.Time`. The API is
nearly identical as `~astropy.time.Time` but supports more time input formats.
You can specify the format, scale, precision, location and other arguments just
as you would do with `~astropy.time.Time`. An example::

    >>> times = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']
    >>> parse_time(times, format='isot', scale='tai')
    <Time object: scale='tai' format='isot' value=['1999-01-01T00:00:00.123' '2010-01-01T00:00:00.000']>

Please be aware that all sunpy functions which require time as an input sanitize the input using `~sunpy.time.parse_time`.

2. Time Ranges
==============

A very standard task in data analysis is to have to deal with pairs of times or time
ranges. This occurs very often with plotting or when searching for data. To deal with
time ranges sunpy provides the `sunpy.time.TimeRange` object. A TimeRange object can be created
very easily by providing it with two time strings, a start time and an end time: ::

    >>> from sunpy.time import TimeRange
    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')

You can also pass the start and end times as a tuple: ::

    >>> time_range = TimeRange(('2010/03/04 00:10', '2010/03/04 00:20'))

This object makes use of parse_time() so it can accept a wide variety of time formats.
A time range object can also be created by providing a start time and a duration.
The duration must be provided as a `~astropy.time.TimeDelta` or
time-equivalent `astropy.units.Quantity` or `datetime.timedelta` object
example: ::

    >>> import astropy.units as u
    >>> time_range = TimeRange('2010/03/04 00:10', 400 * u.second)

or: ::

    >>> import astropy.units as u
    >>> from astropy.time import TimeDelta
    >>> time_range = TimeRange('2010/03/04 00:10', TimeDelta(400 * u.second))

or: ::

    >>> from datetime import timedelta
    >>> time_range = TimeRange('2010/03/04 00:10', timedelta(0, 400))

The time range objects provides a number of useful functions. For example, you can easily
get the time at the center of your interval or the length of your interval in minutes
or days or seconds: ::

    >>> time_range.center
    <Time object: scale='utc' format='isot' value=2010-03-04T00:13:20.000>
    >>> time_range.minutes
    <Quantity 6.66666667 min>
    >>> time_range.days
    <Quantity 0.00462963 d>
    >>> time_range.seconds
    <Quantity 400. s>

It also makes it easy to create new time ranges. The functions next() and previous()
do an inplace update to the object by either adding or subtracting the same time interval
. This could be useful if you need to step through a number of time ranges. For example,
if you needed time ranges that spanned 30 minutes over a period of 4 hours you could do: ::

    >>> for a in range(8):
    ...     print(time_range.next())  # doctest: +IGNORE_OUTPUT
        Start: 2010-03-04 00:16:40
        End:   2010-03-04 00:23:20
        Center:2010-03-04 00:20:00
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:23:20
        End:   2010-03-04 00:30:00
        Center:2010-03-04 00:26:40
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:30:00
        End:   2010-03-04 00:36:40
        Center:2010-03-04 00:33:20
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:36:40
        End:   2010-03-04 00:43:20
        Center:2010-03-04 00:40:00
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:43:20
        End:   2010-03-04 00:50:00
        Center:2010-03-04 00:46:40
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:50:00
        End:   2010-03-04 00:56:40
        Center:2010-03-04 00:53:20
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:56:40
        End:   2010-03-04 01:03:20
        Center:2010-03-04 01:00:00
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 01:03:20
        End:   2010-03-04 01:10:00
        Center:2010-03-04 01:06:40
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>

A time range can also be easily split into sub-intervals of equal length, for example to
split a TimeRange object into two new TimeRange objects: ::

    time_range.split(2)

Check out the code reference for the `sunpy.time.TimeRange` object for more information.
.. _guide:

************
User's Guide
************

Welcome to the user guide for the sunpy core package. sunpy is a community-
developed, free and open-source solar data analysis environment. It is meant to
provide the core functionality and tools to analyze solar data with Python.
This guide provides a walkthrough of the major features in sunpy.
For more details checkout the :ref:`reference`.

.. toctree::
   :maxdepth: 2

   installation
   Brief Tour <tour>
   Data Acquisition <acquiring_data/index>
   Data Types <data_types/index>
   Plotting <plotting>
   Units and Coordinates <units-coordinates>
   Time <time>
   customization
   logger
   ssw
   troubleshooting
.. _plotting:

*****************
Plotting in sunpy
*****************

sunpy makes use of `matplotlib <https://matplotlib.org/>`_ for all of its
plotting - as such, it tries to follow the matplotlib plotting philosophy. It is
therefore useful to go over how matplotlib works as background.

1. Matplotlib Tutorial
**********************

The tutorial provided here is a summary of one that can be found in the `matplotlib
usage documentation <https://matplotlib.org/faq/usage_faq.html>`_.

Matplotlib provides two main pathways for plotting. One is meant for interactive use
(e.g. command-line) and the other for non-interactive use (e.g. modules). It is important
to recognize though that the interactive-use pathway (referred to as pyplot) just
provides shortcuts for doing many of the more advanced non-interactive functions in the
background. It is therefore possible to switch between the two as necessary and
it is possible to use pyplot in a non-interactive way. In this manner pyplot
is just a shortcut to making it quicker to set up plot axes and figures.
In order to get access to the full interactive capabilities of pyplot it is
necessary to turn this feature on.
Pylab is another matplotlib usage scenario but it is essentially just pyplot with the
interactive capabilities turned on and numpy and matplotlib imported into the main
namespace.

2. Pyplot
*********
Here is a simple example of pyplot usage.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    plt.plot(range(10), range(10))
    plt.title("A simple Plot")
    plt.show()

The `~matplotlib.pyplot.show` command opens a plot on the screen and blocks
execution until the plot window is closed. The `~matplotlib.pyplot.show`
command only works once. If you were to call `~matplotlib.pyplot.show` again
after the above code is executed nothing happens. This confusing behavior
is something that the matplotlib devs get complaints about often and so this may change.
A discussion about this can be found `here
<https://stackoverflow.com/questions/5524858/matplotlib-show-doesnt-work-twice>`_.
Don't be confused by another command called `~matplotlib.pyplot.draw`.
This is only used while in interactive mode.

To turn on interactivity for pyplot use the command ::

    >>> plt.ion()   # doctest: +SKIP

In interactive mode, the plot will appear at the first `~matplotlib.pyplot.plot`
command and most commands will update the plot as you call them. Here is some
example code::

    >>> plt.plot(range(10), range(10))   # doctest: +SKIP
    >>> plt.title("Simple Plot")   # doctest: +SKIP

In this example, you'll see that the title appears right on the plot when you call it.
Note that in this case the `~matplotlib.pyplot.show` command is useless as the
plot shows up right when you create it. Also note that some commands will not
automatically update the plot and you have to use the `~matplotlib.pyplot.draw`
command. The following command ::

    >>> plt.ioff()   # doctest: +SKIP

turns off interactivity.

3. Advanced Pyplot
******************
If you need more fine-grained control over plots the recommended path is to use pyplot
and access the figures and axes objects. This is shown in the following example.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np
    x = np.arange(0, 10, 0.2)
    y = np.sin(x)
    fig = plt.figure()
    ax = plt.subplot()
    ax.plot(x, y)
    ax.set_xlabel('x')
    plt.show()

In matplotlib, `~matplotlib.figure.Figure` is the top-level container for all plot elements and
`~matplotlib.axes.Axes` is the top-level container for a particular plot. So the above example,
creates a figure then creates an axes and populates the plot in ``ax``. With this method you
now have your hands on the `~matplotlib.axes.Axes` object so you can do things
like change the labels on the x and y axes or add a legend.
In the previous section, pyplot took care of creating these
objects for you so you didn't have to worry about creating them yourself.

4. sunpy Plotting Convention
****************************

To be consistent with matplotlib, sunpy has developed a standard plotting policy
which supports both simple and advanced matplotlib usage. The following examples
focus on the map object but they should be applicable across all of the data
objects.

4.1 peek()
**********

For quick and easy access to a plot all sunpy base objects (i.e. maps, spectra,
timeseries) define their own `~sunpy.map.mapbase.GenericMap.peek` command which
will create a plot for you and show it without you having to deal with any
matplotlib setup. This is so that it is easy to take a quick look at your data.
For example you can make the following plot.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    smap.peek(draw_limb=True)

This creates a plot window with all axes defined, a plot title, and the image of
the map data defined by the contents of the map. In non-interactive mode the
plot window blocks the command line terminal and must be closed before doing anything else.

4.2 plot()
**********

For more advanced plotting the base sunpy objects also provide a
`~sunpy.map.mapbase.GenericMap.plot` command. This command is similar to the
pyplot `~matplotlib.pyplot.imshow` command in that it will create a figure and
axes object for you if you haven't already.

When you create a plot with `~sunpy.map.GenericMap.peek` or
`~sunpy.map.GenericMap.plot`, sunpy will use `astropy.visualization.wcsaxes` to
represent coordinates on the image accurately, for more information see
:ref:`wcsaxes-plotting`.

Using `~sunpy.map.GenericMap.plot` it is possible to customise the look of the
plot by combining sunpy and matplotlib commands, for example you can over plot
contours on the Map:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import astropy.units as u

    import sunpy.map
    import sunpy.data.sample

    aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia_map.plot()
    aia_map.draw_limb()

    # let's add contours as well
    aia_map.draw_contours([10,20,30,40,50,60,70,80,90] * u.percent)

    plt.colorbar()
    plt.show()


In this example, the `~matplotlib.figure.Figure` and
`~astropy.visualization.wcsaxes.WCSAxes` instances are created explicitly, and
then used to modify the plot:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    import sunpy.map
    import sunpy.data.sample

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure()
    # Provide the Map as a projection, which creates a WCSAxes object
    ax = plt.subplot(projection=smap)

    im = smap.plot()

    # Prevent the image from being re-scaled while overplotting.
    ax.set_autoscale_on(False)

    xc = [0,100,1000] * u.arcsec
    yc = [0,100,1000] * u.arcsec

    coords = SkyCoord(xc, yc, frame=smap.coordinate_frame)

    p = ax.plot_coord(coords, 'o')

    # Set title.
    ax.set_title('Custom plot with WCSAxes')

    plt.colorbar()
    plt.show()

It is possible to create the same plot, explicitly not using `~astropy.visualization.wcsaxes`, however, this will not have the features of `~astropy.visualization.wcsaxes` which include correct representation of rotation and plotting in different coordinate systems.
Please see this example :ref:`sphx_glr_generated_gallery_map_plot_frameless_image.py`.

.. _wcsaxes-plotting:

Plotting Maps with wcsaxes
**************************

By default :ref:`map` uses the `astropy.visualization.wcsaxes` module to improve
the representation of world coordinates, and calling
`~sunpy.map.GenericMap.plot` or `~sunpy.map.GenericMap.peek()` will use wcsaxes
for plotting. Unless a standard `matplotlib.axes.Axes` object is explicitly
created.

To explicitly create a `~astropy.visualization.wcsaxes.WCSAxes` instance do the
following ::

    >>> fig = plt.figure()   # doctest: +SKIP
    >>> ax = plt.subplot(projection=smap)   # doctest: +SKIP

when plotting on an `~astropy.visualization.wcsaxes.WCSAxes` axes, it will by
default plot in pixel coordinates, you can override this behavior and plot in
'world' coordinates by getting the transformation from the axes with
``ax.get_transform('world')``. Note: World coordinates are always in **degrees**
so you will have to convert to degrees.::

    >>> smap.plot()   # doctest: +SKIP
    >>> ax.plot((100*u.arcsec).to(u.deg), (500*u.arcsec).to(u.deg),
    ...         transform=ax.get_transform('world'))   # doctest: +SKIP

Finally, here is a more complex example using sunpy maps, wcsaxes and Astropy
units to plot a AIA image and a zoomed in view of an active region.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    from matplotlib import patches
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    import sunpy.map
    import sunpy.data.sample

    # Define a region of interest
    length = 250 * u.arcsec
    x0 = -100 * u.arcsec
    y0 = -400 * u.arcsec

    # Create a sunpy Map, and a second submap over the region of interest.
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    bottom_left = SkyCoord(x0 - length, y0 - length,
                        frame=smap.coordinate_frame)
    top_right = SkyCoord(x0 + length, y0 + length,
                        frame=smap.coordinate_frame)
    submap = smap.submap(bottom_left, top_right=top_right)

    # Create a new matplotlib figure, larger than default.
    fig = plt.figure(figsize=(5, 12))

    # Add a first Axis, using the WCS from the map.
    ax1 = fig.add_subplot(2, 1, 1, projection=smap)

    # Plot the Map on the axes with default settings.
    smap.plot()

    # Draw a box on the image
    smap.draw_quadrangle(bottom_left, height=length * 2, width=length * 2)

    # Create a second axis on the plot.
    ax2 = fig.add_subplot(2, 1, 2, projection=submap)

    submap.plot()

    # Add a overlay grid.
    submap.draw_grid(grid_spacing=10*u.deg)

    # Change the title.
    ax2.set_title('Zoomed View')

    plt.show()
***********
TimeSeries
***********

Time series data are a fundamental part of many data analysis projects
in heliophysics as well as other areas. sunpy therefore provides a TimeSeries
object to handle this type of data. This new object supersedes the now
deprecated Lightcurve datatype. Once you've read through this guide check out
the :doc:`/code_ref/timeseries` for a more thorough look at sunpy TimeSeries
and to see what data sources it currently supports.

Sections 1 and 2 below describe how to create TimeSeries objects.  Section 3
describes how to inspect the TimeSeries object to examine the data itself and
the metadata.  Section 4 describes how to plot a TimeSeries object, Section 5
describes how a TimeSeries object can be manipulated, and Section 6 describes
in detail the TimeSeries metadata.

.. warning::

   The TimeSeries supersedes the previous LightCurve object but does not
   implement data download methods. To download TimeSeries data files use
   `sunpy.net`.

1. Creating a TimeSeries from an observational data source
==========================================================

A TimeSeries object can be created from local files.  For convenience, sunpy can
download several example timeseries of observational data. These files have names like
``sunpy.data.sample.EVE_TIMESERIES`` and ``sunpy.data.sample.GOES_XRS_TIMESERIES``.
To create the sample `sunpy.timeseries.sources.goes.XRSTimeSeries` type the
following into your interactive Python shell: ::

    >>> import sunpy.timeseries as ts
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> my_timeseries = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')  # doctest: +REMOTE_DATA

.. doctest-skip-all

This is calling the `~sunpy.timeseries.TimeSeries` factory to create a time
series from a GOES XRS FITS file. The TimeSeries factory uses `sunpy.io.fits` to
read the FITS file. Note that if you have not downloaded the data already you
should get an error and some instruction on how to download the sample data.

The variable ``my_timeseries`` is a :ref:`timeseries` object. To create one from
a local GOES/XRS FITS file try the following: ::

    >>> my_timeseries = ts.TimeSeries('/mydirectory/myts.fits', source='XRS')   # doctest: +SKIP

sunpy will attempt to detect automatically the instrument source for most FITS
files. However timeseries data are stored in a variety of file types (FITS, txt,
csv) and formats, and so it is not always possible to detect the source. For
this reason, it is good practice to explicitly state the source for the file.
sunpy ships with a number of known instrumental sources.  If you would like
sunpy to include another instrumental source please follow the
`sunpy contribution guide <https://sunpy.org/contribute>`__.
The `~sunpy.timeseries.TimeSeries` factory has the ability to make a list of
TimeSeries objects using a list of filepaths, a folder or a glob, for example: ::

    >>> my_ts_list = ts.TimeSeries('filepath1', 'filepath2', source='XRS')   # doctest: +SKIP
    >>> my_ts_list = ts.TimeSeries('/goesdirectory/', source='XRS')   # doctest: +SKIP
    >>> my_ts_list = ts.TimeSeries(glob, source='XRS')   # doctest: +SKIP

Note that this functionality will only work with files from the same single
source, generating a source specific child of the `~sunpy.timeseries.GenericTimeSeries`
class such as the `~sunpy.timeseries.sources.goes.XRSTimeSeries` above. For this
reason, all the files should be from that same source for the `~sunpy.timeseries.TimeSeries`
factory to work as intended.

1.1 Creating a Single TimeSeries from Multiple Files
====================================================

You can create a single time series from multiple files for a given source using
the keyword argument ``concatenate=True``, such as:

    >>> my_timeseries = ts.TimeSeries('/mydirectory/myts1.fits', '/mydirectory/myts2.fits', source='XRS', concatenate=True)  # doctest: +SKIP

Note these must all be from the same source if using
`~sunpy.timeseries.GenericTimeSeries.concatenate` from within the TimeSeries
factory. However, if time series `~sunpy.timeseries.GenericTimeSeries.concatenate` method
can be used to make a single time series from multiple TimeSeries from different
sources once they are already in the form of TimeSeries objects.

2. Creating Custom TimeSeries
=============================

Sometimes you will have data that you want to create into a TimeSeries. You can
use the factory to create a `~sunpy.timeseries.GenericTimeSeries`
from a variety of data sources currently including `pandas.DataFrame` and
`astropy.table.Table`.

2.1 Creating a TimeSeries from a Pandas DataFrame
=================================================

A TimeSeries object must be supplied with some data when it is
created.  The data can either be in your current Python session, in a
local file, or in a remote file.  Let's create some data and pass
it into a TimeSeries object: ::

    >>> import numpy as np
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

The first line imports the numpy module used to create and store the data.
The second line creates a basic numpy array of values representing a sine wave.
We can use this array along with a suitable time storing object (such as Astropy
`~astropy.time` or a list of `datetime` objects) to make a Pandas
`~pandas.DataFrame`.  A suitable list of times must contain the same
number of values as the data, this can be created using: ::

    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(24*60, 0, -1)]

The Pandas `~pandas.DataFrame` will use the dates list as the index: ::

    >>> from pandas import DataFrame
    >>> data = DataFrame(intensity, index=times, columns=['intensity'])

This `~pandas.DataFrame` can then be used to construct a TimeSeries: ::

    >>> import sunpy.timeseries as ts
    >>> ts_custom = ts.TimeSeries(data)

Furthermore we could specify the metadata/header and units of this time series
by sending them as arguments to the factory: ::

    >>> from collections import OrderedDict
    >>> import astropy.units as u

    >>> meta = OrderedDict({'key':'value'})
    >>> units = OrderedDict([('intensity', u.W/u.m**2)])
    >>> ts_custom = ts.TimeSeries(data, meta, units)

2.2 Creating Custom TimeSeries from an Astropy Table
====================================================

A Pandas `~pandas.DataFrame` is the underlying object used to store
the data within a TimeSeries, so the above example is the most lightweight to
create a custom TimeSeries, but being scientific data it will often be more
convenient to use an Astropy `~astropy.table.Table` and let the factory
convert this.  An advantage of this method is it allows you to include metadata
and Astropy `~astropy.units.quantity.Quantity` values, which are both supported
in tables, without additional arguments.  For example: ::

    >>> import datetime
    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from astropy.table import Table

    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(24*60, 0, -1)]
    >>> intensity = u.Quantity(np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60)))), u.W/u.m**2)
    >>> tbl_meta = {'t_key':'t_value'}
    >>> table = Table([times, intensity], names=['time', 'intensity'], meta=tbl_meta)
    >>> table.add_index('time')
    >>> ts_table = ts.TimeSeries(table)

Note that due to the properties of the `~astropy.time.Time` object, this will be
a mixin column which since it is a single object, limits the versatility of
the `~astropy.table.Table` a little. For more on mixin columns see the `Astropy
docs <https://docs.astropy.org/en/stable/table/mixin_columns.html>`_.  The units
will be taken from the table quantities for each column, the metadata will
simply be the table.meta dictionary.  You can also explicitly add metadata and
units, these will be added to the relevant dictionaries using the dictionary
update method, with the explicit user-given values taking precedence.

    >>> from sunpy.util.metadata import MetaDict
    >>> from collections import OrderedDict
    >>> import astropy.units as u

    >>> meta = MetaDict({'key':'value'})
    >>> units = OrderedDict([('intensity', u.W/u.m**2)])
    >>> ts_table = ts.TimeSeries(table, meta, units)


3. Inspecting TimeSeries & Getting at the Data
===============================================

A time series holds both data as well as meta data and units data. The meta data
for the time series is accessed by: ::

    >>> header = my_timeseries.meta

This references the `~sunpy.timeseries.TimeSeriesMetaData` object with
the header information as read from the source files. A word of caution: many
data sources provide little to no meta data so this variable might be empty.
The meta data is described in more detail later in this guide. Similarly there
are properties for getting `~sunpy.timeseries.GenericTimeSeries.columns`
as a list of strings, `~sunpy.timeseries.GenericTimeSeries.index`
values and `~sunpy.timeseries.GenericTimeSeries.time_range` of
the data.  The actual data in a sunpy TimeSeries object is accessible through
the `~sunpy.timeseries.GenericTimeSeries.data` attribute.  The
data is implemented as a Pandas `~pandas.DataFrame`, so to get a look at what
data you have available use: ::

    >>> my_timeseries.data  # doctest: +SKIP

You can also get a quick overview of that data using: ::

    >>> my_timeseries.data.info()
    <class 'pandas.core.frame.DataFrame'>
    DatetimeIndex: 42177 entries, 2011-06-06 23:59:59.961999 to 2011-06-07 23:59:57.631999
    Data columns (total 2 columns):
    xrsa    42177 non-null float32
    xrsb    42177 non-null float32
    dtypes: float32(2)
    memory usage: 659.0 KB

Time series are columnar data so to get at a particular datum you need to
first index the column, then the element you want. To get the names of the
available columns: ::

    >>> my_timeseries.data.columns
    Index(['xrsa', 'xrsb'], dtype='object')

You can access the 0th element in the column ``xrsa`` with: ::

    >>> my_timeseries.data['xrsa'][0]
    1e-09

You can also grab all of the data at a particular time: ::

    >>> my_timeseries.data['xrsa']['2011-06-07 00:00:02.008999']
    1e-09

This will return a list of entries with times that match the accuracy of the time
you provide. You can consider the data as x or y values: ::

    >>> x = my_timeseries.data.index
    >>> y = my_timeseries.data.values

You can read more about indexing at the `pandas documentation website
<https://pandas.pydata.org/pandas-docs/stable/>`_.

A TimeSeries can also return an Astropy `~astropy.units.quantity.Quantity` for a
given column using the `~sunpy.timeseries.GenericTimeSeries.quantity`
method, this uses the values stored in the data and units stored in the units
dictionary to determine the `~astropy.units.quantity.Quantity`: ::

    >>> quantity = my_timeseries.quantity('xrsa')

4. Plotting
===========

The sunpy TimeSeries object has its own built-in plot methods so that
it is easy to quickly view your time series. To create a plot just
type:

.. plot::
    :include-source:

    import sunpy.timeseries as ts
    import sunpy.data.sample

    ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
    ts.peek()

This will open a Matplotlib plot on your screen. If you want to save this to a PNG
file you can do so from the Matplotlib GUI.

In addition, to enable users to modify the plot it is possible to use the
`~sunpy.timeseries.GenericTimeSeries.plot` command.  This makes it possible to
use the sunpy plot as the foundation for a more complicated figure:

.. plot::
   :include-source:

   import matplotlib.pyplot as plt

   import sunpy.timeseries as ts
   import sunpy.data.sample

   ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
   fig, ax = plt.subplots()
   ts.plot(axes=ax)
   # Modify the figure here
   fig.savefig('figure.png')


5 Manipulating TimeSeries
=========================

5.1 Modifying the Data
======================

Since the timeseries data is stored as a Pandas `~pandas.DataFrame`
you can easily modify the data directly using all of the usual Pandas methods:
for example, you can modify a single cells value using: ::

    >>> my_timeseries.data['xrsa'][0] = 0.1

Or similarly using a datetime values (as string or datetime object): ::

    >>> my_timeseries.data['xrsa']['2012-06-01 23:59:45.061999'] = 1

You can even change all the values for a given time: ::

    >>> my_timeseries.data['xrsa']['2012-06-01 00:00'] = 1

Note, you will need to be careful to consider units when modifying the
TimeSeries data directly. For further details about editing Pandas DataFames you
can read the `pandas documentation website <https://pandas.pydata.org/pandas-docs/stable/>`_.

Additionally the TimeSeries provides the `~sunpy.timeseries.GenericTimeSeries.add_column`
method which will either add a new column or update a current column if the
colname is already present. This can take numpy array or preferably an Astropy
`~astropy.units.quantity.Quantity` value.  For example: ::

    >>> values = u.Quantity(my_timeseries.data['xrsa'].values[:-2], my_timeseries.units['xrsa']) * 20.5
    >>> my_timeseries.add_column('new col', values)
    <sunpy.timeseries.sources.goes.XRSTimeSeries object at ...>

Note that the values will be converted into the column units if an Astropy
`~astropy.units.quantity.Quantity` is given. Caution should be taken when adding
a new column because this column won't have any associated MetaData entry,
similarly if you use an array of values it won't add an entry into the units
`~collections.OrderedDict`.

5.2 Truncating a TimeSeries
===========================

It is often useful to truncate an existing TimeSeries object to retain a
specific time range.  This is easily achieved by using the `~sunpy.timeseries.GenericTimeSeries.truncate`
method. For example, to trim our GOES data into a period of interest use: ::

    >>> from sunpy.time import TimeRange
    >>> tr = TimeRange('2012-06-01 05:00','2012-06-01 06:30')
    >>> my_timeseries_trunc = my_timeseries.truncate(tr)

This takes a number of different arguments, such as the start and end dates (as
datetime or string objects) or a `~sunpy.time.TimeRange` as used above. Note
that the truncated TimeSeries will have a truncated `~sunpy.timeseries.TimeSeriesMetaData`
object, which may include dropping metadata entries for data totally cut out
from the TimeSeries.  If you want to truncate using slice-like values you can,
for example taking every 2nd value from 0 to 10000 can be done using: ::

    >>> my_timeseries_trunc = my_timeseries.truncate(0,100000,2)

Caution should be used when removing values from the data manually, the
TimeSeries can't guarantee Astropy units are correctly preserved when you
interact with the data directly.

5.3 Down and Up Sampling a TimeSeries Using Pandas
==================================================

Because the data is stored in a Pandas `~pandas.DataFrame` object you
can manipulate it using normal Pandas methods, such as the `~pandas.DataFrame.resample`
method.  To downsample you can use: ::

    >>> downsampled_dataframe = my_timeseries_trunc.data.resample('10T').mean()

Note, here ``10T`` means sample every 10 minutes and 'mean' is the method used
to combine the data. Alternatively the sum method is often used.
You can also upsample, such as: ::

    >>> upsampled_data = my_timeseries_trunc.data.resample('30S').ffill()

Note, here we upsample to 30 second intervals using ``30S`` and use the pandas
fill-forward method. Alternatively the back-fill method could be used.  Caution
should be used when resampling the data, the TimeSeries can't guarantee Astropy
Units are correctly preserved when you interact with the data directly.

5.4 Concatenating TimeSeries
============================

It's common to want to combine a number of TimeSeries together into a single
TimeSeries.  In the simplest scenario this is to combine data from a single
source over several time ranges, for example if you wanted to combine the daily
GOES data to get a week or more of constant data in one TimeSeries.  This can be
performed using the TimeSeries factory with the ``concatenate=True``
keyword argument: ::

    >>> concatenated_timeseries = sunpy.timeseries.TimeSeries(filepath1, filepath2, source='XRS', concatenate=True)  # doctest: +SKIP

Note, you can list any number of files, or a folder or use a glob to select the
input files to be concatenated.  It is possible to concatenate two TimeSeries
after creating them with the factory using the `~sunpy.timeseries.GenericTimeSeries.concatenate`
method.  For example: ::

    >>> concatenated_timeseries = goes_timeseries_1.concatenate(goes_timeseries_2)  # doctest: +SKIP

This will result in a TimeSeries identical to if you used the factory to create
it in one step.  A limitation of the TimeSeries class is that often it is not
easy to determine the source observatory/instrument of a file, generally
because the file formats used vary depending on the scientific working groups,
thus some sources need to be explicitly stated (as a keyword argument) and so it
is not possible to concatenate files from multiple sources with the factory.
To do this you can still use the `~sunpy.timeseries.GenericTimeSeries.concatenate`
method, which will create a new TimeSeries with all the rows and columns of the
source and concatenated TimeSeries in one: ::

    >>> concatenated_timeseries = goes_timeseries.concatenate(eve_timeseries)  # doctest: +SKIP

Note that the more complex `~sunpy.timeseries.TimeSeriesMetaData`
object now has 2 entries and shows details on both: ::

    >>> concatenated_timeseries.meta  # doctest: +SKIP

The metadata object is described in more detail in the next section.


5.5 Creating an Astropy Table from a TimeSeries
===============================================

If you want to take the data from your TimeSeries and use it as a `~astropy.table.Table`
this can be done using the `~sunpy.timeseries.GenericTimeSeries.to_table`
method.  For example: ::

    >>> table = my_timeseries_trunc.to_table()

Note that this `~astropy.table.Table` will contain a mixin column for
containing the Astropy `~astropy.time.Time` object representing the index,
it will also add the relevant units to the columns. One of the most useful
reasons for doing this is that Astropy `~sunpy.timeseries.GenericTimeSeries.to_table`
objects have some very nice options for viewing the data, including the basic
console view: ::

    >>> table
    <Table length=21089>
                 date               xrsa     xrsb
                                   W / m2   W / m2
            datetime64[ns]        float32  float32
    ----------------------------- ------- ----------
    2011-06-06T23:59:59.961999000     0.1 1.8871e-07
    2011-06-07T00:00:04.058999000   1e-09 1.8609e-07
    2011-06-07T00:00:08.151999000   1e-09 1.8609e-07
    2011-06-07T00:00:12.248999000   1e-09 1.8609e-07
    2011-06-07T00:00:16.344999000   1e-09 1.8084e-07
    2011-06-07T00:00:20.441999000   1e-09 1.8084e-07
    2011-06-07T00:00:24.534999000   1e-09 1.8084e-07
    2011-06-07T00:00:28.631999000   1e-09 1.8346e-07
    2011-06-07T00:00:32.728999000   1e-09 1.8346e-07
                              ...     ...        ...
    2011-06-07T23:59:20.768999000   1e-09  1.651e-07
    2011-06-07T23:59:24.864999000   1e-09 1.5985e-07
    2011-06-07T23:59:28.961999000   1e-09 1.5985e-07
    2011-06-07T23:59:33.058999000   1e-09 1.6248e-07
    2011-06-07T23:59:37.151999000   1e-09 1.6248e-07
    2011-06-07T23:59:41.248999000   1e-09 1.5985e-07
    2011-06-07T23:59:45.344999000   1e-09 1.5723e-07
    2011-06-07T23:59:49.441999000   1e-09 1.6248e-07
    2011-06-07T23:59:53.538999000   1e-09 1.5985e-07
    2011-06-07T23:59:57.631999000   1e-09 1.5985e-07

and the more sophisticated browser view using the `~astropy.table.Table.show_in_browser`
method: ::

    >>> table.show_in_browser(jsviewer=True)  # doctest: +SKIP

For further details about editing Astropy tables you can read the `astropy
documentation website <https://docs.astropy.org/en/stable/table/>`_.


6. A Detailed Look at the Metadata
==================================

TimeSeries store metadata in a `~sunpy.timeseries.TimeSeriesMetaData`
object, this object is designed to be able to store multiple basic `~sunpy.util.metadata.MetaDict`
(case-insensitive ordered dictionary) objects and able to identify the relevant
metadata for a given cell in the data. This enables a single TimeSeries to be
created by combining/concatenating multiple TimeSeries source files together
into one and to keep a reliable track of all the metadata relevant to each cell,
column or row.  The metadata can be accessed by: ::

    >>> meta = my_timeseries.meta

You can easily get an overview of the metadata, this will show you a basic
representation of the metadata entries that are relevant to this TimeSeries. ::

    >>> meta
    |-------------------------------------------------------------------------------------------------|
    |TimeRange                  | Columns         | Meta                                              |
    |-------------------------------------------------------------------------------------------------|
    |2011-06-06 23:59:59.961999 | xrsa            | simple: True                                      |
    |            to             | xrsb            | bitpix: 8                                         |
    |2011-06-07 23:59:57.631999 |                 | naxis: 0                                          |
    |                           |                 | extend: True                                      |
    |                           |                 | date: 26/06/2012                                  |
    |                           |                 | numext: 3                                         |
    |                           |                 | telescop: GOES 15                                 |
    |                           |                 | instrume: X-ray Detector                          |
    |                           |                 | object: Sun                                       |
    |                           |                 | origin: SDAC/GSFC                                 |
    |                           |                 | ...                                               |
    |-------------------------------------------------------------------------------------------------|
    <BLANKLINE>

The data within a `~sunpy.timeseries.TimeSeriesMetaData` object is
stored as a list of tuples, each tuple representing the metadata from a source
file or timeseries. The tuple will contain a `~sunpy.time.TimeRange` telling us
which rows the metadata applies to, a list of column name strings for which the
metadata applies to and finally a `~sunpy.util.metadata.MetaDict` object for
storing the key/value pairs of the metadata itself.  Each time a TimeSeries is
concatenated to the original a new set of rows and/or columns will be added to
the `~pandas.DataFrame` and a new entry will be added into the
metadata.  Note that entries are ordered chronologically based on
`~sunpy.time.timerange.TimeRange.start` and generally it's expected that no two
TimeSeries will overlap on both columns and time range.  For example it is not
good practice for alternate row values in a single column to be relevant to
different metadata entries as this would make it impossible to uniquely identify
the metadata relevant to each cell.

If you want the string that's printed then you can use the
`~sunpy.timeseries.TimeSeriesMetaData.to_string` method.  This has the
advantage of having optional keyword arguments that allows you to set the depth
(number of rows for each entry) and width (total number of characters wide)
to better fit your output.  For example: ::

    >>> meta_str = meta.to_string(depth = 20, width=99)

Similar to the TimeSeries, the metadata has some properties for convenient
access to the global metadata details, including
`~sunpy.timeseries.TimeSeriesMetaData.columns` as a list of
strings,  and `~sunpy.timeseries.TimeSeriesMetaData.time_range` of the data.
Beyond this, there are properties to get lists of details for all the entries in
the `~sunpy.timeseries.TimeSeriesMetaData` object, including
`~sunpy.timeseries.TimeSeriesMetaData.timeranges`,
`~sunpy.timeseries.TimeSeriesMetaData.columns` (as a list of string
column names) and `~sunpy.timeseries.TimeSeriesMetaData.metas`.
Similar to TimeSeries objects you can `~sunpy.timeseries.TimeSeriesMetaData.concatenate`
`~sunpy.timeseries.TimeSeriesMetaData`
objects, but generally you won't need to do this as it is done automatically
when actioned on the TimeSeries.
Note that when truncating a `~sunpy.timeseries.TimeSeriesMetaData`
object you will remove any entries outside of the given `~sunpy.time.TimeRange`.
You can also `~sunpy.timeseries.TimeSeriesMetaData.append` a new entry
(as a tuple or list), which will add the entry in the correct chronological
position.  It is frequently necessary to locate the metadata for a given column,
row or cell which can be uniquely identified by both, to do this you can use the
`~sunpy.timeseries.TimeSeriesMetaData.find` method, by adding colname
and/or time/row keyword arguments you get a `~sunpy.timeseries.TimeSeriesMetaData`
object returned which contains only the relevant entries. You can then use the
`~sunpy.timeseries.TimeSeriesMetaData.metas` property to get a list of
just the relevant `~sunpy.util.metadata.MetaDict` objects.  For example: ::

    >>> tsmd_return = my_timeseries.meta.find(colname='xrsa', time='2012-06-01 00:00:33.904999')
    >>> tsmd_return.metas
    []

Note, the colname and time filters are optional, but omitting both filters just
returns an identical `~sunpy.timeseries.TimeSeriesMetaData` object to
the TimeSeries original. A common use case for the metadata is to find out the
instrument/s that gathered the data and in this case you can use the
`~sunpy.timeseries.TimeSeriesMetaData.get` method.  This method takes a
single key string or list of key strings with the optional filters and will
search for any matching values. This method returns another `~sunpy.timeseries.TimeSeriesMetaData`
object, but removes all unwanted key/value pairs.  The result can be converted
into a simple list of strings using the `~sunpy.timeseries.TimeSeriesMetaData.values`
method: ::

    >>> tsmd_return = my_timeseries.meta.get('telescop', colname='xrsa')
    >>> tsmd_return.values()
    ['GOES 15']

Note `~sunpy.timeseries.TimeSeriesMetaData.values` removes duplicate
strings and sorts the returned list.  You can update the values for these
entries efficiently using the `~sunpy.timeseries.TimeSeriesMetaData.update`
method which takes a dictionary argument and updates the values to each of the
dictionaries that match the given colname and time filters, for example: ::

    >>> my_timeseries.meta.update({'telescop': 'G15'}, colname='xrsa', overwrite=True)

Here we have to specify the overwrite=False keyword parameter to allow us to
overwrite values for keys already present in the `~sunpy.util.metadata.MetaDict`
objects, this helps protect the integrity of the original metadata and without
this set (or with it set to False) you can still add new key/value pairs.
Note that the `~sunpy.util.metadata.MetaDict` objects are both case-insensitive
for key strings and have ordered entries, where possible the order is preserved
when updating values.
****
Maps
****

Maps in sunpy are are 2-dimensional data associated with a coordinate system. In
this guide, we will cover some of the basic functionality of maps. Once you've
read through this guide check out :doc:`/code_ref/map` for a more thorough look
at sunpy maps. There you can see what instruments are currently supported or you
can access the code reference for each instrument-specific map subclass.

Creating maps
=============
To make things easy, sunpy can download several example files which are used
throughout the docs. These files have names like
``sunpy.data.sample.AIA_171_IMAGE`` and ``sunpy.data.sample.RHESSI_IMAGE``. To
create a `~sunpy.map.Map` from the the sample AIA image
type the following into your Python shell::

    >>> import sunpy
    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA

    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA

The variable my_map is an `~sunpy.map.sources.AIAMap` object. To create one from a
local FITS file try the following::

    >>> my_map = sunpy.map.Map('/mydirectory/mymap.fits')   # doctest: +SKIP

sunpy should automatically detects the type of file (e.g. FITS), what instrument it is
associated with (e.g. AIA, EIT, LASCO) and will automatically look in the
appropriate places for the FITS keywords it needs to interpret the coordinate
system. If the type of FITS file is not recognized then sunpy will try some
default FITS keywords and return a `~sunpy.map.GenericMap` but results
may vary. sunpy can also create maps from the jpg2000 files from
`helioviewer.org <https://helioviewer.org/>`_.

Creating Custom Maps
====================
It is also possible to create maps using custom data (e.g. from a simulation or an observation
from a data source that is not explicitly supported in sunpy.) To do this you need to provide
`sunpy.map.Map` with both the data array as well as appropriate
meta information. The meta information is important as it informs the `sunpy.map.Map`
of the correct coordinate information associated with the data array. The meta information should be provided to
`sunpy.map.Map` in the form of a header as a `dict` or `~sunpy.util.MetaDict`.

The keys that are required for the header information follows the `FITS standard <https://fits.gsfc.nasa.gov/fits_dictionary.html>`_. sunpy now provides a map header helper function to assist the user in creating a header that contains the correct meta information
to generate a `sunpy.map.Map`.

The helper functionality includes a `~sunpy.map.meta_keywords` function
that will return a `dict` of all the current meta keywords and their descriptions currently used by
`sunpy.map.Map` to make a map::

    >>> from sunpy.map import meta_keywords

    >>> meta_keywords() # doctest: +SKIP
    {'cunit1': 'Units of the coordinate increments along naxis1 e.g. arcsec **required',
     'cunit2': 'Units of the coordinate increments along naxis2 e.g. arcsec **required',
     'crval1': 'Coordinate value at reference point on naxis1 **required'
     ...

There is also functionality also includes a utility function
`~sunpy.map.make_fitswcs_header` that will return a header with the
appropiate FITS keywords once the map data array and an `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames`
is passed. The `astropy.coordinates.SkyCoord` is defined by the user, and contains information on the reference frame,
reference coordinate and observer location. The function returns a `sunpy.util.MetaDict`.
The `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames` must contain an observation time.

The `~sunpy.map.make_fitswcs_header` function also takes optional keywords arguments including ``reference_pixel`` and ``scale`` which describe the pixel coordinate at the reference coordinate (defined by the `~astropy.coordinates.SkyCoord`) and the spatial scale of the pixels, respectively. If neither of these are given their values default to the center of the data array and 1 arcsec, respectively.

Here's an example of creating a header from some generic data and an `astropy.coordinates.SkyCoord`::


    >>> import numpy as np
    >>> import astropy.units as u
    >>> from sunpy.coordinates import frames
    >>> from astropy.coordinates import SkyCoord

    >>> data = np.arange(0,100).reshape(10,10)
    >>> coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime = '2013-10-28', observer = 'earth', frame = frames.Helioprojective)
    >>> header = sunpy.map.make_fitswcs_header(data, coord)
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 5.5
    crpix2: 5.5
    cdelt1: 1.0
    cdelt2: 1.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768


From this we can see now that the function returned a `sunpy.util.MetaDict` that populated
the standard FITS keywords with information provided by the passed `astropy.coordinates.SkyCoord`,
and the data array. Since the ``reference_pixel`` and keywords were not passed in the example above, the
values of ``crpix`` and ``cdelt`` were set to the default values.

These keywords can be passed to the function in the form of an `astropy.units.Quantity` with associated units.
Here's another example of passing ``reference_pixel`` and ``scale`` to the function::

    >>> header = sunpy.map.make_fitswcs_header(data, coord,
    ...                                        reference_pixel=u.Quantity([5, 5]*u.pixel),
    ...                                        scale=u.Quantity([2, 2] *u.arcsec/u.pixel))
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 6.0
    crpix2: 6.0
    cdelt1: 2.0
    cdelt2: 2.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768

As we can see, a list of WCS and observer meta information is contained within the generated headers,
however we may want to include other meta information including the observatory name, the wavelength and
waveunit of the observation. Any of the keywords listed in ``header_helper.meta_keywords`` can be passed
to the `~sunpy.map.make_fitswcs_header` and will then populate the returned MetaDict header.
Furthermore, the following observation keywords can be passed to the `~sunpy.map.make_fitswcs_header`
function and will be translated to the FITS standard: ``observtory``, ``instrument``,``telescope``, ``wavelength``, ``exposure``.

An example of creating a header with these additional keywords::

    >>> header = sunpy.map.make_fitswcs_header(data, coord,
    ...                                        reference_pixel = u.Quantity([5, 5]*u.pixel),
    ...                                        scale = u.Quantity([2, 2] *u.arcsec/u.pixel),
    ...                                        telescope = 'Test case', instrument = 'UV detector',
    ...                                        wavelength = 1000*u.angstrom)
    >>> header  # doctest: +SKIP
    MetaDict([('wcsaxes', 2),
          ('crpix1', 5.0),
          ('crpix2', 5.0),
          ('cdelt1', <Quantity 2. arcsec2 / pix2>),
          ('cdelt2', <Quantity 2. arcsec2 / pix2>),
          ('cunit1', Unit("arcsec")),
          ('cunit2', Unit("arcsec")),
          ('ctype1', 'HPLN-TAN'),
          ('ctype2', 'HPLT-TAN'),
          ('crval1', 0.0),
          ('crval2', 0.0),
          ...
          ('date-obs', '2013-10-28T00:00:00.000'),
          ('hgln_obs', 0.0),
          ('hglt_obs', 4.7711570596394015),
          ('dsun_obs', 148644585949.4918),
          ('rsun_ref', 695700.0),
          ('rsun_obs', 965.3829548285768),
          ('instrume', 'Test case'),
          ('wavelnth', 1000),
          ('detector', 'UV detector'),
          ('waveunit', 'angstrom')])

From these header MetaDict's that are generated, we can now create a custom map::

    >>> my_map = sunpy.map.Map(data, header) # doctest: +SKIP
    >>> my_map.peek() # doctest: +SKIP

Inspecting maps
===============
A map contains a number of data-associated attributes. To get a quick look at
your map simply type::

    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> my_map  # doctest: +REMOTE_DATA
    <sunpy.map.sources.sdo.AIAMap object at ...>
    SunPy Map
    ---------
    Observatory:                 SDO
    Instrument:          AIA 3
    Detector:            AIA
    Measurement:                 171.0 Angstrom
    Wavelength:          171.0 Angstrom
    Observation Date:    2011-06-07 06:33:02
    Exposure Time:               0.234256 s
    Dimension:           [1024. 1024.] pix
    Coordinate System:   helioprojective
    Scale:                       [2.402792 2.402792] arcsec / pix
    Reference Pixel:     [511.5 511.5] pix
    Reference Coord:     [3.22309951 1.38578135] arcsec
    array([[ -95.92475  ,    7.076416 ,   -1.9656711, ..., -127.96519  ,
            -127.96519  , -127.96519  ],
           [ -96.97533  ,   -5.1167884,    0.       , ...,  -98.924576 ,
            -104.04137  , -127.919716 ],
           [ -93.99607  ,    1.0189276,   -4.0757103, ...,   -5.094638 ,
             -37.95505  , -127.87541  ],
           ...,
           [-128.01454  , -128.01454  , -128.01454  , ..., -128.01454  ,
            -128.01454  , -128.01454  ],
           [-127.899666 , -127.899666 , -127.899666 , ..., -127.899666 ,
            -127.899666 , -127.899666 ],
           [-128.03072  , -128.03072  , -128.03072  , ..., -128.03072  ,
            -128.03072  , -128.03072  ]], dtype=float32)

This will show a representation of the data as well as some of its associated
attributes. A number of other attributes are also available, for example the
`~sunpy.map.GenericMap.date`, `~sunpy.map.GenericMap.exposure_time`,
`~sunpy.map.GenericMap.center` and others (see `~sunpy.map.GenericMap`)::

    >>> map_date = my_map.date  # doctest: +REMOTE_DATA
    >>> map_exptime = my_map.exposure_time  # doctest: +REMOTE_DATA
    >>> map_center = my_map.center  # doctest: +REMOTE_DATA

To get a list of all of the attributes check the documentation by typing::

    >>> help(my_map)  # doctest: +SKIP

Many attributes and functions of the map classes accept and return
`~astropy.units.quantity.Quantity` or `~astropy.coordinates.SkyCoord` objects,
please refer to :ref:`units-coordinates-sunpy` for more details.

The meta data for the map is accessed by ::

    >>> header = my_map.meta  # doctest: +REMOTE_DATA

This references the meta data dictionary with the header information as read
from the source file.

Getting at the data
===================
The data in a sunpy Map object is accessible through the
`~sunpy.map.GenericMap.data` attribute.  The data is implemented as a
NumPy `~numpy.ndarray`, so for example, to get
the 0th element in the array ::

    >>> my_map.data[0, 0]  # doctest: +REMOTE_DATA
    -95.92475
    >>> my_map.data[0][0]  # doctest: +REMOTE_DATA
    -95.92475

One important fact to remember is that the first
index is for the y direction while the second index is for the x direction.
For more information about indexing please refer to the
`Numpy documentation <https://docs.scipy.org/doc/numpy-dev/user/quickstart.html#indexing-slicing-and-iterating>`_.

Data attributes like `~numpy.ndarray.dtype` and
`~sunpy.map.GenericMap.dimensions` are accessible through
the SunPyGenericMap object ::

    >>> my_map.dimensions  # doctest: +REMOTE_DATA
    PixelPair(x=<Quantity 1024. pix>, y=<Quantity 1024. pix>)
    >>> my_map.dtype  # doctest: +REMOTE_DATA
    dtype('float32')

Here the dimensions attribute is similar to the `~numpy.ndarray.shape`
attribute, however returning an `~astropy.units.quantity.Quantity`.

If you'd like to use the data in a sunpy `~sunpy.map.GenericMap` object
elsewhere, you can use either of the following::

    >>> var = my_map.data  # doctest: +REMOTE_DATA
    >>> var = my_map.data.copy()  # doctest: +REMOTE_DATA

Python makes use of pointers so if you want to alter the data and keep the
original data in the map intact make sure to copy it.

To create a complete copy of a Map object that is entirely independent of the original,
use the built-in `copy.deepcopy`
method, like so::

    >>> import copy   # doctest: +REMOTE_DATA
    >>> my_map_deepcopy = copy.deepcopy(my_map)   # doctest: +REMOTE_DATA

A deepcopy ensures that any changes in the original Map object are not reflected in the
copied object and vice versa. Note that this is different from simply copying the data of
the Map object - this copies all of the other attributes and methods as well.

Some basic statistical functions on the data array are also passed through to Map
objects::

    >>> my_map.min()  # doctest: +REMOTE_DATA
    -129.78036
    >>> my_map.max()  # doctest: +REMOTE_DATA
    192130.17
    >>> my_map.mean()  # doctest: +REMOTE_DATA
    427.02252

but you can also access all the other `~numpy.ndarray` functions and attributes
by accessing the data array directly. For example::

    >>> my_map.data.std()  # doctest: +REMOTE_DATA
    826.41016

Plotting
========
As is true of all of the sunpy data objects, the sunpy `~sunpy.map.GenericMap`
object (and all of its instrument-specific sub-classes) has its
own built-in plot methods so that it is easy to quickly view your map.
To create a plot just type::

    >>> my_map.peek()   # doctest: +SKIP

This will open a matplotlib plot on your screen.
In addition, to enable users to modify the plot it is possible to grab the
matplotlib axes object by using the `~sunpy.map.GenericMap.plot()` command.
This makes it possible to use the sunpy plot as the foundation for a
more complicated figure. For a bit more information about this and some
examples see :ref:`plotting`.

.. note::

   If the `astropy.visualization.wcsaxes` package is not used (it is used by
   default) the `~sunpy.map.GenericMap.plot()` and
   `~sunpy.map.GenericMap.peek()` methods assume that the data is not rotated,
   i.e. the solar y axis is oriented with the columns of the array. If this
   condition is not met (in the metadata), when the map is plotted a warning
   will be issued. You can create an oriented map by using
   `~sunpy.map.GenericMap.rotate()` before you plot the Map.

Plotting Keywords
-----------------

For Map `~matplotlib.pyplot.imshow` does most of the heavy
lifting in the background while sunpy makes a number of choices for you so that
you don't have to (e.g. colortable, plot title). Changing these defaults
is made possible through two simple interfaces. You can pass any
`~matplotlib.pyplot.imshow` keyword into
the plot command to override the defaults for that particular plot. The following
plot changes the default AIA color table to use an inverse Grey color table.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    fig = plt.figure()
    smap.plot(cmap=plt.cm.Greys_r)
    plt.colorbar()
    plt.show()

You can view or make changes to the default settings through the ``sunpy.map.GenericMap.plot_settings``
dictionary. In the following example we change the title of the plot by changing the
``sunpy.map.GenericMap.plot_settings`` property.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    smap.plot_settings['title'] = "My Second sunpy Plot"
    smap.plot_settings['cmap'] = plt.cm.Blues_r
    fig = plt.figure()
    smap.plot()
    plt.colorbar()
    plt.show()


Colormaps and Normalization
---------------------------

Image data is generally shown in false color in order to better identify it or
to better visualize structures in the image. Matplotlib handles this colormapping
process through the `~matplotlib.colors` module. This process involves two steps:
the data array is first mapped onto the range 0-1 using an instance of
`~matplotlib.colors.Normalize` or a subclass; then this number is mapped to a
color using an instance of a subclass of a `~matplotlib.colors.Colormap`.

sunpy provides the colormaps for each mission as defined by the mission teams.
The Map object chooses the appropriate colormap for you when it is created as
long as it recognizes the instrument. To see what colormaps are available::

    >>> import sunpy.visualization.colormaps as cm
    >>> cm.cmlist.keys()
    dict_keys(['goes-rsuvi94', 'goes-rsuvi131', 'goes-rsuvi171', 'goes-rsuvi195',
    'goes-rsuvi284', 'goes-rsuvi304', 'sdoaia94', 'sdoaia131', 'sdoaia171',
    ...

The sunpy colormaps are registered with matplotlib so you can grab them like
you would any other colormap::

    >>> import matplotlib.pyplot as plt
    >>> import sunpy.visualization.colormaps

You need to import `sunpy.visualization.colormaps` or `sunpy.map` for this to work::

    >>> cmap = plt.get_cmap('sdoaia171')


The following plot shows off all of the colormaps.

.. plot::

    import matplotlib.pyplot as plt
    import sunpy.visualization.colormaps as cm
    cm.show_colormaps()

These can be used with the standard commands to change the colormap. So for
example if you wanted to plot an AIA image but use an EIT colormap, you would
do so as follows.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = plt.get_cmap('sohoeit171')

    fig = plt.figure()
    smap.plot(cmap=cmap)
    plt.colorbar()
    plt.show()

or you can just change the colormap for the map itself as follows::

    >>> smap.plot_settings['cmap'] = plt.get_cmap('sohoeit171')  # doctest: +SKIP

The normalization is also set automatically and is chosen so that all the
data from minimum to maximum is displayed as best as possible for most cases.
This means that it is never necessary to touch the data such as applying a function
such sqrt or log to the data to make your plot look good.
There are many normalizations available from matplotlib such as `~matplotlib.colors.LogNorm`. Other
`more exotic normalizations <https://docs.astropy.org/en/stable/visualization/index.html>`_ are also
made available from Astropy.  Just like the colormap the default normalization
can be changed through the plot_settings dictionary or directly for the individual
plot by passing a keyword argument. The following example shows the difference between
a linear and logarithmic normalization on an AIA image.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure(figsize=(4, 9))

    ax1 = fig.add_subplot(2, 1, 1, projection=smap)
    smap.plot(norm=colors.Normalize(), title='Linear normalization')
    plt.colorbar()

    ax2 = fig.add_subplot(2, 1, 2, projection=smap)
    smap.plot(norm=colors.LogNorm(), title='Logarithmic normalization')
    plt.colorbar()

    plt.show()

Note how the color in the colorbar does not change since these two maps share
the same colormap while the data values associated with each color do because
the normalization is different.

Masking and Clipping Data
=========================
It is often necessary for the purposes of display or otherwise to ignore certain
data in an image. For example, a large data value could be due to
cosmic ray hits and should be ignored. The most straightforward way to ignore
this kind of data in plots without altering the data is to clip it. This can be achieved
very easily by using the ``clip_interval`` keyword. For example::

    >>> import astropy.units as u
    >>> smap.plot(clip_interval=(1, 99.5)*u.percent)  #doctest: +SKIP

This clips out the dimmest 1% of pixels and the brightest 0.5% of pixels.  With those outlier
pixels clipped, the resulting image makes better use of the full range of colors.
If you'd like to see what areas of your images got clipped, you can modify the colormap::

    >>> cmap = map.cmap  # doctest: +SKIP
    >>> cmap.set_over('blue')  # doctest: +SKIP
    >>> cmap.set_under('green')  # doctest: +SKIP

This will color the areas above and below in red and green respectively
(similar to this `example <https://matplotlib.org/examples/pylab_examples/image_masked.html>`_).
You can use the following colorbar command to display these choices::

    >>> plt.colorbar(extend='both')   # doctest: +SKIP

Here is an example of this put to use on an AIA image.

.. plot::
    :include-source:

    import astropy.units as u
    import matplotlib.pyplot as plt

    import sunpy.map
    import sunpy.data.sample

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = smap.cmap.copy()
    cmap.set_over('blue')
    cmap.set_under('green')

    fig = plt.figure(figsize=(12, 4))

    ax1 = fig.add_subplot(1, 2, 1, projection=smap)
    smap.plot(title='Without clipping')
    plt.colorbar()

    ax2 = fig.add_subplot(1, 2, 2, projection=smap)
    smap.plot(clip_interval=(1, 99.5)*u.percent, title='With clipping')
    plt.colorbar(extend='both')

    plt.show()

Another approach to clipping data is to specify explicit values for the minimum and maximum pixel
values using the plotting keywords ``vmin`` and ``vmax``.

Clipping excludes data that has extreme values, but there can be other forms of bad data.
A mask is a boolean
array and so can give you much more fine-grained control over what is not being
displayed.  A `~numpy.ma.MaskedArray`
is a subclass of a numpy array so it has all of the same properties with the
addition of an associated boolean array which holds the mask.
See `this example <https://docs.sunpy.org/en/stable/generated/gallery/computer_vision_techniques/mask_disk.html>`_ in our gallery.

.. the following is a good example which could be fixed and added later
.. The following plot achieves the same goal as above but using a mask instead of clipping.

..    import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    cmap = smap.plot_settings['cmap']
    cmap.set_bad('blue', 1.0)
    smap = sunpy.map.Map('/Users/schriste/Downloads/old downloads/foxsi_ar_data/ssw_cutout_20121030_153001_AIA_94_.fts')
    smap.mask =
    smap.plot()
    plt.colorbar(extend='both')
    plt.show()

.. Hinode XRT image. By inspecting the maximum versus the mean and standard deviation, it is clear that there are some overly bright pixels. This is likely due to cosmic ray hits which is throwing off the default plot making it too dark to see the solar emission.

.. .. plot::

..    import sunpy.map
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map('/Users/schriste/Desktop/sunpy_test_img/XRT20141211_184221.9.fits')
    fig = plt.figure()
    smap.plot()
    txt = r"min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                       max=int(smap.max()),
                                                                       mean=int(smap.mean()),
                                                                       std=int(smap.std()))
    plt.text(-600, 1500, txt, color='white')
    plt.colorbar()
    plt.show()

.. Let's address this by clipping the largest values (in this case everything above 3 sigma). The following plot shows the result of this operation.

.. .. plot::

..     import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    cmap = smap.plot_settings['cmap']
    cmap.set_over('green', 1.0)
    cmap.set_under('purple', 1.0)
    norm = colors.Normalize(vmin=smap.min(), vmax=smap.mean() + 3 *smap.std())
    smap = sunpy.map.Map('/Users/schriste/Desktop/sunpy_test_img/XRT20141211_184221.9.fits')
    smap.plot(norm=norm)
    plt.colorbar(extend='both')
    plt.show()

.. This makes it very visible that there are a number of hot pixels mostly concentrated in the upper half of this image. Now let's address this problem with masking instead of clipping.

.. .. plot::

..     import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy.ma
    smap = sunpy.map.Map('/Users/schriste/Desktop/sunpy_test_img/XRT20141211_184221.9.fits')
    cmap = smap.plot_settings['cmap']
    cmap.set_bad('blue', 1.0)
    smap.data = numpy.ma.masked_greater(smap.data, smap.mean() + 3 *smap.std())
    txt = r"min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                       max=int(smap.max()),
                                                                       mean=int(smap.mean()),
                                                                       std=int(smap.std()))
    plt.text(-600, 1500, txt, color='white')
    norm = colors.Normalize()
    smap.plot(norm = norm)
    plt.colorbar(extend='both')

.. This plot shows a very similar effect to clipping but note that the array properties such as max and min have changed. That's because numpy is now ignoring those masked values. With a masked array
.. (compared to clipping) we can go ahead and make more detailed masking operations so that we are not masking the emission from the bright solar sources. The next plot masks only those bright pixels in the upper area of the plot leaving the bright solar sources which are concentrated in the lower part of the plot intact.

.. .. plot::

..     import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy.ma
    file = '/Users/schriste/Downloads/old downloads/foxsi_ar_data/sXRT20141211_184221.9.fits'
    smap = sunpy.map.Map(file)
    cmap = smap.plot_settings['cmap']
    cmap.set_bad('blue', 1.0)
    smap.data = numpy.ma.masked_greater(smap.data, smap.mean() + 3 *smap.std())
    smap.data.mask[0:250,:] = False
    txt = r"min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                       max=int(smap.max()),
                                                                       mean=int(smap.mean()),
                                                                       std=int(smap.std()))
    plt.text(-600, 1500, txt, color='white')
    norm = colors.Normalize()
    smap.plot(norm = norm)
    plt.colorbar(extend='both')


Composite Maps and Overlaying Maps
==================================

The `~sunpy.map.Map` method described above can also handle a list of maps. If a series of maps
are supplied as inputs, `~sunpy.map.Map` will return a list of maps as the output.  However,
if the 'composite' keyword is set to True, then a `~sunpy.map.CompositeMap` object is
returned.  This is useful if the maps are of a different type (e.g. different
instruments).  For example, to create a simple composite map::

    >>> my_maps = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE, sunpy.data.sample.RHESSI_IMAGE, composite=True)  # doctest: +REMOTE_DATA

A `~sunpy.map.CompositeMap` is different from a regular sunpy `~sunpy.map.GenericMap` object and therefore
different associated methods. To list which maps are part of your composite map use::

    >>> my_maps.list_maps()  # doctest: +REMOTE_DATA
    [<class 'sunpy.map.sources.soho.EITMap'>, <class 'sunpy.map.sources.rhessi.RHESSIMap'>]

The following code adds a new map (which must be instantiated first), sets
its transparency to 25%, turns on contours from 50% to 90% for the second
map, and then plots the result.

.. plot::
    :include-source:

    import sunpy.data.sample
    import sunpy.map
    import matplotlib.pyplot as plt
    my_maps = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE, sunpy.data.sample.RHESSI_IMAGE, composite=True)
    my_maps.add_map(sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE))
    my_maps.set_alpha(2, 0.5)
    my_maps.set_levels(1, [50, 60, 70, 80, 90], percent = True)
    my_maps.plot()
    plt.show()

This is not a particularly pretty plot but it shows what sunpy can do!

Working with your map
=====================
Part of the philosophy of the map object is to provide most of the basic
functionality that a scientist would want therefore a map also contains a number
of map-specific methods such as resizing a map or grabbing a subview. To get
a list of the methods available for a map type::

    >>> help(my_map)  # doctest: +SKIP

and check out the methods section!

MapSequences
============
A `~sunpy.map.MapSequence` is an ordered list of maps.  By default, the maps are ordered by
their observation date, from earlier maps to later maps. A `~sunpy.map.MapSequence` can be
created by supplying multiple existing maps::

    >>> map1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> map2 = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)  # doctest: +REMOTE_DATA
    >>> mc = sunpy.map.Map([map1, map2], sequence=True)  # doctest: +REMOTE_DATA

or by providing a directory full of image files::

    >>> mc = sunpy.map.Map('path/to/my/files/*.fits', sequence=True)   #  doctest: +SKIP

The earliest map in the MapSequence can be accessed by simply indexing the maps
list::

    >>> mc.maps[0]   # doctest: +SKIP

MapSequences can hold maps that have different shapes.  To test if all the
maps in a `~sunpy.map.MapSequence` have the same shape::

    >>> mc.all_maps_same_shape()  # doctest: +REMOTE_DATA
    True

It is often useful to return the image data in a `~sunpy.map.MapSequence` as a single
three dimensional Numpy `~numpy.ndarray`::

    >>> mc.as_array()   # doctest: +SKIP

Note that an array is returned only if all the maps have the same
shape.  If this is not true, an error (ValueError) is returned.  If all the
maps have nx pixels in the x-direction, and ny pixels in the y-direction,
and there are n maps in the MapSequence, the `~numpy.ndarray` array that is
returned has shape (ny, nx, n).  The data of the first map in the `~sunpy.map.MapSequence`
appears in the `~numpy.ndarray` in position ``[:, :, 0]``, the data of second map in
position ``[:, :, 1]``, and so on.  The order of maps in the `~sunpy.map.MapSequence` is
reproduced in the returned `~numpy.ndarray`.

The meta data from each map can be obtained using::

    >>> mc.all_meta()   # doctest: +SKIP

This returns a list of map meta objects that have the same order as
the maps in the `~sunpy.map.MapSequence`.

Coalignment of MapSequences
===========================
A typical data preparation step when dealing with time series of images is to
coalign images taken at different times so that features in different images
remain in the same place.  A common approach to this problem is
to take a representative template that contains the features you are interested
in, and match that to your images.  The location of the best match tells you
where the template is in your image.  The images are then shifted to the
location of the best match.  This aligns your images to the position of the
features in your representative template.

sunpy provides a function to coalign the maps inside the `~sunpy.map.MapSequence`.
The implementation of this functionality requires the installation of the
scikit-image library, a commonly used image processing library.
To coalign a `~sunpy.map.MapSequence`, simply import
the function and apply it to your `~sunpy.map.MapSequence`::

    >>> from sunpy.image.coalignment import mapsequence_coalign_by_match_template
    >>> coaligned = mapsequence_coalign_by_match_template(mc)  # doctest: +REMOTE_DATA

This will return a new `~sunpy.map.MapSequence`, coaligned to a template extracted from the
center of the first map in the `~sunpy.map.MapSequence`, with the map dimensions clipped as
required.  The coalignment algorithm provides many more options for handling
the coalignment of `~sunpy.map.MapSequence` type::

    >>> help(mapsequence_coalign_by_match_template)   # doctest: +SKIP

for a full list of options and functionality.

If you just want to calculate the shifts required to compensate for solar
rotation relative to the first map in the `~sunpy.map.MapSequence` without applying them, use::

    >>> from sunpy.image.coalignment import calculate_match_template_shift
    >>> shifts = calculate_match_template_shift(mc)  # doctest: +REMOTE_DATA

This is the function used to calculate the shifts in `~sunpy.map.MapSequence` coalignment
function above.  Please see `~sunpy.image.coalignment.calculate_match_template_shift` to learn more about its features.
Shifts calculated using calculate_match_template_shift can be passed directly
to the coalignment function.


Compensating for solar rotation in MapSequences
===============================================
Often a set of solar image data consists of fixing the pointing of a
field of view for some time and observing.  Features on the Sun will
rotate according to the Sun's rotation.

A typical data preparation step when dealing with time series of these
types of images is to shift the images so that features do not appear
to move across the field of view.  This requires taking in to account
the rotation of the Sun.  The Sun rotates differentially, depending on
latitude, with features at the equator moving faster than features at
the poles.

sunpy provides a function to shift images in `~sunpy.map.MapSequence` following solar
rotation.  This function shifts an image according to the solar
differential rotation calculated at the latitude of the center of the
field of view.  The image is not *differentially* rotated.  This
function is useful for de-rotating images when the effects of
differential rotation in the `~sunpy.map.MapSequence` can be ignored (for example, if
the spatial extent of the image is small, or when the duration of the
`~sunpy.map.MapSequence` is small; deciding on what 'small' means depends on your
application).

To apply this form of solar derotation to a `~sunpy.map.MapSequence`, simply import the
function and apply it to your `~sunpy.map.MapSequence`::

    >>> from sunpy.physics.solar_rotation import mapsequence_solar_derotate
    >>> derotated = mapsequence_solar_derotate(mc)  # doctest: +SKIP

For more info see `~sunpy.physics.solar_rotation.mapsequence_solar_derotate`.

If you just want to calculate the shifts required to compensate for solar
rotation relative to the first map in the `~sunpy.map.MapSequence` without applying them, use::

    >>> from sunpy.physics.solar_rotation import calculate_solar_rotate_shift
    >>> shifts = calculate_solar_rotate_shift(mc)  # doctest: +SKIP

Please consult the docstring of the `~sunpy.image.coalignment.mapsequence_coalign_by_match_template` function in order to learn about the features of this function.
*******************
Data Types in sunpy
*******************

In this section of the guide we introduce the parts of sunpy used to load and
represent solar data.

.. toctree::
    :maxdepth: 2

    maps
    timeseries
***************************************
Querying and Downloading Data from JSOC
***************************************

Joint Science Operations Center (JSOC) contains data products from the Solar Dynamics Observatory, as well as certain other missions and instruments.
These data are available from the JSOC database, which can be directly accessed by the online `JSOC interface <http://jsoc.stanford.edu/ajax/lookdata.html>`__.

sunpy's JSOC Client provides an easier interface to query for JSOC data and make export requests.
It uses `drms module <https://docs.sunpy.org/projects/drms>`_ as its backend, and exposes a similar API as the VSO Client.

There are two ways of downloading JSOC data. One way is using sunpy's unified search interface, known as ``Fido``.
``Fido`` supplies a single, easy and consistent way to to obtain most forms of solar physics data.
An alternative way to fetch data from JSOC is by using the underlying JSOC Client.
This option can be preferred when you need to separate the staging and downloading steps, which is not supported by Fido.

The JSOC stages data before you can download it, so a JSOC query is a three stage process.
First you query the JSOC for records and a table of these records is returned.
Then you can request these records to be staged for download and then you can download them.
Fido combines the last two stages into a single call to `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`.

Setup
*****

sunpy's Fido module is in `sunpy.net`.
It can be imported as follows::

    >>> from sunpy.net import Fido, attrs as a

The JSOC client handles the particulars of how the data from the data provider is downloaded to your computer.

.. warning::

    You must have an email address registered with JSOC before you are allowed to make a request.
    See `this <http://jsoc.stanford.edu/ajax/register_email.html>`__ to register your email address.

Querying the JSOC
*****************

To search for data in JSOC, your query needs at minimum, a "Series" name and a "PrimeKey".
Different PrimeKeys are supported by different Series, and you can find out the PrimeKeys supported in any Series by:

    >>> import drms
    >>> c = drms.Client()  # doctest: +REMOTE_DATA
    >>> print(c.pkeys('hmi.m_720s'))  # doctest: +REMOTE_DATA
    ['T_REC', 'CAMERA']

The most common PrimeKey, that is supported by every Series is Time, that is denoted by ``T_REC`` or ``T_OBS``.
Hence, Time can always be passed as an attribute while building a query.
Wavelength is another pre-defined attribute which is a PrimeKey.
Other PrimeKeys which need to be passed should be manually passed in
`~sunpy.net.jsoc.attrs.PrimeKey`.
This will be explained later in detail.

Constructing a Basic Query
==========================

Let's start with a very simple query.
We could ask for all ``hmi.v_45s`` series data between January 1st from 00:00 to 01:00, 2014::

    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                   a.jsoc.Series('hmi.v_45s'))  # doctest: +REMOTE_DATA

This returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing information on the available online files which fit the criteria specified by the attrs objects in the above call.
It does not download the files.

To see a summary of results of our query, simply type the name of the variable set to the Fido search, in this case, ``res``::

    >>> res  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    81 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2014.01.01_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:01:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:02:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:04:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:05:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:06:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:06:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:07:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
                        ...      ...        ...      ...     ...
    2014.01.01_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    <BLANKLINE>

Now, let's break down the arguments of ``Fido.search`` to understand better what we've done.
The first argument ``a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00')`` sets the start and end times for the query (any date/time format understood by sunpy's :ref:`parse_time function <parse-time>` can be used to specify dates and time).
The Time attribute takes UTC time, as default.
If you need to pass a Time in some other time scale, such as TAI, pass an `astropy.time.Time` object, like::

    >>> import astropy.time

Then, the Time attribute can be passed as::

    >>> a.Time(astropy.time.Time('2014-01-01T00:00:00', scale='tai'), astropy.time.Time('2014-01-01T01:00:00', scale='tai'))
    <sunpy.net.attrs.Time(2014-01-01 00:00:00.000, 2014-01-01 01:00:00.000)>

The second argument::

    >>> a.jsoc.Series('hmi.v_45s')
    <sunpy.net.jsoc.attrs.Series(hmi.v_45s: Dopplergrams with a cadence of 45 seconds) object ...>

sets the series we are looking for.

So what is going on here?
The notion is that a JSOC query has a set of attribute objects, imported as ``a.jsoc``, that are specified to construct the query.

``a.jsoc.Series()`` is compulsory to be provided in each of the jsoc queries.
Apart from this, at least one PrimeKey must be passed (generally ``a.Time()``).


Querying with other PrimeKeys
=============================

Other than Time, one other PrimeKey is supported with in-built attribute.
In case of AIA series, ``a.Wavelength()`` can be passed as a PrimeKey::

    >>> import astropy.units as u
    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                               a.jsoc.Series('aia.lev1_euv_12s'),
    ...                               a.Wavelength(304*u.AA))  # doctest: +REMOTE_DATA

Note that, only Time and Wavelength are in-built attributes here. If you need to pass any other PrimeKey, it should be passed like this::

    >>> a.jsoc.PrimeKey('HARPNUM', '4864')
    <sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('HARPNUM', '4864')

If 2 or more PrimeKeys need to be passed together::

    >>> a.jsoc.PrimeKey('HARPNUM', '4864') & a.jsoc.PrimeKey('CAMERA', '2')
    <AttrAnd([<sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('HARPNUM', '4864'), <sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('CAMERA', '2')])>

Also, note that the pre-defined primekeys, Time and Wavelength can also be passed as above, but you need to specify the exact keyword for it.
For e.g. by::

    >>> a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.PrimeKey('WAVELNTH', '161')
    (<sunpy.net.attrs.Time(2014-01-01 00:00:00.000, 2014-01-01 01:00:00.000)>, <sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('WAVELNTH', '161'))

If the correct keyword is not specified, or the passed PrimeKey is not supported by the given series, a meaningful error will be thrown, which will give you the PrimeKeys supported by that series.
Hence, by looking at the error, one can easily retry building the query with correct PrimeKeys.

Another important thing to note is that, Wavelength when passed through in-built attribute, should be passed as an astropy quantity.
Specifying spectral units in arguments is necessary or an error will be raised.
For more information on units, see `~astropy.units`.
But, when the same is passed through PrimeKey attribute, it should be passed as a string.
All other PrimeKey values passed through PrimeKey attribute, must be passed as a string.

Manually specifying keyword data to fetch
=========================================

Upon doing ``Fido.search()`` as described above, only a limited set of keywords are returned in the response object.
These default keywords are ``'DATE'``, ``'TELESCOP'``, ``'INSTRUME'``, ``'T_OBS'`` and ``'WAVELNTH'``.

If you want to get a manual set of keywords in the response object, you can pass the set of keywords using :meth:`~sunpy.net.base_client.QueryResponseTable.show` method.

    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                   a.jsoc.Series('hmi.v_45s'))  # doctest: +REMOTE_DATA
    >>> res.show('TELESCOP', 'INSTRUME', 'T_OBS')  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    81 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
    TELESCOP  INSTRUME           T_OBS
    -------- ---------- -----------------------
     SDO/HMI HMI_FRONT2 2014.01.01_00:00:37_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:01:22_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:02:07_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:02:52_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:03:37_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:04:22_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:05:07_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:05:52_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:06:37_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:07:22_TAI
         ...        ...                     ...
     SDO/HMI HMI_FRONT2 2014.01.01_00:53:07_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:53:52_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:54:37_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:55:22_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:56:07_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:56:52_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:57:37_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:58:22_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:59:07_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_00:59:52_TAI
     SDO/HMI HMI_FRONT2 2014.01.01_01:00:37_TAI
    Length = 81 rows
    <BLANKLINE>
    <BLANKLINE>

Passing an incorrect keyword won't throw an error, but the corresponding column in the table will not be displayed.

To display all of the columns, we can use ``show()`` without passing any arguments::

    >>> res.show()  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    81 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
            DATE                DATE__OBS        ... CALVER64
    -------------------- ----------------------- ... --------
    2014-01-05T17:46:02Z 2013-12-31T23:59:39.20Z ...     4370
    2014-01-05T17:47:10Z 2014-01-01T00:00:24.20Z ...     4370
    2014-01-05T17:48:18Z 2014-01-01T00:01:09.20Z ...     4370
    2014-01-05T17:49:25Z 2014-01-01T00:01:54.20Z ...     4370
    2014-01-05T17:50:34Z 2014-01-01T00:02:39.20Z ...     4370
    2014-01-05T17:51:42Z 2014-01-01T00:03:24.20Z ...     4370
    2014-01-05T17:52:50Z 2014-01-01T00:04:09.20Z ...     4370
    2014-01-05T17:53:59Z 2014-01-01T00:04:54.20Z ...     4370
    2014-01-05T17:55:08Z 2014-01-01T00:05:39.20Z ...     4370
    2014-01-05T17:56:16Z 2014-01-01T00:06:24.20Z ...     4370
                     ...                     ... ...      ...
    2014-01-05T19:05:49Z 2014-01-01T00:52:09.20Z ...     4370
    2014-01-05T17:35:43Z 2014-01-01T00:52:54.20Z ...     4370
    2014-01-05T17:36:54Z 2014-01-01T00:53:39.20Z ...     4370
    2014-01-05T17:38:01Z 2014-01-01T00:54:24.20Z ...     4370
    2014-01-05T17:39:09Z 2014-01-01T00:55:09.20Z ...     4370
    2014-01-05T17:40:17Z 2014-01-01T00:55:54.20Z ...     4370
    2014-01-05T17:41:25Z 2014-01-01T00:56:39.20Z ...     4370
    2014-01-05T17:42:33Z 2014-01-01T00:57:24.20Z ...     4370
    2014-01-05T17:43:41Z 2014-01-01T00:58:09.20Z ...     4370
    2014-01-05T17:44:52Z 2014-01-01T00:58:54.20Z ...     4370
    2014-01-05T17:46:03Z 2014-01-01T00:59:39.20Z ...     4370
    Length = 81 rows
    <BLANKLINE>
    <BLANKLINE>

Using Segments
==============
In some cases, more than 1 file are present for the same set of query.
These data are distinguished by what are called Segments.
It is necessary to specify the "Segment" which you need to download.
Providing a segment won't have any affect on the response object returned, but this will be required later, while making an export request.

A list of supported segments of a series, say ``hmi.sharp_720s`` can be obtained by::

    >>> import drms
    >>> c = drms.Client()  # doctest: +REMOTE_DATA
    >>> si = c.info('hmi.sharp_720s')  # doctest: +REMOTE_DATA
    >>> print(si.segments.index.values)  # doctest: +REMOTE_DATA
    ['magnetogram' 'bitmap' 'Dopplergram' 'continuum' 'inclination' 'azimuth'
     'field' 'vlos_mag' 'dop_width' 'eta_0' 'damping' 'src_continuum'
     'src_grad' 'alpha_mag' 'chisq' 'conv_flag' 'info_map' 'confid_map'
     'inclination_err' 'azimuth_err' 'field_err' 'vlos_err' 'alpha_err'
     'field_inclination_err' 'field_az_err' 'inclin_azimuth_err'
     'field_alpha_err' 'inclination_alpha_err' 'azimuth_alpha_err' 'disambig'
     'conf_disambig']

Also, if you provide an incorrect segment name, it will throw a meaningful error, specifying which segment values are supported by the given series::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.sharp_720s'),
    ...             a.jsoc.Segment('image'))  # doctest: +REMOTE_DATA
    Traceback (most recent call last):
    ...
    ValueError: Unexpected Segments were passed. The series hmi.sharp_720s contains the following Segments ['magnetogram', 'bitmap', 'Dopplergram', 'continuum', 'inclination', 'azimuth', 'field', 'vlos_mag', 'dop_width', 'eta_0', 'damping', 'src_continuum', 'src_grad', 'alpha_mag', 'chisq', 'conv_flag', 'info_map', 'confid_map', 'inclination_err', 'azimuth_err', 'field_err', 'vlos_err', 'alpha_err', 'field_inclination_err', 'field_az_err', 'inclin_azimuth_err', 'field_alpha_err', 'inclination_alpha_err', 'azimuth_alpha_err', 'disambig', 'conf_disambig']

To get files for more than 1 segment at the same time, chain ``a.jsoc.Segment()`` using ``AND`` operator::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.sharp_720s'),
    ...             a.jsoc.Segment('continuum') & a.jsoc.Segment('magnetogram'))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    61 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME WAVELNTH CAR_ROT
    ----------------------- -------- --------- -------- -------
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
                        ...      ...       ...      ...     ...
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    Length = 61 rows
    <BLANKLINE>
    <BLANKLINE>

Using Keywords
==============
In some cases, you might want to filter out files based on key metadata, also called keywords.

A list of supported keywords of a series, say ``hmi.sharp_720s`` can be obtained by::


    >>> import drms
    >>> c = drms.Client()  # doctest: +REMOTE_DATA
    >>> keywords = c.keys('hmi.sharp_720s')  # doctest: +REMOTE_DATA
    >>> print(keywords)  # doctest: +REMOTE_DATA
    ['cparms_sg000', 'magnetogram_bzero', 'magnetogram_bscale', 'cparms_sg001', 'bitmap_bzero', 'bitmap_bscale', 'cparms_sg002', 'Dopplergram_bzero', 'Dopplergram_bscale', 'cparms_sg003', 'continuum_bzero', 'continuum_bscale', 'cparms_sg004', 'inclination_bzero', 'inclination_bscale', 'cparms_sg005', 'azimuth_bzero', 'azimuth_bscale', 'cparms_sg006', 'field_bzero', 'field_bscale', 'cparms_sg007', ... 'ERRJHT', 'ERRVF']

Each keyword needs to be compared to a value, e.g., ``a.jsoc.Keyword("bitmap_bzero") == 0`` or ``a.jsoc.Keyword("bitmap_bzero") > 1``.

An of this example is::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.sharp_720s'),a.jsoc.Keyword('bitmap_bzero') == 0) # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    61 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME WAVELNTH CAR_ROT
    ----------------------- -------- --------- -------- -------
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
                        ...      ...       ...      ...     ...
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    Length = 61 rows
    <BLANKLINE>
    <BLANKLINE>


You can pass multiple keywords and they will be chained together inside the query::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Series('hmi.sharp_720s'),
    ...             a.jsoc.Keyword('bitmap_bzero') == 0, a.jsoc.Keyword('continuum_bscale') > 0) # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    61 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME WAVELNTH CAR_ROT
    ----------------------- -------- --------- -------- -------
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
                        ...      ...       ...      ...     ...
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    Length = 61 rows
    <BLANKLINE>
    <BLANKLINE>

If you provide a keyword without a comparison it will raise an error::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.sharp_720s'),
    ...             a.jsoc.Keyword('bitmap_bzero'))  # doctest: +REMOTE_DATA
    Traceback (most recent call last):
    ...
    ValueError: Keyword 'bitmap_bzero' needs to have a comparison to a value.

If you provide an incorrect keyword name it will also raise a error::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.sharp_720s'),
    ...             a.jsoc.Keyword('bac') == 0)  # doctest: +REMOTE_DATA
    Traceback (most recent call last):
    ...
    ValueError: Keyword: 'bac' is not supported by series: hmi.sharp_720s

Using Sample
============
In case you need to query for data, at some interval of time, say every 10 min, you can pass it using `~sunpy.net.attrs.Sample`.
In other words, if you need to query for ``hmi.v_45s`` series data between January 1st from 00:00 to 01:00, 2014, every 10 minutes, you can do::

    >>> import astropy.units as u
    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.v_45s'), a.Sample(10*u.min))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    7 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2014.01.01_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:10:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:20:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:30:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:39:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:49:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    <BLANKLINE>
    <BLANKLINE>

Note that the argument passed in ``a.Sample()`` must be an Astropy quantity, convertible into seconds.

Constructing complex queries
============================

Complex queries can be built using ``OR`` operators.
Let's look for 2 different series data at the same time::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.v_45s') | a.jsoc.Series('aia.lev1_euv_12s'))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    81 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2014.01.01_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:01:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:02:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:04:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:05:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:06:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:06:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:07:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
                        ...      ...        ...      ...     ...
    2014.01.01_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    2107 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
           T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
    -------------------- -------- -------- -------- -------
    2014-01-01T00:00:01Z  SDO/AIA    AIA_4       94    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_1      131    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_2      193    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_2      211    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_1      335    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_4       94    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_1      131    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_3      171    2145
                     ...      ...      ...      ...     ...
    2014-01-01T00:59:49Z  SDO/AIA    AIA_2      211    2145
    2014-01-01T00:59:49Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:59:49Z  SDO/AIA    AIA_1      335    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_4       94    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_1      131    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_2      193    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_2      211    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_1      335    2145
    Length = 2107 rows
    <BLANKLINE>
    <BLANKLINE>

The two series names are joined together by the operator ``|``.
This is the ``OR`` operator.  Think of the above query as setting a set of conditions which get passed to the JSOC.
Let's say you want all the ``hmi.v_45s`` data from two separate days::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00') |
    ...             a.Time('2014-01-02T00:00:00', '2014-01-02T01:00:00'),
    ...             a.jsoc.Series('hmi.v_45s'))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    81 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2014.01.01_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:01:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:02:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:03:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:04:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:05:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:06:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:06:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:07:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
                        ...      ...        ...      ...     ...
    2014.01.01_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    81 Results from the JSOCClient:
    Source: http://jsoc.stanford.edu
    <BLANKLINE>
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2014.01.02_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:01:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:02:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:03:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:03:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:04:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:05:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:06:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:06:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:07:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
                        ...      ...        ...      ...     ...
    2014.01.02_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    <BLANKLINE>

Each of the arguments in this query style can be thought of as setting conditions that the returned records must satisfy.

It should be noted that ``AND`` operator is supported by some of the attributes only.
The attributes which support "&" are `~sunpy.net.jsoc.attrs.PrimeKey` and `~sunpy.net.jsoc.attrs.Segment`.
Using "&" with any other attributes will throw an error.

Downloading data
****************

To download the files located by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`,
you can download them by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.v_45s') | a.jsoc.Series('aia.lev1_euv_12s'),
    ...             a.jsoc.Notify('solar@example.com')  # doctest: +SKIP
    >>> downloaded_files = Fido.fetch(res)  # doctest: +SKIP

To export a request for download, you must have used the `sunpy.net.jsoc.attrs.Notify` attribute at search time to specify your email address.

.. note::

   **Only complete searches can be downloaded from JSOC**, this means that no slicing operations performed on the results object will affect the number of files downloaded.

Using JSOCClient for complex usage
**********************************

Fido interface uses `~sunpy.net.jsoc.JSOCClient` in its backend, and combines the last 2 stages the JSOC process into one.
You can directly use the JSOC client to make queries, instead of the Fido client.
This will allow you to separate the 3 stages of the JSOC process, and perform it individually, hence allowing a greater control over the whole process.

Setup
=====

sunpy's JSOC module is in `~sunpy.net`.  It can be imported as follows::

    >>> from sunpy.net import jsoc
    >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA

This creates your client object.

Making a query
==============

Querying JSOC using the JSOC client is very similar to what we were doing with Fido.
As above, we have to make sure we have an email address registered with JSOC before you are allowed to make a request.
See `this <http://jsoc.stanford.edu/ajax/register_email.html>`__ to register your email address.
We can add an email address to the search query with the `sunpy.net.jsoc.attrs.Notify` attribute.
Please note you can search without this but right now, you can not add the email address after the search::

    >>> from sunpy.net import attrs as a
    >>> res = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                     a.jsoc.Series('hmi.v_45s'),
    ...                     a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA

Apart from the function name, everything is the same. You need to pass the same values in the `~sunpy.net.jsoc.JSOCClient.search` as you did in `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`.
Complex queries can be built in a similar way, and all other things are the same.

Staging the request
===================

JSOC is a 3-stage process, and after getting the query results, we need to stage a request for the data to be downloaded.
Only then, can we download them. The download request can be staged like this::

    >>> requests = client.request_data(res)  # doctest: +SKIP
    >>> print(requests)  # doctest: +SKIP
    <ExportRequest id="JSOC_20170713_1461", status=0>

The function `~sunpy.net.jsoc.JSOCClient.request_data` stages the request.
It returns a `drms.client.ExportRequest` object, which has many attributes.
The most important ones are ``id`` and ``status``. Only when the status is 0, we can move to the third step, i.e. downloading the data.

If you are making more than 1 query at a time, it will return a list of `~drms.client.ExportRequest` objects.
Hence, access the list elements accordingly. You can get the id and status of the request (if it is not a list) by::

    >>> requests.id  # doctest: +SKIP
    JSOC_20170713_1461
    >>> requests.status  # doctest: +SKIP
    0

Downloading data
================

Once the status code is 0 you can download the data using the `~sunpy.net.jsoc.JSOCClient.get_request` method::

    >>> res = client.get_request(requests)  # doctest: +SKIP

This returns a Results instance which can be used to watch the progress of the download::

    >>> res.wait(progress=True)   # doctest: +SKIP
.. The tests in the module are very fragile to change in ordering of VSO
   results. In their current form they are very hard to maintain.

.. Skip all doctests here until the database locked error is fixed.
.. doctest-skip-all

.. _database_guide:

**************************
Using the database package
**************************
.. currentmodule:: sunpy.database

The database package offers the possibility to save retrieved data (e.g.
via the :mod:`sunpy.net` package) onto a local or remote database. The
database may be a single file located on a local hard drive (if a SQLite
database is used) or a local or remote database server (see the SQLAlchemy
documentation for a `list of supported databases
<https://docs.sqlalchemy.org/en/latest/dialects/index.html>`_)
This makes it possible to fetch required data from the local database
instead of downloading it again from a remote server.

The package :mod:`sunpy.database` was developed as part of Google Summer of Code (GSOC) 2013.

1. Connecting and initializing the database
*******************************************
To start a connection to an existing or a new database, create
a :class:`Database` object:

    >>> from sunpy.database import Database
    >>> database = Database('sqlite:///sunpydata.sqlite')

The database object in our example above connects to a new SQLite database with
the file name "sunpydata.sqlite" in the current directory.

The first parameter of :class:`Database` receives one mandatory argument:
a URL which describes how to connect to the database. The supported
format of this URL is described by the documentation of
:func:`sqlalchemy.create_engine`.

Note that a connection is only established when it's really needed,
i.e. if some query is sent to the database to read from it. Transactions
can also be committed explicitly using the :meth:`Database.commit` method.

.. warning::

    If you are using :class:`Database` objects in an interactive Python
    session you must not forget to call the :meth:`Database.commit` method
    on them explicitly before quitting the Python session! Otherwise, all
    changes on the altered databases are lost!

.. note::

    You can set the default database url in the sunpy config file under the
    'database' section. See :ref:`customizing-sunpy` for information on the
    config file. A database section might look like this::

        [database]
        url = sqlite:////home/user/sunpy/my_database.sqlite


2. Adding new entries
*********************
Each entry in a database is an instance of the class
:class:`tables.DatabaseEntry` with the following attributes:

====================== ===================================================
      Attribute                            Description
====================== ===================================================
id                     A unique ID number. By default it is None, but
                       automatically set to the maximum number plus one
                       when an entry is added to the database.
source                 The source is the name of an observatory or the
                       name of a network of observatories.
provider               The name of the server which provides the retrieved
                       data.
physobs                A physical observable identifier used by VSO.
fileid                 The file ID is a string defined by the data
                       provider that should point to a specific data
                       product. The association of fileid to the specific
                       data may change sometime, if the fileid always
                       points to the latest calibrated data.
observation_time_start The date and time when the observation of the data
                       started.
observation_time_end   The date and time when the observation of the data
                       ended.
instrument             The instrument which was used to observe the data.
size                   The size of the data in kilobytes (-1 if unknown).
wavemin                The value of the measured wave length.
wavemax                This is the same value as ``wavemin``. The value is
                       stored twice, because each
                       `sunpy.net.dataretriever.client.QueryResponse` which is
                       used by the vso package contains both these values.
path                   A local file path where the according FITS file is
                       saved.
download_time          The date and time when the files connected to a
                       query have been downloaded. Note: this is not the
                       date and time when this entry has been added to a
                       database!
starred                Entries can be starred to mark them. By default,
                       this value is False.
fits_header_entries    A list of :class:`tables.FitsHeaderEntry` instances.
fits_key_comments      A list of :class:`tables.FitsKeyComment` instances.
tags                   A list of :class:`tables.Tag` instances.
====================== ===================================================

* The ``id`` attribute is automatically set if an entry is added to a database.

* The attributes ``source``, ``provider``, ``physobs``, ``fileid``,
  ``observation_time_start``, ``observation_time_end``, ``instrument``,  ``size``,
  ``wavemin``, and ``wavemax`` are set by methods which use the VSO interface. In
  particular, these are :meth:`Database.add_from_vso_query_result`,
  and possibly :meth:`Database.fetch`.

* The attributes ``path`` and ``download_time`` are set by the method
  ``Database.download`` and also possibly by :meth:`Database.fetch`.

* ``starred`` is set or changed via :meth:`Database.star` or :meth:`Database.unstar`.

* ``tags`` is set via :meth:`Database.tag` or :meth:`Database.remove_tag`.

* The attribute ``fits_header_entries`` is set by the methods
  ``Database.download``, :meth:`Database.add_from_dir`, and
  :meth:`Database.add_from_file`.

2.1 Adding entries from one FITS file
=====================================
The method :meth:`Database.add_from_file` receives one positional argument
(either a string or a file-like object) which is used to add at least one
new entry from the given FITS file to the database. Why "at least one" and
not "exactly one"? The reason is that each database entry does not
represent one file, but one FITS header from a file. That means if you pass
a file which has 5 FITS headers in it, 5 entries will be added to the
database. The file in the following example
(``sunpy.data.sample.AIA_171_IMAGE``) has only one FITS header, so just one
entry is added to the database.

.. note::

  If you are working with Hinode/SOT files you may notice that for
  each file you get two entries, one which refers to the observation and
  another that contains some (useless) telemetry data.

:meth:`Database.add_from_file` saves the value of ``path`` by either simply
passing on the value of the received argument (if it was a string)
or by reading the value of ``file.name`` where ``file`` is the passed argument.
If the path cannot be determined, it stays as ``None`` (the default value).

``wavemin`` and ``wavemax`` are only set if the wavelength unit
of the FITS file can be found out or if the ``default_waveunit`` attribute of
the database object is set. These values are then
used to convert from the used units to nanometers. The rationale behind
this behaviour is that it makes querying for wavelengths more flexible.

``instrument`` is set by looking up the
FITS header key *INSTRUME*. ``observation_time_start`` is set
by searching for the FITS header key *DATE-OBS* or *DATE_OBS*.
Analogously, ``observation_time_end`` is set by searching for *DATE-END* or
*DATE_END*. Finally, the whole FITS header is stored in the attribute
``fits_header_entries`` as a list of :class:`tables.FitsHeaderEntry`
instances, and FITS comments are stored in the attribute
``fits_key_comments`` which is a list of :class:`tables.FitsKeyComment`
instances.

Using the function `len` on a :class:`Database` object returns the
number of saved database entries. To get the first entry of the database,
``database[0]`` is used (the ID number of the entries does not matter,
``database[0]`` always returns the oldest saved entry of the database). If
the database is empty, this expression raises an :exc:`IndexError`.
In section 3, more advanced formats of the slicing syntax are introduced.

    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> database.add_from_file(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    1
    >>> entry = database[0]  # doctest: +REMOTE_DATA
    >>> entry.path == sunpy.data.sample.AIA_171_IMAGE  # doctest: +REMOTE_DATA
    True
    >>> entry.wavemin, entry.wavemax  # doctest: +REMOTE_DATA
    (17.1, 17.1)
    >>> entry.instrument  # doctest: +REMOTE_DATA
    'AIA_3'
    >>> entry.observation_time_start, entry.observation_time_end  # doctest: +REMOTE_DATA
    (datetime.datetime(2011, 6, 7, 6, 33, 2, 770000), None)
    >>> len(entry.fits_header_entries)  # doctest: +REMOTE_DATA
    194
    >>> for fits_header_entry in entry.fits_header_entries[:10]:
    ...     print('{entry.key}\n    {entry.value}'.format(entry=fits_header_entry))  # doctest: +REMOTE_DATA
    SIMPLE
        True
    BITPIX
        -32
    NAXIS
        2
    NAXIS1
        1024
    NAXIS2
        1024
    PCOUNT
        0
    GCOUNT
        1
    XTENSION
        BINTABLE
    BLD_VERS
        V5R12X
    LVL_NUM
        1.5

    >>> for fits_key_comment in entry.fits_key_comments:
    ...     print('{comment.key}\n    {comment.value}'.format(comment=fits_key_comment))  # doctest: +REMOTE_DATA
    SIMPLE
        conforms to FITS standard
    BITPIX
        array data type
    NAXIS
        number of array dimensions
    NAXIS1
        width of table in bytes
    NAXIS2
        number of rows in table
    PCOUNT
        number of group parameters
    GCOUNT
        number of groups
    XTENSION
        binary table extension

2.2 Adding entries from a directory of FITS files
=================================================
Adding all FITS files from a directory works by calling the method
:meth:`Database.add_from_dir` and passing the desired directory to it. By
setting the keyword argument ``ignore_already_added`` to ``True``, no
exception is raised if it is attempted to add an already existing entry

In the following example case, setting this parameter is required because the
file ``sunpy.data.sample.AIA_171_IMAGE`` was already added, which is located
in the directory ``sampledata_dir``

    >>> from sunpy import config  # doctest: +REMOTE_DATA
    >>> sampledata_dir = config.get("downloads", "sample_dir")  # doctest: +REMOTE_DATA
    >>> database.default_waveunit = 'angstrom'  # doctest: +REMOTE_DATA
    >>> database.add_from_dir(sampledata_dir, ignore_already_added=True,
    ...                       time_string_parse_format="%d/%m/%Y")  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    59

2.3 Adding entries using the VSO interface
==========================================

2.3.1 Adding entries from a VSO query result
--------------------------------------------
The number of database entries that will be added from a VSO query result
is equal to the value of ``len(qr)`` in the following code snippet. Note that
:meth:`Database.add_from_vso_query_result` does not download any files,
though. If you want to add new entries using the VSO and also want to
download files at the same time, take a look at the following two
sections.

    >>> from sunpy.net import vso, attrs as a
    >>> client = vso.VSOClient()  # doctest: +REMOTE_DATA

.. testsetup:: *
   >>> from unittest.mock import Mock, patch
   >>> from sunpy.tests.mocks import MockObject
   >>> class MockQR:
   ...    def __init__(self, start_time, instrument, wavelength):
   ...        end_time = '{}{}'.format(start_time[:-1], int(start_time[-1]) + 1)
   ...        self.time = MockObject(start=start_time, end=end_time)
   ...        self.instrument = instrument
   ...        self.wave = MockObject(wavemin=wavelength, wavemax=wavelength, waveunit='nm')
   >>> with patch.object(vso.VSOClient, '__init__', lambda x: None):
   ...    vso.VSOClient().search = Mock(return_value=[MockQR(f"2011050800000{t}", 'AIA', w) for
   ...                                   t, w in zip([0,0,2,3], [17.1, 21.1, 9.4, 33.5])])

After initialising the VSO client:

    >>> qr = client.search(
    ...     a.Time('2011-05-08', '2011-05-08 00:00:05'),
    ...     a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> len(qr)  # doctest: +REMOTE_DATA
    4
    >>> database.add_from_vso_query_result(qr)  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    63

2.3.2 "Clever" Fetching
-----------------------

The method :meth:`Database.fetch` checks if the given query has already been
used once to add entries to the database. Otherwise, the query is used to
download and add new data. The :meth:`Database.fetch` method also accepts an
optional keyword argument ``path`` which is passed as-is to
:meth:`sunpy.net.vso.VSOClient.fetch` and determines the value of the ``path``
attribute of each entry.

Note that the number of entries that will be added depends on the total number
of FITS headers, **not** the number of records in the query.

.. testsetup:: *
   >>> qr1 = [MockQR(f"2012080500000{t}", 'AIA', w) for t,w in zip([1, 1, 2, 2], [9.4, 9.4, 33.5, 33.5])]
   >>> qr2 = [MockQR(f"2013080500000{t}", 'AIA', w) for t,w in zip([1, 1, 2, 2], [9.4, 9.4, 33.5, 33.5])]
   >>> produce_results = [qr1, None, qr2, None]
   >>> add_data = lambda x: database.add_from_vso_query_result(x) if x else None
   >>> database.fetch = Mock(side_effect=map(add_data, produce_results))

In the next code snippet new data is downloaded as this query
has not been downloaded before.

    >>> entries = database.fetch(
    ...     a.Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    67

Any more identical fetch calls don't add any new entries to the database
because they have already been downloaded.

    >>> entries = database.fetch(
    ...     a.Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    67

However, this different fetch call downloads new files because there is a
new date range whose files have not been downloaded yet.

    >>> entries = database.fetch(
    ...     a.Time('2013-08-05', '2013-08-05 00:00:05'),
    ...     a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    71

The caching also ensures that when queries have some results in
common, files for the common results will not be downloaded again. In the
following example, the query is new, but all of it's files have already been
downloaded. This means no new files are downloaded.

    >>> entries = database.fetch(
    ...     a.Time('2012-08-05 00:00:00', '2012-08-05 00:00:01'),
    ...     a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    71

2.4 Adding entries manually
===========================
Although usually not required, it is also possible to add database entries
by specifying the parameters manually. To do so pass the
values as keyword arguments to :class:`tables.DatabaseEntry` as follows:

    >>> from sunpy.database.tables import DatabaseEntry
    >>> entry = DatabaseEntry(instrument='EIT', wavemin=25.0)
    >>> database.add(entry)
    >>> entry in database
    False
    >>> database.commit()
    >>> entry in database
    True
    >>> len(database) # doctest: +REMOTE_DATA
    72

Note that the ``in`` operator only works as expected after the
:meth:`Database.commit` method has been called!

3. Displaying entries in a table
********************************
In the previous code snippets 71 entries have been added,
all of them saving a lot of data. To display the database in a table format
there is a helper function. :func:`tables.display_entries` takes two
arguments: the first one is an iterator of :class:`tables.DatabaseEntry`
instances. Remember that an instance of :class:`Database` yields instances
of :class:`tables.DatabaseEntry` instances, so you can simply pass a
database object. The second argument is an iterable of the resulting
columns in the table to be displayed. Each string in this iterable is used
to access the entry's attribute values. In the following example, the
values of ``entry.id``, ``entry.observation_time_start``,
``entry.observation_time_end``, ``entry.instrument``, ``entry.wavemin``,
and ``entry.wavemax`` are displayed (where ``entry`` stands for the
respective database entry). Note that *N/A* is displayed if the value
cannot be found or is not set.

    >>> from sunpy.database.tables import display_entries
    >>> print(display_entries(database,
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax']))   # doctest:  +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
      3    2011-06-07 06:33:02 ...               17.1               17.1
      4    2011-06-07 06:33:02 ...               21.1               21.1
      5    2011-06-07 06:33:03 ...               33.5               33.5
      6    2011-06-07 06:33:05 ...                9.4                9.4
      7    2011-06-07 06:33:05 ...              160.0              160.0
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
    ...                    ... ...                ...                ...
     62    2011-05-08 00:00:02 ...                9.4                9.4
     63    2011-05-08 00:00:03 ...               33.5               33.5
     64    2012-08-05 00:00:01 ...                9.4                9.4
     65    2012-08-05 00:00:01 ...                9.4                9.4
     66    2012-08-05 00:00:02 ...               33.5               33.5
     67    2012-08-05 00:00:02 ...               33.5               33.5
     68    2013-08-05 00:00:01 ...                9.4                9.4
     69    2013-08-05 00:00:01 ...                9.4                9.4
     70    2013-08-05 00:00:02 ...               33.5               33.5
     71    2013-08-05 00:00:02 ...               33.5               33.5
     72                    N/A ...               25.0                N/A
    Length = 72 rows

The index operator can be used to access certain single database entries like
``database[5]``  to get the 5th entry or  ``database[-2]``
to get the second-latest entry. It is also possible to
use more advanced slicing syntax;
see https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html
for more information.

    >>> print(display_entries(database[9::10],
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax']))   # doctest:  +REMOTE_DATA
     id observation_time_start observation_time_end instrument wavemin wavemax
    --- ---------------------- -------------------- ---------- ------- -------
     10    2011-06-07 06:39:31                  N/A      AIA_2    19.3    19.3
     20    2011-06-06 23:59:55  2011-06-08 00:00:05        GBM     N/A     N/A
     30                    N/A                  N/A        N/A     N/A     N/A
     40                    N/A                  N/A        N/A     N/A     N/A
     50                    N/A                  N/A        N/A     N/A     N/A
     60    2011-05-08 00:00:00  2011-05-08 00:00:01        AIA    17.1    17.1
     70    2013-08-05 00:00:02  2013-08-05 00:00:03        AIA    33.5    33.5

4. Removing entries
*******************
`Database.remove()` can be used to remove database entries from the sunpy
database. It takes a ``tables.DatabaseEntry`` object as argument.

For example, imagine we want to only have database entries which have an
observation time saved. To remove all the entries where the value of both
``observation_time_start`` and ``observation_time_end`` is ``None``,
simply iterate over the database and the :meth:`Database.remove`
method to remove those where there is no time set:

    >>> for database_entry in database:
    ...     if database_entry.observation_time_start is None and database_entry.observation_time_end is None:
    ...         database.remove(database_entry)
    ...
    >>> len(database) # doctest: +REMOTE_DATA
    39
    >>> print(display_entries(database,
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax']))   # doctest:  +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
      3    2011-06-07 06:33:02 ...               17.1               17.1
      4    2011-06-07 06:33:02 ...               21.1               21.1
      5    2011-06-07 06:33:03 ...               33.5               33.5
      6    2011-06-07 06:33:05 ...                9.4                9.4
      7    2011-06-07 06:33:05 ...              160.0              160.0
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
    ...                    ... ...                ...                ...
     61    2011-05-08 00:00:00 ...               21.1               21.1
     62    2011-05-08 00:00:02 ...                9.4                9.4
     63    2011-05-08 00:00:03 ...               33.5               33.5
     64    2012-08-05 00:00:01 ...                9.4                9.4
     65    2012-08-05 00:00:01 ...                9.4                9.4
     66    2012-08-05 00:00:02 ...               33.5               33.5
     67    2012-08-05 00:00:02 ...               33.5               33.5
     68    2013-08-05 00:00:01 ...                9.4                9.4
     69    2013-08-05 00:00:01 ...                9.4                9.4
     70    2013-08-05 00:00:02 ...               33.5               33.5
     71    2013-08-05 00:00:02 ...               33.5               33.5
    Length = 39 rows

5. Editing entries
******************

5.1 Starring and unstarring entries
===================================
The database package supports marking certain entries as "starred" using the
:meth:`Database.star` method. For example, to star
all values that have a wavelength of 20nm or higher:

    >>> for database_entry in database:
    ...     if database_entry.wavemin and database_entry.wavemin > 20:
    ...         database.star(database_entry)
    >>> print(display_entries(
    ...     filter(lambda entry: entry.starred, database),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))   # doctest:  +REMOTE_DATA
     id observation_time_start ...      wavemin           wavemax
    --- ---------------------- ... ----------------- -----------------
      4    2011-06-07 06:33:02 ...              21.1              21.1
      5    2011-06-07 06:33:03 ...              33.5              33.5
      7    2011-06-07 06:33:05 ...             160.0             160.0
     16    2011-06-07 06:32:11 ... 617.3000000000001 617.3000000000001
     61    2011-05-08 00:00:00 ...              21.1              21.1
     63    2011-05-08 00:00:03 ...              33.5              33.5
     66    2012-08-05 00:00:02 ...              33.5              33.5
     67    2012-08-05 00:00:02 ...              33.5              33.5
     70    2013-08-05 00:00:02 ...              33.5              33.5
     71    2013-08-05 00:00:02 ...              33.5              33.5

To remove the star from these entries, the :meth:`Database.unstar` method
works the same way.

5.2 Setting and removing tags
=============================
To add some more information by assigning tags to entries, the database package
also supports *tags*. The ``tags`` property of a database object holds all tags
that are saved in the database. For example, to assign the tag *spring* to all
entries that have been observed in March, April, or May on any day in any
year:

    >>> for database_entry in database:
    ...     if database_entry.observation_time_start.month in [3,4,5]:
    ...         database.tag(database_entry, 'spring')
    >>> database.tags  # doctest: +REMOTE_DATA
    [<Tag(name 'spring')>]
    >>> spring = database.tags[0]  # doctest: +REMOTE_DATA
    >>> print(display_entries(
    ...     filter(lambda entry: spring in entry.tags, database),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))  # doctest: +REMOTE_DATA
     id observation_time_start observation_time_end instrument wavemin wavemax
    --- ---------------------- -------------------- ---------- ------- -------
     17    2014-04-09 06:00:12                  N/A      AIA_3    17.1    17.1
     60    2011-05-08 00:00:00  2011-05-08 00:00:01        AIA    17.1    17.1
     61    2011-05-08 00:00:00  2011-05-08 00:00:01        AIA    21.1    21.1
     62    2011-05-08 00:00:02  2011-05-08 00:00:03        AIA     9.4     9.4
     63    2011-05-08 00:00:03  2011-05-08 00:00:04        AIA    33.5    33.5

5.3 Manually changing attributes
================================
Attributes for entries in the database can be manually edited. The
:meth:`Database.edit` method receives the database entry to be edited and any
number of keyword arguments to describe which values to change and how. For
example, to change the time entires that are ``None`` to the start of the
observation time (because it's possible, not because it is accurate!):

    >>> for database_entry in database:
    ...     if database_entry.observation_time_end is None:
    ...         database.edit(database_entry, observation_time_end=database_entry.observation_time_start)
    ...
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))   # doctest:  +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
      3    2011-06-07 06:33:02 ...               17.1               17.1
      4    2011-06-07 06:33:02 ...               21.1               21.1
      5    2011-06-07 06:33:03 ...               33.5               33.5
      6    2011-06-07 06:33:05 ...                9.4                9.4
      7    2011-06-07 06:33:05 ...              160.0              160.0
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
    ...                    ... ...                ...                ...
     61    2011-05-08 00:00:00 ...               21.1               21.1
     62    2011-05-08 00:00:02 ...                9.4                9.4
     63    2011-05-08 00:00:03 ...               33.5               33.5
     64    2012-08-05 00:00:01 ...                9.4                9.4
     65    2012-08-05 00:00:01 ...                9.4                9.4
     66    2012-08-05 00:00:02 ...               33.5               33.5
     67    2012-08-05 00:00:02 ...               33.5               33.5
     68    2013-08-05 00:00:01 ...                9.4                9.4
     69    2013-08-05 00:00:01 ...                9.4                9.4
     70    2013-08-05 00:00:02 ...               33.5               33.5
     71    2013-08-05 00:00:02 ...               33.5               33.5
    Length = 39 rows

It is possible to use
``database_entry.observation_time_end = database_entry.observation_time_start``",
but it has one major disadvantage: this operation cannot be undone.
See section 6 to see how undoing and redoing works.

6. Undoing and redoing operations
*********************************
A very handy feature of the database package is that every operation that
changes the database in some way can be reverted. The Database class has
the methods :meth:`Database.undo` and :meth:`Database.redo` to undo
and redo the last n commands, respectively.

.. warning::

    The undo and redo history are only saved in-memory! That means in
    particular, that if you work on a Database object in an interactive
    Python session and quit this session, the undo and redo history are
    lost.

In the following snippet, the operations from the sections 5.3, 5.2, and
5.1 are undone. Note the following changes: there are no longer any tags
anymore saved in the database, there is no entry which is starred and
the entries with no end of observation time are back.

    >>> database.undo(4)  # undo the edits from 5.3 (4 records have been affected)   # doctest:  +REMOTE_DATA
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred']))   # doctest:  +REMOTE_DATA
     id observation_time_start observation_time_end ...  tags  starred
    --- ---------------------- -------------------- ... ------ -------
      1    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      2    2011-06-07 06:33:01  2011-06-07 06:33:01 ...    N/A      No
      3    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      4    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A     Yes
      5    2011-06-07 06:33:03  2011-06-07 06:33:03 ...    N/A     Yes
      6    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A      No
      7    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A     Yes
      8    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
      9    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
     10    2011-06-07 06:39:31  2011-06-07 06:39:31 ...    N/A      No
    ...                    ...                  ... ...    ...     ...
     61    2011-05-08 00:00:00  2011-05-08 00:00:01 ... spring     Yes
     62    2011-05-08 00:00:02  2011-05-08 00:00:03 ... spring      No
     63    2011-05-08 00:00:03  2011-05-08 00:00:04 ... spring     Yes
     64    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     65    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     66    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     67    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     68    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     69    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     70    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
     71    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
    Length = 39 rows

The :meth:`~Database.redo` method reverts the last n operations that have been
undone. If not that many operations can be redone (i.e. any number greater than
14 in this example), an exception is raised. The redo call
reverts the original state: the tags have appeared again and all entries with a
wavelength >20nm are starred again, and there are no entries with no
stored end of observation time.

    >>> database.redo(1)  # redo all undone operations   # doctest:  +REMOTE_DATA
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred']))   # doctest:  +REMOTE_DATA
     id observation_time_start observation_time_end ...  tags  starred
    --- ---------------------- -------------------- ... ------ -------
      1    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      2    2011-06-07 06:33:01  2011-06-07 06:33:01 ...    N/A      No
      3    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      4    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A     Yes
      5    2011-06-07 06:33:03  2011-06-07 06:33:03 ...    N/A     Yes
      6    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A      No
      7    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A     Yes
      8    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
      9    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
     10    2011-06-07 06:39:31  2011-06-07 06:39:31 ...    N/A      No
    ...                    ...                  ... ...    ...     ...
     61    2011-05-08 00:00:00  2011-05-08 00:00:01 ... spring     Yes
     62    2011-05-08 00:00:02  2011-05-08 00:00:03 ... spring      No
     63    2011-05-08 00:00:03  2011-05-08 00:00:04 ... spring     Yes
     64    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     65    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     66    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     67    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     68    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     69    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     70    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
     71    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
    Length = 39 rows

7. Querying the database
************************
The API for querying databases is similar to querying the VSO using the
method :meth:`sunpy.net.vso.VSOClient.search`. The :meth:`Database.search`
method accepts any number of ORed query attributes (using \|) and
combines them using AND. It returns a list of matched database entries.
The special thing about querying databases is that all attributes support
the unary operator ``~`` to negate specific attributes. Example: the query
``Instrument.eit`` returns all entries that have *not* been observed
with the EIT.

7.1 Using VSO attributes
========================
Using the attributes from :mod:`sunpy.net.attrs` is quite intuitive:
the simple attributes and the Time attribute work exactly as you expect.

.. note::

  The ``near`` parameter of :class:`sunpy.net.attrs.Time` is ignored, because
  it's behaviour is not documented and it is different depending on the
  server which is requested.

The following query returns the data that was added in section 2.3.2:

    >>> print(display_entries(
    ...     database.search(a.Time('2012-08-05', '2012-08-05 00:00:05'), a.Instrument.aia),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'], sort=True))   # doctest:  +REMOTE_DATA
     id observation_time_start observation_time_end instrument wavemin wavemax
    --- ---------------------- -------------------- ---------- ------- -------
     64    2012-08-05 00:00:01  2012-08-05 00:00:02        AIA     9.4     9.4
     65    2012-08-05 00:00:01  2012-08-05 00:00:02        AIA     9.4     9.4
     66    2012-08-05 00:00:02  2012-08-05 00:00:03        AIA    33.5    33.5
     67    2012-08-05 00:00:02  2012-08-05 00:00:03        AIA    33.5    33.5

.. NOTE the following code does not actually work. There seems to be a bug
    in sunpy.util.unit_conversion.to_angstrom (this is called by the
    __init__ of the Wave attribute)

When using the :class:`sunpy.net.attrs.Wavelength` attribute, you have to
specify a unit using `astropy.units.Quantity`. If not an error is
raised. This also means that there is no default unit that is
used by the class. To know how you can specify a detail using astropy
check `astropy.units`.

    >>> import astropy.units as u
    >>> print(display_entries(
    ...     database.search(a.Wavelength(1.0*u.nm, 2.0*u.nm)),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'], sort=True))   # doctest:  +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
     10    2011-06-07 06:39:31 ...               19.3               19.3
     11    2011-06-07 06:45:55 ...               19.3               19.3
     12    2011-06-07 06:52:19 ...               19.3               19.3
     13    2011-06-07 06:58:43 ...               19.3               19.3
     17    2014-04-09 06:00:12 ...               17.1               17.1
     18    2011-06-07 20:37:52 ...               19.5               19.5
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
      3    2011-06-07 06:33:02 ...               17.1               17.1
     58    2011-06-07 06:33:29 ... 17.400000000000002 17.400000000000002
     60    2011-05-08 00:00:00 ...               17.1               17.1
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3

7.2 Database-specific attributes
================================
There are 5 additional query attributes supported by the database package.
They can be imported from the submodule :mod:`sunpy.database.attrs` and are in
particular:

- Starred
- Tag
- Path
- DownloadTime
- FitsHeaderEntry

The following query searches for all entries that have the tag 'spring' or
(inclusive or!) are starred and have not the FITS header key 'WAVEUNIT'
with the value 'Angstrom':

    >>> import sunpy.database.attrs as dbattrs
    >>> print(display_entries(
    ...     database.search(dbattrs.Tag('spring') | dbattrs.Starred(), ~dbattrs.FitsHeaderEntry('WAVEUNIT', 'Angstrom')),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred'], sort=True))   # doctest:  +REMOTE_DATA
     id observation_time_start observation_time_end ...  tags  starred
    --- ---------------------- -------------------- ... ------ -------
     16    2011-06-07 06:32:11  2011-06-07 06:32:11 ...    N/A     Yes
     17    2014-04-09 06:00:12  2014-04-09 06:00:12 ... spring      No
      4    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A     Yes
      5    2011-06-07 06:33:03  2011-06-07 06:33:03 ...    N/A     Yes
     60    2011-05-08 00:00:00  2011-05-08 00:00:01 ... spring      No
     61    2011-05-08 00:00:00  2011-05-08 00:00:01 ... spring     Yes
     62    2011-05-08 00:00:02  2011-05-08 00:00:03 ... spring      No
     63    2011-05-08 00:00:03  2011-05-08 00:00:04 ... spring     Yes
     66    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     67    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
      7    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A     Yes
     70    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
     71    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes


8. Caching
**********
All entries that are saved in the database are also saved in a cache
in-memory. The type of the cache is determined at the initialization of
the database object and cannot be changed after that. The default type is
:class:`caching.LRUCache` (least-recently used) and the other one which is
supported is :class:`caching.LFUCache` (least-frequently used). By
default the cache size is ``float('inf')``, i.e. infinite. To set the
cache size after the database object has been initialized, use the method
:meth:`Database.set_cache_size`. If the new size is smaller than the
current number of database entries, entries are removed according to the
cache type until the number of entries is equal to the given cache size.
The following call to :meth:`Database.set_cache_size` sets the cache size
to 10 and therefore removes the 5 entries that been used least recently.

    >>> database.set_cache_size(10)
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))   # doctest:  +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
     11    2011-06-07 06:45:55 ...               19.3               19.3
     12    2011-06-07 06:52:19 ...               19.3               19.3
     13    2011-06-07 06:58:43 ...               19.3               19.3
     16    2011-06-07 06:32:11 ...  617.3000000000001  617.3000000000001
     17    2014-04-09 06:00:12 ...               17.1               17.1
     18    2011-06-07 20:37:52 ...               19.5               19.5
     58    2011-06-07 06:33:29 ... 17.400000000000002 17.400000000000002
     59    2011-06-06 00:00:00 ...                N/A                N/A
************************
Using sunpy's HEK module
************************

The Heliophysics Event Knowledgebase (HEK) is a repository of feature
and event information about the Sun.
Entries are generated both by automated algorithms and human observers.
sunpy accesses this information through the `~sunpy.net.hek` module, which was developed through support from the European Space Agency Summer of Code in Space (ESA-SOCIS) 2011.

A simple query
**************

To search the HEK, you need a start time, an end time, and an event type.
Times are specified as strings or Python datetime objects.
Event types are specified as upper case, two letter strings, and are identical to the two letter abbreviations found at the HEK website, http://www.lmsal.com/hek/VOEvent_Spec.html.

    >>> from sunpy.net import attrs as a
    >>> from sunpy.net import Fido
    >>> tstart = '2011/08/09 07:23:56'
    >>> tend = '2011/08/09 12:40:29'
    >>> event_type = 'FL'
    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type))  # doctest: +REMOTE_DATA

``tstart`` and ``tend`` defines the start and end times of the query, and ``event_type`` specifies the event type which in this example we are searching for flares defined as ``FL``.
Line 6 goes out to the web, contacts the HEK, and queries it for the information you have requested.
Event data for ALL flares available in the HEK within the time range 2011/08/09 07:23:56 UT - 2011/08/09 12:40:20 UT will be returned, regardless of which feature recognition method used to detect the flare.

Let's break down the arguments of ``Fido.search``.
The first argument::

    a.Time(tstart,tend)

sets the start and end times for the query.

The second argument::

    a.hek.EventType(event_type)

sets the type of event to look for.
Since we have defined ``event_type = 'FL'``, this sets the query to look for flares.
We could have also set the flare event type using the syntax::

    a.hek.FL

There is more on the attributes of hek in section 3 of this guide.

The result
**********

So, how many flare detections did the query turn up?

The result object returned by ``Fido.search`` is a `~.UnifiedResponse` object which contains all the results from any clients used in the search.
The first thing we need to do is access the results from the HEK client, the only ones for the query we gave.

    >>> len(result['hek'])  # doctest: +REMOTE_DATA
    19

This object is an `astropy.table.Table` object with the columns which correspond to the parameters listed at http://www.lmsal.com/hek/VOEvent_Spec.html.

You can inspect all results very simply:

    >>> result['hek']  # doctest: +SKIP

Remember, the HEK query we made returns all the flares in the time-range stored in the HEK, regardless of the feature recognition method.
The HEK parameter which stores the the feature recognition method is called "frm_name".
We can select just this column:

    >>> result["hek"]["frm_name"]  # doctest: +REMOTE_DATA
    <QueryResponseColumn name='frm_name' dtype='str32' length=19>
                              asainz
                              asainz
                              asainz
                              asainz
                              asainz
                              asainz
                              asainz
                   SSW Latest Events
                                SWPC
    Flare Detective - Trigger Module
    Flare Detective - Trigger Module
                                SWPC
                   SSW Latest Events
    Flare Detective - Trigger Module
    Flare Detective - Trigger Module
    Flare Detective - Trigger Module
    Flare Detective - Trigger Module
    Flare Detective - Trigger Module
    Flare Detective - Trigger Module

It is likely each flare on the Sun was actually detected multiple times by many different methods.

More complex queries
********************

The Fido allows you to make more complex queries.
There are two key features you need to know in order to make use of the full power of Fido.
Firstly, the attribute module - ``attrs.hek`` - describes ALL the parameters stored by the HEK as listed in http://www.lmsal.com/hek/VOEvent_Spec.html, and the HEK client makes these parameters searchable.

To explain this, let's have a closer look at ``attrs.hek``.
The help command is your friend here; scroll down to section DATA you will see:

    >>> help(a.hek) # doctest:+REMOTE_DATA
    Help on module sunpy.net.hek.attrs in sunpy.net.hek:
    <BLANKLINE>
    NAME
        sunpy.net.hek.attrs
    <BLANKLINE>
    DESCRIPTION
        Attributes that can be used to construct HEK queries. They are different to
        the VSO ones in that a lot of them are wrappers that conveniently expose
        the comparisons by overloading Python operators. So, e.g., you are able
        to say AR & AR.NumSpots < 5 to find all active regions with less than 5 spots.
        As with the VSO query, you can use the fundamental logic operators AND and OR
        to construct queries of almost arbitrary complexity. Note that complex queries
        result in multiple requests to the server which might make them less efficient.
    <BLANKLINE>
    CLASSES
    ...

You'll see that one of the attributes is a flare object::

    FL = <sunpy.net.hek.attrs.FL object>

We can replace a.hek.EventType('FL') with a.hek.FL - they do the same thing, setting the query to look for flare events.
Both methods of setting the event type are provided as a convenience

Let's look further at the FRM attribute::

    >>> help(a.hek.FRM) # doctest:+REMOTE_DATA
    Help on FRM in module sunpy.net.hek.attrs object:
    <BLANKLINE>
    class FRM(builtins.object)
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  Contact = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  HumanFlag = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  Identifier = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  Institute = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  Name = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  ParamSet = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  SpecificID = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  URL = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  VersionNumber = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
    <BLANKLINE>

Let's say I am only interested in those flares identified by the SSW Latest Events tool.
I can retrieve those entries only from the HEK with the following command:

    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.FRM.Name == 'SSW Latest Events')  # doctest: +REMOTE_DATA
    >>> len(result[0])  # doctest: +REMOTE_DATA
    2

We can also retrieve all the entries in the time range which were not made by SSW Latest Events with the following command:

    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.FRM.Name != 'SSW Latest Events')  # doctest: +REMOTE_DATA
    >>> len(result[0])  # doctest: +REMOTE_DATA
    19

We are using Python's comparison operators to filter the returns from Fido.
Other comparisons are possible.
For example, let's say I want all the flares that have a peak flux of over 4000.0:

    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.FL.PeakFlux > 4000.0)  # doctest: +REMOTE_DATA
    >>> len(result[0])  # doctest: +REMOTE_DATA
    1

Multiple comparisons can be included.
For example, let's say I want all the flares with a peak flux above 1000 AND west of 800 arcseconds from disk center of the Sun:

    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.Event.Coord1 > 800, a.hek.FL.PeakFlux > 1000.0)  # doctest: +REMOTE_DATA

Multiple comparison operators can be used to filter the results back from the HEK.

The second important feature about the HEK client is that the comparisons we've made above can be combined using Python's logical operators.
This makes complex queries easy to create.
However, some caution is advisable.
Let's say I want all the flares west of 50 arcseconds OR have a peak flux over 1000.0:

    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), (a.hek.Event.Coord1 > 50) or (a.hek.FL.PeakFlux > 1000.0))  # doctest: +REMOTE_DATA

and as a check:

    >>> result["hek"]["fl_peakflux"] # doctest: +REMOTE_DATA
    <QueryResponseColumn name='fl_peakflux' dtype='object' length=17>
       None
       None
       None
       None
       None
       None
       None
    2326.86
    1698.83
       None
       None
    2360.49
    3242.64
    1375.93
    6275.98
    923.984
    1019.83

    >>> result["hek"]["event_coord1"] # doctest: +REMOTE_DATA
    <QueryResponseColumn name='event_coord1' dtype='float64' length=17>
     51.0
     51.0
     51.0
    924.0
    924.0
    924.0
     69.0
    883.2
    883.2
     69.0
     69.0
    883.2
    883.2
    883.2
    883.2
    883.2
    883.2

Note that some of the fluxes are returned as "None".
This is because some feature recognition methods for flares do not report the peak flux.
However, because the location of ``event_coord1`` is greater than 50, the entry from the HEK for that flare detection is returned.

Let's say we want all the flares west of 50 arcseconds AND have a peak flux over 1000.0:

    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), (a.hek.Event.Coord1 > 50) and (a.hek.FL.PeakFlux > 1000.0))  # doctest: +REMOTE_DATA

    >>> result["hek"]["fl_peakflux"] # doctest: +REMOTE_DATA
    <QueryResponseColumn name='fl_peakflux' dtype='float64' length=7>
    2326.86
    1698.83
    2360.49
    3242.64
    1375.93
    6275.98
    1019.83
    >>> result["hek"]["event_coord1"] # doctest: +REMOTE_DATA
    <QueryResponseColumn name='event_coord1' dtype='float64' length=7>
    883.2
    883.2
    883.2
    883.2
    883.2
    883.2
    883.2

In this case none of the peak fluxes are returned with the value `None`.
Since we are using an ```and`` logical operator we need a result from the ``(a.hek.FL.PeakFlux > 1000.0)`` filter.
Flares that have `None` for a peak flux cannot provide this, and so are excluded.
The `None` type in this context effectively means "Don't know"; in such cases the client returns only those results from the HEK that definitely satisfy the criteria passed to it.

Getting data for your event
***************************

The 'hek2vso' module allows you to take an HEK event and acquire VSO records specific to that event and was developed with support from the 2013 Google Summer of Code.

    >>> from sunpy.net import hek2vso
    >>> h2v = hek2vso.H2VClient()  # doctest: +REMOTE_DATA

There are several ways to use this capability.
For example, you can pass in a list of HEK results and get out the corresponding VSO records.
Here are the VSO records returned via the tenth result from the HEK query in Section 2 above:

    >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type))  # doctest: +REMOTE_DATA
    >>> vso_records = h2v.translate_and_query(result[0][10])  # doctest: +REMOTE_DATA
    >>> len(vso_records[0])  # doctest: +REMOTE_DATA
    31

``result[0][10]`` is the HEK entry generated by the "Flare Detective" automated flare detection algorithm running on the AIA 193 angstrom waveband.
The VSO records are for full disk AIA 193 angstrom images between the start and end times of this event.
The 'translate_and_query' function uses exactly that information supplied by the HEK in order to find the relevant data for that event.
Note that the VSO does not generate records for all solar data, so it is possible that an HEK entry corresponds to data that is not accessible via the VSO.

You can also go one step further back, passing in a list of HEK attribute objects to define your search, the results of which are then used to generate their corresponding VSO records:

   >>> q = h2v.full_query((a.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), a.hek.EventType('FL')))  # doctest: +SKIP

The full capabilities of the HEK query module can be used in this function (see above).

Finally, for greater flexibility, it is possible to pass in a list of HEK results and create the corresponding VSO query attributes.

    >>> vso_query = hek2vso.translate_results_to_query(result[0][10])  # doctest: +REMOTE_DATA
    >>> vso_query[0]  # doctest: +REMOTE_DATA
    [<sunpy.net.attrs.Time(2011-08-09 07:22:44.000, 2011-08-09 07:28:56.000)>, <sunpy.net.attrs.Source(SDO: The Solar Dynamics Observatory.) object at ...>, <sunpy.net.attrs.Instrument(AIA: Atmospheric Imaging Assembly) object at ...>, <sunpy.net.attrs.Wavelength(193.0, 193.0, 'Angstrom')>]

This function allows users finer-grained control of VSO queries generated from HEK results.
***********************************
Querying Helioviewer.org with sunpy
***********************************

sunpy can be used to make several basic requests using the The `Helioviewer.org API <https://api.helioviewer.org/docs/v2/>`_ including generating a PNG screenshot and downloading a `JPEG 2000 <https://wiki.helioviewer.org/wiki/JPEG_2000>`_ image.

As you can get JPEG 2000 images, you will need two other pieces of software in order to open them in Python.
The first is OpenJPEG which is an open source library for reading and writing JPEG2000 files.
The other package you will need is `Glymur <https://pypi.python.org/pypi/Glymur/>`_.
Both of these are available as `conda <https://www.anaconda.com/>`_ packages and ideally should be installed in this manner.
Otherwise, please follow the instructions at `the OpenJPEG homepage <http://www.openjpeg.org>`_ and the `Glymur homepage <https://glymur.readthedocs.io/en/latest/>`_.

To interact with the Helioviewer API, users first create a `~sunpy.net.helioviewer.HelioviewerClient` instance.
The client instance can then be used to make various queries against the API using the same parameters one would use when making a web request.
Note that the HelioviewerClient does not currently offer full access to the HelioViewer API.

We provide the follwing functions:

1. Download a JPEG 2000 image for the specified datasource that is the closest match in time to the ``date`` requested.
2. As above but return the full JSON response instead of an image.
3. Return all available datasources.
4. Generates custom screenshots that allow labels and layers of images.
5. Download the header information present in a JPEG2000 image. This includes:
    - FITS header
    - Helioviewer-specific metadata.

Nearly all requests require the user to specify the data they are interested in.
Depending on the function, it consists of passing in either: *observatory*, *instrument*,
*detector* and *measurement* or *source_id* keywords.

To find out what the allowed values are, you can access the `~sunpy.net.helioviewer.HelioviewerClient.data_sources` dictionary that is an attribute of the Helioviewer client.
Let us begin by retrieving the available list of sources that Helioviewer supports by using `~sunpy.net.helioviewer.HelioviewerClient.data_sources`::

    >>> from sunpy.net import helioviewer
    >>> hv = helioviewer.HelioviewerClient()  # doctest: +REMOTE_DATA
    >>> for sourceid, obs in hv.data_sources.items():# doctest: +REMOTE_DATA
    ...     print(f"{sourceid}: {obs}")  # doctest: +REMOTE_DATA
    ('SOHO', 'EIT', None, '171'): 0
    ('SOHO', 'EIT', None, '195'): 1
    ('SOHO', 'EIT', None, '284'): 2
    ('SOHO', 'EIT', None, '304'): 3
    ('SOHO', 'LASCO', 'C2', 'white-light'): 4
    ('SOHO', 'LASCO', 'C3', 'white-light'): 5
    ('SOHO', 'MDI', None, 'magnetogram'): 6
    ('SOHO', 'MDI', None, 'continuum'): 7
    ('SDO', 'AIA', None, '94'): 8
    ('SDO', 'AIA', None, '131'): 9
    ('SDO', 'AIA', None, '171'): 10
    ('SDO', 'AIA', None, '193'): 11
    ('SDO', 'AIA', None, '211'): 12
    ('SDO', 'AIA', None, '304'): 13
    ('SDO', 'AIA', None, '335'): 14
    ('SDO', 'AIA', None, '1600'): 15
    ('SDO', 'AIA', None, '1700'): 16
    ('SDO', 'AIA', None, '4500'): 17
    ('SDO', 'HMI', None, 'continuum'): 18
    ('SDO', 'HMI', None, 'magnetogram'): 19
    ...

Every JPEG 2000 file provided by the Helioviewer Project has been processed to generate an image that
can be used for browsing purposes.
This typically involves following the standard image processing procedure used by each instrument team to convert their science data into an image for a webpage.
The JPEG 2000 image is then scaled between 0 and 255 (byte-scaled).
**Please note that the JPEG 2000 image data is not the same as the original science data.**

Suppose we want to download a JPEG 2000 image of the latest AIA 304 image available on Helioviewer.org.
From the list above, we know that SDO/AIA 304  is ``(('SDO', 'AIA', None, '304'), 13)``.
So ``observatory="SDO"``,``instrument=AIA``, ``detector=None``, ``measurement=304`` and the ``source_id`` is 13.
So we can use the approach as shown in the following example::

   >>> from sunpy.net.helioviewer import HelioviewerClient
   >>> import matplotlib.pyplot as plt
   >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
   >>> file = hv.download_jp2('2012/01/01', observatory="SDO", instrument="AIA",
   ...                        measurement="304")  # doctest: +REMOTE_DATA

Since ``detector=None`` we can ignore this keyword and skip it when we call this function.
As we also have the source_id for AIA 304, which is ``13``, we could make the same request using: ::

   file = hv.download_jp2('2012/01/01', source_id=13)

Since this is a JPEG 2000 image, to plot this image you can either call Glymur directly::

   >>> import glymur # doctest: +SKIP
   >>> im = glymur.Jp2k(file)[:]  # doctest: +SKIP

The better method is to load the image into a sunpy Map object::

   >>> from sunpy.map import Map
   >>> aia = Map(file)  # doctest: +SKIP
   >>> aia.peek()  # doctest: +SKIP

.. image:: helioviewer-1.png

The sunpy Map selects a color table based on the JPEG 2000 image meta data for plotting.
This will be the color table that is used by the Helioviewer Project to display JPEG 2000 images in their own clients.

In this example we will query Helioviewer for the relevant JPEG 2000 file closest to the input time, for a SDO/HMI continuum image and crop to focus on an active region::

   >>> from sunpy.net.helioviewer import HelioviewerClient
   >>> import matplotlib.pyplot as plt
   >>> from astropy.units import Quantity
   >>> from sunpy.map import Map
   >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
   >>> data_sources = hv.get_data_sources()  # doctest: +REMOTE_DATA
   >>> filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO',
   ...                            instrument='HMI', measurement='continuum')  # doctest: +REMOTE_DATA
   >>> hmi = Map(filepath)  # doctest: +SKIP
   >>> xrange = Quantity([200, 550], 'arcsec')  # doctest: +REMOTE_DATA
   >>> yrange = Quantity([-400, 200], 'arcsec')  # doctest: +REMOTE_DATA
   >>> hmi.submap(xrange, yrange).peek()  # doctest: +SKIP

.. image:: helioviewer-2.png

The other main method is `~sunpy.net.helioviewer.HelioviewerClient.download_png`.
This allows more complex images to be created but again these are not the original science data.
The biggest difference is that we do not use the separate keywords but have to pass them as a string of lists.
This is the ``layer`` keyword in this function.

We will recreate the first example using the PNG function::

   >>> from sunpy.net.helioviewer import HelioviewerClient
   >>> import matplotlib.pyplot as plt
   >>> from matplotlib.image import imread
   >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
   >>> file = hv.download_png('2020/01/01', 4.8, "[SDO,AIA,304,1,100]", x0=0, y0=0, width=768, height=768, watermark=True)  # doctest: +REMOTE_DATA
   >>> im = imread(file)  # doctest: +REMOTE_DATA
   >>> plt.imshow(im)  # doctest: +SKIP
   >>> plt.axis('off')  # doctest: +SKIP
   >>> plt.show()  # doctest: +SKIP

.. image:: helioviewer-3.png

Since this is just a PNG, we can use matplotlib directly to plot this image.
Note that the filename of the returned file has the date and time of the request, not of any of the times shown in the image itself.
**This is not a bug.**
The reason for this is that the user may ask for images from multiple sources, and each of them may have a different observation time.
The problem becomes which time is the most appropriate to associate with the resultant image.
Helioviewer.org doesn't choose between the images times, but instead uses the request time to construct the image filename.
This means that the image file names for request times in the future (like in this example) can look a little unusual compared to the times in the image.

After the date string, we have a number (``4.8``) which refers to the image resolution in arcseconds per pixel (larger values mean lower resolution).
The next input is the ``layers`` keyword which is ``"[SDO,AIA,304,1,100]"``.
The first 4 are the observatory, instrument, detector, measurement values from before.
Note that since SDO AIA has no detector value, you can skip this within the list.
The ``1`` and ``100`` in the layer list refer to the visibility and opacity of the datasource.
You can use the ``sourceid`` instead of the keywords, so it would be ``[13,1,100]`` for this example.

Finally, the ``x0`` and ``y0`` are the center points about which to focus and the ``width`` and ``height`` are the pixel values for the image dimensions.
These have defaults set so you do not need to supply these.

In this example we will create a composite PNG image using data from two different SDO AIA wavelengths and LASCO C2 coronagraph data.
The layer string is extended to include the additional data sources, and opacity is throttled down for the second AIA layer so that it does not completely block out the lower layer::

   >>> from sunpy.net.helioviewer import HelioviewerClient
   >>> import matplotlib.pyplot as plt
   >>> from matplotlib.image import imread
   >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
   >>> file = hv.download_png('2012/01/01', 6,
   ...                        "[SDO,AIA,304,1,100],[SDO,AIA,193,1,50],[SOHO,LASCO,C2,white-light,1,100]",
   ...                        x0=0, y0=0, width=768, height=768, watermark=True)  # doctest: +REMOTE_DATA
   >>> im = imread(file)  # doctest: +REMOTE_DATA
   >>> plt.imshow(im)  # doctest: +SKIP
   >>> plt.axis('off')  # doctest: +SKIP
   >>> plt.show()  # doctest: +SKIP

.. image:: helioviewer-4.png

For more information about using querying Helioviewer.org, see the `Helioviewer.org
API documentation <https://api.helioviewer.org/docs/v2/>`_.
.. _fido_guide:

***************************************
Finding and Downloading Data using Fido
***************************************

This guide outlines how to search for and download data using `~sunpy.net.Fido` sunpy's interface for search and download.
`~sunpy.net.Fido` is a unified interface for searching and fetching solar physics data irrespective of the underlying client or webservice through which the data is obtained, e.g. VSO_,JSOC_, etc.
It therefore supplies a single, easy and consistent way to obtain most forms of solar physics data.

Import
******

The `~sunpy.net.Fido` object is in `sunpy.net`.
It can be imported as follows::

    >>> from sunpy.net import Fido, attrs as a

Search Attributes
*****************

To search for data with `~sunpy.net.Fido`, you need to specify attributes to search against.
The range of attributes are found in the `attrs <sunpy.net.attrs>` submodule.
Examples of these attributes are:

- `a.Time <sunpy.net.attrs.Time>`
- `a.Instrument <sunpy.net.attrs.Instrument>`
- `a.Wavelength <sunpy.net.attrs.Wavelength>`

whereas some of these attributes are client specific, and are found under client specific submodules, e.g. `attrs.vso <sunpy.net.vso.attrs>` and `attrs.jsoc <sunpy.net.jsoc.attrs>`.

In to each attribute you have to provide a value to use::

    >>> a.Time('2012/3/4', '2012/3/6'), a.Instrument.lyra
    (<sunpy.net.attrs.Time(2012-03-04 00:00:00.000, 2012-03-06 00:00:00.000)>, <sunpy.net.attrs.Instrument(LYRA: Lyman Alpha Radiometer is the solar UV radiometer on board
    Proba-2.) object at ...>)

For attributes that have no fixed selection of values (``Time`` for example) you will have to provide the range you require.
However, for attributes that have a fixed range of **known** values, it is possible to list all these values.
**Please note that each list is not exhaustive.**

Using ``Instrument`` as the first example, if you print the object::

    >>> print(a.Instrument)
    sunpy.net.attrs.Instrument
    <BLANKLINE>
    Specifies the Instrument name for the search.
    <BLANKLINE>
           Attribute Name          Client          Full Name                                           Description
    --------------------------- ----------- ------------------------ --------------------------------------------------------------------------------
    aia                         VSO         AIA                      Atmospheric Imaging Assembly
    bcs                         VSO         BCS                      Bragg Crystal Spectrometer
    be_continuum                VSO         BE-Continuum             INAF-OACT Barra Equatoriale Continuum Instrument
    be_halpha                   VSO         BE-Halpha                INAF-OACT Barra Equatoriale Hα Instrument
    bigbear                     VSO         Big Bear                 Big Bear Solar Observatory, California TON and GONG+ sites
    caii                        VSO         CAII                     Kanzelhöhe Ca II k Instrument
    cds                         VSO         CDS                      Coronal Diagnostic Spectrometer
    celias                      VSO         CELIAS                   Charge, Element, and Isotope Analysis System
    ...

You get a full list of known values, a description and what "Clients" support those values (if you want to use a specific data source).
This is supported for most attributes including the client specific ones.


For example you can print a list of Series provided by JSOC::

    >>> print(a.jsoc.Series)
    sunpy.net.jsoc.attrs.Series
    <BLANKLINE>
    The JSOC Series to Download.
    <BLANKLINE>
              Attribute Name           Client             Full Name                                                Description
    ---------------------------------- ------ ---------------------------------- --------------------------------------------------------------------------------
    aia_flatfield                      JSOC   aia.flatfield                      AIA flatfield
    aia_lev1                           JSOC   aia.lev1                           AIA Level 1
    aia_lev1_euv_12s                   JSOC   aia.lev1_euv_12s                   AIA Level 1, 12 second cadence
    aia_lev1_uv_24s                    JSOC   aia.lev1_uv_24s                    AIA Level 1, 24 second cadence
    aia_lev1_vis_1h                    JSOC   aia.lev1_vis_1h                    AIA Level 1, 3600 second cadence
    aia_master_pointing3h              JSOC   aia.master_pointing3h              Master Pointing Parameters
    aia_response                       JSOC   aia.response                       AIA instrument response table
    aia_temperature_summary_300s       JSOC   aia.temperature_summary_300s       Temperature Statistics from AIA Housekeeping - Thermal Packet
    hmi_b_135s                         JSOC   hmi.b_135s                         Full-disk Milne-Eddington inversion with the azimuth disambiguation informati...
    ...

Furthermore, you can use tab completion to auto-fill the attribute name, for example by typing ``a.jsoc.aia_f<TAB>``.

Searching for Data Using Fido
*****************************

For example::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.lyra, a.Level.two) # doctest: +REMOTE_DATA

this returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing information on the results which fit the criteria specified by the attrs objects in the above call.
It does not download the files.
For instructions on how to download data using Fido, see :ref:`downloading_data`.

To see a summary of the results of our query, simply type the name of the variable set to the Fido search, in this case, result::

    >>> result  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA ...      ESA     2
    2012-03-05 00:00:00.000 2012-03-05 23:59:59.999       LYRA ...      ESA     2
    2012-03-06 00:00:00.000 2012-03-06 23:59:59.999       LYRA ...      ESA     2
    <BLANKLINE>
    <BLANKLINE>

Queries can be made more flexible or specific by adding more attrs objects to the `~sunpy.net.Fido` search.
Specific passbands can be searched for by supplying an `~astropy.units.Quantity` to the `a.Wavelength <sunpy.net.attrs.Wavelength>` attribute::

    >>> import astropy.units as u
    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.norh,
    ...             a.Wavelength(17*u.GHz))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the NoRHClient:
    Source: https://solar.nro.nao.ac.jp/norh/doc/manuale/node1.html
    <BLANKLINE>
           Start Time               End Time        ... Provider Wavelength
                                                    ...             GHz
    ----------------------- ----------------------- ... -------- ----------
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999 ...      NRO       17.0
    2012-03-05 00:00:00.000 2012-03-05 23:59:59.999 ...      NRO       17.0
    2012-03-06 00:00:00.000 2012-03-06 23:59:59.999 ...      NRO       17.0
    <BLANKLINE>
    <BLANKLINE>

Data of a given cadence can also be specified using the Sample attribute.
To search for data at a given cadence use the `a.Sample <sunpy.net.attrs.Sample>` attribute.

    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.aia,
    ...             a.Wavelength(171*u.angstrom), a.Sample(10*u.minute))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    289 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 19.591 Gbyte
    <BLANKLINE>
           Start Time       ...
                            ...
    ----------------------- ...
    2012-03-04 00:00:00.000 ...
    2012-03-04 00:10:00.000 ...
    2012-03-04 00:20:00.000 ...
    2012-03-04 00:30:00.000 ...
    2012-03-04 00:40:00.000 ...
    2012-03-04 00:50:00.000 ...
    2012-03-04 01:00:00.000 ...
    2012-03-04 01:10:00.000 ...
    2012-03-04 01:20:00.000 ...
    2012-03-04 01:30:00.000 ...
                        ... ...
    2012-03-05 22:30:00.000 ...
    2012-03-05 22:40:00.000 ...
    2012-03-05 22:50:00.000 ...
    2012-03-05 23:00:00.000 ...
    2012-03-05 23:10:00.000 ...
    2012-03-05 23:20:00.000 ...
    2012-03-05 23:30:00.000 ...
    2012-03-05 23:40:00.000 ...
    2012-03-05 23:50:00.000 ...
    2012-03-06 00:00:00.000 ...
    Length = 289 rows
    <BLANKLINE>
    <BLANKLINE>

To search for data from multiple instruments, wavelengths, times etc., use the pipe ``|`` operator.
This joins queries together just as the logical ``OR`` operator would::

    >>> Fido.search(a.Time('2012/3/4', '2012/3/4 02:00'),
    ...             a.Instrument.lyra | a.Instrument.rhessi)  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 3 Providers:
    <BLANKLINE>
    2 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA ...      ESA     2
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA ...      ESA     3
    <BLANKLINE>
    1 Results from the RHESSIClient:
    Source: https://hesperia.gsfc.nasa.gov/hessidata
    <BLANKLINE>
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999     RHESSI ... RHESSI     NASA
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        ...   Size        Info
                                                    ...  Mibyte
    ----------------------- ----------------------- ... -------- --------------
    2012-03-03 22:57:40.000 2012-03-04 00:33:20.000 ... -0.00098 RHESSI level-0
    2012-03-04 00:33:20.000 2012-03-04 01:45:40.000 ... -0.00098 RHESSI level-0
    2012-03-04 01:45:40.000 2012-03-04 02:09:00.000 ... -0.00098 RHESSI level-0
    <BLANKLINE>
    <BLANKLINE>


Working with Search Results
***************************

:meth:`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>` can make multiple queries to multiple clients in one search.
This means that the results of a call to search can contain many sets of records, called responses, from many clients.
The results of a search are represented in a `~sunpy.net.fido_factory.UnifiedResponse` object, which provides access to all the response tables and allows some operations to be performed on all the results at once.
`~sunpy.net.fido_factory.UnifiedResponse` acts both like a two dimensional array, where the first dimension is the response index and the second index is the row index, and a dictionary where you can index the responses by the name of the client.

For example, the following code returns a response containing LYRA data from the `~sunpy.net.dataretriever.LYRAClient`, and EVE data from the `~sunpy.net.vso.VSOClient`::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Level.two,
    ...                       a.Instrument.lyra | a.Instrument.eve)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    2 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2
    <BLANKLINE>
    50 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        ...   Size         Info
                                                    ...  Mibyte
    ----------------------- ----------------------- ... -------- ----------------
    2012-01-01 00:00:00.000 2012-01-01 01:00:00.000 ... -0.00098 L2Lines (merged)
    2012-01-01 00:00:00.000 2012-01-01 01:00:00.000 ... -0.00098 L2Spectra (MEGS)
    2012-01-01 01:00:00.000 2012-01-01 02:00:00.000 ... -0.00098 L2Lines (merged)
                        ...                     ... ...      ...              ...
    2012-01-01 23:00:00.000 2012-01-02 00:00:00.000 ... -0.00098 L2Spectra (MEGS)
    2012-01-02 00:00:00.000 2012-01-02 01:00:00.000 ... -0.00098 L2Lines (merged)
    2012-01-02 00:00:00.000 2012-01-02 01:00:00.000 ... -0.00098 L2Spectra (MEGS)
    Length = 50 rows
    <BLANKLINE>
    <BLANKLINE>


If you then wanted to inspect just the LYRA data for the whole time range specified in the search, you would index this response to see just the results returned by the `~sunpy.net.dataretriever.LYRAClient`::

    >>> results[0, :]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2

Or, equivalently::

    >>> results["lyra"]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2

Normal slicing operations work as with any other Python sequence, e.g. ``results[1,::10]`` to access every tenth file in the result returned by the second client.

Note that the first (response) index is still necessary even if results are only found for a single client.
So in this case the first result would be ``results[0, 0]`` rather than ``results[0]`` (the latter would return all results from the first - and only - client and is therefore the same as ``results``).

Working with Response Tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As we have seen above the `~sunpy.net.fido_factory.UnifiedResponse` object contains many response tables which make up the search results.
Each of the responses are `~sunpy.net.base_client.QueryResponseTable` objects, which are `astropy.table` objects meaning that you can interact with them and filter them like any other tabular data.
This can be used to interact with results which are metadata only, i.e. searches from the HEK, or it can be used to reduce the number of files downloaded by `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`.

For example if we did a query for some AIA and HMI data::

    >>> results = Fido.search(a.Time("2020/01/01", "2020/01/01 00:05"), a.Instrument.aia | a.Instrument.hmi)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    201 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 13.626 Gbyte
    <BLANKLINE>
           Start Time       ...
                            ...
    ----------------------- ...
    2020-01-01 00:00:00.000 ...
    2020-01-01 00:00:04.000 ...
    2020-01-01 00:00:05.000 ...
    2020-01-01 00:00:05.000 ...
    2020-01-01 00:00:06.000 ...
    2020-01-01 00:00:09.000 ...
    2020-01-01 00:00:09.000 ...
    2020-01-01 00:00:11.000 ...
    2020-01-01 00:00:12.000 ...
    2020-01-01 00:00:14.000 ...
                        ... ...
    2020-01-01 00:04:47.000 ...
    2020-01-01 00:04:48.000 ...
    2020-01-01 00:04:52.000 ...
    2020-01-01 00:04:52.000 ...
    2020-01-01 00:04:53.000 ...
    2020-01-01 00:04:54.000 ...
    2020-01-01 00:04:57.000 ...
    2020-01-01 00:04:57.000 ...
    2020-01-01 00:04:59.000 ...
    2020-01-01 00:05:00.000 ...
    Length = 201 rows
    <BLANKLINE>
    21 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        ...            Info
                                                    ...
    ----------------------- ----------------------- ... --------------------------
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000 ... 45sec. Continuum intensity
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000 ...         45sec. Magnetogram
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000 ...         45sec. Dopplergram
    2020-01-01 00:01:07.000 2020-01-01 00:01:08.000 ... 45sec. Continuum intensity
    2020-01-01 00:01:07.000 2020-01-01 00:01:08.000 ...         45sec. Magnetogram
    2020-01-01 00:01:07.000 2020-01-01 00:01:08.000 ...         45sec. Dopplergram
    2020-01-01 00:01:52.000 2020-01-01 00:01:53.000 ... 45sec. Continuum intensity
    2020-01-01 00:01:52.000 2020-01-01 00:01:53.000 ...         45sec. Magnetogram
    2020-01-01 00:01:52.000 2020-01-01 00:01:53.000 ...         45sec. Dopplergram
    2020-01-01 00:02:37.000 2020-01-01 00:02:38.000 ... 45sec. Continuum intensity
    2020-01-01 00:02:37.000 2020-01-01 00:02:38.000 ...         45sec. Magnetogram
    2020-01-01 00:02:37.000 2020-01-01 00:02:38.000 ...         45sec. Dopplergram
    2020-01-01 00:03:22.000 2020-01-01 00:03:23.000 ... 45sec. Continuum intensity
    2020-01-01 00:03:22.000 2020-01-01 00:03:23.000 ...         45sec. Magnetogram
    2020-01-01 00:03:22.000 2020-01-01 00:03:23.000 ...         45sec. Dopplergram
    2020-01-01 00:04:07.000 2020-01-01 00:04:08.000 ... 45sec. Continuum intensity
    2020-01-01 00:04:07.000 2020-01-01 00:04:08.000 ...         45sec. Magnetogram
    2020-01-01 00:04:07.000 2020-01-01 00:04:08.000 ...         45sec. Dopplergram
    2020-01-01 00:04:52.000 2020-01-01 00:04:53.000 ... 45sec. Continuum intensity
    2020-01-01 00:04:52.000 2020-01-01 00:04:53.000 ...         45sec. Magnetogram
    2020-01-01 00:04:52.000 2020-01-01 00:04:53.000 ...         45sec. Dopplergram
    <BLANKLINE>
    <BLANKLINE>

The VSO client returns a lot of information about the records, so the first thing we can do is show only the columns we are interested in.
We can inspect all the available column names in all the responses with the `~.UnifiedResponse.all_colnames` property::

    >>> results.all_colnames  # doctest: +REMOTE_DATA
    ['End Time', 'Extent Length', 'Extent Type', 'Extent Width', 'Info', 'Instrument', 'Physobs', 'Provider', 'Size', 'Source', 'Start Time', 'Wavelength', 'Wavetype', 'fileid']

And then we can pick which ones to see with the :meth:`~.UnifiedResponse.show` method::

    >>> results.show("Start Time", "Instrument", "Physobs", "Wavelength")  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    201 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time       Instrument  Physobs   Wavelength [2]
                                                     Angstrom
    ----------------------- ---------- --------- ----------------
    2020-01-01 00:00:00.000        AIA intensity   335.0 .. 335.0
    2020-01-01 00:00:04.000        AIA intensity   193.0 .. 193.0
    2020-01-01 00:00:05.000        AIA intensity   304.0 .. 304.0
    2020-01-01 00:00:05.000        AIA intensity 4500.0 .. 4500.0
    2020-01-01 00:00:06.000        AIA intensity   131.0 .. 131.0
    2020-01-01 00:00:09.000        AIA intensity   171.0 .. 171.0
    2020-01-01 00:00:09.000        AIA intensity   211.0 .. 211.0
    2020-01-01 00:00:11.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:00:12.000        AIA intensity   335.0 .. 335.0
    2020-01-01 00:00:14.000        AIA intensity 1600.0 .. 1600.0
                        ...        ...       ...              ...
    2020-01-01 00:04:47.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:04:48.000        AIA intensity   335.0 .. 335.0
    2020-01-01 00:04:52.000        AIA intensity 1700.0 .. 1700.0
    2020-01-01 00:04:52.000        AIA intensity   193.0 .. 193.0
    2020-01-01 00:04:53.000        AIA intensity   304.0 .. 304.0
    2020-01-01 00:04:54.000        AIA intensity   131.0 .. 131.0
    2020-01-01 00:04:57.000        AIA intensity   171.0 .. 171.0
    2020-01-01 00:04:57.000        AIA intensity   211.0 .. 211.0
    2020-01-01 00:04:59.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:05:00.000        AIA intensity   335.0 .. 335.0
    Length = 201 rows
    <BLANKLINE>
    21 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time       Instrument      Physobs        Wavelength [2]
                                                              Angstrom
    ----------------------- ---------- ------------------ ----------------
    2020-01-01 00:00:22.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:00:22.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:00:22.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:01:07.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:01:07.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:01:07.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:01:52.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:01:52.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:01:52.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:02:37.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:02:37.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:02:37.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:03:22.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:03:22.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:03:22.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:04:07.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:04:07.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:04:07.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:04:52.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:04:52.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:04:52.000        HMI       LOS_velocity 6173.0 .. 6174.0
    <BLANKLINE>
    <BLANKLINE>

To give an example of filtering post-search, let's only return the rows in the table which are line-of-sight magnetograms from HMI or the 94Å passband from AIA.
You can also always do this filtering with the `a.vso.Physobs <sunpy.net.attrs.Physobs>` and `a.Wavelength <sunpy.net.attrs.Wavelength>` attrs in the search command.

First we split the results in to a table for AIA and a table for HMI::

   >>> aia, hmi = results  # doctest: +REMOTE_DATA

We can use boolean indexing to match the value of the ``"Physobs"`` column::

    >>> hmi_los = hmi[hmi["Physobs"] == "LOS_magnetic_field"]  # doctest: +REMOTE_DATA
    >>> hmi_los.show("Start Time", "Instrument", "Wavelength", "Physobs")  # doctest: +REMOTE_DATA
    <sunpy.net.vso.table_response.VSOQueryResponseTable object at ...>
           Start Time       Instrument  Wavelength [2]       Physobs
                                          Angstrom
    ----------------------- ---------- ---------------- ------------------
    2020-01-01 00:00:22.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:01:07.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:01:52.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:02:37.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:03:22.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:04:07.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:04:52.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field

To match the ``"Wavelength"`` column we need to account for the fact that VSO results return a wavelength range of ``[min, max]`` so we match the min::

    >>> aia_94 = aia[aia["Wavelength"][:, 0] == 94 * u.AA]  # doctest: +REMOTE_DATA
    >>> aia_94.show("Start Time", "Instrument", "Wavelength", "Physobs")  # doctest: +REMOTE_DATA
    <sunpy.net.vso.table_response.VSOQueryResponseTable object at ...>
           Start Time       Instrument Wavelength [2]  Physobs
                                          Angstrom
    ----------------------- ---------- -------------- ---------
    2020-01-01 00:00:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:59.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:59.000        AIA   94.0 .. 94.0 intensity
                        ...        ...            ...       ...
    2020-01-01 00:03:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:59.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:59.000        AIA   94.0 .. 94.0 intensity
    Length = 25 rows

These can then be passed to `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> Fido.fetch(hmi_los, aia_94)  # doctest: +SKIP

.. warning::

   While you can reduce the number of columns and rows in the results, the
   ``fetch()`` method may need certain columns to be present to successfully
   download the files. It is therefore highly recommended that if you are
   planning on downloading data you do not slice out columns, but instead use
   ``.show()`` to only display the ones you are interested in.


.. _downloading_data:

Downloading data
****************
Once you have located your files via a `Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>`, you can download them via `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(results)  # doctest: +SKIP

This downloads the files to the location set in you sunpy config file.
It also returns a `parfive.Results` object ``downloaded_files``, of absolute file paths of where the files have been downloaded to.

You can also specify the path to which you want the data downloaded::

  >>> downloaded_files = Fido.fetch(results, path='/ThisIs/MyPath/to/Data/{file}')  # doctest: +SKIP

This downloads the query results into the directory ``/ThisIs/MyPath/to/Data``, naming each downloaded file with the filename ``{file}`` obtained from the client.
You can also use other properties of the returned query to define the path where the data is saved.
For example, to save the data to a subdirectory named after the instrument, use::

    >>> downloaded_files = Fido.fetch(results, path='./{instrument}/{file}')  # doctest: +SKIP

You can see the list of options that can be specified in path for all the files to be downloaded with ``results.path_format_keys``.

Retrying Downloads
^^^^^^^^^^^^^^^^^^

If any files failed to download, the progress bar will show an incomplete number of files (i.e. 100/150) and the `parfive.Results` object will contain a list of the URLs that failed to transfer and the error associated with them.
This can be accessed with the ``.errors`` attribute or by printing the `~parfive.Results` object::

    >>> print(downloaded_files.errors)  # doctest: +SKIP

The transfer can be retried by passing the `parfive.Results` object back to `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(downloaded_files)  # doctest: +SKIP

doing this will append any newly downloaded file names to the list and replace the ``.errors`` list with any errors that occurred during the second attempt.


.. _VSO: https://sdac.virtualsolar.org/cgi/search
.. _JSOC: http://jsoc.stanford.edu/


Fido Clients
************

`~sunpy.net.Fido` provides access to many sources of data via "clients", these clients can be defined inside sunpy or in other packages.
If you want to see the current list of clients you can do::

    >>> print(Fido)
    sunpy.net.Fido
    <BLANKLINE>
    Fido is a unified data search and retrieval tool.
    <BLANKLINE>
    It provides simultaneous access to a variety of online data sources, some
    cover multiple instruments and data products like the Virtual Solar
    Observatory and some are specific to a single source.
    <BLANKLINE>
    For details of using `~sunpy.net.Fido` see :ref:`fido_guide`.
    <BLANKLINE>
    <BLANKLINE>
          Client                                                    Description
    ----------------- -------------------------------------------------------------------------------------------------------
    CDAWEBClient      Provides access to query and download from the Coordinated Data Analysis Web (CDAWeb).
    EVEClient         Provides access to Level 0C Extreme ultraviolet Variability Experiment (EVE) data.
    GBMClient         Provides access to data from the Gamma-Ray Burst Monitor (GBM) instrument on board the Fermi satellite.
    XRSClient         Provides access to the GOES XRS fits files archive.
    SUVIClient        Provides access to data from the GOES Solar Ultraviolet Imager (SUVI).
    GONGClient        Provides access to the Magnetogram products of NSO-GONG synoptic Maps.
    LYRAClient        Provides access to the LYRA/Proba2 data archive.
    NOAAIndicesClient Provides access to the NOAA solar cycle indices.
    NOAAPredictClient Provides access to the NOAA SWPC predicted sunspot Number and 10.7 cm radio flux values.
    SRSClient         Provides access to the NOAA SWPC solar region summary data.
    NoRHClient        Provides access to the Nobeyama RadioHeliograph (NoRH) averaged correlation time series data.
    RHESSIClient      Provides access to the RHESSI observing summary time series data.
    HEKClient         Provides access to the Heliophysics Event Knowledgebase (HEK).
    HECClient         Provides access to the HELIO webservices.
    JSOCClient        Provides access to the JSOC Data Export service.
    VSOClient         Provides access to query and download from Virtual Solar Observatory (VSO).
.. _sample-data:

***********************
Downloading Sample Data
***********************

sunpy provides a number of sample data files which are referenced in the
documentation and examples. These files are available to download onto your
local machine so that you can try out the code in the documentation. To
download the sample data simply run the following command::

    import sunpy.data.sample

This will download the data to your sample-data directory which can be
customized by editing the sunpyrc file (see :doc:`../customization`).
After running this you can then import the sample data files shortcuts which
are used below.
.. _acquiring_data:

*************************
Acquiring Data with sunpy
*************************

In this section of the guide we will introduce the ways you can obtain different
kind of solar data from different places.

.. toctree::
    :maxdepth: 2

    sample-data
    fido
    jsoc
    hek
    helioviewer
    database
.. _dev_guide:
{% if not is_development %}
.. _newcomers:
.. _remote_data:
{% endif %}

*****************
Developer's Guide
*****************

{% if is_development %}

This article describes the guidelines to be followed by developers working on sunpy.
If you are thinking of contributing to sunpy please read the following carefully.

We currently recommend the :ref:`newcomers` as the place to start.
This goes over the basics and has links to useful tutorials on git.

.. toctree::
   :maxdepth: 2

   contents/newcomers
   contents/code_standards
   contents/documentation
   contents/example_gallery
   contents/tests
   contents/pr_review_procedure
   contents/api
   contents/dependencies
   contents/units_quantities
   contents/new_objects
   contents/extending_fido
   contents/maintainer_workflow
   contents/logger
   contents/remote_data
   contents/config
   contents/funding

{%else%}

Please go `here <https://docs.sunpy.org/en/latest/dev_guide/index.html>`__ for our up to date developer's guide.

{%endif%}
.. _extending_fido:

***************************************
Extending Fido with New Sources of Data
***************************************

Sunpy's data search and retrieval tool (``Fido``) is designed to be extensible, so that new sources of data or metadata can be supported, either inside or outside the sunpy core package.
There are two ways of defining a new client, depending on the complexity of the web service.
A "scraper" client inherits from `~sunpy.net.dataretriever.client.GenericClient` which provides helper methods for downloading from a list of URLs.
If the service you want to add has easily accesible HTTP or FTP URLs that have a well defined folder and filename structure, this is probably the best approach.
If the service you want to add requires making requests to an API with parameters for the search and getting a list of results in return, then you probably want to write a "full" client.

Before writing a new client, ensure you are familiar with how searches are specified by the `sunpy.net.attr` system, including combining them with logical operations.
When choosing a name for your new client it should have the form ``<name>Client`` as sunpy will split the name the name of the class to extract the name of your client.
The main place this is done is when constructing a `~.UnifiedResponse` object, where the name part can be used to index the response object.


.. _new_scraper_client:

Writing a new "scraper" client
==============================

A "scraper" Fido client (also sometimes referred to as a "data retriever" client) is a Fido client which uses the URL `~sunpy.net.scraper.Scraper` to find files on remote servers.
If the data provider you want to integrate does not provide a tree of files with predictable URLs then a "full" client is more likely to provide the functionality you need.

A new "scraper" client inherits from `~sunpy.net.dataretriever.client.GenericClient` and requires a minimum of these three components:

* A class method :meth:`~sunpy.net.base_client.BaseClient.register_values`; this registers the "attrs" that are supported by the client.
  It returns a dictionary where keys are the supported attrs and values are lists of tuples.
  Each ``tuple`` contains the attr value and its description.
* A class attribute ``baseurl``; this is a regular expression which is used to match the URLs supported by the client.
* A class attribute ``pattern``; this is a template used to extract the metadata from URLs matched by ``baseurl``.
  The extraction uses the `~sunpy.extern.parse.parse` format.


For a simple example of a scraper client, we can look at the implementation of `sunpy.net.dataretriever.sources.eve.EVEClient` in sunpy.
A version without documentation strings is reproduced below:

.. code-block:: python

  class EVEClient(GenericClient):
      baseurl = (r'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/'
                 r'L0CS/SpWx/%Y/%Y%m%d_EVE_L0CS_DIODES_1m.txt')
      pattern = '{}/SpWx/{:4d}/{year:4d}{month:2d}{day:2d}_EVE_L{Level:1d}{}'

      @classmethod
      def register_values(cls):
          from sunpy.net import attrs
          adict = {attrs.Instrument: [('EVE', 'Extreme ultraviolet Variability Experiment, which is part of the NASA Solar Dynamics Observatory mission.')],
                   attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
                   attrs.Source: [('SDO', 'The Solar Dynamics Observatory.')],
                   attrs.Provider: [('LASP', 'The Laboratory for Atmospheric and Space Physics.')],
                   attrs.Level: [('0', 'EVE: The specific EVE client can only return Level 0C data. Any other number will use the VSO Client.')]}
          return adict


This client scrapes all the URLs available under the base url ``http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/``.
`~sunpy.net.scraper.Scraper` is primarily focused on URL parsing based on time ranges, so the rest of the ``baseurl`` pattern specifies where in the pattern the time information is located, using `strptime <https://strftime.org/>`__ notation.
The ``pattern`` attribute is used to populate the results table from the URLs matched with the ``baseurl``.
It includes some of the time definitions, as well as names of attrs (in this case "Level").
The supported time keys are: 'year', 'month', 'day', 'hour', 'minute', 'second', 'millisecond'.

The attrs returned in the ``register_values()`` method are used to match your client to a search, as well as adding their values to the attr.
This means that after this client has been imported, running ``print(a.Provider)`` will show that the ``EVEClient`` has registered a provider value of ``LASP``.
In addition to this, a sanitized, lower cased version of the value will be available for tab completing, e.g. ``a.Provider.lasp`` or ``a.Level.zero``.


More Complex Clients
--------------------

Sometimes the attr values may not exist identically in the required URLs, and therefore can not be simply extracted with ``pattern``.
Say, for example, the Wavelength of a file is expressed in the URL as a passband by name; in this case conversion of the `~astropy.units.Quantity` object to the pass band name would be needed.
This is done addressed with the two following methods:

* :meth:`~sunpy.net.dataretriever.client.GenericClient.pre_search_hook` which will convert the passed attrs to their representation in the URL.
* :meth:`~sunpy.net.dataretriever.client.GenericClient.post_search_hook` which converts the retrieved metadata from a URL to the form in which they are desired to be represented in the response table.

A good example of the use of these two methods is the `sunpy.net.dataretriever.sources.norh.NoRHClient` in sunpy.

It may also be possible that the ``baseurl`` property needs to be customised based on attrs other than Time.
Since `~sunpy.net.scraper.Scraper` doesn't currently support generating directories that have non-time variables, the :meth:`~sunpy.net.dataretriever.client.GenericClient.search` needs to be customised.
The search method should in this case, generate a ``baseurl`` dependant on the values of these attrs, and then call ``super().search`` or `~sunpy.net.scraper.Scraper` for each ``baseurl`` generated.
For an example of a complex modification of the ``search()`` method see the implementation of `.SUVIClient.search`.

Examples
--------

Suppose any file of a data archive can be described by this URL ``https://some-domain.com/%Y/%m/%d/satname_{SatelliteNumber}_{Level}_%y%m%d%H%M%S_{any-2-digit-number}.fits``:

``baseurl`` becomes ``r'https://some-domain.com/%Y/%m/%d/satname_(\d){2}_(\d){1}_(\d){12}_(\d){2}\.fits'``.

Note all variables in the filename are converted to regex that will match any possible value for it.
A character enclosed within ``()`` followed by a number enclosed within ``{}`` is used to match the specified number of occurences of that special sequence.
For example, ``%y%m%d%H%M%S`` is a twelve digit variable (with 2 digits for each item) and thus represented by ``r'(\d){12}'``.
Note that ``\`` is used to escape the special character ``.``.

``pattern`` becomes ``'{}/{year:4d}/{month:2d}{day:2d}/satname_{SatelliteNumber:2d}_{Level:1d}_{:6d}{hour:2d}{minute:2d}{second:2d}_{:2d}.fits'``.
Note the sole purpose of ``pattern`` is to extract the information from matched URL, using `~sunpy.extern.parse.parse`.
So the desired key names for returned dictionary should be written in the ``pattern`` within ``{}``, and they should match with the ``attr.__name__``.

``register_values()`` can be written as:

.. code-block:: python

    @classmethod
    def register_values(cls):

        from sunpy.net import attrs
        adict = {
        attrs.Instrument: [("SatName", "The description of Instrument")],
        attrs.Physobs: [('some_physobs', 'Phsyobs description')],
        attrs.Source: [('some_source', 'Source description')],
        attrs.Provider: [('some_provider', 'Provider description')],
        attrs.Level: [("1", "Level 1 data"), ("2", "Level 2 data")],
        attrs.SatelliteNumber: [("16", "Describe it"), ("17", "Describe it")]
        }

        return adict


.. _new_full_client:

Writing a "full" client
=======================

In this section we will describe how to build a "full" Fido client.
You should write a new "full" client if the data you are accessing can not be accessed via a URL template, for instance if you hit a web API with a query to return results for a search.

A new Fido client contains three major components:

* A subclass of `~sunpy.net.base_client.BaseClient` which implements ``search``, ``fetch``, and ``_can_handle_query``.
* Zero or more new `~sunpy.net.attr.Attr` classes to specify search parameters unique to your data source.
* An instance of `~sunpy.net.attr.AttrWalker` which can be used to walk the tree of `~sunpy.net.attr.Attr` instances and convert them into a form useful to your client's search method.

Search Attrs
------------

As described in `~sunpy.net.attr` the attr system allows the construction of complex queries by the user.
To make these complex queries easily processable by the clients the ``AttrWalker`` converts these into a set of queries which can be processed separately.
It does this by converting the input query to a set of queries which are ORed, but are complete queries.
This means the list of queries is an **OR** of **ANDs** (technically called `disjuntive normal form <https://en.wikipedia.org/wiki/Disjunctive_normal_form>`__).

Each query in the list of ORs contains all the information about that query so for example if the user provided a query like::

  a.Time("2020/02/02", "2020/02/03") & (a.Instrument("AIA") | a.Instrument("HMI"))

it would be passed to the client as::

  (a.Time("2020/02/02", "2020/02/03") & a.Instrument("HMI")) | (a.Time("2020/02/02", "2020/02/03") & a.Instrument("AIA"))

So you can process each element of the OR in turn without having to consult any other part of the query.

If the query the user provided contains an OR statement you get passed an instance of `~sunpy.net.attr.AttrOr` and each sub-element of that `~sunpy.net.attr.AttrOr` will be `~sunpy.net.attr.AttrAnd` (or a single other attr class).
If the user query dosen't contain an OR you get a single `~.Attr` instance or an `~.AttrAnd`.

For example you could get any of the following queries (using ``&`` for AND and ``|`` for OR):

* ``(a.Instrument("AIA") & a.Time("2020/02/02", "2020/02/03")) | (a.Instrument("HMI") & a.Time("2020/02/02", "2020/02/03"))``
* ``a.Time("2020/02/02", "2020/02/03")``
* ``a.Instrument("AIA") & a.Time("2020/02/02", "2020/02/03")``
* ``(a.Time(..) & a.Instrument("AIA") & a.Wavelength(30*u.nm, 31*u.nm)) | (a.Time(..) & a.Instrument("AIA") & a.Wavelength(30*u.nm, 31*u.nm))``

but you **would not** be passed queries which look like the following examples, even if that's how the user specified them:

* ``a.Time("2020/02/02", "2020/02/03") & (a.Instrument("AIA") | a.Instrument("HMI"))``
* ``a.Time(..) & (a.Instrument("AIA") | a.Instrument("AIA")) & a.Wavelength(30*u.nm, 31*u.nm))``

The Attr Walker
###############

Given the potential complexity of these combined attrs, converting them into other forms, such as query parameters or JSON etc involves walking the tree and converting each attr to the expected format in a given way.
This parsing and conversion of the query tree is deliberately not done using methods or attributes of the attrs themselves.
The attrs should be independent of any client in their implementation, so they can be shared between the different ``Fido`` clients.

A class is provided to facilitate this conversion, `~sunpy.net.attr.AttrWalker`.
The `~sunpy.net.attr.AttrWalker` class consists of three main components:

* **Creators**: The `~sunpy.net.attr.AttrWalker.create` method is one of two generic functions for which a different function is called for each Attr type.
  The intended use for creators is to return a new object dependant on different attrs.
  It is commonly used to dispatch on `~sunpy.net.attr.AttrAnd` and `~sunpy.net.attr.AttrOr`.

* **Appliers**: The `~sunpy.net.attr.AttrWalker.apply` method is the same as `~sunpy.net.attr.AttrWalker.create` in that it is a generic function.
  The only difference between it and `~sunpy.net.attr.AttrWalker.create` is its intended use.
  Appliers are generally used to modify an object returned by a creator with the values or information contained in other Attrs.

* **Converters**: Adding a converter to the walker adds the function to both the creator and the applier.
  For the VSO client this is used to convert each supported attr into a `~sunpy.net.attr.ValueAttr` which is then later processed by the appliers and creators.
  This pattern can be useful if you would otherwise have to repeat a lot of logic in each of the applier functions for each type of Attr you support.

An Example of ``AttrWalker``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we will write a parser for some simple queries which uses `~sunpy.net.attr.AttrWalker` to convert the query to a `dict` of URL query parameters for a HTTP GET request.
Let's imagine we have a web service which you can do a HTTP GET request to ``https://sfsi.sunpy.org/search`` for some imaginary data from an instrument called SFSI (Sunpy Fake Solar Instrument).
This GET request takes three query parameters ``startTime``, ``endTime`` and ``level``, so a request might look something like: ``https://sfsi.sunpy.org/search?startTime=2020-01-02T00:00:00&endTime=2020-01-02T00:00:00&level=1``.
Which would search for level one data between 2020-01-01 and 2020-01-02.

As `~sunpy.net.attrs` has `~sunpy.net.attrs.Time` and `~sunpy.net.attrs.Level` we don't need to define any of our own attrs for this client.
We do however want to write our own walker to convert them to the form out client's ``search()`` method wants to send them to the server.

The first step is to setup the walker and define a creator method which will return a list of dicts, one for each independent search.

.. code-block:: python

    import sunpy.net.attrs as a
    from sunpy.net.attr import AttrWalker, AttrAnd, AttrOr, DataAttr

    walker = AttrWalker()

    @walker.add_creator(AttrOr)
    def create_or(wlk, tree):
        results = []
        for sub in tree.attrs:
            results.append(wlk.create(sub))

        return results

    @walker.add_creator(AttrAnd, DataAttr)
    def create_and(wlk, tree):
        result = dict()
        wlk.apply(tree, result)
        return [result]


The call ``wlk.apply(...)`` inside the creator will walk any nested attrs and add their values to the dictionary as defined by the applier registered to each attr type.
If we want our client to support searching by ``a.Time`` and ``a.Level`` as in the URL example above, we would need to register an applier for each of these attrs.

.. code-block:: python

    @walker.add_applier(a.Time)
    def _(wlk, attr, params):
        return params.update({'startTime': attr.start.isot,
                              'endTime': attr.end.isot})

    @walker.add_applier(a.Level)
    def _(wlk, attr, params):
        return params.update({'level': attr.value})


This combination of creators and appliers would allow support of any combination of queries consisting of ``a.Time`` and ``a.Level``.
Obviously, most clients would want to support more attrs than these two, and this could be done by adding more applier functions.

Adding "Attrs" to Registry
##########################

Registering of "attrs" ensures discoverability of search attributes supported by the corresponding sunpy Client.
For adding them to the Registry, we need to define a ``classmethod`` :meth:`~sunpy.net.base_client.BaseClient.register_values` that returns a dictionary of registered values.
This dictionary should have `~sunpy.net.attr.Attr` classes as keys and a list of tuples corresponding to that key representing the possible values the key "attr" can take.
Each tuple comprises of two elements.
The first one is a value and the second element contains a brief description of that value.
An example of writing ``register_values()`` for `~sunpy.net.dataretriever.client.GenericClient` is provided above.
Please note that it can be defined in a similar way for full clients too.

An Example of ``register_values()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    @classmethod
    def register_values(cls):

        from sunpy.net import attrs
        adict = {
        attrs.Instrument: [("LASCO", "Large Angle and Spectrometric Coronagraph")],
        attrs.Source: [('SOHO', 'Solar and Heliospheric Observatory')],
        attrs.Provider: [('SDAC', 'Solar Data Analysis Center')],
        attrs.Detector: [('C1', 'Coronograph 1'),
                         ('C2', 'Coronograph 2'),
                         ('C3', 'Coronograph 3')]
        }

        return adict

Writing a Search Method
-----------------------

The ``search()`` method has the job of taking a set of user queries and returning an instance of `.QueryResponseTable` containing the results.

The general flow of a ``search()`` method is:

* Call your instance of an `.AttrWalker` to convert the input into a form expected by your API.
* Make as many requests to your API as needed to fulfill the query. (Generally one per element of the outer `.AttrOr`).
* Process the response from your API into an instance of `.QueryResponseTable`.

To process the query with the `.AttrWalker`, call the :meth:`.AttrWalker.create` method::

  def search(self, query):
    queries = walker.create(query)

Assuming the walker is the one we defined above, queries would be a list of dicts with the attrs processed into query parameters for the API URL.

.. note::

   If you want your search method to be able to be called independently of Fido, then you should accept a variable number of positional arguments (``*args``) and they should have the AND operator applied to them.
   This looks like::

     def search(self, *args):
         query = attr.and_(args)
         queries = walker.create(query)


Once the walker has processed the query into a form designed to be passed to your API, your ``search()`` method then needs to iterate over these parameters, make the requests, and process the results into a table.

In the following example we pretend our client has a method ``_make_search(query_parameters)`` which takes the query parameters and makes a request to our API.
We also pretend that the response is a json object in the form of a Python dictionary, which we want to put into the table.

.. code-block:: python

  def search(self, query):
      queries = walker.create(query)

      results = []
      for query_parameters in queries:
          results.append(self._make_search(query_parameters))

      return QueryResponseTable(results, client=self)


In reality, you probably want to post-process the results from your API before you put them in the table, they should be human readable first, with spaces and capitalisation as appropriate.


Supporting filesize estimates
#############################
The base client has a method for automatically estimating the total size of files in a given query: :meth:`~sunpy.net.base_client.QueryResponseTable.total_size`.
To enable to support for this, make sure the table returned by ``search`` has a column that contains filesizes as astropy quantities convertible to ``u.byte``, and set the ``size_column`` class attribute to the name of this column.


The ``_can_handle_query`` method
---------------------------------

The next required method is ``_can_handle_query``, this method tells ``Fido`` if your client might be able to return results for a given query.
If this method returns `True`, your clients ``search()`` method will be called for that query.
This method gets passed each query (in its independent form), and must either return ``True`` or ``False``.

A simple example, which just checks the type of ``attrs`` and not their values would be::

  @classmethod
  def _can_handle_query(cls, *query):
      query_attrs = set(type(x) for x in query)
      supported_attrs = {a.Time, a.Level}
      return supported_attrs.issuperset(query_attrs)

Note, that this method is a class method, it gets called without instantiating your client to speed up the dispatching.


Writing a Fetch Method
----------------------

The ``fetch()`` method of a Fido client is responsible for converting a set of search results (possibly sliced by the user) into a set of URLs to be downloaded.
Due to the history of clients and how they were implemented in sunpy, some existing clients support use outside of the``Fido`` wrapper, this makes them appear more complex.
In this example we are going to write a ``fetch()`` method which is designed only to be called from ``Fido``.

The parameters for such a method should be::

  def fetch(self, query_results, *, path, downloader, **kwargs):

The parameters here are:

* ``query_results`` which is an instance of `~.QueryResponseTable` or `~sunpy.net.base_client.QueryResponseRow`, these are the results the user wants to download.
* ``path=`` This is the path that the user wants the file to be downloaded to, this can be a template string (i.e. expects to have ``.format()`` called on it).
* ``downloader=`` This is a `parfive.Downloader` object which should be mutated by the ``fetch()`` method.
* ``**kwargs`` It is very important that ``fetch()`` methods accept extra keyword arguments that they don't use, as the user might be passing them to other clients via ``Fido``.


Processing the ``query_results`` Argument
#########################################

The ``query_results`` argument can be of two types `~.QueryResponseTable` or `~sunpy.net.base_client.QueryResponseRow`, as the user can slice the results table down to a single row and then pass that to ``Fido.fetch()``.
If you do not wish to handle a single row any differently to a table, you can place the `~sunpy.net.base_client.convert_row_to_table` decorator on your ``fetch()`` method which will convert the argument to a length one table when it is a single row object.

The primary function of the ``fetch()`` method is for you to convert this results object into a set of URLs for Fido to download.
This logic will be specific to your client.


Formatting the ``path=`` Argument
#################################

The path argument may contain format sections which are processed column names from the response table.
In addition to these it may contain the ``{file}`` format segment which is a placeholder for the filename.
Each row of the results table has a `~sunpy.net.base_client.QueryResponseRow.response_block_map` property which is a dictionary of valid format keys to values for that row.

In addition to the `~sunpy.net.base_client.QueryResponseRow.response_block_map` your fetch method also needs to be able to generate a filename for the file.
The simplest (but unlikely) scenario is that you know the filename for each file you are going to download before you do so, in this situation you would be able to generate the full filepath for each row of the response as follows::

  for row in query_results:
      filename = self._calculate_filename(row)
      filepath = path.format(file=filename, **row.response_block_map)


In the situation where you wish to be told the filename by the webserver you are downloading the file from, it is a little more complex, you need to pass a callback function to :meth:`parfive.Downloader.enqueue_file` which will calculate the full filename in the context of the download, where the headers can be inspected for the filename the webserver provides.

The filename callback passed to :meth:`parfive.Downloader.enqueue_file` accepts two arguments ``resp`` and ``url``.
``resp`` is an `aiohttp.ClientResponse` object which is returned when `parfive` requests the URL.
This response object allows us to inspect the headers of the response before the data is downloaded.
``url`` is the URL that was requested to generate the ``resp`` response.

To combine the formatting of the row with the extraction of the filename from the headers it is common to use `functools.partial` to generate many functions with different fixed parameters.
In the following example we will define a function which takes 4 arguments which we will use to generate the filename for the row.
This function will be called by `parfive` with the ``resp`` and ``url`` arguments.::

  def make_filename(path, row, resp, url):
      # Define a fallback filename based on the information in the search results
      name = f"row['ID'].fits"

      if resp:
        cdheader = resp.headers.get("Content-Disposition", None)
        if cdheader:
          _, params = cgi.parse_header(cdheader)
          name = params.get('filename', "")

      return path.format(file=name, **row.response_block_map)

To reduce this function down to the two arguments expected we pre-specify the first two of these with `~functools.partial` before passing the function to `~parfive.Downloader.enqueue_file` inside the ``fetch()`` method.
Our simple example above now becomes::

  for row in query_results:
      filepath = partial(make_filename, path, row)

Where the ``path`` variable is a `pathlib.Path` object provided as the ``path`` argument to ``fetch()``.


Adding URLs to be Downloaded
############################

For each file you wish for ``Fido`` to download (normally one per row of the ``query_results``) you need to call the :meth:`parfive.Downloader.enqueue_file` of the ``downloader`` argument.
Combining this with the simple example above it may look something like::

  for row in query_results:
      filename = self._calculate_filename(row)
      filepath = path.format(file=filename, **row.response_block_map)

      url = self._calculate_url(row)
      downloader.enqueue_file(url, filename=filepath)


If your filepath is a callback function, pass this to the ``filename=`` argument.

Your fetch method does not need to return anything, as long as ``enqueue_file`` is called for every file you want ``Fido`` to download.

Putting it all Together
-----------------------

An example client class may look something like::

  import cgi

  import sunpy.net.atrrs as a
  from sunpy.net.attr import AttrWalker, AttrAnd, AttrOr, DataAttr
  from sunpy.base_client import QueryResponseTable


  walker = AttrWalker()


  @walker.add_creator(AttrOr)
  def create_or(wlk, tree):
      results = []
      for sub in tree.attrs:
          results.append(wlk.create(sub))

      return results


  @walker.add_creator(AttrAnd, DataAttr)
  def create_and(wlk, tree):
      result = dict()
      wlk.apply(tree, result)
      return [result]


  @walker.add_applier(a.Time)
  def _(wlk, attr, params):
      return params.update({'startTime': attr.start.isot,
                            'endTime': attr.end.isot})


  @walker.add_applier(a.Level)
  def _(wlk, attr, params):
      return params.update({'level': attr.value})


  class ExampleClient(BaseClient):
      size_column = 'Filesize'

      def search(self, query):
          queries = walker.create(query)

          results = []
          for query_parameters in queries:
              results.append(self._make_search(query_parameters))

          return QueryResponseTable(results, client=self)

      def _make_filename(path, row, resp, url):
          # Define a fallback filename based on the information in the search results
          name = f"row['ID'].fits"

          if resp:
            cdheader = resp.headers.get("Content-Disposition", None)
            if cdheader:
              _, params = cgi.parse_header(cdheader)
              name = params.get('filename', "")

          return path.format(file=name, **row.response_block_map)

      @convert_row_to_table
      def fetch(self, query_results, *, path, downloader, **kwargs):
          for row in query_results:
              filepath = partial(self._make_filename, path, row)

              url = f"https://sfsi.sunpy.org/download/{row['ID']}"
              downloader.enqueue_file(url, filename=filepath)

      @classmethod
      def _can_handle_query(cls, *query):
          query_attrs = set(type(x) for x in query)
          supported_attrs = {a.Time, a.Level}
          return supported_attrs.issuperset(query_attrs)
.. _testing:

******************
Testing Guidelines
******************

This section describes the testing framework and format standards for tests in sunpy.
Here we have heavily adapted the `Astropy version <https://docs.astropy.org/en/latest/development/testguide.html>`_, and **it is worth reading that link.**

The testing framework used by sunpy is the `pytest`_ framework, accessed through the ``pytest`` command.

.. _pytest: https://pytest.org/en/latest/

.. note::

    The ``pytest`` project was formerly called ``py.test``, and you may
    see the two spellings used interchangeably.

Dependencies for testing
------------------------

Since the testing dependencies are not actually required to install or use sunpy, they are not included in "install_requires" in "setup.cfg".

Developers who want to run the test suite will need to install the testing packages using pip::

    $ pip install -e .[tests]

If you want to see the current test dependencies, you check "extras_require" in "setup.cfg".

Running Tests
=============

There are currently two different ways to invoke the sunpy tests.
However, we strongly suggest using ``tox`` as the default one.
Each method uses the widely-used ``pytest`` framework and are detailed below.

``tox``
-------

The primary method is to use `tox`_, which is a generic virtualenv management and test command line tool.
We have several environments within our "tox.ini" file and you can list them::

    $ tox -v -l

Then you can run any of them doing::

    $ tox -e <name of env>

This will create a test environment in ".tox" and build, install sunpy and runs the entire test suite.
This is the method that our continuous integration uses.

.. _tox: https://tox.readthedocs.io/en/latest/

``pytest``
----------

The test suite can be run directly from the native ``pytest`` command.
In this case, it is important for developers to be aware that they must manually rebuild any extensions by running ``python setup.py build_ext`` before testing.

To run the entire suite with ``pytest``::

    $ pytest

will use the settings in ``setup.cfg``.

If you want to run one specific test file::

    $ pytest sunpy/map/tests/test_mapbase.py

or one specific test in a test file::

    $ pytest sunpy/map/tests/test_mapbase.py::<test_name>

(This does not work with ``tox`` and is a known issue.)

If a test errors, you can use ``pdb`` to create a debugging session at the moment the test fails::

    $ pytest --pdb

Remote data
-----------
By default, no online tests are selected and so to run the online tests you have to::

    $ tox -e py38-online

or::

    $ pytest --remote-data=any

Figure tests
------------
In order to avoid changes in figures due to different package versions, we recommend using tox to run the figure tests::

    $ tox -e py38-figure

This will ensure that any figures created are checked using the package versions that were used to create the original figure hashes.
Running this will create a folder, "figure_test_images", within your work folder ("<local clone location>/figure_test_images"), which is ignored by git.
Inside this folder will be all the images created, as well as a json file with the hashes of the figures created by the test run.
The current hashes are located within "sunpy/tests/figure_hashes_mpl_<ver>_ft_<ver>_astropy_<ver>.json" and this will be where you will need to update old hashes or create new figure entries if anything changes.
The filenames are the versions of Matplotlib, freetype and astropy used.
If these versions differ to your local setup, the figure tests will not run.
In theory, the Python version does not change the results as we have pinned the packages that cause the hash to vary.

Running tests in parallel
-------------------------

It is possible to speed up sunpy's tests using the `pytest-xdist`_ plugin.
This plugin can be installed using `pip`_::

    pip install pytest-xdist

Once installed, tests can be run in parallel using the ``--parallel`` commandline option.
For example, to use 4 processes::

    $ tox -e <name of environment> -- -n=4

or::

    $ pytest -n 4 ./sunpy

.. _pytest-xdist: https://pypi.python.org/pypi/pytest-xdist
.. _pip: https://pypi.org/project/pip/

Coverage reports
----------------

sunpy can use `pytest-cov`_  generate test coverage reports and settings are stored in ``setup.cfg``.
This plugin can be installed using `pip`_::

    $ pip install pytest-cov

To generate a test coverage report, use::

    $ pytest --cov ./sunpy

This will print to the terminal a report of line coverage of our test suite.
If you want to create a report in html, you can run::

    $ pytest --cov-report xml:cov.xml --cov ./sunpy
    $ coverage html

.. _pytest-cov: https://pypi.org/project/pytest-cov/

Writing tests
=============

``pytest`` has the following `test discovery rules <https://pytest.org/en/latest/goodpractices.html#conventions-for-python-test-discovery>`_::

 * ``test_*.py`` or ``*_test.py`` files
 * ``Test`` prefixed classes (without an ``__init__`` method)
 * ``test_`` prefixed functions and methods

We use the first one for our test files, ``test_*.py`` and we suggest that developers follow this.

A rule of thumb for unit testing is to have at least one unit test per public function.

Simple example
--------------

The following example shows a simple function and a test to test this
function::

    def func(x):
        """Add one to the argument."""
        return x + 1

    def test_answer():
        """Check the return value of func() for an example argument."""
        assert func(3) == 5

If we place this in a ``test.py`` file and then run::

    $ pytest test.py

The result is::

    ============================= test session starts ==============================
    python: platform darwin -- Python 3.8.3 -- pytest-3.2.0
    test object 1: /Users/username/tmp/test.py

    test.py F

    =================================== FAILURES ===================================
    _________________________________ test_answer __________________________________

        def test_answer():
    >       assert func(3) == 5
    E       assert 4 == 5
    E        +  where 4 = func(3)

    test.py:5: AssertionError
    =========================== 1 failed in 0.07 seconds ===========================

Sometimes the output from the test suite will have ``xfail`` meaning a test has passed although it has been marked as ``@pytest.mark.xfail``), or ``skipped`` meaing a test that has been skipped due to not meeting some condition (online and figure tests are the most common).

You need to use the option ``-rs`` for skipped tests and ``-rx`` for xfailed tests, respectively.
Or use ``-rxs`` for detailed information on both skipped and xfailed tests.

Where to put tests
------------------

Each package should include a suite of unit tests, covering as many of the public methods/functions as possible.
These tests should be included inside each package, e.g::

    sunpy/map/tests/

"tests" directories should contain an ``__init__.py`` file so that the tests can be imported.

Online Tests
------------

There are some tests for functions and methods in sunpy that require a working connection to the internet.
``pytest`` is configured in a way that it iterates over all tests that have been marked as ``pytest.mark.remote_data`` and checks if there is an established connection to the internet.
If there is none, the test is skipped, otherwise it is run.

Marking tests is pretty straightforward, use the decorator ``@pytest.mark.remote_data`` to mark a test function as needing an internet connection::

    @pytest.mark.remote_data
    def func(x):
        """Add one to the argument."""
        return x + 1

Tests that create files
-----------------------

Tests may often be run from directories where users do not have write permissions so tests which create files should always do so in temporary directories.
This can be done with the `pytest tmpdir function argument <https://pytest.org/en/latest/tmpdir.html>`_ or with Python's built-in `tempfile module
<https://docs.python.org/3/library/tempfile.html#module-tempfile>`_.

Tests that use test data
------------------------

We store test data in "sunpy/data/test" as long as it is less than about 100 kB.
These data should always be accessed via the :func:`sunpy.data.test.get_test_filepath` and :func:`sunpy.data.test.test_data_filenames` functions.
This way you can use them when you create a test.

You can also use our sample data but this will have to be marked as an online test (see above)::

    import sunpy.data.sample

    @pytest.mark.remote_data
    def func():
        """Returns the file path for the sample data."""
        return sunpy.data.sample.AIA_131_IMAGE

Generally we do not run the tests on our sample data, so only do this if you have a valid reason.

Figure unit tests
-----------------

.. note::
    The figure tests and the hashes they use are only checked on Linux and might be different on other platforms.
    We should suggest if you do not use a Linux, to add a fake hash to the json files and then CircleCi (ran on a PR) will tell you the real hash to use.

You can write sunpy unit tests that test the generation of Matplotlib figures by adding the decorator ``sunpy.tests.helpers.figure_test``.
Here is a simple example::

    import matplotlib.pyplot as plt
    from sunpy.tests.helpers import figure_test

    @figure_test
    def test_simple_plot():
        plt.plot([0,1])

The current figure at the end of the unit test, or an explicitly returned figure, has its hash (currently ``SHA256``) compared against an established hash collection (more on this below).
If the hashes do not match, the figure has changed, and thus the test is considered to have failed.

If you are adding a new figure test you will need to generate a new hash library::

    $ tox -e py38-figure -- --mpl-generate-hash-library=sunpy/tests/figure_hashes_mpl_332_ft_261_astropy_42.json

The filename changes if the version of astropy or Matplotlib or freetype gets updated.
So you might need to adjust this command.
For the development figure tests::

    $ tox -e py38-figure-devdeps -- --mpl-generate-hash-library=sunpy/tests/figure_hashes_mpl_dev_ft_261_astropy_dev.json

This will run the figure test suite and update the hashes stored.

If you want to check what the images look like, you can do::

    $ tox -e py38-figure -- --mpl-generate-path=baseline

The images output from the tests will be stored in a folder called ``.tmp/py38-figure/baseline`` or ``baseline`` in the sunpy folder, so you can double check the test works as you expected.

.. _doctests:

doctests
--------

Code examples in the documentation will also be run as tests and this helps to validate that the documentation is accurate and up to date.
sunpy uses the same system as Astropy, so for information on writing doctests see the astropy `documentation <https://docs.astropy.org/en/latest/development/testguide.html#writing-doctests>`_.

You do not have to do anything extra in order to run any documentation tests.
Within our ``setup.cfg`` file we have set default options for ``pytest``, such that you only need to run::

    $ pytest <file to test>

to run any documentation test.

Bugs discovered
---------------

In addition to writing unit tests new functionality, it is also a good practice to write a unit test each time a bug is found, and submit the unit test along with the fix for the problem.
This way we can ensure that the bug does not re-emerge at a later time.
.. _coding-standards:

****************
Coding Standards
****************

The purpose of the page is to describe the standards that are expected of all the code in the SunPy project.
All potential developers should read and abide by the following standards.
Code which does not follow these standards closely will not be accepted.

We try to closely follow the coding style and conventions proposed by `Astropy <https://docs.astropy.org/en/stable/development/codeguide.html#coding-style-conventions>`_.

Language Standard
=================

* All code must be compatible with Python 3.7 and later.
  Usage of ``six``, ``__future__``, and ``2to3`` is no longer acceptable.

* The new Python 3 formatting style should be used (i.e.
  ``"{0:s}".format("spam")`` instead of ``"%s" % "spam"``).

* The core package and affiliated packages should be importable with no dependencies other than components already in the sunpy core package, the `Python Standard Library <https://docs.python.org/3/library/index.html>`_, and packages already required by the sunpy core package.
  Adding dependencies to sunpy core will be considered but are highly discouraged.
  Such optional dependencies should be recorded in the ``setup.cfg`` file in the ``extras_require`` entry.

Coding Style/Conventions
========================

* The code will follow the standard `PEP8 Style Guide for Python Code <https://www.python.org/dev/peps/pep-0008/>`_.
  In particular, this includes using only 4 spaces for indentation, and never tabs.

* **Follow the existing coding style** within a file and avoid making changes that are purely stylistic.
  Please try to maintain the style when adding or modifying code.

* Following PEP8's recommendation, absolute imports are to be used in general.
  We allow relative imports within a module to avoid circular import chains.

* The ``import numpy as np``, ``import matplotlib as mpl``, and ``import matplotlib.pyplot as plt`` naming conventions should be used wherever relevant.
  ``from packagename import *`` should never be used (except in ``__init__.py``)

* Classes should either use direct variable access, or Python's property mechanism for setting object instance variables.

* Classes should use the builtin `super` function when making calls to methods in their super-class(es) unless there are specific reasons not to.
  `super` should be used consistently in all subclasses since it does not work otherwise.

* Multiple inheritance should be avoided in general without good reason.

* ``__init__.py`` files for modules should not contain any significant implementation code. ``__init__.py`` can contain docstrings and code for organizing the module layout.


Private code
============

It is often useful to designate code as private, which means it is not part of the user facing API, only used internally by sunpy, and can be modified without a deprecation period.
Any classes, functions, or variables that are private should either:

- Have an underscore as the first character of their name, e.g., ``_my_private_function``.
- If you want to do that to entire set of functions in a file, name the file with a underscore as the first character, e.g., ``_my_private_file.py``.

If these might be useful for other packages within the sunpy ecosphere, they should be made public.

Utilities in sunpy
==================

Within ``sunpy``, it might be useful to have a set of utility classes or functions that are used by internally to help with certain tasks or to provide a certain level of abstraction.
These should be placed either:

- ``sunpy.{subpackage}.utils.py``, if it is only used within that sub-package.
- ``sunpy.util`` if it is used across multiple sub-packages.

These can be private (see section above) or public.
The decision is up to the developer, but if these might be useful for other packages within the sunpy ecosphere, they should be made public.

Formatting
==========

We enforce a minimum level of code style with our continuous intergration (the name is ``sunpy.sunpy (python_codestyle [linux]``).
This runs a tool called `pre-commit <https://pre-commit.com/>`__.

The settings and tools we use for the pre-commit can be found in the file :file:`.pre-commit-config.yaml` at the root of the sunpy git repository.
Some of the checks are:
* Checks (but doesn't fix) various PEP8 issues with flake8.
* Sort all imports in any Python files with isort.
* Remove any unused variables or imports with autoflake.

We suggest you use "tox" (which is used to run the sunpy test suite) to run these tools without having to setup anything within your own Python virtual environment::

    $ tox -e codestyle

What you will see is this output (heavily condensed):

.. code-block:: bash

    codestyle create: /home/<USER>/GitHub/sunpy/.tox/codestyle
    codestyle run-test: commands[0] | pre-commit install-hooks
    codestyle run-test: commands[1] | pre-commit run --verbose --all-files --show-diff-on-failure
    flake8...................................................................Passed
    - hook id: flake8
    - duration: 1.35s

    0

    Check for case conflicts.................................................Passed
    - hook id: check-case-conflict
    - duration: 0.08s
    Trim Trailing Whitespace.................................................Failed
    - hook id: trailing-whitespace
    - duration: 0.08s
    - exit code: 1
    - files were modified by this hook

    Fixing docs/dev_guide/code_standards.rst

    pre-commit hook(s) made changes.
    If you are seeing this message in CI, reproduce locally with: `pre-commit run --all-files`.
    To run `pre-commit` as part of git workflow, use `pre-commit install`.
    All changes made by hooks:
    diff --git a/docs/dev_guide/code_standards.rst b/docs/dev_guide/code_standards.rst
    index bed700d90..c6b5df977 100644
    --- a/docs/dev_guide/code_standards.rst
    +++ b/docs/dev_guide/code_standards.rst
    @@ -59,6 +59,8 @@ Instead of installing this, you can use "tox" (which is used to run the sunpy te

        $ tox -e codestyle

    +What you will see
    +
    If you want to setup the pre-commit locally, you can do the following::

        $ pip install pre-commit
    diff --git a/docs/dev_guide/documentation.rst b/docs/dev_guide/documentation.rst
    index 5cd914047..b1017f77a 100644
    --- a/docs/dev_guide/documentation.rst
    +++ b/docs/dev_guide/documentation.rst
    @@ -39,9 +39,9 @@ If there are multiple code elements with the same name (e.g. ``peek()`` is a met

    .. code-block:: rst

    -    `GenericMap.peek` or `CompositeMap.peek`
    +    `.GenericMap.peek` or `.CompositeMap.peek`

    -These will show up as `GenericMap.peek` or `CompositeMap.peek`.
    +These will show up as `.GenericMap.peek` or `.CompositeMap.peek`.
    To still show only the last segment you can add a tilde as prefix:

    ERROR: InvocationError for command /home/nabil/GitHub/sunpy/.tox/codestyle/bin/pre-commit run --verbose --all-files --show-diff-on-failure (exited with code 1)
    ___________________________________________________________________________________________ summary ___________________________________________________________________________________________
    ERROR:   codestyle: commands failed

This will inform you of what checks failed and why, and what changes (if any) the command has made to your code.

If you want to setup the pre-commit locally, you can do the following::

    $ pip install pre-commit

Now you can do::

    $ pre-commit run --all-files

which will run the tools on all files in the sunpy git repository.
The pre-commit tools can change some of the files, but in other cases it will report problems that require manual correction.
If the pre-commit tool changes any files, they will show up as new changes that will need to be committed.

Automate
--------

Instead of running the pre-commit command each time you can install the git hook::

    $ pre-commit install

which installs a command to :file:`.git/hooks/pre-commit` which will run these tools at the time you do ``git commit`` and means you don't have to run the first command each time.
We only suggest doing the install step if you are comfortable with git and the pre-commit tool.

Documentation and Testing
=========================

* American English is the default language for all documentation strings and inline commands.
  Variables names should also be based on English words.

* Documentation strings must be present for all public classes/methods/functions, and must follow the form outlined in the :ref:`docs_guidelines` page.
  Additionally, examples or tutorials in the package documentation are strongly recommended.

* Write usage examples in the docstrings of all classes and functions whenever possible.
  These examples should be short and simple to reproduce–users should be able to copy them verbatim and run them.
  These examples should, whenever possible, be in the :ref:`doctests` format and will be executed as part of the test suite.

* Unit tests should be provided for as many public methods and functions as possible, and should adhere to the standards set in the :ref:`testing` document.

Data and Configuration
======================

* We store test data in ``sunpy/data/test`` as long as it is less than about 100 kB.
  These data should always be accessed via the :func:`sunpy.data.test.get_test_filepath` and :func:`sunpy.data.test.test_data_filenames` functions.

* We store data used for examples in the `sample-data repository <https://github.com/sunpy/sample-data>`_.
  This data should not be used for unit tests but can be within our documentation.

* All persistent configuration should use the :ref:`config` mechanism.
  Such configuration items should be placed at the top of the module or package that makes use of them, and supply a description sufficient for users to understand what the setting
  changes.

Standard output, warnings, and errors
=====================================

The built-in ``print(...)`` function should only be used for output that is explicitly requested by the user, for example ``print_header(...)`` or ``list_catalogs(...)``.
Any other standard output, warnings, and errors should follow these rules:

* For errors/exceptions, one should always use ``raise`` with one of the built-in exception classes, or a custom exception class.
  The nondescript ``Exception`` class should be avoided as much as possible, in favor of more specific exceptions (`IOError`, `ValueError`, etc.).

* For warnings, one should always use the functions in ``sunpy.util.exceptions`` and *not* `warnings.warn`. This ensures we are always raising a sunpy specific warning type.

Including C Code
================

* C extensions are only allowed when they provide a significant performance enhancement over pure Python, or a robust C library already exists to provided the needed functionality.

* The use of `Cython`_ is strongly recommended for C extensions.

* If a C extension has a dependency on an external C library, the source code for the library should be bundled with sunpy, provided the license for the C library is compatible with the sunpy license.
  Additionally, the package must be compatible with using a system-installed library in place of the library included in sunpy.

* In cases where C extensions are needed but `Cython`_ cannot be used, the `PEP 7 Style Guide for C Code <https://www.python.org/dev/peps/pep-0007/>`_ is recommended.

* C extensions (`Cython`_ or otherwise) should provide the necessary information for building the extension.

.. _Cython: https://cython.org/
.. _remote_data:

*******************
Remote Data Manager
*******************

There are functions which require files to be downloaded from the internet before they are executed.
If the files change on the remote server or if they are tampered with after the file is downloaded, the function may return erroneous results.

This Remote Data Manager provides developers with a way to download the files easily and ensure that the files are intact when they are used.
If the file is changed on the remote server, the data manager will raise an error so that the user is made aware that the intended file is not available.
If the file has changed on disk, data manager will redownload the file after warning the user.

Also in `sunpy.util`, there is a `~sunpy.util.hash_file` function.

Usage
=====

The manager has to be imported before it is used::

    from sunpy.data import manager

`~sunpy.data.data_manager.manager.DataManager.require` is a decorator which is used to inform the data manager that the function requires a file for its execution.
`~sunpy.data.data_manager.manager.DataManager.get` function is used to access the file inside the function.
Suppose a function requires a file with url 'http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt' which has a SHA256 hash of '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11'.
The following example will show how this function can be implemented::

    @manager.require('test_file',
                    ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
                    '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11')
    def test_function():
        file_path = manager.get('test_file')
        # do something with file
        return file_path

Cache
-----

Remote files sometimes have to be cached.
This is to ensure that the files are not redownloaded frequently, thus saving users' disk space and well as internet bandwidth.
The cache has an expiry time which is set by the ``sunpyrc`` config file.
To change this please see :ref:`customizing-sunpy`.

Usage
=====

The following example shows how cache can be used::

    from sunpy.data import cache

    def test_function():
        file_path = cache.download('http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt')
        return file_path


    test_function()  # The file will be downloaded with this call
    test_function()  # subsequent calls won't redownload the file


Testing
=======

A pytest fixture is provided for ease of mocking network requests when using cache.
The following example demonstrates the usage of the fixture.::

    @pytest.fixture()
    def local_cache(sunpy_cache):
        sunpy_cache = sunpy_cache('sunpy.test_module.cache')
        sunpy_cache.add('http://example.com/test_file',
                        'test_data_path')

The above snippet creates a pytest fixture called ``local_cache``. This fixture can be used in wherever the files have to be mocked.
An example is given below.::

    def test_test_function(local_cache):
        # inside this function the mocked cache is used

        # test_function uses 'http://example.com/test_file'
        assert test_function() == True
.. _config:

***************
Global Settings
***************

sunpy makes use of a settings file (:file:`sunpyrc`).
This file contains a number of global settings such as where files should be downloaded by default or the default format for displaying times.
When developing new functionality check this file and make use of the default values if appropriate or, if needed, define a new value.
More information can be found in :ref:`customizing-sunpy`.
.. _newcomers:

****************
Newcomers' Guide
****************

Welcome to the SunPy newcomers' guide.
If you have come across this page, you just might be new to the SunPy project.
We aim to be a comprehensive Python package that allows solar physicists to deep dive in the vast amount of solar data available.

Firstly, we want to thank you for your interest in contributing to SunPy!
SunPy is an open project that encourages everyone to contribute in any way possible.

The people who help develop or contribute to SunPy are varied in ability and experience with the vast majority being volunteers who dedicate time each week.
We pride ourselves on being a welcoming community and we would love to have you become a part of our community.

Although this document mainly focuses on how to make contributions to the core sunpy libraries code and documentation, there are other ways to get involved with the SunPy community.

If you have any questions, comments or just want to say hello, we have online chat on `Matrix`_ which requires no registration or a `Google Group`_ which you message.

.. _Matrix: https://openastronomy.element.io/#/room/#sunpy:openastronomy.org
.. _Google Group: https://groups.google.com/forum/#!forum/sunpy

How to Contribute to sunpy
==========================

Not Code
--------

A common misconception (which applies to any package) is that all we really want is some Python code, in fact, we do not require only code contributions!
If you do not have the time or the desire to code, we have severals of areas where we can use help.

Reporting Issues
^^^^^^^^^^^^^^^^

If you use sunpy and stumble upon a problem, the best way to report it is by opening an `issue`_ on our GitHub issue tracker.
This way we can help you work around the problem and hopefully fix the problem!

You will need to sign into `GitHub`_ to report an issue and if you are not already a member of Github, you will have to join.
Joining GitHub will make it easier to report and track issues in the future.

If you do not want to join Github, then another way to report your issue is email the SunPy developers list `sunpy-dev@googlegroups.com`_.

When reporting an issue, please try to provide a short description of the issue with a small code sample, this way we can attempt to reproduce the error.
Also provide any error output generated when you encountered the issue, we can use this information to debug the issue.
For a good example of how to do this see issue `#2879`_.

If there is functionality that is not currently available in sunpy you can make a feature request.
Please write as much information as possible regarding the feature you would like to see in sunpy.

When you go to open either an issue or a feature request, there is a GitHub template that will guide you on the type of information we are seeking.
Please be sure to read over the comments in the GitHub text box.

.. _issue: https://github.com/sunpy/sunpy/issues
.. _sunpy-dev@googlegroups.com: https://groups.google.com/forum/#!forum/sunpy-dev
.. _#2879: https://github.com/sunpy/sunpy/issues/2879

Documentation
^^^^^^^^^^^^^

sunpy has `online documentation`_ and we try to make sure its as comprehensive as possible.
This documentation contains the API of sunpy but also a user guide, an example gallery and developer documents.

However, documentation for any project is a living document.
It is never complete and there are always areas that could be expanded upon or could do with proof reading to check if the text is easy to follow and understandable.
If parts are confusing or difficult to follow, we would love suggestions or improvements!

.. _online documentation: https://docs.sunpy.org/en/latest/index.html

Reviewing a Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^

We at any one time have a variety of `pull requests`_ open and getting reviews is important.
Generally the more people that can look over a pull request the better it will turn out and we encourage everyone to do so.

.. _pull requests: https://github.com/sunpy/sunpy/pulls

Code
----

If you would prefer to code instead, the best way to start is to work on an exisiting and known `issues`_.
We have several repositories you can investigate.
The main one is the sunpy repository with where all the known `issues`_ with sunpy are detailed.
Each issue should have a series of labels that provide information about the nature of the issue.
If you find an issue you'd like to work on, please make sure to add a comment to let people know that you are working on it! This will make it less likely that effort is duplicated.

.. note::

    sunpy is Free and open-source software (FOSS), under the BSD-2 license. By contributing you are stating that you have the right to and agree to have your work distributed under the terms of this license.

    This applies to everyone who wants to contribute during work time no matter who their employer is.
    You should start by checking if there is a Open Source Software Policy (e.g., `Standford's policy <https://otl.stanford.edu/open-source-stanford>`__) for your work place.
    If not, `OSS-Watch <http://oss-watch.ac.uk/resources/contributing>`__ summaries what you will need to check and who to ask, however this resource is aimed at a UK readers.
    As an example, `Standford's guidance <https://otl.stanford.edu/sites/g/files/sbiybj10286/f/otlcopyrightguide.pdf>`__ allows someone to contribute and open source their code.
    If you are unsure if your university or institution allows you to contribute under the BSD-2 license, you should contact the relevant department or administrator that deals with copyright at your institution.

If you are unsure where to start we suggest the `Good First Issue label`_.
These are issues that have been deemed a good way to be eased into sunpy and are achievable with little understanding of the sunpy codebase.
Please be aware that this does not mean the issue is "easy", just that you do not need to be aware of the underlying structure of sunpy.

We also tag issues for specific events such as  `Hacktoberfest`_ under the `Hacktoberfest label`_.
The scope of the issues should be appropriate for that specific event.
We do particpate in several other events but right now we do not have dedicated labels.
So please use the above labels for starting issues!

In addition, we have several other repositories that have open issues and you might find these more interesting than the main repository.

Python:

* `ndcube <https://github.com/sunpy/ndcube>`_
* `drms <https://github.com/sunpy/drms>`_
* `radiospectra <https://github.com/sunpy/radiospectra>`_
* `ablog <https://github.com/sunpy/ablog>`_
* `irispy <https://github.com/sunpy/irispy>`_
* `sunkit-image <https://github.com/sunpy/sunkit-image>`_

CSS/HTML:

* `sunpy-sphinx-theme <https://github.com/sunpy/sunpy-sphinx-theme>`_
* `sunpy.org <https://github.com/sunpy/sunpy.org>`_

.. _issues: https://github.com/sunpy/sunpy/issues
.. _Good First Issue label: https://github.com/sunpy/sunpy/issues?utf8=%E2%9C%93&q=is%3Aissue+is%3Aopen+label%3A%22Good+First+Issue%22
.. _Hacktoberfest: https://hacktoberfest.digitalocean.com/
.. _Hacktoberfest label: https://github.com/sunpy/sunpy/issues?q=is%3Aissue+is%3Aopen+label%3AHacktoberfest

If you already have code that you've developed or already have your own idea of code you'd like to work on please first have a look at the issue list to see if any existing issues are related to your idea.
If you find nothing then create your own issue to stimulate the addition of your code and make sure to let people know about it chat room or by email.
Creating an issue creates a permanent record.
It may be possible your idea may be outside of the scope of the repository you are trying to contribute to and the issue comments are a great place to have that discussion where potential future contributors can also see.

Setting up a development environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The instructions in this following section are based upon three resources:

* `Astropy Dev Workflow <https://docs.astropy.org/en/latest/development/workflow/development_workflow.html>`_
* `Astropy Dev environment <https://docs.astropy.org/en/latest/development/workflow/get_devel_version.html#get-devel>`_
* `Astropy Pull Request Example <https://docs.astropy.org/en/latest/development/workflow/git_edit_workflow_examples.html#astropy-fix-example>`_

**We strongly recommend that you read these links.**
These links are more in-depth than this guide but you will need to replace ``astropy`` with ``sunpy``.

In order to start coding you will need a local Python environment and we would recommend using `Anaconda`_ or `miniconda`_ (shortened to conda from here on).
This method will bypass your operating system Python packages and makes the entire process easier.

The first step is to install the version of conda that corresponds to your operating system and `instructions are here`_.
Next we will want to setup the conda environment and we will need to add the `conda-forge`_ channel as a prerequisite:

.. code:: bash

    $ conda config --add channels conda-forge
    $ conda create -n sunpy-dev pip
    $ conda activate sunpy-dev

This will create a new conda environment called "sunpy-dev" and install the latest version of pip from the conda-forge channel.
The next step is get a developement version of sunpy.
This will require that `git`_ be installed.
If you have a `GitHub`_ account, we suggest that you `fork`_ the `sunpy repository`_ (the fork button is to the top right) and **use that url for the clone step**.
This will make submitting changes easier in the long term for you:

.. warning::

    Do not clone the sunpy repository into ``$HOME/sunpy``. Depending on the operating system this location is used to store downloaded data files.
    This will cause conflicts later on, so the last argument (``sunpy-git``) on the ``git clone`` line will become the local folder name of the cloned repository.
    Otherwise you are free to clone to any other location.

.. code:: bash

    $ git clone https://github.com/<username>/sunpy.git sunpy-git
    $ cd sunpy-git
    $ pip install -e .[dev]

.. note::
    If this does not work, it could be due to a missing C compiler (e.g., ``gcc`` or ``clang``) that is required to build sunpy at install.
    Getting the compiler either from your system package manager, XCode or Anaconda should address this.

Now you have the latest version of sunpy installed and are ready to work on it using your favorite editor!
Ideally, when you start making changes you want to create a git branch:

.. code:: bash

    $ git checkout -b my_fix

You can change ``my_fix`` to anything you prefer.
If you get stuck or want help, just `ask here`_!

.. _Anaconda: https://www.anaconda.com/distribution/
.. _miniconda: https://conda.io/en/latest/miniconda.html
.. _instructions are here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html#installation
.. _conda-forge: https://conda-forge.org/
.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _GitHub: https://github.com/
.. _fork: https://guides.github.com/activities/forking/
.. _sunpy repository: https://github.com/sunpy/sunpy
.. _ask here: https://openastronomy.element.io/#/room/#sunpy:openastronomy.org

Checking the code you have written
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that you have written some code to address an issue.
You will need to check two things:

1. The changes you have made are correct, i.e., it fixes a bug or the feature works.
   This requires you to run the code either manually or by writting/running a test function.
   `pytest`_ is the framework we use for this.

2. The changes you have made follow the correct coding style.
   We follow the `PEP8`_ style for all Python code and depending on your setup, you can use a `linter program <https://realpython.com/python-code-quality/#how-to-improve-python-code-quality>`_ to check your code.
   For documentation, we follow the `numpydoc style <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_.

We provide more more detail about our :ref:`test suite and how to write tests <testing>`, and how to :ref:`create and style documentation <docs_guidelines>`.

.. _pytest: https://docs.pytest.org/en/latest/

Send it back to us
^^^^^^^^^^^^^^^^^^
Once you have some changes you would like to submit, you will need to commit the changes.
This is a three stage process:

1. Use ``git status`` to see that the only changes locally are the right ones.
2. Use ``git add <path to file>`` to add the changes to ``git``.
3. Use ``git commit -m <message>`` to label those changes.
4. Use ``git push`` to update your fork (copy) of sunpy on GitHub.

Here you replace ``<message>`` with some text of the work you have done.
We strongly recommend having a good commit message and this `commit guide`_ is worth reading.

Next step is to open a pull request on GitHub.
If you are new to pull requests, here is the `GitHub guide`_ that is a detailed walkthrough.
Go to the "pull requests" tab on **your fork** and pressing the large green "New pull request" button.
Now on the right side from the box marked "compare" you can select your branch.
Do one final check to make sure the code changes look correct and then press the green "Create pull request" button.

When you open your pull request, we have a GitHub template that will guide you on what to write in the message box.
Please fill this in and title the pull request.
Now the final step is to press the green "Create pull request" button.

As soon as you do this, you will be greeted by a message from the "sunpy bot" as well as several continuous integration checks.
These are explained on our :ref:`Pull Request Review <pr_review>` page.
But what is important to know is that these run a series of tests to make sure that the changes do not cause any new errors.
We strongly recommend that any code changes you have had, follow the `PEP8`_ style and that you have ran the code locally to make sure any changes do not break any existing code.
We provide an overview on how to run the test suite :ref:`here <testing>`.
Now we (the sunpy community) can review the code and offer suggestions and once we are happy, we can merge in the pull request.

If you do not have time to finish what you started on or ran out of time during a sprint and do not want to submit a pull request, you can create a git patch instead:

.. code:: bash

    $ git format-patch main --stdout > my_fix.patch

You can rename ``my_fix`` to something more relevant.
This way, you still get acknowledged for the work you have achieved.
Now you can email this patch to the  `Google Group`_ .

Just remember, if you have any problems get in touch!

.. _commit guide: https://chris.beams.io/posts/git-commit/
.. _GitHub guide: https://guides.github.com/activities/hello-world/
.. _PEP8: https://realpython.com/python-pep8/
.. _Google Group: https://groups.google.com/forum/#!forum/sunpy

Summer of Code(s)
^^^^^^^^^^^^^^^^^

If you are interested in a "Summer of Code" project with sunpy, we have information on our `wiki`_ which has guidelines, advice, application templates and more!
Our projects are located on our umbrella's organization website, `OpenAstronomy`_.

.. _wiki: https://github.com/sunpy/sunpy/wiki#summer-of-codes
.. _OpenAstronomy: https://openastronomy.org/gsoc/
.. _docs_guidelines:

*************************
SunPy Documentation Rules
*************************

Overview
========

All code must be documented and we follow these style conventions described here:

* `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_

We recommend familiarizing yourself with this style.

Referring to other code
-----------------------

To link to other methods, classes, or modules in sunpy you have to use backticks, for example:

.. code-block:: rst

    `sunpy.io.read_file`

generates a link like this: `sunpy.io.read_file`.

We use the sphinx setting ``default_role = 'obj'`` so that you do not nor **SHOULD NOT** use the ``:class:`` qualifier, but ``:func:``, ``:meth:`` are different (more on this below).

Often, you don't want to show the full package and module name.
As long as the target is unambiguous you can simply leave them out:

.. code-block:: rst

    `.read_file`

and the link still works: `.read_file`.

If there are multiple code elements with the same name (e.g. ``peek()`` is a method in multiple classes), you'll have to extend the definition:

.. code-block:: rst

    `.GenericMap.peek` or `.CompositeMap.peek`

These will show up as `.GenericMap.peek` or `.CompositeMap.peek`.
To still show only the last segment you can add a tilde as prefix:

.. code-block:: rst

    `~.GenericMap.peek` or `~.CompositeMap.peek`

will render as `~.GenericMap.peek` or `~.CompositeMap.peek`.

Other packages can also be linked via
`intersphinx <http://www.sphinx-doc.org/en/master/ext/intersphinx.html>`_:

.. code-block:: rst

    `numpy.mean`

will return this link: `numpy.mean`.
This works for Python, Numpy and Astropy (full list is in :file:`docs/conf.py`).

With Sphinx, if you use ``:func:`` or ``:meth:``, it will add closing brackets to the link.
If you get the wrong pre-qualifier, it will break the link, so we suggest that you double check if what you are linking is a method or a function.

.. code-block:: rst

    :class:`numpy.mean()`
    :meth:`numpy.mean()`
    :func:`numpy.mean()`

will return two broken links ("class" and "meth") but "func" will work.

sunpy-Specific Rules
--------------------

* For **all** RST files, we enforce a one sentence per line rule and ignore the line length.
* Standards on docstring length and style are enforced using `docformatter <https://pypi.org/project/docformatter/>`__:

.. code-block:: bash

    $ docformatter -r -i  --pre-summary-newline --make-summary-multi-line

.. _Docs Guidelines for Data Sources:

Documenting Data Sources
----------------------------

Subclasses of `~sunpy.map.GenericMap` or `~sunpy.timeseries.TimeSeries` must provide a detailed docstring providing an overview of the data source that the object represents.
In order to maintain consistency and completeness, the following information must be provided by a data source docstring, if available, and preferably in the following order:

* the name of the mission and instrument and the institution that built it
* short description of the instrument (e.g. Cassegrain reflector, Wolter-1 grazing incidence x-ray, coronagraph) including the type of detector
* description of the platform (e.g. satellite in 28 deg inclined orbit, a telescope on the summit of Mauna Kea in Hawaii)
* description of the primary purpose or science goals of the instrument.
* list of all wavelength(s) or passbands in appropriate units
* description of the emission processes which dominate in those passbands
* appropriate measurement properties such as field of view, angular resolution, time resolution
* description of the operational concept (e.g. operates 24/7, observes from 7 am to 5 pm UT) including mention of unusual operations scenarios (e.g. calibration seasons, eclipse seasons)
* the start and end of the data set

In addition, a reference section must be provided with links to the following resources, if available,

* the mission web page
* the instrument web page
* relevant wikipedia page(s)
* relevant user guide(s)
* the mission paper and instrument paper
* information to interpret metadata keywords such as FITS header reference
* the data archive

An example docstring can be found in the :ref:`Writing a new Instrument Map Class guide <new_maps_ts_etc>`.

Sphinx
======

All of the sunpy documentation (like this page) is built by `Sphinx <https://www.sphinx-doc.org/en/stable/>`_, which is a tool especially well-suited for documenting Python projects.
Sphinx works by parsing files written using a `a Mediawiki-like syntax <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ called `reStructuredText <http://docutils.sourceforge.net/rst.html>`_.
In addition to parsing static files of reStructuredText, Sphinx can also be told to parse code comments.
In fact, in addition to what you are reading right now, the `Python documentation <https://www.python.org/doc/>`_ was also created using Sphinx.

Usage
-----

All of the sunpy documentation is contained in the "docs" folder and code documentation strings.
The examples from the example gallery can be found in the "examples" folder.

In the root directory run::

    $ tox -e build_docs

This will generate HTML documentation for sunpy in the "docs/_build/html" directory.
You can open the "index.html" file to browse the final product.
The gallery examples are located under "docs/_build/html/generated/gallery".
Sphinx builds documentation iteratively, only adding things that have changed.

If you want to build the documentation without building the gallery, i.e. to reduce build times while working on other sections of the documentation you can run::

    $ tox -e build_docs -- -D plot_gallery=False

If you'd like to start from scratch (i.e., remove the tox cache) then run::

    $ tox -e build_docs -- -aE

To build the documentation in your current python environment you must have all the dependencies specified in ``setup.cfg`` installed (``pip install -e .[docs]``).
Then change to the :file:`docs/` directory and run::

    $ make html

For more information on how to use Sphinx, consult the `Sphinx documentation <http://www.sphinx-doc.org/en/stable/contents.html>`_.

Special Sphinx directives
-------------------------

``minigallery`` directive
^^^^^^^^^^^^^^^^^^^^^^^^^

Sphinx will automatically record which functions, classes, etc. are used in each gallery example.
In the documentation, you can insert a mini-gallery of the subset of the gallery examples that uses a particular function, class, etc.
For example, the following RST block::

    .. minigallery:: sunpy.coordinates.RotatedSunFrame

produces this mini-gallery:

.. minigallery:: sunpy.coordinates.RotatedSunFrame

If you want to specify more than one object, separate them by spaces.
This is particularly useful if you need to cover multiple namespaces in which an object may be accessed, e.g.::

    .. minigallery:: sunpy.coordinates.RotatedSunFrame sunpy.coordinates.metaframes.RotatedSunFrame

``generate`` directive
^^^^^^^^^^^^^^^^^^^^^^

In rare circumstances, one may want to insert "raw" HTML directly into the pages written by Sphinx.
For HTML that is statically available (i.e., already written in some form), one can use the `"raw" directive <https://docutils.sourceforge.io/docs/ref/rst/directives.html#raw-data-pass-through>`__.
For HTML that is generated by Python code, sunpy provides the custom directive ``generate``.
Here's an example RST block::

    .. generate:: html
        :html_border:

        import os
        from sunpy.data.sample import file_dict
        print("<table>")
        for key, value in file_dict.items():
            print(f"<tr><th>{key}</th><td>{os.path.basename(value)}</td></tr>")
        print("</table>")

to insert the following HTML table:

.. generate:: html
    :html_border:

    import os
    from sunpy.data.sample import file_dict
    print("<table>")
    for key, value in file_dict.items():
        print(f"<tr><th>{key}</th><td>{os.path.basename(value)}</td></tr>")
    print("</table>")

Troubleshooting
----------------

Sphinx can be very particular about formatting, and the warnings and errors aren't always obvious.

Below are some commonly-encountered warning/error messages along with a human-readable translation:

**WARNING: Duplicate explicit target name: "xxx".**

If you reference the same URL, etc more than once in the same document sphinx will complain.
To avoid, use double-underscores instead of single ones after the URL.

**ERROR: Malformed table. Column span alignment problem at line offset n**

Make sure there is a space before and after each colon in your class and
function docs (e.g. attribute : type, instead of attribute: type).
Also, for some sections (e.g. Attributes) numpydoc seems to complain when a description spans more than one line, particularly if it is the first attribute listed.

**WARNING: Block quote ends without a blank line; unexpected unindent.**

Lists should be indented one level from their parents.

**ERROR: Unknown target name: "xxx"**

In addition to legitimate errors of this type, this error will also occur when variables have a trailing underscore, e.g., ``xxx_``.

**WARNING: Explicit markup ends without a blank line; unexpected unindent.**

This usually occurs when the text following a directive is wrapped to the next line without properly indenting a multi-line text block.

**WARNING: toctree references unknown document '...'** / **WARNING: toctree contains reference to nonexisting document**

This pair of errors is due to the way numpydoc scrapes class members.
.. _dev_logger:

*********************************
Logging, Warnings, and Exceptions
*********************************

Overview
========

sunpy makes use of a logging system to deal with messages (see :ref:`logger`). This provides the users and
developers the ability to decide which messages to show, to capture them, and to optionally also send
them to a file. The logger will log all messages issued directly to it but also warnings issued
through `warnings.warn` as well as exceptions.

The logger is configured as soon as sunpy is imported. You can access it
by importing it explicitly::

    from sunpy import log

Messages can be issued directly to it with the following levels and in the following way::

    log.debug("Detailed information, typically of interest only when diagnosing problems.")

    log.info("A message conveying information about the current task, and confirming that
              things are working as expected.")

    log.warning("An indication that something unexpected happened, or indicative of
                 some problem in the near future (e.g. disk space low).

                 The software is still working as expected.")
    log.error("Due to a more serious problem, the software has not been able to
               perform some function but the task is still continuing.")

    log.critical("A serious error, indicating that the program itself may be unable to
                  continue running. A real error may soon by issued and the task will fail.")

The difference between logging a warning/error/critical compared to issuing a Python warning or raising
an exception are subtle but important.

Use Python `warnings.warn` in library code if the issue is avoidable and the user code should be
modified to eliminate the warning.

Use ``log.warning()`` if there is likely nothing the user can do about the situation, but the event
should still be noted. An example of this might be if the input data are all zeros. This may be unavoidable or
even by design but you may want to let the user know.

True exceptions (not ``log.error()``) should be raised only when there is no way for the function to proceed.

Regardless of the type (log message or warning or exception) messages should be one or two complete sentences
that fully describe the issue and end with a period.

Issuing Warnings
================

sunpy warnings are provided by the `sunpy.util` module. The primary warning which
should be used is `sunpy.util.exceptions.SunpyUserWarning`. For deprecation use `sunpy.util.exceptions.SunpyDeprecationWarning` or
`sunpy.util.exceptions.SunpyPendingDeprecationWarning`.

These three warning types have corresponding functions to raise them::

    >>> from sunpy.util.exceptions import warn_user
    >>> from sunpy.util.exceptions import warn_deprecated
    >>> from sunpy.util.exceptions import warn_metadata

These warning functions must be used to interact correctly with the logging system.
A warning can be issued in the following way::

    >>> from sunpy.util.exceptions import warn_user
    >>> warn_user("You have been warned about something you did not do correctly.")  # doctest: +IGNORE_WARNINGS

See the section above for a discussion about the distinction between ``log.warn()`` and raising a warning.

Raising Exceptions
==================

Raising errors causes the program to halt. Likely the primary error that a sunpy developer will
want to use is

* ValueError: should be raised if the input is not what was expected and could not be used. Functions should not return anything (like None) in that case or issue a warning.

Exceptions are raised simply with::

    >>> raise ValueError("The following error occurred.")  #doctest: +SKIP

For more information on exceptions see the :ref:`python:bltin-exceptions` documentation.
.. _pr_review:

******************************
Pull Requests and GitHub Teams
******************************

This document describes the standards required for a pull request to sunpy and an explanation of our automated tests.

Each pull request **must** meet the following criteria before it is considered for merge:

* The code must be PEP 8 compliant and meet the ill-defined sunpy quality standards.
  We have these in the :ref:`coding-standards` page.

* The PR must contain a changelog entry if it changes the behavior of any code.

* The test coverage should not decrease, and for new features should be at or very close to 100%.

* All code must be properly documented.
  Each function and each class must have an associated documentation string in the correct format.

Review Process
==============

Before the "merge" button is clicked the following criteria must be met:

* All the continuous integration must pass unless there is a known issue.

* At least two members (not the author of the PR) of the "sunpy-developers" group must have approved the PR, one should, ideally, be a relevant subpackage maintainer.

* All comments posted on the thread must be resolved.

It is important that approval for merging the PR is done on the comment thread, as this becomes part of the "permanent record", this includes in during community meetings or in chat.

Continuous Integration
======================

Currently we have a variety of services that respond or activate on an opened pull request.

Comments from bots:

* `pep8speaks <https://github.com/OrkoHunter/pep8speaks>`_: Performs a PEP8 check on any submitted code. This is updated as the code changes.

Checks that appear at the bottom of a pull request:

.. image:: images/checks_pr.png
   :width: 600
   :alt: PR checks

or at the top under the "Checks" tab:

.. image:: images/checks.png
   :width: 600
   :alt: PR checks tab

* `figure-tests (CircleCi) <https://circleci.com/gh/sunpy/sunpy/>`_: Runs two figure tests environments ("py39-figure", "py38-figure-devdeps").

* figure_report/figure_tests (Giles): Show the final results of the figure tests.

* figure_report_devdeps (Giles): Show the final results of the figure tests using development packages.

* changelog: absent | found (Giles): If a changelog is needed, this will check and will pass if a changelog with the correct number is found.

* `docs/readthedocs.org:sunpy (Read the Docs) <https://readthedocs.org/projects/sunpy/>`_: This builds our documentation.
  This primary check is to ensure the documentation has rendered correctly.
  Warnings are not checked on this build but under Azure Pipelines (see below).

* `sunpy.sunpy (Azure Pipelines) <https://dev.azure.com/sunpy/sunpy/_build>`_: Runs our test suite on all three operating systems.
  There are 11 separate checks for this.

* `codecov/patch (CodeCov) <https://codecov.io/gh/sunpy/sunpy/>`_: Checks how many lines of the code lack test coverage for the submitted code in the pull request.

* `codecov/project (CodeCov) <https://codecov.io/gh/sunpy/sunpy/>`_: Checks how many lines of the code lack test coverage in sunpy overall.

* `pre-commit - pr <https://pre-commit.ci>`__: Checks the code style checks have passed. This CI will automatically fix style issues by commenting ``pre-commit.ci autofix`` on its own line in a comment on the PR.

It is common to see some of these checks fail.
This can be happen due to a change that has broken a test (should be fixed) or a remote server has failed (might have to wait for it to come back).
Therefore it is important to check why a task failed and if has a pre-existing issue, it can be safe to ignore a failing check on that pull request.
However, you should try to ensure that as many checks pass before merging.

Understanding Azure Pipelines
-----------------------------

The vast majority of our tests are run on Azure Pipelines and this means you might have to navigate to the results if you want to check why the tests failed.
The tests for Azure Pipelines are split into two phases to reduce the number of builds running at one time.
If your PR fails the the linux offline tests or the "style_check" check the second stage tests will not run.

The Azure checks on GitHub manifest:

.. image:: images/azure_check_pr.png
   :width: 600
   :alt: PR checks tab

This is the main form. There will be one check per Azure job ran, and a summary one called "sunpy.sunpy".
The details text will redirect you to the "Checks" tab.

Doing so will show:

.. image:: images/azure_summary_check.png
   :width: 600
   :alt: Summary of Azure outputs on Checks tab

You get some statistics that you don't need to worry about and then a series of boxes under the "ANNOTATIONS" heading.
Unfortunately, when a Azure step fails you sometimes will get "Bash exited with code '1'." which means you have to go to the page to see what happened.
If the failure is due to a test, you will get a selection of test outputs under this heading.

On the left you should see the entire list of Azure checks.
You can go to a failing check and you will see:

.. image:: images/azure_goto.png
   :width: 600
   :alt: Go to Azure Pipelines

which will take you to the Azure Pipelines website.
This will load up the following:

.. image:: images/azure_steps_in_job.png
   :width: 600
   :alt: Build steps in Azure

Here you can see each step that is undertaken during a job on Azure.
Normally the "Running tox" should be red if the tests have failed.
You will need to click on this so it will load the output from the test suite.

Our test suite is very verbose, so there will be a lot of text outputted.
The important bits of information should be at the bottom as "pytest" prints out a test summary at the end.
For example:

.. code:: bash

    ============================================================================= short test summary info =============================================================================
    SKIPPED [1] d:\a\1\s\.tox\py37\lib\site-packages\pytest_doctestplus\plugin.py:178: unable to import module local('d:\\a\\1\\s\\.tox\\py37\\lib\\site-packages\\sunpy\\io\\setup_package.py')
    SKIPPED [213] d:\a\1\s\.tox\py37\lib\site-packages\pytest_remotedata\plugin.py:87: need --remote-data option to run
    SKIPPED [18] d:\a\1\s\.tox\py37\lib\site-packages\_pytest\doctest.py:387: all tests skipped by +SKIP option
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\map\sources\tests\test_source_type.py:21: Glymur can not be imported.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\map\sources\tests\test_source_type.py:30: Glymur can not be imported.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_ana.py:22: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_ana.py:31: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_ana.py:40: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_ana.py:49: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_ana.py:58: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_ana.py:67: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_filetools.py:54: Glymur can not be imported.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_filetools.py:73: Glymur can not be imported.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_filetools.py:106: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_filetools.py:115: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_filetools.py:122: ANA is not available.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_jp2.py:11: Glymur can not be imported.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_jp2.py:21: Glymur can not be imported.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\io\tests\test_jp2.py:31: Glymur can not be imported.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\net\tests\test_fido.py:298: Windows.
    SKIPPED [1] .tox\py37\lib\site-packages\sunpy\net\tests\test_helioviewer.py:90: Glymur can not be imported.
    FAILED ..\..\.tox\py37\lib\site-packages\sunpy\timeseries\sources\noaa.py::sunpy.timeseries.sources.noaa.NOAAGoesSXRTimeSeries

If you want to find the full test output, you can search the tab for the name of the test out of the ~3 results, one will be that output.

SunPy GitHub Groups
===================

This document has already referred to two SunPy groups, namely "developers" and "maintainers" there is also a third primary SunPy group "owners".

SunPy owners
------------

The SunPy owners group is the group of people who have total control over the SunPy GitHub organization.
The SunPy board have control over who is in this group, it has been decided that generally it will be the Lead Developer and the SunPy board chair and vice-chair.

sunpy Maintainers
-----------------

This is the group of people who have push access to the main sunpy repository.
The membership of this group is at the discretion of the Lead Developer, but shall generally be made up of people who have demonstrated themselves to be trust worthy and active contributors to the project.

This group has `subgroups <https://github.com/orgs/sunpy/teams/sunpy-maintainers/teams>`__ for each section of the repository that has `maintainers <https://sunpy.org/project/#maintainers>`__.
The members of these groups will automatically be requested to review all PRs which change files in that subpackage.

sunpy Developers
----------------

The members of this group have "read" access to the sunpy repository.
As all these repository are open anyway, what this effectively means is that these people can be assigned to issues.
The members of this group are people who are involved in the development of sunpy at a good frequency, they are people who's opinions have been demonstrated to be constructive and informative.
.. _dependency_versions:

*************************
Dependency Support Policy
*************************

.. note::

    This policy is based on `NEP-0029`_.

sunpy has a short list of core dependencies (Python, numpy, astropy, parfive) and a long list of optional dependencies.
The minimum version of these packages that we enforce follows this policy.

* Python: Released in the prior 42 months from the anticipated release date.
* astropy: Released in the prior 12 months from the anticipated release date.
* Everything else: Released in the prior 24 months from the anticipated release date.

Sponsored affiliated packages will support *at least* the sunpy LTS version at the time of their release.

.. _NEP-0029: https://numpy.org/neps/nep-0029-deprecation_policy.html
.. _funding:

.. warning:: HEAVY WIP

*********************
Submodule and funding
*********************

Here we list whatever parts of the SunPy project have been funded by external sources such as grants or Google Summer of Code (GSoC) and the like.

Database
========

Google Summer of Code (2013)
----------------------------

database/attrs.py
database/test_tables.py
database/test_database.py
database/test_commands.py
database/test_caching.py
database/test_attrs.py

Net
---

Google Summer of Code (2014)
----------------------------

fido_factory.py, vso.py, most of the dataretriever submodule

ESA Summer of Code in Space (2011)
----------------------------------

attr.py
tools/hektemplate.py
tools/hek_mkcls.py
.. _new_maps_ts_etc:

************************************************
Creating new sunpy Subclasses (Maps, TimeSeries)
************************************************

Writing a new Instrument Map Class
==================================

All instrument Map classes are subclasses of the generic `~sunpy.map.GenericMap` subclass.
`~sunpy.map.GenericMap` expects to be provided the data array and a header dictionary and will parse the header metadata to the best of its ability based on common standards.
The instrument subclass implements the instrument-specific code to parse the metadata, apply any necessary procedures on the data array, as well as defining other things such what color tables to use.

In practice, the instrument subclass is not directly accessed by users.
The `~sunpy.map.Map` factory is the primary interface for creating Map objects.
Any subclass of `~sunpy.map.GenericMap` which defines a method named ``is_datasource_for()`` will automatically be registered with the `~sunpy.map.Map` factory.
The ``is_datasource_for`` method is used by the `~sunpy.map.Map` factory to check if a file should use a particular instrument Map class.
This function can run any test it needs to determine this.
For example, it might check the value of the ``INSTRUMENT`` key in the metadata dictionary.
The following example shows how this works and includes a sample doc string that is compatible with the :ref:`Docs Guidelines for Data Sources`.

.. code-block:: python

    import sunpy.map
    class NextGenerationTelescopeMap(sunpy.map.GenericMap):
      "NextGenerationTelescope Map.

      The Next Generation Telescope is a type A telescope on board the XYZ mission.
      It operates in low Earth orbit with an altitude of 600 kmn and an inclination of 28.5 degrees.
      It is designed to observe the aliens on the Sun that are responsible for triggering the impulsive release of magnetic energy in the solar corona.
      It observes in the following 3 different passband in visible light, wavelength A, wavelength B, wavelength C.
      The primary emission processes in these passbands are process A and process B.

      The focal plane consists of a MAGIC detector with 2 x 2 pixels.
      The plate scale is 500 arcsec per pixel.
      The field of view is the whole Sun (1000 x 1000 arsec).
      It makes images in each passband every 10 minutes except for when it is in eclipse which occurs every approximately 30 minutes.

      It began operating on 2100 November 1.

      Notes
      -----
      Due to rise of our new insect overlords, the telescope was not operational from 2200 Jan to 2300 Jan.

      References
      ----------
      * List of all required references
      "

        def __init__(self, data, header, **kwargs):

            # will process the header according to common standards
            super(FutureMap, self).__init__(data, header, **kwargs)

            # Any NextGenerationTelescope Instrument-specific manipulation.
            # This typically involves editing the `self.meta` attribute.

        # used by the Map factory to determine if this subclass should be used
        @classmethod
        def is_datasource_for(cls, data, header, **kwargs):
            """Determines if data, header corresponds to a NextGenerationTelescope image"""
            # returns True only if this is data and header from NextGenerationTelescope
            return header.get('instrume', '').startswith('NextGenerationTelescope')


Writing a new Instrument TimeSeries Class
=========================================

To be written.
.. _maintainer-workflow:

************************
Workflow for Maintainers
************************

This page is for maintainers who can merge our own or other peoples' changes into the upstream repository.

Seeing as how you're a maintainer, you should be completely on top of the basic git workflow in :ref:`newcomers` and Astropy's `git workflow`_.

.. _git workflow: https://docs.astropy.org/en/stable/development/workflow/development_workflow.html#development-workflow

Integrating changes via the web interface (recommended)
=======================================================

Whenever possible, merge pull requests automatically via the pull request manager on GitHub.
Merging should only be done manually if there is a really good reason to do this!

Make sure that pull requests do not contain a messy history with merges, etc.
If this is the case, then follow the manual instructions, and make sure the fork is rebased to tidy the history before committing.

Integrating changes manually
============================

First, check out the "sunpy" repository.
Being a maintainer, you've got read-write access.

It's good to have your upstream remote have a scary name, to remind you that it's a read-write remote::

    $ git remote add upstream-rw git@github.com:sunpy/sunpy.git
    $ git fetch upstream-rw

Let's say you have some changes that need to go into trunk (``upstream-rw/main``).

The changes are in some branch that you are currently on.
For example, you are looking at someone's changes like this::

    $ git remote add someone git://github.com/someone/sunpy.git
    $ git fetch someone
    $ git branch cool-feature --track someone/cool-feature
    $ git checkout cool-feature

So now you are on the branch with the changes to be incorporated upstream.
The rest of this section assumes you are on this branch.

If you prefer not to add remotes, you can make git fetch all pull requests opened to Sunpy.
Locate the section for your git remote in the ``.git/config`` file.
It looks like this::

    [remote "upstream"]
            url = git@github.com:sunpy/sunpy.git
            fetch = +refs/heads/*:refs/remotes/upstream/*

Now add the line ``fetch = +refs/pull/*/head:refs/remotes/upstream/pr/*`` to this section.
It ends up looking like this::

    [remote "upstream"]
            url = git@github.com:sunpy/sunpy.git
            fetch = +refs/heads/*:refs/remotes/upstream/*
            fetch = +refs/pull/*/head:refs/remotes/upstream/pr/*

Now fetch all the pull requests::

    $ git fetch upstream
    From github.com:sunpy/sunpy
    * [new ref]         refs/pull/1000/head -> upstream/pr/1000
    * [new ref]         refs/pull/1002/head -> upstream/pr/1002
    * [new ref]         refs/pull/1004/head -> upstream/pr/1004
    * [new ref]         refs/pull/1009/head -> upstream/pr/1009

To check out a particular pull request::

    $ git checkout pr/999
    Branch pr/999 set up to track remote branch pr/999 from upstream.
    Switched to a new branch 'pr/999'

When to remove or combine/squash commits
----------------------------------------

In all cases, be mindful of maintaining a welcoming environment and be helpful with advice, especially for new contributors.
It is expected that a maintainer would offer to help a contributor who is a novice git user do any squashing that that maintainer asks for, or do the squash themselves by directly pushing to the PR branch.

Pull requests **must** be rebased and at least partially squashed (but not necessarily squashed to a single commit) if large (approximately >10KB) non-source code files (e.g. images, data files, etc.) are added and then removed or modified in the PR commit history (The squashing should remove all but the last addition of the file to not use extra space in the repository).

Combining/squashing commits is **encouraged** when the number of commits is excessive for the changes made.
The definition of "excessive" is subjective, but in general one should attempt to have individual commits be units of change, and not include reversions.
As a concrete example, for a change affecting < 50 lines of source code and including a changelog entry, more than a two commits would be excessive.
For a larger pull request adding significant functionality, however, more commits may well be appropriate.

As another guideline, squashing should remove extraneous information but should not be used to remove useful information for how a PR was developed.
For example, 4 commits that are testing changes and have a commit message of just "debug" should be squashed.
But a series of commit messages that are "Implemented feature X", "added test for feature X", "fixed bugs revealed by tests for feature X" are useful information and should not be squashed away without reason.

When squashing, extra care should be taken to keep authorship credit to all individuals who provided substantial contribution to the given PR, e.g. only squash commits made by the same author.

When to rebase
--------------

Pull requests **must** be rebased (but not necessarily squashed to a single commit) if:

* There are commit messages include offensive language or violate the code of conduct (in this case the rebase must also edit the commit messages)

Pull requests **may** be rebased (either manually or with the rebase and merge button) if:

* There are conflicts with main
* There are merge commits from upstream/main in the PR commit history (merge commits from PRs to the user's fork are fine)

Asking contributors who are new to the project or inexperienced with using git is **discouraged**, as is maintainers rebasing these PRs before merge time, as this requires resetting of local git checkouts.


A few commits
-------------

If there are only a few commits, consider rebasing to upstream::

    # Fetch upstream changes
    $ git fetch upstream-rw

    # Rebase
    $ git rebase upstream-rw/main

A long series of commits
------------------------

If there are a longer series of related commits, consider a merge instead::

    $ git fetch upstream-rw
    $ git merge --no-ff upstream-rw/main

Note the ``--no-ff`` above.
This forces git to make a merge commit, rather than doing a fast-forward, so that these set of commits branch off trunk then rejoin the main history with a merge, rather than appearing to have been made directly on top of trunk.

Check the history
-----------------

Now, in either case, you should check that the history is sensible and you have the right commits::

    $ git log --oneline --graph
    $ git log -p upstream-rw/main..

The first line above just shows the history in a compact way, with a text representation of the history graph.
The second line shows the log of commits excluding those that can be reached from trunk (``upstream-rw/main``), and including those that can be reached from current HEAD (implied with the ``..`` at the end).
So, it shows the commits unique to this branch compared to trunk.
The ``-p`` option shows the diff for these commits in patch form.

Push to open pull request
-------------------------

Now you need to push the changes you have made to the code to the open pull request::

    $ git push git@github.com:<username>/sunpy.git HEAD:<name of branch>

You might have to add ``--force`` if you rebased instead of adding new commits.

Using Milestones and Labels
===========================

Current milestone guidelines:

* Only confirmed issues or pull requests that are release critical or for some other reason should be addressed before a release, should have a milestone.
  When in doubt about which milestone to use for an issue, do not use a milestone and ask other the maintainers.

Current labelling guidelines:

* Issues that require fixing in main, but that also are confirmed to apply to supported stable version lines should be marked with a "Affects Release" label.
* All open issues should have a "Priority <level>", "Effort <level>" and "Package <level>", if you are unsure at what level, pick higher ones just to be safe.
  If an issue is more of a question or discussion, you can omit these labels.
* If an issue looks to be straightforward, you should add the "Good first issue" and "Hacktoberfest" label.
* For other labels, you should add them if they fit, like if an issue affects the net submodule, add the "net" label or if it is a feature request etc.

Using Projects
==============

Projects allow us to layout current pull requests and issues in a manner that enables a more "meta" view regarding major releases.
We categorize pull requests and issues into several levels of priorities and whether these can be classed as blockers before a release can be attempted.
Further we can add general notes that someone deems important for a release.

Updating and Maintaining the Changelog
======================================

The changelog will be read by users, so this description should be aimed at sunpy users instead of describing internal changes which are only relevant to the developers.

The current changelog is kept in the file "CHANGELOG.rst" at the root of the repository.
You do not need to update this file as we use `towncrier`_ to update our changelog.
This is built and embedded into our documentation.

Towncrier will automatically reflow your text, so it will work best if you stick to a single paragraph, but multiple sentences and links are OK and encouraged.
You can install towncrier and then run ``towncrier --draft`` if you want to get a preview of how your change will look in the final release notes.

`Instructions on how to write a changelog. <https://github.com/sunpy/sunpy/blob/main/changelog/README.rst>`__.

.. _towncrier: https://pypi.org/project/towncrier/

Releasing sunpy
===============

We have a `step by step checklist`_ on the SunPy Wiki on how to release the sunpy core package.

.. _step by step checklist: https://github.com/sunpy/sunpy/wiki/Home%3A-Release-Checklist
.. doctest-skip-all

.. _units_in_code:

***************************
Use of quantities and units
***************************

Much code perform calculations using physical quantities.
SunPy uses astropy's `quantities and units <https://docs.astropy.org/en/stable/units/index.html>`_ implementation to store, express and convert physical quantities.
New classes and functions should adhere to the SunPy project's `quantity and unit usage guidelines
<https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`_.

This document sets out SunPy's reasons and requirements for the usage of quantities and units.
Briefly, SunPy's `policy <https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`_ is that *all user-facing function/object arguments which accept physical quantities as input **MUST** accept astropy quantities, and **ONLY** astropy quantities*.

Developers should consult the `Astropy Quantities and Units page <https://docs.astropy.org/en/stable/units/index.html>`_ for the latest updates on using quantities and units.  The `astropy tutorial on quantities and units <https://www.astropy.org/astropy-tutorials/Quantities.html>`_ also provides useful examples on their
capabilities.

Astropy provides the decorator `~astropy.units.quantity_input` that checks the units of the input arguments to a function against the expected units of the argument.
We recommend using this decorator to perform function argument unit checks.
The decorator ensures that the units of the input to the function are convertible to that specified by the decorator, for example ::

    >>> import astropy.units as u
    >>> @u.quantity_input
    ... def myfunction(myangle: u.arcsec):
    ...     return myangle**2

This function only accepts arguments that are convertible to arcseconds.
Therefore::

    >>> myfunction(20 * u.degree)
    <Quantity 400. deg2>

returns the expected answer but::

    >>> myfunction(20 * u.km)
    Traceback (most recent call last):
    ...
    astropy.units.core.UnitsError: Argument 'myangle' to function 'myfunction' must be in units convertible to 'arcsec'.

raises an error.

The following is an example of a use-facing function that returns the area of a square, in units that are the square of the input length unit::

    >>> @u.quantity_input
    ... def get_area_of_square(side_length: u.m):
    ...     """
    ...     Compute the area of a square.
    ...
    ...     Parameters
    ...     ----------
    ...     side_length : `~astropy.units.quantity.Quantity`
    ...         Side length of the square
    ...
    ...     Returns
    ...     -------
    ...     area : `~astropy.units.quantity.Quantity`
    ...         Area of the square.
    ...     """
    ...
    ...     return (side_length ** 2)

This more advanced example shows how a private function that does not accept quantities can be wrapped by a function that does::

    >>> @u.quantity_input
    ... def some_function(length: u.m):
    ...     """
    ...     Does something useful.
    ...
    ...     Parameters
    ...     ----------
    ...     length : `~astropy.units.quantity.Quantity`
    ...         A length.
    ...
    ...     Returns
    ...     -------
    ...     length : `~astropy.units.quantity.Quantity`
    ...         Another length
    ...     """
    ...
    ...     # the following function either
    ...     # a] does not accept Quantities
    ...     # b] is slow if using Quantities
    ...     result = _private_wrapper_function(length.convert('meters').value)
    ...
    ...     # now convert back to a quantity
    ...     result = Quantity(result_meters, units_of_the_private_wrapper_function)
    ...
    ...     return result

In this example, the non-user facing function ``_private_wrapper_function`` requires a numerical input in units of meters, and returns a numerical output.
The developer knows that the result of ``_private_wrapper_function`` is in the units ``units_of_the_private_wrapper_function``, and sets the result of ``some_function`` to return the answer in those units.
.. _example_gallery:

***************
Example Gallery
***************

The purpose of the page is to describe the contribution guidelines for the `sunpy Example Gallery <https://docs.sunpy.org/en/stable/generated/gallery/index.html>`_.

All potential contributors to the sunpy Example Gallery should read and abide by the following guidelines.

.. note:: We have an example template located at ``examples/example_template/example_template.py``.

Contribution Guidelines
=======================

* The title of the example should be short yet descriptive and emphasize the goal of the example.
  Try to make the title appeal to a broad audience and avoid referencing a specific instrument, catalog, or anything wavelength dependent.

* Each example should begin with a paragraph that gives a brief overview of the entire example, including relevant astronomy concepts, and motivates the described functionality.

* The examples must be compatible with the versions supported by the last major release of the sunpy core package (i.e., Python >= 3.7).

* All the examples must be fully PEP8 compliant, we recommend using one of the many PEP8 linters that are available (autopep8, flake8 as some examples).

* Wherever possible, the examples should include linked references with links pointing to the appropriate `DOI <https://zenodo.org/record/2551710>`_ or `ADS <https://ui.adsabs.harvard.edu/>`_ entry.

* The example should include links to relevant documentation pages.

* Each example should, where possible, include at least one image, map, or plot to use as the icon in the example gallery.

* The examples should avoid using acronyms without defining them first (e.g. Virtual Solar Observatory, or VSO).
  Similarly complex jargon should be avoided unless clearly explained.

* There should be a good variety of examples for each section (simple and more complex to cater for different levels).
.. _public_api:

******************
sunpy's Public API
******************

Convention in the Python ecosystem is to add an underscore to the start of a function to denote if a function or method is "private" e.g., `~sunpy.coordinates.sun._angular_radius`.
If it is considered to be private, there is no guarantee that the API or behavior will change with a warning, so external use of these functions or methods are strongly discouraged.

sunpy follows this convention but with one extra caveat.
Within each python file, we have a ``__all__`` that defines what is imported into the namespace if you do e.g., ``from sunpy.coordinates.sun import *``.
This is the "public" API of that module.
These functions are the ones listed within our API documentation: :ref:`reference`.
If you do ``import sunpy.coordinates.sun``, you can still access the "private" functions.

This means that all of the public API will follow the deprecation policy detailed below with the exception of `sunpy.util` which is considered to be for internal sunpy use only.

Deprecation Policy and Breaking Changes
=======================================

All public API within the SunPy project (the sunpy core package and stable affiliated packages) will enforce strict standards when it comes to either changing, updating or breaking the API.

.. _deprecation:

Deprecations
------------

If you want to deprecate anything within in sunpy, you should do the following:

.. code-block:: python

    from sunpy.util.decorators import deprecated

    @deprecated(since="3.1", message="We will be moving this", alternative="sunpy.net.Scraper")
    class Scraper:

The deprecation warning has to be in one LTS release before the deprecated code can be removed.
So in the above example, the warning will be in sunpy 3.1 but it can not be removed until sunpy 4.1 after the 4.0 LTS release.

There should be a "deprecation" changelog entry to accompany the deprecation warning.
When the code is actually removed, a "removal" changelog should be added.

The same applies if you want to change the default value of a keyword argument for a function or method, e.g.:

.. code-block:: python

    if response_format is None:
        response_format = "legacy"
        warnings.warn("The default response format from the VSO client will "
                    "be changing to 'table' in version 3.1. "
                    "To remove this warning set response_format='legacy' "
                    "to maintain the old behaviour or response_format='table'"
                    " to use the new behaviour.",
                    SunpyDeprecationWarning,
                    stacklevel=2)

.. note::

    This is a summary of `SEP-0009`_ which is the formal SunPy project deprecation policy.

.. _SEP-0009: https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0009.md#deprecations-and-documentation

.. _breaking:

Breaking Changes
----------------

Every attempt is made to avoid breaking any public API in sunpy but in the case it does happen.

There should be a "breaking" changelog entry to accompany the change with as much detail on how to update a user's code due to this breaking change.
.. _timeseries_code_ref:

Timeseries (`sunpy.timeseries`)
*******************************
One of the core classes in sunpy is a timeseries.
A number of instruments are supported through subclasses of the base `~sunpy.timeseries.GenericTimeSeries` class.
See :ref:`ts-sources` for a list of them.

A timeseries can be created by calling `~sunpy.timeseries.TimeSeries`.

.. _ts-sources:

Instrument TimeSeries Classes
=============================
The generic method to create an instrument-specific TimeSeries is to call `~sunpy.timeseries.TimeSeries` with a file path and the instrument-specific source keyword argument.
In some cases the source can be determined automatically if a FITS file is being loaded.

The following example shows the factory loading a sample file::

    >>> import sunpy.timeseries as ts
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> goes = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')  # doctest: +REMOTE_DATA

The `~sunpy.timeseries.TimeSeries` factory will load the file and create the timeseries instance.
The following instrument classes are supported:

.. automodapi:: sunpy.timeseries
    :no-inheritance-diagram:
    :include-all-objects:

.. automodapi:: sunpy.timeseries.sources

CDF files
=========
`~sunpy.timeseries.GenericTimeSeries` can load a single CDF file, or a list of CDF files if ``concatenate=True`` is passed.

Units
-----
The physical units of different columns in CDF files do not conform to a standard that `astropy.units` understands.
sunpy internally stores a set of common mappings from unit strings to `~astropy.units.Unit`, but you may see a warning about unrecognised unit strings when reading a CDF file.
To register the correct unit definition :func:`astropy.units.add_enabled_units` can be used.
For example, to register 'deg K' as representing Kelvin and '#/cc' as 1/cm^3::

  >>> import astropy.units as u
  >>> _ = u.add_enabled_units([u.def_unit('deg K', represents=u.K), u.def_unit('#/cc', represents=u.cm**-3)])
************
Known Issues
************

This page documents commonly known issues, issues here is defined broadly and refers to oddities or specifics of how sunpy or the Python ecosystem works that could anyone catch out.

Disagreement between `astropy.time.Time` and SSW for pre-1972 dates
===================================================================

Conversions between time scales for pre-1972 dates may not agree between `astropy.time.Time` and SSW.
Prior to 1972, the length of a Universal Time (UT) second was different from the length of an SI second, and `astropy.time.Time` correctly handles this difference.
However, SSW does not (at least in the commonly used routines), and thus SSW conversions for pre-1972 dates between UT and International Atomic Time (TAI, based on SI seconds) can be incorrect by multiple seconds.
The SSW routine ``utc2tai`` notes that its output is invalid for pre-1972 dates, but other SSW routines are missing such a disclaimer.

`Reference issue <https://github.com/sunpy/sunpy/issues/5500>`__.

`~sunpy.net.jsoc.JSOCClient` time filtering
===========================================

The JSOC API filters on ``T_OBS`` for some series, instead of ``T_REC`` which is what is used most of the time.
There is not anything we can do to fix this on our side.
`The specific section of the SDO users guide explaining this <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/sdoguidese4.html#x9-240004.2.4>`__.

`Reference issue <https://github.com/sunpy/sunpy/issues/5447>`__.
Database (`sunpy.database`)
***************************

``sunpy.database`` can be used to provide a local cache of the files and
records retrieved from various remote services. For an introduction to the
database see :ref:`database_guide`.


.. automodapi:: sunpy.database
    :no-inheritance-diagram:

Submodules
==========

.. automodapi:: sunpy.database.tables
    :no-inheritance-diagram:

.. automodapi:: sunpy.database.caching
    :no-inheritance-diagram:

.. automodapi:: sunpy.database.commands
    :no-inheritance-diagram:

.. automodapi:: sunpy.database.attrs
    :no-inheritance-diagram:
Visualization (`sunpy.visualization`)
*************************************

`sunpy.visualization` contains plotting helpers and functions.

.. automodapi:: sunpy.visualization

.. automodapi:: sunpy.visualization.colormaps

.. automodapi:: sunpy.visualization.colormaps.color_tables

.. automodapi:: sunpy.visualization.animator

.. automodapi:: sunpy.visualization.wcsaxes_compat
sunpy
*****

.. automodapi:: sunpy
Physics (`sunpy.physics`)
*************************

``sunpy.physics`` contains routines to calculate various physical parameters
and models for the Sun.

.. automodapi:: sunpy.physics

.. automodapi:: sunpy.physics.differential_rotation

.. automodapi:: sunpy.physics.solar_rotation
Data (`sunpy.data`)
*******************

``sunpy.data`` contains ways to access sample data and small test files for
running the sunpy test suite.

.. automodapi:: sunpy.data
   :include-all-objects:

.. automodapi:: sunpy.data.sample

.. automodapi:: sunpy.data.test

.. automodapi:: sunpy.data.data_manager

.. automodapi:: sunpy.data.data_manager.downloader

.. automodapi:: sunpy.data.data_manager.storage
Utilities (`sunpy.util`)
************************

.. automodapi:: sunpy.util

.. automodapi:: sunpy.util.config

.. automodapi:: sunpy.util.datatype_factory_base

.. automodapi:: sunpy.util.net

.. automodapi:: sunpy.util.xml

.. automodapi:: sunpy.util.sphinx

.. automodapi:: sunpy.util.sphinx.generate
   :headings: ^"

.. automodapi:: sunpy.util.functools

.. automodapi:: sunpy.util.types
Remote data (`sunpy.net`)
*************************

``sunpy.net`` contains a lot of different code for accessing various solar
physics related web services. This submodule contains many layers. Most users
should use `~sunpy.net.Fido`, which is an interface to multiple sources
including all the sources implemented in `~sunpy.net.dataretriever` as well as
`~sunpy.net.vso` and `~sunpy.net.jsoc`. `~sunpy.net.Fido` can be used like so::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument.lyra)  # doctest: +REMOTE_DATA
    >>> files = Fido.fetch(results)  # doctest: +SKIP

.. automodapi:: sunpy.net
   :include-all-objects:
   :headings: =-

.. automodapi:: sunpy.net.attrs
   :headings: ^"

.. automodapi:: sunpy.net.fido_factory
   :headings: ^"

VSO
---

.. automodapi:: sunpy.net.vso
   :headings: ^"

.. automodapi:: sunpy.net.vso.attrs
   :headings: ^"

Dataretriever
-------------

.. automodapi:: sunpy.net.dataretriever
   :headings: ^"

.. automodapi:: sunpy.net.dataretriever.attrs.goes
   :headings: ^"

JSOC
----

.. automodapi:: sunpy.net.jsoc
   :headings: ^"

.. automodapi:: sunpy.net.jsoc.attrs
   :headings: ^"

HEK
---

.. automodapi:: sunpy.net.hek
   :headings: ^"

.. automodapi:: sunpy.net.hek.attrs
   :headings: ^"

.. automodapi:: sunpy.net.hek2vso
   :headings: ^"

CDAWeb
------

.. automodapi:: sunpy.net.cdaweb
   :headings: ^"


HELIO
-----

.. automodapi:: sunpy.net.helio
   :headings: ^"

.. automodapi:: sunpy.net.helio.attrs
   :headings: ^"

Helioviewer
-----------

.. automodapi:: sunpy.net.helioviewer
   :headings: ^"

Internal Classes and Functions
==============================

These classes and functions are designed to be used to help develop new clients
for `sunpy.net.Fido`.

.. automodapi:: sunpy.net.base_client

.. automodapi:: sunpy.net.dataretriever.client

.. automodapi:: sunpy.net.attr

.. automodapi:: sunpy.net.scraper
.. _sun_code_ref:

Solar properties (``sunpy.sun``)
********************************

``sunpy.sun`` contains constants, parameters and models of the Sun.

.. automodapi:: sunpy.sun.constants
   :include-all-objects:

.. automodapi:: sunpy.sun.models
   :include-all-objects:
Current status of sub-packages
******************************

sunpy has variations of stability across different sub-packages.
This document summarizes the current status of the sunpy sub-packages, so that users understand where they might expect changes in future, and which sub-packages they can safely use for production code.
For help or to discuss specific sub-packages, refer to the `maintainer list <https://sunpy.org/project/#maintainers>`__ to find who to contact.

The classification is as follows:

{% raw %}
.. raw:: html
{% endraw %}
   <style>
         .planned:before {
              color: #cbcbcb;
              content: "⬤";
              padding: 5px;
         }
         .dev:before {
              color: #ffad00;
              content: "⬤";
              padding: 5px;
         }
         .stable:before {
              color: #4e72c3;
              content: "⬤";
              padding: 5px;
         }
         .mature:before {
              color: #03a913;
              content: "⬤";
              padding: 5px;
         }
         .pendingdep:before {
              color: #a84b03;
              content: "⬤";
              padding: 5px;
         }
         .deprecated:before {
              color: #ff0000;
              content: "⬤";
              padding: 5px;
         }
    </style>

    <table align='center'>
      <tr>
        <td align='center'><span class="planned"></span></td>
        <td>Planned</td>
      </tr>
      <tr>
        <td align='center'><span class="dev"></span></td>
        <td>Actively developed, be prepared for possible significant changes.</td>
      </tr>
      <tr>
        <td align='center'><span class="stable"></span></td>
        <td>Reasonably stable, any significant changes/additions will generally include backwards-compatiblity.</td>
      </tr>
      <tr>
        <td align='center'><span class="mature"></span></td>
        <td>Mature.  Additions/improvements possible, but no major changes planned. </td>
      </tr>
      <tr>
        <td align='center'><span class="pendingdep"></span></td>
        <td>Pending deprecation.  Might be deprecated in a future version.</td>
      </tr>
      <tr>
        <td align='center'><span class="deprecated"></span></td>
        <td>Deprecated.  Might be removed in a future version.</td>
      </tr>
    </table>

The current planned and existing sub-packages are:

{% raw %}
.. raw:: html
{% endraw %}

    <table border="1" class="docutils stability" align='center'>
        <tr>
            <th class="head">
                Sub-Package
            </th>
            <th class="head">
                &nbsp;
            </th>
            <th class="head">
                Comments
            </th>
        </tr>
    {% for module, prop in sunpy_modules.items() %}
        <tr>
            <td>
                <a href="../code_ref/{{ module }}.html">sunpy.{{ module }}</a>
            </td>
            <td align='center'>
                <span class="{{ prop['status'] }}"></span>
            </td>
            <td>
                {{ prop['comments'] }}
            </td>
        </tr>
    {% endfor %}
    </table>


Taken with love from the `Astropy project. <https://github.com/astropy/astropy/blob/master/LICENSE.rst>`_
Image processing (`sunpy.image`)
********************************

``sunpy.image`` contains routines to process images (i.e. `numpy` arrays). The
routines in this submodule are generally exposed through map-specific functions
in other places.

.. automodapi:: sunpy.image

.. automodapi:: sunpy.image.resample

.. automodapi:: sunpy.image.transform

.. automodapi:: sunpy.image.coalignment
.. _map:

Maps (`sunpy.map`)
******************

.. module:: sunpy.map

.. testsetup::

    import sunpy

.. currentmodule:: sunpy.map

Overview
========
One of core classes in sunpy is a Map. A sunpy Map object is simply a
spatially-aware data array, often an image. In order to make it easy to work
with image data in sunpy, the Map object provides a number of methods for
commonly performed operations.

2D map objects are subclasses of `~sunpy.map.GenericMap` and these objects are
created using the Map factory `~sunpy.map.Map`.

A number of instrument are supported by subclassing this base object. See
:ref:`map-sources` to see a list of all of them. More complex subclasses are also
available. See :ref:`map-classes`.

.. todo :
    1. Map factory and registration
    2. MapBase and Generic Map
    3. MapMeta and the separation from the file io


Creating Map Objects
====================
sunpy Map objects are constructed using the special factory
class `~sunpy.map.Map`: ::

    >>> x = sunpy.map.Map('file.fits')  # doctest: +SKIP

The result of a call to `~sunpy.map.Map` will be either a `~sunpy.map.mapbase.GenericMap` object,
or a subclass of `~sunpy.map.mapbase.GenericMap` which either deals with a specific type of data,
e.g. `~sunpy.map.sources.sdo.AIAMap` or `~sunpy.map.sources.soho.LASCOMap`
(see :ref:`map-classes` to see a list of all of them), or if no
instrument matches, a 2D map `~sunpy.map.mapbase.GenericMap`.

Fixing map metadata
-------------------

If you need to fix the metadata of a fits file before it is handed to `Map`, this can be done as
follows:

    >>> data, header = sunpy.io.fits.read(filepath)[0] # doctest: +SKIP
    >>> header['cunit1'] = 'arcsec' # doctest: +SKIP
    >>> header['cunit2'] = 'arcsec' # doctest: +SKIP
    >>> map = sunpy.map.Map(data, header) # doctest: +SKIP


.. autoclass:: sunpy.map.map_factory.MapFactory

.. _map-classes:

sunpy.map Package
=================

All sunpy Maps are derived from `sunpy.map.GenericMap`, all the methods and attributes are documented in that class.

.. automodapi:: sunpy.map
    :no-main-docstr:
    :inherited-members:
    :include-all-objects:

.. _map-sources:

Instrument Map Classes
======================
Defined in ``sunpy.map.sources`` are a set of `~sunpy.map.GenericMap` subclasses
which convert the specific metadata and other differences in each instruments
data to the standard `~sunpy.map.GenericMap` interface.
These 'sources' also define things like the colormap and default
normalisation for each instrument.
These subclasses also provide a method, which describes to the `~sunpy.map.Map` factory
which data and metadata pairs match its instrument.

.. automodapi:: sunpy.map.sources
    :no-main-docstr:
    :inherited-members:


Writing a new Instrument Map Class
==================================

Any subclass of `~sunpy.map.GenericMap` which defines a method named
``is_datasource_for`` will automatically be registered with
the `~sunpy.map.Map` factory. The ``is_datasource_for`` method describes the form of the
data and metadata for which the `~sunpy.map.GenericMap` subclass is valid. For
example it might check the value of the ``INSTRUMENT`` key in the metadata
dictionary.
This makes it straightforward to define your own
`~sunpy.map.GenericMap` subclass for a new instrument or a custom data source
like simulated data. These classes only have to be imported for this to work, as
demonstrated by the following example.

.. code-block:: python

    import sunpy.map
    class FutureMap(sunpy.map.GenericMap):

        def __init__(self, data, header, **kwargs):

            super(FutureMap, self).__init__(data, header, **kwargs)

            # Any Future Instrument specific keyword manipulation

       # Specify a classmethod that determines if the data-header pair matches
       # the new instrument
       @classmethod
       def is_datasource_for(cls, data, header, **kwargs):
            """Determines if header corresponds to an AIA image"""
            return str(header.get('instrume', '')).startswith('FUTURESCOPE')


This class will now be available through the `~sunpy.map.Map` factory as long as this
class has been defined, i.e. imported into the current session.

If you do not want to create a method named ``is_datasource_for`` you can
manually register your class and matching method using the following method

.. code-block:: python

    import sunpy.map

    sunpy.map.Map.register(FutureMap, FutureMap.some_matching_method)
Time (`sunpy.time`)
*******************

``sunpy.time`` contains helpers for converting strings to `astropy.time.Time`
objects and handling common operations on these objects. As well as this a
`~sunpy.time.TimeRange` object is provided for representing a period of time
and performing operations on that range.

.. automodapi:: sunpy.time
.. _reference:

*************
API Reference
*************

.. toctree::
   :maxdepth: 2

   coordinates/index
   data
   database
   image
   io
   map
   net
   physics
   sun
   sunpy
   time
   timeseries
   util
   visualization
Input/output (`sunpy.io`)
*************************

``sunpy.io`` contains two types of routines, the first reads (data, header)
pairs from files in a way similar to FITS files. The other is special readers
for files that are commonly used in solar physics.

.. warning::

   When reading FITS files, it is strongly recommended that `astropy.io.fits` be used over the tools in `sunpy.io`.
   The sunpy FITS reader is designed to meet the needs of map, and does not represent the structure of the FITS file well.


Unified File Readers
====================

.. automodapi:: sunpy.io

.. automodapi:: sunpy.io.header

.. _iofits:
.. automodapi:: sunpy.io.fits

.. _iojp2:
.. automodapi:: sunpy.io.jp2

.. _ioana:
.. automodapi:: sunpy.io.ana

Special File Readers
====================

.. _iospecialgenx:
.. automodapi:: sunpy.io.special.genx


asdf (Advanced Scientific Data Format)
--------------------------------------

`asdf <https://asdf.readthedocs.io/en/latest/>`__ is a modern file format
designed to meet the needs of the astronomy community. It has deep integration
with Python, SunPy, and Astropy, as well as implementations in other languages.
It can be used to store known Python objects in a portable, well defined file
format. It is primarily useful for storing complex Astropy and SunPy objects in
a way that can be loaded back into the same form as they were saved.

sunpy currently implements support for saving `Map <sunpy.map.GenericMap>` and
`coordinate frame <sunpy.coordinates.frames>` objects into asdf files. As asdf
tightly integrates into Python, saving a map to an asdf file will save the
metadata, data, mask and the shift. The mask and shift are not currently saved
to FITS. The following code shows to to save and load a sunpy Map to an asdf
file

.. doctest-requires:: asdf

  >>> import asdf
  >>> import sunpy.map
  >>> from sunpy.data.sample import AIA_171_IMAGE  # doctest: +REMOTE_DATA
  >>> aiamap = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
  >>> tree = {'amap': aiamap}  # doctest: +REMOTE_DATA
  >>> with asdf.AsdfFile(tree) as asdf_file:  # doctest: +REMOTE_DATA
  ...     asdf_file.write_to("sunpy_map.asdf")  # doctest: +REMOTE_DATA
  >>> input_asdf = asdf.open("sunpy_map.asdf")  # doctest: +REMOTE_DATA
  >>> input_asdf['amap']  # doctest: +REMOTE_DATA
    <sunpy.map.sources.sdo.AIAMap object at ...>
    SunPy Map
    ---------
    Observatory:                 SDO
    Instrument:          AIA 3
    Detector:            AIA
    Measurement:                 171.0 Angstrom
    Wavelength:          171.0 Angstrom
    Observation Date:    2011-06-07 06:33:02
    Exposure Time:               0.234256 s
    Dimension:           [1024. 1024.] pix
    Coordinate System:   helioprojective
    Scale:                       [2.402792 2.402792] arcsec / pix
    Reference Pixel:     [511.5 511.5] pix
    Reference Coord:     [3.22309951 1.38578135] arcsec
    array([[ -95.92475  ,    7.076416 ,   -1.9656711, ..., -127.96519  ,
            -127.96519  , -127.96519  ],
           [ -96.97533  ,   -5.1167884,    0.       , ...,  -98.924576 ,
            -104.04137  , -127.919716 ],
           [ -93.99607  ,    1.0189276,   -4.0757103, ...,   -5.094638 ,
             -37.95505  , -127.87541  ],
           ...,
           [-128.01454  , -128.01454  , -128.01454  , ..., -128.01454  ,
            -128.01454  , -128.01454  ],
           [-127.899666 , -127.899666 , -127.899666 , ..., -127.899666 ,
            -127.899666 , -127.899666 ],
           [-128.03072  , -128.03072  , -128.03072  , ..., -128.03072  ,
            -128.03072  , -128.03072  ]], dtype=float32)
   >>> input_asdf.close()  # doctest: +REMOTE_DATA

CDF (common data format)
------------------------

CDF files are commonly used to store timeseries data observed by instruments
taking in-situ measurements of plasmas throughout the heliosphere.
sunpy provides support for reading in CDF files that conform to the
`Space Physics Guidelines for CDF <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`_.

.. automodapi:: sunpy.io.cdf
.. _sunpy-coordinates-rotatedsunframe:

Differential rotation using coordinate frames
*********************************************

..
  >>> # Due to small differences depending on different processors in numpy 1.22,
  >>> # reduce precision at which results are printed. See https://github.com/matplotlib/matplotlib/pull/21634#issuecomment-1004200517
  >>> # for the likely reason this is needed.
  >>> import numpy as np
  >>> np.set_printoptions(precision=6)

Normally, coordinates refer to a point in inertial space (relative to the barycenter of the solar system).
Transforming to a different observation time does not move the point itself, but if the coordinate frame origin and/or axis change with time, the coordinate representation is updated to account for this change.
In solar physics an example of a frame that changes with time is the `~sunpy.coordinates.frames.HeliographicStonyhurst` (HGS) frame.
Its origin moves with the center of the Sun, and its orientation rotates such that the longitude component of Earth is zero at any given time.
A coordinate in a HGS frame of reference transformed to a HGS frame defined a day later will have a different longitude::

  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> from sunpy.coordinates import HeliographicStonyhurst

  >>> hgs_coord = SkyCoord(0*u.deg, 0*u.deg, radius=1*u.au, frame='heliographic_stonyhurst', obstime="2001-01-01")
  >>> new_frame = HeliographicStonyhurst(obstime="2001-01-02")
  >>> new_hgs_coord = hgs_coord.transform_to(new_frame)
  >>> hgs_coord.lon, new_hgs_coord.lon
  (<Longitude 0. deg>, <Longitude -1.01372559 deg>)

but when transformed to an inertial frame of reference we can see that these two coordinates refer to the same point in space::

  >>> hgs_coord.transform_to('icrs')
  <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, AU)
      (101.79107615, 26.05004621, 0.99601156)>
  >>> new_hgs_coord.transform_to('icrs')
  <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, AU)
      (101.79107615, 26.05004621, 0.99601156)>

To evolve a coordinate in time such that it accounts for the rotational motion of the Sun, one can use the `~sunpy.coordinates.metaframes.RotatedSunFrame` "metaframe" class as described below.
This machinery will take into account the latitude-dependent rotation rate of the solar surface, also known as differential rotation.
Multiple models for differential rotation are supported (see :func:`~sunpy.physics.differential_rotation.diff_rot` for details).

.. note::
   `~sunpy.coordinates.metaframes.RotatedSunFrame` is a powerful metaframe, but can be tricky to use correctly.
   We recommend users to first check if the simpler :func:`~sunpy.coordinates.propagate_with_solar_surface` context manager is sufficient for their needs.

In addition, one may want to account for the translational motion of the Sun as well, and that can be achieved by also using the context manager :func:`~sunpy.coordinates.transform_with_sun_center` for desired coordinate transformations.

Basics of the RotatedSunFrame class
===================================
The `~sunpy.coordinates.metaframes.RotatedSunFrame` class allows one to specify coordinates in a coordinate frame prior to an amount of solar (differential) rotation being applied.
That is, the coordinate will point to a location in inertial space at some time, but will use a coordinate system at a *different* time to refer to that point, while accounting for the differential rotation between those two times.

`~sunpy.coordinates.metaframes.RotatedSunFrame` is not itself a coordinate frame, but is instead a "metaframe".
A new frame class is created on the fly corresponding to each base coordinate frame class.
This tutorial will refer to these new classes as ``RotatedSun*`` frames.

Creating coordinates
--------------------

`~sunpy.coordinates.metaframes.RotatedSunFrame` requires two inputs: the base coordinate frame and the duration of solar rotation.
The base coordinate frame needs to be fully specified, which means a defined ``obstime`` and, if relevant, a defined ``observer``.
Note that the ``RotatedSun*`` frame that is created in this example is appropriately named ``RotatedSunHeliographicStonyhurst``::

  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> from sunpy.coordinates import RotatedSunFrame
  >>> import sunpy.coordinates.frames as f

  >>> base_frame = f.HeliographicStonyhurst(obstime="2001-01-01")
  >>> rs_hgs = RotatedSunFrame(base=base_frame, duration=1*u.day)
  >>> rs_hgs
  <RotatedSunHeliographicStonyhurst Frame (base=<HeliographicStonyhurst Frame (obstime=2001-01-01T00:00:00.000, rsun=695700.0 km)>, duration=1.0 d, rotation_model=howard)>

Once a ``RotatedSun*`` frame is created, it can be used in the same manner as other frames.  Here, we create a `~astropy.coordinates.SkyCoord` using the ``RotatedSun*`` frame::

  >>> rotated_coord = SkyCoord(0*u.deg, 0*u.deg, frame=rs_hgs)
  >>> rotated_coord
  <SkyCoord (RotatedSunHeliographicStonyhurst: base=<HeliographicStonyhurst Frame (obstime=2001-01-01T00:00:00.000, rsun=695700.0 km)>, duration=1.0 d, rotation_model=howard): (lon, lat) in deg
        (0., 0.)>

Transforming this into the original heliographic Stonyhurst frame, we can see that the longitude is equal to the original zero degrees, plus an extra offset to account for one day of differential rotation::

  >>> rotated_coord.transform_to(base_frame).lon
  <Longitude 14.32632838 deg>

Instead of explicitly specifying the duration of solar rotation, one can use the keyword argument ``rotated_time``.
The duration will be automatically calculated from the difference between ``rotated_time`` and the ``obstime`` value of the base coordinate frame.
Here, we also include coordinate data in the supplied base coordinate frame::

  >>> rs_hgc = RotatedSunFrame(base=f.HeliographicCarrington(10*u.deg, 20*u.deg, observer="earth",
  ...                                                        obstime="2020-03-04 00:00"),
  ...                          rotated_time="2020-03-06 12:00")
  >>> rs_hgc
  <RotatedSunHeliographicCarrington Coordinate (base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>)>, duration=2.5 d, rotation_model=howard): (lon, lat) in deg
      (10., 20.)>

A ``RotatedSun*`` frame containing coordinate data can be supplied to ``SkyCoord`` as normal::

  >>> SkyCoord(rs_hgc)
  <SkyCoord (RotatedSunHeliographicCarrington: base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>)>, duration=2.5 d, rotation_model=howard): (lon, lat) in deg
      (10., 20.)>

The above examples used the default differential-rotation model, but any of the models available through :func:`sunpy.physics.differential_rotation.diff_rot` are selectable.
For example, instead of the default ("howard"), one can specify "allen" using the keyword argument ``rotation_model``.
Note the slight difference in the "real" longitude compared to the output above::

  >>> allen = RotatedSunFrame(base=f.HeliographicCarrington(10*u.deg, 20*u.deg, observer="earth",
  ...                                                       obstime="2020-03-04 00:00"),
  ...                         rotated_time="2020-03-06 12:00", rotation_model="allen")
  >>> allen.transform_to(allen.base)
  <HeliographicCarrington Coordinate (obstime=2020-03-04T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, km)
      (45.22266666, 20., 695700.)>

Transforming coordinate arrays
------------------------------
For another transformation example, we define a meridan with a Carrington longitude of 100 degrees, plus 1 day of differential rotation.
Again, the coordinates are already differentially rotated in inertial space; the ``RotatedSun*`` frame allows one to represent the coordinates in a frame *prior* to the differential rotation::

  >>> meridian = RotatedSunFrame([100]*11*u.deg, range(-75, 90, 15)*u.deg,
  ...                            base=f.HeliographicCarrington(observer="earth", obstime="2001-01-01"),
  ...                            duration=1*u.day)
  >>> meridian
  <RotatedSunHeliographicCarrington Coordinate (base=<HeliographicCarrington Frame (obstime=2001-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>)>, duration=1.0 d, rotation_model=howard): (lon, lat) in deg
      [(100., -75.), (100., -60.), (100., -45.), (100., -30.), (100., -15.),
       (100.,   0.), (100.,  15.), (100.,  30.), (100.,  45.), (100.,  60.),
       (100.,  75.)]>

An easy way to "see" the differential rotation is to transform the coordinates to the base coordinate frame.
Note that the points closer to the equator (latitude of 0 degrees) have evolved farther in longitude than the points at high latitudes::

  >>> meridian.transform_to(meridian.base)
  <HeliographicCarrington Coordinate (obstime=2001-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, km)
      [(110.755047, -75., 695700.), (111.706972, -60., 695700.),
       (112.809044, -45., 695700.), (113.682163, -30., 695700.),
       (114.17618 , -15., 695700.), (114.326328,   0., 695700.),
       (114.17618 ,  15., 695700.), (113.682163,  30., 695700.),
       (112.809044,  45., 695700.), (111.706972,  60., 695700.),
       (110.755047,  75., 695700.)]>

.. testsetup::
  # The next test is run with fixed-precision printing to ensure no whitespace appears when tested
  >>> import numpy as np
  >>> old_floatmode = np.get_printoptions()['floatmode']
  >>> np.set_printoptions(floatmode='fixed')

In the specific case of `~sunpy.coordinates.frames.HeliographicCarrington`, this frame rotates with the Sun, but in a non-differential manner.
The Carrington longitude approximately follows the rotation of the Sun.
One can transform to the coordinate frame of 1 day in the future to see the difference between Carrington rotation and differential rotation.
Note that equator rotates slightly faster than the Carrington rotation rate (its longitude is now greater than 100 degrees), but most latitudes rotate slower than the Carrington rotation rate::

  >>> meridian.transform_to(f.HeliographicCarrington(observer="earth", obstime="2001-01-02"))
  <HeliographicCarrington Coordinate (obstime=2001-01-02T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, km)
      [( 96.71777552, -75.1035280, 695509.61226612),
       ( 97.60193088, -60.0954217, 695194.47689542),
       ( 98.68350999, -45.0808511, 694918.44538999),
       ( 99.54760854, -30.0611014, 694697.75301952),
       (100.03737064, -15.0375281, 694544.31380180),
       (100.18622957, -0.01157236, 694467.21969767),
       (100.03737064,  15.0151761, 694471.58239044),
       ( 99.54760854,  30.0410725, 694557.27090716),
       ( 98.68350999,  45.0645144, 694719.82847332),
       ( 97.60193088,  60.0838908, 694951.31065278),
       ( 96.71777552,  75.0975847, 695238.51302901)]>


.. testcleanup::
  >>> np.set_printoptions(floatmode=old_floatmode)

Be aware that transformations with a change in ``obstime`` will also contend with a translation of the center of the Sun.
Note that the ``radius`` component above is no longer precisely on the surface of the Sun.
For precise transformations of solar features, one should also use the context manager :func:`~sunpy.coordinates.transformations.transform_with_sun_center` to account for the translational motion of the Sun.
Using the context manager, the ``radius`` component stays as the solar radius as desired::

  >>> from sunpy.coordinates import transform_with_sun_center
  >>> with transform_with_sun_center():
  ...     print(meridian.transform_to(f.HeliographicCarrington(observer="earth", obstime="2001-01-02")))
  <HeliographicCarrington Coordinate (obstime=2001-01-02T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, km)
      [( 96.570646, -75., 695700.), ( 97.52257 , -60., 695700.),
       ( 98.624643, -45., 695700.), ( 99.497762, -30., 695700.),
       ( 99.991779, -15., 695700.), (100.141927,   0., 695700.),
       ( 99.991779,  15., 695700.), ( 99.497762,  30., 695700.),
       ( 98.624643,  45., 695700.), ( 97.52257 ,  60., 695700.),
       ( 96.570646,  75., 695700.)]>



Transforming multiple durations of rotation
-------------------------------------------

Another common use case for differential rotation is to track a solar feature over a sequence of time steps.
Let's track an active region that starts at `~sunpy.coordinates.frames.Helioprojective` coordinates (-123 arcsec, 456 arcsec), as seen from Earth, and we will look both backwards and forwards in time.
When ``duration`` is an array, the base coordinate will be automatically upgraded to an array if it is a scalar.
We specify a range of durations from -5 days to +5 days, stepping at 1-day increments::

  >>> durations = range(-5, 6, 1)*u.day
  >>> ar_start = f.Helioprojective(-123*u.arcsec, 456*u.arcsec,
  ...                              obstime="2001-01-01", observer="earth")
  >>> ar = RotatedSunFrame(base=ar_start, duration=durations)
  >>> ar
  <RotatedSunHelioprojective Coordinate (base=<Helioprojective Frame (obstime=2001-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>)>, duration=[-5. -4. -3. -2. -1.  0.  1.  2.  3.  4.  5.] d, rotation_model=howard): (Tx, Ty) in arcsec
      [(-123., 456.), (-123., 456.), (-123., 456.), (-123., 456.),
       (-123., 456.), (-123., 456.), (-123., 456.), (-123., 456.),
       (-123., 456.), (-123., 456.), (-123., 456.)]>

Let's convert to the base coordinate frame to reveal the motion of the active region over time::

  >>> ar.transform_to(ar.base)
  <Helioprojective Coordinate (obstime=2001-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      [(-865.549563, 418.102848, 0.982512),
       (-794.67361 , 429.259359, 0.981549),
       (-676.999492, 439.158483, 0.980695),
       (-519.354795, 447.212391, 0.980001),
       (-330.98304 , 452.940564, 0.979507),
       (-123.      , 456.      , 0.979244),
       (  92.27676 , 456.207078, 0.979226),
       ( 302.081349, 453.54936 , 0.979455),
       ( 493.984308, 448.186389, 0.979917),
       ( 656.653862, 440.439434, 0.980585),
       ( 780.541211, 430.770974, 0.981419)]>

Be aware that these coordinates are represented in the `~sunpy.coordinates.frames.Helioprojective` coordinates as seen from Earth at the base time.
Since the Earth moves in its orbit around the Sun, one may be more interested in representing these coordinates as they would been seen by an Earth observer at each time step.
Since the destination frame of the transformation will now have arrays for ``obstime`` and ``observer``, one actually has to construct the initial coordinate with an array for ``obstime`` (and ``observer``) due to a limitation in Astropy.
Note that the active region moves slightly slower across the disk of the Sun because the Earth orbits in the same direction as the Sun rotates, thus reducing the apparent rotation of the Sun::

  >>> ar_start_array = f.Helioprojective([-123]*len(durations)*u.arcsec,
  ...                                    [456]*len(durations)*u.arcsec,
  ...                                    obstime=["2001-01-01"]*len(durations), observer="earth")
  >>> ar_array = RotatedSunFrame(base=ar_start_array, duration=durations)
  >>> earth_hpc = f.Helioprojective(obstime=ar_array.rotated_time, observer="earth")
  >>> ar_array.transform_to(earth_hpc)
  <Helioprojective Coordinate (obstime=['2000-12-27 00:00:00.000' '2000-12-28 00:00:00.000'
   '2000-12-29 00:00:00.000' '2000-12-30 00:00:00.000'
   '2000-12-31 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-02 00:00:00.000' '2001-01-03 00:00:00.000'
   '2001-01-04 00:00:00.000' '2001-01-05 00:00:00.000'
   '2001-01-06 00:00:00.000'], rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      [(-853.35712 , 420.401517, 0.982294),
       (-771.20926 , 429.298481, 0.981392),
       (-650.31062 , 437.85932 , 0.980601),
       (-496.634378, 445.519914, 0.97996 ),
       (-317.863549, 451.731964, 0.9795  ),
       (-123.      , 456.      , 0.979244),
       (  78.103714, 457.916782, 0.979203),
       ( 275.263157, 457.194475, 0.97938 ),
       ( 458.500759, 453.689226, 0.979764),
       ( 618.572111, 447.417202, 0.980336),
       ( 747.448484, 438.560811, 0.981067)]>


Transforming into RotatedSun frames
-----------------------------------

So far, all of the examples show transformations with the ``RotatedSun*`` frame as the starting frame.
The ``RotatedSun*`` frame can also be the destination frame, which can be more intuitive in some situations and even necessary in some others (due to API limitations).
Let's use a coordinate from earlier, which represents the coordinate in a "real" coordinate frame::

  >>> coord = rs_hgc.transform_to(rs_hgc.base)
  >>> coord
  <HeliographicCarrington Coordinate (obstime=2020-03-04T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, km)
      (45.13354448, 20., 695700.)>

If we create a ``RotatedSun*`` frame for a different base time, we can represent that same point using coordinates prior to differential rotation::

  >>> rs_frame = RotatedSunFrame(base=f.HeliographicCarrington(observer="earth",
  ...                                                          obstime=coord.obstime),
  ...                            rotated_time="2020-03-06 12:00")
  >>> rs_frame
  <RotatedSunHeliographicCarrington Frame (base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>)>, duration=2.5 d, rotation_model=howard)>

  >>> new_coord = coord.transform_to(rs_frame)
  >>> new_coord
  <RotatedSunHeliographicCarrington Coordinate (base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>)>, duration=2.5 d, rotation_model=howard): (lon, lat, radius) in (deg, deg, km)
      (10., 20., 695700.)>

There coordinates are stored in the ``RotatedSun*`` frame, but it can be useful to "pop off" this extra layer and retain only the coordinate representation in the base coordinate frame.
There is a convenience method called :meth:`~sunpy.coordinates.metaframes.RotatedSunFrame.as_base()` to do exactly that.
Be aware the resulting coordinate does *not* point to the same location in inertial space, despite the superficial similarity.
Essentially, the component values have been copied from one coordinate frame to a different coordinate frame, and thus this is not merely a transformation between coordinate frames::

  >>> new_coord.as_base()
  <HeliographicCarrington Coordinate (obstime=2020-03-04T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, km)
      (10., 20., 695700.)>

Example uses of RotatedSunFrame
===============================

Here are the examples in our gallery that use `~sunpy.coordinates.metaframes.RotatedSunFrame`:

.. minigallery:: sunpy.coordinates.RotatedSunFrame


..
  >>> # Reset change to default print options
  >>> np.set_printoptions(precision=8)
.. _sunpy-coordinates-velocities:

Coordinates with velocity information
*************************************

Velocity information can be added to any coordinate [#differentials]_.
When the coordinate is transformed to a different coordinate frame, the velocity vector will be transformed as appropriate.
Be aware that the transformation framework does not take any velocity information into account when transforming the position vector.

Creating a SkyCoord with velocity
=================================
Velocity information can be added as keyword arguments to `~astropy.coordinates.SkyCoord`.
For SunPy's frames, the names of the velocities components are the names of the position components prepended by "d\_", e.g.,::

    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> import sunpy.coordinates

    >>> sc = SkyCoord(lon=10*u.deg, lat=20*u.deg, distance=1*u.AU,
    ...               d_lon=3*u.arcmin/u.week, d_lat=4*u.arcmin/u.d, d_distance=5*u.km/u.min,
    ...               frame='heliocentricinertial', obstime='2021-01-01')
    >>> sc
    <SkyCoord (HeliocentricInertial: obstime=2021-01-01T00:00:00.000): (lon, lat, distance) in (deg, deg, AU)
        (10., 20., 1.)
     (d_lon, d_lat, d_distance) in (arcsec / s, arcsec / s, km / s)
        (0.00029762, 0.00277778, 0.08333333)>

See :ref:`astropy-coordinates-velocities` for ways to add velocity information to existing coordinates.

Querying velocity information
=============================
SunPy has functions to query the positions of planets or other objects (e.g., :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`).
For any of these functions, if ``include_velocity=True`` is specified, the returned coordinate will include velocity information, e.g.,::

    >>> from sunpy.coordinates import get_body_heliographic_stonyhurst
    >>> get_body_heliographic_stonyhurst('mercury', '2021-01-01', include_velocity=True)
    <HeliographicStonyhurst Coordinate (obstime=2021-01-01T00:00:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (-156.46460438, -1.38836399, 0.43234904)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (0.09138133, 0.00720229, -7.2513617)>

Transforming velocities
=======================
The transformation of the velocity vector between two coordinate frames takes into account two effects:

* The change in the direction of the velocity vector due to a change in the orientation of the axes between the two frames.
* The "induced" velocity due to the time dependence of the frames themselves.

Orientation change
------------------
To illustrate the orientation change, let's start with the `~astropy.coordinates.SkyCoord` created at the beginning, which was defined in the `~sunpy.coordinates.frames.HeliocentricInertial` frame.
We transform to Astropy's `~astropy.coordinates.HCRS` frame, which is a different inertial [#inertial]_ frame that is also centered at the Sun::

    >>> sc_hcrs = sc.transform_to('hcrs')
    >>> sc_hcrs
    <SkyCoord (HCRS: obstime=2021-01-01T00:00:00.000): (ra, dec, distance) in (deg, deg, AU)
        (80.95428245, 44.31500877, 1.)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (-8828397.4990187, 87659732.13232313, 0.08333338)>

Even though the velocity vectors are oriented very differently in their respective spherical coordinates, their amplitudes are essentially the same::

    >>> sc.velocity.norm()
    <Quantity 2.02654081 km / s>
    >>> sc_hcrs.velocity.norm()
    <Quantity 2.02654083 km / s>

Induced velocity
----------------
To illustrate "induced" velocity, consider the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame, which is defined such that the Earth is always at zero degrees longitude.
That is, this frame rotates around the Sun over time to "follow" the Earth.
Accordingly, the longitude component of Earth's velocity vector will be negligible in this frame::

    >>> from sunpy.coordinates import get_earth
    >>> earth = get_earth('2021-01-01', include_velocity=True)
    >>> earth
    <SkyCoord (HeliographicStonyhurst: obstime=2021-01-01T00:00:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0., -3.02983361, 0.98326486)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (1.82278759e-09, -0.00487486, -0.01720926)>

Transforming this coordinate to the `~sunpy.coordinates.frames.HeliocentricInertial` frame, which does not rotate over time, confirms that the Earth is moving in inertial space at the expected ~1 degree/day in heliographic longitude::

    >>> earth.heliocentricinertial
    <SkyCoord (HeliocentricInertial: obstime=2021-01-01T00:00:00.000): (lon, lat, distance) in (deg, deg, AU)
        (24.55623543, -3.02983361, 0.98326486)
     (d_lon, d_lat, d_distance) in (arcsec / s, arcsec / s, km / s)
        (0.0422321, -0.00487486, -0.01720925)>
    >>> earth.heliocentricinertial.d_lon.to('deg/d')
    <Quantity 1.01357048 deg / d>

Transforming over time
======================
As the transformation framework is currently implemented, transforming between frames with different values of ``obstime`` takes into account any time dependency for the definitions of the frames, but does *not* incorporate any notion of the coordinate itself moving in inertial space.
This behavior does not change even if there is velocity information attached to the coordinate.
For example, if we take the same coordinate created earlier for Earth, and transform it to one day later::

    >>> from sunpy.coordinates import HeliographicStonyhurst
    >>> earth.transform_to(HeliographicStonyhurst(obstime=earth.obstime + 1*u.day))
    <SkyCoord (HeliographicStonyhurst: obstime=2021-01-02T00:00:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (-1.01416251, -3.02979409, 0.98326928)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (-1.19375277e-05, -0.00487485, -0.01743006)>

Note that the location of the Earth in the new frame is ~-1 degree in longitude, as opposed to zero degrees.
That is, this coordinate represents the location of Earth on 2021 January 1 using axes that are defined using the location of Earth on 2021 January 2.

Footnotes
=========

.. [#differentials] Differentials of position with respect to units other than time are also possible, but are not currently well supported.
.. [#inertial] While `~sunpy.coordinates.frames.HeliocentricInertial` and `~astropy.coordinates.HCRS` have fixed axes directions, strictly speaking the small motion of the origin (the Sun) will induce a translational velocity relative to `~astropy.coordinates.ICRS`, but that aspect will cancel out in the transformation.
.. _sunpy-coordinates-carrington:

Calculating Carrington longitude
********************************

Carrington coordinates are a heliographic coordinate system, where the longitude is determined assuming that the Sun has rotated at a sidereal period of ~25.38 days starting at zero longitude at a specific time on 1853 November 9.
For a given observer at a given observation time, the center of the disk of the Sun as seen by that observer (also known as the apparent [#apparent]_ sub-observer point) has a Carrington latitude and longitude.
When the observer is at Earth, the apparent sub-Earth [#subEarth]_ Carrington latitude and longitude are also known as |B0| and |L0|, respectively.

Quick summary
=============
SunPy calculates the sub-observer Carrington longitude in a manner that enables co-alignment of images of the Sun's surface from different observatories at different locations/velocities in the solar system, but that results in a ~20-arcsecond discrepancy with |AA|.

Observer effects
================
An observer will perceive the locations and orientations of solar-system objects different from how they truly are.
There are two primary effects for why the apparent sub-observer point is not the same as the "true" sub-observer point.

Light travel time
-----------------

It takes time for light to travel from the Sun to any observer (e.g., ~500 seconds for an observer at Earth).
Thus, the farther the observer is from the Sun, the farther back in time the Sun will appear to be.
This effect is sometimes called "planetary aberration".

As an additional detail, the observer is closer to the sub-observer point on the Sun's surface than to the center of the Sun (i.e., by the radius of the Sun), which corresponds to ~2.3 lightseconds.

Observer motion
---------------

The motion of an observer will shift the apparent location of objects, and this effect is called "stellar aberration".
A more subtle question is whether observer motion also shifts the apparent sub-observer Carrington longitude.

Stellar aberration causes the apparent ecliptic longitude of an observer to be different from its true ecliptic longitude.
This shift in ecliptic longitude also means a shift in the heliographic longitude of the observer.
For example, for an Earth observer, the apparent heliographic longitude of Earth is ~0.006 degrees (~20 arcseconds) less than its true heliographic longitude.
For some purposes, this shift in the heliographic longitude of the observer should also shift the apparent sub-observer Carrington longitude.

However, shifting the sub-observer Carrington longitude due to observer motion is not always appropriate.
For certain types of observations (e.g., imagery of features on the surface of the Sun), stellar aberration does not shift the positions of solar features **relative** to the center of the disk of the Sun; everything shifts in tandem.
Thus, for the purposes of co-aligning images from different observatories, it is critical to **exclude** any correction for observer motion.

Background information
======================

The IAU definition
------------------

Since the original definition of Carrington coordinates [#Carrington]_, there have been a number of refinements of the definition, but there have been inconsistencies with how those definitions are interpreted.
Possible points of contention include the reference epoch, the rotation period, and the handling of observer effects.

The most recent definition was essentially laid down in 2007 by the IAU [#IAU]_:

* The reference epoch is J2000.0 TDB (2000 January 1 12:00 TDB).
  The difference between Barycentric Dynamical Time (TDB) and Terrestrial Time (TT) is very small, but both are appreciably different from Coordinated Universal Time (UTC).
* The longitude of the prime meridian (|W0|) at the reference epoch is 84.176 degrees.
  This longitude is the "true" longitude, without accounting for any effects that would modify the apparent prime meridian for an observer (e.g., light travel time).
* The sidereal rotation rate is exactly 14.1844 degrees per day.
  This is close to, but not exactly the same as, a sidereal rotation period of 25.38 days or a mean synodic period of 27.2753 days.
* The ICRS coordinates of the Sun's north pole of rotation is given as a right ascension of 286.13 degrees and a declination of 63.87 degrees.

The tuning of |W0|
------------------

Importantly, the IAU parameters (specifically, |W0|) were tuned following an investigation by |AA| [#AA]_.
The prescription of |AA| included both the correction for light travel time and the correction for stellar aberration due to Earth's motion.
When using this prescription, |AA| found that a |W0| value of 84.176 degrees minimized the differences between the modern approach and earlier approaches.
However, for those purposes where one should use sub-observer Carrington longitudes without the stellar-aberration correction, there will be a discrepancy compared to |AA| and to calculations made using earlier approaches.

The approach in SunPy
=====================

SunPy determines apparent sub-observer Carrington longitude by including the correction for the difference in light travel time, but **excluding** any correction for stellar aberration due to observer motion.
The exclusion of a stellar-aberration correction is appropriate for purposes such as image co-alignment, where Carrington longitudes must be consistently associated with features on the Sun's surface.

Comparisons to other sources
============================

Compared to |AA| and older methods of calculation
-------------------------------------------------

|AA| publishes the apparent sub-Earth Carrington longitude (|L0|), and these values include the stellar-aberration correction.
Consequently, SunPy values will be consistently ~0.006 degrees (~20 arcseconds) greater than |AA| values, although these discrepancies are not always apparent at the printed precision (0.01 degrees).

Since |AA| specifically tuned the IAU parameters to minimize the discrepancies with older methods of calulation under their prescription that includes the stellar-aberration correction, SunPy values will also be ~20 arcseconds greater than values calculated using older methods.
Be aware that older methods of calculation may not have accounted for the variable light travel time between the Sun and the Earth, which can cause additional discrepancies of up to ~5 arcseconds.

|AA| does not appear to account for the difference in light travel time between the sub-Earth point on the Sun's surface and the center of the Sun (~2.3 lightseconds), which results in a fixed discrepancy of ~1.4 arcseconds in |L0|.
However, this additional discrepancy is much smaller than the difference in treatment of stellar aberration.

Compared to SPICE
-----------------

In `SPICE <https://naif.jpl.nasa.gov/naif/>`_, the apparent sub-observer Carrington longitude can be calculated using the function `subpnt <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html>`_.
`subpnt`_ allows for two types of aberration correction through the input argument ``abcorr``:

* "LT", which is just the correction for light travel time
* "LT+S", which is the correction for light travel time ("LT") plus the correction for observer motion ("S")

SunPy calculates the apparent sub-observer Carrington longitude in a manner equivalent to specifying "LT" (as opposed to "LT+S").
The discrepancy between SunPy values and SPICE values is no greater than 0.01 arcseconds.

Compared to JPL Horizons
------------------------

In `JPL Horizons <https://ssd.jpl.nasa.gov/?horizons>`_, one can request the "Obs sub-long & sub-lat".
JPL Horizons appears to start from the IAU parameters and to include the correction for light travel time but not the correction for observer motion (i.e., equivalent to specifying "LT" to `subpnt`_ in `SPICE`_).
The discrepancy between SunPy values and JPL Horizons values is no greater than 0.01 arcseconds.

Compared to SunSPICE
--------------------

In `SunSPICE <https://stereo-ssc.nascom.nasa.gov/sunspice.shtml>`_, one can convert to and from the "Carrington" coordinate system using the function `convert_sunspice_coord <https://hesperia.gsfc.nasa.gov/ssw/packages/sunspice/idl/convert_sunspice_coord.pro>`_.
However, these Carrington longitudes are "true" rather than "apparent" because the observer is not specified, so there are no corrections for light travel time or for observer motion.
For example, for an Earth observer, SunPy values will be consistently greater than SunSPICE's coordinate transformation by ~0.82 degrees.

Footnotes
=========

.. [#apparent] The term "apparent" (as opposed to "true") is used to indicate that the observer perceives the locations and orientations of objects different from how they truly are due to effects that include light travel time.
.. [#subEarth] "Sub-Earth" is a special case of "sub-observer", where the observer is at Earth (specifically, Earth center).  That is, the sub-Earth point is the point on the Sun's surface "below" the center of the Earth along the Sun-Earth line.
.. [#Carrington] Carrington (1863), *Observations of the Spots on the Sun*, p. 16 (`<https://archive.org/details/observationsofsp00carr/page/16/mode/1up>`__) and p.27 (`<https://archive.org/details/observationsofsp00carr/page/27/mode/1up>`__)
.. [#IAU] Seidelmann et al. (2007), "Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2006", `<http://dx.doi.org/10.1007/s10569-007-9072-y>`__
.. [#AA] Urban & Kaplan (2007), "Investigation of Change in the Computational Technique of the Sun's Physical Ephemeris in The Astronomical Almanac", `<http://asa.hmnao.com/static/files/sun_rotation_change.pdf>`__

.. |AA| replace:: *The Astronomical Almanac*
.. |B0| replace:: B\ :sub:`0`
.. |L0| replace:: L\ :sub:`0`
.. |W0| replace:: W\ :sub:`0`
.. _sunpy-coordinates-wcs:

Coordinates and WCS
*******************

The `sunpy.coordinates` sub-package provides a mapping between FITS-WCS CTYPE
convention and the coordinate frames as defined in `sunpy.coordinates`. This is
used via the `astropy.wcs.utils.wcs_to_celestial_frame` function, with which the
SunPy frames are registered upon being imported. This list is used by packages
such as ``wcsaxes`` to convert from `astropy.wcs.WCS` objects to coordinate
frames.

The `sunpy.map.GenericMap` class creates `astropy.wcs.WCS` objects as
``amap.wcs``, however, it adds some extra attributes to the `~astropy.wcs.WCS`
object to be able to fully specify the coordinate frame. It adds
``heliographic_observer`` and ``rsun``.

If you want to obtain a un-realized coordinate frame corresponding to a
`~sunpy.map.GenericMap` object you can do the following::

  >>> import sunpy.map
  >>> from sunpy.data.sample import AIA_171_IMAGE  # doctest: +REMOTE_DATA
  >>> amap = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
  >>> amap.observer_coordinate  # doctest: +REMOTE_DATA
    <SkyCoord (HeliographicStonyhurst: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>

which is equivalent to::

  >>> from astropy.wcs.utils import wcs_to_celestial_frame # doctest: +REMOTE_DATA
  >>> wcs_to_celestial_frame(amap.wcs)  # doctest: +REMOTE_DATA
    <Helioprojective Frame (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>)>
.. _sunpy-coordinates:

Coordinates (`sunpy.coordinates`)
*********************************

This sub-package contains:

* A robust framework for working with solar-physics coordinate systems
* Functions to obtain the locations of solar-system bodies (`sunpy.coordinates.ephemeris`)
* Functions to calculate Sun-specific coordinate information (`sunpy.coordinates.sun`)

The SunPy coordinate framework extends the
:ref:`Astropy coordinates framework <astropy:astropy-coordinates>`.

Supported Coordinate Systems
============================

.. list-table::
   :widths: auto
   :header-rows: 1

   * - Coordinate system
     - Abbreviation
     - SunPy/Astropy equivalent
     - Notes
   * - Heliocentric Aries Ecliptic (Mean)
     - HAE (also HEC)
     - Astropy's `~astropy.coordinates.HeliocentricMeanEcliptic`
     -
   * - Heliocentric Cartesian
     - HCC
     - `~sunpy.coordinates.frames.Heliocentric`
     -
   * - Heliocentric Earth Ecliptic
     - HEE
     - `~sunpy.coordinates.frames.HeliocentricEarthEcliptic`
     -
   * - Heliocentric Earth Equatorial
     - HEEQ (also HEQ)
     - `~sunpy.coordinates.frames.HeliographicStonyhurst`
     - Use a Cartesian representation
   * - Heliocentric Inertial
     - HCI
     - `~sunpy.coordinates.frames.HeliocentricInertial`
     -
   * - Heliocentric Radial
     - HCR
     - similar to `~sunpy.coordinates.frames.Heliocentric`
     - Use a cylindrical representation, *but* with a 90-degree offset in ``psi``
   * - Heliocentric/Heliographic Radial-Tangential-Normal
     - HGRTN
     - similar to `~sunpy.coordinates.frames.Heliocentric`
     - The axes are permuted, with HCC X, Y, Z equivalent respectively to HGRTN Y, Z, X
   * - Heliographic Carrington
     - HGC
     - `~sunpy.coordinates.frames.HeliographicCarrington`
     -
   * - Heliographic Stonyhurst
     - HGS
     - `~sunpy.coordinates.frames.HeliographicStonyhurst`
     -
   * - Helioprojective Cartesian
     - HPC
     - `~sunpy.coordinates.frames.Helioprojective`
     -
   * - Geocentric Earth Equatorial (Mean)
     - GEI
     - `~sunpy.coordinates.frames.GeocentricEarthEquatorial`
     -
   * - Geographic
     - GEO
     - Astropy's `~astropy.coordinates.ITRS`
     - The precise geographic definitions may differ
   * - Geocentric Solar Ecliptic
     - GSE
     - `~sunpy.coordinates.frames.GeocentricSolarEcliptic`
     -


For a description of these coordinate systems,
see `Thompson (2006) <https://doi.org/10.1051/0004-6361:20054262>`_
and `Franz & Harper (2002) <https://doi.org/10.1016/S0032-0633(01)00119-2>`_
(and `corrected version <https://www2.mps.mpg.de/homes/fraenz/systems/systems3art/systems3art.html>`_).


Getting Started
===============

The easiest interface to work with coordinates is through the `~astropy.coordinates.SkyCoord` class::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> from sunpy.coordinates import frames

  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame=frames.Helioprojective)
  >>> c = SkyCoord(x=-72241.0*u.km, y=361206.1*u.km, z=589951.4*u.km, frame=frames.Heliocentric)
  >>> c = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
  >>> c
  <SkyCoord (HeliographicStonyhurst: obstime=None, rsun=695700.0 km): (lon, lat) in deg
      (70., -30.)>


It is also possible to use strings to specify the frame but in that case make sure to
explicitly import `sunpy.coordinates` as it registers SunPy's coordinate frames with
the Astropy coordinates framework::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord

  >>> import sunpy.coordinates
  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame='helioprojective', observer='earth')
  >>> c
  <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=earth): (Tx, Ty) in arcsec
      (-100., 500.)>


`~astropy.coordinates.SkyCoord` and all coordinate frames
support array coordinates. These work the same as single-value coordinates,
but they store multiple coordinates in a single object. When you're going to
apply the same operation to many different coordinates, this is a better choice
than a list of `~astropy.coordinates.SkyCoord` objects, because it will be
*much* faster than applying the operation to each
`~astropy.coordinates.SkyCoord` in a ``for`` loop::

   >>> c = SkyCoord([-500, 400]*u.arcsec, [100, 200]*u.arcsec, frame=frames.Helioprojective)
   >>> c
   <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=None): (Tx, Ty) in arcsec
       [(-500.,  100.), ( 400.,  200.)]>
   >>> c[0]
   <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=None): (Tx, Ty) in arcsec
       (-500.,  100.)>


Accessing Coordinates
---------------------

Individual coordinates can be accessed via attributes on the SkyCoord object,
but the names of the components of the coordinates can depend on the the frame and the chosen
representation (e.g., Cartesian versus spherical).

`~sunpy.coordinates.Helioprojective`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the helioprojective frame, the theta_x and theta_y components are accessed as
``Tx`` and ``Ty``, respectively::

  >>> c = SkyCoord(-500*u.arcsec, 100*u.arcsec, frame=frames.Helioprojective)
  >>> c.Tx
  <Longitude -500. arcsec>
  >>> c.Ty
  <Latitude 100. arcsec>

`~sunpy.coordinates.Heliocentric`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Heliocentric is typically used with Cartesian components::

  >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km, frame=frames.Heliocentric)
  >>> c.x
  <Quantity -72241. km>
  >>> c.y
  <Quantity 361206.1 km>
  >>> c.z
  <Quantity 589951.4 km>

`~sunpy.coordinates.HeliographicStonyhurst` and `~sunpy.coordinates.HeliographicCarrington`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both of the heliographic frames have the components of latitude, longitude and radius::

   >>> c = SkyCoord(70*u.deg, -30*u.deg, 1*u.AU, frame=frames.HeliographicStonyhurst)
   >>> c.lat
   <Latitude -30. deg>
   >>> c.lon
   <Longitude 70. deg>
   >>> c.radius
   <Distance 1. AU>

Heliographic Stonyhurst, when used with Cartesian components, is known as Heliocentric
Earth Equatorial (HEEQ).  Here's an example of how to use
`~sunpy.coordinates.frames.HeliographicStonyhurst` for HEEQ coordinates::

  >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km,
  ...              representation_type='cartesian', frame=frames.HeliographicStonyhurst)
  >>> c.x
  <Quantity -72241. km>
  >>> c.y
  <Quantity 361206.1 km>
  >>> c.z
  <Quantity 589951.4 km>

Transforming Between Coordinate Frames
======================================

Both `~astropy.coordinates.SkyCoord` and
`~astropy.coordinates.BaseCoordinateFrame` instances have a
`~astropy.coordinates.SkyCoord.transform_to` method. This can be used to
transform the frame to any other frame, either implemented in SunPy or in
Astropy (see also :ref:`astropy-coordinates-transforming`).
An example of transforming the center of the solar disk to Carrington
coordinates is::

   >>> c = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective, obstime="2017-07-26",
   ...              observer="earth")
   >>> c
   <SkyCoord (Helioprojective: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty) in arcsec
       (0., 0.)>
   >>> c.transform_to(frames.HeliographicCarrington)
   <SkyCoord (HeliographicCarrington: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, AU)
       (283.95956776, 5.31701821, 0.00465047)>

It is also possible to transform to any coordinate system implemented in Astropy. This can be used to find the position of the solar limb in AltAz equatorial coordinates::

    >>> from astropy.coordinates import EarthLocation, AltAz
    >>> time = '2017-07-11 15:00'
    >>> greenbelt = EarthLocation(lat=39.0044*u.deg, lon=-76.8758*u.deg)
    >>> greenbelt_frame = AltAz(obstime=time, location=greenbelt)
    >>> west_limb = SkyCoord(900*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective,
    ...                      observer=greenbelt.get_itrs(greenbelt_frame.obstime), obstime=time)  # doctest: +REMOTE_DATA
    >>> west_limb.transform_to(greenbelt_frame)  # doctest: +REMOTE_DATA
    <SkyCoord (AltAz: obstime=2017-07-11 15:00:00.000, location=(1126916.53031967, -4833386.58391627, 3992696.62211575) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt, distance) in (deg, deg, m)
        (111.40782056, 57.1660434, 1.51859559e+11)>

Observer Location Information
=============================

The `~sunpy.coordinates.frames.Helioprojective`, `~sunpy.coordinates.frames.Heliocentric`
and `~sunpy.coordinates.frames.HeliographicCarrington` frames are defined by the location of
the observer. For example in `~sunpy.coordinates.frames.Helioprojective` the
observer is at the origin of the coordinate system. This information is encoded
in such observer-based frames as the ``observer`` attribute, which is itself an instance of the
`~sunpy.coordinates.frames.HeliographicStonyhurst` frame.
The ``observer`` can be a string for a solar-system body (e.g., "earth" or
"mars"), and
`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` will be used
with the specified ``obstime`` to fully specify the observer location.  If
the observer location is not fully specified, or not present at all, most
transformations cannot be performed.
The location of the observer is automatically populated from meta data when
coordinate frames are created using map.

In the case of `~sunpy.coordinates.frames.HeliographicCarrington`, one can specify ``observer='self'`` to indicate that the coordinate itself should be used as the observer for defining the coordinate frame.

It is possible to convert from a `~sunpy.coordinates.frames.Helioprojective`
frame with one observer location to another
`~sunpy.coordinates.frames.Helioprojective` frame with a different observer
location, by converting through `~sunpy.coordinates.frames.HeliographicStonyhurst`, this
does involve making an assumption of the radius of the Sun to calculate the
position on the solar sphere. The conversion can be performed as follows::

  # Input coordinate
  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth", obstime="2017-07-26", frame=frames.Helioprojective)

  Define a new Helioprojective frame with a different observer.
  >>> import sunpy.coordinates
  >>> hpc_out = sunpy.coordinates.Helioprojective(observer="venus", obstime="2017-07-26")

  Perform the transformation from one to the other.
  >>> hpc2 = hpc1.transform_to(hpc_out)

An example with two maps, named ``aia`` and ``stereo``::

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)  # doctest: +SKIP
  >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)  # doctest: +SKIP


Design of the Coordinates Sub-Package
=====================================

This sub-package works by defining a collection of ``Frames``
(`sunpy.coordinates.frames`), which exists on a transformation graph, where the
transformations between the coordinate frames are then defined and registered
with the transformation graph (`sunpy.coordinates.transformations`). It is also
possible to transform SunPy frames to Astropy frames.

Positions within these ``Frames`` are stored as a ``Representation`` of a
coordinate, a representation being a description of a point in a Cartesian,
spherical or cylindrical system (see :ref:`astropy-coordinates-representations`). A frame
that contains a representation of one or many points is said to have been
'realized'.

For a more in depth look at the design and concepts of the Astropy coordinates
system see :ref:`astropy-coordinates-overview`



Frames and SkyCoord
-------------------

The `~astropy.coordinates.SkyCoord` class is a high level wrapper around the
`astropy.coordinates` sub-package. It provides an easier way to create and transform
coordinates, by using string representations for frames rather than the classes
themselves and some other usability improvements, for more information see the
`~astropy.coordinates.SkyCoord` documentation.

The main advantage provided by `~astropy.coordinates.SkyCoord` is the support it
provides for caching Frame attributes. Frame attributes are extra data specified
with a frame, some examples in `sunpy.coordinates` are ``obstime`` or
``observer`` for observer location. Only the frames where this data is
meaningful have these attributes, i.e. only the Helioprojective frames have
``observer``. However, when you transform into another frame and then back to a
projective frame using `~astropy.coordinates.SkyCoord` it will remember the attributes previously
provided, and repopulate the final frame with them. If you were to do
transformations using the Frames alone this would not happen.

The most important implication for this in `sunpy.coordinates` is the ``rsun``
parameter in the projective frames. If you create a projective frame with a
``rsun`` attribute, if you convert back to a projective frame it will be set
correctly. It should also be noted that, if you create a Heliographic frame and
then transform to a projective frame with an ``rsun`` attribute, it will not
match the ``radius`` coordinate in the Heliographic frame. This is because you may
mean to be describing a point above the defined 'surface' of the Sun.

More Detailed Information
=========================

.. toctree::
   :maxdepth: 1

   carrington
   rotatedsunframe
   velocities
   wcs
   other_api


Reference/API
=============

.. automodapi:: sunpy.coordinates

.. automodapi:: sunpy.coordinates.ephemeris

.. automodapi:: sunpy.coordinates.sun

.. automodapi:: sunpy.coordinates.utils
    :no-inheritance-diagram:


Attribution
===========

Some of this documentation was adapted from Astropy under the terms of the `BSD
License
<https://raw.githubusercontent.com/astropy/astropy/master/LICENSE.rst>`_.

This sub-package was initially developed by Pritish Chakraborty as part of GSOC 2014 and Stuart Mumford.
.. _sunpy-coordinates-other-api:

Reference/API for supporting coordinates modules
************************************************

The parts of the following modules that are useful to a typical user are already imported into the `sunpy.coordinates` namespace.

.. automodapi:: sunpy.coordinates.frames

.. automodapi:: sunpy.coordinates.metaframes

.. automodapi:: sunpy.coordinates.transformations

.. automodapi:: sunpy.coordinates.wcs_utils
Acknowledging or Citing SunPy
=============================

If you use SunPy in your scientific work, we would appreciate citing it in your publications.
The continued growth and development of SunPy is dependent on the community being aware of SunPy.

Citing SunPy in Publications
----------------------------

Please add the following line within your methods, conclusion or acknowledgements sections:

   *This research used version X.Y.Z (software citation) of the SunPy open source
   software package (project citation).*

The project citation should be to the `SunPy paper`_, and the software citation should be the specific `Zenodo DOI`_ for the version used in your work.

.. code:: bibtex

    @ARTICLE{sunpy_community2020,
      doi = {10.3847/1538-4357/ab4f7a},
      url = {https://iopscience.iop.org/article/10.3847/1538-4357/ab4f7a},
      author = {{The SunPy Community} and Barnes, Will T. and Bobra, Monica G. and Christe, Steven D. and Freij, Nabil and Hayes, Laura A. and Ireland, Jack and Mumford, Stuart and Perez-Suarez, David and Ryan, Daniel F. and Shih, Albert Y. and Chanda, Prateek and Glogowski, Kolja and Hewett, Russell and Hughitt, V. Keith and Hill, Andrew and Hiware, Kaustubh and Inglis, Andrew and Kirk, Michael S. F. and Konge, Sudarshan and Mason, James Paul and Maloney, Shane Anthony and Murray, Sophie A. and Panda, Asish and Park, Jongyeob and Pereira, Tiago M. D. and Reardon, Kevin and Savage, Sabrina and Sipőcz, Brigitta M. and Stansby, David and Jain, Yash and Taylor, Garrison and Yadav, Tannmay and Rajul and Dang, Trung Kien},
      title = {The SunPy Project: Open Source Development and Status of the Version 1.0 Core Package},
      journal = {The Astrophysical Journal},
      volume = {890},
      issue = {1},
      pages = {68-},
      publisher = {American Astronomical Society},
      year = {2020}
    }

You can also get this information with ``sunpy.__citation__``.

Other SunPy publications
########################

The SunPy v1.0.8 code was reviewed by `The Journal of Open Source Software (JOSS) <https://joss.theoj.org/papers/10.21105/joss.01832>`__.

.. code:: bibtex

    @ARTICLE{Mumford2020,
      doi = {10.21105/joss.01832},
      url = {https://doi.org/10.21105/joss.01832},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {46},
      pages = {1832},
      author = {Stuart Mumford and Nabil Freij and Steven Christe and Jack Ireland and Florian Mayer and V. Hughitt and Albert Shih and Daniel Ryan and Simon Liedtke and David Pérez-Suárez and Pritish Chakraborty and Vishnunarayan K and Andrew Inglis and Punyaslok Pattnaik and Brigitta Sipőcz and Rishabh Sharma and Andrew Leonard and David Stansby and Russell Hewett and Alex Hamilton and Laura Hayes and Asish Panda and Matt Earnshaw and Nitin Choudhary and Ankit Kumar and Prateek Chanda and Md Haque and Michael Kirk and Michael Mueller and Sudarshan Konge and Rajul Srivastava and Yash Jain and Samuel Bennett and Ankit Baruah and Will Barnes and Michael Charlton and Shane Maloney and Nicky Chorley and Himanshu  and Sanskar Modi and James Mason and Naman9639  and Jose Rozo and Larry Manley and Agneet Chatterjee and John Evans and Michael Malocha and Monica Bobra and Sourav Ghosh and Airmansmith97  and Dominik Stańczak and Ruben De Visscher and Shresth Verma and Ankit Agrawal and Dumindu Buddhika and Swapnil Sharma and Jongyeob Park and Matt Bates and Dhruv Goel and Garrison Taylor and Goran Cetusic and Jacob  and Mateo Inchaurrandieta and Sally Dacie and Sanjeev Dubey and Deepankar Sharma and Erik Bray and Jai Rideout and Serge Zahniy and Tomas Meszaros and Abhigyan Bose and André Chicrala and Ankit  and Chloé Guennou and Daniel D'Avella and Daniel Williams and Jordan Ballew and Nick Murphy and Priyank Lodha and Thomas Robitaille and Yash Krishan and Andrew Hill and Arthur Eigenbrot and Benjamin Mampaey and Bernhard Wiedemann and Carlos Molina and Duygu Keşkek and Ishtyaq Habib and Joseph Letts and Juanjo Bazán and Quinn Arbolante and Reid Gomillion and Yash Kothari and Yash Sharma and Abigail Stevens and Adrian Price-Whelan and Ambar Mehrotra and Arseniy Kustov and Brandon Stone and Trung Dang and Emmanuel Arias and Fionnlagh Dover and Freek Verstringe and Gulshan Kumar and Harsh Mathur and Igor Babuschkin and Jaylen Wimbish and Juan Buitrago-Casas and Kalpesh Krishna and Kaustubh Hiware and Manas Mangaonkar and Matthew Mendero and Mickaël Schoentgen and Norbert Gyenge and Ole Streicher and Rajasekhar Mekala and Rishabh Mishra and Shashank Srikanth and Sarthak Jain and Tannmay Yadav and Tessa Wilkinson and Tiago Pereira and Yudhik Agrawal and Jamescalixto  and Yasintoda  and Sophie Murray},
      title = {SunPy: A Python package for Solar Physics},
      journal = {Journal of Open Source Software}
    }

A paper about version v0.5 of SunPy was published in `Computational Science & Discovery <https://iopscience.iop.org/article/10.1088/1749-4699/8/1/014009>`__.

Acknowledging SunPy in Posters and Talks
----------------------------------------

Please include the `Sunpy logo`_ on the title, conclusion slide, or about page.
For websites please link the image to `sunpy.org`_.
Other versions of the logo are available in the `sunpy-logo repository`_.

.. _SunPy paper: https://doi.org/10.3847/1538-4357/ab4f7a
.. _Sunpy logo: https://github.com/sunpy/sunpy-logo/blob/master/sunpy_logo.svg
.. _sunpy.org: https://sunpy.org/
.. _sunpy-logo repository: https://github.com/sunpy/sunpy-logo/
.. _Zenodo DOI: https://doi.org/10.5281/zenodo.591887
This folder contains color table data. Each color table is in a single .csv files,
with the columns denoting red, green, blue values.
************
sunpy.extern
************

This sub-package contains third-party Python packages/modules that are
required for some of SunPy's core functionality.

In particular, this currently includes for Python:

- `AppDirs`_ (Python 2 and 3 versions): This provides the core config file handling for SunPy's configuration system.
- `distro` _ This provides the functionality of the deprecated `platform.linux_distribution()`, which we need (and more)

To replace any of the other Python modules included in this package, simply remove them and update any imports in SunPy to import the system versions rather than the bundled copies.

.. _AppDirs: https://github.com/ActiveState/appdirs
.. _distro: https://github.com/nir0s/distro.git
