Changelog
=========

Version 0.10.2 (October 16, 2021)
---------------------------------

Small bug-fix release:

- Fix regression in `overlay()` in case no geometries are intersecting (but
  have overlapping total bounds) (#2172).
- Fix regression in `overlay()` with `keep_geom_type=True` in case the
  overlay of two geometries in a GeometryCollection with other geometry types
  (#2177).
- Fix `overlay()` to honor the `keep_geom_type` keyword for the
  `op="differnce"` case (#2164).
- Fix regression in `plot()` with a mapclassify `scheme` in case the
  formatted legend labels have duplicates (#2166).
- Fix a bug in the `explore()` method ignoring the `vmin` and `vmax` keywords
  in case they are set to 0 (#2175).
- Fix `unary_union` to correctly handle a GeoSeries with missing values (#2181).
- Avoid internal deprecation warning in `clip()` (#2179).


Version 0.10.1 (October 8, 2021)
--------------------------------

Small bug-fix release:

- Fix regression in `overlay()` with non-overlapping geometries and a
  non-default `how` (i.e. not "intersection") (#2157).


Version 0.10.0 (October 3, 2021)
--------------------------------

Highlights of this release:

- A new `sjoin_nearest()` method to join based on proximity, with the
  ability to set a maximum search radius (#1865). In addition, the `sindex`
  attribute gained a new method for a "nearest" spatial index query (#1865,
  #2053).
- A new `explore()` method on GeoDataFrame and GeoSeries with native support
  for interactive visualization based on folium / leaflet.js (#1953)
- The `geopandas.sjoin()`/`overlay()`/`clip()` functions are now also
  available as methods on the GeoDataFrame (#2141, #1984, #2150).

New features and improvements:

- Add support for pandas' `value_counts()` method for geometry dtype (#2047).
- The `explode()` method has a new `ignore_index` keyword (consistent with
  pandas' explode method) to reset the index in the result, and a new
  `index_parts` keywords to control whether a cumulative count indexing the
  parts of the exploded multi-geometries should be added (#1871).
- `points_from_xy()` is now available as a GeoSeries method `from_xy` (#1936).
- The `to_file()` method will now attempt to detect the driver (if not
  specified) based on the extension of the provided filename, instead of
  defaulting to ESRI Shapefile (#1609).
- Support for the `storage_options` keyword in `read_parquet()` for
  specifying filesystem-specific options (e.g. for S3) based on fsspec (#2107).
- The read/write functions now support `~` (user home directory) expansion (#1876).
- Support the `convert_dtypes()` method from pandas to preserve the
  GeoDataFrame class (#2115).
- Support WKB values in the hex format in `GeoSeries.from_wkb()` (#2106).
- Update the `estimate_utm_crs()` method to handle crossing the antimeridian
  with pyproj 3.1+ (#2049).
- Improved heuristic to decide how many decimals to show in the repr based on
  whether the CRS is projected or geographic (#1895).
- Switched the default for `geocode()` from GeoCode.Farm to the Photon
  geocoding API (https://photon.komoot.io) (#2007).

Deprecations and compatibility notes:

- The `op=` keyword of `sjoin()` to indicate which spatial predicate to use
  for joining is being deprecated and renamed in favor of a new `predicate=`
  keyword (#1626).
- The `cascaded_union` attribute is deprecated, use `unary_union` instead (#2074).
- Constructing a GeoDataFrame with a duplicated "geometry" column is now
  disallowed. This can also raise an error in the `pd.concat(.., axis=1)`
  function if this results in duplicated active geometry columns (#2046).
- The `explode()` method currently returns a GeoSeries/GeoDataFrame with a
  MultiIndex, with an additional level with indices of the parts of the
  exploded multi-geometries. For consistency with pandas, this will change in
  the future and the new `index_parts` keyword is added to control this.

Bug fixes:

- Fix in the `clip()` function to correctly clip MultiPoints instead of
  leaving them intact when partly outside of the clip bounds (#2148).
- Fix `GeoSeries.isna()` to correctly return a boolean Series in case of an
  empty GeoSeries (#2073).
- Fix the GeoDataFrame constructor to preserve the geometry name when the
  argument is already a GeoDataFrame object (i.e. `GeoDataFrame(gdf)`) (#2138).
- Fix loss of the values' CRS when setting those values as a column
  (`GeoDataFrame.__setitem__`) (#1963)
- Fix in `GeoDataFrame.apply()` to preserve the active geometry column name
  (#1955).
- Fix in `sjoin()` to not ignore the suffixes in case of a right-join
  (`how="right`) (#2065).
- Fix `GeoDataFrame.explode()` with a MultiIndex (#1945).
- Fix the handling of missing values in `to/from_wkb` and `to_from_wkt` (#1891).
- Fix `to_file()` and `to_json()` when DataFrame has duplicate columns to
  raise an error (#1900).
- Fix bug in the colors shown with user-defined classification scheme (#2019).
- Fix handling of the `path_effects` keyword in `plot()` (#2127).
- Fix `GeoDataFrame.explode()` to preserve `attrs` (#1935)

Notes on (optional) dependencies:

- GeoPandas 0.10.0 dropped support for Python 3.6 and pandas 0.24. Further,
  the minimum required versions are numpy 1.18, shapely 1.6, fiona 1.8,
  matplotlib 3.1 and pyproj 2.2.
- Plotting with a classification schema now requires mapclassify version >=
  2.4 (#1737).
- Compatibility fixes for the latest numpy in combination with Shapely 1.7 (#2072)
- Compatibility fixes for the upcoming Shapely 1.8 (#2087).
- Compatibility fixes for the latest PyGEOS (#1872, #2014) and matplotlib
  (colorbar issue, #2066).


Version 0.9.0 (February 28, 2021)
---------------------------------

Many documentation improvements and a restyled and restructured website with
a new logo (#1564, #1579, #1617, #1668, #1731, #1750, #1757, #1759).

New features and improvements:

- The `geopandas.read_file` function now accepts more general
  file-like objects (e.g. `fsspec` open file objects). It will now also
  automatically recognize zipped files (#1535).
- The `GeoDataFrame.plot()` method now provides access to the pandas plotting
  functionality for the non-geometry columns, either using the `kind` keyword
  or the accessor method (e.g. `gdf.plot(kind="bar")` or `gdf.plot.bar()`)
  (#1465).
- New `from_wkt()`, `from_wkb()`, `to_wkt()`, `to_wkb()` methods for
  GeoSeries to construct a GeoSeries from geometries in WKT or WKB
  representation, or to convert a GeoSeries to a pandas Seriew with WKT or WKB
  values (#1710).
- New `GeoSeries.z` attribute to access the z-coordinates of Point geometries
  (similar to the existing `.x` and `.y` attributes) (#1773).
- The `to_crs()` method now handles missing values (#1618).
- Support for pandas' new `.attrs` functionality (#1658).
- The `dissolve()` method now allows dissolving by no column (`by=None`) to
  create a union of all geometries (single-row GeoDataFrame) (#1568).
- New `estimate_utm_crs()` method on GeoSeries/GeoDataFrame to determine the
  UTM CRS based on the bounds (#1646).
- `GeoDataFrame.from_dict()` now accepts `geometry` and `crs` keywords
  (#1619).
- `GeoDataFrame.to_postgis()` and `geopandas.read_postgis()` now supports
  both sqlalchemy engine and connection objects (#1638).
- The `GeoDataFrame.explode()` method now allows exploding based on a
  non-geometry column, using the pandas implementation (#1720).
- Performance improvement in `GeoDataFrame/GeoSeries.explode()` when using
  the PyGEOS backend (#1693).
- The binary operation and predicate methods (eg `intersection()`,
  `intersects()`) have a new `align` keyword which allows optionally not
  aligning on the index before performing the operation with `align=False`
  (#1668).
- The `GeoDataFrame.dissolve()` method now supports all relevant keywords of
  `groupby()`, i.e. the `level`, `sort`, `observed` and `dropna` keywords
  (#1845).
- The `geopandas.overlay()` function now accepts `make_valid=False` to skip
  the step to ensure the input geometries are valid using `buffer(0)` (#1802).
- The `GeoDataFrame.to_json()` method gained a `drop_id` keyword to
  optionally not write the GeoDataFrame's index as the "id" field in the
  resulting JSON (#1637).
- A new `aspect` keyword in the plotting methods to optionally allow retaining
  the original aspect (#1512)
- A new `interval` keyword in the `legend_kwds` group of the `plot()` method
  to control the appearance of the legend labels when using a classification
  scheme (#1605).
- The spatial index of a GeoSeries (accessed with the `sindex` attribute) is
  now stored on the underlying array. This ensures that the spatial index is
  preserved in more operations where possible, and that multiple geometry
  columns of a GeoDataFrame can each have a spatial index (#1444).
- Addition of a `has_sindex` attribute on the GeoSeries/GeoDataFrame to check
  if a spatial index has already been initialized (#1627).
- The `geopandas.testing.assert_geoseries_equal()` and `assert_geodataframe_equal()`
  testing utilities now have a `normalize` keyword (False by default) to
  normalize geometries before comparing for equality (#1826). Those functions
  now also give a more informative error message when failing (#1808).

Deprecations and compatibility notes:

- The `is_ring` attribute currently returns True for Polygons. In the future,
  this will be False (#1631). In addition, start to check it for LineStrings
  and LinearRings (instead of always returning False).
- The deprecated `objects` keyword in the `intersection()` method of the
  `GeoDataFrame/GeoSeries.sindex` spatial index object has been removed
  (#1444).

Bug fixes:

- Fix regression in the `plot()` method raising an error with empty
  geometries (#1702, #1828).
- Fix `geopandas.overlay()` to preserve geometries of the correct type which
  are nested within a GeometryCollection as a result of the overlay
  operation (#1582). In addition, a warning will now be raised if geometries
  of different type are dropped from the result (#1554).
- Fix the repr of an empty GeoSeries to not show spurious warnings (#1673).
- Fix the `.crs` for empty GeoDataFrames (#1560).
- Fix `geopandas.clip` to preserve the correct geometry column name (#1566).
- Fix bug in `plot()` method when using `legend_kwds` with multiple subplots
  (#1583)
- Fix spurious warning with `missing_kwds` keyword of the `plot()` method
  when there are no areas with missing data (#1600).
- Fix the `plot()` method to correctly align values passed to the `column`
  keyword as a pandas Series (#1670).
- Fix bug in plotting MultiPoints when passing values to determine the color
  (#1694)
- The `rename_geometry()` method now raises a more informative error message
  when a duplicate column name is used (#1602).
- Fix `explode()` method to preserve the CRS (#1655)
- Fix the `GeoSeries.apply()` method to again accept the `convert_dtype`
  keyword to be consistent with pandas (#1636).
- Fix `GeoDataFrame.apply()` to preserve the CRS when possible (#1848).
- Fix bug in containment test as `geom in geoseries` (#1753).
- The `shift()` method of a GeoSeries/GeoDataFrame now preserves the CRS
  (#1744).
- The PostGIS IO functionality now quotes table names to ensure it works with
  case-sensitive names (#1825).
- Fix the `GeoSeries` constructor without passing data but only an index (#1798).

Notes on (optional) dependencies:

- GeoPandas 0.9.0 dropped support for Python 3.5. Further, the minimum
  required versions are pandas 0.24, numpy 1.15 and shapely 1.6 and fiona 1.8.
- The `descartes` package is no longer required for plotting polygons. This
  functionality is now included by default in GeoPandas itself, when
  matplotlib is available (#1677).
- Fiona is now only imported when used in `read_file`/`to_file`. This means
  you can now force geopandas to install without fiona installed (although it
  is still a default requirement) (#1775).
- Compatibility with the upcoming Shapely 1.8 (#1659, #1662, #1819).


Version 0.8.2 (January 25, 2021)
--------------------------------

Small bug-fix release for compatibility with PyGEOS 0.9.


Version 0.8.1 (July 15, 2020)
-----------------------------

Small bug-fix release:

- Fix a regression in the `plot()` method when visualizing with a
  JenksCaspallSampled or FisherJenksSampled scheme (#1486).
- Fix spurious warning in `GeoDataFrame.to_postgis` (#1497).
- Fix the un-pickling with `pd.read_pickle` of files written with older
  GeoPandas versions (#1511).


Version 0.8.0 (June 24, 2020)
-----------------------------

**Experimental**: optional use of PyGEOS to speed up spatial operations (#1155).
PyGEOS is a faster alternative for Shapely (being contributed back to a future
version of Shapely), and is used in element-wise spatial operations and for
spatial index in e.g. `sjoin` (#1343, #1401, #1421, #1427, #1428). See the
[installation docs](https://geopandas.readthedocs.io/en/latest/install.html#using-the-optional-pygeos-dependency)
for more info and how to enable it.

New features and improvements:

- IO enhancements:

  - New `GeoDataFrame.to_postgis()` method to write to PostGIS database (#1248).
  - New Apache Parquet and Feather file format support (#1180, #1435)
  - Allow appending to files with `GeoDataFrame.to_file` (#1229).
  - Add support for the `ignore_geometry` keyword in `read_file` to only read
    the attribute data. If set to True, a pandas DataFrame without geometry is
    returned (#1383).
  - `geopandas.read_file` now supports reading from file-like objects (#1329).
  - `GeoDataFrame.to_file` now supports specifying the CRS to write to the file
    (#802). By default it still uses the CRS of the GeoDataFrame.
  - New `chunksize` keyword in `geopandas.read_postgis` to read a query in
    chunks (#1123).

- Improvements related to geometry columns and CRS:

  - Any column of the GeoDataFrame that has a "geometry" dtype is now returned
    as a GeoSeries. This means that when having multiple geometry columns, not
    only the "active" geometry column is returned as a GeoSeries, but also
    accessing another geometry column (`gdf["other_geom_column"]`) gives a
    GeoSeries (#1336).
  - Multiple geometry columns in a GeoDataFrame can now each have a different
    CRS. The global `gdf.crs` attribute continues to returns the CRS of the
    "active" geometry column. The CRS of other geometry columns can be accessed
    from the column itself (eg `gdf["other_geom_column"].crs`) (#1339).
  - New `set_crs()` method on GeoDataFrame/GeoSeries to set the CRS of naive
    geometries (#747).

- Improvements related to plotting:

  - The y-axis is now scaled depending on the center of the plot when using a
    geographic CRS, instead of using an equal aspect ratio (#1290).
  - When passing a column of categorical dtype to the `column=` keyword of the
    GeoDataFrame `plot()`, we now honor all categories and its order (#1483).
    In addition, a new `categories` keyword allows to specify all categories
    and their order otherwise (#1173).
  - For choropleths using a classification scheme (using `scheme=`), the
    `legend_kwds` accept two new keywords to control the formatting of the
    legend: `fmt` with a format string for the bin edges (#1253), and `labels`
    to pass fully custom class labels (#1302).

- New `covers()` and `covered_by()` methods on GeoSeries/GeoDataframe for the
  equivalent spatial predicates (#1460, #1462).
- GeoPandas now warns when using distance-based methods with data in a
  geographic projection (#1378).

Deprecations:

- When constructing a GeoSeries or GeoDataFrame from data that already has a
  CRS, a deprecation warning is raised when both CRS don't match, and in the
  future an error will be raised in such a case. You can use the new `set_crs`
  method to override an existing CRS. See
  [the docs](https://geopandas.readthedocs.io/en/latest/projections.html#projection-for-multiple-geometry-columns).
- The helper functions in the `geopandas.plotting` module are deprecated for
  public usage (#656).
- The `geopandas.io` functions are deprecated, use the top-level `read_file` and
  `to_file` instead (#1407).
- The set operators (`&`, `|`, `^`, `-`) are deprecated, use the
  `intersection()`, `union()`, `symmetric_difference()`, `difference()` methods
  instead (#1255).
- The `sindex` for empty dataframe will in the future return an empty spatial
  index instead of `None` (#1438).
- The `objects` keyword in the `intersection` method of the spatial index
  returned by the `sindex` attribute is deprecated and will be removed in the
  future (#1440).

Bug fixes:

- Fix the `total_bounds()` method to ignore missing and empty geometries (#1312).
- Fix `geopandas.clip` when masking with non-overlapping area resulting in an
  empty GeoDataFrame (#1309, #1365).
- Fix error in `geopandas.sjoin` when joining on an empty geometry column (#1318).
- CRS related fixes: `pandas.concat` preserves CRS when concatenating GeoSeries
  objects (#1340), preserve the CRS in `geopandas.clip` (#1362) and in
  `GeoDataFrame.astype` (#1366).
- Fix bug in `GeoDataFrame.explode()` when 'level_1' is one of the column names
  (#1445).
- Better error message when rtree is not installed (#1425).
- Fix bug in `GeoSeries.equals()` (#1451).
- Fix plotting of multi-part geometries with additional style keywords (#1385).

And we now have a [Code of Conduct](https://github.com/geopandas/geopandas/blob/main/CODE_OF_CONDUCT.md)!

GeoPandas 0.8.0 is the last release to support Python 3.5. The next release
will require Python 3.6, pandas 0.24, numpy 1.15 and shapely 1.6 or higher.


Version 0.7.0 (February 16, 2020)
---------------------------------

Support for Python 2.7 has been dropped. GeoPandas now works with Python >= 3.5.

The important API change of this release is that GeoPandas now requires
PROJ > 6 and pyproj > 2.2, and that the `.crs` attribute of a GeoSeries and
GeoDataFrame no longer stores the CRS information as a proj4 string or dict,
but as a ``pyproj.CRS`` object (#1101).

This gives a better user interface and integrates improvements from pyproj and
PROJ 6, but might also require some changes in your code. Check the
[migration guide](https://geopandas.readthedocs.io/en/latest/projections.html#upgrading-to-geopandas-0-7-with-pyproj-2-2-and-proj-6)
in the documentation.

Other API changes;

- The `GeoDataFrame.to_file` method will now also write the GeoDataFrame index
  to the file, if the index is named and/or non-integer. You can use the
  `index=True/False` keyword to overwrite this default inference (#1059).

New features and improvements:

- A new `geopandas.clip` function to clip a GeoDataFrame to the spatial extent
  of another shape (#1128).
- The `geopandas.overlay` function now works for all geometry types, including
  points and linestrings in addition to polygons (#1110).
- The `plot()` method gained support for missing values (in the column that
  determines the colors). By default it doesn't plot the corresponding
  geometries, but using the new `missing_kwds` argument you can specify how to
  style those geometries (#1156).
- The `plot()` method now also supports plotting GeometryCollection and
  LinearRing objects (#1225).
- Added support for filtering with a geometry or reading a subset of the rows in
  `geopandas.read_file` (#1160).
- Added support for the new nullable integer data type of pandas in
  `GeoDataFrame.to_file` (#1220).

Bug fixes:

- `GeoSeries.reset_index()` now correctly results in a GeoDataFrame instead of DataFrame (#1252).
- Fixed the `geopandas.sjoin` function to handle MultiIndex correctly (#1159).
- Fixed the `geopandas.sjoin` function to preserve the index name of the left GeoDataFrame (#1150).


Version 0.6.3 (February 6, 2020)
---------------------------------

Small bug-fix release:

- Compatibility with Shapely 1.7 and pandas 1.0 (#1244).
- Fix `GeoDataFrame.fillna` to accept non-geometry values again when there are
  no missing values in the geometry column. This should make it easier to fill
  the numerical columns of the GeoDataFrame (#1279).


Version 0.6.2 (November 18, 2019)
---------------------------------

Small bug-fix release fixing a few regressions:

- Fix a regression in passing an array of RRB(A) tuples to the ``.plot()``
  method (#1178, #1211).
- Fix the ``bounds`` and ``total_bounds`` attributes for empty GeoSeries, which
  also fixes the repr of an empty or all-NA GeoSeries (#1184, #1195).
- Fix filtering of a GeoDataFrame to preserve the index type when ending up
  with an empty result (#1190).


Version 0.6.1 (October 12, 2019)
--------------------------------

Small bug-fix release fixing a few regressions:

- Fix `astype` when converting to string with Multi geometries (#1145) or when converting a dataframe without geometries (#1144).
- Fix `GeoSeries.fillna` to accept `np.nan` again (#1149).


Version 0.6.0 (September 27, 2019)
----------------------------------

Important note! This will be the last release to support Python 2.7 (#1031)

API changes:

- A refactor of the internals based on the pandas ExtensionArray interface (#1000). The main user visible changes are:

  - The `.dtype` of a GeoSeries is now a `'geometry'` dtype (and no longer a numpy `object` dtype).
  - The `.values` of a GeoSeries now returns a custom `GeometryArray`, and no longer a numpy array. To get back a numpy array of Shapely scalars, you can convert explicitly using `np.asarray(..)`.

- The `GeoSeries` constructor now raises a warning when passed non-geometry data. Currently the constructor falls back to return a pandas `Series`, but in the future this will raise an error (#1085).
- The missing value handling has been changed to now separate the concepts of missing geometries and empty geometries (#601, 1062). In practice this means that (see [the docs](https://geopandas.readthedocs.io/en/v0.6.0/missing_empty.html) for more details):

  - `GeoSeries.isna` now considers only missing values, and if you want to check for empty geometries, you can use `GeoSeries.is_empty` (`GeoDataFrame.isna` already only looked at missing values).
  - `GeoSeries.dropna` now actually drops missing values (before it didn't drop either missing or empty geometries)
  - `GeoSeries.fillna` only fills missing values (behaviour unchanged).
  - `GeoSeries.align` uses missing values instead of empty geometries by default to fill non-matching index entries.

New features and improvements:

- Addition of a `GeoSeries.affine_transform` method, equivalent of Shapely's function (#1008).
- Addition of a `GeoDataFrame.rename_geometry` method to easily rename the active geometry column (#1053).
- Addition of `geopandas.show_versions()` function, which can be used to give an overview of the installed libraries in bug reports (#899).
- The `legend_kwds` keyword of the `plot()` method can now also be used to specify keywords for the color bar (#1102).
- Performance improvement in the `sjoin()` operation by re-using existing spatial index of the input dataframes, if available (#789).
- Updated documentation to work with latest version of geoplot and contextily (#1044, #1088).
- A new ``geopandas.options`` configuration, with currently a single option to control the display precision of the coordinates (``options.display_precision``). The default is now to show less coordinates (3 for projected and 5 for geographic coordinates), but the default can be overridden with the option.

Bug fixes:

- Also try to use `pysal` instead of `mapclassify` if available (#1082).
- The `GeoDataFrame.astype()` method now correctly returns a `GeoDataFrame` if the geometry column is preserved (#1009).
- The `to_crs` method now uses `always_xy=True` to ensure correct lon/lat order handling for pyproj>=2.2.0 (#1122).
- Fixed passing list-like colors in the `plot()` method in case of "multi" geometries (#1119).
- Fixed the coloring of shapes and colorbar when passing a custom `norm` in the `plot()` method (#1091, #1089).
- Fixed `GeoDataFrame.to_file` to preserve VFS file paths (e.g. when a "s3://" path is specified) (#1124).
- Fixed failing case in ``geopandas.sjoin`` with empty geometries (#1138).


In addition, the minimum required versions of some dependencies have been increased: GeoPandas now requirs pandas >=0.23.4 and matplotlib >=2.0.1 (#1002).


Version 0.5.1 (July 11, 2019)
-----------------------------

- Compatibility with latest mapclassify version 2.1.0 (#1025).

Version 0.5.0 (April 25, 2019)
------------------------------

Improvements:

* Significant performance improvement (around 10x) for `GeoDataFrame.iterfeatures`,
  which also improves `GeoDataFrame.to_file` (#864).
* File IO enhancements based on Fiona 1.8:

    * Support for writing bool dtype (#855) and datetime dtype, if the file format supports it (#728).
    * Support for writing dataframes with multiple geometry types, if the file format allows it (e.g. GeoJSON for all types, or ESRI Shapefile for Polygon+MultiPolygon) (#827, #867, #870).

* Compatibility with pyproj >= 2 (#962).
* A new `geopandas.points_from_xy()` helper function to convert x and y coordinates to Point objects (#896).
* The `buffer` and `interpolate` methods now accept an array-like to specify a variable distance for each geometry (#781).
* Addition of a `relate` method, corresponding to the shapely method that returns the DE-9IM matrix (#853).
* Plotting improvements:

    * Performance improvement in plotting by only flattening the geometries if there are actually 'Multi' geometries (#785).
    * Choropleths: access to all `mapclassify` classification schemes and addition of the `classification_kwds` keyword in the `plot` method to specify options for the scheme (#876).
    * Ability to specify a matplotlib axes object on which to plot the color bar with the `cax` keyword, in order to have more control over the color bar placement (#894).

* Changed the default provider in ``geopandas.tools.geocode`` from Google (now requires an API key) to Geocode.Farm (#907, #975).

Bug fixes:

- Remove the edge in the legend marker (#807).
- Fix the `align` method to preserve the CRS (#829).
- Fix `geopandas.testing.assert_geodataframe_equal` to correctly compare left and right dataframes (#810).
- Fix in choropleth mapping when the values contain missing values (#877).
- Better error message in `sjoin` if the input is not a GeoDataFrame (#842).
- Fix in `read_postgis` to handle nullable (missing) geometries (#856).
- Correctly passing through the `parse_dates` keyword in `read_postgis` to the underlying pandas method (#860).
- Fixed the shape of Antarctica in the included demo dataset 'naturalearth_lowres'
  (by updating to the latest version) (#804).


Version 0.4.1 (March 5, 2019)
-----------------------------

Small bug-fix release for compatibility with the latest Fiona and PySAL
releases:

* Compatibility with Fiona 1.8: fix deprecation warning (#854).
* Compatibility with PySAL 2.0: switched to `mapclassify` instead of `PySAL` as
  dependency for choropleth mapping with the `scheme` keyword (#872).
* Fix for new `overlay` implementation in case the intersection is empty (#800).


Version 0.4.0 (July 15, 2018)
-----------------------------

Improvements:

* Improved `overlay` function (better performance, several incorrect behaviours fixed) (#429)
* Pass keywords to control legend behavior (`legend_kwds`) to `plot` (#434)
* Add basic support for reading remote datasets in `read_file` (#531)
* Pass kwargs for `buffer` operation on GeoSeries (#535)
* Expose all geopy services as options in geocoding (#550)
* Faster write speeds to GeoPackage (#605)
* Permit `read_file` filtering with a bounding box from a GeoDataFrame (#613)
* Set CRS on GeoDataFrame returned by `read_postgis` (#627)
* Permit setting markersize for Point GeoSeries plots with column values (#633)
* Started an example gallery (#463, #690, #717)
* Support for plotting MultiPoints (#683)
* Testing functionality (e.g. `assert_geodataframe_equal`) is now publicly exposed (#707)
* Add `explode` method to GeoDataFrame (similar to the GeoSeries method) (#671)
* Set equal aspect on active axis on multi-axis figures (#718)
* Pass array of values to column argument in `plot` (#770)

Bug fixes:

* Ensure that colorbars are plotted on the correct axis (#523)
* Handle plotting empty GeoDataFrame (#571)
* Save z-dimension when writing files (#652)
* Handle reading empty shapefiles (#653)
* Correct dtype for empty result of spatial operations (#685)
* Fix empty `sjoin` handling for pandas>=0.23 (#762)


Version 0.3.0 (August 29, 2017)
-------------------------------

Improvements:

* Improve plotting performance using ``matplotlib.collections`` (#267)
* Improve default plotting appearance. The defaults now follow the new matplotlib defaults (#318, #502, #510)
* Provide access to x/y coordinates as attributes for Point GeoSeries (#383)
* Make the NYBB dataset available through ``geopandas.datasets`` (#384)
* Enable ``sjoin`` on non-integer-index GeoDataFrames (#422)
* Add ``cx`` indexer to GeoDataFrame (#482)
* ``GeoDataFrame.from_features`` now also accepts a Feature Collection (#225, #507)
* Use index label instead of integer id in output of ``iterfeatures`` and
  ``to_json`` (#421)
* Return empty data frame rather than raising an error when performing a spatial join with non overlapping geodataframes (#335)

Bug fixes:

* Compatibility with shapely 1.6.0 (#512)
* Fix ``fiona.filter`` results when bbox is not None (#372)
* Fix ``dissolve`` to retain CRS (#389)
* Fix ``cx`` behavior when using index of 0 (#478)
* Fix display of lower bin in legend label of choropleth plots using a PySAL scheme (#450)


Version 0.2.0
-------------

Improvements:

* Complete overhaul of the documentation
* Addition of ``overlay`` to perform spatial overlays with polygons (#142)
* Addition of ``sjoin`` to perform spatial joins (#115, #145, #188)
* Addition of ``__geo_interface__`` that returns a python data structure
  to represent the ``GeoSeries`` as a GeoJSON-like ``FeatureCollection`` (#116)
  and ``iterfeatures`` method (#178)
* Addition of the ``explode`` (#146) and ``dissolve`` (#310, #311) methods.
* Addition of the ``sindex`` attribute, a Spatial Index using the optional
  dependency ``rtree`` (``libspatialindex``) that can be used to speed up
  certain operations such as overlays (#140, #141).
* Addition of the ``GeoSeries.cx`` coordinate indexer to slice a GeoSeries based
  on a bounding box of the coordinates (#55).
* Improvements to plotting: ability to specify edge colors (#173), support for
  the ``vmin``, ``vmax``, ``figsize``, ``linewidth`` keywords (#207), legends
  for chloropleth plots (#210), color points by specifying a colormap (#186) or
  a single color (#238).
* Larger flexibility of ``to_crs``, accepting both dicts and proj strings (#289)
* Addition of embedded example data, accessible through
  ``geopandas.datasets.get_path``.

API changes:

* In the ``plot`` method, the ``axes`` keyword is renamed to ``ax`` for
  consistency with pandas, and the ``colormap`` keyword is renamed to ``cmap``
  for consistency with matplotlib (#208, #228, #240).

Bug fixes:

* Properly handle rows with missing geometries (#139, #193).
* Fix ``GeoSeries.to_json`` (#263).
* Correctly serialize metadata when pickling (#199, #206).
* Fix ``merge`` and ``concat`` to return correct GeoDataFrame (#247, #320, #322).
GeoPandas [![Actions Status](https://github.com/geopandas/geopandas/workflows/Tests/badge.svg)](https://github.com/geopandas/geopandas/actions?query=workflow%3ATests) [![Coverage Status](https://codecov.io/gh/geopandas/geopandas/branch/main/graph/badge.svg)](https://codecov.io/gh/geopandas/geopandas) [![Join the chat at https://gitter.im/geopandas/geopandas](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/geopandas/geopandas?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/geopandas/geopandas/main) [![DOI](https://zenodo.org/badge/11002815.svg)](https://zenodo.org/badge/latestdoi/11002815)
=========

Python tools for geographic data

Introduction
------------

GeoPandas is a project to add support for geographic data to
[pandas](http://pandas.pydata.org) objects.  It currently implements
`GeoSeries` and `GeoDataFrame` types which are subclasses of
`pandas.Series` and `pandas.DataFrame` respectively.  GeoPandas
objects can act on [shapely](http://shapely.readthedocs.io/en/latest/)
geometry objects and perform geometric operations.

GeoPandas geometry operations are cartesian.  The coordinate reference
system (crs) can be stored as an attribute on an object, and is
automatically set when loading from a file.  Objects may be
transformed to new coordinate systems with the `to_crs()` method.
There is currently no enforcement of like coordinates for operations,
but that may change in the future.

Documentation is available at [geopandas.org](http://geopandas.org)
(current release) and
[Read the Docs](http://geopandas.readthedocs.io/en/latest/)
(release and development versions).

Install
--------

See the [installation docs](https://geopandas.readthedocs.io/en/latest/install.html)
for all details. GeoPandas depends on the following packages:

- ``pandas``
- ``shapely``
- ``fiona``
- ``pyproj``
- ``packaging``

Further, ``matplotlib`` is an optional dependency, required
for plotting, and [``rtree``](https://github.com/Toblerity/rtree) is an optional
dependency, required for spatial joins. ``rtree`` requires the C library [``libspatialindex``](https://github.com/libspatialindex/libspatialindex).

Those packages depend on several low-level libraries for geospatial analysis, which can be a challenge to install. Therefore, we recommend to install GeoPandas using the [conda package manager](https://conda.io/en/latest/). See the [installation docs](https://geopandas.readthedocs.io/en/latest/install.html) for more details.


Get in touch
------------

- Ask usage questions ("How do I?") on [StackOverflow](https://stackoverflow.com/questions/tagged/geopandas) or [GIS StackExchange](https://gis.stackexchange.com/questions/tagged/geopandas).
- Report bugs, suggest features or view the source code [on GitHub](https://github.com/geopandas/geopandas).
- For a quick question about a bug report or feature request, or Pull Request, head over to the [gitter channel](https://gitter.im/geopandas/geopandas).
- For less well defined questions or ideas, or to announce other projects of interest to GeoPandas users, ... use the [mailing list](https://groups.google.com/forum/#!forum/geopandas).


Examples
--------

    >>> import geopandas
    >>> from shapely.geometry import Polygon
    >>> p1 = Polygon([(0, 0), (1, 0), (1, 1)])
    >>> p2 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    >>> p3 = Polygon([(2, 0), (3, 0), (3, 1), (2, 1)])
    >>> g = geopandas.GeoSeries([p1, p2, p3])
    >>> g
    0         POLYGON ((0 0, 1 0, 1 1, 0 0))
    1    POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))
    2    POLYGON ((2 0, 3 0, 3 1, 2 1, 2 0))
    dtype: geometry

![Example 1](doc/source/gallery/test.png)

Some geographic operations return normal pandas object.  The `area` property of a `GeoSeries` will return a `pandas.Series` containing the area of each item in the `GeoSeries`:

    >>> print(g.area)
    0    0.5
    1    1.0
    2    1.0
    dtype: float64

Other operations return GeoPandas objects:

    >>> g.buffer(0.5)
    0    POLYGON ((-0.3535533905932737 0.35355339059327...
    1    POLYGON ((-0.5 0, -0.5 1, -0.4975923633360985 ...
    2    POLYGON ((1.5 0, 1.5 1, 1.502407636663901 1.04...
    dtype: geometry

![Example 2](doc/source/gallery/test_buffer.png)

GeoPandas objects also know how to plot themselves. GeoPandas uses
[matplotlib](http://matplotlib.org) for plotting. To generate a plot of our
GeoSeries, use:

    >>> g.plot()

GeoPandas also implements alternate constructors that can read any data format recognized by [fiona](http://fiona.readthedocs.io/en/latest/). To read a zip file containing an ESRI shapefile with the [boroughs boundaries of New York City](https://data.cityofnewyork.us/City-Government/Borough-Boundaries/tqmj-j8zm) (GeoPandas includes this as an example dataset):

    >>> nybb_path = geopandas.datasets.get_path('nybb')
    >>> boros = geopandas.read_file(nybb_path)
    >>> boros.set_index('BoroCode', inplace=True)
    >>> boros.sort_index(inplace=True)
    >>> boros
                   BoroName     Shape_Leng    Shape_Area  \
    BoroCode
    1             Manhattan  359299.096471  6.364715e+08
    2                 Bronx  464392.991824  1.186925e+09
    3              Brooklyn  741080.523166  1.937479e+09
    4                Queens  896344.047763  3.045213e+09
    5         Staten Island  330470.010332  1.623820e+09

                                                       geometry
    BoroCode
    1         MULTIPOLYGON (((981219.0557861328 188655.31579...
    2         MULTIPOLYGON (((1012821.805786133 229228.26458...
    3         MULTIPOLYGON (((1021176.479003906 151374.79699...
    4         MULTIPOLYGON (((1029606.076599121 156073.81420...
    5         MULTIPOLYGON (((970217.0223999023 145643.33221...

![New York City boroughs](doc/source/gallery/nyc.png)

    >>> boros['geometry'].convex_hull
    BoroCode
    1    POLYGON ((977855.4451904297 188082.3223876953,...
    2    POLYGON ((1017949.977600098 225426.8845825195,...
    3    POLYGON ((988872.8212280273 146772.0317993164,...
    4    POLYGON ((1000721.531799316 136681.776184082, ...
    5    POLYGON ((915517.6877458114 120121.8812543372,...
    dtype: geometry

![Convex hulls of New York City boroughs](doc/source/gallery/nyc_hull.png)
# GeoPandas Project Code of Conduct

Behind the GeoPandas Project is an engaged and respectful community made up of people
from all over the world and with a wide range of backgrounds.
Naturally, this implies diversity of ideas and perspectives on often complex
problems. Disagreement and healthy discussion of conflicting viewpoints is
welcome: the best solutions to hard problems rarely come from a single angle.
But disagreement is not an excuse for aggression: humans tend to take
disagreement personally and easily drift into behavior that ultimately degrades
a community. This is particularly acute with online communication across
language and cultural gaps, where many cues of human behavior are unavailable.
We are outlining here a set of principles and processes to support a
healthy community in the face of these challenges.

Fundamentally, we are committed to fostering a productive, harassment-free
environment for everyone. Rather than considering this code an exhaustive list
of things that you can’t do, take it in the spirit it is intended - a guide to
make it easier to enrich all of us and the communities in which we participate.

Importantly: as a member of our community, *you are also a steward of these
values*.  Not all problems need to be resolved via formal processes, and often
a quick, friendly but clear word on an online forum or in person can help
resolve a misunderstanding and de-escalate things.

However, sometimes these informal processes may be inadequate: they fail to
work, there is urgency or risk to someone, nobody is intervening publicly and
you don't feel comfortable speaking in public, etc.  For these or other
reasons, structured follow-up may be necessary and here we provide the means
for that: we welcome reports by emailing
[*geopandas-conduct@googlegroups.com*](mailto:geopandas-conduct@googlegroups.com)
or by filling out
[this form](https://docs.google.com/forms/d/e/1FAIpQLSd8Tbi2zNl1i2N9COX0yavHEqTGFIPQ1_cLcy1A3JgVc1OrAQ/viewform).

This code applies equally to founders, developers, mentors and new community
members, in all spaces managed by the GeoPandas Project. This
includes the mailing lists, our GitHub organization, our chat room, in-person
events, and any other forums created by the project team. In addition,
violations of this code outside these spaces may affect a person's ability to
participate within them.

By embracing the following principles, guidelines and actions to follow or
avoid, you will help us make the GeoPandas Project a welcoming and productive community. Feel
free to contact the Code of Conduct Committee at
[*geopandas-conduct@googlegroups.com*](mailto:geopandas-conduct@googlegroups.com)  with any questions.


1. **Be friendly and patient**.

2. **Be welcoming**. We strive to be a community that welcomes and supports
   people of all backgrounds and identities. This includes, but is not limited
   to, members of any race, ethnicity, culture, national origin, color,
   immigration status, social and economic class, educational level, sex, sexual
   orientation, gender identity and expression, age, physical appearance, family
   status, technological or professional choices, academic
   discipline, religion, mental ability, and physical ability.

3. **Be considerate**. Your work will be used by other people, and you in turn
   will depend on the work of others. Any decision you take will affect users
   and colleagues, and you should take those consequences into account when
   making decisions. Remember that we're a world-wide community. You may be
   communicating with someone with a different primary language or cultural
   background.

4. **Be respectful**. Not all of us will agree all the time, but disagreement is
   no excuse for poor behavior or poor manners. We might all experience some
   frustration now and then, but we cannot allow that frustration to turn into a
   personal attack. It’s important to remember that a community where people
   feel uncomfortable or threatened is not a productive one.

5. **Be careful in the words that you choose**. Be kind to others. Do not insult
   or put down other community members. Harassment and other exclusionary
   behavior are not acceptable. This includes, but is not limited to:
   * Violent threats or violent language directed against another person
   * Discriminatory jokes and language
   * Posting sexually explicit or violent material
   * Posting (or threatening to post) other people's personally identifying
     information ("doxing")
   * Personal insults, especially those using racist, sexist, and xenophobic terms
   * Unwelcome sexual attention
   * Advocating for, or encouraging, any of the above behavior
   * Repeated harassment of others. In general, if someone asks you to stop,
     then stop

6. **Moderate your expectations**. Please respect that community members choose
   how they spend their time in the project. A thoughtful question about your
   expectations is preferable to demands for another person's time.

7. **When we disagree, try to understand why**. Disagreements, both social and
   technical, happen all the time and the GeoPandas Project is no exception.  Try to
   understand where others are coming from, as seeing a question from their
   viewpoint may help find a new path forward.  And don’t forget that it is
   human to err: blaming each other doesn’t get us anywhere, while we can learn
   from mistakes to find better solutions.

8. **A simple apology can go a long way**. It can often de-escalate a situation,
   and telling someone that you are sorry is an act of empathy that doesn’t
   automatically imply an admission of guilt.


## Reporting

If you believe someone is violating the code of conduct, please report this in
a timely manner. Code of conduct violations reduce the value of the community
for everyone and we take them seriously.

You can file a report by emailing
[*geopandas-conduct@googlegroups.com*](mailto:geopandas-conduct@googlegroups.com) or by filing out
[this form](https://docs.google.com/forms/d/e/1FAIpQLSd8Tbi2zNl1i2N9COX0yavHEqTGFIPQ1_cLcy1A3JgVc1OrAQ/viewform).

The online form gives you the option to keep your report anonymous or request
that we follow up with you directly. While we cannot follow up on an anonymous
report, we will take appropriate action.

Messages sent to the e-mail address or through the form will be sent
only to the Code of Conduct Committee, which currently consists of:

* Hannah Aizenman
* Joris Van den Bossche
* Martin Fleischmann


## Enforcement

Enforcement procedures within the GeoPandas Project follow Project Jupyter's
[*Enforcement Manual*](https://github.com/jupyter/governance/blob/master/conduct/enforcement.md). For information on enforcement, please view the [original manual](https://github.com/jupyter/governance/blob/master/conduct/enforcement.md).

Original text courtesy of the [*Speak
Up!*](http://web.archive.org/web/20141109123859/http://speakup.io/coc.html),
[*Django*](https://www.djangoproject.com/conduct) and [*Jupyter*](https://github.com/jupyter/governance/blob/master/conduct/code_of_conduct.md) Projects,
modified by the GeoPandas Project. We are grateful to those projects for contributing these materials under open licensing terms for us to easily reuse.

All content on this page is licensed under a [*Creative Commons
Attribution*](http://creativecommons.org/licenses/by/3.0/) license.
Guidelines
==========

Contributions to GeoPandas are very welcome. They are likely to
be accepted more quickly if they follow these guidelines.

At this stage of GeoPandas development, the priorities are to define a
simple, usable, and stable API and to have clean, maintainable,
readable code. Performance matters, but not at the expense of those
goals.

In general, GeoPandas follows the conventions of the pandas project
where applicable. Please read the [contributing
guidelines](https://geopandas.readthedocs.io/en/latest/community/contributing.html).


In particular, when submitting a pull request:

- Install the requirements for the development environment (one can do this
  with either conda, and the environment.yml file, or pip, and the
  requirements-dev.txt file, and can use the pandas contributing guidelines
  as a guide).
- All existing tests should pass. Please make sure that the test
  suite passes, both locally and on
  [GitHub Actions](https://github.com/geopandas/geopandas/actions). Status on
  GHA will be visible on a pull request. GHA are automatically enabled
  on your own fork as well. To trigger a check, make a PR to your own fork.

- New functionality should include tests. Please write reasonable
  tests for your code and make sure that they pass on your pull request.

- Classes, methods, functions, etc. should have docstrings. The first
  line of a docstring should be a standalone summary. Parameters and
  return values should be documented explicitly.

Improving the documentation and testing for code already in GeoPandas
is a great way to get started if you'd like to make a contribution.

Style
-----

- GeoPandas supports Python 3.7+ only. The last version of GeoPandas
  supporting Python 2 is 0.6.

- GeoPandas follows [the PEP 8
  standard](http://www.python.org/dev/peps/pep-0008/) and uses
  [Black](https://black.readthedocs.io/en/stable/) and
  [Flake8](http://flake8.pycqa.org/en/latest/) to ensure a consistent
  code format throughout the project.

- Imports should be grouped with standard library imports first,
  third-party libraries next, and GeoPandas imports third. Within each
  grouping, imports should be alphabetized. Always use absolute
  imports when possible, and explicit relative imports for local
  imports when necessary in tests.

- You can set up [pre-commit hooks](https://pre-commit.com/) to
  automatically run `black` and `flake8` when you make a git
  commit. This can be done by installing `pre-commit`:

    $ python -m pip install pre-commit

  From the root of the geopandas repository, you should then install
  `pre-commit`:

    $ pre-commit install

  Then `black` and `flake8` will be run automatically each time you
  commit changes. You can skip these checks with `git commit
  --no-verify`.
---

name: Bug Report
about: Create a bug report to help us improve geopandas
title: "BUG:"
labels: "bug, needs triage"

---

- [ ] I have checked that this issue has not already been reported.

- [ ] I have confirmed this bug exists on the latest version of geopandas.

- [ ] (optional) I have confirmed this bug exists on the main branch of geopandas.

---

**Note**: Please read [this guide](https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports) detailing how to provide the necessary information for us to reproduce your bug.

#### Code Sample, a copy-pastable example

```python
# Your code here

```

#### Problem description

[this should explain **why** the current behaviour is a problem and why the expected output is a better solution]

#### Expected Output

#### Output of ``geopandas.show_versions()``

<details>

[paste the output of ``geopandas.show_versions()`` here leaving a blank line after the details tag]

</details>
---

name: Submit Question
about: Ask a general question about geopandas
title: "QST:"
labels: "question"

---

- [ ] I have searched the [geopandas] tag on [StackOverflow](https://stackoverflow.com/questions/tagged/geopandas) and [GIS StackExchange](https://gis.stackexchange.com/questions/tagged/geopandas) for similar questions.

- [ ] I have asked my usage related question on [StackOverflow](https://stackoverflow.com) or [GIS StackExhange](https://gis.stackexchange.com).

---

#### Question about geopandas

**Note**: If you'd still like to submit a question, please read [this guide](
https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports) detailing how to provide the necessary information for us to reproduce your question.

```python
# Your code here, if applicable

```
---

name: Feature Request
about: Suggest an idea for geopandas
title: "ENH:"
labels: "enhancement"

---

#### Is your feature request related to a problem?

[this should provide a description of what the problem is, e.g. "I wish I could use geopandas to do [...]"]

#### Describe the solution you'd like

[this should provide a description of the feature request, e.g. "`GeoDataFrame.foo` should get a new parameter `bar` that [...]", try to write a docstring for the desired feature]

#### API breaking implications

[this should provide a description of how this feature will affect the API]

#### Describe alternatives you've considered

[this should provide a description of any alternative solutions or features you've considered]

#### Additional context

[add any other context, code examples, or references to existing implementations about the feature request here]

```python
# Your code here, if applicable

```

---

name: Installation Issue
about: Ask about installing geopandas
title: ""
labels: "installation"

---

- [ ] I have read the [documentation on installation](https://geopandas.org/install.html) and followed the instructions provided.

- [ ] I have looked through [issues labeled "installation"](https://github.com/geopandas/geopandas/labels/installation) in the geopandas repo.

- If your issue is related to installation using the `conda-forge` channel, please open an issue in [geopandas-feedstock repository](https://github.com/conda-forge/geopandas-feedstock) instead.

---

#### System information

[what operating system do you have and what package management system are you
using]

#### Environment details

<details>

[if using conda, paste the output of `conda info` and `conda list`; if using
pip, `pip freeze`]

</details>

# Examples Gallery

Examples are available in the [documentation](https://geopandas.readthedocs.io/en/latest/gallery/index.html). Source Jupyter notebooks are in [`doc/source/gallery`](https://github.com/geopandas/geopandas/tree/main/doc/source/gallery).
# Datasets included with geopandas

- `'nybb'`: Borough boundaries of New York City, data provided Department of City Planning (DCP), https://data.cityofnewyork.us/City-Government/Borough-Boundaries/tqmj-j8zm
- `'naturalearth_cities'`: capital cities, based on http://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-populated-places/
- `'naturalearth_lowres'`: country boundaries, based on http://www.naturalearthdata.com/downloads/110m-cultural-vectors/110m-admin-0-countries/


## Data sources and licenses

### `nybb`

#### Summary

These districts were created by the Department of City Planning to aid city agencies in administering public services.

#### Description

The borough boundaries of New York City clipped to the shoreline at mean high tide.

#### Credits

Department of City Planning.

#### Use limitations

This dataset is being provided by the Department of City Planning (DCP) on DCP’s website for informational purposes only. DCP does not warranty the completeness, accuracy, content, or fitness for any particular purpose or use of the dataset, nor are any such warranties to be implied or inferred with respect to the dataset as furnished on the website. DCP and the City are not liable for any deficiencies in the completeness, accuracy, content, or fitness for any particular purpose or use the dataset, or applications utilizing the dataset, provided by any third party.


# About GeoPandas

```{toctree}
:maxdepth: 2
:caption: About
:hidden:

Team <about/team>
Citing <about/citing>
Logo <about/logo>
```


GeoPandas is an open source project to add support for geographic data to pandas objects. It
currently implements `GeoSeries` and `GeoDataFrame` types which are subclasses of
`pandas.Series` and `pandas.DataFrame` respectively. GeoPandas objects can act on
`shapely` geometry objects and perform geometric operations.

GeoPandas is a community-led project written, used and supported by a wide range of
people from all around of world of a large variety of backgrounds. Want to get involved
in the community? See our [community guidelines](community).

GeoPandas will always be 100% open source software, free for all to use and released
under the liberal terms of the BSD-3-Clause license.

```{container} button

{doc}`Team <about/team>`
{doc}`Citing <about/citing>` {doc}`Logo <about/logo>`
```

## Project history

Kelsey Jordahl founded GeoPandas project in 2013 during the Scipy Conference and
released a version 0.1.0 in July 2014. In 2016, Joris Van den Bossche took the lead and
became the maintainer of the project. Since the beginning, GeoPandas is a BSD-licensed
open-source project supported by a [community of
contributors](https://github.com/geopandas/geopandas/graphs/contributors) from around
the world and is now maintained by a [team](about/team) of core developers.

In 2020 GeoPandas became [NumFOCUS Affiliated
Project](https://numfocus.org/sponsored-projects/affiliated-projects) and received two
[Small Development Grants](https://numfocus.org/programs/sustainability) to support its
development.

## Timeline

- **2013**: Beginning of the development
- **2014**: GeoPandas 0.1.0 released
- **2020**: GeoPandas became [NumFOCUS Affiliated
  Project](https://numfocus.org/sponsored-projects/affiliated-projects)


# Getting Started

```{toctree}
---
maxdepth: 2
caption: Getting Started
hidden:
---

Installation <getting_started/install>
Introduction to GeoPandas <getting_started/introduction>
Examples Gallery <gallery/index>
```

## Installation

GeoPandas is written in pure Python, but has several dependencies written in C 
([GEOS](https://geos.osgeo.org), [GDAL](https://www.gdal.org/), [PROJ](https://proj.org/)).  Those base C libraries can sometimes be a challenge to 
install. Therefore, we advise you to closely follow the recommendations below to avoid 
installation problems.

### Easy way

The best way to install GeoPandas is using ``conda`` and ``conda-forge`` channel:

```
conda install -c conda-forge geopandas
```

### Detailed instructions

Do you prefer ``pip install`` or installation from source? Or specific version? See 
{doc}`detailed instructions <getting_started/install>`.

### What now?

- If you don't have GeoPandas yet, check {doc}`Installation <getting_started/install>`. 
- If you have never used GeoPandas and want to get familiar with it and its core 
  functionality quickly, see {doc}`Getting Started Tutorial <getting_started/introduction>`. 
- Detailed illustration how to work with different parts of GeoPandas, how to make maps,
  manage projections, spatially merge data or geocode are part of our 
  {doc}`User Guide <docs/user_guide>`. 
- And if you are interested in the complete 
  documentation of all classes, functions, method and attributes GeoPandas offers, 
  {doc}`API Reference <docs/reference>` is here for you.


```{container} button

{doc}`Installation <getting_started/install>` {doc}`Introduction <getting_started/introduction>`
{doc}`User Guide <docs/user_guide>` {doc}`API Reference <docs/reference>`
```

## Get in touch

Haven't found what you were looking for?

- Ask usage questions ("How do I?") on [StackOverflow](https://stackoverflow.com/questions/tagged/geopandas) or [GIS StackExchange](https://gis.stackexchange.com/questions/tagged/geopandas).
- Report bugs, suggest features or view the source code on [GitHub](https://github.com/geopandas/geopandas).
- For a quick question about a bug report or feature request, or Pull Request,
  head over to the [gitter channel](https://gitter.im/geopandas/geopandas).
- For less well defined questions or ideas, or to announce other projects of
  interest to GeoPandas users, ... use the [mailing list](https://groups.google.com/forum/#!forum/geopandas).

# GeoPandas logo

GeoPandas project uses a logo derived from [`pandas` logo](https://pandas.pydata.org/about/citing.html), enclosing it in a globe illustrating the geographic nature of our data.

## Versions

We have four versions of our logo:

### Primary logo

The primary logo should be used in a majority of cases. Inverted logo or icon should be used only when necessary.

```{image} ../_static/logo/geopandas_logo.png
:alt: geopandas-logo
:align: center
```

### Inverted colors

If you want to place the GeoPandas logo on a dark background, use the inverted version.

```{image} ../_static/logo/geopandas_logo_green.png
:alt: geopandas-logo-green
:align: center
```

### Icon

Although it is possible to use icon independently, we would prefer using the complete variant above.

```{image} ../_static/logo/geopandas_icon.png
:alt: geopandas-icon
:width: 25%
:align: center
```

### Inverted icon

```{image} ../_static/logo/geopandas_icon_green.png
:alt: geopandas-icon-green
:width: 25%
:align: center
```

## Download

You can download all version in SVG and PNG from [GitHub repository](https://github.com/geopandas/geopandas/tree/main/doc/source/_static/logo).


## Colors

Pink and yellow accent colors are shared with `pandas`.

### Green
```{raw} html
<svg xmlns="http://www.w3.org/2000/svg" width="75" height="75" style="float: left">
    <circle cx="33" cy="33" r="33" fill="#139C5A"></circle>
</svg>
```
**HEX:** #139C5A

**RGB:** (19, 156, 90)

### Yellow
```{raw} html
<svg xmlns="http://www.w3.org/2000/svg" width="75" height="75" style="float: left">
    <circle cx="33" cy="33" r="33" fill="#FFCA00"></circle>
</svg>
```
**HEX:** #FFCA00

**RGB:** (255, 202, 0)

### Pink
```{raw} html
<svg xmlns="http://www.w3.org/2000/svg" width="75" height="75" style="float: left">
    <circle cx="33" cy="33" r="33" fill="#E70488"></circle>
</svg>
```
**HEX:** #E70488

**RGB:** ((31, 4, 136)


# Citing

When citing GeoPandas, you can use [Zenodo DOI](https://zenodo.org/record/3946761#.Xy24LC2ZPOQ) for each release. 
Below is the example of resulting BiBTeX record for GeoPandas 0.8.1 and a reference using APA 6<sup>th</sup> ed.

_Kelsey Jordahl, Joris Van den Bossche, Martin Fleischmann, Jacob Wasserman, James McBride, Jeffrey Gerard, … François Leblanc. (2020, July 15). geopandas/geopandas: v0.8.1 (Version v0.8.1). Zenodo. http://doi.org/10.5281/zenodo.3946761_


```
@software{kelsey_jordahl_2020_3946761,
  author       = {Kelsey Jordahl and
                  Joris Van den Bossche and
                  Martin Fleischmann and
                  Jacob Wasserman and
                  James McBride and
                  Jeffrey Gerard and
                  Jeff Tratner and
                  Matthew Perry and
                  Adrian Garcia Badaracco and
                  Carson Farmer and
                  Geir Arne Hjelle and
                  Alan D. Snow and
                  Micah Cochran and
                  Sean Gillies and
                  Lucas Culbertson and
                  Matt Bartos and
                  Nick Eubank and
                  maxalbert and
                  Aleksey Bilogur and
                  Sergio Rey and
                  Christopher Ren and
                  Dani Arribas-Bel and
                  Leah Wasser and
                  Levi John Wolf and
                  Martin Journois and
                  Joshua Wilson and
                  Adam Greenhall and
                  Chris Holdgraf and
                  Filipe and
                  François Leblanc},
  title        = {geopandas/geopandas: v0.8.1},
  month        = jul,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v0.8.1},
  doi          = {10.5281/zenodo.3946761},
  url          = {https://doi.org/10.5281/zenodo.3946761}
}
```

# Team

## Contributors

GeoPandas is developed by more than [100 volunteer contributors](https://github.com/geopandas/geopandas/graphs/contributors).

## Core developers
- Joris Van den Bossche - **lead maintainer** | [@jorisvandenbossche](https://github.com/jorisvandenbossche)
- Martin Fleischmann | [@martinfleis](https://github.com/martinfleis)
- James McBride | [@jdmcbr](https://github.com/jdmcbr)
- Brendan Ward | [@brendan-ward](https://github.com/brendan-ward)
- Levi Wolf | [@ljwolf](https://github.com/ljwolf)

## Founder

- Kelsey Jordahl | [@kjordahl](https://github.com/kjordahl)

## Alumni developers

- Jacob Wasserman | [@jwass](https://github.com/jwass)
# Roadmap

## Roadmap for GeoPandas 1.0

WIP
# Ecosystem

## GeoPandas dependencies

GeoPandas brings together the full capability of `pandas` and the open-source geospatial
tools `Shapely`, which brings manipulation and analysis of geometric objects backed by
[`GEOS`](https://trac.osgeo.org/geos) library, `Fiona`, allowing us to read and write
geographic data files using [`GDAL`](https://gdal.org), and `pyproj`, a library for
cartographic projections and coordinate transformations, which is a Python interface to
[`PROJ`](https://proj.org).

Furthermore, GeoPandas has several optional dependencies as `rtree`, `pygeos`,
`mapclassify`, or `geopy`.

### Required dependencies

#### [pandas](https://github.com/pandas-dev/pandas)
`pandas` is a Python package that provides fast, flexible, and expressive data
structures designed to make working with structured (tabular, multidimensional,
potentially heterogeneous) and time series data both easy and intuitive. It aims to be
the fundamental high-level building block for doing practical, real world data analysis
in Python. Additionally, it has the broader goal of becoming the most powerful and
flexible open source data analysis / manipulation tool available in any language. It is
already well on its way toward this goal.

#### [Shapely](https://github.com/Toblerity/Shapely)
`Shapely` is a BSD-licensed Python package for manipulation and analysis of planar
geometric objects. It is based on the widely deployed `GEOS` (the engine of PostGIS) and
`JTS` (from which `GEOS` is ported) libraries. `Shapely` is not concerned with data
formats or coordinate systems, but can be readily integrated with packages that are.

#### [Fiona](https://github.com/Toblerity/Fiona)
`Fiona` is `GDAL’s` neat and nimble vector API for Python programmers. Fiona is designed
to be simple and dependable. It focuses on reading and writing data in standard Python
IO style and relies upon familiar Python types and protocols such as files,
dictionaries, mappings, and iterators instead of classes specific to `OGR`. Fiona can
read and write real-world data using multi-layered GIS formats and zipped virtual file
systems and integrates readily with other Python GIS packages such as `pyproj`, `Rtree`,
and `Shapely`.

#### [pyproj](https://github.com/pyproj4/pyproj)
`pyproj` is a Python interface to `PROJ` (cartographic projections and coordinate
transformations library). GeoPandas uses a `pyproj.crs.CRS` object to keep track of the
projection of each `GeoSeries` and its `Transformer` object to manage re-projections.

### Optional dependencies

#### [rtree](https://github.com/Toblerity/rtree)
`Rtree` is a ctypes Python wrapper of `libspatialindex` that provides a number of
advanced spatial indexing features for the spatially curious Python user.

#### [PyGEOS](https://github.com/pygeos/pygeos)
`PyGEOS` is a C/Python library with vectorized geometry functions. The geometry
operations are done in the open-source geometry library `GEOS`. PyGEOS wraps these
operations in `NumPy` ufuncs providing a performance improvement when operating on
arrays of geometries.

#### [mapclassify](https://github.com/pysal/mapclassify)
`mapclassify` provides functionality for Choropleth map classification. Currently,
fifteen different classification schemes are available, including a highly-optimized
implementation of Fisher-Jenks optimal classification. Each scheme inherits a common
structure that ensures computations are scalable and supports applications in streaming
contexts.

#### [geopy](https://github.com/geopy/geopy)
`geopy` is a Python client for several popular geocoding web services. `geopy` makes it
easy for Python developers to locate the coordinates of addresses, cities, countries,
and landmarks across the globe using third-party geocoders and other data sources.

#### [matplotlib](https://github.com/matplotlib/matplotlib)
`Matplotlib` is a comprehensive library for creating static, animated, and interactive
visualizations in Python. Matplotlib produces publication-quality figures in a variety
of hardcopy formats and interactive environments across platforms. Matplotlib can be
used in Python scripts, the Python and IPython shell, web application servers, and
various graphical user interface toolkits.

## GeoPandas ecosystem

Various packages are built on top of GeoPandas addressing specific geospatial data
processing needs, analysis, and visualization. Below is an incomplete list (in no
particular order) of tools which form the GeoPandas-related Python ecosystem.

### Spatial analysis and Machine Learning

#### [PySAL](https://github.com/pysal/pysal)
`PySAL`, the Python spatial analysis library, is an open source cross-platform library
for geospatial data science with an emphasis on geospatial vector data written in
Python. `PySAL` is a family of packages, some of which are listed below.

##### [libpysal](https://github.com/pysal/libpysal)
`libpysal` provides foundational algorithms and data structures that support the rest of
the library. This currently includes the following modules: input/output (`io`), which
provides readers and writers for common geospatial file formats; weights (`weights`),
which provides the main class to store spatial weights matrices, as well as several
utilities to manipulate and operate on them; computational geometry (`cg`), with several
algorithms, such as Voronoi tessellations or alpha shapes that efficiently process
geometric shapes; and an additional module with example data sets (`examples`).

##### [esda](https://github.com/pysal/esda)
`esda` implements methods for the analysis of both global (map-wide) and local (focal)
spatial autocorrelation, for both continuous and binary data. In addition, the package
increasingly offers cutting-edge statistics about boundary strength and measures of
aggregation error in statistical analyses.

##### [segregation](https://github.com/pysal/segregation)
`segregation` package calculates over 40 different segregation indices and provides a
suite of additional features for measurement, visualization, and hypothesis testing that
together represent the state of the art in quantitative segregation analysis.

##### [mgwr](https://github.com/pysal/mgwr)
`mgwr` provides scalable algorithms for estimation, inference, and prediction using
single- and multi-scale geographically weighted regression models in a variety of
generalized linear model frameworks, as well as model diagnostics tools.

##### [tobler](https://github.com/pysal/tobler)
`tobler` provides functionality for areal interpolation and dasymetric mapping.
`tobler` includes functionality for interpolating data using area-weighted approaches,
regression model-based approaches that leverage remotely-sensed raster data as auxiliary
information, and hybrid approaches.

#### [movingpandas](https://github.com/anitagraser/movingpandas)
`MovingPandas` is a package for dealing with movement data. `MovingPandas` implements a
`Trajectory` class and corresponding methods based on GeoPandas. A trajectory has a
time-ordered series of point geometries. These points and associated attributes are
stored in a `GeoDataFrame`. `MovingPandas` implements spatial and temporal data access
and analysis functions as well as plotting functions.

#### [momepy](https://github.com/martinfleis/momepy)
`momepy` is a library for quantitative analysis of urban form - urban morphometrics. It
is built on top of `GeoPandas`, `PySAL` and `networkX`. `momepy` aims to provide a wide
range of tools for a systematic and exhaustive analysis of urban form. It can work with
a wide range of elements, while focused on building footprints and street networks.

#### [geosnap](https://github.com/spatialucr/geosnap)
`geosnap` makes it easier to explore, model, analyze, and visualize the social and
spatial dynamics of neighborhoods. `geosnap` provides a suite of tools for creating
socio-spatial datasets, harmonizing those datasets into consistent set of time-static
boundaries, modeling bespoke neighborhoods and prototypical neighborhood types, and
modeling neighborhood change using classic and spatial statistical methods. It also
provides a set of static and interactive visualization tools to help you display and
understand the critical information at each step of the process.

#### [mesa-geo](https://github.com/Corvince/mesa-geo)
`mesa-geo` implements a GeoSpace that can host GIS-based GeoAgents, which are like
normal Agents, except they have a shape attribute that is a `Shapely` object. You can
use `Shapely` directly to create arbitrary shapes, but in most cases you will want to
import your shapes from a file. Mesa-geo allows you to create GeoAgents from any vector
data file (e.g. shapefiles), valid GeoJSON objects or a GeoPandas `GeoDataFrame`.

#### [Pyspatialml](https://github.com/stevenpawley/Pyspatialml)
`Pyspatialml` is a Python module for applying `scikit-learn` machine learning models to
'stacks' of raster datasets. Pyspatialml includes functions and classes for working with
multiple raster datasets and performing a typical machine learning workflow consisting
of extracting training data and applying the predict or `predict_proba` methods of
`scikit-learn` estimators to a stack of raster datasets. Pyspatialml is built upon the
`rasterio` Python module for all of the heavy lifting, and is also designed for working
with vector data using the `geopandas` module.

#### [PyGMI](https://github.com/Patrick-Cole/pygmi)
`PyGMI` stands for Python Geoscience Modelling and Interpretation. It is a modelling and
interpretation suite aimed at magnetic, gravity and other datasets.

### Visualization

#### [hvPlot](https://hvplot.holoviz.org/user_guide/Geographic_Data.html#Geopandas)
`hvPlot` provides interactive Bokeh-based plotting for GeoPandas
dataframes and series using the same API as the Matplotlib `.plot()`
support that comes with GeoPandas. hvPlot makes it simple to pan and zoom into
your plots, use widgets to explore multidimensional data, and render even the
largest datasets in web browsers using [Datashader](https://datashader.org).

#### [contextily](https://github.com/geopandas/contextily)
`contextily` is a small Python 3 (3.6 and above) package to retrieve tile maps from the
internet. It can add those tiles as basemap to `matplotlib` figures or write tile maps
to disk into geospatial raster files. Bounding boxes can be passed in both WGS84
(EPSG:4326) and Spheric Mercator (EPSG:3857).

#### [cartopy](https://github.com/SciTools/cartopy)
`Cartopy` is a Python package designed to make drawing maps for data analysis and
visualisation easy. It features: object oriented projection definitions; point, line,
polygon and image transformations between projections; integration to expose advanced
mapping in `Matplotlib` with a simple and intuitive interface; powerful vector data
handling by integrating shapefile reading with `Shapely` capabilities.

#### [bokeh](https://github.com/bokeh/bokeh)
`Bokeh` is an interactive visualization library for modern web browsers. It provides
elegant, concise construction of versatile graphics, and affords high-performance
interactivity over large or streaming datasets. `Bokeh` can help anyone who would like
to quickly and easily make interactive plots, dashboards, and data applications.

#### [folium](https://github.com/python-visualization/folium)
`folium` builds on the data wrangling strengths of the Python ecosystem and the mapping
strengths of the `Leaflet.js` library. Manipulate your data in Python, then visualize it
in a `Leaflet` map via `folium`.

#### [kepler.gl](https://github.com/keplergl/kepler.gl)
`Kepler.gl` is a data-agnostic, high-performance web-based application for visual
exploration of large-scale geolocation data sets. Built on top of Mapbox GL and
`deck.gl`, `kepler.gl` can render millions of points representing thousands of trips and
perform spatial aggregations on the fly.

#### [geoplot](https://github.com/ResidentMario/geoplot)
`geoplot` is a high-level Python geospatial plotting library. It's an extension to
`cartopy` and `matplotlib` which makes mapping easy: like `seaborn` for geospatial. It
comes with the high-level plotting API, native projection support and compatibility with
`matplotlib`.

#### [GeoViews](https://github.com/holoviz/geoviews)
`GeoViews` is a Python library that makes it easy to explore and
visualize any data that includes geographic locations, with native
support for GeoPandas dataframes and series objects.  It has
particularly powerful support for multidimensional meteorological and
oceanographic datasets, such as those used in weather, climate, and
remote sensing research, but is useful for almost anything that you
would want to plot on a map!

#### [EarthPy](https://github.com/earthlab/earthpy)
`EarthPy` is a python package that makes it easier to plot and work with spatial raster
and vector data using open source tools. `Earthpy` depends upon `geopandas` which has a
focus on vector data and `rasterio` with facilitates input and output of raster data
files. It also requires `matplotlib` for plotting operations. `EarthPy’s` goal is to
make working with spatial data easier for scientists.

#### [splot](https://github.com/pysal/splot)
`splot` provides statistical visualizations for spatial analysis. It methods for
visualizing global and local spatial autocorrelation (through Moran scatterplots and
cluster maps), temporal analysis of cluster dynamics (through heatmaps and rose
diagrams), and multivariate choropleth mapping (through value-by-alpha maps). A high
level API supports the creation of publication-ready visualizations

#### [legendgram](https://github.com/pysal/legendgram)
`legendgram` is a small package that provides "legendgrams" legends that visualize the
distribution of observations by color in a given map. These distributional
visualizations for map classification schemes assist in analytical cartography and
spatial data visualization.

### Geometry manipulation

#### [TopoJSON](https://github.com/mattijn/topojson)
`topojson` is a library for creating a TopoJSON encoding of nearly any 
geographical object in Python. With topojson it is possible to reduce the size of
your geographical data, typically by orders of magnitude. It is able to do so through
eliminating redundancy through computation of a topology, fixed-precision integer
encoding of coordinates, and simplification and quantization of arcs.

#### [geocube](https://github.com/corteva/geocube)
Tool to convert geopandas vector data into rasterized `xarray` data.

### Data retrieval

#### [OSMnx](https://github.com/gboeing/osmnx)
`OSMnx` is a Python package that lets you download spatial data from OpenStreetMap and
model, project, visualize, and analyze real-world street networks. You can download and
model walkable, drivable, or bikeable urban networks with a single line of Python code
and then easily analyze and visualize them. You can just as easily download and work with
other infrastructure types, amenities/points of interest, building footprints, elevation
data, street bearings/orientations, and speed/travel time.

#### [pyrosm](https://github.com/HTenkanen/pyrosm)
`Pyrosm` is a Python library for reading OpenStreetMap data from Protocolbuffer Binary
Format -files (`*.osm.pbf`) into Geopandas `GeoDataFrames`. Pyrosm makes it easy to
extract various datasets from OpenStreetMap pbf-dumps including e.g. road networks,
buildings, Points of Interest (POI), landuse and natural elements. Also fully customized
queries are supported which makes it possible to parse the data from OSM with more
specific filters.

#### [geobr](https://github.com/ipeaGIT/geobr)
`geobr` is a computational package to download official spatial data sets of Brazil. The
package includes a wide range of geospatial data in geopackage format (like shapefiles
but better), available at various geographic scales and for various years with
harmonized attributes, projection and topology.

#### [cenpy](https://github.com/cenpy-devs/cenpy)
An interface to explore and query the US Census API and return Pandas `Dataframes`. This
package is intended for exploratory data analysis and draws inspiration from
sqlalchemy-like interfaces and `acs.R`. With separate APIs for application developers
and folks who only want to get their data quickly & painlessly, `cenpy` should meet the
needs of most who aim to get US Census Data into Python.

```{admonition} Expand this page
Do know a package which should be here? [Let us
know](https://github.com/geopandas/geopandas/issues) or [add it by
yourself](contributing.rst)!
```
Community
---------

GeoPandas is a community-led project written, used and supported by a wide range of
people from all around of world of a large variety of backgrounds. Everyone is welcome,
each small contribution, no matter if it is a fix of a typo in the documentation, bug
report, an idea, or a question, is valuable. As a member of our community, you should
adhere to the principles presented in the :doc:`Code of Conduct
<community/code_of_conduct>`.

If you'd like to contribute, please read the :doc:`Contributing guide
<community/contributing>`. It will help you to understand the way GeoPandas development
works and us to review your contribution.

GeoPandas is a part of the broader Python :doc:`ecosystem <community/ecosystem>`. It
depends on a range of great tools, and various packages are built on top of GeoPandas
addressing specific needs in geospatial data processing, analysis and visualization.

.. container:: button

    :doc:`Contributing <community/contributing>` :doc:`Code of Conduct <community/code_of_conduct>`
    :doc:`Ecosystem <community/ecosystem>`

.. toctree::
  :maxdepth: 2
  :caption: Community
  :hidden:

  Contributing <community/contributing>
  Code of Conduct <community/code_of_conduct>
  Ecosystem <community/ecosystem>
Documentation
-------------

The documentation of GeoPandas consists of four parts - :doc:`User Guide <docs/user_guide>` with explanation of the basic functionality, :doc:`Advanced Guide <docs/advanced_guide>` covering topics which assume knowledge of basics, :doc:`Examples <gallery/index>`, and :doc:`API reference <docs/reference>` detailing every class, method, function and attribute used implemented by GeoPandas.

.. container:: button

    :doc:`User Guide <docs/user_guide>` :doc:`Advanced Guide <docs/advanced_guide>`
    :doc:`Examples <gallery/index>` :doc:`API reference <docs/reference>`


.. toctree::
  :maxdepth: 2
  :caption: Documentation

  User Guide <docs/user_guide>
  Advanced Guide <docs/advanced_guide>
  API reference <docs/reference>
  Changelog <docs/changelog>GeoPandas |version|
===================

GeoPandas is an open source project to make working with geospatial
data in python easier.  GeoPandas extends the datatypes used by
`pandas`_ to allow spatial operations on geometric types.  Geometric
operations are performed by `shapely`_.  Geopandas further depends on
`fiona`_ for file access and `matplotlib`_ for plotting.

.. _pandas: http://pandas.pydata.org
.. _shapely: https://shapely.readthedocs.io
.. _fiona: https://fiona.readthedocs.io
.. _matplotlib: http://matplotlib.org

Description
-----------

The goal of GeoPandas is to make working with geospatial data in
python easier.  It combines the capabilities of pandas and shapely,
providing geospatial operations in pandas and a high-level interface
to multiple geometries to shapely.  GeoPandas enables you to easily do
operations in python that would otherwise require a spatial database
such as PostGIS.

.. toctree::
   :hidden:

   Home <self>
   About <about>
   Getting started <getting_started>
   Documentation <docs>
   Community <community>

.. container:: button

    :doc:`Getting started <getting_started>` :doc:`Documentation <docs>`
    :doc:`About GeoPandas <about>` :doc:`Community <community>`

Useful links
------------

`Binary Installers (PyPI) <https://pypi.org/project/geopandas/>`_ | `Source Repository (GitHub) <https://github.com/geopandas/geopandas>`_ | `Issues & Ideas <https://github.com/geopandas/geopandas/issues>`_ | `Q&A Support <https://stackoverflow.com/questions/tagged/geopandas>`_


Supported by
------------

.. image:: https://numfocus.org/wp-content/uploads/2017/07/NumFocus_LRG.png
    :alt: numfocus
    :width: 400
    :target: https://numfocus.org

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
{{ fullname | escape | underline}}


.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}

{% if objtype in ['class', 'method', 'function'] %}

.. include:: {{module}}.{{objname}}.examples

.. raw:: html

    <div class="clear"></div>

{% endif %}
{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. automethod:: {{ (module.split('.')[1:] + [objname]) | join('.') }}
Examples Gallery
================

The following examples show off the functionality in GeoPandas. They highlight many of the things you can do with this package, and show off some best-practices.


.. nbgallery::
    :name: nbshpinx-gallery
    :glob:

    ./*Contributing to GeoPandas
=========================

(Contribution guidelines largely copied from `pandas <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_)

Overview
--------

Contributions to GeoPandas are very welcome.  They are likely to
be accepted more quickly if they follow these guidelines.

At this stage of GeoPandas development, the priorities are to define a
simple, usable, and stable API and to have clean, maintainable,
readable code.  Performance matters, but not at the expense of those
goals.

In general, GeoPandas follows the conventions of the pandas project
where applicable.

In particular, when submitting a pull request:

- All existing tests should pass.  Please make sure that the test
  suite passes, both locally and on
  `GitHub Actions <https://github.com/geopandas/geopandas/actions>`_.  Status on
  GHA will be visible on a pull request. GHA are automatically enabled
  on your own fork as well. To trigger a check, make a PR to your own fork.

- New functionality should include tests.  Please write reasonable
  tests for your code and make sure that they pass on your pull request.

- Classes, methods, functions, etc. should have docstrings.  The first
  line of a docstring should be a standalone summary.  Parameters and
  return values should be documented explicitly.

- Follow PEP 8 when possible. We use `Black
  <https://black.readthedocs.io/en/stable/>`_ and `Flake8
  <http://flake8.pycqa.org/en/latest/>`_ to ensure a consistent code
  format throughout the project. For more details see
  :ref:`below <contributing_style>`.

- Imports should be grouped with standard library imports first,
  3rd-party libraries next, and GeoPandas imports third.  Within each
  grouping, imports should be alphabetized.  Always use absolute
  imports when possible, and explicit relative imports for local
  imports when necessary in tests.

- GeoPandas supports Python 3.7+ only. The last version of GeoPandas
  supporting Python 2 is 0.6.


Seven Steps for Contributing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are seven basic steps to contributing to *GeoPandas*:

1) Fork the *GeoPandas* git repository
2) Create a development environment
3) Install *GeoPandas* dependencies
4) Make a ``development`` build of *GeoPandas*
5) Make changes to code and add tests
6) Update the documentation
7) Submit a Pull Request

Each of these 7 steps is detailed below.


1) Forking the *GeoPandas* repository using Git
------------------------------------------------

To the new user, working with Git is one of the more daunting aspects of contributing to *GeoPandas**.
It can very quickly become overwhelming, but sticking to the guidelines below will help keep the process
straightforward and mostly trouble free.  As always, if you are having difficulties please
feel free to ask for help.

The code is hosted on `GitHub <https://github.com/geopandas/geopandas>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* Software Carpentry's `Git Tutorial <http://swcarpentry.github.io/git-novice/>`_
* `Atlassian <https://www.atlassian.com/git/tutorials/what-is-version-control>`_
* the `GitHub help pages <http://help.github.com/>`_.
* Matthew Brett's `Pydagogue <https://matthew-brett.github.io/pydagogue/>`_.

Getting started with Git
~~~~~~~~~~~~~~~~~~~~~~~~

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
~~~~~~~

You will need your own fork to work on the code. Go to the `GeoPandas project
page <https://github.com/geopandas/geopandas>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone git@github.com:your-user-name/geopandas.git geopandas-yourname
    cd geopandas-yourname
    git remote add upstream git://github.com/geopandas/geopandas.git

This creates the directory `geopandas-yourname` and connects your repository to
the upstream (main project) *GeoPandas* repository.

The testing suite will run automatically on GitHub Actions once your pull request is
submitted. The test suite will also automatically run on your branch so you can
check it prior to submitting the pull request.

Creating a branch
~~~~~~~~~~~~~~~~~~

You want your main branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *GeoPandas*. You can have many shiny-new-features
and switch in between them using the git checkout command.

To update this branch, you need to retrieve the changes from the main branch::

    git fetch upstream
    git rebase upstream/main

This will replay your commits on top of the latest GeoPandas git main.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes, you will need to ``stash`` them prior
to updating.  This will effectively store your changes and they can be reapplied
after updating.

.. _contributing.dev_env:

2) Creating a development environment
---------------------------------------
A development environment is a virtual space where you can keep an independent installation of *GeoPandas*.
This makes it easy to keep both a stable version of python in one place you use for work, and a development
version (which you may break while playing with code) in another.

An easy way to create a *GeoPandas* development environment is as follows:

- Install either `Anaconda <http://docs.continuum.io/anaconda/>`_ or
  `miniconda <http://conda.pydata.org/miniconda.html>`_
- Make sure that you have :ref:`cloned the repository <contributing.forking>`
- ``cd`` to the *geopandas** source directory

Using the provided environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*GeoPandas* provides an environment which includes the required dependencies for development.
The environment file is located in the top level of the repo and is named ``environment-dev.yml``.
You can create this environment by navigating to the the *GeoPandas* source directory
and running::

      conda env create -f environment-dev.yml

This will create a new conda environment named ``geopandas_dev``.

Creating the environment manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, it is possible to create a development environment manually.  To do this,
tell conda to create a new environment named ``geopandas_dev``, or any other name you would like
for this environment, by running::

      conda create -n geopandas_dev python

This will create the new environment, and not touch any of your existing environments,
nor any existing python installation.

Working with the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To work in this environment, you need to ``activate`` it. The instructions below
should work for both Windows, Mac and Linux::

      conda activate geopandas_dev

Once your environment is activated, you will see a confirmation message to
indicate you are in the new development environment.

To view your environments::

      conda info -e

To return to you home root environment::

      conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

At this point you can easily do a *development* install, as detailed in the next sections.


3) Installing Dependencies
--------------------------

To run *GeoPandas* in an development environment, you must first install
*GeoPandas*'s dependencies. If you used the provided environment in section 2, skip this
step and continue to section 4. If you created the environment manually, we suggest installing
dependencies using the following commands (executed after your development environment has been activated)::

    conda install -c conda-forge pandas fiona shapely pyproj rtree pytest

This should install all necessary dependencies.


4) Making a development build
-----------------------------

Once dependencies are in place, make an in-place build by navigating to the git
clone of the *GeoPandas* repository and running::

    python setup.py develop


5) Making changes and writing tests
-------------------------------------

*GeoPandas* is serious about testing and strongly encourages contributors to embrace
`test-driven development (TDD) <http://en.wikipedia.org/wiki/Test-driven_development>`_.
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

Adding tests is one of the most common requests after code is pushed to *GeoPandas*.  Therefore,
it is worth getting in the habit of writing tests ahead of time so this is never an issue.

*GeoPandas* uses the `pytest testing system
<http://doc.pytest.org/en/latest/>`_ and the convenient
extensions in `numpy.testing
<http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.

Writing tests
~~~~~~~~~~~~~

All tests should go into the ``tests`` directory. This folder contains many
current examples of tests, and we suggest looking to these for inspiration.

The ``.util`` module has some special ``assert`` functions that
make it easier to make statements about whether GeoSeries or GeoDataFrame
objects are equivalent. The easiest way to verify that your code is correct is to
explicitly construct the result you expect, then compare the actual result to
the expected correct result, using eg the function ``assert_geoseries_equal``.

Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

The tests can then be run directly inside your Git clone (without having to
install *GeoPandas*) by typing::

    pytest

6) Updating the Documentation
-----------------------------

*GeoPandas* documentation resides in the ``doc`` folder. Changes to the docs are made by
modifying the appropriate file in the ``source`` folder within ``doc``. *GeoPandas* docs use
mixture of reStructuredText syntax for ``rst`` files, `which is explained here
<http://www.sphinx-doc.org/en/stable/rest.html#rst-primer>`_ and MyST syntax for ``md``
files `explained here <https://myst-parser.readthedocs.io/en/latest/index.html>`_.
The docstrings follow the `Numpy Docstring standard
<https://github.com/numpy/numpy/blob/main/doc/HOWTO_DOCUMENT.rst.txt>`_. Some pages
and examples are Jupyter notebooks converted to docs using `nbsphinx
<https://nbsphinx.readthedocs.io/>`_. Jupyter notebooks should be stored without the output.

We highly encourage you to follow the `Google developer documentation style guide
<https://developers.google.com/style/highlights>`_ when updating or creating new documentation.

Once you have made your changes, you may try if they render correctly by
building the docs using sphinx. To do so, you can navigate to the `doc` folder::

    cd doc

and type::

    make html

The resulting html pages will be located in ``doc/build/html``.

In case of any errors, you can try to use ``make html`` within a new environment based on
environment.yml specification in the ``doc`` folder. You may need to register Jupyter kernel as
``geopandas_docs``. Using conda::

    cd doc
    conda env create -f environment.yml
    conda activate geopandas_docs
    python -m ipykernel install --user --name geopandas_docs
    make html

For minor updates, you can skip the ``make html`` part as reStructuredText and MyST
syntax are usually quite straightforward.


7) Submitting a Pull Request
------------------------------

Once you've made changes and pushed them to your forked repository, you then
submit a pull request to have them integrated into the *GeoPandas* code base.

You can find a pull request (or PR) tutorial in the `GitHub's Help Docs <https://help.github.com/articles/using-pull-requests/>`_.

.. _contributing_style:

Style Guide & Linting
---------------------

GeoPandas follows the `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ standard
and uses `Black <https://black.readthedocs.io/en/stable/>`_ and
`Flake8 <http://flake8.pycqa.org/en/latest/>`_ to ensure a consistent code
format throughout the project.

Continuous Integration (GitHub Actions) will run those tools and
report any stylistic errors in your code. Therefore, it is helpful before
submitting code to run the check yourself::

   black geopandas
   git diff upstream/main -u -- "*.py" | flake8 --diff

to auto-format your code. Additionally, many editors have plugins that will
apply ``black`` as you edit files.

Optionally (but recommended), you can setup `pre-commit hooks <https://pre-commit.com/>`_
to automatically run ``black`` and ``flake8`` when you make a git commit. If you did not
use the provided development environment in ``environment-dev.yml``, you must first install ``pre-commit``::

   $ python -m pip install pre-commit

From the root of the geopandas repository, you should then install the
``pre-commit`` included in *GeoPandas*::

   $ pre-commit install

Then ``black`` and ``flake8`` will be run automatically
each time you commit changes. You can skip these checks with
``git commit --no-verify``.

Commit message conventions
--------------------------

Commit your changes to your local repository with an explanatory message. GeoPandas
uses the pandas convention for commit message prefixes and layout. Here are
some common prefixes along with general guidelines for when to use them:

* ENH: Enhancement, new functionality
* BUG: Bug fix
* DOC: Additions/updates to documentation
* TST: Additions/updates to tests
* BLD: Updates to the build process/scripts
* PERF: Performance improvement
* TYP: Type annotations
* CLN: Code cleanup

The following defines how a commit message should be structured. Please refer to the
relevant GitHub issues in your commit message using GH1234 or #1234. Either style
is fine, but the former is generally preferred:

* a subject line with `< 80` chars.
* One blank line.
* Optionally, a commit message body.

Now you can commit your changes in your local repository::

    git commit -m
GeoPandas Project Code of Conduct
=================================

Behind the GeoPandas Project is an engaged and respectful community made up of
people from all over the world and with a wide range of backgrounds.
Naturally, this implies diversity of ideas and perspectives on often
complex problems. Disagreement and healthy discussion of conflicting
viewpoints is welcome: the best solutions to hard problems rarely come from a single
angle. But disagreement is not an excuse for aggression: humans tend to take
disagreement personally and easily drift into behavior that ultimately
degrades a community. This is particularly acute with online communication
across language and cultural gaps, where many cues of human behavior are
unavailable. We are outlining here a set of principles and processes to support a
healthy community in the face of these challenges.

Fundamentally, we are committed to fostering a productive, harassment-free
environment for everyone. Rather than considering this code an exhaustive list
of things that you can’t do, take it in the spirit it is intended - a guide to
make it easier to enrich all of us and the communities in which we participate.

Importantly: as a member of our community, *you are also a steward of these
values*. Not all problems need to be resolved via formal processes, and often
a quick, friendly but clear word on an online forum or in person can help
resolve a misunderstanding and de-escalate things.

However, sometimes these informal processes may be inadequate: they fail to
work, there is urgency or risk to someone, nobody is intervening publicly and
you don't feel comfortable speaking in public, etc. For these or other
reasons, structured follow-up may be necessary and here we provide the means
for that: we welcome reports by emailing
`geopandas-conduct@googlegroups.com <mailto:geopandas-conduct@googlegroups.com>`__ or by filling out `this
form <https://docs.google.com/forms/d/e/1FAIpQLSd8Tbi2zNl1i2N9COX0yavHEqTGFIPQ1_cLcy1A3JgVc1OrAQ/viewform>`__.

This code applies equally to founders, developers, mentors and new community
members, in all spaces managed by the GeoPandas Project. This
includes the mailing lists, our GitHub organization, our chat room, in-person
events, and any other forums created by the project team. In addition,
violations of this code outside these spaces may affect a person's ability to
participate within them.

By embracing the following principles, guidelines and actions to follow or
avoid, you will help us make Jupyter a welcoming and productive community. Feel
free to contact the Code of Conduct Committee at
`geopandas-conduct@googlegroups.com <mailto:geopandas-conduct@googlegroups.com>`__ with any questions.

1. **Be friendly and patient**.

2. **Be welcoming**. We strive to be a community that welcomes and supports
   people of all backgrounds and identities. This includes, but is not limited
   to, members of any race, ethnicity, culture, national origin, color,
   immigration status, social and economic class, educational level, sex, sexual
   orientation, gender identity and expression, age, physical appearance, family
   status, technological or professional choices, academic
   discipline, religion, mental ability, and physical ability.

3. **Be considerate**. Your work will be used by other people, and you in turn
   will depend on the work of others. Any decision you take will affect users
   and colleagues, and you should take those consequences into account when
   making decisions. Remember that we're a world-wide community. You may be
   communicating with someone with a different primary language or cultural
   background.

4. **Be respectful**. Not all of us will agree all the time, but disagreement is
   no excuse for poor behavior or poor manners. We might all experience some
   frustration now and then, but we cannot allow that frustration to turn into a
   personal attack. It’s important to remember that a community where people
   feel uncomfortable or threatened is not a productive one.

5. **Be careful in the words that you choose**. Be kind to others. Do not insult
   or put down other community members. Harassment and other exclusionary
   behavior are not acceptable. This includes, but is not limited to:

   -  Violent threats or violent language directed against another person
   -  Discriminatory jokes and language
   -  Posting sexually explicit or violent material
   -  Posting (or threatening to post) other people's personally identifying information ("doxing")
   -  Personal insults, especially those using racist, sexist, and xenophobic terms
   -  Unwelcome sexual attention
   -  Advocating for, or encouraging, any of the above behavior
   -  Repeated harassment of others. In general, if someone asks you to stop, then stop

6. **Moderate your expectations**. Please respect that community members choose
   how they spend their time in the project. A thoughtful question about your
   expectations is preferable to demands for another person's time.

7. **When we disagree, try to understand why**. Disagreements, both social and
   technical, happen all the time and the GeoPandas Project is no exception. Try to
   understand where others are coming from, as seeing a question from their
   viewpoint may help find a new path forward. And don’t forget that it is
   human to err: blaming each other doesn’t get us anywhere, while we can learn
   from mistakes to find better solutions.

8. **A simple apology can go a long way**. It can often de-escalate a situation,
   and telling someone that you are sorry is an act of empathy that doesn’t
   automatically imply an admission of guilt.

Reporting
---------

If you believe someone is violating the code of conduct, please report this in
a timely manner. Code of conduct violations reduce the value of the community
for everyone and we take them seriously.

You can file a report by emailing
`geopandas-conduct@googlegroups.com <mailto:geopandas-conduct@googlegroups.com>`__ or by filing out
`this form <https://docs.google.com/forms/d/e/1FAIpQLSd8Tbi2zNl1i2N9COX0yavHEqTGFIPQ1_cLcy1A3JgVc1OrAQ/viewform>`__.

The online form gives you the option to keep your report anonymous or request
that we follow up with you directly. While we cannot follow up on an anonymous
report, we will take appropriate action.

Messages sent to the e-mail address or through the form will be sent
only to the Code of Conduct Committee, which currently consists of:

- Hannah Aizenman
- Joris Van den Bossche
- Martin Fleischmann

Enforcement
-----------

Enforcement procedures within the GeoPandas Project follow Project Jupyter's `Enforcement
Manual <https://github.com/jupyter/governance/blob/master/conduct/enforcement.md>`__.
For information on enforcement, please view the `original
manual <https://github.com/jupyter/governance/blob/master/conduct/enforcement.md>`__.

Original text courtesy of the `Speak
Up! <http://web.archive.org/web/20141109123859/http://speakup.io/coc.html>`__,
`Django <https://www.djangoproject.com/conduct>`__ and
`Jupyter <https://github.com/jupyter/governance/blob/master/conduct/code_of_conduct.md>`__
Projects, modified by the GeoPandas Project. We are grateful to those projects for
contributing these materials under open licensing terms for us to easily reuse.

All content on this page is licensed under a `Creative Commons
Attribution <http://creativecommons.org/licenses/by/3.0/>`__ license.
Installation
============

GeoPandas depends for its spatial functionality on a large geospatial, open
source stack of libraries (`GEOS`_, `GDAL`_, `PROJ`_). See the
:ref:`dependencies` section below for more details. Those base C
libraries can sometimes be a challenge to install. Therefore, we advise you
to closely follow the recommendations below to avoid installation problems.

.. _install-conda:

Installing with Anaconda / conda
--------------------------------

To install GeoPandas and all its dependencies, we recommend to use the `conda`_
package manager. This can be obtained by installing the
`Anaconda Distribution`_ (a free Python distribution for data science), or
through `miniconda`_ (minimal distribution only containing Python and the
`conda`_ package manager). See also the `installation docs
<https://conda.io/docs/user-guide/install/download.html>`__ for more information
on how to install Anaconda or miniconda locally.

The advantage of using the `conda`_ package manager is that it provides
pre-built binaries for all the required and optional dependencies of GeoPandas
for all platforms (Windows, Mac, Linux).

To install the latest version of GeoPandas, you can then do::

    conda install geopandas


Using the conda-forge channel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`conda-forge`_ is a community effort that provides conda packages for a wide
range of software. It provides the *conda-forge* package channel for conda from
which packages can be installed, in addition to the "*defaults*" channel
provided by Anaconda.
Depending on what other packages you are working with, the *defaults* channel
or *conda-forge* channel may be better for your needs (e.g. some packages are
available on *conda-forge* and not on *defaults*).

GeoPandas and all its dependencies are available on the *conda-forge*
channel, and can be installed as::

    conda install --channel conda-forge geopandas

.. note::

    We strongly recommend to either install everything from the *defaults*
    channel, or everything from the *conda-forge* channel. Ending up with a
    mixture of packages from both channels for the dependencies of GeoPandas
    can lead to import problems.
    See the `conda-forge section on using multiple channels
    <http://conda-forge.org/docs/user/tipsandtricks.html#using-multiple-channels>`__
    for more details.


Creating a new environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating a new environment is not strictly necessary, but given that installing
other geospatial packages from different channels may cause dependency conflicts
(as mentioned in the note above), it can be good practice to install the geospatial
stack in a clean environment starting fresh.

The following commands create a new environment with the name ``geo_env``,
configures it to install packages always from conda-forge, and installs
GeoPandas in it::

    conda create -n geo_env
    conda activate geo_env
    conda config --env --add channels conda-forge
    conda config --env --set channel_priority strict
    conda install python=3 geopandas


.. _install-pip:

Installing with pip
-------------------

GeoPandas can also be installed with pip, if all dependencies can be installed
as well::

    pip install geopandas

.. _install-deps:

.. warning::

    When using pip to install GeoPandas, you need to make sure that all dependencies are
    installed correctly.

    - `fiona`_ provides binary wheels with the dependencies included for Mac and Linux,
      but not for Windows.
    - `pyproj`_, `rtree`_, and `shapely`_ provide binary wheels with dependencies included
      for Mac, Linux, and Windows.
    - Windows wheels for `shapely`, `fiona`, `pyproj` and `rtree`
      can be found at `Christopher Gohlke's website
      <https://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.

    Depending on your platform, you might need to compile and install their
    C dependencies manually. We refer to the individual packages for more
    details on installing those.
    Using conda (see above) avoids the need to compile the dependencies yourself.

Installing from source
----------------------

You may install the latest development version by cloning the
`GitHub` repository and using pip to install from the local directory::

    git clone https://github.com/geopandas/geopandas.git
    cd geopandas
    pip install .

It is also possible to install the latest development version
directly from the GitHub repository with::

    pip install git+git://github.com/geopandas/geopandas.git

For installing GeoPandas from source, the same :ref:`note <install-deps>` on
the need to have all dependencies correctly installed applies. But, those
dependencies can also be installed independently with conda before installing
GeoPandas from source::

    conda install pandas fiona shapely pyproj rtree

See the :ref:`section on conda <install-conda>` above for more details on
getting running with Anaconda.

.. _dependencies:

Dependencies
------------

Required dependencies:

- `numpy`_
- `pandas`_ (version 0.25 or later)
- `shapely`_ (interface to `GEOS`_)
- `fiona`_ (interface to `GDAL`_)
- `pyproj`_ (interface to `PROJ`_; version 2.2.0 or later)
- `packaging`_

Further, optional dependencies are:

- `rtree`_ (optional; spatial index to improve performance and required for
  overlay operations; interface to `libspatialindex`_)
- `psycopg2`_ (optional; for PostGIS connection)
- `GeoAlchemy2`_ (optional; for writing to PostGIS)
- `geopy`_ (optional; for geocoding)


For plotting, these additional packages may be used:

- `matplotlib`_ (>= 3.1.0)
- `mapclassify`_ (>= 2.4.0)


Using the optional PyGEOS dependency
------------------------------------

Work is ongoing to improve the performance of GeoPandas. Currently, the
fast implementations of basic spatial operations live in the `PyGEOS`_
package (but work is under way to contribute those improvements to Shapely).
Starting with GeoPandas 0.8, it is possible to optionally use those
experimental speedups by installing PyGEOS. This can be done with conda
(using the conda-forge channel) or pip::

    # conda
    conda install pygeos --channel conda-forge
    # pip
    pip install pygeos

More specifically, whether the speedups are used or not is determined by:

- If PyGEOS >= 0.8 is installed, it will be used by default (but installing
  GeoPandas will not yet automatically install PyGEOS as dependency, you need
  to do this manually).

- You can still toggle the use of PyGEOS when it is available, by:

  - Setting an environment variable (``USE_PYGEOS=0/1``). Note this variable
    is only checked at first import of GeoPandas.
  - Setting an option: ``geopandas.options.use_pygeos = True/False``. Note,
    although this variable can be set during an interactive session, it will
    only work if the GeoDataFrames you use are created (e.g. reading a file
    with ``read_file``) after changing this value.

.. warning::

    The use of PyGEOS is experimental! Although it is passing all tests,
    there might still be issues and not all functions of GeoPandas will
    already benefit from speedups (one known issue: the `to_crs` coordinate
    transformations lose the z coordinate). But trying this out is very welcome!
    Any issues you encounter (but also reports of successful usage are
    interesting!) can be reported at https://gitter.im/geopandas/geopandas
    or https://github.com/geopandas/geopandas/issues


.. _PyPI: https://pypi.python.org/pypi/geopandas

.. _GitHub: https://github.com/geopandas/geopandas

.. _numpy: http://www.numpy.org

.. _pandas: http://pandas.pydata.org

.. _shapely: https://shapely.readthedocs.io

.. _fiona: https://fiona.readthedocs.io

.. _matplotlib: http://matplotlib.org

.. _geopy: https://github.com/geopy/geopy

.. _psycopg2: https://pypi.python.org/pypi/psycopg2

.. _GeoAlchemy2: https://geoalchemy-2.readthedocs.io/

.. _mapclassify: http://pysal.org/mapclassify

.. _pyproj: https://github.com/pyproj4/pyproj

.. _rtree: https://github.com/Toblerity/rtree

.. _libspatialindex: https://github.com/libspatialindex/libspatialindex

.. _conda: https://conda.io/en/latest/

.. _Anaconda distribution: https://www.anaconda.com/distribution/

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html

.. _conda-forge: https://conda-forge.org/

.. _GDAL: https://www.gdal.org/

.. _GEOS: https://geos.osgeo.org

.. _PROJ: https://proj.org/

.. _PyGEOS: https://github.com/pygeos/pygeos/

.. _packaging: https://packaging.pypa.io/en/latest/.. _reference:

API Reference
=============

The API Reference provides an overview of all public objects, functions and methods implemented in GeoPandas. All classes and function exposed in ``geopandas.*`` namespace plus those listed in the reference are public.

.. warning::
   The ``geopandas.array`` and ``geopandas.base`` modules are private. Stable functionality in such modules is not guaranteed.

.. toctree::
  :maxdepth: 2

  GeoSeries <reference/geoseries>
  GeoDataFrame <reference/geodataframe>
  Input/output <reference/io>
  Tools <reference/tools>
  Spatial index <reference/sindex>
  Testing <reference/testing>

User Guide
==========

The User Guide covers different parts of basic usage of GeoPandas. Each page focuses on a single topic and outlines how it is implemented in GeoPandas, with reproducible examples.

If you don't know anything about GeoPandas, start with the :doc:`Introduction to GeoPandas <../getting_started/introduction>`.

Advanced topics can be found in the :doc:`Advanced Guide <advanced_guide>` and further specification in the :doc:`API Reference <reference>`.

.. toctree::
  :maxdepth: 2

  Data Structures <user_guide/data_structures>
  Reading and Writing Files <user_guide/io>
  Indexing and Selecting Data <user_guide/indexing>
  Making Maps and plots <user_guide/mapping>
  Interactive mapping <user_guide/interactive_mapping>
  Managing Projections <user_guide/projections>
  Geometric Manipulations <user_guide/geometric_manipulations>
  Set Operations with overlay <user_guide/set_operations>
  Aggregation with dissolve <user_guide/aggregation_with_dissolve>
  Merging Data <user_guide/mergingdata>
  Geocoding <user_guide/geocoding>
Advanced Guide
==============

The Advanced Guide covers advanced usage of GeoPandas. Each page focuses on a single
topic and outlines how it is implemented in GeoPandas, with reproducible examples.

If you don't know anything about GeoPandas, start with the :doc:`Introduction to
GeoPandas <../getting_started/introduction>`.

Basic topics can be found in the :doc:`User Guide <user_guide>` and further
specification in the :doc:`API Reference <reference>`.

.. note::
   This section is currently work in progress. See the available pages below.

.. toctree::
  :maxdepth: 2

  user_guide/missing_empty
  user_guide/reproject_fiona
.. include:: ../../../CHANGELOG.md
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas
   import matplotlib.pyplot as plt
   plt.close('all')


Set-Operations with Overlay
============================

When working with multiple spatial datasets -- especially multiple *polygon* or
*line* datasets -- users often wish to create new shapes based on places where
those datasets overlap (or don't overlap). These manipulations are often
referred using the language of sets -- intersections, unions, and differences.
These types of operations are made available in the *geopandas* library through
the :meth:`~geopandas.GeoDataFrame.overlay` method.

The basic idea is demonstrated by the graphic below but keep in mind that
overlays operate at the DataFrame level, not on individual geometries, and the
properties from both are retained. In effect, for every shape in the left
:class:`~geopandas.GeoDataFrame`, this operation is executed against every other shape in the right
:class:`~geopandas.GeoDataFrame`:

.. image:: ../../_static/overlay_operations.png

**Source: QGIS Documentation**

.. note::
   Note to users familiar with the *shapely* library: :meth:`~geopandas.GeoDataFrame.overlay` can be thought
   of as offering versions of the standard *shapely* set-operations that deal with
   the complexities of applying set operations to two *GeoSeries*. The standard
   *shapely* set-operations are also available as :class:`~geopandas.GeoSeries` methods.


The different Overlay operations
--------------------------------

First, we create some example data:

.. ipython:: python

    from shapely.geometry import Polygon
    polys1 = geopandas.GeoSeries([Polygon([(0,0), (2,0), (2,2), (0,2)]),
                                  Polygon([(2,2), (4,2), (4,4), (2,4)])])
    polys2 = geopandas.GeoSeries([Polygon([(1,1), (3,1), (3,3), (1,3)]),
                                  Polygon([(3,3), (5,3), (5,5), (3,5)])])

    df1 = geopandas.GeoDataFrame({'geometry': polys1, 'df1':[1,2]})
    df2 = geopandas.GeoDataFrame({'geometry': polys2, 'df2':[1,2]})

These two GeoDataFrames have some overlapping areas:

.. ipython:: python

    ax = df1.plot(color='red');
    @savefig overlay_example.png width=5in
    df2.plot(ax=ax, color='green', alpha=0.5);

We illustrate the different overlay modes with the above example.
The :meth:`~geopandas.GeoDataFrame.overlay` method will determine the set of all individual geometries
from overlaying the two input GeoDataFrames. This result covers the area covered
by the two input GeoDataFrames, and also preserves all unique regions defined by
the combined boundaries of the two GeoDataFrames.

.. note::
   For historical reasons, the overlay method is also available as a top-level function :func:`overlay`.
   It is recommended to use the method as the function may be deprecated in the future.

When using ``how='union'``, all those possible geometries are returned:

.. ipython:: python

    res_union = df1.overlay(df2, how='union')
    res_union

    ax = res_union.plot(alpha=0.5, cmap='tab10')
    df1.plot(ax=ax, facecolor='none', edgecolor='k');
    @savefig overlay_example_union.png width=5in
    df2.plot(ax=ax, facecolor='none', edgecolor='k');

The other ``how`` operations will return different subsets of those geometries.
With ``how='intersection'``, it returns only those geometries that are contained
by both GeoDataFrames:

.. ipython:: python

    res_intersection = df1.overlay(df2, how='intersection')
    res_intersection

    ax = res_intersection.plot(cmap='tab10')
    df1.plot(ax=ax, facecolor='none', edgecolor='k');
    @savefig overlay_example_intersection.png width=5in
    df2.plot(ax=ax, facecolor='none', edgecolor='k');

``how='symmetric_difference'`` is the opposite of ``'intersection'`` and returns
the geometries that are only part of one of the GeoDataFrames but not of both:

.. ipython:: python

    res_symdiff = df1.overlay(df2, how='symmetric_difference')
    res_symdiff

    ax = res_symdiff.plot(cmap='tab10')
    df1.plot(ax=ax, facecolor='none', edgecolor='k');
    @savefig overlay_example_symdiff.png width=5in
    df2.plot(ax=ax, facecolor='none', edgecolor='k');

To obtain the geometries that are part of ``df1`` but are not contained in
``df2``, you can use ``how='difference'``:

.. ipython:: python

    res_difference = df1.overlay(df2, how='difference')
    res_difference

    ax = res_difference.plot(cmap='tab10')
    df1.plot(ax=ax, facecolor='none', edgecolor='k');
    @savefig overlay_example_difference.png width=5in
    df2.plot(ax=ax, facecolor='none', edgecolor='k');

Finally, with ``how='identity'``, the result consists of the surface of ``df1``,
but with the geometries obtained from overlaying ``df1`` with ``df2``:

.. ipython:: python

    res_identity = df1.overlay(df2, how='identity')
    res_identity

    ax = res_identity.plot(cmap='tab10')
    df1.plot(ax=ax, facecolor='none', edgecolor='k');
    @savefig overlay_example_identity.png width=5in
    df2.plot(ax=ax, facecolor='none', edgecolor='k');


Overlay Countries Example
-------------------------

First, we load the countries and cities example datasets and select :

.. ipython:: python

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    capitals = geopandas.read_file(geopandas.datasets.get_path('naturalearth_cities'))

    # Select South America and some columns
    countries = world[world['continent'] == "South America"]
    countries = countries[['geometry', 'name']]

    # Project to crs that uses meters as distance measure
    countries = countries.to_crs('epsg:3395')
    capitals = capitals.to_crs('epsg:3395')

To illustrate the :meth:`~geopandas.GeoDataFrame.overlay` method, consider the following case in which one
wishes to identify the "core" portion of each country -- defined as areas within
500km of a capital -- using a ``GeoDataFrame`` of countries and a
``GeoDataFrame`` of capitals.

.. ipython:: python

    # Look at countries:
    @savefig world_basic.png width=5in
    countries.plot();

    # Now buffer cities to find area within 500km.
    # Check CRS -- World Mercator, units of meters.
    capitals.crs

    # make 500km buffer
    capitals['geometry']= capitals.buffer(500000)
    @savefig capital_buffers.png width=5in
    capitals.plot();


To select only the portion of countries within 500km of a capital, we specify the ``how`` option to be "intersect", which creates a new set of polygons where these two layers overlap:

.. ipython:: python

   country_cores = countries.overlay(capitals, how='intersection')
   @savefig country_cores.png width=5in
   country_cores.plot(alpha=0.5, edgecolor='k', cmap='tab10');

Changing the "how" option allows for different types of overlay operations. For example, if we were interested in the portions of countries *far* from capitals (the peripheries), we would compute the difference of the two.

.. ipython:: python

   country_peripheries = countries.overlay(capitals, how='difference')
   @savefig country_peripheries.png width=5in
   country_peripheries.plot(alpha=0.5, edgecolor='k', cmap='tab10');


.. ipython:: python
    :suppress:

    import matplotlib.pyplot as plt
    plt.close('all')


keep_geom_type keyword
----------------------

In default settings, :meth:`~geopandas.GeoDataFrame.overlay` returns only geometries of the same geometry type as GeoDataFrame
(left one) has, where Polygon and MultiPolygon is considered as a same type (other types likewise).
You can control this behavior using ``keep_geom_type`` option, which is set to
True by default. Once set to False, ``overlay`` will return all geometry types resulting from
selected set-operation. Different types can result for example from intersection of touching geometries,
where two polygons intersects in a line or a point.


More Examples
-------------

A larger set of examples of the use of :meth:`~geopandas.GeoDataFrame.overlay` can be found `here <https://nbviewer.jupyter.org/github/geopandas/geopandas/blob/main/doc/source/gallery/overlays.ipynb>`_



.. toctree::
   :maxdepth: 2
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas


Merging Data
=========================================

There are two ways to combine datasets in *geopandas* -- attribute joins and spatial joins.

In an attribute join, a :class:`GeoSeries` or :class:`GeoDataFrame` is
combined with a regular :class:`pandas.Series` or :class:`pandas.DataFrame` based on a
common variable. This is analogous to normal merging or joining in *pandas*.

In a Spatial Join, observations from two :class:`GeoSeries` or :class:`GeoDataFrame`
are combined based on their spatial relationship to one another.

In the following examples, we use these datasets:

.. ipython:: python

   world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
   cities = geopandas.read_file(geopandas.datasets.get_path('naturalearth_cities'))

   # For attribute join
   country_shapes = world[['geometry', 'iso_a3']]
   country_names = world[['name', 'iso_a3']]

   # For spatial join
   countries = world[['geometry', 'name']]
   countries = countries.rename(columns={'name':'country'})


Appending
---------

Appending :class:`GeoDataFrame` and :class:`GeoSeries` uses pandas :meth:`~pandas.DataFrame.append` methods.
Keep in mind, that appended geometry columns needs to have the same CRS.

.. ipython:: python

    # Appending GeoSeries
    joined = world.geometry.append(cities.geometry)

    # Appending GeoDataFrames
    europe = world[world.continent == 'Europe']
    asia = world[world.continent == 'Asia']
    eurasia = europe.append(asia)


Attribute Joins
----------------

Attribute joins are accomplished using the :meth:`~pandas.DataFrame.merge` method. In general, it is recommended
to use the ``merge()`` method called from the spatial dataset. With that said, the stand-alone
:func:`pandas.merge` function will work if the :class:`GeoDataFrame` is in the ``left`` argument;
if a :class:`~pandas.DataFrame` is in the ``left`` argument and a :class:`GeoDataFrame`
is in the ``right`` position, the result will no longer be a :class:`GeoDataFrame`.

For example, consider the following merge that adds full names to a :class:`GeoDataFrame`
that initially has only ISO codes for each country by merging it with a :class:`~pandas.DataFrame`.

.. ipython:: python

   # `country_shapes` is GeoDataFrame with country shapes and iso codes
   country_shapes.head()

   # `country_names` is DataFrame with country names and iso codes
   country_names.head()

   # Merge with `merge` method on shared variable (iso codes):
   country_shapes = country_shapes.merge(country_names, on='iso_a3')
   country_shapes.head()


Spatial Joins
----------------

In a Spatial Join, two geometry objects are merged based on their spatial relationship to one another.

.. ipython:: python


   # One GeoDataFrame of countries, one of Cities.
   # Want to merge so we can get each city's country.
   countries.head()
   cities.head()

   # Execute spatial join

   cities_with_country = cities.sjoin(countries, how="inner", predicate='intersects')
   cities_with_country.head()


GeoPandas provides two spatial-join functions:

- :meth:`GeoDataFrame.sjoin`: joins based on binary predicates (intersects, contains, etc.)
- :meth:`GeoDataFrame.sjoin_nearest`: joins based on proximity, with the ability to set a maximum search radius.

.. note::
   For historical reasons, both methods are also available as top-level functions :func:`sjoin` and :func:`sjoin_nearest`.
   It is recommended to use methods as the functions may be deprecated in the future.

Binary Predicate Joins
~~~~~~~~~~~~~~~~~~~~~~

Binary predicate joins are available via :meth:`GeoDataFrame.sjoin`.

:meth:`GeoDataFrame.sjoin` has two core arguments: ``how`` and ``predicate``.

**predicate**

The ``predicate`` argument specifies how ``geopandas`` decides whether or not to join the attributes of one
object to another, based on their geometric relationship.

The values for ``predicate`` correspond to the names of geometric binary predicates and depend on the spatial
index implementation.

The default spatial index in ``geopandas`` currently supports the following values for ``predicate`` which are
defined in the
`Shapely documentation <http://shapely.readthedocs.io/en/latest/manual.html#binary-predicates>`__:

* `intersects`
* `contains`
* `within`
* `touches`
* `crosses`
* `overlaps`

**how**

The `how` argument specifies the type of join that will occur and which geometry is retained in the resultant
:class:`GeoDataFrame`. It accepts the following options:

* ``left``: use the index from the first (or `left_df`) :class:`GeoDataFrame` that you provide
  to :meth:`GeoDataFrame.sjoin`; retain only the `left_df` geometry column
* ``right``: use index from second (or `right_df`); retain only the `right_df` geometry column
* ``inner``: use intersection of index values from both :class:`GeoDataFrame`; retain only the `left_df` geometry column

Note more complicated spatial relationships can be studied by combining geometric operations with spatial join.
To find all polygons within a given distance of a point, for example, one can first use the :meth:`~geopandas.GeoSeries.buffer` method to expand each
point into a circle of appropriate radius, then intersect those buffered circles with the polygons in question.

Nearest Joins
~~~~~~~~~~~~~

Proximity-based joins can be done via :meth:`GeoDataFrame.sjoin_nearest`.

:meth:`GeoDataFrame.sjoin_nearest` shares the ``how`` argument with :meth:`GeoDataFrame.sjoin`, and
includes two additional arguments: ``max_distance`` and ``distance_col``.

**max_distance**

The ``max_distance`` argument specifies a maximum search radius for matching geometries. This can have a considerable performance impact in some cases.
If you can, it is highly recommended that you use this parameter.

**distance_col**

If set, the resultant GeoDataFrame will include a column with this name containing the computed distances between an input geometry and the nearest geometry.
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas


Geocoding
==========

``geopandas`` supports geocoding (i.e., converting place names to
location on Earth) through `geopy`_, an optional dependency of ``geopandas``.
The following example shows how to get the
locations of boroughs in New York City, and plots those locations along
with the detailed borough boundary file included within ``geopandas``.

.. _geopy: http://geopy.readthedocs.io/

.. ipython:: python

    boros = geopandas.read_file(geopandas.datasets.get_path("nybb"))
    boros.BoroName
    boro_locations = geopandas.tools.geocode(boros.BoroName)
    boro_locations

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    boros.to_crs("EPSG:4326").plot(ax=ax, color="white", edgecolor="black");
    @savefig boro_centers_over_bounds.png
    boro_locations.plot(ax=ax, color="red");


By default, the :func:`~geopandas.tools.geocode` function uses the
`Photon geocoding API <https://photon.komoot.io>`__.
But a different geocoding service can be specified with the
``provider`` keyword.

The argument to ``provider`` can either be a string referencing geocoding
services, such as ``'google'``, ``'bing'``, ``'yahoo'``, and
``'openmapquest'``, or an instance of a :mod:`Geocoder <geopy.geocoders>` from :mod:`geopy`. See
``geopy.geocoders.SERVICE_TO_GEOCODER`` for the full list.
For many providers, parameters such as API keys need to be passed as
``**kwargs`` in the :func:`~geopandas.tools.geocode` call.

For example, to use the OpenStreetMap Nominatim geocoder, you need to specify
a user agent:

.. code-block:: python

    geopandas.tools.geocode(boros.BoroName, provider='nominatim', user_agent="my-application")

.. attention::

    Please consult the Terms of Service for the chosen provider. The example
    above uses ``'photon'`` (the default), which expects fair usage
    - extensive usage will be throttled.
    (`Photon's Terms of Use <https://photon.komoot.io>`_).
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas


Indexing and Selecting Data
===========================

GeoPandas inherits the standard pandas_ methods for indexing/selecting data. This includes label based indexing with :attr:`~pandas.DataFrame.loc` and integer position based indexing with :attr:`~pandas.DataFrame.iloc`, which apply to both :class:`GeoSeries` and :class:`GeoDataFrame` objects. For more information on indexing/selecting, see the pandas_ documentation.

.. _pandas: http://pandas.pydata.org/pandas-docs/stable/indexing.html

In addition to the standard pandas_ methods, GeoPandas also provides
coordinate based indexing with the :attr:`~GeoDataFrame.cx` indexer, which slices using a bounding
box. Geometries in the :class:`GeoSeries` or :class:`GeoDataFrame` that intersect the
bounding box will be returned.

Using the ``world`` dataset, we can use this functionality to quickly select all
countries whose boundaries extend into the southern hemisphere.

.. ipython:: python

   world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
   southern_world = world.cx[:, :0]
   @savefig world_southern.png
   southern_world.plot(figsize=(10, 3));
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas


.. _missing-empty:

Missing and empty geometries
============================

GeoPandas supports, just like in pandas, the concept of missing values (NA
or null values). But for geometry values, we have an additional concept of
empty geometries:

- **Empty geometries** are actual geometry objects but that have no coordinates
  (and thus also no area, for example). They can for example originate from
  taking the intersection of two polygons that have no overlap.
  The scalar object (when accessing a single element of a GeoSeries) is still
  a Shapely geometry object.
- **Missing geometries** are unknown values in a GeoSeries. They will typically
  be propagated in operations (for example in calculations of the area or of
  the intersection), or ignored in reductions such as :attr:`~GeoSeries.unary_union`.
  The scalar object (when accessing a single element of a GeoSeries) is the
  Python ``None`` object.

.. warning::

    Starting from GeoPandas v0.6.0, those two concepts are more consistently
    separated. See :ref:`below <missing-empty.changes-0.6.0>` for more details
    on what changed compared to earlier versions.


Consider the following example GeoSeries with one polygon, one missing value
and one empty polygon:

.. ipython:: python

    from shapely.geometry import Polygon
    s = geopandas.GeoSeries([Polygon([(0, 0), (1, 1), (0, 1)]), None, Polygon([])])
    s

In spatial operations, missing geometries will typically propagate (be missing
in the result as well), while empty geometries are treated as a geometry
and the result will depend on the operation:

.. ipython:: python

    s.area
    s.union(Polygon([(0, 0), (0, 1), (1, 1), (1, 0)]))
    s.intersection(Polygon([(0, 0), (0, 1), (1, 1), (1, 0)]))

The :meth:`GeoSeries.isna` method will only check for missing values and not
for empty geometries:

.. ipython:: python
    :okwarning:

    s.isna()

On the other hand, if you want to know which values are empty geometries,
you can use the :attr:`GeoSeries.is_empty` attribute:

.. ipython:: python

    s.is_empty

To get only the actual geometry objects that are neither missing nor empty,
you can use a combination of both:

.. ipython:: python
    :okwarning:

    s.is_empty | s.isna()
    s[~(s.is_empty | s.isna())]


.. _missing-empty.changes-0.6.0:

Changes since GeoPandas v0.6.0
------------------------------

In GeoPandas v0.6.0, the missing data handling was refactored and made more
consistent across the library.

Historically, missing ("NA") values in a GeoSeries could be represented by empty
geometric objects, in addition to standard representations such as ``None`` and
``np.nan``. At least, this was the case in :meth:`GeoSeries.isna` or when a
GeoSeries got aligned in geospatial operations. But, other methods like
:meth:`~GeoSeries.dropna` and :meth:`~GeoSeries.fillna` did not follow this
approach and did not consider empty geometries as missing.

In GeoPandas v0.6.0, the most important change is :meth:`GeoSeries.isna` no
longer treating empty as missing:

* Using the small example from above, the old behaviour treated both the
  empty as missing geometry as "missing":

  .. code-block:: python

    >>> s
    0    POLYGON ((0 0, 1 1, 0 1, 0 0))
    1                              None
    2          GEOMETRYCOLLECTION EMPTY
    dtype: object

    >>> s.isna()
    0    False
    1     True
    2     True
    dtype: bool

* Starting from GeoPandas v0.6.0, it will now only see actual missing values
  as missing:

  .. ipython:: python
    :okwarning:

    s.isna()

  For now, when ``isna()`` is called on a GeoSeries with empty geometries,
  a warning is raised to alert the user of the changed behaviour with an
  indication how to solve this.

Additionally, the behaviour of :meth:`GeoSeries.align` changed to use
missing values instead of empty geometries to fill non-matching indexes.
Consider the following small toy example:

.. ipython:: python

    from shapely.geometry import Point
    s1 = geopandas.GeoSeries([Point(0, 0), Point(1, 1)], index=[0, 1])
    s2 = geopandas.GeoSeries([Point(1, 1), Point(2, 2)], index=[1, 2])
    s1
    s2

* Previously, the ``align`` method would use empty geometries to fill
  values:

  .. code-block:: python

    >>> s1_aligned, s2_aligned = s1.align(s2)

    >>> s1_aligned
    0                 POINT (0 0)
    1                 POINT (1 1)
    2    GEOMETRYCOLLECTION EMPTY
    dtype: object

    >>> s2_aligned
    0    GEOMETRYCOLLECTION EMPTY
    1                 POINT (1 1)
    2                 POINT (2 2)
    dtype: object

  This method is used under the hood when performing spatial operations on
  mis-aligned GeoSeries objects:

  .. code-block:: python

    >>> s1.intersection(s2)
    0    GEOMETRYCOLLECTION EMPTY
    1                 POINT (1 1)
    2    GEOMETRYCOLLECTION EMPTY
    dtype: object

* Starting from GeoPandas v0.6.0, :meth:`GeoSeries.align` will use missing
  values to fill in the non-aligned indices, to be consistent with the
  behaviour in pandas:

  .. ipython:: python

    s1_aligned, s2_aligned = s1.align(s2)
    s1_aligned
    s2_aligned

  This has the consequence that spatial operations will also use missing
  values instead of empty geometries, which can have a different behaviour
  depending on the spatial operation:

  .. ipython:: python
    :okwarning:

    s1.intersection(s2)
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas


Managing Projections
=========================================


Coordinate Reference Systems
-----------------------------

The Coordinate Reference System (CRS) is important because the geometric shapes
in a GeoSeries or GeoDataFrame object are simply a collection of coordinates in
an arbitrary space. A CRS tells Python how those coordinates relate to places on
the Earth.

You can find the codes for most commonly used projections from
`www.spatialreference.org <https://spatialreference.org/>`_.

The same CRS can often be referred to in many ways. For example, one of the most
commonly used CRS is the WGS84 latitude-longitude projection. This can be
referred to using the authority code ``"EPSG:4326"``.

*geopandas* can accept anything accepted by :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`:

- CRS WKT string
- An authority string (i.e. "epsg:4326")
- An EPSG integer code (i.e. 4326)
- A :class:`pyproj.CRS <pyproj.crs.CRS>`
- An object with a to_wkt method.
- PROJ string
- Dictionary of PROJ parameters
- PROJ keyword arguments for parameters
- JSON string with PROJ parameters

For reference, a few very common projections and their EPSG codes:

* WGS84 Latitude/Longitude: ``"EPSG:4326"``
* UTM Zones (North): ``"EPSG:32633"``
* UTM Zones (South): ``"EPSG:32733"``


What is the best format to store the CRS information?
-----------------------------------------------------

Generally, WKT or SRID's are preferred over PROJ strings as they can contain more information about a given CRS.
Conversions between WKT and PROJ strings will in most cases cause a loss of information, potentially leading to erroneous transformations. If possible WKT2 should be used.

For more details, see https://proj.org/faq.html#what-is-the-best-format-for-describing-coordinate-reference-systems


Setting a Projection
----------------------

There are two relevant operations for projections: setting a projection and re-projecting.

Setting a projection may be necessary when for some reason *geopandas* has coordinate data (x-y values), but no information about how those coordinates refer to locations in the real world. Setting a projection is how one tells *geopandas* how to interpret coordinates. If no CRS is set, *geopandas* geometry operations will still work, but coordinate transformations will not be possible and exported files may not be interpreted correctly by other software.

Be aware that **most of the time** you don't have to set a projection. Data loaded from a reputable source (using the :func:`geopandas.read_file()` command) *should* always include projection information. You can see an objects current CRS through the :attr:`GeoSeries.crs` attribute.

From time to time, however, you may get data that does not include a projection. In this situation, you have to set the CRS so *geopandas* knows how to interpret the coordinates.

For example, if you convert a spreadsheet of latitudes and longitudes into a
GeoSeries by hand, you would set the projection by passing the WGS84
latitude-longitude CRS to the :meth:`GeoSeries.set_crs` method (or by setting
the :attr:`GeoSeries.crs` attribute):

.. sourcecode:: python

    my_geoseries = my_geoseries.set_crs("EPSG:4326")
    my_geoseries = my_geoseries.set_crs(epsg=4326)


Re-Projecting
----------------

Re-projecting is the process of changing the representation of locations from one coordinate system to another. All projections of locations on the Earth into a two-dimensional plane `are distortions <https://en.wikipedia.org/wiki/Map_projection#Which_projection_is_best.3F>`_, the projection that is best for your application may be different from the projection associated with the data you import. In these cases, data can be re-projected using the :meth:`GeoDataFrame.to_crs` command:

.. ipython:: python

    # load example data
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    # Check original projection
    # (it's Plate Carrée! x-y are long and lat)
    world.crs

    # Visualize
    ax = world.plot()
    @savefig world_starting.png
    ax.set_title("WGS84 (lat/lon)");

    # Reproject to Mercator (after dropping Antartica)
    world = world[(world.name != "Antarctica") & (world.name != "Fr. S. Antarctic Lands")]
    world = world.to_crs("EPSG:3395") # world.to_crs(epsg=3395) would also work
    ax = world.plot()
    @savefig world_reproj.png
    ax.set_title("Mercator");


Projection for multiple geometry columns
----------------------------------------

GeoPandas 0.8 implements support for different projections assigned to different geometry
columns of the same GeoDataFrame. The projection is now stored together with geometries per column (directly
on the GeometryArray level).

Note that if GeometryArray has assigned projection, it is preferred over the
projection passed to GeoSeries or GeoDataFrame during the creation:

.. code-block:: python

   >>> array.crs
   <Geographic 2D CRS: EPSG:4326>
   Name: WGS 84
   Axis Info [ellipsoidal]:
   - Lat[north]: Geodetic latitude (degree)
   - Lon[east]: Geodetic longitude (degree)
   ...
   >>> GeoSeries(array, crs=3395).crs  # crs=3395 is ignored as array already has CRS
   FutureWarning: CRS mismatch between CRS of the passed geometries and 'crs'. Use 'GeoDataFrame.set_crs(crs, allow_override=True)' to overwrite CRS or 'GeoDataFrame.to_crs(crs)' to reproject geometries. CRS mismatch will raise an error in the future versions of GeoPandas.
       GeoSeries(array, crs=3395).crs

   <Geographic 2D CRS: EPSG:4326>
   Name: WGS 84
   Axis Info [ellipsoidal]:
   - Lat[north]: Geodetic latitude (degree)
   - Lon[east]: Geodetic longitude (degree)
   ...

If you want to overwrite projection, you can then assign it to the GeoSeries
manually or re-project geometries to the target projection using either
``GeoSeries.set_crs(epsg=3395, allow_override=True)`` or
``GeoSeries.to_crs(epsg=3395)``.

All GeometryArray-based operations preserve projection; however, if you loop over a column
containing geometry, this information might be lost.


Upgrading to GeoPandas 0.7 with pyproj > 2.2 and PROJ > 6
---------------------------------------------------------

Starting with GeoPandas 0.7, the `.crs` attribute of a GeoSeries or GeoDataFrame
stores the CRS information as a :class:`pyproj.CRS <pyproj.crs.CRS>`, and no longer as a proj4 string
or dict.

Before, you might have seen this:

.. code-block:: python

   >>> gdf.crs
   {'init': 'epsg:4326'}

while now you will see something like this:

.. code-block:: python

   >>> gdf.crs
   <Geographic 2D CRS: EPSG:4326>
   Name: WGS 84
   Axis Info [ellipsoidal]:
   - Lat[north]: Geodetic latitude (degree)
   - Lon[east]: Geodetic longitude (degree)
   ...
   >>> type(gdf.crs)
   pyproj.crs.CRS

This gives a better user interface and integrates improvements from pyproj and
PROJ 6, but might also require some changes in your code. See `this blogpost
<https://jorisvandenbossche.github.io/blog/2020/02/11/geopandas-pyproj-crs/>`__
for some more background, and the subsections below cover different possible
migration issues.

See the `pyproj docs <https://pyproj4.github.io/pyproj/stable/>`__ for more on
the :class:`pyproj.CRS <pyproj.crs.CRS>` object.

Importing data from files
^^^^^^^^^^^^^^^^^^^^^^^^^

When reading geospatial files with :func:`geopandas.read_file`, things should
mostly work out of the box. For example, reading the example countries dataset
yields a proper CRS:

.. ipython:: python

   df = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
   df.crs

However, in certain cases (with older CRS formats), the resulting CRS object
might not be fully as expected. See the :ref:`section below <unrecognized-crs-reasons>`
for possible reasons and how to solve it.


Manually specifying the CRS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

When specifying the CRS manually in your code (e.g., because your data has not
yet a CRS, or when converting to another CRS), this might require a change in
your code.

**"init" proj4 strings/dicts**

Currently, a lot of people (and also the GeoPandas docs showed that before)
specify the EPSG code using the "init" proj4 string:

.. code-block:: python

   ## OLD
   GeoDataFrame(..., crs={'init': 'epsg:4326'})
   # or
   gdf.crs = {'init': 'epsg:4326'}
   # or
   gdf.to_crs({'init': 'epsg:4326'})

The above will now raise a deprecation warning from pyproj, and instead of the
"init" proj4 string, you should use only the EPSG code itself as follows:

.. code-block:: python

   ## NEW
   GeoDataFrame(..., crs="EPSG:4326")
   # or
   gdf.crs = "EPSG:4326"
   # or
   gdf.to_crs("EPSG:4326")


**proj4 strings/dicts**

Although a full proj4 string is not deprecated (as opposed to the "init" string
above), it is still recommended to change it with an EPSG code if possible.

For example, instead of:

.. code-block:: python

   gdf.crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"

we recommend to do:

.. code-block:: python

   gdf.crs = "EPSG:2163"

*if* you know the EPSG code for the projection you are using.

One possible way to find out the EPSG code is using pyproj for this:

.. code-block:: python

   >>> import pyproj
   >>> crs = pyproj.CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
   >>> crs.to_epsg()
   2163

(you might need to set the ``min_confidence`` keyword of ``to_epsg`` to a lower
value if the match is not perfect)

Further, on websites such as `spatialreference.org <https://spatialreference.org/>`__
and `epsg.io <https://epsg.io/>`__ the descriptions of many CRS can be found
including their EPSG codes and proj4 string definitions.

**Other formats**

Next to the EPSG code mentioned above, there are also other ways to specify the
CRS: an actual :class:`pyproj.CRS <pyproj.crs.CRS>` object, a WKT string, a PROJ JSON string, etc.
Anything that is accepted by :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>` can by specified
to the ``crs`` keyword/attribute in GeoPandas.

Also compatible CRS objects, such as from the :mod:`rasterio` package, can be
passed directly to GeoPandas.


The axis order of a CRS
^^^^^^^^^^^^^^^^^^^^^^^

Starting with PROJ 6 / pyproj 2, the axis order of the official EPSG definition
is honoured. For example, when using geographic coordinates (degrees of longitude
and latitude) in the standard EPSG:4326, the CRS will look like:

.. code-block:: python

   >>> pyproj.CRS(3EPSG:4326")
   <Geographic 2D CRS: EPSG:4326>
   ...
   Axis Info [ellipsoidal]:
   - Lat[north]: Geodetic latitude (degree)
   - Lon[east]: Geodetic longitude (degree)
   ...

This mentions the order as (lat, lon), as that is the official order of coordinates
in EPSG:4326. In GeoPandas, however, the coordinates are always stored as (x, y),
and thus as (lon, lat) order, regardless of the CRS (i.e. the "traditional" order used
in GIS). When reprojecting, GeoPandas and pyproj will under the hood take care of
this difference in axis order, so the user doesn't need to care about this.

.. _unrecognized-crs-reasons:

Why is it not properly recognizing my CRS?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are many file sources and CRS definitions out there "in the wild" that
might have a CRS description that does not fully conform to the new standards of
PROJ > 6 (proj4 strings, older WKT formats, ...). In such cases, you will get a
:class:`pyproj.CRS <pyproj.crs.CRS>` object that might not be fully what you expected (e.g. not equal
to the expected EPSG code). Below we list a few possible cases.

I get a "Bound CRS"?
~~~~~~~~~~~~~~~~~~~~

Some CRS definitions include a *"towgs84" clause*, which can give problems in
recognizing the actual CRS.

For example, both the proj4 and WKT representation for EPSG:31370 (the local
projection used in Belgium) as can be found at `https://spatialreference.org/ref/epsg/31370/ <https://spatialreference.org/ref/epsg/31370/>`__
include this. When taking one of those definitions from that site, and creating
a CRS object:

.. code-block:: python

   >>> import pyproj
   >>> crs = pyproj.CRS("+proj=lcc +lat_1=51.16666723333333 +lat_2=49.8333339 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +towgs84=106.869,-52.2978,103.724,-0.33657,0.456955,-1.84218,1 +units=m +no_defs")
   >>> crs
   <Bound CRS: +proj=lcc +lat_1=51.16666723333333 +lat_2=49.83333 ...>
   Name: unknown
   Axis Info [cartesian]:
   - E[east]: Easting (metre)
   - N[north]: Northing (metre)
   Area of Use:
   - undefined
   Coordinate Operation:
   - name: Transformation from unknown to WGS84
   - method: Position Vector transformation (geog2D domain)
   Datum: Unknown based on International 1909 (Hayford) ellipsoid
   - Ellipsoid: International 1909 (Hayford)
   - Prime Meridian: Greenwich
   Source CRS: unknown

You notice that the above is a not a "Projected CRS" as expected, but a "Bound CRS".
This is because it is "bound" to a conversion to WGS84, and will always use this
when reprojecting instead of letting PROJ determine the best conversion.

To get the actual underlying projected CRS, you can use the ``.source_crs`` attribute:

.. code-block:: python

   >>> crs.source_crs
   <Projected CRS: PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unk ...>
   Name: unknown
   ...

Now we have a "Projected CRS", and now it will also recognize the correct EPSG
number:

.. code-block:: python

   >>> crs.to_epsg()

   >>> crs.source_crs.to_epsg()
   31370

I have a different axis order?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As mentioned above, pyproj now honours the axis order of the EPSG definition.
However, proj4 strings or older WKT versions don't specify this correctly, which
can be a reason that the CRS object is not equal to the expected EPSG code.

Consider the following example of a Canadian projected CRS "EPSG:2953". When
constructing the CRS object from the WKT string as provided on
`https://epsg.io/2953 <https://epsg.io/2953>`__:

.. code-block:: python

   >>> crs = pyproj.CRS("""PROJCS["NAD83(CSRS) / New Brunswick Stereographic",
   ...     GEOGCS["NAD83(CSRS)",
   ...         DATUM["NAD83_Canadian_Spatial_Reference_System",
   ...             SPHEROID["GRS 1980",6378137,298.257222101,
   ...                 AUTHORITY["EPSG","7019"]],
   ...             AUTHORITY["EPSG","6140"]],
   ...         PRIMEM["Greenwich",0,
   ...             AUTHORITY["EPSG","8901"]],
   ...         UNIT["degree",0.0174532925199433,
   ...             AUTHORITY["EPSG","9122"]],
   ...         AUTHORITY["EPSG","4617"]],
   ...     PROJECTION["Oblique_Stereographic"],
   ...     PARAMETER["latitude_of_origin",46.5],
   ...     PARAMETER["central_meridian",-66.5],
   ...     PARAMETER["scale_factor",0.999912],
   ...     PARAMETER["false_easting",2500000],
   ...     PARAMETER["false_northing",7500000],
   ...     UNIT["metre",1,
   ...         AUTHORITY["EPSG","9001"]],
   ...     AUTHORITY["EPSG","2953"]]""")

   >>> crs
   <Projected CRS: PROJCS["NAD83(CSRS) / New Brunswick Stereographic" ...>
   Name: NAD83(CSRS) / New Brunswick Stereographic
   Axis Info [cartesian]:
   - E[east]: Easting (metre)
   - N[north]: Northing (metre)
   ...

Although this is the WKT string as found online for "EPSG:2953", this CRS object
does not evaluate equal to this EPSG code:

.. code-block:: python

   >>> crs == "EPSG:2953"
   False

If we construct the CRS object from the EPSG code (truncated output):

.. code-block:: python

   >>> pyproj.CRS("EPSG:2953")
   <Projected CRS: EPSG:2953>
   Name: NAD83(CSRS) / New Brunswick Stereographic
   Axis Info [cartesian]:
   - N[north]: Northing (metre)
   - E[east]: Easting (metre)
   ...

You can see that the CRS object constructed from the WKT string has a "Easting,
Northing" (i.e. x, y) axis order, while the CRS object constructed from the EPSG
code has a (Northing, Easting) axis order.

Only having this difference in axis order is no problem when using the CRS in
GeoPandas, since GeoPandas always uses a (x, y) order to store the data
regardless of the CRS definition. But, you might still want to verify it is
equivalent to the expected EPSG code. By lowering the `min_confidence`, the axis
order will be ignored:

.. code-block:: python

   >>> crs.to_epsg()

   >>> crs.to_epsg(min_confidence=20)
   2953


The ``.crs`` attribute is no longer a dict or string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you relied on the ``.crs`` object being a dict or a string, such code can
be broken given it is now a :class:`pyproj.CRS <pyproj.crs.CRS>` object. But this object actually
provides a more robust interface to get information about the CRS.

For example, if you used the following code to get the EPSG code:

.. code-block:: python

   gdf.crs['init']

This will no longer work. To get the EPSG code from a ``crs`` object, you can use
the :meth:`~pyproj.crs.CRS.to_epsg` method.

Or to check if a CRS was a certain UTM zone:

.. code-block:: python

   '+proj=utm ' in gdf.crs

could be replaced with the more robust check (requires pyproj 2.6+):

.. code-block:: python

   gdf.crs.utm_zone is not None

And there are many other methods available on the :class:`pyproj.CRS <pyproj.crs.CRS>` class to get
information about the CRS.
Re-projecting using GDAL with Rasterio and Fiona
================================================

The simplest method of re-projecting is :meth:`GeoDataFrame.to_crs`.
It uses ``pyproj`` as the engine and transforms the points within the geometries.

These examples demonstrate how to use ``Fiona`` or ``rasterio`` as the engine to re-project your data.
Fiona and rasterio are powered by GDAL and with algorithms that consider the geometry instead of
just the points the geometry contains. This is particularly useful for antimeridian cutting.
However, this also means the transformation is not as fast.


Fiona Example
-------------

.. code-block:: python

    from functools import partial

    import fiona
    import geopandas
    from fiona.transform import transform_geom
    from packaging import version
    from pyproj import CRS
    from pyproj.enums import WktVersion
    from shapely.geometry import mapping, shape


    # set up Fiona transformer
    def crs_to_fiona(proj_crs):
        proj_crs = CRS.from_user_input(proj_crs)
        if version.parse(fiona.__gdal_version__) < version.parse("3.0.0"):
            fio_crs = proj_crs.to_wkt(WktVersion.WKT1_GDAL)
        else:
            # GDAL 3+ can use WKT2
            fio_crs = proj_crs.to_wkt()
        return fio_crs

    def base_transformer(geom, src_crs, dst_crs):
        return shape(
            transform_geom(
                src_crs=crs_to_fiona(src_crs),
                dst_crs=crs_to_fiona(dst_crs),
                geom=mapping(geom),
                antimeridian_cutting=True,
            )
        )

    # load example data
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    destination_crs = "EPSG:3395"
    forward_transformer = partial(base_transformer, src_crs=world.crs, dst_crs=destination_crs)

    # Reproject to Mercator (after dropping Antartica)
    world = world[(world.name != "Antarctica") & (world.name != "Fr. S. Antarctic Lands")]
    with fiona.Env(OGR_ENABLE_PARTIAL_REPROJECTION="YES"):
        mercator_world = world.set_geometry(world.geometry.apply(forward_transformer), crs=destination_crs)


Rasterio Example
----------------

This example requires rasterio 1.2+ and GDAL 3+.


.. code-block:: python

    import geopandas
    import rasterio.warp
    from shapely.geometry import shape

    # load example data
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    # Reproject to Mercator (after dropping Antartica)
    world = world[(world.name != "Antarctica") & (world.name != "Fr. S. Antarctic Lands")]

    destination_crs = "EPSG:3395"
    geometry = rasterio.warp.transform_geom(
        src_crs=world.crs,
        dst_crs=destination_crs,
        geom=world.geometry.values,
    )
    mercator_world = world.set_geometry(
        [shape(geom) for geom in geometry],
        crs=destination_crs,
    )
.. _geometric_manipulations:

Geometric Manipulations
========================

*geopandas* makes available all the tools for geometric manipulations in the `shapely library <http://shapely.readthedocs.io/en/latest/manual.html>`_.

Note that documentation for all set-theoretic tools for creating new shapes using the relationship between two different spatial datasets -- like creating intersections, or differences -- can be found on the :doc:`set operations <set_operations>` page.

Constructive Methods
~~~~~~~~~~~~~~~~~~~~

.. method:: GeoSeries.buffer(distance, resolution=16)

  Returns a :class:`~geopandas.GeoSeries` of geometries representing all points within a given `distance`
  of each geometric object.

.. attribute:: GeoSeries.boundary

  Returns a :class:`~geopandas.GeoSeries` of lower dimensional objects representing
  each geometries's set-theoretic `boundary`.

.. attribute:: GeoSeries.centroid

  Returns a :class:`~geopandas.GeoSeries` of points for each geometric centroid.

.. attribute:: GeoSeries.convex_hull

  Returns a :class:`~geopandas.GeoSeries` of geometries representing the smallest
  convex `Polygon` containing all the points in each object unless the
  number of points in the object is less than three. For two points,
  the convex hull collapses to a `LineString`; for 1, a `Point`.

.. attribute:: GeoSeries.envelope

  Returns a :class:`~geopandas.GeoSeries` of geometries representing the point or
  smallest rectangular polygon (with sides parallel to the coordinate
  axes) that contains each object.

.. method:: GeoSeries.simplify(tolerance, preserve_topology=True)

  Returns a :class:`~geopandas.GeoSeries` containing a simplified representation of
  each object.

.. attribute:: GeoSeries.unary_union

  Return a geometry containing the union of all geometries in the :class:`~geopandas.GeoSeries`.


Affine transformations
~~~~~~~~~~~~~~~~~~~~~~~~

.. method:: GeoSeries.affine_transform(self, matrix)

  Transform the geometries of the :class:`~geopandas.GeoSeries` using an affine transformation matrix

.. method:: GeoSeries.rotate(self, angle, origin='center', use_radians=False)

  Rotate the coordinates of the :class:`~geopandas.GeoSeries`.

.. method:: GeoSeries.scale(self, xfact=1.0, yfact=1.0, zfact=1.0, origin='center')

 Scale the geometries of the :class:`~geopandas.GeoSeries` along each (x, y, z) dimension.

.. method:: GeoSeries.skew(self, angle, origin='center', use_radians=False)

  Shear/Skew the geometries of the :class:`~geopandas.GeoSeries` by angles along x and y dimensions.

.. method:: GeoSeries.translate(self, xoff=0.0, yoff=0.0, zoff=0.0)

  Shift the coordinates of the :class:`~geopandas.GeoSeries`.



Examples of Geometric Manipulations
------------------------------------

.. sourcecode:: python

    >>> import geopandas
    >>> from geopandas import GeoSeries
    >>> from shapely.geometry import Polygon
    >>> p1 = Polygon([(0, 0), (1, 0), (1, 1)])
    >>> p2 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    >>> p3 = Polygon([(2, 0), (3, 0), (3, 1), (2, 1)])
    >>> g = GeoSeries([p1, p2, p3])
    >>> g
    0         POLYGON ((0 0, 1 0, 1 1, 0 0))
    1    POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))
    2    POLYGON ((2 0, 3 0, 3 1, 2 1, 2 0))
    dtype: geometry

.. image:: ../../_static/test.png

Some geographic operations return normal pandas object.  The :attr:`~geopandas.GeoSeries.area` property of a :class:`~geopandas.GeoSeries` will return a :class:`pandas.Series` containing the area of each item in the :class:`~geopandas.GeoSeries`:

.. sourcecode:: python

    >>> print(g.area)
    0    0.5
    1    1.0
    2    1.0
    dtype: float64

Other operations return GeoPandas objects:

.. sourcecode:: python

    >>> g.buffer(0.5)
    0    POLYGON ((-0.3535533905932737 0.35355339059327...
    1    POLYGON ((-0.5 0, -0.5 1, -0.4975923633360985 ...
    2    POLYGON ((1.5 0, 1.5 1, 1.502407636663901 1.04...
    dtype: geometry

.. image:: ../../_static/test_buffer.png

GeoPandas objects also know how to plot themselves.  GeoPandas uses `matplotlib`_ for plotting. To generate a plot of our GeoSeries, use:

.. sourcecode:: python

    >>> g.plot()

GeoPandas also implements alternate constructors that can read any data format recognized by `fiona`_.  To read a zip file containing an ESRI shapefile with the `borough boundaries of New York City`_ (GeoPandas includes this as an example dataset):

.. sourcecode:: python

    >>> nybb_path = geopandas.datasets.get_path('nybb')
    >>> boros = geopandas.read_file(nybb_path)
    >>> boros.set_index('BoroCode', inplace=True)
    >>> boros.sort_index(inplace=True)
    >>> boros
                   BoroName     Shape_Leng    Shape_Area  \
    BoroCode
    1             Manhattan  359299.096471  6.364715e+08
    2                 Bronx  464392.991824  1.186925e+09
    3              Brooklyn  741080.523166  1.937479e+09
    4                Queens  896344.047763  3.045213e+09
    5         Staten Island  330470.010332  1.623820e+09

                                                       geometry
    BoroCode
    1         MULTIPOLYGON (((981219.0557861328 188655.31579...
    2         MULTIPOLYGON (((1012821.805786133 229228.26458...
    3         MULTIPOLYGON (((1021176.479003906 151374.79699...
    4         MULTIPOLYGON (((1029606.076599121 156073.81420...
    5         MULTIPOLYGON (((970217.0223999023 145643.33221...

.. image:: ../../_static/nyc.png

.. sourcecode:: python

    >>> boros['geometry'].convex_hull
    BoroCode
    1    POLYGON ((977855.4451904297 188082.3223876953,...
    2    POLYGON ((1017949.977600098 225426.8845825195,...
    3    POLYGON ((988872.8212280273 146772.0317993164,...
    4    POLYGON ((1000721.531799316 136681.776184082, ...
    5    POLYGON ((915517.6877458114 120121.8812543372,...
    dtype: geometry

.. image:: ../../_static/nyc_hull.png

To demonstrate a more complex operation, we'll generate a
:class:`~geopandas.GeoSeries` containing 2000 random points:

.. sourcecode:: python

    >>> import numpy as np
    >>> from shapely.geometry import Point
    >>> xmin, xmax, ymin, ymax = 900000, 1080000, 120000, 280000
    >>> xc = (xmax - xmin) * np.random.random(2000) + xmin
    >>> yc = (ymax - ymin) * np.random.random(2000) + ymin
    >>> pts = GeoSeries([Point(x, y) for x, y in zip(xc, yc)])

Now draw a circle with fixed radius around each point:

.. sourcecode:: python

    >>> circles = pts.buffer(2000)

We can collapse these circles into a single :class:`MultiPolygon`
geometry with

.. sourcecode:: python

    >>> mp = circles.unary_union

To extract the part of this geometry contained in each borough, we can
just use:

.. sourcecode:: python

    >>> holes = boros['geometry'].intersection(mp)

.. image:: ../../_static/holes.png

and to get the area outside of the holes:

.. sourcecode:: python

    >>> boros_with_holes = boros['geometry'].difference(mp)

.. image:: ../../_static/boros_with_holes.png

Note that this can be simplified a bit, since ``geometry`` is
available as an attribute on a :class:`~geopandas.GeoDataFrame`, and the
:meth:`~geopandas.GeoSeries.intersection` and :meth:`~geopandas.GeoSeries.difference` methods are implemented with the
"&" and "-" operators, respectively.  For example, the latter could
have been expressed simply as ``boros.geometry - mp``.

It's easy to do things like calculate the fractional area in each
borough that are in the holes:

.. sourcecode:: python

    >>> holes.area / boros.geometry.area
    BoroCode
    1    0.579939
    2    0.586833
    3    0.608174
    4    0.582172
    5    0.558075
    dtype: float64

.. _matplotlib: http://matplotlib.org
.. _fiona: http://fiona.readthedocs.io/en/latest/
.. _geopy: https://github.com/geopy/geopy
.. _geo_interface: https://gist.github.com/sgillies/2217756
.. _borough boundaries of New York City: https://data.cityofnewyork.us/City-Government/Borough-Boundaries/tqmj-j8zm

.. toctree::
   :maxdepth: 2
.. _io:

Reading and Writing Files
=========================

Reading Spatial Data
---------------------

*geopandas* can read almost any vector-based spatial data format including ESRI
shapefile, GeoJSON files and more using the command::

    geopandas.read_file()

which returns a GeoDataFrame object. This is possible because *geopandas* makes
use of the great `fiona <http://fiona.readthedocs.io/en/latest/manual.html>`_
library, which in turn makes use of a massive open-source program called
`GDAL/OGR <http://www.gdal.org/>`_ designed to facilitate spatial data
transformations.

Any arguments passed to :func:`geopandas.read_file` after the file name will be
passed directly to :func:`fiona.open`, which does the actual data importation. In
general, :func:`geopandas.read_file` is pretty smart and should do what you want
without extra arguments, but for more help, type::

    import fiona; help(fiona.open)

Among other things, one can explicitly set the driver (shapefile, GeoJSON) with
the ``driver`` keyword, or pick a single layer from a multi-layered file with
the ``layer`` keyword::

    countries_gdf = geopandas.read_file("package.gpkg", layer='countries')

Where supported in :mod:`fiona`, *geopandas* can also load resources directly from
a web URL, for example for GeoJSON files from `geojson.xyz <http://geojson.xyz/>`_::

    url = "http://d2ad6b4ur7yvpq.cloudfront.net/naturalearth-3.3.0/ne_110m_land.geojson"
    df = geopandas.read_file(url)

You can also load ZIP files that contain your data::

    zipfile = "zip:///Users/name/Downloads/cb_2017_us_state_500k.zip"
    states = geopandas.read_file(zipfile)

If the dataset is in a folder in the ZIP file, you have to append its name::

    zipfile = "zip:///Users/name/Downloads/gadm36_AFG_shp.zip!data"

If there are multiple datasets in a folder in the ZIP file, you also have to
specify the filename::

    zipfile = "zip:///Users/name/Downloads/gadm36_AFG_shp.zip!data/gadm36_AFG_1.shp"

It is also possible to read any file-like objects with a :func:`os.read` method, such
as a file handler (e.g. via built-in :func:`open` function) or :class:`~io.StringIO`::

    filename = "test.geojson"
    file = open(filename)
    df = geopandas.read_file(file)

File-like objects from `fsspec <https://filesystem-spec.readthedocs.io/en/latest>`_
can also be used to read data, allowing for any combination of storage backends and caching
supported by that project::

    path = "simplecache::http://download.geofabrik.de/antarctica-latest-free.shp.zip"
    with fsspec.open(path) as file:
        df = geopandas.read_file(file)

You can also read path objects::

    import pathlib
    path_object = pathlib.path(filename)
    df = geopandas.read_file(path_object)

Reading subsets of the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since geopandas is powered by Fiona, which is powered by GDAL, you can take advantage of
pre-filtering when loading in larger datasets. This can be done geospatially with a geometry
or bounding box. You can also filter rows loaded with a slice. Read more at :func:`geopandas.read_file`.

Geometry Filter
^^^^^^^^^^^^^^^

.. versionadded:: 0.7.0

The geometry filter only loads data that intersects with the geometry.

.. code-block:: python

    gdf_mask = geopandas.read_file(
        geopandas.datasets.get_path("naturalearth_lowres")
    )
    gdf = geopandas.read_file(
        geopandas.datasets.get_path("naturalearth_cities"),
        mask=gdf_mask[gdf_mask.continent=="Africa"],
    )

Bounding Box Filter
^^^^^^^^^^^^^^^^^^^

.. versionadded:: 0.1.0

The bounding box filter only loads data that intersects with the bounding box.

.. code-block:: python

    bbox = (
        1031051.7879884212, 224272.49231459625, 1047224.3104931959, 244317.30894023244
    )
    gdf = geopandas.read_file(
        geopandas.datasets.get_path("nybb"),
        bbox=bbox,
    )

Row Filter
^^^^^^^^^^

.. versionadded:: 0.7.0

Filter the rows loaded in from the file using an integer (for the first n rows)
or a slice object.

.. code-block:: python

    gdf = geopandas.read_file(
        geopandas.datasets.get_path("naturalearth_lowres"),
        rows=10,
    )
    gdf = geopandas.read_file(
        geopandas.datasets.get_path("naturalearth_lowres"),
        rows=slice(10, 20),
    )

Field/Column Filters
^^^^^^^^^^^^^^^^^^^^

Load in a subset of fields from the file:

.. note:: Requires Fiona 1.8+

.. code-block:: python

    gdf = geopandas.read_file(
        geopandas.datasets.get_path("naturalearth_lowres"),
        ignore_fields=["iso_a3", "gdp_md_est"],
    )

Skip loading geometry from the file:

.. note:: Requires Fiona 1.8+
.. note:: Returns :obj:`pandas.DataFrame`

.. code-block:: python

    pdf = geopandas.read_file(
        geopandas.datasets.get_path("naturalearth_lowres"),
        ignore_geometry=True,
    )


Writing Spatial Data
---------------------

GeoDataFrames can be exported to many different standard formats using the
:meth:`geopandas.GeoDataFrame.to_file` method.
For a full list of supported formats, type ``import fiona; fiona.supported_drivers``.

In addition, GeoDataFrames can be uploaded to `PostGIS <https://postgis.net/>`__ database (starting with GeoPandas 0.8)
by using the :meth:`geopandas.GeoDataFrame.to_postgis` method.

.. note::

    GeoDataFrame can contain more field types than supported by most of the file formats. For example tuples or lists
    can be easily stored in the GeoDataFrame, but saving them to e.g. GeoPackage or Shapefile will raise a ValueError.
    Before saving to a file, they need to be converted to a format supported by a selected driver.

**Writing to Shapefile**::

    countries_gdf.to_file("countries.shp")

**Writing to GeoJSON**::

    countries_gdf.to_file("countries.geojson", driver='GeoJSON')

**Writing to GeoPackage**::

    countries_gdf.to_file("package.gpkg", layer='countries', driver="GPKG")
    cities_gdf.to_file("package.gpkg", layer='cities', driver="GPKG")


Spatial databases
-----------------

*geopandas* can also get data from a PostGIS database using the
:func:`geopandas.read_postgis` command.

Writing to PostGIS::

    from sqlalchemy import create_engine
    db_connection_url = "postgresql://myusername:mypassword@myhost:5432/mydatabase";
    engine = create_engine(db_connection_url)
    countries_gdf.to_postgis("countries_table", con=engine)


Apache Parquet and Feather file formats
---------------------------------------

.. versionadded:: 0.8.0

GeoPandas supports writing and reading the Apache Parquet and Feather file
formats.

`Apache Parquet <https://parquet.apache.org/>`__ is an efficient, columnar
storage format (originating from the Hadoop ecosystem). It is a widely used
binary file format for tabular data. The Feather file format is the on-disk
representation of the `Apache Arrow <https://arrow.apache.org/>`__ memory
format, an open standard for in-memory columnar data.

The :func:`geopandas.read_parquet`, :func:`geopandas.read_feather`,
:meth:`GeoDataFrame.to_parquet` and :meth:`GeoDataFrame.to_feather` methods
enable fast roundtrip from GeoPandas to those binary file formats, preserving
the spatial information.

.. warning::

    This is an initial implementation of Parquet file support and
    associated metadata. This is tracking version 0.1.0 of the metadata
    specification at:
    https://github.com/geopandas/geo-arrow-spec

    This metadata specification does not yet make stability promises. As such,
    we do not yet recommend using this in a production setting unless you are
    able to rewrite your Parquet or Feather files.
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas
   import matplotlib
   orig = matplotlib.rcParams['figure.figsize']
   matplotlib.rcParams['figure.figsize'] = [orig[0] * 1.5, orig[1]]
   import matplotlib.pyplot as plt
   plt.close('all')


Mapping and Plotting Tools
=========================================


*geopandas* provides a high-level interface to the matplotlib_ library for making maps. Mapping shapes is as easy as using the :meth:`~GeoDataFrame.plot()` method on a :class:`GeoSeries` or :class:`GeoDataFrame`.

.. _matplotlib: https://matplotlib.org/stable/

Loading some example data:

.. ipython:: python

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    cities = geopandas.read_file(geopandas.datasets.get_path('naturalearth_cities'))

We can now plot those GeoDataFrames:

.. ipython:: python

    # Examine country GeoDataFrame
    world.head()

    # Basic plot, random colors
    @savefig world_randomcolors.png
    world.plot();

Note that in general, any options one can pass to `pyplot <http://matplotlib.org/api/pyplot_api.html>`_ in matplotlib_ (or `style options that work for lines <http://matplotlib.org/api/lines_api.html>`_) can be passed to the :meth:`~GeoDataFrame.plot` method.


Choropleth Maps
-----------------

*geopandas* makes it easy to create Choropleth maps (maps where the color of each shape is based on the value of an associated variable). Simply use the plot command with the ``column`` argument set to the column whose values you want used to assign colors.

.. ipython:: python
   :okwarning:

    # Plot by GDP per capita
    world = world[(world.pop_est>0) & (world.name!="Antarctica")]
    world['gdp_per_cap'] = world.gdp_md_est / world.pop_est
    @savefig world_gdp_per_cap.png
    world.plot(column='gdp_per_cap');


Creating a legend
~~~~~~~~~~~~~~~~~

When plotting a map, one can enable a legend using the ``legend`` argument:

.. ipython:: python

    # Plot population estimates with an accurate legend
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    @savefig world_pop_est.png
    world.plot(column='pop_est', ax=ax, legend=True)

However, the default appearance of the legend and plot axes may not be desirable. One can define the plot axes (with ``ax``) and the legend axes (with ``cax``) and then pass those in to the :meth:`~GeoDataFrame.plot` call. The following example uses ``mpl_toolkits`` to vertically align the plot axes and the legend axes:

.. ipython:: python

    # Plot population estimates with an accurate legend
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig, ax = plt.subplots(1, 1)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    @savefig world_pop_est_fixed_legend_height.png
    world.plot(column='pop_est', ax=ax, legend=True, cax=cax)


And the following example plots the color bar below the map and adds its label using ``legend_kwds``:

.. ipython:: python

    # Plot population estimates with an accurate legend
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    @savefig world_pop_est_horizontal.png
    world.plot(column='pop_est',
               ax=ax,
               legend=True,
               legend_kwds={'label': "Population by Country",
                            'orientation': "horizontal"})


Choosing colors
~~~~~~~~~~~~~~~~

One can also modify the colors used by :meth:`~GeoDataFrame.plot` with the ``cmap`` option (for a full list of colormaps, see the `matplotlib website <http://matplotlib.org/users/colormaps.html>`_):

.. ipython:: python

    @savefig world_gdp_per_cap_red.png
    world.plot(column='gdp_per_cap', cmap='OrRd');


To make the color transparent for when you just want to show the boundary, you have two options. One option is to do ``world.plot(facecolor="none", edgecolor="black")``. However, this can cause a lot of confusion because ``"none"``  and ``None`` are different in the context of using ``facecolor`` and they do opposite things. ``None`` does the "default behavior" based on matplotlib, and if you use it for ``facecolor``, it actually adds a color. The second option is to use ``world.boundary.plot()``. This option is more explicit and clear.:

.. ipython:: python

    @savefig world_gdp_per_cap_transparent.png
    world.boundary.plot();


The way color maps are scaled can also be manipulated with the ``scheme`` option (if you have ``mapclassify`` installed, which can be accomplished via ``conda install -c conda-forge mapclassify``). The ``scheme`` option can be set to any scheme provided by mapclassify (e.g. 'box_plot', 'equal_interval',
'fisher_jenks', 'fisher_jenks_sampled', 'headtail_breaks', 'jenks_caspall', 'jenks_caspall_forced', 'jenks_caspall_sampled', 'max_p_classifier', 'maximum_breaks', 'natural_breaks', 'quantiles', 'percentiles', 'std_mean' or 'user_defined'). Arguments can be passed in classification_kwds dict. See the `mapclassify documentation <https://pysal.org/mapclassify>`_ for further details about these map classification schemes.

.. ipython:: python

    @savefig world_gdp_per_cap_quantiles.png
    world.plot(column='gdp_per_cap', cmap='OrRd', scheme='quantiles');


Missing data
~~~~~~~~~~~~

In some cases one may want to plot data which contains missing values - for some features one simply does not know the value. Geopandas (from the version 0.7) by defaults ignores such features.

.. ipython:: python

    import numpy as np
    world.loc[np.random.choice(world.index, 40), 'pop_est'] = np.nan
    @savefig missing_vals.png
    world.plot(column='pop_est');

However, passing ``missing_kwds`` one can specify the style and label of features containing None or NaN.

.. ipython:: python

    @savefig missing_vals_grey.png
    world.plot(column='pop_est', missing_kwds={'color': 'lightgrey'});

    @savefig missing_vals_hatch.png
    world.plot(
        column="pop_est",
        legend=True,
        scheme="quantiles",
        figsize=(15, 10),
        missing_kwds={
            "color": "lightgrey",
            "edgecolor": "red",
            "hatch": "///",
            "label": "Missing values",
        },
    );

Other map customizations
~~~~~~~~~~~~~~~~~~~~~~~~

Maps usually do not have to have axis labels. You can turn them off using ``set_axis_off()`` or ``axis("off")`` axis methods.

.. ipython:: python

    ax = world.plot()
    @savefig set_axis_off.png
    ax.set_axis_off();

Maps with Layers
-----------------

There are two strategies for making a map with multiple layers -- one more succinct, and one that is a little more flexible.

Before combining maps, however, remember to always ensure they share a common CRS (so they will align).

.. ipython:: python

    # Look at capitals
    # Note use of standard `pyplot` line style options
    @savefig capitals.png
    cities.plot(marker='*', color='green', markersize=5);

    # Check crs
    cities = cities.to_crs(world.crs)

    # Now we can overlay over country outlines
    # And yes, there are lots of island capitals
    # apparently in the middle of the ocean!

**Method 1**

.. ipython:: python

    base = world.plot(color='white', edgecolor='black')
    @savefig capitals_over_countries_1.png
    cities.plot(ax=base, marker='o', color='red', markersize=5);

**Method 2: Using matplotlib objects**

.. ipython:: python

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    # set aspect to equal. This is done automatically
    # when using *geopandas* plot on it's own, but not when
    # working with pyplot directly.
    ax.set_aspect('equal')

    world.plot(ax=ax, color='white', edgecolor='black')
    cities.plot(ax=ax, marker='o', color='red', markersize=5)
    @savefig capitals_over_countries_2.png
    plt.show();

Control the order of multiple layers in a plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When plotting multiple layers, use ``zorder`` to take control of the order of layers being plotted.
The lower the ``zorder`` is, the lower the layer is on the map and vice versa.

Without specified ``zorder``, cities (Points) gets plotted below world (Polygons), following the default order based on geometry types.

.. ipython:: python

    ax = cities.plot(color='k')
    @savefig zorder_default.png
    world.plot(ax=ax);

We can set the ``zorder`` for cities higher than for world to move it of top.

.. ipython:: python

    ax = cities.plot(color='k', zorder=2)
    @savefig zorder_set.png
    world.plot(ax=ax, zorder=1);


Pandas Plots
-----------------

Plotting methods also allow for different plot styles from pandas
along with the default ``geo`` plot. These methods can be accessed using
the ``kind`` keyword argument in :meth:`~GeoDataFrame.plot`, and include:

* ``geo`` for mapping
* ``line`` for line plots
* ``bar`` or ``barh`` for bar plots
* ``hist`` for histogram
* ``box`` for boxplot
* ``kde`` or ``density`` for density plots
* ``area``  for area plots
* ``scatter`` for scatter plots
* ``hexbin`` for hexagonal bin plots
* ``pie`` for pie plots

.. ipython:: python

    gdf = world.head(10)
    @savefig pandas_line_plot.png
    gdf.plot(kind='scatter', x="pop_est", y="gdp_md_est")

You can also create these other plots using the ``GeoDataFrame.plot.<kind>`` accessor methods instead of providing the ``kind`` keyword argument.

.. ipython:: python

    @savefig pandas_bar_plot.png
    gdf.plot.bar()

For more information check out the `pandas documentation <https://pandas.pydata.org/pandas-docs/stable/user_guide/visualization.html>`_.


Other Resources
-----------------
Links to jupyter Notebooks for different mapping tasks:

`Making Heat Maps <http://nbviewer.jupyter.org/gist/perrygeo/c426355e40037c452434>`_


.. ipython:: python
    :suppress:

    matplotlib.rcParams['figure.figsize'] = orig


.. ipython:: python
    :suppress:

    import matplotlib.pyplot as plt
    plt.close('all')
.. currentmodule:: geopandas

.. ipython:: python
   :suppress:

   import geopandas
   import matplotlib
   orig = matplotlib.rcParams['figure.figsize']
   matplotlib.rcParams['figure.figsize'] = [orig[0] * 1.5, orig[1]]



Data Structures
=========================================

GeoPandas implements two main data structures, a :class:`GeoSeries` and a
:class:`GeoDataFrame`.  These are subclasses of :class:`pandas.Series` and
:class:`pandas.DataFrame`, respectively.

GeoSeries
---------

A :class:`GeoSeries` is essentially a vector where each entry in the vector
is a set of shapes corresponding to one observation. An entry may consist
of only one shape (like a single polygon) or multiple shapes that are
meant to be thought of as one observation (like the many polygons that
make up the State of Hawaii or a country like Indonesia).

*geopandas* has three basic classes of geometric objects (which are actually *shapely* objects):

* Points / Multi-Points
* Lines / Multi-Lines
* Polygons / Multi-Polygons

Note that all entries in  a :class:`GeoSeries` need not be of the same geometric type, although certain export operations will fail if this is not the case.

Overview of Attributes and Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`GeoSeries` class implements nearly all of the attributes and
methods of Shapely objects.  When applied to a :class:`GeoSeries`, they
will apply elementwise to all geometries in the series.  Binary
operations can be applied between two :class:`GeoSeries`, in which case the
operation is carried out elementwise.  The two series will be aligned
by matching indices.  Binary operations can also be applied to a
single geometry, in which case the operation is carried out for each
element of the series with that geometry.  In either case, a
:class:`~pandas.Series` or a :class:`GeoSeries` will be returned, as appropriate.

A short summary of a few attributes and methods for GeoSeries is
presented here, and a full list can be found in the :doc:`all attributes and methods page <../reference/geoseries>`.
There is also a family of methods for creating new shapes by expanding
existing shapes or applying set-theoretic operations like "union" described
in :doc:`geometric manipulations <geometric_manipulations>`.

Attributes
^^^^^^^^^^^^^^^
* :attr:`~GeoSeries.area`: shape area (units of projection -- see :doc:`projections <projections>`)
* :attr:`~GeoSeries.bounds`: tuple of max and min coordinates on each axis for each shape
* :attr:`~GeoSeries.total_bounds`: tuple of max and min coordinates on each axis for entire GeoSeries
* :attr:`~GeoSeries.geom_type`: type of geometry.
* :attr:`~GeoSeries.is_valid`: tests if coordinates make a shape that is reasonable geometric shape (`according to this <http://www.opengeospatial.org/standards/sfa>`_).

Basic Methods
^^^^^^^^^^^^^^

* :meth:`~GeoSeries.distance`: returns :class:`~pandas.Series` with minimum distance from each entry to ``other``
* :attr:`~GeoSeries.centroid`: returns :class:`GeoSeries` of centroids
* :meth:`~GeoSeries.representative_point`:  returns :class:`GeoSeries` of points that are guaranteed to be within each geometry. It does **NOT** return centroids.
* :meth:`~GeoSeries.to_crs`: change coordinate reference system. See :doc:`projections <projections>`
* :meth:`~GeoSeries.plot`: plot :class:`GeoSeries`. See :doc:`mapping <mapping>`.

Relationship Tests
^^^^^^^^^^^^^^^^^^^

* :meth:`~GeoSeries.geom_almost_equals`: is shape almost the same as ``other`` (good when floating point precision issues make shapes slightly different)
* :meth:`~GeoSeries.contains`: is shape contained within ``other``
* :meth:`~GeoSeries.intersects`: does shape intersect ``other``


GeoDataFrame
------------

A :class:`GeoDataFrame` is a tabular data structure that contains a :class:`GeoSeries`.

The most important property of a :class:`GeoDataFrame` is that it always has one :class:`GeoSeries` column that holds a special status. This :class:`GeoSeries` is referred to as the :class:`GeoDataFrame`'s "geometry". When a spatial method is applied to a :class:`GeoDataFrame` (or a spatial attribute like ``area`` is called), this commands will always act on the "geometry" column.

The "geometry" column -- no matter its name -- can be accessed through the :attr:`~GeoDataFrame.geometry` attribute (``gdf.geometry``), and the name of the ``geometry`` column can be found by typing ``gdf.geometry.name``.

A :class:`GeoDataFrame` may also contain other columns with geometrical (shapely) objects, but only one column can be the active geometry at a time. To change which column is the active geometry column, use the :meth:`GeoDataFrame.set_geometry` method.

An example using the ``worlds`` GeoDataFrame:

.. ipython:: python

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    world.head()
    #Plot countries
    @savefig world_borders.png
    world.plot();

Currently, the column named "geometry" with country borders is the active
geometry column:

.. ipython:: python

    world.geometry.name

We can also rename this column to "borders":

.. ipython:: python

    world = world.rename(columns={'geometry': 'borders'}).set_geometry('borders')
    world.geometry.name

Now, we create centroids and make it the geometry:

.. ipython:: python
   :okwarning:

    world['centroid_column'] = world.centroid
    world = world.set_geometry('centroid_column')

    @savefig world_centroids.png
    world.plot();


**Note:** A :class:`GeoDataFrame` keeps track of the active column by name, so if you rename the active geometry column, you must also reset the geometry::

    gdf = gdf.rename(columns={'old_name': 'new_name'}).set_geometry('new_name')

**Note 2:** Somewhat confusingly, by default when you use the :func:`~geopandas.read_file` command, the column containing spatial objects from the file is named "geometry" by default, and will be set as the active geometry column. However, despite using the same term for the name of the column and the name of the special attribute that keeps track of the active column, they are distinct. You can easily shift the active geometry column to a different :class:`GeoSeries` with the :meth:`~GeoDataFrame.set_geometry` command. Further, ``gdf.geometry`` will always return the active geometry column, *not* the column named ``geometry``. If you wish to call a column named "geometry", and a different column is the active geometry column, use ``gdf['geometry']``, not ``gdf.geometry``.

Attributes and Methods
~~~~~~~~~~~~~~~~~~~~~~

Any of the attributes calls or methods described for a :class:`GeoSeries` will work on a :class:`GeoDataFrame` -- effectively, they are just applied to the "geometry" :class:`GeoSeries`.

However, :class:`GeoDataFrames <GeoDataFrame>` also have a few extra methods for input and output which are described on the :doc:`Input and Output <io>` page and for geocoding with are described in :doc:`Geocoding <geocoding>`.


.. ipython:: python
    :suppress:

    matplotlib.rcParams['figure.figsize'] = orig


Display options
---------------

GeoPandas has an ``options`` attribute with currently a single configuration
option to control:

.. ipython:: python

    import geopandas
    geopandas.options

The ``geopandas.options.display_precision`` option can control the number of
decimals to show in the display of coordinates in the geometry column.
In the ``world`` example of above, the default is to show 5 decimals for
geographic coordinates:

.. ipython:: python

    world['centroid_column'].head()

If you want to change this, for example to see more decimals, you can do:

.. ipython:: python

    geopandas.options.display_precision = 9
    world['centroid_column'].head()
.. ipython:: python
   :suppress:

   import geopandas
   import matplotlib
   orig = matplotlib.rcParams['figure.figsize']
   matplotlib.rcParams['figure.figsize'] = [orig[0] * 1.5, orig[1]]


Aggregation with dissolve
=============================

Spatial data are often more granular than we need. For example, we might have data on sub-national units, but we're actually interested in studying patterns at the level of countries.

In a non-spatial setting, when all we need are summary statistics of the data, we aggregate our data using the :meth:`~pandas.DataFrame.groupby` function. But for spatial data, we sometimes also need to aggregate geometric features. In the *geopandas* library, we can aggregate geometric features using the :meth:`~geopandas.GeoDataFrame.dissolve` function.

:meth:`~geopandas.GeoDataFrame.dissolve` can be thought of as doing three things:

(a) it dissolves all the geometries within a given group together into a single geometric feature (using the :attr:`~geopandas.GeoSeries.unary_union` method), and
(b) it aggregates all the rows of data in a group using :ref:`groupby.aggregate <groupby.aggregate>`, and
(c) it combines those two results.

:meth:`~geopandas.GeoDataFrame.dissolve` Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we are interested in studying continents, but we only have country-level data like the country dataset included in *geopandas*. We can easily convert this to a continent-level dataset.


First, let's look at the most simple case where we just want continent shapes and names. By default, :meth:`~geopandas.GeoDataFrame.dissolve` will pass ``'first'`` to :ref:`groupby.aggregate <groupby.aggregate>`.

.. ipython:: python

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    world = world[['continent', 'geometry']]
    continents = world.dissolve(by='continent')

    @savefig continents1.png
    continents.plot();

    continents.head()

If we are interested in aggregate populations, however, we can pass different functions to the :meth:`~geopandas.GeoDataFrame.dissolve` method to aggregate populations using the ``aggfunc =`` argument:

.. ipython:: python

   world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
   world = world[['continent', 'geometry', 'pop_est']]
   continents = world.dissolve(by='continent', aggfunc='sum')

   @savefig continents2.png
   continents.plot(column = 'pop_est', scheme='quantiles', cmap='YlOrRd');

   continents.head()


.. ipython:: python
    :suppress:

    matplotlib.rcParams['figure.figsize'] = orig


.. toctree::
   :maxdepth: 2

Dissolve Arguments
~~~~~~~~~~~~~~~~~~

The ``aggfunc =`` argument defaults to 'first' which means that the first row of attributes values found in the dissolve routine will be assigned to the resultant dissolved geodataframe.
However it also accepts other summary statistic options as allowed by :meth:`pandas.groupby <pandas.DataFrame.groupby>` including:

* 'first'
* 'last'
* 'min'
* 'max'
* 'sum'
* 'mean'
* 'median'
* function
* string function name
* list of functions and/or function names, e.g. [np.sum, 'mean']
* dict of axis labels -> functions, function names or list of such.

For example, to get the number of contries on each continent, 
as well as the populations of the largest and smallest country of each,
we can aggregate the ``'name'`` column using ``'count'``,
and the ``'pop_est'`` column using ``'min'`` and ``'max'``:

.. ipython:: python

    world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
    continents = world.dissolve(
        by="continent",
        aggfunc={
            "name": "count",
            "pop_est": ["min", "max"],
        },
    )

   continents.head()=====
Tools
=====
.. currentmodule:: geopandas

.. autosummary::
   :toctree: api/

   sjoin
   sjoin_nearest
   overlay
   clip
   tools.geocode
   tools.reverse_geocode
   tools.collect
   points_from_xy
   datasets.available
   datasets.get_path
=============
Spatial index
=============
.. currentmodule:: geopandas

GeoPandas offers built-in support for spatial indexing using an R-Tree algorithm.
Depending on the ability to import ``pygeos``, GeoPandas will either use
``pygeos.STRtree`` or ``rtree.index.Index``. The main interface for both is the
same and follows the ``pygeos`` model.

``GeoSeries.sindex`` creates a spatial index, which can use the methods and
properties documented below.

Constructor
-----------
.. autosummary::
   :toctree: api/

    GeoSeries.sindex

Spatial Index object
--------------------

The spatial index object returned from :attr:`GeoSeries.sindex` has the following
methods:

.. currentmodule:: geopandas.sindex.SpatialIndex
.. autosummary::
   :toctree: api/

    intersection
    is_empty
    nearest
    query
    query_bulk
    size
    valid_query_predicates

The concrete implementations currently available are
``geopandas.sindex.PyGEOSSTRTreeIndex`` and ``geopandas.sindex.RTreeIndex``.

In addition to the methods listed above, the ``rtree``-based spatial index
(``geopandas.sindex.RTreeIndex``) offers the full capability of
``rtree.index.Index`` - see the full API in the `rtree documentation`_.

Similarly, the ``pygeos``-based spatial index
(``geopandas.sindex.PyGEOSSTRTreeIndex``) offers the full capability of
``pygeos.STRtree``, including nearest-neighbor queries.
See the full API in the `PyGEOS STRTree documentation`_.

.. _rtree documentation: https://rtree.readthedocs.io/en/stable/class.html
.. _PyGEOS STRTree documentation: https://pygeos.readthedocs.io/en/latest/strtree.html
=========
GeoSeries
=========
.. currentmodule:: geopandas

Constructor
-----------
.. autosummary::
   :toctree: api/

   GeoSeries

General methods and attributes
------------------------------

.. autosummary::
   :toctree: api/

   GeoSeries.area
   GeoSeries.boundary
   GeoSeries.bounds
   GeoSeries.total_bounds
   GeoSeries.length
   GeoSeries.geom_type
   GeoSeries.distance
   GeoSeries.representative_point
   GeoSeries.exterior
   GeoSeries.interiors
   GeoSeries.x
   GeoSeries.y
   GeoSeries.z

Unary predicates
----------------

.. autosummary::
   :toctree: api/

   GeoSeries.is_empty
   GeoSeries.is_ring
   GeoSeries.is_simple
   GeoSeries.is_valid
   GeoSeries.has_z


Binary Predicates
-----------------

.. autosummary::
   :toctree: api/

   GeoSeries.contains
   GeoSeries.crosses
   GeoSeries.disjoint
   GeoSeries.geom_equals
   GeoSeries.geom_almost_equals
   GeoSeries.geom_equals_exact
   GeoSeries.intersects
   GeoSeries.overlaps
   GeoSeries.touches
   GeoSeries.within
   GeoSeries.covers
   GeoSeries.covered_by


Set-theoretic Methods
---------------------

.. autosummary::
   :toctree: api/

   GeoSeries.difference
   GeoSeries.intersection
   GeoSeries.symmetric_difference
   GeoSeries.union

Constructive Methods and Attributes
-----------------------------------

.. autosummary::
   :toctree: api/

   GeoSeries.buffer
   GeoSeries.boundary
   GeoSeries.centroid
   GeoSeries.convex_hull
   GeoSeries.envelope
   GeoSeries.simplify

Affine transformations
----------------------

.. autosummary::
   :toctree: api/

   GeoSeries.affine_transform
   GeoSeries.rotate
   GeoSeries.scale
   GeoSeries.skew
   GeoSeries.translate

Aggregating and exploding
-------------------------

.. autosummary::
   :toctree: api/

   GeoSeries.unary_union
   GeoSeries.explode

Serialization / IO / conversion
-------------------------------

.. autosummary::
   :toctree: api/

   GeoSeries.from_file
   GeoSeries.from_wkb
   GeoSeries.from_wkt
   GeoSeries.from_xy
   GeoSeries.to_file
   GeoSeries.to_json
   GeoSeries.to_wkb
   GeoSeries.to_wkt

Projection handling
-------------------

.. autosummary::
   :toctree: api/

   GeoSeries.crs
   GeoSeries.set_crs
   GeoSeries.to_crs
   GeoSeries.estimate_utm_crs

Missing values
--------------

.. autosummary::
   :toctree: api/

   GeoSeries.fillna
   GeoSeries.isna
   GeoSeries.notna

Overlay operations
------------------

.. autosummary::
   :toctree: api/

   GeoSeries.clip

Plotting
--------

.. autosummary::
   :toctree: api/

   GeoSeries.plot
   GeoSeries.explore


Spatial index
-------------

.. autosummary::
   :toctree: api/

   GeoSeries.sindex
   GeoSeries.has_sindex

Indexing
--------

.. autosummary::
   :toctree: api/

   GeoSeries.cx

Interface
---------

.. autosummary::
   :toctree: api/

   GeoSeries.__geo_interface__


Methods of pandas ``Series`` objects are also available, although not
all are applicable to geometric objects and some may return a
``Series`` rather than a ``GeoSeries`` result when appropriate. The methods
``isna()`` and ``fillna()`` have been
implemented specifically for ``GeoSeries`` and are expected to work
correctly.
============
GeoDataFrame
============
.. currentmodule:: geopandas

A ``GeoDataFrame`` is a tabular data structure that contains a column
which contains a ``GeoSeries`` storing geometry.

Constructor
-----------
.. autosummary::
   :toctree: api/

   GeoDataFrame

Serialization / IO / conversion
-------------------------------

.. autosummary::
   :toctree: api/

   GeoDataFrame.from_file
   GeoDataFrame.from_features
   GeoDataFrame.from_postgis
   GeoDataFrame.to_file
   GeoDataFrame.to_json
   GeoDataFrame.to_parquet
   GeoDataFrame.to_feather
   GeoDataFrame.to_postgis
   GeoDataFrame.to_wkb
   GeoDataFrame.to_wkt

Projection handling
-------------------

.. autosummary::
   :toctree: api/

   GeoDataFrame.crs
   GeoDataFrame.set_crs
   GeoDataFrame.to_crs
   GeoDataFrame.estimate_utm_crs

Active geometry handling
------------------------

.. autosummary::
   :toctree: api/

   GeoDataFrame.rename_geometry
   GeoDataFrame.set_geometry

Aggregating and exploding
-------------------------

.. autosummary::
   :toctree: api/

   GeoDataFrame.dissolve
   GeoDataFrame.explode

Spatial joins
-------------

.. autosummary::
   :toctree: api/

   GeoDataFrame.sjoin
   GeoDataFrame.sjoin_nearest

Overlay operations
------------------

.. autosummary::
   :toctree: api/

   GeoDataFrame.clip
   GeoDataFrame.overlay

Plotting
--------

.. autosummary::
   :toctree: api/

   GeoDataFrame.explore


.. autosummary::
   :toctree: api/
   :template: accessor_callable.rst

   GeoDataFrame.plot

Spatial index
-------------

.. autosummary::
   :toctree: api/

   GeoDataFrame.sindex
   GeoDataFrame.has_sindex

Indexing
--------

.. autosummary::
   :toctree: api/

   GeoDataFrame.cx

Interface
---------

.. autosummary::
   :toctree: api/

   GeoDataFrame.__geo_interface__
   GeoDataFrame.iterfeatures

All pandas ``DataFrame`` methods are also available, although they may
not operate in a meaningful way on the ``geometry`` column. All methods
listed in `GeoSeries <geoseries>`__ work directly on an active geometry column of GeoDataFrame.
============
Input/output
============
.. currentmodule:: geopandas

GIS vector files
----------------
.. autosummary::
   :toctree: api/

   read_file
   GeoDataFrame.to_file

PostGIS
-------
.. autosummary::
   :toctree: api/

   read_postgis
   GeoDataFrame.to_postgis


Feather
-------
.. autosummary::
   :toctree: api/

   read_feather
   GeoDataFrame.to_feather

Parquet
-------
.. autosummary::
   :toctree: api/

   read_parquet
   GeoDataFrame.to_parquet
=======
Testing
=======
.. currentmodule:: geopandas

GeoPandas includes specific functions to test its objects.

.. autosummary::
   :toctree: api/

   .. testing.geom_equals
   .. testing.geom_almost_equals
   testing.assert_geoseries_equal
   testing.assert_geodataframe_equal
