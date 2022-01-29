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
