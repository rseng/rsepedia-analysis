# Change Log

## v0.x.x Unreleased
### New features
* Add new convenience function `arviz.extract_dataset` ([1725](https://github.com/arviz-devs/arviz/pull/1725))
* [experimental] Enable dask chunking information to be passed to `InferenceData.from_netcdf` with regex support ([1749](https://github.com/arviz-devs/arviz/pull/1749))
* Allow kwargs to customize appearance of the mean in `plot_lm`

### Maintenance and fixes
* Drop Python 3.6 support ([1430](https://github.com/arviz-devs/arviz/pull/1430))
* Bokeh 3 compatibility. ([1919](https://github.com/arviz-devs/arviz/pull/1919))
* Remove manual setting of 2d KDE limits ([1939](https://github.com/arviz-devs/arviz/pull/1939))
* Pin to bokeh<3 version ([1954](https://github.com/arviz-devs/arviz/pull/1954))
* Fix legend labels in plot_ppc to reflect prior or posterior. ([1967](https://github.com/arviz-devs/arviz/pull/1967))

### Deprecation

### Documentation
* Fixed typo in `Forestplot` documentation
* Restructured contributing section and added several new pages to help contributing to docs ([1903](https://github.com/arviz-devs/arviz/pull/1903))


## v0.11.4 (2021 Oct 3)

### Maintenance and fixes
* Fix standard deviation code in density utils by replacing it with `np.std`. ([1833](https://github.com/arviz-devs/arviz/pull/1833))

## v0.11.3 (2021 Oct 1)
### New features
* Change order of regularization in `psislw` ([1943](https://github.com/arviz-devs/arviz/pull/1943))
* Added `labeller` argument to enable label customization in plots and summary ([1201](https://github.com/arviz-devs/arviz/pull/1201))
* Added `arviz.labels` module with classes and utilities ([1201](https://github.com/arviz-devs/arviz/pull/1201) and [1605](https://github.com/arviz-devs/arviz/pull/1605))
* Added probability estimate within ROPE in `plot_posterior` ([1570](https://github.com/arviz-devs/arviz/pull/1570))
* Added `rope_color` and `ref_val_color` arguments to `plot_posterior` ([1570](https://github.com/arviz-devs/arviz/pull/1570))
* Improved retrieving or pointwise log likelihood in `from_cmdstanpy`, `from_cmdstan` and `from_pystan` ([1579](https://github.com/arviz-devs/arviz/pull/1579) and [1599](https://github.com/arviz-devs/arviz/pull/1599))
* Added interactive legend to bokeh `forestplot` ([1591](https://github.com/arviz-devs/arviz/pull/1591))
* Added interactive legend to bokeh `ppcplot` ([1602](https://github.com/arviz-devs/arviz/pull/1602))
* Add more helpful error message for HDF5 problems reading `InferenceData` from NetCDF ([1637](https://github.com/arviz-devs/arviz/pull/1637))
* Added `data.log_likelihood`, `stats.ic_compare_method` and `plot.density_kind` to `rcParams` ([1611](https://github.com/arviz-devs/arviz/pull/1611))
* Improve error messages in `stats.compare()`, and `var_name` parameter. ([1616](https://github.com/arviz-devs/arviz/pull/1616))
* Added ability to plot HDI contours to `plot_kde` with the new `hdi_probs` parameter. ([1665](https://github.com/arviz-devs/arviz/pull/1665))
* Add dtype parsing and setting in all Stan converters ([1632](https://github.com/arviz-devs/arviz/pull/1632))
* Add option to specify colors for each element in ppc_plot  ([1769](https://github.com/arviz-devs/arviz/pull/1769))

### Maintenance and fixes
* Fix conversion for numpyro models with ImproperUniform latent sites ([1713](https://github.com/arviz-devs/arviz/pull/1713))
* Fixed conversion of Pyro output fit using GPUs ([1659](https://github.com/arviz-devs/arviz/pull/1659))
* Enforced using coordinate values as default labels ([1201](https://github.com/arviz-devs/arviz/pull/1201))
* Integrate `index_origin` with all the library ([1201](https://github.com/arviz-devs/arviz/pull/1201))
* Fix pareto k threshold typo in reloo function ([1580](https://github.com/arviz-devs/arviz/pull/1580))
* Preserve shape from Stan code in `from_cmdstanpy` ([1579](https://github.com/arviz-devs/arviz/pull/1579))
* Updated `from_pystan` converters to follow schema convention ([1585](https://github.com/arviz-devs/arviz/pull/1585)
* Used generator instead of list wherever possible ([1588](https://github.com/arviz-devs/arviz/pull/1588))
* Correctly use chain index when constructing PyMC3 `DefaultTrace` in `from_pymc3` ([1590](https://github.com/arviz-devs/arviz/pull/1590))
* Fix bugs in CmdStanPyConverter ([1595](https://github.com/arviz-devs/arviz/pull/1595) and [1598](https://github.com/arviz-devs/arviz/pull/1598))
* Fix `c` argument in `plot_khat` ([1592](https://github.com/arviz-devs/arviz/pull/1592))
* Fix `ax` argument in `plot_elpd` ([1593](https://github.com/arviz-devs/arviz/pull/1593))
* Remove warning in `stats.py` compare function ([1607](https://github.com/arviz-devs/arviz/pull/1607))
* Fix `ess/rhat` plots in `plot_forest` ([1606](https://github.com/arviz-devs/arviz/pull/1606))
* Fix `from_numpyro` crash when importing model with `thinning=x` for `x > 1` ([1619](https://github.com/arviz-devs/arviz/pull/1619))
* Upload updated mypy.ini in ci if mypy copilot fails ([1624](https://github.com/arviz-devs/arviz/pull/1624))
* Added type checking to raise an error whenever `InferenceData` object is passed using `io_pymc3`'s `trace` argument ([1629](https://github.com/arviz-devs/arviz/pull/1629))
* Fix `xlabels` in `plot_elpd` ([1601](https://github.com/arviz-devs/arviz/pull/1601))
* Renamed `sample` dim to `__sample__` when stacking `chain` and `draw` to avoid dimension collision ([1647](https://github.com/arviz-devs/arviz/pull/1647))
* Removed the `circular` argument in `plot_dist` in favor of `is_circular` ([1681](https://github.com/arviz-devs/arviz/pull/1681))
* Fix `legend` argument in `plot_separation` ([1701](https://github.com/arviz-devs/arviz/pull/1701))
* Removed testing dependency on http download for radon dataset ([1717](https://github.com/arviz-devs/arviz/pull/1717))
* Fixed plot_kde to take labels with kwargs. ([1710](https://github.com/arviz-devs/arviz/pull/1710))
* Fixed xarray related tests. ([1726](https://github.com/arviz-devs/arviz/pull/1726))
* Fix Bokeh deprecation warnings ([1657](https://github.com/arviz-devs/arviz/pull/1657))
* Fix credible inteval percentage in legend in `plot_loo_pit` ([1745](https://github.com/arviz-devs/arviz/pull/1745))
* Arguments `filter_vars` and `filter_groups` now raise `ValueError` if illegal arguments are passed ([1772](https://github.com/arviz-devs/arviz/pull/1772))
* Remove constrained_layout from arviz rcparams ([1764](https://github.com/arviz-devs/arviz/pull/1764))
* Fix plot_elpd for a single outlier ([1787](https://github.com/arviz-devs/arviz/pull/1787))

### Deprecation
* Deprecated `index_origin` and `order` arguments in `az.summary` ([1201](https://github.com/arviz-devs/arviz/pull/1201))

### Documentation
* Language improvements of the first third of the "Label guide" ([1699](https://github.com/arviz-devs/arviz/pull/1699))
* Added "Label guide" page and API section for `arviz.labels` module ([1201](https://github.com/arviz-devs/arviz/pull/1201) and [1635](https://github.com/arviz-devs/arviz/pull/1635))
* Add "Installation guide" page to the documentation ([1551](https://github.com/arviz-devs/arviz/pull/1551))
* Improve documentation on experimental `SamplingWrapper` classes ([1582](https://github.com/arviz-devs/arviz/pull/1582))
* Added example to `plot_hdi` using Inference Data ([1615](https://github.com/arviz-devs/arviz/pull/1615))
* Removed `geweke` diagnostic from `numba` user guide ([1653](https://github.com/arviz-devs/arviz/pull/1653))
* Restructured the documentation sections to improve community and about us information ([1587](https://github.com/arviz-devs/arviz/pull/1587))

## v0.11.2 (2021 Feb 21)
### New features
* Added `to_zarr` and `from_zarr` methods to InferenceData ([1518](https://github.com/arviz-devs/arviz/pull/1518))
* Added confidence interval band to auto-correlation plot ([1535](https://github.com/arviz-devs/arviz/pull/1535))

### Maintenance and fixes
* Updated CmdStanPy converter form compatibility with versions >=0.9.68 ([1558](https://github.com/arviz-devs/arviz/pull/1558) and ([1564](https://github.com/arviz-devs/arviz/pull/1564))
* Updated `from_cmdstanpy`, `from_cmdstan`, `from_numpyro` and `from_pymc3` converters to follow schema convention ([1550](https://github.com/arviz-devs/arviz/pull/1550), [1541](https://github.com/arviz-devs/arviz/pull/1541), [1525](https://github.com/arviz-devs/arviz/pull/1525) and [1555](https://github.com/arviz-devs/arviz/pull/1555))
* Fix calculation of mode as point estimate ([1552](https://github.com/arviz-devs/arviz/pull/1552))
* Remove variable name from legend in posterior predictive plot ([1559](https://github.com/arviz-devs/arviz/pull/1559))
* Added significant digits formatter to round rope values ([1569](https://github.com/arviz-devs/arviz/pull/1569))
* Updated `from_cmdstan`. csv reader, dtype problem fixed and dtype kwarg added for manual dtype casting ([1565](https://github.com/arviz-devs/arviz/pull/1565))

### Deprecation
* Removed Geweke diagnostic ([1545](https://github.com/arviz-devs/arviz/pull/1545))
* Removed credible_interval and include_circ arguments ([1548](https://github.com/arviz-devs/arviz/pull/1548))

### Documentation
* Added an example for converting dataframe to InferenceData ([1556](https://github.com/arviz-devs/arviz/pull/1556))
* Added example for `coords` argument in `plot_posterior` docstring ([1566](https://github.com/arviz-devs/arviz/pull/1566))

## v0.11.1 (2021 Feb 2)
### Maintenance and fixes
* Fixed ovelapping titles and repeating warnings on circular traceplot ([1517](https://github.com/arviz-devs/arviz/pull/1517))
* Removed repetitive variable names from forest plots of multivariate variables ([1527](https://github.com/arviz-devs/arviz/pull/1527))
* Fixed regression in `plot_pair` labels that prevented coord names to be shown when necessary ([1533](https://github.com/arviz-devs/arviz/pull/1533))

### Documentation
* Use tabs in ArviZ example gallery ([1521](https://github.com/arviz-devs/arviz/pull/1521))

## v0.11.0 (2021 Dec 17)
### New features
* Added `to_dataframe` method to InferenceData ([1395](https://github.com/arviz-devs/arviz/pull/1395))
* Added `__getitem__` magic to InferenceData ([1395](https://github.com/arviz-devs/arviz/pull/1395))
* Added group argument to summary ([1408](https://github.com/arviz-devs/arviz/pull/1408))
* Add `ref_line`, `bar`, `vlines` and `marker_vlines` kwargs to `plot_rank` ([1419](https://github.com/arviz-devs/arviz/pull/1419))
* Add observed argument to (un)plot observed data in `plot_ppc` ([1422](https://github.com/arviz-devs/arviz/pull/1422))
* Add support for named dims and coordinates with multivariate observations ([1429](https://github.com/arviz-devs/arviz/pull/1429))
* Add support for discrete variables in rank plots ([1433](https://github.com/arviz-devs/arviz/pull/1433)) and
  `loo_pit` ([1500](https://github.com/arviz-devs/arviz/pull/1500))
* Add `skipna` argument to `plot_posterior` ([1432](https://github.com/arviz-devs/arviz/pull/1432))
* Make stacking the default method to compute weights in `compare` ([1438](https://github.com/arviz-devs/arviz/pull/1438))
* Add `copy()` method to `InferenceData` class. ([1501](https://github.com/arviz-devs/arviz/pull/1501)).


### Maintenance and fixes
* prevent wrapping group names in InferenceData repr_html ([1407](https://github.com/arviz-devs/arviz/pull/1407))
* Updated CmdStanPy interface ([1409](https://github.com/arviz-devs/arviz/pull/1409))
* Remove left out warning about default IC scale in `compare` ([1412](https://github.com/arviz-devs/arviz/pull/1412))
* Fixed a typo found in an error message raised in `distplot.py` ([1414](https://github.com/arviz-devs/arviz/pull/1414))
* Fix typo in `loo_pit` extraction of log likelihood ([1418](https://github.com/arviz-devs/arviz/pull/1418))
* Have `from_pystan` store attrs as strings to allow netCDF storage ([1417](https://github.com/arviz-devs/arviz/pull/1417))
* Remove ticks and spines in `plot_violin`  ([1426 ](https://github.com/arviz-devs/arviz/pull/1426))
* Use circular KDE function and fix tick labels in circular `plot_trace` ([1428](https://github.com/arviz-devs/arviz/pull/1428))
* Fix `pair_plot` for mixed discrete and continuous variables ([1434](https://github.com/arviz-devs/arviz/pull/1434))
* Fix in-sample deviance in `plot_compare` ([1435](https://github.com/arviz-devs/arviz/pull/1435))
* Fix computation of weights in compare ([1438](https://github.com/arviz-devs/arviz/pull/1438))
* Avoid repeated warning in summary ([1442](https://github.com/arviz-devs/arviz/pull/1442))
* Fix hdi failure with boolean array ([1444](https://github.com/arviz-devs/arviz/pull/1444))
* Automatically get the current axes instance for `plt_kde`, `plot_dist` and `plot_hdi` ([1452](https://github.com/arviz-devs/arviz/pull/1452))
* Add grid argument to manually specify the number of rows and columns ([1459](https://github.com/arviz-devs/arviz/pull/1459))
* Switch to `compact=True` by default in our plots ([1468](https://github.com/arviz-devs/arviz/issues/1468))
* `plot_elpd`, avoid modifying the input dict ([1477](https://github.com/arviz-devs/arviz/issues/1477))
* Do not plot divergences in `plot_trace` when `kind=rank_vlines` or `kind=rank_bars` ([1476](https://github.com/arviz-devs/arviz/issues/1476))
* Allow ignoring `observed` argument of `pymc3.DensityDist` in `from_pymc3` ([1495](https://github.com/arviz-devs/arviz/pull/1495))
* Make `from_pymc3` compatible with theano-pymc 1.1.0 ([1495](https://github.com/arviz-devs/arviz/pull/1495))
* Improve typing hints ([1491](https://github.com/arviz-devs/arviz/pull/1491), ([1492](https://github.com/arviz-devs/arviz/pull/1492),
  ([1493](https://github.com/arviz-devs/arviz/pull/1493), ([1494](https://github.com/arviz-devs/arviz/pull/1494) and
  ([1497](https://github.com/arviz-devs/arviz/pull/1497))


### Deprecation
* `plot_khat` deprecate `annotate` argument in favor of `threshold`. The new argument accepts floats ([1478](https://github.com/arviz-devs/arviz/issues/1478))

### Documentation
* Reorganize documentation and change sphinx theme ([1406](https://github.com/arviz-devs/arviz/pull/1406))
* Switch to [MyST](https://myst-parser.readthedocs.io/en/latest/) and [MyST-NB](https://myst-nb.readthedocs.io/en/latest/index.html)
  for markdown/notebook parsing in docs ([1406](https://github.com/arviz-devs/arviz/pull/1406))
* Incorporated `input_core_dims` in `hdi` and `plot_hdi` docstrings ([1410](https://github.com/arviz-devs/arviz/pull/1410))
* Add documentation pages about experimental `SamplingWrapper`s usage ([1373](https://github.com/arviz-devs/arviz/pull/1373))
* Show example titles in gallery page ([1484](https://github.com/arviz-devs/arviz/pull/1484))
* Add `sample_stats` naming convention to the InferenceData schema ([1063](https://github.com/arviz-devs/arviz/pull/1063))
* Extend api documentation about `InferenceData` methods ([1338](https://github.com/arviz-devs/arviz/pull/1338))

### Experimental
* Modified `SamplingWrapper` base API ([1373](https://github.com/arviz-devs/arviz/pull/1373))

## v0.10.0 (2020 Sep 24)
### New features
* Added InferenceData dataset containing circular variables ([1265](https://github.com/arviz-devs/arviz/pull/1265))
* Added `is_circular` argument to `plot_dist` and `plot_kde` allowing for a circular histogram (Matplotlib, Bokeh) or 1D KDE plot (Matplotlib). ([1266](https://github.com/arviz-devs/arviz/pull/1266))
* Added `to_dict` method for InferenceData object ([1223](https://github.com/arviz-devs/arviz/pull/1223))
* Added `circ_var_names` argument to `plot_trace` allowing for circular traceplot (Matplotlib) ([1336](https://github.com/arviz-devs/arviz/pull/1336))
* Ridgeplot is hdi aware. By default displays truncated densities at the specified `hdi_prop` level ([1348](https://github.com/arviz-devs/arviz/pull/1348))
* Added `plot_separation` ([1359](https://github.com/arviz-devs/arviz/pull/1359))
* Extended methods from `xr.Dataset` to `InferenceData` ([1254](https://github.com/arviz-devs/arviz/pull/1254))
* Add `extend` and `add_groups` to `InferenceData` ([1300](https://github.com/arviz-devs/arviz/pull/1300) and [1386](https://github.com/arviz-devs/arviz/pull/1386))
* Added `__iter__` method (`.items`) for InferenceData ([1356](https://github.com/arviz-devs/arviz/pull/1356))
* Add support for discrete variables in `plot_bpv` ([#1379](https://github.com/arviz-devs/arviz/pull/1379))

### Maintenance and fixes
* Automatic conversion of list/tuple to numpy array in distplot ([1277](https://github.com/arviz-devs/arviz/pull/1277))
* `plot_posterior` fix overlap of hdi and rope ([1263](https://github.com/arviz-devs/arviz/pull/1263))
* `plot_dist` bins argument error fixed ([1306](https://github.com/arviz-devs/arviz/pull/1306))
* Improve handling of circular variables in `az.summary` ([1313](https://github.com/arviz-devs/arviz/pull/1313))
* Removed change of default warning in `ELPDData` string representation ([1321](https://github.com/arviz-devs/arviz/pull/1321))
* Update `radon` example dataset to current InferenceData schema specification ([1320](https://github.com/arviz-devs/arviz/pull/1320))
* Update `from_cmdstan` functionality and add warmup groups ([1330](https://github.com/arviz-devs/arviz/pull/1330) and [1351](https://github.com/arviz-devs/arviz/pull/1351))
* Restructure plotting code to be compatible with mpl>=3.3 ([1312](https://github.com/arviz-devs/arviz/pull/1312) and [1352](https://github.com/arviz-devs/arviz/pull/1352))
* Replaced `_fast_kde()` with `kde()` which now also supports circular variables via the argument `circular` ([1284](https://github.com/arviz-devs/arviz/pull/1284)).
* Increased `from_pystan` attrs information content ([1353](https://github.com/arviz-devs/arviz/pull/1353))
* Allow `plot_trace` to return and accept axes ([1361](https://github.com/arviz-devs/arviz/pull/1361))
* Update diagnostics to be on par with posterior package ([1366](https://github.com/arviz-devs/arviz/pull/1366))
* Use method="average" in `scipy.stats.rankdata` ([1380](https://github.com/arviz-devs/arviz/pull/1380))
* Add more `plot_parallel` examples ([1380](https://github.com/arviz-devs/arviz/pull/1380))
* Bump minimum xarray version to 0.16.1 ([1389](https://github.com/arviz-devs/arviz/pull/1389)
* Fix multi rope for `plot_forest` ([1390](https://github.com/arviz-devs/arviz/pull/1390))
* Bump minimum xarray version to 0.16.1 ([1389](https://github.com/arviz-devs/arviz/pull/1389))
* `from_dict` will now store warmup groups even with the main group missing ([1386](https://github.com/arviz-devs/arviz/pull/1386))
* increase robustness for repr_html handling ([1392](https://github.com/arviz-devs/arviz/pull/1392))

## v0.9.0 (2020 June 23)
### New features
* loo-pit plot. The kde is computed over the data interval (this could be shorter than [0, 1]). The HDI is computed analytically ([1215](https://github.com/arviz-devs/arviz/pull/1215))
* Added `html_repr` of InferenceData objects for jupyter notebooks. ([1217](https://github.com/arviz-devs/arviz/pull/1217))
* Added support for PyJAGS via the function `from_pyjags`. ([1219](https://github.com/arviz-devs/arviz/pull/1219) and [1245](https://github.com/arviz-devs/arviz/pull/1245))
* `from_pymc3` can now retrieve `coords` and `dims` from model context ([1228](https://github.com/arviz-devs/arviz/pull/1228), [1240](https://github.com/arviz-devs/arviz/pull/1240) and [1249](https://github.com/arviz-devs/arviz/pull/1249))
* `plot_trace` now supports multiple aesthetics to identify chain and variable
  shape and support matplotlib aliases ([1253](https://github.com/arviz-devs/arviz/pull/1253))
* `plot_hdi` can now take already computed HDI values ([1241](https://github.com/arviz-devs/arviz/pull/1241))
* `plot_bpv`. A new plot for Bayesian p-values ([1222](https://github.com/arviz-devs/arviz/pull/1222))

### Maintenance and fixes
* Include data from `MultiObservedRV` to `observed_data` when using
  `from_pymc3` ([1098](https://github.com/arviz-devs/arviz/pull/1098))
* Added a note on `plot_pair` when trying to use `plot_kde` on `InferenceData`
  objects. ([1218](https://github.com/arviz-devs/arviz/pull/1218))
* Added `log_likelihood` argument to `from_pyro` and a warning if log likelihood cannot be obtained ([1227](https://github.com/arviz-devs/arviz/pull/1227))
* Skip tests on matplotlib animations if ffmpeg is not installed ([1227](https://github.com/arviz-devs/arviz/pull/1227))
* Fix hpd bug where arguments were being ignored ([1236](https://github.com/arviz-devs/arviz/pull/1236))
* Remove false positive warning in `plot_hdi` and fixed matplotlib axes generation ([1241](https://github.com/arviz-devs/arviz/pull/1241))
* Change the default `zorder` of scatter points from `0` to `0.6` in `plot_pair` ([1246](https://github.com/arviz-devs/arviz/pull/1246))
* Update `get_bins` for numpy 1.19 compatibility ([1256](https://github.com/arviz-devs/arviz/pull/1256))
* Fixes to `rug`, `divergences` arguments in `plot_trace` ([1253](https://github.com/arviz-devs/arviz/pull/1253))

### Deprecation
* Using `from_pymc3` without a model context available now raises a
  `FutureWarning` and will be deprecated in a future version ([1227](https://github.com/arviz-devs/arviz/pull/1227))
* In `plot_trace`, `chain_prop` and `compact_prop` as tuples will now raise a
  `FutureWarning` ([1253](https://github.com/arviz-devs/arviz/pull/1253))
* `hdi` with 2d data raises a FutureWarning ([1241](https://github.com/arviz-devs/arviz/pull/1241))

### Documentation
* A section has been added to the documentation at InferenceDataCookbook.ipynb illustrating the use of ArviZ in conjunction with PyJAGS. ([1219](https://github.com/arviz-devs/arviz/pull/1219) and [1245](https://github.com/arviz-devs/arviz/pull/1245))
* Fixed inconsistent capitalization in `plot_hdi` docstring ([1221](https://github.com/arviz-devs/arviz/pull/1221))
* Fixed and extended `InferenceData.map` docs ([1255](https://github.com/arviz-devs/arviz/pull/1255))

## v0.8.3 (2020 May 28)
### Maintenance and fixes
* Restructured internals of `from_pymc3` to handle old pymc3 releases and
  sliced traces and to provide useful warnings ([1211](https://github.com/arviz-devs/arviz/pull/1211))


## v0.8.2 (2020 May 25)
### Maintenance and fixes
* Fixed bug in `from_pymc3` for sliced `pymc3.MultiTrace` input ([1209](https://github.com/arviz-devs/arviz/pull/1209))

## v0.8.1 (2020 May 24)

### Maintenance and fixes
* Fixed bug in `from_pymc3` when used with PyMC3<3.9 ([1203](https://github.com/arviz-devs/arviz/pull/1203))
* Fixed enforcement of rcParam `plot.max_subplots` in `plot_trace` and
  `plot_pair` ([1205](https://github.com/arviz-devs/arviz/pull/1205))
* Removed extra subplot row and column in in `plot_pair` with `marginal=True` ([1205](https://github.com/arviz-devs/arviz/pull/1205))
* Added latest PyMC3 release to CI in addition to using GitHub default branch ([1207](https://github.com/arviz-devs/arviz/pull/1207))

### Documentation
* Use `dev` as version indicator in online documentation ([1204](https://github.com/arviz-devs/arviz/pull/1204))

## v0.8.0 (2020 May 23)

### New features
* Stats and plotting functions that provide `var_names` arg can now filter parameters based on partial naming (`filter="like"`) or regular expressions (`filter="regex"`) (see [1154](https://github.com/arviz-devs/arviz/pull/1154)).
* Add `true_values` argument for `plot_pair`. It allows for a scatter plot showing the true values of the variables ([1140](https://github.com/arviz-devs/arviz/pull/1140))
* Allow xarray.Dataarray input for plots.([1120](https://github.com/arviz-devs/arviz/pull/1120))
* Revamped the `hpd` function to make it work with mutidimensional arrays, InferenceData and xarray objects ([1117](https://github.com/arviz-devs/arviz/pull/1117))
* Skip test for optional/extra dependencies when not installed ([1113](https://github.com/arviz-devs/arviz/pull/1113))
* Add option to display rank plots instead of trace ([1134](https://github.com/arviz-devs/arviz/pull/1134))
* Add out-of-sample groups (`predictions` and `predictions_constant_data`) to `from_dict` ([1125](https://github.com/arviz-devs/arviz/pull/1125))
* Add out-of-sample groups (`predictions` and `predictions_constant_data`) and `constant_data` group to pyro and numpyro translation ([1090](https://github.com/arviz-devs/arviz/pull/1090), [1125](https://github.com/arviz-devs/arviz/pull/1125))
* Add `num_chains` and `pred_dims` arguments to from_pyro and from_numpyro ([1090](https://github.com/arviz-devs/arviz/pull/1090), [1125](https://github.com/arviz-devs/arviz/pull/1125))
* Integrate jointplot into pairplot, add point-estimate and overlay of plot kinds ([1079](https://github.com/arviz-devs/arviz/pull/1079))
* New grayscale style. This also add two new cmaps `cet_grey_r` and `cet_grey_r`. These are perceptually uniform gray scale cmaps from colorcet (linear_grey_10_95_c0) ([1164](https://github.com/arviz-devs/arviz/pull/1164))
* Add warmup groups to InferenceData objects, initial support for PyStan ([1126](https://github.com/arviz-devs/arviz/pull/1126)) and PyMC3 ([1171](https://github.com/arviz-devs/arviz/pull/1171))
* `hdi_prob` will not plot hdi if argument `hide` is passed. Previously `credible_interval` would omit HPD if `None` was passed  ([1176](https://github.com/arviz-devs/arviz/pull/1176))
* Add `stats.ic_pointwise` rcParam ([1173](https://github.com/arviz-devs/arviz/pull/1173))
* Add `var_name` argument to information criterion calculation: `compare`,
  `loo` and `waic` ([1173](https://github.com/arviz-devs/arviz/pull/1173))

### Maintenance and fixes
* Fixed `plot_pair` functionality for two variables with bokeh backend ([1179](https://github.com/arviz-devs/arviz/pull/1179))
* Changed `diagonal` argument for `marginals` and fixed `point_estimate_marker_kwargs` in `plot_pair` ([1167](https://github.com/arviz-devs/arviz/pull/1167))
* Fixed behaviour of `credible_interval=None` in `plot_posterior` ([1115](https://github.com/arviz-devs/arviz/pull/1115))
* Fixed hist kind of `plot_dist` with multidimensional input ([1115](https://github.com/arviz-devs/arviz/pull/1115))
* Fixed `TypeError` in `transform` argument of `plot_density` and `plot_forest` when `InferenceData` is a list or tuple ([1121](https://github.com/arviz-devs/arviz/pull/1121))
* Fixed overlaid pairplots issue ([1135](https://github.com/arviz-devs/arviz/pull/1135))
* Update Docker building steps ([1127](https://github.com/arviz-devs/arviz/pull/1127))
* Updated benchmarks and moved to asv_benchmarks/benchmarks ([1142](https://github.com/arviz-devs/arviz/pull/1142))
* Moved `_fast_kde`, `_fast_kde_2d`, `get_bins` and `_sturges_formula` to `numeric_utils` and `get_coords` to `utils` ([1142](https://github.com/arviz-devs/arviz/pull/1142))
* Rank plot: rename `axes` argument to `ax` ([1144](https://github.com/arviz-devs/arviz/pull/1144))
* Added a warning specifying log scale is now the default in compare/loo/waic functions ([1150](https://github.com/arviz-devs/arviz/pull/1150)).
* Fixed bug in `plot_posterior` with rcParam "plot.matplotlib.show" = True ([1151](https://github.com/arviz-devs/arviz/pull/1151))
* Set `fill_last` argument of `plot_kde` to False by default ([1158](https://github.com/arviz-devs/arviz/pull/1158))
* plot_ppc animation: improve docs and error handling ([1162](https://github.com/arviz-devs/arviz/pull/1162))
* Fix import error when wrapped function docstring is empty ([1192](https://github.com/arviz-devs/arviz/pull/1192))
* Fix passing axes to plot_density with several datasets ([1198](https://github.com/arviz-devs/arviz/pull/1198))

### Deprecation
* `hpd` function deprecated in favor of `hdi`. `credible_interval` argument replaced by `hdi_prob`throughout with exception of `plot_loo_pit` ([1176](https://github.com/arviz-devs/arviz/pull/1176))
* `plot_hpd` function deprecated in favor of `plot_hdi`. ([1190](https://github.com/arviz-devs/arviz/pull/1190))

### Documentation
* Add classifier to `setup.py` including Matplotlib framework ([1133](https://github.com/arviz-devs/arviz/pull/1133))
* Image thumbs generation updated to be Bokeh 2 compatible ([1116](https://github.com/arviz-devs/arviz/pull/1116))
* Add new examples for `plot_pair` ([1110](https://github.com/arviz-devs/arviz/pull/1110))
* Add examples for `psislw` and `r2_score` ([1129](https://github.com/arviz-devs/arviz/pull/1129))
* Add more examples on 2D kde customization ([1158](https://github.com/arviz-devs/arviz/pull/1158))
* Make docs compatible with sphinx3 and configure `intersphinx` for better
  references ([1184](https://github.com/arviz-devs/arviz/pull/1184))
* Extend the developer guide and add it to the website ([1184](https://github.com/arviz-devs/arviz/pull/1184))

## v0.7.0 (2020 Mar 2)

### New features
* Add out-of-sample predictions (`predictions` and  `predictions_constant_data` groups) to pymc3, pystan, cmdstan and cmdstanpy translations ([983](https://github.com/arviz-devs/arviz/pull/983), [1032](https://github.com/arviz-devs/arviz/pull/1032) and [1064](https://github.com/arviz-devs/arviz/pull/1064))
* Started adding pointwise log likelihood storage support ([794](https://github.com/arviz-devs/arviz/pull/794), [1044](https://github.com/arviz-devs/arviz/pull/1044) and [1064](https://github.com/arviz-devs/arviz/pull/1064))
* Add out-of-sample predictions (`predictions` and  `predictions_constant_data` groups) to pymc3 and pystan translations ([983](https://github.com/arviz-devs/arviz/pull/983) and [1032](https://github.com/arviz-devs/arviz/pull/1032))
* Started adding pointwise log likelihood storage support ([794](https://github.com/arviz-devs/arviz/pull/794), [1044](https://github.com/arviz-devs/arviz/pull/1044))
* Violinplot: rug-plot option ([997](https://github.com/arviz-devs/arviz/pull/997))
* Integrated rcParams `plot.point_estimate` ([994](https://github.com/arviz-devs/arviz/pull/994)), `stats.ic_scale` ([993](https://github.com/arviz-devs/arviz/pull/993)) and `stats.credible_interval` ([1017](https://github.com/arviz-devs/arviz/pull/1017))
* Added `group` argument to `plot_ppc` ([1008](https://github.com/arviz-devs/arviz/pull/1008)), `plot_pair` ([1009](https://github.com/arviz-devs/arviz/pull/1009)) and `plot_joint` ([1012](https://github.com/arviz-devs/arviz/pull/1012))
* Added `transform` argument to `plot_trace`, `plot_forest`, `plot_pair`, `plot_posterior`, `plot_rank`, `plot_parallel`,  `plot_violin`,`plot_density`, `plot_joint` ([1036](https://github.com/arviz-devs/arviz/pull/1036))
* Add `skipna` argument to `hpd` and `summary` ([1035](https://github.com/arviz-devs/arviz/pull/1035))
* Added `transform` argument to `plot_trace`, `plot_forest`, `plot_pair`, `plot_posterior`, `plot_rank`, `plot_parallel`,  `plot_violin`,`plot_density`, `plot_joint` ([1036](https://github.com/arviz-devs/arviz/pull/1036))
* Add `marker` functionality to `bokeh_plot_elpd` ([1040](https://github.com/arviz-devs/arviz/pull/1040))
* Add `ridgeplot_quantiles` argument to `plot_forest` ([1047](https://github.com/arviz-devs/arviz/pull/1047))
* Added the functionality [interactive legends](https://docs.bokeh.org/en/1.4.0/docs/user_guide/interaction/legends.html) for bokeh plots of `densityplot`, `energyplot`
  and `essplot` ([1024](https://github.com/arviz-devs/arviz/pull/1024))
* New defaults for cross validation: `loo` (old: waic) and `log` -scale (old: `deviance` -scale) ([1067](https://github.com/arviz-devs/arviz/pull/1067))
* **Experimental Feature**: Added `arviz.wrappers` module to allow ArviZ to refit the models if necessary ([771](https://github.com/arviz-devs/arviz/pull/771))
* **Experimental Feature**: Added `reloo` function to ArviZ ([771](https://github.com/arviz-devs/arviz/pull/771))
* Added new helper function `matplotlib_kwarg_dealiaser` ([1073](https://github.com/arviz-devs/arviz/pull/1073))
* ArviZ version to InferenceData attributes. ([1086](https://github.com/arviz-devs/arviz/pull/1086))
* Add `log_likelihood` argument to `from_pymc3` ([1082](https://github.com/arviz-devs/arviz/pull/1082))
* Integrated rcParams for `plot.bokeh.layout` and `plot.backend`. ([1089](https://github.com/arviz-devs/arviz/pull/1089))
* Add automatic legends in `plot_trace` with compact=True (matplotlib only) ([1070](https://github.com/arviz-devs/arviz/pull/1070))
* Updated hover information for `plot_pair` with bokeh backend ([1074](https://github.com/arviz-devs/arviz/pull/1074))


### Maintenance and fixes
* Fixed bug in density and posterior plot bin computation ([1049](https://github.com/arviz-devs/arviz/pull/1049))
* Fixed bug in density plot ax argument ([1049](https://github.com/arviz-devs/arviz/pull/1049))
* Fixed bug in extracting prior samples for cmdstanpy ([979](https://github.com/arviz-devs/arviz/pull/979))
* Fix erroneous warning in traceplot ([989](https://github.com/arviz-devs/arviz/pull/989))
* Correct bfmi denominator ([991](https://github.com/arviz-devs/arviz/pull/991))
* Removed parallel from jit full ([996](https://github.com/arviz-devs/arviz/pull/996))
* Rename flat_inference_data_to_dict ([1003](https://github.com/arviz-devs/arviz/pull/1003))
* Violinplot: fix histogram ([997](https://github.com/arviz-devs/arviz/pull/997))
* Convert all instances of SyntaxWarning to UserWarning ([1016](https://github.com/arviz-devs/arviz/pull/1016))
* Fix `point_estimate` in `plot_posterior` ([1038](https://github.com/arviz-devs/arviz/pull/1038))
* Fix interpolation `hpd_plot` ([1039](https://github.com/arviz-devs/arviz/pull/1039))
* Fix `io_pymc3.py` to handle models with `potentials` ([1043](https://github.com/arviz-devs/arviz/pull/1043))
* Fix several inconsistencies between schema and `from_pymc3` implementation
  in groups `prior`, `prior_predictive` and `observed_data` ([1045](https://github.com/arviz-devs/arviz/pull/1045))
* Stabilize covariance matrix for `plot_kde_2d` ([1075](https://github.com/arviz-devs/arviz/pull/1075))
* Removed extra dim in `prior` data in `from_pyro` ([1071](https://github.com/arviz-devs/arviz/pull/1071))
* Moved CI and docs (build & deploy) to Azure Pipelines and started using codecov ([1080](https://github.com/arviz-devs/arviz/pull/1080))
* Fixed bug in densityplot when variables differ between models ([1096](https://github.com/arviz-devs/arviz/pull/1096))

### Deprecation
* `from_pymc3` now requires PyMC3>=3.8

### Documentation
* Updated `InferenceData` schema specification (`log_likelihood`,
  `predictions` and `predictions_constant_data` groups)
* Clarify the usage of `plot_joint` ([1001](https://github.com/arviz-devs/arviz/pull/1001))
* Added the API link of function to examples ([1013](https://github.com/arviz-devs/arviz/pull/1013))
* Updated PyStan_schema_example to include example of out-of-sample prediction ([1032](https://github.com/arviz-devs/arviz/pull/1032))
* Added example for `concat` method ([1037](https://github.com/arviz-devs/arviz/pull/1037))


## v0.6.1 (2019 Dec 28)

### New features
* Update for pair_plot divergences can be selected
* Default tools follow global (ArviZ) defaults
* Add interactive legend for a plot, if only two variables are used in pairplot

### Maintenance and fixes
* Change `packaging` import from absolute to relative format, explicitly importing `version` function


## v0.6.0 (2019 Dec 24)

### New features

* Initial bokeh support.
* ArviZ.jl a Julia interface to ArviZ (@sethaxen )

### Maintenance and fixes

* Fully support `numpyro` (@fehiepsi )
* log_likelihood and observed data from pyro
* improve rcparams
* fix `az.concat` functionality (@anzelpwj )

### Documentation
* distplot docstring plotting example (@jscarbor )

## v0.5.1 (2019 Sep 16)

### Maintenance and fixes
* Comment dev requirements in setup.py


## v0.5.0 (2019 Sep 15)

## New features
* Add from_numpyro Integration ([811](https://github.com/arviz-devs/arviz/pull/811))
* Numba Google Summer of Code additions (https://ban-zee.github.io/jekyll/update/2019/08/19/Submission.html)
* Model checking, Inference Data, and Convergence assessments (https://github.com/OriolAbril/gsoc2019/blob/master/final_work_submission.md)


## v0.4.1 (2019 Jun 9)

### New features
* Reorder stats columns ([695](https://github.com/arviz-devs/arviz/pull/695))
* Plot Forest reports ess and rhat by default([685](https://github.com/arviz-devs/arviz/pull/685))
* Add pointwise elpd ([678](https://github.com/arviz-devs/arviz/pull/678))

### Maintenance and fixes
* Fix io_pymc3 bug ([693](https://github.com/arviz-devs/arviz/pull/693))
* Fix io_pymc3 warning ([686](https://github.com/arviz-devs/arviz/pull/686))
* Fix 0 size bug with pystan ([677](https://github.com/arviz-devs/arviz/pull/677))


## v0.4.0 (2019 May 20)

### New features
* Add plot_dist ([592](https://github.com/arviz-devs/arviz/pull/592))
* New rhat and ess ([623](https://github.com/arviz-devs/arviz/pull/623))
* Add plot_hpd ([611](https://github.com/arviz-devs/arviz/pull/611))
* Add plot_rank ([625](https://github.com/arviz-devs/arviz/pull/625))

### Deprecation
* Remove load_data and save_data ([625](https://github.com/arviz-devs/arviz/pull/625))


## v0.3.3 (2019 Feb 23)

### New features
* Plot ppc supports multiple chains ([526](https://github.com/arviz-devs/arviz/pull/526))
* Plot titles now wrap ([441](https://github.com/arviz-devs/arviz/pull/441))
* plot_density uses a grid ([379](https://github.com/arviz-devs/arviz/pull/379))
* emcee reader support ([550](https://github.com/arviz-devs/arviz/pull/550))
* Animations in plot_ppc ([546](https://github.com/arviz-devs/arviz/pull/546))
* Optional dictionary for stat_funcs in summary ([583](https://github.com/arviz-devs/arviz/pull/583))
* Can exclude variables in selections with negated variable names ([574](https://github.com/arviz-devs/arviz/pull/574))

### Maintenance and fixes
* Order maintained with xarray_var_iter ([557](https://github.com/arviz-devs/arviz/pull/557))
* Testing very improved (multiple)
* Fix nan handling in effective sample size ([573](https://github.com/arviz-devs/arviz/pull/573))
* Fix kde scaling ([582](https://github.com/arviz-devs/arviz/pull/582))
* xticks for discrete variables ([586](https://github.com/arviz-devs/arviz/pull/586))
* Empty InferenceData saves consistent with netcdf ([577](https://github.com/arviz-devs/arviz/pull/577))
* Removes numpy pinning ([594](https://github.com/arviz-devs/arviz/pull/594))

### Documentation
* JOSS and Zenodo badges ([537](https://github.com/arviz-devs/arviz/pull/537))
* Gitter badge ([548](https://github.com/arviz-devs/arviz/pull/548))
* Docs for combining InferenceData ([590](https://github.com/arviz-devs/arviz/pull/590))

## v0.3.2 (2019 Jan 15)

### New features

* Support PyStan3 ([464](https://github.com/arviz-devs/arviz/pull/464))
* Add some more information to the inference data of tfp ([447](https://github.com/arviz-devs/arviz/pull/447))
* Use Split R-hat ([477](https://github.com/arviz-devs/arviz/pull/477))
* Normalize from_xyz functions ([490](https://github.com/arviz-devs/arviz/pull/490))
* KDE: Display quantiles ([479](https://github.com/arviz-devs/arviz/pull/479))
* Add multiple rope support to `plot_forest` ([448](https://github.com/arviz-devs/arviz/pull/448))
* Numba jit compilation to speed up some methods ([515](https://github.com/arviz-devs/arviz/pull/515))
* Add `from_dict` for easier creation of az.InferenceData objects  ([524](https://github.com/arviz-devs/arviz/pull/524))
* Add stable logsumexp ([522](https://github.com/arviz-devs/arviz/pull/522))


### Maintenance and fixes

* Fix for `from_pyro` with multiple chains ([463](https://github.com/arviz-devs/arviz/pull/463))
* Check `__version__` for attr ([466](https://github.com/arviz-devs/arviz/pull/466))
* And exception to plot compare ([461](https://github.com/arviz-devs/arviz/pull/461))
* Add Docker Testing to travisCI ([473](https://github.com/arviz-devs/arviz/pull/473))
* fix jointplot warning ([478](https://github.com/arviz-devs/arviz/pull/478))
* Fix tensorflow import bug ([489](https://github.com/arviz-devs/arviz/pull/489))
* Rename N_effective to S_effective ([505](https://github.com/arviz-devs/arviz/pull/505))


### Documentation

* Add docs to plot compare ([461](https://github.com/arviz-devs/arviz/pull/461))
* Add InferenceData tutorial in header ([502](https://github.com/arviz-devs/arviz/pull/502))
* Added figure to InferenceData tutorial ([510](https://github.com/arviz-devs/arviz/pull/510))


## v0.3.1 (2018 Dec 18)

### Maintenance and fixes
* Fix installation problem with release 0.3.0

## v0.3.0 (2018 Dec 14)

* First Beta Release
<img src="https://arviz-devs.github.io/arviz/_static/logo.png" height=100></img>

[![PyPI version](https://badge.fury.io/py/arviz.svg)](https://badge.fury.io/py/arviz)
[![Azure Build Status](https://dev.azure.com/ArviZ/ArviZ/_apis/build/status/arviz-devs.arviz?branchName=main)](https://dev.azure.com/ArviZ/ArviZ/_build/latest?definitionId=1&branchName=main)
[![codecov](https://codecov.io/gh/arviz-devs/arviz/branch/main/graph/badge.svg)](https://codecov.io/gh/arviz-devs/arviz)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/arviz-devs/community)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01143/status.svg)](https://doi.org/10.21105/joss.01143) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2540945.svg)](https://doi.org/10.5281/zenodo.2540945)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)
# ArviZ

ArviZ (pronounced "AR-_vees_") is a Python package for exploratory analysis of Bayesian models.
Includes functions for posterior analysis, data storage, model checking, comparison and diagnostics.

### ArviZ in other languages
ArviZ also has a Julia wrapper available [ArviZ.jl](https://arviz-devs.github.io/ArviZ.jl/stable/).

## Documentation

The ArviZ documentation can be found in the [official docs](https://arviz-devs.github.io/arviz/index.html).
First time users may find the [quickstart](https://arviz-devs.github.io/arviz/getting_started/Introduction.html)
to be helpful. Additional guidance can be found in the
[usage documentation](https://arviz-devs.github.io/arviz/usage.html).


## Installation

### Stable
ArviZ is available for installation from [PyPI](https://pypi.org/project/arviz/).
The latest stable version can be installed using pip:

```
pip install arviz
```

ArviZ is also available through [conda-forge](https://anaconda.org/conda-forge/arviz).

```
conda install -c conda-forge arviz
```

### Development
The latest development version can be installed from the main branch using pip:

```
pip install git+git://github.com/arviz-devs/arviz.git
```

Another option is to clone the repository and install using git and setuptools:

```
git clone https://github.com/arviz-devs/arviz.git
cd arviz
python setup.py install
```

-------------------------------------------------------------------------------
## [Gallery](https://arviz-devs.github.io/arviz/examples/index.html)

<p>
<table>
<tr>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_forest_ridge.html">
  <img alt="Ridge plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_forest_ridge_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_parallel.html">
  <img alt="Parallel plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_parallel_minmax_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_trace.html">
  <img alt="Trace plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_trace_bars_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_density.html">
  <img alt="Density plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_density_thumb.png" />
  </a>
  </td>

  </tr>
  <tr>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_posterior.html">
  <img alt="Posterior plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_posterior_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_joint.html">
  <img alt="Joint plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_joint_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_ppc.html">
  <img alt="Posterior predictive plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_ppc_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_pair.html">
  <img alt="Pair plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_pair_thumb.png" />
  </a>
  </td>

  </tr>
  <tr>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_energy.html">
  <img alt="Energy Plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_pair_hex_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_violin.html">
  <img alt="Violin Plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_violin_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_forest.html">
  <img alt="Forest Plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_forest_thumb.png" />
  </a>
  </td>

  <td>
  <a href="https://arviz-devs.github.io/arviz/examples/plot_autocorr.html">
  <img alt="Autocorrelation Plot"
  src="https://raw.githubusercontent.com/arviz-devs/arviz/gh-pages/_static/plot_autocorr_thumb.png" />
  </a>
  </td>

</tr>
</table>

## Dependencies

ArviZ is tested on Python 3.7, 3.8 and 3.9, and depends on NumPy, SciPy, xarray, and Matplotlib.


## Citation


If you use ArviZ and want to cite it please use [![DOI](http://joss.theoj.org/papers/10.21105/joss.01143/status.svg)](https://doi.org/10.21105/joss.01143)

Here is the citation in BibTeX format

```
@article{arviz_2019,
  doi = {10.21105/joss.01143},
  url = {https://doi.org/10.21105/joss.01143},
  year = {2019},
  publisher = {The Open Journal},
  volume = {4},
  number = {33},
  pages = {1143},
  author = {Ravin Kumar and Colin Carroll and Ari Hartikainen and Osvaldo Martin},
  title = {ArviZ a unified library for exploratory analysis of Bayesian models in Python},
  journal = {Journal of Open Source Software}
}
```


## Contributions
ArviZ is a community project and welcomes contributions.
Additional information can be found in the [Contributing Readme](https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md)


## Code of Conduct
ArviZ wishes to maintain a positive community. Additional details
can be found in the [Code of Conduct](https://github.com/arviz-devs/arviz/blob/main/CODE_OF_CONDUCT.md)

## Donations
ArviZ is a non-profit project under NumFOCUS umbrella. If you want to support ArviZ financially, you can donate [here](https://numfocus.org/donate-to-arviz).

## Sponsors
[![NumFOCUS](https://www.numfocus.org/wp-content/uploads/2017/07/NumFocus_LRG.png)](https://numfocus.org)
# Main Governance Document

## The Project

The ArviZ Project (The Project) is an open source software project
affiliated with the 501c3 NumFOCUS Foundation. The goal of The Project is to
develop open source software and deploy open and public websites and services
for reproducible, exploratory and interactive computing. The Software developed
by The Project is released under the Apache 2 open source license,
developed openly and hosted in public GitHub repositories under the
[GitHub organization](https://github.com/arviz-devs). Examples of
Project Software include the ArviZ code and the Documentation, etc. The Services run by the
Project consist of public websites and web-services that are hosted
at [https://arviz-devs.github.io/arviz/](https://arviz-devs.github.io/arviz/)
The Project is developed by a team of distributed developers, called
Contributors. Contributors are individuals who have contributed code,
documentation, designs or other work to one or more Project repositories.
Anyone can be a Contributor. Contributors can be affiliated with any legal
entity or none. Contributors participate in the project by submitting,
reviewing and discussing GitHub Pull Requests and Issues and participating in
open and public Project discussions on GitHub, Slack, and Gitter chat rooms.
The foundation of Project participation is openness and transparency.

There have been many Contributors to the Project, whose contributions are listed in the
logs of any of the repositories under the ArviZ-devs organization.

The Project Community consists of all Contributors and Users of the Project.
Contributors work on behalf of and are responsible to the larger Project
Community and we strive to keep the barrier between Contributors and Users as
low as possible.

The Project is formally affiliated with the 501c3 NumFOCUS Foundation
([http://numfocus.org](http://numfocus.org)),  NumFOCUS is the
only legal entity that has a formal relationship with the project.

## Governance

This section describes the governance and leadership model of The Project.

The foundations of Project governance are:

-   Openness & Transparency
-   Active Contribution
-   Institutional Neutrality

Traditionally, project leadership was unstructured but primarily driven
by a subset of Core Contributors whose active and consistent
contributions have been recognized by their receiving “commit rights” to the
Project GitHub repositories. In general all Project decisions are made through
consensus among the Core Contributors with input from the Community.

While this approach has served us well, as the Project grows and faces more
legal and financial decisions and interacts with other institutions, we see a
need for a more formal governance model. Moving forward The Project leadership
will consist of a Random Variables Council. We view this governance model as
the formalization of what we are already doing, rather than a change in
direction.


## Community Architecture

* General Contributors
* Core Contributors (of which Council members are also a part of)
* 7 Person steering council (Random Variables Council)

Anyone working with ArviZ has the responsibility to personally uphold
the Code of Conduct. Core Contributors have the additional responsibility
of _enforcing_ the Code of Conduct to maintain a safe community.

## Core Contributors
Core Contributors are those who have provided consistent and meaningful contributions to ArviZ.
These can be, but are not limited to, code contributions, community contributions, tutorial
development etc. Core Contributors will be given the ability to manage the ArviZ GitHub
repository, including code merges to main. This does not necessarily mean Core Contributors
must submit code, but more so signifies trust with the project as a whole.

### Core Contributor Responsibilities
* Enforce code of conduct
* Maintain a check against Council

### Core Contributor Nominations and Confirmation Process
Current Core Contributors can nominate candidates to become Core Contributors by
requesting so in a GitHub issue. If nominated candidates accept their nomination
(explicit comment approving nomination on the issue or "thumbs-up" emoji on the same issue),
then they can be considered by the Council: on the first of the month following a
nomination, the Council will vote on each nominee using the process above.

Voting will be private with results published on the issue ticket. In the case of a
rejection, results must include the reasons behind the decision (e.g. the time since
starting to contribute is deemed too short for now). The candidate
would then have to wait 3 months to be considered again.

### Current Core Contributors
* Oriol Abril-Pla (@OriolAbril)
* Agustina Arroyuelo (@agustinaarroyuelo)
* Alex Andorra (@AlexAndorra)
* Seth Axen (@sethaxen)
* Colin Carroll (@ColCarroll)
* Piyush Gautam (@percygautam)
* Robert P. Goldman (@rpgoldman)
* Ari Hartikainen (@ahartikainen)
* Ravin Kumar (@canyon289)
* Osvaldo Martin (@aloctavodia)
* Mitzi Morris (@mitzimorris)
* Du Phan (@fehiepsi)
* Aki Vehtari (@avehtari)

## Random Variables Council
The Project will have a Steering Council that consists of Core Contributors
who have produced contributions that are substantial in quality and quantity,
and sustained over at least one year. The overall role of the Council is to
ensure, taking input from the Community, the
long-term well-being of the project, both technically and as a community.

During the everyday project activities, council members participate in all
discussions, code review and other project activities as peers with all other
Contributors and the Community. In these everyday activities, Council Members
do not have any special power or privilege through their membership on the
Council. However, it is expected that because of the quality and quantity of
their contributions and their expert knowledge of the Project Software and
Services that Council Members will provide useful guidance, both technical and
in terms of project direction, to potentially less experienced contributors.

Council Members will have the responsibility of
* Removing members, including Council Members, if they are in violation of the Code of Conduct
* Making decisions when regular community discussion does not produce consensus on an issue
  in a reasonable time frame.
* Making decisions about strategic collaborations with other organizations or individuals.
* Making decisions about the overall scope, vision and direction of the project.
* Developing funding sources
* Deciding how to disburse funds with consultation from Core Contributors

The council may choose to delegate these responsibilities to sub-committees. If so, Council members must update this document to make the delegation clear.

Note that an individual council member does not have the power to unilaterally wield these
responsibilities, but the council as a whole must jointly make these decisions.
In other words, Council Members are first and foremost Core Contributors, but only
when needed they can collectively make decisions for the health of the project.

### Current Random Variables Council members
The current RV Council members are:

* Oriol Abril-Pla (@OriolAbril)
* Alex Andorra (@AlexAndorra)
* Seth Axen (@sethaxen)
* Colin Carroll (@ColCarroll)
* Ari Hartikainen (@ahartikainen)
* Ravin Kumar (@canyon289)
* Osvaldo Martin (@aloctavodia)

The election record can be found at [elections/ArviZ_2020.md](https://github.com/arviz-devs/arviz/blob/main/elections/ArviZ_2020.md)

### Council Decision Making Process
By and large we expect the decisions in ArviZ to be made _ad hoc_ and require little formal
coordination and with the community at large. However, for controversial proposals and
new Core Contributors the council may need to intervene to make the final decision
in a group vote.

#### Call for a vote
Core Contributors can call for a vote to resolve a target issue they feel has
been stale for too long and for which informal consensus appears unlikely. For a vote
to be called, the target issue must be at least 2 months old.

To do so, they have to open a proposal issue ticket labeled "Council Vote".
The proposal issue should contain a link to the target issue and a proposal
on how to resolve it. Proposals should include a statement making clear what it
means to "agree" or to "disagree".

Before voting starts, at least 3 days will be left for Core
Contributors to raise doubts about the proposal's _phrasing_, no extra discussion will take
place in the proposal issue. Proposal issues should be locked from creation
to prevent attracting discussion from people not familiar with the decision
process.

#### Voting process

* Each Council Member will vote either "Yes", "No", or "Neutral".
* It is recommended that all Council Members expose their reasons when voting.
  "No" votes, however, _must_ list the reasons for disagreement. Any "No" vote
  with no reason listed will be considered a "Neutral" vote.
* An absence of vote is considered as "Neutral".
* Voting will remain open for at least 3 days.
* For the proposal to pass, at least 60% of the council must vote "Yes", and no more
  than 20% can vote "No".

For decisions about the project the Council will perform it directly on the
proposal issue.
For decisions about people, such as electing or ejecting Core Contributors, the Council
will vote privately. However the decision will be posted publicly in an issue ticket.

### Private communications of the Council

Unless specifically required, all Council discussions and activities will be
between public (GitHub, Gitter), and partially public channels (Slack)
and done in collaboration and discussion with the Core Contributors
and Community. The Council will have a private channel that will be used
sparingly and only when a specific matter requires privacy. When private
communications and decisions are needed, the Council will do its best to
summarize those to the Community after eliding personal/private/sensitive
information that should not be posted to the public internet.

### Conflict of interest

It is expected that Council Members will be employed at a wide
range of companies, universities and non-profit organizations. Because of this,
it is possible that Members will have conflict of interests. Such conflict of
interests include, but are not limited to:

-   Financial interests, such as investments, employment or contracting work,
    outside of The Project that may influence their work on The Project.
-   Access to proprietary information of their employer that could potentially
    leak into their work with the Project.

All members of the Council shall disclose to the rest of the
Council any conflict of interest they may have. Members with a conflict of
interest in a particular issue may participate in Council discussions on that
issue, but must recuse themselves from voting on the issue.

### Council Selection Process

#### Eligibility
* Must be core contributor for at least one year

#### Nominations
* Nominations are taken over a public GitHub issue ticket over the course of 2 weeks.
* Only Core Contributors may nominate folks
* Self Nominations are allowed
* At the conclusion of the 2 weeks, the list of nominations is posted on the ticket and this ticket is closed.


#### Election Process

* Voting occurs over a period of at least 1 week, at the conclusion of the nominations.
  Voting is blind and mediated by either an application or a third party like NumFOCUS.
  Each voter can vote zero or more times, once per each candidate. As this is not about
  ranking but about capabilities, voters vote on a yes/neutral/no basis
  per candidate -- “would I trust this person to lead ArviZ?”.
* Candidates are evaluated independently, each candidate having 60% or more of
  yes votes _and_ less or equal than 20% of no votes is chosen. If the number
  of chosen candidates is >=4 and <=10 all candidates are confirmed and the
  election process stops here.
* In the event that either not enough or too many candidates were confirmed,
  candidates are ranked by interpreting yes=+1, neutral=0 and no=-1. If too
  many candidates were confirmed, the 10 candidates with higher rank are
  elected. If not enough candidates were chosen, the 4 candidates with higher
  rank are elected.
* In the event of a tie there will be a runoff election for the tied candidates. To avoid
  further ties and discriminate more among the tied candidates, this vote will be held
  by [Majority Judgment](https://en.wikipedia.org/wiki/Majority_judgment) (MJ): for each
  candidate, voters judge their suitability for office as either "Excellent", "Very Good",
  "Good", "Acceptable", "Poor", or "Reject". Multiple candidates may be given the same
  grade by a voter. The candidate with the highest median grade is the winner.
* If more than one candidate has the same highest median-grade, the MJ winner is
  discovered by removing (one-by-one) any grades equal in value to the shared median
  grade from each tied candidate's total. This is repeated until only one of the
  previously tied candidates is currently found to have the highest median-grade.
* If ties are still present after this second round, the winner will be chosen at random.
  Each person tied will pick an integer number in the `[1, 100]` interval and send it
  privately to the third party mediating the election. After receiving all the numbers,
  said third party will draw a random integer from random.org. The person with the
  closest circular distance, defined as `min(|a-b|, 100-|a-b|)`, will be selected.
  This process will be repeated as many times as necessary as there may be ties
  resulting from candidates choosing the same number.
* At the conclusion of voting, all the results will be posted. And at least 24
  hours will be left to challenge the election result in case there were
  suspicions of irregularities or the process had not been correctly carried
  out.

#### Length of Tenure and Reverification
* Council members term limits are 4 years, after which point their seat will come up for reelection.
* Each year on April 7th council members will be asked to restate their commitment to being on the council
* Attempts should be made to reach every council member over at least 2 communication media. For example: email, Slack, phone, or GitHub.
* If a council member does not restate their commitment their seat will be vacated.
* Inactivity can be determined by lack of substantial contribution, including votes on council, code or discussion contributions, contributions in the community or otherwise.
* In the event of a vacancy in the council, an election will be held to fill the position.
* There is no limit on the number of terms a Council Member can serve

#### Vote of No Confidence
* In exceptional circumstances, council members as well as core contributors may remove a sitting council member via a vote of no confidence. Core contributors can also call for a vote to remove the entire council -- in which case, Council Members do not vote.
* A no-confidence vote is triggered when a core team member (i.e Council member or Core contributor) calls for one publicly on an appropriate project communication channel, and two other core team members second the proposal. The initial call for a no-confidence vote must specify which type is intended -- whether it is targeting a single member or the council as a whole.
* The vote lasts for two weeks, and the people taking part in it vary:
  * If this is a single-member vote called by Core contributors, both Council members and Core contributors vote, and the vote is deemed successful if at least two thirds of voters express a lack of confidence.
  * If this is a whole-council vote, then it was necessarily called by Core contributors (since Council members can’t remove the whole Council) and only Core contributors vote. The vote is deemed successful if at least two thirds of voters express a lack of confidence.
  * If this is a single-member vote called by Council Members, only Council Members vote, and the vote is deemed successful if at least half the voters express a lack of confidence. Council Members also have the possibility to call for the whole core team to vote (i.e Council members and Core contributors), although this is not the default option. The threshold for successful vote is also at 50% of voters for this option.
* If a single-member vote succeeds, then that member is removed from the council and the resulting vacancy can be handled in the usual way.
* If a whole-council vote succeeds, the council is dissolved and a new council election is triggered immediately.

#### Ejecting Core Contributors
* Core contributors can be ejected through a simple majority vote by the council. Council members vote "Yes" or "No".
* Upon ejecting a core contributor the council must publish an issue ticket, or public document detailing the
  * Violations
  * Evidence if available
  * Remediation plan (if necessary)
  * Signatures majority of council members to validate correctness and accuracy
* Core contributors can also voluntarily leave the project by notifying the community through a public means or by notifying the entire council.

#### Voting Criteria For Future Elections
Voting for first election is restricted to establish stable governance, and to defer major decision to elected leaders
* For the first election only the people registered following the guidelines in `elections/ArviZ_2020.md` can vote
* In the first year, the council must determine voting eligibility for future elections between two criteria:
  * Core contributors
  * The contributing community at large
# ArviZ Community Code of Conduct

ArviZ adopts the NumFOCUS Code of Conduct directly. In other words, we
expect our community to treat others with kindness and understanding.


# THE SHORT VERSION  
Be kind to others. Do not insult or put down others.
Behave professionally. Remember that harassment and sexist, racist,
or exclusionary jokes are not appropriate.

All communication should be appropriate for a professional audience
including people of many different backgrounds. Sexual language and
imagery are not appropriate.

ArviZ is dedicated to providing a harassment-free community for everyone,
regardless of gender, sexual orientation, gender identity, and
expression, disability, physical appearance, body size, race,
or religion. We do not tolerate harassment of community members
in any form.

Thank you for helping make this a welcoming, friendly community for all.


# How to Submit a Report
If you feel that there has been a Code of Conduct violation an anonymous
reporting form is available.
**If you feel your safety is in jeopardy or the situation is an
emergency, we urge you to contact local law enforcement before making
a report. (In the U.S., dial 911.)**

We are committed to promptly addressing any reported issues.
If you have experienced or witnessed behavior that violates this 
Code of Conduct, please complete the form below to
make a report.

**REPORTING FORM:** https://numfocus.typeform.com/to/ynjGdT

Reports are sent to the NumFOCUS Code of Conduct Enforcement Team
(see below).

You can view the Privacy Policy and Terms of Service for TypeForm here.
The NumFOCUS Privacy Policy is here:
https://www.numfocus.org/privacy-policy


# Full Code of Conduct
The full text of the NumFOCUS/ArviZ Code of Conduct can be found on
NumFOCUS's website  
https://numfocus.org/code-of-conduct
# Contributing to ArviZ
This document outlines only the most common contributions.
Please see the [Contributing guide](https://arviz-devs.github.io/arviz/contributing/index.html)
on our documentation for a better view of how can you contribute to ArviZ.
We welcome a wide range of contributions, not only code!

## Reporting issues
If you encounter any bug or incorrect behaviour while using ArviZ,
please report an issue to our [issue tracker](https://github.com/arviz-devs/arviz/issues).
Please include any supporting information, in particular the version of
ArviZ that you are using.
The issue tracker has several templates available to help in writing the issue
and including useful supporting information.

## Contributing code
Thanks for your interest in contributing code to ArviZ!

* If this is your first time contributing to a project on GitHub, please read through our step by step guide to contributing to ArviZ
* If you have contributed to other projects on GitHub you can go straight to our [development workflow]()

### Adding new features
If you are interested in adding a new feature to ArviZ,
first submit an issue using the "Feature Request" label for the community
to discuss its place and implementation within ArviZ.
# Election 2020

## Summary of Process
This is the first election cycle for ArviZ. During this cycle we are aiming
for the following outcomes

* Ratification of GOVERNANCE.md
* Identification and listing of Core Contributors
* Nominations for Council Candidates
* Council Elections. Particularly
  * Size of Council
  * Inaugural ArviZ Council

This being the first election means that the process will have more steps than
the ones listed in GOVERNANCE.md.

The moderators for these tasks live in the PST and CST timezone so the dates are reflective of their localities.

## Process Details
### Identification of Core Contributors (June 7th)
1. A Google Form will be shared over Slack to get core contributor names, emails, and optionally Github username.
    * Emails are needed for Elections steps so NumFOCUS can mail and track votes
    * The Google Form is moderated by Ravin (@canyon289)
2. Ravin will open a Pull Request with all names and update the PR periodically until June 12th 8pm PST, which is also the cutoff time for the form.
    * On June 12th 8pm PST Ravin will make the form's responses public as well, with exception of email
3. The pull request will stay open 24 hours at a minimum. If there is no dissent the pull request will be merged. If there is discussion the timeline below will be shifted.

Note: This is a one time process for identifying core contributors. After this cycle the core contributor acceptance process will performed by the ArviZ Council once elected

### Nominations for Council (June 14th)
1. A github issue labeled *Arviz 2020 Council Elections* will be opened for nominations
2. Core contributors can nominate themselves or other core contributors.
3. Core Contributor must accept acknowledgement with thumbs up on nomination
   * Without an acknowledgement, they will not be a candidate in the election.
   * A nomination can be explicitly declined as well in a comment, if that individual
   feels they will not be able to serve as a Council Member for any reason.
4. Nominations end at June 20th anywhere on Earth. On June 21st nominations will be compiled in the last comment and the issue will be closed.

#### Historical Record
Nominations were completed on June 21st (PST). Conversation is tracked in issue ticket.
https://github.com/arviz-devs/arviz/issues/1243


### First Round Election for Council and Council size (Week of June 21)
The first election will be for council and council size. Elections will take place over
a Google Form owned by NumFOCUS. The tallying will be performed in accordance with `GOVERNANCE.md`

The election period will be exactly 7 days from when NumFOCUS email timestamp for the Google Form.

1. NumFOCUS will email all registered voters a ballot for a vote on
    * Council Size
    * Council Members
2. At conclusion of election period NumFOCUS will remove the voters identifying information and share the results (most likely in Google Sheets)
    * A table of results will be added to this document
    
#### Historical Record
Call for election took place in this issue ticket. https://github.com/arviz-devs/arviz/issues/1260  
NumFOCUS mediated election using [this google form](https://forms.gle/P9P9Yx9tMQ5jbY1L9) as template, and shared the results with us over email. Result verification proof is attached in issue ticket linked above.

The outcome of the election was a vote for a 7 Person Council, comprised of the following core contributors:
* @OriolAbril
* @aloctavodia
* @ColCarroll
* @ahartikainen
* @sethaxen
* @AlexAndorra
* @canyon289

### Tie Votes (1 week)
In the event of a tie there will be a run-off election. The specific timing will depend on the results of the election above.
The tallying will be performed in accordance with `GOVERNANCE.md`

1. NumFOCUS will email all registered voters a ballot for a vote on tied candidates
2. At conclusion of election period NumFOCUS will remove the voters identifying information and share the results (most likely in Google Sheets)
    * A table of results will be added to this document
    
#### Historical Record
A tied election was not needed. The primary election was decisive.
    
    
 
# Elections
Repository for ArviZ elections
---
title: ArviZ a unified library for exploratory analysis of Bayesian models in Python
tags:
- Bayesian statistics
- Visualization
- Probabilistic programming
- Python
authors:
- name: Ravin Kumar
  orcid: 0000-0003-0501-6098
  affiliation: 4
- name: Colin Carroll
  orcid: 0000-0001-6977-0861
  affiliation: 1
- name: Ari Hartikainen
  orcid: 0000-0002-4569-569X
  affiliation: 2
- name: Osvaldo Martin
  orcid: 0000-0001-7419-8978
  affiliation: 3
affiliations:
- name: Freebird Inc., United States
  index: 1
- name: Aalto University, Department of Civil Engineering, Espoo, Finland
  index: 2
- name: Instituto de Matemática Aplicada San Luis, UNSL-CONICET. Ejército de los Andes 950, 5700 San Luis, Argentina
  index: 3
- name: Carbon IT LLC, United States
  index: 4
date: 21 November 2018
bibliography: paper.bib
---

# Summary

While conceptually simple, Bayesian methods can be mathematically and numerically challenging. Probabilistic programming languages (PPLs) implement functions to easily build Bayesian models together with efficient automatic inference methods. This helps separate the model building from the inference, allowing practitioners to focus on their specific problems and leaving PPLs to handle the computational details for them [@zotero-null-295;@bessiere_bayesian_2013;@Ghahramani2015]. The inference process generates a *posterior distribution* — which has a central role in Bayesian statistics — together with other distributions like the *posterior predictive distribution* and the *prior predictive distribution*. The correct visualization, analysis, and interpretation of these distributions is key to properly answer the questions that motivate the inference process.

When working with Bayesian models there are a series of related tasks that need to be addressed besides inference itself:


- Diagnoses of the quality of the inference
- Model criticism, including evaluations of both model assumptions and model predictions
- Comparison of models, including model selection or model averaging
- Preparation of the results for a particular audience

Successfully performing such tasks are central to the iterative and interactive modeling process.
These tasks require both numerical and visual summaries to help statisticians or practitioners
analyze visual summaries. In the words of Persi Diaconis [@DiaconisTheoriesDataAnalysis2011]
"Exploratory data analysis seeks to reveal structure, or simple descriptions in data. We look at
numbers or graphs and try to find patterns. We pursue leads suggested by background information,
imagination, patterns perceived, and experience with other data analyses".

For these reasons we introduce ArviZ, a Python package for exploratory analysis of Bayesian models.
ArviZ aims to be a package that integrates seamlessly with established probabilistic programming 
languages like PyStan [@StanProbabilisticProgramming],
PyMC [@SalvatierProbabilisticProgrammingPython2016],
Edward [@TranDeepProbabilisticProgramming2017; @TranEdwardlibraryprobabilistic2016],
emcee [@emcee], Pyro [@bingham2018pyro], and easily integrated with novel or bespoke Bayesian
analyses.  Where the aim of the probabilistic programming languages is to make it easy to build and
solve Bayesian models, the aim of the ArviZ library is to make it easy to process and analyze the
results from the Bayesian models. We hope ArviZ will become a key Python tool for Bayesian data
analysis by allowing users to focus on problems from their domain knowledge and not on computational details.

Bayesian inference produces naturally high dimensional data. By storing each type of data resulting
from PPLs as an xarray [@HoyerxarrayNDlabeled2017] dataset, ArviZ provides labeled querying of the
data, efficient algorithms, and persistent metadata. These datasets are stored together on disk and
in code using netCDF4 [@56302;@brown_1993] groups, which are themselves built with HDF5, and allows
for well supported serialization. This functionality is implemented in the InferenceData class (see Figure 1).
In addition to the listed benefits of using netCDF and xarray, by using a single data structure all
statistical and visualization functions need to be implemented only once.

![Relationship between netCDF, az.InferenceData, and xarray](https://d2mxuefqeaa7sj.cloudfront.net/s_26E7E0D1516EA1B427269A258102C3AC9090025345CBB4CA6C7DBDA445D6595F_1542830805296_inference_data.png)


In addition to common plots for Bayesian analysis including a trace plot and forest plot,
ArviZ implements other visualizations such as a plot for posterior predictive checks, a pair plot,
and a parallel coordinate plot [@GabryVisualizationBayesianworkflow2017]. Additionally, it supports
a number of statistical checks, such as calculating the effective sample size, the r-hat statistic,
Pareto-smoothed importance sampling leave-one-out cross validation (PSIS-LOO-CV)
[@VehtariPracticalBayesianModel2015], and widely applicable information criterion (WAIC)
[@watanabe_widely_2013].  

# Funding

Work by Osvaldo Martin was supported by CONICET-Argentina and ANPCyT-Argentina (PICT-0218).

# Acknowledgments

We thank the PyMC3 Community — especially Adrian Seyboldt, Junpeng Lao, and Thomas Wiecki — as well
as the Stan community — especially Allen Riddell . We also would like to extend thanks to all the ArviZ contributors, and the contributors of the libraries used to build ArviZ
— particularly xarray, matplotlib, pandas, and numpy.

# Example Plots  
Examples of ArviZ's plotting functionality are shown in Figure 2 through Figure 5. 
Additional examples can be found in the ArviZ documentation.

![Bivariate Hexbin Plot with marginal distributions](plot_joint.png){width=80%}

![2D Kernel Density Estimation](plot_kde_2d.png){width=80%}

![Markov Chain Monte Carlo Trace](plot_trace.png){width=80%}  

![John Kruschke styled Posterior Distribution](plot_posterior.png){width=80%} 

# References

## Description
<!--
Thank you so much for your PR! To help us review your contribution, please consider the following points:

- The PR title should summarize the changes, for example "Add new group argument for the pair plot".
  Avoid non-descriptive titles such as "Addresses issue #348".

- The description should provide at least 1-2 sentences describing the pull request in detail and
  link to any relevant issues.

- Please prefix the title of incomplete contributions with [WIP] (to indicate a work in
  progress). WIPs may be useful to (1) indicate you are working on something to avoid
  duplicated work, (2) request broad review of functionality or API, or (3) seek collaborators. -->

## Checklist
<!-- Feel free to remove check-list items aren't relevant to your change -->

- [ ] Follows [official](https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md#pull-request-checklist) PR format
- [ ] Includes a sample plot to visually illustrate the changes (only for plot-related functions)
- [ ] New features are properly documented (with an example if appropriate)?
- [ ] Includes new or updated tests to cover the new feature
- [ ] Code style  correct (follows pylint and black guidelines)
- [ ] Changes are listed in [changelog](https://github.com/arviz-devs/arviz/blob/main/CHANGELOG.md#v0xx-unreleased)

<!--
Also, please consider reading the contributing guidelines and code of conduct carefully before submitting the PR. They are available at
- https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md
- https://github.com/arviz-devs/arviz/blob/main/CODE_OF_CONDUCT.md

We understand that PRs can sometimes be overwhelming, especially as the
reviews start coming in. Please let us know if the reviews are unclear or
the recommended next step seems overly demanding, if you would like help in
addressing a reviewer's comments, or if you have been waiting too long to hear
back on your PR.
-->
## Description
<!--
Thank you so much for your PR!  To help us review your contribution, please
consider the following points:

- The PR title should summarize the changes, for example "Add new group argument for the
  pair plot".  Avoid non-descriptive titles such as "Addresses issue #348". If your pull
  request addresses an issue, please use the pull request title to describe
  the issue and mention the issue number in the pull request description.

- The description should provide at least 1-2 sentences describing the pull request
  in detail (Why is this change required?  What problem does it solve?) and
  link to any relevant issues. If modifying a plot, render your plot to inspect for changes
  and copy image in the pull request message on Github

- Please prefix the title of incomplete contributions with [WIP] (to indicate a work in
  progress). WIPs may be useful to (1) indicate you are working on something to avoid
  duplicated work, (2) request broad review of functionality or API, or (3) seek collaborators.
-->

## Checklist
<!-- Feel free to remove check-list items aren't relevant to your change -->

- [ ] Does the PR follow [official](https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md#pull-request-checklist)
      PR format?
- [ ] Is the documentation [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) compliant?
- [ ] Is the fix listed in the [Documentation](https://github.com/arviz-devs/arviz/blob/main/CHANGELOG.md#documentation)
      section of the changelog?

<!--
Also, please consider reading the contributing guidelines and code of conduct carefully before submitting the PR. They are available at
- https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md
- https://github.com/arviz-devs/arviz/blob/main/CODE_OF_CONDUCT.md

- If you are contributing fixes to docstrings, please pay attention to
  https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md#docstring-formatting.
  In particular, note the difference between using single backquotes, double backquotes, and
  asterisks in the markup.

We understand that PRs can sometimes be overwhelming, especially as the
reviews start coming in.  Please let us know if the reviews are unclear or
the recommended next step seems overly demanding, if you would like help in
addressing a reviewer's comments, or if you have been waiting too long to hear
back on your PR.
-->
## Description
<!--
Thank you so much for your PR!  To help us review your contribution, please
consider the following points:

- The PR title should summarize the changes, for example "Add new group argument for the
  pair plot".  Avoid non-descriptive titles such as "Addresses issue #348". If your pull
  request addresses an issue, please use the pull request title to describe
  the issue and mention the issue number in the pull request description.

- The description should provide at least 1-2 sentences describing the pull request
  in detail (Why is this change required?  What problem does it solve?) and
  link to any relevant issues. If modifying a plot, render your plot to inspect for changes
  and copy image in the pull request message on Github

- Please prefix the title of incomplete contributions with [WIP] (to indicate a work in
  progress). WIPs may be useful to (1) indicate you are working on something to avoid
  duplicated work, (2) request broad review of functionality or API, or (3) seek collaborators.
-->


## Checklist
<!-- Feel free to remove check-list items aren't relevant to your change -->

- [ ] Does the PR follow [official](https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md#pull-request-checklist)
      PR format?
- [ ] Has included a sample plot to visually illustrate the changes? (only for plot-related functions)
- [ ] Is the new feature properly documented with an example?
- [ ] Does the PR include new or updated tests to cover the new feature (using [pytest fixture pattern](
      https://docs.pytest.org/en/latest/fixture.html#fixture))?
- [ ] Is the code style correct (follows pylint and black guidelines)?
- [ ] Is the new feature listed in the [New features](https://github.com/arviz-devs/arviz/blob/main/CHANGELOG.md#new-features)
      section of the changelog?

<!--
Also, please consider reading the contributing guidelines and code of conduct carefully before submitting the PR. They are available at
- https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md
- https://github.com/arviz-devs/arviz/blob/main/CODE_OF_CONDUCT.md

- If you are contributing fixes to docstrings, please pay attention to
  https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md#docstring-formatting.
  In particular, note the difference between using single backquotes, double backquotes, and
  asterisks in the markup.

We understand that PRs can sometimes be overwhelming, especially as the
reviews start coming in.  Please let us know if the reviews are unclear or
the recommended next step seems overly demanding, if you would like help in
addressing a reviewer's comments, or if you have been waiting too long to hear
back on your PR.
-->
## Description
<!--
Thank you so much for your PR!  To help us review your contribution, please
consider the following points:

- The PR title should summarize the changes, for example "Add new group argument for the
  pair plot".  Avoid non-descriptive titles such as "Addresses issue #348". If your pull
  request addresses an issue, please use the pull request title to describe
  the issue and mention the issue number in the pull request description.

- The description should provide at least 1-2 sentences describing the pull request
  in detail (Why is this change required?  What problem does it solve?) and
  link to any relevant issues. If modifying a plot, render your plot to inspect for changes
  and copy image in the pull request message on Github

- Please prefix the title of incomplete contributions with [WIP] (to indicate a work in
  progress). WIPs may be useful to (1) indicate you are working on something to avoid
  duplicated work, (2) request broad review of functionality or API, or (3) seek collaborators.
-->


## Checklist
<!-- Feel free to remove check-list items aren't relevant to your change -->

- [ ] Does the PR follow [official](https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md#pull-request-checklist)
      PR format?
- [ ] Does the PR include new or updated tests to prevent issue recurrence (using [pytest fixture pattern](
      https://docs.pytest.org/en/latest/fixture.html#fixture))?
- [ ] Is the code style correct (follows pylint and black guidelines)?
- [ ] Is the fix listed in the [Maintenance and fixes](https://github.com/arviz-devs/arviz/blob/main/CHANGELOG.md#maintenance-and-fixes)
      section of the changelog?

<!--
Also, please consider reading the contributing guidelines and code of conduct carefully before submitting the PR. They are available at
- https://github.com/arviz-devs/arviz/blob/main/CONTRIBUTING.md
- https://github.com/arviz-devs/arviz/blob/main/CODE_OF_CONDUCT.md

We understand that PRs can sometimes be overwhelming, especially as the
reviews start coming in.  Please let us know if the reviews are unclear or
the recommended next step seems overly demanding, if you would like help in
addressing a reviewer's comments, or if you have been waiting too long to hear
back on your PR.
-->
---
name: 'Feature Request'
about: Suggest an idea for this project

---

## Tell us about it

The more specific the better. You can assume you have access to any combination of datasets defined in the schema. When considering an implementation, keep in mind that

- models in ArviZ can have one, or lots, of random variables,
- each random variable can have very few or millions of draws,
- each random variable can have 0 or many dimensions.

## Thoughts on implementation

Not required, but if you have thoughts on how to implement the feature, that can be helpful! In case there are academic references, we welcome those here too.
---
name: 'Usage Question'
about: General questions about ArviZ usage

---

## Short Description

Let us know what you're trying to do and we can point you in the right direction. Screenshots/plots or stack traces that include `arviz` are particularly helpful.

## Code Example or link

Please provide a minimal, self-contained, and reproducible example demonstrating what you're trying to do. Ideally it will be a code snippet, link to a notebook, or link to code that can be run on another user's computer.

Also include the ArviZ version and version of any other relevant packages.

## Relevant documentation or public examples

Please provide documentation, public examples, or any additional information which may be relevant to your question
---
name: 'Bug Report'
about: Create a report to help us improve

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior. Ideally a self-contained snippet of code, or link to a notebook or external code. Please include screenshots/images produced with ArviZ here, or the stacktrace including `arviz` code to help.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Additional context**
Versions of `arviz` and other libraries used, operating system used, and anything else that may be useful.
# Scripts
Various scripts that help with ArviZ maintenance

#### container.sh
Script that helps with building docker image, and running jupyter lab,
jupyter notebooks, building docs, or getting a shell or testing

#### container.ps1
Same as container.sh

Example:

`powershell -File container.ps1 --lab`

#### create_testenv.sh
Script that helps install ArviZ on whatever computer is being used

#### lint.sh
Run all linters to check for code quality deviations

#### install_miniconda.sh
Install miniconda code environment managers and creates arviz environment

#### Dockerfile
Dockerfile for ArviZ image
(about_us)=

# About us
ArviZ is an open source project aiming to provide tools for Exploratory
Analysis of Bayesian Models that do not depend on the inference library
used. ArviZ brings together people from several probabilistic programming
libraries and from multiple programming languages. ArviZ will always be
free for everyone to use and released under the
[Apache License](https://github.com/arviz-devs/arviz/blob/main/LICENSE)
(Version 2.0).

ArviZ is developed in the open on [GitHub](https://github.com/arviz-devs/arviz)
with the input of the whole community and consensus of its core contributors.
For more information about our governance approach, read our
[Governance Document](https://github.com/arviz-devs/arviz/blob/main/GOVERNANCE.md).

## Our team
ArviZ is maintained and governed by an international group of core contributors, as defined
by our governance. You can see our current team in the cards below, as well as
in [our governance doc](https://github.com/arviz-devs/arviz/blob/main/GOVERNANCE.md#current-core-contributors).

People interested in joining the project should read
[this section](https://github.com/arviz-devs/arviz/blob/main/GOVERNANCE.md#core-contributors)
of our governance doc.
Those interested in contributing should instead refer to our {ref}`contributing_guide`.

```{include} core_contributors.md
```

## Support ArviZ
If you have found ArviZ useful in your work, research, or company, consider supporting the
project in any of the ways described in this section.

### Donate
If you or your company have the means to do so,
consider a donation to the project commensurate with your resources.
Any amount helps! All donations will be used strictly to fund the development
of ArviZ’s open source software, documentation, and community.

ArviZ is a Sponsored Project of NumFOCUS, a 501(c)(3) nonprofit charity  in the United States.
NumFOCUS provides ArviZ with fiscal, legal, and administrative support to
help ensure the health and sustainability of the project.
Visit [numfocus.org](https://numfocus.org/) for more information.

Donations to ArviZ are managed by NumFOCUS. For donors in the United States,
your gift is tax-deductible to the extent provided by law.
As with any donation, you should consult with your tax advisor about
your particular tax situation.

The ArviZ project will make the decisions on how to best use any funds received
to serve the priorities outlined in our
[Roadmap](https://github.com/arviz-devs/arviz/wiki/ArviZ-2021-roadmap).

::::{div} sd-d-flex-row sd-align-major-center
:::{button-link} https://numfocus.org/donate-to-arviz
:color: primary
:class: btn-lg

Donate to ArviZ
:::
::::

### Cite
If you use ArviZ in your scientific work, you can cite it using [![DOI](http://joss.theoj.org/papers/10.21105/joss.01143/status.svg)](https://doi.org/10.21105/joss.01143)

Here is the citation in BibTeX format

```
@article{arviz_2019,
  doi = {10.21105/joss.01143},
  url = {https://doi.org/10.21105/joss.01143},
  year = {2019},
  publisher = {The Open Journal},
  volume = {4},
  number = {33},
  pages = {1143},
  author = {Ravin Kumar and Colin Carroll and Ari Hartikainen and Osvaldo Martin},
  title = {ArviZ a unified library for exploratory analysis of Bayesian models in Python},
  journal = {Journal of Open Source Software}
}
```

Please also consider citing ArviZ's {ref}`arviz-dependencies` and
the inference library used to build and fit the model.

### ArviZ for enterprise

> ArviZ is now available as part of the Tidelift Subscription!

Tidelift is working with ArviZ and the maintainers of thousands of other open source projects
to deliver commercial support and maintenance for the open source dependencies
you use to build your applications.
Save time, reduce risk, and improve code health, while contributing financially to ArviZ.

::::{div} sd-d-flex-row sd-align-major-center
:::{button-link} https://tidelift.com/subscription/pkg/pypi-arviz?utm_source=pypi-arviz&utm_medium=referral&utm_campaign=enterprise
:color: warning

Learn more
:::
::::

#### Enterprise-ready open source software — managed for you

The Tidelift Subscription is a managed open source subscription for application
dependencies covering millions of open source projects across JavaScript, Python, Java,
PHP, Ruby, .NET, and more. And now, your favorite probabilistic programming language is included in the Tidelift subscription!

Your subscription includes:

* **Security updates**: Tidelift’s security response team coordinates patches for new breaking security vulnerabilities and alerts immediately through a private channel, so your software supply chain is always secure.

* **Licensing verification and indemnification**: Tidelift verifies license information to enable easy policy enforcement and adds intellectual property indemnification to cover creators and users in case something goes wrong. You always have a 100% up-to-date bill of materials for your dependencies to share with your legal team, customers, or partners.

* **Maintenance and code improvement**: Tidelift ensures the software you rely on keeps working as long as you need it to work. Your managed dependencies are actively maintained and Tidelift recruits additional maintainers where required.

* **Package selection and version guidance**: Tidelift helps you choose the best open source packages from the start—and then guides you through updates to stay on the best releases as new issues arise.

* **Roadmap input**: Take a seat at the table with the creators behind the software you use. ArviZ core contributors and other Tidelift’s participating maintainers earn more income as our software is used by more subscribers, so we’re interested in knowing what you need.

* **Tooling and cloud integration**: Tidelift works with GitHub, GitLab, BitBucket, and more. It supports every cloud platform (and other deployment targets, too).

The end result? All of the capabilities you expect from commercial-grade software, for the full breadth of open source you use. That means less time grappling with esoteric open source trivia, and more time building your own applications — and your business.

::::{div} sd-d-flex-row sd-align-major-center
:::{button-link} https://tidelift.com/subscription/pkg/pypi-arviz?utm_source=pypi-arviz&utm_medium=referral&utm_campaign=enterprise
:color: warning

Learn more
:::
::::

## Sponsors
[![NumFOCUS](https://www.numfocus.org/wp-content/uploads/2017/07/NumFocus_LRG.png)](https://numfocus.org)
:::::{grid} 1 2 3 3

::::{grid-item-card}
:img-top: https://raw.githubusercontent.com/pymc-devs/pymc/main/docs/logos/svg/PyMC_square.svg
:link: https://github.com/pymc-devs/pymc-examples

**PyMC example notebooks**

Curated selection of PyMC guides and case studies.
+++
Examples authored by PyMC developers and expert users.
::::


::::{grid-item-card}
:img-top: https://raw.githubusercontent.com/stan-dev/logos/master/logo.png
:link: https://mc-stan.org/users/documentation/case-studies

**Stan case studies**

The case studies on this page are intended to reflect best practices in Bayesian methodology and Stan programming.
+++
Examples authored by Stan developers and expert users.
::::

::::{grid-item-card}
:img-top: _static/images/bambi.png
:link: https://github.com/bambinos/Bambi_resources

**Bambi educational resources**

Bambi is a high-level Bayesian model-building interface written in Python,
designed to make it extremely easy to fit mixed-effects models.
+++
Tutorials and book translations maintained by the Bambi team
::::

::::{grid-item-card}
:img-top: https://oriolabril.github.io/oriol_unraveled/images/nb/labeled_arys.png
:link: https://oriolabril.github.io/oriol_unraveled/categories/#arviz

**Oriol unraveled**

Blog mostly on statistics, diversity and open source. I also use the blog as a playground
ground to test new pages for ArviZ docs.
+++
Oriol Abril Pla is an ArviZ core contributor pursuing his PhD
at Helsinki University
::::
:::::

:::::{grid} 2 3 3 4

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/23738400

Oriol Abril-Pla

[{fab}`github`](https://github.com/oriolabril)
[{fab}`twitter`](https://twitter.com/OriolAbril)
[{fas}`globe`](https://oriolabril.github.io/oriol_unraveled)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/8028618

Agustina Arroyuelo

[{fab}`github`](https://github.com/agustinaarroyuelo)
[{fab}`twitter`](https://twitter.com/AgustinaArroyu1)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/30447180

Alex Andorra

[{fab}`github`](https://github.com/AlexAndorra)
[{fab}`twitter`](https://twitter.com/alex_andorra)
[{fas}`globe`](https://www.learnbayesstats.com/)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/8673634

Seth Axen

[{fab}`github`](https://github.com/sethaxen)
[{fab}`twitter`](https://twitter.com/sethaxen)
[{fas}`globe`](https://sethaxen.com/)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/2295568

Colin Carroll

[{fab}`github`](https://github.com/ColCarroll)
[{fab}`twitter`](https://twitter.com/colindcarroll)
[{fas}`globe`](https://colindcarroll.com/)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/20489158

Piyush Gautam

[{fab}`github`](https://github.com/percygautam)
[{fab}`twitter`](https://twitter.com/percygautam)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/3274

Robert P. Goldman

[{fab}`github`](https://github.com/rpgoldman)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/13161958

Ari Hartikainen

[{fab}`github`](https://github.com/ahartikainen)
[{fab}`twitter`](https://twitter.com/a_hartikainen)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/7213793

Ravin Kumar

[{fab}`github`](https://github.com/canyon289)
[{fab}`twitter`](https://twitter.com/canyon289)
[{fas}`globe`](http://ravinkumar.com)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/1338958

Osvaldo Martin

[{fab}`github`](https://github.com/aloctavodia)
[{fab}`twitter`](https://twitter.com/aloctavodia)
[{fas}`globe`](https://aloctavodia.github.io/)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/5466612

Mitzi Morris

[{fab}`github`](https://github.com/mitzimorris)
[{fab}`twitter`](https://twitter.com/searchbooster_)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/4736342

Du Phan

[{fab}`github`](https://github.com/fehiepsi)
[{fab}`twitter`](https://twitter.com/fehiepsi)
[{fas}`globe`](http://fehiepsi.github.io/)
::::

::::{grid-item-card}
:class-item: p-2
:img-top: https://avatars.githubusercontent.com/u/6705400

Aki Vehtari

[{fab}`github`](https://github.com/avehtari)
[{fab}`twitter`](https://twitter.com/avehtari)
[{fas}`globe`](https://users.aalto.fi/~ave/)
::::

:::::
(community)=
# Community
ArviZ is a community-driven open source project committed to host an open,
inclusive and positive community. Read the
[ArviZ Code of Conduct](https://github.com/arviz-devs/arviz/blob/main/CODE_OF_CONDUCT.md)
for guidance on how to interact with each other and make the community thrive.

Our main goal is to provide backend-agnostic tools for diagnostics and visualizations of Bayesian
inference.

We do this by:

* Creating and maintaining the ArviZ library with functions for posterior analysis, data storage,
  sample diagnostics, model checking, and comparison.
* Creating social infrastructure through a welcoming and diverse community
* Making the right data, tools and best practices more discoverable
* Promoting advocacy for a culture of principled iterative modeling based on the Bayesian Workflow

This page aims to gather resources related to Bayesian modeling and ArviZ from a variety of
sources to work towards the 2nd and 3rd bullet points.
This page covers how can someone participate in the ArviZ community (this someone could be you!),
gives an overview of the ecosystem of Python libraries for Bayesian analysis and serves
as a compilation for community resources around ArviZ.

*But what is the ArviZ community?*
The ArviZ community is a self-identifying group composed of probabilistic programming practitioners who,
together, contribute to the technical and social infrastructure for open and reproducible research.
Our specific focus is on Bayesian modeling and best practices that lower the barriers to working using
a fully fledged Bayesian modeling workflow, as opposed to a rigid toolbox approach.

Community members are people who use, cite and share use cases for ArviZ,
write posts or give talks about using ArviZ,
participate in issues to help define the roadmap and improve the project,
or answer questions in any of the affine forums.

This page is part of the Python ArviZ library and is therefore specific to Python,
but the ArviZ community is not restricted to Python and in fact, aims to act as a bridge
between programming languages and encourage collaboration.

## Participate online
There are many alternatives and channels available to interact with other ArviZ users.

### Twitter
To begin with, you may want to follow us on Twitter! We are [@arviz_devs](https://twitter.com/arviz_devs).
It is the best way to be up to date with the latest developments and events related
to ArviZ.

### Gitter
The chat at [Gitter](https://gitter.im/arviz-devs/community) is a great space
to ask quick questions about ArviZ and to chat with other users and contributors.
For longer questions, the Discourse forums listed below are probably a better platform.

### Affine forums
Many ArviZ contributors are also active in one of [PyMC3 Discourse](https://discourse.pymc.io/)
or [Stan Discourse](https://discourse.mc-stan.org/) (and sometimes even in both!).

## Conferences
* [StanCon](https://mc-stan.org/events/)
* [PyMCon](https://pymcon.com)

# The Bayesian Python ecosystem
In the last years, many libraries for Bayesian data analysis have been created,
and there is a slight tendency towards more modular libraries. ArviZ plays
an important role in this ecosystem both as the go-to library for visualization
and diagnostics in Python for any and all {abbr}`PPL (Probabilistic programming library)`s and as a way to standardize and
share the results of PPL libraries.

The PPLs that integrate with ArviZ are:

* [PyMC3](https://docs.pymc.io)
* [Stan](https://mc-stan.org/users/documentation/)
* [MCX](https://github.com/rlouf/mcx)
* [Pyro](https://pyro.ai/) and [NumPyro](https://pyro.ai/numpyro/)
* [PyJAGS](https://pypi.org/project/pyjags/)
* [TensorFlow Probability](https://www.tensorflow.org/probability)

Moreover, there are other libraries that use ArviZ for visualization and diagnostics
and/or that are compatible with `InfereceData` objects:

* [Bambi](https://bambinos.github.io/bambi/)
* [corner.py](https://corner.readthedocs.io/en/latest/)

# Community educational resources
ArviZ is a transversal and backend agnostic project. One can use ArviZ with _any_ PPL,
but one has to be used in order to have some data for ArviZ to analyze.
This makes writing detailed and complete documentation for ArviZ complicated.
This section aims to act as a compilation of community resources related to ArviZ
that hopefully can bridge the gap. These can be books, case studies in the documentation of
PPLs, personal blogs...

Do you know of resources about ArviZ that are not listed here? Open an issue or a PR and
let us know!

## Books
* Bayesian Modeling and Computation in Python: available soon
* [Bayesian Analysis with Python](https://github.com/aloctavodia/BAP)
* [Bayesian Data Analysis 3](http://www.stat.columbia.edu/~gelman/book/)

## Podcasts
* [Learning Bayesian Statistics](https://www.learnbayesstats.com/)
* [dats'n'stats](https://www.pydata-podcast.com/)

## Blogs and example compilations
If you have a blog that has 2 or more posts tagged as "ArviZ", you can submit
a pull request to add your blog to this list

```{include} external_resources.md
```
(schema)=
# InferenceData schema specification
The `InferenceData` schema approach defines a data structure compatible with [NetCDF](https://www.unidata.ucar.edu/software/netcdf/). Its purpose is to serve the following three goals:
1. Usefulness in the analysis of Bayesian inference results.
2. Reproducibility of Bayesian inference analysis.
3. Interoperability between different inference backends and programming languages.

Currently there are **two beta implementations** of this design:
* [ArviZ](https://arviz-devs.github.io/arviz/) in **Python** which integrates with:
  - [emcee](https://emcee.readthedocs.io/en/stable/)
  - [PyMC3](https://docs.pymc.io)
  - [Pyro](https://pyro.ai/) and [NumPyro](https://pyro.ai/numpyro/)
  - [PyStan](https://pystan.readthedocs.io/en/latest/index.html), [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) and [CmdStanPy](https://cmdstanpy.readthedocs.io/en/latest/index.html)
  - [TensorFlow Probability](https://www.tensorflow.org/probability)
* [ArviZ.jl](https://github.com/arviz-devs/ArviZ.jl) in **Julia** which integrates with:
  - [CmdStan.jl](https://github.com/StanJulia/CmdStan.jl), [StanSample.jl](https://github.com/StanJulia/StanSample.jl) and [Stan.jl](https://github.com/StanJulia/Stan.jl)
  - [Turing.jl](https://turing.ml/dev/) and indirectly any package using [MCMCChains.jl](https://github.com/TuringLang/MCMCChains.jl) to store results

## Terminology
The terminology used in this specification is based on [xarray's terminology](http://xarray.pydata.org/en/stable/terminology.html), however, no xarray knowledge is assumed in this description. There are also some extensions particular to  the {ref}`InferenceData <xarray_for_arviz>` case.

* **Variable**: NetCDF-like variables are multidimensional labeled arrays representing a single quantity. Variables and their dimensions must be named. They can also have attributes describing it. Relevant terms related to `InferenceData` variables are following:
  - *variable_name*
  - *values* (its data)
  - *dimensions*
  - *coordinates*
  - *attributes*
* **Dimension**: The dimensions of an object are its named axes. A variable containing 3D data can have dimensions `[chain, draw, dim0]`, i.e., its `0th`-dimension is `chain`, its `1st`-dimension is `draw`, and so on. Every dimension present in an `InferenceData` variable must share names with a *coordinate*. Given that dimensions must be named, dimension and dimension name are used equivalents.
* **Coordinate**: A named array that labels a dimension. A coordinate named `chain` with values `[0, 1, 2, 3]` would label the `chain` dimension. Coordinate names and values can be loosely thought of as labels and tick labels along a dimension, respectively.
* **Attribute**: An ordered dictionary that can store arbitrary metadata.
* **Group**: Dataset containing one or several variables with a conceptual link between them. Variables inside a group will generally share some dimensions too. For example, the `posterior` group contains a representation of the posterior distribution conditioned on the observations in the `observed_data` group.
* **Matching samples**: Two variables (or groups) will be called to have matching samples if they are generated with the same set of samples. Therefore, they will share dimensions and coordinates corresponding to the sampling process. Sample dimensions (generally `(chain, draw)`) are the ones introduced by the sampling process.
* **Matching variables**: Two groups with matching variables are groups that conceptually share variables, variable dimensions and coordinates of the variable dimensions but do not necessarily share variable names nor sample dimensions. Variable dimensions are the ones intrinsic to the data and model as opposed to sample dimensions which are the ones relative to the sampling process. When talking about specific variables, this same idea is expressed as one variable being the counterpart of the other.

## Current design
`InferenceData` stores all quantities that are relevant to fulfilling its goals in different groups. Different groups generally distinguish conceptually different quantities in Bayesian inference, however, convenience in {ref}`creation <creating_InferenceData>` and {ref}`usage <working_with_InferenceData>` of `InferenceData` objects also plays a role. In general, each quantity (such as posterior distribution or observed data) will be represented by several multidimensional labeled variables.

### Rules
Following are a few rules which should be followed:
* Each group should have one entry per variable and each variable should be named.
* When relevant, the first two dimensions of each variable should be the sample identifier (`chain`, `draw`).
* For groups like `observed_data` or `constant_data`, the two initial dimensions(`chain`, `draw`) are omitted.
* Dimensions must be named and share name with a coordinate specifying the index values, called coordinate values.
* Coordinate values can be repeated and should not necessarily be numerical values.
* Variables must not share names with dimensions.
* Moreover, each group contains the following attributes:
  - `created_at`: the date of creation of the group.
  - `inference_library`: the library used to run the inference.
  - `inference_library_version`: version of the inference library used.

### Relations
`InferenceData` data objects contain any combination of the groups described below. There are also some relations (detailed below) between the variables and dimensions of different groups. Hence, whenever related groups are present they should comply with these relations.

#### `posterior`
Samples from the posterior distribution p(theta|y).

#### `sample_stats`
Information and diagnostics for each `posterior` sample, provided by the inference
backend. It may vary depending on the algorithm used by the backend (i.e. an affine
invariant sampler has no energy associated). Therefore none of these parameters
should be assumed to be present in the `sample_stats` group. The convention
below serves to ensure that if a variable is present with one of these names
it will correspond to the definition given in front of it.

The name convention used for `sample_stats` variables is the following:

* `lp`: The joint log posterior density for the model (up to an additive constant).
* `acceptance_rate`: The average acceptance probabilities of all possible samples in the proposed tree.
* `step_size`: The current integration step size.
* `step_size_nom`: The nominal integration step size. The `step_size` may differ from this for example, if the step size is jittered. It should only be present if `step_size` is also present and it varies between samples (i.e. step size is jittered).
* `tree_depth`: The number of tree doublings in the balanced binary tree.
* `n_steps`: The number of leapfrog steps computed. It is related to `tree_depth` with `n_steps <=
  2^tree_dept`.
* `diverging`: (boolean) Indicates the presence of leapfrog transitions with large energy deviation
  from starting and subsequent termination of the trajectory. "large" is defined as `max_energy_error` going over a threshold.
* `energy`: The value of the Hamiltonian energy for the accepted proposal (up to an
additive constant).
* `energy_error`: The difference in the Hamiltonian energy between the initial point and
the accepted proposal.
* `max_energy_error`: The maximum absolute difference in Hamiltonian energy between the initial point and all possible samples in the proposed tree.
* `int_time`: The total integration time (static HMC sampler)


#### `log_likelihood`
Pointwise log likelihood data. Samples should match with `posterior` ones and its variables
should match `observed_data` variables. The `observed_data` counterpart variable
may have a different name. Moreover, some cases such as a multivariate normal
may require some dimensions or coordinates to be different.

#### `posterior_predictive`
Posterior predictive samples p(y|y) corresponding to the posterior predictive distribution evaluated at the `observed_data`. Samples should match with `posterior` ones and its variables should match `observed_data` variables. The `observed_data` counterpart variable may have a different name.

#### `observed_data`
Observed data on which the `posterior` is conditional. It should only contain data which is modeled as a random variable. Each variable should have a counterpart in `posterior_predictive`, however, the `posterior_predictive` counterpart variable may have a different name.

#### `constant_data`
Model constants, data included in the model which is not modeled as a random variable. It should be the data used to generate samples in all the groups except the `predictions` groups.

#### `prior`
Samples from the prior distribution p(theta). Samples do not need to match `posterior` samples. However, this group will still follow the convention on `chain` and `draw` as first dimensions. It should have matching variables with the `posterior` group.

#### `sample_stats_prior`
Information and diagnostics for the samples in the `prior` group, provided by the inference backend. It may vary depending on the algorithm used by the backend. Variable names follow the same convention defined in `sample_stats`.

#### `prior_predictive`
Samples from the prior predictive distribution. Samples should match `prior` samples and each variable should have a counterpart in `posterior_predictive`/`observed_data`.

#### `predictions`
Out of sample posterior predictive samples p(y'|y). Samples should match `posterior` samples. Its variables should have a counterpart in `posterior_predictive`. However, variables in `predictions` and their counterpart in `posterior_predictive` can have different coordinate values.

#### `predictions_constant_data`
Model constants used to get the `predictions` samples. Its variables should have a counterpart in `constant_data`. However, variables in `predictions_constant_data` and their counterpart in `constant_data` can have different coordinate values.

## Planned features
The `InferenceData` structure is still evolving, with some feature being currently developed. This section aims to describe the roadmap of the specification.

### Sampler parameters
Parameters of the sampling algorithm and sampling backend to be used for analysis reproducibility.

## Examples
In order to clarify the definitions above, an example of `InferenceData` generation for a 1D linear regression is available in several programming languages and probabilistic programming frameworks. This particular inference task has been chosen because it is widely well known while still being useful and it also allows to populate all the fields in the `InferenceData` object.

### Python

```{toctree}
PyMC3 example <PyMC3_schema_example>
PyStan example <PyStan_schema_example>
```
(diataxis_for_arviz)=
# Diátaxis for ArviZ

ArviZ documentation has the {doc}`diataxis:index` as a North Star.

There are however some specificities to ArviZ case and knowing Diátaxis alone
is not enough to be able to understand ArviZ doc structure and correctly place
new content pages within it.
At the same time, having a basic understanding of Diátaxis is needed to understand
this page and the content structure in ArviZ docs.

These are:

* Multiple target audiences and documentation goals
* Strong visual component
* Work in progress

## Multiple target audiences
Diátaxis is a framework based around **user needs** in their exercise of a
**practical craft** (more detail at {doc}`diataxis:introduction` page in Diatáxis website).

The scope of the ArviZ website however is not limited to ArviZ users nor it is limited to practical
crafts like using ArviZ or contributing to it.

After careful consideration, we decided to split the website into three components:
two practical crafts, using and contributing
(each of which follows Diátaxis independently of the other)
and general information about the ArviZ project,
it's people, community, marketing info...
which can nor should follow the Diátaxis framework.

## Strong visual component
ArviZ is a library with a strong visualization component.
This translates for example into some duplicity when it comes to reference content.

Both the {ref}`example_gallery` and the {ref}`api` pages are reference content.
The API reference indexes by name and contains the description of each object,
whereas the Example Gallery indexes by image and to avoid excessive duplication,
links to the respective API page.

This allows users who know the name of the function or who want to check the options
available for the function they are currently using to quickly find it in the
API reference while also allowing users who know how a plot looks but don't know
their name nor the right ArviZ function to find it in the Example Gallery.

## Work in progress
In addition, ArviZ documentation itself and the application of Diátaxis to it
is still an active work in progress. That means that some pages will not
match what is described in {ref}`content_structure` because their updating to
enforce Diátaxis is still pending.
(contributing_guide)=
# Contributing to ArviZ
As a scientific, community-driven, open source software project,
ArviZ welcomes contributions from interested individuals or groups.
These guidelines are provided to give potential contributors information
to make their contributions compliant with the conventions of the ArviZ project,
and maximize the probability of such contributions to be merged as quickly
and efficiently as possible.

Even though contributing code or documentation by sending pull requests on
the main [ArviZ GitHub repository](https://github.com/arviz-devs/arviz) is the most common way of contributing and
therefore this guide dedicates many resources to these, ArviZ welcomes any kind of contribution.
As a clear proof, you can become a core contributor and participate in
[ArviZ governance](https://github.com/arviz-devs/arviz/blob/main/GOVERNANCE.md)
without having contributed code to ArviZ.
All contributions to ArviZ, either on its codebase, its documentation,
or its community are valuable and we really appreciate your help.

:::{tip}
Contact us on [Gitter](https://gitter.im/arviz-devs/community) if you want to
contribute to the project but you are not sure where you can contribute or how to start.
:::

Below we list some examples of contributions that might serve as inspiration.

* [Submitting issues](https://github.com/arviz-devs/arviz/issues/new/choose) related to bugs or desired enhancements. Check the details at {doc}`issue_reports`.
* Fixing outstanding [issues](https://github.com/arviz-devs/arviz/issues) (bugs) with the existing codebase. They range from low-level software bugs to higher-level design problems.
* Contributing or improving the [documentation](https://arviz-devs.github.io/arviz/) (`docs`) or examples (`arviz/examples`).
* Help users submitting issues by making sure their issue is clear and/or providing temporal
  workarounds for them. Check the details at {ref}`issue_triaging`
* Review Pull Requests ensuring contributions are well tested and documented, that documentation
  renders correctly and has no typos or mistakes. Check the details at {ref}`review_prs`
* Adding new or improved functionality to the existing codebase.

Further information on how to contribute to ArviZ can be found on these specific pages:

```{toctree}
:maxdepth: 2
:caption: Contribution types overview

issue_reports
issue_triaging
outreach
review_prs
contributing_prs
```

:::{toctree}
:hidden:
:caption: Tutorials

pr_tutorial
:::

:::{toctree}
:hidden:
:caption: How-to guides

sphinx_doc_build
running_benchmarks
:::

:::{toctree}
:hidden:
:caption: Reference

pr_checklist
architecture
content_structure
docstrings
syntax_guide
developing_in_docker
:::

:::{toctree}
:hidden:
:caption: In depth explanations

doc_toolchain
diataxis_for_arviz
plotting_backends
:::
(pr_checklist)=
# Pull request checklist

We recommend that your contribution complies with the following guidelines before you submit a pull request:

* If your pull request addresses an issue, please use the pull request title to describe the issue and mention the issue number in the pull request description. This will make sure a link back to the original issue is created.

* All public methods must have informative docstrings with sample usage when appropriate.

* Please prefix the title of incomplete contributions with `[WIP]` (to indicate a work in progress). WIPs may be useful to (1) indicate you are working on something to avoid duplicated work, (2) request a broad review of functionality or API, or (3) seek collaborators.

* All other tests pass when everything is rebuilt from scratch.
See {ref}`developing_in_docker` for information on running the test suite locally.

* When adding additional plotting functionality, provide at least one example script in the [arviz/examples/](https://github.com/arviz-devs/arviz/tree/main/examples) folder. Have a look at other examples for reference. Examples should demonstrate why the new functionality is useful in practice and, if possible, compare it to other methods available in ArviZ.

* Added tests follow the [pytest fixture pattern](https://docs.pytest.org/en/latest/fixture.html#fixture).

* Documentation and high-coverage tests are necessary for enhancements to be accepted.

* Documentation follows Numpy style guide.

* Run any of the pre-existing examples in ``docs/source/notebooks`` that contain analyses that would be affected by your changes to ensure that nothing breaks. This is a useful opportunity to not only check your work for bugs that might not be revealed by unit test, but also to show how your contribution improves ArviZ for end users.

* If modifying a plot, render your plot to inspect for changes and copy image in the pull request message on Github.

You can also check for common programming errors with the following
tools:

* Save plots as part of tests. Plots will be saved to a directory named `test_images` by default.

  ```bash
  $ pytest arviz/tests/base_tests/<name of test>.py --save
  ```

* Optionally save plots to a user named directory. This is useful for comparing changes across branches.

  ```bash
  $ pytest arviz/tests/base_tests/<name of test>.py --save user_defined_directory
  ```


* Code coverage **cannot** decrease. Coverage can be checked with **pytest-cov** package:

  ```bash
  $ pip install pytest pytest-cov coverage
  $ pytest --cov=arviz --cov-report=html arviz/tests/
  ```

* Your code has been formatted with [black](https://github.com/ambv/black) with a line length of 100 characters.

  ```bash
  $ pip install black
  $ black arviz/ examples/ asv_benchmarks/
  ```

* Your code passes pylint

  ```bash
  $ pip install pylint
  $ pylint arviz/
  ```

* No code style warnings, check with:

  ```bash
  $ ./scripts/lint.sh
  ```
# Issue reports

We appreciate being notified of problems with the existing ArviZ code.
We prefer that issues be filed on the
[Github Issue Tracker](https://github.com/arviz-devs/arviz/issues),
rather than on social media or by direct email to the developers.

Please verify that your issue is not being currently addressed by other
issues or pull requests by using the GitHub search tool to look for keywords
in the project issue tracker.

When writing an issue, make sure to include any supporting information,
in particular, the version of ArviZ that you are using and how did you install it.
The issue tracker has several templates available to help in writing the issue
and including useful supporting information.

## Minimal reproducible example
If your issue reports a bug or tries to show a specific behaviour,
provide a way to reproduce the issue. Consider the advice on
[stackoverflow](https://stackoverflow.com/help/minimal-reproducible-example)
about writing a minimal reproducible code example.

If even a minimal version of the code is still too long to comfortably share
in the issue text, upload the code snippet as a
[GitHub Gist](https://gist.github.com/) or a standalone Jupyter notebook.

## Sharing error tracebacks
If relevant, also include the complete error messages and tracebacks.
To avoid cluttering the issue body, you can optionally use the `details`
html tag in GitHub to create a simple dropdown, like so:

````
<details><summary>Click to see full traceback</summary>

```
long and barely comprehensible error traceback
```

</details>
````(sphinx_doc_build)=
# Building the documentation

The preferred workflow for contributing to ArviZ is to fork
the [GitHub repository](https://github.com/arviz-devs/arviz/),
clone it to your local machine, and develop on a feature branch. For a detailed description of the recommended development process, see the step-by-step guide {ref}`here <pr_steps>`.

If you are using Docker, see {ref}`building_doc_with_docker`.

## Pull request checks

Every PR has a list of checks that check if your changes follow the rules being followed by ArviZ docs.
You can check why a specific test is failing by clicking the `Details` next to it. It will take you to errors and warning page. This page shows the details of errors, for example in case of docstrings:

`arviz/plots/pairplot.py:127:0: C0301: Line too long (142/100) (line-too-long)`.

It means line 127 of the file `pairplot.py` is too long.
For running tests locally, see the {ref}`pr_checklist`.


(preview_change)=
## Previewing doc changes

There is an easy way to check the preview of docs by opening a PR on GitHub. ArviZ uses `readthedocs` to automatically build the documentation.
For previewing documentation changes, take the following steps:

1. Go to the checks of your PR. Wait for the `docs/readthedocs.org:arviz` to complete.

   ```{note}
   The green tick indicates that the given check has been built successfully.
   ```

2. Click the `Details` button next to it.
3. It will take you to the preview of ArviZ docs of your PR.
4. Go to the webpage of the file you are working on.
5. Check the preview of your changes on that page.

```{note} Note
The preview version of ArviZ docs will have a warning box that says "This page was created from a pull request (#Your PR number)." It shows the PR number whose changes have been implemented.
```

For example, a warning box will look like this:

```{warning}
This page was created from a pull request (#PR Number).
```
(developing_in_docker)=
# Developing in Docker

We have provided a Dockerfile which helps for isolating build problems, and local development.
Install [Docker](https://www.docker.com/) for your operating system, clone [this](https://github.com/arviz-devs/arviz) repo. Docker will generate an environment with your local copy of `arviz` with all the packages in Dockerfile.

## container.sh & container.ps1

Predefined docker commands can be run with a `./scripts/container.sh` (on Linux and macOS)
and with `./scripts/container.ps1`. The scripts enable developers easily to call predefined docker commands.
Users can use one or multiple flags.

They are executed in the following order: clear-cache, build, test, docs, shell, notebook, lab

### For Linux and macOS
    $ ./scripts/container.sh --clear-cache
    $ ./scripts/container.sh --build

    $ ./scripts/container.sh --test
    $ ./scripts/container.sh --docs
    $ ./scripts/container.sh --shell
    $ ./scripts/container.sh --notebook
    $ ./scripts/container.sh --lab

### For Windows
    $ powershell.exe -File ./scripts/container.ps1 --clear-cache
    $ powershell.exe -File ./scripts/container.ps1 --build

    $ powershell.exe -File ./scripts/container.ps1 --test
    $ powershell.exe -File ./scripts/container.ps1 --docs
    $ powershell.exe -File ./scripts/container.ps1 --shell
    $ powershell.exe -File ./scripts/container.ps1 --notebook
    $ powershell.exe -File ./scripts/container.ps1 --lab

## Testing in Docker
Testing the code using docker consists of executing the same file(`./scripts/container.ps1`) 3 times (you may need root privileges to run it).
1. First run `./scripts/container.sh --clear-cache`.
2. Then run `./scripts/container.sh --build`. This starts a local docker image called `arviz`.
3. Finally, run the tests with `./scripts/container.sh --test`.

The `./scripts/container.sh --build` command needs to be executed only once to build the docker image. If you have already built the docker image, and want to test your new changes, then it is a 2 step process.
1. Run `./scripts/container.sh --clear-cache` because the cache generated locally is not valid for docker.
2. Run `./scripts/container.sh --test` to test with docker.

**NOTE**: If you run into errors due to `__pycache__` files (i.e. while testing in
docker after testing locally or installing with pip after testing with
docker), try running `./scripts/container.sh --clear-cache` before the errored
command.

## Using the Docker image interactively
Once the Docker image is built with `./scripts/container.sh --build`, interactive containers can also be run. Therefore, code can be edited and executed using the docker container, but modifying directly the working directory of the host machine.

### For Linux and macOS
To start a bash shell inside Docker, run:

    $ docker run --mount type=bind,source="$(pwd)",target=/opt/arviz/ -it arviz bash

Alternatively, to start a jupyter notebook, there are two steps, first run:

    $ docker run --mount type=bind,source="$(pwd)",target=/opt/arviz/ --name jupyter-dock -it -d -p 8888:8888 arviz
    $ docker exec jupyter-dock bash -c "pip install jupyter"
    $ docker exec -it jupyter-dock bash -c "jupyter notebook --ip 0.0.0.0 --no-browser --allow-root"

### For Windows
For Windows, use %CD% on cmd.exe and $pwd.Path on powershell.

    $ docker run --mount type=bind,source=%CD%,target=/opt/arviz/ -it arviz bash

For running a jupyter notebook, run:

    $ docker run --mount type=bind,source=%CD%,target=/opt/arviz/ --name jupyter-dock -it -d -p 8888:8888 arviz
    $ docker exec jupyter-dock bash -c "pip install jupyter"
    $ docker exec -it jupyter-dock bash -c "jupyter notebook --ip 0.0.0.0 --no-browser --allow-root"

This will output something similar to `http://(<docker container id> or <ip>):8888/?token=<token id>`, and can be accessed at `http://localhost:8888/?token=<token id>`.

(building_doc_with_docker)=
## Building documentation with Docker
The documentation can be built with Docker by running `./scripts/container.sh --docs`. The
`./scripts/container.sh --build` needs to be executed once before running the doc build command.
The docker image contains by default all dependencies needed
for building the documentation. After having build the docs in the Docker
container, they can be checked at `doc/build` and viewed in the browser from `doc/build/index.html`.

To update any file, go to the specific file(`.md`, `.rst`, or `.ipynb`) and make your changes. In order to check the results, you can either check the {ref}`preview_change` or rebuild the docs with Docker. For the latter option, You need to only rebuild the docs, using the same command `./scripts/container.sh --docs`. The resulting docs will be saved in the `doc/build` directory as HTML files. Open your corresponding file and check the results. You wil also need to rebuild the docker image if extensions and dependencies were added or changed.(oureach_contrib)=
# Outreach

:::{caution} Page in construction
:::

You can also contribute to ArviZ by doing outreach, writing blogposts, getting
people to use it at your company or university, helping external libraries integrate with ArviZ (in
any programming language of your choice).

If you have a website or blog where you often post about ArviZ, a library that integrates with ArviZ reach out to us or
open a PR to add your blog to {ref}`community`, we'd love to hear from you!
(issue_triaging)=
# Issue Triaging


As a community-driven, open source collective, we value all user contributions at ArviZ. While most users often tend to look for ways to contribute code via PRs, we welcome and encourage help with issue triaging.

We consider issue triaging an integral part of the library development. It is the main communication
channel between contributors and users of ArviZ. Moreover, it is a task that gives contributors
more flexibility than contributing code.

This page aims to describe how to get started as an issue triager with ArviZ and provide some
guidelines related to the task.

:::{note}
Even if it is often ignored, issue triaging is generally regarded as an important task,
and it is getting more traction lately. A clear example of that is the
[Labelathon](https://ropensci.org/commcalls/apr2021-pkg-community/) recently hosted by
the rOpenSci community. The [Labelathon was recorded](https://ropensci.org/commcalls/apr2021-pkg-community/)
and has also had a couple of [follow-up blogpost](https://ropensci.org/tags/labelathon/) which are great resources related to
issue triaging.
:::

## How to subscribe to ArviZ issues

The first step for you to get started would be to make sure you subscribe and get notifications of any new issue published on ArviZ's GitHub repo. You can then chose to work on it at the moment or come back to it later. Staying subscribed will ensure you are notified whenever issues are created or updated and therefore, see where your help might be needed.

If you are not familiar on how to set up notifications on GitHub, please check the following - [Setting Up Notifications on GitHub](https://docs.github.com/en/github/managing-subscriptions-and-notifications-on-github/setting-up-notifications/configuring-notifications#configuring-your-watch-settings-for-an-individual-repository).
Once you are set up, you can [view your subscriptions](https://docs.github.com/en/github/managing-subscriptions-and-notifications-on-github/managing-subscriptions-for-activity-on-github/viewing-your-subscriptions) and [manage your subscritions](https://docs.github.com/en/github/managing-subscriptions-and-notifications-on-github/managing-subscriptions-for-activity-on-github/managing-your-subscriptions) to ensure you are not being inundated with the volume and are getting notifications only to issues of your interest.

## Triage guidelines and suggestions

Similar to contributing code via PRs, most issue triaging tasks don't require any
specific permissions. Anyone with a GitHub account can help with issue triaging.

:::{important}
The list below provides ideas and examples of issue triaging entails. However, it is not a comprehensive compilation. Often users encounter issues not foreseen or experienced by developers. We encourage users go ahead and take ownership and bring these to the attention of the person who have posted or the contributors working on it.
:::

Make sure the issue contains a [minimal reproducible example](https://stackoverflow.com/help/minimal-reproducible-example), if relevant.
: Sometimes, the issue doesn't contain an example, however, it can still be clear about the problem.
  In that scenario, someone other than the person who posted the issue can generate an example.
  Issues with a reproducible example allow contributors to focus on fixing the bug or testing the
  proposed enhancement directly instead of having to first understand and reproduce the issue. This
  therefore makes things easier for contributors, but at the same time will reduce the time it takes
  to answer or fix an issue, helping issue posters directly too.

Ensure the issue is clear and has references, if needed.
: If an issue is not completely clear, you can comment by asking for clarifications or
  adding extra references so the issue is clear enough for someone else to start working on.
  One example would be [#1694 (comment)](https://github.com/arviz-devs/arviz/issues/1694#issuecomment-840683745)

Suggest fixes or workarounds.
: In some cases, the issue might be a usage question more than an issue (for example
  [#1758](https://github.com/arviz-devs/arviz/issues/1758)), in which case it can be answered directly,
  in some other cases it might be a version mismatch (i.e. users expecting a fresh out of the oven
  feature to be present in the latest version when it's only available on GitHub development
  version like in [#1773](https://github.com/arviz-devs/arviz/issues/1773)).
  Or there may even be a not ideal yet acceptable workaround for users to avoid the issue
  before it's fixed (for example [#1467](https://github.com/arviz-devs/arviz/issues/1467)).

Guide newcomers
: It is also important to introduce people who comment on ArviZ issues for the first time to the
  community, as well as helping people find issues to work on that match their interest
  and abilities.

:::{note}
If you regularly help with issue triaging you'll probably be asked to join the team and be given
some triaging or write permissions on the ArviZ GitHub repo. This would allow you to label and
assign issues in the ArviZ repo.

You can read more about our team and how to join in the {ref}`about_us` page.
:::
(syntax_guide)=
# Written content formatting reference

## Adding targets
Adding custom targets or anchors to the headings is really helpful for cross-referencing. They allow us to link to the heading using a simple syntax.
They are defined using this syntax:

::::{tab-set}

:::{tab-item} rST
:sync: rst

```
.. _mytarget:
```
:::

:::{tab-item} MyST (Markdown)
:sync: myst

```
(mytarget)=
```
:::

::::

They are referred using this syntax:

::::{tab-set}

:::{tab-item} rST
:sync: rst

```
:ref:`mytarget`
```
:::

:::{tab-item} MyST (Markdown)
:sync: myst

```
{ref}`mytarget`
```
:::

::::

For adding anchors in `.ipynb` files, Markdown syntax will be used in the markdown cells of `.ipynb` files.

:::{seealso}
{ref}`sphinx:ref-role`
  Sphinx documentation about the `ref` role.

[Hyperlink targets](https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#hyperlink-targets)
  Docutils documentation on targets (another name for anchors) with extensive detail on how and when
  they can be used.

[MyST targets and cross-referencing docs](https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html#targets-and-cross-referencing)
:::

## Backticks for highlighting code keywords

For highlighting inline code keywords or file names, backticks are used. In Markdown single backticks are used while in rST double backticks are used. For example, for highlighting the file name `conf.py`, we will use this syntax:

:::::{tab-set}
::::{tab-item} rST
```
``conf.py``
```
::::
::::{tab-item} MyST (Markdown)
```
`conf.py`
```
::::
:::::

## Table of content tree
All documentation generated with sphinx is structured via `toctrees`.
Sphinx prints a warning if a file is not part of any toctree unless
that file is marked as an orphan.

The hierarchy defined via toctrees (and not the file hierarchy!) is
the one that defines the nabvar and the contents of the left sidebar.
Keeping the toctrees organized and up to date ensures the sphinx
build works, that the generated documentation is correctly ordered
and can be navigated and that all pages can be reached.
It also enables navigation to the **Previous** and **Next** pages in the footer.

Follow this syntax for adding the table of content:

:::::{tab-set}
::::{tab-item} rST
```rST
.. toctree::
    developer_guide
    doc_guide
```
::::
::::{tab-item} MyST (Markdown)
```markdown
:::{toctree}
developer_guide
doc_guide
:::
```
::::
:::::

Read more about the {ref}`Sphinx toctree directive <sphinx:toctree-directive>`.

(adding_references)=
## Adding references

In ArviZ docs, we use {confval}`sphinx:intersphinx_mapping` to add references to other libraries functions and objects.
The {mod}`~sphinx.ext.intersphinx` ensures us that cross-references to the target project exist.
It is the only way to link for multiversioned docs to link to themselves.
It raises a warning if the target references are changed or removed.
This way we don't need to add the long exact links.
It saves a lot of time and effort.

### Hyperlinks
Complementary functions such as {func}`arviz.compare` and {func}`arviz.plot_compare` should reference
each other in their docstrings using a hyperlink, not only by name. The same
should happen with external functions whose usage is assumed to be known; a
clear example of this situation are docstrings on `kwargs` passed to bokeh or
matplotlib methods. This section covers how to reference functions from any
part of the docstring. Read more about it {doc}`here <sphinx:usage/restructuredtext/roles>`.

(reference_external_libs)=
#### Reference external libraries

Sphinx is configured to ease referencing libraries ArviZ relies heavily on by
using [intersphinx](https://docs.readthedocs.io/en/stable/guides/intersphinx.html).
See guidance on the reference about how to link to objects from external
libraries and the value of intersphinx_mapping in [conf.py](https://github.com/arviz-devs/arviz/blob/main/doc/source/conf.py) for the complete and up to date list of libraries that can be referenced.

In ArviZ docs, you can add references to functions and objects of `matplotlib`, `bokeh`, `xarray`, etc following the simple syntax. Let's try adding a function of few libraries, i.e., {meth}`xarray.Dataset.sel`, {func}`matplotlib.pyplot.subplots` and
{func}`bokeh.plotting.figure`.

:::::{tab-set}
::::{tab-item} rST
```
:meth:`xarray.Dataset.sel`
:func:`matplotlib.pyplot.subplots`
:func:`bokeh.plotting.figure`
```
::::
::::{tab-item} MyST (Markdown)
```
{meth}`xarray.Dataset.sel`
{func}`matplotlib.pyplot.subplots`
{func}`bokeh.plotting.figure`
```
::::
:::::

Note that the `:key:` before
the reference must match the kind of object that is being referenced, it
generally will not be `:ref:` nor `:doc:`. For
example, for functions `:func:` has to be used and for class methods
`:meth:`. The complete list of keys can be found [here](https://github.com/sphinx-doc/sphinx/blob/685e3fdb49c42b464e09ec955e1033e2a8729fff/sphinx/domains/python.py#L845-L881).

The extension [sphobjinv](https://sphobjinv.readthedocs.io/en/latest/) can
also be helpful in order to get the exact type and name of a reference. Below
is an example on getting a reference from matplotlib docs:

```
$ sphobjinv suggest -t 90 -u https://matplotlib.org/objects.inv "axes.plot"

Remote inventory found.

:py:method:`matplotlib.axes.Axes.plot`
:py:method:`matplotlib.axes.Axes.plot_date`
:std:doc:`api/_as_gen/matplotlib.axes.Axes.plot`
:std:doc:`api/_as_gen/matplotlib.axes.Axes.plot_date`
```

We can therefore link to matplotlib docs on `Axes.plot` from any docstring
using:

:::::{tab-set}
::::{tab-item} rST
```
:meth:`mpl:matplotlib.axes.Axes.plot`
```
::::
::::{tab-item} MyST (Markdown)
```
{meth}`mpl:matplotlib.axes.Axes.plot`
```
::::
:::::

The `intersphinx_mappings`
defined for ArviZ can be seen in `conf.py`.
Moreover, the intersphinx key is optional. Thus, the pattern to get sphinx to generate links is:

:::::{tab-set}
::::{tab-item} rST
```
:type_id:`(intersphinx_key:)object_id`
```
::::
::::{tab-item} MyST (Markdown)
```
{type_id}`(intersphinx_key:)object_id`
```
::::
:::::

with the part between brackets being optional. See the docstring on
{meth}`~arviz.InferenceData.to_dataframe` and
[its source](https://arviz-devs.github.io/arviz/_modules/arviz/data/inference_data.html#InferenceData.to_dataframe) for an example.

(reference_arviz_objects)=
#### Referencing ArviZ objects

The same can be done to refer to ArviZ functions, in which case,
``:func:`arviz.loo` `` is enough, there is no need to use intersphinx.
Moreover, using ``:func:`~arviz.loo` `` will only show ``loo`` as link text
due to the preceding ``~``.

## Code tabs

You will find code tabs on every other page in ArviZ docs. As we have two main types of files, i.e, `.rst` and `.md`, we often use two code tabs to show the functionalities in both rST and Markdown.

### Synchronised Tabs
ArviZ docs are using `sphinx-design` extension for adding sync code tabs in [conf.py](https://github.com/arviz-devs/arviz/blob/main/doc/source/conf.py#L61). You can check the syntax and more info about it
at [Synchronised Tabs](https://sphinx-design.readthedocs.io/en/sbt-theme/tabs.html#synchronised-tabs). Using this extension saves us from a lot of raw-html code. Sphinx provides this extension to make our work easy and or code more concise.


## Extensions
Sphinx supports a lot of extensions to improve and customize the doc features. ArviZ makes use of these builtin and external extensions to add extra roles. See the {doc}`Extensions in ArviZ <sphinx-primer:extensions/used_extensions>`.
(content_structure)=
# Documentation content structure

ArviZ documentation has the {doc}`diataxis:index` as a North Star.
This page assumes basic familiarity with Diátaxis.

The ArviZ website is divided into three content blocks:

* ArviZ user documentation
* ArviZ contributor documentation
* ArviZ Project website

## ArviZ user documentation
It is the block formed by the navbar sections
`Getting Started`, `Example Gallery`, `User Guide` and `API Reference`.

It's audience are ArviZ users and it uses Diátaxis to achieve its goal of
teaching how to use ArviZ.

Getting Started
: The Getting Started section should contain mostly tutorials, however, some
  how-to content or overview explanations are also allowed and even encouraged

Example Gallery
: The Example Gallery is a visual reference section that should contain
  **no description** of the functions used. It should instead link to the
  API Reference section.

  Its goal is to serve as visual index for users who know how their desired plot
  looks like but not the name of the ArviZ function that generates it

User Guide
: The User Guide is targeted to people who already have basic ArviZ knowledge
  and should mostly contain how-to guides.
  It can also contain in-depth explanations and some non-python objects reference
  content.

API Reference
: It should index and describe all ArviZ objects that are part of the public API,
  that is, are exposed to ArviZ users. It should be generated from a list of
  public objects only via sphinx extensions, taking the actual content from
  the docstrings in each ArviZ object.

## ArviZ contributor documentation
It is all inside the `Contributing` navbar section.

It's audience are any and all ArviZ contributors, from new contributors to members
of the core team.
It follows Diátaxis quite closely by using multiple captioned toctrees.

Contribution types overview
: Diátaxis explanation content. Contributing to ArviZ can be loosely undestood
  as a single practical craft, but it is more of multiple practical crafts
  coming together.

  This section aims to provide high level overviews explaining different
  contribution types, what skills are needed for them, why are they important
  and guide the reader to their next steps in one of the following sections.

Tutorials
: Diátaxis tutorial content. Tutorials with clearly defined and small scope,
  guiding new contributors in a single "atomic" contribution or even step
  in the contribution process.

How-to guides
: Diátaxis how-to content. Guides that outline the steps for common processes
  in the contributing workflow. Documenting this ensures it is available
  for newcomers but also that all contributors use the same processes
  which can reduce misunderstanding and frictions.

Reference
: Diátaxis reference content.

In depth explanations
: Diátaxis explanation content.

## ArviZ Project website
It is the block formed by the navbar sections
`Community`, `About Us` and by the homepage.

It's audience is anyone interested in the ArviZ project.
This can well be users or contributors but it could also be someone
who is neither. It has also no relation to any practical craft
and does not follow Diátaxis at all.

It should contain all the information related to the ArviZ community,
meeting places (online of physical), shared resources, related
projects and communities; and to the project itself, its members,
its roadmap, its procedures...
# Documentation toolchain

ArviZ documentation is built using a Python documentation tool, [Sphinx](https://www.sphinx-doc.org/en/master/).
Sphinx converts `rst`(restructured text) files into HTML websites.
There are different extensions available for converting other types of files
like markdown, jupyter notebooks, etc. into HTML websites.

ArviZ [docs](https://github.com/arviz-devs/arviz/tree/main/doc/source) consist of `.rst`, `.md` and `.ipynb` files.
It uses `myst-parser` and `myst-nb` for `.md` and `.ipynb` files.
Apart from the content in the `/doc/source` folder, ArviZ documentation also consists of docstrings.
Docstrings are used in the `.py` files to explain the Python objects parameters and return values.

ArviZ docs also uses sphinx extensions for style, layout, navbar and putting code in the documentation.

A more detailed overview of the sphinx toolchain used at ArviZ is available at
{doc}`sphinx-primer:index`
# Running the benchmark suite

To run the **benchmark tests** do the following:

```bash
$ pip install asv
$ cd arviz
$ asv run
```
# Pull request step-by-step

The preferred workflow for contributing to ArviZ is to fork
the [GitHub repository](https://github.com/arviz-devs/arviz/),
clone it to your local machine, and develop on a feature branch.

(pr_steps)=
## Steps

1. Fork the [project repository](https://github.com/arviz-devs/arviz/) by clicking on the 'Fork' button near the top right of the main repository page. This creates a copy of the code under your GitHub user account.

(fork_step_pr)=
2. Clone your fork of the ArviZ repo from your GitHub account to your local disk.

   ::::{tab-set}

   :::{tab-item} SSH
   :sync: ssh

   ```
   $ git clone git@github.com:<your GitHub handle>/arviz.git
   ```
   :::

   :::{tab-item} HTTPS
   :sync: https

   ```
   $ git clone https://github.com/<your GitHub handle>/arviz.git
   ```
   :::

   ::::

3. Navigate to your arviz directory and add the base repository as a remote:

   ::::{tab-set}

   :::{tab-item} SSH
   :sync: ssh

   ```
   $ cd arviz
   $ git remote add upstream git@github.com:arviz-devs/arviz.git
   ```
   :::

   :::{tab-item} HTTPS
   :sync: https

   ```
   $ cd arviz
   $ git remote add upstream https://github.com/arviz-devs/arviz
   ```
   :::

   ::::

(feature_branch_step_pr)=
4. Create a ``feature`` branch to hold your development changes:

   ```bash
   $ git checkout -b my-feature
   ```

   ```{warning}
   Always create a new ``feature`` branch before making any changes. Make your changes
   in the ``feature`` branch. It's good practice to never routinely work on the ``main`` branch of any repository.
   ```

5. Project requirements are in ``requirements.txt``, and libraries used for development are in ``requirements-dev.txt``.  To set up a development environment, you may (probably in a [virtual environment](https://docs.python-guide.org/dev/virtualenvs/)) run:

   ```bash
   $ pip install -r requirements.txt
   $ pip install -r requirements-dev.txt
   $ pip install -r requirements-docs.txt  # to generate docs locally
   ```

   Alternatively, for developing the project in [Docker](https://docs.docker.com/), there is a script to setup the Docker environment for development. See {ref}`developing_in_docker`.

6. Develop the feature on your feature branch. Add your changes using git commands, ``git add`` and then ``git commit``, like:

   ```bash
   $ git add modified_files
   $ git commit -m "commit message here"
   ```

   to record your changes locally.
   After committing, it is a good idea to sync with the base repository in case there have been any changes:
   ```bash
   $ git fetch upstream
   $ git rebase upstream/main
   ```

   Then push the changes to your GitHub account with:

   ```bash
   $ git push -u origin my-feature
   ```

7. Go to the GitHub web page of your fork of the ArviZ repo. Click the 'Pull request' button to send your changes to the project's maintainers for review. This will send an email to the committers.

   :::{tip}
   Now that the PR is ready to submit, check the {ref}`pr_checklist`.
   :::
# Library architecture
ArviZ is organized in modules (the folders in [arviz directory](https://github.com/arviz-devs/arviz/tree/main/arviz)).
The main 3 modules are `data`, `plots` and `stats`.
Then we have 3 more folders. The [tests](https://github.com/arviz-devs/arviz/tree/main/arviz/tests)
folder contains tests for all these 3 modules.

The [static](https://github.com/arviz-devs/arviz/tree/main/arviz/static)
folder is only used to store style and CSS files to get HTML output for `InferenceData`.
Finally we have the [wrappers](https://github.com/arviz-devs/arviz/tree/main/arviz/wrappers)
folder that contains experimental (not tested yet either) features
and interacts closely with both [data](https://github.com/arviz-devs/arviz/tree/main/arviz/data)
and [stats](https://github.com/arviz-devs/arviz/tree/main/arviz/stats) modules.

In addition, there are some files on the higher level directory: `utils.py`, `sel_utils.py`,
`rcparams.py` and `labels.py`.
(review_prs)=
# Review PRs

:::{caution} Page in construction
:::

:::{tip}
Start preparing to be a PR reviewer by looking at how
core contributors approach reviews. It will serve
as a way to see which ways of interacting with contributors work
best while also familiarizing with the developer guide,
coding/documentation conventions and styles
:::

The best place to start is the {ref}`pr_checklist`. Make sure all the
recommendations outlined there are followed.
# Contributing code via pull requests

While issue reporting is valuable, we strongly encourage users who are
inclined to do so to submit patches for new or existing issues via pull
requests. This is particularly the case for simple fixes, such as typos
or tweaks to documentation, which do not require a heavy investment
of time and attention.

Contributors are also encouraged to contribute new code to enhance ArviZ's
functionality, also via pull requests.
Please consult the [ArviZ documentation](https://arviz-devs.github.io/arviz/)
to ensure that any new contribution does not strongly overlap with existing
functionality and open a "Feature Request" issue before starting to work
on the new feature.

## Steps before starting work
Before starting work on a pull request double-check that no one else
is working on the ticket in both issue tickets and pull requests.

ArviZ is a community-driven project and always has multiple people working
on it simultaneously. These guidelines define a set of rules to ensure
that we all make the most of our time and we don't have two contributors
working on the same changes. Let's see what to do when you encounter the following scenarios:

### If an issue ticket exists
If an issue exists check the ticket to ensure no one else has started working on it. If you are first to start work, comment on the ticket to make it evident to others. If the comment looks old or abandoned leave a comment asking if you may start work.

### If an issue ticket doesn't exist
Open an issue ticket for the issue and state that you'll be solving the issue with a pull request. Optionally create a pull request and add `[WIP]` in the title to indicate Work in Progress.

### In the event of a conflict
In the event of two or more people working on the same issue, the general precedence will go to the person who first commented on the issue. If no comments it will go to the first person to submit a PR for review. Each situation will differ though, and the core contributors will make the best judgment call if needed.

### If the issue ticket has someone assigned to it
If the issue is assigned then precedence goes to the assignee. However, if there has been no activity for 2 weeks from the assignment date, the ticket is open for all again and can be unassigned.

(dev_summary)=
## Development process - summary

## Code Formatting
For code generally follow the
[TensorFlow's style guide](https://www.tensorflow.org/community/contribute/code_style)
or the [Google style guide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md).
Both more or less follow PEP 8.

Final formatting is done with [black](https://github.com/ambv/black).

(docstring_formatting)=
# Docstring formatting and type hints

Docstrings should follow the
[numpy docstring guide](https://numpydoc.readthedocs.io/en/latest/format.html).
Extra guidance can also be found in
[pandas docstring guide](https://pandas.pydata.org/pandas-docs/stable/development/contributing_docstring.html).
Please reasonably document any additions or changes to the codebase,
when in doubt, add a docstring.

The different formatting and aim between numpydoc style type description and
[type hints](https://docs.python.org/3/library/typing.html)
should be noted. numpydoc style targets docstrings and aims to be human
readable whereas type hints target function definitions and `.pyi` files and
aim to help third party tools such as type checkers or IDEs. ArviZ does not
require functions to include type hints
however contributions including them are welcome.

## Documentation for user facing methods
If changes are made to a method documented in the {ref}`ArviZ API Guide <api>`
please consider adding inline documentation examples.
You can refer to {func}`az.plot_posterior <arviz.plot_posterior>` for a good example.

### Tests
Section in construction
(docstrings)=
# Docstring style

Docstrings should follow the
[numpy docstring guide](https://numpydoc.readthedocs.io/en/latest/format.html). Read more about it {ref}`docstring_formatting`.
In ArviZ docs, docstrings consit of five main sections, i.e.,
1. Short summary
2. Parameters
3. Returns
4. See Also
5. Examples
Extended summary is strongly encouraged and references is required when relevant in order to cite the papers proposing or explaining the algorithms that are implemented. All other sections can also be used when convenient.

## References
While adding description of parameters, examples, etc, it is important to add references to external libraries and ArviZ.
Docstrings follow the same guide for adding references as the other docs.
For adding references to external libraries functions and objects, see {ref}`reference_external_libs`. For referencing ArviZ objects, follow {ref}`reference_arviz_objects`.

## See Also
In ArviZ docs, we have a lot of interconnected functions both within the library and with external libraries and it can take a lot of time to search for the related functions. It is cruical to add the See Also section to save users time.
For adding the _See Also_ docstring section, you just need to add the function name. Sphinx will
automatically add links to other ArviZ objects and functions listed in the _See Also_
section.

For example, let's add {func}`~arviz.hdi` and {func}`~arviz.plot_ppc` in the _See Also_ section.

```
    See Also
    --------
    hdi : Calculate highest density interval (HDI) of array for given probability.
    plot_ppc : plot for posterior/prior predictive checks.
```

## Kwargs parameters
All the kwargs parameters in {ref}`plots <plot_api>` modules are passed to the matplotlib or bokeh functions. While writing their description, the functions to which they are being passed must be mentioned. In order to check or add those functions, the process is the same for all the kwargs arguments. Let's read the step-by-step guide for `backend_kwargs` as an example, [here](https://github.com/arviz-devs/arviz/wiki/ArviZ-Hacktoberfest-2021).
# Plotting backends
ArviZ supports multiple backends. While adding another backend, please ensure you meet the
following design patterns.

## Code Separation
Each backend should be placed in a different module per the backend.
See `arviz.plots.backends` for examples.

The code in the root level of `arviz.plots` should not contain
any opinion on backend. The idea is that the root level plotting
function performs math and constructs keywords, and the backends
code in `arviz.plots.backends` perform the backend specific
keyword argument defaulting and plot behavior.

The convenience function `get_plotting_function` available in
`arviz.plots.get_plotting_function` should be called to obtain
the correct plotting function from the associated backend. If
adding a new backend follow the pattern provided to programmatically
call the correct backend.

## Test Separation
Tests for each backend should be split into their own module
See [tests.test_plots_matplotlib](https://github.com/arviz-devs/arviz/blob/main/arviz/tests/base_tests/test_plots_matplotlib.py) for an example.

## Gallery Examples
Gallery examples are not required but encouraged. Examples are
compiled into the ArviZ documentation website. The [examples](https://github.com/arviz-devs/arviz/tree/main/examples) directory
can be found in the root of the ArviZ git repository.
---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---

(plots_arguments_guide)=
# Plots' arguments guide

Arviz {ref}`plot <plot_api>` module is used for plotting data. It consists of lot of functions that serves different purposes.
Most of these plotting functions have common arguments. These common arguments are explained below with the examples.

<!--- TODO: use names like centered, non_centered, rugby, radon... -->

```{code-cell}
import arviz as az
data = az.load_arviz_data('centered_eight');
non_centered = az.load_arviz_data('non_centered_eight');
```

:::{warning} Page in construction
:::

(common_var_names)=
## `var_names`

Variables to be plotted, if None all variables are plotted. Prefix the variables by ~ when you want to exclude them from the plot. Let's see the examples.

Plot default autocorrelation

```{code-cell}
az.plot_posterior(data);
```

Plot one variable by setting `var_names=var1`

```{code-cell}
az.plot_posterior(data, var_names=['mu']);
```

Plot subset variables by specifying variable name exactly

```{code-cell}
az.plot_posterior(data, var_names=['mu', 'tau']);
```

Plot variables with regular expressions
```{code-cell}
az.plot_posterior(data, var_names=['mu', '^the'], filter_vars="regex");
```

(common_filter_vars)=
## `filter_vars`
If None (default), interpret `var_names` as the real variables names.

Plot using the default value of `filter_vars` which is `None`

```{code-cell}
az.plot_posterior(data, var_names=['mu']);
```

If “like”, interpret `var_names` as substrings of the real variables names.

Plot using `filter_vars="like"`

```{code-cell}
az.plot_posterior(data, var_names=['mu', '^the'], filter_vars="like");
```

If “regex”, interpret `var_names` as regular expressions on the real variables names. A la `pandas.filter`.

Plot using `filter_vars="regex"`

```{code-cell}
az.plot_posterior(data, var_names=['mu', '^the'], filter_vars="regex");
```

(common_coords)=
## `coords`
Dictionary mapping dimensions to selected coordinates to be plotted. Dimensions without a mapping specified will include all coordinates for that dimension. Defaults to including all coordinates for all dimensions if None.

Using coords argument to plot only a subset of data

```{code-cell}
coords = {"school": ["Choate","Phillips Exeter"]};
az.plot_posterior(data, var_names=["mu", "theta"], coords=coords);
```

(common_hdi_prob)=
## `hdi_prob`

(common_color)=
## `color`

(common_grid)=
## `grid`
Number of rows and columns. Defaults to None, the rows and columns are automatically inferred.

Plot variables in a 4x5 grid

```{code-cell}
az.plot_density([data, non_centered], grid=(4, 5));
```

(common_figsize)=
## `figsize`

`figsize` is short for figure size. If None it will be defined automatically.

(common_textsize)=
## `textsize`

(common_legend)=
## `legend`

(common_ax)=
## `ax`

(common_backend_kwargs)=
## `backend_kwargs`

(common_show)=
## `show`
(wrapper_guide)=
# Sampling wrappers
Sampling wrappers allow ArviZ to call {abbr}`PPL (Probabilistic programming library)`s in order to perform a limited
subset of their capabilities and calculate stats and diagnostics that require
refitting the model on different data.

Their implementation is still experimental and may vary in the future. In fact,
there are currently two possible approaches when creating sampling wrappers.
The first one delegates all calculations to the PPL
whereas the second one externalizes the computation of the pointwise log
likelihood to the user who is expected to write it with xarray/NumPy.

```{toctree}
pystan2_refitting
pystan_refitting
pymc3_refitting
numpyro_refitting
pystan2_refitting_xr_lik
pymc3_refitting_xr_lik
numpyro_refitting_xr_lik
```
# Data structures

```{toctree}
:maxdepth: 2
label_guide
../schema/schema
```

A notebook about working with {class}`~arviz.InferenceData` will be added to
this section soon.
# Computation

```{toctree}
Numba
Dask
```
(plotting_user_guide)=
# Plotting

ArviZ provides multiple plotting functions for different functionalities.
You can browse all the functions in the `plot` module {ref}`by name <plot_api>` or
{ref}`by example output image <example_styles>`.

The pages in this section will cover common functionality between several plotting
functions in tutorial and how-to guide format.

```{toctree}
plots_arguments_guide
```
```{eval-rst}
.. currentmodule:: arviz.labels
```
# Plot utils

(labeller_api)=
## Labellers
See also the {ref}`label_guide`

```{eval-rst}
.. autosummary::
  :toctree: generated/

  BaseLabeller
  DimCoordLabeller
  IdxLabeller
  DimIdxLabeller
  MapLabeller
  NoModelLabeller
  NoVarLabeller
```

## Labeling utils

```{eval-rst}
.. autosummary::
  :toctree: generated/

  mix_labellers
```

## Xarray utils
Low level functions to iterate over xarray objects.

```{eval-rst}
.. currentmodule:: arviz.sel_utils

.. autosummary::
  :toctree: generated/

  xarray_sel_iter
  xarray_var_iter
  xarray_to_ndarray
```
(getting_started)=

# Getting Started

```{toctree}
:maxdepth: 2

Installation
Introduction
XarrayforArviZ
CreatingInferenceData
WorkingWithInferenceData
```
