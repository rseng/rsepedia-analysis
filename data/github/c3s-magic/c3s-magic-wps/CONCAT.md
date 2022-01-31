Contributions are very welcome. Please make sure there is a github issue
associated with with every pull request. Creating an issue is also a good
way to propose new features.This tool allows the user to compute the performance metrics.

    Description of user-changeable settings on webpage':'
    Selection of pressure levels
    Selection of regions (global, tropics 20°N-20°S, northern extratropics 20°-90°N, southern extratropics 20°-90°S)
    Selection of models 
Tool to compute the RMSE between the observed and modelled patterns of variability obtained through classification and their relative relative bias in the percentage of occurrence and the persistence of each mode. Modes of variability are first obtained through EOFs and/or k-mean clustering applied to the anomalies for any season or month, where the user can specify the number of clusters to be computed (e.g. 3 for Arctic sea ice or 4 for North Atlantic sea level pressure). The user also has the option to remove linear or second order trends from the data before obtaining the modes.


    Description of user-changeable settings on webpage':'
    Selection of region.
    Selection of variable.
    Selection of temporal period.
    Selection of RCP scenario.
    Selection of models.
    Selection of detrending method.
    Selection of EOF or spatial data.
EnsClus is a cluster analysis tool in Python, based on the k-means algorithm, for ensembles of climate model simulations. The aim is to group ensemble members according to similar characteristics and to select the most representative member for each cluster.

The user chooses which feature of the data is used to group the ensemble members by the clustering: time mean, maximum, a certain percentile (75% in the examples below), standard deviation and trend over the time period. For each ensemble member this value is computed at each grid point, obtaining N lat-lon maps, where N is the number of ensemble members. The anomaly is computed subtracting the ensemble mean of these maps to each of the single maps. The anomaly is therefore computed with respect to the ensemble members (and not with respect to the time) and the Empirical Orthogonal Function (EOF) analysis is applied to these anomaly maps.

Regarding the EOF analysis, the user can choose either how many Principal Components (PCs) the user wants to retain or the percentage of explained variance the user wants to keep. After reducing dimensionality via EOF analysis, k-means analysis is applied using the desired subset of PCs. The major final outputs are the classification in clusters, i.e. which member belongs to which cluster (in k-means analysis the number k of clusters needs to be defined prior to the analysis) and the most representative member for each cluster, which is the closest member to the cluster centroid.

Other outputs refer to the statistics of clustering: in the PC space, the minimum and the maximum distance between a member in a cluster and the cluster centroid (i.e. the closest and the furthest member), the intra-cluster standard deviation for each cluster (i.e. how much the cluster is compact).

![example output](diagnosticsdata/ensclus/JJApranomaly_hist-min-small.png "Example output")

### Description of user-changeable settings on webpage

1) Selection of models

2) Selection of variable name, as in the NetCDF data file.

3) Selection of total number of ensemble members.

4) Selection of season (DJF, DJFM, NDJFM, JJA so far).

5) Selection of area (EAT, PNA, NH so far).

6) Selection of kind of simulations (historical, scenario).

7) Start year

8) End year

9) Selection of number of clusters.

10) Selection of value to investigate (??th_percentile, mean, maximum, standard deviation, trend).

11) Percentage of variance explained by PCs (perc)

12) Number of PCs (numpcs): if perc is set, the number of PCs is computed according to the required percentage of variance explained by PCs (in the example, the number of PCs that explains al least the 80% of variance is 15: the first 15 PCs explain exactly the 80.73% of variance)

![example output](diagnosticsdata/drought_indicator/drymax.png "Example Output")
**Longest period of consecutive dry days, defined as days with less than 'Dry day limit'.**


![example output](diagnosticsdata/drought_indicator/dryfreq.png "Example Output")
**Number of periods of more than 'Minimum duration limit'.**


The drought indicator currently calculates indicators of meteorological drought by calculating the longest time period of consecutive dry days as well as the number of dry periods of a user defined minimum duration.This python tool is based on the algorithm proposed by [Baldwin and Thompson, 2009], and requires the daily geopotential height field on pressure levels as input. The method is based on an EOF/PC decomposition of the zonally averaged geopotential height, with the leading pattern of variability representative of the (zonal mean) NAM. The calculation is independently repeated at each available pressure level. The daily index can be used to characterize episodic variability of the stratosphere-troposphere connection, while regression on the monthly averaged index is used to quantify the signature of the NAM on the hemispheric climate.

For the diagnostic “Indices of annular modes”, two kind of outputs are produced by the zmnam tool. The first product is the plot of a monthly timeseries of the zonal mean NAM index, from which information on the yearly variability of this pattern can be visualized, and long-term trends can be calculated. When the index is strongly positive, the polar vortex in the stratosphere is concurrently stronger, and positive geopotential height anomalies are found at midlatitudes. The opposite for the other polarity of this index.

![example output](diagnosticsdata/annularmodes/CMIP5_MPI-ESM-MR_amip_r1i1p1_1979-2008_250hPa_mo_ts.png "Example Output")

The monthly averaging smooths fluctuations of the zonal mean index at higher frequencies, which can however be important for the coupling between the troposphere and the stratosphere. For this reason, the daily values of the index are used to build a frequency histogram. From the statistical properties of the resulting probability density function, it is possible to evaluate the realism of a simulated stratospheric variability, by comparing against reanalysis datasets.

![example output](diagnosticsdata/annularmodes/CMIP5_MPI-ESM-MR_amip_r1i1p1_1979-2008_250hPa_da_pdf.png "Example Output")

### Description of user-changeable settings on webpage

1) Selection of model;

2) Selection of period;

3) Selection of supplementary pressure levels

Diagnostic to estimate surge levels along the coast of the North Sea from anomalies in mean sea level pressure and wind components.
![example output](diagnosticsdata/shapefile_selection/Basin_grid2.png "Example Output")

### Description
The user provides a shapefile with one or more polygons. For each polygon, a new timeseries, or CII, is produced with only one time series per polygon. The spatial information is reduced to a representative point ('representative') or as an average of all grid points within the polygon ('mean_inside').     If there are no grid points strictly inside the polygon, 'mean_inside' method defaults to 'representative' for that polygon. Outputs are in the form of a NetCDF file, or as ascii code in csv format.
This namelist implements calculation of the quantile bias to allow evaluation of the precipitation bias based on a user defined quantile in models as compared to a reference dataset. The bias was shown to be typically stronger at higher quantiles (Mehran et al. 2014). The quantile bias (QB) is deﬁned as the ratio of monthly precipitation amounts in each simuation to that of the reference dataset (GPCP observations in the example) above a speciﬁed threshold t (e.g., the 75th percentile of all the local monthly values)

![example output](diagnosticsdata/quantilebias/quantilebias_thumbnail.png "Example Output")

### Description of user-changeable settings on webpage

1) Selection of models

2) Selection of experiment

3) Start year

4) End year

5) Selection of quantile

 This tool computes the diurnal temperature variation indicator as a proxy for energy demand and the wind capacity factor as a proxy for energy supply. The diurnal temperature indicator is the number of days in a season where the daily temperature variation (tasmax - tasmin) exceeds the vulnerability threshold. Here, the vulnerability threshold is based on the mean daily temperature variation for the  reference period (1960 -1990) + 5 degrees. This definition was proposed in a study conducted for the DALKIA energy company. The capacity factor is computed using the manufacturer provided  power  curves that  relate  power  output  to  steady  winds  blowing  at  hub  height.  The capacity  factor can be derived from power curve values easily, by dividing power output by the nominal capacity of the turbine. Here the 100-m wind heights are extrapolated from the surface wind speeds.

    Description of user-changeable settings on webpage':'
    Selection of longitudes and latitudes.
    Selection of projection period.
    Selection of RCP scenario.
    Selection of models.
This tool performs joint timeseries and trend calculation for the HyInt hydroclimatic indices and ETCCDI climate extremes indices over pre-selected continental scale regions or user defined regions. The hydroclimatic indices were selected following Giorgi et al. (2014)':' the simple precipitation intensity index (SDII), the maximum dry spell length (DSL) and wet spell length (WSL), the hydroclimatic intensity index (HY-INT = normalized(SDII) x normalized(DSL)), a measure of the overall behaviour of the hydroclimatic cycle (Giorgi et al., 2011), and the precipitation area (PA), i.e. the area over which at any given day precipitation occurs, (Giorgi et al., 2014). The 27 ETCCDI indices are included in the analysis upon selection from the user. Trends are calculated using the R lm function and significance testing performed with a Student T test on non-null coefficients hypothesis. Trend coefficients are stored together with their statistics which include standard error, t value and Pr(>|t|). The tool produces a variety of types of plots including timeseries with their spread, trend lines and summary plots of trend coefficients.

Figure: output figure type 12 for selected indices and regions calculated for EC-Earth rcp85 over 1976-2100
![example output](diagnosticsdata/hyint_timeseries/hyint_timeseries.png "Example Output")

Figure: output figure type 14 for selected indices and regions calculated for EC-Earth rcp85 over 1976-2100
![example output](diagnosticsdata/hyint_timeseries/hyint_trends1.png "Example Output")
![example output](diagnosticsdata/hyint_timeseries/hyint_trends2.png "Example Output")

### Description of user-changeable settings on webpage:

1)  Selection of model

2)  Selection of period

3)  Selection of reference normalization period to be used for normalized indices

4)  Selection of indices to be plotted from the following list (order-sensitive): "SDII", "DSL", "WSL", "HY-INT", "ABS_INT", "ABS_DSL", "ABS_WSL", "PA", "R95", "altcddETCCDI", "altcsdiETCCDI", "altcwdETCCDI", "altwsdiETCCDI", "cddETCCDI", "csdiETCCDI", "cwdETCCDI", "dtrETCCDI", "fdETCCDI", "gslETCCDI", "idETCCDI", "prcptotETCCDI", "r10mmETCCDI", "r1mmETCCDI", "r20mmETCCDI", "r95pETCCDI", "r99pETCCDI", "rx1dayETCCDI", "rx5dayETCCDI", "sdiiETCCDI", "suETCCDI", "tn10pETCCDI", "tn90pETCCDI", "tnnETCCDI", "tnxETCCDI", "trETCCDI", "tx10pETCCDI", "tx90pETCCDI", "txnETCCDI", "txxETCCDI", "wsdiETCCDI"

5) Select regions for timeseries and maps from the following list: World, World60 (60S/60N), Tropics (30S/30N), South America, Africa, North America, India, Europe, East-Asia, Australia

6) Type of figure: [11] timeseries over required individual region/exp, [12] timeseries over multiple regions/exp, [13] timeseries with multiple models, [14] summary trend coefficients multiple regions, [15] summary trend coefficients multiple models

Tool to compute extreme indices relevant to the insurance industry. These indices are based on the ETCCDI indices, and there are currently 5 available for extreme heati (tx90p), cold (tx10p), wind (wx), drought(cdd) and flooding (rx5day). The individual indices can be combined into a single index with or without weightings for each component. This combined index is roughly analagous to the Actuaries Climate Risk Index.

    Description of user-changeable settings on webpage':'
    Selection of period for defining the baseline thresholds (1960-1990 by default).
    Selection of indice.
    Selection of period for projections.
    Selection of RCP scenario.
    Selection of longitudes and latitudes.
    Selection of models.
This namelist calls RainFARM (https://github.com/jhardenberg/RainFARM.jl). RainFARM is a Julia library and command-line tools implementing the RainFARM stochastic precipitation downscaling method adapted for climate models. The stochastic method is a weather generator which allows to generate fine-scale precipitation fields with a realistic correlation structure, extrapolating to fine scales information  simulated by climate models at regional scales.  RainFARM exploits the nonlinear transformation of a Gaussian random precipitation field obtained extrapolating to small scales the large-scale power spectrum of the fields. It conserves average precipitation at coarse scale (Rebora et al. 2006, D'Onofrio et al. 2014). Description of user-changeable settings on webpage 1) Selection of model; 2) Selection of period; 3) Selection of longitude and latitude.

Figure: original precipitation (mm/day) field (left) and the downscaled field (right) for the EC-Earth model over central Europe.
![example output](diagnosticsdata/rainfarm/RainFARM_example_8x8.png "Example Output")
![example output](diagnosticsdata/rainfarm/RainFARM_example_64x64.png "Example Output")


### Description of user-changeable settings on webpage

1) Selection of model;

2) Selection of period;

3) Selection of longitude and latitude: [lon1,lon2,lat1,lat2], with lon(0/360). Subsetting of region where calculation is to be performed. The selected region needs to have equal number of longitude and latitude grid points;

4) Regridding (preprocessing option): default option false;

5) spatial spectral slope adopted to produce data at local scale (computed from large scales if not specified);

6) number of ensemble members to be calculated;

7) number of subdivisions for downscaling
This functionality is implemented in the namelist MId-Latitude Evaluation System (MiLES), also available as stand-alone package (https://github.com/oloapinivad/MiLES). The tool works on daily 500hPa geopotential height data (with data interpolated on a common 2.5x2.5 grid) and calculates the following diagnostics':'
  -  Z500 Empirical Orthogonal Functions':' Based on CDO "eofs" function. First 4 EOFs for North Atlantic (over the 90W-40E 20N-85N box) and Northern Hemisphere (20N-85N). North Atlantic Oscillation, East Atlantic Pattern, and Arctic Oscillation are thus computed. Figures showing linear regression of PCs on monthly Z500 are provided. PCs and eigenvectors, as well as the variances explained are provided in NetCDF4 Zip format.
  -  North Atlantic Weather Regimes (beta)':' following k-means clustering of 500hPa geopotential height. 4 weather regimes over North Atlantic (80W-40E 30N-87.5N) are evaluted using anomalies from daily seasonal cycle. North Atlantic 4 EOFs are computed to reduce the phase-space dimension and then k-means clustering using Hartigan-Wong algorithm with k=4 is computed. Figures report patterns and frequencies of occurrence. NetCDF4 Zip data are saved. Only 4 regimes and DJF supported so far. 

![example output](diagnosticsdata/teleconnections/EOF1_MPI-ESM-P_r1_1951_2005_DJF.png "NAO EOF1")
![example output](diagnosticsdata/teleconnections/Regime2_MPI-ESM-P_r1_1951_2005_DJF.png "Regimes")

### Description of user-changeable settings on webpage: 

1) Selection of model

2) Selection of reference dataset

3) Selection of period 

4) Selection of season

5) Teleconnection mode (NAO, AO)
This diagnostic contains functions that will provide an integrated view of the historical evolution and projection of the climate system state in terms of the mean state, aggregating the multi-model data into a single signal. This diagnostic will also assess the ensemble variance and the agreement between ensemble members with the mean signal. ensemble mean anomalies for any variable with respect to its mean annual cycle computed on a pre-defined period. ensemble mean variance and agreement over any preselected period and variable of interest.

### Description of user-changeable settings on webpage

- Selection of period for defining the climatology.
- Selection of bias correction method for computing the climatology.
- Selection of either regular, cross-validated or rescaled anomalies.
- Selection of period for projections.
- Selection of longitudes and latitudes.
- Selection of ensemble agreement threshold for stippling.
- Selection of models.
- Selection of temporal averaging (monthly, seasonal, or yearly anomalies)
Tool to combine indices for gridded data from a single or multiple models. The user has the option to include weights, e.g. when combining the components of the insurance risk indice a user may want to give higher weighting to the flood index for a region where flooding is the highest source of insured losses. The user can also compute the area-weighted averages for the indices, or for any other variable, where gridsare assumed to be regular.

    Description of user-changeable settings on webpage':'
    Selection of whether or not to compute area-weighted averages.
    Selection of indices to combine.
    Selection of RCP scenario.
    Selection of longitudes and latitudes.
    Selection of models.
    Selection of weights.

This tool allows the user to compute the single performance metric index.

    Description of user-changeable settings on webpage':'
    Selection of models 
Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna
This functionality is implemented in the namelist MId-Latitude Evaluation System (MiLES), also available as stand-alone package (https://github.com/oloapinivad/MiLES). The tool works on daily 500hPa geopotential height data (with data interpolated on a common 2.5x2.5 grid) and calculates the following diagnostics: 1D Atmospheric Blocking:
  Tibaldi and Molteni (1990) index for Northern Hemisphere. Computed at fixed latitude of 60N, with delta of -5,-2.5,0,2.5,5 deg, fiN=80N and fiS=40N. Full timeseries and climatologies are provided in NetCDF4 Zip format.
  2D Atmospheric blocking: following the index by Davini et al. (2012). It is a 2D version of Tibaldi and Molteni (1990) for Northern Hemisphere atmospheric blocking evaluating meridional gradient reversal at 500hPa. It includes also Meridional Gradient Index and Blocking Intensity index and Rossby wave orientation index, computing both Instantenous Blocking and Blocking Events frequency. Blocking Events allow the estimation of the blocking duration. A supplementary Instantaneous Blocking index with the GHGS2 conditon is also evaluted. Full timeseries and climatologies are provided in NetCDF4 Zip format.

![example output](diagnosticsdata/blocking/TM90_MPI-ESM-P_r1_1951_2005_DJF.png "Example Output")
![example output](diagnosticsdata/blocking/BlockEvents_MPI-ESM-P_r1_1951_2005_DJF.png "Example Output")

### Description of user-changeable settings on webpage: 

1) Selection of model

2) Selection of reference dataset 

3) Selection of period

4) Selection of season
This tool allows the user to compute the total duration of "extreme spells"; the total number of days in a season where the temperature exceeds a threshold over a minimum number of consecutive days. Other tools for computing such events are also available, but the novelty of this tool is that the users can select their own thresholds (based on quantiles) and minimum duration for an event. 

    Description of user-changeable settings on webpage':'
    Selection of longitudes and latitudes.
    Selection of projection period.
    Selection of baseline period.
    Selection of RCP scenario.
    Selection of models.
    Selection of exceedance threshold.
    Selection of minimum duration for an event.
This python tool is based on the algorithm proposed by [Baldwin and Thompson, 2009], and requires the daily geopotential height field on pressure levels as input. The method is based on an EOF/PC decomposition of the zonally averaged geopotential height, with the leading pattern of variability representative of the (zonal mean) NAM. The calculation is independently repeated at each available pressure level. The daily index can be used to characterize episodic variability of the stratosphere-troposphere connection, while regression on the monthly averaged index is used to quantify the signature of the NAM on the hemispheric climate.

To evaluate the modelled strat-trop coupling, the metric is based on the spatial patterns of the zonal mean NAM index. This is obtained by projecting monthly anomalies of the geopotential height field onto the monthly averaged index, then normalized. The well-known annular pattern emerges at upper levels, and it is generally less longitudinally symmetric moving towards the surface.
Having calculated the reanalysis-based spatial patterns, it is possible to compute the difference between these patterns and those reproduced by climate models. The resulting spatial patterns can be used to assess the differences in the strength of this mode of variability and the latitudinal extent.

![example output](diagnosticsdata/stratosphere-troposphere/test250.png "Example Output")

### Description of user-changeable settings on webpage

1) Selection of model;

2) Selection of period;

3) Selection of supplementary pressure levels


Climate explorer allows the user to correlate a time series with other user-defined time series, as well as with system-defined time series of the same temporal resolution. The user is given the option to calculate correlation coefficients using the average of multiple months of climate data. In the example below daily observations are correlated with a timeseries using the correlate function of climate explorer. The result is written as NetCDF in the basket, a link is given which is directly visualizable.
This namelist implements the HyInt tool to calculate a set of 6 indices that allow to evaluate the response of the hydrological cycle to global warming with a joint view of both wet and dry extremes. The indices were selected following Giorgi et al. (2014)':' the simple precipitation intensity index (SDII), the maximum dry spell length (DSL) and wet spell length (WSL), the hydroclimatic intensity index (HY-INT = normalized(SDII) x normalized(DSL)), a measure of the overall behaviour of the hydroclimatic cycle (Giorgi et al., 2011), and the precipitation area (PA), i.e. the area over which at any given day precipitation occurs, (Giorgi et al., 2014). The tool can then produce multi-model and multi-indices maps including the 27 ETCDDI climate extreme indices.

Figure: output figure type 1 for the hyint index calculated for EC-Earth rcp85 multi year mean over 1976-2100 with boxes for user-selected regions

![example output](diagnosticsdata/hyint/hyint_EC-Earth_rcp85_r8i1p1_r320x160_1976_2100_ALL_myear-mean_Globe_map.png "Example Output")

### Description of user-changeable settings on webpage:

1)  Selection of model

2)  Selection of period

3)  Selection of reference normalization period to be used for normalized indices

4)  Selection of indices to be plotted from the following list (order-sensitive): "SDII", "DSL", "WSL", "HY-INT", "ABS_INT", "ABS_DSL", "ABS_WSL", "PA", "R95", "altcddETCCDI", "altcsdiETCCDI", "altcwdETCCDI", "altwsdiETCCDI", "cddETCCDI", "csdiETCCDI", "cwdETCCDI", "dtrETCCDI", "fdETCCDI", "gslETCCDI", "idETCCDI", "prcptotETCCDI", "r10mmETCCDI", "r1mmETCCDI", "r20mmETCCDI", "r95pETCCDI", "r99pETCCDI", "rx1dayETCCDI", "rx5dayETCCDI", "sdiiETCCDI", "suETCCDI", "tn10pETCCDI", "tn90pETCCDI", "tnnETCCDI", "tnxETCCDI", "trETCCDI", "tx10pETCCDI", "tx90pETCCDI", "txnETCCDI", "txxETCCDI", "wsdiETCCDI"

5) Type of figure: 1) lon/lat maps per individual field/exp/multi-year mean, 2) lon/lat maps per individual field exp-ref-diff/multi-year mean, 3) lon/lat maps multi-field/exp-ref-diff/multi-year mean (figure type 3 resembles Giorgi et al. 2011).
CHANGELOG
*********

1.0.0 (2019-11-16)
=====================

* Stable release of c3s-magic-wps. A lot of small fixes were done to all the Processes contained within this Service, and documentation was improved.

1.0.0rc1 (2019-06-06)
=====================

* First release candidate of stable release of c3s-magic-wps. A Web Processing Service for Climate Data Analysis in the MAGIC project. The software in this WPS powers the processes behind the C3S-MAGIC portal. This software was designed to be run on a resource with access to data produced by the CP4CDS project.
C3S MAGIC WPS
=============
.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://opensource.org/licenses/Apache-2.0
    :alt: License

.. image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :target: http://c3s-magic-wps.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.com/c3s-magic/c3s-magic-wps.svg?branch=master
   :target: https://travis-ci.com/c3s-magic/c3s-magic-wps
   :alt: Travis Build

.. image:: https://badge.fury.io/py/c3s-magic-wps.svg
    :target: https://badge.fury.io/py/c3s-magic-wps

.. image:: https://zenodo.org/badge/184254565.svg
   :target: https://zenodo.org/badge/latestdoi/184254565

.. inclusion-marker-start-do-not-remove

Web Processing Service for Climate Data Analysis in the MAGIC project. The software in this WPS powers the processes behind the C3S-MAGIC portal. This software was designed to be run on a resource with access to data produced by the CP4CDS project.

.. inclusion-marker-end-do-not-remove

* Free software: Apache Software License 2.0
* Documentation: https://c3s-magic-wps.readthedocs.io.

Links
-----

* `Climate Data Store`_
* `C3S MAGIC Portal`_
* `CP4CDS Quality Control`_

Credits
-------

This package was created with Cookiecutter_ and the `bird-house/cookiecutter-birdhouse`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`bird-house/cookiecutter-birdhouse`: https://github.com/bird-house/cookiecutter-birdhouse
.. _`Climate Data Store`: https://cds.climate.copernicus.eu
.. _`C3S MAGIC Portal`: http://portal.c3s-magic.eu
.. _`CP4CDS Quality Control`: https://cp4cds-qcapp.ceda.ac.uk
.. _processes:

Processes
=========

.. contents::
    :local:
    :depth: 1


Data Availability
-----------------

A special *Meta* process in the Magic WPS is available to automatically determine the available cmip5 model data by reading the available files from the data folder. As CMIP5 data is structured according to the DRS this can be done automatically. The format used to pass this information is the json output of the `linux tree command`_. If configured correctly, the WPS processes will automatically find the data

.. code-block:: bash

        tree -J -l -d -L 8 /group_workspaces/jasmin2/cp4cds1/data/c3s-cmip5

.. _linux tree command: http://mama.indstate.edu/users/ice/tree/

Metrics
-------

This is a list of metrics available as processes in the MAGIC WPS.

The list of Models, Experiments, and Ensembles is created automatically when the server is started, and shown as "None" here.

.. automodule:: c3s_magic_wps.processes
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
.. _configuration:

Configuration
=============

Command-line options
--------------------

You can overwrite the default `PyWPS`_ configuration by using command-line options.
See the c3s magic wps help which options are available::

    $ c3s_magic_wps start --help
    --hostname HOSTNAME        hostname in PyWPS configuration.
    --port PORT                port in PyWPS configuration.

Start service with different hostname and port::

    $ c3s_magic_wps start --hostname localhost --port 5001

Use a custom configuration file
-------------------------------

You can overwrite the default `PyWPS`_ configuration by providing your own
PyWPS configuration file (just modifiy the options you want to change).
Use one of the existing ``sample-*.cfg`` files as example and copy them to ``etc/custom.cfg``.

For example change the hostname (*demo.org*) and logging level:

.. code-block:: sh

   $ cd c3s_magic_wps
   $ vim etc/custom.cfg
   $ cat etc/custom.cfg
   [server]
   url = http://demo.org:5000/wps
   outputurl = http://demo.org:5000/outputs

   [logging]
   level = DEBUG

   [data]
   archive_root =/cmip5
   obs_root = /obs


Start the service with your custom configuration:

.. code-block:: sh

   # start the service with this configuration
   $ c3s_magic_wps start -c etc/custom.cfg


.. _PyWPS: http://pywps.org/
.. _installation:

Installation
============

.. contents::
    :local:
    :depth: 1

Install from PyPI
------------------

The MAGIC WPS is available as a package in `PyPI <https://pypi.org/>` (c3s-magic-wps). A working installation of `ESMValTool <https://www.esmvaltool.org/>`_ is required, currently version 2.0a2 is supported.

Install from GitHub
-------------------

*Note: These installation instructions assume you have* `Anaconda <https://docs.anaconda.com/anaconda/install/>`_ *installed.*

*Note: these installtion instructions assume you have* `Julia <https://julialang.org/downloads/>`_ *installed*

Check out code from the c3s magic wps GitHub repo and create a conda environment:

.. code-block:: sh

   $ git clone https://github.com/c3s-magic/c3s-magic-wps.git
   $ cd c3s-magic-wps
   $ conda env create -f environment.yml
   $ source activate c3s_magic_wps

*Note: the environment will pull pywps from github via pip, ESMValTool is installed using conda*

Next, to complete the installation of ESMValTool a R and Julia script are required. These are available in the package. Please adjust to match the installation location of conda on your system

.. code-block:: sh

   $ Rscript /opt/conda/envs/c3s-magic-wps/lib/python3.6/site-packages/esmvaltool/install/R/setup.R
   $ julia /opt/conda/envs/c3s-magic-wps/lib/python3.6/site-packages/esmvaltool/install/Julia/setup.jl

Finally install the WPS:

.. code-block:: sh

   $ cd ../c3s-magic-wps
   $ python setup.py develop

Configure c3s magic wps PyWPS service
-------------------------------------

The wps can be configured using the `pywps configuration files <https://pywps.readthedocs.io/en/master/configuration.html>`_. See the etc folder for examples. Create a file called ``.custom.cfg`` to customize settings for your installation. A path to the cmip and obs files is needed to run the metrics.

See the :ref:`installation` section for more info


Start c3s magic wps PyWPS service
---------------------------------

After successful installation you can start the service using the ``c3s_magic_wps`` command-line. An additional environment variable is needed with the location of the model data.

.. code-block:: sh

   $ export CMIP_DATA_ROOT=/path/to/cmip/files


   $ c3s_magic_wps --help # show help
   $ c3s_magic_wps start  # start service with default configuration

   OR

   $ c3s_magic_wps start --daemon # start service as daemon
   loading configuration
   forked process id: 42

*Note: Remember the process ID (PID) so you can stop the service with* ``kill PID``.

The deployed WPS service is by default available on:

http://localhost:5000/wps?service=WPS&version=1.0.0&request=GetCapabilities

You can find which process uses a given port using the following command (here for port 5000):

.. code-block:: sh

   $ netstat -nlp | grep :5000

Check the log files for errors:

.. code-block:: sh

   $ tail -f  pywps.log

Run c3s magic wps as Docker container
-------------------------------------

*Note: These installation instructions assume you have* `Docker <https://docs.docker.com/install/>`_ *installed.*

You can also choose to run c3s magic wps from a Docker container.

Download c3s-magic-wps, build the docker container and run it using docker-compose:  

.. code-block:: sh

   $ git clone https://github.com/c3s-magic/c3s-magic-wps.git
   $ cd c3s-magic-wps
   $ docker-compose build              
   $ docker-compose up

By default the WPS service should be available on port 5000:

 http://localhost:5000/wps?service=wps&request=GetCapabilities

Run docker exec to watch logs:

.. code-block:: sh

   $ docker ps     # find container name
   container_name
   $ docker exec container_name tail -f /opt/wps/pywps.log

Use docker-compose to stop the containers:

.. code-block:: sh

   $ docker-compose down

Use Ansible to deploy c3s magic wps on your System
--------------------------------------------------

Use the `Ansible playbook`_ for PyWPS to deploy c3s magic wps on your system.

.. _Ansible playbook: http://ansible-wps-playbook.readthedocs.io/en/latest/index.html
C3S MAGIC WPS
=============

This is the documentation for the Web Processing Services for Climate Data Analysis developed in the MAGIC project.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   installation
   configuration
   processes
