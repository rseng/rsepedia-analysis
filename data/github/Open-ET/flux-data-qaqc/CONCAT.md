---
title: 'flux-data-qaqc: A Python Package for Energy Balance Closure and Post-Processing of Eddy Flux Data'
tags:
  - Python
  - eddy-covariance
  - energy balance closure
  - post-processing
  - flux tower data
authors:
  - name: John Volk^[corresponding author]
    orcid: 0000-0001-9994-1545
    affiliation: 1
  - name: Justin Huntington
    affiliation: 1
  - name: Richard Allen
    affiliation: 2
  - name: Forrest Melton
    affiliation: "3, 4"
  - name: Martha Anderson
    affiliation: 5
  - name: Ayse Kilic
    affiliation: 6
affiliations:
 - name: Desert Research Institue
   index: 1
 - name: University of Idaho
   index: 2
 - name: NASA Ames Research Center
   index: 3
 - name: California State University Monterey Bay
   index: 4
 - name: USDA Agricultural Research Service
   index: 5
 - name: University of Nebraska Lincoln
   index: 6
date: 27 April 2021
bibliography: paper.bib
---

# Introduction

Eddy Covariance (EC) systems provide continuous direct measurements of turbulent and radiative heat and energy fluxes between the atmosphere and land surface [@baldocchi1988]. Evapotranspiration (ET) estimates from EC systems are highly valuable for many scientific disciplines, water policy decisions, and operational applications such as irrigation management. Although the EC technique is considered accurate under certain atmospheric and spatial conditions, inherent theoretical and practical limitations of the technique are magnified by data uncertainty from post-processing methods. Error can result from the well-documented energy balance closure problem resulting from inherent limitations in the technique but also from data processing errors, instrumentation errors, and other sources [@foken2008; @stoy2013]. Several approaches for post-processing and correcting EC-based ET estimates exist, notably those implemented by major flux measurement networks and datasets like FLUXNET2015 [@pastorello2020]. However, post-processing routines often differ between individual site teams, data networks, and end users, and in some cases post-processing decisions are difficult to reproduce. 

# Statement of need

The ``flux-data-qaqc`` open-source and object-oriented Python package is intended for a wide audience including anyone who needs to process EC flux data to estimate daily or monthly ET that has been corrected for energy balance closure. It is also useful for those who wish to perform visual or flag-based quality control and quality assurance (QA/QC) checks on EC flux or meteorological data to filter out or replace any suspect data points before applying the data for any given use. For example, atmospheric and hydrologic modeling and field scientists may find this software useful for simplifying the post-processing of 30-minute EC flux data into daily ET that has been gap-filled and corrected (or not) for energy imbalance using different approaches. The software will also make calculations of reference ET from EC meteorological data and download reference ET from gridded products [@abatzoglou2013]. Output daily and monthly ET estimates from ``flux-data-qaqc`` can then be used by scientists for model validation and applied science, while knowing that the post-processing steps used to generate them are tested and well documented. 

Post-processing steps of EC data commonly involve at least some of the following steps: data formatting, unit conversions, gap-filling, time aggregation, ET and other atmospheric calculations, energy imbalance assessment and correction, and data QA/QC (automatic or visual). Because ``flux-data-qaqc`` provides a framework for these common EC post-processing steps using Python objects, it also allows scientists or engineers to document and share their processing decisions and analyses with others via the Python scripts they write, which is critical in any serious study. Although the tools provided by ``flux-data-qaqc`` are well documented and "high-level", minimal proficiency in Python is required. 

The interactive  time series and scatter plots (e.g. tools for hover, pan, zoom, selection, and linked axes) produced by ``flux-data-qaqa`` are useful for data QA/QC and for anyone interested in understanding or visually documenting the atmospheric and/or hydrologic conditions at an EC tower. Within the same plot document, daily and monthly surface energy balance variables (average fluxes) as well as all commonly measured meteorological variables (e.g. wind speed, temperature, precipitation, vapor pressure deficit, soil moisture) are displayed. In addition, the plots display computed results useful for EC data QA/QC and many ET-oriented analyses, these results include: energy balance closure (daily and monthly), the energy balance ratio initially and after despiking and smoothing, ET as a fraction of gridded reference ET, as well as daily ET gap information. 

Ultimately, obtaining final ET estimates from EC data involves multiple human-driven processing decisions and this presents an uncertainty issue when comparing data that are processed differently. ``flux-data-qaqc`` does not address issues regarding the deployment configuration  of the EC system or the initial processing of high frequency EC data, for which several software exist. Instead, this software provides a framework and provenance for EC data post-processing (e.g., gap-filling, temporal aggregation), energy balance closure analyses, and ET estimates [@wilson2002].

The [online documentation](https://flux-data-qaqc.readthedocs.io/en/latest/index.html) for ``flux-data-qaqc`` provides more detail regarding these applications, including a comprehensive  tutorial with usage examples. Also included are easy-to-follow installation instructions via the Python Package Index and a virtual environment for automatic third party software dependency management. 

# Design and features

The ``flux-data-qaqc`` package follows a simple hierarchical structure with a few low-level classes that are then inherited in a linear fashion by two classes that are more specific to EC post-processing (\autoref{fig:fig1}). The software is also modular to facilitate incorporation of additional classes and methods. Adding new scientific routines to the package is also facilitated by its standardized and intuitive internal naming schemes and data formats shared by the core classes and their attributes. 

![Diagram showing the relationship between software components, green highlighted boxes indicate software whereas non-highlighted indicate input or output files of the software.\label{fig:fig1}](figure1.pdf)

Input to ``flux-data-qaqc`` is tabular timeseries EC/meterological data and site metadata which are parsed using a configuration file that can be written by the user or using a tool provided with the package, and full [example files](https://flux-data-qaqc.readthedocs.io/en/latest/advanced_config_options.html#example-data-description) are provided. The ``Data`` class parses the data and has tools for unit conversion, quality-based filtering, visualization, and some atmospheric calculations such as dew point temperature and potential solar radiation. From here, the ``QaQc`` class is used for gap-filling, daily and monthly aggregation, energy balance closure correction, writing data using a standard format, and creating interactive diagnostic plots. The ``Data`` and ``QaQc`` classes leverage the Pandas Python library, giving users access to pre- and post-processed data in a Pandas dataframe which facilitates a myriad of other Python-based scientific analyses via Numpy, SciPy and many other packages. The ``Plot`` class contains routines for creating standardized diagnostic plots of pre- and post-processed data, it also provides a robust API for generating interactive (via the [Bokeh package](https://docs.bokeh.org/en/latest/index.html)) scatter and time series plots of arbitrary time series data. The ``Plot`` and ``Convert`` classes along with the ``Util`` module are used or inherited throughout the package; these low-level tools are modular enough to easily be used in other software or custom data pipelines. 

Here are some capabilities and tools provided by ``flux-data-qaqc``:

* ability to read a variety of tabular data formats using a configuration file
* tools for batch processing via writing configuration files programmatically
* automatic unit conversions
* filtering of poor quality data (flag-based and numerical)
* gap-filling of hourly or higher frequency EC data with day and night time treated separately 
* temporal aggregation of high frequency data to daily and monthly
* interactive visual inspection tools at pre-processed data frequency and for daily and monthly post-processed data, for both EC and meteorological data
* management and archiving of pre- and post-processed data and metadata, including standardized naming
* multiple energy balance closure correction routines:  an Energy Balance Ratio method, the Bowen Ratio method, and a multiple linear regression method
* atmospheric calculations, including the hourly and daily forms of the American Society of Civil Engineers standardized reference ET equation [@allen2005] using the [RefET](https://github.com/WSWUP/RefET) python package 
* ability to download gridded weather data and reference ET that corresponds with the EC tower [@abatzoglou2013]
* daily gap-filling of ET using gridded reference ET and site-based fraction of reference ET 

# State of the field 

``flux-data-qaqc`` is most similar to the ONEFlux pipeline [@pastorello2020] in that it includes the same energy balance closure-based corrections with slight modifications. However, there are many differences and improvements for processing of micrometeorological measurements to calculate daily water and energy fluxes. For one, ONEFlux is a command line application as opposed to an object-oriented framework. ONEFlux has many other uses outside of the scope ``flux-data-qaqc``, which is focused on water vapor flux. Similarly, several tools in ``flux-data-qaqc`` are not included in ONEFlux, for example: reading of data in different formats via a configuration file, multiple hourly gap-filling methods, a daily gap-filling method using gridded reference ET, visualization tools, atmospheric calculations such as standardized reference ET, monthly aggregation, and ability for the user to select from multiple alternative energy balance closure correction routines.

Other EC data processing and QA/QC software can be classified into one or more of these groups: 1) processing of high frequency EC data; 2) flux partitioning, e.g. partitioning water vapor flux into transpiration and direct evaporation components; 3) friction velocity (u*) filtering of CO2; 4) post-processed gap-filling and filtering; and 5) flux footprint generation. ``flux-data-qaqc`` is dedicated to tasks in group 4 and has methods and scope unique to other open-source post-processing software known to the authors. Besides ONEFlux, here is a list of all other known and similar EC post-processing software, their purposes and some major differences:

* [MDI Meteo gap-filling tool](http://www.bgc-jena.mpg.de/~MDIwork/meteo/). Online application that performs gap-filling and uses a different methodology.
* [GaFir](https://www.bayceer.uni-bayreuth.de/mm/de/software/software/software_dl.php?id_obj=124194). R package that performs gap-filling and uses a different methodology.
* [ReddyProc](https://www.bgc-jena.mpg.de/bgi/index.php/Services/REddyProcWeb) [@wutzler2018]. R package that performs u* filtering, gap-filling, and flux partitioning. Differences include: CO2 focused, no closure corrections, different gap-filling methods for half-hourly/hourly, no daily gap-filling, no monthly aggregation or plots, strict input formatting, missing ancillary atmospheric calculations. 
* [openeddy](https://github.com/lsigut/openeddy). R package that performs data handling, summarizing, and plotting. Differences include: no closure corrections, different QA/QC procedures, no gap-filling, strict input formatting, missing ancillary atmospheric calculations.
* [PyFluxPro](https://github.com/OzFlux/PyFluxPro). Python GUI application that performs standardized post-processing and QA/QC routines, gap-filling, u* filtering, flux partitioning, and plotting. Differences include: different gap-filling methods, no closure corrections, no daily gap-filling, missing ancillary atmospheric calculations.
* [hesseflux](https://github.com/mcuntz/hesseflux) [@cuntz2020]. Python package similar to ReddyProc with similar differences, namely lack of energy balance closure corrections.

# Research enabled by ``flux-data-qaqc``

This package was designed for and is being used to generate a benchmark ET dataset for ground validation and intercomparison with remote sensing ET models that are part of OpenET [@melton2021]. OpenET is a large collaborative effort to provide satellite-based ET estimates at the field scale across the western  United States, which will greatly improve water management practices among other uses. This package is also being used for processing ET from GRAPEX vineyard EC sites to evaluate spectral-based ET approaches [@kustas2018].

# Acknowledgements

We gratefully acknowledge support for this work from the Walton Family Foundation; Lyda Hill Philanthropies; the S.D. Bechtel, Jr. Foundation; the Gordon and Betty Moore Foundation; the Windward Fund; the Water Funder Initiative; the NASA Applied Science Program and the NASA Western Water Applications Office; the USGS Landsat Science Team; the California State University Agricultural Research Institute; and the Idaho Agricultural Experiment Station and Nebraska Agricultural Experiment Station. 

# References

Change Log
==========

Version 0.1.6
-------------

Add automated tests using GitHb Actions, `see here <https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml>`__ and added in description of how to run tests locally on docs.

Remove ``xlrd`` reader as a dependency due to outdated reading ability as a ``Pandas`` excel reader.

Other minor bug fixes related to ``Plot`` class.

Ass JOSS paper and publish software on Zenodo.

Version 0.1.5
-------------

Add configuration writing function :func:`.util.write_configs` to ``util`` module to facilitate batch processing os similar formatted input files via a station metadata file and data dictionary. 

Update check on energy balance ratio closure correction to also check if the inverse of the energy balance ratio is greater than 0.5, in other words :math:`\frac{1}{EBR} > |0.5|` to avoid closure correction factors that are too small. This check occurs both after step 3 and 6 of the energy balance closure correction routine. 

Version 0.1.4
-------------

Relax default allowance for missing days threshold from 90 (~ 3 days) to 80 % (~ 6 days) in the monthly resample algorithm. In other words if a month has more than 80 % missing daily values, its monthly aggregate will not be resampled, it will be replaced with a null value. The threshold is a keyword argument to the ``util.monthly_resample`` function, but the default is used in any automatic resampling of variables. As a reminder, the number of missing days per month which is tabulated for some variables can be used to fine tune this filter. This change was implemented in version 0.1.4.post1. 

Add daily ASCE standardized reference ET calculation option from the :meth:`.QaQc.daily_ASCE_refET` method. Also added automatic estimation of daily maximum and minimum air temperature from input (e.g. hourly) data and added the input variables to the list of variables that are linearly interpolated before taking daily aggregates in the :obj:`.QaQc` constructor. In other words, the inputs to the daily ASCE reference ET formulation: ea, tmin, tmax, rs, wind speed, are interpolated over daytime and nighttime hourly gaps (2 and 4 default) before taking daily means, mins, maxs, and subsequently used in the daily ASCE calculations. 

Changed default keyword argument ``reference`` to "short" of the :meth:`.Data.hourly_ASCE_refET` method.

Add automatic calculations for high frequency (e.g. hourly or half hourly) data including dew temperature and relative humidity from ea and es if available. The calculations occur when first loading input data, i.e. when :obj:`.Data.df` attribute is accessed. Saturation vapor pressure (es) if calculated at hourly/daily frequency is now saved and added to :obj:`.Data.df` and :obj:`.QaQc.df` properties. 

Require Pandas >= 1.0, changes are not backwards compatible due to internal pandas argument deprecations particularly in the ``pandas.grouper`` object. 

Require Bokeh >= 2.0, changes are not backwards compatible due to legend keyword argument name changes in Bokeh 2.

Minor changes to remove package deprecation warnings from ``Pandas`` and ``Bokeh`` related to their respective large changes. 

Add package dependency ``openpyxl`` package as a fallback for reading in headers of Excel files when ``xlrd`` is unmaintained and failing with previously working tools for reading metadata on Excel files. 

Add a requirements.txt file with package.


Version 0.1.3
-------------

Add option to use gridMET grass reference ET (ETo) and EToF for gap filling daily ET. The default behavior still uses alfalfa reference ET, to use ETo assign the ``refET="ETo"`` keyword argument to :meth:`.QaQc.correct_data` or directly to :meth:`.QaQc._ET_gap_fill`. The ET and ET reference fraction plot labels are updated to show the correct reference ET variable used.

Improve scaling of scatter plots to give equal x and y axis lengths, change return of :meth:`.Plot.scatter_plot` to return tuple of (xmin, xmax, ymin, ymax) for use in plotting one to one lines or limiting axes lengths. 

Version 0.1.2
-------------

Change default functionality of the :meth:`.QaQc.write` method to use the internal variable names (as opposed to the input names) of ``flux-data-qaqc`` in the header files of the output daily and monthly time series CSV files. For example, the column for net radiation is always named and saved as "Rn". This can be reversed to the previous behavior of using the user's input names by setting the new ``use_input_names`` keyword argument to :meth:`.QaQc.write` to ``True``. 

Change the :meth:`.Plot.scatter_plot` underlying call to the ``bokeh`` modules scatter plot as opposed to the set circle glyph plot. This allows the user to change the symbol from circle to others by passing a valid value to the scatter_plot's ``marker`` keyword argument, e.g. ``marker='cross'``.

Version 0.1.1
-------------

Add least squares linear regression method for single or multivariate input; specifically the ``QaQc.lin_regress()`` method. It can be used to correct energy balance components or for any arbitrary time series data loaded in a ``QaQc`` instance. It produces and returns a readable table with regression results (fitted coefficients, root-mean-square-error, etc.) which can be accessed from ``QaQc.lin_regress_results`` after calling the method. The default regression if used to correct energy balance components assumes net radiation is accurate (as the dependent variable):

:math:`Rn = c_0 + c_1 G + c_2 LE + c_3 H`

where :math:`c_0 = 0`.

This regression utilizes the scikit-learn Python module and therefore it was added to the environment and setup files as a dependency.

Version 0.1.0
-------------

Add hourly ASCE standardized reference ET calculation to the ``Data`` class as :meth:`.Data.hourly_ASCE_refET` with options for short and tall (grass and alfalfa) reference ET calculations. If the input data is hourly or higher frequency the input data for the reference ET calculation will automatically be resampled to hourly data. If the input data is hourly then the resulting reference ET time series will be merged with the :attr:`.Data.df` attribute otherwise if the input data is at a temporal frequency > hourly, then the reference ET time series will be return by the :meth:`.Data.hourly_ASCE_refET` method. 

Add methods and options to linearly interpolate energy balance variables based on length of gaps during daytime (:math:`Rn > 0`) and night (:math:`Rn < 0`). These methods are run automatically by the ``QaQc`` constructor if temporal frequency of input is detected as less than daily. New keyword arguments to ``QaQc`` are ``max_interp_hours`` and ``max_interp_hours_night`` respectively.

Other notable changes:

* first release on GitHub
* creation of this file/page (the Change Log)
* add optional return options to plot methods of ``Data`` and ``QaQc`` objects for custimization of default plots or to show/use a subset of them

Version 0.0.9
-------------

Major improvements and notabable changes include:

* add package to PyPI
* change allowable gap percentage for monthly time series to 10 % from 70 %
* add reading of wind direction data, BSD3 license, add package data
* fix bugs related to filtering of subday gaps
* improve plots and other error handling, add feature to hide lines in line plots

Version 0.0.5
-------------

Major improvements and notabable changes include:

* first documentation on `ReadTheDocs <https://flux-data-qaqc.readthedocs.io/en/latest/>`__
* add multiple pages in docs such as installation, config options, basic tutorials, full API reference, etc. 
* improve and streamline config file options
* add vapor pressure and vapor pressure deficit calculations for hourly or lower frequency data in the ``Data.df`` property (upon initial loading of time series into memory
* add automatic unit conversions and checks on select input variables using the ``Convert`` class in the ``util`` module
* add new plots in default plots from ``QaQc`` class, e.g. filtered and raw ETrF
* many rounds of improvements to plots, e.g. hover tooltips, linked axes, style, options for columns, etc. 
* modify Energy Balance Ratio to filter out extreme values of filtered Energy Balance Ratio correction factors
* improve temporal resampling with options to drop days with certain fraction of sub-daily gaps
* track number of gap days in monthly time series of corrected ET 
* add examples of ET gap-filling to docs and change most example data to use Twitchel Island alfalfa site data from AmeriFlux
* add plotting of input data using ``plot`` method of ``Data`` instance which allows for viewing of input data at its initial temporal frequency


Version 0.0.1
-------------

First working version, many changes, milestones included: 

* basic templates and working versions of the ``Data``, ``QaQc``, and ``Plot`` classes 
* versions and improvements to daily and monthly resampling 
* Bowen and Energy Balance Ratio correction routines 
* example Jupyter notebooks including with FLUXNET and USGS data 
* calculation of potential clear sky radiation 
* changing variable naming system to use internal and user names 
* ability to read in multiple soil heat flux and soil moisture measurements and calculate weighted averages 
* make package installable and Conda environment
* add input data filtering using quality control flags (numeric threshold and flags)
* reading of input variables' units
* added the ``util`` submodule with methods for resammpling time series
* ability to take non-weighted averages for any acceptable input variable
* add config file options like date parsing
* removed filtering and smoothing options from Bowen Ratio method and other modifications to it
* add methods for downloading gridMET variables based on location in CONUS
* add routine for gap filling ET based on gridMET ETrF that is smoothed and filtered
* improved ``Plot`` class to contain modular plot methods (line and scatter) for use with arbitrary data
* changed internal variable naming, e.g. etr to ETr
* methods to estimate ET from LE that consider the latent heat of vaporization is affected by air temp.
* other updates to improve code structure and optimization of calculations
.. image:: https://readthedocs.org/projects/flux-data-qaqc/badge/?version=latest
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml/badge.svg
   :target: https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml
   :alt: Automated tests

flux-data-qaqc
================

``flux-data-qaqc`` provides a framework to create reproducible workflows for validation and analysis of eddy covariance data. The package is intended for those who need to post-process flux data, particularly for generating daily and monthly evapotranspiration (ET) timeseries estimates with energy balance closure corrections applied. Applications where this software may be useful include analysis of eddy covariance data, hydrologic or atmospheric model validation, and irrigation and water consumption studies. 

Key functionalities and tools include:

* data validation with methods for quality-based filtering
* time series tools, e.g. gap-filling and temporal aggregation
* energy balance closure algorithms and other meterological calculations
* data provenance, e.g. from metadata management and file structure
* downloading and management of `gridMET <http://www.climatologylab.org/gridmet.html>`__ meterological data
* customizable and interactive visualizations
* built-in unit conversions and batch processing tools

Documentation
-------------

`ReadTheDocs <https://flux-data-qaqc.readthedocs.io/>`_

Installation
------------

Using PIP:

.. code-block:: bash

   pip install fluxdataqaqc

PIP should install the necessary dependencies however it is recommended to use
conda and first install the provided virtual environment. This is useful to
avoid changing your local Python environment. Note, ``flux-data-qaqc`` has been
tested for Python 3.7+, although it may work with versions greater than or
equal to 3.4.

First make sure you have the ``fluxdataqaqc`` environment file, you can download it `here <https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/environment.yml?token=AB3BJKUKL2ELEM7WPLYLXFC45WQOG>`_. Next to install run,

.. code-block:: bash

   conda env create -f environment.yml

To activate the environment before using the ``flux-data-qaqc`` package run,

.. code-block:: bash

   conda activate fluxdataqaqc

Now install using PIP:

.. code-block:: bash

   pip install fluxdataqaqc

Now all package modules and tools should be available in your Python environment PATH and able to be imported. Note if you did not install the Conda virtual environment above, PIP should install dependencies automatically but be sure to be using a version of Python above or equal to 3.4. To test that everything has installed correctly by opening a Python interpretor or IDE and run the following:

.. code-block:: python

   import fluxdataqaqc

and 

.. code-block:: python

   from fluxdataqaqc import Data, QaQc, Plot

If everything has been installed correctly you should get no errors. 


Automated testing with `pytest <https://docs.pytest.org/en/6.2.x/contents.html#>`__
==========================================================================================

Software tests are automatically run each time a change to ``flux-data-qaqc`` is made on the master branch in GitHub using this `GitHub Actions workflow <https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml>`__.  Automated tests help spot potential bugs early so that they be identified and corrected efficiently resulting in an improved user experience.  

Running tests manually
^^^^^^^^^^^^^^^^^^^^^^

``pytest`` is required to run software tests that are provided. You can install ``pytest`` with PIP:

.. code-block:: bash

    pip install pytest

The tests utilize the example flux input data and ``flux-data-qaqc`` configuration files that are provided with the software whether installed from PyPI or GitHub. These files can be found `here <https://github.com/Open-ET/flux-data-qaqc/tree/master/examples>`__.

To run the tests, navivgate to the root directory of the source code (from the command line or shell) and run pytest:

.. code-block:: bash

    pytest
    
This will print out basic test results, usage of ``pytest`` plugins and command line options can be used for getting more information out of the tests.

Contributors
============

OpenET team
^^^^^^^^^^^

Guidance and feedback has been given by the `OpenET <https://etdata.org/>`__
team members,

* Dr. Ayse Kilic
* Christian Dunkerly 
* Forrest Melton
* Dr. Gabriel Senay
* Dr. Joshua Fisher
* Dr. John Volk 
* Dr. Justin Huntington
* Dr. Martha Anderson
* Dr. Richard G. Allen

Acknowledgements
^^^^^^^^^^^^^^^^

The Energy Balance Ratio correction method was modified from the `FLUXNET2015
<https://fluxnet.fluxdata.org/>`__ methodology for correcting daily latent energy and
sensible heat.  Additionally, FLUXNET and `AmeriFlux
<https://ameriflux.lbl.gov/>`__ eddy covariance data was used for code
development and we thank them for mantaining and providing eddy flux data.

We would also like to acknowledge the `RefET <https://github.com/WSWUP/RefET>`__ python package and its developer Charles Morton. ``RefET`` methods are used in ``flux-data-qaqc`` for calculating the American Soceity of Civil Engineers Standardized reference ET using the hourly and daily formulations, as well as the calculation for potential solar radiation.

Contributing
^^^^^^^^^^^^
``flux-data-qaqc`` is an open-source Python package and anyone seriously interested in contributing is encouraged to do so. Look for current issues to get started or offer direct changes with a pull request. 

Report issues or problems with the software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please open a GitHub issue if you have a technical problem with the software. The issues page is `here <https://github.com/Open-ET/flux-data-qaqc/issues>`__. 

Seek support or give feedback
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Any suggestions or feedback that are nontechnical but related to software development or usage please open a `GitHub issue <https://github.com/Open-ET/flux-data-qaqc/issues>`__. Any other feedback should be sent to John.Volk@dri.edu.

.. There are two major functionalities in
   ``flux-data-qaqc``, first, correcting surface energy balance by
   adjusting latent energy and sensible heat fluxes and calculate other
   climatic variables. Second, it serves as a robust way to read in
   different time series data and produce visualizations, e.g. their daily
   and monthly time series.

Configuration Options and Caveats
=================================

This tutorial shows how to use ``flux-data-qaqc`` with climate data of various
formats and generally covers formatting rules of input data and extra options
that can be set in a config file. The major differences when using
``flux-data-qaqc`` for different input data lie in the config file declarations
therefore the entire workflow from :ref:`Tutorial`
will work just the same once your config file is set up correctly. 

Example data description
------------------------

The data used in :ref:`Configuration Options and Caveats`  is provided with
``flux-data-qaqc`` and can be downloaded `here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Config_options>`__.
There are two datasets used in the following examples, the data used for
showing quality based filtering of data based on QC flags is from the FLUXNET
2015 dataset for site "ARM USDA UNL OSU Woodward Switchgrass 1" which contains
switchgrass, more information on this site can be found `here <http://sites.fluxdata.org/US-AR1/>`__. The site that contains multiple soil
heat flux measurements for weighted averaging and soil moisture measurements
for plotting is a subset (shortened for reduced disk space) of the "ARM Southern Great Plains site- Lamont" AmeriFlux site dataset, more information on this site can be found `here <http://ameriflux.lbl.gov/sites/siteinfo/US-ARM>`__.

A reproducible Jupyter Notebook with minor differences of this tutorial can be
found `here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Config_options/advanced_config_options.ipynb>`__.


Setting up a config file
------------------------

``flux-data-qaqc`` starts with the creation of a configuration file, a
text file with extension ".ini" or ".INI" that follows the rules set 
`here <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`__.
A config file for ``flux-data-qaqc`` requires the sections: 1. **METADATA** 2. **DATA**
although you may provide additional for custom uses. 

In **METADATA** you may enter any metadata that you wish so long
as the entrie's *key* is followed by an equal sign and the assigned 
*value*, i.e. 

.. parsed-literal::

    key = value

Here are the *mandatory* metadata entries unique to ``flux-data-qaqc``:

* climate_file_path
* station_elevation
* station_longitude
* station_latitude
* site_id

The “climate_file_path” is the full or relative file path of the input climate
file (excel or CSV, more on formatting this below) containing the climatic data
to be analyzed. The “station_elevation” (meters) and latitude/longitude
(decimal degrees) fields are used to calculate clear sky potential solar
radiation using an ASCE formulation and to locate the nearest gridMET cell
centroid when downloading reference ET. “site_id” is used for saving output
files. 

Optional metadata entries that are used by ``flux-data-qaqc`` include
“missing_data_value”, “qc_threshold”, “qc_flag”, "var_name_delim", "skiprows",
"date_parser", and "gridmet_file_path".  “missing_data_value” is used to
correctly parse missing values in the input climate time series. All other
optional metadata that can be used by ``flux-data-qaqc`` except
"gridmet_file_path" (which is simply the path to a file that is downloaded by
:meth:`.QaQc.download_gridMET`) are explained within this page.

The **DATA** section of the config file is where you specify climate
variables and their units following the same approach explained above. 

Here is a list of all the “expected” climate variable names in the
**DATA** section of a config file, where keys are the keys in the config 
file and values are the internal names used by ``flux-data-qaqc``. This list 
can be accessed at any time from :attr:`.Data.variable_names_dict`:

   >>> from fluxdataqaqc import Data
   >>> Data.variable_names_dict
       {'datestring_col' : 'date' ,
       'net_radiation_col' : 'Rn' ,
       'ground_flux_col' : 'G' ,
       'latent_heat_flux_col' : 'LE' ,
       'latent_heat_flux_corrected_col' : 'LE_user_corr' ,
       'sensible_heat_flux_col' : 'H' ,
       'sensible_heat_flux_corrected_col' : 'H_user_corr' ,
       'shortwave_in_col' : 'sw_in' ,
       'shortwave_out_col' : 'sw_out' ,
       'shortwave_pot_col' : 'sw_pot' ,
       'longwave_in_col' : 'lw_in' ,
       'longwave_out_col' : 'lw_out' ,
       'vap_press_col' : 'vp' ,
       'vap_press_def_col' : 'vpd' ,
       'avg_temp_col' : 't_avg' ,
       'precip_col' : 'ppt' ,
       'wind_spd_col' : 'ws'}

You may view these climate entry keys and values (as found in the config file)
from within Python using the :attr:`.Data.config` property which contains all
information listed in the config file as a :obj:`configparser.ConfigParser`
instance.

    >>> config_path = 'config_for_QC_flag_filtering.ini'
    >>> d = Data(config_path)
    >>> # loop through a list of tuples with keys and values from DATA section
    >>> for each in d.config.items('DATA'):
    >>>     print(each) 
        ('datestring_col', 'date')
        ('net_radiation_col', 'Rn')
        ('net_radiation_units', 'w/m2')
        ('ground_flux_col', 'G')
        ('ground_flux_units', 'w/m2')
        ('latent_heat_flux_col', 'LE')
        ('latent_heat_flux_qc', 'a_qc_value')
        ('latent_heat_flux_units', 'w/m2')
        ('latent_heat_flux_corrected_col', 'LE_corrected')
        ('latent_heat_flux_corrected_units', 'w/m2')
        ('sensible_heat_flux_col', 'H')
        ('sensible_heat_flux_qc', 'a_qc_value')
        ('sensible_heat_flux_units', 'w/m2')
        ('sensible_heat_flux_corrected_col', 'H_corrected')
        ('sensible_heat_flux_corrected_units', 'w/m2')
        ('shortwave_in_col', 'sw_in')
        ('shortwave_in_qc', 'swrad_flag')
        ('shortwave_in_units', 'w/m2')
        ('shortwave_out_col', 'sw_out')
        ('shortwave_out_units', 'w/m2')
        ('shortwave_pot_col', 'sw_pot')
        ('shortwave_pot_units', 'w/m2')
        ('longwave_in_col', 'lw_in')
        ('longwave_in_units', 'w/m2')
        ('longwave_out_col', 'lw_out')
        ('longwave_out_units', 'w/m2')
        ('vap_press_col', 'na')
        ('vap_press_units', 'na')
        ('vap_press_def_col', 'vpd')
        ('vap_press_def_units', 'hPa')
        ('avg_temp_col', 't_avg')
        ('avg_temp_units', 'C')
        ('precip_col', 'ppt')
        ('precip_units', 'mm')
        ('wind_spd_col', 'ws')
        ('wind_spd_units', 'm/s')

You can also access the data from the :attr:`.Data.config` as a dictionary,
for example if your **METADATA** section has an entry for "land_cover", e.g.

.. parsed-literal::
    
    [METADATA]
    land_cover = CROP
    ...

then access this value with :meth:`configparser.ConfigParser.get` which returns 
the value of "land_cover" in the config file's **METADATA** section

   >>> d.config.get('METADATA', 'land_cover')
       CROP

.. tip::
   If you are unsure if your config file's metadata contains a specific entry
   you can pass the ``fallback`` keyword-only argument to the
   :meth:`configparser.ConfigParser.get` method similar to a Python dictionary.

Here is an example,

   >>> d.config.get('METADATA', 'land_cov', fallback='not given')
       "not given"


Input formatting and caveats
----------------------------

Missing data
^^^^^^^^^^^^

For parsing data gaps in input time series assign the
“missing_data_value” to the **METADATA** section of the config file. 
The value should be numeric, e.g.  

.. parsed-literal::

    missing_data_value = -999

If the input time series file does not contain all climate variables that are
expeced by ``flux-data-qaqc``, then specify them as missing (‘na’) in the
config file or simply do not list them in the config. Missing variables will be
ignored for the most part and will not be present in output files/plots,
however if key variables for the energy balance are not present (:math:`LE`,
:math:`H`, :math:`G`, and :math:`Rn`) then you will not be able to run energy
balance closure correction routines.

Data file format
^^^^^^^^^^^^^^^^

``flux-data-qaqc`` accepts Microsoft Excel files (.xlx and .xlsx) and
comma separated value (CSV) text files containing time series input. 
The input file should have a column with combined date and time. Currently there is no
restriction on the temporal frequency of input data however it is
automatically resampled to daily frequency before running correction
routines. Lastly, there should be a single header row containing all
variable names followed by the first entry of climatic variables.

Here is an example of a valid input file’s first 5 rows and 8 columns:

========== ====== ======= ======= ======= ===== === =====
date       t_avg  sw_pot  sw_in   lw_in   vpd   ppt ws
========== ====== ======= ======= ======= ===== === =====
2009-01-01 2.803  186.71  123.108 261.302 1.919 0   3.143
2009-01-02 2.518  187.329 121.842 268.946 0.992 0   2.093
2009-01-03 5.518  188.008 124.241 268.004 2.795 0   4.403
2009-01-04 -3.753 188.742 113.793 246.675 0.892 0   4.336
========== ====== ======= ======= ======= ===== === =====

.. note:: 
   If the the input datas temporal frequency is not recognized
   ``flux-data-qaqc`` will attempt to resample it to daily frequency when it is
   used to create a :obj:`.QaQc` object. Also, if a value is not recognized a
   numeric in any data column it will be forced to a null value.

Data header formatting
^^^^^^^^^^^^^^^^^^^^^^

A common format of some time series data is that the header row may
not start on the first line of the file. If this is the case you must add
an entry to the **METADATA** section of the config file "skiprows" which
states the number of rows to skip before finding the header row. A 
caveat is that if using CSV data files you may have any number of comment
lines before the header so long as they start with a hashtag symbol "#"
(comment), in this case you should not add "skiprows" to **METADATA**. 

Optimize data load time 
^^^^^^^^^^^^^^^^^^^^^^^

``flux-data-qaqc`` utilizes the :mod:`pandas` for most time series data
management, specifically the usage of :obj:`datetime.datetime` objects for
advanced temporal analysis tools. If your file is large you can specify the 
datetime format in the **METADATA** section of the config file to potentially
greatly speedup the loading of data. For example if your date column contains
strings in the format year month day hour minute with no delimiters, e.g. 
201401010000 for 2014 January 1st at midnight, then in the ``flux-data-qaqc``
config file you would enter:

.. parsed-literal::

    date_parser = %Y%m%d%H%M

For more information of the correct date parser string for your date format
see the directives of the :meth:`datetime.datetime.strptime` `here <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior>`__.



--------------

Quality-based data filtering 
----------------------------

Currently ``flux-data-qaqc`` supports filtering out poor quality data
based on user-provided quality control (QC) values (numeric) or flags 
(characters) using the :meth:`.Data.apply_qc_flags` method. This feature 
helps to facilitate manual or semi-manual data filtering which is 
sometimes necessary during data preprocessing.

Flag-based filtering
^^^^^^^^^^^^^^^^^^^^

Let’s say that you have a column in your input data named ‘QC_flag’ that
contains character strings signifying the assigned data quality for a
climate time series. The flag is either ‘g’ meaning a data point is ‘good’
or if the flag is ‘b’ the data point is bad quality and you would like
to filter it. Further let's say that you want the filter to apply to your
latent energy and and sensible heat variables, then in your config file 
you would need to declare the flag for 'bad' data ('b') to be filtered out
in the **METADATA** section:

.. code:: bash

   qc_flag = b

and in the **DATA** section of your config you will state that the
‘QC_flag’ column should be applied to your :math:`LE` and :math:`H` variables:

.. code:: bash

   latent_heat_flux_qc = QC_flag
   sensible_heat_flux_qc = QC_flag

Now, when the :meth:`.Data.apply_qc_flags` method is used the all date entries
of :math:`LE` and :math:`H` that have a "QC_flag" value of 'b' will be forced
to null in the :attr:`.Data.df` property of a :obj:`.Data` instance. 

Threshold-based filtering
^^^^^^^^^^^^^^^^^^^^^^^^^

Another option is to use a numeric quality control *value* that exists
in your input data along with a threshold value which means that when
the quality control value falls below this threshold you would like to
exclude it from the analysis. Let’s assume the column containing the
quality control values is named ‘QC_values’ and it contains values
between 0 and 1 with 0 meaning the poorest quality data and 1 being the
highest and that you would like to remove all data for select variables
with a quality control value below 0.5. Let’s further assume that you
would like this to apply to your incoming solar radiation variable. Then
you would declare the threshold in the **METADATA** section of your
config file:

.. code:: bash

   qc_threshold = 0.5

and in the **DATA** section of your config you will state that the
‘QC_value’ column should be applied to your incoming shortwave radiation
variable:

.. code:: bash

   shortwave_in_qc = QC_value

Now you are all set to use the functionality, note that you may apply
the same quality control value or flag column to multiple climate
variables (as shown in the first example). You may also use both numeric
qualtiy control values and character string flags for the same input
dataset although they cannot both be applied to the same variable. In
other wordsf, if you have a column of quality control numeric values it
cannot also have character strings mixed in. Another option that is used
in the example below is to declare multiple quality control flags that
should be filtered out using a comma separated list. For example in the
provided example config the flags ‘x’ and ‘b’ are used to remove select
days from incoming shorwave radiation,

.. code:: bash

   qc_flag = x, b

There is another option for specifying variables quality control
values/flags. Name the column containing the qualtiy control value/flag
in your input climate file the same as the variable it corresponds to
with the suffix "_QC". For example if your sensible heat column 
was named **sens_h** then your qualtiy control column should be named
**sens_h_QC**. If you use this option you do not need to specify the 
names in your config file. 

Example with flags and thresholds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example uses the provided time series and config files for QC flag filtering on `GitHub <https://github.com/Open-ET/flux-data-qaqc/tree/master/examples/Config_options>`__. The data is from the FLUXNET 2015 site "ARM USDA UNL OSU Woodward Switchgrass 1". This location exhibits switchgrass fields, more information on this site can be found `here <http://ameriflux.lbl.gov/sites/siteinfo/US-AR1>`__.

Because this dataset did not originally contain character Qa/Qc flags they were added for demonstration and applied to shortwave radiation. To view the list of string flags specified in the config file,

    >>> config_path = 'config_for_QC_flag_filtering.ini'
    >>> d = Data(config_path)
    >>> # view or reassign the numeric threshold specified in the config file
    >>> d.qc_threshold
        0.5

And to view the character QC flags if assigned in the config,

    >>> d.qc_flag
        ['x', 'b']

The :attr:`.Data.qc_var_pairs` attribute shows you which variables were found in your input file that have quality control values assigned, it uses the names as found in the input file,

    >>> d.qc_var_pairs
        {'LE': 'a_qc_value', 'H': 'a_qc_value', 'sw_in': 'swrad_flag'}

.. tip::
   The :attr:`.Data.qc_var_pairs` dictionary can be updated in Python to assign
   different columns of QC values to different time series variables.

Now let's apply the QC values.

Note that in this example we mixed both numeric values and threshold
with character flags, the numeric values are being applied to :math:`LE` and :math:`H`
whereas the flags (‘x’ and ‘b’) are applied to incoming shortwave
radiation.

    >>> # make copys of before and after the QC filter is applied
    >>> no_qc = d.df.input_LE.copy()
    >>> no_qc_swrad = d.df.input_sw_in.copy()
    >>> # apply QC flags/values
    >>> d.apply_qc_flags()
    >>> qc_def = d.df.input_LE.copy()
    >>> qc_flag_swrad = d.df.input_sw_in.copy()
        WARNING: renaming column Rn to input_Rn
        WARNING: renaming column G to input_G
        WARNING: renaming column LE to input_LE
        WARNING: renaming column H to input_H
        WARNING: renaming column sw_in to input_sw_in
        WARNING: renaming column sw_out to input_sw_out
        WARNING: renaming column sw_pot to input_sw_pot
        WARNING: renaming column lw_in to input_lw_in
        WARNING: renaming column lw_out to input_lw_out
        WARNING: renaming column vpd to input_vpd
        WARNING: renaming column t_avg to input_t_avg
        WARNING: renaming column ppt to input_ppt
        WARNING: renaming column ws to input_ws


.. note::
   This is a good time to point out that ``flux-data-qaqc`` may change the
   names of your input variables if they exactly match the internal names used
   by the software (see :attr:`.Data.variable_names_dict`, if this is the case
   (as is above) a warning message is printed when reading in the data
   (accessing the ``df`` or ``monthly_df`` properties of :obj:`.Data` or
   :obj:`.QaQc` for the first time) and the names will be modified with a
   prefix of "_input" as shown above.

Here is a plot showing the data before and after applying the filter.

    >>> from bokeh.plotting import ColumnDataSource, figure, show
    >>> from bokeh.models.formatters import DatetimeTickFormatter
    >>> p = figure(x_axis_label='date', y_axis_label='swrad with data removed based on QC value')
    >>> p.line(no_qc_swrad.index, no_qc_swrad, color='red', legend="no flag", line_width=2)
    >>> p.line(no_qc_swrad.index, qc_flag_swrad, color='black', legend="flag = b or x", line_width=2)
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> show(p)


.. raw:: html
    :file: _static/qc_flag1.html

And for :math:`LE`,

    >>> p = figure(x_axis_label='date', y_axis_label='LE with data removed based on QC value')
    >>> p.line(no_qc.index, no_qc, color='red', legend="no QC", line_width=2)
    >>> p.line(no_qc.index, qc_def, color='black', legend="QC=0.5", line_width=2)
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> show(p)


.. raw:: html
    :file: _static/qc_flag2.html


Alternative naming method for QC data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case the climate variables QC columns are named with the same
base name as the climate variables with the ‘\_QC’ suffix. For example
if :math:`LE` is named ‘LE_F_MDS’ in your input files header then the QC column
is named ‘LE_F_MDS_QC’. If your time series data has qualtiy control header names which follow this convention they will automatically be detected and used when you apply them using :meth:`.Data.apply_qc_flags`, i.e. the column names and variables they should be assigned to do not need to be declared in the config.ini file.

    >>> import os
    >>> config_path = os.path.join('..','Basic_usage','fluxnet_config.ini')
    >>> d = Data(config_path)
    >>> # view input files header, note the QC columns 
    >>> d.header
        Index(['TIMESTAMP', 'TA_F', 'TA_F_QC', 'SW_IN_POT', 'SW_IN_F', 'SW_IN_F_QC',
               'LW_IN_F', 'LW_IN_F_QC', 'VPD_F', 'VPD_F_QC', 'PA_F', 'PA_F_QC', 'P_F',
               'P_F_QC', 'WS_F', 'WS_F_QC', 'USTAR', 'USTAR_QC', 'NETRAD', 'NETRAD_QC',
               'PPFD_IN', 'PPFD_IN_QC', 'PPFD_OUT', 'PPFD_OUT_QC', 'SW_OUT',
               'SW_OUT_QC', 'LW_OUT', 'LW_OUT_QC', 'CO2_F_MDS', 'CO2_F_MDS_QC',
               'TS_F_MDS_1', 'TS_F_MDS_1_QC', 'SWC_F_MDS_1', 'SWC_F_MDS_1_QC',
               'G_F_MDS', 'G_F_MDS_QC', 'LE_F_MDS', 'LE_F_MDS_QC', 'LE_CORR',
               'LE_CORR_25', 'LE_CORR_75', 'LE_RANDUNC', 'H_F_MDS', 'H_F_MDS_QC',
               'H_CORR', 'H_CORR_25', 'H_CORR_75', 'H_RANDUNC', 'NEE_VUT_REF',
               'NEE_VUT_REF_QC', 'NEE_VUT_REF_RANDUNC', 'NEE_VUT_25', 'NEE_VUT_50',
               'NEE_VUT_75', 'NEE_VUT_25_QC', 'NEE_VUT_50_QC', 'NEE_VUT_75_QC',
               'RECO_NT_VUT_REF', 'RECO_NT_VUT_25', 'RECO_NT_VUT_50', 'RECO_NT_VUT_75',
               'GPP_NT_VUT_REF', 'GPP_NT_VUT_25', 'GPP_NT_VUT_50', 'GPP_NT_VUT_75',
               'RECO_DT_VUT_REF', 'RECO_DT_VUT_25', 'RECO_DT_VUT_50', 'RECO_DT_VUT_75',
               'GPP_DT_VUT_REF', 'GPP_DT_VUT_25', 'GPP_DT_VUT_50', 'GPP_DT_VUT_75',
               'RECO_SR', 'RECO_SR_N'],
              dtype='object')



Verify that the QC columns have been paired with corresponding climate variables

    >>> d.qc_var_pairs
        {'NETRAD': 'NETRAD_QC',
         'G_F_MDS': 'G_F_MDS_QC',
         'LE_F_MDS': 'LE_F_MDS_QC',
         'H_F_MDS': 'H_F_MDS_QC',
         'SW_IN_F': 'SW_IN_F_QC',
         'SW_OUT': 'SW_OUT_QC',
         'LW_IN_F': 'LW_IN_F_QC',
         'LW_OUT': 'LW_OUT_QC',
         'VPD_F': 'VPD_F_QC',
         'TA_F': 'TA_F_QC',
         'P_F': 'P_F_QC',
         'WS_F': 'WS_F_QC'}


.. note::
   FLUXNET files include their own qualtiy control flags for sensible heat and
   other variables where quality threshold columns are named the same as the
   climate variable they correspond to with the "\_QC" suffix. Therefore they
   do not need to be defined in a config file before applying them. 

For the dataset defined in the example "FLUXNET_config.ini" we did not specify a QC threshold or flag(s) in the config file, therefore we must assign it when calling the :meth:`.Data.apply_qc_flags` method (shown in :ref:`Example of threshold filtering`).

    >>> # view the QC threshold specified in the config file
    >>> print(d.qc_threshold, type(d.qc_threshold))
        None <class 'NoneType'>

Alternatively, you may assign the threshold of flag values at any time directly to a :obj:`.Data` instance:

    >>> d.qc_threshold = .75

Example of threshold filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Be sure to validate QC thresholds or flags before applying them to make sure everything seems correct. Below we see that the lowest QC values correspond with poor quality gap-fill data near the begining of the time series of sensible heat (:math:`H`). 

    >>> from bokeh.models import LinearAxis, Range1d
    >>> p = figure(x_axis_label='date', y_axis_label='sensible heat flux (w/m2)')
    >>> p.extra_y_ranges = {"sec": Range1d(start=-0.1, end=1.1)}
    >>> p.line(d.df.index, d.df['H_F_MDS'], color='red', line_width=1, legend='data')
    >>> p.add_layout(LinearAxis(y_range_name="sec", axis_label='QC value'), 'right')
    >>> p.circle(d.df.index, d.df['H_F_MDS_QC'], line_width=2, y_range_name="sec", legend='QC')
    >>> p.x_range=Range1d(d.df.index[0], d.df.index[365])
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> p.legend.location = "top_left"
    >>> show(p)
        WARNING: Insufficient data to calculate mean for multiple G measurements
        WARNING: Insufficient data to calculate mean for multiple THETA measurements


.. raw:: html
    :file: _static/qc_flag3.html

As a reminder, the routine provided for numeric or theshold filtering removes all data entries that have been assigned to a QC column and have a QC value that falls below some threshold.

    >>> # apply QC numeric threshold filters
    >>> d.apply_qc_flags(threshold=0.5)

Values with QC values < 0.5 are now removed (null) for any variable listed in :attr:`.Data.qc_var_pairs`. 

.. caution::
   The :meth:`.Data.apply_qc_flags()` method applies the filter to all
   variables in the climate file that have a QC column if columns are not
   specified in the config file.

To see all columns (variables) that may have been affected by the previous filter or to constrain them, modify the declarations in the config file or within :attr:`.Data.qc_var_pairs`, i.e.

    >>> d.qc_var_pairs
        {'NETRAD': 'NETRAD_QC',
        'G_F_MDS': 'G_F_MDS_QC',
        'LE_F_MDS': 'LE_F_MDS_QC',
        'H_F_MDS': 'H_F_MDS_QC',
        'SW_IN_F': 'SW_IN_F_QC',
        'SW_OUT': 'SW_OUT_QC',
        'LW_IN_F': 'LW_IN_F_QC',
        'LW_OUT': 'LW_OUT_QC',
        'VPD_F': 'VPD_F_QC',
        'TA_F': 'TA_F_QC',
        'P_F': 'P_F_QC',
        'WS_F': 'WS_F_QC'} 

Now let's view the same sesnible heat flux time series after applying the threshold filter, notice the strange oscillating artifact near the beginning of the time series as been removed:

    >>> p = figure(x_axis_label='date', y_axis_label='sensible heat flux (w/m2)')
    >>> p.extra_y_ranges = {"sec": Range1d(start=-0.1, end=1.1)}
    >>> p.line(d.df.index, d.df['H_F_MDS'], color='red', line_width=1, legend='data')
    >>> p.add_layout(LinearAxis(y_range_name="sec", axis_label='QC value'), 'right')
    >>> p.circle(d.df.index, d.df['H_F_MDS_QC'], line_width=2, y_range_name="sec", legend='QC')
    >>> p.x_range=Range1d(d.df.index[0], d.df.index[365])
    >>> p.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")
    >>> p.legend.location = "top_left"
    >>> show(p)

.. raw:: html
    :file: _static/qc_flag4.html

.. seealso::
   :ref:`Step 0, manual cleaning of poor quality data` for an example that shows   how to filter poor quality data after loading data into a :obj:`.QaQc` 
   object.

--------------

Averaging data from multiple sensors
------------------------------------

Non-weighted averaging
^^^^^^^^^^^^^^^^^^^^^^

If the climate station being analyzed has multiple sensors for the same 
variable (e.g. sensible heat flux) you can easily tell ``flux-data-qaqc``
to use their non-weighted average of for ``flux-data-qaqc`` routines
including energy balance closure corrections or interactive visualizations.
To do so simply list the variable names (as found in the file header) with
a delimiter of your choice and then list the delimiter in the **METADATA**
section. Example, if you have three sensible heat variables named "h_1",
"sens_h_2", and "sensible heat, (w/m2)" then in your config file's 
**METADATA** you would write:

.. parsed-literal::

    var_name_delim = ;

and the sensible heat assignment in the **DATA** section would read:

.. parsed-literal::

    sensible_heat_flux_col = h_1;sens_h_2;sensible heat, (w/m2)

.. caution:: 
   Because there is a comma in the last variable name we cannot use a comma as
   the name delimiter. Also, if you do not state the delimiter of variable
   names in the **METADATA** section of the config file, ``flux-data-qaqc`` will
   look for the single variable name "h_1;sens_h_2;sensible heat, (w/m2)" in
   the header which will not be found.

``flux-data-qaqc`` will name the average in this case as H_mean, in general
it will add the suffix "_mean" to the internal name of the variable used 
by ``flux-data-qaqc`` which can be found in the keys of the :attr:`.Data.variable_names_dict`
dictionary.

.. hint::
   If you use any averaging option for an energy balance component, i.e.
   latent energy, sensible heat, net radiation, or soil heat flux, the average
   will also be used in energy balance closure corrections. 

Weighted averaging
^^^^^^^^^^^^^^^^^^

``flux-data-qaqc`` provides the ability to read in multiple soil heat
flux/moisture variables for a given station location, calculate their
weighted or non weighted average, and write/plot their daily and monthly
time series. *Currently weighted averaging is only provided for 
soil heat flux and soil moisture variables*, using this config option is also
the only way to automatically produce time series plots of these variables
when using :meth:`.QaQc.plot`. This may be useful for comparing/validating multiple soil
heat/moisture probes at varying locations or depths or varying
instrumentation. 

Here is what you need to do to use this functionality:

1. List the multiple soil variable names in your config file's **DATA** section
   following the convention:

-  For multiple soil heat flux variables config names should begin with
   “G\_” or “g\_” followed by an integer starting with 1,2,3,…
   i.e. g_[number]. For example:

.. code:: bash

   g_1 = name_of_my_soil_heat_flux_variable

For soil moisture variables the name of the config variable should follow
“theta_[number]” for example:

.. code:: bash

   theta_1 = name_of_my_soil_moisture_variable

2. List the units of each variable. To specify the units of your soil
   flux/moisture variables add "_units" to the config name you assigned:

.. code:: bash

   g_1_units = w/m2
   theta_1_units = cm

3. To set weights for multiple variables to compute weighted averages
   assign the "_weight" suffix to their names in the config file. For
   example, to set weights for multiple soil heat flux variables:

.. code:: bash

   g_1_weight = 0.25
   g_2_weight = 0.25
   g_3_weight = 0.5

.. hint::
   If weights are not given the arithmetic mean will be calculated. Or if the
   weights do not sum to 1, they will be automatically normalized so that they
   do.

As in the case for non-weighted averaging for any energy balance
component, if you use this option for soil heat flux (:math:`G`), the weighted 
average will also be used in energy balance closure corrections.

Weighted average example
^^^^^^^^^^^^^^^^^^^^^^^^

This example uses time series data recorded from the "ARM Southern Great Plains site- Lamont" AmeriFlux eddy covariance tower, more information on this site can be found `here <http://ameriflux.lbl.gov/sites/siteinfo/US-ARM>`__.

Here is the **DATA** section of the config file that defines the multiple :math:`G` variables in the input data file used for the example below, we put a much higher weight on the :math:`G` sensors "G_2_1_1" and "G_3_1_1",

.. code:: bash

   [DATA]
   g_1 = G_1_1_1
   g_1_units = w/m2
   g_1_weight = 1
   g_2 = G_2_1_1
   g_2_units = w/m2
   g_2_weight = 10
   g_3 = G_3_1_1
   g_3_units = w/m2
   g_3_weight = 10
   g_4 = G_4_1_1
   g_4_units = w/m2
   g_4_weight = 1
   ...

Note, the naming system of these variables (from AmeriFlux conventions) indicates that the multiple :math:`G` sensors are spaced in differing horizontal locations from one another. 

There are many soil moisture sensors at this site, because we are not using these variables within any calculations and simply want them to be loaded in and later plotted we will not assign weights to them and therefore the arithmetic mean will be calculated and added to output plots and time series files. Here is what is listed in the **DATA** section of the config file for multiple soil moisture recordings in this case:

.. code:: bash

   [DATA]
   theta_1 = SWC_1_1_1
   theta_1_units = (%): Soil water content (volumetric), range 0-100
   theta_2 = SWC_2_1_1
   theta_2_units = (%): Soil water content (volumetric), range 0-100
   theta_3 = SWC_1_2_1
   theta_3_units = (%): Soil water content (volumetric), range 0-100
   theta_4 = SWC_2_2_1
   theta_4_units = (%): Soil water content (volumetric), range 0-100
   theta_5 = SWC_3_1_1
   theta_5_units = (%): Soil water content (volumetric), range 0-100
   theta_6 = SWC_4_1_1
   theta_6_units = (%): Soil water content (volumetric), range 0-100
   theta_7 = SWC_3_2_1
   theta_7_units = (%): Soil water content (volumetric), range 0-100
   theta_8 = SWC_4_2_1
   theta_8_units = (%): Soil water content (volumetric), range 0-100
   theta_9 = SWC_1_3_1
   theta_9_units = (%): Soil water content (volumetric), range 0-100
   theta_10 = SWC_1_4_1
   theta_10_units = (%): Soil water content (volumetric), range 0-100
   theta_11 = SWC_1_5_1
   theta_11_units = (%): Soil water content (volumetric), range 0-100
   theta_12 = SWC_1_6_1
   theta_12_units = (%): Soil water content (volumetric), range 0-100
   theta_13 = SWC_2_3_1
   theta_13_units = (%): Soil water content (volumetric), range 0-100
   theta_14 = SWC_2_3_2
   theta_14_units = (%): Soil water content (volumetric), range 0-100
   theta_15 = SWC_2_2_2
   theta_15_units = (%): Soil water content (volumetric), range 0-100
   theta_16 = SWC_2_1_2
   theta_16_units = (%): Soil water content (volumetric), range 0-100
   ...

.. hint:: 
   The units for soil moisture variables will be used in the y-axis daily and
   monthly time series plots when they are created by :meth:`.QaQc.plot`.

Now that the config file has been setup, let's verify that everything was read in correctly,

    >>> # read in the data
    >>> config_path = 'config_for_multiple_soil_vars.ini'
    >>> d = Data(config_path)
    >>> # note the newly added multiple g and theta variables
    >>> d.variables
        {'date': 'TIMESTAMP_START',
        'Rn': 'NETRAD_1_1_1',
        'LE': 'LE_1_1_1',
        'H': 'H_1_1_1',
        'sw_in': 'SW_IN_1_1_1;SW_IN_1_1_2',
        'sw_out': 'SW_OUT_1_1_1',
        'lw_in': 'LW_IN_1_1_1',
        'lw_out': 'LW_OUT_1_1_1',
        'vpd': 'VPD_PI_1_1_1',
        't_avg': 'T_SONIC_1_1_1',
        'ws': 'WS_1_1_1;WS_1_2_1',
        'g_1': 'G_1_1_1',
        'g_2': 'G_2_1_1',
        'g_3': 'G_3_1_1',
        'g_4': 'G_4_1_1',
        'theta_1': 'SWC_1_1_1',
        'theta_2': 'SWC_2_1_1',
        'theta_3': 'SWC_1_2_1',
        'theta_4': 'SWC_2_2_1',
        'theta_5': 'SWC_3_1_1',
        'theta_6': 'SWC_4_1_1',
        'theta_7': 'SWC_3_2_1',
        'theta_8': 'SWC_4_2_1',
        'theta_9': 'SWC_1_3_1',
        'theta_10': 'SWC_1_4_1',
        'theta_11': 'SWC_1_5_1',
        'theta_12': 'SWC_1_6_1',
        'theta_13': 'SWC_2_3_1',
        'theta_14': 'SWC_2_3_2',
        'theta_15': 'SWC_2_2_2',
        'theta_16': 'SWC_2_1_2'}

Note, the windspeed and shortwave incoming radtiation columns were assigned multiple variables as well, these will be used to calculate the non-weighted mean as described in :ref:`Non-weighted averaging`.

Check the units assignment:

    >>> d.units
        {'Rn': 'w/m2',
         'LE': 'w/m2',
         'H': 'w/m2',
         'sw_in': 'w/m2',
         'sw_out': 'w/m2',
         'lw_in': 'w/m2',
         'lw_out': 'w/m2',
         'vpd': 'hPa',
         't_avg': 'C',
         'ws': 'm/s',
         'g_1': 'w/m2',
         'g_2': 'w/m2',
         'g_3': 'w/m2',
         'g_4': 'w/m2',
         'theta_1': '(%): Soil water content (volumetric), range 0-100',
         'theta_2': '(%): Soil water content (volumetric), range 0-100',
         'theta_3': '(%): Soil water content (volumetric), range 0-100',
         'theta_4': '(%): Soil water content (volumetric), range 0-100',
         'theta_5': '(%): Soil water content (volumetric), range 0-100',
         'theta_6': '(%): Soil water content (volumetric), range 0-100',
         'theta_7': '(%): Soil water content (volumetric), range 0-100',
         'theta_8': '(%): Soil water content (volumetric), range 0-100',
         'theta_9': '(%): Soil water content (volumetric), range 0-100',
         'theta_10': '(%): Soil water content (volumetric), range 0-100',
         'theta_11': '(%): Soil water content (volumetric), range 0-100',
         'theta_12': '(%): Soil water content (volumetric), range 0-100',
         'theta_13': '(%): Soil water content (volumetric), range 0-100',
         'theta_14': '(%): Soil water content (volumetric), range 0-100',
         'theta_15': '(%): Soil water content (volumetric), range 0-100',
         'theta_16': '(%): Soil water content (volumetric), range 0-100'}


View these variables and their weights as written in the config file:

    >>> d.soil_var_weight_pairs
        {'g_1': {'name': 'added_G_col', 'weight': '6'},
         'g_2': {'name': 'another_G_var', 'weight': '2'},
         'g_3': {'name': 'G', 'weight': '0.5'},
         'g_4': {'name': 'final_G_var', 'weight': '0.25'},
         'g_5': {'name': 'yet_another_G', 'weight': '0.25'},
         'theta_1': {'name': 'soil_moisture_z1', 'weight': '0.25'},
         'theta_2': {'name': 'soil_moisture_z10', 'weight': '0.75'}}

When the data is first loaded into memory the weighted (and non-weighted) averages are calculated. At this stage weights will be automatically normalized so that they sum to one and the new weights will be printed if this occurs.

    >>> # load daily or monthly dataframe to calculate the weighted averages if they exist
    >>> d.df.head();
        g weights not given or don't sum to one, normalizing
        Here are the new weights:
         G_1_1_1:0.05, G_2_1_1:0.45, G_3_1_1:0.45, G_4_1_1:0.05
        Calculating mean for var: THETA from columns: ['SWC_1_1_1', 'SWC_2_1_1', 'SWC_1_2_1', 'SWC_2_2_1', 'SWC_3_1_1', 'SWC_4_1_1', 'SWC_3_2_1', 'SWC_4_2_1', 'SWC_1_3_1', 'SWC_1_4_1', 'SWC_1_5_1', 'SWC_1_6_1', 'SWC_2_3_1', 'SWC_2_3_2', 'SWC_2_2_2', 'SWC_2_1_2']
        Calculating mean for var: sw_in
         from columns: ['SW_IN_1_1_1', 'SW_IN_1_1_2']
        Calculating mean for var: ws
         from columns: ['WS_1_1_1', 'WS_1_2_1']

In this example, shortwave incoming radiation and windspeed were also averaged (non-weighted) from multiple recordings as described in :ref:`Non-weighted averaging`.

The weights have been changed and updated as we would expect for :math:`G`, you may ignore the weights for soil moisture in this case- because they were not assigned the arithmetic mean is calculated and the weights are not used.

    >>> d.soil_var_weight_pairs
        {'g_1': {'name': 'G_1_1_1', 'weight': 0.045454545454545456},
         'g_2': {'name': 'G_2_1_1', 'weight': 0.45454545454545453},
         'g_3': {'name': 'G_3_1_1', 'weight': 0.45454545454545453},
         'g_4': {'name': 'G_4_1_1', 'weight': 0.045454545454545456},
         'theta_1': {'name': 'SWC_1_1_1', 'weight': 1},
         'theta_2': {'name': 'SWC_2_1_1', 'weight': 1},
         'theta_3': {'name': 'SWC_1_2_1', 'weight': 1},
         'theta_4': {'name': 'SWC_2_2_1', 'weight': 1},
         'theta_5': {'name': 'SWC_3_1_1', 'weight': 1},
         'theta_6': {'name': 'SWC_4_1_1', 'weight': 1},
         'theta_7': {'name': 'SWC_3_2_1', 'weight': 1},
         'theta_8': {'name': 'SWC_4_2_1', 'weight': 1},
         'theta_9': {'name': 'SWC_1_3_1', 'weight': 1},
         'theta_10': {'name': 'SWC_1_4_1', 'weight': 1},
         'theta_11': {'name': 'SWC_1_5_1', 'weight': 1},
         'theta_12': {'name': 'SWC_1_6_1', 'weight': 1},
         'theta_13': {'name': 'SWC_2_3_1', 'weight': 1},
         'theta_14': {'name': 'SWC_2_3_2', 'weight': 1},
         'theta_15': {'name': 'SWC_2_2_2', 'weight': 1},
         'theta_16': {'name': 'SWC_2_1_2', 'weight': 1}}

Now the dataframe also has the weighted means that will be named g_mean and theta_mean,

    >>> d.df.columns
        Index(['input_t_avg', 'input_sw_pot', 'input_sw_in', 'input_lw_in',
               'input_vpd', 'input_ppt', 'input_ws', 'input_Rn', 'input_sw_out',
               'input_lw_out', 'input_G', 'input_LE', 'LE_corrected', 'input_H',
               'H_corrected', 'added_G_col', 'another_G_var', 'final_G_var',
               'yet_another_G', 'soil_moisture_z1', 'soil_moisture_z10', 'a_qc_value',
               'swrad_flag', 'g_mean', 'theta_mean'],
              dtype='object')

.. note:: 
   Even though we did not specify "ground_flux_col" in the config file, the
   weighted average value has now been used to update this variable. Therefore
   the weighted mean will be used in energy balance closure correction routines
   if they are subsequently run.

Check which variable will be used as :math:`G` later if closure corrections are used:

    >>> d.variables.get('G')
        'g_mean'

Now, let's visualize the resulting weighted average of multiple :math:`G` measurements and their individual daily time series,

    >>> # get just G columns for plot arguments
    >>> G_cols = [c for c in d.df.columns if c.startswith(('g_','G_'))]
    >>> G_cols
        ['G_1_1_1', 'G_2_1_1', 'G_3_1_1', 'G_4_1_1', 'g_mean']

The example below creates the time series plot with a short span of data for easier visibility of weighted mean, it also used the plot routines provided by :obj:`.Data` and :obj:`.QaQc` which are inhereted from the :obj:`.Plot` class within ``flux-data-qaqc``. Specifically this example utilizes :meth:`.Plot.add_lines` which makes the time series plotting of multiple variables more efficient and automatically handles the hover tooltips.

    >>> from fluxdataqaqc.plot import ColumnDataSource # for hover tooltips
    >>> # shorter period for visualization
    >>> df = d.df.loc['01/01/2008':'05/01/2008', G_cols]
    >>> plt_vars = G_cols
    >>> colors = ['blue', 'red', 'orange', 'green', 'black']
    >>> x_name = 'date'
    >>> source = ColumnDataSource(df)
    >>> fig = figure(x_axis_label='date', y_axis_label='Soil heat flux (w/m2)')
    >>> Data.add_lines(fig, df, plt_vars, colors, x_name, source, labels=G_cols)
    >>> show(fig)

.. raw:: html
    :file: _static/weighted_g.html

Note, the weighted mean is closer to 'G_2_1_1' and 'G_3_1_1' as we gave them weights of 10 versus 1 to 'G_1_1_1' and 'G_4_1_1'.

Lastly, the code snippets below run the Energy Balance Ratio closure correction and creating the default plots in order to view the daily and monthly time series of multiple soil moisture variables. It also shows how to upload the output plot file into a Jupyter Notebook for viewing.

    >>> # in order to correctly view the output in a Jupyter notebook
    >>> from bokeh.io import output_notebook
    >>> output_notebook()

Within the set of default plots created by the :meth:`.QaQc.plot` method will include interactive daily and monthly time series of multiple :math:`G` and soil moisture variables if they were assigned in the input config file (as in this example), scroll down to view them. 

    >>> from fluxdataqaqc import QaQc
    >>> q = QaQc(d)
    >>> q.correct_data()
    >>> # this will NOT save the plot file, use output_type='save'
    >>> q.plot(output_type='show')


.. raw:: html
    :file: _static/US-ARM_multipe_soilvars_plots.html
 

Tutorial
========

This tutorial demonstrates the most important features of the
``flux-data-qaqc`` Python package for management, analysis, and visualization
of eddy covariance time series data. It is recommended to read the
:ref:`Installation` and :ref:`Configuration Options and Caveats` tutorials before this one. 

A Jupyter Notebook of this tutorial is available `here <https://github.com/Open-ET/flux-data-qaqc/blob/master/examples/Basic_usage/Tutorial.ipynb>`__.

.. Tip:: 
   Currently, the software does not include a command line interface therefore
   to use the software you must use Python, e.g. make your own scripts or use
   an interactive shell. However, you will see that common workflows can be
   accomplished with a few (5-10) lines of code and you can simply follow the
   templates given here to make custom scripts.

Description of example datasets
-------------------------------

The data for this example comes from the “Twitchell Alfalfa” AmeriFlux
eddy covariance flux tower site in California. The site is located in
alfalfa fields and exhibits a mild Mediterranean climate with dry and
hot summers, for more information on this site or to download data click
`here <https://ameriflux.lbl.gov/sites/siteinfo/US-Tw3>`__.

Loading input
-------------

The loading and management of input climatic data and metadata from a
config.ini file is done using the :obj:`fluxdataqaqc.Data` object. In a
nutshell, a :obj:`.Data` object is created from a properly formatted config file
(see :ref:`Setting up a config file`) and has tools for parsing input climate
data, averaging input climate time series, accessing/managing metadata,
flag-based data filtering, and creating interactive visualizations of input
data.

There is only one argument to create a Data object, the path to the
config.ini file:

    >>> # imports for code snippets within tutorial
    >>> import pandas as pd
    >>> from fluxdataqaqc import Data, QaQc, Plot
    >>> from bokeh.plotting import figure, show, ColumnDataSource
    >>> from bokeh.models.formatters import DatetimeTickFormatter
    >>> from bokeh.models import LinearAxis, Range1d
    >>> # create a Data object from the config.ini file    
    >>> config_path = 'US-Tw3_config.ini'
    >>> d = Data(config_path)

Attributes of a Data object
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below are some of the useful attributes of the :obj:`.Data` object and how
they may be used.

The full path to the config.ini file that was used to create the
:obj:`.Data` instance can be accessed, note that it will return a
system-depenedent :obj:`pathlib.Path` object. E.g. on my Linux machine the
path is:

    >>> d.config_file
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/US-Tw3_config.ini')



On a Windows machine the path will have the appropriate backslashes.

Similarly to access the climate time series file:

    >>> d.climate_file
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/AMF_US-Tw3_BASE_HH_5-5.csv')



The :attr:`.Data.config` attribute is a :obj:`configparser.ConfigParser` object,
it allows you to access metadata and data in the config file in multiple
ways and to modify them. In ``flux-data-qaqc`` it is mainly used for
accessing information about the input data.

    >>> # get a list of all entries in the METADATA section of the config.ini
    >>> d.config.items('METADATA') # access the DATA section the same way
        [('climate_file_path', 'AMF_US-Tw3_BASE_HH_5-5.csv'),
         ('station_latitude', '38.1159'),
         ('station_longitude', '-121.6467'),
         ('station_elevation', '-9.0'),
         ('missing_data_value', '-9999'),
         ('skiprows', '2'),
         ('date_parser', '%Y%m%d%H%M'),
         ('site_id', 'US-Tw3'),
         ('country', 'USA'),
         ('doi_contributor_name', 'Dennis Baldocchi'),
         ('doi_contributor_role', 'Author'),
         ('doi_contributor_email', 'baldocchi@berkeley.edu'),
         ('doi_contributor_institution', 'University of California, Berkeley'),
         ('doi_organization', 'California Department of Water Resources'),
         ('doi_organization_role', 'Sponsor'),
         ('flux_measurements_method', 'Eddy Covariance'),
         ('flux_measurements_variable', 'CO2'),
         ('flux_measurements_operations', 'Continuous operation'),
         ('site_name', 'Twitchell Alfalfa'),
         ('igbp', 'CRO'),
         ('igbp_comment',
          'alfalfa is a fast growing leguminous crop raised for animal feed of low stature.  It is planted in rows and typically reaches 60-70 cm in height prior to harvest.'),
         ('land_ownership', 'public'),
         ('network', 'AmeriFlux'),
         ('reference_paper',
          'Baldocchi, D., Penuelas, J. (2018) The Physics And Ecology Of Mining Carbon Dioxide From The Atmosphere By Ecosystems, Global Change Biology, 45(), 9275–9287'),
         ('reference_doi', '10.1111/gcb.14559'),
         ('reference_usage', 'Reference'),
         ('research_topic',
          'The research approach of the University of California, Berkeley Biometeorology Laboratory involves the coordinated use of experimental measurements and theoretical models to understand the physical, biological, and chemical processes that control trace gas fluxes between the biosphere and atmosphere and to quantify their temporal and spatial variations. The research objectives of the Mayberry Wetland, Twitchell Wetland, Sherman Island, Twitchell Island, Twitchell Alfalfa,  and Twitchell Corn sites are as follows: 1) Describe differences in the fluxes of CO2, CH4, H2O, and energy between different land uses, 2) Understand the mechanisms controlling these fluxes, 3) Use ecosystem modeling to understand controls on these mechanisms under different environmental scenarios. These six sites were selected to capture a wide range of inundated conditions within the Sacramento-San Joaquin River Delta. The research focuses on the eddy covariance technique to measure CH4, CO2, H2O, and energy fluxes and works to combine measurements of both net fluxes and partitioned fluxes in order to achieve a mechanistic understanding of the ecological controls on current and future carbon flux in the Delta.'),
         ('terrain', 'Flat'),
         ('aspect', 'FLAT'),
         ('wind_direction', 'W'),
         ('surface_homogeneity', '370.0'),
         ('site_desc',
          "The Twitchell Alfalfa site is an alfalfa field owned by the state of California and leased to third parties for farming. The tower was installed on May 24, 2013. This site and the surrounding region are part of the San Joaquin - Sacramento River Delta drained beginning in the 1850's and subsequently used for agriculture. The field has been alfalfa for X years…., Crop rotation occurs every 5-6 years.  The site is harvested by mowing and bailing several times per year.  The field is fallow typically between November and February. The site is irrigated by periodically-flooded ditches surrounding the field. The site is irrigated by raising, and subsequently lowering the water table??"),
         ('site_funding', 'California Department of Water Resources'),
         ('team_member_name', 'Joe Verfaillie'),
         ('team_member_role', 'Technician'),
         ('team_member_email', 'jverfail@berkeley.edu'),
         ('team_member_institution', 'University of California, Berkeley'),
         ('url_ameriflux', 'http://ameriflux.lbl.gov/sites/siteinfo/US-Tw3'),
         ('utc_offset', '-8'),
         ('mat', '15.6'),
         ('map', '421.0'),
         ('land_owner', 'California Department of Water Resources'),
         ('climate_koeppen', 'Csa'),
         ('doi', '10.17190/AMF/1246149'),
         ('doi_citation',
          'Dennis Baldocchi (2013-) AmeriFlux US-Tw3 Twitchell Alfalfa, 10.17190/AMF/1246149'),
         ('doi_dataproduct', 'AmeriFlux'),
         ('team_member_address',
          'Department of Environmental Science, Policy and Management, 137 Mulford Hall, 345 Hilgard Hall,Berkeley, CA USA 94720-3110'),
         ('url', 'http://nature.berkeley.edu/biometlab/sites.php?tab=US-Tw3'),
         ('dom_dist_mgmt', 'Agriculture'),
         ('site_snow_cover_days', '0.0'),
         ('state', 'CA'),
         ('location_date_start', '20130524.0'),
         ('acknowledgement',
          'Biometeorology Lab, University of California, Berkeley, PI:  Dennis Baldocchi')]


A useful method is the :meth:`configparser.ConfigParser.get` which takes the
section of the config file and the “option” and returns the value:

    >>> d.config.get(section='METADATA', option='site_name')
        'Twitchell Alfalfa'


    >>> # section and option are optional keywords
    >>> d.config.get('METADATA', 'site_name')
        'Twitchell Alfalfa'

.. Tip::
   If you are unsure if an entry or option exists in the config file, use the
   ``fallback`` keyword argument

    >>> # section and option are optional keywords
    >>> d.config.get('METADATA', 'site name', fallback='na')
        'na'

Some metadata entries are added as :obj:`.Data` attributes for easier access
as they are used in multiple ways later, these include:

-  site_id\ :math:`^*`
-  elevation\ :math:`^*`
-  latitude\ :math:`^*`
-  longitude\ :math:`^*`
-  na_val
-  qc_threshold
-  qc_flag

:math:`^*`\ mandatory **METADATA** entries in the config file, see
:ref:`Setting up a Config File` for further explanation.

View all the columns as found in the header row of the input time series
climate file.

    >>> d.header
        array(['TIMESTAMP_START', 'TIMESTAMP_END', 'CO2', 'H2O', 'CH4', 'FC',
               'FCH4', 'FC_SSITC_TEST', 'FCH4_SSITC_TEST', 'G', 'H', 'LE',
               'H_SSITC_TEST', 'LE_SSITC_TEST', 'WD', 'WS', 'USTAR', 'ZL', 'TAU',
               'MO_LENGTH', 'V_SIGMA', 'W_SIGMA', 'TAU_SSITC_TEST', 'PA', 'RH',
               'TA', 'VPD_PI', 'T_SONIC', 'T_SONIC_SIGMA', 'SWC_1_1_1',
               'SWC_1_2_1', 'TS_1_1_1', 'TS_1_2_1', 'TS_1_3_1', 'TS_1_4_1',
               'TS_1_5_1', 'NETRAD', 'PPFD_DIF', 'PPFD_IN', 'PPFD_OUT', 'SW_IN',
               'SW_OUT', 'LW_IN', 'LW_OUT', 'P', 'FC_PI_F', 'RECO_PI_F',
               'GPP_PI_F', 'H_PI_F', 'LE_PI_F'], dtype='<U15')


.. Note::
   All of the header columns will not necessarily be loaded, only those
   specified in the config file. Also, no data other than the header line is
   loaded into memory when creating a :obj:`.Data` object, the time series data
   is only loaded when calling :attr:`.Data.df` for increased efficiency for
   some workflows involving only metadata.

Variable names and units
^^^^^^^^^^^^^^^^^^^^^^^^

In ``flux-data-qaqc`` there are two naming schemes for climate
variables, the names as defined by the column headers in the input time
series file and the internal names for some variables and calculated
variables created by the package. We will refer to these two sets as
“user-defined” and “internal” names hereforth.

The :attr:`.Data.variables` attribute maps the internal to user-defined
variable names:

    >>> d.variables
        {'date': 'TIMESTAMP_START',
         'Rn': 'NETRAD',
         'G': 'G',
         'LE': 'LE_PI_F',
         'H': 'H_PI_F',
         'sw_in': 'SW_IN',
         'sw_out': 'SW_OUT',
         'lw_in': 'LW_IN',
         'lw_out': 'LW_OUT',
         'vpd': 'VPD_PI',
         't_avg': 'T_SONIC',
         'ws': 'WS',
         'theta_1': 'SWC_1_1_1',
         'theta_2': 'SWC_1_2_1'}


And, the :attr:`.Data.inv_map` maps the internal to user-defined names if
they differ, however this is only created once the data is loaded by
calling :attr:`.Data.df`.

    >>> # a similar dictionary attribute for input units
    >>> d.units
        {'Rn': 'w/m2',
         'G': 'w/m2',
         'LE': 'w/m2',
         'H': 'w/m2',
         'sw_in': 'w/m2',
         'sw_out': 'w/m2',
         'lw_in': 'w/m2',
         'lw_out': 'w/m2',
         'vpd': 'hPa',
         't_avg': 'C',
         'ws': 'm/s',
         'theta_1': '(%): Soil water content (volumetric), range 0-100',
         'theta_2': '(%): Soil water content (volumetric), range 0-100'}


Accessing input data
^^^^^^^^^^^^^^^^^^^^

The :py:attr:`Data.df` property gves access to the time series input climate
data for columns specified in the config file as a datetime-indexed
:obj:`pandas.DataFrame` object. This object has numerous powerful built in
tools for time series analysis and visualization.

    >>> # first 5 datetimes that are not gaps
    >>> d.df.dropna().head()

.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>input_G</th>
          <th>WS</th>
          <th>VPD_PI</th>
          <th>T_SONIC</th>
          <th>SWC_1_1_1</th>
          <th>SWC_1_2_1</th>
          <th>NETRAD</th>
          <th>SW_IN</th>
          <th>SW_OUT</th>
          <th>LW_IN</th>
          <th>LW_OUT</th>
          <th>H_PI_F</th>
          <th>LE_PI_F</th>
          <th>theta_mean</th>
        </tr>
        <tr>
          <th>date</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>2013-05-24 12:30:00</th>
          <td>122.194848</td>
          <td>3.352754</td>
          <td>18.853678</td>
          <td>25.682739</td>
          <td>6.6790</td>
          <td>26.1655</td>
          <td>652.648719</td>
          <td>1027.756939</td>
          <td>212.800000</td>
          <td>300.363524</td>
          <td>462.671744</td>
          <td>95.487930</td>
          <td>375.841436</td>
          <td>16.42225</td>
        </tr>
        <tr>
          <th>2013-05-24 13:00:00</th>
          <td>108.054863</td>
          <td>3.882154</td>
          <td>18.560999</td>
          <td>26.057700</td>
          <td>6.7065</td>
          <td>26.1600</td>
          <td>629.990486</td>
          <td>997.749437</td>
          <td>209.933333</td>
          <td>303.269447</td>
          <td>461.095065</td>
          <td>96.584383</td>
          <td>371.619775</td>
          <td>16.43325</td>
        </tr>
        <tr>
          <th>2013-05-24 13:30:00</th>
          <td>79.330662</td>
          <td>4.646089</td>
          <td>18.900260</td>
          <td>26.067374</td>
          <td>6.7120</td>
          <td>26.1545</td>
          <td>595.817687</td>
          <td>954.988747</td>
          <td>206.733333</td>
          <td>303.017852</td>
          <td>455.455579</td>
          <td>84.066406</td>
          <td>358.194935</td>
          <td>16.43325</td>
        </tr>
        <tr>
          <th>2013-05-24 14:00:00</th>
          <td>52.366527</td>
          <td>5.048825</td>
          <td>20.440061</td>
          <td>25.961307</td>
          <td>6.7395</td>
          <td>26.1325</td>
          <td>549.039365</td>
          <td>900.975244</td>
          <td>201.333333</td>
          <td>298.914731</td>
          <td>449.517276</td>
          <td>69.449710</td>
          <td>406.528564</td>
          <td>16.43600</td>
        </tr>
        <tr>
          <th>2013-05-24 14:30:00</th>
          <td>35.658417</td>
          <td>5.302946</td>
          <td>21.064824</td>
          <td>25.954462</td>
          <td>6.7450</td>
          <td>26.1215</td>
          <td>493.519695</td>
          <td>833.458365</td>
          <td>192.066667</td>
          <td>296.791541</td>
          <td>444.663544</td>
          <td>47.774030</td>
          <td>315.295309</td>
          <td>16.43325</td>
        </tr>
      </tbody>
    </table>
    </div>
    <br />



.. Tip:: 
   There are *many* tutorials on how to use the :obj:`pandas.DataFrame` and its
   powerful data analysis tools for multiple purposes online, to get started
   you may want to visit Panda’s own list of tutorials `here
   <https://pandas.pydata.org/pandas-docs/stable/getting_started/tutorials.html#internal-guides>`__.

By default the column names in :attr:`.Data.df` are retained from
user-defined names unless they were named exactly the same as an
internal name. For example the input ground heat flux column in this
dataset is named “G”, therefore it was renamed as “input_g”

    >>> d.df.columns
        Index(['input_G', 'WS', 'VPD_PI', 'T_SONIC', 'SWC_1_1_1', 'SWC_1_2_1',
               'NETRAD', 'SW_IN', 'SW_OUT', 'LW_IN', 'LW_OUT', 'H_PI_F', 'LE_PI_F',
               'theta_mean'],
              dtype='object')


    >>> # the new name was also updated in Data.variables
    >>> d.variables.get('G')
        'input_G'



As stated earlier, :attr:`.Data.inv_map` maps the user-defined names to
internal ``flux-data-qaqc`` names only after loading :attr:`.Data.df`:

    >>> d.inv_map
        {'TIMESTAMP_START': 'date',
         'NETRAD': 'Rn',
         'input_G': 'G',
         'LE_PI_F': 'LE',
         'H_PI_F': 'H',
         'SW_IN': 'sw_in',
         'SW_OUT': 'sw_out',
         'LW_IN': 'lw_in',
         'LW_OUT': 'lw_out',
         'VPD_PI': 'vpd',
         'T_SONIC': 't_avg',
         'WS': 'ws',
         'SWC_1_1_1': 'theta_1',
         'SWC_1_2_1': 'theta_2'}



.. Tip:: 
   The :attr:`.Data.inv_map` is mainly used to rename the dataframe to internal
   names, this can be very useful if you are creating your own custom workflows
   using the ``flux-data-qaqc`` API because it allows you to only know the
   internal names of variables therefore they can be hard coded into your
   workflow and applied to different eddy covariance datasets. For example,
   let’s say we wanted to make HTML tables of basic statistics of just the
   energy balance components for many datasets (that may have different names
   for the same variables) and save the file using the user-defined names:

   >>> d = Data('US-Tw3_config.ini')
   >>> df = d.df.rename(columns=d.inv_map)
   >>> # get some metadata for saving
   >>> site_id = d.site_id
   >>> vars_we_want = ['H', 'LE', 'Rn', 'G']
   >>> # rename variables, calculate basice statistics table and save to HTML
   >>> df[vars_we_want].rename(columns=d.variables).describe().to_html('{}.html'.format(site_id))
       Calculating mean for var: THETA from columns: ['SWC_1_1_1', 'SWC_1_2_1']
       WARNING: renaming column G to input_G
   >>> # which produces the following HTML table with user-defined names:
   >>> from IPython.display import HTML
   >>> HTML(filename='{}.html'.format(site_id))


.. raw:: html
    :file: _static/tutorial/US-Tw3.html

.. raw:: html

       <br />

Another powerful feature of the :attr:`.Data.df` property is that it is
datetime-indexed using the input data’s temporal frequency, view the
date index like so:

    >>> d.df.index
        DatetimeIndex(['2013-01-01 00:00:00', '2013-01-01 00:30:00',
                       '2013-01-01 01:00:00', '2013-01-01 01:30:00',
                       '2013-01-01 02:00:00', '2013-01-01 02:30:00',
                       '2013-01-01 03:00:00', '2013-01-01 03:30:00',
                       '2013-01-01 04:00:00', '2013-01-01 04:30:00',
                       ...
                       '2018-06-04 19:00:00', '2018-06-04 19:30:00',
                       '2018-06-04 20:00:00', '2018-06-04 20:30:00',
                       '2018-06-04 21:00:00', '2018-06-04 21:30:00',
                       '2018-06-04 22:00:00', '2018-06-04 22:30:00',
                       '2018-06-04 23:00:00', '2018-06-04 23:30:00'],
                      dtype='datetime64[ns]', name='date', length=95088, freq=None)
    


Datetime-indexed :obj:`pandas.DataFrame` objects have useful features for
time series analysis like grouping and calculating statistics by time
aggregates. The example below shows how to calculate the day of year
mean for energy balance components, it also demonstrates how to use the
``add_lines`` plotting method available to :obj:`.Data`, :obj:`.QaQc`, and
:obj:`.Plot` objects.

    >>> # convert to internal names, copy dataframe
    >>> df = d.df.rename(columns=d.inv_map)
    >>> # day of year mean of input energy balance components
    >>> vars_we_want = ['H', 'LE', 'Rn', 'G']
    >>> doy_means = df[vars_we_want].groupby(d.df.index.dayofyear).mean()
    >>> # create a Bokeh figure
    >>> fig = figure(x_axis_label='day of year', y_axis_label='day of year mean (w/m2)')
    >>> # arguements needed for creating interactive plots
    >>> plt_vars = vars_we_want
    >>> colors = ['red', 'blue', 'black', 'green']
    >>> x_name = 'date'
    >>> source = ColumnDataSource(doy_means)
    >>> Plot.add_lines(fig, doy_means, plt_vars, colors, x_name, source, labels=vars_we_want,
    >>>     x_axis_type=None) 
    >>> show(fig)

.. raw:: html
    :file: _static/tutorial/doy_mean_example.html
    

.. Note::
   The ``x_axis_type=None`` is a unique argument to :meth:`.Plot.add_lines` and
   :meth:`.Plot.line_plot` that in this case means to not try to force the
   x-axis format to a datetime representation, default is
   ``x_axis_type='date'``.

.. seealso::
   Some routines occur automatically when creating a :obj:`.Data` object,
   including calcuation of weighted and non-weighted averages of soil heat flux
   and soil moisture which is described in :ref:`Averaging data from multiple
   sensors`.

Modifying input data
^^^^^^^^^^^^^^^^^^^^

A last note on the :obj:`.Data` object (same goes for the :obj:`.QaQc` object)
is that :attr:`Data.df` is a class property, in this case that means that it
can be reassigned with a different :obj:`pandas.DataFrame`. This is
critical for manual pre-filtering and validation of data before
proceeding with energy balance closure routines. A simple example is
shown here:

    >>> # add 5 to air temperature, this would effect ET calculations later
    >>> x = d.df
    >>> x['T_SONIC'] += 5
    >>> d.df = x

A realistic use of the reassignability of the :attr:`.Data.df` and
:attr:`.QaQc.df` properties is shown in `manual cleaning of poor quality
data <https://flux-data-qaqc.readthedocs.io/en/latest/closure_explanation.html#step-0-manual-cleaning-of-poor-quality-data>`__.

.. seealso::
   The :meth:`.Data.apply_qc_flags` method allows for reading in quality
   control flags with the input data and filtering specific data out based on
   user-defined numeric or character flags. This routine is specific to
   :obj:`.Data` and includes several attributes that are added to a
   :obj:`.Data` instance, for full explanation and examples see
   :ref:`Quality-based data filtering`.

Visualize input data
--------------------

The :meth:`.Data.plot` method create a series of interactive time series
plots of input data, potential plots inlcude:

-  energy balance components
-  radiation components
-  multiple soil heat flux measurements
-  air temperature
-  vapor pressure and vapor pressure deficit
-  wind speed
-  precipitation
-  latent energy
-  multiple soil moisture measurements

If any of these variables are not found the plot(s) will not be added.

The most useful interactive features of plots created by
``flux-data-qaqc`` are:

-  pan/zoom
-  hover tolltips on var names, values, date
-  linked x-axes on time series plots
-  save plot option (can save specific subplot zoomed in)

Here is an example,

    >>> d.plot(output_type='notebook', plot_width=700)

The output plot is not shown in the online documentation due to 
memory constraints. 


.. hint:: 
   The plot methods of :obj:`.Data` and :obj:`.QaQc` objects have the keyword
   argument ``output_type`` which by default is set to “save”, the other two
   options are “notebook” for showing within a Jupyter Notebook and “show”
   which opens a temporary file in the default web browser.

If you rather save the plot, and maybe you want 2 columns of plots,

    >>> d.plot(ncols=2, plot_width=500) 

After saving a plot without specifying the output file path (keyword
argument ``out_file``), it will be saved to an “output” directory where
the config file is with the file name based on :attr:`.Data.site_id` with the
suffix “\_input_plots”:

    >>> # where the plot file was saved by default
    >>> d.plot_file
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/output/US-Tw3_input_plots.html')

The following plot is not shown due to excessive memory usage needed to build 
online documentation.

    >>> # view outplot plots within Jupyter notebook
    >>> from IPython.display import HTML
    >>> HTML(filename=d.plot_file)


.. hint:: 
   The :meth:`.QaQc.plot` method shown below is similar however it may include
   added plots with calculated and corrected variables (if they exist) and will
   always plot data in daily and monthly temporal frequency because daily
   frequency is required before applying ``flux-data-qaqc`` energy balance
   closure corrections.

Temporal resampling
-------------------

The :obj:`.QaQc` object holds several tools for managing data and eddy
covariance data analysis, but one of it’s primary features is temporal
resampling of input data to daily and monthly frequencies. The
resampling of time series data to daily frequency occurs upon the
creation of a :obj:`.QaQc` instance if the frequency within the preceeding
:obj:`.Data` object is not already daily:

    >>> # the frequency of the input data is 30 minute
    >>> d.df.index[0:5]
        DatetimeIndex(['2013-01-01 00:00:00', '2013-01-01 00:30:00',
                       '2013-01-01 01:00:00', '2013-01-01 01:30:00',
                       '2013-01-01 02:00:00'],
                      dtype='datetime64[ns]', name='date', freq=None)

    >>> # creating a QaQc instance will automatically convert to daily
    >>> q = QaQc(d)
        The input data temporal frequency appears to be less than daily.
        Data is being resampled to daily temporal frequency.
        Filtering days with less then 100.0% or 48/48 sub-daily measurements
        Converting vpd from hpa to kpa


    >>> # first 5 datetime indices are dates now
    >>> q.df.index[0:5]
        DatetimeIndex(['2013-01-01', '2013-01-02', '2013-01-03', '2013-01-04',
                       '2013-01-05'],
                      dtype='datetime64[ns]', name='date', freq=None)



The method used for aggregating different variables, e.g. mean or sum,
when resampling to daily or monthly frequency is defined in the
:attr:`QaQc.agg_dict` class attribute:

    >>> # these are the internal names as keys and temporal aggregation method as values
    >>> QaQc.agg_dict
        {'energy': 'mean',
         'flux': 'mean',
         'flux_corr': 'mean',
         'br': 'mean',
         'ET': 'sum',
         'ET_corr': 'sum',
         'ET_gap': 'sum',
         'ET_fill': 'sum',
         'ET_fill_val': 'sum',
         'ET_user_corr': 'sum',
         'ebr': 'mean',
         'ebr_corr': 'mean',
         'ebr_user_corr': 'mean',
         'ebr_5day_clim': 'mean',
         'gridMET_ETr': 'sum',
         'gridMET_prcp': 'sum',
         'lw_in': 'mean',
         't_avg': 'mean',
         'rso': 'mean',
         'sw_pot': 'mean',
         'sw_in': 'mean',
         'vp': 'mean',
         'vpd': 'mean',
         'ppt': 'sum',
         'ws': 'mean',
         'Rn': 'mean',
         'sw_out': 'mean',
         'lw_out': 'mean',
         'G': 'mean',
         'LE': 'mean',
         'LE_corr': 'mean',
         'LE_user_corr': 'mean',
         'H': 'mean',
         'H_corr': 'mean',
         'H_user_corr': 'mean'}



.. note:: 
   There are several calculated variables above that may not look familiar,
   many are calculated by the energy balance closure correction routines and
   described in :ref:`Closure Methodologies`.  Also, any other variables (not
   found in :attr:`.QaQc.agg_dict` that exist in a :attr:`.QaQc.df` before
   accessing :attr:`.QaQc.monthly_df` the first time will be averaged in the
   monthly time series dataframe (:attr:`.QaQc.monthly_df`).

The :obj:`.QaQc` constructor tries to infer the temporal frequency of the
input time series data, however the method is not always accurate, to
access the inferred initial temporal frequency of the data view the
:attr:`.QaQc.temporal_freq` attribute:

    >>> q.temporal_freq
        '30T'



If the inferred input frequency was accurate you will see a `Pandas datetime
alias
<https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__,
in this case ‘30T’ is thirty minutes. If the temporal frequency is not automatically detected you should be able to rely on the ``n_samples_per_day`` instance attribute that is manually estimated by the ``QaQc`` constructor:

    >>> q.n_samples_per_day
        48

Filter days with sub-daily gaps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``drop_gaps`` and ``daily_frac`` keyword arguments used when creating a :obj:`.QaQc` instance allow you to control how days with sub-daily measurement gaps will or will not be filtered out when resampling to daily frequency.

Sub-daily gaps in energy balance variables :math:`LE`, :math:`H`, :math:`Rn`, and :math:`G` , and daily ASCE standaridized reference ET inputs, e.g. hourly :math:`ea` ("vp"), :math:`rs` ("sw_in"), :math:`t_min`, :math:`t_max`, and :math:`ws`, can be linearly interpolated automatically before daily aggregations. Interpolation is performed over gap lengths measured in hours, with options to control the longest length of gap to interpolate when :math:`Rn \ge 0` controlled by the :obj:`.QaQc` keyword argument ``max_interp_hours`` (default 2 hours) and the longest gap to interpolate when :math:`Rn < 0` set by the ``max_interp_hours_night`` (default 4 hours).:math:`

.. Important:: 

   By default the :obj:`.QaQc` constructor will first linearly interpolate
   energy balance and ASCE ref. ET variables (:math:`LE`, :math:`H`,
   :math:`Rn`, :math:`G`, :math:`ea` ("vp"), :math:`rs` ("sw_in"),
   :math:`t_min`, :math:`t_max`, and :math:`ws`) according to the maximum gap
   lengths (``max_interp_hours`` and ``max_interp_hours_night``) and then count
   sub-daily gaps and drop days (set values to null) for all climate data
   columns (not QC flag or sub-daily gap count columns) where any of the
   sub-daily data are missing because by default ``drop_gaps=True`` and
   ``daily_frac=1.0``. In other words, if you have hourly input data
   for(:math:`LE` and one hour was missing on a given day, by default that hour
   will be linearly interpolated before calculating the daily time series and
   the daily mean will be calculated after. On the other hand, if other climate
   variables had a single hour missing on a given day, e.g. wind direction,
   this day would be filtered out by the :obj:`.QaQc` constructor.  This is
   important because the daily time series is what is used in all energy
   balance closure correction and daily ASCE standardized reference ET
   algorithms.

The percentage of sub-daily samples to require set by the ``daily_frac`` argument and the maximum length of gaps to linearly interpolate set by ``max_interp_hours`` and ``max_interp_hours_night`` complement each other and are used in tandem. For example, if the input data is half-hourly and you only want a maximum of 4 hours to be interpolated on any given day and gap lengths to interpolate should be no more than 2 hours each then you would pass the following parameters to the :obj:`.QaQc` constructor:


    >>> q = QaQc(d, daily_frac=20/24, max_interp_hours=2, max_interp_hours_night=2)
        The input data temporal frequency appears to be less than daily.
        Data is being resampled to daily temporal frequency.
        Linearly interpolating gaps in energy balance components up to 2 hours when Rn < 0 and up to 2 hours when Rn >= 0.
        Filtering days with less then 83.33333333333334% or 40/48 sub-daily measurements    

In this case we set ``daily_frac=20/24`` because we are only allowing a maximum of 4 hours of total gaps in the day in other words we are requiring 40 of the 48 half hourly samples to exist before we filter out a day. Remember, because linear interpolation of gaps is done before counting sub-daily gaps, this could result in retaining days with more than 4 hours of gaps in the original time series of energy balance components. You may also pass the ``daily_frac`` arugment as a decimal fraction, e.g. :math:`0.8333 \approx 20/24`.

To not drop any days and take daily means/sums based on whatever data exists in a given day *without* any interpolation of energy balance variables,

    >>> q = QaQc(d, drop_gaps=False, max_interp_hours=None)
        The input data temporal frequency appears to be less than daily.
        Data is being resampled to daily temporal frequency.


Let’s view a comparison of :math:`Rn` using different options of
filtering days with sub-daily gaps in the working dataset, because it
has several periods of systematic gaps which cause upwards skewing of
daily mean :math:`Rn` if not filtered carefully:

    >>> # make an empty pandas dataframe for Rn series
    >>> Rn_df = pd.DataFrame()
    >>> # recreate multiplt QaQc instances using different sub-day gap filters
    >>> q = QaQc(d, drop_gaps=False, max_interp_hours=None)
    >>> Rn_df['sub_day_gaps'] = q.df.Rn_subday_gaps
    >>> Rn_df['no_filter_no_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, drop_gaps=False)
    >>> Rn_df['no_filter_with_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=0.5) # filter days with less than 50% data
    >>> Rn_df['require_50'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=0.75)
    >>> Rn_df['require_75'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=1, max_interp_hours=24, max_interp_hours_night=24) 
    >>> Rn_df['require_100_with_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> q = QaQc(d, daily_frac=1, max_interp_hours=None) 
    >>> Rn_df['require_100_no_interp'] = q.df.rename(columns=q.inv_map).Rn
    >>> # plot to compare results of day-gap filter
    >>> fig = figure(x_axis_label='date', y_axis_label='mean daily net radiation (w/m2), filtered based on sub-daily gaps')
    >>> # arguments needed for creating interactive line plots
    >>> colors = ['red', 'darkred','orange', 'blue', 'black', 'tan']
    >>> plt_vars = ['no_filter_no_interp', 'no_filter_with_interp', 'require_50', 'require_75', 'require_100_with_interp', 'require_100_no_interp']
    >>> labels = ['no filter wout/interp.', 'no filter w/interp.', 'require > 50% w/interp.', 'require > 75% w/interp.', 'require 100% w/interp.', 'require 100% wout/interp.']
    >>> x_name = 'date'
    >>> source = ColumnDataSource(Rn_df)
    >>> Plot.add_lines(fig, Rn_df, plt_vars, colors, x_name, source, labels=labels) 
    >>> # add daily gap counts to secondary y
    >>> fig.extra_y_ranges['gap_counts'] = Range1d(start=0, end=48)
    >>> fig.add_layout(LinearAxis(y_range_name='gap_counts', axis_label='number of sub-daily gaps'), 'right')
    >>> fig.circle('date', 'sub_day_gaps', legend='n sub-day gaps', y_range_name='gap_counts',
    >>>     color='silver', source=source
    >>> )
    >>> fig.hover[0].tooltips.append(('sub_day_gaps','@{}'.format('sub_day_gaps')))
    >>> fig.legend.location = 'top_right'
    >>> show(fig)


.. raw:: html 
    :file: _static/tutorial/filter_subday.html 

Try zooming in on the gaps filled by the "no filter wout/interp." line to compare which days are retained/filtered by different options, also remove lines by clicking on them in the legend to compare subsets of options.

.. Tip:: 
   For a more fine-grained approach to filtering out days where perhaps
   multiple 2 hour gaps were filled use the newly created daily gap count
   columns: "LE_subday_gaps", "H_subday_gaps", "Rn_subday_gaps", and
   "G_subday_gaps":

      >>> q = QaQc(d)
      >>> df = q.df.rename(columns=q.inv_map)

   For example, you could post-filter out days in any given energy balance
   variable, in this case :math:`Rn` where sub-daily gaps exceed a threshold:

      >>> df.loc[(df.Rn_subday_gaps > 4) & (df.Rn.notna()), ['Rn','Rn_subday_gaps']]

        .. raw:: html

            <div>
            <style scoped>
                .dataframe tbody tr th:only-of-type {
                    vertical-align: middle;
                }
            
                .dataframe tbody tr th {
                    vertical-align: top;
                }
            
                .dataframe thead th {
                    text-align: right;
                }
            </style>
            <table border="1" class="dataframe">
              <thead>
                <tr style="text-align: right;">
                  <th></th>
                  <th>Rn</th>
                  <th>Rn_subday_gaps</th>
                </tr>
                <tr>
                  <th>date</th>
                  <th></th>
                  <th></th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>2015-06-09</td>
                  <td>101.710194</td>
                  <td>5.0</td>
                </tr>
                <tr>
                  <td>2015-11-20</td>
                  <td>47.990988</td>
                  <td>5.0</td>
                </tr>
                <tr>
                  <td>2016-01-15</td>
                  <td>72.495973</td>
                  <td>8.0</td>
                </tr>
                <tr>
                  <td>2018-01-06</td>
                  <td>79.507008</td>
                  <td>7.0</td>
                </tr>
                <tr>
                  <td>2018-05-10</td>
                  <td>160.997332</td>
                  <td>6.0</td>
                </tr>
              </tbody>
            </table>
            </div>
          </br>

Monthly time series
^^^^^^^^^^^^^^^^^^^

The :attr:`.QaQc.monthly_df` property allows for creating the monthly time
series of input anc calculated variables provided by
:meth:`.QaQc.correct_data`. It uses the same temporal aggregation methods as
the daily time series i.e. from :attr:`.QaQc.agg_dict`. Although there are
many similarities there are important differences between :attr:`.QaQc.df`
and :attr:`.QaQc.monthly_df` other than the obvious: when accessing the
:attr:`.QaQc.monthly_df` it will automatically run the default energy balance
closure correction routine provided by :meth:`.QaQc.correct_data` *if* it has
not yet been run. You can check if it has been run at anytime by:

    >>> q.corrected
        False

To show how this works let’s access the monthly data and show the
monthly statistics of the “corrected” evapotranspiration (ET_corr):


    >>> # first note, ET_corr is not in the dataset yet
    >>> 'ET_corr' in q.df.columns
        False

Now access the monthly time series,

    >>> q.monthly_df;
    >>> 'ET_corr' in q.df.columns
        True

By calling the monthly dataframe, the energy balance closure was applied
automatically

    >>> q.monthly_df.ET_corr.describe()
        count     61.000000
        mean      87.858135
        std       49.938287
        min       11.370062
        25%       41.418994
        50%       84.383190
        75%      127.500125
        max      192.033481
        Name: ET_corr, dtype: float64


    >>> q.corrected
        True

.. note:: 
   The :attr:`.QaQc.monthly_df` also filters out months with less than 30% of
   days of the month missing by default. To calculate monthly time series with
   other threshold fractions of days required use the
   :func:`.util.monthly_resample` function and adjust the keyword argument
   ``thresh`` which is the fraction (0-1) of days of the month required to not
   be gaps otherwise the month’s value will be forced to null, e.g. if you
   wanted to caclulate the monthly mean air temperature requiring 30 and 90
   percent of the days in the month to not be gaps:

    >>> from fluxdataqaqc.util import monthly_resample
    >>> # select just t_avg for example
    >>> cols = ['t_avg'] 
    >>> df = q.df.rename(columns=q.inv_map)
    >>> # create temporary df with different monthly resample results
    >>> tmp_df = monthly_resample(df, cols, 'mean', thresh=0.9).rename(
    >>>     columns={'t_avg': 'thresh_90'}
    >>> )
    >>> # join temp dataframe with monthly resample results using different thresh
    >>> monthly_gap_comp = tmp_df.join(monthly_resample(df, cols, 'mean', thresh=0.3).rename(
    >>>     columns={'t_avg': 'thresh_30'})
    >>> )
    >>> # plot to compare results of day-gap filter
    >>> fig = figure(x_axis_label='date', y_axis_label='monthy mean air temperature (C), filtered based on daily gaps')
    >>> # arguments needed for creating interactive line plots
    >>> x = 'date'
    >>> source = ColumnDataSource(monthly_gap_comp)
    >>> # this example also shows how to use other Bokeh plot arguments
    >>> Plot.line_plot(fig,'date','thresh_30',source,'red',label='require > 30%', line_alpha=0.5) 
    >>> Plot.line_plot(fig,'date','thresh_90',source,'black',label='require > 90%',line_dash='dotted', line_width=2) 
    >>> fig.legend.location = 'top_right'
    >>> show(fig)

.. raw:: html 
    :file: _static/tutorial/filter_monthly.html

.. raw:: html 

   <br>

Energy balance corrections
--------------------------

``flux-data-qaqc`` provides routines that adjust surface energy balance fluxes
to improve energy balance closure of eddy covariance flux station data.
These routines ultimately result in a corrected daily and monthly time series
of latent energy, sensible heat, and evapotranspiration with the option to
gap-fill days in corrected ET with ET calculated from gridMET reference ET and
fraction of reference ET.

There are two methods currently implemented: 

* Energy Balance Ratio method (default), modified from the `FLUXNET method <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
* Bowen Ratio approach (forces closure) 
* Multiple least squares regression
  - user defines LHS and RHS from :math:`LE`, :math:`H`, :math:`Rn`, and :math:`G`,

Detailed descriptions of methods including daily ET gap-filling methods
can be found in the online documentation :ref:`Closure Methodologies`
page. A few important notes on the API of these methods and other
hydro-climatic statistical variables that are calculated are shown
below.

The :meth:`.QaQc.correct_data` method is used to run energy balance closure
corrections. Here are a few tips on using them,

    >>> # note above the monthly_df automatically applied the 'ebr' Energy Balance Ratio correction
    >>> q.corr_meth
        'ebr'

    >>> # potential correction options
    >>> q.corr_methods
        ('ebr', 'br', 'lin_regress')


    >>> # to specify the Bowen Raito method:
    >>> q.correct_data(meth='br')

    >>> # the most recently used correction method is now shown
    >>> q.corr_meth
        'br'

.. Tip:: 
   After applying any energy balance closure correction routine all previous
   corrected variables will be overwritten or dropped in :attr:`.QaQc.df`,
   :attr:`.QaQc.monthly_df`, and :attr:`.QaQc.variables`, therefore to make a comparison
   of different methods on the same data make a copy of the ``df`` or
   ``monthly_df`` properties before running the next correction, e.g.

    >>> # make copies of daily results of different correction options
    >>> q.correct_data(meth='ebr')
    >>> ebr_gapfilled = q.df
    >>> q.correct_data(meth='ebr', etr_gap_fill=False)
    >>> ebr_notgapfilled = q.df
    >>> q.correct_data(meth='br')
    >>> br_gapfilled = q.df
    >>> q.correct_data(meth='br', etr_gap_fill=False)
    >>> br_notgapfilled = q.df


ET gap-filling
^^^^^^^^^^^^^^

A few notes on the option that uses reference ET and fraction of daily
reference ET to fill in large gaps in corrected ET, i.e. the keyword
argument ``QaQc.correct_data(etr_gap_fill = True)``.

-  The nearest `gridMET <http://www.climatologylab.org/gridmet.html>`__
   cell’s time series data for precipitation and alfalfa reference ET is
   attempted to be downloaded if it is not found in the
   ``gridmet_file_path`` entry of the config.ini file.

-  If the path to a gridMET file is not found it is re-downloaded, the
   config file will be updated with the new path and resaved.

-  Only the overlapping time period that matches the eddy covariance
   time series data is attempted to be downloaded, i.e. the period in
   ``QaQc.df.index``.

-  When a gridMET file is downloaded it will always be saved in a
   subdirectory where the config file is located called “gridMET_data”
   and named using the :attr:`QaQc.site_id` and gridMET cell centroid
   latitude and longitude.

-  Corrected latent energy (:math:`LE_{corr}`) gaps are also backwards
   filled from gap-filled ET.

.. caution:: 
   `gridMET <http://www.climatologylab.org/gridmet.html>`__ only exists within
   the contiguous United States and from 1979 to present, therefore if your
   station lies outside of this region or you are analyzing eddy flux data
   recorded before 1979 this option will not be ususable and you should always
   run corrections with ``etr_gap_fill=False`` to avoid potential errors.


The Bowen Ratio correction method will produce the ‘br’ variable which
is the Bowen Ratio.

Other calculations
------------------

By default, :meth:`.QaQc.correct_data` also calculates ET from input latent
energy (LE) and air temperature, corrected ET from corrected LE and air
temperature, potential clear sky radiation (ASCE formulation), and the
:obj:`.Data` object attempts to calculate vapor pressure deficit from vapor
pressure and air temperature or vapor pressure from vapor pressure
deficit and air temperature if they exist at hourly or shorter temporal
frequency.

Evapotranspiration
^^^^^^^^^^^^^^^^^^

The evapotranspiration (ET) calculations are described in :ref:`Steps 7 and 8 correct turbulent fluxes, EBR, and ET` of the Energy Balance Ratio correction explanation.

ASCE clear sky radiation
^^^^^^^^^^^^^^^^^^^^^^^^

Daily ASCE potential clear sky radiation (:math:`R_{so}`) is calculated
using equation 19 in the “ASCE Standardized Reference Evapotranspiration
Equation” final report by the Task Committee on Standardization of
Reference Evapotranspiration Environmental and Water Resources Institute
of the American Society of Civil Engineers January, 2005
`here <https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf>`__.
This calculation is a simple method based primarily on elevation and
latitude which results in a theoretical envelope of :math:`R_{so}` as a
function of the day of year,

.. math::  R_{so} = \left(5 + 2 \times 10^{-5} z \right) R_a 

where :math:`z` is elevation in meters and :math:`R_a` is daily
extraterrestrial radiation (radiation with in the absence of an
atmosphere), which itself is a well-behaved function of solar
declination, the day of the year and the solar constant (see equations
21-29 in the `ASCE
report <https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf>`__).

Vapor pressure/deficit
^^^^^^^^^^^^^^^^^^^^^^

The :obj:`.Data` object will attempt to calculate vapor pressure or vapor
pressure deficit if one exists but not the other and average air
temperature time series also exists with the input data at hourly or
shorter temporal frequency. The Magnus equation (eqn. 37 in the ASCE report) 
states that the saturation vapor pressure (:math:`es`) in kPa relates to air temperature,

.. math::  es = 0.6108  e^{\left(\frac{17.27 \cdot T}{(T + 237.3)}\right)} 

where :math:`T` is average hourly air temperature in degrees celcius.
Vapor pressure deficit (:math:`vpd`) is,

.. math::  vpd = es - ea,

where :math:`ea` is actual vapor pressure in kPa. **Note,** The
equations above are defined for hourly measurements however they are used for
hourly or shorter mean variables (:math:`T`, :math:`ea`, or :math:`vpd`)
within ``flux-data-qaqc`` and then converted to daily means, if they are
not present in the input data at hourly or shorter frequencies then they
are not calculated.

These equations can be rearanged to solve for either :math:`es` or
:math:`vpd` given the other variable and air temperature. For example,
if given :math:`T` and :math:`vpd`, then to get actual vapor pressure

.. math::  es = 0.6108  e^{\left(\frac{17.27 \cdot T}{(T + 237.3)}\right)}  

.. math::  ea = es - vpd. 

In ``flux-data-qaqc`` actual vapor pressure is named "vp" not "ea". Also, during these calculations, if relative humidity is not found in the input dataset then it will subsequently be estimated as 

.. math:: rh = 100 \times \frac{ea}{es}.

.. hint:: 
   The same calculations are available at the daily timestep but are not
   automatically applied as the hourly or higher temporal frequency calculation
   is preffered. To apply the estimates of vapor pressure or vapor pressure
   deficit, and saturation vapor pressure, and relative humidity with daily data
   one must call the :meth:`.QaQc._calc_vpd_from_vp` method from a :obj:`.QaQc`
   instance. 
   
Calculated variable reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although variables created by energy balance closure corrections are described
in :ref:`Closure Methodologies` and others are below, here is a reference list 
of all possibly calculated variables created by ``flux-data-qaqc``:
    

   +--------------+-----------------------------------------------------------------+
   |Variable      | Description                                                     |
   +==============+=================================================================+
   |rso           | potential clear sky radiation (ASCE formulation)                |
   +--------------+-----------------------------------------------------------------+
   |flux          | input LE + H                                                    |
   +--------------+-----------------------------------------------------------------+
   |energy        | input Rn - G                                                    |
   +--------------+-----------------------------------------------------------------+
   |ebr_5day_clim | 5 day climatology of the filtered Energy Balance Ratio          |
   +--------------+-----------------------------------------------------------------+
   |LE_corr       | corrected latent energy                                         |
   +--------------+-----------------------------------------------------------------+
   |ebc_cf        | energy balance closure correction factor (inverse of ebr_corr)  |
   +--------------+-----------------------------------------------------------------+
   |ebr_corr      | corrected energy balance ratio                                  |
   +--------------+-----------------------------------------------------------------+
   |flux_corr     | LE_corr + H_corr                                                |
   +--------------+-----------------------------------------------------------------+
   |ebr           | input energy balance ratio                                      |
   +--------------+-----------------------------------------------------------------+
   |br            | bowen ratio                                                     |
   +--------------+-----------------------------------------------------------------+
   |H_corr        | corrected sensible heat                                         |
   +--------------+-----------------------------------------------------------------+
   |ET            | ET calculated from input LE and average air temperature         |
   +--------------+-----------------------------------------------------------------+
   |ET_corr       | ET calculated from LE_corr and avg. air temp.                   |
   +--------------+-----------------------------------------------------------------+
   |gridMET_ETr   | gridMET alfalfa reference ET (nearest cell)                     |
   +--------------+-----------------------------------------------------------------+
   |gridMET_prcp  | gridMET precipitation                                           |
   +--------------+-----------------------------------------------------------------+
   |ETrF          | fraction of reference ET for ET_corr, i.e. ET_corr / gridMET_ETr|
   +--------------+-----------------------------------------------------------------+
   |ETrF_filtered | filtered ETrF                                                   |
   +--------------+-----------------------------------------------------------------+
   |ET_fill       | gridMET_ETr * ETrF_filtered (fills gaps in ET_corr)             |
   +--------------+-----------------------------------------------------------------+
   |ET_gap        | True on gap days in ET_corr, False otherwise                    |  
   +--------------+-----------------------------------------------------------------+
   |ET_fill_val   | value of ET_fill on gap days                                    |  
   +--------------+-----------------------------------------------------------------+


A note on units
---------------

Upon creation of a :obj:`.QaQc` object, variables are checked for valid input
units and converted to required units needed for internal calculations when
running :meth:`.QaQc.correct_data` and for certain default plots (see below).
For a list of valid input units for different variables refer to the
:attr:`.QaQc.allowable_units` attribute:

    >>> q.allowable_units
        {'LE': ['w/m2'],
         'H': ['w/m2'],
         'Rn': ['w/m2'],
         'G': ['w/m2'],
         'lw_in': ['w/m2'],
         'lw_out': ['w/m2'],
         'sw_in': ['w/m2'],
         'sw_out': ['w/m2'],
         'ppt': ['mm', 'in'],
         'vp': ['kpa', 'hpa'],
         'vpd': ['kpa', 'hpa'],
         't_avg': ['c', 'f']}


For each variable above, if given one of the units allowable the units
will automatically be converted to the required units.

To know which variables are required to be in particular units view
:attr:`.Qc.required_units`:

    >>> q.required_units
        {'LE': 'w/m2',
         'H': 'w/m2',
         'Rn': 'w/m2',
         'G': 'w/m2',
         'lw_in': 'w/m2',
         'lw_out': 'w/m2',
         'sw_in': 'w/m2',
         'sw_out': 'w/m2',
         'ppt': 'mm',
         'vp': 'kpa',
         'vpd': 'kpa',
         't_avg': 'c'}



.. note:: 
   The list of allowable units is a work in progress, if your input units are
   not available consider raising an issue on `GitHub
   <https://github.com/Open-ET/flux-data-qaqc/issues>`__ or providing the
   conversion directly with a pull request. Automatic unit conversions are
   handled within the :obj:`.util.Convert` class using the
   :meth:`.util.Convert.convert` class method.

Save resampled and corrected data
---------------------------------

The :meth:`.QaQc.write` method conveniently writes daily and monthly time
series of input and calculated variables to comma separated value (CSV)
files. If the :meth:`.QaQc.correct_data` method has not yet been run it will
be and the monthly time series will also be created using the default
parameters for the correction routine (Energy Balance Ratio method with
ETr-based gap filling).

The default output directory for time series files can be
accessed/changed by the ``out_dir`` attribute, if not changed it will be
located in the same directory of the config.ini file. The daily and
monthly time series file names will begin with the :attr:`.QaQc.site_id`
followed by “daily_data” or “monthly_data” resepctively. For example,

    >>> # new QaQc instance
    >>> q = QaQc(d)
    >>> # a platform dependent pathlib.Path object
    >>> q.out_dir
        PosixPath('/home/john/flux-data-qaqc/examples/Basic_usage/output')

The line below shows that no output files have been written to
:attr:`.QaQc.out_dir` yet,

    >>> # print files in output directory that begin with the site_id
    >>> [f.name for f in q.out_dir.glob('{}*'.format(q.site_id))]
        ['US-Tw3_input_plots.html']

    >>> q.corrected
        False

    >>> # writing files also ran corrections since they were not yet run
    >>> q.write()
    >>> q.corrected
        True


Now the respective daily and monthly time series have been written to
:attr:`.QaQc.out_dir`,

    >>> [f.name for f in q.out_dir.glob('{}*'.format(q.site_id))]
        ['US-Tw3_daily_data.csv', 'US-Tw3_monthly_data.csv', 'US-Tw3_input_plots.html']


.. hint:: 
   You can overwrite the default name of the output directory to save the daily
   and monthly time series using the ``out_dir`` keyword argument to
   :meth:`.QaQc.write`, this option keeps the location within the directory of
   the config file but just changes the name, whereas to change the entire
   output directory path adjust the :attr:`.QaQc.out_dir` attribute directly.
   Also, the naming scheme of output files created will use user-defined names
   for all input variables.

Visualize resampled and corrected data
--------------------------------------

Similar to the :meth:`.Data.plot`, the :meth:`.QaQc.plot` method creates a
series of default time series and scatter plots of input and in this
case calculated variables. The temporal frequency of plots from
:meth:`.QaQc.plot` will always be daily and monthly and additional plots are
created for validation of energy balance closure corrections, otherwise the
same options such as number of subplot columns, super title, subplot
dimensions, output type, output file path, etc. are available. Similar to
:meth:`.QaQc.write` and :attr:`.QaQc.monthly_df`, if the data has not yet been
corrected the ``plot`` method will correct it using the default parameters
before creating the plots. 

Here is an example of the default daily and monthly time series plots produced after running the Energy Balance Ratio closure correction:

    >>> q = QaQc(d)
    >>> q.plot(output_type='notebook', plot_width=700)

.. raw:: html
    :file: _static/tutorial/US-Tw3_plots.html


==========================================================
flux-data-qaqc - Tools for Energy Balance Closure Analysis
==========================================================



|Docs| |Downloads per month| |PyPI version| 

-----------

`View on GitHub <https://github.com/Open-ET/flux-data-qaqc>`_. 

``flux-data-qaqc`` provides a framework to create reproducible workflows for validation and analysis of eddy covariance data. The package is intended for those who need to post-process flux data, particularly for generating daily and monthly evapotransipration (ET) timeseries estimates with energy balance closure corrections applied. Applications where this software may be useful include analysis involving eddy covariance flux tower data, hydrologic or atmospehric model validation, and irrigation and water consumption studies. 

Key functionalities and tools include:

* data validation with methods for quality-based filtering
* time series data tools, e.g. temporal aggregation and resampling
* management of site metadata, data provenance, and file structure
* energy balance closure algorithms and other meterological calculations
* downloading and management of `gridMET <http://www.climatologylab.org/gridmet.html>`__ meterological data
* customizable and interactive data visualizations
* batch processing and unit conversions

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   Configuration Options <advanced_config_options>
   Tutorial <tutorial>
   Closure Algorithms <closure_explanation>
   api
   Tests <software_tests>
   contributors

.. toctree::
   :maxdepth: 1
   
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |Docs| image:: https://readthedocs.org/projects/flux-data-qaqc/badge/?version=latest
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |Downloads per month| image:: https://img.shields.io/pypi/dm/fluxdataqaqc.svg
   :target: https://pypi.python.org/pypi/fluxdataqaqc/

.. |PyPI version| image:: https://img.shields.io/pypi/v/fluxdataqaqc.svg
   :target: https://pypi.python.org/pypi/fluxdataqaqc/
Installation
------------

Using PIP:

.. code-block:: bash

   pip install fluxdataqaqc

PIP should install the necessary dependencies however it is recommended to use
conda and first install the provided virtual environment. This is useful to
avoid changing your local Python environment. Note, ``flux-data-qaqc`` has been
tested for Python 3.7+, although it may work with versions greater than or
equal to 3.4.

First make sure you have the ``fluxdataqaqc`` environment file, you can download it `here <https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/environment.yml?token=AB3BJKUKL2ELEM7WPLYLXFC45WQOG>`_. Next to install run,

.. code-block:: bash

   conda env create -f environment.yml

To activate the environment before using the ``flux-data-qaqc`` package run,

.. code-block:: bash

   conda activate fluxdataqaqc

Now install using PIP:

.. code-block:: bash

   pip install fluxdataqaqc

Now all package modules and tools should be available in your Python environment PATH and able to be imported. Note if you did not install the Conda virtual environment above, PIP should install dependencies automatically but be sure to be using a version of Python above or equal to 3.4. To test that everything has installed correctly by opening a Python interpretor or IDE and run the following:

.. code-block:: python

   import fluxdataqaqc

and 

.. code-block:: python

   from fluxdataqaqc import Data, QaQc, Plot

If everything has been installed correctly you should get no errors. 



Developer mode
^^^^^^^^^^^^^^

If you plan on contributing to ``flux-data-qaqc`` you can install it in
*developer mode* which allows you to easily modify the source code and
test changes within a Python environment.

Clone or download from `GitHub <https://github.com/Open-ET/flux-data-qaqc>`_.  

.. code-block:: bash

   git clone https://github.com/Open-ET/flux-data-qaqc.git

Next move to the root directory and install and activate the conda 
environment:

.. code-block:: bash

   conda env create -f environment.yml

Run the following to install ``flux-data-qaqc`` in developer mode into
your environment,

.. code-block:: bash

   pip install -e .


Closure Methodologies
=====================

``flux-data-qaqc`` currently provides two routines which ultimately adjust
turbulent fluxes in order to improve energy balance closure of eddy covariance
tower data, the Energy Balance Ratio and the Bowen Ratio method. 

Closure methods are assigned as keyword arguments to the 
:meth:`.QaQc.correct_data` method, and for a list of provided 
closure options see :attr:`.QaQc.corr_methods`.
For example, if you would like to run the Bowen Ratio correction routine
assuming you have succesfully created a :obj:`.QaQc` object,

.. code-block:: python

    # q is a QaQc instance
    q.correct_data(meth='br')

The other keyword argument for :meth:`.QaQc.correct_data` allows for gap filling corrected evapotranspiration (:math:`ET`) which is calculated from corrected latent energy (:math:`LE`). By default the gap filling option is set to True, more details on this below in :ref:`Step 9, optionally gap fill corrected ET using gridMET reference ET and reference ET fraction`.

.. Tip:: 
   All interactive visualizations in this page were created using
   :meth:`.Plot.line_plot`, :meth:`.Plot.add_lines`, and
   :meth:`.Plot.scatter_plot` which automatically handle issues with utilizing
   the mouse hover tooltips and other :obj:`bokeh.plotting.figure.Figure`
   features.
    
Data description
----------------

The data for this example comes from the "Twitchell Alfalfa" AmeriFlux eddy
covariance flux tower site in California. The site is located in alfalfa fields and exhibits a mild Mediterranean climate with dry and hot summers, for more information on this site or to download data click `here <https://ameriflux.lbl.gov/sites/siteinfo/US-Tw3>`__. 


Energy Balance Ratio method
---------------------------

The Energy Balance Ratio method (default) is modified from the `FLUXNET
methodology <https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
(step 3 daily heat processing).
The method involves filtering out of extreme values of the daily Energy
Balance Ratio time series, smoothing, and gap filling. Then the inverse of
the filtered and smoothed time series is used as a series of
correction factors for the initial time series of latent energy
(:math:`LE`) and sensible heat (:math:`H`) flux time series.

All steps, abbreviated
^^^^^^^^^^^^^^^^^^^^^^

Below is a step-by-step description of the Energy Balance Ratio
correction routine used by ``flux-data-qaqc``. More details and visual
demonstration of steps are shown below.

**Step 0 (optional):** optionally filter out poor quality data first if quality
control (QC) values or flags are provided with the dataset or other means. For
example, FLUXNET data includes QC values for :math:`H` and :math:`LE`,
e.g. H_F_MDS_QC and LE_F_MDS_QC are QC values for gap filled :math:`H` and
:math:`LE`. This allows for manual pre-QaQc of data.

**Step 1:** calculate the Energy Balance Ratio (EBR =
:math:`\frac{H + LE}{Rn – G}`) daily time series from raw data.

**Step 2:** filter EBR values that are outside 1.5 times the
interquartile range.

**Step 3:** for each day in the daily time series of filtered EBR, a
sliding window of +/- 7 days (15 days) is used to select up to 15
values.

**Step 4:** for each day take a percentile (default 50) of the 15 EBR
values. Check if the inverse of the EBR value is :math:`> |2|` or if the
the inverse of the ratio multiplied by the measured :math:`LE` would
result in a flux greater than 850 or less than -100 :math:`w/m^2`, if so
leave a gap for filling later.

**Step 5:** if less than +/- 5 days exist in the sliding 15 day window,
use the mean EBR for all days in a +/- 5 day (11 day) sliding window.
Apply same criteria for an extreme EBR value as in step 4.

**Step 6:** if no EBR data exist in the +/- 5 sliding window to average,
fill remaining gaps of EBR with the mean from a +/- 5 day sliding window
over the day of year mean for all years on record, i.e. 5 day
climatology. Calculate the 5 day climatology from the filtered and
smoothed EBR as produced from step 5. Apply same criteria for an extreme
EBR value as in steps 4 and 5.

**Step 7:** use the filtered EBR time series from previous steps to
correct :math:`LE` and :math:`H` by multiplying by the energy balance
closure correction factor :math:`{EBC_{CF}} = \frac{1}{EBR}`, where EBR
has been filtered by the previous steps. Use the corrected :math:`LE`
and :math:`H` to calculate the corrected EBR.

**Step 8:** calculate corrected :math:`ET` from corrected :math:`LE`
using average air temperature to adjust the latent heat of vaporization.

**Step 9 (optional):** if desired, fill remaining gaps in the corrected
:math:`ET` time series with :math:`ET` that is calculated by gridMET
reference :math:`ET` (:math:`ETr` or :math:`ETo`) multiplied by the filtered and smoothed
fraction of reference ET (:math:`ETrF` or :math:`EToF`).


Step 0, manual cleaning of poor quality data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below we can see that the daily time series of net radiation (:math:`Rn`) has
some periods of poor quality data. This is a common issue due, e.g. to
instrumentation problems, that cannot always be avoided. In this case the
sensor did not record values at night (or they were not provided with the data)
when :math:`Rn` values are lower for several days (e.g. around 8/26/2014) which
resulted in overestimates of daily mean :math:`Rn` during these periods. Although these days can automatically be filtered out by the :obj:`.QaQc` class, the example below shows a way of manually filtering them because in other cases outliers in the daily data may not be caused by resampling of sub-daily data with systematic measurement gaps. The main point is that manual inspection and potentially pre-filtering of poor quality data before proceeding with energy balance closure corrections is often necessary.

.. raw:: html
    :file: _static/closure_algorithms/step0_NotFiltered.html

There are several ways to conduct manual pre-filtering of poor quality meterological time series data, to filter data based on input quality flags or numeric quality values see :ref:`Quality-based data filtering`. 

``flux-data-qaqc`` also allows for filtering of poor quality data on the fly as shown in this example. In other words, we simply filter out the periods we think have bad data for :math:`Rn` within Python before running the closure correction. After manually determing the date periods with poor quality :math:`Rn`, here is how they were filtered oiut before running the correction:

    >>> import pandas as pd 
    >>> import numpy as np
    >>> from fluxdataqaqc import Data, QaQc
    >>> d = Data('Path/to/config.ini')
    >>> # days with sub daily gaps can be filtered out automatically here, 
    >>> # see "Tip" below the following plot 
    >>> q = QaQc(d, drop_gaps=False) 
    >>> # rename dataframe columns for ease of variable access, adjust
    >>> df = q.df.rename(columns=q.inv_map)

Here were the dates chosen and one way to filter them,

    >>> # make a QC flag column for Rn
    >>> df['Rn_qc'] = 'good'
    >>> df.loc[pd.date_range('2/10/2014','2/10/2014'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('8/25/2014','9/18/2014'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('10/21/2015','10/26/2015'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('10/28/2015','11/1/2015'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('7/23/2016','7/23/2016'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('9/22/2016','9/22/2016'), 'Rn_qc'] = 'bad'
    >>> df.loc[pd.date_range('3/3/2017','3/3/2017'), 'Rn_qc'] = 'bad'
    >>> # filter (make null) based on our QC flag column for Rn
    >>> df.loc[df.Rn_qc == 'bad', 'Rn'] = np.nan
    >>> # reassign to use pre-filtered data for corrections
    >>> q.df = df

The resulting energy balance component plot with :math:`Rn` filtered:

.. raw:: html
    :file: _static/closure_algorithms/step0_Filtered.html

.. tip::
   In this case, the issues with :math:`Rn` were caused by resampling 30 minute
   data with systematic night-time gaps. These sort of issues can be
   automatically handled when creating a :obj:`.QaQc` object; the keyword
   arguments ``drop_gaps`` and ``daily_frac`` to the :obj:`.QaQc` class are
   used to automatically filter out days with measurement gaps of varying size,
   i.e.,
   
   >>> d = Data('path/to/config.ini')
   >>> q = QaQc(d, drop_gaps=True, daily_frac=0.8)
   >>> q.correct_data()

   This would produce very similar energy balance closure results as the manual
   filter above. Another more fine-grained option would have been to flag the
   days with gaps in the sub-daily input time series that you would like to
   filter by :meth:`.Data.apply_qc_flags`.


.. Note::
   The remaining step-by-step explanation in this page uses the pre-filtered
   input time series, however results of the energy balance closure correction
   without pre-filtering outliers of :math:`Rn` are also shown in plots for the
   final steps (8 and 9) for comparison. If you now ran:

   >>> q.df = df
   >>> q.correct_data()
   >>> q.plot(output_type='show')

   This will directly produce the same output of step 9 using the 
   pre-filtered data. 

Steps 1 and 2, filtering outliers of EBR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate daily EBR = :math:`\frac{H + LE}{Rn - G}` time series and
filter out extreme values that are outside 1.5 the interquartile range.
Note, in ``flux-data-qaqc`` this is named as “ebr”.

.. raw:: html
    :file: _static/closure_algorithms/steps1_2_PreFiltered.html

Steps 3, 4, and 5, further filtering of EBR using moving window statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Filter the EBR time series using statistics performed over multiple
moving windows. Specifically, take the median EBR from a +/- 7 day
moving window, if less than 11 days exist in this window take the mean
from a +/- 5 day moving window. In both of these cases check the
resulting value before retaining based on the following criteria:

-  the inverse of the EBR value must be :math:`> |2|`
-  the the inverse of the ratio multiplied by the measured :math:`LE`
   should result in a flux less than 850 and greater than -100
   :math:`w/m^2`

If either of these criteria are not met leave a gap for the day for
filling in later steps.

.. raw:: html
    :file: _static/closure_algorithms/steps3_5_PreFiltered.html

Step 6, calculate the 5 day climatology of EBR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute the 5 day climatology of daily EBR (as adjusted from previous
steps) to fill in remaining gaps of 11 or more days. Specifically,
calculate the the day of year mean of the EBR for all years in record
and then extract the day of year mean using a moving +/- 5 day (11 day)
moving window. The resulting value is also checked against the same
criteria described in steps 3-5:

-  the inverse of the EBR value must be :math:`> |2|`
-  the the inverse of the ratio multiplied by the measured :math:`LE`
   should result in a flux less than 850 and greater than -100
   :math:`w/m^2`

Note, this step is only used for remaining gaps which should be larger
than 11 days in the EBR time series following step 5. This example has a
few time periods that were filled with the 5 day climatology of EBR
which can be seen as the thin blue line in the plot below.

.. raw:: html
    :file: _static/closure_algorithms/step6_PreFiltered.html

``flux-data-qaqc`` also keeps a record of the 5 day climatology of the
Energy Balance Ratio as calculated at this step (shown below), it is 
named by ``flux-data-qaqc`` as ebr_5day_clim.

.. raw:: html
    :file: _static/closure_algorithms/5dayclim_PreFiltered.html

Steps 7 and 8 correct turbulent fluxes, EBR, and ET
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate corrected :math:`LE` and :math:`H` by multiplying by
:math:`\frac{1}{EBR}` where :math:`EBR` is the filtered EBR time series
from previous steps:

.. math:: LE_{corr} = LE \times \frac{1}{EBR}

\ and

.. math:: H_{corr} = H \times \frac{1}{EBR}.

Use corrected LE and H to calculate the corrected EBR,

.. math:: EBR_{corr} = \frac{H_{corr} + LE_{corr}}{Rn - G}.

Calculate ET from LE using average air temperature to adjust the latent
heat of vaporization following the method of Harrison, L.P. 1963,

.. math:: ET_{mm \cdot day^{-1}} = 86400_{sec \cdot day^{-1}} \times \frac{LE_{w \cdot m^{-2}}}{2501000_{MJ \cdot kg^{-1}} - (2361 \cdot T_{C})}, 

where evapotransipiration (:math:`ET`) in :math:`mm \cdot day^{-1}`,
:math:`LE` is latent energy flux in :math:`w \cdot m^{-2}`, and
:math:`T` is air temperature in degrees celcius. The same approach is
used to calculate corrected :math:`ET` (:math:`ET_{corr}`) using
:math:`LE_{corr}`.

The plot below shows the time series of the initial and corrected ET (:math:`ET` and :math:`ET_{corr}`).

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_PreFiltered.html

There were not significant gaps in the energy balance components for this dataset and therefore step 9 was not used, although it is still demonstrated with an artificial gap in the next step. 

The following plot shows the energy balance closure of the initial and corrected data after applying the steps above, including the manual pre-filtering of :math:`Rn`,

.. raw:: html
    :file: _static/closure_algorithms/EBC_scatter_PreFiltered.html

Notice the mean daily corrected energy balance ratio (slope of corrected) is 1 or near perfect closure. However, the same plot below shows the results if we skipped the manual pre-filtering of outlier :math:`Rn` values. In this case the resulting corrected mean closure is only 0.93:

.. raw:: html
    :file: _static/closure_algorithms/EBC_scatter_noPreFilter.html

.. Tip:: 
   These and other interactive visualizations of energy balance closure results 
   are provided by default via the :meth:`.QaQc.plot` method.

In ``flux-data-qaqc`` new variable names from these steps are: LE_corr, H_corr,
ebr, ebr_corr, ebc_cf, ET, ET_corr, ebr_corr, and ebr_5day_clim. The inverse of
the corrected EBR (filtered from previous steps) is named ebc_cf which is short
for energy balance closure correction factor as described by the `FLUXNET
methodology
<https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/data-processing/>`__
(step 3 daily heat processing).

Step 9, optionally gap fill corrected ET using gridMET reference ET and reference ET fraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is done by downloading :math:`ETr` or :math:`ETo` (default is :math:`ETr`)
for the overlapping gridMET cell (site must be in CONUS) and then calculating,

.. math:: ET_{fill} = ETrF \times ET_r,

\ where

.. math:: ETrF = \frac{ET_{corr}}{ET_r}

:math:`ET_{corr}` is the corrected ET produced by step 8 and :math:`ETrF` is
the fraction of reference ET. :math:`ETrF` if first filtered to remove
outliers outside of 1.5 times the interquartile range, it is then smoothed with
a 7 day moving average (minimum of 2 days must exist in window) and lastly it
is linearly interpolated to fill any remaining gaps. 

The same gap filling procedure can easily be done using gridMET grass reference ET (:math:`ETo`) as opposed to alfalfa reference ET (:math:`ETr`).

.. Tip:: 
   The filtered and raw versions of :math:`ETrF`/:math:`EToF`, gridMET
   :math:`ETr`, gridMET :math:`ETo`, gap days, and monthly total number of gap
   filled days are tracked for post-processing and visualized by the
   :meth:`.QaQc.plot` and :meth:`.QaQc.write` methods.

Since the data used in this example does not have gaps, for illustration we have created the following large gap in the measured energy balance components from May through August, 2014:

.. raw:: html
    :file: _static/closure_algorithms/step0_Filtered_withGap.html

The resulting time series of :math:`ET_{corr}` using the optional gap filling method described is shown below.  

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_with_gapFill.html

Note, the gap filled values of :math:`ET` (green line) do not accurately catch the harvesting cycles of alfalfa however the :math:`ET_{corr}` values (blue line) do, this is because the gap filled values are based from gridMET reference ET which is not locally representative. If this is hard to see, try using the box zoom tool on the right of the plot to zoom in on the gap-filled period.

This ET gap-filling step is used by default when running ``flux-data-qaqc``
energy balance closure correction routines, to disable it set the
``etr_gap_fill`` argument of :meth:`QaQc.correct_data` to False, e.g.

.. code-block:: python

    # q is a QaQc instance
    q.correct_data(meth='ebr', etr_gap_fill=False)


In ``flux-data-qaqc`` new variable names from this step are: ETrF,
ETrF_filtered, gridMET_ETr, ET_gap, ET_fill, and ET_fill_val. The
difference between ET_fill and ET_fill_val is that the latter is masked
(null) on days that the fill value was not used to fill gaps in
:math:`ET_{corr}`. Also, ET_gap is a daily series of True and False
values indicating which days (from step 8) of :math:`ET_{corr}` were gaps
that were subsequently filled.

.. Note::
   When using the :math:`ETr`-based gap-filling option, any gap filled days
   will also be used to fill in gaps of :math:`LE_{corr}`, therefore the mean
   closure as found in the daily and monthly closure scatter plot outputs (from
   :meth:`.QaQc.plot`) will be updated to reflect the influence of the
   gap-filled days.

Bowen Ratio method
------------------

The Bowen Ratio energy balance closure correction method implemented
here follows the typical approach where the corrected latent energy
(:math:`LE`) and sensible heat (:math:`H`) fluxes are adjusted the
following way

.. math::  LE_{corr} = \frac{(Rn - G)}{(1 + \beta)}, 

\ and

.. math::  H_{corr} = LE_{corr} \times \beta 

where :math:`\beta` is the Bowen Ratio, the ratio of sensible heat flux
to latent energy flux,

.. math:: \beta = \frac{H}{LE}.

This routine forces energy balance closure for each day in the time
series.

Here is the resulting :math:`ET_{corr}` time series using the pre-filtered (:math:`Rn`) energy balance time series and the Bowen Ratio method:

.. raw:: html
    :file: _static/closure_algorithms/ET_ts_BR_PreFiltered.html

And here is the energy balance closure scatter plot which shows the forced closure of the method:

.. raw:: html
    :file: _static/closure_algorithms/EBC_BR_scatter_PreFiltered.html

New variables produced by ``flux-data-qaqc`` by this method include: br
(Bowen Ratio), ebr, ebr_corr, LE_corr, H_corr, ET, ET_corr, energy, flux, and
flux_corr.


.. include:: ../../CHANGELOG.rst
API Reference
=============

This page documents objects and functions provided by ``flux-data-qaqc``.



Data
----

.. autoclass:: fluxdataqaqc.Data
    :members:
    :undoc-members:
    :show-inheritance:

QaQc
----

.. autoclass:: fluxdataqaqc.QaQc
    :members:
    :undoc-members:
    :show-inheritance:

Plot
----

.. autoclass:: fluxdataqaqc.Plot
    :members:
    :undoc-members:
    :show-inheritance:

utility classes and functions
-----------------------------

.. automodule:: fluxdataqaqc.util
    :members:
    :undoc-members:
    :show-inheritance:

