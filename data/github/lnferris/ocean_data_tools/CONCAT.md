# ocean_data_tools: a MATLAB toolbox for interacting with bulk freely-available oceanographic data

<img src="https://user-images.githubusercontent.com/24570061/85356569-4664ba80-b4dd-11ea-9ec7-8ec26df76dcf.png" width="500">

[![GitHub license](https://img.shields.io/github/license/lnferris/ocean_data_tools)](https://github.com/lnferris/ocean_data_tools/blob/master/LICENSE) [![GitHub stars](https://img.shields.io/github/stars/lnferris/ocean_data_tools)](https://github.com/lnferris/ocean_data_tools/stargazers) [![GitHub forks](https://img.shields.io/github/forks/lnferris/ocean_data_tools)](https://github.com/lnferris/ocean_data_tools/network) 
[![View ocean_data_tools on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/80047-ocean_data_tools) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02497/status.svg)](https://doi.org/10.21105/joss.02497)

**Copyright (c) 2020 lnferris** 

ocean_data_tools simplifies the process of extracting, formatting, and visualizing freely-available oceanographic data. While a wealth of oceanographic data is accessible online, some end-users may be dissuaded from utilizing this data due to the overhead associated with obtaining and formatting it into usable data structures. ocean_data_tools solves this problem by allowing the user to transform common oceanographic data sources into uniform structs, call generalized functions on these structs, easily perform custom calculations, and make graphics.

Find a bug, have a question, or want to chat about contributing? Open an issue or email lnferris@alum.mit.edu.

### [Getting Started](#getting-started-1)

### [Dependencies](#dependencies-1)

### [Accessing Help](#accessing-help-1)

### [How to Contribute](#how-to-contribute-1)

### [Contents](#contents-1)

### [Finding Data](#finding-data-1)

### [Citing ODT](#citing-odt-1)

## Getting Started

1. Download [bathymetry](#bathymetry).
2. Download [nctoolbox](https://github.com/nctoolbox/nctoolbox). You will need to run the command ``setup_nctoolbox`` at the beginning of each MATLAB session.
3. Add ocean_data_tools and nctoolbox to the path. Specifically, the following folders must be added to the [path](https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html):

- ocean_data_tools/ocean_data_tools
- ocean_data_tools/ocean_data_tools/utilities
- nctoolbox/

4. Run each demonstration in **demos/demos.m**, which contains example usages for all functions. All required test data is included in **data/**.

Functions are named using a two-part system. The prefix (``argo_``, ``bathymetry_``, ``general_``, etc.) indicates the appropriate data source, while the suffix (``\_build``, ``\_profiles``, ``\_section``, etc.) indicates the action performed. Functions with the ``\_build`` suffix load raw data into uniform structs (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``). Uniform structs created by ``\_build`` functions are compatable with any ``general_`` function.

Data sources currently supported:
| Data Source | DOI, Product Code, or Link    |
|:--  |:--|
| Argo floats | [doi:10.17882/42182](https://doi.org/10.17882/42182) |
| Smith & Sandwell bathymetry | [doi:10.1126/science.277.5334.1956](https://doi.org/10.1126/science.277.5334.1956) |
| IOOS Glider DAC | https://gliders.ioos.us/ |
| MOCHA Climatology | [doi:10.7282/T3XW4N4M](https://doi.org/10.7282/T3XW4N4M) |
| HYbrid Coordinate Ocean Model | https://hycom.org |
| CMEMS Global Ocean 1/12° Physics Analysis and Forecast | GLOBAL_ANALYSIS_FORECAST_ PHY_001_024 |
| CMEMS Global Ocean Waves Multi Year | GLOBAL_REANALYSIS_WAV_001_032 |
| GO-SHIP hydrographic cruises | https://www.go-ship.org/ |
| World Ocean Atlas 2018 | https://www.ncei.noaa.gov/products/world-ocean-atlas |
| World Ocean Database | https://www.ncei.noaa.gov/products/world-ocean-database |

Main functions are located in **ocean_data_tools/**. Demonstrations are located in **demos/**. Test datas are located in **data/**. Shell scripts for batch downloading data are located in **shell_scripts/**. While shell scripts can be run directly in a macOS Terminal, running them in Windows requires [Cygwin](https://www.cygwin.com/) (and perhaps slight modification of commands). Python syntax examples are located in **python/**, which may be grow to become a module in the future.

## Dependencies

The only true dependency is [nctoolbox](https://github.com/nctoolbox/nctoolbox).

It is recommended to also download [Gibbs-SeaWater (GSW) Oceanographic Toolbox](http://www.teos-10.org/software.htm#1). A benefit of ocean_data_tools is that neatly packs data into uniform structs; at which point a user can easily apply custom calculations or functions from other toolboxes such as GSW. See an [example](docs/gsw_example.md).

## Accessing Help

To access help, run the command ``doc ocean_data_tools``.

## How to Contribute

* Want to make changes or add a new function? (1) Fork the repository (make your own separate copy), (2) make changes, and (3) open a 'pull request'. Once approved, it can be merged into the master branch. If you wish to chat beforehand about your contribution, open an issue or email lnferris@alum.mit.edu. 
* Don't use git often and don't want to remember all the terminal commands? Download [GitHub Desktop](https://desktop.github.com/).
* Find a bug in the code? Open an 'issue' to notify contributors and create an official record.

Before contributing, please see [Contents](#contents-1) and consider how your function fits into ocean_data_tools and its ethos of structure arrays. At a minimum, functions must be well-documented and address a specific freely-available oceanographic data source which can be accessed by anyone online.

Adding a new function isn't the only way to contribute. Python, Julia, etc. translations of existing Matlab functions are also welcomed!

If you are interested in becoming a formal collaborator (e.g. have direct access and co-manage this repository), please reach out.

## Contents

#### [Building uniform structs from data sources](#building-uniform-structs-from-data-sources-1)

#### [General functions for subsetting and plotting uniform structs](#general-functions-for-subsetting-and-plotting-uniform-structs-1)

#### [Plotting gridded data without building structs](#plotting-gridded-data-without-building-structs-1)

#### [Adding bathymetry to existing plots](#adding-bathymetry-to-existing-plots-1)

#### [Additional functions for inspecting Argo data](#additional-functions-for-inspecting-argo-data-1)

#### [Miscellaneous utilities](#miscellaneous-utilities-1)


### Building uniform structs from data sources

**[argo_build](docs/argo_build.md)** searches the locally-stored Argo profiles matching the specified region & time period and builds a uniform struct

**[glider_build](docs/glider_build.md)** loads an archived glider survey (downloaded from gliders.ioos.us/erddap) and builds a uniform struct

**[mocha_build_profiles](docs/mocha_build_profiles.md)** builds a uniform struct of profiles from the MOCHA Mid-Atlantic Bight climatology

**[model_build_profiles](docs/model_build_profiles.md)**  builds a uniform struct of profiles from HYCOM or Operational Mercator CMEMS GLOBAL_ANALYSIS_FORECAST_PHY_001_024

<img src="https://user-images.githubusercontent.com/24570061/88250150-ac776580-cc74-11ea-8c72-cea7cc50b4d9.png" width="700">

**waves_build** builds a uniform struct of timeseries from CMEMS Global Ocean Waves Multi Year product GLOBAL_REANALYSIS_WAV_001_032

**[whp_cruise_build](docs/whp_cruise_build.md)** builds a uniform struct of profiles from GO-SHIP cruise data in WHP-Exchange Format

**[woa_build_profiles](docs/woa_build_profiles.md)** builds a uniform struct of profiles from World Ocean Atlas 2018 Statistical Mean for All Decades, Objectively Analyzed Mean Fields

**[wod_build](docs/wod_build.md)** builds a uniform struct of profiles from World Ocean Database data

*Don't see a function yet for your preferred data source? Email lnferris@alum.mit.edu to request or [contribute](#how-to-contribute-1).*

### General functions for subsetting and plotting uniform structs

**[general_depth_subset](docs/general_depth_subset.md)** subsets a uniform struct by depth

**[general_map](docs/general_map.md)** plots coordinate locations in a uniform struct, with optional bathymetry contours

**[general_profiles](docs/general_profiles.md)** plots vertical profiles in a uniform struct

**[general_region_subset](docs/general_region_subset.md)** subsets a uniform struct by polygon region

<img src="https://user-images.githubusercontent.com/24570061/88250944-358f9c00-cc77-11ea-9b0d-2d582ad186dd.png" width="700">

**[general_remove_duplicates](docs/general_remove_duplicates.md)** removes spatially (or spatially and temporally) non-unique profiles from a uniform struct

**[general_section](docs/general_section.md)** plots a data section from a uniform struct


### Plotting gridded data without building structs

**[mocha_domain_plot](docs/mocha_domain_plot.md)** plots a 3-D domain from the MOCHA Mid-Atlantic Bight climatology

**[mocha_simple_plot](docs/mocha_simple_plot.md)** plots a 2-D layer from the MOCHA Mid-Atlantic Bight climatology

**[model_domain_plot](docs/model_domain_plot.md)** plots a 3-D domain from HYCOM or Operational Mercator CMEMS GLOBAL_ANALYSIS_FORECAST_PHY_001_024

**[model_simple_plot](docs/model_simple_plot.md)** plots a 2-D layer from HYCOM or Operational Mercator CMEMS GLOBAL_ANALYSIS_FORECAST_PHY_001_024

<img src="https://user-images.githubusercontent.com/24570061/88250403-8900ea80-cc75-11ea-8a5d-8a474d2e5c3f.png" width="700">

**[woa_domain_plot](docs/woa_domain_plot.md)** plots a 3-D domain from World Ocean Atlas 2018 Statistical Mean for All Decades, Objectively Analyzed Mean Fields

**[woa_simple_plot](docs/woa_simple_plot.md)** plots a 2-D layer from World Ocean Atlas 2018 Statistical Mean for All Decades, Objectively Analyzed Mean Fields

### Adding bathymetry to existing plots

**[bathymetry_extract](docs/bathymetry_extract.md)** extracts a region of Smith & Sandwell Global Topography and outputs as arrays

**[bathymetry_plot](docs/bathymetry_plot.md)** adds bathymetry to 2-D (latitude vs. longitude) or 3-D (latitude vs. longitude vs. depth) data plots

<img src="https://user-images.githubusercontent.com/24570061/88251161-ed24ae00-cc77-11ea-87d6-0e3b4484764d.jpg" width="700">

**[bounding_region](docs/bounding_region.md)** finds the rectangular region around a uniform struct and/or list of coordinates to pass as an argument for other bathymetry functions

**[bathymetry_section](docs/bathymetry_section.md)** adds Smith & Sandwell Global Topography to a section from plot using bathymetry data nearest to specified coordinates

<img src="https://user-images.githubusercontent.com/24570061/88250660-3d027580-cc76-11ea-808c-f51d5105e420.png" width="700">

### Additional functions for inspecting Argo data

**[argo_platform_map](docs/argo_platform_map.md)** plots locations of Argo profiles in a uniform struct, coloring markers by platform (individual Argo float)

<img src="https://user-images.githubusercontent.com/24570061/88250439-a2099b80-cc75-11ea-9516-ad3d1f65fdf9.jpg" width="700">

**[argo_platform_subset](docs/argo_platform_subset.md)** subsets a uniform struct of Argo data to one platform (individual Argo float)

**[argo_profiles_map](docs/argo_profiles_map.md)** plots coordinate locations of Argo profiles in uniform struct argo, using colors corresponding to argo_profiles called on the same struct

**[argo_profiles](docs/argo_profiles.md)** plots vertical Argo profiles in uniform struct argo, using colors corresponding to argo_profiles_map called on the same struct


### Miscellaneous utilities

**[region_select](docs/region_select.md)** creates coordinate list (which represents vertices of a polygon region) by clicking stations on a plot

**[transect_select](docs/transect_select.md)** creates a coordinate list (which represents a virtual transect) by clicking stations on a plot

<img src="https://user-images.githubusercontent.com/24570061/88250639-2b20d280-cc76-11ea-9c94-3ce16300f735.png" width="700">


## Finding Data

There two types of datasets: those that need to be downloaded manually<sup>1</sup> and those that can be accessed remotely<sup>2</sup> through OpenDAP (e.g. the data can be accessed directly on the the internet using a url). 

#### argo<sup>1</sup>

Download [Argo data](https://argo.ucsd.edu/) directly from GDAC FTP servers using either the [Coriolis selection tool](http://www.argodatamgt.org/Access-to-data/Argo-data-selection), or the [US GDAC](https://nrlgodae1.nrlmry.navy.mil/cgi-bin/argo_select.pl). See the [Argo User's Manual](http://www.argodatamgt.org/Documentation) for more information.

Alternatively run **shell_scripts/download_argo** to download data via File Transfer Protocol.

#### bathymetry<sup>1</sup>

To get bathymetry data (for ``bathymetry_dir``), download Smith & Sandwell under [Global Topography V19.1](https://topex.ucsd.edu/marine_topo/) in netcdf form (topo_20.1.nc).

#### glider<sup>1</sup>

Vist [gliders.ioos.us/erddap](https://gliders.ioos.us/erddap/index.html). Click "View a List of All 779 Datasets" or use the "Advanced Search". After choosing a dataset, navigate to the [Data Access Form](https://gliders.ioos.us/erddap/tabledap/ce_311-20170725T1930-delayed.html). To get started, select these variables:

<img src="https://user-images.githubusercontent.com/24570061/94058620-419af580-fdaf-11ea-859a-616c8b5b1433.png" width="700">

Scroll to "File type:". In the drop-down menu, select ".nc". Click "Submit".

#### mocha<sup>2</sup>

The url for MOCHA Mid-Atlantic Bight climatology is embedded. See [Rutgers Marine catalog](http://tds.marine.rutgers.edu/thredds/catalog.html).

#### model<sup>1,2</sup>

HYCOM data may be accessed remotely using OpenDAP. Get the data url by visiting the [HYCOM website](https://www.hycom.org/dataserver/gofs-3pt1/analysis). For example, click Access Data Here -> GLBv0.08/expt_57.7 (Jun-01-2017 to Sep-30-2017)/ -> Hindcast Data: Jun-01-2017 to Sep-30-2017. Click on the OpenDAP link. Copy the url as and use this as the ``source`` in ``model_build_profiles``.

Alteratively, download subsetted HYCOM data using NCSS. Get the data url by visiting the [HYCOM website](https://www.hycom.org/dataserver/gofs-3pt1/analysis). For example, click Access Data Here -> GLBv0.08/expt_57.7 (Jun-01-2017 to Sep-30-2017)/ -> Hindcast Data: Jun-01-2017 to Sep-30-2017. Click on the NetcdfSubset link. Set constraints and copy the NCSS Request URL at the bottom of the page. Run **shell_scripts/download_hycom_lite**. To download multiple months or years, run **shell_scripts/download_hycom_bulk_daily** (partition files by day) or **shell_scripts/download_hycom_bulk_monthly** (partition files by month). Please use responsibly.

For Mercator, download Copernicus Marine data directly from FTP servers. First make a [Copernicus account](http://marine.copernicus.eu/services-portfolio/access-to-products/). Use the selection tool to download GLOBAL_ANALYSIS_FORECAST_PHY_001_024. Alternatively run **shell_scripts/download_mercator**. Before running the script, follow the instructions for modifying your ~/.netrc file in the comments of the script.

#### waves<sup>1</sup>

First make a [Copernicus account](http://marine.copernicus.eu/services-portfolio/access-to-products/). Use the selection tool to download CMEMS Global Ocean Waves Multi Year product GLOBAL_REANALYSIS_WAV_001_032.

#### whp_cruise<sup>1</sup>

For [GO-SHIP data](https://usgoship.ucsd.edu/hydromap/), get CTD data (for ``ctdo_dir``) by choosing a [GO-SHIP cruise](https://cchdo.ucsd.edu/search?q=GO-SHIP) and downloading the CTD data in whp_netcdf format. More information about whp_netcdf parameters is available [here](https://exchange-format.readthedocs.io/en/latest/index.html#). Get LADCP data (for ``uv_dir``, ``wke_dir``) [here](https://currents.soest.hawaii.edu/go-ship/ladcp/). There is information about LACDP processing [here](https://www.ldeo.columbia.edu/~ant/LADCP.html).

#### woa<sup>2</sup>

Functions build the World Ocean Atlas url at maximum resolution based on arguments, but coarser resolutions and seasonal climatologies are available at the [NODC website](https://www.nodc.noaa.gov/OC5/woa18/woa18data.html). Note NCEI is scheduled to update data urls in the near future. Functions will be updated as such.

#### wod<sup>1</sup>

Search the [World Ocean Database](https://www.nodc.noaa.gov/OC5/SELECT/dbsearch/dbsearch.html) and select products.

## Citing ODT

Ferris, L., (2020).  ocean_data_tools:  A MATLAB toolbox for interacting with bulk freely-available oceanographic data. Journal of Open Source Software, 5(54), 2497. https://doi.org/10.21105/joss.02497

---
title: 'ocean_data_tools: A MATLAB toolbox for interacting with bulk freely-available oceanographic data'
tags:
  - MATLAB
  - oceanography
authors:
  - name: Laur Ferris
    orcid: 0000-0001-6446-9340
    affiliation: 1
affiliations:
 - name: Virginia Institute of Marine Science
   index: 1
date: 03 July 2020
bibliography: paper.bib
---

# Statement of Need

``ocean_data_tools`` simplifies the process of extracting, formatting, and 
visualizing freely-available oceanographic data. A wealth of oceanographic 
data (from research cruises, autonomous floats, global ocean models, etc.)
is accessible online. However, many oceanographers and environmental 
scientists (particularly those from subdisciplines not accustomed to working
with large datasets) can be dissuaded from utilizing this data because of the
overhead associated with determining how to batch download data and 
format it into easily-manipulable data structures. ``ocean_data_tools``
solves this problem by allowing the user to transform common oceanographic 
data sources into uniform structure arrays, call general functions on these structure arrays, 
perform custom calculations, and make graphics. 

![Building a virtual cruise from the Operational Mercator global ocean
analysis and forecast system at 1/12 degree with 3D bathymetry [@Smith:1997]. 
Showing (a) a 3D velocity plot created using ``model_domain_plot``, (b) 
virtual cruise selection using ``transect_select``, and ``model_build_profiles``, 
(c) coordinates of the resulting uniform structure array, and (d) a temperature section 
plotted using ``general_section`` with ``bathymetry_section``. Three of the 
subplots use colormaps from cmocean [@Thyng:2016]. \label{fig:1}](figure.png)

# Summary

Structure arrays, the common currency of ``ocean_data_tools``, are more user-friendly than the native data storage underlying many of the datasets because they allow the user to neatly group related data of any type or size into containers called fields. Both the structure array and its fields are mutable, and data is directly visible and accessible in the Matlab workspace (unlike NetCDF which requires a function call to read variables).
Matlab was chosen as the language of choice for this toolbox because it is already extensively used within the oceanographic community.
It is also a primary language for much of the community, which is important because this toolbox aims to lower the barrier to entry for using the growing variety of freely-available field- and model-derived oceanographic datasets.

The workflow of ``ocean_data_tools`` is to build uniform structure arrays (e.g. ``argo``,
``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``) from raw datasets and 
call general functions on these structure arrays to map, subset, or plot. Functions with 
the ``\_build`` suffix load raw data into uniform structure arrays. Structure arrays are 
compatible with all ``general_`` functions, and serve to neatly contain the data for use with 
custom user-defined calculations or other toolboxes such as the commonly-used
Gibbs-SeaWater (GSW) Oceanographic Toolbox [@McDougall:2011]. One application of the ``\_build`` 
feature is to create virtual cruises from model output \autoref{fig:1}. The user
draws transects on a map (or passes coordinates as an argument) to build vertical profiles 
from model data. This may be used as a cruise planning tool, to facilitate 
comparison of observations with model output, or to support decision-making in underwater glider 
piloting (using model forecasts to inform ballasting or adjust flight for ocean currents).  Some ``ocean_data_tools`` functions
employ ``nctoolbox`` [@nctoolbox].

| Data Source | DOI, Product Code, or Link    |
|:--  |:--|
| Argo floats | [doi:10.17882/42182](https://doi.org/10.17882/42182) |
| Smith & Sandwell bathymetry | [doi:10.1126/science.277.5334.1956](https://doi.org/10.1126/science.277.5334.1956) |
| IOOS Glider DAC | https://gliders.ioos.us/ |
| MOCHA Climatology | [doi:10.7282/T3XW4N4M](https://doi.org/10.7282/T3XW4N4M) |
| HYbrid Coordinate Ocean Model | https://hycom.org |
| CMEMS Global Ocean 1/12° Physics Analysis and Forecast | GLOBAL_ANALYSIS_FORECAST_ PHY_001_024 |
| GO-SHIP hydrographic cruises | https://www.go-ship.org/ |
| World Ocean Atlas 2018 | https://www.ncei.noaa.gov/products/world-ocean-atlas |
| World Ocean Database | https://www.ncei.noaa.gov/products/world-ocean-database |

: Current ocean_data_tools data sources. \label{table:1}

There are several high-quality ocean and/or climate related Matlab toolboxes such as Climate Data Toolbox for Matlab [@Greene:2019], those part of [SEA-MAT: Matlab Tools for Oceanographic Analysis](https://sea-mat.github.io/sea-mat/), and Gibbs-SeaWater (GSW) Oceanographic Toolbox [@McDougall:2011]. However, there are no other documented and designed-to-be-shared toolboxes filling the same data exploration niche as this one. ``ocean_data_tools`` is unique in encouraging the user to invoke a variety of freely-available data into their exploration and does not expect the user to provide privately-collected measurements or privately-generated model output. It connects users to specific, well-documented data sources Table \ref{table:1}. ``ocean_data_tools`` has already been used for data exploration in support of scientific publications [@Bemis:2020; @Crear:2020]. This toolbox is built for extensibility; the objective is to welcome contributors and continuously add support for additional datasets such as [Remote Sensing 
Systems](http://www.remss.com/) products and European Centre for Medium-Range 
Weather Forecasts (ECMWF) products. The source code for ``ocean_data_tools`` has
been archived to Zenodo with the linked DOI: [@Ferris:2020].

# Acknowledgements

The Virginia Institute of Marine Science (VIMS) provided financial support for this project.
I am grateful to Donglai Gong for ongoing mentorship. I thank the many organizations providing freely-available
data to the oceanography community including (but not limited to) Argo, the HYCOM 
consortium, the Copernicus Programme, the International Global Ship-based Hydrographic
Investigations Program (GO-SHIP), and the National Oceanic and Atmospheric 
Administration (NOAA). I also thank the two reviewers for helpful feedback, especially 
Kelly Kearney for her insightful suggestions. The figure was generated using E.U. Copernicus Marine Service Information. This paper is Contribution No.3960 of the Virginia Institute of Marine Science, William & Mary.

# References
### whp_cruise_build

#### Syntax

```Matlab
[cruise] = whp_cruise_build(ctdo_dir,uv_dir,wvke_dir,variable_list)
```
#### Description

``[cruise] = whp_cruise_build(ctdo_dir,uv_dir,wvke_dir,variable_list)`` searches pathways ``ctdo_dir``, ``uv_dir``, ``wvke_dir`` for CTD data in whp_netcdf format, horizontal LADCP data in netcdf format, and vertical LACDP data in netcdf format respectively. Variable lists for LADCP are fixed, while the CTD variable list is specified using ``variable_list`` (station, woce_date, longitude, latitude, and pressure are included automatically.) Lat/lon information (metadata) is pulled from the CTD files by default. If CTD is not found, metadata from LACDP files are used instead.

The paths used as arguments should point to data from the *same oceanographic cruise*.

``ctdo_dir`` is a character array search path with wildcards. The search path should be the path to the CTD netcdf files (in whp_netcdf format) themselves, not their directory.

``variable_list`` is a cell array where each element is the string name of a variable to be included from CTD files.

``uv_dir`` is a character array search path with wildcards. The search path should be the path to the horizontal LADCP data netcdf files themselves, not their directory.

``wvke_dir`` is a character array path to all files in the directory. 

Example paths: 

```Matlab
ctdo_dir = '/Users/lnferris/Documents/S14/ctd/*.nc'; 
uv_dir = '/Users/lnferris/Documents/S14/whp_cruise/uv/*.nc';
wvke_dir = '/Users/lnferris/Documents/S14/whp_cruise/wvke/';
```

#### Example 1


```Matlab

% Paths to data:

ctdo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/whp_cruise/ctd/*.nc'; % included
uv_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/whp_cruise/uv/*.nc'; % included
wvke_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/whp_cruise/wvke/'; % included
bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc'; % need to download

% Get information about available CTD+ variables:

listing = dir(ctdo_dir);
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

variable_list = {'salinity','temperature','oxygen'};

% Build a uniform struct of cruise data:

[cruise] = whp_cruise_build(ctdo_dir,uv_dir,wvke_dir,variable_list); % Use a dummy path (e.g. uv_dir ='null') if missing data. 

% Map cruise stations:

general_map(cruise,bathymetry_dir,'2Dcontour')

```
<img src="https://user-images.githubusercontent.com/24570061/88341972-89989000-cd0c-11ea-8961-70c76fc639d8.png" width="700">


```Matlab
% Plot a salinity section:

variable = 'salinity'; % See cruise for options.
xref = 'lon'; % See cruise for options.
zref = 'pressure'; % See cruise for options.
general_section(cruise,variable,xref,zref)

```
<img src="https://user-images.githubusercontent.com/24570061/88346894-7939e280-cd17-11ea-96a2-735db5a527d5.png" width="800">


[Back](https://github.com/lnferris/ocean_data_tools#building-uniform-structs-from-data-sources-1)

### general_map

#### Syntax

```Matlab
general_map(object)
general_map(object,bathymetry_dir)
general_map(object,bathymetry_dir,ptype)
```
#### Description

``general_map(object)`` plots coordinate locations (``object.lon`` and ``object.lat``); where ``object`` is a struct created by any of the ``_build`` functions in ocean_data_tools (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``). 

``general_map(object,bathymetry_dir)`` adds bathymetry contours from Smith & Sandwell Global Topography with path ``bathymetry_dir``.

``general_map(object,bathymetry_dir,ptype)`` allows the user to modify plot type from the default contours e.g. ``ptype = '2Dscatter'`` or ``'2Dcontour'``

``bathymetry_dir`` is the character array path to the Smith & Sandwell Global Topography file "topo_20.1.nc"

#### Example 1


```Matlab

% Get variable information:

argo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/argo/*profiles*.nc';
listing = dir(argo_dir); 
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '01-Nov-2015 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Make a map:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
ptype = '2Dcontour'; % '2Dscatter' '2Dcontour'
general_map(argo,bathymetry_dir,ptype) % bathymetry_dir, ptype optional

```
<img src="https://user-images.githubusercontent.com/24570061/88301724-fd1dab80-ccd2-11ea-9ea7-7badf1424865.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#general-functions-for-subsetting-and-plotting-uniform-structs-1)

### general_profiles

#### Syntax

```Matlab
general_profiles(object,variable,zref)
```
#### Description

``general_profiles(object,variable,zref)`` plots vertical profiles of the specified ``variable`` in struct ``object`` as a function of the depth field specified by ``zref``; where ``object`` is a struct created by any of the ``_build`` functions in ocean_data_tools (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``) and ``variable`` is a field name.

``variable`` is the string name of the field (of ``object``) to be plotted as the x-variable of the vertical profile

``zref`` is the string name of the field (of ``object``) to be plotted as the depth variable of the vertical profile

#### Example 1


```Matlab

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '01-Nov-2015 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Plot profiles:

object = argo;  % argo, cruise, hycom, mercator, woa, wod
variable = 'TEMP_ADJUSTED'; % see particular struct for options
zref = 'depth'; % see particular struct for options
general_profiles(object,variable,zref)

```
<img src="https://user-images.githubusercontent.com/24570061/88301788-11fa3f00-ccd3-11ea-9cdf-1622f701bfe9.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#general-functions-for-subsetting-and-plotting-uniform-structs-1)

### bathymetry_section

#### Syntax

```Matlab
[bath_section,lon_section,lat_section,time_section] = bathymetry_section(bathy,xcoords,ycoords,xref)
[bath_section,lon_section,lat_section,time_section] = bathymetry_section(bathy,xcoords,ycoords,xref,filled)
```
#### Description

``[bath_section,lon_section,lat_section] = bathymetry_section(bathy,xcoords,ycoords,xref)`` makes a section plot from ``bathy``, where ``bathy`` is a struct of Smith & Sandwell Global Topography created using ``bathymetry_extract``. ``xcoords`` (longitude) and ``ycoords`` (latitude) are densified to a 1/60-deg grid before bathymetry is interpolated. The bathymetry section is plotted against ``xref``; where ``xref = 'lon'``, ``'lat'``,``'km'``, or a time vector of length(xcoords). The extracted data is output ``bath_section``, ``lon_section``, ``lat_section``, and ``time_section``; output vectors are sorted by the selected reference axis (longitude, latitude, or time).
 
``[bath_section,lon_section,lat_section,time_section] = bathymetry_section(bathy,xcoords,ycoords,xref,filled)`` allows the bathymetry to be filled in black down to the x-axis (instead of a simple line). Set ``filled=1`` to turn on, ``filled=0`` to turn off.

``xcoords`` and ``ycoords`` are vectors of coordinates. Rows or columns are fine, and both -180/180 or 0/360 notation are fine.

When ``xref`` is a time vector, it must be of ``length(xcoords)`` and elements of the vector must be datenums. Otherwise set ``xref = 'lon'`` or  ``xref = 'lat'``. Alteratively pass ``xref = 'km'`` to plot in along-track distance, assuming spherical earth.  

#### Example 1

```Matlab

% Add bathymetry to a temperature section plot from the list of coordinates stored in struct cruise:

xref = 'lon'; 
general_section(cruise,'temperature',xref,'pressure') % plot temperature section
xcoords = cruise.lon; 
ycoords = cruise.lat;
filled = 1;
[bathy] = bathymetry_extract(bathymetry_dir,bounding_region(cruise));
bathymetry_section(bathy,xcoords,ycoords,xref,filled)
```
<img src="https://user-images.githubusercontent.com/24570061/88436173-b8c50500-cdd1-11ea-8270-22930d42843c.png" width="800">

#### Example 2
```Matlab
% Plot bathymetry nearest to a list of coordinates. Use latitude as the x-axis:

xref = 'lat'; 
xcoords = [60 60.1 60.4 60.2 59.9]; 
ycoords = [10 20.1 15.0 16.1 13.7]; 
[bathy] = bathymetry_extract(bathymetry_dir,bounding_region([],xcoords,ycoords));
figure
bathymetry_section(bathy,xcoords,ycoords,xref)
```
#### Example 3
```Matlab
% Plot bathymetry nearest to a list of coordinates. Use a time as the x-axis:

xref = [737009 737010 737011 737012 737013]; 
xcoords = [60 60.1 60.4 60.2 59.9]; 
ycoords = [10 20.1 15.0 16.1 13.7]; 
[bathy] = bathymetry_extract(bathymetry_dir,bounding_region([],xcoords,ycoords));
figure
bathymetry_section(bathy,xcoords,ycoords,xref)
```

[Back](https://github.com/lnferris/ocean_data_tools#adding-bathymetry-to-existing-plots-1)
### general_region_subset

#### Syntax

```Matlab
subobject] = general_region_subset(object,xcoords,ycoords)
```
#### Description

``[subobject] = general_region_subset(object,xcoords,ycoords)`` subsets object into a polygon region specified by ``xcoords`` and ``ycoords`` (vertices of the polygon); where ``object`` is a struct created by any of the ``_build`` functions in ocean_data_tools (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``). 

``xcoords`` and ``ycoords`` are vectors of coordinates. Rows or columns are fine. -180/180 or 0/360 notation should match that of ``object``

``subobject`` is a struct which is structurally identical to ``object`` but contains only data within the polygon region

#### Example 1

```Matlab

% Get variable information:

argo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/argo/*profiles*.nc';
listing = dir(argo_dir); 
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '01-Nov-2015 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Choose a region for subsetting the uniform struct:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
general_map(argo,bathymetry_dir)
[xcoords,ycoords] = region_select(); % click desired  region on the figure

```
<img src="https://user-images.githubusercontent.com/24570061/88415528-a684a000-cdac-11ea-8f79-189818c2351c.png" width="600">

```Matlab
% Subset the struct and remap:

[subargo] = general_region_subset(argo,xcoords,ycoords); 
general_map(subargo,bathymetry_dir)

```
<img src="https://user-images.githubusercontent.com/24570061/88416080-8e615080-cdad-11ea-8014-e245935da2cf.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#general-functions-for-subsetting-and-plotting-uniform-structs-1)
### mocha_build_profiles

#### Syntax

```Matlab
[mocha] = mocha_build_profiles(month,xcoords,ycoords)
[mocha] = mocha_build_profiles(month,xcoords,ycoords,zgrid)
```
#### Description

``[mocha] = mocha_build_profiles(month,xcoords,ycoords)`` builds a unform struct, ``mocha`` of profiles from the MOCHA Mid-Atlantic Bight climatology, pulling profiles nearest to coordinates specified by ``xcoords`` and ``ycoords``. The calendar month is specified by ``month``.

``[mocha] = mocha_build_profiles(month,xcoords,ycoords,zgrid)`` depth-interpolates the profiles to a vertical grid of ``zgrid``, in meters. ``zgrid=2`` would produce profiles interpolated to 2 meter vertical grid.

``xcoords`` and ``ycoords`` are vectors of coordinates. Rows or columns are fine, and both -180/180 or 0/360 notation are fine.

``month`` is an integer between 1 (January) and 12 (December).

#### Example 1

```Matlab

% Setup nctoolbox:

setup_nctoolbox

% Plot surface temperature:

month = 10; % Month (1 through 12).
depth = 0;
variable = 'temperature'; %  'temperature' 'salinity'
region = [34 42  -80 -70]; % [30 48 -80 -58]
mocha_simple_plot(month,depth,variable,region)

% Click stations on the plot to create a coordinate list:

[xcoords,ycoords] = transect_select('densify',10); % click desired transect on the figure, densify selection by 10x 

```

<img src="https://user-images.githubusercontent.com/24570061/88334226-73d09e00-ccff-11ea-867d-860d64744dc0.png" width="600">

```Matlab

% Build a uniform struct of profiles:

zgrid = 1; % vertical grid for linear interpolation in meters
[mocha] = mocha_build_profiles(month,xcoords,ycoords,zgrid); % zgrid optional, no interpolation if unspecified

% Make plots:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
general_map(mocha,bathymetry_dir,'2Dcontour')
general_section(mocha,'temperature','stn','depth',1,1)
```
<img src="https://user-images.githubusercontent.com/24570061/88334243-78955200-ccff-11ea-8196-4db6298cdec5.png" width="600">
<img src="https://user-images.githubusercontent.com/24570061/88334248-79c67f00-ccff-11ea-926b-a713efbb94d0.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#building-uniform-structs-from-data-sources-1)

### wod_build

#### Syntax

```Matlab
[wod] = wod_build(wod_dir,variable_list)
```
#### Description

``wod_build(wod_dir,variable_list)`` loads profiles in path ``wod_dir`` into the struct ``wod`` with all variables specified in ``variable_list``. Variables lon, lat, date, z are included automatically.

``wod_dir`` is a character array search path with wildcards. The search path should be the path to the netcdf files themselves, not their directory. 

``wod`` is a uniform struct containing data from profiles in the path. Some data is included automatically while some must be specified. The variables lon, lat, date, and z are included automatically. Additional variables must be specified in ``variable_list``, a cell array where each element is the string name of a variable.

#### Example 1


% wod_build

```Matlab

% Get variable information:

wod_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/wod/*.nc'; % included
listing = dir(wod_dir));
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load data in path:

variable_list = {'Temperature','Salinity'}; % Variables to read (besides lon, lat, date, z).
[wod] = wod_build(wod_dir,variable_list);

% Plot profiles:

general_profiles(wod,'Temperature','depth')

```
<img src="https://user-images.githubusercontent.com/24570061/88361566-748d2280-cd47-11ea-82a7-0458d6e2c8dc.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#building-uniform-structs-from-data-sources-1)

### argo_platform_map

#### Syntax

```Matlab
argo_platform_map(argo)
argo_platform_map(argo,annotate)
```
#### Description

``argo_platform_map(argo)`` plots locations of Argo profiles in ``argo``, coloring markers by the specific Argo platform which made the measurements; where ``argo`` is a struct created by ``argo_build``.

``argo_platform_map(argo,annotate)`` adds number annotations to the markers which correspond to a longer Argo platform ID in the legend. By default ``annotate=0``. Set ``annotate=1`` to turn on annotation.


#### Example 1

```Matlab

% Get variable information:

argo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/argo/*profiles*.nc';
listing = dir(argo_dir); 
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '01-Nov-2015 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Plot locations of Argo profiles, coloring by unique platform:

annotate = 1; 
argo_platform_map(argo,annotate) % annotate optional,  1=on 0=off
bathymetry_plot(bathymetry_extract(bathymetry_dir,bounding_region(argo)),'2Dcontour')

```
<img src="https://user-images.githubusercontent.com/24570061/88316847-6955da80-cce6-11ea-8bb0-d9d0523a3a29.png" width="700">

[Back](https://github.com/lnferris/ocean_data_tools#additional-functions-for-inspecting-argo-data-1)

### bathymetry_plot

#### Syntax

```Matlab
bathymetry_plot(bathy,ptype)
```
#### Description

``bathymetry_plot(bathy,ptype)`` makes a 2D (latitude vs. longitude) or 3D (latitude vs. longitude vs. depth) plot from ``bathy``, where ``bathy`` is a struct of Smith & Sandwell Global Topography created using ``bathymetry_extract``. ``type = '2Dscatter'`` or ``'2Dcontour'`` or ``'3Dsurf'`` specifies the plot type.
                     
#### Example 1

```Matlab

% Setup nctoolbox: 

setup_nctoolbox

% Plot a 3-D velocity domain from Operational Mercator:

model_type = 'mercator'; % 'hycom' 'mercator'
source = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/mercator/global-analysis-forecast-phy-001-024_1593408360353.nc'; 
date = '18-Mar-2020 00:00:00';   
variable = 'thetao'; 
region = [60.0, 70.0 ,-80, -60];      % [-90 90 -180 180]
variable = 'velocity'; 
model_domain_plot(model_type,source,date,variable,region)

% Add bathymetry:

[bathy] = bathymetry_extract(bathymetry_dir,region);
bathymetry_plot(bathy,'3Dsurf')
caxis([0 1])

```
<img src="https://user-images.githubusercontent.com/24570061/88409944-ab912180-cda3-11ea-84bc-f848a4f795bc.png" width="700">

[Back](https://github.com/lnferris/ocean_data_tools#adding-bathymetry-to-existing-plots-1)

### mocha_simple_plot

#### Syntax

```Matlab
mocha_simple_plot(month,depth,variable,region)
```
#### Description

``mocha_simple_plot(month,depth,variable,region)`` plots the nearest available depth-level to ``depth``. ``variable`` specifies the parameter to be plotted and ``region`` is the rectangular region to be plotted. The calendar month is specified by ``month``.

``month`` is an integer between 1 (January) and 12 (December).

``depth`` is (a single, double, integer) indicates negative meters below the surface.

``variable`` is a string or character array and is the name of the parameter to be plotted 

``region`` is a vector containing the bounds [S N W E] of the region to be plotted, -180/180 or 0/360 longtitude format is fine.  Limits may cross the dateline e.g. [35 45 170 -130] but this is a Mid-Atlantic product so there is no reason to try that.

#### Example 1

```Matlab

% Setup nctoolbox:

setup_nctoolbox

% Plot surface temperature:

month = 10; % Month (1 through 12).
depth = 0;
variable = 'temperature'; %  'temperature' 'salinity'
region = [34 42  -80 -70]; % [30 48 -80 -58]
mocha_simple_plot(month,depth,variable,region)
```

<img src="https://user-images.githubusercontent.com/24570061/88336317-c52e5c80-cd02-11ea-8b69-1e671bf6536b.png" width="600">


[Back](https://github.com/lnferris/ocean_data_tools#plotting-gridded-data-without-building-structs-1)

### bounding_region

#### Syntax

```Matlab
[region] = bounding_region(object)
[region] = bounding_region(object,xcoords,ycoords)
[region] = bounding_region([],xcoords,ycoords)
```
#### Description

``[region] = bounding_region(object)`` finds a rectangular ``region`` = [S N W E]  around a struct ``object``; where ``object`` is a struct created by any of the ``_build`` functions in ocean_data_tools (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``). 

``[region] = bounding_region(object,xcoords,ycoords)`` ensures that the ``region`` bounding the above struct also ecompasses the points specified  by ``xcoords`` (longitude) and ``ycoords`` (latitude). This is useful for extracting bathymetry only once before using ``bathymetry_plot`` and a ``bathymetry_section``.

``[region] = bounding_region([],xcoords,ycoords)`` finds a rectangular ``region``  = [S N W E] around the points specified  by ``xcoords`` (longitude) and ``ycoords`` (latitude).

``xcoords`` and ``ycoords`` are vectors of coordinates. Rows or columns are fine, and both -180/180 or 0/360 notation are fine when using this function with ``bathymetry_extract``.

#### Example 1

```Matlab

% Map Argo profiles, coloring by platform:

argo_platform_map(argo,1)

% Extract relevant bathymetry around struct argo:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
region = bounding_region(argo);
[bathy] = bathymetry_extract(bathymetry_dir,region);

% Add bathymetry contours:

bathymetry_plot(bathy,region,'2Dcontour')
```
<img src="https://user-images.githubusercontent.com/24570061/88435475-430c6980-cdd0-11ea-9fa8-417bf9b71583.png" width="600">

#### Example 2

```Matlab
% Find the region around a list of coordinates (to be later used with bathymetry_section):

xcoords = -75:1/48:-74;
ycoords = 65:1/48:66;
region = bounding_region([],xcoords,ycoords);
```

[Back](https://github.com/lnferris/ocean_data_tools#adding-bathymetry-to-existing-plots-1)
### mocha_domain_plot

#### Syntax

```Matlab
mocha_domain_plot(month,variable,region)
```
#### Description

``mocha_domain_plot(month,variable,region)`` plots all depth-levels of ``variable`` over the specified ``region``. The calendar month is specified by ``month``.

``month`` is an integer between 1 (January) and 12 (December).

``variable`` is a string or character array and is the name of the parameter to be plotted.

``region`` is a vector containing the bounds [S N W E] of the region to be  plotted, -180/180 or 0/360 longtitude format is fine.  Limits may cross the dateline e.g. [35 45 170 -130] but this is a Mid-Atlantic product so there is no reason to try that.

#### Example 1


```Matlab

% Setup nctoolbox:

setup_nctoolbox

% Plot a 3-D temperature domain:

month = 10; % Month (1 through 12).
variable = 'temperature'; %  'temperature' 'salinity'
region = [34 42  -80 -70]; % [30 48 -80 -58]
mocha_domain_plot(month,variable,region)

% Add bathymetry:
bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
bathymetry_plot(bathymetry_dir,region,'3Dsurf')
caxis([27 38])

```

<img src="https://user-images.githubusercontent.com/24570061/88339261-ba29fb00-cd07-11ea-97a3-bef5868a1fa4.png" width="600">


[Back](https://github.com/lnferris/ocean_data_tools#plotting-gridded-data-without-building-structs-1)

### transect_select

#### Syntax

```Matlab
[xcoords,ycoords] = transect_select() 
[xcoords,ycoords] = transect_select(scheme,value)
```
#### Description

``[xcoords,ycoords] = transect_select()`` creates a list of x and y coordinates selected by clicking stations on an existing (latitude vs. longitude) plot, returning them as ``xcoords`` and ``ycoords``

``[xcoords,ycoords] = transect_select(scheme,value)`` auto-generates additional stations with based on the scheme chosen. ``scheme = 'densify'`` auto-generates additional stations with the multiplier ``value``; ``value=10`` would fill in 10 stations for every station clicked using linear interpolation of complex coordinates. ``scheme = 'spacing'`` auto-generates additional stations with the specified spacing ``value``, where ``value`` is the longitude or latitude spacing in degrees; ``value=0.5`` would fill in stations such that stations are 0.5 degrees apart. 

If ``scheme = 'densify'``, ``value`` (no units) should be an integer. If it is not an integer it will be rounded to an integer.

If ``scheme = 'spacing'``, ``value`` (in degrees) should be single, double, or integer and represents the prescribed grid spacing between auto-generated stations in degrees. The spacing criterion is applied to each dimension separately and is not the diagonal displacement.

``xcoords`` and ``ycoords`` are vectors of coordinates representing a polygonal chain. -180/180 or 0/360 notation will match that of the existing plot.

#### Example 1


```Matlab

% Plot HYCOM surface salinity:

model_type = 'hycom'; 
source = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7';
date = '28-Aug-2017 00:00:00';  
variable = 'salinity';                
region = [-5.0, 45.0 ,160,-150 ];      
depth = -150;                                                   
model_simple_plot(model_type,source,date,variable,region,depth);

% Click stations on the plot to create a coordinate list:

[xcoords,ycoords] = transect_select('densify',10); % click desired transect on the figure, densify selection by 10x

```
<img src="https://user-images.githubusercontent.com/24570061/88406388-9f569580-cd9e-11ea-9871-e4d55941d7c4.png" width="500">

```Matlab

% Build a uniform struct from the coordinates:

variable_list = {'water_temp','salinity'}; % 'water_u' 'water_v' 'water_temp' 'salinity'
[hycom] =  model_build_profiles(source,date,variable_list,xcoords,ycoords);

% Map stations:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
general_map(hycom,bathymetry_dir,'2Dcontour')

```
<img src="https://user-images.githubusercontent.com/24570061/88406404-a67da380-cd9e-11ea-8d49-bd4db591c282.png" width="500">


[Back](https://github.com/lnferris/ocean_data_tools#miscellaneous-utilities-1)

### region_select

#### Syntax

```Matlab
[xcoords,ycoords] = region_select()
```
#### Description

``[xcoords,ycoords] = region_select()`` creates a list of x and y coordinates (which represent vertices of a polygon) selected by clicking stations on an existing (latitude vs. longitude) plot, returning them as ``xcoords`` and ``ycoords``.

``xcoords`` and ``ycoords`` are vectors of coordinates representing vertices of a polygon. -180/180 or 0/360 notation will match that of the existing plot.

#### Example 1

```Matlab

% Get variable information:

argo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/argo/*profiles*.nc';
listing = dir(argo_dir); 
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '01-Nov-2015 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Choose a region for subsetting the uniform struct:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
general_map(argo,bathymetry_dir)
[xcoords,ycoords] = region_select(); % click desired  region on the figure

```
<img src="https://user-images.githubusercontent.com/24570061/88415528-a684a000-cdac-11ea-8f79-189818c2351c.png" width="600">

```Matlab
% Subset the struct:

[subargo] = general_region_subset(argo,xcoords,ycoords); 

```


[Back](https://github.com/lnferris/ocean_data_tools#miscellaneous-utilities-1)

### general_depth_subset

#### Syntax

```Matlab
[subobject] = general_depth_subset(object,zrange)
[subobject] = general_depth_subset(object,zrange,depth_list)
```
#### Description

``[subobject] =  general_depth_subset(object,zrange)`` subsets ``object`` by depth-range ``zrange``; where ``object`` is a struct created by any of the ``_build`` functions in ocean_data_tools (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``). The default depth-variable used to subset is ``'depth'``. ``zrange`` is a 2-element vector e.g. ``zrange=[0 200]`` in meters or dbar. Order does not matter, but the sign convention should be the same as the depth variable in ``object``.

``[subobject] =  general_depth_subset(object,zrange,depth_list)`` enables the user to specify one or more depth variables (instead of using default ``'depth'``) e.g. ``depth_list = {'pressure'}`` or ``depth_list = {'pressure','z','depth','depth_vke'}``.

``subobject`` is a struct which is structurally identical to ``object`` but contains only data within the specified depth range. In other words, profiles within ``object`` have been truncated.

#### Example 1


```Matlab
% Build a uniform struct from HYCOM and plot a temperature section:

[hycom] =  model_build_profiles(source,date,variable_list,xcoords,ycoords,zgrid);
general_section(hycom,'water_temp','lat','depth',1,1)

```
<img src="https://user-images.githubusercontent.com/24570061/88417509-0892d480-cdb0-11ea-9685-0da4b82f99f4.png" width="600">

```Matlab
% Subset to upper 450 meters and replot the temperature section:

[hycom] =  general_depth_subset(hycom,[-450 0]);
general_section(hycom,'water_temp','lat','depth',1,1)

```
<img src="https://user-images.githubusercontent.com/24570061/88417524-10527900-cdb0-11ea-8a5c-b8c6607a2386.png" width="600">


[Back](https://github.com/lnferris/ocean_data_tools#general-functions-for-subsetting-and-plotting-uniform-structs-1)

### woa_domain_plot

#### Syntax

```Matlab
woa_domain_plot(variable,time,region)
```
#### Description

``woa_domain_plot(variable,time,region)`` plots all depth levels of  World Ocean Atlas 2018 Statistical Mean for All Decades, Objectively Analyzed Mean Fields at Standard Depth Levels over the specified ``region``. ``variable`` specifies the parameter to be plotted and ``region`` is the rectangular region to be plotted. ``time`` specifies monthly or annual climatology; ``time = '00'`` for annual climatology and ``'01'`` ``'10'`` etc. for monthly climatology. The function builds the url, extracting the maximum resolution available (typically 0.25-deg or 1.00-degree grid). 

Available variables are:

``'temperature'`` (degrees Celsius)    
``'salinity'`` (psu)                    
``'oxygen'`` (umol/kg)                 
``'o2sat'`` (%)

``'AOU'`` (umol/kg)                
``'silicate'`` (umol/kg)          
``'phosphate'`` (umol/kg)   
``'nitrate'`` (umol/kg)     

``time`` is a string or character array. ``'00'`` is annual climatology, while other codes e.g. ``'02'`` (February) or ``'11'`` (November) indicate monthly climatology.

``variable`` is a string or character array and is the name of the parameter to be plotted.

``region`` is a vector containing the bounds [S N W E] of the region to be plotted, -180/180 or 0/360 longtitude format is fine.  Limits may cross the dateline e.g. [35 45 170 -130].

#### Example 1


```Matlab

% Setup nctoolbox:

setup_nctoolbox

% Plot a 3-D nitrate domain:

variable = 'nitrate';
time = '03';
region = [-5.0, 45.0 ,-120, -150]; 
woa_domain_plot(variable,time,region)

```
<img src="https://user-images.githubusercontent.com/24570061/88359635-5e7c6380-cd41-11ea-8f39-62beb8ecf912.png" width="900">

[Back](https://github.com/lnferris/ocean_data_tools#plotting-gridded-data-without-building-structs-1)
### glider_build

#### Syntax

```Matlab
[glider] = glider_build(glider_dir)
[glider] = glider_build(glider_dir,variable_list)
```
#### Description

``[glider] = glider_build(glider_dir,variable_list)`` loads data (a single netCDF file downloaded from [gliders.ioos.us/erddap](https://gliders.ioos.us/erddap/index.html)) from ``glider_dir`` into struct array ``glider``. Glider profiles are loaded with all variables specified in ``variable_list``. 

The only required argument is ``glider_dir``. By default, the following variables are loaded:
```
'profile_id'
'time'
'latitude'
'longitude'
'precise_time'
'depth'
'pressure'
'temperature'
'conductivity'
'salinity'
'density'
'precise_lat'
'precise_lon'
'time_uv'
'lat_uv',
'lon_uv'
'u'
'v' 
```

This default list can be overridden by passing a user-defined ``variable_list``, a cell array where each element is the string name of a variable. 

``glider_dir`` is a character array search path to a single netcdf file downloaded from gliders.ioos.us/erddap). The path should be the path to the netcdf file itself, not its directory. 

#### Example 1


```Matlab

% Get variable information:

glider_dir = '/Users/lnferris/Desktop/ce_311-20170725T1930.nc';  % included
ncdisp(glider_dir) % Peek at netCDF header info to inform choice of variable_list.

% Load glider data:

[glider] = glider_build(glider_dir); 

% Make plots:

figure
general_map(glider,bathymetry_dir)

figure
general_section(glider,'salinity','km','pressure')

```
<img src="https://user-images.githubusercontent.com/24570061/94057510-b66d3000-fdad-11ea-8261-c2dddbf72439.png" width="600">
<img src="https://user-images.githubusercontent.com/24570061/94057516-b79e5d00-fdad-11ea-85e0-3e274feff845.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#building-uniform-structs-from-data-sources-1)

### model_domain_plot

#### Syntax

```Matlab
model_domain_plot(model_type,source,date,variable,region)
```
#### Description

``model_domain_plot(model,source,date,variable,region)`` plots all depth levels of HYCOM or Operational Mercator GLOBAL_ANALYSIS_FORECAST_PHY_001_024 over a particular rectangular ``region``. ``variable`` specifies the parameter to be plotted.  In addition to innate variables, each model has the additional derived variable ``variable='velocity'``. ``model_type ='hycom'`` or ``model_type ='mercator'`` specifies the model used. ``source`` is the url or local path of the relevant dataset.

``source`` (a character array) is the path to either a local netcdf file or an OpenDAP url.

``date`` is a date string in format 'dd-mmm-yyyy HH:MM:SS'. 

``variable`` is a string or character array and is the name of the parameter to be plotted.

``region`` is a vector containing the bounds [S N W E] of the region to be plotted, -180/180 or 0/360 longtitude format is fine.  Limits may cross the dateline e.g. [35 45 170 -130].

HYCOM variables: 
```Matlab
'water_u' 
'water_v' 
'water_temp' 
'salinity' 
'velocity' 
```
Mercator variables: 
```Matlab
'thetao' 
'so' 
'uo' 
'vo' 
'velocity'
```
                     
#### Example 1


```Matlab

% Setup nctoolbox: 

setup_nctoolbox

% Plot a 3-D velocity domain from Operational Mercator:

model_type = 'mercator'; % 'hycom' 'mercator'
source = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/mercator/global-analysis-forecast-phy-001-024_1593408360353.nc'; 
date = '18-Mar-2020 00:00:00';   
variable = 'thetao'; 
region = [60.0, 70.0 ,-80, -60];      % [-90 90 -180 180]
variable = 'velocity'; 
model_domain_plot(model,source,date,variable,region)

% Add bathymetry:

bathymetry_plot(bathymetry_dir,region,'3Dsurf')
caxis([0 1])

```
<img src="https://user-images.githubusercontent.com/24570061/88409944-ab912180-cda3-11ea-84bc-f848a4f795bc.png" width="700">

[Back](https://github.com/lnferris/ocean_data_tools#plotting-gridded-data-without-building-structs-1)
### bathymetry_extract

#### Syntax

```Matlab
[bathy] =  bathymetry_extract(bathymetry_dir,region)
```
#### Description

``[bathy] =  bathymetry_extract(bathymetry_dir,region)`` extracts Smith & Sandwell Global Topography in path ``bathymetry_dir`` over the specified rectangular ``region``. Struct ``bathy`` has fields ``z``, ``lat``, and ``lon``. Struct bathy has fields z, lat, and lon and contains only the data from the specified region.

``bathymetry_dir`` is the character array path to the Smith & Sandwell Global Topography file "topo_20.1.nc"

``region`` is a vector containing the bounds [S N W E] of the search region, with limits [-90 90 -180 180]. Limits may cross the dateline e.g. [35 45 170 -130].

#### Example 1

```Matlab
% Extract relevant bathymetry over a region:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
region = [-60.0 -50.0 150.0 160.0];
[bathy] = bathymetry_extract(bathymetry_dir,region);
```

#### Example 2

```Matlab
% Extract relevant bathymetry around struct argo:

bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
region = bounding_region(argo);
[bathy] = bathymetry_extract(bathymetry_dir,region);

```

[Back](https://github.com/lnferris/ocean_data_tools#adding-bathymetry-to-existing-plots-1)

### argo_build

#### Syntax

```Matlab
[argo,matching_files] = argo_build(argo_dir)
[argo,matching_files] = argo_build(argo_dir,region)
[argo,matching_files] = argo_build(argo_dir,region,start_date)
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date)
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list)
```
#### Description

``[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list)`` searches pathway ``argo_dir`` for profiles meeting the search criteria ``region``, ``start_date``, and ``end_date``. Profiles are loaded into the struct array ``argo`` with all variables specified in ``variable_list``. Files containing matching profiles are listed in ``matching_files``.

The only required argument is ``argo_dir``. The default state is to load all profiles in path ``argo_dir``, writing variables TEMP_ADJUSTED and PSAL_ADJUSTED into the uniform struct ``argo``.

``argo_dir`` is a character array search path with wildcards. The search path should be the path to the netcdf files themselves, not their directory. 

``region`` is a vector containing the bounds [S N W E] of the search region, with limits [-90 90 -180 180]. Limits may cross the dateline e.g. [35 45 170 -130].

``start_date`` and ``end_date`` are date strings in format ``'dd-mmm-yyyy HH:MM:SS'``.

``argo`` is a uniform struct containing data from the profiles matching the region and date criteria. Some data is included automatically while some must be specificed. The variables PLATFORM_NUMBER, LONGITUDE, LATITUDE, JULD, and PRES_ADJUSTED are included automatically. Additional variables must be specified in ``variable_list``, a cell array where each element is the string name of a variable.

``matching_files`` is a string array where each string is the full path to a file which contained a profile matching the region and date criteria.

#### Example 1


```Matlab

% Get variable information:

argo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/argo/*profiles*.nc';
listing = dir(argo_dir); 
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '01-Nov-2015 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Make plots:

general_profiles(argo,'TEMP_ADJUSTED','depth')
general_map(argo,bathymetry_dir,'2Dcontour')

```
<img src="https://user-images.githubusercontent.com/24570061/88301724-fd1dab80-ccd2-11ea-9ea7-7badf1424865.png" width="600">
<img src="https://user-images.githubusercontent.com/24570061/88301788-11fa3f00-ccd3-11ea-9cdf-1622f701bfe9.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#building-uniform-structs-from-data-sources-1)

### GSW Example

```Matlab

% Using ocean_data_tools, build a uniform struct from HYCOM and subset to upper 450 meters:

[hycom] =  model_build_profiles(source,date,variable_list,xcoords,ycoords,zgrid);
[hycom] =  general_depth_subset(hycom,[0 450]);

% Using ocean_data_tools, plot temperature and salinity sections:

general_section(hycom,'water_temp','lat','depth',1,1)
general_section(hycom,'salinity','lat','depth',1,1)

```
<img src="https://user-images.githubusercontent.com/24570061/88403123-16d5f600-cd9a-11ea-9ca6-c1e0403a44da.png" width="700">
<img src="https://user-images.githubusercontent.com/24570061/88403126-176e8c80-cd9a-11ea-846f-8e97e80f3805.png" width="700">

```Matlab

% Using GSW, append the struct with absolute salinity, conservative temperature, and density:

[hycom.SA, ~] = gsw_SA_from_SP(hycom.salinity,-hycom.depth,hycom.lon,hycom.lat);
hycom.CT = gsw_CT_from_t(hycom.SA,hycom.water_temp,-hycom.depth);
hycom.rho = gsw_rho(hycom.SA,hycom.CT,-hycom.depth);

% Using ocean_data_tools, plot a density section:

general_section(hycom,'rho','lat','depth',1,1)

```
<img src="https://user-images.githubusercontent.com/24570061/88403129-19d0e680-cd9a-11ea-9b4e-4e733cdb3c5b.png" width="700">

[Back](https://github.com/lnferris/ocean_data_tools#dependencies-1)

### argo_profiles_map

#### Syntax

```Matlab
argo_profiles_map(argo)
argo_profiles_map(argo,annotate)
```
#### Description

``argo_profiles_map(argo)`` plots locations of Argo profiles in struct ``argo``, coloring markers by profile; where ``argo`` is a struct created by ``argo_build``. The colors of profiles corresponds to those of ``argo_profiles`` called on the same struct.

``argo_profiles_map(argo,annotate)`` adds number annotations to the markers. By default ``annotate=0``. Set ``annotate=1`` to turn on annotation. The annotations of profiles correspond to those of ``argo_profiles`` called on the same struct.

#### Example 1


```Matlab

% Get variable information:

argo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/argo/*profiles*.nc';
listing = dir(argo_dir); 
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '28-Dec-2016 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Plot profiles with annotations:

variable = 'TEMP_ADJUSTED'; % See object for options.
annotate = 1; 
argo_profiles(argo,variable,annotate) % annotate optional,  1=on 0=off

```
<img src="https://user-images.githubusercontent.com/24570061/88327048-b04acc80-ccf4-11ea-8e1d-4ae00634953c.png" width="600">

```Matlab
% Map profiles with annotations:

annotate = 1; 
argo_profiles_map(argo,annotate) % annotate optional,  1=on 0=off
bathymetry_plot(bathymetry_extract(bathymetry_dir,bounding_region(argo)),'2Dcontour') % add bathymetry contours
```

<img src="https://user-images.githubusercontent.com/24570061/88327053-b2149000-ccf4-11ea-8d1d-c0493bf9ef43.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#additional-functions-for-inspecting-argo-data-1)
### argo_profiles

#### Syntax

```Matlab
argo_profiles(argo,variable) 
argo_profiles(argo,variable,annotate)
```
#### Description

``argo_profiles(argo,variable)`` plots vertical profiles of the specified variable in struct ``argo`` as a function of depth (PRES_ADJUSTED); where ``argo`` is a struct created by ``argo_build`` and ``variable`` is a field name.
 
``argo_profiles(argo,variable,annotate)`` adds number annotations to the markers by default ``annotate=0``. Set ``annotate=1`` to turn on annotation. The annotations of profiles correspond to those of ``argo_profiles_map`` called on the same struct.

``variable`` is the string name of the field (of ``argo``) to be plotted as a vertical profile

#### Example 1


```Matlab

% Get variable information:

argo_dir = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/argo/*profiles*.nc';
listing = dir(argo_dir); 
ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.

% Load Argo data from west of New Zealand:

region = [-60.0 -50.0 150.0 160.0]; %  Search region [-90 90 -180 180]
start_date = '28-Dec-2016 00:00:00';
end_date = '01-Jan-2017 00:00:00';
variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
[argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);

% Plot profiles with annotations:

variable = 'TEMP_ADJUSTED'; % See object for options.
annotate = 1; 
argo_profiles(argo,variable,annotate) % annotate optional,  1=on 0=off

```
<img src="https://user-images.githubusercontent.com/24570061/88327048-b04acc80-ccf4-11ea-8e1d-4ae00634953c.png" width="600">

```Matlab
% Map profiles with annotations:

annotate = 1; 
argo_profiles_map(argo,annotate) % annotate optional,  1=on 0=off
bathymetry_plot(bathymetry_extract(bathymetry_dir,bounding_region(argo)),'2Dcontour') % add bathymetry contours
```

<img src="https://user-images.githubusercontent.com/24570061/88327053-b2149000-ccf4-11ea-8d1d-c0493bf9ef43.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#additional-functions-for-inspecting-argo-data-1)

### model_simple_plot

#### Syntax

```Matlab
[data,lat,lon] = model_simple_plot(model_type,source,date,variable,region,depth)
[data,lat,lon] = model_simple_plot(model_type,source,date,variable,region,depth,arrows)
```
#### Description

``[data,lat,lon] = model_simple_plot(model_type,source,date,variable,region,depth)`` plots the nearest available depth-level of HYCOM or Operational Mercator GLOBAL_ANALYSIS_FORECAST_PHY_001_024. ``variable`` specifies the parameter to be plotted and ``region`` is the rectangular region to be plotted. ``model_type ='hycom'`` or ``model_type ='mercator'`` specifies the model used. ``source`` is the url or local path of the relevant dataset. ``data``, ``lat``, and ``lon`` from the plotted layer are available outputs.

``[data,lat,lon] = model_simple_plot(model_type,source,date,variable,region,depth,arrows)`` adds directional arrows if it is a velocity magnitude plot. ``arrows=1`` is on, ``arrows=0`` is off.

``source`` (a character array) is the path to either a local netcdf file or an OpenDAP url.

``date`` is a date string in format 'dd-mmm-yyyy HH:MM:SS'. 

``variable`` is a string or character array and is the name of the parameter to be plotted.

``depth`` is (a single, double, integer) indicates negative meters below the surface.

``region`` is a vector containing the bounds [S N W E] of the region to be plotted, -180/180 or 0/360 longtitude format is fine.  Limits may cross the dateline e.g. [35 45 170 -130].

``data``, ``lon``, and ``lat`` are double arrays containing the plotted data layer. As such, this function can be used to extract data layers from HYCOM or Operational Mercator GLOBAL_ANALYSIS_FORECAST_PHY_001_024.

HYCOM variables: 
```Matlab
'water_u' 
'water_v' 
'water_temp' 
'salinity' 
'velocity' 
'surf_el' 
'water_u_bottom' 
'water_v_bottom' 
'water_temp_bottom' 
'salinity_bottom'
```
Mercator variables: 
```Matlab
'thetao' 
'so' 
'uo' 
'vo' 
'velocity'
'mlotst' 
'siconc'
'usi' 
'vsi' 
'sithick'
'bottomT' 
'zos'
```
                     

#### Example 1


```Matlab

% Setup nctoolbox: 

setup_nctoolbox

% Plot temperature at the depth level closest to 150m:

model_type = 'mercator'; % 'hycom' 'mercator'
source = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/mercator/global-analysis-forecast-phy-001-024_1593408360353.nc'; 
date = '18-Mar-2020 00:00:00';   
variable = 'thetao'; 
region = [60.0, 70.0 ,-80, -60];      % [-90 90 -180 180]
depth = -150;                
arrows = 0;  
model_simple_plot(model_type,source,date,variable,region,depth,arrows)

```
<img src="https://user-images.githubusercontent.com/24570061/88408553-b8147a80-cda1-11ea-98bf-e1d9b45aa53c.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#plotting-gridded-data-without-building-structs-1)

### general_section

#### Syntax

```Matlab
general_section(object,variable,xref,zref)
general_section(object,variable,xref,zref,interpolate)
general_section(object,variable,xref,zref,interpolate,contours)
```
#### Description

``general_section(object,variable,xref,zref)`` creates a section plot from ``object``; where ``object`` is a struct created by any of the ``_build`` functions in ocean_data_tools (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``). The color field is specified by ``variable``. ``xref`` and ``zref`` specify fields to use for the x-axis and z-axis. Alteratively pass ``xref = 'km'`` to plot in along-track distance, assuming spherical earth.

``general_section(object,variable,xref,zref,interpolate)`` interpolates the plot using the MATLAB ``shading`` function. ``interpolate=1`` for on, ``interpolate=0`` for off.

``general_section(object,variable,xref,zref,interpolate,contours)`` adds contours to the section plot. ``contours=1`` for on, ``contours=0`` for off.

``variable`` is the string name of the field (of ``object``) to be plotted as the color variable of the section plot

``zref`` is the string name of the field (of ``object``) to be plotted as the depth variable of the section plot

``xref`` is the string name of the field (of ``object``) to be plotted as the horizontal variable of the section plot, usually ``'stn'``, ``'lat'``, or ``'lon'``. Alteratively pass ``xref = 'km'`` to plot in along-track distance, assuming spherical earth.


#### Example 1


```Matlab

% Setup nctoolbox:

setup_nctoolbox

% Built a uniform struct from MOCHA climatology:

month = 10; % Month (1 through 12).
depth = 0;
variable = 'temperature'; %  'temperature' 'salinity'
region = [34 42  -80 -70]; % [30 48 -80 -58]
mocha_simple_plot(month,depth,variable,region)
[xcoords,ycoords] = transect_select('densify',10); % click desired transect on the figure, densify selection by 10x 
zgrid = 1; % vertical grid for linear interpolation in meters
[mocha] = mocha_build_profiles(month,xcoords,ycoords,zgrid); % zgrid optional, no interpolation if unspecified
```
<img src="https://user-images.githubusercontent.com/24570061/88334226-73d09e00-ccff-11ea-867d-860d64744dc0.png" width="600">

```Matlab

% Make a temperature section:

object = mocha; % argo, cruise, hycom, mercator, woa, wod
variable = 'temperature'; % see particular struct for options
xref = 'stn';  % 'lat' 'lon' 'stn';
zref = 'depth;  See particular struct for options
interpolate = 1; % 1=on 0=off
contours = 1; % 1=on 0=off
general_section(object,variable,xref,zref,interpolate,contours)
```
<img src="https://user-images.githubusercontent.com/24570061/88334248-79c67f00-ccff-11ea-926b-a713efbb94d0.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#general-functions-for-subsetting-and-plotting-uniform-structs-1)


### woa_build_profiles

#### Syntax

```Matlab
[woa] =   woa_build_profiles(variable_list,time,xcoords,ycoords)
[woa] =   woa_build_profiles(variable_list,time,xcoords,ycoords,zgrid)
```
#### Description

``[woa] =   woa_build_profiles(variable_list,time,xcoords,ycoords)`` builds a struct of profiles from World Ocean Atlas 2018 Statistical Mean for All Decades Objectively Analyzed Mean Fields at Standard Depth Levels, pulling profiles nearest to coordinates specified by ``xcoords`` and ``ycoords``. ``time`` specifies monthly or annual climatology; ``time ='00'`` for annual climatology and ``'01'`` ``'10'`` etc. for monthly climatology. Profiles are loaded into the struct array ``woa`` with all variables specified in ``variable_list``. The function builds the url, extracting the maximum resolution available (typically 0.25-deg or 1.00-degree grid). Resolution depends on the variable. If the user requests only 0.25-degree variables in variable_list, data will be returned in  0.25-degree resolution. If any requested variable is coarser (1-degree) all variables will be returned in 1-degree resolution.

``[woa] =  woa_build_profiles(variable_list,time,xcoords,ycoords,zgrid)`` depth-interpolates the profiles to a vertical grid of ``zgrid``, in meters. ``zgrid=2`` would produce profiles interpolated to 2 meter vertical grid.

``time`` is a string or character array. ``'00'`` is annual climatology, while other codes e.g. ``'02'`` (February) or ``'11'`` (November) indicate monthly climatology.

``variable_list`` is a cell array where each element is the string name of a variable to be read and included in struct ``woa``.

``xcoords`` and ``ycoords`` are vectors of coordinates. Rows or columns are fine, and both -180/180 or 0/360 notation are fine.

Available variables are:

``'temperature'`` (degrees Celsius)    
``'salinity'`` (psu)                    
``'oxygen'`` (umol/kg)                 
``'o2sat'`` (%)

``'AOU'`` (umol/kg)                
``'silicate'`` (umol/kg)          
``'phosphate'`` (umol/kg)   
``'nitrate'`` (umol/kg)                      

#### Example 1


```Matlab

% Setup nctoolbox:

setup_nctoolbox

% Plot surface nitrate:

variable = 'nitrate';
time = '03';
region = [-5.0, 45.0 ,-120, -150]; 
depth = -0; 
woa_simple_plot(variable,time,region,depth)

% Click stations on the plot to create a coordinate list:

[xcoords,ycoords] = transect_select('densify',10); % click desired transect on the figure, densify selection by 10x 

```
<img src="https://user-images.githubusercontent.com/24570061/88359631-5c1a0980-cd41-11ea-8e2b-c331e28e1e09.png" width="900">

```Matlab

% Build a uniform struct of profiles:

variable_list = {'temperature','salinity','oxygen'}; % 'temperature' 'salinity' 'oxygen' 'o2sat' 'AOU' 'silicate' 'phosphate' 'nitrate'
time = '00'; % '00' for annual climatology '01' '10' etc. for monthly climatology
zgrid = 1; % vertical grid for linear interpolation in meters
[woa] =  woa_build_profiles(variable_list,time,xcoords,ycoords,zgrid); % zgrid optional, no interpolation if unspecified
[woa] = general_remove_duplicates(woa); % thin struct to gridding of source (optional)

% Make a section plot:

general_section(woa,'salinity','lon','depth',1,1)
```

<img src="https://user-images.githubusercontent.com/24570061/88359633-5d4b3680-cd41-11ea-84bd-ef62b5f42fbb.png" width="900">

[Back](https://github.com/lnferris/ocean_data_tools#building-uniform-structs-from-data-sources-1)

### argo_platform_subset

#### Syntax

```Matlab
[subargo] = argo_platform_subset(argo,platform_id) 
```
#### Description

``[subargo] = argo_platform_subset(argo,platform_id)`` subsets ``argo`` by Argo platform ID into struct ``subargo``; where ``argo`` is a struct created by ``argo_build`` and platform_id is the integer ID

``platform_id`` must be an integer corresponding to the ID of an Argo float.

``subargo`` is a struct which is structurally identical to ``argo`` but contains data only from the chosen Argo float

#### Example 1

```Matlab
% Plot locations of Argo profiles:
argo_platform_map(argo,1)
bathymetry_plot(bathymetry_extract(bathymetry_dir,bounding_region(argo)),'2Dcontour')
```

<img src="https://user-images.githubusercontent.com/24570061/88316847-6955da80-cce6-11ea-8bb0-d9d0523a3a29.png" width="700">

```Matlab
% Subset by platform:

platform_id = 5904421;
[subargo] = argo_platform_subset(argo,platform_id);

% Make a new plot:

argo_platform_map(subargo,1)
bathymetry_plot(bathymetry_extract(bathymetry_dir,bounding_region(argo)),'2Dcontour')

```

<img src="https://user-images.githubusercontent.com/24570061/88324607-ec306280-ccf1-11ea-8f9a-81320046ccf4.png" width="700">

[Back](https://github.com/lnferris/ocean_data_tools#additional-functions-for-inspecting-argo-data-1)

### model_build_profiles

#### Syntax

```Matlab
[model] = model_build_profiles(source,date,variable_list,xcoords,ycoords)
[model] = model_build_profiles(source,date,variable_list,xcoords,ycoords,zgrid)
```
#### Description

``[model] = model_build_profiles(source,date,variable_list,xcoords,ycoords)`` builds a struct of profiles from HYCOM or Operational Mercator GLOBAL_ANALYSIS_FORECAST_PHY_001_024, pulling profiles nearest to coordinates specified by ``xcoords`` and ``ycoords``. Profiles are loaded into the struct array ``model`` with all variables specified in ``variable_list``.

``[model] = model_build_profiles(source,date,variable_list,xcoords,ycoords,zgrid)`` depth-interpolates the profiles to a vertical grid of ``zgrid``, in meters. ``zgrid=2`` would produce profiles interpolated to 2 meter vertical grid.

``source`` (a character array) is the path to either a local netcdf file or an OpenDAP url.

``date`` is a date string in format 'dd-mmm-yyyy HH:MM:SS'. 

``variable_list`` is a cell array where each element is the string name of a variable to be read and included in struct ``model``.

``xcoords`` and ``ycoords`` are vectors of coordinates. Rows or columns are fine, and both -180/180 or 0/360 notation are fine.

HYCOM variables: 
```Matlab
'water_u' 
'water_v' 
'water_temp' 
'salinity' 
```
Mercator variables: 
```Matlab
'thetao' 
'so' 
'uo' 
'vo' 
```
                     
#### Example 1


```Matlab

% Setup nctoolbox: 

setup_nctoolbox

% Plot temperature at the depth level closest to 150m:

model_type = 'mercator'; % 'hycom' 'mercator'
source = '/Users/lnferris/Documents/GitHub/ocean_data_tools/data/mercator/global-analysis-forecast-phy-001-024_1593408360353.nc'; 
date = '18-Mar-2020 00:00:00';   
variable = 'thetao'; 
region = [60.0, 70.0 ,-80, -60];      % [-90 90 -180 180]
depth = -150;                
arrows = 0;  
model_simple_plot(model_type,source,date,variable,region,depth,arrows)

% Click stations on the plot to create a coordinate list:

[xcoords,ycoords] = transect_select('densify',10); % click desired transect on the figure, densify selection by 10x 

```
<img src="https://user-images.githubusercontent.com/24570061/88411026-3c1c3180-cda5-11ea-81d3-d3b656315464.png" width="600">

```Matlab

% Build a uniform struct of profiles:

variable_list = {'thetao','so','uo'}; % thetao' 'so' 'uo' 'vo'
zgrid = 1; % vertical grid for linear interpolation in meters
[mercator] =  model_build_profiles(source,date,variable_list,xcoords,ycoords,zgrid); % zgrid optional, no interpolation if unspecified

% Make plots:

general_map(mercator,bathymetry_dir,'2Dcontour')
general_section(mercator,'thetao','stn','depth',1,1)

```
<img src="https://user-images.githubusercontent.com/24570061/88411140-6b32a300-cda5-11ea-922e-48bf06df90b3.png" width="600">
<img src="https://user-images.githubusercontent.com/24570061/88411172-7685ce80-cda5-11ea-9ae9-0989c763bef9.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#building-uniform-structs-from-data-sources-1)

### woa_simple_plot

#### Syntax

```Matlab
[data,lat,lon] = woa_simple_plot(variable,time,region,depth)
```
#### Description

``[data,lat,lon] = woa_simple_plot(variable,time,region,depth)`` plots the nearest available depth-level to ``depth``. ``variable`` specifies the parameter to be plotted and ``region`` is the rectangular region to be plotted. ``time``specifies monthly or annual climatology; ``time='00'`` for annual climatology and ``'01'`` ``'10'`` etc. for monthly climatology. The function builds the url, extracting the maximum resolution available (typically 0.25-deg or 1.00-degree grid). ``data``, ``lat``, and ``lon`` from the plotted layer are available outputs.

Available variables are:

``'temperature'`` (degrees Celsius)    
``'salinity'`` (psu)                    
``'oxygen'`` (umol/kg)                 
``'o2sat'`` (%)

``'AOU'`` (umol/kg)                
``'silicate'`` (umol/kg)          
``'phosphate'`` (umol/kg)   
``'nitrate'`` (umol/kg)

``time`` is a string or character array. ``'00'`` is annual climatology, while other codes e.g. ``'02'`` (February) or ``'11'`` (November) indicate monthly climatology.

``variable`` is a string or character array and is the name of the parameter to be plotted.

``depth`` is (a single, double, integer) indicates negative meters below the surface.

``region`` is a vector containing the bounds [S N W E] of the region to be plotted, -180/180 or 0/360 longtitude format is fine.  Limits may cross the dateline e.g. [35 45 170 -130].

``data``, ``lon``, and ``lat`` are double arrays containing the plotted data layer. As such, this function can be used to extract data layers from World Ocean Atlas 2018.

#### Example 1


```Matlab

% Setup nctoolbox:

setup_nctoolbox

% Plot surface nitrate:

variable = 'nitrate';
time = '03';
region = [-5.0, 45.0 ,-120, -150]; 
depth = -0; 
woa_simple_plot(variable,time,region,depth)

```
<img src="https://user-images.githubusercontent.com/24570061/88360391-7a810480-cd43-11ea-9472-2d18a306389b.png" width="900">

[Back](https://github.com/lnferris/ocean_data_tools#plotting-gridded-data-without-building-structs-1)

### general_remove_duplicates

#### Syntax

```Matlab
[subobject] = general_remove_duplicates(object)
[subobject] = general_remove_duplicates(object,var3)
```
#### Description

``[subobject] = general_remove_duplicates(object)`` removes profiles without unique coordinate locations; where ``object`` is a struct created by any of the ``_build`` functions in ocean_data_tools (e.g. ``argo``, ``cruise``, ``hycom``, ``mercator``, ``woa``, ``wod``).  This function can be used to remove duplicates when the user accidentally built an ``object`` using a coordinate list (``xcoords``, ``ycoords``) that exceeded the spatial resolution of the model or other raw data source itself.

``[subobject] = general_remove_duplicates(object,var3)`` uses a third field ``var3`` as a uniqueness criterion, usually date. This avoids removing profiles with the same coordinate location but unique dates. ``var3`` should be a fieldname used to filter the data e.g. ``var3='date'``.

``subobject`` is a struct which is structurally identical to ``object`` but contains data only spatially-unique (or also var3-unique) profiles

#### Example 1


```Matlab

% Build a uniform struct from HYCOM, requesting 49 profiles at a higher resolution than the model itself:

source = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7';
date = '28-Aug-2017 00:00:00';
xcoords = -75:1/48:-74;
ycoords = 65:1/48:66;
variable_list = {'water_temp','salinity'}; 
[hycom] = model_build_profiles(source,date,variable_list,xcoords,ycoords);

% Remove duplicate profiles, resulting in 39 unique profiles:

object = hycom;
[subobject] = general_remove_duplicates(object);

```
<img src="https://user-images.githubusercontent.com/24570061/88433788-bd3aef00-cdcc-11ea-99c4-653bca43d1d0.png" width="600">

[Back](https://github.com/lnferris/ocean_data_tools#general-functions-for-subsetting-and-plotting-uniform-structs-1)


