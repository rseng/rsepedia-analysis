# Welcome to leafmap

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)
[![image](https://img.shields.io/pypi/v/leafmap.svg)](https://pypi.python.org/pypi/leafmap)
[![image](https://img.shields.io/conda/vn/conda-forge/leafmap.svg)](https://anaconda.org/conda-forge/leafmap)
[![image](https://pepy.tech/badge/leafmap)](https://pepy.tech/project/leafmap)
[![image](https://github.com/giswqs/leafmap/workflows/docs/badge.svg)](https://leafmap.org)
[![image](https://github.com/giswqs/leafmap/workflows/Linux%20build/badge.svg)](https://github.com/giswqs/leafmap/actions)
[![image](https://img.shields.io/lgtm/grade/python/g/giswqs/leafmap.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/giswqs/leafmap/context:python)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![image](https://img.shields.io/badge/YouTube-Channel-red)](https://www.youtube.com/c/QiushengWu)
[![image](https://img.shields.io/twitter/follow/giswqs?style=social)](https://twitter.com/giswqs)
[![status](https://joss.theoj.org/papers/10.21105/joss.03414/status.svg)](https://doi.org/10.21105/joss.03414)

**A Python package for geospatial analysis and interactive mapping in a Jupyter environment.**

-   GitHub repo: <https://github.com/giswqs/leafmap>
-   Documentation: <https://leafmap.org>
-   PyPI: <https://pypi.org/project/leafmap>
-   Conda-forge: <https://anaconda.org/conda-forge/leafmap>
-   Leafmap tutorials on YouTube: <https://www.youtube.com/c/QiushengWu>
-   Free software: [MIT license](https://opensource.org/licenses/MIT)

## Introduction

**Leafmap** is a Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment. It is a spin-off project of the [geemap](https://geemap.org) Python package, which was designed specifically to work with [Google Earth Engine](https://earthengine.google.com) (GEE). However, not everyone in the geospatial community has access to the GEE cloud computing platform. Leafmap is designed to fill this gap for non-GEE users. It is a free and open-source Python package that enables users to analyze and visualize geospatial data with minimal coding in a Jupyter environment, such as Google Colab, Jupyter Notebook, and JupyterLab. Leafmap is built upon several open-source packages, such as [folium](https://github.com/python-visualization/folium) and [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) (for creating interactive maps), [WhiteboxTools](https://github.com/jblindsay/whitebox-tools) and [whiteboxgui](https://github.com/giswqs/whiteboxgui) (for analyzing geospatial data), and [ipywidgets](https://github.com/jupyter-widgets/ipywidgets) (for designing interactive graphical user interface [GUI]). Leafmap has a toolset with various interactive tools that allow users to load vector and raster data onto the map without coding. In addition, users can use the powerful analytical backend (i.e., WhiteboxTools) to perform geospatial analysis directly within the leafmap user interface without writing a single line of code. The WhiteboxTools library currently contains **468** tools for advanced geospatial analysis, such as [GIS Analysis](https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html), [Geomorphometric Analysis](https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html), [Hydrological Analysis](https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html), [LiDAR Data Analysis](https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html), [Mathematical and Statistical Analysis](https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html), and [Stream Network Analysis](https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html).

## Statement of Need

There are a plethora of Python packages for geospatial analysis, such as [geopandas](https://github.com/geopandas/geopandas) for vector data analysis and [xarray](https://github.com/pydata/xarray) for raster data analysis. However, few Python packages provide interactive GUIs for loading geospatial data in a Jupyter environment. It might take many lines to code to load and display geospatial data with various file formats on an interactive map, which can be a challenging task for novice users with limited coding skills. There are also some notable Python packages for visualizing geospatial data in a Jupyter environment, such as [plotly](https://github.com/plotly/plotly.py) and [kepler.gl](https://docs.kepler.gl/docs/keplergl-jupyter). However, plotly is designed for displaying static data, which lacks bidirectional communication between the front-end and the backend. Kepler.gl provides unique 3D functionality for visualizing large-scale geospatial datasets, but it lacks tools for performing geospatial analysis, such as hydrological analysis and LiDAR data analysis. In contrast, leafmap provides many convenient functions for loading and visualizing geospatial datasets with only one line of code. Users can also use the interactive GUI to load geospatial datasets without coding. Leafmap is intended for anyone who would like to analyze and visualize geospatial data interactively in a Jupyter environment. It is particularly suited for novice users with limited programming skills. Advanced programmers can also find leafmap a useful tool for analyzing geospatial data and building interactive web apps.

Launch the interactive notebook tutorial for the **leafmap** Python package with Google Colab or Binder now:

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)

Check out this excellent article on Medium - [Leafmap a new Python Package for Geospatial data science](https://link.medium.com/HRRKDcynYgb)

To learn more about leafmap, check out the leafmap documentation website - <https://leafmap.org>

![](https://i.imgur.com/abd8pTH.gif)

## Key Features

Below is a partial list of features available for the leafmap package. Please check the [examples](https://github.com/giswqs/leafmap/tree/master/examples) page for notebook examples, GIF animations, and video tutorials.

-   Create an interactive map with only one-line of code.
-   Select from a variety of basemaps interactively without coding.
-   Add XYZ, WMS, and vector tile services to the map.
-   Convert CSV to points and display points as a marker cluster.
-   Add local vector data (e.g., shapefile, GeoJSON, KML) to the map without coding.
-   Add local raster data (e.g., GeoTIFF) to the map without coding.
-   Add Cloud Optimized GeoTIFF (COG) and SpatialTemporal Asset Catalog (STAC) to the map.
-   Add OpenStreetMap data to the map with a single line of code.
-   Add a GeoPandas GeoDataFrame to the map with a single line of code.
-   Add a point layer with popup attributes to the map.
-   Add data from a PostGIS database to the map.
-   Add custom legends and colorbars to the map.
-   Perform geospatial analysis using WhiteboxTools and whiteboxgui.
-   Create split-panel map and linked maps.
-   Publish interactive maps with a single line of code.
-   Download and display OpenStreetMap data with a single line of code.

## Citations

If you find **leafmap** useful in your research, please consider citing the following paper to support my work. Thank you for your support.

-   Wu, Q. (2021). Leafmap: A Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment. _Journal of Open Source Software_, 6(63), 3414. <https://doi.org/10.21105/joss.03414>

## Demo

![](https://wetlands.io/file/images/leafmap_demo.gif)

## YouTube Channel

I have created a [YouTube Channel](https://www.youtube.com/c/QiushengWu) for sharing geospatial tutorials. You can subscribe to my channel for regular updates. If there is any specific tutorial you would like to see, please submit a feature request [here](https://github.com/giswqs/leafmap/issues).

[![Earth Engine Tutorials on YouTube](https://wetlands.io/file/images/youtube.png)](https://www.youtube.com/c/QiushengWu)
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
giswqs@gmail.com.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
---
title: "Leafmap: A Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment"
tags:
    - Python
    - geospatial
    - ipyleaflet
    - mapping
    - Jupyter
authors:
    - name: Qiusheng Wu
      orcid: 0000-0001-5437-4073
      affiliation: "1"
affiliations:
    - name: Department of Geography, University of Tennessee, Knoxville, TN 37996, United States
      index: 1
date: 22 June 2021
bibliography: paper.bib
---

# Summary

**Leafmap** is a Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment. It is a spin-off project of the [geemap](https://geemap.org) Python package [@Wu2020], which was designed specifically to work with [Google Earth Engine](https://earthengine.google.com) (GEE) [@Gorelick2017]. However, not everyone in the geospatial community has access to the GEE cloud computing platform. Leafmap is designed to fill this gap for non-GEE users. It is a free and open-source Python package that enables users to analyze and visualize geospatial data with minimal coding in a Jupyter environment, such as Google Colab, Jupyter Notebook, and JupyterLab. Leafmap is built upon several open-source packages, such as [folium](https://github.com/python-visualization/folium) [@Filipe2021] and [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) [@Renou2021] (for creating interactive maps), [WhiteboxTools](https://github.com/jblindsay/whitebox-tools) [@Lindsay2018] and [whiteboxgui](https://github.com/giswqs/whiteboxgui) (for analyzing geospatial data), and [ipywidgets](https://github.com/jupyter-widgets/ipywidgets) [@Grout2021] (for designing interactive graphical user interface [GUI]). Leafmap has a toolset with various interactive tools that allow users to load vector and raster data onto the map without coding. In addition, users can use the powerful analytical backend (i.e., WhiteboxTools) to perform geospatial analysis directly within the leafmap user interface without writing a single line of code. The WhiteboxTools library currently contains **468** tools for advanced geospatial analysis, such as [GIS Analysis](https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html), [Geomorphometric Analysis](https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html), [Hydrological Analysis](https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html), [LiDAR Data Analysis](https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html), [Mathematical and Statistical Analysis](https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html), and [Stream Network Analysis](https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html).

# Statement of Need

There is a plethora of Python packages for geospatial analysis, such as [geopandas](https://github.com/geopandas/geopandas) for vector data analysis [@Jordahl2021] and [xarray](https://github.com/pydata/xarray) for raster data analysis [@Hoyer2017]. However, few Python packages provide interactive GUIs for loading geospatial data in a Jupyter environment. It might take many lines to code to load and display geospatial data with various file formats on an interactive map, which can be a challenging task for novice users with limited coding skills. There are also some notable Python packages for visualizing geospatial data in a Jupyter environment, such as [plotly](https://github.com/plotly/plotly.py) [@Mease2021] and [kepler.gl](https://docs.kepler.gl/docs/keplergl-jupyter) [@He2021]. However, plotly is designed for displaying static data, which lacks bidirectional communication between the front-end and the backend. Kepler.gl provides unique 3D functionality for visualizing large-scale geospatial datasets, but it lacks tools for performing geospatial analysis, such as hydrological analysis and LiDAR data analysis. In contrast, leafmap provides many convenient functions for loading and visualizing geospatial datasets with only one line of code. Users can also use the interactive GUI to load geospatial datasets without coding. Leafmap is intended for anyone who would like to analyze and visualize geospatial data interactively in a Jupyter environment. It is particularly suited for novice users with limited programming skills. Advanced programmers can also find leafmap a useful tool for analyzing geospatial data and building interactive web apps.

# Leafmap Plotting Backends

**Leafmap** has three plotting backends, including [folium](https://github.com/python-visualization/folium), [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet), and [here-map-widget-for-jupyter](https://github.com/heremaps/here-map-widget-for-jupyter) [@Kharude2021]. An interactive map created using one of the plotting backends can be displayed in a Jupyter environment, such as Google Colab, Jupyter Notebook, and JupyterLab. By default, `import leafmap` in [Jupyter Notebook](https://gishub.org/leafmap-binder) and [JupyterLab](https://gishub.org/leafmap-binder) will use the `ipyleaflet` plotting backend, whereas `import leafmap` in [Google Colab](https://gishub.org/leafmap-colab) will use the `folium` plotting backend. Note that Google Colab does not yet support custom widgets, such as `ipyleaflet` and `heremap widget` ([source](https://github.com/googlecolab/colabtools/issues/498#issuecomment-695335421)). Therefore, interactive maps created using the `ipyleaflet` and `heremap widget` backends won't show up in Google Colab, even though the code might run successfully without any errors.

The three plotting backends do not offer equal functionality. The `ipyleaflet` plotting backend provides the richest interactive functionality, including the custom toolset for loading, analyzing, and visualizing geospatial data interactively without coding. For example, users can add vector data (e.g., GeoJSON, Shapefile, KML, GeoDataFrame) and raster data (e.g., GeoTIFF, Cloud Optimized GeoTIFF [COG]) to the map with a few clicks (see Figure 1). Users can also perform geospatial analysis using the WhiteboxTools GUI with 468 geoprocessing tools directly within the map interface (see Figure 2). Other interactive functionality (e.g., split-panel map, linked map, time slider, time-series inspector) can also be useful for visualizing geospatial data. The `ipyleaflet` package is built upon `ipywidgets` and allows bidirectional communication between the front-end and the backend enabling the use of the map to capture user input ([source](https://blog.jupyter.org/interactive-gis-in-jupyter-with-ipyleaflet-52f9657fa7a)). In contrast, `folium` has relatively limited interactive functionality. It is meant for displaying static data only. The `folium` plotting backend is included in this package to support using `leafmap` in Google Colab. Note that the aforementioned custom toolset and interactive functionality are not available for the `folium` plotting backend. Compared with `ipyleaflet` and `folium`, the `heremap widget` plotting backend provides some unique [3D functionality](https://github.com/heremaps/here-map-widget-for-jupyter#use-ipywidgets-controls-to-build-an-interactive-gui) for visualizing geospatial data. An [API key](https://developer.here.com/documentation/identity-access-management/dev_guide/topics/dev-apikey.html) from the [Here Developer Portal](https://developer.here.com/) is required to use `heremap`.

![](https://i.imgur.com/pe7CoC7.png)
**Figure 1.** The leafmap user interface built upon ipyleaflet and ipywidgets.

![](https://i.imgur.com/5GzDG3W.png)
**Figure 2.** The WhiteboxTools graphical user interface integrated into leafmap.

# Leafmap Modules

The key functionality of the leafmap Python package is organized into nine modules as shown in the table below.

| Module    | Description                                                                   |
| --------- | ----------------------------------------------------------------------------- |
| basemaps  | A collection of XYZ and WMS tile layers to be used as basemaps                |
| colormaps | Commonly used colormaps and palettes for visualizing geospatial data          |
| common    | Functions being used by multiple plotting backends to process geospatial data |
| foliumap  | A plotting backend built upon the folium Python package                       |
| heremap   | A plotting backend built upon the here-map-widget-for-jupyter                 |
| leafmap   | The default plotting backend built upon the ipyleaflet Python package         |
| legends   | Built-in legends for commonly used geospatial datasets                        |
| osm       | Functions for extracting and downloading OpenStreetMap data                   |
| toolbar   | A custom toolset with interactive tools built upon ipywidgets and ipyleaflet  |

# Leafmap Tutorials

Comprehensive documentation and API reference of the leafmap package is available at https://geemap.org. A list of notebook examples and video tutorials for using leafmap can be found at https://leafmap.org/tutorials. Users can also try out leafmap using Google Colab (https://gishub.org/leafmap-colab) and Binder (https://gishub.org/leafmap-binder) using an Internet browser without having to set up the Python environment and install leafmap on their computer.

# Acknowledgements

The author would like to thank the developers of ipyleaflet and ipywidgets, which empower the interactive mapping functionality of leafmap, including [Martin Renou](https://github.com/martinRenou), [David Brochart](https://github.com/davidbrochart), and [Sylvain Corlay](https://github.com/SylvainCorlay). The authors would also like to express thanks to [John Lindsay](https://github.com/jblindsay) for developing the WhiteboxTools library, which serves as the geospatial analysis backend of leafmap. Special thanks go to all leafmap contributors, especially [Sachin Kharude](https://github.com/sackh) for contributing the heremap plotting backend to leafmap. Last but not least, the author would like to thank [Tomas Beuzen](https://github.com/TomasBeuzen) and [Martin Fleischmann](https://github.com/martinfleis) for reviewing this paper and the `leafmap` package. Their constructive comments greatly improved the quality of the source code and documentation of the `leafmap` package as well as this paper.

# References
---
name: Bug Report
about: Create a bug report to help us improve
labels: bug
---

<!-- Please search existing issues to avoid creating duplicates. -->

### Environment Information

-   leafmap version:
-   Python version:
-   Operating System:

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
---
name: Feature Request
about: Submit a feature request to help us improve
labels: Feature Request
---

<!-- Please search existing issues to avoid creating duplicates. -->

### Description

Describe the feature (e.g., new functions/tutorials) you would like to propose.
Tell us what can be achieved with this new feature and what's the expected outcome.

### Source code

```
Paste your source code here if have sample code to share.
```
# leafmap

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)
[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)
[![image](https://img.shields.io/pypi/v/leafmap.svg)](https://pypi.python.org/pypi/leafmap)
[![image](https://img.shields.io/conda/vn/conda-forge/leafmap.svg)](https://anaconda.org/conda-forge/leafmap)
[![image](https://pepy.tech/badge/leafmap)](https://pepy.tech/project/leafmap)
[![image](https://github.com/giswqs/leafmap/workflows/docs/badge.svg)](https://leafmap.org)
[![image](https://github.com/giswqs/leafmap/workflows/Linux%20build/badge.svg)](https://github.com/giswqs/leafmap/actions)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![image](https://img.shields.io/badge/YouTube-Channel-red)](https://www.youtube.com/c/QiushengWu)
[![image](https://img.shields.io/twitter/follow/giswqs?style=social)](https://twitter.com/giswqs)

## Tutorials

1. Introducing the leafmap Python package for interactive mapping ([video](https://youtu.be/-UPt7x3Gn60) | [gif](https://i.imgur.com/2pRxunR.gif) | [notebook](https://leafmap.org/notebooks/01_leafmap_intro))
2. Using basemaps in leafmap ([video](https://youtu.be/uylpjbDZesY) | [gif](https://youtu.be/-lOo-vxjrDM) | [notebook](https://leafmap.org/notebooks/02_using_basemaps))
3. Using Cloud Optimized GeoTIFF (COG) and SpatioTemporal Asset Catalog (STAC) ([notebook](https://leafmap.org/notebooks/03_cog_stac))
4. Creating a virtual mosaic of Cloud Optimized GeoTIFFs (COG) ([notebook](https://leafmap.org/notebooks/04_cog_mosaic))
5. Loading local raster datasets with leafmap ([notebook](https://leafmap.org/notebooks/05_load_raster))
6. Adding custom legends to the map ([notebook](https://leafmap.org/notebooks/06_legend))
7. Adding custom colorbars to the map ([notebook](https://leafmap.org/notebooks/07_colorbar))
8. Using WhiteboxTools with leafmap ([notebook](https://leafmap.org/notebooks/08_whitebox))
9. Converting CSV to points ([notebook](https://leafmap.org/notebooks/09_csv_to_points))
10. Adding local vector data (e.g., shp, geojson, kml) to the map ([notebook](https://leafmap.org/notebooks/10_add_vector))
11. Creating linked maps for visualizing multiple maps simultaneously ([notebook](https://leafmap.org/notebooks/11_linked_maps))
12. Creating a split-panel map with a single line of code ([notebook](https://leafmap.org/notebooks/12_split_map))
13. Adding a GeoPandas GeoDataFrame to the map with a single line of code ([notebook](https://leafmap.org/notebooks/13_geopandas))
14. Adding data from a PostGIS database to the map ([notebook](https://leafmap.org/notebooks/14_postgis))
15. Downloading OpenStreetMap data to the map with a single line of code ([notebook](https://leafmap.org/notebooks/15_add_osm))
16. Use [HERE Map Widget for Jupyter](https://github.com/heremaps/here-map-widget-for-jupyter) as plotting backend ([notebook](https://leafmap.org/notebooks/16_heremap))
17. Adding vector tile layers to the map ([notebook](https://leafmap.org/notebooks/17_vector_tile_layer))
18. Adding a point layer with popup attributes to the map ([notebook](https://leafmap.org/notebooks/18_point_layer))
19. Saving maps as a html file ([notebook](https://leafmap.org/notebooks/19_map_to_html))
20. Adding Planet global monthly and quarterly mosaic ([notebook](https://leafmap.org/notebooks/20_planet_imagery))
21. Using timeseries inspector with one click ([notebook](https://leafmap.org/notebooks/21_ts_inspector))
22. Using time slider for visualizing timeseries images ([notebook](https://leafmap.org/notebooks/22_time_slider))
23. Creating colormaps with a single line of code ([notebook](https://leafmap.org/notebooks/23_colormaps))
24. Creating heat map from csv ([notebook](https://leafmap.org/notebooks/24_heatmap))
25. Creating a population heat map with a colorbar and map title ([notebook](https://leafmap.org/notebooks/25_map_title))
26. Creating an interactive map using the kepler.gl plotting backend ([notebook](https://leafmap.org/notebooks/26_kepler_gl))
27. Creating a basemap gallery ([notebook](https://leafmap.org/notebooks/27_basemap_gallery))
28. Publishing maps with a single line of code ([notebook](https://leafmap.org/notebooks/28_publish_map))
29. Using the pydeck plotting backend ([notebook](https://leafmap.org/notebooks/29_pydeck))
30. Using U.S. Census data ([notebook](https://leafmap.org/notebooks/30_census_data))
31. Searching basemaps with xyzservices ([notebook](https://leafmap.org/notebooks/31_search_basemaps))
32. Loading local raster datasets and Cloud Optimized GeoTIFF (COG) ([notebook](https://leafmap.org/notebooks/32_local_tile))
33. Adding image overlay to the map ([notebook](https://leafmap.org/notebooks/33_image_overlay))
34. Adding points from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://leafmap.org/notebooks/34_add_points_from_xy))
35. Adding circle markers from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://leafmap.org/notebooks/35_circle_markers))
36. Adding labels to the map ([notebook](https://leafmap.org/notebooks/36_add_labels))
37. Adding Planetary Computer STAC item to the map ([notebook](https://leafmap.org/notebooks/37_planetary_computer))
38. Using the plotly plotting backend ([notebook](https://leafmap.org/notebooks/38_plotly))
39. Getting pixel values using the Inspector tool ([notebook](https://leafmap.org/notebooks/39_inspector_tool))
40. Using the interactive plotly toolbar GUI ([notebook](https://leafmap.org/notebooks/40_plotly_gui))
41. Loading COG/STAC items using the raster GUI ([notebook](https://leafmap.org/notebooks/41_raster_gui))
42. Creating Cloud Optimized GeoTIFF (COG) ([notebook](https://leafmap.org/notebooks/42_create_cog))
43. Searching for locations and features in vector data ([notebook](https://leafmap.org/notebooks/43_search_control))
44. Opening vector data attribute table without coding ([notebook](https://leafmap.org/notebooks/44_attribute_table))
45. Creating vector data interactively without coding ([notebook](https://leafmap.org/notebooks/45_create_vector))
46. Editing existing vector data interactively without coding ([notebook](https://leafmap.org/notebooks/46_edit_vector))

## Demo

![](https://wetlands.io/file/images/leafmap_demo.gif)

## YouTube Channel

I have created a [YouTube Channel](https://www.youtube.com/c/QiushengWu) for sharing geospaital tutorials. You can subscribe to my channel for regular updates. If there is any specific tutorial you would like to see, please submit a feature request [here](https://github.com/giswqs/leafmap/issues).

[![Earth Engine Tutorials on YouTube](https://wetlands.io/file/images/youtube.png)](https://www.youtube.com/c/QiushengWu)
## Data Sources

Some datasets contained in the repository were downloaded from various sources. Credits to the original owners of the datasets.

The following datasets are downloaded from https://github.com/keplergl/kepler.gl/tree/master/bindings/kepler.gl-jupyter/notebooks

-   hex_data.csv
-   hex_config.json
-   sf_zip_geo.json
# toolbar module

::: leafmap.toolbar
# Usage

You can try out leafmap by using Goolge Colab ([![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)) or Binder ([![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)) without having to install anything on your computer.

## Launch Jupyter notebook

```bash
conda activate geo
jupyter notebook
```

## Use ipyleaflet plotting backend

```python
import leafmap
```

## Use folium plotting backend

```python
import leafmap.foliumap as leafmap
```

## Create an interactive map

```python
m = leafmap.Map(center=(40, -100), zoom=4)
m
```

## Customize map height

```python
m = leafmap.Map(height="450px")
m
```

## Set control visibility

```python
m = leafmap.Map(draw_control=False, measure_control=False, fullscreen_control=False, attribution_control=True)
m
```

## Change basemaps

```python
m = leafmap.Map(google_map="TERRAIN")
m.add_basemap("HYBRID")
m
```

## Add XYZ tile layer

```python
m = leafmap.Map()
m.add_tile_layer(url="https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}", name="Google Satellite", attribution="Google")
m
```

## Add WMS tile layer

```python
m = leafmap.Map()
naip_url = 'https://services.nationalmap.gov/arcgis/services/USGSNAIPImagery/ImageServer/WMSServer?'
m.add_wms_layer(url=naip_url, layers='0', name='NAIP Imagery', format='image/png', shown=True)
m
```

## Use HERE Map Widget for Jupyter plotting backend

### Prerequisites

-   A HERE developer account, free and available under [HERE Developer Portal](https://developer.here.com)
-   An [API key](https://developer.here.com/documentation/identity-access-management/dev_guide/topics/dev-apikey.html) from the [HERE Developer Portal](https://developer.here.com)
-   Export API key into environment variable `HEREMAPS_API_KEY`

```bash
export HEREMAPS_API_KEY=YOUR-ACTUAL-API-KEY
```

```python
import leafmap.heremap as leafmap
```

### Create an interactive map

```python
import os
api_key = os.environ.get("HEREMAPS_API_KEY") # read api_key from environment variable.
m = leafmap.Map(api_key=api_key, center=(40, -100), zoom=4)
m
```
# plotlymap module

::: leafmap.plotlymap
# Changelog

## v0.7.5 - January 27, 2022

**New Features**:

-   Added vector creation GUI [#179](https://github.com/giswqs/leafmap/issues/179) [#194](https://github.com/giswqs/leafmap/pull/194)

## v0.7.4 - January 24, 2022

**New Features**:

-   Added attribute table GUI [#179](https://github.com/giswqs/leafmap/issues/179)

**Improvement**

-   Improved add_labels function [#188](https://github.com/giswqs/leafmap/discussions/188)
-   Improved GitHub workflows [#192](https://github.com/giswqs/leafmap/pull/192)
-   Improved add_raster function [#191](https://github.com/giswqs/leafmap/pull/191)
-   Removed nominatim URL from Search Control [#182](https://github.com/giswqs/leafmap/issues/182)
-   Fixed search control bug [#183](https://github.com/giswqs/leafmap/pull/183)

## v0.7.3 - January 21, 2022

**New Features**:

-   Added search control GUI [#182](https://github.com/giswqs/leafmap/issues/182) [#183](https://github.com/giswqs/leafmap/pull/183)
-   Added COG creation [#176](https://github.com/giswqs/leafmap/issues/176)

**Improvement**

-   Removed COG mosaic function #180
-   Updated binder env

## v0.7.2 - January 11, 2022

**New Features**:

-   Added GUI for loading COG/STAC [#164](https://github.com/giswqs/leafmap/issues/164)
-   Added ROI to GeoJSON function [#170](https://github.com/giswqs/leafmap/issues/170)
-   Added add_geojson for plotly [#163](https://github.com/giswqs/leafmap/issues/163) [#167](https://github.com/giswqs/leafmap/pull/167)

## v0.7.1 - January 3, 2022

**New Features**:

-   Added plotly toolbar GUI [#160](https://github.com/giswqs/leafmap/issues/160)
-   Added layer control [#160](https://github.com/giswqs/leafmap/issues/160)
-   Added Inspector support for local tile [#162](https://github.com/giswqs/leafmap/issues/162)
-   Added add_gdf for plotly [#163](https://github.com/giswqs/leafmap/issues/163)

**Improvement**

-   Improved COG visualization [#161](https://github.com/giswqs/leafmap/issues/161)
-   Fixed citation bug [#165](https://github.com/giswqs/leafmap/pull/165)

## v0.7.0 - December 29, 2021

**New Features**:

-   Added Planetary Computer STAC support [#137](https://github.com/giswqs/leafmap/issues/137)
-   Added plotly backend [#109](https://github.com/giswqs/leafmap/issues/109)
-   Added Inspector tool [#158](https://github.com/giswqs/leafmap/pull/158)
-   Added plotly COG STAC support [#109](https://github.com/giswqs/leafmap/issues/109)
-   Added plotly planet imagery support [#109](https://github.com/giswqs/leafmap/issues/109)
-   Added plotly toolbar [#160](https://github.com/giswqs/leafmap/issues/160)
-   Added geojson_to_df and geom_type functions

**Improvement**

-   Removed pangeo broken binder links
-   Improved kepler config options [#150](https://github.com/giswqs/leafmap/discussions/150)
-   Improved stac tile function [#137](https://github.com/giswqs/leafmap/issues/156)
-   Updated STAC notebook example [#156](https://github.com/giswqs/leafmap/issues/156)

## v0.6.1 - December 23, 2021

**New Features**:

-   Added image overlay functionality [#136](https://github.com/giswqs/leafmap/issues/136)
-   Added marker cluster function [#138](https://github.com/giswqs/leafmap/issues/138)
-   Added locate control to folium
-   Added cesium_to_streamlit function [#139](https://github.com/giswqs/leafmap/issues/139)
-   Added add_points_from_xy function [#138](https://github.com/giswqs/leafmap/issues/138)
-   Added circle markers function [#140](https://github.com/giswqs/leafmap/issues/143)

**Improvement**

-   Added localtileserver to env.yml
-   Fixed gdf style callback bug [#119](https://github.com/giswqs/leafmap/issues/119)
-   Added ts_inspector docstring [#147](https://github.com/giswqs/leafmap/discussions/147)
-   Improved streamlit download button

## v0.6.0 - November 27, 2021

**New Features**:

-   Added add_marker function
-   Added save_data function
-   Added support for local tile [#129](https://github.com/giswqs/leafmap/issues/129)
-   Added open raster GUI [#129](https://github.com/giswqs/leafmap/issues/129)
-   Added zoom to tile [#129](https://github.com/giswqs/leafmap/issues/129)

## v0.5.5 - November 9, 2021

**New Features**:

-   Added YouthMappers workshop [notebook](https://leafmap.org/workshops/YouthMappers_2021/)

**Improvement**

-   Fixed `add_legend` bug
-   Changed default `max_zoom` to 24

## v0.5.4 - November 2, 2021

**New Features**:

-   Added search basemaps GUI [#93](https://github.com/giswqs/leafmap/issues/93)
-   Added get wms layers function
-   Made streamlit map width reponsive [#126](https://github.com/giswqs/leafmap/issues/126)
-   Added function read file from url
-   Added streamlit download button
-   Added SIGSPATIAL workshop notebook

**Improvement**

-   Fixed layer attribution error [#93](https://github.com/giswqs/leafmap/issues/93)
-   Fixed open vector bug [#124](https://github.com/giswqs/leafmap/discussions/124)
-   Improved streamlit support

## v0.5.3 - October 17, 2021

**New Features**:

-   Added support for US Census data with hundreds of WMS layers [#123](https://github.com/giswqs/leafmap/issues/123)

## v0.5.2 - October 17, 2021

**Improvement**

-   Fixed pydeck import error

## v0.5.1 - October 17, 2021

**New Features**:

-   Added support for pydeck [#122](https://github.com/giswqs/leafmap/issues/122)
-   Added streamlit support for heremap [#118](https://github.com/giswqs/leafmap/issues/118)
-   Added create_colormap function

**Improvement**

-   Added optional postgis port param [#144](https://github.com/giswqs/leafmap/pull/114)
-   Added STAC time slider example to notebook [#177](https://github.com/giswqs/leafmap/pull/117)
-   Fixed geojson style callback bug [#119](https://github.com/giswqs/leafmap/issues/119)
-   Updated foss4g notebook
-   Fixed planet imagery bug
-   Improved vector to geojson
-   Added streamlit app link to docs

## v0.4.3 - September 17, 2021

**New Features**:

-   Added `sandbox_path` option allowing users to restrict Voila app access to system directories [#113](https://github.com/giswqs/leafmap/issues/113)

## v0.4.2 - September 10, 2021

**New Features**:

-   Changed default plotting backend on Colab from folium to ipyleaflet [#112](https://github.com/giswqs/leafmap/issues/112)
-   Added streamlit support [#96](https://github.com/giswqs/leafmap/issues/96)
-   Added support for xyzservices provider [#92](https://github.com/giswqs/leafmap/issues/92)
-   Added a basemap gallery [#91](https://github.com/giswqs/leafmap/issues/91)

**Improvement**

-   Fixed linked maps bug
-   Improved folium basemaps [#91](https://github.com/giswqs/leafmap/issues/91)

## v0.4.1 - August 4, 2021

**New Features**:

-   Added 200+ basemaps from xyzservices [#91](https://github.com/giswqs/leafmap/issues/91)

**Improvement**

-   Fixed typo [#90](https://github.com/giswqs/leafmap/pull/90)
-   Added kepler module to mkdocs
-   Removed support for Python 3.6 due to xyzservices

## v0.4.0 - July 28, 2021

**New Features**:

-   Added kepler.gl plotting backend [#88](https://github.com/giswqs/leafmap/issues/88)
-   Added keplergl sample data [#88](https://github.com/giswqs/leafmap/issues/88)
-   Added keplergl sample html [#88](https://github.com/giswqs/leafmap/issues/88)

**Improvement**

-   Added CITATIONS.cff

## v0.3.5 - July 26, 2021

**New Features**:

-   Added kepler.gl plotting backend [#88](https://github.com/giswqs/leafmap/issues/88)

**Improvement**

-   Added unittest for toolbar module [#83](https://github.com/giswqs/leafmap/issues/83)
-   Updated paper.md

## v0.3.4 - July 21, 2021

**New Features**:

-   Added map title function [#84](https://github.com/giswqs/leafmap/issues/84)

**Improvement**

-   Improved add_ahp and add_kml for http
-   Added codespell to docs.yml
-   Made XYZ tiles attribution required [#83](https://github.com/giswqs/leafmap/issues/83)
-   Changed some functions to be private [#83](https://github.com/giswqs/leafmap/issues/83)
-   Added more info about plotting backends [#83](https://github.com/giswqs/leafmap/issues/83)
-   Added text description to notebooks [#83](https://github.com/giswqs/leafmap/issues/83)
-   Added NotImplementedError for foliumap [#83](https://github.com/giswqs/leafmap/issues/83)
-   Fixed typos using codespell [#83](https://github.com/giswqs/leafmap/issues/83)
-   Added Code of Conduct [#83](https://github.com/giswqs/leafmap/issues/83)
-   Made usage page interactive [#83](https://github.com/giswqs/leafmap/issues/83)
-   Added key features notebook [#83](https://github.com/giswqs/leafmap/issues/83)
-   Added plotting backend comparison [#83](https://github.com/giswqs/leafmap/issues/83)
-   Added leafmap and foliumap unittest [#83](https://github.com/giswqs/leafmap/issues/83)
-   Improved JOSS paper [#83](https://github.com/giswqs/leafmap/issues/83)

## v0.3.3 - July 8, 2021

**New Features**:

-   Added troubleshooting section [#76](https://github.com/giswqs/leafmap/issues/76)
-   Added pandas_to_geojson function [#75](https://github.com/giswqs/leafmap/issues/75)
-   Added creating heat map from csv [#64](https://github.com/giswqs/leafmap/issues/64)
-   Added cog mosaic from file [#61](https://github.com/giswqs/leafmap/issues/61)
-   Added colormap notebook [#60](https://github.com/giswqs/leafmap/issues/60)

**Improvement**

-   Changed COG and STAC function names [#61](https://github.com/giswqs/leafmap/issues/61)
-   Updated colormap example [#60](https://github.com/giswqs/leafmap/issues/60)

## v0.3.2 - June 22, 2021

**New Features**:

-   Added time slider [#42](https://github.com/giswqs/leafmap/issues/42)
-   Added JOSS manuscript
-   Added unittests

## v0.3.1 - June 20, 2021

**New Features**:

-   Added GUI for loading COG [#50](https://github.com/giswqs/leafmap/issues/50)
-   Added methods to add vector data on heremap [#43 ](https://github.com/giswqs/leafmap/pull/43)
-   Added Planet imagery GUI [#9](https://github.com/giswqs/leafmap/commit/2bea287e08886b8d20b96a80364d898237b425bd)

**Improvement**

-   Improved support for folium styles [#47](https://github.com/giswqs/leafmap/discussions/47)
-   Improved save map to image [#37](https://github.com/giswqs/leafmap/issues/37)
-   Updated toolbar icons [#9](https://github.com/giswqs/leafmap/issues/9)
-   Added LGTM
-   Updated installation docs

## v0.3.0 - June 14, 2021

**New Features**:

-   Added Planet basemaps GUI [#9](https://github.com/giswqs/leafmap/issues/9)
-   Added open point layer GUI [#29](https://github.com/giswqs/leafmap/issues/29)
-   Improved GUI for opening vector data from http [#33](https://github.com/giswqs/leafmap/issues/33)
-   Added map to html function [#32](https://github.com/giswqs/leafmap/issues/32)
-   Added point layer with popup [#27](https://github.com/giswqs/leafmap/issues/27)
-   Added vector tile layer support [#26](https://github.com/giswqs/leafmap/pull/26)
-   Added HERE map plotting backend [#20](https://github.com/giswqs/leafmap/pull/20)

**Improvement**

-   Allow json file in open data widget
-   Added five notebook tutorials
-   Fixed folium map custom size bug [#21](https://github.com/giswqs/leafmap/issues/21)

## v0.2.0 - June 5, 2021

**New Features**:

-   Added handle-draw function [#2](https://github.com/giswqs/leafmap/issues/2)
-   Added split-panel map [#7](https://github.com/giswqs/leafmap/issues/7)
-   Added GeoPandas support [#16](https://github.com/giswqs/leafmap/issues/16)
-   Added support for PostGIS [#15](https://github.com/giswqs/leafmap/issues/15)
-   Added support for downloading OpenStreetMap data [#10](https://github.com/giswqs/leafmap/issues/10) [#12](https://github.com/giswqs/leafmap/issues/12)

**Improvement**

-   Fixed basemap bug [#5](https://github.com/giswqs/leafmap/discussions/5)
-   Fixed output scroll bug [#11](https://github.com/giswqs/leafmap/issues/11)
-   Changed COG and STAC functions to snake_case
-   Added binder badge to notebooks
-   Added binder env
-   Added 15 tutorials
-   Added domain name leafmap.org

## v0.1.0 - May 25, 2021

**New Features**:

-   Create an interactive map with only one-line of code.
-   Select from a variety of basemaps interactively without coding.
-   Add XYZ and WMS tile services to the map.
-   Convert CSV to points and display points as a marker cluster.
-   Add local vector data (e.g., shapefile, GeoJSON, KML) to the map without coding.
-   Add local raster data (e.g., GeoTIFF) to the map without coding.
-   Add Cloud Optimized GeoTIFF (COG) and SpatialTemporal Asset Catalog (STAC) to the map.
-   Add custom legends and colorbars to the map.
-   Perform geospatial analysis using WhiteboxTools and whiteboxgui.
-   Publish interactive maps with only one line of code.
# kepler module

::: leafmap.kepler
# Welcome to leafmap

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)
[![image](https://img.shields.io/pypi/v/leafmap.svg)](https://pypi.python.org/pypi/leafmap)
[![image](https://img.shields.io/conda/vn/conda-forge/leafmap.svg)](https://anaconda.org/conda-forge/leafmap)
[![image](https://pepy.tech/badge/leafmap)](https://pepy.tech/project/leafmap)
[![image](https://github.com/giswqs/leafmap/workflows/docs/badge.svg)](https://leafmap.org)
[![image](https://github.com/giswqs/leafmap/workflows/Linux%20build/badge.svg)](https://github.com/giswqs/leafmap/actions)
[![image](https://img.shields.io/lgtm/grade/python/g/giswqs/leafmap.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/giswqs/leafmap/context:python)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![image](https://img.shields.io/badge/YouTube-Channel-red)](https://www.youtube.com/c/QiushengWu)
[![image](https://img.shields.io/twitter/follow/giswqs?style=social)](https://twitter.com/giswqs)
[![status](https://joss.theoj.org/papers/10.21105/joss.03414/status.svg)](https://doi.org/10.21105/joss.03414)

**A Python package for geospatial analysis and interactive mapping in a Jupyter environment.**

-   GitHub repo: <https://github.com/giswqs/leafmap>
-   Documentation: <https://leafmap.org>
-   PyPI: <https://pypi.org/project/leafmap>
-   Conda-forge: <https://anaconda.org/conda-forge/leafmap>
-   Leafmap tutorials on YouTube: <https://www.youtube.com/c/QiushengWu>
-   Free software: [MIT license](https://opensource.org/licenses/MIT)

## Introduction

**Leafmap** is a Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment. It is a spin-off project of the [geemap](https://geemap.org) Python package, which was designed specifically to work with [Google Earth Engine](https://earthengine.google.com) (GEE). However, not everyone in the geospatial community has access to the GEE cloud computing platform. Leafmap is designed to fill this gap for non-GEE users. It is a free and open-source Python package that enables users to analyze and visualize geospatial data with minimal coding in a Jupyter environment, such as Google Colab, Jupyter Notebook, and JupyterLab. Leafmap is built upon several open-source packages, such as [folium](https://github.com/python-visualization/folium) and [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) (for creating interactive maps), [WhiteboxTools](https://github.com/jblindsay/whitebox-tools) and [whiteboxgui](https://github.com/giswqs/whiteboxgui) (for analyzing geospatial data), and [ipywidgets](https://github.com/jupyter-widgets/ipywidgets) (for designing interactive graphical user interface [GUI]). Leafmap has a toolset with various interactive tools that allow users to load vector and raster data onto the map without coding. In addition, users can use the powerful analytical backend (i.e., WhiteboxTools) to perform geospatial analysis directly within the leafmap user interface without writing a single line of code. The WhiteboxTools library currently contains **468** tools for advanced geospatial analysis, such as [GIS Analysis](https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html), [Geomorphometric Analysis](https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html), [Hydrological Analysis](https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html), [LiDAR Data Analysis](https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html), [Mathematical and Statistical Analysis](https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html), and [Stream Network Analysis](https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html).

## Statement of Need

There are a plethora of Python packages for geospatial analysis, such as [geopandas](https://github.com/geopandas/geopandas) for vector data analysis and [xarray](https://github.com/pydata/xarray) for raster data analysis. However, few Python packages provide interactive GUIs for loading geospatial data in a Jupyter environment. It might take many lines to code to load and display geospatial data with various file formats on an interactive map, which can be a challenging task for novice users with limited coding skills. There are also some notable Python packages for visualizing geospatial data in a Jupyter environment, such as [plotly](https://github.com/plotly/plotly.py) and [kepler.gl](https://docs.kepler.gl/docs/keplergl-jupyter). However, plotly is designed for displaying static data, which lacks bidirectional communication between the front-end and the backend. Kepler.gl provides unique 3D functionality for visualizing large-scale geospatial datasets, but it lacks tools for performing geospatial analysis, such as hydrological analysis and LiDAR data analysis. In contrast, leafmap provides many convenient functions for loading and visualizing geospatial datasets with only one line of code. Users can also use the interactive GUI to load geospatial datasets without coding. Leafmap is intended for anyone who would like to analyze and visualize geospatial data interactively in a Jupyter environment. It is particularly suited for novice users with limited programming skills. Advanced programmers can also find leafmap a useful tool for analyzing geospatial data and building interactive web apps.

Launch the interactive notebook tutorial for the **leafmap** Python package with Google Colab or Binder now:

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)

Check out this excellent article on Medium - [Leafmap a new Python Package for Geospatial data science](https://link.medium.com/HRRKDcynYgb)

## Key Features

Below is a partial list of features available for the leafmap package. Please check the [examples](https://github.com/giswqs/leafmap/tree/master/examples) page for notebook examples, GIF animations, and video tutorials.

-   Create an interactive map with only one-line of code.
-   Select from a variety of basemaps interactively without coding.
-   Add XYZ, WMS, and vector tile services to the map.
-   Convert CSV to points and display points as a marker cluster.
-   Add local vector data (e.g., shapefile, GeoJSON, KML) to the map without coding.
-   Add local raster data (e.g., GeoTIFF) to the map without coding.
-   Add Cloud Optimized GeoTIFF (COG) and SpatialTemporal Asset Catalog (STAC) to the map.
-   Add OpenStreetMap data to the map with a single line of code.
-   Add a GeoPandas GeoDataFrame to the map with a single line of code.
-   Add a point layer with popup attributes to the map.
-   Add data from a PostGIS database to the map.
-   Add custom legends and colorbars to the map.
-   Perform geospatial analysis using WhiteboxTools and whiteboxgui.
-   Create split-panel map and linked maps.
-   Publish interactive maps with a single line of code.
-   Download and display OpenStreetMap data with a single line of code.

## Citations

If you find **leafmap** useful in your research, please consider citing the following paper to support my work. Thank you for your support.

-   Wu, Q. (2021). Leafmap: A Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment. _Journal of Open Source Software_, 6(63), 3414. <https://doi.org/10.21105/joss.03414>

## Demo

![](data/leafmap_demo.gif)

## YouTube Channel

I have created a [YouTube Channel](https://www.youtube.com/c/QiushengWu) for sharing geospatial tutorials. You can subscribe to my channel for regular updates. If there is any specific tutorial you would like to see, please submit a feature request [here](https://github.com/giswqs/leafmap/issues).

[![Leafmap Tutorials on YouTube](https://wetlands.io/file/images/youtube.png)](https://www.youtube.com/c/QiushengWu)
# Installation

## Install from PyPI

**leafmap** is available on [PyPI](https://pypi.org/project/leafmap/). To install **leafmap**, run this command in your terminal:

```bash
    pip install leafmap
```

## Install from conda-forge

**leafmap** is also available on [conda-forge](https://anaconda.org/conda-forge/leafmap). If you have
[Anaconda](https://www.anaconda.com/distribution/#download-section) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your computer, you can install leafmap using the following command:

```bash
    conda install leafmap -c conda-forge
```

The leafmap package has some optional dependencies (e.g., [geopandas](https://geopandas.org/) and [localtileserver](https://github.com/banesullivan/localtileserver)), which can be challenging to install on some computers, especially Windows. It is highly recommended that you create a fresh conda environment to install geopandas and leafmap. Follow the commands below to set up a conda env and isntall [geopandas](https://geopandas.org), [localtileserver](https://github.com/banesullivan/localtileserver), [keplergl](https://docs.kepler.gl/docs/keplergl-jupyter), [pydeck](https://deckgl.readthedocs.io/), and leafmap.

```bash
    conda create -n geo python=3.9
    conda activate geo
    conda install geopandas
    conda install mamba -c conda-forge
    mamba install localtileserver keplergl pydeck leafmap -c conda-forge
```

Optionally, you can install some [Jupyter notebook extensions](https://github.com/ipython-contrib/jupyter_contrib_nbextensions), which can improve your productivity in the notebook environment. Some useful extensions include Table of Contents, Gist-it, Autopep8, Variable Inspector, etc. See this [post](https://towardsdatascience.com/jupyter-notebook-extensions-517fa69d2231) for more information.

```bash
    conda install jupyter_contrib_nbextensions -c conda-forge
```

## Install from GitHub

To install the development version from GitHub using [Git](https://git-scm.com/), run the following command in your terminal:

```bash
    pip install git+https://github.com/giswqs/leafmap
```

## Upgrade leafmap

If you have installed **leafmap** before and want to upgrade to the latest version, you can run the following command in your terminal:

```bash
    pip install -U leafmap
```

If you use conda, you can update leafmap to the latest version by running the following command in your terminal:

```bash
    conda update -c conda-forge leafmap
```

To install the development version from GitHub directly within Jupyter notebook without using Git, run the following code:

```python
    import leafmap
    leafmap.update_package()
```

## Troubleshooting

If the interactive map does not show up on Jupyter Notebook and JupyterLab, it is probably because the [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) extentsion is not installed properly.
For example, you might receive an error message saying `Error displaying widget: model not found`. This a well-known issue related to ipyleaflet. See some relevant issues below.

-   [How to display map object using ipyleaflet in jupyter notebook or jupyter Lab](https://github.com/jupyter-widgets/ipyleaflet/issues/739)
-   [ipyleaflet does not work in jupyter lab - "Error displaying widget: model not found"](https://github.com/jupyter-widgets/ipyleaflet/issues/418)
-   [Error displaying widget: model not found](https://github.com/jupyter-widgets/ipyleaflet/issues/504)

Try some of the options below to resolve the issue. If the issue persists after trying these steps, you can open an issue on the [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet/issues) repository.

For Jupyter notebook, try running the following two commands within your leafmap conda environment:

```
jupyter nbextension install --py --symlink --sys-prefix ipyleaflet
jupyter nbextension enable --py --sys-prefix ipyleaflet
```

For JupyterLab, try running the following command within your leafmap conda environment:

```
jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-leaflet

```

Alternatively, you can run leafmap directly using binder:

-   <https://gishub.org/leafmap-binder>
-   <https://gishub.org/leafmap-binder>
# common module

::: leafmap.common
# basemaps module

::: leafmap.basemaps
# pydeck module

::: leafmap.deck
# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given. You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at <https://github.com/giswqs/leafmap/issues>.

If you are reporting a bug, please include:

-   Your operating system name and version.
-   Any details about your local setup that might be helpful in troubleshooting.
-   Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

### Write Documentation

leafmap could always use more documentation, whether as part of the official leafmap docs, in docstrings, or even on the web in blog posts, articles, and such.

### Submit Feedback

The best way to send feedback is to file an issue at <https://github.com/giswqs/leafmap/issues>.

If you are proposing a feature:

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible, to make it easier to implement.
-   Remember that this is a volunteer-driven project, and that contributions are welcome :)

## Get Started

Ready to contribute? Here's how to set up _leafmap_ for local development.

1. Fork the [leafmap](https://github.com/giswqs/leafmap) repo on GitHub.

2. Clone your fork locally:

    ```
    git clone git@github.com:your_name_here/leafmap.git
    ```

3. Install your local copy into a conda env. Assuming you have conda installed, this is how you set up your fork for local development:

    ```
    conda create -n leafmap-test python
    ```

    ```
    conda activate leafmap-test
    ```

    ```
    cd leafmap/
    ```

    ```
    pip install -e .
    ```

4. Create a branch for local development:

    ```
    git checkout -b name-of-your-bugfix-or-feature
    ```

    Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox:

    ```
    flake8 leafmap tests
    ```

    ```
    python setup.py test or pytest
    ```

    To get flake8 and tox, just pip install them into your conda env.

6. Commit your changes and push your branch to GitHub:

    ```
    git add .
    ```

    ```
    git commit -m "Your detailed description of your changes."
    ```

    ```
    git push origin name-of-your-bugfix-or-feature
    ```

7. Submit a pull request through the GitHub website.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put your new functionality into a function with a docstring, and add the feature to the list in README.rst.
3. The pull request should work for Python 3.6, 3.7 and 3.8, and for PyPy. Check <https://github.com/giswqs/leafmap/actions> and make sure that the tests pass for all supported Python versions.

## Code of Conduct

### Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

### Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

-   Demonstrating empathy and kindness toward other people
-   Being respectful of differing opinions, viewpoints, and experiences
-   Giving and gracefully accepting constructive feedback
-   Accepting responsibility and apologizing to those affected by our mistakes,
    and learning from the experience
-   Focusing on what is best not just for us as individuals, but for the
    overall community

Examples of unacceptable behavior include:

-   The use of sexualized language or imagery, and sexual attention or
    advances of any kind
-   Trolling, insulting or derogatory comments, and personal or political attacks
-   Public or private harassment
-   Publishing others' private information, such as a physical or email
    address, without their explicit permission
-   Other conduct which could reasonably be considered inappropriate in a
    professional setting

### Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

### Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
giswqs@gmail.com.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

### Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

#### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

#### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

#### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

#### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# legends module

::: leafmap.legends
# foliumap module

::: leafmap.foliumap
# heremap module

::: leafmap.heremap
# Get Started

This Get Started guide is intended as a quick way to start programming with **leafmap**. You can try out leafmap by using Goolge Colab ([![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)) or Binder ([![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)) without having to install anything on your computer.

## Important Note

**Leafmap** has three plotting backends, including [folium](https://github.com/python-visualization/folium), [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet), and [here-map-widget-for-jupyter](https://github.com/heremaps/here-map-widget-for-jupyter). An interactive map created using one of the plotting backends can be displayed in a Jupyter environment, such as Google Colab, Jupyter Notebook, and JupyterLab. By default, `import leafmap` in [Jupyter Notebook](https://gishub.org/leafmap-binder) and [JupyterLab](https://gishub.org/leafmap-binder) will use the `ipyleaflet` plotting backend, whereas `import leafmap` in [Google Colab](https://gishub.org/leafmap-colab) will use the `folium` plotting backend. Note that Google Colab does not yet support custom widgets, such as `ipyleaflet` and `heremap widget` ([source](https://github.com/googlecolab/colabtools/issues/498#issuecomment-695335421)). Therefore, interactive maps created using the `ipyleaflet` and `heremap widget` backends won't show up in Google Colab, even though the code might run successfully without any errors.

The three plotting backends do not offer equal functionality. The `ipyleaflet` plotting backend provides the richest interactive functionality, including the custom toolset for loading, analyzing, and visualizing geospatial data interactively without coding. For example, users can add vector data (e.g., GeoJSON, Shapefile, KML, GeoDataFrame) and raster data (e.g., GeoTIFF, Cloud Optimized GeoTIFF [COG]) to the map with a few clicks (see Figure 1). Users can also perform geospatial analysis using the WhiteboxTools GUI with 468 geoprocessing tools directly within the map interface (see Figure 2). Other interactive functionality (e.g., split-panel map, linked map, time slider, time-series inspector) can also be useful for visualizing geospatial data. The `ipyleaflet` package is built upon `ipywidgets` and allows bidirectional communication between the front-end and the backend enabling the use of the map to capture user input ([source](https://blog.jupyter.org/interactive-gis-in-jupyter-with-ipyleaflet-52f9657fa7a)). In contrast, `folium` has relatively limited interactive functionality. It is meant for displaying static data only. The `folium` plotting backend is included in this package to support using `leafmap` in Google Colab. Note that the aforementioned custom toolset and interactive functionality are not available for the `folium` plotting backend. Compared with `ipyleaflet` and `folium`, the `heremap widget` plotting backend provides some unique [3D functionality](https://github.com/heremaps/here-map-widget-for-jupyter#use-ipywidgets-controls-to-build-an-interactive-gui) for visualizing geospatial data. An [API key](https://developer.here.com/documentation/identity-access-management/dev_guide/topics/dev-apikey.html) from the [Here Developer Portal](https://developer.here.com/) is required to use `heremap`.

![](https://i.imgur.com/pe7CoC7.png)
**Figure 1.** The leafmap user interface built upon ipyleaflet and ipywidgets.

![](https://i.imgur.com/5GzDG3W.png)
**Figure 2.** The WhiteboxTools graphical user interface integrated into leafmap.

To use a specific plotting backend, use one of the following:

-   `import leafmap.leafmap as leafmap`
-   `import leafmap.foliumap as leafmap`
-   `import leafmap.heremap as leafmap`

## Leafmap Modules

The key functionality of the leafmap Python package is organized into nine modules as shown in the table below.

| Module    | Description                                                                   |
| --------- | ----------------------------------------------------------------------------- |
| basemaps  | A collection of XYZ and WMS tile layers to be used as basemaps                |
| colormaps | Commonly used colormaps and palettes for visualizing geospatial data          |
| common    | Functions being used by multiple plotting backends to process geospatial data |
| foliumap  | A plotting backend built upon the folium Python package                       |
| heremap   | A plotting backend built upon the here-map-widget-for-jupyter                 |
| leafmap   | The default plotting backend built upon the ipyleaflet Python package         |
| legends   | Built-in legends for commonly used geospatial datasets                        |
| osm       | Functions for extracting and downloading OpenStreetMap data                   |
| toolbar   | A custom toolset with interactive tools built upon ipywidgets and ipyleaflet  |

## Launch Jupyter notebook

```bash
conda activate env_name
jupyter notebook
```

## Use ipyleaflet plotting backend

```python
import leafmap
m = leafmap.Map(center=(40, -100), zoom=4)
m
```

![](https://i.imgur.com/CUtzD19.png)

## Use folium plotting backend

```python
import leafmap.foliumap as leafmap
m = leafmap.Map(center=(40, -100), zoom=4)
m
```

![](https://i.imgur.com/cwdskMb.png)

## Use heremap plotting backend

### Prerequisites

-   A HERE developer account, free and available under [HERE Developer Portal](https://developer.here.com)
-   An [API key](https://developer.here.com/documentation/identity-access-management/dev_guide/topics/dev-apikey.html) from the [HERE Developer Portal](https://developer.here.com)
-   Export API key into environment variable `HEREMAPS_API_KEY`

```bash
export HEREMAPS_API_KEY=YOUR-ACTUAL-API-KEY
```

### Create an interactive map

```python
import os
import leafmap.heremap as leafmap
os.environ["HEREMAPS_API_KEY"] = "YOUR_HEREMAPS_API_KEY"
api_key = os.environ.get("HEREMAPS_API_KEY") # read api_key from environment variable.
m = leafmap.Map(api_key=api_key, center=(40, -100), zoom=4)
m
```

![](https://i.imgur.com/TWfgHsV.png)

## Leafmap demo with ipyleaflet backend

![](data/leafmap_demo.gif)
# osm module

::: leafmap.osm
 
# leafmap module

::: leafmap.leafmap# Interactive maps

This page demonstrates some interactive maps created using the [kepler.gl](https://kepler.gl/) plotting backend.

## Create an interactive map

You can specify various parameters to initialize the map, such as `center`, `zoom`, `height`, and `widescreen`.

```python
import leafmap.kepler as leafmap
m = leafmap.Map(center=[40, -100], zoom=2, height=600, widescreen=False)
m
```

<iframe width=760 height=500 frameBorder=0 src="../html/kepler.html"></iframe>

## Add a CSV

```python
m = leafmap.Map(center=[37.7621, -122.4143], zoom=12)
in_csv = 'https://raw.githubusercontent.com/giswqs/leafmap/master/examples/data/hex_data.csv'
config = 'https://raw.githubusercontent.com/giswqs/leafmap/master/examples/data/hex_config.json'
m.add_csv(in_csv, layer_name="hex_data", config=config)
m
```

<iframe width=760 height=500 frameBorder=0 src="../html/kepler_hex.html"></iframe>

## Add a GeoJSON

```python
m = leafmap.Map(center=[20, 0], zoom=1)
lines = 'https://raw.githubusercontent.com/giswqs/leafmap/master/examples/data/cable-geo.geojson'
m.add_geojson(lines, layer_name="Cable lines")
m
```

<iframe width=760 height=500 frameBorder=0 src="../html/kepler_lines.html"></iframe>

Add a GeoJSON with US state boundaries to the map.

```python
m = leafmap.Map(center=[50, -110], zoom=2)
polygons = 'https://raw.githubusercontent.com/giswqs/leafmap/master/examples/data/us-states.json'
m.add_geojson(polygons, layer_name="Countries")
m
```

<iframe width=760 height=500 frameBorder=0 src="../html/kepler_states.html"></iframe>

## Add a shapefile

```python
m = leafmap.Map(center=[20, 0], zoom=1)
in_shp = "https://github.com/giswqs/leafmap/raw/master/examples/data/countries.zip"
m.add_shp(in_shp, "Countries")
m
```

<iframe width=760 height=500 frameBorder=0 src="../html/kepler_countries.html"></iframe>

## Add a GeoDataFrame

```python
import geopandas as gpd
gdf = gpd.read_file("https://raw.githubusercontent.com/giswqs/leafmap/master/examples/data/world_cities.geojson")
m = leafmap.Map(center=[20, 0], zoom=1)
m.add_gdf(gdf, "World cities")
m
```

<iframe width=760 height=500 frameBorder=0 src="../html/kepler_cities.html"></iframe>
# colormaps module

::: leafmap.colormaps
# pc module

::: leafmap.pc
# Tutorials

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/leafmap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/leafmap-binder)

1. Introducing the leafmap Python package for interactive mapping ([video](https://youtu.be/-UPt7x3Gn60) | [gif](https://i.imgur.com/2pRxunR.gif) | [notebook](https://leafmap.org/notebooks/01_leafmap_intro))
2. Using basemaps in leafmap ([video](https://youtu.be/uylpjbDZesY) | [gif](https://youtu.be/-lOo-vxjrDM) | [notebook](https://leafmap.org/notebooks/02_using_basemaps))
3. Using Cloud Optimized GeoTIFF (COG) and SpatioTemporal Asset Catalog (STAC) ([notebook](https://leafmap.org/notebooks/03_cog_stac))
4. Creating a virtual mosaic of Cloud Optimized GeoTIFFs (COG) ([notebook](https://leafmap.org/notebooks/04_cog_mosaic))
5. Loading local raster datasets with leafmap ([notebook](https://leafmap.org/notebooks/05_load_raster))
6. Adding custom legends to the map ([notebook](https://leafmap.org/notebooks/06_legend))
7. Adding custom colorbars to the map ([notebook](https://leafmap.org/notebooks/07_colorbar))
8. Using WhiteboxTools with leafmap ([notebook](https://leafmap.org/notebooks/08_whitebox))
9. Converting CSV to points ([notebook](https://leafmap.org/notebooks/09_csv_to_points))
10. Adding local vector data (e.g., shp, geojson, kml) to the map ([notebook](https://leafmap.org/notebooks/10_add_vector))
11. Creating linked maps for visualizing multiple maps simultaneously ([notebook](https://leafmap.org/notebooks/11_linked_maps))
12. Creating a split-panel map with only one line of code ([notebook](https://leafmap.org/notebooks/12_split_map))
13. Adding a GeoPandas GeoDataFrame to the map with a single line of code ([notebook](https://leafmap.org/notebooks/13_geopandas))
14. Adding data from a PostGIS database to the map ([notebook](https://leafmap.org/notebooks/14_postgis))
15. Downloading OpenStreetMap data with a single line of code ([notebook](https://leafmap.org/notebooks/15_openstreetmap))
16. Using [HERE Map Widget for Jupyter](https://github.com/heremaps/here-map-widget-for-jupyter) as a plotting backend ([notebook](https://leafmap.org/notebooks/16_heremap))
17. Adding vector tile layers to the map ([notebook](https://leafmap.org/notebooks/17_vector_tile_layer))
18. Adding a point layer with popup attributes to the map ([notebook](https://leafmap.org/notebooks/18_point_layer))
19. Saving maps as a html file ([notebook](https://leafmap.org/notebooks/19_map_to_html))
20. Adding Planet global monthly and quarterly mosaic ([notebook](https://leafmap.org/notebooks/20_planet_imagery))
21. Using timeseries inspector with a single click ([notebook](https://leafmap.org/notebooks/21_ts_inspector))
22. Using time slider for visualizing timeseries images ([notebook](https://leafmap.org/notebooks/22_time_slider))
23. Creating colormaps with a single line of code ([notebook](https://leafmap.org/notebooks/23_colormaps))
24. Creating heat map from csv ([notebook](https://leafmap.org/notebooks/24_heatmap))
25. Creating a population heat map with a colorbar and map title ([notebook](https://leafmap.org/notebooks/25_map_title))
26. Creating an interactive map using the kepler.gl plotting backend ([notebook](https://leafmap.org/notebooks/26_kepler_gl))
27. Creating a basemap gallery ([notebook](https://leafmap.org/notebooks/27_basemap_gallery))
28. Publishing maps with a single line of code ([notebook](https://leafmap.org/notebooks/28_publish_map))
29. Using the pydeck plotting backend ([notebook](https://leafmap.org/notebooks/29_pydeck))
30. Using U.S. Census data ([notebook](https://leafmap.org/notebooks/30_census_data))
31. Searching basemaps with xyzservices ([notebook](https://leafmap.org/notebooks/31_search_basemaps))
32. Loading local raster datasets and Cloud Optimized GeoTIFF (COG) ([notebook](https://leafmap.org/notebooks/32_local_tile))
33. Adding image overlay to the map ([notebook](https://leafmap.org/notebooks/33_image_overlay))
34. Adding points from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://leafmap.org/notebooks/34_add_points_from_xy))
35. Adding circle markers from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://leafmap.org/notebooks/35_circle_markers))
36. Adding labels to the map ([notebook](https://leafmap.org/notebooks/36_add_labels))
37. Adding Planetary Computer STAC item to the map ([notebook](https://leafmap.org/notebooks/37_planetary_computer))
38. Using the plotly plotting backend ([notebook](https://leafmap.org/notebooks/38_plotly))
39. Getting pixel values using the Inspector tool ([notebook](https://leafmap.org/notebooks/39_inspector_tool))
40. Using the interactive plotly toolbar GUI ([notebook](https://leafmap.org/notebooks/40_plotly_gui))
41. Loading COG/STAC items using the raster GUI ([notebook](https://leafmap.org/notebooks/41_raster_gui))
42. Creating Cloud Optimized GeoTIFF (COG) ([notebook](https://leafmap.org/notebooks/42_create_cog))
43. Searching for locations and features in vector data ([notebook](https://leafmap.org/notebooks/43_search_control))
44. Opening vector data attribute table without coding ([notebook](https://leafmap.org/notebooks/44_attribute_table))
45. Creating vector data interactively without coding ([notebook](https://leafmap.org/notebooks/45_create_vector))
46. Editing existing vector data interactively without coding ([notebook](https://leafmap.org/notebooks/46_edit_vector))

## Demo

![](https://wetlands.io/file/images/leafmap_demo.gif)
# FAQ

## How do I report an issue or make a feature request

Please go to <https://github.com/giswqs/leafmap/issues>.

## What's the difference between folium and ipyleaflet

A key difference between [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) and [folium](https://github.com/python-visualization/folium) is that ipyleaflet is built upon ipywidgets and allows bidirectional communication between the front-end and the backend enabling the use of the map to capture user input, while folium is meant for displaying static data only ([source](https://blog.jupyter.org/interactive-gis-in-jupyter-with-ipyleaflet-52f9657fa7a)). Note that [Google Colab](https://colab.research.google.com/) currently does not support ipyleaflet ([source](https://github.com/googlecolab/colabtools/issues/498#issuecomment-695335421)). Therefore, if you are using leafmap
with Google Colab, `import leafmap` will automatically use the `folium` plotting backend. If you are using leafmap with Jupyter installed locally, `import leafmap` will automatically use the `ipyleaflet', which provides more functionalities for capturing user input (e.g., mouse-clicking and moving).

## How to use a specific plotting backend

`leafmap` has three plotting backends: [folium](https://github.com/python-visualization/folium), [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet), and [here-map-widget-for-jupyter](https://github.com/heremaps/here-map-widget-for-jupyter). If you are using `leafmap` with Jupyter installed locally, `import leafmap` will use the `ipyleaflet` plotting backend by default. If you are using `leafmap` with [Google Colab](https://githubtocolab.com/giswqs/leafmap/blob/master/examples/notebooks/01_leafmap_intro.ipynb), `import leafmap` will use the `folium` plotting backend by default. Note that Google Colab does not yet support `ipyleaflet` ([source](https://github.com/googlecolab/colabtools/issues/498#issuecomment-695335421)). Therefore, you won't be able to access the `leafmap` toolbar in Colab. Note that the backends do not offer equal functionality. Some interactive functionality in `ipyleaflet` might not be available in `folium` or `heremap`. To use a specific plotting backend, use one of the following:

-   `import leafmap.leafmap as leafmap`
-   `import leafmap.foliumap as leafmap`
-   `import leafmap.heremap as leafmap`

## Why the interactive map does not show up

If the interactive map does not show up on Jupyter Notebook and JupyterLab, it is probably because the [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) extension is not installed properly.
For example, you might receive an error message saying `Error displaying widget: model not found`. This a well-known issue related to ipyleaflet. See some relevant issues below.

-   [How to display map object using ipyleaflet in jupyter notebook or jupyter Lab](https://github.com/jupyter-widgets/ipyleaflet/issues/739)
-   [ipyleaflet does not work in jupyter lab - "Error displaying widget: model not found"](https://github.com/jupyter-widgets/ipyleaflet/issues/418)
-   [Error displaying widget: model not found](https://github.com/jupyter-widgets/ipyleaflet/issues/504)

Try some of the options below to resolve the issue. If the issue persists after trying these steps, you can open an issue on the [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet/issues) repository.

For Jupyter notebook, try running the following two commands within your leafmap conda environment:

```
jupyter nbextension install --py --symlink --sys-prefix ipyleaflet
jupyter nbextension enable --py --sys-prefix ipyleaflet
```

For JupyterLab, try running the following command within your leafmap conda environment:

```
jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-leaflet

```

Alternatively, you can run leafmap directly using binder:

-   <https://gishub.org/leafmap-binder>
-   <https://gishub.org/leafmap-binder>

## How to use leafmap in countries where Google Services are blocked

If you are trying to use leafmap in countries where Google Services are blocked (e.g., China), you will need a VPN. Use `leafmap.set_proxy(port=your-port-number)` to connect to Google servers. Otherwise, you might encounter a connection timeout issue.

```python
import leafmap
leafmap.set_proxy(port=your-port-number)
m = leafmap.Map()
m
```
## Data Sources

Some datasets contained in the repository were downloaded from various sources. Credits to the original owners of the datasets.

The following datasets are downloaded from https://github.com/keplergl/kepler.gl/tree/master/bindings/kepler.gl-jupyter/notebooks

-   hex_data.csv
-   hex_config.json
-   sf_zip_geo.json
