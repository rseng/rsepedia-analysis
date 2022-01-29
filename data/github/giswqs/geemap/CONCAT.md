# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at qwu18@utk.edu. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
---
title: "geemap: A Python package for interactive mapping with Google Earth Engine"
tags:
    - Python
    - Google Earth Engine
    - ipyleaflet
    - mapping
    - Jupyter notebook
authors:
    - name: Qiusheng Wu
      orcid: 0000-0001-5437-4073
      affiliation: "1"
affiliations:
    - name: Department of Geography, University of Tennessee, Knoxville, TN 37996, United States
      index: 1
date: 21 May 2020
bibliography: paper.bib
---

# Summary

**geemap** is a Python package for interactive mapping with [Google Earth Engine](https://earthengine.google.com/) (GEE),
which is a cloud computing platform with a [multi-petabyte catalog](https://developers.google.com/earth-engine/datasets/)
of satellite imagery and geospatial datasets (e.g., Landsat, Sentinel, MODIS, NAIP) [@Gorelick2017]. During the past few years,
GEE has become very popular in the geospatial community and it has empowered numerous environmental applications at local, regional,
and global scales. Some of the notable environmental applications include mapping global forest change [@Hansen2013],
global urban change [@Liu2020], global surface water change [@Pekel2016], wetland inundation dynamics [@Wu2019], vegetation
phenology [@Li2019], and time series analysis [@Kennedy2018].

GEE provides both JavaScript and Python APIs for making computational requests to the Earth Engine servers.
Compared with the comprehensive [documentation](https://earthengine.google.com/) and interactive IDE
(i.e., [GEE JavaScript Code Editor](https://code.earthengine.google.com/)) of the GEE JavaScript API, the
GEE Python API lacks good documentation and lacks functionality for visualizing results interactively.
The **geemap** Python package is created to fill this gap. It is built upon
[ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) and [ipywidgets](https://github.com/jupyter-widgets/ipywidgets),
enabling GEE users to analyze and visualize Earth Engine datasets interactively with Jupyter notebooks.

# geemap Audience

**geemap** is intended for students and researchers who would like to utilize the Python ecosystem of diverse libraries and
tools to explore Google Earth Engine. It is also designed for existing GEE users who would like to transition from the GEE
JavaScript API to a Python API. The automated JavaScript-to-Python [conversion module](https://github.com/giswqs/geemap/blob/master/geemap/conversion.py)
of the **geemap** package can greatly reduce the time needed to convert existing GEE JavaScripts to Python scripts and Jupyter notebooks.

# geemap Functionality

The interactive mapping functionality of the **geemap** package is built upon [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) and [folium](https://github.com/python-visualization/folium), both of which rely on Jupyter notebooks for creating interactive maps. A key difference between ipyleaflet and folium is that ipyleaflet is built upon ipywidgets and allows bidirectional communication between the frontend and the backend, enabling the use of the map to capture user input, while folium is meant for displaying static data only [@QuantStack2019]. It should be noted that [Google Colab](https://colab.research.google.com/) currently does not support ipyleaflet. Therefore, if one wants to use **geemap** on Google Colab, one should `import geemap.foliumap as geemap`, which provides limited interactive mapping functionality. To utilize the full interactive mapping functionality of **geemap**, one should `import geemap` on a local computer or secured server with Jupyter notebook installed.

The key functionality of **geemap** is organized into several modules:

-   [geemap](https://geemap.readthedocs.io/en/latest/source/geemap.html#module-geemap.geemap): the main module for interactive mapping with Google Earth Engine, ipyleaflet, and ipywidgets.

-   [foliumap](https://geemap.readthedocs.io/en/latest/source/geemap.html#module-geemap.foliumap): a module for interactive mapping with Earth Engine and folium. It is designed for users to run geemap with Google Colab.

-   [conversion](https://geemap.readthedocs.io/en/latest/source/geemap.html#module-geemap.conversion): utilities for automatically converting Earth Engine JavaScripts to Python scripts and Jupyter notebooks.

-   [basemaps](https://geemap.readthedocs.io/en/latest/source/geemap.html#module-geemap.basemaps): a module for adding various XYZ and WMS tiled basemaps.

-   [legends](https://geemap.readthedocs.io/en/latest/source/geemap.html#module-geemap.legends): a module for adding customized legends to interactive maps.

# geemap Tutorials

Various tutorials and documentation are available for using **geemap**, including:

-   [20+ video tutorials with corresponding notebook examples](https://github.com/giswqs/geemap/tree/master/examples)
-   [360+ Jupyter notebook examples for using Google Earth Engine](https://github.com/giswqs/earthengine-py-notebooks)
-   [Complete documentation on geemap modules and methods](https://geemap.readthedocs.io/en/latest/source/geemap.html)

# Acknowledgements

The author would like to thank the developers of ipyleaflet and ipywidgets, which empower the interactive mapping functionality of **geemap**, including [Martin Renou](https://github.com/martinRenou), [Sylvain Corlay](https://github.com/SylvainCorlay), and [David Brochart](https://github.com/davidbrochart). The author would also like to acknowledge source code contributions from [Justin Braaten](https://github.com/jdbcode), [Cesar Aybar](https://github.com/csaybar), [Oliver Burdekin](https://github.com/Ojaybee), [Diego Garcia Diaz](https://github.com/Digdgeo), and [Stephan Büttig](https://twitter.com/stephan_buettig).

# References
---
name: Bug Report
about: Create a bug report to help us improve
labels: bug
---

<!-- Please search existing issues to avoid creating duplicates. -->

### Environment Information

-   geemap version:
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
# geemap

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/geemap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/geemap-binder)
[![image](https://img.shields.io/pypi/v/geemap.svg)](https://pypi.python.org/pypi/geemap)
[![image](https://img.shields.io/conda/vn/conda-forge/geemap.svg)](https://anaconda.org/conda-forge/geemap)
[![image](https://pepy.tech/badge/geemap)](https://pepy.tech/project/geemap)
[![image](https://readthedocs.org/projects/geemap/badge/?version=latest)](https://geemap.org/geemap)
[![image](https://img.shields.io/badge/YouTube-GEE%20Tutorials-red)](https://gishub.org/geemap)
[![image](https://img.shields.io/twitter/follow/giswqs?style=social)](https://twitter.com/giswqs)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## YouTube Channel

More video tutorials for geemap and Earth Engine are available on my [YouTube channel](https://www.youtube.com/c/QiushengWu). If you can't access YouTube in your country, you can try [西瓜视频](http://gishub.org/xigua) or [哔哩哔哩](https://space.bilibili.com/527404442)。

[![Earth Engine Tutorials on YouTube](https://wetlands.io/file/images/youtube.png)](https://www.youtube.com/c/QiushengWu)

## Tutorials

1. [Introducing the geemap Python package for interactive mapping with Google Earth Engine](#1-introducing-the-geemap-python-package-for-interactive-mapping-with-google-earth-engine) ([video](https://youtu.be/h0pz3S6Tvx0) | [gif](https://i.imgur.com/pI39k7v.gif) | [notebook](https://geemap.org/notebooks/01_geemap_intro))
2. [Using basemaps in geemap and ipyleaflet for interactive mapping with Google Earth Engine](#2-using-basemaps-in-geemap-and-ipyleaflet-for-interactive-mapping-with-google-earth-engine) ([video](https://youtu.be/6J5ZCIUPXfI) | [gif](https://i.imgur.com/P5B2f7p.gif) | [notebook](https://geemap.org/notebooks/02_using_basemaps))
3. [Introducing the Inspector tool for Earth Engine Python API](#3-introducing-the-inspector-tool-for-earth-engine-python-api) ([video](https://youtu.be/k477ksjkaXw) | [gif](https://i.imgur.com/8d77gtI.gif) | [notebook](https://geemap.org/notebooks/03_inspector_tool))
4. [Creating a split-panel map for visualizing Earth Engine data](#4-creating-a-split-panel-map-for-visualizing-earth-engine-data) ([video](https://youtu.be/9EUTX8j-YVM) | [gif](https://i.imgur.com/kql7pC3.gif) | [notebook](https://geemap.org/notebooks/04_split_panel_map))
5. [Using drawing tools to interact with Earth Engine data](#5-using-drawing-tools-to-interact-with-earth-engine-data) ([video](https://youtu.be/N7rK2aV1R4c) | [gif](https://i.imgur.com/Lm5pDUr.gif) | [notebook](https://geemap.org/notebooks/05_drawing_tools))
6. [Creating an interactive map with a marker cluster](#6-creating-an-interactive-map-with-a-marker-cluster) ([video](https://youtu.be/4HycJPrwpuo) | [gif](https://i.imgur.com/GF4cOqh.gif) | [notebook](https://geemap.org/notebooks/06_marker_cluster))
7. [Converting data formats between GeoJSON and Earth Engine](#7-converting-data-formats-between-geojson-and-earth-engine) ([video](https://youtu.be/DbK_SRgrCHw) | [gif](https://i.imgur.com/hVPmUG1.gif) | [notebook](https://geemap.org/notebooks/07_geojson))
8. [Automated conversion from Earth Engine JavaScripts to Python scripts and Jupyter notebooks](#8-automated-conversion-from-earth-engine-javascripts-to-python-scripts-and-jupyter-notebooks) ([video](https://youtu.be/RpIaalFk4H8) | [gif](https://i.imgur.com/BW0zJnN.gif) | [notebook](https://geemap.org/notebooks/08_ee_js_to_ipynb))
9. [Interactive plotting of Earth Engine data with minimal coding](#9-interactive-plotting-of-earth-engine-data-with-minimal-coding) ([video](https://youtu.be/PDab8mkAFL0) | [gif](https://i.imgur.com/iGMRnRb.gif) | [notebook](https://geemap.org/notebooks/09_plotting))
10. [Using shapefiles with Earth Engine without having to upload data to GEE](#10-using-shapefiles-with-earth-engine-without-having-to-upload-data-to-gee) ([video](https://youtu.be/OlNlqfj4uHo) | [gif](https://i.imgur.com/W3vNSdX.gif) | [notebook](https://geemap.org/notebooks/10_shapefiles))
11. [Exporting Earth Engine Image and ImageCollection as GeoTIFF and Numpy array](#11-exporting-earth-engine-image-and-imagecollection-as-geotiff-and-numpy-array) ([video](https://youtu.be/_6JOA-iiEGU) | [gif](https://i.imgur.com/RonLr0j.gif) | [notebook](https://geemap.org/notebooks/11_export_image))
12. [Computing zonal statistics with Earth Engine and exporting results as CSV or shapefile](#12-computing-zonal-statistics-with-earth-engine-and-exporting-results-as-csv-or-shapefile) ([video](https://youtu.be/ou-Xm3CLitM) | [gif](https://i.imgur.com/8xmUitW.gif) | [notebook](https://geemap.org/notebooks/12_zonal_statistics))
13. [Calculating zonal statistics by group (e.g., analyzing land cover composition of each country/state)](#13-calculating-zonal-statistics-by-group-eg-analyzing-land-cover-composition-of-each-countrystate) ([video](https://youtu.be/cORcGGH03gg) | [gif](https://i.imgur.com/LxD2em9.gif) | [notebook](https://geemap.org/notebooks/13_zonal_statistics_by_group))
14. [Adding a customized legend for Earth Engine data](#14-adding-a-customized-legend-for-earth-engine-data) ([video](https://youtu.be/NwnW_qOkNRw) | [gif](https://i.imgur.com/idkZHQp.gif) | [notebook](https://geemap.org/notebooks/14_legends))
15. [Converting Earth Engine JavaScripts to Python code directly within Jupyter notebook](#15-converting-earth-engine-javascripts-to-python-code-directly-within-jupyter-notebook) ([video](https://youtu.be/nAzZjKKd4w0) | [gif](https://i.imgur.com/aGCBWSV.gif) | [notebook](https://geemap.org/notebooks/15_convert_js_to_py))
16. [Adding animated text to GIF images generated from Earth Engine data](#16-adding-animated-text-to-gif-images-generated-from-earth-engine-data) ([video](https://youtu.be/fDnDVuM_Ke4) | [gif](https://i.imgur.com/MSde1om.gif) | [notebook](https://geemap.org/notebooks/16_add_animated_text))
17. [Adding colorbar and images to GIF animations generated from Earth Engine data](#17-adding-colorbar-and-images-to-gif-animations-generated-from-earth-engine-data) ([video](https://youtu.be/CpT3LQPNKJs) | [gif](https://i.imgur.com/13tFMSI.gif) | [notebook](https://geemap.org/notebooks/17_add_colorbar_to_gif))
18. [Creating Landsat timelapse animations with animated text using Earth Engine](#18-creating-landsat-timelapse-animations-with-animated-text-using-earth-engine) ([video](https://youtu.be/OwjSJnGWKJs) | [gif](https://i.imgur.com/XOHOeXk.gif) | [notebook](https://geemap.org/notebooks/18_create_landsat_timelapse))
19. [How to search and import datasets from Earth Engine Data Catalog](#19-how-to-search-and-import-datasets-from-earth-engine-data-catalog) ([video](https://youtu.be/lwtgzrHrXj8) | [gif](https://i.imgur.com/E09p64F.gif) | [notebook](https://geemap.org/notebooks/19_search_places_and_datasets))
20. [Using timeseries inspector to visualize landscape changes over time](#20-using-timeseries-inspector-to-visualize-landscape-changes-over-time) ([video](https://youtu.be/0CZ7Aj8hCyo) | [gif](https://i.imgur.com/61wbRjK.gif) | [notebook](https://geemap.org/notebooks/20_timeseries_inspector))
21. [Exporting Earth Engine maps as HTML files and PNG images](#21-exporting-earth-engine-maps-as-html-files-and-png-images) ([video](https://youtu.be/GWMvaNQz3kY) | [gif](https://i.imgur.com/rJuXH4a.gif) | [notebook](https://geemap.org/notebooks/21_export_map_to_html_png))
22. [How to import Earth Engine Python scripts into Jupyter notebook?](#22-how-to-import-earth-engine-python-scripts-into-jupyter-notebook) ([video](https://youtu.be/V7CbB9W41w8) | [gif](https://i.imgur.com/WwJoBHF.gif) | [notebook](https://geemap.org/notebooks/22_import_scripts))
23. [How to search Earth Engine API and import assets from GEE personal account?](#23-how-to-search-earth-engine-api-and-import-assets-from-gee-personal-account) ([video](https://youtu.be/c9VJ_uRYSkw) | [gif](https://i.imgur.com/b1auzkr.gif) | [notebook](https://geemap.org/notebooks/22_import_assets))
24. [How to publish interactive Earth Engine maps?](#24-how-to-publish-interactive-earth-engine-maps) ([video](https://youtu.be/NNrrLBIqroY) | [gif](https://i.imgur.com/Hpfzazk.gif) | [notebook](https://geemap.org/notebooks/24_publish_maps))
25. [How to load local raster datasets with geemap?](#25-how-to-load-local-raster-datasets-with-geemap) ([video](https://youtu.be/6XIehAnoazk) | [gif](https://i.imgur.com/nsqEt2O.gif) | [notebook](https://geemap.org/notebooks/25_load_rasters))
26. [How to create and deploy Earth Engine Apps using Python?](https://i.imgur.com/Hpfzazk.gif) ([video](https://youtu.be/nsIjfD83ggA) | [gif](https://i.imgur.com/Hpfzazk.gif) | [notebook](https://geemap.org/notebooks/26_heroku))
27. [How to create an interactive Earth Engine App for creating Landsat timelapse?](https://i.imgur.com/doHfnKp.gif) ([video](https://youtu.be/whIXudC6r_s) | [gif](https://i.imgur.com/doHfnKp.gif) | [notebook](https://geemap.org/notebooks/27_timelapse_app))
28. [How to use your local computer as a web server for hosting Earth Engine Apps?](https://i.imgur.com/q0sJSyi.gif) ([video](https://youtu.be/eRDZBVJcNCk) | [gif](https://i.imgur.com/q0sJSyi.gif) | [notebook](https://geemap.org/notebooks/28_voila))
29. [How to use pydeck for rendering Earth Engine data](https://i.imgur.com/HjFB95l.gif) ([video](https://youtu.be/EIkEH4okFF4) | [gif](https://i.imgur.com/HjFB95l.gif) | [notebook](https://geemap.org/notebooks/29_pydeck))
30. [How to get image basic properties and descriptive statistics](https://i.imgur.com/3B6YhkI.gif) ([video](https://youtu.be/eixBPPWgWs8) | [gif](https://i.imgur.com/3B6YhkI.gif) | [notebook](https://geemap.org/notebooks/30_image_props_stats))
31. [Machine Learning with Earth Engine - Unsupervised Classification](https://i.imgur.com/uNQfrFx.gif) ([video](https://youtu.be/k9MEy2awVJQ) | [gif](https://i.imgur.com/uNQfrFx.gif) | [notebook](https://geemap.org/notebooks/31_unsupervised_classification))
32. [Machine Learning with Earth Engine - Supervised Classification](https://i.imgur.com/jJ2Xiu6.gif) ([video](https://youtu.be/qWaEfgWi21o) | [gif](https://i.imgur.com/jJ2Xiu6.gif) | [notebook](https://geemap.org/notebooks/32_supervised_classification))
33. [Machine Learning with Earth Engine - Performing Accuracy Assessment for Image Classification](https://i.imgur.com/1JkIrF3.gif) ([video](https://youtu.be/JYptiw-I8dc) | [gif](https://i.imgur.com/1JkIrF3.gif) | [notebook](https://geemap.org/notebooks/33_accuracy_assessment))
34. [Interactive extraction of pixel values and interactive region reduction](https://i.imgur.com/LXRqSTu.gif) ([video](https://t.co/D0NC63KgF3) | [gif](https://i.imgur.com/LXRqSTu.gif) | [notebook](https://geemap.org/notebooks/34_extract_values))
35. How to use geemap and Earth Engine in Google Colab ([video](https://youtu.be/fG6kx9vq7hs) | [gif](https://i.imgur.com/OJCasMe.gif) | [notebook](https://geemap.org/notebooks/35_geemap_colab))
36. How to find out the greenest day of the year ([video](https://youtu.be/9KEaW4Ks5fQ) | [gif](https://i.imgur.com/eLDeb4t.gif) | [notebook](https://geemap.org/notebooks/36_quality_mosaic))
37. How to use Earth Engine with pydeck for 3D terrain visualization ([video](https://youtu.be/4E3zOP3-md8) | [gif](https://i.imgur.com/Gx7Y015.gif) | [notebook](https://geemap.org/notebooks/37_pydeck_3d))
38. How to use Cloud Optimized GeoTIFF with Earth Engine ([video](https://youtu.be/2P2PGSMj-wM) | [gif](https://i.imgur.com/z2mfrrZ.gif) | [notebook](https://geemap.org/notebooks/38_cloud_geotiff))
39. How to create Landsat timelapse animations without coding ([video](https://youtu.be/ab0oUhnd_7U) | [gif](https://i.imgur.com/7eyMcZQ.gif) | [notebook](https://geemap.org/notebooks/39_timelapse))
40. How to add interactive widgets to the map ([video](https://youtu.be/KsIxGq6cHtw) | [gif](https://i.imgur.com/peRZZjj.gif) | [notebook](https://geemap.org/notebooks/40_ipywidgets))
41. How to develop an Earth Engine app for mapping surface water dynamics ([video](https://youtu.be/fHdwV3LEMYo) | [gif](https://i.imgur.com/GUWSVZs.gif) | [notebook](https://geemap.org/notebooks/41_water_app))
42. How to upload data to Earth Engine Apps using ipywidgets ([video](https://youtu.be/4-WeaiObj84) | [gif](https://i.imgur.com/INLzqdw.gif) | [notebook](https://geemap.org/notebooks/42_upload_data))
43. How to extract pixel values from an Earth Engine image using a point shapefile ([video](https://youtu.be/UbQ8jyc4VP4) | [gif](https://i.imgur.com/pbt6neQ.gif) | [notebook](https://geemap.org/notebooks/43_extract_values_to_points))
44. How to use Cloud Optimized GeoTIFF (COG) and SpatioTemporal Asset Catalog (STAC) ([video](https://youtu.be/yLlYoy01RxA) | [gif](https://i.imgur.com/XjG3zYq.gif) | [notebook](https://geemap.org/notebooks/44_cog_stac))
45. How to load a virtual mosaic of Cloud Optimized GeoTIFFs (COG) ([video](https://youtu.be/jDUaopr0Dhg) | [gif](https://i.imgur.com/My8Ksh7.gif) | [notebook](https://geemap.org/notebooks/45_cog_mosaic))
46. How to use locally trained machine learning models with Earth Engine ([video](https://youtu.be/nq_Ro7E0b6E) | [gif](https://i.imgur.com/muwDfkC.gif) | [notebook](https://geemap.org/notebooks/46_local_rf_training))
47. How to download image thumbnails from Earth Engine ([video](https://youtu.be/qwXZDSbfyE8) | [gif](https://i.imgur.com/gqr7CNz.gif) | [notebook](https://geemap.org/notebooks/47_image_thumbnails))
48. How to add a draggable legend to folium maps ([video](https://youtu.be/-rO1MztlLMo) | [gif](https://i.imgur.com/i2Bye9X.gif) | [notebook](https://geemap.org/notebooks/48_folium_legend))
49. How to add a colorbar to the map ([video](https://youtu.be/qiKns09X1Ao) | [gif](https://i.imgur.com/VpMq8M9.gif) | [notebook](https://geemap.org/notebooks/49_colorbar))
50. How to create publication quality maps using cartoee ([video](https://youtu.be/t24_lpYA1ko) | [gif](https://i.imgur.com/fwCzZTi.gif) | [notebook](https://geemap.org/notebooks/50_cartoee_quickstart))
51. How to create publication quality maps with custom projections ([video](https://youtu.be/3dS2EkAuAxM) | [gif](https://i.imgur.com/vvvF94j.gif) | [notebook](https://geemap.org/notebooks/51_cartoee_projections))
52. How to create timelapse animations with custom projection, scale bar, and north arrow ([video](https://youtu.be/ejuugljSut4) | [gif](https://i.imgur.com/MVQFyHN.gif) | [notebook](https://geemap.org/notebooks/52_cartoee_gif))
53. How to change layer visualization interactively with a GUI ([video](https://youtu.be/4E7gg6yaHBg) | [gif](https://i.imgur.com/VqqlMSK.gif) | [notebook](https://geemap.org/notebooks/53_layer_vis))
54. Visualizing Earth Engine vector data interactively with a GUI ([video](https://youtu.be/SIMnvbn8d-4) | [gif](https://youtu.be/F7xa5OaweY0) | [notebook](https://geemap.org/notebooks/54_vector_vis))
55. Visualizing Earth Engine raster data interactively with a GUI ([video](https://youtu.be/2R3933NFIa0) | [gif](https://youtu.be/6HFGvyXOXJM) | [notebook](https://geemap.org/notebooks/55_raster_vis))
56. Loading local vector and raster data to geemap without coding ([video](https://youtu.be/Zhwz0uS4Xi0) | [gif](https://youtu.be/SWLpnYnsqMw) | [notebook](https://geemap.org/notebooks/56_local_data))
57. Creating publication-quality maps with multiple Earth Engine layers ([video](https://youtu.be/v-FWj9dAMJ8) | [gif](https://youtu.be/85Cu3cVLmOY) | [notebook](https://geemap.org/notebooks/57_cartoee_blend))
58. Loading vector data (e.g., shp, kml, geojson) to the map without coding ([video](https://youtu.be/10KA7uhEWUM) | [gif](https://youtu.be/UsmigaIDNpE) | [notebook](https://geemap.org/notebooks/58_add_vector))
59. Using whitebox with geemap ([video](https://youtu.be/n8ODeZpuyCE) | gif | [notebook](https://geemap.org/notebooks/59_whitebox))
60. Visualizing Earth Engine data with over 200 colormaps through dot notation ([video](https://youtu.be/RBCf7wgK3Cg) | gif | [notebook](https://geemap.org/notebooks/60_colormaps))
61. Adding a scale bar to a cartoee map (video | gif | [notebook](https://geemap.org/notebooks/61_cartoee_scalebar))
62. Using the time slider for visualizing Earth Engine time-series images ([video](https://youtu.be/w_nWkNz8fyI) | gif | [notebook](https://geemap.org/notebooks/62_time_slider))
63. Creating interactive charts from Earth Engine data (video | [gif](https://youtu.be/e-GTdUUc8N8) | [notebook](https://geemap.org/notebooks/63_charts))
64. Accessing the Earth Engine Data Catalog via dot notation with autocompletion (video | [gif](https://youtu.be/hGbs2cl7otk) | [notebook](https://geemap.org/notebooks/64_data_catalog))
65. Styling Earth Engine vector data (video | gif | [notebook](https://geemap.org/notebooks/65_vector_styling))
66. Adding a legend to publication quality maps using cartoee (video | gif | [notebook](https://geemap.org/notebooks/66_cartoee_legend))
67. Creating training samples for machine learning and supervised image classification (video | [gif](https://youtu.be/VWh5PxXPZw0) | [notebook](https://geemap.org/notebooks/67_training_samples))
68. Converting NetCDF to Earth Engine Image (video | gif | [notebook](https://geemap.org/notebooks/68_netcdf_to_ee))
69. Plotting Earth Engine vector data with cartoee (video | [gif](https://youtu.be/Gr6GBuBWnnk) | [notebook](https://geemap.org/notebooks/69_cartoee_vector))
70. Creating linked maps with a few lines of code (video | [gif](https://youtu.be/AFUGje3VWM8) | [notebook](https://geemap.org/notebooks/70_linked_maps))
71. Creating Landsat timelapse animations with a few clicks (video | [gif](https://youtu.be/mA21Us_3m28) | [notebook](https://geemap.org/notebooks/71_timelapse))
72. Creating time-series cloud-free composites with a few clicks (video | [gif](https://youtu.be/kEltQkNia6o) | [notebook](https://geemap.org/notebooks/72_time_slider_gui))
73. Generating transects along lines with Earth Engine without coding (video | [gif](https://youtu.be/0TNXSs6fwrg) | [notebook](https://geemap.org/notebooks/73_transect))
74. Creating points from CSV without coding (video | gif | [notebook](https://geemap.org/notebooks/74_csv_to_points))
75. Visualizing land cover change with inteactive Sankey diagrams (video | gif | [notebook](https://geemap.org/notebooks/75_sankee))
76. Downloading and visualizing OpenStreetMap data (video | gif | [notebook](https://geemap.org/notebooks/76_osm_to_ee))
77. Adding Planet global monthly and quarterly mosaic (video | gif | [notebook](https://geemap.org/notebooks/77_planet_imagery))
78. Using timeseries inspector with one click (video | [gif](https://i.imgur.com/s1GoEOV.gif) | [notebook](https://geemap.org/notebooks/78_ts_inspector))
79. Creating histograms using the geemap chart module (video | gif | [notebook](https://geemap.org/notebooks/79_chart_histogram/))
80. Adding a point layer with popup attributes (video | gif | [notebook](https://geemap.org/notebooks/80_point_layer))
81. Creating timelapse animations from GEOS weather satellites (video | gif | [notebook](https://geemap.org/notebooks/81_goes_timelapse))
82. Creating elevation contours for any location around the globe (video | gif | [notebook](https://geemap.org/notebooks/82_contours))
83. Loading local raster datasets and Cloud Optimized GeoTIFF (COG) ([notebook](https://geemap.org/notebooks/83_local_tile))
84. Downloading OpenStreetMap data with a single line of code ([notebook](https://geemap.org/notebooks/84_openstreetmap))
85. Converting PostGIS data to ee.FeatureCollection ([notebook](https://geemap.org/notebooks/85_postgis))
86. Adding image overlay to the map ([notebook](https://geemap.org/notebooks/86_image_overlay))
87. Adding points from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://geemap.org/notebooks/87_add_points_from_xy))
88. Adding circle markers from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://geemap.org/notebooks/88_circle_markers))
89. Labeling Earth Engine FeatureCollection on the map ([notebook](https://geemap.org/notebooks/89_add_labels))
90. Creating 1-m resolution NAIP imagery timelapse ([notebook](https://geemap.org/notebooks/90_naip_timelapse))
91. Adding Planetary Computer STAC item to the map ([notebook](https://geemap.org/notebooks/91_planetary_computer))
92. Using plotly with Earth Engine ([notebook](https://geemap.org/notebooks/92_plotly))
93. Getting pixel values from COG/STAC using the Inspector tool ([notebook](https://leafmap.org/notebooks/93_cog_inspector))
94. Using heremap with Earth Engine ([notebook](https://geemap.org/notebooks/94_heremap))
95. Creating Cloud Optimized GeoTIFF (COG) ([notebook](https://geemap.org/notebooks/95_create_cog))

### 1. Introducing the geemap Python package for interactive mapping with Google Earth Engine

![Intro geemap](https://i.imgur.com/pI39k7v.gif)

### 2. Using basemaps in geemap and ipyleaflet for interactive mapping with Google Earth Engine

![basemaps](https://i.imgur.com/P5B2f7p.gif)

### 3. Introducing the Inspector tool for Earth Engine Python API

![Inspector tool](https://i.imgur.com/8d77gtI.gif)

### 4. Creating a split-panel map for visualizing Earth Engine data

![Split panel](https://i.imgur.com/kql7pC3.gif)

### 5. Using drawing tools to interact with Earth Engine data

![Drawing tools](https://i.imgur.com/Lm5pDUr.gif)

### 6. Creating an interactive map with a marker cluster

![Marker cluster](https://i.imgur.com/GF4cOqh.gif)

### 7. Converting data formats between GeoJSON and Earth Engine

![GeoJSON](https://i.imgur.com/hVPmUG1.gif)

### 8. Automated conversion from Earth Engine JavaScripts to Python scripts and Jupyter notebooks

![Conversion](https://i.imgur.com/BW0zJnN.gif)

### 9. Interactive plotting of Earth Engine data with minimal coding

![Plotting](https://i.imgur.com/iGMRnRb.gif)

### 10. Using shapefiles with Earth Engine without having to upload data to GEE

![shapefile](https://i.imgur.com/W3vNSdX.gif)

### 11. Exporting Earth Engine Image and ImageCollection as GeoTIFF and Numpy array

![exporting](https://i.imgur.com/RonLr0j.gif)

### 12. Computing zonal statistics with Earth Engine and exporting results as CSV or shapefile

![zonal](https://i.imgur.com/8xmUitW.gif)

### 13. Calculating zonal statistics by group (e.g., analyzing land cover composition of each country/state)

![zonal by group](https://i.imgur.com/LxD2em9.gif)

### 14. Adding a customized legend for Earth Engine data

![legend](https://i.imgur.com/idkZHQp.gif)

### 15. Converting Earth Engine JavaScripts to Python code directly within Jupyter notebook

![js-py](https://i.imgur.com/aGCBWSV.gif)

### 16. Adding animated text to GIF images generated from Earth Engine data

![animated text](https://i.imgur.com/MSde1om.gif)

### 17. Adding colorbar and images to GIF animations generated from Earth Engine data

![logo](https://i.imgur.com/13tFMSI.gif)

### 18. Creating Landsat timelapse animations with animated text using Earth Engine

![timelapse](https://i.imgur.com/XOHOeXk.gif)

### 19. How to search and import datasets from Earth Engine Data Catalog

![search](https://i.imgur.com/E09p64F.gif)

### 20. Using timeseries inspector to visualize landscape changes over time

![ts inspector](https://i.imgur.com/61wbRjK.gif)

### 21. Exporting Earth Engine maps as HTML files and PNG images

![export html](https://i.imgur.com/rJuXH4a.gif)

### 22. How to import Earth Engine Python scripts into Jupyter notebook?

![import scripts](https://i.imgur.com/WwJoBHF.gif)

### 23. How to search Earth Engine API and import assets from GEE personal account?

![import assets](https://i.imgur.com/b1auzkr.gif)

### 24. How to publish interactive Earth Engine maps?

![publish maps](https://i.imgur.com/Hpfzazk.gif)

### 25. How to load local raster datasets with geemap?

![load rasters](https://i.imgur.com/nsqEt2O.gif)
# geemap-heroku

Python scripts for deploying Earth Engine Apps to heroku, try it out: <https://geemap-demo.herokuapp.com/>

## How to deploy your own Earth Engine Apps?

- [Sign up](https://signup.heroku.com/) for a free heroku account.
- Follow the [instructions](https://devcenter.heroku.com/articles/getting-started-with-python#set-up) to install [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and Heroku Command Line Interface (CLI).
- Authenticate heroku using the `heroku login` command.
- Clone this repository: <https://github.com/giswqs/geemap-heroku>
- Create your own Earth Engine notebook and put it under the `notebooks` directory.
- Add Python dependencies in the `requirements.txt` file if needed.
- Edit the `Procfile` file by replacing `notebooks/geemap.ipynb` with the path to your own notebook.
- Commit changes to the repository by using `git add . && git commit -am "message"`.
- Create a heroku app: `heroku create`
- Run the `config_vars.py` script to extract Earth Engine token from your computer and set it as an environment variable on heroku: `python config_vars.py`
- Deploy your code to heroku: `git push heroku master`
- Open your heroku app: `heroku open`

## Optional steps

- To specify a name for your app, use `heroku apps:create example`
- To preview your app locally, use `heroku local web`
- To hide code cells from your app, you can edit the `Procfile` file and set `--strip_sources=True`
- To periodically check for idle kernels, you can edit the `Procfile` file and set `--MappingKernelManager.cull_interval=60 --MappingKernelManager.cull_idle_timeout=120`
- To view information about your running app, use `heroku logs --tail`
- To set an environment variable on heroku, use `heroku config:set NAME=VALUE`
- To view environment variables for your app, use `heroku config`

## Credits

The instructions above on how to deploy a voila application on heroku are adapted from [voila-dashboards/voila-heroku](https://github.com/voila-dashboards/voila-heroku).
# earthengine-py-documentation
Unofficial Google Earth Engine Python Documentation
# toolbar module

::: geemap.toolbar
# conversion module

::: geemap.conversion# Usage

Below is a list of some commonly used functions available in the **geemap** Python package. Please check the [API Reference](https://geemap.org/geemap) for a complete list of all available functions.

To create an ipyleaflet-based interactive map:

    import geemap
    Map = geemap.Map(center=[40,-100], zoom=4)
    Map

To create a folium-based interactive map:

    import geemap.foliumap as geemap
    Map = geemap.Map(center=[40,-100], zoom=4)
    Map

To add an Earth Engine data layer to the Map:

    Map.addLayer(ee_object, vis_params, name, shown, opacity)

To center the map view at a given coordinates with the given zoom level:

    Map.setCenter(lon, lat, zoom)

To center the map view around an Earth Engine object:

    Map.centerObject(ee_object, zoom)

To add LayerControl to a folium-based Map:

    Map.addLayerControl()

To add a minimap (overview) to an ipyleaflet-based Map:

    Map.add_minimap()

To add additional basemaps to the Map:

    Map.add_basemap('Esri Ocean')
    Map.add_basemap('Esri National Geographic')

To add an XYZ tile layer to the Map:

    url = 'https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}'
    Map.add_tile_layer(url, name='Google Map', attribution='Google')

To add a WMS layer to the Map:

    naip_url = 'https://services.nationalmap.gov/arcgis/services/USGSNAIPImagery/ImageServer/WMSServer?'
    Map.add_wms_layer(url=naip_url, layers='0', name='NAIP Imagery', format='image/png', shown=True)

To convert a shapefile to Earth Engine object and add it to the Map:

    ee_object = geemap.shp_to_ee(shp_file_path)
    Map.addLayer(ee_object, {}, 'Layer name')

To convert a GeoJSON file to Earth Engine object and add it to the Map:

    ee_object = geemap.geojson_to_ee(geojson_file_path)
    Map.addLayer(ee_object, {}, 'Layer name')

To download an ee.FeatureCollection as a shapefile:

    geemap.ee_to_csv(ee_object, filename, selectors)

To export an ee.FeatureCollection to other formats, including shp, csv,
json, kml, and kmz:

    geemap.ee_export_vector(ee_object, filename, selectors)

To export an ee.Image as a GeoTIFF file:

    geemap.ee_export_image(ee_object, filename, scale, crs, region, file_per_band)

To export an ee.ImageCollection as GeoTIFF files:

    geemap.ee_export_image_collection(ee_object, output, scale, crs, region, file_per_band)

To extract pixels from an ee.Image into a 3D numpy array:

    geemap.ee_to_numpy(ee_object, bands, region, properties, default_value)

To calculate zonal statistics:

    geemap.zonal_statistics(in_value_raster, in_zone_vector, out_file_path, statistics_type='MEAN')

To calculate zonal statistics by group:

    geemap.zonal_statistics_by_group(in_value_raster, in_zone_vector, out_file_path, statistics_type='SUM')

To create a split-panel Map:

    Map.split_map(left_layer='HYBRID', right_layer='ESRI')

To add a marker cluster to the Map:

    Map.marker_cluster()
    feature_collection = ee.FeatureCollection(Map.ee_markers)

To add a customized legend to the Map:

    legend_dict = {
        'one': (0, 0, 0),
        'two': (255,255,0),
        'three': (127, 0, 127)
    }
    Map.add_legend(legend_title='Legend', legend_dict=legend_dict, position='bottomright')
    Map.add_legend(builtin_legend='NLCD')

To download a GIF from an Earth Engine ImageCollection:

    geemap.download_ee_video(tempCol, videoArgs, saved_gif)

To add animated text to an existing GIF image:

    geemap.add_text_to_gif(in_gif, out_gif, xy=('5%', '5%'), text_sequence=1984, font_size=30, font_color='#0000ff', duration=100)

To create a colorbar for an Earth Engine image:

    palette = ['blue', 'purple', 'cyan', 'green', 'yellow', 'red']
    create_colorbar(width=250, height=30, palette=palette, vertical=False,add_labels=True, font_size=20, labels=[-40, 35])

To create a Landsat timelapse animation and add it to the Map:

    Map.add_landsat_ts_gif(label='Place name', start_year=1985, bands=['NIR', 'Red', 'Green'], frames_per_second=5)

To convert all GEE JavaScripts in a folder recursively to Python
scripts:

    from geemap.conversion import *
    js_to_python_dir(in_dir, out_dir)

To convert all GEE Python scripts in a folder recursively to Jupyter
notebooks:

    from geemap.conversion import *
    template_file = get_nb_template()
    py_to_ipynb_dir(in_dir, template_file, out_dir)

To execute all Jupyter notebooks in a folder recursively and save output
cells:

    from geemap.conversion import *
    execute_notebook_dir(in_dir)

To search Earth Engine API documentation with Jupyter notebooks:

    import geemap
    geemap.ee_search()

To publish an interactive GEE map with Jupyter notebooks:

    Map.publish(name, headline, visibility)

To add a local raster dataset to the map:

    Map.add_raster(image, bands, colormap, layer_name)

To get image basic properties:

    geemap.image_props(image).getInfo()

To get image descriptive statistics:

    geemap.image_stats(image, region, scale)

To remove all user-drawn geometries:

    geemap.remove_drawn_features()

To extract pixel values based on user-drawn geometries:

    geemap.extract_values_to_points(out_shp)
# plotlymap module

::: geemap.plotlymap# Changelog

## v0.11.1 - January 24, 2022

**New Features**:

-   Added ee_extra to algorithms [#868](https://github.com/giswqs/geemap/pull/868)
-   Added COG creation [bc83fdf](https://github.com/giswqs/geemap/commit/bc83fdf959dd91bdeb40de6af41fea79933c57c2)
-   Added heremap plotting backend [#382](https://github.com/giswqs/geemap/issues/382)
-   Added COG Inspector GUI [#841](https://github.com/giswqs/geemap/issues/841)

**Improvement**:

-   Improved GitHub workflows [#879](https://github.com/giswqs/geemap/pull/879)
-   Fixed ee_stac_list bug [#873](https://github.com/giswqs/geemap/issues/873)
-   Fixed js py conversion [ce7fee0](https://github.com/giswqs/geemap/commit/ce7fee0e87ab2a35069e45d141fb2c84333f392b)
-   Updated notebook 07 [#871](https://github.com/giswqs/geemap/pull/871)
-   Added IR band to goes_timelapse [#870](https://github.com/giswqs/geemap/pull/870)
-   Updated ee_basemaps [#869](https://github.com/giswqs/geemap/pull/869)
-   Removed COG mosaic
-   Fixed cartoee legend bug
-   Updated installation instructions

## v0.11.0 - January 7, 2022

**New Features**:

-   Added support for plotly [#842](https://github.com/giswqs/geemap/issues/842)
-   Added colorbar to timelapse [#846](https://github.com/giswqs/geemap/issues/846)
-   Added save_colorbar function [#846](https://github.com/giswqs/geemap/issues/846)
-   Added ocean color timelapse [#845](https://github.com/giswqs/geemap/issues/845)
-   Added support for xyzservices basemaps [#795](https://github.com/giswqs/geemap/issues/795)
-   Added labeling gdf shp geojson [#815](https://github.com/giswqs/geemap/issues/815)
-   Added remove_labels [#815](https://github.com/giswqs/geemap/issues/815)
-   Added Planetary Computer STAC support
-   Added bbox_to_gdf function

**Improvement**:

-   Fixed cartoee projection bug [#843](https://github.com/giswqs/geemap/discussions/843)
-   Improved COG visualization [#844](https://github.com/giswqs/geemap/issues/844)
-   Updated STAC notebook example [#841](https://github.com/giswqs/geemap/issues/841)
-   Improved stac tile functions [#839](https://github.com/giswqs/geemap/pull/839)
-   Removed pangeo broken binder links

## v0.10.2 - December 23, 2021

**New Features**:

-   Add locate control to folium [#809](https://github.com/giswqs/geemap/issues/809)
-   Added add_points_from_xy function [#812](https://github.com/giswqs/geemap/issues/812)
-   Added heatmap function
-   Added add_labels function [#815](https://github.com/giswqs/geemap/issues/815)
-   Added NAIP timelapse [#789](https://github.com/giswqs/geemap/issues/789)

**Improvement**:

-   Improved js_to_py function [#805](https://github.com/giswqs/geemap/discussions/805)
-   Renamed popups to popup [#812](https://github.com/giswqs/geemap/issues/812)
-   Changed default map view [#821](https://github.com/giswqs/geemap/issues/821)
-   Fixed centerObject bug [#823](https://github.com/giswqs/geemap/issues/823)
-   Fixed typo [#824](https://github.com/giswqs/geemap/pull/824)

## v0.10.1 - December 6, 2021

**Improvement**:

-   A temporary fix for ipyleaflet basemap error [#795](https://github.com/giswqs/geemap/issues/795)

## v0.10.0 - November 28, 2021

**New Features**:

-   Added remove_legend function [#761](https://github.com/giswqs/geemap/issues/761)
-   Added add_marker function [#765](https://github.com/giswqs/geemap/pull/765)
-   Added support for local tile and raster GUI [#758](https://github.com/giswqs/geemap/issues/758), [#769](https://github.com/giswqs/geemap/pull/769)
-   Added a new osm module [#770](https://github.com/giswqs/geemap/issues/770) [#772](https://github.com/giswqs/geemap/pull/772)
-   Added support for PostGIS [#771](https://github.com/giswqs/geemap/issues/771) [#772](https://github.com/giswqs/geemap/pull/772)
-   Added ImageOverlay from local files [#773](https://github.com/giswqs/geemap/issues/773)

## v0.9.5 - November 22, 2021

**New Features**:

-   Added timelapse module
-   Added quarter and monthly timelapse [#746](https://github.com/giswqs/geemap/issues/746)
-   Improved create timeseries [#736](https://github.com/giswqs/geemap/issues/736)
-   Added Sentinel-2 timelapse [#733](https://github.com/giswqs/geemap/issues/733) [#736](https://github.com/giswqs/geemap/issues/736)
-   Added MODIS NDVI timelapse [#728](https://github.com/giswqs/geemap/issues/728)
-   Added GOES timelapse [#717](https://github.com/giswqs/geemap/issues/717)
-   Added time slider opacity param [#720](https://github.com/giswqs/geemap/discussions/720)
-   Added contour function [#688](https://github.com/giswqs/geemap/issues/688)
-   Added more gif functions
-   Added make_gif and gif_to_mp4 functions
-   Improved date sequence
-   Added Alibaba font type
-   Added ESA Land Cover legend
-   Added zoom to bounds function
-   Added streamlit download button

**Improvement**:

-   Fixed encoding bug [#747](https://github.com/giswqs/geemap/issues/747)

## v0.9.4 - October 23, 2021

**New Features**:

-   Made streamlit map width responsive [#713](https://github.com/giswqs/geemap/issues/713)
-   Added function read file from url

**Improvement**:

-   Fixed map width bug [#712](https://github.com/giswqs/geemap/issues/712)
-   Fixed algorithms module bug
-   Updated environment.yml

## v0.9.3 - October 23, 2021

**New Features**:

-   Added streamlit support [#697](https://github.com/giswqs/geemap/issues/697)
-   Added point layer function [#702](https://github.com/giswqs/geemap/issues/702)
-   Added river width module [#682](https://github.com/giswqs/geemap/issues/682)
-   Added census data and xyzservices
-   Added nlcd notebook
-   Added river width module notebook
-   Added GEE workshop notebook

**Improvement**:

-   Fixed geojson style callback bug [#692](https://github.com/giswqs/geemap/issues/692)
-   Fixed open vector bug [#124](https://github.com/giswqs/geemap/issues/124)
-   Removed py36 due to xyzservices

## v0.9.2 - October 1, 2021

**New Features**:

-   Added RivWidthCloud module [#682](https://github.com/giswqs/geemap/issues/682)
-   Added RivWidthCloud notebook [#682](https://github.com/giswqs/geemap/issues/682)
-   Added [NLCD notebook](https://geemap.org/notebooks/nlcd_app/)
-   Added a close button to timeseries inspector

**Improvement**:

-   Fixed hover countries notebook [#686](https://github.com/giswqs/geemap/pull/686)
-   Improved cartoee colorbar with custom label size [#681](https://github.com/giswqs/geemap/discussions/681)

## v0.9.1 - September 17, 2021

**New Features**:

-   Added `sandbox_path` option allowing users to restrict Voila app access to system directories [#673](https://github.com/giswqs/geemap/issues/673)

## v0.9.0 - September 10, 2021

**New Features**:

-   Get current device latlon [#618](https://github.com/giswqs/geemap/issues/618)

**Improvement**:

-   Improved Colab support [#661](https://github.com/giswqs/geemap/issues/661)
-   Improved folium colorbar [#586](https://github.com/giswqs/geemap/issues/586)
-   Fixed broken link [#653](https://github.com/giswqs/geemap/issues/653)
-   Fixed extract pixel values bug [#610](https://github.com/giswqs/geemap/issues/610)
-   Fixed color palette bug [#605](https://github.com/giswqs/geemap/pull/605)
-   Fixed typos [#589](https://github.com/giswqs/geemap/pull/589)

## v0.8.18 - July 8, 2021

**New Features**:

-   Added pandas_to_geojson [#557](https://github.com/giswqs/geemap/discussions/557)
-   Added feature_histogram function to chart module [#553](https://github.com/giswqs/geemap/pull/553)
-   Added feature_groups function to chart module [#539](https://github.com/giswqs/geemap/pull/539)
-   Added random forest probability output [#550](https://github.com/giswqs/geemap/pull/550)

**Improvement**:

-   Renamed eefolium module to foliumap
-   Changed COG and STAC to lowercase
-   Changed .format() to fstring [#561](https://github.com/giswqs/geemap/pull/561)
-   Fixed random forest string to label bug [#545](https://github.com/giswqs/geemap/pull/545)
-   Improved split-panel map [#543](https://github.com/giswqs/geemap/discussions/543)
-   Updated otsu example [#535](https://github.com/giswqs/geemap/discussions/535)

## v0.8.17 - June 20, 2021

**New Features**:

-   Added Planet global mosaic [#527](https://github.com/giswqs/geemap/issues/527)
-   Add LCMS dataset option for sankee [#517](https://github.com/giswqs/geemap/issues/517)
-   Added add_osm function [#503](https://github.com/giswqs/geemap/discussions/503)

**Improvement**:

-   Added otsu example [#535](https://github.com/giswqs/geemap/discussions/535)
-   Fixed timeseries plotting bug [#513](https://github.com/giswqs/geemap/discussions/513)
-   Fixed shp deletion bug [#509](https://github.com/giswqs/geemap/discussions/509)
-   Fixed csv_to_points bug [#490](https://github.com/giswqs/geemap/discussions/490)
-   Improved ee_to_geojson [#486](https://github.com/giswqs/geemap/pull/486)
-   Improved random sampling notebook [#479](https://github.com/giswqs/geemap/discussions/479)
-   Fixed link bug [#480](https://github.com/giswqs/geemap/issues/480)
-   Improved sankee notebook [#471](https://github.com/giswqs/geemap/issues/471)
-   Updated installation docs
-   Added binder env

## v0.8.16 - May 10, 2021

**New Features**:

-   Added csv_to_points GUI [#461](https://github.com/giswqs/geemap/issues/461)
-   Added GUI for creating transects [#454](https://github.com/giswqs/geemap/issues/454)
-   Added csv_to_ee and csv_to_makers [#461](https://github.com/giswqs/geemap/issues/461)
-   Added geopandas support [#455](https://github.com/giswqs/geemap/issues/455)

**Improvement**:

-   Improved geojson style [#459](https://github.com/giswqs/geemap/issues/459) [#460](https://github.com/giswqs/geemap/issues/460)
-   Improved vector support [#455](https://github.com/giswqs/geemap/issues/455)
-   Improved add_colorbar function [#450](https://github.com/giswqs/geemap/issues/450)
-   Improved add_raster function [#449](https://github.com/giswqs/geemap/pull/449)
-   Updated notebooks

## v0.8.15 - April 28, 2021

**Improvement**:

-   Improved shp_to_geojson function [#430](https://github.com/giswqs/geemap/discussions/430)
-   Improved add_styled_vector function [#432](https://github.com/giswqs/geemap/discussions/432)
-   Fixed map publish bug [#445](https://github.com/giswqs/geemap/issues/445)
-   Improved add_colorbar function [dc7e548](https://github.com/giswqs/geemap/commit/dc7e54856694a1994b6d4f4044385babe04bd086)

## v0.8.14 - April 20, 2021

**New Features**:

-   Added timelapse GUI [#359](https://github.com/giswqs/geemap/issues/359)
-   Added timeslider GUI [#359](https://github.com/giswqs/geemap/issues/359) [#387](https://github.com/giswqs/geemap/issues/387)

**Improvement**:

-   Improved add_geojson function [731e59e](https://github.com/giswqs/geemap/commit/731e59efc4a1f629db13f6b6cc4e9ef6b06cbe8f)
-   Added GeoPython workshop notebook [6efd5e](https://geemap.org/workshops/GeoPython_2021)
-   Improved cartoee colorbar [#413](https://github.com/giswqs/geemap/discussions/413)
-   Improved cartoee add_layer function [#368](https://github.com/giswqs/geemap/issues/368)

## v0.8.13 - March 22, 2021

**New Features**:

-   Added linked maps [#375](https://github.com/giswqs/geemap/issues/375)
-   Added cartoee legend [#343](https://github.com/giswqs/geemap/issues/343)
-   Added chart by feature property [#339](https://github.com/giswqs/geemap/issues/339)
-   Added tool gui template [#239](https://github.com/giswqs/geemap/issues/239)
-   Added GEE Toolbox GUI [#362](https://github.com/giswqs/geemap/issues/362)
-   Added support for multiple legends [#365](https://github.com/giswqs/geemap/discussions/365)

**Improvement**:

-   Improved dataset module to use GEE STAC [#346](https://github.com/giswqs/geemap/issues/346)
-   Improved training sample tool [#326](https://github.com/giswqs/geemap/issues/326)
-   Added netcdf_to_ee example [#285](https://github.com/giswqs/geemap/issues/285)
-   Improved to_html function [#361](https://github.com/giswqs/geemap/discussions/361)
-   Changed colorbar plotting backend [#372](https://github.com/giswqs/geemap/issues/372)
-   Improved get_colorbar function [#372](https://github.com/giswqs/geemap/issues/372)
-   Added vector styling example
-   Improved zonal statistics

## v0.8.12 - March 8, 2021

**New Features**:

-   Added a dataset module for accessing the Earth Engine Data Catalog via dot notation [#345](https://github.com/giswqs/geemap/issues/345)
-   Added a chart module for creating interactive charts for Earth Engine data [#343](https://github.com/giswqs/geemap/issues/343)
-   Added a time slider for visualizing Earth Engine time-series images [#335 ](https://github.com/giswqs/geemap/issues/335) [#344](https://github.com/giswqs/geemap/issues/344)
-   Added a `netcdf_to_ee` function [#342](https://github.com/giswqs/geemap/pull/342)
-   Added a `numpy_to_ee` function [#337](https://github.com/giswqs/geemap/pull/337)
-   Added vertical colorbar support [#322](https://github.com/giswqs/geemap/issues/322)
-   Added GUI for creating training samples [#326](https://github.com/giswqs/geemap/issues/326)

**Improvement**:

-   Added layer control by default to folium map [#323](https://github.com/giswqs/geemap/issues/323)
-   Added geemap matplotlib example [#319](https://github.com/giswqs/geemap/discussions/319)
-   Added lgtm continuous integration [#330](https://github.com/giswqs/geemap/issues/330)
-   Fixed layer palette bug [#334](https://github.com/giswqs/geemap/issues/334)
-   Fixed minimap zoom parameter [#329](https://github.com/giswqs/geemap/pull/329)
-   Fixed centerObject bug

## v0.8.11 - February 23, 2021

**New Features**:

-   Added a colormap module [#302](https://github.com/giswqs/geemap/issues/302)
-   Added a new cartoee scale bar function [#313](https://github.com/giswqs/geemap/pull/313)
-   Added extract pixel values function [#315](https://github.com/giswqs/geemap/issues/315)
-   Visualizing Earth Engine image with >200 matplotlib colormaps via dot notation ([example](https://geemap.org/notebooks/60_colormaps/))

**Improvement**:

-   Improved the basemap module accessible via dot notation [#302](https://github.com/giswqs/geemap/issues/302)
-   Added googledrivedownloader and python-box to requirements [#310](https://github.com/giswqs/geemap/discussions/310)
-   Fixed folium layer name bug [#314](https://github.com/giswqs/geemap/issues/314)

## v0.8.10 - February 16, 2021

**New Features**:

-   Added default basemap options when creating the Map [#293](https://github.com/giswqs/geemap/issues/293)
-   Added GUI for change basemaps [#294](https://github.com/giswqs/geemap/issues/294)
-   Added GUI for js2py conversion [#296](https://github.com/giswqs/geemap/issues/296)
-   Added geemap cheat sheet [#276](https://github.com/giswqs/geemap/issues/276)
-   Added `Map.zoomToObject()` method [#303](https://github.com/giswqs/geemap/issues/303)

**Improvement**:

-   Improved `Map.centerObject()` method [#303](https://github.com/giswqs/geemap/issues/303)

## v0.8.9 - February 4, 2021

**New Features**:

-   Added [whiteboxgui](https://github.com/giswqs/whiteboxgui) with 477 geoprocessing tools [#254](https://github.com/giswqs/geemap/issues/254)

**Improvement**:

-   Fixed file open encoding bug

## v0.8.8 - January 17, 2021

**New Features**:

-   Added support for converting Pandas/GeoPandas DataFrame to ee.FeatureCollection and vice versa [#268](https://github.com/giswqs/geemap/issues/268)
-   Added KML/KMZ support [#247](https://github.com/giswqs/geemap/issues/247)
-   Added Code of Conduct

**Improvement**:

-   Fixed CSV encoding bug [#267](https://github.com/giswqs/geemap/issues/267)
-   Improved downloading shp support [#263](https://github.com/giswqs/geemap/issues/263)
-   Fixed WMS bug [#250](https://github.com/giswqs/geemap/discussions/250)
-   Added cartoee subplots example [#238](https://github.com/giswqs/geemap/discussions/238)
-   Reformatted code using black formatter
-   Improved support for shp and geojson [#244](https://github.com/giswqs/geemap/issues/244)
-   Fixed layer control bug
-   Added cartoee blend tutorial [#241](https://github.com/giswqs/geemap/issues/241)
-   Improved drawing tools [#240](https://github.com/giswqs/geemap/issues/240)
-   Improved Inspector tool

## v0.8.7 - December 27, 2020

**New Features**:

-   Added toolbar GUI [#215](https://github.com/giswqs/geemap/issues/215)
-   Added layer vis [#215](https://github.com/giswqs/geemap/issues/215)
-   Added raster/vector colormap [#215](https://github.com/giswqs/geemap/issues/215)
-   Added support for linking legend with layer [#234](https://github.com/giswqs/geemap/issues/234)
-   Added styled vector function [#235](https://github.com/giswqs/geemap/issues/235)
-   Added mouse click observe to toolbar [#215](https://github.com/giswqs/geemap/issues/215)
-   Added new tool for opening local data [#239](https://github.com/giswqs/geemap/issues/239)

**Improvement**:

-   Fixed COG mosaic bug [#236](https://github.com/giswqs/geemap/issues/236) and [#237](https://github.com/giswqs/geemap/issues/237)

## v0.8.6 - December 22, 2020

**New Features**:

-   Added GUI for changing layer visualization interactively [#215](https://github.com/giswqs/geemap/issues/215)
-   Added a toolbar [#215](https://github.com/giswqs/geemap/issues/215)
-   Added color bar support [#223](https://github.com/giswqs/geemap/issues/223)
-   Added draggable legend to folium maps [#224](https://github.com/giswqs/geemap/issues/224)
-   Added `get_image_collection_gif()` function [#225](https://github.com/giswqs/geemap/issues/225)
-   Added `image_dates()` function [#216](https://github.com/giswqs/geemap/issues/216)

**Improvement**:

-   Added `max_zoom` parameter to `add_tile_layer()` [#227](https://github.com/giswqs/geemap/issues/227)
-   Added mouse latlon to insepctor tool [#229](https://github.com/giswqs/geemap/discussions/229)
-   Added download icon to notebooks [#202](https://github.com/giswqs/geemap/issues/202)
-   Added GitHub issue template [#202](https://github.com/giswqs/geemap/issues/202)
-   Added more tutorials (cartoee gif, legend, color bar, vis GUI, etc.)
-   Fixed remove control bug [#218](https://github.com/giswqs/geemap/discussions/218)
-   Fixed split-panel map bug
-   Improved Exception handling

## v0.8.5 - December 12, 2020

**New Features**:

-   Add toolbar [#6](https://github.com/giswqs/geemap/issues/6)
-   Add functions for downloading imgae thumbnails [#214](https://github.com/giswqs/geemap/issues/214)
-   Add func for getting image collection dates [#216](https://github.com/giswqs/geemap/issues/216)
-   Add cartoee scale bar and north arrow [#191](https://github.com/giswqs/geemap/issues/191)
-   Add support for COG mosaic [#200](https://github.com/giswqs/geemap/issues/200)

**Improvement**:

-   Improve support for locally trained models [#210](https://github.com/giswqs/geemap/issues/210)
-   Add verbose option of downloading functions [#197](https://github.com/giswqs/geemap/pull/197)
-   Improve Inspector tool for point geometry [#198](https://github.com/giswqs/geemap/issues/198)
-   Add tutorials (COG, STAC, local RF, image thumbnails)

## v0.8.4 - December 6, 2020

**New Features:**

-   Add support for Cloud Optimized GeoTIFF (COG) and SpatioTemporal Asset Catalog (STAC) #192
-   Add [Map.add_cog_layer()](https://geemap.org/geemap/#geemap.geemap.Map.add_cog_layer) and [Map.add_stac_layer()](https://geemap.org/geemap/#geemap.geemap.Map.add_stac_layer)
-   Add new COG functions, e.g., `cog_tile()`, `cog_bounds()`, `cog_center()`, `cog_bands()`
-   Add new STAC functions, e.g., `stac_tile()`, `stac_bounds()`, `stac_center()`, `stac_bands()`

**Improvements:**

-   Improve Google Colab support #193. Use `import geemap` rather than `import geemap.foliumap as geemap`
-   Add `Open in Colab` button to notebooks #194

## v0.8.3 - December 2, 2020

**New Features:**

-   Add button for removing user-drawn features #182
-   Add function for moving drawn layer to top
-   Add remove_last_drawn() function #130
-   Add support for QGIS Layer Style File #174
-   Add mouse click get coordinates example #173
-   Add cartoee colab example #157
-   Add notebooks to mkdocs

**Improvements:**

-   Improve ee_Initialize() #189 #190
-   Fix cartoee map orientation bug #177 #183
-   Fix problematic Date field in shapefile #176
-   Fix Windows unzip bug

## v0.8.2 - November 6, 2020

**Improvements**

-   Reorganize modules
-   Add a new module common.py
-   Add new domain geemap.org
-   Format code using black
-   Add more init options for Map class

## v0.8.1 - October 27, 2020

**New Features:**

-   Add machine learning module #124 #156
-   Add cartoee module #157 #161
-   Add more tutorials (e.g., timelapse, water app, ipywidgets)

**Improvements:**

-   Make ee_Initialize() optional for Map class

BIG THANK YOU to [Kel Markert](https://github.com/kmarkert) for adding the cartoee and ml modules!!

## v0.8.0 - October 10, 2020

**Improvements**

-   Add support for loading Cloud Optimized GeoTIFFs as ee.Image and ee.ImageCollection
-   Make fmask optional when creating Landsat timelapse
-   Add support for creating timelapse of spectral indices (e.g., NDWI, NDVI)
-   Add geemap Colab tutorial
-   Add timelapse download option for voila
-   Add pydeck tutorial for visualizing 3D terrain data
-   Add qualityMosaic() tutorial

**Fixes**

-   Fix Windows zipfile bug

## v0.7.13 - September 15, 2020

**Improvements**

-   Improve ee authentication in Colab #145
-   Improve non-interactive mode #138
-   Add Colab notebook example

**Fixes**

-   Fix automated testing error
-   Fix Windows ee_search() bug

## v0.7.12 - September 1, 2020

**Improvements**

-   Rebuild docs using mkdocs-material
-   Add Internet proxy function
-   Add support for exporting shp and geojson #63

**Fixes**

-   Fix heroko config bug
-   Fix landsat timelapse bug #99 #134
-   Fix js_py conversion bug #136

## v0.7.11 - August 16, 2020

**Improvements:**

-   Add function for removing drawn features #130
-   Add function for extracting pixel values #131
-   Add function for interactive region reduction #35
-   Add machine learning tutorials

**Fixes:**

-   [Fix js_py conversion bug](https://github.com/giswqs/geemap/commit/6c0ebe4006d60f9ebb4390d0914400fc276e2c7d)
-   Fix typos

## v0.7.10 - August 5, 2020

**Improvements:**

-   Add function for getting image properties
-   Add function for calculating descriptive statistics (i.e., min, max, mean, std, sum)
-   Add more utils functions

## v0.7.7 - August 5, 2020

**Improvements:**

-   Add support for publishing maps #109
-   Add `find_layer()` function
-   Add `layer_opacity()` function
-   Update Readthedocs

**Fixes:**

-   Fix duplicate layer bug

## v0.7.0 - May 22, 2020

## v0.6.0 - April 5, 2020

## v0.5.0 - March 23, 2020

## v0.4.0 - March 19, 2020

## v0.3.0 - March 18, 2020

## v0.2.0 - March 17, 2020

## v0.1.0 - March 8, 2020
# Welcome to geemap

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/geemap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/geemap-binder)
[![image](https://img.shields.io/pypi/v/geemap.svg)](https://pypi.python.org/pypi/geemap)
[![image](https://img.shields.io/conda/vn/conda-forge/geemap.svg)](https://anaconda.org/conda-forge/geemap)
[![image](https://pepy.tech/badge/geemap)](https://pepy.tech/project/geemap)
[![image](https://github.com/giswqs/geemap/workflows/docs/badge.svg)](https://geemap.org)
[![image](https://github.com/giswqs/geemap/workflows/build/badge.svg)](https://github.com/giswqs/geemap/actions?query=workflow%3Abuild)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/giswqs/geemap.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/giswqs/geemap/context:python)
[![image](https://img.shields.io/badge/YouTube-Channel-red)](https://www.youtube.com/c/QiushengWu)
[![image](https://img.shields.io/twitter/follow/giswqs?style=social)](https://twitter.com/giswqs)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![image](https://joss.theoj.org/papers/10.21105/joss.02305/status.svg)](https://joss.theoj.org/papers/10.21105/joss.02305)

**A Python package for interactive mapping with Google Earth Engine, ipyleaflet, and ipywidgets.**

-   GitHub repo: <https://github.com/giswqs/geemap>
-   Documentation: <https://geemap.org>
-   PyPI: <https://pypi.org/project/geemap>
-   Conda-forge: <https://anaconda.org/conda-forge/geemap>
-   360+ GEE notebook examples: <https://github.com/giswqs/earthengine-py-notebooks>
-   GEE Tutorials on YouTube: <https://www.youtube.com/c/QiushengWu>
-   Free software: [MIT license](https://opensource.org/licenses/MIT)

## Introduction

**geemap** is a Python package for interactive mapping with [Google Earth Engine](https://earthengine.google.com/) (GEE), which is a cloud computing platform with a [multi-petabyte catalog](https://developers.google.com/earth-engine/datasets/) of satellite imagery and geospatial datasets. During the past few years, GEE has become very popular in the geospatial community and it has empowered numerous environmental applications at local, regional, and global scales. GEE provides both JavaScript and Python APIs for making computational requests to the Earth Engine servers. Compared with the comprehensive [documentation](https://developers.google.com/earth-engine) and interactive IDE (i.e., [GEE JavaScript Code Editor](https://code.earthengine.google.com/)) of the GEE JavaScript API, the GEE Python API has relatively little documentation and limited functionality for visualizing results interactively. The geemap Python package was created to fill this gap. It is built upon [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) and [ipywidgets](https://github.com/jupyter-widgets/ipywidgets), and enables users to analyze and visualize Earth Engine datasets interactively within a Jupyter-based environment.

**geemap** is intended for students and researchers, who would like to utilize the Python ecosystem of diverse libraries and tools to explore Google Earth Engine. It is also designed for existing GEE users who would like to transition from the GEE JavaScript API to Python API. The automated JavaScript-to-Python [conversion module](https://github.com/giswqs/geemap/blob/master/geemap/conversion.py) of the geemap package can greatly reduce the time needed to convert existing GEE JavaScripts to Python scripts and Jupyter notebooks.

For video tutorials and notebook examples, please visit the [examples page](https://github.com/giswqs/geemap/tree/master/examples). For complete documentation on geemap modules and methods, please visit the [API Reference](https://geemap.org/geemap/).

If you find geemap useful in your research, please consider citing the following papers to support my work. Thank you for your support.

-   Wu, Q., (2020). geemap: A Python package for interactive mapping with Google Earth Engine. The Journal of Open Source Software, 5(51), 2305. <https://doi.org/10.21105/joss.02305>
-   Wu, Q., Lane, C. R., Li, X., Zhao, K., Zhou, Y., Clinton, N., DeVries, B., Golden, H. E., & Lang, M. W. (2019). Integrating LiDAR data and multi-temporal aerial imagery to map wetland inundation dynamics using Google Earth Engine. Remote Sensing of Environment, 228, 1-13. <https://doi.org/10.1016/j.rse.2019.04.015> ([pdf](https://gishub.org/2019_rse) | [source code](https://doi.org/10.6084/m9.figshare.8864921))

Check out the geemap workshop I presented at the GeoPython Conference 2021. This workshop gives a comprehensive introduction to the key features of geemap.

[![geemap workship](https://img.youtube.com/vi/wGjpjh9IQ5I/0.jpg)](https://www.youtube.com/watch?v=wGjpjh9IQ5I)

## Key Features

Below is a partial list of features available for the geemap package. Please check the [examples](https://github.com/giswqs/geemap/tree/master/examples) page for notebook examples, GIF animations, and video tutorials.

-   Convert Earth Engine JavaScripts to Python scripts and Jupyter notebooks.
-   Display Earth Engine data layers for interactive mapping.
-   Support Earth Engine JavaScript API-styled functions in Python, such as `Map.addLayer()`, `Map.setCenter()`, `Map.centerObject()`, `Map.setOptions()`.
-   Create split-panel maps with Earth Engine data.
-   Retrieve Earth Engine data interactively using the Inspector Tool.
-   Interactive plotting of Earth Engine data by simply clicking on the map.
-   Convert data format between GeoJSON and Earth Engine.
-   Use drawing tools to interact with Earth Engine data.
-   Use shapefiles with Earth Engine without having to upload data to one's GEE account.
-   Export Earth Engine FeatureCollection to other formats (i.e., shp, csv, json, kml, kmz).
-   Export Earth Engine Image and ImageCollection as GeoTIFF.
-   Extract pixels from an Earth Engine Image into a 3D numpy array.
-   Calculate zonal statistics by group.
-   Add a customized legend for Earth Engine data.
-   Convert Earth Engine JavaScripts to Python code directly within Jupyter notebook.
-   Add animated text to GIF images generated from Earth Engine data.
-   Add colorbar and images to GIF animations generated from Earth Engine data.
-   Create Landsat timelapse animations with animated text using Earth Engine.
-   Search places and datasets from Earth Engine Data Catalog.
-   Use timeseries inspector to visualize landscape changes over time.
-   Export Earth Engine maps as HTML files and PNG images.
-   Search Earth Engine API documentation within Jupyter notebooks.
-   Import Earth Engine assets from personal account.
-   Publish interactive GEE maps directly within Jupyter notebook.
-   Add local raster datasets (e.g., GeoTIFF) to the map.
-   Perform image classification and accuracy assessment.
-   Extract pixel values interactively and export as shapefile and csv.

## YouTube Channel

I have created a [YouTube Channel](https://www.youtube.com/c/QiushengWu) for sharing **geemap** tutorials. You can subscribe to my channel for regular updates. If there is any specific tutorial you would like to see, please submit a feature request [here](https://github.com/giswqs/geemap/issues).

[![Earth Engine Tutorials on YouTube](https://wetlands.io/file/images/youtube.png)](https://www.youtube.com/c/QiushengWu)
# Installation

## Earth Engine Account

To use **geemap**, you must first [sign up](https://earthengine.google.com/signup/) for a [Google Earth Engine](https://earthengine.google.com/) account.
You cannot use Google Earth Engine unless your application has been approved. Once you receive the application approval email, you can log in to
the [Earth Engine Code Editor](https://code.earthengine.google.com/) to get familiar with the JavaScript API.

![signup](https://i.imgur.com/ng0FzUT.png)

## Install from PyPI

**geemap** is available on [PyPI](https://pypi.org/project/geemap/). To install **geemap**, run this command in your terminal:

```bash
    pip install geemap
```

## Install from conda-forge

**geemap** is also available on [conda-forge](https://anaconda.org/conda-forge/geemap). If you have
[Anaconda](https://www.anaconda.com/distribution/#download-section) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your computer, you can install geemap using the following command:

```bash
    conda install geemap -c conda-forge
```

The geemap package has an optional dependency - [geopandas](https://geopandas.org/), which can be challenging to install on some computers, especially Windows. It is highly recommended that you create a fresh conda environment to install geopandas and geemap. Follow the commands below to set up a conda env and install geopandas, xarray_leaflet, and geemap.

```bash
    conda create -n gee python=3.9
    conda activate gee
    conda install geopandas
    conda install mamba -c conda-forge
    mamba install geemap localtileserver -c conda-forge
```

Optionally, you can install some [Jupyter notebook extensions](https://github.com/ipython-contrib/jupyter_contrib_nbextensions), which can improve your productivity in the notebook environment. Some useful extensions include Table of Contents, Gist-it, Autopep8, Variable Inspector, etc. See this [post](https://towardsdatascience.com/jupyter-notebook-extensions-517fa69d2231) for more information.

```bash
    conda install jupyter_contrib_nbextensions -c conda-forge
```

Check the **YouTube** video below on how to install geemap using conda.

[![geemap](http://img.youtube.com/vi/h0pz3S6Tvx0/0.jpg)](http://www.youtube.com/watch?v=h0pz3S6Tvx0 "Install geemap")

## Install from GitHub

To install the development version from GitHub using [Git](https://git-scm.com/), run the following command in your terminal:

```bash
    pip install git+https://github.com/giswqs/geemap
```

## Upgrade geemap

If you have installed **geemap** before and want to upgrade to the latest version, you can run the following command in your terminal:

```bash
    pip install -U geemap
```

If you use conda, you can update geemap to the latest version by running the following command in your terminal:

```bash
    mamba update -c conda-forge geemap
```

To install the development version from GitHub directly within Jupyter notebook without using Git, run the following code:

```python
    import geemap
    geemap.update_package()
```

## Use Docker

To use geemap in a Docker container, check out the following docker containers with geemap installed.

-   [gee-community/ee-jupyter-contrib](https://github.com/gee-community/ee-jupyter-contrib/tree/master/docker/gcp_ai_deep_learning_platform)
-   [bkavlak/geemap](https://hub.docker.com/r/bkavlak/geemap)
-   [giswqs/geemap](https://hub.docker.com/r/giswqs/geemap)
# common module

::: geemap.common# basemaps module

::: geemap.basemaps# timelapse module

::: geemap.timelapse
# chart module

::: geemap.chart
# ml module

::: geemap.ml# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given. You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at <https://github.com/giswqs/geemap/issues>.

If you are reporting a bug, please include:

- Your operating system name and version.
- Any details about your local setup that might be helpful in troubleshooting.
- Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

### Write Documentation

geemap could always use more documentation, whether as part of the official geemap docs, in docstrings, or even on the web in blog posts, articles, and such.

### Submit Feedback

The best way to send feedback is to file an issue at <https://github.com/giswqs/geemap/issues>.

If you are proposing a feature:

- Explain in detail how it would work.
- Keep the scope as narrow as possible, to make it easier to implement.
- Remember that this is a volunteer-driven project, and that contributions are welcome :)

## Get Started

Ready to contribute? Here's how to set up _geemap_ for local development.

1. Fork the [geemap](https://github.com/giswqs/geemap) repo on GitHub.

2. Clone your fork locally:

    ```
    git clone git@github.com:your_name_here/geemap.git
    ```

3. Install your local copy into a conda env. Assuming you have conda installed, this is how you set up your fork for local development:

    ```
    conda create -n geemap-test python
    ```

    ```
    conda activate geemap-test
    ```
    
    ```
    cd geemap/
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
    flake8 geemap tests
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
3. The pull request should work for Python 3.6, 3.7 and 3.8, and for PyPy. Check <https://github.com/giswqs/geemap/actions> and make sure that the tests pass for all supported Python versions.
# river module

::: geemap.algorithms.river
# legends module

::: geemap.legends# foliumap module

::: geemap.foliumap
# heremap module

::: geemap.heremap
# Get Started

This Get Started guide is intended as a quick way to start programming with **geemap** and the Earth Engine Python API.

## Launch Jupyter notebook

    conda activate gee
    jupyter notebook

## Import libraries

    import ee
    import geemap

## Create an interactive map

    Map = geemap.Map(center=(40, -100), zoom=4)
    Map

## Add Earth Engine data

    # Add Earth Engine dataset
    dem = ee.Image('USGS/SRTMGL1_003')
    landcover = ee.Image("ESA/GLOBCOVER_L4_200901_200912_V2_3").select('landcover')
    landsat7 = ee.Image('LE7_TOA_5YEAR/1999_2003')
    states = ee.FeatureCollection("TIGER/2018/States")

## Set visualization parameters

    dem_vis = {
    'min': 0,
    'max': 4000,
    'palette': ['006633', 'E5FFCC', '662A00', 'D8D8D8', 'F5F5F5']}

    landsat_vis = {
        'min': 20,
        'max': 200,
        'bands': ['B4', 'B3', 'B2']
    }

## Display data on the map

    Map.addLayer(dem, dem_vis, 'SRTM DEM', True, 0.5)
    Map.addLayer(landcover, {}, 'Land cover')
    Map.addLayer(landsat7, landsat_vis, 'Landsat 7')
    Map.addLayer(states, {}, "US States")

## Interact with the map

Once data are added to the map, you can interact with data using various tools, such as the drawing tools, inspector tool, plotting tool.
Check the video below on how to use the Inspector tool to query Earth Engine interactively.

[![geemap](http://img.youtube.com/vi/k477ksjkaXw/0.jpg)](http://www.youtube.com/watch?v=k477ksjkaXw "inspector")
# geemap cheat sheet

## Installation

### Install from PyPI

```console
pip install geemap
```

### Install from conda-forge

```console
conda install geemap -c conda-forge
```

### Create a new conda env

```console
conda create -n gee python=3.8
conda activate gee
conda install mamba -c conda-forge
mamba install geemap -c conda-forge
```

## Upgrade

### Upgrade from PyPI

```console
pip install -U geemap
```

### Upgrade from conda-forge

```console
conda update geemap -c conda-forge
```

### Upgrade from GitHub

```console
import geemap
geemap.update_package()
```

## Map

### Create an interactive map

```python
Map = geemap.Map(center=(lon, lat), zoom=4)
Map
```

### Change the default basemap

```python
Map = geemap.Map(basemap='HYBRID')
```

### Add basemaps

```python
Map.add_basemap('OpenTopoMap')
```

### Add XYZ layers

```python
url = 'https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}'
Map.add_tile_layer(url, name='Google Satellite', attribution='Google')
```

### Add WMS layers

```python
url = 'https://services.nationalmap.gov/arcgis/services/USGSNAIPImagery/ImageServer/WMSServer?'
Map.add_wms_layer(url=url, layers='0', name='NAIP Imagery', format='image/png', shown=True)
```

### Add Earth Engine layers

```python
image = ee.Image('USGS/SRTMGL1_003')
vis_params = {
  'min': 0,
  'max': 4000,
  'palette': ['006633', 'E5FFCC', '662A00', 'D8D8D8', 'F5F5F5']
  }
Map.addLayer(image, vis_params, 'SRTM DEM', True, 0.5)
```

### Set map center

```python
Map.setCenter(lon, lat, zoom)
```

### Center map around an object

```python
Map.centerObject(ee_object, zoom)
```

### Add built-in legends

```python
Map.add_legend(builtin_legend='NLCD')
```

### Add custom legends

```python
Map.add_legend(legend_title, legend_dict, layer_name)
```

## Export data

### Export vector to local

```python
geemap.ee_to_shp(ee_object, filename)
geemap.ee_export_geojson(ee_object, filename)
geemap.ee_export_vector(ee_object, filename)
```

### Export vector to Google Drive

```python
ee_export_vector_to_drive(ee_object, description, folder, file_format='shp', selectors=None)
```

### Export image to local

```python
ee_export_image(ee_object, filename, scale=None, crs=None, region=None, file_per_band=False)
```

### Export image collection to local

```python
ee_export_image_collection(ee_object, out_dir, scale=None, crs=None, region=None, file_per_band=False)
```

### Export image to Google Drive

```python
ee_export_image_to_drive(ee_object, description, folder=None, region=None, scale=None, crs=None, file_format='GeoTIFF')
```

### Export image collection to Google Drive

```python
ee_export_image_collection_to_drive(ee_object, descriptions=None, folder=None, region=None, scale=None, crs=None, file_format='GeoTIFF')
```
# osm module

::: geemap.osm
# colormaps module

::: geemap.colormaps# cartoee module

::: geemap.cartoee# Tutorials

## YouTube Channel

More video tutorials for geemap and Earth Engine are available on my [YouTube channel](https://www.youtube.com/c/QiushengWu). If you can't access YouTube in your country, you can try [西瓜视频](http://gishub.org/xigua) or [哔哩哔哩](https://space.bilibili.com/527404442)。

[![Earth Engine Tutorials on YouTube](https://wetlands.io/file/images/youtube.png)](https://www.youtube.com/c/QiushengWu)

## geemap Tutorials

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://gishub.org/geemap-colab)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/geemap-binder)
[![image](https://mybinder.org/badge_logo.svg)](https://gishub.org/geemap-binder)

1. [Introducing the geemap Python package for interactive mapping with Google Earth Engine](#1-introducing-the-geemap-python-package-for-interactive-mapping-with-google-earth-engine) ([video](https://youtu.be/h0pz3S6Tvx0) | [gif](https://i.imgur.com/pI39k7v.gif) | [notebook](https://geemap.org/notebooks/01_geemap_intro))
2. [Using basemaps in geemap and ipyleaflet for interactive mapping with Google Earth Engine](#2-using-basemaps-in-geemap-and-ipyleaflet-for-interactive-mapping-with-google-earth-engine) ([video](https://youtu.be/6J5ZCIUPXfI) | [gif](https://i.imgur.com/P5B2f7p.gif) | [notebook](https://geemap.org/notebooks/02_using_basemaps))
3. [Introducing the Inspector tool for Earth Engine Python API](#3-introducing-the-inspector-tool-for-earth-engine-python-api) ([video](https://youtu.be/k477ksjkaXw) | [gif](https://i.imgur.com/8d77gtI.gif) | [notebook](https://geemap.org/notebooks/03_inspector_tool))
4. [Creating a split-panel map for visualizing Earth Engine data](#4-creating-a-split-panel-map-for-visualizing-earth-engine-data) ([video](https://youtu.be/9EUTX8j-YVM) | [gif](https://i.imgur.com/kql7pC3.gif) | [notebook](https://geemap.org/notebooks/04_split_panel_map))
5. [Using drawing tools to interact with Earth Engine data](#5-using-drawing-tools-to-interact-with-earth-engine-data) ([video](https://youtu.be/N7rK2aV1R4c) | [gif](https://i.imgur.com/Lm5pDUr.gif) | [notebook](https://geemap.org/notebooks/05_drawing_tools))
6. [Creating an interactive map with a marker cluster](#6-creating-an-interactive-map-with-a-marker-cluster) ([video](https://youtu.be/4HycJPrwpuo) | [gif](https://i.imgur.com/GF4cOqh.gif) | [notebook](https://geemap.org/notebooks/06_marker_cluster))
7. [Converting data formats between GeoJSON and Earth Engine](#7-converting-data-formats-between-geojson-and-earth-engine) ([video](https://youtu.be/DbK_SRgrCHw) | [gif](https://i.imgur.com/hVPmUG1.gif) | [notebook](https://geemap.org/notebooks/07_geojson))
8. [Automated conversion from Earth Engine JavaScripts to Python scripts and Jupyter notebooks](#8-automated-conversion-from-earth-engine-javascripts-to-python-scripts-and-jupyter-notebooks) ([video](https://youtu.be/RpIaalFk4H8) | [gif](https://i.imgur.com/BW0zJnN.gif) | [notebook](https://geemap.org/notebooks/08_ee_js_to_ipynb))
9. [Interactive plotting of Earth Engine data with minimal coding](#9-interactive-plotting-of-earth-engine-data-with-minimal-coding) ([video](https://youtu.be/PDab8mkAFL0) | [gif](https://i.imgur.com/iGMRnRb.gif) | [notebook](https://geemap.org/notebooks/09_plotting))
10. [Using shapefiles with Earth Engine without having to upload data to GEE](#10-using-shapefiles-with-earth-engine-without-having-to-upload-data-to-gee) ([video](https://youtu.be/OlNlqfj4uHo) | [gif](https://i.imgur.com/W3vNSdX.gif) | [notebook](https://geemap.org/notebooks/10_shapefiles))
11. [Exporting Earth Engine Image and ImageCollection as GeoTIFF and Numpy array](#11-exporting-earth-engine-image-and-imagecollection-as-geotiff-and-numpy-array) ([video](https://youtu.be/_6JOA-iiEGU) | [gif](https://i.imgur.com/RonLr0j.gif) | [notebook](https://geemap.org/notebooks/11_export_image))
12. [Computing zonal statistics with Earth Engine and exporting results as CSV or shapefile](#12-computing-zonal-statistics-with-earth-engine-and-exporting-results-as-csv-or-shapefile) ([video](https://youtu.be/ou-Xm3CLitM) | [gif](https://i.imgur.com/8xmUitW.gif) | [notebook](https://geemap.org/notebooks/12_zonal_statistics))
13. [Calculating zonal statistics by group (e.g., analyzing land cover composition of each country/state)](#13-calculating-zonal-statistics-by-group-eg-analyzing-land-cover-composition-of-each-countrystate) ([video](https://youtu.be/cORcGGH03gg) | [gif](https://i.imgur.com/LxD2em9.gif) | [notebook](https://geemap.org/notebooks/13_zonal_statistics_by_group))
14. [Adding a customized legend for Earth Engine data](#14-adding-a-customized-legend-for-earth-engine-data) ([video](https://youtu.be/NwnW_qOkNRw) | [gif](https://i.imgur.com/idkZHQp.gif) | [notebook](https://geemap.org/notebooks/14_legends))
15. [Converting Earth Engine JavaScripts to Python code directly within Jupyter notebook](#15-converting-earth-engine-javascripts-to-python-code-directly-within-jupyter-notebook) ([video](https://youtu.be/nAzZjKKd4w0) | [gif](https://i.imgur.com/aGCBWSV.gif) | [notebook](https://geemap.org/notebooks/15_convert_js_to_py))
16. [Adding animated text to GIF images generated from Earth Engine data](#16-adding-animated-text-to-gif-images-generated-from-earth-engine-data) ([video](https://youtu.be/fDnDVuM_Ke4) | [gif](https://i.imgur.com/MSde1om.gif) | [notebook](https://geemap.org/notebooks/16_add_animated_text))
17. [Adding colorbar and images to GIF animations generated from Earth Engine data](#17-adding-colorbar-and-images-to-gif-animations-generated-from-earth-engine-data) ([video](https://youtu.be/CpT3LQPNKJs) | [gif](https://i.imgur.com/13tFMSI.gif) | [notebook](https://geemap.org/notebooks/17_add_colorbar_to_gif))
18. [Creating Landsat timelapse animations with animated text using Earth Engine](#18-creating-landsat-timelapse-animations-with-animated-text-using-earth-engine) ([video](https://youtu.be/OwjSJnGWKJs) | [gif](https://i.imgur.com/XOHOeXk.gif) | [notebook](https://geemap.org/notebooks/18_create_landsat_timelapse))
19. [How to search and import datasets from Earth Engine Data Catalog](#19-how-to-search-and-import-datasets-from-earth-engine-data-catalog) ([video](https://youtu.be/lwtgzrHrXj8) | [gif](https://i.imgur.com/E09p64F.gif) | [notebook](https://geemap.org/notebooks/19_search_places_and_datasets))
20. [Using timeseries inspector to visualize landscape changes over time](#20-using-timeseries-inspector-to-visualize-landscape-changes-over-time) ([video](https://youtu.be/0CZ7Aj8hCyo) | [gif](https://i.imgur.com/61wbRjK.gif) | [notebook](https://geemap.org/notebooks/20_timeseries_inspector))
21. [Exporting Earth Engine maps as HTML files and PNG images](#21-exporting-earth-engine-maps-as-html-files-and-png-images) ([video](https://youtu.be/GWMvaNQz3kY) | [gif](https://i.imgur.com/rJuXH4a.gif) | [notebook](https://geemap.org/notebooks/21_export_map_to_html_png))
22. [How to import Earth Engine Python scripts into Jupyter notebook?](#22-how-to-import-earth-engine-python-scripts-into-jupyter-notebook) ([video](https://youtu.be/V7CbB9W41w8) | [gif](https://i.imgur.com/WwJoBHF.gif) | [notebook](https://geemap.org/notebooks/22_import_scripts))
23. [How to search Earth Engine API and import assets from GEE personal account?](#23-how-to-search-earth-engine-api-and-import-assets-from-gee-personal-account) ([video](https://youtu.be/c9VJ_uRYSkw) | [gif](https://i.imgur.com/b1auzkr.gif) | [notebook](https://geemap.org/notebooks/22_import_assets))
24. [How to publish interactive Earth Engine maps?](#24-how-to-publish-interactive-earth-engine-maps) ([video](https://youtu.be/NNrrLBIqroY) | [gif](https://i.imgur.com/Hpfzazk.gif) | [notebook](https://geemap.org/notebooks/24_publish_maps))
25. [How to load local raster datasets with geemap?](#25-how-to-load-local-raster-datasets-with-geemap) ([video](https://youtu.be/6XIehAnoazk) | [gif](https://i.imgur.com/nsqEt2O.gif) | [notebook](https://geemap.org/notebooks/25_load_rasters))
26. [How to create and deploy Earth Engine Apps using Python?](https://i.imgur.com/Hpfzazk.gif) ([video](https://youtu.be/nsIjfD83ggA) | [gif](https://i.imgur.com/Hpfzazk.gif) | [notebook](https://geemap.org/notebooks/26_heroku))
27. [How to create an interactive Earth Engine App for creating Landsat timelapse?](https://i.imgur.com/doHfnKp.gif) ([video](https://youtu.be/whIXudC6r_s) | [gif](https://i.imgur.com/doHfnKp.gif) | [notebook](https://geemap.org/notebooks/27_timelapse_app))
28. [How to use your local computer as a web server for hosting Earth Engine Apps?](https://i.imgur.com/q0sJSyi.gif) ([video](https://youtu.be/eRDZBVJcNCk) | [gif](https://i.imgur.com/q0sJSyi.gif) | [notebook](https://geemap.org/notebooks/28_voila))
29. [How to use pydeck for rendering Earth Engine data](https://i.imgur.com/HjFB95l.gif) ([video](https://youtu.be/EIkEH4okFF4) | [gif](https://i.imgur.com/HjFB95l.gif) | [notebook](https://geemap.org/notebooks/29_pydeck))
30. [How to get image basic properties and descriptive statistics](https://i.imgur.com/3B6YhkI.gif) ([video](https://youtu.be/eixBPPWgWs8) | [gif](https://i.imgur.com/3B6YhkI.gif) | [notebook](https://geemap.org/notebooks/30_image_props_stats))
31. [Machine Learning with Earth Engine - Unsupervised Classification](https://i.imgur.com/uNQfrFx.gif) ([video](https://youtu.be/k9MEy2awVJQ) | [gif](https://i.imgur.com/uNQfrFx.gif) | [notebook](https://geemap.org/notebooks/31_unsupervised_classification))
32. [Machine Learning with Earth Engine - Supervised Classification](https://i.imgur.com/jJ2Xiu6.gif) ([video](https://youtu.be/qWaEfgWi21o) | [gif](https://i.imgur.com/jJ2Xiu6.gif) | [notebook](https://geemap.org/notebooks/32_supervised_classification))
33. [Machine Learning with Earth Engine - Performing Accuracy Assessment for Image Classification](https://i.imgur.com/1JkIrF3.gif) ([video](https://youtu.be/JYptiw-I8dc) | [gif](https://i.imgur.com/1JkIrF3.gif) | [notebook](https://geemap.org/notebooks/33_accuracy_assessment))
34. [Interactive extraction of pixel values and interactive region reduction](https://i.imgur.com/LXRqSTu.gif) ([video](https://t.co/D0NC63KgF3) | [gif](https://i.imgur.com/LXRqSTu.gif) | [notebook](https://geemap.org/notebooks/34_extract_values))
35. How to use geemap and Earth Engine in Google Colab ([video](https://youtu.be/fG6kx9vq7hs) | [gif](https://i.imgur.com/OJCasMe.gif) | [notebook](https://geemap.org/notebooks/35_geemap_colab))
36. How to find out the greenest day of the year ([video](https://youtu.be/9KEaW4Ks5fQ) | [gif](https://i.imgur.com/eLDeb4t.gif) | [notebook](https://geemap.org/notebooks/36_quality_mosaic))
37. How to use Earth Engine with pydeck for 3D terrain visualization ([video](https://youtu.be/4E3zOP3-md8) | [gif](https://i.imgur.com/Gx7Y015.gif) | [notebook](https://geemap.org/notebooks/37_pydeck_3d))
38. How to use Cloud Optimized GeoTIFF with Earth Engine ([video](https://youtu.be/2P2PGSMj-wM) | [gif](https://i.imgur.com/z2mfrrZ.gif) | [notebook](https://geemap.org/notebooks/38_cloud_geotiff))
39. How to create Landsat timelapse animations without coding ([video](https://youtu.be/ab0oUhnd_7U) | [gif](https://i.imgur.com/7eyMcZQ.gif) | [notebook](https://geemap.org/notebooks/39_timelapse))
40. How to add interactive widgets to the map ([video](https://youtu.be/KsIxGq6cHtw) | [gif](https://i.imgur.com/peRZZjj.gif) | [notebook](https://geemap.org/notebooks/40_ipywidgets))
41. How to develop an Earth Engine app for mapping surface water dynamics ([video](https://youtu.be/fHdwV3LEMYo) | [gif](https://i.imgur.com/GUWSVZs.gif) | [notebook](https://geemap.org/notebooks/41_water_app))
42. How to upload data to Earth Engine Apps using ipywidgets ([video](https://youtu.be/4-WeaiObj84) | [gif](https://i.imgur.com/INLzqdw.gif) | [notebook](https://geemap.org/notebooks/42_upload_data))
43. How to extract pixel values from an Earth Engine image using a point shapefile ([video](https://youtu.be/UbQ8jyc4VP4) | [gif](https://i.imgur.com/pbt6neQ.gif) | [notebook](https://geemap.org/notebooks/43_extract_values_to_points))
44. How to use Cloud Optimized GeoTIFF (COG) and SpatioTemporal Asset Catalog (STAC) ([video](https://youtu.be/yLlYoy01RxA) | [gif](https://i.imgur.com/XjG3zYq.gif) | [notebook](https://geemap.org/notebooks/44_cog_stac))
45. How to load a virtual mosaic of Cloud Optimized GeoTIFFs (COG) ([video](https://youtu.be/jDUaopr0Dhg) | [gif](https://i.imgur.com/My8Ksh7.gif) | [notebook](https://geemap.org/notebooks/45_cog_mosaic))
46. How to use locally trained machine learning models with Earth Engine ([video](https://youtu.be/nq_Ro7E0b6E) | [gif](https://i.imgur.com/muwDfkC.gif) | [notebook](https://geemap.org/notebooks/46_local_rf_training))
47. How to download image thumbnails from Earth Engine ([video](https://youtu.be/qwXZDSbfyE8) | [gif](https://i.imgur.com/gqr7CNz.gif) | [notebook](https://geemap.org/notebooks/47_image_thumbnails))
48. How to add a draggable legend to folium maps ([video](https://youtu.be/-rO1MztlLMo) | [gif](https://i.imgur.com/i2Bye9X.gif) | [notebook](https://geemap.org/notebooks/48_folium_legend))
49. How to add a colorbar to the map ([video](https://youtu.be/qiKns09X1Ao) | [gif](https://i.imgur.com/VpMq8M9.gif) | [notebook](https://geemap.org/notebooks/49_colorbar))
50. How to create publication quality maps using cartoee ([video](https://youtu.be/t24_lpYA1ko) | [gif](https://i.imgur.com/fwCzZTi.gif) | [notebook](https://geemap.org/notebooks/50_cartoee_quickstart))
51. How to create publication quality maps with custom projections ([video](https://youtu.be/3dS2EkAuAxM) | [gif](https://i.imgur.com/vvvF94j.gif) | [notebook](https://geemap.org/notebooks/51_cartoee_projections))
52. How to create timelapse animations with custom projection, scale bar, and north arrow ([video](https://youtu.be/ejuugljSut4) | [gif](https://i.imgur.com/MVQFyHN.gif) | [notebook](https://geemap.org/notebooks/52_cartoee_gif))
53. How to change layer visualization interactively with a GUI ([video](https://youtu.be/4E7gg6yaHBg) | [gif](https://i.imgur.com/VqqlMSK.gif) | [notebook](https://geemap.org/notebooks/53_layer_vis))
54. Visualizing Earth Engine vector data interactively with a GUI ([video](https://youtu.be/SIMnvbn8d-4) | [gif](https://youtu.be/F7xa5OaweY0) | [notebook](https://geemap.org/notebooks/54_vector_vis))
55. Visualizing Earth Engine raster data interactively with a GUI ([video](https://youtu.be/2R3933NFIa0) | [gif](https://youtu.be/6HFGvyXOXJM) | [notebook](https://geemap.org/notebooks/55_raster_vis))
56. Loading local vector and raster data to geemap without coding ([video](https://youtu.be/Zhwz0uS4Xi0) | [gif](https://youtu.be/SWLpnYnsqMw) | [notebook](https://geemap.org/notebooks/56_local_data))
57. Creating publication-quality maps with multiple Earth Engine layers ([video](https://youtu.be/v-FWj9dAMJ8) | [gif](https://youtu.be/85Cu3cVLmOY) | [notebook](https://geemap.org/notebooks/57_cartoee_blend))
58. Loading vector data (e.g., shp, kml, geojson) to the map without coding ([video](https://youtu.be/10KA7uhEWUM) | [gif](https://youtu.be/UsmigaIDNpE) | [notebook](https://geemap.org/notebooks/58_add_vector))
59. Using whitebox with geemap ([video](https://youtu.be/n8ODeZpuyCE) | gif | [notebook](https://geemap.org/notebooks/59_whitebox))
60. Visualizing Earth Engine data with over 200 colormaps through dot notation ([video](https://youtu.be/RBCf7wgK3Cg) | gif | [notebook](https://geemap.org/notebooks/60_colormaps))
61. Adding a scale bar to a cartoee map (video | gif | [notebook](https://geemap.org/notebooks/61_cartoee_scalebar))
62. Using the time slider for visualizing Earth Engine time-series images ([video](https://youtu.be/w_nWkNz8fyI) | gif | [notebook](https://geemap.org/notebooks/62_time_slider))
63. Creating interactive charts from Earth Engine data (video | [gif](https://youtu.be/e-GTdUUc8N8) | [notebook](https://geemap.org/notebooks/63_charts))
64. Accessing the Earth Engine Data Catalog via dot notation with autocompletion (video | [gif](https://youtu.be/hGbs2cl7otk) | [notebook](https://geemap.org/notebooks/64_data_catalog))
65. Styling Earth Engine vector data (video | gif | [notebook](https://geemap.org/notebooks/65_vector_styling))
66. Adding a legend to publication quality maps using cartoee (video | gif | [notebook](https://geemap.org/notebooks/66_cartoee_legend))
67. Creating training samples for machine learning and supervised image classification (video | [gif](https://youtu.be/VWh5PxXPZw0) | [notebook](https://geemap.org/notebooks/67_training_samples))
68. Converting NetCDF to Earth Engine Image (video | gif | [notebook](https://geemap.org/notebooks/68_netcdf_to_ee))
69. Plotting Earth Engine vector data with cartoee (video | [gif](https://youtu.be/Gr6GBuBWnnk) | [notebook](https://geemap.org/notebooks/69_cartoee_vector))
70. Creating linked maps with a few lines of code (video | [gif](https://youtu.be/AFUGje3VWM8) | [notebook](https://geemap.org/notebooks/70_linked_maps))
71. Creating Landsat timelapse animations with a few clicks (video | [gif](https://youtu.be/mA21Us_3m28) | [notebook](https://geemap.org/notebooks/71_timelapse))
72. Creating time-series cloud-free composites with a few clicks (video | [gif](https://youtu.be/kEltQkNia6o) | [notebook](https://geemap.org/notebooks/72_time_slider_gui))
73. Generating transects along lines with Earth Engine without coding (video | [gif](https://youtu.be/0TNXSs6fwrg) | [notebook](https://geemap.org/notebooks/73_transect))
74. Creating points from CSV without coding (video | gif | [notebook](https://geemap.org/notebooks/74_csv_to_points))
75. Visualizing land cover change with inteactive Sankey diagrams (video | gif | [notebook](https://geemap.org/notebooks/75_sankee))
76. Downloading and visualizing OpenStreetMap data (video | gif | [notebook](https://geemap.org/notebooks/76_osm_to_ee))
77. Adding Planet global monthly and quarterly mosaic (video | gif | [notebook](https://geemap.org/notebooks/77_planet_imagery))
78. Using timeseries inspector with one click (video | [gif](https://i.imgur.com/s1GoEOV.gif) | [notebook](https://geemap.org/notebooks/78_ts_inspector))
79. Creating histograms using the geemap chart module (video | gif | [notebook](https://geemap.org/notebooks/79_chart_histogram/))
80. Adding a point layer with popup attributes (video | gif | [notebook](https://geemap.org/notebooks/80_point_layer))
81. Creating timelapse animations from GEOS weather satellites (video | gif | [notebook](https://geemap.org/notebooks/81_goes_timelapse))
82. Creating elevation contours for any location around the globe (video | gif | [notebook](https://geemap.org/notebooks/82_contours))
83. Loading local raster datasets and Cloud Optimized GeoTIFF (COG) ([notebook](https://geemap.org/notebooks/83_local_tile))
84. Downloading OpenStreetMap data with a single line of code ([notebook](https://geemap.org/notebooks/84_openstreetmap))
85. Converting PostGIS data to ee.FeatureCollection ([notebook](https://geemap.org/notebooks/85_postgis))
86. Adding image overlay to the map ([notebook](https://geemap.org/notebooks/86_image_overlay))
87. Adding points from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://geemap.org/notebooks/87_add_points_from_xy))
88. Adding circle markers from xy data (e.g., CSV, Pandas DataFrame) ([notebook](https://geemap.org/notebooks/88_circle_markers))
89. Labeling Earth Engine FeatureCollection on the map ([notebook](https://geemap.org/notebooks/89_add_labels))
90. Creating 1-m resolution NAIP imagery timelapse ([notebook](https://geemap.org/notebooks/90_naip_timelapse))
91. Adding Planetary Computer STAC item to the map ([notebook](https://geemap.org/notebooks/91_planetary_computer))
92. Using plotly with Earth Engine ([notebook](https://geemap.org/notebooks/92_plotly))
93. Getting pixel values from COG/STAC using the Inspector tool ([notebook](https://leafmap.org/notebooks/93_cog_inspector))
94. Using heremap with Earth Engine ([notebook](https://geemap.org/notebooks/94_heremap))
95. Creating Cloud Optimized GeoTIFF (COG) ([notebook](https://geemap.org/notebooks/95_create_cog))
# FAQ

## How do I report an issue or make a feature request

Please go to <https://github.com/giswqs/geemap/issues>.

## How do I cite geemap in publications

Wu, Q., (2020). geemap: A Python package for interactive mapping with Google Earth Engine. _The Journal of Open Source Software_, 5(51), 2305. <https://doi.org/10.21105/joss.02305>

```
Bibtex:
@article{wu2020geemap,
    title={geemap: A Python package for interactive mapping with Google Earth Engine},
    author={Wu, Qiusheng},
    journal={Journal of Open Source Software},
    volume={5},
    number={51},
    pages={2305},
    year={2020}
}
```

## Why does geemap use two plotting backends: folium and ipyleaflet

A key difference between [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) and [folium](https://github.com/python-visualization/folium) is that ipyleaflet is built upon ipywidgets and allows bidirectional communication between the front-end and the backend enabling the use of the map to capture user input, while folium is meant for displaying static data only ([source](https://blog.jupyter.org/interactive-gis-in-jupyter-with-ipyleaflet-52f9657fa7a)). Note that [Google Colab](https://colab.research.google.com/) currently does not support ipyleaflet ([source](https://github.com/googlecolab/colabtools/issues/60#issuecomment-596225619)). Therefore, if you are using geemap
with Google Colab, you should use `import geemap.foliumap as geemap`. If you are using geemap with a local Jupyter notebook server, you can
use `import geemap`, which provides more functionalities for capturing user input (e.g., mouse-clicking and moving).

## Why the interactive map does not show up

If the interactive map does not show up on Jupyter notebook and JupyterLab, it is probably because the ipyleaflet extentsion is not installed properly.
For Jupyter notebook, try running the following two commands within your geemap conda environment:

```
jupyter nbextension install --py --symlink --sys-prefix ipyleaflet
jupyter nbextension enable --py --sys-prefix ipyleaflet
```

For JupyterLab, try running the following command within your geemap conda environment:

```
jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-leaflet

```

## How to use geemap in countries where Google Services are blocked

If you are trying to use geemap in coutries where Gooogle Services are blocked (e.g., China), you will need a VPN. Use `geemap.set_proxy(port=your-port-number)` to connect to Earth Engine servers. Otherwise, you might encounter a connection timeout issue.

```
import geemap
geemap.set_proxy(port=your-port-number)
Map = geemap.Map()
Map
```

![](https://i.imgur.com/AHY9rT2.png)
# datasets module

::: geemap.datasets
# geemap module

::: geemap.geemap# Citations

If you find **geemap** useful in your research, please consider citing the following papers to support my work. Thank you for your support.

* Wu, Q., (2020). geemap: A Python package for interactive mapping with Google Earth Engine. The Journal of Open Source Software, 5(51), 2305. <https://doi.org/10.21105/joss.02305>
* Wu, Q., Lane, C. R., Li, X., Zhao, K., Zhou, Y., Clinton, N., DeVries, B., Golden, H. E., & Lang, M. W. (2019). Integrating LiDAR data and multi-temporal aerial imagery to map wetland inundation dynamics using Google Earth Engine. Remote Sensing of Environment, 228, 1-13. <https://doi.org/10.1016/j.rse.2019.04.015> ([pdf](https://gishub.org/2019_rse) | [source code](https://doi.org/10.6084/m9.figshare.8864921))