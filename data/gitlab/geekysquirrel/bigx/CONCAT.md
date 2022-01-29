[![GPL-3.0](https://img.shields.io/badge/license-GPL3-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html)&nbsp;
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02537/status.svg)](https://doi.org/10.21105/joss.02537)
[![Branch coverage status](https://geekysquirrel.gitlab.io/bigx/reports/coverage/branch_coverage.svg)](https://geekysquirrel.gitlab.io/bigx/reports/coverage/lcov-report/)&nbsp;
[![Statement coverage status](https://geekysquirrel.gitlab.io/bigx/reports/coverage/statement_coverage.svg)](https://geekysquirrel.gitlab.io/bigx/reports/coverage/lcov-report/)&nbsp;
[![Quality Gate status](https://sonarcloud.io/api/project_badges/measure?project=geekysquirrel_bigx&metric=alert_status)](https://sonarcloud.io/dashboard?id=geekysquirrel_bigx)&nbsp;
[![Log in to see pipeline status](https://gitlab.com/geekysquirrel/bigx/badges/dev/pipeline.svg)](https://gitlab.com/geekysquirrel/bigx/-/pipelines/latest)

# Welcome to BigX

![](src/img/bigx.png)

[![Demo](doc/img/demo.png)](https://geekysquirrel.gitlab.io/bigx/)
[![Setup](doc/img/setup.png)](doc/setup.md)
[![Usage](doc/img/usage.png)](doc/layer-reference.md)

# About

BigX is a tool to create custom maps using open data (GeoJSON and CSV).

Its main features are:

* Can be configured to work completely offline
* Layer configuration files make it easy to load different datasets
  - Definition works in a purely declarative way - no code needs to be written!
  - all layer definitions are text-based so configuration can be done without access to a GUI
    and can easily be scripted, e.g. for use in containerised or cloud environments
* Load any geographical data in (Geo)JSON or CSV format and combine your data with other datasets
* Dynamically load data instead of pre-loading everything to reduce waiting time

**Who this is for:**

* People with some computer literacy and basic understanding of algorithms and data structures but not necessarily any coding skills
* People who have some data and want to quickly visualise it and combine it without having to write code
* People who want to share interactive maps with a restricted audience while controlling the physical location of the data
* The resulting maps, if built/deployed with bundled data by somebody who *knows what they're doing*™ are for anybody who can browse a map.

**Who this is *not* for:**

* Complete beginners who have never heard of things such as map data structures, JSON or CSV
* People aiming to do their entire data analysis using this tool.
  - While it can be used to manipulate and filter, the more of it you do, the slower it will get.
* People aiming to visualise "Big Data" (whatever you imagine when you hear that term).
  - If you would like to get good performance, it's a good idea to preprocess your data using
    a more suitable tool such as Python. Have a look at the [sample script](./scripts/preprocessUtil.py) to get some ideas.

**Alternatives:**

* You could build your own custom map from scratch using Python and its many libraries.
* If you would like a GUI desktop-based approach to visualise your data, you could try [QGIS](https://www.qgis.org/) or [uDig](http://udig.refractions.net/).
* Paid hosted services for interactive maps are available from places such as [ArcGIS](https://www.arcgis.com/index.html) or [Carto](https://carto.com/pricing/).
* Maps can be coded using [R](https://www.r-project.org/)/[Shiny](https://shiny.rstudio.com/) and compiled down into Javascript.


Before you start please note:

BigX does not contain data by default with the exception of some example data.
To use it as it was intended, you will have to acquire the data yourself and responsibility for using the proper credits lies **entirely with you**! I will not take responsibility if you use data without permission or fail to give credit to the data owners. I do not own any of the data shown in this README but an effort has been made to credit where necessary. If you're a data owner and you think I've used your data inappropriately, please get in touch.

# Demo

You can find a running demo instance [here](https://geekysquirrel.gitlab.io/bigx/).

This is not supposed to be a "production" version of BigX, just a demo. If you like it please download it and deploy it yourself.

**Why is it slow?**
As described in the "Architecture" section below, BigX runs in the browser on your computer ("client-side"). This means, the data needs to first get to you from the server where it is stored. If you run BigX on your local machine as described in the [the setup instructions](doc/setup.md), these two machines are the same, so the data gets sent to your browser *really fast™*. That's why ultimately for serious use, such as exploratory data analysis, you should download BigX.

# BigX in action

Different types of layers can be combined into one map. The map key is dynamically updated depending on which layers are shown. *(© OpenStreetMap contributors under [ODbL](http://www.openstreetmap.org/copyright). Contains OS data © Crown copyright and database right 2018)*:

![Different layers on one map](src/img/demo/dynamic-legend.gif)

CSV data can be added on top of geographic data to generate complex maps *(© Stadia Maps, OpenMapTiles. Data by OpenStreetMap, under ODbL. Contains OS data © Crown copyright and database right 2018, election results © UK Parliament 2017)*:

![CSV data attached to geographical data](src/img/demo/csv.png)

All layers can be combined to give context to the data and create meaningful maps. The colour scheme can be generated using the features' data. *(© Stadia Maps, OpenMapTiles. Data by OpenStreetMap, under ODbL. Contains Ordnance Survey, Office of National Statistics, National Records Scotland and LPS Intellectual Property data © Crown copyright and database right 2016. Data licensed under the terms of the Open Government Licence. Ordnance Survey data covered by OS OpenData Licence.)*:

![Combined map with generated colour spectrum](src/img/demo/generated-colour.png)


Use satellite images to explore what's going on. *(© Esri, DigitalGlobe, GeoEye, i-cubed, USDA FSA, USGS, AEX, Getmapping, Aerogrid, IGN, IGP, swisstopo, and the GIS User Community. Data from Cambodian Mine Action and Victim Assistance Authority (CMAA) and Office for the Coordination of Humanitarian Affairs (OCHA) to Humanitarian Data Exchange (HDX).)*:

![Satellite](src/img/demo/satellite.png)

Generate overlays from points to visually filter the data *(© OpenStreetMap contributors under [ODbL](http://www.openstreetmap.org/copyright). Contains OS data © Crown copyright and database right 2018. © 2009-2012 Atos and RSP)*:

!["Overlay" layers for visual filtering](src/img/demo/overlay.png)

Load large datasets thanks to clustering. Clusters contain content-aware icons *(© OpenStreetMap contributors under [ODbL](http://ww.openstreetmap.org/copyright). Contains OS data © Crown copyright and database right 2018)*:

![Clustered markers with dynamic preview](src/img/demo/clustering.png)

Easily generate heat maps from points-based layers using different styles *(© Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL. © 2018 Food Standards Agency and ONS Postcode Directory, Open Government license)*:

![Heat map](src/img/demo/heat.png)

Now with dark mode! Switch UI elements and base map - be kind to your eyes *(© Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL. © 2018 Food Standards Agency and ONS Postcode Directory, Open Government license)*:

![Hexbin map](src/img/demo/hexbins.png)

Credit where credit is due: all copyright information is dynamically generated depending on what layers are shown:

![Dynamically generated credits](src/img/demo/credits.gif)

Plan data-driven trips using the built-in measurement tool *(© OpenStreetMap contributors under [ODbL](http://www.openstreetmap.org/copyright). © 2018 OpenStreetMap contributors, Contains OS data © Crown copyright and database right 2018, © 2009-2012 Atos and RSP, © 2003-2019 Walking Englishman)*:

![Route planning tools](src/img/demo/measurement.png)

Scrape live web-content to display on the map by combining with geographic data, using e.g. postcodes to map to coordinates *(© Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL. © Property data by Rightmove.com, Contains OS data © Crown copyright and database right. Contains Royal Mail data © Royal Mail copyright and database right © 2020 Source: Office for National Statistics licensed under the Open Government Licence v.3.0 © Geonames.org)*:

![Live layers](src/img/demo/livelayers.png)

Create multiple views that logically group layers into different, themed maps. Easily switch between them while keeping your layers control tidy:

![Views](src/img/demo/views.png)

Bring your GPX files (or any other temporal data) to life by animating them on the map *(© Map tiles by Stamen Design, under CC BY 3.0. Data by OpenStreetMap, under CC BY SA. Contains OS data © Crown copyright and database right 2018. GPX data © 2018 Stefanie Wiegand)*:

![Animated](src/img/demo/animated.gif)

Show movement by importing a directed graph (here: bike trips between rental stations in London); potentially by combining multiple datasets for locations and movement *(© Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL. Powered by TfL Open Data. Contains OS data © Crown copyright and database rights 2016 and Geomni UK Map data © and database rights 2019)*:

![Migration](src/img/demo/migration.gif)

# Architecture

A typical setup of BigX looks like this:

![The typical setup](src/img/demo/arch.png)

BigX is running locally and all data (geographic data, map tiles and other data) is stored on the same machine. The internet conection is optional if you want to use online base maps but none of your data gets transferred anywhere outside your machine.

Other setups are possible, and can be used in scenarios where it's not practical for end-users to store the data, such as the [demo](https://geekysquirrel.gitlab.io/bigx/), where data is stored on the server that hosts BigX and gets transferred to your browser as required.

# Installation, Configuration and Usage

Check out the [docs](doc/README.md) for all documentation.

# Contributing

To suggest new features, report bugs or get help you can raise an issue on the [issue tracker](https://gitlab.com/geekysquirrel/bigx/-/issues). Merge requests are also welcome, provided all tests pass and coverage has not decreased. For more information on how to run/test BigX, please read the [setup instractions](doc/setup.md).

# Copyright and licenses

[BigX](https://gitlab.com/geekysquirrel/bigx) is licensed under [GPL3](https://www.gnu.org/licenses/gpl-3.0.md).
It uses the following libraries and plugins:

* [Leaflet](https://github.com/Leaflet/Leaflet) ([BSD2](https://github.com/Leaflet/Leaflet/blob/master/LICENSE))
  - [Leaflet.GroupedLayerControl](https://github.com/NHellFire/leaflet-groupedlayercontrol) ([MIT](https://github.com/NHellFire/leaflet-groupedlayercontrol/blob/gh-pages/LICENSE.txt))
  - [Leaflet.LinearMeasurement](https://github.com/NLTGit/Leaflet.LinearMeasurement) ([BSD2](https://github.com/NLTGit/Leaflet.LinearMeasurement/blob/master/LICENSE))
  - [L.Control.Zoomslider](https://github.com/kartena/Leaflet.zoomslider) ([BSD2](https://github.com/kartena/Leaflet.zoomslider/blob/master/LICENSE))
  - [Leaflet-SVGIcon](https://github.com/iatkin/leaflet-svgicon) ([MIT](https://github.com/iatkin/leaflet-svgicon/blob/master/LICENSE))
  - [Leaflet MaskCanvas](https://github.com/domoritz/leaflet-maskcanvas) ([MIT](https://github.com/domoritz/leaflet-maskcanvas/blob/master/LICENSE))
  - [Leaflet.markercluster](https://github.com/Leaflet/Leaflet.markercluster) ([MIT](https://github.com/Leaflet/Leaflet.markercluster/blob/master/MIT-LICENCE.txt))
  - [Leaflet.migrationLayer](https://github.com/lit-forest/leaflet.migrationLayer) ([MIT](https://github.com/lit-forest/leaflet.migrationLayer/blob/master/LICENSE))
  - [Leaflet-Slider](https://github.com/Eclipse1979/leaflet-slider)
  - [sidebar-v2](https://github.com/Turbo87/sidebar-v2) ([MIT](https://github.com/Turbo87/sidebar-v2/blob/master/LICENSE))
  - [Leaflet.heat](https://github.com/Leaflet/Leaflet.heat) ([BSD2](https://github.com/Leaflet/Leaflet.heat/blob/gh-pages/LICENSE))
  - [Leaflet.draw](https://github.com/Leaflet/Leaflet.draw) ([MIT](https://github.com/Leaflet/Leaflet.draw/blob/develop/MIT-LICENSE.md))
  - [leaflet-geotiff-2](https://github.com/danwild/leaflet-geotiff-2) ([MIT](https://github.com/danwild/leaflet-geotiff-2/blob/master/LICENSE))
  - [Leaflet.EasyButton](https://github.com/cliffcloud/Leaflet.EasyButton) ([MIT](https://github.com/CliffCloud/Leaflet.EasyButton/blob/master/LICENSE))
  - [leaflet-d3](https://github.com/Asymmetrik/leaflet-d3) ([MIT](https://github.com/Asymmetrik/leaflet-d3/blob/master/LICENSE))
  - [Leaflet.SvgShapeMarkers](https://github.com/rowanwins/Leaflet.SvgShapeMarkers) ([MIT](https://github.com/rowanwins/Leaflet.SvgShapeMarkers/blob/gh-pages/LICENSE))
* Other
  - [jQuery](https://github.com/jquery/jquery) ([MIT](https://github.com/jquery/jquery/blob/master/LICENSE.txt))
  - [TopoJSON](https://github.com/topojson/topojson) ([BSD](https://github.com/topojson/topojson/blob/master/LICENSE.md))
  - [GPXParser.js](https://github.com/Luuka/GPXParser.js) ([MIT](https://github.com/Luuka/GPXParser.js/blob/master/LICENSE))
  - [PapaParse](https://papaparse.com) ([MIT](https://github.com/mholt/PapaParse/blob/master/LICENSE))
  - [html2canvas](https://html2canvas.hertzen.com/) ([MIT](https://github.com/niklasvh/html2canvas/blob/master/LICENSE))
  - [merge-images](https://github.com/lukechilds/merge-images) ([MIT](https://github.com/lukechilds/merge-images/blob/master/LICENSE))
  - [geotiff.js](https://github.com/geotiffjs/geotiff.js) ([MIT](https://github.com/geotiffjs/geotiff.js/blob/master/LICENSE))
  - [plotty](https://github.com/santilland/plotty) ([MIT](https://github.com/santilland/plotty/blob/master/LICENSE))
  - [d3](https://github.com/d3/d3) ([BSD](https://github.com/d3/d3/blob/master/LICENSE))
  - [d3-hexbin](https://github.com/d3/d3-hexbin) ([BSD](https://github.com/d3/d3-hexbin/blob/master/LICENSE))
  - [RainbowVis-JS](https://github.com/anomal/RainbowVis-JS) ([EPL](https://github.com/anomal/RainbowVis-JS/blob/master/license.md))
![](img/bigx.png)

# Welcome to BigX!

Using this tool should be straightforward, but here are the main points:

## Layers

Layers are the building blocks that make up a map in BigX.
Each layer has a type, such as "point" or "heat" which influence what it looks like.

You can select layers by hovering over the layer button in the top right corner of the map:

![Layer icon](img/manual/layer-icon.png)

In the menu you can select any number of layers to be displayed on the map. Deselect a layer to remove it from the map.

![Layer menu](img/manual/layer-menu.png)

Some layers take a while to load if the data they show is large or your connection is slow.
They will display a spinner while they're still loading.

## Views

A view is a group of layers that "belong together". They will likely relate to the same place.

You can select a view using the "Views" icon on the left hand side. Only one view can be selected at the same time.

A view defines a starting point and zoom level and a default base map but you can change that by selecting a different one.

![Views](img/demo/views.png)

## Tools

### Measuring distance

You can use the tape measure to measure distances on the map.
Click once to add a new stop and double click to finish measuring.

![Views](img/manual/measurement.gif)
*(© OpenStreetMap contributors under [ODbL](http://ww.openstreetmap.org/copyright).)*

You can remove the measurement by clicking the "x" on the final distance or click the measurement button.

### Export map

This allows you to select a part of the map to export as an image which you can download.

This is currently experimental and may or may not work.

## Copyright

If you're planning to capture images or videos please note that both the map tiles and the data displayed may be subject to copyright.

To see who owns the copyright you can check the credits section in the bottom left corner.

![Dynamically generated credits](img/manual/credits.gif)

This is dynamically updated depending on what layers and base map you've selected and will only show copyright notices for the currently displayed map.
# Here be dragons!

These scripts are not really meant for end users and tracked in this repository more for myself.
However, they helped me work with the data so other people might find them useful.

Apart from in-code comments, there is no documentation for these scripts and there are currently no plans of adding one. Feel free to contribute one yourself if you think this project could benefit from it.

Ony use these scripts if you know what you're doing and in any case **always take a backup**, if you lose your data, it's your own fault.

![](../src/img/bigx.png)

# Setup

## Using docker

Assuming you've got [docker](https://www.docker.com/) installed, you can run BigX locally using the command below.

```bash
docker run -d --name bigx -p 8080:80 -v $PWD/src/:/usr/local/apache2/htdocs/ httpd:alpine
```

This does the following things:

* ```bash
  docker run -d --name bigx
  ```

  runs BigX in the background, creating a new container with the name "bigx"

* ```bash
  -p 8080:80
  ```

  maps port 80 on the container to port 8080 on your machine.
  You can change it to something else if the port is already taken.

* ```bash
  -v $PWD/src/:/usr/local/apache2/htdocs/
  ```

  mounts the BigX src directory (wherever you put it)
  as the main public directory of the docker container.

* ```bash
  httpd:alpine
  ```

  specifies the docker image that we're using - a slim version of the apache web server.
  
You can then add more data, create and modify views and layers and generally use BigX.
Note that your browser probably caches files, so if you've made a change and refreshed the page but the change doesn't show, you might have to clear the cache. On Chrome, you can easily do that by pressing F12 to enable the console, then right-click the refresh button and choose "Empty Cache and Hard Reload"; then close the console again using F12.
  
To stop the container, you can type

```bash
docker stop bigx
```

You can the restart the container at any time using 

```bash
docker start bigx
```

or you can fully remove the docker container from your system, using

```bash
docker rm bigx
```

although this won't actually remove any of the BigX files as the container only uses them via a virtual volume.

## Using Electron

You can use [Electron](https://www.electronjs.org/) to run BigX "standalone", i.e. without a web server.

To do that, you need to [install Node.js](https://www.electronjs.org/docs/tutorial/development-environment) and [yarn](https://yarnpkg.com/getting-started/install) and then run

```bash
yarn
yarn start
```

to start BigX.
For copyright reasons only a very small number of datasets is contained by default where the copyright holder's terms allowed it. If you would like to explore the demo data, you need to download it first before running the app. Each layer file contains a comment with instructions on how to do this.

## Running the application on an existing webserver

Copy/unpack the files inside /src into your webserver's directory.
For example with Apache2 on Ubuntu, that's
`/var/www/html` by default. Make sure all permissions are set correctly.
Consult your webserver's manual for other options.

# Configuration

## Get GeoJSON files

Download GeoJSON files or convert shapefiles, e.g. using [OGR2OGR](https://www.gdal.org/ogr2ogr.html).

Example: Convert UK National Grid shapefile to WGS84 GeoJSON:

```bash
ogr2ogr -f GeoJSON -s_srs epsg:27700 -t_srs epsg:4326 county_region.json county_region.shp
```

Interesting datasets are e.g.

* Ordnance Survey's [OpenData](https://www.ordnancesurvey.co.uk/opendatadownload/products.html)
* OSM Shapefiles from [Geofabrik](http://download.geofabrik.de/)
  *(make sure to download *.shp.zip)*

Put the files in the data directory.

Since GeoJSON files can be quite large, BigX also supports TopoJSON.
You can convert (and optionally simplify) your files at [MapShaper.org](https://mapshaper.org/) to reduce the file size significantly.

## Get additional data files

Data can be acquired from many sources, e.g. https://data.gov.uk but could also be generated.
Interesting datasets are e.g.

* the RailDeliveryGroup's [Network Rail INSPIRE data](https://data.gov.uk/dataset/1ed6306a-48a4-4186-b5cb-369eafd8adaa/network-rail-inspire-data)
* the UK [House Price Index](https://data.gov.uk/dataset/d3fd3c42-2dab-4ca0-98ad-36a77561dd6c/house-price-index)
* Ofcom [Connected Nations](https://data.gov.uk/dataset/e218662f-2bdb-4e16-b4a8-8c16fdfd52bd/ofcom-connected-nations-previously-called-infrastructure-report-uk-internet-speeds-and-coverage-broadband-wifi-and-mobile)
* the UK [2017 General Election results](https://data.gov.uk/dataset/b77fcedb-4792-4de4-935f-4f344ed4c2c6/general-election-results-2017)

Some of the files might come in different formats and will need to be converted to CSV using a custom parser. You will need some programming skills to do this.

It is also possible to crawl data from publicly accessible websites, such as https://www.walkingenglishman.com/. Make sure to read the terms and conditions for every website you crawl and if unsure ask the owner's permission. Be a nice person and **credit every data source** properly, even if it's free.

Put the files in the data directory.

## Prepare your data

You can use the [preprocessUtil.py](../src/scripts/preprocessUtil.py) script to preprocess your data files:

* Prune large CSV files by removing data
* Group point features in CSV files and obtain their mean coordinates
* Optionally convert CSV files containing points to GeoJSON

You need to adapt the script to match your input files.

## Get map tiles

You *could* use the OpenStreetMap tile server, but if you use this tool heavily, you might want to use your own tiles so you don't steal their resources. See their [tile usage policy](https://operations.osmfoundation.org/policies/tiles/) if you decide to use them anyway. There are free alternatives, see https://wiki.openstreetmap.org/wiki/Tile_servers.

To set up your own tile server, head to https://switch2osm.org/serving-tiles/ to find out how.

Once you've set it up, you could just use it as it is. However, if you would like the application to work offline and load tiles much faster, you could pre-generate your tiles.

To do that, edit the [osm-tile-crawler.py](../src/scripts/osm-tile-crawler.py) script:

* Set the min/max zoom levels for the tiles you want to generate.
  Note that if you generate large maps, you might run into problems regarding disk space and filesystem-specific issues, such as support for symbolic links. Compatible file systems are e.g. EXT2/3/4.
* Define your bounding box so you only generate the tiles you need
* Point the script to your tile server URL
* Set the output directory name for your generated tiles
* Run the script

On an 8 core laptop with 4GB dedicated to the tile server database, it took about 20 hours to generate a map for the UK at zoom levels 1-16.

**🛑 UNDER NO CIRCUMSTANCES CRAWL SOMEBODY ELSE'S TILE SERVER!!! 🛑**

A word on crawled tiles:

Be aware that map size can be a problem. Tiles covering the whole of the UK at zoom levels 6 - 16 is about 30GB in size. While this may fit on a 32GB USB stick, it will not work on a file system such as FAT or NTSF or even EXT4, because the number of files is too large (only up to 2^32 are supported). You could either use lower zoom levels (configure your application accordingly) or a file system that supports the number of nodes you need (Apple's APFS for instance supports up to 2^63 or you could use ZFS or Btrfs).

## Configure the app

Edit the [js/config.js](../src/js/config.js) file by setting up your base layer(s). If you're using a local tile server, it will look something like this:

```js
{
    name: "OSM local server",
    url: "https://127.0.0.1:8080/tileserver/{z}/{y}/{x}.png",
    credits: {
        text: "OpenStreetMap contributors",
        link: "https://www.openstreetmap.org/copyright",
        date: "2018"
    }
},
```

However, if you have generated the tiles locally, all you need to do is point the map source to your tiles directory, relative to the main application directory:

```js
{
    name: "OSM local tiles",
    url: "atlas/{z}/{x}/{y}.png",
    credits: {
        text: "OpenStreetMap contributors",
        link: "https://www.openstreetmap.org/copyright",
        date: "2018"
    }
},
```

## Set up your views

To view any maps, you need to create views containing layers.

The [js/views/example.js](../src/js/views/example.js) file contains documentation describing how to create a view.

Views live inside the [js/views](../src/js/views) folder. The [js/views/default.js](../src/js/views/default.js) file is the default view. It's empty by default but you can add your own layers to it.

To define multiple views, you can create multiple view files and configure BigX to use them in the [js/config.js](../src/js/config.js) file. The views section looks like this:

```js
defaultView: 'awesome.js',
views: [
    "default.js",
    "awesome.js",
    "library.js",
    ...
]
```

You can add your own views here or change the order. If the section is missing or empty, the default layer will be loaded.

To set a particular view as the default active view, use the `defaultView` property.

## Set up your data layers

Please refer to the [layer reference document](./layer-reference.md).

The `js/views/layers/` directory contains some sample layers with instructions on how to acquire/preprocess the data and the example layer at [js/views/layers/example.js](../src/js/views/layers/example.js) shows the different options.

You will probably want to center your map and define the default/min/max zoom levels.

Make sure to fill in the credits section for each layer. This corresponds to the credits section as explained in the [js/views/layers/example.js](../src/js/views/layers/example.js) file although a base layer only has one credit object. Some free to use online tile servers are contained in the config file.

# Usage

Once installed and configured, open the [index.html](../src/index.html) page in a browser. Do this by navigating to

```bash
http://localhost/BigX
```

if you've copied BigX to your webserver's root directory or to

```bash
http://localhost:8080/BigX
```

if you're running BigX using docker.

# Development

You can easily extend BigX yourself.
The enviromnent in which it is developed uses node.js.
The following scripts have been defined in the project's [package.json](../package.json) and are executed using yarn:

* ```bash
  yarn lint
  ```

  Runs [ESLint](https://eslint.org/) to detect and fix problems in the code. When using CI as
  configured in [.gitlab-ci.yml](../.gitlab-ci.yml), the pipeline will fail if the linting stage
  fails to ensure only clean code is deployed.<br /><br />

* ```bash
  yarn jshint
  ```

  Runs [JSHint](https://jshint.com/) to make sure the code is to a good standard.<br /><br />

* ```bash
  yarn test
  ```

  Runs all unit tests for the project. Like the linting stage, only code where all tests pass
  will be deployed.<br /><br />

* ```bash
  yarn test1 [your.test.js]
  ```

  Runs a single unit test in isolation.<br /><br />

* ```bash
  yarn badges
  ```

  Creates badges for test coverage of the most recent test run.<br /><br />

* ```bash
  yarn stryker
  ```

  Runs mutation testing for the project. This can take a long time (1h+).<br /><br />

* ```bash
  yarn start
  ```

  Runs BigX using [Electron](https://www.electronjs.org/), i.e. in an isolated stand-alone
  environment instead of a webserver/browser setup.<br /><br />
# BigX documentation

![](src/img/bigx.png)

Welcome to the BigX docs for advanced users.

For "end users" (i.e. people who just browse an already curated map), the [manual](../src/manual.md) should be sufficient. These documents aim at advanced users who either set up BigX or use it for exploratory data analysis.

If you're overwhelmed by the terminology, I suggest you check out [this excellent overview](https://mapschool.io/) to get started.

* [**Setup**](./setup.md)\
  Instructions for installing/running BigX.
* **Usage**\
  Documentation for people who set up BigX by adding new data.
  - [**Layer reference**](./layer-reference.md)\
    A comprehensive guide to curating new layers for BigX.
  - **Views reference**\
    Coming soon
# Layer reference

This document describes the anatomy of a layer.

* [Layer properties](#layer-properties)
* [`mainData`](#maindata)
  - [`file`](#file)
  - [`type`](#type)
  - [*`latField`/`lonField`*](#latfieldlonfield)
  - [*`dataFields`*](#datafields)
  - [*`filters`*](#filters)
  - [*`rawcontent`*](#rawcontent)
* [*`extraData`*](#extradata)
  - [`mainID`/`dataID`](#mainiddataid)
  - [*`matchingMode`*](#matchingmode)
  - [*`caseSensitive`*](#casesensitive)
  - [*`groupname`*](#groupname)
* [`options`](#options)
  - [*`credits`*](#credits)
  - [*`whitelist`/`blacklist`*](#whitelistblacklist)
  - [*`minZoom`/`maxZoom`*](#minzoommaxzoom)
  - [*`featureTypes`*](#featuretypes)
  - [*`allowMultipleFeatureTypes`*](#allowmultiplefeaturetypes)
  - [*`hideTypes`*](#hidetypes)
  - [*`hover`*](#hover)
  - [*`tooltip`*](#tooltip)
  - [*`dynamicDataLoading`*](#dynamicdataloading)
* [Layer type-specific options](#layertype-specific-options)
  - [*`point`*](#point)
  - [*`line`*](#line)
  - [*`area`*](#area)
  - [*`overlay`*](#overlay)
  - [*`heat`*](#heat)
  - [*`migration`*](#migration)
  - [*`raster`*](#raster)
  - [*`timeline`*](#timeline)
* [General](#general)
  - [Comparison operators](#comparison-operators)

Each layer contains markup that describes what it looks like.
The idea is that layer can be defined without any coding knowledge.
Because each layer is only loaded when it is selected, the layer file is a `.js` file instead of a `.json` file.
This also has the advantage that `.js` files are more forgiving when it comes to syntax and don't require the property names to be in quotes.
This mechanism treats the layer definitions as Javascript objects that are imported into BigX, so it requires them to be exported in their definition file.

The following snippet shows the outline of a layer:

```js
export default {
    name: 'My custom layer',
    type: 'area',
    mainData: {},
    extraData: [],
    options: {
        credits: [],
        whitelist: {},
        blacklist: {},
        minZoom: 10,
        maxZoom: 12,
        featureTypes: {},
        allowMultipleFeatureTypes: true,
        hideTypes: [ 'bus', 'giant', undefined ],
        hover: {},
        tooltip: [],
        dynamicDataLoading: false,
        
        point: {},
        line: {},
        area: {},
        overlay: {},
        heat: {},
        migration: {},
        raster: {},
        timeline: {}
    }
}
```
The following sections will explain each part of the layer definition where *italic* properties are optional

## Layer properties

### `name`

Each layer must provide a human-readable layer name.

```js
name: 'My custom layer'
```

This name is used in the map key and the layer selection.

![Layer name appearing in the map key](img/name-key.png)

![Layer name appearing in the layer selection](img/name-layers.png)

Furthermore the names need to be unique within each category and layer type so there can't be two layers of the same name and type in the same category.

### `type`

The type defines what kind of layer it is.

```js
type: "area"
```

The options are
* point
* line
* area
* overlay
* heat
* migration
* raster
* timeline

![Layer type appearing in the layer selection](img/type.png)

Please refer to the individual layer options for more details.

## `mainData`

This section defines the "main" data you want to visualise.
It contains the individual features to be shown on the map.

```js
mainData: {
    file: 'data/csv/my_data.csv',
    type: 'point',
    latField: 'lat',
    lonField: 'lon',
    dataFields: ['layer_name', {'some_data': 'newfieldname'}, 'more_data' ],
    filters: [
        { field: 'Movement', largerThan: 0.001 }
    ],
    rawcontent: {}
}
```

Originally this was called "geoData" as it is usually geographic data but sometimes the "main" data is not geographic but contains a reference to additional data which contains geographic locations.

### `file`

The name of the file including the file extension.

```js
file: 'data/csv/my_data.csv'
```

If you use a local file, it should be in the data directory and this
path is relative to the BigX source directory (`/src`).
You can also use files from online sources, but with one caveat:
- the server from which you load the file needs to have CORS enabled or
- you need to use a CORS proxy and prepend it to your URL like so:
  `https://my-cors-proxy.org/https://example.com/data/route.gpx`

### *`type`*

```js
type: 'point'
```

Using this optional parameter you can choose to force a default feature type if the data does not explicitly contain types. The options are:
- point: parse file as points if possible
- line: merge all points into one line
- area: merge all points into one area

This allows you to reinterpret your existing data depending on how you would like to present it.

### *`latField`/`lonField`*

```js
latField: 'lat',
lonField: 'lon'
```

These fields let you specify the name of the columns containing latitude/longitude if your file is a `.csv` file. This will be used to generate a GeoJSON layer "on the fly", which is the default data format in Leaflet. Currently, only `.csv` files that contain points are supported.

These fields are optional and will only be used if a feature does not have coordinates yet.

### *`dataFields`*

```js
dataFields: ['layer_name', {'some_data': 'newfieldname'}, 'more_data' ]
```

You can specify the names/headings of the data fields/columns in your file that should be read. If undefined, all data found in the file will be attached to the feature, if empty, none will be attached. This allows loading only part of a file into BigX, preserving memory by leaving out data that is not needed without having to modify the file.

You can also use key/value pairs to define a new name for a field.
The key is the identifier in the data file and the value is the new name of the field in your feature.

### *`filters`*

```js
filters: [
    { field: 'Age', smallerOrEqual: 12 }
]
```

Similarly to the `dataFields` option, filters allow you to only use part of a file. While dataFields filter out unwanted columns or properties, filters allow you to leave out features that don't conform to a rule.
`field` refers either to a column in a `.csv` file or a property of a feature.

All [comparison operators](#comparison-operators) can be used.

### *`rawcontent`*

This option lets you scrape content live off the internet.
It is an advanced feature and currently only documented in the [/src/js/views/layers/example.js](../src/js/views/layers/example.js) file.

## *`extraData`*

You can choose to combine your single main dataset with one or more additional datasets.

This is particularly useful if

* your main data does not contain geographic information or
* your main data contains only geographic information but no data properties or
* the data you want to visualise exists in several different files but combining them is not feasible because duplicate information would inflate the file size and make loading slow

```js
extraData: [
        {
            file: 'data/path/to/file.csv',
            mainID: 'ID',
            dataID: 'feature_id',
            matchingMode: 'equal',
            caseSensitive: false,
            groupname: 'propertyGroup'
        }
    ],
```

Just like for `mainData`, the following fields apply:

* `file`
* *`type`*
* *`latField`/`lonField`*
* *`dataFields`*
* *`filters`*


There are several additional fields that can be specified per extra data file.
### `mainID`/`dataID`

To combine multiple datasts into one, you need an identifier for your features which is used to map/join your datasets.

If there's only one row for each feature from `mainData` in this `extraData`, all the column names can be mapped directly into the feature's properties:

```json
{
    "type": "FeatureCollection", "crs": {...},
    "features": [
        {
            "type": "Feature",
            "properties": {"id": "1"},
            "geometry": {...}
        },
        {
            "type": "Feature",
            "properties": {"id": "2"},
            "geometry": {...}
        },
        ...
    ]
}
```

Now we attach the following CSV:

```csv
identifier,name
1,Foo
2,Bar
...
```

The config snipped to do this looks like this:

```js
mainID: 'id',
dataID: 'identifier'
```

The resulting features as loaded into BigX will look like this:

```json
{
    "type": "FeatureCollection", "crs": {...},
    "features": [
        {
            "type": "Feature",
            "properties": {"id": "1", "name": "Foo"},
            "geometry": {...}
        },
        {
            "type": "Feature",
            "properties": {"id": "2", "name": "Bar"},
            "geometry": {...}
        },
        ...
    ]
}
```

This works out of the box for CSV files where the identifier is unique. If your `extraData` provides multiple data points for each feature, see the [*groupname*](#groupname) option for how to set this up.

### `matchingMode`

Data can be mapped based on a multitude of conditions.
All [comparison operators](#comparison-operators) can be used.
The default is `equal`.

```js
matchingMode: 'contains'
```

### *`caseSensitive`*

Matching of `mainID` to `dataID` is case-sensitive by default. This can be disabled here.

```js
caseSensitive: false
```

### *`groupName`*

This optional field lets you group properties for features.

This is necessary if there is more than one row in the `extraData` that matches the feature's identifier, the values in each row need to be kept together to avoid losing the context. In that case, a groupname is used to access that data later.

If for example you would like to visualise data for european electoral regions, you'd start with the GeoJSON for these regions *(Contains OS data © Crown copyright and database right 2018)*:

```json
{
    "type": "FeatureCollection",
    "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
    "features": [
        {
            "type": "Feature",
            "properties": {"EER13CD": "E15000001", "EER13NM": "North East"},
            "geometry": {"type": "Point", "coordinates": [-1.647849865331362, 55.50051475212853]}
        },
        {
            "type": "Feature",
            "properties": {"EER13CD": "E15000002", "EER13NM": "North West"},
            "geometry": {"type": "Point", "coordinates": [-3.056281529493948, 54.1808916482116]}
        },
        ...
    ]
}
```

There is not extra data in these features apart from their code and name.
However, other datasets exist that don't contain the locations but do use the code to refer to the features, e.g. the traffic statistics by the Department for Transport *(Source: Office for National Statistics licensed under the Open Government Licence v.3.0. Contains OS data © Crown copyright and database right 2018. © 2019 Department for Transport)*:

```csv
ons_code,year,all_motor_vehicles
E15000001,1993,10262258321
E15000001,1994,10448184265
E15000001,1995,10624571491
...
```

The CSV contains multiple rows with that identifier; this is why we've chosen the GeoJSON file as our main file, so that there will only be one feature for each european electoral region.
By linking the identifiers using the `mainID` and `dataID` fields, you can tell BigX that the field `EER13CD` in the GeoJSON file corresponds to the `ons_code` field in the CSV file but with every new line from the CSV that contains an already existing identifier, the data value(s) would be overwritten.
This can be avoided by providing a `groupname` to tell BigX that we're expecting multiple data points per feature and want to keep them all.


The layer config snippet looks like this:

```js
mainID: 'EER13CD',
dataID: 'ons_code',
groupname: 'traffic'
```

The resulting features will be as follows:

```json
{
    "type": "FeatureCollection", "crs": {...},
    "features": [
        {
            "type": "Feature",
            "geometry": {...},
            "properties": {
                "EER13CD": "E15000001", "EER13NM": "North East",
                "traffic": [
                    { "year": 1993, "all_motor_vehicles": 10262258321 },
                    { "year": 1994, "all_motor_vehicles": 10448184265 },
                    { "year": 1995, "all_motor_vehicles": 10624571491 },
                    ...
                ]
            }
        },
        ...
    ]
}
```

## `options`

This is the configuration part of a layer. Apart from layer type-specific configuration, all fields are optional and if not given will either not apply or use default values.

### *`credits`*

Credits are technically optional but your data may require you to credit its provider.
You can add multiple credits; each one will be visible while the layer is selected.

Within each credit, all fields are optional but you need at least either a text or a link.
If no copyright symbol is given, it will be added automatically.

```js
credits: [
    {
        text: 'Super Data Corporation',
        link: 'https://example.com',
        date: '2012-2019'
    },
    {
        text: 'Other Co.',
        link: 'https://example.com'
    },
    {
        link: 'https://data.gov.uk'
    },
    {
        text: '&copy; Robert Smith, all rights reserved'
    }
]
```

The configuration above would look like this:

![Generated credits](img/creditsconfig.png)

### *`whitelist`/`blacklist`*

Add filters here as black- or whitelists. Note that whitelists have priority
over blacklists and only one can be used at the same time.

For example to only load features where the `LEGEND` property is either `Viewpoint` or `Nature Reserve` write:

```js
whitelist: {
    'LEGEND': ['Viewpoint', 'Nature Reserve']
}
```

To load all features except those where the `TYPE` property is `Trail` write:

```js
blacklist: {
    'TYPE': 'Trail'
}
```

Values for can be given as string or array as explained in more detail in [comparison operators](#comparison-operators).
Filters are executed *before* any extra data is loaded. This means,
that they can only filter by properties, which are present in the original
dataset. To exclude data using additional properties use type filters.

### *`minZoom`/`maxZoom`*

Define the zoom level range, in which this layer should be displayed.
This is useful for large layers, that don't make sense to display when
the map is all zoomed out because there are too many of them, or layers
that show an overview, such as a heatmap, that doesn't make much sense
when the map is all zoomed in.

Both are optional and if not specified the layer will be displayed at any zoom level.

```js
minZoom: 10,
maxZoom: 12
```

### *`featureTypes`*

Here you can define (mutually exclusive) feature types for your data.
The labelling is done by matching a property against a value and assigning the `featureType` property with the `featureType` key as an extra property for your feature if it matched the rule.

Types can then be used to assign different markers, styles or tooltips for features of each type.
Matching is done in order of definition and is greedy by default.
Any feature with a type not defined here will have the type `undefined`.
All [comparison operators](#comparison-operators) can be used.

Here are some examples for feature types:
```js
featureTypes: {
    'standard': {key: 'LEGEND', equal: 'Standard Gauge'},
    'rural': {key: 'location', notEqual: ['city','urban']},
    'male': {key: 'title', startsWith: 'Mr.'},
    'academic': {key: 'title', endsWith: ['MD', 'PhD']},
    'vehicle': {key: 'label', contains: ['car', 'van', 'vehic']},
    'tiny': {key: 'radius', smallerThan: 5},
    'small': {key: 'radius', smallerOrEqual: 5},
    'large': {key: 'diameter', largerThan: 30, smallerOrEqual: 30},
    'giant': {key: 'diameter', largerOrEqual: 30}
}
```

### *`allowMultipleFeatureTypes`*

The default behaviour for assigning a feature type is to use use the first matching type. To make the matching non-greedy, you can allow multiple feature types, which will turn
the `featureType` property into an array containing all matched feature types or an
empty array if non were found.

```js
allowMultipleFeatureTypes: true
```

### *`hideTypes`*

Once feature types have been defined, they can be excluded from being displayed:

```js
hideTypes: [ 'bus', 'giant', undefined ]
```

### *`hover`*

This defines the behaviour when you hover over a feature. The factors are applied to the feature's initial style as defined in the line/area options and can be positive or negative.

* Brightness is increased/decreased based on the percentage value given.
* Opacity is capped at 0 and 1 respectively and markers are optional.

Currently, start and end markers are only supported for lines, not areas.
All values that are not set will fall back to the default hover behaviour shown below:

```js
hover: {
    lineWeightFactor: 2,
    lineOpacityFactor: 3,
    lineBrightnessPercentage: 20,
    fillOpacityFactor: 3,
    fillBrightnessPercentage: 20,
    shadow: true,
    startMarker: ["#f66", "S", "#000"],
    endMarker: ["#6c0", "F", "#000"]
}
```

Example hover images:

![Hover effect for a line feature](img/hover-line.png)\
*(© OpenStreetMap contributors. Contains Ordnance Survey, Office of National Statistics, National Records Scotland and LPS Intellectual Property data © Crown copyright and database right 2016. Data licensed under the terms of the Open Government Licence. Ordnance Survey data covered by OS OpenData Licence.)*

![Hover effect for an area feature](img/hover-area.png)\
*(© Map tiles by Stamen Design, under CC BY 3.0. Data by OpenStreetMap, under ODbL. Contains OS data © Crown copyright and database right 2018)*

![Start and end markers for line hover effect](img/hover-markers.png)\
*(© OpenStreetMap contributors, © 2003-2020 Walking Englishman)*

### *`tooltip`*

While hover effects are applied to all features, tooltips can depend on properties of individual features. A tooltip is a concatenation of information and hence represented as an array.

Apart from simple strings and HTML inside a string, the following types are supported:

- `["var","property name"]`\
  Prints the property of the feature.
- `["eval", "expression"]`\
  Appends the result of the expression. You can access everything in the scope of the feature.
- `["nvar", "extraProps", 1]`\
  Prints a nested variable up to a depth (optional). The resulting HTML uses unordered lists.
- `["if", {key: "prop", operator: "value"}, "Available", "Not available"]`\
  Parses a condition and prints the first argument if it's true and the second argument if it's false. See [comparison operators](#comparison-operators) for supported operators.

Since all elements are evaluated recursively, arrays can be used, as long as they don't use one of the supported keywords. The order of elements in an array is retained.

```js
tooltip: [
    "<b>", ["var", "FERRY_FROM"], " - ", ["var", "FERRY_TO"], "</b><br />", 
    "Type: ", ["var", "FERRY_TYPE"], "<br />",
    ["eval", "2+(feature.properties.trains!=undefined?" +           
                "Object.keys(feature.properties.trains).length:0)"],
    ["nvar", "destinations", 3],
    ["if", {key: "available", equal: "1"}, "Available", "Not available"],
    ["if", {key: "open", equal: true}, ["Open", ["var", "open"]], "Closed"]
]
```

![Hover effect for a line feature](img/tooltip.png)\
*(© OpenStreetMap contributors, © 2009-2012 Atos and RSP)*

### *`dynamicDataLoading`*

This defines whether the features in a layer should be loaded dynamically
whenever the map is oved or zoomed. Depending on the data, this can affect
the loading speed and responsiveness of a layer. The default is "true".

```js
dynamicDataLoading: false,
```

## Layer type-specific options


### *`point`*


### *`line`*


### *`area`*


### *`overlay`*


### *`heat`*


### *`migration`*


### *`raster`*


### *`timeline`*

## General

### Comparison operators

These operators are used to compare two values. Typically a comparison looks something like this:

```js
{ field: "Temperature", largerThan: 10 }
```

This will compare all entities for which the `Temperature` field is > 10.

The right-hand side of a comparison can either be a single value or a set of values:

```js
{ field: "Name", startsWith: ['Ba', 'Wi', 'De'] }
```

Operators for numeric values can be combined:

```js
{ field: "Temperature", largerThan: 10, smallerOrEqual: 20 }
```

The available operators are:

- `equal`: a number or string or an array, where any item is matched (OR)
- `notEqual`: a number or string or an array where none of the items must be matched (AND)
- `startsWith`: a string that starts with any of the given values (OR)
- `endsWith`: a string that ends with any of the given values (OR)
- `contains`: a string or array, at least one of which is contained in the value (OR)
- `smallerThan`: a single number, representing a maximum exclusive value
- `smallerOrEqual`: a single number, representing a maximum inclusive value
- `largerThan`: a single number, representing a minimum exclusive value
- `largerOrEqual`: a single number, representing a minimum inclusive value
---
title: 'BigX: A geographical dataset visualisation tool'
tags:
  - JavaScript
  - Leaflet.js
  - data visualisation
authors:
  - name: Stefanie Wiegand
    orcid: 0000-0003-3030-9853
    affiliation: 1
affiliations:
 - name:  Bioengineering Sciences Research Group, School of Engineering, Faculty of Engineering and Physical Sciences, University of Southampton 
   index: 1
date: 18 June 2020
bibliography: paper.bib
---

# Statement of Need

Geographic open data is widely available today not only to researchers but to the general public. Globally operating organisations such as [OpenStreetMap](https://www.openstreetmap.org/), and more localised portals such as the [Open Geography Portal](https://geoportal.statistics.gov.uk/) and [Open Development Cambodia](https://opendevelopmentcambodia.net/) provide geographic data in a multitude of formats to download and analyse.

There are a variety of tools available to suit users with different objectives, requirements and technical abilities, all of which have different advantages and disadvantages regarding their ease of use, interoperability, extensibility, licence, cost, and other factors.

Technically inclined users can harness programming languages such as [Python](https://www.python.org/) or [R](https://www.r-project.org/) to process and analyse their data and plot it using tools and libraries such as [MatLab](https://www.mathworks.com/products/matlab.html), [Matplotlib](https://matplotlib.org/) or geoplot [@geoplot] and the Javascript-based library Leaflet [@leaflet] has made the visualisation of geographic data more accessible to software developers outside the scientific community.

Closed-source commercial desktop-based tools such as [ArcGIS](https://www.arcgis.com/index.html) are popular, but open-source alternatives such as [QGIS](https://www.qgis.org/) and [uDig](http://udig.refractions.net/) exist. For paying subscribers, online services such as ArcGIS Online, and [Carto](https://carto.com/pricing/) are available. While free services such as [QGIS Cloud](https://qgiscloud.com/) exist, private maps are typically premium features, which raises questions of data privacy.
In order to securely share interactive maps with a small number of trusted users online while retaining control over the physical location of the data in accordance with data protection laws, users can publish their own maps using tools such as qgis2web [@qgis2web] to convert their maps to static web-based maps although "many elements cannot be reproduced" [@qgis2web] and a new static version has to be produced after each change to the original map.

Statistical data can readily be obtained from places, such as [data.gov.uk](https://data.gov.uk/), but it often relates to geographic locations using identifiers. While these are typically defined in publicly available datasets, they often come in formats such as CSV and need to be linked to the geographic data before the data can be visualised.

Without the technical skills to perform preliminary data processing, cleaning, and mapping, or the money to commission such work, certain groups are left out:

* Researchers from non-technical fields such as social sciences miss out on the benefits data visualisation can give them or lose time on manually processing large datasets.
* Users and organisations with an interest in analysing their own data but without the budget to pay for professional software or services often have no other option than handing over their data to online services, but
* in places where internet connectivity is bad (such as Lower- and Middle Income Countries), transferring large datasets is not always possible and
* sensitive data, as used in medical-related fields, often comes with restrictions on how and where it can be stored, leaving users at the mercy of providers of such services and their data protection policies.

# Summary

BigX addresses these issues as follows:

* By using an open-source license, it is freely available to anyone.
* All data processing is done client-side, which makes it easily configurable to work offline or semi-offline (using an online base-map but local data). This means that:
  - no data ever leaves the machine, so it cannot be intercepted or otherwise leaked,
  - because data does not need to be transferred, response times are much faster (subject to the machine it runs on).
* BigX can be deployed using simple access control measures such as .htaccess files and allows making updates to maps in a simple text editor, thus supporting environments where GUI access may not be present.
* The markup-like definition of layers does not require the user to have programming skills --- a basic understanding of the data is sufficient to generate a data layer. \autoref{fig:map} shows an example for a layer markup and the result. There are various layer types and options available to suit different types of data, all of which are explained in the [manual](https://gitlab.com/geekysquirrel/bigx/-/blob/dev/src/manual.md).

![The layer markup (left) and the generated layer on the map (right). Contains Ordnance Survey, Office of National Statistics, National Records Scotland and LPS Intellectual Property data © Crown copyright and database right 2016. Data licensed under the terms of the Open Government Licence. Ordnance Survey data covered by OS OpenData Licence. Contains HM Land Registry data © Crown copyright and database right 2018. This data is licensed under the Open Government Licence v3.0.\label{fig:map}](code_and_map.png)

So far, BigX has been used to map population movements in China during the COVID-19 pandemic [@popmap] and is being used in current research to map travel distances for prosthetics patients in Cambodia as part of a wider project [@lmic].

# Acknowledgements

While BigX was not funded directly, the author would like to thank Dr. Alex Dickinson and Prof. Andy Tatem for the opportunity to apply it to real-world research.

# References
