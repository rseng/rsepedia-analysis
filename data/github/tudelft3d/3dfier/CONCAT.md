[![GitHub license](https://img.shields.io/github/license/tudelft3d/3dfier)](https://github.com/tudelft3d/3dfier/blob/master/LICENSE)
[![docs](https://img.shields.io/badge/docs-http://tudelft3d.github.io/3dfier-brightgreen)](http://tudelft3d.github.io/3dfier)
[![GitHub issues](https://img.shields.io/github/issues/tudelft3d/3dfier)](https://github.com/tudelft3d/3dfier/issues)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02866/status.svg)](https://doi.org/10.21105/joss.02866)


## 3dfier
<img src="docs/images/3dfierLogo.png" width="300">

Takes 2D GIS datasets (e.g. topographical datasets) and "3dfies" them (as in "making them three-dimensional") by lifting every polygon to 3D.
The elevation is obtained from a point cloud (we support LAS/LAZ at this moment), and the semantics of every polygon is used to perform the lifting.
That is, water polygons are extruded to horizontal polygons, buildings to LOD1 blocks, roads as smooth surfaces, etc.
Every polygon is triangulated (constrained Delaunay triangulation) and the lifted polygons are "stitched" together so that one digital surface model (DSM) is constructed.
Our aim is to obtain one DSM that is error-free, i.e. no intersecting triangles, no holes (the surface is watertight), where buildings are integrated in the surface, etc.
This surface will then be used as input in simulation software for instance.

![](docs/images/leiden3dfier.png)

<a href="https://vimeo.com/181421237">This video</a> illustrates the process and what 3dfier is about.

The lifting options can be configured in the [YAML](https://yaml.org/) file provided, an example is provided in `/resources/config_files/myconfig.yml`.
Any 2D input (which should be a planar partition) can be used as input, and each class must be mapped to one of the following:

  1. Building
  1. Terrain
  1. Road
  1. Water
  1. Forest
  1. Bridge
  1. Separation (used for walls and fences)

It is possible to define new classes, although that would require a bit of programming.

Output is in the following formats: OBJ, CityGML, CityJSON, CSV (for buildings only, i.e. their ID and height (ground+roof) are output in a tabular format), PostGIS, and STL.
The ID of each polygon is preserved, and there is a 1-to-1 mapping between the input and the output. 

If you use it, feedback is very much appreciated.


## Documentation
The [3dfier documentation](http://tudelft3d.github.io/3dfier) has extensive information on the installation, usage, and how 3dfier works.


## If you use 3dfier in a scientific context, please cite this article:

Ledoux H, Biljecki F, Dukai B, Kumar K, Peters R, Stoter J, and Commandeur T (2021). 3dfier: automatic reconstruction of 3D city models. *Journal of Open Source Software*, 6(57), 2866. 

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02866/status.svg)](https://doi.org/10.21105/joss.02866)

```
@article{3dfier,
  author = {Ledoux, Hugo and Biljecki, Filip and Dukai, Balázs and Kumar, Kavisha and Peters, Ravi and Stoter, Jantien and Commandeur, Tom},
  doi = {10.21105/joss.02866},
  journal = {Journal of Open Source Software},
  number = {57},
  pages = {2866},
  title = {3dfier: automatic reconstruction of 3D city models},
  volume = {6},
  year = {2021}
}
```

## LAS/LAZ Pointcloud

We expect the LAS/LAZ to be classified according to the ASPRS Standard LIDAR Point Classes v1.4 (Table 4.9 of this [PDF](http://www.asprs.org/wp-content/uploads/2010/12/LAS_1-4_R6.pdf)), and at a minimum these should be defined:

  - 0-1: Created, never classified and/or unclassified
  - 2: Ground
  - 3-5: Vegetation

If the vegetation is not classified or not filtered out, then buildings might be taller and there might be artefacts in the terrain.

## Binary releases for Windows and Mac OS X

In order to make easy use of 3dfier we created pre-build binaries which can be downloaded from the [releases page](https://github.com/tudelft3d/3dfier/releases). 

Download the latest release and unzip the archive in a easy to find location, not in the download folder of your browser. 

To be able to quickly test 3dfier one can download the [example dataset](https://github.com/tudelft3d/3dfier/releases/tag/example_data) and unzip the archive in the folder of 3dfier. 

## Test data

In the folder `example_data` (download [example dataset](https://github.com/tudelft3d/3dfier/releases/tag/example_data)) there is a small part of the [BGT datasets](http://www.kadaster.nl/web/Themas/Registraties/BGT.htm) (2D 1:1k topographic datasets of the Netherlands), and a part of the [AHN3 LIDAR dataset](https://www.pdok.nl/nl/ahn3-downloads) that can be used for testing. 
The resulting model (in OBJ) can be found in `example_data/output/test_area.obj`

Further, there is an [open data website](https://3d.bk.tudelft.nl/opendata/3dfier/) that contains 3D models of a few Dutch cities, generated with 3dfier.

## Validate config file
The configuration is stored in [YAML format](http://docs.ansible.com/ansible/latest/YAMLSyntax.html) and needs to be valid for the parser to read the file. 
Config files can be schema validated using [YAML Lint](http://www.yamllint.com)

## Run 3dfier:
**Windows** 
Open a command line (click start and type `command` or `cmd`). Using the command line browse to the folder where you extracted the example files and run:
`3dfier myconfig.yml -o output.ext`

**Mac OS X and Linux**
Open a console. Using the console browse to the folder where you extracted the example files and run:
`$ ./3dfier myconfig.yml --OBJ output.obj`

**Docker**

3dfier offers a alpine base image which tries to give you as much freedom for your vector data source as possible. Vector data is read by GDAL/OGR.

To run 3dfier over Docker simply execute:

    $ docker run --rm --name 3dfier -v <local path where your files are>:/data tudelft3d/3dfier:<tag> 3dfier <name of config file> <... 3dfier parameters>

All your input data needs to be in `<local path where your files are>` and in the config file you need to reference your input data relative to `<local path where your files are>`. To achieve this either move your data and config into `<local path where your files are>` (and subdirectories), or set `<local path where your files are>` to the lowest common ancestor that contains all the data and config files you need.

**Keep in mind that `<local path where your files are>` need to be writable by any user, otherwise your output won't be saved.**

For instance to run it on the example data set (on Linux):

    $ cd 3dfier/example_data
    $ docker run --rm -it -v 3dfier/example_data:/data tudelft3d/3dfier:latest 3dfier testarea_config_unix.yml --OBJ test.obj

There is also a [tutorial](https://github.com/tudelft3d/3dfier/wiki/General-3dfier-tutorial-to-generate-LOD1-models) on how to generate a 3D model with 3dfier.


## Prepare BGT data
For preparing BGT data as input for 3dfier look at [resources/BGT_prepare/ReadMe.md](https://github.com/tudelft3d/3dfier/blob/master/resources/BGT_prepare/ReadMe.md)




# Contributing

+ The GitHub Issues section is the primary communication channel, please post issues, questions, feature requests there.

+ The repository follows the [Gitflow branching model](http://nvie.com/posts/a-successful-git-branching-model/).

+ To contribute **documentation** consider extending, fixing or adjusting the `/docs` and submit a pull request.

+ To contribute **code**, please submit a pull request.
---
title: '3dfier: automatic reconstruction of 3D city models'
tags:
  - GIS
  - 3D city modelling
authors:
  - name: Hugo Ledoux^[Corresponding author]
    orcid: 0000-0002-1251-8654
    affiliation: 1 
  - name: Filip Biljecki
    orcid: 0000-0002-6229-7749
    affiliation: 2
  - name: Balázs Dukai
    orcid: 0000-0003-0193-5863
    affiliation: 1
  - name: Kavisha Kumar
    affiliation: 1
  - name: Ravi Peters
    affiliation: 1
  - name: Jantien Stoter    
    affiliation: 1
  - name: Tom Commandeur
    affiliation: 1
affiliations:
  - name: Delft University of Technology, the Netherlands
    index: 1
  - name: National University of Singapore, Singapore
    index: 2
date: 16 October 2020
bibliography: paper.bib
---

# Summary

Three-dimensional city models are essential to assess the impact that environmental factors will have on citizens, because they are the input to several simulation and prediction software.
Examples of such environmental factors are noise [@Stoter08], wind [@GarciaSanchez14], air pollution [@Ujang13], and temperature [@Lee13; @Hsieh11].

However, those 3D models, which typically contain buildings and other man-made objects such as roads, overpasses, bridges, and trees, are in practice complex to obtain, and it is very time-consuming and tedious to reconstruct them manually.

The software *3dfier* addresses this issue by automating the 3D reconstruction process.
It takes 2D geographical datasets (e.g., topographic datasets) that consist of polygons and "3dfies" them (as in "making them three-dimensional"). 
The elevation is obtained from an aerial point cloud dataset, and the semantics of the polygons is used to perform the lifting to the third dimension, so that it is realistic.
The resulting 3D dataset is semantically decomposed/labelled based on the input polygons, and together they form one(many) surface(s) that aim(s) to be error-free: no self-intersections, no gaps, etc.
Several output formats are supported (including the international standards), and the 3D city models are optimised for use in different software. 


# Statement of need

The 3D city models needed as input in environmental simulations have specific requirements that go beyond the typical 3D models used for visualisation: they require semantic information (i.e., an object, modelled with one or more surfaces, "knows" what it is, for instance a window or a roof surface) and they should be free of geometric errors.
It is known that practitioners and researchers wanting to perform some simulations or analysis can spend a significant part of their time constructing and repairing the input 3D models; @McKenney98 estimates this to as much as 70\% of their time.
Furthermore, the formats required by the different software and/or the agencies (for instance [the international standard CityGML](https://www.ogc.org/standards/citygml)), are complex to generate [@Ledoux19].

The software *3dfier* automates the reconstruction step, it enriches the data with semantics, and it supports several output formats (used in different fields).

It builds upon previous work done for reconstructing the whole country of the Netherlands (with 10M+ buildings) [@OudeElberink13], and provides the following improvements: use of recent and maintained libraries (e.g., [CGAL](https://www.cgal.org/), [GDAL](https://gdal.org/) and [Boost](https://www.boost.org/)), a clear open-source license, recent formats and international standards are supported, no geometric errors in output.

There exist different commercial software that can perform extrusion (e.g., [Safe FME](https://www.safe.com/fme/) or [ArcGIS](https://www.arcgis.com/)), but these extrude each classes of objects separately (usually only buildings), without ensuring that adjacent objects should be "stitched" together.
As a result, the resulting 3D dataset is often unsuitable as input for other spatial analysis software.


# Overview of the reconstruction steps

![Overview of 3dfier.\label{fig:overview}](extrusion.png){ width=90% }

As shown in \autoref{fig:overview}, as input we use geographical datasets that are readily available for an area (often as open data too):

  1. 2D polygons representing buildings, lakes, roads, parts, etc. ([OpenStreetMap](https://www.openstreetmap.org) is one option);
  2. elevation points, usually acquired with a laser-scanner and available in [LAS format](https://www.asprs.org/wp-content/uploads/2010/12/LAS_1_4_r13.pdf), or derived from aerial images.

Each of the classes in the input 2D polygons is mapped to a specific class: *Terrain*, *Forest*, *Water*, *Road*, *Building*, *Bridge/Overpass*, and *Separation* (walls and fences).

![1D visualisation of the reconstruction process.\label{fig:steps}](steps.pdf)

The semantics of every input 2D polygon is used to perform the lifting to the third dimension.
For example, water polygons are extruded to horizontal polygons, buildings to prismatic blocks, roads as smooth surfaces, etc. 
Every polygon is triangulated and in a next step the lifted polygons are "stitched" together so that one surface is reconstructed. 
In this step, priority is given to "hard" objects such as roads, i.e., vegetation polygons are moved to be aligned with the road polygons.
The output of the software is one watertight surface with no intersecting triangles and no holes, and features such as buildings and trees can be added or omitted.  
Triangles are grouped and labelled with the class and the attributes that were in the input 2D polygons they decompose.
This surface can be used directly as input in several urban applications, such as simulations of noise or wind.


# Use of the software

The software is a command-line interface (CLI), and uses a configuration file as input (a [YAML file](https://yaml.org/)).
This file allows the user to control the mapping between the input and the *3dfier* classes, to specify which [LAS](https://www.asprs.org/wp-content/uploads/2010/12/LAS_1_4_r13.pdf) classes to use/ignore, to control how the lifting is performed for the different classes, etc.

The software, being modular, is also extensible for other use cases or for use in different countries.
As an example, new topographic classes (for instance trees) could be added by simply creating a new C++ class that inherits from the parent class, and the output for the different formats supported must be added.

Great care was taken to keep the software as efficient as possible and make it suitable for reconstructing very large areas. For instance *3dfier* is the software that enables the Dutch national mapping agency (Kadaster) to create the [3D base registration for the Netherlands](https://www.pdok.nl/3d-basisvoorziening).

![An example of the output of 3dfier, for the city of Leiden in the Netherlands.\label{fig:results}](results.png)


# Acknowlegements

This work was funded by the [Netherlands Kadaster](https://www.kadaster.nl/) and received funding from the European Research Council (ERC) under the European Unions Horizon2020 Research & Innovation Programme (grant agreement no. 677312 UMnD: Urban modelling in higher dimensions).

# References


## To run:

`3dfier testarea_config.yml --OBJ out/myoutput.obj`

if you want CityJSON output, change command to `--CityJSON` and run:

`3dfier testarea_config.yml --CityJSON out/myoutput.json`


## 2D input data: the BGT

The files in the folder `bgt` are a crop of the [BGT datasets](http://www.kadaster.nl/web/Themas/Registraties/BGT.htm) in Delft. There were created with the script in `resources/BGT_prepare/..` 

## LIDAR point cloud input: AHN3

The 2 files in the folder `ahn3` are part of the [AHN3](https://www.pdok.nl/nl/ahn3-downloads) cropped for the same area. The 2 tiles used are `c_37en1.laz` and `c_37en2.laz`.

## Example output 

What you should get when you run 3dfier is in `output/testarea.obj` and `output/testarea.json`.
If you use [MeshLab](http://meshlab.sourceforge.net) to visualise the OBJ file, the colours for each class can be activated in the menu `Render/Color/Per Face`. It looks like that:

![](output/testarea.png)


For CityJSON, use [ninja](https://ninja.cityjson.org).
# 3dfier2sim
                                            
Small utility program to convert the output of 3dfier to a mesh that can be used directly (in most cases! we can't guarantee it) as input in [ANSYS](http://www.ansys.com/) software to perform simulation.

Main operations performed:

  1. the model, which is a surface, is closed to form a volume. This is done by expanding the edges of the dataset by roughly 200m and adding vertices at the elevation zero, and then adding new faces to close the volume.
  2. the surface of the model is automatically repaired 
    - by filling holes in it with new surfaces (which are also triangulated)
    - by flipping the orientation of the surfaces so that the model is a valid 2-manifold
  3. the new surfaces added to close the volume are all triangulated.

## Compilation

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ ./3dfier2sim


## Usage

To convert the test-dataset of 3dfier (in `/example_data/`):

  1. remove the bridges in the yml file
  2. run 3dfier with OBJ output, let us call it `testarea.obj`
  3. it's a good idea to translate the area to (minx, miny) since the Dutch coordinates are rather large and can affect precision, use `/resources/translate2min-obj.py' and create the file `testarea_t.obj`
  4. convert that file to an [OFF file](https://en.wikipedia.org/wiki/OFF_(file_format)), eg `testarea.off` in our case
  5. `3dfier2sim testarea.off` will give a summary of the operations done and output a file `out.off`
  6. voilà, you should be able to import this file without spending too much time (semi-)manually fixing the errors.


When preparing the BGT GML files as input for 3dfier, the BGT_conversion.bat batch script can be used. 
The script requires GDAL >2.0 to perform the filtering and conversion.
OGR is used to only select polygonal geometries and stroke CurvePolygons, it also filters the history of objects by selecting those of which the 'eindregistratie' is not set.

Usage:
- Extract the BGT GML files in a folder
- Copy all GFS files from "BGT gfs files" into the same folder
- Run BGT_conversion.bat (or the commands within if on Linux/Mac)
- Use generated GPKG files as input to 3dfier

Known issues:
* Sqlite doesn't support dashes in column names (IMGeo atrributes) thus for IMGeo output we changed to conversion to GeoPackage
* ogr2ogr fails with message on 'eindregistratie no such attribute'; remove the 'eindregistratie == NULL' from the command. This happens when the GML file does not contain any historical objects with an 'eindregistratie' attribute.
* ogr2ogr fails with 'ERROR 1: Did not get at least 3 values or invalid number of set of coordinates'; in bgt_begroeidterreindeel.gml and bgt_onbegroeidterreindeel.gml replace all '<gml:posList>' strings with '<gml:posList srsDimension="2">'. Not all 'kruinlijn' geometries are writting with the srsDimension argument which ogr does not accept.
---
title: The open-source tool for creation of 3D models
keywords: 3dfier homepage
summary: 
sidebar: 3dfier_sidebar
permalink: index.html
---

3dfier tries to fill the gap for simply creating 3D models. It takes 2D GIS datasets (e.g. topographical datasets) and "3dfies" them (as in "making them three-dimensional"). The elevation is obtained from a point cloud (we support LAS/LAZ at the moment), and the semantics of every polygon is used to perform the lifting. After lifting, elevation gaps between the polygons are removed by "stitching" together the polygons based on rules so that a watertight digital surface model (DSM) fused with 3D objects is constructed. A rule based stucture is used to extrude water as horizontal polygons, create LOD1 blocks for buildings, smooth road surfaces and construct bridges (3D polygonal surface). This software is developed by the [3D Geoinformation group](https://3d.bk.tudelft.nl) at the **Delft University of Technology**.

Our aim is to obtain one model that is error-free, so no intersecting objects, no holes (the surface is watertight) and buildings are integrated in the surface.

{% include imagezoom.html file="leiden3dfier.png" alt="3dfier result of Delft" %}---
title: "Page Not Found"
search: exclude
---  

Sorry, but the page you were trying to view does not exist. Try searching for it or looking at the URL to see if it looks correct.
---
title: Model output flow
keywords: output model flow
sidebar: 3dfier_sidebar
permalink: output_flow.html
---

## Writing a model including database and GDAL flow
To output a model all output formats defined in the command line options are iterated and the corresponding format is created. Below the steps for format decision making is shown including writing to GDAL formats and a PostGIS database. The code for writing to PostGIS is depending on GDAL too. Therefore the code for creating features, geometry and attributes is reused. The difference is where the GDAL formats create a dataset (file) per layer where the PostGIS database creates a single dataset (database) with multiple layers. For all other layers lookup the corresponding flow diagrams in the following sections.

{% include imagezoom.html file="flows/3dfier_writing_model.png" alt="Flow diagram for writing models" %}

## CityJSON
Creating CityJSON output is similar to [OBJ](#obj) with some small differences. CityJSON supports Solid objects and attributes where OBJ does not. For Building objects it makes LoD1 Solid objects and includes floors with `floor: true` configured. All other objects are created as MultiSurface objects.

For buildings it writes the `min height surface` attribute with the floor height and `measuredHeight` attribute with the absolute roof height (roof-floor).
{% include imagezoom.html file="flows/3dfier_writing_model_CityJSON.png" alt="Flow diagram for writing models in CityJSON" %}

## OBJ
OBJ output supports only writing objects as MultiSurface. Buildings are supported in both LoD0 and LoD1. When configured LoD0 produces a building with the roof polygon at floor height. This resembles the foundation of the house and can be used in case there are existing building models.

Using a work-around objects in an OBJ model contain an ID and are semantically labelled by defining a material.
~~~
o b222b6a92-00b5-11e6-b420-2bdcc4ab5d7f
usemtl Terrain
f 2596 2597 2598
f 2598 2599 2600
f 2597 2601 2602
~~~
{% include imagezoom.html file="flows/3dfier_writing_model_OBJ.png" alt="Flow diagram for writing models in OBJ" %}

## CityGML
The output for CityGML can be influenced by the configuration in two ways. Setting `triangulation: true` will create a triangulated MultiSurface Building object instead of a LoD1 Solid. `floor: true` will result in the model to contain the floor polygons in buildings. 

For buildings it writes the `min height surface` attribute with the floor height and `measuredHeight` attribute with the absolute roof height (roof-floor).
{% include imagezoom.html file="flows/3dfier_writing_model_CityGML.png" alt="Flow diagram for writing models in CityGML" %}

## IMGeo
IMGeo is an Application Domain Extension of [CityGML](#citygml). This means the core of the data model is equal to that of CityGML but some classes are extended to describe specific (meta)data. This data format is used by the government in The Netherlands and many of the datasets are in this data format. The flow is comparable to that of CityGML except for two things; 1) LoD0 floor and roof polygons are not stored and 2) IMGeo has specific attributes that are copied from the input to the correct output attributes.
{% include imagezoom.html file="flows/3dfier_writing_model_IMGeo.png" alt="Flow diagram for writing models in IMGeo" %}

## CSV
The model output in CSV format is a fantastic tool for building statistics. The data can be output in three different versions:
1. Single percentile - Write `object id,roof height at percentile,floor height at percentile`
2. Multiple percentile - Write `object id,roof height at percentiles,floor height at percentiles` for floor percentiles; 0, 10, 20, 30, 40, 50 and for roof percentiles 0, 10, 25, 50, 75, 90, 95, 99.
3. All elevation values - Write `object id,elevation value,elevation value,elevation value,...` 

{% include imagezoom.html file="flows/3dfier_writing_model_CSV.png" alt="Flow diagram for writing models in CSV" %}
---
title: Main data flow
keywords: main flow
sidebar: 3dfier_sidebar
permalink: main_flow.html
---

## Guide for reading
The architecture of 3dfier is described mainly using a collection of flow diagrams. All flow diagrams are connected and originate from the general flow. When looking for a specific topic look through the general flow and look up the corresponding name for the next flow diagram containing the information for the topic. The [Threedfy flow]({{site.baseurl}}/threedfy_flow.html) and the [Output flow]({{site.baseurl}}/output_flow.html) are located in separate pages in the Architecture menu. They are split from the general flow to keep the overview.

## Data flow through the program
The flow diagrams use specific shapes and colours to explain the use of a block. The legend explains which block describes what step. The flow diagrams only have a directed flow which means there is only one way to walk through the diagram. It can contain loops and references to other flow diagrams but it will always resume execution after finishing that block. 

Nested flow diagrams come in two different types, a completely separate flow diagram in the case of a large flow or nested within the same diagram. When a flow is nested within the same diagram it is outlined by a grey block.

Some of the flow diagrams have a rounded orange blocks. These blocks are there to explain the combined use of actions contained by the block.

{% include imagezoom.html file="flows/legend.png" alt="Legend for flow diagrams of 3dfier" %}

### General flow
When executing 3dfier it is launched with command line options and a configuration file. The command line options point to the correct [configuration](#configuration) and directs the output to the wanted file format(s). The configuration contains all settings needed for the algorithm execute. The [settings are validated](#configuration) first, the passed values are checked to be of the expected format like string, boolean or numerical format. Settings that are part of the configuration that do not exist as an option in 3dfier are ignored. If the configuration is deemed valid the input files for both polygons and points are opened to check [if they exist](#data-integrity-check). This is done to overcome exceptions later in the process due to misconfiguration.

In case all files exist and are readable the [polygon files are read](#reading-polygons) and all polygons are added to an R-tree for fast intersection queries. From the R-tree bounds the bounding box of the polygons extent is calculated and stored for when the point files are read.

After the polygon files are added the [point files are read](#reading-points). This can be either in LAS or LAZ format. The header of the file is read to acquire the bounding box. The bounding box is intersected with the polygons extent stored in the previous step. The point file is discarded when it doesn't intersect with the polygons.

When all input data is read the magic happens within the [Threedfy](#magic) step. It does the lifting, stitching and creates vertical walls. More about the core of the algorithm can be read in the [Threedfy flow page]({{site.baseurl}}/threedfy_flow.html).

Last is the output of the model. There are various implementations where each is specifically made for a data format. The main flow and the implementation for database and GDAL formats is shown in [Writing model](#writing-model). Other outputs are detailed in the [Output flow page]({{site.baseurl}}/output_flow.html).

{% include imagezoom.html file="flows/3dfier_general_flow.png" alt="General flow diagram 3dfier" %}

### Data types
3dfier makes use of some data types specifically designed for the algorithm. The sections below describe these data types for better readability of the architecture pages. 

#### Map3D
The main object of 3dfier is the Map3D. It is the object that stores all configurations and objects used in the reconstruction process. The software only creates a single Map3D during the reconstruction. After reconstruction writing the model is done by feeding it the Map3D. The Map3D also contains all TopoFeatures.

#### TopoFeature
TopoFeature is short for Topological Feature. This object contains the 2D geometry of a polygon, a list of heights per vertex that is collected from the process that reads 3D points, the NodeColumn's and the final 3D object created during the reconstruction process. The TopoFeature class contains overloads from three subclasses (Flat, Boundary3D and TIN) that contain the seven subclasses that implement the lifting classes.

{% include imagezoom.html file="flows/3dfier_classes.png" alt="Diagram of 3dfier classes" %}

#### NodeColumn
The NodeColumn is a nested vector of integers containing the various heights on a given x,y-location. For each vertex that is processed by the algorithm an entry is created in the NodeColumn. If an object has that vertex the height calculated from the statistics is pushed to the NodeColumn. Finally the NodeColumn is used when making the decision at what height to put the final object.

### Configuration
Reading and validation of the configuration is done using the [YAML-CPP library](https://github.com/jbeder/yaml-cpp). It parses the file into an object and allows for easy access to all values. The configuration file is read and validated first. Value types are tested and options with limited options are verified. If validation fails an error message will be printed and the program execution will stop.

When validation passes, all values are used to set the corresponding parameters in the software. After this initial read the configuration file is not used anymore since all values are now stored in memory.

### Data integrity check
For all input files configured to use a check is done if they exist and readable. The polygon files are opened using an OGR reader and closed directly after. Same is done for the LAS/LAZ files, they are opened with LASlib and directly closed. This prevents failure of the program after partly executing. Also the output filepaths are tested to be existing before execution of the algorithm.

*NOTE: Output folder must exist, it's not created automatically*

### Reading polygons
If all sanity checks are passed the program turns to reading the polygon files. It takes into account the extent setting, only polygons which bounding box intersects with the given extent will be added. 

A file is opened for reading and the id and height_field attributes are loaded. If these attributes do not exist the program will stop execution with an exception. The total number of features in the layer is logged and counters are set for the amount of multipolygons. Next each feature is read, the attributes are stored in a new TopoFeature. A TopoFeature is the internal storage class used by the algorithm to store attributes, 2D and 3D geometries and functionality to lift the objects. The geometry is extracted and depending on the polygon type it is pre-processed. For a multipolygon each separate polygon is added to a new TopoFeature in which the attributes are copied. The id of the feature is altered by adding a trailing dash and counter, e.g. 'id-0, id-1'. Curvepolygons are automatically stoked into small straight line segments. 

*Important: When reading the geometry the boost::geometry library is used to read WKT and both functions boost::geometry::unique and boost::geometry::correct are used for removing duplicate vertices and correcting ring orientation.*

During the creation of a TopoFeature and storing its geometry, all used buffers and vectors are resized to correspond to the amount of points in the polygon outer and inner rings. This preallocates most of the needed memory for storing the elevation information acquired when reading the point clouds.

For a later stage a spatial index is created that stores the bounding box of the geometries with a reference to the TopoFeature. This speeds up filtering the TopoFeatures within range for the point to polygon distance calculations done with the points from the point cloud. The spatial index is a [Boost R-tree](https://www.boost.org/doc/libs/1_72_0/libs/geometry/doc/html/geometry/reference/spatial_indexes/boost__geometry__index__rtree.html). There are two R-trees used, one for Buildings and one for all other objects combined. The reason for this is the two vertex radius settings which can differ, [radius_vertex_elevation]({{site.baseurl}}/output_options.html#radius_vertex_elevation) and [building_radius_vertex_elevation]({{site.baseurl}}/output_options.html#building_radius_vertex_elevation).

The bounding box of all polygons is calculated from the two R-trees combined plus the vertex radius. Since the R-trees contain the bounding box of the object geometries the total bounding box can be larger then the actual area.

{% include imagezoom.html file="flows/3dfier_reading_polygons.png" alt="Flow diagram polygon reading" %}

### Reading points
The elevation information is read from LAS or LAZ files. The file is opened and the header is read. The bounding box from the header is intersected with the bounding box of all polygons. When the file doesn't overlap it is skipped. When the bounds overlap the point count is logged together with thinning setting if configured. For each point read the following things are checked before going into adding a point to a TopoFeature:
- is point to be used according to *i % thinning == 0*
- is classification not within the list of *omit_LAS_classes*
- is point within the bounding box of the polygons

When all checks pass the points is added to the map. The next step is to query both R-trees and return the objects that intersect with this point plus vertex radius in all directions. Each returned object is iterated and the point classification is checked based on the object types. The classification allowed refers to [use_LAS_classes]({{site.baseurl}}/lifting_options.html#use_las_classes) for each object class type configured in the [Lifting options]({{site.baseurl}}/lifting_options.html#). When the point with a certain class needs to be within the polygon, as configured in [use_LAS_classes_within]({{site.baseurl}}/lifting_options.html#use_las_classes_within), an additional boolean is set.

The elevation point is added to the object if its classification is configured for use. Then followed by a check if the point needs to be within the polygon, or just within range. The height is then assigned to each vertex within range of the point. The height is multiplied by 100 (for centimetre accuracy), rounded to an integer (to keep the memory footprint as low as possible) and then added to the elevation vector of the vertex.

Each object class has its own additional set of rules for a point to be used.

{% include imagezoom.html file="flows/3dfier_reading_points.png" alt="Flow diagram point reading" %}

#### Terrain and Forest (TIN)
These classes have the [simplification settings]({{site.baseurl}}/lifting_options.html#simplification) that adds points randomly. Also it stores all elevation points that are within the polygon respecting the [innerbuffer setting]({{site.baseurl}}/lifting_options.html#innerbuffer) in a vector. This vector is used to add additional points inside the polygons for a more precise ground model. 

#### Building and Water (Flat)
For these classes there is a difference to the others. Since these features are modelled as flat surfaces the elevation vector is not maintained per vector but for the TopoFeature as a whole. This lowers the memory footprint and calculation time while maintaining the same accuracy. All points within the polygon or within range of the vertices are stored as described in [radius_vertex_elevation]({{site.baseurl}}/output_options.html#radius_vertex_elevation). For buildings there are two elevation vectors, one for the ground and one for the roof.

#### Road, Separation and Bridge (Boundary3D)
These classes do not have any specific rules for assigning elevation information.

### Magic
The magic of the algorithm happens inside the threedfy [3-D-fy] function. This code is what makes the input data to be merged into a 3D model. Yellow boxes are explained in separate flow diagrams in the [Threedfy flow section]({{site.baseurl}}/threedfy_flow.html).

{% include imagezoom.html file="flows/threedfy.png" alt="Flow diagram for threedfy" %}

### Writing model
Besides creation of the data an important part is to write the created model into a data standard applicable for future use. There are several options available from a 3D modelling, computer graphics or statistics point of view. Most code for writing the model to a file is specifically created for that output format. The flow diagram contains part of all available flows. A more detailed description can be found in [Output flows]({{site.baseurl}}/output_flow).
{% include imagezoom.html file="flows/3dfier_writing_model.png" alt="Flow diagram for writing 3D models" %}---
title: Threedfy flow
keywords: threedfy magic flow
sidebar: 3dfier_sidebar
permalink: threedfy_flow.html
---

## Threedfy flow
All steps done by the Threedfy algorithm are described in the flow diagram below. The yellow blocks refer to nested flow diagrams described in the other sections on this page.

{% include imagezoom.html file="flows/threedfy.png" alt="Flow diagram for threedfying data" %}

## Lifting
Lifting means that all vertices imported from the 2D polygonal geometries are lifted to 3D by statistical analysis of the height points read from the point clouds. All 3D points in a given radius around the vertex are used to extract an appropriate height for the objects vertex.
{% include imagezoom.html file="flows/threedfy_lifting.png" alt="Flow diagram for lifting of polygons" %}

## Stitching
Stitching of polygons is filling holes created by the lifting step in the threedfy algorithm. Topologically connected objects are stitched (like sewing) vertex wise by looking at their type, height an connectedness. Stitching is the most complex part of the algorithm that defines the rules for the final model.

{% include imagezoom.html file="flows/threedfy_stitching.png" alt="Flow diagram for stitching of polygons" %}

## Fix bow ties
A bow tie is a location where due to the stitching the heights of two adjacent objects intersect like a bow tie. An example is two vertices A and B from two adjacent objects. Vertex A from object 1 has height 0 and vertex A from object 2 has height 10 meters. Vertex B in object 1 has a height of 5 meters and vertex B in object 2 has a height of 0 meters. This will result in the creation of a bow tie and therefore not a closed watertight surface. When fixing these issues object types are checked like with the stitching and one of the two objects is stitched to the other in a way the bow tie will be removed.

{% include imagezoom.html file="flows/threedfy_fix_bowties.png" alt="Flow diagram for fixing bow ties" %}

## Create vertical walls
Vertical walls are assigned to objects during the stitching process. When the height difference between the objects is larger then the configured height jump threshold a boolean is set in the object to register it contains vertical walls. The vertical walls are created in this final step. The walls are reconstructed as a collection of triangles in the shape of a fan. Using this the algorithm can make sure all possible heights are topologically connected to the neighbouring objects.

{% include imagezoom.html file="flows/threedfy_vertical_walls.png" alt="Flow diagram for creation of vertical walls" %}---
title: Do's & Don'ts
keywords: do dont good bad practice 
sidebar: 3dfier_sidebar
permalink: do_donts.html
---

## DO
- Do validate you YAML configuration using [www.yamllint.com](http://www.yamllint.com)
- Do make sure the input is topologically correct if you expect a watertight output
- Do [*omit_LAS_classes*] of points you do not use, this improves speed a lot and can help overcome configuration mistakes
- Do make sure all [*input_polygons: datasets: lifting*] class has a corresponding [*lifting_options*] class (e.g. [*input_polygons: datasets: lifting: Building*] & [*lifting_options: Building*])

## DON'T
- Don't expect magic
- Don't combine [*simplification*] and [*simplification_tinsimp*] settings, the latter is always preferred
- Don't use lidar [*thinning*] setting for other then testing to improve speed since this is a simple skip amount while reading points
- Don't forget to configure a value in [*use_LAS_classes*] when using [*use_LAS_classes_within*], it defaults to all classes---
title: Output options
keywords: settings config configuration
sidebar: 3dfier_sidebar
permalink: output_options.html
---

## options
~~~ yaml
building_radius_vertex_elevation: 3.0  # Radius in meters used for point-vertex distance between 3D points and vertices of building polygons, radius_vertex_elevation used when not specified
radius_vertex_elevation: 1.0           # Radius in meters used for point-vertex distance between 3D points and vertices of polygons
threshold_jump_edges: 0.5              # Threshold in meters for stitching adjacent objects, when the height difference is larger then the threshold a vertical wall is created 
threshold_bridge_jump_edges: 0.5       # Threshold in meters for stitching bridges to adjacent objects, if not specified it falls back to threshold_jump_edges
max_angle_curvepolygon: 0.0            # The largest allowed angle along the stroked arc of a curved polygon. Use zero for the default setting. (https://gdal.org/doxygen/ogr__api_8h.html#a87f8bce40c82b3513e36109ea051dff2)
extent: xmin, ymin, xmax, ymax         # Filter the input polygons to this extent
~~~

### radius_vertex_elevation
*Default value: 1.0m.* 
Describes the maximum distance between a vertex and a height points. If the point is within radius of the vertex it will be stored in the vertex height list.

From the image below it is shown what is meant with the radius vertex elevation. The distance to the vertex on the edge of the building in blue is calcuated for each height point and if the point is within the red buffer distance the point is used.
{% include imagezoom.html file="/settings/radius_vertex_elevation.png" alt="" %}

### building_radius_vertex_elevation
*Default value: 3.0m.*
Same as [radius_vertex_elevation](#radius_vertex_elevation) but specific to vertices of buildings. Radius for buildings is larger for finding appropriate ground points due to shading effects of the building itself. In locations where streets are narrow and buildings close to each other it can be hard to find ground points. By setting a larger radius this can be solved. In cases with lots of relief, a larger radius can introduce an incorrect ground height since too many points are taken into account.

### threshold_jump_edges
*Default value: 0.5m.*
When stitching objects their height difference is taken into account. This threshold sets the minimum height difference between to objects to result in a vertical wall to be created. Two objects with a height difference smaller then the threshold will be stitched by adjusting their heights based on set rules.

### threshold_bridge_jump_edges
*Default value: 0.5m.*
Same as [threshold_jump_edges](#threshold_jump_edges) but specific for bridges. In the case of bridges the threshold must sometimes be set to a larger value due to a lower accuracy in height data of the terrain around a bridge. Mostly the case when using Dense Matching point clouds.

### max_angle_curvepolygon
*Default value: 4 degrees.*
Use this setting when using curved polygons as input. Creation of a 3D object from an arc is not possible. Therefore the arcs need to be stroked to lines. The OGR algorithm strokes the arcs identical in both directions of the arc. Because of this the topology is maintained and assured for the stroked lines. Please refer to the OGR API documentation for the function [OGR_G_ApproximateArcAngles()](https://gdal.org/doxygen/ogr__api_8h.html#a87f8bce40c82b3513e36109ea051dff2) that is used to stroke an arc to a line string.

### extent
*Download [YAML]({{site.baseurl}}/assets/configs/extent_green.yml) and [OBJ]({{site.baseurl}}/assets/configs/extent_green.obj)*
As one can see from the examples below, all objects of which its bounding box intersects with the extent is added to the output. This can be used to clip an area without the need to change input configuration.

{% include imagezoom.html file="/settings/settings_extent.png" alt="" %}

There are several areas clipped to show the results. In the following image the used extents are shown and their colors represent that of the results below.
{% include imagezoom.html file="/settings/extents.png" alt="" %}

**Green extent**
*Download [YAML]({{site.baseurl}}/assets/configs/extent_green.yml) and [OBJ]({{site.baseurl}}/assets/configs/extent_green.obj)*
{% include imagezoom.html file="/settings/extent_green.png" alt="" %}
**Blue extent**
*Download [YAML]({{site.baseurl}}/assets/configs/extent_blue.yml) and [OBJ]({{site.baseurl}}/assets/configs/extent_blue.obj)*
{% include imagezoom.html file="/settings/extent_blue.png" alt="" %}
**Purple extent**
*Download [YAML]({{site.baseurl}}/assets/configs/extent_purple.yml) and [OBJ]({{site.baseurl}}/assets/configs/extent_purple.obj)*
{% include imagezoom.html file="/settings/extent_purple.png" alt="" %}
**Red extent**
*Download [YAML]({{site.baseurl}}/assets/configs/extent_red.yml) and [OBJ]({{site.baseurl}}/assets/configs/extent_red.obj)*
{% include imagezoom.html file="/settings/extent_red.png" alt="" %}---
title: Input options
keywords: settings config configuration
sidebar: 3dfier_sidebar
permalink: input_options.html
---

In the input_options we define the settings used for reading the input files for both polygons and point clouds.

If a class is not defined in the YAML the default values will be used. Default values can be found in [resources/config_files/myconfig_DEFAULTS.yml](https://github.com/{{site.repository}}/raw/master/resources/config_files/myconfig_DEFAULTS.yml).

## input_polygons
One of the most important parts of the configuration is setting what vector data to read and to map it to the correct class for lifting. The files, or layers, are all mapped to a single lifting class. With this mapping it defines what rules to use for lifting and stitching each polygon.

### datasets
There are two ways to map datasets to the lifting class depending on the input files. For file formats that have a flat structure and only contain a single layer per file we have the flat configuration. File formats like GML supporting multiple layers can be configured using a nested configuration.

1. Flat config
~~~ yaml
- datasets:
  - D:\data\campus\partof.shp
  - D:\data\campus\otherpart.shp
  uniqueid: FACE_ID
  lifting: Building
- datasets:
  - D:\data\campus\terrain.shp
  uniqueid: TERRAIN_ID
  lifting: Terrain
~~~
The simplest configuration that is mainly used it a flat configuration. It describes a dataset configuration for a single class. First create a new level with a leading dash per *datasets* and add a new level with a leading dash per filename. Make sure the spacing is correct; the filename needs to be at the same level as *uniqueid* and *lifting* and a spacing deeper then *datasets*. There can be multiple files configured for this class in a single dataset configuration. In this case both files must contain the same attribute configured at *uniqueid*.

2. Nested config
~~~ yaml
- datasets:
  - D:\data\campus\polygons.gml
	uniqueid: fid
	lifting_per_layer:
		Buildingfootprint: Building
		Terrainobjects: Terrain
		Waterways: Water
~~~
As an alternative a file with multiple layers can be used for input polygons. In this example the polygons.gml file consists of three layers: *Buildingfootprint*, *Terrainobjects* and *Waterways*. These layers are mapped using the `lifting_per_layer` option by defining the class per layer name.

### uniqueid
This is the name of the attribute in the input file that contains a unique value per feature. In case nothing is put here it uses *fid* by default. **CHECK THIS**
When identifiers from different datasets conflict with each other, many of the output formats will be invalid according to their standard. This is because each object in the final model has to have a unique id, even when originating from different sources. **Note: Multipolygons are split so the identifier is changed by adding a trailing counter**.

### lifting
This is a class name corresponding to the list of the lifting options defined in [Lifting options]({{site.baseurl}}/lifting_options). In total there are 7 different possibilities to choose from;
- Building
- Terrain
- Forest
- Water
- Road
- Separation
- Bridge/Overpass

{% include imagezoom.html file="/settings/settings_input_polygon_class_lifting_class.png" alt="" %}

### height_field
This is the attribute name of an integer attribute. This integer attribute describes the different relative height levels of polygons overlapping in 3D. This describes the vertical relationship between the objects. In Dutch BGT all objects with height_field = 0 describe the ground level. All objects above are counted in positive direction (a bridge would have level 1) and objects below ground level are counted in negative direction (a tunnel would have level -1). [More information on how this is applied in BGT](http://imgeo.geostandaarden.nl/def/imgeo-object/overbruggingsdeel/inwinningsregel-imgeo/toelichting-relatieve-hoogte)

### handle_multiple_heights
When the input polygons use the [height_field](#height_field) this boolean setting describes if these multiple vertical overlapping objects need to be used (height_field attribute with any value) or to just use the ground level polygons (height_field attribute with value of 0).

## input_elevation
Without point clouds the algorithm cannot perform *'3dfying'* of the polygons. Therefore we need to configure the point cloud files and what to read. 

### datasets
1. Single file
~~~ yaml
- datasets:                         # List of data set with specific settings
  - D:\data\top10nl\schie\ahn3.laz  # Definition for one or multiple LAS/LAZ files using the same parameters
~~~
To read a LAS or LAZ file with points the filenames have to be defined in *datasets* section of the *input_elevation* section. First create a new level with a leading dash per *datasets* and add a new level with a leading dash per filename. Make sure the spacing is correct; the filename needs to be a spacing deeper then *datasets*.

2. Directory
~~~ yaml
- datasets:                         # List of data set with specific settings
  - D:\data\top10nl\schie\*         # Definition for reading a full directory with LAS/LAZ files using the same parameters
~~~
Alternatively all files in a directory can be read at once. Instead of writing a filename an *\** is used to define to use all files in that directory. There is no check for file extension so do this for a directory with LAS/LAZ files only.

### omit_LAS_classes
~~~ yaml
omit_LAS_classes:      # Option to omit classes defined in the LAS/LAZ files
  - 1 # unclassified   # ASPRS Standard LiDAR Point Classes classification value
  - 3 # vegetation
  - 4 # vegetation
  - 5 # vegetation
~~~
Many point clouds contain classes that are not used for the creation of the 3D models. Classes like unclassified or vegetation points. These classes can be completely omitted during the reading process. When defined in this section all points with these classifications are ignored. 

### thinning
~~~ yaml
thinning: 10           # Thinning factor for points, this is the amount of points skipped during read, a value of 10 would result in points number 10, 20, 30, 40 being used
~~~
Thinning is used to skip a certain amount of points while reading. When using point clouds with a large amount of points one can use this for getting results quicker. When testing a configuration or area one can set this value to skip points and make the process faster. There is nothing smart done here, the reader simply skips the set amount of points during reading. The equation used here is `point number % thinning value == 0` so only points that have value modulo equal to zero are used. When setting this value to 5 points the reader will use point 5, 10, 15, 20 and so on.---
title: Lifting options
keywords: settings config configuration
sidebar: 3dfier_sidebar
permalink: lifting_options.html
---

In the lifting_options we define the settings used for lifting polygons of the specified class. In total there are 7 different possibilities to choose from;
- Building
- Terrain
- Forest
- Water
- Road
- Separation
- Bridge/Overpass

If a class is not defined the default values will be used. Default values can be found in [resources/config_files/myconfig_DEFAULTS.yml](https://github.com/{{site.repository}}/raw/master/resources/config_files/myconfig_DEFAULTS.yml).

### Building
~~~ yaml
roof:
 height: percentile-90     # Percentile of points within radius of building vertices for roof height lifting, building_radius_vertex_elevation is defined in options
 use_LAS_classes:          # LAS classes to be used for this class. If empty all classes are used
   - 6
ground:
 height: percentile-10     # Percentile of points within radius of building vertices for floor height lifting
 use_LAS_classes:
   - 2
   - 9
lod: 1                      # Define the Level Of Detail for Buildings (0 and 1 possible)
floor: true                 # Set if floors should be created
inner_walls: true           # Set if the walls between to adjacent buildings within a block should be created
triangulate: false          # Set if the output should be triangulated, only works for non-triangular output formats like CityGML/IMGeo
~~~

#### roof and ground settings
Settings within the roof section only apply to the calculation of the roof height where settings in the ground section apply to the calculation of the floor height.

#### height: percentile-xx
To define the height of an object calculated from the points close to it we use percentiles. The percentile defined here is a number between 0-100 written at the location of the **xx** e.g. `percentile-10`. The percentile is calculated from the heights of the points corresponding to the vertex of the object as described in [output_options:radius_vertex_elevation]({{site.baseurl}}/output_options.html#radius_vertex_elevation). All heights are ordered and the value at the configured percentile is calculated.

Lets say we have a list of 10 height values and we configure a percentile of 90. We order the 10 height values from low to height and we will take the 9th value in the list. The 9th value because: `10 values * 90th percentile / 100 = 9`.

Example;
Height values: `1-1-5-7-6-9-6-3-4-2`
Configuration: `percentile-50`
Ordered height values: `1-1-2-3-4-5-6-6-7-9`
Height calculated: `4`

#### use_LAS_classes
This defines what LAS classes to use for the height calculation of its related object class. Only points classified with the LAS classes configured here are allowed to be added to the height lists of the vertices. 

When points are read the algorithm checks:
1. the classification of the point according to the [input_elevation:omit_LAS_classes]({{site.baseurl}}/input_options.html#omit_las_classes)
2. the classification of the points according to this list
3. if the point is within distance set by [output_options:radius_vertex_elevation]({{site.baseurl}}/output_options.html#radius_vertex_elevation)

Lets say this setting corresponds to the class `Road` and the LAS classes configured are `2 (ground) and 11 (road surface)`. All points read that are within distance of the vertex of the road object and have LAS classification of 2 or 11 are added to the heights list of that vertex of the road object.

**If no classes are defined here, all points read are used!**

**TODO: Add image**

#### use_LAS_classes_within
Additional to the [use_LAS_classes](#use_las_classes) there is a more specific version that defines additionally that points have to lay within a polygon to be allowed to be added to the heights list. All 3 steps described in the previous setting are taken plus an additional one that checks if the point falls within the polygon. This can be used for example for objects that are influenced by surrounding objects and their heights like bridges crossing over bridges.

**If use_LAS_classes are not defined but only use_LAS_classes_within, points within distance of a vertex but not within the polygon are not filtered by class and all classes are used!**

**TODO: Add image of vertex buffer including within and image of excluding use_LAS_classes**

#### lod: 0
*Download [YAML]({{site.baseurl}}/assets/configs/lod0.yml) and [OBJ]({{site.baseurl}}/assets/configs/lod0.obj)*

Create a building as a flat surface at ground level. The calculation of the floor height is done identical to LoD 1 and reconstructed with vertical walls from the terrain to the floor surface of the building. The vertical walls are added to make the model water tight and close the terrain surface.

{% include imagezoom.html file="/settings/settings_lod0.png" alt="" %}
{% include imagezoom.html file="/settings/lod0.png" alt="" %}

#### lod: 1
*Download [YAML]({{site.baseurl}}/assets/configs/lod1.yml) and [OBJ]({{site.baseurl}}/assets/configs/lod1.obj)*

Create a building as an extrusion from the ground level to the roof level. The calculation of the floor height is done using the `Building:ground:height` percentile and the roof height using `Building:roof:height` percentile. When the heights are calculated the two surfaces, roof and floor, are built and calculation of the vertical walls is started. The walls are topologically noded to neighbouring buildings. If a building is taller then its neighbour an additional vertex is added at the roof height of its neighbour. Same holds for a buildings floor being lower then that of a neighbour. 

For example we have two buildings with their floors at 0 meters and roofs at 10 and respectively 9 meters height. The building with the roof of 10 meters will contain three vertices at the connecting wall segment, at 0, 9 and 10 meters. The building with the roof of 9 meters will contain two vertices at the connecting wall segment, at 0 and 9 meters. This way we can assure the model to be both watertight and topologically correct.

{% include imagezoom.html file="/settings/settings_lod1.png" alt="" %}
{% include imagezoom.html file="/settings/lod1.png" alt="" %}

#### floor: true
*Download [YAML]({{site.baseurl}}/assets/configs/floor_true.yml) and [OBJ]({{site.baseurl}}/assets/configs/floor_true.obj)*

Reconstruct polygons of the floor of buildings. When set to `true` (and configuration set to `lod: 1`) the model will create a watertight 3D solid object per building.

{% include imagezoom.html file="/settings/settings_floor_true.png" alt="" %}
{% include imagezoom.html file="/settings/floor_true.png" alt="" %}

#### floor: false
*Download [YAML]({{site.baseurl}}/assets/configs/floor_false.yml) and [OBJ]({{site.baseurl}}/assets/configs/floor_false.obj)*

Do not reconstruct polygons of the floor of buildings. The buildings are still connected to the surrounding objects through the walls but do not contain polygons for the floors.

{% include imagezoom.html file="/settings/settings_floor_false.png" alt="" %}
{% include imagezoom.html file="/settings/floor_false.png" alt="" %}

#### inner_walls: true
*Download [YAML]({{site.baseurl}}/assets/configs/inner_walls_true.yml) and [OBJ]({{site.baseurl}}/assets/configs/inner_walls_true.obj)*

Reconstruct polygons for walls in between connected buildings. When set to `true` (and configuration set to `lod: 1`) the model will create walls in between houses that are connected with at least one side.

{% include imagezoom.html file="/settings/settings_inner_walls_true.png" alt="" %}
{% include imagezoom.html file="/settings/inner_walls_true.png" alt="" %}

#### inner_walls: false
*Download [YAML]({{site.baseurl}}/assets/configs/inner_walls_false.yml) and [OBJ]({{site.baseurl}}/assets/configs/inner_walls_false.obj)*

Do not reconstruct polygons for walls in between connected buildings.

{% include imagezoom.html file="/settings/settings_inner_walls_false.png" alt="" %}
{% include imagezoom.html file="/settings/inner_walls_false.png" alt="" %}

#### triangulate: true
Create multisurface geometry for buildings by triangulating polygons. This is the default for output formats that don't support solids or only support multipolygon geometries.

#### triangulate: false
Create solid geometry for buildings, only works when an output format supports solid geometries.

### Water
~~~ yaml
height: percentile-10     # Percentile of points within radius of water vertices for lifting, radius_vertex_elevation is defined in options
use_LAS_classes:          # LAS classes to be used for this class, if empty all classes are used
  - 2 
use_LAS_classes_within:   # LAS classes to be used for this class, but only if points fall within the polygon and the range of the vertex
  - 9
~~~

Class for objects that describe water surfaces. Water objects will be raised to a single height since water is expected to be flat.


### Road
~~~ yaml
height: percentile-50
filter_outliers: true     # Filter outliers by iterative Least Squares fitting of 3D quadric surface, only replace heights of detected outliers
flatten: true             # Filter outliers by iterative Least Squares fitting of 3D quadric surface, replace all heights of polygon with the fitted plane which results in smoother roads
use_LAS_classes:          # LAS classes to be used for this class, if empty all classes are used
  - 2
  - 11
~~~

Class for objects that describe road surfaces. For road object each vertex will be raised to the height calculated by the points that correspond to the vertex. Extra options are available to remove outliers within road objects to overcome unwanted features. 

#### filter_outliers: true
Option to clean road polygons that contain spikes due to used height points. The filter implemented is an iterative Least Squares fitting of a 3D quadratic surface. The input are the polygon vertices with their calculated heights based on the percentile. A minimum of 6 vertices is needed for the filter to work. The filter fits a surface and finds the largest outlier. When the outlier is smaller then 2 sigma of the standard deviation iteration will stop. 

All heights of vertices which are marked as an outlier are replaced with the value represented by the fitted surface. See [flatten](#road-flatten-true) for further replacement of height values.

**Important: not all height points are used for performance reasons, only vertices. Therefore this filter is extremely sensitive to the arrangement of vertices of the road polygon**

**TODO: add image of outlier filter**

#### filter_outliers: false
Disable filtering of outliers using iterative Least Squares fitting of a 3D quadratic surface.

#### flatten: true {#road-flatten-true}
Option to even further fit the road objects towards the quadratic surface calculated by [filter_outliers](#filter_outliers). By enabling this all heights of the polygon vertices are replaced by the values represented by the fitted surface.

**TODO: add image of flatten filter**

#### flatten: false {#road-flatten-false}
Disable height replacement of vertices not marked as outliers by the Least Squares 3D quadratic surface fitting in [filter_outliers](#filter_outliers).

### Separation
~~~ yaml
height: percentile-80
~~~

Class that describes objects inspired by objects contained in the Dutch BGT. These objects generally describe walls or wall-like objects with a flat top.

### Bridge/Overpass
~~~ yaml
height: percentile-50
flatten: true             # Filter outliers by iterative Least Squares fitting of 3D quadric surface, replace all heights of polygon with the fitted plane which results in smoother bridges
~~~

Settings for reconstruction of bridges is handled in this class. Some of the sides of a bridge must be floating with no connection to neighbouring polygons and some of the sides need to be stitched to the neighbours. A rather complex method is used to identify what vertices to stitch to adjacent objects and what vertices to skip so they *'float'*.

The Least Squares 3D quadratic surface fitting as described in [road:filter_outliers](#filter_outliers-true) is always applied.

**TODO: add image of bridge**

#### flatten: true
See [road:flatten: true](#road-flatten-true) for more information. This adjusts all vertices of the bridge to the height of the fitted quadratic surface.

#### flatten: false
See [road:flatten: false](#road-flatten-false) for more information.

### Terrain
~~~ yaml
simplification: 0            # Simplification factor for points added within terrain polygons, points are added random
simplification_tinsimp: 0.1  # Simplification threshold for points added within terrain polygons, points are removed from triangulation until specified error threshold value is reached
innerbuffer: 0.0             # Inner buffer in meters where no additional points will be added within boundary of the terrain polygon
~~~
Since it is a special case the *Terrain* class is adding raw points from the point cloud to the interior of the polygon. Since a terrain polygon can be rather large and contains information about relief within its boundaries more information is needed.

Point clouds normally contain a lot of information. Since we don't want the 3D model to contain lots of unnecessary information, we introduce simplification algorithms. These algorithms make sure only a selection of points from the point cloud are added to the interior of terrain polygons.

#### simplification
Random filtering using a [uniform integer distribution](https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution) between 1 and supplied *simplification* value. A value of 6 generated equal changes compared to throwing a 6-sided dice.

Usage of [simplification_tinsimp](#simplification_tinsimp) is preferred, this is just a cheap alternative.

**Don't use simplification and simplification_tinsimp in the same configuration!**

#### simplification_tinsimp
For simplification of terrain there is an algorithm implemented that describes the terrain with as little triangles as possible. It minimises the amount of triangles while it makes sure the set maximum error is not exceeded. Points that do not add more detail to the terrain when added are not used. The minimum detail it needs to add must be higher then the configured simplification distance. The algorithm is based on the paper of [Heckbert, P. S., & Garland, M. (1997). Survey of polygonal surface simplification algorithms](https://people.eecs.berkeley.edu/~jrs/meshpapers/GarlandHeckbert.pdf). It uses Greedy Insertion to add points in a specific order to the triangle mesh so the point with the largest impact on the terrain is processed first. The points are added iteratively up-to the point the calculated error is less then the configured threshold.

Below are three examples of the impact of this simplification setting. 
1. The first is 0.0 meters which results in all height points within the terrain polygon to be added. It results in a lot of triangles in the yellow terrain polygons. 
2. The second is 0.1 meters where it means that the terrain TIN is allowed to have a maximum deviation of 0.1 meters with the point cloud. The result shows a lot less triangles compared to the first example. In area's with a bit more relief close to the buildings there are quite some triangles left.
3. The final example is 0.5 meters an it results in only points added at locations that greatly influence the height of terrain since the example area is quite flat.

*Download [YAML]({{site.baseurl}}/assets/configs/tinsimp_0.yml) and [OBJ]({{site.baseurl}}/assets/configs/tinsimp_0.obj) of simplification_tinsimp:0*
{% include imagezoom.html file="/settings/settings_tinsimp_0.png" alt="" %}
{% include imagezoom.html file="/settings/tinsimp_0.png" alt="" %}

*Download [YAML]({{site.baseurl}}/assets/configs/tinsimp_01.yml) and [OBJ]({{site.baseurl}}/assets/configs/tinsimp_01.obj) of simplification_tinsimp:0.1*
{% include imagezoom.html file="/settings/settings_tinsimp_01.png" alt="" %}
{% include imagezoom.html file="/settings/tinsimp_01.png" alt="" %}

*Download [YAML]({{site.baseurl}}/assets/configs/tinsimp_05.yml) and [OBJ]({{site.baseurl}}/assets/configs/tinsimp_05.obj) of simplification_tinsimp:0.5*
{% include imagezoom.html file="/settings/settings_tinsimp_05.png" alt="" %}
{% include imagezoom.html file="/settings/tinsimp_05.png" alt="" %}

#### innerbuffer
In case there is a need to prevent points to be added close to polygon boundaries there is the *innerbuffer* setting. When this is used the interior points are only added when the distance to the boundary is a minimum of the supplied value in meters. Below are three examples that show what happens when setting and increasing the *innerbuffer* value.

*Download [YAML]({{site.baseurl}}/assets/configs/innerbuffer_05.yml) and [OBJ]({{site.baseurl}}/assets/configs/innerbuffer_05.obj)*
{% include imagezoom.html file="/settings/settings_innerbuffer_05.png" alt="" %}
{% include imagezoom.html file="/settings/innerbuffer_05.png" alt="" %}

*Download [YAML]({{site.baseurl}}/assets/configs/innerbuffer_1.yml) and [OBJ]({{site.baseurl}}/assets/configs/innerbuffer_1.obj)*
{% include imagezoom.html file="/settings/settings_innerbuffer_1.png" alt="" %}
{% include imagezoom.html file="/settings/innerbuffer_1.png" alt="" %}

*Download [YAML]({{site.baseurl}}/assets/configs/innerbuffer_3.yml) and [OBJ]({{site.baseurl}}/assets/configs/innerbuffer_3.obj)*
{% include imagezoom.html file="/settings/settings_innerbuffer_3.png" alt="" %}
{% include imagezoom.html file="/settings/innerbuffer_3.png" alt="" %}

### Forest
~~~ yaml
simplification: 0            # Simplification factor for points added within forest polygons, points are added random
simplification_tinsimp: 0.1  # Simplification threshold for points added within forest polygons, points are removed from triangulation until specified error threshold value is reached
innerbuffer: 0.0             # Inner buffer in meters where no additional points will be added within boundary of the forest polygon
~~~

All options available for [Terrain](#terrain) are also available for the *Forest* class.

There is no proper solution yet on how to nicely represent vegetation. Therefore in general configurations with only ground points are used for the *Forest* class.---
title: Useful scripts
keywords: resources scripts
sidebar: 3dfier_sidebar
permalink: useful_scripts.html
---

In the [resources folder](https://github.com/{{site.repository}}/raw/master/resources/) in the repository you can find example configuration files, default settings file and helpful scripts for data visualisation and conversion

## Resources folder in repository
Description of the contents of the resources repository

Name | Folder | Description
-----|--------|-------------
3dfier.mtl | | Material description file for coloring an OBJ file
3dfier2sim | x | Source code for program to convert 3dfier OBJ output to OFF file with a volume, see README.md in folder
BGT_prepare | x | Script and GFS files to stroke arcs and convert BGT GML to GeoPackage, see ReadMe.md in folder
BGT_prepare\BGT_gfs_files | x | GFS files describing the content of different BGT GML files
build_ubuntu | x | Scripts to build 3dfier program and dependencies for Ubuntu
build_ubuntu1604.sh | | Script to build 3dfier program and dependencies for Ubuntu 16.04
config_files | x | Example config files
config_files\myconfig.yml | | Example config file
config_files\myconfig_DEFAULTS.yml | | Config file describing programmed defaults
config_files\myconfig_README.yml | | Config file with complete description of all settings options
create_vegetation_off.py | | See script
Example_data | x | Zip file with example input data and configuration for running 3dfier
flowdiagrams | x | Flow diagrams of 3dfier used within the documentation
obj2off.py | | Script to convert an OBJ to OFF file format
splitobj | x | Script for splitting a single OBJ file into several different OBJ files by class
splitobj_cpp | x | C++ program for splitting a single OBJ file into several different OBJ files by class
translate2min-obj.py | | Script for shifting origin of OBJ to 0,0
---
title: Installation
keywords: 3dfier installation ubuntu docker windows compile
summary: These instructions will help you to install 3dfier on various operating systems. For Windows please use the binary files and do not compile from
sidebar: 3dfier_sidebar
permalink: installation.html
---

## Install on Windows using binaries
Binary releases exist only for Windows users. Others will have to follow one of the other installation guides for [Linux](#ubuntu-1604) or [Docker](#docker)
There exists a ready-to-use version of [3dfier for Windows 64-bit](https://github.com/{{site.repository}}/releases/latest). Download and extract the files to any given folder and follow the instructions in the [Get started guide]({{site.baseurl}}/index).

### Release binaries content
Description of files in the zip file

All dll files distributed with 3dfier belong to GDAL or other packages used in the GDAL drivers. Other depencencies used are statically built within the executable. 

Filename | Package
---------|--------
3dfier.exe | 3dfier
expat.dll | GDAL
freexl.dll | GDAL
gdal204.dll | GDAL
geos.dll | GDAL
geos_c.dll | GDAL
hdf5.dll | GDAL
hdf5_hl.dll | GDAL
iconv-2.dll | GDAL
iconv.dll | GDAL
jpeg.dll | GDAL
libcurl.dll | GDAL
libeay32.dll | GDAL
libgmp-10.dll | GDAL
liblzma.dll | GDAL
libmpfr-4.dll | GDAL
libmysql.dll | GDAL
libpng16.dll | GDAL
libpq.dll | GDAL
libxml2.dll | GDAL
lwgeom.dll | GDAL
netcdf.dll | GDAL
ogdi.dll | GDAL
openjp2.dllv | GDAL
proj_5_2.dll | GDAL
spatialite.dll | GDAL
sqlite3.dll | GDAL
ssleay32.dll | GDAL
szip.dll | GDAL
xerces-c_3_2.dll | GDAL
zlib1.dll | GDAL
zstd.dll | GDAL


## macOS

You need to install the following free libraries:

  1. [CMake](http://www.cmake.org)
  1. [CGAL v5.0+](http://www.cgal.org) 
  1. [GDAL](https://gdal.org/)
  1. [yaml-cpp](https://github.com/jbeder/yaml-cpp)
  1. [LASlib](https://github.com/LAStools/LAStools/tree/master/LASlib)
  1. [LASzip](https://github.com/LAStools/LAStools/tree/master/LASzip)

We suggest using [Homebrew](http://brew.sh/) for the first 4:

    $ brew install cgal
    $ brew install gdal
    $ brew install cmake

For LASlib/LASzip, follow the instruction and use the `CMakeLists.txt`.

To compile 3dfier:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make install
    $ 3dfier


## Ubuntu 20.04
### 1. Adding *ubuntugis-unstable* repository
To install *GDAL* on Ubuntu 20.04 LTS it is probably the easiest to add the [*ubuntugis-unstable*](https://launchpad.net/~ubuntugis/+archive/ubuntu/ubuntugis-unstable?field.series_filter=xenial). It contains *GDAL* >= 2.1 under (`libgdal-dev`) package.

*Note: ubuntugis-stable repository doesn't contain any Ubuntu 20.04 packages yet.*

Add the *ubuntugis-unstable* repository by running:
```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
```

### 2. Install dependencies
*CGAL* (`libcgal-dev`), *Boost* (`libboost-all-dev`) and *yaml-cpp* (`libyaml-cpp-dev`) are part of the *Ubuntu Universe* repository.

Once you have all the repos added, you can use a package manager, e.g. `apt` or *Synaptic* to install them. 

E.g. using apt-get:
```
sudo apt-get install libcgal-dev libboost-all-dev libyaml-cpp-dev libgdal-dev
```

### 3. Compile LAStools
*LAStools* contains both the *LASlib* and *LASzip*.

This step requires *CMake* and *UnZip* packages. Some Linux distributions don't have them preinstalled:
```
sudo apt-get install -y unzip cmake
```
Download and compile *LAStools* with the following:
```
wget  http://lastools.github.io/download/LAStools.zip
unzip LAStools.zip
cd LAStools; mv LASlib/src/LASlib-config.cmake LASlib/src/laslib-config.cmake
mkdir build; cd build
cmake ..
sudo make install
```

### 4. Compile 3dfier

Download and compile 3dfier:
```
git clone https://github.com/tudelft3d/3dfier.git
cd 3dfier; mkdir build; cd build
cmake ..
sudo make install
```
Test the installation by trying out the [first run]({{site.baseurl}}/first_run).

## Windows
This guide will talk you through the compilation of 3dfier on Windows 10 64-bit using Visual Studio (steps are identical for Visual Studio 2017 and 2019).

There are some steps to be taken to prepare the build environment for 3dfier. Most important is installing software to compile and downloading the libraries 3dfier is depending on.

*Note 1: Versions used in this guide are the versions used at time of writing. Future versions of libraries could be supported but usage can change and the Visual Studio Solution might need changing for them to work.*

*Note 2: In this guide we build 3dfier in 64-bit (x64) and all dependencies must also be built in 64-bit. For project files created with CMake the `-A x64` switch should explicitly be used.*

*Note 3: Since CGAL 5.0 the library is header only. Building the library is not needed anymore. The current Visual Studio project file in the repository is made for version 4.xx, when using verions >5.0 one should remove the library includes from the project file.*

### 1. Running installers
First you will need to download and install in their default directorties:
1. [Visual Studio Community (2017 or later)](https://www.visualstudio.com/downloads/). Install at least the C++ part.
1. [CMake (3.15 or later)](https://cmake.org/download/), download and install `Windows win64-x64 Installer`. Add variable to the PATH during installation.
1. [Boost precompiled binaries (1.71 or later)](https://sourceforge.net/projects/boost/files/boost-binaries). Pick the latest version that is built for your [MSVC++ compiler version](https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B). Install boost using the installer.
1. [OSGeo4W (with GDAL 2.3.0 or later)](https://trac.osgeo.org/osgeo4w), download the [64-bit installer](http://download.osgeo.org/osgeo4w/osgeo4w-setup-x86_64.exe). From this package you will need to install at least the GDAL package.
1. [CGAL (4.12 or later)](https://github.com/CGAL/cgal/releases), download `CGAL-4.12-Setup.exe` (or newer) and install. Select *GMP and MPFR precompiled libs*, *Environment variables to set CGAL_DIR* and *Add CGAL/auxilary/gmp/lib to the PATH* during setup.

### 2. Compilation of dependencies
Next, we need to download and compile Yaml-cpp and LAStools. 

#### Yaml-cpp
Download [yaml-cpp (0.5.3 or later)](https://github.com/jbeder/yaml-cpp/releases) and extract to e.g. `C:\dev\yaml-cpp`. There are two options of getting the Visual Studio project files using CMake:

1. using CMake GUI ([tutorial here](https://cmake.org/runningcmake/)).

1. using command line. 
Open a Command prompt (press windows button+R, type cmd and press enter) and navigate to the yaml-cpp directory:
```
cd C:\dev\yaml-cpp
```
Generate the Visual Studio build files with
```
mkdir vs_build
cd vs_build
cmake .. -G "Visual Studio 15 2017 Win64"
```

After generation open the Visual Studio file `YAML_CPP.sln`. Set the solution configuration to `Release` in the main toolbar. From the menu bar select Build and click `Build Solution`.

#### LAStools
Download [LAStools](https://rapidlasso.com/lastools/) and extract to e.g. `C:\dev\lastools`.

Use CMake as explained previous for [Yaml-cpp](#yaml-cpp) to generate the Visual Studio solution files.

After generation open the Visual Studio file `LASlib.sln`. Set the solution configuration to `Release` in the main toolbar. From the menu bar select Build and click `Build Solution`.

### 3. Set environment variables
Go to `Control Panel > System > Advanced system settings > Environment Variables` and add the following user variables. Note that the version numbers and the installation paths may be different!
* `BOOST_ROOT`=`C:\boost_1_71_0`
* `BOOST_LIBRARYDIR`=`C:\boost_1_71_0\lib64-msvc-14.0`
* `CGAL_DIR`=`C:\dev\CGAL-4.12`
* `GDAL_ROOT`=`C:\OSGeo4W64`
* `LASLIB_ROOT`=`C:\dev\lastools\LASlib`
* `LASZIP_ROOT`=`C:\dev\lastools\LASzip`
* `OSGEO4W_ROOT`=`C:\OSGeo4W64`
* `YAML-CPP_DIR`=`C:\dev\yaml-cpp`

Go to `Control Panel > System > Advanced system settings > Environment Variables` and add the following directory to Path.
* `C:\OSGeo4W64\bin`

### 4. Compile 3dfier
Download and extract the source code from the menu on the left or fork directly from GitHub. Browse to the vs_build folder and open the Visual Studio file `3dfier.sln`.

If in any case the Visual Studio solution is not working its possible to generate them from the CMake files directly as explained previous for [Yaml-cpp](#yaml-cpp).

### 5. Run 3dfier!
If all is good you should now be able to run 3dfier! Go to the [First run]({{site.baseurl}}/first_run) and start producing models.

* * * 
#### Help, Visual Studio complains that some file can not be found!
Check whether the directories and files specified in the environment variables are correct. Also check these places in Visual Studio:
* the include folders in `Project > Properties > C/C++ > General > Additional Include Directories`
* the library folders in `Project > Properties > Linker > General > Additional Library Directories`
* the libraries files in `Project > Properties > Linker > Input > Additional Dependencies`

Make sure each directory or file exists on your drive. For example: you may need to change a version number somewhere.

## Docker
We offer built docker images from the `master`, `development` branches and each release. You'll find the images and instructions on using them at [Docker Hub](https://hub.docker.com/r/tudelft3d/3dfier).
---
title: Supported file formats
keywords: file formats gdal shapefile geopackage obj cityjson citygml postgis
sidebar: 3dfier_sidebar
permalink: supported_file_formats.html
---

## Input formats
For 3dfier to do its magic you need to input at minimum a set of polygons and a point cloud. This page describes the file formats allowed for both. 

### 1. Polygons
For reading polygon input GDAL is used. This opens up the possibility of reading lots of different file types. Take into account that even though the file format can be read by 3dfier using GDAL, the data doesn't automatically in the right format. Some file formats like GML allow for a feature to have multiple geometry types. When the reader is not configured correctly the outcome can be unexpected. 

For all file formats supported for reading you can look the [OGR Vector drivers list](https://gdal.org/drivers/vector/index.html).

Multipolygons are supported. They will be split into single polygons. The unique identifier will be duplicated and a trailing counter will be added ``id-1, id-2, etc``.

Curvedpolygons are supported starting v1.3. When the reader comes across a curvedpolygon this is discretized using the default OGR settings. This algorithm makes sure the discretization for a curve into a linestring is identical for traversing the line from start->end as from end->start.

### 2. Point cloud
Point cloud reading is implemented using the [LASLib library](https://github.com/LAStools/LAStools/tree/master/LASlib) that is the OpenSource reading library of [LAStools](https://rapidlasso.com/lastools/). This library supports reading of LAS and LAZ point clouds.

## Output formats
The output format is defined using the command-line parameter e.g. `--CityJSON`. The output format is case sensitive.

### CityJSON
[CityJSON](http://www.cityjson.org)

CityJSON is a JSON-based encoding for storing 3D city models, also called digital maquettes or digital twins. CityJSON used the information model of CityGML.

### OBJ
[Wikipedia OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file)

OBJ is a geometry definition file format first developed by Wavefront Technologies. The file format is open and has been adopted by many 3D graphics application vendors. The OBJ file format is a simple data-format that represents 3D geometry alone.

When the goal is to visualise the geometries in a viewer like [MeshLab](http://www.meshlab.net/) this is the best choice as an output format. If you put the material definition file [3dfier.mtl](https://github.com/{{site.repository}}/raw/master/resources/3dfier.mtl) within the same folder of the .obj file, you will even have objects coloured according to their class.

### STL
[Wikipedia STL](https://en.wikipedia.org/wiki/STL_%28file_format%29)

3dfier currently exports to the text version of STL (binary support coming soon).

### CityGML
[CityGML.org](http://www.citygml.org/)

CityGML is a common information model and XML-based encoding for the representation, storage, and exchange of virtual 3D city and landscape models.

Writing of CityGML is supported. The output created consists of a valid CityGML schema as can be tested using this [CityGML schema validator](http://geovalidation.bk.tudelft.nl/schemacitygml/).

To write a separate file for each class there is the format specifier `CityGML-Multifile`. This will add the layer name behind the output filename followed by the proper extension like `filename + layername + .gml`. Take this into account when defining the output filename for example `testdata-` will result in `testdata-Buildings.gml`.

### IMGeo
[IMGeo catalog](https://www.geonovum.nl/geo-standaarden/bgt-imgeo/gegevenscatalogus-imgeo-versie-211)
[BGT-IMGeo](https://www.geonovum.nl/geo-standaarden/bgt-imgeo)

The Dutch basic registration BGT are stored in a CityGML ADE format called IMGeo. This `CityGML-IMGeo` output format creates a valid schema according to the [Geonovum IMGeo 2.1.1 GML Application Schema validator](http://validatie.geostandaarden.nl/etf-webapp/testruns/create-direct?testProjectId=a6a9ddd2-9ab6-3f87-98bc-bbdeb274d679)

To write a separate file for each class there is the format specifier `CityGML-IMGeo-Multifile`. This will add the layer name behind the output filename followed by the proper extension like `filename + layername + .gml`. Take this into account when defining the output filename for example `testdata-` will result in `testdata-Buildings.gml`.

### CSV
For the purpose of calculation building statistics some CSV outputs are implemented. The output writes a CSV with unique id and height statistics per building.

`CSV-BUILDINGS`, `CSV-BUILDINGS-MULTIPLE`, `CSV-BUILDINGS-ALL-Z`

The available options for CSV output are:

**`CSV-BUILDINGS`**

Returns the xx-th percentile of the z-values of the ground and roof points within a building polygon. Values are in **cm**. The percentiles are provided as:

```
lifting_options: 
  Building:
    ground:
      height: percentile-10
    roof:
      height: percentile-90
```
Output:

| id                                    | roof | floor |
|---------------------------------------|------|-------|
| b31e1feb1-00ba-11e6-b420-2bdcc4ab5d7f | 1248 | 20    |

**`CSV-BUILDINGS-MULTIPLE`**

Returns multiple percentiles of the z-values of the ground and roof points within a building polygon. Values are in **m**.

Output:

| id                                    | ground-0.00 | ground-0.10 | ground-0.20 | ground-0.30 | ground-0.40 | ground-0.50 | roof-0.00 | roof-0.10 | roof-0.25 | roof-0.50 | roof-0.75 | roof-0.90 | roof-0.95 | roof-0.99 |
|---------------------------------------|-------------|-------------|-------------|-------------|-------------|-------------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| b31bdd44c-00ba-11e6-b420-2bdcc4ab5d7f | -0.03       | 0.01        | 0.02        | 0.02        | 0.03        | 0.04        | -0.03     | 0.01      | 0.02      | 0.04      | 0.11      | 2.36      | 2.38      | 2.44      |

**`CSV-BUILDINGS-ALL-Z`**

Returns all z-values of the points within a building polygon. Values are in **cm**.

Output:

| id                                    | allzvalues |
|---------------------------------------|------------|
| b31bdfb7a-00ba-11e6-b420-2bdcc4ab5d7f | -3\|-3\|-2\|-2\|-2\|... |




### Shapefile
Writing of Shapefiles is supported in a simplified way; it writes MultiPolygonZ with vertical faces. Not all software has support for vertical faces and will show errors in the output. When using QGIS the file opens as expected.

To write a separate file for each class there is the format specifier `Shapefile-Multifile`. This will add the layer name behind the output filename followed by the proper extension like `filename + layername + .gml`. Take this into account when defining the output filename for example `testdata-` will result in `testdata-Buildings.gml`.

### PostGIS
Output of a PostGIS database is supported. Instead of the file name you have to supply the PostGIS connection string as used for GDAL too:
`3dfier example_data\testarea_config.yml --PostGIS "PG:dbname='3dfier' host='localhost' port='5432' user='username' password='password'"`.
Please make sure to double quote the connection string so its passed as a single parameter.

`PostGIS` specifier creates a database with a table per layer, a column per attribute and a geometry column with type of MultiPolygonZ.

`PostGIS-PDOK` outputs a PostGIS database as described before with an additional column XML that contains the IMGeo XML string of the object

`PostGIS-PDOK-CityGML` outputs a PostGIS database as described before with an additional column XML that contains the CityGML XML string of the object

### GDAL
Besides the list of previous described output formats 3dfier also implements a generic GDAL output driver. Same as for the file reading the output format might not support the geometric output as created. Nevertheless one can try to use this driver at its own discretion. For this driver to work you need to create a new section in the configuration file in which to configure the GDAL driver to use. The driver used needs to support creation of geometries and MultiPolygonZ.

~~~ yaml
output:
  gdal_driver: "ESRI Shapefile"
~~~---
title: Visualisation of various data formats
keywords: visualisation
sidebar: 3dfier_sidebar
permalink: data_visualisation.html
---

### CityJSON
- [CityJSON supported software](https://www.cityjson.org/software/)

### OBJ
- [Meshlab](http://www.meshlab.net/)
- [CloudCompare](https://www.danielgm.net/cc/)
- [Blender](https://www.blender.org/)
- [FME](https://www.safe.com/fme/)

### CityGML
- [Schema Validation](http://geovalidation.bk.tudelft.nl/schemacitygml)
- [CityGML-tools](https://github.com/citygml4j/citygml-tools)
- [Azul](https://github.com/tudelft3d/azul)
- [FME](https://www.safe.com/fme/)

### IMGeo (CityGML ADE)
- [Schema Validation](http://validatie.geostandaarden.nl/etf-webapp/testprojects?testdomain=IMGeo%20GML)
- [Azul](https://github.com/tudelft3d/azul)

### CSV
Using polygon files and table merge in GIS---
title: Your first time running
keywords: first use
summary: These instructions will help you to get started with 3dfier. Other sections will provide additional information on the use of the software.
sidebar: 3dfier_sidebar
permalink: first_run.html
---

3dfier is a command-line utility that does not have a graphical user interface (GUI). This means you need to run it from [Command Prompt](https://www.lifewire.com/how-to-open-command-prompt-2618089).

## Verify installation
First thing to do is to test if the [Installation]({{site.baseurl}}/installation) was successful.
You can do this by opening the Command Prompt (press windows button+R, type cmd and press enter) and dragging the 3dfier program into the Command Prompt window. If you now press enter the 3dfier program should output:

```
3dfier Copyright (C) 2015-2019  3D geoinformation research group, TU Delft
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions; for details run 3dfier with the '--license' option.

ERROR: one YAML config file must be specified.

Allowed options:
  --help                        View all options
  --version                     View version
  --license                     View license
  --OBJ arg                     Output
  --OBJ-NoID arg                Output
  --CityGML arg                 Output
  --CityGML-Multifile arg       Output
  --CityGML-IMGeo arg           Output
  --CityGML-IMGeo-Multifile arg Output
  --CityJSON arg                Output
  --CSV-BUILDINGS arg           Output
  --CSV-BUILDINGS-MULTIPLE arg  Output
  --CSV-BUILDINGS-ALL-Z arg     Output
  --Shapefile arg               Output
  --Shapefile-Multifile arg     Output
  --PostGIS arg                 Output
  --PostGIS-PDOK arg            Output
  --PostGIS-PDOK-CityGML arg    Output
  --GDAL arg                    Output
```

## Minimum system requirements
The software is built with the idea to make creating 3D models as easy as possible. Therefore the software doesn't need supercomputers or special hardware. As long as you are not reconstructing extremely large areas the software allows to use any normal hardware. Most of the developments and testing have been done using a laptop with an Intel Core i7-6560U CPU, 16Gb RAM and SSD for data storage. A simple model can be created even using a laptop with less RAM and slower CPU.

## Command line options
When executing 3dfier without any command line options the text output prints all possible options as shown in [Verify installation](#verify-installation). The `--help` option prints the same text output as when executed without any options. Using `--version` one can see the version of 3dfier and that of some of the libraries used. The `--license` option shows the license message of 3dfier and list the licenses of libraries used.

All other options (marked with Output) are model output formats that can be used to write the output. The option is the file format name followed by the arguments needed for the format. 

In case of `OBJ`, `OBJ-NoID`, `CityGML`, `CityGML-IMGeo`,  `CityJSON`, `CSV-BUILDINGS`, `CSV-BUILDINGS-MULTIPLE`, `CSV-BUILDINGS-ALL-Z` and `Shapefile` the argument is the file name of the output.

In case of `CityGML-Multifile`, `CityGML-IMGeo-Multifile` and`Shapefile-Multifile` the argument is the first part of the file name that will be followed by the input layer name and the file extension. If `arg` is `filename_` the resulting format is `filename_layername.ext`.

In case of `PostGIS`, `PostGIS-PDOK` and `PostGIS-PDOK-CityGML` the argument is a [PostGIS connection string](https://gdal.org/drivers/vector/pg.html) in the format used by GDAL. The string must be surrounded by single quotes so 3dfier understands it as a single option. Example: `'PG:"dbname='databasename' host='addr' port='5432' user='x' password='y'"'`.

## Prepare example data
For this example we use [BGT_Delft_Example.zip](https://github.com/{{site.repository}}/raw/master/resources/Example_data/BGT_Delft_Example.zip) from the GitHub repository located in `3dfier/resources/Example_data/`. Create a folder with 3dfier and the depencency dll's by following the [Installation]({{site.baseurl}}/installation) instructions and add the `example_data folder`.

{% include imagezoom.html file="example_folder.png" caption="Folder layout" alt="Folder layout" %}

## Running with example data
Opening the Command Prompt (press windows button+R, type cmd and press enter) and navigate to the folder where 3dfier is located.

Now go into the example data folder using `cd example_data` and run `..\3dfier.exe testarea_config.yml --OBJ output\myfirstmodel.obj`. Now 3dfier will start processing and when finished it produced your first 3D model!

The output file is written to `example_data\output\myfirstmodel.obj` and console output should be as follows.

```
3dfier Copyright (C) 2015-2019  3D geoinformation research group, TU Delft
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions; for details run 3dfier with the '--license' option.

Config file is valid.
Reading input dataset: bgt\bgt_waterdeel.sqlite
        Layer: waterdeel
        (3 features --> Water)
Reading input dataset: bgt\bgt_ondersteunendwaterdeel.sqlite
        Layer: ondersteunendwaterdeel
        (0 features --> Water)
Reading input dataset: bgt\bgt_onbegroeidterreindeel.sqlite
        Layer: onbegroeidterreindeel
        (81 features --> Terrain)
Reading input dataset: bgt\bgt_wegdeel.sqlite
        Layer: trafficarea
        (151 features --> Road)
Reading input dataset: bgt\bgt_ondersteunendwegdeel.sqlite
        Layer: auxiliarytrafficarea
        (0 features --> Road)
Reading input dataset: bgt\bgt_pand.sqlite
        Layer: buildingpart
        (160 features --> Building)
Reading input dataset: bgt\bgt_begroeidterreindeel.sqlite
        Layer: plantcover
        (126 features --> Forest)
Reading input dataset: bgt\bgt_scheiding.sqlite
        Layer: scheiding
        (57 features --> Separation)
Reading input dataset: bgt\bgt_kunstwerkdeel.sqlite
        Layer: kunstwerkdeel
        (1 features --> Separation)
Reading input dataset: bgt\bgt_overigbouwwerk.sqlite
        Layer: overigbouwwerk
        (0 features --> Separation)
Reading input dataset: bgt\bgt_overbruggingsdeel.sqlite
        Layer: bridgeconstructionelement
        (3 features --> Bridge/Overpass)

Total # of polygons: 570
Constructing the R-tree... done.
Spatial extent: (84,616.468, 447,422.999) (85,140.839, 447,750.636)
Reading LAS/LAZ file: ahn3\ahn3_cropped_1.laz
        (640,510 points in the file)
        (all points used, no skipping)
        (omitting LAS classes: 0 1 )
[==================================================] 100%
Reading LAS/LAZ file: ahn3\ahn3_cropped_2.laz
        (208,432 points in the file)
        (all points used, no skipping)
        (omitting LAS classes: 0 1 )
[==================================================] 100%
All points read in 1 seconds || 00:00:01
3dfying all input polygons...
===== /LIFTING =====
===== LIFTING/ =====
=====  /ADJACENT FEATURES =====
=====  ADJACENT FEATURES/ =====
=====  /STITCHING =====
=====  STITCHING/ =====
=====  /BOWTIES =====
=====  BOWTIES/ =====
=====  /VERTICAL WALLS =====
=====  VERTICAL WALLS/ =====
Lifting, stitching and vertical walls done in 0 seconds || 00:00:00
=====  /CDT =====
=====  CDT/ =====
CDT created in 0 seconds || 00:00:00
...3dfying done.
OBJ output: output\myfirstmodel.obj
Features written in 0 seconds || 00:00:00
Successfully terminated in 2 seconds || 00:00:02
```

## Data that is generated
Only 3D data is generated for 2D polygons supplied in the input. If no polygon is input for terrain it is not automagically generated. There is no option (yet) for the creation of terrain within the polygons extent or point file extent. A [Generate terrain polygon tutorial]({{site.baseurl}}/generate_terrain) is written to create a topologically terrain polygon with building footprints cut out.

## Known issues
### Missing DLL's
If you see an error message about missing MSCVRxxx.dll instead of this text you may need to do one of the following:

* MSCVR100.dll is missing: Install the [Microsoft Visual C++ 2010 Redistributable](https://www.microsoft.com/en-us/download/details.aspx?id=14632)
* MSCVR120.dll is missing: Install the [Microsoft Visual C++ 2013 Redistributable](https://www.microsoft.com/en-us/download/details.aspx?id=40784)
* MSCVR140.dll is missing: Install the [Microsoft Visual C++ 2015 Redistributable](https://www.microsoft.com/en-us/download/details.aspx?id=48145)
---
title: What does it do?
keywords: 3dfier algorithm details
summary: 
sidebar: 3dfier_sidebar
permalink: what_does_it_do.html
---

What the 3dfier algorithm actually do to make 3D objects out of 2D polygons and point clouds? This page explains the steps implemented to create 3D models.

The input of the algorithm is a set of topologically connected polygons. An example is shown below from a horizontal point of view. All object belong to one of the 7 object classes supported by 3dfier;
- Terrain
- Forest
- Water
- Road
- Separation
- Building
- Bridge / Overpass

{% include imagezoom.html file="algorithm/objects_horizontal_view.png" alt="" %}

The objects when viewed in a GIS software using the same colors as in the horizontal view.
{% include imagezoom.html file="algorithm/objects_overview_color.png" alt="" %}


On top of the objects we drape the points read from a point cloud.
{% include imagezoom.html file="algorithm/objects_horizontal_view_with_points.png" alt="" %}

Using the points read from the point cloud we select all points within the 2D buffer distance of a vertex. The height is stored into a list of heights per vertex. In the configuration the lifting percentile per object type is stored. The type of the object is requested and using the percentile and all heights stored the height of the vertex for this object is calculated. This process is repeated for each vertex in the imported 2D polygons. To describe objects in the `terrain` and `forest` classes in a more realistic way additional height points are added within the 2D polygons. These additional height points represent better the height differences. This process is called `lifting`.
{% include imagezoom.html file="algorithm/objects_horizontal_view_lifted.png" alt="" %}

After lifting the polygons from 2D to 3D many connected polygons will have vertical height jumps. These height jumps are shown in the left side of the picture below where on the right side these are solved and vertical walls (red) are added. Since the goal of the algorithm is to create a watertight model we have to solve for these holes in the model. This is explained in the next step.
{% include imagezoom.html file="algorithm/3d_objects_with_holes_filled.png" alt="" %}

Now all polygons are lifted to height we will solve for the holes we created. Based on type of object we say it is `soft` or `hard`, the best description would be nature and man-made. Heights of soft features are allowed to be moved towards that of a hard feature. Also `water` objects are forced to stay flat and its heights will not be moved. The result of this process is a 3D topologically connected surface containing all 2D polygons as 3D objects as shown below. This process is called `stitching`. 
{% include imagezoom.html file="algorithm/objects_horizontal_view_stitched.png" alt="" %}

From the created model there are two output options for buildings;
1. Trianglated MultiSurfaces
2. Solids

In the image below one can see MultiSurface output where the walls and roof of the building are stored in the object as a collection of triangles.
{% include imagezoom.html file="algorithm/objects_horizontal_view_multisurface.png" alt="" %}

In the image below one can see buildings output as Solid objects where the walls and roof are stored as rectangular surfaces.
{% include imagezoom.html file="algorithm/objects_horizontal_view_solid.png" alt="" %}---
title: Different LoD of buildings
keywords: LoD buildings examples
sidebar: 3dfier_sidebar
permalink: lod_buildings.html
---

## What's the difference between LoD0 and LoD1 buildings?
The most simple way to explain the difference between LoD0 and LoD1 is to use the way we describe a building in real life. The floor, walls and roof of a building together make a building. When reconstructing an LoD0 building only the floor is created as a flat surface within the surrounding terrain. A LoD1 building adds the walls and roof to a building. In the LoD1 version the roof is describe as a flat surface.


**IMPORTANT: LoD0 is currently implemented for OBJ output only**

## LoD0 reconstruction
*Download [YAML]({{site.baseurl}}/assets/configs/lod0.yml) and [OBJ]({{site.baseurl}}/assets/configs/lod0.obj)*

**Point classes used: ground points**

Buildings created in LoD0 are flat surfaces. The surface height is calculated according to the ground points detected around the footprint. The percentile statistics determine the final height. All objects around the building are stitched to these heights or in case of water an additional vertical wall is created. Result is that terrain around the building and the buildings flat surface reconstruction is identical to that of the LoD1 reconstruction.

The LoD0 reconstruction is particularly useful when not in need of buildings in the model. E.g. if there is an existing LoD2 model and a watertight terrain model is to be created.

{% include imagezoom.html file="/lod0-delft.png" alt="" %}
{% include imagezoom.html file="/lod0-delft-zoom.png" alt="" %}

Use the following settings for the example dataset to reconstruct LoD0 buildings:
~~~ yaml
input_polygons:
  - datasets: 
      - bgt\bgt_pand.sqlite
    uniqueid: gml_id
    lifting: Building
    height_field: relatievehoogteligging
~~~
~~~ yaml
lifting_options: 
  Building:
    lod: 0
    floor: true
    triangulate: false
    ground:
      height: percentile-10
      use_LAS_classes:
        - 2
        - 9
    roof:
      height: percentile-90
      use_LAS_classes:
        - 6
~~~
~~~ yaml
input_elevation:
  - datasets:
      - ahn3\ahn3_cropped_1.laz
      - ahn3\ahn3_cropped_2.laz
    omit_LAS_classes:
      - 0 # Never classified
      - 1 # Unclassified
    thinning: 0

options:
  building_radius_vertex_elevation: 3.0
  radius_vertex_elevation: 1.0
  threshold_jump_edges: 0.5
~~~

## LoD1 reconstruction
*Download [YAML]({{site.baseurl}}/assets/configs/lod1.yml) and [OBJ]({{site.baseurl}}/assets/configs/lod1.obj)*

**Point classes used: ground and non-ground points**

The reconstruction of buildings in LoD1 result in cubes with flat roof surfaces and extruded walls between floor and roof surfaces. Height for the floor is calculated like explained in [LoD0 reconstruction](#lod0-reconstruction). Roof height is based on all points in distance of the vertices and points within the footprint. Points are filtered by class set with the *use_LAS_classes* setting. The percentile statistics determine the final height.

{% include imagezoom.html file="/lod1-delft.png" alt="" %}

Use the following settings for the example dataset to reconstruct LoD1 buildings:
~~~ yaml
input_polygons:
  - datasets: 
      - bgt\bgt_pand.sqlite
    uniqueid: gml_id
    lifting: Building
    height_field: relatievehoogteligging
~~~
~~~ yaml
lifting_options: 
  Building:
    lod: 1
    floor: true
    inner_walls: true
    triangulate: false
    ground:
      height: percentile-10
      use_LAS_classes:
        - 2
        - 9
    roof:
      height: percentile-90
      use_LAS_classes:
        - 6
~~~
~~~ yaml
input_elevation:
  - datasets:
      - ahn3\ahn3_cropped_1.laz
      - ahn3\ahn3_cropped_2.laz
    omit_LAS_classes:
      - 0 # Never classified
      - 1 # Unclassified
    thinning: 0

options:
  building_radius_vertex_elevation: 3.0
  radius_vertex_elevation: 1.0
  threshold_jump_edges: 0.5
~~~---
title: Minimal data requirements
keywords: 3dfier homepage
sidebar: 3dfier_sidebar
permalink: minimal_data_requirements.html
---

## What to input?
Using the software can have different reasons and end-goals. These examples describe what the minimal data requirements are to create a 3D model using 3dfier. The first thing needed is polygons describing the objects to make 3D. Next is the height information to be added to these polygons in the form of a point cloud. Both input types must have a file format from the list of [Supported file formats]({{site.baseurl}}/supported_file_formats.html).

Besides the data files the software needs a configuration file formatted as YAML. This configuration contains all settings except the output format. All possibilities are explained in [Settings]({{site.baseurl}}/input_options.html). In the following sections various possibilities and their minimum data requirements are covered.

{% include imagezoom.html file="extrusion.png" alt="Data input and extrusion steps" %}

## Buildings only
For the creation of a 3D model consisting of only LoD1 buildings there is the need for the following data:
- Polygons of building footprints
- Point cloud with minimum classification:
	- Ground
	- Non-ground

*Why the need for classification in ground/non-ground for points?* While developing the algorithms thorough testing with non classified point clouds have been done. The results showed that for the algorithm to work it cannot do without differentiating between ground and non-ground. The main reason is due to vegetation. The noise introduced by vegetation cannot be overcome by the statistical analysis.

In this example configuration the ground height for the building is calculated by using only the points classified as ground (class 2). The height of the roof of buildings is calculated by using all point classes. This will also take into account ground points for the roof. However using a high percentile makes the result less prone to picking a ground points as the roof height.

*Note: having more LAS classes gives better results since there is less noise within the configured classes. E.g. there is (almost) no vegetation points in the building classes so building heights are not interfered by trees.*

~~~yaml
input_polygons:
  - datasets: 
      - bgt\bgt_pand.sqlite
    uniqueid: gml_id
    lifting: Building
    height_field: relatievehoogteligging

lifting_options: 
  Building:
    ground:
      height: percentile-10
      use_LAS_classes:
        - 2
    roof:
      height: percentile-90

input_elevation:
  - datasets:
      - ahn3\ahn3_cropped_1.laz
      - ahn3\ahn3_cropped_2.laz
~~~

## All objects
For the creation of a 3D model containing all objects in an area there is the need for the following data:
- Topologically connected polygons
- Point cloud with minimum classification:
	- Ground
	- Non-ground

In this example we configured all classes to used the ground/non-ground classification. All lifting_options of objects related to ground-like objects (terrain, forest, water and road) are set to only use ground (class 2) and other objects to use all points. 

~~~yaml
input_polygons:
  - datasets: 
      - bgt\bgt_waterdeel.sqlite
      - bgt\bgt_ondersteunendwaterdeel.sqlite
    uniqueid: gml_id
    lifting: Water
    height_field: relatievehoogteligging
  - datasets: 
      - bgt\bgt_onbegroeidterreindeel.sqlite
    uniqueid: gml_id
    lifting: Terrain
    height_field: relatievehoogteligging
  - datasets: 
      - bgt\bgt_wegdeel.sqlite
      - bgt\bgt_ondersteunendwegdeel.sqlite
    uniqueid: gml_id
    lifting: Road
    height_field: relatievehoogteligging
  - datasets: 
      - bgt\bgt_pand.sqlite
    uniqueid: gml_id
    lifting: Building
    height_field: relatievehoogteligging
  - datasets: 
      - bgt\bgt_begroeidterreindeel.sqlite
    uniqueid: gml_id
    lifting: Forest
    height_field: relatievehoogteligging
  - datasets: 
      - bgt\bgt_scheiding.sqlite
      - bgt\bgt_kunstwerkdeel.sqlite
      - bgt\bgt_overigbouwwerk.sqlite
    uniqueid: gml_id
    lifting: Separation
    height_field: relatievehoogteligging
  - datasets: 
      - bgt\bgt_overbruggingsdeel.sqlite
    uniqueid: gml_id
    lifting: Bridge/Overpass
    height_field: relatievehoogteligging

lifting_options: 
  Building:
    ground:
      height: percentile-10
      use_LAS_classes:
        - 2
    roof:
      height: percentile-90
  Terrain:
    simplification_tinsimp: 0.1
    use_LAS_classes:
      - 2
  Forest:
    simplification_tinsimp: 0.1
    use_LAS_classes:
      - 2
  Water:
    height: percentile-10
  Road:
    height: percentile-50
    use_LAS_classes:
      - 2
  Separation:
    height: percentile-80
  Bridge\Overpass:
    height: percentile-50

input_elevation:
  - datasets:
      - ahn3\ahn3_cropped_1.laz
      - ahn3\ahn3_cropped_2.laz
~~~---
title: Simplification of terrain
keywords: examples terrain forest simplification
sidebar: 3dfier_sidebar
permalink: terrain_simplification.html
---

## Why simplify terrain and forest?
Modelling in 3D creates a tremendous amount of data. Maintaining and using large datasets presents a challenge. As long as the data that is represented by a 3D model is of added value overcoming these challenges are not a problem. When creating terrain models the algorithm uses Triangulated Irregular Networks (TIN). This is a collection of connected 3D triangles that form a closed surface. The amount of triangles, and therefore vertices, is key in the total storage size of a terrain model.

As described in [What does it do?]({{site.baseurl}}/what_does_it_do.html) Terrain and Forest class get added additional height points during the reconstruction. These points are the height points that are read from the point files. If the algorithm would add all points the resulting density would be extremely high while the added value would be low for the height information. The software contains several options to filter these additional points. These options are described in the lifting options for [Terrain]({{site.baseurl}}/lifting_options.html#terrain) and [Forest]({{site.baseurl}}/lifting_options.html#forest). The simple version is random simplification opposed to TIN simplification that is a smart error minimisation algorithm.

## Random simplification
When configured to use random simplification the algorithm uses a random number generator to decide if a point is used or not. Random filtering using a [uniform integer distribution](https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution) between 1 and configured *simplification* value. A value of 6 generated equal changes compared to throwing a 6-sided dice.

It is not advised to use this random simplification for production stage. It creates a different outcome each run and does not create a model that describes the terrain in the best way possible. For testing and research purposes this filtering works a lot faster.

## TIN simplification
A more robust algorithm for simplification of terrain and forest classes is implemented. It minimises the amount of triangles while it makes sure the set maximum error is not exceeded. Points that do not add more detail to the terrain when added are not used. The minimum detail it needs to add must be higher then the configured simplification distance. The algorithm is based on the paper of [Heckbert, P. S., & Garland, M. (1997). Survey of polygonal surface simplification algorithms](https://people.eecs.berkeley.edu/~jrs/meshpapers/GarlandHeckbert.pdf). It uses Greedy Insertion to add points in a specific order to the triangle mesh so the point with the largest impact on the terrain is processed first. The points are added iteratively up-to the point the calculated error is less then the configured threshold.

In the image the black dashed line called true ground surface is the surface described by the point clouds. The error defined by the TIN and the point cloud is the smallest distance from the point to the TIN.

{% include imagezoom.html file="/tinsimp_error_specification.png" alt="Terrain simplification error definition" %}

In the example below there are three cases that show the impact of the `simplification_tinsimp` setting. On the left is an example of using all points for the Terrain without simplification. The middle permits a maximum error of 5 centimeter (`simplification_tinsimp: 0.05`) and the right a maximum error of 15 centimeter (`simplification_tinsimp: 0.15`).
{% include imagezoom.html file="/terrain_simplification_tinsimp.png" alt="" %}---
title: Building footprints from OpenStreetMap
keywords: footprints OSM
sidebar: 3dfier_sidebar
permalink: building_footprints_from_openstreetmap.html
---

This guide is about extracting building footprints as a Shapefile of polygons from the [OpenStreetMap](https://www.openstreetmap.org) dataset, in order to use them as input for 3dfier.

## Download OSM data 

First, you need to download the OSM data for your area of interest:

1. Through your browser, visit the [OpenStreetMap](https://www.openstreetmap.org) website

1. Zoom at the area you want to work on

1. Press the `Export` button on the top

1. From the left panel, you may select to further specify the area you want to download (via the `Manually select a different area` option)

1. When finished, press the `Export` button from the left panel.

You should be asked to download an `.osm` file. Just store it somewhere on your computer.

## Extract buildings through QGIS 

Originally, the data from OpenStreetMap are just geometries with key-value pairs assigned to them. You can easily filter the buildings from all geometries inside an `.osm` file.

1. Open QGIS.

1. From the menu, select `Layers`->`Add Layer`->`Add Vector layer...`.

1. Select the `.osm` file with the area downloaded.

1. When prompted about the layer you want to add, you can only select the `multipolygons` layer (it's OK if you add all of them, but buildings are on this layer).

1. From the `Layers Panel`, right-click on the multipolygons layer you've just added and select `Filter...`.

1. On the dialog, provide the following expression: `"building" is not null` (essentially, what that means is that all polygons without a _building_ key, will be filtered out). You should now see only the buildings on the map.

1. From the `Layers Panel`, right-click on the multipolygons layer, again, and now select `Save as...`.

1. On the dialog, select the `ESRI Shapefile` format and provide the output file. You may also want to reproject the geometries on another CRS, now, if the elevation data you are going to use on 3dfier are not on WGS 84 (EPSG:4326).

When you save, you should have a shapefile with footprint of the buildings for this area, at the CRS you specified. You can use this, now, as an input for 3dfier.---
title: Generate terrain polygon
keywords: terrain
sidebar: 3dfier_sidebar
permalink: generate_terrain.html
---

Without a 2D polygon describing terrain extent, 3dfier is not able to create a terrain from a point cloud without. The solution is to create a polygon of the area of interest save this as a shapefile.

We make sure the terrain and buildings are stitched nicely and the buildings do not intersect with the terrain. We do this by using symetrical difference between the terrain polygon and the buildings. For nice tutorial about this see section *D. Symmetrical Difference* at [GrindGIS](http://grindgis.com/software/qgis/basic-editing-tools-in-qgis).

### Terrain
{% include imagezoom.html file="https://user-images.githubusercontent.com/30265851/32222599-30d38b70-be3a-11e7-90ad-000218305924.png" external=true %}

### 2D terrain + buildings
{% include imagezoom.html file="https://user-images.githubusercontent.com/30265851/32222677-7acc1a3a-be3a-11e7-8737-c0637005a6ba.png" external=true %}

### Symmetrical Difference operation result
{% include imagezoom.html file="https://user-images.githubusercontent.com/30265851/32222645-62cf1964-be3a-11e7-9a70-746222d0c53a.png" external=true %}

### 3D model as OBJ in Meshlab
{% include imagezoom.html file="https://user-images.githubusercontent.com/30265851/32222728-b87a5112-be3a-11e7-892e-32e4afee01ca.png" external=true %}

### Viewpoint below the 3D model
{% include imagezoom.html file="https://user-images.githubusercontent.com/30265851/32222770-e95b2a40-be3a-11e7-88c2-c026412ebc9e.png" external=true %}

Thanks to [@antoinebio](https://github.com/antoinebio) for the images supplied in [issue #48](https://github.com/tudelft3d/3dfier/issues/48).---
title: Generate LoD1 models
keywords: terrain
sidebar: 3dfier_sidebar
permalink: generate_lod1.html
---

Goal: In this tutorial you will learn how to generate 3D city models from open data using extrusion and the open-source software [3dfier](https://github.com/tudelft3d/3dfier).

## Introduction to extrusion

One popular way to generate 3D city models is extrusion: features from a 2D dataset such as a cadastral database, are lifted to a single height, creating volumetric 3D city models.
The heights are usually obtained from laser scanning (e.g. the average elevation of all points within a footprint), from a cadastral database, or from volunteered geoinformation (e.g. using the number of floors).
The first case is illustrated:

{% include imagezoom.html file="extrusion.png" alt="Data input and extrusion steps" %}

This is method is simple and straightforward.
The resulting 3D city models, while simple (they have flat tops; the so-called LOD1 models), offer much advantantage over 2D datasets.
For example, they may be used for shadow analyses and line of sight predictions.

However, there are some challenges to be aware of, e.g. errors in the 2D data propagate to the generated 3D model, while simple this method is not followed by many implementations, and because point clouds are usually large the calculation of height of each feature may be quite slow.

Here at the [3D geoinformation group at TU Delft](https://3d.bk.tudelft.nl) we have developed <a href="https://github.com/tudelft3d/3dfier">3dfier</a> for creating 3D models.
In this tutorial we will briefly demonstrate how to generate a 3D model with it, using open data.

## Installation of the software

The software is command-line, that is, it doesn't have a graphical interface.
However, it is still very simple to use.
The first step is to install it following the [installation instructions]({{site.baseurl}}/installation). 
To test whether you installed it correctly, just run `./3dfier` and you should get something like this:

```
$ ./3dfier 
3dfier Copyright (C) 2015-2019  3D geoinformation research group, TU Delft
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions; for details run 3dfier with the '--license' option.

ERROR: one YAML config file must be specified.

Allowed options:
  --help                        View all options
  --version                     View version
  --license                     View license
  --OBJ arg                     Output
  --OBJ-NoID arg                Output
  --CityGML arg                 Output
  --CityGML-Multifile arg       Output
  --CityGML-IMGeo arg           Output
  --CityGML-IMGeo-Multifile arg Output
  --CityJSON arg                Output
  --CSV-BUILDINGS arg           Output
  --CSV-BUILDINGS-MULTIPLE arg  Output
  --CSV-BUILDINGS-ALL-Z arg     Output
  --Shapefile arg               Output
  --Shapefile-Multifile arg     Output
  --PostGIS arg                 Output
  --PostGIS-PDOK arg            Output
  --PostGIS-PDOK-CityGML arg    Output
  --GDAL arg                    Output
```

## Short overview of 3dfier

3dfier requires one or more 2D datasets and one or more elevation datasets as input.
The 2D datasets can be of any OGR format, such as SHP or GML. The elevation datasets will be in LAS/LAZ.

Besides defining input datasets, 3dfier enables defining certain parameters such as point cloud thinning.

All these are defined in a single text file with the `.yml` extension.
For example, if you have building footprints in a separate 2D file this is how the input will look like:

```
input_polygons:
    datasets: 
      - bgt/bgt_pand.sqlite
    uniqueid: gml_id
    lifting: Building
```

Don't worry, 3dfier comes with a prepared sample config file, which you can edit to adapt to your case.


## Generating your own dataset

After downloading the software you may have noticed that there is a directory called `example_data`. It contains everything you need to create a sample 3D city model which we will use in this short tutorial.

The configuration file `testarea_config.yml` is prepared with all the required information.

First familiarise yourself with the input datasets.
In this example for 2D data we will use [BGT](https://www.kadaster.nl/bgt), that is, the Dutch large-scale topographic map. BGT is open data courtesy of [Kadaster](https://www.kadaster.nl), the national mapping agency of the Netherlands.
The folder `bgt` contains the 2D dataset in multiple files. The area we will work on is the centre of Delft (you can use the free [QGIS](http://www.qgis.org/en/site/) to view the files):

{% include imagezoom.html file="delft-bgt.png" alt="" %}

So the 2D dataset contains not only buildings, but also other features such as water, roads, and vegetation.

For the elevation we will use the National height model of the Netherlands ([AHN](https://www.pdok.nl/nl/ahn3-downloads)). It is available as open data as well.
The area comes in two files, both stored in the folder `ahn3`. Here is how the point cloud looks like (you can use the free software [CloudCompare](http://www.danielgm.net/cc/) for this):

{% include imagezoom.html file="delft-ahn3.png" alt="" %}

The information about the input point clouds is noted in the configuration file as well:

```
input_elevation:
  - datasets:
      - ahn3/ahn3_cropped_1.laz
      - ahn3/ahn3_cropped_2.laz
    omit_LAS_classes:
      - 1 # unclassified
    thinning: 0
```

Notice that in the configuration file you can also specify omitting certain classes in the point cloud (such as vegetation) and thinning the points to speed up the process.
These two lidar files are small, hence we will take into account all lidar points.

Now that you have checked the input datasets, let's have a look at some other options of 3dfier.
One important option is to specify how the elevation of the top of each building is determined.
In our case:

```
lifting_options:
  Building:
    roof:
        height: percentile-90
    ground:
        height: percentile-10
```

the elevation of the 'roof' is at the 90th percentile of the elevation of all points that are within the building footprint.
This should roughly correspond to the elevation of the top of the building (a value of 90 percentile is given to filter out outliers and features such as chimneys).
The elevation of the bottom (that is, the ground plate) is at the 10th percentile.
You can play with these values and determine what works best for you.
Some people prefer to rather use `height_roof: percentile-50` to get the height of the top at the median of all points.

Another important option is the format of the resulting 3D city model. The options are shown above at [Installation of the software](#installation-of-the-software). The option is passed as a command line argument.

3dfier offers [CityJSON](https://www.cityjson.org), [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file), and many more formats.
OBJ is widely supported by 3D computer graphics software, so you can create a nice render of the 3D model:

{% include imagezoom.html file="ams-dof.png" alt="" %}

On the other hand, CityJSON is a powerful 3D GIS format, enabling spatial analyses and structuring of objects.

Now we are ready to generate the 3D model, both in CityJSON and OBJ.
Generating the 3D model requires only one simple command:

```
3dfier testarea_config.yml --OBJ output/testarea.obj --CityJSON output/testarea.json
```

3dfier will report on the process of the 3D generation, but overall for this example it should not take more than half a minute.
If everything went well with the input data, the file should be available in the directory where you specified it (in our case in `/output`).

If you generated an OBJ you can view it with the free software [MeshLab](http://meshlab.sourceforge.net).
If you opted for CityJSON you can visualise it with our [CityJSON web-viewer](https://tudelft3d.github.io/CityJSON-viewer/).

An OBJ is composed of triangles, so the result will look like this:

{% include imagezoom.html file="delft-meshlab.png" alt="" %}

3dfier also comes with a material file, so if you switch the corresponding options in MeshLab you can visualise different semantic classes:

{% include imagezoom.html file="meshlab-options.png" alt="" %}

{% include imagezoom.html file="delft-meshlab-colors.png" alt="" %}

CityJSON looks similar:

{% include imagezoom.html file="delft-webviewer.png" alt="" %}

Congratulations, you have created a 3D city model!

## Try it with your own data

Now that you have familiarised yourself with 3dfier, you can try to generate a 3D city model with your own data.
Let us know if you generate a nice dataset so we can showcase it on our website.