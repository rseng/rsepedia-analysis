# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased
## Changed:
- update documentation with table of features
- drop python 3.5 due to problematic lazperf dependency, includes CI on python 3.7 and 3.8 

## 0.4.2 - 2020-09-18
## Added:
- faster implementation of spatial selection

## Changed:
- all elements in WKT and shapefiles containing multiple polygons are considered (not only the first element)
- the validity conditions for multipolygons have been relaxed: valid adjacent polygons are now accepted

## Fixed:
- bug in copy_point_cloud with masks and zero-filled feature-arrays

## 0.4.1 - 2020-05-20
## Added:
- select_equal filter accepts list of values to compare to the points' attributes
- also the attribute-based filter functions optionally return a mask to allow filter combinations

## Fixed:
- bug in writing/reading 'None' as parameter in the PLY comments

## 0.4.0 - 2020-05-13
## Added:
- build_volume module
- the most relevant functions can be now imported directly from laserchicken
- reading/writing of binary PLY and LAZ files, with optional writing of selected attributes
- utility function to merge point-cloud data
- extra log tasks implemented: point-cloud log entries are introduced upon point-cloud loading, filtering, normalizing, merging and assigning to targets.
- select_polygon now supports multi-polygons and optionally return a mask for the selected points

## Changed:
- compute neighborhoods returns generator with neighborhoods instead of nested neighborhoods like it did before (breaking change!)
- Some of the existing modules have been renamed/restructured (breaking changes!):
    - `normalization` --> `normalize`
    - `feature_extraction` created (functions moved from `feature_extractor/__init__.py`)
    - `select` and `spatial_selection` merged into `filter`, with the function `select_polygon` allowing to deal with all the spatial selection functionalities
    - format-specific `read_*` and `write_*` modules replaced by `load` and `export`
- Dependency on `laspy` replaced by `pylas` + `lazperf` (easier reading/writing of LAS/LAZ files)

## 0.3.2 - 2019-12-12
## Added
- Features added:
    - max_intensity
    - min_intensity
    - range_intensity
    - mean_intensity
    - std_intensity
    - coeff_var_intensity

## Changed
- Features renamed:
    - max_norm_z --> min_normalized_height
    - min_norm_z --> max_normalized_height
    - range_norm_z --> range_normalized_height
    - mean_norm_z --> mean_normalized_height
    - std_norm_z --> std_normalized_height
    - coeff_var_norm_z --> coeff_var_normalized_height
    - density_absolute_mean_norm_z --> density_absolute_mean_normalized_height
    - entropy_norm_z --> entropy_normalized_height
    - kurto_norm_z --> kurto_normalized_height
    - skew_norm_z --> median_normalized_height
    - var_z --> skew_normalized_height
    - perc_1_norm_z --> var_normalized_height
    - perc_100_norm_z --> perc_1_normalized_height

## 0.3.1 - 2019-09-25
## Added
- Percentiles 1-100
- Percentiles normalized height 1-100
- Band ratio feature extractor 
- Function to list available feature extractors
- Tutorial notebook

## Changed
- Echo ratio no longer gives percentage (removed factor 100)

## Fixed
- Bug in reading ply file with comments in unexpected format
- Bug in normal vector and slope

## Removed

## 0.3.0 - 2019-02-27
## Added
- Normalization module
- General tests that all current and future feature extractors will be checked against.
- Possibility to have a randomly subsampled (fixed) number of neighbors (eg for faster feature calculation)
- Added feature extractors based on normalized height

## Changed
- Many feature calculations are done in a vectorized way
- z_entropy feature renamed to entropy_z for consistency
- density_absolute_mean feature renamed to density_absolute_mean_z for consistency
- perc_10 features renamed to perc_10_z for consistency

## Fixed
- Corner cases for many feature extractors (e.g. zero points)
- Bug in slope feature calculation
- Bug in fix density absolute mean feature calculation

## Removed

## [0.1.0] - 2018-04-17
# Point Cloud File Formats -work in progress
Handeling huge point clouds requires a spatial data structure to efficiently access each point.
There is a multitude of File Formats that can be used with LAS and LAZ files. A summary of [File Formats](http://www.cloudcompare.org/doc/wiki/index.php?title=FILE_I/O), out of which we would focus on: 

* [CSV](https://docs.python.org/3/library/csv.html)

* [PCD](http://pointclouds.org/documentation/tutorials/pcd_file_format.php) - the Point Cloud Library format [(PCL)](http://pointclouds.org/), python binding [NLeSc/python-pcl](https://github.com/NLeSC/python-pcl)

* [PLY](http://paulbourke.net/dataformats/ply/) - Polygon File Format or the [Stanford Triangle Format](http://www.graphics.stanford.edu/data/3Dscanrep) 

## CSV:

The CSV (Comma Separated Values) file format is the most common import and export format for spreadsheets and databases. Easy to use within python but once dealing with huge point clouds with extra attributes the file will become too big to process easily. And it is completely non standard.

## PCD:

The [PCD](http://pointclouds.org/documentation/tutorials/pcd_file_format.php) is used as a file format to support 3D point cloud data. Please refer to [pointclouds.org](http://pointclouds.org/documentation/tutorials/pcd_file_format.php) for reference and more information.
Each PCD file contains a header (ASCII) that identifies and declares certain properties of the point cloud data stored in the file.
This is PCL's format which is very complicated to install. 

**HEADER:** 

The header entries must be specified precisely in the following order:
```
VERSION - PCD file version
FIELDS -  name of each dimension/field that a point can have, for example: xyz, rgb (colors), surface normals, moment invariants (j1-3) and more...
* FIELDS x y z                                
* FIELDS x y z rgb                            
* FIELDS x y z normal_x normal_y normal_z     
* FIELDS j1 j2 j3                             
SIZE - size of each dimension in bytes
TYPE - type of each dimension as a char.
* I - represents signed types int8 (char), int16 (short), and int32 (int)
* U - represents unsigned types uint8 (unsigned char), uint16 (unsigned short), uint3(unsigned int)
* F - represents float types
COUNT - specifies how many elements does each dimension have. For example, x data usually has 1 element, but a feature descriptor like the VFH has 308. Basically this is a way to introduce n-D histogram descriptors at each point, and treating them as a single contiguous block of memory. Default: count=1.
WIDTH - width of the point cloud dataset in the number of points, one of the 2 meanings: 
* total number of points in the cloud (equal with POINTS see below) for unorganized datasets
* total number of points in a row of an organized point cloud dataset.
HEIGHT - height of the point cloud dataset in the number of points, one of the 2 meanings:
* HEIGHT=1 for unorganized datasets (thus used to check whether a dataset is organized or not).
* specify the height (total number of rows) of an organized point cloud dataset
VIEWPOINT - specifies an acquisition viewpoint for the points in the dataset. The viewpoint information is specified as a translation (tx ty tz) + quaternion (qw qx qy qz). Default: VIEWPOINT 0 0 0 1 0 0 0
POINTS - total number of points in the cloud. 
DATA - data type that the point cloud data is stored in: ascii/binary.
```
**Attributes**

PCL comes with a variety of pre-defined point types, ranging from SSE-aligned structures for XYZ data, to more complex n-dimensional histogram representations such as PFH (Point Feature Histograms). 
A list of all the [PCD pre-defined point types](https://github.com/PointCloudLibrary/pcl/blob/master/common/include/pcl/impl/point_types.hpp).


* You can [add custom point type](http://pointclouds.org/documentation/tutorials/adding_custom_ptype.php) to the point cloud.

## PLY:

Format for storing graphical objects that are described as a collection of polygons. 
PLY is composed of an header (in ASCII format) followed by a list of vertices and a list of polygons. The header specifies the elements of a mesh and their types and states what properties are associated with each vertex, such as (x,y,z) coordinates, normals and color. The polygon faces are simply lists of indices into the vertex list, and each face begins with a count of the number of elements in each list. The header is followed by the list of elements.

Please refer to [paulbourke- PLY file format](http://paulbourke.net/dataformats/ply/) for reference and more information.

**HEADER**

```
ply
format ascii 1.0
comment Mars model by Paul Bourke
element vertex 259200
property float x
property float y
property float z
element face 516960
property list uchar int vertex_indices
end_header
```
* ply - file format
* format - ascii 1.0, binary_little_endian 1.0, binary_big_endian 1.0
* comment - a comment which is ignored
* element - a description of how some particular data element is stored and how many of them there are. Hence, in a file where there are 12 vertices, each represented as a floating point (X,Y,Z) triple, one would expect to see:
```
element vertex 12
property float x
property float y
property float z
```
Other properties may indicate any data items that are stored at each vertex (for example return number) and the data type of that information. Data types can be one of two variants, depending on the source of the ply file. The type can be specified with one of  [char uchar short ushort int uint float double], or one of  [int8 uint8 int16 uint16 int32 uint32 float32 float64].

At the end of the header, there must always be the line:
```
end_header
```
A sample file addapted from [paulbourke](http://paulbourke.net/dataformats/ply/example1.ply):

```
ply
format ascii 1.0
comment Mars model by Paul Bourke
element vertex 259200
property float x
property float y
property float z
element face 516960
property list uchar int vertex_indices
end_header
15081.5 -3.45644e+06 65.8061
15081 -3.45659e+06 197.422
15078.2 -3.45648e+06 329.009
15075.4 -3.45663e+06 460.597
15071.2 -3.4567e+06 592.148
15065.6 -3.45674e+06 723.653
15059.9 -3.457e+06 855.16
15050.7 -3.45674e+06 986.473

     lots of vertices follow

14541.2 3.33642e+06 -698.464
14547.7 3.33663e+06 -571.58
14551.5 3.33649e+06 -444.589
14552.7 3.336e+06 -317.541
14556.9 3.33645e+06 -190.56
14558.7 3.33661e+06 -63.5247
3 0 721 1
3 721 0 720
3 1 722 2
3 722 1 721
3 2 723 3
3 723 2 722

     lots of triangular facets follow
```
**Attributes**

A variety of properties can be stored, including: xyz, color, surface normals, texture coordinates and data confidence values.

The structure of a typical PLY file:
```
  Header
  Vertex List
  Face List
  (lists of other elements)
```
* Adding new attributes - applications can create new properties that are attached to elements of an object. The format for defining a new element is exactly the same as for vertices, faces and edges. 

**Log**

A log or comments can be placed in the header by using the word comment at the start of the line. Everything from there until the end of the line is ignored.
```
  comment this is ignored
```
**Point Cloud libraries**

* PLY can be used in [CloudCompare](http://www.cloudcompare.org/) which can run on Windows, MacOS and Linux. 

* [python-plyfile](https://github.com/dranjan/python-plyfile) is a NumPy-based text/binary PLY file reader/writer for Python. Dependencies python2 >= 2.7 or python3. The test suite was designed to test python2.7 and 3.4 but works for python3.5. 

Installing: ```pip install plyfile ```

## Size comparison

ODM: [OPALS](http://geo.tuwien.ac.at/opals/html/index.html) - 23.3 MB

LAS: 7.15 MB

CSV/TXT: 9.24 MB

PCD: 3.36 MB

PLY: 5.04 MB

## Notes:

* [NLeSc load and query](https://github.com/NLeSC/pointcloud-benchmark): NLeSc has a pointcloud Python package which can be used to deal with point clouds in various Point Cloud Data Management Systems (PCDMS's), such as [Postgresql](https://www.postgresql.org/), MonetDB, Oracle and a combination of LAStools tools (lassort, lasindex, lasmerge and lasclip).

# Internal data structure

To stay close to the chosen file format the python data structure will look like this:
```
{'log': [{'time': '2018-01-18 16:01', 'module': 'load', 'parameters': [], 'version': '0.9.2'},
         {'time': '2018-01-18 16:01', 'module': 'filter', 'parameters': [('z', 'gt', '1.5')], 'version': '0.9.2'}],
   'pointcloud':
       {'offset': {'type': 'double', 'data': 12.1}},
   'vertex':
       {'x': {'type': 'float', 'data': np.array([0.1, 0.2, 0.3])},
        'y': {'type': 'float', 'data': np.array([0.1, 0.2, 0.3])},
        'z': {'type': 'float', 'data': np.array([0.1, 0.2, 0.3])},
        'return': {'type': 'int', 'data': np.array([1, 1, 2])}}}
        }
 }
 ```
 After we compute some features we enrich the data structure with extra attributes:
 ```
 {'log': [{'time': '2018-01-18 16:01', 'module': 'load', 'parameters': [], 'version': '0.9.2'},
         {'time': '2018-01-18 16:01', 'module': 'filter', 'parameters': [('z', 'gt', '1.5')], 'version': '0.9.2'}],
   'pointcloud':
       {'offset': {'type': 'double', 'data': 12.1}},
   'vertex':
       {'x': {'type': 'float', 'data': np.array([0.1, 0.2, 0.3])},
        'y': {'type': 'float', 'data': np.array([0.1, 0.2, 0.3])},
        'z': {'type': 'float', 'data': np.array([0.1, 0.2, 0.3])},
        'return': {'type': 'int', 'data': np.array([1, 1, 2])}}},
        'eigen_val_1': {'type': 'float', 'data': [0.1, 0.5,  0.25 ])},
        'eigen_val_2': {'type': 'float', 'data': [0.02, 0.05,  0.025 ])},
        ...
        'echo_ratio': {'type': 'float', 'data': np.array([0.05, 0.04, 0.36])}
  }     
    
 ```
 
 This gives us the three data types that we want to store in the memory (and file):
 * Logging
 * Points with unlimited custom attributes
 * Point cloud wide attributes (offset for instance).
 The logging will be stored in COMMENT lines in the PLY file. Example:
 ```
 ply
format ascii 1.0
comment [
comment {'time': '2018-01-18 16:01', 'module': 'load', 'parameters': [], 'version': '0.9.2'}
comment {'time': '2018-01-18 16:01', 'module': 'filter', 'parameters': [('z', 'gt', '1.5')], 'version': '0.9.2'}
comment ]
element vertex 3
property float x
property float y
property float z
property int return
element pointcloud 1
property double offset
end_header
0.1 0.1 0.1 1
0.2 0.2 0.2 1
0.3 0.3 0.3 2
12.1
```
Please cite the software if you are using it in your scientific publication.
<p align="left">
  <img src="https://raw.githubusercontent.com/eEcoLiDAR/laserchicken/master/laserchicken_logo.png" width="500"/>
</p>

[![Build Status](https://travis-ci.org/eEcoLiDAR/laserchicken.svg?branch=master)](https://travis-ci.org/eEcoLiDAR/laserchicken)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/6e3836750fe14f34ba85e26956e8ef10)](https://www.codacy.com/app/c-meijer/eEcoLiDAR?utm_source=www.github.com&amp;utm_medium=referral&amp;utm_content=eEcoLiDAR/eEcoLiDAR&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/eEcoLiDAR/eEcoLiDAR/badge.svg)](https://coveralls.io/github/eEcoLiDAR/eEcoLiDAR)
[![DOI](https://zenodo.org/badge/95649056.svg)](https://zenodo.org/badge/latestdoi/95649056)
[![Documentation Status](https://readthedocs.org/projects/laserchicken/badge/?version=latest)](https://laserchicken.readthedocs.io/en/latest/)

Toolkit for handling point clouds created using airborne laser scanning (ALS). Find neighboring points in your point cloud and describe them as feature values. Read our [user manual](https://laserchicken.readthedocs.io/) and our (very modest) [tutorial](https://github.com/eEcoLiDAR/laserchicken/blob/master/tutorial.ipynb).

# Installation
Prerequisites:
- Python 3.6 or higher
- pip
```
pip install laserchicken
```

#### Necessary steps for making a new release
* Check citation.cff using general DOI for all version (option: create file via 'cffinit')
* Create .zenodo.json file from CITATION.cff (using cffconvert)  
```cffconvert --validate```  
```cffconvert --ignore-suspect-keys --outputformat zenodo --outfile .zenodo.json```
* Set new version number in laserchicken/_version.txt
* Check that documentation uses the correct version
* Edit Changelog (based on commits in https://github.com/eecolidar/laserchicken/compare/v0.3.2...master)
* Test if package can be installed with pip (`pip install .`)
* Create Github release
* Upload to pypi:  
```python setup.py sdist bdist_wheel```  
```python -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/*```  
(or ```python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*``` to test first)
* Check doi on zenodo


## Feature testing

All features were tested for the following general conditions:
- Output consistent point clouds and don't crash with artificial data, real data, all zero data (x, y or z), data without points, data with very low number of neighbors (0, 1, 2)
- Input should not be changed by the feature extractor

The specific features were tested as follows.

*Echo ratio*

A test was written with artificial data to check the correctness of the calculation with manually calculated ratio. Also tested on real data to make sure it doesn't crash, without checking for correctness. We could add a test for correctness with real data but we would need both that data and a verified ground truth.

*Eigenvalues*

Only sanity tests (l1>l2>l3) on real data and corner cases but no actual test for correctness. The code is very simple though and mainly calls numpy.linalg.eig.

*Height statistics (max_z','min_z','mean_z','median_z','std_z','var_z','coeff_var_z','skew_z','kurto_z)*

Tested on real data for correctness. It is however unclear where the ground truths come from. Code is mainly calling numpy methods that do all the work already. Only calculations in our code are:

```
range_z = max_z - min_z
coeff_var_z = np.std(z) / np.mean(z)
```
   
I don't know about any packages that could provide an out of the box coefficient of variance. This is probably because the calculation is so simple.

*Pulse penetration ratio*

Tested for correctness using artificial data against manually calculated values. No comparison was made with other implementations.

*Sigma_z*

Tested for correctness using artificial data against manually calculated values. No comparison was made with other implementations.

*Percentiles*

Tested for correctness using a simple case with artificial data against manually calculated values.

*point_density*

Tested for correctness on artificial data.



# Contributing

## Pull Request Submission Guidelines

Before you submit your pull request consider the following guidelines: 
* Please communicate with us up front about any new feature you would like to add, to avoid disappointment later. You can do this by creating an [issue](https://github.com/eEcoLiDAR/eEcoLiDAR/issues)
* Fork the repository to your own github account if you don't have write access
* Clone the repository on your local machine
* Make your changes in a new git branch:  
`git checkout -b my-fix-branch master`
* Install the development environment:  
`python setup.py develop`
* Make your changes and add tests demonstrating that you fixed the bug or covering the new feature you added
* Order your imports:  
`isort -w 120 your_changed_file.py`
* Format your code according the the project standard:  
`yapf -i your_changed_file.py`
* Check that your code is clean and fix any issues:  
`prospector your_changed_file.py`
* Run tests and make sure they pass:  
`python setup.py test`
* Commit your changes and upload:
  ```
  git add changed_file_1.py changed_file_2.py
  git commit -m 'Your commit message'
  git push
  ```
* Create a [pull request](https://github.com/eEcoLiDAR/eEcoLiDAR/pulls)
# Spatial selection using the `filter` module

Filter is a module which provides functionality to run a range spatial selection over a set of points. The range should be specified as a **Polygon**, stored either in a WKT string or file, or in ESRI shapefile while points are defined as a *point_cloud* (**pc**).

## Two step spatial selection

Range selections are defined in two steps: filtering and refinement. The filtering step uses the polygon's bounding-box and the extent of the *point-cloud* to verify if they overlap. If they overlap, then in the refinement step **contains** from *shapely* is used to extract all points within the Polygon boundaries. 

## Examples

### Example with a WKT string.
```
from laserchicken import load
from laserchicken.filter import select_polygon

pc_in = load("testdata/AHN2.las")
wkt_string = "POLYGON(( 243590.0 572110.0, 243640.0 572160.0, 243700.0 572110.0, 243640.0 572060.0, 243590.0 572110.0 ))"
pc_out = select_polygon(pc_in, wkt_string)
```

### Example with a WKT file.
```
#At the moment only the first polygon will be used for the range selection.
from laserchicken import load
from laserchicken.filter import select_polygon

pc_in = load("testdata/AHN2.las")
filename = "testdata/ahn2_geometries_wkt/ahn2_polygon.wkt"
pc_out = select_polygon(pc_in, filename, read_from_file=True)
```

### Example with a ESRI shapefile.
```
#At the moment only the first polygon will be used for the range selection.
from laserchicken import load
from laserchicken.filter import select_polygon

pc_in = load("testdata/AHN2.las")
filename = "testdata/ahn2_geometries_shp/ahn2_polygon.shp"
pc_out = select_polygon(pc_in, filename, read_from_file=True)
```
.. Laserchicken documentation master file, created by
   sphinx-quickstart on Thu Oct  3 15:28:52 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Laserchicken's documentation!
========================================
Laserchicken is a user-extendable, cross-platform Python tool for
extracting statistical properties (features in machine learning jargon)
of flexibly defined subsets of point cloud data.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

User manual
===========

Laserchicken processes point clouds from Airborne Laser Scanning in LAS/LAZ format, and normalizes, and filters the points, finds neighbors for targeted points and calculates features.

.. image:: figures/workflow.png
  :alt: Laserchicken workflow

The figure shows the default workflow for which Laserchicken was intended. In short, a point cloud is loaded from a LAS or LAZ or PLY file. After this, points can be filtered by various criteria, and the height can be normalized. A point cloud of targets are loaded. See below a description of what is the concept of a target point. For every target point, neighbors will be computed. Based on the list of neighbors, features are extracted that effectively describe the neighborhood of each target point.

Environment and target point cloud
==================================

In Laserchicken, the LiDAR dataset is referred to as the environment point cloud (EPC), and the subsets of points over which a metric is to be calculated are referred to as neighborhoods. Each neighborhood is defined by a target volume and a target point (e.g. a cube of a certain size and its centroid, respectively), with all points enclosed in the volume constituting the neighborhood.

See below an illustration of a target point cloud (green points) representing the centroids of a regular grid cell, and the neighborhoods (red points) defined by a square infinite cell target volume (red columns). Features are calculated over the neighborhood of each target point and then associated with the target point, thus forming the enriched target point cloud (eTPC). Points that are not included in the neighborhoods are shown in in black.

.. image:: figures/targets.png
  :width: 450
  :alt: targets and environment point cloud

Four volume definitions are implemented: an infinite square cell, an infinite cylinder, a cube and a sphere.

.. image:: figures/cell.png
  :width: 300
  :alt: cell
.. image:: figures/cylinder.png
  :width: 300
  :alt: cylinder
.. image:: figures/sphere.png
  :width: 300
  :alt: sphere
.. image:: figures/voxel.png
  :width: 300
  :alt: voxel

All target points together form the target point cloud (TPC). The TPC can be freely defined by the user, and can be for instance identical to the environment point cloud or alternatively a regular grid as illustrated above.
Features are calculated over the list of neighborhoods, with the feature values being associated with each neighborhood's defining target point, thus forming the enriched target point cloud (eTPC).

Modules
=======

Each module from the workflow is described below and an example of its usage is given. If you need more examples, be sure to have a look of the unit tests in the Laserchicken's source code.

Load
----

The load module provides functionality to load  point cloud datasets provided in ASPRS LAS/LAZ, or in PLY format, and is used for both input point clouds. In conjunction with
the PDAL library (https://pdal.io/), this provides access to a comprehensive range of point cloud data formats.

Example from the tutorial notebook::

   from laserchicken import load
   point_cloud = load('testdata/AHN3.las')

Normalize
---------

A number of features require the normalized height above ground as input. Laserchicken provides the option of internally constructing a digital terrain model (DTM) and deriving this quantity. To this end, the EPC is divided into small cells 1m or 2.5m squared). The lowest point in each cell is taken as the height of the DTM. Each point in the cell is then assigned a normalized height with respect to the derived DTM height. This results in strictly positive heights and smooths variations in elevation on scales larger than the cell size. The normalized EPC can be used directly in further analysis, or serialized to disk.

Example from the tutorial notebook::

   from laserchicken.normalize import normalize
   normalize(point_cloud)

Filter
------
Laserchicken provides the option of filtering the EPC prior to extracting features. Points may be filtered on the value of a single attribute relative to a specified threshold (e.g. above a certain normalized height above ground), or on specific values of their attributes (e.g. LAS standard classification). It is also possible to filter with (geo-)spatial layers such as polygons (e.g. regions of interest, land cover types), i.e. selectively including or excluding points.

Example of spatial filtering from the tutorial notebook::

   from laserchicken.filter import select_polygon
   polygon = "POLYGON(( 131963.984125 549718.375000," + \
                      " 132000.000125 549718.375000," + \
                      " 132000.000125 549797.063000," + \
                      " 131963.984125 549797.063000," + \
                      " 131963.984125 549718.375000))"
   point_cloud = select_polygon(point_cloud, polygon)

Example of applying a filter on the theshold of an attribute::

   from laserchicken.filter import select_above, select_below
   points_below_1_meter = select_below(point_cloud, 'normalized_height', 1)
   points_above_1_meter = select_above(point_cloud, 'normalized_height', 1)


Compute neighbors
-----------------

The Compute neighbors module constructs the neighborhoods as defined by the TPC
and target volume by identifying the points in the EPC which reside in the specified volume
centered on the target points, returning each as a list of indices to the EPC. This essential step of computing neighboring points for large samples of points is computationally expensive. Laserchicken uses the optimized ckDtree class (kdTrees are a space-partitioning data structure) provided by the scipy library to organize both the EPC
and TPC in kdTrees in an initial step prior to the computation of neighbors, subsequently accelerating the process of computing neighbors by using the indices of the points with respect to the kDtrees.

Example from the tutorial notebook::

   from laserchicken import compute_neighborhoods
   from laserchicken import build_volume
   targets = point_cloud
   volume = build_volume('sphere', radius=5)
   neighborhoods = compute_neighborhoods(point_cloud, targets, volume)

Note that in the above example, ``neighborhoods`` is a generator and can only be iterated once. If you would want to do multiple calculations without recalculating the neighbors, you can copy the neighborhoods to a list. This is not done by default because neighborhoods can quickly grow quite large so that available RAM unnecessarily becomes the bottle neck.

Features
--------

Feature extraction requires the EPC, the TPC, the computed list of neighborhoods, and a list of requested features as input. For each target point it selects the points of the associated neighborhood and calculates a vector of the requested features over these. This feature vector is appended to the target point, thus defining the eTPC.

Currently, a number of features are implemented, including percentiles of the height distribution and eigenvectors. Computationally expensive calculations requiring multi-dimensional linear algebraic operations (e.g. eigenvectors and eigenvalues) have been vectorized using the einsum function of the numpy library to optimize performance. Their implementation can serve as a
template for new features requiring similar operations.

Example from the tutorial notebook::

   from laserchicken import compute_features
   compute_features(point_cloud, neighborhoods, targets, ['std_z','mean_z','slope'], volume)

Features can be parameterized. If you need different parameters than their defaults you need to register them with these prior to using them.

Example of adding a few parameterized band ratio features on different attributes::

   from laserchicken import register_new_feature_extractor
   from laserchicken.feature_extractor.band_ratio_feature_extractor import BandRatioFeatureExtractor
   register_new_feature_extractor(BandRatioFeatureExtractor(None,1,data_key='normalized_height'))
   register_new_feature_extractor(BandRatioFeatureExtractor(1,2,data_key='normalized_height'))
   register_new_feature_extractor(BandRatioFeatureExtractor(2,None,data_key='normalized_height'))
   register_new_feature_extractor(BandRatioFeatureExtractor(None,0,data_key='z'))

The currently registered features can be listed as follows::

   from laserchicken.feature_extractor import list_feature_names
   sorted(list_feature_names())

Which outputs something like::

   ['band_ratio_1<normalized_height<2',
    'band_ratio_2<normalized_height',
    'band_ratio_2<normalized_height<3',
    'band_ratio_3<normalized_height',
    'band_ratio_normalized_height<1',
    'band_ratio_z<0',
    'coeff_var_norm_z',
    'coeff_var_z',
    'density_absolute_mean_norm_z',
    'density_absolute_mean_z',
    'echo_ratio',
    'eigenv_1',
    'eigenv_2',
    'eigenv_3',
    'entropy_norm_z',
    'entropy_z',
    'kurto_norm_z',
    'kurto_z',
    'max_norm_z',
    'max_z',
    'mean_norm_z',
    'mean_z',
    'median_norm_z',
    'median_z',
    'min_norm_z',
    'min_z',
    'normal_vector_1',
    'normal_vector_2',
    'normal_vector_3',
    'perc_100_normalized_height',
    'perc_100_z',
    'perc_10_normalized_height',
    'perc_10_z',
    ...
    'perc_99_normalized_height',
    'perc_99_z',
    'perc_9_normalized_height',
    'perc_9_z',
    'point_density',
    'pulse_penetration_ratio',
    'range_norm_z',
    'range_z',
    'sigma_z',
    'skew_norm_z',
    'skew_z',
    'slope',
    'std_norm_z',
    'std_z',
    'var_norm_z',
    'var_z']

The following table includes the list of all implemented features:

.. table::
   :widths: 35 30 25 10

   ============================================================  ========================================================================================================================================  ===========================================================  ==================================
    Feature name                                                  Formal description                                                                                                                        Example of use                                               Refs.
   ============================================================  ========================================================================================================================================  ===========================================================  ==================================
   Point density (``point_density``)                             :math:`N/A` where :math:`S` is the neighborhood target volume (area) for finite (infinite) cells                                          Point cloud spatial distribution                             |
   Pulse penetration ratio (``pulse_penetration_ratio``)         :math:`N_{\mathrm{ground}}/N_{\mathrm{tot}}`                                                                                              Tree species classification                                  :cite:`yu2014`
   Echo ratio (``echo_ratio``)                                   :math:`N_{\mathrm{sphere}}/N_{\mathrm{cylinder}}`                                                                                         Roof detection                                               :cite:`car2009`
   Skewness (``skew_z``) [a]_                                    :math:`1/\sigma^3 \cdot \sum{(Z_i - \bar{Z})^3/N}`                                                                                        Vegetation, ground, and roof classification and detection    :cite:`Crosilla2013`
   Kurtosis (``kurto_z``) [a]_                                   :math:`1/\sigma^4 \cdot \sum{(Z_i - \bar{Z})^4/N}`                                                                                        Vegetation, ground, and roof classification and detection    :cite:`Crosilla2013`
   Standard deviation (``std_z``) [a]_ [b]_                      :math:`\sqrt{\sum{(Z_i - \bar{Z})^2/(N - 1)}}`                                                                                            Classification of reed within wetland                        :cite:`zlinszky2012`
   Variance (``var_z``) [a]_                                     :math:`\sum{(Z_i - \bar{Z})^2/(N - 1)}`                                                                                                   Classification of reed within wetland                        :cite:`zlinszky2012`
   Sigma Z (``sigma_z``) [a]_                                    :math:`\sqrt{\sum{(R_i - \bar{R})^2/(N - 1)}}` where :math:`R_i` is  the residual after plane fitting                                     |                                                            :cite:`zlinszky2012`
   Minimum Z (``min_z``) [a]_ [b]_                               :math:`Z_{\mathrm{min}}`                                                                                                                  Simple digital terrain model in wetlands                     :cite:`zlinszky2012`
   Maximum Z (``max_z``) [a]_ [b]_                               :math:`Z_{\mathrm{max}}`                                                                                                                  Height and structure of forests                              :cite:`naesset2002`
   Mean Z (``mean_z``) [a]_ [b]_                                 :math:`\sum{Z_{i}}/N`                                                                                                                     Height and structure of forests                              :cite:`naesset2002`
   Median Z (``median_z``) [a]_                                  :math:`Z_{\mathrm{median}}`                                                                                                               Height and structure of forests                              :cite:`naesset2002`
   Range Z (``range_z``) [a]_ [b]_                               :math:`|Z_{\mathrm{max}} - Z_{\mathrm{min}}|`                                                                                             Height and structure of forests                              :cite:`naesset2002`
   Percentiles Z (``perc_X_z`` with ``X`` in (0:100]) [a]_       Height of every :math:`10^{\mathrm{th}}` percentile.                                                                                      Height and structure of forests                              :cite:`naesset2002`
   Eigenvalues (``eigenv_X``, with ``X`` in (1,2,3))             :math:`\lambda_1, \lambda_2, \lambda_3 ` , with :math:`|\lambda_1| \ge |\lambda_2| \ge |\lambda_3|`                                       Classification of urban objects                              :cite:`weinmann2017`
   Normal vector (``normal_vector_X``, with ``X`` in (1,2,3))    eigen vector :math:`\vec{v}_3`                                                                                                            Roof detection                                               :cite:`Dorninger2008`
   Slope (``slope``)                                             :math:`\tan(\mathrm{arccos}(\vec{v}_3\cdot\vec{k}))` , where :math:`\vec{k} = [0,0,1]^T`                                                  Planar surface detection                                     :cite:`doi:10.1002/esp.3606`
   Entropy Z  (``entropy_z``) [a]_                               :math:`-\sum_{i}{P_i \cdot \mathrm{log}_2{P_i}}`, with :math:`P_i = N_i/\sum_{j}{N_j}`  and :math:`N_i` points in bin :math:`i`           Foliage height diversity                                     :cite:`Bae2014`
   Coefficient variance Z (``coeff_var_z``) [a]_ [b]_            :math:`\frac{1}{\bar{Z}} \cdot \sqrt{\sum{\frac{(Z_i - \bar{Z})^2}{N - 1}}}`                                                              Urban tree species classification                            :cite:`koma2016urban`
   Density absolute mean (``density_absolute_mean_z``) [a]_      :math:`100 \cdot \sum [Z_i > \bar{Z}]/N`                                                                                                  Urban tree species classification                            :cite:`koma2016urban`
   Band ratio (``band_ratio_Z1<z<Z2``) [c]_                      :math:`N_{Z_1<Z<Z_2}/N_{\mathrm{tot}}` where :math:`Z_1` and :math:`Z_2` are user-provided values                                         Height and structure of forests                              |
   ============================================================  ========================================================================================================================================  ===========================================================  ==================================

.. [a] Also available for the normalized height (e.g. ``mean_normalized_height``)
.. [b] Also available for the intensity (e.g. ``mean_intensity``)
.. [c] Fully customizable in variable and range

Below is an example. The figure visualizes the slope feature for a small neighborhood size. We used the same target point cloud as the environment point cloud. The image was generated using mayavi plotting software (https://docs.enthought.com/mayavi/mayavi/).

.. image:: figures/slope.png
  :width: 450
  :alt: slope on small point cloud

Export
------

Laserchicken can write to PLY or LAS/LAZ format for further analysis with the user's choice of software. The PLY format is preferred, as it is flexibly extendable and is the only format Laserchicken will write provenance data to. It is also a widely supported point cloud format.

Example from the tutorial notebook::

   from laserchicken import export
   export(point_cloud, 'my_output.ply')

.. bibliography:: bibliography.bib
   :cited:
   :style: unsrt


