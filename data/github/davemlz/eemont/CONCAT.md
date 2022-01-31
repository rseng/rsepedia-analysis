<p align="center">
  <a href="https://github.com/davemlz/eemont"><img src="https://raw.githubusercontent.com/davemlz/eemont/master/docs/_static/header2.png" alt="header"></a>
</p>
<p align="center">
    <em>A python package that extends Google Earth Engine</em>
</p>
<p align="center">
<a href='https://pypi.python.org/pypi/eemont'>
    <img src='https://img.shields.io/pypi/v/eemont.svg' alt='PyPI' />
</a>
<a href='https://anaconda.org/conda-forge/eemont'>
    <img src='https://img.shields.io/conda/vn/conda-forge/eemont.svg' alt='conda-forge' />
</a>
<a href="https://opensource.org/licenses/MIT" target="_blank">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License">
</a>
<a href='https://eemont.readthedocs.io/en/latest/?badge=latest'>
    <img src='https://readthedocs.org/projects/eemont/badge/?version=latest' alt='Documentation Status' />
</a>
<a href="https://github.com/davemlz/eemont/actions/workflows/tests.yml" target="_blank">
    <img src="https://github.com/davemlz/eemont/actions/workflows/tests.yml/badge.svg" alt="Tests">
</a>
<a href="https://github.com/sponsors/davemlz" target="_blank">
    <img src="https://img.shields.io/badge/GitHub%20Sponsors-Donate-ff69b4.svg" alt="GitHub Sponsors">
</a>
<a href="https://www.buymeacoffee.com/davemlz" target="_blank">
    <img src="https://img.shields.io/badge/Buy%20me%20a%20coffee-Donate-ff69b4.svg" alt="Buy me a coffee">
</a>
<a href="https://ko-fi.com/davemlz" target="_blank">
    <img src="https://img.shields.io/badge/kofi-Donate-ff69b4.svg" alt="Ko-fi">
</a>
<a href="https://developers.google.com/earth-engine/tutorials/community/developer-resources" target="_blank">
    <img src="https://img.shields.io/badge/GEE%20Community-Developer%20Resources-00b6ff.svg" alt="GEE Community">
</a>
<a href="https://twitter.com/dmlmont" target="_blank">
    <img src="https://img.shields.io/twitter/follow/dmlmont?style=social" alt="Twitter">
</a>
<a href='https://joss.theoj.org/papers/34696c5b62c50898b4129cbbe8befb0a'>
    <img src='https://joss.theoj.org/papers/34696c5b62c50898b4129cbbe8befb0a/status.svg' alt='JOSS' />
</a>
<a href="https://github.com/psf/black" target="_blank">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black">
</a>
<a href="https://pycqa.github.io/isort/" target="_blank">
    <img src="https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336" alt="isort">
</a>
</p>

---

**GitHub**: [https://github.com/davemlz/eemont](https://github.com/davemlz/eemont)

**Documentation**: [https://eemont.readthedocs.io/](https://eemont.readthedocs.io/)

**PyPI**: [https://pypi.org/project/eemont/](https://pypi.org/project/eemont/)

**Conda-forge**: [https://anaconda.org/conda-forge/eemont](https://anaconda.org/conda-forge/eemont)

**Tutorials**: [https://github.com/davemlz/eemont/tree/master/docs/tutorials](https://github.com/davemlz/eemont/tree/master/docs/tutorials)

**Paper**: [https://joss.theoj.org/papers/10.21105/joss.03168](https://joss.theoj.org/papers/10.21105/joss.03168)

---

## Overview

[Google Earth Engine](https://earthengine.google.com/) is a cloud-based service for 
geospatial processing of vector and raster data. The Earth Engine platform has a 
[JavaScript and a Python API](https://developers.google.com/earth-engine/guides) with 
different methods to process geospatial objects. Google Earth Engine also provides a 
[HUGE PETABYTE-SCALE CATALOG](https://developers.google.com/earth-engine/datasets/) of 
raster and vector data that users can process online (e.g. Landsat Missions Image 
Collections, Sentinel Missions Image Collections, MODIS Products Image Collections, World 
Database of Protected Areas, etc.). The eemont package extends the 
[Google Earth Engine Python API](https://developers.google.com/earth-engine/guides/python_install) 
with pre-processing and processing tools for the most used satellite platforms by adding 
utility methods for different 
[Earth Engine Objects](https://developers.google.com/earth-engine/guides/objects_methods_overview) 
that are friendly with the Python method chaining.


## Google Earth Engine Community: Developer Resources

The eemont Python package can be found in the 
[Earth Engine Community: Developer Resources](https://developers.google.com/earth-engine/tutorials/community/developer-resources) 
together with other awesome resources such as [geemap](https://geemap.org/) and 
[rgee](https://github.com/r-spatial/rgee).


## How does it work?

The eemont python package extends the following Earth Engine classes:

- [ee.Feature](https://developers.google.com/earth-engine/guides/features)
- [ee.FeatureCollection](http://developers.google.com/earth-engine/guides/feature_collections)
- [ee.Geometry](https://developers.google.com/earth-engine/guides/geometries)
- [ee.Image](https://developers.google.com/earth-engine/guides/image_overview)
- [ee.ImageCollection](https://developers.google.com/earth-engine/guides/ic_creating)
- [ee.List](https://developers.google.com/earth-engine/guides/objects_methods_overview)
- [ee.Number](https://developers.google.com/earth-engine/guides/objects_methods_overview)

New utility methods and constructors are added to above-mentioned classes in order
to create a more fluid code by being friendly with the Python method chaining. These
methods are mandatory for some pre-processing and processing tasks (e.g. clouds masking,
shadows masking, image scaling, spectral indices computation, etc.), and they are
presented as simple functions that give researchers, students and analysts the chance to
analyze data with far fewer lines of code.

Look at this simple example where a
[Sentinel-2 Surface Reflectance Image Collection](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR)
is pre-processed and processed in just one step:

```python
import ee, eemont
   
ee.Authenticate()
ee.Initialize()

point = ee.Geometry.PointFromQuery(
    'Cali, Colombia',
    user_agent = 'eemont-example'
) # Extended constructor

S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
    .filterBounds(point)
    .closest('2020-10-15') # Extended (pre-processing)
    .maskClouds(prob = 70) # Extended (pre-processing)
    .scaleAndOffset() # Extended (pre-processing)
    .spectralIndices(['NDVI','NDWI','BAIS2'])) # Extended (processing)
```

And just like that, the collection was pre-processed, processed and ready to be analyzed!

## Installation

Install the latest version from PyPI:

```
pip install eemont
```

Upgrade `eemont` by running:

```
pip install -U eemont
```

Install the latest version from conda-forge:

```
conda install -c conda-forge eemont
```

Install the latest dev version from GitHub by running:

```
pip install git+https://github.com/davemlz/eemont
```

## Features

Let's see some of the main features of eemont and how simple they are compared to the GEE
Python API original methods:

### Overloaded Operators

The following operators are overloaded: +, -, \*\, /, //, %, \**\ , <<, >>, &, \|\, <, <=,
==, !=, >, >=, -, ~. (and you can avoid the `ee.Image.expression()` method!)

<table>

<tr>
<th> GEE Python API </th>
<th> eemont-style </th>
</tr>

<tr>
<td>
  
``` python
ds = 'COPERNICUS/S2_SR'
          
S2 = (ee.ImageCollection(ds)
.first())

def scaleImage(img):
    scaling = img.select('B.*')
    x = scaling.multiply(0.0001)
    scaling = img.select(['AOT','WVP'])
    scaling = scaling.multiply(0.001)
    x = x.addBands(scaling)
    notScaling = img.select([
        'SCL',
        'TCI.*',
        'MSK.*',
        'QA.*'
    ]))
    return x.addBands(notScaling)
    
S2 = scaleImage(S2)

exp = '2.5*(N-R)/(N+(6*R)-(7.5*B)+1)'

imgDict = {
'N': S2.select('B8'),
'R': S2.select('B4'),
'B': S2.select('B2')
}

EVI = S2.expression(exp,imgDict)
```

</td>
<td>

``` python
ds = 'COPERNICUS/S2_SR'
          
S2 = (ee.ImageCollection(ds)
.first()
.scale())

N = S2.select('B8')
R = S2.select('B4')
B = S2.select('B2')

EVI = 2.5*(N-R)/(N+(6*R)-(7.5*B)+1)
```
</td>
</tr>

</table>

### Clouds and Shadows Masking

Masking clouds and shadows can be done using eemont with just one method: `maskClouds()`!

<table>

<tr>
<th> GEE Python API </th>
<th> eemont-style </th>
</tr>

<tr>
<td>
  
``` python
ds = 'LANDSAT/LC08/C01/T1_SR'
          
def maskCloudsShadows(img):
    c = (1 << 3)
    s = (1 << 5)
    qa = 'pixel_qa'
    qa = img.select(qa)
    cm = qa.bitwiseAnd(c).eq(0)
    sm = qa.bitwiseAnd(s).eq(0)
    mask = cm.And(sm)
    return img.updateMask(mask)
    
(ee.ImageCollection(ds)
    .map(maskCloudsShadows))
```

</td>
<td>

``` python
ds = 'LANDSAT/LC08/C01/T1_SR'
          
(ee.ImageCollection(ds)
    .maskClouds())
```
</td>
</tr>

</table>

### Image Scaling and Offsetting

Scaling and offsetting can also be done using eemont with just one method: `scale()`!

<table>

<tr>
<th> GEE Python API </th>
<th> eemont-style </th>
</tr>

<tr>
<td>
  
``` python
def scaleBands(img):
    scaling = img.select([
    'NDVI',
    'EVI',
    'sur.*'
    ])
    x = scaling.multiply(0.0001)
    scaling = img.select('.*th')
    scaling = scaling.multiply(0.01)
    x = x.addBands(scaling)
    notScaling = img.select([
    'DetailedQA',
    'DayOfYear',
    'SummaryQA'
    ])
    return x.addBands(notScaling)              

ds = 'MODIS/006/MOD13Q1'

(ee.ImageCollection(ds)
    .map(scaleBands))
```

</td>
<td>

``` python
ds = 'MODIS/006/MOD13Q1'
          
(ee.ImageCollection(ds)
    .scaleAndOffset())
```
</td>
</tr>

</table>


### Complete Preprocessing

The complete preprocessing workflow (Masking clouds and shadows, and image scaling and
offsetting) can be done using eemont with just one method: `preprocess()`!


<table>

<tr>
<th> GEE Python API </th>
<th> eemont-style </th>
</tr>

<tr>
<td>
  
``` python
ds = 'LANDSAT/LC08/C01/T1_SR'
          
def maskCloudsShadows(img):
    c = (1 << 3)
    s = (1 << 5)
    qa = 'pixel_qa'
    qa = img.select(qa)
    cm = qa.bitwiseAnd(c).eq(0)
    sm = qa.bitwiseAnd(s).eq(0)
    mask = cm.And(sm)
    return img.updateMask(mask)
    
def scaleBands(img):
    scaling = img.select('B[1-7]')
    x = scaling.multiply(0.0001)
    scaling = img.select([
    'B10','B11'
    ])
    scaling = scaling.multiply(0.1)
    x = x.addBands(scaling)
    notScaling = img.select([
    'sr_aerosol',
    'pixel_qa',
    'radsat_qa'
    ])
    return x.addBands(notScaling)
    
(ee.ImageCollection(ds)
    .map(maskCloudsShadows)
    .map(scaleBands))
```

</td>
<td>

``` python
ds = 'LANDSAT/LC08/C01/T1_SR'
          
(ee.ImageCollection(ds)
    .preprocess())
```
</td>
</tr>

</table>


### Spectral Indices

Do you need to compute several spectral indices? Use the `spectralIndices()` method! The
indices are taken from [Awesome Spectral Indices](https://github.com/davemlz/awesome-spectral-indices).

<table>

<tr>
<th> GEE Python API </th>
<th> eemont-style </th>
</tr>

<tr>
<td>
  
``` python
ds = 'LANDSAT/LC08/C01/T1_SR'
          
def scaleImage(img):
    scaling = img.select('B[1-7]')
    x = scaling.multiply(0.0001)
    scaling = img.select(['B10','B11'])
    scaling = scaling.multiply(0.1)
    x = x.addBands(scaling)
    notScaling = img.select([
        'sr_aerosol',
        'pixel_qa',
        'radsat_qa'
    ]))
    return x.addBands(notScaling)

def addIndices(img):
    x = ['B5','B4']
    a = img.normalizedDifference(x)
    a = a.rename('NDVI')
    x = ['B5','B3']
    b = img.normalizedDifference(x)
    b = b.rename('GNDVI')
    x = ['B3','B6']
    c = img.normalizedDifference(x)
    c = b.rename('NDSI')
    return img.addBands([a,b,c])                    

(ee.ImageCollection(ds)
    .map(scaleImage)
    .map(addIndices))
```

</td>
<td>

``` python
ds = 'LANDSAT/LC08/C01/T1_SR'
          
(ee.ImageCollection(ds)
    .scaleAndOffset()
    .spectralIndices([
        'NDVI',
        'GNDVI',
        'NDSI'])
)
```
</td>
</tr>

</table>

The list of available indices can be retrieved by running:

``` python 
eemont.listIndices()
```

Information about the indices can also be checked:

``` python 
indices = eemont.indices() 
indices.BAIS2.formula
indices.BAIS2.reference
```

### Closest Image to a Specific Date

Struggling to get the closest image to a specific date? Here is the solution: the
`closest()` method!

<table>

<tr>
<th> GEE Python API </th>
<th> eemont-style </th>
</tr>

<tr>
<td>
  
``` python
ds = 'COPERNICUS/S5P/OFFL/L3_NO2'
          
xy = [-76.21, 3.45]
poi = ee.Geometry.Point(xy)

date = ee.Date('2020-10-15')
date = date.millis()

def setTimeDelta(img):              
    prop = 'system:time_start'
    prop = img.get(prop)
    prop = ee.Number(prop)              
    delta = prop.subtract(date)
    delta = delta.abs()              
    return img.set(
    'dateDist',
    delta)                     

(ee.ImageCollection(ds)
    .filterBounds(poi)
    .map(setTimeDelta)
    .sort('dateDist')
    .first())
```

</td>
<td>

``` python
ds = 'COPERNICUS/S5P/OFFL/L3_NO2'
          
xy = [-76.21, 3.45]
poi = ee.Geometry.Point(xy)

(ee.ImageCollection(ds)
    .filterBounds(poi)
    .closest('2020-10-15'))
```
</td>
</tr>

</table>


### Time Series By Regions

The JavaScript API has a method for time series extraction (included in the `ui.Chart`
module), but this method is missing in the Python API... so, here it is!

PD: Actually, there are two methods that you can use: `getTimeSeriesByRegion()` and
`getTimeSeriesByRegions()`!

``` python
f1 = ee.Feature(ee.Geometry.Point([3.984770,48.767221]).buffer(50),{'ID':'A'})
f2 = ee.Feature(ee.Geometry.Point([4.101367,48.748076]).buffer(50),{'ID':'B'})
fc = ee.FeatureCollection([f1,f2])

S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
    .filterBounds(fc)
    .filterDate('2020-01-01','2021-01-01')
    .maskClouds()
    .scaleAndOffset()
    .spectralIndices(['EVI','NDVI']))

# By Region
ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                geometry = fc,
                                bands = ['EVI','NDVI'],
                                scale = 10)

# By Regions
ts = S2.getTimeSeriesByRegions(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                collection = fc,
                                bands = ['EVI','NDVI'],
                                scale = 10)
```

### Constructors by Queries

Don't you have the coordinates of a place? You can construct them by using queries!

``` python
usr = 'my-eemont-query-example'
   
seattle_bbox = ee.Geometry.BBoxFromQuery('Seattle',user_agent = usr)
cali_coords = ee.Feature.PointFromQuery('Cali, Colombia',user_agent = usr)
amazonas_river = ee.FeatureCollection.MultiPointFromQuery('Río Amazonas',user_agent = usr)
```

### JavaScript Modules

This is perhaps the most important feature in `eeExtra`! What if you could use a
JavaScript module (originally just useful for the Code Editor) in python or R? Well,
wait no more for it!

<table>

<tr>
<th> JS (Code Editor) </th>
<th> Python (eemont) </th>
<th> R (rgee+) </th>
</tr>

<tr>
<td>
  
``` javascript
var usr = 'users/sofiaermida/'
var rep = 'landsat_smw_lst:'
var fld = 'modules/'
var fle = 'Landsat_LST.js'
var pth = usr+rep+fld+fle
var mod = require(pth)
var LST = mod.collection(
    ee.Geometry.Rectangle([
        -8.91,
        40.0,
        -8.3,
        40.4
    ]),
    'L8',
    '2018-05-15',
    '2018-05-31',
    true
)
```

</td>
<td>
  
``` python
import ee, eemont
ee.Initialize()
usr = 'users/sofiaermida/'
rep = 'landsat_smw_lst:'
fld = 'modules/'
fle = 'Landsat_LST.js'
pth = usr+rep+fld+fle
ee.install(pth)
mod = ee.require(pth)
LST = mod.collection(
    ee.Geometry.Rectangle([
        -8.91,
        40.0,
        -8.3,
        40.4
    ]),
    'L8',
    '2018-05-15',
    '2018-05-31',
    True
)
```

</td>
<td>

``` r
library(rgee)
library(rgeeExtra)
ee_Initialize()
usr <- 'users/sofiaermida/'
rep <- 'landsat_smw_lst:'
fld <- 'modules/'
fle <- 'Landsat_LST.js'
pth <- paste0(usr,rep,fld,fle)
mod <- ee$require(pth)
LST = mod$collection(
    ee$Geometry$Rectangle(c(
        -8.91,
        40.0,
        -8.3,
        40.4
    )),
    'L8',
    '2018-05-15',
    '2018-05-31',
    TRUE
)
```
</td>
</tr>

</table>

## License

The project is licensed under the MIT license.

## How to cite

Do you like using eemont and think it is useful? Share the love by citing it!:

```
Montero, D., (2021). eemont: A Python package that extends Google Earth Engine. 
Journal of Open Source Software, 6(62), 3168, https://doi.org/10.21105/joss.03168
```
   
If required, here is the BibTex!:

```
@article{Montero2021,
    doi = {10.21105/joss.03168},
    url = {https://doi.org/10.21105/joss.03168},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {62},
    pages = {3168},
    author = {David Montero},
    title = {eemont: A Python package that extends Google Earth Engine},
    journal = {Journal of Open Source Software}
}
```

## Artists

- [David Montero Loaiza](https://github.com/davemlz): Lead Developer of eemont and eeExtra.
- [César Aybar](https://github.com/csaybar): Lead Developer of rgee and eeExtra.
- [Aaron Zuspan](https://github.com/aazuspan): Plus Codes Constructors and Methods, Panchromatic Sharpening and Histogram Matching Developer.

## Credits

Special thanks to [Justin Braaten](https://github.com/jdbcode) for featuring eemont in
tutorials and the GEE Community: Developer Resources Page, to
[César Aybar](https://github.com/csaybar) for the formidable help with Awesome Spectral
Indices and to the JOSS Review Team ([Katy Barnhart](https://github.com/kbarnhart),
[Jayaram Hariharan](https://github.com/elbeejay), [Qiusheng Wu](https://github.com/giswqs)
and [Patrick Gray](https://github.com/patrickcgray)) for the comments, suggestions and contributions!---
title: 'eemont: A Python package that extends Google Earth Engine'
tags:
  - Python
  - Google Earth Engine
  - Remote Sensing
  - GIS
  - QGIS
  - R
authors:
  - name: David Montero
    orcid: 0000-0002-9010-3286
    affiliation: "1"
affiliations:
 - name: Independent Researcher
   index: 1
date: 26 March 2021
bibliography: paper.bib
---

# Summary

[Google Earth Engine](https://earthengine.google.com/) (GEE) is a cloud-based service for geospatial processing of vector and raster data [@Gorelick2017].
The GEE platform has [JavaScript and Python APIs](https://developers.google.com/earth-engine/guides) with different methods to process geospatial objects.
GEE also provides a [multi-petabyte data catalog](https://developers.google.com/earth-engine/datasets/) of geospatial datasets and satellite imagery (e.g., Landsat, Sentinel, MODIS).
The eemont package extends the [GEE Python API](https://developers.google.com/earth-engine/guides/python_install) with pre-processing and processing tools
for the commonly used satellite imagery (e.g., Landsat, Sentinel, MODIS) by adding new methods for different
[Earth Engine Objects](https://developers.google.com/earth-engine/guides/objects_methods_overview) that are friendly with the Python method chaining. The package
can be used with [geemap](https://geemap.org/) [@Wu2020] for interactive visualization in Jupyter Notebooks, with the
[GEE Plugin for QGIS](https://gee-community.github.io/qgis-earthengine-plugin/) for processing and visualization inside [QGIS](https://www.qgis.org/es/site/) [@QGIS_software], and it can be used in R with [rgee](https://github.com/r-spatial/rgee) [@Aybar2020].

# Statement of need

The typical pre-processing and processing steps of satellite imagery are long and complex, making it challenging for analysts to move from data curation
to analysis. These steps have been simplified through eemont with simple and clearer pythonic methods by extending the main Earth Engine objects.

The eemont python package extends the following Earth Engine classes:

- [ee.Feature](https://developers.google.com/earth-engine/guides/features)
- [ee.FeatureCollection](https://developers.google.com/earth-engine/guides/feature_collections)
- [ee.Geometry](https://developers.google.com/earth-engine/guides/geometries)
- [ee.Image](https://developers.google.com/earth-engine/guides/image_overview)
- [ee.ImageCollection](https://developers.google.com/earth-engine/guides/ic_creating)

New utility methods and constructors are added to above-mentioned classes in order to create a more fluid code by being friendly with the Python method chaining.
The added methods are useful for pre-processing and processing tasks (e.g., clouds masking, shadows masking, image scaling, spectral indices computation, time series, etc.),
and they are presented as simple functions that give researchers, students and analysts the chance to process a large number of geospatial datasets with a few lines of code, 
improving code-writing and producing analysis-ready geospatial datasets.

The following script shows an example of the required code to pre-process and process the Landsat-8
Surface Reflectance Image Collection using the standard GEE Python API:

```python
import ee

ee.Authenticate()
ee.Initialize()

point = ee.Geometry.Point([-74.0592,11.3172])

def maskClouds(img):
    cloudShadowBitMask = (1 << 3)
    cloudsBitMask = (1 << 5)
    qa = img.select('pixel_qa')
    cloudShadowMask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
    cloudMask = qa.bitwiseAnd(cloudsBitMask).eq(0)
    mask = cloudShadowMask.And(cloudMask)
    return image.updateMask(mask)

def scaleImage(img):
    sc = img.select('B[1-7]').multiply(0.0001)
    sc = sc.addBands(img.select(['B10','B11']).multiply(0.1))
    sc = sc.addBands(img.select(['sr_aerosol','pixel_qa','radsat_qa']))
    return sc.copyProperties(img,img.propertyNames())

def addIndices(img):
    NDVI = img.normalizedDifference(['B5','B4']).rename('NDVI')
    imgDict = {
        'N': img.select('B5'),
        'R': img.select('B4'),
        'B': img.select('B2')
    }
    formula = '2.5 * (N - R) / (N + 6.0 * R - 7.5 * B + 1.0)'
    EVI = img.expression(formula,imgDict).rename('EVI')
    GNDVI = img.normalizedDifference(['B5','B3']).rename('GNDVI')
    return img.addBands([NDVI,EVI,GNDVI])

L8 = (ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
    .filterBounds(point)
    .map(maskClouds)
    .map(scaleImage)
    .map(addIndices)
```

The above 39 lines of code can be simplified as 9 lines of code using eemont, which supports method chaining:

```python
import ee, eemont

ee.Authenticate()
ee.Initialize()

L8 = (ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
    .maskClouds()
    .scale()
    .index(['NDVI','EVI','GNDVI']))
```

These methods support multiple datasets from the [GEE STAC](https://earthengine-stac.storage.googleapis.com/catalog/catalog.json)
and can be coupled with additional packages built on top of the GEE Python API such as [geemap](https://geemap.org/) [@Wu2020],
[geetools](https://github.com/gee-community/gee_tools) [@Principe2021] and [restee](https://kmarkert.github.io/restee/) [@Markert2021], extending and simplifying the
use of GEE.

# Google Earth Engine Community: Developer Resources

The eemont package is featured on GEE's official website ([GEE Community: Developer Resources](https://developers.google.com/earth-engine/tutorials/community/developer-resources))
together with a curated list of resources developed by the GEE developer community, and works as one of the modules that simplify workflows,
extending the GEE Python API for the community.

# Compatibility with QGIS

The eemont python package can be used in QGIS with the [GEE Plugin for QGIS](https://gee-community.github.io/qgis-earthengine-plugin/) by installing the package
using the OSGeo4W Shell:

```
py3_env
python -m pip install eemont
```

After installation, eemont can be used in the Python console inside QGIS:

```python
import ee, eemont
from ee_plugin import Map

S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
    .maskClouds()
    .scale()
    .index(['NDVI','EVI','GNDVI'])
    .first())

Map.addLayer(S2,{'min':0,'max':1,'bands':'NDVI'},'NDVI',True)
```

# Compatibility with R

The eemont python package can be used in R with [rgee](https://github.com/r-spatial/rgee) by using the [reticulate](https://rstudio.github.io/reticulate/) package [@Ushey2021].
This compatibility increases the number of researchers who can use the eemont functionalities by using a different programming language.
The following chunk shows the eemont configuration and usage for R:

```r
library(rgee)
library(reticulate)

ee_Initialize()

py_install("eemont",pip = TRUE)

eemont <- import("eemont")

point <- ee$Geometry$Point(c(-74.0592,11.3172))
L8 <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$filterBounds(point)
L8 <- L8$maskClouds()$scale()$index("NDWI")
```

# Acknowledgements

The author would like to thank the Google Earth Engine team, Justin Braaten, Qiusheng Wu, César Aybar and Gennadii Donchyts for their big contribution to the GEE community.

# References
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Please paste here the code and the error you got.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Setup (please complete the following information):**
 - OS: [e.g. iOS]
 - python version [e.g. 3.9]
 - eemont version [e.g. 0.2.5]
 - earthengine-api version [e.g. 0.1.274]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is this a Feature Request for eemont or for eeExtra?**
If your request involves using just the earthengine-api and numpy modules, please submit your feature request to eeExtra. Otherwise, if your request involves using more modules, submit your feature request using this template here in eemont.

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Do you want to work on it?**
Contributions to eemont are welcome! If you want to work on this feature request, please, let us know!

**Additional context**
Add any other context or screenshots about the feature request here.
---
title: "Using eemont and geemap in R with rgee"
output: html_notebook
---

# Using eemont and geemap in R with rgee

_Tutorial created by *David Montero Loaiza*_: [GitHub](https://github.com/davemlz) | [Twitter](https://twitter.com/dmlmont)

- GitHub Repo: [https://github.com/davemlz/eemont](https://github.com/davemlz/eemont)
- PyPI link: [https://pypi.org/project/eemont/](https://pypi.org/project/eemont/)
- Conda-forge: [https://anaconda.org/conda-forge/eemont](https://anaconda.org/conda-forge/eemont)
- Documentation: [https://eemont.readthedocs.io/](https://eemont.readthedocs.io/)
- More tutorials: [https://github.com/davemlz/eemont/tree/master/docs/tutorials](https://github.com/davemlz/eemont/tree/master/docs/tutorials)

## Let's start!

Import the required packages:

```{r}
library(rgee)
library(reticulate)
```

Initialize GEE:

```{r}
ee_Initialize()
```

If required, please uncomment:

```{r}
#py_install("eemont")
#py_install("geemap")
```

Import *eemont* and *geemap.colormaps*:

```{r}
eemont <- import("eemont")
cm <- import("geemap.colormaps")
```

Point of interest:

```{r}
point <- ee$Geometry$Point(c(-74.0592,11.3172))
```

Get and filter the Landsat 8 SR collection:

```{r}
L8 <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$filterBounds(point)
```

## Using *eemont*:

The same Python methods can be used here. Remember to change "." for "$":

```{r}
L8 <- L8$maskClouds()$scaleAndOffset()$spectralIndices("NDWI")
```

## Visualization (let's use the colormaps from *geemap*)

Let's define the visualization:

```{r}
visParamsNDWI <- list(
  min = 0,
  max = 1,
  bands = "NDWI",
  palette = cm$get_palette("Spectral_r") # Matplotlib colormap from geemap
)
```

Display results:

```{r}
Map$centerObject(point,8)
Map$addLayer(L8$first(),visParamsNDWI,"L8 NDWI")
```
Contributing
===================

Contributions to eemont are welcome! Here you will find how to do it:

- **Bugs:** If you find a bug, please report it by opening an issue. if possible, please attach the error, code, version, and other details. 

- **Fixing Issues:** If you want to contributte by fixing an issue, please   check the eemont issues: contributions are welcome for open issues with labels :code:`bug` and :code:`help wanted`.

- **Enhancement:** New features and modules are welcome! You can check the eemont issues: contributions are welcome for open issues with labels :code:`enhancement` and :code:`help wanted`.

- **Documentation:** You can add examples, notes and references to the eemont documentation by using the NumPy Docstrings of the eemont documentation, or by creating blogs, tutorials or papers.

Contribution Steps
---------------------

First, fork the `eemont <https://github.com/davemlz/eemont>`_ repository and clone it to your local machine. Then, create a development branch::

   git checkout -b name-of-dev-branch
   
eemont is divided according to Earth Engine classes, and you will find a module for each class (e.g. :code:`imagecollection.py`). Look for the required class as follows:

- ee.Feature: :code:`feature.py`
- ee.FeatureCollection: :code:`featurecollection.py`
- ee.Geometry: :code:`geometry.py`
- ee.Image: :code:`image.py`
- ee.ImageCollection: :code:`imagecollection.py`

The :code:`common.py` is used for methods that can be used for more than one Earth Engine class.

When creating new features, please start with the :code:`self` argument and add the corresponding decorator (
:code:`@extend()` from the :code:`extending` module). Check this example:

.. code-block:: python

   from .extending import extend
   
   @extend(ee.image.Image, static = False)
   def my_new_method(self,other):
        '''Returns the addition of and image and a float.
    
        Parameters
        ----------    
        self : ee.Image [this]
            Image to add.
        other : float
            Float to add.

        Returns
        -------    
        ee.Image
            Addition of an ee.Image and a float.

        Examples
        --------
        >>> import ee, eemont
        >>> ee.Initialize()
        >>> img = ee.Image(0).my_new_method(other = 3.14)
        '''
        return self.add(other)
        
By using the :code:`@extend()` decorator, the :code:`my_new_method()` method is added to the :code:`ee.Image` class. If you want to add a static method, please set the :code:`static` argument to :code:`False`. Look for the required class as follows:

- ee.Feature: :code:`ee.feature.Feature`
- ee.FeatureCollection: :code:`ee.featurecollection.FeatureCollection`
- ee.Geometry: :code:`ee.geometry.Geometry`
- ee.Image: :code:`ee.image.Image`
- ee.ImageCollection: :code:`ee.imagecollection.ImageCollection`
- ee.List: :code:`ee.ee_list.List`
- ee.Number: :code:`ee.ee_number.Number`

Remember to use `Black <https://github.com/psf/black>`_!

In order to test additions, you can use :code:`pytest` over the :code:`tests` folder::

   pytest tests
   
This will autmatically test all modules for the available satellite platforms through eemont. If you have added a new feature, please include it in the tests.

To test across different Python versions, please use :code:`tox`.

Now it's time to commit your changes and push your development branch::

   git add .
   git commit -m "Description of your work"
   git push origin name-of-dev-branch
  
And finally, submit a pull request.Contributors
============

.. list-table::
   :widths: 50 50
   
   * - .. image:: _static/davemlz.jpg
          :width: 400
          :alt: David Montero Loaiza          
             
       **David Montero Loaiza**
       
       Lead Developer of `eemont <https://github.com/davemlz/eemont>`_ and `eeExtra <https://github.com/r-earthengine/ee_extra>`_.
       
       - **GitHub**: `davemlz <https://github.com/davemlz>`_
       - **Twitter**: `dmlmont <https://twitter.com/dmlmont>`_         
       
       **Bio**: Remote Sensing Analyst. Research Assistant at `RSC4Earth <https://rsc4earth.de/authors/dmontero/>`_. PhD Student at the `University of Leipzig <https://www.physgeo.uni-leipzig.de/en/institute-of-geophysics-and-geology/research/remote-sensing-centre-for-earth-system-research-rsc4earth/>`_.
       
       **Other projects**: 
       
           - `eeExtra <https://github.com/r-earthengine/ee_extra>`_: A ninja Python package behind rgee, rgeeExtra and eemont.
           - `Awesome EE Spectral Indices <https://github.com/davemlz/awesome-ee-spectral-indices>`_: A ready-to-use curated list of spectral indices for Google Earth Engine.
           - `spectral <https://github.com/davemlz/spectral>`_: Awesome Spectral Indices for the Google Earth Engine JavaScript API (Code Editor).
       
     - .. image:: _static/csaybar.jpg
          :width: 400
          :alt: Cesar Aybar       
             
       **César Aybar**
       
       Lead Developer of `rgee <https://github.com/r-spatial/rgee>`_ and `eeExtra <https://github.com/r-earthengine/ee_extra>`_.
       
       - **GitHub**: `csaybar <https://github.com/csaybar>`_
       - **Twitter**: `csaybar <https://twitter.com/csaybar>`_         
       
       **Bio**: Master Student at the `Copernicus Master in Digital Earth <https://github.com/r-spatial/rgee>`_.
       
       **Other projects**: 
       
           - `eeExtra <https://github.com/r-earthengine/ee_extra>`_: A ninja Python package behind rgee, rgeeExtra and eemont.
           - `rgee <https://github.com/r-spatial/rgee>`_: Google Earth Engine for R.
           - `rgeeExtra <https://github.com/r-earthengine/rgeeExtra>`_: High-level functions to process spatial and simple Earth Engine objects.
   * - .. image:: _static/logo.png
          :width: 400
          :alt: Aaron Zuspan         
             
       **Aaron Zuspan**
       
       Pan-sharpening, histogram matching, and Plus Codes Developer.
       
       - **GitHub**: `aazuspan <https://github.com/aazuspan>`_
       - **Twitter**: `aazuspan <https://twitter.com/aazuspan>`_         
       
       **Bio**: Geospatial analyst and software developer, interested in using code and remote sensing to study forests and fire.
       
       **Other projects**: 
       
           - `sankee <https://github.com/aazuspan/sankee>`_: Interactive Sankey plots of land cover changes in Earth Engine.
           - `eexarray <https://github.com/aazuspan/eexarray>`_: A Python interface between Earth Engine and xarray.
           - `geeSharp.js <https://github.com/aazuspan/geeSharp.js>`_: Pan-sharpening in the Earth Engine Code Editor.
       
     -Tutorials
============

Here you can find a collection of eemont tutorials in Jupyter Notebooks and RMarkdown files:

1. Masking Clouds and Shadows in Sentinel-2 (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/001-Clouds-and-Shadows-Masking-Sentinel-2.ipynb>`_) 
2. Scaling a Sentinel-2 Image Collection (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/002-Sentinel-2-Image-Collection-Scaling.ipynb>`_) 
3. Getting the Closest Image to a Specific Date (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/003-Closest-Image-to-Date-MOD16A2.ipynb>`_) 
4. Computing Spectral Indices on Landsat 8 (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/004-Computing-Spectral-Indices-Landsat-8.ipynb>`_) 
5. Computing the EVI with Overloaded Operators (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/005-EVI-with-Overloaded-Operators-Sentinel-2.ipynb>`_) 
6. Computing NDSI and Snow Cover using Overloaded Operators and Rich Comparisons (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/006-NDSI-and-Snow-Cover-Sentinel-2-MOD10A2.ipynb>`_) 
7. Masking Clouds in Sentinel-3 (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/007-Clouds-Masking-Sentinel-3.ipynb>`_) 
8. Cloudless MOD09Q1 Median Composite (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/008-Cloudless-MOD09Q1-Median-Composite.ipynb>`_) 
9. Using eemont and geemap in R with rgee (`RMarkdown <https://github.com/davemlz/eemont/blob/master/tutorials/009-eemont-And-geemap-In-R-With-rgee.Rmd>`_) 
10. Creating Points From Queries (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/010-Creating-Points-From-Queries.ipynb>`_) 
11. Creating a Bounding Box From a Query (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/011-Creating-A-Bounding-Box-From-Query.ipynb>`_) 
12. Computing Spectral Indices on MODIS (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/012-Spectral-Indices-MODIS-MOD09GA.ipynb>`_) 
13. Time Series By Region and Conversion to Pandas (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/013-Time-Series-By-Region-Pandas.ipynb>`_) 
14. Time Series By Regions and Conversion to Pandas (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/014-Time-Series-By-Regions-Pandas.ipynb>`_)
15. Scaling and Offseting ANY RASTER DATASET From the GEE STAC (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/015-Scaling-ANY-Raster-From-GEE-STAC.ipynb>`_) 
16. Spectral Indices From the Awesome Spectral Indices for GEE (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/016-Spectral-Indices-From-Awesome-Spectral-Indices-List.ipynb>`_) 
17. Masking Clouds and Shadows in VIIRS Products (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/017-VIIRS-Products-Clouds-Masking.ipynb>`_) 
18. Complete Preprocessing (Clouds Masking, Shadows Masking, Scaling and Offsetting) With Just One Method (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/018-Complete-Preprocessing-One-Method.ipynb>`_) 
19. Checking the STAC Info of ANY RASTER DATASET in the GEE STAC (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/019-Checking-STAC-Info.ipynb>`_) 
20. Overloaded Operators for the ee.Number Object Class (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/020-Overloaded-Operators-Number.ipynb>`_) 
21. Citation Tools for ANY RASTER DATASET in the GEE STAC (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/021-Citation-Tools.ipynb>`_) 
22. Overloaded Operators for the ee.List Object Class (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/022-Overloaded-Operators-List.ipynb>`_)
23. Creating Geometries from Plus Codes (And Viceversa) (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/023-Creating-Geometries-From-Plus-Codes.ipynb>`_)
24. Container Emulation Methods for ee.Image and ee.ImageCollection (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/024-Container-Image-ImageCollection.ipynb>`_)
25. Landsat 8 - Collection 2 (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/025-Landsat-8-Collection-2.ipynb>`_)
26. Histogram Matching (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/026-Histogram-Matching.ipynb>`_)
27. Panchromatic Sharpening (`Jupyter Notebook <https://github.com/davemlz/eemont/blob/master/tutorials/027-Panchromatic-Sharpening.ipynb>`_)Tutorials
============

.. toctree::
    :caption: Clouds and Shadows Masking
    :maxdepth: 3

    tutorials/001-Clouds-and-Shadows-Masking-Sentinel-2.ipynb
    tutorials/008-Cloudless-MOD09Q1-Median-Composite.ipynb
    tutorials/007-Clouds-Masking-Sentinel-3.ipynb
    tutorials/017-VIIRS-Products-Clouds-Masking.ipynb
    tutorials/018-Complete-Preprocessing-One-Method.ipynb    

.. toctree::
    :caption: Scaling and Offsetting
    :maxdepth: 3

    tutorials/002-Sentinel-2-Image-Collection-Scaling.ipynb
    tutorials/015-Scaling-ANY-Raster-From-GEE-STAC.ipynb
    tutorials/025-Landsat-8-Collection-2.ipynb
    tutorials/029-Landsat457-Collection-2.ipynb

.. toctree::
    :caption: Spectral Transformations
    :maxdepth: 3

    tutorials/004-Computing-Spectral-Indices-Landsat-8.ipynb
    tutorials/012-Spectral-Indices-MODIS-MOD09GA.ipynb
    tutorials/016-Spectral-Indices-From-Awesome-Spectral-Indices-List.ipynb
    tutorials/028-Tasseled-Cap.ipynb
    tutorials/030-Awesome-Spectral-Indices-v003.ipynb
    tutorials/032-Combining-eemont-wxee.ipynb

.. toctree::
    :caption: Time Series
    :maxdepth: 3

    tutorials/003-Closest-Image-to-Date-MOD16A2.ipynb
    tutorials/013-Time-Series-By-Region-Pandas.ipynb
    tutorials/014-Time-Series-By-Regions-Pandas.ipynb

.. toctree::
    :caption: Histogram Matching and Pan Sharpening
    :maxdepth: 3

    tutorials/026-Histogram-Matching.ipynb
    tutorials/027-Panchromatic-Sharpening.ipynb

.. toctree::
    :caption: Overloaded Operators and Container Methods
    :maxdepth: 3

    tutorials/005-EVI-with-Overloaded-Operators-Sentinel-2.ipynb
    tutorials/006-NDSI-and-Snow-Cover-Sentinel-2-MOD10A2.ipynb
    tutorials/020-Overloaded-Operators-Number.ipynb
    tutorials/022-Overloaded-Operators-List.ipynb
    tutorials/024-Container-Image-ImageCollection.ipynb

.. toctree::
    :caption: Constructors
    :maxdepth: 3

    tutorials/010-Creating-Points-From-Queries.ipynb
    tutorials/011-Creating-A-Bounding-Box-From-Query.ipynb
    tutorials/023-Creating-Geometries-From-Plus-Codes.ipynb

.. toctree::
    :caption: Earth Engine Apps
    :maxdepth: 3

    tutorials/031-App-Manager.ipynb
    tutorials/033-Earth-Engine-Apps.ipynb

.. toctree::
    :caption: Other STAC Methods
    :maxdepth: 3

    tutorials/019-Checking-STAC-Info.ipynb    
    tutorials/021-Citation-Tools.ipynbWelcome to eemont!
==================

.. raw:: html

   <embed>
     <p align="center">
       <a href="https://github.com/davemlz/eemont"><img src="https://raw.githubusercontent.com/davemlz/davemlz/main/eemont.png" height="200px"/></a>       
       <br>
       <b>A Python package that extends <a href="https://earthengine.google.com/">Google Earth Engine</a></b>
     </p>
   </embed>

.. image:: https://img.shields.io/pypi/v/eemont.svg
        :target: https://pypi.python.org/pypi/eemont
        :alt: PyPI Version
        
.. image:: https://img.shields.io/conda/vn/conda-forge/eemont.svg
        :target: https://anaconda.org/conda-forge/eemont
        :alt: conda-forge Version
        
.. image:: https://img.shields.io/badge/License-MIT-blue.svg
        :target: https://opensource.org/licenses/MIT
        :alt: License
        
.. image:: https://readthedocs.org/projects/eemont/badge/?version=latest
        :target: https://eemont.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://github.com/davemlz/eemont/actions/workflows/tests.yml/badge.svg
        :target: https://github.com/davemlz/eemont/actions/workflows/tests.yml
        :alt: Tests Status  

.. image:: https://img.shields.io/badge/GitHub%20Sponsors-Donate-ff69b4.svg
        :target: https://github.com/sponsors/davemlz
        :alt: Donation GitHub Sponsors

.. image:: https://img.shields.io/badge/Buy%20me%20a%20coffee-Donate-ff69b4.svg
        :target: https://www.buymeacoffee.com/davemlz
        :alt: Donation Buy me a coffee
        
.. image:: https://img.shields.io/badge/kofi-Donate-ff69b4.svg
        :target: https://ko-fi.com/davemlz
        :alt: Donation kofi
        
.. image:: https://static.pepy.tech/personalized-badge/eemont?period=total&units=international_system&left_color=grey&right_color=green&left_text=Downloads
        :target: https://pepy.tech/project/eemont
        :alt: Downloads
        
.. image:: https://img.shields.io/badge/GEE%20Community-Developer%20Resources-00b6ff.svg
        :target: https://developers.google.com/earth-engine/tutorials/community/developer-resources
        :alt: GEE Community
        
.. image:: https://img.shields.io/twitter/follow/dmlmont?style=social
        :target: https://twitter.com/dmlmont
        :alt: Twitter
        
.. image:: https://joss.theoj.org/papers/34696c5b62c50898b4129cbbe8befb0a/status.svg
        :target: https://joss.theoj.org/papers/34696c5b62c50898b4129cbbe8befb0a
        :alt: JOSS Paper
        
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/psf/black
        :alt: Black

.. image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
        :target: https://pycqa.github.io/isort/
        :alt: isort


.. toctree::   
   :maxdepth: 2
   :caption: Extended Classes 
   :hidden:      
   
   classes/index
   
.. toctree::
   :maxdepth: 2
   :caption: User Guide      
   :hidden:
      
   guide/index
   
.. toctree::
   :maxdepth: 2
   :caption: Extensions    
   :hidden:
   
   guide/eemontR
   
.. toctree::
   :maxdepth: 2
   :caption: What's new? 
   :hidden:
      
   changelog
   
.. toctree::
   :maxdepth: 2
   :caption: Tutorials 
   :hidden:
      
   tutorials
   
.. toctree::
   :maxdepth: 2
   :caption: Contributors 
   :hidden:
      
   contributors

- GitHub: `https://github.com/davemlz/eemont <https://github.com/davemlz/eemont>`_
- Documentation: `https://eemont.readthedocs.io/ <https://eemont.readthedocs.io/>`_
- PyPI: `https://pypi.org/project/eemont/ <https://pypi.org/project/eemont/>`_
- Conda-Forge: `https://anaconda.org/conda-forge/eemont <https://anaconda.org/conda-forge/eemont>`_
- Tutorials: `https://github.com/davemlz/eemont/tree/master/docs/tutorials <https://github.com/davemlz/eemont/tree/master/docs/tutorials>`_

**Table of Contents**

- `Overview`_
- `Google Earth Engine Community: Developer Resources`_
- `Additional Resources`_
- `How does it work?`_
- `Installation`_
- `Features`_
- `Supported Platforms`_
- `License`_
- `Contributing`_
- `How to cite`_
- `Artists`_
- `Credits`_

Overview
-------------------

`Google Earth Engine <https://earthengine.google.com/>`_ is a cloud-based service for geospatial processing of vector and raster data. The Earth Engine platform has a `JavaScript and a Python API <https://developers.google.com/earth-engine/guides>`_ with different methods to process geospatial objects. Google Earth Engine also provides a `HUGE PETABYTE-SCALE CATALOG <https://developers.google.com/earth-engine/datasets/>`_ of raster and vector data that users can process online (e.g. Landsat Missions Image Collections, Sentinel Missions Image Collections, MODIS Products Image Collections, World Database of Protected Areas, etc.). The eemont package extends the `Google Earth Engine Python API <https://developers.google.com/earth-engine/guides/python_install>`_ with pre-processing and processing tools for the most used satellite platforms by adding utility methods for different `Earth Engine Objects <https://developers.google.com/earth-engine/guides/objects_methods_overview>`_ that are friendly with the Python method chaining.

Google Earth Engine Community: Developer Resources
-----------------------------------------------------

The eemont Python package can be found in the `Earth Engine Community: Developer Resources <https://developers.google.com/earth-engine/tutorials/community/developer-resources>`_ together with other awesome resources such as `geemap <https://geemap.org/>`_ and `rgee <https://github.com/r-spatial/rgee>`_.

Additional Resources
--------------------

If you like eemont, you might be interested in...

- `Awesome Spectral Indices for GEE <https://github.com/davemlz/awesome-ee-spectral-indices>`_: A ready-to-use curated list of spectral indices for Google Earth Engine.
- `spectral <https://github.com/davemlz/spectral>`_: Awesome Spectral Indices for the Google Earth Engine JavaScript API (Code Editor).
- `eeExtra <https://github.com/r-earthengine/ee_extra>`_: A ninja Python package behind rgee, rgeeExtra and eemont.
- `rgeeExtra <https://github.com/r-earthengine/rgeeExtra>`_: High-level functions to process spatial and simple Earth Engine objects.

How does it work?
-------------------

The eemont python package extends the following Earth Engine classes:

- `ee.Feature <https://developers.google.com/earth-engine/guides/features>`_
- `ee.FeatureCollection <https://developers.google.com/earth-engine/guides/feature_collections>`_
- `ee.Geometry <https://developers.google.com/earth-engine/guides/geometries>`_
- `ee.Image <https://developers.google.com/earth-engine/guides/image_overview>`_
- `ee.ImageCollection <https://developers.google.com/earth-engine/guides/ic_creating>`_
- `ee.List <https://developers.google.com/earth-engine/guides/objects_methods_overview>`_
- `ee.Number <https://developers.google.com/earth-engine/guides/objects_methods_overview>`_

New utility methods and constructors are added to above-mentioned classes in order to create a more fluid code by being friendly with the Python method chaining. These methods are mandatory for some pre-processing and processing tasks (e.g. clouds masking, shadows masking, image scaling, spectral indices computation, etc.), and they are presented as simple functions that give researchers, students and analysts the chance to analyze data with far fewer lines of code.

Look at this simple example where a `Sentinel-2 Surface Reflectance Image Collection <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR>`_ is pre-processed and processed in just one step:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()
   
   point = ee.Geometry.PointFromQuery('Cali, Colombia',user_agent = 'eemont-example') # Extended constructor
   
   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
       .filterBounds(point)
       .closest('2020-10-15') # Extended (pre-processing)
       .maskClouds(prob = 70) # Extended (pre-processing)
       .scaleAndOffset() # Extended (pre-processing)
       .spectralIndices(['NDVI','NDWI','BAIS2'])) # Extended (processing)

And just like that, the collection was pre-processed, processed and ready to be analyzed!

Installation
------------

Install the latest eemont version from PyPI by running:

.. code-block::   
      
   pip install eemont

Upgrade eemont by running:

.. code-block::   
      
   pip install -U eemont

Install the development version from GitHub by running:

.. code-block::   
      
   pip install git+https://github.com/davemlz/eemont
   
Install the latest eemont version from conda-forge by running:

.. code-block::   
      
   conda install -c conda-forge eemont

Features
--------

Let's see some of the main features of eemont and how simple they are compared to the GEE Python API original methods:

Overloaded Operators
~~~~~~~~~~~~~~~~~~~~~~~

The following operators are overloaded: +, -, \*\, /, //, %, \**\ , <<, >>, &, \|\, <, <=, ==, !=, >, >=, -, ~. (and you can avoid the :code:`ee.Image.expression()` method!)

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - GEE Python API
     - eemont-style     
   * - .. code-block:: python             
          
          ds = 'COPERNICUS/S2_SR'
          
          S2 = (ee.ImageCollection(ds)
            .first())
          
          def scaleImage(img):
              scaling = img.select('B.*')
              x = scaling.multiply(0.0001)
              scaling = img.select(['AOT','WVP'])
              scaling = scaling.multiply(0.001)
              x = x.addBands(scaling)
              notScaling = img.select([
                  'SCL',
                  'TCI.*',
                  'MSK.*',
                  'QA.*'
              ]))
              return x.addBands(notScaling)
              
          S2 = scaleImage(S2)
          
          exp = '2.5*(N-R)/(N+(6*R)-(7.5*B)+1)'
          
          imgDict = {
            'N': S2.select('B8'),
            'R': S2.select('B4'),
            'B': S2.select('B2')
          }
   
          EVI = S2.expression(exp,imgDict)
     - .. code-block:: python                     
   
          ds = 'COPERNICUS/S2_SR'
          
          S2 = (ee.ImageCollection(ds)
            .first()
            .scale())

          N = S2.select('B8')
          R = S2.select('B4')
          B = S2.select('B2')

          EVI = 2.5*(N-R)/(N+(6*R)-(7.5*B)+1)

Clouds and Shadows Masking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Masking clouds and shadows can be done using eemont with just one method: :code:`maskClouds()`!

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - GEE Python API
     - eemont-style     
   * - .. code-block:: python   
          :linenos:
          
          ds = 'LANDSAT/LC08/C01/T1_SR'
          
          def maskCloudsShadows(img):
              c = (1 << 3)
              s = (1 << 5)
              qa = 'pixel_qa'
              qa = img.select(qa)
              cm = qa.bitwiseAnd(c).eq(0)
              sm = qa.bitwiseAnd(s).eq(0)
              mask = cm.And(sm)
              return img.updateMask(mask)
              
          (ee.ImageCollection(ds)
            .map(maskCloudsShadows))
     - .. code-block:: python 
          :linenos:          
   
          ds = 'LANDSAT/LC08/C01/T1_SR'
          
          (ee.ImageCollection(ds)
            .maskClouds())

Image Scaling and Offsetting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scaling and offsetting can also be done using eemont with just one method: :code:`scale()`!

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - GEE Python API
     - eemont-style     
   * - .. code-block:: python
          :linenos:          
   
          def scaleBands(img):
              scaling = img.select([
                'NDVI',
                'EVI',
                'sur.*'
              ])
              x = scaling.multiply(0.0001)
              scaling = img.select('.*th')
              scaling = scaling.multiply(0.01)
              x = x.addBands(scaling)
              notScaling = img.select([
                'DetailedQA',
                'DayOfYear',
                'SummaryQA'
              ])
              return x.addBands(notScaling)              
          
          ds = 'MODIS/006/MOD13Q1'
          
          (ee.ImageCollection(ds)
            .map(scaleBands))
     - .. code-block:: python
          :linenos:          
   
          ds = 'MODIS/006/MOD13Q1'
          
          (ee.ImageCollection(ds)
            .scale())
            
Complete Preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The complete preprocessing workflow (Masking clouds and shadows, and image scaling and offsetting) can be done using eemont with just one method: :code:`preprocess()`!

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - GEE Python API
     - eemont-style     
   * - .. code-block:: python   
          :linenos:
          
          ds = 'LANDSAT/LC08/C01/T1_SR'
          
          def maskCloudsShadows(img):
              c = (1 << 3)
              s = (1 << 5)
              qa = 'pixel_qa'
              qa = img.select(qa)
              cm = qa.bitwiseAnd(c).eq(0)
              sm = qa.bitwiseAnd(s).eq(0)
              mask = cm.And(sm)
              return img.updateMask(mask)
              
          def scaleBands(img):
              scaling = img.select('B[1-7]')
              x = scaling.multiply(0.0001)
              scaling = img.select([
                'B10','B11'
              ])
              scaling = scaling.multiply(0.1)
              x = x.addBands(scaling)
              notScaling = img.select([
                'sr_aerosol',
                'pixel_qa',
                'radsat_qa'
              ])
              return x.addBands(notScaling)
              
          (ee.ImageCollection(ds)
            .map(maskCloudsShadows)
            .map(scaleBands))
     - .. code-block:: python 
          :linenos:          
   
          ds = 'LANDSAT/LC08/C01/T1_SR'
          
          (ee.ImageCollection(ds)
            .preprocess())

Spectral Indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Do you need to compute several spectral indices? Use the :code:`index()` method! A lot of built-in vegetation, burn, water, snow, drought and kernel indices can be computed:

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - GEE Python API
     - eemont-style     
   * - .. code-block:: python                    
   
          ds = 'LANDSAT/LC08/C01/T1_SR'
          
          def scaleImage(img):
              scaling = img.select('B[1-7]')
              x = scaling.multiply(0.0001)
              scaling = img.select(['B10','B11'])
              scaling = scaling.multiply(0.1)
              x = x.addBands(scaling)
              notScaling = img.select([
                  'sr_aerosol',
                  'pixel_qa',
                  'radsat_qa'
              ]))
              return x.addBands(notScaling)
          
          def addIndices(img):
              x = ['B5','B4']
              a = img.normalizedDifference(x)
              a = a.rename('NDVI')
              x = ['B5','B3']
              b = img.normalizedDifference(x)
              b = b.rename('GNDVI')
              x = ['B3','B6']
              c = img.normalizedDifference(x)
              c = b.rename('NDSI')
              return img.addBands([a,b,c])                    
          
          (ee.ImageCollection(ds)
            .map(scaleImage)
            .map(addIndices))
          
     - .. code-block:: python                 
   
          ds = 'LANDSAT/LC08/C01/T1_SR'
          
          (ee.ImageCollection(ds)
            .scale()
            .index(['NDVI','GNDVI','NDSI']))

The list of available indices can be retrieved by running:

.. code-block:: python  
   
   eemont.listIndices()

Information about the indices can also be checked:

.. code-block:: python   
       
   indices = eemont.indices() 
   indices.BAIS2.formula
   indices.BAIS2.reference

Closest Image to a Specific Date
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Struggling to get the closest image to a specific date? Here is the solution: the :code:`closest()` method!

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - GEE Python API
     - eemont-style     
   * - .. code-block:: python
          :linenos:          
   
          ds = 'COPERNICUS/S5P/OFFL/L3_NO2'
          
          xy = [-76.21, 3.45]
          poi = ee.Geometry.Point(xy)
          
          date = ee.Date('2020-10-15')
          date = date.millis()
          
          def setTimeDelta(img):              
              prop = 'system:time_start'
              prop = img.get(prop)
              prop = ee.Number(prop)              
              delta = prop.subtract(date)
              delta = delta.abs()              
              return img.set(
                'dateDist',
                delta)                     
          
          (ee.ImageCollection(ds)
            .filterBounds(poi)
            .map(setTimeDelta)
            .sort('dateDist')
            .first())
          
     - .. code-block:: python
          :linenos:          
   
          ds = 'COPERNICUS/S5P/OFFL/L3_NO2'
          
          xy = [-76.21, 3.45]
          poi = ee.Geometry.Point(xy)
          
          (ee.ImageCollection(ds)
            .filterBounds(poi)
            .closest('2020-10-15'))
       
Time Series By Regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The JavaScript API has a method for time series extraction (included in the ui.Chart module), but this method is missing in the Python API... so, here it is!

PD: Actually, there are two methods that you can use: :code:`getTimeSeriesByRegion()` and :code:`getTimeSeriesByRegions()`!

.. code-block:: python

   f1 = ee.Feature(ee.Geometry.Point([3.984770,48.767221]).buffer(50),{'ID':'A'})
   f2 = ee.Feature(ee.Geometry.Point([4.101367,48.748076]).buffer(50),{'ID':'B'})
   fc = ee.FeatureCollection([f1,f2])

   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
      .filterBounds(fc)
      .filterDate('2020-01-01','2021-01-01')
      .maskClouds()
      .scale()
      .index(['EVI','NDVI']))

   # By Region
   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10)
   
   # By Regions
   ts = S2.getTimeSeriesByRegions(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                  collection = fc,
                                  bands = ['EVI','NDVI'],
                                  scale = 10)
                                  
Constructors by Queries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Don't you have the coordinates of a place? You can construct them by using queries!

.. code-block:: python

   usr = 'my-eemont-query-example'
   
   seattle_bbox = ee.Geometry.BBoxFromQuery('Seattle',user_agent = usr)
   cali_coords = ee.Feature.PointFromQuery('Cali, Colombia',user_agent = usr)
   amazonas_river = ee.FeatureCollection.MultiPointFromQuery('Río Amazonas',user_agent = usr)

Supported Platforms
------------------------

The Supported Platforms for each method can be found in the eemont documentation.

- Masking clouds and shadows supports Sentinel Missions (Sentinel-2 SR and Sentinel-3), Landsat Missions (SR products) and some MODIS Products. Check all details in User Guide > Masking Clouds and Shadows > Supported Platforms.
- Image scaling supports Sentinel Missions (Sentinel-2 and Sentinel-3), Landsat Missions and most MODIS Products. Check all details in User Guide > Image Scaling > Supported Platforms.
- Spectral indices computation supports Sentinel-2 and Landsat Missions. Check all details in User Guide > Spectral Indices > Supported Platforms.
- Getting the closest image to a specific date and time series supports all image collections with the :code:`system:time_start` property.

License
-------

The project is licensed under the MIT license.

Contributing
------------------

Contributions to eemont are welcome! Here you will find how to do it:

- **Bugs:** If you find a bug, please report it by opening an issue. if possible, please attach the error, code, version, and other details. 

- **Fixing Issues:** If you want to contributte by fixing an issue, please   check the eemont issues: contributions are welcome for open issues with labels :code:`bug` and :code:`help wanted`.

- **Enhancement:** New features and modules are welcome! You can check the eemont issues: contributions are welcome for open issues with labels :code:`enhancement` and :code:`help wanted`.

- **Documentation:** You can add examples, notes and references to the eemont documentation by using the NumPy Docstrings of the eemont documentation, or by creating blogs, tutorials or papers.

Contribution Steps
~~~~~~~~~~~~~~~~~~~~~~~~

First, fork the `eemont <https://github.com/davemlz/eemont>`_ repository and clone it to your local machine. Then, create a development branch::

   git checkout -b name-of-dev-branch
   
eemont is divided according to Earth Engine classes, and you will find a module for each class (e.g. :code:`imagecollection.py`). Look for the required class as follows:

- ee.Feature: :code:`feature.py`
- ee.FeatureCollection: :code:`featurecollection.py`
- ee.Geometry: :code:`geometry.py`
- ee.Image: :code:`image.py`
- ee.ImageCollection: :code:`imagecollection.py`

The :code:`common.py` is used for methods that can be used for more than one Earth Engine class.

When creating new features, please start with the :code:`self` argument and add the corresponding decorator (
:code:`@extend()` from the :code:`extending` module). Check this example:

.. code-block:: python

   from .extending import extend
   
   @extend(ee.image.Image, static = False)
   def my_new_method(self,other):
        '''Returns the addition of and image and a float.
    
        Parameters
        ----------    
        self : ee.Image [this]
            Image to add.
        other : float
            Float to add.

        Returns
        -------    
        ee.Image
            Addition of an ee.Image and a float.

        Examples
        --------
        >>> import ee, eemont
        >>> ee.Initialize()
        >>> img = ee.Image(0).my_new_method(other = 3.14)
        '''
        return self.add(other)
        
By using the :code:`@extend()` decorator, the :code:`my_new_method()` method is added to the :code:`ee.Image` class. If you want to add a static method, please set the :code:`static` argument to :code:`False`. Look for the required class as follows:

- ee.Feature: :code:`ee.feature.Feature`
- ee.FeatureCollection: :code:`ee.featurecollection.FeatureCollection`
- ee.Geometry: :code:`ee.geometry.Geometry`
- ee.Image: :code:`ee.image.Image`
- ee.ImageCollection: :code:`ee.imagecollection.ImageCollection`
- ee.List: :code:`ee.ee_list.List`
- ee.Number: :code:`ee.ee_number.Number`

Remember to use `Black <https://github.com/psf/black>`_!

In order to test additions, you can use :code:`pytest` over the :code:`tests` folder::

   pytest tests
   
This will autmatically test all modules for the available satellite platforms through eemont. If you have added a new feature, please include it in the tests.

To test across different Python versions, please use :code:`tox`.

Now it's time to commit your changes and push your development branch::

   git add .
   git commit -m "Description of your work"
   git push origin name-of-dev-branch
  
And finally, submit a pull request.

How to cite
-----------

Do you like using eemont and think it is useful? Share the love by citing it!::

   Montero, D., (2021). eemont: A Python package that extends Google Earth Engine. Journal of Open Source Software, 6(62), 3168, https://doi.org/10.21105/joss.03168
   
If required, here is the BibTex!::

   @article{Montero2021,
     doi = {10.21105/joss.03168},
     url = {https://doi.org/10.21105/joss.03168},
     year = {2021},
     publisher = {The Open Journal},
     volume = {6},
     number = {62},
     pages = {3168},
     author = {David Montero},
     title = {eemont: A Python package that extends Google Earth Engine},
     journal = {Journal of Open Source Software}
   }

Artists
-------

- `David Montero Loaiza <https://github.com/davemlz>`_: Lead Developer of eemont and eeExtra.
- `César Aybar <https://github.com/csaybar>`_: Lead Developer of rgee and eeExtra.
- `Aaron Zuspan <https://github.com/aazuspan>`_: Plus Codes Constructors and Methods, Panchromatic Sharpening and Histogram Matching Developer.

Credits
-------

Special thanks to `Justin Braaten <https://github.com/jdbcode>`_ for featuring eemont in tutorials and the GEE Community: Developer Resources Page, to `César Aybar <https://github.com/csaybar>`_ for the formidable help with Awesome Spectral Indices for GEE and to the JOSS Review Team (`Katy Barnhart <https://github.com/kbarnhart>`_, `Jayaram Hariharan <https://github.com/elbeejay>`_, `Qiusheng Wu <https://github.com/giswqs>`_ and `Patrick Gray <https://github.com/patrickcgray>`_) for the comments, suggestions and contributions!Changelog
============

v0.3.0
--------------

New Features
~~~~~~~~~~~~~~~~~~~~~~

..
   - The :code:`require()` extended method for ee was created.
   - The :code:`install()` extended method for ee was created.
   - The :code:`uninstall()` extended method for ee was created.
- The :code:`tasseledCap()` extended method for ee.Image and ee.ImageCollection was created.
- The :code:`ee.App` class was created.
- The :code:`listApps()` extended method for ee was created.
- The :code:`listDatasets()` extended method for ee was created.

Improvements
~~~~~~~~~~~~~~~~~~~~~~

- The Awesome Spectral Indices list was updated to v0.0.3.
- Errors (`#43 <https://github.com/davemlz/eemont/issues/43>`_) of the :code:`getTimeSeriesByRegion` method of the *ee.ImageCollection* module were solved.
- Errors (`#41 <https://github.com/davemlz/eemont/issues/41>`_) of the :code:`getTimeSeriesByRegions` method of the *ee.ImageCollection* module were solved.
- The :code:`maskClouds()`, :code:`spectralIndices()`, :code:`scaleAndOffset()` and :code:`preprocess()` extended methods for ee.Image and ee.ImageCollection classes now support the following platforms:
   
   - `USGS Landsat 5 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C02_T1_L2>`_
   - `USGS Landsat 4 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2>`_

Deprecation
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`index` method was deprecated for ee.Image and ee.ImageCollection.
- The :code:`scale` method was deprecated for ee.Image and ee.ImageCollection.

Contributors
~~~~~~~~~~~~~~~~~~~~~~

- `César Aybar <https://github.com/csaybar>`_
- `Aaron Zuspan <https://github.com/aazuspan>`_

v0.2.5
--------------

New Features
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`panSharpen()` extended method for ee.Image and ee.ImageCollection classes was created.
- The :code:`matchHistogram()` extended method for ee.Image classes was created.

Contributors
~~~~~~~~~~~~~~~~~~~~~~

- `Aaron Zuspan <https://github.com/aazuspan>`_

v0.2.4
--------------

Improvements
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`maskClouds()`, :code:`spectralIndices()`, :code:`scaleAndOffset()` and :code:`preprocess()` extended methods for ee.Image and ee.ImageCollection classes now support the following platforms:
   
   - `USGS Landsat 8 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2>`_
   - `USGS Landsat 7 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2>`_

v0.2.3
--------------

Improvements
~~~~~~~~~~~~~~~~~~~~~~

- Conflicts of the :code:`__len__` Container Emulation method of the *ee.Dictionary* module were solved.

v0.2.2
--------------

New Modules
~~~~~~~~~~~~~~~~~~~~~~

- The *ee.Dictionary* module was created.

New Features
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`__contains__` container emulation method was overloaded for ee.Dictionary and ee.List objects.
- The :code:`__len__` container emulation method was overloaded for ee.Dictionary, ee.List, ee.FeatureCollection and ee.ImageCollection objects.
- The :code:`__getitem__` container emulation method was overloaded for ee.Dictionary, ee.List, ee.Feature, ee.FeatureCollection, ee.Image and ee.ImageCollection objects.

v0.2.1
--------------

New Features
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`LinearRingFromPlusCodes()` extended constructor for ee.Gemoetry classes was created.
- The :code:`LineStringFromPlusCodes()` extended constructor for ee.Gemoetry classes was created.
- The :code:`MultiLineStringFromPlusCodes()` extended constructor for ee.Gemoetry classes was created.
- The :code:`MultiPointFromPlusCodes()` extended constructor for ee.Gemoetry classes was created.
- The :code:`MultiPolygonFromPlusCodes()` extended constructor for ee.Gemoetry classes was created.
- The :code:`PointFromPlusCode()` extended constructor for ee.Gemoetry classes was created.
- The :code:`PolygonFromPlusCodes()` extended constructor for ee.Gemoetry classes was created.
- The :code:`RectangleFromPlusCodes()` extended constructor for ee.Gemoetry classes was created.
- The :code:`plusCodes()` extended method for ee.Gemoetry and ee.Feature classes was created.

Contributors
~~~~~~~~~~~~~~~~~~~~~~

- `Aaron Zuspan <https://github.com/aazuspan>`_

v0.2.0
--------------

New Modules
~~~~~~~~~~~~~~~~~~~~~~

- The *ee.Number* module was created.
- The *ee.List* module was created.
- The *extending* module was created.

New Features
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`getOffsetParams()` extended method for ee.Image and ee.ImageCollection classes was created.
- The :code:`getScaleParams()` extended method for ee.Image and ee.ImageCollection classes was created.
- The :code:`scaleAndOffset()` extended method for ee.Image and ee.ImageCollection classes was created and will replace the :code:`scale()` method.
- The :code:`spectralIndices()` extended method for ee.Image and ee.ImageCollection classes was created and will replace the :code:`index()` method.
- The :code:`preprocess()` extended method for ee.Image and ee.ImageCollection classes was created.
- The :code:`getDOI()` extended method for ee.Image and ee.ImageCollection classes was created.
- The :code:`getCitation()` extended method for ee.Image and ee.ImageCollection classes was created.
- The :code:`getSTAC()` extended method for ee.Image and ee.ImageCollection classes was created.
- The binary operators (+, -, \*\, /, //, %, \**\ , <<, >>, &, |) were overloaded for ee.Number objects.
- The binary operators (+, \*\) were overloaded for ee.List objects.
- The rich comparisons (<, <=, ==, !=, >, >=) were overloaded for ee.Number objects.
- The unary operators (-, ~) were overloaded for ee.Number objects.

Improvements
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`maskClouds()` extended method for ee.Image and ee.ImageCollection classes now supports the following platforms:
   
   - `VNP09GA: VIIRS Surface Reflectance Daily 500m and 1km <https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_001_VNP09GA?hl=en>`_
   - `VNP13A1: VIIRS Vegetation Indices 16-Day 500m <https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_001_VNP13A1?hl=en>`_
- The :code:`scaleAndOffset()` extended method for ee.Image and ee.ImageCollection classes now supports ALL raster datasets from the `Google Earth Engine STAC Catalog <https://developers.google.com/earth-engine/datasets>`_.
- The :code:`spectralIndices()` extended method for ee.Image and ee.ImageCollection classes now supports ALL indices from the `Awesome List of Spectral Indices for Google Earth Engine <https://github.com/davemlz/awesome-ee-spectral-indices>`_.
   
Pending Deprecation
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`index()` method for ee.Image and ee.ImageCollection classes will be deprecated in future versions.
- The :code:`scale()` method for ee.Image and ee.ImageCollection classes will be deprecated in future versions.

v0.1.9
--------------

Improvements
~~~~~~~~~~~~~~~~~~~~~~

- :code:`kernel`, :code:`sigma`, :code:`p` and :code:`c` parameters were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection classes.
- The following vegetation indices were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection:

   - 'GARI' : Green Atmospherically Resistant Vegetation Index.
   - 'GEMI' : Global Environment Monitoring Index.
   - 'GLI' : Green Leaf Index.
   - 'GVMI' : Global Vegetation Moisture Index.
   - 'VARI' : Visible Atmospherically Resistant Index.
- The following drought indices were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection:

   - 'NDDI' : Normalized Difference Drought Index.
- The following kernel indices were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection:

   - 'kEVI' : Kernel Enhanced Vegetation Index.
   - 'kNDVI' : Kernel Normalized Difference Vegetation Index.
   - 'kRVI' : Kernel Ratio Vegetation Index.
   - 'kVARI' : Kernel Visible Atmospherically Resistant Index.

v0.1.8
--------------

New Modules
~~~~~~~~~~~~~~~~~~~~~~

- The *ee.Feature* module was created.
- The *ee.FeatureCollection* module was created.
- The *ee.Geometry* module was created.

New Features
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`getTimeSeriesByRegion()` extended method for ee.ImageCollection classes was created.
- The :code:`getTimeSeriesByRegions()` extended method for ee.ImageCollection classes was created.
- The :code:`indices()` function was created.
- The :code:`listIndices()` function was created.
- The :code:`BBoxFromQuery()` extended constructor for ee.Geometry and ee.Feature classes was created.
- The :code:`PointFromQuery()` extended constructor for ee.Geometry and ee.Feature classes was created.
- The :code:`MultiPointFromQuery()` extended constructor for ee.Geometry and ee.FeatureCollection classes was created.


Improvements
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`index()` extended method for ee.Image and ee.ImageCollection classes now supports the following platforms:
   
   - `MCD43A4.006 MODIS Nadir BRDF-Adjusted Reflectance Daily 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD43A4>`_
   - `MOD09GQ.006 Terra Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GQ>`_
   - `MOD09GA.006 Terra Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GA>`_
   - `MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09Q1>`_
   - `MOD09A1.006 Terra Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1>`_
   - `MYD09GQ.006 Aqua Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GQ>`_
   - `MYD09GA.006 Aqua Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GA>`_
   - `MYD09Q1.006 Aqua Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09Q1>`_
   - `MYD09A1.006 Aqua Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09A1>`_
- The :code:`maskClouds()` extended method for ee.Image and ee.ImageCollection classes now supports the following platforms:
   
   - `MYD09GA.006 Aqua Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GA>`_
   - `MYD09Q1.006 Aqua Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09Q1>`_
   - `MYD09A1.006 Aqua Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09A1>`_   
   - `MYD17A2H.006: Aqua Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD17A2H>`_   
   - `MYD13Q1.006 Aqua Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13Q1>`_
   - `MYD13A1.006 Aqua Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A1>`_
   - `MYD13A2.006 Aqua Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A2>`_
- The :code:`scale()` extended method for ee.Image and ee.ImageCollection classes now supports the following platforms:
   
   - `MYD09GQ.006 Aqua Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GQ>`_
   - `MYD09GA.006 Aqua Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GA>`_
   - `MYD09Q1.006 Aqua Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09Q1>`_
   - `MYD09A1.006 Aqua Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09A1>`_
   - `MYD10A1.006 Aqua Snow Cover Daily Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD10A1>`_
   - `MYD11A1.006 Aqua Land Surface Temperature and Emissivity Daily Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD11A1>`_
   - `MYD11A2.006 Aqua Land Surface Temperature and Emissivity 8-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD11A2>`_
   - `MYDOCGA.006 Aqua Ocean Reflectance Daily Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYDOCGA>`_
   - `MYD14A1.006: Aqua Thermal Anomalies & Fire Daily Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD14A1>`_   
   - `MYD17A2H.006: Aqua Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD17A2H>`_
   - `MYD17A3HGF.006: Aqua Net Primary Production Gap-Filled Yearly Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD17A3HGF>`_   
   - `MYD13Q1.006 Aqua Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13Q1>`_
   - `MYD13A1.006 Aqua Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A1>`_
   - `MYD13A2.006 Aqua Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A2>`_
   - `MYD08_M3.061 Aqua Atmosphere Monthly Global Product <https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD08_M3>`_
- The following vegetation indices were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection:

   - 'EVI2' : Two-Band Enhanced Vegetation Index.
   
- The following burn indices were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection:

   - 'CSIT' : Char Soil Index Thermal.
   - 'NBRT' : Normalized Burn Ratio Thermal.
   - 'NDVIT' : Normalized Difference Vegetation Index Thermal
   - 'SAVIT' : Soil-Adjusted Vegetation Index Thermal.

v0.1.7
--------------

New Modules
~~~~~~~~~~~~~~~~~~~~~~

- The *pd.DataFrame* module was created.
- The *common* module was created (it feeds the :code:`index()`, :code:`scale()` and :code:`maskClouds()` methods for both ee.Image and ee.ImageCollection).

New Features
~~~~~~~~~~~~~~~~~~~~~~

- The :code:`toEEFeatureCollection()` extended method for pd.DataFrame classes was created.
- The binary operators (+, -, \*\, /, //, %, \**\ , <<, >>, &, |) were overloaded for ee.Image objects.
- The rich comparisons (<, <=, ==, !=, >, >=) were overloaded for ee.Image objects.
- The unary operators (-, ~) were overloaded for ee.Image objects.

Improvements
~~~~~~~~~~~~~~~~~~~~~~

- *Exceptions* and *Warnings* were added to most methods.
- Conflicts between the Gain factor and the Green band in the :code:`index()` method were solved.
- :code:`tolerance` and :code:`unit` parameters were added to the :code:`closest()` extended method for ee.ImageCollection classes.
- The :code:`maskClouds()` extended method for ee.Image and ee.ImageCollection classes now supports the following platforms:

   - `Sentinel-3 OLCI EFR: Ocean and Land Color Instrument Earth Observation Full Resolution <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S3_OLCI>`_
   - `MOD09GA.006 Terra Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GA>`_
   - `MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09Q1>`_
   - `MOD09A1.006 Terra Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1>`_
   - `MCD15A3H.006 MODIS Leaf Area Index/FPAR 4-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD15A3H>`_
   - `MOD17A2H.006: Terra Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD17A2H>`_
   - `MOD16A2.006: Terra Net Evapotranspiration 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2>`_
   - `MOD13Q1.006 Terra Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13Q1>`_
   - `MOD13A1.006 Terra Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A1>`_
   - `MOD13A2.006 Terra Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A2>`_
- The :code:`scale()` extended method for ee.Image and ee.ImageCollection classes now supports the following platforms:

   - `Sentinel-3 OLCI EFR: Ocean and Land Color Instrument Earth Observation Full Resolution <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S3_OLCI>`_
   - `MCD43A4.006 MODIS Nadir BRDF-Adjusted Reflectance Daily 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD43A4>`_
   - `MCD43A3.006 MODIS Albedo Daily 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD43A3>`_
   - `MOD09GQ.006 Terra Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GQ>`_
   - `MOD09GA.006 Terra Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GA>`_
   - `MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09Q1>`_
   - `MOD09A1.006 Terra Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1>`_
   - `MOD10A1.006 Terra Snow Cover Daily Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD10A1>`_
   - `MOD11A1.006 Terra Land Surface Temperature and Emissivity Daily Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD11A1>`_
   - `MOD11A2.006 Terra Land Surface Temperature and Emissivity 8-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD11A2>`_
   - `MODOCGA.006 Terra Ocean Reflectance Daily Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MODOCGA>`_
   - `MOD14A1.006: Terra Thermal Anomalies & Fire Daily Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD14A1>`_
   - `MCD43A1.006 MODIS BRDF-Albedo Model Parameters Daily 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD43A1>`_
   - `MCD15A3H.006 MODIS Leaf Area Index/FPAR 4-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD15A3H>`_
   - `MOD17A2H.006: Terra Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD17A2H>`_
   - `MOD17A3HGF.006: Terra Net Primary Production Gap-Filled Yearly Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD17A3HGF>`_
   - `MOD16A2.006: Terra Net Evapotranspiration 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2>`_
   - `MOD13Q1.006 Terra Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13Q1>`_
   - `MOD13A1.006 Terra Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A1>`_
   - `MOD13A2.006 Terra Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A2>`_
   - `MOD08_M3.061 Terra Atmosphere Monthly Global Product <https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD08_M3>`_
- The following vegetation indices were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection:

   - 'GBNDVI' : Green-Blue Normalized Difference Vegetation Index.
   - 'GRNDVI' : Green-Red Normalized Difference Vegetation Index.
   - 'MNDVI' : Modified Normalized Difference Vegetation Index.
- The following snow indices were added to the :code:`index()` extended method for ee.Image and ee.ImageCollection:

   - 'NDSI' : Normalized Difference Snow Index.
- The 'SR' vegetation index was replaced by 'RVI' in the :code:`index()` extended method for ee.Image and ee.ImageCollection.Masking Clouds and Shadows
====================================

Masking clouds and shadows may seem hard, but it isn't! Let's take a look on it!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the following methods:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

   maskClouds
   preprocess
      
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   maskClouds
   preprocess

Supported Platforms
----------------------

This method automatically masks clouds and shadows on the following supported satellite platforms:

Sentinel Missions
~~~~~~~~~~~~~~~~~~~

- `Sentinel-3 OLCI EFR: Ocean and Land Color Instrument Earth Observation Full Resolution <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S3_OLCI>`_
- `Sentinel-2 MSI: MultiSpectral Instrument, Level-2A <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR?hl=en>`_

Landsat Missions
~~~~~~~~~~~~~~~~~~~

- `USGS Landsat 8 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR>`_
- `USGS Landsat 8 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2>`_
- `USGS Landsat 7 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_SR>`_
- `USGS Landsat 7 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2>`_
- `USGS Landsat 5 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1_SR>`_
- `USGS Landsat 4 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C01_T1_SR>`_

MODIS Products (Terra + Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MCD15A3H.006 MODIS Leaf Area Index/FPAR 4-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD15A3H>`_

MODIS Products (Terra)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MOD09GA.006 Terra Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GA>`_
- `MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09Q1>`_
- `MOD09A1.006 Terra Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1>`_
- `MOD17A2H.006: Terra Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD17A2H>`_
- `MOD16A2.006: Terra Net Evapotranspiration 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2>`_
- `MOD13Q1.006 Terra Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13Q1>`_
- `MOD13A1.006 Terra Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A1>`_
- `MOD13A2.006 Terra Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A2>`_

MODIS Products (Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MYD09GA.006 Aqua Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GA>`_
- `MYD09Q1.006 Aqua Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09Q1>`_
- `MYD09A1.006 Aqua Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09A1>`_   
- `MYD17A2H.006: Aqua Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD17A2H>`_   
- `MYD13Q1.006 Aqua Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13Q1>`_
- `MYD13A1.006 Aqua Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A1>`_
- `MYD13A2.006 Aqua Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A2>`_

VIIRS Products
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `VNP09GA: VIIRS Surface Reflectance Daily 500m and 1km <https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_001_VNP09GA?hl=en>`_
- `VNP13A1: VIIRS Vegetation Indices 16-Day 500m <https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_001_VNP13A1?hl=en>`_

.. warning::
   Not supported satellite platforms will raise a *Warning*.   

QA Method
----------------------

By default, the :code:`maskClouds()` uses the QA band of each paltform to compute the clouds and shadows masks (except for Sentinel-2, where the default method is Cloud Probability). The following table shows the band and the bits used for each platform (The value in parentheses is the valid value of the bitmask):

.. list-table:: QA bits used for clouds/shadows masking
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Platform
     - QA Band
     - Cloud Bitmask
     - Cirrus Bitmask
     - Shadow Bitmask
   * - Sentinel-3
     - quality_flags
     - 27 (0)
     -
     -
   * - Sentinel-2
     - QA60
     - 10 (0)
     - 11 (0)
     -
   * - Landsat Series
     - pixel_qa
     - 5 (0)
     - 
     - 3 (0)
   * - Landsat Collection 2
     - QA_PIXEL
     - 3 (0)
     - 2 (0)
     - 4 (0)
   * - MOD09GA
     - state_1km
     - 0 (0)
     - 8 (0)
     - 2 (0)
   * - MOD09Q1
     - State
     - 0 (0)
     - 8 (0)
     - 2 (0)
   * - MOD09A1
     - StateQA
     - 0 (0)
     - 8 (0)
     - 2 (0)
   * - MCD15A3H
     - FparExtra_QC
     - 5 (0)
     - 4 (0)
     - 6 (0)
   * - MOD17A2H
     - Psn_QC
     - 3 (0)
     - 
     - 
   * - MOD16A2
     - ET_QC
     - 3 (0)
     - 
     - 
   * - MOD13Q1
     - SummaryQA
     - 0 (0)
     - 
     - 
   * - MOD13A1
     - SummaryQA
     - 0 (0)
     - 
     -
   * - VNP09GA
     - QF1, QF2
     - 2 (0)
     - 6 (0), 7 (0)
     - 3 (0)

Usage
-----

Let's check how to use the :code:`maskClouds()` method for different platforms:

Sentinel-3
~~~~~~~~~~~~~~

On Sentinel 3, clouds are masked according to the bright pixels in the quality_flags band of the `Sentinel-3 OLCI EFR: Ocean and Land Color Instrument Earth Observation Full Resolution <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S3_OLCI>`_.

.. warning::
   This method may mask water as well on Sentinel-3 images. 

Let's take the Sentinel-3 image collection as example:

.. code-block:: python

   S3 = ee.ImageCollection('COPERNICUS/S3/OLCI')
   
There is no need to specify any arguments, since they're ignored.

.. code-block:: python

   S3.maskClouds()
  
This method can also be applied to a single image:

.. code-block:: python

   S3.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   S3.scale().maskClouds()

Sentinel-2
~~~~~~~~~~~~~~~~

On Sentinel 2, clouds can be masked using two methods: *QA* and *Cloud Probability*. The *QA* method uses the QA60 band in the `Surface Reflectance Product <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR?hl=en>`_ to mask clouds, while the *Cloud Probability* method uses the
`COPERNICUS/S2_CLOUD_PROBABILITY <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_CLOUD_PROBABILITY?hl=en>`_ collection to do it.

Shadows are masked based on the clouds mask, where shadows are searched within a range from clouds edges in the shadows direction.

.. seealso::
   For more info on masking shadows, please visit
   `‘Braaten, J. 2020. Sentinel-2 Cloud Masking with s2cloudless. Google Earth Engine, Community Tutorials’ <https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless>`_.

First, let's take the Sentinel-2 image collection:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR')

In order to use the *QA* method, it must be specified using the :code:`method` parameter:

.. code-block:: python

   S2.maskClouds(method = 'qa')
   
This line maps the *QA* masking method over the whole collection, but the method can also be applied to a single image:

.. code-block:: python

   S2.first().maskClouds(method = 'qa')
   
The *QA* method gives us the option to avoid masking cirrus clouds, but it must be specified using the :code:`maskCirrus` parameter:

.. code-block:: python

   S2.maskClouds(method = 'qa', maskCirrus = False)

And we can also avoid masking shadows by specifying the :code:`maskShadows` parameter:

.. code-block:: python

   S2.maskClouds(method = 'qa', maskShadows = False)
   
Now, in order to use the *Cloud Probability* method, we can specify it in the :code:`method` parameter:

.. code-block:: python

   S2.maskClouds(method = 'cloud_prob')
   
But, it is the default method, so you can just let the extended method with no additional parameters:

.. code-block:: python

   S2.maskClouds()
   
The *Cloud Probability* method uses a probability threshold to mask clouds, by default, the threshold is set to 60, but it can be modified using the :code:`prob` parameter:

.. code-block:: python

   S2.maskClouds(prob = 70)
   
If your image or collection is scaled, the :code:`scaledImage` parameter must be set to :code:`True`:

.. code-block:: python

   S2.scale().maskClouds(scaledImage = True)
   
In order to search for shadows, portental shadow pixels must be specified. Pixels with a NIR reflectance below 0.15 are considered potential shadow pixels, but this can be modified using the
:code:`dark` parameter:

.. code-block:: python

   S2.maskClouds(dark = 0.2)
   
Shadows are searched whitin a maximum range of 1000 m in the shadow direction from cloud edges, but this range can be modified using the :code:`cloudDist` parameter:

.. code-block:: python

   S2.maskClouds(cloudDist = 1500)
   
After finding all clouds and shadows, the mask can be dilated to avoid border effects. By default, clouds and shadows are dilated by 250 m, but this can be modified using the :code:`buffer` parameter:

.. code-block:: python

   S2.maskClouds(buffer = 100)
   
Finally, in order to avoid confusion between clouds and bright surface objects, the Cloud Displacement Index (CDI) can be used. By default, the CDI is not used, but it can be modified
using the :code:`cdi` parameter:

.. code-block:: python

   S2.maskClouds(cdi = -0.5)
   
.. seealso::
   For more info on CDI, please visit
   `‘Frantz, D., HaS, E., Uhl, A., Stoffels, J., Hill, J. 2018. Improvement of the Fmask algorithm for Sentinel-2 images: 
   Separating clouds from bright surfaces based on parallax effects. Remote Sensing of Environment 2015: 471-481’ 
   <https://www.sciencedirect.com/science/article/pii/S0034425718302037#:~:text=In%20this%20paper%2C%20we%20present,separated%20from%20bright%20ground%20objects.>`_.

Landsat Series
~~~~~~~~~~~~~~~~

On Landsat Series, both clouds and shadows are masked based on the pixel_qa band in the `Surface Reflectance Products <https://developers.google.com/earth-engine/datasets/catalog/landsat>`_.

Let's take the Landsat 8 image collection as example:

.. code-block:: python

   L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
   
There is no need to specify most of the arguments showed for Sentinel-2, since they're ignored.

.. code-block:: python

   L8.maskClouds()
   
Shadows are masked by default, but if required, the :code:`maskShadows` parameter can be modified.

.. code-block:: python

   L8.maskClouds(maskShadows = False)
   
This method can also be applied to a single image:

.. code-block:: python

   L8.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   L8.scale().maskClouds()
   
MODIS Products
~~~~~~~~~~~~~~~~

On MODIS Products, clouds and shadows are masked according to the specific QA band.

Let's take the MOD13Q1 image collection as example:

.. code-block:: python

   MOD13Q1 = ee.ImageCollection('MODIS/006/MOD13Q1')
   
There is no need to specify most of the arguments showed for Sentinel-2, since they're ignored.

.. code-block:: python

   MOD13Q1.maskClouds()
   
MOD13Q1, MOD13A1, MOD17A2H and MOD16A2 products don't have cirrus and shadow bitmasks, therefore, the arguments :code:`maskShadows` and :code:`maskCirrus` are ignored. MOD09GA, MOD09Q1, MOD09A1 and MCD15A3H products have cirrus and shadows bitmasks, and by default, they are set to *True*. If required, they can be set to *False*:

.. code-block:: python

   MOD09GA = ee.ImageCollection('MODIS/006/MOD09GA').maskClouds(maskShadows = False, maskCirrus = False)
   
This method can also be applied to a single image:

.. code-block:: python

   MOD09GA.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   MOD09GA.scale().maskClouds()
   
MOD13A2 doesn't have a bitmask QA band, instead, it has a Class QA band, where a value of zero means that the pixel has good data.

.. code-block:: python

   MOD13A2 = ee.ImageCollection('MODIS/006/MOD13A2').maskClouds()
   
VIIRS Products
~~~~~~~~~~~~~~~~

On VIIRS Products, clouds and shadows are masked according to the specific QA band.

Let's take the VNP09GA image collection as example:

.. code-block:: python

   VNP09GA = ee.ImageCollection('NOAA/VIIRS/001/VNP13A1')
   
There is no need to specify most of the arguments showed for Sentinel-2, since they're ignored.

.. code-block:: python

   VNP09GA.maskClouds()
   
If required, the arguments :code:`maskShadows` and :code:`maskCirrus` can be set to *False*:

.. code-block:: python

   VNP09GA.maskClouds(maskShadows = False, maskCirrus = False)
   
This method can also be applied to a single image:

.. code-block:: python

   VNP09GA.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   VNP09GA.scale().maskClouds()
   
VNP13A1 doesn't have a bitmask QA band, instead, it has a Class QA band, where a value of 9 means that the pixel is a cloud, while a value of 7 means that the pixel is a cloud shadows.

.. code-block:: python

   VNP13A1 = ee.ImageCollection('NOAA/VIIRS/001/VNP13A1').maskClouds()
   
If required, the argument :code:`maskShadows` can be set to *False*:

.. code-block:: python

   VNP13A1 = ee.ImageCollection('NOAA/VIIRS/001/VNP13A1').maskClouds(maskShadows = False)
   
Preprocessing
~~~~~~~~~~~~~~~
   
Additionally, if required, the :code:`preprocess()` method can be used to scale and mask images and image collections:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').preprocess()
   
The arguments for the :code:`maskClouds()` method can be used inside the :code:`preprocess()` method:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').preprocess(prob = 55,maskCirrus = False,buffer = 300,cdi = -0.5)Image Scaling
====================================

Image scaling now is A LOT EASIER with eemont! Let's see how!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the following methods:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

   getOffsetParams
   getScaleParams
   preprocess
   scaleAndOffset
      
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   getOffsetParams
   getScaleParams
   preprocess
   scaleAndOffset

Supported Platforms
----------------------

This method automatically scales ALL images from the `Google Earth Engine STAC Catalog <https://developers.google.com/earth-engine/datasets>`_ by using the `List of Scale and Offset Paramaters from the GEE STAC Catalog Repository <https://github.com/davemlz/ee-catalog-scale-offset-params>`_.

Usage
------------------

The :code:`scaleAndOffset()` method scales each image according to the *Scale* and *Offset* parameters for each band of the image. 

Let's take the Sentinel-2 SR image collection as example:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR')
   
The spectral bands from Sentinel-2 have values close to the (0, 10000) range, but they're unscaled. In order to get the real values, each spectral band must be multiplied by 0.0001, while AOT and WVP bands must be multiplied by 0.001. This scaling is automatically done by the :code:`scaleAndOffset()` method, without any additional parameters:

.. code-block:: python

   S2.scaleAndOffset()

The *Scale* and *Offset* parameters vary according to the satellite platform and the :code:`scaleAndOffset()` method detects the platform and do the scaling according to its parameters.

Let's take now the MOD11A2 product from MODIS. The LST_Day_1km and LST_Night_1km bands must be multiplied by 0.02, the Day_view_time and Night_view_time bands must be multiplied by 0.1,
the Emis_31 and Emis_32 bands must be multiplied by 0.002 and added by 0.49, while the Day_view_angl and Night_view_angl bands must be added by -65. All of this scaling
is simply done by the :code:`scaleAndOffset()` method:

.. code-block:: python

   MOD11A2scaled = ee.ImageCollection('MODIS/006/MOD11A2').scaleAndOffset()
   
The :code:`scale()` method can be applied to single images as well:

.. code-block:: python

   MOD11A2scaled = ee.ImageCollection('MODIS/006/MOD11A2').first().scaleAndOffset()
   
Preprocessing
~~~~~~~~~~~~~~~~~
   
Additionally, if required, the :code:`preprocess()` method can be used to scale and mask images and image collections:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').preprocess()Constructors
====================================

Let's see how to use the extended constructors available through eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Geometry, ee.Feature and ee.FeatureCollection classes with the following constructors:

ee.Geometry
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.geometry

.. autosummary::

   BBoxFromQuery
   LinearRingFromPlusCodes
   LineStringFromPlusCodes
   MultiLineStringFromPlusCodes
   MultiPointFromPlusCodes
   MultiPointFromQuery
   MultiPolygonFromPlusCodes
   PointFromPlusCode
   PointFromQuery
   PolygonFromPlusCodes
   RectangleFromPlusCodes   

ee.Feature
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.feature

.. autosummary::

   BBoxFromQuery
   PointFromQuery

ee.FeatureCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.featurecollection

.. autosummary::

   MultiPointFromQuery
   
Usage
------------------

Constructors By Query
~~~~~~~~~~~~~~~~~~~~~~~~

The constructors by query use the `geopy package <https://pypi.org/project/geopy/>`_ to construct geometries, features and feature collections.

The simplest geometry to construct is the Point, and can be constructed for ee.Geometry and ee.Feature classes;

.. code-block:: python

   geometry = ee.Geometry.PointFromQuery('Cali, Colombia',user_agent = 'eemont-user-guide-constructors')
   feature = ee.Feature.PointFromQuery('Cali, Colombia',user_agent = 'eemont-user-guide-constructors')
   
It has to be noted that the :code:`user_agent` argument is mandatory and it can be set to any string representing the app or the GEE username.

.. warning::
   If the :code:`user_agent` argument is not defined, the constructor will raise an Exception.

By default, to geocode the query, the Nominatim geocoder is used. If required, this parameter can be modified:

.. code-block:: python

   geometry = ee.Geometry.PointFromQuery('Cali, Colombia',geocoder = 'arcgis',user_agent = 'eemont-user-guide-constructors')
                                 
The second geometry to construct is the MultiPoint, and can be constructed for ee.Geometry and ee.FeatureCollection classes:

.. code-block:: python

   geometry = ee.Geometry.MultiPointFromQuery('Colombia',user_agent = 'eemont-user-guide-constructors')
   feature_collection = ee.FeatureCollection.MultiPointFromQuery('Colombia',user_agent = 'eemont-user-guide-constructors')

.. note::
   When a query is submitted, a set of locations is retrieved. The MultiPoint constructors create a class taking all locations into account. The Point constructors just take the first one.

The last geometry to construct is the Bounding Box, and can be constructed for ee.Geometry and ee.Feature classes:

.. code-block:: python

   geometry = ee.Geometry.BBoxFromQuery('Europe',user_agent = 'eemont-user-guide-constructors')
   feature = ee.Feature.BBoxFromQuery('Europe',user_agent = 'eemont-user-guide-constructors')
   
.. note::
   When using constructors for ee.Feature and ee.FeatureCollection classes, the raw properties of the location, or locations, are set for the feature or feature collection.
   
Constructors By Plus Codes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
   The awesome Plus Codes constructors and methods were created by `Aaron Zuspan <https://github.com/aazuspan>`_ (creator of `sankee <https://github.com/aazuspan/sankee>`_).

`Plus Codes <https://maps.google.com/pluscodes/>`_ are street addresses that represent an area based on longitude and latitude coordinates (e.g. 
"89MH+PW").

.. warning::
   In order to use Plus Codes constructors, it is required to install :code:`openlocationcode`. Please install it by running :code:`pip install openlocationcode`.
   
There are two ways to use the Plus Codes constructors.

1. By using a full Plus Code (e.g. "85FPQXGV+XH").
2. By using a short Plus Code (e.g. "QXGV+XH Denver, CO, USA").

When using full Plus Codes, just use it as the principal argument in the constructor:

.. code-block:: python

   point = ee.Geometry.PointFromPlusCode("85FPQXGV+XH")
   
When using a short Plus Code, it is required to use a geocoder (just like the constructors by queries), and therefore, an user agent must be declared:

.. code-block:: python

   point = ee.Geometry.PointFromPlusCode("QXGV+XH Denver, CO, USA",user_agent = 'eemont-user-guide-constructors')
   
More complex geometries can be constructed using a list of Plus Codes or a nested list of Plus Codes:

.. code-block:: python

   codes = ['85FQ2222+22', '85FR2222+22', '85GR2222+22']
   
   multipoint = ee.Geometry.MultiPointFromPlusCodes(codes)
   linestring = ee.Geometry.LineStringFromPlusCodes(codes)
   polygon = ee.Geometry.PolygonFromPlusCodes(codes)
   
   nestedCodes = [
        ['85FQ2222+22', '85FR2222+22', '85GR2222+22'], 
        ['85FP8PC2+G2', '85FPJF23+G4', '85FPMW2R+RP'],
        ['85FPRJ5W+82', '85GP3M67+CP', '85GQ2R7C+38'],
    ]
    
    multilinestring = ee.Geometry.MultiLineStringFromPlusCodes(nestedCodes)
    multipolygon = ee.Geometry.MultiPolygonFromPlusCodes(nestedCodes)Tasseled Cap
====================================
Guide by `Aaron Zuspan <https://github.com/aazuspan>`_


Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the method :code:`tasseledCap()`:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

    tasseledCap

ee.ImageCollection
~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

    tasseledCap


Supported Platforms
----------------------

Tasseled cap transformation coefficients are published for a number of sensors and processing levels:

* `Sentinel-2 MSI Level 1C <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2>`_ [1]_
* `Landsat 8 OLI TOA <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_TOA>`_ [2]_
* `Landsat 7 ETM+ TOA <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_TOA>`_ [3]_
* `Landsat 5 TM Raw DN <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1>`_ [4]_
* `Landsat 4 TM Raw DN <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C01_T1>`_ [5]_
* `Landsat 4 TM Surface Reflectance <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2>`_ [6]_
* `MODIS NBAR <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD43A4>`_ [7]_


Usage
------------------
Tasseled cap transformation calculates brightness, greenness, and wetness components from spectral bands. These data transformations can be helpful for 
both visualization and dimensionality reduction. 

To calculate tasseled cap components, we'll load a supported image (Sentinel-2) and run the :code:`tasseledCap` method.

.. code-block:: python

    img = ee.Image("COPERNICUS/S2/20160111T112432_20160111T113311_T28PDT")
    img = img.tasseledCap()

Now if we check the image band names, we'll see tasseled cap brightness (TCB), greenness (TCG), and wetness (TCW) component bands 
have been added to our image.

.. code-block:: python

    >>> img.bandNames().getInfo()
    ['B1', 'B2', 'B3' ... 'TCB', 'TCG', 'TCW']

.. warning::

    :code:`tasseledCap` requires a specific set of bands depending on the dataset. If your image does not contain the correct bands, 
    you will receive an error like :code:`Image.select: Pattern 'B1' did not match any bands`. To avoid that, it is safest to run :code:`tasseledCap`
    before subsetting your image bands.

For convenience, we could also run the :code:`tasseledCap` method on an :code:`ee.ImageCollection`, which would add TCB, TCG, and TCW bands
to each image.

.. code-block:: python

    col = ee.ImageCollection("COPERNICUS/S2")
    col = col.tasseledCap()


References
----------

.. [1] Shi, T., & Xu, H. (2019). Derivation of Tasseled Cap Transformation
    Coefficients for Sentinel-2 MSI At-Sensor Reflectance Data. IEEE Journal
    of Selected Topics in Applied Earth Observations and Remote Sensing, 1–11.
    doi:10.1109/jstars.2019.2938388
.. [2] Baig, M.H.A., Zhang, L., Shuai, T. and Tong, Q., 2014. Derivation of a
    tasselled cap transformation based on Landsat 8 at-satellite reflectance.
    Remote Sensing Letters, 5(5), pp.423-431.
.. [3] Huang, C., Wylie, B., Yang, L., Homer, C. and Zylstra, G., 2002.
    Derivation of a tasselled cap transformation based on Landsat 7 at-satellite
    reflectance. International journal of remote sensing, 23(8), pp.1741-1748.
.. [4] Crist, E.P., Laurin, R. and Cicone, R.C., 1986, September. Vegetation and
    soils information contained in transformed Thematic Mapper data. In
    Proceedings of IGARSS’86 symposium (pp. 1465-1470). Paris: European Space
    Agency Publications Division.
.. [5] Crist, E.P. and Cicone, R.C., 1984. A physically-based transformation of
    Thematic Mapper data---The TM Tasseled Cap. IEEE Transactions on Geoscience
    and Remote sensing, (3), pp.256-263.
.. [6] Crist, E.P., 1985. A TM tasseled cap equivalent transformation for
    reflectance factor data. Remote sensing of Environment, 17(3), pp.301-306.
.. [7] Lobser, S.E. and Cohen, W.B., 2007. MODIS tasselled cap: land cover
    characteristics expressed through transformed MODIS data. International
    Journal of Remote Sensing, 28(22), pp.5079-5101.
Closest Image to a Specific Date
====================================

Let's see how to get the closest image (or set of images) to a specific date.

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.ImageCollection class with the method :code:`closest()`:

.. currentmodule:: eemont.imagecollection

.. autosummary::

   closest

This method automatically filters any image collection to get the closest image to a specific date.

.. warning::
   This method uses the :code:`system:time_start` property, therefore, make sure your image collection has it!   

Usage
------------------

The :code:`closest()` method works on any image colection that has a :code:`system:time_start` property.

First, let's take the Sentinel-2 image collection as example:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR')
   
Now, we have to filter the collection to our ROI. The result of the :code:`closest()` method will vary depending on this.
Let's assume a single point is our ROI.

.. code-block:: python

   ROI = ee.Geometry.Point([-76.45, 4.32])
   S2 = S2.filterBounds(ROI)

Now, the :code:`closest()` method has just one parameter, :code:`date`, and this parameter can be a string...

.. code-block:: python

   S2.closest('2020-10-15')
   
Or an ee.Date class:

.. code-block:: python

   dateOfInterest = ee.Date('2020-10-15')
   S2.closest(dateOfInterest)
   
Both chunks will give you the same result here: an ee.ImageCollection of size 1. The result has just one image since our ROI intersects just one scene.
To get that image as a single image, we can use the :code:`first()` method.

.. code-block:: python

   S2.closest('2020-10-15').first()

By default, the image collection is filtered according to +/- 1 month from the :code:`date` parameter (:code:`tolerance = 1` and :code:`unit = 'month'`). This is done to speed up the searching process, but if required (if there are not images in that range), the :code:`tolerance` and :code:`unit` parameters can be modified:

.. code-block:: python

   S2.closest('2020-10-15', tolerance = 2, unit = 'year')

Now, let's assume that our ROI is larger, in this case, a whole department (state) of Colombia:

.. code-block:: python

   ROI = (ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level1')
       .filter(ee.Filter.eq('ADM1_NAME','Valle Del Cauca')))
   S2 = ee.ImageCollection('COPERNICUS/S2_SR').filterBounds(ROI).closest('2020-10-15')

You'll note that the size of the resulting ee.ImageCollection here is greater than 1. This result has more than one image since our ROI now intersects more than one scene.
To get those images together as a single image, you can mosaic them or use an ee.Reducer, for example :code:`median()`.

.. code-block:: python
   
   S2.median()Overloaded Operators
====================================

Let's see how to use overloaded operators in eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.Number classes with the binary and unary operators (including rich comparisons).

ee.Image 
-------------------

Binary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of binary operators that are overloaded:

.. list-table:: Binary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Overloaded Operator
   * - Addition
     - Image1.add(Image2)
     - Image1 + Image2 
   * - Subtraction
     - Image1.subtract(Image2)
     - Image1 - Image2
   * - Multiplication
     - Image1.multiply(Image2)
     - Image1 * Image2
   * - Division
     - Image1.divide(Image2)
     - Image1 / Image2
   * - Floor Division
     - Image1.divide(Image2).floor()
     - Image1 // Image2
   * - Modulo
     - Image1.mod(Image2)
     - Image1 % Image2
   * - Power
     - Image1.pow(Image2)
     - Image1 ** Image2
   * - Left Shift
     - Image1.leftShift(Image2)
     - Image1 << Image2
   * - Right Shift
     - Image1.rightShift(Image2)
     - Image1 >> Image2
   * - And
     - Image1.And(Image2)
     - Image1 & Image2
   * - Or
     - Image1.Or(Image2)
     - Image1 | Image2
          
Rich Comparisons
~~~~~~~~~~~~~~~~~~~

The following table shows the list of rich comparisons that are overloaded:

.. list-table:: Rich comparisons.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Lower Than
     - Image1.lt(Image2)
     - Image1 < Image2 
   * - Lower Than or Equal
     - Image1.lte(Image2)
     - Image1 <= Image2
   * - Equal
     - Image1.eq(Image2)
     - Image1 == Image2
   * - Not Equal
     - Image1.neq(Image2)    
     - Image1 != Image2
   * - Greater Than
     - Image1.gt(Image2)
     - Image1 > Image2 
   * - Greater Than or Equal
     - Image1.gte(Image2)
     - Image1 >= Image2
   * - Equal
     - Image1.eq(Image2)
     - Image1 == Image2
     
Unary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of unary operators that are overloaded:

.. list-table:: Unary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Negation
     - Image.multiply(-1)
     - \-\ Image
   * - Invert
     - Image.Not()
     - ~ Image
     
ee.Number 
-------------------

Binary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of binary operators that are overloaded:

.. list-table:: Binary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Overloaded Operator
   * - Addition
     - Number1.add(Number2)
     - Number1 + Number2 
   * - Subtraction
     - Number1.subtract(Number2)
     - Number1 - Number2
   * - Multiplication
     - Number1.multiply(Number2)
     - Number1 * Number2
   * - Division
     - Number1.divide(Number2)
     - Image1 / Image2
   * - Floor Division
     - Number1.divide(Number2).floor()
     - Number1 // Number2
   * - Modulo
     - Number1.mod(Number2)
     - Number1 % Number2
   * - Power
     - Number1.pow(Number2)
     - Number1 ** Number2
   * - Left Shift
     - Number1.leftShift(Number2)
     - Number1 << Number2
   * - Right Shift
     - Number1.rightShift(Number2)
     - Number1 >> Number2
   * - And
     - Number1.And(Number2)
     - Number1 & Number2
   * - Or
     - Number1.Or(Number2)
     - Number1 | Number2
          
Rich Comparisons
~~~~~~~~~~~~~~~~~~~

The following table shows the list of rich comparisons that are overloaded:

.. list-table:: Rich comparisons.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Lower Than
     - Number1.lt(Number2)
     - Number1 < Number2 
   * - Lower Than or Equal
     - Number1.lte(Number2)
     - Number1 <= Number2
   * - Equal
     - Number1.eq(Number2)
     - Number1 == Number2
   * - Not Equal
     - Number1.neq(Number2)    
     - Number1 != Number2
   * - Greater Than
     - Number1.gt(Number2)
     - Number1 > Number2 
   * - Greater Than or Equal
     - Number1.gte(Number2)
     - Number1 >= Number2
   * - Equal
     - Number1.eq(Number2)
     - Number1 == Number2
     
Unary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of unary operators that are overloaded:

.. list-table:: Unary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Negation
     - Number.multiply(-1)
     - \-\ Number
   * - Invert
     - Number.Not()
     - ~ Number
     
ee.List 
-------------------

Binary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of binary operators that are overloaded:

.. list-table:: Binary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Overloaded Operator
   * - Concatenation
     - List1.cat(List2)
     - List1 + List2
   * - Repeat
     - ee.List.repeat(List,Value)
     - List * Value

Usage
------------------

Overloaded operators can be used on any ee.Image or ee.Number object. Let's see how to compute the EVI using overloaded operators!

Let's take the Sentinel-2 SR image collection as example (remember to scale your image or image collection!):

.. code-block:: python

   point = ee.Geometry.Point([-76.0269,2.92846])
   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
      .filterBounds(point)
      .sort('CLOUDY_PIXEL_PERCENTAGE')
      .first()
      .maskClouds()
      .scale())

Now, let's take apart the bands that we need (it is not necessary, but it's easier to use :code:`N` instead of :code:`S2.select('B8')`):

.. code-block:: python

   N = S2.select('B8')
   R = S2.select('B4')
   B = S2.select('B2')
   
And finally, let's compute the EVI using overloaded operators:

.. code-block:: python

   EVI = 2.5 * (N - R) / (N + 6.0 * R - 7.5 * B + 1.0)

Let's see another example, but using rich comparisons. We are going to compute a snow cover mask!

First, compute the NDSI:

.. code-block:: python

   S2 = S2.index('NDSI')   
   
And now, let's take apart the bands that we need:

.. code-block:: python

   NDSI = S2.select('NDSI')
   N = S2.select('B8')
   G = S2.select('B3')
   
Finally, compute the snow cover mask `(Hall et al., 2001) <https://modis.gsfc.nasa.gov/data/atbd/atbd_mod10.pdf>`_:

.. code-block:: python

   snowPixels = (NDSI > 0.4) & (N >= 0.1) & (G > 0.11)

And update the mask (if required):

.. code-block:: python

   S2 = S2.updateMask(snowPixels)Extensions
====================================

Did you know that you can use eemont inside QGIS or R? Let's see how!

QGIS
-----------

In order to use eemont inside QGIS, please follow these steps:

First, make sure that you have successfully installed the `Google Earth Engine Plugin for QGIS <https://gee-community.github.io/qgis-earthengine-plugin/>`_.

Then, open the OSGeo4W shell and run the following line:

.. code-block::

   py3_env
   
This will set the Python 3 environment. Afterwards, you can install eemont by running:

.. code-block::

   python -m pip install eemont
   
After installation, eemont can be used in the Python console inside QGIS:

.. code-block:: python

   import ee, eemont
   from ee_plugin import Map

   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
       .maskClouds()
       .scale()
       .index(['NDVI','EVI','GNDVI'])
       .first())

   Map.addLayer(S2,{'min':0,'max':1,'bands':'NDVI'},'NDVI',True)


R
------

In order to use eemont inside R, please follow these steps:

First, make sure that you have successfully installed the `rgee <https://github.com/r-spatial/rgee>`_ and `reticulate <https://rstudio.github.io/reticulate/>`_.

Then, open a new R script and run the following chunk:

.. code-block:: r

   library(rgee)
   library(reticulate)
   
   ee_Initialize()

Now, we are ready to go!

First, we have to install :code:`eemont` (if required):

.. code-block:: r

   py_install("eemont",pip = TRUE)
   
Then, :code:`eemont` can be imported!

.. code-block:: r

   eemont <- import("eemont")
   
All python methods are available here, let's take a look!

Define a point of interest:

.. code-block:: r

   point <- ee$Geometry$Point(c(-74.0592,11.3172))
   
Get and filter the Landsat 8 SR collection:

.. code-block:: r

   L8 <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$filterBounds(point)
   
And use :code:`eemont` as you wish!

.. code-block:: r

   L8 <- L8$maskClouds()$scale()$index("NDWI")Data Conversion
====================================

Let's see how to convert non-Earth Engine classes to Earth Engine classes.

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   import pandas as pd
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the pd.DataFrame classes with the method :code:`toEEFeatureCollection()`:

pd.DataFrame
~~~~~~~~

.. currentmodule:: eemont.dataframe

.. autosummary::

   toEEFeatureCollection

Methods
-----------

A table of availabe conversion options is shown below:

.. list-table:: Available options
   :widths: 30 30 40
   :header-rows: 1

   * - From
     - To
     - Method
   * - pd.DataFrame
     - ee.FeatureCollection
     - :code:`toEEFeatureCollection()`

Usage
------------------

Let's create a pandas data frame:

.. code-block:: python

   df = pd.DataFrame()
   df['lat'] = [2.92846, 4.8927]
   df['lon'] = [-76.0269, -75.3188]
   df['name'] = ['Nevado del Huila', 'Nevado del Ruiz']
   
This data frame can be easily converted into a ee.FeatureCollection (with no geometries) using the :code:`toEEFeatureCollection()` method for pd.DataFrame classes:

.. code-block:: python

   fcWithNoGeometries = df.toEEFeatureCollection()

If the data frame has latitude and longitude columns, these can be specified in the :code:`latitude` and :code:`longitude` parameters:

.. code-block:: python

   fcWithGeometries = df.toEEFeatureCollection(latitude = 'lat',longitude = 'lon')Time Series By Regions
====================================

Let's see how to extract time series by region (or regions) with eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont, geemap
   import pandas as pd
   import numpy as np
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.ImageCollection classes with the methods :code:`getTimeSeriesByRegion()` and :code:`getTimeSeriesByRegions()`:
    
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   getTimeSeriesByRegion
   getTimeSeriesByRegions
   
Usage
------------------

The :code:`getTimeSeriesByRegion()` and :code:`getTimeSeriesByRegions()` methods extract time series by region (or regions) according to the specified ee.Geometry, ee.Feature or ee.FeatureCollection.

Let's create an ee.FeatureCollection as example (two crops in France with two identifiers: A and B):

.. code-block:: python

   f1 = ee.Feature(ee.Geometry.Point([3.984770,48.767221]).buffer(50),{'ID':'A'})
   f2 = ee.Feature(ee.Geometry.Point([4.101367,48.748076]).buffer(50),{'ID':'B'})
   fc = ee.FeatureCollection([f1,f2])

Let's take the Sentinel-2 SR image collection as example (compute NDVI and EVI to extract their values):

.. code-block:: python

   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
      .filterBounds(fc)
      .filterDate('2020-01-01','2021-01-01')
      .maskClouds()
      .scale()
      .index(['EVI','NDVI']))

Time Series By Region
~~~~~~~~~~~~~~~~~~~~~~~~

Let's assume that we want the mean of all pixels inside our collection as a single geometry. In that case, we'll use the :code:`getTimeSeriesByRegion()`:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = ee.Reducer.mean(),
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10)
                                 
Here, we are extracting the EVI and NDVI time series from the S2 collection by the geometry of our feature collection. If required, we can use more than one reducer:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10)
                                 
Now we are not extracting just the mean values, but also the median values. A column named :code:`reducer` is created indicating the corresponding reducer.

We can also add arguments that are valid for the :code:`reduceRegion()` method:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 bestEffort = True,
                                 maxPixels = 1e13,
                                 tileScale = 2)

By default, the output date column is named :code:`reducer`, but it can be modified:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 dateColumn = 'my_date_column')
                                 
The date value is by defult retrieved in the Standard ISO format, but it can be changed by :code:`ms` (milliseconds) or any other format:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 dateFormat = 'YYYYMMdd')

When the region is masked and the reducer doesn't retrieve any value, a NA value of :code:`-9999` is used, but if required, it can be modified: 

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 naValue = -9999999)
                                 
Time Series By Regions
~~~~~~~~~~~~~~~~~~~~~~~~

Now, if we want a time series by each feature in our feature collection, we require the :code:`getTimeSeriesByRegions()` method:

.. code-block:: python

   ts = S2.getTimeSeriesByRegions(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                  collection = fc,
                                  bands = ['EVI','NDVI'],
                                  scale = 10)

The same parameters of the :code:`getTimeSeriesByRegion()` method can be used here, except for :code:`bestEffort` and :code:`maxPixels`:

.. code-block:: python

   ts = S2.getTimeSeriesByRegions(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                  collection = fc,
                                  bands = ['EVI','NDVI'],
                                  scale = 10,
                                  naValue = -99999999,
                                  dateColumn = 'my_date_colum',
                                  dateFormat = 'ms')
                                  
Conversion to Pandas
~~~~~~~~~~~~~~~~~~~~~~~~

The time series is always retrieved as an ee.FeatureCollection. To convert the collection to a pandas data frame, we'll use the geemap package:

.. code-block:: python

   tsPandas = geemap.ee_to_pandas(ts)
   
Then we can convert the NA value to a real NA and the date column to a datetime class:

.. code-block:: python

   tsPandas[tsPandas == -9999] = np.nan
   tsPandas['date'] = pd.to_datetime(tsPandas['date'],infer_datetime_format = True)Panchromatic Sharpening
====================================
Guide by `Aaron Zuspan <https://github.com/aazuspan>`_


Panchromatic sharpening is simple in eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the method :code:`panSharpen()`:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

   panSharpen

ee.ImageCollection
~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::
   
      panSharpen

Supported Platforms
----------------------

Pansharpening is supported for the following platforms:

Landsat Missions
~~~~~~~~~~~~~~~~~~~

- `USGS Landsat 8 Collection 1 Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1>`_
- `USGS Landsat 8 Collection 1 Tier 1 OLI <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LO08_C01_T1>`_
- `USGS Landsat 8 Collection 1 Tier 1 and Real-Time data <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_RT>`_
- `USGS Landsat 8 Collection 1 Tier 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T2>`_
- `USGS Landsat 8 Collection 1 Tier 2 OLI <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LO08_C01_T2>`_
- `USGS Landsat 8 Collection 1 Tier 1 TOA Reflectance <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_TOA>`_
- `USGS Landsat 8 Collection 1 Tier 1 and Real-Time data TOA Reflectance <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_RT_TOA>`_
- `USGS Landsat 8 Collection 1 Tier 2 TOA Reflectance <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T2_TOA>`_


- `USGS Landsat 7 Collection 1 Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1>`_
- `USGS Landsat 7 Collection 1 Tier 1 and Real-Time data <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_RT>`_
- `USGS Landsat 7 Collection 1 Tier 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T2>`_
- `USGS Landsat 7 Collection 1 Tier 1 TOA Reflectance <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_TOA>`_
- `USGS Landsat 7 Collection 1 Tier 1 and Real-Time data TOA Reflectance <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_RT_TOA>`_
- `USGS Landsat 7 Collection 1 Tier 2 TOA Reflectance <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T2_TOA>`_

.. warning::
    Surface Reflectance products do not contain panchromatic bands and do not support pan-sharpening.

Algorithms
-----------

The :code:`panSharpen` method can be run using a variety of different algorithms by setting the :code:`method` argument. The following sharpening algorithms are supported by eemont:

.. list-table:: Available algorithm methods.
   :widths: 25 65 15
   :header-rows: 1

   * - Method
     - Name
     - kwargs
   * - :code:`SFIM` (default)
     - Smoothing-Filter Based Intensity Modulation
     - No
   * - :code:`HPFA`
     - High-Pass Filter Addition
     - No
   * - :code:`PCS`
     - Principal Component Substitution
     - Yes
   * - :code:`SM`
     - Simple Mean
     - No


.. seealso::
     Some algorithms take additional keyword arguments (kwargs) such as :code:`maxPixels`, :code:`geometry`, or 
     :code:`bestEffort`. These are passed to :code:`ee.Image.reduceRegion`. More information on how to set these arguments
     can be found `here <https://developers.google.com/earth-engine/guides/reducers_reduce_region>`_.
 
Usage
------------------
Let's load a supported image from Landsat 8:

.. code-block:: python

    img = ee.Image("LANDSAT/LC08/C01/T1_TOA/LC08_047027_20160819")

And sharpen it using the :code:`panSharpen` method with the default :code:`SFIM` algorithm.

.. code-block:: python

    sharp = img.panSharpen()

Easy as that!

We can also try sharpening with a different algorithm. Remember that some algorirthms take additional keyword arguments.
Here, we'll specify :code:`maxPixels` to avoid pixel limitations with the :code:`PCS` algorithm.

.. code-block:: python

    sharp = img.panSharpen(method="PCS", maxPixels=1e13)

Pan-sharpening image collections is identical to sharpening images:

.. code-block:: python

    imgCollection = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA")
    sharpCollection = imgCollection.panSharpen()

Quality Assessment
~~~~~~~~~~~~~~~~
Quality assessment (QA) metrics can be calculated by the :code:`panSharpen` method to measure how much spectral distortion 
was introduced by sharpening and to compare different sharpening algorithms. eemont supports the following QA metrics:

.. list-table:: Available QA metrics.
   :widths: 20 65 25 25
   :header-rows: 1

   * - QA
     - Name
     - Ideal Value
     - Mode
   * - :code:`MSE`
     - Mean Squared Error
     - 0
     - Band
   * - :code:`RMSE`
     - Root-Mean Squared Error
     - 0
     - Band
   * - :code:`RASE`
     - Relative Average Spectral Error 
     - 0
     - Image
   * - :code:`ERGAS`
     - Dimensionless Global Relative Error of Synthesis
     - 0
     - Image
   * - :code:`DIV`
     - Difference in Variance
     - 0
     - Band
   * - :code:`bias`
     - Bias
     - 0
     - Band
   * - :code:`CC`
     - Correlation Coefficient
     - 1
     - Band
   * - :code:`CML`
     - Change in mean luminance
     - 1
     - Band
   * - :code:`CMC`
     - Change in mean contrast
     - 1
     - Band
   * - :code:`UIQI`
     - Universal Image Quality Index
     - 1
     - Band

.. note::
    Some metrics are calculated image-wise and others are calculated band-wise. Image-wise metrics return one value for
    the entire image while band-wise metrics return one value for each band.

QA metrics are calculated by passing a list of one or more metrics to the :code:`qa` argument of the :code:`panSharpen`
method. Below, we'll calculate :code:`RASE` and :code:`UIQI` while sharpening an image.

.. code-block:: python

    metrics = ["RASE", "UIQI"]
    img = ee.Image("LANDSAT/LC08/C01/T1_TOA/LC08_047027_20160819")
    sharp = img.panSharpen(qa=metrics, maxPixels=1e13)

.. seealso::
    All QA metrics take additional keyword arguments (kwargs) such as :code:`maxPixels`, :code:`geometry`, or 
    :code:`bestEffort`. These are passed to :code:`ee.Image.reduceRegion`. More information can be found `here <https://developers.google.com/earth-engine/guides/reducers_reduce_region>`_.

Calculated QA metrics are set as properties of the sharpened image and can be retrieved with the :code:`get` method.
QA property names use the format :code:`eemont:{QA}`.

.. code-block:: python

    sharp.get("eemont:RASE").getInfo()

Histogram Matching
====================================
Guide by `Aaron Zuspan <https://github.com/aazuspan>`_


Overview
-----------

The eemont package extends the ee.Image class with the method :code:`matchHistogram()`:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

    matchHistogram

Usage
------------------
Histogram matching performs band-wise adjustments to match the spectral response of one image to a target image. 
Let's look at a few examples of how histogram matching can be used for different platform combinations in eemont.

.. seealso::
    For more info on histogram matching and details on eemont's implementation, please visit `Braaten, J. 2021. Histogram 
    Matching. Google Earth Engine, Community Tutorials <https://developers.google.com/earth-engine/tutorials/community/histogram-matching>`_.

Landsat 5 and 8
~~~~~~~~~~~~~~~
In this example, we'll match the histograms of a Landsat 5 image to a Landsat 8 image. This can be helpful for performing 
consistent time series analysis using different sensors.

First, we'll load the source image and the target image.

.. code-block:: python

    source = ee.Image("LANDSAT/LT05/C01/T1_SR/LT05_195028_20110208")
    target = ee.Image("LANDSAT/LC08/C02/T1_L2/LC08_196027_20130714")

.. warning::
    Images must overlap. Cloud cover or major changes in land cover may cause inaccurate results.

Now, we'll specify which bands should be matched. We do this using a dictionary with the source bands as keys and target
bands as values.

.. code-block:: python

    bands = {
        # Red
        "B3": "SR_B4",
        # Green
        "B2": "SR_B3",
        # Blue
        "B1": "SR_B2"
    }

.. note::
    Any bands that aren't listed will be removed from the matched image.

Finally, we'll call the :code:`matchHistogram` method to create a new image representing the source image matched to the
target image.

.. code-block:: python

    matched = source.matchHistogram(target, bands)


Sentinel-2 and MODIS
~~~~~~~~~~~~~~~~~~~~
In this example, we'll match the histogram of a Sentinel-2 image to a MODIS image.

.. code-block:: python
    
    source = ee.Image("COPERNICUS/S2/20180923T081641_20180923T083023_T35PQQ")
    target = ee.Image("MODIS/006/MOD09A1/2018_08_05")

Specify the matching bands between images:

.. code-block:: python

    bands = {
        # Red
        "B4": "sur_refl_b01",
        # Green
        "B3": "sur_refl_b04",
        # Blue
        "B2": "sur_refl_b03"
    }

And match the Sentinel-2 image to the MODIS image!

.. code-block:: python

    matched = source.match(target, bands)


You can adjust the quality and speed of histogram matching using the :code:`maxBuckets` argument. Fewer buckets will run
faster but produce less accurate matching. By default, 256 buckets are used.

.. code-block:: python

    matched = source.matchHistogram(target, bands, maxBuckets=64)

.. note::

    :code:`maxBuckets` are automatically adjusted to the nearest power of 2, so :code:`maxBuckets=50` is the same as :code:`maxBuckets=64.`Spectral Indices Computation
====================================

Let's see how to compute built-in spectral indices with eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the method :code:`spectralIndices()`:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

   spectralIndices
      
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   spectralIndices

Supported Platforms
----------------------

This method automatically computes spectral indices for the following supported satellite platforms:

Sentinel Missions
~~~~~~~~~~~~~~~~~~~

- `Sentinel-2 MSI: MultiSpectral Instrument, Level-2A <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR?hl=en>`_
- `Sentinel-2 MSI: MultiSpectral Instrument, Level-1C <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2#bands>`_

Landsat Missions
~~~~~~~~~~~~~~~~~~~

- `USGS Landsat 8 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR>`_
- `USGS Landsat 8 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2>`_
- `USGS Landsat 7 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_SR>`_
- `USGS Landsat 7 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2>`_
- `USGS Landsat 5 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1_SR>`_
- `USGS Landsat 4 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C01_T1_SR>`_

MODIS Products (Terra + Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MCD43A4.006 MODIS Nadir BRDF-Adjusted Reflectance Daily 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD43A4>`_

MODIS Products (Terra)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MOD09GQ.006 Terra Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GQ>`_
- `MOD09GA.006 Terra Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GA>`_
- `MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09Q1>`_
- `MOD09A1.006 Terra Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1>`_

MODIS Products (Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MYD09GQ.006 Aqua Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GQ>`_
- `MYD09GA.006 Aqua Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GA>`_
- `MYD09Q1.006 Aqua Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09Q1>`_
- `MYD09A1.006 Aqua Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09A1>`_

.. important::
   It is highly recommended to scale the image (or image collection) before computing spectral indices. See the :code:`scaleAndOffset()` method for more info.  

List of Indices
----------------------

The list of indices is retrieved from the `Awesome Spectral Indices <https://github.com/davemlz/awesome-spectral-indices>`_ list.

Vegetation Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in vegetation indices:

.. list-table:: Built-in vegetation indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - BNDVI
     - Blue Normalized Difference Vegetation Index
     - `Index DataBase BNDVI <https://www.indexdatabase.de/db/i-single.php?id=135>`_  
   * - CIG
     - Chlorophyll Index - Green
     - `Index DataBase CIG <https://www.indexdatabase.de/db/i-single.php?id=128>`_
   * - CVI
     - Chlorophyll Vegetation Index
     - `Index DataBase CVI <https://www.indexdatabase.de/db/i-single.php?id=391>`_
   * - EVI
     - Enhanced Vegetation Index
     - `Index DataBase EVI <https://www.indexdatabase.de/db/i-single.php?id=16>`_
   * - EVI2
     - Two-Band Enhanced Vegetation Index
     - `(Jiang et al., 2008) <https://doi.org/10.1016/j.rse.2008.06.006>`_
   * - GARI
     - Green Atmospherically Resistant Vegetation Index
     - `Index DataBase GARI <https://www.indexdatabase.de/db/i-single.php?id=363>`_
   * - GBNDVI
     - Green-Blue Normalized Difference Vegetation Index
     - `Index DataBase GBNDVI <https://www.indexdatabase.de/db/i-single.php?id=186>`_
   * - GEMI
     - Global Environment Monitoring Index
     - `Index DataBase GEMI <https://www.indexdatabase.de/db/i-single.php?id=25>`_
   * - GLI
     - Green Leaf Index
     - `Index DataBase GLI <https://www.indexdatabase.de/db/i-single.php?id=375>`_
   * - GNDVI
     - Green Normalized Difference Vegetation Index
     - `Index DataBase GNDVI <https://www.indexdatabase.de/db/i-single.php?id=401>`_
   * - GRNDVI
     - Green-Red Normalized Difference Vegetation Index
     - `Index DataBase GRNDVI <https://www.indexdatabase.de/db/i-single.php?id=185>`_
   * - GVMI
     - Global Vegetation Moisture Index
     - `Index DataBase GVMI <https://www.indexdatabase.de/db/i-single.php?id=372>`_
   * - MNDVI
     - Modified Normalized Difference Vegetation Index
     - `Index DataBase MNDVI <https://www.indexdatabase.de/db/i-single.php?id=245>`_
   * - NDVI
     - Normalized Difference Vegetation Index
     - `Index DataBase NDVI <https://www.indexdatabase.de/db/i-single.php?id=58>`_
   * - NGRDI
     - Normalized Green Red Difference Index
     - `Index DataBase NGRDI <https://www.indexdatabase.de/db/i-single.php?id=390>`_
   * - RVI
     - Ratio Vegetation Index
     - `Index DataBase RVI <https://www.indexdatabase.de/db/i-single.php?id=72>`_
   * - SAVI
     - Soil-Adjusted Vegetation Index
     - `Index DataBase SAVI <https://www.indexdatabase.de/db/i-single.php?id=87>`_
   * - VARI
     - Visible Atmospherically Resistant Index
     - `Index DataBase VARI <https://www.indexdatabase.de/db/i-single.php?id=356>`_     
     
Burn Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in burn indices:

.. list-table:: Built-in burn indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - BAI
     - Burned Area Index
     - `(Martín, 1998) [spanish] <https://digital.csic.es/bitstream/10261/6426/1/Martin_Isabel_Serie_Geografica.pdf>`_ 
   * - BAIS2
     - Burned Area Index for Sentinel 2
     - `(Filipponi, 2018) <https://doi.org/10.3390/ecrs-2-05177>`_
   * - CSIT
     - Char Soil Index Thermal
     - `(Smith et al., 2007) <https://doi.org/10.1080/01431160600954704>`_
   * - NBR
     - Normalized Burn Ratio
     - `Index DataBase NBR <https://www.indexdatabase.de/db/i-single.php?id=53>`_
   * - NBRT
     - Normalized Burn Ratio Thermal
     - `(Holden et al., 2005) <https://doi.org/10.1080/01431160500239008>`_
   * - NDVIT
     - Normalized Difference Vegetation Index Thermal
     - `(Smith et al., 2007) <https://doi.org/10.1080/01431160600954704>`_
   * - SAVIT
     - Soil-Adjusted Vegetation Index Thermal
     - `(Smith et al., 2007) <https://doi.org/10.1080/01431160600954704>`_
     
Water Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in water indices:

.. list-table:: Built-in water indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - MNDWI
     - Modified Normalized Difference Water Index
     - `(Xu, 2006) <https://doi.org/10.1080/01431160600589179>`_  
   * - NDWI
     - Normalized Difference Water Index
     - `(McFeeters, 1996) <https://doi.org/10.1080/01431169608948714>`_
     
Snow Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in snow indices:

.. list-table:: Built-in snow indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - NDSI
     - Normalized Difference Snow Index
     - `(Riggs et al., 1994) <https://doi.org/10.1109/IGARSS.1994.399618>`_
     
Drought Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in drought indices:

.. list-table:: Built-in snow indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - NDDI
     - Normalized Difference Drought Index
     - `(Gu et al., 2007) <https://doi.org/10.1029/2006GL029127>`_
     
Generalized Kernel Indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in kernel indices:

.. list-table:: Built-in kernel indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - kEVI
     - Kernel Enhanced Vegetation Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_
   * - kNDVI
     - Kernel Normalized Difference Vegetation Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_
   * - kRVI
     - Kernel Ratio Vegetation Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_
   * - kVARI
     - Kernel Visible Atmospherically Resistant Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_

Kernels
----------------------

In the case of generalized kernel indices, the following kernels are available:

Linear Kernel
~~~~~~~~~~~~~~~~~~~

The linear kernel for generalized kernel indices can be selected by setting :code:`kernel = 'linear'`.

.. math::

   k(a,b) = ab
   
RBF Kernel
~~~~~~~~~~~~~~~~~~~

The Radial Basis Function (RBF) kernel for generalized kernel indices can be selected by setting :code:`kernel = 'RBF'`.

.. math::

   k(a,b) = exp(- \frac{(a - b) ^ 2}{2 \sigma ^ 2})
   
Where :math:`\sigma` is a free length-scale parameter.
   
Polynomial Kernel
~~~~~~~~~~~~~~~~~~~

The polynomial kernel for generalized kernel indices can be selected by setting :code:`kernel = 'poly'`.

.. math::

   k(a,b) = (ab + c) ^ p

Where :math:`c` is a free parameter that trades off the influence of higher-order versus lower-order terms and :math:`p` is the kernel degree.

List of Bands
----------------------

The following table shows the list of bands used for spectral indices computation:

.. list-table:: Bands used for spectral indices computation.
   :widths: 18 18 16 16 16 16
   :header-rows: 1

   * - Description
     - Name     
     - Sentinel-2
     - Landsat 8
     - Landsat 457
     - MODIS     
   * - Aerosols
     - A
     - B1
     - B1
     -
     -     
   * - Blue
     - B
     - B2
     - B2
     - B1
     - B3 
   * - Green
     - G
     - B3
     - B3
     - B2
     - B4    
   * - Red
     - R
     - B4
     - B4
     - B3
     - B1
   * - Red Edge 1
     - RE1
     - B5
     - 
     -
     -     
   * - Red Edge 2
     - RE2
     - B6
     - 
     -
     -     
   * - Red Edge 3
     - RE3
     - B7
     - 
     -
     -     
   * - Red Edge 4
     - RE4
     - B8A
     - 
     -
     -     
   * - NIR
     - N
     - B8
     - B5
     - B4
     - B2
   * - SWIR 1
     - S1
     - B11
     - B6
     - B5
     - B6     
   * - SWIR 2
     - S2
     - B12
     - B7
     - B7
     - B7   
   * - Thermal 1
     - T1
     - 
     - B10
     - B6
     -     
   * - Thermal 2
     - T2
     - 
     - B11
     - 
     -     

.. warning::
   If the satellite platform doesn't have the required bands for computing an index, it won't be computed.

Usage
------------------

The :code:`spectralIndices()` method computes the specified spectral index and adds it as a new band.

Let's take the Sentinel-2 SR image collection as example (remember to scale your image or image collection!):

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').scaleAndOffset()
   
By default, the :code:`spectralIndices()` method computes the NDVI:

.. code-block:: python

   S2withIndices = S2.spectralIndices()
   S2withIndices.select('NDVI')
   
If required, any of the above-mentioned indices can be computed by modifying the :code:`index` parameter:

.. code-block:: python

   S2withIndices = S2.spectralIndices(index = 'EVI')
   S2withIndices.select('EVI')
   
Specific index-parameters can be changed, for example, the canopy background adjustment L is set to 1.0 for EVI, but for SAVI it can be changed to 0.5:

.. code-block:: python

   S2withIndices = S2.spectralIndices('SAVI',L = 0.5)
   S2withIndices.select('SAVI')
   
If more than one index is required, a list of indices can be used:

.. code-block:: python

   S2withIndices = S2.spectralIndices(['CIG','NBR','NDWI'])
   S2withIndices.select('CIG')
   S2withIndices.select('NBR')
   S2withIndices.select('NDWI')
   
Indices can also be computed for single images:

.. code-block:: python

   S2withIndices = S2.first().spectralIndices(['GBNDVI','MNDVI','EVI'])
   S2withIndices.select('GBNDVI')
   S2withIndices.select('MNDVI')
   S2withIndices.select('EVI')
   
All vegetation indices can be computed by setting :code:`index = vegetation`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('vegetation')
   S2withIndices.select('NDVI')
   S2withIndices.select('GNDVI')
   S2withIndices.select('RVI')
   # ...
   
All burn indices can be computed by setting :code:`index = burn`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('burn')
   S2withIndices.select('BAI')
   S2withIndices.select('BAIS2')
   S2withIndices.select('NBR')
   
All water indices can be computed by setting :code:`index = water`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('water')
   S2withIndices.select('NDWI')
   S2withIndices.select('MNDWI')
   
All snow indices can be computed by setting :code:`index = snow`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('snow')
   S2withIndices.select('NDSI')
   
If you want to compute all available indices, you can set :code:`index = all`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('all')
   S2withIndices.select('NDVI')
   S2withIndices.select('BAI')
   S2withIndices.select('NDWI')
   S2withIndices.select('NDSI')
   # ...

Generalized Kernel Indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generalized kernel indices are availabe through eemont (e.g. kNDVI):

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI')
   S2withIndices.select('kNDVI')
   
By default, the RBF kernel is used and the :code:`sigma` parameter is :code:`0.5 * (a + b)` (this means, that for :code:`k(N,R)`, :code:`sigma = '0.5 * (N + R)'`). If required, :code:`sigma` can be modified by another expression (using :code:`a` and :code:`b`) or a float:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI',sigma = 1)
   S2withIndices.select('kNDVI')
   
The kernel can be modified by modifying the :code:`kernel` parameter:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI',kernel = 'poly')
   S2withIndices.select('kNDVI')
   
For the polynomial kernel, the :code:`p` and :code:`c` parameters can be modified:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI',kernel = 'poly',p = 4,c = 0)
   S2withIndices.select('kNDVI')
   
All kernel indices can be computed by setting :code:`index = kernel`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kernel')
   S2withIndices.select('kEVI')
   S2withIndices.select('kNDVI')
   S2withIndices.select('kRVI')
   S2withIndices.select('kVARI')
   
.. seealso::
   For more info on generalized kernel indices, please visit
   `‘Camps-Valls, G., et al. 2021. A unified vegetation index for quantifying the terrestrial biosphere. Science Advances 7 (9): eabc7447. Doi: 10.1126/sciadv.abc7447’ <https://doi.org/10.1126/sciadv.abc7447>`_.User Guide
===========

Here you can find the user guide for the extended methods:

.. toctree::
   :caption: Table of Contents
   :maxdepth: 3

   closestImage
   constructors
   dataConversion
   containers
   histogramMatching
   imageScaling
   maskingClouds
   overloadedOperators
   panSharpening
   spectralIndices   
   timeSeries
   tasseledCapEmulating Containers
====================================

Let's see how to use container emulation methods in eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image, ee.ImageCollection, ee.Feature, ee.FeatureCollection, ee.List and ee.Dictionary classes with container emulation methods.

ee.Image
~~~~~~~~~~~~~~~~~~~

The following table shows the list of container emulation methods that are overloaded:

.. list-table:: Container emulation methods.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Container Emulation Method
   * - Get Item By Key
     - Image.select(bands)
     - Image[bands]
   * - Get Item By Index
     - Image.select(index)
     - Image[index]
   * - Get Item By Slice
     - Image.slice(start,stop)
     - Image[start:stop]
     
ee.ImageCollection 
~~~~~~~~~~~~~~~~~~~

The following table shows the list of container emulation methods that are overloaded:

.. list-table:: Container emulation methods.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Container Emulation Method
   * - Get Item By Key
     - ImageCollection.select(bands)
     - ImageCollection[bands]
   * - Get Item By Index
     - ImageCollection.select(index)
     - ImageCollection[index]
   * - Get Item By Slice
     - ImageCollection.map(lambda x: x.slice(start,stop))
     - ImageCollection[start:stop]
   * - Length
     - ImageCollection.size().getInfo()
     - len(ImageCollection)

ee.Feature 
~~~~~~~~~~~~~~~~~~~

The following table shows the list of container emulation methods that are overloaded:

.. list-table:: Container emulation methods.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Container Emulation Method
   * - Get Item By Key
     - Feature.select(properties)
     - Feature[properties]

ee.FeatureCollection 
~~~~~~~~~~~~~~~~~~~~~~

The following table shows the list of container emulation methods that are overloaded:

.. list-table:: Container emulation methods.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Container Emulation Method
   * - Get Item By Key
     - FeatureCollection.select(properties)
     - FeatureCollection[properties]
   * - Length
     - FeatureCollection.size().getInfo()
     - len(FeatureCollection)
     
ee.List 
~~~~~~~~~~~~~~~~~~~

The following table shows the list of container emulation methods that are overloaded:

.. list-table:: Container emulation methods.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Container Emulation Method
   * - Get Item By Index
     - List.select(index)
     - List[index]
   * - Get Item By Slice
     - List.slice(start,stop)
     - List[start:stop]
   * - Length
     - List.length().getInfo()
     - len(List)
   * - Contains
     - List.contains(value).getInfo()
     - value in List
     
ee.Dictionary 
~~~~~~~~~~~~~~~~~~~

The following table shows the list of container emulation methods that are overloaded:

.. list-table:: Container emulation methods.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Container Emulation Method
   * - Get Item By Key
     - Dictionary.get(key)
     - Dictionary[key]
   * - Contains
     - Dictionary.contains(key).getInfo()
     - key in Dictionary

Usage
------------------

Container emulation methods can be used on any of the Earth Engine objects mentioned above. Let's see how to use them!

Raster Types
~~~~~~~~~~~~~~~~~~~

Let's take the Sentinel-2 SR image collection as example:

.. code-block:: python

   point = ee.Geometry.Point([-76.0269,2.92846])
   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
      .filterBounds(point)
      .sort('CLOUDY_PIXEL_PERCENTAGE')
      .first()
      .maskClouds()
      .scale())

Now, if we want to select a specific band, we can do it as follows:

.. code-block:: python

   NIR = S2['B8']
   
Or multiple bands:

.. code-block:: python

   NIRRED = S2[['B8','B4']]
   
We can also use regex!:

.. code-block:: python

   bands = S2['B.*']
   
Or an index:

.. code-block:: python

   BLUE = S2[1]
   
Or even better, a slice:

.. code-block:: python

   RGB = S2[1:4]
   
All of these methods can also be done for ee.ImageCollection objects:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').filterBounds(point)['B.*']
   
And, additionally, we can check the size of the image collection by using :code:`len()`:

.. code-block:: python

   len(S2)
 
Vector Types
~~~~~~~~~~~~~~~~~~~

Let's see another example, but using features:

.. code-block:: python

   WDPA = ee.FeatureCollection("WCMC/WDPA/current/polygons") 
   
And now, let's take some properties:

.. code-block:: python

   WDPA = WDPA[['WDPAID','NAME','REP_AREA']]
   
Now, let's check the size of the feature collection:

.. code-block:: python

   len(WDPA)
   
For the ee.List objects, we can also use container emulaion methods!

Lists
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   l = ee.List([100,120,230,310,450])

You can get an item by using its index:

.. code-block:: python

   l[0]
   
Or multiple items by using an slice:

.. code-block:: python

   l[1:4]
   
You can also check if an item is in the list:

.. code-block:: python

   370 in l
   
And get the length of the list:

.. code-block:: python

   len(l)

Dictionaries
~~~~~~~~~~~~~~~~~~~

Things work in a similar way for ee.Dictionary classes:

.. code-block:: python

   d = ee.Dictionary({'ID': 1,'Name': 'Natural Park','Area': 3240})
   
We can get a value by using its key:

.. code-block:: python

   d['Name']
   
We can also check if a key is in a dictionary:

.. code-block:: python

   'Area' in dExtensions
====================================

Did you know that you can use eemont inside QGIS or R? Let's see how!

QGIS
-----------

In order to use eemont inside QGIS, please follow these steps:

First, make sure that you have successfully installed the `Google Earth Engine Plugin for QGIS <https://gee-community.github.io/qgis-earthengine-plugin/>`_.

Then, open the OSGeo4W shell and run the following line:

.. code-block::

   py3_env
   
This will set the Python 3 environment. Afterwards, you can install eemont by running:

.. code-block::

   python -m pip install eemont
   
After installation, eemont can be used in the Python console inside QGIS:

.. code-block:: python

   import ee, eemont
   from ee_plugin import Map

   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
       .maskClouds()
       .scale()
       .index(['NDVI','EVI','GNDVI'])
       .first())

   Map.addLayer(S2,{'min':0,'max':1,'bands':'NDVI'},'NDVI',True)


R
------

In order to use eemont inside R, please follow these steps:

First, make sure that you have successfully installed the `rgee <https://github.com/r-spatial/rgee>`_ and `reticulate <https://rstudio.github.io/reticulate/>`_.

Then, open a new R script and run the following chunk:

.. code-block:: r

   library(rgee)
   library(reticulate)
   
   ee_Initialize()

Now, we are ready to go!

First, we have to install :code:`eemont` (if required):

.. code-block:: r

   py_install("eemont",pip = TRUE)
   
Then, :code:`eemont` can be imported!

.. code-block:: r

   eemont <- import("eemont")
   
All python methods are available here, let's take a look!

Define a point of interest:

.. code-block:: r

   point <- ee$Geometry$Point(c(-74.0592,11.3172))
   
Get and filter the Landsat 8 SR collection:

.. code-block:: r

   L8 <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$filterBounds(point)
   
And use :code:`eemont` as you wish!

.. code-block:: r

   L8 <- L8$maskClouds()$scale()$index("NDWI")Time Series By Regions
====================================

Let's see how to extract time series by region (or regions) with eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont, geemap
   import pandas as pd
   import numpy as np
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.ImageCollection classes with the methods :code:`getTimeSeriesByRegion()` and :code:`getTimeSeriesByRegions()`:
    
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   getTimeSeriesByRegion
   getTimeSeriesByRegions
   
Usage
------------------

The :code:`getTimeSeriesByRegion()` and :code:`getTimeSeriesByRegions()` methods extract time series by region (or regions) according to the specified ee.Geometry, ee.Feature or ee.FeatureCollection.

Let's create an ee.FeatureCollection as example (two crops in France with two identifiers: A and B):

.. code-block:: python

   f1 = ee.Feature(ee.Geometry.Point([3.984770,48.767221]).buffer(50),{'ID':'A'})
   f2 = ee.Feature(ee.Geometry.Point([4.101367,48.748076]).buffer(50),{'ID':'B'})
   fc = ee.FeatureCollection([f1,f2])

Let's take the Sentinel-2 SR image collection as example (compute NDVI and EVI to extract their values):

.. code-block:: python

   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
      .filterBounds(fc)
      .filterDate('2020-01-01','2021-01-01')
      .maskClouds()
      .scale()
      .index(['EVI','NDVI']))

Time Series By Region
~~~~~~~~~~~~~~~~~~~~~~~~

Let's assume that we want the mean of all pixels inside our collection as a single geometry. In that case, we'll use the :code:`getTimeSeriesByRegion()`:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = ee.Reducer.mean(),
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10)
                                 
Here, we are extracting the EVI and NDVI time series from the S2 collection by the geometry of our feature collection. If required, we can use more than one reducer:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10)
                                 
Now we are not extracting just the mean values, but also the median values. A column named :code:`reducer` is created indicating the corresponding reducer.

We can also add arguments that are valid for the :code:`reduceRegion()` method:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 bestEffort = True,
                                 maxPixels = 1e13,
                                 tileScale = 2)

By default, the output date column is named :code:`reducer`, but it can be modified:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 dateColumn = 'my_date_column')
                                 
The date value is by defult retrieved in the Standard ISO format, but it can be changed by :code:`ms` (milliseconds) or any other format:

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 dateFormat = 'YYYYMMdd')

When the region is masked and the reducer doesn't retrieve any value, a NA value of :code:`-9999` is used, but if required, it can be modified: 

.. code-block:: python

   ts = S2.getTimeSeriesByRegion(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                 geometry = fc,
                                 bands = ['EVI','NDVI'],
                                 scale = 10,
                                 naValue = -9999999)
                                 
Time Series By Regions
~~~~~~~~~~~~~~~~~~~~~~~~

Now, if we want a time series by each feature in our feature collection, we require the :code:`getTimeSeriesByRegions()` method:

.. code-block:: python

   ts = S2.getTimeSeriesByRegions(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                  collection = fc,
                                  bands = ['EVI','NDVI'],
                                  scale = 10)

The same parameters of the :code:`getTimeSeriesByRegion()` method can be used here, except for :code:`bestEffort` and :code:`maxPixels`:

.. code-block:: python

   ts = S2.getTimeSeriesByRegions(reducer = [ee.Reducer.mean(),ee.Reducer.median()],
                                  collection = fc,
                                  bands = ['EVI','NDVI'],
                                  scale = 10,
                                  naValue = -99999999,
                                  dateColumn = 'my_date_colum',
                                  dateFormat = 'ms')
                                  
Conversion to Pandas
~~~~~~~~~~~~~~~~~~~~~~~~

The time series is always retrieved as an ee.FeatureCollection. To convert the collection to a pandas data frame, we'll use the geemap package:

.. code-block:: python

   tsPandas = geemap.ee_to_pandas(ts)
   
Then we can convert the NA value to a real NA and the date column to a datetime class:

.. code-block:: python

   tsPandas[tsPandas == -9999] = np.nan
   tsPandas['date'] = pd.to_datetime(tsPandas['date'],infer_datetime_format = True)Closest Image to a Specific Date
====================================

Let's see how to get the closest image (or set of images) to a specific date.

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.ImageCollection class with the method :code:`closest()`:

.. currentmodule:: eemont.imagecollection

.. autosummary::

   closest

This method automatically filters any image collection to get the closest image to a specific date.

.. warning::
   This method uses the :code:`system:time_start` property, therefore, make sure your image collection has it!   

Usage
------------------

The :code:`closest()` method works on any image colection that has a :code:`system:time_start` property.

First, let's take the Sentinel-2 image collection as example:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR')
   
Now, we have to filter the collection to our ROI. The result of the :code:`closest()` method will vary depending on this.
Let's assume a single point is our ROI.

.. code-block:: python

   ROI = ee.Geometry.Point([-76.45, 4.32])
   S2 = S2.filterBounds(ROI)

Now, the :code:`closest()` method has just one parameter, :code:`date`, and this parameter can be a string...

.. code-block:: python

   S2.closest('2020-10-15')
   
Or an ee.Date class:

.. code-block:: python

   dateOfInterest = ee.Date('2020-10-15')
   S2.closest(dateOfInterest)
   
Both chunks will give you the same result here: an ee.ImageCollection of size 1. The result has just one image since our ROI intersects just one scene.
To get that image as a single image, we can use the :code:`first()` method.

.. code-block:: python

   S2.closest('2020-10-15').first()

By default, the image collection is filtered according to +/- 1 month from the :code:`date` parameter (:code:`tolerance = 1` and :code:`unit = 'month'`). This is done to speed up the searching process, but if required (if there are not images in that range), the :code:`tolerance` and :code:`unit` parameters can be modified:

.. code-block:: python

   S2.closest('2020-10-15', tolerance = 2, unit = 'year')

Now, let's assume that our ROI is larger, in this case, a whole department (state) of Colombia:

.. code-block:: python

   ROI = (ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level1')
       .filter(ee.Filter.eq('ADM1_NAME','Valle Del Cauca')))
   S2 = ee.ImageCollection('COPERNICUS/S2_SR').filterBounds(ROI).closest('2020-10-15')

You'll note that the size of the resulting ee.ImageCollection here is greater than 1. This result has more than one image since our ROI now intersects more than one scene.
To get those images together as a single image, you can mosaic them or use an ee.Reducer, for example :code:`median()`.

.. code-block:: python
   
   S2.median()Image Scaling
====================================

Image scaling now is A LOT EASIER with eemont! Let's see how!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the following methods:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

   getOffsetParams
   getScaleParams
   preprocess
   scaleAndOffset
      
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   getOffsetParams
   getScaleParams
   preprocess
   scaleAndOffset

Supported Platforms
----------------------

This method automatically scales ALL images from the `Google Earth Engine STAC Catalog <https://developers.google.com/earth-engine/datasets>`_ by using the `List of Scale and Offset Paramaters from the GEE STAC Catalog Repository <https://github.com/davemlz/ee-catalog-scale-offset-params>`_.

Usage
------------------

The :code:`scaleAndOffset()` method scales each image according to the *Scale* and *Offset* parameters for each band of the image. 

Let's take the Sentinel-2 SR image collection as example:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR')
   
The spectral bands from Sentinel-2 have values close to the (0, 10000) range, but they're unscaled. In order to get the real values, each spectral band must be multiplied by 0.0001, while AOT and WVP bands must be multiplied by 0.001. This scaling is automatically done by the :code:`scaleAndOffset()` method, without any additional parameters:

.. code-block:: python

   S2.scaleAndOffset()

The *Scale* and *Offset* parameters vary according to the satellite platform and the :code:`scaleAndOffset()` method detects the platform and do the scaling according to its parameters.

Let's take now the MOD11A2 product from MODIS. The LST_Day_1km and LST_Night_1km bands must be multiplied by 0.02, the Day_view_time and Night_view_time bands must be multiplied by 0.1,
the Emis_31 and Emis_32 bands must be multiplied by 0.002 and added by 0.49, while the Day_view_angl and Night_view_angl bands must be added by -65. All of this scaling
is simply done by the :code:`scaleAndOffset()` method:

.. code-block:: python

   MOD11A2scaled = ee.ImageCollection('MODIS/006/MOD11A2').scaleAndOffset()
   
The :code:`scale()` method can be applied to single images as well:

.. code-block:: python

   MOD11A2scaled = ee.ImageCollection('MODIS/006/MOD11A2').first().scaleAndOffset()
   
Preprocessing
~~~~~~~~~~~~~~~~~
   
Additionally, if required, the :code:`preprocess()` method can be used to scale and mask images and image collections:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').preprocess()Data Conversion
====================================

Let's see how to convert non-Earth Engine classes to Earth Engine classes.

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   import pandas as pd
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the pd.DataFrame classes with the method :code:`toEEFeatureCollection()`:

pd.DataFrame
~~~~~~~~

.. currentmodule:: eemont.dataframe

.. autosummary::

   toEEFeatureCollection

Methods
-----------

A table of availabe conversion options is shown below:

.. list-table:: Available options
   :widths: 30 30 40
   :header-rows: 1

   * - From
     - To
     - Method
   * - pd.DataFrame
     - ee.FeatureCollection
     - :code:`toEEFeatureCollection()`

Usage
------------------

Let's create a pandas data frame:

.. code-block:: python

   df = pd.DataFrame()
   df['lat'] = [2.92846, 4.8927]
   df['lon'] = [-76.0269, -75.3188]
   df['name'] = ['Nevado del Huila', 'Nevado del Ruiz']
   
This data frame can be easily converted into a ee.FeatureCollection (with no geometries) using the :code:`toEEFeatureCollection()` method for pd.DataFrame classes:

.. code-block:: python

   fcWithNoGeometries = df.toEEFeatureCollection()

If the data frame has latitude and longitude columns, these can be specified in the :code:`latitude` and :code:`longitude` parameters:

.. code-block:: python

   fcWithGeometries = df.toEEFeatureCollection(latitude = 'lat',longitude = 'lon')User Guide
===========

Here you can find the user guide for the extended methods:

.. toctree::
   :caption: Table of Contents
   :maxdepth: 3

   closestImage
   constructors
   dataConversion
   containers
   histogramMatching
   imageScaling
   maskingClouds
   overloadedOperators
   panSharpening
   spectralIndices   
   timeSeriesSpectral Indices Computation
====================================

Let's see how to compute built-in spectral indices with eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the method :code:`spectralIndices()`:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

   spectralIndices
      
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   spectralIndices

Supported Platforms
----------------------

This method automatically computes spectral indices for the following supported satellite platforms:

Sentinel Missions
~~~~~~~~~~~~~~~~~~~

- `Sentinel-2 MSI: MultiSpectral Instrument, Level-2A <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR?hl=en>`_
- `Sentinel-2 MSI: MultiSpectral Instrument, Level-1C <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2#bands>`_

Landsat Missions
~~~~~~~~~~~~~~~~~~~

- `USGS Landsat 8 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR>`_
- `USGS Landsat 8 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2>`_
- `USGS Landsat 7 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_SR>`_
- `USGS Landsat 7 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2>`_
- `USGS Landsat 5 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1_SR>`_
- `USGS Landsat 4 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C01_T1_SR>`_

MODIS Products (Terra + Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MCD43A4.006 MODIS Nadir BRDF-Adjusted Reflectance Daily 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD43A4>`_

MODIS Products (Terra)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MOD09GQ.006 Terra Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GQ>`_
- `MOD09GA.006 Terra Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GA>`_
- `MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09Q1>`_
- `MOD09A1.006 Terra Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1>`_

MODIS Products (Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MYD09GQ.006 Aqua Surface Reflectance Daily Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GQ>`_
- `MYD09GA.006 Aqua Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GA>`_
- `MYD09Q1.006 Aqua Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09Q1>`_
- `MYD09A1.006 Aqua Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09A1>`_

.. important::
   It is highly recommended to scale the image (or image collection) before computing spectral indices. See the :code:`scaleAndOffset()` method for more info.  

List of Indices
----------------------

The list of indices is retrieved from the `Awesome List of Spectral Indices for Google Earth Engine <https://github.com/davemlz/awesome-ee-spectral-indices>`_

Vegetation Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in vegetation indices:

.. list-table:: Built-in vegetation indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - BNDVI
     - Blue Normalized Difference Vegetation Index
     - `Index DataBase BNDVI <https://www.indexdatabase.de/db/i-single.php?id=135>`_  
   * - CIG
     - Chlorophyll Index - Green
     - `Index DataBase CIG <https://www.indexdatabase.de/db/i-single.php?id=128>`_
   * - CVI
     - Chlorophyll Vegetation Index
     - `Index DataBase CVI <https://www.indexdatabase.de/db/i-single.php?id=391>`_
   * - EVI
     - Enhanced Vegetation Index
     - `Index DataBase EVI <https://www.indexdatabase.de/db/i-single.php?id=16>`_
   * - EVI2
     - Two-Band Enhanced Vegetation Index
     - `(Jiang et al., 2008) <https://doi.org/10.1016/j.rse.2008.06.006>`_
   * - GARI
     - Green Atmospherically Resistant Vegetation Index
     - `Index DataBase GARI <https://www.indexdatabase.de/db/i-single.php?id=363>`_
   * - GBNDVI
     - Green-Blue Normalized Difference Vegetation Index
     - `Index DataBase GBNDVI <https://www.indexdatabase.de/db/i-single.php?id=186>`_
   * - GEMI
     - Global Environment Monitoring Index
     - `Index DataBase GEMI <https://www.indexdatabase.de/db/i-single.php?id=25>`_
   * - GLI
     - Green Leaf Index
     - `Index DataBase GLI <https://www.indexdatabase.de/db/i-single.php?id=375>`_
   * - GNDVI
     - Green Normalized Difference Vegetation Index
     - `Index DataBase GNDVI <https://www.indexdatabase.de/db/i-single.php?id=401>`_
   * - GRNDVI
     - Green-Red Normalized Difference Vegetation Index
     - `Index DataBase GRNDVI <https://www.indexdatabase.de/db/i-single.php?id=185>`_
   * - GVMI
     - Global Vegetation Moisture Index
     - `Index DataBase GVMI <https://www.indexdatabase.de/db/i-single.php?id=372>`_
   * - MNDVI
     - Modified Normalized Difference Vegetation Index
     - `Index DataBase MNDVI <https://www.indexdatabase.de/db/i-single.php?id=245>`_
   * - NDVI
     - Normalized Difference Vegetation Index
     - `Index DataBase NDVI <https://www.indexdatabase.de/db/i-single.php?id=58>`_
   * - NGRDI
     - Normalized Green Red Difference Index
     - `Index DataBase NGRDI <https://www.indexdatabase.de/db/i-single.php?id=390>`_
   * - RVI
     - Ratio Vegetation Index
     - `Index DataBase RVI <https://www.indexdatabase.de/db/i-single.php?id=72>`_
   * - SAVI
     - Soil-Adjusted Vegetation Index
     - `Index DataBase SAVI <https://www.indexdatabase.de/db/i-single.php?id=87>`_
   * - VARI
     - Visible Atmospherically Resistant Index
     - `Index DataBase VARI <https://www.indexdatabase.de/db/i-single.php?id=356>`_     
     
Burn Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in burn indices:

.. list-table:: Built-in burn indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - BAI
     - Burned Area Index
     - `(Martín, 1998) [spanish] <https://digital.csic.es/bitstream/10261/6426/1/Martin_Isabel_Serie_Geografica.pdf>`_ 
   * - BAIS2
     - Burned Area Index for Sentinel 2
     - `(Filipponi, 2018) <https://doi.org/10.3390/ecrs-2-05177>`_
   * - CSIT
     - Char Soil Index Thermal
     - `(Smith et al., 2007) <https://doi.org/10.1080/01431160600954704>`_
   * - NBR
     - Normalized Burn Ratio
     - `Index DataBase NBR <https://www.indexdatabase.de/db/i-single.php?id=53>`_
   * - NBRT
     - Normalized Burn Ratio Thermal
     - `(Holden et al., 2005) <https://doi.org/10.1080/01431160500239008>`_
   * - NDVIT
     - Normalized Difference Vegetation Index Thermal
     - `(Smith et al., 2007) <https://doi.org/10.1080/01431160600954704>`_
   * - SAVIT
     - Soil-Adjusted Vegetation Index Thermal
     - `(Smith et al., 2007) <https://doi.org/10.1080/01431160600954704>`_
     
Water Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in water indices:

.. list-table:: Built-in water indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - MNDWI
     - Modified Normalized Difference Water Index
     - `(Xu, 2006) <https://doi.org/10.1080/01431160600589179>`_  
   * - NDWI
     - Normalized Difference Water Index
     - `(McFeeters, 1996) <https://doi.org/10.1080/01431169608948714>`_
     
Snow Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in snow indices:

.. list-table:: Built-in snow indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - NDSI
     - Normalized Difference Snow Index
     - `(Riggs et al., 1994) <https://doi.org/10.1109/IGARSS.1994.399618>`_
     
Drought Indices
~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in drought indices:

.. list-table:: Built-in snow indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - NDDI
     - Normalized Difference Drought Index
     - `(Gu et al., 2007) <https://doi.org/10.1029/2006GL029127>`_
     
Generalized Kernel Indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table shows the list of built-in kernel indices:

.. list-table:: Built-in kernel indices.
   :widths: 20 50 30
   :header-rows: 1

   * - Index
     - Description     
     - Reference
   * - kEVI
     - Kernel Enhanced Vegetation Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_
   * - kNDVI
     - Kernel Normalized Difference Vegetation Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_
   * - kRVI
     - Kernel Ratio Vegetation Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_
   * - kVARI
     - Kernel Visible Atmospherically Resistant Index
     - `(Camps-Valls et al., 2021) <https://doi.org/10.1126/sciadv.abc7447>`_

Kernels
----------------------

In the case of generalized kernel indices, the following kernels are available:

Linear Kernel
~~~~~~~~~~~~~~~~~~~

The linear kernel for generalized kernel indices can be selected by setting :code:`kernel = 'linear'`.

.. math::

   k(a,b) = ab
   
RBF Kernel
~~~~~~~~~~~~~~~~~~~

The Radial Basis Function (RBF) kernel for generalized kernel indices can be selected by setting :code:`kernel = 'RBF'`.

.. math::

   k(a,b) = exp(- \frac{(a - b) ^ 2}{2 \sigma ^ 2})
   
Where :math:`\sigma` is a free length-scale parameter.
   
Polynomial Kernel
~~~~~~~~~~~~~~~~~~~

The polynomial kernel for generalized kernel indices can be selected by setting :code:`kernel = 'poly'`.

.. math::

   k(a,b) = (ab + c) ^ p

Where :math:`c` is a free parameter that trades off the influence of higher-order versus lower-order terms and :math:`p` is the kernel degree.

List of Bands
----------------------

The following table shows the list of bands used for spectral indices computation:

.. list-table:: Bands used for spectral indices computation.
   :widths: 18 18 16 16 16 16
   :header-rows: 1

   * - Description
     - Name     
     - Sentinel-2
     - Landsat 8
     - Landsat 457
     - MODIS     
   * - Aerosols
     - A
     - B1
     - B1
     -
     -     
   * - Blue
     - B
     - B2
     - B2
     - B1
     - B3 
   * - Green
     - G
     - B3
     - B3
     - B2
     - B4    
   * - Red
     - R
     - B4
     - B4
     - B3
     - B1
   * - Red Edge 1
     - RE1
     - B5
     - 
     -
     -     
   * - Red Edge 2
     - RE2
     - B6
     - 
     -
     -     
   * - Red Edge 3
     - RE3
     - B7
     - 
     -
     -     
   * - Red Edge 4
     - RE4
     - B8A
     - 
     -
     -     
   * - NIR
     - N
     - B8
     - B5
     - B4
     - B2
   * - SWIR 1
     - S1
     - B11
     - B6
     - B5
     - B6     
   * - SWIR 2
     - S2
     - B12
     - B7
     - B7
     - B7   
   * - Thermal 1
     - T1
     - 
     - B10
     - B6
     -     
   * - Thermal 2
     - T2
     - 
     - B11
     - 
     -     

.. warning::
   If the satellite platform doesn't have the required bands for computing an index, it won't be computed.

Usage
------------------

The :code:`spectralIndices()` method computes the specified spectral index and adds it as a new band.

Let's take the Sentinel-2 SR image collection as example (remember to scale your image or image collection!):

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').scaleAndOffset()
   
By default, the :code:`spectralIndices()` method computes the NDVI:

.. code-block:: python

   S2withIndices = S2.spectralIndices()
   S2withIndices.select('NDVI')
   
If required, any of the above-mentioned indices can be computed by modifying the :code:`index` parameter:

.. code-block:: python

   S2withIndices = S2.spectralIndices(index = 'EVI')
   S2withIndices.select('EVI')
   
Specific index-parameters can be changed, for example, the canopy background adjustment L is set to 1.0 for EVI, but for SAVI it can be changed to 0.5:

.. code-block:: python

   S2withIndices = S2.spectralIndices('SAVI',L = 0.5)
   S2withIndices.select('SAVI')
   
If more than one index is required, a list of indices can be used:

.. code-block:: python

   S2withIndices = S2.spectralIndices(['CIG','NBR','NDWI'])
   S2withIndices.select('CIG')
   S2withIndices.select('NBR')
   S2withIndices.select('NDWI')
   
Indices can also be computed for single images:

.. code-block:: python

   S2withIndices = S2.first().spectralIndices(['GBNDVI','MNDVI','EVI'])
   S2withIndices.select('GBNDVI')
   S2withIndices.select('MNDVI')
   S2withIndices.select('EVI')
   
All vegetation indices can be computed by setting :code:`index = vegetation`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('vegetation')
   S2withIndices.select('NDVI')
   S2withIndices.select('GNDVI')
   S2withIndices.select('RVI')
   # ...
   
All burn indices can be computed by setting :code:`index = burn`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('burn')
   S2withIndices.select('BAI')
   S2withIndices.select('BAIS2')
   S2withIndices.select('NBR')
   
All water indices can be computed by setting :code:`index = water`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('water')
   S2withIndices.select('NDWI')
   S2withIndices.select('MNDWI')
   
All snow indices can be computed by setting :code:`index = snow`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('snow')
   S2withIndices.select('NDSI')
   
If you want to compute all available indices, you can set :code:`index = all`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('all')
   S2withIndices.select('NDVI')
   S2withIndices.select('BAI')
   S2withIndices.select('NDWI')
   S2withIndices.select('NDSI')
   # ...

Generalized Kernel Indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generalized kernel indices are availabe through eemont (e.g. kNDVI):

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI')
   S2withIndices.select('kNDVI')
   
By default, the RBF kernel is used and the :code:`sigma` parameter is :code:`0.5 * (a + b)` (this means, that for :code:`k(N,R)`, :code:`sigma = '0.5 * (N + R)'`). If required, :code:`sigma` can be modified by another expression (using :code:`a` and :code:`b`) or a float:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI',sigma = 1)
   S2withIndices.select('kNDVI')
   
The kernel can be modified by modifying the :code:`kernel` parameter:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI',kernel = 'poly')
   S2withIndices.select('kNDVI')
   
For the polynomial kernel, the :code:`p` and :code:`c` parameters can be modified:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kNDVI',kernel = 'poly',p = 4,c = 0)
   S2withIndices.select('kNDVI')
   
All kernel indices can be computed by setting :code:`index = kernel`:

.. code-block:: python

   S2withIndices = S2.spectralIndices('kernel')
   S2withIndices.select('kEVI')
   S2withIndices.select('kNDVI')
   S2withIndices.select('kRVI')
   S2withIndices.select('kVARI')
   
.. seealso::
   For more info on generalized kernel indices, please visit
   `‘Camps-Valls, G., et al. 2021. A unified vegetation index for quantifying the terrestrial biosphere. Science Advances 7 (9): eabc7447. Doi: 10.1126/sciadv.abc7447’ <https://doi.org/10.1126/sciadv.abc7447>`_.Masking Clouds and Shadows
====================================

Masking clouds and shadows may seem hard, but it isn't! Let's take a look on it!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.ImageCollection classes with the following methods:

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::

   maskClouds
   preprocess
      
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection

.. autosummary::

   maskClouds
   preprocess

Supported Platforms
----------------------

This method automatically masks clouds and shadows on the following supported satellite platforms:

Sentinel Missions
~~~~~~~~~~~~~~~~~~~

- `Sentinel-3 OLCI EFR: Ocean and Land Color Instrument Earth Observation Full Resolution <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S3_OLCI>`_
- `Sentinel-2 MSI: MultiSpectral Instrument, Level-2A <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR?hl=en>`_

Landsat Missions
~~~~~~~~~~~~~~~~~~~

- `USGS Landsat 8 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR>`_
- `USGS Landsat 8 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2>`_
- `USGS Landsat 7 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_SR>`_
- `USGS Landsat 7 Level 2, Collection 2, Tier 1 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2>`_
- `USGS Landsat 5 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1_SR>`_
- `USGS Landsat 4 Surface Reflectance Tier 1 and 2 <https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C01_T1_SR>`_

MODIS Products (Terra + Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MCD15A3H.006 MODIS Leaf Area Index/FPAR 4-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD15A3H>`_

MODIS Products (Terra)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MOD09GA.006 Terra Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09GA>`_
- `MOD09Q1.006 Terra Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09Q1>`_
- `MOD09A1.006 Terra Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1>`_
- `MOD17A2H.006: Terra Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD17A2H>`_
- `MOD16A2.006: Terra Net Evapotranspiration 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2>`_
- `MOD13Q1.006 Terra Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13Q1>`_
- `MOD13A1.006 Terra Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A1>`_
- `MOD13A2.006 Terra Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A2>`_

MODIS Products (Aqua)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MYD09GA.006 Aqua Surface Reflectance Daily Global 1km and 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09GA>`_
- `MYD09Q1.006 Aqua Surface Reflectance 8-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09Q1>`_
- `MYD09A1.006 Aqua Surface Reflectance 8-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD09A1>`_   
- `MYD17A2H.006: Aqua Gross Primary Productivity 8-Day Global 500M 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD17A2H>`_   
- `MYD13Q1.006 Aqua Vegetation Indices 16-Day Global 250m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13Q1>`_
- `MYD13A1.006 Aqua Vegetation Indices 16-Day Global 500m <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A1>`_
- `MYD13A2.006 Aqua Vegetation Indices 16-Day Global 1km <https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MYD13A2>`_

VIIRS Products
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `VNP09GA: VIIRS Surface Reflectance Daily 500m and 1km <https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_001_VNP09GA?hl=en>`_
- `VNP13A1: VIIRS Vegetation Indices 16-Day 500m <https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_001_VNP13A1?hl=en>`_

.. warning::
   Not supported satellite platforms will raise a *Warning*.   

QA Method
----------------------

By default, the :code:`maskClouds()` uses the QA band of each paltform to compute the clouds and shadows masks (except for Sentinel-2, where the default method is Cloud Probability). The following table shows the band and the bits used for each platform (The value in parentheses is the valid value of the bitmask):

.. list-table:: QA bits used for clouds/shadows masking
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Platform
     - QA Band
     - Cloud Bitmask
     - Cirrus Bitmask
     - Shadow Bitmask
   * - Sentinel-3
     - quality_flags
     - 27 (0)
     -
     -
   * - Sentinel-2
     - QA60
     - 10 (0)
     - 11 (0)
     -
   * - Landsat Series
     - pixel_qa
     - 5 (0)
     - 
     - 3 (0)
   * - Landsat Collection 2
     - QA_PIXEL
     - 3 (0)
     - 2 (0)
     - 4 (0)
   * - MOD09GA
     - state_1km
     - 0 (0)
     - 8 (0)
     - 2 (0)
   * - MOD09Q1
     - State
     - 0 (0)
     - 8 (0)
     - 2 (0)
   * - MOD09A1
     - StateQA
     - 0 (0)
     - 8 (0)
     - 2 (0)
   * - MCD15A3H
     - FparExtra_QC
     - 5 (0)
     - 4 (0)
     - 6 (0)
   * - MOD17A2H
     - Psn_QC
     - 3 (0)
     - 
     - 
   * - MOD16A2
     - ET_QC
     - 3 (0)
     - 
     - 
   * - MOD13Q1
     - SummaryQA
     - 0 (0)
     - 
     - 
   * - MOD13A1
     - SummaryQA
     - 0 (0)
     - 
     -
   * - VNP09GA
     - QF1, QF2
     - 2 (0)
     - 6 (0), 7 (0)
     - 3 (0)

Usage
-----

Let's check how to use the :code:`maskClouds()` method for different platforms:

Sentinel-3
~~~~~~~~~~~~~~

On Sentinel 3, clouds are masked according to the bright pixels in the quality_flags band of the `Sentinel-3 OLCI EFR: Ocean and Land Color Instrument Earth Observation Full Resolution <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S3_OLCI>`_.

.. warning::
   This method may mask water as well on Sentinel-3 images. 

Let's take the Sentinel-3 image collection as example:

.. code-block:: python

   S3 = ee.ImageCollection('COPERNICUS/S3/OLCI')
   
There is no need to specify any arguments, since they're ignored.

.. code-block:: python

   S3.maskClouds()
  
This method can also be applied to a single image:

.. code-block:: python

   S3.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   S3.scale().maskClouds()

Sentinel-2
~~~~~~~~~~~~~~~~

On Sentinel 2, clouds can be masked using two methods: *QA* and *Cloud Probability*. The *QA* method uses the QA60 band in the `Surface Reflectance Product <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR?hl=en>`_ to mask clouds, while the *Cloud Probability* method uses the
`COPERNICUS/S2_CLOUD_PROBABILITY <https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_CLOUD_PROBABILITY?hl=en>`_ collection to do it.

Shadows are masked based on the clouds mask, where shadows are searched within a range from clouds edges in the shadows direction.

.. seealso::
   For more info on masking shadows, please visit
   `‘Braaten, J. 2020. Sentinel-2 Cloud Masking with s2cloudless. Google Earth Engine, Community Tutorials’ <https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless>`_.

First, let's take the Sentinel-2 image collection:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR')

In order to use the *QA* method, it must be specified using the :code:`method` parameter:

.. code-block:: python

   S2.maskClouds(method = 'qa')
   
This line maps the *QA* masking method over the whole collection, but the method can also be applied to a single image:

.. code-block:: python

   S2.first().maskClouds(method = 'qa')
   
The *QA* method gives us the option to avoid masking cirrus clouds, but it must be specified using the :code:`maskCirrus` parameter:

.. code-block:: python

   S2.maskClouds(method = 'qa', maskCirrus = False)

And we can also avoid masking shadows by specifying the :code:`maskShadows` parameter:

.. code-block:: python

   S2.maskClouds(method = 'qa', maskShadows = False)
   
Now, in order to use the *Cloud Probability* method, we can specify it in the :code:`method` parameter:

.. code-block:: python

   S2.maskClouds(method = 'cloud_prob')
   
But, it is the default method, so you can just let the extended method with no additional parameters:

.. code-block:: python

   S2.maskClouds()
   
The *Cloud Probability* method uses a probability threshold to mask clouds, by default, the threshold is set to 60, but it can be modified using the :code:`prob` parameter:

.. code-block:: python

   S2.maskClouds(prob = 70)
   
If your image or collection is scaled, the :code:`scaledImage` parameter must be set to :code:`True`:

.. code-block:: python

   S2.scale().maskClouds(scaledImage = True)
   
In order to search for shadows, portental shadow pixels must be specified. Pixels with a NIR reflectance below 0.15 are considered potential shadow pixels, but this can be modified using the
:code:`dark` parameter:

.. code-block:: python

   S2.maskClouds(dark = 0.2)
   
Shadows are searched whitin a maximum range of 1000 m in the shadow direction from cloud edges, but this range can be modified using the :code:`cloudDist` parameter:

.. code-block:: python

   S2.maskClouds(cloudDist = 1500)
   
After finding all clouds and shadows, the mask can be dilated to avoid border effects. By default, clouds and shadows are dilated by 250 m, but this can be modified using the :code:`buffer` parameter:

.. code-block:: python

   S2.maskClouds(buffer = 100)
   
Finally, in order to avoid confusion between clouds and bright surface objects, the Cloud Displacement Index (CDI) can be used. By default, the CDI is not used, but it can be modified
using the :code:`cdi` parameter:

.. code-block:: python

   S2.maskClouds(cdi = -0.5)
   
.. seealso::
   For more info on CDI, please visit
   `‘Frantz, D., HaS, E., Uhl, A., Stoffels, J., Hill, J. 2018. Improvement of the Fmask algorithm for Sentinel-2 images: 
   Separating clouds from bright surfaces based on parallax effects. Remote Sensing of Environment 2015: 471-481’ 
   <https://www.sciencedirect.com/science/article/pii/S0034425718302037#:~:text=In%20this%20paper%2C%20we%20present,separated%20from%20bright%20ground%20objects.>`_.

Landsat Series
~~~~~~~~~~~~~~~~

On Landsat Series, both clouds and shadows are masked based on the pixel_qa band in the `Surface Reflectance Products <https://developers.google.com/earth-engine/datasets/catalog/landsat>`_.

Let's take the Landsat 8 image collection as example:

.. code-block:: python

   L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
   
There is no need to specify most of the arguments showed for Sentinel-2, since they're ignored.

.. code-block:: python

   L8.maskClouds()
   
Shadows are masked by default, but if required, the :code:`maskShadows` parameter can be modified.

.. code-block:: python

   L8.maskClouds(maskShadows = False)
   
This method can also be applied to a single image:

.. code-block:: python

   L8.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   L8.scale().maskClouds()
   
MODIS Products
~~~~~~~~~~~~~~~~

On MODIS Products, clouds and shadows are masked according to the specific QA band.

Let's take the MOD13Q1 image collection as example:

.. code-block:: python

   MOD13Q1 = ee.ImageCollection('MODIS/006/MOD13Q1')
   
There is no need to specify most of the arguments showed for Sentinel-2, since they're ignored.

.. code-block:: python

   MOD13Q1.maskClouds()
   
MOD13Q1, MOD13A1, MOD17A2H and MOD16A2 products don't have cirrus and shadow bitmasks, therefore, the arguments :code:`maskShadows` and :code:`maskCirrus` are ignored. MOD09GA, MOD09Q1, MOD09A1 and MCD15A3H products have cirrus and shadows bitmasks, and by default, they are set to *True*. If required, they can be set to *False*:

.. code-block:: python

   MOD09GA = ee.ImageCollection('MODIS/006/MOD09GA').maskClouds(maskShadows = False, maskCirrus = False)
   
This method can also be applied to a single image:

.. code-block:: python

   MOD09GA.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   MOD09GA.scale().maskClouds()
   
MOD13A2 doesn't have a bitmask QA band, instead, it has a Class QA band, where a value of zero means that the pixel has good data.

.. code-block:: python

   MOD13A2 = ee.ImageCollection('MODIS/006/MOD13A2').maskClouds()
   
VIIRS Products
~~~~~~~~~~~~~~~~

On VIIRS Products, clouds and shadows are masked according to the specific QA band.

Let's take the VNP09GA image collection as example:

.. code-block:: python

   VNP09GA = ee.ImageCollection('NOAA/VIIRS/001/VNP13A1')
   
There is no need to specify most of the arguments showed for Sentinel-2, since they're ignored.

.. code-block:: python

   VNP09GA.maskClouds()
   
If required, the arguments :code:`maskShadows` and :code:`maskCirrus` can be set to *False*:

.. code-block:: python

   VNP09GA.maskClouds(maskShadows = False, maskCirrus = False)
   
This method can also be applied to a single image:

.. code-block:: python

   VNP09GA.first().maskClouds()
   
And can be used for scaled images without specifying it:

.. code-block:: python

   VNP09GA.scale().maskClouds()
   
VNP13A1 doesn't have a bitmask QA band, instead, it has a Class QA band, where a value of 9 means that the pixel is a cloud, while a value of 7 means that the pixel is a cloud shadows.

.. code-block:: python

   VNP13A1 = ee.ImageCollection('NOAA/VIIRS/001/VNP13A1').maskClouds()
   
If required, the argument :code:`maskShadows` can be set to *False*:

.. code-block:: python

   VNP13A1 = ee.ImageCollection('NOAA/VIIRS/001/VNP13A1').maskClouds(maskShadows = False)
   
Preprocessing
~~~~~~~~~~~~~~~
   
Additionally, if required, the :code:`preprocess()` method can be used to scale and mask images and image collections:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').preprocess()
   
The arguments for the :code:`maskClouds()` method can be used inside the :code:`preprocess()` method:

.. code-block:: python

   S2 = ee.ImageCollection('COPERNICUS/S2_SR').preprocess(prob = 55,maskCirrus = False,buffer = 300,cdi = -0.5)Overloaded Operators
====================================

Let's see how to use overloaded operators in eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Image and ee.Number classes with the binary and unary operators (including rich comparisons).

ee.Image 
-------------------

Binary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of binary operators that are overloaded:

.. list-table:: Binary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Overloaded Operator
   * - Addition
     - Image1.add(Image2)
     - Image1 + Image2 
   * - Subtraction
     - Image1.subtract(Image2)
     - Image1 - Image2
   * - Multiplication
     - Image1.multiply(Image2)
     - Image1 * Image2
   * - Division
     - Image1.divide(Image2)
     - Image1 / Image2
   * - Floor Division
     - Image1.divide(Image2).floor()
     - Image1 // Image2
   * - Modulo
     - Image1.mod(Image2)
     - Image1 % Image2
   * - Power
     - Image1.pow(Image2)
     - Image1 ** Image2
   * - Left Shift
     - Image1.leftShift(Image2)
     - Image1 << Image2
   * - Right Shift
     - Image1.rightShift(Image2)
     - Image1 >> Image2
   * - And
     - Image1.And(Image2)
     - Image1 & Image2
   * - Or
     - Image1.Or(Image2)
     - Image1 | Image2
          
Rich Comparisons
~~~~~~~~~~~~~~~~~~~

The following table shows the list of rich comparisons that are overloaded:

.. list-table:: Rich comparisons.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Lower Than
     - Image1.lt(Image2)
     - Image1 < Image2 
   * - Lower Than or Equal
     - Image1.lte(Image2)
     - Image1 <= Image2
   * - Equal
     - Image1.eq(Image2)
     - Image1 == Image2
   * - Not Equal
     - Image1.neq(Image2)    
     - Image1 != Image2
   * - Greater Than
     - Image1.gt(Image2)
     - Image1 > Image2 
   * - Greater Than or Equal
     - Image1.gte(Image2)
     - Image1 >= Image2
   * - Equal
     - Image1.eq(Image2)
     - Image1 == Image2
     
Unary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of unary operators that are overloaded:

.. list-table:: Unary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Negation
     - Image.multiply(-1)
     - \-\ Image
   * - Invert
     - Image.Not()
     - ~ Image
     
ee.Number 
-------------------

Binary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of binary operators that are overloaded:

.. list-table:: Binary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Overloaded Operator
   * - Addition
     - Number1.add(Number2)
     - Number1 + Number2 
   * - Subtraction
     - Number1.subtract(Number2)
     - Number1 - Number2
   * - Multiplication
     - Number1.multiply(Number2)
     - Number1 * Number2
   * - Division
     - Number1.divide(Number2)
     - Image1 / Image2
   * - Floor Division
     - Number1.divide(Number2).floor()
     - Number1 // Number2
   * - Modulo
     - Number1.mod(Number2)
     - Number1 % Number2
   * - Power
     - Number1.pow(Number2)
     - Number1 ** Number2
   * - Left Shift
     - Number1.leftShift(Number2)
     - Number1 << Number2
   * - Right Shift
     - Number1.rightShift(Number2)
     - Number1 >> Number2
   * - And
     - Number1.And(Number2)
     - Number1 & Number2
   * - Or
     - Number1.Or(Number2)
     - Number1 | Number2
          
Rich Comparisons
~~~~~~~~~~~~~~~~~~~

The following table shows the list of rich comparisons that are overloaded:

.. list-table:: Rich comparisons.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Lower Than
     - Number1.lt(Number2)
     - Number1 < Number2 
   * - Lower Than or Equal
     - Number1.lte(Number2)
     - Number1 <= Number2
   * - Equal
     - Number1.eq(Number2)
     - Number1 == Number2
   * - Not Equal
     - Number1.neq(Number2)    
     - Number1 != Number2
   * - Greater Than
     - Number1.gt(Number2)
     - Number1 > Number2 
   * - Greater Than or Equal
     - Number1.gte(Number2)
     - Number1 >= Number2
   * - Equal
     - Number1.eq(Number2)
     - Number1 == Number2
     
Unary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of unary operators that are overloaded:

.. list-table:: Unary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method          
     - Overloaded Operator
   * - Negation
     - Number.multiply(-1)
     - \-\ Number
   * - Invert
     - Number.Not()
     - ~ Number
     
ee.List 
-------------------

Binary Operators
~~~~~~~~~~~~~~~~~~~

The following table shows the list of binary operators that are overloaded:

.. list-table:: Binary operators.
   :widths: 20 40 40
   :header-rows: 1

   * - Operation
     - GEE Python method     
     - Overloaded Operator
   * - Concatenation
     - List1.cat(List2)
     - List1 + List2
   * - Repeat
     - ee.List.repeat(List,Value)
     - List * Value

Usage
------------------

Overloaded operators can be used on any ee.Image or ee.Number object. Let's see how to compute the EVI using overloaded operators!

Let's take the Sentinel-2 SR image collection as example (remember to scale your image or image collection!):

.. code-block:: python

   point = ee.Geometry.Point([-76.0269,2.92846])
   S2 = (ee.ImageCollection('COPERNICUS/S2_SR')
      .filterBounds(point)
      .sort('CLOUDY_PIXEL_PERCENTAGE')
      .first()
      .maskClouds()
      .scale())

Now, let's take apart the bands that we need (it is not necessary, but it's easier to use :code:`N` instead of :code:`S2.select('B8')`):

.. code-block:: python

   N = S2.select('B8')
   R = S2.select('B4')
   B = S2.select('B2')
   
And finally, let's compute the EVI using overloaded operators:

.. code-block:: python

   EVI = 2.5 * (N - R) / (N + 6.0 * R - 7.5 * B + 1.0)

Let's see another example, but using rich comparisons. We are going to compute a snow cover mask!

First, compute the NDSI:

.. code-block:: python

   S2 = S2.index('NDSI')   
   
And now, let's take apart the bands that we need:

.. code-block:: python

   NDSI = S2.select('NDSI')
   N = S2.select('B8')
   G = S2.select('B3')
   
Finally, compute the snow cover mask `(Hall et al., 2001) <https://modis.gsfc.nasa.gov/data/atbd/atbd_mod10.pdf>`_:

.. code-block:: python

   snowPixels = (NDSI > 0.4) & (N >= 0.1) & (G > 0.11)

And update the mask (if required):

.. code-block:: python

   S2 = S2.updateMask(snowPixels)Constructors
====================================

Let's see how to use the extended constructors available through eemont!

Before anything, let's import our modules and authenticate in Google Earth Engine:

.. code-block:: python

   import ee, eemont
   
   ee.Authenticate()
   ee.Initialize()

Now, we are ready to go!

Overview
-----------

The eemont package extends the ee.Geometry, ee.Feature and ee.FeatureCollection classes with the following constructors:

ee.Geometry
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.geometry

.. autosummary::

   BBoxFromQuery
   LinearRingFromPlusCodes
   LineStringFromPlusCodes
   MultiLineStringFromPlusCodes
   MultiPointFromPlusCodes
   MultiPointFromQuery
   MultiPolygonFromPlusCodes
   PointFromPlusCode
   PointFromQuery
   PolygonFromPlusCodes
   RectangleFromPlusCodes   

ee.Feature
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.feature

.. autosummary::

   BBoxFromQuery
   PointFromQuery

ee.FeatureCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.featurecollection

.. autosummary::

   MultiPointFromQuery
   
Usage
------------------

Constructors By Query
~~~~~~~~~~~~~~~~~~~~~~~~

The constructors by query use the `geopy package <https://pypi.org/project/geopy/>`_ to construct geometries, features and feature collections.

The simplest geometry to construct is the Point, and can be constructed for ee.Geometry and ee.Feature classes;

.. code-block:: python

   geometry = ee.Geometry.PointFromQuery('Cali, Colombia',user_agent = 'eemont-user-guide-constructors')
   feature = ee.Feature.PointFromQuery('Cali, Colombia',user_agent = 'eemont-user-guide-constructors')
   
It has to be noted that the :code:`user_agent` argument is mandatory and it can be set to any string representing the app or the GEE username.

.. warning::
   If the :code:`user_agent` argument is not defined, the constructor will raise an Exception.

By default, to geocode the query, the Nominatim geocoder is used. If required, this parameter can be modified:

.. code-block:: python

   geometry = ee.Geometry.PointFromQuery('Cali, Colombia',geocoder = 'arcgis',user_agent = 'eemont-user-guide-constructors')
                                 
The second geometry to construct is the MultiPoint, and can be constructed for ee.Geometry and ee.FeatureCollection classes:

.. code-block:: python

   geometry = ee.Geometry.MultiPointFromQuery('Colombia',user_agent = 'eemont-user-guide-constructors')
   feature_collection = ee.FeatureCollection.MultiPointFromQuery('Colombia',user_agent = 'eemont-user-guide-constructors')

.. note::
   When a query is submitted, a set of locations is retrieved. The MultiPoint constructors create a class taking all locations into account. The Point constructors just take the first one.

The last geometry to construct is the Bounding Box, and can be constructed for ee.Geometry and ee.Feature classes:

.. code-block:: python

   geometry = ee.Geometry.BBoxFromQuery('Europe',user_agent = 'eemont-user-guide-constructors')
   feature = ee.Feature.BBoxFromQuery('Europe',user_agent = 'eemont-user-guide-constructors')
   
.. note::
   When using constructors for ee.Feature and ee.FeatureCollection classes, the raw properties of the location, or locations, are set for the feature or feature collection.
   
Constructors By Plus Codes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
   The awesome Plus Codes constructors and methods were created by `Aaron Zuspan <https://github.com/aazuspan>`_ (creator of `sankee <https://github.com/aazuspan/sankee>`_).

`Plus Codes <https://maps.google.com/pluscodes/>`_ are street addresses that represent an area based on longitude and latitude coordinates (e.g. 
"89MH+PW").

.. warning::
   In order to use Plus Codes constructors, it is required to install :code:`openlocationcode`. Please install it by running :code:`pip install openlocationcode`.
   
There are two ways to use the Plus Codes constructors.

1. By using a full Plus Code (e.g. "85FPQXGV+XH").
2. By using a short Plus Code (e.g. "QXGV+XH Denver, CO, USA").

When using full Plus Codes, just use it as the principal argument in the constructor:

.. code-block:: python

   point = ee.Geometry.PointFromPlusCode("85FPQXGV+XH")
   
When using a short Plus Code, it is required to use a geocoder (just like the constructors by queries), and therefore, an user agent must be declared:

.. code-block:: python

   point = ee.Geometry.PointFromPlusCode("QXGV+XH Denver, CO, USA",user_agent = 'eemont-user-guide-constructors')
   
More complex geometries can be constructed using a list of Plus Codes or a nested list of Plus Codes:

.. code-block:: python

   codes = ['85FQ2222+22', '85FR2222+22', '85GR2222+22']
   
   multipoint = ee.Geometry.MultiPointFromPlusCodes(codes)
   linestring = ee.Geometry.LineStringFromPlusCodes(codes)
   polygon = ee.Geometry.PolygonFromPlusCodes(codes)
   
   nestedCodes = [
        ['85FQ2222+22', '85FR2222+22', '85GR2222+22'], 
        ['85FP8PC2+G2', '85FPJF23+G4', '85FPMW2R+RP'],
        ['85FPRJ5W+82', '85GP3M67+CP', '85GQ2R7C+38'],
    ]
    
    multilinestring = ee.Geometry.MultiLineStringFromPlusCodes(nestedCodes)
    multipolygon = ee.Geometry.MultiPolygonFromPlusCodes(nestedCodes)ee.List
=======

Extended methods for the ee.List class:

.. currentmodule:: eemont.eeList

.. autosummary::
   :toctree: stubs

   __add__
   __contains__
   __getitem__
   __len__
   __mul__
   __radd__
   __rmul__ee.Feature
===========

Extended methods for the ee.Feature class:

.. currentmodule:: eemont.feature

.. autosummary::
   :toctree: stubs

   BBoxFromQuery
   PointFromQuery
   plusCodesee
==

Extended methods for the ee module:

..
   .. currentmodule:: eemont.extra

   .. autosummary::
      :toctree: stubs

      install
      require
      uninstall

.. currentmodule:: eemont.app

.. autosummary::
   :toctree: stubs

   listAppsee.Geometry
=============

Extended methods for the ee.Geometry class:

.. currentmodule:: eemont.geometry

.. autosummary::
   :toctree: stubs

   BBoxFromQuery
   LinearRingFromPlusCodes
   LineStringFromPlusCodes
   MultiLineStringFromPlusCodes
   MultiPointFromPlusCodes
   MultiPointFromQuery
   MultiPolygonFromPlusCodes
   PointFromPlusCode
   PointFromQuery
   PolygonFromPlusCodes
   RectangleFromPlusCodes
   plusCodesee.App
===========

.. autoclass:: eemont.app.App
   :members:

   .. autoattribute:: ee.Image
==========

Extended methods for the ee.Image class:

.. currentmodule:: eemont.image

.. autosummary::
   :toctree: stubs
   
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   index
   maskClouds
   matchHistogram
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndices
   tasseledCapCommon
===========

Functions of the common module:

.. currentmodule:: eemont.common

.. autosummary::
   :toctree: stubs

   indices
   listIndicesAPI Reference
================

.. toctree::
   :caption: Earth Engine Module
   :maxdepth: 2
   :hidden:
   
   ee

.. toctree::
   :caption: Earth Engine Classes
   :maxdepth: 2
   :hidden:
   
   eedictionary
   eefeature
   eefeaturecollection
   eegeometry
   eeimage
   eeimagecollection
   eelist

.. toctree::
   :caption: New Classes
   :maxdepth: 2
   :hidden:
   
   eeapp
   
.. toctree::
   :caption: Non-Earth Engine Classes
   :maxdepth: 2
   :hidden:
   
   pddataframe
   
.. toctree::
   :caption: Additional Modules
   :maxdepth: 2
   :hidden:
   
   common

Extended Methods for the Earth Engine Module
--------------------------------------------

ee
~~

..
   .. currentmodule:: eemont.extra

   .. autosummary::

      install
      require
      uninstall

.. currentmodule:: eemont.app

.. autosummary::

   listApps

.. currentmodule:: eemont.common

.. autosummary::

   listDatasets

Extended Earth Engine Object Classes
------------------------------------------

Here you can find the reference of the new methods for each one of the Earth Engine classes:

ee.Dictionary
~~~~~~~~~~~~~

.. currentmodule:: eemont.eeDictionary

.. autosummary::

   __contains__
   __getitem__

ee.Feature
~~~~~~~~~~~~~

.. currentmodule:: eemont.feature

.. autosummary::

   BBoxFromQuery
   PointFromQuery
   plusCodes

ee.FeatureCollection
~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.featurecollection

.. autosummary::
   
   MultiPointFromQuery
   
ee.Geometry
~~~~~~~~~~~~~

.. currentmodule:: eemont.geometry

.. autosummary::

   BBoxFromQuery
   LinearRingFromPlusCodes
   LineStringFromPlusCodes
   MultiLineStringFromPlusCodes
   MultiPointFromPlusCodes
   MultiPointFromQuery
   MultiPolygonFromPlusCodes
   PointFromPlusCode
   PointFromQuery
   PolygonFromPlusCodes
   RectangleFromPlusCodes
   plusCodes

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::
   
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   index
   maskClouds
   matchHistogram
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndices
   tasseledCap
   
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection
.. autosummary::

   closest
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   getTimeSeriesByRegion
   getTimeSeriesByRegions
   index
   maskClouds
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndices
   tasseledCap

ee.List
~~~~~~~
   
.. currentmodule:: eemont.eeList

.. autosummary::

   __add__
   __contains__
   __getitem__
   __len__
   __mul__
   __radd__
   __rmul__

Extended Non-Earth Engine Object Classes
------------------------------------------

Non-Earth Engine classes such as pd.DataFrame are also extended:

pd.DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.dataframe
.. autosummary::

   toEEFeatureCollection
   
Additional Modules
------------------------------------------

Functions of additional modules:

Common
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.common
.. autosummary::

   indices
   listDatasets
   listIndicesee.Dictionary
=============

Extended methods for the ee.Dictionary class:

.. currentmodule:: eemont.eeDictionary

.. autosummary::
   :toctree: stubs

   __contains__
   __getitem__pd.DataFrame
==================

Extended methods for the pd.DataFrame class:

.. currentmodule:: eemont.dataframe

.. autosummary::
   :toctree: stubs

   toEEFeatureCollectionee.ImageCollection
====================

Extended methods for the ee.ImageCollection class:

.. currentmodule:: eemont.imagecollection

.. autosummary::
   :toctree: stubs

   closest
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   getTimeSeriesByRegion
   getTimeSeriesByRegions
   index
   maskClouds
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndices
   tasseledCapee.FeatureCollection
=======================

Extended methods for the ee.FeatureCollection class:

.. currentmodule:: eemont.featurecollection

.. autosummary::
   :toctree: stubs
   
   MultiPointFromQueryee.Feature
===========

Extended methods for the ee.Feature class:

.. currentmodule:: eemont.feature

.. autosummary::
   :toctree: stubs

   BBoxFromQuery
   PointFromQuery
   plusCodesee.ImageCollection
====================

Extended methods for the ee.ImageCollection class:

.. currentmodule:: eemont.imagecollection

.. autosummary::
   :toctree: stubs

   closest
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   getTimeSeriesByRegion
   getTimeSeriesByRegions
   index
   maskClouds
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndicespd.DataFrame
==================

Extended methods for the pd.DataFrame class:

.. currentmodule:: eemont.dataframe

.. autosummary::
   :toctree: stubs

   toEEFeatureCollectionAPI Reference
================

.. toctree::
   :caption: Earth Engine Classes
   :maxdepth: 2
   :hidden:
   
   eefeature
   eefeaturecollection
   eegeometry
   eeimage
   eeimagecollection
   
.. toctree::
   :caption: Non-Earth Engine Classes
   :maxdepth: 2
   :hidden:
   
   pddataframe
   
.. toctree::
   :caption: Additional Modules
   :maxdepth: 2
   :hidden:
   
   common

Extended Earth Engine Object Classes
------------------------------------------

Here you can find the reference of the new methods for each one of the Earth Engine classes:

ee.Feature
~~~~~~~~~~~~~

.. currentmodule:: eemont.feature

.. autosummary::

   BBoxFromQuery
   PointFromQuery
   plusCodes

ee.FeatureCollection
~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.featurecollection

.. autosummary::
   
   MultiPointFromQuery
   
ee.Geometry
~~~~~~~~~~~~~

.. currentmodule:: eemont.geometry

.. autosummary::

   BBoxFromQuery
   LinearRingFromPlusCodes
   LineStringFromPlusCodes
   MultiLineStringFromPlusCodes
   MultiPointFromPlusCodes
   MultiPointFromQuery
   MultiPolygonFromPlusCodes
   PointFromPlusCode
   PointFromQuery
   PolygonFromPlusCodes
   RectangleFromPlusCodes
   plusCodes

ee.Image
~~~~~~~~

.. currentmodule:: eemont.image

.. autosummary::
   
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   index
   maskClouds
   matchHistogram
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndices
   
ee.ImageCollection
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.imagecollection
.. autosummary::

   closest
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   getTimeSeriesByRegion
   getTimeSeriesByRegions
   index
   maskClouds
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndices

Extended Non-Earth Engine Object Classes
------------------------------------------

Non-Earth Engine classes such as pd.DataFrame are also extended:

pd.DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.dataframe
.. autosummary::

   toEEFeatureCollection
   
Additional Modules
------------------------------------------

Functions of additional modules:

Common
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: eemont.common
.. autosummary::

   indices
   listIndicesee.FeatureCollection
=======================

Extended methods for the ee.FeatureCollection class:

.. currentmodule:: eemont.featurecollection

.. autosummary::
   :toctree: stubs
   
   MultiPointFromQueryee.Image
==========

Extended methods for the ee.Image class:

.. currentmodule:: eemont.image

.. autosummary::
   :toctree: stubs
   
   getCitation
   getDOI
   getOffsetParams
   getScaleParams
   getSTAC
   index
   maskClouds
   matchHistogram
   panSharpen
   preprocess
   scale
   scaleAndOffset
   spectralIndicesCommon
===========

Functions of the common module:

.. currentmodule:: eemont.common

.. autosummary::
   :toctree: stubs

   indices
   listIndicesee.Geometry
=============

Extended methods for the ee.Geometry class:

.. currentmodule:: eemont.geometry

.. autosummary::
   :toctree: stubs

   BBoxFromQuery
   LinearRingFromPlusCodes
   LineStringFromPlusCodes
   MultiLineStringFromPlusCodes
   MultiPointFromPlusCodes
   MultiPointFromQuery
   MultiPolygonFromPlusCodes
   PointFromPlusCode
   PointFromQuery
   PolygonFromPlusCodes
   RectangleFromPlusCodes
   plusCodeseemont.imagecollection.spectralIndices
======================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: spectralIndiceseemont.geometry.BBoxFromQuery
=============================

.. currentmodule:: eemont.geometry

.. autofunction:: BBoxFromQueryeemont.feature.PointFromQuery
=============================

.. currentmodule:: eemont.feature

.. autofunction:: PointFromQueryeemont.image.scaleAndOffset
===========================

.. currentmodule:: eemont.image

.. autofunction:: scaleAndOffseteemont.image.preprocess
=======================

.. currentmodule:: eemont.image

.. autofunction:: preprocesseemont.geometry.PointFromQuery
==============================

.. currentmodule:: eemont.geometry

.. autofunction:: PointFromQueryeemont.geometry.MultiPointFromQuery
===================================

.. currentmodule:: eemont.geometry

.. autofunction:: MultiPointFromQueryeemont.dataframe.toEEFeatureCollection
======================================

.. currentmodule:: eemont.dataframe

.. autofunction:: toEEFeatureCollectioneemont.image.getScaleParams
===========================

.. currentmodule:: eemont.image

.. autofunction:: getScaleParamseemont.common.listIndices
=========================

.. currentmodule:: eemont.common

.. autofunction:: listIndiceseemont.imagecollection.closest
==============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: closesteemont.image.\_\_radd\_\_
=========================

.. currentmodule:: eemont.image

.. autofunction:: __radd__eemont.imagecollection.getOffsetParams
======================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getOffsetParamseemont.featurecollection.MultiPointFromQuery
============================================

.. currentmodule:: eemont.featurecollection

.. autofunction:: MultiPointFromQueryeemont.imagecollection.scaleAndOffset
=====================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: scaleAndOffseteemont.image.scale
==================

.. currentmodule:: eemont.image

.. autofunction:: scaleeemont.image.maskClouds
=======================

.. currentmodule:: eemont.image

.. autofunction:: maskCloudseemont.image.index
==================

.. currentmodule:: eemont.image

.. autofunction:: indexeemont.imagecollection.getTimeSeriesByRegions
=============================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getTimeSeriesByRegionseemont.imagecollection.getTimeSeriesByRegion
============================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getTimeSeriesByRegioneemont.imagecollection.preprocess
=================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: preprocesseemont.feature.BBoxFromQuery
============================

.. currentmodule:: eemont.feature

.. autofunction:: BBoxFromQueryeemont.image.getOffsetParams
============================

.. currentmodule:: eemont.image

.. autofunction:: getOffsetParamseemont.imagecollection.index
============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: indexeemont.image.\_\_add\_\_
========================

.. currentmodule:: eemont.image

.. autofunction:: __add__eemont.common.indices
=====================

.. currentmodule:: eemont.common

.. autofunction:: indiceseemont.imagecollection.scale
============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: scaleeemont.imagecollection.maskClouds
=================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: maskCloudseemont.image.spectralIndices
============================

.. currentmodule:: eemont.image

.. autofunction:: spectralIndiceseemont.imagecollection.getScaleParams
=====================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getScaleParamseemont.imagecollection.spectralIndices
======================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: spectralIndices﻿eemont.imagecollection.tasseledCap
==================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: tasseledCapeemont.imagecollection.getCitation
==================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getCitation﻿eemont.eeList.\_\_contains\_\_
==============================

.. currentmodule:: eemont.eeList

.. autofunction:: __contains__eemont.feature.plusCodes
========================

.. currentmodule:: eemont.feature

.. autofunction:: plusCodeseemont.geometry.BBoxFromQuery
=============================

.. currentmodule:: eemont.geometry

.. autofunction:: BBoxFromQueryeemont.feature.PointFromQuery
=============================

.. currentmodule:: eemont.feature

.. autofunction:: PointFromQueryeemont.imagecollection.getDOI
=============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getDOIeemont.image.scaleAndOffset
===========================

.. currentmodule:: eemont.image

.. autofunction:: scaleAndOffseteemont.image.getDOI
===================

.. currentmodule:: eemont.image

.. autofunction:: getDOIeemont.geometry.LineStringFromPlusCodes
=======================================

.. currentmodule:: eemont.geometry

.. autofunction:: LineStringFromPlusCodeseemont.image.preprocess
=======================

.. currentmodule:: eemont.image

.. autofunction:: preprocesseemont.geometry.LinearRingFromPlusCodes
=======================================

.. currentmodule:: eemont.geometry

.. autofunction:: LinearRingFromPlusCodeseemont.geometry.PointFromQuery
==============================

.. currentmodule:: eemont.geometry

.. autofunction:: PointFromQuery﻿eemont.image.tasseledCap
========================

.. currentmodule:: eemont.image

.. autofunction:: tasseledCap﻿eemont.image.panSharpen
=======================

.. currentmodule:: eemont.image

.. autofunction:: panSharpeneemont.geometry.MultiPointFromQuery
===================================

.. currentmodule:: eemont.geometry

.. autofunction:: MultiPointFromQueryeemont.dataframe.toEEFeatureCollection
======================================

.. currentmodule:: eemont.dataframe

.. autofunction:: toEEFeatureCollectioneemont.image.getScaleParams
===========================

.. currentmodule:: eemont.image

.. autofunction:: getScaleParams﻿eemont.eeList.\_\_rmul\_\_
==========================

.. currentmodule:: eemont.eeList

.. autofunction:: __rmul__﻿eemont.eeDictionary.\_\_getitem\_\_
===================================

.. currentmodule:: eemont.eeDictionary

.. autofunction:: __getitem__eemont.common.listIndices
=========================

.. currentmodule:: eemont.common

.. autofunction:: listIndices﻿eemont.eeDictionary.\_\_contains\_\_
====================================

.. currentmodule:: eemont.eeDictionary

.. autofunction:: __contains__﻿eemont.extra.install
====================

.. currentmodule:: eemont.extra

.. autofunction:: installeemont.imagecollection.closest
==============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: closest﻿eemont.eeList.\_\_add\_\_
=========================

.. currentmodule:: eemont.eeList

.. autofunction:: __add__eemont.image.\_\_radd\_\_
=========================

.. currentmodule:: eemont.image

.. autofunction:: __radd__eemont.geometry.plusCodes
=========================

.. currentmodule:: eemont.geometry

.. autofunction:: plusCodeseemont.geometry.PolygonFromPlusCodes
====================================

.. currentmodule:: eemont.geometry

.. autofunction:: PolygonFromPlusCodeseemont.imagecollection.getOffsetParams
======================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getOffsetParams﻿eemont.eeList.\_\_getitem\_\_
=============================

.. currentmodule:: eemont.eeList

.. autofunction:: __getitem__eemont.featurecollection.MultiPointFromQuery
============================================

.. currentmodule:: eemont.featurecollection

.. autofunction:: MultiPointFromQueryeemont.geometry.MultiLineStringFromPlusCodes
============================================

.. currentmodule:: eemont.geometry

.. autofunction:: MultiLineStringFromPlusCodeseemont.imagecollection.scaleAndOffset
=====================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: scaleAndOffseteemont.image.scale
==================

.. currentmodule:: eemont.image

.. autofunction:: scaleeemont.image.maskClouds
=======================

.. currentmodule:: eemont.image

.. autofunction:: maskCloudseemont.geometry.RectangleFromPlusCodes
======================================

.. currentmodule:: eemont.geometry

.. autofunction:: RectangleFromPlusCodeseemont.image.getSTAC
====================

.. currentmodule:: eemont.image

.. autofunction:: getSTACeemont.image.index
==================

.. currentmodule:: eemont.image

.. autofunction:: indexeemont.imagecollection.getTimeSeriesByRegions
=============================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getTimeSeriesByRegions﻿eemont.eeList.\_\_len\_\_
=========================

.. currentmodule:: eemont.eeList

.. autofunction:: __len__﻿eemont.extra.uninstall
======================

.. currentmodule:: eemont.extra

.. autofunction:: uninstalleemont.imagecollection.getTimeSeriesByRegion
============================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getTimeSeriesByRegioneemont.imagecollection.getSTAC
==============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getSTAC﻿eemont.eeList.\_\_radd\_\_
==========================

.. currentmodule:: eemont.eeList

.. autofunction:: __radd__eemont.imagecollection.preprocess
=================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: preprocess﻿eemont.extra.require
====================

.. currentmodule:: eemont.extra

.. autofunction:: requireeemont.geometry.MultiPolygonFromPlusCodes
=========================================

.. currentmodule:: eemont.geometry

.. autofunction:: MultiPolygonFromPlusCodeseemont.feature.BBoxFromQuery
============================

.. currentmodule:: eemont.feature

.. autofunction:: BBoxFromQuery﻿eemont.image.matchHistogram
===========================

.. currentmodule:: eemont.image

.. autofunction:: matchHistogrameemont.image.getOffsetParams
============================

.. currentmodule:: eemont.image

.. autofunction:: getOffsetParamseemont.geometry.MultiPointFromPlusCodes
=======================================

.. currentmodule:: eemont.geometry

.. autofunction:: MultiPointFromPlusCodeseemont.imagecollection.index
============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: indexeemont.image.\_\_add\_\_
========================

.. currentmodule:: eemont.image

.. autofunction:: __add__eemont.common.indices
=====================

.. currentmodule:: eemont.common

.. autofunction:: indiceseemont.image.getCitation
========================

.. currentmodule:: eemont.image

.. autofunction:: getCitationeemont.imagecollection.scale
============================

.. currentmodule:: eemont.imagecollection

.. autofunction:: scaleeemont.geometry.PointFromPlusCode
=================================

.. currentmodule:: eemont.geometry

.. autofunction:: PointFromPlusCode﻿eemont.eeList.\_\_mul\_\_
=========================

.. currentmodule:: eemont.eeList

.. autofunction:: __mul__﻿eemont.imagecollection.panSharpen
=================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: panSharpeneemont.imagecollection.maskClouds
=================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: maskCloudseemont.image.spectralIndices
============================

.. currentmodule:: eemont.image

.. autofunction:: spectralIndiceseemont.imagecollection.getScaleParams
=====================================

.. currentmodule:: eemont.imagecollection

.. autofunction:: getScaleParamseemont.common.indices
=====================

.. currentmodule:: eemont.common

.. autofunction:: indiceseemont.dataframe.toEEFeatureCollection
======================================

.. currentmodule:: eemont.dataframe

.. autofunction:: toEEFeatureCollectioncommon
======

Functions of the common module:

.. automodule:: eemont.common
   :members:
   :undoc-members:
   :show-inheritance:Other Modules
===============

Additional eemont modules are listed below:

.. toctree::
   :caption: Other Modules
   :maxdepth: 2

   commonOther Modules
===============

Additional eemont modules are listed below:

.. toctree::
   :caption: Other Modules
   :maxdepth: 2

   commoncommon
======

Functions of the common module:

.. automodule:: eemont.common
   :members:
   :undoc-members:
   :show-inheritance: