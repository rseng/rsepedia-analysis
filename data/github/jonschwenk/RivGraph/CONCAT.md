---
title: 'RivGraph: Automatic extraction and analysis of river and delta channel network topology'
tags:
  - Python
  - rivers
  - deltas
  - image processing
  - networks
  - channel network extraction
  - fluvial geomorphology
authors:
  - name: Jon Schwenk^[Corresponding author]
    orcid: 0000-0001-5803-9686
    affiliation: "1"
  - name: Jayaram Hariharan
    orcid: 0000-0002-1343-193X
    affiliation: "2"
affiliations:
 - name: Los Alamos National Laboratory, Division of Earth and Environmental Sciences
   index: 1
 - name: Department of Civil, Architectural and Environmental Engineering, The University of Texas at Austin
   index: 2
date: 01 January 2021
bibliography: paper.bib
---

# Summary
River networks sustain life and landscapes by carrying and distributing water, sediment, and nutrients throughout ecosystems and communities. At the largest scale, river networks drain continents through tree-like *tributary* networks. At typically smaller scales, river deltas and braided rivers form loopy, complex *distributary* river networks via avulsions and bifurcations. In order to model flows through these networks or analyze network structure, the topology, or connectivity, of the network must be resolved. Additionally, morphologic properties of each river channel as well as the direction of flow through the channel inform how fluxes travel through the network's channels.

`RivGraph` is a Python package that automates the extraction and characterization of river channel networks from a user-provided binary image, or mask, of a channel network (Fig. 1). Masks may be derived from (typically remotely-sensed) imagery, simulations, or even hand-drawn. `RivGraph` will create explicit representations of the channel network by resolving river centerlines as links, and junctions as nodes. Flow directions are solved for each link of the network without using auxiliary data, e.g., a digital elevation model (DEM). Morphologic properties are computed as well, including link lengths, widths, sinuosities, branching angles, and braiding indices. If provided, `RivGraph` will preserve georeferencing information of the mask and will export results as ESRI shapefiles, GeoJSONs, and GeoTIFFs for easy import into GIS software. `RivGraph` can also return extracted networks as `networkx` objects for convenient interfacing with the full-featured `networkx` package [@hagberg2008]. Finally, `RivGraph` offers a suite of topologic metrics that were specifically designed for river channel network analysis [@tejedor2015].

![The core functionality of RivGraph for a delta channel network.](./examples/images/rivgraph_overview_white.PNG)

# Statement of need

Satellite and aerial photography have provided unprecedented opportunities to study the structure and dynamics of rivers and their networks. As both the quantity and quality of these remotely-sensed observations grow, the need for tools that automatically map and measure river channel network properties has grown in turn. The genesis of `RivGraph` is rooted in the work of [@tejedor2015, @tejedor2015a, and @tejedor2017] in a revitalized effort to see river channel networks through the lenses of their network structure. The authors were relegated to time-consuming hand-delineations of the delta channel networks they analyzed.  `RivGraph` was thus born from a need to transform binary masks of river channel networks into their graphical representations accurately, objectively, and efficiently.

`RivGraph` has already been instrumental in a number of investigations. The development of the flow directions algorithms itself provided insights into the nature of river channel network structure in braided rivers and deltas [@schwenk2020]. For deltas specifically, `RivGraph`-extracted networks have been used to study how water and sediment are partitioned at bifurcations [@dong2020], to determine how distance to the channel network plays a controlling role on Arctic delta lake dynamics [@vulis2020], and to construct a network-based model of nitrate removal across the Wax Lake Delta [@knights2020]. For braided rivers, `RivGraph` was used to extract channel networks from hydrodynamic simulations in order to develop the novel "entropic braiding index" [eBI, @tejedor2019], and a function for computing the eBI (as well as the classic braiding index) for braided rivers is provided in `RivGraph`. The work of @marra2014 represented an effort to understand braided rivers through their topologies, although their networks were apparently extracted manually. Ongoing, yet-unpublished work is using `RivGraph` to study river dynamics, delta loopiness, and nutrient transport through Arctic deltas.

We are aware of one other package that extracts network topology from channel network masks. The  `Orinoco` Python package [@marshak2020] uses a fast marching method to resolve the channel network in contrast to `RivGraph`'s skeletonization approach. `Orinoco` uses only a shortest-path approach for setting flow directions rather than `RivGraph`'s exploitation of many morphologic features (including shortest path) to set flow directions. If a DEM of the channel network is available, the Lowpath [@hiatt2020] add-on to the [Topological Tools for Geomorphological Analysis](https://github.com/tue-alga/ttga) package may be of interest. `RivGraph`'s along-river mesh generation for braided rivers (Fig. 2) was inspired by RivMAP [@schwenk2017b].

# Functionality

`RivGraph` requires the user to provide a binary mask of a channel network. If the provided mask is georeferenced (e.g., a GeoTIFF), `RivGraph` will export results in the same coordinate reference system (CRS) for easy analysis with a Geographical Information System (GIS) software. Otherwise, a "dummy" CRS is applied, and calculated physical quantities (e.g., length and width) will be in units of pixels.  The channel mask is the basis for all `RivGraph` processing, so the user should consider carefully the features to include and the desired level of smoothing. `RivGraph` respects the connectivity of the channel mask such that all groups of pixels connected in the mask will be connected in the vectorized representation as well. The user may therefore wish to preprocess their mask to fill small or unwanted islands or smooth channel boundaries. `RivGraph`'s `im_utils()` module contains a number of functions for achieving these tasks, including island-filling and morphological operators. Detailed information about how to create and prepare masks is provided in [the documentation](https://jonschwenk.github.io/RivGraph/maskmaking/index.html).

### Basic Functionality

`RivGraph` was designed with an emphasis on user-friendliness and accessibility, guided by the idea that even novice Python users should be able to make use of its functionality. Anticipated common workflows are gathered into two classes that manage georeferencing conversions, path management, and I/O with simple, clearly-named methods. Beginning users will want to instantiate either a `delta` or a (braided) `river` class and apply the relevant methods, which are as follows:

- `skeletonize()` : skeletonizes the mask; minor conditioning of the skeleton is performed to simplify the topology. For example, if a "+" pattern with the center pixel "off" appears in the skeleton, the center pixel will be added to the skeleton to reduce the number of branchpoints from four to one. 
- `compute_network()` : walks along the skeleton to resolve the links and nodes. All pixels in the unpruned skeleton will be visited and therefore represented in the vectorized output.
- `prune_network()` : removes portions of the network that do not contribute meaningfully to its topology. The skeletonization process often results in many "dangling links," or links connected to the network at only one end. During pruning, all dangling links are removed except those connected to inlet or outlet nodes. For the `delta` class, user-provided shoreline and inlet nodes files are required so that `RivGraph` can prune the network to the shoreline and identify the inlet and outlet nodes. Additionally, bridging links, or links whose removal results in two subnetworks, are removed if one of the resulting subnetworks contains no inlet or outlet nodes. The corresponding subnetwork without inlet or outlet nodes is also removed.
- `compute_link_width_and_length()` : adds width and length attributes to each link. Width is computed via sampling a distance transform image along the link (centerline) coordinates and multiplying by two. Length is the sum of the Euclidean distance between each pair of pixels along a link.
- `assign_flow_directions()` : uses a suite of algorithms to set the flow direction of each link in the network. Separate "recipes" are provided for the `delta` and `river` classes, but users may also create their own. Rationale of the various algorithms and details of recipe construction are given in @schwenk2020.

Additional methods are available for plotting, exporting GeoTIFFs and geovectors, saving/loading the network, converting to adjacency matrices, computing junction angles, and finding islands.

Braided rivers should be analyzed with the `river` class, which instead of a user-provided shoreline requires a two-character string denoting the *exit sides* of the river with respect to the mask, e.g., 'NS' for a river whose upstream terminus is at the top of the image and downstream at the bottom. `RivGraph` exploits the general direction of the braided river's channel belt to set flow directions and generate an along-river mesh (Fig. 2) that can be used for characterizing downstream changes. In addition to the methods above, the `river` class also features:

- `compute_centerline()` : computes the centerline of the holes-filled river mask (not individual channels)
- `compute_mesh()` : creates a mesh of evenly-spaced transects that are approximately perpendicular to the centerline. The user can specify the mesh spacing, transect width, and degree of smoothing.

![Figure 2. A RivGraph-generated mesh for a mask of the Indus River.](./examples/images/indus_mesh_paper.JPG){width=40%}

### Advanced Functionality

`RivGraph` is organized into a set of modules such that users can find particular functions based on their general class. Customized workflows can be created by calling appropriate functions from these modules, which include

- `classes` : contains the `river` and `delta` classes and associated methods
- `directionality` : algorithms for setting flow directions that are not specific to deltas or braided rivers
- `geo_utils` : functions for handling geospatial data
- `im_utils` : image processing utilities, including morphologic operators
- `io_utils` : functions for reading and writing data and results
- `ln_utils` : functions for building and manipulating the links and nodes of the network
- `mask_to_graph` : the algorithm for converting the mask to a set of links and nodes
- `walk` : functions for walking along the skeleton and identifying branchpoints
- `deltas/delta_directionality` : delta-specific algorithms for setting flow directions
- `deltas/delta_metrics` : functions for computing topologic metrics
- `deltas/delta_utils` : algorithm for pruning deltas and clipping the delta network by the shoreline
- `rivers/river_directionality` : river-specific algorithms for setting flow directions
- `rivers/river_utils` : algorithms for pruning rivers and generating along-river meshes

# Dependencies

`RivGraph` relies on functionality from the following Python packages: GDAL [@gdal2020], NumPy [@harris2020], Matplotlib [@hunter2007], GeoPandas [@jordahl2020], Shapely [@gillies2007], Fiona [@gillies2011], pyproj [@snow2020], scikit-image [@vanderwalt2014], OpenCV [@bradski2000], networkx [@hagberg2008], and fastdtw [@slaypni2020].

# Acknowledgements

We thank Efi Foufoula-Georgiou, Alejandro Tejedor, Anthony Longjas, Lawrence Vulius, Kensuke Naito, and Deon Knights for providing test cases and feedback for `RivGraph`'s development. We are also grateful to Anastasia Piliouras and Joel Rowland for providing valuable insights and subsequent testing of `RivGraph`'s flow directionality algorithms.

`RivGraph` has received financial support from NSF under EAR-1719670, the United States Department of Energy, and Los Alamos National Laboratory's Lab Directed Research and Development (LDRD) program. Special thanks are due to Dr. Efi Foufoula-Georgiou for providing support during the nascent phase of `RivGraph`'s development.

# References

[![build](https://github.com/jonschwenk/RivGraph/actions/workflows/build.yml/badge.svg)](https://github.com/jonschwenk/RivGraph/actions/workflows/build.yml)
[![Coverage Status](https://coveralls.io/repos/github/jonschwenk/RivGraph/badge.svg?branch=master)](https://coveralls.io/github/jonschwenk/RivGraph?branch=master)
![docs](https://github.com/jonschwenk/RivGraph/workflows/docs/badge.svg)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02952/status.svg)](https://doi.org/10.21105/joss.02952)
<br />

[![RivGraph logo](https://github.com/jonschwenk/RivGraph/blob/master/docs/logos/rg_logo_full.png)](https://jonschwenk.github.io/RivGraph/ "Go to documentation.")

About
-----

RivGraph is a Python package that provides tools for converting a binary mask of a channel network into a directed, weighted graph (i.e. a set of connected links and nodes).

![Core functionality of RivGraph.\label{fig:corefunctions}](https://github.com/jonschwenk/RivGraph/blob/master/examples/images/rivgraph_overview_white.PNG)

The figure above demonstrates the core components of RivGraph, but many other features are provided, including:

- Morphologic metrics (lengths, widths, branching angles, braiding indices)
- Algebraic representations of the channel network graph
- Topologic metrics (both topologic and dynamic such as alternative paths, flux sharing, entropies, mutual information, etc.)
- Tools for cleaning and preparing your binary channel network mask
- Island detection, metrics, and filtering
- Mesh generation for characterizing along-river characteristics
- (beta) Tools for centerline migration analysis

All of RivGraph's functionality maintains and respects georeferencing information. If you start with a georeferenced mask (e.g. a GeoTIFF), RivGraph exports your results in the CRS (coordinate reference system) of your mask for convenient mapping, analysis, and fusion with other datasets in a GIS.

You can see some description of RivGraph's functionality via this [AGU poster](https://www.researchgate.net/publication/329845073_Automatic_Extraction_of_Channel_Network_Topology_RivGraph), and the flow directionality logic and validation is described in our [ESurf Dynamics paper](https://www.earth-surf-dynam.net/8/87/2020/esurf-8-87-2020.html). Examples demonstrating the basic RivGraph features are available for a [delta channel network](https://github.com/jonschwenk/RivGraph/blob/master/examples/delta_example.ipynb) and a [braided river](https://github.com/jonschwenk/RivGraph/blob/master/examples/braided_river_example.ipynb).

Installing
-----
RivGraph v0.4 is hosted on the anaconda channel [jschwenk](https://anaconda.org/jschwenk/rivgraph). We recommend installing into a fresh conda environment to minimize the risk of dependency clashes. The easiest way to do this is to download the [environment.yml](https://github.com/jonschwenk/RivGraph/blob/master/environment.yml) file, then open Terminal (Mac/Unix) or Anaconda Prompt (Windows) and type:

<pre><code>conda env create --file /path/to/environment.yml  # the environment name will be 'rivgraph', but you can change the environment file to name it anything</code></pre>

You may then want to install Spyder or your preferred IDE. Conda should fetch all the required dependencies and handle versioning.

If you want to install RivGraph into an already-existing environment, you can run <pre><code>conda activate myenv
conda install rivgraph -c jschwenk</code></pre>

You may also [install RivGraph from this Github repo](https://jonschwenk.github.io/RivGraph/install/index.html#installation-from-source).

Instructions for testing your installation are available [here](https://jonschwenk.github.io/RivGraph/install/index.html#installation-from-source).

How to use?
-----
Please see the [documentation](https://jonschwenk.github.io/RivGraph/) for more detailed instructions.

RivGraph requires that you provide a binary mask of your network. [This page](https://jonschwenk.github.io/RivGraph/maskmaking/index.html) provides some help, hints, and tools for finding or creating your mask.

To see what RivGraph does and how to operate it, you can work through the [Colville Delta example](https://github.com/jonschwenk/RivGraph/blob/master/examples/delta_example.ipynb) or the [Brahmaputra River example](https://github.com/jonschwenk/RivGraph/blob/master/examples/braided_river_example.ipynb). Both examples include sample masks.

RivGraph contains two primary classes (`delta` and `river`) that provide convenient methods for creating a processing workflow for a channel network. As the examples demonstrate, you can instantiate a delta or river class, then apply associated methods for each. After looking at the examples, take a look at [classes.py](https://github.com/jonschwenk/RivGraph/blob/master/rivgraph/classes.py) to understand what methods are available.

**Note**: there are many functions under the hood that may be useful to you. Check out the [im_utils script](https://github.com/jonschwenk/RivGraph/blob/master/rivgraph/im_utils.py) (image utilities) in particular for functions to help whip your mask into shape!


Contributing
------------
If you think you're not skilled or experienced enough to contribute, think again! We agree wholeheartedly with the sentiments expressed by this [Imposter syndrome disclaimer](https://github.com/Unidata/MetPy#contributing). We welcome all forms of user contributions including feature requests, bug reports, code, documentation requests, and code. Simply open an issue in the [tracker](https://github.com/jonschwenk/RivGraph/issues). For code development contributions, please contact us via email to be added to our slack channel where we can hash out a plan for your contribution.

Citing RivGraph
------------

Citations help us justify the effort that goes into building and maintaining this project. If you used RivGraph for your research, please consider citing us.

If you use RivGraph's flow directionality algorithms, please cite our [ESurf Dynamics paper](https://www.earth-surf-dynam.net/8/87/2020/esurf-8-87-2020.html). Additionally, if you publish work wherein RivGraph was used to process your data, please cite our [JOSS Paper](https://joss.theoj.org/papers/10.21105/joss.02952).

Contacting us
-------------

The best way to get in touch is to [open an issue](https://github.com/jonschwenk/rivgraph/issues/new) or comment on any open issue or pull request. Otherwise, send an email to j.........k@gmail.com


License
------------

This is free software: you can redistribute it and/or modify it under the terms of the **BSD 3-clause License**. A copy of this license is provided in [LICENSE.txt](https://github.com/jonschwenk/RivGraph/blob/master/LICENSE.txt).

RivGraph has been assigned number C19049 by the Feynman Center for Innovation.
Welcome to RivGraph's documentation!
====================================

RivGraph provides two classes that you can instantiate with your binary channel mask: *delta* and (braided) *river*. These classes contain methods that may be applied to execute network extractions, perform analyses, and write (georeferenced) results. While most of RivGraph's processing is automated, deltaic channel networks require the user to also create and provide a shoreline and inlet nodes shapefiles. Take a look at the `FAQ  <todo>`_for quick answers.

If the documentation you seek is not herein, please open a request via the `Issue Tracker <https://github.com/jonschwenk/RivGraph/issues>`_. `Detailed examples <https://github.com/jonschwenk/RivGraph/tree/master/examples>`_ for use are available, and most functions have been documented following a modified `numpy docstring format <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

Documentation
-------------

.. toctree::
   :maxdepth: 1

   quickstart/index
   background/index
   install/index
   examples/index
   maskmaking/index
   shoreline/index
   linksnodes/index
   issues/index
   featuredevelopment/index
   contributing/index
   gallery/index
   apiref/index
.. _issues:

************
Known issues 
************

There are several known issues that new, and experienced, users may encounter. These may arise when some of the assumptions of RivGraph are violated, and this will serve as a resource to troubleshoot what is wrong. Many of these isues have been identified in the `issues page <https://github.com/jonschwenk/RivGraph/issues>`_ and found during regular use of RivGraph. If you encounter more issues, please add them to the issues page.

Shoreline clipping issues
=========================
See the `shoreline documentation <https://jonschwenk.github.io/RivGraph/shoreline/index.html>`_ further information on methods to generate a shoreline. 

Left id failure
---------------

Shoreline features should not contain any attributes named left_fid, LEFT_ID, etc. which popular GIS software such as ArcGIS & QGIS can add. It is recommended to enter a shoreline with a single attribute, id, with a null value inside to avoid an `error <https://github.com/jonschwenk/RivGraph/issues/9>`_.

Non-point intersection
----------------------

An assumption of `clip_by_shoreline` is that the channel skeleton and the shoreline intersect at a point, not at a line segment. If a shoreline was extracted using an automatic method such as the `Opening Angle Method <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2008GL033963>`_, a line intersection can arise. It's recommended to smooth the shoreline, e.g. using `spline interpolation <https://gis.stackexchange.com/questions/24827/smoothing-polygons-in-contour-map>`_ or manually adjust the shoreline in that case to ensure point intersection. 

Mask issues
===========

The basics of maskmaking and mask assumptions are covered in the `maskmaking documentation <https://jonschwenk.github.io/RivGraph/maskmaking/index.html>`_, and this section details what happens when something goes awry. 

Skeleton around the image forming a cyclic graph
------------------------------------------------

As noted in the `maskmaking docs <https://jonschwenk.github.io/RivGraph/maskmaking/index.html>`_, there should be no no-data pixels in the mask. This is because no-data pixels are treated as 1 (i.e. water), and `this can result <https://github.com/jonschwenk/RivGraph/issues/34>`_ in a cyclical network which loops back to the delta apex. Here's an example of the resulting skeleton where the teal is no data, blue is water, grey is land, yellow is the shoreline, and orange is the channel network skeleton after pruning. This issue can arise when projecting a mask from WGS84 (EPSG:4326) imagery to a local UTM zone, which results in the generation of triangular bands of no-data on the edge of the image. The NA values in the projected mask should be set as land (0) prior to running channel network skeletonization to prevent any issues. 

.. image:: https://user-images.githubusercontent.com/18738680/103107918-c4725d00-45f7-11eb-990f-c6b49bebeba9.png

Multiple inlet nodes after assigning directionality
---------------------------------------------------

After successfully extracting the channel network and assigning flow directionality, an `error <https://github.com/jonschwenk/RivGraph/issues/52>`_ that may arise when computing the delta metrics is an error that there are multiple inlet nodes found. This arises from an assumption built into RivGraph as of v0.4: the skeletonization of the mask will fill any islands 4 pixels or less, in a 4-neighbor sense, in the network. If these islands exist in the mask the skeleton can intersect one of these islands, which may in some instances cause an error. Here is an example of where this occurs and an issue arises. In the image below the link intersecting the island will have an `wid_adj` of zero, as the adjusted vector `wid_pix` will have a value of essentially zero. 

.. image:: https://user-images.githubusercontent.com/14874485/118341935-709cdd80-b4de-11eb-87a6-bbbc24fc8aec.png

Without this hole-filling assumption the following skeleton would be generated:

.. image:: https://user-images.githubusercontent.com/14874485/118342093-236d3b80-b4df-11eb-81b1-ec032d841e0e.png

Depending on the analysis you are performing, it may or may not be appropriate to have the skeleton shown in the second case. This issue can arise in both the ``delta`` and the ``river`` skeletonization methods.

As changing this assumption may lead to unknown downstream changes in RivGraph's functionality, there are no plans at this time to remove it from the skeletonization procedure. 

Therefore there are several options available to treat this issue if it arises. If the island is relevant to the problem at hand, i.e. represents a significant feature in the network being analyzed: 

1) Comment out the following three lines in ``skeletonize_mask`` for your relevant class:

``Iskel = imu.fill_holes(Iskel, maxholesize=4)``

``Iskel = morphology.skeletonize(Iskel)``

``Iskel = simplify_skel(Iskel)``

2. Downscale the mask such that the island becomes larger than 4 pixels. This will increase processing time, but for relatively small masks may not be significant. 

If the island doesn't represent a significant feature in the network and could be removed: 

3. Remove islands than 4 pixels during mask processing. The function `im_utils:fill_holes` provides this functionality.
.. _maskmaking:

==========
Maskmaking
==========

.. image:: ../../images/lena_mask.PNG
  :alt: The Lena River Delta
  :align: center
.. centered::
  The Lena River Delta as seen by Bing Virtual Earth and its mask. Note that this mask is Landsat-derived so it does not correspond perfectly with the image on the left.

-----------------------------------
Maskmaker, Maskmaker Make Me a Mask
-----------------------------------

*RivGraph* requires that you provide a mask of your channel network. In this document, we'll cover the following:


 - :ref:`whatismask`
 - :ref:`maskconstraints`
 - :ref:`maskcapture`
 - :ref:`wheretoget`
 - :ref:`howtoprep`
 - :ref:`georef`
 - :ref:`nonriver`
 - :ref:`supportedfiletypes`

.. important::
  One thing to keep in mind: although *RivGraph* contains functions for pruning and otherwise modifying your channel network, it will always honor the mask you provide. You may need to iterate between altering your mask and *RivGraph*-ing it to achieve your desired results.
  "Garbage in, garbage out."

.. _whatismask:

---------------
What is a mask?
---------------
A mask is simply a binary image (only ones and zeros) where pixels belonging to the channel network are ones, like the right panel of the Lena Delta above. Before processing your mask with *RivGraph*, you should ensure that your mask contains `no no-data <https://github.com/jonschwenk/RivGraph/issues/34>`_. One way to ensure this is to convert your mask to a boolean datatype, using for example numpy:

:code:`Mask_binary = np.array(Mask, dtype=np.bool)`

The mask is the cornerstone for using *RivGraph*. You should always ensure that it contains the features you want and none of the ones you don't.

.. tip:: Make sure to remove all the objects (connected "on" pixels) in your mask that you do not want to analyze. Leaving in other objects can cause unexpected behavior or `errors <https://github.com/jonschwenk/RivGraph/issues/32>`_. Often, a quick way to achieve this is via the :obj:`rivgraph.im_utils.largest_blobs()` function, which will keep only the largest connected component.


.. _maskconstraints:

-----------------------------------------------
What does RivGraph expect my mask to look like?
-----------------------------------------------
*RivGraph* can handle two types of masks: deltas and braided (or single-threaded) rivers. 

**All masks should contain a single connected component** (or blob). Before doing your analysis, use the :obj:`rivgraph.im_utils.largest_blobs()` to ensure you have only one connected component. This will remove any isolated portions of the mask. For example,

.. code-block:: python3

   from rivgraph import im_utils
   Mask_one_blob = im_utils.largest_blobs(Mask, nlargest=1, action='keep')

This will create an image wherein only the largest connected component remains (hopefully your channel network). 

Delta masks often have a large waterbody (ocean or lake) connecting all the outlets. This is completely fine, as the waterbody will be removed in the pruning stage. 


.. _maskcapture:

---------------------------
What should a mask capture?
---------------------------

Above, we defined a mask as all the pixels "belonging to the channel network." But which pixels are part of the channel network? The obvious starting point is to consider all *surface water* pixels as defining the channel network. Is it important that your mask shows bankfull channels or includes small streams that may only be active during flood conditions?


.. image:: ../../images/brahma_masks_2004.PNG
 :align: center

As an example, four masks are shown above of a portion of the Brahmaputra River at four different days of 2004. The hydrograph of the river is shown in the left panel. As expected, you can see that the mask changes as the river floods and recedes. You can also see that the river channel network has been rearranged in some places by the flooding. For example, **A** and **D** show the river at roughly the same discharge, but the network in **D** is a result of reworking by the flood.

The bottom line here is that your mask should reflect your analysis goals. For example, if you're only interested in counting the number of loops in a delta channel network, then you might only care about ensuring all the channels are represented. If you want to route fluxes through your extracted network, then you probably should try to obtain a mask that captures the full width of the channels at some representative discharge.

.. _wheretoget:

----------------------
Where do I get a mask?
----------------------
Masks can come from a variety of sources, but in my experience there are three primary methods for mask generation:

  - automatically generated from satellite imagery
  - manually drawn by hand
  - model/simulation outputs

There are *many* methods available for creating masks automatically from remotely-sensed imagery. We won't get into the details of those here, but note that machine learning has proved a very valuable tool for maskmaking. There are also simple, proven techniques available as well. The Brahmaputra masks above were created by thresholding the Landsat-derived NDVI (`Normalized Difference Vegetation Index <https://www.usgs.gov/core-science-systems/nli/landsat/landsat-normalized-difference-vegetation-index>`_
), which is a simple ratio of band values.

Drawing a mask by hand is often not an ideal choice, but might be the most efficient way to move forward. In these cases, I would typically use QGIS to draw polygons that cover the channel network, then use the `Rasterize  <https://docs.qgis.org/2.8/en/docs/user_manual/processing_algs/gdalogr/gdal_conversion/rasterize.html>`_
tool to convert the polygons to a binary raster (image). If you go this route, be sure to specify an appropriate coordinate reference system for your polygons in order to preserve the georeferencing information (don't use EPSG:4326). You will also need to specify a pixel resolution for your mask upon conversion.

If you're analyzing the output of a simulation, it is unlikely that the simulation will provide binary channel masks as an output. In these cases, you will need to develop a way to identify the channel network from the available simulation results. For example, while developing the entropic Braided Index (`eBI <https://ui.adsabs.harvard.edu/abs/2019AGUFMEP51E2163T/abstract>`_
), we used Delft3D simulations to test hypotheses about how the eBI changes under various sedimentation schemes. To make masks, we developed a combined depth + discharge threshold to identify which pixels were part of the "active river channel."

Here are some resources that either provide masks or tools for you to make your own.

- Published masks:

  - `Arctic deltas <https://data.ess-dive.lbl.gov/view/doi:10.15485/1505624>`_, made with eCognition and Landsat imagery.
  - `Indus and Brahmaputra Rivers <https://esurf.copernicus.org/articles/8/87/2020/#section6>`_, clipped from GRWL dataset.
  - `Global mask <https://zenodo.org/record/1297434>`_ of Landsat-derived rivers at "mean annual discharge." Has some issues at tile boundaries, and can be "feathery" along braided rivers, but not a bad global mask.
  - `Global Surface Water Dataset <https://global-surface-water.appspot.com/>`_ - provides all water pixels in the Landsat archive as monthly global images and as integrated-through-time images. For example, can threshold on the "Occurrence" product to make a mask. Use `Google Earth Engine <https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_2_GlobalSurfaceWater>`_ to access and create your masks.
  - If you know of more, please mention them in the `Issue Tracker <https://github.com/jonschwenk/RivGraph/issues>`_!

.. image:: ../../images/jrc_mackenzie.PNG
 :align: center

.. centered::
  The Global Surface Water's *Occurrence* map shows the fraction of time an observable Landsat pixel was water.


- You can relatively quickly train and apply ML models using `Google Earth Engine <https://earthengine.google.com/>`_, although the learning curve may be a little steep if you haven't used it before.

- `DeepWaterMap  <https://github.com/isikdogan/deepwatermap>`_ is a trained deep convolutional neural network that you can apply to Landsat/Sentinel multispectral imagery to create your own masks. You can also improve DeepWaterMap's base model by adding more training data. Requires some knowledge of Tensorflow.




.. _howtoprep:

-----------------------------
How do I edit my mask?
-----------------------------

As a mask is simply a single-band image, any pixel-based image editing software can be used for hand-editing (Photoshop, GIMP, MSPaint, etc.). However, there are a few issues with using these tools:

- These softwares will generally not preserve georeferencing information of your source image. You will have to add it back to the edited image.
- The softwares may have difficulty opening/editing a single-band image as opposed to the more standard RGB (3 band).
- Filetypes are sometimes not compatible between Python-exported images and these softwares and will thus require extra attention.

I have found three effective ways to edit georeferenced masks. The one you choose depends on the quantity and quality of editing you need to achieve.

1) Edit your mask directly in QGIS.

   a) `Serval  <https://plugins.qgis.org/plugins/Serval/>`_ plugin for QGIS allows for single-pixel manipulations. Good if you only need to edit a handful of pixels.

   b) `ThRaSe  <https://plugins.qgis.org/plugins/ThRasE/>`_ plugin for QGIS appears to have more sophisticated raster-editing capabilities, but I haven't used it.

2)  `Paint.NET <https://www.getpaint.net/download.html>`_ is an image-editing software that preserves georeferencing information. It's fairly basic and easy to use. If you have a significant amount of hand-editing to do, look into it.

3) Use image processing tools in *RivGraph* to edit your mask. There are morphological operators like :obj:`rivgraph.im_utils.dilate()` and :obj:`rivgraph.im_utils.erode()`, :obj:`rivgraph.im_utils.regionprops()` for filtering objects based on their properties (areas, lengths, perimeters, etc.), and :obj:`rivgraph.im_utils.largest_blobs()` for keeping/removing the largest connected components in the mask. There is also a :obj:`rivgraph.im_utils.hand_clean()` utility that allows you to draw polygons one-at-a-time and specify their pixel values. I usually find these tools sufficient for cleaning a mask, regardless of the amount of editing required.


.. _georef:

--------------------------------------
Does my mask need to be georeferenced?
--------------------------------------

Most masks are already produced in a GIS context and are already geographically referenced. However, *RivGraph* does not require that your mask image be georeferenced (e.g. a GeoTIFF). If you provide a mask without any georeference information, *RivGraph* will assign it a "dummy" projection in order to proceed. This has no effect on the network extraction. However, it is strongly advised that you provide a georeferenced mask. There are three primary reasons for this:

1) The coordinate reference system (CRS) of your mask will be carried through all your analysis, meaning that shapefiles and GeoTIFFs you export using *RivGraph* will align perfectly with your mask. Additionally, your results will be easily importable into a GIS for further analysis or fusion with other geospatial data.

2) *RivGraph* computes morphologic metrics (length and width) using pixel coordinates. A georeferenced mask contains information about the units of the mask, and thus any metrics of physical distance will inherit these units. If your CRS is meters-based, your results will be in meters.

3) Some of *RivGraph*'s functionality under the hood requires some heuristic thresholds or parameters. While these were designed to be as CRS-agnostic as possible, these functions will undoubtedly perform better when pixels have known units. As an example, generating a mesh along a braided river corridor requires some parameters defining the size and smoothness of the mesh. Having a mask with physically-meaningful units makes this parameterization much simpler and more intuitive.

.. warning::
  You should **avoid** degree-based CRSs (like EPSG:4326). This is because the length of a degree is not uniform, but varies with latitude. For example, at the equator, a degree of longitude is roughly 111 km. In Anchorage, Alaska, a degree of longitude is approximately 55 km. Effectively, degrees are meaningless units of physical measurements. A more prudent approach would be to first project your mask into a meters-based CRS (e.g. the appropriate `UTM zone <https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>`_) before analysis with *RivGraph*.

.. _nonriver:

---------------------------------------------------
Can my mask represent something that isn't a river?
---------------------------------------------------

Perhaps you'd like to vectorize a road network or a vascular system. This is possible to do with *RivGraph*. However, you will not be able to instantiate the convenient *delta* or *river* classes as they are designed only for river channel networks. Instead, you will need to poke around the API to figure out which functions will work for you. A good starting point is to skeletonize your mask with :obj:`rivgraph.mask_to_graph.skeletonize_mask()` then run :obj:`rivgraph.mask_to_graph.skel_to_graph()` to convert the skeleton to a set of links and nodes. If you have an interesting non-river use-case, please send an email to j........k@gmail.com and we can add it as an example.

.. _supportedfiletypes:

-----------------------------------------
What filetypes are supported for my mask?
-----------------------------------------

Any `gdal-readable filetype <https://gdal.org/drivers/raster/index.html>`_ should be fine. GeoTIFF is most common and recommended if possible.
.. _shoreline:

============================
Shoreline creation
============================

Every delta analyzed by *RivGraph* requires that the user create a shoreline as well. This shoreline is used to determine the location of the outlet nodes of the network. Here, guidance is provided for how to create this shoreline for your delta mask.

 - :ref:`whyshoreline`
 - :ref:`howshoreline`


.. _whyshoreline:

------------------------
Purpose of the shoreline
------------------------

Consider the following mask and its skeleton:

.. image:: ../../images/colville_mask_skel.PNG

How can we identify the "ends" of the network--i.e. the outlet locations? There are two options; first, we could manually specify each outlet point individually. This is a viable option, but it is also tedious and does not lend itself well to automation. Instead, *RivGraph* takes a different approach that still requires manual input, but is robust, less tedious, and has the potential for automation.

Instead of placing nodes individually, you will provide a shoreline shapefile (or any geopandas-readable PolyLine). *RivGraph* will intersect your shoreline with the skeleton (actually the vectorized links and nodes that make up the skeleton), place nodes at the connected intersection points, and trim away all the skeleton that lies outside the river network.

.. _howshoreline:

--------------------------
How do I make a shoreline?
--------------------------
There are numerous tools available to generate your shoreline, such as the `Opening Angle Method <http://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2008GL033963>`_ or you may automate a procedure. Here, we describe shoreline creation with QGIS, although other GIS software may also be used.

1. Generate the skeleton of your mask.

.. code-block:: python3

   from rivgraph.classes import delta
   mydelta = delta(Mask, path_to_results)
   mydelta.skeletonize()
   mydelta.to_geotiff('skeleton')

*RivGraph* will write your georeferenced skeleton at `path_to_results`.

2. Drag your mask and skeleton into QGIS. You will see something like the above figure.
3. Create a new layer in the ``Layer -> Create Layer -> Create Shapefile Layer`` dropdown.

    * Make sure to select ``Line`` for ``Geometry Type``.

    * Make sure to set the CRS of your new layer to be the same as your mask.

    * Specify a filepath, preferably the same as `path_to_results`, but it can be anywhere.

4. Turn on editing for this layer by clicking the ``Toggle Editing`` icon.
5. Create a new line feature using the ``Add Line Feature`` button.

Now we're at the point of actually drawing the shoreline. Our goal is to intersect all the skeleton links by the shoreline at locations where we'd like outlet nodes to be placed. Identify the first outlet of the mask, and begin drawing shoreline segments across your channels. Here are some tips:

    * The only thing that matters is where your shoreline intersects the outlet links. Don't worry if your shoreline doesn't actually follow the shore.
    * I find it helpful to connect the ends of islands to define where to create shoreline nodes; this typically ensures you're cutting across the correct part of the channel. See the figure below.
    * Make sure your final shoreline cuts the skeleton into two disconnected components.
    * Again, don't worry about intersecting portions of the skeleton other than the outlet links.
    * It gets easier with a little practice, and you may have to iterate a time or two to achieve your desired result.

.. image:: ../../images/shoreline_howto1.PNG
  :align: center

After we run ``mydelta.prune_network()`` (and specifying the proper paths for the shoreline and inlet nodes) with the shoreline above, we get the following result:

.. image:: ../../images/shoreline_howto2.PNG
  :align: center
.. centered::
  The pruned network is in blue; the outlet nodes are yellow.

Notice that all the spurious skeleton portions have been trimmed, as have all the links in the ocean. We also see outlet nodes have been placed exactly where the shoreline intersected the skeleton.
.. _background:

==========
Background
==========

---------------------
The Birth of RivGraph
---------------------

In 2016, I walked into the postdoc office at the Saint Anthony Falls Laboratory to see deltas on display as 8.5 x 11 sheets of paper 
taped together with hand-drawn markings all over. My colleagues Alejandro Tejedor and Anthony Longjas were trying to whip these deltas
into shape for a series of papers (`1 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014WR016577%4010.1002/%28ISSN%291944-7973.CONART1>`_, `2 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014WR016604>`_) introducing a new way to think about deltas: through the lens of their channel networks. As I 
looked at their handiwork, I couldn't help thinking that there must be a better way. My own research at the time had led me to develop 
`RivMAP <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016EA000196>`_, a Matlab toolbox for analyzing the morphodynamics of meandering rivers using binary masks. And so RivGraph was born as a small set of Matlab scripts with very limited functionality.

I finished my PhD with RivGraph in a barely-formed state, and followed my adviser, `Efi <http://efi.eng.uci.edu/>`_, to UC-Irvine for a short postdoc. Alex and Anthony were pushing the deltas work even further, and I kept adding scripts to RivGraph. The inevitable idea was hatched to analyze dozens of deltas worldwide, and eventually became a `successful proposal <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1812019>`_. In the meantime, I took a postdoc position at Los Alamos
National Laboratory with an eye toward global river morphodynamics. My arrival at LANL marked the beginning of my Python journey, and
soon I had converted all the RivGraph scripts to Python. My research needs led to a fuller development of RivGraph and brought braided
rivers into the mix, and sometime in 2017, RivGraph officially became RivGraph.

I cut my Python teeth while developing RivGraph, so there is a degree of inherited clunkiness and inefficiency baked in. Initially, there
were also some implements (such as automatically adding artificial nodes to parallel edges) that were tailored to somewhat particular use
cases. However, as the user base continues to grow, improvements and upgrades have been implemented to meet their needs. Please add your requests
to the mix and report any bugs you find using Github's `issue tracker <https://github.com/jonschwenk/RivGraph/issues>`_.

--------------------------------
What's different about RivGraph?
--------------------------------

RivGraph fills a void in coding space by providing a full-featured package for working with masks of river channels. 
In my opinion, most of RivGraph is just a convenient collection of already-existing functionality like CRS handling, path management, 
etc. That said, there are three somewhat novel components of RivGraph that are not available elsewhere, or at least not available in
Python. 

The first of these is the walking algorithm that breaks a skeleton into its constituent links and nodes. While there is similar
functionality available in GIS packages, RivGraph ensures that the resulting skeleton is parsimonious--i.e. contains as few nodes as possible
while fully preserving the topology of the mask. This is achieved with help from a double-convolution, where the first convolution identifies possible
branchpoints and the second reduces branchpoint clusters to a minimum-required set. 

The second novelty RivGraph offers is its ability
to automatically set flow directions. This is trivial for some networks (like Wax Lake Delta), but many channel networks are wild beasts. RivGraph's
solution is correspondingly complicated, but does a pretty good job. We `published  <https://esurf.copernicus.org/articles/8/87/2020/esurf-8-87-2020.html>`_ the method, its validation, and its implications.

The third novel component of RivGraph is its abilty to generate an along-channel mesh that approximately follows a river's centerline
while transecting the centerline approximately perpendicularly. This was a surprisingly tricky function to get right, and I'm not even
sure it's *there* yet. The first iteration of this appeared in RivMAP, and there have probably been 3-4 method changes before settling on the current version, which uses Dynamic Time Warping (thanks `Zoltan et. al <https://pubs.geoscienceworld.org/gsa/geology/article/47/3/263/568705/High-curvatures-drive-river-meandering>`_ for introducing this to me, although I still contend that it's inappropriate to use for measuring channel migration rates) to iteratively map vertices on buffered centerlines away from the original.

--------------------------------------
Is RivGraph only for channel networks?
--------------------------------------

While RivGraph is designed around channel networks, it contains a smattering of tools that can be useful across a broad range of analyses.
For instance, the `mask_to_graph.py <https://github.com/jonschwenk/RivGraph/blob/master/rivgraph/mask_to_graph.py>`_ script contains tools that will convert *any* binary mask to a vectorized skeleton, not just river
channel networks. There are a number of image processing tools in `im_utils.py <https://github.com/jonschwenk/RivGraph/blob/master/rivgraph/im_utils.py>`_ that I use frequently in other projects, like a Matlab-like
implementation of *regionprops()* for measuring blob properties of a binary image. I have personally found one of the most broadly useful tools in RivGraph is `write_geotiff <https://github.com/jonschwenk/RivGraph/blob/f2284de77a79b8f8812d04c579d52852f584de1d/rivgraph/io_utils.py#L281>`_, which does what it says.
 .. _contributing:

============
Contributing
============

If you think you're not skilled or experienced enough to contribute, think
again! We agree wholeheartedly with the sentiments expressed by this
`Imposter syndrome disclaimer <https://github.com/Unidata/MetPy#contributing>`_.

We welcome all forms of user contributions including feature requests, bug
reports, code, documentation requests, and code. Simply open an
`issue <https://github.com/jonschwenk/RivGraph/issues>`_. For code development
contributions, please contact us via email to be added to our slack channel
where we can hash out a plan for your contribution.
.. _featuredevelopment:

Feature Development
===================

There’s always more to be done! The features listed below are
suggestions that require some significant amount of development. If
you’re interested in tackling, or helping tackle, any of these, please
email j………k@gmail.com and we’ll add you to the RivGraph slack channel to
start discussions.

Automatic Shoreline Extraction
------------------------------

Currently, users must provide their own shoreline file for each delta
they wish to analyze. While the strategy for creating shorelines is
`documented <https://jonschwenk.github.io/RivGraph/shoreline/index.html>`__,
a preferable option would have RivGraph automatically generate
shorelines. There have been a number of published softwares that attempt
to solve this problem (e.g. opening-angle method), but in practice these
have been found to be too slow and/or need finegaling to interface with
RivGraph.

Lake/wetland connectivity
-------------------------

Many deltas, especially Arctic ones, are lake-dense, and the
connectivity of these lakes has implications for biogeochemical cycling
and transport times of water and sediment. We have spent a few weeks
developing a ``rivgraph-lakes`` branch that attempts to resolve lake
connectivity. However, it’s a difficult problem and we didn’t quite
cross the finish line.

.. figure:: ../../images/lakemask_example.PNG
   :alt: image-20220107142427137

   An example mask with labeled channel network (gray) and lakes (white).

In this formulation, a lake mask would be provided in conjunction with
the river mask, or lakes would be labeled distinctly in a river mask, as
shown by the above mask. The difficulties here lie in the number of ways
lakes can be connected; the figure shown is a rather simple case, but
sometimes lakes themselves intersect both the channel network *and* the
ocean/draining body. The many possibilities make simply resolving their
position within the channel network a formidable task.

The second major issue arises when trying to set flow directions. Are
lakes sources, sinks, or both? The flow directionality algorithms need
to know! And they need to be adapted accordingly. We have made
significant progress on this feature, but it is not ready for
Production.

Another point to mention here is that while we are focusing on lakes
here, the concept may be more generally applicable. For example, if a
user also had a mask of wetlands–or any objects that are connected to
the channel network, this framework could handle those cases.

Machine Learned Flow Directions
-------------------------------

The current scheme for `setting flow
directions <https://esurf.copernicus.org/articles/8/87/2020/>`__ in the
links of a channel network is quite complicated. It is
“physically-based” in the sense that the rules defining flow directions
have physical bases. There are two “hard” constraints–1) there can be no
interior sources or sinks and 2) there can be no cycles (although this
can be violated if RivGraph cannot find a cycle-less solutions).

The paper (and attached dataset) linked above contains plenty of
“training data” that could be used to take a machine learning approach
to set flow directions. Many of the algorithms in place, for example
determining if a link is part of a “main channel,” or
synthetically-generated slopes, can be used to generate the features of
an AI approach. It’s not clear to me how you would enforce the
constraints listed above, but I’ll bet it could be done. Another
difficulty might be that in the “physically-based” scheme, directions
are set iteratively starting with the most certain. By setting a link,
that information is useful in setting its neighbors’ directions, and
so-on. Continuity can only be exploited when only one link (in a group
of adjacent links) has an unknown direction. I would therefore guess
that a ML model would have to follow a similar iterative path.

However, who knows? An AI approach could be faster, especially for the
larger deltas. Post-processing corrections would also be useful (for
either approach, really).

Flow Direction Uncertainty
--------------------------

Many links’ flow directions are not certain and may in fact be
bidirectional. One approached proposed to handle this is: rather than
export a single adjacency matrix with RivGraph’s “best guess,” we could
export a family of adjacency matrices. Alej Tejedor has done some recent
work on the concept of “effective resolution”–i.e. the resolution of the
underlying mask at which graph-based flux partitioning is significantly
impacted (as we coarsen the mask). This concept could be combined with
the family of adjacency matrices to prevent the family from being huge.
In other words, “effective resolution” could prune links that don’t
contribute much to flux routing–and these links are typically the most
uncertain.

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "gallery\knights2020.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_gallery_knights2020.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_gallery_knights2020.py:


Nitrate Removal Across Ecogeomorphic Zones in Wax Lake Delta, Louisiana (USA)
=============================================================================
*Deon Knights, Audrey H. Sawyer, Rebecca T. Barnes, Anastasia Piliouras, 
Jon Schwenk, Douglas A. Edmonds, and Alexander M. Brown*

`This publication <https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2019WR026867>`_
made use of RivGraph's flux modeling to estimate nitrate uptake rates
in channels on the Wax Lake Delta. 

.. image:: ../gallery_source/images/knights_et_al_2020.PNG

.. GENERATED FROM PYTHON SOURCE LINES 13-16

.. code-block:: default

    # (left) Removal rate calculated using submerged‐delta approach. (right) Removal rate calculated using 
    # nutrient spiraling approach as a percentage of amount of nitrate entering each link.



.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  0.000 seconds)


.. _sphx_glr_download_gallery_knights2020.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: knights2020.py <knights2020.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: knights2020.ipynb <knights2020.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
:orphan:



.. _sphx_glr_gallery:

RivGraph in the wild
====================

Some selected publications that have used RivGraph to answer some cool science questions.

`Let us know <https://github.com/jonschwenk/RivGraph/issues>`_ if you'd like to add an application to the gallery!




.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="`This publication &lt;https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2019WR026867&gt;`_ made...">

.. only:: html

 .. figure:: /gallery/images/thumb/sphx_glr_knights2020_thumb.png
     :alt: Nitrate Removal Across Ecogeomorphic Zones in Wax Lake Delta, Louisiana (USA)

     :ref:`sphx_glr_gallery_knights2020.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /gallery/knights2020
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
.. _examples:

========
Examples
========

Delta Examples
--------------

.. toctree::
   :maxdepth: 1

   delta_example/delta_example

River Examples
--------------

.. toctree::
   :maxdepth: 1

   braided_river_example/braided_river_example
Let’s demo RivGraph on the Colville Delta!
------------------------------------------

This demo shows some of the core functionality and convenient plotting
and exporting features provided by RivGraph for analyzing delta channel
networks. The basic steps of RivGraph include:

1. Instantiate delta class
2. Skeletonize the binary mask
3. Compute the network (links and nodes)
4. Prune the network (requires user-created shoreline and input nodes
   for deltas)
5. Compute morphologic metrics (lengths, widths)
6. Assign flow directions for each link.
7. Compute some topologic metrics.

Along the way, we’ll export some geotiffs and GeoJSONs (or shapefiles if
you prefer) for inspection in QGIS. RivGraph requires a **binary mask of
the channel network**, preferably georeferenced (i.e., a GeoTiff). For
deltas, you will also need to create two shapefiles/GeoJSONs: one of the
**shoreline**, and one of the **inlet nodes**. See section 4 for
guidance on how to create these required geovector files.

1. Instantiate delta class
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    from rivgraph.classes import delta
    import matplotlib.pyplot as plt

    # Define the path to the georeferenced binary image.
    mask_path = "./data/Colville_Delta/Colville_mask.tif"

    # Results will be saved with this name
    name = 'Colville'

    # A folder called Colville will be created within this path for storing outputs
    results_folder = './data/Colville_Delta/Results'

    # Boot up the delta class! We set verbose=True to see progress of processing.
    colville = delta(name, mask_path, results_folder=results_folder, verbose=True)

    # The mask has been re-binarized and stored as an attribute of colville:
    plt.imshow(colville.Imask)




.. parsed-literal::

    <matplotlib.image.AxesImage at 0x273d8f8e3a0>




.. image:: output_2_1.png


2. Skeletonize the binary mask
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Simply use the skeletonize() method.
    colville.skeletonize()

    # After running, colville has a new attribute: Iskel. Let's take a look.
    plt.imshow(colville.Iskel)


.. parsed-literal::

    Skeletonizing mask...done.




.. parsed-literal::

    <matplotlib.image.AxesImage at 0x273d790afd0>




.. image:: output_4_2.png


The skeleton is hard to see; perhaps we’d like to look at it closer? One
option is to save it as a geotiff and pull it up in a GIS (like QGIS).

.. code:: ipython3

    # We use the write_geotiff() method with the "skeleton" option.
    colville.to_geotiff('skeleton')


.. parsed-literal::

    Geotiff written to data\Colville_Delta\Results\Colville_skel.tif.


The georeferenced Colville skeleton has been written to disk, so we can
pull it up in QGIS along with the georeferenced mask:

.. figure:: images/colville_qgis_mask_skel_large.png
   :alt: colville_qgis_mask_skel_large.PNG

   colville_qgis_mask_skel_large.PNG

Or a bit zoomed-in:

.. figure:: images/colville_qgis_mask_skel_zoom.png
   :alt: colville_qgis_mask_skel_zoom.PNG

   colville_qgis_mask_skel_zoom.PNG

3. Compute the network (links and nodes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Simply use the compute_network() method.
    colville.compute_network()


.. parsed-literal::

    Resolving links and nodes...done.


.. code:: ipython3

    # Now we can see that the "links" and "nodes" dictionaries have been added as colville attributes:
    links = colville.links
    nodes = colville.nodes
    print('links: {}'.format(links.keys()))
    print('nodes: {}'.format(nodes.keys()))


.. parsed-literal::

    links: dict_keys(['idx', 'conn', 'id', 'n_networks'])
    nodes: dict_keys(['idx', 'conn', 'id'])


The *links* dictionary currently contains four keys: - idx: a list of
all the pixel indices that make up the link (indices created with input
mask shape and np.ravel_multi_index) - conn : a two-element list
containing the node *id*\ s of the link’s endpoints - id: each link has
a unique *id*; the ordering is irrelevant - n_networks: the number of
disconnected networks (==1 if the input mask contains a single connected
blob)

The *nodes* dictionary currently contains three keys: - idx: the index
of the node’s position within the original image
(i.e. np.ravel_multi_index()) - conn: an N-element list containing the N
link *id*\ s of the links connected to this node. - id: each node has a
unique *id*; the ordering is irrelevant

We can visualze the network in a couple of ways. First, we can plot with
matplotlib:

.. code:: ipython3

    colville.plot('network')



.. image:: output_12_0.png


Nodes and links are labeled with their ids. Kind of hard to see, so we
can zoom in OR we can export the network to geovectors and pull ’em into
QGIS:

.. code:: ipython3

    colville.to_geovectors('network', ftype='json') # ftype can be either 'shp' or 'json'

    # Let's see where the network geovector files were written:
    print(colville.paths['links'])
    print(colville.paths['nodes'])


.. parsed-literal::

    data\Colville_Delta\Results\Colville_links.json
    data\Colville_Delta\Results\Colville_nodes.json


And dragging these into QGIS: |colville_network_unpruned.PNG| You can
query different links and nodes using the Identify tool. Note that their
properties (‘conn’ and ‘id’) are appended.

.. |colville_network_unpruned.PNG| image:: images/colville_network_unpruned.png

4. Pruning the network
~~~~~~~~~~~~~~~~~~~~~~

You notice in the above image that there are many superfluous links
along the shoreline. This is a result of skeletonizing such a massive,
connected waterbody (i.e. the ocean in this case). Additionally, the
network contains a number of “dangling” links, or those that are
connected only at one end. We want to keep the inlet and outlet dangling
links, but not the others! RivGraph will automatically prune the
network, but it requires (for deltas) two additional pieces of
information: the location of the inlet nodes, and a delineation of the
shoreline. We can create both of these in QGIS:

.. figure:: images/colville_shoreline_inlet_outlet.png
   :alt: colville_shoreline_inlet_outlet.png

   colville_shoreline_inlet_outlet.png

Shoreline: Create a polyline vector layer. The shoreline should be drawn
to intersect all the outlet links. It should separate all the unwanted
ocean links from the actual links of the delta channel network. If you
get errors, you may need to adjust your shoreline a little–try to ensure
it does not intersect any nodes!

Inlet nodes: Create a point vector layer. Simply place points at nodes
that represent the inlets to the network. The placement does not need to
be exact; RivGraph will find the closest node to the one(s) you create.
These will be marked as inlet nodes and won’t be removed during pruning.

Saving: For convencience, these files should be saved in the Results
folder that you initialized the class. Save as
results_folder/Colville_shoreline.shp and
results_folder/Colville_inlet_nodes.shp. However, this is not mandatory
as you can also point to the files during pruning.

Now that we have identified the shoreline and inlet/outlet nodes, let’s
prune the network!

.. code:: ipython3

    colville.prune_network()
    # Note that we can also specify the location of the shoreline and inlet nodes:
    # colville.prune_network(path_shoreline='/path/to/shoreline/file', path_inletnodes='/path/to/inletnodes/file')

    # Now that we've pruned, we should re-export the network:
    colville.to_geovectors()
    # Note that this time we didn't specify the arguments; by default 'network' will be exported as type 'json'.



.. parsed-literal::

    [917, 919, 923, 926, 927, 930, 931, 933, 935, 938, 939, 941, 944, 946, 948, 950, 951, 954, 955, 958, 959, 962, 963]


Let’s see how the pruned version compares to the unpruned:

.. figure:: images/colville_shoreline_inlet_outlet_pruned.png
   :alt: colville_shoreline_inlet_outlet_pruned.png

   colville_shoreline_inlet_outlet_pruned.png

Wow, we really clipped off a lot of links! We also added some new nodes
at the shoreline–notice how each link that intersects the shoreline was
truncated, and outlet nodes were placed there (RivGraph remembers which
nodes are outlet nodes). You may be concerned that some of the dangling
links or subnetworks were pruned–this is by design, and if you want to
retain any dangling links, you need to mark their upstream-most nodes as
inlet nodes in your shapefile.

Compare with the figure above this one; the set of nodes was also
reduced. As links were removed from the network, some nodes were no
longer needed as they only connected two links.

5. Compute morphologic metrics (lengths, widths)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that the network is resolved and pruned, we can compute some link
metrics.

.. code:: ipython3

    # Compute link widths and lengths
    colville.compute_link_width_and_length()

    # Lets look at histograms of link widths and lengths:
    trash = plt.hist(colville.links['len_adj'], bins=50)
    plt.ylabel('count')
    plt.xlabel('link length (m)')
    plt.title('Histogram of link lengths')


.. parsed-literal::

    Computing link widths and lengths...done.




.. parsed-literal::

    Text(0.5, 1.0, 'Histogram of link lengths')




.. image:: output_22_2.png


In the above figure, we see that almost all the links are 1 km or
shorter, with three being much longer. This histogram will be different
for each delta, and can depend on the resolution of your input binary
mask.

Note: the lengths are reported in meters because that is the unit of the
original geotiff CRS. You can check this unit with
``print(colville.unit)``. It is highly unadvisable to use degrees
(EPSG:4326 and others) to compute distances.

.. code:: ipython3

    print(colville.unit)


.. parsed-literal::

    meter


Note: we used the ‘len_adj’ field rather than the ‘len’ field. The
difference is addressed in a separate Jupyter notebook called XXX.

We can do the same for the widths:

.. code:: ipython3

    trash = plt.hist(colville.links['wid_adj'], bins=50)
    plt.ylabel('count')
    plt.xlabel('link width (m)')
    plt.title('Histogram of link widths')




.. parsed-literal::

    Text(0.5, 1.0, 'Histogram of link widths')




.. image:: output_26_1.png


6. Assign flow directions for each link.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we wish to determine the long-term, steady-state flow direction in
each link. The algorithms used here are described in `this
paper <https://www.earth-surf-dynam.net/8/87/2020/esurf-8-87-2020.html>`__.

.. code:: ipython3

    colville.assign_flow_directions()


.. parsed-literal::

    A file has been created for manually setting link directions at data\Colville_Delta\Results\Colville_fixlinks.csv.
    No cycles were found in network.


If RivGraph has any problems assigning link directions, it will let us
know. Here, we see no error messages, and a message indicating no cycles
were found in the graph. Great!

We also notice that RivGraph mentiones that a .csv file was created for
us to manually set flow directions. If we inspect the flow directions
and find some that are incorrect, these can be fixed by entering the
link ID and the appropriate upstream node in this .csv, and running
``assign_flow_directions()`` again. See the `braided river
example <https://github.com/jonschwenk/RivGraph/blob/master/examples/braided_river_example.ipynb>`__,
section 7 for more details. Note that any links entered into this .csv
will be forced to have the upstream node as indicated. RivGraph sets
links’ directions iteratively, so if you find a problematic area in the
link directions (i.e. a number of links whose directions are wrong), you
can usually fix it by setting a few key links without needing to flip
all of them manually.

Let’s look at some plots.

.. code:: ipython3

    # Plot the links with the directionality marked
    colville.plot('directions')



.. image:: output_30_0.png


Links are colored such that upstream is cyan and downstream is purple.
Similar to the skeleton, we can export the link directions as a geotiff
for inspection in a GIS:

.. code:: ipython3

    colville.to_geotiff('directions')


.. parsed-literal::

    Geotiff written to data\Colville_Delta\Results\Colville_link_directions.tif.


Pulling this into QGIS and applying a similar color ramp, we see

.. figure:: images/colville_link_directions.PNG
   :alt: colville_link_directions.PNG

   colville_link_directions.PNG

The pixel values along each link have been rescaled from 0 (upstream) to
1 (downstream).

Now that flow directions have been computed, we can also compute
junction angles at each node.

.. code:: ipython3

    # As of 3/4/2020, this method only computes junction angles at nodes that have exactly three connecting links.
    colville.compute_junction_angles(weight=None) # See XXX for a description and meaning of the weight options.

    # If we check the the nodes dictionary, we should see that three new fields exist: 'int_ang', 'jtype', and 'width_ratio'.
    # 'int_ang' is the junction angle. 'jtype' is either 'b' (bifurcation), 'c' (confluence), or -1 for nodes for which the
    # junction angles cannot be computed. 'width_ratio' refers to the ratio between the larger and smaller links.
    print(colville.nodes.keys())


.. parsed-literal::

    dict_keys(['idx', 'conn', 'id', 'inlets', 'outlets', 'int_ang', 'jtype', 'width_ratio'])


7. Compute topologic metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RivGraph will compute a number of topologic metrics for your delta
channel network. These metrics are explained and demonstrated in Tejedor
et. al 2015a (doi.org/10.1002/2014WR016577)
and
2015b (doi.org/10.1002/2014WR016604).
Note that some pre-processing is done to the topology to compute these
metrics; it is highly recommended that you understand these
preprocessing steps and/or compute the metrics yourself.

.. code:: ipython3

    colville.compute_topologic_metrics() # You may get an overflow warning

    # The metrics are stored in an attribute dictionary:
    print(colville.topo_metrics.keys())


.. parsed-literal::

    dict_keys(['nonlin_entropy_rate', 'nER_prob_exceedence', 'nER_randomized', 'top_mutual_info', 'top_conditional_entropy', 'top_link_sharing_idx', 'n_alt_paths', 'resistance_distance', 'top_pairwise_dependence', 'flux_sharing_idx', 'leakage_idx', 'dyn_pairwise_dependence', 'dyn_mutual_info', 'dyn_conditional_entropy'])


.. code:: ipython3

    # Query different metrics by accessing the dictionary by key.
    print(colville.topo_metrics['nonlin_entropy_rate'])


.. parsed-literal::

    0.7623661979554095


.. code:: ipython3

    # Most metrics are computed for each outlet node
    print(colville.topo_metrics['top_mutual_info']) # The first column are node IDs, the second are the topological mutual information values.


.. parsed-literal::

    [[154.           4.48446645]
     [155.           4.38451158]
     [156.           4.4265762 ]
     [157.           4.87862286]
     [158.           4.3490749 ]
     [159.           4.83200719]
     [160.           4.27096245]
     [161.           4.33538991]
     [162.           4.33403524]
     [163.           3.71888165]
     [164.           4.3345187 ]
     [165.           3.69583589]
     [166.           4.4236153 ]
     [167.           4.3345187 ]
     [168.           4.39590258]
     [169.           3.79725696]
     [170.           4.40327503]
     [171.           4.41122671]
     [172.           4.42854618]
     [173.           4.46143106]
     [174.           4.46675082]
     [175.           4.32908791]
     [176.           4.33634811]]


If you wish to compute your own metrics or perform topological analyses,
you’ll probably need an adjacency matrix. RivGraph will provide this
with the following method:

.. code:: ipython3

    # Unweighted, unnormalized adjacency matrix
    adj = colville.adjacency_matrix()
    print(adj)


.. parsed-literal::

    [[0. 0. 0. ... 0. 0. 0.]
     [1. 0. 0. ... 0. 0. 0.]
     [0. 0. 0. ... 0. 0. 0.]
     ...
     [0. 0. 0. ... 0. 0. 0.]
     [0. 0. 0. ... 0. 0. 0.]
     [0. 0. 0. ... 0. 0. 0.]]


.. code:: ipython3

    # You may also want an adjacency matrix weighted by link width.
    adj_w = colville.adjacency_matrix(weight='wid_adj') # Can also weight by 'len_adj' or provide a vector of your own weights.
    print(adj_w)


.. parsed-literal::

    [[ 0.          0.          0.         ...  0.          0.
       0.        ]
     [95.23358675  0.          0.         ...  0.          0.
       0.        ]
     [ 0.          0.          0.         ...  0.          0.
       0.        ]
     ...
     [ 0.          0.          0.         ...  0.          0.
       0.        ]
     [ 0.          0.          0.         ...  0.          0.
       0.        ]
     [ 0.          0.          0.         ...  0.          0.
       0.        ]]


.. code:: ipython3

    # And you may want this adjacency matrix normalized.
    adj_w_n = colville.adjacency_matrix(weight='wid_adj', normalized=True)
    print(adj_w_n) # Each row sums to 1


.. parsed-literal::

    [[0. 0. 0. ... 0. 0. 0.]
     [1. 0. 0. ... 0. 0. 0.]
     [0. 0. 0. ... 0. 0. 0.]
     ...
     [0. 0. 0. ... 0. 0. 0.]
     [0. 0. 0. ... 0. 0. 0.]
     [0. 0. 0. ... 0. 0. 0.]]
Let’s demo RivGraph on the Brahmaputra River!
---------------------------------------------

This demo shows some of the core functionality and convenient plotting
and exporting features provided by RivGraph for a braided river network.
The basic steps of RivGraph-ing a braided river include:

1. Instantiate river class
2. Skeletonize the binary mask
3. Compute the network (links and nodes)
4. Prune the network
5. Compute morphologic metrics (lengths, widths)
6. Compute a mesh (6.1 - Adjust mesh parameters)
7. Assign flow directions for each link
8. A note on topologic metrics

Along the way, we’ll export some geotiffs and GeoJSONs (or shapefiles if
you prefer) for inspection in QGIS. RivGraph requires a **binary mask of
the channel network**, preferably georeferenced (i.e., a GeoTiff) in a
projected coordinate reference system.

1. Instantiate river class
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    from rivgraph.classes import river
    import matplotlib.pyplot as plt
    import os
    
    # Define the path to the georeferenced binary image.
    mask_path = './data/Brahmaputra_Braided_River/Brahmaputra_mask.tif'
    
    # Results will be saved with this name
    name = 'Brahma' 
    
    # Where to store RivGraph-generated geotiff and geovector files.
    results_folder = './data/Brahmaputra_Braided_River/Results' 
    
    # Set the exit sides of the river relative to the image. In this case, the
    # Brahmaputra is "entering" the image from the North and "exiting" the 
    # image from the South.
    es = 'NS' # The first character is the upstream side
    
    # Boot up the river class! We set verbose=True to see progress of processing.
    brahma = river(name, mask_path, results_folder, exit_sides=es, verbose=True) 
    
    # The mask has been re-binarized and stored as an attribute of brahma:
    plt.imshow(brahma.Imask)




.. parsed-literal::

    <matplotlib.image.AxesImage at 0x1eb4f279e50>




.. image:: output_2_1.png


2. Skeletonize the binary mask
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Simply use the skeletonize() method.
    brahma.skeletonize()
    
    # The skeletonized image is stored as an attribute to the brahm class. Let's take a look.
    plt.imshow(brahma.Iskel)




.. parsed-literal::

    <matplotlib.image.AxesImage at 0x1eb4f7e90a0>




.. image:: output_4_1.png


The skeleton is hard to see; perhaps we’d like to look at it closer? One
option is to save it as a geotiff and pull it up in a GIS (like QGIS).

.. code:: ipython3

    # We use the write_geotiff() method with the 'skeleton' option.
    brahma.to_geotiff('skeleton')


.. parsed-literal::

    Geotiff written to data\Brahmaputra_Braided_River\Results\Brahma_skel.tif.


The georeferenced Brahmaputra skeleton has been written to disk, so we
can pull it up in QGIS along with the georeferenced mask:

.. figure:: images/brahma_qgis_mask_skel.png
   :alt: brahma_qgis_mask_skel.PNG

   brahma_qgis_mask_skel.PNG

3. Compute the network (links and nodes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Simply use the compute_network() method.
    brahma.compute_network()


.. parsed-literal::

    Resolving links and nodes...done.


.. code:: ipython3

    # Now we can see that the "links" and "nodes" dictionaries ahve been added 
    # as attributes to the brahma class:
    links = brahma.links
    nodes = brahma.nodes
    print('links: {}'.format(links.keys()))
    print('nodes: {}'.format(nodes.keys()))


.. parsed-literal::

    links: dict_keys(['idx', 'conn', 'id', 'n_networks'])
    nodes: dict_keys(['idx', 'conn', 'id'])


The *links* dictionary currently contains four keys: - idx: a list of
all the pixel indices that make up the link (indices created with input
mask shape and np.ravel_multi_index) - conn : a two-element list
containing the node *id*\ s of the link’s endpoints - id: each link has
a unique *id*; the ordering is irrelevant - n_networks: the number of
disconnected networks (==1 if the input mask contains a single connected
blob)

The *nodes* dictionary currently contains three keys: - idx: the index
of the node’s position within the original image
(i.e. np.ravel_multi_index()) - conn: an N-element list containing the N
link *id*\ s of the links connected to this node. - id: each node has a
unique *id*; the ordering is irrelevant

We can visualze the network in a couple of ways. First, we can plot with
matplotlib:

.. code:: ipython3

    brahma.plot('network')



.. image:: output_12_0.png


Nodes and links are labeled with their ids. We can zoom in if plotting
in an interactive matplotlib window, *or* we can export the network
links and nodes as geovector files and pull ’em into QGIS:

.. code:: ipython3

    brahma.to_geovectors('network', ftype='json')
    
    # Let's see where the network geovector files were written:
    print(brahma.paths['links'])
    print(brahma.paths['nodes'])


.. parsed-literal::

    data\Brahmaputra_Braided_River\Results\Brahma_links.json
    data\Brahmaputra_Braided_River\Results\Brahma_nodes.json


And dragging these into QGIS: |brahma_qgis_network_unpruned.PNG| You can
query different links and nodes using the Identify tool. Note that their
properties (‘conn’ and ‘id’) are appended.

.. |brahma_qgis_network_unpruned.PNG| image:: images/brahma_qgis_network_unpruned.png

If you pan around a bit, you’ll notice that there are many “dangling”
links, or links connected only at one end. These links play no role in
the connectivity of the network, and we therefore may want to remove
them. We don’t want to remove the inlet and outlet links, however.
RivGraph will prune these dangling links but preserve the inlet and
outlet links by exploiting the *exit_sides* we supplied earlier. Let’s
prune the network.

4. Prune the network
~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    brahma.prune_network()
    
    # We see that 'inlets' and 'outlets' have been added to the nodes dictionary:
    print(brahma.nodes.keys())
    
    # We can get the node ids of the inlets and outlets
    print('inlets:', brahma.nodes['inlets'])
    print('outlets:', brahma.nodes['outlets'])


.. parsed-literal::

    dict_keys(['idx', 'conn', 'id', 'inlets', 'outlets'])
    inlets: [247]
    outlets: [1919]


5. Compute morphologic metrics (lengths, widths)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that the network is resolved and pruned, we can compute some link
metrics.

.. code:: ipython3

    brahma.compute_link_width_and_length()
    
    # Let's look at histograms of link widths and lengths:
    trash = plt.hist(brahma.links['len_adj'], bins=50)
    plt.ylabel('count')
    plt.xlabel('link length (m)')
    plt.title('Histogram of link lengths')


.. parsed-literal::

    Computing distance transform...done.
    Computing link widths and lengths...done.




.. parsed-literal::

    Text(0.5, 1.0, 'Histogram of link lengths')




.. image:: output_20_2.png


In the above figure, we see that almost all the links are shorter than 5
km. This histogram will be different for each braided river, and can
depend on the resolution of your input binary mask. Resolving smaller
channels generally, though not always, produces smaller average link
lengths as longer links are broken to connect to the smaller ones.

**Note**: the lengths are reported in **meters** bcause that is the unit
of the provided mask’s CRS (coordinate reference system). You can check
this unit:

.. code:: ipython3

    print(brahma.unit)


.. parsed-literal::

    meter


6. Compute a mesh
~~~~~~~~~~~~~~~~~

In preparation for setting flow directions of each link, we will compute
an “along-valley” mesh. This mesh is created based on the overall
morphology of the braided river as opposed to individual channels within
the network. The objective of the mesh is to create an along-channel
grid that contains transects that are roughly perpendicular to the river
“centerline” and approximately evenly-spaced.

While this mesh is necessary for setting flow directions, it is also
useful for measuring along-river characteristics (like width, erosion
rates, braiding index, etc.). The generation of this mesh requires a few
parameters that ideally would be user-defined each time. However, for
simplicity, RivGraph will make estimates of these parameters based
primarily off the average link width. These parameters are described
after we generate a default mesh.

.. code:: ipython3

    # Note that we provide no arguments to the compute_mesh() function.
    brahma.compute_mesh()


.. parsed-literal::

    Computing centerline...done.
    Computing link widths and lengths...done.
    Generating mesh...done.


Before we play with the mesh-generation parameters, let’s take a look at
what we’ve generated with the ``compute_mesh()`` function.

First, we see that a centerline was computed. We can access this
centerline:

.. code:: ipython3

    print(brahma.centerline)


.. parsed-literal::

    (array([764625., 764595., 764565., ..., 785475., 785505., 785535.]), array([2753535., 2753505., 2753475., ..., 2633625., 2633655., 2633655.]))


We get a numpy array of two arrays of columns, rows of the centerline.
This isn’t very interpretable as-is, but we can export the centerline as
a geovector file:

.. code:: ipython3

    brahma.to_geovectors('centerline', ftype='json')

The centerline is exported as a shapefile or GeoJSON, depending on the
filetype you specify. (GeoJSON is default.) Let’s take a look at this
centerline in QGIS:

.. figure:: images/brahma_qgis_centerline.png
   :alt: brahma_qgis_centerline.png

   brahma_qgis_centerline.png

A little jagged, but it’s fine for our purposes. The centerline is
computed by filling in all the islands of the network, then using a
distance transform to find the “centermost” pixels.

Now let’s take a look at the mesh we generated. We first need to export
it:

.. code:: ipython3

    brahma.to_geovectors('mesh', ftype='json')

The mesh consists of two files: ``meshlines``, or transects, and
``meshpolys``. If we want to see where these files were generated, we
can check the paths dictionary:

.. code:: ipython3

    print(brahma.paths['meshlines'])
    print(brahma.paths['meshpolys'])


.. parsed-literal::

    data\Brahmaputra_Braided_River\Results\Brahma_meshlines.json
    data\Brahmaputra_Braided_River\Results\Brahma_meshpolys.json


Let’s see what the mesh looks like by dragging these GeoJSONs into QGIS:

.. figure:: images/brahma_qgis_initial_meshlines.png
   :alt: brahma_qgis_initial_meshlines.png

   brahma_qgis_initial_meshlines.png

The mesh does a pretty good job of meeting our two criteria:
approximately perpendicular to the centerline, and evenly-spaced along
the centerline. It’s not perfect, but we can play with the parameters to
adjust it.

6.1 Adjust mesh parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~

We may want to alter certain features of the mesh, like the spacing of
transects or the overall width of the mesh. When we call
``compute_mesh()``, there are three optional keyword arguments we can
provide to control mesh properties. These are:

``grid_spacing`` : The along-centerline distance between each transect.
The default value is the average link width, found by
``np.mean(brahm.links['wid_adj'])``.

``smoothing`` : The degree of smoothing to perform on the centerline
before creating its offsets. This parameter is expressed in terms of the
fraction of the total centerline length, so the default
``smoothing=0.1`` will use a moving-average with a window size of 10% of
the length of the centerline.

``buf_halfwidth`` : The distance from the centerline to each edgeline of
the buffer. The default is 10% wider than the maximum width of the mask,
which ensures that the mask is fully covered by the mesh.

We can check what values of each of these were used above:

.. code:: ipython3

    # grid_spacing
    print('grid_spacing: {}'.format(brahma.avg_chan_width))
    # buf_halfwidth
    print('buf_halfwidth: {}'.format(brahma.max_valley_width_pixels * brahma.pixlen * 1.1))
    # smoothing by default was 0.1
    print('smoothing: {}'.format(0.1))



.. parsed-literal::

    grid_spacing: 623.1819853141467
    buf_halfwidth: 15631.83213830036
    smoothing: 0.1


You may have noticed that the mesh transects near the beginning of the
reach (top of image) aren’t quite a perpendicular to the centerline as
we might like. Let’s try smoothing the centerline a bit more and
reducing the buffer width to make these more perpendicular. The grid
spacing will not affect the overall mesh strucutre, so we’ll leave it
the same for comparison purposes.

.. code:: ipython3

    brahma.compute_mesh(buf_halfwidth=10000, smoothing=0.25)


.. parsed-literal::

    Generating mesh...done.


Here’s a side-by-side comparison with the default mesh:

.. figure:: images/brahma_mesh_comparison.png
   :alt: brahma_mesh_comparison.png

   brahma_mesh_comparison.png

There are two primary differences between our custom mesh and the
default one. First, the custom mesh is narrower as a direct result of us
reducing ``buf_halfwidth`` from 15631 to 10000. Secondly, the custom
mesh transects are significantly more parallel to the centerline than
the default mesh. This resulted from increasing ``smoothing`` from 0.1
to 0.25. These two parameters, plus ``grid_spacing``, can be used to
design a mesh to suit your needs.

7. Assign flow directions to each link
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we want to determine the long-term, steady-state flow direction in
each link. The algorithms used here are described in `this
paper <https://www.earth-surf-dynam.net/8/87/2020/esurf-8-87-2020.html>`__.

.. code:: ipython3

    brahma.assign_flow_directions()


.. parsed-literal::

    Setting link directionality...Using data\Brahmaputra_Braided_River\Results\Brahma_fixlinks.csv to manually set flow directions.
    Attempting to fix 2 cycles.
    All cycles were resolved.
    done.


A statement appears that tells us that a .csv file was created for us to
set flow directions manually. This file is only created if it does not
already exist. So if we re-run ``assign_flow_directions()``, we should
not see this message again.

.. code:: ipython3

    brahma.assign_flow_directions()


.. parsed-literal::

    Setting link directionality...Using data\Brahmaputra_Braided_River\Results\Brahma_fixlinks.csv to manually set flow directions.
    Attempting to fix 2 cycles.
    All cycles were resolved.
    done.


Now we see a message stating that RivGraph will use the .csv file to set
flow directions manually. However, we haven’t populated the file yet so
no there is no information for RivGraph to use. We’ll revisit this at
the end of this section.

We also see a message stating that RivGraph is
``Attempting to fix 2 cycles.`` A cycle is a set of links within a graph
that feeds into itself. RivGraph’s directionality algorithms enforce the
condition that no cycles should be present in the resulting graph. If
any cycles are found after setting all link directions, RivGraph
attempts to fix them. This does not mean that all the links’ directions
are correct, but rather that the resulting graph contains no cycles. In
this case, RivGraph was able to resolve both cycles.

Let’s inspect the resulting link directions. We again have two options:
we can plot the links via matplotlib within Python:

.. code:: ipython3

    brahma.plot('directions')



.. image:: output_44_0.png


Cyan represents the upstream portion of a link, and magenta the
downstream. As before, this is difficult to see without zoooming. So we
can export a geotiff that contains directionality information for
inspection in QGIS:

.. code:: ipython3

    brahma.to_geotiff('directions')


.. parsed-literal::

    Geotiff written to data\Brahmaputra_Braided_River\Results\Brahma_link_directions.tif.


Dragging into QGIS (and adding a colorbar legend), we see

.. figure:: images/brahma_qgis_initial_directions.png
   :alt: brahma_qgis_initial_directions.png

   brahma_qgis_initial_directions.png

I’ve circled two short links in yellow, and noted that their flow
direction was set as going right-to-left. Ideally, this junction would
be comprised of a single node replacing these links, but RivGraph does
not have the capability yet to simplify the network (`a feature request
has been made for
this) <https://github.com/jonschwenk/RivGraph/issues/11>`__. From
cursory inspection, we could make the argument that flow should instead
go left-to-right for these two links. Let’s force flow the opposite
direction as an example of how to manually set links.

First, we need to identify the link ID’s of each of these links, as well
as the desired upstream nodes. Using the Identify tool in QGIS with the
links and nodes GeoJSON layers turned on, this is easy:

.. figure:: images/brahma_qgis_identify_links_for_reversal.png
   :alt: brahma_qgis_identify_links_for_reversal.png

   brahma_qgis_identify_links_for_reversal.png

We see that the red-highlighted link’s ID is ``2280``, and it’s
(upstream, downstream) nodes are ``1631, 1650``. We want to reverse this
order so that the upstream node ID is ``1650``. Open the
``Brahma_fixlinks`` csv and enter this information. I have also done
this for the link to the right of this one.

.. figure:: images/brahma_fixlinks_csv.png
   :alt: brahma_fixlinks_csv.png

   brahma_fixlinks_csv.png

Simply save the .csv and re-run the ``assign_flow_directions()`` method:

.. code:: ipython3

    brahma.assign_flow_directions()


.. parsed-literal::

    Setting link directionality...Using data\Brahmaputra_Braided_River\Results\Brahma_fixlinks.csv to manually set flow directions.
    Attempting to fix 2 cycles.
    All cycles were resolved.
    done.


After exporting the ``directions`` geotiff again, we can plot it against
the original:

.. figure:: images/brahma_flow_direction_reversed_manually.png
   :alt: brahma_flow_direction_reversed_manually.png

   brahma_flow_direction_reversed_manually.png

We see that RivGraph used the information in our .csv to manually set
the flow directions for both these links in the opposite direction of
the original solution. It is important to note that RivGraph guarantees
that whatever link flow directions are specified in the .csv will be
preserved throughout the direction-setting process. RivGraph uses an
iterative approach to set link directionalities, where knowing the
direction of a nearby link can be used to set links of unknown
directions. This means that if you incorrectly set a link’s direction
manually, you could “infect” the nearby links with wrong directions!
Indeed, we see from the above comparison that the upper-right link’s
direction has also been reversed, even though we didn’t specify it.
Whether or not this is correct is determined by the user.

8. A note on topologic metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you’ve looked through the `delta
example <https://github.com/jonschwenk/RivGraph/blob/master/examples/delta_example.ipynb>`__,
you’ll see the final section covers computing topolgic metrics. In order
to compute these metrics, some additional finagling of the network is
required. We have not yet implemented the required pre-processing for
braided rivers. However, many of the functions in the `delta metrics
script <https://github.com/jonschwenk/RivGraph/blob/master/rivgraph/deltas/delta_metrics.py>`__
can be used on braided rivers, provided you first pre-process your
braided river network properly.
RivGraph in the wild
====================

Some selected publications that have used RivGraph to answer some cool science questions.

`Let us know <https://github.com/jonschwenk/RivGraph/issues>`_ if you'd like to add an application to the gallery!

.. _rivgraph:

========
rivgraph
========

.. automodule:: rivgraph.classes
   :members:
   :special-members:

.. automodule:: rivgraph.directionality
   :members:
   :special-members:

.. automodule:: rivgraph.geo_utils
   :members:
   :special-members:

.. automodule:: rivgraph.im_utils
   :members:
   :special-members:

.. automodule:: rivgraph.io_utils
   :members:
   :special-members:

.. automodule:: rivgraph.ln_utils
   :members:
   :special-members:

.. automodule:: rivgraph.mask_to_graph
   :members:
   :special-members:

.. automodule:: rivgraph.mask_utils
   :members:
   :special-members:

.. automodule:: rivgraph.walk
   :members:
   :special-members:
.. _rivers:

========
rivers
========

.. automodule:: rivgraph.rivers.river_directionality
   :members:
   :special-members:

.. automodule:: rivgraph.rivers.river_utils
   :members:
   :special-members:
.. _apiref:

=============
API Reference
=============

.. toctree::
   :maxdepth: 3
   :caption: Module Contents:

   rivgraph
   deltas
   rivers 
.. _deltas:

========
deltas
========

.. automodule:: rivgraph.deltas.delta_directionality
   :members:
   :special-members:

.. automodule:: rivgraph.deltas.delta_metrics
   :members:
   :special-members:

.. automodule:: rivgraph.deltas.delta_utils
   :members:
   :special-members:
.. _install:

=========================
Installation Instructions
=========================

.. note::
   *RivGraph* requires the installation of common geospatial Python packages such as `GDAL <https://gdal.org/>`_.
   These packages can be difficult to install properly and often create dependency errors.
   Because of this, we recommend using `Anaconda <https://www.anaconda.com/products/individual>`_ to create a virtual environment for *RivGraph*, and to manage the installation of Python libraries as it will handle package versions and dependencies for you.

Installation via *conda*
--------------------------

The latest 'stable' version of *RivGraph* can be installed via `conda`.
We recommend installing *RivGraph* into a fresh conda environment to minimize the risk of dependency clashes.
The easiest way to do this is by first downloading the `environment.yml <https://github.com/jonschwenk/RivGraph/blob/master/environment.yml>`_ (go to link, click "Raw", then copy the contents into a text editor and save as 'environment.yml'), opening Terminal (Mac/Unix) or Anaconda Prompt (Windows) and typing:
::

   $ conda env create --file /path/to/environment.yml

.. tip::
   The default environment name is 'rivgraph' (per the `environment.yml` file), but you can change the environment file to name it anything you like.

If you would rather install *RivGraph* into a pre-existing environment "myenv", you can use the following commands:
::

   conda activate myenv
   conda install rivgraph -c jschwenk

.. warning::

 *RivGraph* dependencies may be pinned to specific versions of packages (e.g. geopandas 0.7) that may not mesh with your existing environment.
 Check the `environment file <https://github.com/jonschwenk/RivGraph/blob/master/environment.yml>`_ for these cases.

Installation from source
------------------------

If you would prefer to install the *RivGraph* package from source, then follow these steps:

.. warning::

   *Rivgraph* uses many geospatial dependencies (e.g. GDAL) that can be
   difficult to install. Note that so long as the continuous integration
   workflows on GitHub are working (denoted by a green check next to the latest
   commit message on the source repository home page), the latest version of
   the source code is stable and installation from source will work if the
   dependencies are correctly installed.

1. Clone the repository
::

   $ git clone https://github.com/jonschwenk/RivGraph.git

2. Install dependencies; note these can be installed via conda from the
`environment.yml <https://github.com/jonschwenk/RivGraph/blob/master/environment.yml>`_ file, however a list is also
provided below with links to the homepage for each dependency.

**RivGraph Dependencies:**
   - `python <https://www.python.org/>`_ (tested on v3.6)
   - `GDAL <https://gdal.org/>`_
   - `NumPy <https://numpy.org/>`_
   - `GeoPandas <https://geopandas.org/>`_ (v0.7.0)
   - `scikit-image <https://scikit-image.org/>`_
   - `OpenCV <https://github.com/skvark/opencv-python>`_
   - `NetworkX <https://networkx.org/>`_
   - `Matplotlib <https://matplotlib.org/>`_
   - `pyproj <https://pyproj4.github.io/pyproj/stable/>`_
   - `Shapely <https://shapely.readthedocs.io/en/latest/>`_
   - `Fiona <https://fiona.readthedocs.io/en/latest/>`_
   - `FastDTW <https://github.com/slaypni/fastdtw>`_

3. From the cloned folder, run the following in the command line:
::

   $ python setup.py install

to install the *RivGraph* package.

.. note::
   If you run into issues installing *RivGraph* at this stage, please check
   to see whether you've installed all of the required dependencies.

4. To test your installation, you need to install the `pytest <https://docs.pytest.org/en/stable/index.html>`_ package.
Then from the cloned folder you can run the unit tests with the following command:
::

   $ pytest
.. _linksnodes:

===========================
Link and Node Dictionaries
===========================

After defining a network in *RivGraph* with :obj:`rivgraph.classes.rivnetwork.compute_network()`, two dictionaries will be created:
**links** and **nodes**.
This page of the documentation will describe the *keys* present in these dictionaries.

.. note::
   Not all dictionary key:value pairs are the same length as the number of links or nodes. Some dictionary key:value pairs contain meta-information about the network, or information that only applies to a subset of links/nodes. When the number of key:value pairs matches the number of links or nodes, then they are aligned such that the i'th index of any key:value pair refers to the same link or node.

- :ref:`links`
- :ref:`nodes`

.. _links:

--------------------
The Links Dictionary
--------------------

Links Key Values
----------------

**N** represents the number of links.

.. csv-table:: Generic Link Keys
   :file: links_generic.csv
   :header-rows: 1

.. csv-table:: Directionality-specific Keys
   :file: links_dir.csv
   :header-rows: 1

.. csv-table:: Braided River Exclusive Keys
   :file: links_river.csv
   :header-rows: 1

.. csv-table:: Delta Exclusive Keys
   :file: links_delta.csv
   :header-rows: 1

Accessing Link Values
---------------------

For example, if you know the *link_id* of the link you are interested in, you can get its index with :code:`links['id'].index(link_id)`.

.. _nodes:

--------------------
The Nodes Dictionary
--------------------

Nodes Key Values
----------------

**M** represents the number of nodes.

.. csv-table:: Generic Node Keys
   :file: nodes_generic.csv
   :header-rows: 1

Accessing Node Values
---------------------
For example, if you wish to find the links connected to :code:`node_id == 66`, you can use :code:`nodes['conn'][nodes['id'].index(66)]`.
.. _quickstart:

==========
Quickstart
==========

Install RivGraph
----------------

Follow one of the installation methods provided in :doc:`Installation Methods <../install/index>`.

Run Script
----------

.. todo:: create some simple test case/script to run
