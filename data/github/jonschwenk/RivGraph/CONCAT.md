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
