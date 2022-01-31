<img src="https://user-images.githubusercontent.com/520137/94442680-873b3200-01a4-11eb-8ec0-cb7ce5a56cc9.png" alt="MorphoPy logo" width="450"/>

# A Python package for feature extraction of neural morphologies
(c) 2020 by Sophie Laturnus, Adam von Daranyi, Ziwei Huang and Philipp Berens.


********

## <a name="content">Content</a> 
- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Uninstallation](#uninstallation)
- [Usage](#usage)
- [Module description](#module)

## <a name="overview">Overview</a> 

MorphoPy is a Python3 package that uses networkX to compute and show information about neurites.
The input can be passed with single swc-files or it can handle whole directories with multiple files at once.
You can use MorphoPy imported in Python, as well as from the command line as a batch-tool.

The current working build:
 
 **version 0.7.1**
 
All builds are tested on Linux (Debian and CentOS) and Windows 10.

You can find all working builds at [pypi](https://pypi.org/project/morphopy/).

And a full documentation of the source code is available at [readthedocs](https://morphopy.readthedocs.io/en/latest/).

[back to start](#content)

## <a name="requirements">Software Requirements</a>

In the following, all commands written in boxes need to be entered into your terminal.

### Mac:

 - homebrew (to install latest version of python3)
 - python >3.4: without homebrew you can find python [here](https://www.python.org/downloads/mac-osx/)
 - pip
 
**1)** If you want the latest python3 version you need to install homebrew,
       for that just run this command in your terminal:

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    
**2)** Then you can easily install the newest version of python3:

	brew install python

**3)** Now you should have python3 and pip3 installed. You can check the version like this:

	python3 --version
	pip3 -V

If both are enough you can continue with installing MorphoPy, see below.

### Linux:

 - python >3.4
 - pip 3 

Install Python3 on **Debian/Ubuntu** systems:

	apt-get install python3 python3-pip

Install Python3 on **RedHat/CentOS/Fedora systems** (perhaps python version differs):

	yum install python36 python3-pip

That's it. All requirements are met and you can continue with the MorphoPy installation. See below. 

### Windows:

 - python >3.5 (x64): download from [here](https://www.python.org/downloads/windows/)
 - pip : be sure you selected it during installation with the base python package
 - Microsoft Build Tools >14.0 you can download them [here](https://visualstudio.microsoft.com/visual-cpp-build-tools/) 

**1)** Install python with pip by executing the downloaded installation file and
       be sure to check the option to add python paths to enviroment variable at the first step!
       
**2)** Install Microsoft Build Tools with the automatic installation tool from the website above.
       Select C/C++ Compiler Tools at the selection page

All requirements are met now, see below for continue with MorphoPy installation.

[back to start](#content)

## <a name="installation">Installation:</a>

###  Install the MorphoPy package (all platforms):

Install the python package with pip3 and get the latest build:

	pip3 install morphopy

## <a name="uninstallation">Uninstallation:</a>

You can simply uninstall the package with pip3:

	pip3 uninstall morphopy

	
[back to start](#content)

## <a name="usage">Usage</a> 
After the installation you can simply call:

	morphopy <options>
Help:

	Usage: morphopy -c <compute_feature> -i <swc_file>|<directory> [-o <output directory>]
	                   [--long | --func <persistence_function> | --conf <config_file>] [--help]
	
	Options:
    -h, --help                  show this help.
    
	-c, --compute               parameter for selecting the computing feature:
	                            persistence: compute persistence data
                                stats      : compute morphometric statistics
                                density    : create density maps
       statistics options:
       --long                   you can change your output format, in long
                                format you get all values in a separate row.
                                (default: all values in one row)   
       persistence options:
       --func                   if persistence is selected as feature, you can
                                specify with this option a method function you
                                want to use at computing the persistence.
                                (default: radial distance function)
       density map options:
       --conf                   if density map is selected, you can pass a
                                config file with more parameters for creating
                                the density maps. (optional)
    -i, --input                 specifies a swc-file or a directory as input
                                for morphopy (mandatory)
    -o, --output                specifies the output directory for saving the
                                results in. (default: same as source)


Available functions for persistence at the moment are:
 - radial_distance (default function)
 - height
 - path_length
 - branch_order

A sample config file for density maps looks like this (stored in a text file):

	[global]
	# specific distance for resampling nodes:
    distance: 1
    # width of each bin in microns across all dimensions
    #bin_size: 20
    # number of bins for each dimension (only if you don't use bin_size)
    n_bins_x: 20
    n_bins_y: 20
    n_bins_z: 20
    # if true: probabilty density is returned, count histogram otherwise
    density: True
    # smoothing the density data
    smooth: True
    # sigma used at smoothing
    sigma: 2
    # normalization bounds for density map:
    [norm_bound]
    r_max_x: 238.85
    r_max_y: 140.95
    r_max_z: 285.97
    r_min_x: -236.17
    r_min_y: -24.2
    r_min_z: -173.72

[back to start](#content)

## <a name="contributing">Contributing to MorphoPy </a>

We tested MorphoPy to the best of our knowledge and abilities in the scope of several projects. If you still find a bug
or you are missing a feature, please do not hesitate to contact us via [GitHub issues](https://github.com/berenslab/MorphoPy/issues).
Please try to provide a minimal example that reproduces the bug you want to be fixed.
If you want to develop the code base further, you can work with git pull requests. Please make sure that you document
the code and add tests and examples of how to use your code.


## <a name="citation"> Citing MorphoPy </a>

If you are using MorphoPy for your research or your work project please make sure to cite us and this repository:
```
@article{Laturnus2020,
  doi = {10.21105/joss.02339},
  url = {https://doi.org/10.21105/joss.02339},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2339},
  author = {Sophie Laturnus and Adam von Daranyi and Ziwei Huang and Philipp Berens},
  title = {MorphoPy: A python package for feature extraction of neural morphologies.},
  journal = {Journal of Open Source Software}
}
```
or
```
@misc{morphopy,
  author = {Laturnus, Sophie and von Daranyi, Adam and Huang, Ziwei and Berens, Philipp},
  title = {MorphoPy},
  publisher = {GitHub},
  journal = {GitHub repository},
  url = {https://github.com/berenslab/MorphoPy}
}
```
## <a name="module">Module description</a> 

**Important:** MorphoPy requires the soma to be one single point. If more than 3 somatic points are present in the
reconstruction file they will be automatically collapsed to the centroid of their convex hull on file loading. If the 
soma is described by 2 to 3 points they will be automatically collapsed to their mean (also see utils.get_standardized_swc).

MorphoPy currently only supports neurites that connect back to the soma. This means, axons that emerge from dendritic
structures can not be handled.

A neuron is represented as a directed acyclic graph with node attributes _id, x-, y-, z- position, radius_ and _type_id_ (soma: 1, axon: 2, dendrite: 3, apical dendrite: 4), and with edge attributes _path_length_ and _euclidean_dist_. Positions, radius and length mesaures are assumed to be given in microns. 

![Node and edge attributes](https://user-images.githubusercontent.com/520137/80974465-0d836000-8e21-11ea-87b3-0bc9fdb41a3f.png)

*Fig. 1: Node and edge attributes associated with each neuron graph.*

All data is stored in the [tidy data format](http://vita.had.co.nz/papers/tidy-data.pdf).

Please also refer to our [tutorial](https://github.com/berenslab/MorphoPy/blob/master/notebooks/MORPHOPY%20Tutorial.ipynb)
and the [documentation](https://morphopy.readthedocs.io/en/latest/).


### Density maps
Density maps are marginal histograms over the neural mass. MorphoPy allows you to create density maps of different projections through the function compute_density_maps(). Per default it computes x, y, z, xy, xz and yz density maps from the point cloud of the original reconstruction. The point cloud is constructed through resampling along all neurites with a default distance of 1 micron. The resulting point cloud is then binned into bins of 20 microns and smoothed using Gaussian smoothing with std of 1.

However, you can customize all these parameters by passing a config file to the function (see [above](#usage)).

### Morphometric statistics

`MorphoPy` offers a default selection of 28 single-valued morphometric statistics, namely:
- number of branch points
- width (x-extent), depth (y-extent), height (z-extent)
- number of tips
- number of neurites extending from the soma directly (stems)
- the total path length (in microns)
- average and maximal radius thickness (with the soma excluded)
- total surface and volume
- maximal path distance to the soma
- maximal branch order
- maximal, min and median path angle
- average soma exit angle
- maximal path length of a segment
- median intermediate and median terminal segment length
- log of max, min and median tortuosity across all edges (= path length/euclidean length)
- max, min and average branch angle
- maximal branching degree (with soma excluded)
- _weighted proportional sum of absolute deviations_ as a measure of tree asymmetry (for more details see https://www.sciencedirect.com/science/article/pii/0165027086901196)

![Morphometric statistics that can be queried.](https://user-images.githubusercontent.com/520137/80974473-0f4d2380-8e21-11ea-8ce2-acb8153cece4.png)

*Fig. 2: Explanatory schematic of the morphometric statistics that can be computed on all nodes. Left: distance measures, Right: angles.*

### Morphometric distributions

Morphometric distributions are not (yet) available via the command line tool.
Frequency histograms or Gaussian kernel density estimates can be queried for all listed key statistics using the
methods `get_histogram(key)` or `get_kde_distribution(key)`. If you provide a distance measure (e.g. branch order,
path distance from soma or radial distance) the returned distribution will be two-dimensional and allows to investigate
a features' development across space.
Additionally, it is possible to compute [Sholl intersection profiles](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1244622/) using the function `get_sholl_intersection_profile()`. 

Key statistics are
- branch orders
- Strahler order
- branch angles
- path angles
- root angles
- thickness
- segment lengths
- path length to soma
- radial distance

### Persistence

Persistence diagrams are a concept from topology. They have been introduced as descriptors of neural morphologies by [Kanari et al.](https://link.springer.com/article/10.1007/s12021-017-9341-1) and [Li et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182184) 

<img src="https://user-images.githubusercontent.com/520137/80973456-b0d37580-8e1f-11ea-92f4-2dfa5d9729d9.png" alt="Schematic of persistence diagrams" width="450"/>

*Fig. 3: Schematic of how a persistence diagram is generated. The longest branch with the longest 'lifetime' is marked in red. Taken from Kanari et al. 2018.*


The recorded birth and death times in the figure above are based on a certain distance (or lifetime) function. `MorphoPy` implements four different distance functions to choose from: radial distance (default), height, path length and branch order. 
They all compute the distance of a point with respect to the soma. In the command line tool you can switch between them using the `--func` keyword (see [above](#usage)). 
To provide your own distance function, add its code and its keyword to the `persistence_functions.py` file, but make sure 
that the distance functions interface fits the specification `custom_distance(networkx.DiGraph,node_id_end, node_id_start)` 
(see Fig. 4 and the [tutorial](https://github.com/berenslab/MorphoPy/blob/master/notebooks/MORPHOPY%20Tutorial.ipynb) for an example). 

<img src="https://user-images.githubusercontent.com/520137/80983512-eaf74400-8e2c-11ea-94cc-040275f6aeca.png" alt="How to add a custom distance function" width="751"/>

*Fig. 4: How to add a custom persistence distance function. To be able to call it from the command line you need to add it to the functions list.*

If you are using the API you can simply pass a function to the `get_persistence()`-method (see the [tutorial](https://github.com/berenslab/MorphoPy/blob/master/notebooks/MORPHOPY%20Tutorial.ipynb) for an example).


### Not enough? ###

You want to compute your own features? Go for it! We recommend you to check out `networkx` and `shapely` for more options.


[back to start](#content)

---
title: 'MorphoPy: A python package for feature extraction of neural morphologies.'
tags:
    - Python
    - neurons
    - morphologies
    - anatomy
authors:
    - name: Sophie Laturnus
      orcid: 0000-0001-9532-788X
      affiliation: "1,3,4"
    - name: Adam von Daranyi
      affiliation: 4
    - name: Ziwei Huang
      affiliation: "1,3,4"
    - name: Philipp Berens
      orcid: 0000-0002-0199-4727
      affiliation: "1,2,3,4"
affiliations:
    - name: Institute for Ophthalmic Research, University of Tübingen, Germany
      index: 1
    - name: Institute for Bioinformatics and Medical Informatics, University of Tübingen, Germany
      index: 2
    - name: Bernstein Center for Computational Neuroscience, University of Tübingen, Germany
      index: 3
    - name: Center for Integrative Neuroscience, University of Tübingen, Germany
      index: 4
date: 26 May 2020
bibliography: morphopy.bib
---

# Summary

For a long time, the anatomy of a neuron has been considered a defining feature of neural cell types.
However, computational analysis of neuronal morphologies persists to be a challenging problem. It usually begins with choosing a
feature representation in order to make individual morphologies amenable to statistical data analysis or to processing by
machine learning methods. Over the years, many different feature representations have been suggested in the literature,
such as density maps [@jefferis:2007], single valued summary statistics (morphometrics) [@scorcioni:2008; @neurom] or,
more recently, persistence images [@li:2017; @kanari:2018].
Unfortunately, current software packages for extracting them from morphological data are often focused on solely one such
representation and implemented across various programming languages.

Our software package `MorphoPy` provides straightforward access to different feature representations from neural morphologies for
downstream statistical analysis. It bundles common representations such as density maps, morphometrics, morphometric distributions
and persistence images in one simple open source framework implemented in Python to make them accessible to a larger community.
`MorphoPy` can be used either as a standalone command line tool or as a package within an interactive computing workflow.

![Neural reconstructions are represented as direct acyclic graphs with node and edge attributes.\label{fig:attributes}](./figures/Fig1_attributes_small.png)

`MorphoPy` builds on the functionality of the networkx package [@hagberg:2008] and represents each neuron as a directed
acyclic tree graph with node and edge attributes \autoref{fig:attributes}. The package supports to read in files in the common
[swc-format](http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html), offers functions to
compute various feature representations and provides 2D plotting routines for data exploration (\autoref{fig:plot}).
```python

from morphopy.neurontree import NeuronTree as nt
from morphopy.computation import file_manager as fm

N = fm.load_swc_file("../data/EC3-80604.CNG.swc")
Dendrites = N.get_dendritic_tree()

from morphopy.neurontree.plotting import show_threeview
fig = plt.figure(figsize=(10,10))
show_threeview(N, fig)
```

![Plotting reconstructions in 2D. \label{fig:plot}](./figures/threeview_dendrites.png)

As shown in the code snippet above, it is also possible to split the reconstruction into its different parts (axon or dendrites only)
and operate on each neurite type separately. 

Currently, MorphoPy supports the following feature representations:

_Density maps_ are computed on the basis of a configuration file (or a dictionary) that controls parameters such as bin size
and binning ranges. Additionally, users can specify whether and to which degree
they want to smooth each density map (see \autoref{fig:dms}).

![XY-density map of the dendrite plotted above with different degrees of Gaussian smoothing. \label{fig:dms}](./figures/density_map_smoothing.png)

A variety of _morphometric statistics_ can be computed on the nodes and edges of each reconstruction.
The `get_morphometric_statistics()`-method offers a precompiled single valued selection of these statistics including e.g.
min/max branch angles, maximal branch order, and maximal path length to the soma (see \autoref{fig:morphometrics}), but in principle, they can be adjusted to the user's personal preference.

![Node and edge related morphometric statistics. \label{fig:morphometrics}](./figures/fig_morphometrics.png)

Additionally, it is possible to query the entire _morphometric distribution_ of each statistic either in form of a histogram or as a
Gaussian kernel density estimate (kde). \autoref{fig:morphdist}, for example, shows the kde of radial distances, branch angles and their
combination for the dendrites shown in \autoref{fig:plot}.

![Kernel density estimate of branch angles as a function of the radial distance from the soma. \label{fig:morphdist}](./figures/2D_morph_dist.png)

Furthermore, `MorphoPy` supports the generation of _2D persistence diagrams_. Persistence diagrams describe the branching
of neural morphologies [@li:2017; @kanari:2018] with respect to a specified distance function. By default, `MorphoPy`
computes a persistence diagram based on the radial distance from the soma, but users can choose from four different
pre-implemented distance functions (radial distance, path length, height or branch order) or provide their own.

```python
from morphopy.computation.feature_presentation import get_persistence

import numpy as np
def custom_distance(G, u, v):
    """
    Returns a distance between nodes u and v, 
    which both are part of the graph given in G.
    """
    n = G.node[u]['pos']
    r = G.node[v]['pos']
    return np.dot(n, r)

df = get_persistence(Dendrites.get_topological_minor())
df_custom = get_persistence(Dendrites.get_topological_minor(), f=custom_distance)
```

In addition to working as a package in interactive computing environments, `MorphoPy` can be called from the command line
to operate on single files or entire batches.
```bash
MorphoPy.py -c [density|persistence|stats] 
            -i ['path_to_file'|'path_to_folder']
```
For a full documentation of `MorphoPy`'s functionality please refer to our [documentation](https://morphopy.readthedocs.io/en/latest/) and tutorial on our [GitHub page](https://github.com/berenslab/MorphoPy).

`MorphoPy` has been developed in the context of a benchmarking study for cortical interneuron cell type classification
based on their morphology [@laturnus:2019]. It has already been used in a series of scientific publications that relate
transcriptome, electrophysiology and morphology of cortical interneurons in V1 and M1 [@scala:2019; @scala:2020].

# Acknowledgements

This work was funded by the German Ministry of Education and Research(FKZ 01GQ1601), the German Research Foundation (DFG)
under Germany’s Excellence Strategy (EXC2064/1 – 390727645; BE5601/4-1, SFB 1233 “Robust Vision”, Project number 276693517),
and the National Institute of Mental Health (U19MH114830).

# References
Computation
-----------
.. toctree::
   :maxdepth: 2

   feature_presentation.rst
   file_manager.rst
   persistence_functions.rstCommand line script
-------------------

.. toctree::
   :maxdepth: 2

   morphopy.rstUtils
=====

.. automodule:: morphopy.neurontree.utils
    :members:
Feature Presentation
====================

.. automodule:: morphopy.computation.feature_presentation
    :members:
NeuronTree
----------
.. toctree::
   :maxdepth: 2

   neurontree.rst
   plotting.rst
   utils.rst
NeuronTree
==========

.. autoclass:: morphopy.neurontree.NeuronTree.NeuronTree
    :members:
    
    .. automethod:: __init__
Welcome to Morphopy's documentation!
====================================

MorphoPy is a Python3 package that uses networkX to compute and show information
about neurites. The input can be passed with single swc-files or it can handle
whole directories with multiple files at once. You can use MorphoPy imported in
Python or from command line as batch-tool as well.

You can find MorphoPy with more information `here... <https://github.com/berenslab/MorphoPy>`_

Contents:
=========

.. toctree::
   :maxdepth: 2

   command.rst
   computation.rst
   neuron.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Morphopy
==========

.. automodule:: morphopy.MorphoPy
    :members:
File Manager
============

.. automodule:: morphopy.computation.file_manager
    :members:
Plotting
========

.. automodule:: morphopy.neurontree.plotting
    :members:
Persistence Functions
=====================

.. automodule:: morphopy.computation.persistence_functions
    :members:
