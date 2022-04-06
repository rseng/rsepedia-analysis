# bmm: Bayesian Map-Matching

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03651/status.svg)](https://doi.org/10.21105/joss.03651)

Map-matching using particle smoothing methods.

[Docs](https://bmm.readthedocs.io/en/latest/) and [methodology](https://arxiv.org/abs/2012.04602).

Contributing guidelines can be found in the repo's `CONTRIBUTING.md` file.

## Install
```
pip install bmm
```

## Load graph and convert to UTM
UTM (Universal Transverse Mercator) is a commonly used projection of spherical longitude-latitude
coordinates into square x-y coordinates.
```python
import numpy as np
import pandas as pd
import osmnx as ox
import json

import bmm

graph = ox.graph_from_place('Porto, Portugal')
graph = ox.project_graph(graph)
```

## Load polyline and convert to UTM
```python
data_path = 'simulations/porto/test_route.csv'
polyline_longlat = json.loads(pd.read_csv(data_path)['POLYLINE'][0])
polyline_utm = bmm.long_lat_to_utm(polyline_longlat, graph)
```
or generate fake data
```python
fake_route, fake_polyline_utm = bmm.sample_route(graph, timestamps=15, num_obs=25)
```

## Offline map-matching
```python
matched_particles = bmm.offline_map_match(graph, polyline=polyline_utm, n_samps=100, timestamps=15)
```

## Online map-matching
```python
# Initiate with first observation
matched_particles = bmm.initiate_particles(graph, first_observation=polyline_utm[0], n_samps=100)

# Update when new observation comes in
matched_particles = bmm.update_particles(graph, matched_particles, new_observation=polyline_utm[1], time_interval=15)
```

## Plot
```python
bmm.plot(graph, particles=matched_particles, polyline=polyline_utm)
```
![porto_mm](simulations/porto/test_route.png?raw=true "Map-matched route - Porto")




---
title: 'bmm: Bayesian Map-matching'
tags:
  - python
  - map-matching
  - GPS
authors:
  - name: Samuel Duffield
    affiliation: 1
affiliations:
 - name: University of Cambridge
   index: 1
date: 21 July 2021
bibliography: paper.bib
   
---

# Summary

`bmm` is a Python package providing probabilistic map-matching with uncertainty quantification.
Map-matching is the task of converting a polyline (series of noisy location observations - e.g. GPS data)
and a graph (collection of edges and nodes) into a continuous route trajectory restricted to the graph. Here a
continuous route is represented by series of connected edges as well as positions along said edges at
observation time. `bmm` uses Bayesian particle smoothing methods to produce a collection of particles, each of which
representing a continuous, plausible route along edges in the graph.

`bmm` is built on top of `osmnx` [@Boeing2017] - a python package assisting with the retrieval and processing
of OpenStreetMap data [@OpenStreetMap]. Although, `bmm` is applicable to be used on any suitably
labelled NetworkX graph [@Hagberg2008].

In addition, `bmm` utilises `numpy` [@Harris2020] and `numba` [@Lam2015] for fast scientific calculations,
`pandas` [@Reback2020] and `geopandas` [@Jordahl2020] for spatial data storage and manipulation
as well as `matplotlib` [@Hunter2007] for visualisation.

Documentation for `bmm` can be found at [bmm.readthedocs.io](https://bmm.readthedocs.io/en/latest/).


# Statement of need

Map-matching is a vital task for data driven inference involving GPS data.
Map-matching is often non-trivial, i.e. when the graph is dense, the observation noise is significant
and/or the time between observations is large. In these cases there may be multiple routes
that could have feasibly generated the observed polyline and returning a single trajectory is suboptimal.
Indeed, of 500 routes successfully map-matched using `bmm` from the Porto taxi dataset [@taxidata], 467 exhibited
multi-modality. This uncertainty over the inferred route would not be captured in the single trajectory
approach that is adopted by the most prominent map-matching software @Luxen2011 and @Yang2018, which adapt a Viterbi
algorithm - first applied to map-matching in @Newson2009. The code for @Luxen2011 is found as part of
the [OSRM project](https://github.com/Project-OSRM/osrm-backend) and represents an efficient C++ implementation
although is not easily accessible through Python. The software package accompanying @Yang2018 is found
at [fmm](https://github.com/cyang-kth/fmm) and provides extremely fast map-matching but without the convenience and
accessibility of working directly with an `osmnx` graph.

`bmm` adopts a state-space model approach as described in @Duffield2020
and produces a particle approximation that duly represents probabilistic
uncertainty in both the route taken and the positions at observation times. Additionally, `bmm` offers
support for both offline and online computation.


# Core Functionality

`bmm` can be used to convert a polyline (ordered series of GPS coordinates) into a collection of possible routes
along edges within a graph.

We assume that the graph is stored as a NetworkX [@Hagberg2008] object (which can easily be
achieved for a given region using `osmnx` [@Boeing2017]) and that the polyline is stored as an array or list of
two-dimensional coordinates in the same coordinate system as the graph. A common choice for coordinate system
is UTM (Universal Transverse Mercator) which as a square coordinate system (with unit metres) is less
cumbersome than the spherical longitude-latitude coordinates system (with unit degrees). `bmm` can convert
longitude-latitude to UTM using the `bmm.long_lat_to_utm` function.

### Offline Map-matching

Given a suitable graph and polyline `bmm` can be easily used to map-match
```python
matched_particles = bmm.offline_map_match(graph, polyline=polyline_utm,
                                          n_samps=100, timestamps=15)
```
Here the `n_samps` parameter represents the number of particles/trajectories to output and `timestamps` is the
number of seconds between polyline observations - this can be a float if all observation times are equally spaced,
an array of length one less than that of the polyline representing the unequal times between observations or an 
array of length equal to the polyline representing UNIX timestamps for the observation times.

The output of `bmm.offline_map_match` is a `bmm.MMParticles` object that contains a `particles` attributes listing
the possible trajectories the algorithm has managed to fit to the polyline - full details can be found at
[bmm.readthedocs.io](https://bmm.readthedocs.io/en/latest/).

### Online Map-matching

`bmm` can also map-match data that arrives in an online or sequential manner. Initiate a `bmm.MMParticles`
with the first observation
```python
matched_particles = bmm.initiate_particles(graph,
                                           first_observation=polyline_utm[0],
                                           n_samps=100)
```
and then update as new data comes in
```python
matched_particles = bmm.update_particles(graph,
                                         matched_particles,
                                         new_observation=polyline_utm[1],
                                         time_interval=15)
```

### Parameter Tuning

The statistical model described in @Duffield2020 has various parameters which can be adjusted to fit the features
of the graph and time interval setup. This can be done by adjusting a `bmm.MapMatchingModel` argument or its
default `bmm.ExponetialMapMatchingModel` which is taken as an optional `mm_model` argument in the above
map-matching functions. In addition, these parameters can be learnt from a series of polylines using `bmm.offline_em`. 


### Plotting

Once a polyline has been succesfully map-matched, it can be visualised using `bmm`
```python
bmm.plot(graph, particles=matched_particles, polyline=polyline_utm)
```
![](simulations/porto/test_route.png)



# Acknowledgements

Samuel Duffield acknowledges support from the EPSRC.


# References# Contributing to `bmm`
`bmm` would love contributions! Whether that's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## We Develop with GitHub
`bmm` uses GitHub to host code, track issues and feature requests, as well as accepting pull requests.

## We Use [GitHub Flow](https://guides.github.com/introduction/flow/index.html), So All Code Changes Happen Through Pull Requests
Pull requests are the best way to propose changes to the codebase (we use [GitHub Flow](https://guides.github.com/introduction/flow/index.html)). We actively welcome your pull requests:

1. Fork the repo and create your branch from `master`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. Issue that pull request!

## Any contributions you make will be under the MIT Software License
In short, when you submit code changes, your submissions will be understood to be under the same [MIT License](https://github.com/SamDuffield/bmm/blob/master/LICENSE.txt) that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using GitHub's [issues](https://github.com/SamDuffield/bmm/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](https://github.com/SamDuffield/bmm/issues/new); it's that easy!

## Write bug reports with detail, background, and sample code
**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

People *love* thorough bug reports. I'm not even kidding.

## License
By contributing, you agree that your contributions will be licensed under its [MIT License](https://github.com/SamDuffield/bmm/blob/master/LICENSE.txt).

## References
This document was adapted from [briandk](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62) and [Facebook's Draft](https://github.com/facebook/draft-js/blob/master/CONTRIBUTING.md).
Index
=====
Welcome to ``bmm``'s documentation!
====================================

``bmm`` provides map-matching with uncertainty quantification for both online and offline inference!

Map-matching converts a series of noisy GPS coordinates into a continuous trajectory that is restricted to a graph (i.e. road network) or in the case of ``bmm`` a collection of continuous trajectories representing multiple plausible routes!

``bmm`` is built on top of ``osmnx``, an `awesome package for retrieving and processing OpenStreetMap data <https://github.com/gboeing/osmnx>`_.

The probabilistic model and particle smoothing methodology behind ``bmm`` can be found on `arXiv <https://arxiv.org/abs/2012.04602>`_
and the source code on `GitHub <https://github.com/SamDuffield/bmm>`_.


Docs
=====================

.. toctree::
    :maxdepth: 3

    functions
    classes
    genindex


Install
=========
``pip install bmm``


Quickstart
===========
Load graph and convert to UTM (Universal Transverse Mercator), a commonly used projection of spherical longtitude-latitude
coordinates into square x-y coordinates::

    import numpy as np
    import pandas as pd
    import osmnx as ox
    import json
    import bmm

    graph = ox.graph_from_place('Porto, Portugal')
    graph = ox.project_graph(graph)

Beware that downloading graphs using ``osmnx`` can take a few minutes, especially for large cities.

Load polyline and convert to UTM::

    data_path = 'simulations/porto/test_route.csv'
    polyline_longlat = json.loads(pd.read_csv(data_path)['POLYLINE'][0])
    polyline_utm = bmm.long_lat_to_utm(polyline_longlat, graph)

Offline map-matching
^^^^^^^^^^^^^^^^^^^^^
::

    matched_particles = bmm.offline_map_match(graph, polyline=polyline_utm, n_samps=100, timestamps=15)

Online map-matching
^^^^^^^^^^^^^^^^^^^^
Initiate with first observation::

    matched_particles = bmm.initiate_particles(graph, first_observation=polyline_utm[0], n_samps=100)

Update when new observation comes in ::

    matched_particles = bmm.update_particles(graph, matched_particles, new_observation=polyline_utm[1], time_interval=15)



Sanity Check
=============
You can manually test that ``bmm`` is working sensibly for a given graph by generating synthetic data::

    graph = ox.graph_from_place('London, UK')
    graph = ox.project_graph(graph)
    generated_route, generated_polyline = bmm.sample_route(graph, timestamps=15, num_obs=20)

Note that the London graph takes some time (~10mins) to download and for testing on synthetic data it may be worth considering a smaller region
(although not so small that the ``sample_route`` function consistently terminates early due to reaching the edge of the graph).

Run map-matching on the generated polyline::

    matched_particles = bmm.offline_map_match(graph, generated_polyline, n_samps=100, timestamps=15)

Plot true generated route::

    bmm.plot(graph, generated_route, generated_polyline, particles_color='green')

.. image:: sanity_check_truth.png

Plot map-matched particles::

    bmm.plot(graph, matched_particles, generated_polyline)

.. image:: sanity_check_mm.png

Classes
========

.. autoclass:: bmm.MMParticles
    :members:

.. autoclass:: bmm.MapMatchingModel
    :members:

.. autoclass:: bmm.ExponentialMapMatchingModel
    :members:

Functions
==========

.. autofunction:: bmm.offline_map_match

.. autofunction:: bmm.initiate_particles
.. autofunction:: bmm.update_particles
.. autofunction:: bmm._offline_map_match_fl

.. autofunction:: bmm.sample_route
.. autofunction:: bmm.random_positions

.. autofunction:: bmm.offline_em

.. autofunction:: bmm.plot

.. autofunction:: bmm.get_possible_routes

.. autofunction:: bmm.cartesianise_path
.. autofunction:: bmm.get_geometry
.. autofunction:: bmm.discretise_edge
.. autofunction:: bmm.observation_time_indices
.. autofunction:: bmm.observation_time_rows
.. autofunction:: bmm.long_lat_to_utm


