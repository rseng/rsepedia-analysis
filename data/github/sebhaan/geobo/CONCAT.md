![](https://github.com/sebhaan/geobo/blob/master/docs/Logo_GeoBO.png?raw=True)

# GeoBO: A Python package for Multi-Objective Bayesian Optimisation and Joint Inversion in Geosciences

``GeoBO`` is build upon a probabilistic framework using Gaussian Process (GP) priors to jointly solve multi-linear forward models. This software generates multi-output 3D cubes of geophysical properties (e.g. density, magnetic susceptibility, mineral concentrations) and their uncertainties from 2D survey data (e.g. magnetics and gravity) and any pre-existing drillcore measurements. The reconstructed 3D model is then used to query the next most promising measurement location given an expensive cost function (e.g. for drillcores). A ranked list of new measurements is proposed based on user-defined objectives as defined in the acquisition function which typically aims to optimize exploration (reducing global model uncertainty) and exploitation (focusing on highly promising regions) while minimizing costs. 

![GeoBO Framework](https://github.com/sebhaan/geobo/blob/master/docs/Overview_illustration.png?raw=True)

## Table of Contents
- [Definitions](#definitions)
- [Functionality](#functionality)
- [Installation And Requirements](#installation-and-requirements)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Documentation](#documentation)
- [Usage and Settings](#usage-and-settings) 
- [Examples and Tests](#examples)
  - [Synthetic Models](#synthetic-models)
  - [Drillcore Test Example](#drillcore-test-example)
  - [Results and Output Files](#results-and-output-files)
- [Options and Customization](#options-and-customization)
  - [Custom Linear Forward Models](#custom-linear-forward-models)
  - [Gaussian Process Kernel Functions](#gaussian-process-kernel-functions)
- [Literature](#literature)
- [Related Software](#related-software)
- [Attribution and Acknowledgments](#attribution-and-acknowledgements)
  - [Project Contributors](#project-contributors)
- [License](#license)


## Definitions

Bayesian Optimisation (BO) is a powerful framework for finding the extrema of objective functions that are noisy, expensive
to evaluate, do not have a closed-form (e.g. black-boxfunctions), or have no accessible derivatives. The model used for approximating the objective function is called surrogate model, which is typically based on a [Gaussian Process models](https://en.wikipedia.org/wiki/Gaussian_process) for tractability. Gaussian Processes define a prior over functions (typically given by a kernel function) and is used to propose points in the search space where sampling is likely to yield an improvement. The specific set of objectives for the improvement are defined in an acquisition function, which guides the search for a user-defined optimum. An example use case scenario is described in a nutshell in [OPTIMIZATION_FOR_ACTIVE_SENSORFUSION_IN_A_NUTSHELL.pdf](https://github.com/sebhaan/geobo/blob/master/docs/OPTIMIZATION_FOR_ACTIVE_SENSORFUSION_IN_A_NUTSHELL.pdf).

### Acquisition function 

The key of BO is the acquisition function, which typically has to balance between: 

a) exploration, i.e., querying points that maximise the information gain and minimize the uncertainty of a model
b) exploitation, i.e. querying points that maximise the reward (e.g. concentrating search in the
vicinity locations with high value such as minerals)
c) minimize the number of samples given an expensive cost function for any new measurement.


### Forward Models and Joint Inversion
In geology and geophysics, inversion problems occur whenever the goal is to reconstruct the geological conditions, i.e. the 3D distribution of physical rock properties, that give rise to a set of (2D) geophysical observations. Since the number of possible geological configurations is typically greater than the number of observational constraints, the problem is nearly always under-determined.
Forward models transform the localized measurement of a remote sensor grid into a 3D representation of geophysical properties of a region. The most common geophysical linear forward model are gravity and magnetic forward models, which are computed using Li’s tractable approximation. Joint inversion is  simultaneously interpreting  multiple (distinct) sensor measurements using a single model to provide a better constrained joint solution rather than taking individual solutions that only satisfy their aspect of data on their own. 



## Functionality

GeoBO's probabilistic framework includes all steps from  prior selection, data fusion and inversion, to sensor optimisation and real world model output. The main functionalities of GeoBO are summarised in the following:

 - Joint probabilistic inversion tool by solving simultaneously multi-linear forward models (e.g. gravity, magnetics) using cross-variances between geophysical properties (cross-variance terms can be specified by user). 
 - Output 1: Generation of cubes and computation of complete posterior distribution for all geophysical properties (described by their mean and variance value at each location (cubecell aka voxel). 
 - Output 2: Generation of ranked proposal list for new most promising drillcores based on global optimisation of acquisition function
 - Templates for acquisition function to use in Bayesian Optimisation
 - Flexible parameter settings for exploration-exploitation trade-off and inclusion of local 3D cost function in acquisition function 


Other features are:

 - Generation of simulated geophysical data with a choice of three different models
 - Package includes geological survey/drillcore sample as well as synthetic data and functions for synthetic data generation
 - Generation of 2D/3D visualisation plots of reconstructed cubes and survey data
 - 3D Cube export in VTK format (for subsequent analysis, e.g., in Python or with ParaView)
 - Options to include any pre-existing drillcore data 
 - Included linear forward models: density-to-gravity and magnetic susceptibility-to-magnetic field; custom linear forward models can be added (see [Custom Linear Forward Models](#custom-linear-forward-models))
 - Library of Gaussian Process (GP) kernels including sparse GP kernels
 - Flexible settings for any cube geometry and resolution
 - (Optional) Optimization of GP hyperparameters and cross-correlation coefficients via computation of marginal GP likelihood

Example outputs can be found in the directory `examples/results/`.

## Installation And Requirements


### Installation

To install GeoBO locally using setuptools: 

```sh
python setup.py build
python setup.py install
```

or using pip:

```sh
pip3 install geobo
```

The installation can be tested by running the example with included synthetic data and default settings:

```sh
cd geobo/
python main.py tests/settings_example1.yaml
```

### Requirements

- python >=3.6
- numpy
- matplotlib
- scikit_image
- scipy
- rasterio
- pandas
- pyvista
- skimage
- PyYAML


### Documentation 

Documentation conversion is generated using pandoc. The README markdown file can be converted to PDF:

```bash
pandoc -V geometry:margin=1.0in README.md -o README.pdf
```

A complete API documentation for all modules can be found here:

- [`run_geobo.py`](http://htmlpreview.github.io/?https://github.com/sebhaan/geobo/blob/release/docs/APIdocs/geobo/run_geobo.html)
- [`inversion.py`](http://htmlpreview.github.io/?https://github.com/sebhaan/geobo/blob/release/docs/APIdocs/geobo/inversion.html)
- [`kernels.py`](http://htmlpreview.github.io/?https://github.com/sebhaan/geobo/blob/release/docs/APIdocs/geobo/kernels.html)
- [`cubeshow.py`](http://htmlpreview.github.io/?https://github.com/sebhaan/geobo/blob/release/docs/APIdocs/geobo/cubeshow.html)
- [`sensormodel.py`](http://htmlpreview.github.io/?https://github.com/sebhaan/geobo/blob/release/docs/APIdocs/geobo/sensormodel.html)
- [`simcube.py`](http://htmlpreview.github.io/?https://github.com/sebhaan/geobo/blob/release/docs/APIdocs/geobo/simcube.html)
- [`utils.py`](http://htmlpreview.github.io/?https://github.com/sebhaan/geobo/blob/release/docs/APIdocs/geobo/utils.html)


## Usage and Settings

1) Change the main settings such as filenames and parameters in `settings.yaml`. These settings specify:

- directory, filenames, and geophysical drillcore properties
- the generated cube's geometry, size, and resolution
- Gaussian Process settings (lengthscale, input data uncertainty, correlation coefficients, kernel function)
- local Earth's magnetic field vector
- Bayesian Optimisation Settings (vertical/non-vertical drillcores, the exploration/exploitation and cost weighting)
- plotting settings
- optional generation of simulated data 

2) Then run geobo
```sh
cd geobo/
python main.py settings.yaml
```

The main functions for the acquisition function can be found in [`run_geobo.py`](docs/APIdocs/geobo/run_geobo.html); visualisation functions and VTK export are defined in [`cubeshow.py`](docs/APIdocs/geobo/cubeshow.html); inversion functions are defined in [`inversion.py`](docs/APIdocs/geobo/inversion.html). 


## Examples and Tests


### Synthetic Models

Synthetic geophysical models can be created by setting switching on `gen_simulation` in the settings yaml file. 
Three different models are so far implemented:

- two-cylindric dipping bodies (`modelname: 'cylinders'`) 
- two layer model (`modelname: 'layers2'`)
- three layer model (`modelname: 'layers3'`)
For each model a 3D voxel cube with geological structures is generated with density and magnetic susceptibility properties, plus
the corresponding 2D gravity and magnetic remote sensor measurements. Other custom models can be included by adding a new model in function `create_syncube()` in `simcube.py`

Result examples of the synthetic models are stored in the subfolder `examples/testdata/synthetic/`.
An example settings file is given in `settings_example1.yaml` and can be run by

```sh
cd geobo/
python main.py tests/settings_example1.yaml 
```


### Drillcore Test Example

Another examples includes drillcore and gravity/magnetic survey data (`examples/testdata/sample/`). This example can be run with

```sh
cd geobo/
python main.py tests/settings_example2.yaml
```
and creates the reconstructed density and magnetic susceptibility cubes, uncertainty cubes


### Results and Output Files

The output results include the generated reconstructed density and magnetic susceptibility cubes and their corresponding uncertainty cubes, visualisations of original survey data and reconstructed properties, and list of new measurement proposals.

List of Joint Inversion Result Outputs: 

- 3D Cube files in vtk format (to use, e.g., with PyVista or ParaView): Output of cross-correlated reconstructed properties (density: **cube_density.vtk**, magnetic susceptibility: **cube_magsus.vtk**, property from drill: **cube_drill.vtk**) and their uncertainity cubes in terms of standard deviation (**cube_density_variance.vtk, cube_drill_variance.vtk, cube_drill_variance.vtk**). In case the cross-covariance terms (`gp_coeff` in the settings file) are all set to zero, the solutions for each forward model are independent from each other.
- Optional (Default optiion: plot=True in function `read_surveydata()`): Images of the input gravitational field (**gravfield.png**) and magnetic field (**magfield.png**) and their corresponding downsampled images (**gravfield_downsampled.png, magfield_downsampled.png**) 
- Optional (if `plot3d`:True in settings): 3D Preview of reconstructed properties: density (**density-mesh3D.png**), magnetic suscpetibility (**magsus-mesh3D.png**), and drill property (**drill-mesh3D.png**). However, it is recommended to use PyVista or Paraview instead. Future updates for 3D visualisations will be based on PyVista. 
- Optional (if `plot_vertical`:True in settings): 2D maps of vertically (along z-axis) mean value of cube properties (**dens_rec2D_loc2.png, magsus_rec2D_loc2.png, drill_rec2D_loc2.png**)

List of Bayesian Optimisation Output:

- List of all new measurement proposals (here for drillcores) ranked from maximum (hightest gain) to minimum of optimisation function. The results are saved as csv file (**newdrill_proposals_non-vertical.csv** or **newdrill_proposals_vertical.csv**) and include for vertical drillcores the position (Easting, Northing) and for non-vertical drillcores position and drill angles (Azimuth, Dip).
- Points of proposed measurement positions on top of reconstructed drill property image (mean projection along z-axis): The output figure (non-vertical drillcores:**newdrill_proposals.png** or vertical drillcores:**newdrill_vertical_proposals.png**) shows the location of the already existing drills as given by input measurements (black points), the new proposed drill positions (white), and the best (maximum of optimsation function) new drill location (red). 

![Example image of new measurement proposals (black_ existing, white: new proposed, red: best proposal) on top of reconstructed property (mean value projection)](https://github.com/sebhaan/geobo/blob/master/docs/newdrill_proposals.png?raw=True).

## Options and Customization


### Custom Linear Forward Models

The relationship between a physical system (or its model parameters) ***P*** and the observed sensor data **y** is described by a linear forward model

**y** = **G** ***P***,  

where **G** is the transformation operator or matrix. The gravitational and magnetic forward model can be determined analytically by using Li's tractable approximation (see Li and Oldenburg 1998) for a 3D field of prisms of constant susceptibility and density, and GeoBO applies this prism shape model to compute the corresponding sensor sensitivity for gravity and anomalous magnetic field related to each prism cell.

The current implementation includes magnetic and gravity forward models, which are defined in the module `sensormodel.py` by the functions `A_sens()`,`grav_func()`, and `magn_func()`. The easiest way to add custom models is to create a new forward model function similar to the included functions `grav_func()` or `magn_func` and to compute the forward model matrix with `A_sens()`, if possible. The custom function need to describe the sensitivity or relationship for a particular point relative to the sensor origin (see, e.g., `grav_func()`).

In general any linear forward model can be added by changing accordingly the forward model matrix as computed by `A_sens()` as long as this function returns the matrix **G** that satisfies the linear relation **y** = **G** ***P***.


### Gaussian Process Kernel Functions

Gaussian Processes (GPs) are a flexible, probabilistic approach using kernel machines and can propagate consistently uncertainties from input to output space under the Bayesian formalism. Another advantage of GPs is that their marginal likelihood function is well defined by the values of their hyper-parameters, and can thus be optimized.
The choice for an appropriate covariance (kernel) function is important and there are many stationary (invariant to translation in input space) and non-stationary covariance functions available (for an overview, see, e.g., Rasmussen and Williams 2006). To handle the computational problem of inverting a large covariance matrix, GeoBO uses by default an intrinsically sparse covariance function (Melkumyan et, al. 2009). However, other standard kernel functions are available (see module `kernels.py`), which includes the squared exponential 
and Matern32 function and their their corresponding multi-kernel covariance functions (see Melkumyan et. al. 2011). 

The settings yaml file allows you to choose the kernel function by configuring the parameter `kernelfunc`,  which can be set either to 'sparse' (Default), 'exp' (squared exponential) or 'matern32'. New custom kernels can be a added in the module `kernels.py`, which requires to write their covariance function (see as example `gpkernel()`) and cross-covariance function (see as example `gpkernel_sparse()`), and then to add their function name to settings.yaml and to `create_cov()` in `kernels.py`.  

The hyperparameters of the GP kernel can be configured in the settings yaml file (see Gaussian Process Settings) and are given by the lengthscale (`gp_lengthscale`), the noise regularization terms (`gp_err`) per forward model, and the cross-covariance amplitude terms which (`w1,w2,w3`) that coorrepond to the correlation coefficients between the model properties (e.g., rock density ,magnetic susceptibility, and drillcore measurements). The mathematical details for construction of the Multi-Kernel Covariance Functions are described in Haan et al 2020.


### Bayesian Optimisation Options

To find the optimal new sampling point, GeoBO maximises the objective function (acquisition function) which is defined by the Upper Confidence Bound (UCB)

UCB(x) = *m*(x) + *k* *sigma*(x) - *b* *c*(x)

with the mean value for the prediction *m*(x), the variance *sigma*<sup>2</sup>(x), and a cost function *c*(x), which is defined by the cost of obtaining a measurement at the sample point x. The parameter *k* and *b* define the trade-off in exploration-to-exploitation and gain-to-cost, respectively. For example, maximizing the mean value can be beneficial if the goal is to sample new data at locations with high density or mineral content, and not only where the uncertainty is high. 
The parameters *k* and *b* can be accordingly specified by the user in the settings yaml file.  Moreover, the settings allow the user to choose between vertical and non-vertical drillcore; in the latter case GeoBO is optimising also dip and azimuthal angle of the  drillcore in addition to drillcore position.


## Literature

Sebastian Haan, Fabio Ramos, Dietmar Muller, "Multi-Objective Bayesian Optimisation and Joint Inversion for Active Sensor Fusion", Geophysics, 86(1), pp.1-78. [arXiv Preprint](https://arxiv.org/abs/2010.05386)

Carl Edward Rasmussen and Christopher KI Williams, Gaussian process for machine learning, MIT press, 2006.

Li Yaoguo and Douglas W Oldenburg, “3d-inversion of gravity data,” Geophysics, vol. 63, no. 1, pp. 109–119, 1998.

Arman Melkumyan and Fabio Ramos, “A sparse covariance function for exact gaussian process inference in large datasets.,” in IJCAI, 2009, vol. 9, pp. 1936–1942

Armon Melkuyman and Fabio Ramos, “Multi-kernel gaussian processes,” in IJCAI, 2011, vol. 22, p. 1408

Reid, A., O. Simon Timothy, E. V. Bonilla, L. McCalman, T. Rawling, and F. Ramos, 2013, Bayesian joint inversions for the exploration of earth resources.: IJCAI, 2877

Eric Brochu, Vlad M Cora, and Nando De Freitas, “A tutorial on bayesian optimization of expensive cost functions, with application to active user modeling and hierarchical reinforcement learning,” arXiv preprint arXiv:1012.2599, 2010.


## Related Software

For the inversion part, GeoBO uses a direct inversion method via transformation of Gaussian Process priors, which enables joint inversion but is limited to linear forward models (e.g. gravity, magnetics, drillcores). For solving more complex non-linear forward models (e.g., seismic, or prior geological knowledge), the following bayesian inversion methods can potentially be applied to generate 3D geophysical surrogate models or to further refine GeoBo's 3D posterior model:

- hIPPYlib: an Extensible Software Framework for Large-scale Deterministic and Bayesian Inverse Problems. [Publication Link](https://www.theoj.org/joss-papers/joss.00940/10.21105.joss.00940.pdf); the software code is available at [hippylib.github.io](https://hippylib.github.io/)

- Obsidian: a flexible software platform for MCMC sampling of 3-D multi-modal geophysical models on distributed computing clusters. [Publication Link](https://gmd.copernicus.org/articles/12/2941/2019/); the code for version 0.1.2 of Obsidian is available at [https://doi.org/10.5281/zenodo.2580422](https://doi.org/10.5281/zenodo.2580422)

- GemPy: open-source stochastic geological modeling and inversion; geoscientific model development. See [gempy.org](https://www.gempy.org/)



## Attribution and Acknowledgments

Acknowledgments are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub.

If you make use of this code for your research project, please include the following acknowledgment:

“This research was supported by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney.”

### Project Contributors

Key project contributors to the GeoBO project are:

- Dr. Sebastian Haan (USYD, Sydney Informatics Hub): Expert in machine learning and physics, main contributor and software development of GeoBO.
- Prof. Fabian Ramos (USYD): Computational scientist and research expert in machine learning and bayesian computational techniques.
- Prof. Dietmar Muller (USYD, School of Geoscience): Research expert in geophysics and geoscience applications.
- Dr. Ben Mather (USYD, Sydney Informatics Hub/School of Geoscience ): Computational Geophysicist, GeoBO testing.


## License

Copyright 2020 Sebastian Haan

GeoBO is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License (AGPL version 3) as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program (see LICENSE.md). If not, see [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).


# Contributing 

We welcome all contributions to GeoBO. Contributions can be in the form of new code functions, improvements to the documentation, or by pointing out a bug or potential improvement.  


### Bugs reports and feature requests

For bugs and suggestions, please submit a bug report or feature request use the [Github issue tracker](https://github.com/sebhaan/geobo/issues). Search for existing and closed issues first. If your problem or idea is not yet addressed, please open a new issue. Github
allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.


### Submitting code

If you would like to contribute, please submit a pull request. (See the [Github Hello World](https://guides.github.com/activities/hello-world/) example, if you are new to Github). Please document your functions clearly.

By contributing to the repository you state you own the copyright to those contributions and agree to include your contributions as part of this project under the APGL license.

## Future Work

We welcome all contributions to GeoBO, in particular we would appreciate all pull requests pertaining to the following areas:

- Implementing new forward models in `sensormodel.py`
- Expanding the example studies for other use cases


### Testing

(Currently in development).
Continuous integration and testing is performed with `pytest`, which means running this in the source directory:

```bash
python setup.py test
```

The existing tests should be passing before you start coding (raise an issue if that is not the case) and any new functions should also have tests that we can use to verify the code. 

---
title: 'GeoBO: Python package for Multi-Objective Bayesian Optimisation and Joint Inversion in Geosciences'
tags:
  - Python
  - Geoscience
  - Bayesian Optimisation
  - Inversion
  - Gaussian Processes
authors:
  - name: Sebastian Haan
    orcid: 0000-0002-5994-5637
    affiliation: "1"
affiliations:
  - name: Sydney Informatics Hub, The University of Sydney, Australia
    index: 1
date: 14 November 2020
bibliography: paper.bib
---
<!-- pandoc -V geometry:margin=1in -V fontsize:11pt --filter pandoc-citeproc  -o paper.pdf paper.md -->
<!--add "author: Sebastian Haan" to meta for standard pandoc conversion --> 


# Summary

A critical decision process in geophysical data collection is how to efficiently combine a variety of sensor types and where to collect new observations, in particular if measurements are very costly or resources are limited. This may include drill-core placements or the layout of a geophysical survey. Bayesian Optimisation (BO) solves this decision making problem by finding global optimal solutions given multiple, often distinct measurements and model uncertainties. One of the major challenges for applying BO in geoscience still lies in constructing a 3D model from sparse, typically 2D sensor data that can be computationally efficiently evaluated and that transforms the combined data and uncertainties from multiple sensor measurements coherently into a posterior function approximating the true model; this is also known as an inversion problem.

``GeoBO`` is a Python package that solves both Multi-Objective Optimisation and Joint Inversion within one probabilistic framework: from prior selection and input, data fusion and inversion, to the final sensor optimisation and real world model output. 
<!-- While the current implementation generates 3D geophysical properties based on gravity and magnetic inversion and searching for optimal new drillcore measurements, in principle, the same model can be applied to a wide range of sensor fusion problems and allocation problems, such as: How to position or activate sensors for quasi-linear inverse problems (e.g. seismic, tomography) or if the model state is dynamic?  Where to sample if the cost function is incomplete? -->
Fig. 1 shows a graphical overview model of ``GeoBO``; the two main process steps are summarized in the following. 

First, ``GeoBO`` jointly solves multi-linear forward models of 2D survey data (e.g., magnetics and gravity) to 3D-geophysical properties using Gaussian Process kernels [@Rasmussen:2006; @Melkumyan:2009], without the requirement of prior geological knowledge. In addition, sparse drillhole measurements can be taken into account. One of the key features is to consider correlations between geophysical properties by solving simultaneously multiple forward models [@Melkumyan:2011; @Reid:2013]. Thus, this approach provides a better constrained combined solution of multiple, distinct sensor types, rather than taking their individual solutions. Another practical advantage of this probabilistic approach is that predictions are described by posterior distributions for each location (voxel or cube cell), which quantify the uncertainty in the predictions and their credible intervals. Most other work in the field of joint inversion only report point estimates of the quantities of interest [@Zeyen:1993; @Gallardo:2003; @Moorkamp:2011].

While this inversion method is limited to linear forward models, the reconstructed posterior distribution of the geophysical properties can be further refined using site-specific knowledge or more complex inversion methods for sampling non-linear forward models, such as hIPPYlib [@Hippylib], Obsidian [@Obsidian], or GemPy [@GemPy].

In a second step, the BO in `GeoBO` allows the user to define objectives in an acquisition function [@Mockus:1975; @Brochu:2010; @Haan:2020] which typically has to balance between a) exploration, i.e., querying points that maximise the information gain and minimize the uncertainty of a geophysical site model; b) exploitation, i.e., querying points that maximise the reward (e.g., concentrating search in the vicinity locations with high value such as minerals); and c) minimize the number of samples given a cost function for any new measurement (e.g., costly drill-cores). The reconstructed 3D model output of the joint inversion model is then used to query for the next most promising point based on the aquisition function, which guides the search for a user-defined optimum.


``GeoBO`` allows applications to specify priors as additional options, such as a) the typical correlation lengthscale of the geological structure via GP hyperparameters, b) the gain/cost parameters in the BO acquisition function, and c) the correlation coefficients between different geophysical properties. This package includes forward models for gravity and magnetics surveys [@Li:1998], drillcores, and test-sets of geological models and data. 

While it has been tested on solving the problem of allocating iteratively new expensive drill-core samples [@Haan:2020], the same method can be applied to a wide range of optimisation and sensor fusion problems. 


![Overview of the probabilistic framework for multi-objective optimisation and joint inversion; the process stages are: 1) Joint inversion via Gaussian Processes based on multiple sensor data fusion and forward model,  2) Posterior model generation of 3D multi-geophysical properties. 3) Maximisation of acquisition function to allocate optimal new sensor location and type. 4)  New acquired data is combined with existing data; process repeats until maximum number of iterations is achieved. For more details see @Haan:2020](graphmodel2.png)

## Installation, Dependencies and Usage
``GeoBO`` can be installed using setuptools or pip and requires the following python packages: numpy, matplotlib, scikit_image, scipy, rasterio, pandas, pyvista, skimage, PyYAML (for details see requirements.txt). The installation can be tested by running the script with included synthetic data and default settings.

To use ``GeoBO``, change the main settings such as filenames and parameters in settings.yaml (see description there), and run 
```sh
cd geobo/
python main.py settings.yaml
```

Core features of ``GeoBO``:

 - Joint probabilistic inversion tool by solving simultaneously multi-linear forward models (e.g. gravity, magnetics) using cross-variances between geophysical properties (cross-variance terms can be specified by use in settings.yaml)
 - Computation of complete posterior distribution for all geophysical properties (described by their mean and variance value at each location/voxel) 
 - Generation of ranked proposal list for new most promising drillcores (optional: non-vertical drillcores) based on global optimisation of acquisition function (see settings.yaml)
 - Templates for acquisition function to use in Bayesian Optimisation (see 'futility' functions in run_geobo.py)
 - Flexible parameter settings (settings.yaml) for exploration-exploitation trade-off and inclusion of local 3D cost function in acquisition function (see function 'create_costcube' in run_geobo.py)

Other features are:

 - Package includes geological survey/drillcore sample as well as synthetic data and functions for synthetic data generation
 - Generation of 2D/3D visualisation plots of reconstructed cubes and survey data (see plot settings in settings.yaml)
 - 3D Cube export in VTK format (for subsequent analysis, e.g., with ParaView)
 - Options to include any pre-existing drillcore data (see settings.yaml for feature names as well as format of sample drilldata)
 - Included linear forward models: density-to-gravity and magnetic susceptibility-to-magnetic field; custom linear forward models can be added (sensormodel.py)
 - Library of Gaussian Process (GP) kernels including sparse GP kernels (kernels.py)
 - Flexible settings for any cube geometry and resolution (cube length/width/height and voxel resolutions)
 - Optional  marginal GP likelihood for  optimization of GP hyperparameters and inversion process


## Results and Output

First, the joint inversion generates the reconstructed properties and their uncertainties as 3D cubes, such as density and magnetic susceptibility. Other 3D physical properties can be obtained by adding custom forward models (see section *Custom Linear Forward Models* in README). Secondly, Bayesian Optimisation uses these probabilistic inversion results to query new potential measurements based on a specific set of objectives as defined in an acquisition function, which guides the search for a user-defined optimum. Finally a list of all new measurement proposals ranked from maximum to minimum improvement is generated including a map of the most promising new measurement locations.


# Acknowledgements
Sebastian Haan would like to acknowledge Dietmar Muller, Fabio Ramos, and Ben Mather for their very valuable contributions to this project. This research was supported by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney.


# References
::: {role="main"}
::: {#section-intro .section}
Copyright 2020 Sebastian Haan

This script creates reconstructed Cubes with mean subtracted properties
of density, magnetic susceptibility, and drill core properties plus
their predicted variance cubes (see inversion.py for more details)

See settings.yaml for specifications.

ToDo: - Implement GP hyperparameter optimisation at adavanced user level

Expand source code

    """
    Copyright 2020 Sebastian Haan

    This script creates reconstructed Cubes with mean subtracted properties of density, magnetic susceptibility, 
    and  drill core properties plus their predicted variance cubes (see inversion.py for more details)

    See settings.yaml for specifications.

    ToDo: 
    - Implement GP hyperparameter optimisation at adavanced user level

    """
    import numpy as np
    import pandas as pd
    import rasterio
    from scipy import interpolate
    from scipy.ndimage.interpolation import zoom
    from scipy.optimize import minimize, least_squares, shgo
    from scipy.stats import norm, pareto
    from scipy import interpolate
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm

    # Import custom functions
    from config_loader import * # loads settings
    import cubeshow as cs 
    from utils import *
    import inversion
    import simcube



    def read_surveydata(plot = True):
            """
            Reading gravity data and magnetic data, including cropping and downsampling to Cube size

            PARAMETER
            :param plot: boolean, if True plots survey data and downsampled

            RETURN
            gravity data
            magnetic data
            sensor locations
            """
            if FNAME_gravsurvey is not None:
                    gravimg = rasterio.open(inpath + FNAME_gravsurvey)
                    grav = gravimg.read(1)
            else:
                    grav = None
            if FNAME_magsurvey is not None:
                    magimg = rasterio.open(inpath + FNAME_magsurvey)
                    mag = magimg.read(1)
            else:
                    mag = None
            # downsample to Cube voxelsize
            zoomfac = xNcube *1./grav.shape[1]
            grav2 = zoom(grav, zoomfac)
            assert grav2.shape == (yNcube, xNcube)
            zoomfac = xNcube *1./mag.shape[1]
            mag2 = zoom(mag, zoomfac)
            assert mag2.shape == (yNcube, xNcube)
            # Create survey sensor coordinates:
            x_s = np.linspace(0.5, xNcube - 0.5, xNcube)  * xvoxsize
            y_s = np.linspace(0.5, yNcube - 0.5, yNcube)  * yvoxsize
            z_s = zmax + zoff
            x_sensor, y_sensor, z_sensor = np.meshgrid(x_s,y_s,z_s)
            locations = np.asarray([x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten()]).T
            if not os.path.exists(outpath):
                    os.makedirs(outpath)
            if plot:
                    extent=[xmin,xmax,ymin,ymax]
                    plt.imshow(grav, aspect = 'equal', cmap = 'viridis', extent=extent)
                    plt.colorbar()
                    plt.savefig(outpath + 'gravfield.png')
                    plt.clf()
                    plt.imshow(mag, aspect = 'equal', cmap = 'viridis', extent=extent)
                    plt.colorbar()
                    plt.savefig(outpath + 'magfield.png')
                    plt.clf()
                    plt.imshow(grav2, aspect = 'equal', cmap = 'viridis', extent=extent)
                    plt.colorbar()
                    plt.savefig(outpath + 'gravfield_downsampled.png')
                    plt.clf()
                    plt.imshow(mag2, aspect = 'equal', cmap = 'viridis', extent=extent)
                    plt.colorbar()
                    plt.savefig(outpath + 'magfield_downsampled.png')
                    plt.close("all")
            return grav2.flatten(), mag2.flatten(), locations


    def read_drilldata(features):
            """
            Reading drill data

            PARAMETER
            :param features: List of drill features of interest (see settings.yaml)

            RETURN
            Drill data
            x,y,z coordinates of drilldata
            x,y,z min,max arrays of drilldata
            """
            if FNAME_drilldata is not None:
                    drill = pd.read_csv(inpath + FNAME_drilldata)
                    drill = drill[(drill.x >= xmin) & (drill.x <= xmax) & (drill.y >= ymin) & (drill.y <= ymax) & (drill.z <= zmax) & (drill.z >= zmin)].copy()
            else: 
                    drill = None
            # Select data only within extent
            # Convert to local coordinate system with origin X,Y = 0,0
            drill['x'] = drill['x'] - xmin
            drill['y'] = drill['y'] - ymin
            # Set up voxel coordinate grid
            xdrill = drill['x'].values
            ydrill = drill['y'].values
            zdrill = drill['z'].values 
            try:
                    drillfirst = drill.groupby('SiteID').first()
                    drilllast = drill.groupby('SiteID').last()
                    xdrillminmax = np.asarray([drillfirst.x.values, drilllast.x.values]).T
                    ydrillminmax = np.asarray([drillfirst.y.values, drilllast.y.values]).T
                    zdrillminmax = np.asarray([drillfirst.z.values, drilllast.z.values]).T
            except:
                    xdrillminmax = 0.#np.asarray([xdrill.min(), xdrill.min()])
                    ydrillminmax = 0.#np.asarray([ydrill.min(), ydrill.min()])
                    zdrillminmax = 0.#np.asarray([zdrill.min(), zdrill.min()])
            coord = np.vstack([xdrill, ydrill, zdrill]).T
            drilldata = []
            for feature in features:
                    data = drill[feature]
                    drilldata.append(align_drill(coord, data))
            return np.asarray(drilldata), coord, (xdrillminmax, ydrillminmax, zdrillminmax)


    def align_drill(coord, data):
            """
            Convert drill-core data in Model Cube shape with iteration over all voxels

            PARAMETER
            param coord: xyz Coordinates of drillcore, shape (N_drill, 3)
            param data: 1 dim data for drillcore, shape (N_drill)

            RETURN
            drilldata in voxel shape
            """
            dx = xvoxsize
            dy = yvoxsize
            dz = zvoxsize
            data = np.asarray(data)
            res = np.zeros_like(xxx)
            for ix in range(xxx.shape[0]):
                    for iy in range(xxx.shape[1]):
                            for iz in range(xxx.shape[2]):
                                    sel = np.where(((xxx[ix,iy,iz] - dx) <= coord[:, 0]) & (coord[:, 0] < (xxx[ix,iy,iz] + dx))
                                            & ((yyy[ix,iy,iz] - dy) <= coord[:, 1]) & (coord[:, 1] < (yyy[ix,iy,iz] + dy))
                                            & ((zzz[ix,iy,iz] - dz) <= coord[:, 2]) & (coord[:, 2] < (zzz[ix,iy,iz] + dz)))
                                    if np.size(data[sel]) > 0:
                                            #print(sel)
                                            m = np.nanmean(data[sel])
                                            if np.isfinite(m):
                                                    res[ix,iy,iz] = m
            return res



    """
    Below: defineition of acquisition function functions for BO

    The key of BO is the acquisition function, which typically has to balance between 
    a) exploration, i.e., querying points that maximise the information gain and minimize the uncertainty of a model, 
    b) exploitation, i.e. querying points that maximise the reward 
    (e.g. concentrating search in the vicinity locations with high value such as minerals), and 
    c) minimize the number of samples given an expensive cost function for any new measurement. 

    Exploration-exploitation and cost trade-off parameters can be set in settings.yaml
    """

    def futility_vertical(params, costs = None):
            """ 
            Utility/Aquisition function for bayesian optimisation assuming vertical drillcores

            PARAMETER
            params: parameters for drillhole (x_drill,y_drill) 
            param costs: cube of costs with same shape as reconstructed cube

            RETURN
            Output of utility function (scalar)
            """
            if costs is None:
                    costs = drill_rec * 0.
            params = np.asarray(params)
            xmaxvox = drill_rec.shape[0] - 1 
            ymaxvox = drill_rec.shape[1] - 1
            if np.isfinite(params).all():
                    xd = int(np.round(params[0]))
                    yd = int(np.round(params[1]))
                    if (xd > 0) & (xd < xmaxvox) & (yd > 0) & (yd < ymaxvox):
                            func = np.sum(drill_rec[xd, yd, :]) + kappa * np.sqrt(np.sum(drill_var[xd, yd, :])) - beta * np.sum(costs[xd,yd, :])
                    else:
                            func = -np.inf
            else:
                    func = -np.inf
            return -func # (negative if optimiser use minimize)


    def futility_drill(params, length_newdrill = zLcube, costs = None): 
            """
            Calculates utility function for proposed non-vertical drillcore, 
            which takes into account the azimuth and dip in addition to location and drill-core length

            PARAMETER
            param params: [x0, y0, azimuth, dip] x0 and y0 in meters; azimuth and dip in degree 
            param length_newdrill: length of drillcore in meters
            param costs: cube of costs with same shape as reconstructed cube

            RETURN
            Output of utility function (scalar)
            """
            if costs is None:
                    costs = drill_rec * 0.
            x0, y0, azimuth, dip = params
            nstep = int(2 * length_newdrill/np.min([xvoxsize, yvoxsize, zvoxsize]))
            rladder = np.linspace(0,length_newdrill, nstep)
            x0 = rladder * 0 + x0
            y0 = rladder * 0 + y0
            z0 = rladder * 0 + zmax
            azimuth = rladder * 0 + azimuth
            dip = rladder * 0 + dip
            try:
                    xdrillnew, ydrillnew, zdrillnew = spherical2cartes(x0, y0, z0, azimuth * np.pi/180., (180 - dip)* np.pi/180., rladder)
                    xnew = (xdrillnew/xvoxsize).astype(int)
                    ynew = (ydrillnew/yvoxsize).astype(int)
                    znew = (-zdrillnew/zvoxsize).astype(int)
                    funct = np.sum(drill_rec[xnew,ynew,znew])  + kappa * np.sqrt(np.sum(drill_var[xnew,ynew,znew])) - beta * np.sum(costs[xnew,ynew,znew])
            except:
                    funct = 0.
            return -funct

    """

    Chose in settings.yaml sepcifications and whether vertical or non-vertical drillcores

    The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
    This has the advantage of proposing and sorting multiple local maxima

    """

    def bayesopt_vert(drillcoord = None):
            """
            New Drillcore proposal based on mean, uncertainty, and costs as defined in acquisitian function futility_vertical

            The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
            This has the advantage of proposing and sorting multiple local maxima

            PARAMETER
            param drillcoord: exisiting x, y, x drill coordinates

            RETURN
            Saves list of proposal of new drillcore coordinates in file newdrill_vertical_proposals.csv

            """
            print('Calculating propoals list of new vertical drillholes...')
            # Define boundaries for optimisation:
            bp0 = np.asarray([yNcube/2, xNcube/2])
            blb = np.asarray([1, yNcube-1])
            bub = np.asarray([1, xNcube-1])
            #bopt_res = minimize(futility_vertical, bp0, bounds = (blb, bub),  method='SLSQP', options={'maxiter': 500}) #tol =1e-6, method='SLSQP'
            bopt_res = shgo(futility_vertical, bounds = ((1,yNcube - 1), (1,xNcube- 1)), n=20, iters=20, sampling_method='sobol' ) #tol =1e-6, method='SLSQP'
            if not bopt_res.success:
                    print('WARNING: ' + bopt_res.message) 
            print('New Drillcore Proposal:')
            print('_______________________')
            print('EASTING [meters]: ', np.round(bopt_res.x[1] * xvoxsize) + xmin + 0.5 * xvoxsize)
            print('NORTHING [meters]: ', np.round(bopt_res.x[0] * yvoxsize) + ymin + 0.5 * yvoxsize)

            # Write new drillcore proposals as csv file
            print('Saving all new drillhole proposals in file: newdrill_proposals.csv')
            newdf_values = bopt_res.xl 
            newdf_head = ['NORTHING','EASTING']
            df_newdrill = pd.DataFrame(np.round(newdf_values,2), columns = newdf_head)
            df_newdrill['EASTING'] = np.round(df_newdrill['EASTING']) * xvoxsize + xmin + 0.5 * xvoxsize
            df_newdrill['NORTHING'] = np.round(df_newdrill['NORTHING']) * yvoxsize + ymin + 0.5 * yvoxsize
            df_newdrill['BO_GAIN'] = -np.round(bopt_res.funl,4)
            df_newdrill.to_csv(outpath+ 'newdrill_vertical_proposals.csv', index=False)

            # Create image of porposed drillcore positions
            densimg = density_rec.mean(axis = 2)
            magsusimg = magsus_rec.mean(axis = 2)
            drillimg = drill_rec.mean(axis = 2)
            drillvarimg = drill_var.mean(axis = 2)
            #extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
            extent=[xmin, xmax , ymin , ymax]
            plt.clf()
            plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
            plt.xlabel('EASTING')
            plt.ylabel('NORTHING')
            if drillcoord is not None:
                    xdrill, ydrill, zdrill = drillcoord.T
                    plt.scatter(xdrill+xmin, ydrill+ymin, color='k')
            plt.scatter(df_newdrill['EASTING'].values, df_newdrill['NORTHING'].values, color='white')
            plt.scatter(df_newdrill.loc[0,'EASTING'], df_newdrill.loc[0,'NORTHING'], color='red')
            plt.title('Proposed Vertical Drillcores')
            plt.tight_layout()
            plt.savefig(outpath + 'newdrill_vertical_proposals.png')
            plt.close("all")


    def bayesopt_nonvert(drillcoord = None):
            """
            New Drillcore proposal based on mean, uncertainty, and costs as defined in acquisitian function futility_vertical

            The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
            This has the advantage of proposing and sorting multiple local maxima

            PARAMETER
            param drillcoord: exisiting x, y, x drill coordinates

            RETURN
            Saves list of proposal of new drillcore coordinates in file newdrill_vertical_proposals.csv
            """

            bnds = ((yvoxsize,yLcube - yvoxsize), (xvoxsize,xLcube- xvoxsize), (0,360),(30,90))
            bopt_res = shgo(futility_drill, bnds, n=10, iters=500, sampling_method='sobol')
            print('New Drillcore Proposal:')
            print('_______________________')
            print('EASTING [meters]: ', np.round(bopt_res.x[1] + ymin))
            print('NORTHING [meters]: ', np.round(bopt_res.x[0] + xmin))
            print('Azimuth Angle [degree]: ', np.round(bopt_res.x[2],1))
            print('Dip Angle [degree]: ', np.round(bopt_res.x[3],1))

            # Write new drillcore proposals as csv file
            print('Saving all new drillhole proposals in file: newdrill_proposals.csv')
            newdf_values = bopt_res.xl
            newdf_head = ['NORTHING', 'EASTING', 'AZIMUTH', 'DIP']
            df_newdrill = pd.DataFrame(np.round(newdf_values,2), columns = newdf_head)
            df_newdrill['EASTING'] = np.round(df_newdrill['EASTING'] + xmin,1)
            df_newdrill['NORTHING'] = np.round(df_newdrill['NORTHING'] + ymin,1)
            df_newdrill['BO_GAIN'] = -np.round(bopt_res.funl,4)
            df_newdrill.to_csv(outpath+ 'newdrill_proposals.csv', index=False)

            # Create image of prpposed drillcore positions
            densimg = density_rec.mean(axis = 2)
            magsusimg = magsus_rec.mean(axis = 2)
            drillimg = drill_rec.mean(axis = 2)
            drillvarimg = drill_var.mean(axis = 2)
            #extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
            extent=[xmin, xmax, ymin, ymax]
            plt.clf()
            plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
            plt.xlabel('EASTING')
            plt.ylabel('NORTHING')
            if drillcoord is not None:
                    xdrill, ydrill, zdrill = drillcoord.T
                    plt.scatter(xdrill+xmin, ydrill+ymin, color='k')
            plt.scatter(df_newdrill['EASTING'].values, df_newdrill['NORTHING'].values, color='white')
            plt.scatter(df_newdrill.loc[0,'EASTING'], df_newdrill.loc[0,'NORTHING'], color='red')
            plt.title('Proposed Drillcores')
            plt.tight_layout()
            plt.savefig(outpath + 'newdrill_proposals.png')
            plt.close('all')


    def create_costcube(cubeshape = (xNcube, yNcube, zNcube)):
            """
            User function to create costcube for drilling. Need to be same cube shape as reconstructed cube.
            By default set to zero costs.

            INPUT
            param cubeshape: shape of cube, by default: (xNcube, yNcube, zNcube)

            RETURN
            costcube
            """
            # Here costs are set to zero, change below
            costcube = np.zeros(cubeshape)

            return costcube



    ### Main computation ###

    if gen_simulation:
            # Create first simulated datacube, see settings.yaml
            simcube.create_simdata(modelname)  

    # Initiate inversion  class
    inv = inversion.Inversion()

    # Create new cube geometry
    voxelpos = inv.create_cubegeometry()
    xxx, yyy, zzz = voxelpos
    xxx = inv.xxx = xxx.reshape(xNcube, yNcube, zNcube)
    yyy = inv.yyy = yyy.reshape(xNcube, yNcube, zNcube)
    zzz = inv.zzz = zzz.reshape(xNcube, yNcube, zNcube)

    # Read in survey data
    gravfield, magfield, sensor_locations = read_surveydata()

    # # Read in existing drillcore data
    drilldata, drillcoord, drillminmax = read_drilldata(drill_features)
    drilldata0 = drilldata[ifeature]
    drillfield = drilldata0[drilldata0 != 0]
    xdrillminmax, ydrillminmax, zdrillminmax = drillminmax

    # Joint Inversion and reconmstyrcuting of cube with geophysical properties:
    density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var = inv.cubing(gravfield, magfield, drillfield, sensor_locations, drilldata0)

    ### Create VTK cubes or reconstructed cubes:  
    origin = (voxelpos[0].min(), voxelpos[1].min(), voxelpos[2].min())
    voxelsize = (xvoxsize, yvoxsize,zvoxsize)
    cs.create_vtkcube(density_rec, origin, voxelsize, fname = outpath + 'cube_density.vtk')
    cs.create_vtkcube(magsus_rec, origin, voxelsize, fname = outpath + 'cube_magsus.vtk')
    cs.create_vtkcube(drill_rec, origin, voxelsize, fname = outpath + 'cube_drill.vtk')
    cs.create_vtkcube(density_var, origin, voxelsize, fname = outpath + 'cube_density_variance.vtk')
    cs.create_vtkcube(magsus_var, origin, voxelsize, fname = outpath + 'cube_magsus_variance.vtk')
    cs.create_vtkcube(drill_var, origin, voxelsize, fname = outpath + 'cube_drill_variance.vtk')


    # Create plots of 2D maps of vertically integrated cube properties:
    if plot_vertical:       
            densimg = density_rec.mean(axis = 2)
            magsusimg = magsus_rec.mean(axis = 2)
            drillimg = drill_rec.mean(axis = 2)
            extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
            plt.clf()
            plt.imshow(densimg, aspect = 'equal', cmap = 'viridis', extent=extent)
            plt.colorbar()
            plt.savefig(outpath + 'dens_rec2D_loc2.png')
            plt.clf()
            plt.imshow(magsusimg, aspect = 'equal', cmap = 'viridis', extent=extent)
            plt.colorbar()
            plt.savefig(outpath + 'magsus_rec2D_loc2.png')
            plt.clf()
            plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent)
            plt.colorbar()
            plt.savefig(outpath + 'drill_rec2D_loc2.png')
            plt.close('all')


    # Create 3D Cube plot if needed (same data as in VTK cube)
    if plot3d:
                    plt.clf()
                    #cs.skplot3(density_rec, Nsize = (xNcube, yNcube, zNcube), drill = (yyydrill/xvoxsize,xxxdrill/yvoxsize,-zzzdrill/zvoxsize), sensor = (x_sensor/xvoxsize, y_sensor/yvoxsize, x_sensor * 0.), show = False, path_out = outpath, filename = 'density-drill-mesh.png')
                    cs.skplot3(density_rec, Nsize = (yNcube, xNcube, zNcube), drill = (ydrillminmax/xvoxsize, xdrillminmax/yvoxsize,-zdrillminmax/zvoxsize), sensor = (sensor_locations[1]/xvoxsize,sensor_locations[0]/yvoxsize, sensor_locations[2] * 0. + zmax), show = False, path_out = outpath, filename = 'density-mesh3D.png')
                    plt.clf()
                    cs.skplot3(drill_rec, Nsize = (yNcube, xNcube, zNcube), drill = (ydrillminmax/xvoxsize, xdrillminmax/yvoxsize,-zdrillminmax/zvoxsize), sensor = (sensor_locations[1]/xvoxsize,sensor_locations[0]/yvoxsize, sensor_locations[2] * 0. + zmax), show = False, path_out = outpath, filename = 'drill-mesh3D.png')
                    plt.close('all')


    # Default: no costs, change in function create_costcube to specify costs 
    costcube = create_costcube

    # Propose new drill-core, chose in settings.yaml sepcifications and whether vertical or non-vertical drillcores
    if bayesopt_vertical:
            bayesopt_vert(drillcoord)

    if bayesopt_nonvertical:
            bayesopt_nonvert(drillcoord)
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Functions {#header-functions .section-title}
---------

` def align_drill(coord, data)`{.name .flex}

:   ::: {.section .desc}
    Convert drill-core data in Model Cube shape with iteration over all
    voxels

    PARAMETER param coord: xyz Coordinates of drillcore, shape
    (N\_drill, 3) param data: 1 dim data for drillcore, shape (N\_drill)

    RETURN drilldata in voxel shape
    :::

    Expand source code

        def align_drill(coord, data):
                """
                Convert drill-core data in Model Cube shape with iteration over all voxels

                PARAMETER
                param coord: xyz Coordinates of drillcore, shape (N_drill, 3)
                param data: 1 dim data for drillcore, shape (N_drill)

                RETURN
                drilldata in voxel shape
                """
                dx = xvoxsize
                dy = yvoxsize
                dz = zvoxsize
                data = np.asarray(data)
                res = np.zeros_like(xxx)
                for ix in range(xxx.shape[0]):
                        for iy in range(xxx.shape[1]):
                                for iz in range(xxx.shape[2]):
                                        sel = np.where(((xxx[ix,iy,iz] - dx) <= coord[:, 0]) & (coord[:, 0] < (xxx[ix,iy,iz] + dx))
                                                & ((yyy[ix,iy,iz] - dy) <= coord[:, 1]) & (coord[:, 1] < (yyy[ix,iy,iz] + dy))
                                                & ((zzz[ix,iy,iz] - dz) <= coord[:, 2]) & (coord[:, 2] < (zzz[ix,iy,iz] + dz)))
                                        if np.size(data[sel]) > 0:
                                                #print(sel)
                                                m = np.nanmean(data[sel])
                                                if np.isfinite(m):
                                                        res[ix,iy,iz] = m
                return res

` def bayesopt_nonvert(drillcoord=None)`{.name .flex}

:   ::: {.section .desc}
    New Drillcore proposal based on mean, uncertainty, and costs as
    defined in acquisitian function futility\_vertical

    The optimisation of the acquisitioan function is performed with a
    global optimiser: SHGO with sobol sampling This has the advantage of
    proposing and sorting multiple local maxima

    PARAMETER param drillcoord: exisiting x, y, x drill coordinates

    RETURN Saves list of proposal of new drillcore coordinates in file
    newdrill\_vertical\_proposals.csv
    :::

    Expand source code

        def bayesopt_nonvert(drillcoord = None):
                """
                New Drillcore proposal based on mean, uncertainty, and costs as defined in acquisitian function futility_vertical

                The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
                This has the advantage of proposing and sorting multiple local maxima

                PARAMETER
                param drillcoord: exisiting x, y, x drill coordinates

                RETURN
                Saves list of proposal of new drillcore coordinates in file newdrill_vertical_proposals.csv
                """

                bnds = ((yvoxsize,yLcube - yvoxsize), (xvoxsize,xLcube- xvoxsize), (0,360),(30,90))
                bopt_res = shgo(futility_drill, bnds, n=10, iters=500, sampling_method='sobol')
                print('New Drillcore Proposal:')
                print('_______________________')
                print('EASTING [meters]: ', np.round(bopt_res.x[1] + ymin))
                print('NORTHING [meters]: ', np.round(bopt_res.x[0] + xmin))
                print('Azimuth Angle [degree]: ', np.round(bopt_res.x[2],1))
                print('Dip Angle [degree]: ', np.round(bopt_res.x[3],1))

                # Write new drillcore proposals as csv file
                print('Saving all new drillhole proposals in file: newdrill_proposals.csv')
                newdf_values = bopt_res.xl
                newdf_head = ['NORTHING', 'EASTING', 'AZIMUTH', 'DIP']
                df_newdrill = pd.DataFrame(np.round(newdf_values,2), columns = newdf_head)
                df_newdrill['EASTING'] = np.round(df_newdrill['EASTING'] + xmin,1)
                df_newdrill['NORTHING'] = np.round(df_newdrill['NORTHING'] + ymin,1)
                df_newdrill['BO_GAIN'] = -np.round(bopt_res.funl,4)
                df_newdrill.to_csv(outpath+ 'newdrill_proposals.csv', index=False)

                # Create image of prpposed drillcore positions
                densimg = density_rec.mean(axis = 2)
                magsusimg = magsus_rec.mean(axis = 2)
                drillimg = drill_rec.mean(axis = 2)
                drillvarimg = drill_var.mean(axis = 2)
                #extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
                extent=[xmin, xmax, ymin, ymax]
                plt.clf()
                plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
                plt.xlabel('EASTING')
                plt.ylabel('NORTHING')
                if drillcoord is not None:
                        xdrill, ydrill, zdrill = drillcoord.T
                        plt.scatter(xdrill+xmin, ydrill+ymin, color='k')
                plt.scatter(df_newdrill['EASTING'].values, df_newdrill['NORTHING'].values, color='white')
                plt.scatter(df_newdrill.loc[0,'EASTING'], df_newdrill.loc[0,'NORTHING'], color='red')
                plt.title('Proposed Drillcores')
                plt.tight_layout()
                plt.savefig(outpath + 'newdrill_proposals.png')
                plt.close('all')

` def bayesopt_vert(drillcoord=None)`{.name .flex}

:   ::: {.section .desc}
    New Drillcore proposal based on mean, uncertainty, and costs as
    defined in acquisitian function futility\_vertical

    The optimisation of the acquisitioan function is performed with a
    global optimiser: SHGO with sobol sampling This has the advantage of
    proposing and sorting multiple local maxima

    PARAMETER param drillcoord: exisiting x, y, x drill coordinates

    RETURN Saves list of proposal of new drillcore coordinates in file
    newdrill\_vertical\_proposals.csv
    :::

    Expand source code

        def bayesopt_vert(drillcoord = None):
                """
                New Drillcore proposal based on mean, uncertainty, and costs as defined in acquisitian function futility_vertical

                The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
                This has the advantage of proposing and sorting multiple local maxima

                PARAMETER
                param drillcoord: exisiting x, y, x drill coordinates

                RETURN
                Saves list of proposal of new drillcore coordinates in file newdrill_vertical_proposals.csv

                """
                print('Calculating propoals list of new vertical drillholes...')
                # Define boundaries for optimisation:
                bp0 = np.asarray([yNcube/2, xNcube/2])
                blb = np.asarray([1, yNcube-1])
                bub = np.asarray([1, xNcube-1])
                #bopt_res = minimize(futility_vertical, bp0, bounds = (blb, bub),  method='SLSQP', options={'maxiter': 500}) #tol =1e-6, method='SLSQP'
                bopt_res = shgo(futility_vertical, bounds = ((1,yNcube - 1), (1,xNcube- 1)), n=20, iters=20, sampling_method='sobol' ) #tol =1e-6, method='SLSQP'
                if not bopt_res.success:
                        print('WARNING: ' + bopt_res.message) 
                print('New Drillcore Proposal:')
                print('_______________________')
                print('EASTING [meters]: ', np.round(bopt_res.x[1] * xvoxsize) + xmin + 0.5 * xvoxsize)
                print('NORTHING [meters]: ', np.round(bopt_res.x[0] * yvoxsize) + ymin + 0.5 * yvoxsize)

                # Write new drillcore proposals as csv file
                print('Saving all new drillhole proposals in file: newdrill_proposals.csv')
                newdf_values = bopt_res.xl 
                newdf_head = ['NORTHING','EASTING']
                df_newdrill = pd.DataFrame(np.round(newdf_values,2), columns = newdf_head)
                df_newdrill['EASTING'] = np.round(df_newdrill['EASTING']) * xvoxsize + xmin + 0.5 * xvoxsize
                df_newdrill['NORTHING'] = np.round(df_newdrill['NORTHING']) * yvoxsize + ymin + 0.5 * yvoxsize
                df_newdrill['BO_GAIN'] = -np.round(bopt_res.funl,4)
                df_newdrill.to_csv(outpath+ 'newdrill_vertical_proposals.csv', index=False)

                # Create image of porposed drillcore positions
                densimg = density_rec.mean(axis = 2)
                magsusimg = magsus_rec.mean(axis = 2)
                drillimg = drill_rec.mean(axis = 2)
                drillvarimg = drill_var.mean(axis = 2)
                #extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
                extent=[xmin, xmax , ymin , ymax]
                plt.clf()
                plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
                plt.xlabel('EASTING')
                plt.ylabel('NORTHING')
                if drillcoord is not None:
                        xdrill, ydrill, zdrill = drillcoord.T
                        plt.scatter(xdrill+xmin, ydrill+ymin, color='k')
                plt.scatter(df_newdrill['EASTING'].values, df_newdrill['NORTHING'].values, color='white')
                plt.scatter(df_newdrill.loc[0,'EASTING'], df_newdrill.loc[0,'NORTHING'], color='red')
                plt.title('Proposed Vertical Drillcores')
                plt.tight_layout()
                plt.savefig(outpath + 'newdrill_vertical_proposals.png')
                plt.close("all")

` def costcube(cubeshape=(25, 16, 16))`{.name .flex}

:   ::: {.section .desc}
    User function to create costcube for drilling. Need to be same cube
    shape as reconstructed cube. By default set to zero costs.

    INPUT param cubeshape: shape of cube, by default: (xNcube, yNcube,
    zNcube)

    RETURN costcube
    :::

    Expand source code

        def create_costcube(cubeshape = (xNcube, yNcube, zNcube)):
                """
                User function to create costcube for drilling. Need to be same cube shape as reconstructed cube.
                By default set to zero costs.

                INPUT
                param cubeshape: shape of cube, by default: (xNcube, yNcube, zNcube)

                RETURN
                costcube
                """
                # Here costs are set to zero, change below
                costcube = np.zeros(cubeshape)

                return costcube

` def create_costcube(cubeshape=(25, 16, 16))`{.name .flex}

:   ::: {.section .desc}
    User function to create costcube for drilling. Need to be same cube
    shape as reconstructed cube. By default set to zero costs.

    INPUT param cubeshape: shape of cube, by default: (xNcube, yNcube,
    zNcube)

    RETURN costcube
    :::

    Expand source code

        def create_costcube(cubeshape = (xNcube, yNcube, zNcube)):
                """
                User function to create costcube for drilling. Need to be same cube shape as reconstructed cube.
                By default set to zero costs.

                INPUT
                param cubeshape: shape of cube, by default: (xNcube, yNcube, zNcube)

                RETURN
                costcube
                """
                # Here costs are set to zero, change below
                costcube = np.zeros(cubeshape)

                return costcube

` def futility_drill(params, length_newdrill=800.0, costs=None)`{.name .flex}

:   ::: {.section .desc}
    Calculates utility function for proposed non-vertical drillcore,
    which takes into account the azimuth and dip in addition to location
    and drill-core length

    PARAMETER param params: \[x0, y0, azimuth, dip\] x0 and y0 in
    meters; azimuth and dip in degree param length\_newdrill: length of
    drillcore in meters param costs: cube of costs with same shape as
    reconstructed cube

    RETURN Output of utility function (scalar)
    :::

    Expand source code

        def futility_drill(params, length_newdrill = zLcube, costs = None): 
                """
                Calculates utility function for proposed non-vertical drillcore, 
                which takes into account the azimuth and dip in addition to location and drill-core length

                PARAMETER
                param params: [x0, y0, azimuth, dip] x0 and y0 in meters; azimuth and dip in degree 
                param length_newdrill: length of drillcore in meters
                param costs: cube of costs with same shape as reconstructed cube

                RETURN
                Output of utility function (scalar)
                """
                if costs is None:
                        costs = drill_rec * 0.
                x0, y0, azimuth, dip = params
                nstep = int(2 * length_newdrill/np.min([xvoxsize, yvoxsize, zvoxsize]))
                rladder = np.linspace(0,length_newdrill, nstep)
                x0 = rladder * 0 + x0
                y0 = rladder * 0 + y0
                z0 = rladder * 0 + zmax
                azimuth = rladder * 0 + azimuth
                dip = rladder * 0 + dip
                try:
                        xdrillnew, ydrillnew, zdrillnew = spherical2cartes(x0, y0, z0, azimuth * np.pi/180., (180 - dip)* np.pi/180., rladder)
                        xnew = (xdrillnew/xvoxsize).astype(int)
                        ynew = (ydrillnew/yvoxsize).astype(int)
                        znew = (-zdrillnew/zvoxsize).astype(int)
                        funct = np.sum(drill_rec[xnew,ynew,znew])  + kappa * np.sqrt(np.sum(drill_var[xnew,ynew,znew])) - beta * np.sum(costs[xnew,ynew,znew])
                except:
                        funct = 0.
                return -funct

` def futility_vertical(params, costs=None)`{.name .flex}

:   ::: {.section .desc}
    Utility/Aquisition function for bayesian optimisation assuming
    vertical drillcores

    PARAMETER params: parameters for drillhole (x\_drill,y\_drill) param
    costs: cube of costs with same shape as reconstructed cube

    RETURN Output of utility function (scalar)
    :::

    Expand source code

        def futility_vertical(params, costs = None):
                """ 
                Utility/Aquisition function for bayesian optimisation assuming vertical drillcores

                PARAMETER
                params: parameters for drillhole (x_drill,y_drill) 
                param costs: cube of costs with same shape as reconstructed cube

                RETURN
                Output of utility function (scalar)
                """
                if costs is None:
                        costs = drill_rec * 0.
                params = np.asarray(params)
                xmaxvox = drill_rec.shape[0] - 1 
                ymaxvox = drill_rec.shape[1] - 1
                if np.isfinite(params).all():
                        xd = int(np.round(params[0]))
                        yd = int(np.round(params[1]))
                        if (xd > 0) & (xd < xmaxvox) & (yd > 0) & (yd < ymaxvox):
                                func = np.sum(drill_rec[xd, yd, :]) + kappa * np.sqrt(np.sum(drill_var[xd, yd, :])) - beta * np.sum(costs[xd,yd, :])
                        else:
                                func = -np.inf
                else:
                        func = -np.inf
                return -func # (negative if optimiser use minimize)

` def read_drilldata(features)`{.name .flex}

:   ::: {.section .desc}
    Reading drill data

    PARAMETER :param features: List of drill features of interest (see
    settings.yaml)

    RETURN Drill data x,y,z coordinates of drilldata x,y,z min,max
    arrays of drilldata
    :::

    Expand source code

        def read_drilldata(features):
                """
                Reading drill data

                PARAMETER
                :param features: List of drill features of interest (see settings.yaml)

                RETURN
                Drill data
                x,y,z coordinates of drilldata
                x,y,z min,max arrays of drilldata
                """
                if FNAME_drilldata is not None:
                        drill = pd.read_csv(inpath + FNAME_drilldata)
                        drill = drill[(drill.x >= xmin) & (drill.x <= xmax) & (drill.y >= ymin) & (drill.y <= ymax) & (drill.z <= zmax) & (drill.z >= zmin)].copy()
                else: 
                        drill = None
                # Select data only within extent
                # Convert to local coordinate system with origin X,Y = 0,0
                drill['x'] = drill['x'] - xmin
                drill['y'] = drill['y'] - ymin
                # Set up voxel coordinate grid
                xdrill = drill['x'].values
                ydrill = drill['y'].values
                zdrill = drill['z'].values 
                try:
                        drillfirst = drill.groupby('SiteID').first()
                        drilllast = drill.groupby('SiteID').last()
                        xdrillminmax = np.asarray([drillfirst.x.values, drilllast.x.values]).T
                        ydrillminmax = np.asarray([drillfirst.y.values, drilllast.y.values]).T
                        zdrillminmax = np.asarray([drillfirst.z.values, drilllast.z.values]).T
                except:
                        xdrillminmax = 0.#np.asarray([xdrill.min(), xdrill.min()])
                        ydrillminmax = 0.#np.asarray([ydrill.min(), ydrill.min()])
                        zdrillminmax = 0.#np.asarray([zdrill.min(), zdrill.min()])
                coord = np.vstack([xdrill, ydrill, zdrill]).T
                drilldata = []
                for feature in features:
                        data = drill[feature]
                        drilldata.append(align_drill(coord, data))
                return np.asarray(drilldata), coord, (xdrillminmax, ydrillminmax, zdrillminmax)

` def read_surveydata(plot=True)`{.name .flex}

:   ::: {.section .desc}
    Reading gravity data and magnetic data, including cropping and
    downsampling to Cube size

    PARAMETER :param plot: boolean, if True plots survey data and
    downsampled

    RETURN gravity data magnetic data sensor locations
    :::

    Expand source code

        def read_surveydata(plot = True):
                """
                Reading gravity data and magnetic data, including cropping and downsampling to Cube size

                PARAMETER
                :param plot: boolean, if True plots survey data and downsampled

                RETURN
                gravity data
                magnetic data
                sensor locations
                """
                if FNAME_gravsurvey is not None:
                        gravimg = rasterio.open(inpath + FNAME_gravsurvey)
                        grav = gravimg.read(1)
                else:
                        grav = None
                if FNAME_magsurvey is not None:
                        magimg = rasterio.open(inpath + FNAME_magsurvey)
                        mag = magimg.read(1)
                else:
                        mag = None
                # downsample to Cube voxelsize
                zoomfac = xNcube *1./grav.shape[1]
                grav2 = zoom(grav, zoomfac)
                assert grav2.shape == (yNcube, xNcube)
                zoomfac = xNcube *1./mag.shape[1]
                mag2 = zoom(mag, zoomfac)
                assert mag2.shape == (yNcube, xNcube)
                # Create survey sensor coordinates:
                x_s = np.linspace(0.5, xNcube - 0.5, xNcube)  * xvoxsize
                y_s = np.linspace(0.5, yNcube - 0.5, yNcube)  * yvoxsize
                z_s = zmax + zoff
                x_sensor, y_sensor, z_sensor = np.meshgrid(x_s,y_s,z_s)
                locations = np.asarray([x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten()]).T
                if not os.path.exists(outpath):
                        os.makedirs(outpath)
                if plot:
                        extent=[xmin,xmax,ymin,ymax]
                        plt.imshow(grav, aspect = 'equal', cmap = 'viridis', extent=extent)
                        plt.colorbar()
                        plt.savefig(outpath + 'gravfield.png')
                        plt.clf()
                        plt.imshow(mag, aspect = 'equal', cmap = 'viridis', extent=extent)
                        plt.colorbar()
                        plt.savefig(outpath + 'magfield.png')
                        plt.clf()
                        plt.imshow(grav2, aspect = 'equal', cmap = 'viridis', extent=extent)
                        plt.colorbar()
                        plt.savefig(outpath + 'gravfield_downsampled.png')
                        plt.clf()
                        plt.imshow(mag2, aspect = 'equal', cmap = 'viridis', extent=extent)
                        plt.colorbar()
                        plt.savefig(outpath + 'magfield_downsampled.png')
                        plt.close("all")
                return grav2.flatten(), mag2.flatten(), locations
:::

::: {.section}
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Functions](#header-functions) {#functions}

    -   `align_drill`
    -   `bayesopt_nonvert`
    -   `bayesopt_vert`
    -   `costcube`
    -   `create_costcube`
    -   `futility_drill`
    -   `futility_vertical`
    -   `read_drilldata`
    -   `read_surveydata`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
::: {role="main"}
::: {#section-intro .section}
Script for running inversion and reconstructing 3D cubes from 2D sensor
data using Gaussina processes

Author: Sebastian Haan

Expand source code

    """
    Script for running inversion and reconstructing 3D cubes from 2D sensor data using Gaussina processes

    Author: Sebastian Haan
    """

    import numpy as np
    import sys
    from scipy.linalg import pinv, solve, cholesky, solve_triangular
    from config_loader import *
    import kernels as kernel #for Gaussian Process prior
    import sensormodel

    class Inversion:
            """
            Class for Inversion and reconstructon of 3D cubes from 2D sensor data.

            To solve the inverse problem for linear systems, we deploy a Bayesian framework with Gaussian process priors 
            over the physical properties of interest. Gaussian processes (GP) are a flexible, probabilistic approach 
            using kernel machines with non-parametric priors (see e.g. Rasmussen 2006).  An important advantage of 
            the Bayesian method is that it generates a predictive distribution with a mean and variance, 
            which are prerequisites for Bayesian optimisation with respect to, e.g., information gain from a new measurement. 
            Another advantage is that the GP marginal likelihood function is well defined by the values of their hyper-parameters, 
            which allows it to optimise them exactly. This reasoning about functions under uncertainty and their well-tuned 
            interpolation character allows GPs to work extremely well for sparse data as well as for big data 
            (using sparse covariance matrix, see e.g. Melkumyan 2000).

            To take fully into account cross-covariances between multiple model parameters 
            (e.g., rock density and magnetic susceptibility), we construct sparse cross-covariance terms 
            between all kernel pairs following Melkumyan 2011. The full covariance matrix is then constructed 
            by setting sparse covariance functions on all diagonal block elements and sparse-sparse cross-covariance functions 
            on all non-diagonal block elements. Cross-covariance amplitude terms are given by the linear correlation coefficients 
            between the corresponding geophysical properties.

            Call function cubing() for running inversion.
            """
            def __init__(self):
                    # Set Lengthscale and amplitude for GP
                    self.gp_length = gp_lengthscale * np.asarray([xvoxsize, xvoxsize, xvoxsize]) # 2*voxelsize seems to have max logl
                    self.gp_sigma = np.asarray(gp_err)
                    self.coeffm = np.asarray(gp_coeff) # coefficients for GP kernel mix


            def create_cubegeometry(self):
                    """
                    Create Cube geometry and coordinates

                    See settings.yaml for cube specifications
                    """
                    # Create voxel Edge coordinates:
                    xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
                    yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
                    zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
                    #zedge = np.linspace(0.,Ncube/zfactor, Ncube/zfactor+1) * voxsize
                    xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
                    self.Edges = np.asarray([xEdges, yEdges, -zEdges])
                    # Create voxel coordinates
                    xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
                    ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
                    znew = zmax - np.arange(zvoxsize/2., zLcube + zvoxsize/2., zvoxsize)
                    xx ,yy = np.meshgrid(xnew,ynew)
                    self.xxx, self.yyy, self.zzz = np.meshgrid(xnew, ynew, znew)
                    self.voxelpos = np.vstack([self.xxx.flatten(), self.yyy.flatten(), self.zzz.flatten()])
                    return self.voxelpos


            def predict3(self, calclogl = False):
                    """
                    Calculate mean, variance, and likelihood likelihood for GP with 3x3 kernel block matrix 
                
                PARAMETERS

                :param calclogl: if True calculate marginal GP loglikelihood, if False logl return is set to 0

                RETURN: 
                mean
                covariance
                log-likelihood
                    """
                    # Calculate kernel 3x3 block matrix:
                    self.datastd = np.mean([np.nanstd(self.gravfield), np.nanstd(self.magfield), np.nanstd(self.drillfield)])
                    self.kcov = kernel.create_cov(self.D2, self.gp_length, self.coeffm, fkernel = kernelfunc)
                    #Ak = np.dot(self.Asens3, kcov) # shape (3*Nsensor, 3*Ncube^3)
                    yerr = np.hstack((self.gravfield * 0. + self.gp_sigma[0], self.magfield * 0. + self.gp_sigma[1], self.drillfield *0. + self.gp_sigma[2]))
                    #try:
                    AkA = np.dot(self.Asens3, np.dot(self.kcov, self.Asens3.T)) + np.diag(yerr**2)
                    #AkA = np.dot(Ak, self.Asens3.T) + np.diag(yerr**2) # shape(2*Nsensor, 2*Nsensor)
                    # Cholesky decomposition
                    try:
                            AkA_chol = cholesky(AkA, lower=True)
                    except:
                            print("Cholesky decompostion failed, AkA matrix i likely not positive semitive.")
                            print("Change GP parameter settings")
                            sys.exit(1)
                    usolve = solve_triangular(AkA_chol, self.Fs3, lower=True) #shape(2*Nsensor)
                    # Calculate Likelihood if option is set
                    if calclogl:
                            log_det_AkA = np.log(np.diag(AkA_chol)**2).sum()
                            n_log_2pi = xNcube * yNcube * zNcube * np.log(2 * np.pi)
                            logl = -0.5 * (np.dot(usolve, usolve) + log_det_AkA + n_log_2pi)
                    else:
                            logl = 0.
                    # predicted mean
                    self._V = solve_triangular(AkA_chol, np.dot(self.Asens3, self.kcov), lower=True) #shape(2*Nsensor, 2* Ncube^3)
                    mu = np.dot(self._V.T, usolve) #* _fstd + _fmean
                    # covariance
                    covar = self.kcov - np.dot(self._V.T, self._V)
                    #var = np.diagonal(covar) #* _fstd**2
                    #except:
                    #       print("Warning: GP Inversion failed!")
                    #       mu, covar, logl = np.zeros(3*xNcube*yNcube*zNcube), kcov * 1e9, np.inf 
                    return mu, covar, logl

            
            def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0):
                    """
                    Joint Inversion and Cubing of sensor data

                    PARAMETERS

                    :param gravfield: 2D gravitational survey data
                    :param magfield: 2D magnetic survey data
                    :param drillfield: drilldata

                    RETURN

                    density_rec: cube with mean density 
                    magsus_rec: cube with mean magnetic susceptibility
                    drill_rec: cube with mean drill-core property  
                    density_var: varicance cube with density 
                    magsus_var: variance cube with magnetic susceptibility  
                    drill_var: varianbce cube with drill-core property  


                    """
                    self.gravfield = gravfield
                    self.magfield = magfield
                    self.drillfield = drillfield
                    self.sensor_locations = sensor_locations
                    self.drilldata0 = drilldata0
                    # Normalize data
                    grav_mean, grav_std = self.gravfield.mean(), self.gravfield.std()
                    gravfield_norm = (self.gravfield - grav_mean) / grav_std
                    magn_mean, magn_std = self.magfield.mean(), self.magfield.std()
                    magfield_norm = (self.magfield - magn_mean) / magn_std
                    drill_mean, drill_std = self.drillfield.mean(), self.drillfield.std()
                    drillfield_norm = (self.drillfield - drill_mean) / drill_std
                    # Create kernel distance matrix
                    self.points3D = kernel.calcGridPoints3D((xNcube, yNcube, zNcube), (xvoxsize, yvoxsize, zvoxsize))
                    self.D2 = kernel.calcDistanceMatrix(self.points3D)
                    # Combine Data and Forward models
                    voxelpos_drill = np.vstack([self.xxx[self.drilldata0 != 0], self.yyy[self.drilldata0 != 0], self.zzz[self.drilldata0 != 0]])
                    # Combine data:
                    self.Fs3 = np.hstack((gravfield_norm, magfield_norm, drillfield_norm))
                    # Calculate Sensitivity matrix
                    Asens_grav, _ = A_sens(magneticField * 0., self.sensor_locations, self.Edges, 'grav')
                    Asens_mag, _ = A_sens(magneticField, self.sensor_locations, self.Edges, 'magn')
                    Asens_drill = A_drill(voxelpos_drill, self.voxelpos)
                    # Calculate total transformation matrix
                    Asens_g= np.vstack([Asens_grav, np.zeros_like(Asens_mag), np.zeros_like(Asens_drill)])
                    Asens_m= np.vstack([np.zeros_like(Asens_grav), Asens_mag, np.zeros_like(Asens_drill)])
                    Asens_d= np.vstack([np.zeros_like(Asens_grav), np.zeros_like(Asens_mag), Asens_drill])
                    self.Asens3 = np.hstack([Asens_g, Asens_m, Asens_d])
                    # Run actual GP inversion
                    self.mu_rec, self.cov_rec, self.logl = self.predict3(calclogl = True)
                    # Reshape results and multiply with stadnard deviatio 
                    results_rec = self.mu_rec.reshape(3, yNcube, xNcube, zNcube)
                    results_var = np.diag(self.cov_rec).reshape(3, yNcube, xNcube, zNcube)
                    #Cut off outer voxel rows
                    #results_rec = results_rec[:, 1:yNcube-1, 1:xNcube-1, :]
                    #results_var = results_var[:, 1:yNcube-1, 1:xNcube-1, :]
                    density_rec = results_rec[0] * grav_std  # Model represents deviation from the mean
                    density_var = results_var[0] * grav_std**2
                    magsus_rec = results_rec[1] * magn_std # Model represents deviation from the mean
                    magsus_var = results_var[1] * magn_std**2
                    drill_rec = results_rec[2] * drill_std  # Model represents deviation from the mean
                    drill_var = results_var[2] * drill_std**2
                    return density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Classes {#header-classes .section-title}
-------

` class Inversion`{.flex .name .class}

:   ::: {.section .desc}
    Class for Inversion and reconstructon of 3D cubes from 2D sensor
    data.

    To solve the inverse problem for linear systems, we deploy a
    Bayesian framework with Gaussian process priors over the physical
    properties of interest. Gaussian processes (GP) are a flexible,
    probabilistic approach using kernel machines with non-parametric
    priors (see e.g. Rasmussen 2006). An important advantage of the
    Bayesian method is that it generates a predictive distribution with
    a mean and variance, which are prerequisites for Bayesian
    optimisation with respect to, e.g., information gain from a new
    measurement. Another advantage is that the GP marginal likelihood
    function is well defined by the values of their hyper-parameters,
    which allows it to optimise them exactly. This reasoning about
    functions under uncertainty and their well-tuned interpolation
    character allows GPs to work extremely well for sparse data as well
    as for big data (using sparse covariance matrix, see e.g. Melkumyan
    2000).

    To take fully into account cross-covariances between multiple model
    parameters (e.g., rock density and magnetic susceptibility), we
    construct sparse cross-covariance terms between all kernel pairs
    following Melkumyan 2011. The full covariance matrix is then
    constructed by setting sparse covariance functions on all diagonal
    block elements and sparse-sparse cross-covariance functions on all
    non-diagonal block elements. Cross-covariance amplitude terms are
    given by the linear correlation coefficients between the
    corresponding geophysical properties.

    Call function cubing() for running inversion.
    :::

    Expand source code

        class Inversion:
                """
                Class for Inversion and reconstructon of 3D cubes from 2D sensor data.

                To solve the inverse problem for linear systems, we deploy a Bayesian framework with Gaussian process priors 
                over the physical properties of interest. Gaussian processes (GP) are a flexible, probabilistic approach 
                using kernel machines with non-parametric priors (see e.g. Rasmussen 2006).  An important advantage of 
                the Bayesian method is that it generates a predictive distribution with a mean and variance, 
                which are prerequisites for Bayesian optimisation with respect to, e.g., information gain from a new measurement. 
                Another advantage is that the GP marginal likelihood function is well defined by the values of their hyper-parameters, 
                which allows it to optimise them exactly. This reasoning about functions under uncertainty and their well-tuned 
                interpolation character allows GPs to work extremely well for sparse data as well as for big data 
                (using sparse covariance matrix, see e.g. Melkumyan 2000).

                To take fully into account cross-covariances between multiple model parameters 
                (e.g., rock density and magnetic susceptibility), we construct sparse cross-covariance terms 
                between all kernel pairs following Melkumyan 2011. The full covariance matrix is then constructed 
                by setting sparse covariance functions on all diagonal block elements and sparse-sparse cross-covariance functions 
                on all non-diagonal block elements. Cross-covariance amplitude terms are given by the linear correlation coefficients 
                between the corresponding geophysical properties.

                Call function cubing() for running inversion.
                """
                def __init__(self):
                        # Set Lengthscale and amplitude for GP
                        self.gp_length = gp_lengthscale * np.asarray([xvoxsize, xvoxsize, xvoxsize]) # 2*voxelsize seems to have max logl
                        self.gp_sigma = np.asarray(gp_err)
                        self.coeffm = np.asarray(gp_coeff) # coefficients for GP kernel mix


                def create_cubegeometry(self):
                        """
                        Create Cube geometry and coordinates

                        See settings.yaml for cube specifications
                        """
                        # Create voxel Edge coordinates:
                        xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
                        yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
                        zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
                        #zedge = np.linspace(0.,Ncube/zfactor, Ncube/zfactor+1) * voxsize
                        xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
                        self.Edges = np.asarray([xEdges, yEdges, -zEdges])
                        # Create voxel coordinates
                        xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
                        ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
                        znew = zmax - np.arange(zvoxsize/2., zLcube + zvoxsize/2., zvoxsize)
                        xx ,yy = np.meshgrid(xnew,ynew)
                        self.xxx, self.yyy, self.zzz = np.meshgrid(xnew, ynew, znew)
                        self.voxelpos = np.vstack([self.xxx.flatten(), self.yyy.flatten(), self.zzz.flatten()])
                        return self.voxelpos


                def predict3(self, calclogl = False):
                        """
                        Calculate mean, variance, and likelihood likelihood for GP with 3x3 kernel block matrix 
                    
                    PARAMETERS

                    :param calclogl: if True calculate marginal GP loglikelihood, if False logl return is set to 0

                    RETURN: 
                    mean
                    covariance
                    log-likelihood
                        """
                        # Calculate kernel 3x3 block matrix:
                        self.datastd = np.mean([np.nanstd(self.gravfield), np.nanstd(self.magfield), np.nanstd(self.drillfield)])
                        self.kcov = kernel.create_cov(self.D2, self.gp_length, self.coeffm, fkernel = kernelfunc)
                        #Ak = np.dot(self.Asens3, kcov) # shape (3*Nsensor, 3*Ncube^3)
                        yerr = np.hstack((self.gravfield * 0. + self.gp_sigma[0], self.magfield * 0. + self.gp_sigma[1], self.drillfield *0. + self.gp_sigma[2]))
                        #try:
                        AkA = np.dot(self.Asens3, np.dot(self.kcov, self.Asens3.T)) + np.diag(yerr**2)
                        #AkA = np.dot(Ak, self.Asens3.T) + np.diag(yerr**2) # shape(2*Nsensor, 2*Nsensor)
                        # Cholesky decomposition
                        try:
                                AkA_chol = cholesky(AkA, lower=True)
                        except:
                                print("Cholesky decompostion failed, AkA matrix i likely not positive semitive.")
                                print("Change GP parameter settings")
                                sys.exit(1)
                        usolve = solve_triangular(AkA_chol, self.Fs3, lower=True) #shape(2*Nsensor)
                        # Calculate Likelihood if option is set
                        if calclogl:
                                log_det_AkA = np.log(np.diag(AkA_chol)**2).sum()
                                n_log_2pi = xNcube * yNcube * zNcube * np.log(2 * np.pi)
                                logl = -0.5 * (np.dot(usolve, usolve) + log_det_AkA + n_log_2pi)
                        else:
                                logl = 0.
                        # predicted mean
                        self._V = solve_triangular(AkA_chol, np.dot(self.Asens3, self.kcov), lower=True) #shape(2*Nsensor, 2* Ncube^3)
                        mu = np.dot(self._V.T, usolve) #* _fstd + _fmean
                        # covariance
                        covar = self.kcov - np.dot(self._V.T, self._V)
                        #var = np.diagonal(covar) #* _fstd**2
                        #except:
                        #       print("Warning: GP Inversion failed!")
                        #       mu, covar, logl = np.zeros(3*xNcube*yNcube*zNcube), kcov * 1e9, np.inf 
                        return mu, covar, logl

                
                def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0):
                        """
                        Joint Inversion and Cubing of sensor data

                        PARAMETERS

                        :param gravfield: 2D gravitational survey data
                        :param magfield: 2D magnetic survey data
                        :param drillfield: drilldata

                        RETURN

                        density_rec: cube with mean density 
                        magsus_rec: cube with mean magnetic susceptibility
                        drill_rec: cube with mean drill-core property  
                        density_var: varicance cube with density 
                        magsus_var: variance cube with magnetic susceptibility  
                        drill_var: varianbce cube with drill-core property  


                        """
                        self.gravfield = gravfield
                        self.magfield = magfield
                        self.drillfield = drillfield
                        self.sensor_locations = sensor_locations
                        self.drilldata0 = drilldata0
                        # Normalize data
                        grav_mean, grav_std = self.gravfield.mean(), self.gravfield.std()
                        gravfield_norm = (self.gravfield - grav_mean) / grav_std
                        magn_mean, magn_std = self.magfield.mean(), self.magfield.std()
                        magfield_norm = (self.magfield - magn_mean) / magn_std
                        drill_mean, drill_std = self.drillfield.mean(), self.drillfield.std()
                        drillfield_norm = (self.drillfield - drill_mean) / drill_std
                        # Create kernel distance matrix
                        self.points3D = kernel.calcGridPoints3D((xNcube, yNcube, zNcube), (xvoxsize, yvoxsize, zvoxsize))
                        self.D2 = kernel.calcDistanceMatrix(self.points3D)
                        # Combine Data and Forward models
                        voxelpos_drill = np.vstack([self.xxx[self.drilldata0 != 0], self.yyy[self.drilldata0 != 0], self.zzz[self.drilldata0 != 0]])
                        # Combine data:
                        self.Fs3 = np.hstack((gravfield_norm, magfield_norm, drillfield_norm))
                        # Calculate Sensitivity matrix
                        Asens_grav, _ = A_sens(magneticField * 0., self.sensor_locations, self.Edges, 'grav')
                        Asens_mag, _ = A_sens(magneticField, self.sensor_locations, self.Edges, 'magn')
                        Asens_drill = A_drill(voxelpos_drill, self.voxelpos)
                        # Calculate total transformation matrix
                        Asens_g= np.vstack([Asens_grav, np.zeros_like(Asens_mag), np.zeros_like(Asens_drill)])
                        Asens_m= np.vstack([np.zeros_like(Asens_grav), Asens_mag, np.zeros_like(Asens_drill)])
                        Asens_d= np.vstack([np.zeros_like(Asens_grav), np.zeros_like(Asens_mag), Asens_drill])
                        self.Asens3 = np.hstack([Asens_g, Asens_m, Asens_d])
                        # Run actual GP inversion
                        self.mu_rec, self.cov_rec, self.logl = self.predict3(calclogl = True)
                        # Reshape results and multiply with stadnard deviatio 
                        results_rec = self.mu_rec.reshape(3, yNcube, xNcube, zNcube)
                        results_var = np.diag(self.cov_rec).reshape(3, yNcube, xNcube, zNcube)
                        #Cut off outer voxel rows
                        #results_rec = results_rec[:, 1:yNcube-1, 1:xNcube-1, :]
                        #results_var = results_var[:, 1:yNcube-1, 1:xNcube-1, :]
                        density_rec = results_rec[0] * grav_std  # Model represents deviation from the mean
                        density_var = results_var[0] * grav_std**2
                        magsus_rec = results_rec[1] * magn_std # Model represents deviation from the mean
                        magsus_var = results_var[1] * magn_std**2
                        drill_rec = results_rec[2] * drill_std  # Model represents deviation from the mean
                        drill_var = results_var[2] * drill_std**2
                        return density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var

    ### Methods

    ` def create_cubegeometry(self)`{.name .flex}

    :   ::: {.section .desc}
        Create Cube geometry and coordinates

        See settings.yaml for cube specifications
        :::

        Expand source code

            def create_cubegeometry(self):
                    """
                    Create Cube geometry and coordinates

                    See settings.yaml for cube specifications
                    """
                    # Create voxel Edge coordinates:
                    xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
                    yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
                    zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
                    #zedge = np.linspace(0.,Ncube/zfactor, Ncube/zfactor+1) * voxsize
                    xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
                    self.Edges = np.asarray([xEdges, yEdges, -zEdges])
                    # Create voxel coordinates
                    xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
                    ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
                    znew = zmax - np.arange(zvoxsize/2., zLcube + zvoxsize/2., zvoxsize)
                    xx ,yy = np.meshgrid(xnew,ynew)
                    self.xxx, self.yyy, self.zzz = np.meshgrid(xnew, ynew, znew)
                    self.voxelpos = np.vstack([self.xxx.flatten(), self.yyy.flatten(), self.zzz.flatten()])
                    return self.voxelpos

    ` def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0)`{.name .flex}

    :   ::: {.section .desc}
        Joint Inversion and Cubing of sensor data

        PARAMETERS

        :param gravfield: 2D gravitational survey data :param magfield:
        2D magnetic survey data :param drillfield: drilldata

        RETURN

        density\_rec: cube with mean density magsus\_rec: cube with mean
        magnetic susceptibility drill\_rec: cube with mean drill-core
        property\
        density\_var: varicance cube with density magsus\_var: variance
        cube with magnetic susceptibility\
        drill\_var: varianbce cube with drill-core property
        :::

        Expand source code

            def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0):
                    """
                    Joint Inversion and Cubing of sensor data

                    PARAMETERS

                    :param gravfield: 2D gravitational survey data
                    :param magfield: 2D magnetic survey data
                    :param drillfield: drilldata

                    RETURN

                    density_rec: cube with mean density 
                    magsus_rec: cube with mean magnetic susceptibility
                    drill_rec: cube with mean drill-core property  
                    density_var: varicance cube with density 
                    magsus_var: variance cube with magnetic susceptibility  
                    drill_var: varianbce cube with drill-core property  


                    """
                    self.gravfield = gravfield
                    self.magfield = magfield
                    self.drillfield = drillfield
                    self.sensor_locations = sensor_locations
                    self.drilldata0 = drilldata0
                    # Normalize data
                    grav_mean, grav_std = self.gravfield.mean(), self.gravfield.std()
                    gravfield_norm = (self.gravfield - grav_mean) / grav_std
                    magn_mean, magn_std = self.magfield.mean(), self.magfield.std()
                    magfield_norm = (self.magfield - magn_mean) / magn_std
                    drill_mean, drill_std = self.drillfield.mean(), self.drillfield.std()
                    drillfield_norm = (self.drillfield - drill_mean) / drill_std
                    # Create kernel distance matrix
                    self.points3D = kernel.calcGridPoints3D((xNcube, yNcube, zNcube), (xvoxsize, yvoxsize, zvoxsize))
                    self.D2 = kernel.calcDistanceMatrix(self.points3D)
                    # Combine Data and Forward models
                    voxelpos_drill = np.vstack([self.xxx[self.drilldata0 != 0], self.yyy[self.drilldata0 != 0], self.zzz[self.drilldata0 != 0]])
                    # Combine data:
                    self.Fs3 = np.hstack((gravfield_norm, magfield_norm, drillfield_norm))
                    # Calculate Sensitivity matrix
                    Asens_grav, _ = A_sens(magneticField * 0., self.sensor_locations, self.Edges, 'grav')
                    Asens_mag, _ = A_sens(magneticField, self.sensor_locations, self.Edges, 'magn')
                    Asens_drill = A_drill(voxelpos_drill, self.voxelpos)
                    # Calculate total transformation matrix
                    Asens_g= np.vstack([Asens_grav, np.zeros_like(Asens_mag), np.zeros_like(Asens_drill)])
                    Asens_m= np.vstack([np.zeros_like(Asens_grav), Asens_mag, np.zeros_like(Asens_drill)])
                    Asens_d= np.vstack([np.zeros_like(Asens_grav), np.zeros_like(Asens_mag), Asens_drill])
                    self.Asens3 = np.hstack([Asens_g, Asens_m, Asens_d])
                    # Run actual GP inversion
                    self.mu_rec, self.cov_rec, self.logl = self.predict3(calclogl = True)
                    # Reshape results and multiply with stadnard deviatio 
                    results_rec = self.mu_rec.reshape(3, yNcube, xNcube, zNcube)
                    results_var = np.diag(self.cov_rec).reshape(3, yNcube, xNcube, zNcube)
                    #Cut off outer voxel rows
                    #results_rec = results_rec[:, 1:yNcube-1, 1:xNcube-1, :]
                    #results_var = results_var[:, 1:yNcube-1, 1:xNcube-1, :]
                    density_rec = results_rec[0] * grav_std  # Model represents deviation from the mean
                    density_var = results_var[0] * grav_std**2
                    magsus_rec = results_rec[1] * magn_std # Model represents deviation from the mean
                    magsus_var = results_var[1] * magn_std**2
                    drill_rec = results_rec[2] * drill_std  # Model represents deviation from the mean
                    drill_var = results_var[2] * drill_std**2
                    return density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var

    ` def predict3(self, calclogl=False)`{.name .flex}

    :   ::: {.section .desc}
        Calculate mean, variance, and likelihood likelihood for GP with
        3x3 kernel block matrix

        PARAMETERS

        :param calclogl: if True calculate marginal GP loglikelihood, if
        False logl return is set to 0

        RETURN: mean covariance log-likelihood
        :::

        Expand source code

            def predict3(self, calclogl = False):
                    """
                    Calculate mean, variance, and likelihood likelihood for GP with 3x3 kernel block matrix 
                
                PARAMETERS

                :param calclogl: if True calculate marginal GP loglikelihood, if False logl return is set to 0

                RETURN: 
                mean
                covariance
                log-likelihood
                    """
                    # Calculate kernel 3x3 block matrix:
                    self.datastd = np.mean([np.nanstd(self.gravfield), np.nanstd(self.magfield), np.nanstd(self.drillfield)])
                    self.kcov = kernel.create_cov(self.D2, self.gp_length, self.coeffm, fkernel = kernelfunc)
                    #Ak = np.dot(self.Asens3, kcov) # shape (3*Nsensor, 3*Ncube^3)
                    yerr = np.hstack((self.gravfield * 0. + self.gp_sigma[0], self.magfield * 0. + self.gp_sigma[1], self.drillfield *0. + self.gp_sigma[2]))
                    #try:
                    AkA = np.dot(self.Asens3, np.dot(self.kcov, self.Asens3.T)) + np.diag(yerr**2)
                    #AkA = np.dot(Ak, self.Asens3.T) + np.diag(yerr**2) # shape(2*Nsensor, 2*Nsensor)
                    # Cholesky decomposition
                    try:
                            AkA_chol = cholesky(AkA, lower=True)
                    except:
                            print("Cholesky decompostion failed, AkA matrix i likely not positive semitive.")
                            print("Change GP parameter settings")
                            sys.exit(1)
                    usolve = solve_triangular(AkA_chol, self.Fs3, lower=True) #shape(2*Nsensor)
                    # Calculate Likelihood if option is set
                    if calclogl:
                            log_det_AkA = np.log(np.diag(AkA_chol)**2).sum()
                            n_log_2pi = xNcube * yNcube * zNcube * np.log(2 * np.pi)
                            logl = -0.5 * (np.dot(usolve, usolve) + log_det_AkA + n_log_2pi)
                    else:
                            logl = 0.
                    # predicted mean
                    self._V = solve_triangular(AkA_chol, np.dot(self.Asens3, self.kcov), lower=True) #shape(2*Nsensor, 2* Ncube^3)
                    mu = np.dot(self._V.T, usolve) #* _fstd + _fmean
                    # covariance
                    covar = self.kcov - np.dot(self._V.T, self._V)
                    #var = np.diagonal(covar) #* _fstd**2
                    #except:
                    #       print("Warning: GP Inversion failed!")
                    #       mu, covar, logl = np.zeros(3*xNcube*yNcube*zNcube), kcov * 1e9, np.inf 
                    return mu, covar, logl
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Classes](#header-classes) {#classes}

    -   #### `Inversion`

        -   `create_cubegeometry`
        -   `cubing`
        -   `predict3`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
::: {role="main"}
::: {#section-intro .section}
Framework for plotting 3D volumetric data (3D cube) as voxel grid.

Includes export of datacube as vtk file

Author: Sebastian Haan

Expand source code

    """
    Framework for plotting 3D volumetric data (3D cube) as voxel grid.

    Includes export of datacube as vtk file

    Author: Sebastian Haan

    """

    from __future__ import print_function
    import numpy as np
    import os
    import sys
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import LightSource
    import pyvista as pv # helper module for the Visualization Toolkit (VTK)
    #import plotly.plotly as py
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage import measure



    def skplot2(density, Nsize, drill = [], sensor = [], savefig = True, show = True, path_out = '', filename = 'density-drill-mesh.png'):
        """
        Makes volumetric 3D mesh plot at different density levels
        and plots vertical drillholes as well
        :params drill: x and y drillpositions
        """
        if savefig:
            if not os.path.exists(path_out):
                raise ValueError("Path not found: " + path_out)
            else:
                outfile = path_out + filename
        Nsize = np.asarray(Nsize)
        if len(drill) > 0:
            drill = np.asarray(drill)
            xdrill = drill[0]
            ydrill = drill[1]
        if len(sensor) > 0:
            sensor = np.asarray(sensor)
            x_sensor = sensor[0]
            y_sensor = sensor[1]
            z_sensor = sensor[2]
        colorscheme = cm.viridis
        fstd = density.std()
        fmean = density.mean()
        fmax = density.max()  
        if fmean < 0:
            sfaces = np.linspace(0+0.25*fstd, fmax*0.99, 5)
        else:
            sfaces = np.linspace(fmean+0.25*fstd, fmax*0.99, 5)
        cfaces = colorscheme(sfaces/sfaces.max())
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        # plot density
        sfact = sfaces.min()/sfaces.max()
        cfaces[:,3] = 0.2*sfaces/sfaces.max() # introduced in new library
        for i in range(len(sfaces)):
            if sfaces[i] > 0.:
                verts, faces, normals, values = measure.marching_cubes_lewiner(density, sfaces[i])
                mesh = Poly3DCollection(verts[faces])
                mesh.set_facecolor(cfaces[i])
                #mesh.set_alpha(0.2*sfaces[i]/sfaces.max())
                mesh.set_edgecolor(cfaces[i])
                ax.add_collection3d(mesh)
        # plot drillhole
        if len(drill) > 0:
            for i in range(len(xdrill)):
                ax.plot([xdrill[i],xdrill[i]],[ydrill[i],ydrill[i]],[0,density.shape[2]-1], 'k')
        # plot grav and mag sensors
        if len(sensor) > 0:
            ax.scatter(x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten(), c='r', marker='o')

        ax.set_xlabel("Y")
        ax.set_ylabel("X")
        ax.set_zlabel("-Z")

        ax.set_xlim(1, Nsize[0])  # a = 6 (times two for 2nd ellipsoid)
        ax.set_ylim(1, Nsize[1])  # b = 10
        ax.set_zlim(0, Nsize[2]-1)  # c = 16
        ax.invert_zaxis()
        ax.set_title("Drillholes: " + str(len(xdrill)))
        m = cm.ScalarMappable(cmap = colorscheme)
        m.set_array(sfaces)
        #cbar = plt.colorbar(m)
        plt.tight_layout()
        if savefig:
            plt.savefig(outfile)
        if show:
            plt.show()


    def skplot3(density, Nsize, drill = [], sensor = [], savefig = True, show = True, path_out = '', filename = 'density-drill-mesh.png'):
        """
        Makes columetric 3D mesh plot at different density levels
        and plots non-vertical drillholes
        :params drill: x, y, z drillpositions
        """
        if savefig:
            if not os.path.exists(path_out):
                raise ValueError("Path not found: " + path_out)
            else:
                outfile = path_out + filename
        Nsize = np.asarray(Nsize)
        if len(drill) > 0:
            drill = np.asarray(drill)
            xdrill = drill[0]
            ydrill = drill[1]
            zdrill = drill[2]
        if len(sensor) > 0:
            sensor = np.asarray(sensor)
            x_sensor = sensor[0]
            y_sensor = sensor[1]
            z_sensor = sensor[2]
        colorscheme = cm.viridis
        fstd = density.std()
        fmean = density.mean()
        fmax = np.percentile(density, 99)
        #fmax = density.max() 
        fmin = np.percentile(density, 1)
        if fmean < 0:
            sfaces = np.linspace(0+0.25*fstd, fmax, 5)
        else:
            sfaces = np.linspace(fmin+0.25*fstd, fmax, 5)
        cfaces = colorscheme((sfaces - fmin)/(sfaces - fmin).max())
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        # plot density
        sfact = sfaces.min()/sfaces.max()
        cfaces[:,3] = 0.2*sfaces/sfaces.max() # introduced in new library
        for i in range(len(sfaces)):
            if sfaces[i] > 0.:
                verts, faces, normals, values = measure.marching_cubes_lewiner(density, sfaces[i])
                mesh = Poly3DCollection(verts[faces])
                mesh.set_facecolor(cfaces[i])
                #mesh.set_alpha(0.2*sfaces[i]/sfaces.max())
                mesh.set_edgecolor(cfaces[i])
                ax.add_collection3d(mesh)
        # plot drillhole
        if len(drill) > 0:
            #ax.scatter(xdrill, ydrill, zdrill, c = 'k', marker = 'o')
            for i in range(len(xdrill)):
                ax.plot([xdrill[i, 0],xdrill[i, 1]],[ydrill[i, 0],ydrill[i, 1]],[zdrill[i, 0],zdrill[i, 1]], 'k')
        # plot grav and mag sensors
        if len(sensor) > 0:
            ax.scatter(x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten(), c='r', marker='o')

        ax.set_xlabel("Y")
        ax.set_ylabel("X")
        ax.set_zlabel("-Z")

        ax.set_xlim(1, Nsize[0])  # a = 6 (times two for 2nd ellipsoid)
        ax.set_ylim(1, Nsize[1])  # b = 10
        ax.set_zlim(0, Nsize[2]-1)  # c = 16
        ax.invert_zaxis()
        ax.set_title("Drillholes: " + str(len(xdrill)))
        m = cm.ScalarMappable(cmap = colorscheme)
        m.set_array(sfaces)
        #cbar = plt.colorbar(m)
        plt.tight_layout()
        if savefig:
            plt.savefig(outfile)
        if show:
            plt.show()


    def create_vtkcube(density, origin, voxelsize, fname):
        """
        Export Cube as VTK file (can be used in e.g. ParaView)
        and create a range of 3D cube plots with pyvista
        :param density: 3D cube in shape (xdim, ydim, zdim)
        :param origin: origin cooridnates of cube
        :param voxelsize: voxel sizes in (xsize, ysize, zsize)
        :param fname: path + filename for files
        """
        grid = pv.UniformGrid()
        grid.dimensions = np.array(density.shape) + 1
        grid.origin = origin
        grid.spacing = voxelsize
        grid.cell_arrays["values"] = density.flatten(order="F")
        grid.save(fname)
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Functions {#header-functions .section-title}
---------

` def create_vtkcube(density, origin, voxelsize, fname)`{.name .flex}

:   ::: {.section .desc}
    Export Cube as VTK file (can be used in e.g. ParaView) and create a
    range of 3D cube plots with pyvista :param density: 3D cube in shape
    (xdim, ydim, zdim) :param origin: origin cooridnates of cube :param
    voxelsize: voxel sizes in (xsize, ysize, zsize) :param fname: path +
    filename for files
    :::

    Expand source code

        def create_vtkcube(density, origin, voxelsize, fname):
            """
            Export Cube as VTK file (can be used in e.g. ParaView)
            and create a range of 3D cube plots with pyvista
            :param density: 3D cube in shape (xdim, ydim, zdim)
            :param origin: origin cooridnates of cube
            :param voxelsize: voxel sizes in (xsize, ysize, zsize)
            :param fname: path + filename for files
            """
            grid = pv.UniformGrid()
            grid.dimensions = np.array(density.shape) + 1
            grid.origin = origin
            grid.spacing = voxelsize
            grid.cell_arrays["values"] = density.flatten(order="F")
            grid.save(fname)

` def skplot2(density, Nsize, drill=[], sensor=[], savefig=True, show=True, path_out='', filename='density-drill-mesh.png')`{.name .flex}

:   ::: {.section .desc}
    Makes volumetric 3D mesh plot at different density levels and plots
    vertical drillholes as well :params drill: x and y drillpositions
    :::

    Expand source code

        def skplot2(density, Nsize, drill = [], sensor = [], savefig = True, show = True, path_out = '', filename = 'density-drill-mesh.png'):
            """
            Makes volumetric 3D mesh plot at different density levels
            and plots vertical drillholes as well
            :params drill: x and y drillpositions
            """
            if savefig:
                if not os.path.exists(path_out):
                    raise ValueError("Path not found: " + path_out)
                else:
                    outfile = path_out + filename
            Nsize = np.asarray(Nsize)
            if len(drill) > 0:
                drill = np.asarray(drill)
                xdrill = drill[0]
                ydrill = drill[1]
            if len(sensor) > 0:
                sensor = np.asarray(sensor)
                x_sensor = sensor[0]
                y_sensor = sensor[1]
                z_sensor = sensor[2]
            colorscheme = cm.viridis
            fstd = density.std()
            fmean = density.mean()
            fmax = density.max()  
            if fmean < 0:
                sfaces = np.linspace(0+0.25*fstd, fmax*0.99, 5)
            else:
                sfaces = np.linspace(fmean+0.25*fstd, fmax*0.99, 5)
            cfaces = colorscheme(sfaces/sfaces.max())
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            # plot density
            sfact = sfaces.min()/sfaces.max()
            cfaces[:,3] = 0.2*sfaces/sfaces.max() # introduced in new library
            for i in range(len(sfaces)):
                if sfaces[i] > 0.:
                    verts, faces, normals, values = measure.marching_cubes_lewiner(density, sfaces[i])
                    mesh = Poly3DCollection(verts[faces])
                    mesh.set_facecolor(cfaces[i])
                    #mesh.set_alpha(0.2*sfaces[i]/sfaces.max())
                    mesh.set_edgecolor(cfaces[i])
                    ax.add_collection3d(mesh)
            # plot drillhole
            if len(drill) > 0:
                for i in range(len(xdrill)):
                    ax.plot([xdrill[i],xdrill[i]],[ydrill[i],ydrill[i]],[0,density.shape[2]-1], 'k')
            # plot grav and mag sensors
            if len(sensor) > 0:
                ax.scatter(x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten(), c='r', marker='o')

            ax.set_xlabel("Y")
            ax.set_ylabel("X")
            ax.set_zlabel("-Z")

            ax.set_xlim(1, Nsize[0])  # a = 6 (times two for 2nd ellipsoid)
            ax.set_ylim(1, Nsize[1])  # b = 10
            ax.set_zlim(0, Nsize[2]-1)  # c = 16
            ax.invert_zaxis()
            ax.set_title("Drillholes: " + str(len(xdrill)))
            m = cm.ScalarMappable(cmap = colorscheme)
            m.set_array(sfaces)
            #cbar = plt.colorbar(m)
            plt.tight_layout()
            if savefig:
                plt.savefig(outfile)
            if show:
                plt.show()

` def skplot3(density, Nsize, drill=[], sensor=[], savefig=True, show=True, path_out='', filename='density-drill-mesh.png')`{.name .flex}

:   ::: {.section .desc}
    Makes columetric 3D mesh plot at different density levels and plots
    non-vertical drillholes :params drill: x, y, z drillpositions
    :::

    Expand source code

        def skplot3(density, Nsize, drill = [], sensor = [], savefig = True, show = True, path_out = '', filename = 'density-drill-mesh.png'):
            """
            Makes columetric 3D mesh plot at different density levels
            and plots non-vertical drillholes
            :params drill: x, y, z drillpositions
            """
            if savefig:
                if not os.path.exists(path_out):
                    raise ValueError("Path not found: " + path_out)
                else:
                    outfile = path_out + filename
            Nsize = np.asarray(Nsize)
            if len(drill) > 0:
                drill = np.asarray(drill)
                xdrill = drill[0]
                ydrill = drill[1]
                zdrill = drill[2]
            if len(sensor) > 0:
                sensor = np.asarray(sensor)
                x_sensor = sensor[0]
                y_sensor = sensor[1]
                z_sensor = sensor[2]
            colorscheme = cm.viridis
            fstd = density.std()
            fmean = density.mean()
            fmax = np.percentile(density, 99)
            #fmax = density.max() 
            fmin = np.percentile(density, 1)
            if fmean < 0:
                sfaces = np.linspace(0+0.25*fstd, fmax, 5)
            else:
                sfaces = np.linspace(fmin+0.25*fstd, fmax, 5)
            cfaces = colorscheme((sfaces - fmin)/(sfaces - fmin).max())
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            # plot density
            sfact = sfaces.min()/sfaces.max()
            cfaces[:,3] = 0.2*sfaces/sfaces.max() # introduced in new library
            for i in range(len(sfaces)):
                if sfaces[i] > 0.:
                    verts, faces, normals, values = measure.marching_cubes_lewiner(density, sfaces[i])
                    mesh = Poly3DCollection(verts[faces])
                    mesh.set_facecolor(cfaces[i])
                    #mesh.set_alpha(0.2*sfaces[i]/sfaces.max())
                    mesh.set_edgecolor(cfaces[i])
                    ax.add_collection3d(mesh)
            # plot drillhole
            if len(drill) > 0:
                #ax.scatter(xdrill, ydrill, zdrill, c = 'k', marker = 'o')
                for i in range(len(xdrill)):
                    ax.plot([xdrill[i, 0],xdrill[i, 1]],[ydrill[i, 0],ydrill[i, 1]],[zdrill[i, 0],zdrill[i, 1]], 'k')
            # plot grav and mag sensors
            if len(sensor) > 0:
                ax.scatter(x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten(), c='r', marker='o')

            ax.set_xlabel("Y")
            ax.set_ylabel("X")
            ax.set_zlabel("-Z")

            ax.set_xlim(1, Nsize[0])  # a = 6 (times two for 2nd ellipsoid)
            ax.set_ylim(1, Nsize[1])  # b = 10
            ax.set_zlim(0, Nsize[2]-1)  # c = 16
            ax.invert_zaxis()
            ax.set_title("Drillholes: " + str(len(xdrill)))
            m = cm.ScalarMappable(cmap = colorscheme)
            m.set_array(sfaces)
            #cbar = plt.colorbar(m)
            plt.tight_layout()
            if savefig:
                plt.savefig(outfile)
            if show:
                plt.show()
:::

::: {.section}
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Functions](#header-functions) {#functions}

    -   `create_vtkcube`
    -   `skplot2`
    -   `skplot3`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
::: {role="main"}
::: {#section-intro .section}
Expand source code

    # Utility functions for coordinate conversion and reading ion drillcore and survey data

    import numpy as np
    import pandas as pd
    import rasterio
    import os
    from scipy.ndimage.interpolation import zoom
    from config_loader import *

    def spherical2cartes(x0, y0, z0, phi, theta, r):
            """ Conversion from spherical coordinates to cartesian

            PARAMETER
            param x0, y0,z0: coordinates of origin
            param phi: azimuthal angle
            param theta: polar angle
            param r: radial length

            RETURN
            x ,y, z coordinates
            """
            x = x0 + r * np.sin(theta) * np.cos(phi)
            y = y0 + r * np.sin(theta) * np.sin(phi)
            z = z0 + r * np.cos(theta)
            return x, y, z


    def cartes2spherial(x0, y0, z0, x1, y1, z1):
            """ Conversion from cartesian coordinates to spehrical

            PARAMETER
            param x0, y0,z0: coordinates of origin
            param x1, y1,z1: coordinates of end point

            RETURN
            radius, polar angle, azimuthal angle
            """
            r = np.sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)
            theta = np.arccos((z1-z0)/r)
            phi = np.arctan2((y1-y0),(x1-x0))
            return r, theta, phi


    def align_drill2(coord, data):
            """Convert drill-core data in Model Cube shape with iteratinon over all voxels

            PARAMETER

            param coord: xyz Coordinates of drillcore, shape (N_drill, 3)
            param data: 1 dim data for drillcore, shape (N_drill)

            RETURN

            drillcore voxel data cube
            """
            dx = xvoxsize
            dy = yvoxsize
            dz = zvoxsize
            data = np.asarray(data)
            res = np.zeros_like(xxx)
            for ix in range(xxx.shape[0]):
                    for iy in range(xxx.shape[1]):
                            for iz in range(xxx.shape[2]):
                                    sel = np.where(((xxx[ix,iy,iz] - dx) <= coord[:, 0]) & (coord[:, 0] < (xxx[ix,iy,iz] + dx))
                                            & ((yyy[ix,iy,iz] - dy) <= coord[:, 1]) & (coord[:, 1] < (yyy[ix,iy,iz] + dy))
                                            & ((zzz[ix,iy,iz] - dz) <= coord[:, 2]) & (coord[:, 2] < (zzz[ix,iy,iz] + dz)))
                                    if np.size(sel) > 0:
                                            #print(sel)
                                            m = np.nanmean(data[sel])
                                            if np.isfinite(m):
                                                    res[ix,iy,iz] = m
            return res


    def normalize(x):
            """ Normalise x
            :param x: input 1D array to normalise
            
            Return
            normalized array
            """
            if abs(x.max() - x.min()) > 0:
                    norm = (x- x.min()) / (x.max() - x.min())  
            else:
                    norm = x / x.max()
            return norm


    def readcsv_drill(fname, prop_names, pos_names, cartesian = True):
            """ Read csv file with drill-core data and ectract data for density, magsus, and mineral content
            Drill-core data is converted in casrtesian xyz coordinates if only xy position, length and azimuth and dip angle available

            PARAMETER

            :param fname: Path and filename of csv file
            :param prop_names: string array of property names in header of file, e.g ['density', 'magsus', 'mineral']
            :param pos_names: string array of position names in header ['x', 'y', 'z'] 
                    or ['East', 'North', 'Elev', DepthFrom', 'DepthTo', 'Azimuth', 'Dip'] with Easting and Northings in meter and Azimuth and Drip in degrees
            :param cartesian: if 'True' (Default), the drill positions are provided in xyz fromat in meters, 
                    if 'False' the drill positions are cualcuated from drill depth, azimuth and dip


            RETURN

            returns: pandas array with positions and drill properties
            """

            names = prop_names + pos_names
            drill = pd.read_csv(fname_drilldata, names = names)
            if cartesian:
                    xdrill, ydrill, zdrill = drill[pos_names[0]].values, drill[pos_names[1]].values, drill[pos_names[2]].values
            else:
                    xdrill,ydrill,zdrill = spherical2cartes(drill.East.values, drill.North.values, drill.Elev.values, 
                            drill.Azimuth.values * np.pi/180., (90 - drill.Dip.values)* np.pi/180., 0.5*(drill.DepthFrom.values + drill.DepthTo.values))
            data = pd.DataFrame()
            data['x'] = xdrill
            data['y'] = ydrill
            data['z'] = zdrill
            data['density'] = drill[prop_names[0]].values
            data['magsus'] = drill[prop_names[1]].values
            data['mineral'] = drill[prop_names[2]].values
            return data


    def readraster_survey(fname, pixres = None, clipext = None):
            """ Read geotif rasterfile with one band per property

            PARAMETER
            :param fname: Path and filename of csv file
            :param pixres: desired pixel resolution in meters (assuming same resolution for x and y) 
            :param clipext: provide boundary box [xmin,ymin,xmax,ymax] to crop data to extent


            RETURN
            returns: numpy array with survey data and x,y pixelcoordinates
            """

            img = rasterio.open(fname)
            bounds = img.bounds #format: [xmin, ymin, xmax, ymax]
            img_shape = img.shape
            # Native pixel resolutiojn:
            ypixres, xpixres = img.res
            # if ypixsize == xpixsize:
            #       print("WARNING, Input raster pixelsize (resolution) is not the same in vertical and horizontal direction")
            # Read in data into numpy array:
            array = img.read(1)
            # Regrid array to desired 
            if pixres is not None:
                    assert ypixres == xpixres, " Raster pixelsize (resolution) is not the same in vertical and horizontal direction"
                    array  = zoom(array, ypixsize / pixres)
                    xpixres = ypixres = pixres
            # Define pixel coordinates of original image
            xx, yy = np.meshgrid(np.linspace(bounds[0] + 0.5*xpixres, bounds[2] - 0.5*xpixres, array.shape[1]), 
                                    np.linspace(bounds[1] + 0.5*ypixres, bounds[3] - 0.5*ypixres, array.shape[0]))
            if clipext is not None:
                    clipext = np.asarray(clipext)
                    # check if clip extent is within boudning box';
                    if ((bounds[0]<= clipext[0] <= bounds[2]) &  
                            (bounds[1]<= clipext[1] <= bounds[3]) &
                            (bounds[0]<= clipext[2] <= bounds[2]) &
                            (bounds[1]<= clipext[3] <= bounds[3])):
                            # Clip array:
                            array = array[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
                             & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
                            xx2 = xx[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
                             & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
                            yy2 = yy[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
                             & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
                            array = array.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
                            xx2 = xx2.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
                            yy2 = yy2.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
                    else:
                            print('WARNING: Clip extent exceeds image boundary!... No clipping is applied')
            else: 
                    xx2, yy2 = xx, yy
            return array, np.asarray([xx2, yy2])
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Functions {#header-functions .section-title}
---------

` def align_drill2(coord, data)`{.name .flex}

:   ::: {.section .desc}
    Convert drill-core data in Model Cube shape with iteratinon over all
    voxels

    PARAMETER

    param coord: xyz Coordinates of drillcore, shape (N\_drill, 3) param
    data: 1 dim data for drillcore, shape (N\_drill)

    RETURN

    drillcore voxel data cube
    :::

    Expand source code

        def align_drill2(coord, data):
                """Convert drill-core data in Model Cube shape with iteratinon over all voxels

                PARAMETER

                param coord: xyz Coordinates of drillcore, shape (N_drill, 3)
                param data: 1 dim data for drillcore, shape (N_drill)

                RETURN

                drillcore voxel data cube
                """
                dx = xvoxsize
                dy = yvoxsize
                dz = zvoxsize
                data = np.asarray(data)
                res = np.zeros_like(xxx)
                for ix in range(xxx.shape[0]):
                        for iy in range(xxx.shape[1]):
                                for iz in range(xxx.shape[2]):
                                        sel = np.where(((xxx[ix,iy,iz] - dx) <= coord[:, 0]) & (coord[:, 0] < (xxx[ix,iy,iz] + dx))
                                                & ((yyy[ix,iy,iz] - dy) <= coord[:, 1]) & (coord[:, 1] < (yyy[ix,iy,iz] + dy))
                                                & ((zzz[ix,iy,iz] - dz) <= coord[:, 2]) & (coord[:, 2] < (zzz[ix,iy,iz] + dz)))
                                        if np.size(sel) > 0:
                                                #print(sel)
                                                m = np.nanmean(data[sel])
                                                if np.isfinite(m):
                                                        res[ix,iy,iz] = m
                return res

` def cartes2spherial(x0, y0, z0, x1, y1, z1)`{.name .flex}

:   ::: {.section .desc}
    Conversion from cartesian coordinates to spehrical

    PARAMETER param x0, y0,z0: coordinates of origin param x1, y1,z1:
    coordinates of end point

    RETURN radius, polar angle, azimuthal angle
    :::

    Expand source code

        def cartes2spherial(x0, y0, z0, x1, y1, z1):
                """ Conversion from cartesian coordinates to spehrical

                PARAMETER
                param x0, y0,z0: coordinates of origin
                param x1, y1,z1: coordinates of end point

                RETURN
                radius, polar angle, azimuthal angle
                """
                r = np.sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)
                theta = np.arccos((z1-z0)/r)
                phi = np.arctan2((y1-y0),(x1-x0))
                return r, theta, phi

` def normalize(x)`{.name .flex}

:   ::: {.section .desc}
    Normalise x :param x: input 1D array to normalise

    Return normalized array
    :::

    Expand source code

        def normalize(x):
                """ Normalise x
                :param x: input 1D array to normalise
                
                Return
                normalized array
                """
                if abs(x.max() - x.min()) > 0:
                        norm = (x- x.min()) / (x.max() - x.min())  
                else:
                        norm = x / x.max()
                return norm

` def readcsv_drill(fname, prop_names, pos_names, cartesian=True)`{.name .flex}

:   ::: {.section .desc}
    Read csv file with drill-core data and ectract data for density,
    magsus, and mineral content Drill-core data is converted in
    casrtesian xyz coordinates if only xy position, length and azimuth
    and dip angle available

    PARAMETER

    :param fname: Path and filename of csv file :param prop\_names:
    string array of property names in header of file, e.g \[\'density\',
    \'magsus\', \'mineral\'\] :param pos\_names: string array of
    position names in header \[\'x\', \'y\', \'z\'\] or \[\'East\',
    \'North\', \'Elev\', DepthFrom\', \'DepthTo\', \'Azimuth\',
    \'Dip\'\] with Easting and Northings in meter and Azimuth and Drip
    in degrees :param cartesian: if \'True\' (Default), the drill
    positions are provided in xyz fromat in meters, if \'False\' the
    drill positions are cualcuated from drill depth, azimuth and dip

    RETURN

    returns: pandas array with positions and drill properties
    :::

    Expand source code

        def readcsv_drill(fname, prop_names, pos_names, cartesian = True):
                """ Read csv file with drill-core data and ectract data for density, magsus, and mineral content
                Drill-core data is converted in casrtesian xyz coordinates if only xy position, length and azimuth and dip angle available

                PARAMETER

                :param fname: Path and filename of csv file
                :param prop_names: string array of property names in header of file, e.g ['density', 'magsus', 'mineral']
                :param pos_names: string array of position names in header ['x', 'y', 'z'] 
                        or ['East', 'North', 'Elev', DepthFrom', 'DepthTo', 'Azimuth', 'Dip'] with Easting and Northings in meter and Azimuth and Drip in degrees
                :param cartesian: if 'True' (Default), the drill positions are provided in xyz fromat in meters, 
                        if 'False' the drill positions are cualcuated from drill depth, azimuth and dip


                RETURN

                returns: pandas array with positions and drill properties
                """

                names = prop_names + pos_names
                drill = pd.read_csv(fname_drilldata, names = names)
                if cartesian:
                        xdrill, ydrill, zdrill = drill[pos_names[0]].values, drill[pos_names[1]].values, drill[pos_names[2]].values
                else:
                        xdrill,ydrill,zdrill = spherical2cartes(drill.East.values, drill.North.values, drill.Elev.values, 
                                drill.Azimuth.values * np.pi/180., (90 - drill.Dip.values)* np.pi/180., 0.5*(drill.DepthFrom.values + drill.DepthTo.values))
                data = pd.DataFrame()
                data['x'] = xdrill
                data['y'] = ydrill
                data['z'] = zdrill
                data['density'] = drill[prop_names[0]].values
                data['magsus'] = drill[prop_names[1]].values
                data['mineral'] = drill[prop_names[2]].values
                return data

` def readraster_survey(fname, pixres=None, clipext=None)`{.name .flex}

:   ::: {.section .desc}
    Read geotif rasterfile with one band per property

    PARAMETER :param fname: Path and filename of csv file :param pixres:
    desired pixel resolution in meters (assuming same resolution for x
    and y) :param clipext: provide boundary box \[xmin,ymin,xmax,ymax\]
    to crop data to extent

    RETURN returns: numpy array with survey data and x,y
    pixelcoordinates
    :::

    Expand source code

        def readraster_survey(fname, pixres = None, clipext = None):
                """ Read geotif rasterfile with one band per property

                PARAMETER
                :param fname: Path and filename of csv file
                :param pixres: desired pixel resolution in meters (assuming same resolution for x and y) 
                :param clipext: provide boundary box [xmin,ymin,xmax,ymax] to crop data to extent


                RETURN
                returns: numpy array with survey data and x,y pixelcoordinates
                """

                img = rasterio.open(fname)
                bounds = img.bounds #format: [xmin, ymin, xmax, ymax]
                img_shape = img.shape
                # Native pixel resolutiojn:
                ypixres, xpixres = img.res
                # if ypixsize == xpixsize:
                #       print("WARNING, Input raster pixelsize (resolution) is not the same in vertical and horizontal direction")
                # Read in data into numpy array:
                array = img.read(1)
                # Regrid array to desired 
                if pixres is not None:
                        assert ypixres == xpixres, " Raster pixelsize (resolution) is not the same in vertical and horizontal direction"
                        array  = zoom(array, ypixsize / pixres)
                        xpixres = ypixres = pixres
                # Define pixel coordinates of original image
                xx, yy = np.meshgrid(np.linspace(bounds[0] + 0.5*xpixres, bounds[2] - 0.5*xpixres, array.shape[1]), 
                                        np.linspace(bounds[1] + 0.5*ypixres, bounds[3] - 0.5*ypixres, array.shape[0]))
                if clipext is not None:
                        clipext = np.asarray(clipext)
                        # check if clip extent is within boudning box';
                        if ((bounds[0]<= clipext[0] <= bounds[2]) &  
                                (bounds[1]<= clipext[1] <= bounds[3]) &
                                (bounds[0]<= clipext[2] <= bounds[2]) &
                                (bounds[1]<= clipext[3] <= bounds[3])):
                                # Clip array:
                                array = array[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
                                 & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
                                xx2 = xx[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
                                 & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
                                yy2 = yy[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
                                 & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
                                array = array.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
                                xx2 = xx2.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
                                yy2 = yy2.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
                        else:
                                print('WARNING: Clip extent exceeds image boundary!... No clipping is applied')
                else: 
                        xx2, yy2 = xx, yy
                return array, np.asarray([xx2, yy2])

` def spherical2cartes(x0, y0, z0, phi, theta, r)`{.name .flex}

:   ::: {.section .desc}
    Conversion from spherical coordinates to cartesian

    PARAMETER param x0, y0,z0: coordinates of origin param phi:
    azimuthal angle param theta: polar angle param r: radial length

    RETURN x ,y, z coordinates
    :::

    Expand source code

        def spherical2cartes(x0, y0, z0, phi, theta, r):
                """ Conversion from spherical coordinates to cartesian

                PARAMETER
                param x0, y0,z0: coordinates of origin
                param phi: azimuthal angle
                param theta: polar angle
                param r: radial length

                RETURN
                x ,y, z coordinates
                """
                x = x0 + r * np.sin(theta) * np.cos(phi)
                y = y0 + r * np.sin(theta) * np.sin(phi)
                z = z0 + r * np.cos(theta)
                return x, y, z
:::

::: {.section}
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Functions](#header-functions) {#functions}

    -   `align_drill2`
    -   `cartes2spherial`
    -   `normalize`
    -   `readcsv_drill`
    -   `readraster_survey`
    -   `spherical2cartes`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
::: {role="main"}
::: {#section-intro .section}
Calcuation of forward models of gravity and magnetic sensitivity as well
as drill core transfromation.

The forward models transform the localized measurement of a remote
sensor grid into a 3D representation of geophysical properties of a
region, here gravity and magnetic forward models: The gravity forward
model is defined by using Li\'s tractable approximation for a 3-D field
of constant density prisms ( Li and Oldenbur, 3D-inversion of gravity
data, 1998}) and can be determined analytically. The induced magnetic
field calculation uses Li\'s tractable approximation for a 3-D field of
prisms of constant magnetic susceptibility, which depends on the
magnetic mineral content below the surface and is measured by the
response of their magnetic dipoles induced by the Earth\'s magnetic
field. The joint GP inversion takes into account a covariance that
exists between density and magnetic susceptibility.

Author: Sebastian Haan

Expand source code

    """
    Calcuation of forward models of gravity and magnetic sensitivity as well as drill core transfromation.

    The forward models transform the localized measurement of a remote sensor grid into a 3D representation 
    of geophysical properties of a region, here  gravity and magnetic forward models: 
    The gravity forward model is defined by using Li's tractable approximation for a 3-D field 
    of constant density prisms ( Li and Oldenbur, 3D-inversion of gravity data, 1998}) 
    and can be determined analytically. The induced magnetic field calculation uses Li's tractable approximation 
    for a 3-D field of prisms of constant magnetic susceptibility, 
    which depends on the magnetic mineral content below the surface and is measured 
    by the response of their magnetic dipoles induced by the Earth's magnetic field. 
    The joint GP inversion takes into account a covariance that exists between density and magnetic susceptibility.

    Author: Sebastian Haan
    """

    import numpy as np
    from config_loader import *



    def A_sens(magneticField, locations, Edges, func):
            """
            calculate gravity and magnetic forward model matrix

            PARAMETER

            :param magneicField: ambient magnetic field
            :param func: 'grav' for gravity, 'magn' for magnetic

            RETURNS

            Forward model
            eZ
            """
            Edges = np.asarray(Edges)
            xEdges = Edges[0]
            yEdges = Edges[1]
            zEdges = Edges[2]
            bx = magneticField[0]
            by = magneticField[1]
            bz = magneticField[2]
            nprism = xNcube* yNcube * zNcube
            nedge = (xNcube+1)* (yNcube+1) * (zNcube+1)
            sens = np.zeros((1 * xNcube * yNcube, nprism))
            result_ez = np.zeros((1 * xNcube * yNcube, nedge))

            # For each sensor location
            for n in range(1 * xNcube * yNcube):
                    x0 = xEdges - locations[n, 0]
                    y0 = yEdges - locations[n, 1]
                    z0 = zEdges - locations[n, 2] # z-axis already inverted

                    # Lazy edge padding for both grav and mag
                    aLongWay = 1e6 # Metres, chosen as in Obsedian
                    x0[0] -= aLongWay
                    y0[0] -= aLongWay
                    x0[-1] += aLongWay
                    y0[-1] += aLongWay

                    #Precompute eZ for each point
                    if func == 'grav':
                            eZ = grav_func(x0, y0, z0) # shape (Ncube+1,Ncube+1,Ncube+1)
                    elif func == 'magn':
                            eZ = magn_func(x0, y0, z0, bx, by, bz)
                    else:
                            print('function not supported')
                    result_ez[n,:] = eZ.reshape((xNcube+1)*(yNcube+1)*(zNcube+1))

                    # Compute the sensitivities     
                    idx = 0
                    for i in range(xEdges.shape[0]-1):
                        for j in range(yEdges.shape[1]-1):
                            for k in range(zEdges.shape[2]-1):
                                    sens[n, idx] = -((eZ[i + 1, j + 1, k + 1] - eZ[i + 1, j + 1, k] - eZ[i + 1, j, k + 1] + eZ[i + 1, j, k]) 
                                            - (eZ[i, j + 1, k + 1] - eZ[i, j + 1, k] - eZ[i, j, k + 1] + eZ[i, j, k]))
                                    idx += 1

            if func == 'grav':
                    sens = c_MILLIGALS_UNITS * sens / fcor_grav
            if func == 'magn':
                    sens = sens / fcor_mag

            return sens, result_ez


    def grav_func(x, y, z):
            """
            Compute the vertical component of gravity
            Computes the sensitivity for a particular point in the gravity

            PARAMETER

            :param x: x coordinate of the position.
        :param y: y coordinate of the position.
        :param z: z coordinate of the position.
            """
            eps = 1e-9 
            r = np.sqrt(x**2 + y**2 + z**2)
            func =  x * np.log(y + r) + y * np.log(x + r) - z * np.arctan((x * y) / (z * r + eps))
            return func


    def magn_func(x,y,z,bx,by,bz):
            """
            Compute magnetic forward model
            Calculating the magnetic sensitivity at a particular position relative to the origin.

            PARAMETER

        :param x: x coordinate of the position.
        :param y: y coordinate of the position.
        :param z: z coordinate of the position.
        :param bx: The magnetic field in x-direction at this position
        :param by: The magnetic field in y-direction at this position
        :param bz: The magnetic field in z-direction at this position
            """
            r = np.sqrt(x**2 + y**2 + z**2)
            #Compute the normalisation factor for the magnetic field
            normB = np.sqrt(bx * bx + by * by + bz * bz)
            func = 1./ normB * ((2. * by * bz * np.log(x + r)) + (2. * bz * bx * np.log(y + r)) + (2. * by * bx * np.log(z + r))
                + (bz * bz - by * by) * np.arctan((x * z) / (y * r)) + (bz * bz - bx * bx) * np.arctan((y * z) / (x * r))) 
            res = -func
            return res


    def A_drill(loc, voxelpos):
            """
            Transform the voxel cube into filter matrix for drill hole with sensitivity 1

            PARAMETER
            
            :param loc: x,y,z drillcore coordinates in shape (Ndrillcore, 3)
            :param voxelpos: x,y,z coordinates in shape (3, Nvoxel)
            """
            x = voxelpos[0].flatten() #- 1
            y = voxelpos[1].flatten() #- 1
            z = voxelpos[2].flatten() #- 1
            sens = np.zeros((loc.shape[1], xNcube * yNcube * zNcube))
            for i in range(loc.shape[1]):
                    coord = loc[:,i]#.astype(int)
                    sel = np.where((x == coord[0]) & (y == coord[1]) & (z == coord[2]))
                    sens[i, sel] = 1 
            return sens
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Functions {#header-functions .section-title}
---------

` def A_drill(loc, voxelpos)`{.name .flex}

:   ::: {.section .desc}
    Transform the voxel cube into filter matrix for drill hole with
    sensitivity 1

    PARAMETER

    :param loc: x,y,z drillcore coordinates in shape (Ndrillcore, 3)
    :param voxelpos: x,y,z coordinates in shape (3, Nvoxel)
    :::

    Expand source code

        def A_drill(loc, voxelpos):
                """
                Transform the voxel cube into filter matrix for drill hole with sensitivity 1

                PARAMETER
                
                :param loc: x,y,z drillcore coordinates in shape (Ndrillcore, 3)
                :param voxelpos: x,y,z coordinates in shape (3, Nvoxel)
                """
                x = voxelpos[0].flatten() #- 1
                y = voxelpos[1].flatten() #- 1
                z = voxelpos[2].flatten() #- 1
                sens = np.zeros((loc.shape[1], xNcube * yNcube * zNcube))
                for i in range(loc.shape[1]):
                        coord = loc[:,i]#.astype(int)
                        sel = np.where((x == coord[0]) & (y == coord[1]) & (z == coord[2]))
                        sens[i, sel] = 1 
                return sens

` def A_sens(magneticField, locations, Edges, func)`{.name .flex}

:   ::: {.section .desc}
    calculate gravity and magnetic forward model matrix

    PARAMETER

    :param magneicField: ambient magnetic field :param func: \'grav\'
    for gravity, \'magn\' for magnetic

    RETURNS

    Forward model eZ
    :::

    Expand source code

        def A_sens(magneticField, locations, Edges, func):
                """
                calculate gravity and magnetic forward model matrix

                PARAMETER

                :param magneicField: ambient magnetic field
                :param func: 'grav' for gravity, 'magn' for magnetic

                RETURNS

                Forward model
                eZ
                """
                Edges = np.asarray(Edges)
                xEdges = Edges[0]
                yEdges = Edges[1]
                zEdges = Edges[2]
                bx = magneticField[0]
                by = magneticField[1]
                bz = magneticField[2]
                nprism = xNcube* yNcube * zNcube
                nedge = (xNcube+1)* (yNcube+1) * (zNcube+1)
                sens = np.zeros((1 * xNcube * yNcube, nprism))
                result_ez = np.zeros((1 * xNcube * yNcube, nedge))

                # For each sensor location
                for n in range(1 * xNcube * yNcube):
                        x0 = xEdges - locations[n, 0]
                        y0 = yEdges - locations[n, 1]
                        z0 = zEdges - locations[n, 2] # z-axis already inverted

                        # Lazy edge padding for both grav and mag
                        aLongWay = 1e6 # Metres, chosen as in Obsedian
                        x0[0] -= aLongWay
                        y0[0] -= aLongWay
                        x0[-1] += aLongWay
                        y0[-1] += aLongWay

                        #Precompute eZ for each point
                        if func == 'grav':
                                eZ = grav_func(x0, y0, z0) # shape (Ncube+1,Ncube+1,Ncube+1)
                        elif func == 'magn':
                                eZ = magn_func(x0, y0, z0, bx, by, bz)
                        else:
                                print('function not supported')
                        result_ez[n,:] = eZ.reshape((xNcube+1)*(yNcube+1)*(zNcube+1))

                        # Compute the sensitivities     
                        idx = 0
                        for i in range(xEdges.shape[0]-1):
                            for j in range(yEdges.shape[1]-1):
                                for k in range(zEdges.shape[2]-1):
                                        sens[n, idx] = -((eZ[i + 1, j + 1, k + 1] - eZ[i + 1, j + 1, k] - eZ[i + 1, j, k + 1] + eZ[i + 1, j, k]) 
                                                - (eZ[i, j + 1, k + 1] - eZ[i, j + 1, k] - eZ[i, j, k + 1] + eZ[i, j, k]))
                                        idx += 1

                if func == 'grav':
                        sens = c_MILLIGALS_UNITS * sens / fcor_grav
                if func == 'magn':
                        sens = sens / fcor_mag

                return sens, result_ez

` def grav_func(x, y, z)`{.name .flex}

:   ::: {.section .desc}
    Compute the vertical component of gravity Computes the sensitivity
    for a particular point in the gravity

        PARAMETER

        :param x: x coordinate of the position.

    :param y: y coordinate of the position. :param z: z coordinate of
    the position.
    :::

    Expand source code

        def grav_func(x, y, z):
                """
                Compute the vertical component of gravity
                Computes the sensitivity for a particular point in the gravity

                PARAMETER

                :param x: x coordinate of the position.
            :param y: y coordinate of the position.
            :param z: z coordinate of the position.
                """
                eps = 1e-9 
                r = np.sqrt(x**2 + y**2 + z**2)
                func =  x * np.log(y + r) + y * np.log(x + r) - z * np.arctan((x * y) / (z * r + eps))
                return func

` def magn_func(x, y, z, bx, by, bz)`{.name .flex}

:   ::: {.section .desc}
    Compute magnetic forward model Calculating the magnetic sensitivity
    at a particular position relative to the origin.

        PARAMETER

    :param x: x coordinate of the position. :param y: y coordinate of
    the position. :param z: z coordinate of the position. :param bx: The
    magnetic field in x-direction at this position :param by: The
    magnetic field in y-direction at this position :param bz: The
    magnetic field in z-direction at this position
    :::

    Expand source code

        def magn_func(x,y,z,bx,by,bz):
                """
                Compute magnetic forward model
                Calculating the magnetic sensitivity at a particular position relative to the origin.

                PARAMETER

            :param x: x coordinate of the position.
            :param y: y coordinate of the position.
            :param z: z coordinate of the position.
            :param bx: The magnetic field in x-direction at this position
            :param by: The magnetic field in y-direction at this position
            :param bz: The magnetic field in z-direction at this position
                """
                r = np.sqrt(x**2 + y**2 + z**2)
                #Compute the normalisation factor for the magnetic field
                normB = np.sqrt(bx * bx + by * by + bz * bz)
                func = 1./ normB * ((2. * by * bz * np.log(x + r)) + (2. * bz * bx * np.log(y + r)) + (2. * by * bx * np.log(z + r))
                    + (bz * bz - by * by) * np.arctan((x * z) / (y * r)) + (bz * bz - bx * bx) * np.arctan((y * z) / (x * r))) 
                res = -func
                return res
:::

::: {.section}
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Functions](#header-functions) {#functions}

    -   `A_drill`
    -   `A_sens`
    -   `grav_func`
    -   `magn_func`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
::: {role="main"}
::: {#section-intro .section}
Script for generating simulated cubes and sensors

Set of multiple distinct synthetic geophysical models can be created,
including two-dipping body (\"cylinders\") and layered models. For each
model a 3D voxel cube with geological structures is generated given by
their density and magnetic susceptibility as well as the corresponding
2D gravity and magnetic remote sensor measurements.

Other custom models can be included by adding a new model in function
create\_syncube()

Author: Sebastian Haan

Expand source code

    """ 
    Script for generating simulated cubes and sensors

    Set of multiple distinct synthetic geophysical models can be created, 
    including two-dipping body ("cylinders") and  layered models. 
    For each model  a 3D voxel cube with geological structures is generated 
    given by their density and magnetic susceptibility 
    as well as the corresponding 2D gravity and magnetic remote sensor measurements. 

    Other custom models can be included by adding a new model in function create_syncube()

    Author: Sebastian Haan
    """

    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import pandas as pd
    import random
    import rasterio
    from config_loader import *
    import inversion
    from sensormodel import * 
    import cubeshow as cs


    def create_syncube(modelname, voxelpos):
            """Creates synthetic cube for density and magnetic susceptibility
            
            Generates two output files, one vtk cube and one csv file

            PARAMETER

            param modelname: String, options: "layers_2", "layers_3", "cylinders" 
            param voxelpos: voxel positions (x, y, z)

            RETURN

            density cube
            susceptibility cube 
            """
            print("Creating simulated cube data ...")
            xxx, yyy, zzz = voxelpos
            x3 = xxx.reshape(yNcube, xNcube, zNcube)
            y3 = yyy.reshape(yNcube, xNcube, zNcube)
            z3 = zzz.reshape(yNcube, xNcube, zNcube)
            if modelname == 'layers_2':
                    zshift = zLcube/8. * 1. / (1 + np.exp(2.*(-y3 + zLcube/2)))
                    layer1 = 4. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.3 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.325 + zshift))))
                    cut1 = np.percentile(layer1,90)
                    layer1[layer1 < cut1] = 0.
                    layer1[layer1 >= cut1] = layer1.max()
                    layer2 = 8. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.25 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.275 + zshift))))
                    cut2 = np.percentile(layer2,90)
                    layer2[layer2 < cut2] = 0.
                    layer2[layer2 >= cut2] = layer2.max()
                    density = 0.5 + layer1 + layer2 #
                    # assume simple correlation between magnetism and density:
                    magsus = gp_coeff[1] * density
            if modelname == 'layers_3':
                    zshift = zLcube/8. * 1. / (1 + np.exp(2.*(-y3 + yLcube/2.)))
                    layer3 = 6. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.35 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.375 + zshift))))
                    cut3 = np.percentile(layer3,90)
                    layer3[layer3 < cut3] = 0.
                    layer3[layer3 >= cut3] = layer3.max()
                    layer1 = 4. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.3 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.325 + zshift))))
                    cut1 = np.percentile(layer1,90)
                    layer1[layer1 < cut1] = 0.
                    layer1[layer1 >= cut1] = layer1.max()
                    layer2 = 8. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.25 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.275 + zshift))))
                    cut2 = np.percentile(layer2,90)
                    layer2[layer2 < cut2] = 0.
                    layer2[layer2 >= cut2] = layer2.max()
                    density = 0.5 + layer1 + layer2 + layer3
                    magsus = gp_coeff[1] * density
            if modelname == 'cylinders':
                    rad = yLcube/18.
                    rc1 = ((y3-yLcube/1.3 - rad)**2 ) + ((z3 + zLcube/(4)  - rad)**2)
                    rc2 = ((y3-yLcube/4. - rad)**2 ) + ((z3 + zLcube/(4)  - rad)**2)
                    density = x3 * 0. + 0.1
                    density[rc2 <= rad**2] = 1.
                    density[rc1 <= rad**2] = 1.
                    density[(x3<xLcube/5.) | (x3>xLcube * 4./5.)] = 0.1
                    #magnetic = density
                    magsus = gp_coeff[1] * density

            # Create simulated VTK cube
            origin = (voxelpos[0].min(), voxelpos[1].min(), voxelpos[2].min())
            voxelsize = (xvoxsize, yvoxsize,zvoxsize)
            cs.create_vtkcube(density, origin, voxelsize, fname = inpath + 'simcube_' + modelname + '.vtk')

            # Save simulated data as csv file:
            newdf_head = ['x','y','z', 'DENSITY', 'MAGSUS']
            data = np.asarray([x3.flatten(), y3.flatten(), z3.flatten(), density.flatten(), magsus.flatten()])
            df= pd.DataFrame(data.T, columns = newdf_head)
            df.to_csv(inpath + 'simcube_' + modelname + '.csv', index = False)

            # Create simulated drill data from 4 drillcores:
            #select four random x,y coordinates
            dfdrill = df.copy()
            xdrill = [random.randint(2,xNcube -2) for p in range(2)] 
            ydrill = [random.randrange(2,yNcube -2) for p in range(2)]
            xdrill = np.asarray(xdrill) * xvoxsize + 0.5 * xvoxsize
            ydrill = np.asarray(ydrill) * yvoxsize+ 0.5 * yvoxsize
            dfdrill = dfdrill.loc[(dfdrill['x'].isin(xdrill)) & (dfdrill['y'].isin(ydrill))]
            dfdrill['SiteID'] = 'SiteID_' + dfdrill['x'].astype(str) + dfdrill['y'].astype(str)
            dfdrill.to_csv(inpath + 'simdrill_' + modelname + '.csv', index = False)

            return density, magsus


    def create_synsurvey(modelname, density, magsus):
            """
            Generates synthetic gravity and magntics survey data. Sensors are positiones on top of cube data

            PARAMETER

            param density: density cube data
            param magsus: magnetic susceptibiity cube data

            RETURN

            gravity 2D array
            magnetic 2D array
            """
            print("Creating simulated sensor data...")
            # Create voxel edges
            xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
            yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
            zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
            xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
            Edges = np.asarray([xEdges, yEdges, -zEdges])
            # Create sensor positions
            xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
            ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
            xx ,yy = np.meshgrid(xnew,ynew)
            zz = xx * 0. + zoff
            sensor_locations = np.asarray([xx.flatten(), yy.flatten(), zz.flatten()]).T
            # Calculate sensitivities
            gravsens, _ = A_sens(magneticField * 0., sensor_locations, Edges, 'grav')
            magsens, _ = A_sens(magneticField, sensor_locations, Edges, 'magn')
            gravfield = np.dot(gravsens, density.flatten()) # shape Nsensor, equivalent to (gravsens * properties).sum(axis = 1)
            magfield =  np.dot(magsens, magsus.flatten())
            grav2D = gravfield.reshape(xx.shape)
            magn2D = magfield.reshape(xx.shape)
            # Write csv file
            newdf_head = ['X','Y', 'GRAVITY', 'MAGNETIC']
            data = np.asarray([xx.flatten(), yy.flatten(), gravfield, magfield])
            df= pd.DataFrame(data.T, columns = newdf_head)
            df.to_csv(inpath + 'simsurveydata_' + modelname + '.csv', index = False)

            return grav2D, magn2D


    def create_simdata(modelname = "cylinders", plot = True):
            """ 
            Generates two simulated 3D cubes (density and magnetic susceptibility) 
            plus their corresponding gravity and magnetics 2D sensor data above surface.

            Sensor data is calculated via forward models as specified in sensormodel.py.
            Other settings such as cube geometry, earth's magnetic field, and sensor height 
            above ground are specified in settings. Sensor x,y positions are by default centered 
            as grid on top of voxel centers.

            PARAMETER
            param modelname: String ["layers_2", "layers_3", "cylinders"], defaults to "cylinders" 
            param plot: boolean, if True, 2D plots of sensor and simulated data are created

            RETURN
            The simulated data is saved in output directory (settings) as csv files (one for cube and one for sensor data) 
            The Two cubes are aslo saved in addition as VTK format files.
            Gravity and magentic survey data as tif file.
            Plots of sensor and vertically integrated cube data are saved in output directory.
            """
            inv = inversion.Inversion()

            voxelpos = inv.create_cubegeometry()

            # Create cubes
            density, magsus = create_syncube(modelname, voxelpos)

            # Create sensor data
            grav2D, magn2D = create_synsurvey(modelname, density, magsus)

            # Create tif survey files:
            with rasterio.open(inpath + 'gravity_simdata_' + modelname +'.tif', 'w',driver = 'GTiff', width = grav2D.shape[1], height=grav2D.shape[0], count=1, dtype='float32') as dst: 
                    dst.write(grav2D.astype(rasterio.float32), 1) 
            with rasterio.open(inpath + 'magnetic_simdata_' + modelname +'.tif', 'w',driver = 'GTiff', width = magn2D.shape[1], height=magn2D.shape[0], count=1, dtype='float32') as dst: 
                    dst.write(magn2D.astype(rasterio.float32), 1) 

            # Plot sensor data and integrated cube data along vertical axis
            if plot:
                    print("Creating plots of simulated data ...")
                    extent = np.asarray([0, xLcube, 0, yLcube])
                    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(10, 8))
                    axs[0,0].imshow(grav2D, extent = extent)
                    axs[0,0].set_title('Gravity Measurements')
                    axs[0,0].grid(True)
                    axs[0,1].imshow(magn2D, extent = extent)
                    axs[0,1].set_title('Magnetic Meauerements')
                    axs[0,1].grid(True)
                    axs[1,0].imshow(np.sum(density, axis =2), extent = extent)
                    axs[1,0].set_title('Vertical Sum Density')
                    axs[1,0].grid(True)
                    axs[1,1].imshow(np.sum(magsus, axis =2), extent = extent)
                    axs[1,1].set_title('Vertical Sum Magnetic Susceptibility')
                    axs[1,1].grid(True)
                    plt.tight_layout()
                    plt.savefig(inpath + 'figure_simdata_' + modelname +'.png', dpi = 300)
                    plt.close()
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Functions {#header-functions .section-title}
---------

` def create_simdata(modelname='cylinders', plot=True)`{.name .flex}

:   ::: {.section .desc}
    Generates two simulated 3D cubes (density and magnetic
    susceptibility) plus their corresponding gravity and magnetics 2D
    sensor data above surface.

    Sensor data is calculated via forward models as specified in
    sensormodel.py. Other settings such as cube geometry, earth\'s
    magnetic field, and sensor height above ground are specified in
    settings. Sensor x,y positions are by default centered as grid on
    top of voxel centers.

    PARAMETER param modelname: String \[\"layers\_2\", \"layers\_3\",
    \"cylinders\"\], defaults to \"cylinders\" param plot: boolean, if
    True, 2D plots of sensor and simulated data are created

    RETURN The simulated data is saved in output directory (settings) as
    csv files (one for cube and one for sensor data) The Two cubes are
    aslo saved in addition as VTK format files. Gravity and magentic
    survey data as tif file. Plots of sensor and vertically integrated
    cube data are saved in output directory.
    :::

    Expand source code

        def create_simdata(modelname = "cylinders", plot = True):
                """ 
                Generates two simulated 3D cubes (density and magnetic susceptibility) 
                plus their corresponding gravity and magnetics 2D sensor data above surface.

                Sensor data is calculated via forward models as specified in sensormodel.py.
                Other settings such as cube geometry, earth's magnetic field, and sensor height 
                above ground are specified in settings. Sensor x,y positions are by default centered 
                as grid on top of voxel centers.

                PARAMETER
                param modelname: String ["layers_2", "layers_3", "cylinders"], defaults to "cylinders" 
                param plot: boolean, if True, 2D plots of sensor and simulated data are created

                RETURN
                The simulated data is saved in output directory (settings) as csv files (one for cube and one for sensor data) 
                The Two cubes are aslo saved in addition as VTK format files.
                Gravity and magentic survey data as tif file.
                Plots of sensor and vertically integrated cube data are saved in output directory.
                """
                inv = inversion.Inversion()

                voxelpos = inv.create_cubegeometry()

                # Create cubes
                density, magsus = create_syncube(modelname, voxelpos)

                # Create sensor data
                grav2D, magn2D = create_synsurvey(modelname, density, magsus)

                # Create tif survey files:
                with rasterio.open(inpath + 'gravity_simdata_' + modelname +'.tif', 'w',driver = 'GTiff', width = grav2D.shape[1], height=grav2D.shape[0], count=1, dtype='float32') as dst: 
                        dst.write(grav2D.astype(rasterio.float32), 1) 
                with rasterio.open(inpath + 'magnetic_simdata_' + modelname +'.tif', 'w',driver = 'GTiff', width = magn2D.shape[1], height=magn2D.shape[0], count=1, dtype='float32') as dst: 
                        dst.write(magn2D.astype(rasterio.float32), 1) 

                # Plot sensor data and integrated cube data along vertical axis
                if plot:
                        print("Creating plots of simulated data ...")
                        extent = np.asarray([0, xLcube, 0, yLcube])
                        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(10, 8))
                        axs[0,0].imshow(grav2D, extent = extent)
                        axs[0,0].set_title('Gravity Measurements')
                        axs[0,0].grid(True)
                        axs[0,1].imshow(magn2D, extent = extent)
                        axs[0,1].set_title('Magnetic Meauerements')
                        axs[0,1].grid(True)
                        axs[1,0].imshow(np.sum(density, axis =2), extent = extent)
                        axs[1,0].set_title('Vertical Sum Density')
                        axs[1,0].grid(True)
                        axs[1,1].imshow(np.sum(magsus, axis =2), extent = extent)
                        axs[1,1].set_title('Vertical Sum Magnetic Susceptibility')
                        axs[1,1].grid(True)
                        plt.tight_layout()
                        plt.savefig(inpath + 'figure_simdata_' + modelname +'.png', dpi = 300)
                        plt.close()

` def create_syncube(modelname, voxelpos)`{.name .flex}

:   ::: {.section .desc}
    Creates synthetic cube for density and magnetic susceptibility

    Generates two output files, one vtk cube and one csv file

    PARAMETER

    param modelname: String, options: \"layers\_2\", \"layers\_3\",
    \"cylinders\" param voxelpos: voxel positions (x, y, z)

    RETURN

    density cube susceptibility cube
    :::

    Expand source code

        def create_syncube(modelname, voxelpos):
                """Creates synthetic cube for density and magnetic susceptibility
                
                Generates two output files, one vtk cube and one csv file

                PARAMETER

                param modelname: String, options: "layers_2", "layers_3", "cylinders" 
                param voxelpos: voxel positions (x, y, z)

                RETURN

                density cube
                susceptibility cube 
                """
                print("Creating simulated cube data ...")
                xxx, yyy, zzz = voxelpos
                x3 = xxx.reshape(yNcube, xNcube, zNcube)
                y3 = yyy.reshape(yNcube, xNcube, zNcube)
                z3 = zzz.reshape(yNcube, xNcube, zNcube)
                if modelname == 'layers_2':
                        zshift = zLcube/8. * 1. / (1 + np.exp(2.*(-y3 + zLcube/2)))
                        layer1 = 4. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.3 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.325 + zshift))))
                        cut1 = np.percentile(layer1,90)
                        layer1[layer1 < cut1] = 0.
                        layer1[layer1 >= cut1] = layer1.max()
                        layer2 = 8. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.25 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.275 + zshift))))
                        cut2 = np.percentile(layer2,90)
                        layer2[layer2 < cut2] = 0.
                        layer2[layer2 >= cut2] = layer2.max()
                        density = 0.5 + layer1 + layer2 #
                        # assume simple correlation between magnetism and density:
                        magsus = gp_coeff[1] * density
                if modelname == 'layers_3':
                        zshift = zLcube/8. * 1. / (1 + np.exp(2.*(-y3 + yLcube/2.)))
                        layer3 = 6. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.35 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.375 + zshift))))
                        cut3 = np.percentile(layer3,90)
                        layer3[layer3 < cut3] = 0.
                        layer3[layer3 >= cut3] = layer3.max()
                        layer1 = 4. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.3 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.325 + zshift))))
                        cut1 = np.percentile(layer1,90)
                        layer1[layer1 < cut1] = 0.
                        layer1[layer1 >= cut1] = layer1.max()
                        layer2 = 8. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.25 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.275 + zshift))))
                        cut2 = np.percentile(layer2,90)
                        layer2[layer2 < cut2] = 0.
                        layer2[layer2 >= cut2] = layer2.max()
                        density = 0.5 + layer1 + layer2 + layer3
                        magsus = gp_coeff[1] * density
                if modelname == 'cylinders':
                        rad = yLcube/18.
                        rc1 = ((y3-yLcube/1.3 - rad)**2 ) + ((z3 + zLcube/(4)  - rad)**2)
                        rc2 = ((y3-yLcube/4. - rad)**2 ) + ((z3 + zLcube/(4)  - rad)**2)
                        density = x3 * 0. + 0.1
                        density[rc2 <= rad**2] = 1.
                        density[rc1 <= rad**2] = 1.
                        density[(x3<xLcube/5.) | (x3>xLcube * 4./5.)] = 0.1
                        #magnetic = density
                        magsus = gp_coeff[1] * density

                # Create simulated VTK cube
                origin = (voxelpos[0].min(), voxelpos[1].min(), voxelpos[2].min())
                voxelsize = (xvoxsize, yvoxsize,zvoxsize)
                cs.create_vtkcube(density, origin, voxelsize, fname = inpath + 'simcube_' + modelname + '.vtk')

                # Save simulated data as csv file:
                newdf_head = ['x','y','z', 'DENSITY', 'MAGSUS']
                data = np.asarray([x3.flatten(), y3.flatten(), z3.flatten(), density.flatten(), magsus.flatten()])
                df= pd.DataFrame(data.T, columns = newdf_head)
                df.to_csv(inpath + 'simcube_' + modelname + '.csv', index = False)

                # Create simulated drill data from 4 drillcores:
                #select four random x,y coordinates
                dfdrill = df.copy()
                xdrill = [random.randint(2,xNcube -2) for p in range(2)] 
                ydrill = [random.randrange(2,yNcube -2) for p in range(2)]
                xdrill = np.asarray(xdrill) * xvoxsize + 0.5 * xvoxsize
                ydrill = np.asarray(ydrill) * yvoxsize+ 0.5 * yvoxsize
                dfdrill = dfdrill.loc[(dfdrill['x'].isin(xdrill)) & (dfdrill['y'].isin(ydrill))]
                dfdrill['SiteID'] = 'SiteID_' + dfdrill['x'].astype(str) + dfdrill['y'].astype(str)
                dfdrill.to_csv(inpath + 'simdrill_' + modelname + '.csv', index = False)

                return density, magsus

` def create_synsurvey(modelname, density, magsus)`{.name .flex}

:   ::: {.section .desc}
    Generates synthetic gravity and magntics survey data. Sensors are
    positiones on top of cube data

    PARAMETER

    param density: density cube data param magsus: magnetic
    susceptibiity cube data

    RETURN

    gravity 2D array magnetic 2D array
    :::

    Expand source code

        def create_synsurvey(modelname, density, magsus):
                """
                Generates synthetic gravity and magntics survey data. Sensors are positiones on top of cube data

                PARAMETER

                param density: density cube data
                param magsus: magnetic susceptibiity cube data

                RETURN

                gravity 2D array
                magnetic 2D array
                """
                print("Creating simulated sensor data...")
                # Create voxel edges
                xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
                yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
                zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
                xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
                Edges = np.asarray([xEdges, yEdges, -zEdges])
                # Create sensor positions
                xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
                ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
                xx ,yy = np.meshgrid(xnew,ynew)
                zz = xx * 0. + zoff
                sensor_locations = np.asarray([xx.flatten(), yy.flatten(), zz.flatten()]).T
                # Calculate sensitivities
                gravsens, _ = A_sens(magneticField * 0., sensor_locations, Edges, 'grav')
                magsens, _ = A_sens(magneticField, sensor_locations, Edges, 'magn')
                gravfield = np.dot(gravsens, density.flatten()) # shape Nsensor, equivalent to (gravsens * properties).sum(axis = 1)
                magfield =  np.dot(magsens, magsus.flatten())
                grav2D = gravfield.reshape(xx.shape)
                magn2D = magfield.reshape(xx.shape)
                # Write csv file
                newdf_head = ['X','Y', 'GRAVITY', 'MAGNETIC']
                data = np.asarray([xx.flatten(), yy.flatten(), gravfield, magfield])
                df= pd.DataFrame(data.T, columns = newdf_head)
                df.to_csv(inpath + 'simsurveydata_' + modelname + '.csv', index = False)

                return grav2D, magn2D
:::

::: {.section}
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Functions](#header-functions) {#functions}

    -   `create_simdata`
    -   `create_syncube`
    -   `create_synsurvey`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
::: {role="main"}
::: {#section-intro .section}
Kernel library for Gaussian Processes including sparse kernels and
cross-covariance terms

The choice for an appropriate covariance function is important, as the
GP\'s output directly depends on it. These parameters of the covariance
function are referred to as the hyperparameters of the GP, which can be
either given by a fixed covariance scale and noise, or learned from data
by optimising the marginal likelihood. To handle the computational
problem of inverting a large covariance matrix, sparse covariance
function are included here as well. To take fully into account
cross-covariances between multiple model parameters (e.g., rock density
and magnetic susceptibility), we construct cross-covariance terms (here
labeled with edning \"2\") between all kernel pairs. One important
requirement for constructing cross-covariance terms is that they must be
defined to be both positive semi-definite and informative; for an
overview how to construct such as matrix in detail see Melkumyan 2011.

Author: Sebastian Haan

Expand source code

    """
    Kernel library for Gaussian Processes including sparse kernels and cross-covariance terms

    The choice for an appropriate covariance function is important, 
    as the GP's output directly depends on it. These parameters of the covariance function are 
    referred to as the hyperparameters of the GP, which can be either given by a fixed covariance scale 
    and noise, or learned from data by optimising the marginal likelihood. To handle the computational problem 
    of inverting a large covariance matrix, sparse covariance function are included here as well.
    To take fully into account cross-covariances between multiple model parameters 
    (e.g., rock density and magnetic susceptibility), we construct cross-covariance terms (here labeled with edning "2")
    between all kernel pairs. One important requirement for constructing cross-covariance terms is that they must be defined 
    to be both positive semi-definite and informative; for an overview how to construct such as matrix in detail see Melkumyan 2011.

    Author: Sebastian Haan
    """
    from scipy import reshape, sqrt, identity
    import numpy as np


    def calcGridPoints3D(Lpix, pixscale):
        """
        returns grid points for distance matrix calculation.
        :param Lpix: number of pixels in each dimension as array (xLpix, yLpix, zLpix)
        :param pixscale: pixelscale in each dimension as array (xpixscale, ypixscale, zpixscale)
        """
        Lpix = np.asarray(Lpix)
        pixscale = np.asarray(pixscale)
        xLpix, yLpix, zLpix = Lpix[0], Lpix[1], Lpix[2]
        xpixscale, ypixscale, zpixscale = pixscale[0], pixscale[1], pixscale[2]
        xrange = np.arange(1, xLpix+1) * xpixscale
        yrange = np.arange(1, yLpix+1) * ypixscale
        zrange = np.arange(1, zLpix+1) * zpixscale
        _xg, _yg, _zg = np.meshgrid(xrange, yrange, zrange)
        xr, yr, zr = _xg.ravel(), _yg.ravel(), _zg.ravel()
        return np.asarray([xr, yr, zr]).T


    def calcDistanceMatrix(nDimPoints, 
                           distFunc=lambda deltaPoint: np.sum(deltaPoint[d]**2 for d in range(len(deltaPoint)))):
        """ Returns the matrix of squared distances from one coordinate to any other
        :param nDimPoints: list of n-dim tuples
        :param distFunc: calculates the distance based on the differences
        """
        nDimPoints = np.array(nDimPoints)
        dim = len(nDimPoints[0])
        delta = [None]*dim
        for d in range(dim):
            data = nDimPoints[:,d]
            delta[d] = data - np.reshape(data,(len(data),1)) # computes all possible combinations

        dist = distFunc(delta)
        #dist = dist + np.identity(len(data))*dist.max() # eliminate self matching
        # returns  squared distance:
        return dist 


    def calc_square_distances2D(Lpix, pixscale):
        """
        Initialize (squared) distance matrix for stationary kernel.
        """
        Lpix = np.asarray(Lpix)
        pixscale = np.asarray(pixscale)
        xLpix, yLpix = Lpix[0], Lpix[1]
        xpixscale, ypixscale = pixscale[0], pixscale[1]
        xrange = (np.arange(0, xLpix) - xLpix/2.0) * xpixscale
        yrange = (np.arange(0, yLpix) - xLpix/2.0) * ypixscale
        _xg, _yg = np.meshgrid(xrange, yrange)
        xr, yr = _xg.ravel(), _yg.ravel()
        Dx = xr[:, np.newaxis] - xr[np.newaxis,:]
        Dy = yr[:, np.newaxis] - yr[np.newaxis,:]
        return Dx**2 + Dy**2


    def gpkernel(D2, gamma):
        """2-D round  RBF kernel, with length scale = standard deviation of
        the PSF of a Gaussian process scene drawn from this kernel.
        Squared Exponential kernel
        :param D2: pairwise square distances
        :param gamma: kernel length scale
        """
        return np.exp(-0.5 * D2/gamma**2)

    def gpkernel2(D2, gammas):
        """exp squared x qxp squared kernel
        the PSF of a Gaussian process scene drawn from this kernel.
        Squared Exponential kernel
        :param D2: pairwise square distances
        :param gamma: kernel length scales (gamma1, gamma2)
        """
        l1 = gammas[0]
        l2 =gammas[1]
        return np.sqrt(2.*l1 * l2/(l1**2 + l2**2)) * np.exp(- D2/(l1**2 + l2**2))

    def gpkernel_sparse(D2, gamma):
        """2-D round sparse RBF kernel, defined in Melkumyan and Ramos, 2009
        lengthscale is roughly equivlanet to 4 times the lengthcale of squared exponential
        Same as in gpcubesolve_r1.py
        :param D2: pairwise square distances
        :param gamma: kernel length scale
        """
        D2 = np.sqrt(D2)
        #gamma = 4 * gamma
        res = np.zeros_like(D2)
        res[D2 < gamma] = (2 + np.cos(2*np.pi * D2[D2 < gamma] /gamma))/3.*(1-D2[D2 < gamma] /gamma) + 1/(2.*np.pi) * np.sin(2*np.pi*D2[D2 < gamma] /gamma)
        # Remove floating errors
        res[res<0.] = 0.
        return res

    def gpkernel_sparse2(D2, gammas):
        """ sparse matrix x sparse matrix
        :param D2: pairwise square distances
        :param gamma: kernel length scales (gamma1, gamma2)
        """
        D2 = np.sqrt(D2)
        l1 = gammas[0]
        l2 = gammas[1]
        # offset for avoiding non-zero terms
        if l1 == l2:
            l2 += 1e-3 * l2
        lmean = np.mean([l1, l2])
        lmin = np.min([l1,l2])
        lmax = np.max([l1,l2])
        #if D2 >= lmean:
        res = np.zeros_like(D2)
        #if D2 <= abs(l2 - l1) / 2.:
        res[D2 <= abs(l2 - l1) / 2.] = 2./(3*np.sqrt(l1*l2)) * (lmin + 1/np.pi*lmax**3/(lmax**2 - lmin**2) * np.sin(np.pi*lmin/lmax*np.cos(2*np.pi*D2[D2 <= abs(l2 - l1) / 2.]/lmax)))
        #if (D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.):
        res[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)] = 2./(3*np.sqrt(l1*l2)) * (lmean - D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)] + l1**3*np.sin(np.pi*(l2-2.*D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)])/l1)/(2*np.pi*(l1**2-l2**2)) - l2**3*np.sin(np.pi*(l1-2.*D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)])/l2)/(2*np.pi*(l1**2-l2**2)))
        # Remove floating errors
        res[res<0.] = 0.
        return res

    def gpkernel_matern32(D2, gamma):
        ''' Matern3/2 kernel
        :param D2: pairwise square distances
        :param gamma: kernel length scale
        '''
        nu = np.sqrt(3) * np.sqrt(D2) / gamma    
        return (1 + nu) * np.exp(-nu)

    def gpkernel_matern32_2(D2, gammas):
        ''' Matern3/2 x Matern3/2 kernel
        :param D2: pairwise square distances
        :param gammas: kernel length scales (gamma1, gamma2)
        '''
        l1 = gammas[0]
        l2 = gammas[1]
        norm = 2*np.sqrt(l1 * l2)/(l1**2 - l2**2)
        return  norm * (l1 * np.exp(- np.sqrt(3*D2)/l1) - l2 * np.exp(- np.sqrt(3*D2)/l2))

    def create_cov(D2, gplength, crossweights = [1,1,1], fkernel = 'sparse'):
        """
        Compute cross-correlation covariance matrix
        :param D2: distance matrix
        :param gplengths: scalar or array with lengthscales for each of the properties (up to three)
        :param crossweights: array with cross-correlation coefficents correspdoning to the output property (w1,w2,w3)
        :param fkernel: kernel function, which can be either 'sparse' (Default), 'exp' (squared exponential) or 'matern32'

        Cross weights (w1,w2,w3) defined as 
        w1: density - drill correltion
        w2: magentic - drill correlation 
        w3: density - magnetic correlation

        Return: cross-correlation matrix
        """
        # first calculate kernel
        params= np.asarray(gplength)
        if params[1] == params[0]:
            params[1] = 1.01 * params[0]
        if params[2] == params[0]:
            params[1] = 1.02 * params[0]
        if params[2] == params[1]:
            params[2] = 1.01 * params[1]
        w1, w2, w3 = np.asarray(crossweights)
        #kcov = np.asarray([gpkernel(D2, params[0]), gpkernel2(D2, params[0:2]), gpkernel2(D2, params[0:2]), gpkernel(D2, params[1])]).T
        if fkernel == 'matern32':
            kcov1 = np.vstack([gpkernel_matern32(D2, params[0]), w3 * gpkernel_matern32_2(D2, params[[0,1]]), w1 * gpkernel_matern32_2(D2, params[[0,2]])])
            kcov2 = np.vstack([w3 * gpkernel_matern32_2(D2, params[[1,0]]), gpkernel_matern32(D2, params[1]), w2 * gpkernel_matern32_2(D2, params[[1,2]])])
            kcov3 = np.vstack([w1* gpkernel_matern32_2(D2, params[[2,0]]), w2 * gpkernel_matern32_2(D2, params[[2,1]]), gpkernel_matern32(D2, params[2])])
        if fkernel == 'sparse':
            kcov1 = np.vstack([gpkernel_sparse(D2, params[0]), w3 * gpkernel_sparse2(D2, params[[0,1]]), w1 * gpkernel_sparse2(D2, params[[0,2]])])
            kcov2 = np.vstack([w3 * gpkernel_sparse2(D2, params[[1,0]]), gpkernel_sparse(D2, params[1]), w2 * gpkernel_sparse2(D2, params[[1,2]])])
            kcov3 = np.vstack([w1 * gpkernel_sparse2(D2, params[[2,0]]), w2 * gpkernel_sparse2(D2, params[[2,1]]), gpkernel_sparse(D2, params[2])])
        if fkernel == 'exp':
            kcov1 = np.vstack([gpkernel(D2, params[0]), w3 * gpkernel2(D2, params[[0,1]]), w1 * gpkernel2(D2, params[[0,2]])])
            kcov2 = np.vstack([w3 * gpkernel2(D2, params[[1,0]]), gpkernel(D2, params[1]), w2 * gpkernel2(D2, params[[1,2]])])
            kcov3 = np.vstack([w1 * gpkernel2(D2, params[[2,0]]), w2 * gpkernel2(D2, params[[2,1]]), gpkernel(D2, params[2])])
        return np.hstack([kcov1, kcov2, kcov3])
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Functions {#header-functions .section-title}
---------

` def calcDistanceMatrix(nDimPoints, distFunc=<function <lambda>>)`{.name .flex}

:   ::: {.section .desc}
    Returns the matrix of squared distances from one coordinate to any
    other :param nDimPoints: list of n-dim tuples :param distFunc:
    calculates the distance based on the differences
    :::

    Expand source code

        def calcDistanceMatrix(nDimPoints, 
                               distFunc=lambda deltaPoint: np.sum(deltaPoint[d]**2 for d in range(len(deltaPoint)))):
            """ Returns the matrix of squared distances from one coordinate to any other
            :param nDimPoints: list of n-dim tuples
            :param distFunc: calculates the distance based on the differences
            """
            nDimPoints = np.array(nDimPoints)
            dim = len(nDimPoints[0])
            delta = [None]*dim
            for d in range(dim):
                data = nDimPoints[:,d]
                delta[d] = data - np.reshape(data,(len(data),1)) # computes all possible combinations

            dist = distFunc(delta)
            #dist = dist + np.identity(len(data))*dist.max() # eliminate self matching
            # returns  squared distance:
            return dist 

` def calcGridPoints3D(Lpix, pixscale)`{.name .flex}

:   ::: {.section .desc}
    returns grid points for distance matrix calculation. :param Lpix:
    number of pixels in each dimension as array (xLpix, yLpix, zLpix)
    :param pixscale: pixelscale in each dimension as array (xpixscale,
    ypixscale, zpixscale)
    :::

    Expand source code

        def calcGridPoints3D(Lpix, pixscale):
            """
            returns grid points for distance matrix calculation.
            :param Lpix: number of pixels in each dimension as array (xLpix, yLpix, zLpix)
            :param pixscale: pixelscale in each dimension as array (xpixscale, ypixscale, zpixscale)
            """
            Lpix = np.asarray(Lpix)
            pixscale = np.asarray(pixscale)
            xLpix, yLpix, zLpix = Lpix[0], Lpix[1], Lpix[2]
            xpixscale, ypixscale, zpixscale = pixscale[0], pixscale[1], pixscale[2]
            xrange = np.arange(1, xLpix+1) * xpixscale
            yrange = np.arange(1, yLpix+1) * ypixscale
            zrange = np.arange(1, zLpix+1) * zpixscale
            _xg, _yg, _zg = np.meshgrid(xrange, yrange, zrange)
            xr, yr, zr = _xg.ravel(), _yg.ravel(), _zg.ravel()
            return np.asarray([xr, yr, zr]).T

` def calc_square_distances2D(Lpix, pixscale)`{.name .flex}

:   ::: {.section .desc}
    Initialize (squared) distance matrix for stationary kernel.
    :::

    Expand source code

        def calc_square_distances2D(Lpix, pixscale):
            """
            Initialize (squared) distance matrix for stationary kernel.
            """
            Lpix = np.asarray(Lpix)
            pixscale = np.asarray(pixscale)
            xLpix, yLpix = Lpix[0], Lpix[1]
            xpixscale, ypixscale = pixscale[0], pixscale[1]
            xrange = (np.arange(0, xLpix) - xLpix/2.0) * xpixscale
            yrange = (np.arange(0, yLpix) - xLpix/2.0) * ypixscale
            _xg, _yg = np.meshgrid(xrange, yrange)
            xr, yr = _xg.ravel(), _yg.ravel()
            Dx = xr[:, np.newaxis] - xr[np.newaxis,:]
            Dy = yr[:, np.newaxis] - yr[np.newaxis,:]
            return Dx**2 + Dy**2

` def create_cov(D2, gplength, crossweights=[1, 1, 1], fkernel='sparse')`{.name .flex}

:   ::: {.section .desc}
    Compute cross-correlation covariance matrix :param D2: distance
    matrix :param gplengths: scalar or array with lengthscales for each
    of the properties (up to three) :param crossweights: array with
    cross-correlation coefficents correspdoning to the output property
    (w1,w2,w3) :param fkernel: kernel function, which can be either
    \'sparse\' (Default), \'exp\' (squared exponential) or \'matern32\'

    Cross weights (w1,w2,w3) defined as w1: density - drill correltion
    w2: magentic - drill correlation w3: density - magnetic correlation

    Return: cross-correlation matrix
    :::

    Expand source code

        def create_cov(D2, gplength, crossweights = [1,1,1], fkernel = 'sparse'):
            """
            Compute cross-correlation covariance matrix
            :param D2: distance matrix
            :param gplengths: scalar or array with lengthscales for each of the properties (up to three)
            :param crossweights: array with cross-correlation coefficents correspdoning to the output property (w1,w2,w3)
            :param fkernel: kernel function, which can be either 'sparse' (Default), 'exp' (squared exponential) or 'matern32'

            Cross weights (w1,w2,w3) defined as 
            w1: density - drill correltion
            w2: magentic - drill correlation 
            w3: density - magnetic correlation

            Return: cross-correlation matrix
            """
            # first calculate kernel
            params= np.asarray(gplength)
            if params[1] == params[0]:
                params[1] = 1.01 * params[0]
            if params[2] == params[0]:
                params[1] = 1.02 * params[0]
            if params[2] == params[1]:
                params[2] = 1.01 * params[1]
            w1, w2, w3 = np.asarray(crossweights)
            #kcov = np.asarray([gpkernel(D2, params[0]), gpkernel2(D2, params[0:2]), gpkernel2(D2, params[0:2]), gpkernel(D2, params[1])]).T
            if fkernel == 'matern32':
                kcov1 = np.vstack([gpkernel_matern32(D2, params[0]), w3 * gpkernel_matern32_2(D2, params[[0,1]]), w1 * gpkernel_matern32_2(D2, params[[0,2]])])
                kcov2 = np.vstack([w3 * gpkernel_matern32_2(D2, params[[1,0]]), gpkernel_matern32(D2, params[1]), w2 * gpkernel_matern32_2(D2, params[[1,2]])])
                kcov3 = np.vstack([w1* gpkernel_matern32_2(D2, params[[2,0]]), w2 * gpkernel_matern32_2(D2, params[[2,1]]), gpkernel_matern32(D2, params[2])])
            if fkernel == 'sparse':
                kcov1 = np.vstack([gpkernel_sparse(D2, params[0]), w3 * gpkernel_sparse2(D2, params[[0,1]]), w1 * gpkernel_sparse2(D2, params[[0,2]])])
                kcov2 = np.vstack([w3 * gpkernel_sparse2(D2, params[[1,0]]), gpkernel_sparse(D2, params[1]), w2 * gpkernel_sparse2(D2, params[[1,2]])])
                kcov3 = np.vstack([w1 * gpkernel_sparse2(D2, params[[2,0]]), w2 * gpkernel_sparse2(D2, params[[2,1]]), gpkernel_sparse(D2, params[2])])
            if fkernel == 'exp':
                kcov1 = np.vstack([gpkernel(D2, params[0]), w3 * gpkernel2(D2, params[[0,1]]), w1 * gpkernel2(D2, params[[0,2]])])
                kcov2 = np.vstack([w3 * gpkernel2(D2, params[[1,0]]), gpkernel(D2, params[1]), w2 * gpkernel2(D2, params[[1,2]])])
                kcov3 = np.vstack([w1 * gpkernel2(D2, params[[2,0]]), w2 * gpkernel2(D2, params[[2,1]]), gpkernel(D2, params[2])])
            return np.hstack([kcov1, kcov2, kcov3])

` def gpkernel(D2, gamma)`{.name .flex}

:   ::: {.section .desc}
    2-D round RBF kernel, with length scale = standard deviation of the
    PSF of a Gaussian process scene drawn from this kernel. Squared
    Exponential kernel :param D2: pairwise square distances :param
    gamma: kernel length scale
    :::

    Expand source code

        def gpkernel(D2, gamma):
            """2-D round  RBF kernel, with length scale = standard deviation of
            the PSF of a Gaussian process scene drawn from this kernel.
            Squared Exponential kernel
            :param D2: pairwise square distances
            :param gamma: kernel length scale
            """
            return np.exp(-0.5 * D2/gamma**2)

` def gpkernel2(D2, gammas)`{.name .flex}

:   ::: {.section .desc}
    exp squared x qxp squared kernel the PSF of a Gaussian process scene
    drawn from this kernel. Squared Exponential kernel :param D2:
    pairwise square distances :param gamma: kernel length scales
    (gamma1, gamma2)
    :::

    Expand source code

        def gpkernel2(D2, gammas):
            """exp squared x qxp squared kernel
            the PSF of a Gaussian process scene drawn from this kernel.
            Squared Exponential kernel
            :param D2: pairwise square distances
            :param gamma: kernel length scales (gamma1, gamma2)
            """
            l1 = gammas[0]
            l2 =gammas[1]
            return np.sqrt(2.*l1 * l2/(l1**2 + l2**2)) * np.exp(- D2/(l1**2 + l2**2))

` def gpkernel_matern32(D2, gamma)`{.name .flex}

:   ::: {.section .desc}
    Matern3/2 kernel :param D2: pairwise square distances :param gamma:
    kernel length scale
    :::

    Expand source code

        def gpkernel_matern32(D2, gamma):
            ''' Matern3/2 kernel
            :param D2: pairwise square distances
            :param gamma: kernel length scale
            '''
            nu = np.sqrt(3) * np.sqrt(D2) / gamma    
            return (1 + nu) * np.exp(-nu)

` def gpkernel_matern32_2(D2, gammas)`{.name .flex}

:   ::: {.section .desc}
    Matern3/2 x Matern3/2 kernel :param D2: pairwise square distances
    :param gammas: kernel length scales (gamma1, gamma2)
    :::

    Expand source code

        def gpkernel_matern32_2(D2, gammas):
            ''' Matern3/2 x Matern3/2 kernel
            :param D2: pairwise square distances
            :param gammas: kernel length scales (gamma1, gamma2)
            '''
            l1 = gammas[0]
            l2 = gammas[1]
            norm = 2*np.sqrt(l1 * l2)/(l1**2 - l2**2)
            return  norm * (l1 * np.exp(- np.sqrt(3*D2)/l1) - l2 * np.exp(- np.sqrt(3*D2)/l2))

` def gpkernel_sparse(D2, gamma)`{.name .flex}

:   ::: {.section .desc}
    2-D round sparse RBF kernel, defined in Melkumyan and Ramos, 2009
    lengthscale is roughly equivlanet to 4 times the lengthcale of
    squared exponential Same as in gpcubesolve\_r1.py :param D2:
    pairwise square distances :param gamma: kernel length scale
    :::

    Expand source code

        def gpkernel_sparse(D2, gamma):
            """2-D round sparse RBF kernel, defined in Melkumyan and Ramos, 2009
            lengthscale is roughly equivlanet to 4 times the lengthcale of squared exponential
            Same as in gpcubesolve_r1.py
            :param D2: pairwise square distances
            :param gamma: kernel length scale
            """
            D2 = np.sqrt(D2)
            #gamma = 4 * gamma
            res = np.zeros_like(D2)
            res[D2 < gamma] = (2 + np.cos(2*np.pi * D2[D2 < gamma] /gamma))/3.*(1-D2[D2 < gamma] /gamma) + 1/(2.*np.pi) * np.sin(2*np.pi*D2[D2 < gamma] /gamma)
            # Remove floating errors
            res[res<0.] = 0.
            return res

` def gpkernel_sparse2(D2, gammas)`{.name .flex}

:   ::: {.section .desc}
    sparse matrix x sparse matrix :param D2: pairwise square distances
    :param gamma: kernel length scales (gamma1, gamma2)
    :::

    Expand source code

        def gpkernel_sparse2(D2, gammas):
            """ sparse matrix x sparse matrix
            :param D2: pairwise square distances
            :param gamma: kernel length scales (gamma1, gamma2)
            """
            D2 = np.sqrt(D2)
            l1 = gammas[0]
            l2 = gammas[1]
            # offset for avoiding non-zero terms
            if l1 == l2:
                l2 += 1e-3 * l2
            lmean = np.mean([l1, l2])
            lmin = np.min([l1,l2])
            lmax = np.max([l1,l2])
            #if D2 >= lmean:
            res = np.zeros_like(D2)
            #if D2 <= abs(l2 - l1) / 2.:
            res[D2 <= abs(l2 - l1) / 2.] = 2./(3*np.sqrt(l1*l2)) * (lmin + 1/np.pi*lmax**3/(lmax**2 - lmin**2) * np.sin(np.pi*lmin/lmax*np.cos(2*np.pi*D2[D2 <= abs(l2 - l1) / 2.]/lmax)))
            #if (D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.):
            res[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)] = 2./(3*np.sqrt(l1*l2)) * (lmean - D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)] + l1**3*np.sin(np.pi*(l2-2.*D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)])/l1)/(2*np.pi*(l1**2-l2**2)) - l2**3*np.sin(np.pi*(l1-2.*D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)])/l2)/(2*np.pi*(l1**2-l2**2)))
            # Remove floating errors
            res[res<0.] = 0.
            return res
:::

::: {.section}
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Functions](#header-functions) {#functions}

    -   `calcDistanceMatrix`
    -   `calcGridPoints3D`
    -   `calc_square_distances2D`
    -   `create_cov`
    -   `gpkernel`
    -   `gpkernel2`
    -   `gpkernel_matern32`
    -   `gpkernel_matern32_2`
    -   `gpkernel_sparse`
    -   `gpkernel_sparse2`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
