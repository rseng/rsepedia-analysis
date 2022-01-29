---
title: 'flowTorch - a Python library for analysis and reduced-order modeling of fluid flows'
tags:
  - Python
  - PyTorch
  - fluid flows
  - reduced-order modeling
  - modal analysis
authors:
  - name: Andre Weiner
    orcid: 0000-0001-5617-1560
    affiliation: 1
  - name: Richard Semaan
    orcid: 0000-0002-3219-0545
    affiliation: 1
affiliations:
 - name: Technische Universität Braunschweig, Institute of Fluid Mechanics, Flow Modeling and Control Group
   index: 1
date: 05 August 2021
bibliography: paper.bib
---

# Summary

The `flowTorch` library enables researchers to access, analyze, and model
fluid flow data from experiments or numerical simulations. Instead of a black-box end-to-end solution,
`flowTorch` provides modular components allowing to assemble transparent and reproducible workflows with ease. Popular
data formats for fluid flows like [OpenFOAM](https://www.openfoam.com/), [VTK](https://vtk.org/), or
[DaVis](https://www.lavision.de/en/products/davis-software/) may be accessed via a common interface in
a few lines of Python code. Internally, the data are organized as [PyTorch](https://pytorch.org/) tensors.
Relying on PyTorch tensors as primary data structure enables fast array operations, parallel processing on
CPU and GPU, and exploration of novel deep learning-based analysis and modeling approaches.
The `flowTorch` packages also includes a collection of Jupyter notebooks demonstrating  how to apply the
library components in a variety of different use cases, e.g., finding coherent flow structures with modal analysis or creating
reduced-order models.

# Statement of need

Thanks to the increased processing power of modern hardware, fluid flow experiments as well as numerical simulations
are producing vast amounts of highly resolved, complex data. Those data offer great opportunities to optimize industrial processes or to understand natural phenomena. As modern datasets continue to grow, post-processing pipelines will be increasingly important for synthesizing different data formats and facilitating complex data analysis. While most researchers prefer simple text-encoded comma-separated value (CSV) files, big datasets require
special binary formats, such as [HDF5](https://www.hdfgroup.org/solutions/hdf5/) or [NetCDF](https://en.wikipedia.org/wiki/NetCDF).
If the data are associated with a structured or an unstructured mesh, VTK files are a popular choice. Other simulation libraries for
fluid flows, like OpenFOAM, organize mesh and field data in custom folder and file structures. On the experimental side, software packages
like DaVis allow exporting particle image velocimetry (PIV) snapshots as CSV files. Reading CSV files can be a daunting task, too. A sequence of
snapshots might be organized in one or multiple files. If the data are stored in a single file, the file must be read first and then the individual snapshots must be extracted following some initially unknown pattern. If the data are spread out over multiple files, the time might be encoded in the file name, but it could be also the case that
the files are located in individual folders whose names encode the time. The latter structure is typical for OpenFOAM run time post-processing data.
Moreover, different software packages will create different file headers, which may have to be parsed or sometimes ignored. CSV, VTK, or OpenFOAM data may come
as binary or text-encoded files. This list is by no means comprehensive in terms of available formats and presents only the tip of the iceberg.

A common research task may be to compare and combine different data sources of the same fluid flow problem for cross-validation
or to leverage each source's strengths in different kinds of analysis. A typical example would be to compare or combine PIV data with sampled
planes extracted from a numerical simulation. The simulation offers greater details and additional field information, while the PIV experiment is more trustworthy since it
is closer to the real application. The PIV data may have to be processed and cleaned before using it in consecutive analysis steps. Often, significant research time is spent on accessing, converting, and processing the data with different tools and different formats to finally analyze the data in yet another tool.
Text-encoded file format might be convenient at first when exchanging data between tools, but for large datasets the additional conversion is unsuitable.

`flowTorch` aims to simplify access to data by providing a unified interface to various data formats via the subpackage `flowtorch.data`. Accessing data from a
distributed OpenFOAM simulation is as easy as loading VTK or PIV data and requires only a few lines of Python code. All field data
are converted internally to PyTorch tensors [@paszke2015]. Once the data are available as PyTorch tensors, further processing steps like scaling, clipping, masking, splitting, or merging are readily available as single function calls. The same is true for computing the mean, the standard deviation, histograms, or quantiles. Modal analysis techniques, like dynamic mode decomposition (DMD)[@schmid2010; @kutz2016] and proper orthogonal decomposition (POD)[@brunton2019; @semaan2020], are available via the subpackage `flowtorch.analysis`. The third subpackage, `flowtorch.rom`, enables adding reduced-order models (ROMs), like cluster-based network modeling (CNM)[@fernex2021], to the post-processing pipeline. Computationally intensive tasks may be offloaded to the GPU if needed, which greatly accelerates parameter studies. The entire analysis workflow described in the previous section can be performed in a single ecosystem sketched in \autoref{fig:ft_structure}. Moreover, re-using an analysis pipeline in a different problem setting is straightforward.

![Components of flowTorch and library dependencies.\label{fig:ft_structure}](media/flowtorch_components_plain.pdf){ width=90% }

Besides the subpackages already available in `flowTorch`, the library also integrates nicely with related software packages like [ParaView](https://www.paraview.org/) or [VisIt](https://visit-dav.github.io/visit-website/index.html) for mesh-based post-processing as well as specialized analysis and modeling packages like PyDMD [@demo2018], PySINDy [@desilva2020], or [modred](https://github.com/belson17/modred). Rather than re-implementing functionality already existing in other established libraries, `flowTorch` wraps around them to simplify their usage and streamline the overall post-processing pipeline. For example, we use ParaView's [vtk](https://pypi.org/project/vtk/) package to access various types of VTK files in Python. Gathering point coordinates, write times, or snapshots from several VTK files requires very different steps than when dealing with OpenFOAM or DaVis data.  However, due to the common interface to data sources in `flowTorch`, these tasks appear to be exactly the same for the user. In contrast to `flowTorch`, PyDMD offers a wide range of DMD variants but does not provide access to data. If an advanced DMD algorithm is required, our library can be used to access and pre-process a dataset, before PyDMD is used to perform the modal decomposition.

Another more general issue we want to address is the reproducibility of research outcomes. Popular algorithms, like POD or DMD, may be relatively easy to
implement with libraries like NumPy, SciPy, or PyTorch. However,
applying these algorithms to real datasets typically requires several pre-processing steps, like cropping, clipping, or normalizing the
data, and careful tuning of the algorithms' free parameters (hyperparameters). Therefore, it is often unclear which exact steps were
taken to produce the reported results and how robust the results are to changes in the free parameters or the data. Even if the authors are willing to provide
more details, essential information may not be accessible due to black-box (closed-source) analysis tools used somewhere in the process.

With `flowTorch`, we attempt to make analysis and modeling workflows accessible, streamlined, and transparent in several ways:

- we provide Jupyter notebooks with start-to-end workflows, including short explanations for each step taken in the process; the notebooks' content varies from
toy examples through common benchmark problems to the analysis of real turbulent flow data; the datasets used in the notebooks are also part of
the library
- the library is modular and often wraps around other libraries to make them easier to use; a few lines of Python code are sufficient to implement
a basic workflow; the modular structure and the rich documentation of the source code simplify writing extensions and enable quick automated experimentation

Ultimately, our goal is to reduce redundant work as much as possible and enable users to focus on what matters - understanding and modeling flow dynamics.

# Examples

In this section, we demonstrate two applications of `flowTorch`. In the first example, DMD is employed to identify relevant modes in a transonic flow displaying shock-boundary-layer interactions. In the second example, a ROM of the flow past a circular cylinder [@noack2003] is constructed employing CNM [@fernex2021]. Both examples are also available as Jupyter notebooks and in the `flowTorch` documentation.

## DMD analysis of airfoil surface data

For this example, we need only a handful of `flowTorch` components.
```
import torch as pt
from flowtorch import DATASETS
from flowtorch.data import CSVDataloader, mask_box
from flowtorch.analysis import DMD
```
`DATASETS` is a dictionary holding names and paths of all available datasets. The `CSVDataloader` provides easy access to the data, and the `mask_box` function allows selecting only a spatial subset of the raw data. As the name suggests, the `DMD` class enables us to perform a DMD analysis.

The dataset we use here consists of surface pressure coefficient distributions sampled over a NACA-0012 airfoil in transonic flow conditions. The OpenFOAM configuration files to produce the dataset are available in a separate [GitHub repository](https://github.com/AndreWeiner/naca0012_shock_buffet). At a Reynolds number of $Re=10^6$, a Mach number of $Ma=0.75$ and $\alpha = 4^\circ$ angle of attack, the flow displays a so-called shock buffet on the upper side of the airfoil. The shock buffet is a self-sustained unsteady interaction between the shock and the boundary layer separation. Our aim is to extract flow structures (modes) associated with the buffet phenomenon.

A code snippet to read the data, mask part of it, and build the data matrix reads:
```
...
path = DATASETS["csv_naca0012_alpha4_surface"]
loader = CSVDataloader.from_foam_surface(
    path, "total(p)_coeff_airfoil.raw", "cp")
vertices = loader.vertices
vertices /= (vertices[:, 0].max() - vertices[:, 0].min())
mask = mask_box(vertices,
    lower=[-1.0, 0.0, -1.0], upper=[0.9999, 1.0, 1.0])
points_upper = mask.sum().item()
data_matrix = pt.zeros((points_upper, len(times)), dtype=pt.float32)
for i, time in enumerate(times):
    snapshot = loader.load_snapshot("cp", time)
    data_matrix[:, i] = pt.masked_select(snapshot, mask)
```
The `CSVDataloader` has a class method designed to read raw sample data created by OpenFOAM simulations. Every `Dataloader` implementation provides access to one or multiple snapshots of one or multiple fields and the associated vertices. The airfoil coordinates are typically normalized with the chord length, which is the difference between largest and smallest value of the $x$-coordinate in the present example. The data contain pressure coefficients from both upper and lower surfaces, so we create a spatial mask to extract values from the upper surface. The DMD expects a data matrix as input whose columns are individual snapshots. Therefore, we allocate a new 2D tensor with as many rows as selected points and as many columns as selected snapshots (the data loader also provides access to the available write time - not shown here). Finally, we loop over the snapshot times and fill the data matrix.

Creating a new `DMD` instance automatically performs the mode decomposition based on the provided input. We can analyze
the obtained spectrum and the associated modes. The modes have real and imaginary parts, which are equally important for the
reconstruction of the flow field. It is usually enough to visualize either real or imaginary part for the physical interpretation of modes.
```
dmd = DMD(data_matrix, dt, rank=200)
amplitudes = dmd.amplitudes
frequencies = dmd.frequency
modes_real = dmd.modes.real
```
In contrast to POD, the DMD modes are not sorted by their variance, but rather form a spectrum.
\autoref{fig:dmd} presents the real part of three spatial modes with the largest amplitudes. Also shown is their corresponding frequency.

![Real part of three dominant DMD modes over the upper surface of a NACA-0012 airfoil. The modes are normalized to the range $[0,1]$. The coordinates are normalized with the chord $c$. The shock is located at $x/c\approx 0.25$. Modes 8 and 18 are harmonics. The motion of the shock front is correlated with changes in the pressure values close to the trailing edge. This effect can be nicely observed via the mode animations in the documentation and indicates the existence of a physical link between both effects.\label{fig:dmd}](media/dmd_modes_airfoil_cp.pdf){ width=90% }

## CNM of the flow past a circular cylinder

This example demonstrates how to model a flow using the CNM algorithm [@fernex2021]. Compared to the original CNM implementation available on [GitHub](https://github.com/fernexda/cnm), the version in `flowTorch` is refactored, more user-friendly, and extendible. In `flowTorch`, creating a ROM always consists of three steps: i) encoding/reduction, ii) time evolution, and iii) decoding/reconstruction. In the code snippet below, we use an encoder based on the singular value decomposition (SVD) to reduce the dimensionality of the original snapshot sequence, and then predict the temporal evolution and reconstruct the flow over the period of $1s$.

```
...
from flowtorch.rom import CNM, SVDEncoder
# load data
...
encoder = SVDEncoder(rank=20)
info = encoder.train(data_matrix)
reduced_state = encoder.encode(data_matrix)
cnm = CNM(reduced_state, encoder, dt, n_clusters=20, model_order=4)
prediction = cnm.predict(data_matrix[:, :5], end_time=1.0, step_size=dt)
```
The `predict` function computes the temporal evolution in the reduced state space and automatically performs the reconstruction. If we are only interested in the phase space, we can use `predict_reduced` instead, and reconstruct selected states using the encoder's `decode` method. The temporal evolution in the phase-space is displayed in \autoref{fig:cnm}.

![Phase-space representation of data clustering (large dots) and trajectory; the numbering reflects the sequence in which the centroids are visited; the smaller dots mark interpolated time steps between the centroids and are colored by their cluster affiliation (only for visualization).\label{fig:cnm}](media/cnm_cluster_transition.pdf){ width=70% }

# Acknowledgements

The authors gratefully acknowledge financial support by the German Research Foundation (DFG) received within the research unit [FOR 2895](https://www.for2895.uni-stuttgart.de/en/) *Unsteady flow and interaction phenomena at high speed stall conditions*.

# References![FOR2895Logo](media/for2895_logo.png)

# flowTorch

[![status](https://joss.theoj.org/papers/57b32d31997c90a40b3f4bdc20782e55/status.svg)](https://joss.theoj.org/papers/57b32d31997c90a40b3f4bdc20782e55)

**flowTorch** - a Python library for analysis and reduced order modeling of fluid flows

*flowTorch* is developed primarily by [@AndreWeiner](https://github.com/AndreWeiner) in the [Flow Modeling and Control group](https://www.tu-braunschweig.de/en/ism/research-workgroups/flow-modelling-and-control) led by [Richard Semaan](https://www.tu-braunschweig.de/en/ism/research/flow-modelling-and-control/staff/semaan). The development is financed by the German Research Foundation (DFG) within the research program [FOR 2895](https://www.for2895.uni-stuttgart.de/)


> unsteady flow and interaction phenomena at high speed stall conditions

with the primary goal to investigate flow conditions that lead to [buffeting](https://en.wikipedia.org/wiki/Aeroelasticity#Buffeting) at airfoils in the transonic flow regime. The animation below shows the shock buffet on a NACA-0012 airfoil at *Re=10^7*, *Ma=0.75*, and 4 degrees angle of attack. The simulation was conducted in OpenFOAM; follow [this link](https://github.com/AndreWeiner/naca0012_shock_buffet) for more information about the setup.

https://user-images.githubusercontent.com/8482575/120886182-f2b78800-c5ec-11eb-9b93-efb9a139c431.mp4

## Why *flowTorch*?

The *flowTorch* project was started to make the analysis and modeling of fluid data **easy** and **accessible** to everyone. The library design intends to strike a balance between **usability** and **flexibility**. Instead of a monolithic, black-box analysis tool, the library offers modular components that allow assembling custom analysis and modeling workflows with ease. *flowTorch* helps to fuse data from a wide range of file formats typical for fluid flow data, for example, to compare experiments simulations. The available analysis and modeling tools are rigorously tested and demonstrated on a variety of different fluid flow datasets. Moreover, one can significantly accelerate the entire process of accessing, cleaning, analysing, and modeling fluid flow data by starting with one of the pipelines available in the *flowTorch* [documentation](https://flowmodelingcontrol.github.io/flowtorch-docs/1.0/index.html).

To get a first impression of how working with *flowTorch* looks like, the code snippet below shows part of a pipeline for performing a dynamic mode decomposition (DMD) of a transient *OpenFOAM* simulation.

```
import torch as pt
from flowtorch import DATASETS
from flowtorch.data import FOAMDataloader, mask_box
from flowtorch.analysis.dmd import DMD

path = DATASETS["of_cylinder2D_binary"]
loader = FOAMDataloader(path)

# select a subset of the available snapshots
times = loader.write_times
window_times = [time for time in times if float(time) >= 4.0]

# load vertices, discard z-coordinate, and create a mask
vertices = loader.vertices[:, :2]
mask = mask_box(vertices, lower=[0.1, -1], upper=[0.75, 1])

# assemble the data matrix
data_matrix = pt.zeros((mask.sum().item(), len(window_times)), dtype=pt.float32)
for i, time in enumerate(window_times):
    # load the vorticity vector field, take the z-component [:, 2], and apply the mask
    data_matrix[:, i] = pt.masked_select(loader.load_snapshot("vorticity", time)[:, 2], mask)

# perform DMD
dmd = DMD(data_matrix, rank=19)
# analyse dmd.modes or dmd.eigvals
# ...
```

Currently, the following sub-packages are under active development. Note that some of the components are not yet available in the public release because further developments and testing are required:

| package | content |
| :------ | :-------|
|flowtorch.data | data loading, domain reduction (masked selection) |
| flowtorch.analysis | algorithms for dimensionality reduction, including *proper orthogonal decomposition* (POD), *dynamic mode decomposition* (DMD), autoencoders, and variants thereof |
| flowtorch.rom | reduced-order modeling using [cluster-based network models (CNM)](https://github.com/fernexda/cnm) |

*flowTorch* uses the [PyTorch](https://github.com/pytorch/pytorch) library as a backend for data structures, data types, and linear algebra operations on CPU and GPU. Some cool features of *flowTorch* include:

- data accessors return PyTorch tensors, which can be used directly within your favorite machine learning library, e.g., *PyTorch*, *SkLearn* or *Tensorflow*
- most algorithms run on CPU as well as on GPU
- mixed-precision operations (single/double); switching to single precision makes your life significantly easier when dealing with large datasets
- user-friendly Python library that integrates easily with popular tools and libraries like *Jupyterlab*, *Matplotlib*, *Pandas*, or *Numpy*
- a rich tutorial collection to help you getting started
- interfaces to common data formats like [OpenFOAM](https://www.openfoam.com/), [VTK](https://vtk.org/) (for Flexi and SU2), [TAU](https://www.dlr.de/as/desktopdefault.aspx/tabid-395/526_read-694/), [iPSP](https://www.dlr.de/as/en/desktopdefault.aspx/tabid-183/251_read-13334/), CSV (for DaVis PIV data and raw OpenFOAM output)

*flowTorch* can be also used easily in combination with existing Python packages for analysis and reduced-order modeling thanks to the interoperability between PyTorch and NumPy. Great examples are (by no means a comprehensive list):

- [PyDMD](https://github.com/mathLab/PyDMD) - Python Dynamic Mode Decomposition
- [PySINDy](https://github.com/dynamicslab/pysindy) - sparse identification of nonlinear dynamical systems from data

## Getting started

The easiest way to install *flowTorch* is as follows:
```
# install via pip
pip3 install git+https://github.com/FlowModelingControl/flowtorch
# to uninstall flowTorch, run
pip3 uninstall flowtorch
```
Alternatively, you can also clone the repository manually by running
```
git clone git@github.com:FlowModelingControl/flowtorch.git
```
and install the dependencies listed in *requirements.txt*:
```
pip3 install -r requirements.txt
```

To get an overview of what *flowTorch* can do for you, have a look at the [online documentation](https://flowmodelingcontrol.github.io/flowtorch-docs/1.0/index.html). The examples presented in the online documentation are also contained in this repository. In fact, the documentation is a static version of several [Jupyter labs](https://jupyter.org/) with start-to-end analyses. If you are interested in an interactive version of one particular example, navigate to `./docs/source/notebooks` and run `jupyter lab`. Note that to execute some of the notebooks, the **corresponding datasets are required**. The datasets can be downloaded [here](https://cloudstorage.tu-braunschweig.de/getlink/fiQUyeDFx3sg2T6LLHBQoCCx/datasets_29_10_2021.tar.gz) (~1.4GB). If the data are only required for unit testing, a reduced dataset may be downloaded [here](https://cloudstorage.tu-braunschweig.de/getlink/fiFZaHCgTWYeq1aZVg3hAui1/datasets_minimal_29_10_2021.tar.gz) (~384MB). Download the data into a directory of your choice and navigate into that directory. To extract the archive, run:
```
# full dataset
tar xzf datasets_29_10_2021.tar.gz
# reduced dataset
tar xzf datasets_minimal_29_10_2021.tar.gz
```
To tell *flowTorch* where the datasets are located, define the `FLOWTORCH_DATASETS` environment variable:
```
# add export statement to bashrc; assumes that the extracted 'datasets' or 'datasets_minimal'
# folder is located in the current directory
# full dataset
echo "export FLOWTORCH_DATASETS=\"$(pwd)/datasets/\"" >> ~/.bashrc
# reduced dataset
echo "export FLOWTORCH_DATASETS=\"$(pwd)/datasets_minimal/\"" >> ~/.bashrc
# reload bashrc
. ~/.bashrc
```

## Development
### Documentation

To build the flowTorch documentation, the following additional packages are required:
```
pip3 install sphinx sphinx_rtd_theme nbsphinx recommonmark
```
To build the HTML version of the API documentation, navigate to `./docs` and run:
```
make html
```

### Unit testing
All sub-packages contain unit tests, which require the installation of PyTest:
```
pip3 install pytest
```
Moreover, the flowTorch datasets must be downloaded and referenced as described in the previous section.
To run all unit tests of all sub-packages, execute:
```
pytest flowtorch
```
You can also execute all tests in a sub-package, e.g., data
```
pytest flowtorch/data
```
or run individual test modules, e.g.,
```
pytest flowtorch/data/test_FOAMDataloader.py
```

## Getting help

If you encounter any issues using *flowTorch* or if you have any questions regarding current and future development plans, please use the repository's [issue tracker](https://github.com/FlowModelingControl/flowtorch/issues). Consider the following steps before and when opening a new issue:

0. Have you searched for similar issues that may have been already reported? The issue tracker has a *filter* function to search for keywords in open issues.
1. Click on the green *New issue* button in the upper right corner and describe your problem as detailed as possible. The issue should state what **the problem** is, what the **expected behavior** should be, and, maybe, suggest a **solution**. Note that you can also attach files or images to the issue.
2. Select a suitable label from the drop-down menu called *Labels*.
3. Click on the green *Submit new issue* button and wait for a reply.

## Reference

If *flowTorch* aids your work, you may support our work by referencing the following software article:
```
@article{Weiner2021,
  doi = {10.21105/joss.03860},
  url = {https://doi.org/10.21105/joss.03860},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3860},
  author = {Andre Weiner and Richard Semaan},
  title = {flowTorch - a Python library for analysis and reduced-order modeling of fluid flows},
  journal = {Journal of Open Source Software}
}
```

## License

*flowTorch* is [GPLv3](https://en.wikipedia.org/wiki/GNU_General_Public_License)-licensed; refer to the [LICENSE](https://github.com/FlowModelingControl/flowtorch/blob/main/LICENSE) file for more information.

# How to contribute to *flowTorch*

We appreciate all efforts contributing to the *flowTorch* project, may it be bug-fixes, feature contributions, feature suggestions, additional examples, or other kinds of improvements. If you would like to contribute, you may consider the following steps.

## 0. Open a new issue

It is always useful to open a new issue as a first step. The issue helps the developers to plan and organize developments and to provide quick feedback on potential problems or existing solutions. For example, it might be that the bug you are reporting has already been fixed on a development branch or that someone is already working on a similar feature to the one you are suggesting. *flowTorch* is still a rather small project, so there are typically few open issues. Nonetheless, you should give your issue a suitable label to follow common best practices (e.g., feature, bug, documentation, ...).

## 1. Fork the *flowTorch* repository and create a new branch

The typical workflow of forking and branching is very well described in the [GitHub documentation](https://docs.github.com/en/get-started/quickstart/fork-a-repo).

## 2. Ensure code quality

*flowTorch* uses the [PyTorch library](https://pytorch.org/docs/stable/index.html) as backend for array-like data structures (tensors) and operations thereon. When implementing new features, try to rely as much as possible on the functionality offered by PyTorch instead of using NumPy, SciPy or similar libraries.

Most of the library contains [type hints](https://docs.python.org/3/library/typing.html). Type hints are not strictly necessary to run the code, but they make the lives of everybody much easier, so please use type hint in all parts of your code.

Python is a language that allows implementing operations with enormous complexity in a single line of code. Therefore, it is extremely important to provide a detailed documentation of new functionality containing all considerations the developer had in mind and also potential references or resources that were used as basis. *flowTorch* generates the documentation using [Sphinx](https://www.sphinx-doc.org/en/master/), and therefore, doc-stings should be formatted as [reStructuredText](https://docutils.sourceforge.io/rst.html).

We use PyLint to ensure proper code formatting. If VSCode is your editor of choice, have a look at the [documentation for linting](https://code.visualstudio.com/docs/python/linting) to setup automated code formatting.

## 3. Provide unit tests

If new features are added, accompanying unit tests should be provided. We use [PyTest](https://docs.pytest.org/en/6.2.x/) for testing. If the tests require additional datasets, please make sure that you have the permission to share the data such that the new data can be added to the *flowTorch* datasets in the next release. It might be necessary in some cases to create fake data (data that behave the same way as real data but that might be smaller and not protected).

## 4. Push changes and create a pull-request

To complete your contribution, create a new [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) for your feature or bug-fix. The changes will then be tested by someone with write access to the main repository before they are merged. If datasets are required, please provide a download link such that unit tests or examples can be executed.

**Thank you for considering to contribute to the *flowTorch* project!**# NACA 0012-34 airfoil simulation

## Geometry generation

The script *generate_airfoil.py* creates two files:
- **naca0012-34.stl**: geometry serving a input STL for *snappyHexMesh*
- **NACA-0012-34.pdf**: a simple plot comparing the profile data generated by the script against external profile data stored in *naca0012-34.csv*

## Simulation conditions

- Reynolds number $Re=3.0\times 10^6$
- angle of attack $\alpha = 3.5^\circ$
- Mach number $Ma = 0.73$

The Mach number is defined as the ratio between free-stream velocity $U_{in}$ and the speed of sound $c$:

$$ Ma = U_{in}/c \text{ with } c=\sqrt{\kappa R_s T_{in}} $$

Some useful expression are:
- specific gas  constant: $R_s = R/M = 8.3145 J/(mol K) / 28.9 g/mol  \cdot 1000 g/kg = 287.7 J/(kg K)$
- head capacity ratio of dry air at room temperature: $\kappa = 1.4$
- speed of sound: $c=\sqrt{1.4 \cdot 287.7 \cdot 293} = 343.5$
- transonic inlet velocity: $U_{in} = Ma \cdot c = 0.73\cdot 343.5 = 250.8$
- free-stream density: $\rho = p_{inf}/(R_s T_{in}) = 10^5 / (287.7 \cdot 293) = 1.186 kg/m^3$
- kinematic viscosity based on free-stream values: $\nu = \mu / \rho = 1.82\times 10^{-5} / 1.186 = 1.5342\times 10^{-5} m^2/s$
- characteristic length (chord length of airfoil): $L = Re \nu / U_{in} = 3\times 10^6 \cdot 1.5342\times 10^{-5} / 250.8 = 0.18352$