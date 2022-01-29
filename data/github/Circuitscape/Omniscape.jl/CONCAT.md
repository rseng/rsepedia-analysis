<img src='docs/src/assets/logo.png' width=300/>

| **Docs** | **Chat** | **Cite** | **Status** |
|:-----------------------------------------------------:|:------------------------------------:|:-----------:|:-------:|
| [![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.circuitscape.org/Omniscape.jl/stable) [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://docs.circuitscape.org/Omniscape.jl/latest) | [![gitter](https://badges.gitter.im/Circuitscape/Omniscape.jl.png)](https://gitter.im/Circuitscape/Omniscape.jl) | [![DOI](https://joss.theoj.org/papers/10.21105/joss.02829/status.svg)](https://doi.org/10.21105/joss.02829) <br> [![Omniscape Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/Omniscape)](https://pkgs.genieframework.com?packages=Omniscape) | [![Build Status](https://github.com/Circuitscape/Omniscape.jl/workflows/CI/badge.svg)](https://github.com/Circuitscape/Omniscape.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/Circuitscape/Omniscape.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Circuitscape/Omniscape.jl) | 

Omniscape.jl is built on [Circuitscape.jl](https://github.com/Circuitscape/Circuitscape.jl) and implements the Omniscape connectivity modeling algorithm to map omni-directional habitat connectivity. The Omniscape algorithm was developed by [McRae and colleagues](https://www.researchgate.net/publication/304842896_Conserving_Nature's_Stage_Mapping_Omnidirectional_Connectivity_for_Resilient_Terrestrial_Landscapes_in_the_Pacific_Northwest) in 2016. **Check out [the docs](https://circuitscape.github.io/Omniscape.jl/stable) for additional information.**

## Installation

**The latest version of Omniscape.jl requires Julia version 1.5.4 or greater**. You can install Julia [here](https://julialang.org/downloads/). Once installation is complete, open a Julia terminal and run the following code to install Omniscape.jl.
```julia
using Pkg; Pkg.add("Omniscape")
```
If you want to install the latest (unreleased) development version of Omniscape, you can get it by running:
```julia
using Pkg; Pkg.add(PackageSpec(name = "Omniscape", rev = "main"))
```

## Citing Omniscape.jl 

Please cite [Landau et al. (2021)](https://doi.org/10.21105/joss.02829) when using Omniscape.jl.
> Landau, V.A., V.B. Shah, R. Anantharaman, and K.R. Hall. 2021. Omniscape.jl: Software to compute omnidirectional landscape connectivity. *Journal of Open Source Software*, *6*(57), 2829.

Here's a bibtex entry:
```
@article{Landau2021,
  doi = {10.21105/joss.02829},
  url = {https://doi.org/10.21105/joss.02829},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {57},
  pages = {2829},
  author = {Vincent A. Landau and Viral B. Shah and Ranjan Anantharaman and Kimberly R. Hall},
  title = {Omniscape.jl: Software to compute omnidirectional landscape connectivity},
  journal = {Journal of Open Source Software}
}

```

Please be sure to also cite the [original work](https://www.researchgate.net/publication/304842896_Conserving_Nature's_Stage_Mapping_Omnidirectional_Connectivity_for_Resilient_Terrestrial_Landscapes_in_the_Pacific_Northwest) where the Omniscape algorithm was first described:
> McRae, B. H., K. Popper, A. Jones, M. Schindel, S. Buttrick, K. R. Hall, R. S. Unnasch, and J. Platt. 2016. Conserving Nature’s Stage: Mapping Omnidirectional Connectivity for Resilient Terrestrial Landscapes in the Pacific Northwest. *The Nature Conservancy*, Portland, Oregon.

## Contributing
Contributions in the form of pull requests are always welcome and appreciated. To report a bug or make a feature request, please [file an issue](https://github.com/Circuitscape/Omniscape.jl/issues/new). For general discussions and questions about usage, start a conversation on [gitter](https://gitter.im/Circuitscape/Omniscape.jl).

## Acknowledgments
Development of this software package was made possible by funding from [NASA's Ecological Forecasting program](https://appliedsciences.nasa.gov/what-we-do/ecological-forecasting) and the [Wilburforce Foundation](http://www.wilburforce.org/) through a project led by Kim Hall at The Nature Conservancy. This software package would not have been possible without Brad McRae (1966-2017), the visionary behind Circuitscape, the Omniscape algorithm, and several other software tools for assessing connectivity. Omniscape.jl is built on [Circuitscape.jl](https://github.com/Circuitscape/Circuitscape.jl), which was authored by Ranjan Anantharaman and Viral Shah, both of whom have been incredibly helpful in steering and guiding the development of Omniscape.jl. Kim Hall, Aaron Jones, Carrie Schloss, Melissa Clark, Jim Platt, and early Omniscape.jl users helped steer software development by providing valuable feedback and insight.
---
title: 'Omniscape.jl: Software to compute omnidirectional landscape connectivity'
tags:
  - julia
  - ecology
  - circuit theory
  - habitat connectivity
  - omniscape
  - circuitscape
authors:
  - name: Vincent A. Landau
    orcid: 0000-0001-9290-9438
    affiliation: 1
  - name: Viral B. Shah
    orcid: 0000-0001-9602-4012
    affiliation: 2
  - name: Ranjan Anantharaman
    ordic: 0000-0002-4409-3937
    affiliation: 3
  - name: Kimberly R. Hall
    orcid: 0000-0002-7802-3558
    affiliation: 4
affiliations:
 - name: Conservation Science Partners, Inc., Fort Collins, Colorado, United States
   index: 1
 - name: Julia Computing Inc., Cambridge, Massachusetts, United States
   index: 2
 - name: Massachusetts Institute of Technology, Cambridge, Massachusetts, United States
   index: 3
 - name: The Nature Conservancy, Lansing, Michigan, United States
   index: 4
date: 9 October 2020
bibliography: paper.bib
---

# Summary

Omniscape.jl is a software package that implements the Omniscape algorithm [@mcrae2016] to compute landscape connectivity. It is written in the Julia programming language [@bezanson2017] to be fast, scalable, and easy-to-use. Circuitscape.jl [@anantharaman2020], the package on which Omniscape.jl builds and expands, abstracts landscapes as two-dimensional electrical networks and solves for current flow. The current flow that results represents landscape connectivity. Omniscape.jl is novel in that it produces maps of "omni-directional" connectivity, which provide a spatial representation of connectivity between every possible pair of start and endpoints in the landscape. These maps can be used by researchers and landscape managers to understand and predict how ecological processes (e.g., animal movement, disease transmission, gene flow, and fire behavior) are likely to manifest in geographic space. Omniscape.jl makes use of Julia's native multi-threading, making it readily scalable and deployable to high performance compute nodes. More information on the broader Circuitscape project, which is home to Circuitscape.jl and Omniscape.jl, can be found at [circuitscape.org](https://circuitscape.org).

# Motivation

Modeling where and how ecological processes are connected provides valuable information for landscape management and researchers. A common output from connectivity modeling efforts is a map that provides a spatial representation of connectivity by showing likely paths of flow for ecological processes. A common method for identifying flow or movement corridors uses graph-theory to identify least-cost paths [@bunn2000]. The circuit-theoretic approach to connectivity modeling was a more recent innovation, and it has been gaining popularity over the past decade [@mcrae2006; @mcrae2008; @dickson2019]. In the circuit-theoretic approach, the landscape is abstracted as a network of current sources, grounds, and resistors. The resulting current flow through the electrical network is then related to the movement or flow intensity of the ecological process of interest. These models were first implemented in the Circuitscape software package [@shah2008], which was recently updated and rewritten as Circuitscape.jl [@anantharaman2020] in the Julia programming language. 

Circuitscape.jl is most often run in "pairwise" mode, where current flow is calculated between pairs of user-defined "cores," which are usually habitat patches. Results from this method can be highly sensitive to the location of cores. This can be problematic in cases where core location is arbitrary, or when there is uncertainty about where cores should be placed. The Omniscape algorithm [@mcrae2016] offers an alternative, "coreless" approach to pairwise Circuitscape.jl and computes omni-directional landscape connectivity by implementing Circuitscape.jl iteratively in a moving window. By precluding the need to identify and delineate discrete cores, the Omniscape algorithm also allows for a more detailed evaluation of connectivity within natural areas that may otherwise be defined as cores in Circuitscape. Code for Python was developed to implement the Omniscape algorithm in @mcrae2016, but a user-friendly software package was not available. To fill this need, we developed Omniscape.jl, an easy-to-use software package written in Julia. Omniscape.jl is useful for modeling connectivity in landscapes that do not have discrete cores, for example landscapes that are a combination of natural, semi-natural, and human-modified lands, or in cases where understanding connectivity within natural areas is of interest.

# The Omniscape Algorithm

Omniscape.jl works by applying Circuitscape.jl iteratively through the landscape in a circular moving window with a user-specified radius (\autoref{fig:window}). Omniscape.jl requires two basic spatial data inputs: a resistance raster, and a source strength raster. The resistance raster defines the traversal cost (measured in ohms) for every pixel in the landscape, that is, the relative cost for the ecological process of interest to move through each pixel. The source strength raster defines for every pixel the relative amount of current (measured in amperes) to be injected into that pixel. In the case of modeling animal movement, a pixel with a high source strength corresponds to relatively more individuals originating from that pixel.

![A diagram of the moving window used in Omniscape.jl, adapted with permission from @mcrae2016.\label{fig:window}](fig1.png)

The algorithm works as follows:

1. The circular window centers on a pixel in the source strength surface that has a source strength greater than 0 (or a user-specified threshold). This is referred to as the target pixel.
2. The source strength and resistance rasters are clipped to the circular window centered on the target pixel.
3. Every source strength pixel within the search radius that also has a source strength greater than 0 is identified. These are referred to as the source pixels.
4. Circuitscape.jl is run using the clipped resistance raster in “advanced" mode, where the target pixel is set to ground, and the source pixels are set as current sources. The total amount of current injected is equal to the source strength of the target pixel, and is divvied up among the source pixels in proportion to their source strengths.

Steps 1-4 are repeated for every potential target pixel. The resulting current maps from each moving window iteration are summed to get a final map of cumulative current flow. Individual moving window iterations can be run independently. Omniscape.jl makes use of Julia's multi-threaded parallel processing to solve individual moving windows in parallel.

In addition to cumulative current, Omniscape.jl also optionally provides two additional outputs: flow potential, and normalized cumulative current (\autoref{fig:outputs}). Flow potential represents current flow under "null" resistance conditions and demonstrates what current flow would look like when unconstrained by resistance and barriers. Flow potential is calculated exactly as cumulative current flow, but with resistance set to one for the entire landscape. Normalized cumulative current flow is calculated by dividing cumulative current flow by flow potential. Normalized current helps identify areas where current is impeded or channelized (e.g., more or less current than expected under null resistance conditions). High values mean current flow is channelized, and low values mean current is impeded.

![An example of the three different Omniscape.jl outputs. Outputs are from the Maryland forest connectivity example in the software documentation. Cumulative current flow shows the total current for each landscape pixel. Flow potential shows predicted current under resistance-free conditions. Normalized current shows the degree to which a pixel has more or less current than expected under resistance-free conditions (cumulative current flow divided by flow potential). Each layer is visualized using a quantile stretch with 30 breaks. Axes show easting and northing coordinates for reference.\label{fig:outputs}](fig2.png)

# Usage

Omniscape.jl is run from the Julia REPL. It offers a single user-facing function, `run_omniscape`, which has two methods. The first method accepts a single argument specifying the path to an [INI file](https://en.wikipedia.org/wiki/INI_file) that contains input data file paths and run options. Spatial data inputs can be in either ASCII or GeoTIFF raster formats, and outputs can also be written in either format. The second method of `run_omniscape` accepts arrays representing resistance and other spatial data inputs, and a dictionary of arguments specifying algorithm options. A complete user guide for Omniscape.jl, including installation instructions, function documentation, examples, and a complete list of options and their defaults, can be found in the [Omniscape.jl documentation](https://docs.circuitscape.org/Omniscape.jl/latest/).


# Acknowledgments
Development of this software package was made possible by funding from NASA's Ecological Forecasting program (grant NNX17AF58G) and the Wilburforce Foundation. This software package would not have been possible without Brad McRae (1966-2017), the visionary behind Circuitscape, the Omniscape algorithm, and several other software tools for assessing connectivity. Aaron Jones developed the diagram in \autoref{fig:window}. Aaron Jones, Carrie Schloss, Melissa Clark, Jim Platt, and early Omniscape.jl users helped steer software development by providing valuable feedback and insight.

# References
# User Guide

## Installation
The latest version of Omniscape.jl **requires Julia version 1.5.4 or greater**. You can install Julia [here](https://julialang.org/downloads/). Once installation is complete, open a Julia terminal and run the following code to install Omniscape.jl.
```julia
using Pkg; Pkg.add("Omniscape")
```
If you want to install the latest (unreleased) development version of Omniscape.jl, you can get it by running:
```julia
using Pkg; Pkg.add(PackageSpec(name = "Omniscape", rev = "main"))
```

Don't forget to "watch" the Omniscape.jl repository on GitHub to be notified of new releases! 

```@raw html
<td>
	<figure>
		<img src='../figs/watch_releases.gif' alt='How to be notified of new releases'><br>
		<figcaption><em></em></figcaption>
	</figure>
</td>
```

## Running Omniscape

For detailed examples of how to use Omniscape, check out the [Examples](@ref) section. Omniscape.jl provides a single user-facing function: [`run_omniscape`](@ref).

`run_omniscape()` offers two methods. The first, shown above, accepts the path to an [INI file](https://en.wikipedia.org/wiki/INI_file) that specifies file paths for raster inputs and other user-specified options. An INI file can be created using any text editor (e.g. notepad) and saved with the .ini file extension. The following code block shows an example INI file. The headings in square brackets are not required. They are there for organization purposes and are ignored by `run_omniscape()`.

```
[Required]
resistance_file = resistance_surface.tif
radius = 100
block_size = 5
project_name = output/example

[General options]
source_from_resistance = true
r_cutoff = 50
calc_normalized_current = true

parallelize = true
parallel_batch_size = 20

[Output options]
write_raw_currmap = true
```

The second method of `run_omniscape` accepts in-memory objects representing resistance and other spatial data inputs, and a dictionary specifying Omniscape settings and options.

The full suite of settings that are supported are described in detail in [Settings and Options](@ref) below.

### Parallel Processing

Omniscape uses parallel processing by default, but currently, Julia requires that the number of parallel threads to use be specified by an environment variable called `JULIA_NUM_THREADS`. This environment variable needs to be defined prior to launching Julia. The following examples demonstrate how to set `JULIA_NUM_THREADS` and start up Julia to use 4 threads from your terminal.

On Linux/Mac
```bash
export JULIA_NUM_THREADS=4
julia
```
On Windows
```bash
set JULIA_NUM_THREADS=4
julia
```
For Julia versions 1.5 and greater, you can now specify the number of threads with the `-t` flag when running Julia, e.g.:
```bash
julia -t 4
```

## Settings and Options

### Required

**`resistance_file`:** The path to the resistance layer input. This file can be in ASCII (.asc) or GeoTiff (.tif) format. If the file is in .asc format, Omniscape will also detect and use any associated .prj file in the same directory to determine the projection of the input file. The same applies for all other inputs described below that may be in .asc  format.

**`radius`:** A positive integer specifying the radius *in pixels* of the moving window used to identify sources to connect to each target.

**`project_name`:** The name of the project. Omniscape will create a directory called `project_name` in the directory from which you run Julia, and write any outputs to that directory. Supports the use of full path specification (e.g. path/to/directory).

*If `source_from_resistance` (described below) is false:*

**`source_file`:** The path to the source layer input. Must be an ASCII or GeoTIFF. This raster must have the same dimensions as `resistance_file`, and it is recommended that they have the exact same projection to ensure proper alignment. NoData values will be assigned a source strength of 0.  Does not need to be provided if `source_from_resistance` = true.

### Optional
#### General Options

**`block_size`:** An odd integer. Defaults to 1. An odd, positive integer specifying the side length for source layer blocking in target generation. The block size option coarsens the source strength surface for the purposes of identifying target pixels and assigning source strength values. The figure below shows two source strength grids. On the left is the case when `block_size = 1`. In this scenario, every pixel in the grid with a source strength greater than 0 (or, if specified, `source_threshold`, described below) will be a target pixel, and there will be a Circuitscape solve for each one of these pixels. The figure below and to the right represents the case when `block_size = 3`. In this case, the source strength grid is broken up into chunks of 9 pixels (3x3 blocks), each shown with a thick black outline. Only the centers of these 3x3 blocks (the white pixels) will be considered as potential target pixels. The maximum number of Circuitscape solves for this source grid is reduced from 81 when `block_size = 1` to just 9 when `block_size = 3`. In the `block_size = 3` case, the white pixels will be assigned a new source strength equal to the sum of the 9 pixels in the 3x3 block after setting any pixels with a source strength less than `source_threshold` to 0. This ensures that the total amount of current injected will be the same regardless of the value of `block_size`. Using a `block_size` > 1 can significantly reduce compute times and result in only negligable differences in the cumulative current map output.

```@raw html
<table border="0"><tr>
<td>
	<figure>
		<img src='../figs/sources_block_of_1.png' alt='missing'><br>
		<figcaption><em>Block size of 1</em></figcaption>
	</figure>
</td>
<td>
	<figure>
		<img src='../figs/sources_block_of_3.png' alt='missing'><br>
		<figcaption><em>Block size of 3</em></figcaption>
	</figure>
</td>
</tr></table>
```

**`source_from_resistance`**: One of true, false. Should a source layer be derived using the resistance layer? If true, sources are calculated as the inverse of the resistance layer, and therefore it is not recommended that your resistance layer contain values less than 1. Sources will be set to 0 for all cells with a resistance greater than `r_cutoff` (described below). Defaults to false.

**`resistance_is_conductance`:** One of true, false. Defaults to false. Specify whether the file specified by `resistance_file` is a conductance (rather than resistance) surface. Conductance is the inverse of resistance. Note that `r_cutoff` (an optional setting described below) must be in units of resistance even if a conductance file is supplied as input.

**`r_cutoff`**: The maximum resistance value a cell can have to be included as a source. Only applies when `source_from_resistance` = true. Defaults to Inf (which allows all cells to be considered as sources regardless of the resistance value).

**`buffer`**: A positive integer. Defaults to 0. Specifies an additional buffer distance beyond `radius` to clip the resistance and source layers to for each moving window iteration. Any source pixels beyond the `radius` but within the buffered area are set to 0. If 0, resistance and source layers will be clipped to a circle of size `radius` for each moving window iteration.

**`source_threshold`:** Positive number. Defaults to 0. Only pixels in the source layer greater than `source_threshold` will be included as sources.

**`calc_normalized_current`:** One of true, false. Defaults to false. Specify whether to calculate normalized current flow. Normalized current is calculated as raw current divided by flow potential. If true, a normalized current flow raster called "`normalized_cum_currmap`" (with either .tif or .asc extension, see `write_as_tif` below) will be written to the `project_name` directory.

**`calc_flow_potential`:** One of true, false. Defaults to false. Specify whether to calculate flow potential. Flow potential calculates current flow in "null" conditions, where the resistance of the entire landscape is 1. If true, a flow potential raster called "`flow_potential`" (with either .tif or .asc extension, see `write_as_tif` below) written to the `project_name` directory. This can still be set to false even if `calc_normalized_current` = true if you want to avoid writing the flow potential raster to disk.

**`allow_different_projections`:** One of true, false. Defaults to false. If true, warnings about non-matching projections are suppressed.

**`connect_four_neighbors_only`:** One of true, false. Defaults to false. Circuitscape creates a graph (network) by connecting cells to their four or eight immediate neighbors. The default is eight (four cardinal and four diagonal neighbors). Set `connect_four_neighbors_only` to true if you want to connect cells to their four cardinal neighbors only.

**`solver`:** One of "cg+amg" or "cholmod". Defaults to "cg+amg". The linear solver method to use in Circuitscape. See the Circuitscape.jl paper, [Anantharaman et al. (2019)](https://proceedings.juliacon.org/papers/10.21105/jcon.00058), for more information.

#### Resistance Reclassification
Omniscape.jl allows you to reclassify categorical resistance surfaces internally based on a user-provided reclass table. This allows the user to avoid reclassifying rasters manually in a GIS, and can streamline your workflow.

!!! note
    If instead of a resistance raster, you provide Omniscape a conductance raster, then conductance is what Omniscape will reclassify based on the provided reclass table.

**`reclassify_resistance`**: One of true, false. Defaults to false. Do you want Omniscape to reclassify your resistance/conductance raster using a reclass table that you provide?

**`reclass_table`**: If `reclassify_resistance = true`, the file path to the reclass table you wish to use. The reclass table is a two column, tab-separated .txt file. The first column contains the original resistance values in the resistance surface, and the second column specifies what those values should be changed to. You can reclassify values to `missing`
to replace them with infinite resistance (NoData). Note that you don't need to include an entry for every value in your original raster. If you only want to reclassify two specific resistance values, then only include entries for those two values.

Example reclass_table.txt:

```
1	3
2	5
3	1
4	2
5	missing
```

**`write_reclassified_resistance`**: One of true, false. Defaults to false. Should the reclassified resistance/conductance raster be saved to the output folder?


#### Processing Options
**`parallelize`:** One of true, false. Defaults to true. Specify whether to use parallel processing.

**`parallel_batch_size`:** Integer. Defaults to 10. The batch size (number of jobs) to send to each parallel worker. Particularly in cases where single solves are very fast, setting this to a larger number can reduce I/O overhead when scheduling/sending jobs to parallel workers. If set too high, then you will not be fully utilizing parallel workers.

**`precision`:** One of single, double. Defaults to double. Single precision uses less memory, but is less accurate than double precision. In certain cases (e.g. with extremelely large resistance values and/or extremely small source strengths), computations with single precision may be subject to [numerical underflow](https://en.wikipedia.org/wiki/Arithmetic_underflow), resulting in incorrect results. Use single precision with caution.

#### Output Options

**`write_raw_currmap`:** One of true, false. Defaults to true. Save the raw cumulative current map to disk? Should always be set to true unless `calc_flow_potential`, `calc_normalized_current`, or both are true and you do not need the raw current output. If true, the cumulative current map is saved to disk as a raster called "flow_potential" with either a .tif or .asc extension (see `write_as_tif` below) in the `project_name` directory.

**`mask_nodata`:** One of true, false. Defaults to true. Specify whether to mask current flow outputs according to NoData values in resistance surface. (i.e. pixels in current flow outputs that line up with NoData values in resistance are set to no data if `mask_nodata` = true).

**`write_as_tif`:** One of true, false. Defaults to true. Should outputs be written in tif format? If false, outputs are written in .asc format.


#### Conditional Connectivity Options

**`conditional`:** One of true, false. Defaults to false. Should conditional source/target matching be used? That is, should a given target only be connected to sources that meet similarity conditions to the target? If false, _none_ of the options described below are needed. If true, then gridded data with values for each pixel are used to compare targets and sources and determine which pairs should be connected according to user-specified criteria.

**`n_conditions`:** One of 1, 2. Defaults to 1. The number of conditions to use for conditional source/target matching.

*If `n_conditions` = 1:*

**`condition1_file`:** The file path to the data representing condition one in present day. Only needed if `conditional` = true. Must be an ASCII or GeoTIFF. This raster must have the same dimensions as `resistance_file`, and it is recommended that it also has the exact same projection to ensure proper alignment. Ensure that every pixel in the source strength raster has a corresponding value (not NoData) in `condition1_file`.

**`comparison1`:** One of within or equal. Defaults to within. How should conditions be compared when determining whether to connect a source/target pair. If within, then the value of condition 1 for the source must be within the following range, where target is the value at the target pixel or block: (target + `condition1_lower`, target + `condition1_upper`).  `condition1_lower` and `condition1_upper` are explained further below. If equal, then the value at the source pixel must be equal to the value at the target pixel.

**`condition1_lower`:** Number. Only required if `comparison1` = within. If `condition1_lower` = -1, then a source may have a condition 1 value up to 1 unit smaller than the target's value to be connected.

**`condition1_upper`:** Number. Only required if `comparison1` = within. If `condition1_upper` = 1, then a source may have a condition 1 value up to 1 unit larger than the target's value and it will still be connected.

*If `n_conditions` = 2:*

**`condition2_file`:** The file path to the data representing condition two in present day. Only needed if `conditional` = true and `n_conditions` = 2. Must be an ASCII or GeoTIFF. This raster must have the same dimensions as `resistance_file`, and it is recommended that it also has the exact same projection to ensure proper alignment. Ensure that every pixel in the source strength raster has a corresponding value (not NoData) in `condition2_file`.

**`comparison2`:** One of within or equal. Defaults to within. Only applies if `n_conditions` = 2. How should conditions be compared when determining whether to connect a source/target pair. If within, then the value of condition 2 for the source must be within the following range, where target is the value at the target pixel or block: (target + `condition2_lower`, target + `condition2_upper`).  `condition2_lower` and `condition2_upper` are explained further below. If equal, then the value at the source pixel must be equal to the value at the target pixel.

**`condition2_lower`:** Number. Only required if `n_conditions` = 2 and `comparison1` = within. If `condition2_lower` = -1, then a source may have a condition 2 value up to 1 unit smaller than the target's value and it will still be connected.

**`condition2_upper`:** Number. Only required if `n_conditions` = 2 and `comparison1` = within. If `condition2_lower` = 1, then a source may have a condition 2 value up to 1 unit larger than the target's value and it will still be connected.

*Using future conditions:*

**`compare_to_future`:** One of none, 1, 2, or both. Which condition(s) should compare the future condition in targets with present-day conditions in sources when determining which pairs to connect? For any condition(s) specified in this option, two data layers are needed: one with future condition values for all pixels in the study area, and one for present day condition values for all pixels in the study area. Defaults to "none".

**`condition1_future_file`:** The file path to the data representing condition one in the future. Only needed if `compare_to_future` = 1 or `compare_to_future` = both. Must be an ASCII or GeoTIFF. This raster must have the same dimensions as `resistance_file`, and it is recommended that they have the exact same projection to ensure proper alignment. Ensure that every pixel in the source strength raster has a corresponding value (not NoData) in `condition1_future_file`.

**`condition2_future_file`:** The file path to the data representing condition two in the future. Only needed if `n_conditions` = 2 *and* `compare_to_future` = 2 or `compare_to_future` = both. Must be an ASCII or GeoTIFF. This raster must have the same dimensions as `resistance_file`, and it is recommended that they have the exact same projection to ensure proper alignment. Ensure that every pixel in the source strength raster has a corresponding value (not NoData) in `condition2_future_file`.

## Omniscape in Docker

A Docker image with the latest version of Omniscape is [available on Docker Hub](https://hub.docker.com/r/vlandau/omniscape). To pull the image and start the Docker container from your terminal, navigate to the directory containing your Omniscape input files via `cd` and run the following code (set `JULIA_NUM_THREADS` to the number of threads you want to use for parallel processing):

On Linux/Mac:
```
docker run -it --rm \
	-v $(pwd):/home/omniscape \
	-w /home/omniscape \
	-e JULIA_NUM_THREADS=2 \
	vlandau/omniscape:latest
```

On Windows (via Windows Command Line):
```
docker run -it --rm^
 -v %cd%:/home/omniscape^
 -w /home/omniscape^
 -e JULIA_NUM_THREADS=2^
 vlandau/omniscape:latest
```
The `-v` flag and subsequent code will mount the files in your current working directory and make them available to the Docker container (which is why you need to run the code above from the directory that contains your input files). Once you're in Julia in the Docker container, you're ready to go! Make sure that the file paths in your .ini file are relative to the working directory from which you ran Docker.
# Omniscape.jl

Package repository: [https://github.com/Circuitscape/Omniscape.jl](https://github.com/Circuitscape/Omniscape.jl)

!!! note

    Before proceeding, it is strongly recommended that you familiarize yourself with the circuit theoretic approach to modeling landscape connectivity. See [McRae 2006](https://circuitscape.org/pubs/McRae_2006_IBR_Evolution.pdf) and [McRae et al. 2008](https://circuitscape.org/pubs/McRae_et_al_2008_Ecology.pdf) to learn more. See [Anantharaman et al. 2020](https://proceedings.juliacon.org/papers/10.21105/jcon.00058) for more on [Circuitscape.jl](https://github.com/Circuitscape/Omniscape.jl).

### Table of Contents

```@contents
Pages = ["index.md", "algorithm.md","usage.md", "examples.md", "apidocs.md"]
Depth = 2
```

## About Omniscape.jl

Omniscape.jl implements the Omniscape connectivity algorithm developed by [McRae et al. (2016)](https://www.researchgate.net/publication/304842896_Conserving_Nature's_Stage_Mapping_Omnidirectional_Connectivity_for_Resilient_Terrestrial_Landscapes_in_the_Pacific_Northwest). This software package can be used to produce maps of omni-directional habitat connectivity useful for scientific research as well as landscape management and conservation. Omniscape.jl is built on [Circuitscape.jl](https://github.com/Circuitscape/Circuitscape.jl). It offers a unique approach to connectivity modeling, particularly among circuit theoretic methods, by allowing the sources, destinations, and intensity of animal movement or ecological flow (modeled as electrical current) to be informed by continuous spatial data (such as a habitat suitability map). This information is combined with other spatial information on landscape resistance to movement or flow to produce models of habitat connectivity. To learn about how the algorithm works, see [The Omniscape Algorithm](@ref). Check out the [Examples](@ref) section for a step-by-step demonstration of how to use Omniscape.jl.

### Outputs

Omniscape.jl provides three different outputs.
1. **Cumulative current flow**: the total current flowing through the landscape -- the result of the Omniscape algorithm described above.
2. **Flow potential** (optional): current flow under "null" resistance conditions. Flow potential demonstrates what movement/flow would look like when movement is unconstrained by resistance and barriers. Flow potential is calculated exactly as cumulative current flow is, but with resistance set to 1 for the entire landscape.
3. **Normalized current flow** (optional): calculated as cumulative current flow divided by flow potential. Normalized current helps identify areas where current is impeded or channelized (e.g. more or less current than expected under null resistance conditions). High values mean current flow is channelized, low values mean current is impeded.

### Climate Connectivity

Climate connectivity can be modeled using the conditional connectivity options in Omniscape. These options options allow the user to impose extra constraints on source and target identification and matching. For example the present day climate of the source pixels might be required to be similar to the projected future climate for the target pixel. Info on constraints is provided to Omniscape via raster layers. See the documentation on [Conditional Connectivity Options](@ref) for more info on how to implement this feature.

## Citing Omniscape.jl

Please cite [Landau et al. (2021)](https://doi.org/10.21105/joss.02829) when using Omniscape.jl.
> Landau, V.A., V.B. Shah, R. Anantharaman, and K.R. Hall. 2021. Omniscape.jl: Software to compute omnidirectional landscape connectivity. *Journal of Open Source Software*, *6*(57), 2829.

Here's a bibtex entry:
```
@article{Landau2021,
  doi = {10.21105/joss.02829},
  url = {https://doi.org/10.21105/joss.02829},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {57},
  pages = {2829},
  author = {Vincent A. Landau and Viral B. Shah and Ranjan Anantharaman and Kimberly R. Hall},
  title = {Omniscape.jl: Software to compute omnidirectional landscape connectivity},
  journal = {Journal of Open Source Software}
}

```

Be sure to also cite the [original work](https://www.researchgate.net/publication/304842896_Conserving_Nature's_Stage_Mapping_Omnidirectional_Connectivity_for_Resilient_Terrestrial_Landscapes_in_the_Pacific_Northwest) where the Omniscape algorithm was first described:
> McRae, B. H., K. Popper, A. Jones, M. Schindel, S. Buttrick, K. R. Hall, R. S. Unnasch, and J. Platt. 2016. Conserving Nature’s Stage: Mapping Omnidirectional Connectivity for Resilient Terrestrial Landscapes in the Pacific Northwest. *The Nature Conservancy*, Portland, Oregon.

## Acknowledgments
Development of this software package was made possible by funding from NASA's Ecological Forecasting program and the Wilburforce Foundation through a project led by Kim Hall at The Nature Conservancy. This software package would not have been possible without Brad McRae (1966-2017), the visionary behind Circuitscape, the Omniscape algorithm, and several other software tools for assessing connectivity. Omniscape.jl is built on [Circuitscape.jl](https://github.com/Circuitscape/Circuitscape.jl), which was authored by Ranjan Anantharaman and Viral Shah, both of whom have been incredibly helpful in steering and guiding the development of Omniscape.jl. Kim Hall, Aaron Jones, Carrie Schloss, Melissa Clark, Jim Platt, and early Omniscape.jl users helped steer software development by providing valuable feedback and insight.


## References

Anantharaman, R., Hall, K., Shah, V., & Edelman, A. (2020). Circuitscape in Julia: Circuitscape in Julia: High Performance Connectivity Modelling to Support Conservation Decisions. *Proceedings of the JuliaCon Conferences*. DOI: 10.21105/jcon.00058.

McRae, B. H. (2006). Isolation by resistance. *Evolution*, 60(8), 1551-1561.

McRae, B. H., Dickson, B. G., Keitt, T. H., & Shah, V. B. (2008). Using circuit theory to model connectivity in ecology, evolution, and conservation. *Ecology*, 89(10), 2712-2724.

McRae, B. H., Popper, K., Jones, A., Schindel, M., Buttrick, S., Hall, K., Unnasch, B. & Platt, J. (2016). Conserving nature’s stage: mapping omnidirectional connectivity for resilient terrestrial landscapes in the Pacific Northwest. *The Nature Conservancy*, Portland, Oregon.
# API Documentation
```@docs
run_omniscape

missingarray

missingarray_to_array
```
# How It Works

## The Omniscape Algorithm

The Omniscape algorithm works by applying Circuitscape iteratively through the landscape in a moving window with a user-specified radius. Omniscape requires two basic inputs: a resistance raster, and a source strength raster. The resistance raster defines the traversal cost for every pixel in the landscape. The source strength raster defines for every pixel the relative amount of current to be injected into that pixel. A diagram of the moving window, adapted and borrowed from McRae et al. 2016, is shown in figure 1 below.

```@raw html
<img src='../figs/moving_window.png' width=350)> <br><em><b>Figure 1</b>: An illustration of a moving window iteration in the Omniscape algorithm.</em><br><br>
```

The algorithm works as follows:
1. The window centers on a pixel in the source strength surface that has a source strength greater than 0 (or a user specified threshold). This is referred to as the target pixel.
2. The source strength and resistance rasters are clipped to the circular window.
3. Every source pixel within the search radius that also has a source strength greater than 0 is identified. These are referred to as the source pixels.
4. Circuitscape is run using the clipped resistance raster in “advanced” mode, where the target pixel is set to ground, and the source pixels are set as current sources. The total amount of current injected is equal to the source strength of the target pixel, and is divvied up among the source pixels in proportion to their source strengths.
5. Steps 1-4 are repeated for every potential target pixel. The resulting current maps are summed to get a map of cumulative current flow.

The Omniscape algorithm evaluates connectivity between every possible pair of pixels in the landscape that are a) valid sources (i.e. have a source strength greater than 0 or other user-specified threshold) and b) no further apart than the moving window radius.
# Examples

## Forest connectivity in central Maryland

Land cover datasets are commonly used to parameterize resistance for connectivity modeling. This example uses the [National Land Cover Dataset](https://www.usgs.gov/centers/eros/science/national-land-cover-database) for the United States to model forest connectivity in central Maryland. Each value in the categorical land cover dataset is assigned a resistance score. We can have Omniscape.jl assign these values internally by providing a reclassification table (see [Resistance Reclassification](@ref)).

First, install the necessary packages and import them:

```julia
using Pkg; Pkg.add(["Omniscape", "Rasters", "Plots"])
using Omniscape, Rasters, Plots
```
```@setup mdforest
using Pkg; Pkg.add(["Omniscape", "Rasters", "Plots"])
using Omniscape, Rasters, Plots
url_base = "https://raw.githubusercontent.com/Circuitscape/datasets/main/"
download(string(url_base, "data/nlcd_2016_frederick_md.tif"),
         "nlcd_2016_frederick_md.tif")
nothing
```

Next, download the landcover data we'll use in this example, and plot it:

```julia
url_base = "https://raw.githubusercontent.com/Circuitscape/datasets/main/"
# Download the NLCD tile used to create the resistance surface and load it
download(string(url_base, "data/nlcd_2016_frederick_md.tif"),
         "nlcd_2016_frederick_md.tif")

# Plot the landcover data
values = [11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95]
palette = ["#476BA0", "#DDC9C9", "#D89382", "#ED0000", "#AA0000",
           "#b2b2b2", "#68AA63", "#1C6330", "#B5C98E", "#CCBA7C",
           "#E2E2C1", "#DBD83D", "#AA7028", "#BAD8EA", "#70A3BA"]

plot(Raster("nlcd_2016_frederick_md.tif"),
     title = "Land Cover Type", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(palette, (values .- 12) ./ 84, categorical = true),
     size = (700, 640))
```
```@raw html
<img src='../figs/mdlc.png' width=500><br>
```

Now, load the array using Omniscape's internal `read_raster()` function or a function from a GIS Julia package of your choice. `read_raster()` returns a tuple with the data array, a wkt string containing geographic projection info, and an array containing geotransform values. We'll use the wkt and geotransform later.

```@example mdforest
land_cover, wkt, transform = Omniscape.read_raster("nlcd_2016_frederick_md.tif", Float64)
```

The next step is to create a resistance reclassification table that defines a resistance value for each land cover value. Land cover values go in the left column, and resistance values go in the right column. In this case, we are modeling forest connectivity, so forest classes receive the lowest resistance score of one. Other "natural" land cover types are assigned moderate values, and human-developed land cover types are assigned higher values. Medium- to high-intensity development are given a value of `missing`, which denotes infinite resistance (absolute barriers to movement).

```@example mdforest
# Create the reclassification table used to translate land cover into resistance
reclass_table = [
    11.	100; # Water
    21	500; # Developed, open space
    22	1000; # Developed, low intensity
    23	missing; # Developed, medium intensity
    24	missing; # Developed, high intensity
    31	100; # Barren land
    41	1; # Deciduous forest
    42	1; # Evergreen forest
    43	1; # Mixed forest
    52	20; # Shrub/scrub
    71	30; # Grassland/herbaceous
    81	200; # Pasture/hay
    82	300; # Cultivated crops
    90	20; # Woody wetlands
    95	30; # Emergent herbaceous wetlands
]
```

Next, we define the configuration options for this model run. See the [Settings and Options](@ref) section in the [User Guide](@ref) for more information about each of the configuration options.

```@example mdforest
# Specify the configuration options
config = Dict{String, String}(
    "radius" => "100",
    "block_size" => "21",
    "project_name" => "md_nlcd_omniscape_output",
    "source_from_resistance" => "true",
    "r_cutoff" => "1", # Only forest pixels should be sources
    "reclassify_resistance" => "true",
    "calc_normalized_current" => "true",
    "calc_flow_potential" => "true"
)
```

Finally, compute connectivity using `run_omniscape()`, feeding in the configuration dictionary, the resistance array, the reclass table, as well as the wkt and geotransform information loaded earlier. Passing in the wkt and geotransform, along with `true` for the `write_outputs` argument, will allow Omniscape to write the outputs as properly projected rasters. `run_omniscape` will print some information to the console and show progress, along with an ETA, in the form of a progress bar.

```@example mdforest
currmap, flow_pot, norm_current = run_omniscape(config,
                                                land_cover,
                                                reclass_table = reclass_table,
                                                wkt = wkt,
                                                geotransform = transform,
                                                write_outputs = true)
```

You'll see that outputs are written to a new folder called "md\_nlcd\_omniscape\_output". This is specified by the "project\_name" value in `config` above. The cumulative current map will always be called "cum\_currmap.tif", and it will be located in the output folder. We also specified in the run configuration that flow potential and normalized current should be computed as well. These are called "flow\_potential.tif" and "normalized\_cum\_currmap.tif", respectively. See [Outputs](@ref) for a description of each of these outputs.

Now, plot the outputs. Load the outputs into Julia as spatial data and plot them.

First, the cumulative current map:

```julia
current = Raster("md_nlcd_omniscape_output/cum_currmap.tif")
plot(current,
     title = "Cumulative Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (600, 550))
```
```@raw html
<img src='../figs/md-curmap.png' width=500> <br><em>Cumulative current flow representing forest connectivity. Note that areas in white correspond to built up areas (NLCD values of 23 and 24) that act as absolute barriers to movement.</em><br><br>
```

Next, flow potential. This map shows what connectivity looks like under "null" conditions (resistance equals 1 for the whole landscape).

```julia
fp = Raster("md_nlcd_omniscape_output/flow_potential.tif")
plot(fp,
     title = "Flow Potential", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))
```
```@raw html
<img src='../figs/md-fp.png' width=500> <br><em>Flow potential, which shows what connectivity would look like in the absence of barriers to movement. The blocking that you can see is an artifact of setting a large block_size to make the example run faster. Set a smaller block_size to reduce/remove this issue.</em><br><br>
```

Finally, map normalized current flow, which is calculated as flow potential divided by cumulative current.

```julia
normalized_current = Raster("md_nlcd_omniscape_output/normalized_cum_currmap.tif")
plot(normalized_current,
     title = "Normalized Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))
```
```@raw html
<img src='../figs/md-norm-cur.png' width=500> <br><em>Normalized cumulative current. Values greater than one indicate areas with channelized/bottlenecked flow. Values around 1 (cumulative current ≈ flow potential) indicate diffuse current. Values less than 1 indicate impeded flow.</em><br><br>
```
