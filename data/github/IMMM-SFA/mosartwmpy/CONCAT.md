[![build](https://github.com/IMMM-SFA/mosartwmpy/actions/workflows/build.yml/badge.svg)](https://github.com/IMMM-SFA/mosartwmpy/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/IMMM-SFA/mosartwmpy/branch/main/graph/badge.svg?token=IPOY8984MB)](https://codecov.io/gh/IMMM-SFA/mosartwmpy)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03221/status.svg)](https://doi.org/10.21105/joss.03221)
[![DOI](https://zenodo.org/badge/312114600.svg)](https://zenodo.org/badge/latestdoi/312114600)

## mosartwmpy

`mosartwmpy` is a python translation of MOSART-WM, a model for water routing and reservoir management written in Fortran. The original code can be found at [IWMM](https://github.com/IMMM-SFA/iwmm) and [E3SM](https://github.com/E3SM-Project/E3SM), in which MOSART is the river routing component of a larger suite of earth-science models. The motivation for rewriting is largely for developer convenience -- running, debugging, and adding new capabilities were becoming increasingly difficult due to the complexity of the codebase and lack of familiarity with Fortran. This version aims to be intuitive, lightweight, and well documented, while still being highly interoperable.
For a quick start, check out the [Jupyter notebook tutorial](https://github.com/IMMM-SFA/mosartwmpy/blob/main/notebooks/tutorial.ipynb)!

## getting started

Ensure you have Python >= 3.7 available (consider using a [virtual environment](https://github.com/pyenv/pyenv), see the docs [here](https://mosartwmpy.readthedocs.io/en/latest/virtualenv.html) for a brief tutorial), then install `mosartwmpy` with:
```shell
pip install mosartwmpy
```

Alternatively, install via conda with:
```shell
conda install -c conda-forge mosartwmpy
```

Download a sample input dataset spanning May 1981 by running the following and selecting option `1` for "tutorial". This will download and unpack the inputs to your current directory. Optionally specify a path to download and extract to instead of the current directory.

```shell
python -m mosartwmpy.download
```

Settings are defined by the merger of the `mosartwmpy/config_defaults.yaml` and a user specified file which can override any of the default settings. Create a `config.yaml` file that defines your simulation (if you chose an alternate download directory in the step above, you will need to update the paths to point at your data):

> `config.yaml`
> ```yaml
> simulation:
>   name: tutorial
>   start_date: 1981-05-24
>   end_date: 1981-05-26
>
> grid:
>   path: ./input/domains/mosart.nc
>   land:
>     path: ./input/domains/land.nc
>
> runoff:
>   read_from_file: true
>   path: ./input/runoff/runoff_1981_05.nc
>
> water_management:
>   enabled: true
>   demand:
>     read_from_file: true
>     path: ./input/demand/demand_1981_05.nc
>   reservoirs:
>     enable_istarf: true
>     parameters:
>       path: ./input/reservoirs/reservoirs.nc
>     dependencies:
>       path: ./input/reservoirs/dependency_database.parquet
>     streamflow:
>       path: ./input/reservoirs/mean_monthly_reservoir_flow.parquet
>     demand:
>       path: ./input/reservoirs/mean_monthly_reservoir_demand.parquet
> ```

`mosartwmpy` implements the [Basic Model Interface](https://csdms.colorado.edu/wiki/BMI) defined by the CSDMS, so driving it should be familiar to those accustomed to the BMI. To launch the simulation, open a python shell and run the following:

```python
from mosartwmpy import Model

# path to the configuration yaml file
config_file = 'config.yaml'

# initialize the model
mosart_wm = Model()
mosart_wm.initialize(config_file)

# advance the model one timestep
mosart_wm.update()

# advance until the `simulation.end_date` specified in config.yaml
mosart_wm.update_until(mosart_wm.get_end_time())
```

## model input

Input for `mosartwmpy` consists of many files defining the characteristics of the discrete grid, the river network, surface and subsurface runoff, water demand, and dams/reservoirs.
Currently, the gridded data is expected to be provided at the same spatial resolution.
Runoff input can be provided at any time resolution; each timestep will select the runoff at the closest time in the past.
Currently, demand input is read monthly but will also pad to the closest time in the past.
Efforts are under way for more robust demand handling.
Dams/reservoirs require four different input files: the physical characteristics, the average monthly flow expected during the simulation period, the average monthly demand expected during the simulation period, and a database mapping each GRanD ID to grid cell IDs allowed to extract water from it.
These dam/reservoir input files can be generated from raw GRanD data, raw elevation data, and raw ISTARF data using the [provided utility](mosartwmpy/utilities/CREATE_GRAND_PARAMETERS.md).
The best way to understand the expected format of the input files is to examine the sample inputs provided by the download utility: `python -m mosartwmpy.download`.

## model output

By default, key model variables are output on a monthly basis at a daily averaged resolution to `./output/<simulation name>/<simulation name>_<year>_<month>.nc`. See the configuration file for examples of how to modify the outputs, and the `./mosartwmpy/state/state.py` file for state variable names.

Alternatively, certain model outputs deemed most important can be accessed using the BMI interface methods. For example:
```python
from mosartwmpy import Model

mosart_wm = Model()
mosart_wm.initialize()

# get a list of model output variables
mosart_wm.get_output_var_names()

# get the flattened numpy.ndarray of values for an output variable
supply = mosart_wm.get_value_ptr('supply_water_amount')
```

## visualization

`Model` instances can plot the current value of certain input and output variables (those available from `Model.get_output_var_name` and `Model.get_input_var_names`):

```python
from mosartwmpy import Model
config_file = 'config.yaml'
mosart_wm = Model()
mosart_wm.initialize(config_file)
for _ in range(8):
    mosart_wm.update()

mosart_wm.plot_variable('outgoing_water_volume_transport_along_river_channel', log_scale=True)
```
![River transport](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/river_transport.png)

Using provided utility functions, the output of a simulation can be plotted as well.

Plot the storage, inflow, and outflow of a particular GRanD dam:
```python
from mosartwmpy import Model
from mosartwmpy.plotting.plot import plot_reservoir
config_file = 'config.yaml'
mosart_wm = Model()
mosart_wm.initialize(config_file)
mosart_wm.update_until()

plot_reservoir(
    model=mosart_wm,
    grand_id=310,
    start='1981-05-01',
    end='1981-05-31',
)
```
![Grand Coulee](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/grand_coulee_1981_05.png)

Plot a particular output variable (as defined in `config.yaml`) over time:
```python
from mosartwmpy import Model
from mosartwmpy.plotting.plot import plot_variable
config_file = 'config.yaml'
mosart_wm = Model()
mosart_wm.initialize(config_file)
mosart_wm.update_until()

plot_variable(
    model=mosart_wm,
    variable='RIVER_DISCHARGE_OVER_LAND_LIQ',
    start='1981-05-01',
    end='1981-05-31',
    log_scale=True,
    cmap='winter_r',
)
```
![River network no tiles](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/river_without_tiles_1981_05.png)

If `cartopy`, `scipy`, and `geoviews` are installed, tiles can be displayed along with the plot:
```python
plot_variable(
    model=mosart_wm,
    variable='RIVER_DISCHARGE_OVER_LAND_LIQ',
    start='1981-05-01',
    end='1981-05-31',
    log_scale=True,
    cmap='winter_r',
    tiles='StamenWatercolor'
)
```
![River network with tiles](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/river_with_tiles_1981_05.png)

## model coupling

A common use case for `mosartwmpy` is to run coupled with output from the Community Land Model (CLM). To see an example of how to drive `mosartwmpy` with runoff from a coupled model, check out the [Jupyter notebook tutorial](https://github.com/IMMM-SFA/mosartwmpy/blob/main/notebooks/tutorial.ipynb)!

## testing and validation

Before running the tests or validation, make sure to download the "sample_input" and "validation" datasets using the download utility `python -m mosartwmpy.download`.

To execute the tests, run `./test.sh` or `python -m unittest discover mosartwmpy/tests` from the repository root.

To execute the validation, run a model simulation that includes the years 1981 - 1982, note your output directory, and then run `python -m mosartwmpy.validate` from the repository root. This will ask you for the simulation output directory, think for a moment, and then open a figure with several plots representing the NMAE (Normalized Mean Absolute Error) as a percentage and the spatial sums of several key variables compared between your simulation and the validation scenario. Use these plots to assist you in determining if the changes you have made to the code have caused unintended deviation from the validation scenario. The NMAE should be 0% across time if you have caused no deviations. A non-zero NMAE indicates numerical difference between your simulation and the validation scenario. This might be caused by changes you have made to the code, or alternatively by running a simulation with different configuration or parameters (i.e. larger timestep, fewer iterations, etc). The plots of the spatial sums can assist you in determining what changed and the overall magnitude of the changes.

If you wish to merge code changes that intentionally cause significant deviation from the validation scenario, please work with the maintainers to create a new validation dataset.
[Issues]: https://github.com/IMMM-SFA/mosartwmpy/issues
[new Issue]: https://github.com/IMMM-SFA/mosartwmpy/issues/new/choose
[new Pull Request]: https://github.com/IMMM-SFA/mosartwmpy/compare
[IM3]: https://im3.pnnl.gov/

## How to contribute
___

__I found a bug!__
* Nice sleuthing!
* Check if the bug has already been reported by searching the existing GitHub [Issues]. If you find a match, consider adding additional details to the existing ticket.
* Open a [new Issue], being sure to include a clear title and description along with as much detail as possible; code samples or log messages demonstrating the bug are quite helpful.

__I fixed a bug!__
* First of all, thanks!
* Open a [new Pull Request] with the fix. Ensure the description clearly outlines the bug and the solution. Include the Issue number if applicable.

__I created a new feature!__
* You're the best!
* Consider opening a [new Issue] to describe use cases for the new feature. This will offer a platform for discussion and critique.
* Then, open a [new Pull Request] with clear documentation of the methodology. Be sure to include new unit tests if appropriate.

The [IM3] team truly appreciates and encourages community involvement. We're all in this together! ---
title: 'mosartwmpy: A Python implementation of the MOSART-WM coupled hydrologic routing and water management model'
tags:
  - Python
  - hydrology
  - water management
  - multisector dynamics
  - reservoir modeling
authors:
  - name: Travis Thurber
    orcid: 0000-0002-4370-9971
    affiliation: 1
  - name: Chris R. Vernon
    orcid: 0000-0002-3406-6214
    affiliation: 1
  - name: Ning Sun
    orcid: 0000-0002-4094-4482
    affiliation: 1
  - name: Sean W. D. Turner
    orcid: 0000-0003-4400-9800
    affiliation: 1
  - name: Jim Yoon
    orcid: 0000-0002-8025-2587
    affiliation: 1
  - name: Nathalie Voisin
    orcid: 0000-0002-6848-449X
    affiliation: 1
affiliations:
 - name: Pacific Northwest National Laboratory, Richland, WA., USA
   index: 1
date: 22 March 2021
bibliography: paper.bib
---

# Summary
`mosartwmpy` is a Python implementation of the Model for Scale Adaptive River Transport with Water Management (MOSART-WM). This new version retains the functionality of the legacy model (written in FORTRAN) while providing new features to enhance user experience and extensibility. MOSART is a large-scale river-routing model used to study riverine dynamics of water, energy, and biogeochemistry cycles across local, regional, and global scales [@li2013physically]. The WM component introduced by @voisin2013improved represents river regulation through reservoir storage and release operations, diversions from reservoir releases, and allocation to sectoral water demands. Each reservoir release is independently calibrated using long-term mean monthly inflow into the reservoir, long-term mean monthly demand associated with this reservoir, and reservoir goals (flood control, irrigation, recreation, etc.). Generic monthly pre-release rules and storage targets are set up for individual reservoirs; however, those releases are updated annually for inter-annual variability (dry or wet year) and daily for environmental constraints such as flow minimum release and minimum/maximum storage levels. The WM model allows an evaluation of the impact of water management over multiple river basins at once (global, continental scales) and with consistent representation of human operations over the full domain.

# Statement of Need
MOSART-WM is often utilized as the hydrological component in a larger suite of earth-science models, such as in @doecode_10475. In this context, MOSART-WM is quite efficient and streamlined when running on a supported High-Performance Computing (HPC) cluster. However, learning how to use, extend, and test a complex codebase written with domain knowledge implied in a lower-level programming language may greatly increase the turnaround time and error rate for setting up and executing novel experiments. Broadening the code’s accessibility using a programming language such as Python, which in 2020 was the second most utilized language on GitHub [@octoverse_2020], provides a more accessible option to learn and contribute to the science on most computational platforms.

`mosartwmpy` was designed to bridge the gap between the domain scientist who wants a performant software that can still be extended for future research needs, and the new user who may not have expertise within the hydrologic sciences but wishes to integrate the process into their own workflow or quickly become capable in conducting hands-on experimentation for educational purposes. A refactor of MOSART-WM in Python ameliorates the steep learning curve of the FORTRAN version by providing an easy to learn, use, and modify interface. `mosartwmpy` was also built with interoperability in mind by implementing the Community Surface Dynamics Modeling System (CSDMS) Basic Model Interface (BMI) [@peckham2013component; @hutton2020basic], which offers a familiar set of controls for operating the model. By leveraging the BMI, `mosartwmpy` can be readily coupled with other earth system models to perform cross-domain experiments.

The target audience for `mosartwmpy` is the data scientist or hydrologist who wishes to rapidly prototype, test, and develop novel modeling capabilities relating to reservoirs and water demand. `mosartwmpy` can be accessed on GitHub (https://github.com/IMMM-SFA/mosartwmpy), and a walkthrough of key functionality and use can be found here: [tutorial](https://mosartwmpy.readthedocs.io/en/latest/). Model results have been validated against historical simulations [@https://doi.org/10.1029/2020WR027902], and a validation utility is included with the code.

# State of the field
In addition to `mosartwmpy`'s ancestor MOSART-WM, several other models are commonly used in hydrologic modeling, each excelling at different processes. The Community Land Model [CLM; @https://doi.org/10.1029/2018MS001583] focuses on the traditional water cycle (i.e. precipitation, evaporation) as well as plant and soil chemistry; runoff output from CLM can be used as input for `mosartwmpy`. StateMod [@statemod] focuses on water allocation based on legal constraints (water rights) as well as supply and demand. MODFLOW [@modflow] focuses on solving for complex groundwater flow properties on three-dimensional grids. GLOFRIM [@hoch2017glofrim] focuses on providing a framework to couple models implementing the BMI across multi-scale grids, for instance by coupling a meteorological model to a water cycle model to river routing model to hydrodynamic model. In this context, `mosartwmpy` focuses on the interactions between river routing, reservoir management, and water demand.

# Functionality
Model input for `mosartwmpy` consists of channel and reservoir geometry, groundwater and subsurface runoff (i.e. rain and ice melt), and water demand. The runoff populates the local streams and river channels, which generally flow toward the oceans (as in \autoref{fig:flow}). In many locations along the main channels, reservoirs collect and store water, releasing portions of the storage downstream over time for various societal uses (i.e. hydroelectricity, drinking water, flood control, and irrigation). Water is also consumed directly from the channels and tributaries for irrigation or non-irrigation uses.

![River basin flow over the continental United States as output from `mosartwmpy`.\label{fig:flow}](figure_1.png){ width=75% }

# Ongoing Research
As open-source software, `mosartwmpy` promotes collaborative and community model development to meet scientific objectives. Two examples of the potential to extend this codebase are currently underway as a part of the U.S. DOE’s Department of Energy’s Integrated Multisector Multiscale Modeling (IM3) basic research project: an Agent Based Model (ABM) for more accurately simulating irrigation demand based on the economics of various crop types [@abm; @Yoone2020431118], and an improved reservoir operations module for more accurately simulating reservoir release based on data-driven harmonic functions [@reservoirs]. Future planned experiments will study the uncertainty characterization and quantification [@hess-19-3239-2015] of the model based on perturbations in key parameters.

# Acknowledgements
This research was supported by the U.S. Department of Energy, Office of Science, as part of research in MultiSector Dynamics, Earth and Environmental System Modeling Program. The Pacific Northwest National Laboratory is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830. The views and opinions expressed in this paper are those of the authors alone.

# References
Validation files will be downloaded and placed in this directory -- they will be ignored by git.
Output files will be placed into subfolders in this directory -- they will be ignored by git.
For convenience, input files can be placed in this directory -- they will be ignored by git.
## create_grand_parameters.py

This utility method generates the four dam/reservoir related input files expected by `mosartwmpy`:
* `grand_reservoir_parameters.nc` - dam/reservoir physical and behavioral parameters
* `grand_average_monthly_flow.parquet` - mean monthly flow across the reservoir during the expected simulation period
* `grand_average_monthly_demand.parquet` - mean monthly demand on the reservoir's water during the expected simulation period
* `grand_dependency_database.parquet` - mapping between GRanD ID and grid cell IDs allowed to extract water

Several datasets are required to perform this operation:
* GRanD reservoir shapefiles, i.e. [GRanD_v1_3](http://globaldamwatch.org/grand/)
* Output data from a previous `mosartwmpy` simulation, with water management disabled, covering the desired simulation time period
* Monthly demand input that will be used for the simulation
* The `mosartwmpy` grid (river network) file
* Elevation data covering the `mosartwmpy` domain in parquet format (this data can be upscaled within the utility method if necessary); for instance from [HYDROSHEDS](https://www.hydrosheds.org/downloads) void-filled elevation
* ISTARF dataset with the data-driven reservoir coefficients, i.e. [ISTARF v0.0.1](https://zenodo.org/record/4602277)

Note that the reservoir parameters and dependency database provided in the tutorial are reasonably robust for a 1/8 degree grid --
for most use cases it would be sufficient to simply update the mean flow and demand files as appropriate to your simulation.

Once the necessary data has been collected, run the utility with the `create_grand_parameters` that is installed along with `mosartwmpy`.
This script will ask for the locations of the datasets and desired output locations.
For more complete control of the script, examine the [method signature](create_grand_parameters.py) and invoke directly from python.[![build](https://github.com/IMMM-SFA/mosartwmpy/actions/workflows/build.yml/badge.svg)](https://github.com/IMMM-SFA/mosartwmpy/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/IMMM-SFA/mosartwmpy/branch/main/graph/badge.svg?token=IPOY8984MB)](https://codecov.io/gh/IMMM-SFA/mosartwmpy)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03221/status.svg)](https://doi.org/10.21105/joss.03221)
[![DOI](https://zenodo.org/badge/312114600.svg)](https://zenodo.org/badge/latestdoi/312114600)

## mosartwmpy

`mosartwmpy` is a python translation of MOSART-WM, a model for water routing and reservoir management written in Fortran. The original code can be found at [IWMM](https://github.com/IMMM-SFA/iwmm) and [E3SM](https://github.com/E3SM-Project/E3SM), in which MOSART is the river routing component of a larger suite of earth-science models. The motivation for rewriting is largely for developer convenience -- running, debugging, and adding new capabilities were becoming increasingly difficult due to the complexity of the codebase and lack of familiarity with Fortran. This version aims to be intuitive, lightweight, and well documented, while still being highly interoperable.
For a quick start, check out the [Jupyter notebook tutorial](https://github.com/IMMM-SFA/mosartwmpy/blob/main/notebooks/tutorial.ipynb)!

## getting started

Ensure you have Python >= 3.7 available (consider using a [virtual environment](https://github.com/pyenv/pyenv), see the docs [here](https://mosartwmpy.readthedocs.io/en/latest/virtualenv.html) for a brief tutorial), then install `mosartwmpy` with:
```shell
pip install mosartwmpy
```

Alternatively, install via conda with:
```shell
conda install -c conda-forge mosartwmpy
```

Download a sample input dataset spanning May 1981 by running the following and selecting option `1` for "tutorial". This will download and unpack the inputs to your current directory. Optionally specify a path to download and extract to instead of the current directory.

```shell
python -m mosartwmpy.download
```

Settings are defined by the merger of the `mosartwmpy/config_defaults.yaml` and a user specified file which can override any of the default settings. Create a `config.yaml` file that defines your simulation (if you chose an alternate download directory in the step above, you will need to update the paths to point at your data):

> `config.yaml`
> ```yaml
> simulation:
>   name: tutorial
>   start_date: 1981-05-24
>   end_date: 1981-05-26
>
> grid:
>   path: ./input/domains/mosart.nc
>   land:
>     path: ./input/domains/land.nc
>
> runoff:
>   read_from_file: true
>   path: ./input/runoff/runoff_1981_05.nc
>
> water_management:
>   enabled: true
>   demand:
>     read_from_file: true
>     path: ./input/demand/demand_1981_05.nc
>   reservoirs:
>     enable_istarf: true
>     parameters:
>       path: ./input/reservoirs/reservoirs.nc
>     dependencies:
>       path: ./input/reservoirs/dependency_database.parquet
>     streamflow:
>       path: ./input/reservoirs/mean_monthly_reservoir_flow.parquet
>     demand:
>       path: ./input/reservoirs/mean_monthly_reservoir_demand.parquet
> ```

`mosartwmpy` implements the [Basic Model Interface](https://csdms.colorado.edu/wiki/BMI) defined by the CSDMS, so driving it should be familiar to those accustomed to the BMI. To launch the simulation, open a python shell and run the following:

```python
from mosartwmpy import Model

# path to the configuration yaml file
config_file = 'config.yaml'

# initialize the model
mosart_wm = Model()
mosart_wm.initialize(config_file)

# advance the model one timestep
mosart_wm.update()

# advance until the `simulation.end_date` specified in config.yaml
mosart_wm.update_until(mosart_wm.get_end_time())
```

## model input

Input for `mosartwmpy` consists of many files defining the characteristics of the discrete grid, the river network, surface and subsurface runoff, water demand, and dams/reservoirs.
Currently, the gridded data is expected to be provided at the same spatial resolution.
Runoff input can be provided at any time resolution; each timestep will select the runoff at the closest time in the past.
Currently, demand input is read monthly but will also pad to the closest time in the past.
Efforts are under way for more robust demand handling.
Dams/reservoirs require four different input files: the physical characteristics, the average monthly flow expected during the simulation period, the average monthly demand expected during the simulation period, and a database mapping each GRanD ID to grid cell IDs allowed to extract water from it.
These dam/reservoir input files can be generated from raw GRanD data, raw elevation data, and raw ISTARF data using the [provided utility](mosartwmpy/utilities/CREATE_GRAND_PARAMETERS.md).
The best way to understand the expected format of the input files is to examine the sample inputs provided by the download utility: `python -m mosartwmpy.download`.

## model output

By default, key model variables are output on a monthly basis at a daily averaged resolution to `./output/<simulation name>/<simulation name>_<year>_<month>.nc`. See the configuration file for examples of how to modify the outputs, and the `./mosartwmpy/state/state.py` file for state variable names.

Alternatively, certain model outputs deemed most important can be accessed using the BMI interface methods. For example:
```python
from mosartwmpy import Model

mosart_wm = Model()
mosart_wm.initialize()

# get a list of model output variables
mosart_wm.get_output_var_names()

# get the flattened numpy.ndarray of values for an output variable
supply = mosart_wm.get_value_ptr('supply_water_amount')
```

## visualization

`Model` instances can plot the current value of certain input and output variables (those available from `Model.get_output_var_name` and `Model.get_input_var_names`):

```python
from mosartwmpy import Model
config_file = 'config.yaml'
mosart_wm = Model()
mosart_wm.initialize(config_file)
for _ in range(8):
    mosart_wm.update()

mosart_wm.plot_variable('outgoing_water_volume_transport_along_river_channel', log_scale=True)
```
![River transport](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/river_transport.png)

Using provided utility functions, the output of a simulation can be plotted as well.

Plot the storage, inflow, and outflow of a particular GRanD dam:
```python
from mosartwmpy import Model
from mosartwmpy.plotting.plot import plot_reservoir
config_file = 'config.yaml'
mosart_wm = Model()
mosart_wm.initialize(config_file)
mosart_wm.update_until()

plot_reservoir(
    model=mosart_wm,
    grand_id=310,
    start='1981-05-01',
    end='1981-05-31',
)
```
![Grand Coulee](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/grand_coulee_1981_05.png)

Plot a particular output variable (as defined in `config.yaml`) over time:
```python
from mosartwmpy import Model
from mosartwmpy.plotting.plot import plot_variable
config_file = 'config.yaml'
mosart_wm = Model()
mosart_wm.initialize(config_file)
mosart_wm.update_until()

plot_variable(
    model=mosart_wm,
    variable='RIVER_DISCHARGE_OVER_LAND_LIQ',
    start='1981-05-01',
    end='1981-05-31',
    log_scale=True,
    cmap='winter_r',
)
```
![River network no tiles](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/river_without_tiles_1981_05.png)

If `cartopy`, `scipy`, and `geoviews` are installed, tiles can be displayed along with the plot:
```python
plot_variable(
    model=mosart_wm,
    variable='RIVER_DISCHARGE_OVER_LAND_LIQ',
    start='1981-05-01',
    end='1981-05-31',
    log_scale=True,
    cmap='winter_r',
    tiles='StamenWatercolor'
)
```
![River network with tiles](https://github.com/IMMM-SFA/mosartwmpy/raw/main/docs/_static/river_with_tiles_1981_05.png)

## model coupling

A common use case for `mosartwmpy` is to run coupled with output from the Community Land Model (CLM). To see an example of how to drive `mosartwmpy` with runoff from a coupled model, check out the [Jupyter notebook tutorial](https://github.com/IMMM-SFA/mosartwmpy/blob/main/notebooks/tutorial.ipynb)!

## testing and validation

Before running the tests or validation, make sure to download the "sample_input" and "validation" datasets using the download utility `python -m mosartwmpy.download`.

To execute the tests, run `./test.sh` or `python -m unittest discover mosartwmpy/tests` from the repository root.

To execute the validation, run a model simulation that includes the years 1981 - 1982, note your output directory, and then run `python -m mosartwmpy.validate` from the repository root. This will ask you for the simulation output directory, think for a moment, and then open a figure with several plots representing the NMAE (Normalized Mean Absolute Error) as a percentage and the spatial sums of several key variables compared between your simulation and the validation scenario. Use these plots to assist you in determining if the changes you have made to the code have caused unintended deviation from the validation scenario. The NMAE should be 0% across time if you have caused no deviations. A non-zero NMAE indicates numerical difference between your simulation and the validation scenario. This might be caused by changes you have made to the code, or alternatively by running a simulation with different configuration or parameters (i.e. larger timestep, fewer iterations, etc). The plots of the spatial sums can assist you in determining what changed and the overall magnitude of the changes.

If you wish to merge code changes that intentionally cause significant deviation from the validation scenario, please work with the maintainers to create a new validation dataset.
Python Virtual Environments
===========================

Maintaining different versions of Python can be a chore. Thankfully, there are many tools for managing Python environments; here are a few recommendations:

* `PyCharm <https://www.jetbrains.com/pycharm/>`_ IDE -- great for developing code in Python, and can automatically create virtual environments for a codebase by detecting versions and dependencies from the ``setup.py`` or ``setup.cfg``.
* `Conda <https://docs.conda.io>`_ package manager -- a Python package manager focused on scientific computing that can also manage virtual environments.
* `pyenv <https://github.com/pyenv/pyenv>`_ CLI -- a shell based tool for installing and switching between different versions of Python and dependencies. I will give a brief tutorial of using ``pyenv`` below, but recognize that the instructions may change over time so the ``pyenv`` documentation is the best place to look.

To create a Python 3.9 virtual environment, try the following steps:

* Install pyenv:

  - if on Mac, use `brew <https://brew.sh/>`_: brew install pyenv
  - if on a linux system, try `pyenv-installer <https://github.com/pyenv/pyenv-installer>`_
  - if on Windows, try `pyenv-win <https://github.com/pyenv-win/pyenv-win>`_

* Install Python 3.9:

  - in a shell, run ``pyenv install 3.9.1``

* Activate Python 3.9 in the current shell

  - in the shell, run ``pyenv shell 3.9.1``

* Proceed with the install of mosartwmpy:

  - in the same shell, run ``pip install mosartwmpy``

* Now you can interact with ``mosartwmpy`` in this current shell session

  - if you start a new shell session you will need to run ``pyenv shell 3.9.1`` again before proceeding
  - this new shell session should maintain all previously pip installed modules for Python 3.9.1mosartwmpy Python API
=====================

mosartwmpy.model module
-----------------------

.. automodule:: mosartwmpy.model
   :members:
   :undoc-members:
   :show-inheritance::notoc:

.. mosartwmpy documentation master file, created by

.. module:: mosartwmpy

************************
mosartwmpy documentation
************************

**Date**: |today| **Version**: |version|

**Useful links**:
`Source Repository <https://github.com/immm-sfa/mosartwmpy>`__ |
`Issues & Ideas <https://github.com/immm-sfa/mosartwmpy/issues>`__

`mosartwmpy` is a Python translation of MOSART-WM, a water routing and reservoir management model written in Fortran.


.. panels::
    :card: + intro-card text-center
    :column: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex

    >>>
    :img-top: _static/cognitive.svg

    Getting started
    ^^^

    Get to know the *mosartwmpy* model.

    +++
    .. link-button:: README
            :type: ref
            :text: Getting started
            :classes: btn-block btn-secondary stretched-link

    >>>
    :img-top: _static/education.svg

    Tutorial
    ^^^

    Follow along with this Jupyter notebook to learn the ropes of *mosartwmpy*.

    +++
    .. link-button:: tutorial
            :type: ref
            :text: Tutorial
            :classes: btn-block btn-secondary stretched-link

    >>>
    :img-top: _static/soccer.svg

    Tips & tricks
    ^^^

    Learn about ways to manage Python virtual environments.

    +++
    .. link-button:: virtualenv
            :type: ref
            :text: Virtual environments
            :classes: btn-block btn-secondary stretched-link

    >>>
    :img-top: _static/api.svg

    API reference
    ^^^

    A detailed description of the *mosartwmpy* API.

    +++
    .. link-button:: mosartwmpy
            :type: ref
            :text: API
            :classes: btn-block btn-secondary stretched-link


.. toctree::
    :maxdepth: 1
    :hidden:
    :titlesonly:

    README.md
    tutorial.ipynb
    virtualenv
    mosartwmpy
