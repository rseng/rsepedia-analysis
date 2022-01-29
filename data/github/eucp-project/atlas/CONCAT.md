# Atlas of constrained climate projections- EUCP WP2

[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=eucp-project_atlas&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=eucp-project_atlas)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5654741.svg)](https://doi.org/10.5281/zenodo.5654741)
[![Research Software Directory Badge](https://img.shields.io/badge/rsd-eucp_atlas-00a3e3.svg)](https://www.research-software.nl/software/eucp-atlas)

This repository contains the source code and content for an [Atlas of
constrained climate projections](https://eucp-project.github.io/atlas/). The
atlas demonstrates outputs from the probabilistic projection methods developed
or assessed in the [European Climate Projection
system](https://www.eucp-project.eu/) (EUCP) Horizon2020 project. For more info,
see the [Atlas about page](https://eucp-project.github.io/atlas/about).

## Citation

To cite this repository, use the information available at [CITATION.cff](CITATION.cff),
and to cite the content of the Atlas, see the [Atlas about page](https://eucp-project.github.io/atlas/about).

## License

The source code is licensed under [Apache 2.0](./LICENSE), whereas the content
e.g. maps in the `static` directory of this repository, are licensed under
[CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Maintainers

Current maintainers of the Atlas are Research software engineers from the
[Netherlands eScience Center](https://www.esciencecenter.nl/). If you have any
questions or concerns, please submit an
[issue](https://github.com/eucp-project/atlas/issues). Maintainers will do their
best to help you.

## Contributions

For information on how to contribute to this Atlas, please check the
[contributing guidelines](CONTRIBUTING.md).

## Acknowledgements

FIXME
# Contributing guidelines

## In general

Contributions are very welcome. Please make sure there is a GitHub issue
associated with with every pull request. [Creating an issue](https://github.com/eucp-project/atlas/issues/new) is also a good way to
propose new features.

## Add maps to the Atlas

The maps shown in Atlas are created using the data available at [Zenodo](https://zenodo.org/record/5645154),
and are stored in `.png` format in the [assets](./assets) directory of this
repository. You can use the scripts and notebooks as examples to develop a new
script for your own data.

### Processing data and preparing plots

Content for the Atlas page is generated in two steps:

1. Preprocess data to clean NetCDF files
2. Create maps in a uniform way based on the preprocessed data files

The scripts that perform these actions are stored in the [python](./python)
directory. For more instructions, see [README](./python/README.md).

When adding a new script or notebook, please place a short description in the
[README](./python/README.md) of this directory.

## Edit the Atlas pages

The atlas is created with
[Nuxt.js](https://nuxtjs.org/docs/get-started/installation).

### Create a local build

To locally render the Atlas, run the following:

```bash
# install dependencies
$ npm install

# serve with hot reload at localhost:3000
$ npm run dev

# build for production and launch server
$ npm run build
$ npm run start

# generate static project
$ npm run generate
```

For detailed explanations on how things work, check out [Nuxt.js
docs](https://nuxtjs.org).

### Edit the pages

The [pages](./pages) directory contains your Application Views and Routes. The
framework reads all the `*.vue` files inside this directory and creates the
router of your application.

More information about the usage of this directory in [the
documentation](https://nuxtjs.org/guide/routing).

## Making a release

### Author information

Ensure all authors are present in:

- `CITATION.cff`

### Confirm release info

Ensure the right date and upcoming version number is set in:

- `CITATION.cff`
- `package.json`

### Release on GitHub

Open [releases](https://github.com/eucp-project/atlas/releases) and draft a new
release.

Tag the release according to semantic versioning guidelines, preceded with a `v`
(e.g.: v1.0.0). The release title is the tag and the release date together
(e.g.: v1.0.0 (2019-07-25)). Tick the pre-release box in case the release is a
candidate release, and amend the version tag with `rc` and the candidate number.

### Release on Zenodo

Confirm the new release on [Zenodo](https://zenodo.org/record/5654741).

### Release on the Research Software Directory

Wait a few hours, then confirm the addition of a new release on the
[RSD](https://www.research-software.nl/software/eucp-atlas).
# PLUGINS

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains Javascript plugins that you want to run before mounting the root Vue.js application.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/plugins).
---
question: Which method is best?
sort: 4
---
The question of which constraint method is best depends on what the projection
is being used for. Some methods have particular strengths over others. For
example, some methods provide relative weights of available simulations that
could be used to select models to drive downstream climate impact analysis, and
some methods can provide physically consistent projections for a number of
variables.
![summary of skills](summary_available_info.png)

Method performance can be different in different areas or with different
variables. One method, for instance, could bring significant improvements over
unconstrained projections in northern Europe but not in the Mediterranean
region.
![summary of skills](summary_skill.png)

All methods show improvements in some regions for temperature projections, whereas there is less evidence that they improve precipitation projections.  One
method, CALL, produces projections for the Mediterranean that are markedly worse
than just using the raw climate model projections.
---
question: Can I use these figures?
sort: 6
---

Yes. The images are licences under a [CC-BY
licence](https://creativecommons.org/licenses/by/4.0/), which means you are free
to use them, but please acknowledge the source.

FIXME: add citation info for the atlas itself once published.
---
question: What do these maps show?
sort: 1
---

This Atlas shows several different model projections of the European climate in
2050, illustrating the effect of constraining the projections. Each map shows
differences between 2041-2060 mean conditions with respect to the 1995-2014
baseline. Not all climate simulations are equally likely, however
this is not routinely addressed in current climate projections. Methods of
isolating the more likely projections include comparing a simulation of the past
to observed climate data recorded at the time. This can suggest how well the
simulation may project the future. This process is known as 'constraint' and can
affect the results given by an ensemble of climate simulations by ruling out the
poorer simulations.
---
question: I'd like to cite the Atlas in my paper, how can I do this?
sort: 7
---

Great! You can cite the Atlas using its DOI: 10.5281/zenodo.5654741, and
citation information on the [Zenodo
repository](https://doi.org/10.5281/zenodo.5654741).

If you want to cite a particular figure or piece of content shown in this Atlas,
the information for its citation can be found under "Where can I find more
information?".---
question: Which maps are available?
sort: 9
---

![available combinations](placeholder.png)

_The ASK unconstrained maps:_
  
  The ASK methodology does not use the information on
  a prior or unconstrained range, instead using the information from a model-based
  estimate of the human fingerprint in the observed climate, to scale the
  multi-model mean future climate response.
---
question: What do the percentiles mean?
sort: 2
---

A percentile is the level below which a given percentage of simulations fall
when they are put in order according to a variable such as temperature. For
instance, the 50th temperature percentile shows the highest temperature in the
lower 50% of the distribution of temperature changes. The 10th percentile of
precipitation projections may show significant drying, while the 90th percentile
might show much wetter conditions.
---
question: What are the differences between the methods?
sort: 3
---

There are several differences between the constraint methods used in this analysis, including the kind of observational data they use and how the results are handled statistically. For instance, the UKCP and ClimWIP methods both use multiple observed variables (e.g. temperature, precipitation) to constrain their projections, rather than just using the target variable being reported.

Some methods constrain their variables at different spatial scales, from global to local, before deriving constrained local projections. Different methods also handle uncertainty in the projections differently.

The following table includes the key characteristics of the different methods and is reprinted from Brunner et al. (2020), which contains more information on the differences between these methods.

![methods table](methods.png)
---
question: Is the data available as well?
sort: 8
---

Data the underpins the Climate Projection Atlas is being made available for
public download via Zenodo:
[doi:10.5281/zenodo.5645153](https://doi.org/10.5281/zenodo.5645153).

If you use the data in your work, please cite it following the citation
information provided on the Zenodo page.
---
question: What if I want to apply these methods myself?
sort: 10
---

For three of the methods, software implementation are available:

* ClimWIP
    * [ESMValTool implementation](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html)
    * [ClimWIP development version](https://github.com/lukasbrunner/ClimWIP/)
    * [ClimWIP version used for Brunner 2020 paper](https://doi.org/10.5281/zenodo.4073039)
* KCC
    * [R package](https://doi.org/10.5281/zenodo.5233947)
    * [Shiny app](https://saidqasmi.shinyapps.io/KCC-shinyapp/)
    * [Jupyter Notebook](https://gitlab.com/saidqasmi/kcc_notebook)
* [REA software implementation](http://doi.org/10.5281/zenodo.3890966)

For the other methods, please refer to their reference papers listed under
"more information" below.
---
question: Where can I find more information?
sort: 11
---

Please view the references section below for more on the background research
behind this Atlas. For more on the EUCP project, please head to [our
website](https://www.eucp-project.eu/). This Atlas has been produced by Work
Package 2 of EUCP, which is described in more detail
[here](https://www.eucp-project.eu/the-project/work-packages-wps/wp-2/).

The following paper provides an overview including all methods:

* [Comparing Methods to Constrain Future European Climate Projections Using a
  Consistent Framework](https://doi.org/10.1175/JCLI-D-19-0953.1)

For the individual methods, the following references provide more information:

* REA (CNRM/CNRS):
    * [Making climate projections conditional on historical observation](https://doi.org/10.1126/sciadv.abc0671)
    * [Reducing uncertainty in local climate projection](https://doi.org/10.21203/rs.3.rs-364943/v1)
* ClimWIP (ETHZ):
    * [Quantifying uncertainty in European climate projections using combined performance-independence weighting](https://doi.org/10.1088/1748-9326/ab492f)
    * [Reduced global warming from CMIP6 projections when weighting models by performance and independence. Earth System Dynamics](https://doi.org/10.5194/esd-11-995-2020)
    * [An investigation of weighting schemes suitable for incorporating large ensembles into multi-model ensembles](https://doi.org/10.5194/esd-11-807-2020)
    * [Partitioning climate projection uncertainty with multiple large ensembles and CMIP5/6](https://doi.org/10.5194/esd-11-491-2020)
* CALL (UOxf):
    * [Calibrating large-ensemble European climate projections using observational data](https://doi.org/10.5194/esd-11-1033-2020)
* ASK (UEdin):
    * [Towards consistent observational constraints in climate predictions and projection](doi:10.3389/fclim.2021.678109)
* UKCP (UKMO):
    * [UKCP18 Land Projections: Science Report](https://www.metoffice.gov.uk/pub/data/weather/uk/ukcp18/science-reports/UKCP18-Land-report.pdf)
    * [UKCP18 Science Overview Report](https://www.metoffice.gov.uk/pub/data/weather/uk/ukcp18/science-reports/UKCP18-Overview-report.pdf)
---
question: Who are these projections for?
sort: 5
---

Climate projections can be useful for many people, from communities to
policymakers to businesses. We have prepared some example use cases to
illustrate the utility of the projections in this Atlas.
See the
<NuxtLink :to="`/examples`" class="hover:text-blue-400 underline">EXAMPLES</NuxtLink>.
# Processing data with Python (notebooks)

In this folder, all scripts and notebooks are stored that are used to pre-process
data and to generate the maps in the Atlas.

## Requirements

The dependencies of the notebooks and scripts can be installed in a Conda environment with

```shell
# From this directory
conda install mamba -n base -c conda-forge -y
mamba env create --file environment.yml
conda activate atlas
```

## Preprocess NetCDF data

We provide some notebooks that check the original/raw data, fix/add the metadata
using
[CF-conventions](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.9/cf-conventions.html)
and save data in a NetCDF format. As the output of a method (i.e.
original/raw data) is provided by a specific institute, there is one notebook
per each `institute`-`method`:

- [Preprocess CNRM-KCC data](cleanup_CNRM_KCC_atlas_netcdf.ipynb)
- [Preprocess EdinU-ASK data](cleanup_EdinU_ASK_atlas_netcdf.ipynb)
- [Preprocess ETHZ-ClimWIP data](cleanup_ETHZ_ClimWIP_atlas_netcdf.ipynb)
- [Preprocess ICTP-REA data](cleanup_ICTP_REA_atlas_netcdf.ipynb)
- [Preprocess UKMO-UKCP data](cleanup_UKMO_UKCP_atlas_netcdf.ipynb)
- [Preprocess UOxf-CALL data](cleanup_UOxf_CALL_atlas_netcdf.ipynb)

To run a notebook, you only need to specify the path to raw data as `datapath`
and a path to store the output as `output_path`. Defaults are:

- `datapath = "./AtlasData/raw"`

- `output_path = "./AtlasData/preprocess"`

The pre-processed data follows the following standards:

### coordinates

- climatology_bounds (climatology_bounds) datetime64[ns] ['2050-06-01', '2050-09-01', '2050-12-01', '2051-03-01']
- time (time) (datetime64[ns]) [2050-07-16 2051-01-16] # "JJA", "DJF"
- latitude (lat) (float64) [30, ..., 75]
- longitude (lon) (float64) [-10, ..., 40]
- percentile (percentile) (int64) [10, 25, 50, 75, 90]

### variables

- tas (time, latitude, longitude, percentile) (float64)
- pr (time, latitude, longitude, percentile) (float64)

### attributes

The attributes of variables and coordinates are defined as:

- "tas": {
    "description": "Change in Air Temperature",
    "standard_name": "Change in Air Temperature",
    "long_name": "Change in Near-Surface Air Temperature",
    "units": "K", 
    "cell_methods": "time: mean changes over 20 years 2041-2060 vs 1995-2014",
},
- "pr": {
    "description": "Relative precipitation",
    "standard_name": "Relative precipitation",
    "long_name": "Relative precipitation",
    "units": "%",  
    "cell_methods": "time: mean changes over 20 years 2041-2060 vs 1995-2014",
},
- "latitude": {"units": "degrees_north", "long_name": "latitude", "axis": "Y"},
- "longitude": {"units": "degrees_east", "long_name": "longitude", "axis": "X"},
- "time": {
    "climatology": "climatology_bounds",
    "long_name": "time",
    "axis": "T",
    "climatology_bounds": ["2050-6-1", "2050-9-1", "2050-12-1", "2051-3-1"],
    "description": "mean changes over 20 years 2041-2060 vs 1995-2014. The mid point 2050 is chosen as the representative time.",
},
- "percentile": {"units": "%", "long_name": "percentile", "axis": "Z"},

The attributes of the data is defined as:

- "description": "Contains modified `institute` `method` data used for Atlas in EUCP project.",
- "history": "original `institute` `method` data files ...",

### output file names

output_file_name = `prefix_activity_institution-id_source_method_sub-method_cmor-var`

> example: atlas_EUCP_CNRM_CMIP6_KCC_cons_tas.nc


## Create maps

> make sure that the conda environment `atlas` is activated.

Maps are created using pre-processed data and `maps_creator_atlas_data.py` script:

```shell
python ./atlas/python/maps_creator_atlas_data.py --inputdir "./AtlasData/preprocess" --outputdir "./atlas/assets/processed_figures"
```

## Add new maps

If you want to add new maps to the Atlas, please publish the pre-processed data on
Zenodo and add references to [FIXME].

## Additional notebooks

These notebooks can be used to tweak plot settings and preview maps using raw
model data.

- [Preview map of relative precipitation](maps_prototype_prec.ipynb)
- [Preview map of temperature](maps_prototype_tas.ipynb)
# MIDDLEWARE

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains your application middleware.
Middleware let you define custom functions that can be run before rendering either a page or a group of pages.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/routing#middleware).
# STATIC

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains your static files.
Each file inside this directory is mapped to `/`.
Thus you'd want to delete this README.md before deploying to production.

Example: `/static/robots.txt` is mapped as `/robots.txt`.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/assets#static).
# STORE

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains your Vuex Store files.
Vuex Store option is implemented in the Nuxt.js framework.

Creating a file in this directory automatically activates the option in the framework.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/vuex-store).
# LAYOUTS

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains your Application Layouts.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/views#layouts).
