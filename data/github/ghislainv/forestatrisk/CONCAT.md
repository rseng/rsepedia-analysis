---
title: "forestatrisk: a Python package for modelling and forecasting deforestation in the tropics"
tags:
  - Python
  - land use change
  - spatial modelling
  - spatial analysis
  - forecasting
  - spatial autocorrelation
  - tropical forests
  - roads
  - protected areas
  - biodiversity scenario
  - ipbes
  - co2 emissions
  - ipcc
authors:
  - name: Ghislain Vieilledent
    orcid: 0000-0002-1685-4997
    affiliation: "1, 2, 3, 4"
affiliations:
  - name: CIRAD, UMR AMAP, F--34398 Montpellier, France
    index: 1
  - name: CIRAD, Forêts et Sociétés, F--34398 Montpellier, France.
    index: 2
  - name: AMAP, Univ Montpellier, CIRAD, CNRS, INRAE, IRD, Montpellier, France.
    index: 3
  - name: European Commission, Joint Research Centre (JRC), I--21027 Ispra (VA), Italy.
    index: 4
date: 6 December 2020
output:
  bookdown::pdf_document2:
    keep_md: yes
    keep_tex: yes
    citation_package: "natbib"
bibliography: paper.bib
#bibliography: /home/ghislain/Documents/Bibliography/biblio.bib
link-citations: yes
---

# Summary

The `forestatrisk` Python package can be used to model the spatial probability of deforestation and predict future forest cover in the tropics. The spatial data used to model deforestation come from georeferenced raster files, which can be very large (several gigabytes). The functions available in the `forestatrisk` package process large rasters by blocks of data, making calculations fast and efficient. This allows deforestation to be modeled over large geographic areas (e.g., at the scale of a country) and at high spatial resolution (e.g., $\leq$ 30 m). The `forestatrisk` package offers the possibility of using logistic regression with auto-correlated spatial random effects to model the deforestation process. The spatial random effects make possible to structure the residual spatial variability of the deforestation process, not explained by the variables of the model and often very large. In addition to these new features, the `forestatrisk` Python package is open source (GPLv3 license), cross-platform, scriptable (via Python), user-friendly (functions provided with full documentation and examples), and easily extendable (with additional statistical models for example). The `forestatrisk` Python package has been used to model deforestation and predict future forest cover by 2100 across the humid tropics.

<!-- The JOSS paper should only include:

- A list of the authors of the software and their affiliations.
- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.
- A clear statement of need that illustrates the purpose of the software.
- A description of how this software compares to other commonly-used packages in this research area.
- Mentions (if applicable) of any ongoing research projects using the software or recent scholarly publications enabled by it.
- A list of key references including a link to the software archive.

Compile with the following command:
docker run --rm \
    --volume $PWD/paper:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/paperdraft

-->

# Statement of Need

Commonly called the "Jewels of the Earth", tropical forests shelter 30 million species of plants and animals representing half of the Earth's wildlife and at least two-thirds of its plant species [@Gibson2011]. Through photosynthesis and carbon sequestration, tropical forests play an important role in the global carbon cycle, and in regulating the global climate [@Baccini2017]. Despite the many ecosystem services they provide, tropical forests are disappearing at an alarming rate [@Keenan2015; @Vancutsem2020], mostly because of human activities. Currently, around 8 Mha (twice the size of Switzerland) of tropical forest are disappearing each year [@Keenan2015]. Spatial modelling of deforestation allows identifying the main factors that determine the spatial risk of deforestation and quantifying their relative effects. Forecasting forest cover change is paramount as it allows anticipating the consequences of deforestation (in terms of carbon emissions or biodiversity loss) under various technological, political, and socioeconomic scenarios, and informs decision makers accordingly [@Clark2001]. Because both biodiversity and carbon vary greatly in space [@Allnutt2008; @Baccini2017], it is necessary to provide spatial forecasts of forest cover change to properly quantify biodiversity loss and carbon emissions associated with future deforestation.

The `forestatrisk` Python package can be used to model tropical deforestation spatially, predict the spatial risk of deforestation, and forecast future forest cover in the tropics (\autoref{fig:prob}). Several other software tools allow modeling and forecasting of forest cover change [@Mas2014]. The most famous land cover change software tools include [Dinamica-EGO](https://csr.ufmg.br/dinamica/) [@Soares-Filho2002], [Land Change Modeller](https://clarklabs.org/terrset/land-change-modeler/) [@Eastman2017], and [CLUE](http://www.ivm.vu.nl/en/Organisation/departments/spatial-analysis-decision-support/Clue/) [@Verburg2009]. Despite the many functionalities they provide, these software tools are not open source and might not all be cross-platform, scriptable, and completely user-friendly. Moreover, the statistical approaches they propose to model land cover change do not take into account the residual spatial variability in the deforestation process that is not explained by the model's variables, and which is often very large. The more sophisticated algorithms they use (genetic algorithms, artificial neural networks, or machine learning algorithms) might also have the tendency to over-fit the data [@Mas2014]. Finally, application of these software tools to large spatial scales (e.g., at the country or continental scale) with high resolution data (e.g., $\leq$ 30 m) has not yet been demonstrated [but see @Soares-Filho2006 for a study in the Amazon at 1 km resolution]. The `forestatrisk` Python package aims to fill some of these gaps and to enlarge the range of software available to model and forecast tropical deforestation.

# Main functionalities

## A set of functions for modelling and forecasting deforestation

The `forestatrisk` Python package includes functions to (i) compute the forest cover change raster and the rasters of explanatory variables for a given country from several global datasets (such as OpenStreetMap or the SRTM Digital Elevation Database v4.1 for example) (ii) efficiently sample forest cover change observations and retrieve information on spatial explanatory variables for each observation, (iii) estimate the parameters of various statistical deforestation models, (iv) predict the spatial probability of deforestation, (v) forecast the likely forest cover in the future, (vi) validate the models and the projected maps of forest cover change, (vii) estimate carbon emissions associated with future deforestation, and (viii) plot the results. The `forestatrisk` package includes a hierarchical Bayesian logistic regression model with autocorrelated spatial random effects, which is well suited for modeling deforestation (see below). Any statistical model class with a `.predict()` method can potentially be used together with the function `forestatrisk.predict_raster()` to predict the spatial risk of deforestation. This allows a wide variety of additional statistical models from other Python packages to be used, such as `scikit-learn` [@Pedregosa2011] for example, for deforestation modelling and forecasting.

## Ability to process large raster data

Spatially-distributed forest cover change and explanatory variables are commonly available as georeferenced raster data. Raster data consist of rows and columns of cells (or pixels), with each cell storing a single value. The resolution of the raster dataset is its pixel width in ground units. Depending on the number of pixels (which is a function of the raster's geographical extent and resolution), raster files might occupy a space of several gigabytes on disk. Processing such large rasters in memory can be prohibitively intensive. Functions in the `forestatrisk` package process large rasters by blocks of pixels representing subsets of the raster data. This makes computation efficient, with low memory usage. Reading and writing subsets of raster data is done by using two methods from the GDAL Python bindings [@GDAL2020]: `gdal.Dataset.ReadAsArray()` and `gdal.Band.WriteArray()`. Numerical computations on arrays are performed with the `NumPy` Python package, whose core is mostly made of optimized and compiled C code that runs quickly [@Harris2020]. This allows the `forestatrisk` Python package to model and forecast forest cover change on large spatial scales (e.g., at the country or continental scale) using high resolution data (e.g., $\leq$ 30 m), even on personal computers with average performance hardware. For example, the `forestatrisk` Python package has been used on a personal computer to model and forecast the forest cover change at 30-m resolution for the Democratic Republic of the Congo (total area of 2,345 million km$^2$), processing large raster files of 71,205 $\times$ 70,280 cells without issues.

## Statistical model with autocorrelated spatial random effects

The `forestatrisk` Python package includes a function called `.model_binomial_iCAR()` to estimate the parameters of a logistic regression model including auto-correlated spatial random effects. The model considers the random variable $y_i$ which takes value 1 if a forest pixel $i$ is deforested in a given period of time, and 0 if it is not. The model assumes that $y_i$ follows a Bernoulli distribution of parameter $\theta_i$ (\autoref{eq:icar}). $\theta_i$ represents the spatial relative probability of deforestation for pixel $i$ and is linked, through a logit function, to a linear combination of the explanatory variables $X_i \beta$, where $X_i$ is the vector of explanatory variables for pixel $i$, and $\beta$ is the vector of effects $[\beta_1, \ldots, \beta_n]$ associated with the $n$ variables. The model can include (or not) an intercept $\alpha$. To account for the residual spatial variation in the deforestation process, the model includes additional random effects $\rho_{j(i)}$ for the cells of a spatial grid covering the study-area. The spatial grid resolution has to be chosen in order to have a reasonable balance between a good representation of the spatial variability and a limited number of parameters to estimate. Each observation $i$ is associated with one spatial cell $j(i)$. Random effects $\rho_j$ are assumed to be spatially autocorrelated through an intrinsic conditional autoregressive (iCAR) model [@Besag1991]. In an iCAR model, the random effect $\rho_j$ associated with cell $j$ depends on the values of the random effects $\rho_{j^{\prime}}$ associated with neighboring cells $j^{\prime}$. The variance of the spatial random effects $\rho_j$ is denoted by $V_{\rho}$. The number of neighbouring cells for cell $j$ (which might vary) is denoted by $n_j$. Spatial random effects $\rho_j$ account for unmeasured or unmeasurable variables [@Clark2005], which explain a part of the residual spatial variation in the deforestation process that is not explained by the fixed (i.e., explanatory) variables ($X_i$). The parameter inference is done in a hierarchical Bayesian framework. The `.model_binomial_iCAR()` function calls an adaptive Metropolis-within-Gibbs algorithm [@Rosenthal2011] written in C for maximum computation speed.

\begin{equation}
\begin{split}
  y_i \sim \mathcal{B}ernoulli(\theta_i)\\
  \text{logit}(\theta_i) = \alpha + X_i \beta + \rho_{j(i)}\\
  \rho_{j(i)} \sim \mathcal{N}ormal(\sum_{j^{\prime}} \rho_{j^{\prime}} / n_j,V_{\rho} / n_j)
\end{split}
\label{eq:icar}
\end{equation}

# Applications and perspectives

The Python package `forestatrisk` was recently used to model the spatial probability of deforestation and predict forest cover change by 2100 across the humid tropics (<https://forestatrisk.cirad.fr>). Future developments of the package will focus on expanding documentation, case studies, statistical models, and validation tools. We are convinced that the `forestatrisk` package could be of great help in obtaining estimates of carbon emissions and biodiversity loss under various scenarios of deforestation in the tropics. Such scenarios should help decision-makers take initiatives to tackle climate change and the biodiversity crisis. The results from the `forestatrisk` package could contribute to future IPCC and IPBES reports [@IPCC2014; @IPBES2019], or help implement the REDD$+$ mechanism of the [Paris Agreement](https://unfccc.int/process-and-meetings/the-paris-agreement/the-paris-agreement).

# Figures

![**Map of the spatial probability of deforestation in the Guadeloupe archipelago.** This map was produced with the `forestatrisk` Python package. Colored pixels represent the extent of the natural old-growth tropical moist forest in 2020. The original map has a 30-m resolution. A relative spatial probability of deforestation was computed for each forest pixel. Probability of deforestation is a function of several explanatory variables describing: topography (altitude and slope), accessibility (distances to nearest road, town, and river), forest landscape (distance to forest edge), deforestation history (distance to past deforestation), and land conservation status (presence of a protected area). This map can be reproduced with the /Get started/ tutorial at <https://ecology.ghislainv.fr/forestatrisk>.\label{fig:prob}](prob_joss.png)

# Acknowledgements

I am grateful to Clovis Grinand, Romuald Vaudry, and Matthieu Tiberghien who gave me the opportunity to work on deforestation modeling when we were leading forest conservation projects in Madagascar. I also warmly thank Frédéric Achard and all the members of the IFORCE group for their invaluable support during the first phase of development of the package, while I was seconded to the JRC in Ispra. I would also like to thank Chris Garrard for writing the book _"Geoprocessing with Python"_ [@Garrard2016], which has been of great help in the development of the `forestatrisk` package. This work benefited from funding by FRB-FFEM (BioSceneMada project, AAP-SCEN-2013 I), the European Commission (Roadless Forest project), and CNRT (RELIQUES project).

# References
---
title: ""
output:
  bookdown::html_document2:
    highlight: tango
    number_sections: no
    toc: yes
    toc_float: no
---

```{r setup, echo=FALSE, message=FALSE}
library(knitr)
opts_chunk$set(
  fig.align="center",
  echo=TRUE,
  message=FALSE,
  warning=FALSE,
  cache=TRUE,
  cache.lazy=FALSE 
)
```

### Prerequisites {-}

#### Download the Rmd file

The `.Rmd` file for this R notebook is available [here](https://github.com/ghislainv/forestatrisk/blob/master/docsrc/notebooks_R/get_started_R.Rmd) for download.

#### Install Python and the `reticulate` R package {-}

To use the `forestatrisk` Python package within R, you first need to install Python for your system. The recommended way is to install the latest version of `miniconda3` (which includes Python) following the recommendations [here](https://docs.conda.io/en/latest/miniconda.html).

Second, you need to install the [`reticulate`](https://rstudio.github.io/reticulate/) R package which provides a comprehensive set of tools for interoperability between Python and R.

```{r R_libraries, cache=FALSE}
# Set Miniconda path if necessary
# conda_path <- "/home/ghislain/.pyenv/versions/miniconda3-latest"
# Sys.setenv(RETICULATE_MINICONDA_PATH=conda_path)

# Install library
if (!require("reticulate", character.only=TRUE)) {
  install.packages("reticulate")
  require("reticulate", character.only=TRUE)
}
```

#### Create a Python virtual environment {-}

Using functions from the `reticulate` R package, install the Python virtual environment with a selection of Python packages. Use Python 3.7 and the `conda-forge` channel to install the packages.

```{r condaenv, cache=FALSE}
if (!("conda-far" %in% conda_list()$name)) {
	# Vector of packages
	conda_pkg <- c("python=3.7", "gdal", "numpy", "matplotlib",
				   "pandas", "patsy", "pip", "statsmodels",
				   "earthengine-api")
	# Create conda virtual environment and install packages
	conda_create("conda-far", packages=conda_pkg, forge=TRUE)
	# Install packages not available with conda using pip
	pip_pkg <- c("pywdpa", "sklearn")
	conda_install("conda-far", pip_pkg, pip=TRUE,
				  pip_ignore_installed=FALSE)
	# Install the forestatrisk package
	conda_install("conda-far", "forestatrisk", pip=TRUE,
				  pip_ignore_installed=FALSE)
}
```

Note: for a manual installation of the virtual environment, execute the following command lines in a miniconda terminal.

```
conda create --name conda-far -c conda-forge python=3.7 gdal numpy matplotlib pandas patsy pip statsmodels earthengine-api --yes
conda activate conda-far
pip install pywdpa sklearn
pip install forestatrisk
```

#### Use the Python virtual environment in R {-}

Specify the Python version to use with `reticulate` and check that it is the right one. Function `py_config` will inform you on the Python and Numpy versions detected. If everything is fine, you should obtain an output similar to the following one indicating that the `conda-far` virtual environment is being used.

```{r python_env, cache=FALSE}
# Specify the Python environment to use
use_condaenv("conda-far", required=TRUE)

# Check python configuration
py_config()
```

You are now ready to import the `forestatrisk` Python module in R.

```{r importpy, cache=FALSE}
far <- import("forestatrisk")
```

We will follow the steps of the `Get started` tutorial in Python [here](https://ecology.ghislainv.fr/forestatrisk/notebooks/get_started.html#) and adapt the commands to use the functions from the `forestatrisk` Python package in R. See more in details the [documentation](https://rstudio.github.io/reticulate/) of the `reticulate` R package to see how to adapt Python commands and objects to R.

We first create a directory to hold the outputs.

```{r output-dir}
# Make output directory
dir.create("output")
```

### 1. Data

#### 1.1 Import and unzip the data

We use the [Guadeloupe](https://en.wikipedia.org/wiki/Guadeloupe) archipelago as a case study.

```{r import-data}
# Source of the data
url <- "https://github.com/ghislainv/forestatrisk/raw/master/docsrc/notebooks/data_GLP.zip"

if (!file.exists("data_GLP.zip")) {
  download.file(url, "data_GLP.zip")
  unzip("data_GLP.zip", exdir="data")
}
```

#### 1.2 Files

The `data` folder includes, among other files:

- The forest cover change data for the period 2010--2020 as a GeoTiff raster file (`data/fcc23.tif`).
- Spatial variables as GeoTiff raster files (`.tif` extension, eg. `data/dist_edge.tif `for distance to forest edge).

#### 1.3 Sampling the observations

```{r sampling}
# Sample points
dataset <- far$sample(nsamp=10000L, adapt=TRUE, seed=1234L, csize=10L,
                     var_dir="data",
                     input_forest_raster="fcc23.tif",
                     output_file="output/sample.txt",
                     blk_rows=0L)
```

```{r dataset}
# Remove NA from data-set (otherwise scale() and
# model_binomial_iCAR do not work)
dataset <- dataset[complete.cases(dataset), ]
# Set number of trials to one for far.model_binomial_iCAR()
dataset$trial <- 1
# Print the first five rows
head(dataset)
```

### 2. Model

#### 2.1 Model preparation

```{r model-prep}
# Neighborhood for spatial-autocorrelation
neighborhood <- far$cellneigh(raster="data/fcc23.tif", csize=10L, rank=1L)
nneigh <- neighborhood[[1]]
adj <- neighborhood[[2]]

# List of variables
variables <- c("scale(altitude)", "scale(slope)",
             "scale(dist_defor)", "scale(dist_edge)", "scale(dist_road)",
             "scale(dist_town)", "scale(dist_river)")

# Formula
right_part <- paste0(" + ", paste(variables, collapse=" + "), " + cell")
left_part <- "I(1-fcc23) + trial ~ "
formula <- paste0(left_part, right_part)

# Starting values
beta_start <- -99  # Simple GLM estimates

# Priors
priorVrho <- -1  # -1="1/Gamma"
```

#### 2.2 iCAR model

```{r icar-model}
# Run the model
mod_binomial_iCAR <- far$model_binomial_iCAR(
    # Observations
    suitability_formula=formula, data=dataset,
    # Spatial structure
    n_neighbors=np_array(nneigh,dtype="int32"),
    neighbors=np_array(adj,dtype="int32"),
    # Environment
    eval_env=-1L,
    # Priors
    priorVrho=priorVrho,
    # Chains
    burnin=1000L, mcmc=1000L, thin=1L,
    # Starting values
    beta_start=beta_start)
```

#### 2.3 Model summary

```{r mod-summary}
# Predictions
pred_icar <- mod_binomial_iCAR$theta_pred

# Summary
print(mod_binomial_iCAR)
```

```{r save-summary, results='hide'}
# Write summary in file
sink("output/summary_icar.txt")
mod_binomial_iCAR
sink()
```


### 3 Predict

#### 3.1 Interpolate spatial random effects

```{r interpolate}
# Spatial random effects
rho <- mod_binomial_iCAR$rho

# Interpolate
far$interpolate_rho(rho=rho, input_raster="data/fcc23.tif",
                    output_file="output/rho.tif",
                    csize_orig=10L, csize_new=1L)
```

#### 3.2 Predict deforestation probability

```{r predict, results='hide'}
# Update dist_edge and dist_defor at t3
file.rename("data/dist_edge.tif", "data/dist_edge.tif.bak")
file.rename("data/dist_defor.tif", "data/dist_defor.tif.bak")
file.copy("data/forecast/dist_edge_forecast.tif", "data/dist_edge.tif")
file.copy("data/forecast/dist_defor_forecast.tif", "data/dist_defor.tif")

# Compute predictions
far$predict_raster_binomial_iCAR(
    mod_binomial_iCAR, var_dir="data",
    input_cell_raster="output/rho.tif",
    input_forest_raster="data/forest/forest_t3.tif",
    output_file="output/prob.tif",
    blk_rows=10L  # Reduced number of lines to avoid memory problems
)

# Reinitialize data
file.remove("data/dist_edge.tif")
file.remove("data/dist_defor.tif")
file.rename("data/dist_edge.tif.bak", "data/dist_edge.tif")
file.rename("data/dist_defor.tif.bak", "data/dist_defor.tif")
```

### 4. Project future forest cover change

```{r annual-defor}
# Forest cover
fc <- vector()
dates <- c("t2", "t3")
ndates <- length(dates)
for (i in 1:ndates) {
  rast <- paste0("data/forest/forest_", dates[i], ".tif")
  val <- far$countpix(input_raster=rast, value=1)
  fc <- c(fc, val$area)  # area in ha
}
# Save results to disk
fileConn <- file("output/forest_cover.txt")
writeLines(as.character(fc), fileConn)
close(fileConn)
# Annual deforestation
T <- 10.0
annual_defor <- (fc[1] - fc[2]) / T
cat(paste0("Mean annual deforested area during the period 2010-2020: ",
		   annual_defor,
		   " ha/yr"))
```

```{r project}
# Projected deforestation (ha) during 2020-2050
defor <- annual_defor * 30

# Compute future forest cover in 2050
stats <- far$deforest(
    input_raster="output/prob.tif",
    hectares=defor,
    output_file="output/fcc_2050.tif",
    blk_rows=128L)
```

### 5. Figures

#### 5.1 Historical forest cover change

Forest cover change for the period 2000-2010-2020

```{r plot-fcc123}
# Plot forest
fig_fcc123 <- far$plot$fcc123(
    input_fcc_raster="data/forest/fcc123.tif",
    maxpixels=1e8,
    output_file="output/fcc123.png",
    borders="data/ctry_PROJ.shp",
    linewidth=0.2,
    figsize=c(5, 4), dpi=800)

knitr::include_graphics("output/fcc123.png")
```

#### 5.2 Spatial random effects

```{r plot-rho-orig}
# Original spatial random effects
fig_rho_orig <- far$plot$rho("output/rho_orig.tif",
                            borders="data/ctry_PROJ.shp",
                            linewidth=0.5,
                            output_file="output/rho_orig.png",
                            figsize=c(9,5), dpi=150)

knitr::include_graphics("output/rho_orig.png")
```

```{r plot-rho}
# Interpolated spatial random effects
fig_rho <- far$plot$rho("output/rho.tif",
                       borders="data/ctry_PROJ.shp",
                       linewidth=0.5,
                       output_file="output/rho.png",
                       figsize=c(9,5), dpi=150)

knitr::include_graphics("output/rho.png")
```

#### 5.3 Spatial probability of deforestation

```{r plot-prob}
# Spatial probability of deforestation
fig_prob <- far$plot$prob("output/prob.tif",
                         maxpixels=1e8,
                         borders="data/ctry_PROJ.shp",
                         linewidth=0.2,
                         legend=TRUE,
                         output_file="output/prob.png",
                         figsize=c(5, 4), dpi=800)

knitr::include_graphics("output/prob.png")
```

#### 5.4 Future forest cover

```{r plot-proj}
# Projected forest cover change (2020-2050)
fcc_2050 <- far$plot$fcc("output/fcc_2050.tif",
                        maxpixels=1e8,
                        borders="data/ctry_PROJ.shp",
                        linewidth=0.2,
                        output_file="output/fcc_2050.png",
                        figsize=c(5, 4), dpi=800)

knitr::include_graphics("output/fcc_2050.png")
```


Community guidelines
====================

The ``forestatrisk`` Python package is Open Source and released under
the `GNU GPL version 3 license <license.html>`__. Anybody who is
interested can contribute to the package development. There are many
ways to contribute, such as writing tutorials, examples or tests,
improving documentation, submitting bug reports and feature requests,
or writing code to provide new functionalities which can be
incorporated into future versions of the package. Every contributor
must agree to follow the project's `Code Of Conduct <code_of_conduct.html>`__.

Report an issue
+++++++++++++++

If you want to report a bug, request a feature, or discuss an issue,
please open an `issue
<https://github.com/ghislainv/forestatrisk/issues>`__ on the GitHub
project page.

Contribute to code
++++++++++++++++++

Changes to the source code and documentation should be made via GitHub
pull requests (PR).

You can learn how to do this from this *free* video series `How to
Contribute to an Open Source Project on GitHub
<https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github>`__,
Aaron Meurer's `tutorial on the git workflow
<https://www.asmeurer.com/git-workflow/>`__, or the guide `How to
Contribute to Open Source
<https://opensource.guide/how-to-contribute/>`__.

The important steps to follow are:

1. Start by creating a fork of the ``forestatrisk`` repository
2. Make changes to the source code on a development branch, not the default *master* branch
3. Keep your fork's master and development branches up to date with changes in the ``forestatrisk`` repository
4. Commit the changes you made. Chris Beams has written a `guide <https://chris.beams.io/posts/git-commit/>`__
   on how to write good commit messages.
5. Push to your fork and submit a pull request.
Changelog
=========

forestatrisk 1.0
++++++++++++++++

* Version associated to the publication in JOSS.
* Adding a Contributing section.
* Adding Community guidelines and a Code of conduct.
  
forestatrisk 0.2
++++++++++++++++

* New module organization.
* Adding a logo for the package.
* Package website at `<https://ecology.ghislainv.fr/forestatrisk/>`_\ .
* Update docstring in Python functions.
* New documentation with Sphinx.
* Continuous Integration with GitHub Actions.
* CI: Automated tests with pytest
* CI: Wheel build for PyPI.
  
forestatrisk 0.1.1
++++++++++++++++++

* New ``ee_jrc`` function to compute forest cover with GEE.
* Use of ``rclone`` to interact with GoogleDrive.
* Updated dependencies.
* Use of ``pywdpa`` to download protected areas.
* New tutorials.
* Tests have been added.

forestatrisk 0.1
++++++++++++++++

* First release of the package (previously called ``deforestprob``).
  
Code of conduct
===============

Our Pledge
----------

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our project
and our community a harassment-free experience for everyone, regardless
of age, body size, disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic
status, nationality, personal appearance, race, religion, or sexual
identity and orientation.

Our Standards
-------------

Examples of behavior that contributes to creating a positive environment
include:

-  Using welcoming and inclusive language
-  Being respectful of differing viewpoints and experiences
-  Gracefully accepting constructive criticism
-  Focusing on what is best for the community
-  Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

-  The use of sexualized language or imagery and unwelcome sexual
   attention or advances
-  Trolling, insulting/derogatory comments, and personal or political
   attacks
-  Public or private harassment
-  Publishing others' private information, such as a physical or
   electronic address, without explicit permission
-  Other conduct which could reasonably be considered inappropriate in a
   professional setting

Our Responsibilities
--------------------

Project maintainers are responsible for clarifying the standards of
acceptable behavior and are expected to take appropriate and fair
corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that they
deem inappropriate, threatening, offensive, or harmful.

Scope
-----

This Code of Conduct applies within all project spaces, and it also
applies when an individual is representing the project or its community
in public spaces. Examples of representing a project or community
include using an official project e-mail address, posting via an
official social media account, or acting as an appointed representative
at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior
may be reported by contacting the project team at *ghislain [dot]
vieilledent [at] cirad [dot] fr*.  All complaints will be reviewed and
investigated and will result in a response that is deemed necessary
and appropriate to the circumstances.  The project team is obligated
to maintain confidentiality with regard to the reporter of an
incident. Further details of specific enforcement policies may be
posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in
good faith may face temporary or permanent repercussions as determined
by other members of the project's leadership.

Attribution
-----------

This Code of Conduct is adapted from the `Contributor
Covenant <https://www.contributor-covenant.org>`__, version 1.4,
available at
https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
..
   # ==============================================================================
   # author          :Ghislain Vieilledent
   # email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
   # web             :https://ecology.ghislainv.fr
   # license         :GPLv3
   # ==============================================================================

.. image:: https://ecology.ghislainv.fr/forestatrisk/_static/logo-far.svg
   :align: right
   :target: https://ecology.ghislainv.fr/forestatrisk
   :alt: Logo forestatrisk
   :width: 140px

``forestatrisk`` Python package
*******************************


|Python version| |PyPI version| |GitHub Actions| |License| |Zenodo| |JOSS|


Overview
========

The ``forestatrisk`` Python package can be used to **model** the
tropical deforestation spatially, **predict** the spatial risk of
deforestation, and **forecast** the future forest cover in the
tropics. It provides functions to estimate the spatial probability of
deforestation as a function of various spatial explanatory variables.

Spatial explanatory variables can be derived from topography
(altitude, slope, and aspect), accessibility (distance to roads,
towns, and forest edge), deforestation history (distance to previous
deforestation), or land conservation status (eg. protected area) for
example.

.. image:: https://ecology.ghislainv.fr/forestatrisk/_static/forestatrisk.png
   :align: center
   :target: https://ecology.ghislainv.fr/forestatrisk
   :alt: prob_AFR
   :width: 800px

Scientific publication
======================

**Vieilledent G.** 2021. ``forestatrisk``: a Python package for
modelling and forecasting deforestation in the tropics.
*Journal of Open Source Software*. 6(59): 2975.
[doi: `10.21105/joss.02975 <https://doi.org/10.21105/joss.02975>`__]. |pdf|
	   
Statement of Need
=================

Spatial modelling of the deforestation allows identifying the main
factors determining the spatial risk of deforestation and quantifying
their relative effects. Forecasting forest cover change is paramount
as it allows anticipating the consequences of deforestation (in terms
of carbon emissions or biodiversity loss) under various technological,
political and socio-economic scenarios, and informs decision makers
accordingly. Because both biodiversity and carbon vary greatly in
space, it is necessary to provide spatial forecasts of forest cover
change to properly quantify biodiversity loss and carbon emissions
associated with future deforestation.

The ``forestatrisk`` Python package can be used to model the tropical
deforestation spatially, predict the spatial risk of deforestation,
and forecast the future forest cover in the tropics. The spatial data
used to model deforestation come from georeferenced raster files,
which can be very large (several gigabytes). The functions available
in the ``forestatrisk`` package process large rasters by blocks of
data, making calculations fast and efficient. This allows
deforestation to be modeled over large geographic areas (e.g. at the
scale of a country) and at high spatial resolution
(eg. ≤ 30 m). The ``forestatrisk`` package offers the possibility
of using logistic regression with auto-correlated spatial random
effects to model the deforestation process. The spatial random effects
make possible to structure the residual spatial variability of the
deforestation process, not explained by the variables of the model and
often very large. In addition to these new features, the
``forestatrisk`` Python package is open source (GPLv3 license),
cross-platform, scriptable (via Python), user-friendly (functions
provided with full documentation and examples), and easily extendable
(with additional statistical models for example). The ``forestatrisk``
Python package has been used to model deforestation and predict future
forest cover by 2100 across the humid tropics
(`<https://forestatrisk.cirad.fr>`__).

Installation
============

You will need several dependencies to run the ``forestatrisk`` Python
package. The best way to install the package is to create a Python
virtual environment, either through ``conda`` (recommended) or ``virtualenv``.

Using ``conda`` (recommended)
+++++++++++++++++++++++++++++

You first need to have ``miniconda3`` installed (see `here
<https://docs.conda.io/en/latest/miniconda.html>`__).

Then, create a conda environment (details `here
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__)
and install the ``forestatrisk`` package with the following commands:

.. code-block:: shell
		
   conda create --name conda-far -c conda-forge python=3.7 gdal numpy matplotlib pandas patsy pip statsmodels earthengine-api --yes
   conda activate conda-far
   pip install pywdpa sklearn # Packages not available with conda
   pip install forestatrisk # For PyPI version
   # pip install https://github.com/ghislainv/forestatrisk/archive/master.zip # For GitHub dev version
   # conda install -c conda-forge python-dotenv rclone --yes  # Potentially interesting libraries

To deactivate and delete the conda environment:

.. code-block:: shell
		
   conda deactivate
   conda env remove --name conda-far

Using ``virtualenv``
++++++++++++++++++++

You first need to have the ``virtualenv`` package installed (see `here <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__).

Then, create a virtual environment and install the ``forestatrisk``
package with the following commands:

.. code-block:: shell

   cd ~
   mkdir venvs # Directory for virtual environments
   cd venvs
   virtualenv --python=/usr/bin/python3 venv-far
   source ~/venvs/venv-far/bin/activate
   # Install numpy first
   pip install numpy
   # Install gdal (the correct version) 
   pip install --global-option=build_ext --global-option="-I/usr/include/gdal" gdal==$(gdal-config --version)
   pip install forestatrisk # For PyPI version, this will install all other dependencies
   # pip install https://github.com/ghislainv/forestatrisk/archive/master.zip # For GitHub dev version
   pip install statsmodels # Optional additional packages

To deactivate and delete the virtual environment:

.. code-block:: shell
		
   deactivate
   rm -R ~/venvs/venv-far # Just remove the repository

Installation testing
++++++++++++++++++++

You can test that the package has been correctly installed using the
command ``forestatrisk`` in a terminal:

.. code-block:: shell

  forestatrisk

This should return a short description of the ``forestatrisk`` package
and the version number:

.. code-block:: shell

  # forestatrisk: modelling and forecasting deforestation in the tropics.
  # https://ecology.ghislainv.fr/forestatrisk/
  # forestatrisk version x.x.

You can also test the package executing the commands in the `Get
started
<https://ecology.ghislainv.fr/forestatrisk/notebooks/get_started.html>`__
tutorial.
   
Main functionalities
====================

Sample
++++++

Function ``.sample()`` sample observations points from a forest cover
change map. The sample is balanced and stratified between deforested
and non-deforested pixels. The function also retrieves information
from explanatory variables for each sampled point. Sampling is done by
block to allow computation on large study areas (e.g. country or
continental scale) with a high spatial resolution (e.g. 30m).

Model
+++++

Function ``.model_binomial_iCAR()`` can be used to fit the
deforestation model. A linear Binomial logistic regression model is
used in this case. The model includes an intrinsic Conditional
Autoregressive (iCAR) process to account for the spatial
autocorrelation of the observations. Parameter inference is done in a
hierarchical Bayesian framework. The function calls a Gibbs sampler
with a Metropolis algorithm written in pure C code to reduce
computation time.

Other models (such as a simple GLM or a Random Forest model) can also
be used.

Predict and project
+++++++++++++++++++

Function ``.predict()`` allows predicting the deforestation
probability on the whole study area using the deforestation model
fitted with ``.model_*()`` functions. The prediction is done by block
to allow the computation on large study areas (e.g. country or
continental scale) with a high spatial resolution (e.g. 30m).

Function ``.deforest()`` predicts the future forest cover map based on a
raster of probability of deforestation (rescaled from 1 to 65535),
which is obtained from function ``.predict()``, and an area (in
hectares) to be deforested.

Validate
++++++++

A set of functions (eg. ``.cross_validation()`` or
``.map_accuracy()``\ ) is also provided to perform model and map
validation.

Contributing
============

The ``forestatrisk`` Python package is Open Source and released under
the `GNU GPL version 3 license
<https://ecology.ghislainv.fr/forestatrisk/license.html>`__. Anybody
who is interested can contribute to the package development following
our `Community guidelines
<https://ecology.ghislainv.fr/forestatrisk/contributing.html>`__. Every
contributor must agree to follow the project's `Code of conduct
<https://ecology.ghislainv.fr/forestatrisk/code_of_conduct.html>`__.


.. |Python version| image:: https://img.shields.io/pypi/pyversions/forestatrisk?logo=python&logoColor=ffd43b&color=306998
   :target: https://pypi.org/project/forestatrisk
   :alt: Python version

.. |PyPI version| image:: https://img.shields.io/pypi/v/forestatrisk
   :target: https://pypi.org/project/forestatrisk
   :alt: PyPI version

.. |GitHub Actions| image:: https://github.com/ghislainv/forestatrisk/workflows/PyPkg/badge.svg
   :target: https://github.com/ghislainv/forestatrisk/actions
   :alt: GitHub Actions
	 
.. |License| image:: https://img.shields.io/badge/licence-GPLv3-8f10cb.svg
   :target: https://www.gnu.org/licenses/gpl-3.0.html
   :alt: License GPLv3	 

.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.996337.svg
   :target: https://doi.org/10.5281/zenodo.996337
   :alt: Zenodo

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.02975/status.svg
   :target: https://doi.org/10.21105/joss.02975
   :alt: JOSS

.. |pdf| image:: https://ecology.ghislainv.fr/forestatrisk/_static/logo-pdf.png
   :target: https://www.theoj.org/joss-papers/joss.02975/10.21105.joss.02975.pdf
   :alt: pdf
:orphan:
   
.. include:: ../CONTRIBUTING.rst
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Python API
==========

Supackages and submodules of the ``forestatrisk`` Python package.

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   subpackages/forestatrisk.build_data
   subpackages/forestatrisk.model
   subpackages/forestatrisk.predict
   subpackages/forestatrisk.project
   subpackages/forestatrisk.validate
   subpackages/forestatrisk.plot
   subpackages/forestatrisk.misc

Submodules
----------

.. toctree::
   :maxdepth: 4
	      
   submodules/forestatrisk.forestatrisk
   submodules/forestatrisk.hbm

Citation
========

Please cite `Vieilledent (2021)
<https://doi.org/10.21105/joss.02975>`__ if you use the
``forestatrisk`` Python package for a scientific publication. You can
use the BibTeX entry available below.

**Vieilledent G.** 2021. ``forestatrisk``: a Python package for
modelling and forecasting deforestation in the tropics.
*Journal of Open Source Software*. 6(59): 2975.
[doi: `10.21105/joss.02975 <https://doi.org/10.21105/joss.02975>`__]. |pdf|

.. code-block:: bibtex

    @article{Vieilledent2021,
      author    = {Ghislain Vieilledent},
      title     = {forestatrisk: a {P}ython package for modelling and forecasting deforestation in the tropics},
      journal   = {{Journal of Open Source Software}},
      year      = {2021},
      volume    = {6},
      number    = {59},
      pages     = {2975},
      doi       = {10.21105/joss.02975},
      url       = {https://doi.org/10.21105/joss.02975},
      publisher = {The Open Journal},
    }

    
.. |pdf| image:: https://ecology.ghislainv.fr/forestatrisk/_static/logo-pdf.png
   :target: https://www.theoj.org/joss-papers/joss.02975/10.21105.joss.02975.pdf
   :alt: pdf
.. include:: ../README.rst

Table of contents
================= 

.. toctree::
   :maxdepth: 2

   Home <self>	      
   notebooks/get_started
   articles
   package_contents      
   indices
   citation
   changelog
   license
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

:orphan:

.. include:: ../CODE_OF_CONDUCT.rst
.. include:: ../CHANGELOG.rst
Articles
========

Python notebooks
----------------

.. toctree::
   :maxdepth: 1

   notebooks/get_started.ipynb
   notebooks/country_data/country_data
   notebooks/far_tropics.ipynb
   notebooks/new_caledonia/new_caledonia

   
R notebooks
-----------

.. toctree::
   :maxdepth: 1

   notebooks_R/get_started_R

   
forestatrisk.hbm module
-----------------------

.. automodule:: forestatrisk.hbm
   :members:
   :undoc-members:
   :show-inheritance: 
forestatrisk.forestatrisk module
--------------------------------

.. automodule:: forestatrisk.forestatrisk
   :members:
   :undoc-members:
   :show-inheritance: 
forestatrisk.build\_data package
================================

forestatrisk.build\_data.data module
------------------------------------

.. automodule:: forestatrisk.build_data.data
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.build\_data.ee\_gfc module
---------------------------------------

.. automodule:: forestatrisk.build_data.ee_gfc
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.build\_data.ee\_jrc module
---------------------------------------

.. automodule:: forestatrisk.build_data.ee_jrc
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.build\_data.sample module
--------------------------------------

.. automodule:: forestatrisk.build_data.sample
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.predict package
============================

forestatrisk.predict.interpolate\_rho module
--------------------------------------------

.. automodule:: forestatrisk.predict.interpolate_rho
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.predict.predict\_raster module
-------------------------------------------

.. automodule:: forestatrisk.predict.predict_raster
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.predict.predict\_raster\_binomial\_iCAR module
-----------------------------------------------------------

.. automodule:: forestatrisk.predict.predict_raster_binomial_iCAR
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.predict.wrast\_rho module
--------------------------------------

.. automodule:: forestatrisk.predict.wrast_rho
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.project package
============================

forestatrisk.project.deforest module
------------------------------------

.. automodule:: forestatrisk.project.deforest
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.project.deforest\_diffusion module
-----------------------------------------------

.. automodule:: forestatrisk.project.deforest_diffusion
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.project.emissions module
-------------------------------------

.. automodule:: forestatrisk.project.emissions
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.validate package
=============================

forestatrisk.validate.diffproj module
-------------------------------------

.. automodule:: forestatrisk.validate.diffproj
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.validate.map\_accuracy module
------------------------------------------

.. automodule:: forestatrisk.validate.map_accuracy
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.validate.map\_validation module
--------------------------------------------

.. automodule:: forestatrisk.validate.map_validation
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.validate.model\_validation module
----------------------------------------------

.. automodule:: forestatrisk.validate.model_validation
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.validate.resample\_sum module
------------------------------------------

.. automodule:: forestatrisk.validate.resample_sum
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.plot package
=========================

forestatrisk.plot.plot module
-----------------------------

.. automodule:: forestatrisk.plot.plot
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.misc package
=========================

forestatrisk.misc.countpix module
---------------------------------

.. automodule:: forestatrisk.misc.countpix
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.misc.miscellaneous module
--------------------------------------

.. automodule:: forestatrisk.misc.miscellaneous
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.model package
==========================

forestatrisk.model.cellneigh module
-----------------------------------

.. automodule:: forestatrisk.model.cellneigh
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.model.model\_binomial\_iCAR module
-----------------------------------------------

.. automodule:: forestatrisk.model.model_binomial_iCAR
   :members:
   :undoc-members:
   :show-inheritance:

forestatrisk.model.model\_random\_forest module
-----------------------------------------------

.. automodule:: forestatrisk.model.model_random_forest
   :members:
   :undoc-members:
   :show-inheritance:

