[![docs](https://github.com/JGCRI/plutus/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/JGCRI/plutus/actions/workflows/pkgdown.yaml)
[![build](https://github.com/JGCRI/plutus/actions/workflows/rcmd.yml/badge.svg)](https://github.com/JGCRI/plutus/actions/workflows/rcmd.yml)
[![codecov](https://codecov.io/gh/JGCRI/plutus/branch/main/graph/badge.svg?token=1PK34KIHKE)](https://codecov.io/gh/JGCRI/plutus)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03212/status.svg)](https://doi.org/10.21105/joss.03212)

# plutus
Plutus is designed for GCAM v5.3 (excluding GCAM-USA).
<br />

<!-------------------------->
<!-------------------------->
# <a name="Contents"></a>Contents
<!-------------------------->
<!-------------------------->

- [Key Links](#KeyLinks)
- [Introduction](#Introduction)
- [Citation](#Citation)
- [Installation Guide](#InstallGuides)
- [How-to Guides](#How-toGuides) 

<br />

<!-------------------------->
<!-------------------------->
# <a name="KeyLinks"></a>Key Links
<!-------------------------->
<!-------------------------->

- Github: https://github.com/JGCRI/plutus
- Webpage: https://jgcri.github.io/plutus/

[Back to Contents](#Contents)

<br />

<!-------------------------->
<!-------------------------->
# <a name="Introduction"></a>Introduction
<!-------------------------->
<!-------------------------->

`plutus` post-processes outputs from the Global Change Analysis Model (GCAM) to calculate the electricity investment costs and stranded asset costs associated with GCAM projections of future power sector energy generation by technology.


[Back to Contents](#Contents)

<br />

<!-------------------------->
<!-------------------------->
# <a name="Citation"></a>Citation
<!-------------------------->
<!-------------------------->

Zhao, M., Binsted, M., Wild, T.B., Khan, Z., Yarlagadda, B., Iyer, G., Vernon, C., Patel, P., Santos da Silva, S.R., Calvin, K.V., (2021). plutus - An R package to calculate electricity investments and stranded assets from the Global Change Analysis Model (GCAM). Journal of Open Source Software, 6(65), 3212, https://doi.org/10.21105/joss.03212


[Back to Contents](#Contents)

<br />


<!-------------------------->
<!-------------------------->
# <a name="InstallationGuides"></a>Installation Guides
<!-------------------------->
<!-------------------------->

1. Download and install:

    - R (https://www.r-project.org/)
    - R studio (https://www.rstudio.com/)

2. For Linux users, install following libraries:

```
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
sudo apt-get install libxml2-dev
```
    
3. Open R studio:

```
install.packages('devtools')
devtools::install_github('JGCRI/rgcam')
devtools::install_github('JGCRI/plutus')
```

4. `Metis` installation

`Metis` provides functions to visualize the outputs from `plutus`. The installation guide for `Metis` can be accessed at [Metis Github Page](https://github.com/JGCRI/metis).

[Back to Contents](#Contents)

<br />


<!-------------------------->
<!-------------------------->
# <a name="How-toGuides"></a>How-to Guides
<!-------------------------->
<!-------------------------->
`plutus::gcamInvest` provides all-in-one workflow that reads gcamdata, processes queries, and estimates stranded assets and capital investments. Please visit the followings for detailed instructions.

- [Instruction on `plutus::gcamInvest`](https://jgcri.github.io/plutus/articles/gcamInvest.html)
- [Case tutorial](https://jgcri.github.io/plutus/articles/CaseTutorial.html)

[Back to Contents](#Contents)

<br />
<!-- ------------------------>
<!-- ------------------------>
# plutus 0.1.0 (Under development)
<p align="center"> <img src="READMEfigs/plutus.PNG"></p>
<!-- ------------------------>
<!-- ------------------------>

## New features

## Bug Fixes
# How to Contribute to `plutus`

Thank you for taking the time to contribute and help us advance the science and architecture of `plutus`. We provide few guidelines that we ask contributors to follow. The guidelines aim to ease the maintainers' organizational and logistical duties, while encouraging development by others.

Before you start:

* Make sure you have a [GitHub account](https://github.com/signup/free).
* Trivial changes to comments or documentation do not require creating a new issue.

## Did you find a bug?

* Make sure the bug was not already reported in the Github [Issues](https://github.com/JGCRI/plutus/issues).
* [Open an issue](https://github.com/JGCRI/plutus/issues/new) and clearly describe the issue with as much information as possible. A code sample or an executable test case are recommended.
  
## Did you plan to write a patch that fixes a bug?

  * [Open an issue](https://github.com/JGCRI/plutus/issues/new) and clearly describes the problem and discuss how your solution will affect `plutus`.
  * Fork the repository on GitHub to work on the patch.
  * Interact with the project maintainers to refine/change/prioritize your issue.

## Making changes

* Start your work on your fork of the repository.
* Check for unnecessary whitespace with `git diff --check` and format code.
* Make sure your commit messages are descriptive but succinct, describing what was changed and why, and **reference the relevant issue number**. Make commits of logical units.
* Make sure you have added the necessary tests for your changes.
* Run **all** the tests to assure nothing else was accidentally broken.

## Submitting changes

* Submit a pull request with clear documentation of the methodology to the main `plutus` repository.
* **Your pull request should include one of the following two statements**:
   * You own the copyright on the code being contributed, and you hereby grant PNNL unlimited license to use this code in this version or any future version of `plutus`. You reserve all other rights to the code.
   * Somebody else owns the copyright on the code being contributed (e.g., your employer because you did it as part of your work for them); you are authorized by that owner to grant PNNL an unlimited license to use this code in this version or any future version of `plutus`, and you hereby do so. All other rights to the code are reserved by the copyright owner.
* The core team looks at Pull Requests, and will respond as soon as possible.

# Additional Resources

* [General GitHub documentation](http://help.github.com/)
* [GitHub pull request documentation](http://help.github.com/send-pull-requests/)
---
title: "Case Tutorial"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Case Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


<!-------------------------->
<!-------------------------->
## Installation
<!-------------------------->
<!-------------------------->
<p align="center"> <img src="vignettesFigs/divider.png"></p>

1. Download and install:

    - R (https://www.r-project.org/)
    
    - R studio (https://www.rstudio.com/)

2. For Linux users, install following libraries:

```
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
sudo apt-get install libxml2-dev
```
    
3. Open R studio:

```
install.packages('devtools')
devtools::install_github('JGCRI/rgcam')
devtools::install_github('JGCRI/plutus')
```


<br />

<!-------------------------->
<!-------------------------->
## Example
<!-------------------------->
<!-------------------------->
<p align="center"> <img src="vignettesFigs/divider.png"></p>

This example will give you a step-by-step instruction on using the core functionality of `plutus`.

#### 1. Download example dataset

  - Download GCAM sample XML database at http://doi.org/10.5281/zenodo.4641500
  - Unzip dataset to your desired location. Make sure the path is the one that holds the dataset (.basex files). In this example, the path to the dataset is: E:/plutus_example/IDBNexus_gcam5p3_HadGEM2-ES_rcp8p5

#### 2. Prepare your own input data and assumptions

The default data and assumptions used in `plutus` can be found in GCAM v5.3 (See [Table 1](#table1)). However, if any of those files listed in [Table 1](#table1) have been updated in your GCAM run, it is necessary to use the updated values in `plutus` as well. `plutus` uses argument `gcamdataFile` to specify the path to ALL of the updated files.

<a name="table1"></a>
**Table 1:** Data and assumption files. "~\\" represents the path to your GCAM v5.3 home folder.

| Data or Assumption | Technology | Region | Data File in GCAM v5.3|
|---|---|---|---|
| Overnight capital costs | Electricity generation technologies | Global | ~\\input\\gcamdata\\outputs\\L2233.GlobalIntTechCapital_elec.csv <br /> ~\\input\\gcamdata\\outputs\\L2233.GlobalTechCapital_elecPassthru.csv |
| Overnight capital costs | Cooling technologies | Global | ~\\input\\gcamdata\\outputs\\L2233.GlobalIntTechCapital_elec_cool.csv <br /> ~\\input\\gcamdata\\outputs\\L2233.GlobalTechCapital_elec_cool.csv |
| Capacity factors | Electricity generation technologies | Global | ~\\input\\gcamdata\\outputs\\L223.GlobalTechCapFac_elec.csv |
| Capacity factors | Intermittant technologies | Global | ~\\input\\gcamdata\\outputs\\L223.GlobalIntTechCapFac_elec.csv |
| Capacity factors | Intermittant technologies | Regional | ~\\input\\gcamdata\\outputs\\L223.StubTechCapFactor_elec.csv |
| Lifetime and steepness| Electricity generation technologies | Global | ~\\input\\gcamdata\\inst\\extdata\\energy\\A23.globaltech_retirement.csv |

*Note*: In this example, our GCAM sample dataset is based on the default data and assumptions, so we do not need to update those files in `plutus`.

#### 3.  Use `plutus::gcamInvest` to estimate stranded assets and electricity investments

Here are the core arguments we used in this example:

  - `gcamdatabase`: FULL path to the GCAM output XML database folder.
  - `reReadData`: If TRUE, `plutus` will read the GCAM database and create a queryData.proj file under the output directory.
  - `scenOrigNames`: Choose scenarios to read in. Default = 'All' will read all the scenarios.
  - `regionsSelect`: Choose regions to read in. Defaul = NULL will read all the regions.
  - `saveData`: If TRUE, `plutus` will save data by different aggregated classes.


```{r eval=F, results='hide'}
# Load required packages
library(plutus)
library(dplyr)

# Set your directory based on the example dataset location
workdir <- 'E:/plutus_example/'

# provide path to the desired GCAM database folder that holds .basex files.
path_to_gcamdatabase <- file.path(workdir, 'IDBNexus_gcam5p3_HadGEM2-ES_rcp8p5')
# Specify the path to data file folder if any data files has been updated
# path_to_gcamdataFile <- 'E:/plutus_example/gcamdataFile'

# Use plutus::gcamInvest to calculate stranded assets and electricity investments
invest <- plutus::gcamInvest(gcamdatabase = path_to_gcamdatabase,
                             dataProjFile = file.path(workdir, 'outputs', 'dataProj.proj'),
                             dirOutputs = file.path(workdir, 'outputs'),
                             reReadData = T,
                             # gcamdataFile = path_to_gcamdataFile, # Use this argument if any data files has been updated
                             scenOrigNames = c('Reference'),
                             regionsSelect = c('USA', 'China'),
                             saveData = T)
```

```{r eval=T, echo=F, results='hide', include=FALSE}
library(plutus)
library(dplyr)

invest <- plutus::gcamInvest(gcamdatabase = NULL,
                             dataProjFile = plutus::exampleGCAMproj,
                             scenOrigNames = 'Reference',
                             regionsSelect = c('USA', 'China'),
                             saveData = F)
```

<br />

#### 4. Check outputs

`plutus::gcamInvest` returns a list containing:

- `data`: a dataframe with the post-processed GCAM output showing stranded assets and electricity investments (see [Table 2](#table2)) by scenario, region, technology, and time period 
- `dataAggParam`: a dataframe with the data aggregated to the parameter
- `dataAggClass1`: a dataframe with the data aggregated to class 1
- `dataAggclass2`: a dataframe with the data aggregated to class 2
- `scenarios`: A list of the scenarios
- `queries`: A list of the queries used to extract the data

`data` contains the most detailed information. Check the structure and content of `data` and select the most relevant columns including scenario, region, param, class1 (= fuel types), x (= year), units, and value.

```{r eval=T}
# Get dataframe with post-processed outputs for stranded assets and electricity investment
df_invest <- invest$data
# Check dataframe structure
str(df_invest)

# Select key variables
df_invest_sub <- df_invest %>% 
  dplyr::select(scenario, region, param, class1, x, units, value) %>% 
  dplyr::rename(fuel = class1,
                year = x)
head(df_invest_sub)

# Check out 8 parameters
unique(df_invest$param)

```


There are **eight** different parameters in the output showing stranded assets or electricity investments in terms of monetary value (Billion 2010 USD) or installed capacity (GW). The detailed descriptions are in [Table 2](#table2). GCAM v5.3 is operated in five-years time step with 2015 as the final calibration year (base year). Annual value is the value calculated at each time step. Cumulative value represents the cumulative total over five-year time step.

<a name="table2"></a>
**Table 2:** Descriptions of output parameters for stranded assets and electricity investments.

| Parameter | Description | Unit |
|---|---|---|
| elecNewCapCost | Annual electricity capacity installations | Billion 2010 USD |
| elecNewCapGW | Annual electricity capacity installations | Gigawatts |
| elecAnnualRetPrematureCost | Annual premature retirements | Billion 2010 USD |
| elecAnnualRetPrematureGW | Annual premature retirements | Gigawatts |
| elecCumCapCost | Cumulative electricity capacity installations | Billion 2010 USD |
| elecCumCapGW | Cumulative electricity capacity installations | Gigawatts |
| elecCumRetPrematureCost | Cumulative premature retirements | Billion 2010 USD |
| elecCumRetPrematureGW | Cumulative premature retirements | Gigawatts |

[Back to Top](#installation)
---
title: "gcamInvest"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{gcamInvest}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<!-------------------------->
<!-------------------------->
# Structure
<!-------------------------->
<!-------------------------->

<p align="center"> <img src="vignettesFigs/divider.png"></p>

`plutus::gcamInvest` integrates the functionality of reading CGAM database and calculating stranded assets and electricity investments in the power sector. This function will return a list containing:

- data: a dataframe with the post-processed GCAM output showing stranded assets and electricity investment (see [Table 1](#table1)) by scenario, region, technology, and time period 
- dataAggParam: a dataframe with the data aggregated to the parameter
- dataAggClass1: a dataframe with the data aggregated to class 1
- dataAggclass2: a dataframe with the data aggregated to class 2
- scenarios: A list of the scenarios
- queries: A list of the queries used to extract the data

<br />

<a name="table1"></a>
**Table 1:** Descriptions of output parameters for stranded assets and electricity investments.

| Parameter | Description | Unit |
|---|---|---|
| elecNewCapCost | Annual electricity capacity installations | Billion 2010 USD |
| elecNewCapGW | Annual electricity capacity installations | Gigawatts |
| elecAnnualRetPrematureCost | Annual premature retirements | Billion 2010 USD |
| elecAnnualRetPrematureGW | Annual premature retirements | Gigawatts |
| elecCumCapCost | Cumulative electricity capacity installations | Billion 2010 USD |
| elecCumCapGW | Cumulative electricity capacity installations | Gigawatts |
| elecCumRetPrematureCost | Cumulative premature retirements | Billion 2010 USD |
| elecCumRetPrematureGW | Cumulative premature retirements | Gigawatts |


The details of arguments and their default values in `plutus::gcamInvest` can be found in the [gcamInvest reference page](https://jgcri.github.io/plutus/reference/gcamInvest.html). The following sections provide step-by-step instructions on using `plutus::gcamInvest`.

```{r eval=F}
# Default argument values
plutus::gcamInvest(gcamdatabase = NULL,
                   queryFile = NULL,
                   reReadData = T,
                   dataProjFile = paste(getwd(), "/outputs/dataProj.proj", sep = ""),
                   gcamdataFile = NULL,
                   scenOrigNames = 'All',
                   scenNewNames = NULL,
                   regionsSelect = NULL,
                   dirOutputs = paste(getwd(), "/outputs", sep = ""),
                   folderName = NULL,
                   nameAppend = "",
                   saveData = T)
```

<br />

<!-------------------------->
<!-------------------------->
# Read GCAM Data
<!-------------------------->
<!-------------------------->
<p align="center"> <img src="vignettesFigs/divider.png"></p>

## Read .proj File

`plutus::gcamInvest` is able to read rgcam-based .proj file from GCAM output by providing the path to the .proj file. `plutus` also includes an example .proj dataset `plutus::exampleGCAMproj`.

```{r eval=F}
library(plutus)

invest <- plutus::gcamInvest(# dataProjFile = path_to_projfile,
                             dataProjFile = plutus::exampleGCAMproj)

# Explore the list returned from plutus::gcamInvest
df <- invest$data; df
dfParam <- invest$dataAggParam; dfParam
dfClass1 <- invest$dataAggClass1; dfClass1
dfScenario <- invest$scenarios; dfScenario
dfQuery <- invest$queries; dfQuery

```

<br />

## Read GCAM XML Database

`plutus::gcamInvest` can directly read GCAM output XML database. Assign the path of GCAM database to `gcamdatabase` argument.


```{r eval=F}
library(plutus)

# provide path to the desired GCAM database folder.
path_to_gcamdatabase <- 'E:/gcam-core-gcam-v5.3/output/databse_basexdb'
invest <- plutus::gcamInvest(gcamdatabase = path_to_gcamdatabase)

```

<br />

<!-------------------------->
<!-------------------------->
## Subset GCAM Data
<!-------------------------->
<!-------------------------->

`plutus::gcamInvest` provides options to subset GCAM data by scenario and region.


```{r eval=F}
library(plutus)

invest <- plutus::gcamInvest(dataProjFile = plutus::exampleGCAMproj,
                             scenOrigNames = c('Reference', 'Impacts', 'Policy'),
                             scenNewNames = c('Reference', 'Climate Impacts', 'Climate Policy'),
                             regionsSelect = c('USA', 'Argentina'))

df <- invest$data; df
dfParam <- invest$dataAggParam; dfParam
dfClass1 <- invest$dataAggClass1; dfClass1

```

<br />

<!-------------------------->
<!-------------------------->
# Input Data and Assumptions
<!-------------------------->
<!-------------------------->
<p align="center"> <img src="vignettesFigs/divider.png"></p>

`plutus` provides default data and assumptions files from GCAM v5.3 that are collected from Default-GCAM-5.3-Folder/input/gcamdata/. All data files (CSV) associated with different data categories and assumptions are listed in [Table 2](#table2).

<a name="table2"></a>
**Table 2:** Data and assumption files.

| Data or Assumption | Technology | Region | Data File |
|---|---|---|---|
| Overnight capital costs | Electricity generation technologies | Global | L2233.GlobalIntTechCapital_elec.csv <br /> L2233.GlobalTechCapital_elecPassthru.csv |
| Overnight capital costs | Cooling technologies | Global | L2233.GlobalIntTechCapital_elec_cool.csv <br /> L2233.GlobalTechCapital_elec_cool.csv |
| Capacity factors | Electricity generation technologies | Global | L223.GlobalTechCapFac_elec.csv |
| Capacity factors | Intermittant technologies | Global | L223.GlobalIntTechCapFac_elec.csv |
| Capacity factors | Intermittant technologies | Regional | L223.StubTechCapFactor_elec.csv |
| Lifetime and steepness| Electricity generation technologies | Global | A23.globaltech_retirement.csv |

However, users may run GCAM v5.3 with updated/customized values associated with capital costs, capacity factors, and lifetime by technology and region. In this case,  `plutus::gcamInvest` allows users to input their own data associated with their GCAM runs. All you need to do is to provide the path that includes **ALL** CSV files listed in [Table 2](#table2). Here are two options:

- Users can update CSV files listed in [Table 2](#table2) directly in Your-GCAM-5.3-Folder/input/gcamdata. Then, specify argument to `gcamdataFile = Your-GCAM-5.3-Folder/input/gcamdata`.
- Users can copy and paste all CSV files listed in  [Table 2](#table2) from Your-GCAM-5.3-Folder/input/gcamdata to a new folder and update the values according to your GCAM run. Then, specify argument to `gcamdataFile = Path-to-Your-New-Folder`.

While updating the values in those CSV files, please keep the file name and the format of each updated file unchanged.

```{r eval=F}
library(plutus)

invest <- plutus::gcamInvest(dataProjFile = plutus::exampleGCAMproj,
                             gcamdataFile = 'E:/gcam-core-gcam-v5.3/input/gcamdata')
```

<br />

<!-------------------------->
<!-------------------------->
# Output Options
<!-------------------------->
<!-------------------------->
<p align="center"> <img src="vignettesFigs/divider.png"></p>

`plutus::gcamInvest` automatically saves outputs, such as .proj files, queries, and output data tables in default path based on the working directory. This section will introduce arguments in `plutus::gcamInvest` that control different output options.

## Re-read Data

`plutus::gcamInvest` automatically save a copy of .proj file (default name: dataProj.proj) in the output directory after extracting GCAM output data if provided a GCAM database folder using `gcamdatabase` argument. Reload the same (subsetted) GCAM dataset can be much faster by using the automatically saved dataProj.proj file. Users can choose to turnoff auto save function by specifying argument to `reReadData = F`.

```{r eval=F}
library(plutus)

path_to_gcamdatabase <- 'E:/gcam-core-gcam-v5.3/output/databse_basexdb'
invest <- plutus::gcamInvest(gcamdatabase = path_to_gcamdatabase,
                             reReadData = F) # Default is reReadData = T
```

<br />

## Output Directory

It is easy to change output folder names and names appended to the output data tables.

  - `folderName` creates a folder with the specified name under output directory 
  - `nameAppend` appends the specified string at the end of each output data table file 

```{r eval=F}
library(plutus)

invest <- plutus::gcamInvest(dataProjFile = plutus::exampleGCAMproj,
                             # dirOutputs = Your-desired-output-path, # Default is paste(getwd(), "/outputs", sep = "")
                             folderName = 'USA', # Default is fodlerName = NULL
                             nameAppend = '_Invest') # Default is nameAppend = ''
```

You can also choose not to save any output by setting argument to `saveData = F`.
```{r eval=F}
library(plutus)

invest <- plutus::gcamInvest(dataProjFile = plutus::exampleGCAMproj,
                             saveData = F) # Default is saveData = T
```

<br />

<!-------------------------->
<!-------------------------->
# Visualization
<!-------------------------->
<!-------------------------->
<p align="center"> <img src="vignettesFigs/divider.png"></p>

`plutus` is designed to integrate with `metis` (Khan et al., 2020), an R package to harmonize and analyze multi-sectoral data and linkages at variable spatial scales. `plutus::gcamInvest` generates a data structure that can be directly used by `metis`. More details on metis can be accessed via metis repository at https://github.com/JGCRI/metis.

## Charts

Below is an example of using `metis.chartsProcess` to generate charts of stranded assets and electricity investments in the USA and Argentina.

```{r eval=F}
library(plutus)
library(metis)

invest <- plutus::gcamInvest(dataProjFile = plutus::exampleGCAMproj,
                             scenOrigNames = c('Reference', 'Impacts', 'Policy'),
                             scenNewNames = c('Reference', 'Climate Impacts', 'Climate Policy'),
                             regionsSelect = c('USA', 'Argentina'))
rTable_i <- invest$data

# Plot charts for each scenario and comparison between scenarios and regions across time series
metis.chartsProcess(rTable = rTable_i,
                    paramsSelect = 'All',
                    regionsSelect = c('USA', 'Argentina'),
                    xCompare = c("2030", "2040", "2050"),
                    scenRef = "Reference",
                    dirOutputs = paste(getwd(), "/outputs", sep = ""),
                    regionCompare = 1,
                    scenarioCompareOnly = 1,
                    multiPlotOn = T,
                    multiPlotFigsOnly = T,
                    folderName = "Charts_USA_Argentina",
                    xRange = c(2025, 2030, 2035, 2040, 2045, 2050),
                    colOrderName1 = "scenario",
                    pdfpng = 'pdf')
```

<p align="center" style="font-size:18px;"> *New Electricity Capacity Installations for USA and Argentina* </p>
<p align="center"> <img src="vignettesFigs/elecNewCapCost_figBar_USA_Argentina_Reference.png"></p>

<br />

## Maps

Users can also use `metis.mapsProcess` to produce spatial maps by scenario, region, technology, and time period. Here is an example to map wind capacity installations for 32 geopolitical regions.

```{r eval=F}
library(plutus)
library(metis)

invest <- plutus::gcamInvest(dataProjFile = plutus::exampleGCAMproj,
                             scenOrigNames = 'Reference',
                             regionsSelect = NULL,
                             saveData = F)

rTable_i <- invest$data
# Filter data to only show new installations for wind
polygonTable_i <- rTable_i %>% 
  dplyr::filter(class1 %in% 'Wind', param %in% 'elecNewCapGW') %>% 
  dplyr::select(subRegion, value, x, class1, scenario, param) %>% 
  dplyr::rename(class = class1)

# Plot wind new installation maps for all 32 geopolitical regions across time
metis.mapsProcess(polygonTable = polygonTable_i,
                  subRegCol = 'subRegion',
                  subRegType = 'subRegion',
                  folderName = 'Maps_Wind_eleNewCapGW',
                  facetCols = 2,
                  animateOn = F)
```

<p align="center" style="font-size:18px;"> *Regional Electricity Capacity Installations for Wind Technology* </p>
<p align="center"> <img src="vignettesFigs/map_GCAMReg32_elecNewCapGW_Reference_KMEANS.png"></p>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_capac_fac_region}
\alias{data_capac_fac_region}
\title{data_capac_fac_region}
\format{
.csv
}
\source{
paste(rawDataFolder,"L223.StubTechCapFactor_elec.csv", sep="")
}
\usage{
data_capac_fac_region
}
\description{
data_capac_fac_region
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_capac_fac_region
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elecInvest.R
\name{elecInvest}
\alias{elecInvest}
\title{elecInvest}
\usage{
elecInvest(
  elec_gen_vintage,
  gcamdataFile,
  world_regions,
  start_year = 2015,
  end_year = 2050
)
}
\arguments{
\item{elec_gen_vintage}{Electricity vintage query result}

\item{gcamdataFile}{Default = NULL. Optional. For example, gcamdataFile = "~/gcam-core-gcam-v5.3/input/gcamdata".}

\item{world_regions}{GCAM regions for which to collect data}

\item{start_year}{Start year of time frame of interest for analysis}

\item{end_year}{end_year of time frame of interest for analysis}
}
\value{
Returns data in a form required by plutus::gcamInvest.R
}
\description{
Function that calculates electricity subsector investment requirements from a given GCAM run.
}
\keyword{infrastructure}
\keyword{investments,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_cap_cost_int_tech}
\alias{data_cap_cost_int_tech}
\title{data_cap_cost_int_tech}
\format{
.csv
}
\source{
paste(rawDataFolder,"L2233.GlobalIntTechCapital_elec.csv", sep="")
}
\usage{
data_cap_cost_int_tech
}
\description{
data_cap_cost_int_tech
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_cap_cost_int_tech
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hydroInvest.R
\name{hydroInvest}
\alias{hydroInvest}
\title{hydroInvest}
\usage{
hydroInvest(addition_costs, start_year = 2010)
}
\arguments{
\item{addition_costs}{list formatted as produced by elecInvest.R function}

\item{start_year}{start year for analysis}
}
\value{
Returns data in a form required by gcamInvest.R
}
\description{
Function that calculates electricity subsector investment requirements from a given GCAM run.
}
\keyword{infrastructure}
\keyword{investments,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcamInvest.R
\name{gcamInvest}
\alias{gcamInvest}
\title{gcamInvest}
\usage{
gcamInvest(
  gcamdatabase = NULL,
  queryFile = NULL,
  reReadData = T,
  dataProjFile = paste(getwd(), "/outputs/dataProj.proj", sep = ""),
  gcamdataFile = NULL,
  scenOrigNames = "All",
  scenNewNames = NULL,
  regionsSelect = NULL,
  dirOutputs = paste(getwd(), "/outputs", sep = ""),
  folderName = NULL,
  nameAppend = "",
  saveData = T
)
}
\arguments{
\item{gcamdatabase}{Default = NULL. Full path to GCAM database folder.}

\item{queryFile}{Defualt = NULL. When NULL plutus loads pre-saved xml file plutus::xmlElecQueries}

\item{reReadData}{Default = TRUE. If TRUE will read the GCAM data base and create a queryData.proj file
under the output directory. If FALSE will load a '.proj' file if a file
with full path is provided otherwise it will search for a dataProj.proj file in the existing
folder which may have been created from an old run.}

\item{dataProjFile}{Default = paste(getwd(), "/outputs/dataProj.proj", sep = ""). Optional. A default 'dataProj.proj' is produced if no .Proj file is specified.}

\item{gcamdataFile}{Default = NULL. Optional. For example, gcamdataFile = "~/gcam-core-gcam-v5.3/input/gcamdata".
Use full path to GCAM 'gcamdata' folder that contains costs and capacity data. Data files including:

(1) A23.globaltech_retirement.csv;
(2) L223.StubTechCapFactor_elec.csv;
(3) L223.GlobalTechCapFac_elec.csv;
(4) L2233.GlobalIntTechCapital_elec.csv;
(5) L2233.GlobalIntTechCapital_elec_cool.csv;
(6) L2233.GlobalTechCapital_elec_cool.csv;
(7) L2233.GlobalTechCapital_elecPassthru.csv}

\item{scenOrigNames}{Default = "All". Original Scenarios names in GCAM database in a string vector.
For example c('scenario1','scenario2).}

\item{scenNewNames}{Default = NULL. New Names which may be shorter and more useful for figures etc.
Default will use Original Names. For example c('scenario1','scenario2)}

\item{regionsSelect}{Default = NULL. The regions to analyze in a vector. Example c('Colombia','Argentina'). Full list:

USA, Africa_Eastern, Africa_Northern, Africa_Southern, Africa_Western, Australia_NZ, Brazil, Canada
Central America and Caribbean, Central Asia, China, EU-12, EU-15, Europe_Eastern, Europe_Non_EU,
European Free Trade Association, India, Indonesia, Japan, Mexico, Middle East, Pakistan, Russia,
South Africa, South America_Northern, South America_Southern, South Asia, South Korea, Southeast Asia,
Taiwan, Argentina, Colombia, Uruguay}

\item{dirOutputs}{Default = paste(getwd(), "/outputs", sep = ""). Full path to directory for outputs}

\item{folderName}{Default = NULL}

\item{nameAppend}{Default = "". Name to append to saved files.}

\item{saveData}{Default = TRUE. Set to F if do not want to save any data to file.}
}
\value{
Returns (1) annual and cumulative costs and power generation of
premature retired electricity infrustructure (including hydropower);
(2) annual and cumulative costs and power generation of new capacity from new electricity infrustructures (including hydropower).
}
\description{
This function connects and reads a gcamdatabase or a query data proj file and calculates
electricity and hydropower investments.
}
\keyword{database,}
\keyword{gcam}
\keyword{gcam,}
\keyword{query}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assumptions.R
\name{assumptions}
\alias{assumptions}
\title{assumptions}
\usage{
assumptions(name = NULL)
}
\arguments{
\item{name}{Default=NULL. Name of assumption object.
List of Assumptions
\itemize{
\item GCAMbaseYear,
\item convEJ2MTOE,
\item convEJ2TWh,
\item convEJ2GW,
\item convEJ2GWh,
\item convGW_kW,
\item conv_C_CO2,
\item conv_MT_GT,
\item hydro_cap_fact,
\item hydro_cost_GW,
\item convUSD_1975_2010,
\item convUSD_1975_2015,
\item conv1975USDperGJ22017USDperMWh,
\item conv1975USDperGJ22017USDperMBTU}}
}
\value{
A list of assumptions
}
\description{
This function loads holds the different assumptions used throughout the plutus package.
}
\examples{
library(plutus)
a <- plutus::assumptions("GCAMbaseYear")
}
\keyword{assumptions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_A23.globaltech_retirement}
\alias{data_A23.globaltech_retirement}
\title{data_A23.globaltech_retirement}
\format{
.csv
}
\source{
paste(rawDataFolder,"A23.globaltech_retirement.csv",sep="")
}
\usage{
data_A23.globaltech_retirement
}
\description{
data_A23.globaltech_retirement
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_A23.globaltech_retirement
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_cap_cost_int_cool}
\alias{data_cap_cost_int_cool}
\title{data_cap_cost_int_cool}
\format{
.csv
}
\source{
paste(rawDataFolder,"L2233.GlobalIntTechCapital_elec_cool.csv", sep="")
}
\usage{
data_cap_cost_int_cool
}
\description{
data_cap_cost_int_cool
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_cap_cost_int_cool
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_cap_cost_tech}
\alias{data_cap_cost_tech}
\title{data_cap_cost_tech}
\format{
.csv
}
\source{
paste(rawDataFolder,"L2233.GlobalTechCapital_elecPassthru.csv", sep="")
}
\usage{
data_cap_cost_tech
}
\description{
data_cap_cost_tech
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_cap_cost_tech
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{xmlElecQueries}
\alias{xmlElecQueries}
\title{elecQueries xml file}
\format{
.xml
}
\source{
gcam query
}
\usage{
xmlElecQueries
}
\description{
elecQueries xml file
}
\examples{
\dontrun{
 library(plutus); library(XML)
 plutus::xmlElecQueries
 # Can save xml
 XML::saveXML(plutus::xmlElecQueries, file=paste(getwd(), "/ElecQueries.xml", sep = ""))
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plutus.R
\docType{package}
\name{plutus}
\alias{plutus}
\title{plutus}
\description{
\itemize{
\item Github: https://github.com/JGCRI/plutus
\item Webpage: https://jgcri.github.io/plutus/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_cap_cost_cool}
\alias{data_cap_cost_cool}
\title{data_cap_cost_cool}
\format{
.csv
}
\source{
paste(rawDataFolder,"L2233.GlobalTechCapital_elec_cool.csv", sep="")
}
\usage{
data_cap_cost_cool
}
\description{
data_cap_cost_cool
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_cap_cost_cool
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{exampleGCAMproj}
\alias{exampleGCAMproj}
\title{Example GCAM .proj file}
\format{
R table or .csv
}
\source{
metis.readgcam() run outputs saved
}
\usage{
exampleGCAMproj
}
\description{
Example GCAM .proj file
}
\examples{
\dontrun{
 library(plutus);
 plutus::exampleGCAMproj
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_capac_fac_int}
\alias{data_capac_fac_int}
\title{data_capac_fac_int}
\format{
.csv
}
\source{
paste(rawDataFolder,"L223.GlobalIntTechCapFac_elec.csv", sep="")
}
\usage{
data_capac_fac_int
}
\description{
data_capac_fac_int
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_capac_fac_int
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_tech_mapping}
\alias{data_tech_mapping}
\title{data_tech_mapping}
\format{
.csv
}
\source{
paste(rawDataFolder,"agg_tech_mapping.csv", sep="")
}
\usage{
data_tech_mapping
}
\description{
data_tech_mapping
}
\examples{
\dontrun{
 library(plutus);
 plutus::tech_mapping
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_capac_fac}
\alias{data_capac_fac}
\title{data_capac_fac}
\format{
.csv
}
\source{
paste(rawDataFolder,"L223.GlobalTechCapFac_elec.csv", sep="")
}
\usage{
data_capac_fac
}
\description{
data_capac_fac
}
\examples{
\dontrun{
 library(plutus);
 plutus::data_capac_fac
}
}
\keyword{datasets}
