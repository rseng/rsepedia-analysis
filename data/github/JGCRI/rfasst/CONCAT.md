[![build](https://github.com/JGCRI/rfasst/actions/workflows/rcmd.yml/badge.svg?branch=main)](https://github.com/JGCRI/rfasst/actions/workflows/rcmd.yml)
[![docs](https://github.com/JGCRI/rfasst/actions/workflows/pkgdown.yaml/badge.svg?branch=main)](https://github.com/JGCRI/rfasst/actions/workflows/pkgdown.yaml)
[![codecov](https://codecov.io/gh/JGCRI/rfasst/branch/main/graph/badge.svg?token=2IBODRZKVF)](https://codecov.io/gh/JGCRI/rfasst)
[![DOI](https://zenodo.org/badge/344924589.svg)](https://zenodo.org/badge/latestdoi/344924589)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03820/status.svg)](https://doi.org/10.21105/joss.03820)


<!-- ------------------------>
<!-- ------------------------>
# <a name="Contents"></a>Contents
<!-- ------------------------>
<!-- ------------------------>

- [Key Links](#KeyLinks)
- [Introduction](#Introduction)
- [Citation](#Citation)
- [Installation Guide](#InstallGuide)
- [How-to guide](#howto) 
- [Publications](#Publications)

<!-- ------------------------>
<!-- ------------------------>
# <a name="KeyLinks"></a>Key Links
<!-- ------------------------>
<!-- ------------------------>

[Back to Contents](#Contents)

- Github: https://github.com/JGCRI/rfasst
- Webpage: https://jgcri.github.io/rfasst/

<!-- ------------------------>
<!-- ------------------------>
# <a name="Introduction"></a>Introduction
<!-- ------------------------>
<!-- ------------------------>

[Back to Contents](#Contents)

`rfasst` reports a consistent range of adverse health and agricultural effects attributable to air pollution for any scenario run by the Global Change Analysis Model ([GCAM](http://www.globalchange.umd.edu/gcam/)), by replicating the calculations of the air quality reduced-form model [TM5-FASST]( https://ec.europa.eu/jrc/en/publication/tm5-fasst-global-atmospheric-source-receptor-model-rapid-impact-analysis-emission-changes-air).


<!-- ------------------------>
<!-- ------------------------>
# <a name="Citation"></a>Citation
<!-- ------------------------>
<!-- ------------------------>

[Back to Contents](#Contents)


<!-- ------------------------>
<!-- ------------------------>
# <a name="InstallGuide"></a>Installation Guide
<!-- ------------------------>
<!-- ------------------------>

[Back to Contents](#Contents)

1. Download and install:
    - R (https://www.r-project.org/)
    - R studio (https://www.rstudio.com/)  
    - (For cloning the repo) Git (https://git-scm.com/downloads) 
    
    
2. Open R studio:

```r
install.packages("devtools")
devtools::install_github("JGCRI/rfasst")
```

(Optional) 

To clone the repository to the local machine: Git bash in the working directory (right click "Git Bash Here") -> In the Git console type:  

```r
git clone https://github.com/JGCRI/rfasst.git
```

Then, open the Rproject (rfasst.Rproj): In the Rstudio menu, click "Build -> Install and restart" (Ctrl+Shift+B)
  

<!-- ------------------------>
<!-- ------------------------>
# <a name="keyfunctions"></a> How to guides
<!-- ------------------------>
<!-- ------------------------>

[Back to Contents](#Contents)

The package consists of a set of functions divided in four different modules:
- Module 1. Emissions re-scaling: Process emissions by GCAM-region and re-scale them to TM5-FASST regions, and make some additional pollutant-related adjustments. More details in the [Module1 emissions](https://jgcri.github.io/rfasst/articles/Module1_emissions.html) page. 
- Module 2. Concentration: Estimate fine particulate matter (PM2.5) and ozone (O3) concentration levels (measured by different indicators) for each region. More details in the [Module2 concentration](https://jgcri.github.io/rfasst/articles/Module2_concentration.html) page. 
- Module 3. Health: Report adverse health effects attributable to exposure to fine particulate matter (PM2.5) and ozone (O3; M6M). More details in the [Module3 health](https://jgcri.github.io/rfasst/articles/Module3_health.html) page. 
- Module 4. Agriculture: Estimate adverse agricultural impacts associated to ozone exposure, including relative yield losses (RYLs) and production and revenue losses. More details in the [Module4 agriculture](https://jgcri.github.io/rfasst/articles/Module4_agriculture.html) page. 

In addition, the package includes some default mapping files and default values, that are read by the different functions. These can be changed by the user. Some of these constants include:
- Years to be analyzed: In the `all_years` vector, the user can select the years to be included in the analysis. All the avialble years are '2005','2010','2020','2030','2040','2050','2060','2070','2080','2090','2100'.It is no possible to add any other year, but they can be reduced if desired (for example to reduce computation time).
- GCAM crop categories to be included in the analysis
- Shares to allocate emissions between Russia Eastern (RUE) and Russia Western (RUS)
- Coefficients and/or counterfactual values for exposure-response functions applied to estimate adverse health and agricultural impacts.
- Baseline moratlity rates.
- Median values for the health impact economic assessment (Value of Statistical Life)
- Other


<!-- ------------------------>
<!-- ------------------------>
# <a name="Publications"></a>Publications
<!-- ------------------------>
<!-- ------------------------>

[Back to Contents](#Contents)

Previous to the development of this package, different studies have combined the use of GCAM and TM5-FASST:

- Sampedro, J., Smith, S.J., Arto, I., González-Eguino, M., Markandya, A., Mulvaney, K.M., Pizarro-Irizar, C. and Van Dingenen, R., 2020. Health co-benefits and mitigation costs as per the Paris Agreement under different technological pathways for energy supply. Environment international, 136, p.105513.

- Markandya, A., Sampedro, J., Smith, S.J., Van Dingenen, R., Pizarro-Irizar, C., Arto, I. and González-Eguino, M., 2018. Health co-benefits from air pollution and mitigation costs of the Paris Agreement: a modelling study. The Lancet Planetary Health, 2(3), pp.e126-e133.

- Sampedro, J., Waldhoff, S.T., Van de Ven, D.J., Pardo, G., Van Dingenen, R., Arto, I., Del Prado, A. and Sanz, M.J., 2020. Future impacts of ozone driven damages on agricultural systems. Atmospheric Environment, 231, p.117538.
<!-- ------------------------>
<!-- ------------------------>
# rfasst 1.0.0
<p align="center"> <img src="READMEfigs/metisHeaderThick.PNG"></p>
<!-- ------------------------>
<!-- ------------------------>

## New features

## Bug Fixes## How to contribute
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
title: 'rfasst: An R tool to estimate air pollution impacts on health and agriculture'
tags:
  - R
  - emissions
authors:
  - name: Jon Sampedro
    orcid: 0000-0002-2277-1530
    affiliation: 1
  - name: Zarrar Khan
    orcid: 0000-0002-8147-8553
    affiliation: 1
  - name: Chris R. Vernon
    orcid: 0000-0002-3406-6214
    affiliation: 1  
  - name: Steven J. Smith
    orcid: 0000-0003-3248-5607
    affiliation: 1
  - name: Stephanie Waldhoff
    orcid: 0000-0002-8073-0868
    affiliation: 1
  - name: Rita Van Dingenen
    orcid: 0000-0003-2521-4972
    affiliation: 2
affiliations:
 - name: Joint Global Change Research Institute, Pacific Northwest National Laboratory, College Park, MD, USA
   index: 1
 - name: European Commission, Joint Research Centre (JRC), Ispra, Italy
   index: 2
date: 17 April 2021
bibliography: paper.bib
---
# Summary
Existing scientific literature shows that health and agricultural impacts attributable to air pollution are significant and should be considered in the integrated analysis of human and Earth-system interactions.
The implementation of policies that affect the power sector, the composition of the vehicle fleet or the investments and deployment of different energy sources will in turn result in different levels of air pollution. 
Even though the various methodologies for estimating the impacts of air pollution, such as exposure-response functions, are extensively applied by the scientific community, 
they are normally not included in integrated assessment modeling outputs.

`rfasst` is an R package designed to estimate future human-health and agricultural damages attributable to air pollution using scenario outputs from the Global Change Analysis Model (GCAM), namely emission pathways and agricultural production and prices. 
The package combines these with the calculations from the TM5-FASST air quality model to estimate the associated adverse health and agricultural impacts. 
The structure of the `rfasst` package is summarized in Figure 1.

![Structure of the `rfasst` package](figure_rfasst.png)

`rfasst` can be accessed via the web at the public domain https://github.com/JGCRI/rfasst. The following code is a simplified example that shows how to run the package. In addition, we provide an R vignette step-by-step tutorial for users to get started with `rfasst` which is accessible here: [Tutorial](https://jgcri.github.io/rfasst/).

```r
install.packages("devtools")
library(devtools)
devtools::install_github("JGCRI/rfasst")
library(rfasst)

db_path<-"path/to/my-gcam-database"
query_path<-"path/to/queries_rfasst.xml"
db_name<-"my-gcam-database"
prj_name<-"my-project.dat"
scen_name<-"name of the GCAM scenario"
queries<-queries_rfasst.xml 

# Using the m1_emissions_rescale function as a representative example. 
# All the functions are listed in the documentation vignettes.
m1_emissions_rescale(db_path,query_path,db_name,prj_name,
		     scen_name,queries,saveOutput=T, map=T)  

```


# Statement of need

According to the World Health Organization's (WHO) [Ambient air quality database](https://www.who.int/data/gho/data/themes/topics/topic-details/GHO/ambient-air-pollution), more than 90% of people breathe unhealthy air at a global level.
Therefore, premature mortality associated with air pollution is one of the biggest threats for human health, accounting for more than 8 million deaths globally per year [@burnett2018global], but heavily concentrated in developing Asia.
Likewise, air pollution leads to a significant decrease of crop yields. 
Ozone, which is formed by the reaction of air pollutants with solar radiation, is considered the most hazardous pollutant for crop yields [@emberson2018ozone]. 
Current high ozone concentration levels entail substantial economic damages and would increase pressures on several measures associated with food security [@van2009global]. 
The integration of these effects into integrated assessment models, such as GCAM, can provide valuable insights for scenario analysis.

The GCAM model [@calvin2019gcam], developed at the Joint Global Change Research Institute (JGCRI), is an integrated assessment multi-sector model designed to explore human and Earth-system dynamics. 
For each scenario representing an alternate future, GCAM reports a full suite of emissions of greenhouse gases and air pollutants, by region and time period through 2100. 
GCAM outputs also include regional agricultural production projections for a range of crops, detailed in online [documentation](https://github.com/JGCRI/gcam-doc/blob/gh-pages/aglu.md). However, GCAM does not include the atmospheric 
and meteorological information required to translate the greenhouse gas and air pollutant emissions into particulate matter ($PM_{2.5}$) and ozone ($O_{3}$) concentration levels. 
This transformation from emissions to concentration is addressed by full chemistry models or by simplified air quality emulators, such as TM5-FASST [@van2018tm5].
These concentration levels are the inputs for the exposure-response functions that are normally used to calculate adverse human-health and agricultural effects associated with exposure to $PM_{2.5}$ and $O_{3}$.  

Therefore, the combined use of these models, which is the essence of `rfasst`, is a powerful methodology to estimate a consistent range of health and agricultural damages and the co-benefits associated with different strategies or climate policies, that complements existing models and methods. 
Historically, the assessment of air pollution and the subsequent impacts has been developed with global chemical transport models (CTMs) such as GEOS-Chem [@bey2001global] or the The Community Multiscale Air Quality Modeling System (CMAQ; @appel2021community). 
These models combine atmospheric, meteorological, and land-use data with anthropogenic emissions, which can be extracted from global emissions inventories, such as  The Emissions Database for Global Atmospheric Research (EDGAR; @crippa2016forty) or The Community Emissions Data System (CEDS; @hoesly2018historical), 
or generated by additional tools, namely SMOKE [@smoke] or EmissV [@schuch2018emissv]. 
While these CTMs provide detailed outputs, they are also computationally expensive, so they have limited application for either rapid assessment or evaluation of multiple scenarios.
A number of air quality tools and emulators, such as TM5-FASST [@van2018tm5], BenMap [@sacks2018environmental], SHERPA [@clappier2015new] or AirQ+ [@world2018airq], 
have been developed to allow faster estimation of air pollution levels and associated impacts for global or regional changes in emission levels. 
`rfasst` contributes to this set of tools by making the already widely used TM5-FASST tool available as an easily accessible `R` package that can be readily used for air pollution and impact assessment of alternative “what-if” global scenarios.

Prior to the development of this package, we have used GCAM and TM5-FASST to analyze these co-benefits in different studies. We showed that health co-benefits attributable to air pollution are larger than mitigation costs 
for different technological scenarios consistent with the 2°C target of the Paris Agreement [@sampedro_2020a]. 
Previously, we demonstrated that these health co-benefits outweigh mitigation costs in multiple decarbonization scenarios based on different emissions abatement efforts across regions [@markandya2018health]. 
In addition, we have applied this methodology to show how high $O_{3}$ levels generate substantial crop losses and, subsequently, negative economic impacts in the agricultural sector [@sampedro_2020b].
Taking all these results into consideration, we understand that a tool that systematically addresses air pollution driven human-health and agricultural damages within an integrated assessment modeling framework, 
is a significant contribution to this community, and of interest for a range of stakeholders, particularly for the designers of alternative transition strategies. 



# Functionality
The package includes several functions that have been classified in four different modules. 
Note that all the functions are listed in the [Tutorial](https://jgcri.github.io/rfasst/), which includes individual documentation pages for each of these modules.

+ Module 1: Static downscaling of GCAM emissions to country-level and re-aggregation into a new regional distribution (consistent with TM5-FASST), and some additional pollutant-related adjustments (e.g., organic carbon to organic matter).
+ Module 2: Calculation of regional fine particulate matter ($PM_{2.5}$) and ozone ($O_{3}$) concentration levels using different indicators.
+ Module 3: Estimation of health impacts attributable to $PM_{2.5}$ and $O_{3}$ exposure. The package reports both physical damages, such as premature mortality, years of life lost (YLLs), and disability adjusted life years (DALYs),
and the associated monetized damages based on the Value of Statistical Life (VSL).
+ Module 4: Estimation of agricultural damages attributable to $O_{3}$ exposure, including relative yield losses (RYLs) and losses in agricultural production and revenue ($Revenue=Prod \cdot Price$).

The package also includes additional input information, namely constant values and mapping files, that need to be read in for running the different functions and can be modifiable by the user.
The [Tutorial](https://jgcri.github.io/rfasst/) explains which values can be changed within each module. These include the time horizon (from 2010 to 2100 in 10-year periods, +2005), 
the crop categories to be included in the analysis (see @kyle2011gcam for a detailed mapping of GCAM crop categories), the coefficients or counterfactual values for the exposure-response functions (both for health and agricultural damages),
the base Value of Statistical Life (VSL) or Value of Statistical Life Year (VSLY), and additional ancillary information.

The outputs generated by the package consist of both comma-separated values (CSV) files and maps (as Portable Network Graphic files) that can be controlled by the user. If the parameter `saveOutput` is set to `TRUE`, the function writes a CSV file with the selected outcome in the corresponding sub-directory. 
In addition, if `map` is set to `TRUE`, the function generates a suite of maps and animations for the corresponding output. We note that these maps are generated using the [rmap](https://github.com/JGCRI/rmap) package, documented in the following [website](https://jgcri.github.io/rmap).
As an example, the following Figure 2 shows the average $PM_{2.5}$ concentration levels per region, for a GCAM-v5.3 reference scenario.


![$PM_{2.5}$ concentration per country and period in a reference scenario (ug/m3)](figure_conc.png){ height=200% }


Finally, the package is continually being developed to address science objectives and some additional features are scheduled for future releases. For example, an alternative dynamic GDP-based downscaling technique
for re-scaling GCAM emissions in Module 1 (as developed in @gidden2019global), additional age-specific functions for the health impact assessment, as well as a more flexible structure, to allow users to be able to read in emission pathways from different models. 


# Acknowledgements
The research described in this paper was conducted under the Laboratory Directed Research and Development Program at Pacific Northwest National Laboratory, a multiprogram national laboratory operated by Battelle for the U.S. Department of Energy. 
The views and opinions expressed in this paper are those of the authors alone.

# References
---
title: "Module3_health"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Module3_health}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The functions in this module estimate adverse health effects attributable to fine particulate matter (PM2.5) and ozone (O3; M6M) exposure, measured as premature mortalities, years of life lost (YLLs), and disability adjusted life years (DALYs). The following code shows as an example, the functions to estimate premature mortalities attributable to PM2.5 and the DALYs associated with ozone exposure (M6M).All the functions associated to this module can be found in the [References](https://jgcri.github.io/rfasst/reference/index.html) page. 


```{r setup_mort, eval = F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario",
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the packae, "queries_rfasst.xml")


  # To write the csv files into the output folder:

     m3_get_mort_pm25(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 

     m3_get_daly_o3(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 

  # To save this data by TM5-FASST region in year 2050 as dataframes:

     pm25.mort.2050<-dplyr::bind_rows(m3_get_mort_pm25(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F, map=F)) %>% dplyr::filter(year==2050)
     head(pm25.mort.2050)

     o3.daly.2050<-dplyr::bind_rows(m3_get_daly_o3(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F, map=F)) %>% dplyr::filter(year==2050)
     head(o3.daly.2050)

```

The calculation of premature mortalities is consistent with the Global Burden of Disease 2016 ([Lim et al](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(12)61766-8/fulltext)). Mortalities are estimated for five causes, namely stroke, ischemic heart disease (IHD), chronic obstructive pulmonary disease (COPD), acute lower respiratory illness diseases (ALRI) and lung cancer (LC). As explained in the TM5-FASST documentation paper ([Van Dingenen et al 2018](https://acp.copernicus.org/articles/18/16173/2018/acp-18-16173-2018-discussion.html)), cause-specific deaths are estimated following the following equation:
$$Mort_{t,r,c,j}=mo_{r,c,j} \cdot ((RR_{c,j}-1)/RR_{c,j})\cdot Pop_{t,r}$$
So premature mortality in period $t$, region $r$, for cause $c$, associated with exposure to pollutant $j$ is calculated as the product between the baseline mortality rate, the change in the RR relative risk of death attributable to a change in population-weighted mean pollutant concentration, and the population exposed.

* Cause-specific baseline mortality rates are taken from the World Health Organization projections [WHO 2013](https://www.who.int/healthinfo/global_burden_disease/cod_2008_sources_methods.pdf). 

* For PM2.5, relative risk are estimated based on the Integrated Exposure-Response functions (ERFs) from [Burnett et al 2014](https://pubmed.ncbi.nlm.nih.gov/24518036/). For O3, relative risk is based on the ERFs from [Jerret et al 2009](https://www.nejm.org/doi/full/10.1056/nejmoa0803894).

* Population exposed is cause-specific. We assume that for all the causes adults are exposed (>30 Years) while for ALRI, children under 5 years. Population fractions are calculated from the from the [SSP database](https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=welcome).

Years of life lost and DALYs are estimated by transforming premature mortalities into either YLLs or DALYs using constant ratios. YLL-to-Mortalities ratios are based on a reduction in statistical life expectancy (months) per 10µg/m³ increase in anthropogenic PM2.5 (TM5-FASST calculations). DALY-to-Mortalities ratios are estimated using the latest available data on mortailities and DALYs from the Institute for Health Metrics and Evaluation ([IHME](http://www.healthdata.org/)). We note that these ratios (both for YLLs and for DALYs) are calculated regionally so they differ across regions.

Economic damages associated to the health impacts are calculated multiplying either the premature mortalities or the years of life lost by the corresponding Value of Statistical Life (VSL) or Value of Statistical Life Year (VSLY). Both the VSL and the VSLY are based on the widely accepted OECD values for 2005. According to the literature, this value ranges between US$1.8 and $4.5 million for VSL. The calculations for all regions are based on the  “unit value transfer approach” which adjusts the VSL according to their GDP and GDP growth rates, as detialed in [Narain and Sall 2016](https://openknowledge.worldbank.org/handle/10986/24440). In this version, we assume the income elasticity of both the VSL and the VSLY are equal to 0.8, with no adjustments across regions. A potential regional adjustment of these income elasticities will be explored in further developments of the package. We note that the central OECD values could be adjusted by the user. Finally, the damages reported by the functions use the median values of VSL and VSLY for simplicity. The reporting of Lower and Upper Bounds can be selected by the user, but requires some code changes. The systematic reporting of LB and UB is planned for the next package development.

The following code shows some examples for the monetized premature deaths associated to exposure to PM2.5 and monetized years of life lost attributable to O3 exposure:

```{r setup_ecoloss, eval = F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-" Name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"Name of the GCAM scenario"
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the packae, "queries_rfasst.xml")

  # To write the csv files into the output folder:

     m3_get_mort_o3_ecoloss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 

     m3_get_yll_o3_Ecoloss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 

  # To save this data by TM5-FASST region in year 2050 as dataframes:

     pm25.mort.ecoloss.2050<-dplyr::bind_rows(m3_get_mort_o3_ecoloss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050)
     head(pm25.mort.ecoloss.2050)

     o3.yll.ecoloss.2050<-dplyr::bind_rows(m3_get_yll_o3_Ecoloss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050)
     head(o3.yll.ecoloss.2050)

```

The complete list of outputs generated by the suite of functions that form this module is listed below:

* `PM25_MORT_[scenario]_[year].csv`
* `PM25_YLL_[scenario]_[year].csv`
* `PM25_DALY_[scenario]_[year].csv`
* `PM25_MORT_ECOLOSS_[scenario]_[year].csv`
* `PM25_YLL_ECOLOSS__[scenario]_[year].csv`
* `O3_MORT_[scenario]_[year].csv`
* `O3_YLL_[scenario]_[year].csv`
* `O3_DALY_[scenario]_[year].csv`
* `O3_MORT_ECOLOSS_[scenario]_[year].csv`
* `O3_YLL_ECOLOSS__[scenario]_[year].csv`

As in Module 2, for all these functions, the package allows to produce different figures and/or animations, generated using the [rmap](https://github.com/JGCRI/rmap) package documented in the following [page](https://jgcri.github.io/rmap/). To generate these maps, the user needs to include the `map=T` parameter, and they will be generated and stored in the corresponding output sub-directory. As an example for this module, the following map shows the premature mortalities attributable to PM2.5 exposure by cause in 2050.

<!-------------------------->
<!-------------------------->
<p align="center" style="font-size:18px;"> *Premature Mortalities attributable to PM2.5 concentration by cause in 2050 (#)* </p>
<p align="center"> <img src="https://raw.githubusercontent.com/JGCRI/rfasst/main/vignettes/vignetteFigs/PM2.5_mort_2050.png"></p>
<!-------------------------->
<!-------------------------->














---
title: "Module1_emissions"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Module1_emissions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
The first function of the package is `m1_emissions_rescale`. It re-scales emissions from GCAM regions to TM5-FASST regions which is the purpose of Module 1.

The function processes emissions outputs from GCAM to make them match with the regional disaggregation of TM5-FASST. Concretely, GCAM emissions are downscaled to country level, using a static downscaling approach based on the pollutant-specific country-level emissions in year 2000 reported in the [RCP database](https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=welcome). This approach is planned to be revised in the future, for example with the implementation of a dynamic GDP-based downscaling procedure linked with the scenario-specific socioeconomic data in GCAM. 

Apart from the regional adjustments, this function makes some additional pollutant-related changes such as the transformation of organic carbon into organic matter.

Inputs:

* NonCO2 emissions by sector: By default, they will be automatically queried from the GCAM-database
* NonCO2 emissions from international aviation and international shipping: By default, queried from the GCAM database.
* Mapping of pollutants that are reported by GCAM to more aggregated categories which are input in TM5-FASST (my_pol)
* Mapping of GCAM regions to countries and ISO3 codes for downscaling GCAM results to country-level (GCAM_Reg_Adj)
* Mapping of TM5-FASST regions to ISO3 country codes for re-grouping country-level downscaled emissions (fasst_reg)
* Shares to distribute emissions of different species between Russia Eastern (RUE) and Western (RUS) (adj_rus)

Outputs:

The function produces re-scaled emissions of the main pollutants for each period and TM5-FASST region. If `saveOutput` is set to `TRUE`, the function writes the following csv files in the `output/m1` sub-directory: `[scenario]_[year].csv`

In addition, by setting `map` to `TRUE`, the function generates air pollutant emission maps by year and specie, using the [rmap](https://github.com/JGCRI/rmap) package documented in the following [page](https://jgcri.github.io/rmap/). The function also 
generates animations and individual figures for each pollutant by modifying the `mapIndivPol`and `anim` parameters.

The outputs from this module generated with the function `m1_emissions_rescale` are going to read by the set of functions in Module 2, which calculate  PM2.5 and O3 average concentration levels based on these emission sets.


```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup, eval=F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario"
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the package, "queries_rfasst.xml")

   #To write the re-scaled emissions (csv files) for all years into the output folder:
     m1_emissions_rescale(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=T,mapIndivPol=F, anim=T )

   #To save as data frame emisions of main pollutants by TM5-FASST region in 2050:

     em.2050<-dplyr::bind_rows(m1_emissions_rescale(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050) 
                         
     head(m1_emissions_rescale)

```

<!-------------------------->
<!-------------------------->
<p align="center" style="font-size:18px;"> *Air pollutant emissions by specie in 2050 (Gg)* </p>
<p align="center"> <img src="https://raw.githubusercontent.com/JGCRI/rfasst/main/vignettes/vignetteFigs/emissions_2050.png"></p>

<!-------------------------->
<!-------------------------->
---
title: "Module4_agriculture"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Module4_agriculture}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The functions that form this module estimate adverse agricultural impacts attributable to ozone exposure. These functions calculate relative yield losses (RYLs) and the subsequent production and revenue losses damages, combining O3 concentration data (see [Module2 concentration](file:///C:/Users/samp699/Desktop/rfasst/docs/articles/[Module2 concentration].html)), with agricultural production and price projections from GCAM.

RYLs are estimated for four main crops (maize, rice, soybeans, and wheat), using exposure-response functions (ERFs), and are calculated using two different indicators for O3 exposure: accumulated daytime hourly O3 concentration above a threshold of 40 ppbV (AOT40) and the seasonal mean daytime O3 concentration (M7 for the 7-hour mean and M12 for the 12-hour mean). These damges are then extended to the rest of the commodities based on their carbon fixation pathway, as explained below. For AOT40, RYLs are calculated with a linear function described in [Mills et al (2007)](https://www.sciencedirect.com/science/article/pii/S1352231006011356?via%3Dihub):

$$RYL_{t,r,j}=\alpha_j \cdot AOT40_{t,r,j} $$

For Mi, the exposure-response function is detailed in [Wang and Mauzerall (2004)](https://www.sciencedirect.com/science/article/pii/S1352231004004182?via%3Dihub) and follows a Weibull distribution (for $M_i$>$c$, 0 otherwise):

$$RYL_{t,r,j}=1-\frac{exp~[-(\frac{M_{t,r,j}}{a_j})^{b_j}]}{exp~[-(\frac{c_j}{a_j})^{b_j}]}$$

More information on the methodology used for the estimation of agricultural damages can be found in the TM5-FASST documentation paper ([Van Dingenen et al (2018)](https://acp.copernicus.org/articles/18/16173/2018/acp-18-16173-2018-discussion.html)) and in [Van Dingenen et al (2009)](https://www.sciencedirect.com/science/article/pii/S1352231008009424?via%3Dihub). The combined use of the models for estimating agricultural damages is explained in more detail in [Sampedro et al (2020)](https://www.sciencedirect.com/science/article/pii/S1352231020302739).

Production losses are calculated combining the RYLs with projected agricultural production levels of the analyzed GCAM scenario. In addition, multiplying the projected price by the production losses (quantities), we estimate revenue losses for period $t$, region $i$, and crop $j$ as:

   $$Damage_{t,r,j}=Prod_{t,r,j} \cdot Price_{t,r,j} \cdot Pulse.RYL_{t,r,j}$$ 
As explained above, the module combines O3 concentration data with agricultural production and price projections. However, in order to complete all the calculations, the package includes some additional information:

* 2010 data on harvested area for downscalling the damages to country level in order to re-scale them to the corresponding regional disaggregation  (`d.ha`;`area_harvest.csv`)
* RYLs based on the different exposure-response functions are calculated for four main crops (maize, rice, soybeans, and wheat), so the damages need to be expanded to the rest of the commodities. The package includes a commodity mapping that expands the damages to all the crops based on their carbon fixation pathway (C3 or C4 categorization) (`d.gcam.commod.o3`; `GCAM_commod_map_o3.csv`)

We note that these two files could be easily modified by the user if any other mapping/downscalling technique is preferred. The files are stored in the `/inst/extada/mapping` folder.

Finally, all the functions that form this module provide the following list of outputs:

* `RYL_AOT40_[scenario]_[year].csv`
* `RYL_Mi_[scenario]_[year].csv`
* `PROD_LOSS_[scenario]_[year].csv`-- These files include production losses using both the AOT40 and Mi indicators
* `REV_LOSS_[scenario]_[year].csv`-- These files include production losses using both the AOT40 and Mi indicators

The code below shows some examples that show how to produce some of these outputs:

```{r setup, eval = F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario"
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the package, "queries_rfasst.xml")


  # To write the some outputs as csv files into the corresponding output folder:

    # Relative yield losses.
       Using AOT40 as O3 exposure indicator: m4_get_ryl_aot40(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 
       Using Mi as O3 exposure indicator: m4_get_ryl_mi(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 

    # Production losses (includes losses using both AOT40 and Mi).
       m4_get_prod_loss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 

    # Revenue losses (includes losses using both AOT40 and Mi).
       m4_get_rev_loss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F) 

#----------------------------------------------

  # To save this data by TM5-FASST region in year 2050 as dataframes:
    
    # Relative yield losses.
      # Using AOT40 as O3 exposure indicator:
         ryl.aot40.2050<-dplyr::bind_rows(m4_get_ryl_aot40(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F, map=F)) %>% dplyr::filter(year==2050)
                                                              
      # Using Mi as O3 exposure indicator:
         ryl.mi.2050<-dplyr::bind_rows(m4_get_ryl_mi(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F, map=F)) %>% dplyr::filter(year==2050)

    # Production losses (includes losses using both AOT40 and Mi).
       prod.loss.2050<-dplyr::bind_rows(m4_get_prod_loss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F)) %>% dplyr::filter(year==2050)

    # Revenue losses (includes losses using both AOT40 and Mi).
       rev.loss.2050<-dplyr::bind_rows(m4_get_rev_loss(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T, map=F)) %>% dplyr::filter(year==2050)


```

As in Modules 2 and 3, for all these functions, the package allows to produce different figures and/or animations, generated using the [rmap](https://github.com/JGCRI/rmap) package documented in the following [page](jgcri.github.io/rmap/). To generate these maps, the user needs to include the `map=T` parameter, and they will be generated and stored in the corresponding output sub-directory. As an example for this module, the following map shows the  production losses attributable to O3 exposure in 2050.

<!-------------------------->
<!-------------------------->
<p align="center" style="font-size:18px;"> *Agricultural production losses attributable to O3 exposure by commodity in 2050 (Mt)* </p>
<p align="center"> <img src="https://raw.githubusercontent.com/JGCRI/rfasst/main/vignettes/vignetteFigs/Prod_Loss_2050.png"></p>
<!-------------------------->
<!-------------------------->





---
title: "Module2_concentration"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Module2_concentration}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
The functions that form this module take the emission sets per TM5-FASST region and year generated in Module 1 (by `m1_emissions_rescale`), and estimate the fine particulate matter (PM2.5) and ozone (O3) concentration levels for each period and TM5-FASST region. In particular the package reports the following range of PM2.5 and O3 indicators:

Fine particulate matter (PM2.5) concentration levels (`m2_get_conc_pm25`)

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup_pm5, eval = F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario"
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the packae, "queries_rfasst.xml")
 saveOutput: Writes the files.By default=T
 map: Produce the maps. By default=F

  # To write the csv files into the output folder:

     m2_get_conc_pm25(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T,map=F) 

  # To save as data frame pm2.5 concentration levels by TM5-FASST region in year 2050:

     pm25.conc.2050<-dplyr::bind_rows(m2_get_conc_pm25(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050) 
                         
     head(pm25.conc.2050)

```

* Ozone concentration levels (`m2_get_conc_o3`)

```{r setup_o3, eval =F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario",
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the packae, "queries_rfasst.xml")
 saveOutput: Writes the files.By default=T
 map: Produce the maps. By default=F

  # To write the csv files into the output folder:

     m2_get_conc_o3(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T,map=F) 

  # To save as data frame ozone (O3) concentration levels by TM5-FASST region in year 2050:

     o3.conc.2050<-dplyr::bind_rows(m2_get_conc_o3(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050) 
                         
     head(o3.conc.2050)

```

* Maximum 6-monthly running average of daily maximum hourly O3 (M6M, also known as 6mDMA1) (`m2_get_conc_m6m`)

```{r setup_m6m, eval = F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario",
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the packae, "queries_rfasst.xml")
 saveOutput: Writes the files.By default=T
 map: Produce the maps. By default=F

  # To write the csv files into the output folder:

     m2_get_conc_m6m(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T,map=F) 

  # To save as data frame ozone-M6M (6mDMAI) concentration levels by TM5-FASST region in year 2050:

     m6m.conc.2050<-dplyr::bind_rows(m2_get_conc_m6m(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050) 
                         
     head(m6m.conc.2050)

```

* Accumulated daytime hourly O3 concentration above a threshold of 40 ppbV (AOT40) (`m2_get_conc_aot40`)

```{r setup_aot40, eval = F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario"
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the packae, "queries_rfasst.xml")
 saveOutput: Writes the files.By default=T
 map: Produce the maps. By default=F

  # To write the csv files into the output folder:

     m2_get_conc_aot40(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T,map=F) 

  # To save as data frame ozone-AOT40 concentration levels by TM5-FASST region in year 2050:

     aot40.conc.2050<-dplyr::bind_rows(m2_get_conc_aot40(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050) 
                         
     head(aot40.conc.2050)

```

* Seasonal mean daytime O3 concentration (`m2_get_conc_mi`)

```{r setup_mi, eval = F}
library(rfasst)
library(magrittr)

 db_path<-"path_to_your_gcam_database"
 query_path<-"path_to_your_gcam_queries_file"
 db_name<-"name of the database"
 prj_name<-"name for a Project to add extracted results to" # (any name should work, avoid spaces just in case) 
 scen_name<-"name of the GCAM scenario",
 queries<-"Name of the query file" # (the package includes a default query file that includes all the queries required in every function in the packae, "queries_rfasst.xml")
 saveOutput: Writes the files.By default=T
 map: Produce the maps. By default=F

  # To write the csv files into the output folder:

     m2_get_conc_mi(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=T,map=F) 

  # To save as data frame ozone-Mi concentration levels by TM5-FASST region in year 2050:

     mi.conc.2050<-dplyr::bind_rows(m2_get_conc_mi(db_path,query_path,db_name,prj_name,scen_name,queries,saveOutput=F)) %>% dplyr::filter(year==2050) 
                         
     head(mi.conc.2050)

```

In addition, for all these functions, the package allows to produce different figures and/or animations, generated using the [rmap](https://github.com/JGCRI/rmap) package documented in the following [page](jgcri.github.io/rmap/). To generate these maps, the user needs to include the `map=T` parameter, and they will be generated and stored in the corresponding output sub-directory. As an example for this module, the following map shows the PM2.5 concentrations in 2050:
<!-------------------------->
<!-------------------------->
<p align="center" style="font-size:18px;"> *PM2.5 concntration by region in 2050 (ug/m3)* </p>
<p align="center"> <img src="https://raw.githubusercontent.com/JGCRI/rfasst/main/vignettes/vignetteFigs/PM2.5_concentration_2050.png"></p>
<!-------------------------->
<!-------------------------->

As indicated in the documentation of TM5-FASST [Van Dingenen et al 2018](https://acp.copernicus.org/articles/18/16173/2018/acp-18-16173-2018-discussion.html),estimates of PM2.5 and O3 concentration levels in a receptor region driven by the emissions of different precursors in different sources are based on parametrizations of meteorology and atmospheric chemistry drawn from the more complex TM5 model. In summary, concentration of a pollutant $j$, in region $y$, from all the precursors ($i$) emitted in all regions ($x_k$), is calculated as: 

$$C_j (y)=C_{j,base} (y)+∑_{k=1}^{n_x}∑_{i=1}^{n_i}SRC_{i,j} [x_k,y]\cdot[E_i (x_k)-E_{i,base} (x_k )]$$	

 * $C_{j,base} (y)$  is the base-run concentration level of pollutant j in region y, pre-computed with TM5.
 
 * $E_{i,base} (x_k)$ are the base-run precursor i emissions in region $x_k$ 
 
 * $SRC_{i,j} [x_k,y]$ = the $i$–to-$j$ source-receptor coefficient for source region $x_k$ and receptor region y, pre-computed from a 20% emission reduction of component i  in region $x_k$ relative to the base-run
 
 * $E_i (x_k)$ are the emissions of precursor i in region $x_k$ in the analyzed scenario.
 
Following this equation, base-run emissions and concentrations, and source-receptor coefficient matrixes need to be combined with emissions pathways for the analyzed scenario. The package includes the following input information:

* Source-receptor coefficient matrixes (SRC):
  + PM2.5 Source-receptor matrixes:
    - SO4: SO2, NOx and NH3
    - NO3: NOx, SO2, and NH3
    - NH4: NH3, NOx and SO2
    - BC: BC_POP
    - POM: POM_POP
    
  + O3 Source-receptor matrixes:
    - For O3 exposure (O3): NOx, NMVOC and SO2
    - For health calculations: M6M (6mDMA1): NOx, NMVOC, SO2 and CH4
    - For agricultural damages:
      + AOT40: NOx, NMVOC, SO2 and CH4
      + Mi (M7 and M12): NOx, NMVOC, SO2 and CH4

The package also includes the base emission and concentration levels for all these indicators.

In addition, primary PM2.5 emissions (BC and POM) are assumed to have a more direct influence in urban (more dense) areas, so the emission-concentration relation for these two pollutants is modified using adjustment coefficients that are included in the package.



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_aot40_nmvoc}
\alias{src.maize_aot40_nmvoc}
\title{src.maize_aot40_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_aot40_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and AOT40 (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_aot40_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.AOT_SOY}
\alias{coef.AOT_SOY}
\title{coef.AOT_SOY}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Mills, G., Buse, A., Gimeno, B., Bermejo, V., Holland, M., Emberson, L. and Pleijel, H., 2007. A synthesis of AOT40-based response functions and critical levels of ozone for agricultural and horticultural crops. Atmospheric Environment, 41(12), pp.2630-2643.
\dontrun{
 library(rfasst);
 rfasst::coef.AOT_SOY
}
}
\usage{
coef.AOT_SOY
}
\description{
Coefficient for AOT40-Soy
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.o3_nmvoc}
\alias{src.o3_nmvoc}
\title{src.o3_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.o3_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and O3 concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.o3_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.bc}
\alias{src.bc}
\title{SRC BC}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.bc
}
\description{
Source-receptor coefficients (SRC) between BC emissions and BC concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.bc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.base_em}
\alias{raw.base_em}
\title{Base emissions}
\format{
.csv
}
\source{
RCP database: https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=welcome
}
\usage{
raw.base_em
}
\description{
Emissions in the base year (2000) per TM5-FASST region of the main precursors for Particulate Matter (PM2.5), and Ozone (O3 and M6M) formation. These include Black Carbon (BC), Carbon Dioxide (CO2), Methane (CH4), Nitrogen Dioxide (N2O), Organic Matter (POM), Nitrogen Oxides (NOx), sulphur dioxide (SO2), ammonia (NH3) and non-methane volatile organic compounds (NMVOC or VOC)
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.base_em
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.urb_incr}
\alias{raw.urb_incr}
\title{Urban Increment}
\format{
.csv
}
\source{
Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\usage{
raw.urb_incr
}
\description{
Adjustment factor for primary PM2.5 concentrations (BC and POM), which are assumed to be concentrated in urban areas.
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.urb_incr
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m4_agriculture.R
\name{m4_get_ryl_aot40}
\alias{m4_get_ryl_aot40}
\title{m4_get_ryl_aot40}
\usage{
m4_get_ryl_aot40(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
RYLs for each TM5-FASST regions for all years (%). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce Relative Yield Losses (RYLs) based on the AOT40 indicator for O3 exposure
}
\keyword{AOT40}
\keyword{RYLS,}
\keyword{agriculture,}
\keyword{module_4,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{d.ha}
\alias{d.ha}
\title{d.ha}
\format{
An object of class \code{data.frame} with 1109 rows and 4 columns.
}
\source{
GFDL-NOAA
\dontrun{
 library(rfasst);
 rfasst::d.ha
}
}
\usage{
d.ha
}
\description{
Harvested area by crop for weights to map O3 crops to GCAM commodities
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gdp_eu_2005}
\alias{gdp_eu_2005}
\title{gdp_eu_2005}
\format{
An object of class \code{numeric} of length 1.
}
\source{
OECD
\dontrun{
 library(rfasst);
 rfasst::gdp_eu_2005
}
}
\usage{
gdp_eu_2005
}
\description{
Base GDP for EU in 2005
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_mi_so2}
\alias{src.wheat_mi_so2}
\title{src.wheat_mi_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_mi_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and Mi (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_mi_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{adj_rus}
\alias{adj_rus}
\title{adj_rus}
\format{
An object of class \code{data.frame} with 20 rows and 3 columns.
}
\source{
TM5-FASST
}
\usage{
adj_rus
}
\description{
Shares to distribute emissions of different species between Russia Eastern (RUE) and Western (RUS)
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::adj_rus
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{d.gcam.commod.o3}
\alias{d.gcam.commod.o3}
\title{d.gcam.commod.o3}
\format{
An object of class \code{data.frame} with 28 rows and 2 columns.
}
\source{
Own assumptions
\dontrun{
 library(rfasst);
 rfasst::d.gcam.commod.o3
}
}
\usage{
d.gcam.commod.o3
}
\description{
O3 to GCAM commodities (based on their carbon fixation pathways; C3 and C4 categories)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m2_conc.R
\name{m2_get_conc_aot40}
\alias{m2_get_conc_aot40}
\title{m2_get_conc_aot40}
\usage{
m2_get_conc_aot40(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Produce AOT40 levels for each TM5-FASST regions for all years (ppm.h). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce AOT40 concentration levels based on re-scaled emission pathways. AOT40 is the accumulated daytime hourly O3 concentration above a threshold of 40 ppbV (AOT40)
}
\keyword{AOT40}
\keyword{concentration,}
\keyword{module_2,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_aot40_so2}
\alias{src.maize_aot40_so2}
\title{src.maize_aot40_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_aot40_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and AOT40 (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_aot40_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONV_1975_2010}
\alias{CONV_1975_2010}
\title{CONV_1975_2010}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
CONV_1975_2010
}
\description{
1975-2010 deflator
\dontrun{
 library(rfasst);
 rfasst::CONV_1975_2010
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GCAM_reg}
\alias{GCAM_reg}
\title{GCAM_reg}
\format{
An object of class \code{data.frame} with 234 rows and 3 columns.
}
\source{
GCAM
}
\usage{
GCAM_reg
}
\description{
Mapping of countries to GCAM regions
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::GCAM_reg
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.gdp}
\alias{raw.gdp}
\title{GDP-SSP database}
\format{
.csv
}
\source{
https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=welcome
}
\usage{
raw.gdp
}
\description{
Filtered GDP data per SSP (SSP_database_v9.csv). To be consistent we make use of the IIASA Model/scenarios
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.gdp
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.ssp.data}
\alias{raw.ssp.data}
\title{SSP data}
\format{
.csv
}
\source{
https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=welcome
}
\usage{
raw.ssp.data
}
\description{
Country level population and GDP data per SSP (SSP_database_v9.csv). To be consistent we make use of the IIASA-WIC Model/scenarios
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.ssp.data
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_daly_o3}
\alias{m3_get_daly_o3}
\title{m3_get_daly_o3}
\source{
Institute for Health Metrics and Evaluation (http://www.healthdata.org/)
}
\usage{
m3_get_daly_o3(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Disability Adjusted Life Years (DALYs) attributable to O3 exposure for each TM5-FASST regions for all years (# DALYs). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce Disability Adjusted Life Years (DALYs) attributable to O3 (M6M) exposure. See calc_daly_o3 for detials on DALY-to-Mortality ratios.
}
\keyword{DALY,}
\keyword{O3,}
\keyword{module_3,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_mi_nox}
\alias{src.wheat_mi_nox}
\title{src.wheat_mi_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_mi_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and Mi (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_mi_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_daly_pm25}
\alias{m3_get_daly_pm25}
\title{m3_get_daly_pm25}
\source{
Institute for Health Metrics and Evaluation (http://www.healthdata.org/)
}
\usage{
m3_get_daly_pm25(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Disability Adjusted Life Years (DALYs) attributable to PM2.5 exposure for each TM5-FASST regions for all years (# DALYs). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce Disability Adjusted Life Years (DALYs) attributable to PM2.5 exposure. See calc_daly_pm for detials on DALY-to-Mortality ratios.
}
\keyword{DALY,}
\keyword{PM2.5,}
\keyword{module_3,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_aot40_nox}
\alias{src.maize_aot40_nox}
\title{src.maize_aot40_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_aot40_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and AOT40 (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_aot40_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m4_agriculture.R
\name{m4_get_rev_loss}
\alias{m4_get_rev_loss}
\title{m4_get_rev_loss}
\usage{
m4_get_rev_loss(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Revenue losses attributable to O3 for each GCAM region for all years (Billion$2010). The countries that form each GCAM region are listed in the documentation repository: https://github.com/JGCRI/gcam-doc/blob/gh-pages/overview.md The list of commodities within each category can be found in: Kyle, G.P., Luckow, P., Calvin, K.V., Emanuel, W.R., Nathan, M. and Zhou, Y., 2011. GCAM 3.0 agriculture and land use: data sources and methods (No. PNNL-21025). Pacific Northwest National Lab.(PNNL), Richland, WA (United States), and in Sampedro, J., Waldhoff, S.T., Van de Ven, D.J., Pardo, G., Van Dingenen, R., Arto, I., del Prado, A. and Sanz, M.J., 2020. Future impacts of ozone driven damages on agricultural systems. Atmospheric Environment, 231, p.117538, Table S2.
}
\description{
Produce agricultural revenue losses attributable to ozone exposure for all GCAM crops. Losses have been calculated using two ozone exposure indicators: AOT40 and Mi.
}
\keyword{O3,}
\keyword{agriculture,}
\keyword{losses}
\keyword{module_4,}
\keyword{revenue}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m4_agriculture.R
\name{m4_get_ryl_mi}
\alias{m4_get_ryl_mi}
\title{m4_get_ryl_mi}
\usage{
m4_get_ryl_mi(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
RYLs for each TM5-FASST regions for all years (%).  The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce Relative Yield Losses (RYLs) based on the Mi indicator for O3 exposure
}
\keyword{Mi}
\keyword{RYLS,}
\keyword{agriculture,}
\keyword{module_4,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rfasst.R
\docType{package}
\name{rfasst}
\alias{rfasst}
\title{rfasst}
\description{
\itemize{
\item Github: https://github.com/JGCRI/rfasst
\item Webpage: https://jgcri.github.io/rfasst/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_mort_o3}
\alias{m3_get_mort_o3}
\title{m3_get_mort_o3}
\source{
Jerrett, M., Burnett, R.T., Pope III, C.A., Ito, K., Thurston, G., Krewski, D., Shi, Y., Calle, E. and Thun, M., 2009. Long-term ozone exposure and mortality. New England Journal of Medicine, 360(11), pp.1085-1095.
}
\usage{
m3_get_mort_o3(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Premature mortality attributable to O3 exposure for  TM5-FASST regions for all years (# mortalties). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce premature mortality attributable to O3 exposure (measured by the M6M indicator) based on the integrated exposure-response functions (IER) from Jerret et al (2009), consistent with the GBD 2016 study.
}
\keyword{O3}
\keyword{module_3,}
\keyword{mortality,}
\keyword{premature}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_mi_nmvoc}
\alias{src.soy_mi_nmvoc}
\title{src.soy_mi_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_mi_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and Mi (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_mi_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_mi_nmvoc}
\alias{src.maize_mi_nmvoc}
\title{src.maize_mi_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_mi_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and Mi (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_mi_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{b.Mi_SOY}
\alias{b.Mi_SOY}
\title{b.Mi_SOY}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::b.Mi_SOY
}
}
\usage{
b.Mi_SOY
}
\description{
"b" coefficient for Mi-Soy
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.Mi_SOY}
\alias{coef.Mi_SOY}
\title{coef.Mi_SOY}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::coef.Mi_SOY
}
}
\usage{
coef.Mi_SOY
}
\description{
Coefficient for Mi-Soy
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.AOT_WHEAT}
\alias{coef.AOT_WHEAT}
\title{coef.AOT_WHEAT}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Mills, G., Buse, A., Gimeno, B., Bermejo, V., Holland, M., Emberson, L. and Pleijel, H., 2007. A synthesis of AOT40-based response functions and critical levels of ozone for agricultural and horticultural crops. Atmospheric Environment, 41(12), pp.2630-2643.
\dontrun{
 library(rfasst);
 rfasst::coef.AOT_WHEAT
}
}
\usage{
coef.AOT_WHEAT
}
\description{
Coefficient for AOT40-Wheat
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_mort_pm25_ecoloss}
\alias{m3_get_mort_pm25_ecoloss}
\title{m3_get_mort_pm25_ecoloss}
\source{
Narain, U. and Sall, C., 2016. Methodology for Valuing the Health Impacts of Air Pollution//// Markandya, A., Sampedro, J., Smith, S.J., Van Dingenen, R., Pizarro-Irizar, C., Arto, I. and González-Eguino, M., 2018. Health co-benefits from air pollution and mitigation costs of the Paris Agreement: a modelling study. The Lancet Planetary Health, 2(3), pp.e126-e133.
}
\usage{
m3_get_mort_pm25_ecoloss(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  ssp = "SSP2",
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Economic damages associated with mortality attributable to PM2.5 exposure for each TM5-FASST regions for all years (Million$2015). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce economic damages associated with premature mortality attributable to PM2.5 exposure based on the IER functions from Burnett et al (2014), consistent with the GBD 2016 study. The economic valuation takes as a base value the widely accepted Value of Statistical Life (VSL) of the OECD for 2005. This value, according to the literature ranges between US$1.8 and $4.5 million. The calculations for all regions are based on the  “unit value transfer approach” which adjusts the VSL according to their GDP and GDP growth rates. (Markandya et al 2018)
}
\keyword{,premature}
\keyword{PM2.5}
\keyword{VSL}
\keyword{module_3,}
\keyword{mortality,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.no3_nh3}
\alias{src.no3_nh3}
\title{no3_nh3}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.no3_nh3
}
\description{
Source-receptor coefficients (SRC) between NH3 emissions and NO3 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.no3_nh3
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.m6m_nox}
\alias{src.m6m_nox}
\title{src.m6m_nox}
\format{
.csv

.csv
}
\source{
Results from TM5
}
\usage{
src.m6m_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and M6M (O3) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.m6m_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_aot40_ch4}
\alias{src.wheat_aot40_ch4}
\title{src.wheat_aot40_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_aot40_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and AOT40 (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_aot40_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.mort.rates}
\alias{raw.mort.rates}
\title{raw.mort.rates}
\format{
.csv
}
\source{
https://www.who.int/healthinfo/global_burden_disease/cod_2008_sources_methods.pdf
}
\usage{
raw.mort.rates
}
\description{
cause-specific baseline mortalities from stroke, ischemic heart disease (IHD), chronic obstructive pulmonary disease (COPD), acute lower respiratory illness diseases (ALRI) and lung cancer (LC).
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.mort.rates
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{inc_elas_vsl}
\alias{inc_elas_vsl}
\title{inc_elas_vsl}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Own assumptions
\dontrun{
 library(rfasst);
 rfasst::inc_elas_vsl
}
}
\usage{
inc_elas_vsl
}
\description{
Income elasticity for the Value of Statistical Life (VSL)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_aot40_so2}
\alias{src.rice_aot40_so2}
\title{src.rice_aot40_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_aot40_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and AOT40 (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_aot40_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vsl_eu_2005_lb}
\alias{vsl_eu_2005_lb}
\title{vsl_eu_2005_lb}
\format{
An object of class \code{numeric} of length 1.
}
\source{
OECD
\dontrun{
 library(rfasst);
 rfasst::vsl_eu_2005_lb
}
}
\usage{
vsl_eu_2005_lb
}
\description{
Lower bound for the Value of Statistical Life (VSL)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{country_iso}
\alias{country_iso}
\title{country_iso}
\format{
An object of class \code{data.frame} with 249 rows and 2 columns.
}
\usage{
country_iso
}
\description{
Mapping of countries to iso3 codes
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::country_iso
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{perc_pop_rue}
\alias{perc_pop_rue}
\title{perc_pop_rue}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
perc_pop_rue
}
\description{
Percentages to divide population between Russia and Russia Eastern
\dontrun{
 library(rfasst);
 rfasst::perc_pop_rue
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.nh4_nh3}
\alias{src.nh4_nh3}
\title{nh4_nh3}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.nh4_nh3
}
\description{
Source-receptor coefficients (SRC) between NH3 emissions and NH4 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.nh4_nh3
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{a.Mi_WHEAT}
\alias{a.Mi_WHEAT}
\title{a.Mi_WHEAT}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::a.Mi_WHEAT
}
}
\usage{
a.Mi_WHEAT
}
\description{
"a" coefficient for Mi-Wheat
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_mi_nox}
\alias{src.maize_mi_nox}
\title{src.maize_mi_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_mi_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and Mi (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_mi_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.Mi_WHEAT}
\alias{coef.Mi_WHEAT}
\title{coef.Mi_WHEAT}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::coef.Mi_WHEAT
}
}
\usage{
coef.Mi_WHEAT
}
\description{
Coefficient for Mi-Wheat
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m4_agriculture.R
\name{calc_rev_gcam}
\alias{calc_rev_gcam}
\title{calc_rev_gcam}
\source{
GCAM
}
\usage{
calc_rev_gcam(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Future agricultural revenue levels from GCAM for each GCAM region for all years (Billion$2010). The countries that form each GCAM region are listed in the documentation repository: https://github.com/JGCRI/gcam-doc/blob/gh-pages/overview.md The list of commodities within each category can be found in: Kyle, G.P., Luckow, P., Calvin, K.V., Emanuel, W.R., Nathan, M. and Zhou, Y., 2011. GCAM 3.0 agriculture and land use: data sources and methods (No. PNNL-21025). Pacific Northwest National Lab.(PNNL), Richland, WA (United States), and in Sampedro, J., Waldhoff, S.T., Van de Ven, D.J., Pardo, G., Van Dingenen, R., Arto, I., del Prado, A. and Sanz, M.J., 2020. Future impacts of ozone driven damages on agricultural systems. Atmospheric Environment, 231, p.117538, Table S2.
}
\description{
Extract future agricultural revenues from GCAM
}
\keyword{agriculture,}
\keyword{module_4,}
\keyword{revenues}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.o3_ch4}
\alias{src.o3_ch4}
\title{src.o3_ch4}
\format{
.csv
}
\source{
Fiore, A.M., West, J.J., Horowitz, L.W., Naik, V. and Schwarzkopf, M.D., 2008. Characterizing the tropospheric ozone response to methane emission controls and the benefits to climate and air quality. Journal of Geophysical Research: Atmospheres, 113(D8).
}
\usage{
src.o3_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and O3 concentration level, normalized from existing literature (not pre-computed with TM5-FASST)
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.o3_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_mi_so2}
\alias{src.soy_mi_so2}
\title{src.soy_mi_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_mi_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and Mi (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_mi_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.o3_nox}
\alias{src.o3_nox}
\title{src.o3_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.o3_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and O3 concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.o3_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vsl_eu_2005_ub}
\alias{vsl_eu_2005_ub}
\title{vsl_eu_2005_ub}
\format{
An object of class \code{numeric} of length 1.
}
\source{
OECD
\dontrun{
 library(rfasst);
 rfasst::vsl_eu_2005_lb
}
}
\usage{
vsl_eu_2005_ub
}
\description{
Upper bound for the Value of Statistical Life (VSL)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m4_agriculture.R
\name{m4_get_prod_loss}
\alias{m4_get_prod_loss}
\title{m4_get_prod_loss}
\usage{
m4_get_prod_loss(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Ag losses attributable to O3 for each GCAM region for all years (Mt).The countries that form each GCAM region are listed in the documentation repository: https://github.com/JGCRI/gcam-doc/blob/gh-pages/overview.md The list of commodities within each category can be found in: Kyle, G.P., Luckow, P., Calvin, K.V., Emanuel, W.R., Nathan, M. and Zhou, Y., 2011. GCAM 3.0 agriculture and land use: data sources and methods (No. PNNL-21025). Pacific Northwest National Lab.(PNNL), Richland, WA (United States), and in Sampedro, J., Waldhoff, S.T., Van de Ven, D.J., Pardo, G., Van Dingenen, R., Arto, I., del Prado, A. and Sanz, M.J., 2020. Future impacts of ozone driven damages on agricultural systems. Atmospheric Environment, 231, p.117538, Table S2.
}
\description{
Produce agricultural production losses attributable to ozone exposure for all GCAM crops. Losses have been calculated using two ozone exposure indicators: AOT40 and Mi.
}
\keyword{O3,production}
\keyword{agriculture,}
\keyword{losses}
\keyword{module_4,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONV_OC_POM}
\alias{CONV_OC_POM}
\title{CONV_OC_POM}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Kanakidou, M., Seinfeld, J.H., Pandis, S.N., Barnes, I., Dentener, F.J., Facchini, M.C., Dingenen, R.V., Ervens, B., Nenes, A.N.C.J.S.E., Nielsen, C.J. and Swietlicki, E., 2005. Organic aerosol and global climate modelling: a review. Atmospheric Chemistry and Physics, 5(4), pp.1053-1123.
\dontrun{
 library(rfasst);
 rfasst::CONV_OC_POM
}
}
\usage{
CONV_OC_POM
}
\description{
Transform OC to POM
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{all_years}
\alias{all_years}
\title{all_years}
\format{
An object of class \code{character} of length 11.
}
\usage{
all_years
}
\description{
Years to be analyzed: c('2005','2010','2020','2030','2040','2050','2060','2070','2080','2090','2100')
\dontrun{
 library(rfasst);
 rfasst::all_years
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_mi_so2}
\alias{src.maize_mi_so2}
\title{src.maize_mi_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_mi_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and Mi (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_mi_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ancillary_functions.R
\name{calc_gdp_pc}
\alias{calc_gdp_pc}
\title{calc_gdp_pc}
\source{
https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=welcome
}
\usage{
calc_gdp_pc(ssp = "SSP2")
}
\arguments{
\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}
}
\value{
GDP_pc for TM5-FASST regions for all years
}
\description{
Get GDP_pc from the SSP database (SSP_database_v9) for the economic assessment of the health impacts. To be consistent we make use of the IIASA Model/scenarios.
}
\keyword{GDP}
\keyword{socioeconomics,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Percen}
\alias{Percen}
\title{#' Percen}
\format{
An object of class \code{data.frame} with 17490 rows and 6 columns.
}
\source{
RCP database: https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=welcome
\dontrun{
 library(rfasst);
 rfasst::Percen
}
}
\usage{
Percen
}
\description{
Percentages to downscale GCAM emissions to country-level
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.m6m_ch4}
\alias{src.m6m_ch4}
\title{src.m6m_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.m6m_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and M6M (O3) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.m6m_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.twn.pop}
\alias{raw.twn.pop}
\title{Socioeconomics Taiwan}
\format{
.csv
}
\source{
https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=welcome
}
\usage{
raw.twn.pop
}
\description{
Population and GDP data per SSP (SSP_database_v9.csv). Given that IIASA-WIC does not report data for Taiwan, we use data from "OECD_Env-Growth"
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.twn.pop
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_aot40_ch4}
\alias{src.maize_aot40_ch4}
\title{src.maize_aot40_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_aot40_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and AOT40 (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_aot40_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.maize_mi_ch4}
\alias{src.maize_mi_ch4}
\title{src.maize_mi_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.maize_mi_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and Mi (O3) concentration level for maize, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.maize_mi_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rr_resp_o3}
\alias{rr_resp_o3}
\title{rr_resp_o3}
\format{
An object of class \code{numeric} of length 1.
}
\source{
TM5-FASST
\dontrun{
 library(rfasst);
 rfasst::rr_resp_o3
}
}
\usage{
rr_resp_o3
}
\description{
Relative risk for respiratory disease associated to ozone exposure
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_aot40_nmvoc}
\alias{src.soy_aot40_nmvoc}
\title{src.soy_aot40_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_aot40_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and AOT40 (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_aot40_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MTC_MTCO2}
\alias{MTC_MTCO2}
\title{MTC_MTCO2}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
MTC_MTCO2
}
\description{
Transform MTC to MTCO2
\dontrun{
 library(rfasst);
 rfasst::MTC_MTCO2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_mi_nmvoc}
\alias{src.wheat_mi_nmvoc}
\title{src.wheat_mi_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_mi_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and Mi (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_mi_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_mi_nmvoc}
\alias{src.rice_mi_nmvoc}
\title{src.rice_mi_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_mi_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and Mi (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_mi_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{d.iso}
\alias{d.iso}
\title{d.iso}
\format{
An object of class \code{data.frame} with 238 rows and 4 columns.
}
\usage{
d.iso
}
\description{
Countries to GCAM regions
\dontrun{
 library(rfasst);
 rfasst::d.iso
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONV_TG_T}
\alias{CONV_TG_T}
\title{CONV_TG_T}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
CONV_TG_T
}
\description{
Unit converter: Teragram to tonne
\dontrun{
 library(rfasst);
 rfasst::CONV_TG_T
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.daly}
\alias{raw.daly}
\title{raw.daly}
\format{
An object of class \code{data.frame} with 85680 rows and 18 columns.
}
\source{
Institute for Health Metrics and Evaluation (http://www.healthdata.org/)
}
\usage{
raw.daly
}
\description{
Data on Disability Adjusted Life Years (DALYs). Used the latest available data (2019)
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.daly
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.pom}
\alias{src.pom}
\title{SRC POM}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.pom
}
\description{
Source-receptor coefficients (SRC) between POM emissions and POM concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.pom
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ch4_htap_pert}
\alias{ch4_htap_pert}
\title{ch4_htap_pert}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Fiore, A.M., West, J.J., Horowitz, L.W., Naik, V. and Schwarzkopf, M.D., 2008. Characterizing the tropospheric ozone response to methane emission controls and the benefits to climate and air quality. Journal of Geophysical Research: Atmospheres, 113(D8).
\dontrun{
 library(rfasst);
 rfasst::ch4_htap_pert
}
}
\usage{
ch4_htap_pert
}
\description{
Normalized CH4-O3 relation from Fiore et al (2008)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m2_conc.R
\name{m2_get_conc_pm25}
\alias{m2_get_conc_pm25}
\title{m2_get_conc_pm25}
\usage{
m2_get_conc_pm25(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Particulate matter (PM2.5) concentration levels for each TM5-FASST regions for all years (ug/m3).The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce fine particulate matter (PM2.5) concentration levels for TM5-FASST regions based on re-scaled emission pathways.
}
\keyword{PM2.5}
\keyword{concentration,}
\keyword{module_2,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vsly_eu_2005}
\alias{vsly_eu_2005}
\title{vsly_eu_2005}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Schlander, M., Schaefer, R. and Schwarz, O., 2017. Empirical studies on the economic value of a Statistical Life Year (VSLY) in Europe: what do they tell us?. Value in Health, 20(9), p.A666.
\dontrun{
 library(rfasst);
 rfasst::vsly_eu_2005
}
}
\usage{
vsly_eu_2005
}
\description{
Base Value of Statistical Life Year (VSLY) in $2005
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{calc_daly_pm25}
\alias{calc_daly_pm25}
\title{calc_daly_pm25}
\source{
Institute for Health Metrics and Evaluation (http://www.healthdata.org/)
}
\usage{
calc_daly_pm25()
}
\value{
DALY-to-Mortality ratios for TM5-FASST regions for all years and PM2.5-related causes (ALRI, COPD, LC, IHD, STROKE).The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Get the DALY-to-Mortality ratios used to estimate the Disability Adjusted Life Years attributable to fine particulate matter (PM2.5) exposure
}
\keyword{DALYs,}
\keyword{PM2.5}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_yll_o3_ecoloss}
\alias{m3_get_yll_o3_ecoloss}
\title{m3_get_yll_o3_ecoloss}
\source{
Schlander, M., Schaefer, R. and Schwarz, O., 2017. Empirical studies on the economic value of a Statistical Life Year (VSLY) in Europe: what do they tell us?. Value in Health, 20(9), p.A666.
}
\usage{
m3_get_yll_o3_ecoloss(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Economic damages associated with YLLs attributable to O3 (M6M) exposure for each TM5-FASST regions for all years (Thous$2015). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce Economic damages associated with YLLs attributable to O3 (M6M).The economic valuation takes as a base value the Value of Statistical Life Year (VSLY) for EU from Schlander et al (2017) and expands the value to other regions based on the“unit value transfer approach” which adjusts the VSLY according to their GDP and GDP growth rates. YLL-to-Mortalities ratios are based on TM5-FASST calculations. Premature mortalities are  based on the integrated exposure-response functions (IER) from Burnett et al (2014), consistent with the GBD 2016 study.
}
\keyword{O3,}
\keyword{VSLY}
\keyword{YLL,}
\keyword{module_3,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.m6m_nmvoc}
\alias{src.m6m_nmvoc}
\title{src.m6m_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.m6m_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and M6M (O3) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.m6m_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{map_pol}
\alias{map_pol}
\title{map_pol}
\format{
An object of class \code{character} of length 6.
}
\usage{
map_pol
}
\description{
Indicate the pollutants whose emissions are mapped (if map=T in m1_emissions_rescale)
\dontrun{
 library(rfasst);
 rfasst::map_pol
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_mi_nox}
\alias{src.soy_mi_nox}
\title{src.soy_mi_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_mi_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and Mi (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_mi_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_aot40_ch4}
\alias{src.rice_aot40_ch4}
\title{src.rice_aot40_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_aot40_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and AOT40 (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_aot40_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{d.weight.gcam}
\alias{d.weight.gcam}
\title{d.weight.gcam}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 872 rows and 4 columns.
}
\source{
Own assumptions
\dontrun{
 library(rfasst);
 rfasst::d.weight.gcam
}
}
\usage{
d.weight.gcam
}
\description{
O3 to GCAM commodities (based on their carbon fixation pathways; C3 and C4 categories)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vsl_eu_2005}
\alias{vsl_eu_2005}
\title{vsl_eu_2005}
\format{
An object of class \code{numeric} of length 1.
}
\source{
OECD
\dontrun{
 library(rfasst);
 rfasst::vsl_eu_2005
}
}
\usage{
vsl_eu_2005
}
\description{
Median bound for the Value of Statistical Life (VSL)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ancillary_functions.R
\name{interpLinear}
\alias{interpLinear}
\title{interpLinear}
\usage{
interpLinear(d, y_start, y_end)
}
\arguments{
\item{d}{Data frame to be interpolated}

\item{y_start}{Starting year (start of the decade)}

\item{y_end}{End year (End of the decade)}
}
\description{
Function to interpolate annual values using decade-averages
}
\keyword{interpolate}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m2_conc.R
\name{m2_get_conc_mi}
\alias{m2_get_conc_mi}
\title{m2_get_conc_mi}
\usage{
m2_get_conc_mi(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Produce Mi levels for each TM5-FASST regions for all years (ppb).The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce Mi concentration levels based on re-scaled emission pathways from module 1. Mi is the the seasonal mean daytime O3 concentration (M7 for the 7-hour mean and M12 for the 12-hour mean)
}
\keyword{(M7}
\keyword{M12)}
\keyword{Mi}
\keyword{and}
\keyword{concentration,}
\keyword{module_2,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_mi_ch4}
\alias{src.soy_mi_ch4}
\title{src.soy_mi_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_mi_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and Mi (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_mi_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ancillary_functions.R
\name{calc_pop}
\alias{calc_pop}
\title{calc_pop}
\source{
https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=welcome
}
\usage{
calc_pop(ssp = "SSP2")
}
\arguments{
\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}
}
\value{
Population and population shares (<5Y; >30Y) for TM5-FASST regions for all years
}
\description{
Get population data and shares of population under 5 Years and above 30 Years from the SSP database (SSP_database_v9).To be consistent we make use of the IIASA-WIC Model/scenarios. Given that IIASA-WIC does not report data for Taiwan, we use data from "OECD_Env-Growth" for this region.
}
\keyword{population}
\keyword{socioeconomics,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.yll.pm25}
\alias{raw.yll.pm25}
\title{raw.yll.pm25}
\format{
An object of class \code{data.frame} with 56 rows and 8 columns.
}
\source{
TM5-FASST
}
\usage{
raw.yll.pm25
}
\description{
Years of Life Lost (YLLs) to Mortality ratios attributable to PM2.5 exposure
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.yll.pm25
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{b.Mi_RICE}
\alias{b.Mi_RICE}
\title{b.Mi_RICE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::b.Mi_RICE
}
}
\usage{
b.Mi_RICE
}
\description{
"b" coefficient for Mi-Rice
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.so4_nh3}
\alias{src.so4_nh3}
\title{so4_nh3}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.so4_nh3
}
\description{
Source-receptor coefficients (SRC) between NH3 emissions and SO4 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.so4_nh3
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CROP_ANALYSIS}
\alias{CROP_ANALYSIS}
\title{CROP_ANALYSIS}
\format{
An object of class \code{character} of length 12.
}
\usage{
CROP_ANALYSIS
}
\description{
Crops that are included in the analysis
\dontrun{
 library(rfasst);
 rfasst::CROP_ANALYSIS
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.m6m_so2}
\alias{src.m6m_so2}
\title{src.m6m_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.m6m_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and M6M (O3) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.m6m_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.nh4_so2}
\alias{src.nh4_so2}
\title{nh4_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.nh4_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and NH4 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.nh4_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_mi_so2}
\alias{src.rice_mi_so2}
\title{src.rice_mi_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_mi_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and Mi (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_mi_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.nh4_nox}
\alias{src.nh4_nox}
\title{nh4_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.nh4_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and NH4 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.nh4_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Regions}
\alias{Regions}
\title{Regions}
\format{
An object of class \code{data.frame} with 249 rows and 4 columns.
}
\usage{
Regions
}
\description{
Combined regions
\dontrun{
 library(rfasst);
 rfasst::Regions
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{my_pol}
\alias{my_pol}
\title{#' my_pol}
\format{
An object of class \code{data.frame} with 42 rows and 2 columns.
}
\source{
GCAM
\dontrun{
 library(rfasst);
 rfasst::my_pol
}
}
\usage{
my_pol
}
\description{
Information about GCAM and TM5-FASST regions and pollutants and their equivalences.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_aot40_nox}
\alias{src.wheat_aot40_nox}
\title{src.wheat_aot40_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_aot40_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and AOT40 (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_aot40_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m2_conc.R
\name{m2_get_conc_m6m}
\alias{m2_get_conc_m6m}
\title{m2_get_conc_m6m}
\usage{
m2_get_conc_m6m(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
M6M levels for each TM5-FASST regions for all years. The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce M6M concentration levels based on re-scaled emission pathways. M6M is maximum 6-monthly running average of daily maximum hourly O3 (ppb) (6mDMA1).
}
\keyword{M6M}
\keyword{concentration,}
\keyword{module_2,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_aot40_nox}
\alias{src.rice_aot40_nox}
\title{src.rice_aot40_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_aot40_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and AOT40 (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_aot40_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_aot40_so2}
\alias{src.wheat_aot40_so2}
\title{src.wheat_aot40_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_aot40_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and AOT40 (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_aot40_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.AOT_MAIZE}
\alias{coef.AOT_MAIZE}
\title{coef.AOT_MAIZE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Mills, G., Buse, A., Gimeno, B., Bermejo, V., Holland, M., Emberson, L. and Pleijel, H., 2007. A synthesis of AOT40-based response functions and critical levels of ozone for agricultural and horticultural crops. Atmospheric Environment, 41(12), pp.2630-2643.
\dontrun{
 library(rfasst);
 rfasst::coef.AOT_MAIZE
}
}
\usage{
coef.AOT_MAIZE
}
\description{
Coefficient for AOT40-Maize
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_mi_nox}
\alias{src.rice_mi_nox}
\title{src.rice_mi_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_mi_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and Mi (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_mi_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.AOT_RICE}
\alias{coef.AOT_RICE}
\title{coef.AOT_RICE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Mills, G., Buse, A., Gimeno, B., Bermejo, V., Holland, M., Emberson, L. and Pleijel, H., 2007. A synthesis of AOT40-based response functions and critical levels of ozone for agricultural and horticultural crops. Atmospheric Environment, 41(12), pp.2630-2643.
\dontrun{
 library(rfasst);
 rfasst::coef.AOT_RICE
}
}
\usage{
coef.AOT_RICE
}
\description{
Coefficient for AOT40-Rice
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONV_OCawb_POM}
\alias{CONV_OCawb_POM}
\title{CONV_OCawb_POM}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Kanakidou, M., Seinfeld, J.H., Pandis, S.N., Barnes, I., Dentener, F.J., Facchini, M.C., Dingenen, R.V., Ervens, B., Nenes, A.N.C.J.S.E., Nielsen, C.J. and Swietlicki, E., 2005. Organic aerosol and global climate modelling: a review. Atmospheric Chemistry and Physics, 5(4), pp.1053-1123.
\dontrun{
 library(rfasst);
 rfasst::CONV_OCawb_POM
}
}
\usage{
CONV_OCawb_POM
}
\description{
Transform OC from biogenic sources to POM
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{a.Mi_RICE}
\alias{a.Mi_RICE}
\title{a.Mi_RICE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::a.Mi_RICE
}
}
\usage{
a.Mi_RICE
}
\description{
"a" coefficient for Mi-Rice
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONV_BIL}
\alias{CONV_BIL}
\title{CONV_BIL}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
CONV_BIL
}
\description{
Transform $Billion to $
\dontrun{
 library(rfasst);
 rfasst::CONV_BIL
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m4_agriculture.R
\name{calc_prod_gcam}
\alias{calc_prod_gcam}
\title{calc_prod_gcam}
\source{
GCAM
}
\usage{
calc_prod_gcam(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Future agricultural production levels from GCAM for each GCAM region for all years (Mt). The countries that form each GCAM region are listed in the documentation repository: https://github.com/JGCRI/gcam-doc/blob/gh-pages/overview.md The list of commodities within each category can be found in: Kyle, G.P., Luckow, P., Calvin, K.V., Emanuel, W.R., Nathan, M. and Zhou, Y., 2011. GCAM 3.0 agriculture and land use: data sources and methods (No. PNNL-21025). Pacific Northwest National Lab.(PNNL), Richland, WA (United States), and in Sampedro, J., Waldhoff, S.T., Van de Ven, D.J., Pardo, G., Van Dingenen, R., Arto, I., del Prado, A. and Sanz, M.J., 2020. Future impacts of ozone driven damages on agricultural systems. Atmospheric Environment, 231, p.117538, Table S2.
}
\description{
Extract future agricultural production levels from GCAM
}
\keyword{agriculture,}
\keyword{module_4,}
\keyword{production}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.so4_nox}
\alias{src.so4_nox}
\title{so4_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.so4_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and SO4 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.so4_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.base_aot}
\alias{raw.base_aot}
\title{Base AOT40 conc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
raw.base_aot
}
\description{
O3 concentration above a threshold of 40 ppbV (AOT40), in the base year (2000) per TM5-FASST region
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.base_aot
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{b.Mi_MAIZE}
\alias{b.Mi_MAIZE}
\title{b.Mi_MAIZE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::b.Mi_MAIZE
}
}
\usage{
b.Mi_MAIZE
}
\description{
"b" coefficient for Mi-Maize
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.Mi_RICE}
\alias{coef.Mi_RICE}
\title{coef.Mi_RICE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::coef.Mi_RICE
}
}
\usage{
coef.Mi_RICE
}
\description{
Coefficient for Mi-Rice
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fasst_reg}
\alias{fasst_reg}
\title{fasst_reg}
\format{
An object of class \code{data.frame} with 249 rows and 2 columns.
}
\source{
TM5-FASST
}
\usage{
fasst_reg
}
\description{
Mapping of countries to TM5-FASST regions
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::fasst_reg
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_yll_o3}
\alias{m3_get_yll_o3}
\title{m3_get_yll_o3}
\usage{
m3_get_yll_o3(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
YLLs attributable to O3 exposure for each TM5-FASST regions for all years (# YLLs). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce YLLs attributable to O3 (M6M) exposure. YLL-to-Mortalities ratios are based on TM5-FASST calculations. Premature mortalities are based on the IER functions from Jerret et al (2009), consistent with the GBD 2016 study.
}
\keyword{O3}
\keyword{YLL,}
\keyword{module_3,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_yll_pm25_ecoloss}
\alias{m3_get_yll_pm25_ecoloss}
\title{m3_get_yll_pm25_ecoloss}
\source{
Schlander, M., Schaefer, R. and Schwarz, O., 2017. Empirical studies on the economic value of a Statistical Life Year (VSLY) in Europe: what do they tell us?. Value in Health, 20(9), p.A666.
}
\usage{
m3_get_yll_pm25_ecoloss(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Economic damages associated with YLLs attributable to PM2.5 exposure for each TM5-FASST regions for all years (Thous$2015).The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce Economic damages associated with YLLs attributable to PM2.5.The economic valuation takes as a base value the Value of Statistical Life Year (VSLY) for EU from Schlander et al (2017) and expands the value to other regions based on the“unit value transfer approach” which adjusts the VSLY according to their GDP and GDP growth rates. . .YLL-to-Mortalities ratios are based on TM5-FASST calculations. Premature mortalities are  based on the integrated exposure-response functions (IER) from Burnett et al (2014), consistent with the GBD 2016 study.
}
\keyword{PM2.5,}
\keyword{VSLY}
\keyword{YLL,}
\keyword{module_3,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_aot40_nmvoc}
\alias{src.rice_aot40_nmvoc}
\title{src.rice_aot40_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_aot40_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and AOT40 (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_aot40_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{perc_pop_rus}
\alias{perc_pop_rus}
\title{perc_pop_rus}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
perc_pop_rus
}
\description{
Percentages to divide population between Russia and Russia Eastern
\dontrun{
 library(rfasst);
 rfasst::perc_pop_rus
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TG_KG}
\alias{TG_KG}
\title{TG_KG}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
TG_KG
}
\description{
Transform Tg to Kg
\dontrun{
 library(rfasst);
 rfasst::TG_KG
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.rr}
\alias{raw.rr}
\title{raw.rr}
\format{
.csv
}
\source{
Lim, S.S., Vos, T., Flaxman, A.D., Danaei, G., Shibuya, K., Adair-Rohani, H., AlMazroa, M.A., Amann, M., Anderson, H.R., Andrews, K.G. and Aryee, M., 2012. A comparative risk assessment of burden of disease and injury attributable to 67 risk factors and risk factor clusters in 21 regions, 1990–2010: a systematic analysis for the Global Burden of Disease Study 2010. The lancet, 380(9859), pp.2224-2260.
}
\usage{
raw.rr
}
\description{
Relative risk of death attributable to a change in population-weighted mean pollutant concentration. From Van Dingenen et al (2018): "RR for PM2:5 exposure is calculated from the integrated exposure-response (IER) functions developed by Burnett et al. (2014) and first applied in the GBD study (Lim et al., 2012)
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.rr
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_aot40_nox}
\alias{src.soy_aot40_nox}
\title{src.soy_aot40_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_aot40_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and AOT40 (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_aot40_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_aot40_ch4}
\alias{src.soy_aot40_ch4}
\title{src.soy_aot40_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_aot40_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and AOT40 (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_aot40_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{a.Mi_MAIZE}
\alias{a.Mi_MAIZE}
\title{a.Mi_MAIZE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::a.Mi_MAIZE
}
}
\usage{
a.Mi_MAIZE
}
\description{
"a" coefficient for Mi-Maize
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_yll_pm25}
\alias{m3_get_yll_pm25}
\title{m3_get_yll_pm25}
\usage{
m3_get_yll_pm25(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
YLLs attributable to PM2.5 exposure for each TM5-FASST regions for all years (# YLLs). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce YLLs attributable to PM2.5 exposure. YLL-to-Mortalities ratios are based on TM5-FASST calculations. Premature mortalities are  based on the integrated exposure-response functions (IER) from Burnett et al (2014), consistent with the GBD 2016 study.
}
\keyword{PM2.5}
\keyword{YLL,}
\keyword{module_3,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{calc_daly_o3}
\alias{calc_daly_o3}
\title{calc_daly_o3}
\source{
Institute for Health Metrics and Evaluation (http://www.healthdata.org/)
}
\usage{
calc_daly_o3()
}
\value{
DALY-to-Mortality ratios for TM5-FASST regions for all years and O3-related causes (respiratory disease).The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Get the DALY-to-Mortality ratios used to estimate the Disability Adjusted Life Years (DALYs) attributable to ozone (O3) exposure.
}
\keyword{DALYs,}
\keyword{O3}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{calc_mort_rates}
\alias{calc_mort_rates}
\title{calc_mort_rates}
\source{
https://www.who.int/healthinfo/global_burden_disease/cod_2008_sources_methods.pdf
}
\usage{
calc_mort_rates()
}
\value{
Baseline mortality rates for TM5-FASST regions for all years and causes (ALRI, COPD, LC, IHD, STROKE). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Get cause-specific baseline mortalities from stroke, ischemic heart disease (IHD), chronic obstructive pulmonary disease (COPD), acute lower respiratory illness diseases (ALRI) and lung cancer (LC).
}
\keyword{Baseline}
\keyword{mortality}
\keyword{rates}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vsly_eu_2014}
\alias{vsly_eu_2014}
\title{vsly_eu_2014}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Schlander, M., Schaefer, R. and Schwarz, O., 2017. Empirical studies on the economic value of a Statistical Life Year (VSLY) in Europe: what do they tell us?. Value in Health, 20(9), p.A666.
\dontrun{
 library(rfasst);
 rfasst::vsly_eu_2014
}
}
\usage{
vsly_eu_2014
}
\description{
Base Value of Statistical Life Year (VSLY) in $2014
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m2_conc.R
\name{m2_get_conc_o3}
\alias{m2_get_conc_o3}
\title{m2_get_conc_o3}
\usage{
m2_get_conc_o3(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  ch4_o3 = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{ch4_o3}{Includes the CH4 effect on O3 based on Fiore et al (2008).By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Produce ozone (O3) levels for each TM5-FASST regions for all years (ppb). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce ozone (O3) concentration levels based on re-scaled emission pathways.
}
\keyword{O3}
\keyword{concentration,}
\keyword{module_2,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.yll.o3}
\alias{raw.yll.o3}
\title{raw.yll.o3}
\format{
An object of class \code{data.frame} with 56 rows and 2 columns.
}
\source{
TM5-FASST
}
\usage{
raw.yll.o3
}
\description{
Years of Life Lost (YLLs) to Mortality ratios attributable to O3 exposure
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.yll.o3
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dis}
\alias{dis}
\title{dis}
\format{
An object of class \code{character} of length 5.
}
\source{
TM5-FASST
\dontrun{
 library(rfasst);
 rfasst::dis
}
}
\usage{
dis
}
\description{
List of diseases for reading relative risk
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_mort_o3_ecoloss}
\alias{m3_get_mort_o3_ecoloss}
\title{m3_get_mort_o3_ecoloss}
\source{
Jerrett, M., Burnett, R.T., Pope III, C.A., Ito, K., Thurston, G., Krewski, D., Shi, Y., Calle, E. and Thun, M., 2009. Long-term ozone exposure and mortality. New England Journal of Medicine, 360(11), pp.1085-1095.//// Markandya, A., Sampedro, J., Smith, S.J., Van Dingenen, R., Pizarro-Irizar, C., Arto, I. and González-Eguino, M., 2018. Health co-benefits from air pollution and mitigation costs of the Paris Agreement: a modelling study. The Lancet Planetary Health, 2(3), pp.e126-e133.
}
\usage{
m3_get_mort_o3_ecoloss(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  ssp = "SSP2",
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Economic damages associated with mortality attributable to O3 (M6M) exposure for each TM5-FASST regions for all years (Million$2015). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Produce economic damages associated with premature mortality attributable to O3 (M6M) exposure based on the IER functions from Jerret et al (2009), consistent with the GBD 2016 study. The economic valuation takes as a base value the widely accepted Value of Statistical Life (VSL) of the OECD for 2005. This value, according to the literature ranges between US$1.8 and $4.5 million. The calculations for all regions are based on the  “unit value transfer approach” which adjusts the VSL according to their GDP and GDP growth rates. (Markandya et al 2018)
}
\keyword{,premature}
\keyword{O3}
\keyword{VSL}
\keyword{module_3,}
\keyword{mortality,}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{m3_get_mort_pm25}
\alias{m3_get_mort_pm25}
\title{m3_get_mort_pm25}
\source{
Burnett, R.T., Pope III, C.A., Ezzati, M., Olives, C., Lim, S.S., Mehta, S., Shin, H.H., Singh, G., Hubbell, B., Brauer, M. and Anderson, H.R., 2014. An integrated risk function for estimating the global burden of disease attributable to ambient fine particulate matter exposure. Environmental health perspectives, 122(4), pp.397-403.
}
\usage{
m3_get_mort_pm25(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  ssp = "SSP2",
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{ssp}{Set the ssp narrative associated to the GCAM scenario. c("SSP1","SSP2","SSP3","SSP4","SSP5"). By default is SSP2}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Premature mortality attributable to PM2.5 exposure for each TM5-FASST regions for all years (# mortalities).
}
\description{
Produce premature mortality attributable to PM2.5 exposure based on the integrated exposure-response functions (IER) from Burnett et al (2014), consistent with the GBD 2016 study.
}
\keyword{PM2.5}
\keyword{module_3,}
\keyword{mortality,}
\keyword{premature}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_aot40_nmvoc}
\alias{src.wheat_aot40_nmvoc}
\title{src.wheat_aot40_nmvoc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_aot40_nmvoc
}
\description{
Source-receptor coefficients (SRC) between NMVOC emissions and AOT40 (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_aot40_nmvoc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.o3_so2}
\alias{src.o3_so2}
\title{src.o3_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.o3_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and O3 concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.o3_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cf_o3}
\alias{cf_o3}
\title{cf_o3}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Jerrett, M., Burnett, R.T., Pope III, C.A., Ito, K., Thurston, G., Krewski, D., Shi, Y., Calle, E. and Thun, M., 2009. Long-term ozone exposure and mortality. New England Journal of Medicine, 360(11), pp.1085-1095.
\dontrun{
 library(rfasst);
 rfasst::cf_o3
}
}
\usage{
cf_o3
}
\description{
Counterfactual threshold for ozone (Jerret et al 2009)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coef.Mi_MAIZE}
\alias{coef.Mi_MAIZE}
\title{coef.Mi_MAIZE}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::coef.Mi_MAIZE
}
}
\usage{
coef.Mi_MAIZE
}
\description{
Coefficient for Mi-Maize
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONV_MIL}
\alias{CONV_MIL}
\title{CONV_MIL}
\format{
An object of class \code{numeric} of length 1.
}
\usage{
CONV_MIL
}
\description{
Transform $Million to $
\dontrun{
 library(rfasst);
 rfasst::CONV_MIL
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{b.Mi_WHEAT}
\alias{b.Mi_WHEAT}
\title{b.Mi_WHEAT}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::b.Mi_WHEAT
}
}
\usage{
b.Mi_WHEAT
}
\description{
"b" coefficient for Mi-Maize
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.no3_so2}
\alias{src.no3_so2}
\title{no3_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.no3_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and NO3 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.no3_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{a.Mi_SOY}
\alias{a.Mi_SOY}
\title{a.Mi_SOY}
\format{
An object of class \code{numeric} of length 1.
}
\source{
Wang, X. and Mauzerall, D.L., 2004. Characterizing distributions of surface ozone and its impact on grain production in China, Japan and South Korea: 1990 and 2020. Atmospheric Environment, 38(26), pp.4383-4402.
\dontrun{
 library(rfasst);
 rfasst::a.Mi_SOY
}
}
\usage{
a.Mi_SOY
}
\description{
"a" coefficient for Mi-Soy
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.so4_so2}
\alias{src.so4_so2}
\title{so4_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.so4_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and SO4 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.so4_so2
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m4_agriculture.R
\name{calc_price_gcam}
\alias{calc_price_gcam}
\title{calc_price_gcam}
\source{
GCAM
}
\usage{
calc_price_gcam(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files.By default=T}

\item{map}{Produce the maps. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Future agricultural price levels from GCAM for each GCAM region for all years ($1975/kg). The countries that form each GCAM region are listed in the documentation repository: https://github.com/JGCRI/gcam-doc/blob/gh-pages/overview.md The list of commodities within each category can be found in: Kyle, G.P., Luckow, P., Calvin, K.V., Emanuel, W.R., Nathan, M. and Zhou, Y., 2011. GCAM 3.0 agriculture and land use: data sources and methods (No. PNNL-21025). Pacific Northwest National Lab.(PNNL), Richland, WA (United States), and in Sampedro, J., Waldhoff, S.T., Van de Ven, D.J., Pardo, G., Van Dingenen, R., Arto, I., del Prado, A. and Sanz, M.J., 2020. Future impacts of ozone driven damages on agricultural systems. Atmospheric Environment, 231, p.117538, Table S2.
}
\description{
Extract future agricultural price levels from GCAM
}
\keyword{agriculture,}
\keyword{module_4,}
\keyword{price}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.wheat_mi_ch4}
\alias{src.wheat_mi_ch4}
\title{src.wheat_mi_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.wheat_mi_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and Mi (O3) concentration level for wheat, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.wheat_mi_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.base_conc}
\alias{raw.base_conc}
\title{Base concentration}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
raw.base_conc
}
\description{
Fine Particulate Matter (PM2.5), and Ozone (O3 and M6M) concentration levels in the base year (2000) per TM5-FASST region
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.base_conc
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.rice_mi_ch4}
\alias{src.rice_mi_ch4}
\title{src.rice_mi_ch4}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.rice_mi_ch4
}
\description{
Source-receptor coefficients (SRC) between CH4 emissions and Mi (O3) concentration level for rice, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.rice_mi_ch4
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.no3_nox}
\alias{src.no3_nox}
\title{no3_nox}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.no3_nox
}
\description{
Source-receptor coefficients (SRC) between NOx emissions and NO3 (PM2.5) concentration level, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.no3_nox
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m3_health.R
\name{calc_daly_tot}
\alias{calc_daly_tot}
\title{calc_daly_tot}
\source{
Institute for Health Metrics and Evaluation (http://www.healthdata.org/)
}
\usage{
calc_daly_tot()
}
\value{
DALY-to-Mortality ratios for TM5-FASST regions for all years and causes (ALRI, COPD, LC, IHD, STROKE).The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Get the DALY-to-Mortality ratios used to estimate the Disability Adjusted Life Years (DALYs).
}
\keyword{DALYs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw.base_mi}
\alias{raw.base_mi}
\title{Base Mi conc}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
raw.base_mi
}
\description{
Seasonal mean daytime O3 concentration (M7 for the 7-hour mean and M12 for the 12-hour mean), in the base year (2000) per TM5-FASST region
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::raw.base_mi
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m1_emissions_rescale.R
\name{m1_emissions_rescale}
\alias{m1_emissions_rescale}
\title{m1_emissions_rescale}
\usage{
m1_emissions_rescale(
  db_path,
  query_path,
  db_name,
  prj_name,
  scen_name,
  queries,
  saveOutput = T,
  map = F,
  mapIndivPol = F,
  anim = T
)
}
\arguments{
\item{db_path}{Path to the GCAM database}

\item{query_path}{Path to the query file}

\item{db_name}{Name of the GCAM database}

\item{prj_name}{Name of the rgcam project. This can be an existing project, or, if not, this will be the name}

\item{scen_name}{Name of the GCAM scenario to be processed}

\item{queries}{Name of the GCAM query file. The file by default includes the queries required to run rfasst}

\item{saveOutput}{Writes the emission files. By default=T}

\item{map}{Produce the maps. By default=F}

\item{mapIndivPol}{If set to T, it produces the maps for individual pollutants. By default=F}

\item{anim}{If set to T, produces multi-year animations. By default=T}
}
\value{
Emissions per TM5-FASST region (Kg) for all the selected years (all_years). The list of countries that form each region and the full name of the region can be found in Table S2.2 in the TM5-FASST documentation paper: Van Dingenen, R., Dentener, F., Crippa, M., Leitao, J., Marmer, E., Rao, S., Solazzo, E. and Valentini, L., 2018. TM5-FASST: a global atmospheric source-receptor model for rapid impact analysis of emission changes on air quality and short-lived climate pollutants. Atmospheric Chemistry and Physics, 18(21), pp.16173-16211.
}
\description{
Re-scale emissions from GCAM regions to TM5-FASST regions. The function also completes some additional transformation required such as pollutant transformation (e.g. OC to POM) or unit changes (Tg to Kg).
}
\keyword{emissions}
\keyword{module_1,}
\keyword{re-scale}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{src.soy_aot40_so2}
\alias{src.soy_aot40_so2}
\title{src.soy_aot40_so2}
\format{
.csv
}
\source{
Results from TM5
}
\usage{
src.soy_aot40_so2
}
\description{
Source-receptor coefficients (SRC) between SO2 emissions and AOT40 (O3) concentration level for soy, pre-computed from a 20% emission reduction of component i  in region xk relative to the base run
}
\examples{
\dontrun{
 library(rfasst);
 rfasst::src.soy_aot40_so2
}
}
\keyword{datasets}
