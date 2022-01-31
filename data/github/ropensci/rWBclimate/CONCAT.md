`rWBclimate`
==========

[![Build Status](https://api.travis-ci.org/ropensci/rWBclimate.png)](https://travis-ci.org/ropensci/rWBclimate)
[![Build status](https://ci.appveyor.com/api/projects/status/28njj5uw980frlpu/branch/master)](https://ci.appveyor.com/project/sckott/rwbclimate/branch/master)
[![codecov.io](https://codecov.io/github/ropensci/rWBclimate/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rWBclimate?branch=master)

Introduction
========================================================
rWBclimate is an R interface for the World Bank climate data used in the World Bank [climate knowledge portal](http://sdwebx.worldbank.org/climateportal/index.cfm).  

Installation
------
Right now the package is only installable from github with [devtools](http://cran.r-project.org/web/packages/devtools/index.html):

```R
require(rWBclimate)
```
Package description
----

`rWBclimate` provides access to three different classes of climate data at two different spatial scales.  The three different classes of data are GCM model , ensemble and historical data.  Each data class will let you download two different four different types for two different variables.  The two variables are either precipitation expressed in millimeters or temperature as degrees celcius, and for each variable you can download your data in one of four types. The data is fully described below along with examples.

Data classes
---
*__Model data__*

Almost all model data in the Climate Data API are derived from 15 global circulation models (GCMs) used by the Intergovernmental Panel on Climate Change (IPCC) 4th Assessment Reports. The models simulate the response of the global climate system to increasing greenhouse gas concentrations. The data in the Climate Data API have been aggregated to both the country and basin levels, as explained below. Note these data are modeled estimates of temperature and precipitation changes in different time periods under different GCMs and scenarios. They include changes for future time periods and also as “backcasting” (model representations of the past) set for past time periods. The latter should not be confused with any instrumental or observed data. There is a specific dataset with historical measured climate data as well.

**Data types**

|Type|Description|
|----|----|
|Monthly average|The monthly average for all 12 months for a given time period|
|Annual average|a single average for a given time period|
|Monthly anomaly|Average monthly change (anomaly).  The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|
|Annual anomaly|Average annual change (anomaly). The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|

**Data time scales**

Climate model data is only available as averages for 20 year chunks.  The package will automatically convert any start and end data into valid API calls and will return all data between the given start and end date.  The following time periods are available


|Past|     |Future|    |
|----|-----|------|----|
|*start*|  *end*  |  *start*  |*end*|
|1920  | 1939|  2020 | 2039 |
|1940 |  1959|  2040 |2059|
|1960 |   1979|  2060 |2079 |
|1980  | 1999| 2080 | 2099 |

**Data spatial scales**
Data is available at two spatial scales.  The first is country level. You can download data for any country in the world using a valid [ISO 3 letter country code](http://userpage.chemie.fu-berlin.de/diverse/doc/ISO_3166.html). Alternatively you can download data at the basin network for a slightly higher resolution represented as a number 1-468.  This is based on level 2 boundaries from [waterbase.org](http://www.waterbase.org), or you can view a [full map of all the available basins.](http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)


**_Downloading GCM Model Data_**

Model data is downloaded for two different scenarios, the [A2 and B1](http://en.wikipedia.org/wiki/Special_Report_on_Emissions_Scenarios). Generally A2 scenarios are where there is little difference between the future and now and B1 is a more ecologically friendly world with greater decrease in emissions. Both the A2 and the B1 will be downloaded for 15 different GCM models listed in the table below:


|Name in output|Model name|
|--------------|----------|
|bccr_bcm2_0 |[BCM 2.0](http://www-pcmdi.llnl.gov/ipcc/model_documentation/BCCR_BCM2.0.htm)|
|csiro_mk3_5|[CSIRO Mark 3.5](http://www.cawcr.gov.au/publications/technicalreports/CTR_021.pdf)|
|ingv_echam4|[ECHAM 4.6](http://www.bo.ingv.it/)|
|cccma_cgcm3_1|[CGCM 3.1 (T47)](http://www.ec.gc.ca/ccmac-cccma/default.asp?lang=En)|
|cnrm_cm3|[CNRM CM3](http://www.cnrm.meteo.fr/scenario2004/indexenglish.html)|
|gfdl_cm2_0|[GFDL CM2.0](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|gfdl_cm2_1|[GFDL CM2.1](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|ipsl_cm4|[IPSL-CM4](http://mc2.ipsl.jussieu.fr/simules.html)|
|microc3_2_medres|[MIROC 3.2 (medres)](https://esg.llnl.gov:8443/metadata/browseCatalog.do?uri=http://esgcet.llnl.gov/metadata/pcmdi/ipcc/thredds/miroc3_2_medres.sresb1/pcmdi.ipcc4.miroc3_2_medres.sresb1.thredds)|
|miub_echo_g|[ECHO-G](http://www-pcmdi.llnl.gov/projects/modeldoc/cmip/echo-g_tbls.html)|
|mpi_echam5|[ECHAM5/MPI-OM](http://www.mpimet.mpg.de/en/science/models/echam.html)|
|mri_cgcm2_3_2a|[MRI-CGCM2.3.2](http://www.mri-jma.go.jp/Welcome.html)|
|inmcm3_0|[INMCM3.0](http://www.ipcc-data.org/ar4/model-INM-CM3.html)|
|ukmo_hadcm3|[UKMO HadCM3](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadcm3)|
|ukmo_hadgem1|[UKMO HadGEM1](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadgem1)|

The model data can be downloaded with two main functions:
```R
get_model_temp()  ## Get model temperature data
get_model_precip() ## Get model precipitation data
```
*Example 1: Plotting monthly data from different GCM's*
Say you want to compare temperature from two different models in the USA to see how they vary.  You can download data for the USA and then subset it to the specific models you're interested in and then plot them.





```r
usa.dat <- get_model_temp("USA", "mavg", 2080, 2100)
usa.dat.bcc <- usa.dat[usa.dat$gcm == "bccr_bcm2_0", ]
usa.dat.had <- usa.dat[usa.dat$gcm == "ukmo_hadcm3", ]
## Add a unique ID to each for easier plotting
usa.dat.bcc$ID <- paste(usa.dat.bcc$scenario, usa.dat.bcc$gcm, sep = "-")
usa.dat.had$ID <- paste(usa.dat.had$scenario, usa.dat.had$gcm, sep = "-")
plot.df <- rbind(usa.dat.bcc, usa.dat.had)
ggplot(plot.df, aes(x = as.factor(month), y = data, group = ID, colour = gcm,
    linetype = scenario)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") +
    xlab("Month") + theme_bw()
```

![plot of chunk getmodeldata](figure/getmodeldata.png)


Subsetting all the data can be a bit tedious.  You could also compare all the models but just for one scenario, the A2.


```r
ggplot(usa.dat[usa.dat$scenario == "a2", ], aes(x = month, y = data, group = gcm,
    colour = gcm)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") +
    xlab("Month") + theme_bw()
```

![plot of chunk plotallmodeldata](figure/plotallmodeldata.png)


*Example 2: Plotting annual data for different countries*

Data can be extracted from countries or basins submitted as vectors. Here we will plot the expected temperature anomaly for each 20 year period over a baseline control period of 1961-2000.  These countries chosen span the north to south pole.  It's clear from the plot that the northern most countries (US and Canada) have the biggest anomaly, and Belize, the most equatorial country, has the smallest anomaly.

```r
country.list <- c("CAN", "USA", "MEX", "BLZ", "COL", "PER", "BOL", "ARG")
country.dat <- get_model_temp(country.list, "annualanom", 2010, 2100)
# Subset data
country.dat.bcc <- country.dat[country.dat$gcm == "bccr_bcm2_0", ]
## Exclude A2 scenario
country.dat.bcc <- subset(country.dat.bcc, country.dat.bcc$scenario != "a2")
ggplot(country.dat.bcc, aes(x = fromYear, y = data, group = locator, colour = locator)) +
    geom_point() + geom_path() + ylab("Temperature anomaly over baseline") +
    theme_bw()
```

![plot of chunk annualdata](figure/annualdata.png)



**_Downloading Ensemble Model Data_**

Getting raw model data is really useful for comparing scenarios or different GCM's.  However many users might be more interested in aggregated measurements. If so you can download ensemble climate data.  This will give access to all 15 GCM's combined and queries will return the 10th, 50th and 90th quantiles for the ensemble of all GCM's.  All queries are constructed the same as raw model data above with the same data types, spatial scales and time scales.

*Example 1: Comparing model quantiles*

Let's look at monthly precipitation predictions for Indonesia for the period of 2080-2100.  We'll create a plot of the precipitation expected in the future for the two different scenarios and plot them with different colours and dashed lines to indicate quantiles.


```r
idn.dat <- get_ensemble_precip("IDN", "mavg", 2080, 2100)
# Set line types
ltype <- rep(1, dim(idn.dat)[1])
ltype[idn.dat$percentile != 50] <- 2
idn.dat$ltype <- ltype
# Create uniqueIDs
idn.dat$uid <- paste(idn.dat$scenario, idn.dat$percentile, sep = "-")
ggplot(idn.dat, aes(x = as.factor(month), y = data, group = uid, colour = scenario,
    linetype = as.factor(ltype))) + geom_point() + geom_path() + xlab("Month") +
    ylab("Rain in mm") + theme_bw()
```

![plot of chunk comparing quantiles](figure/comparing_quantiles.png)


*Example 2: Ensemble statistics*

You can also download 13 different ensemble statistics about basins or countries aside from the raw data as above.  These derived statistics are often given relative to a control period, 1961-2000.  The time periods however are only for 2046-2065, or 2081-2100, not the standard 20 year intervals available for raw model data.  Below is a full table of available statistics.

|Parameter name|Description|Units|
|---------------|-------|-------|
|*tmin_means*|Average daily minimum temperature|degrees Celsius|
|*tmax_means*|Average daily maximum temperature|degrees Celsius|
|*tmax_days90th*|Number of days with maximum temperature above the control period’s 90th percentile (hot days)|days|
|*tmin_days90th*|Number of days with minimum  temperature above the control period’s 90th percentile (warm nights)|days|
|*tmax_days10th*|Number of days with maximum temperature below the control period’s 10th percentile (cool days)|days|
|*tmin_days10th*|Number of days with minimum  temperature below the control period’s 10th percentile (cold nights)|days|
|*tmin_days0*|Number of days with minimum  temperature below 0 degrees Celsius|days|
|*ppt_days*|Number of days with precipitation greater than 0.2 mm|days|
|*ppt_days2*|Number of days with precipitation greater than 2 mm|days|
|*ppt_days10*|Number of days with precipitation greater than 10 mm|days|
|*ppt_days90th*|Number of days with precipitation greater than the control period's 90th percentile|days|
|*ppt_dryspell*|Average number of days between precipitation events|days|
|*ppt_means*|Average daily precipitation|mm|

Similar to our previous example where we looked at temperature anomaly along a latitudinal gradient, we can examine more complex ensemble statistics.  Here we examine the average minimum daily temperature in Scavavadian countries.  Also given that only two time periods are available, there is no need to enter in a specific year range

```r
country.list <- c("ISL", "FIN", "NOR", "SWE")
country.dat <- get_ensemble_stats(country.list, "mavg", "tmin_means")
####### Subset data Exclude A2 scenario
country.dat.b1 <- subset(country.dat, country.dat$scenario == "b1")
# choose just one percentile
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$percentile == 50)
# get just one year period
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$fromYear == 2081)


ggplot(country.dat.b1, aes(x = month, y = data, group = locator, colour = locator)) +
    geom_point() + geom_path() + ylab("Average daily minimum temperature") +
    theme_bw() + xlab("Month")
```

![plot of chunk enesmble data](figure/enesmble_data.png)



**_Downloading Historical Data_**

It is possible to extract historical data from GCM queries, but this is acutally model backcasting output, not true historical records.  The climate data api can download data at the same spatial scales as all other requests (country or basin), but the time scale differs. You can download data at monthly, yearly or decadanal time scales.  Monthly data is actually the mean monthly temperature or precipitation value for each month from 1901-2009 for countries, or 1960-2009 for basins.  You can indicate the type of data you want in the call by setting the `time_scale` parameter to *month*,*year* or *decade*.

*Example 1: Downloading monthly data*

You can download historical precipitation for Belize, Colombia, Peru and Bolivia.  



```r
country.list <- c("BLZ", "COL", "PER", "BOL")
country.dat <- get_historical_precip(country.list, "month")

ggplot(country.dat, aes(x = month, y = data, group = locator, colour = locator)) +
    geom_point() + geom_path() + ylab("Average historical precipitation (mm)") +
    theme_bw() + xlab("Month")
```

![plot of chunk historicalmonth](figure/historicalmonth.png)


*Example 2: Downloading annual data*

Another use of historical data is to look at increases in temperature over time.

```r
country.list <- c("USA", "MEX", "CAN", "BLZ")
country.dat <- get_historical_temp(country.list, "year")

ggplot(country.dat, aes(x = year, y = data, group = locator)) + geom_point() +
    geom_path() + ylab("Average annual temperature of Canada") + theme_bw() +
    xlab("Year") + stat_smooth(se = F, colour = "black") + facet_wrap(~locator,
    scale = "free")
```

![plot of chunk historicalyear](figure/historicalyear.png)


**_Mapping climate Data_**

Data can be mapped in [ggplot2](http://docs.ggplot2.org/current/) by downloading KML files from the climate database, creating dataframes and plotting using ggplot2. Maps are available at the basin or country spatial scale.  Because KML files can be large files are downloaded and locally cached in a directory you must set using the kmlpath option.  You can set to any directory as follows: `options(kmlpath = <yourpath>)`.  KML files will be stored there and only downloaded again if they aren't found locally.  If you wish to get rid of them, you will need to manually delete them.  `rWBclimate` allows you to download maps and plot them as ggplot maps, or it can automatically create maps of climate data.  For your convenience we've created basin and country vectors for all six continents *(Antartica has no data in the database)*.  These data files are automatically loaded with the package and can be accessed for basins as *NoAm_basin*, *SoAm_basin*,*Asia_basin*,*Eur_basin*,*Africa_basin*,and *Oceana_basin*. Countries can be similarly accessed as *NoAm_country*, *SoAm_country*,*Asia_country*,*Eur_country*,*Africa_country*,and *Oceana_country*.  If you are plotting by continent, it can take some time to first download all the KML files.

*Example 1: Creating map dataframe and plotting*

Creating map data frames is straightforward, simply provide a list of valid country codes to the function `create_map_df`

```r
# Set the kmlpath option
options(kmlpath = "~/kmltemp")
## Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)
```

```r
ggplot(af_basin, aes(x = long, y = lat, group = group)) + geom_polygon() + theme_bw()
```

![plot of chunk map_plot](figure/map_plot.png)


*Example 2: Mapping climate data*

In order to map climate data you need to have a single point of data for each spatial polygon(country or basin) in your map dataframe.  We've already created the Africa basin map dataframe, so now you need to get some data for each basin.  


```r
af_basin_dat <- get_ensemble_temp(Africa_basin, "annualanom", 2080, 2100)
## Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat, af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat, af_basin_dat$percentile == 50)
```


Now that we have both the map dataframe and the data we can bind them together with the function `climate_map()`.  The function has two options, it can return a `ggplot2` map or it can return a dataframe that you can plot easily with `ggplot2`.  This is useful because it means that you could bind together multiple dataframes for bigger plots, or perhaps add new data yourself


```r

af_map <- climate_map(af_basin, af_basin_dat, return_map = T)
af_map + scale_fill_continuous("Temperature \n anomaly", low = "yellow", high = "red") +
    theme_bw()
```

![plot of chunk climatemap](figure/climatemap.png)



*Example 3: Creating a temperature map of the world*


```r
options(kmlpath = "~/kmltemp")
### Combine all country vectors

world <- c(NoAm_country, SoAm_country, Eur_country, Asia_country, Africa_country,
    Oceana_country)
world_map_df <- create_map_df(world)
```

```r
world_dat <- get_ensemble_temp(world, "annualavg", 2080, 2100)
## Subset data to just one scenario, and one percentile
world_dat <- subset(world_dat, world_dat$scenario == "a2")
world_dat <- subset(world_dat, world_dat$percentile == 50)

world_map <- climate_map(world_map_df, world_dat, return_map = T)
world_map + scale_fill_continuous("Temperature \n at 2100", low = "yellow",
    high = "red") + theme_bw()
```

![plot of chunk worldmap](figure/worldmap.png)


To cite package ‘rWBclimate’ in publications use:

```coffee

  Edmund Hart (). rWBclimate: A package for accessing World Bank climate data. R package version 0.1.3
  http://github.com/ropensci/rWBclimate

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {rWBclimate: A package for accessing World Bank climate data},
    author = {Edmund Hart},
    note = {R package version 0.1.3},
    url = {http://github.com/ropensci/rWBclimate},
  }
```


[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rWBclimate}
-->
Introduction 
========================================================
rWBclimate is an R interface for the World Bank climate data used in the World Bank [climate knowledge portal](http://sdwebx.worldbank.org/climateportal/index.cfm).  

Installation
------
Right now the package is only installable from github with [devtools](http://cran.r-project.org/web/packages/devtools/index.html):

```R
require(rWBclimate)
```
Package description
----

`rWBclimate` provides access to three different classes of climate data at two different spatial scales.  The three different classes of data are GCM model , ensemble and historical data.  Each data class will let you download two different four different types for two different variables.  The two variables are either precipitation expressed in millimeters or temperature as degrees celcius, and for each variable you can download your data in one of four types. The data is fully described below along with examples.

Data classes
---
*__Model data__*

Almost all model data in the Climate Data API are derived from 15 global circulation models (GCMs) used by the Intergovernmental Panel on Climate Change (IPCC) 4th Assessment Reports. The models simulate the response of the global climate system to increasing greenhouse gas concentrations. The data in the Climate Data API have been aggregated to both the country and basin levels, as explained below. Note these data are modeled estimates of temperature and precipitation changes in different time periods under different GCMs and scenarios. They include changes for future time periods and also as "backcasting" (model representations of the past) set for past time periods. The latter should not be confused with any instrumental or observed data. There is a specific dataset with historical measured climate data as well.

**Data types**

|Type|Description|
|----|----|
|Monthly average|The monthly average for all 12 months for a given time period|
|Annual average|a single average for a given time period|
|Monthly anomaly|Average monthly change (anomaly).  The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|
|Annual anomaly|Average annual change (anomaly). The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|

**Data time scales**

Climate model data is only available as averages for 20 year chunks.  The package will automatically convert any start and end data into valid API calls and will return all data between the given start and end date.  The following time periods are available


|Past|     |Future|    |
|----|-----|------|----|
|*start*|  *end*  |  *start*  |*end*|
|1920  | 1939|  2020 | 2039 |
|1940 |  1959|  2040 |2059|
|1960 |   1979|  2060 |2079 |
|1980  | 1999| 2080 | 2099 |

**Data spatial scales**
Data is available at two spatial scales.  The first is country level. You can download data for any country in the world using a valid [ISO 3 letter country code](http://userpage.chemie.fu-berlin.de/diverse/doc/ISO_3166.html). Alternatively you can download data at the basin network for a slightly higher resolution represented as a number 1-468.  This is based on level 2 boundaries from [waterbase.org](http://www.waterbase.org), or you can view a [full map of all the available basins.](http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)


**_Downloading GCM Model Data_**

Model data is downloaded for two different scenarios, the [A2 and B1](http://en.wikipedia.org/wiki/Special_Report_on_Emissions_Scenarios). Generally A2 scenarios are where there is little difference between the future and now and B1 is a more ecologically friendly world with greater decrease in emissions. Both the A2 and the B1 will be downloaded for 15 different GCM models listed in the table below:


|Name in output|Model name|
|--------------|----------|
|bccr_bcm2_0 |[BCM 2.0](http://www-pcmdi.llnl.gov/ipcc/model_documentation/BCCR_BCM2.0.htm)|
|csiro_mk3_5|[CSIRO Mark 3.5](http://www.cawcr.gov.au/publications/technicalreports/CTR_021.pdf)|
|ingv_echam4|[ECHAM 4.6](http://www.bo.ingv.it/)|
|cccma_cgcm3_1|[CGCM 3.1 (T47)](http://www.ec.gc.ca/ccmac-cccma/default.asp?lang=En)|
|cnrm_cm3|[CNRM CM3](http://www.cnrm.meteo.fr/scenario2004/indexenglish.html)|
|gfdl_cm2_0|[GFDL CM2.0](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|gfdl_cm2_1|[GFDL CM2.1](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|ipsl_cm4|[IPSL-CM4](http://mc2.ipsl.jussieu.fr/simules.html)|
|microc3_2_medres|[MIROC 3.2 (medres)](https://esg.llnl.gov:8443/metadata/browseCatalog.do?uri=http://esgcet.llnl.gov/metadata/pcmdi/ipcc/thredds/miroc3_2_medres.sresb1/pcmdi.ipcc4.miroc3_2_medres.sresb1.thredds)|
|miub_echo_g|[ECHO-G](http://www-pcmdi.llnl.gov/projects/modeldoc/cmip/echo-g_tbls.html)|
|mpi_echam5|[ECHAM5/MPI-OM](http://www.mpimet.mpg.de/en/science/models/echam.html)|
|mri_cgcm2_3_2a|[MRI-CGCM2.3.2](http://www.mri-jma.go.jp/Welcome.html)|
|inmcm3_0|[INMCM3.0](http://www.ipcc-data.org/ar4/model-INM-CM3.html)|
|ukmo_hadcm3|[UKMO HadCM3](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadcm3)|
|ukmo_hadgem1|[UKMO HadGEM1](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadgem1)|

The model data can be downloaded with two main functions:
```R
get_model_temp()  ## Get model temperature data
get_model_precip() ## Get model precipitation data
```
*Example 1: Plotting monthly data from different GCM's* 
Say you want to compare temperature from two different models in the USA to see how they vary.  You can download data for the USA and then subset it to the specific models you're interested in and then plot them.





```r
usa.dat <- get_model_temp("USA", "mavg", 2080, 2100)
usa.dat.bcc <- usa.dat[usa.dat$gcm == "bccr_bcm2_0", ]
usa.dat.had <- usa.dat[usa.dat$gcm == "ukmo_hadcm3", ]
## Add a unique ID to each for easier plotting
usa.dat.bcc$ID <- paste(usa.dat.bcc$scenario, usa.dat.bcc$gcm, sep = "-")
usa.dat.had$ID <- paste(usa.dat.had$scenario, usa.dat.had$gcm, sep = "-")
plot.df <- rbind(usa.dat.bcc, usa.dat.had)
ggplot(plot.df, aes(x = as.factor(month), y = data, group = ID, colour = gcm, 
    linetype = scenario)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") + 
    xlab("Month") + theme_bw()
```

![plot of chunk getmodeldata](figure/getmodeldata.png) 


Subsetting all the data can be a bit tedious.  You could also compare all the models but just for one scenario, the A2.


```r
ggplot(usa.dat[usa.dat$scenario == "a2", ], aes(x = month, y = data, group = gcm, 
    colour = gcm)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") + 
    xlab("Month") + theme_bw()
```

![plot of chunk plotallmodeldata](figure/plotallmodeldata.png) 


*Example 2: Plotting annual data for different countries*

Data can be extracted from countries or basins submitted as vectors. Here we will plot the expected temperature anomaly for each 20 year period over a baseline control period of 1961-2000.  These countries chosen span the north to south pole.  It's clear from the plot that the northern most countries (US and Canada) have the biggest anomaly, and Belize, the most equatorial country, has the smallest anomaly.

```r
country.list <- c("CAN", "USA", "MEX", "BLZ", "COL", "PER", "BOL", "ARG")
country.dat <- get_model_temp(country.list, "annualanom", 2010, 2100)
# Subset data
country.dat.bcc <- country.dat[country.dat$gcm == "bccr_bcm2_0", ]
## Exclude A2 scenario
country.dat.bcc <- subset(country.dat.bcc, country.dat.bcc$scenario != "a2")
ggplot(country.dat.bcc, aes(x = fromYear, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Temperature anomaly over baseline") + 
    theme_bw()
```

![plot of chunk annualdata](figure/annualdata.png) 



**_Downloading Ensemble Model Data_**

Getting raw model data is really useful for comparing scenarios or different GCM's.  However many users might be more interested in aggregated measurements. If so you can download ensemble climate data.  This will give access to all 15 GCM's combined and queries will return the 10th, 50th and 90th quantiles for the ensemble of all GCM's.  All queries are constructed the same as raw model data above with the same data types, spatial scales and time scales.

*Example 1: Comparing model quantiles*

Let's look at monthly precipitation predictions for Indonesia for the period of 2080-2100.  We'll create a plot of the precipitation expected in the future for the two different scenarios and plot them with different colours and dashed lines to indicate quantiles.


```r
idn.dat <- get_ensemble_precip("IDN", "mavg", 2080, 2100)
# Set line types
ltype <- rep(1, dim(idn.dat)[1])
ltype[idn.dat$percentile != 50] <- 2
idn.dat$ltype <- ltype
# Create uniqueIDs
idn.dat$uid <- paste(idn.dat$scenario, idn.dat$percentile, sep = "-")
ggplot(idn.dat, aes(x = as.factor(month), y = data, group = uid, colour = scenario, 
    linetype = as.factor(ltype))) + geom_point() + geom_path() + xlab("Month") + 
    ylab("Rain in mm") + theme_bw()
```

![plot of chunk comparing quantiles](figure/comparing_quantiles.png) 


*Example 2: Ensemble statistics*

You can also download 13 different ensemble statistics about basins or countries aside from the raw data as above.  These derived statistics are often given relative to a control period, 1961-2000.  The time periods however are only for 2046-2065, or 2081-2100, not the standard 20 year intervals available for raw model data.  Below is a full table of available statistics.

|Parameter name|Description|Units|
|---------------|-------|-------|
|*tmin_means*|Average daily minimum temperature|degrees Celsius|
|*tmax_means*|Average daily maximum temperature|degrees Celsius|
|*tmax_days90th*|Number of days with maximum temperature above the control period's 90th percentile (hot days)|days|
|*tmin_days90th*|Number of days with minimum  temperature above the control period's 90th percentile (warm nights)|days|
|*tmax_days10th*|Number of days with maximum temperature below the control period's 10th percentile (cool days)|days|
|*tmin_days10th*|Number of days with minimum  temperature below the control period's 10th percentile (cold nights)|days|
|*tmin_days0*|Number of days with minimum  temperature below 0 degrees Celsius|days|
|*ppt_days*|Number of days with precipitation greater than 0.2 mm|days|
|*ppt_days2*|Number of days with precipitation greater than 2 mm|days|
|*ppt_days10*|Number of days with precipitation greater than 10 mm|days|
|*ppt_days90th*|Number of days with precipitation greater than the control period's 90th percentile|days|
|*ppt_dryspell*|Average number of days between precipitation events|days|
|*ppt_means*|Average daily precipitation|mm|

Similar to our previous example where we looked at temperature anomaly along a latitudinal gradient, we can examine more complex ensemble statistics.  Here we examine the average minimum daily temperature in Scavavadian countries.  Also given that only two time periods are available, there is no need to enter in a specific year range

```r
country.list <- c("ISL", "FIN", "NOR", "SWE")
country.dat <- get_ensemble_stats(country.list, "mavg", "tmin_means")
####### Subset data Exclude A2 scenario
country.dat.b1 <- subset(country.dat, country.dat$scenario == "b1")
# choose just one percentile
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$percentile == 50)
# get just one year period
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$fromYear == 2081)


ggplot(country.dat.b1, aes(x = month, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Average daily minimum temperature") + 
    theme_bw() + xlab("Month")
```

![plot of chunk enesmble data](figure/enesmble_data.png) 



**_Downloading Historical Data_**

It is possible to extract historical data from GCM queries, but this is acutally model backcasting output, not true historical records.  The climate data api can download data at the same spatial scales as all other requests (country or basin), but the time scale differs. You can download data at monthly, yearly or decadanal time scales.  Monthly data is actually the mean monthly temperature or precipitation value for each month from 1901-2009 for countries, or 1960-2009 for basins.  You can indicate the type of data you want in the call by setting the `time_scale` parameter to *month*,*year* or *decade*.

*Example 1: Downloading monthly data*

You can download historical precipitation for Belize, Colombia, Peru and Bolivia.  



```r
country.list <- c("BLZ", "COL", "PER", "BOL")
country.dat <- get_historical_precip(country.list, "month")

ggplot(country.dat, aes(x = month, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Average historical precipitation (mm)") + 
    theme_bw() + xlab("Month")
```

![plot of chunk historicalmonth](figure/historicalmonth.png) 


*Example 2: Downloading annual data*

Another use of historical data is to look at increases in temperature over time.

```r
country.list <- c("USA", "MEX", "CAN", "BLZ")
country.dat <- get_historical_temp(country.list, "year")

ggplot(country.dat, aes(x = year, y = data, group = locator)) + geom_point() + 
    geom_path() + ylab("Average annual temperature of Canada") + theme_bw() + 
    xlab("Year") + stat_smooth(se = F, colour = "black") + facet_wrap(~locator, 
    scale = "free")
```

![plot of chunk historicalyear](figure/historicalyear.png) 


**_Mapping climate Data_**

Data can be mapped in [ggplot2](http://docs.ggplot2.org/current/) by downloading KML files from the climate database, creating dataframes and plotting using ggplot2. Maps are available at the basin or country spatial scale.  Because KML files can be large files are downloaded and locally cached in a directory you must set using the kmlpath option.  You can set to any directory as follows: `options(kmlpath = <yourpath>)`.  KML files will be stored there and only downloaded again if they aren't found locally.  If you wish to get rid of them, you will need to manually delete them.  `rWBclimate` allows you to download maps and plot them as ggplot maps, or it can automatically create maps of climate data.  For your convenience we've created basin and country vectors for all six continents *(Antartica has no data in the database)*.  These data files are automatically loaded with the package and can be accessed for basins as *NoAm_basin*, *SoAm_basin*,*Asia_basin*,*Eur_basin*,*Africa_basin*,and *Oceana_basin*. Countries can be similarly accessed as *NoAm_country*, *SoAm_country*,*Asia_country*,*Eur_country*,*Africa_country*,and *Oceana_country*.  If you are plotting by continent, it can take some time to first download all the KML files.

*Example 1: Creating map dataframe and plotting*

Creating map data frames is straightforward, simply provide a list of valid country codes to the function `create_map_df`

```r
# Set the kmlpath option
options(kmlpath = "~/kmltemp")
## Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)
```

```
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=                                                                |   1%
  |                                                                       
  |=                                                                |   2%
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |===                                                              |   4%
  |                                                                       
  |====                                                             |   6%
  |                                                                       
  |====                                                             |   7%
  |                                                                       
  |=====                                                            |   8%
  |                                                                       
  |======                                                           |   9%
  |                                                                       
  |======                                                           |  10%
  |                                                                       
  |=======                                                          |  11%
  |                                                                       
  |========                                                         |  12%
  |                                                                       
  |=========                                                        |  13%
  |                                                                       
  |=========                                                        |  14%
  |                                                                       
  |==========                                                       |  16%
  |                                                                       
  |===========                                                      |  17%
  |                                                                       
  |============                                                     |  18%
  |                                                                       
  |============                                                     |  19%
  |                                                                       
  |=============                                                    |  20%
  |                                                                       
  |==============                                                   |  21%
  |                                                                       
  |==============                                                   |  22%
  |                                                                       
  |===============                                                  |  23%
  |                                                                       
  |================                                                 |  24%
  |                                                                       
  |=================                                                |  26%
  |                                                                       
  |=================                                                |  27%
  |                                                                       
  |==================                                               |  28%
  |                                                                       
  |===================                                              |  29%
  |                                                                       
  |====================                                             |  30%
  |                                                                       
  |====================                                             |  31%
  |                                                                       
  |=====================                                            |  32%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |======================                                           |  34%
  |                                                                       
  |=======================                                          |  36%
  |                                                                       
  |========================                                         |  37%
  |                                                                       
  |=========================                                        |  38%
  |                                                                       
  |=========================                                        |  39%
  |                                                                       
  |==========================                                       |  40%
  |                                                                       
  |===========================                                      |  41%
  |                                                                       
  |===========================                                      |  42%
  |                                                                       
  |============================                                     |  43%
  |                                                                       
  |=============================                                    |  44%
  |                                                                       
  |==============================                                   |  46%
  |                                                                       
  |==============================                                   |  47%
  |                                                                       
  |===============================                                  |  48%
  |                                                                       
  |================================                                 |  49%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================                                |  51%
  |                                                                       
  |==================================                               |  52%
  |                                                                       
  |===================================                              |  53%
  |                                                                       
  |===================================                              |  54%
  |                                                                       
  |====================================                             |  56%
  |                                                                       
  |=====================================                            |  57%
  |                                                                       
  |======================================                           |  58%
  |                                                                       
  |======================================                           |  59%
  |                                                                       
  |=======================================                          |  60%
  |                                                                       
  |========================================                         |  61%
  |                                                                       
  |========================================                         |  62%
  |                                                                       
  |=========================================                        |  63%
  |                                                                       
  |==========================================                       |  64%
  |                                                                       
  |===========================================                      |  66%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |============================================                     |  68%
  |                                                                       
  |=============================================                    |  69%
  |                                                                       
  |==============================================                   |  70%
  |                                                                       
  |==============================================                   |  71%
  |                                                                       
  |===============================================                  |  72%
  |                                                                       
  |================================================                 |  73%
  |                                                                       
  |================================================                 |  74%
  |                                                                       
  |=================================================                |  76%
  |                                                                       
  |==================================================               |  77%
  |                                                                       
  |===================================================              |  78%
  |                                                                       
  |===================================================              |  79%
  |                                                                       
  |====================================================             |  80%
  |                                                                       
  |=====================================================            |  81%
  |                                                                       
  |=====================================================            |  82%
  |                                                                       
  |======================================================           |  83%
  |                                                                       
  |=======================================================          |  84%
  |                                                                       
  |========================================================         |  86%
  |                                                                       
  |========================================================         |  87%
  |                                                                       
  |=========================================================        |  88%
  |                                                                       
  |==========================================================       |  89%
  |                                                                       
  |==========================================================       |  90%
  |                                                                       
  |===========================================================      |  91%
  |                                                                       
  |============================================================     |  92%
  |                                                                       
  |=============================================================    |  93%
  |                                                                       
  |=============================================================    |  94%
  |                                                                       
  |==============================================================   |  96%
  |                                                                       
  |===============================================================  |  97%
  |                                                                       
  |================================================================ |  98%
  |                                                                       
  |================================================================ |  99%
  |                                                                       
  |=================================================================| 100%
```

```r
ggplot(af_basin, aes(x = long, y = lat, group = group)) + geom_polygon() + theme_bw()
```

![plot of chunk map_plot](figure/map_plot.png) 


*Example 2: Mapping climate data*

In order to map climate data you need to have a single point of data for each spatial polygon(country or basin) in your map dataframe.  We've already created the Africa basin map dataframe, so now you need to get some data for each basin.  


```r
af_basin_dat <- get_ensemble_temp(Africa_basin, "annualanom", 2080, 2100)
## Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat, af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat, af_basin_dat$percentile == 50)
```


Now that we have both the map dataframe and the data we can bind them together with the function `climate_map()`.  The function has two options, it can return a `ggplot2` map or it can return a dataframe that you can plot easily with `ggplot2`.  This is useful because it means that you could bind together multiple dataframes for bigger plots, or perhaps add new data yourself


```r

af_map <- climate_map(af_basin, af_basin_dat, return_map = T)
af_map + scale_fill_continuous("Temperature \n anomaly", low = "yellow", high = "red") + 
    theme_bw()
```

![plot of chunk climatemap](figure/climatemap.png) 



*Example 3: Creating a temperature map of the world*


```r
options(kmlpath = "~/kmltemp")
### Combine all country vectors

world <- c(NoAm_country, SoAm_country, Eur_country, Asia_country, Africa_country, 
    Oceana_country)
world_map_df <- create_map_df(world)
```

```
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=                                                                |   1%
  |                                                                       
  |=                                                                |   2%
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |==                                                               |   4%
  |                                                                       
  |===                                                              |   4%
  |                                                                       
  |===                                                              |   5%
  |                                                                       
  |====                                                             |   6%
  |                                                                       
  |====                                                             |   7%
  |                                                                       
  |=====                                                            |   7%
  |                                                                       
  |=====                                                            |   8%
  |                                                                       
  |======                                                           |   8%
  |                                                                       
  |======                                                           |   9%
  |                                                                       
  |======                                                           |  10%
  |                                                                       
  |=======                                                          |  10%
  |                                                                       
  |=======                                                          |  11%
  |                                                                       
  |========                                                         |  12%
  |                                                                       
  |========                                                         |  13%
  |                                                                       
  |=========                                                        |  13%
  |                                                                       
  |=========                                                        |  14%
  |                                                                       
  |==========                                                       |  15%
  |                                                                       
  |==========                                                       |  16%
  |                                                                       
  |===========                                                      |  17%
  |                                                                       
  |============                                                     |  18%
  |                                                                       
  |============                                                     |  19%
  |                                                                       
  |=============                                                    |  19%
  |                                                                       
  |=============                                                    |  20%
  |                                                                       
  |=============                                                    |  21%
  |                                                                       
  |==============                                                   |  21%
  |                                                                       
  |==============                                                   |  22%
  |                                                                       
  |===============                                                  |  22%
  |                                                                       
  |===============                                                  |  23%
  |                                                                       
  |===============                                                  |  24%
  |                                                                       
  |================                                                 |  24%
  |                                                                       
  |================                                                 |  25%
  |                                                                       
  |=================                                                |  25%
  |                                                                       
  |=================                                                |  26%
  |                                                                       
  |=================                                                |  27%
  |                                                                       
  |==================                                               |  27%
  |                                                                       
  |==================                                               |  28%
  |                                                                       
  |===================                                              |  29%
  |                                                                       
  |===================                                              |  30%
  |                                                                       
  |====================                                             |  30%
  |                                                                       
  |====================                                             |  31%
  |                                                                       
  |=====================                                            |  32%
  |                                                                       
  |=====================                                            |  33%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |======================                                           |  34%
  |                                                                       
  |=======================                                          |  35%
  |                                                                       
  |=======================                                          |  36%
  |                                                                       
  |========================                                         |  36%
  |                                                                       
  |========================                                         |  37%
  |                                                                       
  |=========================                                        |  38%
  |                                                                       
  |=========================                                        |  39%
  |                                                                       
  |==========================                                       |  39%
  |                                                                       
  |==========================                                       |  40%
  |                                                                       
  |==========================                                       |  41%
  |                                                                       
  |===========================                                      |  41%
  |                                                                       
  |===========================                                      |  42%
  |                                                                       
  |============================                                     |  42%
  |                                                                       
  |============================                                     |  43%
  |                                                                       
  |============================                                     |  44%
  |                                                                       
  |=============================                                    |  44%
  |                                                                       
  |=============================                                    |  45%
  |                                                                       
  |==============================                                   |  46%
  |                                                                       
  |==============================                                   |  47%
  |                                                                       
  |===============================                                  |  47%
  |                                                                       
  |===============================                                  |  48%
  |                                                                       
  |================================                                 |  49%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================                                |  50%
  |                                                                       
  |=================================                                |  51%
  |                                                                       
  |==================================                               |  52%
  |                                                                       
  |==================================                               |  53%
  |                                                                       
  |===================================                              |  53%
  |                                                                       
  |===================================                              |  54%
  |                                                                       
  |====================================                             |  55%
  |                                                                       
  |====================================                             |  56%
  |                                                                       
  |=====================================                            |  56%
  |                                                                       
  |=====================================                            |  57%
  |                                                                       
  |=====================================                            |  58%
  |                                                                       
  |======================================                           |  58%
  |                                                                       
  |======================================                           |  59%
  |                                                                       
  |=======================================                          |  59%
  |                                                                       
  |=======================================                          |  60%
  |                                                                       
  |=======================================                          |  61%
  |                                                                       
  |========================================                         |  61%
  |                                                                       
  |========================================                         |  62%
  |                                                                       
  |=========================================                        |  63%
  |                                                                       
  |=========================================                        |  64%
  |                                                                       
  |==========================================                       |  64%
  |                                                                       
  |==========================================                       |  65%
  |                                                                       
  |===========================================                      |  66%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |============================================                     |  67%
  |                                                                       
  |============================================                     |  68%
  |                                                                       
  |=============================================                    |  69%
  |                                                                       
  |=============================================                    |  70%
  |                                                                       
  |==============================================                   |  70%
  |                                                                       
  |==============================================                   |  71%
  |                                                                       
  |===============================================                  |  72%
  |                                                                       
  |===============================================                  |  73%
  |                                                                       
  |================================================                 |  73%
  |                                                                       
  |================================================                 |  74%
  |                                                                       
  |================================================                 |  75%
  |                                                                       
  |=================================================                |  75%
  |                                                                       
  |=================================================                |  76%
  |                                                                       
  |==================================================               |  76%
  |                                                                       
  |==================================================               |  77%
  |                                                                       
  |==================================================               |  78%
  |                                                                       
  |===================================================              |  78%
  |                                                                       
  |===================================================              |  79%
  |                                                                       
  |====================================================             |  79%
  |                                                                       
  |====================================================             |  80%
  |                                                                       
  |====================================================             |  81%
  |                                                                       
  |=====================================================            |  81%
  |                                                                       
  |=====================================================            |  82%
  |                                                                       
  |======================================================           |  83%
  |                                                                       
  |=======================================================          |  84%
  |                                                                       
  |=======================================================          |  85%
  |                                                                       
  |========================================================         |  86%
  |                                                                       
  |========================================================         |  87%
  |                                                                       
  |=========================================================        |  87%
  |                                                                       
  |=========================================================        |  88%
  |                                                                       
  |==========================================================       |  89%
  |                                                                       
  |==========================================================       |  90%
  |                                                                       
  |===========================================================      |  90%
  |                                                                       
  |===========================================================      |  91%
  |                                                                       
  |===========================================================      |  92%
  |                                                                       
  |============================================================     |  92%
  |                                                                       
  |============================================================     |  93%
  |                                                                       
  |=============================================================    |  93%
  |                                                                       
  |=============================================================    |  94%
  |                                                                       
  |==============================================================   |  95%
  |                                                                       
  |==============================================================   |  96%
  |                                                                       
  |===============================================================  |  96%
  |                                                                       
  |===============================================================  |  97%
  |                                                                       
  |================================================================ |  98%
  |                                                                       
  |================================================================ |  99%
  |                                                                       
  |=================================================================| 100%
```

```r
world_dat <- get_ensemble_temp(world, "annualavg", 2080, 2100)
## Subset data to just one scenario, and one percentile
world_dat <- subset(world_dat, world_dat$scenario == "a2")
world_dat <- subset(world_dat, world_dat$percentile == 50)

world_map <- climate_map(world_map_df, world_dat, return_map = T)
world_map + scale_fill_continuous("Temperature \n at 2100", low = "yellow", 
    high = "red") + theme_bw()
```

![plot of chunk worldmap](figure/worldmap.png) 

<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rWBclimate}
-->
Introduction 
========================================================
rWBclimate is an R interface for the World Bank climate data used in the World Bank [climate knowledge portal](http://sdwebx.worldbank.org/climateportal/index.cfm).  

Installation
------
Right now the package is only installable from github with [devtools](http://cran.r-project.org/web/packages/devtools/index.html):

```R
require(rWBclimate)
```
Package description
----

`rWBclimate` provides access to three different classes of climate data at two different spatial scales.  The three different classes of data are GCM model , ensemble and historical data.  Each data class will let you download two different four different types for two different variables.  The two variables are either precipitation expressed in millimeters or temperature as degrees celcius, and for each variable you can download your data in one of four types. The data is fully described below along with examples.

Data classes
---
*__Model data__*

Almost all model data in the Climate Data API are derived from 15 global circulation models (GCMs) used by the Intergovernmental Panel on Climate Change (IPCC) 4th Assessment Reports. The models simulate the response of the global climate system to increasing greenhouse gas concentrations. The data in the Climate Data API have been aggregated to both the country and basin levels, as explained below. Note these data are modeled estimates of temperature and precipitation changes in different time periods under different GCMs and scenarios. They include changes for future time periods and also as "backcasting" (model representations of the past) set for past time periods. The latter should not be confused with any instrumental or observed data. There is a specific dataset with historical measured climate data as well.

**Data types**

|Type|Description|
|----|----|
|Monthly average|The monthly average for all 12 months for a given time period|
|Annual average|a single average for a given time period|
|Monthly anomaly|Average monthly change (anomaly).  The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|
|Annual anomaly|Average annual change (anomaly). The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|

**Data time scales**

Climate model data is only available as averages for 20 year chunks.  The package will automatically convert any start and end data into valid API calls and will return all data between the given start and end date.  The following time periods are available


|Past|     |Future|    |
|----|-----|------|----|
|*start*|  *end*  |  *start*  |*end*|
|1920  | 1939|  2020 | 2039 |
|1940 |  1959|  2040 |2059|
|1960 |   1979|  2060 |2079 |
|1980  | 1999| 2080 | 2099 |

**Data spatial scales**
Data is available at two spatial scales.  The first is country level. You can download data for any country in the world using a valid [ISO 3 letter country code](http://userpage.chemie.fu-berlin.de/diverse/doc/ISO_3166.html). Alternatively you can download data at the basin network for a slightly higher resolution represented as a number 1-468.  This is based on level 2 boundaries from [waterbase.org](http://www.waterbase.org), or you can view a [full map of all the available basins.](http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)


**_Downloading GCM Model Data_**

Model data is downloaded for two different scenarios, the [A2 and B1](http://en.wikipedia.org/wiki/Special_Report_on_Emissions_Scenarios). Generally A2 scenarios are where there is little difference between the future and now and B1 is a more ecologically friendly world with greater decrease in emissions. Both the A2 and the B1 will be downloaded for 15 different GCM models listed in the table below:


|Name in output|Model name|
|--------------|----------|
|bccr_bcm2_0 |[BCM 2.0](http://www-pcmdi.llnl.gov/ipcc/model_documentation/BCCR_BCM2.0.htm)|
|csiro_mk3_5|[CSIRO Mark 3.5](http://www.cawcr.gov.au/publications/technicalreports/CTR_021.pdf)|
|ingv_echam4|[ECHAM 4.6](http://www.bo.ingv.it/)|
|cccma_cgcm3_1|[CGCM 3.1 (T47)](http://www.ec.gc.ca/ccmac-cccma/default.asp?lang=En)|
|cnrm_cm3|[CNRM CM3](http://www.cnrm.meteo.fr/scenario2004/indexenglish.html)|
|gfdl_cm2_0|[GFDL CM2.0](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|gfdl_cm2_1|[GFDL CM2.1](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|ipsl_cm4|[IPSL-CM4](http://mc2.ipsl.jussieu.fr/simules.html)|
|microc3_2_medres|[MIROC 3.2 (medres)](https://esg.llnl.gov:8443/metadata/browseCatalog.do?uri=http://esgcet.llnl.gov/metadata/pcmdi/ipcc/thredds/miroc3_2_medres.sresb1/pcmdi.ipcc4.miroc3_2_medres.sresb1.thredds)|
|miub_echo_g|[ECHO-G](http://www-pcmdi.llnl.gov/projects/modeldoc/cmip/echo-g_tbls.html)|
|mpi_echam5|[ECHAM5/MPI-OM](http://www.mpimet.mpg.de/en/science/models/echam.html)|
|mri_cgcm2_3_2a|[MRI-CGCM2.3.2](http://www.mri-jma.go.jp/Welcome.html)|
|inmcm3_0|[INMCM3.0](http://www.ipcc-data.org/ar4/model-INM-CM3.html)|
|ukmo_hadcm3|[UKMO HadCM3](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadcm3)|
|ukmo_hadgem1|[UKMO HadGEM1](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadgem1)|

The model data can be downloaded with two main functions:
```R
get_model_temp()  ## Get model temperature data
get_model_precip() ## Get model precipitation data
```
*Example 1: Plotting monthly data from different GCM's* 
Say you want to compare temperature from two different models in the USA to see how they vary.  You can download data for the USA and then subset it to the specific models you're interested in and then plot them.





```r
library(rWBclimate)
library(ggplot2)
usa.dat <- get_model_temp("USA", "mavg", 2080, 2100)
usa.dat.bcc <- usa.dat[usa.dat$gcm == "bccr_bcm2_0", ]
usa.dat.had <- usa.dat[usa.dat$gcm == "ukmo_hadcm3", ]
## Add a unique ID to each for easier plotting
usa.dat.bcc$ID <- paste(usa.dat.bcc$scenario, usa.dat.bcc$gcm, sep = "-")
usa.dat.had$ID <- paste(usa.dat.had$scenario, usa.dat.had$gcm, sep = "-")
plot.df <- rbind(usa.dat.bcc, usa.dat.had)
ggplot(plot.df, aes(x = as.factor(month), y = data, group = ID, colour = gcm, 
    linetype = scenario)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") + 
    xlab("Month") + theme_bw()
```

![plot of chunk getmodeldata](figure/getmodeldata.png)


Subsetting all the data can be a bit tedious.  You could also compare all the models but just for one scenario, the A2.


```r
ggplot(usa.dat[usa.dat$scenario == "a2", ], aes(x = month, y = data, group = gcm, 
    colour = gcm)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") + 
    xlab("Month") + theme_bw()
```

![plot of chunk plotallmodeldata](figure/plotallmodeldata.png) 


*Example 2: Plotting annual data for different countries*

Data can be extracted from countries or basins submitted as vectors. Here we will plot the expected temperature anomaly for each 20 year period over a baseline control period of 1961-2000.  These countries chosen span the north to south pole.  It's clear from the plot that the northern most countries (US and Canada) have the biggest anomaly, and Belize, the most equatorial country, has the smallest anomaly.

```r
country.list <- c("CAN", "USA", "MEX", "BLZ", "COL", "PER", "BOL", "ARG")
country.dat <- get_model_temp(country.list, "annualanom", 2010, 2100)
# Subset data
country.dat.bcc <- country.dat[country.dat$gcm == "bccr_bcm2_0", ]
## Exclude A2 scenario
country.dat.bcc <- subset(country.dat.bcc, country.dat.bcc$scenario != "a2")
ggplot(country.dat.bcc, aes(x = fromYear, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Temperature anomaly over baseline") + 
    theme_bw()
```

![plot of chunk annualdata](figure/annualdata.png) 



**_Downloading Ensemble Model Data_**

Getting raw model data is really useful for comparing scenarios or different GCM's.  However many users might be more interested in aggregated measurements. If so you can download ensemble climate data.  This will give access to all 15 GCM's combined and queries will return the 10th, 50th and 90th quantiles for the ensemble of all GCM's.  All queries are constructed the same as raw model data above with the same data types, spatial scales and time scales.

*Example 1: Comparing model quantiles*

Let's look at monthly precipitation predictions for Indonesia for the period of 2080-2100.  We'll create a plot of the precipitation expected in the future for the two different scenarios and plot them with different colours and dashed lines to indicate quantiles.


```r
idn.dat <- get_ensemble_precip("IDN", "mavg", 2080, 2100)
# Set line types
ltype <- rep(1, dim(idn.dat)[1])
ltype[idn.dat$percentile != 50] <- 2
idn.dat$ltype <- ltype
# Create uniqueIDs
idn.dat$uid <- paste(idn.dat$scenario, idn.dat$percentile, sep = "-")
ggplot(idn.dat, aes(x = as.factor(month), y = data, group = uid, colour = scenario, 
    linetype = as.factor(ltype))) + geom_point() + geom_path() + xlab("Month") + 
    ylab("Rain in mm") + theme_bw()
```

![plot of chunk comparing quantiles](figure/comparing_quantiles.png) 


*Example 2: Ensemble statistics*

You can also download 13 different ensemble statistics about basins or countries aside from the raw data as above.  These derived statistics are often given relative to a control period, 1961-2000.  The time periods however are only for 2046-2065, or 2081-2100, not the standard 20 year intervals available for raw model data.  Below is a full table of available statistics.

|Parameter name|Description|Units|
|---------------|-------|-------|
|*tmin_means*|Average daily minimum temperature|degrees Celsius|
|*tmax_means*|Average daily maximum temperature|degrees Celsius|
|*tmax_days90th*|Number of days with maximum temperature above the control period's 90th percentile (hot days)|days|
|*tmin_days90th*|Number of days with minimum  temperature above the control period's 90th percentile (warm nights)|days|
|*tmax_days10th*|Number of days with maximum temperature below the control period's 10th percentile (cool days)|days|
|*tmin_days10th*|Number of days with minimum  temperature below the control period's 10th percentile (cold nights)|days|
|*tmin_days0*|Number of days with minimum  temperature below 0 degrees Celsius|days|
|*ppt_days*|Number of days with precipitation greater than 0.2 mm|days|
|*ppt_days2*|Number of days with precipitation greater than 2 mm|days|
|*ppt_days10*|Number of days with precipitation greater than 10 mm|days|
|*ppt_days90th*|Number of days with precipitation greater than the control period's 90th percentile|days|
|*ppt_dryspell*|Average number of days between precipitation events|days|
|*ppt_means*|Average daily precipitation|mm|

Similar to our previous example where we looked at temperature anomaly along a latitudinal gradient, we can examine more complex ensemble statistics.  Here we examine the average minimum daily temperature in Scavavadian countries.  Also given that only two time periods are available, there is no need to enter in a specific year range

```r
country.list <- c("ISL", "FIN", "NOR", "SWE")
country.dat <- get_ensemble_stats(country.list, "mavg", "tmin_means")
####### Subset data Exclude A2 scenario
country.dat.b1 <- subset(country.dat, country.dat$scenario == "b1")
# choose just one percentile
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$percentile == 50)
# get just one year period
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$fromYear == 2081)


ggplot(country.dat.b1, aes(x = month, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Average daily minimum temperature") + 
    theme_bw() + xlab("Month")
```

![plot of chunk enesmble data](figure/enesmble_data.png) 



**_Downloading Historical Data_**

It is possible to extract historical data from GCM queries, but this is acutally model backcasting output, not true historical records.  The climate data api can download data at the same spatial scales as all other requests (country or basin), but the time scale differs. You can download data at monthly, yearly or decadanal time scales.  Monthly data is actually the mean monthly temperature or precipitation value for each month from 1901-2009 for countries, or 1960-2009 for basins.  You can indicate the type of data you want in the call by setting the `time_scale` parameter to *month*,*year* or *decade*.

*Example 1: Downloading monthly data*

You can download historical precipitation for Belize, Colombia, Peru and Bolivia.  



```r
country.list <- c("BLZ", "COL", "PER", "BOL")
country.dat <- get_historical_precip(country.list, "month")

ggplot(country.dat, aes(x = month, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Average historical precipitation (mm)") + 
    theme_bw() + xlab("Month")
```

![plot of chunk historicalmonth](figure/historicalmonth.png) 


*Example 2: Downloading annual data*

Another use of historical data is to look at increases in temperature over time.

```r
country.list <- c("USA", "MEX", "CAN", "BLZ")
country.dat <- get_historical_temp(country.list, "year")

ggplot(country.dat, aes(x = year, y = data, group = locator)) + geom_point() + 
    geom_path() + ylab("Average annual temperature of Canada") + theme_bw() + 
    xlab("Year") + stat_smooth(se = F, colour = "black") + facet_wrap(~locator, 
    scale = "free")
```

![plot of chunk historicalyear](figure/historicalyear.png) 


**_Mapping climate Data_**

Data can be mapped in [ggplot2](http://docs.ggplot2.org/current/) by downloading KML files from the climate database, creating dataframes and plotting using ggplot2. Maps are available at the basin or country spatial scale.  Because KML files can be large files are downloaded and locally cached in a directory you must set using the kmlpath option.  You can set to any directory as follows: `options(kmlpath = <yourpath>)`.  KML files will be stored there and only downloaded again if they aren't found locally.  If you wish to get rid of them, you will need to manually delete them.  `rWBclimate` allows you to download maps and plot them as ggplot maps, or it can automatically create maps of climate data.  For your convenience we've created basin and country vectors for all six continents *(Antartica has no data in the database)*.  These data files are automatically loaded with the package and can be accessed for basins as *NoAm_basin*, *SoAm_basin*,*Asia_basin*,*Eur_basin*,*Africa_basin*,and *Oceana_basin*. Countries can be similarly accessed as *NoAm_country*, *SoAm_country*,*Asia_country*,*Eur_country*,*Africa_country*,and *Oceana_country*.  If you are plotting by continent, it can take some time to first download all the KML files.

*Example 1: Creating map dataframe and plotting*

Creating map data frames is straightforward, simply provide a list of valid country codes to the function `create_map_df`

```r
# Set the kmlpath option
options(kmlpath = "~/kmltemp")
## Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)
```

```
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=                                                                |   1%
  |                                                                       
  |=                                                                |   2%
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |===                                                              |   4%
  |                                                                       
  |====                                                             |   6%
  |                                                                       
  |====                                                             |   7%
  |                                                                       
  |=====                                                            |   8%
  |                                                                       
  |======                                                           |   9%
  |                                                                       
  |======                                                           |  10%
  |                                                                       
  |=======                                                          |  11%
  |                                                                       
  |========                                                         |  12%
  |                                                                       
  |=========                                                        |  13%
  |                                                                       
  |=========                                                        |  14%
  |                                                                       
  |==========                                                       |  16%
  |                                                                       
  |===========                                                      |  17%
  |                                                                       
  |============                                                     |  18%
  |                                                                       
  |============                                                     |  19%
  |                                                                       
  |=============                                                    |  20%
  |                                                                       
  |==============                                                   |  21%
  |                                                                       
  |==============                                                   |  22%
  |                                                                       
  |===============                                                  |  23%
  |                                                                       
  |================                                                 |  24%
  |                                                                       
  |=================                                                |  26%
  |                                                                       
  |=================                                                |  27%
  |                                                                       
  |==================                                               |  28%
  |                                                                       
  |===================                                              |  29%
  |                                                                       
  |====================                                             |  30%
  |                                                                       
  |====================                                             |  31%
  |                                                                       
  |=====================                                            |  32%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |======================                                           |  34%
  |                                                                       
  |=======================                                          |  36%
  |                                                                       
  |========================                                         |  37%
  |                                                                       
  |=========================                                        |  38%
  |                                                                       
  |=========================                                        |  39%
  |                                                                       
  |==========================                                       |  40%
  |                                                                       
  |===========================                                      |  41%
  |                                                                       
  |===========================                                      |  42%
  |                                                                       
  |============================                                     |  43%
  |                                                                       
  |=============================                                    |  44%
  |                                                                       
  |==============================                                   |  46%
  |                                                                       
  |==============================                                   |  47%
  |                                                                       
  |===============================                                  |  48%
  |                                                                       
  |================================                                 |  49%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================                                |  51%
  |                                                                       
  |==================================                               |  52%
  |                                                                       
  |===================================                              |  53%
  |                                                                       
  |===================================                              |  54%
  |                                                                       
  |====================================                             |  56%
  |                                                                       
  |=====================================                            |  57%
  |                                                                       
  |======================================                           |  58%
  |                                                                       
  |======================================                           |  59%
  |                                                                       
  |=======================================                          |  60%
  |                                                                       
  |========================================                         |  61%
  |                                                                       
  |========================================                         |  62%
  |                                                                       
  |=========================================                        |  63%
  |                                                                       
  |==========================================                       |  64%
  |                                                                       
  |===========================================                      |  66%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |============================================                     |  68%
  |                                                                       
  |=============================================                    |  69%
  |                                                                       
  |==============================================                   |  70%
  |                                                                       
  |==============================================                   |  71%
  |                                                                       
  |===============================================                  |  72%
  |                                                                       
  |================================================                 |  73%
  |                                                                       
  |================================================                 |  74%
  |                                                                       
  |=================================================                |  76%
  |                                                                       
  |==================================================               |  77%
  |                                                                       
  |===================================================              |  78%
  |                                                                       
  |===================================================              |  79%
  |                                                                       
  |====================================================             |  80%
  |                                                                       
  |=====================================================            |  81%
  |                                                                       
  |=====================================================            |  82%
  |                                                                       
  |======================================================           |  83%
  |                                                                       
  |=======================================================          |  84%
  |                                                                       
  |========================================================         |  86%
  |                                                                       
  |========================================================         |  87%
  |                                                                       
  |=========================================================        |  88%
  |                                                                       
  |==========================================================       |  89%
  |                                                                       
  |==========================================================       |  90%
  |                                                                       
  |===========================================================      |  91%
  |                                                                       
  |============================================================     |  92%
  |                                                                       
  |=============================================================    |  93%
  |                                                                       
  |=============================================================    |  94%
  |                                                                       
  |==============================================================   |  96%
  |                                                                       
  |===============================================================  |  97%
  |                                                                       
  |================================================================ |  98%
  |                                                                       
  |================================================================ |  99%
  |                                                                       
  |=================================================================| 100%
```

```r
ggplot(af_basin, aes(x = long, y = lat, group = group)) + geom_polygon() + theme_bw()
```

![plot of chunk map_plot](figure/map_plot.png) 


*Example 2: Mapping climate data*

In order to map climate data you need to have a single point of data for each spatial polygon(country or basin) in your map dataframe.  We've already created the Africa basin map dataframe, so now you need to get some data for each basin.  


```r
af_basin_dat <- get_ensemble_temp(Africa_basin, "annualanom", 2080, 2100)
## Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat, af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat, af_basin_dat$percentile == 50)
```


Now that we have both the map dataframe and the data we can bind them together with the function `climate_map()`.  The function has two options, it can return a `ggplot2` map or it can return a dataframe that you can plot easily with `ggplot2`.  This is useful because it means that you could bind together multiple dataframes for bigger plots, or perhaps add new data yourself


```r

af_map <- climate_map(af_basin, af_basin_dat, return_map = T)
af_map + scale_fill_continuous("Temperature \n anomaly", low = "yellow", high = "red") + 
    theme_bw()
```

![plot of chunk climatemap](figure/climatemap.png) 



*Example 3: Creating a temperature map of the world*


```r
options(kmlpath = "~/kmltemp")
### Combine all country vectors

world <- c(NoAm_country, SoAm_country, Eur_country, Asia_country, Africa_country, 
    Oceana_country)
world_map_df <- create_map_df(world)
```

```
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=                                                                |   1%
  |                                                                       
  |=                                                                |   2%
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |==                                                               |   4%
  |                                                                       
  |===                                                              |   4%
  |                                                                       
  |===                                                              |   5%
  |                                                                       
  |====                                                             |   6%
  |                                                                       
  |====                                                             |   7%
  |                                                                       
  |=====                                                            |   7%
  |                                                                       
  |=====                                                            |   8%
  |                                                                       
  |======                                                           |   8%
  |                                                                       
  |======                                                           |   9%
  |                                                                       
  |======                                                           |  10%
  |                                                                       
  |=======                                                          |  10%
  |                                                                       
  |=======                                                          |  11%
  |                                                                       
  |========                                                         |  12%
  |                                                                       
  |========                                                         |  13%
  |                                                                       
  |=========                                                        |  13%
  |                                                                       
  |=========                                                        |  14%
  |                                                                       
  |==========                                                       |  15%
  |                                                                       
  |==========                                                       |  16%
  |                                                                       
  |===========                                                      |  17%
  |                                                                       
  |============                                                     |  18%
  |                                                                       
  |============                                                     |  19%
  |                                                                       
  |=============                                                    |  19%
  |                                                                       
  |=============                                                    |  20%
  |                                                                       
  |=============                                                    |  21%
  |                                                                       
  |==============                                                   |  21%
  |                                                                       
  |==============                                                   |  22%
  |                                                                       
  |===============                                                  |  22%
  |                                                                       
  |===============                                                  |  23%
  |                                                                       
  |===============                                                  |  24%
  |                                                                       
  |================                                                 |  24%
  |                                                                       
  |================                                                 |  25%
  |                                                                       
  |=================                                                |  25%
  |                                                                       
  |=================                                                |  26%
  |                                                                       
  |=================                                                |  27%
  |                                                                       
  |==================                                               |  27%
  |                                                                       
  |==================                                               |  28%
  |                                                                       
  |===================                                              |  29%
  |                                                                       
  |===================                                              |  30%
  |                                                                       
  |====================                                             |  30%
  |                                                                       
  |====================                                             |  31%
  |                                                                       
  |=====================                                            |  32%
  |                                                                       
  |=====================                                            |  33%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |======================                                           |  34%
  |                                                                       
  |=======================                                          |  35%
  |                                                                       
  |=======================                                          |  36%
  |                                                                       
  |========================                                         |  36%
  |                                                                       
  |========================                                         |  37%
  |                                                                       
  |=========================                                        |  38%
  |                                                                       
  |=========================                                        |  39%
  |                                                                       
  |==========================                                       |  39%
  |                                                                       
  |==========================                                       |  40%
  |                                                                       
  |==========================                                       |  41%
  |                                                                       
  |===========================                                      |  41%
  |                                                                       
  |===========================                                      |  42%
  |                                                                       
  |============================                                     |  42%
  |                                                                       
  |============================                                     |  43%
  |                                                                       
  |============================                                     |  44%
  |                                                                       
  |=============================                                    |  44%
  |                                                                       
  |=============================                                    |  45%
  |                                                                       
  |==============================                                   |  46%
  |                                                                       
  |==============================                                   |  47%
  |                                                                       
  |===============================                                  |  47%
  |                                                                       
  |===============================                                  |  48%
  |                                                                       
  |================================                                 |  49%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================                                |  50%
  |                                                                       
  |=================================                                |  51%
  |                                                                       
  |==================================                               |  52%
  |                                                                       
  |==================================                               |  53%
  |                                                                       
  |===================================                              |  53%
  |                                                                       
  |===================================                              |  54%
  |                                                                       
  |====================================                             |  55%
  |                                                                       
  |====================================                             |  56%
  |                                                                       
  |=====================================                            |  56%
  |                                                                       
  |=====================================                            |  57%
  |                                                                       
  |=====================================                            |  58%
  |                                                                       
  |======================================                           |  58%
  |                                                                       
  |======================================                           |  59%
  |                                                                       
  |=======================================                          |  59%
  |                                                                       
  |=======================================                          |  60%
  |                                                                       
  |=======================================                          |  61%
  |                                                                       
  |========================================                         |  61%
  |                                                                       
  |========================================                         |  62%
  |                                                                       
  |=========================================                        |  63%
  |                                                                       
  |=========================================                        |  64%
  |                                                                       
  |==========================================                       |  64%
  |                                                                       
  |==========================================                       |  65%
  |                                                                       
  |===========================================                      |  66%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |============================================                     |  67%
  |                                                                       
  |============================================                     |  68%
  |                                                                       
  |=============================================                    |  69%
  |                                                                       
  |=============================================                    |  70%
  |                                                                       
  |==============================================                   |  70%
  |                                                                       
  |==============================================                   |  71%
  |                                                                       
  |===============================================                  |  72%
  |                                                                       
  |===============================================                  |  73%
  |                                                                       
  |================================================                 |  73%
  |                                                                       
  |================================================                 |  74%
  |                                                                       
  |================================================                 |  75%
  |                                                                       
  |=================================================                |  75%
  |                                                                       
  |=================================================                |  76%
  |                                                                       
  |==================================================               |  76%
  |                                                                       
  |==================================================               |  77%
  |                                                                       
  |==================================================               |  78%
  |                                                                       
  |===================================================              |  78%
  |                                                                       
  |===================================================              |  79%
  |                                                                       
  |====================================================             |  79%
  |                                                                       
  |====================================================             |  80%
  |                                                                       
  |====================================================             |  81%
  |                                                                       
  |=====================================================            |  81%
  |                                                                       
  |=====================================================            |  82%
  |                                                                       
  |======================================================           |  83%
  |                                                                       
  |=======================================================          |  84%
  |                                                                       
  |=======================================================          |  85%
  |                                                                       
  |========================================================         |  86%
  |                                                                       
  |========================================================         |  87%
  |                                                                       
  |=========================================================        |  87%
  |                                                                       
  |=========================================================        |  88%
  |                                                                       
  |==========================================================       |  89%
  |                                                                       
  |==========================================================       |  90%
  |                                                                       
  |===========================================================      |  90%
  |                                                                       
  |===========================================================      |  91%
  |                                                                       
  |===========================================================      |  92%
  |                                                                       
  |============================================================     |  92%
  |                                                                       
  |============================================================     |  93%
  |                                                                       
  |=============================================================    |  93%
  |                                                                       
  |=============================================================    |  94%
  |                                                                       
  |==============================================================   |  95%
  |                                                                       
  |==============================================================   |  96%
  |                                                                       
  |===============================================================  |  96%
  |                                                                       
  |===============================================================  |  97%
  |                                                                       
  |================================================================ |  98%
  |                                                                       
  |================================================================ |  99%
  |                                                                       
  |=================================================================| 100%
```

```r
world_dat <- get_ensemble_temp(world, "annualavg", 2080, 2100)
## Subset data to just one scenario, and one percentile
world_dat <- subset(world_dat, world_dat$scenario == "a2")
world_dat <- subset(world_dat, world_dat$percentile == 50)

world_map <- climate_map(world_map_df, world_dat, return_map = T)
world_map + scale_fill_continuous("Temperature \n at 2100", low = "yellow", 
    high = "red") + theme_bw()
```

![plot of chunk worldmap](figure/worldmap.png) 

`rWBclimate`
==========

[![Build Status](https://api.travis-ci.org/ropensci/rWBclimate.png)](https://travis-ci.org/ropensci/rWBclimate)
[![Build status](https://ci.appveyor.com/api/projects/status/28njj5uw980frlpu/branch/master)](https://ci.appveyor.com/project/sckott/rwbclimate/branch/master)
[![codecov.io](https://codecov.io/github/ropensci/rWBclimate/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rWBclimate?branch=master)

Introduction
========================================================
rWBclimate is an R interface for the World Bank climate data used in the World Bank [climate knowledge portal](http://sdwebx.worldbank.org/climateportal/index.cfm).

Installation
------
Right now the package is only installable from github with [devtools](http://cran.r-project.org/web/packages/devtools/index.html):

```R
require(rWBclimate)
```
Package description
----

`rWBclimate` provides access to three different classes of climate data at two different spatial scales.  The three different classes of data are GCM model , ensemble and historical data.  Each data class will let you download two different four different types for two different variables.  The two variables are either precipitation expressed in millimeters or temperature as degrees celcius, and for each variable you can download your data in one of four types. The data is fully described below along with examples.

Data classes
---
*__Model data__*

Almost all model data in the Climate Data API are derived from 15 global circulation models (GCMs) used by the Intergovernmental Panel on Climate Change (IPCC) 4th Assessment Reports. The models simulate the response of the global climate system to increasing greenhouse gas concentrations. The data in the Climate Data API have been aggregated to both the country and basin levels, as explained below. Note these data are modeled estimates of temperature and precipitation changes in different time periods under different GCMs and scenarios. They include changes for future time periods and also as “backcasting” (model representations of the past) set for past time periods. The latter should not be confused with any instrumental or observed data. There is a specific dataset with historical measured climate data as well.

**Data types**

|Type|Description|
|----|----|
|Monthly average|The monthly average for all 12 months for a given time period|
|Annual average|a single average for a given time period|
|Monthly anomaly|Average monthly change (anomaly).  The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|
|Annual anomaly|Average annual change (anomaly). The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|

**Data time scales**

Climate model data is only available as averages for 20 year chunks.  The package will automatically convert any start and end data into valid API calls and will return all data between the given start and end date.  The following time periods are available


|Past|     |Future|    |
|----|-----|------|----|
|*start*|  *end*  |  *start*  |*end*|
|1920  | 1939|  2020 | 2039 |
|1940 |  1959|  2040 |2059|
|1960 |   1979|  2060 |2079 |
|1980  | 1999| 2080 | 2099 |

**Data spatial scales**
Data is available at two spatial scales.  The first is country level. You can download data for any country in the world using a valid [ISO 3 letter country code](http://userpage.chemie.fu-berlin.de/diverse/doc/ISO_3166.html). Alternatively you can download data at the basin network for a slightly higher resolution represented as a number 1-468.  This is based on level 2 boundaries from [waterbase.org](http://www.waterbase.org), or you can view a [full map of all the available basins.](http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)


**_Downloading GCM Model Data_**

Model data is downloaded for two different scenarios, the [A2 and B1](http://en.wikipedia.org/wiki/Special_Report_on_Emissions_Scenarios). Generally A2 scenarios are where there is little difference between the future and now and B1 is a more ecologically friendly world with greater decrease in emissions. Both the A2 and the B1 will be downloaded for 15 different GCM models listed in the table below:


|Name in output|Model name|
|--------------|----------|
|bccr_bcm2_0 |[BCM 2.0](http://www-pcmdi.llnl.gov/ipcc/model_documentation/BCCR_BCM2.0.htm)|
|csiro_mk3_5|[CSIRO Mark 3.5](http://www.cawcr.gov.au/publications/technicalreports/CTR_021.pdf)|
|ingv_echam4|[ECHAM 4.6](http://www.bo.ingv.it/)|
|cccma_cgcm3_1|[CGCM 3.1 (T47)](http://www.ec.gc.ca/ccmac-cccma/default.asp?lang=En)|
|cnrm_cm3|[CNRM CM3](http://www.cnrm.meteo.fr/scenario2004/indexenglish.html)|
|gfdl_cm2_0|[GFDL CM2.0](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|gfdl_cm2_1|[GFDL CM2.1](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|ipsl_cm4|[IPSL-CM4](http://mc2.ipsl.jussieu.fr/simules.html)|
|microc3_2_medres|[MIROC 3.2 (medres)](https://esg.llnl.gov:8443/metadata/browseCatalog.do?uri=http://esgcet.llnl.gov/metadata/pcmdi/ipcc/thredds/miroc3_2_medres.sresb1/pcmdi.ipcc4.miroc3_2_medres.sresb1.thredds)|
|miub_echo_g|[ECHO-G](http://www-pcmdi.llnl.gov/projects/modeldoc/cmip/echo-g_tbls.html)|
|mpi_echam5|[ECHAM5/MPI-OM](http://www.mpimet.mpg.de/en/science/models/echam.html)|
|mri_cgcm2_3_2a|[MRI-CGCM2.3.2](http://www.mri-jma.go.jp/Welcome.html)|
|inmcm3_0|[INMCM3.0](http://www.ipcc-data.org/ar4/model-INM-CM3.html)|
|ukmo_hadcm3|[UKMO HadCM3](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadcm3)|
|ukmo_hadgem1|[UKMO HadGEM1](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadgem1)|

The model data can be downloaded with two main functions:
```R
get_model_temp()  ## Get model temperature data
get_model_precip() ## Get model precipitation data
```
*Example 1: Plotting monthly data from different GCM's*
Say you want to compare temperature from two different models in the USA to see how they vary.  You can download data for the USA and then subset it to the specific models you're interested in and then plot them.

```{r loadlib, echo=FALSE,warning=FALSE,message=FALSE,results='hide'}
library(ggplot2)
library(rWBclimate)
```

```{r getmodeldata,message=FALSE}
usa.dat <- get_model_temp("USA","mavg",2080,2100)
usa.dat.bcc <- usa.dat[usa.dat$gcm == "bccr_bcm2_0", ]
usa.dat.had <- usa.dat[usa.dat$gcm=="ukmo_hadcm3",]
## Add a unique ID to each for easier plotting
usa.dat.bcc$ID <- paste(usa.dat.bcc$scenario,usa.dat.bcc$gcm,sep="-")
usa.dat.had$ID <- paste(usa.dat.had$scenario,usa.dat.had$gcm,sep="-")
plot.df <- rbind(usa.dat.bcc,usa.dat.had)
ggplot(plot.df,aes(x=as.factor(month),y=data,group=ID,colour=gcm,linetype=scenario))+geom_point()+geom_path()+ylab("Average temperature in degrees C \n between 2080 and 2100") + xlab("Month")+theme_bw()
```

Subsetting all the data can be a bit tedious.  You could also compare all the models but just for one scenario, the A2.

```{r plotallmodeldata}
ggplot(usa.dat[usa.dat$scenario=="a2",],aes(x=month,y=data,group=gcm,colour=gcm))+geom_point()+geom_path()+ylab("Average temperature in degrees C \n between 2080 and 2100") + xlab("Month")+theme_bw()
```

*Example 2: Plotting annual data for different countries*

Data can be extracted from countries or basins submitted as vectors. Here we will plot the expected temperature anomaly for each 20 year period over a baseline control period of 1961-2000.  These countries chosen span the north to south pole.  It's clear from the plot that the northern most countries (US and Canada) have the biggest anomaly, and Belize, the most equatorial country, has the smallest anomaly.
```{r annualdata}
country.list <- c("CAN","USA","MEX","BLZ","COL","PER","BOL","ARG")
country.dat <- get_model_temp(country.list,"annualanom",2010,2100)
#Subset data
country.dat.bcc <- country.dat[country.dat$gcm == "bccr_bcm2_0", ]
##Exclude A2 scenario
country.dat.bcc <- subset(country.dat.bcc,country.dat.bcc$scenario != "a2")
ggplot(country.dat.bcc,aes(x=fromYear,y=data,group=locator,colour=locator))+geom_point()+geom_path()+ylab("Temperature anomaly over baseline")+theme_bw()
```


**_Downloading Ensemble Model Data_**

Getting raw model data is really useful for comparing scenarios or different GCM's.  However many users might be more interested in aggregated measurements. If so you can download ensemble climate data.  This will give access to all 15 GCM's combined and queries will return the 10th, 50th and 90th quantiles for the ensemble of all GCM's.  All queries are constructed the same as raw model data above with the same data types, spatial scales and time scales.

*Example 1: Comparing model quantiles*

Let's look at monthly precipitation predictions for Indonesia for the period of 2080-2100.  We'll create a plot of the precipitation expected in the future for the two different scenarios and plot them with different colours and dashed lines to indicate quantiles.

```{r comparing quantiles}
idn.dat <- get_ensemble_precip("IDN","mavg",2080,2100)
#Set line types
ltype <- rep(1,dim(idn.dat)[1])
ltype[idn.dat$percentile != 50] <- 2
idn.dat$ltype <- ltype
#Create uniqueIDs
idn.dat$uid <- paste(idn.dat$scenario,idn.dat$percentile,sep="-")
ggplot(idn.dat,aes(x=as.factor(month),y=data,group=uid,colour=scenario,linetype=as.factor(ltype)))+geom_point()+geom_path()+xlab("Month")+ylab("Rain in mm")+theme_bw()
```

*Example 2: Ensemble statistics*

You can also download 13 different ensemble statistics about basins or countries aside from the raw data as above.  These derived statistics are often given relative to a control period, 1961-2000.  The time periods however are only for 2046-2065, or 2081-2100, not the standard 20 year intervals available for raw model data.  Below is a full table of available statistics.

|Parameter name|Description|Units|
|---------------|-------|-------|
|*tmin_means*|Average daily minimum temperature|degrees Celsius|
|*tmax_means*|Average daily maximum temperature|degrees Celsius|
|*tmax_days90th*|Number of days with maximum temperature above the control period’s 90th percentile (hot days)|days|
|*tmin_days90th*|Number of days with minimum  temperature above the control period’s 90th percentile (warm nights)|days|
|*tmax_days10th*|Number of days with maximum temperature below the control period’s 10th percentile (cool days)|days|
|*tmin_days10th*|Number of days with minimum  temperature below the control period’s 10th percentile (cold nights)|days|
|*tmin_days0*|Number of days with minimum  temperature below 0 degrees Celsius|days|
|*ppt_days*|Number of days with precipitation greater than 0.2 mm|days|
|*ppt_days2*|Number of days with precipitation greater than 2 mm|days|
|*ppt_days10*|Number of days with precipitation greater than 10 mm|days|
|*ppt_days90th*|Number of days with precipitation greater than the control period's 90th percentile|days|
|*ppt_dryspell*|Average number of days between precipitation events|days|
|*ppt_means*|Average daily precipitation|mm|

Similar to our previous example where we looked at temperature anomaly along a latitudinal gradient, we can examine more complex ensemble statistics.  Here we examine the average minimum daily temperature in Scavavadian countries.  Also given that only two time periods are available, there is no need to enter in a specific year range
```{r enesmble data}
country.list <- c("ISL","FIN","NOR","SWE")
country.dat <- get_ensemble_stats(country.list,"mavg","tmin_means")
####### Subset data
## Exclude A2 scenario
country.dat.b1 <- subset(country.dat,country.dat$scenario == "b1")
# choose just one percentile
country.dat.b1 <- subset(country.dat.b1,country.dat.b1$percentile== 50)
# get just one year period
country.dat.b1 <- subset(country.dat.b1,country.dat.b1$fromYear== 2081)


ggplot(country.dat.b1,aes(x=month,y=data,group=locator,colour=locator))+geom_point()+geom_path()+ylab("Average daily minimum temperature")+theme_bw() + xlab("Month")
```


**_Downloading Historical Data_**

It is possible to extract historical data from GCM queries, but this is acutally model backcasting output, not true historical records.  The climate data api can download data at the same spatial scales as all other requests (country or basin), but the time scale differs. You can download data at monthly, yearly or decadanal time scales.  Monthly data is actually the mean monthly temperature or precipitation value for each month from 1901-2009 for countries, or 1960-2009 for basins.  You can indicate the type of data you want in the call by setting the `time_scale` parameter to *month*,*year* or *decade*.

*Example 1: Downloading monthly data*

You can download historical precipitation for Belize, Colombia, Peru and Bolivia.


```{r historicalmonth}
country.list <- c("BLZ","COL","PER","BOL")
country.dat <- get_historical_precip(country.list,"month")

ggplot(country.dat,aes(x=month,y=data,group=locator,colour=locator))+geom_point()+geom_path()+ylab("Average historical precipitation (mm)")+theme_bw()+xlab("Month")

```

*Example 2: Downloading annual data*

Another use of historical data is to look at increases in temperature over time.
```{r historicalyear, message=FALSE,warning=FALSE}
country.list <- c("USA","MEX","CAN","BLZ")
country.dat <- get_historical_temp(country.list,"year")

ggplot(country.dat,aes(x=year,y=data,group=locator))+geom_point()+geom_path()+ylab("Average annual temperature of Canada")+theme_bw()+xlab("Year")+stat_smooth(se=F,colour="black")+facet_wrap(~locator,scale="free")

```

**_Mapping climate Data_**

Data can be mapped in [ggplot2](http://docs.ggplot2.org/current/) by downloading KML files from the climate database, creating dataframes and plotting using ggplot2. Maps are available at the basin or country spatial scale.  Because KML files can be large files are downloaded and locally cached in a directory you must set using the kmlpath option.  You can set to any directory as follows: `options(kmlpath = <yourpath>)`.  KML files will be stored there and only downloaded again if they aren't found locally.  If you wish to get rid of them, you will need to manually delete them.  `rWBclimate` allows you to download maps and plot them as ggplot maps, or it can automatically create maps of climate data.  For your convenience we've created basin and country vectors for all six continents *(Antartica has no data in the database)*.  These data files are automatically loaded with the package and can be accessed for basins as *NoAm_basin*, *SoAm_basin*,*Asia_basin*,*Eur_basin*,*Africa_basin*,and *Oceana_basin*. Countries can be similarly accessed as *NoAm_country*, *SoAm_country*,*Asia_country*,*Eur_country*,*Africa_country*,and *Oceana_country*.  If you are plotting by continent, it can take some time to first download all the KML files.

*Example 1: Creating map dataframe and plotting*

Creating map data frames is straightforward, simply provide a list of valid country codes to the function `create_map_df`
```{r map_plot,fig.width=5.5, fig.height=7,warning=FALSE,message=FALSE}
#Set the kmlpath option
options(kmlpath = "~/kmltemp")
##Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)
ggplot(af_basin, aes(x=long, y=lat,group=group))+ geom_polygon() + theme_bw()
```

*Example 2: Mapping climate data*

In order to map climate data you need to have a single point of data for each spatial polygon(country or basin) in your map dataframe.  We've already created the Africa basin map dataframe, so now you need to get some data for each basin.

```{r data_to_map}
af_basin_dat <- get_ensemble_temp(Africa_basin,"annualanom",2080,2100)
##  Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat,af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat,af_basin_dat$percentile == 50)
```

Now that we have both the map dataframe and the data we can bind them together with the function `climate_map()`.  The function has two options, it can return a `ggplot2` map or it can return a dataframe that you can plot easily with `ggplot2`.  This is useful because it means that you could bind together multiple dataframes for bigger plots, or perhaps add new data yourself

```{r climatemap,fig.width=5.5, fig.height=7,message=FALSE}

af_map <- climate_map(af_basin,af_basin_dat,return_map = T)
af_map + scale_fill_continuous("Temperature \n anomaly",low="yellow",high = "red") + theme_bw()

```


*Example 3: Creating a temperature map of the world*

```{r worldmap,fig.width=8, fig.height=6,warning=FALSE,message=FALSE}
options(kmlpath = "~/kmltemp")
### Combine all country vectors

world <- c(NoAm_country,SoAm_country,Eur_country,Asia_country,Africa_country,Oceana_country)
world_map_df <- create_map_df(world)
world_dat <- get_ensemble_temp(world,"annualavg",2080,2100)
  ##  Subset data to just one scenario, and one percentile
world_dat <- subset(world_dat,world_dat$scenario == "a2")
world_dat <- subset(world_dat,world_dat$percentile == 50)

world_map <- climate_map(world_map_df,world_dat,return_map = T)
world_map + scale_fill_continuous("Temperature \n at 2100",low="yellow",high = "red") + theme_bw()

```

To cite package ‘rWBclimate’ in publications use:

```coffee

  Edmund Hart (). rWBclimate: A package for accessing World Bank climate data. R package version 0.1.3
  http://github.com/ropensci/rWBclimate

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {rWBclimate: A package for accessing World Bank climate data},
    author = {Edmund Hart},
    note = {R package version 0.1.3},
    url = {http://github.com/ropensci/rWBclimate},
  }
```


[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rWBclimate}
-->
Introduction 
========================================================
rWBclimate is an R interface for the World Bank climate data used in the World Bank [climate knowledge portal](http://sdwebx.worldbank.org/climateportal/index.cfm).  

Installation
------
Right now the package is only installable from github with [devtools](http://cran.r-project.org/web/packages/devtools/index.html):

```R
require(rWBclimate)
```
Package description
----

`rWBclimate` provides access to three different classes of climate data at two different spatial scales.  The three different classes of data are GCM model , ensemble and historical data.  Each data class will let you download two different four different types for two different variables.  The two variables are either precipitation expressed in millimeters or temperature as degrees celcius, and for each variable you can download your data in one of four types. The data is fully described below along with examples.

Data classes
---
*__Model data__*

Almost all model data in the Climate Data API are derived from 15 global circulation models (GCMs) used by the Intergovernmental Panel on Climate Change (IPCC) 4th Assessment Reports. The models simulate the response of the global climate system to increasing greenhouse gas concentrations. The data in the Climate Data API have been aggregated to both the country and basin levels, as explained below. Note these data are modeled estimates of temperature and precipitation changes in different time periods under different GCMs and scenarios. They include changes for future time periods and also as " backcasting" (model representations of the past) set for past time periods. The latter should not be confused with any instrumental or observed data. There is a specific dataset with historical measured climate data as well.

**Data types**

|Type|Description|
|----|----|
|Monthly average|The monthly average for all 12 months for a given time period|
|Annual average|a single average for a given time period|
|Monthly anomaly|Average monthly change (anomaly).  The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|
|Annual anomaly|Average annual change (anomaly). The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|

**Data time scales**

Climate model data is only available as averages for 20 year chunks.  The package will automatically convert any start and end data into valid API calls and will return all data between the given start and end date.  The following time periods are available


|Past|     |Future|    |
|----|-----|------|----|
|*start*|  *end*  |  *start*  |*end*|
|1920  | 1939|  2020 | 2039 |
|1940 |  1959|  2040 |2059|
|1960 |   1979|  2060 |2079 |
|1980  | 1999| 2080 | 2099 |

**Data spatial scales**
Data is available at two spatial scales.  The first is country level. You can download data for any country in the world using a valid [ISO 3 letter country code](http://userpage.chemie.fu-berlin.de/diverse/doc/ISO_3166.html). Alternatively you can download data at the basin network for a slightly higher resolution represented as a number 1-468.  This is based on level 2 boundaries from [waterbase.org](http://www.waterbase.org), or you can view a [full map of all the available basins.](http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)


**_Downloading GCM Model Data_**

Model data is downloaded for two different scenarios, the [A2 and B1](http://en.wikipedia.org/wiki/Special_Report_on_Emissions_Scenarios). Generally A2 scenarios are where there is little difference between the future and now and B1 is a more ecologically friendly world with greater decrease in emissions. Both the A2 and the B1 will be downloaded for 15 different GCM models listed in the table below:


|Name in output|Model name|
|--------------|----------|
|bccr_bcm2_0 |[BCM 2.0](http://www-pcmdi.llnl.gov/ipcc/model_documentation/BCCR_BCM2.0.htm)|
|csiro_mk3_5|[CSIRO Mark 3.5](http://www.cawcr.gov.au/publications/technicalreports/CTR_021.pdf)|
|ingv_echam4|[ECHAM 4.6](http://www.bo.ingv.it/)|
|cccma_cgcm3_1|[CGCM 3.1 (T47)](http://www.ec.gc.ca/ccmac-cccma/default.asp?lang=En)|
|cnrm_cm3|[CNRM CM3](http://www.cnrm.meteo.fr/scenario2004/indexenglish.html)|
|gfdl_cm2_0|[GFDL CM2.0](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|gfdl_cm2_1|[GFDL CM2.1](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|ipsl_cm4|[IPSL-CM4](http://mc2.ipsl.jussieu.fr/simules.html)|
|microc3_2_medres|[MIROC 3.2 (medres)](https://esg.llnl.gov:8443/metadata/browseCatalog.do?uri=http://esgcet.llnl.gov/metadata/pcmdi/ipcc/thredds/miroc3_2_medres.sresb1/pcmdi.ipcc4.miroc3_2_medres.sresb1.thredds)|
|miub_echo_g|[ECHO-G](http://www-pcmdi.llnl.gov/projects/modeldoc/cmip/echo-g_tbls.html)|
|mpi_echam5|[ECHAM5/MPI-OM](http://www.mpimet.mpg.de/en/science/models/echam.html)|
|mri_cgcm2_3_2a|[MRI-CGCM2.3.2](http://www.mri-jma.go.jp/Welcome.html)|
|inmcm3_0|[INMCM3.0](http://www.ipcc-data.org/ar4/model-INM-CM3.html)|
|ukmo_hadcm3|[UKMO HadCM3](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadcm3)|
|ukmo_hadgem1|[UKMO HadGEM1](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadgem1)|

The model data can be downloaded with two main functions:
```R
get_model_temp()  ## Get model temperature data
get_model_precip() ## Get model precipitation data
```
*Example 1: Plotting monthly data from different GCM's* 
Say you want to compare temperature from two different models in the USA to see how they vary.  You can download data for the USA and then subset it to the specific models you're interested in and then plot them.

```{r loadlib, echo=FALSE,warning=FALSE,message=FALSE,results='hide'}
library(ggplot2)
library(rWBclimate)
```

```{r getmodeldata,message=FALSE}
usa.dat <- get_model_temp("USA","mavg",2080,2100)
usa.dat.bcc <- usa.dat[usa.dat$gcm == "bccr_bcm2_0", ]
usa.dat.had <- usa.dat[usa.dat$gcm=="ukmo_hadcm3",]
## Add a unique ID to each for easier plotting
usa.dat.bcc$ID <- paste(usa.dat.bcc$scenario,usa.dat.bcc$gcm,sep="-")
usa.dat.had$ID <- paste(usa.dat.had$scenario,usa.dat.had$gcm,sep="-")
plot.df <- rbind(usa.dat.bcc,usa.dat.had)
ggplot(plot.df,aes(x=as.factor(month),y=data,group=ID,colour=gcm,linetype=scenario))+geom_point()+geom_path()+ylab("Average temperature in degrees C \n between 2080 and 2100") + xlab("Month")+theme_bw()
```

Subsetting all the data can be a bit tedious.  You could also compare all the models but just for one scenario, the A2.

```{r plotallmodeldata}
ggplot(usa.dat[usa.dat$scenario=="a2",],aes(x=month,y=data,group=gcm,colour=gcm))+geom_point()+geom_path()+ylab("Average temperature in degrees C \n between 2080 and 2100") + xlab("Month")+theme_bw()
```

*Example 2: Plotting annual data for different countries*

Data can be extracted from countries or basins submitted as vectors. Here we will plot the expected temperature anomaly for each 20 year period over a baseline control period of 1961-2000.  These countries chosen span the north to south pole.  It's clear from the plot that the northern most countries (US and Canada) have the biggest anomaly, and Belize, the most equatorial country, has the smallest anomaly.
```{r annualdata}
country.list <- c("CAN","USA","MEX","BLZ","COL","PER","BOL","ARG")
country.dat <- get_model_temp(country.list,"annualanom",2010,2100)
#Subset data
country.dat.bcc <- country.dat[country.dat$gcm == "bccr_bcm2_0", ]
##Exclude A2 scenario
country.dat.bcc <- subset(country.dat.bcc,country.dat.bcc$scenario != "a2")
ggplot(country.dat.bcc,aes(x=fromYear,y=data,group=locator,colour=locator))+geom_point()+geom_path()+ylab("Temperature anomaly over baseline")+theme_bw()
```


**_Downloading Ensemble Model Data_**

Getting raw model data is really useful for comparing scenarios or different GCM's.  However many users might be more interested in aggregated measurements. If so you can download ensemble climate data.  This will give access to all 15 GCM's combined and queries will return the 10th, 50th and 90th quantiles for the ensemble of all GCM's.  All queries are constructed the same as raw model data above with the same data types, spatial scales and time scales.

*Example 1: Comparing model quantiles*

Let's look at monthly precipitation predictions for Indonesia for the period of 2080-2100.  We'll create a plot of the precipitation expected in the future for the two different scenarios and plot them with different colours and dashed lines to indicate quantiles.

```{r comparing quantiles}
idn.dat <- get_ensemble_precip("IDN","mavg",2080,2100)
#Set line types
ltype <- rep(1,dim(idn.dat)[1])
ltype[idn.dat$percentile != 50] <- 2 
idn.dat$ltype <- ltype
#Create uniqueIDs
idn.dat$uid <- paste(idn.dat$scenario,idn.dat$percentile,sep="-")
ggplot(idn.dat,aes(x=as.factor(month),y=data,group=uid,colour=scenario,linetype=as.factor(ltype)))+geom_point()+geom_path()+xlab("Month")+ylab("Rain in mm")+theme_bw()
```

*Example 2: Ensemble statistics*

You can also download 13 different ensemble statistics about basins or countries aside from the raw data as above.  These derived statistics are often given relative to a control period, 1961-2000.  The time periods however are only for 2046-2065, or 2081-2100, not the standard 20 year intervals available for raw model data.  Below is a full table of available statistics.

|Parameter name|Description|Units|
|---------------|-------|-------|
|*tmin_means*|Average daily minimum temperature|degrees Celsius|
|*tmax_means*|Average daily maximum temperature|degrees Celsius|
|*tmax_days90th*|Number of days with maximum temperature above the control period's 90th percentile (hot days)|days|
|*tmin_days90th*|Number of days with minimum  temperature above the control period's 90th percentile (warm nights)|days|
|*tmax_days10th*|Number of days with maximum temperature below the control period's 10th percentile (cool days)|days|
|*tmin_days10th*|Number of days with minimum  temperature below the control period's 10th percentile (cold nights)|days|
|*tmin_days0*|Number of days with minimum  temperature below 0 degrees Celsius|days|
|*ppt_days*|Number of days with precipitation greater than 0.2 mm|days|
|*ppt_days2*|Number of days with precipitation greater than 2 mm|days|
|*ppt_days10*|Number of days with precipitation greater than 10 mm|days|
|*ppt_days90th*|Number of days with precipitation greater than the control period's 90th percentile|days|
|*ppt_dryspell*|Average number of days between precipitation events|days|
|*ppt_means*|Average daily precipitation|mm|

Similar to our previous example where we looked at temperature anomaly along a latitudinal gradient, we can examine more complex ensemble statistics.  Here we examine the average minimum daily temperature in Scavavadian countries.  Also given that only two time periods are available, there is no need to enter in a specific year range
```{r enesmble data}
country.list <- c("ISL","FIN","NOR","SWE")
country.dat <- get_ensemble_stats(country.list,"mavg","tmin_means")
####### Subset data
## Exclude A2 scenario
country.dat.b1 <- subset(country.dat,country.dat$scenario == "b1")
# choose just one percentile
country.dat.b1 <- subset(country.dat.b1,country.dat.b1$percentile== 50)
# get just one year period
country.dat.b1 <- subset(country.dat.b1,country.dat.b1$fromYear== 2081)


ggplot(country.dat.b1,aes(x=month,y=data,group=locator,colour=locator))+geom_point()+geom_path()+ylab("Average daily minimum temperature")+theme_bw() + xlab("Month")
```


**_Downloading Historical Data_**

It is possible to extract historical data from GCM queries, but this is acutally model backcasting output, not true historical records.  The climate data api can download data at the same spatial scales as all other requests (country or basin), but the time scale differs. You can download data at monthly, yearly or decadanal time scales.  Monthly data is actually the mean monthly temperature or precipitation value for each month from 1901-2009 for countries, or 1960-2009 for basins.  You can indicate the type of data you want in the call by setting the `time_scale` parameter to *month*,*year* or *decade*.

*Example 1: Downloading monthly data*

You can download historical precipitation for Belize, Colombia, Peru and Bolivia.  


```{r historicalmonth}
country.list <- c("BLZ","COL","PER","BOL")
country.dat <- get_historical_precip(country.list,"month")

ggplot(country.dat,aes(x=month,y=data,group=locator,colour=locator))+geom_point()+geom_path()+ylab("Average historical precipitation (mm)")+theme_bw()+xlab("Month")

```

*Example 2: Downloading annual data*

Another use of historical data is to look at increases in temperature over time.
```{r historicalyear, message=FALSE,warning=FALSE}
country.list <- c("USA","MEX","CAN","BLZ")
country.dat <- get_historical_temp(country.list,"year")

ggplot(country.dat,aes(x=year,y=data,group=locator))+geom_point()+geom_path()+ylab("Average annual temperature of Canada")+theme_bw()+xlab("Year")+stat_smooth(se=F,colour="black")+facet_wrap(~locator,scale="free")

```

**_Mapping climate Data_**

Data can be mapped in [ggplot2](http://docs.ggplot2.org/current/) by downloading KML files from the climate database, creating dataframes and plotting using ggplot2. Maps are available at the basin or country spatial scale.  Because KML files can be large files are downloaded and locally cached in a directory you must set using the kmlpath option.  You can set to any directory as follows: `options(kmlpath = <yourpath>)`.  KML files will be stored there and only downloaded again if they aren't found locally.  If you wish to get rid of them, you will need to manually delete them.  `rWBclimate` allows you to download maps and plot them as ggplot maps, or it can automatically create maps of climate data.  For your convenience we've created basin and country vectors for all six continents *(Antartica has no data in the database)*.  These data files are automatically loaded with the package and can be accessed for basins as *NoAm_basin*, *SoAm_basin*,*Asia_basin*,*Eur_basin*,*Africa_basin*,and *Oceana_basin*. Countries can be similarly accessed as *NoAm_country*, *SoAm_country*,*Asia_country*,*Eur_country*,*Africa_country*,and *Oceana_country*.  If you are plotting by continent, it can take some time to first download all the KML files.

*Example 1: Creating map dataframe and plotting*

Creating map data frames is straightforward, simply provide a list of valid country codes to the function `create_map_df`
```{r map_plot,fig.width=5.5, fig.height=7,warning=FALSE,message=FALSE}
#Set the kmlpath option
options(kmlpath = "~/kmltemp")
##Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)
ggplot(af_basin, aes(x=long, y=lat,group=group))+ geom_polygon() + theme_bw()
```

*Example 2: Mapping climate data*

In order to map climate data you need to have a single point of data for each spatial polygon(country or basin) in your map dataframe.  We've already created the Africa basin map dataframe, so now you need to get some data for each basin.  

```{r data_to_map}
af_basin_dat <- get_ensemble_temp(Africa_basin,"annualanom",2080,2100)
##  Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat,af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat,af_basin_dat$percentile == 50)
```

Now that we have both the map dataframe and the data we can bind them together with the function `climate_map()`.  The function has two options, it can return a `ggplot2` map or it can return a dataframe that you can plot easily with `ggplot2`.  This is useful because it means that you could bind together multiple dataframes for bigger plots, or perhaps add new data yourself

```{r climatemap,fig.width=5.5, fig.height=7,message=FALSE}

af_map <- climate_map(af_basin,af_basin_dat,return_map = T)
af_map + scale_fill_continuous("Temperature \n anomaly",low="yellow",high = "red") + theme_bw()

```


*Example 3: Creating a temperature map of the world*

```{r worldmap,fig.width=8, fig.height=6,warning=FALSE,message=FALSE}
options(kmlpath = "~/kmltemp")
### Combine all country vectors

world <- c(NoAm_country,SoAm_country,Eur_country,Asia_country,Africa_country,Oceana_country)
world_map_df <- create_map_df(world)
world_dat <- get_ensemble_temp(world,"annualavg",2080,2100)
  ##  Subset data to just one scenario, and one percentile
world_dat <- subset(world_dat,world_dat$scenario == "a2")
world_dat <- subset(world_dat,world_dat$percentile == 50)

world_map <- climate_map(world_map_df,world_dat,return_map = T)
world_map + scale_fill_continuous("Temperature \n at 2100",low="yellow",high = "red") + theme_bw()

```
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rWBclimate}
-->
Introduction 
========================================================
rWBclimate is an R interface for the World Bank climate data used in the World Bank [climate knowledge portal](http://sdwebx.worldbank.org/climateportal/index.cfm).  

Installation
------
Right now the package is only installable from github with [devtools](http://cran.r-project.org/web/packages/devtools/index.html):

```R
require(rWBclimate)
```
Package description
----

`rWBclimate` provides access to three different classes of climate data at two different spatial scales.  The three different classes of data are GCM model , ensemble and historical data.  Each data class will let you download two different four different types for two different variables.  The two variables are either precipitation expressed in millimeters or temperature as degrees celcius, and for each variable you can download your data in one of four types. The data is fully described below along with examples.

Data classes
---
*__Model data__*

Almost all model data in the Climate Data API are derived from 15 global circulation models (GCMs) used by the Intergovernmental Panel on Climate Change (IPCC) 4th Assessment Reports. The models simulate the response of the global climate system to increasing greenhouse gas concentrations. The data in the Climate Data API have been aggregated to both the country and basin levels, as explained below. Note these data are modeled estimates of temperature and precipitation changes in different time periods under different GCMs and scenarios. They include changes for future time periods and also as "backcasting" (model representations of the past) set for past time periods. The latter should not be confused with any instrumental or observed data. There is a specific dataset with historical measured climate data as well.

**Data types**

|Type|Description|
|----|----|
|Monthly average|The monthly average for all 12 months for a given time period|
|Annual average|a single average for a given time period|
|Monthly anomaly|Average monthly change (anomaly).  The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|
|Annual anomaly|Average annual change (anomaly). The control period is 1961-1999 for temperature and precipitation variables, and 1961-2000 for derived statistics.|

**Data time scales**

Climate model data is only available as averages for 20 year chunks.  The package will automatically convert any start and end data into valid API calls and will return all data between the given start and end date.  The following time periods are available


|Past|     |Future|    |
|----|-----|------|----|
|*start*|  *end*  |  *start*  |*end*|
|1920  | 1939|  2020 | 2039 |
|1940 |  1959|  2040 |2059|
|1960 |   1979|  2060 |2079 |
|1980  | 1999| 2080 | 2099 |

**Data spatial scales**
Data is available at two spatial scales.  The first is country level. You can download data for any country in the world using a valid [ISO 3 letter country code](http://userpage.chemie.fu-berlin.de/diverse/doc/ISO_3166.html). Alternatively you can download data at the basin network for a slightly higher resolution represented as a number 1-468.  This is based on level 2 boundaries from [waterbase.org](http://www.waterbase.org), or you can view a [full map of all the available basins.](http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)


**_Downloading GCM Model Data_**

Model data is downloaded for two different scenarios, the [A2 and B1](http://en.wikipedia.org/wiki/Special_Report_on_Emissions_Scenarios). Generally A2 scenarios are where there is little difference between the future and now and B1 is a more ecologically friendly world with greater decrease in emissions. Both the A2 and the B1 will be downloaded for 15 different GCM models listed in the table below:


|Name in output|Model name|
|--------------|----------|
|bccr_bcm2_0 |[BCM 2.0](http://www-pcmdi.llnl.gov/ipcc/model_documentation/BCCR_BCM2.0.htm)|
|csiro_mk3_5|[CSIRO Mark 3.5](http://www.cawcr.gov.au/publications/technicalreports/CTR_021.pdf)|
|ingv_echam4|[ECHAM 4.6](http://www.bo.ingv.it/)|
|cccma_cgcm3_1|[CGCM 3.1 (T47)](http://www.ec.gc.ca/ccmac-cccma/default.asp?lang=En)|
|cnrm_cm3|[CNRM CM3](http://www.cnrm.meteo.fr/scenario2004/indexenglish.html)|
|gfdl_cm2_0|[GFDL CM2.0](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|gfdl_cm2_1|[GFDL CM2.1](http://data1.gfdl.noaa.gov/nomads/forms/deccen/CM2.X)|
|ipsl_cm4|[IPSL-CM4](http://mc2.ipsl.jussieu.fr/simules.html)|
|microc3_2_medres|[MIROC 3.2 (medres)](https://esg.llnl.gov:8443/metadata/browseCatalog.do?uri=http://esgcet.llnl.gov/metadata/pcmdi/ipcc/thredds/miroc3_2_medres.sresb1/pcmdi.ipcc4.miroc3_2_medres.sresb1.thredds)|
|miub_echo_g|[ECHO-G](http://www-pcmdi.llnl.gov/projects/modeldoc/cmip/echo-g_tbls.html)|
|mpi_echam5|[ECHAM5/MPI-OM](http://www.mpimet.mpg.de/en/science/models/echam.html)|
|mri_cgcm2_3_2a|[MRI-CGCM2.3.2](http://www.mri-jma.go.jp/Welcome.html)|
|inmcm3_0|[INMCM3.0](http://www.ipcc-data.org/ar4/model-INM-CM3.html)|
|ukmo_hadcm3|[UKMO HadCM3](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadcm3)|
|ukmo_hadgem1|[UKMO HadGEM1](http://www.metoffice.gov.uk/research/modelling-systems/unified-model/climate-models/hadgem1)|

The model data can be downloaded with two main functions:
```R
get_model_temp()  ## Get model temperature data
get_model_precip() ## Get model precipitation data
```
*Example 1: Plotting monthly data from different GCM's* 
Say you want to compare temperature from two different models in the USA to see how they vary.  You can download data for the USA and then subset it to the specific models you're interested in and then plot them.





```r
library(rWBclimate)
library(ggplot2)
usa.dat <- get_model_temp("USA", "mavg", 2080, 2100)
usa.dat.bcc <- usa.dat[usa.dat$gcm == "bccr_bcm2_0", ]
usa.dat.had <- usa.dat[usa.dat$gcm == "ukmo_hadcm3", ]
## Add a unique ID to each for easier plotting
usa.dat.bcc$ID <- paste(usa.dat.bcc$scenario, usa.dat.bcc$gcm, sep = "-")
usa.dat.had$ID <- paste(usa.dat.had$scenario, usa.dat.had$gcm, sep = "-")
plot.df <- rbind(usa.dat.bcc, usa.dat.had)
ggplot(plot.df, aes(x = as.factor(month), y = data, group = ID, colour = gcm, 
    linetype = scenario)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") + 
    xlab("Month") + theme_bw()
```

![plot of chunk getmodeldata](figure/getmodeldata.png)


Subsetting all the data can be a bit tedious.  You could also compare all the models but just for one scenario, the A2.


```r
ggplot(usa.dat[usa.dat$scenario == "a2", ], aes(x = month, y = data, group = gcm, 
    colour = gcm)) + geom_point() + geom_path() + ylab("Average temperature in degrees C \n between 2080 and 2100") + 
    xlab("Month") + theme_bw()
```

![plot of chunk plotallmodeldata](figure/plotallmodeldata.png) 


*Example 2: Plotting annual data for different countries*

Data can be extracted from countries or basins submitted as vectors. Here we will plot the expected temperature anomaly for each 20 year period over a baseline control period of 1961-2000.  These countries chosen span the north to south pole.  It's clear from the plot that the northern most countries (US and Canada) have the biggest anomaly, and Belize, the most equatorial country, has the smallest anomaly.

```r
country.list <- c("CAN", "USA", "MEX", "BLZ", "COL", "PER", "BOL", "ARG")
country.dat <- get_model_temp(country.list, "annualanom", 2010, 2100)
# Subset data
country.dat.bcc <- country.dat[country.dat$gcm == "bccr_bcm2_0", ]
## Exclude A2 scenario
country.dat.bcc <- subset(country.dat.bcc, country.dat.bcc$scenario != "a2")
ggplot(country.dat.bcc, aes(x = fromYear, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Temperature anomaly over baseline") + 
    theme_bw()
```

![plot of chunk annualdata](figure/annualdata.png) 



**_Downloading Ensemble Model Data_**

Getting raw model data is really useful for comparing scenarios or different GCM's.  However many users might be more interested in aggregated measurements. If so you can download ensemble climate data.  This will give access to all 15 GCM's combined and queries will return the 10th, 50th and 90th quantiles for the ensemble of all GCM's.  All queries are constructed the same as raw model data above with the same data types, spatial scales and time scales.

*Example 1: Comparing model quantiles*

Let's look at monthly precipitation predictions for Indonesia for the period of 2080-2100.  We'll create a plot of the precipitation expected in the future for the two different scenarios and plot them with different colours and dashed lines to indicate quantiles.


```r
idn.dat <- get_ensemble_precip("IDN", "mavg", 2080, 2100)
# Set line types
ltype <- rep(1, dim(idn.dat)[1])
ltype[idn.dat$percentile != 50] <- 2
idn.dat$ltype <- ltype
# Create uniqueIDs
idn.dat$uid <- paste(idn.dat$scenario, idn.dat$percentile, sep = "-")
ggplot(idn.dat, aes(x = as.factor(month), y = data, group = uid, colour = scenario, 
    linetype = as.factor(ltype))) + geom_point() + geom_path() + xlab("Month") + 
    ylab("Rain in mm") + theme_bw()
```

![plot of chunk comparing quantiles](figure/comparing_quantiles.png) 


*Example 2: Ensemble statistics*

You can also download 13 different ensemble statistics about basins or countries aside from the raw data as above.  These derived statistics are often given relative to a control period, 1961-2000.  The time periods however are only for 2046-2065, or 2081-2100, not the standard 20 year intervals available for raw model data.  Below is a full table of available statistics.

|Parameter name|Description|Units|
|---------------|-------|-------|
|*tmin_means*|Average daily minimum temperature|degrees Celsius|
|*tmax_means*|Average daily maximum temperature|degrees Celsius|
|*tmax_days90th*|Number of days with maximum temperature above the control period's 90th percentile (hot days)|days|
|*tmin_days90th*|Number of days with minimum  temperature above the control period's 90th percentile (warm nights)|days|
|*tmax_days10th*|Number of days with maximum temperature below the control period's 10th percentile (cool days)|days|
|*tmin_days10th*|Number of days with minimum  temperature below the control period's 10th percentile (cold nights)|days|
|*tmin_days0*|Number of days with minimum  temperature below 0 degrees Celsius|days|
|*ppt_days*|Number of days with precipitation greater than 0.2 mm|days|
|*ppt_days2*|Number of days with precipitation greater than 2 mm|days|
|*ppt_days10*|Number of days with precipitation greater than 10 mm|days|
|*ppt_days90th*|Number of days with precipitation greater than the control period's 90th percentile|days|
|*ppt_dryspell*|Average number of days between precipitation events|days|
|*ppt_means*|Average daily precipitation|mm|

Similar to our previous example where we looked at temperature anomaly along a latitudinal gradient, we can examine more complex ensemble statistics.  Here we examine the average minimum daily temperature in Scavavadian countries.  Also given that only two time periods are available, there is no need to enter in a specific year range

```r
country.list <- c("ISL", "FIN", "NOR", "SWE")
country.dat <- get_ensemble_stats(country.list, "mavg", "tmin_means")
####### Subset data Exclude A2 scenario
country.dat.b1 <- subset(country.dat, country.dat$scenario == "b1")
# choose just one percentile
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$percentile == 50)
# get just one year period
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$fromYear == 2081)


ggplot(country.dat.b1, aes(x = month, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Average daily minimum temperature") + 
    theme_bw() + xlab("Month")
```

![plot of chunk enesmble data](figure/enesmble_data.png) 



**_Downloading Historical Data_**

It is possible to extract historical data from GCM queries, but this is acutally model backcasting output, not true historical records.  The climate data api can download data at the same spatial scales as all other requests (country or basin), but the time scale differs. You can download data at monthly, yearly or decadanal time scales.  Monthly data is actually the mean monthly temperature or precipitation value for each month from 1901-2009 for countries, or 1960-2009 for basins.  You can indicate the type of data you want in the call by setting the `time_scale` parameter to *month*,*year* or *decade*.

*Example 1: Downloading monthly data*

You can download historical precipitation for Belize, Colombia, Peru and Bolivia.  



```r
country.list <- c("BLZ", "COL", "PER", "BOL")
country.dat <- get_historical_precip(country.list, "month")

ggplot(country.dat, aes(x = month, y = data, group = locator, colour = locator)) + 
    geom_point() + geom_path() + ylab("Average historical precipitation (mm)") + 
    theme_bw() + xlab("Month")
```

![plot of chunk historicalmonth](figure/historicalmonth.png) 


*Example 2: Downloading annual data*

Another use of historical data is to look at increases in temperature over time.

```r
country.list <- c("USA", "MEX", "CAN", "BLZ")
country.dat <- get_historical_temp(country.list, "year")

ggplot(country.dat, aes(x = year, y = data, group = locator)) + geom_point() + 
    geom_path() + ylab("Average annual temperature of Canada") + theme_bw() + 
    xlab("Year") + stat_smooth(se = F, colour = "black") + facet_wrap(~locator, 
    scale = "free")
```

![plot of chunk historicalyear](figure/historicalyear.png) 


**_Mapping climate Data_**

Data can be mapped in [ggplot2](http://docs.ggplot2.org/current/) by downloading KML files from the climate database, creating dataframes and plotting using ggplot2. Maps are available at the basin or country spatial scale.  Because KML files can be large files are downloaded and locally cached in a directory you must set using the kmlpath option.  You can set to any directory as follows: `options(kmlpath = <yourpath>)`.  KML files will be stored there and only downloaded again if they aren't found locally.  If you wish to get rid of them, you will need to manually delete them.  `rWBclimate` allows you to download maps and plot them as ggplot maps, or it can automatically create maps of climate data.  For your convenience we've created basin and country vectors for all six continents *(Antartica has no data in the database)*.  These data files are automatically loaded with the package and can be accessed for basins as *NoAm_basin*, *SoAm_basin*,*Asia_basin*,*Eur_basin*,*Africa_basin*,and *Oceana_basin*. Countries can be similarly accessed as *NoAm_country*, *SoAm_country*,*Asia_country*,*Eur_country*,*Africa_country*,and *Oceana_country*.  If you are plotting by continent, it can take some time to first download all the KML files.

*Example 1: Creating map dataframe and plotting*

Creating map data frames is straightforward, simply provide a list of valid country codes to the function `create_map_df`

```r
# Set the kmlpath option
options(kmlpath = "~/kmltemp")
## Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)
```

```
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=                                                                |   1%
  |                                                                       
  |=                                                                |   2%
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |===                                                              |   4%
  |                                                                       
  |====                                                             |   6%
  |                                                                       
  |====                                                             |   7%
  |                                                                       
  |=====                                                            |   8%
  |                                                                       
  |======                                                           |   9%
  |                                                                       
  |======                                                           |  10%
  |                                                                       
  |=======                                                          |  11%
  |                                                                       
  |========                                                         |  12%
  |                                                                       
  |=========                                                        |  13%
  |                                                                       
  |=========                                                        |  14%
  |                                                                       
  |==========                                                       |  16%
  |                                                                       
  |===========                                                      |  17%
  |                                                                       
  |============                                                     |  18%
  |                                                                       
  |============                                                     |  19%
  |                                                                       
  |=============                                                    |  20%
  |                                                                       
  |==============                                                   |  21%
  |                                                                       
  |==============                                                   |  22%
  |                                                                       
  |===============                                                  |  23%
  |                                                                       
  |================                                                 |  24%
  |                                                                       
  |=================                                                |  26%
  |                                                                       
  |=================                                                |  27%
  |                                                                       
  |==================                                               |  28%
  |                                                                       
  |===================                                              |  29%
  |                                                                       
  |====================                                             |  30%
  |                                                                       
  |====================                                             |  31%
  |                                                                       
  |=====================                                            |  32%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |======================                                           |  34%
  |                                                                       
  |=======================                                          |  36%
  |                                                                       
  |========================                                         |  37%
  |                                                                       
  |=========================                                        |  38%
  |                                                                       
  |=========================                                        |  39%
  |                                                                       
  |==========================                                       |  40%
  |                                                                       
  |===========================                                      |  41%
  |                                                                       
  |===========================                                      |  42%
  |                                                                       
  |============================                                     |  43%
  |                                                                       
  |=============================                                    |  44%
  |                                                                       
  |==============================                                   |  46%
  |                                                                       
  |==============================                                   |  47%
  |                                                                       
  |===============================                                  |  48%
  |                                                                       
  |================================                                 |  49%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================                                |  51%
  |                                                                       
  |==================================                               |  52%
  |                                                                       
  |===================================                              |  53%
  |                                                                       
  |===================================                              |  54%
  |                                                                       
  |====================================                             |  56%
  |                                                                       
  |=====================================                            |  57%
  |                                                                       
  |======================================                           |  58%
  |                                                                       
  |======================================                           |  59%
  |                                                                       
  |=======================================                          |  60%
  |                                                                       
  |========================================                         |  61%
  |                                                                       
  |========================================                         |  62%
  |                                                                       
  |=========================================                        |  63%
  |                                                                       
  |==========================================                       |  64%
  |                                                                       
  |===========================================                      |  66%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |============================================                     |  68%
  |                                                                       
  |=============================================                    |  69%
  |                                                                       
  |==============================================                   |  70%
  |                                                                       
  |==============================================                   |  71%
  |                                                                       
  |===============================================                  |  72%
  |                                                                       
  |================================================                 |  73%
  |                                                                       
  |================================================                 |  74%
  |                                                                       
  |=================================================                |  76%
  |                                                                       
  |==================================================               |  77%
  |                                                                       
  |===================================================              |  78%
  |                                                                       
  |===================================================              |  79%
  |                                                                       
  |====================================================             |  80%
  |                                                                       
  |=====================================================            |  81%
  |                                                                       
  |=====================================================            |  82%
  |                                                                       
  |======================================================           |  83%
  |                                                                       
  |=======================================================          |  84%
  |                                                                       
  |========================================================         |  86%
  |                                                                       
  |========================================================         |  87%
  |                                                                       
  |=========================================================        |  88%
  |                                                                       
  |==========================================================       |  89%
  |                                                                       
  |==========================================================       |  90%
  |                                                                       
  |===========================================================      |  91%
  |                                                                       
  |============================================================     |  92%
  |                                                                       
  |=============================================================    |  93%
  |                                                                       
  |=============================================================    |  94%
  |                                                                       
  |==============================================================   |  96%
  |                                                                       
  |===============================================================  |  97%
  |                                                                       
  |================================================================ |  98%
  |                                                                       
  |================================================================ |  99%
  |                                                                       
  |=================================================================| 100%
```

```r
ggplot(af_basin, aes(x = long, y = lat, group = group)) + geom_polygon() + theme_bw()
```

![plot of chunk map_plot](figure/map_plot.png) 


*Example 2: Mapping climate data*

In order to map climate data you need to have a single point of data for each spatial polygon(country or basin) in your map dataframe.  We've already created the Africa basin map dataframe, so now you need to get some data for each basin.  


```r
af_basin_dat <- get_ensemble_temp(Africa_basin, "annualanom", 2080, 2100)
## Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat, af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat, af_basin_dat$percentile == 50)
```


Now that we have both the map dataframe and the data we can bind them together with the function `climate_map()`.  The function has two options, it can return a `ggplot2` map or it can return a dataframe that you can plot easily with `ggplot2`.  This is useful because it means that you could bind together multiple dataframes for bigger plots, or perhaps add new data yourself


```r

af_map <- climate_map(af_basin, af_basin_dat, return_map = T)
af_map + scale_fill_continuous("Temperature \n anomaly", low = "yellow", high = "red") + 
    theme_bw()
```

![plot of chunk climatemap](figure/climatemap.png) 



*Example 3: Creating a temperature map of the world*


```r
options(kmlpath = "~/kmltemp")
### Combine all country vectors

world <- c(NoAm_country, SoAm_country, Eur_country, Asia_country, Africa_country, 
    Oceana_country)
world_map_df <- create_map_df(world)
```

```
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=                                                                |   1%
  |                                                                       
  |=                                                                |   2%
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |==                                                               |   4%
  |                                                                       
  |===                                                              |   4%
  |                                                                       
  |===                                                              |   5%
  |                                                                       
  |====                                                             |   6%
  |                                                                       
  |====                                                             |   7%
  |                                                                       
  |=====                                                            |   7%
  |                                                                       
  |=====                                                            |   8%
  |                                                                       
  |======                                                           |   8%
  |                                                                       
  |======                                                           |   9%
  |                                                                       
  |======                                                           |  10%
  |                                                                       
  |=======                                                          |  10%
  |                                                                       
  |=======                                                          |  11%
  |                                                                       
  |========                                                         |  12%
  |                                                                       
  |========                                                         |  13%
  |                                                                       
  |=========                                                        |  13%
  |                                                                       
  |=========                                                        |  14%
  |                                                                       
  |==========                                                       |  15%
  |                                                                       
  |==========                                                       |  16%
  |                                                                       
  |===========                                                      |  17%
  |                                                                       
  |============                                                     |  18%
  |                                                                       
  |============                                                     |  19%
  |                                                                       
  |=============                                                    |  19%
  |                                                                       
  |=============                                                    |  20%
  |                                                                       
  |=============                                                    |  21%
  |                                                                       
  |==============                                                   |  21%
  |                                                                       
  |==============                                                   |  22%
  |                                                                       
  |===============                                                  |  22%
  |                                                                       
  |===============                                                  |  23%
  |                                                                       
  |===============                                                  |  24%
  |                                                                       
  |================                                                 |  24%
  |                                                                       
  |================                                                 |  25%
  |                                                                       
  |=================                                                |  25%
  |                                                                       
  |=================                                                |  26%
  |                                                                       
  |=================                                                |  27%
  |                                                                       
  |==================                                               |  27%
  |                                                                       
  |==================                                               |  28%
  |                                                                       
  |===================                                              |  29%
  |                                                                       
  |===================                                              |  30%
  |                                                                       
  |====================                                             |  30%
  |                                                                       
  |====================                                             |  31%
  |                                                                       
  |=====================                                            |  32%
  |                                                                       
  |=====================                                            |  33%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |======================                                           |  34%
  |                                                                       
  |=======================                                          |  35%
  |                                                                       
  |=======================                                          |  36%
  |                                                                       
  |========================                                         |  36%
  |                                                                       
  |========================                                         |  37%
  |                                                                       
  |=========================                                        |  38%
  |                                                                       
  |=========================                                        |  39%
  |                                                                       
  |==========================                                       |  39%
  |                                                                       
  |==========================                                       |  40%
  |                                                                       
  |==========================                                       |  41%
  |                                                                       
  |===========================                                      |  41%
  |                                                                       
  |===========================                                      |  42%
  |                                                                       
  |============================                                     |  42%
  |                                                                       
  |============================                                     |  43%
  |                                                                       
  |============================                                     |  44%
  |                                                                       
  |=============================                                    |  44%
  |                                                                       
  |=============================                                    |  45%
  |                                                                       
  |==============================                                   |  46%
  |                                                                       
  |==============================                                   |  47%
  |                                                                       
  |===============================                                  |  47%
  |                                                                       
  |===============================                                  |  48%
  |                                                                       
  |================================                                 |  49%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================                                |  50%
  |                                                                       
  |=================================                                |  51%
  |                                                                       
  |==================================                               |  52%
  |                                                                       
  |==================================                               |  53%
  |                                                                       
  |===================================                              |  53%
  |                                                                       
  |===================================                              |  54%
  |                                                                       
  |====================================                             |  55%
  |                                                                       
  |====================================                             |  56%
  |                                                                       
  |=====================================                            |  56%
  |                                                                       
  |=====================================                            |  57%
  |                                                                       
  |=====================================                            |  58%
  |                                                                       
  |======================================                           |  58%
  |                                                                       
  |======================================                           |  59%
  |                                                                       
  |=======================================                          |  59%
  |                                                                       
  |=======================================                          |  60%
  |                                                                       
  |=======================================                          |  61%
  |                                                                       
  |========================================                         |  61%
  |                                                                       
  |========================================                         |  62%
  |                                                                       
  |=========================================                        |  63%
  |                                                                       
  |=========================================                        |  64%
  |                                                                       
  |==========================================                       |  64%
  |                                                                       
  |==========================================                       |  65%
  |                                                                       
  |===========================================                      |  66%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |============================================                     |  67%
  |                                                                       
  |============================================                     |  68%
  |                                                                       
  |=============================================                    |  69%
  |                                                                       
  |=============================================                    |  70%
  |                                                                       
  |==============================================                   |  70%
  |                                                                       
  |==============================================                   |  71%
  |                                                                       
  |===============================================                  |  72%
  |                                                                       
  |===============================================                  |  73%
  |                                                                       
  |================================================                 |  73%
  |                                                                       
  |================================================                 |  74%
  |                                                                       
  |================================================                 |  75%
  |                                                                       
  |=================================================                |  75%
  |                                                                       
  |=================================================                |  76%
  |                                                                       
  |==================================================               |  76%
  |                                                                       
  |==================================================               |  77%
  |                                                                       
  |==================================================               |  78%
  |                                                                       
  |===================================================              |  78%
  |                                                                       
  |===================================================              |  79%
  |                                                                       
  |====================================================             |  79%
  |                                                                       
  |====================================================             |  80%
  |                                                                       
  |====================================================             |  81%
  |                                                                       
  |=====================================================            |  81%
  |                                                                       
  |=====================================================            |  82%
  |                                                                       
  |======================================================           |  83%
  |                                                                       
  |=======================================================          |  84%
  |                                                                       
  |=======================================================          |  85%
  |                                                                       
  |========================================================         |  86%
  |                                                                       
  |========================================================         |  87%
  |                                                                       
  |=========================================================        |  87%
  |                                                                       
  |=========================================================        |  88%
  |                                                                       
  |==========================================================       |  89%
  |                                                                       
  |==========================================================       |  90%
  |                                                                       
  |===========================================================      |  90%
  |                                                                       
  |===========================================================      |  91%
  |                                                                       
  |===========================================================      |  92%
  |                                                                       
  |============================================================     |  92%
  |                                                                       
  |============================================================     |  93%
  |                                                                       
  |=============================================================    |  93%
  |                                                                       
  |=============================================================    |  94%
  |                                                                       
  |==============================================================   |  95%
  |                                                                       
  |==============================================================   |  96%
  |                                                                       
  |===============================================================  |  96%
  |                                                                       
  |===============================================================  |  97%
  |                                                                       
  |================================================================ |  98%
  |                                                                       
  |================================================================ |  99%
  |                                                                       
  |=================================================================| 100%
```

```r
world_dat <- get_ensemble_temp(world, "annualavg", 2080, 2100)
## Subset data to just one scenario, and one percentile
world_dat <- subset(world_dat, world_dat$scenario == "a2")
world_dat <- subset(world_dat, world_dat$percentile == 50)

world_map <- climate_map(world_map_df, world_dat, return_map = T)
world_map + scale_fill_continuous("Temperature \n at 2100", low = "yellow", 
    high = "red") + theme_bw()
```

![plot of chunk worldmap](figure/worldmap.png) 

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Oceana_basin}
\alias{Oceana_basin}
\title{Basin codes for Oceana, used in downloading maps}
\description{
Basin codes for Oceana, used in downloading maps
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ensemble_climate_data.R
\name{get_ensemble_climate_data}
\alias{get_ensemble_climate_data}
\title{Download ensemble climate data}
\usage{
get_ensemble_climate_data(locator, geo_type, type, cvar, start, end)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about.
(http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]}

\item{geo_type}{basin or country depending on the locator type}

\item{type}{the type of data you want "mavg" for monthly averages, "annualavg"}

\item{cvar}{The variable you're interested in. "pr" for precipitation, "tas" for
temperature in celcius.}

\item{start}{The starting year you want data for, can be in the past or the future.
Must conform to the periods outlined in the world bank API.  If given values don't
conform to dates, the fuction will automatically round them.}

\item{end}{The ending year you want data for, can be in the past or the future.
Similar to the start date, dates will be rounded to the nearest end dat.}
}
\description{
Download ensemble data for all models, returns the 10th, 50th and 90th
percentile of all models (15 for A1, 13 for B2).  Ensemble requets can be for
countries or basins.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_model_temp.R
\name{get_model_temp}
\alias{get_model_temp}
\title{Download GCM temperature data}
\usage{
get_model_temp(locator, type, start, end)
}
\arguments{
\item{locator}{A vector of either watershed basin ID's from http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf
It can be just a single basin id, or a vector of ids.  ids should be strings.}

\item{type}{the type of data to retrieve, must be "mavg" for monthly averages,
"annualavg" for annual averages, "manom" for monthly anomaly, and "annualanom" 
 for annual anomaly.}

\item{start}{the start year to gather data for.}

\item{end}{the end year to gather data to.}
}
\value{
a dataframe with temperature predictions in degrees C for all scenarios, gcms, for each time period.
}
\description{
Function wraps get_climate_data() and returns temperature
by basin or country in degrees C as output from all 15 models, for the a1 and b2 scenarios.
}
\details{
start and end year can be any years, but all years will be coerced
        into periods outlined by the API (http://data.worldbank.org/developers/climate  -data-api)
        anomaly periods are only valid for future scenarios and based on a 
        reference period of 1969 - 1999, see API for full details.
}
\examples{
\dontrun{
# Get data for 2 basins, annual average temperature for all valid time periods
# then subset them, and plot
temp_dat <- get_model_temp(c("2","231"),"annualavg",1900,3000)
temp_dat <- subset(temp_dat,temp_dat$gcm=="ukmo_hadcm3")
temp_dat <- subset(temp_dat,temp_dat$scenario!="b1")
ggplot(temp_dat,aes(x=fromYear,y=data,group=locator,
colour=locator))+geom_path()

### Get data for 4 countries with monthly tempitation values
temp_dat <- get_model_temp(c("USA","BRA","CAN","YEM"),"mavg",2020,2030)
temp_dat <- subset(temp_dat,temp_dat$gcm=="ukmo_hadcm3")
temp_dat <- subset(temp_dat,temp_dat$scenario!="b1")
ggplot(temp_dat,aes(x=as.factor(month),y=data,group=locator,
colour=locator))+geom_path()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_kml.R
\name{download_kml}
\alias{download_kml}
\title{Download kml files}
\usage{
download_kml(locator)
}
\arguments{
\item{locator}{The a vector of ISO3 country code's that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID's [1-468] (http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)}
}
\description{
Downloads map data from in kml format and writes it to a temporary directory.  You must specify a temporary directory to write files to in your options.
}
\details{
kml files can be quite large making downloading them every time you want to make a map time consuming.  To 
reduce this time it's easiest to download kml files and store them.  To set the directory use a line like this: \code{options(kmlpath="/Users/emh/kmltemp")}  The option must be called "kmlpath".
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_climate_data.R
\name{get_climate_data}
\alias{get_climate_data}
\title{get_climate_data}
\usage{
get_climate_data(locator, geo_type, type, cvar, start, end)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]}

\item{geo_type}{basin or country depending on the locator type}

\item{type}{the type of data you want "mavg" for monthly averages, "annualavg"}

\item{cvar}{The variable you're interested in. "pr" for precipitation, "tas" for temperature in celcius.}

\item{start}{The starting year you want data for, can be in the past or the future. Must conform to the periods outlined in the world bank API.  If given values don't conform to dates, the fuction will automatically round them.}

\item{end}{The ending year you want data for, can be in the past or the future.  Similar to the start date, dates will be rounded to the nearest end dat.}
}
\description{
Download monthly average climate data from the world bank climate 
            data api. Ideally you'll want to use the wrapper functions that call this.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_model_precip.R
\name{get_model_precip}
\alias{get_model_precip}
\title{Download GCM precipitation data}
\usage{
get_model_precip(locator, type, start, end)
}
\arguments{
\item{locator}{A vector of either watershed basin ID's from http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf
It can be just a single basin id, or a vector of ids.  ids should be strings.}

\item{type}{the type of data to retrieve, must be "mavg" for monthly averages,
"annualavg" for annual averages, "manom" for monthly anomaly, and "annualanom" for 
annual anomaly.}

\item{start}{the start year to gather data for.}

\item{end}{the end year to gather data to.}
}
\value{
a dataframe with precipitation predictions in mm for all scenarios, gcms, for each time period.
}
\description{
Function wraps get_climate_data() and returns precipitation
by basin or country in mm as output from all 15 models, for the a1 and b2 scenarios.
}
\details{
start and end year can be any years, but all years will be coerced
        into periods outlined by the API (http://data.worldbank.org/developers/climate  -data-api)
        anomaly periods are only valid for future scenarios and based on a 
        reference period of 1969 - 1999, see API for full details.
}
\examples{
\dontrun{
# Get data for 2 basins, annual average precipitation for all valid time periods
# then subset them, and plot
precip_dat <- get_model_precip(c("2","231"),"annualavg",1900,3000)
precip_dat <- subset(precip_dat,precip_dat$gcm=="ukmo_hadcm3")
precip_dat <- subset(precip_dat,precip_dat$scenario!="b1")
ggplot(precip_dat,aes(x=fromYear,y=annualData,group=locator,colour=locator))+geom_path()

### Get data for 4 countries with monthly precipitation values
precip_dat <- get_model_precip(c("USA","BRA","CAN","YEM"),"mavg",2020,2030)
precip_dat <- subset(precip_dat,precip_dat$gcm=="ukmo_hadcm3")
precip_dat <- subset(precip_dat,precip_dat$scenario!="b1")
ggplot(precip_dat,aes(x=as.factor(month),y=monthVals,group=locator,colour=locator))+geom_path()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_historical_precip.R
\name{get_historical_precip}
\alias{get_historical_precip}
\title{Download historical precipitation data}
\usage{
get_historical_precip(locator, time_scale)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]. This can be a vector of all basins or all countries.}

\item{time_scale}{The time scale you want to return values on.  Must be "\emph{month}", "\emph{year}" or "\emph{decade}"}
}
\value{
a dataframe with historical precipitation data
}
\description{
The Climate Data API provides access to historical precipitation data. These data are separate from the outputs of the GCMs,
and they are based on gridded climatologies from the Climate Research Unit.
}
\details{
The historical period for country is 1901 - 2009, and 1960 - 2009 for basin. The time_scale parameter returns a different number of variables depending on the input timescale. \emph{Month} will return 12 values, a historical average for that month across all years.  \emph{Year} will return yearly averages for each year, and \emph{decade} will return decade averages.
}
\examples{
\dontrun{
## Plot annual historical data for USA, Brazil and Australia
hist_dat <- get_historical_precip(c("USA","BRA","AUS"),"year")
ggplot(hist_dat,aes(x = year,y = data, group = locator,
colour = locator)) + geom_point() + geom_path() + ylab("Mean annual precipitaion")

## Plot monthly historical data
hist_mo_dat <- get_historical_precip(c("USA","AUS","BRA","IDN"),time_scale="month")
ggplot(hist_mo_dat,aes(x = month,y = data, group = locator,
colour = locator)) + geom_point() + geom_path() + ylab("Mean monthly precipitaion")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{codes}
\alias{codes}
\title{isocodes data}
\description{
isocodes data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ensemble_temp.R
\name{get_ensemble_temp}
\alias{get_ensemble_temp}
\title{Download ensemble temperature data}
\usage{
get_ensemble_temp(locator, type, start, end)
}
\arguments{
\item{locator}{A vector of either watershed basin ID's from
http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf
It can be just a single basin id, or a vector of ids.  ids should be strings.}

\item{type}{the type of data to retrieve, must be "mavg" for monthly averages,
"annualavg" for annual averages, "manom" for monthly anomaly, and "annualanom" for
annual anomaly.}

\item{start}{the start year to gather data for.}

\item{end}{the end year to gather data to.}
}
\value{
a dataframe with precipitation predictions in mm for all scenarios, gcms,
for each time period.
}
\description{
Function wraps \code{\link{get_ensemble_climate_data}} and returns
precipitation by basin or country in mm.  Output is the 10th 50th and 90th percentile
for all gcm's for the a1 and b2 scenarios.
}
\details{
start and end year can be any years, but all years will be coerced
into periods outlined by the API (http://data.worldbank.org/developers/climate-data-api)
anomaly periods are only valid for future scenarios and based on a
reference period of 1969 - 1999, see API for full details.
}
\examples{
\dontrun{
# Get data for 2 basins, annual average precipitation for all valid time periods
# then subset them, and plot
temp_dat <- get_ensemble_temp(locator=c(2,231), type="annualavg", start=1900, end=3000)
temp_dat <- subset(temp_dat,temp_dat$scenario!="b1")
temp_dat$uniqueGroup <- paste(temp_dat$percentile,temp_dat$locator,sep="-")
ggplot(temp_dat, aes(x=fromYear, y=data, group=uniqueGroup,
       colour=as.factor(locator), linetype=as.factor(percentile))) +
   geom_path()

### Get data for 2 countries with monthly precipitation values
temp_dat <- get_ensemble_temp(locator = c("USA","BRA"), type = "mavg", start = 2020, end = 2030)
temp_dat <- subset(temp_dat,temp_dat$scenario!="b1")
temp_dat$uniqueGroup <- paste(temp_dat$percentile, temp_dat$locator,sep="-")
ggplot(temp_dat, aes(x=as.factor(month), y=data, group=uniqueGroup, colour=locator)) +
   geom_path()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{package}
\name{rWBclimate}
\alias{rWBclimate}
\title{rWBclimate}
\description{
rWBclimate
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_map_df.R
\name{create_map_df}
\alias{create_map_df}
\title{Create mapable dataframe}
\usage{
create_map_df(locator)
}
\arguments{
\item{locator}{The a vector of ISO3 country code's that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID's [1-468] (http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf)}
}
\description{
A function that will download maps for a vector of basins or country codes and return a data frame that has the kml output processed such that it can be plotted with ggplot2 and other mapping functions:
}
\details{
kml files can be quite large (100k-600k per country) making downloading them every time you want to make a map time consuming.  To 
reduce this time it's easiest to download kml files and store them.  To set the directory use a line like this: \code{options(kmlpath="/Users/emh/kmltemp")}  The option must be called "kmlpath".  These files will be persistent until you delete them.
}
\examples{
\dontrun{
to_map <- create_map_df(c("USA","MEX","CAN"))
ggplot(to_map, aes(x=long, y=lat,group=group))+ geom_polygon()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_historical_data.R
\name{get_historical_data}
\alias{get_historical_data}
\title{Download historical climate data}
\usage{
get_historical_data(locator, cvar, time_scale)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468].
The historical period for country is 1901 - 2009, and 1960 - 2009 for basin}

\item{cvar}{The climate variable you're interested in. "\emph{pr}" for precipitation, "\emph{tas}" for temperature in celcius.}

\item{time_scale}{The time scale you want to return values on.  Must be "\emph{month}", "\emph{year}" or "\emph{decade}"}
}
\value{
a dataframe with historical climate data
}
\description{
The Climate Data API provides access to historical temperature
and precipitation data. These data are separate from the outputs of the GCMs,
and they are based on gridded climatologies from the Climate Research Unit.
}
\details{
The time_scale parameter returns a different number of variables depending on the input timescale. \emph{Month} will return 12 values, a historical average for that month across all years.  \emph{Year} will return yearly averages for each year, and \emph{decade} will return decade averages.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_historical_data_recursive.R
\name{get_historical_data_recursive}
\alias{get_historical_data_recursive}
\title{Download historical climate data recursively}
\usage{
get_historical_data_recursive(locator, cvar, time_scale)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468].
The historical period for country is 1901 - 2009, and 1960 - 2009 for basin}

\item{cvar}{The climate variable you're interested in. "\emph{pr}" for precipitation, "\emph{tas}" for temperature in celcius.}

\item{time_scale}{The time scale you want to return values on.  Must be "\emph{month}", "\emph{year}" or "\emph{decade}"}
}
\value{
a dataframe with historical climate data
}
\description{
Recursively get historical data
}
\details{
The time_scale parameter returns a different number of variables depending on the input timescale. \emph{Month} will return 12 values, a historical average for that month across all years.  \emph{Year} will return yearly averages for each year, and \emph{decade} will return decade averages.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ensemble_precip.R
\name{get_ensemble_precip}
\alias{get_ensemble_precip}
\title{Download ensemble precipitation data}
\usage{
get_ensemble_precip(locator, type, start, end)
}
\arguments{
\item{locator}{A vector of either watershed basin ID's from http://data.worldbank.org/sites/default/files/climate_data_api_basins.pdf
It can be just a single basin id, or a vector of ids.  ids should be strings.}

\item{type}{the type of data to retrieve, must be "mavg" for monthly averages,
"annualavg" for annual averages, "manom" for monthly anomaly, and "annualanom" for 
annual anomaly.}

\item{start}{the start year to gather data for.}

\item{end}{the end year to gather data to.}
}
\value{
a dataframe with precipitation predictions in mm for all scenarios, gcms, for each time period.
}
\description{
Function wraps get_ensemble_climate_data() and returns precipitation
by basin or country in mm.  Output is the 10th 50th and 90th percentile for all 
gcm's for the a1 and b2 scenarios.
}
\details{
start and end year can be any years, but all years will be coerced
        into periods outlined by the API (http://data.worldbank.org/developers/climate  -data-api)
        anomaly periods are only valid for future scenarios and based on a 
        reference period of 1969 - 1999, see API for full details.
}
\examples{
\dontrun{
# Get data for 2 basins, annual average precipitation for all valid time periods
# then subset them, and plot
precip_dat <- get_ensemble_precip(c("2","231"),"annualavg",1900,3000)
precip_dat <- subset(precip_dat,precip_dat$scenario!="b1")
precip_dat$uniqueGroup <- paste(precip_dat$percentile,precip_dat$locator,sep="-")
ggplot(precip_dat,aes(x=fromYear,y=annualVal,group=uniqueGroup,colour=as.factor(locator),
linetype=as.factor(percentile)))+ geom_path()

### Get data for 2 countries with monthly precipitation values
precip_dat <- get_ensemble_precip(c("USA","BRA"),"mavg",2020,2030)
precip_dat <- subset(precip_dat,precip_dat$scenario!="b1")
precip_dat$uniqueGroup <- paste(precip_dat$percentile,precip_dat$locator,sep="-")
ggplot(precip_dat,aes(x=as.factor(month),y=monthVals,group=uniqueGroup,
colour=locator))+geom_path()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{NoAm_basin}
\alias{NoAm_basin}
\title{Basin codes for NoAm, used in downloading maps}
\description{
Basin codes for NoAm, used in downloading maps
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Asia_country}
\alias{Asia_country}
\title{Country codes for all of Asia}
\description{
Country codes for all of Asia
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{SoAm_basin}
\alias{SoAm_basin}
\title{Basin codes for SoAm, used in downloading maps}
\description{
Basin codes for SoAm, used in downloading maps
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_ISO_code.R
\name{check_ISO_code}
\alias{check_ISO_code}
\title{check country codes}
\usage{
check_ISO_code(iso)
}
\arguments{
\item{iso}{The 3 letter country code based on ISO3 Country abbreviations (http://unstats.un.org/unsd/methods/m49/m49alpha.htm)}
}
\value{
TRUE if a valid code, otherwise an error is returned
}
\description{
Checks if the country code entered is a valid country code that data exists for
}
\examples{
\dontrun{
check_ISO_code("USA")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate_map.R
\name{climate_map}
\alias{climate_map}
\title{Map climate data}
\usage{
climate_map(map_df, data_df, return_map = TRUE)
}
\arguments{
\item{map_df}{a map dataframe generated from create_map_df()}

\item{data_df}{a climate dataframe with one piece of data to be mapped to each unique spatial polygon.}

\item{return_map}{True returns a ggplot2 object, False returns a dataframe where data items are matched to their polygon that you can plot later on.}
}
\value{
Either a ggplot2 map or a dataframe depending on the parameter return_map
}
\description{
Create maps of climate data.  It requires two data inputs, a map dataframe, and a climate dataframe.  The climate data must have one data point per spatial mapping point,e.g. 1 data point per country or basin being mapped.
}
\examples{
\dontrun{
#Set the kmlpath option
options(kmlpath = "~/kmltemp2")
##Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)
af_basin_dat <- get_ensemble_temp(Africa_basin,"annualanom",2080,2100)
##  Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat,af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat,af_basin_dat$percentile == 50)
af_map <- climate_map(map_df = af_basin, data_df = af_basin_dat, return_map = TRUE)
af_map + scale_fill_continuous("Temperature \n anomaly",low="yellow",high = "red") + theme_bw()

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Eur_country}
\alias{Eur_country}
\title{Country codes for all of Eur}
\description{
Country codes for all of Eur
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_locator.R
\name{check_locator}
\alias{check_locator}
\title{Checks for what kind of locator a user input}
\usage{
check_locator(locator)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]}
}
\value{
geo_ref a string indicating what kind of geography to use in the api
}
\description{
Checks for what kind of locator a user input
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Africa_country}
\alias{Africa_country}
\title{Country codes for all of Africa}
\description{
Country codes for all of Africa
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data_recursive.R
\name{get_data_recursive}
\alias{get_data_recursive}
\title{wratpper for get_climate_data()}
\usage{
get_data_recursive(locator, geo_type, type, cvar, start, end)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]}

\item{geo_type}{basin or country depending on the locator type}

\item{type}{the type of data you want "mavg" for monthly averages, "annualavg"}

\item{cvar}{The variable you're interested in. "pr" for precipitation, "tas" for temperature in celcius.}

\item{start}{The starting year you want data for, can be in the past or the future. Must conform to the periods outlined in the world bank API.  If given values don't conform to dates, the fuction will automatically round them.}

\item{end}{The ending year you want data for, can be in the past or the future.  Similar to the start date, dates will be rounded to the nearest end dat.}
}
\description{
Function to recursively call the get_climate_data().  Handles a vector of basins
or countries as well as multiple dates.
}
\examples{
\dontrun{
 get_ensemble_data_recursive(c("1","2"),"basin","mavg","pr",1920,1940)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ensemble_data_recursive.R
\name{get_ensemble_data_recursive}
\alias{get_ensemble_data_recursive}
\title{Wrapper for get_ensemble_climate_data()}
\usage{
get_ensemble_data_recursive(locator, geo_type, type, cvar, start, end)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]}

\item{geo_type}{basin or country depending on the locator type}

\item{type}{the type of data you want "mavg" for monthly averages, "annualavg"}

\item{cvar}{The variable you're interested in. "pr" for precipitation, "tas" for temperature in celcius.}

\item{start}{The starting year you want data for, can be in the past or the future. Must conform to the periods outlined in the world bank API.  If given values don't conform to dates, the fuction will automatically round them.}

\item{end}{The ending year you want data for, can be in the past or the future.  Similar to the start date, dates will be rounded to the nearest end dat.}
}
\description{
Function to recursively call the get_ensemble_climate_data().  Handles a vector of basins
or countries as well as multiple dates.
}
\examples{
\dontrun{
 get_ensemble_data_recursive(c("1","2"),"basin","mavg","pr",1920,1940)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ensemble_stats.R
\name{get_ensemble_stats}
\alias{get_ensemble_stats}
\title{Download ensemble statistics}
\usage{
get_ensemble_stats(locator, type, stat)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]}

\item{type}{the type of data you want "mavg" for monthly averages, "annualavg"}

\item{stat}{The statistics of interest, must be one of the ones listed above.}
}
\description{
Statistics can be from either two time periods: 2046 - 2065 and 2081 - 2100 and are all given in units relative to a control period:  1961 - 2000.  Derived statistics can be any of the following:

\tabular{lll}{
\strong{Statistic} \tab \strong{Description} \tab \strong{Units} \cr
\emph{tmin_means} \tab Average daily minimum  temperature \tab degrees Celsius \cr 
\emph{tmax_means} \tab Average daily maximum temperature \tab degrees Celsius \cr 
\emph{tmax_days90th} \tab Number of days with maximum temperature above the control period 90th percentile (hot days) \tab days \cr
\emph{tmin_days90th} \tab Number of days with minimum  temperature above the control period 90th percentile (warm nights) \tab days \cr
\emph{tmax_days10th} \tab Number of days with maximum temperature below the control period 10th percentile (cool days) \tab days \cr
\emph{tmin_days10th} \tab Number of days with minimum  temperature below the control period 10th percentile (cold nights) \tab days \cr
\emph{tmin_days0} \tab Number of days with minimum  temperature below 0 degrees Celsius \tab days \cr 
\emph{ppt_days} \tab Number of days with precipitation greater than 0.2 mm \tab days \cr 
\emph{ppt_days2} \tab Number of days with precipitation greater than 2 mm \tab days \cr 
\emph{ppt_days10} \tab Number of days with precipitation greater than 10 mm \tab days \cr 
\emph{ppt_days90th} \tab Number of days with precipitation greater than the control periods 90th percentile \tab days \cr 
\emph{ppt_dryspell} \tab Average number of days between precipitation events \tab days \cr 
\emph{ppt_means} \tab Average daily precipitation \tab mm 
}
}
\examples{
\dontrun{
 ### Request data on the US for days of rain over 2 mm
 ens_dat <- get_ensemble_stats("USA","mavg","ppt_days2")
 # subset to the 50th percentile and just until the year 2100
 ens_dat <- subset(ens_dat, ens_dat$percentile == 50)
 ens_dat <- subset(ens_dat,ens_dat$toYear == 2100)
 ggplot(ens_dat,aes(x = as.factor(month), y= monthVals, group=scenario, 
 colour=scenario)) + geom_point() + geom_line()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Oceana_country}
\alias{Oceana_country}
\title{Country codes for all of Oceana}
\description{
Country codes for all of Oceana
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kml_to_sp.R
\name{kml_to_sp}
\alias{kml_to_sp}
\title{Convert kml to polygon}
\usage{
kml_to_sp(map_df, df = NULL, crs_string = "+proj=longlat +datum=WGS84")
}
\arguments{
\item{map_df}{a map dataframe generated from create_map_df()}

\item{df}{a climate dataframe with one piece of data to be mapped to each unique spatial polygon.}

\item{crs_string}{Coordinate reference string to use in spatial projection.  Default is WSGS84.}
}
\value{
a SpatialPolygon object
}
\description{
Create an sp SpatialPolygon or SpatialPolygonDataFrame object from a downloaded KML file and data file
}
\details{
If a dataframe is included, a spatial polygon dataframe object is created.  The dataframe must have one unique piece of information per polygon, otherwise an error will be thrown.  However just a basic spatial polygon will be created if no dataframe is included.
}
\examples{
\dontrun{
sa_map <- create_map_df(locator=SoAm_country)
sa_dat <- get_ensemble_temp(SoAm_country,"annualanom",2080,2100)
sa_dat <- subset(sa_dat,sa_dat$scenario == "a2")
sa_dat <- subset(sa_dat,sa_dat$percentile == 50)
sa_poly <- kml_to_sp(sa_map,df = sa_dat)
### colors are a bit off, but just to verify that data is 
spplot(sa_poly,"data")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Eur_basin}
\alias{Eur_basin}
\title{Basin codes for Eur, used in downloading maps}
\description{
Basin codes for Eur, used in downloading maps
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{NoAm_country}
\alias{NoAm_country}
\title{Country codes for all of NoAm}
\description{
Country codes for all of NoAm
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Asia_basin}
\alias{Asia_basin}
\title{Basin codes for Asia, used in downloading maps}
\description{
Basin codes for Asia, used in downloading maps
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{SoAm_country}
\alias{SoAm_country}
\title{Country codes for all of SoAm}
\description{
Country codes for all of SoAm
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rWBclimate-package.r
\docType{data}
\name{Africa_basin}
\alias{Africa_basin}
\title{Basin codes for Africa, used in downloading maps}
\description{
Basin codes for Africa, used in downloading maps
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_historical_temp.R
\name{get_historical_temp}
\alias{get_historical_temp}
\title{Download historical temperature data}
\usage{
get_historical_temp(locator, time_scale)
}
\arguments{
\item{locator}{The ISO3 country code that you want data about. (http://unstats.un.org/unsd/methods/m49/m49alpha.htm) or the basin ID [1-468]. This can be a vector of all basins or all countries.}

\item{time_scale}{The time scale you want to return values on.  Must be "\emph{month}", "\emph{year}" or "\emph{decade}"}
}
\value{
a dataframe with historical temperature data
}
\description{
The Climate Data API provides access to historical precipitation data. These data are separate from the outputs of the GCMs,
and they are based on gridded climatologies from the Climate Research Unit.
}
\details{
The historical period for country is 1901 - 2009, and 1960 - 2009 for basin. The time_scale parameter returns a different number of variables depending on the input timescale. \emph{Month} will return 12 values, a historical average for that month across all years.  \emph{Year} will return yearly averages for each year, and \emph{decade} will return decade averages.
}
\examples{
\dontrun{
## Plot annual historical data for USA, Brazil and Australia
hist_dat <- get_historical_precip(c("USA","BRA","AUS"),"year")
ggplot(hist_dat,aes(x = year,y = data, group = locator,
colour = locator)) + geom_point() + geom_path() + ylab("Mean annual temperature")

## Plot monthly historical data
hist_mo_dat <- get_historical_precip(c("USA","AUS","BRA","IDN"),time_scale="month")
ggplot(hist_mo_dat,aes(x = month,y = data, group = locator,
colour = locator)) + geom_point() + geom_path() + ylab("Mean monthly temperature")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/date_correct.R
\name{date_correct}
\alias{date_correct}
\title{correct data values}
\usage{
date_correct(start, end)
}
\arguments{
\item{start}{The start year}

\item{end}{The end year}
}
\value{
a 2xM matrix where M in the number of periods in the data api
}
\description{
Round start and end dates to conform with data api standards.  See api documentation (http://data.worldbank.org/developers/climate-data-api)
for full details
}
\examples{
\dontrun{
date_correct(1921,1957)
}
}
