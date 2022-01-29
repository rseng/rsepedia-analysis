---
title: 'pollucheck v1.0: A package to explore open-source air pollution data'
authors:
- affiliation: '1'
  name: Adithi R. Upadhya
  orcid: 0000-0002-1764-1379
- affiliation: '1'
  name: Pratyush Agrawal
- affiliation: '2'
  name: Sreekanth Vakacherla
  orcid: 0000-0003-0400-6584
- affiliation: '1'
  name: Meenakshi Kushwaha
date: "12 June 2021"
bibliography: paper.bib
tags:
- R
- open-source
- air quality
- shiny
affiliations:
- index: '1'
  name: ILK Labs, Bengaluru, India
- index: '2'
  name: Center for Study of Science, Technology and Policy, Bengaluru, India
---

# Summary

Air pollution impacts human health, quality of living, climate, and the economy [@Hystad:2020]. To assess its impact and facilitate mitigation actions, quantification of air pollution is vital. Measurements are the most accurate way of quantifying air pollution. Many countries conduct regulatory measurements of various air pollutants (e.g., fine and respirable particulate matter, nitrogen dioxide, sulfur dioxide, and surface ozone) and make the data available publicly.

Air pollution data sets typically span several seasons or years and real-time data are recorded typically every hour or at a higher frequency. With the ever increasing amount of data and number of data providers, there is a clear need for tools to handle, analyze, and visualize large data sets. The current Shiny app `pollucheck` aims at a simple workflow to generate a suite of statistical plots and summary statistics [@Shiny:2021].  Users do  not need any programming background to analyze time series data and generate a variety of plots.

`pollucheck` can handle real-time pollution and co-located meteorological data (if available) from the three most popular open-source air pollution databases: [OpenAQ](openaq.org), [AirNow](airnow.gov), and [Indian Central Pollution Control Board (CPCB) dashboard](app.cpcbccr.com). While CPCB data are specific to Indian regulatory monitoring stations, OpenAQ hosts the global open-source pollution databases and AirNow hosts the global PM~2.5~ (mass concentration of particulate matter with an aerodynamic diameter less than or equal to 2.5 microns) data, collected under the United States Embassy and Consulates' air quality monitoring programmes.

The output of `pollucheck` is displayed in seven tabs. Different packages used for building `pollucheck` include `tidyverse`, `openair`, `shiny`, `bslib`, `forecast`, `biwavelet`, `readxl`, `DT`, `data.table`, `nortest`, and `zoo`.

# Statement of Need

Pollution data from these sources are typically in different file formats and templates that require customised codes or programmes for analysis. Also, a rigorous quality check of the data is preferred before visualization (plotting) and reporting. `pollucheck` offers a single-stop solution for

(i) handling the pollution data from the open-source databases,
(ii) applying a suite of quality check options,
(iii) generating a variety of summary statistics at various averaging intervals,
(iv) performing time series analysis,
(v) generating a bunch of temporal and statistical plots, and 
(vi) comparing data from two input files.

To our knowledge, currently there is no application that can generate utilisable summary statistics and plots using the data from the pollution databases. However, there are a few Shiny apps that deal with data cleaning and visualization of pollution data collected from single/multiple air quality instruments [@Salmon:2017; @Upadhya:2020].

# App Display

i) The `File` tab is used to upload the input file and to specify the source and time resolution of the input data. The default time zone is set to *Asia/Kolkata*. For OpenAQ and AirNow data sets, appropriate time zones need to be selected based on the input file. For the CPCB data set, the time zone option is default and inactive. A set of quality check options for the (a) removal of negative values, (b) removal of consecutive duplicate values, and (c) detection of outliers are provided. Data completeness criteria (minimum percentage of data required) for computing daily mean values can be specified. If the input file contains simultaneous PM~2.5~ and PM~10~ (mass concentration of particulate matter with an aerodynamic diameter less than or equal to 10 microns) data, the app computes the PM~2.5~/PM~10~ ratio, a useful metric in the air pollution field to identify sources of PM and to estimate PM~2.5~ when only PM~10~ is available [@Chan:2008; @Chu:2015; @Spandana:2021]. The selected quality check or completeness criteria will be applied to all the parameters of the input file. Hourly or daily mean values of all the parameters can be displayed and downloaded (as `.csv`) from this tab.

ii) The `Summary` tab provides various statistics (central tendencies, percentiles, minimum, maximum, standard deviation, interquartile range, etc.) for all the parameters in the input file at three different averaging intervals. The averaging intervals can be selected using the drop-down menu. The displayed statistics can be downloaded.

iii) The `Summary Plots` tab generates (a) a data availability plot for all the parameters (based on daily mean values), (b) a time series plot, (c) box and whisker plots, (d) a vertical bar plot, and (e) diurnal variability plots. Except for the data availability plot, the parameter of interest to plot needs to be selected from the drop-down menu. Plots can be generated using hourly or daily mean data. The diurnal variability plots can be plotted either by aggregating the whole data in the input file or month wise. Considering the general log-normal nature of the pollution data, an option is provided for the diurnal variability plots to be plotted using mean and standard deviations or median and interquartile ranges. The title and y-axis labels of the plots are editable.

iv) The `Statistical Plots` tab can be used to conduct normality tests (Anderson-Darling and Shapiro-Wilk), generate density and quantile-quantile (QQ) plots, generate autocorrelogram, and conduct trends and periodicity analysis on the parameter selected. While autocorrelogram is generated based on monthly mean values, trend (the Mann-Kendall test) and periodicity (wavelet analysis) analyses are conducted on daily mean values of the selected parameter. For trend and periodicity analyses and generating autocorrelogram, the missing daily mean values are imputed using the `forecast` package [@Hyndman:2008].

v) The `Linear Regression` tab allows a user to perform univariable and multiple linear regression analyses among the parameters of choice. For univariable linear regression, a scatter plot will be generated with least-squares linear fit. For multiple linear regression, multiple independent parameters can be selected. A scatter plot between the dependent variable and fitted data (using regression coefficients) will be generated. Relevant statistical coefficients are provided along with the plots.

vi) The `Compare` tab allows users to upload a second data file to compare data between the selected parameters from the two input files. The selected quality check criteria conditions applied on the parameters of the first input file will be automatically applied to the parameters in the second input file. Time-series, scatter, and diurnal variability plots of the two parameters of interest will be generated.

vii) Some features of the widely used `openair` package [@Carslaw:2012] are integrated into  `pollucheck` with permission. Calendar and time variation plots of the selected parameter are generated in this tab. Daily data will be used for Calendar plots and hourly data will be used for time variation plots.

An extensive list of frequently asked questions (FAQs) is provided as a separate tab for a better understanding of the `pollucheck` functioning, detailed features of the plots,  and analysis and the various packages used to build `pollucheck`.

# Limitations

1)  `pollucheck` does not download data automatically from the cloud. Downloaded files need to be provided as input.
2)  Multiple files cannot be uploaded to `pollucheck` at a given time.
3)  The current version of `pollucheck` is limited to accepting real-time data files from only three data sources.
4)  Some analyses (e.g., periodicity analysis) can be performed using daily mean values only.
5)  Caution needs to be exercised when using the averaged wind direction data, since wind direction is a vector quantity; hence it needs to be processed in a different way which has not been implemented here. Only wind direction data at one hour resolution is processed correctly.
6)  Any manipulation or alteration to the downloaded file before giving it as input to the app can lead to erroneous results.

# Installation

`pollucheck` is hosted online on *shinyapps.io* and can be installed to serve locally from [GitHub](https://github.com/).

Load and run `pollucheck` as follows:

``` {.r}
install.packages("devtools")
devtools::install_github("adithirgis/pollucheck")
pollucheck::pollucheck_run()
```

`pollucheck` is furnished with a preloaded data set for a quick user tour of the analysis, plotting options, and the functions available. In the `Compare` tab, the preloaded data set acts as the second input file if no second file is uploaded.

# Case Study

For better understanding of the major functionalities of `pollucheck`, we present a case study based on 18 months of pollution data set. This data is downloaded from the Central Pollution Control Board dashboard for the monitoring station located at Hebbal, Bengaluru, India at a time resolution of 60 minutes. Only plots related to PM~2.5~ data generated through the app are shown here. Figure 1 depicts the efficiency of the app in detecting the outliers. The top panel of Figure 1 shows the hourly time series of the raw PM~2.5~ (few outliers were synthetically added to the data), while the bottom panel depicts the quality checked data. Almost all the sporadically high values were detected by the app as outliers and removed.



![](inst/shiny/WWW/figure1a.JPG){ width=100% }


![Hourly time series of raw (top panel) and cleaned (bottom panel) PM~2.5~.](inst/shiny/WWW/figure1b.JPG){ width=100% }
  


Figure 2 depicts the difference between **Month and year box plot** and **Monthly box plot**. These plots are highly useful if the dataset length is more than a year. **Monthly box plot** (bottom panel) partitions all the data points into the calendar month bins irrespective of the year. While **Month and year box plot** (top panel) accounts for the entire timeline, i.e., including the year.

![](inst/shiny/WWW/figure2a.JPG){ width=100% }


![Box plots depicting the monthly variations in hourly PM~2.5~.](inst/shiny/WWW/figure2b.JPG){ width=100% }


Diurnal variation in PM~2.5~ based on mean (and standard deviation) and median (and interquartile range) are shown in the top and bottom panels of Figure 3, respectively. The choice between mean and median is useful when the distribution of the data deviates from normal. In the top panel, the line depicts the mean and the vertical bars depict standard deviation.  In the bottom panel, the line depicts the median and the vertical bars depict the interquartile range.


![](inst/shiny/WWW/figure3a.JPG){ width=100% }


![Diurnal variations in PM~2.5~.](inst/shiny/WWW/figure3b.JPG){ width=100% }

In Figure 4, a linear regression is shown between PM~2.5~ and the PM~2.5~/PM~10~ ratio.  The app computes the ratio using the individual PM~2.5~ and PM~10~ data sets. The blue line depicts the least square linear fit. The R-square and the equation of the linear fit are also provided on the panel.

![Linear regression analysis.](inst/shiny/WWW/figure4.JPG){ width=100% }


The periodicity in PM~2.5~ is shown as a wavelet periodogram (Figure 5). Wavelet analysis is useful in analysing non-stationary time series data. Only daily averaged data will be used for this analysis and missing data is imputed to perform the wavelet analysis. 

![PM~2.5~ periodicity analysis based on wavelet transform.](inst/shiny/WWW/figure5.JPG){ width=100% }


# Acknowledgements

We wish to thank Prof. Julian D. Marshall (University of Washington, Seattle), Prof. Joshua Apte (University of California, Berkeley), Dr. Jai Asundi (Center for Study of Science, Technology and Policy, Bengaluru), Dr Saumya Singh (University of California, Berkeley), and the R community for their help and support.

# References
# pollucheck: Open-Source Air Quality App!


[![R-CMD-check](https://github.com/adithirgis/pollucheck/workflows/R-CMD-check/badge.svg)](https://github.com/adithirgis/pollucheck/actions) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5128607.svg)](https://doi.org/10.5281/zenodo.5128607)[![DOI](https://joss.theoj.org/papers/10.21105/joss.03435/status.svg)](https://doi.org/10.21105/joss.03435)

pollucheck helps exploring the open-source air quality data.

-   pollucheck allows users to handle open-source air quality data sets available from [OpenAQ](https://openaq.org/#/countries), Central Pollution Control Board [CPCB](https://app.cpcbccr.com/ccr/#/caaqm-dashboard-all/caaqm-landing), and [AirNow](https://www.airnow.gov/international/us-embassies-and-consulates/).
-   Users can visualize data, analyze data, perform basic statistical operations, and generate a variety of publication-ready plots.
-   We have also included the popular [openair](https://cran.r-project.org/web/packages/openair/index.html) package in this app.
-   We have hosted this app here - https://aruapps.shinyapps.io/OpenSourceAirQualityApp/ 

### We are in this together!

A walk through to use this app for everyone -

#### Download data - CPCB website

-   Example - Where do you live in India?

-   Find the nearest [CPCB station](https://app.cpcbccr.com/ccr/#/caaqm-dashboard-all/caaqm-landing) to download data from a regulatory air quality monitor.

-   Visit [CPCB website](https://app.cpcbccr.com/ccr/#/caaqm-dashboard-all/caaqm-landing) to access the Central/State Pollution Control Board Data.

![](inst/shiny/WWW/CPCB_data_down_S.jpeg)

-   Select the Indian state from the "State Name" drop-down.

![](inst/shiny/WWW/CPCB_Station.jpeg)

-   Now select the city for which the data needs to be downloaded using the "City Name" drop-down menu.

![](inst/shiny/WWW/CPCB_Station_city.jpeg)

-   Now from the "Station Name" drop-down select the desired station.

-   Select the Parameters. Note- Multiple parameters can be selected at a time.

![](inst/shiny/WWW/CPCB_Station_parameters.jpeg)

-   Report Format- To use the pollucheck app, Please keep the format as "tabular".

-   Criteria- This drop-down will help you to select between different time averaging of data. Note- pollucheck app only supports 15 min, 30 min, and 60 min average data.

-   Select the Start Date and End date of the data and click on "Submit".

-   Download that data (15, 30, 60 min resolution would be good).

![](inst/shiny/WWW/CPCB_Station_TA.jpeg)

#### Download data - OpenAQ 

- Click on [OpenAQ](https://openaq.org/#/countries/IN?_k=5ecycz) in the app to download the OpenAQ data set.

<img src="inst/shiny/WWW/oaq_1.jpg" width="286"/>

- Now, click on the *Download* option on your browser.

<img src="inst/shiny/WWW/oaq_2.jpg" width="517"/>

- In data download, you can download the data by *Locations* or by *Datasets*.

<img src="inst/shiny/WWW/oaq_3.jpg" width="484"/>

- Now select the desired, *Country*, *City/Region* and *Location* using the drop-down.

<img src="inst/shiny/WWW/oaq_4.jpg" width="306"/><img src="inst/shiny/WWW/oaq_5.jpg" width="103"/><img src="inst/shiny/WWW/oaq_6.jpg" width="215"/>

- To proceed further select the *Start Date* and the *End Date* for the data.

<img src="inst/shiny/WWW/oaq_7.jpg" width="525"/>

- Now select which type of sensor is available at that location. Is it a *Low-cost Sensor* or *Reference Grade* monitor?

![](inst/shiny/WWW/oaq_8.jpg)

- Finally, select the parameters from the list of *Core* and *Additional Parameters* which you wish to download from that particular location and click on *Download Selection*. The file will be saved in your local disk in .csv format.

<img src="inst/shiny/WWW/oaq_9.jpg" width="519"/>

#### Download data - AirNow 

- Click on [AirNow-US Embassies](https://www.airnow.gov/international/us-embassies-and-consulates/) to visit official website to download the data.

<img src="inst/shiny/WWW/a_1.JPG" width="557"/>

- Select the desired city and the parameters.

<img src="inst/shiny/WWW/a_2.JPG" width="554"/>

- Now go to the *Historical* sub-menu on the Homepage and download the desired file.

<img src="inst/shiny/WWW/a_3.JPG" width="557"/>

### App usage

- Select the source from where the data was downloaded.

- Now select the time resolution at which the data was downloaded.

<img src="inst/shiny/WWW/image_1.JPG" width="315"/>

- Select the check box according to your need.

    -   Remove Negative values- Negative values do not represent concentration,they represent missing values, so it is always advised to remove them.This option helps you to remove all the negative values from your entire data set.

    -   Remove duplicate consecutive values- Sometimes when the instrument breaks down, it tends to show exactly same consecutive values, it is advised to remove these as well. This feature removes consecutive repetitive values in your data set.

    -   Specify a multiple (X) to remove outliers based on Mean and SD- If you want to clean your data set based on outliers, not usually necessary, use only if you want to remove outliers based on Mean and Standard Deviation values.

    -   Specify % of data completeness for computing daily mean values- If you are looking for entire/complete data set to be present for analysis and not less, you can use this to select the desired level of completeness in a day using the scroll bar.

    -   Remove PM2.5 and PM10 above- Usually, values above 9999 are incorrect, also because the instruments usually measure only to 999 values in PM instruments. This can be removed using this filter option.

<img src="inst/shiny/WWW/image_2.JPG" width="410"/>

- Output aggregation- The uploaded data can be converted into daily or hourly mean values.

- "Download as csv" or click on "Show Data" to see the data in the app.

![](inst/shiny/WWW/image_3.JPG)

- Look at the time series of pollutant concentrations in the **Plots** tab (time series are plots with x axis representing time). Do you see patterns? Are there times of the month or times of the day where concentrations are particularly higher or lower? Are there particular months in a year that are more polluted than others?
- Think about sources in the particular location: traffic, industries, garbage burning, etc.
- What more do you want to learn? Talk to the Humans of ILK.

### App tabs

##### File tab 

- Displays the data after cleaning process.

![](inst/shiny/WWW/image_4.JPG)

##### Summary tab 

- Displays the summary statistics for daily, monthly or for the entire data set.

![](inst/shiny/WWW/image_5.JPG)

##### Summary Plots tab 

- generates time series, box plot, and diurnal plot of the selected parameter.

- Data availability plot of all the pollutants after the cleaning process can be generated.

![](inst/shiny/WWW/DA_plot.JPG)

- The parameter to plot and the data aggregation options are available.

![](inst/shiny/WWW/image_6.JPG)

- Options to edit the Title and axis labels are available.

- Time-series plot

![](inst/shiny/WWW/image_7.JPG)

![](inst/shiny/WWW/image_8.JPG)

- Month and year box plot

![](inst/shiny/WWW/image_9.JPG)

- Monthly box plot

![](inst/shiny/WWW/image_10.JPG)

- Vertical bar plot

![](inst/shiny/WWW/image_11.JPG)

- Diurnal pot using hourly values - has two types using all data or distributed month wise. There is an option to plot point and bars as Median and IQR respectively or Mean and Standard Deviation. The data used for plotting can be downloaded as csv file.

![](inst/shiny/WWW/image_12.JPG)

![](inst/shiny/WWW/image_13.JPG)

##### Statistical Plots tab 

- Tests for normality, pattern and generates density plot, qq plot of the selected parameter.

- Using a selected parameter and aggregation methods, normality test using the Anderson Darling test (for N \> 500) or Shapiro-Wilk test can be conducted.

- Density plot

![](inst/shiny/WWW/image_14.JPG)

![](inst/shiny/WWW/image_15.JPG)

- Q-Q plot

![](inst/shiny/WWW/image_16.JPG)

![](inst/shiny/WWW/image_17.JPG)

- Trend Analysis is also available for daily values. For trend analysis using Mann-Kendall test we use [mk.test](https://www.rdocumentation.org/packages/trend/versions/1.1.4/topics/mk.test). For imputing values in the discontinuous data set we use [forecast package](https://cran.r-project.org/web/packages/forecast/forecast.pdf). For continuous wavelet transform we use [biwavelet package](https://cran.r-project.org/web/packages/biwavelet/biwavelet.pdf). In periodicity analysis, the contours covered by black lines represent the significant periodicity at 95% significant 519 level.


![](inst/shiny/WWW/image_30.JPG)

<img src="inst/shiny/WWW/image_31.JPG" width="590"/>

<img src="inst/shiny/WWW/image_18.JPG" width="590"/>

##### Linear Regression tab 

- There is an option of plotting linear regression plots between various parameters available.

![](inst/shiny/WWW/image_20.JPG)

![](inst/shiny/WWW/image_21.JPG)

- Also multi linear regression can be performed.

![](inst/shiny/WWW/image_22.JPG)

![](inst/shiny/WWW/image_23.JPG)

##### Compare tab 

- Allows user to upload data from another site for comparison and generate time series and a scatter plot between parameters selected from different sites.

- There are options to generate time series, scatter plot / linear regression and diurnal plots for both the sites.

![](inst/shiny/WWW/image_24.JPG)

<img src="inst/shiny/WWW/image_32.JPG" width="397"/>

![](inst/shiny/WWW/image_26.JPG)

##### `openair` tab 

- allows users use the package's widely used functions for the selected parameter.

- Calendar plot

![](inst/shiny/WWW/image_27.JPG)

![](inst/shiny/WWW/image_28.JPG)

- Time variation plot

![](inst/shiny/WWW/image_29.JPG)

## Installation

`pollucheck` is hosted online on *shinyapps.io* and can be installed to serve locally from [GitHub](https://github.com/).

Load and run `pollucheck` as follows:

``` {.r}
install.packages("devtools")
devtools::install_github("adithirgis/pollucheck")
pollucheck::pollucheck_run()
```

`pollucheck` is furnished with a preloaded data set for a quick user tour of the analysis, plotting options and the functions available. In the `Compare` tab, the preloaded data set acts as the second input file if no second file is uploaded.

## Community guidelines

1. [Contribute to the software](.github/CONTRIBUTING.md)

2. Report issues or problems with the software / Seek Support

- Please open an issue in the [issue tracker of the project.](https://github.com/adithirgis/pollucheck/issues)

3. Contributors must adhere to the [Code of Conduct](.github/CODE_OF_CONDUCT.md).

4. We are happy to incorporate more features based one what users need. Write to us at [contact\@ilklabs.com](mailto:contact@ilklabs.com).

## Author credit statement

Adithi R. Upadhya created the package, Meenakshi Kushwaha supervised and will maintain. Pratyush Agrawal contributed to designing, testing the app and data collection while Sreekanth Vakacherla designed and supervised. All the authors contributed to the manuscript. 
# pollucheck 1.0                     

## Major changes 

## Bug fixes# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
# Contributing to pollucheck

This outlines how to propose a change to pollucheck. 
For more detailed info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib). 

## Fixing typos

You can refer to [this](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork) for proposing changes with pull request from a fork. 
You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 

If you find any typos in an `.Rd` file, then please make changes in the corresponding `.R` file, as `.Rd` files are automatically generated by [roxygen2](https://roxygen2.r-lib.org/articles/roxygen2.html) and should not be edited by hand.

## Bigger changes

The first step here for you to make any changes is to install `devtools` using `install.packages("devtools")`.
If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. We recommend using `usethis::create_from_github("adithirgis/pollucheck", fork = TRUE)`. For details on `usethis` click [here](https://cran.r-project.org/web/packages/usethis/index.html).

*   Install all development dependencies with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. It should pass without errors and warnings. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the pollucheck project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

# Frequently Asked Questions

### File

**What does the app do?**  
*This app helps to analyze and visualize open source air quality data
available. An example dataset is already loaded to help you walk through
all the features of the app. If you need help with downloading your own,
go to the next question. This app can process all parameters except Wind
Direction which when downloaded at 1 hour is processed correctly, any
other time-resolution is used then the app will not be process wind
direction correctly.*

**What is ratio\_PM?**  
*ratio\_PM is = PM2.5 / PM10 or ratio of PM2.5 and PM10.*

**What parameters are applied to pre-loaded data?**  
*The pre-loaded data, downloaded from
[CPCB](https://app.cpcbccr.com/ccr/#/caaqm-dashboard-all/caaqm-landing)
has all the conditions applied ie: removed negative values and repeated
consecutive measurements, removed outliers using mean and 3 times std
dev and also checked for 75% of completeness of data in a day for all
parameters.*

**Where can we download the air quality e data from?**  
*Air Quality data from India can be downloaded from 3 sources -
[CPCB](https://app.cpcbccr.com/ccr/#/caaqm-dashboard-all/caaqm-landing),
[OpenAQ](https://openaq.org/#/countries/IN?_k=5ecycz), and
[AirNow](https://www.airnow.gov/international/us-embassies-and-consulates/#India).*

**Can I use this app to analyze app work for data from countries other
than India?**  
*Yes, two of the sources mentioned above OpenAQ, and AirNow have data
sets from many other countries. CPCB is specific to India. Make sure you
specify the correct timezone in the app when analyzing data from
different countries.*

**How to use the app?**  
*Please check [here](https://github.com/adithirgis/pollucheck) for a
walk through.*

**Will this application read data from a reference monitor which is not
from the sources mentioned in the app?**  
*No, this app can only do post-processing of the data downloaded from
the mentioned sources.*

**What are the different file formats which can be uploaded to this
app?**  
*The data can be uploaded only in .csv or .xlsx format.*

**What if I do not use any quality check options?**  
*Data from sources specified here, sometimes contains fill values which
can be due to instrument malfunction or unavailability of data. The data
will be as it is except that it will remove fill values from the data.*

**What is the naming convention of downloaded files?**  
The File tab downloaded data is the hourly or daily average so
\_average.csv, while the summary statistics file has an extension of
\_summary.csv, and the diurnal table download will have a suffix of
\_diurnal.csv.

**What are the different packages used in the app?**  
*The different packages we use here are - tidyverse, ggplot2, openair,
lubridate, shiny, bslib, forecast, biwavelet, readxl, DT, data.table,
nortest, janitor, zoo.*

**Can the input file be manipulated/edited?**  
*No. This app is currently designed to read the preset format of the
input file. Any alterations in the input file may lead to unpredictable
results.*

**What is the standard outlier condition?**  
*The standard outlier detection condition is to calculate the daily mean
and standard deviation and then removing values that are lower than Mean
subtracted by x times SD and greater than Mean added to x times SD,
where x is any real number.*

**Why should we remove repeated consecutive measurements?**  
*Sometimes when the instrument breaks down, it tends to show exactly the
same consecutive values (often this value is the last measured value),
it is suggested that such consecutive values be removed from the
dataset*

**Why should we remove values above 9999 from the data?**  
*Usually, values above 9999 are incorrect, also because the instruments
usually measure only up to 999 values in PM instruments. This can be
removed using this filter option, along with this 9999 could also be a
fill value.*

**What are the available time resolutions for plots and tables?**  
*Only hourly and daily average plots can be generated in this
application, although the input data can be of 60, 30, and 15 minutes
time resolution.*

**How to export plots?**  
*To save the image right click on the plot and select save image as. You
can save it on your local computer in .png or .jpeg format.*

**Can we upload data sets of different timezone?**  
*Yes, for data from OpenAQ and AirNow, different time formats are
supported, just select the right timezone.*

**Are the plots customizable?**  
*Only the axis labels and Plot title can be customized.*

**Whom to contact for doubts or suggestions?**  
*We are happy to incorporate more features based on what users need.
Write to us at <contact@ilklabs.com>. For reporting any bug/issue click
[here](https://github.com/adithirgis/pollucheck/issues).*

**Can the plots generated using pollucheck be used for publications?**  
*Depends on the user. If you find the resolution and everything usable,
please go ahead. Remember to cite us!*

**How to cite pollucheck?**  
*We will update this soon!*

**Is the theme of the application customizable?**  
*No, the theme of the app is not customizable, if you have any
suggestions and resources, please Pull Request in the github repo link
provided above.*

### Summary

**What is the naming convention of the file downloaded?**  
*The data downloaded in this tab will have a suffix of \_average.csv.*

### Summary Plots

**What does data availability plot indicate?**  
*The data availability plot shows the available percentage of daily data
points for each parameter.*

**How to read the box plots?**  
*A boxplot is a way to show a [five-number
summary](https://www.statisticshowto.com/statistics-basics/how-to-find-a-five-number-summary-in-statistics/)
in a chart. The main part of the chart (the “box”) shows where the
middle portion of the data is: the interquartile range. At the ends of
the box, you’ll find the first
[quartile](https://www.statisticshowto.com/what-are-quartiles/)(the 25%
mark) and the third quartile (the 75% mark). The whiskers on the far
most of the box plots show the minimum and maximum values in the data
sets. Finally, the median is represented by a vertical bar in the center
of the box.
[source](https://www.statisticshowto.com/probability-and-statistics/descriptive-statistics/box-plot/)*

**What do the bars in the vertical bar plots indicate?**  
*The vertical bars represent mean + sd and mean - sd.*

**What is the difference between month-yearly plot and monthly plot?**  
*The month-yearly plot displays box plot based on each month and year,
while monthly plot displays the plot for all months in multiple years.*

**What are the bars in the diurnal plot and why are there two
options?**  
*Diurnal plots represent the time of the day variability using mean (of
that hour) and sd or median (of that hour) and IQR (InterQuartile
Range). So when you choose mean and sd the bars represent sd, while when
you choose median and IQR, the bars represent IQR.*

**Can we plot ratios of different pollutants on a time series?**  
*No, this feature is not implemented yet.*

### Statistical Plots

**What is the difference between the Shapiro-Wilk normality test and the
Anderson- Darling normality test?**  
*Shapiro-Wilk test is used to check for normality when the sample size
is below 5000, while Anderson- Darling test is used when sample size is
greater than 5000 if you have data points above 5000, please check
Anderson-Darling test.*

**How is Q-Q Plot different from a Scatter Plot?**  
*A Q-Q plot is a scatterplot created by plotting two sets of quantiles
against one another. If both sets of quantiles come from the same
distribution, we should see the points forming a line that’s roughly
straight.
([Source](https://data.library.virginia.edu/understanding-q-q-plots/#:~:text=A%20Q%2DQ%20plot%20is%20a,truly%20come%20from%20Normal%20distributions.))*

**How is the data imputed for trend analysis?**  
*For imputing values in the discontinuous data set we use [forecast
package](https://cran.r-project.org/web/packages/forecast/forecast.pdf),
please check the resource for more details.*

**What packages were used for Trend Analysis?**  
*Trend Analysis is available for daily values alone. For trend analysis
using Mann-Kendall test we use
[mk.test](https://www.rdocumentation.org/packages/trend/versions/1.1.4/topics/mk.test).
For imputing values in the discontinuous data set we use [forecast
package](https://cran.r-project.org/web/packages/forecast/forecast.pdf).
For continuous wavelet transform we use [biwavelet
package](https://cran.r-project.org/web/packages/biwavelet/biwavelet.pdf).
In periodicity analysis, the contours covered by black lines represent
the significant periodicity at 95% significant 519 level.*

**What are the time resolutions used for Trend Analysis?**  
*Daily mean and monthly mean is used to perform Trend Analysis by
Mann-Kendall test.*

**What is Autocorrelogram plot?**  
*Autocorrelations or lagged correlations are used to assess whether a
time series is dependent on its past. These are generated using monthly
mean values. To have continuous data we imputed using the same method
discussed above.
[source](https://www.rdocumentation.org/packages/stats/versions/3.3.1/topics/acf)*

**What are the blue dashed lines in the Autocorrelogram plot?**  
*The lines give the values beyond which the autocorrelations are
(statistically) significantly different from zero.
[source](https://stats.stackexchange.com/questions/49571/understanding-the-blue-dotted-lines-in-an-acf-from-r#:~:text=The%20lines%20give%20the%20values,ACF%20seems%20to%20indicate%20seasonality.)*

**What are the unit of the Lag in the Autocorrelogram plot?**  
*The lag is effectively a time lag and is dependent of the frequency you
choose to build the time series. Each bar on the graph correspond to one
time unit. The unit here is month for Lag in x axis.*

**What are the time resolutions used for Periodicity Analysis?**  
*Daily mean data is used to generate the biwavelet. We have used the
default options in the biwavelet package to generate the plot.*

### Linear Regression

**What are Linear and Multi linear regression?**  
*Linear regression is which consist of single independent variable while
Multi linear regression consist of three or more variables.*

**What is the difference between multiple R- square and adjusted R-
square?**  
*Multiple R squared is simply a measure of Rsquared for models that have
multiple predictor variables. Therefore it measures the amount of
variation in the response variable that can be explained by the
predictor variables. The fundamental point is that when you add
predictors to your model, the multiple Rsquared will always increase, as
a predictor will always explain some portion of the variance. Adjusted
Rsquared controls against this increase, and adds penalties for the
number of predictors in the model. Therefore it shows a balance between
the most parsimonious model, and the best fitting model. Generally, if
you have a large difference between your multiple and your adjusted
Rsquared that indicates you may have overfit your model.
[source](https://stats.stackexchange.com/questions/241283/what-is-the-main-difference-between-multiple-r-squared-and-adjusted-r-squared/241298)*

**How to check if the data is significant or not?**  
*The p-value in the last column tells you the significance of the
regression coefficient for a given parameter. If the p-value is small
enough to claim statistical significance, that means there is strong
evidence that the coefficient is different from 0.
[source](https://stats.stackexchange.com/questions/37912/how-to-determine-which-variables-are-statistically-significant-in-multiple-regre)*

**What are residuals?**  
*The data points usually don’t fall exactly on the regression equation
line; they are scattered around. A residual is the vertical distance
between a data point and the regression line.
[source](https://www.statisticshowto.com/residual/)*

**What are degrees of freedom?**  
*The degrees of freedom indicate the number of independent values that
can vary in an analysis without breaking any constraints.
[source](https://statisticsbyjim.com/hypothesis-testing/degrees-freedom-statistics/)*

### Compare

**What conditions are applied on the data for comparison?**  
*The same conditions applied for the original data uploaded in **File**
tab will be applied example - remove negative values, or averaging
period of hourly or daily.*

**Can I compare two data sets of different timezone?**  
*No, right now the app allows for comparison with data from the same
time zone.*

**Can I compare two data sets of different time periods?**  
*Yes, you can compare the time series of both the data sets, but the
scatter plot will not be generated as they belong to different time
series.*

**Can I compare data sets from different downloading sources (different
websites)?**  
*Yes, you can definitely compare from different available sources.*

**How many data sets can be compared at a time in the app?**  
*As of now, only one single file can be compared with the uploaded
dataset.*

**Is it possible to compare two data sets with different time
resolutions?**  
*Yes! It is possible to compare data sets with different time
resolutions since the application averages to 1 hour or daily average.*

### openair

**What is openair?**  
*openair is an R package developed for the purpose of analysing air
quality data - or more generally atmospheric composition data. We use
functions from openair to plot. The openair package usage and other
details can be found [here](https://github.com/davidcarslaw/openair).*

**What is the difference between “median and quantiles” and “mean and
95% confidence intervals” in Time Variation Plot?**  
*Since the Time Variation plot shows the variation in parameters based
on different time resolutions. It is the same option as time variation
plots to plot mean and sd or median and IQR.*

**How many different plots can be created using openair?**  
*Right now, this app supports only two plots - calendar plot and
temporal variation plot.*

**Where can I find the documentation of openair?**  
*Please click
[here](https://cran.r-project.org/web/packages/openair/openair.pdf).
Also, remember to cite them!*

**What are the ambient levels in India?**  
*Please check the
[link](https://app.cpcbccr.com/ccr_docs/FINAL-REPORT_AQI_.pdf) for
details.*
