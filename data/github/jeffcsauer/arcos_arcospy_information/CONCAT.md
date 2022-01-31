---
title: 'arcos and arcospy: R and Python packages for accessing the DEA ARCOS database from 2006 - 2014'
tags:
  - R
  - Python
  - API
  - Open science
  - Health
authors:
  - name: Steven Rich
    affiliation: 1
  - name: Andrew Ba Tran
    affiliation: 1
  - name: Aaron Williams
    affiliation: 1
  - name: Jason Holt
    affiliation: 1    
  - name: Jeffery Sauer  
    orcid: 0000-0002-8581-6919
    affiliation: 2
  - name: Taylor M. Oshan
    affiliation: 2
affiliations:
 - name: The Washington Post
   index: 1
 - name: University of Maryland, College Park, Department of Geographical Sciences
   index: 2
date: 09 July 2020
bibliography: paper.bib
---

# Summary

In the early 2000s governmental agencies across the United States began to observe increases in the number of all-cause opioid-involved deaths [@rudd]. The Centers for Disease Control (CDC) describe this Opioid Overdose Epidemic as occuring in three waves, with the first wave attributed to the widespread distribution of prescription opioids, reaching an opioid prescribing rate as high as 81.3 per 100 persons in 2012 [@cdcweb]. While the more recent wave two and wave three of the Opioid Overdose Epidemic are largely attributed to illicit opioids like heroin and fentanyl, prescribing rates and prescription misuse remain high [@Guy2017;@ko2020]. 

Researchers, journalists, and government agencies are still actively investigating the myriad impacts of prescription opioids on the Opioid Overdose Epidemic. One powerful tool for understanding prescription opioid distribution is the Drug Enforcement Administration's (DEA) Automation of Reports and Consolidated Orders System (ARCOS). ARCOS tracks the commercial distribution of controlled substances in the United States, including opioid analgesics. ARCOS data is highly detailed, tracking commercial origin, pharmacy order frequency, point-of-sale distribution, and more. For a variety of reasons ranging from patient confidentiality to protecting trade secrets, access to sub-state ARCOS data is available only for approved requests (e.g. research or litigation) [@arcosprivacy]. Recent litigation efforts by *The Washington Post*, HD Media, and local journalists allowed for the public release of an anonymized, large portion of the ARCOS database from 2006 to 2012, with additional data for 2013 and 2014 now also available. `arcos` and `arcospy` are open-source API wrappers in `R` and `Python`, respectively, that allow researchers and interested citizens to easily access this newly available portion of the ARCOS database.

## Statement of Need

Previously, researchers wanting to use ARCOS data relied on what was made available by the DEA, typically in the form of state-level estimates, or submitted special access requests to the DEA [@Kenan2012;@reisman_2009]. While alternative data on prescription records are offered by the Centers for Medicare & Medicaid Services in the [Medicare Provider Utilization and Payment Datasets](https://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/Medicare-Provider-Charge-Data), this data pertains to a specific sample of the population and spans a different set of years (2011 to 2017). In addition, the level of detail in the ARCOS data offers substantial opportunity for commodity analysis about prescription opioids including market dynamics, product demand, supply chain flow, and more. This relationship between commercial distribution and the Opioid Overdose Epidemic is an important area for future study to define and recognize the warning signs of potentially problematic prescribing practices [@VanZee2009]. The release of national, longitudinal, sub-state ARCOS data is a major contribution for researchers interested in the distribution of prescription opioids and the subsequent sociomedical impacts.

In raw format, the ARCOS database is more than 130 gigabytes and includes several hundred columns. Thus, the purpose of `arcos` and `arcospy` are meant to:

-   Simplify access to an open, large, robust prescription opioid database
-   Provide measures of prescription opioid distribution relevant to both the medical and social sciences 
-   Promote analytical flexibility and reproducibility through mirrored functionality across `R` and `Python` 

## API Structure

The `arcos` and `arcospy` API is [publically available](https://arcos-api.ext.nile.works/__swagger__/) and hosted using the [OpenAPI specification](https://swagger.io/specification/). The primary maintainers of API database are members of the Data Reporting Team at *The Washington Post*. A key is required to use the API. The standard key is `WaPo` and additional keys may be sourced from [the Github repository](https://github.com/wpinvestigative/arcos). Guidelines on using the API are available from [*The Washington Post*](https://www.washingtonpost.com/national/2019/07/18/how-download-use-dea-pain-pills-database/).

All commands share the same name between `arcos` and `arcospy`. This allows users to easily switch between languages if the need arises. Outputs from all of the functions are delivered in popular formats - `data.frames` in `R` and `pandas.DataFrame` in `python` - to enable statistical, spatial, network, or other types of analysis.

Both `arcos` and `arcospy` use parameter delivery - `urltools` in `R` and `requests` in `Python` - to build the API query. Checks are in place to ensure that invalid inputs are not passed to the API. For example, a series of integers cannot be passed as a county name. Corrective warning messages are returned to users who provide invalid inputs.

## Data Availability and Basic Usage

 Data can be gathered at the pharmacy, distributor, county, or state as the geographic unit of analysis. Depending on the geographic level, there may be raw, summarized, or supplemental data available. For example, the `county_raw()` command returns each individual ARCOS record for a given county from 2006 to 2014. The following code chunk demonstrates this function in `R`:

 ```{.r}
library(arcos)
# Gather all ARCOS records for Hill County, Montana
HillRaw <- county_raw(county = "Hill", 
                      state = "MT", 
                      key = "WaPo")
head(HillRaw)
```

```
| REPORTER_DEA_NO | REPORTER_BUS_ACT | REPORTER_NAME        | ... | dos_str |
|-----------------|------------------|----------------------|-----|---------|
| PM0023046       | DISTRIBUTOR      | MCKESSON CORPORATION | ... | 5.0     |
| PM0023046       | DISTRIBUTOR      | MCKESSON CORPORATION | ... | 5.0     |
| PM0023046       | DISTRIBUTOR      | MCKESSON CORPORATION | ... | 5.0     |
| PM0023046       | DISTRIBUTOR      | MCKESSON CORPORATION | ... | 5.0     |
| PM0023046       | DISTRIBUTOR      | MCKESSON CORPORATION | ... | 7.5     |
| PM0023046       | DISTRIBUTOR      | MCKESSON CORPORATION | ... | 20.0    |
```

 However, the `summarized_county_annual()` command returns the annual summarized totals for a given county for each year of 2006 to 2012. The following code chunk demonstrates this function in `Python`:

 ```{.python}
from arcos import summarized_county_annual
# Gather summarized ARCOS records for Hill County, Montana
HillSummarized = summarized_county_annual(county = "Hill", 
                                          state = "MT", 
                                          key = "WaPo")
HillSummarized.head()
```
```
 	BUYER_COUNTY 	BUYER_STATE  year 	count 	DOSAGE_UNIT 	countyfips
0 	HILL 	        MT 	         2006 	1516 	594700 	        30041
1 	HILL 	        MT 	         2007 	1710 	505430 	        30041
2 	HILL 	        MT 	         2008 	2467 	715560 	        30041
3 	HILL 	        MT 	         2009 	3200 	851560 	        30041
4 	HILL 	        MT 	         2010 	3290 	803760 	        30041

 ```

 Given the number of records, users should anticipate that commands querying for raw data will take longer than commands querying for summarized data. Full documentation on how data is collected by the DEA is available in the [ARCOS Registrant Handbook](https://www.deadiversion.usdoj.gov/arcos/handbook/full.pdf). `arcos` and `arcospy` also include supplemental commands that return relevant auxiliary data - such as county population - gathered from the American Community Survey. A description of each of the functions currently offered, as well as examples in `R` and `Python` demonstrating functionality, are available on [the shared `arcos` and `arcospy` repository](https://github.com/jeffcsauer/arcos_arcospy_information).

There are several ways to conceptualize the unit of analysis for opioids from the present data. These include the total number of records, the total number of all opioid pills, the total number of specific opioid pills (i.e. oxycodone versus hydrocodone), or the total amount (in weight) of all or specific opioid pills. Other common units of analysis that may be of interest include morphine milligram equivalents (MMEs) or prescription counts [@stopka]. Users should choose a unit of analysis that has precedent in their discipline and take appropriate steps to standardize the data (e.g. by population or another stratum) when necessary.

## Conclusion

`arcos` and `arcospy` allows access to a substantial amount of previously unavailable data on prescription opioid distribution in the United States during the years leading up to the present Opioid Crisis. Data from the DEA ARCOS system has been used in scientific publications, primarily at the intersection of health and criminology, to investigate trends in analgesic use and potential abuse [@gilson; @joranson]. Additionally, the data made available by `arcos` has been used extensively by journalists at *The Washington Post* and local news outlets to report on trends in prescription opioid distribution [@diez_2019; @top_2020]. ARCOS data can be merged (non-spatially or spatially) with other United States statistical products through packages like ``tidycensus`` in ``R`` and ``cenpy`` in ``Python``, opening numerous doors for research and teaching exercises. Examples of these possibilities are available on the Github repositories for [``arcos``](https://github.com/wpinvestigative/arcos) and [``arcospy``](https://github.com/wpinvestigative/arcos). The flexibility to query the ARCOS DEA database using commands of the same name enhances reproducibility across languages and ease of access. Expanding the ways in which researchers and journalists can analyze robust datasets like ARCOS is an important step towards understanding how the United States arrived at the present Opioid Overdose Epidemic.

# Availability

``arcos`` is available on [CRAN](https://cran.r-project.org/web/packages/arcos/index.html) as well as [Github](https://github.com/wpinvestigative/arcos).

``arcospy`` is available on [PyPI](https://pypi.org/project/arcospy/) as a pip installable package as well as [Github](https://github.com/jeffcsauer/arcospy).

The repository for this article and additional information is stored on the shared the shared `arcos` and `arcospy` repository on [Github](https://github.com/jeffcsauer/arcos_arcospy_information).

# Citations
# Contributing

Contributions are welcome to both `arcos` and `arcospy`! For major contributions, please fork the master `arcos` or `arcospy` branch and open a pull request with the suggested changes. The appropriate team will then review the changes and determine if there is a generalizable solution to both `R` and `Python`. If there is no generalizable solution we will strive to make your contribution visible on the appropriate Github page.

Issues posted to the specific `arcos` and `arcospy` repositories will be cross-posted and centralized on the `arcos_arcospy_information` so that users in both languages can be made aware of happenings in the software. These issues will be adapted into issue formatted Github templates by the `arcos` and `arcospy` maintainers. Users may also post issues directly on the `arcos_arcospy_information` repository.

## Documentation

Improvements to the documentation for `arcos` and `arcospy` are welcome. As stated above, we will try to generalize all contributions to both packages. For example, if the description or working for a specific command is unclear we can improve in the documentation in both packages. Additionally, if there are features of a command that you think should be included in the documentation, please let us know so we can improve the user experience.

## Contributing to the core dataset

Presently, the core functionality of the data and API is maintained by the Data Reporting Team at *The Washington Post*. Questions and suggestions regarding the underlying data may be posted on either Github repository, although additional communication with members of the Data Reporting Team at *The Washington Post* may be necessary. 

## Contributing to the API functionality

There is ample room for users to suggest functions that can be added to `arcos` and `arcospy`. As stated above, the API itself is maintained by the Data Reporting Team at *The Washington Post*. Well-described functions with clear use cases can be posted to Github and then considered for integration by the Data Reporting Team at *The Washington Post*. 

It should be noted that existing functions in the `arcos` and `arcospy` functions have flexibility to be used in other functions and loops. Please take this into consideration before requesting a new API command. For example, users interested in downloading data for the entire nation at the county level might consider downloading all of the data (using the `raw_data()` command) and then processing it to own their needs locally. Alternatively, users could also loop the `county_summarized_annual()` over a nationwide FIPS list.

## Issues

Trackable issue pages are available on both the `arcos` and `arcospy` Github pages. Issues may be related to anything from malfunctioning commands to inconsistent data. We encourage an active discussion and hope to readily address any errors. We recommend that when submitting an issue users provide specific tags and examples of the aberrant behavior. If you are able to solve an issue with the code independently, please open a pull request with the corrected code and a short explanation for the bug fix.
## arcos (R) and arcospy (Python): twin software packages to access the DEA ARCOS database from 2006 to 2014
Welcome to an informational landing page for the R `arcos` and Python `arcospy` software. This landing page acts as a home base for background information on the software and a high-level overview of its features. Specific information on the R or Python versions can be found at their Github repositories.

[**arcos** repository](https://github.com/wpinvestigative/arcos)

[**arcospy** repository](https://github.com/jeffcsauer/arcospy)

**We recommend reading [the vignette](https://github.com/jeffcsauer/arcos_arcospy_information/tree/master/docs/R_Python_arcos_vignette.html) that provides an overview of the API, a quick start guide, and demonstration of R/Python dual functionality!**

## Motivation
The ongoing Opioid Crisis in the United States poses serious public health issues. Understanding *how* the United States arrived at the present crisis is crucial to avoid future crises. One powerful tool for understanding trends in prescription opioid distribution is the Drug Enforcement Agency's (DEA) Automation of Reports and Consolidated Orders System (ARCOS). In raw format, the available portion of the ARCOS database is more than 150 gigabytes and includes several hundred columns. Thus, `arcos` and `arcospy` are meant to simplify access to the available data so researchers and interested citizens can rapidly gather relevant measures of prescription opioid distribution. The data available in `arcos` and `arcospy` spans the years 2006 to 2014. These are the years leading up to the now widely recognized Opioid Crisis. The data has been featured in a series of investigative reporting articles by both national and regional newspapers. There are numerous potential applications of this data to the study of health sciences, criminology, medical sociology, and more.  

## Development history
`arcos` is the result of efforts by the data reporting team at *The Washington Post* to make a large portion of the DEA ARCOS more easily accessible to the public. `arcos` was released on CRAN in August of 2019 (note: `arcos` has been temporarily removed from CRAN due to COVID-19, there is intention to get the package back on CRAN as soon as possible - it remains `devtools` installable). Shortly thereafer, health geographers working at the University of Maryland translated the package to python under the name `arcospy`. Due to this development timeline and the fact that `arcospy` is a translation of `arcos`, at times we defer to the built-out documentation of `arcos` to avoid redundancy.

## Build status
Build status is evaluated via CRAN for `arcos` and travis for `arcospy`.

**arcos**: ![CRAN Badge](http://www.r-pkg.org/badges/version/arcos)

**arcospy**: ![Build Status](https://travis-ci.com/jeffcsauer/arcospy.svg?token=sRx5dHJBVzwnJnFuh3p9&branch=master)

## Dependencies

**arcos**: R (≥ 3.3.0), stringr, magrittr, jsonlite, dplyr, urltools, vroom

**arcospy**: Python (≥ 3.0.0), pandas, urltools

## Features

**arcos** and **arcospy** offers nearly 30 functions for users to access DEA ARCOS data using pharmacies, distributors, counties, or states as the unit of analysis, as well as useful supplementary information. These functions have the exact same name in both `R` and `Python`, allowing users to rapidly switch between languages if the need arises (i.e. if a certain type of analysis is available only in R or Python).


__Functions and the available datasets (Read the [reference page](https://wpinvestigative.github.io/arcos/reference/index.html) for more info):__

| Function                                                                  | What                                                                            | Type         | Years       | Drugs                   | Buyers                                             |
|---------------------------------------------------------------------------|---------------------------------------------------------------------------------|--------------|-------------|-------------------------|----------------------------------------------------|
| [buyer_addresses()](https://wpinvestigative.github.io/arcos/reference/buyer_addresses.html)                       | Get DEA designated addresses for each pharmacy                                  | Raw          |             | All                     | All                                                |
| [buyer_details()](https://wpinvestigative.github.io/arcos/reference/buyer_details.html)                           | Get monthly summarized pill totals by county                                    | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [buyer_list()](https://wpinvestigative.github.io/arcos/reference/buyer_list.html)                                 | Get list of business types listed in the BUYER_BUS_ACT in the ARCOS database    | Raw          | 2006 - 2014 |                         | All                                                |
| [combined_buyer_annual()](https://wpinvestigative.github.io/arcos/reference/combined_buyer_annual.html)           | Get annual total pills for each buyer (pharmacy, etc) in a county               | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [combined_buyer_monthly()](https://wpinvestigative.github.io/arcos/reference/combined_buyer_monthly.html)         | Get annual total pills for each buyer (pharmacy, etc) in a county               | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [county_list()](https://wpinvestigative.github.io/arcos/reference/county_list.html)                               | Get list of counties and states and fips codes represented in the ARCOS data    | Raw          | 2006 - 2014 |                         |                                                    |
| [county_population()](https://wpinvestigative.github.io/arcos/reference/county_population.html)                   | Get annual population for counties between 2006 and 2012                        | Supplemental | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [county_raw()](https://wpinvestigative.github.io/arcos/reference/county_raw.html)                                 | Download raw prescription data for specified county (by state and county names) | Raw          |             | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [county_raw_fips()](https://wpinvestigative.github.io/arcos/reference/county_raw_fips.html)                       | Download raw prescription data for specified county (by county FIPS code)       | Raw          |             | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [drug_county_biz()](https://wpinvestigative.github.io/arcos/reference/drug_county_biz.html)                       | Raw data by county and individual drug and business type                        | Raw          | 2006 - 2014 | All                     | All                                                |
| [drug_county_raw()](https://wpinvestigative.github.io/arcos/reference/drug_county_raw.html)                       | Raw data by county and individual drug and business type via fips code          | Raw          | 2006 - 2014 | All                     | All                                                |
| [drug_list()](https://wpinvestigative.github.io/arcos/reference/drug_list.html)                                   | Get list of the 14 drugs tracked in the ARCOS data                              | Raw          | 2006 - 2014 | All                     |                                                    |
| [not_pharmacies()](https://wpinvestigative.github.io/arcos/reference/not_pharmacies.html)                         | Get list of misidentified pharmacies by BUYER_DEA_NOs                           | Supplemental | 2006 - 2012 |                         | Retail Pharmacy, Chain Pharmacy                    |
| [pharm_cbsa()](https://wpinvestigative.github.io/arcos/reference/pharm_cbsa.html)                                 | Get the core-based statistical area GEOID for each pharmacy                     | Supplemental | 2006 - 2014 |                         | Retail Pharmacy, Chain Pharmacy                    |
| [pharm_counties()](https://wpinvestigative.github.io/arcos/reference/pharm_counties.html)                         | Get county GEOID for each pharmacy                                              | Supplemental | 2006 - 2014 |                         | Retail Pharmacy, Chain Pharmacy                    |
| [pharm_latlon()](https://wpinvestigative.github.io/arcos/reference/pharm_latlon.html)                             | Get latitude and longitude data for each pharmacy                               | Supplemental | 2006 - 2014 |                         | Retail Pharmacy, Chain Pharmacy                    |
| [pharm_tracts()](https://wpinvestigative.github.io/arcos/reference/pharm_tracts.html)                             | Get census tract GEOID for each pharmacy                                        | Supplemental | 2006 - 2014 |                         | Retail Pharmacy, Chain Pharmacy                    |
| [pharmacy_raw()](https://wpinvestigative.github.io/arcos/reference/pharmacy_raw.html)                             | Download raw prescription data for specified pharmacy into R                    | Raw          | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy                    |
| [raw_data()](https://wpinvestigative.github.io/arcos/reference/raw_data.html)                                     | Download raw ARCOS data                                                         | Raw          | 2006 - 2014 | All                     | All                                                |
| [reporter_addresses()](https://wpinvestigative.github.io/arcos/reference/reporter_addresses.html)                 | Get DEA designated addresses for each Reporter                                  | Raw          | 2006 - 2014 | All                     | All                                                |
| [state_population()](https://wpinvestigative.github.io/arcos/reference/state_population.html)                     | Get annual population for states between 2006 and 2014                          | Supplemental | 2006 - 2014 |                         |                                                    |
| [summarized_county_annual()](https://wpinvestigative.github.io/arcos/reference/summarized_county_annual.html)     | Get annual summarized pill totals by county                                     | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [summarized_county_monthly()](https://wpinvestigative.github.io/arcos/reference/summarized_county_monthly.html)   | Get monthly summarized pill totals by county                                    | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [total_distributors_county()](https://wpinvestigative.github.io/arcos/reference/total_distributors_county.html)   | Get total pills for each distributor in a county                                | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone |                                                    |
| [total_distributors_state()](https://wpinvestigative.github.io/arcos/reference/total_distributors_state.html)     | Get total pills for each distributor in a state                                 | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone |                                                    |
| [total_manufacturers_county()](https://wpinvestigative.github.io/arcos/reference/total_manufacturers_county.html) | Get total pills for each manufacturer in a county                               | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone |                                                    |
| [total_manufacturers_state()](https://wpinvestigative.github.io/arcos/reference/total_manufacturers_state.html)   | Get total pills for each manufacturer in a state                                | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone |                                                    |
| [total_pharmacies_county()](https://wpinvestigative.github.io/arcos/reference/total_pharmacies_county.html)       | Get total pills for each pharmacy in a county                                   | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |
| [total_pharmacies_state()](https://wpinvestigative.github.io/arcos/reference/total_pharmacies_state.html)         | Get total pills for each pharmacy in a state                                    | Summarized   | 2006 - 2014 | Oxycodone & Hydrocodone | Retail Pharmacy, Chain Pharmacy, and Practitioners |

## Examples

Several examples of using both `arcos` and `arcospy` are available on their respective Github pages. These examples include querying the data, combining the data with other census products at different geographic levels, and making the data spatial. Additional examples are provided to show how users can make use of base commands in loops to expand the functionality of the software.

[**arcos** examples](https://github.com/wpinvestigative/arcos/tree/master/vignettes)

[**arcospy** examples](https://github.com/jeffcsauer/arcospy/tree/master/demos)

## Installation

**R: installing arcos**

```R
#Get the latest stable release from CRAN:
install.packages("arcos")

#...or via devtools:

# install.packages("devtools")
devtools::install_github('wpinvestigative/arcos')
```

**Python: installing arcospy**

`arcospy` is available via PyPI:

```python
pip install arcospy
```

## Contribute

Contributing guidelines are outlined in [contributing.md](https://github.com/jeffcsauer/arcos_arcospy_information/blob/master/contributing.md).

## Credits

**arcos**: Steven Rich (The Washington Post), Andrew Ba Tran (The Washington Post), Aaron Williams (The Washington Post), Jason Hold (The Washington Post)

**arcospy**: Jeffery Sauer (University of Maryland, College Park), Dr. Taylor Oshan (University of Maryland, College Park)

## License

MIT (2019); The Washington Post and The Charleston Gazette-Mail (2019)
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

Thank you for taking the time to report a bug for the `arcos` and `arcospy` software! Please fill out as much of the following template as possible. Feel free to delete excess text after you have provided the relevant information.

**Software information**

- Software: are you using R or Python?

- If using Python:
```python
>>> import os; print(os.name, os.sys.platform)
>>> import sys; print(sys.version)
>>> import arcospy; print(arcospy.__version__)
```

- If using R:
```r
sessionInfo()
packageVersion("arcos")
```

**Describe the bug**
A clear and concise description of what the bug is.

**Code to Reproduce**
Please paste a reproducible example of the bug you are encountering.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

Thank you for taking the time to suggest a new feature for `arcos` and `arcospy`! Please see the guidelines outlined in the contributing.md document before filling out the following form.

**Is this a feature request related to the underlying ARCOS data or arcos/arcospy API?**
Please identify if your feature request seeks something in the ARCOS data or is a new function for the arcos/arcospy API.

**Would this feature help you obtain a measurement or data product that has precedent in the appropriate literature?**
If so, please link or identify the relevant literature. If not, please briefly explain why this measurement or data product would be useful to a wider audience.

**Is your feature request related to an existing problem? Please describe.**
A clear and concise description of what the problem is. 
Ex. "The current function does X, which is a problem because..."

**Describe the solution you'd like**
A clear and concise description of what you want to happen with the feature. This may inlcude:
- New field(s)
- Resulting data type(s)
- Geographic level(s) of the data (if applicable)
- Timespan of the data (if applicable)

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.
Ex. Have you tried creating the function locally using existing arcos/arcospy commands?

**Additional context**
Add any other context or screenshots about the feature request here.
Hi @edonnachie,

Thank you so much for the detailed comments! They were very helpful for improving the paper. I have updated the `paper.md` with several changes and provide responses to each of your comments below. 

---

*Author responses to reviewer comments:*

> Comment #1: It would be helpful to explain what the Opioid Overdose Crisis is, perhaps with link to the CDC. It's not clear to me whether the overdoses are accidental (due to overuse of this drug therapy) or suicide-related (due to the controls being too lax).

This is an excellent point. I have done a few things to address your comment. Firstly, I rewrote the opening of the introduction to contextualize the U.S. Opioid Overdose Epidemic. Secondly, I added the phrase "all-cause" when describing the increases in opioid-involved mortality as the Opioid Overdose Epidemic is concerned with both unintentional, intentional, and homicidal drug-involved deaths. The system-level *cause* (e.g. 'overuse of this drug therapy' or 'controls being too lax') remains an active area of research.  

> Comment #2: You write that data was released until 2012. You should mention here that the data for 2013-2014 is now also available.

Clarified that additional data for 2013-2014 is now also available. 

> Comment #3: You write that patient confidentiality is one reason that access to ARCOS is usually restricted. If this were true, then providing wide access to this data would be legally and ethically problematic (at least from a European perspective). However, I would tend to the view that the risk of revealing information about patients is negligible. Either way, this needs some explaining.

After consultation with co-authors I have reworked this sentence. Firstly, I replaced the word 'restricted' with 'only for approved requests (e.g. research or litigation)'. I have also added a citation for the Drug Enforcement Administration (DEA) Privacy Impact Assessment, which offers extremely detailed information on potential privacy concerns for the ARCOS data. I have also added the word 'anonymized' in the following sentence to clarify that no patient-level identifying information (e.g. name) exists in the present ARCOS data. The present ARCOS data represents *shipments* of controlled opioid medications from manufacturers to distributors, so you are correct in that  revealing information about patients is negligible. The original wording of the sentence was unclear and we hope that the new phrasing is helpful.

> Comment #4: Both from a technical and a data protection perspective, it would seem important to explain what is behind the API: How and where are the data stored? Who controls and maintains this database? Are requests logged, and if so, are personally identifiable data kept (e.g. IP address)?

I have added additional information on the API in the **API Structure** section including API specification and guidelines on use. I have reached out to members of the Data Reporting Team at *The Washington Post* for exact details on your comments regarding if and what personal identifiable information is stored.

> Comment #5: You should explain why the ARCOS database is useful. How does this compare with other sources of prescribing data (pharmacy data centres, medicaid etc)? I suspect that the main advantage lies in the analysis of how market dynamics facilitated the epidemic.

This is a great point. We feel that the first part of the comment ("How does this compare with other sources of prescribing data (pharmacy data centres, medicaid etc)?") is addressed by the sentences in the **Statement of Need** section. However, we have added additional sentences relating to the second point of the comment ("I suspect that the main advantage lies in the analysis of how market dynamics facilitated the epidemic.") as this is indeed a very important aspect of the Opioid Overdose Epidemic that has ample room left for exploration and quantification. 

> Comment #6: Repeated word: "In raw format, the the ARCOS database "

Deleted repeated word.

> Comment #7: Should "from 2006 to 2012" not be "from 2006 to 2014"?

Changed 2012 to 2014.

> Comment #8: The statement "while still maintaining individual medical privacy " is confusing, see above.

This is a great point. I have removed the quoted phrase as it is confusing and unnecessary.

> Comment #9: What you don't mention is that the data provide a different perspective to the medical claims data that would usually be used for such analyses. If I understand correctly, you're not looking at who received an opioid, but at how these substances traveled across the supply chain (and in what quantity). I think this is an important and interesting distinction to be made.

See author response to Comment #5 above. We agree and have added additional sentences relating to potential analyses of opioid prescription market dynamics.---
title: "Accessing ARCOS data in R and Python"
author: "Authors: arcos and arcospy developers"
date: "Last updated: August 7, 2020"
output:
  pdf_document: default
  html_document:
    code_folding: show
    theme: readable
    toc: yes
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(reticulate)
use_python("C:/Users/USER/AppData/Local/Programs/Python/Python38/python.exe")
```

## Installation

The original `arcos` software was written in R by members of the data reporting team (Steven Rich, Andrew Ba Tran, Aaron Williams, Jason Holt) at *The Washington Post*. Interested researchers at the University of Maryland (Jeffery Sauer, Dr. Taylor M. Oshan) translated the software to python and called it `arcospy`. Both software offer the same functionality as an API wrapper to query the ARCOS data (hosted [here](https://arcos-api.ext.nile.works/__swagger__/)). 

### R

The latest release of `arcos` is available via CRAN, while the development version is available from GitHub.

```{r, eval=FALSE, warning=FALSE}

# R

# CRAN

install.packages("arcos")

# Github

# install.packages("devtools")
devtools::install_github('wpinvestigative/arcos')

```

The R dependencies are `string`, `stringr`, `magrittr`, `jsonlite`, `dplyr`, `urltools`, and `vroom`.

### Python

Python is available on PyPI and can be installed via `pip`.

```{python, eval = FALSE, warning = FALSE}

# Python

pip install arcospy

```

The Python dependencies are `pandas` and `requests`. 

## Quick start

All commands share the name name between `arcos` and `arcospy`. Users can choose between commands that return raw data, commands that return summarized data, and commands that return potentially relevant supplemental data. Note that commands that return raw data may take significantly longer to run due to large size of the underlying datasets.

Once users have chosen a command that suits their needs, the general format of the API calls for the geographic areas of interest (such as county and state) as well as an API key. The default key is "WaPo" (additional keys available [here](https://github.com/wpinvestigative/arcos-api/blob/master/keys/keys.txt)). Outputs from all of the functions are delivered in popular formats, specifically `data.frames` in R and `pandas.DataFrame` in Python. We now demonstrate generic functionality with the `state_population` command to get annual population values for states between 2006 and 2012. 

### R API

```{r, warning = FALSE}

# R

library(arcos)
statepop <- state_population(state = "WV", key = "WaPo")
head(statepop)

```

Help for the R `arcos` package can be retrieved via `?arcos`. Additional documentation on the R API is available via [CRAN](https://cran.r-project.org/web/packages/arcos/arcos.pdf).

### Python API

```{python, warning = FALSE}

# Python

import arcospy
statepop = arcospy.state_population(state = "WV", key = "WaPo")
statepop.head()

```

Help for the Python `arcospy` package can be retrieved via `help(arcospy)`. Additional documentation on the Python API is available via [Github](https://github.com/jeffcsauer/arcospy).

## Example use: county-level analysis over time 

### Gathering data

One of the most straightforward analyses of the ARCOS data is examining trends in opioid orders at the county level over time. This analysis can be used to highlight relative increases and decreases of orders during the ramp-up phase of the Opioid Crisis (mid-2000s) towards its peak (early 2010s). To do so, we make use of the `summarized_county_annual` function to get the total number of opioid dosage units provided to each county of a given state from 2006 to 2014. For the purposes of this tutorial, we limit the analysis to Massachusetts. 

First, gather the total number of dosage units for each county of Massachusetts.

```{r, warning = FALSE, message = FALSE}

# R

MA_orders <- summarized_county_annual(state = "MA", key = "WaPo")

head(MA_orders)

```

```{python, warning = FALSE, message = FALSE}

# Python

MA_orders = arcospy.summarized_county_annual(state = "MA", key = "WaPo")

MA_orders.head()

```

### Adjusting by population

However, the above orders are unadjusted That is, they do not account for the underlying population wherein counties with more people are more likely to order more opioids. `arcos` has a built-in function to gather the population data called `county_population`. We now gather the population data, summarize the population data at the county level, join it to the opioid orders data, and create an average pills per person per year metric. 

```{r, warning = FALSE, message = FALSE}

# R

# Gather the population data

MA_pop <- county_population(state = "MA", key = "WaPo")

# Summarize the data at the county level

library(dplyr)
MA_pop <- MA_pop %>% 
  group_by(BUYER_COUNTY, BUYER_STATE, countyfips) %>% 
  summarize(avg_pop=mean(population, na.rm=T)) %>% 
  rename(buyer_county=BUYER_COUNTY, buyer_state=BUYER_STATE)

# Join to the opioid orders data

MA_merge <- left_join(MA_orders, MA_pop)

# Create an adjusted average pills per person per year metric

MA_merge <- MA_merge %>% 
  mutate(pills_per_person_avg=DOSAGE_UNIT/avg_pop/9)

head(MA_merge)
```

```{python, warning = FALSE, message = FALSE}

# Python

# Gather the population data

MA_pop = arcospy.county_population(state="MA", key="WaPo")

# Summarize the data at the county level

MA_pop = MA_pop.groupby(['BUYER_COUNTY', 'BUYER_STATE', 'countyfips']).mean().reset_index()

# Join to the opioid orders data

import pandas as pd
MA_merge = pd.merge(MA_orders, MA_pop, on=['BUYER_COUNTY', 'BUYER_STATE', 'countyfips'])

# Clean up duplicated year column 
del MA_merge['year_y']
MA_merge = MA_merge.rename(columns={"year_x": "year"})

# Create an adjusted average pills per person per year metric

MA_merge['pills_per_person'] = MA_merge['DOSAGE_UNIT']/MA_merge['population']/9

# Sort padnas.df to match R output from above
MA_merge = MA_merge.sort_values('BUYER_COUNTY')
MA_merge.head()

```

### Rendering graphs

Now that we have an adjusted measure of the number of pills ordered to a given county in a given year, we can create some simple visualizations to see county-level trends. For the purposes of this tutorial, we will render some straightforward line graphs, although additional documentation is available for [R users](https://wpinvestigative.github.io/arcos/articles/annual-maps.html) and [Python users](https://github.com/jeffcsauer/arcospy/blob/master/docs/Intro%20to%20arcospy%20-%20downloading%20data%20and%20making%20it%20spatial.ipynb) on how to make maps using the ARCOS data.

```{r, warning=FALSE, message=FALSE}

# R

library(ggplot2)

# Render plot

MA_merge %>%
  ggplot(aes(x=year, 
             y=pills_per_person_avg,
             group=BUYER_COUNTY,
             color=BUYER_COUNTY)) +
    geom_line() +
  ylab('Average pills per person')+
  xlab('Year') +
  ggtitle('County-level trends in opioid pills in Massachusets from 2006 to 2014. \nData source: ARCOS/WP.')


```

```{python, warning=FALSE, message=FALSE, error=FALSE}

# Python

import seaborn as sns
sns.set(font_scale=0.8)
import matplotlib.pyplot as plt

# Render plot
# Note: we need to include a few extra lines to render the 
# legend outside the plot. See: 
# https://stackoverflow.com/questions/30490740/move-legend-outside-figure-in-seaborn-tsplot

fig, ax1 = plt.subplots(1,1)

g = sns.lineplot(x="year", y="pills_per_person", hue="BUYER_COUNTY", data=MA_merge, ax=ax1)
g.set(xlabel="Year", ylabel="Average pills per person")
g.set_title("County-level trends in opioid pills in Massachusets from 2006 to 2014. \nData source: ARCOS/WP.")

# Automatically adjust bounding box of figure

box = g.get_position()
g.set_position([box.x0, box.y0, box.width * 0.85, box.height]) 

# Place legend

g.legend(loc='center right', bbox_to_anchor=(1.30, 0.5), ncol=1)

plt.show()

```

The above spaghetti graphs help to identify the range of opioid availability in a given county in Massachusetts. While the numbers may seem relatively low (e.g. only 1-5 pills), remember that this is a **population-level** measure. There is also clear variation among counties in a given state, even after adjusting for population. 

Users should note that the above visualization methods may not be suitable for all states. Massachusetts has only 10 counties, whereas most states have more than 20 counties. So many lines would likely render the above plot uninterpretable, so we encourage users to pursue additional options. 

## Important considerations for use

As preivously hinted at with the adjustment of population, there are a few important considerations users should make when using the ARCOS data. 

First and foremost, adjusting by population is almost always a safe adjustment to be made across raw, summarized, and pharmacy-level data. While we provide one method of adjustment above, users may also adjust by yearly-specific measures of population as well.

Secondly, users might consider looking at specific subsets of drug distributors. Specifically, users should examine the variable `BUYER_BUS_ACT` returned by several functions. This variable describes the type of distributor (chain, retail pharmacy, private, medical, etc). Depending on the type of interest of the user, certain types of distributors may or may not need to be excluded. 

Lastly, as discussed elsewhere in the [`arcos` and `arcospy` documentation](https://github.com/jeffcsauer/arcos_arcospy_information/blob/master/paper.md) there are numerous ways to conceptualize the unit of analysis for the distribution of opioids. In the data currently available via `arcos` and `arcospy`, these units are primarily:

- total number of records for a geographic unit (adjusted or unadjusted)
- total number of all opioid pills (adjusted or unadjusted)
- total number of specific opioid pills (adjusted or unadjusted)
- total amount (in weight) of all or specific opioid pills (adjusted or unadjusted)

However, there are other common units of analysis that appear in the literature, such as morphine milligram equivalents (MMEs) or prescription counts, although these are not directly observable in the present data. It is the responsibility of the user to select an appropriate unit of analysis that has precedent in their field or take appropriate steps to standardize the data.

## Additional resources

- [The Opioid Files: a series of investigative articles by *The Washington Post*](https://www.washingtonpost.com/national/2019/07/20/opioid-files/?arc404=true)

- [Drug Enforcement Agency presentation on the structure and purpose of ARCOS](https://www.deadiversion.usdoj.gov/mtgs/distributor/conf_2016/weisman.pdf#search=ARCOS%20Registrant%20Handbook)

- [Example academic publication using earlier versions of ARCOS data](https://www.jpsmjournal.com/article/S0885-3924(04)00182-4/fulltext)
