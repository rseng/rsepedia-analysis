# awardFindR
<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/awardFindR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/awardFindR/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/awardFindR/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/awardFindR?branch=master)
[![](https://badges.ropensci.org/432_status.svg)](https://github.com/ropensci/software-review/issues/432)
<!-- badges: end -->

`awardFindR` is a framework that scrapes and searches a variety of grant and award databases for specific keywords. These include major US federal agencies like NSF and NIH as well as private organizations like the Bill & Melinda Gates Foundation. Results from searching each of these databases are collected and made available to users. The package is designed to be modular and extensible, supporting any number of APIs and other web-based sources that provide award data available online.

Packages that have provided similar functionality include [awardsBot](https://github.com/NCEAS/awards-bot).

## Installation
Install `awardFindR` directly from this repository using the `remotes` package
```
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("ropensci/awardFindR")
```

## Supported sources

- [Arnold Ventures](https://www.arnoldventures.org/grants-search) (`arnold`)
- [Carnegie Corporation of New York](https://www.carnegie.org/grants/grants-database/) (`carnegie`)
- [Federal Reporter](https://federalreporter.nih.gov/) (`fedreporter`)
- [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/) (`gates`)
- [MacArthur Foundation](https://www.macfound.org/grants/) (`macarthur`)
- [Mellon Foundation](https://mellon.org/grants/grants-database/advanced-search/) (`mellon`)
- [National Endowment for the Humanities](https://securegrants.neh.gov/publicquery/main.aspx) (`neh`)
- [NIH RePORTER](https://projectreporter.nih.gov/reporter.cfm) (`nih`)
- [National Science Foundation](https://nsf.gov/awardsearch/) (`nsf`)
- [Open Philanthropy](https://www.openphilanthropy.org/giving/grants) (`ophil`)
- [Open Society Foundations](https://www.opensocietyfoundations.org/grants/past) (`osociety`)
- [Rockefeller Foundation](https://www.rockefellerfoundation.org/) (`rockefeller`)
- [Russell Sage Foundation](https://www.russellsage.org) (`rsf`)
- [Robert Wood Johnson Foundation](https://www.rwjf.org/en/how-we-work/grants-explorer.html) (`rwjf`)
- [Alfred P. Sloan Foundation](https://sloan.org/grants-database) (`sloan`)
- [Social Science Research Council](https://www.ssrc.org) (`ssrc`)
- [John Templeton Foundation](https://www.templeton.org/grants/grant-database) (`templeton`)
- [USASpending.gov](https://www.usaspending.gov/search) (`usaspend`)

## How the package works

`search_awards` has parameters to change keywords, sources and dates as search criteria, which are passed on to source routines. Dates are interpreted with varying degrees of precision based on the data available from each source. See included help on individual sources to understand their respective limitations.

Individual sources have functions of their own. `get_nsf` specifically fetches results from NSF, and `get_nih` fetches results from NIH, for example. These are meant to provide end users with higher levels of detail if they're interested in a specific source.

### Quick introduction

Search for all awards matching a keyword since the default cutoff of Jan 1, 2019 to today

```
awardFindR::search_awards(keywords="ethnography")
```

See the included vignette for additional examples.

For those interested in the results from a specific source, each source has its own function. For example, someone interested in NSF results for "ethnography" between 2018 and 2020 could run the following:

```
awardFindR::get_nsf("ethnography", "2018-01-01", "2020-01-01")
```

Similar functions exist for each supported source. See included help for further details, as the arguments differ slightly between each.

### Dependencies

This package depends on `rvest`, `xml2` and `httr`.

### Funding

This package was developed with support from [a grant from the Sloan Foundation](https://sloan.org/grant-detail/8725).
# Version 1.0.1

- Instead of extracting the 2020s csv from NEH, we use the actual API

# awardFindR (1.0)

- Transferred repo to ropensci
- Added NEWS.md file
- Changed status from "WIP" to "Active"
- Added `verbose` option to all functions that pull HTTP resources, so users can see the status of individual HTTP requests
- Changed all function names to a standardized verb_object format
- `sources` argument in `search_awards()` is now validated properly with `match.args()`
- Various bug fixes
# Contributing
`awardFindR` is built to be easily extendable to include additional data sources and we welcome contributions adding support for additional databases of research funding. 

Adding a new source involves three steps:
1. Adding a `sourcename.R` file that parse the database using its API or webscraping as discussed
2. Adding the source to default `sources` argument of the `search_awards()` function in `main.R`
3. Adding relevant tests

## Adding a new source 

New sources are entirely self-contained in a single `.R` file.

Sources need to have a minimum of two functions:

1.  a `get_*` function that scrapes raw data from the source and returns it as-is,
2.  and a hidden `.standardize_*` (note the preceding period) function that harmonizes the data to fit the data.frame output from `search_awards()`.

The file containing the source should have an easily identifiable filename, which should be reflected in the naming of the `get_*` and `standardize_*` functions. For example, import from the "Megafunds Foundation" should be in `megafunds.R`, containing `get_megafunds()` and `.standardize_megafunds()`.

### `get_*`

The `get_*` routine is exported, so end users can call it directly if they have interest in a specific data source. It should be as faithful a reproduction of the original data source as possible, though variables should still be R-friendly (for example, award amounts as strings like "\$300,000" aren't useful. An integer of 300000 is much more useful.)

These routines will usually use HTML or JSON-based HTTP resources. `utils.R` in this package has a `request()` function that handles HTTP POST and GET requests with `httr` and returns xml2 objects (which can be used with `rvest` or `xml2`) or json as a base R `list()` object automatically. In order to centralize HTTP error handling and standardize message output to the user, please employ this function as much as possible for HTTP requests.

### Style and dependencies

This function can be largely tailored to the needs of an API, but should follow some basic style guidelines. Three arguments are required:

1. Keyword
2. From date
3. To date

The function should accept at least one keyword and two date terms as arguments. _keyword_ is a string, or, if the source can handle more than one keyword at a time, _keywords_ is a vector of strings. Dates are typically exact days or years, depending on the capabilities of the search function available from the source. Exact dates should be date objects named _from_date_ and _to_date_, while years should be integers named _from_year_ and _to_year_. 

Additional source-specific variables can be added to the function, but should have default values specified. `*_get` routines should function as expected with these three variables alone.

If at all possible, searching should be done server-side to reduce HTTP traffic. Downloading the whole grant database and searching with `grep()` or a similar function should be a last resort. This also minimizes CPU load.

JSON sources shouldn't need any additional dependencies, since the `request()` function handles encoding and decoding from/to lists automatically. HTML/XML sources, however, should use `rvest` or `xml2`, respectively. 

If the source does not use HTML, XML or JSON, it may be necessary to employ additional dependencies as necessary. This is an extreme case and should be done with care.

### `.standardize_*`

The `awardFindR` calls the  internal`.*_standardize` function, which should in turn call the `get_*` function described above. All `.standardize_*` functions need to accept the exact same input. This includes:

1. *keywords* - A vector of keywords to search 
2. *from_date* - A beginning date object (e.g. "2015-01-01")
3. *end_date* - And ending date object

The date objects should be translated into whatever the source-specific requirements are for search terms. If an API can only delineate searches by year, get the year. If an API can only handle one keyword at a time, loop the function through multiple keywords with `lapply()`.

This function needs to return a `data.frame` with the following columns:

1.  *institution* - Institution name
2.  *pi* -  Principal investigator name
3.  *year* - Award year
4.  *start* - Award start date
5.  *end* - Award end date
6.  *program* - Program or department name
7.  *amount* - Award amount (in US dollars)
8.  *id* - Award ID number (or any way of uniquely identifying duplicates from one source)
9.  *title* - Project title or description
10. *keyword* - The searched keyword that returned this result
11. *source* - The source the result came from

## Including the Source in `search_awards()`

Source routines are included in the `search_awards()` function through the _sources_ argument. Include the name of your new source (the section before "_standardize") in the default value of the _sources_ for `search_awards()` in `main.R`. This will include your new source routine in the default functionality, and expose the name in the documentation where users are likely to discover it.

## Adding Tests

The main tests in `test-full.R` attempt to exercise both the successful and no results code branches. The latter always works. Unfortunately, the default search for the former (the terms "qualitative analysis" and "ethnography" in 2018) does not actually return results for some smaller sources. For this and other reasons, it's a good idea to create a source-specific test for the `.standardize_*` function.

Tests should actually return results, but should also minimally tax the API. Less than 10 results would be ideal. There is one exception to this: when a source needs to loop through multiple pages. In this case, the smallest number of results that triggers the looping is ideal, in order to ensure maximum test coverage.

One test should verify reproducibility; i.e. we get the same results again for a query in a specified date range. When used with a real-world HTTP resource, this should alert us if there is a change in how the resource is provided.

Finally, HTTP requests for tests should be cached using the ROpenSci package [`vcr`](https://docs.ropensci.org/vcr/). Individual sources should have their own cassettes. This limits stress on the APIs from frequent testing of internal logic, especially by the continuous integration system. Failure to do this could potentially result in service-level bans if we're not careful.
