ghql
====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/ghql)](https://cranchecks.info/pkgs/ghql)
[![R-check](https://github.com/ropensci/ghql/workflows/R-check/badge.svg)](https://github.com/ropensci/ghql/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/ghql/coverage.svg?branch=master)](https://codecov.io/github/ropensci/ghql?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ghql)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/ghql)](https://cran.r-project.org/package=ghql)

`ghql` - a GraphQL client for R

GraphQL - <https://graphql.org>

Examples of GraphQL APIs:

* GitHub: <https://docs.github.com/en/graphql/guides/introduction-to-graphql>
* Opentargets: <https://genetics-docs.opentargets.org/technical-pipeline/graphql-api>
* Countries GraphQL API: <https://github.com/trevorblades/countries>

Other GraphQL R packages:

* [graphql][] - GraphQL query parser
* [gqlr][] - GraphQL server and query methods

## Install

CRAN version


```r
install.packages("ghql")
```

Development version


```r
remotes::install_github("ropensci/ghql")
```


```r
library("ghql")
library("jsonlite")
library("dplyr")
```

## Package Documentation

<https://docs.ropensci.org/ghql/>

## Meta

* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[gqlr]: https://github.com/schloerke/gqlr
[graphql]: https://github.com/ropensci/graphql
[libgraphqlparser]: https://github.com/graphql/libgraphqlparser
ghql 0.1.0
==========

### NEW FEATURES

* Released to CRAN
## Test environments

* local OS X install, R 3.6.2 patched
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

New submission

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.

---

This is a new release. I have read and agree to the the 
CRAN policies at https://cran.r-project.org/web/packages/policies.html

Thanks!
Scott Chamberlain
<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/ghql/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/ghql.git`
* Make sure to track progress upstream (i.e., on our version of `ghql` at `ropensci/ghql`) by doing `git remote add upstream https://github.com/ropensci/ghql.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/ghql`

### Check out our [discussion forum](https://discuss.ropensci.org)
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
ghql
====

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/ghql)](https://cranchecks.info/pkgs/ghql)
[![R-check](https://github.com/ropensci/ghql/workflows/R-check/badge.svg)](https://github.com/ropensci/ghql/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/ghql/coverage.svg?branch=master)](https://codecov.io/github/ropensci/ghql?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ghql)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/ghql)](https://cran.r-project.org/package=ghql)

`ghql` - a GraphQL client for R

GraphQL - <https://graphql.org>

Examples of GraphQL APIs:

* GitHub: <https://docs.github.com/en/graphql/guides/introduction-to-graphql>
* Opentargets: <https://genetics-docs.opentargets.org/technical-pipeline/graphql-api>
* Countries GraphQL API: <https://github.com/trevorblades/countries>

Other GraphQL R packages:

* [graphql][] - GraphQL query parser
* [gqlr][] - GraphQL server and query methods

## Install

CRAN version

```{r eval=FALSE}
install.packages("ghql")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/ghql")
```

```{r}
library("ghql")
library("jsonlite")
library("dplyr")
```

## Package Documentation

<https://docs.ropensci.org/ghql/>

## Meta

* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[gqlr]: https://github.com/schloerke/gqlr
[graphql]: https://github.com/ropensci/graphql
[libgraphqlparser]: https://github.com/graphql/libgraphqlparser
---
title: Introduction to ghql
author: Scott Chamberlain
date: "2021-01-25"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



`ghql` - a GraphQL client for R

## GitHub Authentication

Note: To be clear, this R package isn't just for the GitHub GraphQL API, but it
is the most public GraphQL API we can think of, so is used in examples
throughout here.

See <https://docs.github.com/en/graphql/guides/forming-calls-with-graphql#authenticating-with-graphql> for getting an OAuth token.

Store the token in a env var called `GITHUB_GRAPHQL_TOKEN`

## Install

CRAN version


```r
install.packages("ghql")
```

Development version


```r
remotes::install_github("ropensci/ghql")
```


```r
library("ghql")
library("jsonlite")
library("dplyr")
```

## initialize client


```r
token <- Sys.getenv("GITHUB_GRAPHQL_TOKEN")
con <- GraphqlClient$new(
  url = "https://api.github.com/graphql",
  headers = list(Authorization = paste0("Bearer ", token))
)
```

## load schema

Since not every GraphQL server has a schema at the base URL, have to manually
load the schema in this case


```r
con$load_schema()
```


## Queries

Make a `Query` class object


```r
qry <- Query$new()
```

When you construct queries we check that they are properly formatted using the 
[graphql][] package that leverages the [libgraphqlparser][] C++ parser. If the query
is malformed, we return a message as to why the query is malformed.

Get some stargazer counts


```r
qry$query('mydata', '{
  repositoryOwner(login:"sckott") {
    repositories(first: 5, orderBy: {field:PUSHED_AT,direction:DESC}, isFork:false) {
      edges {
        node {
          name
          stargazers {
            totalCount
          }
        }
      }
    }
  }
}')
qry
#> <ghql: query>
#>   queries:
#>     mydata
qry$queries$mydata
#>  
#>  {
#>   repositoryOwner(login:"sckott") {
#>     repositories(first: 5, orderBy: {field:PUSHED_AT,direction:DESC}, isFork:false) {
#>       edges {
#>         node {
#>           name
#>           stargazers {
#>             totalCount
#>           }
#>         }
#>       }
#>     }
#>   }
#> }
```


```r
# returns json
(x <- con$exec(qry$queries$mydata))
#> [1] "{\"data\":{\"repositoryOwner\":{\"repositories\":{\"edges\":[{\"node\":{\"name\":\"cranchecksdocs\",\"stargazers\":{\"totalCount\":6}}},{\"node\":{\"name\":\"badges\",\"stargazers\":{\"totalCount\":0}}},{\"node\":{\"name\":\"mutant-proposal\",\"stargazers\":{\"totalCount\":0}}},{\"node\":{\"name\":\"extcite\",\"stargazers\":{\"totalCount\":6}}},{\"node\":{\"name\":\"Headstart\",\"stargazers\":{\"totalCount\":140}}}]}}}}\n"
# parse to an R list
jsonlite::fromJSON(x)
#> $data
#> $data$repositoryOwner
#> $data$repositoryOwner$repositories
#> $data$repositoryOwner$repositories$edges
#>         node.name node.totalCount
#> 1  cranchecksdocs               6
#> 2          badges               0
#> 3 mutant-proposal               0
#> 4         extcite               6
#> 5       Headstart             140
```

## Parameterize a query by a variable

Define a query


```r
qry <- Query$new()
qry$query('getgeninfo', 'query getGeneInfo($genId: String!){
  geneInfo(geneId: $genId) {
    id
    symbol
    chromosome
    start
    end
    bioType
    __typename
  }
}')
```

Define a variable as a named list


```r
variables <- list(genId = 'ENSG00000137033')
```

Creat a clint and make a request, passing in the query and then the variables


```r
con <- GraphqlClient$new('https://genetics-api.opentargets.io/graphql')
res <- con$exec(qry$queries$getgeninfo, variables)
jsonlite::fromJSON(res)
#> $data
#> $data$geneInfo
#> $data$geneInfo$id
#> [1] "ENSG00000137033"
#> 
#> $data$geneInfo$symbol
#> [1] "IL33"
#> 
#> $data$geneInfo$chromosome
#> [1] "9"
#> 
#> $data$geneInfo$start
#> [1] 6215786
#> 
#> $data$geneInfo$end
#> [1] 6257983
#> 
#> $data$geneInfo$bioType
#> [1] "protein_coding"
#> 
#> $data$geneInfo$`__typename`
#> [1] "Gene"
```

## Example: Datacite

[Datacite](https://datacite.org/) provides DOIs for research data. Check out the 
[Datacite GraphQL docs](https://support.datacite.org/docs/datacite-graphql-api-guide)
to get started. A minimal example:


```r
con <- GraphqlClient$new("https://api.datacite.org/graphql")
qry <- Query$new()
qry$query('dc', '{
  publications(query: "climate") {
    totalCount

    nodes {
      id
      titles {
        title
      }
      descriptions {
        description
      }
      creators {
        name
        familyName
      }
      fundingReferences {
        funderIdentifier
        funderName
        awardTitle
        awardNumber
      }
    }
  }
}')
res <- con$exec(qry$queries$dc)
head(jsonlite::fromJSON(res)$data$publications$nodes)
#>                                     id                                titles
#> 1 https://doi.org/10.4122/1.1000000046 Single Cell Protein from Landfill Gas
#> 2 https://doi.org/10.4122/1.1000000046 Single Cell Protein from Landfill Gas
#> 3 https://doi.org/10.4122/1.1000000047 Single Cell Protein from Landfill Gas
#> 4 https://doi.org/10.4122/1.1000000047 Single Cell Protein from Landfill Gas
#> 5 https://doi.org/10.4122/1.1000000048 Single Cell Protein from Landfill Gas
#> 6 https://doi.org/10.4122/1.1000000048 Single Cell Protein from Landfill Gas
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       descriptions
#> 1 Municipal solid waste (MSW) landfills are one of the largest human-generated sources of methane emissions in the United States and other countries globally. Methane is believed to be a very potent greenhouse gas that is a key contributor to global climate change, over 21 times stronger than CO2. Methane also has a short (10-year) atmospheric life. Because methane is both potent and short-lived, reducing methane emissions from MSW landfills is one of the best ways to achieve a near-term beneficial impact in mitigating global climate change. The United States Environmental Protection Agency estimates that a landfill gas (LFG) project will capture roughly 60-90% of the methane emitted from the landfill, depending on system design and effectiveness. The captured methane can be then purified and used for industrial applications, as in this case the production of SCP. Utilizing methane in this way decreases its demand from fossil fuels which is its traditional source.
#> 2 Municipal solid waste (MSW) landfills are one of the largest human-generated sources of methane emissions in the United States and other countries globally. Methane is believed to be a very potent greenhouse gas that is a key contributor to global climate change, over 21 times stronger than CO2. Methane also has a short (10-year) atmospheric life. Because methane is both potent and short-lived, reducing methane emissions from MSW landfills is one of the best ways to achieve a near-term beneficial impact in mitigating global climate change. The United States Environmental Protection Agency estimates that a landfill gas (LFG) project will capture roughly 60-90% of the methane emitted from the landfill, depending on system design and effectiveness. The captured methane can be then purified and used for industrial applications, as in this case the production of SCP. Utilizing methane in this way decreases its demand from fossil fuels which is its traditional source.
#> 3 Municipal solid waste (MSW) landfills are one of the largest human-generated sources of methane emissions in the United States and other countries globally. Methane is believed to be a very potent greenhouse gas that is a key contributor to global climate change, over 21 times stronger than CO2. Methane also has a short (10-year) atmospheric life. Because methane is both potent and short-lived, reducing methane emissions from MSW landfills is one of the best ways to achieve a near-term beneficial impact in mitigating global climate change. The United States Environmental Protection Agency estimates that a landfill gas (LFG) project will capture roughly 60-90% of the methane emitted from the landfill, depending on system design and effectiveness. The captured methane can be then purified and used for industrial applications, as in this case the production of SCP. Utilizing methane in this way decreases its demand from fossil fuels which is its traditional source.
#> 4 Municipal solid waste (MSW) landfills are one of the largest human-generated sources of methane emissions in the United States and other countries globally. Methane is believed to be a very potent greenhouse gas that is a key contributor to global climate change, over 21 times stronger than CO2. Methane also has a short (10-year) atmospheric life. Because methane is both potent and short-lived, reducing methane emissions from MSW landfills is one of the best ways to achieve a near-term beneficial impact in mitigating global climate change. The United States Environmental Protection Agency estimates that a landfill gas (LFG) project will capture roughly 60-90% of the methane emitted from the landfill, depending on system design and effectiveness. The captured methane can be then purified and used for industrial applications, as in this case the production of SCP. Utilizing methane in this way decreases its demand from fossil fuels which is its traditional source.
#> 5 Municipal solid waste (MSW) landfills are one of the largest human-generated sources of methane emissions in the United States and other countries globally. Methane is believed to be a very potent greenhouse gas that is a key contributor to global climate change, over 21 times stronger than CO2. Methane also has a short (10-year) atmospheric life. Because methane is both potent and short-lived, reducing methane emissions from MSW landfills is one of the best ways to achieve a near-term beneficial impact in mitigating global climate change. The United States Environmental Protection Agency estimates that a landfill gas (LFG) project will capture roughly 60-90% of the methane emitted from the landfill, depending on system design and effectiveness. The captured methane can be then purified and used for industrial applications, as in this case the production of SCP. Utilizing methane in this way decreases its demand from fossil fuels which is its traditional source.
#> 6 Municipal solid waste (MSW) landfills are one of the largest human-generated sources of methane emissions in the United States and other countries globally. Methane is believed to be a very potent greenhouse gas that is a key contributor to global climate change, over 21 times stronger than CO2. Methane also has a short (10-year) atmospheric life. Because methane is both potent and short-lived, reducing methane emissions from MSW landfills is one of the best ways to achieve a near-term beneficial impact in mitigating global climate change. The United States Environmental Protection Agency estimates that a landfill gas (LFG) project will capture roughly 60-90% of the methane emitted from the landfill, depending on system design and effectiveness. The captured methane can be then purified and used for industrial applications, as in this case the production of SCP. Utilizing methane in this way decreases its demand from fossil fuels which is its traditional source.
#>                                                            creators
#> 1 Babi, Deenesh, Price, Jason, Woodley, Prof. John, Babi, Price, NA
#> 2 Babi, Deenesh, Price, Jason, Woodley, Prof. John, Babi, Price, NA
#> 3 Babi, Deenesh, Price, Jason, Woodley, Prof. John, Babi, Price, NA
#> 4 Babi, Deenesh, Price, Jason, Woodley, Prof. John, Babi, Price, NA
#> 5 Babi, Deenesh, Price, Jason, Woodley, Prof. John, Babi, Price, NA
#> 6 Babi, Deenesh, Price, Jason, Woodley, Prof. John, Babi, Price, NA
#>   fundingReferences
#> 1              NULL
#> 2              NULL
#> 3              NULL
#> 4              NULL
#> 5              NULL
#> 6              NULL
```

## Example: Countries Data
A public GraphQL API for information about countries, continents, and languages. This project uses Countries List and provinces as data sources, so the schema follows the shape of that data, with a few exceptions:

Link to the GraphQL schema api

```r
link <- 'https://countries.trevorblades.com/'
```

Create a new graphqlClient object 

```r
con <- GraphqlClient$new(url = link)
```

Define a Graphql Query

```r
query <- '
query($code: ID!){
  country(code: $code){
    name
    native
    capital
    currency
    phone
    languages{
      code
      name
    }
  }
}'
```

The `ghql` query class and define query in a character string

```r
new <- Query$new()$query('link', query)
```

Inspecting the schema

```r
new$link
#>  
#>  
#> query($code: ID!){
#>   country(code: $code){
#>     name
#>     native
#>     capital
#>     currency
#>     phone
#>     languages{
#>       code
#>       name
#>     }
#>   }
#> }
```

define a variable as a named list

```r
variable <- list(
  code = "DE"
)
```

Making a request, passing in the query and then the variables. Then you convert the raw object to a structured json object

```r
result <- con$exec(new$link, variables = variable) %>% 
  fromJSON(flatten = FALSE)
result
#> $data
#> $data$country
#> $data$country$name
#> [1] "Germany"
#> 
#> $data$country$native
#> [1] "Deutschland"
#> 
#> $data$country$capital
#> [1] "Berlin"
#> 
#> $data$country$currency
#> [1] "EUR"
#> 
#> $data$country$phone
#> [1] "49"
#> 
#> $data$country$languages
#>   code   name
#> 1   de German
```

Convert the json data into a tibble object

```r
country_data <- result$data$country %>% 
  as_tibble()
country_data
#> # A tibble: 1 x 6
#>   name    native      capital currency phone languages$code $name 
#>   <chr>   <chr>       <chr>   <chr>    <chr> <chr>          <chr> 
#> 1 Germany Deutschland Berlin  EUR      49    de             German
```

## run a local GraphQL server

* Copy the `server.js` file from this package located at `inst/server.js` somewhere on your machine. Can locate it on your machine like `system.file("js/server.js", package = "ghql")`. Or you can run the file from where it's at, up to you.
* Make sure node is installed. If not, see <https://nodejs.org>
* Run `node server.js`
* Navigate to your browser - go to http://localhost:4000/graphql
* Back in R, user that URL to connect


```r
(con <- GraphqlClient$new("http://localhost:4000/graphql"))
#> <ghql client>
#>   url: http://localhost:4000/graphql
```


```r
xxx <- Query$new()
xxx$query('query', '{
  __schema {
    queryType {
      name, 
      fields {
        name,
        description
      }
    }
  }
}')
```



```r
con$exec(xxx$queries$query)
#> $data
#> $data$`__schema`
#> $data$`__schema`$queryType
#> $data$`__schema`$queryType$name
#> [1] "Query"
#> 
#> $data$`__schema`$queryType$fields
#>    name description
#> 1 hello            
#> 2  name 
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragment-class.R
\name{Fragment}
\alias{Fragment}
\title{Fragment}
\value{
a `Fragment` class (R6 class)
}
\description{
ghql fragment class
}
\examples{
# make a fragment class
frag <- Fragment$new()

# define a fragment
frag$fragment('Watchers', '
  fragment on Repository {
    watchers(first: 3) {
      edges {
        node {
          name
       }
    }
  }
}')

# define another fragment
frag$fragment('Stargazers', '
  fragment on Repository {
    stargazers(first: 3) {
      edges {
        node {
          name
       }
    }
  }
}')
frag
frag$fragments
frag$fragments$Watchers
frag$fragments$Stargazers
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{fragments}}{(list) list of fragments}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{Fragment$print()}}
\item \href{#method-fragment}{\code{Fragment$fragment()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the `Fragment` class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Fragment$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragment"></a>}}
\if{latex}{\out{\hypertarget{method-fragment}{}}}
\subsection{Method \code{fragment()}}{
create a fragment by name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Fragment$fragment(name, x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{(character) fragment name}

\item{\code{x}}{(character) the fragment}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets fragments internally
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client.R
\name{GraphqlClient}
\alias{GraphqlClient}
\title{GraphqlClient}
\value{
a `GraphqlClient` class (R6 class)
}
\description{
R6 class for constructing GraphQL queries
}
\examples{
x <- GraphqlClient$new()
x

\dontrun{
# make a client
token <- Sys.getenv("GITHUB_GRAPHQL_TOKEN")
cli <- GraphqlClient$new(
  url = "https://api.github.com/graphql",
  headers = list(Authorization = paste0("Bearer ", token))
)

# if the GraphQL server has a schema, you can load it
cli$load_schema()

# dump schema to local file
f <- tempfile(fileext = ".json")
cli$dump_schema(file = f)
readLines(f)
jsonlite::fromJSON(readLines(f))

# after dumping to file, you can later read schema from file for faster loading
rm(cli)
cli <- GraphqlClient$new(
  url = "https://api.github.com/graphql",
  headers = list(Authorization = paste0("Bearer ", token))
)
cli$load_schema(schema_file = f)

# variables
cli$url
cli$schema
cli$schema$data
cli$schema$data$`__schema`
cli$schema$data$`__schema`$queryType
cli$schema$data$`__schema`$mutationType
cli$schema$data$`__schema`$subscriptionType
head(cli$schema$data$`__schema`$types)
cli$schema$data$`__schema`$directives


# methods
## ping - hopefully you get TRUE
cli$ping()

## dump schema
cli$schema2json()


## define query
### creat a query class first
qry <- Query$new()
## another
qry$query('repos', '{
  viewer {
    repositories(last: 10, isFork: false, privacy: PUBLIC) {
      edges {
        node {
          isPrivate
          id
          name
        }
      }
    }
  }
}')
qry
qry$queries
qry$queries$repos
### execute the query
cli$exec(qry$queries$repos)


# query with a fragment
### define query without fragment, but referring to it
qry <- Query$new()
qry$query('queryfrag', '{
  ropensci: repositoryOwner(login:"ropensci") {
    repositories(first: 3) {
      edges {
        node {
          ...Watchers
        }
      }
    }
  }
  ropenscilabs: repositoryOwner(login:"ropenscilabs") {
    repositories(first: 3) {
      edges {
        node {
          ...Watchers
        }
      }
    }
  }
}')

### define a fragment
frag <- Fragment$new()
frag$fragment('Watchers', '
  fragment on Repository {
    watchers(first: 3) {
      edges {
        node {
          name
       }
    }
  }
}')
frag$fragments
frag$fragments$Watchers

### add the fragment to the query 'queryfrag'
qry$add_fragment('queryfrag', frag$fragments$Watchers)
qry
qry$queries$queryfrag

### execute query: we'll hook together the query and your fragment internally
cli$exec(qry$queries$queryfrag)
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{url}}{(character) list of fragments}

\item{\code{headers}}{list of named headers}

\item{\code{schema}}{holds schema}

\item{\code{result}}{holds result from http request}

\item{\code{fragments}}{(list) list of fragments}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{GraphqlClient$new()}}
\item \href{#method-print}{\code{GraphqlClient$print()}}
\item \href{#method-ping}{\code{GraphqlClient$ping()}}
\item \href{#method-load_schema}{\code{GraphqlClient$load_schema()}}
\item \href{#method-dump_schema}{\code{GraphqlClient$dump_schema()}}
\item \href{#method-schema2json}{\code{GraphqlClient$schema2json()}}
\item \href{#method-fragment}{\code{GraphqlClient$fragment()}}
\item \href{#method-exec}{\code{GraphqlClient$exec()}}
\item \href{#method-prep_query}{\code{GraphqlClient$prep_query()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new `GraphqlClient` object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$new(url, headers)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{url}}{(character) URL for the GraphQL schema}

\item{\code{headers}}{Any acceptable headers, a named list. See examples}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new `GraphqlClient` object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the `GraphqlClient` class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ping"></a>}}
\if{latex}{\out{\hypertarget{method-ping}{}}}
\subsection{Method \code{ping()}}{
ping the GraphQL server
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$ping(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{curl options passed on to [crul::verb-HEAD]}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
`TRUE` if successful response, `FALSE` otherwise
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-load_schema"></a>}}
\if{latex}{\out{\hypertarget{method-load_schema}{}}}
\subsection{Method \code{load_schema()}}{
load schema, from URL or local file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$load_schema(schema_url = NULL, schema_file = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{schema_url}}{(character) url for a schema file}

\item{\code{schema_file}}{(character) path to a schema file}

\item{\code{...}}{curl options passed on to [crul::verb-GET]}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing, loads schema into `$schema` slot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dump_schema"></a>}}
\if{latex}{\out{\hypertarget{method-dump_schema}{}}}
\subsection{Method \code{dump_schema()}}{
dump schema to a local file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$dump_schema(file)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file}}{(character) path to a file}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing, writes schema to `file`
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-schema2json"></a>}}
\if{latex}{\out{\hypertarget{method-schema2json}{}}}
\subsection{Method \code{schema2json()}}{
convert schema to JSON
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$schema2json(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{options passed on to [jsonlite::toJSON()]}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
json
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragment"></a>}}
\if{latex}{\out{\hypertarget{method-fragment}{}}}
\subsection{Method \code{fragment()}}{
load schema, from URL or local file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$fragment(name, x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{(character) fragment name}

\item{\code{x}}{(character) the fragment}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets fragments internally
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exec"></a>}}
\if{latex}{\out{\hypertarget{method-exec}{}}}
\subsection{Method \code{exec()}}{
execute the query
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$exec(query, variables, encoding = "UTF-8", ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{query}}{(character) a query, of class `query` or `fragment`}

\item{\code{variables}}{(list) named list with query variables values}

\item{\code{encoding}}{(character) encoding to use to parse the response. passed
on to [crul::HttpResponse] `$parse()` method. default: "UTF-8"}

\item{\code{...}}{curl options passed on to [crul::verb-POST]}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character string of response, if successful
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-prep_query"></a>}}
\if{latex}{\out{\hypertarget{method-prep_query}{}}}
\subsection{Method \code{prep_query()}}{
not used right now
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GraphqlClient$prep_query(query)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{query}}{(character) a query, of class `query` or `fragment`}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ghql-package.R
\docType{package}
\name{ghql-package}
\alias{ghql-package}
\alias{ghql}
\title{ghql}
\description{
General purpose GraphQL client
}
\section{ghql API}{

The main interface in this package is [GraphqlClient], which produces
a client (R6 class) with various methods for interacting with a
GraphQL server. [GraphqlClient] also accepts various input parameters
to set a base URL, and any headers required, which is usually the required
set of things needed to connect to a GraphQL service.

[Query] is an interface to creating GraphQL queries,
which works together with [GraphqlClient]

[Fragment] is an interface to creating GraphQL fragments,
which works together with [GraphqlClient]
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query-class.R
\name{Query}
\alias{Query}
\title{Query}
\value{
a `Query` class (R6 class)
}
\description{
ghql query class
}
\note{
we run an internal method `check_query()` that runs the public
method `parse2json()` - if the query doesn't pass the libgraphqlparser
parser, we return the error message
}
\examples{
# make a client
qry <- Query$new()

## define query
qry$query('query2', '{
  viewer {
    repositories(last: 10, isFork: false, privacy: PUBLIC) {
      edges {
        node {
          isPrivate
          id
          name
        }
      }
    }
  }
}')
qry
qry$queries
qry$queries$query2

# fragments
## by hand
qry$query('querywithfrag', '{
  ropensci: repositoryOwner(login:"ropensci") {
    repositories(first: 3) {
      edges {
        node {
          ...Watchers
        }
      }
    }
  }
  ropenscilabs: repositoryOwner(login:"ropenscilabs") {
    repositories(first: 3) {
      edges {
        node {
          ...Watchers
        }
      }
    }
  }
}
fragment Watchers on Repository {
  watchers(first: 3) {
    edges {
      node {
        name
      }
    }
  }
}')
qry
qry$queries
qry$queries$querywithfrag


\dontrun{
token <- Sys.getenv("GITHUB_GRAPHQL_TOKEN")
con <- GraphqlClient$new(
  url = "https://api.github.com/graphql",
  headers = list(Authorization = paste0("Bearer ", token))
)
jsonlite::fromJSON(con$exec(qry$queries$querywithfrag))

## use Fragment class fragments generator
### define query without fragment, but referring to it
qry$query('queryfrag', '{
  ropensci: repositoryOwner(login:"ropensci") {
    repositories(first: 3) {
      edges {
        node {
          ...Watchers
        }
      }
    }
  }
  ropenscilabs: repositoryOwner(login:"ropenscilabs") {
    repositories(first: 3) {
      edges {
        node {
          ...Watchers
        }
      }
    }
  }
}')

### define a fragment, and use it later
frag <- Fragment$new()
frag$fragment('Watchers', '
  fragment on Repository {
    watchers(first: 3) {
      edges {
        node {
          name
       }
    }
  }
}')
frag$fragments
frag$fragments$Watchers

### add the fragment to the query 'queryfrag'
qry$add_fragment('queryfrag', frag$fragments$Watchers)
qry
qry$queries
qry$queries$queryfrag
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{queries}}{(list) list of queries}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{Query$print()}}
\item \href{#method-query}{\code{Query$query()}}
\item \href{#method-add_fragment}{\code{Query$add_fragment()}}
\item \href{#method-parse2json}{\code{Query$parse2json()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the `Query` class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Query$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-query"></a>}}
\if{latex}{\out{\hypertarget{method-query}{}}}
\subsection{Method \code{query()}}{
define query in a character string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Query$query(name, x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{(character) name of the query}

\item{\code{x}}{(character) the query}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets query with `name` internally
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add_fragment"></a>}}
\if{latex}{\out{\hypertarget{method-add_fragment}{}}}
\subsection{Method \code{add_fragment()}}{
add a fragment to a query
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Query$add_fragment(query_name, fragment)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{query_name}}{(character) the query name to add the fragment to}

\item{\code{fragment}}{(character) the fragment itself}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets the fragment with the query
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse2json"></a>}}
\if{latex}{\out{\hypertarget{method-parse2json}{}}}
\subsection{Method \code{parse2json()}}{
parse query string with libgraphqlparser and get back JSON
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Query$parse2json(query, parse_schema = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{query}}{(character) a query to parse}

\item{\code{parse_schema}}{(logical) enable schema definition parsing?
default: `FAlSE`}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
adf
}
}
}
