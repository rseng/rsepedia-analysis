# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
[![Build Status](https://travis-ci.org/ropensci/chromer.svg?branch=master)](https://travis-ci.org/ropensci/chromer)   
[![Build status](https://ci.appveyor.com/api/projects/status/b1xjatd4i1gx1o6n?svg=true)](https://ci.appveyor.com/project/karthik/chromer)

![chromer logo](https://github.com/ropensci/chromer/raw/master/extra/logo.png)

This package provides programmatic access to the [Chromosome Counts Database (CCDB)](http://ccdb.tau.ac.il/home/) [API](http://ccdb.tau.ac.il/services/). The CCDB is a community resource for plant chromosome numbers. For more details on the database, see the associated publication by [Rice et al.](http://onlinelibrary.wiley.com/doi/10.1111/nph.13191/full) in *New Phytologist*. 

This package is maintained by [Paula Andrea Martinez](https://twitter.com/orchid00) and formerly maintained by [Matthew Pennell](http://mwpennell.github.io/) who are not affiliated with the CCDB group. The URL for Chromer docs is [https://docs.ropensci.org/chromer/](https://docs.ropensci.org/chromer/). 

## Installing
The package can be installed directly from CRAN, but it is currently outdated -- PLEASE install directly from GitHub

```r
install.packages("chromer")
```

or, for the latest version, you can install directly from GitHub using [devtools](http://github.com/hadley/devtools)

```r
## install.packages("devtools")
devtools::install_github("ropensci/chromer")
```

## Querying the CCDB

It is possible to query the database in three ways: by `species`, `genus`, `family`, and `majorGroup`. For example, if we are interested in the genus *Solanum* (Solanaceae), which contains the potato, tomato, and eggplant, we would query the database as follows

```r
library(chromer)
sol_gen <- chrom_counts(taxa = "Solanum", rank = "genus")
head(sol_gen)
nrow(sol_gen)
```

There are over 3000 records for Solanum alone! If we are interested in a particular species, such as tomatoes, we can search for the species directly. 

```r
sol_tom <- chrom_counts(taxa = "Solanum_lycopersicum", rank = "species")
head(sol_tom)
```

Note that `taxa="Solanum lycopersicum"` (including a space between the genus and species name) will also work here.

If we wanted to get data on the whole family, we simply type

```r
sol_fam <- chrom_counts(taxa = "Solanaceae", rank = "family")
head(sol_fam)
```

Or, expand the scope much further and get all Angiosperms (this will take some time)

```r
ang <- chrom_counts(taxa = "Angiosperms", rank = "majorGroup")
head(ang)
```

There are two options for returning data. The first (default) is to only return the species name information (including taxonomic resolutions made by [Taxonome](http://taxonome.bitbucket.org/)) and the haploid and diploid counts. Setting the argument 
`full=TRUE`

```r
sol_gen_full <- chrom_counts("Solanum", rank = "genus", full = TRUE)
```

returns a bunch more info on the records.

```r
head(sol_gen_full)
```

## Summarizing the data

The Chromosome Counts Database is a fantastic resource but as it is a compilation of a large number of resources and studies, the data is somewhat messy and challenging to work with. We have written a little function that does some post-processing to make it easier to handle. The function `summarize_counts()` does the following:

1. Aggregates multiple records for the same species

2. Infers the gametophytic (haploid) number of chromosomes when only the sporophytic (diploid) counts are available. 

3. Parses the records for numeric values. In some cases chromosomal counts also include text characters (e.g., #-#; c.#; #,#,#; and many other varieties). As there are many possible ways that chromosomal counts may be listed in the database, the function takes the naive approach and simply searches the strings for integers. In most cases, this is sensible but may produces weird results on occasion. **Some degree of manual curation will probably be necessary and the output of the summary should be used with caution in downstream analyses**.

To summarize and clean the count data obtained from `chrom_counts()` simply use
```
summarize_counts(sol_gen)
``` 

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/chromer/issues).
* License: CC0
* Get citation information for `chromer` in R doing `citation(package = "chromer")`
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean-data.R
\name{summarize_counts}
\alias{summarize_counts}
\title{Summarize chromosome counts from API call}
\usage{
summarize_counts(counts)
}
\arguments{
\item{counts}{A 'chrom.counts' object inherited from \code{\link{chrom_counts}}.}
}
\value{
A \code{data.frame} containing the resolved binomial, the count type (gametophytic or sporophytic), the counts, the inferred gametophytic count (for sporophytic records) and the number of records supporting each count.
}
\description{
This function processes and cleans the data returned from the API call for use in downstream analysis.
}
\details{
The results from the API call are a bit messy and difficult to use for downstream analyses. This function cleans up the data in three ways. First, it combines aggregates and summarizes all records from each species. Second, many of the counts are combined with text characters (e.g., \code{"#-#"}, \code{"c.#"}, and \code{"#, #, #"}. This function uses regular expressions to pull out all and any numeric values from these strings. Third, some of the records are gametophytic (n) counts and others are from sporophytes (2n); the function simply divides the sporophytic counts in half so that all measurements are on a common scale.

IMPORTANT: Use this function with caution. Parsing the counts programmatically may be useful but it may generate erroneous results in some cases if input is in an odd format. For example, if the count is \code{"#+-#"}, the function will return both the first and second \code{#} as valid counts . Given the creativity(?) of researchers in entering data, it is hard to predict all possible ways that the counts may be represented. Therefore, some manual checking will probably be necessary.
}
\examples{
\dontrun{

## Get all counts for genus Castilleja
res <- chrom_counts("Castilleja", "genus")

## summarize the results
summarize_counts(res)

}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counts.R
\name{chrom_counts}
\alias{chrom_counts}
\title{Returns chromosome counts from Chromosome Counts Database API}
\usage{
chrom_counts(taxa, rank = c("species", "genus", "family", "majorGroup"),
  full = FALSE, foptions = list())
}
\arguments{
\item{taxa}{Taxonomic name(s) to query. Can be either a single name, a vector of multiple names or a list. If supplying multiple names, these must all be of the same \code{rank}.}

\item{rank}{Rank to query.}

\item{full}{Logical. Whether to return full records. Defaults to \code{FALSE} which returns only partial records. Partial records includes the resolved name as well as the gametophytic (n) and sporophytic (2n) counts.}

\item{foptions}{additional options to be passed to \code{httr::GET}}
}
\value{
A \code{data.frame} containing all records matched by query
}
\description{
This function calls the Chromsome Counts Database (CCDB) API and returns all counts for specified higher taxa.
}
\details{
When using the API to query for species, both matched names and resolved names are searched. This means that all records for potential synonyms will be returned as well. Currently species binomials must be specified by either 'genus species' (i.e., space between genus and species) or 'genus_species'.

To search for subspecies (subsp.) or varieties (var.), you can use search terms like:

\code{"Solanum acaule var. albicans"}.

Searching for \code{"Solanum acaule"} will return all subspecies and varieties.

Currently the only acceptable search terms when specifying \code{"majorGroup"} are \code{"Angiosperms"}, \code{"Gymnosperms"}, \code{"Pteridophytes"}, or \code{"Bryophytes"}.
}
\examples{
\dontrun{

## Get all counts for genus Castilleja
chrom_counts("Castilleja", "genus")

## Get all counts for both Castilleja and Lachemilla
chrom_counts(c("Castilleja", "Lachemilla"), "genus")

## Get all counts for Castilleja miniata
chrom_counts("Castilleja miniata", "species")

## Get all counts for only Castilleja miniata subsp. elata
chrom_counts("Castilleja miniata subsp. elata", "species")

## Note that searching for "Castilleja miniata" will return
## all subspecies and varieties

## Get all counts for the Orobanchaceae
chrom_counts("Orobanchaceae", "family")

}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chromer-package.R
\docType{package}
\name{chromer}
\alias{chromer}
\alias{chromer-chromer}
\alias{chromer-package}
\title{Chromer: Programmatic access to Chromosome Counts Database}
\description{
This package is a R based interface to the Chromosome Counts Database (CCDB). Currently, this consists of just one function \code{\link{chrom_counts}}, which queries the database using the taxa supplied. The function \code{\link{summarize_counts}} cleans and processes the data.
}
\references{
Rice, A., L. Glick, S. Abadi, M. Einhorn, N.M. Kopelman, A. Salman-Minkov, J. Mayzel, O. Chay, and I. Mayrose. 2014. The Chromosome Counts Database (CCDB) -- a community resource of plant chromosome numbers. New Phytologist doi:10.1111/nph.13191.
}

