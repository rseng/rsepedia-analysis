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
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/ropensci/rentrez.png)](https://travis-ci.org/ropensci/rentrez) [![Build status](https://ci.appveyor.com/api/projects/status/y8mq2v4mpgou8rhp/branch/master)](https://ci.appveyor.com/project/sckott/rentrez/branch/master) [![Coverage Status](https://coveralls.io/repos/ropensci/rentrez/badge.svg?branch=master)](https://coveralls.io/r/ropensci/rentrez?branch=master) [![CRAN](http://cranlogs.r-pkg.org/badges/rentrez)](http://cran.rstudio.com/web/packages/rentrez/index.html) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32420.svg)](http://dx.doi.org/10.5281/zenodo.32420)

rentrez
=======

`rentrez` provides functions that work with the [NCBI Eutils](http://www.ncbi.nlm.nih.gov/books/NBK25500/) API to search, download data from, and otherwise interact with NCBI databases.

Install
-------

`rentrez` is on CRAN, so you can get the latest stable release with `install.packages("rentrez")`. This repository will sometimes be a little ahead of the CRAN version, if you want the latest (and possibly greatest) version you can install the current github version using Hadley Wickham's [devtools](https://github.com/hadley/devtools).

    library(devtools)
    install_github("ropensci/rentrez")

The EUtils API
--------------

Each of the functions exported by `rentrez` is documented, and this README and the package vignette provide examples of how to use the functions together as part of a workflow. The API itself is [well-documented](http://www.ncbi.nlm.nih.gov/books/NBK25500/). Be sure to read the official documentation to get the most out of API. In particular, be aware of the NCBI's usage policies and try to limit very large requests to off peak (USA) times (`rentrez` takes care of limiting the number of requests per second, and setting the appropriate entrez tool name in each request).

Hopefully this README, and the package's vignette and in-line documentation, provide you with enough information to get started with `rentrez`. If you need more help, or if you discover a bug in `rentrez` please let us know, either through one of the [contact methods described here](http://ropensci.org/contact.html), or [by filing an issue](https://github.com/ropensci/rentrez/issues)

Examples
--------

In many cases, doing something interesting with `EUtils` will take multiple calls. Here are a few examples of how the functions work together (check out the package vignette for others).

### Getting data from that great paper you've just read

Let's say I've just read a paper on the evolution of Hox genes, [Di-Poi *et al*. (2010)](dx.doi.org/10.1038/nature08789), and I want to get the data required to replicate their results. First, I need the unique ID for this paper in pubmed (the PMID). Unfortunately, many journals don't give PMIDS for their papers, but we can use `entrez_search` to find the paper using the doi field:

``` r
library(rentrez)
hox_paper <- entrez_search(db="pubmed", term="10.1038/nature08789[doi]")
hox_paper$ids
```

    # [1] "20203609"

Now, what sorts of data are available from other NCBI database for this paper?

``` r
hox_data <- entrez_link(db="all", id=hox_paper$ids, dbfrom="pubmed")
hox_data
```

    # elink object with contents:
    #  $links: IDs for linked records from NCBI
    # 

In this case all the data is in the `links` element:

``` r
hox_data$links
```

    # elink result with information from 14 databases:
    #  [1] pubmed_medgen              pubmed_pmc_refs           
    #  [3] pubmed_pubmed              pubmed_pubmed_alsoviewed  
    #  [5] pubmed_pubmed_citedin      pubmed_pubmed_combined    
    #  [7] pubmed_pubmed_five         pubmed_pubmed_reviews     
    #  [9] pubmed_pubmed_reviews_five pubmed_mesh_major         
    # [11] pubmed_nuccore             pubmed_nucleotide         
    # [13] pubmed_protein             pubmed_taxonomy_entrez

Each of the character vectors in this object contain unique IDs for records in the named databases. These functions try to make the most useful bits of the returned files available to users, but they also return the original file in case you want to dive into the XML yourself.

In this case we'll get the protein sequences as fasta files, using ' `entrez_fetch`:

``` r
hox_proteins <- entrez_fetch(db="protein", id=hox_data$links$pubmed_protein, rettype="fasta")
```

    # No encoding supplied: defaulting to UTF-8.

``` r
cat(substr(hox_proteins, 1, 237))
```

    # >gi|290760438|gb|ADD54588.1| HOXA10, partial [Saiphos equalis]
    # MACSESPAANSFLVDSLISSASVRGEGGGGGGGGGGAGGGGGEGGGGGGGVYYPNNSSVYLPQTSELSYG
    # LPSYGLFPVLSKRNEGPSQSMVPASHTYMSGMEVWLDPPRSCRLEDPESPQATSCSFTPNIKEENSYCLY
    # DSDKGPKEATATDLSTFPRLTSEVCSMNNV

### Retrieving datasets associated a particular organism.

I like spiders. So let's say I want to learn a little more about New Zealand's endemic "black widow" the katipo. Specifically, in the past the katipo has been split into two species, can we make a phylogeny to test this idea?

The first step here is to use the function `entrez_search` to find datasets that include katipo sequences. The `popset` database has sequences arising from phylogenetic or population-level studies, so let's start there.

``` r
library(rentrez)
katipo_search <- entrez_search(db="popset", term="Latrodectus katipo[Organism]")
katipo_search$count
```

    # [1] 6

In this search `count` is the total number of hits returned for the search term. We can use `entrez_summary` to learn a little about these datasets. `rentrez` will parse this xml into a list of `esummary` records, with each list entry corresponding to one of the IDs it is passed. In this case we get six records, and we see what each one contains like so:

``` r
katipo_summs <- entrez_summary(db="popset", id=katipo_search$ids)
katipo_summs
```

    # List of  6 esummary records. First record:
    # 
    #  $`167843272`
    # esummary result with 17 items:
    #  [1] uid        caption    title      extra      gi         settype   
    #  [7] createdate updatedate flags      taxid      authors    article   
    # [13] journal    strain     statistics properties oslt

An we can extract specific elements from list of summary records with `extract_from_esummary`:

``` r
titles <- extract_from_esummary(katipo_summs, "title")
unname(titles)
```

    # [1] "Latrodectus katipo 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence."
    # [2] "Latrodectus katipo cytochrome oxidase subunit 1 (COI) gene, partial cds; mitochondrial."                                                                                                                                 
    # [3] "Latrodectus 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence."       
    # [4] "Latrodectus cytochrome oxidase subunit 1 (COI) gene, partial cds; mitochondrial."                                                                                                                                        
    # [5] "Latrodectus tRNA-Leu (trnL) gene, partial sequence; and NADH dehydrogenase subunit 1 (ND1) gene, partial cds; mitochondrial."                                                                                            
    # [6] "Theridiidae cytochrome oxidase subunit I (COI) gene, partial cds; mitochondrial."

Let's just get the two mitochondrial loci (COI and trnL), using `entrez_fetch`:

``` r
COI_ids <- katipo_search$ids[c(2,6)]
trnL_ids <- katipo_search$ids[5]
COI <- entrez_fetch(db="popset", id=COI_ids, rettype="fasta")
```

    # No encoding supplied: defaulting to UTF-8.

``` r
trnL <- entrez_fetch(db="popset", id=trnL_ids, rettype="fasta")
```

    # No encoding supplied: defaulting to UTF-8.

The "fetched" results are fasta formatted characters, which can be written to disk easily:

``` r
write(COI, "Test/COI.fasta")
write(trnL, "Test/trnL.fasta")
```

Once you've got the sequences you can do what you want with them, but I wanted a phylogeny and we can do that entirly within R. To get a nice tree with legible tip labels I'm gong to use `stringr` to extract just the species names and `ape` to built and root and neighbor joining tree:

``` r
library(ape)
tf <- tempfile()
write(COI, tf)
coi <- read.dna(tf, format="fasta")
coi_aligned <- muscle(coi)
tree <- nj(dist.dna(coi_aligned))
tree$tip.label <- stringr::str_extract(tree$tip.label, "Steatoda [a-z]+|Latrodectus [a-z]+")
plot( root(tree, outgroup="Steatoda grossa" ), cex=0.8)
```

![](http://i.imgur.com/8n9UeIi.png)

### web\_history and big queries

The NCBI provides search history features, which can be useful for dealing with large lists of IDs or repeated searches.

As an example, imagine you wanted to learn something about all of the SNPs in the non-recombing portion of the Y chromsome in humans. You could first find these SNPs using `entrez_search`, using the "CHR" (chromosome) and "CPOS" (position in chromosome) to specify the region of interest. (The syntax for these search terms is described in the vignette and the documentation for `entrez_search`):

``` r
snp_search <- entrez_search(db="snp", 
                            term="(Y[CHR] AND Homo[ORGN]) NOT 10001:2781479[CPOS]")
snp_search
```

    # Entrez search result with 234154 hits (object contains 20 IDs and no web_history object)
    #  Search term (as translated):  (Y[CHR] AND "Homo"[Organism]) NOT 10001[CHRPOS] :  ...

When I wrote this that was a little over 200 000 SNPs. It's probably not a good idea to set `retmax` to 200 000 and just download all of those identifiers. Instead, we could store this list of IDs on the NCBI's server and refer to them in later calles to functions like `entrez_link` and `entrez_fetch` that accept a web history object.

``` r
snp_search <- entrez_search(db="snp", 
                            term="(Y[CHR] AND Homo[ORGN]) NOT 10001:2781479[CPOS]", 
                            use_history = TRUE)
snp_search
```

    # Entrez search result with 234154 hits (object contains 20 IDs and a web_history object)
    #  Search term (as translated):  (Y[CHR] AND "Homo"[Organism]) NOT 10001[CHRPOS] :  ...

As you can see, the result of the search now includes a `web_history` object. We can use that object to refer to these IDs in later calls. Heree we will just fetch complete records of the first 5 SNPs.

``` r
recs <- entrez_fetch(db="snp", web_history=snp_search$web_history, retmax=5, rettype="xml", parsed=TRUE)
class(recs)
```

    # [1] "XMLInternalDocument" "XMLAbstractDocument"

The records come to us as parsed XML objects, which you could futher process with the `XML` library or write to disk for later use.

### Getting information about NCBI databases

Most of the examples above required some background information about what databases NCBI has to offer, and how they can be searched. `rentrez` provides a set of functions with names starting `entrez_db` that help you to discover this information in an interactive session.

First up, `entrez_dbs()` gives you a list of database names

``` r
entrez_dbs()
```

    #  [1] "pubmed"          "protein"         "nuccore"        
    #  [4] "nucleotide"      "nucgss"          "nucest"         
    #  [7] "structure"       "genome"          "annotinfo"      
    # [10] "assembly"        "bioproject"      "biosample"      
    # [13] "blastdbinfo"     "books"           "cdd"            
    # [16] "clinvar"         "clone"           "gap"            
    # [19] "gapplus"         "grasp"           "dbvar"          
    # [22] "epigenomics"     "gene"            "gds"            
    # [25] "geoprofiles"     "homologene"      "medgen"         
    # [28] "mesh"            "ncbisearch"      "nlmcatalog"     
    # [31] "omim"            "orgtrack"        "pmc"            
    # [34] "popset"          "probe"           "proteinclusters"
    # [37] "pcassay"         "biosystems"      "pccompound"     
    # [40] "pcsubstance"     "pubmedhealth"    "seqannot"       
    # [43] "snp"             "sra"             "taxonomy"       
    # [46] "unigene"         "gencoll"         "gtr"

Some of the names are a little opaque, so you can get some more descriptive information about each with `entrez_db_summary()`

``` r
entrez_db_summary("cdd")
```

    #  DbName: cdd
    #  MenuName: Conserved Domains
    #  Description: Conserved Domain Database
    #  DbBuild: Build150814-1106.1
    #  Count: 50648
    #  LastUpdate: 2015/08/14 18:42

`entrez_db_searchable()` lets you discover the fields available for search terms for a given database. You get back a named-list, with names are fields. Each element has additional information about each named search field (you can also use `as.data.frame` to create a dataframe, with one search-field per row):

``` r
search_fields <- entrez_db_searchable("pmc")
search_fields$GRNT
```

    #  Name: GRNT
    #  FullName: Grant Number
    #  Description: NIH Grant Numbers
    #  TermCount: 2272841
    #  IsDate: N
    #  IsNumerical: N
    #  SingleToken: Y
    #  Hierarchy: N
    #  IsHidden: N

Finally, `entrez_db_links` takes a database name, and returns a list of other NCBI databases which might contain linked-records.

``` r
entrez_db_links("omim")
```

    # Databases with linked records for database 'omim'
    #  [1] biosample   biosystems  books       clinvar     dbvar      
    #  [6] gene        genetests   geoprofiles gtr         homologene 
    # [11] mapview     medgen      medgen      nuccore     nucest     
    # [16] nucgss      omim        pcassay     pccompound  pcsubstance
    # [21] pmc         protein     pubmed      pubmed      sra        
    # [26] structure   unigene

### Trendy topics in genetics

This is one is a little more trivial, but you can also use entrez to search pubmed and the EUtils API allows you to limit searches by the year in which the paper was published. That gives is a chance to find the trendiest -omics going around (this has quite a lot of repeated searching, so it you want to run your own version be sure to do it in off peak times).

Let's start by making a function that finds the number of records matching a given search term for each of several years (using the `mindate` and `maxdate` terms from the Eutils API):

    library(rentrez)
    papers_by_year <- function(years, search_term){
        return(sapply(years, function(y) entrez_search(db="pubmed",term=search_term, mindate=y, maxdate=y, retmax=0)$count))
    }

With that we can fetch the data for each term and, by searching with no term, find the total number of papers published in each year:

``` r
years <- 1990:2015
total_papers <- papers_by_year(years, "")
omics <- c("genomic", "epigenomic", "metagenomic", "proteomic", "transcriptomic", "pharmacogenomic", "connectomic" )
trend_data <- sapply(omics, function(t) papers_by_year(years, t))
trend_props <- trend_data/total_papers
```

That's the data, let's plot it:

``` r
library(reshape)
library(ggplot2)
trend_df <- melt(data.frame(years, trend_props), id.vars="years")
p <- ggplot(trend_df, aes(years, value, colour=variable))
p + geom_line(size=1) + scale_y_log10("number of papers")
```

Giving us... well this:

![](http://i.imgur.com/oSYuWqz.png)

------------------------------------------------------------------------

This package is part of a richer suite called [fulltext](https://github.com/ropensci/fulltext), along with several other packages, that provides the ability to search for and retrieve full text of open access scholarly articles.

------------------------------------------------------------------------

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# CONTRIBUTING 

## Please contribute!
We love collaboration.

## Bugs?

* Submit an issue on the Issues page [here](https://github.com/ropensci/rentrez/issues)

## Code/Documentation contributions

###Scope 

We would love to have your help in developing `rentrez`. Bear in mind, `rentrez`
is intended to be a general-purpose wrapper for the `Eutils` API. Functions that
make use of `rentrez` to perform more specific tasks should be added to existing
packages, or form the basis of new ones. 

###Perferred way to contribute code

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rentrez.git`
* Make sure to track progress upstream (i.e., on our version of `rentrez` at `ropensci/rentrez`) by doing `git remote add upstream https://github.com/ropensci/rentrez.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rentrez`

### Questions? 

Get in touch: [david.winter@gmail.com](mailto:david.winter@gmail.com).

Please note that this project is released with a [Contributor Code of
Conduct](CONDUCT.md). By participating in this project you agree to abide by its
terms.

 Thanks for contibuting!

---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

Sometimes `rentrez` will throw obscure error messages because the NCBI's webservers are down or otherwise not behaving as they are meant to. Before reporting issues about errors thrown while contacting an NCBI database (search, fetch, summary,	 link...) you can test the NCBI is responding to requests by pasting the following URLS into a web browser.	

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=6060535&retmote=rsr	
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=rentrez	

They should give you a plain text file and a small XML file respectively. If you receive errors when trying to aces the URLs it is likely they the NCBI is	having an intermittent problem, try your request again in a few minutes. If those pages load as expected then please file your issue. We appreciate user	
feedback and will do our best to help.
---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#"
)
knitr::opts_knit$set(upload.fun = knitr::imgur_upload, base.url = NULL)  

```

[![Build Status](https://travis-ci.org/ropensci/rentrez.png)](https://travis-ci.org/ropensci/rentrez)
[![Build status](https://ci.appveyor.com/api/projects/status/y8mq2v4mpgou8rhp/branch/master)](https://ci.appveyor.com/project/sckott/rentrez/branch/master)
[![Coverage Status](https://coveralls.io/repos/ropensci/rentrez/badge.svg?branch=master)](https://coveralls.io/r/ropensci/rentrez?branch=master)
[![CRAN](http://cranlogs.r-pkg.org/badges/rentrez)](http://cran.rstudio.com/web/packages/rentrez/index.html)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32420.svg)](http://dx.doi.org/10.5281/zenodo.32420)


#rentrez

`rentrez` provides functions that work with the [NCBI Eutils](http://www.ncbi.nlm.nih.gov/books/NBK25500/)
API to search, download data from, and otherwise interact with NCBI databases.


##Install

`rentrez` is on CRAN, so you can get the latest stable release with `install.packages("rentrez")`. This repository will sometimes be a little ahead of the CRAN version, if you want the latest (and possibly greatest) version you can install
the current github version using Hadley Wickham's [devtools](https://github.com/hadley/devtools).

```
library(devtools)
install_github("ropensci/rentrez")
```



##The EUtils API

Each of the functions exported by `rentrez` is documented, and this README and
the package vignette provide examples of how to use the functions together as part
of a workflow. The API itself is [well-documented](http://www.ncbi.nlm.nih.gov/books/NBK25500/).
Be sure to read the official documentation to get the most out of API. In particular, be aware of the NCBI's usage
policies and try to limit very large requests to off peak (USA) times (`rentrez`
takes care of limiting the number of requests per second, and setting the
appropriate entrez tool name in each request).

Hopefully this README, and the package's vignette and in-line documentation,
provide you with enough information to get started with `rentrez`. If you need
more help, or if you discover a bug in `rentrez` please let us know, either through
one of the [contact methods described here](http://ropensci.org/contact.html),
or [by filing an issue](https://github.com/ropensci/rentrez/issues)


##Examples

In many cases, doing something interesting with `EUtils` will take multiple
calls. Here are a few examples of how the functions work together (check out the
package vignette for others).

###Getting data from that great paper you've just read

Let's say I've just read a paper on the evolution of Hox genes,
[Di-Poi _et al_. (2010)](dx.doi.org/10.1038/nature08789), and I want to get the
data required to replicate their results. First, I need the unique ID for this
paper in pubmed (the PMID). Unfortunately, many journals don't give PMIDS for their
papers, but we can use `entrez_search` to find the paper using the doi field:

```{r doi}
library(rentrez)
hox_paper <- entrez_search(db="pubmed", term="10.1038/nature08789[doi]")
hox_paper$ids
```

Now, what sorts of data are available from other NCBI database for this paper?

```{r links}
hox_data <- entrez_link(db="all", id=hox_paper$ids, dbfrom="pubmed")
hox_data
```
In this case all the data is in the `links` element:

```{r showlinks}
hox_data$links
```

Each of the character vectors in this object contain unique IDs for records in
the named databases. These functions try to make the most useful bits of the
returned files available to users, but they also return the original file in case
you want to dive into the XML yourself.

In this case we'll get the protein sequences as fasta files, using '
`entrez_fetch`:

```{r proteins}
hox_proteins <- entrez_fetch(db="protein", id=hox_data$links$pubmed_protein, rettype="fasta")
cat(substr(hox_proteins, 1, 237))
```

###Retrieving datasets associated a particular organism.

I like spiders. So let's say I want to learn a little more about New Zealand's
endemic "black widow" the katipo. Specifically, in the past the katipo has
been split into two species, can we make a phylogeny to test this idea?

The first step here is to use the function `entrez_search` to find datasets
that include katipo sequences. The `popset` database has sequences arising from
phylogenetic or population-level studies, so let's start there.

```{r katipo}
library(rentrez)
katipo_search <- entrez_search(db="popset", term="Latrodectus katipo[Organism]")
katipo_search$count
```

In this search `count` is the total number of hits returned for the search term.
We can use `entrez_summary` to learn a little about these datasets. `rentrez`
will parse this xml into a list of `esummary` records, with each list entry
corresponding to one of the IDs it is passed. In this case we get six records,
and we see what each one contains like so:


```{r summ}
katipo_summs <- entrez_summary(db="popset", id=katipo_search$ids)
katipo_summs
```
An  we can extract specific elements from list of summary records with
`extract_from_esummary`:

```{r extract}
titles <- extract_from_esummary(katipo_summs, "title")
unname(titles)
```

Let's just get the two mitochondrial loci (COI and trnL), using `entrez_fetch`:

```{r fetch}
COI_ids <- katipo_search$ids[c(2,6)]
trnL_ids <- katipo_search$ids[5]
COI <- entrez_fetch(db="popset", id=COI_ids, rettype="fasta")
trnL <- entrez_fetch(db="popset", id=trnL_ids, rettype="fasta")
```

The "fetched" results are fasta formatted characters, which can be written
to disk easily:

```r
write(COI, "Test/COI.fasta")
write(trnL, "Test/trnL.fasta")
```

Once you've got the sequences you can do what you want with them, but I wanted
a phylogeny and we can do that entirly within R. 
To get a nice tree with legible tip labels I'm gong to
use `stringr` to extract just the species names and `ape` to built and root and
neighbor joining tree:

```r
library(ape)
tf <- tempfile()
write(COI, tf)
coi <- read.dna(tf, format="fasta")
coi_aligned <- muscle(coi)
tree <- nj(dist.dna(coi_aligned))
tree$tip.label <- stringr::str_extract(tree$tip.label, "Steatoda [a-z]+|Latrodectus [a-z]+")
plot( root(tree, outgroup="Steatoda grossa" ), cex=0.8)
```

![](http://i.imgur.com/8n9UeIi.png)


### web_history and big queries

The NCBI provides search history features, which can be useful for dealing with 
large lists of IDs or repeated searches. 

As an example, imagine you wanted to learn something about all of the SNPs in 
the non-recombing portion of the Y chromsome in humans. 
You could first find these SNPs using `entrez_search`, using the "CHR"
(chromosome) and "CPOS" (position in chromosome) to specify the region of
interest. (The syntax for these search terms is described in the vignette and
the documentation for `entrez_search`):


```{r snp_search}
snp_search <- entrez_search(db="snp", 
                            term="(Y[CHR] AND Homo[ORGN]) NOT 10001:2781479[CPOS]")
snp_search
```

When I wrote this that was a little over 200 000 SNPs. It's probably not a good
idea to set `retmax` to 200 000 and just download all of those identifiers.
Instead, we could store this list of IDs on the NCBI's server and refer to them
in later calles to functions like `entrez_link` and `entrez_fetch` that accept
a web history object. 



```{r snp_history}
snp_search <- entrez_search(db="snp", 
                            term="(Y[CHR] AND Homo[ORGN]) NOT 10001:2781479[CPOS]", 
                            use_history = TRUE)
snp_search
```

As you can see, the result of the search now includes a `web_history` object. We can
use that object to refer to these IDs in later calls. Heree we will just fetch 
complete records of the first 5 SNPs.

```{r snp_fetch}
recs <- entrez_fetch(db="snp", web_history=snp_search$web_history, retmax=5, rettype="xml", parsed=TRUE)
class(recs)
```

The records come to us as parsed XML objects, which you could futher process
with the `XML` library or write to disk for later use.


###Getting information about NCBI databases

Most of the examples above required some background information about what
databases NCBI has to offer, and how they can be searched. `rentrez` provides
a set of functions with names starting `entrez_db` that help you to discover
this information in an interactive session.

First up, `entrez_dbs()` gives you a list of database names


```{r dbs}
entrez_dbs()
```

Some of the names are a little opaque, so you can get some more descriptive
information about each with `entrez_db_summary()`

```{r summary}
entrez_db_summary("cdd")
```

`entrez_db_searchable()` lets you discover the fields available for search terms
for a given database. You get back a named-list, with names are fields. Each
element has additional information about each named search field (you can also
use `as.data.frame` to create a dataframe, with one search-field per row):

```{r fields}
search_fields <- entrez_db_searchable("pmc")
search_fields$GRNT
```

Finally, `entrez_db_links` takes a database name, and returns a list of other
NCBI databases which might contain linked-records.

```{r elinks}
entrez_db_links("omim")
```

###Trendy topics in genetics

This is one is a little more trivial, but you can also use entrez to search pubmed and
the EUtils API allows you to limit searches by the year in which the paper was published.
That gives is a chance to find the trendiest -omics going around (this has quite a lot
of repeated searching, so it you want to run your own version be sure to do it
in off peak times).

Let's start by making a function that finds the number of records matching a given
search term for each of several years (using the `mindate` and `maxdate` terms from
the Eutils API):

```
library(rentrez)
papers_by_year <- function(years, search_term){
    return(sapply(years, function(y) entrez_search(db="pubmed",term=search_term, mindate=y, maxdate=y, retmax=0)$count))
}
```

With that we can fetch the data for each term and, by searching with no term,
find the total number of papers published in each year:


```r
years <- 1990:2015
total_papers <- papers_by_year(years, "")
omics <- c("genomic", "epigenomic", "metagenomic", "proteomic", "transcriptomic", "pharmacogenomic", "connectomic" )
trend_data <- sapply(omics, function(t) papers_by_year(years, t))
trend_props <- trend_data/total_papers
```

That's the data, let's plot it:

```r
library(reshape)
library(ggplot2)
trend_df <- melt(data.frame(years, trend_props), id.vars="years")
p <- ggplot(trend_df, aes(years, value, colour=variable))
p + geom_line(size=1) + scale_y_log10("number of papers")
```


Giving us... well this:

![](http://i.imgur.com/oSYuWqz.png)



---

This package is part of a richer suite called [fulltext](https://github.com/ropensci/fulltext), along with several other packages, that provides the ability to search for and retrieve full text of open access scholarly articles.

---

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: Rentrez Tutorial
author: "David winter"
date: "`r Sys.Date()`"
output: 
    rmarkdown::html_vignette:
        toc: true
vignette: >
  %\VignetteIndexEntry{Rentrez Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, count_recs, echo=FALSE}
#are we on CRAN?
library(testthat)
build_this <- identical(Sys.getenv("BUILD_RENTREZ_VIGNETTE"), "true")

knitr::opts_chunk$set(eval = build_this) 


library(rentrez)
count_recs <- function(db, denom) {
    Sys.sleep(0.3)
    nrecs <-  rentrez::entrez_db_summary(db)["Count"]
    round(as.integer(nrecs)/denom, 1)
}
```
## Introduction: The NCBI, entrez and `rentrez`.

The NCBI shares a _lot_ of data. At the time this document was compiled, there
were 31.7 million papers in [PubMed](https://pubmed.ncbi.nlm.nih.gov/),
including 6.6 million full-text records available in [PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/).
[The NCBI Nucleotide Database](https://www.ncbi.nlm.nih.gov/nuccore) (which includes GenBank) has data for 432
million different sequences, and [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) describes
702 million different genetic variants. All of these
records can be cross-referenced with the  1.86 million
species in the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) or 27 thousand disease-associated records
in [OMIM](https://www.ncbi.nlm.nih.gov/omim).


The NCBI makes this data available through a [web interface](https://www.ncbi.nlm.nih.gov/),
an [FTP server](ftp://ftp.ncbi.nlm.nih.gov/) and through a REST API called the
[Entrez Utilities](https://www.ncbi.nlm.nih.gov/books/NBK25500/) (`Eutils` for
short). This package provides functions to use that API, allowing users to
gather and combine data from multiple NCBI databases in the comfort of an R
session or script.

## Getting started with the rentrez

To make the most of all the data the NCBI shares you need to know a little about
their databases, the records they contain and the ways you can find those
records. The [NCBI provides extensive documentation for each of their
databases](https://www.ncbi.nlm.nih.gov/home/documentation/) and for the
[EUtils API that `rentrez` takes advantage of](https://www.ncbi.nlm.nih.gov/books/NBK25501/).
There are also some helper functions in `rentrez` that help users learn their
way around the NCBI's databases. 

First, you can use `entrez_dbs()` to find the list of available databases:

```{r, dbs}
entrez_dbs()
```
There is a set of functions with names starting `entrez_db_` that can be used to
gather more information about each of these databases:

**Functions that help you learn about NCBI databases**

| Function name            | Return                                               |
|--------------------------|------------------------------------------------------|
| `entrez_db_summary()`    | Brief description of what the database is            |
| `entrez_db_searchable()` | Set of search terms that can used with this database |
| `entrez_db_links() `     | Set of databases that might contain linked records   |

For instance, we can get a description of the somewhat cryptically named
database 'cdd'...

```{r, cdd}
entrez_db_summary("cdd")
```

... or find out which search terms can be used with the Sequence Read Archive (SRA)
database (which contains raw data from sequencing projects):

```{r, sra_eg}
entrez_db_searchable("sra")
```

Just how these 'helper' functions might be useful will become clearer once
you've started using `rentrez`, so let's get started.

## Searching databases: `entrez_search()`

Very often, the first thing you'll want to do with `rentrez` is search a given
NCBI database to find records that match some keywords. You can do this using
the function `entrez_search()`. In the simplest case you just need to provide a 
database name (`db`) and a search term (`term`) so let's search PubMed for
articles about the `R language`:


```{r eg_search}
r_search <- entrez_search(db="pubmed", term="R Language")
```
The object returned by a search acts like a list, and you can get a summary of 
its contents by printing it.

```{r print_search}
r_search
```

There are a few things to note here. First, the NCBI's server has worked out
that we meant R as a programming language, and so included the
['MeSH' term](https://www.ncbi.nlm.nih.gov/mesh) term associated with programming
languages. We'll worry about MeSH terms and other special queries later, for now
just note that you can use this feature to check that your search term was interpreted in the way
you intended. Second, there are many more 'hits' for this search than there
are unique IDs contained in this object. That's because the optional argument
`retmax`, which controls the maximum number of returned values has a default
value of 20.

The IDs are the most important thing returned here. They
allow us to fetch records matching those IDs, gather summary data about them or find
cross-referenced records in other databases. We access the IDs as a vector using the
`$` operator:


```{r search_ids}
r_search$ids
```

If we want to get more than 20 IDs we can do so by increasing the `ret_max` argument.

```{r searchids_2}
another_r_search <- entrez_search(db="pubmed", term="R Language", retmax=40)
another_r_search
```

If we want to get IDs for all of the thousands of records that match this
search, we can use the NCBI's web history feature [described below](#web_history).


### Building search terms

The EUtils API uses a special syntax to build search terms. You can search a
database against a specific term using the format `query[SEARCH FIELD]`, and
combine multiple such searches using the boolean operators `AND`, `OR` and `NOT`.

For instance, we can find next generation sequence datasets for the (amazing...) ciliate
_Tetrahymena thermophila_ by using the organism ('ORGN') search field:


```{r, Tt}
entrez_search(db="sra",
              term="Tetrahymena thermophila[ORGN]",
              retmax=0)
```

We can narrow our focus to only those records that have been added recently (using the colon to 
specify a range of values):


```{r, Tt2}
entrez_search(db="sra",
              term="Tetrahymena thermophila[ORGN] AND 2013:2015[PDAT]",
              retmax=0)
```

Or include recent records for either _T. thermophila_ or it's close relative _T.
borealis_ (using parentheses to make ANDs and ORs explicit).


```{r, Tt3}
entrez_search(db="sra",
              term="(Tetrahymena thermophila[ORGN] OR Tetrahymena borealis[ORGN]) AND 2013:2015[PDAT]",
              retmax=0)
```

The set of search terms available varies between databases. You can get a list
of available terms or any given data base with `entrez_db_searchable()`

```{r, sra_searchable}
entrez_db_searchable("sra")
```

### Using the Filter field

"Filter" is a special field that, as the names suggests, allows you to limit 
records returned by a search to set of filtering criteria. There is no programmatic 
way to find the particular terms that can be used with the Filter field. 
However, the NCBI's website provides an "advanced search" tool for some 
databases that can be used to discover these terms. 


For example, to find the list of possible to find all of the terms that can be
used to filter searches to the nucleotide database using the 
[advanced search for that databse](https://www.ncbi.nlm.nih.gov/nuccore/advanced).
On that page selecting "Filter" from the first drop-down box then clicking 
"Show index list" will allow the user to scroll through possible filtering
terms.

###Precise queries using MeSH terms

In addition to the search terms described above, the NCBI allows searches using
[Medical Subject Heading (MeSH)](https://www.ncbi.nlm.nih.gov/mesh) terms. These
terms create a 'controlled vocabulary',  and allow users to make very finely
controlled queries of databases.

For instance, if you were interested in reviewing studies on how a class of
anti-malarial drugs called Folic Acid Antagonists work against _Plasmodium vivax_ (a
particular species of malarial parasite), you could use this search:

```{r, mesh}
entrez_search(db   = "pubmed",
              term = "(vivax malaria[MeSH]) AND (folic acid antagonists[MeSH])")
```

The complete set of MeSH terms is available as a database from the NCBI. That
means it is possible to download detailed information about each term and find 
the ways in which terms relate to each other using `rentrez`. You can search 
for specific terms with `entrez_search(db="mesh", term =...)` and learn about the
results of your search using the tools described below.

### Advanced counting

As you can see above, the  object returned by `entrez_search()` includes the
number of records matching a given search. This means you can learn a little
about the composition of, or trends in, the records stored in the NCBI's
databases using only the search utility. For instance, let's track the rise of
the scientific buzzword "connectome" in PubMed, programmatically creating 
search terms for the `PDAT` field:

```{r, connectome, fig.width=5, fig.height=4, fig.align='center'}
search_year <- function(year, term){
    query <- paste(term, "AND (", year, "[PDAT])")
    entrez_search(db="pubmed", term=query, retmax=0)$count
}

year <- 2008:2014
papers <- sapply(year, search_year, term="Connectome", USE.NAMES=FALSE)

plot(year, papers, type='b', main="The Rise of the Connectome")
```

## Finding cross-references : `entrez_link()`:


One of the strengths of the NCBI databases is the degree to which records of one
type are connected to  other records within the NCBI or to external data
sources. The function `entrez_link()` allows users to discover these links
between records.

### My god, it's full of links

To get an idea of the degree to which records in the NCBI are cross-linked we
can find all NCBI data associated with a single gene (in this case the 
Amyloid Beta Precursor gene, the product of which is associated with the 
plaques that form in the brains of  Alzheimer's Disease patients).

The function `entrez_link()` can be used to find cross-referenced records. In
the most basic case we need to provide an ID (`id`), the database from which this
ID comes (`dbfrom`) and the name of a database in which to find linked records (`db`).
If we set this last argument to 'all' we can find links in multiple databases:

```{r elink0}
all_the_links <- entrez_link(dbfrom='gene', id=351, db='all')
all_the_links
```
Just as with `entrez_search` the returned object behaves like a list, and we can
learn a little about its contents by printing it. In the case, all of the
information is in `links` (and there's a lot of them!):


```{r elink_link}
all_the_links$links
```
The names of the list elements are in the format `[source_database]_[linked_database]` 
and the elements themselves contain a vector of linked-IDs. So, if we want to
find open access publications associated with this gene we could get linked records
in PubMed Central:

```{r, elink_pmc}
all_the_links$links$gene_pmc[1:10]
```

Or if were interested in this gene's role in diseases we could find links to clinVar:

```{r, elink_omim}
all_the_links$links$gene_clinvar

```

### Narrowing our focus

If we know beforehand what sort of links we'd like to find , we can
to use the `db` argument to narrow the focus of a call to `entrez_link`.

For instance, say we are interested in knowing about all of the 
RNA transcripts associated with the Amyloid Beta Precursor gene in humans. 
Transcript sequences are stored in the nucleotide database (referred
to as `nuccore` in EUtils), so to find transcripts associated with a given gene
we need to set `dbfrom=gene` and `db=nuccore`.

```{r, elink1}
nuc_links <- entrez_link(dbfrom='gene', id=351, db='nuccore')
nuc_links
nuc_links$links
```
The object we get back contains links to the nucleotide database generally, but
also to special subsets of that database like [refseq](https://www.ncbi.nlm.nih.gov/refseq/). 
We can take advantage of this narrower set of links to find IDs that match unique
transcripts from our gene of interest.

```{r, elinik_refseqs}
nuc_links$links$gene_nuccore_refseqrna
```
We can use these ids in calls to `entrez_fetch()` or `entrez_summary()` to learn
more about the transcripts they represent. 

### External links

In addition to finding data within the NCBI, `entrez_link` can turn up
connections to external databases. Perhaps the most interesting example is
finding links to the full text of papers in PubMed. For example, when I wrote
this document the first paper linked to Amyloid Beta Precursor  had a unique ID of
`25500142`. We can find links to the full text of that paper with `entrez_link` 
by setting the `cmd` argument to 'llinks':

```{r, outlinks}
paper_links <- entrez_link(dbfrom="pubmed", id=25500142, cmd="llinks")
paper_links
```

Each element of the `linkouts` object contains information about an external
source of data on this paper:

```{r, urls}
paper_links$linkouts
```

Each of those linkout objects contains quite a lot of information, but the URL 
is probably the most useful. For that reason, `rentrez` provides the
function `linkout_urls` to make extracting just the URL simple:

```{r just_urls}
linkout_urls(paper_links)
```

The full list of options for the `cmd` argument are given in in-line
documentation (`?entrez_link`). If you are interested in finding full text
records for a large number of articles checkout the package 
[fulltext](https://github.com/ropensci/fulltext) which makes use of multiple
sources (including the NCBI) to discover the full text articles.

### Using more than one ID

It is possible to pass more than one ID to `entrez_link()`. By default, doing so
will give you a single elink object containing the complete set of links for 
_all_ of the IDs that you specified. So, if you were looking for protein IDs
related to specific genes you could do:

```{r, multi_default}
all_links_together  <- entrez_link(db="protein", dbfrom="gene", id=c("93100", "223646"))
all_links_together
all_links_together$links$gene_protein
```

Although this behaviour might sometimes be useful, it means we've lost track of
which `protein` ID is linked to which `gene` ID. To retain that information we
can set `by_id` to `TRUE`. This gives us a list of elink objects, each once
containing links from a single `gene` ID:

```{r, multi_byid}
all_links_sep  <- entrez_link(db="protein", dbfrom="gene", id=c("93100", "223646"), by_id=TRUE)
all_links_sep
lapply(all_links_sep, function(x) x$links$gene_protein)
```

## Getting summary data: `entrez_summary()`

Having found the unique IDs for some records via `entrez_search` or `entrez_link()`, you are
probably going to want to learn something about them. The `Eutils` API has two
ways to get information about a record. `entrez_fetch()` returns 'full' records
in varying formats and `entrez_summary()` returns less information about each
record, but in relatively simple format. Very often the summary records have the information
you are after, so `rentrez` provides functions to parse and summarise summary
records.


### The summary record

`entrez_summary()` takes a vector of unique IDs for the samples you want to get
summary information from. Let's start by finding out something about the paper
describing [Taxize](https://github.com/ropensci/taxize), using its PubMed ID:


```{r, Summ_1}
taxize_summ <- entrez_summary(db="pubmed", id=24555091)
taxize_summ
```

Once again, the object returned by `entrez_summary` behaves like a list, so you can extract
elements using `$`. For instance, we could convert our PubMed ID to another
article identifier...

```{r, Summ_2}
taxize_summ$articleids
```
...or see how many times the article has been cited in PubMed Central papers

```{r, Summ_3}
taxize_summ$pmcrefcount
```

### Dealing with many records

If you give `entrez_summary()` a vector with more than one ID you'll get a
list of summary records back. Let's get those _Plasmodium vivax_ papers we found
in the `entrez_search()` section back, and fetch some summary data on each paper:

```{r, multi_summ}
vivax_search <- entrez_search(db = "pubmed",
                              term = "(vivax malaria[MeSH]) AND (folic acid antagonists[MeSH])")
multi_summs <- entrez_summary(db="pubmed", id=vivax_search$ids)
```

`rentrez` provides a helper function, `extract_from_esummary()` that takes one
or more elements from every summary record in one of these lists. Here it is
working with one...

```{r, multi_summ2}
extract_from_esummary(multi_summs, "fulljournalname")
```
... and several elements:

```{r, multi_summ3}
date_and_cite <- extract_from_esummary(multi_summs, c("pubdate", "pmcrefcount",  "title"))
knitr::kable(head(t(date_and_cite)), row.names=FALSE)
```

## Fetching full records: `entrez_fetch()`

As useful as the summary records are, sometimes they just don't have the
information that you need. If you want a complete representation of a record you
can use `entrez_fetch`, using the argument `rettype` to specify the format you'd
like the record in.

### Fetch DNA sequences in fasta format

Let's extend the example given in the `entrez_link()` section about finding
transcript for a given gene. This time we will fetch cDNA sequences of those
transcripts.We can start by repeating the steps in the earlier example
to get nucleotide IDs for refseq transcripts of two genes:

```{r, transcript_ids}
gene_ids <- c(351, 11647)
linked_seq_ids <- entrez_link(dbfrom="gene", id=gene_ids, db="nuccore")
linked_transripts <- linked_seq_ids$links$gene_nuccore_refseqrna
head(linked_transripts)
```

Now we can get our sequences with `entrez_fetch`, setting `rettype` to "fasta"
(the list of formats available for [each database is give in this table](https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/)):

```{r fetch_fasta}
all_recs <- entrez_fetch(db="nuccore", id=linked_transripts, rettype="fasta")
class(all_recs)
nchar(all_recs)
```

Congratulations, now you have a really huge character vector! Rather than
printing all those thousands of bases we can take a peak at the top of the file:

```{r, peak}
cat(strwrap(substr(all_recs, 1, 500)), sep="\n")
```

If we wanted to use these sequences in some other application we could write them 
to file:

```r
write(all_recs, file="my_transcripts.fasta")
```

Alternatively, if you want to use them  within an R session  
we could write them to a  temporary file then read that. In this case I'm using `read.dna()` from the
pylogenetics package ape (but not executing the code block in this vignette, so
you don't have to install that package):

```r
temp <- tempfile()
write(all_recs, temp)
parsed_recs <- ape::read.dna(all_recs, temp)
```

### Fetch a parsed XML document

Most of the NCBI's databases can return records in XML format. In additional to
downloading the text-representation of these files, `entrez_fetch()` can return
objects parsed by the `XML` package. As an example, we can check out the  Taxonomy
database's record for (did I mention they are amazing....) _Tetrahymena
thermophila_, specifying we want the result to be parsed by setting
`parsed=TRUE`:

```{r, Tt_tax}
Tt <- entrez_search(db="taxonomy", term="(Tetrahymena thermophila[ORGN]) AND Species[RANK]")
tax_rec <- entrez_fetch(db="taxonomy", id=Tt$ids, rettype="xml", parsed=TRUE)
class(tax_rec)
```

The package XML (which you have if you have installed `rentrez`) provides
functions to get information from these files. For relatively simple records
like this one you can use `XML::xmlToList`:

```{r, Tt_list}
tax_list <- XML::xmlToList(tax_rec)
tax_list$Taxon$GeneticCode
```

For more complex records, which generate deeply-nested lists, you can use
[XPath expressions](https://en.wikipedia.org/wiki/XPath) along with the function 
`XML::xpathSApply` or the extraction operatord `[` and `[[` to extract specific 
parts of the file. For instance, we can get the scientific name of each taxon 
in _T. thermophila_'s lineage by specifying a path through the XML

```{r, Tt_path}
tt_lineage <- tax_rec["//LineageEx/Taxon/ScientificName"]
tt_lineage[1:4]
```

As the name suggests, `XML::xpathSApply()` is a counterpart of base R's
`sapply`, and can be used to apply a function to
nodes in an XML object. A particularly useful function to apply is `XML::xmlValue`, 
which returns the content of the node:

```{r, Tt_apply}
XML::xpathSApply(tax_rec, "//LineageEx/Taxon/ScientificName", XML::xmlValue)
```
There are a few more complex examples of using `XPath` [on the rentrez wiki](https://github.com/ropensci/rentrez/wiki)

<a name="web_history"></a>

## Using NCBI's Web History features

When you are dealing with very large queries it can be time  consuming to pass
long vectors of unique IDs to and from the NCBI. To avoid this problem, the NCBI
provides a feature called "web history" which allows users to store IDs on the
NCBI servers then refer to them in future calls.

###Post a set of IDs to the NCBI for later use: `entrez_post()`

If you have a list of many NCBI IDs that you want to use later on, you can post
them to the NCBI's severs. In order to provide a brief example, I'm going to post just one
ID, the `omim` identifier for asthma:

```{r, asthma}
upload <- entrez_post(db="omim", id=600807)
upload
```
The NCBI sends you back some information you can use to refer to the posted IDs. 
In `rentrez`, that information is represented as a `web_history` object. 

Note that if you have a very long list of IDs you may receive a 414 error when
you try to upload them. If you have such a list (and they come from an external
sources rather than a search that can be save to a `web_history` object), you
may have to 'chunk' the IDs into smaller sets that can processed. 

### Get a `web_history` object from `entrez_search` or `entrez_link()`

In addition to directly uploading IDs to the NCBI, you can use the web history
features with `entrez_search` and `entrez_link`. For instance, imagine you wanted to
find all of the sequences of the widely-studied gene COI from all snails
(which are members of the taxonomic group Gastropoda):

```{r, snail_search}
entrez_search(db="nuccore", term="COI[Gene] AND Gastropoda[ORGN]")
```

That's a lot of sequences! If you really wanted to download all of these it
would be a good idea to save all those IDs to the server by setting
`use_history` to `TRUE` (note you now get a `web_history` object along with your
normal search result):

```{r, snail_history}
snail_coi <- entrez_search(db="nuccore", term="COI[Gene] AND Gastropoda[ORGN]", use_history=TRUE)
snail_coi
snail_coi$web_history
```

Similarity, `entrez_link()` can return `web_history` objects by using the `cmd`
`neighbor_history`. Let's find genetic variants (from the clinvar database)
associated with asthma (using the same OMIM ID we identified earlier):


```{r, asthma_links}
asthma_clinvar <- entrez_link(dbfrom="omim", db="clinvar", cmd="neighbor_history", id=600807)
asthma_clinvar$web_histories
```

As you can see, instead of returning lists of IDs for each linked database (as
it would be default), `entrez_link()` now  returns a list of web_histories.

### Use a `web_history` object

Once you have those IDs stored on the NCBI's servers, you are going to want to
do something with them. The functions `entrez_fetch()` `entrez_summary()` and
`entrez_link()` can all use `web_history` objects in exactly the same way they
use IDs. 

So, we could repeat the last example (finding variants linked to asthma), but this
time using the ID we uploaded earlier

```{r, asthma_links_upload}
asthma_variants <- entrez_link(dbfrom="omim", db="clinvar", cmd="neighbor_history", web_history=upload)
asthma_variants
```

... if we want to get some genetic information about these variants we need to
map our clinvar IDs to SNP IDs:
                           

```{r, links}
snp_links <- entrez_link(dbfrom="clinvar", db="snp", 
                         web_history=asthma_variants$web_histories$omim_clinvar,
                         cmd="neighbor_history")
snp_summ <- entrez_summary(db="snp", web_history=snp_links$web_histories$clinvar_snp)
knitr::kable(extract_from_esummary(snp_summ, c("chr", "fxn_class", "global_maf")))
```

If you really wanted to you could also use `web_history` objects to download all those thousands of COI sequences.
When downloading large sets of data, it is a good idea to take advantage of the
arguments `retmax` and `restart` to split the request up into smaller chunks.
For instance, we could get the first 200 sequences in 50-sequence chunks:

(note: this code block is not executed as part of the vignette to save time and bandwidth):


```r
for( seq_start in seq(1,200,50)){
    recs <- entrez_fetch(db="nuccore", web_history=snail_coi$web_history,
                         rettype="fasta", retmax=50, retstart=seq_start)
    cat(recs, file="snail_coi.fasta", append=TRUE)
    cat(seq_start+49, "sequences downloaded\r")
}
```

## Rate-limiting and API Keys

By default, the NCBI limits users to making only 3 requests per second (and
`rentrez` enforces that limit). Users who register for an "API key" are able to
make up to ten requests per second. Getting one of these keys is simple, you
just need to [register for "my ncbi"
account](https://www.ncbi.nlm.nih.gov/account/) then click on a button in the
[account settings page](https://www.ncbi.nlm.nih.gov/account/settings/).


Once you have an API key, rentrez will allow you to take advantage of it. For 
one-off cases, this is as simple as adding the `api_key` argument to given
function call. (Note these examples are not executed as part of this document, 
as the  API key used here not a real one).

```r
entrez_link(db="protein", dbfrom="gene", id=93100, api_key ="ABCD123")
```

It most cases you will want to use your API for each of several calls to the
NCBI. `rentrez` makes this easy by allowing you to set an environment variable
,`ENTREZ_KEY`. Once this value is set to your key `rentrez` will use it for all
requests to the NCBI. To set the value for a single R session you can use the
function `set_entrez_key()`. Here we set the value and confirm it is available.

```{r, set_key}
set_entrez_key("ABCD123")
Sys.getenv("ENTREZ_KEY")
```

If you use `rentrez` often you should edit your `.Renviron` file (see `r
help(Startup)` for description of this file)  to include your key. Doing so will
mean all requests you send will take advantage of your API key.

```ini
ENTREZ_KEY=ABCD123
```

As long as an API key is set by one of these methods, `rentrez` will allow you
to make up to ten requests per second.

### Slowing rentrez down when you hit the rate-limit

`rentrez` won't let you _send_ requests to the NCBI at a rate higher than the
rate-limit, but it is sometimes possible that they will _arrive_ too close
together an produce errors. If you are using `rentrez` functions in a for loop
and find rate-limiting errors are occuring, you may consider adding a call to
`Sys.sleep(0.1)` before each message sent to the NCBI. This will ensure you stay
beloe the rate limit.


## What next ?

This tutorial has introduced you to the core functions of `rentrez`, there are
almost limitless ways that you could put them together. [Check out the wiki](https://github.com/ropensci/rentrez/wiki)
for more specific examples, and be sure to read the inline-documentation for
each function. If you run into problem with rentrez, or just need help with the
package and `Eutils` please contact us by opening an issue at the [github
repository](https://github.com/ropensci/rentrez/issues)
 


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_summary.r
\name{entrez_summary}
\alias{entrez_summary}
\title{Get summaries of objects in NCBI datasets from a unique ID}
\usage{
entrez_summary(
  db,
  id = NULL,
  web_history = NULL,
  version = c("2.0", "1.0"),
  always_return_list = FALSE,
  retmode = NULL,
  config = NULL,
  ...
)
}
\arguments{
\item{db}{character Name of the database to search for}

\item{id}{vector with unique ID(s) for records in database \code{db}.
In the case of sequence databases these IDs can take form of an
NCBI accession followed by a version number (eg AF123456.1 or AF123456.2)}

\item{web_history}{A web_history object}

\item{version}{either 1.0 or 2.0 see above for description}

\item{always_return_list}{logical, return a list  of esummary objects even
when only one ID is provided (see description for a note about this option)}

\item{retmode}{either "xml" or "json". By default, xml will be used for
version 1.0 records, json for version 2.0.}

\item{config}{vector configuration options passed to \code{httr::GET}}

\item{\dots}{character Additional terms to add to the request, see NCBI
documentation linked to in references for a complete list}
}
\value{
A list of esummary records (if multiple IDs are passed and
always_return_list if FALSE) or a single record.

file XMLInternalDocument xml file containing the entire record
returned by the NCBI.
}
\description{
The NCBI offer two distinct formats for summary documents.
Version 1.0 is a relatively limited summary of a database record based on a 
shared Document Type Definition. Version 1.0 summaries are only available as
XML and are not available for some newer databases
Version 2.0 summaries generally contain more information about a given
record, but each database has its own distinct format. 2.0 summaries are 
available for records in all databases and as JSON and XML files. 
As of version 0.4, rentrez fetches version 2.0 summaries by default and
uses JSON as the exchange format (as JSON object can be more easily converted
into native R types). Existing scripts which relied on the structure and
naming of the "Version 1.0" summary files can be updated by setting the new
\code{version} argument to "1.0".
}
\details{
By default, entrez_summary returns a single record when only one ID is
passed and a list of such records when multiple IDs are passed. This can lead
to unexpected behaviour when the results of a variable number of IDs (perhaps the
result of \code{entrez_search}) are processed with an apply family function
or in a for-loop. If you use this function as part of a function or script that
generates a variably-sized vector of IDs setting \code{always_return_list} to 
\code{TRUE} will avoid these problems. The function
\code{extract_from_esummary} is  provided for the specific case of extracting
named elements from a list of esummary objects, and is designed to work on
single objects as well as lists.
}
\examples{
\dontrun{
 pop_ids = c("307082412", "307075396", "307075338", "307075274")
 pop_summ <- entrez_summary(db="popset", id=pop_ids)
 extract_from_esummary(pop_summ, "title")
 
 # clinvar example
 res <- entrez_search(db = "clinvar", term = "BRCA1", retmax=10)
 cv <- entrez_summary(db="clinvar", id=res$ids)
 cv
 extract_from_esummary(cv, "title", simplify=FALSE)
 extract_from_esummary(cv, "trait_set")[1:2] 
 extract_from_esummary(cv, "gene_sort") 
}
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESummary_}
}
\seealso{
\code{\link[httr]{config}} for available configs

\code{\link{extract_from_esummary}} which can be used to extract
elements from a list of esummary records
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_info.r
\name{entrez_db_links}
\alias{entrez_db_links}
\title{List available links for records from a given NCBI database}
\usage{
entrez_db_links(db, config = NULL)
}
\arguments{
\item{db}{character, name of database to search}

\item{config}{config vector passed to \code{httr::GET}}
}
\value{
An eInfoLink object (sub-classed from list) summarizing linked-databases.
Can be coerced to a data-frame with \code{as.data.frame}. Printing the object
the name of each element (which is the correct name for \code{entrez_link},
and can be used to get (a little) more information about each linked database
(see example below).
}
\description{
For a given database, fetch a list of other databases that contain
cross-referenced records. The names of these records can be used as the
\code{db} argument in \code{\link{entrez_link}}
}
\examples{
\dontrun{
taxid <- entrez_search(db="taxonomy", term="Osmeriformes")$ids
tax_links <- entrez_db_links("taxonomy")
tax_links
entrez_link(dbfrom="taxonomy", db="pmc", id=taxid)

sra_links <- entrez_db_links("sra")
as.data.frame(sra_links)
}
}
\seealso{
\code{\link{entrez_link}}

Other einfo: 
\code{\link{entrez_db_searchable}()},
\code{\link{entrez_db_summary}()},
\code{\link{entrez_dbs}()},
\code{\link{entrez_info}()}
}
\concept{einfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_info.r
\name{entrez_db_searchable}
\alias{entrez_db_searchable}
\title{List available search fields for a given database}
\usage{
entrez_db_searchable(db, config = NULL)
}
\arguments{
\item{db}{character, name of database to get search field from}

\item{config}{config vector passed to \code{httr::GET}}
}
\value{
An eInfoSearch object (subclassed from list) summarizing linked-databases. 
Can be coerced to a data-frame with \code{as.data.frame}. Printing the object
shows only the names of each available search field.
}
\description{
Fetch a list of search fields that can be used with a given database. Fields
can be used as part of the \code{term} argument to \code{\link{entrez_search}}
}
\examples{
\dontrun{
pmc_fields <- entrez_db_searchable("pmc")
pmc_fields[["AFFL"]]
entrez_search(db="pmc", term="Otago[AFFL]", retmax=0)
entrez_search(db="pmc", term="Auckland[AFFL]", retmax=0)

sra_fields <- entrez_db_searchable("sra")
as.data.frame(sra_fields)
}
}
\seealso{
\code{\link{entrez_search}}

Other einfo: 
\code{\link{entrez_db_links}()},
\code{\link{entrez_db_summary}()},
\code{\link{entrez_dbs}()},
\code{\link{entrez_info}()}
}
\concept{einfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_pubmed_xml.r
\name{parse_pubmed_xml}
\alias{parse_pubmed_xml}
\title{Summarize an XML record from pubmed.}
\usage{
parse_pubmed_xml(record)
}
\arguments{
\item{record}{Either and XMLInternalDocument or character the record to be
parsed ( expected to come from \code{\link{entrez_fetch}})}
}
\value{
Either a single pubmed_record object, or a list of several
}
\description{
Note: this function assumes all records are of the type "PubmedArticle"
and will return an empty record for any other type (including books).
}
\examples{
\donttest{
hox_paper <- entrez_search(db="pubmed", term="10.1038/nature08789[doi]")
hox_rel <- entrez_link(db="pubmed", dbfrom="pubmed", id=hox_paper$ids)
recs <- entrez_fetch(db="pubmed", 
                       id=hox_rel$links$pubmed_pubmed[1:3], 
                       rettype="xml")
parse_pubmed_xml(recs)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_citmatch.r
\name{entrez_citmatch}
\alias{entrez_citmatch}
\title{Fetch pubmed ids matching specially formatted citation strings}
\usage{
entrez_citmatch(bdata, db = "pubmed", retmode = "xml", config = NULL)
}
\arguments{
\item{bdata}{character, containing citation data. 
Each citation must be represented in a pipe-delimited format
 journal_title|year|volume|first_page|author_name|your_key|
The final field "your_key" is arbitrary, and can used as you see
fit. Fields can be left empty, but be sure to keep 6 pipes.}

\item{db}{character, the database to search. Defaults to pubmed,
the only database currently available}

\item{retmode}{character, file format to retrieve. Defaults to xml, as 
per the API documentation, though note the API only returns plain text}

\item{config}{vector configuration options passed to httr::GET}
}
\value{
A character vector containing PMIDs
}
\description{
Fetch pubmed ids matching specially formatted citation strings
}
\examples{
\dontrun{
ex_cites <- c("proc natl acad sci u s a|1991|88|3248|mann bj|test1|",
              "science|1987|235|182|palmenberg ac|test2|")
entrez_citmatch(ex_cites)
}
}
\seealso{
\code{\link[httr]{config}} for available configs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_fetch.r
\name{entrez_fetch}
\alias{entrez_fetch}
\title{Download data from NCBI databases}
\usage{
entrez_fetch(
  db,
  id = NULL,
  web_history = NULL,
  rettype,
  retmode = "",
  parsed = FALSE,
  config = NULL,
  ...
)
}
\arguments{
\item{db}{character, name of the database to use}

\item{id}{vector (numeric or character), unique ID(s) for records in database
\code{db}. In the case of sequence databases these IDs can take form of an
NCBI accession followed by a version number (eg AF123456.1 or AF123456.2).}

\item{web_history, }{a web_history object}

\item{rettype}{character, format in which to get data (eg, fasta, xml...)}

\item{retmode}{character, mode in which to receive data, defaults to an empty
string (corresponding to the default mode for rettype).}

\item{parsed}{boolean should entrez_fetch attempt to parse the resulting 
file. Only works with xml records (including those with rettypes other than
"xml") at present}

\item{config}{vector, httr configuration options passed to httr::GET}

\item{\dots}{character, additional terms to add to the request, see NCBI
documentation linked to in references for a complete list}
}
\value{
character string containing the file created

XMLInternalDocument a parsed XML document if parsed=TRUE and
rettype is a flavour of XML.
}
\description{
Pass unique identifiers to an NCBI database and receive data files in a
variety of formats.
A set of unique identifiers mustbe specified with either the \code{db}
argument (which directly specifies the IDs as a numeric or character vector)
or a \code{web_history} object as returned by 
\code{\link{entrez_link}}, \code{\link{entrez_search}} or 
\code{\link{entrez_post}}.
}
\details{
The format for returned records is set by that arguments \code{rettype} (for
a particular format) and \code{retmode} for a general format (JSON, XML text
etc). See  \href{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/}{Table 1} 
in the linked reference for the set of 
formats available for each database. In particular, note that sequence
databases (nuccore, protein and their relatives) use specific format names
(eg "native", "ipg") for different flavours of xml.

For the most part, this function returns a character vector containing the 
fetched records. For XML records (including 'native', 'ipg', 'gbc' sequence
records), setting \code{parsed} to \code{TRUE} will return an
\code{XMLInternalDocument},
}
\examples{
\dontrun{
katipo <- "Latrodectus katipo[Organism]"
katipo_search <- entrez_search(db="nuccore", term=katipo)
kaitpo_seqs <- entrez_fetch(db="nuccore", id=katipo_search$ids, rettype="fasta")
#xml
kaitpo_seqs <- entrez_fetch(db="nuccore", id=katipo_search$ids, rettype="native")
}
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_}
}
\seealso{
\code{\link[httr]{config}} for available '\code{httr}` configs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_info.r
\name{entrez_info}
\alias{entrez_info}
\title{Get information about EUtils databases}
\usage{
entrez_info(db = NULL, config = NULL)
}
\arguments{
\item{db}{character database about which to retrieve information (optional)}

\item{config}{config vector passed on to \code{httr::GET}}
}
\value{
XMLInternalDocument with information describing either all the
databases available in Eutils (if db is not set) or one particular database
(set by 'db')
}
\description{
Gather information about EUtils generally, or a given Eutils database.
Note: The most common uses-cases for the einfo util are finding the list of
search fields available for a given database or the other NCBI databases to
which records in a given database might be linked. Both these use cases
are implemented in higher-level functions that return just this information
(\code{entrez_db_searchable} and \code{entrez_db_links} respectively).
Consequently most users will not have a reason to use this function (though
it is exported by \code{rentrez} for the sake of completeness.
}
\examples{
\dontrun{
all_the_data <- entrez_info()
XML::xpathSApply(all_the_data, "//DbName", xmlValue)
entrez_dbs()
}
}
\seealso{
\code{\link[httr]{config}} for available httr configurations

Other einfo: 
\code{\link{entrez_db_links}()},
\code{\link{entrez_db_searchable}()},
\code{\link{entrez_db_summary}()},
\code{\link{entrez_dbs}()}
}
\concept{einfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_link.r
\name{linkout_urls}
\alias{linkout_urls}
\title{Extract URLs from an elink object}
\usage{
linkout_urls(elink)
}
\arguments{
\item{elink}{elink object (returned by entrez_link) containing Urls}
}
\value{
list of character vectors, one per ID each containing of URLs for that
ID.
}
\description{
Extract URLs from an elink object
}
\seealso{
entrez_link
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_summary.r
\name{extract_from_esummary}
\alias{extract_from_esummary}
\title{Extract elements from a list of esummary records}
\usage{
extract_from_esummary(esummaries, elements, simplify = TRUE)
}
\arguments{
\item{esummaries}{Either an esummary or an esummary_list (as returned by
entrez_summary).}

\item{elements}{the names of the element to extract}

\item{simplify}{logical, if possible return a vector}
}
\value{
List or vector containing requested elements
}
\description{
Extract elements from a list of esummary records
}
\seealso{
\code{\link{entrez_summary}} for examples of this function in action.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_info.r
\name{entrez_dbs}
\alias{entrez_dbs}
\title{List databases available from the NCBI}
\usage{
entrez_dbs(config = NULL)
}
\arguments{
\item{config}{config vector passed to \code{httr::GET}}
}
\value{
character vector listing available dbs
}
\description{
Retrieves the names of  databases available through the EUtils API
}
\examples{
\dontrun{
entrez_dbs()
}
}
\seealso{
Other einfo: 
\code{\link{entrez_db_links}()},
\code{\link{entrez_db_searchable}()},
\code{\link{entrez_db_summary}()},
\code{\link{entrez_info}()}
}
\concept{einfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_search.r
\name{entrez_search}
\alias{entrez_search}
\title{Search the NCBI databases using EUtils}
\usage{
entrez_search(
  db,
  term,
  config = NULL,
  retmode = "xml",
  use_history = FALSE,
  ...
)
}
\arguments{
\item{db}{character, name of the database to search for.}

\item{term}{character, the search term. The syntax used in making these
searches is described in the Details of this help message, the package
vignette and reference given below.}

\item{config}{vector configuration options passed to httr::GET}

\item{retmode}{character, one of json (default) or xml. This will make no
difference in most cases.}

\item{use_history}{logical. If TRUE return a web_history object for use in 
later calls to the NCBI}

\item{\dots}{character, additional terms to add to the request, see NCBI
documentation linked to in references for a complete list}
}
\value{
ids integer Unique IDS returned by the search

count integer Total number of hits for the search

retmax integer Maximum number of hits returned by the search

web_history A web_history object for use in subsequent calls to NCBI

QueryTranslation character, search term as the NCBI interpreted it

file either and XMLInternalDocument xml file resulting from search, parsed with
\code{\link[XML]{xmlTreeParse}} or, if \code{retmode} was set to json a list
resulting from the returned JSON file being parsed with
\code{\link[jsonlite]{fromJSON}}.
}
\description{
Search a given NCBI database with a particular query.
}
\details{
The NCBI uses a search term syntax where search terms can be associated with 
a specific search field with square brackets. So, for instance ``Homo[ORGN]''
denotes a search for Homo in the ``Organism'' field. The names and
definitions of these fields can be identified using
\code{\link{entrez_db_searchable}}.

Searches can make use of several fields by combining them via the boolean
operators AND, OR and NOT. So, using the search term``((Homo[ORGN] AND APP[GENE]) NOT
Review[PTYP])'' in PubMed would identify articles matching the gene APP in
humans, and exclude review articles. More examples of the use of these search
terms, and the more specific MeSH terms for precise searching, 
is given in the package vignette. \code{rentrez} handles special characters
and URL encoding (e.g. replacing spaces with plus signs) on the client side,
so there is no need to include these in search term

The\code{rentrez} tutorial provides some tips on how to make the most of 
searches to the NCBI. In particular, the sections on uses of the "Filter"
field and MeSH terms may in formulating precise searches.
}
\examples{
\dontrun{
   query <- "Gastropoda[Organism] AND COI[Gene]"
   web_env_search <- entrez_search(db="nuccore", query, use_history=TRUE)
   cookie <- web_env_search$WebEnv
   qk <- web_env_search$QueryKey 
   snail_coi <- entrez_fetch(db = "nuccore", WebEnv = cookie, query_key = qk,
                             file_format = "fasta", retmax = 10)
}
\dontrun{

fly_id <- entrez_search(db="taxonomy", term="Drosophila")
#Oh, right. There is a genus and a subgenus name Drosophila...
#how can we limit this search
(tax_fields <- entrez_db_searchable("taxonomy"))
#"RANK" loots promising
tax_fields$RANK
entrez_search(db="taxonomy", term="Drosophila & Genus[RANK]")
}
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_}
}
\seealso{
\code{\link[httr]{config}} for available httr configurations

\code{\link{entrez_db_searchable}} to get a set of search fields that
can be used in \code{term} for any database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_info.r
\name{entrez_db_summary}
\alias{entrez_db_summary}
\title{Retrieve summary information about an NCBI database}
\usage{
entrez_db_summary(db, config = NULL)
}
\arguments{
\item{db}{character, name of database to summaries}

\item{config}{config vector passed to \code{httr::GET}}
}
\value{
Character vector with the following data

DbName Name of database

Description Brief description of the database

Count Number of records contained in the database

MenuName Name in web-interface to EUtils

DbBuild Unique ID for current build of database

LastUpdate Date of most recent update to database
}
\description{
Retrieve summary information about an NCBI database
}
\examples{
\dontrun{
entrez_db_summary("pubmed")
}
}
\seealso{
Other einfo: 
\code{\link{entrez_db_links}()},
\code{\link{entrez_db_searchable}()},
\code{\link{entrez_dbs}()},
\code{\link{entrez_info}()}
}
\concept{einfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_keys.r
\name{set_entrez_key}
\alias{set_entrez_key}
\title{Set the ENTREZ_KEY variable to be used by all rentrez functions}
\usage{
set_entrez_key(key)
}
\arguments{
\item{key}{character. Value to set ENTREZ_KEY to (i.e. your API key).}
}
\value{
A logical of length one, TRUE is the value was set FALSE if not.
value is returned inside invisible(), i.e. it is not printed to screen 
when the function is called.
}
\description{
The NCBI allows users to access more records (10 per second) if they 
register for and use an API key. This function allows users to set this key
for all calls to rentrez functions during a particular R session. See the
vignette section "Using API keys" for a detailed description.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_post.r
\name{entrez_post}
\alias{entrez_post}
\title{Post IDs to Eutils for later use}
\usage{
entrez_post(db, id = NULL, web_history = NULL, config = NULL, ...)
}
\arguments{
\item{db}{character Name of the database from which the IDs were taken}

\item{id}{vector with unique ID(s) for records in database \code{db}.}

\item{web_history}{A web_history object. Can be used to add to additional
identifiers to an existing web environment on the NCBI}

\item{config}{vector of configuration options passed to httr::GET}

\item{\dots}{character Additional terms to add to the request, see NCBI
documentation linked to in references for a complete list}
}
\description{
Post IDs to Eutils for later use
}
\examples{
\dontrun{  
so_many_snails <- entrez_search(db="nuccore", 
                      "Gastropoda[Organism] AND COI[Gene]", retmax=200)
upload <- entrez_post(db="nuccore", id=so_many_snails$ids)
first <- entrez_fetch(db="nuccore", rettype="fasta", web_history=upload,
                      retmax=10)
second <- entrez_fetch(db="nuccore", file_format="fasta", web_history=upload,
                       retstart=10, retmax=10)
}
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EPost_}
}
\seealso{
\code{\link[httr]{config}} for available httr configurations
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_global_query.r
\name{entrez_global_query}
\alias{entrez_global_query}
\title{Find the number of records that match a given term across all NCBI Entrez databases}
\usage{
entrez_global_query(term, config = NULL, ...)
}
\arguments{
\item{term}{the search term to use}

\item{config}{vector configuration options passed to httr::GET}

\item{...}{additional arguments to add to the query}
}
\value{
a named vector with counts for each a database
}
\description{
Find the number of records that match a given term across all NCBI Entrez databases
}
\examples{
\dontrun{ 
NCBI_data_on_best_butterflies_ever <- entrez_global_query(term="Heliconius")
}
}
\seealso{
\code{\link[httr]{config}} for available configs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help.r
\docType{package}
\name{rentrez}
\alias{rentrez}
\alias{rentrez-package}
\title{rentrez}
\description{
rentrez provides functions to search for, discover and download data from
the NCBI's databases using their EUtils function.
}
\details{
Users are expected to know a little bit about the EUtils API, which is well
documented: \url{https://www.ncbi.nlm.nih.gov/books/NBK25500/}

The NCBI will ban IPs that don't use EUtils within their \href{https://www.ncbi.nlm.nih.gov/corehtml/query/static/eutils_help.html}{user guidelines}. In particular
/enumerated{
 /item  Don't send more than three request per second (rentrez enforces this limit)
 /item  If you plan on sending a sequence of more than ~100 requests, do so outside of peak times for the US
 /item  For large requests use the web history method (see examples for \code{\link{entrez_search}} or use \code{\link{entrez_post}} to upload IDs)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_link.r
\name{entrez_link}
\alias{entrez_link}
\title{Get links to datasets related to records from an NCBI database}
\usage{
entrez_link(
  dbfrom,
  web_history = NULL,
  id = NULL,
  db = NULL,
  cmd = "neighbor",
  by_id = FALSE,
  config = NULL,
  ...
)
}
\arguments{
\item{dbfrom}{character Name of database from which the Id(s) originate}

\item{web_history}{a web_history object}

\item{id}{vector with unique ID(s) for records in database \code{db}.}

\item{db}{character Name of the database to search for links (or use "all" to 
search all databases available for \code{db}. \code{entrez_db_links} allows you
to discover databases that might have linked information (see examples).}

\item{cmd}{link function to use. Allowed values include
\itemize{
  \item neighbor (default). Returns a set of IDs in \code{db} linked to the
  input IDs in \code{dbfrom}.
  \item neighbor_score. As `neighbor'', but additionally returns similarity scores.
  \item neighbor_history. As `neighbor', but returns web history objects.
  \item acheck. Returns a list of linked databases available from NCBI for a set of IDs.
  \item ncheck. Checks for the existence of links within a single database.
  \item lcheck. Checks for external (i.e. outside NCBI) links.
  \item llinks. Returns a list of external links for each ID, excluding links
  provided by libraries.
  \item llinkslib. As 'llinks' but additionally includes links provided by
  libraries.
  \item prlinks. As 'llinks' but returns only the primary external link for
  each ID.
}}

\item{by_id}{logical If FALSE (default) return a single 
\code{elink} objects containing links for all of the provided \code{id}s. 
Alternatively, if TRUE return a list of \code{elink} objects, one for each 
ID in \code{id}.}

\item{config}{vector configuration options passed to httr::GET}

\item{\dots}{character Additional terms to add to the request, see NCBI
documentation linked to in references for a complete list}
}
\value{
An elink object containing the data defined by the \code{cmd} argument
(if by_id=FALSE) or a list of such object (if by_id=TRUE).

file XMLInternalDocument xml file resulting from search, parsed with
\code{\link{xmlTreeParse}}
}
\description{
Discover records related to a set of unique identifiers from
an NCBI database. The object returned by this function depends on the value
set for the \code{cmd} argument. Printing the returned object lists the names
, and provides a brief description, of the elements included in the object.
}
\examples{
\dontrun{
 pubmed_search <- entrez_search(db = "pubmed", term ="10.1016/j.ympev.2010.07.013[doi]")
 linked_dbs <- entrez_db_links("pubmed")
 linked_dbs
 nucleotide_data <- entrez_link(dbfrom = "pubmed", id = pubmed_search$ids, db ="nuccore")
 #Sources for the full text of the paper 
 res <- entrez_link(dbfrom="pubmed", db="", cmd="llinks", id=pubmed_search$ids)
 linkout_urls(res)
}

}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ELink_}
}
\seealso{
\code{\link[httr]{config}} for available configs

\code{entrez_db_links}
}
