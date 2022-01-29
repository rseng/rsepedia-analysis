solrium
=======



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/solrium)](https://cranchecks.info/pkgs/solrium)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/solrium?color=2ED968)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/solrium)](https://cran.r-project.org/package=solrium)

**A general purpose R interface to [Solr](https://solr.apache.org/)**

Development is now following Solr v7 and greater - which introduced many changes, which means many functions here may not work with your Solr installation older than v7.

Be aware that currently some functions will only work in certain Solr modes, e.g, `collection_create()` won't work when you are not in Solrcloud mode. But, you should get an error message stating that you aren't.

Currently developing against Solr `v8.2.0`

## Package API and ways of using the package

The first thing to look at is `SolrClient` to instantiate a client connection
to your Solr instance. `ping` and `schema` are helpful functions to look
at after instantiating your client.

There are two ways to use `solrium`:

1. Call functions on the `SolrClient` object
2. Pass the `SolrClient` object to functions

For example, if we instantiate a client like `conn <- SolrClient$new()`, then
to use the first way we can do `conn$search(...)`, and the second way by doing
`solr_search(conn, ...)`. These two ways of using the package hopefully
make the package more user friendly for more people, those that prefer a more
object oriented approach, and those that prefer more of a functional approach.

**Collections**

Functions that start with `collection` work with Solr collections when in
cloud mode. Note that these functions won't work when in Solr standard mode

**Cores**

Functions that start with `core` work with Solr cores when in standard Solr
mode. Note that these functions won't work when in Solr cloud mode

**Documents**

The following functions work with documents in Solr

```
#>  - add
#>  - delete_by_id
#>  - delete_by_query
#>  - update_atomic_json
#>  - update_atomic_xml
#>  - update_csv
#>  - update_json
#>  - update_xml
```

**Search**

Search functions, including `solr_parse` for parsing results from different
functions appropriately

```
#>  - solr_all
#>  - solr_facet
#>  - solr_get
#>  - solr_group
#>  - solr_highlight
#>  - solr_mlt
#>  - solr_parse
#>  - solr_search
#>  - solr_stats
```


## Install

Stable version from CRAN


```r
install.packages("solrium")
```

Or development version from GitHub


```r
remotes::install_github("ropensci/solrium")
```


```r
library("solrium")
```

## Setup

Use `SolrClient$new()` to initialize your connection. These examples use a remote Solr server, but work on any local Solr server.


```r
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))
#> <Solr Client>
#>   host: api.plos.org
#>   path: search
#>   port: 
#>   scheme: http
#>   errors: simple
#>   proxy:
```

You can also set whether you want simple or detailed error messages (via `errors`), and whether you want URLs used in each function call or not (via `verbose`), and your proxy settings (via `proxy`) if needed. For example:


```r
SolrClient$new(errors = "complete")
```

Your settings are printed in the print method for the connection object


```r
cli
#> <Solr Client>
#>   host: api.plos.org
#>   path: search
#>   port: 
#>   scheme: http
#>   errors: simple
#>   proxy:
```

For local Solr server setup:

```
bin/solr start -e cloud -noprompt
bin/post -c gettingstarted example/exampledocs/*.xml
```


## Search


```r
(res <- cli$search(params = list(q='*:*', rows=2, fl='id')))
#> # A tibble: 2 x 1
#>   id                                   
#>   <chr>                                
#> 1 10.1371/journal.pbio.1000146/title   
#> 2 10.1371/journal.pbio.1000146/abstract
```

And you can get search metadata from the attributes:


```r
attributes(res)
#> $names
#> [1] "id"
#> 
#> $row.names
#> [1] 1 2
#> 
#> $class
#> [1] "tbl_df"     "tbl"        "data.frame"
#> 
#> $numFound
#> [1] 2542432
#> 
#> $start
#> [1] 0
```


### Search grouped data

Most recent publication by journal


```r
cli$group(params = list(q='*:*', group.field='journal', rows=5, group.limit=1,
                        group.sort='publication_date desc',
                        fl='publication_date, score'))
#>                   groupValue numFound start     publication_date score
#> 1               plos biology    45430     0 2021-05-18T00:00:00Z     1
#> 2 plos computational biology    68336     0 2021-05-18T00:00:00Z     1
#> 3              plos genetics    78511     0 2021-05-18T00:00:00Z     1
#> 4              plos medicine    33148     0 2021-05-18T00:00:00Z     1
#> 5                       none    57571     0 2012-10-23T00:00:00Z     1
```

First publication by journal


```r
cli$group(params = list(q = '*:*', group.field = 'journal', group.limit = 1,
                        group.sort = 'publication_date asc',
                        fl = c('publication_date', 'score'),
                        fq = "publication_date:[1900-01-01T00:00:00Z TO *]"))
#>                          groupValue numFound start     publication_date score
#> 1                      plos biology    45430     0 2003-08-18T00:00:00Z     1
#> 2        plos computational biology    68336     0 2005-06-24T00:00:00Z     1
#> 3                     plos genetics    78511     0 2005-06-17T00:00:00Z     1
#> 4                     plos medicine    33148     0 2004-09-07T00:00:00Z     1
#> 5                              none    57571     0 2005-08-23T00:00:00Z     1
#> 6              plos clinical trials      521     0 2006-04-21T00:00:00Z     1
#> 7  plos neglected tropical diseases    75150     0 2007-08-30T00:00:00Z     1
#> 8                    plos pathogens    73595     0 2005-07-22T00:00:00Z     1
#> 9                          plos one  2110161     0 2006-12-20T00:00:00Z     1
#> 10                     plos medicin        9     0 2012-04-17T00:00:00Z     1
```

Search group query : Last 3 publications of 2013.


```r
gq <- 'publication_date:[2013-01-01T00:00:00Z TO 2013-12-31T00:00:00Z]'
cli$group(
  params = list(q='*:*', group.query = gq,
                group.limit = 3, group.sort = 'publication_date desc',
                fl = 'publication_date'))
#>   numFound start     publication_date
#> 1   307446     0 2013-12-31T00:00:00Z
#> 2   307446     0 2013-12-31T00:00:00Z
#> 3   307446     0 2013-12-31T00:00:00Z
```

Search group with format simple


```r
cli$group(params = list(q='*:*', group.field='journal', rows=5,
                        group.limit=3, group.sort='publication_date desc',
                        group.format='simple', fl='journal, publication_date'))
#>   numFound start                    journal     publication_date
#> 1  2542432     0               PLOS Biology 2021-05-18T00:00:00Z
#> 2  2542432     0               PLOS Biology 2021-05-18T00:00:00Z
#> 3  2542432     0               PLOS Biology 2021-05-18T00:00:00Z
#> 4  2542432     0 PLOS Computational Biology 2021-05-18T00:00:00Z
#> 5  2542432     0 PLOS Computational Biology 2021-05-18T00:00:00Z
```

### Facet


```r
cli$facet(params = list(q='*:*', facet.field='journal', facet.query=c('cell', 'bird')))
#> $facet_queries
#> # A tibble: 2 x 2
#>   term   value
#>   <chr>  <int>
#> 1 cell  199069
#> 2 bird   21680
#> 
#> $facet_fields
#> $facet_fields$journal
#> # A tibble: 9 x 2
#>   term                             value  
#>   <chr>                            <chr>  
#> 1 plos one                         2110161
#> 2 plos genetics                    78511  
#> 3 plos neglected tropical diseases 75150  
#> 4 plos pathogens                   73595  
#> 5 plos computational biology       68336  
#> 6 plos biology                     45430  
#> 7 plos medicine                    33148  
#> 8 plos clinical trials             521    
#> 9 plos medicin                     9      
#> 
#> 
#> $facet_pivot
#> NULL
#> 
#> $facet_dates
#> NULL
#> 
#> $facet_ranges
#> NULL
```

### Highlight


```r
cli$highlight(params = list(q='alcohol', hl.fl = 'abstract', rows=2))
#> # A tibble: 2 x 2
#>   names                 abstract                                                
#>   <chr>                 <chr>                                                   
#> 1 10.1371/journal.pone… Background: Binge drinking, an increasingly common form…
#> 2 10.1371/journal.pone… Background and Aim: Harmful <em>alcohol</em> consumptio…
```

### Stats


```r
out <- cli$stats(params = list(q='ecology', stats.field=c('counter_total_all','alm_twitterCount'), stats.facet='journal'))
```


```r
out$data
#>                   min     max count missing       sum sumOfSquares        mean
#> counter_total_all   0 2629673 57016       0 372950117 2.403729e+13 6541.148397
#> alm_twitterCount    0    3804 57016       0    322545 8.667376e+07    5.657096
#>                        stddev
#> counter_total_all 19463.00439
#> alm_twitterCount     38.57705
```

### More like this

`solr_mlt` is a function to return similar documents to the one


```r
out <- cli$mlt(params = list(q='title:"ecology" AND body:"cell"', mlt.fl='title', mlt.mindf=1, mlt.mintf=1, fl='counter_total_all', rows=5))
```


```r
out$docs
#> # A tibble: 5 x 2
#>   id                           counter_total_all
#>   <chr>                                    <int>
#> 1 10.1371/journal.pbio.1001805             26190
#> 2 10.1371/journal.pbio.1002559             15937
#> 3 10.1371/journal.pbio.0020440             26740
#> 4 10.1371/journal.pone.0072451              4734
#> 5 10.1371/journal.pone.0087217             23421
```


```r
out$mlt
#> $`10.1371/journal.pbio.1001805`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     5450     0 10.1371/journal.pone.0098876              4455
#> 2     5450     0 10.1371/journal.pone.0082578              3750
#> 3     5450     0 10.1371/journal.pcbi.1007811               927
#> 4     5450     0 10.1371/journal.pone.0193049              3532
#> 5     5450     0 10.1371/journal.pone.0102159              2889
#> 
#> $`10.1371/journal.pbio.1002559`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     6857     0 10.1371/journal.pone.0155028              5121
#> 2     6857     0 10.1371/journal.pone.0041684             32685
#> 3     6857     0 10.1371/journal.pone.0023086             10853
#> 4     6857     0 10.1371/journal.pone.0155989              4195
#> 5     6857     0 10.1371/journal.pone.0223982              1233
#> 
#> $`10.1371/journal.pbio.0020440`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     1567     0 10.1371/journal.pone.0162651              4238
#> 2     1567     0 10.1371/journal.pone.0003259              3615
#> 3     1567     0 10.1371/journal.pone.0102679              5924
#> 4     1567     0 10.1371/journal.pone.0068814             10451
#> 5     1567     0 10.1371/journal.pntd.0003377              5000
#> 
#> $`10.1371/journal.pone.0072451`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1    30732     0 10.1371/journal.pntd.0004689              8298
#> 2    30732     0 10.1371/journal.pone.0000461             20728
#> 3    30732     0 10.1371/journal.pone.0006532             26214
#> 4    30732     0 10.1371/journal.ppat.0020122             10449
#> 5    30732     0 10.1371/journal.pone.0106526              3821
#> 
#> $`10.1371/journal.pone.0087217`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     6320     0 10.1371/journal.pone.0175497              2712
#> 2     6320     0 10.1371/journal.pone.0204743               558
#> 3     6320     0 10.1371/journal.pone.0159131              7088
#> 4     6320     0 10.1371/journal.pone.0220409              2326
#> 5     6320     0 10.1371/journal.pone.0123774              2728
```

### Parsing

`solr_parse` is a general purpose parser function with extension methods `solr_parse.sr_search`, `solr_parse.sr_facet`, and `solr_parse.sr_high`, for parsing `solr_search`, `solr_facet`, and `solr_highlight` function output, respectively. `solr_parse` is used internally within those three functions (`solr_search`, `solr_facet`, `solr_highlight`) to do parsing. You can optionally get back raw `json` or `xml` from `solr_search`, `solr_facet`, and `solr_highlight` setting parameter `raw=TRUE`, and then parsing after the fact with `solr_parse`. All you need to know is `solr_parse` can parse

For example:


```r
(out <- cli$highlight(params = list(q='alcohol', hl.fl = 'abstract', rows=2),
                      raw=TRUE))
#> [1] "{\n  \"response\":{\"numFound\":36140,\"start\":0,\"maxScore\":4.629626,\"docs\":[\n      {\n        \"id\":\"10.1371/journal.pone.0218147\",\n        \"journal\":\"PLOS ONE\",\n        \"eissn\":\"1932-6203\",\n        \"publication_date\":\"2019-12-10T00:00:00Z\",\n        \"article_type\":\"Research Article\",\n        \"author_display\":[\"Victor M. Jimenez Jr.\",\n          \"Erik W. Settles\",\n          \"Bart J. Currie\",\n          \"Paul S. Keim\",\n          \"Fernando P. Monroy\"],\n        \"abstract\":[\"Background: Binge drinking, an increasingly common form of alcohol use disorder, is associated with substantial morbidity and mortality; yet, its effects on the immune system’s ability to defend against infectious agents are poorly understood. Burkholderia pseudomallei, the causative agent of melioidosis can occur in healthy humans, yet binge alcohol intoxication is increasingly being recognized as a major risk factor. Although our previous studies demonstrated that binge alcohol exposure increased B. pseudomallei near-neighbor virulence in vivo and increased paracellular diffusion and intracellular invasion, no experimental studies have examined the extent to which bacterial and alcohol dosage play a role in disease progression. In addition, the temporal effects of a single binge alcohol dose prior to infection has not been examined in vivo. Principal findings: In this study, we used B. thailandensis E264 a close genetic relative of B. pseudomallei, as useful BSL-2 model system. Eight-week-old female C57BL/6 mice were utilized in three distinct animal models to address the effects of 1) bacterial dosage, 2) alcohol dosage, and 3) the temporal effects, of a single binge alcohol episode. Alcohol was administered comparable to human binge drinking (≤ 4.4 g/kg) or PBS intraperitoneally before a non-lethal intranasal infection. Bacterial colonization of lung and spleen was increased in mice administered alcohol even after bacterial dose was decreased 10-fold. Lung and not spleen tissue were colonized even after alcohol dosage was decreased 20 times below the U.S legal limit. Temporally, a single binge alcohol episode affected lung bacterial colonization for more than 24 h after alcohol was no longer detected in the blood. Pulmonary and splenic cytokine expression (TNF-α, GM-CSF) remained suppressed, while IL-12/p40 increased in mice administered alcohol 6 or 24 h prior to infection. Increased lung and not intestinal bacterial invasion was observed in human and murine non-phagocytic epithelial cells exposed to 0.2% v/v alcohol in vitro. Conclusions: Our results indicate that the effects of a single binge alcohol episode are tissue specific. A single binge alcohol intoxication event increases bacterial colonization in mouse lung tissue even after very low BACs and decreases the dose required to colonize the lungs with less virulent B. thailandensis. Additionally, the temporal effects of binge alcohol alters lung and spleen cytokine expression for at least 24 h after alcohol is detected in the blood. Delayed recovery in lung and not spleen tissue may provide a means for B. pseudomallei and near-neighbors to successfully colonize lung tissue through increased intracellular invasion of non-phagocytic cells in patients with hazardous alcohol intake. \"],\n        \"title_display\":\"Persistence of <i>Burkholderia thailandensis</i> E264 in lung tissue after a single binge alcohol episode\",\n        \"score\":4.629626},\n      {\n        \"id\":\"10.1371/journal.pone.0138021\",\n        \"journal\":\"PLOS ONE\",\n        \"eissn\":\"1932-6203\",\n        \"publication_date\":\"2015-09-16T00:00:00Z\",\n        \"article_type\":\"Research Article\",\n        \"author_display\":[\"Pavel Grigoriev\",\n          \"Evgeny M. Andreev\"],\n        \"abstract\":[\"Background and Aim: Harmful alcohol consumption has long been recognized as being the major determinant of male premature mortality in the European countries of the former USSR. Our focus here is on Belarus and Russia, two Slavic countries which continue to suffer enormously from the burden of the harmful consumption of alcohol. However, after a long period of deterioration, mortality trends in these countries have been improving over the past decade. We aim to investigate to what extent the recent declines in adult mortality in Belarus and Russia are attributable to the anti-alcohol measures introduced in these two countries in the 2000s. Data and Methods: We rely on the detailed cause-specific mortality series for the period 1980–2013. Our analysis focuses on the male population, and considers only a limited number of causes of death which we label as being alcohol-related: accidental poisoning by alcohol, liver cirrhosis, ischemic heart diseases, stroke, transportation accidents, and other external causes. For each of these causes we computed age-standardized death rates. The life table decomposition method was used to determine the age groups and the causes of death responsible for changes in life expectancy over time. Conclusion: Our results do not lead us to conclude that the schedule of anti-alcohol measures corresponds to the schedule of mortality changes. The continuous reduction in adult male mortality seen in Belarus and Russia cannot be fully explained by the anti-alcohol policies implemented in these countries, although these policies likely contributed to the large mortality reductions observed in Belarus and Russia in 2005–2006 and in Belarus in 2012. Thus, the effects of these policies appear to have been modest. We argue that the anti-alcohol measures implemented in Belarus and Russia simply coincided with fluctuations in alcohol-related mortality which originated in the past. If these trends had not been underway already, these huge mortality effects would not have occurred. \"],\n        \"title_display\":\"The Huge Reduction in Adult Male Mortality in Belarus and Russia: Is It Attributable to Anti-Alcohol Measures?\",\n        \"score\":4.627285}]\n  },\n  \"highlighting\":{\n    \"10.1371/journal.pone.0218147\":{\n      \"abstract\":[\"Background: Binge drinking, an increasingly common form of <em>alcohol</em> use disorder, is associated\"]},\n    \"10.1371/journal.pone.0138021\":{\n      \"abstract\":[\"Background and Aim: Harmful <em>alcohol</em> consumption has long been recognized as being the major\"]}}}\n"
#> attr(,"class")
#> [1] "sr_high"
#> attr(,"wt")
#> [1] "json"
```

Then parse


```r
solr_parse(out, 'df')
#> # A tibble: 2 x 2
#>   names                 abstract                                                
#>   <chr>                 <chr>                                                   
#> 1 10.1371/journal.pone… Background: Binge drinking, an increasingly common form…
#> 2 10.1371/journal.pone… Background and Aim: Harmful <em>alcohol</em> consumptio…
```

### Progress bars

only supported in the core search methods: `search`, `facet`, `group`, `mlt`, `stats`, `high`, `all`


```r
library(httr)
invisible(cli$search(params = list(q='*:*', rows=100, fl='id'), progress = httr::progress()))
|==============================================| 100%
```

### Advanced: Function Queries

Function Queries allow you to query on actual numeric fields in the SOLR database, and do addition, multiplication, etc on one or many fields to sort results. For example, here, we search on the product of counter_total_all and alm_twitterCount, using a new temporary field "_val_"


```r
cli$search(params = list(q='_val_:"product(counter_total_all,alm_twitterCount)"',
  rows=5, fl='id,title', fq='doc_type:full'))
#> # A tibble: 5 x 2
#>   id                     title                                                  
#>   <chr>                  <chr>                                                  
#> 1 10.1371/journal.pmed.… Why Most Published Research Findings Are False         
#> 2 10.1371/journal.pcbi.… Ten simple rules for structuring papers                
#> 3 10.1371/journal.pone.… A Multi-Level Bayesian Analysis of Racial Bias in Poli…
#> 4 10.1371/journal.pone.… Long-Term Follow-Up of Transsexual Persons Undergoing …
#> 5 10.1371/journal.pone.… More than 75 percent decline over 27 years in total fl…
```

Here, we search for the papers with the most citations


```r
cli$search(params = list(q='_val_:"max(counter_total_all)"',
    rows=5, fl='id,counter_total_all', fq='doc_type:full'))
#> # A tibble: 5 x 2
#>   id                           counter_total_all
#>   <chr>                                    <int>
#> 1 10.1371/journal.pmed.0020124           3256592
#> 2 10.1371/journal.pone.0133079           2629673
#> 3 10.1371/journal.pcbi.1003149           1794503
#> 4 10.1371/journal.pmed.1000376           1249497
#> 5 10.1371/journal.pmed.1000097           1012313
```

Or with the most tweets


```r
cli$search(params = list(q='_val_:"max(alm_twitterCount)"',
    rows=5, fl='id,alm_twitterCount', fq='doc_type:full'))
#> # A tibble: 5 x 2
#>   id                           alm_twitterCount
#>   <chr>                                   <int>
#> 1 10.1371/journal.pcbi.1005619             4935
#> 2 10.1371/journal.pmed.0020124             3890
#> 3 10.1371/journal.pone.0141854             3804
#> 4 10.1371/journal.pone.0115069             3083
#> 5 10.1371/journal.pmed.1001953             2825
```

### Using specific data sources

__USGS BISON service__

The occurrences service


```r
conn <- SolrClient$new(scheme = "https", host = "bison.usgs.gov", path = "solr/occurrences/select", port = NULL)
conn$search(params = list(q = '*:*', fl = c('decimalLatitude','decimalLongitude','scientificName'), rows = 2))
#> # A tibble: 2 x 3
#>   decimalLongitude scientificName     decimalLatitude
#>              <dbl> <chr>                        <dbl>
#> 1            -95.7 Oreothlypis celata            30.1
#> 2            -75.9 Oreothlypis celata            45.4
```

The species names service


```r
conn <- SolrClient$new(scheme = "https", host = "bison.usgs.gov", path = "solr/scientificName/select", port = NULL)
conn$search(params = list(q = '*:*'))
#> # A tibble: 10 x 2
#>    scientificName              `_version_`
#>    <chr>                             <dbl>
#>  1 Epuraea ambigua                 1.68e18
#>  2 Dictyopteris polypodioides      1.68e18
#>  3 Lonicera iberica                1.68e18
#>  4 Pseudopomala brachyptera        1.68e18
#>  5 Oceanococcus                    1.68e18
#>  6 Mactra alata                    1.68e18
#>  7 Reithrodontomys wetmorei        1.68e18
#>  8 Cristellaria orelliana          1.68e18
#>  9 Syringopora rara                1.68e18
#> 10 Aster cordifolius alvearius     1.68e18
```

__PLOS Search API__

Most of the examples above use the PLOS search API... :)

## Solr server management

This isn't as complete as searching functions show above, but we're getting there.

### Cores


```r
conn <- SolrClient$new()
```

Many functions, e.g.:

* `core_create()`
* `core_rename()`
* `core_status()`
* ...

Create a core


```r
conn$core_create(name = "foo_bar")
```

### Collections

Many functions, e.g.:

* `collection_create()`
* `collection_list()`
* `collection_addrole()`
* ...

Create a collection


```r
conn$collection_create(name = "hello_world")
```

### Add documents

Add documents, supports adding from files (json, xml, or csv format), and from R objects (including `data.frame` and `list` types so far)


```r
df <- data.frame(id = c(67, 68), price = c(1000, 500000000))
conn$add(df, name = "books")
```

Delete documents, by id


```r
conn$delete_by_id(name = "books", ids = c(3, 4))
```

Or by query


```r
conn$delete_by_query(name = "books", query = "manu:bank")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/solrium/issues)
* License: MIT
* Get citation information for `solrium` in R doing `citation(package = 'solrium')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
solrium 1.2.0
=============

### MINOR IMPROVEMENTS

* fix vignette titles (#123)
* vignette dependency fix


solrium 1.1.4
=============

### BUG FIXES

* fixed typo in code that made the `delete_by_query()`/`$delete_by_query()` method not work correctly (#121) thanks @abhik1368 for the report


solrium 1.1.0
=============

### MINOR IMPROVEMENTS

* all `data_frame` and `as_data_frame` usage converted to `as_tibble` (#119)
* change to markdown format docs
* fix some examples and update some broken URLs

### BUG FIXES

* group queries were failing because when there were no group results AND when response metadata was available it lead to a bug because you can't set attributes on `NULL`  (#118)


solrium 1.0.2
=============

### NEW FEATURES

* the major search methods on `SolrClient` and their function equivalents gain parameter `progress` that supports for now only `httr::progress()` (#115)
* gains new method `$json_request()` on `SolrClient` and new function `solr_json_request()` for working with the JSON request API (#117)

### MINOR IMPROVEMENTS

* now returning `responseHeader` and `nextCursorMark` when available as attributes on the returned object (#114)
* fixes for upcoming version of `tibble` (#116)


solrium 1.0.0
=============

This is v1, indicating breaking changes from the previous version!

### NEW FEATURES

* Package has been reworked to allow control over what parameters are sent
as query parameters and which as body. If only query parameters given, we do a
`GET` request, but if any body parameters given (even if query params given)
we do a `POST` request.  This means that all `solr_*` functions have more or
less the same parameters, and you now pass query parameters to `params` and
body parameters to `body`. This definitely breaks previous code, apologies
for that, but the bump in major version is a big indicator of the breakage.
* As part of overhaual, moved to using an `R6` setup for the Solr connection
object. The connection object deals with connection details, and you can call
all methods on the object created. Additionally, you can simply
pass the connection object to standalone methods. This change means
you can create connection objects to >1 Solr instance, so you can use many
Solr instances in one R session. (#100)
* gains new functions `update_atomic_json` and `update_atomic_xml` for doing
atomic updates (#97) thanks @yinghaoh
* `solr_search` and `solr_all` gain attributes that include `numFound`,
`start`, and `maxScore` (#94)
* `solr_search`/`solr_all`/`solr_mlt` gain new feature where we automically
check for and adjust `rows` parameter for you if you allow us to.
You can toggle this behavior and you can set a minimum number for rows
to be optimized with `minOptimizedRows`. See (#102) (#104) (#105) for
discussion. Thanks @1havran

### MINOR IMPROVEMENTS

* Replaced `httr` with `crul`. Should only be noticeable with respect
to specifying curl options (#98)
* Added more tests (#56)
* `optimize` renamed to `solr_optimize` (#107)
* now `solr_facet` fails better when no `facet.*` fields given (#103)

### BUG FIXES

* Fixed `solr_highlight` parsing to data.frame bug (#109)


solrium 0.4.0
=============

### MINOR IMPROVEMENTS

* Change `dplyr::rbind_all()` (deprecated) to `dplyr::bind_rows()` (#90)
* Added additional examples of using pivot facetting to `solr_facet()` (#91)
* Fix to `solr_group()` (#92)
* Replaced dependency `XML` with `xml2` (#57)
* Added examples and tests for a few more public Solr instances (#30)
* Now using `tibble` to give back compact data.frame's
* namespace all base package calls
* Many changes to internal parsers to use `xml2` instead of `XML`, and
improvements

solrium 0.3.0
=============

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.0.5 patched
* ubuntu 14.04 (on GitHub Actions), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 2 downstream dependencies.
  (Summary at <https://github.com/ropensci/solrium/blob/master/revdep/README.md>), with no problems found.

-----

This version fixes the rmarkdown/markdown dependency issue for vignettes that Kurt emailed maintainers about.

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

## Do not send the maintainer an email

Please open an issue instead of e-mailing the maintainer. E-mails will be a very low priority to answer.

## Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/solrium/issues) - be sure to include R session information and a reproducible example.

## Code contributions

### Broad overview of contributing workflow

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone git@github.com:<yourgithubusername>/solrium.git`
* Make sure to track progress upstream (i.e., on our version of `solrium` at `ropensci/solrium`) by doing `git remote add upstream git@github.com:ropensci/solrium.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs (see Tests below)
* Push up to your account
* Submit a pull request to home base at `ropensci/solrium`

### Tests

To add tests, go to the folder `tests/testthat/`. Tests are generally organized as individual files for each exported function from the package (that is, listed as an export in the `NAMESPACE` file). If you are adding a new exported function, add a new test file. If you are changing an existing function, work in the tests file for that function, unless it doesn't have tests, in which case make a new test file.

The book R packages book provides [a chapter on testing in general](http://r-pkgs.had.co.nz/tests.html). Do consult that first if you aren't familiar with testing in R.

The easiest set up to run tests is from within an R session:

```r
library(devtools)
library(testthat)
# loads the package
load_all()
```

To test an individual test file

```r
test_file("tests/testthat/test-foobar.R")
```

To run all tests

```r
devtools::test()
```

Or you can run from the CLI like `Rscript -e "devtools::test()"`

If you are running tests that have `skip_on_cran()` in them, set `Sys.setenv(NOT_CRAN = "true")` prior to running tests.


### Making changes

In addition to changing the code, do make sure to udpate the documentation if applicable. The R packages book book has a [chapter on documentation](http://r-pkgs.had.co.nz/man.html) you should read if you aren't familiar.

After code and documentation has been changed, update documentation by running either `devtools::document()` or `roxygen2::roxygenise()`.

Make sure if you change what packages or even functions within packages are imported, most likely add the package to Imports in the DESCRIPTION file and list what functions are imported in the `solrium-package.R` file.

Be very conservative about adding new dependencies.

### Style

* Make sure code, documentation, and comments are no more than 80 characters in width.
* Use `<-` instead of `=` for assignment
* Always use `snake_case` (and all lowercase) instead of `camelCase`

## Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) 

In addition, if this concerns a Solr installation you manage or have admin access to, please include the Solr version.
-->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Cores/collections management}
%\VignetteEncoding{UTF-8}
-->



Cores/collections management
============================

## Installation

Stable version from CRAN


```r
install.packages("solrium")
```

Or the development version from GitHub


```r
install.packages("devtools")
devtools::install_github("ropensci/solrium")
```

Load


```r
library("solrium")
```

Initialize connection


```r
(conn <- SolrClient$new())
```

```
#> <Solr Client>
#>   host: 127.0.0.1
#>   path: 
#>   port: 8983
#>   scheme: http
#>   errors: simple
#>   proxy:
```

## Cores

There are many operations you can do on cores, including:

* `core_create()` - create a core
* `core_exists()` - check if a core exists
* `core_mergeindexes()` - merge indexes
* `core_reload()` - reload a core
* `core_rename()` - rename a core
* `core_requeststatus()` - check request status
* `core_split()` - split a core
* `core_status()` - check core status
* `core_swap()` - core swap
* `core_unload()` - delete a core

### Create a core


```r
conn$core_create()
```

### Delete a core


```r
conn$core_unload()
```

## Collections

There are many operations you can do on collections, including:

* `collection_addreplica()`
* `collection_addreplicaprop()`
* `collection_addrole()`
* `collection_balanceshardunique()`
* `collection_clusterprop()`
* `collection_clusterstatus()`
* `collection_create()`
* `collection_createalias()`
* `collection_createshard()`
* `collection_delete()`
* `collection_deletealias()`
* `collection_deletereplica()`
* `collection_deletereplicaprop()`
* `collection_deleteshard()`
* `collection_list()`
* `collection_migrate()`
* `collection_overseerstatus()`
* `collection_rebalanceleaders()`
* `collection_reload()`
* `collection_removerole()`
* `collection_requeststatus()`
* `collection_splitshard()`

### Create a collection


```r
conn$collection_create()
```

### Delete a collection


```r
conn$collection_delete()
```
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Document management}
%\VignetteEncoding{UTF-8}
-->



Document management
===================

## Installation

Stable version from CRAN


```r
install.packages("solrium")
```

Or the development version from GitHub


```r
install.packages("devtools")
devtools::install_github("ropensci/solrium")
```

Load


```r
library("solrium")
```

Initialize connection. By default, you connect to `http://localhost:8983`


```r
(conn <- SolrClient$new())
```

```
#> <Solr Client>
#>   host: 127.0.0.1
#>   path: 
#>   port: 8983
#>   scheme: http
#>   errors: simple
#>   proxy:
```

## Create documents from R objects

For now, only lists and data.frame's supported.

### data.frame


```r
df <- data.frame(id = c(67, 68), price = c(1000, 500000000))
conn$add(df, "books")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 319
```

### list




```r
ss <- list(list(id = 1, price = 100), list(id = 2, price = 500))
conn$add(ss, "books")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 47
```

## Delete documents

### By id

Add some documents first




```r
docs <- list(list(id = 1, price = 100, name = "brown"),
             list(id = 2, price = 500, name = "blue"),
             list(id = 3, price = 2000L, name = "pink"))
conn$add(docs, "gettingstarted")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 206
```

And the documents are now in your Solr database


```r
conn$search(name = "gettingstarted", params = list(q = "*:*", rows = 3))
```

```
#> # A tibble: 3 x 3
#>   id    title   `_version_`
#>   <chr> <chr>         <dbl>
#> 1 10    adfadsf     1.65e18
#> 2 12    though      1.65e18
#> 3 14    animals     1.65e18
```

Now delete those documents just added


```r
conn$delete_by_id(ids = c(1, 2, 3), "gettingstarted")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 151
```

And now they are gone


```r
conn$search("gettingstarted", params = list(q = "*:*", rows = 4))
```

```
#> # A tibble: 3 x 3
#>   id    title   `_version_`
#>   <chr> <chr>         <dbl>
#> 1 10    adfadsf     1.65e18
#> 2 12    though      1.65e18
#> 3 14    animals     1.65e18
```

### By query

Add some documents first


```r
conn$add(docs, "gettingstarted")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 40
```

And the documents are now in your Solr database


```r
conn$search("gettingstarted", params = list(q = "*:*", rows = 5))
```

```
#> # A tibble: 5 x 5
#>   id    title   `_version_` price name 
#>   <chr> <chr>         <dbl> <int> <chr>
#> 1 10    adfadsf     1.65e18    NA <NA> 
#> 2 12    though      1.65e18    NA <NA> 
#> 3 14    animals     1.65e18    NA <NA> 
#> 4 1     <NA>        1.65e18   100 brown
#> 5 2     <NA>        1.65e18   500 blue
```

Now delete those documents just added


```r
conn$delete_by_query(query = "(name:blue OR name:pink)", "gettingstarted")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 74
```

And now they are gone


```r
conn$search("gettingstarted", params = list(q = "*:*", rows = 5))
```

```
#> # A tibble: 4 x 5
#>   id    title   `_version_` price name 
#>   <chr> <chr>         <dbl> <int> <chr>
#> 1 10    adfadsf     1.65e18    NA <NA> 
#> 2 12    though      1.65e18    NA <NA> 
#> 3 14    animals     1.65e18    NA <NA> 
#> 4 1     <NA>        1.65e18   100 brown
```

## Update documents from files

This approach is best if you have many different things you want to do at once, e.g., delete and add files and set any additional options. The functions are:

* `update_xml()`
* `update_json()`
* `update_csv()`

There are separate functions for each of the data types as they take slightly different parameters - and to make it more clear that those are the three input options for data types.

### JSON


```r
file <- system.file("examples", "books.json", package = "solrium")
conn$update_json(file, "books")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 61
```

### Add and delete in the same file

Add a document first, that we can later delete


```r
ss <- list(list(id = 456, name = "cat"))
conn$add(ss, "books")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 43
```

Now add a new document, and delete the one we just made


```r
file <- system.file("examples", "add_delete.xml", package = "solrium")
cat(readLines(file), sep = "\n")
```

```
#> <update>
#> 	<add>
#> 	  <doc>
#> 	    <field name="id">978-0641723445</field>
#> 	    <field name="cat">book,hardcover</field>
#> 	    <field name="name">The Lightning Thief</field>
#> 	    <field name="author">Rick Riordan</field>
#> 	    <field name="series_t">Percy Jackson and the Olympians</field>
#> 	    <field name="sequence_i">1</field>
#> 	    <field name="genre_s">fantasy</field>
#> 	    <field name="inStock">TRUE</field>
#> 	    <field name="pages_i">384</field>
#> 	  </doc>
#> 	</add>
#> 	<delete>
#> 		<id>456</id>
#> 	</delete>
#> </update>
```

```r
conn$update_xml(file, "books")
```

```
#> $responseHeader
#> $responseHeader$rf
#> [1] 1
#> 
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 95
```

### Notes

Note that `update_xml()` and `update_json()` have exactly the same parameters, but simply use different data input formats. `update_csv()` is different in that you can't provide document or field level boosts or other modifications. In addition `update_csv()` can accept not just csv, but tsv and other types of separators.

<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Solr search}
%\VignetteEncoding{UTF-8}
-->



Solr search
===========

**A general purpose R interface to [Apache Solr](https://lucene.apache.org/solr/)**

## Solr info

+ [Solr home page](https://lucene.apache.org/solr/)
+ [Highlighting help](https://lucene.apache.org/solr/guide/8_2/highlighting.html)
+ [Faceting help](https://lucene.apache.org/solr/guide/8_2/faceting.html)
+ [Install and Setup SOLR in OSX, including running Solr](http://risnandar.wordpress.com/2013/09/08/how-to-install-and-setup-apache-lucene-solr-in-osx/)

## Installation

Stable version from CRAN


```r
install.packages("solrium")
```

Or the development version from GitHub


```r
install.packages("devtools")
devtools::install_github("ropensci/solrium")
```

Load


```r
library("solrium")
```

## Setup connection

You can setup for a remote Solr instance or on your local machine.


```r
(conn <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))
#> <Solr Client>
#>   host: api.plos.org
#>   path: search
#>   port: 
#>   scheme: http
#>   errors: simple
#>   proxy:
```

## Rundown

`solr_search()` only returns the `docs` element of a Solr response body. If `docs` is
all you need, then this function will do the job. If you need facet data only, or mlt
data only, see the appropriate functions for each of those below. Another function,
`solr_all()` has a similar interface in terms of parameter as `solr_search()`, but
returns all parts of the response body, including, facets, mlt, groups, stats, etc.
as long as you request those.

## Search docs

`solr_search()` returns only docs. A basic search:


```r
conn$search(params = list(q = '*:*', rows = 2, fl = 'id'))
#> # A tibble: 2 x 1
#>   id                                                
#>   <chr>                                             
#> 1 10.1371/journal.pone.0058099/materials_and_methods
#> 2 10.1371/journal.pone.0030394/introduction
```

__Search in specific fields with `:`__

Search for word ecology in title and cell in the body


```r
conn$search(params = list(q = 'title:"ecology" AND body:"cell"', fl = 'title', rows = 5))
#> # A tibble: 5 x 1
#>   title                                                    
#>   <chr>                                                    
#> 1 The Ecology of Collective Behavior                       
#> 2 Ecology's Big, Hot Idea                                  
#> 3 Chasing Ecological Interactions                          
#> 4 Spatial Ecology of Bacteria at the Microscale in Soil    
#> 5 Biofilm Formation As a Response to Ecological Competition
```

__Wildcards__

Search for word that starts with "cell" in the title field


```r
conn$search(params = list(q = 'title:"cell*"', fl = 'title', rows = 5))
#> # A tibble: 5 x 1
#>   title                                                                    
#>   <chr>                                                                    
#> 1 Cancer Stem Cell-Like Side Population Cells in Clear Cell Renal Cell Car…
#> 2 Tumor Cell Recognition Efficiency by T Cells                             
#> 3 Enhancement of Chemotactic Cell Aggregation by Haptotactic Cell-To-Cell …
#> 4 Cell-Cell Adhesions and Cell Contractility Are Upregulated upon Desmosom…
#> 5 Dcas Supports Cell Polarization and Cell-Cell Adhesion Complexes in Deve…
```

__Proximity search__

Search for words "sports" and "alcohol" within four words of each other


```r
conn$search(params = list(q = 'everything:"stem cell"~7', fl = 'title', rows = 3))
#> # A tibble: 3 x 1
#>   title                                                                    
#>   <chr>                                                                    
#> 1 Effect of Dedifferentiation on Time to Mutation Acquisition in Stem Cell…
#> 2 A Mathematical Model of Cancer Stem Cell Driven Tumor Initiation: Implic…
#> 3 Phenotypic Evolutionary Models in Stem Cell Biology: Replacement, Quiesc…
```

__Range searches__

Search for articles with Twitter count between 5 and 10


```r
conn$search(params = list(q = '*:*', fl = c('alm_twitterCount', 'id'), fq = 'alm_twitterCount:[5 TO 50]', rows = 10))
#> # A tibble: 10 x 2
#>    id                                                  alm_twitterCount
#>    <chr>                                                          <int>
#>  1 10.1371/journal.pbio.0030378/title                                 8
#>  2 10.1371/journal.pbio.0030378/abstract                              8
#>  3 10.1371/journal.pbio.0030378/references                            8
#>  4 10.1371/journal.pone.0184491                                      10
#>  5 10.1371/journal.pone.0184491/title                                10
#>  6 10.1371/journal.pone.0184491/abstract                             10
#>  7 10.1371/journal.pone.0184491/references                           10
#>  8 10.1371/journal.pone.0184491/body                                 10
#>  9 10.1371/journal.pone.0184491/introduction                         10
#> 10 10.1371/journal.pone.0184491/results_and_discussion               10
```

__Boosts__

Assign higher boost to title matches than to body matches (compare the two calls)


```r
conn$search(params = list(q = 'title:"cell" abstract:"science"', fl = 'title', rows = 3))
#> # A tibble: 3 x 1
#>   title                                                                    
#>   <chr>                                                                    
#> 1 I Want More and Better Cells! – An Outreach Project about Stem Cells and…
#> 2 Globalization of Stem Cell Science: An Examination of Current and Past C…
#> 3 Virtual Reconstruction and Three-Dimensional Printing of Blood Cells as …
```


```r
conn$search(params = list(q = 'title:"cell"^1.5 AND abstract:"science"', fl = 'title', rows = 3))
#> # A tibble: 3 x 1
#>   title                                                                    
#>   <chr>                                                                    
#> 1 I Want More and Better Cells! – An Outreach Project about Stem Cells and…
#> 2 Virtual Reconstruction and Three-Dimensional Printing of Blood Cells as …
#> 3 Globalization of Stem Cell Science: An Examination of Current and Past C…
```

## Search all

`solr_all()` differs from `solr_search()` in that it allows specifying facets, mlt, groups,
stats, etc, and returns all of those. It defaults to `parsetype = "list"` and `wt="json"`,
whereas `solr_search()` defaults to `parsetype = "df"` and `wt="csv"`. `solr_all()` returns
by default a list, whereas `solr_search()` by default returns a data.frame.

A basic search, just docs output


```r
conn$all(params = list(q = '*:*', rows = 2, fl = 'id'))
#> $search
#> # A tibble: 2 x 1
#>   id                                                
#>   <chr>                                             
#> 1 10.1371/journal.pone.0058099/materials_and_methods
#> 2 10.1371/journal.pone.0030394/introduction         
#> 
#> $facet
#> list()
#> 
#> $high
#> # A tibble: 0 x 0
#> 
#> $mlt
#> $mlt$docs
#> # A tibble: 2 x 1
#>   id                                                
#>   <chr>                                             
#> 1 10.1371/journal.pone.0058099/materials_and_methods
#> 2 10.1371/journal.pone.0030394/introduction         
#> 
#> $mlt$mlt
#> list()
#> 
#> 
#> $group
#>   numFound start                                                 id
#> 1  2263584     0 10.1371/journal.pone.0058099/materials_and_methods
#> 2  2263584     0          10.1371/journal.pone.0030394/introduction
#> 
#> $stats
#> NULL
```

Get docs, mlt, and stats output


```r
conn$all(params = list(q = 'ecology', rows = 2, fl = 'id', mlt = 'true', mlt.count = 2, mlt.fl = 'abstract', stats = 'true', stats.field = 'counter_total_all'))
#> $search
#> # A tibble: 2 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pone.0001248
#> 2 10.1371/journal.pone.0059813
#> 
#> $facet
#> list()
#> 
#> $high
#> # A tibble: 0 x 0
#> 
#> $mlt
#> $mlt$docs
#> # A tibble: 2 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pone.0001248
#> 2 10.1371/journal.pone.0059813
#> 
#> $mlt$mlt
#> $mlt$mlt$`10.1371/journal.pone.0001248`
#> # A tibble: 2 x 3
#>   numFound start id                          
#>      <int> <int> <chr>                       
#> 1   236603     0 10.1371/journal.pbio.1002448
#> 2   236603     0 10.1371/journal.pone.0155843
#> 
#> $mlt$mlt$`10.1371/journal.pone.0059813`
#> # A tibble: 2 x 3
#>   numFound start id                          
#>      <int> <int> <chr>                       
#> 1   228703     0 10.1371/journal.pone.0204749
#> 2   228703     0 10.1371/journal.pone.0175014
#> 
#> 
#> 
#> $group
#>   numFound start                           id
#> 1    49638     0 10.1371/journal.pone.0001248
#> 2    49638     0 10.1371/journal.pone.0059813
#> 
#> $stats
#> $stats$data
#>                   min     max count missing       sum sumOfSquares    mean
#> counter_total_all   0 1322780 49638       0 264206214 1.119659e+13 5322.66
#>                     stddev
#> counter_total_all 14044.15
#> 
#> $stats$facet
#> NULL
```


## Facet


```r
conn$facet(params = list(q = '*:*', facet.field = 'journal', facet.query = c('cell', 'bird')))
#> $facet_queries
#> # A tibble: 2 x 2
#>   term   value
#>   <chr>  <int>
#> 1 cell  181404
#> 2 bird   19370
#> 
#> $facet_fields
#> $facet_fields$journal
#> # A tibble: 9 x 2
#>   term                             value  
#>   <fct>                            <fct>  
#> 1 plos one                         1878564
#> 2 plos genetics                    69743  
#> 3 plos pathogens                   62807  
#> 4 plos neglected tropical diseases 61216  
#> 5 plos computational biology       56361  
#> 6 plos biology                     39732  
#> 7 plos medicine                    27839  
#> 8 plos clinical trials             521    
#> 9 plos medicin                     9      
#> 
#> 
#> $facet_pivot
#> NULL
#> 
#> $facet_dates
#> NULL
#> 
#> $facet_ranges
#> NULL
```

## Highlight


```r
conn$highlight(params = list(q = 'alcohol', hl.fl = 'abstract', rows = 2))
#> # A tibble: 2 x 2
#>   names                 abstract                                           
#>   <chr>                 <chr>                                              
#> 1 10.1371/journal.pone… "\nAcute <em>alcohol</em> administration can lead …
#> 2 10.1371/journal.pone… Objectives: <em>Alcohol</em>-related morbidity and…
```

## Stats


```r
out <- conn$stats(params = list(q = 'ecology', stats.field = c('counter_total_all', 'alm_twitterCount'), stats.facet = c('journal', 'volume')))
```


```r
out$data
#>                   min     max count missing       sum sumOfSquares
#> counter_total_all   0 1322780 49638       0 264206214 1.119659e+13
#> alm_twitterCount    0    3438 49638       0    304236 8.148536e+07
#>                          mean     stddev
#> counter_total_all 5322.660341 14044.1498
#> alm_twitterCount     6.129095    40.0507
```


```r
out$facet
#> $counter_total_all
#> $counter_total_all$volume
#> # A tibble: 17 x 9
#>    volume   min     max count missing      sum  sumOfSquares   mean stddev
#>    <chr>  <dbl>   <dbl> <int>   <int>    <dbl>         <dbl>  <dbl>  <dbl>
#>  1 11         0  287793  5264       0 18895052  343471377354  3589.  7237.
#>  2 12         0  516798  5076       0 14049766  457686530584  2768.  9084.
#>  3 13         0  153970  4715       0  5936439   86585551129  1259.  4097.
#>  4 14         0  141744  3293       0  2631930   53269058704   799.  3942.
#>  5 15         0   50449   380       0  1602416   24738162692  4217.  6888.
#>  6 16         0   36261   154       0   829924    8698444092  5389.  5255.
#>  7 17         0   42811   108       0   228662    2794636940  2117.  4647.
#>  8 1       2127  331494    81       0  1989864  220314081736 24566. 46291.
#>  9 2       2141  155085   482       0  7419080  277345781732 15392. 18417.
#> 10 3       1621  138038   741       0  9765264  305018983590 13178. 15436.
#> 11 4        866  404016  1010       0 12543775  550954124849 12420. 19790.
#> 12 5        125  248082  1539       0 16437056  451150556196 10680. 13386.
#> 13 6         95  396324  2948       0 25087216  729164435884  8510. 13228.
#> 14 7         62  270034  4825       0 36102767  872710683283  7482. 11176.
#> 15 8         34  611601  6360       0 42214728 1329098761134  6638. 12843.
#> 16 9         57 1322780  6620       0 40933803 3662802686107  6183. 22697.
#> 17 10       428  887162  6042       0 27538472 1820785780330  4558. 16752.
#> 
#> $counter_total_all$journal
#> # A tibble: 9 x 9
#>   journal   min     max count missing       sum  sumOfSquares   mean stddev
#>   <chr>   <dbl>   <dbl> <int>   <int>     <dbl>         <dbl>  <dbl>  <dbl>
#> 1 1           0  391184  1272       0  22099069 1251395524029 17373. 26125.
#> 2 2           0  270034   335       0   6608325  400062632503 19726. 28417.
#> 3 3        1354  221758  1249       0   9088713  227814882055  7277. 11382.
#> 4 4        9348   18512     2       0     27860     430079248 13930   6480.
#> 5 5           0  887162 41465       0 187628028 6833430419632  4525. 12014.
#> 6 6           0  138374   906       0   7870667  145153611759  8687.  9211.
#> 7 7           0  153848  1101       0  10283498  199515141164  9340.  9698.
#> 8 8           0  323752  2290       0  12009181  244148093179  5244.  8897.
#> 9 9           0 1322780  1018       0   8590873 1894639252767  8439. 42328.
#> 
#> 
#> $alm_twitterCount
#> $alm_twitterCount$volume
#> # A tibble: 17 x 9
#>    volume   min   max count missing   sum sumOfSquares  mean stddev
#>    <chr>  <dbl> <dbl> <int>   <int> <dbl>        <dbl> <dbl>  <dbl>
#>  1 11         0  2142  5264       0 52248     11395704  9.93  45.5 
#>  2 12         0  1877  5076       0 38615     10903883  7.61  45.7 
#>  3 13         0   578  4715       0 15920      2404482  3.38  22.3 
#>  4 14         0   984  3293       0 14247      3956965  4.33  34.4 
#>  5 15         0   453   380       0  5640       906810 14.8   46.6 
#>  6 16         0   456   154       0  2830       414458 18.4   48.7 
#>  7 17         0    44   108       0   170         4882  1.57   6.57
#>  8 1          0    47    81       0   208         6306  2.57   8.49
#>  9 2          0   125   482       0  1116        63658  2.32  11.3 
#> 10 3          0   504   741       0  1407       271861  1.90  19.1 
#> 11 4          0   313  1010       0  1655       155545  1.64  12.3 
#> 12 5          0   165  1539       0  2694       142070  1.75   9.45
#> 13 6          0   968  2948       0  5709      1631337  1.94  23.4 
#> 14 7          0   860  4825       0 21860      2578512  4.53  22.7 
#> 15 8          0  2029  6360       0 40719     10854695  6.40  40.8 
#> 16 9          0  1880  6620       0 55292     17030904  8.35  50.0 
#> 17 10         0  3438  6042       0 43906     18763286  7.27  55.3 
#> 
#> $alm_twitterCount$journal
#> # A tibble: 9 x 9
#>   journal   min   max count missing    sum sumOfSquares  mean stddev
#>   <chr>   <dbl> <dbl> <int>   <int>  <dbl>        <dbl> <dbl>  <dbl>
#> 1 1           0  2142  1272       0  39045     15091795 30.7  105.  
#> 2 2           0   831   335       0   5826      1217690 17.4   57.8 
#> 3 3           0   455  1249       0   7368       576290  5.90  20.7 
#> 4 4           0     3     2       0      3            9  1.5    2.12
#> 5 5           0  3438 41465       0 211665     60964301  5.10  38.0 
#> 6 6           0   250   906       0   8476       414748  9.36  19.3 
#> 7 7           0   230  1101       0   9237       461117  8.39  18.7 
#> 8 8           0   968  2290       0  11804      1555638  5.15  25.6 
#> 9 9           0   578  1018       0  10812      1203770 10.6   32.7
```

## More like this

`solr_mlt` is a function to return similar documents to the one


```r
out <- conn$mlt(params = list(q = 'title:"ecology" AND body:"cell"', mlt.fl = 'title', mlt.mindf = 1, mlt.mintf = 1, fl = 'counter_total_all', rows = 5))
out$docs
#> # A tibble: 5 x 2
#>   id                           counter_total_all
#>   <chr>                                    <int>
#> 1 10.1371/journal.pbio.1001805             23958
#> 2 10.1371/journal.pbio.0020440             26090
#> 3 10.1371/journal.pbio.1002559             11628
#> 4 10.1371/journal.pone.0087217             16196
#> 5 10.1371/journal.pbio.1002191             27371
```


```r
out$mlt
#> $`10.1371/journal.pbio.1001805`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     4678     0 10.1371/journal.pone.0098876              4047
#> 2     4678     0 10.1371/journal.pone.0082578              3244
#> 3     4678     0 10.1371/journal.pone.0102159              2434
#> 4     4678     0 10.1371/journal.pone.0193049              1274
#> 5     4678     0 10.1371/journal.pcbi.1003408             11685
#> 
#> $`10.1371/journal.pbio.0020440`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     1375     0 10.1371/journal.pone.0162651              3463
#> 2     1375     0 10.1371/journal.pone.0003259              3417
#> 3     1375     0 10.1371/journal.pntd.0003377              4613
#> 4     1375     0 10.1371/journal.pone.0068814              9701
#> 5     1375     0 10.1371/journal.pone.0101568              6017
#> 
#> $`10.1371/journal.pbio.1002559`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     6288     0 10.1371/journal.pone.0155028              2881
#> 2     6288     0 10.1371/journal.pone.0023086              9361
#> 3     6288     0 10.1371/journal.pone.0041684             26571
#> 4     6288     0 10.1371/journal.pone.0155989              2519
#> 5     6288     0 10.1371/journal.pone.0129394              2111
#> 
#> $`10.1371/journal.pone.0087217`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     5565     0 10.1371/journal.pone.0204743               103
#> 2     5565     0 10.1371/journal.pone.0175497              1088
#> 3     5565     0 10.1371/journal.pone.0159131              4937
#> 4     5565     0 10.1371/journal.pcbi.0020092             26453
#> 5     5565     0 10.1371/journal.pone.0133941              1336
#> 
#> $`10.1371/journal.pbio.1002191`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1    14595     0 10.1371/journal.pbio.1002232              3055
#> 2    14595     0 10.1371/journal.pone.0191705              1040
#> 3    14595     0 10.1371/journal.pone.0070448              2497
#> 4    14595     0 10.1371/journal.pone.0131700              3353
#> 5    14595     0 10.1371/journal.pone.0121680              4980
```

## Groups

`solr_groups()` is a function to return similar documents to the one


```r
conn$group(params = list(q = 'ecology', group.field = 'journal', group.limit = 1, fl = c('id', 'alm_twitterCount')))
#>                         groupValue numFound start
#> 1                         plos one    41465     0
#> 2       plos computational biology     1018     0
#> 3                     plos biology     1272     0
#> 4 plos neglected tropical diseases     2290     0
#> 5                   plos pathogens      906     0
#> 6                    plos genetics     1101     0
#> 7                             none     1249     0
#> 8                    plos medicine      335     0
#> 9             plos clinical trials        2     0
#>                             id alm_twitterCount
#> 1 10.1371/journal.pone.0001248                0
#> 2 10.1371/journal.pcbi.1003594               21
#> 3 10.1371/journal.pbio.0060300                0
#> 4 10.1371/journal.pntd.0004689               13
#> 5 10.1371/journal.ppat.1005780               19
#> 6 10.1371/journal.pgen.1005860              135
#> 7 10.1371/journal.pone.0043894                0
#> 8 10.1371/journal.pmed.1000303                1
#> 9 10.1371/journal.pctr.0020010                0
```

## Parsing

`solr_parse()` is a general purpose parser function with extension methods for parsing outputs from functions in `solr`. `solr_parse()` is used internally within functions to do parsing after retrieving data from the server. You can optionally get back raw `json`, `xml`, or `csv` with the `raw=TRUE`, and then parse afterwards with `solr_parse()`.

For example:


```r
(out <- conn$highlight(params = list(q = 'alcohol', hl.fl = 'abstract', rows = 2), raw = TRUE))
#> [1] "{\n  \"response\":{\"numFound\":31528,\"start\":0,\"maxScore\":4.6573215,\"docs\":[\n      {\n        \"id\":\"10.1371/journal.pone.0201042\",\n        \"journal\":\"PLOS ONE\",\n        \"eissn\":\"1932-6203\",\n        \"publication_date\":\"2018-07-26T00:00:00Z\",\n        \"article_type\":\"Research Article\",\n        \"author_display\":[\"Graeme Knibb\",\n          \"Carl. A. Roberts\",\n          \"Eric Robinson\",\n          \"Abi Rose\",\n          \"Paul Christiansen\"],\n        \"abstract\":[\"\\nAcute alcohol administration can lead to a loss of control over drinking. Several models argue that this ‘alcohol priming effect’ is mediated by the effect of alcohol on inhibitory control. Alternatively, beliefs about how alcohol affects behavioural regulation may also underlie alcohol priming and alcohol-induced inhibitory impairments. Here two studies examine the extent to which the alcohol priming effect and inhibitory impairments are moderated by beliefs regarding the effects of alcohol on the ability to control behaviour. In study 1, following a priming drink (placebo or .5g/kg of alcohol), participants were provided with bogus feedback regarding their performance on a measure of inhibitory control (stop-signal task; SST) suggesting that they had high or average self-control. However, the bogus feedback manipulation was not successful. In study 2, before a SST, participants were exposed to a neutral or experimental message suggesting acute doses of alcohol reduce the urge to drink and consumed a priming drink and this manipulation was successful. In both studies craving was assessed throughout and a bogus taste test which measured ad libitum drinking was completed. Results suggest no effect of beliefs on craving or ad lib consumption within either study. However, within study 2, participants exposed to the experimental message displayed evidence of alcohol-induced impairments of inhibitory control, while those exposed to the neutral message did not. These findings do not suggest beliefs about the effects of alcohol moderate the alcohol priming effect but do suggest beliefs may, in part, underlie the effect of alcohol on inhibitory control.\\n\"],\n        \"title_display\":\"The effect of beliefs about alcohol’s acute effects on alcohol priming and alcohol-induced impairments of inhibitory control\",\n        \"score\":4.6573215},\n      {\n        \"id\":\"10.1371/journal.pone.0185457\",\n        \"journal\":\"PLOS ONE\",\n        \"eissn\":\"1932-6203\",\n        \"publication_date\":\"2017-09-28T00:00:00Z\",\n        \"article_type\":\"Research Article\",\n        \"author_display\":[\"Jacqueline Willmore\",\n          \"Terry-Lynne Marko\",\n          \"Darcie Taing\",\n          \"Hugues Sampasa-Kanyinga\"],\n        \"abstract\":[\"Objectives: Alcohol-related morbidity and mortality are significant public health issues. The purpose of this study was to describe the prevalence and trends over time of alcohol consumption and alcohol-related morbidity and mortality; and public attitudes of alcohol use impacts on families and the community in Ottawa, Canada. Methods: Prevalence (2013–2014) and trends (2000–2001 to 2013–2014) of alcohol use were obtained from the Canadian Community Health Survey. Data on paramedic responses (2015), emergency department (ED) visits (2013–2015), hospitalizations (2013–2015) and deaths (2007–2011) were used to quantify the acute and chronic health effects of alcohol in Ottawa. Qualitative data were obtained from the “Have Your Say” alcohol survey, an online survey of public attitudes on alcohol conducted in 2016. Results: In 2013–2014, an estimated 595,300 (83%) Ottawa adults 19 years and older drank alcohol, 42% reported binge drinking in the past year. Heavy drinking increased from 15% in 2000–2001 to 20% in 2013–2014. In 2015, the Ottawa Paramedic Service responded to 2,060 calls directly attributable to alcohol. Between 2013 and 2015, there were an average of 6,100 ED visits and 1,270 hospitalizations per year due to alcohol. Annually, alcohol use results in at least 140 deaths in Ottawa. Men have higher rates of alcohol-attributable paramedic responses, ED visits, hospitalizations and deaths than women, and young adults have higher rates of alcohol-attributable paramedic responses. Qualitative data of public attitudes indicate that alcohol misuse has greater repercussions not only on those who drink, but also on the family and community. Conclusions: Results highlight the need for healthy public policy intended to encourage a culture of drinking in moderation in Ottawa to support lower risk alcohol use, particularly among men and young adults. \"],\n        \"title_display\":\"The burden of alcohol-related morbidity and mortality in Ottawa, Canada\",\n        \"score\":4.65702}]\n  },\n  \"highlighting\":{\n    \"10.1371/journal.pone.0201042\":{\n      \"abstract\":[\"\\nAcute <em>alcohol</em> administration can lead to a loss of control over drinking. Several models argue\"]},\n    \"10.1371/journal.pone.0185457\":{\n      \"abstract\":[\"Objectives: <em>Alcohol</em>-related morbidity and mortality are significant public health issues\"]}}}\n"
#> attr(,"class")
#> [1] "sr_high"
#> attr(,"wt")
#> [1] "json"
```

Then parse


```r
solr_parse(out, 'df')
#> # A tibble: 2 x 2
#>   names                 abstract                                           
#>   <chr>                 <chr>                                              
#> 1 10.1371/journal.pone… "\nAcute <em>alcohol</em> administration can lead …
#> 2 10.1371/journal.pone… Objectives: <em>Alcohol</em>-related morbidity and…
```

[Please report any issues or bugs](https://github.com/ropensci/solrium/issues).
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Local Solr setup}
%\VignetteEncoding{UTF-8}
-->

Local Solr setup
======

The Solr version you are working with my differ from below. Don't worry, just go with the version you want to use.

### OSX

__Based on https://lucene.apache.org/solr/guide/8_2/solr-tutorial.html__

1. Download most recent version from an Apache mirror https://lucene.apache.org/solr/
2. Unzip/untar the file. Move to your desired location. Now you have Solr `v.7.0.0`
3. Go into the directory you just created: `cd solr-7.0.0`
4. Launch Solr: `bin/solr start -e cloud -noprompt` - Sets up SolrCloud mode, rather
than Standalone mode. As far as I can tell, SolrCloud mode seems more common.
5. Once Step 4 completes, you can go to `http://localhost:8983/solr/` now, which is
the admin interface for Solr.
6. Load some documents: `bin/post -c gettingstarted docs/`
7. Once Step 6 is complete (will take a few minutes), navigate in your browser to `http://localhost:8983/solr/gettingstarted/select?q=*:*&wt=json` and you should see a
bunch of documents


### Linux

> You should be able to use the above instructions for OSX on a Linux machine.

#### Linuxbrew

Linuxbrew (http://linuxbrew.sh/) is a port of Mac OS homebrew to linux.  Operation is essentially the same as for homebrew.  Follow the installation instructions for linuxbrew and then the instructions for using homebrew (above) should work without modification.

### Windows

You should be able to use the above instructions for OSX on a Windows machine, but with some slight differences. For example, the `bin/post` tool for OSX and Linux doesn't work on Windows, but see https://lucene.apache.org/solr/guide/8_2/post-tool.html#PostTool-Windows for an equivalent.

### `solrium` usage

First, setup a connection object


```r
(conn <- SolrClient$new())
```

```
## <Solr Client>
##   host: 127.0.0.1
##   path: 
##   port: 8983
##   scheme: http
##   errors: simple
##   proxy:
```

And we can now use the `solrium` R package to query the Solr database to get raw JSON data:


```r
conn$search("gettingstarted", params = list(q = '*:*', rows = 3), raw = TRUE)
#> [1] "{\"responseHeader\":{\"status\":0,\"QTime\":8,\"params\":{\"q\":\"*:*\",\"rows\":\"3\",\"wt\":\"json\"}},\"response\":{\"numFound\":3577,\"start\":0,\"maxScore\":1.0,\"docs\":[{\"id\":\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmenter.html\",\"stream_size\":[9016],\"date\":[\"2015-06-10T00:00:00Z\"],\"x_parsed_by\":[\"org.apache.tika.parser.DefaultParser\",\"org.apache.tika.parser.html.HtmlParser\"],\"stream_content_type\":[\"text/html\"],\"dc_title\":[\"Uses of Interface org.apache.solr.highlight.SolrFragmenter (Solr 5.2.1 API)\"],\"content_encoding\":[\"UTF-8\"],\"resourcename\":[\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmenter.html\"],\"title\":[\"Uses of Interface org.apache.solr.highlight.SolrFragmenter (Solr 5.2.1 API)\"],\"content_type\":[\"text/html\"],\"_version_\":1507965023127863296},{\"id\":\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmentsBuilder.html\",\"stream_size\":[10336],\"date\":[\"2015-06-10T00:00:00Z\"],\"x_parsed_by\":[\"org.apache.tika.parser.DefaultParser\",\"org.apache.tika.parser.html.HtmlParser\"],\"stream_content_type\":[\"text/html\"],\"dc_title\":[\"Uses of Class org.apache.solr.highlight.SolrFragmentsBuilder (Solr 5.2.1 API)\"],\"content_encoding\":[\"UTF-8\"],\"resourcename\":[\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmentsBuilder.html\"],\"title\":[\"Uses of Class org.apache.solr.highlight.SolrFragmentsBuilder (Solr 5.2.1 API)\"],\"content_type\":[\"text/html\"],\"_version_\":1507965023153029120},{\"id\":\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/internal/csv/CSVParser.html\",\"stream_size\":[32427],\"date\":[\"2015-06-10T00:00:00Z\"],\"x_parsed_by\":[\"org.apache.tika.parser.DefaultParser\",\"org.apache.tika.parser.html.HtmlParser\"],\"stream_content_type\":[\"text/html\"],\"dc_title\":[\"CSVParser (Solr 5.2.1 API)\"],\"content_encoding\":[\"UTF-8\"],\"resourcename\":[\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/internal/csv/CSVParser.html\"],\"title\":[\"CSVParser (Solr 5.2.1 API)\"],\"content_type\":[\"text/html\"],\"_version_\":1507965023221186560}]}}\n"
#> attr(,"class")
#> [1] "sr_search"
#> attr(,"wt")
#> [1] "json"
```

Or parsed data to a data.frame (just looking at a few columns for brevity):


```r
conn$search("gettingstarted", params = list(q = '*:*', fl = c('date', 'title')))
#> Source: local data frame [10 x 2]
#>
#>                    date                                                                         title
#> 1  2015-06-10T00:00:00Z   Uses of Interface org.apache.solr.highlight.SolrFragmenter (Solr 5.2.1 API)
#> 2  2015-06-10T00:00:00Z Uses of Class org.apache.solr.highlight.SolrFragmentsBuilder (Solr 5.2.1 API)
#> 3  2015-06-10T00:00:00Z                                                    CSVParser (Solr 5.2.1 API)
#> 4  2015-06-10T00:00:00Z                                                     CSVUtils (Solr 5.2.1 API)
#> 5  2015-06-10T00:00:00Z                                 org.apache.solr.internal.csv (Solr 5.2.1 API)
#> 6  2015-06-10T00:00:00Z                 org.apache.solr.internal.csv Class Hierarchy (Solr 5.2.1 API)
#> 7  2015-06-10T00:00:00Z       Uses of Class org.apache.solr.internal.csv.CSVStrategy (Solr 5.2.1 API)
#> 8  2015-06-10T00:00:00Z          Uses of Class org.apache.solr.internal.csv.CSVUtils (Solr 5.2.1 API)
#> 9  2015-06-10T00:00:00Z                                                    CSVConfig (Solr 5.2.1 API)
#> 10 2015-06-10T00:00:00Z                                             CSVConfigGuesser (Solr 5.2.1 API)
```

## Other Vignettes

See the other vignettes for more thorough examples:

* `Document management`
* `Cores/collections management`
* `Solr Search`
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.5 Patched (2021-03-31 r80136) |
|os       |macOS Big Sur 10.16                         |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-05-13                                  |

# Dependencies

|package    |old   |new    |Δ  |
|:----------|:-----|:------|:--|
|solrium    |1.1.4 |1.2.0  |*  |
|cli        |NA    |2.5.0  |*  |
|crul       |NA    |1.1.0  |*  |
|curl       |NA    |4.3.1  |*  |
|dplyr      |NA    |1.0.6  |*  |
|ellipsis   |NA    |0.3.2  |*  |
|lifecycle  |NA    |1.0.0  |*  |
|pillar     |NA    |1.6.0  |*  |
|rlang      |NA    |0.4.11 |*  |
|tibble     |NA    |3.1.1  |*  |
|tidyselect |NA    |1.1.1  |*  |
|utf8       |NA    |1.2.1  |*  |
|vctrs      |NA    |0.3.8  |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*