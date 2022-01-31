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

*Wow, no problems at all. :)**Wow, no problems at all. :)*solrium
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("solrium")
```

Or development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/solrium")
```

```{r}
library("solrium")
```

## Setup

Use `SolrClient$new()` to initialize your connection. These examples use a remote Solr server, but work on any local Solr server.

```{r}
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))
```

You can also set whether you want simple or detailed error messages (via `errors`), and whether you want URLs used in each function call or not (via `verbose`), and your proxy settings (via `proxy`) if needed. For example:

```{r eval=FALSE}
SolrClient$new(errors = "complete")
```

Your settings are printed in the print method for the connection object

```{r}
cli
```

For local Solr server setup:

```
bin/solr start -e cloud -noprompt
bin/post -c gettingstarted example/exampledocs/*.xml
```


## Search

```{r}
(res <- cli$search(params = list(q='*:*', rows=2, fl='id')))
```

And you can get search metadata from the attributes:

```{r}
attributes(res)
```


### Search grouped data

Most recent publication by journal

```{r}
cli$group(params = list(q='*:*', group.field='journal', rows=5, group.limit=1,
                        group.sort='publication_date desc',
                        fl='publication_date, score'))
```

First publication by journal

```{r}
cli$group(params = list(q = '*:*', group.field = 'journal', group.limit = 1,
                        group.sort = 'publication_date asc',
                        fl = c('publication_date', 'score'),
                        fq = "publication_date:[1900-01-01T00:00:00Z TO *]"))
```

Search group query : Last 3 publications of 2013.

```{r}
gq <- 'publication_date:[2013-01-01T00:00:00Z TO 2013-12-31T00:00:00Z]'
cli$group(
  params = list(q='*:*', group.query = gq,
                group.limit = 3, group.sort = 'publication_date desc',
                fl = 'publication_date'))
```

Search group with format simple

```{r}
cli$group(params = list(q='*:*', group.field='journal', rows=5,
                        group.limit=3, group.sort='publication_date desc',
                        group.format='simple', fl='journal, publication_date'))
```

### Facet

```{r}
cli$facet(params = list(q='*:*', facet.field='journal', facet.query=c('cell', 'bird')))
```

### Highlight

```{r}
cli$highlight(params = list(q='alcohol', hl.fl = 'abstract', rows=2))
```

### Stats

```{r}
out <- cli$stats(params = list(q='ecology', stats.field=c('counter_total_all','alm_twitterCount'), stats.facet='journal'))
```

```{r}
out$data
```

### More like this

`solr_mlt` is a function to return similar documents to the one

```{r}
out <- cli$mlt(params = list(q='title:"ecology" AND body:"cell"', mlt.fl='title', mlt.mindf=1, mlt.mintf=1, fl='counter_total_all', rows=5))
```

```{r}
out$docs
```

```{r}
out$mlt
```

### Parsing

`solr_parse` is a general purpose parser function with extension methods `solr_parse.sr_search`, `solr_parse.sr_facet`, and `solr_parse.sr_high`, for parsing `solr_search`, `solr_facet`, and `solr_highlight` function output, respectively. `solr_parse` is used internally within those three functions (`solr_search`, `solr_facet`, `solr_highlight`) to do parsing. You can optionally get back raw `json` or `xml` from `solr_search`, `solr_facet`, and `solr_highlight` setting parameter `raw=TRUE`, and then parsing after the fact with `solr_parse`. All you need to know is `solr_parse` can parse

For example:

```{r}
(out <- cli$highlight(params = list(q='alcohol', hl.fl = 'abstract', rows=2),
                      raw=TRUE))
```

Then parse

```{r}
solr_parse(out, 'df')
```

### Progress bars

only supported in the core search methods: `search`, `facet`, `group`, `mlt`, `stats`, `high`, `all`

```{r eval = FALSE}
library(httr)
invisible(cli$search(params = list(q='*:*', rows=100, fl='id'), progress = httr::progress()))
|==============================================| 100%
```

### Advanced: Function Queries

Function Queries allow you to query on actual numeric fields in the SOLR database, and do addition, multiplication, etc on one or many fields to sort results. For example, here, we search on the product of counter_total_all and alm_twitterCount, using a new temporary field "_val_"

```{r}
cli$search(params = list(q='_val_:"product(counter_total_all,alm_twitterCount)"',
  rows=5, fl='id,title', fq='doc_type:full'))
```

Here, we search for the papers with the most citations

```{r}
cli$search(params = list(q='_val_:"max(counter_total_all)"',
    rows=5, fl='id,counter_total_all', fq='doc_type:full'))
```

Or with the most tweets

```{r}
cli$search(params = list(q='_val_:"max(alm_twitterCount)"',
    rows=5, fl='id,alm_twitterCount', fq='doc_type:full'))
```

### Using specific data sources

__USGS BISON service__

The occurrences service

```{r}
conn <- SolrClient$new(scheme = "https", host = "bison.usgs.gov", path = "solr/occurrences/select", port = NULL)
conn$search(params = list(q = '*:*', fl = c('decimalLatitude','decimalLongitude','scientificName'), rows = 2))
```

The species names service

```{r}
conn <- SolrClient$new(scheme = "https", host = "bison.usgs.gov", path = "solr/scientificName/select", port = NULL)
conn$search(params = list(q = '*:*'))
```

__PLOS Search API__

Most of the examples above use the PLOS search API... :)

## Solr server management

This isn't as complete as searching functions show above, but we're getting there.

### Cores

```{r eval=FALSE}
conn <- SolrClient$new()
```

Many functions, e.g.:

* `core_create()`
* `core_rename()`
* `core_status()`
* ...

Create a core

```{r eval=FALSE}
conn$core_create(name = "foo_bar")
```

### Collections

Many functions, e.g.:

* `collection_create()`
* `collection_list()`
* `collection_addrole()`
* ...

Create a collection

```{r eval=FALSE}
conn$collection_create(name = "hello_world")
```

### Add documents

Add documents, supports adding from files (json, xml, or csv format), and from R objects (including `data.frame` and `list` types so far)

```{r eval=FALSE}
df <- data.frame(id = c(67, 68), price = c(1000, 500000000))
conn$add(df, name = "books")
```

Delete documents, by id

```{r eval=FALSE}
conn$delete_by_id(name = "books", ids = c(3, 4))
```

Or by query

```{r eval=FALSE}
conn$delete_by_query(name = "books", query = "manu:bank")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/solrium/issues)
* License: MIT
* Get citation information for `solrium` in R doing `citation(package = 'solrium')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
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

```{r}
(conn <- SolrClient$new())
```

And we can now use the `solrium` R package to query the Solr database to get raw JSON data:

```{r eval=FALSE}
conn$search("gettingstarted", params = list(q = '*:*', rows = 3), raw = TRUE)
#> [1] "{\"responseHeader\":{\"status\":0,\"QTime\":8,\"params\":{\"q\":\"*:*\",\"rows\":\"3\",\"wt\":\"json\"}},\"response\":{\"numFound\":3577,\"start\":0,\"maxScore\":1.0,\"docs\":[{\"id\":\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmenter.html\",\"stream_size\":[9016],\"date\":[\"2015-06-10T00:00:00Z\"],\"x_parsed_by\":[\"org.apache.tika.parser.DefaultParser\",\"org.apache.tika.parser.html.HtmlParser\"],\"stream_content_type\":[\"text/html\"],\"dc_title\":[\"Uses of Interface org.apache.solr.highlight.SolrFragmenter (Solr 5.2.1 API)\"],\"content_encoding\":[\"UTF-8\"],\"resourcename\":[\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmenter.html\"],\"title\":[\"Uses of Interface org.apache.solr.highlight.SolrFragmenter (Solr 5.2.1 API)\"],\"content_type\":[\"text/html\"],\"_version_\":1507965023127863296},{\"id\":\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmentsBuilder.html\",\"stream_size\":[10336],\"date\":[\"2015-06-10T00:00:00Z\"],\"x_parsed_by\":[\"org.apache.tika.parser.DefaultParser\",\"org.apache.tika.parser.html.HtmlParser\"],\"stream_content_type\":[\"text/html\"],\"dc_title\":[\"Uses of Class org.apache.solr.highlight.SolrFragmentsBuilder (Solr 5.2.1 API)\"],\"content_encoding\":[\"UTF-8\"],\"resourcename\":[\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/highlight/class-use/SolrFragmentsBuilder.html\"],\"title\":[\"Uses of Class org.apache.solr.highlight.SolrFragmentsBuilder (Solr 5.2.1 API)\"],\"content_type\":[\"text/html\"],\"_version_\":1507965023153029120},{\"id\":\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/internal/csv/CSVParser.html\",\"stream_size\":[32427],\"date\":[\"2015-06-10T00:00:00Z\"],\"x_parsed_by\":[\"org.apache.tika.parser.DefaultParser\",\"org.apache.tika.parser.html.HtmlParser\"],\"stream_content_type\":[\"text/html\"],\"dc_title\":[\"CSVParser (Solr 5.2.1 API)\"],\"content_encoding\":[\"UTF-8\"],\"resourcename\":[\"/Users/sacmac/solr-5.2.1/docs/solr-core/org/apache/solr/internal/csv/CSVParser.html\"],\"title\":[\"CSVParser (Solr 5.2.1 API)\"],\"content_type\":[\"text/html\"],\"_version_\":1507965023221186560}]}}\n"
#> attr(,"class")
#> [1] "sr_search"
#> attr(,"wt")
#> [1] "json"
```

Or parsed data to a data.frame (just looking at a few columns for brevity):

```{r eval=FALSE}
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
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Cores/collections management}
%\VignetteEncoding{UTF-8}
-->

```{r, echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

Cores/collections management
============================

## Installation

Stable version from CRAN

```{r eval=FALSE}
install.packages("solrium")
```

Or the development version from GitHub

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/solrium")
```

Load

```{r}
library("solrium")
```

Initialize connection

```{r}
(conn <- SolrClient$new())
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

```{r eval=FALSE}
conn$core_create()
```

### Delete a core

```{r eval=FALSE}
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

```{r eval=FALSE}
conn$collection_create()
```

### Delete a collection

```{r eval=FALSE}
conn$collection_delete()
```
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Document management}
%\VignetteEncoding{UTF-8}
-->

```{r, echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

Document management
===================

## Installation

Stable version from CRAN

```{r eval=FALSE}
install.packages("solrium")
```

Or the development version from GitHub

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/solrium")
```

Load

```{r}
library("solrium")
```

Initialize connection. By default, you connect to `http://localhost:8983`

```{r}
(conn <- SolrClient$new())
```

## Create documents from R objects

For now, only lists and data.frame's supported.

### data.frame

```{r}
df <- data.frame(id = c(67, 68), price = c(1000, 500000000))
conn$add(df, "books")
```

### list

```{r echo=FALSE, results='hide'}
conn$delete_by_id(1:2, "books")
```

```{r}
ss <- list(list(id = 1, price = 100), list(id = 2, price = 500))
conn$add(ss, "books")
```

## Delete documents

### By id

Add some documents first

```{r echo=FALSE, results='hide'}
conn$delete_by_id(1:3, "gettingstarted")
```

```{r}
docs <- list(list(id = 1, price = 100, name = "brown"),
             list(id = 2, price = 500, name = "blue"),
             list(id = 3, price = 2000L, name = "pink"))
conn$add(docs, "gettingstarted")
```

And the documents are now in your Solr database

```{r}
conn$search(name = "gettingstarted", params = list(q = "*:*", rows = 3))
```

Now delete those documents just added

```{r deleteid}
conn$delete_by_id(ids = c(1, 2, 3), "gettingstarted")
```

And now they are gone

```{r}
conn$search("gettingstarted", params = list(q = "*:*", rows = 4))
```

### By query

Add some documents first

```{r}
conn$add(docs, "gettingstarted")
```

And the documents are now in your Solr database

```{r}
conn$search("gettingstarted", params = list(q = "*:*", rows = 5))
```

Now delete those documents just added

```{r deletequery}
conn$delete_by_query(query = "(name:blue OR name:pink)", "gettingstarted")
```

And now they are gone

```{r}
conn$search("gettingstarted", params = list(q = "*:*", rows = 5))
```

## Update documents from files

This approach is best if you have many different things you want to do at once, e.g., delete and add files and set any additional options. The functions are:

* `update_xml()`
* `update_json()`
* `update_csv()`

There are separate functions for each of the data types as they take slightly different parameters - and to make it more clear that those are the three input options for data types.

### JSON

```{r}
file <- system.file("examples", "books.json", package = "solrium")
conn$update_json(file, "books")
```

### Add and delete in the same file

Add a document first, that we can later delete

```{r}
ss <- list(list(id = 456, name = "cat"))
conn$add(ss, "books")
```

Now add a new document, and delete the one we just made

```{r}
file <- system.file("examples", "add_delete.xml", package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_xml(file, "books")
```

### Notes

Note that `update_xml()` and `update_json()` have exactly the same parameters, but simply use different data input formats. `update_csv()` is different in that you can't provide document or field level boosts or other modifications. In addition `update_csv()` can accept not just csv, but tsv and other types of separators.

<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Solr search}
%\VignetteEncoding{UTF-8}
-->

```{r, echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("solrium")
```

Or the development version from GitHub

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/solrium")
```

Load

```{r}
library("solrium")
```

## Setup connection

You can setup for a remote Solr instance or on your local machine.

```{r}
(conn <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))
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

```{r}
conn$search(params = list(q = '*:*', rows = 2, fl = 'id'))
```

__Search in specific fields with `:`__

Search for word ecology in title and cell in the body

```{r}
conn$search(params = list(q = 'title:"ecology" AND body:"cell"', fl = 'title', rows = 5))
```

__Wildcards__

Search for word that starts with "cell" in the title field

```{r}
conn$search(params = list(q = 'title:"cell*"', fl = 'title', rows = 5))
```

__Proximity search__

Search for words "sports" and "alcohol" within four words of each other

```{r}
conn$search(params = list(q = 'everything:"stem cell"~7', fl = 'title', rows = 3))
```

__Range searches__

Search for articles with Twitter count between 5 and 10

```{r}
conn$search(params = list(q = '*:*', fl = c('alm_twitterCount', 'id'), fq = 'alm_twitterCount:[5 TO 50]', rows = 10))
```

__Boosts__

Assign higher boost to title matches than to body matches (compare the two calls)

```{r}
conn$search(params = list(q = 'title:"cell" abstract:"science"', fl = 'title', rows = 3))
```

```{r}
conn$search(params = list(q = 'title:"cell"^1.5 AND abstract:"science"', fl = 'title', rows = 3))
```

## Search all

`solr_all()` differs from `solr_search()` in that it allows specifying facets, mlt, groups,
stats, etc, and returns all of those. It defaults to `parsetype = "list"` and `wt="json"`,
whereas `solr_search()` defaults to `parsetype = "df"` and `wt="csv"`. `solr_all()` returns
by default a list, whereas `solr_search()` by default returns a data.frame.

A basic search, just docs output

```{r}
conn$all(params = list(q = '*:*', rows = 2, fl = 'id'))
```

Get docs, mlt, and stats output

```{r}
conn$all(params = list(q = 'ecology', rows = 2, fl = 'id', mlt = 'true', mlt.count = 2, mlt.fl = 'abstract', stats = 'true', stats.field = 'counter_total_all'))
```


## Facet

```{r}
conn$facet(params = list(q = '*:*', facet.field = 'journal', facet.query = c('cell', 'bird')))
```

## Highlight

```{r}
conn$highlight(params = list(q = 'alcohol', hl.fl = 'abstract', rows = 2))
```

## Stats

```{r}
out <- conn$stats(params = list(q = 'ecology', stats.field = c('counter_total_all', 'alm_twitterCount'), stats.facet = c('journal', 'volume')))
```

```{r}
out$data
```

```{r}
out$facet
```

## More like this

`solr_mlt` is a function to return similar documents to the one

```{r}
out <- conn$mlt(params = list(q = 'title:"ecology" AND body:"cell"', mlt.fl = 'title', mlt.mindf = 1, mlt.mintf = 1, fl = 'counter_total_all', rows = 5))
out$docs
```

```{r}
out$mlt
```

## Groups

`solr_groups()` is a function to return similar documents to the one

```{r}
conn$group(params = list(q = 'ecology', group.field = 'journal', group.limit = 1, fl = c('id', 'alm_twitterCount')))
```

## Parsing

`solr_parse()` is a general purpose parser function with extension methods for parsing outputs from functions in `solr`. `solr_parse()` is used internally within functions to do parsing after retrieving data from the server. You can optionally get back raw `json`, `xml`, or `csv` with the `raw=TRUE`, and then parse afterwards with `solr_parse()`.

For example:

```{r}
(out <- conn$highlight(params = list(q = 'alcohol', hl.fl = 'abstract', rows = 2), raw = TRUE))
```

Then parse

```{r}
solr_parse(out, 'df')
```

[Please report any issues or bugs](https://github.com/ropensci/solrium/issues).
---
title: Local Solr setup
author: Scott Chamberlain
date: "2020-04-22"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Local Solr setup}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

<!--
to run this vignette:
docker run -d -p 8983:8983 --name my_solr solr:latest solr-precreate gettingstarted
docker exec  -it <container id> /bin/bash # go into container
bin/post -c gettingstarted example/exampledocs/ 
-->

The Solr version you are working with my differ from below. Don't worry, just go with the version you want to use.

### OSX

1. Download most recent version from an Apache mirror https://solr.apache.org/
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

### Windows

You should be able to use the above instructions for OSX on a Windows machine, but with some slight differences. For example, the `bin/post` tool for OSX and Linux doesn't work on Windows, but see https://solr.apache.org/guide/8_2/post-tool.html#PostTool-Windows for an equivalent.

### `solrium` usage

First, setup a connection object


```r
library(solrium)
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
```

```
## [1] "{\n  \"responseHeader\":{\n    \"status\":0,\n    \"QTime\":5,\n    \"params\":{\n      \"q\":\"*:*\",\n      \"rows\":\"3\",\n      \"wt\":\"json\"}},\n  \"response\":{\"numFound\":48,\"start\":0,\"docs\":[\n      {\n        \"id\":\"TWINX2048-3200PRO\",\n        \"name\":[\"CORSAIR  XMS 2GB (2 x 1GB) 184-Pin DDR SDRAM Unbuffered DDR 400 (PC 3200) Dual Channel Kit System Memory - Retail\"],\n        \"manu\":[\"Corsair Microsystems Inc.\"],\n        \"manu_id_s\":\"corsair\",\n        \"cat\":[\"electronics\",\n          \"memory\"],\n        \"features\":[\"CAS latency 2,  2-3-3-6 timing, 2.75v, unbuffered, heat-spreader\"],\n        \"price\":[185.0],\n        \"popularity\":[5],\n        \"inStock\":[true],\n        \"store\":[\"37.7752,-122.4232\"],\n        \"manufacturedate_dt\":\"2006-02-13T15:26:37Z\",\n        \"payloads\":[\"electronics|6.0 memory|3.0\"],\n        \"_version_\":1664688417970061312},\n      {\n        \"id\":\"VS1GB400C3\",\n        \"name\":[\"CORSAIR ValueSelect 1GB 184-Pin DDR SDRAM Unbuffered DDR 400 (PC 3200) System Memory - Retail\"],\n        \"manu\":[\"Corsair Microsystems Inc.\"],\n        \"manu_id_s\":\"corsair\",\n        \"cat\":[\"electronics\",\n          \"memory\"],\n        \"price\":[74.99],\n        \"popularity\":[7],\n        \"inStock\":[true],\n        \"store\":[\"37.7752,-100.0232\"],\n        \"manufacturedate_dt\":\"2006-02-13T15:26:37Z\",\n        \"payloads\":[\"electronics|4.0 memory|2.0\"],\n        \"_version_\":1664688418109521920},\n      {\n        \"id\":\"VDBDB1A16\",\n        \"name\":[\"A-DATA V-Series 1GB 184-Pin DDR SDRAM Unbuffered DDR 400 (PC 3200) System Memory - OEM\"],\n        \"manu\":[\"A-DATA Technology Inc.\"],\n        \"manu_id_s\":\"corsair\",\n        \"cat\":[\"electronics\",\n          \"memory\"],\n        \"features\":[\"CAS latency 3,   2.7v\"],\n        \"popularity\":[0],\n        \"inStock\":[true],\n        \"store\":[\"45.18414,-93.88141\"],\n        \"manufacturedate_dt\":\"2006-02-13T15:26:37Z\",\n        \"payloads\":[\"electronics|0.9 memory|0.1\"],\n        \"_version_\":1664688418113716224}]\n  }}\n"
## attr(,"class")
## [1] "sr_search"
## attr(,"wt")
## [1] "json"
```

Or parsed data to a data.frame (just looking at a few columns for brevity):


```r
conn$search("gettingstarted", params = list(q = '*:*', fl = c('id', 'address_s')))
```

```
## # A tibble: 10 x 2
##    id                address_s                                                 
##    <chr>             <chr>                                                     
##  1 TWINX2048-3200PRO <NA>                                                      
##  2 VS1GB400C3        <NA>                                                      
##  3 VDBDB1A16         <NA>                                                      
##  4 SOLR1000          <NA>                                                      
##  5 adata             46221 Landing Parkway Fremont, CA 94538                   
##  6 apple             1 Infinite Way, Cupertino CA                              
##  7 asus              800 Corporate Way Fremont, CA 94539                       
##  8 ati               33 Commerce Valley Drive East Thornhill, ON L3T 7N6 Canada
##  9 belkin            12045 E. Waterfront Drive Playa Vista, CA 90094           
## 10 canon             One Canon Plaza Lake Success, NY 11042
```

## Other Vignettes

See the other vignettes for more thorough examples:

* `Document management`
* `Cores/collections management`
* `Solr Search`
---
title: Cores/collections management
author: Scott Chamberlain
date: "2020-04-22"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Cores/collections management}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---




## Installation

Stable version from CRAN


```r
install.packages("solrium")
```

Or the development version from GitHub


```r
remotes::install_github("ropensci/solrium")
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
---
title: Document management
author: Scott Chamberlain
date: "2020-04-22"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Document management}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



<!--
to run this vignette:
cd vignettes/docker
docker-compose up -d
-->

## Installation

Stable version from CRAN


```r
install.packages("solrium")
```

Or the development version from GitHub


```r
remotes::install_github("ropensci/solrium")
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
if (!collection_exists(conn, "books")) {
  collection_create(conn, name = "books", numShards = 1)
}
```

```
#> $responseHeader
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 6513
#> 
#> 
#> $success
#> $success$`172.20.0.5:8983_solr`
#> $success$`172.20.0.5:8983_solr`$responseHeader
#> $success$`172.20.0.5:8983_solr`$responseHeader$status
#> [1] 0
#> 
#> $success$`172.20.0.5:8983_solr`$responseHeader$QTime
#> [1] 5024
#> 
#> 
#> $success$`172.20.0.5:8983_solr`$core
#> [1] "books_shard1_replica_n1"
#> 
#> 
#> 
#> $warning
#> [1] "Using _default configset. Data driven schema functionality is enabled by default, which is NOT RECOMMENDED for production use. To turn it off: curl http://{host:port}/solr/books/config -d '{\"set-user-property\": {\"update.autoCreateFields\":\"false\"}}'"
```

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
#> [1] 987
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
#> [1] 61
```

## Delete documents

### By id

Create collection if it doesn't exist yet


```r
if (!collection_exists(conn, "gettingstarted")) {
  collection_create(conn, name = "gettingstarted", numShards = 1)
}
```

```
#> $responseHeader
#> $responseHeader$status
#> [1] 0
#> 
#> $responseHeader$QTime
#> [1] 6446
#> 
#> 
#> $success
#> $success$`172.20.0.7:8983_solr`
#> $success$`172.20.0.7:8983_solr`$responseHeader
#> $success$`172.20.0.7:8983_solr`$responseHeader$status
#> [1] 0
#> 
#> $success$`172.20.0.7:8983_solr`$responseHeader$QTime
#> [1] 5112
#> 
#> 
#> $success$`172.20.0.7:8983_solr`$core
#> [1] "gettingstarted_shard1_replica_n1"
#> 
#> 
#> 
#> $warning
#> [1] "Using _default configset. Data driven schema functionality is enabled by default, which is NOT RECOMMENDED for production use. To turn it off: curl http://{host:port}/solr/gettingstarted/config -d '{\"set-user-property\": {\"update.autoCreateFields\":\"false\"}}'"
```

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
#> [1] 1108
```

And the documents are now in your Solr database


```r
conn$search(name = "gettingstarted", params = list(q = "*:*", rows = 3))
```

```
#> # A tibble: 3 x 4
#>   id    price name  `_version_`
#>   <chr> <int> <chr>       <dbl>
#> 1 1       100 brown     1.66e18
#> 2 2       500 blue      1.66e18
#> 3 3      2000 pink      1.66e18
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
#> [1] 48
```

And now they are gone


```r
conn$search("gettingstarted", params = list(q = "*:*", rows = 4))
```

```
#> # A tibble: 0 x 0
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
#> [1] 72
```

And the documents are now in your Solr database


```r
conn$search("gettingstarted", params = list(q = "*:*", rows = 5))
```

```
#> # A tibble: 3 x 4
#>   id    price name  `_version_`
#>   <chr> <int> <chr>       <dbl>
#> 1 1       100 brown     1.66e18
#> 2 2       500 blue      1.66e18
#> 3 3      2000 pink      1.66e18
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
#> [1] 122
```

And now they are gone


```r
conn$search("gettingstarted", params = list(q = "*:*", rows = 5))
```

```
#> # A tibble: 1 x 4
#>   id    price name  `_version_`
#>   <chr> <int> <chr>       <dbl>
#> 1 1       100 brown     1.66e18
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
#> [1] 782
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
#> [1] 100
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
#> [1] 279
```

### Notes

Note that `update_xml()` and `update_json()` have exactly the same parameters, but simply use different data input formats. `update_csv()` is different in that you can't provide document or field level boosts or other modifications. In addition `update_csv()` can accept not just csv, but tsv and other types of separators.

---
title: solr search
author: Scott Chamberlain
date: "2020-04-22"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{solr search}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



**A general purpose R interface to [Apache Solr](https://solr.apache.org/)**


## Installation

Stable version from CRAN


```r
install.packages("solrium")
```

Or the development version from GitHub


```r
remotes::install_github("ropensci/solrium")
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
#> 1 10.1371/journal.pone.0020843
#> 2 10.1371/journal.pone.0022257
```

__Search in specific fields with `:`__

Search for word ecology in title and cell in the body


```r
conn$search(params = list(q = 'title:"ecology" AND body:"cell"', fl = 'title', rows = 5))
#> # A tibble: 5 x 1
#>   title                                                    
#>   <chr>                                                    
#> 1 The Ecology of Collective Behavior                       
#> 2 Chasing Ecological Interactions                          
#> 3 Ecology's Big, Hot Idea                                  
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
#> 1 Altered Development of NKT Cells, γδ T Cells, CD8 T Cells and NK Cells in a P…
#> 2 Cancer Stem Cell-Like Side Population Cells in Clear Cell Renal Cell Carcinom…
#> 3 Cell-Cell Contact Preserves Cell Viability via Plakoglobin                    
#> 4 Transgelin-2 in B-Cells Controls T-Cell Activation by Stabilizing T Cell - B …
#> 5 Tetherin Can Restrict Cell-Free and Cell-Cell Transmission of HIV from Primar…
```

__Proximity search__

Search for words "sports" and "alcohol" within four words of each other


```r
conn$search(params = list(q = 'everything:"stem cell"~7', fl = 'title', rows = 3))
#> # A tibble: 3 x 1
#>   title                                                                         
#>   <chr>                                                                         
#> 1 Effect of Dedifferentiation on Time to Mutation Acquisition in Stem Cell-Driv…
#> 2 Symmetric vs. Asymmetric Stem Cell Divisions: An Adaptation against Cancer?   
#> 3 Phenotypic Evolutionary Models in Stem Cell Biology: Replacement, Quiescence,…
```

__Range searches__

Search for articles with Twitter count between 5 and 10


```r
conn$search(params = list(q = '*:*', fl = c('alm_twitterCount', 'id'), fq = 'alm_twitterCount:[5 TO 50]', rows = 10))
#> # A tibble: 10 x 2
#>    id                           alm_twitterCount
#>    <chr>                                   <int>
#>  1 10.1371/journal.pone.0020913               14
#>  2 10.1371/journal.pone.0023176               12
#>  3 10.1371/journal.pone.0023286               16
#>  4 10.1371/journal.pone.0025390               32
#>  5 10.1371/journal.pone.0013199               24
#>  6 10.1371/journal.pone.0013196                5
#>  7 10.1371/journal.pone.0009019               11
#>  8 10.1371/journal.pone.0008993                6
#>  9 10.1371/journal.pone.0021336               18
#> 10 10.1371/journal.pone.0009042               18
```

__Boosts__

Assign higher boost to title matches than to body matches (compare the two calls)


```r
conn$search(params = list(q = 'title:"cell" abstract:"science"', fl = 'title', rows = 3))
#> # A tibble: 3 x 1
#>   title                                                                         
#>   <chr>                                                                         
#> 1 I Want More and Better Cells! – An Outreach Project about Stem Cells and Its …
#> 2 Virtual Reconstruction and Three-Dimensional Printing of Blood Cells as a Too…
#> 3 Globalization of Stem Cell Science: An Examination of Current and Past Collab…
```


```r
conn$search(params = list(q = 'title:"cell"^1.5 AND abstract:"science"', fl = 'title', rows = 3))
#> # A tibble: 3 x 1
#>   title                                                                         
#>   <chr>                                                                         
#> 1 I Want More and Better Cells! – An Outreach Project about Stem Cells and Its …
#> 2 Virtual Reconstruction and Three-Dimensional Printing of Blood Cells as a Too…
#> 3 Globalization of Stem Cell Science: An Examination of Current and Past Collab…
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
#> 1 10.1371/journal.pone.0020843
#> 2 10.1371/journal.pone.0022257
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
#> 1 10.1371/journal.pone.0020843
#> 2 10.1371/journal.pone.0022257
#> 
#> $mlt$mlt
#> list()
#> 
#> 
#> $group
#>   numFound start                           id
#> 1  2349539     0 10.1371/journal.pone.0020843
#> 2  2349539     0 10.1371/journal.pone.0022257
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
#> 1   245661     0 10.1371/journal.pone.0001275
#> 2   245661     0 10.1371/journal.pbio.1002448
#> 
#> $mlt$mlt$`10.1371/journal.pone.0059813`
#> # A tibble: 2 x 3
#>   numFound start id                          
#>      <int> <int> <chr>                       
#> 1   237548     0 10.1371/journal.pone.0210707
#> 2   237548     0 10.1371/journal.pone.0204749
#> 
#> 
#> 
#> $group
#>   numFound start                           id
#> 1    51895     0 10.1371/journal.pone.0001248
#> 2    51895     0 10.1371/journal.pone.0059813
#> 
#> $stats
#> $stats$data
#>                   min     max count missing       sum sumOfSquares     mean
#> counter_total_all   0 1532460 51895       0 324414287 1.533069e+13 6251.359
#>                     stddev
#> counter_total_all 16010.71
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
#> 1 cell  186815
#> 2 bird   20072
#> 
#> $facet_fields
#> $facet_fields$journal
#> # A tibble: 9 x 2
#>   term                             value  
#>   <fct>                            <fct>  
#> 1 plos one                         1946213
#> 2 plos genetics                    72617  
#> 3 plos pathogens                   66100  
#> 4 plos neglected tropical diseases 65251  
#> 5 plos computational biology       59985  
#> 6 plos biology                     42060  
#> 7 plos medicine                    29915  
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
#> 1 10.1371/journal.pone… Background: Binge drinking, an increasingly common form…
#> 2 10.1371/journal.pone… Background and Aim: Harmful <em>alcohol</em> consumptio…
```

## Stats


```r
out <- conn$stats(params = list(q = 'ecology', stats.field = c('counter_total_all', 'alm_twitterCount'), stats.facet = c('journal', 'volume')))
```


```r
out$data
#>                   min     max count missing       sum sumOfSquares        mean
#> counter_total_all   0 1532460 51895       0 324414287 1.533069e+13 6251.359225
#> alm_twitterCount    0    3439 51895       0    308886 8.194604e+07    5.952134
#>                        stddev
#> counter_total_all 16010.71374
#> alm_twitterCount     39.28964
```


```r
out$facet
#> $counter_total_all
#> $counter_total_all$volume
#> # A tibble: 18 x 9
#>    volume   min     max count missing      sum  sumOfSquares   mean stddev
#>    <chr>  <dbl>   <dbl> <int>   <int>    <dbl>         <dbl>  <dbl>  <dbl>
#>  1 11         0  300326  5264       0 26131455  498787135805  4964.  8374.
#>  2 12         0  637916  5078       0 23045890  754959012682  4538. 11318.
#>  3 13         0  187650  4772       0 14469852  185727513548  3032.  5453.
#>  4 14         0  202041  4117       0  8965786  135038722582  2178.  5298.
#>  5 15         0   56955  1558       0  3703520   44446597378  2377.  4785.
#>  6 16         0   51981   278       0  1639028   23503754454  5896.  7069.
#>  7 17         0   64350   145       0   957264   12461694744  6602.  6531.
#>  8 18        10   10293    25       0    87761     452935251  3510.  2457.
#>  9 1       2201  394613    81       0  2136280  281140741696 26374. 53009.
#> 10 2          0  160743   483       0  7874872  316645622198 16304. 19763.
#> 11 3          0  140263   743       0 10322036  340678112084 13892. 16306.
#> 12 4          0  421769  1010       0 13320977  619882873519 13189. 20982.
#> 13 5          0  261787  1539       0 17578591  518230843883 11422. 14367.
#> 14 6          0  421511  2949       0 27083643  844617111025  9184. 14217.
#> 15 7          0  337925  4826       0 39377491 1060307559235  8159. 12376.
#> 16 8          0  662892  6363       0 46269702 1567973982796  7272. 13913.
#> 17 9          0 1532460  6622       0 45978003 4777504680537  6943. 25949.
#> 18 10         0 1245813  6042       0 35472136 3348332631936  5871. 22799.
#> 
#> $counter_total_all$journal
#> # A tibble: 9 x 9
#>   journal   min     max count missing       sum  sumOfSquares   mean stddev
#>   <chr>   <dbl>   <dbl> <int>   <int>     <dbl>         <dbl>  <dbl>  <dbl>
#> 1 1        1443  227764  1249       0   9827787  251663612483  7869. 11819.
#> 2 2           0  429864  1329       0  24811193 1527110853031 18669. 28304.
#> 3 3           0  337925   355       0   7588181  509085551899 21375. 31303.
#> 4 4        9559   19071     2       0     28630     455077522 14315   6726.
#> 5 5           0 1245813 43312       0 236610221 9778305707321  5463. 13997.
#> 6 6           0  181861   957       0   9440164  202232267528  9864. 10683.
#> 7 7           0  163569  1153       0  11587785  234170621039 10050. 10108.
#> 8 8           0  343133  2440       0  14360560  292773478918  5885.  9240.
#> 9 9           0 1532460  1098       0  10159766 2534894355612  9253. 47170.
#> 
#> 
#> $alm_twitterCount
#> $alm_twitterCount$volume
#> # A tibble: 18 x 9
#>    volume   min   max count missing   sum sumOfSquares  mean stddev
#>    <chr>  <dbl> <dbl> <int>   <int> <dbl>        <dbl> <dbl>  <dbl>
#>  1 11         0  2142  5264       0 52454     11452818  9.96 45.6  
#>  2 12         0  1890  5078       0 39451     11059295  7.77 46.0  
#>  3 13         0   581  4772       0 17551      2521623  3.68 22.7  
#>  4 14         0   984  4117       0 15341      4027537  3.73 31.1  
#>  5 15         0   453  1558       0  6235       912739  4.00 23.9  
#>  6 16         0   456   278       0  2853       414803 10.3  37.3  
#>  7 17         0    44   145       0   224         5628  1.54  6.06 
#>  8 18         0     3    25       0     5           11  0.2   0.645
#>  9 1          0    47    81       0   208         6306  2.57  8.49 
#> 10 2          0   125   483       0  1119        63677  2.32 11.3  
#> 11 3          0   504   743       0  1408       271870  1.90 19.0  
#> 12 4          0   313  1010       0  1658       155910  1.64 12.3  
#> 13 5          0   165  1539       0  2698       142098  1.75  9.45 
#> 14 6          0   972  2949       0  5733      1647077  1.94 23.6  
#> 15 7          0   864  4826       0 21873      2587889  4.53 22.7  
#> 16 8          0  2029  6363       0 40730     10857672  6.40 40.8  
#> 17 9          0  1880  6622       0 55312     17034646  8.35 50.0  
#> 18 10         0  3439  6042       0 44033     18784443  7.29 55.3  
#> 
#> $alm_twitterCount$journal
#> # A tibble: 9 x 9
#>   journal   min   max count missing    sum sumOfSquares  mean stddev
#>   <chr>   <dbl> <dbl> <int>   <int>  <dbl>        <dbl> <dbl>  <dbl>
#> 1 1           0   456  1249       0   7370       577214  5.90  20.7 
#> 2 2           0  2142  1329       0  39113     15093463 29.4  102.  
#> 3 3           0   832   355       0   5830      1219358 16.4   56.3 
#> 4 4           0     3     2       0      3            9  1.5    2.12
#> 5 5           0  3439 43312       0 216075     61407349  4.99  37.3 
#> 6 6           0   250   957       0   8508       415112  8.89  18.8 
#> 7 7           0   230  1153       0   9273       461699  8.04  18.3 
#> 8 8           0   972  2440       0  11854      1563704  4.86  24.8 
#> 9 9           0   581  1098       0  10860      1208134  9.89  31.7
```

## More like this

`solr_mlt` is a function to return similar documents to the one


```r
out <- conn$mlt(params = list(q = 'title:"ecology" AND body:"cell"', mlt.fl = 'title', mlt.mindf = 1, mlt.mintf = 1, fl = 'counter_total_all', rows = 5))
out$docs
#> # A tibble: 5 x 2
#>   id                           counter_total_all
#>   <chr>                                    <int>
#> 1 10.1371/journal.pbio.1001805             25190
#> 2 10.1371/journal.pbio.1002559             13578
#> 3 10.1371/journal.pbio.0020440             26486
#> 4 10.1371/journal.pone.0087217             19956
#> 5 10.1371/journal.pbio.1002191             30074
```


```r
out$mlt
#> $`10.1371/journal.pbio.1001805`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     4900     0 10.1371/journal.pone.0098876              4302
#> 2     4900     0 10.1371/journal.pone.0082578              3542
#> 3     4900     0 10.1371/journal.pone.0193049              2661
#> 4     4900     0 10.1371/journal.pone.0102159              2687
#> 5     4900     0 10.1371/journal.pcbi.1002652              4687
#> 
#> $`10.1371/journal.pbio.1002559`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     6461     0 10.1371/journal.pone.0155028              4235
#> 2     6461     0 10.1371/journal.pone.0041684             29479
#> 3     6461     0 10.1371/journal.pone.0023086             10177
#> 4     6461     0 10.1371/journal.pone.0155989              3893
#> 5     6461     0 10.1371/journal.pone.0223982               728
#> 
#> $`10.1371/journal.pbio.0020440`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     1428     0 10.1371/journal.pone.0162651              3981
#> 2     1428     0 10.1371/journal.pone.0003259              3543
#> 3     1428     0 10.1371/journal.pone.0102679              5638
#> 4     1428     0 10.1371/journal.pone.0068814             10143
#> 5     1428     0 10.1371/journal.pntd.0003377              4835
#> 
#> $`10.1371/journal.pone.0087217`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1     5802     0 10.1371/journal.pone.0175497              2448
#> 2     5802     0 10.1371/journal.pone.0204743               377
#> 3     5802     0 10.1371/journal.pone.0159131              6664
#> 4     5802     0 10.1371/journal.pone.0220409              1224
#> 5     5802     0 10.1371/journal.pone.0123774              2454
#> 
#> $`10.1371/journal.pbio.1002191`
#> # A tibble: 5 x 4
#>   numFound start id                           counter_total_all
#>      <int> <int> <chr>                                    <int>
#> 1    15002     0 10.1371/journal.pbio.1002232                 0
#> 2    15002     0 10.1371/journal.pone.0131700              3824
#> 3    15002     0 10.1371/journal.pone.0070448              2713
#> 4    15002     0 10.1371/journal.pone.0191705              1971
#> 5    15002     0 10.1371/journal.pone.0160798              4387
```

## Groups

`solr_groups()` is a function to return similar documents to the one


```r
conn$group(params = list(q = 'ecology', group.field = 'journal', group.limit = 1, fl = c('id', 'alm_twitterCount')))
#>                         groupValue numFound start                           id
#> 1                         plos one    43312     0 10.1371/journal.pone.0001248
#> 2       plos computational biology     1098     0 10.1371/journal.pcbi.1003594
#> 3                     plos biology     1329     0 10.1371/journal.pbio.0060300
#> 4                    plos genetics     1153     0 10.1371/journal.pgen.1005860
#> 5 plos neglected tropical diseases     2440     0 10.1371/journal.pntd.0004689
#> 6                   plos pathogens      957     0 10.1371/journal.ppat.1005780
#> 7                             none     1249     0 10.1371/journal.pone.0046761
#> 8                    plos medicine      355     0 10.1371/journal.pmed.1000303
#> 9             plos clinical trials        2     0 10.1371/journal.pctr.0020010
#>   alm_twitterCount
#> 1                0
#> 2               21
#> 3                0
#> 4              135
#> 5               13
#> 6               19
#> 7                0
#> 8                1
#> 9                0
```

## Parsing

`solr_parse()` is a general purpose parser function with extension methods for parsing outputs from functions in `solr`. `solr_parse()` is used internally within functions to do parsing after retrieving data from the server. You can optionally get back raw `json`, `xml`, or `csv` with the `raw=TRUE`, and then parse afterwards with `solr_parse()`.

For example:


```r
(out <- conn$highlight(params = list(q = 'alcohol', hl.fl = 'abstract', rows = 2), raw = TRUE))
#> [1] "{\n  \"response\":{\"numFound\":32898,\"start\":0,\"maxScore\":4.737952,\"docs\":[\n      {\n        \"id\":\"10.1371/journal.pone.0218147\",\n        \"journal\":\"PLOS ONE\",\n        \"eissn\":\"1932-6203\",\n        \"publication_date\":\"2019-12-10T00:00:00Z\",\n        \"article_type\":\"Research Article\",\n        \"author_display\":[\"Victor M. Jimenez Jr.\",\n          \"Erik W. Settles\",\n          \"Bart J. Currie\",\n          \"Paul S. Keim\",\n          \"Fernando P. Monroy\"],\n        \"abstract\":[\"Background: Binge drinking, an increasingly common form of alcohol use disorder, is associated with substantial morbidity and mortality; yet, its effects on the immune system’s ability to defend against infectious agents are poorly understood. Burkholderia pseudomallei, the causative agent of melioidosis can occur in healthy humans, yet binge alcohol intoxication is increasingly being recognized as a major risk factor. Although our previous studies demonstrated that binge alcohol exposure increased B. pseudomallei near-neighbor virulence in vivo and increased paracellular diffusion and intracellular invasion, no experimental studies have examined the extent to which bacterial and alcohol dosage play a role in disease progression. In addition, the temporal effects of a single binge alcohol dose prior to infection has not been examined in vivo. Principal findings: In this study, we used B. thailandensis E264 a close genetic relative of B. pseudomallei, as useful BSL-2 model system. Eight-week-old female C57BL/6 mice were utilized in three distinct animal models to address the effects of 1) bacterial dosage, 2) alcohol dosage, and 3) the temporal effects, of a single binge alcohol episode. Alcohol was administered comparable to human binge drinking (≤ 4.4 g/kg) or PBS intraperitoneally before a non-lethal intranasal infection. Bacterial colonization of lung and spleen was increased in mice administered alcohol even after bacterial dose was decreased 10-fold. Lung and not spleen tissue were colonized even after alcohol dosage was decreased 20 times below the U.S legal limit. Temporally, a single binge alcohol episode affected lung bacterial colonization for more than 24 h after alcohol was no longer detected in the blood. Pulmonary and splenic cytokine expression (TNF-α, GM-CSF) remained suppressed, while IL-12/p40 increased in mice administered alcohol 6 or 24 h prior to infection. Increased lung and not intestinal bacterial invasion was observed in human and murine non-phagocytic epithelial cells exposed to 0.2% v/v alcohol in vitro. Conclusions: Our results indicate that the effects of a single binge alcohol episode are tissue specific. A single binge alcohol intoxication event increases bacterial colonization in mouse lung tissue even after very low BACs and decreases the dose required to colonize the lungs with less virulent B. thailandensis. Additionally, the temporal effects of binge alcohol alters lung and spleen cytokine expression for at least 24 h after alcohol is detected in the blood. Delayed recovery in lung and not spleen tissue may provide a means for B. pseudomallei and near-neighbors to successfully colonize lung tissue through increased intracellular invasion of non-phagocytic cells in patients with hazardous alcohol intake. \"],\n        \"title_display\":\"Persistence of <i>Burkholderia thailandensis</i> E264 in lung tissue after a single binge alcohol episode\",\n        \"score\":4.737952},\n      {\n        \"id\":\"10.1371/journal.pone.0138021\",\n        \"journal\":\"PLOS ONE\",\n        \"eissn\":\"1932-6203\",\n        \"publication_date\":\"2015-09-16T00:00:00Z\",\n        \"article_type\":\"Research Article\",\n        \"author_display\":[\"Pavel Grigoriev\",\n          \"Evgeny M. Andreev\"],\n        \"abstract\":[\"Background and Aim: Harmful alcohol consumption has long been recognized as being the major determinant of male premature mortality in the European countries of the former USSR. Our focus here is on Belarus and Russia, two Slavic countries which continue to suffer enormously from the burden of the harmful consumption of alcohol. However, after a long period of deterioration, mortality trends in these countries have been improving over the past decade. We aim to investigate to what extent the recent declines in adult mortality in Belarus and Russia are attributable to the anti-alcohol measures introduced in these two countries in the 2000s. Data and Methods: We rely on the detailed cause-specific mortality series for the period 1980–2013. Our analysis focuses on the male population, and considers only a limited number of causes of death which we label as being alcohol-related: accidental poisoning by alcohol, liver cirrhosis, ischemic heart diseases, stroke, transportation accidents, and other external causes. For each of these causes we computed age-standardized death rates. The life table decomposition method was used to determine the age groups and the causes of death responsible for changes in life expectancy over time. Conclusion: Our results do not lead us to conclude that the schedule of anti-alcohol measures corresponds to the schedule of mortality changes. The continuous reduction in adult male mortality seen in Belarus and Russia cannot be fully explained by the anti-alcohol policies implemented in these countries, although these policies likely contributed to the large mortality reductions observed in Belarus and Russia in 2005–2006 and in Belarus in 2012. Thus, the effects of these policies appear to have been modest. We argue that the anti-alcohol measures implemented in Belarus and Russia simply coincided with fluctuations in alcohol-related mortality which originated in the past. If these trends had not been underway already, these huge mortality effects would not have occurred. \"],\n        \"title_display\":\"The Huge Reduction in Adult Male Mortality in Belarus and Russia: Is It Attributable to Anti-Alcohol Measures?\",\n        \"score\":4.7355423}]\n  },\n  \"highlighting\":{\n    \"10.1371/journal.pone.0218147\":{\n      \"abstract\":[\"Background: Binge drinking, an increasingly common form of <em>alcohol</em> use disorder, is associated\"]},\n    \"10.1371/journal.pone.0138021\":{\n      \"abstract\":[\"Background and Aim: Harmful <em>alcohol</em> consumption has long been recognized as being the major\"]}}}\n"
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

[Please report any issues or bugs](https://github.com/ropensci/solrium/issues).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_get.R
\name{config_get}
\alias{config_get}
\title{Get Solr configuration details}
\usage{
config_get(conn, name, what = NULL, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core. If not given, all cores.}

\item{what}{(character) What you want to look at. One of solrconfig or
schema. Default: solrconfig}

\item{wt}{(character) One of json (default) or xml. Data type returned.
If json, uses \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses
\code{\link[xml2:read_xml]{xml2::read_xml()}} to parse.}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list, \code{xml_document}, or character
}
\description{
Get Solr configuration details
}
\details{
Note that if \code{raw=TRUE}, \code{what} is ignored. That is,
you get all the data when \code{raw=TRUE}.
}
\examples{
\dontrun{
# start Solr with Cloud mode via the schemaless eg: bin/solr -e cloud
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# all config settings
conn$config_get("gettingstarted")

# just znodeVersion
conn$config_get("gettingstarted", "znodeVersion")

# just znodeVersion
conn$config_get("gettingstarted", "luceneMatchVersion")

# just updateHandler
conn$config_get("gettingstarted", "updateHandler")

# just updateHandler
conn$config_get("gettingstarted", "requestHandler")

## Get XML
conn$config_get("gettingstarted", wt = "xml")
conn$config_get("gettingstarted", "updateHandler", wt = "xml")
conn$config_get("gettingstarted", "requestHandler", wt = "xml")

## Raw data - what param ignored when raw=TRUE
conn$config_get("gettingstarted", raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_reload.R
\name{core_reload}
\alias{core_reload}
\title{Reload a core}
\usage{
core_reload(conn, name, raw = FALSE, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Reload a core
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#  bin/solr start -e schemaless
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# Status of particular cores
if (conn$core_exists("gettingstarted")) {
  conn$core_reload("gettingstarted")
  conn$core_status("gettingstarted")
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_mlt.r
\name{solr_mlt}
\alias{solr_mlt}
\title{"more like this" search}
\usage{
solr_mlt(
  conn,
  name = NULL,
  params = NULL,
  body = NULL,
  callopts = list(),
  raw = FALSE,
  parsetype = "df",
  concat = ",",
  optimizeMaxRows = TRUE,
  minOptimizedRows = 50000L,
  progress = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if not needed.}

\item{params}{(list) a named list of parameters, results in a GET request
as long as no body parameters given}

\item{body}{(list) a named list of parameters, if given a POST request
will be performed}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{raw}{(logical) If TRUE, returns raw data in format specified by wt param}

\item{parsetype}{(character) One of 'list' or 'df'}

\item{concat}{(character) Character to concatenate elements of longer than length 1.
Note that this only works reliably when data format is json (wt='json'). The parsing
is more complicated in XML format, but you can do that on your own.}

\item{optimizeMaxRows}{(logical) If \code{TRUE}, then rows parameter will be
adjusted to the number of returned results by the same constraints.
It will only be applied if rows parameter is higher
than \code{minOptimizedRows}. Default: \code{TRUE}}

\item{minOptimizedRows}{(numeric) used by \code{optimizedMaxRows} parameter,
the minimum optimized rows. Default: 50000}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}. only supports
\code{httr::progress} for now. See the README for an example.}

\item{...}{Further args to be combined into query}
}
\value{
XML, JSON, a list, or data.frame
}
\description{
Returns only more like this items
}
\section{More like this parameters}{

\itemize{
\item q Query terms, defaults to '\emph{:}', or everything.
\item fq Filter query, this does not affect the search, only what gets returned
\item mlt.count The number of similar documents to return for each result. Default is 5.
\item mlt.fl The fields to use for similarity. NOTE: if possible these should have a stored
TermVector DEFAULT_FIELD_NAMES = new String[] {"contents"}
\item mlt.mintf Minimum Term Frequency - the frequency below which terms will be ignored in
the source doc. DEFAULT_MIN_TERM_FREQ = 2
\item mlt.mindf Minimum Document Frequency - the frequency at which words will be ignored which
do not occur in at least this many docs. DEFAULT_MIN_DOC_FREQ = 5
\item mlt.minwl minimum word length below which words will be ignored.
DEFAULT_MIN_WORD_LENGTH = 0
\item mlt.maxwl maximum word length above which words will be ignored.
DEFAULT_MAX_WORD_LENGTH = 0
\item mlt.maxqt maximum number of query terms that will be included in any generated query.
DEFAULT_MAX_QUERY_TERMS = 25
\item mlt.maxntp maximum number of tokens to parse in each example doc field that is not stored
with TermVector support. DEFAULT_MAX_NUM_TOKENS_PARSED = 5000
\item mlt.boost (true/false) set if the query will be boosted by the interesting term relevance.
DEFAULT_BOOST = false
\item mlt.qf Query fields and their boosts using the same format as that used in
DisMaxQParserPlugin. These fields must also be specified in mlt.fl.
\item fl Fields to return. We force 'id' to be returned so that there is a unique identifier
with each record.
\item wt (character) Data type returned, defaults to 'json'. One of json or xml. If json,
uses \code{\link[jsonlite]{fromJSON}} to parse. If xml, uses \code{\link[XML]{xmlParse}} to
parse. csv is only supported in \code{\link{solr_search}} and \code{\link{solr_all}}.
\item start Record to start at, default to beginning.
\item rows Number of records to return. Defaults to 10.
\item key API key, if needed.
}
}

\examples{
\dontrun{
# connect
(conn <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))

# more like this search
conn$mlt(params = list(q='*:*', mlt.count=2, mlt.fl='abstract', fl='score',
  fq="doc_type:full"))
conn$mlt(params = list(q='*:*', rows=2, mlt.fl='title', mlt.mindf=1,
  mlt.mintf=1, fl='alm_twitterCount'))
conn$mlt(params = list(q='title:"ecology" AND body:"cell"', mlt.fl='title',
  mlt.mindf=1, mlt.mintf=1, fl='counter_total_all', rows=5))
conn$mlt(params = list(q='ecology', mlt.fl='abstract', fl='title', rows=5))
solr_mlt(conn, params = list(q='ecology', mlt.fl='abstract',
  fl=c('score','eissn'), rows=5))
solr_mlt(conn, params = list(q='ecology', mlt.fl='abstract',
  fl=c('score','eissn'), rows=5, wt = "xml"))

# get raw data, and parse later if needed
out <- solr_mlt(conn, params=list(q='ecology', mlt.fl='abstract', fl='title',
 rows=2), raw=TRUE)
solr_parse(out, "df")
}
}
\references{
See https://lucene.apache.org/solr/guide/8_2/morelikethis.html
for more information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_exists.R
\name{collection_exists}
\alias{collection_exists}
\title{Check if a collection exists}
\usage{
collection_exists(conn, name, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core. If not given, all cores.}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A single boolean, \code{TRUE} or \code{FALSE}
}
\description{
Check if a collection exists
}
\details{
Simply calls \code{\link[=collection_list]{collection_list()}} internally
}
\examples{
\dontrun{
# start Solr with Cloud mode via the schemaless eg: bin/solr -e cloud
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below
(conn <- SolrClient$new())

# exists
conn$collection_exists("gettingstarted")

# doesn't exist
conn$collection_exists("hhhhhh")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_create.R
\name{collection_create}
\alias{collection_create}
\title{Add a collection}
\usage{
collection_create(
  conn,
  name,
  numShards = 1,
  maxShardsPerNode = 1,
  createNodeSet = NULL,
  collection.configName = NULL,
  replicationFactor = 1,
  router.name = NULL,
  shards = NULL,
  createNodeSet.shuffle = TRUE,
  router.field = NULL,
  autoAddReplicas = FALSE,
  async = NULL,
  raw = FALSE,
  callopts = list(),
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{numShards}{(integer) The number of shards to be created as part of the
collection. This is a required parameter when using the 'compositeId' router.}

\item{maxShardsPerNode}{(integer) When creating collections, the shards and/or replicas
are spread across all available (i.e., live) nodes, and two replicas of the same shard
will never be on the same node. If a node is not live when the CREATE operation is called,
it will not get any parts of the new collection, which could lead to too many replicas
being created on a single live node. Defining maxShardsPerNode sets a limit on the number
of replicas CREATE will spread to each node. If the entire collection can not be fit into
the live nodes, no collection will be created at all. Default: 1}

\item{createNodeSet}{(logical) Allows defining the nodes to spread the new collection
across. If not provided, the CREATE operation will create shard-replica spread across all
live Solr nodes. The format is a comma-separated list of node_names, such as
localhost:8983_solr, localhost:8984_solr, localhost:8985_solr. Default: \code{NULL}}

\item{collection.configName}{(character) Defines the name of the configurations (which
must already be stored in ZooKeeper) to use for this collection. If not provided, Solr
will default to the collection name as the configuration name. Default: \code{compositeId}}

\item{replicationFactor}{(integer) The number of replicas to be created for each shard.
Default: 1}

\item{router.name}{(character) The router name that will be used. The router defines
how documents will be distributed among the shards. The value can be either \code{implicit},
which uses an internal default hash, or \code{compositeId}, which allows defining the specific
shard to assign documents to. When using the 'implicit' router, the shards parameter is
required. When using the 'compositeId' router, the numShards parameter is required.
For more information, see also the section Document Routing. Default: \code{compositeId}}

\item{shards}{(character) A comma separated list of shard names, e.g.,
shard-x,shard-y,shard-z . This is a required parameter when using the 'implicit' router.}

\item{createNodeSet.shuffle}{(logical)    Controls whether or not the shard-replicas created
for this collection will be assigned to the nodes specified by the createNodeSet in a
sequential manner, or if the list of nodes should be shuffled prior to creating individual
replicas.  A 'false' value makes the results of a collection creation predictible and
gives more exact control over the location of the individual shard-replicas, but 'true'
can be a better choice for ensuring replicas are distributed evenly across nodes. Ignored
if createNodeSet is not also specified. Default: \code{TRUE}}

\item{router.field}{(character) If this field is specified, the router will look at the
value of the field in an input document to compute the hash and identify a shard instead of
looking at the uniqueKey field. If the field specified is null in the document, the document
will be rejected. Please note that RealTime Get or retrieval by id would also require the
parameter \emph{route} (or shard.keys) to avoid a distributed search.}

\item{autoAddReplicas}{(logical)    When set to true, enables auto addition of replicas on
shared file systems. See the section autoAddReplicas Settings for more details on settings
and overrides. Default: \code{FALSE}}

\item{async}{(character) Request ID to track this action which will be processed
asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Add a collection
}
\examples{
\dontrun{
# connect
(conn <- SolrClient$new())

if (!conn$collection_exists("helloWorld")) {
  conn$collection_create(name = "helloWorld")
}
if (!conn$collection_exists("tablesChairs")) {
  conn$collection_create(name = "tablesChairs")
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_requeststatus.R
\name{core_requeststatus}
\alias{core_requeststatus}
\title{Request status of asynchronous CoreAdmin API call}
\usage{
core_requeststatus(conn, requestid, raw = FALSE, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{requestid}{The name of one of the cores to be removed. Required}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Request status of asynchronous CoreAdmin API call
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#   bin/solr start -e schemaless

# FIXME: not tested yet...
# (conn <- SolrClient$new())
# conn$core_requeststatus(requestid = 1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_addreplicaprop.R
\name{collection_addreplicaprop}
\alias{collection_addreplicaprop}
\title{Add a replica property}
\usage{
collection_addreplicaprop(
  conn,
  name,
  shard,
  replica,
  property,
  property.value,
  shardUnique = FALSE,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{shard}{(character) Required. The name of the shard the replica
belongs to}

\item{replica}{(character) Required. The replica, e.g. core_node1.}

\item{property}{(character) Required. The property to add. Note: this will
have the literal 'property.' prepended to distinguish it from
system-maintained properties. So these two forms are equivalent:
\code{property=special} and \code{property=property.special}}

\item{property.value}{(character) Required. The value to assign to
the property}

\item{shardUnique}{(logical) If \code{TRUE}, then setting this property in one
replica will (1) remove the property from all other replicas in that shard
Default: \code{FALSE}}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Assign an arbitrary property to a particular replica and give it
the value specified. If the property already exists, it will be overwritten
with the new value.
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("addrep")) {
  conn$collection_create(name = "addrep", numShards = 1)
  # OR bin/solr create -c addrep
}

# status
conn$collection_clusterstatus()$cluster$collections$addrep$shards

# add the value world to the property hello
conn$collection_addreplicaprop(name = "addrep", shard = "shard1",
  replica = "core_node1", property = "hello", property.value = "world")

# check status
conn$collection_clusterstatus()$cluster$collections$addrep$shards
conn$collection_clusterstatus()$cluster$collections$addrep$shards$shard1$replicas$core_node1
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/delete.R
\name{delete}
\alias{delete}
\alias{delete_by_id}
\alias{delete_by_query}
\title{Delete documents by ID or query}
\usage{
delete_by_id(
  conn,
  ids,
  name,
  commit = TRUE,
  commit_within = NULL,
  overwrite = TRUE,
  boost = NULL,
  wt = "json",
  raw = FALSE,
  ...
)

delete_by_query(
  conn,
  query,
  name,
  commit = TRUE,
  commit_within = NULL,
  overwrite = TRUE,
  boost = NULL,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{ids}{Document IDs, one or more in a vector or list}

\item{name}{(character) A collection or core name. Required.}

\item{commit}{(logical) If \code{TRUE}, documents immediately searchable.
Deafult: \code{TRUE}}

\item{commit_within}{(numeric) Milliseconds to commit the change, the
document will be added within that time. Default: \code{NULL}}

\item{overwrite}{(logical) Overwrite documents with matching keys.
Default: \code{TRUE}}

\item{boost}{(numeric) Boost factor. Default: \code{NULL}}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to
parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{query}{Query to use to delete documents}
}
\description{
Delete documents by ID or query
}
\details{
We use json internally as data interchange format for this function.
}
\examples{
\dontrun{
(cli <- SolrClient$new())

# add some documents first
ss <- list(list(id = 1, price = 100), list(id = 2, price = 500))
cli$add(ss, name = "gettingstarted")

# Now, delete them
# Delete by ID
cli$delete_by_id(ids = 1, "gettingstarted")
## Many IDs
cli$delete_by_id(ids = c(3, 4), "gettingstarted")

# Delete by query 
cli$search("gettingstarted", params=list(q="*:*")) # apple is there
cli$delete_by_query(query = 'id:apple', "gettingstarted") # delete it
cli$search("gettingstarted", params=list(q='id:apple')) # apple is now gone
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_requeststatus.R
\name{collection_requeststatus}
\alias{collection_requeststatus}
\title{Get request status}
\usage{
collection_requeststatus(conn, requestid, raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{requestid}{(character) Required. The user defined request-id for the
request. This can be used to track the status of the submitted asynchronous
task. \code{-1} is a special request id which is used to cleanup the stored
states for all of the already completed/failed tasks.}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Request the status of an already submitted Asynchronous
Collection API call. This call is also used to clear up the stored statuses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_createalias.R
\name{collection_createalias}
\alias{collection_createalias}
\title{Create an alias for a collection}
\usage{
collection_createalias(
  conn,
  alias,
  collections,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{alias}{(character) Required. The alias name to be created}

\item{collections}{(character) Required. A character vector of collections
to be aliased}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \code{\link[crul]{HttpClient}}}
}
\description{
Create a new alias pointing to one or more collections. If an
alias by the same name already exists, this action will replace the existing
alias, effectively acting like an atomic "MOVE" command.
}
\examples{
\dontrun{
(conn <- SolrClient$new())

if (!conn$collection_exists("thingsstuff")) {
  conn$collection_create(name = "thingsstuff")
}

conn$collection_createalias("tstuff", "thingsstuff")
conn$collection_clusterstatus()$cluster$collections$thingsstuff$aliases
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_xml.R
\name{update_xml}
\alias{update_xml}
\title{Update documents with XML data}
\usage{
update_xml(
  conn,
  files,
  name,
  commit = TRUE,
  optimize = FALSE,
  max_segments = 1,
  expunge_deletes = FALSE,
  wait_searcher = TRUE,
  soft_commit = FALSE,
  prepare_commit = NULL,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{files}{Path to a single file to load into Solr}

\item{name}{(character) Name of the core or collection}

\item{commit}{(logical) If \code{TRUE}, documents immediately searchable.
Deafult: \code{TRUE}}

\item{optimize}{Should index optimization be performed before the method returns.
Default: \code{FALSE}}

\item{max_segments}{optimizes down to at most this number of segments. Default: 1}

\item{expunge_deletes}{merge segments with deletes away. Default: \code{FALSE}}

\item{wait_searcher}{block until a new searcher is opened and registered as the
main query searcher, making the changes visible. Default: \code{TRUE}}

\item{soft_commit}{perform a soft commit - this will refresh the 'view' of the
index in a more performant manner, but without "on-disk" guarantees.
Default: \code{FALSE}}

\item{prepare_commit}{The prepareCommit command is an expert-level API that
calls Lucene's IndexWriter.prepareCommit(). Not passed by default}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite]{fromJSON}} to parse. If xml, uses
\code{\link[xml2]{read_xml}} to parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \code{\link[crul]{HttpClient}}}
}
\description{
Update documents with XML data
}
\details{
You likely may not be able to run this function against many
public Solr services, but should work locally.
}
\examples{
\dontrun{
# start Solr: bin/solr start -f -c -p 8983

# connect
(conn <- SolrClient$new())

# create a collection
if (!conn$collection_exists("books")) {
  conn$collection_create(name = "books", numShards = 2)
}

# Add documents
file <- system.file("examples", "books.xml", package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_xml(file, "books")

# Update commands - can include many varying commands
## Add files
file <- system.file("examples", "books2_delete.xml", package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_xml(file, "books")

## Delete files
file <- system.file("examples", "updatecommands_delete.xml",
package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_xml(file, "books")

## Add and delete in the same document
## Add a document first, that we can later delete
ss <- list(list(id = 456, name = "cat"))
conn$add(ss, "books")
## Now add a new document, and delete the one we just made
file <- system.file("examples", "add_delete.xml", package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_xml(file, "books")
}
}
\seealso{
Other update: 
\code{\link{update_csv}()},
\code{\link{update_json}()}
}
\concept{update}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_status.R
\name{core_status}
\alias{core_status}
\title{Get core status}
\usage{
core_status(
  conn,
  name = NULL,
  indexInfo = TRUE,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{indexInfo}{(logical)}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get core status
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#   bin/solr start -e schemaless
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# Status of all cores
conn$core_status()

# Status of particular cores
conn$core_status("gettingstarted")

# Get index info or not
## Default: TRUE
conn$core_status("gettingstarted", indexInfo = TRUE)
conn$core_status("gettingstarted", indexInfo = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_overseerstatus.R
\name{collection_overseerstatus}
\alias{collection_overseerstatus}
\title{Get overseer status}
\usage{
collection_overseerstatus(conn, raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Returns the current status of the overseer, performance
statistics of various overseer APIs as well as last 10 failures per
operation type.
}
\examples{
\dontrun{
(conn <- SolrClient$new())
conn$collection_overseerstatus()
res <- conn$collection_overseerstatus()
res$responseHeader
res$leader
res$overseer_queue_size
res$overseer_work_queue_size
res$overseer_operations
res$collection_operations
res$overseer_queue
res$overseer_internal_queue
res$collection_queue
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.r
\name{pivot_flatten_tabular}
\alias{pivot_flatten_tabular}
\title{Flatten facet.pivot responses}
\usage{
pivot_flatten_tabular(df_w_pivot)
}
\arguments{
\item{df_w_pivot}{a \code{data.frame} with another
\code{data.frame} nested inside representing a
pivot reponse}
}
\value{
a \code{data.frame}
}
\description{
Convert a nested hierarchy of facet.pivot elements
to tabular data (rows and columns)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_balanceshardunique.R
\name{collection_balanceshardunique}
\alias{collection_balanceshardunique}
\title{Balance a property}
\usage{
collection_balanceshardunique(
  conn,
  name,
  property,
  onlyactivenodes = TRUE,
  shardUnique = NULL,
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{property}{(character) Required. The property to balance. The literal
"property." is prepended to this property if not specified explicitly.}

\item{onlyactivenodes}{(logical) Normally, the property is instantiated
on active nodes only. If \code{FALSE}, then inactive nodes are also included
for distribution. Default: \code{TRUE}}

\item{shardUnique}{(logical) Something of a safety valve. There is one
pre-defined property (preferredLeader) that defaults this value to \code{TRUE}.
For all other properties that are balanced, this must be set to \code{TRUE} or
an error message is returned}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Insures that a particular property is distributed evenly
amongst the physical nodes that make up a collection. If the property
already exists on a replica, every effort is made to leave it there. If the
property is not on any replica on a shard one is chosen and the property
is added.
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("addrep")) {
  conn$collection_create(name = "mycollection")
  # OR: bin/solr create -c mycollection
}

# balance preferredLeader property
conn$collection_balanceshardunique("mycollection", property = "preferredLeader")

# examine cluster status
conn$collection_clusterstatus()$cluster$collections$mycollection
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_params.R
\name{config_params}
\alias{config_params}
\title{Set Solr configuration params}
\usage{
config_params(
  conn,
  name,
  param = NULL,
  set = NULL,
  unset = NULL,
  update = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core. If not given, all cores.}

\item{param}{(character) Name of a parameter}

\item{set}{(list) List of key:value pairs of what to set. Create or
overwrite a parameter set map. Default: NULL (nothing passed)}

\item{unset}{(list) One or more character strings of keys to unset.
Default: \code{NULL} (nothing passed)}

\item{update}{(list) List of key:value pairs of what to update. Updates
a parameter set map. This essentially overwrites the old parameter set,
so all parameters must be sent in each update request.}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list with response from server
}
\description{
Set Solr configuration params
}
\details{
The Request Parameters API allows creating parameter sets that can
override or take the place of parameters defined in solrconfig.xml. It is
really another endpoint of the Config API instead of a separate API, and
has distinct commands. It does not replace or modify any sections of
solrconfig.xml, but instead provides another approach to handling parameters
used in requests. It behaves in the same way as the Config API, by storing
parameters in another file that will be used at runtime. In this case,
the parameters are stored in a file named params.json. This file is kept in
ZooKeeper or in the conf directory of a standalone Solr instance.
}
\examples{
\dontrun{
# start Solr in standard or Cloud mode
# connect
(conn <- SolrClient$new())

# set a parameter set
myFacets <- list(myFacets = list(facet = TRUE, facet.limit = 5))
config_params(conn, "gettingstarted", set = myFacets)

# check a parameter
config_params(conn, "gettingstarted", param = "myFacets")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_group.r
\name{solr_group}
\alias{solr_group}
\title{Grouped search}
\usage{
solr_group(
  conn,
  name = NULL,
  params = NULL,
  body = NULL,
  callopts = list(),
  raw = FALSE,
  parsetype = "df",
  concat = ",",
  progress = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if not needed.}

\item{params}{(list) a named list of parameters, results in a GET request
as long as no body parameters given}

\item{body}{(list) a named list of parameters, if given a POST request
will be performed}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{raw}{(logical) If TRUE, returns raw data in format specified by wt param}

\item{parsetype}{(character) One of 'list' or 'df'}

\item{concat}{(character) Character to concatenate elements of longer than length 1.
Note that this only works reliably when data format is json (wt='json'). The parsing
is more complicated in XML format, but you can do that on your own.}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}. only supports
\code{httr::progress} for now. See the README for an example.}

\item{...}{Further args to be combined into query}
}
\value{
XML, JSON, a list, or data.frame
}
\description{
Returns only group items
}
\section{Group parameters}{

\itemize{
\item q Query terms, defaults to '\emph{:}', or everything.
\item fq Filter query, this does not affect the search, only what gets returned
\item fl Fields to return
\item wt (character) Data type returned, defaults to 'json'. One of json or xml. If json,
uses \code{\link[jsonlite]{fromJSON}} to parse. If xml, uses \code{\link[XML]{xmlParse}} to
parse. csv is only supported in \code{\link{solr_search}} and \code{\link{solr_all}}.
\item key API key, if needed.
\item group.field (fieldname) Group based on the unique values of a field. The
field must currently be single-valued and must be either indexed, or be another
field type that has a value source and works in a function query - such as
ExternalFileField. Note: for Solr 3.x versions the field must by a string like
field such as StrField or TextField, otherwise a http status 400 is returned.
\item group.func (function query) Group based on the unique values of a function
query. Solr4.0 This parameter only is supported on 4.0
\item group.query (query) Return a single group of documents that also match the
given query.
\item rows (number) The number of groups to return. Defaults to 10.
\item start (number) The offset into the list of groups.
\item group.limit (number) The number of results (documents) to return for each
group. Defaults to 1.
\item group.offset (number) The offset into the document list of each group.
\item sort How to sort the groups relative to each other. For example,
sort=popularity desc will cause the groups to be sorted according to the highest
popularity doc in each group. Defaults to "score desc".
\item group.sort How to sort documents within a single group. Defaults
to the same value as the sort parameter.
\item group.format One of grouped or simple. If simple, the grouped documents are
presented in a single flat list. The start and rows parameters refer to numbers of
documents instead of numbers of groups.
\item group.main (logical) If true, the result of the last field grouping command
is used as the main result list in the response, using group.format=simple
\item group.ngroups (logical) If true, includes the number of groups that have
matched the query. Default is false. Solr4.1 WARNING: If this parameter is set
to true on a sharded environment, all the documents that belong to the same group
have to be located in the same shard, otherwise the count will be incorrect. If you
are using SolrCloud, consider using "custom hashing"
\item group.cache.percent (0-100) If > 0 enables grouping cache. Grouping is executed
actual two searches. This option caches the second search. A value of 0 disables
grouping caching. Default is 0. Tests have shown that this cache only improves search
time with boolean queries, wildcard queries and fuzzy queries. For simple queries like
a term query or a match all query this cache has a negative impact on performance
}
}

\examples{
\dontrun{
# connect
(cli <- SolrClient$new())

# by default we do a GET request
cli$group("gettingstarted",
  params = list(q='*:*', group.field='compName_s'))
# OR
solr_group(cli, "gettingstarted",
  params = list(q='*:*', group.field='compName_s'))

# connect
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))

# Basic group query
solr_group(cli, params = list(q='ecology', group.field='journal',
  group.limit=3, fl=c('id','score')))
solr_group(cli, params = list(q='ecology', group.field='journal',
  group.limit=3, fl='article_type'))

# Different ways to sort (notice diff btw sort of group.sort)
# note that you can only sort on a field if you return that field
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
   fl=c('id','score')))
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
   fl=c('id','score','alm_twitterCount'), group.sort='alm_twitterCount desc'))
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
   fl=c('id','score','alm_twitterCount'), sort='score asc',
   group.sort='alm_twitterCount desc'))

# Two group.field values
out <- solr_group(cli, params = list(q='ecology', group.field=c('journal','article_type'),
  group.limit=3, fl='id'), raw=TRUE)
solr_parse(out)
solr_parse(out, 'df')

# Get two groups, one with alm_twitterCount of 0-10, and another group
# with 10 to infinity
solr_group(cli, params = list(q='ecology', group.limit=3, fl=c('id','alm_twitterCount'),
 group.query=c('alm_twitterCount:[0 TO 10]','alm_twitterCount:[10 TO *]')))

# Use of group.format and group.simple.
## The raw data structure of these two calls are slightly different, but
## the parsing inside the function outputs the same results. You can
## of course set raw=TRUE to get back what the data actually look like
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl=c('id','score'), group.format='simple'))
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl=c('id','score'), group.format='grouped'))
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl=c('id','score'), group.format='grouped', group.main='true'))

# xml back
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl=c('id','score'), wt = "xml"))
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl=c('id','score'), wt = "xml"), parsetype = "list")
res <- solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl=c('id','score'), wt = "xml"), raw = TRUE)
library("xml2")
xml2::read_xml(unclass(res))

solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl='article_type', wt = "xml"))
solr_group(cli, params = list(q='ecology', group.field='journal', group.limit=3,
  fl='article_type', wt = "xml"), parsetype = "list")
}
}
\references{
See
https://lucene.apache.org/solr/guide/8_2/collapse-and-expand-results.html
for more information.
}
\seealso{
\code{\link[=solr_highlight]{solr_highlight()}}, \code{\link[=solr_facet]{solr_facet()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ping.R
\name{ping}
\alias{ping}
\title{Ping a Solr instance}
\usage{
ping(conn, name, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) Name of a collection or core. Required.}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses [xml2::read_xml)] to parse

[xml2::read_xml)]: R:xml2::read_xml)}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
if \code{wt="xml"} an object of class \code{xml_document}, if
\code{wt="json"} an object of class \code{list}
}
\description{
Ping a Solr instance
}
\details{
You likely may not be able to run this function against many public
Solr services as they hopefully don't expose their admin interface to the
public, but works locally.
}
\examples{
\dontrun{
# start Solr, in your CLI, run: `bin/solr start -e cloud -noprompt`
# after that, if you haven't run `bin/post -c gettingstarted docs/` yet,
# do so

# connect: by default we connect to localhost, port 8983
(cli <- SolrClient$new())

# ping the gettingstarted index
cli$ping("gettingstarted")
ping(cli, "gettingstarted")
ping(cli, "gettingstarted", wt = "xml")
ping(cli, "gettingstarted", verbose = FALSE)
ping(cli, "gettingstarted", raw = TRUE)

ping(cli, "gettingstarted", wt="xml", verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_clusterprop.R
\name{collection_clusterprop}
\alias{collection_clusterprop}
\title{Add, edit, delete a cluster-wide property}
\usage{
collection_clusterprop(conn, name, val, raw = FALSE, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) Name of the core or collection}

\item{val}{(character) Required. The value of the property. If the value is
empty or null, the property is unset.}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Important: whether add, edit, or delete is used is determined
by the value passed to the \code{val} parameter. If the property name is
new, it will be added. If the property name exists, and the value is
different, it will be edited. If the property name exists, and the value
is \code{NULL} or empty the property is deleted (unset).
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# add the value https to the property urlScheme
collection_clusterprop(conn, name = "urlScheme", val = "https")

# status again
collection_clusterstatus(conn)$cluster$properties

# delete the property urlScheme by setting val to NULL or a 0 length string
collection_clusterprop(conn, name = "urlScheme", val = "")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_highlight.r
\name{solr_highlight}
\alias{solr_highlight}
\title{Highlighting search}
\usage{
solr_highlight(
  conn,
  name = NULL,
  params = NULL,
  body = NULL,
  callopts = list(),
  raw = FALSE,
  parsetype = "df",
  progress = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if not needed.}

\item{params}{(list) a named list of parameters, results in a GET request
as long as no body parameters given}

\item{body}{(list) a named list of parameters, if given a POST request
will be performed}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{raw}{(logical) If TRUE (default) raw json or xml returned. If FALSE,
parsed data returned.}

\item{parsetype}{One of list of df (data.frame)}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}. only supports
\code{httr::progress} for now. See the README for an example.}

\item{...}{Further args to be combined into query}
}
\value{
XML, JSON, a list, or data.frame
}
\description{
Returns only highlight items
}
\section{Facet parameters}{

\itemize{
\item q Query terms. See examples.
\item hl.fl A comma-separated list of fields for which to generate highlighted snippets.
If left blank, the fields highlighted for the LuceneQParser are the defaultSearchField
(or the df param if used) and for the DisMax parser the qf fields are used. A '\emph{' can
be used to match field globs, e.g. 'text_}' or even '\emph{' to highlight on all fields where
highlighting is possible. When using '}', consider adding hl.requireFieldMatch=TRUE.
\item hl.snippets Max no. of highlighted snippets to generate per field. Note:
it is possible for any number of snippets from zero to this value to be generated.
This parameter accepts per-field overrides. Default: 1.
\item hl.fragsize The size, in characters, of the snippets (aka fragments) created by
the highlighter. In the original Highlighter, "0" indicates that the whole field value
should be used with no fragmenting.
\item hl.q Set a query request to be highlighted. It overrides q parameter for
highlighting. Solr query syntax is acceptable for this parameter.
\item hl.mergeContiguous Collapse contiguous fragments into a single fragment. "true"
indicates contiguous fragments will be collapsed into single fragment. This parameter
accepts per-field overrides. This parameter makes sense for the original Highlighter
only. Default: FALSE.
\item hl.requireFieldMatch If TRUE, then a field will only be highlighted if the
query matched in this particular field (normally, terms are highlighted in all
requested fields regardless of which field matched the query). This only takes effect
if "hl.usePhraseHighlighter" is TRUE. Default: FALSE.
\item hl.maxAnalyzedChars How many characters into a document to look for suitable
snippets. This parameter makes sense for the original Highlighter only. Default: 51200.
You can assign a large value to this parameter and use hl.fragsize=0 to return
highlighting in large fields that have size greater than 51200 characters.
\item hl.alternateField If a snippet cannot be generated (due to no terms matching),
you can specify a field to use as the fallback. This parameter accepts per-field overrides.
\item hl.maxAlternateFieldLength If hl.alternateField is specified, this parameter
specifies the maximum number of characters of the field to return. Any value less than or
equal to 0 means unlimited. Default: unlimited.
\item hl.preserveMulti Preserve order of values in a multiValued list. Default: FALSE.
\item hl.maxMultiValuedToExamine When highlighting a multiValued field, stop examining
the individual entries after looking at this many of them. Will potentially return 0
snippets if this limit is reached before any snippets are found. If maxMultiValuedToMatch
is also specified, whichever limit is hit first will terminate looking for more.
Default: Integer.MAX_VALUE
\item hl.maxMultiValuedToMatch When highlighting a multiValued field, stop examining
the individual entries after looking at this many matches are found. If
maxMultiValuedToExamine is also specified, whichever limit is hit first will terminate
looking for more. Default: Integer.MAX_VALUE
\item hl.formatter Specify a formatter for the highlight output. Currently the only
legal value is "simple", which surrounds a highlighted term with a customizable pre- and
post text snippet. This parameter accepts per-field overrides. This parameter makes
sense for the original Highlighter only.
\item hl.simple.pre The text which appears before and after a highlighted term when using
the simple formatter. This parameter accepts per-field overrides. The default values are
\code{<em>} and \code{</em>} This parameter makes sense for the original Highlighter only. Use
hl.tag.pre and hl.tag.post for FastVectorHighlighter (see example under hl.fragmentsBuilder)
\item hl.simple.post The text which appears before and after a highlighted term when using
the simple formatter. This parameter accepts per-field overrides. The default values are
\code{<em>} and \code{</em>} This parameter makes sense for the original Highlighter only. Use
hl.tag.pre and hl.tag.post for FastVectorHighlighter (see example under hl.fragmentsBuilder)
\item hl.fragmenter Specify a text snippet generator for highlighted text. The standard
fragmenter is gap (which is so called because it creates fixed-sized fragments with gaps
for multi-valued fields). Another option is regex, which tries to create fragments that
"look like" a certain regular expression. This parameter accepts per-field overrides.
Default: "gap"
\item hl.fragListBuilder Specify the name of SolrFragListBuilder.  This parameter
makes sense for FastVectorHighlighter only. To create a fragSize=0 with the
FastVectorHighlighter, use the SingleFragListBuilder. This field supports per-field
overrides.
\item hl.fragmentsBuilder Specify the name of SolrFragmentsBuilder. This parameter makes
sense for FastVectorHighlighter only.
\item hl.boundaryScanner Configures how the boundaries of fragments are determined. By
default, boundaries will split at the character level, creating a fragment such as "uick
brown fox jumps over the la". Valid entries are breakIterator or simple, with breakIterator
being the most commonly used. This parameter makes sense for FastVectorHighlighter only.
\item hl.bs.maxScan Specify the length of characters to be scanned by SimpleBoundaryScanner.
Default: 10.  This parameter makes sense for FastVectorHighlighter only.
\item hl.bs.chars Specify the boundary characters, used by SimpleBoundaryScanner.
This parameter makes sense for FastVectorHighlighter only.
\item hl.bs.type Specify one of CHARACTER, WORD, SENTENCE and LINE, used by
BreakIteratorBoundaryScanner. Default: WORD. This parameter makes sense for
FastVectorHighlighter only.
\item hl.bs.language Specify the language for Locale that is used by
BreakIteratorBoundaryScanner. This parameter makes sense for FastVectorHighlighter only.
Valid entries take the form of ISO 639-1 strings.
\item hl.bs.country Specify the country for Locale that is used by
BreakIteratorBoundaryScanner. This parameter makes sense for FastVectorHighlighter only.
Valid entries take the form of ISO 3166-1 alpha-2 strings.
\item hl.useFastVectorHighlighter Use FastVectorHighlighter. FastVectorHighlighter
requires the field is termVectors=on, termPositions=on and termOffsets=on. This
parameter accepts per-field overrides. Default: FALSE
\item hl.usePhraseHighlighter Use SpanScorer to highlight phrase terms only when
they appear within the query phrase in the document. Default: TRUE.
\item hl.highlightMultiTerm If the SpanScorer is also being used, enables highlighting
for range/wildcard/fuzzy/prefix queries. Default: FALSE. This parameter makes sense
for the original Highlighter only.
\item hl.regex.slop Factor by which the regex fragmenter can stray from the ideal
fragment size (given by hl.fragsize) to accomodate the regular expression. For
instance, a slop of 0.2 with fragsize of 100 should yield fragments between 80
and 120 characters in length. It is usually good to provide a slightly smaller
fragsize when using the regex fragmenter. Default: .6. This parameter makes sense
for the original Highlighter only.
\item hl.regex.pattern The regular expression for fragmenting. This could be
used to extract sentences (see example solrconfig.xml) This parameter makes sense
for the original Highlighter only.
\item hl.regex.maxAnalyzedChars Only analyze this many characters from a field
when using the regex fragmenter (after which, the fragmenter produces fixed-sized
fragments). Applying a complicated regex to a huge field is expensive.
Default: 10000. This parameter makes sense for the original Highlighter only.
\item start Record to start at, default to beginning.
\item rows Number of records to return.
\item wt (character) Data type returned, defaults to 'json'. One of json or xml. If json,
uses \code{\link[jsonlite]{fromJSON}} to parse. If xml, uses \code{\link[XML]{xmlParse}} to
parse. csv is only supported in \code{\link{solr_search}} and \code{\link{solr_all}}.
\item fl Fields to return
\item fq Filter query, this does not affect the search, only what gets returned
}
}

\examples{
\dontrun{
# connect
(conn <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))

# highlight search
solr_highlight(conn, params = list(q='alcohol', hl.fl = 'abstract', rows=10),
  parsetype = "list")
solr_highlight(conn, params = list(q='alcohol', hl.fl = c('abstract','title'),
  rows=3), parsetype = "list")

# Raw data back
## json
solr_highlight(conn, params = list(q='alcohol', hl.fl = 'abstract', rows=10),
   raw=TRUE)
## xml
solr_highlight(conn, params = list(q='alcohol', hl.fl = 'abstract', rows=10,
   wt='xml'), raw=TRUE)
## parse after getting data back
out <- solr_highlight(conn, params = list(q='theoretical math',
   hl.fl = c('abstract','title'), hl.fragsize=30, rows=10, wt='xml'),
   raw=TRUE)
solr_parse(out, parsetype='list')
}
}
\references{
See https://lucene.apache.org/solr/guide/8_2/highlighting.html
for more information on highlighting.
}
\seealso{
\code{\link[=solr_search]{solr_search()}}, \code{\link[=solr_facet]{solr_facet()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.r
\name{collapse_pivot_names}
\alias{collapse_pivot_names}
\title{Collapse Pivot Field and Value Columns}
\usage{
collapse_pivot_names(data)
}
\arguments{
\item{data}{a \code{data.frame} with every 2 columns
representing a field and value and the final representing
a count}
}
\value{
a \code{data.frame}
}
\description{
Convert a table consisting of columns in sets of 3
into 2 columns assuming that the first column of every set of 3
(field) is duplicated throughout all rows and should be removed.
This type of structure is usually returned by facet.pivot responses.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_deletealias.R
\name{collection_deletealias}
\alias{collection_deletealias}
\title{Delete a collection alias}
\usage{
collection_deletealias(conn, alias, raw = FALSE, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{alias}{(character) Required. The alias name to be created}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Delete a collection alias
}
\examples{
\dontrun{
(conn <- SolrClient$new())

if (!conn$collection_exists("thingsstuff")) {
  conn$collection_create(name = "thingsstuff")
}

conn$collection_createalias("tstuff", "thingsstuff")
conn$collection_clusterstatus()$cluster$collections$thingsstuff$aliases # new alias
conn$collection_deletealias("tstuff")
conn$collection_clusterstatus()$cluster$collections$thingsstuff$aliases # gone
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.r
\name{solr_parse}
\alias{solr_parse}
\alias{solr_parse.sr_high}
\alias{solr_parse.sr_search}
\alias{solr_parse.sr_all}
\alias{solr_parse.sr_mlt}
\alias{solr_parse.sr_stats}
\alias{solr_parse.sr_group}
\title{Parse raw data from solr_search, solr_facet, or solr_highlight.}
\usage{
solr_parse(input, parsetype = NULL, concat)

\method{solr_parse}{sr_high}(input, parsetype = "list", concat = ",")

\method{solr_parse}{sr_search}(input, parsetype = "list", concat = ",")

\method{solr_parse}{sr_all}(input, parsetype = "list", concat = ",")

\method{solr_parse}{sr_mlt}(input, parsetype = "list", concat = ",")

\method{solr_parse}{sr_stats}(input, parsetype = "list", concat = ",")

\method{solr_parse}{sr_group}(input, parsetype = "list", concat = ",")
}
\arguments{
\item{input}{Output from solr_facet}

\item{parsetype}{One of 'list' or 'df' (data.frame)}

\item{concat}{Character to conactenate strings by, e.g,. ',' (character).
Used in solr_parse.sr_search only.}
}
\description{
Parse raw data from solr_search, solr_facet, or solr_highlight.
}
\details{
This is the parser used internally in solr_facet, but if you
output raw data from solr_facet using raw=TRUE, then you can use this
function to parse that data (a sr_facet S3 object) after the fact to a
list of data.frame's for easier consumption. The data format type is
detected from the attribute "wt" on the sr_facet object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_overlay.R
\name{config_overlay}
\alias{config_overlay}
\title{Get Solr configuration overlay}
\usage{
config_overlay(conn, name, omitHeader = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core. If not given, all cores.}

\item{omitHeader}{(logical) If \code{TRUE}, omit header. Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list with response from server
}
\description{
Get Solr configuration overlay
}
\examples{
\dontrun{
# start Solr with Cloud mode via the schemaless eg: bin/solr -e cloud
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# get config overlay
conn$config_overlay("gettingstarted")

# without header
conn$config_overlay("gettingstarted", omitHeader = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_json_request.R
\name{solr_json_request}
\alias{solr_json_request}
\title{Solr json request}
\usage{
solr_json_request(
  conn,
  name = NULL,
  body = NULL,
  callopts = list(),
  progress = NULL
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if not needed.}

\item{body}{(list) a named list, or a valid JSON character string}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}.
only supports \code{httr::progress} for now. See the README for an example.}
}
\value{
JSON character string
}
\description{
search using the JSON request API
}
\note{
SOLR v7.1 was first version to support this. See
\url{https://issues.apache.org/jira/browse/SOLR-11244}

POST request only, no GET request available
}
\examples{
\dontrun{
# Connect to a local Solr instance
(conn <- SolrClient$new())

## body as JSON 
a <- conn$json_request("gettingstarted", body = '{"query":"*:*"}')
jsonlite::fromJSON(a)
## body as named list
b <- conn$json_request("gettingstarted", body = list(query = "*:*"))
jsonlite::fromJSON(b)

## body as JSON 
a <- solr_json_request(conn, "gettingstarted", body = '{"query":"*:*"}')
jsonlite::fromJSON(a)
## body as named list
b <- solr_json_request(conn, "gettingstarted", body = list(query = "*:*"))
jsonlite::fromJSON(b)
}
}
\references{
See https://lucene.apache.org/solr/guide/7_6/json-request-api.html
for more information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.r
\name{is.sr_facet}
\alias{is.sr_facet}
\alias{is.sr_high}
\alias{is.sr_search}
\title{Test for sr_facet class}
\usage{
is.sr_facet(x)

is.sr_high(x)

is.sr_search(x)
}
\arguments{
\item{x}{Input}
}
\description{
Test for sr_facet class

Test for sr_high class

Test for sr_search class
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_facet.r
\name{solr_facet}
\alias{solr_facet}
\title{Faceted search}
\usage{
solr_facet(
  conn,
  name = NULL,
  params = list(q = "*:*"),
  body = NULL,
  callopts = list(),
  raw = FALSE,
  parsetype = "df",
  concat = ",",
  progress = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if not needed.}

\item{params}{(list) a named list of parameters, results in a GET request
as long as no body parameters given}

\item{body}{(list) a named list of parameters, if given a POST request
will be performed}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{raw}{(logical) If TRUE (default) raw json or xml returned. If FALSE,
parsed data returned.}

\item{parsetype}{(character) One of 'list' or 'df'}

\item{concat}{(character) Character to concatenate elements of longer than length 1.
Note that this only works reliably when data format is json (wt='json'). The parsing
is more complicated in XML format, but you can do that on your own.}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}. only supports
\code{httr::progress} for now. See the README for an example.}

\item{...}{Further args, usually per field arguments for faceting.}
}
\value{
Raw json or xml, or a list of length 4 parsed elements
(usually data.frame's).
}
\description{
Returns only facet items
}
\details{
A number of fields can be specified multiple times, in which case you can separate
them by commas, like \code{facet.field='journal,subject'}. Those fields are:
\itemize{
\item facet.field
\item facet.query
\item facet.date
\item facet.date.other
\item facet.date.include
\item facet.range
\item facet.range.other
\item facet.range.include
\item facet.pivot
}

\strong{Options for some parameters}:

\strong{facet.sort}: This param determines the ordering of the facet field constraints.
\itemize{
\item {count} sort the constraints by count (highest count first)
\item {index} to return the constraints sorted in their index order (lexicographic
by indexed term). For terms in the ascii range, this will be alphabetically sorted.
}
The default is count if facet.limit is greater than 0, index otherwise. This
parameter can be specified on a per field basis.

\strong{facet.method}:
This parameter indicates what type of algorithm/method to use when faceting a field.
\itemize{
\item {enum} Enumerates all terms in a field, calculating the set intersection of
documents that match the term with documents that match the query. This was the
default (and only) method for faceting multi-valued fields prior to Solr 1.4.
\item {fc} (Field Cache) The facet counts are calculated by iterating over documents
that match the query and summing the terms that appear in each document. This was
the default method for single valued fields prior to Solr 1.4.
\item {fcs} (Field Cache per Segment) works the same as fc except the underlying
cache data structure is built for each segment of the index individually
}
The default value is fc (except for BoolField which uses enum) since it tends to use
less memory and is faster then the enumeration method when a field has many unique
terms in the index. For indexes that are changing rapidly in NRT situations, fcs may
be a better choice because it reduces the overhead of building the cache structures
on the first request and/or warming queries when opening a new searcher -- but tends
to be somewhat slower then fc for subsequent requests against the same searcher. This
parameter can be specified on a per field basis.

\strong{facet.date.other}: This param indicates that in addition to the counts for each date
range constraint between facet.date.start and facet.date.end, counts should also be
computed for...
\itemize{
\item {before} All records with field values lower then lower bound of the first
range
\item {after} All records with field values greater then the upper bound of the
last range
\item {between} All records with field values between the start and end bounds
of all ranges
\item {none} Compute none of this information
\item {all} Shortcut for before, between, and after
}
This parameter can be specified on a per field basis. In addition to the all option,
this parameter can be specified multiple times to indicate multiple choices -- but
none will override all other options.

\strong{facet.date.include}: By default, the ranges used to compute date faceting between
facet.date.start and facet.date.end are all inclusive of both endpoints, while
the "before" and "after" ranges are not inclusive. This behavior can be modified
by the facet.date.include param, which can be any combination of the following
options...
\itemize{
\item{lower} All gap based ranges include their lower bound
\item{upper} All gap based ranges include their upper bound
\item{edge} The first and last gap ranges include their edge bounds (ie: lower
for the first one, upper for the last one) even if the corresponding upper/lower
option is not specified
\item{outer} The "before" and "after" ranges will be inclusive of their bounds,
even if the first or last ranges already include those boundaries.
\item{all} Shorthand for lower, upper, edge, outer
}
This parameter can be specified on a per field basis. This parameter can be specified
multiple times to indicate multiple choices.

\strong{facet.date.include}: This param indicates that in addition to the counts for each range
constraint between facet.range.start and facet.range.end, counts should also be
computed for...
\itemize{
\item{before} All records with field values lower then lower bound of the first
range
\item{after} All records with field values greater then the upper bound of the
last range
\item{between} All records with field values between the start and end bounds
of all ranges
\item{none} Compute none of this information
\item{all} Shortcut for before, between, and after
}
This parameter can be specified on a per field basis. In addition to the all option,
this parameter can be specified multiple times to indicate multiple choices -- but
none will override all other options.

\strong{facet.range.include}: By default, the ranges used to compute range faceting between
facet.range.start and facet.range.end are inclusive of their lower bounds and
exclusive of the upper bounds. The "before" range is exclusive and the "after"
range is inclusive. This default, equivalent to lower below, will not result in
double counting at the boundaries. This behavior can be modified by the
facet.range.include param, which can be any combination of the following options...
\itemize{
\item{lower} All gap based ranges include their lower bound
\item{upper} All gap based ranges include their upper bound
\item{edge} The first and last gap ranges include their edge bounds (ie: lower
for the first one, upper for the last one) even if the corresponding upper/lower
option is not specified
\item{outer} The "before" and "after" ranges will be inclusive of their bounds,
even if the first or last ranges already include those boundaries.
\item{all} Shorthand for lower, upper, edge, outer
}
Can be specified on a per field basis. Can be specified multiple times to indicate
multiple choices. If you want to ensure you don't double-count, don't choose both
lower & upper, don't choose outer, and don't choose all.
}
\section{Facet parameters}{

\itemize{
\item name Name of a collection or core. Or leave as \code{NULL} if not needed.
\item q Query terms. See examples.
\item facet.query This param allows you to specify an arbitrary query in the
Lucene default syntax to generate a facet count. By default, faceting returns
a count of the unique terms for a "field", while facet.query allows you to
determine counts for arbitrary terms or expressions. This parameter can be
specified multiple times to indicate that multiple queries should be used as
separate facet constraints. It can be particularly useful for numeric range
based facets, or prefix based facets -- see example below (i.e.
\code{price:[* TO 500]} and \code{price:[501 TO *]}).
\item facet.field This param allows you to specify a field which should be
treated as a facet. It will iterate over each Term in the field and generate a
facet count using that Term as the constraint. This parameter can be specified
multiple times to indicate multiple facet fields. None of the other params in
this section will have any effect without specifying at least one field name
using this param.
\item facet.prefix Limits the terms on which to facet to those starting with
the given string prefix. Note that unlike fq, this does not change the search
results -- it merely reduces the facet values returned to those beginning with
the specified prefix. This parameter can be specified on a per field basis.
\item facet.sort See Details.
\item facet.limit This param indicates the maximum number of constraint counts
that should be returned for the facet fields. A negative value means unlimited.
Default: 100. Can be specified on a per field basis.
\item facet.offset This param indicates an offset into the list of constraints
to allow paging. Default: 0. This parameter can be specified on a per field basis.
\item facet.mincount This param indicates the minimum counts for facet fields
should be included in the response. Default: 0. This parameter can be specified
on a per field basis.
\item facet.missing Set to "true" this param indicates that in addition to the
Term based constraints of a facet field, a count of all matching results which
have no value for the field should be computed. Default: FALSE. This parameter
can be specified on a per field basis.
\item facet.method See Details.
\item facet.enum.cache.minDf This param indicates the minimum document frequency
(number of documents matching a term) for which the filterCache should be used
when determining the constraint count for that term. This is only used when
facet.method=enum method of faceting. A value greater than zero will decrease
memory usage of the filterCache, but increase the query time. When faceting on
a field with a very large number of terms, and you wish to decrease memory usage,
try a low value of 25 to 50 first. Default: 0, causing the filterCache to be used
for all terms in the field. This parameter can be specified on a per field basis.
\item facet.threads This param will cause loading the underlying fields used in
faceting to be executed in parallel with the number of threads specified. Specify
as facet.threads=# where # is the maximum number of threads used. Omitting this
parameter or specifying the thread count as 0 will not spawn any threads just as
before. Specifying a negative number of threads will spin up to Integer.MAX_VALUE
threads. Currently this is limited to the fields, range and query facets are not
yet supported. In at least one case this has reduced warmup times from 20 seconds
to under 5 seconds.
\item facet.date Specify names of fields (of type DateField) which should be
treated as date facets. Can be specified multiple times to indicate multiple
date facet fields.
\item facet.date.start The lower bound for the first date range for all Date
Faceting on this field. This should be a single date expression which may use
the DateMathParser syntax. Can be specified on a per field basis.
\item facet.date.end The minimum upper bound for the last date range for all
Date Faceting on this field (see facet.date.hardend for an explanation of what
the actual end value may be greater). This should be a single date expression
which may use the DateMathParser syntax. Can be specified on a per field basis.
\item facet.date.gap The size of each date range expressed as an interval to
be added to the lower bound using the DateMathParser syntax. Eg:
facet.date.gap=+1DAY. Can be specified on a per field basis.
\item facet.date.hardend A Boolean parameter instructing Solr what to do in the
event that facet.date.gap does not divide evenly between facet.date.start and
facet.date.end. If this is true, the last date range constraint will have an
upper bound of facet.date.end; if false, the last date range will have the smallest
possible upper bound greater then facet.date.end such that the range is exactly
facet.date.gap wide. Default: FALSE. This parameter can be specified on a per
field basis.
\item facet.date.other See Details.
\item facet.date.include See Details.
\item facet.range Indicates what field to create range facets for. Example:
facet.range=price&facet.range=age
\item facet.range.start The lower bound of the ranges. Can be specified on a
per field basis. Example: f.price.facet.range.start=0.0&f.age.facet.range.start=10
\item facet.range.end The upper bound of the ranges. Can be specified on a per
field basis. Example: f.price.facet.range.end=1000.0&f.age.facet.range.start=99
\item facet.range.gap The size of each range expressed as a value to be added
to the lower bound. For date fields, this should be expressed using the
DateMathParser syntax. (ie: facet.range.gap=+1DAY). Can be specified
on a per field basis. Example: f.price.facet.range.gap=100&f.age.facet.range.gap=10
\item facet.range.hardend A Boolean parameter instructing Solr what to do in the
event that facet.range.gap does not divide evenly between facet.range.start and
facet.range.end. If this is true, the last range constraint will have an upper
bound of facet.range.end; if false, the last range will have the smallest possible
upper bound greater then facet.range.end such that the range is exactly
facet.range.gap wide. Default: FALSE. This parameter can be specified on a
per field basis.
\item facet.range.other See Details.
\item facet.range.include See Details.
\item facet.pivot This param allows you to specify a single comma-separated string
of fields to allow you to facet within the results of the parent facet to return
counts in the format of SQL group by operation
\item facet.pivot.mincount This param indicates the minimum counts for facet fields
to be included in the response. Default: 0. This parameter should only be specified
once.
\item start Record to start at, default to beginning.
\item rows Number of records to return.
\item key API key, if needed.
\item wt (character) Data type returned, defaults to 'json'. One of json or xml. If json,
uses \code{\link[jsonlite]{fromJSON}} to parse. If xml, uses \code{\link[XML]{xmlParse}} to
parse. csv is only supported in \code{\link{solr_search}} and \code{\link{solr_all}}.
}
}

\examples{
\dontrun{
# connect - local Solr instance
(cli <- SolrClient$new())
cli$facet("gettingstarted", params = list(q="*:*", facet.field='name'))
cli$facet("gettingstarted", params = list(q="*:*", facet.field='name'),
  callopts = list(verbose = TRUE))
cli$facet("gettingstarted", body = list(q="*:*", facet.field='name'),
  callopts = list(verbose = TRUE))

# Facet on a single field
solr_facet(cli, "gettingstarted", params = list(q='*:*', facet.field='name'))

# Remote instance
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))

# Facet on multiple fields
solr_facet(cli, params = list(q='alcohol',
  facet.field = c('journal','subject')))

# Using mincount
solr_facet(cli, params = list(q='alcohol', facet.field='journal',
  facet.mincount='500'))

# Using facet.query to get counts
solr_facet(cli, params = list(q='*:*', facet.field='journal',
  facet.query=c('cell','bird')))

# Using facet.pivot to simulate SQL group by counts
solr_facet(cli, params = list(q='alcohol', facet.pivot='journal,subject',
             facet.pivot.mincount=10))
## two or more fields are required - you can pass in as a single
## character string
solr_facet(cli, params = list(q='*:*', facet.pivot = "journal,subject",
  facet.limit =  3))
## Or, pass in as a vector of length 2 or greater
solr_facet(cli, params = list(q='*:*', facet.pivot = c("journal", "subject"),
  facet.limit =  3))

# Date faceting
solr_facet(cli, params = list(q='*:*', facet.date='publication_date',
  facet.date.start='NOW/DAY-5DAYS', facet.date.end='NOW',
  facet.date.gap='+1DAY'))
## two variables
solr_facet(cli, params = list(q='*:*',
  facet.date=c('publication_date', 'timestamp'),
  facet.date.start='NOW/DAY-5DAYS', facet.date.end='NOW',
  facet.date.gap='+1DAY'))

# Range faceting
solr_facet(cli, params = list(q='*:*', facet.range='counter_total_all',
  facet.range.start=5, facet.range.end=1000, facet.range.gap=10))

# Range faceting with > 1 field, same settings
solr_facet(cli, params = list(q='*:*',
  facet.range=c('counter_total_all','alm_twitterCount'),
  facet.range.start=5, facet.range.end=1000, facet.range.gap=10))

# Range faceting with > 1 field, different settings
solr_facet(cli, params = list(q='*:*',
  facet.range=c('counter_total_all','alm_twitterCount'),
  f.counter_total_all.facet.range.start=5,
  f.counter_total_all.facet.range.end=1000,
  f.counter_total_all.facet.range.gap=10,
  f.alm_twitterCount.facet.range.start=5,
  f.alm_twitterCount.facet.range.end=1000,
  f.alm_twitterCount.facet.range.gap=10))

# Get raw json or xml
## json
solr_facet(cli, params = list(q='*:*', facet.field='journal'), raw=TRUE)
## xml
solr_facet(cli, params = list(q='*:*', facet.field='journal', wt='xml'),
  raw=TRUE)

# Get raw data back, and parse later, same as what goes on internally if
# raw=FALSE (Default)
out <- solr_facet(cli, params = list(q='*:*', facet.field='journal'),
  raw=TRUE)
solr_parse(out)
out <- solr_facet(cli, params = list(q='*:*', facet.field='journal',
  wt = 'xml'), raw=TRUE)
solr_parse(out)

# Using the USGS BISON API (https://bison.usgs.gov/#solr)
## The occurrence endpoint
(cli <- SolrClient$new(host = "bison.usgs.gov", scheme = "https",
  path = "solr/occurrences/select", port = NULL))
solr_facet(cli, params = list(q='*:*', facet.field='year'))
solr_facet(cli, params = list(q='*:*', facet.field='computedStateFips'))

# using a proxy
# cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL,
#   proxy = list(url = "http://54.195.48.153:8888"))
# solr_facet(cli, params = list(facet.field='journal'),
#   callopts=list(verbose=TRUE))
}
}
\references{
See https://lucene.apache.org/solr/guide/8_2/faceting.html for
more information on faceting.
}
\seealso{
\code{\link[=solr_search]{solr_search()}}, \code{\link[=solr_highlight]{solr_highlight()}}, \code{\link[=solr_parse]{solr_parse()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_csv.R
\name{update_csv}
\alias{update_csv}
\title{Update documents with CSV data}
\usage{
update_csv(
  conn,
  files,
  name,
  separator = ",",
  header = TRUE,
  fieldnames = NULL,
  skip = NULL,
  skipLines = 0,
  trim = FALSE,
  encapsulator = NULL,
  escape = NULL,
  keepEmpty = FALSE,
  literal = NULL,
  map = NULL,
  split = NULL,
  rowid = NULL,
  rowidOffset = NULL,
  overwrite = NULL,
  commit = NULL,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{files}{Path to a single file to load into Solr}

\item{name}{(character) Name of the core or collection}

\item{separator}{Specifies the character to act as the field separator. Default: ','}

\item{header}{TRUE if the first line of the CSV input contains field or column names.
Default: \code{TRUE}. If the fieldnames parameter is absent, these field names
will be used when adding documents to the index.}

\item{fieldnames}{Specifies a comma separated list of field names to use when adding
documents to the Solr index. If the CSV input already has a header, the names
specified by this parameter will override them. Example: fieldnames=id,name,category}

\item{skip}{A comma separated list of field names to skip in the input. An alternate
way to skip a field is to specify it's name as a zero length string in fieldnames.
For example, \code{fieldnames=id,name,category&skip=name} skips the name field,
and is equivalent to \code{fieldnames=id,,category}}

\item{skipLines}{Specifies the number of lines in the input stream to discard
before the CSV data starts (including the header, if present). Default: \code{0}}

\item{trim}{If true remove leading and trailing whitespace from values. CSV parsing
already ignores leading whitespace by default, but there may be trailing whitespace,
or there may be leading whitespace that is encapsulated by quotes and is thus not
removed. This may be specified globally, or on a per-field basis.
Default: \code{FALSE}}

\item{encapsulator}{The character optionally used to surround values to preserve
characters such as the CSV separator or whitespace. This standard CSV format handles
the encapsulator itself appearing in an encapsulated value by doubling the
encapsulator.}

\item{escape}{The character used for escaping CSV separators or other reserved
characters. If an escape is specified, the encapsulator is not used unless also
explicitly specified since most formats use either encapsulation or escaping, not both.}

\item{keepEmpty}{Keep and index empty (zero length) field values. This may be specified
globally, or on a per-field basis. Default: \code{FALSE}}

\item{literal}{Adds fixed field name/value to all documents. Example: Adds a "datasource"
field with value equal to "products" for every document indexed from the CSV
\code{literal.datasource=products}}

\item{map}{Specifies a mapping between one value and another. The string on the LHS of
the colon will be replaced with the string on the RHS. This parameter can be specified
globally or on a per-field basis. Example: replaces "Absolutely" with "true" in every
field \code{map=Absolutely:true}. Example: removes any values of "RemoveMe" in the
field "foo" \code{f.foo.map=RemoveMe:&f.foo.keepEmpty=false }}

\item{split}{If TRUE, the field value is split into multiple values by another
CSV parser. The CSV parsing rules such as separator and encapsulator may be specified
as field parameters.}

\item{rowid}{If not null, add a new field to the document where the passed in parameter
name is the field name to be added and the current line/rowid is the value. This is
useful if your CSV doesn't have a unique id already in it and you want to use the line
number as one. Also useful if you simply want to index where exactly in the original
CSV file the row came from}

\item{rowidOffset}{In conjunction with the rowid parameter, this integer value will be
added to the rowid before adding it the field.}

\item{overwrite}{If true (the default), check for and overwrite duplicate documents,
based on the uniqueKey field declared in the solr schema. If you know the documents you
are indexing do not contain any duplicates then you may see a considerable speed up
with &overwrite=false.}

\item{commit}{Commit changes after all records in this request have been indexed. The
default is commit=false to avoid the potential performance impact of frequent commits.}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Update documents with CSV data
}
\note{
SOLR v1.2 was first version to support csv. See
https://issues.apache.org/jira/browse/SOLR-66
}
\examples{
\dontrun{
# start Solr: bin/solr start -f -c -p 8983

# connect
(conn <- SolrClient$new())

if (!conn$collection_exists("helloWorld")) {
  conn$collection_create(name = "helloWorld", numShards = 2)
}

df <- data.frame(id=1:3, name=c('red', 'blue', 'green'))
write.csv(df, file="df.csv", row.names=FALSE, quote = FALSE)
conn$update_csv("df.csv", "helloWorld", verbose = TRUE)

# give back raw xml
conn$update_csv("df.csv", "helloWorld", wt = "xml")
## raw json
conn$update_csv("df.csv", "helloWorld", wt = "json", raw = TRUE)
}
}
\seealso{
Other update: 
\code{\link{update_json}()},
\code{\link{update_xml}()}
}
\concept{update}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_split.R
\name{core_split}
\alias{core_split}
\title{Split a core}
\usage{
core_split(
  conn,
  name,
  path = NULL,
  targetCore = NULL,
  ranges = NULL,
  split.key = NULL,
  async = NULL,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{path}{(character) Two or more target directory paths in which a piece
of the index will be written}

\item{targetCore}{(character) Two or more target Solr cores to which a piece
of the index will be merged}

\item{ranges}{(character) A list of number ranges, or hash ranges in
hexadecimal format. If numbers, they get converted to hexidecimal format
before being passed to your Solr server.}

\item{split.key}{(character) The key to be used for splitting the index}

\item{async}{(character) Request ID to track this action which will be
processed asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
SPLIT splits an index into two or more indexes. The index being
split can continue to handle requests. The split pieces can be placed into
a specified directory on the server's filesystem or it can be merged into
running Solr cores.
}
\details{
The core index will be split into as many pieces as the number of
\code{path} or \code{targetCore} parameters.

Either \code{path} or \code{targetCore} parameter must be specified but not
both. The \code{ranges} and \code{split.key} parameters are optional and only one of
the two should be specified, if at all required.
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg: bin/solr start -e schemaless
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# Swap a core
## First, create two cores
# conn$core_split("splitcoretest0") # or create in the CLI: bin/solr create -c splitcoretest0
# conn$core_split("splitcoretest1") # or create in the CLI: bin/solr create -c splitcoretest1
# conn$core_split("splitcoretest2") # or create in the CLI: bin/solr create -c splitcoretest2

## check status
conn$core_status("splitcoretest0", FALSE)
conn$core_status("splitcoretest1", FALSE)
conn$core_status("splitcoretest2", FALSE)

## split core using targetCore parameter
conn$core_split("splitcoretest0", targetCore = c("splitcoretest1", "splitcoretest2"))

## split core using split.key parameter
### Here all documents having the same route key as the split.key i.e. 'A!'
### will be split from the core index and written to the targetCore
conn$core_split("splitcoretest0", targetCore = "splitcoretest1", split.key = "A!")

## split core using ranges parameter
### Solr expects hash ranges in hexidecimal, but since we're in R,
### let's not make our lives any harder, so you can pass in numbers
### but you can still pass in hexidecimal if you want.
rgs <- c('0-1f4', '1f5-3e8')
conn$core_split("splitcoretest0",
  targetCore = c("splitcoretest1", "splitcoretest2"), ranges = rgs)
rgs <- list(c(0, 500), c(501, 1000))
conn$core_split("splitcoretest0",
  targetCore = c("splitcoretest1", "splitcoretest2"), ranges = rgs)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SolrClient.R
\name{SolrClient}
\alias{SolrClient}
\title{Solr connection client}
\arguments{
\item{host}{(character) Host url. Deafault: 127.0.0.1}

\item{path}{(character) url path.}

\item{port}{(character/numeric) Port. Default: 8389}

\item{scheme}{(character) http scheme, one of http or https. Default: http}

\item{proxy}{List of arguments for a proxy connection, including one or
more of: url, port, username, password, and auth. See
\link[crul:proxies]{crul::proxy} for  help, which is used to construct the
proxy connection.}

\item{errors}{(character) One of \code{"simple"} or \code{"complete"}. Simple gives
http code and  error message on an error, while complete gives both http
code and error message, and stack trace, if available.}

\item{auth}{an object of class \code{auth}, created by calling \code{\link[crul:auth]{crul::auth()}}}
}
\value{
Various output, see help files for each grouping of methods.
}
\description{
Solr connection client
}
\details{
\code{SolrClient} creates a R6 class object. The object is
not cloneable and is portable, so it can be inherited across packages
without complication.

\code{SolrClient} is used to initialize a client that knows about your
Solr instance, with options for setting host, port, http scheme,
and simple vs. complete error reporting
}
\section{SolrClient methods}{


Each of these methods also has a matching standalone exported
function that you can use by passing in the connection object made
by calling \code{SolrClient$new()}. Also, see the docs for each method for
parameter definitions and their default values.
\itemize{
\item \code{ping(name, wt = 'json', raw = FALSE, ...)}
\item \code{schema(name, what = '', raw = FALSE, ...)}
\item \code{commit(name, expunge_deletes = FALSE, wait_searcher = TRUE, soft_commit = FALSE, wt = 'json', raw = FALSE, ...)}
\item \code{optimize(name, max_segments = 1, wait_searcher = TRUE, soft_commit = FALSE, wt = 'json', raw = FALSE, ...)}
\item \code{config_get(name, what = NULL, wt = "json", raw = FALSE, ...)}
\item \code{config_params(name, param = NULL, set = NULL, unset = NULL, update = NULL, ...)}
\item \code{config_overlay(name, omitHeader = FALSE, ...)}
\item \code{config_set(name, set = NULL, unset = NULL, ...)}
\item \code{collection_exists(name, ...)}
\item \code{collection_list(raw = FALSE, ...)}
\item \code{collection_create(name, numShards = 1, maxShardsPerNode = 1, createNodeSet = NULL, collection.configName = NULL, replicationFactor = 1, router.name = NULL, shards = NULL, createNodeSet.shuffle = TRUE, router.field = NULL, autoAddReplicas = FALSE, async = NULL, raw = FALSE, callopts=list(), ...)}
\item \code{collection_addreplica(name, shard = NULL, route = NULL, node = NULL, instanceDir = NULL, dataDir = NULL, async = NULL, raw = FALSE, callopts=list(), ...)}
\item \code{collection_addreplicaprop(name, shard, replica, property, property.value, shardUnique = FALSE, raw = FALSE, callopts=list())}
\item \code{collection_addrole(role = "overseer", node, raw = FALSE, ...)}
\item \code{collection_balanceshardunique(name, property, onlyactivenodes = TRUE, shardUnique = NULL, raw = FALSE, ...)}
\item \code{collection_clusterprop(name, val, raw = FALSE, callopts=list())}
\item \code{collection_clusterstatus(name = NULL, shard = NULL, raw = FALSE, ...)}
\item \code{collection_createalias(alias, collections, raw = FALSE, ...)}
\item \code{collection_createshard(name, shard, createNodeSet = NULL, raw = FALSE, ...)}
\item \code{collection_delete(name, raw = FALSE, ...)}
\item \code{collection_deletealias(alias, raw = FALSE, ...)}
\item \code{collection_deletereplica(name, shard = NULL, replica = NULL, onlyIfDown = FALSE, raw = FALSE, callopts=list(), ...)}
\item \code{collection_deletereplicaprop(name, shard, replica, property, raw = FALSE, callopts=list())}
\item \code{collection_deleteshard(name, shard, raw = FALSE, ...)}
\item \code{collection_migrate(name, target.collection, split.key, forward.timeout = NULL, async = NULL, raw = FALSE, ...)}
\item \code{collection_overseerstatus(raw = FALSE, ...)}
\item \code{collection_rebalanceleaders(name, maxAtOnce = NULL, maxWaitSeconds = NULL, raw = FALSE, ...)}
\item \code{collection_reload(name, raw = FALSE, ...)}
\item \code{collection_removerole(role = "overseer", node, raw = FALSE, ...)}
\item \code{collection_requeststatus(requestid, raw = FALSE, ...)}
\item \code{collection_splitshard(name, shard, ranges = NULL, split.key = NULL, async = NULL, raw = FALSE, ...)}
\item \code{core_status(name = NULL, indexInfo = TRUE, raw = FALSE, callopts=list())}
\item \code{core_exists(name, callopts = list())}
\item \code{core_create(name, instanceDir = NULL, config = NULL, schema = NULL, dataDir = NULL, configSet = NULL, collection = NULL, shard = NULL, async=NULL, raw = FALSE, callopts=list(), ...)}
\item \code{core_unload(name, deleteIndex = FALSE, deleteDataDir = FALSE, deleteInstanceDir = FALSE, async = NULL, raw = FALSE, callopts = list())}
\item \code{core_rename(name, other, async = NULL, raw = FALSE, callopts=list())}
\item \code{core_reload(name, raw = FALSE, callopts=list())}
\item \code{core_swap(name, other, async = NULL, raw = FALSE, callopts=list())}
\item \code{core_mergeindexes(name, indexDir = NULL, srcCore = NULL, async = NULL, raw = FALSE, callopts = list())}
\item \code{core_requeststatus(requestid, raw = FALSE, callopts = list())}
\item \code{core_split(name, path = NULL, targetCore = NULL, ranges = NULL, split.key = NULL, async = NULL, raw = FALSE, callopts=list())}
\item \code{search(name = NULL, params = NULL, body = NULL, callopts = list(), raw = FALSE,  parsetype = 'df', concat = ',', optimizeMaxRows = TRUE, minOptimizedRows = 50000L, progress = NULL, ...)}
\item \code{facet(name = NULL, params = NULL, body = NULL, callopts = list(), raw = FALSE,  parsetype = 'df', concat = ',', progress = NULL, ...)}
\item \code{stats(name = NULL, params = list(q = '*:*', stats.field = NULL, stats.facet = NULL), body = NULL, callopts=list(), raw = FALSE, parsetype = 'df', progress = NULL, ...)}
\item \code{highlight(name = NULL, params = NULL, body = NULL, callopts=list(), raw = FALSE, parsetype = 'df', progress = NULL, ...)}
\item \code{group(name = NULL, params = NULL, body = NULL, callopts=list(), raw=FALSE, parsetype='df', concat=',', progress = NULL, ...)}
\item \code{mlt(name = NULL, params = NULL, body = NULL, callopts=list(), raw=FALSE, parsetype='df', concat=',', optimizeMaxRows = TRUE, minOptimizedRows = 50000L, progress = NULL, ...)}
\item \code{all(name = NULL, params = NULL, body = NULL, callopts=list(), raw=FALSE, parsetype='df', concat=',', optimizeMaxRows = TRUE, minOptimizedRows = 50000L, progress = NULL, ...)}
\item \code{json_request(name = NULL, body = NULL, callopts=list(),  progress = NULL)}
\item \code{get(ids, name, fl = NULL, wt = 'json', raw = FALSE, ...)}
\item \code{add(x, name, commit = TRUE, commit_within = NULL, overwrite = TRUE, boost = NULL, wt = 'json', raw = FALSE, ...)}
\item \code{delete_by_id(ids, name, commit = TRUE, commit_within = NULL, overwrite = TRUE, boost = NULL, wt = 'json', raw = FALSE, ...)}
\item \code{delete_by_query(query, name, commit = TRUE, commit_within = NULL, overwrite = TRUE, boost = NULL, wt = 'json', raw = FALSE, ...)}
\item \code{update_json(files, name, commit = TRUE, optimize = FALSE, max_segments = 1, expunge_deletes = FALSE, wait_searcher = TRUE, soft_commit = FALSE, prepare_commit = NULL, wt = 'json', raw = FALSE, ...)}
\item \code{update_xml(files, name, commit = TRUE, optimize = FALSE, max_segments = 1, expunge_deletes = FALSE, wait_searcher = TRUE, soft_commit = FALSE, prepare_commit = NULL, wt = 'json', raw = FALSE, ...)}
\item \code{update_csv(files, name, separator = ',', header = TRUE, fieldnames = NULL, skip = NULL, skipLines = 0, trim = FALSE, encapsulator = NULL, escape = NULL, keepEmpty = FALSE, literal = NULL, map = NULL, split = NULL, rowid = NULL, rowidOffset = NULL, overwrite = NULL, commit = NULL, wt = 'json', raw = FALSE, ...)}
\item \code{update_atomic_json(body, name, wt = 'json', raw = FALSE, ...)}
\item \code{update_atomic_xml(body, name, wt = 'json', raw = FALSE, ...)}
}
}

\section{number of results}{

When the \verb{$search()} method returns a data.frame, metadata doesn't fit
into the output data.frame itself. You can access number of results
(\code{numFound}) in the attributes of the results. For example,
\code{attr(x, "numFound")} for number of results, and \code{attr(x, "start")}
for the offset value (if one was given). Or you can get all
attributes like \code{attributes(x)}. These metadata are not in the
attributes when requesting raw xml or json though as those metadata
are in the payload (unless \code{wt="csv"}).
}

\examples{
\dontrun{
# make a client
(cli <- SolrClient$new())

# variables
cli$host
cli$port
cli$path
cli$scheme

# ping
## ping to make sure it's up
cli$ping("gettingstarted")

# version
## get Solr version information
cli$schema("gettingstarted")
cli$schema("gettingstarted", "fields")
cli$schema("gettingstarted", "name")
cli$schema("gettingstarted", "version")$version

# Search
cli$search("gettingstarted", params = list(q = "*:*"))
cli$search("gettingstarted", body = list(query = "*:*"))

# set a different host
SolrClient$new(host = 'stuff.com')

# set a different port
SolrClient$new(host = 3456)

# set a different http scheme
SolrClient$new(scheme = 'https')

# set a proxy
SolrClient$new(proxy = list(url = "187.62.207.130:3128"))

prox <- list(url = "187.62.207.130:3128", user = "foo", pwd = "bar")
cli <- SolrClient$new(proxy = prox)
cli$proxy

# set simple authentication details
SolrClient$new(auth = crul::auth(user = "hello", pwd = "world"))

# A remote Solr instance to which you don't have admin access
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))
res <- cli$search(params = list(q = "memory"))
res
attr(res, "numFound")
attr(res, "start")
attr(res, "maxScore")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_json.R
\name{update_json}
\alias{update_json}
\title{Update documents with JSON data}
\usage{
update_json(
  conn,
  files,
  name,
  commit = TRUE,
  optimize = FALSE,
  max_segments = 1,
  expunge_deletes = FALSE,
  wait_searcher = TRUE,
  soft_commit = FALSE,
  prepare_commit = NULL,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{files}{Path to a single file to load into Solr}

\item{name}{(character) Name of the core or collection}

\item{commit}{(logical) If \code{TRUE}, documents immediately searchable.
Deafult: \code{TRUE}}

\item{optimize}{Should index optimization be performed before the method returns.
Default: \code{FALSE}}

\item{max_segments}{optimizes down to at most this number of segments. Default: 1}

\item{expunge_deletes}{merge segments with deletes away. Default: \code{FALSE}}

\item{wait_searcher}{block until a new searcher is opened and registered as the
main query searcher, making the changes visible. Default: \code{TRUE}}

\item{soft_commit}{perform a soft commit - this will refresh the 'view' of the
index in a more performant manner, but without "on-disk" guarantees.
Default: \code{FALSE}}

\item{prepare_commit}{The prepareCommit command is an expert-level API that
calls Lucene's IndexWriter.prepareCommit(). Not passed by default}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite]{fromJSON}} to parse. If xml, uses
\code{\link[xml2]{read_xml}} to parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \code{\link[crul]{HttpClient}}}
}
\description{
Update documents with JSON data
}
\details{
You likely may not be able to run this function against many
public Solr services, but should work locally.
}
\examples{
\dontrun{
# start Solr: bin/solr start -f -c -p 8983

# connect
(conn <- SolrClient$new())

# Add documents
file <- system.file("examples", "books2.json", package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_json(files = file, name = "books")
update_json(conn, files = file, name = "books")

# Update commands - can include many varying commands
## Add file
file <- system.file("examples", "updatecommands_add.json",
  package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_json(file, "books")

## Delete file
file <- system.file("examples", "updatecommands_delete.json",
  package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_json(file, "books")

# Add and delete in the same document
## Add a document first, that we can later delete
ss <- list(list(id = 456, name = "cat"))
conn$add(ss, "books")
}
}
\seealso{
Other update: 
\code{\link{update_csv}()},
\code{\link{update_xml}()}
}
\concept{update}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_stats.r
\name{solr_stats}
\alias{solr_stats}
\title{Solr stats}
\usage{
solr_stats(
  conn,
  name = NULL,
  params = list(q = "*:*", stats.field = NULL, stats.facet = NULL),
  body = NULL,
  callopts = list(),
  raw = FALSE,
  parsetype = "df",
  progress = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if
not needed.}

\item{params}{(list) a named list of parameters, results in a GET request
as long as no body parameters given}

\item{body}{(list) a named list of parameters, if given a POST request
will be performed}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{raw}{(logical) If TRUE, returns raw data in format specified by
wt param}

\item{parsetype}{(character) One of 'list' or 'df'}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}. only supports
\code{httr::progress} for now. See the README for an example.}

\item{...}{Further args to be combined into query}
}
\value{
XML, JSON, a list, or data.frame
}
\description{
Returns only stat items
}
\section{Stats parameters}{

\itemize{
\item q Query terms, defaults to '\emph{:}', or everything.
\item stats.field The number of similar documents to return for each result.
\item stats.facet You can not facet on multi-valued fields.
\item wt (character) Data type returned, defaults to 'json'. One of json
or xml. If json, uses \code{\link[jsonlite]{fromJSON}} to parse. If xml,
uses \code{\link[XML]{xmlParse}} to parse. csv is only supported in
\code{\link{solr_search}} and \code{\link{solr_all}}.
\item start Record to start at, default to beginning.
\item rows Number of records to return. Defaults to 10.
\item key API key, if needed.
}
}

\examples{
\dontrun{
# connect
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))

# get stats
solr_stats(cli, params = list(q='science', stats.field='counter_total_all'),
  raw=TRUE)
solr_stats(cli, params = list(q='title:"ecology" AND body:"cell"',
   stats.field=c('counter_total_all','alm_twitterCount')))
solr_stats(cli, params = list(q='ecology',
  stats.field=c('counter_total_all','alm_twitterCount'),
  stats.facet='journal'))
solr_stats(cli, params = list(q='ecology',
  stats.field=c('counter_total_all','alm_twitterCount'),
  stats.facet=c('journal','volume')))

# Get raw data, then parse later if you feel like it
## json
out <- solr_stats(cli, params = list(q='ecology',
  stats.field=c('counter_total_all','alm_twitterCount'),
  stats.facet=c('journal','volume')), raw=TRUE)
library("jsonlite")
jsonlite::fromJSON(out)
solr_parse(out) # list
solr_parse(out, 'df') # data.frame

## xml
out <- solr_stats(cli, params = list(q='ecology',
  stats.field=c('counter_total_all','alm_twitterCount'),
  stats.facet=c('journal','volume'), wt="xml"), raw=TRUE)
library("xml2")
xml2::read_xml(unclass(out))
solr_parse(out) # list
solr_parse(out, 'df') # data.frame

# Get verbose http call information
solr_stats(cli, params = list(q='ecology', stats.field='alm_twitterCount'),
   callopts=list(verbose=TRUE))
}
}
\references{
See
https://lucene.apache.org/solr/guide/8_2/the-stats-component.html for
more information on Solr stats.
}
\seealso{
\code{\link[=solr_highlight]{solr_highlight()}}, \code{\link[=solr_facet]{solr_facet()}}, \code{\link[=solr_search]{solr_search()}}, \code{\link[=solr_mlt]{solr_mlt()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_createshard.R
\name{collection_createshard}
\alias{collection_createshard}
\title{Create a shard}
\usage{
collection_createshard(
  conn,
  name,
  shard,
  createNodeSet = NULL,
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{shard}{(character) Required. The name of the shard to be created.}

\item{createNodeSet}{(character) Allows defining the nodes to spread the new
collection across. If not provided, the CREATE operation will create
shard-replica spread across all live Solr nodes. The format is a
comma-separated list of node_names, such as localhost:8983_solr,
localhost:8984_s olr, localhost:8985_solr.}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Create a shard
}
\examples{
\dontrun{
(conn <- SolrClient$new())
## FIXME - doesn't work right now
# conn$collection_create(name = "trees")
# conn$collection_createshard(name = "trees", shard = "newshard")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_splitshard.R
\name{collection_splitshard}
\alias{collection_splitshard}
\title{Create a shard}
\usage{
collection_splitshard(
  conn,
  name,
  shard,
  ranges = NULL,
  split.key = NULL,
  async = NULL,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{shard}{(character) Required. The name of the shard to be split}

\item{ranges}{(character) A comma-separated list of hash ranges in
hexadecimal e.g. ranges=0-1f4,1f5-3e8,3e9-5dc}

\item{split.key}{(character) The key to use for splitting the index}

\item{async}{(character) Request ID to track this action which will be
processed asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Create a shard
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("trees")) {
  conn$collection_create("trees")
}

# find shard names
names(conn$collection_clusterstatus()$cluster$collections$trees$shards)

# split a shard by name
conn$collection_splitshard(name = "trees", shard = "shard1")

# now we have three shards
names(conn$collection_clusterstatus()$cluster$collections$trees$shards)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_unload.R
\name{core_unload}
\alias{core_unload}
\title{Unload (delete) a core}
\usage{
core_unload(
  conn,
  name,
  deleteIndex = FALSE,
  deleteDataDir = FALSE,
  deleteInstanceDir = FALSE,
  async = NULL,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{deleteIndex}{(logical) If \code{TRUE}, will remove the index when
unloading the core. Default: \code{FALSE}}

\item{deleteDataDir}{(logical)    If \code{TRUE}, removes the data directory
and all sub-directories. Default: \code{FALSE}}

\item{deleteInstanceDir}{(logical)    If \code{TRUE}, removes everything related to
the core, including the index directory, configuration files and other
related files. Default: \code{FALSE}}

\item{async}{(character) Request ID to track this action which will be
processed asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Unload (delete) a core
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#   bin/solr start -e schemaless

# connect
(conn <- SolrClient$new())

# Create a core
conn$core_create(name = "books")

# Unload a core
conn$core_unload(name = "books")
## not found
# conn$core_unload(name = "books")
# > Error: 400 - Cannot unload non-existent core [books]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_set.R
\name{config_set}
\alias{config_set}
\title{Set Solr configuration details}
\usage{
config_set(conn, name, set = NULL, unset = NULL, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core. If not given, all cores.}

\item{set}{(list) List of key:value pairs of what to set.
Default: \code{NULL} (nothing passed)}

\item{unset}{(list) One or more character strings of keys to unset.
Default: \code{NULL} (nothing passed)}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list with response from server
}
\description{
Set Solr configuration details
}
\examples{
\dontrun{
# start Solr with Cloud mode via the schemaless eg: bin/solr -e cloud
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# set a property
conn$config_set("gettingstarted", 
  set = list(query.filterCache.autowarmCount = 1000))

# unset a property
conn$config_set("gettingstarted", unset = "query.filterCache.size", 
  verbose = TRUE)

# many properties
conn$config_set("gettingstarted", set = list(
   query.filterCache.autowarmCount = 1000,
   query.commitWithin.softCommit = 'false'
 )
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_mergeindexes.R
\name{core_mergeindexes}
\alias{core_mergeindexes}
\title{Merge indexes (cores)}
\usage{
core_mergeindexes(
  conn,
  name,
  indexDir = NULL,
  srcCore = NULL,
  async = NULL,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{indexDir}{(character)    Multi-valued, directories that would be merged.}

\item{srcCore}{(character)    Multi-valued, source cores that would be merged.}

\item{async}{(character) Request ID to track this action which will be
processed asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Merges one or more indexes to another index. The indexes must
have completed commits, and should be locked against writes until the merge
is complete or the resulting merged index may become corrupted. The target
core index must already exist and have a compatible schema with the one or
more indexes that will be merged to it.
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#  bin/solr start -e schemaless

# connect
(conn <- SolrClient$new())

## FIXME: not tested yet

# use indexDir parameter
# conn$core_mergeindexes(core="new_core_name",
#    indexDir = c("/solr_home/core1/data/index",
#    "/solr_home/core2/data/index"))

# use srcCore parameter
# conn$core_mergeindexes(name = "new_core_name", srcCore = c('core1', 'core2'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_removerole.R
\name{collection_removerole}
\alias{collection_removerole}
\title{Remove a role from a node}
\usage{
collection_removerole(conn, role = "overseer", node, raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{role}{(character) Required. The name of the role. The only supported
role as of now is overseer (set as default).}

\item{node}{(character) Required. The name of the node.}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Remove an assigned role. This API is used to undo the roles
assigned using \code{\link{collection_addrole}}
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# get list of nodes
nodes <- conn$collection_clusterstatus()$cluster$live_nodes
conn$collection_addrole(node = nodes[1])
conn$collection_removerole(node = nodes[1])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_search.r
\name{solr_search}
\alias{solr_search}
\title{Solr search}
\usage{
solr_search(
  conn,
  name = NULL,
  params = list(q = "*:*"),
  body = NULL,
  callopts = list(),
  raw = FALSE,
  parsetype = "df",
  concat = ",",
  optimizeMaxRows = TRUE,
  minOptimizedRows = 50000L,
  progress = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if not needed.}

\item{params}{(list) a named list of parameters, results in a GET request
as long as no body parameters given}

\item{body}{(list) a named list of parameters, if given a POST request
will be performed}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{raw}{(logical) If TRUE, returns raw data in format specified by wt param}

\item{parsetype}{(character) One of 'list' or 'df'}

\item{concat}{(character) Character to concatenate elements of longer than length 1.
Note that this only works reliably when data format is json (wt='json'). The parsing
is more complicated in XML format, but you can do that on your own.}

\item{optimizeMaxRows}{(logical) If \code{TRUE}, then rows parameter will be
adjusted to the number of returned results by the same constraints.
It will only be applied if rows parameter is higher
than \code{minOptimizedRows}. Default: \code{TRUE}}

\item{minOptimizedRows}{(numeric) used by \code{optimizedMaxRows} parameter,
the minimum optimized rows. Default: 50000}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}. only supports
\code{httr::progress} for now. See the README for an example.}

\item{...}{Further args to be combined into query}
}
\value{
XML, JSON, a list, or data.frame
}
\description{
Returns only matched documents, and doesn't return other items,
including facets, groups, mlt, stats, and highlights.
}
\note{
SOLR v1.2 was first version to support csv. See
https://issues.apache.org/jira/browse/SOLR-66
}
\section{Parameters}{

\itemize{
\item q Query terms, defaults to '\emph{:}', or everything.
\item sort Field to sort on. You can specify ascending (e.g., score desc) or
descending (e.g., score asc), sort by two fields (e.g., score desc, price asc),
or sort by a function (e.g., sum(x_f, y_f) desc, which sorts by the sum of
x_f and y_f in a descending order).
\item start Record to start at, default to beginning.
\item rows Number of records to return. Default: 10.
\item pageDoc If you expect to be paging deeply into the results (say beyond page 10,
assuming rows=10) and you are sorting by score, you may wish to add the pageDoc
and pageScore parameters to your request. These two parameters tell Solr (and Lucene)
what the last result (Lucene internal docid and score) of the previous page was,
so that when scoring the query for the next set of pages, it can ignore any results
that occur higher than that item. To get the Lucene internal doc id, you will need
to add \code{docid} to the &fl list.
\item pageScore See pageDoc notes.
\item fq Filter query, this does not affect the search, only what gets returned.
This parameter can accept multiple items in a lis or vector. You can't pass more than
one parameter of the same name, so we get around it by passing multiple queries
and we parse internally
\item fl Fields to return, can be a character vector like \code{c('id', 'title')},
or a single character vector with one or more comma separated names, like
\code{'id,title'}
\item defType Specify the query parser to use with this request.
\item timeAllowed The time allowed for a search to finish. This value only applies
to the search and not to requests in general. Time is in milliseconds. Values \code{<= 0}
mean no time restriction. Partial results may be returned (if there are any).
\item qt Which query handler used. Options: dismax, others?
\item NOW Set a fixed time for evaluating Date based expresions
\item TZ Time zone, you can override the default.
\item echoHandler If \code{TRUE}, Solr places the name of the handle used in the
response to the client for debugging purposes. Default:
\item echoParams The echoParams parameter tells Solr what kinds of Request
parameters should be included in the response for debugging purposes, legal values
include:
\itemize{
\item none - don't include any request parameters for debugging
\item explicit - include the parameters explicitly specified by the client in the request
\item all - include all parameters involved in this request, either specified explicitly
by the client, or implicit because of the request handler configuration.
}
\item wt (character) One of json, xml, or csv. Data type returned, defaults
to 'csv'. If json, uses \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml,
uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to parse. If csv, uses \code{\link[=read.table]{read.table()}} to parse.
\code{wt=csv} gives the fastest performance at least in all the cases we have
tested in, thus it's the default value for \code{wt}
}
}

\section{number of results}{

Because \code{solr_search()} returns a data.frame, metadata doesn't fit into the
output data.frame itself. You can access number of results (\code{numFound}) in
the attributes of the results. For example, \code{attr(x, "numFound")} for
number of results, and \code{attr(x, "start")} for the offset value (if one
was given). Or you can get all attributes like \code{attributes(x)}. These
metadata are not in the attributes when \code{raw=TRUE} as those metadata
are in the payload (unless \code{wt="csv"}).
}

\examples{
\dontrun{
# Connect to a local Solr instance
(cli <- SolrClient$new())
cli$search("gettingstarted", params = list(q = "features:notes"))

solr_search(cli, "gettingstarted")
solr_search(cli, "gettingstarted", params = list(q = "features:notes"))
solr_search(cli, "gettingstarted", body = list(query = "features:notes"))

(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))
cli$search(params = list(q = "*:*"))
cli$search(params = list(q = "title:golgi", fl = c('id', 'title')))

cli$search(params = list(q = "*:*", facet = "true"))


# search
solr_search(cli, params = list(q='*:*', rows=2, fl='id'))

# search and return all rows
solr_search(cli, params = list(q='*:*', rows=-1, fl='id'))

# Search for word ecology in title and cell in the body
solr_search(cli, params = list(q='title:"ecology" AND body:"cell"',
  fl='title', rows=5))

# Search for word "cell" and not "body" in the title field
solr_search(cli, params = list(q='title:"cell" -title:"lines"', fl='title',
  rows=5))

# Wildcards
## Search for word that starts with "cell" in the title field
solr_search(cli, params = list(q='title:"cell*"', fl='title', rows=5))

# Proximity searching
## Search for words "sports" and "alcohol" within four words of each other
solr_search(cli, params = list(q='everything:"sports alcohol"~7',
  fl='abstract', rows=3))

# Range searches
## Search for articles with Twitter count between 5 and 10
solr_search(cli, params = list(q='*:*', fl=c('alm_twitterCount','id'),
  fq='alm_twitterCount:[5 TO 50]', rows=10))

# Boosts
## Assign higher boost to title matches than to body matches
## (compare the two calls)
solr_search(cli, params = list(q='title:"cell" abstract:"science"',
  fl='title', rows=3))
solr_search(cli, params = list(q='title:"cell"^1.5 AND abstract:"science"',
  fl='title', rows=3))

# FunctionQuery queries
## This kind of query allows you to use the actual values of fields to
## calculate relevancy scores for returned documents

## Here, we search on the product of counter_total_all and alm_twitterCount
## metrics for articles in PLOS Journals
solr_search(cli, params = list(q="{!func}product($v1,$v2)",
  v1 = 'sqrt(counter_total_all)',
  v2 = 'log(alm_twitterCount)', rows=5, fl=c('id','title'),
  fq='doc_type:full'))

## here, search on the product of counter_total_all and alm_twitterCount,
## using a new temporary field "_val_"
solr_search(cli,
  params = list(q='_val_:"product(counter_total_all,alm_twitterCount)"',
  rows=5, fl=c('id','title'), fq='doc_type:full'))

## papers with most citations
solr_search(cli, params = list(q='_val_:"max(counter_total_all)"',
   rows=5, fl=c('id','counter_total_all'), fq='doc_type:full'))

## papers with most tweets
solr_search(cli, params = list(q='_val_:"max(alm_twitterCount)"',
   rows=5, fl=c('id','alm_twitterCount'), fq='doc_type:full'))

## many fq values
solr_search(cli, params = list(q="*:*", fl=c('id','alm_twitterCount'),
   fq=list('doc_type:full','subject:"Social networks"',
           'alm_twitterCount:[100 TO 10000]'),
   sort='counter_total_month desc'))

## using wt = csv
solr_search(cli, params = list(q='*:*', rows=50, fl=c('id','score'),
  fq='doc_type:full', wt="csv"))
solr_search(cli, params = list(q='*:*', rows=50, fl=c('id','score'),
  fq='doc_type:full'))

# using a proxy
# cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL,
#   proxy = list(url = "http://186.249.1.146:80"))
# solr_search(cli, q='*:*', rows=2, fl='id', callopts=list(verbose=TRUE))

# Pass on curl options to modify request
## verbose
solr_search(cli, params = list(q='*:*', rows=2, fl='id'),
  callopts = list(verbose=TRUE))

# using a cursor for deep paging
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))
## json, raw data
res <- solr_search(cli, params = list(q = '*:*', rows = 100, sort = "id asc", cursorMark = "*"), 
  parsetype = "json", raw = TRUE, callopts=list(verbose=TRUE))
res
## data.frame
res <- solr_search(cli, params = list(q = '*:*', rows = 100, sort = "id asc", cursorMark = "*"))
res
attributes(res)
attr(res, "nextCursorMark")
## list
res <- solr_search(cli, params = list(q = '*:*', rows = 100, sort = "id asc", cursorMark = "*"),
  parsetype = "list")
res
attributes(res)
attr(res, "nextCursorMark")
}
}
\references{
See https://lucene.apache.org/solr/guide/8_2/searching.html
for more information.
}
\seealso{
\code{\link[=solr_highlight]{solr_highlight()}}, \code{\link[=solr_facet]{solr_facet()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_clusterstatus.R
\name{collection_clusterstatus}
\alias{collection_clusterstatus}
\title{Get cluster status}
\usage{
collection_clusterstatus(conn, name = NULL, shard = NULL, raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{shard}{(character) The shard(s) for which information is requested.
Multiple shard names can be specified as a character vector.}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Fetch the cluster status including collections, shards,
replicas, configuration name as well as collection aliases and cluster
properties.
}
\examples{
\dontrun{
(conn <- SolrClient$new())
conn$collection_clusterstatus()
res <- conn$collection_clusterstatus()
res$responseHeader
res$cluster
res$cluster$collections
res$cluster$collections$gettingstarted
res$cluster$live_nodes
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_deletereplicaprop.R
\name{collection_deletereplicaprop}
\alias{collection_deletereplicaprop}
\title{Delete a replica property}
\usage{
collection_deletereplicaprop(
  conn,
  name,
  shard,
  replica,
  property,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{shard}{(character) Required. The name of the shard the replica
belongs to.}

\item{replica}{(character) Required. The replica, e.g. core_node1.}

\item{property}{(character) Required. The property to delete. Note: this
will have the literal 'property.' prepended to distinguish it from
system-maintained properties. So these two forms are equivalent:
\code{property=special} and  \code{property=property.special}}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Deletes an arbitrary property from a particular replica.
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("deleterep")) {
  conn$collection_create(name = "deleterep")
  # OR bin/solr create -c deleterep
}

# status
conn$collection_clusterstatus()$cluster$collections$deleterep$shards

# add the value bar to the property foo
conn$collection_addreplicaprop(name = "deleterep", shard = "shard1",
  replica = "core_node1", property = "foo", property.value = "bar")

# check status
conn$collection_clusterstatus()$cluster$collections$deleterep$shards
conn$collection_clusterstatus()$cluster$collections$deleterep$shards$shard1$replicas$core_node1

# delete replica property
conn$collection_deletereplicaprop(name = "deleterep", shard = "shard1",
   replica = "core_node1", property = "foo")

# check status - foo should be gone
conn$collection_clusterstatus()$cluster$collections$deleterep$shards$shard1$replicas$core_node1
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_addreplica.R
\name{collection_addreplica}
\alias{collection_addreplica}
\title{Add a replica}
\usage{
collection_addreplica(
  conn,
  name,
  shard = NULL,
  route = NULL,
  node = NULL,
  instanceDir = NULL,
  dataDir = NULL,
  async = NULL,
  raw = FALSE,
  callopts = list(),
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{shard}{(character) The name of the shard to which replica is to be added.
If \code{shard} is not given, then \code{route} must be.}

\item{route}{(character) If the exact shard name is not known, users may pass
the \code{route} value and the system would identify the name of the shard.
Ignored if the \code{shard} param is also given}

\item{node}{(character) The name of the node where the replica should be created}

\item{instanceDir}{(character) The instanceDir for the core that will be created}

\item{dataDir}{(character)    The directory in which the core should be created}

\item{async}{(character) Request ID to track this action which will be processed
asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for details on
supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Add a replica to a shard in a collection. The node name can be
specified if the replica is to be created in a specific node
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("foobar")) {
  conn$collection_create(name = "foobar", numShards = 2)
  # OR bin/solr create -c foobar
}

# status
conn$collection_clusterstatus()$cluster$collections$foobar

# add replica
if (!conn$collection_exists("foobar")) {
  conn$collection_addreplica(name = "foobar", shard = "shard1")
}

# status again
conn$collection_clusterstatus()$cluster$collections$foobar
conn$collection_clusterstatus()$cluster$collections$foobar$shards
conn$collection_clusterstatus()$cluster$collections$foobar$shards$shard1
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize.R
\name{solr_optimize}
\alias{solr_optimize}
\title{Optimize}
\usage{
solr_optimize(
  conn,
  name,
  max_segments = 1,
  wait_searcher = TRUE,
  soft_commit = FALSE,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) A collection or core name. Required.}

\item{max_segments}{optimizes down to at most this number of segments.
Default: 1}

\item{wait_searcher}{block until a new searcher is opened and registered
as the main query searcher, making the changes visible. Default: \code{TRUE}}

\item{soft_commit}{perform a soft commit - this will refresh the 'view'
of the index in a more performant manner, but without "on-disk" guarantees.
Default: \code{FALSE}}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to
parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Optimize
}
\examples{
\dontrun{
(conn <- SolrClient$new())

solr_optimize(conn, "gettingstarted")
solr_optimize(conn, "gettingstarted", max_segments = 2)
solr_optimize(conn, "gettingstarted", wait_searcher = FALSE)

# get xml back
solr_optimize(conn, "gettingstarted", wt = "xml")
## raw xml
solr_optimize(conn, "gettingstarted", wt = "xml", raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solrium-package.R
\docType{package}
\name{solrium-package}
\alias{solrium-package}
\alias{solrium}
\title{General purpose R interface to Solr.}
\description{
This package has support for all the search endpoints, as well as a suite
of functions for managing a Solr database, including adding and deleting
documents.
}
\section{Important search functions}{

\itemize{
\item \code{\link[=solr_search]{solr_search()}} - General search, only returns documents
\item \code{\link[=solr_all]{solr_all()}} - General search, including all non-documents
in addition to documents: facets, highlights, groups, mlt, stats.
\item \code{\link[=solr_facet]{solr_facet()}} - Faceting only (w/o general search)
\item \code{\link[=solr_highlight]{solr_highlight()}} - Highlighting only (w/o general search)
\item \code{\link[=solr_mlt]{solr_mlt()}} - More like this (w/o general search)
\item \code{\link[=solr_group]{solr_group()}} - Group search (w/o general search)
\item \code{\link[=solr_stats]{solr_stats()}} - Stats search (w/o general search)
}
}

\section{Important Solr management functions}{

\itemize{
\item \code{\link[=update_json]{update_json()}} - Add or delete documents using json in a file
\item \code{\link[=add]{add()}} - Add documents via an R list or data.frame
\item \code{\link[=delete_by_id]{delete_by_id()}} - Delete documents by ID
\item \code{\link[=delete_by_query]{delete_by_query()}} - Delete documents by query
}
}

\section{Vignettes}{


See the vignettes for help \code{browseVignettes(package = "solrium")}
}

\section{Performance}{


\code{v0.2} and above of this package will have \code{wt=csv} as the default.
This  should give significant performance improvement over the previous
default of \code{wt=json}, which pulled down json, parsed to an R list,
then to a data.frame. With \code{wt=csv}, we pull down csv, and read that
in directly to a data.frame.

The http library we use, \pkg{crul}, sets gzip compression header by
default. As long as compression is used server side, you're good to go on
compression, which should be a good peformance boost. See
https://wiki.apache.org/solr/SolrPerformanceFactors#Query_Response_Compression
for notes on how to enable compression.

There are other notes about Solr performance at
https://wiki.apache.org/solr/SolrPerformanceFactors that can be
used server side/in your Solr config, but aren't things to tune here in
this R client.

Let us know if there's any further performance improvements we can make.
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commit.R
\name{commit}
\alias{commit}
\title{Commit}
\usage{
commit(
  conn,
  name,
  expunge_deletes = FALSE,
  wait_searcher = TRUE,
  soft_commit = FALSE,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) A collection or core name. Required.}

\item{expunge_deletes}{merge segments with deletes away. Default: \code{FALSE}}

\item{wait_searcher}{block until a new searcher is opened and registered as
the main query searcher, making the changes visible. Default: \code{TRUE}}

\item{soft_commit}{perform a soft commit - this will refresh the 'view' of
the index in a more performant manner, but without "on-disk" guarantees.
Default: \code{FALSE}}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Commit
}
\examples{
\dontrun{
(conn <- SolrClient$new())

conn$commit("gettingstarted")
conn$commit("gettingstarted", wait_searcher = FALSE)

# get xml back
conn$commit("gettingstarted", wt = "xml")
## raw xml
conn$commit("gettingstarted", wt = "xml", raw = TRUE)
}
}
\references{
<>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_rename.R
\name{core_rename}
\alias{core_rename}
\title{Rename a core}
\usage{
core_rename(conn, name, other, async = NULL, raw = FALSE, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{other}{(character) The new name of the core. Required.}

\item{async}{(character) Request ID to track this action which will be
processed asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Rename a core
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#   bin/solr start -e schemaless
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# Status of particular cores
path <- "~/solr-8.2.0/server/solr/testcore/conf"
dir.create(path, recursive = TRUE)
files <- list.files(
"~/solr-8.2.0/server/solr/configsets/sample_techproducts_configs/conf/",
full.names = TRUE)
invisible(file.copy(files, path, recursive = TRUE))
conn$core_create("testcore") # or create in CLI: bin/solr create -c testcore

# rename
conn$core_rename("testcore", "newtestcore")
## status
conn$core_status("testcore") # core missing
conn$core_status("newtestcore", FALSE) # not missing

# cleanup
conn$core_unload("newtestcore")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_delete.R
\name{collection_delete}
\alias{collection_delete}
\title{Add a collection}
\usage{
collection_delete(conn, name, raw = FALSE, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \code{\link[crul]{HttpClient}}}
}
\description{
Add a collection
}
\examples{
\dontrun{
(conn <- SolrClient$new())

if (!conn$collection_exists("helloWorld")) {
  conn$collection_create(name = "helloWorld")
}

collection_delete(conn, name = "helloWorld")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_exists.R
\name{core_exists}
\alias{core_exists}
\title{Check if a core exists}
\usage{
core_exists(conn, name, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A single boolean, \code{TRUE} or \code{FALSE}
}
\description{
Check if a core exists
}
\details{
Simply calls \code{\link[=core_status]{core_status()}} internally
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#   bin/solr start -e schemaless
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or create as below

# connect
(conn <- SolrClient$new())

# exists
conn$core_exists("gettingstarted")

# doesn't exist
conn$core_exists("hhhhhh")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_list.R
\name{collection_list}
\alias{collection_list}
\title{List collections}
\usage{
collection_list(conn, raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
List collections
}
\examples{
\dontrun{
(conn <- SolrClient$new())

conn$collection_list()
conn$collection_list()$collections
collection_list(conn)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_create.R
\name{core_create}
\alias{core_create}
\title{Create a core}
\usage{
core_create(
  conn,
  name,
  instanceDir = NULL,
  config = NULL,
  schema = NULL,
  dataDir = NULL,
  configSet = NULL,
  collection = NULL,
  shard = NULL,
  async = NULL,
  raw = FALSE,
  callopts = list(),
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{instanceDir}{(character) Path to instance directory}

\item{config}{(character) Path to config file}

\item{schema}{(character) Path to schema file}

\item{dataDir}{(character) Name of the data directory relative to
instanceDir.}

\item{configSet}{(character) Name of the configset to use for this core.
For more information, see
https://lucene.apache.org/solr/guide/6_6/config-sets.html}

\item{collection}{(character) The name of the collection to which this core
belongs. The default is the name of the \verb{core.collection.<param>=<value>}
causes a property of \verb{<param>=<value>} to be set if a new collection is being
created. Use \verb{collection.configName=<configname>} to point to the
configuration for a new collection.}

\item{shard}{(character) The shard id this core represents. Normally you
want to be auto-assigned a shard id.}

\item{async}{(character) Request ID to track this action which will be
processed asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/6_6/defining-core-properties.html)}
}
\description{
Create a core
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#   bin/solr start -e schemaless
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or create as below

# connect
(conn <- SolrClient$new())

# Create a core
path <- "~/solr-8.2.0/server/solr/newcore/conf"
dir.create(path, recursive = TRUE)
files <- list.files("~/solr-8.2.0/server/solr/configsets/sample_techproducts_configs/conf/",
full.names = TRUE)
invisible(file.copy(files, path, recursive = TRUE))
conn$core_create(name = "newcore", instanceDir = "newcore",
  configSet = "sample_techproducts_configs")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_atomic_json.R
\name{update_atomic_json}
\alias{update_atomic_json}
\title{Atomic updates with JSON data}
\usage{
update_atomic_json(conn, body, name, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{body}{(character) JSON as a character string}

\item{name}{(character) Name of the core or collection}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Atomic updates to parts of Solr documents
}
\examples{
\dontrun{
# start Solr in Cloud mode: bin/solr start -e cloud -noprompt

# connect
(conn <- SolrClient$new())

# create a collection
if (!conn$collection_exists("books")) {
  conn$collection_delete("books")
  conn$collection_create("books")
}

# Add documents
file <- system.file("examples", "books2.json", package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_json(file, "books")

# get a document
conn$get(ids = 343334534545, "books")

# atomic update
body <- '[{
 "id": "343334534545",
 "genre_s": {"set": "mystery" },
 "pages_i": {"inc": 1 }
}]'
conn$update_atomic_json(body, "books")

# get the document again
conn$get(ids = 343334534545, "books")

# another atomic update
body <- '[{
 "id": "343334534545",
 "price": {"remove": "12.5" }
}]'
conn$update_atomic_json(body, "books")

# get the document again
conn$get(ids = 343334534545, "books")
}
}
\references{
https://lucene.apache.org/solr/guide/7_0/updating-parts-of-documents.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_get.R
\name{solr_get}
\alias{solr_get}
\title{Real time get}
\usage{
solr_get(conn, ids, name, fl = NULL, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{ids}{Document IDs, one or more in a vector or list}

\item{name}{(character) A collection or core name. Required.}

\item{fl}{Fields to return, can be a character vector like
\code{c('id', 'title')}, or a single character vector with one or more
comma separated names, like \code{'id,title'}}

\item{wt}{(character) One of json (default) or xml. Data type returned.
If json, uses \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses
\code{\link[xml2:read_xml]{xml2::read_xml()}} to parse.}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get documents by id
}
\details{
We use json internally as data interchange format for this function.
}
\examples{
\dontrun{
(cli <- SolrClient$new())

# add some documents first
ss <- list(list(id = 1, price = 100), list(id = 2, price = 500))
add(cli, ss, name = "gettingstarted")

# Now, get documents by id
solr_get(cli, ids = 1, "gettingstarted")
solr_get(cli, ids = 2, "gettingstarted")
solr_get(cli, ids = c(1, 2), "gettingstarted")
solr_get(cli, ids = "1,2", "gettingstarted")

# Get raw JSON
solr_get(cli, ids = 1, "gettingstarted", raw = TRUE, wt = "json")
solr_get(cli, ids = 1, "gettingstarted", raw = TRUE, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_deletereplica.R
\name{collection_deletereplica}
\alias{collection_deletereplica}
\title{Delete a replica}
\usage{
collection_deletereplica(
  conn,
  name,
  shard = NULL,
  replica = NULL,
  onlyIfDown = FALSE,
  raw = FALSE,
  callopts = list(),
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) Required. The name of the collection.}

\item{shard}{(character) Required. The name of the shard that includes the replica to
be removed.}

\item{replica}{(character) Required. The name of the replica to remove.}

\item{onlyIfDown}{(logical) When \code{TRUE} will not take any action if the replica
is active. Default: \code{FALSE}}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for details on
supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Delete a replica from a given collection and shard. If the
corresponding core is up and running the core is unloaded and the entry is
removed from the clusterstate. If the node/core is down , the entry is taken
off the clusterstate and if the core comes up later it is automatically
unregistered.
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("foobar2")) {
  conn$collection_create(name = "foobar2", maxShardsPerNode = 2)
}

# status
conn$collection_clusterstatus()$cluster$collections$foobar2$shards$shard1

# add replica
conn$collection_addreplica(name = "foobar2", shard = "shard1")

# delete replica
## get replica name
nms <- names(conn$collection_clusterstatus()$cluster$collections$foobar2$shards$shard1$replicas)
conn$collection_deletereplica(name = "foobar2", shard = "shard1", replica = nms[1])

# status again
conn$collection_clusterstatus()$cluster$collections$foobar2$shards$shard1
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{makemultiargs}
\alias{makemultiargs}
\title{Function to make make multiple args of the same name from a
single input with length > 1}
\usage{
makemultiargs(x)
}
\arguments{
\item{x}{Value}
}
\description{
Function to make make multiple args of the same name from a
single input with length > 1
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_migrate.R
\name{collection_migrate}
\alias{collection_migrate}
\title{Migrate documents to another collection}
\usage{
collection_migrate(
  conn,
  name,
  target.collection,
  split.key,
  forward.timeout = NULL,
  async = NULL,
  raw = FALSE,
  callopts = list()
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{target.collection}{(character) Required. The name of the target collection
to which documents will be migrated}

\item{split.key}{(character) Required. The routing key prefix. For example, if
uniqueKey is a!123, then you would use split.key=a!}

\item{forward.timeout}{(integer) The timeout (seconds), until which write requests
made to the source collection for the given \code{split.key} will be forwarded to the
target shard. Default: 60}

\item{async}{(character) Request ID to track this action which will be processed
asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Migrate documents to another collection
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("migrate_from")) {
  conn$collection_create(name = "migrate_from")
  # OR: bin/solr create -c migrate_from
}

# create another collection
if (!conn$collection_exists("migrate_to")) {
  conn$collection_create(name = "migrate_to")
  # OR bin/solr create -c migrate_to
}

# add some documents
file <- system.file("examples", "books.csv", package = "solrium")
x <- read.csv(file, stringsAsFactors = FALSE)
conn$add(x, "migrate_from")

# migrate some documents from one collection to the other
## FIXME - not sure if this is actually working....
# conn$collection_migrate("migrate_from", "migrate_to", split.key = "05535")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_swap.R
\name{core_swap}
\alias{core_swap}
\title{Swap a core}
\usage{
core_swap(conn, name, other, async = NULL, raw = FALSE, callopts = list())
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{other}{(character) The name of one of the cores to be swapped.
Required.}

\item{async}{(character) Request ID to track this action which will be
processed asynchronously}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
SWAP atomically swaps the names used to access two existing
Solr cores. This can be used to swap new content into production. The
prior core remains available and can be swapped back, if necessary. Each
core will be known by the name of the other, after the swap
}
\details{
Do not use \code{core_swap} with a SolrCloud node. It is not
supported and can result in the core being unusable. We'll try to stop
you if you try.
}
\examples{
\dontrun{
# start Solr with Schemaless mode via the schemaless eg:
#   bin/solr start -e schemaless
# you can create a new core like: bin/solr create -c corename
# where <corename> is the name for your core - or creaate as below

# connect
(conn <- SolrClient$new())

# Swap a core
## First, create two cores
conn$core_create("swapcoretest1")
# - or create on CLI: bin/solr create -c swapcoretest1
conn$core_create("swapcoretest2")
# - or create on CLI: bin/solr create -c swapcoretest2

## check status
conn$core_status("swapcoretest1", FALSE)
conn$core_status("swapcoretest2", FALSE)

## swap core
conn$core_swap("swapcoretest1", "swapcoretest2")

## check status again
conn$core_status("swapcoretest1", FALSE)
conn$core_status("swapcoretest2", FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/schema.R
\name{schema}
\alias{schema}
\title{Get the schema for a collection or core}
\usage{
schema(conn, name, what = "", raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) Name of a collection or core. Required.}

\item{what}{(character) What to retrieve. By default, we retrieve the entire
schema. Options include: fields, dynamicfields, fieldtypes, copyfields, name,
version, uniquekey, similarity, "solrqueryparser/defaultoperator"}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get the schema for a collection or core
}
\examples{
\dontrun{
# start Solr, in your CLI, run: `bin/solr start -e cloud -noprompt`
# after that, if you haven't run `bin/post -c gettingstarted docs/` yet, do so

# connect: by default we connect to localhost, port 8983
(cli <- SolrClient$new())

# get the schema for the gettingstarted index
schema(cli, name = "gettingstarted")

# Get parts of the schema
schema(cli, name = "gettingstarted", "fields")
schema(cli, name = "gettingstarted", "dynamicfields")
schema(cli, name = "gettingstarted", "fieldtypes")
schema(cli, name = "gettingstarted", "copyfields")
schema(cli, name = "gettingstarted", "name")
schema(cli, name = "gettingstarted", "version")
schema(cli, name = "gettingstarted", "uniquekey")
schema(cli, name = "gettingstarted", "similarity")
schema(cli, name = "gettingstarted", "solrqueryparser/defaultoperator")

# get raw data
schema(cli, name = "gettingstarted", "similarity", raw = TRUE)
schema(cli, name = "gettingstarted", "uniquekey", raw = TRUE)

# start Solr in Schemaless mode: bin/solr start -e schemaless
# schema(cli, "gettingstarted")

# start Solr in Standalone mode: bin/solr start
# then add a core: bin/solr create -c helloWorld
# schema(cli, "helloWorld")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add}
\alias{add}
\title{Add documents from R objects}
\usage{
add(
  x,
  conn,
  name,
  commit = TRUE,
  commit_within = NULL,
  overwrite = TRUE,
  boost = NULL,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{x}{Documents, either as rows in a data.frame, or a list.}

\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) A collection or core name. Required.}

\item{commit}{(logical) If \code{TRUE}, documents immediately searchable.
Default: \code{TRUE}}

\item{commit_within}{(numeric) Milliseconds to commit the change, the
document will be added within that time. Default: NULL}

\item{overwrite}{(logical) Overwrite documents with matching keys.
Default: \code{TRUE}}

\item{boost}{(numeric) Boost factor. Default: NULL}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to
parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Add documents from R objects
}
\details{
Works for Collections as well as Cores (in SolrCloud and Standalone
modes, respectively)
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create the boooks collection
if (!collection_exists(conn, "books")) {
  collection_create(conn, name = "books", numShards = 1)
}

# Documents in a list
ss <- list(list(id = 1, price = 100), list(id = 2, price = 500))
add(ss, conn, name = "books")
conn$get(c(1, 2), "books")

# Documents in a data.frame
## Simple example
df <- data.frame(id = c(67, 68), price = c(1000, 500000000))
add(df, conn, "books")
df <- data.frame(id = c(77, 78), price = c(1, 2.40))
add(df, conn, "books")

## More complex example, get file from package examples
# start Solr in Schemaless mode first: bin/solr start -e schemaless
file <- system.file("examples", "books.csv", package = "solrium")
x <- read.csv(file, stringsAsFactors = FALSE)
class(x)
head(x)
if (!collection_exists(conn, "mybooks")) {
  collection_create(conn, name = "mybooks", numShards = 2)
}
add(x, conn, "mybooks")

# Use modifiers
add(x, conn, "mybooks", commit_within = 5000)

# Get back XML instead of a list
ss <- list(list(id = 1, price = 100), list(id = 2, price = 500))
# parsed XML
add(ss, conn, name = "books", wt = "xml")
# raw XML
add(ss, conn, name = "books", wt = "xml", raw = TRUE)
}
}
\seealso{
\code{\link{update_json}}, \code{\link{update_xml}},
\code{\link{update_csv}} for adding documents from files
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_addrole.R
\name{collection_addrole}
\alias{collection_addrole}
\title{Add a role to a node}
\usage{
collection_addrole(conn, role = "overseer", node, raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{role}{(character) Required. The name of the role. The only supported role
as of now is overseer (set as default).}

\item{node}{(character) Required. The name of the node. It is possible to assign a
role even before that node is started.}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Assign a role to a given node in the cluster. The only supported role
as of 4.7 is 'overseer' . Use this API to dedicate a particular node as Overseer.
Invoke it multiple times to add more nodes. This is useful in large clusters where
an Overseer is likely to get overloaded . If available, one among the list of
nodes which are assigned the 'overseer' role would become the overseer. The
system would assign the role to any other node if none of the designated nodes
are up and running
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# get list of nodes
nodes <- conn$collection_clusterstatus()$cluster$live_nodes
collection_addrole(conn, node = nodes[1])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr_all.r
\name{solr_all}
\alias{solr_all}
\title{All purpose search}
\usage{
solr_all(
  conn,
  name = NULL,
  params = NULL,
  body = NULL,
  callopts = list(),
  raw = FALSE,
  parsetype = "df",
  concat = ",",
  optimizeMaxRows = TRUE,
  minOptimizedRows = 50000L,
  progress = NULL,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{Name of a collection or core. Or leave as \code{NULL} if not needed.}

\item{params}{(list) a named list of parameters, results in a GET request
as long as no body parameters given}

\item{body}{(list) a named list of parameters, if given a POST request
will be performed}

\item{callopts}{Call options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{raw}{(logical) If TRUE, returns raw data in format specified by wt param}

\item{parsetype}{(character) One of 'list' or 'df'}

\item{concat}{(character) Character to concatenate elements of longer than length 1.
Note that this only works reliably when data format is json (wt='json'). The parsing
is more complicated in XML format, but you can do that on your own.}

\item{optimizeMaxRows}{(logical) If \code{TRUE}, then rows parameter will be
adjusted to the number of returned results by the same constraints.
It will only be applied if rows parameter is higher
than \code{minOptimizedRows}. Default: \code{TRUE}}

\item{minOptimizedRows}{(numeric) used by \code{optimizedMaxRows} parameter,
the minimum optimized rows. Default: 50000}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}. only supports
\code{httr::progress} for now. See the README for an example.}

\item{...}{Further args to be combined into query}
}
\value{
XML, JSON, a list, or data.frame
}
\description{
Includes documents, facets, groups, mlt, stats, and highlights
}
\section{Parameters}{

\itemize{
\item q Query terms, defaults to '\emph{:}', or everything.
\item sort Field to sort on. You can specify ascending (e.g., score desc) or
descending (e.g., score asc), sort by two fields (e.g., score desc, price asc),
or sort by a function (e.g., sum(x_f, y_f) desc, which sorts by the sum of
x_f and y_f in a descending order).
\item start Record to start at, default to beginning.
\item rows Number of records to return. Default: 10.
\item pageDoc If you expect to be paging deeply into the results (say beyond page 10,
assuming rows=10) and you are sorting by score, you may wish to add the pageDoc
and pageScore parameters to your request. These two parameters tell Solr (and Lucene)
what the last result (Lucene internal docid and score) of the previous page was,
so that when scoring the query for the next set of pages, it can ignore any results
that occur higher than that item. To get the Lucene internal doc id, you will need
to add \code{docid} to the &fl list.
\item pageScore See pageDoc notes.
\item fq Filter query, this does not affect the search, only what gets returned.
This parameter can accept multiple items in a lis or vector. You can't pass more than
one parameter of the same name, so we get around it by passing multiple queries
and we parse internally
\item fl Fields to return, can be a character vector like \code{c('id', 'title')},
or a single character vector with one or more comma separated names, like
\code{'id,title'}
\item defType Specify the query parser to use with this request.
\item timeAllowed The time allowed for a search to finish. This value only applies
to the search and not to requests in general. Time is in milliseconds. Values \code{<= 0}
mean no time restriction. Partial results may be returned (if there are any).
\item qt Which query handler used. Options: dismax, others?
\item NOW Set a fixed time for evaluating Date based expresions
\item TZ Time zone, you can override the default.
\item echoHandler If \code{TRUE}, Solr places the name of the handle used in the
response to the client for debugging purposes. Default:
\item echoParams The echoParams parameter tells Solr what kinds of Request
parameters should be included in the response for debugging purposes, legal values
include:
\itemize{
\item none - don't include any request parameters for debugging
\item explicit - include the parameters explicitly specified by the client in the request
\item all - include all parameters involved in this request, either specified explicitly
by the client, or implicit because of the request handler configuration.
}
\item wt (character) One of json, xml, or csv. Data type returned, defaults
to 'csv'. If json, uses \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml,
uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to parse. If csv, uses \code{\link[=read.table]{read.table()}} to parse.
\code{wt=csv} gives the fastest performance at least in all the cases we have
tested in, thus it's the default value for \code{wt}
}
}

\examples{
\dontrun{
# connect
(cli <- SolrClient$new(host = "api.plos.org", path = "search", port = NULL))

solr_all(cli, params = list(q='*:*', rows=2, fl='id'))

# facets
solr_all(cli, params = list(q='*:*', rows=2, fl='id', facet="true",
  facet.field="journal"))

# mlt
solr_all(cli, params = list(q='ecology', rows=2, fl='id', mlt='true',
  mlt.count=2, mlt.fl='abstract'))

# facets and mlt
solr_all(cli, params = list(q='ecology', rows=2, fl='id', facet="true",
  facet.field="journal", mlt='true', mlt.count=2, mlt.fl='abstract'))

# stats
solr_all(cli, params = list(q='ecology', rows=2, fl='id', stats='true',
  stats.field='counter_total_all'))

# facets, mlt, and stats
solr_all(cli, params = list(q='ecology', rows=2, fl='id', facet="true",
  facet.field="journal", mlt='true', mlt.count=2, mlt.fl='abstract',
  stats='true', stats.field='counter_total_all'))

# group
solr_all(cli, params = list(q='ecology', rows=2, fl='id', group='true',
 group.field='journal', group.limit=3))

# facets, mlt, stats, and groups
solr_all(cli, params = list(q='ecology', rows=2, fl='id', facet="true",
 facet.field="journal", mlt='true', mlt.count=2, mlt.fl='abstract',
 stats='true', stats.field='counter_total_all', group='true',
 group.field='journal', group.limit=3))

# using wt = xml
solr_all(cli, params = list(q='*:*', rows=50, fl=c('id','score'),
  fq='doc_type:full', wt="xml"), raw=TRUE)
}
}
\references{
See https://lucene.apache.org/solr/guide/8_2/searching.html for
more information.
}
\seealso{
\code{\link[=solr_highlight]{solr_highlight()}}, \code{\link[=solr_facet]{solr_facet()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_atomic_xml.R
\name{update_atomic_xml}
\alias{update_atomic_xml}
\title{Atomic updates with XML data}
\usage{
update_atomic_xml(conn, body, name, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{body}{(character) XML as a character string}

\item{name}{(character) Name of the core or collection}

\item{wt}{(character) One of json (default) or xml. If json, uses
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} to parse. If xml, uses \code{\link[xml2:read_xml]{xml2::read_xml()}} to parse}

\item{raw}{(logical) If \code{TRUE}, returns raw data in format specified by
\code{wt} param}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Atomic updates to parts of Solr documents
}
\examples{
\dontrun{
# start Solr in Cloud mode: bin/solr start -e cloud -noprompt

# connect
(conn <- SolrClient$new())

# create a collection
if (!conn$collection_exists("books")) {
  conn$collection_delete("books")
  conn$collection_create("books")
}

# Add documents
file <- system.file("examples", "books.xml", package = "solrium")
cat(readLines(file), sep = "\n")
conn$update_xml(file, "books")

# get a document
conn$get(ids = '978-0641723445', "books", wt = "xml")

# atomic update
body <- '
<add>
 <doc>
   <field name="id">978-0641723445</field>
   <field name="genre_s" update="set">mystery</field>
   <field name="pages_i" update="inc">1</field>
 </doc>
</add>'
conn$update_atomic_xml(body, name="books")

# get the document again
conn$get(ids = '978-0641723445', "books", wt = "xml")

# another atomic update
body <- '
<add>
 <doc>
   <field name="id">978-0641723445</field>
   <field name="price" update="remove">12.5</field>
 </doc>
</add>'
conn$update_atomic_xml(body, "books")

# get the document again
conn$get(ids = '978-0641723445', "books", wt = "xml")
}
}
\references{
https://lucene.apache.org/solr/guide/7_0/updating-parts-of-documents.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_reload.R
\name{collection_reload}
\alias{collection_reload}
\title{Reload a collection}
\usage{
collection_reload(conn, name, raw = FALSE, callopts)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{callopts}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Reload a collection
}
\examples{
\dontrun{
(conn <- SolrClient$new())

if (!conn$collection_exists("helloWorld")) {
  conn$collection_create(name = "helloWorld")
}

conn$collection_reload(name = "helloWorld")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collections.R
\name{collections}
\alias{collections}
\alias{cores}
\title{List collections or cores}
\usage{
collections(conn, ...)

cores(conn, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
character vector
}
\description{
List collections or cores
}
\details{
Calls \code{\link[=collection_list]{collection_list()}} or \code{\link[=core_status]{core_status()}} internally,
and parses out names for you.
}
\examples{
\dontrun{
# connect
(conn <- SolrClient$new())

# list collections
conn$collection_list()
collections(conn)

# list cores
conn$core_status()
cores(conn)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_rebalanceleaders.R
\name{collection_rebalanceleaders}
\alias{collection_rebalanceleaders}
\title{Rebalance leaders}
\usage{
collection_rebalanceleaders(
  conn,
  name,
  maxAtOnce = NULL,
  maxWaitSeconds = NULL,
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) The name of the core to be created. Required}

\item{maxAtOnce}{(integer) The maximum number of reassignments to have queue
up at once. Values <=0 are use the default value Integer.MAX_VALUE. When
this number is reached, the process waits for one or more leaders to be
successfully assigned before adding more to the queue.}

\item{maxWaitSeconds}{(integer) Timeout value when waiting for leaders to
be reassigned. NOTE: if maxAtOnce is less than the number of reassignments
that will take place, this is the maximum interval that any single wait for
at least one reassignment. For example, if 10 reassignments are to take
place and maxAtOnce is 1 and maxWaitSeconds is 60, the upper bound on the
time that the command may wait is 10 minutes. Default: 60}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{You can pass in parameters like \code{property.name=value}    to set
core property name to value. See the section Defining core.properties for
details on supported properties and values.
(https://lucene.apache.org/solr/guide/8_2/defining-core-properties.html)}
}
\description{
Reassign leaders in a collection according to the preferredLeader
property across active nodes
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("mycollection2")) {
  conn$collection_create(name = "mycollection2")
  # OR: bin/solr create -c mycollection2
}

# balance preferredLeader property
conn$collection_balanceshardunique("mycollection2", property = "preferredLeader")

# balance preferredLeader property
conn$collection_rebalanceleaders("mycollection2")

# examine cluster status
conn$collection_clusterstatus()$cluster$collections$mycollection2
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_deleteshard.R
\name{collection_deleteshard}
\alias{collection_deleteshard}
\title{Delete a shard}
\usage{
collection_deleteshard(conn, name, shard, raw = FALSE, ...)
}
\arguments{
\item{conn}{A solrium connection object, see \link{SolrClient}}

\item{name}{(character) Required. The name of the collection that includes the shard
to be deleted}

\item{shard}{(character) Required. The name of the shard to be deleted}

\item{raw}{(logical) If \code{TRUE}, returns raw data}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Deleting a shard will unload all replicas of the shard and remove
them from clusterstate.json. It will only remove shards that are inactive, or
which have no range given for custom sharding.
}
\examples{
\dontrun{
(conn <- SolrClient$new())

# create collection
if (!conn$collection_exists("buffalo")) {
  conn$collection_create(name = "buffalo")
  # OR: bin/solr create -c buffalo
}

# find shard names
names(conn$collection_clusterstatus()$cluster$collections$buffalo$shards)

# split a shard by name
collection_splitshard(conn, name = "buffalo", shard = "shard1")

# now we have three shards
names(conn$collection_clusterstatus()$cluster$collections$buffalo$shards)

# delete shard
conn$collection_deleteshard(name = "buffalo", shard = "shard1_1")
}
}
