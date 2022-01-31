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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

[![Build
Status](https://travis-ci.org/ropensci/tidypmc.svg?branch=master)](https://travis-ci.org/ropensci/tidypmc)
[![Coverage
status](https://codecov.io/gh/ropensci/tidypmc/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tidypmc?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/tidypmc)](https://cran.r-project.org/package=tidypmc)
[![Downloads](https://cranlogs.r-pkg.org/badges/tidypmc)](https://CRAN.R-project.org/package=tidypmc)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/tidypmc?color=orange)](https://CRAN.R-project.org/package=tidypmc)

# tidypmc

The [Open Access subset](https://europepmc.org/downloads/openaccess) of
[Pubmed Central](https://europepmc.org) (PMC) includes 2.5 million
articles from biomedical and life sciences journals. The full text XML
files are freely available for text mining from the [REST
service](https://europepmc.org/RestfulWebService) or [FTP
site](https://europepmc.org/ftp/oa/) but can be challenging to parse.
For example, section tags are nested to arbitrary depths, formulas and
tables may return incomprehensible text blobs and superscripted
references are pasted at the end of words. The functions in the
`tidypmc` package are intended to return readable text and maintain the
document structure, so gene names and other terms can be associated with
specific sections, paragraphs, sentences or table rows.

## Installation

Use [remotes](https://github.com/r-lib/remotes) to install the package.

``` r
remotes::install_github("ropensci/tidypmc")
```

## Load XML

Download a single XML document like
[PMC2231364](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2231364/fullTextXML)
from the [REST service](https://europepmc.org/RestfulWebService) using
the `pmc_xml` function.

``` r
library(tidypmc)
library(tidyverse)
doc <- pmc_xml("PMC2231364")
doc
#  {xml_document}
#  <article article-type="research-article" xmlns:xlink="http://www.w3.org/1999/xlink">
#  [1] <front>\n  <journal-meta>\n    <journal-id journal-id-type="nlm-ta">BMC Microbiol</journal-id ...
#  [2] <body>\n  <sec>\n    <title>Background</title>\n    <p><italic>Yersinia pestis </italic>is th ...
#  [3] <back>\n  <ack>\n    <sec>\n      <title>Acknowledgements</title>\n      <p>We thank Dr. Chen ...
```

The [europepmc](https://github.com/ropensci/europepmc) package includes
additional functions to search PMC and download full text. Be sure to
include the `OPEN_ACCESS` field in the search since these are the only
articles with full text XML available.

``` r
library(europepmc)
yp <- epmc_search("title:(Yersinia pestis virulence) OPEN_ACCESS:Y")
#  19 records found, returning 19
select(yp, pmcid, pubYear, title) %>%
  print(n=5)
#  # A tibble: 19 x 3
#    pmcid      pubYear title                                                                          
#    <chr>      <chr>   <chr>                                                                          
#  1 PMC5505154 2017    Crystal structure of Yersinia pestis virulence factor YfeA reveals two polyspe…
#  2 PMC3521224 2012    Omics strategies for revealing Yersinia pestis virulence.                      
#  3 PMC2704395 2009    Involvement of the post-transcriptional regulator Hfq in Yersinia pestis virul…
#  4 PMC2736372 2009    The NlpD lipoprotein is a novel Yersinia pestis virulence factor essential for…
#  5 PMC3109262 2011    A comprehensive study on the role of the Yersinia pestis virulence markers in …
#  # … with 14 more rows
```

Save all 19 results to a list of XML documents using the `epmc_ftxt` or
`pmc_xml` function.

``` r
docs <- map(yp$pmcid, epmc_ftxt)
```

See the [PMC FTP
vignette](https://github.com/ropensci/tidypmc/blob/master/vignettes/pmcftp.md)
for details on parsing the large XML files on the [FTP
site](https://europepmc.org/ftp/oa/) with 10,000 articles each.

## Parse XML

The package includes five functions to parse the
`xml_document`.

| R function      | Description                                                                 |
| :-------------- | :-------------------------------------------------------------------------- |
| `pmc_text`      | Split section paragraphs into sentences with full path to subsection titles |
| `pmc_caption`   | Split figure, table and supplementary material captions into sentences      |
| `pmc_table`     | Convert table nodes into a list of tibbles                                  |
| `pmc_reference` | Format references cited into a tibble                                       |
| `pmc_metadata`  | List journal and article metadata in front node                             |

The `pmc_text` function uses the
[tokenizers](https://lincolnmullen.com/software/tokenizers/) package to
split section paragraphs into sentences. The function also removes any
tables, figures or formulas that are nested within paragraph tags,
replaces superscripted references with brackets, adds carets and
underscores to other superscripts and subscripts and includes the full
path to the subsection title.

``` r
txt <- pmc_text(doc)
#  Note: removing disp-formula nested in sec/p tag
txt
#  # A tibble: 194 x 4
#     section    paragraph sentence text                                                                         
#     <chr>          <int>    <int> <chr>                                                                        
#   1 Title              1        1 Comparative transcriptomics in Yersinia pestis: a global view of environment…
#   2 Abstract           1        1 Environmental modulation of gene expression in Yersinia pestis is critical f…
#   3 Abstract           1        2 Using cDNA microarray technology, we have analyzed the global gene expressio…
#   4 Abstract           2        1 To provide us with a comprehensive view of environmental modulation of globa…
#   5 Abstract           2        2 Almost all known virulence genes of Y. pestis were differentially regulated …
#   6 Abstract           2        3 Clustering enabled us to functionally classify co-expressed genes, including…
#   7 Abstract           2        4 Collections of operons were predicted from the microarray data, and some of …
#   8 Abstract           2        5 Several regulatory DNA motifs, probably recognized by the regulatory protein…
#   9 Abstract           3        1 The comparative transcriptomics analysis we present here not only benefits o…
#  10 Background         1        1 Yersinia pestis is the etiological agent of plague, alternatively growing in…
#  # … with 184 more rows
count(txt, section, sort=TRUE)
#  # A tibble: 21 x 2
#     section                                                                                                   n
#     <chr>                                                                                                 <int>
#   1 Results and Discussion; Clustering analysis and functional classification of co-expressed gene clust…    22
#   2 Background                                                                                               20
#   3 Results and Discussion; Virulence genes in response to multiple environmental stresses                   20
#   4 Methods; Collection of microarray expression data                                                        17
#   5 Results and Discussion; Computational discovery of regulatory DNA motifs                                 16
#   6 Methods; Gel mobility shift analysis of Fur binding                                                      13
#   7 Results and Discussion; Verification of predicted operons by RT-PCR                                      10
#   8 Abstract                                                                                                  8
#   9 Methods; Discovery of regulatory DNA motifs                                                               8
#  10 Methods; Clustering analysis                                                                              7
#  # … with 11 more rows
```

Load the [tidytext](https://www.tidytextmining.com/) package for further
text processing.

``` r
library(tidytext)
x1 <- unnest_tokens(txt, word, text) %>%
  anti_join(stop_words) %>%
  filter(!word %in% 1:100)
#  Joining, by = "word"
filter(x1, str_detect(section, "^Results"))
#  # A tibble: 1,269 x 4
#     section                paragraph sentence word         
#     <chr>                      <int>    <int> <chr>        
#   1 Results and Discussion         1        1 comprehensive
#   2 Results and Discussion         1        1 analysis     
#   3 Results and Discussion         1        1 sets         
#   4 Results and Discussion         1        1 microarray   
#   5 Results and Discussion         1        1 expression   
#   6 Results and Discussion         1        1 data         
#   7 Results and Discussion         1        1 dissect      
#   8 Results and Discussion         1        1 bacterial    
#   9 Results and Discussion         1        1 adaptation   
#  10 Results and Discussion         1        1 environments 
#  # … with 1,259 more rows
filter(x1, str_detect(section, "^Results")) %>%
  count(word, sort = TRUE)
#  # A tibble: 595 x 2
#     word           n
#     <chr>      <int>
#   1 genes         45
#   2 cluster       24
#   3 expression    21
#   4 pestis        21
#   5 data          19
#   6 dna           15
#   7 gene          15
#   8 figure        13
#   9 fur           12
#  10 operons       12
#  # … with 585 more rows
```

The `pmc_table` function formats tables by collapsing multiline headers,
expanding rowspan and colspan attributes and adding subheadings into a
new column.

``` r
tbls <- pmc_table(doc)
#  Parsing 4 tables
#  Adding footnotes to Table 1
map_int(tbls, nrow)
#  Table 1 Table 2 Table 3 Table 4 
#       39      23       4      34
tbls[[1]]
#  # A tibble: 39 x 5
#     subheading              `Potential operon (r va… `Gene ID`   `Putative or predicted functi… `Reference (s)`
#     <chr>                   <chr>                    <chr>       <chr>                          <chr>          
#   1 Iron uptake or heme sy… yfeABCD operon* (r > 0.… YPO2439-24… Transport/binding chelated ir… yfeABCD [54]   
#   2 Iron uptake or heme sy… hmuRSTUV operon (r > 0.… YPO0279-02… Transport/binding hemin        hmuRSTUV [55]  
#   3 Iron uptake or heme sy… ysuJIHG* (r > 0.95)      YPO1529-15… Iron uptake                    -              
#   4 Iron uptake or heme sy… sufABCDS* (r > 0.90)     YPO2400-24… Iron-regulated Fe-S cluster a… -              
#   5 Iron uptake or heme sy… YPO1854-1856* (r > 0.97) YPO1854-18… Iron uptake or heme synthesis? -              
#   6 Sulfur metabolism       tauABCD operon (r > 0.9… YPO0182-01… Transport/binding taurine      tauABCD [56]   
#   7 Sulfur metabolism       ssuEADCB operon (r > 0.… YPO3623-36… Sulphur metabolism             ssu operon [57]
#   8 Sulfur metabolism       cys operon (r > 0.92)    YPO3010-30… Cysteine synthesis             -              
#   9 Sulfur metabolism       YPO1317-1319 (r > 0.97)  YPO1317-13… Sulfur metabolism?             -              
#  10 Sulfur metabolism       YPO4109-4111 (r > 0.90)  YPO4109-41… Sulfur metabolism?             -              
#  # … with 29 more rows
```

Use `collapse_rows` to join column names and cell values in a semi-colon
delimited string (and then search using functions in the next section).

``` r
collapse_rows(tbls, na.string="-")
#  # A tibble: 100 x 3
#     table     row text                                                                                         
#     <chr>   <int> <chr>                                                                                        
#   1 Table 1     1 subheading=Iron uptake or heme synthesis; Potential operon (r value)=yfeABCD operon* (r > 0.…
#   2 Table 1     2 subheading=Iron uptake or heme synthesis; Potential operon (r value)=hmuRSTUV operon (r > 0.…
#   3 Table 1     3 subheading=Iron uptake or heme synthesis; Potential operon (r value)=ysuJIHG* (r > 0.95); Ge…
#   4 Table 1     4 subheading=Iron uptake or heme synthesis; Potential operon (r value)=sufABCDS* (r > 0.90); G…
#   5 Table 1     5 subheading=Iron uptake or heme synthesis; Potential operon (r value)=YPO1854-1856* (r > 0.97…
#   6 Table 1     6 subheading=Sulfur metabolism; Potential operon (r value)=tauABCD operon (r > 0.90); Gene ID=…
#   7 Table 1     7 subheading=Sulfur metabolism; Potential operon (r value)=ssuEADCB operon (r > 0.97); Gene ID…
#   8 Table 1     8 subheading=Sulfur metabolism; Potential operon (r value)=cys operon (r > 0.92); Gene ID=YPO3…
#   9 Table 1     9 subheading=Sulfur metabolism; Potential operon (r value)=YPO1317-1319 (r > 0.97); Gene ID=YP…
#  10 Table 1    10 subheading=Sulfur metabolism; Potential operon (r value)=YPO4109-4111 (r > 0.90); Gene ID=YP…
#  # … with 90 more rows
```

The other three `pmc` functions are described in the package
[vignette](https://github.com/ropensci/tidypmc/blob/master/vignettes/tidypmc.md).

## Searching text

There are a few functions to search within the `pmc_text` or collapsed
`pmc_table` output. `separate_text` uses the
[stringr](https://stringr.tidyverse.org/) package to extract any regular
expression or vector of words.

``` r
separate_text(txt, "[ATCGN]{5,}")
#  # A tibble: 9 x 5
#    match        section                         paragraph sentence text                                        
#    <chr>        <chr>                               <int>    <int> <chr>                                       
#  1 ACGCAATCGTT… Results and Discussion; Comput…         2        3 A 16 basepair (bp) box (5'-ACGCAATCGTTTTCNT…
#  2 AAACGTTTNCGT Results and Discussion; Comput…         2        4 It is very similar to the E. coli PurR box …
#  3 TGATAATGATT… Results and Discussion; Comput…         2        5 A 21 bp box (5'-TGATAATGATTATCATTATCA-3') w…
#  4 GATAATGATAA… Results and Discussion; Comput…         2        6 It is a 10-1-10 inverted repeat that resemb…
#  5 TGANNNNNNTC… Results and Discussion; Comput…         2        7 A 15 bp box (5'-TGANNNNNNTCAA-3') was found…
#  6 TTGATN       Results and Discussion; Comput…         2        8 It is a part of the E. coli Fnr box (5'-AAW…
#  7 NATCAA       Results and Discussion; Comput…         2        8 It is a part of the E. coli Fnr box (5'-AAW…
#  8 GTTAATTAA    Results and Discussion; Comput…         3        4 The ArcA regulator can recognize a relative…
#  9 GTTAATTAATGT Results and Discussion; Comput…         3        5 An ArcA-box-like sequence (5'-GTTAATTAATGT-…
```

A few wrappers search pre-defined patterns and add an extra step to
expand matched ranges. `separate_refs` matches references within
brackets using `\\[[0-9, -]+\\]` and expands ranges like `[7-9]`.

``` r
separate_refs(txt)
#  # A tibble: 93 x 6
#        id match section   paragraph sentence text                                                              
#     <dbl> <chr> <chr>         <int>    <int> <chr>                                                             
#   1     1 [1]   Backgrou…         1        1 Yersinia pestis is the etiological agent of plague, alternatively…
#   2     2 [2]   Backgrou…         1        3 To produce a transmissible infection, Y. pestis colonizes the fle…
#   3     3 [3]   Backgrou…         1        9 However, a few bacilli are taken up by tissue macrophages, provid…
#   4     4 [4,5] Backgrou…         1       10 Residence in this niche also facilitates the bacteria's resistanc…
#   5     5 [4,5] Backgrou…         1       10 Residence in this niche also facilitates the bacteria's resistanc…
#   6     6 [6]   Backgrou…         2        1 A DNA microarray is able to determine simultaneous changes in all…
#   7     7 [7-9] Backgrou…         2        2 We and others have measured the gene expression profiles of Y. pe…
#   8     8 [7-9] Backgrou…         2        2 We and others have measured the gene expression profiles of Y. pe…
#   9     9 [7-9] Backgrou…         2        2 We and others have measured the gene expression profiles of Y. pe…
#  10    10 [10]  Backgrou…         2        2 We and others have measured the gene expression profiles of Y. pe…
#  # … with 83 more rows
```

`separate_genes` will find microbial genes like tauD (with a capitalized
4th letter) and expand operons like `tauABCD` into four genes.
`separate_tags` will find and expand locus tag ranges below.

``` r
collapse_rows(tbls, na="-") %>%
  separate_tags("YPO") %>%
  filter(id == "YPO1855")
#  # A tibble: 3 x 5
#    id      match        table    row text                                                                      
#    <chr>   <chr>        <chr>  <int> <chr>                                                                     
#  1 YPO1855 YPO1854-1856 Table…     5 subheading=Iron uptake or heme synthesis; Potential operon (r value)=YPO1…
#  2 YPO1855 YPO1854-1856 Table…    21 subheading=Category C: Hypothetical; Gene ID=YPO1854-1856; Description=Pu…
#  3 YPO1855 YPO1854-YPO… Table…     2 Cluster=Cluster II; Genes or operons for motif discovery=hmuRSTUV, YPO068…
```

See the
[vignette](https://github.com/ropensci/tidypmc/blob/master/vignettes/tidypmc.md)
for more details including code to parse XML documents using the
[xml2](https://github.com/r-lib/xml2) package. The [PMC FTP
vignette](https://github.com/ropensci/tidypmc/blob/master/vignettes/pmcftp.md)
has details on parsing XML files at the Europe PMC [FTP
site](https://europepmc.org/ftp/oa/).

### Community Guidelines

This project is released with a [Contributor Code of
Conduct](CONDUCT.md). By participating in this project you agree to
abide by its terms. Feedback, bug reports, and feature requests are
welcome [here](https://github.com/ropensci/tidypmc/issues).
tidypmc 1.8 (dev)
=========================

### DOCUMENTATION FIXES

  * Added a NEWS.md file (#2)


tidypmc 1.7 (2019-08-01)
=========================

### NEW FEATURES

  * released to CRAN
Introduction to tidypmc
================
Chris Stubben
August 6, 2019

The `tidypmc` package parses XML documents in the Open Access subset of
[Pubmed Central](https://europepmc.org). Download the full text using
`pmc_xml`.

``` r
library(tidypmc)
doc <- pmc_xml("PMC2231364")
doc
#  {xml_document}
#  <article article-type="research-article" xmlns:xlink="http://www.w3.org/1999/xlink">
#  [1] <front>\n  <journal-meta>\n    <journal-id journal-id-type="nlm-ta"> ...
#  [2] <body>\n  <sec>\n    <title>Background</title>\n    <p><italic>Yersi ...
#  [3] <back>\n  <ack>\n    <sec>\n      <title>Acknowledgements</title>\n  ...
```

The package includes five functions to parse the
`xml_document`.

| R function      | Description                                                                 |
| :-------------- | :-------------------------------------------------------------------------- |
| `pmc_text`      | Split section paragraphs into sentences with full path to subsection titles |
| `pmc_caption`   | Split figure, table and supplementary material captions into sentences      |
| `pmc_table`     | Convert table nodes into a list of tibbles                                  |
| `pmc_reference` | Format references cited into a tibble                                       |
| `pmc_metadata`  | List journal and article metadata in front node                             |

`pmc_text` splits paragraphs into sentences and removes any tables,
figures or formulas that are nested within paragraph tags, replaces
superscripted references with brackets, adds carets and underscores to
other superscripts and subscripts and includes the full path to the
subsection title.

``` r
library(dplyr)
txt <- pmc_text(doc)
txt
#  # A tibble: 194 x 4
#     section    paragraph sentence text                                                               
#     <chr>          <int>    <int> <chr>                                                              
#   1 Title              1        1 Comparative transcriptomics in Yersinia pestis: a global view of e…
#   2 Abstract           1        1 Environmental modulation of gene expression in Yersinia pestis is …
#   3 Abstract           1        2 Using cDNA microarray technology, we have analyzed the global gene…
#   4 Abstract           2        1 To provide us with a comprehensive view of environmental modulatio…
#   5 Abstract           2        2 Almost all known virulence genes of Y. pestis were differentially …
#   6 Abstract           2        3 Clustering enabled us to functionally classify co-expressed genes,…
#   7 Abstract           2        4 Collections of operons were predicted from the microarray data, an…
#   8 Abstract           2        5 Several regulatory DNA motifs, probably recognized by the regulato…
#   9 Abstract           3        1 The comparative transcriptomics analysis we present here not only …
#  10 Background         1        1 Yersinia pestis is the etiological agent of plague, alternatively …
#  # … with 184 more rows
count(txt, section)
#  # A tibble: 21 x 2
#     section                                                  n
#     <chr>                                                <int>
#   1 Abstract                                                 8
#   2 Authors' contributions                                   6
#   3 Background                                              20
#   4 Conclusion                                               3
#   5 Methods; Clustering analysis                             7
#   6 Methods; Collection of microarray expression data       17
#   7 Methods; Discovery of regulatory DNA motifs              8
#   8 Methods; Gel mobility shift analysis of Fur binding     13
#   9 Methods; Operon prediction                               5
#  10 Methods; Verification of predicted operons by RT-PCR     7
#  # … with 11 more rows
```

`pmc_caption` splits figure, table and supplementary material captions
into sentences.

``` r
cap1 <- pmc_caption(doc)
#  Found 5 figures
#  Found 4 tables
#  Found 3 supplements
filter(cap1, sentence == 1)
#  # A tibble: 12 x 4
#     tag      label               sentence text                                                       
#     <chr>    <chr>                  <int> <chr>                                                      
#   1 figure   Figure 1                   1 Environmental modulation of expression of virulence genes. 
#   2 figure   Figure 2                   1 RT-PCR analysis of potential operons.                      
#   3 figure   Figure 3                   1 Schematic representation of the clustered microarray data. 
#   4 figure   Figure 4                   1 Graphical representation of the consensus patterns by moti…
#   5 figure   Figure 5                   1 EMSA analysis of the binding of Fur protein to promoter DN…
#   6 table    Table 1                    1 Stress-responsive operons in Y. pestis predicted from micr…
#   7 table    Table 2                    1 Classification of the gene members of the cluster II in Fi…
#   8 table    Table 3                    1 Motif discovery for the clustering genes                   
#   9 table    Table 4                    1 Designs for expression profiling of Y. pestis              
#  10 supplem… Additional file 1 …        1 Growth curves of Y. pestis strain 201 under different cond…
#  11 supplem… Additional file 2 …        1 All the transcriptional changes of 4005 genes of Y. pestis…
#  12 supplem… Additional file 3 …        1 List of oligonucleotide primers used in this study.
```

`pmc_table` formats tables by collapsing multiline headers, expanding
rowspan and colspan attributes and adding subheadings into a new column.

``` r
tab1 <- pmc_table(doc)
#  Parsing 4 tables
#  Adding footnotes to Table 1
sapply(tab1, nrow)
#  Table 1 Table 2 Table 3 Table 4 
#       39      23       4      34
tab1[[1]]
#  # A tibble: 39 x 5
#     subheading           `Potential operon (r … `Gene ID`  `Putative or predicted fu… `Reference (s)`
#     <chr>                <chr>                  <chr>      <chr>                      <chr>          
#   1 Iron uptake or heme… yfeABCD operon* (r > … YPO2439-2… Transport/binding chelate… yfeABCD [54]   
#   2 Iron uptake or heme… hmuRSTUV operon (r > … YPO0279-0… Transport/binding hemin    hmuRSTUV [55]  
#   3 Iron uptake or heme… ysuJIHG* (r > 0.95)    YPO1529-1… Iron uptake                -              
#   4 Iron uptake or heme… sufABCDS* (r > 0.90)   YPO2400-2… Iron-regulated Fe-S clust… -              
#   5 Iron uptake or heme… YPO1854-1856* (r > 0.… YPO1854-1… Iron uptake or heme synth… -              
#   6 Sulfur metabolism    tauABCD operon (r > 0… YPO0182-0… Transport/binding taurine  tauABCD [56]   
#   7 Sulfur metabolism    ssuEADCB operon (r > … YPO3623-3… Sulphur metabolism         ssu operon [57]
#   8 Sulfur metabolism    cys operon (r > 0.92)  YPO3010-3… Cysteine synthesis         -              
#   9 Sulfur metabolism    YPO1317-1319 (r > 0.9… YPO1317-1… Sulfur metabolism?         -              
#  10 Sulfur metabolism    YPO4109-4111 (r > 0.9… YPO4109-4… Sulfur metabolism?         -              
#  # … with 29 more rows
```

Captions and footnotes are added as attributes.

``` r
attributes(tab1[[1]])
#  $names
#  [1] "subheading"                     "Potential operon (r value)"    
#  [3] "Gene ID"                        "Putative or predicted function"
#  [5] "Reference (s)"                 
#  
#  $row.names
#   [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
#  [33] 33 34 35 36 37 38 39
#  
#  $class
#  [1] "tbl_df"     "tbl"        "data.frame"
#  
#  $caption
#  [1] "Stress-responsive operons in Y. pestis predicted from microarray expression data"
#  
#  $footnotes
#  [1] "'r' represents the correlation coefficient of adjacent genes; '*' represent the defined operon has the similar expression pattern in two other published microarray datasets [7, 21]; '?' inferred functions of uncharacterized genes; '-' means the corresponding operons have not been experimentally validated in other bacteria."
```

Use `collapse_rows` to join column names and cell values in a semi-colon
delimited string (and then search using functions in the next section).

``` r
collapse_rows(tab1, na.string="-")
#  # A tibble: 100 x 3
#     table     row text                                                                               
#     <chr>   <int> <chr>                                                                              
#   1 Table 1     1 subheading=Iron uptake or heme synthesis; Potential operon (r value)=yfeABCD opero…
#   2 Table 1     2 subheading=Iron uptake or heme synthesis; Potential operon (r value)=hmuRSTUV oper…
#   3 Table 1     3 subheading=Iron uptake or heme synthesis; Potential operon (r value)=ysuJIHG* (r >…
#   4 Table 1     4 subheading=Iron uptake or heme synthesis; Potential operon (r value)=sufABCDS* (r …
#   5 Table 1     5 subheading=Iron uptake or heme synthesis; Potential operon (r value)=YPO1854-1856*…
#   6 Table 1     6 subheading=Sulfur metabolism; Potential operon (r value)=tauABCD operon (r > 0.90)…
#   7 Table 1     7 subheading=Sulfur metabolism; Potential operon (r value)=ssuEADCB operon (r > 0.97…
#   8 Table 1     8 subheading=Sulfur metabolism; Potential operon (r value)=cys operon (r > 0.92); Ge…
#   9 Table 1     9 subheading=Sulfur metabolism; Potential operon (r value)=YPO1317-1319 (r > 0.97); …
#  10 Table 1    10 subheading=Sulfur metabolism; Potential operon (r value)=YPO4109-4111 (r > 0.90); …
#  # … with 90 more rows
```

`pmc_reference` extracts the id, pmid, authors, year, title, journal,
volume, pages, and DOIs from reference tags.

``` r
ref1 <- pmc_reference(doc)
#  Found 76 citation tags
ref1
#  # A tibble: 76 x 9
#        id pmid   authors                  year title                 journal   volume pages doi      
#     <int> <chr>  <chr>                   <int> <chr>                 <chr>     <chr>  <chr> <chr>    
#   1     1 89938… Perry RD, Fetherston JD  1997 Yersinia pestis--eti… Clin Mic… 10     35-66 <NA>     
#   2     2 16053… Hinnebusch BJ            2005 The evolution of fle… Curr Iss… 7      197-… <NA>     
#   3     3 64693… Straley SC, Harmon PA    1984 Yersinia pestis grow… Infect I… 45     655-… <NA>     
#   4     4 15557… Huang XZ, Lindler LE     2004 The pH 6 antigen is … Infect I… 72     7212… 10.1128/…
#   5     5 15721… Pujol C, Bliska JB       2005 Turning Yersinia pat… Clin Imm… 114    216-… 10.1016/…
#   6     6 12732… Rhodius VA, LaRossa RA   2003 Uses and pitfalls of… Curr Opi… 6      114-… 10.1016/…
#   7     7 15342… Motin VL, Georgescu AM…  2004 Temporal global chan… J Bacter… 186    6298… 10.1128/…
#   8     8 15557… Han Y, Zhou D, Pang X,…  2004 Microarray analysis … Microbio… 48     791-… <NA>     
#   9     9 15777… Han Y, Zhou D, Pang X,…  2005 DNA microarray analy… Microbes… 7      335-… 10.1016/…
#  10    10 15808… Han Y, Zhou D, Pang X,…  2005 Comparative transcri… Res Micr… 156    403-… 10.1016/…
#  # … with 66 more rows
```

Finally, `pmc_metadata` saves journal and article metadata to a list.

``` r
pmc_metadata(doc)
#  $PMCID
#  [1] "PMC2231364"
#  
#  $Title
#  [1] "Comparative transcriptomics in Yersinia pestis: a global view of environmental modulation of gene expression"
#  
#  $Authors
#  [1] "Yanping Han, Jingfu Qiu, Zhaobiao Guo, He Gao, Yajun Song, Dongsheng Zhou, Ruifu Yang"
#  
#  $Year
#  [1] 2007
#  
#  $Journal
#  [1] "BMC Microbiology"
#  
#  $Volume
#  [1] "7"
#  
#  $Pages
#  [1] "96"
#  
#  $`Published online`
#  [1] "2007-10-29"
#  
#  $`Date received`
#  [1] "2007-6-2"
#  
#  $DOI
#  [1] "10.1186/1471-2180-7-96"
#  
#  $Publisher
#  [1] "BioMed Central"
```

## Searching text

There are a few functions to search within the `pmc_text` or collapsed
`pmc_table` output. `separate_text` uses the
[stringr](https://stringr.tidyverse.org/) package to extract any
matching regular expression.

``` r
separate_text(txt, "[ATCGN]{5,}")
#  # A tibble: 9 x 5
#    match       section                       paragraph sentence text                                 
#    <chr>       <chr>                             <int>    <int> <chr>                                
#  1 ACGCAATCGT… Results and Discussion; Comp…         2        3 A 16 basepair (bp) box (5'-ACGCAATCG…
#  2 AAACGTTTNC… Results and Discussion; Comp…         2        4 It is very similar to the E. coli Pu…
#  3 TGATAATGAT… Results and Discussion; Comp…         2        5 A 21 bp box (5'-TGATAATGATTATCATTATC…
#  4 GATAATGATA… Results and Discussion; Comp…         2        6 It is a 10-1-10 inverted repeat that…
#  5 TGANNNNNNT… Results and Discussion; Comp…         2        7 A 15 bp box (5'-TGANNNNNNTCAA-3') wa…
#  6 TTGATN      Results and Discussion; Comp…         2        8 It is a part of the E. coli Fnr box …
#  7 NATCAA      Results and Discussion; Comp…         2        8 It is a part of the E. coli Fnr box …
#  8 GTTAATTAA   Results and Discussion; Comp…         3        4 The ArcA regulator can recognize a r…
#  9 GTTAATTAAT… Results and Discussion; Comp…         3        5 An ArcA-box-like sequence (5'-GTTAAT…
```

A few wrappers search pre-defined patterns and add an extra step to
expand matched ranges. `separate_refs` matches references within
brackets using `\\[[0-9, -]+\\]` and expands ranges like `[7-9]`.

``` r
x <- separate_refs(txt)
x
#  # A tibble: 93 x 6
#        id match section   paragraph sentence text                                                    
#     <dbl> <chr> <chr>         <int>    <int> <chr>                                                   
#   1     1 [1]   Backgrou…         1        1 Yersinia pestis is the etiological agent of plague, alt…
#   2     2 [2]   Backgrou…         1        3 To produce a transmissible infection, Y. pestis coloniz…
#   3     3 [3]   Backgrou…         1        9 However, a few bacilli are taken up by tissue macrophag…
#   4     4 [4,5] Backgrou…         1       10 Residence in this niche also facilitates the bacteria's…
#   5     5 [4,5] Backgrou…         1       10 Residence in this niche also facilitates the bacteria's…
#   6     6 [6]   Backgrou…         2        1 A DNA microarray is able to determine simultaneous chan…
#   7     7 [7-9] Backgrou…         2        2 We and others have measured the gene expression profile…
#   8     8 [7-9] Backgrou…         2        2 We and others have measured the gene expression profile…
#   9     9 [7-9] Backgrou…         2        2 We and others have measured the gene expression profile…
#  10    10 [10]  Backgrou…         2        2 We and others have measured the gene expression profile…
#  # … with 83 more rows
filter(x, id == 8)
#  # A tibble: 5 x 6
#       id match    section                         paragraph sentence text                            
#    <dbl> <chr>    <chr>                               <int>    <int> <chr>                           
#  1     8 [7-9]    Background                              2        2 We and others have measured the…
#  2     8 [8-13,1… Background                              2        4 In order to acquire more regula…
#  3     8 [7-13,1… Results and Discussion                  2        1 Recently, many signature expres…
#  4     8 [7-9]    Results and Discussion; Virule…         3        1 As described previously, expres…
#  5     8 [8-10]   Methods; Collection of microar…         1        6 The genome-wide transcriptional…
```

`separate_genes` expands microbial gene operons like `hmsHFRS` into four
separate genes.

``` r
separate_genes(txt)
#  # A tibble: 103 x 6
#     gene  match  section                         paragraph sentence text                             
#     <chr> <chr>  <chr>                               <int>    <int> <chr>                            
#   1 purR  PurR   Abstract                                2        5 Several regulatory DNA motifs, p…
#   2 phoP  PhoP   Background                              2        3 We also identified the regulons …
#   3 ompR  OmpR   Background                              2        3 We also identified the regulons …
#   4 oxyR  OxyR   Background                              2        3 We also identified the regulons …
#   5 csrA  CsrA   Results and Discussion                  1        3 After the determination of the C…
#   6 slyA  SlyA   Results and Discussion                  1        3 After the determination of the C…
#   7 phoPQ PhoPQ  Results and Discussion                  1        3 After the determination of the C…
#   8 hmsH  hmsHF… Results and Discussion; Virule…         3        3 For example, the hemin storage l…
#   9 hmsF  hmsHF… Results and Discussion; Virule…         3        3 For example, the hemin storage l…
#  10 hmsR  hmsHF… Results and Discussion; Virule…         3        3 For example, the hemin storage l…
#  # … with 93 more rows
```

Finally, `separate_tags` expands locus tag ranges.

``` r
collapse_rows(tab1, na="-") %>%
  separate_tags("YPO")
#  # A tibble: 270 x 5
#     id      match      table    row text                                                             
#     <chr>   <chr>      <chr>  <int> <chr>                                                            
#   1 YPO2439 YPO2439-2… Table…     1 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   2 YPO2440 YPO2439-2… Table…     1 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   3 YPO2441 YPO2439-2… Table…     1 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   4 YPO2442 YPO2439-2… Table…     1 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   5 YPO0279 YPO0279-0… Table…     2 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   6 YPO0280 YPO0279-0… Table…     2 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   7 YPO0281 YPO0279-0… Table…     2 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   8 YPO0282 YPO0279-0… Table…     2 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#   9 YPO0283 YPO0279-0… Table…     2 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#  10 YPO1529 YPO1529-1… Table…     3 subheading=Iron uptake or heme synthesis; Potential operon (r va…
#  # … with 260 more rows
```

### Using `xml2`

The `pmc_*` functions use the [xml2](https://github.com/r-lib/xml2)
package for parsing and may fail in some situations, so it helps to know
how to parse `xml_documents`. Use `cat` and `as.character` to view nodes
returned by `xml_find_all`.

``` r
library(xml2)
refs <- xml_find_all(doc, "//ref")
refs[1]
#  {xml_nodeset (1)}
#  [1] <ref id="B1">\n  <citation citation-type="journal">\n    <person-group person-group-type="aut ...
cat(as.character(refs[1]))
#  <ref id="B1">
#    <citation citation-type="journal">
#      <person-group person-group-type="author">
#        <name>
#          <surname>Perry</surname>
#          <given-names>RD</given-names>
#        </name>
#        <name>
#          <surname>Fetherston</surname>
#          <given-names>JD</given-names>
#        </name>
#      </person-group>
#      <article-title>Yersinia pestis--etiologic agent of plague</article-title>
#      <source>Clin Microbiol Rev</source>
#      <year>1997</year>
#      <volume>10</volume>
#      <fpage>35</fpage>
#      <lpage>66</lpage>
#      <pub-id pub-id-type="pmid">8993858</pub-id>
#    </citation>
#  </ref>
```

Many journals use superscripts for references cited so they usually
appear after words like `results9` below.

``` r
# doc1 <- pmc_xml("PMC6385181")
doc1 <- read_xml(system.file("extdata/PMC6385181.xml", package = "tidypmc"))
gsub(".*\\. ", "", xml_text( xml_find_all(doc1, "//sec/p"))[2])
#  [1] "RNA-seq identifies the most relevant genes and RT-qPCR validates its results9, especially in the field of environmental and host adaptation10,11 and antimicrobial response12."
```

Find the tags using `xml_find_all` and then update the nodes by adding
brackets or other text.

``` r
bib <- xml_find_all(doc1, "//xref[@ref-type='bibr']")
bib[1]
#  {xml_nodeset (1)}
#  [1] <xref ref-type="bibr" rid="CR1">1</xref>
xml_text(bib) <- paste0(" [", xml_text(bib), "]")
bib[1]
#  {xml_nodeset (1)}
#  [1] <xref ref-type="bibr" rid="CR1"> [1]</xref>
```

The text is now separated from the reference. Note the `pmc_text`
function adds the brackets by default.

``` r
gsub(".*\\. ", "", xml_text( xml_find_all(doc1, "//sec/p"))[2])
#  [1] "RNA-seq identifies the most relevant genes and RT-qPCR validates its results [9], especially in the field of environmental and host adaptation [10], [11] and antimicrobial response [12]."
```

Genes, species and many other terms are often included within italic
tags. You can mark these nodes using the same code above or simply list
all the names in italics and search text or tables for matches, for
example three letter gene names in text below.

``` r
library(tibble)
x <- xml_name(xml_find_all(doc, "//*"))
tibble(tag=x) %>%
  count(tag, sort=TRUE)
#  # A tibble: 84 x 2
#     tag               n
#     <chr>         <int>
#   1 td              398
#   2 given-names     388
#   3 name            388
#   4 surname         388
#   5 italic          235
#   6 pub-id          129
#   7 tr              117
#   8 xref            108
#   9 year             80
#  10 article-title    77
#  # … with 74 more rows
it <- xml_text(xml_find_all(doc, "//sec//p//italic"), trim=TRUE)
it2 <- tibble(italic=it) %>%
  count(italic, sort=TRUE)
it2
#  # A tibble: 53 x 2
#     italic        n
#     <chr>     <int>
#   1 Y. pestis    46
#   2 in vitro      5
#   3 E. coli       4
#   4 psaEFABC      3
#   5 r             3
#   6 cis           2
#   7 fur           2
#   8 n             2
#   9 nrdHIEF       2
#  10 sufABCDSE     2
#  # … with 43 more rows
filter(it2, nchar(italic) == 3)
#  # A tibble: 8 x 2
#    italic     n
#    <chr>  <int>
#  1 cis        2
#  2 fur        2
#  3 cys        1
#  4 hmu        1
#  5 ybt        1
#  6 yfe        1
#  7 yfu        1
#  8 ymt        1
separate_text(txt, c("fur", "cys", "hmu", "ybt", "yfe", "yfu", "ymt"))
#  # A tibble: 9 x 5
#    match section                               paragraph sentence text                               
#    <chr> <chr>                                     <int>    <int> <chr>                              
#  1 ymt   Results and Discussion; Virulence ge…         3        4 The ymt gene encoding Yersinia mur…
#  2 fur   Results and Discussion; Clustering a…         3        2 It is noticeable that almost all o…
#  3 yfe   Results and Discussion; Clustering a…         3        4 Genes in category A (yfe, hmu, yfu…
#  4 hmu   Results and Discussion; Clustering a…         3        4 Genes in category A (yfe, hmu, yfu…
#  5 yfu   Results and Discussion; Clustering a…         3        4 Genes in category A (yfe, hmu, yfu…
#  6 ybt   Results and Discussion; Clustering a…         3        4 Genes in category A (yfe, hmu, yfu…
#  7 cys   Results and Discussion; Clustering a…         4        2 Genes responsible for sulfur uptak…
#  8 cys   Results and Discussion; Clustering a…         4        3 Cluster III contains members of th…
#  9 fur   Methods; Gel mobility shift analysis…         1        1 The entire coding region of the fu…
```
Parsing Europe PMC FTP files
================
Chris Stubben
August 6, 2019

The [Europe PMC FTP](https://europepmc.org/ftp/oa/) includes 2.5 million
open access articles separated into files with 10K articles each.
Download and unzip a recent series of PMC ids and load into R using the
`readr` package. A sample file with the first 10 articles is included in
the `tidypmc` package.

``` r
library(readr)
pmcfile <- system.file("extdata/PMC6358576_PMC6358589.xml", package = "tidypmc")
pmc <- read_lines(pmcfile)
```

Find the start of the article nodes.

``` r
a1 <- grep("^<article ", pmc)
head(a1)
#  [1]  2 30 38 52 62 69
n <- length(a1)
n
#  [1] 10
```

Read a single article by collapsing the lines into a new line separated
string.

``` r
library(xml2)
x1 <- paste(pmc[2:29], collapse="\n")
doc <- read_xml(x1)
doc
#  {xml_document}
#  <article article-type="case-report" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:mml="http://www.w3.org/1998/Math/MathML">
#  [1] <front>\n  <journal-meta>\n    <journal-id journal-id-type="nlm-ta">ACG Case Rep J</journal-i ...
#  [2] <body>\n  <sec sec-type="intro" id="sec1">\n    <title>Introduction</title>\n    <p>Bezoars a ...
#  [3] <back>\n  <ref-list>\n    <title>References</title>\n    <ref id="B1">\n      <label>1.</labe ...
```

Loop through the articles and save the metadata and text below. All 10K
articles takes about 10 minutes to run on a Mac laptop and returns 1.7M
sentences.

``` r
library(tidypmc)
a1 <- c(a1, length(pmc))
met1 <- vector("list", n)
txt1 <- vector("list", n)
for(i in seq_len(n)){
  doc <- read_xml(paste(pmc[a1[i]:(a1[i+1]-1)], collapse="\n"))
  m1 <- pmc_metadata(doc)
  id <- m1$PMCID
  message("Parsing ", i, ". ", id)
  met1[[i]] <- m1
  txt1[[i]] <- pmc_text(doc)
}
#  Parsing 1. PMC6358576
#  Parsing 2. PMC6358577
#  Parsing 3. PMC6358578
#  Parsing 4. PMC6358579
#  Parsing 5. PMC6358580
#  Parsing 6. PMC6358581
#  Parsing 7. PMC6358585
#  Note: removing table-wrap nested in sec/p tag
#  Note: removing fig nested in sec/p tag
#  Parsing 8. PMC6358587
#  Note: removing table-wrap nested in sec/p tag
#  Note: removing fig nested in sec/p tag
#  Parsing 9. PMC6358588
#  Note: removing fig nested in sec/p tag
#  Parsing 10. PMC6358589
#  Note: removing table-wrap nested in sec/p tag
#  Note: removing fig nested in sec/p tag
```

Combine the list of metadata and text into tables.

``` r
library(dplyr)
met <- bind_rows(met1)
names(txt1) <- met$PMCID
txt <- bind_rows(txt1, .id="PMCID")
met
#  # A tibble: 10 x 12
#     PMCID Title Authors  Year Journal Volume Pages `Published onli… `Date received` DOI   Publisher
#     <chr> <chr> <chr>   <int> <chr>   <chr>  <chr> <chr>            <chr>           <chr> <chr>    
#   1 PMC6… Endo… Dana B…  2018 ACG Ca… 5      e87   2018-12-5        2018-7-8        10.1… American…
#   2 PMC6… Chro… Scott …  2018 ACG Ca… 5      e94   2018-12-5        2018-5-5        10.1… American…
#   3 PMC6… Bile… Steffi…  2018 ACG Ca… 5      e88   2018-12-5        2018-5-7        10.1… American…
#   4 PMC6… New … Gordon…  2018 ACG Ca… 5      e92   2018-12-5        2018-3-3        10.1… American…
#   5 PMC6… Bile… Michae…  2018 ACG Ca… 5      e89   2018-12-5        2017-11-3       10.1… American…
#   6 PMC6… Fuso… Akshay…  2018 ACG Ca… 5      e99   2018-12-19       2018-3-8        10.1… American…
#   7 PMC6… Chor… Marcia…  2019 Genes … 20     56-68 2018-1-24        2017-9-1        10.1… Nature P…
#   8 PMC6… The … Tao Zh…  2019 Spinal… 57     141-… 2018-8-8         2017-12-19      10.1… Nature P…
#   9 PMC6… Natu… Marjol…  2019 Molecu… 20     115-… 2018-12-16       2018-10-22      10.1… Elsevier 
#  10 PMC6… Pred… Yury O…  2019 Molecu… 20     63-78 2018-11-16       2018-9-10       10.1… Elsevier 
#  # … with 1 more variable: Issue <chr>
txt
#  # A tibble: 1,083 x 5
#     PMCID    section    paragraph sentence text                                                      
#     <chr>    <chr>          <int>    <int> <chr>                                                     
#   1 PMC6358… Title              1        1 Endoscopic versus Surgical Intervention for Jejunal Bezoa…
#   2 PMC6358… Abstract           1        1 Bezoar-induced small bowel obstruction is a rare entity, …
#   3 PMC6358… Abstract           1        2 The cornerstone of treatment for intestinal bezoars has b…
#   4 PMC6358… Abstract           1        3 We present a patient with obstructive jejunal phytobezoar…
#   5 PMC6358… Introduct…         1        1 Bezoars are aggregates of undigested foreign material tha…
#   6 PMC6358… Introduct…         1        2 There are currently four classifications of bezoars: phyt…
#   7 PMC6358… Introduct…         1        3 Endoscopic treatment of bezoars causing intestinal obstru…
#   8 PMC6358… Case Repo…         1        1 A 60-year old diabetic woman with a past cholecystectomy …
#   9 PMC6358… Case Repo…         1        2 Physical examination revealed mild diffuse abdominal tend…
#  10 PMC6358… Case Repo…         1        3 Computed tomography (CT) of the abdomen and pelvis reveal…
#  # … with 1,073 more rows
```
---
output: github_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "# "
)
```

[![Build Status](https://travis-ci.org/ropensci/tidypmc.svg?branch=master)](https://travis-ci.org/ropensci/tidypmc)
[![Coverage status](https://codecov.io/gh/ropensci/tidypmc/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tidypmc?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tidypmc)](https://cran.r-project.org/package=tidypmc)
[![Downloads](https://cranlogs.r-pkg.org/badges/tidypmc)](https://CRAN.R-project.org/package=tidypmc)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/tidypmc?color=orange)](https://CRAN.R-project.org/package=tidypmc)

# tidypmc

The [Open Access subset] of [Pubmed Central] (PMC) includes 2.5 million articles
from biomedical and life sciences journals.  The full text XML files are freely
available for text mining from the [REST service] or [FTP site] but can be
challenging to parse. For example, section tags are nested to arbitrary depths,
formulas and tables may return incomprehensible text blobs and superscripted
references are pasted at the end of words.  The functions in the `tidypmc`
package are intended to return readable text and maintain the document
structure, so gene names and other terms can be associated with specific
sections, paragraphs, sentences or table rows.


## Installation

Use [remotes] to install the package.

```{r install, eval=FALSE}
remotes::install_github("ropensci/tidypmc")
```

## Load XML

Download a single XML document like [PMC2231364] from the [REST service] using
the `pmc_xml` function.

```{r pmc_xml, message=FALSE, echo=-1}
options(width=100)
library(tidypmc)
library(tidyverse)
doc <- pmc_xml("PMC2231364")
doc
```

The [europepmc] package includes additional functions to search PMC
and download full text.  Be sure to include the `OPEN_ACCESS` field in
the search since these are the only articles with full text XML available.

```{r epmc, echo=-1}
options(width=100)
library(europepmc)
yp <- epmc_search("title:(Yersinia pestis virulence) OPEN_ACCESS:Y")
select(yp, pmcid, pubYear, title) %>%
  print(n=5)
```


Save all `r nrow(yp)` results to a list of XML documents using the `epmc_ftxt` or `pmc_xml` function.

```{r purrr, eval=FALSE}
docs <- map(yp$pmcid, epmc_ftxt)
```


See the [PMC FTP vignette] for details on parsing the large XML files on the [FTP site]
with 10,000 articles each.


## Parse XML


The package includes five functions to parse the `xml_document`.


|R function     |Description                                                                |
|:--------------|:--------------------------------------------------------------------------|
|`pmc_text`     |Split section paragraphs into sentences with full path to subsection titles|
|`pmc_caption`  |Split figure, table and supplementary material captions into sentences     |
|`pmc_table`    |Convert table nodes into a list of tibbles                                 |
|`pmc_reference`|Format references cited into a tibble                                      |
|`pmc_metadata` |List journal and article metadata in front node                            |


The `pmc_text` function uses the [tokenizers] package to split section paragraphs into
sentences.  The function also removes any tables, figures or formulas that are nested
within paragraph tags, replaces superscripted references with brackets, adds carets and
underscores to other superscripts and subscripts and includes the full path to the
subsection title.

```{r pmc_text, echo=-1}
options(width=110)
txt <- pmc_text(doc)
txt
count(txt, section, sort=TRUE)
```


Load the [tidytext] package for further text processing.

```{r tidytext, echo=-1}
options(width=110)
library(tidytext)
x1 <- unnest_tokens(txt, word, text) %>%
  anti_join(stop_words) %>%
  filter(!word %in% 1:100)
filter(x1, str_detect(section, "^Results"))
filter(x1, str_detect(section, "^Results")) %>%
  count(word, sort = TRUE)
```



The `pmc_table` function formats tables by collapsing multiline headers,
expanding rowspan and colspan attributes and adding subheadings into a new column.

```{r pmc_table, echo=-1}
options(width=110)
tbls <- pmc_table(doc)
map_int(tbls, nrow)
tbls[[1]]
```

Use `collapse_rows` to join column names and cell values in a semi-colon delimited string (and
then search using functions in the next section).

```{r collapserows, echo=-1}
options(width=110)
collapse_rows(tbls, na.string="-")
```

The other three `pmc` functions are described in the package [vignette].


## Searching text

There are a few functions to search within the `pmc_text` or collapsed
`pmc_table` output.  `separate_text` uses the [stringr] package to extract any
regular expression or vector of words.


```{r separate_text, echo=-1}
options(width=110)
separate_text(txt, "[ATCGN]{5,}")
```

A few wrappers search pre-defined patterns and add an extra step to expand
matched ranges. `separate_refs` matches references within brackets using
`\\[[0-9, -]+\\]` and expands ranges like `[7-9]`.

```{r separate_refs, echo=-1}
options(width=110)
separate_refs(txt)
```

`separate_genes` will find microbial genes like tauD (with a
capitalized 4th letter)  and expand operons like `tauABCD` into
four genes.  `separate_tags` will find and expand locus tag ranges below.


```{r locus_tags, echo=-1}
options(width=110)
collapse_rows(tbls, na="-") %>%
  separate_tags("YPO") %>%
  filter(id == "YPO1855")
```


See the [vignette] for more details including code to parse
XML documents using the [xml2] package.  The [PMC FTP vignette]
has details on parsing XML files at the Europe PMC [FTP site].


### Community Guidelines

This project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms. Feedback, bug
reports, and feature requests are welcome
[here](https://github.com/ropensci/tidypmc/issues).


[remotes]: https://github.com/r-lib/remotes
[PMC2231364]: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2231364/fullTextXML
[Open Access subset]: https://europepmc.org/downloads/openaccess
[REST service]: https://europepmc.org/RestfulWebService
[FTP site]: https://europepmc.org/ftp/oa/
[tidytext]: https://www.tidytextmining.com/
[stringr]: https://stringr.tidyverse.org/
[vignette]: https://github.com/ropensci/tidypmc/blob/master/vignettes/tidypmc.md
[PMC FTP vignette]: https://github.com/ropensci/tidypmc/blob/master/vignettes/pmcftp.md
[tokenizers]: https://lincolnmullen.com/software/tokenizers/
[xml2]: https://github.com/r-lib/xml2
[europepmc]: https://github.com/ropensci/europepmc
[Pubmed Central]: https://europepmc.org
---
title: "Parsing Europe PMC FTP files"
author: "Chris Stubben"
date: '`r gsub("  ", " ", format(Sys.time(), "%B %e, %Y"))`'
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Parse PMC FTP files}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "# "
)
```


The [Europe PMC FTP] includes 2.5 million open access articles separated into
files with 10K articles each.  Download and unzip a recent series of PMC ids
and load into R using the `readr` package.   A sample file with the first 10
articles is included in the `tidypmc` package.

```{r load}
library(readr)
pmcfile <- system.file("extdata/PMC6358576_PMC6358589.xml", package = "tidypmc")
pmc <- read_lines(pmcfile)
```


Find the start of the article nodes.

```{r startnode}
a1 <- grep("^<article ", pmc)
head(a1)
n <- length(a1)
n
```

Read a single article by collapsing the lines into a new line separated string.


```{r read1, echo=-1}
options(width=100)
library(xml2)
x1 <- paste(pmc[2:29], collapse="\n")
doc <- read_xml(x1)
doc
```


Loop through the articles and save the metadata and text below.
All 10K articles takes about 10 minutes to run on a Mac laptop and returns 1.7M
sentences.


```{r loop}
library(tidypmc)
a1 <- c(a1, length(pmc))
met1 <- vector("list", n)
txt1 <- vector("list", n)
for(i in seq_len(n)){
  doc <- read_xml(paste(pmc[a1[i]:(a1[i+1]-1)], collapse="\n"))
  m1 <- pmc_metadata(doc)
  id <- m1$PMCID
  message("Parsing ", i, ". ", id)
  met1[[i]] <- m1
  txt1[[i]] <- pmc_text(doc)
}
```


Combine the list of metadata and text into tables.


```{r combine, echo=-1, message=FALSE}
options(width=100)
library(dplyr)
met <- bind_rows(met1)
names(txt1) <- met$PMCID
txt <- bind_rows(txt1, .id="PMCID")
met
txt
```




[Europe PMC FTP]: https://europepmc.org/ftp/oa/
---
title: "Introduction to tidypmc"
author: "Chris Stubben"
date: '`r gsub("  ", " ", format(Sys.time(), "%B %e, %Y"))`'
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Introduction to tidypmc}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "# "
)
```

The `tidypmc` package parses XML documents in the Open Access subset of [Pubmed Central].
Download the full text using `pmc_xml`.

```{r epmc_ftxt}
library(tidypmc)
doc <- pmc_xml("PMC2231364")
doc
```

The package includes five functions to parse the `xml_document`.


|R function     |Description                                                                |
|:--------------|:--------------------------------------------------------------------------|
|`pmc_text`     |Split section paragraphs into sentences with full path to subsection titles|
|`pmc_caption`  |Split figure, table and supplementary material captions into sentences     |
|`pmc_table`    |Convert table nodes into a list of tibbles                                 |
|`pmc_reference`|Format references cited into a tibble                                      |
|`pmc_metadata` |List journal and article metadata in front node                            |



`pmc_text` splits paragraphs into sentences and  removes any tables, figures or
formulas that are nested within paragraph tags, replaces superscripted
references with brackets, adds carets and underscores to other superscripts and
subscripts and includes the full path to the subsection title.

```{r pmc_text, message=FALSE, echo=-1}
options(width=100)
library(dplyr)
txt <- pmc_text(doc)
txt
count(txt, section)
```

`pmc_caption` splits figure, table and supplementary material captions into sentences.


```{r pmc_caption, echo=-1}
options(width=100)
cap1 <- pmc_caption(doc)
filter(cap1, sentence == 1)
```

`pmc_table` formats tables by collapsing multiline headers, expanding rowspan and
colspan attributes and adding subheadings into a new column.

```{r pmc_table, echo=-1}
options(width=100)
tab1 <- pmc_table(doc)
sapply(tab1, nrow)
tab1[[1]]
```

Captions and footnotes are added as attributes.

```{r attributes}
attributes(tab1[[1]])
```


Use `collapse_rows` to join column names and cell values in a semi-colon delimited string (and
then search using functions in the next section).

```{r collapserows, echo=-1}
options(width=100)
collapse_rows(tab1, na.string="-")
```


`pmc_reference` extracts the id, pmid, authors, year, title, journal, volume, pages,
and DOIs from reference tags.


```{r pmc_ref, echo=-1}
options(width=100)
ref1 <- pmc_reference(doc)
ref1
```


Finally, `pmc_metadata` saves journal and article metadata to a list.

```{r pmc_metadata}
pmc_metadata(doc)
```


## Searching text

There are a few functions to search within the `pmc_text` or collapsed `pmc_table` output.
`separate_text` uses the [stringr]  package to extract any matching regular expression.


```{r separate_text, echo=-1}
options(width=100)
separate_text(txt, "[ATCGN]{5,}")
```

A few wrappers search pre-defined patterns and add an extra step to expand matched ranges. `separate_refs`
matches references within brackets using `\\[[0-9, -]+\\]` and expands ranges like `[7-9]`.

```{r separate_refs, echo=-1}
options(width=100)
x <- separate_refs(txt)
x
filter(x, id == 8)
```

`separate_genes` expands microbial gene operons like `hmsHFRS` into four separate genes.

```{r separate_genes, echo=-1}
options(width=100)
separate_genes(txt)
```

Finally, `separate_tags` expands locus tag ranges.


```{r locus_tags, echo=-1}
options(width=100)
collapse_rows(tab1, na="-") %>%
  separate_tags("YPO")
```


### Using `xml2`

The `pmc_*` functions use the [xml2] package for parsing and may fail in some situations, so
it helps to know how to parse `xml_documents`.  Use `cat` and `as.character` to view nodes
returned by `xml_find_all`.

```{r catchar}
library(xml2)
refs <- xml_find_all(doc, "//ref")
refs[1]
cat(as.character(refs[1]))
```


Many journals use superscripts for references cited so they usually
appear after words like `results9` below.

```{r pmcdoc1, message=FALSE}
# doc1 <- pmc_xml("PMC6385181")
doc1 <- read_xml(system.file("extdata/PMC6385181.xml", package = "tidypmc"))
gsub(".*\\. ", "", xml_text( xml_find_all(doc1, "//sec/p"))[2])
```

Find the tags using `xml_find_all` and then update the nodes by adding brackets
or other text.

```{r bib}
bib <- xml_find_all(doc1, "//xref[@ref-type='bibr']")
bib[1]
xml_text(bib) <- paste0(" [", xml_text(bib), "]")
bib[1]
```

The text is now separated from the reference.  Note the `pmc_text` function adds the brackets by default.

```{r pmc_text2, message=FALSE}
gsub(".*\\. ", "", xml_text( xml_find_all(doc1, "//sec/p"))[2])
```


Genes, species and many other terms are often included within italic tags.  You
can mark these nodes using the same code above or simply list all the names
in italics and search text or tables for matches, for example three letter gene
names in text below.


```{r italicgenes}
library(tibble)
x <- xml_name(xml_find_all(doc, "//*"))
tibble(tag=x) %>%
  count(tag, sort=TRUE)
it <- xml_text(xml_find_all(doc, "//sec//p//italic"), trim=TRUE)
it2 <- tibble(italic=it) %>%
  count(italic, sort=TRUE)
it2
filter(it2, nchar(italic) == 3)
separate_text(txt, c("fur", "cys", "hmu", "ybt", "yfe", "yfu", "ymt"))
```




[stringr]: https://stringr.tidyverse.org/
[xml2]: https://github.com/r-lib/xml2
[europepmc]: https://github.com/ropensci/europepmc
[Pubmed Central]: https://europepmc.org
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmc_table.R
\name{pmc_table}
\alias{pmc_table}
\title{Convert table nodes to tibbles}
\usage{
pmc_table(doc)
}
\arguments{
\item{doc}{\code{xml_document} from PubMed Central}
}
\value{
a list of tibbles
}
\description{
Convert PubMed Central table nodes into a list of tibbles
}
\note{
Saves the caption and footnotes as attributes and collapses multiline
headers, expands all rowspan and colspan attributes and adds
subheadings to column one.
}
\examples{
# doc <- pmc_xml("PMC2231364")
doc <- xml2::read_xml(system.file("extdata/PMC2231364.xml",
  package = "tidypmc"
))
x <- pmc_table(doc)
sapply(x, dim)
x
attributes(x[[1]])
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapse_rows.R
\name{collapse_rows}
\alias{collapse_rows}
\title{Collapse a list of PubMed Central tables}
\usage{
collapse_rows(pmc, na.string)
}
\arguments{
\item{pmc}{a list of tables, usually from \code{\link{pmc_table}}}

\item{na.string}{additional cell values to skip, default is NA and ""}
}
\value{
A tibble with table and row number and collapsed text
}
\description{
Collapse rows into a semi-colon delimited list with column names and cell
values
}
\examples{
x <- data.frame(
  genes = c("aroB", "glnP", "ndhA", "pyrF"),
  fold_change = c(2.5, 1.7, -3.1, -2.6)
)
collapse_rows(list(`Table 1` = x))
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmc_caption.R
\name{pmc_caption}
\alias{pmc_caption}
\title{Split captions into sentences}
\usage{
pmc_caption(doc)
}
\arguments{
\item{doc}{\code{xml_document} from PubMed Central}
}
\value{
a tibble with tag, label, sentence number and text
}
\description{
Split figure, table and supplementary material captions into sentences
}
\examples{
# doc <- pmc_xml("PMC2231364") # OR
doc <- xml2::read_xml(system.file("extdata/PMC2231364.xml",
  package = "tidypmc"
))
x <- pmc_caption(doc)
x
dplyr::filter(x, sentence == 1)
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidypmc-package.R
\docType{package}
\name{tidypmc}
\alias{tidypmc}
\alias{tidypmc-package}
\title{\code{tidypmc} package}
\description{
Parse full text XML documents from PubMed Central
}
\details{
See the Github page for details at \url{https://github.com/ropensci/tidypmc}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmc_reference.R
\name{pmc_reference}
\alias{pmc_reference}
\title{Format references cited}
\usage{
pmc_reference(doc)
}
\arguments{
\item{doc}{\code{xml_document} from PubMed Central}
}
\value{
a tibble with id, pmid, authors, year, title, journal, volume, pages,
and doi.
}
\description{
Format references cited
}
\note{
Mixed citations without any child tags are added to the author column.
}
\examples{
# doc <- pmc_xml("PMC2231364")
doc <- xml2::read_xml(system.file("extdata/PMC2231364.xml",
  package = "tidypmc"
))
x <- pmc_reference(doc)
x
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmc_xml.R
\name{pmc_xml}
\alias{pmc_xml}
\title{Download XML from PubMed Central}
\source{
\url{https://europepmc.org/RestfulWebService}
}
\usage{
pmc_xml(id)
}
\arguments{
\item{id}{a PMC id starting with 'PMC'}
}
\value{
\code{xml_document}
}
\description{
Download XML from PubMed Central
}
\examples{
\dontrun{
doc <- pmc_xml("PMC2231364")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmc_metadata.R
\name{pmc_metadata}
\alias{pmc_metadata}
\title{Get article metadata}
\usage{
pmc_metadata(doc)
}
\arguments{
\item{doc}{\code{xml_document} from PubMed Central}
}
\value{
a list
}
\description{
Get a list of journal and article metadata in /front tag
}
\examples{
# doc <- pmc_xml("PMC2231364") # OR
doc <- xml2::read_xml(system.file("extdata/PMC2231364.xml",
  package = "tidypmc"
))
pmc_metadata(doc)
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separate_text.R
\name{separate_text}
\alias{separate_text}
\title{Separate all matching text into multiple rows}
\usage{
separate_text(txt, pattern, column = "text")
}
\arguments{
\item{txt}{a tibble, usually results from \code{pmc_text}}

\item{pattern}{either a regular expression or a vector of words to find in
text}

\item{column}{column name, default "text"}
}
\value{
a tibble
}
\description{
Separate all matching text into multiple rows
}
\note{
passed to \code{grepl} and \code{str_extract_all}
}
\examples{
# doc <- pmc_xml("PMC2231364")
doc <- xml2::read_xml(system.file("extdata/PMC2231364.xml",
        package = "tidypmc"))
txt <- pmc_text(doc)
separate_text(txt, "[ATCGN]{5,}")
separate_text(txt, "\\\\([A-Z]{3,6}s?\\\\)")
# pattern can be a vector of words
separate_text(txt, c("hmu", "ybt", "yfe", "yfu"))
# wrappers for separate_text with extra step to expand matched ranges
separate_refs(txt)
separate_genes(txt)
separate_tags(txt, "YPO")

}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmc_text.R
\name{pmc_text}
\alias{pmc_text}
\title{Split section paragraphs into sentences}
\usage{
pmc_text(doc)
}
\arguments{
\item{doc}{\code{xml_document} from PubMed Central}
}
\value{
a tibble with section, paragraph and sentence number and text
}
\description{
Split section paragraph tags into a table with subsection titles and
sentences using \code{tokenize_sentences}
}
\note{
Subsections may be nested to arbitrary depths and this function will
return the entire path to the subsection title as a delimited string like
"Results; Predicted functions; Pathogenicity".  Tables, figures and
formulas that are nested in section paragraphs are removed, superscripted
references are replaced with brackets, and any other superscripts or
subscripts are separared with ^ and _.
}
\examples{
# doc <- pmc_xml("PMC2231364")
doc <- xml2::read_xml(system.file("extdata/PMC2231364.xml",
  package = "tidypmc"
))
txt <- pmc_text(doc)
txt
dplyr::count(txt, section, sort = TRUE)
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separate_tags.R
\name{separate_tags}
\alias{separate_tags}
\title{Separate locus tag into multiple rows}
\usage{
separate_tags(txt, pattern, column = "text")
}
\arguments{
\item{txt}{a table}

\item{pattern}{regular expression to match locus tags like YPO[0-9-]+ or
the locus tag prefix like YPO.}

\item{column}{column name to search, default "text"}
}
\value{
a tibble with locus tag, matching text and rows.
}
\description{
Separates locus tags mentioned in full text and expands ranges like
YPO1970-74 into new rows
}
\examples{
x <- data.frame(row = 1, text = "some genes like YPO1002 and YPO1970-74")
separate_tags(x, "YPO")
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separate_refs.R
\name{separate_refs}
\alias{separate_refs}
\title{Separate references cited into multiple rows}
\usage{
separate_refs(txt, column = "text")
}
\arguments{
\item{txt}{a table}

\item{column}{column name, default "text"}
}
\value{
a tibble
}
\description{
Separates references cited in brackets or parentheses into multiple rows and
splits the comma-delimited numeric strings and expands ranges like 7-9 into
new rows
}
\examples{
x <- data.frame(row = 1, text = "some important studies [7-9,15]")
separate_refs(x)
}
\author{
Chris Stubben
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separate_genes.R
\name{separate_genes}
\alias{separate_genes}
\title{Separate genes and operons into multiple rows}
\usage{
separate_genes(txt, pattern = "\\\\b[A-Za-z][a-z]{2}[A-Z0-9]+\\\\b",
  genes, operon = 6, column = "text")
}
\arguments{
\item{txt}{a table}

\item{pattern}{regular expression to match genes, default is to match
microbial genes like AbcD, default [A-Za-z][a-z]{2}[A-Z0-9]+}

\item{genes}{an optional vector of genes, set pattern to NA to only match
this list.}

\item{operon}{operon length, default 6. Split genes with 6 or more letters
into separate genes, for example AbcDEF is split into abcD, abcE and abcF.}

\item{column}{column name to search, default "text"}
}
\value{
a tibble with gene name, matching text and rows.
}
\description{
Separate genes and operons mentioned in full text into multiple rows
}
\note{
Check for genes in italics using \code{xml_text(xml_find_all(doc,
"//sec//p//italic"))} and update the pattern or add additional genes as an
optional vector if needed
}
\examples{
x <- data.frame(row = 1, text = "Genes like YacK, hmu and sufABC")
separate_genes(x)
separate_genes(x, genes = "hmu")
}
\author{
Chris Stubben
}
