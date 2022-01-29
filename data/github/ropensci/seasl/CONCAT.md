seasl
=====



[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
![R build status](https://github.com/ropensci/seasl/workflows/R-CMD-check/badge.svg)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/seasl?color=C9A115)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/seasl)](https://cran.r-project.org/package=seasl)
<!-- [![Build Status](https://travis-ci.org/ropensci/seasl.svg?branch=master)](https://travis-ci.org/ropensci/seasl) -->

`seasl` is an R client for exploring CSL styles.

This package is inspired by the Ruby gem `csl`: https://github.com/inukshuk/csl-ruby

The Citation Style Language 1.0.1 specification: http://docs.citationstyles.org/en/1.0.1/specification.html

Package API:

 - `csl_locales`
 - `as.location`
 - `csl_styles`
 - `csl_locale_exists`
 - `csl_cache`
 - `csl_fetch_styles`
 - `csl_style_find`
 - `csl_style_xml`
 - `csl_style_load`
 - `csl_style_exists`
 - `csl_locale_load`
 - `csl_fetch_locales`

## Install


```r
install.packages("seasl")
```

or


```r
remotes::install_github("ropensci/seasl")
```


```r
library("seasl")
```

## Download styles and locales

First, you may want to download style and locale files. `csl_fetch_styles()` and `csl_fetch_locales()`
download the files to your machine. See `?csl_cache` for caching information, including
how to change the cache location.

Styles retrieved from the Github repo at https://github.com/citation-style-language/styles-distribution


```r
csl_fetch_styles()
#>
#> Done! Files put in /Users/sckott/Library/Caches/R/seasl/styles
```

Locales retrieved from the Github repo at https://github.com/citation-style-language/locales


```r
csl_fetch_locales()
#>
#> Done! Files put in /Users/sckott/Library/Caches/R/seasl/locales
```

## File paths to CSL styles and locales

calling `csl_styles` without inputs gives all styles, with separate lists for 
dependent and independent styles


```r
csl_styles()
#> $independent
#>    [1] "academy-of-management-review"                                                                                        
#>    [2] "accident-analysis-and-prevention"                                                                                    
#>    [3] "aci-materials-journal"                                                                                               
#>    [4] "acm-sig-proceedings-long-author-list"                                                                                
#>    [5] "acm-sig-proceedings"                                                                                                 
#>    [6] "acm-sigchi-proceedings-extended-abstract-format"                                                                     
#>    [7] "acm-sigchi-proceedings"                                                                                              
#>    [8] "acm-siggraph"                                                                                                        
#>    [9] "acme-an-international-journal-for-critical-geographies"                                                              
...
```

calling `csl_styles` with an input gives the path to that style, if found


```r
csl_styles("apa")
#> [1] "/Users/sckott/Library/Caches/R/seasl/styles/apa.csl"
csl_styles("archeologie-medievale")
#> [1] "/Users/sckott/Library/Caches/R/seasl/styles/archeologie-medievale.csl"
```

Same patterns go for locales (note that there are far fewer locales than styles)


```r
# just locale names
csl_locales()
#>  [1] "locales-af-ZA" "locales-ar"    "locales-bg-BG" "locales-ca-AD"
#>  [5] "locales-cs-CZ" "locales-cy-GB" "locales-da-DK" "locales-de-AT"
#>  [9] "locales-de-CH" "locales-de-DE" "locales-el-GR" "locales-en-GB"
#> [13] "locales-en-US" "locales-es-CL" "locales-es-ES" "locales-es-MX"
#> [17] "locales-et-EE" "locales-eu"    "locales-fa-IR" "locales-fi-FI"
#> [21] "locales-fr-CA" "locales-fr-FR" "locales-he-IL" "locales-hr-HR"
#> [25] "locales-hu-HU" "locales-id-ID" "locales-is-IS" "locales-it-IT"
#> [29] "locales-ja-JP" "locales-km-KH" "locales-ko-KR" "locales-la"   
#> [33] "locales-lt-LT" "locales-lv-LV" "locales-mn-MN" "locales-nb-NO"
#> [37] "locales-nl-NL" "locales-nn-NO" "locales-pl-PL" "locales-pt-BR"
#> [41] "locales-pt-PT" "locales-ro-RO" "locales-ru-RU" "locales-sk-SK"
#> [45] "locales-sl-SI" "locales-sr-RS" "locales-sv-SE" "locales-th-TH"
#> [49] "locales-tr-TR" "locales-uk-UA" "locales-vi-VN" "locales-zh-CN"
#> [53] "locales-zh-TW"
```


```r
# when locale given, you get the full path
csl_locales("fr-FR")
#> [1] "/Users/sckott/Library/Caches/R/seasl/locales/locales-fr-FR.xml"
```

Alternatively, you can try to find a style by using `csl_style_find()`


```r
# single match
csl_style_find(x = "American Journal of Epidemiology")
#> [1] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-epidemiology.csl"
```


```r
# many matches
csl_style_find(x = "American Journal")
#>  [1] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-agricultural-economics.csl"                                     
#>  [2] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-archaeology.csl"                                                
#>  [3] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-botany.csl"                                                     
#>  [4] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-climate-change.csl"                                             
#>  [5] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-clinical-pathology.csl"                                         
#>  [6] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-enology-and-viticulture.csl"                                    
#>  [7] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-epidemiology.csl"                                               
#>  [8] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-health-behavior.csl"                                            
#>  [9] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-hypertension.csl"                                               
#> [10] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-medical-genetics.csl"                                           
...
```

## Load CSL style from a URL


```r
jps <- csl_style_load('http://www.zotero.org/styles/american-journal-of-political-science')
```

## Query style information


```r
jps$info
#> $title
#> [1] "American Journal of Political Science"
#> 
#> $`title-short`
#> [1] "AJPS"
#> 
#> $id
#> [1] "http://www.zotero.org/styles/american-journal-of-political-science"
#> 
#> $contributor
...
```


```r
jps$info$title
#> [1] "American Journal of Political Science"
```


```r
jps$macros
#> [[1]]
#> [[1]]$name
#> [1] "editor"
#> 
#> [[1]][[2]]
#> [[1]][[2]]$names
#> [[1]][[2]]$names$variable
#> [1] "editor"
#> 
#> [[1]][[2]]$names$delimiter
...
```


```r
jps$citation
#> $sort
#> $sort$key
#> list()
#> attr(,"macro")
#> [1] "author-short"
#> 
#> $sort$key
#> list()
#> attr(,"macro")
#> [1] "year-date"
...
```


```r
jps$bibliography
#> $attributes
#> $attributes$`hanging-indent`
#> [1] "true"
#> 
#> $attributes$`et-al-min`
#> [1] "4"
#> 
#> $attributes$`et-al-use-first`
#> [1] "1"
#> 
...
```

## Get raw XML


```r
csl_style_xml('http://www.zotero.org/styles/american-journal-of-political-science')
#> {xml_document}
#> <style class="in-text" version="1.0" demote-non-dropping-particle="sort-only" default-locale="en-US" xmlns="http://purl.org/net/xbiblio/csl">
#>  [1] <info>\n  <title>American Journal of Political Science</title>\n  <title ...
#>  [2] <macro name="editor">\n  <names variable="editor" delimiter=", ">\n    < ...
#>  [3] <macro name="author">\n  <names variable="author">\n    <name name-as-so ...
#>  [4] <macro name="author-short">\n  <names variable="author">\n    <name form ...
#>  [5] <macro name="access">\n  <choose>\n    <if type="legal_case" match="none ...
#>  [6] <macro name="title">\n  <choose>\n    <if type="bill book graphic legal_ ...
#>  [7] <macro name="legal_case">\n  <group prefix=" " delimiter=" ">\n    <text ...
#>  [8] <macro name="publisher">\n  <choose>\n    <if type="thesis" match="none" ...
...
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/seasl/issues).
* License: MIT
* Citation: execute `citation(package = 'seasl')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/seasl/blob/master/CODE_OF_CONDUCT.md
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

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

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that seasl is released with a [Contributor Code of Conduct][coc]. By contributing to this project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html) for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 

[coc]: https://github.com/ropensci/seasl/blob/master/CODE_OF_CONDUCT.md
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
