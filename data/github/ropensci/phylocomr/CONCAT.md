phylocomr
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/phylocomr/actions/workflows/R-check.yml/badge.svg)](https://github.com/ropensci/phylocomr/actions/workflows/R-check.yml)
[![cran checks](https://cranchecks.info/badges/worst/phylocomr)](https://cranchecks.info/pkgs/phylocomr)
[![codecov](https://codecov.io/gh/ropensci/phylocomr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/phylocomr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/phylocomr)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/phylocomr)](https://cran.r-project.org/package=phylocomr)

`phylocomr` gives you access to the [Phylocom C library](https://github.com/phylocom/phylocom/), licensed under [BSD 2-clause](http://www.opensource.org/licenses/bsd-license.php)

## Package API

* `ecovolve`/`ph_ecovolve` - interface to `ecovolve` executable, and a higher
level interface
* `phylomatic`/`ph_phylomatic` - interface to `phylomatic` executable, and a higher
level interface
* `phylocom` - interface to `phylocom` executable
* `ph_aot` - higher level interface to `aot`
* `ph_bladj` - higher level interface to `bladj`
* `ph_comdist`/`ph_comdistnt` - higher level interface to comdist
* `ph_comstruct` - higher level interface to comstruct
* `ph_comtrait` - higher level interface to comtrait
* `ph_pd` - higher level interface to Faith's phylogenetic diversity

## A note about files

As a convienence you can pass ages, sample and trait data.frame's, and phylogenies as strings, to `phylocomr` functions. However, `phylocomr` has to write these data.frame's/strings to disk (your computer's file system) to be able to run the Phylocom code on them. Internally, `phylocomr` is writing to a temporary file to run Phylocom code, and then the file is removed.

In addition, you can pass in files instead of data.frame's/strings. These are not themselves used. Instead, we read and write those files to temporary files. We do this for two reasons. First, Phylocom expects the files its using to be in the same directory, so if we control the file paths that becomes easier. Second, Phylocom is case sensitive, so we simply standardize all taxon names by lower casing all of them. We do this case manipulation on the temporary files so that your original data files are not modified.

## Installation

Stable version:


```r
install.packages("phylocomr")
```

Development version:


```r
remotes::install_github("ropensci/phylocomr")
```


```r
library("phylocomr")
library("ape")
```

## ecovolve


```r
ph_ecovolve(speciation = 0.05, extinction = 0.005, time_units = 50)
```

## phylomatic


```r
taxa_file <- system.file("examples/taxa", package = "phylocomr")
phylo_file <- system.file("examples/phylo", package = "phylocomr")
(taxa_str <- readLines(taxa_file))
```

```
#> [1] "campanulaceae/lobelia/lobelia_conferta"          
#> [2] "cyperaceae/mapania/mapania_africana"             
#> [3] "amaryllidaceae/narcissus/narcissus_cuatrecasasii"
```

```r
(phylo_str <- readLines(phylo_file))
```

```
#> [1] "(((((eliea_articulata,homalanthus_populneus)malpighiales,rosa_willmottiae),((macrocentrum_neblinae,qualea_clavata),hibiscus_pohlii)malvids),(((lobelia_conferta,((millotia_depauperata,(layia_chrysanthemoides,layia_pentachaeta)layia),senecio_flanaganii)asteraceae)asterales,schwenkia_americana),tapinanthus_buntingii)),(narcissus_cuatrecasasii,mapania_africana))poales_to_asterales;"
```

```r
ph_phylomatic(taxa = taxa_str, phylo = phylo_str)
```

```
#> [1] "(lobelia_conferta:5.000000,(mapania_africana:1.000000,narcissus_cuatrecasasii:1.000000):1.000000)poales_to_asterales:1.000000;\n"
#> attr(,"taxa_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/taxa_2c8c72514e28"
#> attr(,"phylo_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/phylo_2c8c37e25eba"
```

use various different references trees


```r
library(brranching)
library(ape)

r2 <- ape::read.tree(text=brranching::phylomatic_trees[['R20120829']])
smith2011 <- ape::read.tree(text=brranching::phylomatic_trees[['smith2011']])
zanne2014 <- ape::read.tree(text=brranching::phylomatic_trees[['zanne2014']])

# R20120829 tree
taxa_str <- c(
  "asteraceae/bidens/bidens_alba",
  "asteraceae/cirsium/cirsium_arvense",
  "fabaceae/lupinus/lupinus_albus"
)
ph_phylomatic(taxa = taxa_str, phylo = r2)
```

```
#> [1] "(((bidens_alba:13.000000,cirsium_arvense:13.000000):19.000000,lupinus_albus:27.000000):12.000000)euphyllophyte:1.000000;\n"
#> attr(,"taxa_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/taxa_2c8c3c672523"
#> attr(,"phylo_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/phylo_2c8c17b73bd0"
```

```r
# zanne2014 tree
taxa_str <- c(
  "zamiaceae/dioon/dioon_edule",
  "zamiaceae/encephalartos/encephalartos_dyerianus",
  "piperaceae/piper/piper_arboricola"
)
ph_phylomatic(taxa = taxa_str, phylo = zanne2014)
```

```
#> [1] "(((dioon_edule:121.744843,encephalartos_dyerianus:121.744850)zamiaceae:230.489838,piper_arboricola:352.234711)spermatophyta:88.058670):0.000000;\n"
#> attr(,"taxa_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/taxa_2c8c7eb7e0d9"
#> attr(,"phylo_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/phylo_2c8c5a0b078d"
```

```r
# zanne2014 subtree
zanne2014_subtr <- ape::extract.clade(zanne2014, node='Loganiaceae')
zanne_subtree_file <- tempfile(fileext = ".txt")
ape::write.tree(zanne2014_subtr, file = zanne_subtree_file)
taxa_str <- c(
  "loganiaceae/neuburgia/neuburgia_corynocarpum",
  "loganiaceae/geniostoma/geniostoma_borbonicum",
  "loganiaceae/strychnos/strychnos_darienensis"
)
ph_phylomatic(taxa = taxa_str, phylo = zanne2014_subtr)
```

```
#> [1] "((neuburgia_corynocarpum:32.807743,(geniostoma_borbonicum:32.036335,strychnos_darienensis:32.036335):0.771406):1.635496)loganiaceae:0.000000;\n"
#> attr(,"taxa_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/taxa_2c8c6fce2295"
#> attr(,"phylo_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/phylo_2c8c46409d8b"
```

```r
ph_phylomatic(taxa = taxa_str, phylo = zanne_subtree_file)
```

```
#> [1] "((neuburgia_corynocarpum:32.807743,(geniostoma_borbonicum:32.036335,strychnos_darienensis:32.036335):0.771406):1.635496)loganiaceae:0.000000;\n"
#> attr(,"taxa_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/taxa_2c8c3c2734b5"
#> attr(,"phylo_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/phylo_2c8c31f575ed"
```

## aot


```r
traits_file <- system.file("examples/traits_aot", package = "phylocomr")
phylo_file <- system.file("examples/phylo_aot", package = "phylocomr")
traitsdf_file <- system.file("examples/traits_aot_df", package = "phylocomr")
traits <- read.table(text = readLines(traitsdf_file), header = TRUE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(phylo_file)
ph_aot(traits = traits, phylo = phylo_str)
```

```
#> $trait_conservatism
#> # A tibble: 124 x 28
#>    trait trait.name  node name    age ntaxa n.nodes tip.mn tmn.ranklow
#>    <int> <chr>      <int> <chr> <dbl> <int>   <int>  <dbl>       <int>
#>  1     1 traitA         0 a         5    32       2   1.75        1000
#>  2     1 traitA         1 b         4    16       2   1.75         640
#>  3     1 traitA         2 c         3     8       2   1.75         679
#>  4     1 traitA         3 d         2     4       2   1.5          251
#>  5     1 traitA         4 e         1     2       2   1             59
#>  6     1 traitA         7 f         1     2       2   2           1000
#>  7     1 traitA        10 g         2     4       2   2           1000
#>  8     1 traitA        11 h         1     2       2   2           1000
#>  9     1 traitA        14 i         1     2       2   2           1000
#> 10     1 traitA        17 j         3     8       2   1.75         664
#> # … with 114 more rows, and 19 more variables: tmn.rankhi <int>, tip.sd <dbl>,
#> #   tsd.ranklow <int>, tsd.rankhi <int>, node.mn <dbl>, nmn.ranklow <int>,
#> #   nmn.rankhi <int>, nod.sd <dbl>, nsd.ranklow <int>, nsd.rankhi <int>,
#> #   sstipsroot <dbl>, sstips <dbl>, percvaramongnodes <dbl>,
#> #   percvaratnode <dbl>, contributionindex <dbl>, sstipvnoderoot <dbl>,
#> #   sstipvnode <dbl>, ssamongnodes <dbl>, sswithinnodes <dbl>
#> 
#> $independent_contrasts
#> # A tibble: 31 x 17
#>     node name    age n.nodes contrast1 contrast2 contrast3 contrast4 contrastsd
#>    <int> <chr> <dbl>   <int>     <dbl>     <dbl>     <dbl>     <dbl>      <dbl>
#>  1     0 a         5       2     0         0         0         0.254       1.97
#>  2     1 b         4       2     0         1.03      0         0.516       1.94
#>  3     2 c         3       2     0.267     0.535     0         0           1.87
#>  4     3 d         2       2     0.577     0         1.15      0           1.73
#>  5     4 e         1       2     0         0         0.707     0           1.41
#>  6     7 f         1       2     0         0         0.707     0           1.41
#>  7    10 g         2       2     0         0         1.15      0           1.73
#>  8    11 h         1       2     0         0         0.707     0           1.41
#>  9    14 i         1       2     0         0         0.707     0           1.41
#> 10    17 j         3       2     0.267     0.535     0         0           1.87
#> # … with 21 more rows, and 8 more variables: lowval1 <dbl>, hival1 <dbl>,
#> #   lowval2 <dbl>, hival2 <dbl>, lowval3 <dbl>, hival3 <dbl>, lowval4 <dbl>,
#> #   hival4 <dbl>
#> 
#> $phylogenetic_signal
#> # A tibble: 4 x 5
#>   trait  ntaxa varcontr varcn.ranklow varcn.rankhi
#>   <chr>  <int>    <dbl>         <int>        <int>
#> 1 traitA    32    0.054             1         1000
#> 2 traitB    32    0.109             1         1000
#> 3 traitC    32    0.622            74          927
#> 4 traitD    32    0.011             1         1000
#> 
#> $ind_contrast_corr
#> # A tibble: 3 x 6
#>   xtrait ytrait ntaxa  picr  npos ncont
#>   <chr>  <chr>  <int> <dbl> <dbl> <int>
#> 1 traitA traitB    32 0.248  18.5    31
#> 2 traitA traitC    32 0.485  27.5    31
#> 3 traitA traitD    32 0      16.5    31
```

## bladj


```r
ages_file <- system.file("examples/ages", package = "phylocomr")
phylo_file <- system.file("examples/phylo_bladj", package = "phylocomr")
ages_df <- data.frame(
  a = c('malpighiales','salicaceae','fabaceae','rosales','oleaceae',
        'gentianales','apocynaceae','rubiaceae'),
  b = c(81,20,56,76,47,71,18,56)
)
phylo_str <- readLines(phylo_file)
(res <- ph_bladj(ages = ages_df, phylo = phylo_str))
```

```
#> [1] "((((((lomatium_concinnum:20.250000,campanula_vandesii:20.250000):20.250000,(((veronica_candidissima:10.125000,penstemon_paniculatus:10.125000)plantaginaceae:10.125000,justicia_oblonga:20.250000):10.125000,marsdenia_gilgiana:30.375000):10.125000):10.125000,epacris_alba-compacta:50.625000)ericales_to_asterales:10.125000,((daphne_anhuiensis:20.250000,syzygium_cumini:20.250000)malvids:20.250000,ditaxis_clariana:40.500000):20.250000):10.125000,thalictrum_setulosum:70.875000)eudicots:10.125000,((dendrocalamus_giganteus:27.000000,guzmania_densiflora:27.000000)poales:27.000000,warczewiczella_digitata:54.000000):27.000000)malpighiales:1.000000;\n"
#> attr(,"ages_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/ages"
#> attr(,"phylo_file")
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpFwkBRm/phylo_2c8c76ed1b96"
```

```r
plot(ape::read.tree(text = res))
```

![plot of chunk unnamed-chunk-9](man/figures/unnamed-chunk-9-1.png)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/phylocomr/issues).
* License: MIT
* Get citation information for `phylocomr` in R doing `citation(package = 'phylocomr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
phylocomr 0.3.2
===============

### MINOR IMPROVEMENTS

* move readme image into man/figures (#30)

### BUG FIXES

* fix for gcc `-fno-common` (#29) (#31)

phylocomr 0.3.0
===============

### NEW FEATURES

* via (#26) (see below) - we now no longer use file paths passed in directly to functions, but instead write to temporary files to run with Phylocom so that we do not alter at all the users files. We note this in the README and package level manual file `?phylocomr`
* package gains new manual file `?phylocomr-inputs` that details the four types of inputs to functions and what format they are expected in, including how they differ for passing in data.frame's vs. file paths (#28)

### BUG FIXES

* for all data.frame traits inputs to fxns, check that the first column is called `name` (Phylocom doesn't accept anything else) (#27)
* fix was originally for `ph_aot()`, but realized this touches almost all functions: Phylocom is case sensitive. We were already making sure all taxon names in phylogenies (tips and nodes) were lowercased, and were lowercasing names in tables passed in, but were not fixing case in file paths passed in by the user. Now across all functions we make sure case is all lowercase for taxon names in any user inputs, so case problems should no longer be an issue. (#26) via @Jez-R

phylocomr 0.2.0
===============

### BUG FIXES

* two fixes for `ph_bladj()`: 1) now we lowercase the taxon name column in the ages data.frame before writing the data.frame to disk to avoid any mismatch due to case (we write the phylogeny to disk with lowercased names); 2) bladj expects the root node name from the phyologeny to be in the ages file; we now check that (#25)


phylocomr 0.1.4
===============

### BUG FIXES

* fix examples (#22)
* improve class checks in internal code, swap `inherits` for `class` (#23)
* small fix to use of fread in C lib; check that fread worked, and if not if it was an EOF error or other error (#24)


phylocomr 0.1.2
===============

### MINOR IMPROVEMENTS

* fixes for failed checks on Solaris

### BUG FIXES

* fix to internals of all functions that handle a phylogeny. `ph_phylomatic` was working fine with very simple trees in all lowercase. we now inernally lowercase all node and tip labels, on any phylogeny inputs (phylo object, newick string, file path (read, then re-write back to disk)). phylomatic wasn't working with any uppercase labels. 


phylocomr 0.1.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 3.6.2 patched
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

No errors were found with the one reverse dependency.

--------

This release includes a fix for installation failures with gcc trunk, and moves README image to man/figures.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/phylocomr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/phylocomr.git`
* Make sure to track progress upstream (i.e., on our version of `phylocomr` at `ropensci/phylocomr`) by doing `git remote add upstream https://github.com/ropensci/phylocomr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/phylocomr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Get in touch: [myrmecocystus@gmail.com](mailto:myrmecocystus@gmail.com)

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. In addition, add a reproducible example if possible. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
phylocomr
=========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.path = "man/figures/"
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/phylocomr/actions/workflows/R-check.yml/badge.svg)](https://github.com/ropensci/phylocomr/actions/workflows/R-check.yml)
[![cran checks](https://cranchecks.info/badges/worst/phylocomr)](https://cranchecks.info/pkgs/phylocomr)
[![codecov](https://codecov.io/gh/ropensci/phylocomr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/phylocomr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/phylocomr)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/phylocomr)](https://cran.r-project.org/package=phylocomr)

`phylocomr` gives you access to the [Phylocom C library](https://github.com/phylocom/phylocom/), licensed under [BSD 2-clause](http://www.opensource.org/licenses/bsd-license.php)

## Package API

* `ecovolve`/`ph_ecovolve` - interface to `ecovolve` executable, and a higher
level interface
* `phylomatic`/`ph_phylomatic` - interface to `phylomatic` executable, and a higher
level interface
* `phylocom` - interface to `phylocom` executable
* `ph_aot` - higher level interface to `aot`
* `ph_bladj` - higher level interface to `bladj`
* `ph_comdist`/`ph_comdistnt` - higher level interface to comdist
* `ph_comstruct` - higher level interface to comstruct
* `ph_comtrait` - higher level interface to comtrait
* `ph_pd` - higher level interface to Faith's phylogenetic diversity

## A note about files

As a convienence you can pass ages, sample and trait data.frame's, and phylogenies as strings, to `phylocomr` functions. However, `phylocomr` has to write these data.frame's/strings to disk (your computer's file system) to be able to run the Phylocom code on them. Internally, `phylocomr` is writing to a temporary file to run Phylocom code, and then the file is removed.

In addition, you can pass in files instead of data.frame's/strings. These are not themselves used. Instead, we read and write those files to temporary files. We do this for two reasons. First, Phylocom expects the files its using to be in the same directory, so if we control the file paths that becomes easier. Second, Phylocom is case sensitive, so we simply standardize all taxon names by lower casing all of them. We do this case manipulation on the temporary files so that your original data files are not modified.

## Installation

Stable version:

```{r eval=FALSE}
install.packages("phylocomr")
```

Development version:

```{r eval=FALSE}
remotes::install_github("ropensci/phylocomr")
```

```{r}
library("phylocomr")
library("ape")
```

## ecovolve

```{r eval=FALSE}
ph_ecovolve(speciation = 0.05, extinction = 0.005, time_units = 50)
```

## phylomatic

```{r}
taxa_file <- system.file("examples/taxa", package = "phylocomr")
phylo_file <- system.file("examples/phylo", package = "phylocomr")
(taxa_str <- readLines(taxa_file))
(phylo_str <- readLines(phylo_file))
ph_phylomatic(taxa = taxa_str, phylo = phylo_str)
```

use various different references trees

```{r}
library(brranching)
library(ape)

r2 <- ape::read.tree(text=brranching::phylomatic_trees[['R20120829']])
smith2011 <- ape::read.tree(text=brranching::phylomatic_trees[['smith2011']])
zanne2014 <- ape::read.tree(text=brranching::phylomatic_trees[['zanne2014']])

# R20120829 tree
taxa_str <- c(
  "asteraceae/bidens/bidens_alba",
  "asteraceae/cirsium/cirsium_arvense",
  "fabaceae/lupinus/lupinus_albus"
)
ph_phylomatic(taxa = taxa_str, phylo = r2)

# zanne2014 tree
taxa_str <- c(
  "zamiaceae/dioon/dioon_edule",
  "zamiaceae/encephalartos/encephalartos_dyerianus",
  "piperaceae/piper/piper_arboricola"
)
ph_phylomatic(taxa = taxa_str, phylo = zanne2014)

# zanne2014 subtree
zanne2014_subtr <- ape::extract.clade(zanne2014, node='Loganiaceae')
zanne_subtree_file <- tempfile(fileext = ".txt")
ape::write.tree(zanne2014_subtr, file = zanne_subtree_file)
taxa_str <- c(
  "loganiaceae/neuburgia/neuburgia_corynocarpum",
  "loganiaceae/geniostoma/geniostoma_borbonicum",
  "loganiaceae/strychnos/strychnos_darienensis"
)
ph_phylomatic(taxa = taxa_str, phylo = zanne2014_subtr)
ph_phylomatic(taxa = taxa_str, phylo = zanne_subtree_file)
```

## aot

```{r}
traits_file <- system.file("examples/traits_aot", package = "phylocomr")
phylo_file <- system.file("examples/phylo_aot", package = "phylocomr")
traitsdf_file <- system.file("examples/traits_aot_df", package = "phylocomr")
traits <- read.table(text = readLines(traitsdf_file), header = TRUE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(phylo_file)
ph_aot(traits = traits, phylo = phylo_str)
```

## bladj

```{r}
ages_file <- system.file("examples/ages", package = "phylocomr")
phylo_file <- system.file("examples/phylo_bladj", package = "phylocomr")
ages_df <- data.frame(
  a = c('malpighiales','salicaceae','fabaceae','rosales','oleaceae',
        'gentianales','apocynaceae','rubiaceae'),
  b = c(81,20,56,76,47,71,18,56)
)
phylo_str <- readLines(phylo_file)
(res <- ph_bladj(ages = ages_df, phylo = phylo_str))
plot(ape::read.tree(text = res))
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/phylocomr/issues).
* License: MIT
* Get citation information for `phylocomr` in R doing `citation(package = 'phylocomr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Introduction to the phylocomr package"
author: "Scott Chamberlain and Jeroen Ooms"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to the phylocomr package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

`phylocomr` is an R client for Phylocom - a C library for Analysis of Phylogenetic
Community Structure and Character Evolution.

Phylocom home page is at <http://phylodiversity.net/phylocom/>. The source code
for Phylocom is at <https://github.com/phylocom/phylocom/>.

Phylocom is usually used either on the command line or through the R package
`picante`, which has duplicated some of the Phylocom functionality.

Phylocom has been cited nearly 1000 times according to Google Scholar,
so is clearly a very widely used piece of software. The goal with this package
is to make it even easier to use - and in particular, to incorporate its use into
a reproducible workflow entirely in R instead of going to the shell/command line
for Phylocom usage. (Yes, some of Phylocom functionality is in `picante`, but
not all.)

In terms of performance, some functionality will be faster here than in `picante`,
but the maintainers of `picante` have re-written some Phylocom functionality
in C/C++, so performance should be similar in those cases.

## Install

Install `ape` for below examples:

```{r eval=FALSE}
install.packages('ape')
```

Stable `phylocomr` version from CRAN

```{r eval=FALSE}
install.packages("phylocomr")
```

Or, the development version from Github

```{r eval=FALSE}
devtools::install_github("ropensci/phylocomr")
```

```{r}
library("phylocomr")
```

## phylomatic

```r
taxa_file <- system.file("examples/taxa", package = "phylocomr")
phylo_file <- system.file("examples/phylo", package = "phylocomr")
(taxa_str <- readLines(taxa_file))
(phylo_str <- readLines(phylo_file))
ph_phylomatic(taxa = taxa_str, phylo = phylo_str)
```

## aot

```r
traits_file <- system.file("examples/traits_aot", package = "phylocomr")
phylo_file <- system.file("examples/phylo_aot", package = "phylocomr")
traitsdf_file <- system.file("examples/traits_aot_df", package = "phylocomr")
traits <- read.table(text = readLines(traitsdf_file), header = TRUE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(phylo_file)
ph_aot(traits = traits, phylo = phylo_str)
```

## bladj

```r
ages_file <- system.file("examples/ages", package = "phylocomr")
phylo_file <- system.file("examples/phylo_bladj", package = "phylocomr")
ages_df <- data.frame(
  a = c('malpighiales','salicaceae','fabaceae','rosales','oleaceae',
        'gentianales','apocynaceae','rubiaceae'),
  b = c(81,20,56,76,47,71,18,56)
)
phylo_str <- readLines(phylo_file)
res <- ph_bladj(ages = ages_df, phylo = phylo_str)

if (requireNamespace("ape")) {
  library(ape)
  plot(read.tree(text = res))
}
```

![plot](img/blad_tree.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bladj.R
\name{ph_bladj}
\alias{ph_bladj}
\title{bladj}
\usage{
ph_bladj(ages, phylo)
}
\arguments{
\item{ages}{(data.frame/character) ages data.frame, or path to an ages
file. required. column names do not matter, and are discarded anyway.
the first column must be the node names, and the second column the node
ages. See Details.}

\item{phylo}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}
}
\value{
newick string with attributes for where ages and phylo files
used are stored
}
\description{
Bladj take a phylogeny and fixes the root node at a specified age,
and fixes other nodes you might have age estimates for. It then sets all
other branch lengths by placing the nodes evenly between dated nodes,
and between dated nodes and terminals (beginning with the longest
'chains').
}
\details{
See \link{phylocomr-inputs} for expected input formats
}
\section{Common Errors}{

A few issues to be aware of:
\itemize{
\item the ages table must have a row for the root node with an age estimate.
bladj will not work without that. We attempt to check this but can only
check it if you pass in a phylo object; there's no easy way to parse a
newick string without requiring ape
\item bladj is case sensitive. internally we lowercase all tip and node labels
and taxon names in your ages file to avoid any case sensitivity problems
}
}

\examples{
\dontrun{
ages_file <- system.file("examples/ages", package = "phylocomr")
phylo_file <- system.file("examples/phylo_bladj", package = "phylocomr")

# from data.frame
ages_df <- data.frame(
  a = c('malpighiales','eudicots','ericales_to_asterales','plantaginaceae',
        'malvids', 'poales'),
  b = c(81, 20, 56, 76, 47, 71)
)
phylo_str <- readLines(phylo_file)
(res <- ph_bladj(ages = ages_df, phylo = phylo_str))
if (requireNamespace("ape")) {
  library(ape)
  plot(read.tree(text = res))
}

# from files
ages_file2 <- file.path(tempdir(), "ages")
write.table(ages_df, file = ages_file2, row.names = FALSE,
  col.names = FALSE, quote = FALSE)
phylo_file2 <- tempfile()
cat(phylo_str, file = phylo_file2, sep = '\n')
(res <- ph_bladj(ages_file2, phylo_file2))
if (requireNamespace("ape")) {
  library(ape)
  plot(read.tree(text = res))
}

# using a ape phylo phylogeny object
x <- read.tree(text = phylo_str)
if (requireNamespace("ape")) {
  library(ape)
  plot(x)
}

(res <- ph_bladj(ages_file2, x))
if (requireNamespace("ape")) {
  library(ape)
  tree <- read.tree(text = res)
  plot(tree)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylomatic.R
\name{ph_phylomatic}
\alias{ph_phylomatic}
\title{phylomatic}
\usage{
ph_phylomatic(taxa, phylo, tabular = FALSE, lowercase = FALSE, nodes = FALSE)
}
\arguments{
\item{taxa}{(character) all taxa as a character vector (will be written to
a temp file if provided) - OR a path to taxa file. Required. See Details.}

\item{phylo}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}

\item{tabular}{(logical) Output a tabular representation of phylogeny.
Default: \code{FALSE}}

\item{lowercase}{(logical) Convert all chars in taxa file to lowercase.
Default: \code{FALSE}}

\item{nodes}{(logical) label all nodes with default names.
Default: \code{FALSE}}
}
\description{
Phylomatic is a tool for extracting a phylogeny from a master
phylogeny using only a user-supplied list of taxa.
}
\details{
The \code{taxa} character vector must have each element of the
form \code{family/genus/genus_epithet}. If a file is passed in, each
line should have a \code{family/genus/genus_epithet} string - make sure
only one per line, and a newline (i.e., press ENTER) at the end of
each line
}
\examples{
\dontrun{
taxa_file <- system.file("examples/taxa", package = "phylocomr")
phylo_file <- system.file("examples/phylo", package = "phylocomr")

# from strings
(taxa_str <- readLines(taxa_file))
(phylo_str <- readLines(phylo_file))
(tree <- ph_phylomatic(taxa = taxa_str, phylo = phylo_str))

# from files
taxa_file2 <- tempfile()
cat(taxa_str, file = taxa_file2, sep = '\n')
phylo_file2 <- tempfile()
cat(phylo_str, file = phylo_file2, sep = '\n')
(tree <- ph_phylomatic(taxa = taxa_file2, phylo = phylo_file2))

if (requireNamespace("ape")) {
  library(ape)
  plot(read.tree(text = tree))
}
}
}
\references{
Phylomatic is also available as a web service
(https://github.com/camwebb/phylomatic-ws) - but is based on a different
code base (https://github.com/camwebb/phylomatic-ws)
See \href{https://doi.org/10.1111/j.1471-8286.2004.00829.x}{ Webb and Donoghue (2005)}
for more information on the goals of Phylomatic.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comstruct.R
\name{ph_comstruct}
\alias{ph_comstruct}
\title{comstruct}
\usage{
ph_comstruct(
  sample,
  phylo,
  null_model = 0,
  randomizations = 999,
  swaps = 1000,
  abundance = TRUE
)
}
\arguments{
\item{sample}{(data.frame/character) sample data.frame or path to a
sample file}

\item{phylo}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}

\item{null_model}{(integer) which null model to use. See Details.}

\item{randomizations}{(numeric) number of randomizations. Default: 999}

\item{swaps}{(numeric) number of swaps. Default: 1000}

\item{abundance}{(logical) If \code{TRUE} (default) computed accounting
for abundance. Otherwise, uses presence-absence.}
}
\value{
data.frame
}
\description{
Calculates mean phylogenetic distance (MPD) and mean nearest
phylogenetic taxon distance (MNTD; aka MNND) for each sample, and
compares them to MPD/MNTD values for randomly generated samples
(null communities) or phylogenies.
}
\details{
See \link{phylocomr-inputs} for expected input formats
}
\section{Null models}{

\itemize{
\item 0 - Phylogeny shuffle: This null model shuffles species labels across
the entire phylogeny. This randomizes phylogenetic relationships among
species.
\item 1 - Species in each sample become random draws from sample pool:
This null model maintains the species richness of each sample, but the
identities of the species occurring in each sample are randomized. For
each sample, species are drawn without replacement from the list of all
species actually occurring in at least one sample. Thus, species in the
phylogeny that are not actually observed to occur in a sample will not
be included in the null communities
\item 2 - Species in each sample become random draws from phylogeny pool:
This null model maintains the species richness of each sample, but the
identities of the species occurring in each sample are randomized. For
each sample, species are drawn without replacement from the list of all
species in the phylogeny pool. All species in the phylogeny will have
equal probability of being included in the null communities. By changing
the phylogeny, different species pools can be simulated. For example, the
phylogeny could include the species present in some larger region.
\item 3 - Independent swap: The independent swap algorithm (Gotelli and
Entsminger, 2003); also known as ‘SIM9’ (Gotelli, 2000) creates swapped
versions of the sample/species matrix.
}
}

\section{Taxon name case}{

In the \code{sample} table, if you're passing in a file, the names in the
third column must be all lowercase; if not, we'll lowercase them for you.
If you pass in a data.frame, we'll lowercase them for your. All phylo
tip/node labels are also lowercased to avoid any casing problems
}

\examples{
sfile <- system.file("examples/sample_comstruct", package = "phylocomr")
pfile <- system.file("examples/phylo_comstruct", package = "phylocomr")

# from data.frame
sampledf <- read.table(sfile, header = FALSE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(pfile)
(res <- ph_comstruct(sample = sampledf, phylo = phylo_str))

# from files
sample_str <- paste0(readLines(sfile), collapse = "\n")
sfile2 <- tempfile()
cat(sample_str, file = sfile2, sep = '\n')
pfile2 <- tempfile()
cat(phylo_str, file = pfile2, sep = '\n')
(res <- ph_comstruct(sample = sfile2, phylo = pfile2))

# different null models
ph_comstruct(sample = sfile2, phylo = pfile2, null_model = 0)
ph_comstruct(sample = sfile2, phylo = pfile2, null_model = 1)
ph_comstruct(sample = sfile2, phylo = pfile2, null_model = 2)
# ph_comstruct(sample = sfile2, phylo = pfile2, null_model = 3)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comtrait.R
\name{ph_comtrait}
\alias{ph_comtrait}
\title{comtrait}
\usage{
ph_comtrait(
  sample,
  traits,
  binary = NULL,
  metric = "variance",
  null_model = 0,
  randomizations = 999,
  abundance = TRUE
)
}
\arguments{
\item{sample}{(data.frame/character) sample data.frame or path to a
sample file}

\item{traits}{(data.frame/character) traits data.frame or path to a
traits file}

\item{binary}{(logical) a logical vector indicating what columns are to
be treated as binary characters - all others are treated as continuous}

\item{metric}{(integer) metric to calculate. One of variance, mpd, mntd,
or range (converted to phylocom integer format internally)}

\item{null_model}{(integer) which null model to use. See Details.}

\item{randomizations}{(numeric) number of randomizations. Default: 999}

\item{abundance}{(logical) If \code{TRUE} (default) computed accounting
for abundance. Otherwise, uses presence-absence.}
}
\value{
data.frame of the form:
\itemize{
\item trait - Trait name
\item sample - Sample name
\item ntaxa - Number of taxa in sample
\item mean - Mean value of trait in sample
\item metric - Observed metric in sample
\item meanrndmetric - Mean value of metric in null models
\item sdrndmetric - Standard deviation of metric in null models
\item sesmetric - Standardized effect size of metric
\item ranklow - Number of randomizations with metric lower than observed
\item rankhigh - Number of randomizations with metric higher than observed
\item runs - Number of randomizations
}
}
\description{
Calculate measures of trait dispersion within each community, and
compare observed patterns to those expected under a null model.
}
\details{
See \link{phylocomr-inputs} for expected input formats

If you give a data.frame to \code{traits} parameter it expects data.frame like
\itemize{
\item species - the taxon labels matching the sample data to \code{sample}
parameter
\item col1,col2,col3,etc. - any number of trait columns - column names do
not matter
}

When giving a data.frame to \code{traits} make sure to pass in a binary
vector for what traits are to be treated as binary.
}
\section{Null models}{

\itemize{
\item 0 - This null model shuffles trait values across species.
\item 1 - Species in each sample become random draws from sample pool.
This null model maintains the species richness of each sample, but
the identities of the species occurring in each sample are randomized.
For each sample, species are drawn without replacement from the list of
all species actually occurring in at least one sample
\item 2 - Species in each sample become random draws from traits data.
This null model maintains the species richness of each sample, but the
identities of the species occurring in each sample are randomized. For
each sample, species are drawn without replacement from the list of all
species with trait values. This function is redundant since by definition
the sample and trait species must match, but is included for consistency
with the comstruct function.
\item 3 - Independent swap: Same as for \link{ph_comdist} and
\link{ph_comstruct}
}
}

\section{Taxon name case}{

In the \code{sample} and \code{trait} tables, if you're passing in a file, the names
in the third and first columns, respectively, must be all lowercase; if not,
we'll lowercase them for you. If you pass in a data.frame, we'll lowercase
them for your. All phylo tip/node labels are also lowercased to avoid
any casing problems
}

\examples{
\dontrun{
sfile <- system.file("examples/sample_comstruct", package = "phylocomr")
tfile <- system.file("examples/traits_aot", package = "phylocomr")

# from files
sample_str <- paste0(readLines(sfile), collapse = "\n")
sfile2 <- tempfile()
cat(sample_str, file = sfile2, sep = '\n')

traits_str <- paste0(readLines(tfile), collapse = "\n")
tfile2 <- tempfile()
cat(traits_str, file = tfile2, sep = '\n')

ph_comtrait(sample = sfile2, traits = tfile2)

# from data.frame
sampledf <- read.table(sfile, header = FALSE,
  stringsAsFactors = FALSE)
traitsdf_file <- system.file("examples/traits_aot_df",
  package = "phylocomr")
traitsdf <- read.table(text = readLines(traitsdf_file), header = TRUE,
  stringsAsFactors = FALSE)
ph_comtrait(sample = sampledf, traits = traitsdf,
  binary = c(FALSE, FALSE, FALSE, TRUE))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylocomr-inputs.R
\name{phylocomr-inputs}
\alias{phylocomr-inputs}
\title{Expected inputs}
\description{
Expected inputs
}
\section{Ages data}{

A file or data.frame, with 2 columns:
\itemize{
\item (character) node taxonomic name
\item (numeric) age estimate for the node
}

If a file path is used, the table must not have headers

Applies to:
\itemize{
\item \code{\link[=ph_bladj]{ph_bladj()}}
}
}

\section{Sample data}{

A file or data.frame, with 3 columns, sorted by column 1, one row per taxon:
\itemize{
\item (character) sample plot/quadrat/trap/etc. name (no spaces, must begin with
a letter, not a number or symbol)
\item (integer) abundance (leave as 1 for presence/absence data)
\item (character) species code (same as in the phylogeny, must begin with
a letter, not a number or symbol)
}

If a file path is used, the table must not have headers, and must be
tab-delimited

Applies to:
\itemize{
\item \code{\link[=ph_comtrait]{ph_comtrait()}}
\item \code{\link[=ph_comstruct]{ph_comstruct()}}
\item \code{\link[=ph_comdist]{ph_comdist()}}
\item \code{\link[=ph_pd]{ph_pd()}}
\item \code{\link[=ph_rao]{ph_rao()}}
}
}

\section{Traits data}{

A tab-delimited file with the first line as
\verb{type<TAB>n<TAB>n<TAB>... [up to the number of traits]}, for example:
\verb{type<TAB>3<TAB>3<TAB>3<TAB>0}

where n indicates the type of trait in each of the four columns. Types:
\itemize{
\item \code{0}: binary (only one binary trait may be included, and it must be in
the first column) 1 for unordered multistate (no algorithms currently
implemented)
\item \code{2}: ordered multistate (currently treated as continuous)
\item \code{3}: continuous
}

Optional: The second line can start with the word name (lower case only)
and then list the names of the traits in order. These will appear in
the full output file

Subsequent lines should have the taxon name, which must be identical to
its appearance in phylo, and the data columns separated by tabs. For
example: \verb{sp1<TAB>1<TAB>1<TAB>1<TAB>0}
\itemize{
\item OR -
}

A data.frame, where the first column called \code{name} has each taxon, and
any number of columns withh traits, with each with the column name of
the trait. The first column name must be \code{name}.

Applies to:
\itemize{
\item \code{\link[=ph_comtrait]{ph_comtrait()}}
\item \code{\link[=ph_aot]{ph_aot()}}
}
}

\section{Phylogenies}{

Phylocom expects trees in Newick format. The basic Newick format used by
Phylocom is: \verb{((A,B),C);}. See the Phylocom manual
(http://phylodiversity.net/phylocom/) for more details on what they expect.

Applies to: all functions except \code{\link[=ph_phylomatic]{ph_phylomatic()}} and \code{\link[=ph_ecovolve]{ph_ecovolve()}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comdist.R
\name{ph_comdist}
\alias{ph_comdist}
\alias{ph_comdistnt}
\title{comdist}
\usage{
ph_comdist(
  sample,
  phylo,
  rand_test = FALSE,
  null_model = 0,
  randomizations = 999,
  abundance = TRUE
)

ph_comdistnt(
  sample,
  phylo,
  rand_test = FALSE,
  null_model = 0,
  randomizations = 999,
  abundance = TRUE
)
}
\arguments{
\item{sample}{(data.frame/character) sample data.frame or path to a
sample file}

\item{phylo}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}

\item{rand_test}{(logical) do you want to use null models?
Default: \code{FALSE}}

\item{null_model}{(integer) which null model to use. See Details.}

\item{randomizations}{(numeric) number of randomizations. Default: 999}

\item{abundance}{(logical) If \code{TRUE} (default) computed accounting
for abundance. Otherwise, uses presence-absence.}
}
\value{
data.frame or a list of data.frame's if use null models
}
\description{
Outputs the phylogenetic distance between samples, based on phylogenetic
distances of taxa in one sample to the taxa in the other
}
\details{
See \link{phylocomr-inputs} for expected input formats
}
\section{Null models}{

\itemize{
\item 0 - Phylogeny shuffle: This null model shuffles species labels across
the entire phylogeny. This randomizes phylogenetic relationships among
species.
\item 1 - Species in each sample become random draws from sample pool:
This null model maintains the species richness of each sample, but the
identities of the species occurring in each sample are randomized. For
each sample, species are drawn without replacement from the list of all
species actually occurring in at least one sample. Thus, species in the
phylogeny that are not actually observed to occur in a sample will not
be included in the null communities
\item 2 - Species in each sample become random draws from phylogeny pool:
This null model maintains the species richness of each sample, but the
identities of the species occurring in each sample are randomized. For
each sample, species are drawn without replacement from the list of all
species in the phylogeny pool. All species in the phylogeny will have
equal probability of being included in the null communities. By changing
the phylogeny, different species pools can be simulated. For example, the
phylogeny could include the species present in some larger region.
\item 3 - Independent swap: The independent swap algorithm (Gotelli and
Entsminger, 2003); also known as ‘SIM9’ (Gotelli, 2000) creates swapped
versions of the sample/species matrix.
}
}

\section{Taxon name case}{

In the \code{sample} table, if you're passing in a file, the names in the
third column must be all lowercase; if not, we'll lowercase them for you.
If you pass in a data.frame, we'll lowercase them for your. All phylo
tip/node labels are also lowercased to avoid any casing problems
}

\examples{
sfile <- system.file("examples/sample_comstruct", package = "phylocomr")
pfile <- system.file("examples/phylo_comstruct", package = "phylocomr")

# from data.frame
sampledf <- read.table(sfile, header = FALSE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(pfile)
ph_comdist(sample = sampledf, phylo = phylo_str)
ph_comdistnt(sample = sampledf, phylo = phylo_str)
ph_comdist(sample = sampledf, phylo = phylo_str, rand_test = TRUE)
ph_comdistnt(sample = sampledf, phylo = phylo_str, rand_test = TRUE)

# from files
sample_str <- paste0(readLines(sfile), collapse = "\n")
sfile2 <- tempfile()
cat(sample_str, file = sfile2, sep = '\n')
pfile2 <- tempfile()
cat(phylo_str, file = pfile2, sep = '\n')
ph_comdist(sample = sfile2, phylo = pfile2)
ph_comdistnt(sample = sfile2, phylo = pfile2)
ph_comdist(sample = sfile2, phylo = pfile2, rand_test = TRUE)
ph_comdistnt(sample = sfile2, phylo = pfile2, rand_test = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aot.R
\name{ph_aot}
\alias{ph_aot}
\title{aot}
\usage{
ph_aot(
  traits,
  phylo,
  randomizations = 999,
  trait_contrasts = 1,
  ebl_unstconst = FALSE
)
}
\arguments{
\item{traits}{(data.frame/character) trait data.frame or path to
traits file. required. See Details.}

\item{phylo}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}

\item{randomizations}{(numeric) number of randomizations. Default: 999}

\item{trait_contrasts}{(numeric) Specify which trait should be used as 'x'
variable for contrasts. Default: 1}

\item{ebl_unstconst}{(logical) Use equal branch lengths and unstandardized
contrasts. Default: \code{FALSE}}
}
\value{
a list of data.frames
}
\description{
AOT conducts univariate and bivariate tests of phylogenetic signal and
trait correlations, respectively, and node-level analyses of trait
means and diversification.
}
\details{
See \link{phylocomr-inputs} for expected input formats
}
\section{Taxon name case}{

In the \code{traits} table, if you're passing in a file, the names in the
first column must be all lowercase; if not, we'll lowercase them for you.
If you pass in a data.frame, we'll
lowercase them for your. All phylo tip/node labels are also lowercased
to avoid any casing problems
}

\examples{
\dontrun{
traits_file <- system.file("examples/traits_aot", package = "phylocomr")
phylo_file <- system.file("examples/phylo_aot", package = "phylocomr")

# from data.frame
traitsdf_file <- system.file("examples/traits_aot_df",
  package = "phylocomr")
traits <- read.table(text = readLines(traitsdf_file), header = TRUE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(phylo_file)
(res <- ph_aot(traits, phylo = phylo_str))

# from files
traits_str <- paste0(readLines(traits_file), collapse = "\n")
traits_file2 <- tempfile()
cat(traits_str, file = traits_file2, sep = '\n')
phylo_file2 <- tempfile()
cat(phylo_str, file = phylo_file2, sep = '\n')
(res <- ph_aot(traits_file2, phylo_file2))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/executables.R
\name{executables}
\alias{executables}
\alias{ecovolve}
\alias{phylocom}
\alias{phylomatic}
\title{Executables}
\usage{
ecovolve(args = "--help", intern = FALSE)

phylocom(args = "help", intern = FALSE)

phylomatic(args = "--help", intern = FALSE)
}
\arguments{
\item{args}{a character vector of arguments to command.}

\item{intern}{capture output as character vector. Default: \code{FALSE}}
}
\description{
Executables
}
\examples{
ecovolve()
phylocom()
phylomatic()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rao.R
\name{ph_rao}
\alias{ph_rao}
\title{rao - Rao's quadratic entropy}
\usage{
ph_rao(sample, phylo)
}
\arguments{
\item{sample}{(data.frame/character) sample data.frame or path to a
sample file}

\item{phylo}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}
}
\value{
A list of 6 data.frame's:
\strong{Diversity components}:
\itemize{
\item overall alpha (within-site)
\item beta (among-site)
\item total diversity
\item Fst statistic of differentiation for diversity and phylogenetic
diversity
}

\strong{Within-community diversity}:
\itemize{
\item Plot - Plot name
\item NSpp - Number of species
\item NIndiv - Number of individuals
\item PropIndiv - Proportion of all individuals found in this plot
\item D - Diversity (= Simpson’s diversity)
\item Dp - Phylogenetic diversity (= Diversity weighted by interspecific
phylogenetic distances)
}

The remaining 4 tables compare each community pairwise:
\itemize{
\item among_community_diversity_d - Among-community diversities
\item among_community_diversity_h - Among-community diversities excluding
within-community diversity
\item among_community_phylogenetic_diversity_dp - Among-community
phylogenetic diversities
\item among_community_phylogenetic_diversity_hp - Among-community
phylogenetic diversities excluding within-community diversity
}
}
\description{
A measure of within- and among-community diversity taking species
dissimilarity (phylogenetic dissimilarity) into account
}
\details{
See \link{phylocomr-inputs} for expected input formats
}
\section{Taxon name case}{

In the \code{sample} table, if you're passing in a file, the names
in the third column must be all lowercase; if not,
we'll lowercase them for you. If you pass in a data.frame, we'll lowercase
them for your. All phylo tip/node labels are also lowercased to avoid
any casing problems
}

\examples{
sfile <- system.file("examples/sample_comstruct", package = "phylocomr")
pfile <- system.file("examples/phylo_comstruct", package = "phylocomr")

# sample from data.frame, phylogeny from a string
sampledf <- read.table(sfile, header = FALSE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(pfile)

ph_rao(sample = sampledf, phylo = phylo_str)

# both from files
sample_str <- paste0(readLines(sfile), collapse = "\n")
sfile2 <- tempfile()
cat(sample_str, file = sfile2, sep = '\n')
pfile2 <- tempfile()
phylo_str <- readLines(pfile)
cat(phylo_str, file = pfile2, sep = '\n')

ph_rao(sample = sfile2, phylo = pfile2)
}
\seealso{
Other phylogenetic-diversity: 
\code{\link{ph_pd}()}
}
\concept{phylogenetic-diversity}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecovolve.R
\name{ph_ecovolve}
\alias{ph_ecovolve}
\title{ecovolve}
\usage{
ph_ecovolve(
  speciation = 0.05,
  extinction = 0.01,
  time_units = 100,
  out_mode = 3,
  prob_env = "3211000000",
  extant_lineages = FALSE,
  only_extant = FALSE,
  taper_change = NULL,
  competition = FALSE
)
}
\arguments{
\item{speciation}{(numeric) Probability of speciation per unit time.
Default: 0.05}

\item{extinction}{(numeric) Probability of extinction per unit time.
Default: 0.01}

\item{time_units}{(integer) Time units to simulate over. Default: 100}

\item{out_mode}{(integer) Output mode (2 = LTT; 3 = newick). Default: 3}

\item{prob_env}{(character) Probability envelope for character change.
must be a string of 10 integers. Default: 3211000000}

\item{extant_lineages}{(logical) Stop simulation after this number of
extant lineages. Default: \code{FALSE}}

\item{only_extant}{(logical) Output phylogeny pruned only for extant taxa.
Default: \code{FALSE}}

\item{taper_change}{(numeric/integer) Taper character change by
\code{e^(-time/F)}. This produces more conservatism in traits
(see Kraft et al., 2007). Default: \code{NULL}, not passed}

\item{competition}{(logical) Simulate competition, with trait proximity
increasing extinction. Default: \code{FALSE}}
}
\value{
a list with three elements:
\itemize{
\item phylogeny - a phylogeny as a newick string. In the case of
\code{out_mode = 2} gives a Lineage Through Time data.frame instead of a
newick phylogeny
\item sample - a data.frame with three columns, "sample" (all "alive"),
"abundance" (all 1's), "name" (the species code). In the case of
\code{out_mode = 2} gives an empty data.frame
\item traits - a data.frame with first column with spcies code ("name"),
then 5 randomly evolved and independent traits. In the case of
\code{out_mode = 2} gives an empty data.frame
}
}
\description{
Ecovolve generates a phylogeny via a random birth and death process,
generates a traits file with five randomly evolving, in-dependent traits,
and a sample file with a single sample unit (‘alive’) containing all
extant members of the phylogeny.
}
\section{Clean up}{

Two files, "ecovolve.sample" and "ecovolve.traits" are written to the
current working directory when this function runs - we read these files
in, then delete the files via \link{unlink}
}

\section{Failure behavior}{

Function occasionally fails with error "call to 'ecovolve' failed
with status 8. only 1 taxon; > 1 required" - this just means that only
1 taxon was created in the random process, so the function can't proceed
}

\examples{
\dontrun{
# ph_ecovolve(speciation = 0.05)
# ph_ecovolve(speciation = 0.1)
# ph_ecovolve(extinction = 0.005)
# ph_ecovolve(time_units = 50)
# ph_ecovolve(out_mode = 2)
# ph_ecovolve(extant_lineages = TRUE)
# ph_ecovolve(extant_lineages = FALSE)
# ph_ecovolve(only_extant = FALSE)
# ph_ecovolve(only_extant = TRUE, speciation = 0.1)
# ph_ecovolve(taper_change = 2)
# ph_ecovolve(taper_change = 10)
# ph_ecovolve(taper_change = 500)

if (requireNamespace("ape")) {
  # library(ape)
  # x <- ph_ecovolve(speciation = 0.05)
  # plot(read.tree(text = x$phylogeny))
}

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylocomr-package.R
\docType{package}
\name{phylocomr-package}
\alias{phylocomr-package}
\alias{phylocomr}
\title{Phylocom interface}
\description{
\code{phylocomr} gives you access to Phylocom, specifically the
Phylocom C library (https://github.com/phylocom/phylocom/),
licensed under BSD 2-clause
(http://www.opensource.org/licenses/bsd-license.php)
}
\details{
This package isn't doing system calls to a separately installed Phylocom
instance - but actually includes Phylocom itself in the package.

Phylocom is usually used either on the command line or through the
R package \pkg{picante}, which has duplicated some of the Phylocom
functionality.

In terms of performance, some functionality will be faster here than
in \code{picante}, but the maintainers of \code{picante} have re-written some
Phylocom functionality in C/C++, so performance should be similar in
those cases.
}
\section{A note about files}{

As a convienence you can pass ages, sample and trait data.frame's, and
phylogenies as strings, to \code{phylocomr} functions. However, \code{phylocomr}
has to write these data.frame's/strings to disk (your computer's
file system) to be able to run the Phylocom code on them. Internally,
\code{phylocomr} is writing to a temporary file to run Phylocom code, and
then the file is removed.

In addition, you can pass in files instead of data.frame's/strings.
These are not themselves used. Instead, we read and write those
files to temporary files. We do this for two reasons. First,
Phylocom expects the files its using to be in the same directory,
so if we control the file paths that becomes easier. Second,
Phylocom is case sensitive, so we simply standardize all taxon
names by lower casing all of them. We do this case manipulation
on the temporary files so that your original data files are
not modified.
}

\section{Package API}{

\itemize{
\item \code{\link[=ecovolve]{ecovolve()}}/\code{\link[=ph_ecovolve]{ph_ecovolve()}} - interface to \code{ecovolve} executable,
and a higher level interface
\item \code{\link[=phylomatic]{phylomatic()}}/\code{\link[=ph_phylomatic]{ph_phylomatic()}} - interface to \code{phylomatic}
executable, and a higher level interface
\item \code{\link[=phylocom]{phylocom()}} - interface to \code{phylocom} executable
\item \code{\link[=ph_aot]{ph_aot()}} - higher level interface to \code{aot}
\item \code{\link[=ph_bladj]{ph_bladj()}} - higher level interface to \code{bladj}
\item \code{\link[=ph_comdist]{ph_comdist()}}/\code{\link[=ph_comdistnt]{ph_comdistnt()}} - higher level interface to comdist
\item \code{\link[=ph_comstruct]{ph_comstruct()}} - higher level interface to comstruct
\item \code{\link[=ph_comtrait]{ph_comtrait()}} - higher level interface to comtrait
\item \code{\link[=ph_pd]{ph_pd()}} - higher level interface to Faith's phylogenetic diversity
}
}

\author{
Scott Chamberlain

Jeroen Ooms
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd.R
\name{ph_pd}
\alias{ph_pd}
\title{pd - Faith's index of phylogenetic diversity}
\usage{
ph_pd(sample, phylo)
}
\arguments{
\item{sample}{(data.frame/character) sample data.frame or path to a
sample file. required}

\item{phylo}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}
}
\value{
A single data.frame, with the colums:
\itemize{
\item sample - community name/label
\item ntaxa - number of taxa
\item pd - Faith's phylogenetic diversity
\item treebl - tree BL
\item proptreebl - proportion tree BL
}
}
\description{
Calculates Faith’s (1992) index of phylogenetic diversity (PD) for
each sample in the phylo.
}
\details{
See \link{phylocomr-inputs} for expected input formats
}
\section{Taxon name case}{

In the \code{sample} table, if you're passing in a file, the names
in the third column must be all lowercase; if not,
we'll lowercase them for you. If you pass in a data.frame, we'll lowercase
them for your. All phylo tip/node labels are also lowercased to avoid
any casing problems
}

\examples{
sfile <- system.file("examples/sample_comstruct", package = "phylocomr")
pfile <- system.file("examples/phylo_comstruct", package = "phylocomr")

# from data.frame
sampledf <- read.table(sfile, header = FALSE,
  stringsAsFactors = FALSE)
phylo_str <- readLines(pfile)
ph_pd(sample = sampledf, phylo = phylo_str)

# from files
sample_str <- paste0(readLines(sfile), collapse = "\n")
sfile2 <- tempfile()
cat(sample_str, file = sfile2, sep = '\n')
pfile2 <- tempfile()
phylo_str <- readLines(pfile)
cat(phylo_str, file = pfile2, sep = '\n')

ph_pd(sample = sfile2, phylo = pfile2)
}
\seealso{
Other phylogenetic-diversity: 
\code{\link{ph_rao}()}
}
\concept{phylogenetic-diversity}
