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
[Bradley Jones](https://github.com/brj1)
-------------
+ beast1 tree support
  - <https://github.com/GuangchuangYu/treeio/pull/4>

[Casey Dunn](https://github.com/caseywdunn)
----------
+ `drop.tip` method
  - <https://github.com/GuangchuangYu/ggtree/pull/90>
  - <https://github.com/GuangchuangYu/ggtree/pull/97>
  - <https://github.com/GuangchuangYu/treeio/pull/2>
+ add tip number in NHX annotation data
  - <https://github.com/GuangchuangYu/ggtree/pull/83>
+ use `textConnection` to pass string to the `file` parameter.
  - <https://github.com/GuangchuangYu/ggtree/pull/55>

[Shuangbin Xu](https://github.com/xiangpin)
----------
+ parse phyloxml file
  - <https://github.com/YuLab-SMU/treeio/pull/27>
  - <https://github.com/YuLab-SMU/treeio/pull/28>
+ parse `n` and `nm` in jplace when they coexists
  - <https://github.com/YuLab-SMU/treeio/pull/19>
+ jplace version 1 support
  - <https://github.com/YuLab-SMU/treeio/pull/25>

[Tyler Bradley](https://github.com/tbradley1013)
------------
+ `tree_subset` methods for `phylo` and `treedata` objects
  - <https://github.com/GuangchuangYu/treeio/pull/8>


[Konstantinos Geles](https://github.com/ConYel)
------------
+ prototype of `as.treedata` method for `pvclust` object


<!-- README.md is generated from README.Rmd. Please edit that file -->

#  treeio: Base classes and functions for phylogenetic tree input and output <a href="https://yulab-smu.top/treedata-book/"><img src="man/figures/logo.png" align="right" height="139" /></a>

[![](https://badges.ropensci.org/179_status.svg)](https://github.com/ropensci/onboarding/issues/179)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/treeio.svg)](https://www.bioconductor.org/packages/devel/bioc/html/treeio.html#since)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![platform](http://www.bioconductor.org/shields/availability/devel/treeio.svg)](https://www.bioconductor.org/packages/devel/bioc/html/treeio.html#archives)
[![codecov](https://codecov.io/gh/GuangchuangYu/treeio/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/treeio)

[![](https://img.shields.io/badge/release%20version-1.14.3-green.svg)](https://www.bioconductor.org/packages/treeio)
[![](https://img.shields.io/badge/devel%20version-1.15.2-green.svg)](https://github.com/guangchuangyu/treeio)
[![Linux Travis Build
Status](https://img.shields.io/travis/GuangchuangYu/treeio/master.svg?label=Linux)](https://travis-ci.org/GuangchuangYu/treeio)
[![AppVeyor Build
Status](https://img.shields.io/appveyor/ci/Guangchuangyu/treeio/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/treeio)

[![](https://img.shields.io/badge/download-68943/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/treeio)
[![](https://img.shields.io/badge/download-2590/month-blue.svg)](https://bioconductor.org/packages/stats/bioc/treeio)
[![download](http://www.bioconductor.org/shields/downloads/release/treeio.svg)](https://bioconductor.org/packages/stats/bioc/treeio)

‘treeio’ is an R package to make it easier to import and store
phylogenetic tree with associated data; and to link external data from
different sources to phylogeny. It also supports exporting phylogenetic
tree with heterogeneous associated data to a single tree file and can be
served as a platform for merging tree with associated data and
converting file formats.

Visit
<a href="https://yulab-smu.top/treedata-book/" class="uri">https://yulab-smu.top/treedata-book/</a>
for details.

[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=treeio&url=http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628/abstract&screen_name=guangchuangyu)
[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/GuangchuangYu)
[![](https://img.shields.io/badge/follow%20me%20on-WeChat-green.svg)](https://yulab-smu.top/images/biobabble.jpg)

:writing\_hand: Authors
-----------------------

Guangchuang YU

School of Basic Medical Sciences, Southern Medical University

<a href="https://yulab-smu.top" class="uri">https://yulab-smu.top</a>

If you use [treeio](http://bioconductor.org/packages/treeio) in
published research, please cite:

-   LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
    Jones, T Bradley, H Zhu, Y Guan, Y Jiang, **G Yu**<sup>\*</sup>.
    treeio: an R package for phylogenetic tree input and output with
    richly annotated and associated data. ***Molecular Biology and
    Evolution***. 2020, 37(2):599-603. doi:
    [10.1093/molbev/msz240](http://dx.doi.org/10.1093/molbev/msz240).

:arrow\_double\_down: Installation
----------------------------------

Get the released version from Bioconductor:

    ## try http:// if https:// URLs are not supported
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    ## BiocManager::install("BiocUpgrade") ## you may need this
    BiocManager::install("treeio")

Or the development version from github:

    ## install.packages("devtools")
    devtools::install_github("GuangchuangYu/treeio")

:sparkling\_heart: Contributing
-------------------------------

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).

:houses: Package Affiliations
-----------------------------

The `treeio` package is a part of the Bioconductor and rOpenSci
projects.

| [![bioconductor\_footer](http://bioconductor.org/images/logo_bioconductor.gif)](http://bioconductor.org) | [![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org) |
|:--------------------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
# treeio 1.19.1

+ bug fixed in `groupClade.treedata` to return a `treedata` object instead of `phylo` (2021-11-12, Fri)

# treeio 1.18.0

+ Bioconductor 3.14 release

# treeio 1.17.2

+ allow additional parameter to pass to `drop.tip` methods (2021-06-23, Wed, @xiangpin, #62)
+ `as.phylo` and `as.treedata` for `data.frame` (2021-06-12, Sat)
+ `as.ultrametric` method to force a tree to be ultrametric (2021-06-09, Wed)
+ introduce `force.ultrametric` parameter in `read.mcmctree` 

# treeio 1.17.1

+ `read.mcmctree` for PAML MCMCTree result (2021-06-04, Fri)

# treeio 1.16.0

+ Bioconductor 3.13 release

# treeio 1.15.6

+ optimized `read.nhx` for large tree file (2021-03-12, Fri)
- <https://github.com/YuLab-SMU/treeio/pull/51>

# treeio 1.15.5

+ `read.beast.newick` and `write.beast.newick` for importing and exporting newick text with metadata in BEAST style (2021-03-11, Thu)
  - <https://github.com/YuLab-SMU/treeio/pull/50>

# treeio 1.15.4

+ support parsing tree qza file from qiime2 (2020-03-01, Mon)
  - <https://github.com/YuLab-SMU/treeio/pull/46/files>

# treeio 1.15.3

+ support parsing phyloxml (2021-02-04, Thu)
  - <https://github.com/YuLab-SMU/treeio/pull/44>

# treeio 1.15.2

+ bug fixed of parsing nhx, now compatible with missing nhx tag (2020-11-19, Thu)
  - <https://github.com/YuLab-SMU/treeio/pull/40>

# treeio 1.15.1

+ remove magrittr::`%<>%` as it throw error of 'Error: could not find function "%>%<-"' (2020-11-19, Thu) 

# treeio 1.14.0

+ Bioconductor 3.12 release (2020-10-28, Wed)

# treeio 1.13.1

+ `as_tibble` for `pvclust` (2020-06-22, Mon)
+ `as.phylo` and `as.treedata` methods for `pvclust` object (2020-06-21, Sun)

# treeio 1.12.0

+ Bioconductor 3.11 release

# treeio 1.11.3

+ change according to dplyr (v=1.0.0) (2020-04-09, Thu)
  - remove mutate_, rename_, select_ and group_by_
+ remove data_frame for it was deprecated in tibble (v=3.0.0)

# treeio 1.11.2

+ update citation (2020-02-18, Tue)
+ phyloxml parser `read.phyloxml` (2019-12-05, Thu)
  
# treeio 1.11.1

+ support jplace version 1 (2019-11-25, Mon)
  - <https://github.com/YuLab-SMU/treeio/pull/25>
+ `offspring` return `integer(0)` instead of throw error if input `.node` is a tip (2019-11-21, Thu)

# treeio 1.10.0

+ Bioconductor 3.10 release

# treeio 1.9.3

+ add citation information (2019-10-05, Sta)
+ rename `phyPML` to `as.treedata.pml` (2019-10-01, Tue)
+ `as.phylo` method for `igraph` (only work with tree graph) (2019-09-28, Sat)

# treeio 1.9.2

+ `nodeid` and `nodelab` methods for converting between node number and labels (2019-08-09, Fri)
+ `parent`, 'ancestor`, `child`, `offspring` and `rootnode` methods for `treedata` (2019-08-07, Wed)
+ `read.mega_tabular` to parse MEGA Tabular output (2019-07-16, Tue)
+ `read.mega` to parse MEGA NEXUS (actually BEAST compatible)

# treeio 1.9.1

+ `rename_taxa` now use 1st column as key and 2nd column as value by default (2019-05-28, Tue)
+ enable `tree_subset` to specify `group_name` and enable to incorporate `root.edge` by setting `root_edge = TRUE` (2019-05-27, Mon)
+ `full_join` method for phylo object (2019-05-22, Wed)
+ redefined root method to wrape `ape::root.phylo` for compatibility (2019-05-20, Mon)
  - <https://github.com/GuangchuangYu/treeio/issues/18>
  
# treeio 1.8.0

+ Bioconductor 3.9 release

# treeio 1.7.4

+ update test according to the change of default RNG in the comming R-3.6 (2019-04-02, Tue)

# treeio 1.7.3

+ `rescale_tree` from `ggtree` (2019-01-11, Fri)

# treeio 1.7.2

+ `MRCA` methods for `phylo` and `treedata`  (2019-01-10, Thu)
+ mv vignettes to [treedata-book](https://yulab-smu.github.io/treedata-book/)
+ `root` method from `ggtree::reroot` (2018-12-28, Fri)
  - rename to `root` for importing `ape::root` generic

# treeio 1.7.1

+ compatible with `tibble` v=2.0.0 (2018-11-29, Thu)

# treeio 1.6.0

+ Bioconductor 3.8 release

# treeio 1.5.3

+ `read.jplace` compatible with output of [TIPars](https://github.com/id-bioinfo/TIPars) (2018-08-07, Tue)

# treeio 1.5.2

+ bug fixed of `as.phylo.ggtree` and `as.treedata.ggtree` (2018-07-19, Thu)
+ fixed R check for `tree_subset` by using `rlang::quo` and import `utils::head`
  and `utils::tail` (2018-05-24, Thu)
+ `tree_subset` methods contributed by [@tbradley1013](https://github.com/tbradley1013)
+ `drop.tip` works with `tree@extraInfo` (2018-05-23, Wed)
  - <https://github.com/GuangchuangYu/tidytree/pull/6#issuecomment-390259901>

# treeio 1.5.1

+ bug fixed of `groupOTU.treedata` (2018-05-23, Wed)
  - <https://github.com/GuangchuangYu/treeio/issues/7>

# treeio 1.4.0

+ Bioconductor 3.7 release

# treeio 1.3.15

+ Supports convert edge list (matrix, data.frame or tibble) to `phylo` and `treedata` object, now
  `ggtree` can be used to visualize all tree-like graph. (2018-04-23, Mon)

# treeio 1.3.14

+ rename_taxa (2018-04-19, Thu)
  - <https://guangchuangyu.github.io/2018/04/rename-phylogeny-tip-labels-in-treeio/>
+ read.astral (2018-04-17, Tue)
+ read.iqtree

# treeio 1.3.13

+ mv project website to <https://guangchuangyu.github.io/software/treeio>
+ update for rOpenSci acceptance
  - <https://github.com/ropensci/onboarding/issues/179#issuecomment-372127781>


# treeio 1.3.12

+ read.beast now compatible with taxa label contains ', " and space (2018-02-27,
  Wed)
+ update according to rOpenSci comments (2018-02-26, Mon)
  - <https://github.com/ropensci/onboarding/issues/179#issuecomment-365144565>
  - <https://github.com/ropensci/onboarding/issues/179#issuecomment-366800716>

# treeio 1.3.11

+ deprecate read.phyloT as read.tree in ape v5 now supports phyloT newick text
  <2018-01-11, Thu>
+ fixed goodpractice check <2018-01-10, Wed>
    - <https://github.com/ropensci/onboarding/issues/179#event-1416196637>
    - avoid using = for assignment
    - avoid code line > 80 characters
    - avoid sapply, instead using vapply and lapply
    - avoid using 1:length, 1:nrow and 1:ncol, use `seq_len` and `seq_along`
    - more unit tests

# treeio 1.3.10

* added 'Parsing jtree format' session in Importer vignette <2017-12-20, Wed>
* added 'Exporting tree data to JSON format' in Exporter vignette
* `read.jtree` and `write.jtree` functions
* added 'Combining tree with external data' and 'Merging tree data from
  different sources' sessions in Exporter vignette
* added 'Combining tree data' and 'Manipulating tree data using tidytree' sessions in Importer vignette
* full_join method for treedata object and added 'Linking external data to phylogeny' session in Importer vignette <2017-12-15, Fri>

# treeio 1.3.9

* move treedata class, show, get.fields methods to tidytree <2017-12-14, Thu>
* Exporter.Rmd vignette <2017-12-13, Wed>

# treeio 1.3.8

* mv treeio.Rmd vignette to Importer.Rmd and update the contents <2017-12-13, Wed>
* write.beast for treedata object <2017-12-12, Tue>
* add "connect" parameter in groupOTU <2017-12-12, Tue>
   + <https://groups.google.com/forum/#!msg/bioc-ggtree/Q4LnwoTf1DM/yEe95OFfCwAJ>

# treeio 1.3.7

* export groupClade.phylo method <2017-12-11, Mon>

# treeio 1.3.6

* re-defined groupOTU and groupClade generic using S3 <2017-12-11, Mon>

# treeio 1.3.5

* parent, ancestor, child, offspring, rootnode and sibling generic and method for phylo <2017-12-11, Mon>
* update mask and merge_tree function according to the treedata object <2017-12-11, Mon>

# treeio 1.3.4

* support tbl_tree object defined in tidytree <2017-12-08, Fri>

# treeio 1.3.3

* read.codeml output treedata, remove codeml class and clean up code <2017-12-07, Thu>

# treeio 1.3.2

* read.codeml_mlc output treedata object and remove codeml_mlc class <2017-12-06, Wed>
* read.paml_rst output treedata and remove paml_rst class <2017-12-06, Wed>
* read.phylip.tree and read.phylip.seq
* read.phylip output treedata object and phylip class definition was removed
* read.hyphy output treedata object; hyphy class definition was removed
* remove r8s class, read.r8s now output multiPhylo object
* jplace class inherits treedata <2017-12-05, Tue>
* using treedata object to store beast and mrbayes tree
* export read.mrbayes

# treeio 1.3.1

* compatible to parse beast output that only contains HPD range <2017-11-01, Wed>
   + https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/bioc-ggtree/RF2Ly52U_gc/jEP97nNPAwAJ

# treeio 1.2.0

* BioC 3.6 release <2017-11-01, Wed>

# treeio 1.1.2

* new project site using blogdown <2017-09-28, Thu>

# treeio 1.1.1

* parse mlc file without dNdS <2017-08-31, Thu>
   + https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!topic/bioc-ggtree/hTRj-uldgAg
* better implementation of merge_tree <2017-08-31, Thu>

# treeio 0.99.11

* bug fixed in get.fields method for paml_rst <2017-03-20, Mon>
* fixed raxml2nwk for using treedata as output of read.raxml <2017-03-17, Fri>
* taxa_rename function <2017-03-15, Wed>
* phyPML method moved from ggtree <2017-03-06, Mon>

# treeio 0.99.10

* remove raxml class, now read.raxml output treedata object <2017-02-28, Tue>
* bug fixed of read.beast <2017-02-27, Mon>

# treeio 0.99.9

* read.newick for parsing node.label as support values <2017-01-03, Tue>
* read.beast support MrBayes output <2016-12-30, Fri>
* export as.phylo.ggtree <2016-12-30, Fri>

# treeio 0.99.8

* as.treedata.ggtree <2016-12-30, Fri>
* as.treedata.phylo4 & as.treedata.phylo4d <2016-12-28, Wed>

# treeio 0.99.7

* groupOTU, groupClade, gzoom methods from ggtree <2016-12-21, Wed>

# treeio 0.99.6

* add unit test of NHX (move from ggtree) <2016-12-14, Wed>

# treeio 0.99.3

* fixed BiocCheck by adding examples <2016-12-07, Wed>

# treeio 0.99.1

* fixed link in DESCRIPTION <2016-12-06, Tue>

# treeio 0.99.0

* add vignette <2016-12-06, Tue>
* move parser functions from ggtree <2016-12-06, Tue>

# treeio 0.0.1

* `read.nhx` from ggtree <2016-12-06, Tue>
* `as.phylo.treedata` to access `phylo` from `treedata` object <2016-12-06, Tue>
* `as.treedata.phylo` to convert `phylo` to `treedata` object <2016-12-06, Tue>
* `treedata` class definition <2016-12-06, Tue>
# TODO LIST

+ [ ] improve read.phyloxml
+ [ ] re-write read.beast to optimize parsing large file
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/GuangchuangYu/treeio/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/{repo}.git`
* Make sure to track progress upstream (i.e., on our version of `treeio` at `GuangchuangYu/treeio`) by doing `git remote add upstream https://github.com/GuangchuangYu/treeio.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `GuangchuangYu/treeio`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email?

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
### Prerequisites

+ [ ] Have you read [Feedback](https://guangchuangyu.github.io/software/treeio/#feedback) and follow the [guide](https://guangchuangyu.github.io/2016/07/how-to-bug-author/)?
	* [ ] make sure your are using the latest release version
	* [ ] read the [documents](https://guangchuangyu.github.io/software/treeio/documentation/)
	* [ ] google your quesion/issue

### Describe you issue

* [ ] Make a reproducible example (*e.g.* [1](https://gist.github.com/talonsensei/e1fad082657054207f249ec98f0920eb))
* [ ] your code should contain comments to describe the problem (*e.g.* what expected and actually happened?)


### Ask in right place

* [ ] for bugs or feature requests, post here (github issue)
* [ ] for questions, please post to [google group](https://groups.google.com/forum/#!forum/bioc-ggtree)


<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>

---
output:
  md_document:
    variant: gfm
html_preview: false
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


#  treeio: Base classes and functions for phylogenetic tree input and output <a href="https://yulab-smu.top/treedata-book/"><img src="man/figures/logo.png" align="right" height="139" /></a>

```{r echo=FALSE, results="hide", message=FALSE}
#library("txtplot")
library("badger")
library("ypages")

Biocpkg <- function (pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}
```

[![](https://badges.ropensci.org/179_status.svg)](https://github.com/ropensci/onboarding/issues/179)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/treeio.svg)](https://www.bioconductor.org/packages/devel/bioc/html/treeio.html#since)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![platform](http://www.bioconductor.org/shields/availability/devel/treeio.svg)](https://www.bioconductor.org/packages/devel/bioc/html/treeio.html#archives)
[![codecov](https://codecov.io/gh/GuangchuangYu/treeio/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/treeio)


`r badge_bioc_release("treeio", "green")`
`r badge_devel("guangchuangyu/treeio", "green")`
[![Linux Travis Build Status](https://img.shields.io/travis/GuangchuangYu/treeio/master.svg?label=Linux)](https://travis-ci.org/GuangchuangYu/treeio)
[![AppVeyor Build Status](https://img.shields.io/appveyor/ci/Guangchuangyu/treeio/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/treeio)


`r badge_bioc_download("treeio", "total", "blue")`
`r badge_bioc_download("treeio", "month", "blue")`
`r badge_bioc_download_rank("treeio")`



```{r comment="", echo=FALSE, results='asis'}
cat(packageDescription('treeio')$Description)
```

Visit <https://yulab-smu.top/treedata-book/> for details.


[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=treeio&url=http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628/abstract&screen_name=guangchuangyu)
[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/GuangchuangYu)
`r badger::badge_custom("follow me on", "WeChat", "green", "https://yulab-smu.top/images/biobabble.jpg")`




## :writing_hand: Authors

Guangchuang YU 

School of Basic Medical Sciences, Southern Medical University

<https://yulab-smu.top>

If you use `r Biocpkg('treeio')` in published research, please cite:

+ LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR Jones, T Bradley, H Zhu, Y Guan, Y Jiang, __G Yu__^\*^. treeio: an R package for phylogenetic tree input and output with richly annotated and associated data. __*Molecular Biology and Evolution*__. 2020, 37(2):599-603. doi: [10.1093/molbev/msz240](http://dx.doi.org/10.1093/molbev/msz240).


## :arrow_double_down: Installation

Get the released version from Bioconductor:

```r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("treeio")
```

Or the development version from github:

```r
## install.packages("devtools")
devtools::install_github("GuangchuangYu/treeio")
```



## :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree to
abide by the terms outlined in the [Contributor Code of Conduct](CONDUCT.md).


## :houses: Package Affiliations

The `treeio` package is a part of the Bioconductor and rOpenSci projects.

| [![bioconductor_footer](http://bioconductor.org/images/logo_bioconductor.gif)](http://bioconductor.org) | [![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org) |
|:-------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------:|

---
title: "treeio: Base Classes and Functions for Phylogenetic Tree Input and Output"
author: "Guangchuang Yu\\

        School of Basic Medical Sciences, Southern Medical University"
date: "`r Sys.Date()`"
output:
  prettydoc::html_pretty:
    toc: true
    theme: cayman
    highlight: github
  pdf_document:
    toc: true
vignette: >
  %\VignetteIndexEntry{treeio}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r style, echo=FALSE, results="asis", message=FALSE}
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)
```




Please go to <https://yulab-smu.github.io/treedata-book/> (the first three chapters) for the full vignette.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-reroot.R
\name{root.treedata}
\alias{root.treedata}
\title{root}
\usage{
\method{root}{treedata}(phy, outgroup, node = NULL, edgelabel = TRUE, ...)
}
\arguments{
\item{phy}{tree object}

\item{outgroup}{a vector of mode numeric or character specifying the new outgroup}

\item{node}{node to reroot}

\item{edgelabel}{a logical value specifying whether to treat node labels as
edge labels and thus eventually switching them so that they are associated
with the correct edges.}

\item{...}{additional parameters passed to ape::root.phylo}
}
\value{
rerooted treedata
}
\description{
re-root a tree
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloxml.R
\name{read.phyloxml}
\alias{read.phyloxml}
\title{read.phyloxml}
\usage{
read.phyloxml(file)
}
\arguments{
\item{file}{phyloxml file}
}
\value{
treedata class or treedataList class
}
\description{
read.phyloxml
}
\examples{
xmlfile1 <- system.file("extdata/phyloxml", "test_x2.xml", package="treeio")
px1 <- read.phyloxml(xmlfile1)
px1
xmlfile2 <- system.file("extdata/phyloxml", "phyloxml_examples.xml", package="treeio")
px2 <- read.phyloxml(xmlfile2)
px2
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.treedataList}
\alias{print.treedataList}
\title{print}
\usage{
\method{print}{treedataList}(x, ...)
}
\arguments{
\item{x}{a list of treedata objects}

\item{...}{no used}
}
\value{
message
}
\description{
print information of a list of treedata objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rescale_tree.R
\name{rescale_tree}
\alias{rescale_tree}
\title{rescale_tree}
\usage{
rescale_tree(tree_object, branch.length)
}
\arguments{
\item{tree_object}{tree object}

\item{branch.length}{numerical features (e.g. dN/dS)}
}
\value{
update tree object
}
\description{
rescale branch length of tree object
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-as-phylo.R
\name{get.tree}
\alias{get.tree}
\title{get.tree}
\usage{
get.tree(x, ...)
}
\arguments{
\item{x}{tree object}

\item{...}{additional parameters}
}
\value{
phylo object
}
\description{
access phylo slot
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_tree.R
\name{merge_tree}
\alias{merge_tree}
\title{merge_tree}
\usage{
merge_tree(obj1, obj2)
}
\arguments{
\item{obj1}{tree object 1}

\item{obj2}{tree object 2}
}
\value{
tree object
}
\description{
merge two tree object
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/jplace.R
\name{get.placements}
\alias{get.placements}
\alias{get.placements.jplace}
\title{get.placements}
\usage{
get.placements(tree, ...)

\method{get.placements}{jplace}(tree, by = "best", ...)
}
\arguments{
\item{tree}{tree object}

\item{...}{additional parameters}

\item{by}{one of 'best' and 'all'}
}
\value{
placement tibble
}
\description{
access placement information
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequence-utilities.R
\name{read.fasta}
\alias{read.fasta}
\title{read.fasta}
\usage{
read.fasta(fasta, type = "auto")
}
\arguments{
\item{fasta}{fasta file}

\item{type}{sequence type of the input file, one of 'NT' or 'AA'.
Default is 'auto' and guess the sequence type automatically}
}
\value{
DNAbin or AAbin object
}
\description{
read FASTA file
}
\details{
This function supports both DNA or AA sequences
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylip.R
\name{read.phylip.seq}
\alias{read.phylip.seq}
\title{read.phylip.seq}
\usage{
read.phylip.seq(file)
}
\arguments{
\item{file}{phylip file, currently only sequential format is supported}
}
\value{
DNAbin object
}
\description{
read aligned sequences from phylip format
}
\references{
\url{http://evolution.genetics.washington.edu/phylip/doc/sequence.html}
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newick.R
\name{read.newick}
\alias{read.newick}
\title{read.newick}
\usage{
read.newick(file, node.label = "label", ...)
}
\arguments{
\item{file}{newick file}

\item{node.label}{parse node label as 'label' or 'support' value}

\item{...}{additional parameter, passed to 'read.tree'}
}
\value{
phylo or treedata object
}
\description{
read newick tree
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{jplace-class}
\alias{jplace-class}
\title{Class "jplace"
This class stores phylogenetic placements}
\description{
Class "jplace"
This class stores phylogenetic placements
}
\section{Slots}{

\describe{
\item{\code{phylo}}{phylo object for tree structure}

\item{\code{treetext}}{newick tree string}

\item{\code{data}}{associated data}

\item{\code{extraInfo}}{extra information, reserve for merge_tree}

\item{\code{file}}{tree file}

\item{\code{placements}}{reserve for jplace file to store placement information}

\item{\code{info}}{extra information, e.g. metadata, software version etc.}
}}

\author{
Guangchuang Yu \url{https://guangchuangyu.github.io}
}
\keyword{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jtree.R
\name{read.jtree}
\alias{read.jtree}
\title{read.jtree}
\usage{
read.jtree(file)
}
\arguments{
\item{file}{tree file}
}
\value{
treedata object
}
\description{
Import tree data from jtree file, which is JSON-based text and probably output by write.jtree
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylip.R
\name{read.phylip.tree}
\alias{read.phylip.tree}
\title{read.phylip.tree}
\usage{
read.phylip.tree(file)
}
\arguments{
\item{file}{phylip file}
}
\value{
phylo or multiPhylo object
}
\description{
parse tree from phylip file
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RAxML.R
\name{read.raxml}
\alias{read.raxml}
\title{read.raxml}
\usage{
read.raxml(file)
}
\arguments{
\item{file}{RAxML bootstrapping analysis output}
}
\value{
treedata object
}
\description{
parse RAxML bootstrapping analysis output
}
\examples{
raxml_file <- system.file("extdata/RAxML", "RAxML_bipartitionsBranchLabels.H3", package="treeio")
read.raxml(raxml_file)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mask.R
\name{mask}
\alias{mask}
\title{mask}
\usage{
mask(tree_object, field, site, mask_site = FALSE)
}
\arguments{
\item{tree_object}{tree object}

\item{field}{selected field}

\item{site}{site}

\item{mask_site}{if TRUE, site will be masked.
if FALSE, selected site will not be masked, while other sites will be masked.}
}
\value{
updated tree object
}
\description{
site mask
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-as-ultrametric.R
\name{as.ultrametric}
\alias{as.ultrametric}
\title{as.ultrametric}
\usage{
as.ultrametric(tree, ...)
}
\arguments{
\item{tree}{tree object}

\item{...}{additional parameters}
}
\value{
treedata or phylo object
}
\description{
as.ultrametric
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-utilities.R
\name{label_branch_paml}
\alias{label_branch_paml}
\title{label_branch_paml}
\usage{
label_branch_paml(tree, node, label)
}
\arguments{
\item{tree}{phylo object}

\item{node}{node number}

\item{label}{label of branch, e.g. #1}
}
\value{
updated phylo object
}
\description{
label branch for PAML to infer selection pressure using branch model
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nhx.R
\name{read.nhx}
\alias{read.nhx}
\title{read.nhx}
\usage{
read.nhx(file)
}
\arguments{
\item{file}{nhx file}
}
\value{
nhx object
}
\description{
read nhx tree file
}
\examples{
nhxfile <- system.file("extdata/NHX", "ADH.nhx", package="treeio")
read.nhx(nhxfile)
}
\author{
Guangchuang Yu \url{https://guangchuangyu.github.io}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-utilities.R
\name{getNodeNum}
\alias{getNodeNum}
\alias{Nnode2}
\title{getNodeNum}
\usage{
getNodeNum(tree)

Nnode2(tree)
}
\arguments{
\item{tree}{tree object}
}
\value{
number
}
\description{
calculate total number of nodes
}
\examples{
getNodeNum(rtree(30))
Nnode2(rtree(30))
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jplace.R
\name{read.jplace}
\alias{read.jplace}
\title{read.jplace}
\usage{
read.jplace(file)
}
\arguments{
\item{file}{jplace file}
}
\value{
\code{jplace} instance
}
\description{
read jplace file
}
\examples{
jp <- system.file("extdata", "sample.jplace", package="treeio")
read.jplace(jp)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/method-drop-tip.R
\docType{methods}
\name{drop.tip}
\alias{drop.tip}
\alias{drop.tip,treedata-method}
\alias{drop.tip,treedata}
\alias{drop.tip,phylo-method}
\alias{drop.tip,phylo}
\title{drop.tip method}
\source{
drop.tip for phylo object is a wrapper method of ape::drop.tip
from the ape package. The documentation you should
read for the drop.tip function can be found here: \link[ape]{drop.tip}
}
\usage{
drop.tip(object, tip, ...)

drop.tip(object, tip, ...)

\S4method{drop.tip}{phylo}(object, tip, ...)
}
\arguments{
\item{object}{A treedata or phylo object}

\item{tip}{a vector of mode numeric or character specifying the tips to delete}

\item{...}{additional parameters}
}
\value{
updated object
}
\description{
drop.tip method

drop.tip method
}
\examples{
nhxfile <- system.file("extdata/NHX", "ADH.nhx", package="treeio")
nhx <- read.nhx(nhxfile)
drop.tip(nhx, c("ADH2", "ADH1"))
}
\seealso{
\link[ape]{drop.tip}
}
\author{
Casey Dunn \url{http://dunnlab.org}  and Guangchuang Yu \url{https://guangchuangyu.github.io}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paml_rst.R
\name{read.paml_rst}
\alias{read.paml_rst}
\title{read.paml_rst}
\usage{
read.paml_rst(rstfile, type = "Joint")
}
\arguments{
\item{rstfile}{rst file}

\item{type}{one of 'Marginal' or 'Joint'}
}
\value{
A \code{treedata} object
}
\description{
read rst file from paml (both baseml and codeml) output
}
\examples{
rstfile <- system.file("extdata/PAML_Baseml", "rst", package="treeio")
read.paml_rst(rstfile)
}
\author{
Guangchuang Yu \url{https://guangchuangyu.github.io}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylip.R
\name{read.phylip}
\alias{read.phylip}
\title{read.phylip}
\usage{
read.phylip(file)
}
\arguments{
\item{file}{phylip file}
}
\value{
an instance of 'phylip'
}
\description{
parsing phylip tree format
}
\examples{
phyfile <- system.file("extdata", "sample.phy", package="treeio")
read.phylip(phyfile)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codeml.R
\name{read.codeml}
\alias{read.codeml}
\title{read.codeml}
\usage{
read.codeml(rstfile, mlcfile, tree = "mlc", type = "Joint")
}
\arguments{
\item{rstfile}{rst file}

\item{mlcfile}{mlc file}

\item{tree}{one of 'mlc' or 'rst'}

\item{type}{one of 'Marginal' or 'Joint'}
}
\value{
A \code{treedata} object
}
\description{
read baseml output
}
\examples{
rstfile <- system.file("extdata/PAML_Codeml", "rst", package="treeio")
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="treeio")
read.codeml(rstfile, mlcfile)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r8s.R
\name{read.r8s}
\alias{read.r8s}
\title{read.r8s}
\usage{
read.r8s(file)
}
\arguments{
\item{file}{r8s output log file}
}
\value{
multiPhylo object
}
\description{
parse output from r8s
}
\examples{
read.r8s(system.file("extdata/r8s", "H3_r8s_output.log", package="treeio"))
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCTree.R
\name{read.mcmctree}
\alias{read.mcmctree}
\title{read.mcmctree}
\usage{
read.mcmctree(file, force.ultrametric = FALSE)
}
\arguments{
\item{file}{the output tree file of MCMCTree}

\item{force.ultrametric}{logical whether convert the tree 
to be ultrametric, if it is not ultrametric, default is FALSE.
When the tree is ultrametric, branch times will be calculated 
automatically.}
}
\value{
treedata object
}
\description{
read MCMCTree output Tree
}
\examples{
file <- system.file("extdata/MCMCTree", "mcmctree_output.tree", package="treeio")
tr <- read.mcmctree(file)
tr
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raxml2nwk.R
\name{raxml2nwk}
\alias{raxml2nwk}
\title{raxml2nwk}
\usage{
raxml2nwk(infile, outfile = "raxml.tree")
}
\arguments{
\item{infile}{input file}

\item{outfile}{output file}
}
\value{
newick file
}
\description{
convert raxml bootstrap tree to newick format
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyphy.R
\name{read.hyphy.seq}
\alias{read.hyphy.seq}
\title{read.hyphy.seq}
\usage{
read.hyphy.seq(file)
}
\arguments{
\item{file}{output of hyphy ancestral sequence inference; nexus format}
}
\value{
DNAbin object
}
\description{
parse sequences from hyphy output
}
\examples{
ancseq <- system.file("extdata/HYPHY", "ancseq.nex", package="treeio")
read.hyphy.seq(ancseq)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASTRAL.R
\name{read.astral}
\alias{read.astral}
\title{read.astral}
\usage{
read.astral(file)
}
\arguments{
\item{file}{ASTRAL Newick file}
}
\value{
treedata object
}
\description{
parse ASTRAL output newick text
}
\examples{
tt <- paste0(
  "((species1,(species2,species3)'[pp1=0.75;pp2=0.24;pp3=0.01]':",
  "1.2003685744180805)'[pp1=0.98;pp2=0.02;pp3=0]':0.9679599282730038,",
  "((species4,species5)'[pp1=0.88;pp2=0.11;pp3=0.01]':1.2454851536484994))"
)
read.astral(textConnection(tt))
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iqtree.R
\name{read.iqtree}
\alias{read.iqtree}
\title{read.iqtree}
\usage{
read.iqtree(file)
}
\arguments{
\item{file}{IQ-TREE Newick text}
}
\value{
treedata object
}
\description{
parse IQ-TREE output
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename.R
\name{rename_taxa}
\alias{rename_taxa}
\title{rename_taxa}
\usage{
rename_taxa(tree, data, key = 1, value = 2)
}
\arguments{
\item{tree}{tree object, either treedata or phylo}

\item{data}{data frame}

\item{key}{column in data that match tip label (use 1st column by default)}

\item{value}{column in data for rename tip label (use 2nd column by default)}
}
\value{
tree object
}
\description{
rename tip label of phylogenetic tree
}
\examples{
tree <- rtree(3)
d <- data.frame(old = paste0('t', 1:3), new = LETTERS[1:3])
rename_taxa(tree, d)
rename_taxa(tree, d, old, new)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beast.R, R/mega.R
\name{read.beast}
\alias{read.beast}
\alias{read.mrbayes}
\alias{read.beast.newick}
\alias{read.mega}
\title{read.beast}
\usage{
read.beast(file)

read.mrbayes(file)

read.beast.newick(file)

read.mega(file)
}
\arguments{
\item{file}{newick file}
}
\value{
treedata object

treedata object
}
\description{
read beast/mrbayes/mega Nexus output

read beast/mrbayes/mega newick file format
}
\examples{
file <- system.file("extdata/BEAST", "beast_mcc.tree", package="treeio")
read.beast(file)
file <- system.file("extdata/MrBayes", "Gq_nxs.tre", package="treeio")
read.mrbayes(file)
tree <- read.beast.newick(textConnection('(a[&rate=1]:2,(b[&rate=1.1]:1,c[&rate=0.9]:1)[&rate=1]:1);'))
}
\author{
Guangchuang Yu \url{https://guangchuangyu.github.io}

Bradley R Jones
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeqza.R
\name{read.treeqza}
\alias{read.treeqza}
\title{read.treeqza}
\usage{
read.treeqza(treeqza, node.label = "label", ...)
}
\arguments{
\item{treeqza}{the qiime2 output file contained tree file.}

\item{node.label}{parse node label as 'label' or 'support' value.}

\item{...}{additional parameter, passed to 'read.tree'.}
}
\value{
phylo tree object or treedata object when node.label was parsed 'support'.
}
\description{
read.treeqza
}
\examples{
qzafile1 <- system.file("extdata/qiime2treeqza", "fasttree-tree.qza", package="treeio")
qzafile2 <- system.file("extdata/qiime2treeqza", "iqt-tree.qza", package="treeio")
qzafile3 <- system.file("extdata/qiime2treeqza", "raxml-cat-tree.qza", package="treeio")
tr1 <- read.treeqza(qzafile1)
tr1
tr2 <- read.treeqza(qzafile2)
tr2
tr3 <- read.treeqza(qzafile3)
tr3
# parse node label as 'support' value.
qzafile4 <- system.file("extdata/qiime2treeqza", "raxml-cat-bootstrap-tree.qza", package="treeio")
tr4 <- read.treeqza(qzafile4, node.label="support")
tr4
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/ape.R, R/method-as-phylo.R,
%   R/reexport.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{treedata}
\alias{read.tree}
\alias{read.nexus}
\alias{rtree}
\alias{write.tree}
\alias{write.nexus}
\alias{Nnode}
\alias{Ntip}
\alias{is.rooted}
\alias{root}
\alias{as.phylo}
\alias{\%>\%}
\alias{get.fields}
\alias{get.data}
\alias{as.treedata}
\alias{ancestor}
\alias{parent}
\alias{child}
\alias{offspring}
\alias{rootnode}
\alias{nodeid}
\alias{nodelab}
\alias{MRCA}
\alias{full_join}
\alias{as_tibble}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{ape}{\code{\link[ape:summary.phylo]{Nnode}}, \code{\link[ape:summary.phylo]{Ntip}}, \code{\link[ape]{as.phylo}}, \code{\link[ape:root]{is.rooted}}, \code{\link[ape]{read.nexus}}, \code{\link[ape]{read.tree}}, \code{\link[ape]{root}}, \code{\link[ape]{rtree}}, \code{\link[ape]{write.nexus}}, \code{\link[ape]{write.tree}}}

  \item{dplyr}{\code{\link[dplyr:mutate-joins]{full_join}}}

  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}

  \item{tibble}{\code{\link[tibble]{as_tibble}}}

  \item{tidytree}{\code{\link[tidytree]{MRCA}}, \code{\link[tidytree]{ancestor}}, \code{\link[tidytree]{as.treedata}}, \code{\link[tidytree]{child}}, \code{\link[tidytree:get.data-methods]{get.data}}, \code{\link[tidytree:get.fields-methods]{get.fields}}, \code{\link[tidytree]{nodeid}}, \code{\link[tidytree]{nodelab}}, \code{\link[tidytree]{offspring}}, \code{\link[tidytree]{parent}}, \code{\link[tidytree]{rootnode}}, \code{\link[tidytree]{treedata}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-as-treedata.R, R/phangorn.R
\name{as.treedata.phylo}
\alias{as.treedata.phylo}
\alias{as.treedata.pml}
\title{as.treedata}
\usage{
\method{as.treedata}{phylo}(tree, boot = NULL, ...)

\method{as.treedata}{pml}(tree, type = "ml", ...)
}
\arguments{
\item{tree}{input tree, a \code{phylo} object}

\item{boot}{optional, can be bootstrap value from ape::boot.phylo}

\item{...}{additional parameters}

\item{type}{one of 'ml' and 'bayes' for inferring ancestral sequences}
}
\description{
convert phylo to treedata
}
\details{
converting phylo object to treedata object
}
\author{
Guangchuang Yu

Yu Guangchuang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mega.R
\name{read.mega_tabular}
\alias{read.mega_tabular}
\title{read.mega_tabular}
\usage{
read.mega_tabular(file)
}
\arguments{
\item{file}{MEGA tabular file}
}
\value{
treedata object
}
\description{
parse tabular output of MEGA
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/isTip.R
\name{isTip}
\alias{isTip}
\alias{isTip.tbl_tree}
\alias{isTip.phylo}
\alias{isTip.treedata}
\title{isTip}
\usage{
isTip(.data, .node, ...)

\method{isTip}{tbl_tree}(.data, .node, ...)

\method{isTip}{phylo}(.data, .node, ...)

\method{isTip}{treedata}(.data, .node, ...)
}
\arguments{
\item{.data}{phylo, treedata or tbl_tree object}

\item{.node}{node number}

\item{...}{additional parameters}
}
\value{
logical value
}
\description{
whether the node is a tip
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-utilities.R
\name{is.ggtree}
\alias{is.ggtree}
\title{is.ggtree}
\usage{
is.ggtree(x)
}
\arguments{
\item{x}{object}
}
\value{
TRUE or FALSE
}
\description{
test whether input object is produced by ggtree function
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write-beast.R
\name{write.beast}
\alias{write.beast}
\title{write.beast}
\usage{
write.beast(treedata, file = "", translate = TRUE, tree.name = "UNTITLED")
}
\arguments{
\item{treedata}{\code{treedata} object}

\item{file}{output file. If file = "", print the output content on screen}

\item{translate}{whether translate taxa labels}

\item{tree.name}{name of the tree}
}
\value{
output file or file content on screen
}
\description{
Export \code{treedata} object to BEAST NEXUS file. This function was adopted and modified from ape::write.nexus
}
\examples{
nhxfile <- system.file("extdata/NHX", "phyldog.nhx", package="treeio")
nhx <- read.nhx(nhxfile)
write.beast(nhx)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/method-get-treetext.R
\docType{methods}
\name{get.treetext}
\alias{get.treetext}
\alias{get.treetext,treedata-method}
\title{get.treetext method}
\usage{
get.treetext(object, ...)

\S4method{get.treetext}{treedata}(object)
}
\arguments{
\item{object}{treedata object}

\item{...}{additional parameter}
}
\value{
phylo object
}
\description{
access tree text (newick text) from tree object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write-beast.R
\name{write.beast.newick}
\alias{write.beast.newick}
\title{write.beast.newick}
\usage{
write.beast.newick(
  treedata,
  file = "",
  append = FALSE,
  digits = 10,
  tree.prefix = ""
)
}
\arguments{
\item{treedata}{\code{treedata} object}

\item{file}{output file. If file = "", print the output content on screen}

\item{append}{logical. Only used if the argument 'file' is the name of file
(and not a connection or "|cmd").  If 'TRUE' output will be appended to 
'file'; otherwise, it will overwrite the contents of file.}

\item{digits}{integer, the indicating the number of decimal places, default is 10.}

\item{tree.prefix, }{character the tree prefix, default is "".}
}
\value{
output file or file content on screen
}
\description{
Export \code{treedata} object to BEAST Newick file. This is useful for making BEAST starting trees with metadata
}
\examples{
nhxfile <- system.file("extdata/NHX", "phyldog.nhx", package="treeio")
nhx <- read.nhx(nhxfile)
write.beast.newick(nhx)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jtree.R
\name{write.jtree}
\alias{write.jtree}
\title{write.jtree}
\usage{
write.jtree(treedata, file = "")
}
\arguments{
\item{treedata}{\code{treedata} object}

\item{file}{output file. If file = "", print the output content on screen}
}
\value{
output file or file content on screen
}
\description{
Export \code{treedata} object to json tree file
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyphy.R
\name{read.hyphy}
\alias{read.hyphy}
\title{read.hyphy}
\usage{
read.hyphy(nwk, ancseq, tip.fasfile = NULL)
}
\arguments{
\item{nwk}{tree file in nwk format, one of hyphy output}

\item{ancseq}{ancestral sequence file in nexus format,
one of hyphy output}

\item{tip.fasfile}{tip sequence file}
}
\value{
A hyphy object
}
\description{
read HYPHY output
}
\examples{
nwk <- system.file("extdata/HYPHY", "labelledtree.tree", package="treeio")
ancseq <- system.file("extdata/HYPHY", "ancseq.nex", package="treeio")
read.hyphy(nwk, ancseq)
}
\author{
Guangchuang Yu \url{https://guangchuangyu.github.io}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codeml_mlc.R
\name{read.codeml_mlc}
\alias{read.codeml_mlc}
\title{read.codeml_mlc}
\usage{
read.codeml_mlc(mlcfile)
}
\arguments{
\item{mlcfile}{mlc file}
}
\value{
A \code{codeml_mlc} object
}
\description{
read mlc file of codeml output
}
\examples{
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="treeio")
read.codeml_mlc(mlcfile)
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ape.R
\name{Nnode.treedata}
\alias{Nnode.treedata}
\title{Nnode}
\usage{
\method{Nnode}{treedata}(phy, internal.only = TRUE, ...)
}
\arguments{
\item{phy}{treedata object}

\item{internal.only}{whether only count internal nodes}

\item{...}{additional parameters}
}
\value{
number of nodes
}
\description{
number of nodes
}
\examples{
Nnode(rtree(30))
}
\author{
Guangchuang Yu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeio-package.R
\docType{package}
\name{treeio-package}
\alias{treeio}
\alias{treeio-package}
\title{treeio: Base Classes and Functions for Phylogenetic Tree Input and Output}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

'treeio' is an R package to make it easier to import and store phylogenetic tree with associated data; and to link external data from different sources to phylogeny. It also supports exporting phylogenetic tree with heterogeneous associated data to a single tree file and can be served as a platform for merging tree with associated data and converting file formats.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/YuLab-SMU/treeio (devel)}
  \item \url{https://docs.ropensci.org/treeio/ (docs)}
  \item \url{https://yulab-smu.top/treedata-book/ (book)}
  \item \url{https://doi.org/10.1093/molbev/msz240 (paper)}
  \item Report bugs at \url{https://github.com/YuLab-SMU/treeio/issues}
}

}
\author{
\strong{Maintainer}: Guangchuang Yu \email{guangchuangyu@gmail.com} (\href{https://orcid.org/0000-0002-6485-8781}{ORCID})

Other contributors:
\itemize{
  \item Tommy Tsan-Yuk Lam \email{tylam.tommy@gmail.com} [contributor, thesis advisor]
  \item Shuangbin Xu \email{xshuangbin@163.com} (\href{https://orcid.org/0000-0003-3513-5362}{ORCID}) [contributor]
  \item Bradley Jones \email{brj1@sfu.ca} [contributor]
  \item Casey Dunn \email{casey_dunn@brown.edu} [contributor]
  \item Tyler Bradley \email{tcb85@drexel.edu} [contributor]
  \item Konstantinos Geles \email{konstantinos.geles@studenti.unicz.it} [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-subset.R
\name{tree_subset}
\alias{tree_subset}
\alias{tree_subset.phylo}
\alias{tree_subset.treedata}
\title{Subset tree objects by related nodes}
\usage{
tree_subset(
  tree,
  node,
  levels_back = 5,
  group_node = TRUE,
  group_name = "group",
  root_edge = TRUE
)

\method{tree_subset}{phylo}(
  tree,
  node,
  levels_back = 5,
  group_node = TRUE,
  group_name = "group",
  root_edge = TRUE
)

\method{tree_subset}{treedata}(
  tree,
  node,
  levels_back = 5,
  group_node = TRUE,
  group_name = "group",
  root_edge = TRUE
)
}
\arguments{
\item{tree}{a tree object of class phylo}

\item{node}{either a tip label or a node number for the given
tree that will be the focus of the subsetted tree}

\item{levels_back}{a number specifying how many nodes back from
the selected node the subsetted tree should include}

\item{group_node}{whether add grouping information of selected node}

\item{group_name}{group name (default 'group') for storing grouping information if group_node = TRUE}

\item{root_edge}{If TRUE (by default), set root.edge to path length of orginal root to the root of subset tree}
}
\description{
This function allows for a tree object to be subset by specifying a
node and returns all related nodes within a selected number of
levels
}
\details{
This function will take a tree and a specified node from
that tree and subset the tree showing all relatives back to a specified
number of nodes. This function allows for a combination of
\code{ancestor} and \code{offspring} to return a subsetted
tree that is of class phylo. This allows for easy graphing of the tree
with \code{ggtree}
}
\examples{
\dontrun{
  nwk <- system.file("extdata", "sample.nwk", package="treeio")
  tree <- read.tree(nwk)

  sub_tree <- tree_subset(tree, node = "A", levels_back = 3)
  ggtree(sub_tree) + geom_tiplab() + geom_nodelab()
}

\dontrun{
  nwk <- system.file("extdata", "sample.nwk", package="treeio")
  tree <- read.tree(nwk)

  sub_tree <- tree_subset(tree, node = "A", levels_back = 3)
  ggtree(sub_tree) + geom_tiplab() + geom_nodelab()
}

}
