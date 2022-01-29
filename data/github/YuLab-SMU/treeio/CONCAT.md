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

