
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our community a harassment-free experience for everyone, regardless of age, body size, visible or invisible disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming, diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes, and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of acceptable behavior and will take appropriate and fair corrective action in response to any behavior that they deem inappropriate, threatening, offensive, or harmful.

Community leaders have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, and will communicate reasons for moderation decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when an individual is officially representing the community in public spaces. Examples of representing our community include using an official e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by
opening an issue or contacting one or more of the project maintainers.

All community leaders are obligated to respect the privacy and security of the reporter of any incident.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 2.0,
available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct enforcement ladder](https://github.com/mozilla/diversity).

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://www.contributor-covenant.org/translations.

# taxa

[![Project Status: WIP - Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/)
[![cran
version](http://www.r-pkg.org/badges/version/taxa)](https://cran.r-project.org/package=taxa)

This is an R package that provides classes to store and manipulate
taxonomic data. Most of the classes can be used like base R vectors.
This project is a partial rewrite of the previous version of `taxa` and
is currently under development.

**A note about recent changes:**

This is beginning of a complete rewrite of the previous `taxa` package
to make the more basic component classes more like base R vectors. The
`taxmap` class is not yet reimplemented, but will be similar to the
class in the previous versions of taxa. The old version of `taxa` has
been incorporated into the `metacoder` package until this version of
taxa is mature, at which time `metacoder` will also use this version.

## Contributors

-   [Zachary Foster](https://github.com/zachary-foster)
-   [Scott Chamberlain](https://github.com/sckott)

## Comments and contributions

We welcome comments, criticisms, and especially contributions! GitHub
issues are the preferred way to report bugs, ask questions, or request
new features. You can submit issues here:

<https://github.com/ropensci/taxa/issues>

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/taxa/issues).
-   License: MIT
-   Get citation information for `taxa` in R doing
    `citation(package = 'taxa')`
-   Please note that this project is released with a Contributor Code of
    Conduct (see CONDUCT.md). By participating in this project you agree
    to abide by its terms.
taxa 0.4.0
==========

LARGE CHANGES:

The beginning of a complete rewrite of the `taxa` package to make the more basic component classes more like base R vectors.
The `taxmap` class is not yet reimplemented, but will be similar to the class in the previous versions of taxa.
The old version of `taxa` has been incorperated into the `metacoder` package unitl this version of taxa is mature, at which time `metacoder` will also use this version.

taxa 0.3.4
==========

### Bug fixes

* Fixed bug in `n_obs` that would cause an error when used on an object with tables with columns named by numbers.
* Various minor bug fixes

taxa 0.3.3
==========

### Bug fixes

* Numeric column names in tables in `taxmap` are now supported
* Various small bug fixes

taxa 0.3.2
==========

### Bug fixes

* Parsers now correclty handle zero-length inputs ([issue #185](https://github.com/ropensci/taxa/issues/185)).
* `taxonomy_table` option `add_id_col` now works ([issue #191](https://github.com/ropensci/taxa/issues/191)).

### Improvements

* The `parse_tax_data` option `class_col` now accepts negative column indexes, meaning "all other columns".

taxa 0.3.1
==========

### New features

* Added the `taxonomy_table` function that converts the information in a `taxmap` or `taxonomy` object into a table with taxa as rows and ranks as columns.
* Added the `print_tree` function that prints text-based trees of `taxmap` or `taxonomy` objects ([issue #173](https://github.com/ropensci/taxa/issues/173)).
* Added `get_dataset` function to get a single data set from `taxmap` objects. Useful for piping with `%>%`.
* `filter_taxa` and `filter_obs` can now subset anything that has names, length, and can be subset, not just tables, lists, and vectors. For example, `DNAbin` objects from the `ape` package can now be used in `taxmap` objects ([issue #178](https://github.com/ropensci/taxa/issues/178)).

### Improvements

* Parsers are somewhat faster and use less RAM ([issue #177](https://github.com/ropensci/taxa/issues/177)).
* `taxmap` and `taxonomy` parsers now treat taxa with the same name and same place in the taxonomy, but different ranks, database IDs, or authorities, as different taxa.
* `filter_obs` can now filter multiple datasets at once if they are the same length ([issue #179](https://github.com/ropensci/taxa/issues/179)).
* `select_obs` and `arrange_obs` can now work on multiple datasets at once.

### Bug fixes

* Made the `"taxon_rank"` value for the `class_key` options work with `extract_tax_data`. 
* Fixed bug in `taxmap` print method when printing tables with only a taxon ID column ([issue #181](https://github.com/ropensci/taxa/issues/181)).

### Changes

* Option `target` in many functions renamed to `data` to make it more intuitive. 

taxa 0.2.1
==========

### Improvements

* `parse_tax_data` can now incorporate rank information which can be accessed by `result$taxon_ranks()` ([issue #113](https://github.com/ropensci/taxa/issues/113)).
* `taxmap` print methods now have more information and color ([issue #124](https://github.com/ropensci/taxa/issues/124)).
* Added `leaves_apply` function that works like `subtaxa_apply`, but on leaves ([issue #126](https://github.com/ropensci/taxa/issues/126)).
* Functions with a `value` option now return named taxon indexes by default, instead of unnamed taxon indexes ([issue #128](https://github.com/ropensci/taxa/issues/128)).
* `lookup_tax_data` and `extract_tax_data` can now use "fuzzy" matching when looking up taxon names, so taxon names can be misspelled and still be founds.
* `lookup_tax_data` and `extract_tax_data` now only look up unique sequence IDs, improving download speed.
* `filter_obs` now can filter out observations in non-target data sets that are associated with taxa that are removed when `drop_taxa = TRUE` ([issue #143](https://github.com/ropensci/taxa/issues/143)). This is done using `filter_taxa`, so the `supertaxa`, `subtaxa`, and `reassign_obs` options are now available to `filter_obs` to control how taxon removal is done.
* `lookup_tax_data` and `extract_tax_data` now have progress bars instead of printing lots of text when downloading information.
* `mutate_obs` now creates new vector/tables if the data set specified does not exist ([issue #121](https://github.com/ropensci/taxa/issues/124)).
* Add `filter_taxa` option `keep_order` that preserves input taxon order. It is `TRUE` by default, which changes how it used to work. Set to `FALSE` for old behavior.
* Using NSE with an ambiguous name (appears in multiple datasets) now produces a warning ([issue #153](https://github.com/ropensci/taxa/issues/153)).

### Changes

* The `simplify` option in many functions is now always handled the same way: If all vectors in a list are names, then unique key-value pairs are returned. Otherwise, names are ignored and unique values are returned.
* The `leaves` option now behaves like `subtaxa`, returning all leaves for each taxon. The old behavior can be replicated by setting the new `simplify` option to `TRUE` ([issue #127](https://github.com/ropensci/taxa/issues/127)).

### Bug fixes

* `filter_taxa` now has better error messages for invalid inputs ([issue #117](https://github.com/ropensci/taxa/issues/117)).
* Fix a bug that caused an error in `filter_taxa` when no taxa pass filter ([issue #116](https://github.com/ropensci/taxa/issues/116)).
* Fixed a bug in `parse_tax_data` when `class_key` was not named ([issue #131](https://github.com/ropensci/taxa/issues/131)).
* Fixed bug in `hierarchy` print method with `taxon_id` class was not used ([issue #138](https://github.com/ropensci/taxa/issues/138)).
* Fixed bug in `parse_tax_data` when all classification data is `NA.`
* Fixed bug in `taxmap` print method when printing zero-length lists and vectors ([issue #148](https://github.com/ropensci/taxa/issues/148)).

taxa 0.2.0
==========

### Bug fixes

* Fixed a few problems with using duplicated inputs to `subset` ([issue #88](https://github.com/ropensci/taxa/issues/85), [issue #89](https://github.com/ropensci/taxa/issues/85))
* Fixed a bug that caused an error when using unnamed vectors ([issue #86](https://github.com/ropensci/taxa/issues/86))
* Fixed a bug that prevents using sequence accession numbers ([issue #85](https://github.com/ropensci/taxa/issues/85))
* Fixed bug in `lookup_tax_data` and `extract_tax_data` that caused an error when one of the queries failed too download.
* Fixed bug that caused "data" argument of `obs_apply` to not work when passed as a variable ([issue #97](https://github.com/ropensci/taxa/issues/97))

### Improvements

* Added `map_data_` for mapping without using NSE.
* Make default dataset for `n_obs` and `n_obs_1` and make them available for NSE ([issue #91](https://github.com/ropensci/taxa/issues/91)
* `parse_tax_data`/`extract_tax_data` can now parse things like `phylum;Nitrosopumilales;order;Nitrosopumilaceae;family;` and split out the rank and taxon names by using multiple matches to the `class_regex` when `class_sep` is NULL. 
* `extract_tax_data` now gives warnings if a regex does not match.
* Added `n_supertaxa_1` function to get number of immediate supertaxa (always 1 or 0).
* Added `branches` function to go with `roots`, `leaves`, and `stems`. ([issue #56](https://github.com/ropensci/taxa/issues/56))
* Added `internodes` and `is_internode` functions to go with `roots`, `leaves`, `branches`, and `stems`. Useful for removing uninformative taxonomic ranks/taxa.
* Started to incorporate ability for `taxon`, `taxon_name`, `taxon_id`, `taxon_rank`, and `taxa` to handle `NULL` inputs as first class citizens to handle cases when you have essentially a blank taxon (use case comes from `taxize` package) [#95](https://github.com/ropensci/taxa/issues/95) [#107](https://github.com/ropensci/taxa/issues/107)
* data parsers: Put long, often unused columns last ([issue #93](https://github.com/ropensci/taxa/issues/93))
* When parsing classifications that have per-taxon info add input id column ([issue #92](https://github.com/ropensci/taxa/issues/92))
* New function `classification` as an abstraction to get either hierarchy of taxon indexes, names, or ids ([issue #57](https://github.com/ropensci/taxa/issues/57))
* New function `get_data_frame` for both `Taxonomy` and `Taxmap` objects that wraps around `get_data` to coerce into a `data.frame`. ([issue #58](https://github.com/ropensci/taxa/issues/58)) ([PR #105](https://github.com/ropensci/taxa/issues/105))

### Changes

* In the output of the taxmap parsing functions like `parse_tax_data`, I moved "taxon_id" and "input_index" columns to front and "input" to rear. Also "tax_data" now comes before "class_data".

taxa 0.1.0
==========

### NEW FEATURES

* Released to CRAN.
## Test environments and check results

### Local computer: Pop!_OS 20.04 LTS, R version 4.0.3

0 errors | 0 warnings | 0 notes

### win-builder (devel and release)

0 errors | 0 warnings | 1 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Zachary Foster <zacharyfoster1989@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  Vectorized (6:42)

### Rhub (Fedora Linux, R-devel and	Ubuntu Linux 20.04.1 LTS, R-release)

0 errors | 0 warnings | 1 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Zachary Foster <zacharyfoster1989@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  Vectorized (6:42)

## Reverse dependencies

### taxlist

The only code used is a function to convert a class used in `taxa` to one used in `taxlist`.
This function will break, but no other functionality in `taxlist` should be affected.
I have alerted the maintainer to this issue 2 weeks ago (Jun 24).
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our community a harassment-free experience for everyone, regardless of age, body size, visible or invisible disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming, diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes, and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of acceptable behavior and will take appropriate and fair corrective action in response to any behavior that they deem inappropriate, threatening, offensive, or harmful.

Community leaders have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, and will communicate reasons for moderation decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when an individual is officially representing the community in public spaces. Examples of representing our community include using an official e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

All community leaders are obligated to respect the privacy and security of the reporter of any incident.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 2.0,
available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct enforcement ladder](https://github.com/mozilla/diversity).

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://www.contributor-covenant.org/translations.

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
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/taxa/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/taxa.git`
* Make sure to track progress upstream (i.e., on our version of `taxa` at `ropensci/taxa`) by doing `git remote add upstream https://github.com/ropensci/taxa.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/taxa`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

```

```
---
output: md_document
---

taxa
====

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/)
[![cran version](http://www.r-pkg.org/badges/version/taxa)](https://cran.r-project.org/package=taxa)

This is an R package that provides classes to store and manipulate taxonomic data.
Most of the classes can be used like base R vectors.
This project is a partial rewrite of the previous version of `taxa` and is currently under development.

**A note about recent changes:**

This is beginning of a complete rewrite of the previous `taxa` package to make the more basic component classes more like base R vectors.
The `taxmap` class is not yet reimplemented, but will be similar to the class in the previous versions of taxa.
The old version of `taxa` has been incorporated into the `metacoder` package until this version of taxa is mature, at which time `metacoder` will also use this version.

## Contributors

* [Zachary Foster](https://github.com/zachary-foster)
* [Scott Chamberlain](https://github.com/sckott)

## Comments and contributions

We welcome comments, criticisms, and especially contributions!
GitHub issues are the preferred way to report bugs, ask questions, or request new features.
You can submit issues here:

https://github.com/ropensci/taxa/issues

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/taxa/issues).
* License: MIT
* Get citation information for `taxa` in R doing `citation(package = 'taxa')`
* Please note that this project is released with a Contributor Code of Conduct (see CONDUCT.md). By participating in this project you agree to abide by its terms.---
title: "taxa: An R package for taxonomic data"
author: Zachary Foster, Scott Chamberlain, and Niklaus Grunwald
output:
  beamer_presentation:
    theme: "default"
    colortheme: "beaver"
    fonttheme: "structurebold"
    fig_caption: false
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
# options(crayon.enabled = TRUE)
options(width = 100)
```

## The challenges of taxonomic data

\begin{itemize}
  \setlength\itemsep{1em}
  \item Taxonomic data is hierarchical
  \item It is often associated with other data
  \item "Taxa" can be names, classifications of names, or IDs
  \item Each source of taxonomic data formats things differently
\end{itemize}

\vspace{5mm}

![](database_format_table.png)

## What taxa provides

\begin{columns}
\begin{column}{0.45\textwidth}

\begin{itemize}
  \setlength\itemsep{1em}
  \item Classes to hold taxa, taxonomies, and associated data
  \item Flexible parsers to convert raw data to these classes
  \item Dplyr-inspired functions to manipulate these classes
  \item A flexible base for other packages to use
\end{itemize}

\end{column}
\begin{column}{0.8\textwidth}  %%<--- here
    \begin{center}
     \includegraphics[width=0.8\textwidth]{class_diagram.png}
     \end{center}
\end{column}
\end{columns}


## `taxmap`: user-defined data mapped to a taxonomy

\begin{figure}
  \makebox[\textwidth][c]{\includegraphics[width=0.9\textwidth]{taxmap_printed.png}}%
\end{figure}

## Reading data from diverse formats

\begin{figure}
  \makebox[\textwidth][c]{\includegraphics[width=1.2\textwidth]{parsing_guide.png}}%
\end{figure}

## `Dplyr`-like manipulation of taxonomic data

Subset taxonomy and data to one taxon:

```{r echo=TRUE, eval=FALSE}
filter_taxa(x, taxon_names == "Plantae", subtaxa = TRUE)
```

Subset taxonomy to one rank:

```{r echo=TRUE, eval=FALSE}
filter_taxa(x, taxon_ranks == "genus", supertaxa = TRUE)
```

Subset data and remove any taxa not in subset:

```{r echo=TRUE, eval=FALSE}
filter_obs(x, "info", n_legs == 4, drop_taxa = TRUE)
```

Add a column to a dataset:

```{r echo=TRUE, eval=FALSE}
mutate_obs(x, "info", bipedal = n_legs == 2)
```

---
title: "Taxa and metacoder: R packages for parsing, visualization, and manipulation of taxonomic data"
author: Zachary Foster, Scott Chamberlain, Thomas Sharpton, and Niklaus Grünwald
header-includes:
   - \usepackage{subfig}
   - \usepackage{graphicx}
output:
  beamer_presentation:
    theme: "default"
    colortheme: "beaver"
    fonttheme: "structurebold"
    fig_caption: false
    keep_tex: yes
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, comment = "", warning = FALSE, cache = TRUE, autodep = TRUE)
# options(crayon.enabled = TRUE)
options(width = 80)
library(taxa)
library(metacoder)
library(cowplot)
```

## The challenges of taxonomic data


\begin{columns}

\begin{column}{0.3\textwidth}
  \includegraphics[height=8cm]{images/taxonomy.png}
\end{column}

\begin{column}{0.8\textwidth}  %%<--- here
  \begin{itemize}
    \setlength\itemsep{1em}
    \item Taxonomic data is hierarchical
    \item Associated with tabular data
    \item Can be names, classifications, or IDs
    \item Many different taxonomic systems
    \item Many different data formats
    \item Hierarchical visualization is difficult 
  \end{itemize}
\end{column}
  
\end{columns}



## Sources of taxonomic data

**DNA sequence databases**

\begin{center}
\resizebox{.99\textwidth}{!}{%
\includegraphics[height=1cm]{images/ncbi_icon.png}%
\quad
\includegraphics[height=1cm]{images/unite_icon.png}%
\quad
\includegraphics[height=1cm]{images/greengenes_icon.png}%
}
\end{center}

**Species occurrence databases**

\begin{center}
\resizebox{.99\textwidth}{!}{%
\includegraphics[height=1cm]{images/gbif_icon.png}%
\quad
\includegraphics[height=1cm]{images/inaturalist_icon.png}%
\quad
\includegraphics[height=1cm]{images/atlas_icon.jpg}%
}
\end{center}

**Museum records**

\begin{center}
\resizebox{.99\textwidth}{!}{%
\includegraphics[height=1cm]{images/smith_icon.jpg}%
\quad
\includegraphics[height=1cm]{images/bmnh_icon.jpg}%
\quad
\includegraphics[height=1cm]{images/nhm_icon.jpg}%
}
\end{center}



## Sources of taxonomic data: DNA sequences {.fragile}

\scriptsize

\textbf{NCBI GenBank}

\begin{verbatim}
>AC073210.8 Homo sapiens BAC clone RP11-460N20 from 7, complete sequence
AACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTC...
\end{verbatim}

\textbf{UNITE}

\begin{verbatim}
>SH099456.05FU_FJ357315_refs k__Fungi;p__Ascomycota;c__Dothideomycetes
;o__Pleosporales;f__Pleosporaceae;g__Embellisia;s__Embellisia_planifunda  
GCTGGCGGCGTGCCTAACACATGTAAGTCGAACGGGACTGGGGGCAACTCCAGT...
\end{verbatim}

\textbf{RDP}

\begin{verbatim}
>S000448483 Sparassis crispa; MBUH-PIRJO&ILKKA94-1587/ss5
Lineage=Root;rootrank;Fungi;domain;Basidiomycota;phylum;Agaricomycetes;
class;Polyporales;order;Sparassidaceae;family;Sparassis;genus
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGAATGCTTAACACATGAAAC...
\end{verbatim}

\textbf{SILVA}

\begin{verbatim}
>GCVF01000431.1.2369 Bacteria;Proteobacteria;Gammaproteobacteria;
Oceanospirillales;Alcanivoraceae;Alcanivorax;Thalassiosira rotula
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGTATGCTTAACACATGCAAG...
\end{verbatim}


## Sources of taxonomic data: DNA sequences {.fragile}

\scriptsize

\textbf{NCBI GenBank}

\begin{Verbatim}[commandchars=\\\{\}]
>\underline{AC073210.8} \underline{Homo sapiens} BAC clone RP11-460N20 from 7, complete sequence
AACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTC...
\end{Verbatim}

\textbf{UNITE}

\begin{Verbatim}[commandchars=\\\{\}]
>SH099456.05FU_FJ357315_refs \underline{k__Fungi;p__Ascomycota;c__Dothideomycetes}
\underline{;o__Pleosporales;f__Pleosporaceae;g__Embellisia;s__Embellisia_planifunda}  
GCTGGCGGCGTGCCTAACACATGTAAGTCGAACGGGACTGGGGGCAACTCCAGT...
\end{Verbatim}

\textbf{RDP}

\begin{Verbatim}[commandchars=\\\{\}]
>S000448483 \underline{Sparassis crispa}; MBUH-PIRJO&ILKKA94-1587/ss5
Lineage=\underline{Root;rootrank;Fungi;domain;Basidiomycota;phylum;Agaricomycetes;}
\underline{class;Polyporales;order;Sparassidaceae;family;Sparassis;genus}
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGAATGCTTAACACATGAAAC...
\end{Verbatim}

\textbf{SILVA}

\begin{Verbatim}[commandchars=\\\{\}]
>GCVF01000431.1.2369 \underline{Bacteria;Proteobacteria;Gammaproteobacteria;}
\underline{Oceanospirillales;Alcanivoraceae;Alcanivorax;Thalassiosira rotula}
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGTATGCTTAACACATGCAAG...
\end{Verbatim}


## Sources of taxonomic data: Occurrence records

**Global Biodiversity Information Facility : Archaea database**

\scriptsize
```{r echo = TRUE, message=FALSE}
readr::read_tsv("datasets/gbif_archea.csv")[4:8]
```


## Sources of taxonomic data: Museum records

**Smithsonian Museum of Natural History: Mammal database**

\scriptsize
```{r echo = TRUE, message=FALSE}
readr::read_csv("datasets/SNMNH.csv")[9]
```


## The `taxa` package

![](images/taxa_header.png)

\vfill

The `taxa` package is designed to be a solid foundation for using taxonomic data in R. 

\vfill

\begin{itemize}
  \setlength\itemsep{1em}
  \item R6 classes to hold taxa, taxonomies, and associated data
  \item Flexible parsers to convert raw data to these classes
  \item Dplyr-inspired functions to manipulate these classes
  \item Functions to get data associated with each taxon in a taxonomy
\end{itemize}


## The `metacoder` package: visualization of taxon data

```{r include=FALSE}
set.seed(6)
```


\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE, warning=FALSE}
readr::read_tsv("datasets/gbif_archea.csv") %>%
  parse_tax_data(class_cols = 4:8) %>%
  filter_taxa(taxon_names != "") %>%
  heat_tree(node_label = taxon_names, node_color = n_obs,
            node_size = n_obs, layout = "da")
```

## Classes

\begin{center}
\Huge
Classes defined by taxa
\end{center}

## Classes defined by `taxa`: Relationships

\begin{center}
\resizebox{.8\textwidth}{!}{%
\includegraphics[height=1cm]{images/class_diagram.png}%
}
\end{center}


## Classes defined by `taxa`: Relationships

\begin{center}
\resizebox{.8\textwidth}{!}{%
\includegraphics[height=1cm]{images/class_diagram_selected.png}%
}
\end{center}


## Classes defined by `taxa`: The `taxmap` class

\begin{center}
\hspace*{-0.6cm}
\resizebox{1.1\textwidth}{!}{%
\includegraphics{images/taxmap_printed.png}%
}
\end{center}


## Parsing
  
\begin{center}
\Huge
Parsing
\end{center}


## Parsing

\begin{center}
\vspace*{-1cm}\hspace*{-1cm}
\resizebox{1.18\textwidth}{!}{%
\includegraphics{images/parsing_guide.png}%
}
\end{center}


## Parsing: vectors of classifications

\scriptsize
```{r echo=TRUE, message=FALSE, cache=TRUE}
x <- c("Mammalia;Theria;Metatheria;Diprotodontia;Macropodiformes",
       "Mammalia;Theria;Eutheria;Primates;Haplorrhini;Simiiformes") 

parse_tax_data(x, class_sep = ";")
```


## Parsing: vectors of names

\scriptsize
```{r echo=TRUE, message=FALSE, cache=TRUE}
x <- c("Homo sapiens", "Macropus", "Chordata")

lookup_tax_data(x, type = "taxon_name", database = "ncbi")
```


## Parsing: vectors of taxon or sequence IDs

\footnotesize

**Taxon IDs**

```{r echo=TRUE, message=FALSE, cache=TRUE, eval=FALSE}
x <- c("9606", "207598", "7711") # NCBI taxon IDs
lookup_tax_data(x, type = "taxon_id", database = "ncbi")
```

**Sequence IDs**

```{r echo=TRUE, message=FALSE, cache=TRUE, eval=FALSE}
x <- c("AC073210", "MG014608", "AE006468") # NCBI sequence IDs
lookup_tax_data(x, type = "seq_id", database = "ncbi")
```


## Parsing: tables

**Global Biodiversity Information Facility : Archaea database**

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
readr::read_tsv("datasets/gbif_archea.csv")[4:8]
```


## Parsing: tables

\scriptsize
```{r echo = TRUE, message=FALSE, warning=FALSE, cache=TRUE}
x = readr::read_tsv("datasets/gbif_archea.csv")

parse_tax_data(x, class_cols = 4:8)
```


## Parsing: complex strings (NCBI Genbank)

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
x = c("AC073210.8 Homo sapiens BAC clone RP11-460N20 from 7, complete sequence",
      "AE006468.2 Salmonella enterica subsp. enterica serovar Typhimurium",
      "MG014608.1 Macropus fuliginosus Csf1r gene, enhancer")

extract_tax_data(x, database = "ncbi", regex = "([A-Z0-9.]+) (.+)",
                 key = c(my_ncbi_id = "seq_id", my_desc = "info"))
```


## Taxon attributes

\begin{center}
\Huge
Taxon attributes
\end{center}


## Taxon attributes: Taxonomy terminology

```{r include=FALSE, cache=TRUE}
hmp_otus$lineage <- sub(hmp_otus$lineage, 
                        pattern = "^r__Root;",
                        replacement = "r__Bacteria;")
hmp_otus$lineage <- paste0("r__Life;", hmp_otus$lineage)
x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                   class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                   class_regex = "^(.+)__(.+)$")
x = filter_taxa(x, taxon_names == "Proteobacteria", subtaxa = TRUE, supertaxa = TRUE, reassign_obs = F)

plot_one <- function(...) {
  set.seed(2)
  x %>%
    heat_tree(
      node_label = ifelse(is_root, "Root", ""),
      node_size = n_obs,
      node_label_size_range = c(0.05, 0.05),
      edge_color = "grey",
      layout = "da",
      initial_layout = "re",
      make_node_legend = FALSE,
      title_size = 0.09,
      ...)
  
}

selected = "red"
unselected = "grey"
target = "black"

my_color <- ifelse(x$taxon_indexes() %in% subtaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T),
                   selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
subtaxa_recursive <- plot_one(title = "Subtaxa (recursive = T)", node_color = my_color)

my_color <- ifelse(x$taxon_indexes() %in% subtaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T, recursive = FALSE), selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
subtaxa_immediate <- plot_one(title = "Subtaxa (recursive = F)", node_color = my_color)

my_color <- ifelse(x$taxon_indexes() %in% supertaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T), selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
supertaxa_recursive <- plot_one(title = "Supertaxa (recursive = T)", node_color = my_color)

my_color <- ifelse(x$taxon_indexes() %in% supertaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T, recursive = FALSE), selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
supertaxa_immediate <- plot_one(title = "Supertaxa (recursive = F)", node_color = my_color)

leaf_plot <- plot_one(title = "Leaves", node_color = ifelse(is_leaf, selected, unselected))

root_plot <- plot_one(title = "Roots", node_color = ifelse(is_root, selected, unselected))

stem_plot <- plot_one(title = "Stems", node_color = ifelse(is_stem, selected, unselected))

internode_plot <- plot_one(title = "Internodes", node_color = ifelse(is_internode, selected, unselected))

branch_plot <- plot_one(title = "Branches", node_color = ifelse(is_branch, selected, unselected))
```

\centering
```{r warning=FALSE, cache=TRUE}
set.seed(2)
x %>%
    heat_tree(
      node_label = taxon_names,
      node_size = n_obs,
      edge_color = "grey",
      node_color = "grey",
      layout = "da",
      initial_layout = "re",
      make_node_legend = FALSE)
```


## Taxon attributes: subtaxa and supertaxa

\centering
```{r warning=FALSE, cache=TRUE}
cowplot::plot_grid(plotlist = list(subtaxa_recursive, subtaxa_immediate, supertaxa_recursive, supertaxa_immediate))
```


## Taxon attributes: parts of a tree

\centering
```{r warning=FALSE, cache=TRUE}
cowplot::plot_grid(plotlist = list(leaf_plot, root_plot, stem_plot, internode_plot))
```


## Taxon attributes: functions

**Ranks, names, and IDs**

`taxon_names`, `taxon_ranks`, `taxon_ids`

**Parts of the tree**

`branches`, `internodes`, `leaves`, `roots`, `stems`, `supertaxa`, `subtaxa`

**Filtering helpers**

`is_branch`, `is_internode`, `is_leaf`, `is_root`, `is_stem`

**Numbers of supertaxa/subtaxa/data**

`n_supertaxa`, `n_subtaxa`, `n_obs`, `n_supertaxa_1`, `n_subtaxa_1`, `n_obs_1` 
 

## Taxon attributes: Ranks, names, and IDs

These are derived from the list of `taxon` objects.

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
taxon_names(ex_taxmap) %>% head
taxon_ranks(ex_taxmap) %>% head
taxon_ids(ex_taxmap) %>% head
```


## Taxon attributes: `subtaxa`

These return a list of vectors named by taxon IDs.

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
subtaxa(ex_taxmap, value = "taxon_names")[1:3]
```


## Taxon attributes: `subtaxa`

These return a list of vectors named by taxon IDs.

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
subtaxa(ex_taxmap, value = "taxon_names", recursive = FALSE)[1:3]
```



## Taxon attributes: counts

These return counts of things per taxon.

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
n_subtaxa(ex_taxmap)
n_supertaxa(ex_taxmap)
n_obs(ex_taxmap, "info")
n_obs(ex_taxmap, "abund")
```


## Manipulating

\begin{center}
\Huge
Manipulating
\end{center}


## Manipulating: example data

Here is the example object that will be used:

```{r include=FALSE}
obj <- ex_taxmap$clone(deep = TRUE)
obj$data$abund <- NULL
```


\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
print(obj)
```


## Manipulating: example data

Here is the example object that will be used:

```{r include=FALSE}
set.seed(1)
```


\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
heat_tree(obj, node_label = taxon_names, layout = "da")
```


## Manipulating: Subsetting the taxonomy

Subset taxonomy and user-defined data to one taxon:

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
filter_taxa(obj, taxon_names == "Plantae", subtaxa = TRUE)
```


## Manipulating: Subsetting the taxonomy

Subset taxonomy and user-defined data to one taxon:

```{r include=FALSE}
set.seed(1)
```

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE, fig.width=3}
filter_taxa(obj, taxon_names == "Plantae", subtaxa = TRUE) %>%
  heat_tree(node_label = taxon_names, layout = "da")
```


## Manipulating: Subsetting the taxonomy

Subset taxonomy to one rank:

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
filter_taxa(obj, taxon_ranks == "family", supertaxa = TRUE)
```


## Manipulating: Subsetting the taxonomy

Subset taxonomy to one rank:

```{r include=FALSE}
set.seed(8)
```

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE, fig.width=3}
filter_taxa(obj, taxon_ranks == "family", supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names, layout = "da")
```


## Manipulating: Subsetting user data

Subset user-defined data and remove any taxa not in subset:

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
filter_obs(obj, "info", n_legs == 4, drop_taxa = TRUE)
```


## Manipulating: Subsetting user data

Subset data and remove any taxa not in subset:

```{r include=FALSE}
set.seed(7)
```

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE, fig.width=2.5}
filter_obs(obj, "info", n_legs == 4, drop_taxa = TRUE) %>%
  heat_tree(node_label = taxon_names, layout = "da")
```


## Manipulating: Adding user data

Add a column to a dataset:

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
mutate_obs(obj, "info", bipedal = n_legs == 2)
```


## Acknowledgements

\begin{center}
\resizebox{.99\textwidth}{!}{%
\includegraphics[height=0.5cm]{images/ropensci_icon.png}%
\quad
\includegraphics[height=0.5cm]{images/r_icon.jpg}%
\quad
\includegraphics[height=0.5cm]{images/user_icon.png}%
}
\end{center}


\begin{center}
\resizebox{.99\textwidth}{!}{%
\includegraphics[height=0.5cm]{images/osu_icon.png}%
\quad
\includegraphics[height=0.5cm]{images/ars_icon.png}%
}
\end{center}
  
\begin{center}
\includegraphics[height=7cm]{images/grunwaldlab.jpg}%
\end{center}


## Manipulating: values accessible to NSE

The following can be used in manipulation functions as if they were independent variables using Non-Standard Evaluation (NSE):

\begin{itemize}
  \setlength\itemsep{1em}
  \item Functions that return per-taxon information
  \item User-defined table columns
  \item User-defined vectors and lists
  \item User-defined functions
\end{itemize}

\scriptsize
```{r echo = TRUE, message=FALSE, cache=TRUE}
unname(all_names(obj))
```


## Metacoder: visualization of taxonomic data

\begin{center}
\vspace*{-0.3cm}\hspace*{-1cm}
\resizebox{0.9\textwidth}{!}{%
\includegraphics{images/heat_tree_matrix.png}%
}
\end{center}


## Metacoder: visualization of taxonomic data

\begin{center}
\vspace*{-0.7cm}\hspace*{-0.6cm}
\resizebox{1.15\textwidth}{!}{%
\includegraphics{images/metacoder_multiroot.png}%
}
\end{center}


## Thanks for listening!


\begin{center}
\Huge
Questions?
\end{center}


# Talk outline

## Introduction (0.5 min)


## Taxonomic information in the wild (1 min)

* sequence databases
* BLAST results
* Sepcies occurance databases
* museum records


## The purpose of the taxa package (0.5 min)


## The classes (3 min total)

### The basic classes (1 min)

### Taxonomy and Taxmap (2 min)


## Taxon attributes (4.5 min total)

### Ranks, names, and IDs (0.5 min)

* taxon_names, taxon_ranks, taxon_ids

### Parts of the tree (3 min)

* branches, internodes, leaves, roots, stems, supertaxa, subtaxa

### Numbers of supertaxa/subtaxa/data (1 min)

* n_supertaxa, n_subtaxa, n_obs, n_supertaxa_1, n_subtaxa_1, n_obs_1 


## Parsing (5 min total)

### Vectors (1 min)

### Tables (2 min)

### Complex strings (2 min)


## Manipulation (4 min total)

### Subsetting the taxonomy (2 min)

### Subsetting user data (1.5 min)

### Adding user data (0.5 min)


## Compatibility with other packages (3 min total)

### Metacoder (2.5 min)

### Taxize (0.5 min)



## Conclusion


## Genbank FASTA

```
>AC073210.8 Homo sapiens BAC clone RP11-460N20 from 7, complete sequence
GAATTCTTCAGGTAGCTTCCTAGGGTTTCCAAGGCAATACAAGAAGAATTTTGATAGGCAGGAAAATGCA...
```

## UNITE QIIME

```
SH099456.05FU_FJ357315_refs    k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__Pleosporaceae;g__Embellisia;s__Embellisia_planifunda
```

## Greengenes

```
1111875	k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Limnohabitans; s__
```

## RDP 

```
>S000448483 Sparassis crispa; MBUH-PIRJO&ILKKA94-1587/ss5	Lineage=Root;rootrank;Fungi;domain;Basidiomycota;phylum;Agaricomycetes;class;Polyporales;order;Sparassidaceae;family;Sparassis;genus
```

## SILVA

```
>GCVF01000431.1.2369 Bacteria;Proteobacteria;Gammaproteobacteria;Oceanospirillales;Alcanivoraceae;Alcanivorax;Thalassiosira rotula
CGUGCACGGUGGAUGCCUUGGCAGCCAGAGGCGAUGAAGGACGUUGUAGCCUGCGAUAAGCUCCGGUUAGGUGGCAAACA
```

## Smithsonian Museum of Natural History: Mammals

```{r}
smith_data <- readr::read_csv("datasets/SNMNH.csv")
```

## GBIF: archea

```{r}
gbif_data <- readr::read_tsv("datasets/gbif_archea.csv")
```


## Atlas of Living Australia

```
```# Issue 009) Adapt `taxmap` utility functions to `taxonomy`

## Thu Sep 22 14:18:34 2016

I noticed that the edge list in `taxonomy` does not have an edge for each taxon. 
This is different than how it was structured in `metacoder` and this is making adapting functions like `supertaxa` more difficult than I thought. 
Also, the algolrythm used to construct the edge list involves an all-vs-all comparison of `taxon` objects. 
I think I had originally done something like this in `metacoder`, but eventully opted for a recursive strategy to make it faster on large datasets.
However, I dont think I have ever taken into account situations like:

a > b > c
b > c

Something to look into...


## Wed Dec 7 18:25:57 2016

I think I have a non-recursive way to make an edge list from hierarchies. 

```
mammalia <- taxon(
  name = taxon_name("Mammalia"),
  rank = taxon_rank("class"),
  id = taxon_id(9681)
)
plantae <- taxon(
  name = taxon_name("Plantae"),
  rank = taxon_rank("kingdom"),
  id = taxon_id(33090)
)
unidentified <- taxon(
  name = taxon_name("unidentified"),
  rank = taxon_rank("species"),
  id = taxon_id(0)
)
x <- taxonomy(unidentified_plant, unidentified_animal)
> x
<Taxonomy>
  3 taxa: 1. Plantae, 2. unidentified, 3. Mammalia
  2 edges: 1->2, 3->2
```

# Notes on taxize integration

```{r}
library('taxize')
```


## General thoughts

* Add a `simplify` option for all functions that whose focus is to return information that could be encoded in `taxon_name`, `taxon_id`, `taxon`, or `taxon_rank` vectors, but all return other associated information in tables. If `simplify = TRUE` the just the vector is returned, if `simplify = FALSE` then the vector is the first column in a table containing the other information.


## gnr_resolve

```{r}
temp <- gnr_resolve(names = c("Helianthos annus", "Homo saapiens"))
head(temp)
```

The `matched_name` could be a `taxon_name` or `taxon` vector and the `data_source_title` could be encoded in its `taxon_id` field.


## get_* functions

```{r}
mynames <- c("Helianthus annuus ssp. jaegeri", "Helianthus annuus ssp. lenticularis", "Helianthus annuus ssp. texanus")
(tsn <- get_tsn(mynames, accepted = FALSE))
```

These could return `taxon_id` vectors, which would store the database in their `taxon_db` fields instead of using a class for the database.
The other attributes could remain the same or be put returned as columns in a table.

  

## accepted names functions

```{r}
itis_acceptname(searchtsn = tsn)
```

The `acceptedname`, `acceptedtsn`, and `author` columns could be converted to a `taxon_name` vector column.
Alternativly, a vector of `taxon` objects could be returned with the accepted as the first synonym and the query as the second.

## classification

```{r}
specieslist <- c("Abies procera","Pinus contorta")
classification(specieslist, db = 'itis')
```

This could return a `classification` vector.


## as.*id functions

Would probably be replaced by the `taxon_id` constructor to replicate the `check = FALSE` behavior.
For the `check = TRUE` behavior, perhaps `taxon_id` class should have a `validated` boolean field and then these functions could be replaced with a `check_taxon_id` function that populates that field.
This would allow a vector of IDs from different databases to be validated.


## downstream


```{r}
downstream("Dasypodidae", downto = "Species", db = "col")
downstream("Dasypodidae", downto = "Species", db = "col", intermediate = TRUE)
```

When `intermediate = FALSE`, the return value could be a list of `taxon_name` vectors or the columns with ids, names, ranks, and name status could be replaved with a `taxon_name` vector.

When `intermediate = TRUE`, the same as `intermediate = FALSE` as described above could be used, or a list of `taxonomy` could be returned.


## children

```{r}
children("Pinaceae", db = "ncbi")
```

Can also be list of `taxon_name` or tables with `taxon_name` as a column.


## genbank2uid


```{r}
genbank2uid(id = 'AJ748748')
```

Can return `taxon_id` vector.


## apgOrders and apgFamilies

```{r}
x = apgOrders()
head(x)
```

Can return a `taxon` vector so synonyms can be included.


## apg_lookup

```{r}
apg_lookup(taxa = "Hyacinthaceae", rank = "family")
apg_lookup(taxa = "Asteraceae", rank = "family")
```


Can return a `taxon_name` with the the db set to `apg`


## as.data.frame.*id

```{r}
as.data.frame(as.boldid(c("1973","101009","98597"), check=FALSE))
```

This functionality could mostly be handled by `as.data.frame.taxon_id`, except for the `multiple_matches` `pattern_match` columns.
Perhaps these functions can call  `as.data.frame.taxon_id` and add on those columns


## bold_search

```{r}
bold_search(name=c("Apis","Puma concolor"))
bold_search(id=c(88899, 88899), dataTypes="all")
bold_search(id=c(88899, 88899), dataTypes="stats")
bold_search(id=c(88899, 88899), dataTypes="geo")
```

Can replace `taxid`, `taxon`, and  `tax_rank` columns with a `taxon_name` column


## class2tree

Perhaps this could be modified to take `taxonomy` class as input, since that could be the output of `classifications`.
If so, it might make more sense to be in the `taxa` package, or another focused on visluizing taxonomic info.


## comm2sci

```{r}
comm2sci(commnames='american black bear')
comm2sci(commnames='american black bear')
comm2sci(commnames='black bear', db='itis', simplify = FALSE)
```

Could return a list of `taxon_name` vectors


## eol_search

```{r}
x <- eol_search(terms='Homo')
x
```

pageid, name, link columns could be replaced by `taxon_name` vector.


## eubon_search

```{r}
eubon_search("Prionus")
```

A lot of these columns could be made into a `taxon_name` vector or a `taxon` vector if we want the synonyms to be incorperated. 


## fg_name_search

```{r}
fg_name_search(q = "Gymnopus")
```

A lot of these columns could be made into a `taxon_name` vector or a `taxon` vector if we want the synonyms to be incorperated. 


## id2name

```{r}
id2name(19322, db = "itis")
```

Can return `taxon_name` vector


```{r echo = FALSE}
options(max.print = 100, width = 100)
```


```{r eval = FALSE}
# devtools::install_github("ropensci/taxa", ref = "vectorize")
```

```{r}
library(taxa)
```


I used the `vctrs` package to try out making S3 vector-like versions of classes to hold taxon databases and taxon id.
If this works out, I would like to do the same for taxon names, ranks, and taxa as a whole, which would include the ids, names, and ranks.
The ids, names, and ranks each can hold their own database ID.

## taxon_db 

The `taxon_db` class is just a character vector of database names with extra validity checks to make sure that names are valid.

```{r}
taxon_db()
taxon_db('ncbi')
taxon_db(rep('itis', 100))
taxon_db(rep('itis', 2000))
```

The acceptable databases are stored in `database_definitions` in the same way knitr chunk options are stored and changed.

```{r}
database_definitions$get()
```

All names in `taxon_db` objects must be in this list or errors are thrown:

```{r error=TRUE}
taxon_db(c('ncbi', 'itis', 'custom', 'tps'))
```

This also happens on assigning values:

```{r error=TRUE}
x = taxon_db(rep('itis', 10))
x
x[2:3] <- "ncbi"
x
x[2:3] <- "custom"
```

However, you can add new databases to the list:

```{r}
database_definitions$set(
  name = "custom",
  url = "http://www.my_tax_database.com",
  description = "I just made this up",
  id_regex = ".*"
)
```

And then you can use that name too:

```{r}
x[2:3] <- "custom"
x
```


## taxon_id

The `taxon_id` class stores taxon IDs and their associated database (as a `taxon_db` vector).

```{r}
taxon_id()
taxon_id(1:10)
taxon_id(1:10, db = 'ncbi')
taxon_id(1:10000, db = 'itis')
```

> Run code in an R console to see colored output. I forgot how to make it work in Rmarkdown.

Since a `taxon_db` object is used, database names must be valid:

```{r error=TRUE}
taxon_id(1:10, db = 'custom2')
```

The database values can be accessed and set with `taxon_db`:

```{r}
x <- taxon_id(1:10, db = 'ncbi')
x
taxon_db(x)
taxon_db(x) <- 'itis'
x
taxon_db(x)[1:3] <- 'tps'
x
```

... if they are valid:

```{r error=TRUE}
taxon_db(x)[1:3] <- 'custom2'
```

If a database is defined, then the ids must match the ID regex for that database in `database_definitions$get()`:

```{r error=TRUE}
taxon_id(letters[1:3], c('ncbi', 'itis', 'itis'))
taxon_id(letters[1:3])
```

You can set taxon id values with numbers or characters, but only `taxon_id` inputs can supply a defined database.

```{r error=TRUE}
x <- taxon_id(1:10, db = 'ncbi')
x
x[3:4] <- 1
x[6:7] <- c("a", "b")
x[8:9] <- taxon_id(1, db = c('itis', 'tps'))
x
```


## Use in tables

The `taxon_id` prints well in `data.frame`s and `tibble`s, using color when possible:

```{r}
data.frame(x = taxon_id(1:10, "ncbi"), y = taxon_db(rep("itis", 10)), z = 1:10)
tibble::tibble(x = taxon_id(1:10, "ncbi"), y = taxon_db(rep("itis", 10)), z = 1:10)
str(taxon_id(1:10, 'ncbi'))
```


## Plans for taxon name and rank

`taxon_name` and `taxon_rank` would be very similar to `taxon_id`.
`taxon_rank` would either be based on a factor, or be factor-like and would have database-based validity checks for rank names and order.


## Plans for taxon

`taxon` would contain at least `taxon_name`, `taxon_rank` and `taxon_id` vectors.
It would probably also have authority information.
It could also have field for arbitrary information, but this might not be needed since it could be in a table that could contain that info in other columns.


## Plans for taxonomy, taxmap, and hierarcy

Probably leave as R6 for now.






```{r}
library(metacoder)
hmp_otus$lineage <- sub(hmp_otus$lineage, 
                        pattern = "^r__Root;",
                        replacement = "r__Bacteria;")
hmp_otus$lineage <- paste0("r__Life;", hmp_otus$lineage)
x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                   class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                   class_regex = "^(.+)__(.+)$")
x = filter_taxa(x, taxon_names == "Proteobacteria", subtaxa = TRUE, supertaxa = TRUE, reassign_obs = F)
```


```{r}
plot_one <- function(...) {
  set.seed(2)
  x %>%
    heat_tree(
      node_label = ifelse(is_root, "Root", ""),
      node_size = n_obs,
      # edge_color = "grey",
      layout = "da",
      initial_layout = "re",
      make_node_legend = FALSE,
      title_size = 0.07,
      ...)
  
}
```


* subtaxa, recursive
* subtaxa, immediate
* supertaxa, recusrive
* supertaxa, immediate
* leafs
* roots
* stems
* internodes
* branches

```{r}
selected = "red"
unselected = "grey"
target = "black"

my_color <- ifelse(x$taxon_indexes() %in% subtaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T),
                   selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
subtaxa_recursive <- plot_one(title = "Subtaxa (recursive = T)", node_color = my_color)

my_color <- ifelse(x$taxon_indexes() %in% subtaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T, recursive = FALSE), selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
subtaxa_immediate <- plot_one(title = "Subtaxa (recursive = F) ", node_color = my_color)

my_color <- ifelse(x$taxon_indexes() %in% supertaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T), selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
supertaxa_recursive <- plot_one(title = "Supertaxa (recursive = T)", node_color = my_color)

my_color <- ifelse(x$taxon_indexes() %in% supertaxa(x, subset = taxon_names == "Betaproteobacteria", simplify = T, recursive = FALSE), selected, unselected)
my_color[x$taxon_names() == "Betaproteobacteria"] <- target
supertaxa_immediate <- plot_one(title = "Supertaxa (recursive = F)", node_color = my_color)

leaf_plot <- plot_one(title = "Leaves", node_color = ifelse(is_leaf, selected, unselected))

root_plot <- plot_one(title = "Roots", node_color = ifelse(is_root, selected, unselected))

stem_plot <- plot_one(title = "Stems", node_color = ifelse(is_stem, selected, unselected))

internode_plot <- plot_one(title = "Internodes", node_color = ifelse(is_internode, selected, unselected))

branch_plot <- plot_one(title = "Branches", node_color = ifelse(is_branch, selected, unselected))
```


```{r}
library(cowplot)
cowplot::plot_grid(plotlist = list(subtaxa_recursive, subtaxa_immediate, supertaxa_recursive, supertaxa_immediate,
                                   leaf_plot, root_plot, stem_plot, internode_plot, branch_plot))
```





# Adding getters and setters

Make consistnet getters/setters for all classes with the following characteristics:

* The source values are private and prefixed with `my_`
* The active bindings return objects, but check the object type when setting and can accept vectors when setting. They basically pretend to be the source values.
* objects returned by active bindings can be converted to characters by `as.character`
* Remove the ability to store vectors instead of objects. Vectors given to active bindings are converted to objects automatically.
* Use `taxon_names`, `taxon_ranks`, and `taxon_ids` S3 functions for all classes as vector-based getter/setters. Make them usable as setters( e.g. `taxon_names<-`). It might be a bit odd to use `taxon_names` to get/set of `taxon` objects, but having only one function name for everything will make it easy to use.
* Make `taxonomy`/`taxmap` and `hierarcy` inherit `taxa`, so its getters/setters work for them. 
* S3 setters pass by value, R6 setters pass by reference
* All vailidity checks are in the active binding setters when possible. Any ones that rely on combinations of arguments are in the class initialization after the active binding checks

# Goals of the taxa package

## Encode taxonomic information from divserse sources

Taxonomic data can come from many different sources including:

* Sequence data annotations
* Occurance data annotations
* Taxonomy databases
* Museum records
* Species checklists
* Herbarium records

Each of these might encode taxonomic data in different ways.
Sometimes the taxonomic data itself is the focus and sometimes it is only one of many ways to annotate the data.
Therefore `taxa` needs to be intuitve enough so that users not interested in the intirciaes of taxonomy can ignore them, while sill allowing for more adavnce usages when needed.

## Combine taxonomic data from multiple sources

Combining taxonomic data is difficult since there any many conflicting taxonomies, so `taxa` needs to support conflicting taxonomies.
Often users will have conflicting taxonomic data when combining sources, but they will want to combine it in such a way that a consensus taxonomy can be reached.
Therefore, tools will be needed to identify and resolve conflicts

## Filtering taxonomic data

Hierarchical data is difficult to filter, so `taxa` needs special filtering functions.
These need to operate fast enough to work on datasets with 10,000 taxa or more.


## Manipulating taxonomic data

There should be tools to modify taxonomies, including:

* Add taxa
* Remove taxa
* Combine taxa as synonyms
* Spilt taxa with synonyms into different taxa
* Change/add ranks, names, authorities, etc


## Exploring taxonomic data

The `taxa` pacakge will need ways to explore taxonomic hierarchies interactivly.
Plotting taxonomic data is complex so is best left to other packages, but text-based tree exploration should be supported. 
Combined with effective filtering techniques, this should be sufficient.
---
output:
  pdf_document: default
---

```{r, warning = FALSE, fig.width=12, cache=TRUE}
# Look up plant occurrence data for Oregon
library(rgbif)
occ <- rgbif::occ_data(stateProvince = "Oregon",
                       scientificName = "Plantae")

# Parse data with taxa
library(taxa)
obj <- parse_tax_data(occ$data,
                      class_cols = c(22:26, 28),
                      named_by_rank = TRUE)

# Plot number of occurrences for each taxon
library(metacoder)
obj %>% 
  filter_obs("tax_data",
             basisOfRecord == "PRESERVED_SPECIMEN",
             drop_taxa = TRUE) %>%
  filter_taxa(taxon_ranks != "specificEpithet") %>%
  filter_taxa(! is.na(taxon_names)) %>%
  filter_taxa(taxon_names == "Tracheophyta",
              subtaxa = TRUE) %>%
  filter_taxa(taxon_ranks == "order", n_subtaxa > 10,
              subtaxa = TRUE, supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_color = n_obs,
            node_size = n_obs,
            node_color_axis_label = "# occurrences",
            output_file = c("example_analysis_figure.eps",
                            "example_analysis_figure.pdf",
                            "example_analysis_figure.tiff"))
```

---
title: "Introduction to the taxa package"
author: "Scott Chamberlain and Zachary Foster"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
  toc: true
  vignette: >
    %\VignetteIndexEntry{Introduction to the taxa package}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
comment = "#>",
collapse = TRUE,
warning = FALSE,
message = FALSE,
eval = FALSE
)
```

`taxa` defines taxonomic classes and functions to manipulate them. The goal is to use these classes as low
level fundamental taxonomic classes that other R packages can build on and supply robust manipulation functions (e.g. subsetting) that are broadly useful.

There are two distinct types of classes in `taxa`:

* Classes that are concerned only with taxonomic information: `taxon`, `taxonomy`, `hierarchy`, etc.
* A class called `taxmap` that is concerned with combining taxonomic data with
user-defined data of any type (e.g. molecular sequences, abundance counts etc.)

Diagram of class concepts for `taxa` classes:

```{r results='asis', echo = FALSE}
flowchart_path <- "class_diagram.png"
width <- 718
if (knitr:::child_mode()) { # if run as a child
  flowchart_path <- file.path("vignettes", flowchart_path)
}
cat(paste0('<img src="', flowchart_path, '" title="taxa classes diagram" width="',
       width, '">'))
```

Relationship between classes implemented in the taxa package. Diamond-tipped arrows indicate that objects of one class are used in another class. For example, a database object can stored in the taxon_rank, taxon_name, or taxon_id objects. A standard arrow indicates inheritance. For example, the taxmap class inherits the taxonomy class. `*` means that the object (e.g. a database object) can be replaced by a simple character vector. `?` means that the data is optional (Note: being able to replace objects with characters might be going away soon).

## Install

For the latest "stable" release, use the CRAN version:

```{r eval=FALSE}
install.packages("taxa")
```

For all the latest improvements, bug fixes, and bugs, you can download the development version:

```{r eval=FALSE}
devtools::install_github("ropensci/taxa")
```

```{r}
library("taxa")
```

## The classes

### Minor component classes

There are a few optional classes used to store information in other classes.
These will probably mostly be of interest to developers rather than users.

#### database

Taxonomic data usually comes from a database.
A common example is the [NCBI Taxonomy Database](https://www.ncbi.nlm.nih.gov/taxonomy) used to provide taxonomic classifications to sequences deposited in [other NCBI databases](https://www.ncbi.nlm.nih.gov/guide/all/).
The `database` class stores the name of the database and associated information:

```{r}
(ncbi <- taxon_database(
  name = "ncbi",
  url = "http://www.ncbi.nlm.nih.gov/taxonomy",
  description = "NCBI Taxonomy Database",
  id_regex = "*"
))
ncbi$name
ncbi$url
```

To save on memory, a selection of common databases is provided with the package (`database_list`) and any in this list can be used by name instead of making a new database object (e.g. `"ncbi"` instead of the `ncbi` above).

```{r}
database_list
```


#### rank

Taxa might have defined ranks (e.g. species, family, etc.), ambiguous ranks (e.g. "unranked", "unknown"), or no rank information at all.
The particular selection and format of valid ranks varies with database, so the database can be optionally defined.
If no database is defined, any ranks in any order are allowed.

```{r}
taxon_rank(name = "species", database = "ncbi")
```

#### `taxon_name`

The taxon name can be defined in the same way as rank.

```{r}
taxon_name("Poa", database = "ncbi")
```

#### taxon_id

Each database has its set of unique taxon IDs.
These IDs are better than using the taxon name directly because they are guaranteed to be unique, whereas there are often duplicates of taxon names (e.g. *Orestias elegans* is the name of both an orchid and a fish).

```{r}
taxon_id(12345, database = "ncbi")
```

### The "taxon" class

The `taxon` class combines the classes containing the name, rank, and ID for the taxon.
There is also a place to define an authority of the taxon.

```{r}
(x <- taxon(
  name = taxon_name("Poa annua"),
  rank = taxon_rank("species"),
  id = taxon_id(93036),
  authority = "Linnaeus"
))
```

Instead of the name, rank, and ID classes, simple character vectors can be supplied.
These will be converted to objects automatically.

```{r}
(x <- taxon(
  name = "Poa annua",
  rank = "species",
  id = 93036,
  authority = "Linnaeus"
))
```

The `taxa` class is just a list of `taxon` classes.
It is meant to store an arbitrary list of `taxon` objects.

```{r}
grass <- taxon(
  name = taxon_name("Poa annua"),
  rank = taxon_rank("species"),
  id = taxon_id(93036)
)
mammalia <- taxon(
  name = taxon_name("Mammalia"),
  rank = taxon_rank("class"),
  id = taxon_id(9681)
)
plantae <- taxon(
  name = taxon_name("Plantae"),
  rank = taxon_rank("kingdom"),
  id = taxon_id(33090)
)

taxa(grass, mammalia, plantae)
```

### The "hierarchy" class

[Taxonomic classifications](https://en.wikipedia.org/wiki/Taxonomy_(biology)#Classifying_organisms) are an ordered set of taxa, each at a different rank.
The `hierarchy` class stores a list of `taxon` classes like `taxa`, but `hierarchy` is meant to store all of the taxa in a classification in the correct order.

```{r}
x <- taxon(
  name = taxon_name("Poaceae"),
  rank = taxon_rank("family"),
  id = taxon_id(4479)
)

y <- taxon(
  name = taxon_name("Poa"),
  rank = taxon_rank("genus"),
  id = taxon_id(4544)
)

z <- taxon(
  name = taxon_name("Poa annua"),
  rank = taxon_rank("species"),
  id = taxon_id(93036)
)

(hier1 <- hierarchy(z, y, x))
```

Multiple `hierarchy` classes are stored in the `hierarchies` class, similar to how multiple `taxon` are stored in `taxa`.

```{r}
a <- taxon(
  name = taxon_name("Felidae"),
  rank = taxon_rank("family"),
  id = taxon_id(9681)
)
b <- taxon(
  name = taxon_name("Puma"),
  rank = taxon_rank("genus"),
  id = taxon_id(146712)
)
c <- taxon(
  name = taxon_name("Puma concolor"),
  rank = taxon_rank("species"),
  id = taxon_id(9696)
)
(hier2 <- hierarchy(c, b, a))
```

```{r}
hierarchies(hier1, hier2)
```

### The "taxonomy" class

The `taxonomy` class stores unique `taxon` objects in a tree structure.
Usually this kind of complex information would be the output of a file parsing function, but the code below shows how to construct a `taxonomy` object from scratch (you would not normally do this).

```{r}
# define taxa
notoryctidae <- taxon(name = "Notoryctidae", rank = "family", id = 4479)
notoryctes <- taxon(name = "Notoryctes", rank = "genus", id = 4544)
typhlops <- taxon(name = "typhlops", rank = "species", id = 93036)
mammalia <- taxon(name = "Mammalia", rank = "class", id = 9681)
felidae <- taxon(name = "Felidae", rank = "family", id = 9681)
felis <- taxon(name = "Felis", rank = "genus", id = 9682)
catus <- taxon(name = "catus", rank = "species", id = 9685)
panthera <- taxon(name = "Panthera", rank = "genus", id = 146712)
tigris <- taxon(name = "tigris", rank = "species", id = 9696)
plantae <- taxon(name = "Plantae", rank = "kingdom", id = 33090)
solanaceae <- taxon(name = "Solanaceae", rank = "family", id = 4070)
solanum <- taxon(name = "Solanum", rank = "genus", id = 4107)
lycopersicum <- taxon(name = "lycopersicum", rank = "species", id = 49274)
tuberosum <- taxon(name = "tuberosum", rank = "species", id = 4113)
homo <- taxon(name = "homo", rank = "genus", id = 9605)
sapiens <- taxon(name = "sapiens", rank = "species", id = 9606)
hominidae <- taxon(name = "Hominidae", rank = "family", id = 9604)

# define hierarchies
tiger <- hierarchy(mammalia, felidae, panthera, tigris)
cat <- hierarchy(mammalia, felidae, felis, catus)
human <- hierarchy(mammalia, hominidae, homo, sapiens)
mole <- hierarchy(mammalia, notoryctidae, notoryctes, typhlops)
tomato <- hierarchy(plantae, solanaceae, solanum, lycopersicum)
potato <- hierarchy(plantae, solanaceae, solanum, tuberosum)

# make taxonomy
(tax <- taxonomy(tiger, cat, human, tomato, potato))
```

Unlike the `hierarchies` class, each unique `taxon` object is only represented once in the `taxonomy` object.
Each taxon has a corresponding entry in an [edge list](https://en.wikipedia.org/wiki/Adjacency_list) that encode how it is related to other taxa.
This makes `taxonomy` more compact, but harder to manipulate using standard indexing.
To make manipulation easier, there are functions like `filter_taxa` and `subtaxa` that will be covered later.
In general, the `taxonomy` and `taxmap` objects (covered later) would be instantiated using a parser like `parse_tax_data`. 
This is covered in detail in the parsing vignette.

#### supertaxa

A "supertaxon" is a taxon of a coarser rank that encompasses the taxon of interest (e.g. "Homo" is a supertaxon of "sapiens").
The `supertaxa` function returns the supertaxa of all or a subset of the taxa in a `taxonomy` object.

```{r}
supertaxa(tax)
```

By default, the taxon IDs for the supertaxa of all taxa are returned in the same order they appear in the edge list.
Taxon IDs (character) or edge list indexes (integer) can be supplied to the `subset` option to only return information for some taxa.

```{r}
supertaxa(tax, subset = "m")
```

What is returned can be modified with the `value` option:

```{r}
supertaxa(tax, subset = "m", value = "taxon_names")
```

```{r}
supertaxa(tax, subset = "m", value = "taxon_ranks")
```

You can also subset based on a logical test:

```{r}
supertaxa(tax, subset = taxon_ranks == "genus", value = "taxon_ranks")
```


The `subset` and `value` work the same for most of the following functions as well.
See `all_names(tax)` for what can be used with `value` and `subset`.
Note how `value` takes a character vector (`"taxon_ranks"`), but `subset` can use the same value (`taxon_ranks`) as a part of an expression.
`taxon_ranks` is actually a function that is run automatically when its name is used this way:

```{r}
taxon_ranks(tax)
```

This is an example of [Non-standard evaluation](http://adv-r.had.co.nz/Computing-on-the-language.html) (NSE).
NSE makes codes easier to read an write.
The call to `supertaxa` could also have been written without NSE like so:

```{r}
supertaxa(tax, subset = taxon_ranks(tax) == "genus", value = "taxon_ranks")
```


#### subtaxa

The "subtaxa" of a taxon are all those of a finer rank encompassed by that taxon.
For example, *sapiens* is a subtaxon of *Homo*.
The `subtaxa` function returns all subtaxa for each taxon in a `taxonomy` object.

```{r}
subtaxa(tax, value = "taxon_names")
```

This and the following functions behaves much like `supertaxa`, so we will not go into the same details here.

#### roots

We call taxa that have no supertaxa "roots".
The `roots` function returns these taxa.

```{r}
roots(tax, value = "taxon_names")
```


#### leaves

We call taxa without any subtaxa "leaves".
The `leaves` function returns these taxa.

```{r}
leaves(tax, value = "taxon_names")
```

#### other functions

There are many other functions to interact with `taxonomy` object, such as `stems` and `n_subtaxa`, but these will not be described here for now.

### The "taxmap" class

The `taxmap` class is used to store any number of tables, lists, or vectors associated with taxa.
It is basically the same as the `taxonomy` class, but with the following additions:

* A list called `data` that stores arbitrary user data associated with taxa
* A list called `funcs` that stores user defined functions

All the functions described above for the `taxonomy` class can be used with the `taxmap` class.

```{r}
info <- data.frame(name = c("tiger", "cat", "mole", "human", "tomato", "potato"),
                   n_legs = c(4, 4, 4, 2, 0, 0),
                   dangerous = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))

phylopic_ids <- c("e148eabb-f138-43c6-b1e4-5cda2180485a",
                  "12899ba0-9923-4feb-a7f9-758c3c7d5e13",
                  "11b783d5-af1c-4f4e-8ab5-a51470652b47",
                  "9fae30cd-fb59-4a81-a39c-e1826a35f612",
                  "b6400f39-345a-4711-ab4f-92fd4e22cb1a",
                  "63604565-0406-460b-8cb8-1abe954b3f3a")

foods <- list(c("mammals", "birds"),
              c("cat food", "mice"),
              c("insects"),
              c("Most things, but especially anything rare or expensive"),
              c("light", "dirt"),
              c("light", "dirt"))

reaction <- function(x) {
  ifelse(x$data$info$dangerous,
         paste0("Watch out! That ", x$data$info$name, " might attack!"),
         paste0("No worries; its just a ", x$data$info$name, "."))
}

my_taxmap <- taxmap(tiger, cat, mole, human, tomato, potato,
                    data = list(info = info,
                                phylopic_ids = phylopic_ids,
                                foods = foods),
                    funcs = list(reaction = reaction))
```

In most functions that work with taxmap objects, the names of list/vector data sets, table columns, or functions can be used as if they were separate variables on their own (i.e. NSE).
In the case of functions, instead of returning the function itself, the results of the functions are returned.
To see what variables can be used this way, use `all_names`.

```{r}
all_names(my_taxmap)
```

For example using `my_taxmap$data$info$n_legs` or `n_legs` will have the same effect inside manipulation functions like `filter_taxa` described below.
This is similar to how `taxon_ranks` was used in `supertaxa` in a previous section.
To get the values of these variables, use `get_data`.

```{r}
get_data(my_taxmap)
```

Note how "taxon_names" and "dangerous" are used below.

#### Filtering

In addition to all of the functions like `subtaxa` that work with `taxonomy`, `taxmap` has a set of functions to manipulate data in a taxonomic context using functions based on **dplyr**.
Like many operations on `taxmap` objects, there are a pair of functions that modify the taxa as well as the associated data, which we call "observations".
The `filter_taxa` and `filter_obs` functions are an example of such a pair that can filter taxa and observations respectively.
For example, we can use `filter_taxa` to subset all taxa with a name starting with "t":

```{r}
filter_taxa(my_taxmap, startsWith(taxon_names, "t"))
```

There can be any number of filters that resolve to TRUE/FALSE vectors, taxon ids, or edge list indexes.
For example, below is a combination of a TRUE/FALSE vectors and taxon id filter:

```{r eval=FALSE}
filter_taxa(my_taxmap, startsWith(taxon_names, "t"), c("b", "r", "o"))
```

There are many options for `filter_taxa` that make it very flexible.
For example, the `supertaxa` option can make all the supertaxa of selected taxa be preserved.

```{r}
filter_taxa(my_taxmap, startsWith(taxon_names, "t"), supertaxa = TRUE)
```

The `filter_obs` function works in a similar way, but subsets observations in `my_taxmap$data`.

```{r}
filter_obs(my_taxmap, "info", dangerous == TRUE)
```

You can choose to filter out taxa whose observations did not pass the filter as well:

```{r}
filter_obs(my_taxmap, "info", dangerous == TRUE, drop_taxa = TRUE)
```

Note how both the taxonomy and the associated data sets were filtered. 
The `drop_obs` option can be used to specify which non-target (i.e. not `"info"`) data sets are filtered when taxa are removed.

#### Sampling

The functions `sample_n_obs` and `sample_n_taxa` are similar to `filter_obs` and `filter_taxa`, except taxa/observations are chosen  randomly.
All of the options of the "filter_" functions are available to the "sample_" functions

```{r}
set.seed(1)
sample_n_taxa(my_taxmap, 3) # "3" here is a taxon index in the edge list
set.seed(1)
sample_n_taxa(my_taxmap, 3, supertaxa = TRUE)
```


#### Adding columns

Adding columns to tabular data sets is done using `mutate_obs`.

```{r}
mutate_obs(my_taxmap, "info",
           new_col = "Im new",
           newer_col = paste0(new_col, "er!"))
```

Note how you can use newly created columns in the same call.

#### Subsetting columns

Subsetting columns in tabular data sets is done using `select_obs`.

```{r}
# Selecting a column by name
select_obs(my_taxmap, "info", dangerous)

# Selecting a column by index
select_obs(my_taxmap, "info", 3)

# Selecting a column by regular expressions (i.e. TRUE/FALSE)
select_obs(my_taxmap, "info", matches("^dange"))
```

#### Sorting

Sorting the edge list and observations is done using `arrage_taxa` and `arrange_obs`.

```{r}
arrange_taxa(my_taxmap, taxon_names)
arrange_obs(my_taxmap, "info", name)
```


#### Parsing data

The `taxmap` class has the ability to contain and manipulate very complex data.
However, this can make it difficult to parse the data into a `taxmap` object.
For this reason, there are three functions to help creating `taxmap` objects from nearly any kind of data that a taxonomy can be associated with or derived from.
The figure below shows simplified versions of how to create `taxmap` objects from different types of data in different formats.

```{r results='asis', echo = FALSE}
fig_path <- "parsing_guide.png"
width <- 800
if (knitr:::child_mode()) { # if run as a child
  fig_path <- file.path("vignettes", fig_path)
}
cat(paste0('<img src="', fig_path, '" title="parsing diagram" width="',
       width, '">'))
```


The `parse_tax_data` and `lookup_tax_data` have, in addition to the functionality above, the ability to include additional data sets that are somehow associated with the source data sets (e.g. share a common identifier). Elements in these data sets will be assigned the taxa defined in the source data, so functions like `filter_taxa` and `filter_obs` will work on all of the data set at once.
See the parsing vignette for more information.

## Parsing Hierarchy and hierarchies objects

A set of functions are available for parsing objects of class `Hierarchy` and
`hierarchies`. These functions are being ported from the CRAN package `binomen`.

The functions below are "taxonomically aware" so that you can use for example
`>` and `<` operators to filter your taxonomic names data.

### pick

`pick()` - Pick out specific taxa, while others are dropped

```{r}
ex_hierarchy1
# specific ranks by rank name
pick(ex_hierarchy1, ranks("family"))
# two elements by taxonomic name
pick(ex_hierarchy1, nms("Poaceae", "Poa"))
# two elements by taxonomic identifier
pick(ex_hierarchy1, ids(4479, 4544))
# combine types
pick(ex_hierarchy1, ranks("family"), ids(4544))
```

### pop

`pop()` - Pop out taxa, that is, drop them

```{r}
ex_hierarchy1
# specific ranks by rank name
pop(ex_hierarchy1, ranks("family"))
# two elements by taxonomic name
pop(ex_hierarchy1, nms("Poaceae", "Poa"))
# two elements by taxonomic identifier
pop(ex_hierarchy1, ids(4479, 4544))
# combine types
pop(ex_hierarchy1, ranks("family"), ids(4544))
```

### span

`span()` - Select a range of taxa, either by two names, or relational operators

```{r}
ex_hierarchy1
# keep all taxa between family and genus
# - by rank name, taxonomic name or ID
span(ex_hierarchy1, nms("Poaceae", "Poa"))

# keep all taxa greater than genus
span(ex_hierarchy1, ranks("> genus"))

# keep all taxa greater than or equal to genus
span(ex_hierarchy1, ranks(">= genus"))

# keep all taxa less than Felidae
span(ex_hierarchy2, nms("< Felidae"))

## Multiple operator statements - useful with larger classifications
ex_hierarchy3
span(ex_hierarchy3, ranks("> genus"), ranks("< phylum"))
```


## For more information

This vignette is meant to be just an outline of what `taxa` can do.
In the future, we plan to release additional, in-depth vignettes for specific topics.
More information for specific functions and examples can be found on their man pages by typing the name of the function prefixed by a `?` in an R session. For example, `?filter_taxa` will  pull up the help page for `filter_taxa`.


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R, R/documentation.R, R/taxon.R,
%   R/taxon_authority.R, R/taxon_db.R, R/taxon_db_def.R, R/taxon_id.R,
%   R/taxon_rank.R, R/taxon_rank_level.R, R/taxonomy.R
\name{vec_cast.taxa_classification}
\alias{vec_cast.taxa_classification}
\alias{vec_cast.taxa_classification.default}
\alias{vec_cast.taxa_classification.taxa_classification}
\alias{vec_cast.character.taxa_classification}
\alias{vec_cast.factor.taxa_classification}
\alias{vec_cast.integer.taxa_classification}
\alias{taxa_casting_funcs}
\alias{vec_cast.taxa_taxon}
\alias{vec_cast.taxa_taxon.default}
\alias{vec_cast.taxa_taxon.taxa_taxon}
\alias{vec_cast.taxa_taxon.character}
\alias{vec_cast.character.taxa_taxon}
\alias{vec_cast.taxa_taxon.factor}
\alias{vec_cast.factor.taxa_taxon}
\alias{vec_cast.taxa_taxon_authority}
\alias{vec_cast.taxa_taxon_authority.default}
\alias{vec_cast.taxa_taxon_authority.taxa_taxon_authority}
\alias{vec_cast.taxa_taxon_authority.character}
\alias{vec_cast.character.taxa_taxon_authority}
\alias{vec_cast.taxa_taxon_authority.factor}
\alias{vec_cast.factor.taxa_taxon_authority}
\alias{vec_cast.taxa_taxon_authority.integer}
\alias{vec_cast.data.frame.taxa_taxon_authority}
\alias{vec_cast.taxa_taxon_db}
\alias{vec_cast.taxa_taxon_db.default}
\alias{vec_cast.taxa_taxon_db.taxa_taxon_db}
\alias{vec_cast.taxa_taxon_db.character}
\alias{vec_cast.character.taxa_taxon_db}
\alias{vec_cast.taxa_taxon_db.factor}
\alias{vec_cast.factor.taxa_taxon_db}
\alias{vec_cast.taxa_taxon_db_def}
\alias{vec_cast.taxa_taxon_db_def.default}
\alias{vec_cast.taxa_taxon_db_def.taxa_taxon_db_def}
\alias{vec_cast.taxa_taxon_id}
\alias{vec_cast.taxa_taxon_id.default}
\alias{vec_cast.taxa_taxon_id.taxa_taxon_id}
\alias{vec_cast.taxa_taxon_id.character}
\alias{vec_cast.character.taxa_taxon_id}
\alias{vec_cast.taxa_taxon_id.factor}
\alias{vec_cast.factor.taxa_taxon_id}
\alias{vec_cast.taxa_taxon_id.integer}
\alias{vec_cast.integer.taxa_taxon_id}
\alias{vec_cast.taxa_taxon_id.double}
\alias{vec_cast.double.taxa_taxon_id}
\alias{vec_cast.data.frame.taxa_taxon_id}
\alias{vec_cast.taxa_taxon_rank}
\alias{vec_cast.taxa_taxon_rank.default}
\alias{vec_cast.taxa_taxon_rank.taxa_taxon_rank}
\alias{vec_cast.taxa_taxon_rank.character}
\alias{vec_cast.character.taxa_taxon_rank}
\alias{vec_cast.taxa_taxon_rank.factor}
\alias{vec_cast.factor.taxa_taxon_rank}
\alias{vec_cast.taxa_taxon_rank.double}
\alias{vec_cast.double.taxa_taxon_rank}
\alias{vec_cast.data.frame.taxa_taxon_rank}
\alias{vec_cast.taxa_taxon_rank_level}
\alias{vec_cast.taxa_taxon_rank_level.default}
\alias{vec_cast.taxa_taxon_rank_level.taxa_taxon_rank_level}
\alias{vec_cast.taxa_taxon_rank_level.character}
\alias{vec_cast.character.taxa_taxon_rank_level}
\alias{vec_cast.taxa_taxonomy}
\alias{vec_cast.taxa_taxonomy.default}
\alias{vec_cast.taxa_taxonomy.taxa_taxonomy}
\alias{vec_cast.taxa_taxonomy.character}
\alias{vec_cast.character.taxa_taxonomy}
\alias{vec_cast.taxa_taxonomy.factor}
\alias{vec_cast.factor.taxa_taxonomy}
\alias{vec_cast.taxa_taxonomy.taxa_taxon}
\alias{vec_cast.taxa_taxon.taxa_taxonomy}
\title{taxa casting functions}
\usage{
\method{vec_cast}{taxa_classification}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_classification}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_classification}{taxa_classification}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_classification}(x, to, ..., x_arg, to_arg)

\method{vec_cast.factor}{taxa_classification}(x, to, ..., x_arg, to_arg)

\method{vec_cast.integer}{taxa_classification}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxon}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon}{taxa_taxon}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon}{character}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_taxon}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon}{factor}(x, to, ..., x_arg, to_arg)

\method{vec_cast.factor}{taxa_taxon}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxon_authority}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_authority}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_authority}{taxa_taxon_authority}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_authority}{character}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_taxon_authority}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_authority}{factor}(x, to, ..., x_arg, to_arg)

\method{vec_cast.factor}{taxa_taxon_authority}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_authority}{integer}(x, to, ..., x_arg, to_arg)

\method{vec_cast.data.frame}{taxa_taxon_authority}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxon_db}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_db}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_db}{taxa_taxon_db}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_db}{character}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_taxon_db}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_db}{factor}(x, to, ..., x_arg, to_arg)

\method{vec_cast.factor}{taxa_taxon_db}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxon_db_def}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_db_def}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_db_def}{taxa_taxon_db_def}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxon_id}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_id}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_id}{taxa_taxon_id}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_id}{character}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_taxon_id}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_id}{factor}(x, to, ..., x_arg, to_arg)

\method{vec_cast.factor}{taxa_taxon_id}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_id}{integer}(x, to, ..., x_arg, to_arg)

\method{vec_cast.integer}{taxa_taxon_id}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_id}{double}(x, to, ..., x_arg, to_arg)

\method{vec_cast.double}{taxa_taxon_id}(x, to, ..., x_arg, to_arg)

\method{vec_cast.data.frame}{taxa_taxon_id}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxon_rank}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank}{taxa_taxon_rank}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank}{character}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_taxon_rank}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank}{factor}(x, to, ..., x_arg, to_arg)

\method{vec_cast.factor}{taxa_taxon_rank}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank}{double}(x, to, ..., x_arg, to_arg)

\method{vec_cast.double}{taxa_taxon_rank}(x, to, ..., x_arg, to_arg)

\method{vec_cast.data.frame}{taxa_taxon_rank}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxon_rank_level}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank_level}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank_level}{taxa_taxon_rank_level}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon_rank_level}{character}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_taxon_rank_level}(x, to, ..., x_arg, to_arg)

\method{vec_cast}{taxa_taxonomy}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxonomy}{default}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxonomy}{taxa_taxonomy}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxonomy}{character}(x, to, ..., x_arg, to_arg)

\method{vec_cast.character}{taxa_taxonomy}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxonomy}{factor}(x, to, ..., x_arg, to_arg)

\method{vec_cast.factor}{taxa_taxonomy}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxonomy}{taxa_taxon}(x, to, ..., x_arg, to_arg)

\method{vec_cast.taxa_taxon}{taxa_taxonomy}(x, to, ..., x_arg, to_arg)
}
\description{
Functions used internally for casting taxon objects to other types. They have to be exported to
work, but they are not intended to be directly used by most users.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_authority.R
\name{taxon_authority}
\alias{taxon_authority}
\title{Taxon authority class}
\usage{
taxon_authority(
  author = character(),
  date = NA,
  citation = NA,
  .names = NA,
  extract_date = TRUE
)
}
\arguments{
\item{author}{Zero or more author names.}

\item{date}{Zero or more dates.}

\item{citation}{Zero or more literature citations.}

\item{.names}{The names of the vector.}

\item{extract_date}{If \code{TRUE} (the default), then if a date is detected in the \code{author} input and
no \code{date} input is given, then the date is separated from the author input.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_authority}
}
\description{
\Sexpr[results=rd, stage=render]{taxa:::lifecycle("maturing")} Used to store information on taxon authorities, such as author names, date, and citation.
}
\examples{

# Making new objects
x <- taxon_authority(c('A', 'B', 'C'))
x <- taxon_authority(c('Cham. & Schldl.', 'L.'),
                     date = c('1827', '1753'))

# Manipulating objects
as.character(x)
x[2]
x[2] <- 'ABC'
names(x) <- c('a', 'b')
x['b'] <- 'David Bowie'
tax_author(x)[1] <- tolower(tax_author(x)[1])
tax_author(x)
tax_date(x) <- c('2000', '1234')
tax_date(x)
tax_cite(x)[2] <- c('Linnaeus, C. (1771). Mantissa plantarum altera generum.')
tax_cite(x)

# Using as columns in tables
tibble::tibble(x = x, y = 1:2)
data.frame(x = x, y = 1:2)

# Converting to tables
tibble::as_tibble(x)
as_data_frame(x)

}
\seealso{
Other classes: 
\code{\link{[.taxa_classification}()},
\code{\link{classification}()},
\code{\link{taxon_db}()},
\code{\link{taxon_id}()},
\code{\link{taxon_rank}()},
\code{\link{taxon}()}
}
\concept{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{is_leaf}
\alias{is_leaf}
\title{Check if taxa are leaves}
\usage{
is_leaf(x)
}
\arguments{
\item{x}{The object to get leaves for, such as a \link{taxonomy} object}
}
\description{
Check if each taxon is a leaf. A leaf is a taxon with no subtaxa.
subtaxa.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
is_leaf(x)

}
\seealso{
Other leaf functions: 
\code{\link{leaves}()},
\code{\link{n_leaves}()}
}
\concept{leaf functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{is_internode}
\alias{is_internode}
\title{Check if taxa are internodes}
\usage{
is_internode(x)
}
\arguments{
\item{x}{The object to get internodes for, such as a \link{taxonomy} object.}
}
\description{
Check if each taxon is an internode. An internode is a taxon with exactly one
supertaxon and one subtaxon. These taxa can be removed without losing
information on the relationships of the remaining taxa.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
is_internode(x)

}
\seealso{
Other internode functions: 
\code{\link{internodes}()}
}
\concept{internode functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_authority.R
\name{is_taxon_authority}
\alias{is_taxon_authority}
\title{Check if is a \link{taxon_authority}}
\usage{
is_taxon_authority(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is of the \link{taxon_authority} class
}
\examples{
x <- taxon_authority(c('Cham. & Schldl.', 'L.'),
                     date = c('1827', '1753'))
is_taxon_authority(x)
is_taxon_authority(1:3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R, R/documentation.R,
%   R/taxon.R, R/taxon_authority.R, R/taxon_db.R, R/taxon_id.R, R/taxon_rank.R,
%   R/taxon_rank_level.R, R/taxonomy.R
\name{vec_ptype2.taxa_classification}
\alias{vec_ptype2.taxa_classification}
\alias{vec_ptype2.taxa_classification.default}
\alias{vec_ptype2.taxa_classification.vctrs_unspecified}
\alias{vec_ptype2.taxa_classification.taxa_classification}
\alias{vec_ptype2.taxa_classification.character}
\alias{vec_ptype2.character.taxa_classification}
\alias{vec_ptype2.taxa_classification.factor}
\alias{vec_ptype2.factor.taxa_classification}
\alias{taxa_coercion_funcs}
\alias{vec_ptype2.taxa_taxon}
\alias{vec_ptype2.taxa_taxon.default}
\alias{vec_ptype2.taxa_taxon.vctrs_unspecified}
\alias{vec_ptype2.taxa_taxon.taxa_taxon}
\alias{vec_ptype2.taxa_taxon.character}
\alias{vec_ptype2.character.taxa_taxon}
\alias{vec_ptype2.taxa_taxon.factor}
\alias{vec_ptype2.factor.taxa_taxon}
\alias{vec_ptype2.taxa_taxon_authority}
\alias{vec_ptype2.taxa_taxon_authority.default}
\alias{vec_ptype2.taxa_taxon_authority.vctrs_unspecified}
\alias{vec_ptype2.taxa_taxon_authority.taxa_taxon_authority}
\alias{vec_ptype2.taxa_taxon_authority.character}
\alias{vec_ptype2.character.taxa_taxon_authority}
\alias{vec_ptype2.taxa_taxon_authority.factor}
\alias{vec_ptype2.factor.taxa_taxon_authority}
\alias{vec_ptype2.taxa_taxon_db}
\alias{vec_ptype2.taxa_taxon_db.default}
\alias{vec_ptype2.taxa_taxon_db.vctrs_unspecified}
\alias{vec_ptype2.taxa_taxon_db.taxa_taxon_db}
\alias{vec_ptype2.taxa_taxon_db.character}
\alias{vec_ptype2.character.taxa_taxon_db}
\alias{vec_ptype2.taxa_taxon_db.factor}
\alias{vec_ptype2.factor.taxa_taxon_db}
\alias{vec_ptype2.taxa_taxon_id}
\alias{vec_ptype2.taxa_taxon_id.default}
\alias{vec_ptype2.taxa_taxon_id.vctrs_unspecified}
\alias{vec_ptype2.taxa_taxon_id.taxa_taxon_id}
\alias{vec_ptype2.taxa_taxon_id.character}
\alias{vec_ptype2.character.taxa_taxon_id}
\alias{vec_ptype2.taxa_taxon_id.factor}
\alias{vec_ptype2.factor.taxa_taxon_id}
\alias{vec_ptype2.taxa_taxon_rank}
\alias{vec_ptype2.taxa_taxon_rank.default}
\alias{vec_ptype2.taxa_taxon_rank.vctrs_unspecified}
\alias{vec_ptype2.taxa_taxon_rank.taxa_taxon_rank}
\alias{vec_ptype2.taxa_taxon_rank.character}
\alias{vec_ptype2.character.taxa_taxon_rank}
\alias{vec_ptype2.taxa_taxon_rank.factor}
\alias{vec_ptype2.factor.taxa_taxon_rank}
\alias{vec_ptype2.taxa_taxon_rank_level}
\alias{vec_ptype2.taxa_taxon_rank_level.default}
\alias{vec_ptype2.taxa_taxon_rank_level.vctrs_unspecified}
\alias{vec_ptype2.taxa_taxon_rank_level.taxa_taxon_rank_level}
\alias{vec_ptype2.taxa_taxon_rank_level.character}
\alias{vec_ptype2.character.taxa_taxon_rank_level}
\alias{vec_ptype2.taxa_taxonomy}
\alias{vec_ptype2.taxa_taxonomy.default}
\alias{vec_ptype2.taxa_taxonomy.vctrs_unspecified}
\alias{vec_ptype2.taxa_taxonomy.taxa_taxonomy}
\alias{vec_ptype2.taxa_taxonomy.character}
\alias{vec_ptype2.character.taxa_taxonomy}
\alias{vec_ptype2.taxa_taxonomy.factor}
\alias{vec_ptype2.factor.taxa_taxonomy}
\title{taxa coercion functions}
\usage{
\method{vec_ptype2}{taxa_classification}(x, y, ...)

\method{vec_ptype2.taxa_classification}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_classification}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_classification}{taxa_classification}(x, y, ...)

\method{vec_ptype2.taxa_classification}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_classification}(x, y, ...)

\method{vec_ptype2.taxa_classification}{factor}(x, y, ...)

\method{vec_ptype2.factor}{taxa_classification}(x, y, ...)

\method{vec_ptype2}{taxa_taxon}(x, y, ...)

\method{vec_ptype2.taxa_taxon}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_taxon}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_taxon}{taxa_taxon}(x, y, ...)

\method{vec_ptype2.taxa_taxon}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_taxon}(x, y, ...)

\method{vec_ptype2.taxa_taxon}{factor}(x, y, ...)

\method{vec_ptype2.factor}{taxa_taxon}(x, y, ...)

\method{vec_ptype2}{taxa_taxon_authority}(x, y, ...)

\method{vec_ptype2.taxa_taxon_authority}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_taxon_authority}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_taxon_authority}{taxa_taxon_authority}(x, y, ...)

\method{vec_ptype2.taxa_taxon_authority}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_taxon_authority}(x, y, ...)

\method{vec_ptype2.taxa_taxon_authority}{factor}(x, y, ...)

\method{vec_ptype2.factor}{taxa_taxon_authority}(x, y, ...)

\method{vec_ptype2}{taxa_taxon_db}(x, y, ...)

\method{vec_ptype2.taxa_taxon_db}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_taxon_db}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_taxon_db}{taxa_taxon_db}(x, y, ...)

\method{vec_ptype2.taxa_taxon_db}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_taxon_db}(x, y, ...)

\method{vec_ptype2.taxa_taxon_db}{factor}(x, y, ...)

\method{vec_ptype2.factor}{taxa_taxon_db}(x, y, ...)

\method{vec_ptype2}{taxa_taxon_id}(x, y, ...)

\method{vec_ptype2.taxa_taxon_id}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_taxon_id}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_taxon_id}{taxa_taxon_id}(x, y, ...)

\method{vec_ptype2.taxa_taxon_id}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_taxon_id}(x, y, ...)

\method{vec_ptype2.taxa_taxon_id}{factor}(x, y, ...)

\method{vec_ptype2.factor}{taxa_taxon_id}(x, y, ...)

\method{vec_ptype2}{taxa_taxon_rank}(x, y, ...)

\method{vec_ptype2.taxa_taxon_rank}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_taxon_rank}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_taxon_rank}{taxa_taxon_rank}(x, y, ...)

\method{vec_ptype2.taxa_taxon_rank}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_taxon_rank}(x, y, ...)

\method{vec_ptype2.taxa_taxon_rank}{factor}(x, y, ...)

\method{vec_ptype2.factor}{taxa_taxon_rank}(x, y, ...)

\method{vec_ptype2}{taxa_taxon_rank_level}(x, y, ...)

\method{vec_ptype2.taxa_taxon_rank_level}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_taxon_rank_level}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_taxon_rank_level}{taxa_taxon_rank_level}(x, y, ...)

\method{vec_ptype2.taxa_taxon_rank_level}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_taxon_rank_level}(x, y, ...)

\method{vec_ptype2}{taxa_taxonomy}(x, y, ...)

\method{vec_ptype2.taxa_taxonomy}{default}(x, y, ..., x_arg = "", y_arg = "")

\method{vec_ptype2.taxa_taxonomy}{vctrs_unspecified}(x, y, ...)

\method{vec_ptype2.taxa_taxonomy}{taxa_taxonomy}(x, y, ...)

\method{vec_ptype2.taxa_taxonomy}{character}(x, y, ...)

\method{vec_ptype2.character}{taxa_taxonomy}(x, y, ...)

\method{vec_ptype2.taxa_taxonomy}{factor}(x, y, ...)

\method{vec_ptype2.factor}{taxa_taxonomy}(x, y, ...)
}
\description{
Functions used internally for coercing taxon objects and other objects to common data types.
They have to be exported to work, but they are not intended to be directly used by most users.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{is_taxonomy}
\alias{is_taxonomy}
\title{Check if something is a \link{taxonomy}}
\usage{
is_taxonomy(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is of the \link{taxonomy} class
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
is_taxonomy(x)
is_taxonomy(1:2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{n_subtaxa}
\alias{n_subtaxa}
\title{Number of subtaxa per taxon}
\usage{
n_subtaxa(x, subset = NULL, max_depth = NULL, include = FALSE)
}
\arguments{
\item{x}{The object to get subtaxa for, such as a \link{taxonomy} object.}

\item{subset}{The subset of the tree to search. Can be indexes or names.}

\item{max_depth}{The number of ranks to traverse. For example, \code{max_depth = 1} returns only immediate subtaxa. By default (NULL) information for all
subtaxa is returned (i.e. subtaxa of subtaxa, etc).}

\item{include}{If \code{TRUE}, include information for each taxon in the output.}
}
\description{
Get the number of subtaxa per taxon.
}
\examples{
# Generate example data
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))

# Find number of subtaxa (including subtaxa of subtaxa, etc)
n_subtaxa(x)

# Find the number of subtaxa one rank below each taxon
n_subtaxa(x, max_depth = 1)

# Only return data for some taxa (faster than subsetting the whole result)
n_subtaxa(x, subset = 1:3)

}
\seealso{
Other subtaxa functions: 
\code{\link{subtaxa}()}
}
\concept{subtaxa functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{stems}
\alias{stems}
\title{Get stems}
\usage{
stems(x, value = NULL, ...)
}
\arguments{
\item{x}{An object with taxonomic relationships, like \link{taxonomy} objects.}

\item{value}{Something to return instead of indexes. Must be the same length as the number of taxa.}

\item{...}{Additional arguments.}
}
\description{
Get stem indexes for each taxon or another per-taxon value.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris'),
              supertaxa = c(NA, 1, 2, 3, 3))
x <- c(x, x)
stems(x)
stems(x, value = tax_name(x))

}
\seealso{
Other taxonomy functions: 
\code{\link{internodes}()},
\code{\link{leaves}()},
\code{\link{roots}()},
\code{\link{subtaxa}()},
\code{\link{supertaxa}()}

Other stem functions: 
\code{\link{is_stem}()}
}
\concept{stem functions}
\concept{taxonomy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{new_taxonomy}
\alias{new_taxonomy}
\title{Minimal taxonomy constructor}
\usage{
new_taxonomy(taxa = taxon(), supertaxa = integer())
}
\arguments{
\item{taxa}{A \link{taxon} vector.}

\item{supertaxa}{The indexes of \code{taxa} for each taxon's supertaxon.}
}
\value{
An \code{S3} object of class \code{taxa_taxon}
}
\description{
Minimal taxonomy constructor for internal use. Only use when the input is known to be valid since
few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printing_style.R
\name{font_tax_name}
\alias{font_tax_name}
\title{Taxon name formatting in print methods}
\usage{
font_tax_name(text)
}
\arguments{
\item{text}{What to print}
}
\description{
A simple wrapper to make changing the formatting of text printed easier.
This is used for taxon names.
}
\seealso{
Other printer fonts: 
\code{\link{font_default}()},
\code{\link{font_na}()},
\code{\link{font_punct}()},
\code{\link{font_secondary}()}
}
\concept{printer fonts}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_id.R
\name{taxa_taxon_id-class}
\alias{taxa_taxon_id-class}
\alias{taxa_taxon_id}
\title{Taxon ID class}
\description{
Taxon ID class. See \link{taxon_id} for more information
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{is_classification}
\alias{is_classification}
\title{Check if is a classification}
\usage{
is_classification(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is the classification class
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_id.R
\name{is_taxon_id}
\alias{is_taxon_id}
\title{Check if something is a \link{taxon_id} object}
\usage{
is_taxon_id(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is of the \link{taxon_id} class
}
\examples{
x <- taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi')
is_taxon_id(x)
is_taxon_id(1:3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{n_supertaxa}
\alias{n_supertaxa}
\title{Number of supertaxa per taxon}
\usage{
n_supertaxa(x, subset = NULL, max_depth = NULL, include = FALSE)
}
\arguments{
\item{x}{The object to get supertaxa for, such as a \link{taxonomy} object.}

\item{subset}{The subset of the tree to search for roots to that subset. Can
be indexes or names.}

\item{max_depth}{The number of levels to traverse. For example, \code{max_depth = 1} returns only immediate supertaxa. By default (NULL) information for all
supertaxa is returned.}

\item{include}{If \code{TRUE}, include information for each taxon in the output.}
}
\description{
Get the number of supertaxa each taxon is contained in.
}
\examples{
# Generate example data
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))

# Find number of supertaxa each taxon is contained in
n_supertaxa(x)

# Only return data for some taxa (faster than subsetting the whole result)
n_supertaxa(x, subset = 1:3)

}
\seealso{
Other supertaxa functions: 
\code{\link{supertaxa}()}
}
\concept{supertaxa functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank_level.R
\name{impute_order_na}
\alias{impute_order_na}
\title{Fill in NA values in sequence}
\usage{
impute_order_na(order, inc = 1)
}
\arguments{
\item{order}{An ascending sequences, possibly with NAs}

\item{inc}{The increment size to use for values in NA blocks at the start and end of the sequence.}
}
\description{
Fill in the NA values in a ascending sequence based on nearby non-NA values.
Used to guess the order values for unknown ranks based on the values of known ranks.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{roots}
\alias{roots}
\title{Get root taxa}
\usage{
roots(x, subset = NULL)
}
\arguments{
\item{x}{An object containing taxonomic relationships, such as \link{taxonomy} objects.}

\item{subset}{The subset of the tree to search for roots to that subset. Can be indexes or names.}
}
\description{
Get the indexes of root taxa in a taxonomy.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
roots(x)
roots(x, subset = 2:8)

}
\seealso{
Other taxonomy functions: 
\code{\link{internodes}()},
\code{\link{leaves}()},
\code{\link{stems}()},
\code{\link{subtaxa}()},
\code{\link{supertaxa}()}

Other root functions: 
\code{\link{is_root}()}
}
\concept{root functions}
\concept{taxonomy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_generics.R
\name{as_data_frame}
\alias{as_data_frame}
\title{Convert a taxa object to a \code{data.frame}}
\usage{
as_data_frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  stringsAsFactors = default.stringsAsFactors()
)
}
\arguments{
\item{x}{An object defined by taxa, such as \link{taxon} or \link{taxon_id}}

\item{row.names}{\code{NULL} or a character vector giving the row
    names for the data frame.  Missing values are not allowed.}

\item{optional}{logical. If \code{TRUE}, setting row names and
    converting column names (to syntactic names: see
    \code{\link[base]{make.names}}) is optional.  Note that all of \R's
    \pkg{base} package \code{as.data.frame()} methods use
    \code{optional} only for column names treatment, basically with the
    meaning of \code{\link[base]{data.frame}(*, check.names = !optional)}.
    See also the \code{make.names} argument of the \code{matrix} method.}

\item{...}{additional arguments to be passed to or from methods.}

\item{stringsAsFactors}{logical: should the character vector be converted
    to a factor?}
}
\description{
Convert the information in a taxa object to a \code{data.frame} using base R
vectors as columns. Use \link{as_tibble} to convert to tibbles.
}
\examples{
x <- taxon(name = c('Homo sapiens', 'Bacillus', 'Ascomycota', 'Ericaceae'),
           rank = c('species', 'genus', 'phylum', 'family'),
           id = taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi'),
           auth = c('Linnaeus, 1758', 'Cohn 1872', NA, 'Juss., 1789'))
as_data_frame(x)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_id.R
\name{printed_taxon_id}
\alias{printed_taxon_id}
\title{Prepare taxon_id for printing}
\usage{
printed_taxon_id(x, color = FALSE)
}
\arguments{
\item{color}{Use color?}
}
\value{
character
}
\description{
Prepare taxon_id for printing. Makes color optional.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{is_stem}
\alias{is_stem}
\title{Check if taxa are stems}
\usage{
is_stem(x)
}
\arguments{
\item{x}{An object with taxonomic relationships, like \link{taxonomy} objects.}
}
\description{
Check if each taxon is a stem. A stem is any taxa from a root to the first taxon with multiple subtaxa.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris'),
              supertaxa = c(NA, 1, 2, 3, 3))
is_stem(x)

}
\seealso{
Other stem functions: 
\code{\link{stems}()}
}
\concept{stem functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printing_style.R
\name{font_secondary}
\alias{font_secondary}
\title{Font for secondary data}
\usage{
font_secondary(text)
}
\arguments{
\item{text}{What to print.}
}
\description{
A wrapper to make changing the formatting of text printed easier.
This is used for print data associated with other data.
}
\seealso{
Other printer fonts: 
\code{\link{font_default}()},
\code{\link{font_na}()},
\code{\link{font_punct}()},
\code{\link{font_tax_name}()}
}
\concept{printer fonts}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.R
\name{taxon}
\alias{taxon}
\title{Taxon class}
\usage{
taxon(name = character(0), rank = NA, id = NA, auth = NA, .names = NA, ...)
}
\arguments{
\item{name}{The names of taxa. Inputs with be coerced into a \link{character} vector if anything else
is given.}

\item{rank}{The ranks of taxa. Inputs with be coerced into a \link{taxon_rank} vector if anything else
is given.}

\item{id}{The ids of taxa. These should be unique identifier and are usually associated with a
database. Inputs with be coerced into a \link{taxon_id} vector if anything else is given.}

\item{auth}{The authority of the taxon. Inputs with be coerced into a \link{taxon_authority} vector if
anything else is given.}

\item{.names}{The names of the vector.}

\item{...}{Additional arguments.}
}
\value{
An \code{S3} object of class \code{taxa_taxon}
}
\description{
\Sexpr[results=rd, stage=render]{taxa:::lifecycle("maturing")}
Used to store information about taxa, such as names, ranks, and IDs.
}
\examples{

# Create taxon name vector
x <- taxon(c('A', 'B', 'C'))
x <- taxon(name = c('Homo sapiens', 'Bacillus', 'Ascomycota', 'Ericaceae'),
           rank = c('species', 'genus', 'phylum', 'family'),
           id = taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi'),
           auth = c('Linnaeus, 1758', 'Cohn 1872', NA, 'Juss., 1789'))
names(x) <- c('a', 'b', 'c', 'd')

# Get parts of the taxon name vector
tax_name(x)
tax_rank(x)
tax_id(x)
tax_db(x)
tax_auth(x)
tax_author(x)
tax_date(x)
tax_cite(x)

# Set parts of the taxon name vector
tax_name(x) <- tolower(tax_name(x))
tax_rank(x)[1] <- NA
tax_name(x)['b'] <- 'Billy'
tax_id(x) <- '9999'
tax_db(x) <- 'itis'
tax_auth(x) <- NA
tax_author(x)[2:3] <- c('Joe', 'Billy')
tax_date(x) <- c('1999', '2013', '1796', '1899')
tax_cite(x)[1] <- 'Linnaeus, C. (1771). Mantissa plantarum altera generum.'

# Manipulate taxon name vectors
x[1:3]
x[tax_rank(x) > 'family']
x['b'] <- NA
x[c('c', 'd')] <- 'unknown'
is.na(x)

# Use as columns in tables
tibble::tibble(x = x, y = 1:4)
data.frame(x = x, y = 1:4)

# Converting to tables
tibble::as_tibble(x)
as_data_frame(x)

}
\seealso{
Other classes: 
\code{\link{[.taxa_classification}()},
\code{\link{classification}()},
\code{\link{taxon_authority}()},
\code{\link{taxon_db}()},
\code{\link{taxon_id}()},
\code{\link{taxon_rank}()}
}
\concept{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_db_def.R
\name{new_taxon_db_def}
\alias{new_taxon_db_def}
\title{Minimal taxon_db_def constructor}
\usage{
new_taxon_db_def(
  name = character(),
  url = character(),
  desc = character(),
  id_regex = character(),
  rank_levels = list()
)
}
\arguments{
\item{name}{Name of the database in lower case. Inputs will be transformed to a \code{character} vector.}

\item{url}{URL of the database website. Inputs will be transformed to a \code{character} vector.}

\item{desc}{Description of the database. Inputs will be transformed to a \code{character} vector.}

\item{id_regex}{A regular expression for taxon IDs of the database. Inputs will be transformed to a \code{character} vector.}

\item{rank_levels}{Valid taxonomic ranks for the database. Should be a list of \code{numeric} vectors named by taxonomic ranks.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_db_def}
}
\description{
Minimal taxon_db_def constructor for internal use. Only use when the input is known to be valid since
few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_id.R
\name{new_taxon_id}
\alias{new_taxon_id}
\title{Minimal taxon_id constructor}
\usage{
new_taxon_id(.names = NULL, id = character(), db = taxon_db())
}
\arguments{
\item{.names}{The names to apply to the vector}

\item{id}{Zero or more taxonomic ids. Inputs will be transformed to a \code{character} vector.}

\item{db}{The name(s) of the database(s) associated with the IDs. If not \code{NA} (the
default), the input must consist of names of databases in \link{database_ref}. The length must be
0, 1, or equal to the number of IDs.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_id}
}
\description{
Minimal taxon_id constructor for internal use. Only use when the input is known to be valid since
few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.R
\name{taxa_taxon-class}
\alias{taxa_taxon-class}
\alias{taxa_taxon}
\title{Taxon class}
\description{
Taxon class. See \link{taxon} for more information
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{n_leaves}
\alias{n_leaves}
\title{Number of leaves per taxon}
\usage{
n_leaves(x)
}
\arguments{
\item{x}{The object to get leaves for, such as a \link{taxonomy} object}
}
\description{
Get the number of leaves per taxon. A leaf is a taxon with no subtaxa.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
n_leaves(x)

}
\seealso{
Other leaf functions: 
\code{\link{is_leaf}()},
\code{\link{leaves}()}
}
\concept{leaf functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{delete_unused_class_taxa}
\alias{delete_unused_class_taxa}
\title{Removes taxa from the taxonomy of a classification that are not used by any of the instances}
\usage{
delete_unused_class_taxa(x)
}
\description{
Removes taxa from the taxonomy of a classification that are not used by any of the instances
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{printed_taxonomy}
\alias{printed_taxonomy}
\title{Prepare taxonomy for printing}
\usage{
printed_taxonomy(x, color = FALSE)
}
\arguments{
\item{color}{Use color?}
}
\value{
character
}
\description{
Prepare taxonomy for printing. Makes color optional.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank.R
\name{new_taxon_rank}
\alias{new_taxon_rank}
\title{Minimal taxon_rank constructor}
\usage{
new_taxon_rank(rank = character(), levels = taxon_rank_level())
}
\arguments{
\item{rank}{Zero or more taxonomic rank names. Inputs will be transformed to
a \code{character} vector.}

\item{levels}{A named numeric vector indicating the names and orders of
possible taxonomic ranks. Higher numbers indicate for fine-scale groupings.
Ranks of unknown order can be indicated with \code{NA} instead of a number.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_rank}
}
\description{
Minimal taxon_rank constructor for internal use. Only use when the input is
known to be valid since few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.R
\name{new_taxon}
\alias{new_taxon}
\title{Minimal taxon constructor}
\usage{
new_taxon(
  .names = NULL,
  name = character(),
  rank = taxon_rank(),
  id = taxon_id(),
  auth = taxon_authority(),
  ...
)
}
\arguments{
\item{.names}{The names of the vector.}

\item{name}{The names of taxa as a \link{character} vector.}

\item{rank}{The ranks of taxa as a \link{taxon_rank} vector.}

\item{id}{The ids of taxa as a \link{taxon_id} vector.}

\item{auth}{The authority of the taxon as a \link{taxon_authority} vector.}
}
\value{
An \code{S3} object of class \code{taxa_taxon}
}
\description{
Minimal taxon constructor for internal use. Only use when the input is known to be valid since
few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank.R
\name{rank_level_color_funcs}
\alias{rank_level_color_funcs}
\title{Get font color for levels}
\usage{
rank_level_color_funcs(levels)
}
\description{
Make list of crayon style functions to print taxon rank levels in color.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank.R
\name{printed_taxon_rank}
\alias{printed_taxon_rank}
\title{Prepare taxon_rank for printing}
\usage{
printed_taxon_rank(x, color = FALSE)
}
\arguments{
\item{color}{Use color?}
}
\value{
character
}
\description{
Prepare taxon_rank for printing. Makes color optional.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank_level.R
\name{is_taxon_rank_level}
\alias{is_taxon_rank_level}
\title{Check if is a taxon id}
\usage{
is_taxon_rank_level(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is the taxon id class
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{internodes}
\alias{internodes}
\title{Get internodes}
\usage{
internodes(x)
}
\arguments{
\item{x}{The object to get internodes for, such as a \link{taxonomy} object.}
}
\description{
Get internodes indexes for each taxon or another per-taxon value. An
internode is a taxon with exactly one supertaxon and one subtaxon. These taxa
can be removed without losing information on the relationships of the
remaining taxa.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
internodes(x)

}
\seealso{
Other taxonomy functions: 
\code{\link{leaves}()},
\code{\link{roots}()},
\code{\link{stems}()},
\code{\link{subtaxa}()},
\code{\link{supertaxa}()}

Other internode functions: 
\code{\link{is_internode}()}
}
\concept{internode functions}
\concept{taxonomy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation.R, R/taxon.R, R/taxon_id.R,
%   R/taxon_rank.R, R/taxonomy.R
\name{taxa_comparison_funcs}
\alias{taxa_comparison_funcs}
\alias{vec_proxy_equal.taxa_taxon}
\alias{vec_proxy_equal.taxa_taxon_id}
\alias{vec_proxy_compare.taxa_taxon_rank}
\alias{vec_proxy_equal.taxa_taxon_rank}
\alias{Ops.taxa_taxon_rank}
\alias{vec_proxy_equal.taxa_taxonomy}
\title{taxa comparison functions}
\usage{
\method{vec_proxy_equal}{taxa_taxon}(x, ...)

\method{vec_proxy_equal}{taxa_taxon_id}(x, ...)

\method{vec_proxy_compare}{taxa_taxon_rank}(x, ...)

\method{vec_proxy_equal}{taxa_taxon_rank}(x, ...)

\method{Ops}{taxa_taxon_rank}(e1, e2)

\method{vec_proxy_equal}{taxa_taxonomy}(x, ...)
}
\description{
Functions used internally for casting taxon objects to other types. They have to be exported to
work, but they are not intended to be directly used by most users.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{printed_classification}
\alias{printed_classification}
\title{Prepare classification for printing}
\usage{
printed_classification(x, color = FALSE)
}
\arguments{
\item{color}{Use color?}
}
\value{
character
}
\description{
Prepare classification for printing. Makes color optional.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxon_id.R, R/taxonomy.R
\name{tax_db.taxa_classification}
\alias{tax_db.taxa_classification}
\alias{tax_db<-.taxa_classification}
\alias{tax_db}
\alias{tax_db<-}
\alias{tax_db.taxa_taxon}
\alias{tax_db<-.taxa_taxon}
\alias{tax_db.taxa_taxon_id}
\alias{tax_db<-.taxa_taxon_id}
\alias{tax_db.taxa_taxonomy}
\alias{tax_db<-.taxa_taxonomy}
\title{Set and get taxon ID databases}
\usage{
\method{tax_db}{taxa_classification}(x)

\method{tax_db}{taxa_classification}(x) <- value

tax_db(x)

tax_db(x) <- value

\method{tax_db}{taxa_taxon}(x)

\method{tax_db}{taxa_taxon}(x) <- value

\method{tax_db}{taxa_taxon_id}(x)

\method{tax_db}{taxa_taxon_id}(x) <- value

\method{tax_db}{taxa_taxonomy}(x)

\method{tax_db}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon authority dates.}

\item{value}{The taxon citations to set. Inputs will be coerced into a \link{taxon_db} vector.}
}
\description{
Set and get the taxon ID databases in objects that have them, such as \link{taxon_id} objects.
}
\examples{
x <- taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi')
tax_db(x)
tax_db(x) <- 'nbn'
tax_db(x)[2] <- 'itis'

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{subtaxa}
\alias{subtaxa}
\title{Get subtaxa}
\usage{
subtaxa(x, subset = NULL, max_depth = NULL, include = FALSE, value = NULL, ...)
}
\arguments{
\item{x}{The object to get subtaxa for, such as a \link{taxonomy} object.}

\item{subset}{The subset of the tree to search. Can be indexes or names.}

\item{max_depth}{The number of ranks to traverse. For example, \code{max_depth = 1} returns only immediate subtaxa. By default (NULL) information for all
subtaxa is returned (i.e. subtaxa of subtaxa, etc).}

\item{include}{If \code{TRUE}, include information for each taxon in the output.}

\item{value}{Something to return instead of indexes. Must be the same length
as the number of taxa.}

\item{...}{Additional arguments.}
}
\description{
Get subtaxa indexes for each taxon or another per-taxon value. Subtaxa are
taxa contained within a taxon.
}
\examples{
# Generate example data
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))

# The indexes of all subtaxa (with subtaxa of subtaxa, etc) for each taxon
subtaxa(x)

# The indexes of immediate subtaxa (without subtaxa of subtaxa, etc) for each taxon
subtaxa(x, max_depth = 1)

# Return something other than index
subtaxa(x, value = tax_name(x))

# Include each taxon with its subtaxa
subtaxa(x, value = tax_name(x), include = TRUE)

# Only return data for some taxa (faster than subsetting the whole result)
subtaxa(x, subset = 3)

}
\seealso{
Other taxonomy functions: 
\code{\link{internodes}()},
\code{\link{leaves}()},
\code{\link{roots}()},
\code{\link{stems}()},
\code{\link{supertaxa}()}

Other subtaxa functions: 
\code{\link{n_subtaxa}()}
}
\concept{subtaxa functions}
\concept{taxonomy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{num_range}
\alias{num_range}
\title{dplyr select_helpers}
\description{
dplyr select_helpers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_authority.R
\name{taxa_taxon_authority-class}
\alias{taxa_taxon_authority-class}
\alias{taxa_taxon_authority}
\title{Taxon authority class}
\description{
Taxon authority class. See \link{taxon_authority} for more information
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{ends_with}
\alias{ends_with}
\title{dplyr select_helpers}
\description{
dplyr select_helpers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{leaves}
\alias{leaves}
\title{Get leaves}
\usage{
leaves(x, value = NULL, ...)
}
\arguments{
\item{x}{The object to get leaves for, such as a \link{taxonomy} object}

\item{value}{Something to return instead of indexes. Must be the same length as the number of taxa.}

\item{...}{Additional arguments.}
}
\description{
Get leaves indexes for each taxon or another per-taxon value. Leaves are taxa with no subtaxa.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
leaves(x)
leaves(x, value = tax_name(x))

}
\seealso{
Other taxonomy functions: 
\code{\link{internodes}()},
\code{\link{roots}()},
\code{\link{stems}()},
\code{\link{subtaxa}()},
\code{\link{supertaxa}()}

Other leaf functions: 
\code{\link{is_leaf}()},
\code{\link{n_leaves}()}
}
\concept{leaf functions}
\concept{taxonomy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank.R
\name{taxa_taxon_rank-class}
\alias{taxa_taxon_rank-class}
\alias{taxa_taxon_rank}
\title{Taxon rank class}
\description{
Taxon rank class. See \link{taxon_rank} for more information
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{data}
\name{db_ref}
\alias{db_ref}
\title{Valid taxonomy databases}
\format{
An object of class \code{list} of length 3.
}
\usage{
db_ref
}
\description{
This defines the valid taxonomic databases that can be used in \link{taxon_db}
objects and objects that use \link{taxon_db} objects, such as \link{taxon_id} and
\link{taxon}. \code{db_ref$get} can be used to see information for the databases. Users
can add their own custom databases to the list using \code{db_ref$set}. For each
database the following information is included:
\itemize{
\item The URL for the website associated with the database
\item A short description
\item The regular expression that defines valid taxon IDs
\item The ranks used in the database if specified
}
}
\section{Attribution}{


This code is based on the code handling options in \link{knitr}.
}

\examples{

# List all database definitions
db_ref$get()

# Get a specific database definition
db_ref$get('ncbi')

# Add or overwrite a database definition
db_ref$set(
  name = "my_new_database",
  url = "http://www.my_tax_database.com",
  desc = "I just made this up",
  id_regex = ".*"
)

# Reset definitions to default values
db_ref$reset()

}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{duplicated_index_taxonomy}
\alias{duplicated_index_taxonomy}
\title{Returns the index of the first occurrence of each unique taxon in a taxonomy}
\usage{
duplicated_index_taxonomy(x, proxy_func = NULL)
}
\description{
Returns the index of the first occurrence of each unique taxon in a taxonomy
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.R
\name{as_taxon}
\alias{as_taxon}
\title{Convert to a \link{taxon} vector}
\usage{
as_taxon(x, ...)
}
\arguments{
\item{x}{An object to be converted to a taxon vector}

\item{...}{Additional parameters.}
}
\description{
Convert other objects to \link{taxon} vectors. Compatible base R vectors can also
be converted using the \link[=taxon]{taxon constructor}.
}
\examples{

# Convert a taxonomy object to a taxon vector
x <- taxonomy(taxon(name = c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                             'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
                    rank = c('order', 'family', 'genus', 'species',
                             'species', 'family', 'genus', 'species'),
                    id = taxon_id(c('33554', '9681', '9688', '9689',
                                    '9694', '9632', '9639', '9644'),
                                  db = 'ncbi'),
                    auth = c('Bowdich, 1821', 'Fischer de Waldheim, 1817', 'Oken, 1816', 'L., 1758',
                             'L., 1758', 'Fischer de Waldheim, 1817', 'L., 1758', 'L., 1758')),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
names(x) <- letters[1:8]
as_taxon(x)

# Convert base R vectors
as_taxon(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo'))
as_taxon(factor(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo')))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_db_def.R
\name{is_valid_regex}
\alias{is_valid_regex}
\title{Check regex validity}
\usage{
is_valid_regex(text)
}
\arguments{
\item{text}{The putative regex to check.}
}
\description{
Check if a regular expression is valid
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{contains}
\alias{contains}
\title{dplyr select_helpers}
\description{
dplyr select_helpers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxon_authority.R, R/taxonomy.R
\name{tax_author.taxa_classification}
\alias{tax_author.taxa_classification}
\alias{tax_author<-.taxa_classification}
\alias{tax_author}
\alias{tax_author<-}
\alias{tax_author.taxa_taxon}
\alias{tax_author<-.taxa_taxon}
\alias{tax_author<-.taxa_taxon_authority}
\alias{tax_author.taxa_taxon_authority}
\alias{tax_author.taxa_taxonomy}
\alias{tax_author<-.taxa_taxonomy}
\title{Set and get taxon authors}
\usage{
\method{tax_author}{taxa_classification}(x)

\method{tax_author}{taxa_classification}(x) <- value

tax_author(x)

tax_author(x) <- value

\method{tax_author}{taxa_taxon}(x)

\method{tax_author}{taxa_taxon}(x) <- value

\method{tax_author}{taxa_taxon_authority}(x) <- value

\method{tax_author}{taxa_taxon_authority}(x)

\method{tax_author}{taxa_taxonomy}(x)

\method{tax_author}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon authors.}

\item{value}{The taxon authors to set. Inputs will be coerced into a \link{character} vector.}
}
\description{
Set and get taxon authors in objects that have them, such as \link{taxon_authority} objects.
}
\examples{
x <- taxon_authority(c('Cham. & Schldl.', 'L.'),
                     date = c('1827', '1753'))
tax_author(x)
tax_author(x)[1] <- "Billy"
tax_author(x) <- tolower(tax_author(x))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.R
\name{is_taxon}
\alias{is_taxon}
\title{Check if something is a \link{taxon} object}
\usage{
is_taxon(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is of the \link{taxon} class
}
\examples{
x <- taxon(c('A', 'B', 'C'))
is_taxon(x)
is_taxon(1:2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_db.R
\name{new_taxon_db}
\alias{new_taxon_db}
\title{Minimal taxon_db constructor}
\usage{
new_taxon_db(db = character(), ...)
}
\arguments{
\item{db}{Zero or more taxonomic database names. Should be a name contained in
\code{names(db_ref)}. Inputs will be transformed to a \code{character} vector.}

\item{...}{Additional arguments.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_db}
}
\description{
Minimal taxon_db constructor for internal use. Only use when the input is known to be valid since
few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank.R
\name{is_taxon_rank}
\alias{is_taxon_rank}
\title{Check if something is a \link{taxon_rank}}
\usage{
is_taxon_rank(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is of the \link{taxon_rank} class
}
\examples{
x <- taxon_rank(c('species', 'species', 'phylum', 'family'))
is_taxon_rank(x)
is_taxon_rank(1:3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_authority.R
\name{new_taxon_authority}
\alias{new_taxon_authority}
\title{Minimal taxon_authority constructor}
\usage{
new_taxon_authority(
  .names = NULL,
  author = character(),
  date = character(),
  citation = character()
)
}
\arguments{
\item{.names}{The names of the vector.}

\item{author}{Zero or more author names.}

\item{date}{Zero or more dates.}

\item{citation}{Zero or more literature citations.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_authority}
}
\description{
Minimal taxon_authority constructor for internal use. Only use when the input is known to be valid since
few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.R
\name{printed_taxon}
\alias{printed_taxon}
\title{Prepare taxon for printing}
\usage{
printed_taxon(x, color = FALSE)
}
\arguments{
\item{color}{Use color?}
}
\value{
character
}
\description{
Prepare taxon for printing. Makes color optional.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R, R/taxonomy.R
\name{[.taxa_classification}
\alias{[.taxa_classification}
\alias{[[.taxa_classification}
\alias{taxonomy}
\alias{names.taxa_taxonomy}
\alias{names<-.taxa_taxonomy}
\alias{[.taxa_taxonomy}
\alias{[[.taxa_taxonomy}
\title{Taxonomy class}
\usage{
\method{[}{taxa_classification}(x, ...)

\method{[[}{taxa_classification}(x, i)

taxonomy(taxa = taxon(), supertaxa = NA, .names = NULL)

\method{names}{taxa_taxonomy}(x)

\method{names}{taxa_taxonomy}(x) <- value

\method{[}{taxa_taxonomy}(x, ..., subtaxa = TRUE, supertaxa = FALSE, invert = FALSE)

\method{[[}{taxa_taxonomy}(x, i, ..., subtaxa = TRUE, supertaxa = FALSE, invert = FALSE)
}
\arguments{
\item{taxa}{A \link{taxon} vector or something that can be converted to a \link{taxon} vector.}

\item{supertaxa}{The indexes of \code{taxa} for each taxon's supertaxon.}

\item{.names}{The names of the vector (not the names of taxa).}
}
\value{
An \code{S3} object of class \code{taxa_taxon}
}
\description{
\Sexpr[results=rd, stage=render]{taxa:::lifecycle("experimental")}
Used to store information about a set of taxa forming a taxonomic tree.
}
\examples{

x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))

x <- taxonomy(taxon(name = c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                             'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
                    rank = c('order', 'family', 'genus', 'species',
                             'species', 'family', 'genus', 'species'),
                    id = taxon_id(c('33554', '9681', '9688', '9689',
                                    '9694', '9632', '9639', '9644'),
                                  db = 'ncbi'),
                    auth = c('Bowdich, 1821', 'Fischer de Waldheim, 1817', 'Oken, 1816', 'L., 1758',
                             'L., 1758', 'Fischer de Waldheim, 1817', 'L., 1758', 'L., 1758')),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
names(x) <- letters[1:8]

# Subset taxonomy vector
x[2] # By default, all subtaxa are included
x['b'] # Names can also be used
x[2:3, subtaxa = FALSE] # Disable subtaxa
x[3, supertaxa = TRUE] # include supertaxa
x[is_leaf(x)] # Subset by logical vector

# Get parts of the taxonomy vector
tax_name(x)
tax_rank(x)
tax_id(x)
tax_db(x)
tax_auth(x)
tax_author(x)
tax_date(x)
tax_cite(x)

# Set parts of the taxonomy vector
tax_name(x) <- tolower(tax_name(x))
tax_rank(x)[1] <- NA
tax_id(x) <- '9999'
tax_db(x) <- 'itis'
tax_auth(x) <- NA
tax_author(x)[2:3] <- c('Joe', 'Billy')
tax_date(x) <- c('1999', '2013', '1796', '1899',
                 '1997', '2003', '1996', '1859')
tax_cite(x)['c'] <- 'Linnaeus, C. (1771). Mantissa plantarum altera generum.'

# Convert to table
tibble::as_tibble(x)
as_data_frame(x)

# Get taxonomy attributes
subtaxa(x)
subtaxa(x, value = tax_name(x))
subtaxa(x, value = as_taxon(x))
n_subtaxa(x)
supertaxa(x)
n_supertaxa(x)
leaves(x)
n_leaves(x)
is_leaf(x)
stems(x)
is_stem(x)
roots(x)
is_root(x)
internodes(x)
is_internode(x)

}
\seealso{
Other classes: 
\code{\link{classification}()},
\code{\link{taxon_authority}()},
\code{\link{taxon_db}()},
\code{\link{taxon_id}()},
\code{\link{taxon_rank}()},
\code{\link{taxon}()}
}
\concept{classes}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{duplicated_index}
\alias{duplicated_index}
\title{Returns the index of the first occurrence of each unique element}
\usage{
duplicated_index(x, proxy_func = NULL)
}
\description{
Returns the index of the first occurrence of each unique element
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_authority.R
\name{printed_taxon_authority}
\alias{printed_taxon_authority}
\title{Prepare taxon_authority for printing}
\usage{
printed_taxon_authority(x, color = FALSE)
}
\arguments{
\item{color}{Use color?}
}
\value{
character
}
\description{
Prepare taxon_authority for printing. Makes color optional.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{one_of}
\alias{one_of}
\title{dplyr select_helpers}
\description{
dplyr select_helpers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{data}
\name{database_ref}
\alias{database_ref}
\title{Database list}
\format{
An object of class \code{taxa_taxon_db_def} (inherits from \code{vctrs_rcrd}, \code{vctrs_vctr}) of length 8.
}
\usage{
database_ref
}
\description{
The list of known databases.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{limited_print}
\alias{limited_print}
\title{Print a subset of a character vector}
\usage{
limited_print(
  chars,
  prefix = "",
  sep = ", ",
  mid = " ... ",
  trunc_char = "[truncated]",
  max_chars = getOption("width") - nchar(prefix) - 5,
  type = "message"
)
}
\arguments{
\item{chars}{(\code{character}) What to print.}

\item{prefix}{(\code{character} of length 1) What to print before
\code{chars}, on the same line.}

\item{sep}{What to put between consecutive values}

\item{mid}{What is used to indicate omitted values}

\item{max_chars}{(\code{numeric} of length 1) The maximum number of
characters to print.}

\item{type}{(\code{"error"}, \code{"warning"}, \code{"message"}, \code{"cat"}, \code{"print"}, \code{"silent"}, \code{"plain"})}

\item{trunc}{What is appended onto truncated values}
}
\value{
\code{NULL}
}
\description{
Prints the start and end values for a character vector. The number of values
printed depend on the width of the screen by default.
}
\examples{
taxa:::limited_print(1:100)
taxa:::limited_print(1:10000)
taxa:::limited_print(1:10000, prefix = "stuff:")

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{is_root}
\alias{is_root}
\title{Test if taxa are roots}
\usage{
is_root(x, subset = NULL)
}
\arguments{
\item{x}{An object containing taxonomic relationships, such as \link{taxonomy} objects.}

\item{subset}{The subset of the tree to search for roots to that subset. Can be indexes or names.}
}
\description{
Check if each taxon is a root. A root is a taxon with no supertaxon.
}
\examples{
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))
is_root(x)
is_root(x, subset = 2:8)

}
\seealso{
Other root functions: 
\code{\link{roots}()}
}
\concept{root functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{everything}
\alias{everything}
\title{dplyr select_helpers}
\description{
dplyr select_helpers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank_level.R
\name{taxon_rank_level}
\alias{taxon_rank_level}
\title{Taxon rank level}
\usage{
taxon_rank_level(
  level = character(),
  order = NULL,
  guess_order = TRUE,
  impute_na = FALSE
)
}
\arguments{
\item{level}{Zero or more taxonomic rank names. If a named numeric is
applied, the names are used for levels and the numeric values are used
for the order. Inputs will be transformed to a \code{character} vector.}

\item{order}{Integers that determine the relative order of taxonomic levels.
Inputs will be transformed to a \code{integer} vector. \code{NA}s can be used to
indicate that the order is not known.}

\item{guess_order}{If \code{TRUE} and no order is given, try to guess order based on rank names.}

\item{impute_na}{If \code{TRUE}, fill in NAs based on nearby values (assumed in ascending order).}
}
\value{
An \code{S3} object of class \code{taxa_taxon_rank_level}
}
\description{
Used to store taxon rank level information. This is used in \code{\link[=taxon_rank]{taxon_rank()}} objects.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printing_style.R
\name{font_punct}
\alias{font_punct}
\title{Punctuation formatting in print methods}
\usage{
font_punct(text)
}
\arguments{
\item{text}{What to print}
}
\description{
A simple wrapper to make changing the formatting of text printed easier.
This is used for non-data, formatting characters
}
\seealso{
Other printer fonts: 
\code{\link{font_default}()},
\code{\link{font_na}()},
\code{\link{font_secondary}()},
\code{\link{font_tax_name}()}
}
\concept{printer fonts}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R, R/documentation.R, R/taxon.R,
%   R/taxon_authority.R, R/taxon_db.R, R/taxon_db_def.R, R/taxon_id.R,
%   R/taxon_rank.R, R/taxon_rank_level.R, R/taxonomy.R
\name{format.taxa_classification}
\alias{format.taxa_classification}
\alias{obj_print_data.taxa_classification}
\alias{obj_print_footer.taxa_classification}
\alias{vec_ptype_abbr.taxa_classification}
\alias{vec_ptype_full.taxa_classification}
\alias{pillar_shaft.taxa_classification}
\alias{taxa_printing_funcs}
\alias{format.taxa_taxon}
\alias{obj_print_data.taxa_taxon}
\alias{obj_print_footer.taxa_taxon}
\alias{vec_ptype_abbr.taxa_taxon}
\alias{vec_ptype_full.taxa_taxon}
\alias{pillar_shaft.taxa_taxon}
\alias{format.taxa_taxon_authority}
\alias{obj_print_data.taxa_taxon_authority}
\alias{vec_ptype_abbr.taxa_taxon_authority}
\alias{vec_ptype_full.taxa_taxon_authority}
\alias{pillar_shaft.taxa_taxon_authority}
\alias{format.taxa_taxon_db}
\alias{vec_ptype_abbr.taxa_taxon_db}
\alias{vec_ptype_full.taxa_taxon_db}
\alias{obj_print_data.taxa_taxon_db_def}
\alias{vec_ptype_full.taxa_taxon_db_def}
\alias{format.taxa_taxon_id}
\alias{obj_print_data.taxa_taxon_id}
\alias{vec_ptype_abbr.taxa_taxon_id}
\alias{vec_ptype_full.taxa_taxon_id}
\alias{pillar_shaft.taxa_taxon_id}
\alias{format.taxa_taxon_rank}
\alias{obj_print_data.taxa_taxon_rank}
\alias{obj_print_footer.taxa_taxon_rank}
\alias{vec_ptype_abbr.taxa_taxon_rank}
\alias{vec_ptype_full.taxa_taxon_rank}
\alias{pillar_shaft.taxa_taxon_rank}
\alias{format.taxa_taxon_rank_level}
\alias{obj_print_data.taxa_taxon_rank_level}
\alias{vec_ptype_abbr.taxa_taxon_rank_level}
\alias{vec_ptype_full.taxa_taxon_rank_level}
\alias{toString.taxa_taxon_rank_level}
\alias{format.taxa_taxonomy}
\alias{obj_print_data.taxa_taxonomy}
\alias{obj_print_footer.taxa_taxonomy}
\alias{vec_ptype_abbr.taxa_taxonomy}
\alias{vec_ptype_full.taxa_taxonomy}
\alias{pillar_shaft.taxa_taxonomy}
\title{taxa printing functions}
\usage{
\method{format}{taxa_classification}(x, ...)

\method{obj_print_data}{taxa_classification}(x, ...)

\method{obj_print_footer}{taxa_classification}(x, ...)

\method{vec_ptype_abbr}{taxa_classification}(x, ...)

\method{vec_ptype_full}{taxa_classification}(x, ...)

\method{pillar_shaft}{taxa_classification}(x, ...)

\method{format}{taxa_taxon}(x, ...)

\method{obj_print_data}{taxa_taxon}(x, ...)

\method{obj_print_footer}{taxa_taxon}(x, ...)

\method{vec_ptype_abbr}{taxa_taxon}(x, ...)

\method{vec_ptype_full}{taxa_taxon}(x, ...)

\method{pillar_shaft}{taxa_taxon}(x, ...)

\method{format}{taxa_taxon_authority}(x, ...)

\method{obj_print_data}{taxa_taxon_authority}(x, ...)

\method{vec_ptype_abbr}{taxa_taxon_authority}(x, ...)

\method{vec_ptype_full}{taxa_taxon_authority}(x, ...)

\method{pillar_shaft}{taxa_taxon_authority}(x, ...)

\method{format}{taxa_taxon_db}(x, ...)

\method{vec_ptype_abbr}{taxa_taxon_db}(x, ...)

\method{vec_ptype_full}{taxa_taxon_db}(x, ...)

\method{obj_print_data}{taxa_taxon_db_def}(x, ...)

\method{vec_ptype_full}{taxa_taxon_db_def}(x, ...)

\method{format}{taxa_taxon_id}(x, ...)

\method{obj_print_data}{taxa_taxon_id}(x, ...)

\method{vec_ptype_abbr}{taxa_taxon_id}(x, ...)

\method{vec_ptype_full}{taxa_taxon_id}(x, ...)

\method{pillar_shaft}{taxa_taxon_id}(x, ...)

\method{format}{taxa_taxon_rank}(x, ...)

\method{obj_print_data}{taxa_taxon_rank}(x, ...)

\method{obj_print_footer}{taxa_taxon_rank}(x, ...)

\method{vec_ptype_abbr}{taxa_taxon_rank}(x, ...)

\method{vec_ptype_full}{taxa_taxon_rank}(x, ...)

\method{pillar_shaft}{taxa_taxon_rank}(x, ...)

\method{format}{taxa_taxon_rank_level}(x, ...)

\method{obj_print_data}{taxa_taxon_rank_level}(x, ...)

\method{vec_ptype_abbr}{taxa_taxon_rank_level}(x, ...)

\method{vec_ptype_full}{taxa_taxon_rank_level}(x, ...)

\method{toString}{taxa_taxon_rank_level}(x, ...)

\method{format}{taxa_taxonomy}(x, ...)

\method{obj_print_data}{taxa_taxonomy}(x, ...)

\method{obj_print_footer}{taxa_taxonomy}(x, ...)

\method{vec_ptype_abbr}{taxa_taxonomy}(x, ...)

\method{vec_ptype_full}{taxa_taxonomy}(x, ...)

\method{pillar_shaft}{taxa_taxonomy}(x, ...)
}
\description{
Functions used internally for printing information in taxon objects. They have to be exported to
work, but they are not intended to be directly used by most users.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{matches}
\alias{matches}
\title{dplyr select_helpers}
\description{
dplyr select_helpers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_db.R
\name{is_taxon_db}
\alias{is_taxon_db}
\title{Check if something is a \link{taxon_db}}
\usage{
is_taxon_db(x)
}
\arguments{
\item{x}{An object to test}
}
\description{
Check if an object is of the \link{taxon_db} class
}
\examples{
x <- taxon_db(c('ncbi', 'ncbi', 'itis'))
is_taxon_db(x)
is_taxon_db(1:3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxonomy.R
\name{tax_id.taxa_classification}
\alias{tax_id.taxa_classification}
\alias{tax_id<-.taxa_classification}
\alias{tax_id}
\alias{tax_id<-}
\alias{tax_id.taxa_taxon}
\alias{tax_id<-.taxa_taxon}
\alias{tax_id.taxa_taxonomy}
\alias{tax_id<-.taxa_taxonomy}
\title{Set and get taxon IDs}
\usage{
\method{tax_id}{taxa_classification}(x)

\method{tax_id}{taxa_classification}(x) <- value

tax_id(x)

tax_id(x) <- value

\method{tax_id}{taxa_taxon}(x)

\method{tax_id}{taxa_taxon}(x) <- value

\method{tax_id}{taxa_taxonomy}(x)

\method{tax_id}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon IDs.}

\item{value}{The taxon IDs to set. Inputs will be coerced into a \link{taxon_id} vector.}
}
\description{
Set and get the taxon IDs in objects that have them, such as \link{taxon} objects.
}
\examples{
x <- taxon(name = c('Homo sapiens', 'Bacillus', 'Ascomycota', 'Ericaceae'),
           rank = c('species', 'genus', 'phylum', 'family'),
           id = taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi'),
           auth = c('Linnaeus, 1758', 'Cohn 1872', NA, 'Juss., 1789'))

tax_id(x)
tax_id(x) <- paste0('00', tax_id(x))
tax_id(x)[1] <- '00000'

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{new_classification}
\alias{new_classification}
\title{Minimal classfication constructor}
\usage{
new_classification(taxonomy = taxonomy(), instances = integer())
}
\arguments{
\item{taxonomy}{A \code{\link[=taxonomy]{taxonomy()}} object.}

\item{instances}{The indexes of each instance of a taxon in the taxonomy. Can be any length.}
}
\value{
An \code{S3} object of class \code{taxa_classification}
}
\description{
Minimal classfication constructor for internal use. Only use when the input is known to be valid
since few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxonomy.R
\name{tax_name.taxa_classification}
\alias{tax_name.taxa_classification}
\alias{tax_name<-.taxa_classification}
\alias{tax_name}
\alias{tax_name<-}
\alias{tax_name.taxa_taxon}
\alias{tax_name<-.taxa_taxon}
\alias{tax_name.taxa_taxonomy}
\alias{tax_name<-.taxa_taxonomy}
\title{Set and get taxon names}
\usage{
\method{tax_name}{taxa_classification}(x)

\method{tax_name}{taxa_classification}(x) <- value

tax_name(x)

tax_name(x) <- value

\method{tax_name}{taxa_taxon}(x)

\method{tax_name}{taxa_taxon}(x) <- value

\method{tax_name}{taxa_taxonomy}(x)

\method{tax_name}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon names.}

\item{value}{The taxon names to set. Inputs will be coerced into a \link{character} vector.}
}
\description{
Set and get the taxon names in objects that have them, such as \link{taxon} objects.
Note that this is not the same as adding vector names with \link{names}.
}
\examples{
x <- taxon(name = c('Homo sapiens', 'Bacillus', 'Ascomycota', 'Ericaceae'),
           rank = c('species', 'genus', 'phylum', 'family'),
           id = taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi'),
           auth = c('Linnaeus, 1758', 'Cohn 1872', NA, 'Juss., 1789'))

tax_name(x)
tax_name(x) <- tolower(tax_name(x))
tax_name(x)[1] <- 'Billy'

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank_level.R
\name{printed_taxon_rank_level}
\alias{printed_taxon_rank_level}
\title{Prepare taxon_rank_level for printing}
\usage{
printed_taxon_rank_level(x, color = FALSE)
}
\arguments{
\item{color}{Use color?}
}
\value{
character
}
\description{
Prepare taxon_rank_level for printing. Makes color optional.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printing_style.R
\name{font_default}
\alias{font_default}
\title{Default font}
\usage{
font_default(text)
}
\arguments{
\item{text}{What to print.}
}
\description{
A wrapper to make changing the formatting of text printed easier.
}
\seealso{
Other printer fonts: 
\code{\link{font_na}()},
\code{\link{font_punct}()},
\code{\link{font_secondary}()},
\code{\link{font_tax_name}()}
}
\concept{printer fonts}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_db.R
\name{taxa_taxon_db-class}
\alias{taxa_taxon_db-class}
\alias{taxa_taxon_db}
\title{Taxon database class}
\description{
Taxon database class. See \link{taxon_db} for more information
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_db_def.R
\name{taxon_db_def}
\alias{taxon_db_def}
\title{Taxon database definition class}
\usage{
taxon_db_def(
  name = character(),
  url = NA_character_,
  desc = NA_character_,
  id_regex = NA_character_,
  rank_levels = rep(list(NULL), length(name))
)
}
\arguments{
\item{name}{Name of the database in lower case. Inputs will be transformed to a \code{character} vector.}

\item{url}{URL of the database website. Inputs will be transformed to a \code{character} vector.}

\item{desc}{Description of the database. Inputs will be transformed to a \code{character} vector.}

\item{id_regex}{A regular expression for taxon IDs of the database. Inputs will be transformed to a \code{character} vector.}

\item{rank_levels}{Valid taxonomic ranks for the database. Should be a list of \code{numeric} vectors named by taxonomic ranks.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_db_def}
}
\description{
Used to store information on taxonomic databases that is used to validate information in other classes.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank_level.R
\name{new_taxon_rank_level}
\alias{new_taxon_rank_level}
\title{Minimal taxon_rank_level constructor}
\usage{
new_taxon_rank_level(level = character(), order = numeric())
}
\arguments{
\item{level}{Zero or more taxonomic rank names. If a named numeric is
applied, the names are used for levels and the numeric values are used
for the order. Inputs will be transformed to a \code{character} vector.}

\item{order}{Integers that determine the relative order of taxonomic levels.
Inputs will be transformed to a \code{integer} vector. \code{NA}s can be used to
indicate that the order is not known.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_rank_level}
}
\description{
Minimal taxon_rank_level constructor for internal use. Only use when the
input is known to be valid since few validity checks are done.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxon_authority.R, R/taxonomy.R
\name{tax_cite.taxa_classification}
\alias{tax_cite.taxa_classification}
\alias{tax_cite<-.taxa_classification}
\alias{tax_cite}
\alias{tax_cite<-}
\alias{tax_cite.taxa_taxon}
\alias{tax_cite<-.taxa_taxon}
\alias{tax_cite<-.taxa_taxon_authority}
\alias{tax_cite.taxa_taxon_authority}
\alias{tax_cite.taxa_taxonomy}
\alias{tax_cite<-.taxa_taxonomy}
\title{Set and get taxon authority citations}
\usage{
\method{tax_cite}{taxa_classification}(x)

\method{tax_cite}{taxa_classification}(x) <- value

tax_cite(x)

tax_cite(x) <- value

\method{tax_cite}{taxa_taxon}(x)

\method{tax_cite}{taxa_taxon}(x) <- value

\method{tax_cite}{taxa_taxon_authority}(x) <- value

\method{tax_cite}{taxa_taxon_authority}(x)

\method{tax_cite}{taxa_taxonomy}(x)

\method{tax_cite}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon authority dates.}

\item{value}{The taxon citations to set. Inputs will be coerced into a \link{taxon_authority} vector.}
}
\description{
Set and get the taxon authority citations in objects that have them, such as \link{taxon_authority} objects.
}
\examples{
x <- taxon_authority(c('Cham. & Schldl.', 'L.'),
                     date = c('1827', '1753'),
                     citation = c(NA, 'Species Plantarum'))
tax_cite(x)
tax_cite(x)[1] <- "Cham. et al 1984"

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{classification}
\alias{classification}
\title{Taxon class}
\usage{
classification(x = NULL, taxonomy = NULL, .names = NULL)
}
\arguments{
\item{x}{One of:
\itemize{
\item A list where each item represents a series of nested taxa. The contents of
the list can be in any form that can be converted to a \link{taxon} vector.
\item The indexes/names of each instance of a taxon in a \link{taxonomy} object specified by the \code{taxonomy} option. Can
be any length, but must consist of valid indexes for taxa in the \code{taxonomy}
object.
}}

\item{taxonomy}{A \link{taxonomy} object. Only needed if taxon indexes are supplied as the first argument.}

\item{.names}{The names of the vector.}
}
\value{
An \code{S3} object of class \code{taxa_classification}
}
\description{
\Sexpr[results=rd, stage=render]{taxa:::lifecycle("experimental")} Used to
store classifications in reference to a taxonomic tree.
}
\examples{

# Create classification vector with a list
x <- classification(list(
  c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo'),
  c('Carnivora', 'Felidae', 'Panthera', 'Panthera tigris'),
  c('Carnivora', 'Ursidae', 'Ursus', 'Ursus arctos'),
  c('Carnivora', 'Ursidae', 'Ursus', 'Ursus arctos'),
  c('Carnivora', 'Felidae', 'Panthera', 'Panthera tigris')
))


# Create classification vector with indexes and a taxonomy
x <- classification(c(3, 4, 4, 5, 5, 6, 8, 8, 2, 5, 6, 2),
                    taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                               'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
                             supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7)))

x <- classification(c(3, 4, 4, 5, 5, 6, 8, 8, 2, 5, 6, 2),
                    taxonomy(taxon(name = c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                                            'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
                                   rank = c('order', 'family', 'genus', 'species',
                                            'species', 'family', 'genus', 'species'),
                                   id = taxon_id(c('33554', '9681', '9688', '9689',
                                                   '9694', '9632', '9639', '9644'),
                                                 db = 'ncbi'),
                                   auth = c('Bowdich, 1821', 'Fischer, 1817',
                                            'Oken, 1816', 'L., 1758',
                                            'L., 1758', 'Fischer, 1817',
                                            'L., 1758', 'L., 1758')),
                             supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7)))
names(x) <- letters[1:12]


# Get parts of the classification vector
tax_name(x)
tax_rank(x)
tax_id(x)
tax_db(x)
tax_auth(x)
tax_author(x)
tax_date(x)
tax_cite(x)

# Manipulate classification vectors
x[1:3]
x[tax_rank(x) > 'family']
# c(x, x)
# x['b'] <- NA
is.na(x)
# as.data.frame(x)
# tibble::as_tibble(x)

# Use as columns in tables
tibble::tibble(x = x, y = 1:12)
data.frame(x = x, y = 1:12)

}
\seealso{
Other classes: 
\code{\link{[.taxa_classification}()},
\code{\link{taxon_authority}()},
\code{\link{taxon_db}()},
\code{\link{taxon_id}()},
\code{\link{taxon_rank}()},
\code{\link{taxon}()}
}
\concept{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printing_style.R
\name{font_na}
\alias{font_na}
\title{Font for NAs in print methods}
\usage{
font_na(text)
}
\arguments{
\item{text}{What to print}
}
\description{
A simple wrapper to make changing the formatting of text printed easier.
This is used for \code{NA}s
}
\seealso{
Other printer fonts: 
\code{\link{font_default}()},
\code{\link{font_punct}()},
\code{\link{font_secondary}()},
\code{\link{font_tax_name}()}
}
\concept{printer fonts}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{data}
\name{rank_ref}
\alias{rank_ref}
\title{All known taxonomic ranks}
\format{
An object of class \code{numeric} of length 52.
}
\usage{
rank_ref
}
\description{
A list of taxonomic ranks from all databases used combined into a single
vector to make it easier to maintain the relative order of ranks when data
from multiple databases are combined.
}
\section{Attribution}{


This list was adapted from a similar one in the \code{taxize} package.
}

\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxon_authority.R, R/taxonomy.R
\name{tax_date.taxa_classification}
\alias{tax_date.taxa_classification}
\alias{tax_date<-.taxa_classification}
\alias{tax_date}
\alias{tax_date<-}
\alias{tax_date.taxa_taxon}
\alias{tax_date<-.taxa_taxon}
\alias{tax_date<-.taxa_taxon_authority}
\alias{tax_date.taxa_taxon_authority}
\alias{tax_date.taxa_taxonomy}
\alias{tax_date<-.taxa_taxonomy}
\title{Set and get taxon authority dates}
\usage{
\method{tax_date}{taxa_classification}(x)

\method{tax_date}{taxa_classification}(x) <- value

tax_date(x)

tax_date(x) <- value

\method{tax_date}{taxa_taxon}(x)

\method{tax_date}{taxa_taxon}(x) <- value

\method{tax_date}{taxa_taxon_authority}(x) <- value

\method{tax_date}{taxa_taxon_authority}(x)

\method{tax_date}{taxa_taxonomy}(x)

\method{tax_date}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon authority dates.}

\item{value}{The taxon authority dates to set. Inputs will be coerced into a \link{character} vector.}
}
\description{
Set and get the taxon authority dates in objects that have them, such as \link{taxon_authority} objects.
}
\examples{
x <- taxon_authority(c('Cham. & Schldl.', 'L.'),
                     date = c('1827', '1753'))
tax_date(x)
tax_date(x)[1] <- "1984"
tax_date(x) <- c(NA, '1800')

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_id.R
\name{taxon_id}
\alias{taxon_id}
\title{Taxon ID class}
\usage{
taxon_id(id = character(), db = NA, .names = NULL)
}
\arguments{
\item{id}{Zero or more taxonomic ids. Inputs will be transformed to a \link{character} vector if
possible.}

\item{db}{The name(s) of the database(s) associated with the IDs. If not \code{NA} (the default), the
input must consist of names of databases in \link[=db_ref]{db_ref$get()}.}

\item{.names}{The names that will be applied to the vector.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_id}
}
\description{
\Sexpr[results=rd, stage=render]{taxa:::lifecycle("maturing")}
Used to store the ID corresponding to taxa, either arbitrary or from a
particular taxonomy database. This is typically used to store taxon IDs in
\link{taxon} objects.
}
\examples{

# Making new objects
x <- taxon_id(c('A', 'B', 'C'))
x <- taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi')
x <- taxon_id(c('9606', '1386', '4890', '4345'),
              db = c('ncbi', 'ncbi', 'itis', 'itis'))
names(x) <- c('a', 'b', 'c', 'd')

# Manipulating objects
as.character(x)
x[2:3]
x[2:3] <- 'ABC'
x[c('a', 'c')] <- '123'
x[['b']] <- taxon_id('123423', db = 'ncbi')
tax_db(x)
tax_db(x) <- 'nbn'
c(x, x)

# Using as columns in tables
tibble::tibble(x = x, y = 1:4)
data.frame(x = x, y = 1:4)

# Convert to tables
tibble::as_tibble(x)
as_data_frame(x)

# Trying to use an invalid ID with a specified database causes an error
#taxon_id('NOLETTERS', db = 'ncbi')

}
\seealso{
Other classes: 
\code{\link{[.taxa_classification}()},
\code{\link{classification}()},
\code{\link{taxon_authority}()},
\code{\link{taxon_db}()},
\code{\link{taxon_rank}()},
\code{\link{taxon}()}
}
\concept{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{unname_fields}
\alias{unname_fields}
\title{Remove names from fields in a vctrs rcrd}
\usage{
unname_fields(x)
}
\arguments{
\item{x}{a vctrs rcrd}
}
\description{
Remove names from fields in a vctrs rcrd
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxa_taxonomy-class}
\alias{taxa_taxonomy-class}
\alias{taxa_taxonomy}
\title{Taxonomy class}
\description{
Taxonomy class. See \link{taxonomy} for more information
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{\%>\%}
\alias{\%>\%}
\title{magrittr forward-pipe operator}
\description{
magrittr forward-pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printing_style.R
\name{interleave}
\alias{interleave}
\title{Interleaves two vectors}
\usage{
interleave(v1, v2)
}
\description{
Taken from "http://r.789695.n4.nabble.com/Interleaving-elements-of-two-vectors-td795123.html"
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank_level.R
\name{check_taxon_rank_order}
\alias{check_taxon_rank_order}
\title{Check that order is ascending}
\usage{
check_taxon_rank_order(level, order, warn = FALSE)
}
\arguments{
\item{level}{Zero or more taxonomic rank names.}

\item{order}{Integers that determine the relative order of taxonomic levels.}

\item{warn}{If \code{TRUE}, issue a warning when not in ascending order.}
}
\description{
Check that order is ascending and reorder the orders and their levels if needed.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxonomy.R
\name{tax_auth.taxa_classification}
\alias{tax_auth.taxa_classification}
\alias{tax_auth<-.taxa_classification}
\alias{tax_auth}
\alias{tax_auth<-}
\alias{tax_auth.taxa_taxon}
\alias{tax_auth<-.taxa_taxon}
\alias{tax_auth.taxa_taxonomy}
\alias{tax_auth<-.taxa_taxonomy}
\title{Set and get taxon authorities}
\usage{
\method{tax_auth}{taxa_classification}(x)

\method{tax_auth}{taxa_classification}(x) <- value

tax_auth(x)

tax_auth(x) <- value

\method{tax_auth}{taxa_taxon}(x)

\method{tax_auth}{taxa_taxon}(x) <- value

\method{tax_auth}{taxa_taxonomy}(x)

\method{tax_auth}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon authorities.}

\item{value}{The taxon IDs to set. Inputs will be coerced into a \link{taxon_id} vector.}
}
\description{
Set and get the taxon authorities in objects that have them, such as \link{taxon} objects.
Note that this sets all the authority information, such as author name, date, and citations.
To set or get just one of part of the authorities, use \link{tax_author}, \link{tax_date}, or \link{tax_cite} instead.
}
\examples{
x <- taxon(name = c('Homo sapiens', 'Bacillus', 'Ascomycota', 'Ericaceae'),
           rank = c('species', 'genus', 'phylum', 'family'),
           id = taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi'),
           auth = c('Linnaeus, 1758', 'Cohn 1872', NA, 'Juss., 1789'))

tax_auth(x)
tax_auth(x) <- tolower(tax_auth(x))
tax_auth(x)[1] <- 'Billy'

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_generics.R
\name{\%in\%}
\alias{\%in\%}
\title{Value matching for taxa package}
\usage{
x \%in\% table
}
\arguments{
\item{x}{vector or \code{NULL}: the values to be matched.
    \link[base]{Long vectors} are supported.}

\item{table}{vector or \code{NULL}: the values to be matched against.
    \link[base]{Long vectors} are not supported.}
}
\description{
A wrapper for the base value matching \%in\% that is used to take into consideration features of the taxa package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_rank.R
\name{taxon_rank}
\alias{taxon_rank}
\title{Taxon rank class}
\usage{
taxon_rank(
  rank = character(),
  .names = NULL,
  levels = NULL,
  guess_order = TRUE
)
}
\arguments{
\item{rank}{Zero or more taxonomic rank names. Inputs will be transformed to a \link{character}
vector.}

\item{.names}{The names of the vector}

\item{levels}{A named numeric vector indicating the names and orders of possible taxonomic ranks.
Higher numbers indicate for fine-scale groupings. Ranks of unknown order can be indicated with
\code{NA} instead of a number.}

\item{guess_order}{If \code{TRUE} and no rank order is given using numbers, try to guess order based
on rank names.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_rank}
}
\description{
\Sexpr[results=rd, stage=render]{taxa:::lifecycle("maturing")}
Used to store taxon ranks, possibly associated with a taxonomy database. This is typically used to
store taxon ranks in \link{taxon} objects.
}
\examples{

# Making new objects
x <- taxon_rank(c('species', 'species', 'phylum', 'family'))

# Specifiying level order
taxon_rank(c('A', 'B', 'C', 'D', 'A', 'D', 'D'),
           levels = c('D', 'C', 'B', 'A'))
taxon_rank(c('A', 'B', 'C', 'D', 'A', 'D', 'D'),
           levels = c(D = NA, A = 10, B = 20, C = 30))
names(x) <- c('a', 'b', 'c', 'd')

# Manipulating objects
as.character(x)
as.factor(x)
as.ordered(x)
x[2:3]
x[x > 'family'] <- taxon_rank('unknown')
x[1] <- taxon_rank('order')
x['b']
x['b'] <- 'order'

# Using as columns in tables
tibble::tibble(x = x, y = 1:4)
data.frame(x = x, y = 1:4)

# Converting to tables
tibble::as_tibble(x)
as_data_frame(x)

# Trying to add an unknown level as a character causes an error
#x[2] <- 'superkingdom'

# But you can add a new level using taxon_rank objects
x[2] <- taxon_rank('superkingdom')

}
\seealso{
Other classes: 
\code{\link{[.taxa_classification}()},
\code{\link{classification}()},
\code{\link{taxon_authority}()},
\code{\link{taxon_db}()},
\code{\link{taxon_id}()},
\code{\link{taxon}()}
}
\concept{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_db.R
\name{taxon_db}
\alias{taxon_db}
\title{Taxon database class}
\usage{
taxon_db(db = character(), .names = NULL, ...)
}
\arguments{
\item{db}{Zero or more taxonomic database names. Should be a name contained in
\link{db_ref}. Inputs will be transformed to a \link{character} vector if possible.}

\item{.names}{The names of the vector.}

\item{...}{Additional arguments.}
}
\value{
An \code{S3} object of class \code{taxa_taxon_db}
}
\description{
\Sexpr[results=rd, stage=render]{taxa:::lifecycle("maturing")}
Used to store the names of taxon databases defined in \link{db_ref}. Primarily
used in other classes like \link{taxon_id} to define databases for each item.
}
\examples{

# Making new objects
x <- taxon_db(c('ncbi', 'ncbi', 'itis'))
x

# Manipulating objects
as.character(x)
x[2:3]
x[2:3] <- 'nbn'
names(x) <- c('a', 'b', 'c')
x['b']
x['b'] <- 'nbn'
x[x == 'itis'] <- 'gbif'

# Using as columns in tables
tibble::tibble(x = x, y = 1:3)
data.frame(x = x, y = 1:3)

# Converting to tables
tibble::as_tibble(x)
as_data_frame(x)

# Trying to use an invalid database generates an error
# x <- taxon_db(c('ncbi', 'ncbi', 'my_custom_db'))
# x[x == 'itis'] <- 'my_custom_db'

# Listing known databases and their properties
db_ref$get()

# Adding and using a new database
db_ref$set(name = 'my_custom_db', desc = 'I just made this up')
db_ref$get()
x <- taxon_db(c('ncbi', 'ncbi', 'my_custom_db'))

}
\seealso{
Other classes: 
\code{\link{[.taxa_classification}()},
\code{\link{classification}()},
\code{\link{taxon_authority}()},
\code{\link{taxon_id}()},
\code{\link{taxon_rank}()},
\code{\link{taxon}()}
}
\concept{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{default_db_ref}
\alias{default_db_ref}
\title{Defines valid taxonomic databases}
\usage{
default_db_ref(defaults = list())
}
\arguments{
\item{name}{(character) name of the database}

\item{url}{(character) url for the database}

\item{desc}{(character) description of the database}

\item{id_regex}{(character) id regex}
}
\description{
Defines valid taxonomic databases
}
\section{Attribution}{


This code is copied from the code handling options in \link{knitr}.
}

\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{starts_with}
\alias{starts_with}
\title{dplyr select_helpers}
\description{
dplyr select_helpers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R,
%   R/getter_setter_s3_generics.R, R/taxon.R, R/taxonomy.R
\name{tax_rank.taxa_classification}
\alias{tax_rank.taxa_classification}
\alias{tax_rank<-.taxa_classification}
\alias{tax_rank}
\alias{tax_rank<-}
\alias{tax_rank.taxa_taxon}
\alias{tax_rank<-.taxa_taxon}
\alias{tax_rank.taxa_taxonomy}
\alias{tax_rank<-.taxa_taxonomy}
\title{Set and get taxon ranks}
\usage{
\method{tax_rank}{taxa_classification}(x)

\method{tax_rank}{taxa_classification}(x) <- value

tax_rank(x)

tax_rank(x) <- value

\method{tax_rank}{taxa_taxon}(x)

\method{tax_rank}{taxa_taxon}(x) <- value

\method{tax_rank}{taxa_taxonomy}(x)

\method{tax_rank}{taxa_taxonomy}(x) <- value
}
\arguments{
\item{x}{An object with taxon ranks.}

\item{value}{The taxon ranks to set. Inputs will be coerced into a \link{taxon_rank} vector.}
}
\description{
Set and get the taxon ranks in objects that have them, such as \link{taxon} objects.
}
\examples{
x <- taxon(name = c('Homo sapiens', 'Bacillus', 'Ascomycota', 'Ericaceae'),
           rank = c('species', 'genus', 'phylum', 'family'),
           id = taxon_id(c('9606', '1386', '4890', '4345'), db = 'ncbi'),
           auth = c('Linnaeus, 1758', 'Cohn 1872', NA, 'Juss., 1789'))

tax_rank(x)
tax_rank(x) <- 'species'
tax_rank(x)[1] <- taxon_rank('family')

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{supertaxa}
\alias{supertaxa}
\title{Get supertaxa}
\usage{
supertaxa(
  x,
  subset = NULL,
  max_depth = NULL,
  include = FALSE,
  value = NULL,
  use_na = FALSE,
  ...
)
}
\arguments{
\item{x}{The object to get supertaxa for, such as a \link{taxonomy} object.}

\item{subset}{The subset of the tree to search for roots to that subset. Can
be indexes or names.}

\item{max_depth}{The number of levels to traverse. For example, \code{max_depth = 1} returns only immediate supertaxa. By default (NULL) information for all
supertaxa is returned.}

\item{include}{If \code{TRUE}, include information for each taxon in the output.}

\item{value}{Something to return instead of indexes. Must be the same length
as the number of taxa.}

\item{use_na}{Add a NA to represent the root of the taxonomy (i.e. no
supertaxon)}

\item{...}{Additional arguments.}
}
\description{
Get supertaxa indexes for each taxon or another per-taxon value. Supertaxa
are taxa a taxon is contained in.
}
\examples{
# Generate example data
x <- taxonomy(c('Carnivora', 'Felidae', 'Panthera', 'Panthera leo',
                'Panthera tigris', 'Ursidae', 'Ursus', 'Ursus arctos'),
              supertaxa = c(NA, 1, 2, 3, 3, 1, 6, 7))

# The indexes of all supertaxa (with supertaxa of supertaxa, etc) for each taxon
supertaxa(x)

# Return something other than index
supertaxa(x, value = tax_name(x))

# Include each taxon with its supertaxa
supertaxa(x, value = tax_name(x), include = TRUE)

# Only return data for some taxa (faster than subsetting the whole result)
supertaxa(x, subset = 3)

}
\seealso{
Other taxonomy functions: 
\code{\link{internodes}()},
\code{\link{leaves}()},
\code{\link{roots}()},
\code{\link{stems}()},
\code{\link{subtaxa}()}

Other supertaxa functions: 
\code{\link{n_supertaxa}()}
}
\concept{supertaxa functions}
\concept{taxonomy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printing_style.R
\name{print_with_color}
\alias{print_with_color}
\title{Print that works with color}
\usage{
print_with_color(x, original_length = length(x), ...)
}
\arguments{
\item{x}{What to print, typically a character vector}

\item{original_length}{The length of the full vector if only part was given.}

\item{...}{Passed to \code{print}}
}
\description{
The same as the \code{print} function, but can print colored text. Its a bit of a hack, but the only
way I found to replicate the behavior of \code{print} without rewriting the entire \code{print} function.
}
\keyword{internal}
