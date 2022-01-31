[![CRAN status](https://www.r-pkg.org/badges/version/iconr)](https://CRAN.R-project.org/package=iconr)
[![status](https://joss.theoj.org/papers/e68e041e66a613918f76bf43db3f8b02/status.svg)](https://joss.theoj.org/papers/e68e041e66a613918f76bf43db3f8b02)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4767529.svg)](https://doi.org/10.5281/zenodo.4767529)
[![R-CMD-check](https://github.com/zoometh/iconr/workflows/R-CMD-check/badge.svg)](https://github.com/zoometh/iconr/actions) [![Build Status](https://travis-ci.org/zoometh/iconr.svg?branch=master)](https://travis-ci.org/zoometh/iconr)
          
# ***iconr*** package <br> modeling Prehistoric iconography <img src="logo/iconr_logo.png" width='100px' align="right"/>
> Created by [Thomas Huet](mailto:thomashuet7@gmail.com), [Jose M Pozo](mailto:josmpozo@gmail.com), [Craig Alexander](mailto:craiga304@gmail.com)
  
The R package ***iconr*** is grounded in graph theory, spatial analysis ([composition analysis](#iconr-stable-version-the-analysis-of-compositions)), and shape analysis ([geometric morphometric](#iconr-development-version-the-analysis-of-compositions-and-geometric-morphometry)) to offer concepts and functions for modeling Prehistoric iconographic compositions and preparing for further analysis (clustering, typology tree, Harris diagram, etc.). The package purpose is to contribute to cross-cultural comparison through a greater normalization of quantitative analysis.

The theoretical background is as follows: some objects can have a decoration, and a decoration is composed of graphical units (GUs). The whole decoration is considered as a graph and can be analyzed with Graph Theory and compared to other graphs (ie, decorations). As favored [GIS entry](https://zoometh.github.io/iconr/articles/gis.html), GIS indexes can also be used. 

See: how to [contribute](.github/CONTRIBUTING.md) to the next package release, and how to [report an issue](https://github.com/zoometh/iconr/issues) using the [issue template](.github/ISSUE_TEMPLATE.md).
  
## ***iconr*** stable version: the Analysis of Compositions

The ***iconr*** v. 0.1.0 stable version can be installed from the CRAN

```
install.packages("iconr")
```

The v. 0.1.0 allows the analysis of compositions

<center>
  
![](doc/img/solana_voronoi.png){width=800px}
  
</center>
  

Photograph of Solana de Cabañas (Extremadura, Spain) Late Bronze Age stele[^1] [**1**]. Graphical units (GUs) drawing [**2**]. Each GU is recorded with a vertex (POINT) [**3**], and each contiguous vertex (Voronoi cell) is linked with an edge [**4**, **5**]. 

### Overview of the functions for the analysis of compositions

The ***iconr*** v. 0.1.0 functions' descriptions and examples are available on this [website](https://zoometh.github.io/iconr/articles/index.html). 

#### Plot a decoration

Set data folder, select the decoration to be plotted, read files of nodes, edges, and images and plot a decoration

```
dataDir <- system.file("extdata", package = "iconr")
site <- "Brozas" ; decor <- "Brozas"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
imgs <- read.table(paste0(dataDir, "/imgs.tsv"),
                   sep="\t", stringsAsFactors = FALSE)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor,
              dir = dataDir,
              lbl.size = 0.4,
              nd.var = "type")
```
  
  
<img src="doc/img/plot_dec_graph.png" align="center"/>
  
#### Common edges

Common edges between pairs of decorations allow to measure the similarities in their composition. Common edges are more accurate than common nodes (see also, [nds_compar()](https://zoometh.github.io/iconr/reference/list_compar.html) and [eds_compar()](https://zoometh.github.io/iconr/reference/list_compar.html))
For example, we plot common edges shared by the three first decorations of the training dataset with the [plot_compar()](https://zoometh.github.io/iconr/reference/plot_compar.html) function 

```
imgs <- read.table(system.file("extdata", "imgs.tsv", package = "iconr"),
                   sep="\t", stringsAsFactors = FALSE)
nodes <- read.table(system.file("extdata", "nodes.tsv", package = "iconr"),
                    sep="\t", stringsAsFactors = FALSE)
edges <- read.table(system.file("extdata", "edges.tsv", package = "iconr"),
                    sep="\t", stringsAsFactors = FALSE)
lgrph <- list_dec(imgs, nodes, edges)
g.compar <- list_compar(lgrph, nd.var="type")
plot_compar(g.compar, c(1, 2, 3), 
            focus = "edges",
            nd.size = c(0.5, 1.5),
            ed.width = c(1, 2.5),
            dir = dataDir,
            img.format = "png")
```

  
<img src="doc/img/_compar_eds_1_2.png" align="center"/>
<img src="doc/img/_compar_eds_1_3.png" align="center"/>
<img src="doc/img/_compar_eds_2_3.png" align="center"/>

The same result, but in the form of a coincidence matrix, can be obtained with the function [same_elements()](https://zoometh.github.io/iconr/reference/same_elements.html)

## ***iconr*** development version: the Analysis of Compositions and Geometric Morphometry

The ***iconr*** latest version, or development version v. 0.1.1, handle graphical units with POLYGON geometries to perform Geometric Morphometry measurements (GMM). This development version can be downloaded from GitHub

```
devtools::install_github("zoometh/iconr")
```

### Geometric Morphometry measurements' functions overview

The ***iconr*** v. 0.1.1 functions are named `morph_*` (morphology) and `conv_*` (conversions). Find their documentation directly on R (eg, `?morph_nds_compar`). They are performed with the [R package Momocs](https://momx.github.io/Momocs/)

#### Compare GUs' shapes

The [sample dataset](https://github.com/zoometh/iconr/tree/master/doc/datasets/PPN) is composed on 5 decorated objects, belonging from 4 sites of the Near-East Pre-Pottery Neolithic:


<font size="2" align="left">
<table style="width:100%">
	<tr align="center">
		<td>
				<img src="doc/datasets/PPN/Ain Ghazal/Ain Ghazal.stat_2_gref.jpg" width='150px' />
		</td><td>
				<img src="doc/datasets/PPN/Ain Ghazal/Ain Ghazal.stat_5_gref.jpg" width='150px' />
		</td><td>
				<img src="doc/datasets/PPN/Jericho/Jericho.tete_afe_gref.jpg" width='150px' />
		</td><td>
				<img src="doc/datasets/PPN/Kfar Hahoresh/Kfar Hahoresh.crane_afg.jpg" width='150px' />
		</td><td>
				<img src="doc/datasets/PPN/Qarassa/Qarassa.figurine__wx.jpg" width='150px' />
		</td>
	</tr><tr>
			  <th style="padding:5px">Ain Ghazal, statue 2, cache 2</th>
		    <th style="padding:5px">Ain Ghazal, statue 5, cache 2</th>
		    <th style="padding:5px">Jericho, statue A, cache 195</th>
		    <th style="padding:5px">Kfar Hahoresh, modelled skull</th>
		    <th style="padding:5px">Qarassa, bone wand</th>
	</tr>
	</table>
</font>

The graphical units 'faces' ('*visages*'), 'eyes' ('*oeils*'), and 'mouths' ('*bouches*') are drawn in a GIS

<p align="center">
  <img alt="img-name" src="doc/img/gis_gmm.png" width="700">
</p>


#### Resume the GUs geometries

After downloading the PPNB dataset, set 'PPN' as the current working directory (`setwd("*my_path*/PPN"")`), read and convert the 'nodes.csv' Well-Known Text geometries to JPG, and resume information

```
nodes <- read.csv2("*my_path*/PPN/_out/nodes.csv")
conv_wkt_to_jpg(nodes = nodes)
morph_resume(dataDir = "*my_path*/PPN",
             nodes = nodes)
```
  
<img src="doc/img/1_resume.png" align="center"/>

#### Compare the different types of GUs

Stack the countours of 'faces' ('*visages*'), 'eyes' ('*oeils*'), and 'mouths' ('*bouches*')

```
conv_wkt_to_jpg(nodes = nodes)
nodes <- read.csv2("*my_path*/PPN/_out/nodes.csv")
conv_wkt_to_jpg(nodes = nodes)
morph_resume(dataDir = "*my_path*/PPN",
             nodes = nodes)
```

<p align="center">
  <img alt="img-name" src="doc/img/visage_compar_stack.png" width="350">
  <br>
  <img alt="img-name" src="doc/img/oeil_compar_stack.png" width="350">
  <br>
  <img alt="img-name" src="doc/img/bouche_compar_stack.png" width="350">
  <br>
</p>


## Citation

Use the canonical form to cite the package (`citation("iconr")`):
```
@Manual{Huet21pckg,
  title = {iconr: Graphical and Spatial Analysis for Prehistoric Iconography},
  author = {Thomas Huet and Jose Pozo},
  year = {2021},
  note = {R package version 0.1.0},
  url = {https://CRAN.R-project.org/package=iconr},
}
```

The ***iconr*** v. 0.1.0 package has also been published in the [Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.03191) under this BibTex reference:

```
@article{Huet21joss,
  doi = {10.21105/joss.03191},
  url = {https://doi.org/10.21105/joss.03191},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {61},
  pages = {3191},
  author = {Thomas Huet and Jose M. Pozo and Craig Alexander},
  title = {Analysis of Prehistoric Iconography with the R package iconr},
  journal = {Journal of Open Source Software}
}
```


## Next release

### Typology of GUs

***iconr*** aims to use a hierarchical _thesaurus_ (tree-like) with controlled vocabularies for GUs' typology . Identity between GU name and value must be unique (URL). See for example the [whole typological tree](https://zoometh.github.io/iconr/articles/img/typo_gu_ug.html):

<center>
   
[![](doc/img/typology_gu.png)](https://zoometh.github.io/iconr/articles/img/typo_gu_ug.html)
  
</center>
  
  
Or these different subtrees: [geometric](https://zoometh.github.io/iconr/articles/img/typo_gu_geometrique.html), [figurative](https://zoometh.github.io/iconr/articles/img/typo_gu_figuratif.html), [zoomorphic](https://zoometh.github.io/iconr/articles/img/typo_gu_zoomorphe.html), 
[technomorphic](https://zoometh.github.io/iconr/articles/img/typo_gu_technomorphe.html), or [anthropomorphic](https://zoometh.github.io/iconr/articles/img/typo_gu_anthropomorphe.html). Such a structure should also be used for other fields than the GU type (eg, 'technique'). Multi-linguism equivalences -- starting with English --, metadata insertion (EXIF) and standardization of the vocabulary (Dublin Core, CIDOC-CRM) is needed

### Superimpositions

The *diachronic* edge `->-` allows to register the superimposition. The next ***iconr*** will integrate an on-the-fly function allowing to create Harris matrices of GUs when such an edge exists. For example here, the Ibahernando stele shows a Latin writing overlaping a spear and a shield representations

<p align="center">
  <img alt="img-name" src="doc/img/ibahernando.png" width="500">
</p>

### Magic wand

The selection of a colored continuous range can be done from a POINT coordinates (x, y) overlapping this colored range (ie. a GU displayed in black on white background). The next ***iconr*** will integrate a function allowing to extract automatically the shape *behind* the POINTS

<p align="center">
  <img alt="img-name" src="doc/img/solana_magic_wand.png" width="250">
</p>

[^1]: credits Museo Arqueologico de Madrid
# Version 0.1.1 (17 May 2021)
First release.

# Version 0.1.0 (12 February 2021)
First official CRAN release.
## Resubmission 3

This is the third resubmission, responding to feedback from Gregor Seyer after the second resubmission

* Rd-tags. Missing Rd-tags:
  + labels_shadow.Rd: \value
  + side_plot.Rd: \value
  
Done

* Rd-tags. More about the structure of the output (class) and also what the output means

Done

* Please make sure that you do not change the user's options, `par()` or working directory
  
    `oldpar <- par(no.readonly = TRUE)`  
    `on.exit(par(oldpar))`

Done. These changes have been added before all calls of the `par()` function. The current package version does not contains any changes of user's working directory (`setwd()`)


## Resubmission 2

This is the second resubmission, responding to feedback from Uwe Ligges after the first resubmission

* Write the url in the form <https://hal.archives-ouvertes.fr/hal-02913656>

Done

## Resubmission 1

This is the first resubmission, responding to feedback from Uwe Ligges after the first submission

* Not more than 5 MB for a CRAN package, please.

Done

* Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>?

Done

## Test environments
* local Windows 8 install, R 4.0.3
* win-builder (devel and release)
* ubuntu 16.04 (on travis-ci), R 4.0.3
* macOS 10.15.7, R 3.6.3
* Rhub

## R CMD check results

0 errors | 0 warnings | 0 note

## Downstream dependencies

There are currently no downstream dependencies for this package
# Table of Contents
1. [_iconr_ 0.1.0](#v1)
3. [Collaborations](#col)
4. [Collaborators](#collab)


# 1. _iconr_ 0.1.0

## publishing

### briefing notes

  - [ ] [Bulletin de la Société Préhistorique Française - Actualités (*scientific news*)](https://www.jstor.org/journal/bullsociprehfran)

  - [ ] [INORA newsletter](https://www.icomos.org/en/resources/publicationall/165-articles-en-francais/centre-de-documentation/557-inora-international-newsletter-on-rock-art) (NOT SURE IT STILL PUBLISH)
  
### communications in congress

  - [x] 2021 [CAA congress](https://2021.caaconference.org/sessions/)
    - [x] [S17. Tools for the Revolution: developing packages for scientific programming in archaeology (Standard)](https://2021.caaconference.org/sessions/#17)
  - [x] 2021 [EAA congress](https://eaa.klinkhamergroup.com/eaa2021/)
    - [x] Between Variability and Singularity: Crossing Theoretical, Qualitative and Computer-based Approaches to Types and Typologies in Archaeology [PaM]
    - [x] Human Visual Archives, Globally. Materials, Forms and Meanings of Human Representations in Ancient Times
  - [ ] 2022 [WAC congress](https://www.wac-9.org/sessions/) (TOO EXPENSIVE FEES)
  - [ ] 2023 [AURA Congress](http://www.ifrao.com/the-aura-congress/)


### Papers in scientific reviews

#### IT and Methods

  - [x] 2021 [JOSS](https://github.com/zoometh/iconr/blob/master/paper.md) 
  - [ ] [JSS](https://raw.githubusercontent.com/zoometh/jss_iconr/main/article_rvTH17.pdf)
    - [draft](https://raw.githubusercontent.com/zoometh/jss_iconr/main/article_rvTH17.pdf) 
  
#### Prehistorical issues

  - [ ] [JAMT](https://www.springer.com/journal/10816) 

  
# 2. Collaborations

[Users/collaborators](#collab) should be specialists on any consensual group of decorations (Paleolithic rock-art, Iron Age stelae, etc.). We will: 

a. [inform](#info) them  
b. propose a [master class](#col.mc)  
c. propose different [training(s)](#col.mc)  
d. propose a [congress session/round-table](#col.rt)  
e. work to get a [publication](#pub.col)  

## A. Information {#info}

contact mail and schedule a [master class](#col.mc)

***
**mail subject**:  

R package _iconr_ - An open source computer-based method to study ancient iconography 

**mail body**:  

Dear Colleague,  

We have the pleasure to introduce you the first version of the R package  [***iconr***](https://github.com/zoometh/iconr#readme). 

The ***iconr*** package is dedicated to Prehistoric iconography modeling and analysis. Grounded on graph theory and spatial analysis, it aims to offer concepts and functions for a greater normalization of quantitative analysis, to facilitate cross-cultural comparisons. The main principle of the package is to consider any iconographic composition (here, 'decoration') as a geometric graph of graphical units. Geometric graphs -- also known as planar graphs or spatialized graphs -- allow to model the neighborhood of these graphical units which are the fundamental relationships of visual semiotics.

<p align="center">
  <img alt="img-name" src="docs/man/figures/edges_compar.png" width="600">
  <br>
    <em>same edges identification</em>
</p>


The first version of the package has been recently uploaded to the [CRAN]((https://cran.r-project.org/web/packages/iconr/index.html)). A brief description of the package has been published on the [_Journal of Open Statistical Software_](https://joss.theoj.org/papers/10.21105/joss.03191) (attached here), and online documentation is already available:

  + [package description](https://zoometh.github.io/iconr/articles/index.html)
  + [GIS use for data entry](https://zoometh.github.io/iconr/articles/gis.html)
  + [interactive examples](https://zoometh.github.io/iconr/articles/shiny.html)
  + [training datasets](https://zoometh.github.io/iconr/articles/examples.html)

The next ***iconr*** release will integrate shape analysis of graphical units, tree-like structures for graphical units thesaurus (multi-linguism & shared vocabularies), use of directed acyclic graphs (DAG) to model the graphical units' superimpositions (ie, Harris matrix).
  
To promote the package utilization, we have schedule a first presentation on ZOOM (duration: ca. 45 minutes), the **xx/xx/xx at 16:00 UTC**. During the presentation, We will explain how it works and what are the expected outcomes (40 minutes) and respond to the audience questions (5 minutes). We will be very happy to meet you at this moment. If you are interested to participate, please conserve the following information

```
Topic: iconr R package - masterclass
Time: xxx

Join Zoom Meeting
https://xxx

Meeting ID: xxx
Passcode: xxx
```
  
If you cannot be present, but you are interested by the package or the presentation, thank you for letting us know. The presentation will be recorded and available on a video platform.  

Best regards,

Thomas Huet, LabEx ARCHIMEDE ANR-11-LABX-0032-01   
Jose M Pozo, Independant researcher  
Craig Alexander, Independant researcher  

***

## B. Masterclass

Online masterclass to present the _iconr_ package and schedule [training(s)](#col.mc) (open dates with Doodle)

1. Geometric graph heuristic

2. Graph analysis indexes (degrees, same edges, etc.)

3. Case studies

## C. Training

### C.1. R and RStudio basic knowledge

Tutorial on R and RStudio install and basic functions

### C.2. Data entry through a GIS

1. Create graph decorations on GIS

### C.3. _iconr_ package

Datasets presentations/Training to use the _iconr_ package:

1. Read, plot and compare graph decorations

2. Further analysis

* schedule [congress session/round-table](#col.rt) (open dates with Doodle)

* open a GitHub forum for FAQs and scientific exchanges

## D. round-table

After the [training](#col.tr), people are supposed to have used *iconr* on their favorite dataset (i.e. a selection of decorations). We will propose them to present their results in an oral communication in the frame of a round-table. Depending on the size of the decoration they have proceed, each used could choose between:

1. **3 + decorations**: qualitative analysis
2. **7 +decorations**: semi-quantitative analysis (ie. rank-based, non-parametric tests: Mann-Withney, Spearman, etc.)
3. **30 + decorations**: quantitative analysis (ie. mean-based, parametric tests: Shapiro-Wilk, Student, etc.)

The aim of the round-table is to find the appropriate way to compare heterogeneous decorative contents (by periods, families, themes, techniques, etc.)

<p align="center">
  <img alt="img-name" src="doc/img/famille_thm.png" width="500">
  <br>
    <em>families x themes</em>
</p>

see: [families examples](https://zoometh.github.io/iconr/articles/shiny.html)


## E. publication {#pub.col}

  - [ ] [JAMT](https://www.springer.com/journal/10816) 
  - [ ] [Adoranten](https://www.rockartscandinavia.com/adoranten-vv4.php)

# Collaborators {#collab}

Identify who will be our next collaborators

## By families {#fam}

  - [ ] Palaeolithic cave rock-art: 
    - France: 
      - [ ] [Julien Monney](http://theses.fr/184904846)
    - Spain: 
      - [ ] [Aitor Ruiz-Redondo](https://southampton.academia.edu/AitorRuizRedondo)
    
  - [ ] Palaeolithic portable art: ...
  
  - [ ] Azilian painted peebles: ...
  
  - [ ] Near East PPNB painting: ...
  
  - [ ] Early Neolithic anthromorph potteries:
    - [ ] [Johanna Recchia](https://www.theses.fr/236657178)
    
  - [ ] Sahara rock-art
    - [ ] [Fatima Zohra Khaled](http://www.theses.fr/2015MON30070)
    
  - [ ] Levantine/Macro-Schematic rock-art: 
    - [ ] Ines Domingo Sanz
    - [ ] Esther Lopez-Montalvo
    
  - [ ] Chalcolithic stelae (Rouergue, Provence, Languedoc): 
    - [ ] [Jules Masson Mourey](http://www.theses.fr/s163490)
    
  - [ ] Schematic rock-art: 
    - France: 
      - [ ] [Claudia Desfrasne](https://lampea.cnrs.fr/spip.php?article3640)
    - Spain: ...
    
  - [ ] Cups-and-rings: 
    - Bretagne: 
      - [ ] Serge Cassen
    - Galicia (`PENA TU`, etc.): ...
    - Great Britain: 
      - [ ] Guillaume Robin
      - [ ] Aron Mazel
      - [ ] Marta Diaz-Guardamino 
    
  - [ ] Scandinavian rock-art:
    - [ ] [Johan Ling](https://www.ancientportsantiques.com/wp-content/uploads/Documents/PLACES/UK-EUNorth/CopperScandinavia-Ling2015.pdf)
    - [ ] [Nimura, Courtney](https://www.researchgate.net/publication/305391767_Rock_art_and_coastal_change_in_Bronze_Age_Scandinavia/figures?lo=1)
    - [ ] [Gjerde J.](https://munin.uit.no/bitstream/handle/10037/2741/401-506gjerde-thesis-5.pdf?sequence=8&isAllowed=y)
  
  - [ ] Mycenian potteries with figurative decorations: ...
  
  - [ ] Warrior Stelae: 
    - [ ] Pierre-Yves Milcent
    - [ ] Marta Diaz-Guardamino
  
  - [ ] Mailhac potteries with figurative decorations:
    - [ ] Gomes de Soto
  
  - [ ] First Iron Age potteries with figurative decorations (Sopron, Darslup, etc.)
    - [ ] Christian Huth
  
  - [ ] Second Iron Age "bas aragon" stelae (`STELE BAS ARAGON`, etc.): ///
    
## By sites

  - [ ] Mount Bego: 
    - [ ] Nicoletta Bianchi
  - [ ] Valcamonica (`VALCAMONICA`):
    - [ ] Alberto Marretta
    - [ ] Andrea Arca
    - [ ] Paolo Rondini
  - [ ] Morro du Chapéu: ...
    
## By geographical areas
 
  - [ ] Australian arborigen rock-art:
    - [ ] Ines Domingo Sanz
  - [ ] South African San rock-art: 
    - [ ] Aron Mazel
  - [ ] Australian churingas: ...
  - [ ] Tihuanacu potteries with figurative decorations: ...
  - [ ] Aztec codex: ...
  - [ ] Native Americans sustaining memory paintings (e.g. Dakota Bible): ...
  
## By themes

 - [ ] "Mappe di pietra":
   - [ ] Andrea Arca
   
---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
---

Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://stackoverflow.com/> or <https://community.rstudio.com/>.

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](http://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

Brief description of the problem

```r
# insert reprex here
```
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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
# Contributing to iconr

This outlines how to propose a change to iconr.

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation (`.Rd` files) directly using the GitHub web interface.  

## Bigger changes

To develop the versatility of the package, we encourage contributions directly related to:

* [Shape analysis](https://zoometh.github.io/iconr/articles/next.html#shape-analysis-1): in incor 0.1.0, graphical units are considered for their centroid, a `POINT`. For the next releases of the package, we would like the possibility to deal with `POLYGON` geometries, in order to process the graphical units with shape analysis (Procrustes analysis, etc.).  

* [Typology](https://zoometh.github.io/iconr/articles/next.html#typology-1): we aim to create a hierarchical vocabulary (ie, structured vocabulary), in the form of a directed graph, to describe both graphical units typology and technology, mostly for Prehistoric and Protohistoric iconography (like [this](https://raw.githubusercontent.com/zoometh/iconr/master/doc/img/typology_gu.png)). To get this thesaurus shared between different users, we need it in in different languages, easily editable, etc.  

* [Harris Matrix](https://zoometh.github.io/iconr/articles/next.html#harris-matrix-1): incor 0.1.0 has a [*diachronic* edge](https://zoometh.github.io/iconr/articles/index.html#contemporaneous-elements) which is the first step to model overlaping/superimposition/relative chronology of the iconographical content. We want to go further, for example to recreate dynamically a Harrix Matrix for each graphical item.

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it's needed. If you've found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("zoometh/iconr", fork = TRUE)`.

*   Make sure the package passes R CMD check by running `devtools::check()`. If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of [`NEWS.md`](https://github.com/zoometh/iconr/blob/master/NEWS.md).


## Code of Conduct

Please note that the iconr project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
---
title: 'Analysis of Prehistoric Iconography with the R package iconr'
tags:
  - Iconography
  - Semiotic
  - Prehistory
  - Archaeology
  - Graph Theory
  - GIS
authors:
  - name: Thomas Huet
    orcid: 0000-0002-1112-6122
    affiliation: 1
  - name: Jose M Pozo
    orcid: 0000-0002-0759-3510
    affiliation: 2
  - name: Craig Alexander
    orcid: 0000-0001-7539-6415
    affiliation: 2
affiliations:
  - name: LabEx ARCHIMEDE, ANR-11-LABX-0032-01
    index: 1
  - name: Independent Researcher
    index: 2
date: 01 April 2021
bibliography: paper.bib
---

# Background

By definition, prehistorical societies are characterized by the absence of a writing system. During, the largest part of human history, and everywhere in the world, symbolic expressions belong mostly to illiterate societies which express themselves with rock-art paintings, pottery decorations, figurines, statuary, etc., and a lot of now disappeared carved woods, textile design, etc. These graphical expressions are the most significant remaining part of humankind's symbolism. At the composition level, the presence of recurrent patterns of signs (i.e., graphical syntax) in meaningful associations indicates the existence of social conventions in the way to display and to read these expressions. Well-established and shared methods to record and study these graphical contents would open the possibility of cross-cultural comparisons at a large scale and over the long-term.

# Statement of need

Ancient iconography is often perceived as different from other '*current*' archaeological remains [lithics, potteries, settlements, etc., @Chenorkian95]. Indeed, the inherent variability of ancient iconography has led to considerable problems in its study, drastically limiting the possibility to draw a synthesis of graphic expressions at a large scale and over the long-term:

 + Spatial proximities between the graphic units are not precisely quantified. Graphical units are attached to sub-areas of the support (e.g. upper part of a rock, neck of a pottery, centre of a stele).
 + Groupings -- like graphical units grouped into *figures*, *figures* grouped into *patterns*, *patterns* grouped into *motives*, etc. -- are not self-explanatory and introduce a tedious number of groups and hinder their systematic analysis.
 + Relationships and similarities between these groups are often not self-explanatory and unquantified.
 + Descriptive vocabularies and methods of analysis are site-dependent or period-dependent.

Even the reevaluation of semiotics paradigms following the scientific trends  -- *structuralist turn* during the *Processual archaeology* period, ca 1960-1980 [@Saussure89; @Binford62], *iconic turn* during the *Post-processual archaeology* period, ca 1980-2010 [@Gell98; @Hodder82], did not led to the development of efficient tools for ancient iconography studies, such as common descriptive variables, or common interpretation grids.

# Core functionality

The R package `iconr` is designed to offer a greater normalization of quantitative indexes for iconography studies [@Alexander08; @Huet18a]. It is grounded in graph theory and spatial analysis to offer concepts and functions for modeling prehistoric iconographic compositions and preparing them for further analysis: clustering, typology tree, Harris diagram [i.e. temporal succession of archaeological contexts, @Harris14], etc. The main principle of the `iconr` package is to consider any iconographic composition (here, 'decoration') as a geometric graph of graphical units. Geometric graphs, also known as *planar graphs* or *spatialized graphs*, allow to model the neighborhood of these graphical unit which are the fundamental relationships of visual semiotics [@SaintMartin11]. Graphical units are decorated surfaces (`POLYGONS`) modeled as nodes (`POINTS`) and tagged with semantic content (type, color, orientation, etc.). Separable graphical units showing a main graphical content (e.g., type = anthropomorphic figure) are considered as *main* nodes. Graphical units showing a specification of a *main* node (e.g. a sword handed by this anthropomorphic figure) are considered as *attribute* nodes. Each pair of *main* nodes thought to be contemporary that share a border (binary topological relationship: *touches*) of their Voronoi cells, are connected by an undirected edge (`LINES`).
  
  
<center>

![GIS view. The Late Bronze Age stele from Solana de Cabañas (Extremadura, Spain). 1. Original photograph (credits: Museo Arqueológico Nacional, Madrid); 2. Archaeological drawing of engraved parts [credits: @DiazGuardamino10]; 3. Digitalization/Polygonization of engraved parts (i.e., graphical units) and calculation of their their centroids (red points); 4. Voronoi diagram of each graphical unit (*seed*) and dual graph of the Voronoi diagram (i.e., Delaunay triangulation); 5. Identification of graphical units' types](https://raw.githubusercontent.com/zoometh/iconr/master/doc/img/solana_voronoi.png)

</center> 

# Overview

The `iconr` package takes in charge of the geometric graphs management (step 5 in the previous figure). Steps 1 to 4 do not need to be included in the package since efficient implementations already exist: graph elements can be drawn directly on the decorated support drawing or photograph, preferably inside a GIS to make easier the calculation of nodes and edges coordinates. The `iconr` package allows the user to i) read data structures of nodes and edges (.tsv, .csv, .shp) and images (.jpg, .png, .tif, .gif, etc.), ii) plot nodes and edges separately, or together (geometric graph), over the decoration picture, iii) compare different decorations depending on common nodes or common edges. The package stable version is on the CRAN [@iconr]; the latest development version is available from GitHub (https://github.com/zoometh/iconr); the package documentation is available at https://zoometh.github.io/iconr/.

# Examples

## Read

Read the nodes of the Cerro Muriano 1 stele (Andalusia, Spain) with the function [`read_nds()`](https://zoometh.github.io/iconr/reference/read_nds.html).

```r
library(iconr)
dataDir <- system.file("extdata", package = "iconr")
site <- "Cerro Muriano"
decor <- "Cerro Muriano 1"
read_nds(site, decor, dataDir)
```
```
##            site           decor id          type        x         y
## 1 Cerro Muriano Cerro Muriano 1  1    personnage 349.8148 -298.3244
## 2 Cerro Muriano Cerro Muriano 1  2        casque 349.8148 -243.9851
## 3 Cerro Muriano Cerro Muriano 1  3         lance 238.4637 -298.3244
## 4 Cerro Muriano Cerro Muriano 1  4      bouclier 446.0222 -381.1697
## 5 Cerro Muriano Cerro Muriano 1  5        peigne 283.0041 -358.0086
## 6 Cerro Muriano Cerro Muriano 1  7 sexe_masculin 342.6884 -427.4917
## 7 Cerro Muriano Cerro Muriano 1  8    lingot_pdb 451.1489 -237.4782
```

## Plot

Plot the Cerro Muriano 1 stele decoration graph with the function [`plot_dec_grph()`](https://zoometh.github.io/iconr/reference/plot_dec_grph.html).

```r
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
imgs <- read.table(paste0(dataDir, "/imgs.tsv"),
                   sep="\t", stringsAsFactors = FALSE)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir)
```

<center>

![R view. Cerro Muriano 1 decoration graph. Between two *main* nodes, *normal* edges are shown as plain lines. Between *main* nodes and *attribute* nodes, *attribute* edges are shown as dotted lines drawing [credits: @DiazGuardamino10]](https://raw.githubusercontent.com/zoometh/iconr/master/doc/img/cm1.png){width=350px}

</center> 

## Compare

Compare and classify the `iconr` decoration training dataset according to pairwise comparisons between decorations based on their common nodes and common edges; functions [`list_dec()`](https://zoometh.github.io/iconr/reference/list_dec.html) and [`same_elements()`](https://zoometh.github.io/iconr/reference/same_elements.html).

```r
imgs <- read.table(file.path(dataDir, "imgs.csv"), sep=";")
nodes <- read.table(file.path(dataDir, "nodes.csv"), sep=";")
edges <- read.table(file.path(dataDir, "edges.csv"), sep=";")
lgrph <- list_dec(imgs, nodes, edges)
df.same_edges <- same_elements(lgrph, "type", "edges")
df.same_nodes<- same_elements(lgrph, "type", "nodes")
dist.nodes <- dist(df.same_nodes, method = "euclidean")
dist.edges <- dist(df.same_edges, method = "euclidean")
hc.nds <- hclust(dist.nodes, method = "ward.D")
hc.eds <- hclust(dist.edges, method = "ward.D") 
par(mfrow=c(1, 2))
plot(hc.nds, main = "Common nodes", cex = .8)
plot(hc.eds, main = "Common edges", cex = .8)
```
<center>

![Results of the hierarchical clustering on the `iconr` decoration training dataset (five Late Bronze Age stelae) on common nodes (left) and common edges (right)](https://raw.githubusercontent.com/zoometh/iconr/master/doc/img/hc.png)

</center> 

# Acknowledgements

This  project was partly  supported  by  the LabEx  ARCHIMEDE  from  “Investissement  d’Avenir”  program  ANR-11-LABX-0032-01.

# References
---
title: "Interactive examples on training datasets"
---

```{r, echo=FALSE}
# print(getwd())
# print(file.path("../logo", "iconr_logo.png"))
htmltools::img(src = knitr::image_uri(file.path("../man/figures", "iconr_logo.png")),
               # htmltools::img(src = knitr::image_uri(file.path("img", "iconr_logo.png")), 
               # htmltools::img(src = knitr::image_uri(file.path(R.home("doc"), "html", "logo.jpg")), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;')
# setwd("C:/Users/supernova/Dropbox/My PC (supernova-pc)/Documents/iconr/doc/datasets")
```

<style>
.shiny-app-frame {
  position: fixed;
  left: 0;
  top: 250px;
  bottom: 200px;
  right: 0;
}
.shiny-app-frame iframe {
  width: 100%;
  height: 100%;
  border: none;
}
</style>

<div>
  <p>This RShiny app show examples on how to record the GUs and their proximity links with the <b>iconr</b> package 
<a href="https://zoometh.github.io/iconr/reference/plot_dec_grph.html">plot_dec()</a> function</p>
</div>
  
<div class="shiny-app-frame"> 
<iframe src="https://epispat.shinyapps.io/decoration_graph">
</iframe>
</div>
---
title: "GIS interface entry"
author: "Thomas Huet, Jose M Pozo, Craig Alexander"
bibliography: references.bib
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
---

```{r, echo=FALSE}
# print(getwd())
# print(file.path("../logo", "iconr_logo.png"))
htmltools::img(src = knitr::image_uri(file.path("../man/figures", "iconr_logo.png")),
# htmltools::img(src = knitr::image_uri(file.path("img", "iconr_logo.png")), 
# htmltools::img(src = knitr::image_uri(file.path(R.home("doc"), "html", "logo.jpg")), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;')
```

---

<style>
.figure {
   margin-top: 0px;
   margin-bottom: 40px;
}
table {
    margin-top: 0px;
    margin-bottom: 24px;
}
</style>


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(kableExtra)
library(dplyr)

solanas.vor.path <- "../doc/tutorial/img/all_process.gif"
```
  

Any iconographic contents can be modeled with a geometric graph where nodes, also called **graphical units (GUs)**, linked together with edges and then analyzed with the **Graph Theory** and spatial analysis at the support scale.

# Forewords

The ***iconr*** R package grounds concepts and provides normalized tools to manage iconographic contents. This modeling is particularly interesting for compositions coming from Paleolithic, Neolithic and Protohistoric times [@Alexander08; @HuetAlexander15; @Huet18a]. To record large series of iconographic contents, the [GIS interface](#gis) appears to be the most appropriate one for users. This demo explains how to construct the graph _before_ using the package, offering also tipping points to facilitate the recording process.

<center>

![GIS view. A Late Bronze Age stelae (Solana de Cabañas, Cáceres, Spain).Stelae photograph (photograph credits: Museo Arqueológico Nacional, Madrid); Georeferencing of the steale drawing over its photograph (drawing credits: Diaz-Guardamino 2010); Binarization and polygonization/vectorization of the graphical content (graphical units, GUs) of the stelae (now `POLYGONS`); Calcul of their centroid (red points); Calcul of their Voronoi cells; Binary topological relationships (*birel*) for each pairwise of Voronoi cells: the ones that share a border (*touches* = TRUE) will share a link (ie, edge, in red) between their centroids (ie, nodes); Identification of the different types (nodes' column `type`)](`r solanas.vor.path`)

</center>


# A GIS ? Yes but why ? {#gis}

The *tenet* of the ***iconr*** is to always keep the user connected with the iconographic content -- his primary data source -- and emphasize the significance of the spatial dimension for any graphical content. Geographical Information Systems (GIS) offer multiple tools and options to facilitate the data entry.  Use of GIS offers a graphic interface and ensures the correctness of spatial relationships between GUs. It forms a permanent interface between the image of he decoration and the database. Obviously, the main GIS facility is the presence of **attributes tables** which allow to record, filter and sort GUs on many information: types, techniques, orientations, lengths, etc. 

The other most important GIS facilities for the recording process are:

* edition tools
  + snapping tools  
  
* scale tools (see [Absolute scale](#scale.abs))
  + measurement tools
  + georeferencing tools

From far, our software preference goes to [QGIS](https://www.qgis.org/fr/site/), because it is open source, offers a large rank of database connections facilities (with PostGRES/GIS for example), has a large user community, but also because the source file is a XML (.qgs, .qgz) structure that can be parsed, modified, copied and moved with scripting languages like R and Python

# Always start with an **image**

The image will be the reference space of the graph. So, before anything, start by opening the image decoration into the GIS. In this tutorial we will take the example of the South Western Iberian **Abela stela** dated to the Middle Bronze Age. The original drawing can be download [here](https://github.com/zoometh/iconr/blob/master/docs/tutorial/img/Abela.jpg) [@DiazGuardamino10]


<center>
  
![Drag and drop the decoration image before anything](../doc/tutorial/img/openImageGIS.gif){width=70%}
</center>
  

## Scales

Use of GIS easy the scaling process. The creation of a spatialized graph permit to combine network analysis with spatial analysis of the graphical content

### Relative scale {#scale.rel}

The image extent is measured in pixels with a top-left corner origin (0,0). The coordinates system is irrelevant: image, nodes and edges are measured on the pixel grid


<center>

![Coordinates system in pixels](../doc/tutorial/img/coordsQGIS.gif){width=70%}
</center>

### Absolute scale {#scale.abs}

To retrieve to true scale of the decoration, one can create a scale bar and apply a simple rule of three to convert pixels into centimeters, or meters. For example, if the scale belongs to another drawing, you can import it and 'georeferenced' it on the original drawing with the [*Freehand raster georeferencer* plugin](https://gvellut.github.io/FreehandRasterGeoreferencer/), and then create the scale bar

<center>

![Importing information on the scale and creation of a scale bar](../doc/tutorial/img/scaleAll.gif){width=70%}

</center>

Parallely, the dimensions of each GU can be measured with the QGIS **Measure Line** tool. At first, only the maximum length of the GU is important. It has also to be noted that if a Polygonization is done on the GUs, the maximum length -- between all other type of shape analysis indexes -- do not have to be calculated 

<center>

![Measure and record the dimensions of each GUs in pixels](../doc/tutorial/img/scaleRuleOfThree.gif){width=70%}

</center>

To retrieve the real dimensions of each GUs, first get the pixel sizes with the Measure line tools, then apply a simple rule of three with:

1. the size of the scale in pixels (native QGIS function `$length`)
2. the real size of the scale in cm (here, 100 cm)
3. the size of each GUs in cm


# **Nodes** and **edges** tables

[Nodes](#nodes) and [edges](#edges) are the graph elements. At first, we have to create attribute tables for each of them. For example the nodes shapefiles:


<center>
  
![Create the **nodes** shapefile with its attribute table](../doc/tutorial/img/createNodes.gif){width=70%}
  
</center>

In the GIS, [add a node](#graph.enodes.add) for each GUs and [add an edge](#graph.edges.add) between two contiguous GUs. Always start with the nodes

## Nodes {#nodes}

**Nodes** represent the basic information of the graphical content. For example, it would be easy to distinguish a decoration with aurochs (`type = auroch`) from a decoration with swords (`type = sword`). The former iconographic content should be probably related to the Late Paleolithic/Early Neolithic while the second one is more probably related to a period starting with the Bronze Age. Nodes are created as a shapefile of `POINTS`. The attribute table of the nodes has at least four (4) fields:

1. `site` (Text): name of the site
2. `decor` (Text): name of the decoration
3. `id` (Integer): node identifier
4. `type`(Text): one of the relevant characteristics of each node

The nodes are created near the centroids of each different graphical units (GUs).

<center>
  
![Adding **nodes** for each GU](../doc/tutorial/img/addNodes.gif){width=70%}
  
</center>

In this example, beside the *main* nodes sword (`epee`), anciform (`anciforme`) and halberd (`hallebarde`), we can also observe that the sword is connected to a belt and the anciform is worn as necklace. We probably would like to register this piece of graphical content as we also probably would like to characterize the types of blades for the sword and halberd, point out the presence of rivets on the sword depiction, etc. To do so, create 'Attribute edges' (`-+-`)

## Edges {#edges}

**Edges** types (field `type`) give information on nodes relative locations and on the nature of these nodes (main node *vs* attribute node, overlapping *vs* overlapped node, etc.). Edges are created as a shapefile of `LINES`. Edges attribute table has at least five (5) fields:

1. `site` (Text): name of the site
2. `decor` (Text): name of the decoration
3. `a` (Integer): *starting* node
4. `b` (Integer) *ending* node
5. `type` (Text): values `=`, `+` or `>`; textual notation `-=-`, `-+-` or `->-`

Theoretically, between two *main* nodes, edges exist when their Voronoi cells are contiguous. In practice, if you consider that two GUs are neighbors, you can create an edge between their two nodes: GIS snapping tool help !

<center>
  
![Adding edges between contiguous GU](../doc/tutorial/img/addEdges.gif){width=70%}
  
</center>


## Summary

For the Abela decoration, we have created three (3) nodes (`1`,`2`,`3`) and two (2) *normal* edges (`1-=-2`,`1-=-2`). We named the nodes shapefile `nodes.shp` and the edges shapefile `edges.shp` because this is their default name in the ***incor*** package

<center>
  
![](../doc/tutorial/img/inputSummary.gif){width=70%}
  
</center>

# The **table of decorations**

The **table of decorations** is a correspondence table which records joins between nodes and edges dataframes. In the package, the default name of this table is `imgs.tsv`, a tabulate separated-values (but it also can be a`.csv`, comma separated-values)

<center>
  
![](../doc/tutorial/img/createTableOfDecorations.gif){width=70%}
  
</center>

The **table of decorations** has four (4) mandatory fields:

1. `idf` (Integer): short name of the decoration
1. `site` (Text): name of the site
2. `decor` (Text): name of the decoration
4. `img` (Text): name of the drawing/photograph/...


# References

---
title: "Training datasets"
---

```{r, echo=FALSE}
# print(getwd())
# print(file.path("../logo", "iconr_logo.png"))
htmltools::img(src = knitr::image_uri(file.path("../man/figures", "iconr_logo.png")),
               # htmltools::img(src = knitr::image_uri(file.path("img", "iconr_logo.png")), 
               # htmltools::img(src = knitr::image_uri(file.path(R.home("doc"), "html", "logo.jpg")), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;')
# setwd("C:/Users/supernova/Dropbox/My PC (supernova-pc)/Documents/iconr/doc/datasets")
```

---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(shiny)
library(iconr)
library(dplyr)
library(knitr)
library(kableExtra)

downgit.root <- "https://downgit.github.io/#/home?url="
gh.root <- paste0(downgit.root, "https://github.com/zoometh/iconr/")
dir.root <- paste0(gh.root, "tree/master/doc/datasets/")
url.df <- c(paste0(dir.root, "stele%20bas%20aragon"),
            paste0(dir.root, "Valcamonica"),
            paste0(dir.root, "Pena%20Tu"))
url.ref <- c("https://doi.org/10.3989/gladius.2013.0002",
             "https://scholar.google.fr/scholar?hl=fr&as_sdt=0%2C5&q=Le+Pietre+degli+Dei%2C+menhir+e+stele+dell%27Eta+del+Rame+in+Valcamonica+e+Valtellina&btnG=",
             "https://doi.org/10.1016/j.anthro.2005.09.009")
df.dset <- data.frame("download" = c("stele bas aragon",
                                     "Valcamonica",
                                     "Pena Tu"),
                      "name" = c("Estelas Ibericas con Lanzas",
                                 "Ossimo",
                                 "Peña Tu style"),
                      "nb.decor" = c("5",
                                     "3",
                                     "5"),
                      "description" = c("Iron Age stelae with ranges of spears, writings, etc.",
                                        "Chalcolithic rock-art from Ossimo, Valcamonica",
                                        "Bell-Beaker stelae and rock-art"),
                      "drawings.references" = c("Vargas 2013",
                                                "Fedele 1994",
                                                "Ramírez 2005"),
                      stringsAsFactors = F)
```


These training datasets can be downloaded from the below dataframe. We also provide ready-to-use QGIS projects (.qgz) for each decoration:

```{r, echo = F}
# row.names(df.dset) <- df.dset$family
df.dset %>% 
  mutate(download = cell_spec(download, "html", link = url.df)) %>%
  mutate(drawings.references = cell_spec(drawings.references, "html", link = url.ref)) %>%
  kable("html", escape = FALSE, row.names = F) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F)
```
---
title: "Analysis of Prehistoric Iconography with R"
author: "Thomas Huet, Jose M Pozo, Craig Alexander"
email: "thomashuet7@gmail.com"
bibliography: references.bib
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Analysis of Prehistoric Iconography with R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style>
.figure {
   margin-top: 0px;
   margin-bottom: 20px;
}
table {
    margin-top: 0px;
    margin-bottom: 24px;
}
</style>

```{r, include = FALSE}
library(knitr)
library(igraph)
library(dplyr)
library(kableExtra)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">",
  fig.pos = 'H'
)
# ibahernando.path <- paste0(getwd(),"/img/ibahernando_256colours.png")
# brozas.path <- paste0(getwd(),"/img/brozas_256colours.png")
# solanas.path <- paste0(getwd(),"/img/solana_detail_256colours.png")
# solanas.vor.path <- paste0(getwd(),"/img/solana_voronoi_256colours.png")
# ibahernando.path <- "img/ibahernando_256colours.png"
# brozas.path <- "img/brozas_256colours.png"
# solanas.path <- "img/solana_detail_256colours.png"
# solanas.vor.path <- "img/solana_voronoi_256colours.png"
ibahernando.path <- "../man/figures/ibahernando_256colours.png"
brozas.path <- "../man/figures/brozas_256colours.png"
solanas.path <- "../man/figures/solana_detail_256colours.png"
solanas.vor.path <- "../man/figures/solana_voronoi_256colours.png"
```
  
The R package `iconr` is grounded in graph theory and spatial analysis. It offers concepts and functions for modeling Prehistoric iconographic compositions and for their preparation for further analysis (clustering, Harris diagram, etc.) in order to contribute to cross-cultural iconography comparison studies through a greater normalization of quantitative analysis [@Alexander08; @HuetAlexander15; @Huet18a].  
  
The flexibility of graph theory and tools available for the GIS database make the `iconr` package useful in managing, plotting and comparing (potentially large) sets of iconographic content: Atlantic rock-art, Scandinavian rock art, Late Bronze Age stelae, Mycenean figurative pottery, etc. 
  
# Decoration graphs  

The main principle of the `iconr` package is to consider any iconographic composition (here, 'decoration') as a geometric graph of graphical units (GUs). This geometric graph is also known as a planar graph or spatialized graph. The GUs are decorated surfaces (`POLYGONS`) modeled as nodes (`POINTS`). When these GUs are *main* nodes, and not *attribute* nodes, they share edges (`LINES`) with one another when their Voronoi cells share a border (*birel*: touches).  

<center>

![GIS view. The Solana de Cabanas stelae: from its photograph to the modeling of its graphical content](`r solanas.vor.path`){width=700px}

</center>  
  
  
Graph theory offers a conceptual framework and indices (global at the entire graph scale, local at the vertex scale) to deal with notions of networks, relationships and neighbourhoods. The geometric graph is commonly built within a GIS interface. Indeed, use of GIS allows one to create a spatial database of the decoration's iconographic contents and facilitates data recording and visualization. For example, snapping options can connect GUs (nodes) with lines (edges) and we can exploit tools such as feature symbology, layer transparency, etc.  
  
The latest development version of the `iconr` package and its vignette can be downloaded from GitHub

```{r down,eval=FALSE, echo=TRUE}
devtools::install_github("zoometh/iconr", build_vignettes=TRUE)
```

The R package `iconr` is composed of [functions](#functions) and a [example dataset](#data). The main R packages used by the `iconr` package are:

* [magick-image](https://CRAN.R-project.org/package=magick): for drawing/image management 
* [igraph](https://CRAN.R-project.org/package=igraph): for graph management 
* [rgdal](https://CRAN.R-project.org/package=rgdal): for shapefile management

Load the package `iconr`

```{r load, echo=TRUE}
library(iconr)
``` 

# Dataset {#data}

The input dataset is expected to include decoration images and corresponding node and edge data in a single data folder. This folder should include the following files:

* **[Table of decorations](#decorations):** A tabular file storing the set of decoration identifiers and corresponding image filenames.
* **[Images](#drawings):** An image file for each decoration.
* **[Node data](#nd):** A single file storing the data of each node for all decorations. 
* **[Edge data](#ed):** A single file storing the data of each edge for all decorations.

The `iconr` package includes an example dataset with the input files in several alternative formats. The path for the example dataset is the package *extdata* folder. This folder is also -- by default -- the output folder. In order to differentiate input data from output data, the output data filenames always start with a digit or an underscore (ie, a punctuation). So, the input data are :

```{r ls_ext_data}
dataDir <- system.file("extdata", package = "iconr")
input.files <- list.files(dataDir)
cat(input.files, sep="\n")
```

The table of decorations is given in two formats: comma- or semicolon-separated values (`imgs.csv`) and tab-separated value (`imgs.tsv`).

```{r paths.imgs, echo=TRUE}
imgs_path <- paste0(dataDir, "/imgs.csv")
imgs <- read.table(imgs_path, sep=";", stringsAsFactors = FALSE)
```

Each decoration is identified by its name (column `decor`) and the name of the site (column `site`) to which it belongs. In the example dataset, this is transparent in the name of each decoration image, included in jpg format. Any other image format supported by the R package `magick` (jpg, png, tiff, pdf, etc.) is suitable.

As we have stated, a GIS interface is often the most practical way to record graph nodes and graph edges with `POINTS` and `LINES` geometries, respectively. This is typically saved in shapefile (shp) format, which is composed of at least 3 files with extensions `.shp` (geometries), `.shx` (indices), and `dbf` (attribute data). The example dataset includes them for nodes and edges separately, with obvious names:

```{r paths, echo=TRUE}
nodes_path <- paste0(dataDir, "/nodes.shp")
nodes.shp <- rgdal::readOGR(dsn = nodes_path, verbose = FALSE)
nodes <- as.data.frame(nodes.shp)
edges_path <- paste0(dataDir, "/edges.shp")
edges.shp <- rgdal::readOGR(dsn = edges_path, verbose = FALSE)
edges <- as.data.frame(edges.shp)
```

Nodes and edges can also be recorded in tabular format: `.csv` or `.tsv`.  

```{r paths.1, echo=TRUE}
nodes_path <- paste0(dataDir, "/nodes.tsv")
nodes <- read.table(nodes_path, sep="\t", stringsAsFactors = FALSE)
edges_path <- paste0(dataDir, "/edges.tsv")
edges <- read.table(edges_path, sep="\t", stringsAsFactors = FALSE)
```

The list of graph decorations is created with the `list_dec()` function. These graphs are `igraph` objects

```{r graph.clss}
lgrph <- list_dec(imgs, nodes, edges)
g <- lgrph[[1]]
as.character(class(g))
```

By default, the `plot.igraph()` function (ie, `igraph::plot()`) spatialization (`layout`) is based on `x` and `y` columns, when these exist. This is appropriate in our case since we are working with geometric graphs. If this spatialization (node coordinates) is ignored, then very different layouts are possible for the same graph (*graph drawing*):

```{r igraph.1, warning=FALSE, fig.align="center", fig.width=6.5, fig.asp=0.58}
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))           
par(mar=c(1, 0, 2, 0), mfrow=c(1, 2), cex.main = 0.9, font.main = 1)
coords <- layout.fruchterman.reingold(lgrph[[1]])
plot(g,
     vertex.size = 15,
     vertex.frame.color="white",
     vertex.label.family = "sans",
     vertex.label.cex = .8,
     main = "Graph drawing based on x, y coordinates"
)
plot(g,
     layout = layout.fruchterman.reingold(g),
     vertex.size = 5 + degree(g)*10,
     vertex.frame.color="white",
     vertex.label.family = "sans",
     vertex.label.cex = .8,
     main = "Force-directed graph drawing,\nwith degree-dependent node size."
)
mtext(g$decor, cex = 1, side = 1, line = -1, outer = TRUE)
```
  
## Table of decorations {#decorations}

In any compatible file format readable by R, such as comma-separated values (csv) or tabular-separated values (tsv). Decorations identifiers and decorations image filenames are stored in a dataframe, by default the `imgs` dataframe:

```{r imgs,fig.width=6, fig.height=6, fig.align="center",warning=FALSE, fig.cap="\\label{fig:figs}imgs.tsv"}
imgs_path <- paste0(dataDir, "/imgs.tsv")
imgs <- read.table(imgs_path, sep="\t", stringsAsFactors = FALSE)
knitr::kable(imgs, "html") %>% 
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12)
```
  
The decoration's unique identifiers are the concatenation of the site name and the decoration name. For example, the name of the Cerrano Muriano 1 decoration is: ``r imgs[1,"img"]``. 

## Images {#drawings}

Image or drawing (eg, bitmap images, rasters, grids) formats accepted are the common ones (jpg, png, jpeg, tiff, pdf, etc.). The images in the current example dataset come from a PhD thesis, published by M. Diaz-Guardamino [@DiazGuardamino10]. 

For a given decoration, its image is the reference space of the graph: nodes and edges inherit their coordinates from this image. However, different image or graphic systems use different coordinate conventions. The following table illustrates some of the most relevant ones:

```{r xy_coords,out.width="50%", fig.align="center",echo=FALSE,warning=FALSE}
df.equi <- data.frame(
  "Device/Package" = c("*R graphics*", "*R raster*", "*R magick*", "***GIS interface***"),
  "Unit of measure" = c("number of pixels", "number of pixels", "number of pixels", "**number of pixels**"),
  "Origin" = c("bottom-left corner", "top-left corner", "top-left corner", "**top-left corner**"),
  "x-axis orientation" = c("rightward", "downward", "rightward", "**rightward**"),
  "y-axis orientation" = c("upward", "rightward", "downward", "**upward**"),
  check.names = FALSE)
knitr::kable(df.equi) %>%
  kable_styling(full_width = F)
```

The package `iconr` follows the coordinate conventions of the GIS interface. This is in agreement with our advice above to use the GIS as the preferred interface to extract nodes and edges from the decoration image. Observe that the y-axis is oriented upwards from the top-left corner, which implies that y-coordinates are always negative in the image (ie, downwards). The following image illustrates a decoration with the `iconr` (GIS) coordinates for the four corners. Note the contrast with the coordinates in the code generating it with the standard graphics package.

```{r drawing, out.width="50%", fig.width=6, fig.asp=750/666, fig.align="center", warning=FALSE, echo=TRUE, message=FALSE, fig.cap="\\label{fig:figs} `iconr` (GIS) coordinate convention: decoration `Cerro_Muriano.Cerro_Muriano_1.jpg` with the coordinates of its corners."}
library(magick)
library(graphics)
dataDir <- system.file("extdata", package = "iconr")
imgs_path <- paste0(dataDir, "/imgs.csv")
imgs <- read.table(imgs_path, sep=";", stringsAsFactors = FALSE)
cm1 <- image_read(paste0(dataDir, "/", imgs$img[1]))
W <- image_info(cm1)$width
H <- image_info(cm1)$height
oldpar <- par(no.readonly = TRUE)   
on.exit(par(oldpar))            
par(mar = c(0, 0, 0, 0))
plot(cm1)
box(lwd = 2)
text(0, H, paste0(0, ",", 0), cex = 2, adj = c(0, 1.1))
text(W, H, paste0(W, ",", 0), cex = 2, adj = c(1, 1.1))
text(0, 0, paste0(0, ",", -H), cex = 2, adj = c(0, -0.2))
text(W, 0, paste0(W, ",", -H), cex = 2, adj = c(1, -0.2))
```

## Node data {#nd}

Nodes can be stored in any of the following formats: comma-separated values (`.csv`), tab-separated values (`.tsv`), or shapefile (`.shp`). Any of these formats is consistenly read by the function `iconr::read_nds()` generating a data.frame with the node data. For the generic csv and tsv formats, at least the following five columns are required to be included: 

* **site**: decoration site 
* **decor**: decoration name 
* **id**: id uniquelly identifying each node
* **x**: x-coordinate of nodes
* **y**: y-coordinate of nodes 

Additional columns can be included providing relevant characteristics of each node. The example dataset includes the column:

* **type**: type of GU (ie, what it represents)

These node characteristics can then be used for image annotations, decoration comparisons and decoration analysis. Nodes attributes coming from the GIS interface:

```{r nodes.df, warning=FALSE,fig.align="center",warning=FALSE}
nds.df <- read_nds(site = "Cerro Muriano", decor = "Cerro Muriano 1", dir = dataDir) 
knitr::kable(nds.df, "html") %>% 
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12)
```
  
In principle, the nodes are defined to be located at the exact centroids of each GU. However, for practical reasons, they can also be manually located close to these centroids.  

### Nodes types {#nd.type}

Nodes can be *main* nodes or *attribute* nodes depending on the type of edges they share with other nodes (see [edge types](#ed.type)). For example, considering the Cerro Muriano 1 stelae:

* **main nodes**: the anthropomorphic figure, the spear, the shield, the comb and the ingot can be considered as *main* nodes because they are separated and they show different graphical contents.

* **attribute nodes**: the helmet and the male sex of the antropomorphic figure can be considered as *attribute* nodes of the anthropomorphic figure because they are characteristic features of this latter.

## Edge data {#ed}

Edges can be stored in the same three formats as nodes: csv, tsv, or shp. At least the following five columns or variables are required to be included:

* **site**: decoration site 
* **decor**: decoration name
* **a**: *id* of the first [node](#nd) of each edge
* **b**: *id* of the second [node](#nd) of each edge
* **type**: [edge types](#ed.type)

These columns can be seen in the example dataset:

```{r edges.df, warning=FALSE}
edges <- read.table(edges_path, sep = "\t", stringsAsFactors = FALSE)
knitr::kable(head(edges), "html", align = "llccc") %>%
  kableExtra::kable_styling(full_width = FALSE, position = "center",
                            font_size=12) %>%
  gsub("\\+", "$+$", .)
```
  
Analogously to nodes, any of these formats is consistently read by the function `iconr::read_eds()` generating a data.frame with the edges data. Additionally, this function includes columns with the coordinates of the edge nodes:

* **xa**: x-coordinate of the first node of each edge
* **ya**: y-coordinate of the first node of each edge 
* **xb**: x-coordinate of the second node of each edge
* **yb**: y-coordinate of the second node of each edge 

Those coordinates are imported from the `x` and `y` coordinates in the corresponding [node](#nd) file, identified by the node `id` as stated in the edge fields `a` and `b`, respectively.

```{r edges.df.1, warning=FALSE}
eds.df <- read_eds(site = "Cerro Muriano", decor = "Cerro Muriano 1", dir = dataDir) 
knitr::kable(head(eds.df), "html", align = "llcccrrrr") %>% 
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12) %>%
  gsub("\\+", "$+$", .)
```
  
**column names** 

* `type` : [edges types](#ed.type)  

  * `=` : *normal* edges between contiguous and contemporaneous *main* nodes (undirected edge)  

  * `+` : *attribute* edges, between contemporaneous nodes where the *attribute* node `b` is an attribute of *main* node `a` (directed edge)

  * `>` : *diachronic* edges, between non-contemporaneous nodes where the node `a` overlaps node `b`, or node `a` is more ancient than node `b` (directed edge)

* `xa, ya`: coordinates of the *starting* node, or *main* node, or *overlapping* node, or *more recent* node  (`a`) 

* `xb, yb`: coordinates of the *ending* node, or *attribute* node, or *overlapped* node, , or *more ancient* node (`b`) 
  
  
### Edge types {#ed.type}

Graph theory says that edges can be undirected or directed. In the `iconr` package, by default: 

* all contemporaneous nodes have [*normal* edges](#ed.type.norm) or [*attribute* edges](#ed.type.attrib)  edges displayed in <span style="color:orange"><b>orange</b></span>

* all non-contemporaneous nodes have [*diachronic* edges](#ed.type.over) displayed in <span style="color:blue"><b>blue</b></span> edges (see [contemporaneous nodes](#contemp))

The `named_elements()` function allows one to display the textual notation of the different types of edges (`-=-`, `-+-` or `->-`)

```{r count4, warning=FALSE}
named_elements(lgrph[[1]], focus = "edges", nd.var="type")[1]      
```

When there are nodes with the same `nd.var`, this function adds the suffix `#` to the `nd.var` in order to disambiguate the node list. This is the case, for example, for the `chariot_char-+-cheval` (x2) and `chariot_char-+-roue` (x2) edges of the Zarza de Montanchez stelae (decoration 4)

```{r count1, warning=FALSE}
named_elements(lgrph[[4]], focus = "edges", nd.var="type")
```

Employed with the basic R functions for loop (`lapply()`), count (`table()`) and order (`order()`), and removing the the suffix `#`, this function can be used to count the different types of edge. Here, we enumerate the most represented types of the example dataset:

```{r count.all, warning=FALSE}
all.edges <- unlist(lapply(lgrph, named_elements, 
                           focus = "edges", nd.var="type", disamb.marker=""))
edges.count <- as.data.frame(table(all.edges))
edges.order <- order(edges.count$Freq, decreasing = TRUE)
edges.count <- edges.count[edges.order, ] 
knitr::kable(head(edges.count), row.names = FALSE) %>%
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12)
```
  
#### ***Normal* edges** {#ed.type.norm}

Different *main* nodes considered to be contemporaneous and close to one another, may share an edge with the  value `=` for their `type`. By convention, these edges are called *normal* and displayed as a plain line. Its textual notation is `-=-`. The *normal* edges are undirected:  `1-=-2` is equal to `2-=-1`, where node 1 and node 2 are two different *main* nodes.  

#### ***Attribute* edges** {#ed.type.attrib}

When a node is an attribute of another, edges are identified with a `+` and displayed with a dashed line. For example, at the bottom of the Zarza de Montanchez stelae (decoration 4), the main node 7 (`chariot`) is connected with four (4) attribute nodes:

- two horses (`cheval`): 8 and 9  

- two wheels (`roue`): 10 and 11 

```{r graph.attribute.plot.type, out.width="60%", fig.width=6, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Zarza De Montanchez stelae (decoration 4) showing *normal* and *attribute*  edges"}
site <- "Zarza de Montanchez" 
decor <- "Zarza De Montanchez"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              ed.lwd = 1, ed.color = c("darkorange"),
              lbl.size = 0.7)
```
  
The textual notation of an *attribute* edge is `-+-`:

```{r count_att, warning=FALSE}
sort(named_elements(lgrph[[4]], focus = "edges", nd.var = "type"))
```
  
The *attribute* edges are directed: `1-+-2` is not equal to `2-+-1`, `1-+-2` means that node 1 is the *main* node and node 2 is one of its *attribute* nodes. 
 
#### ***Diachronic* edges** {#ed.type.over}

When a node overlaps another or is more recent than another, edges are identified with a `>` and displayed with a blue plain line. For example, the Ibahernando stelae has a latin inscription (`ecriture`) overlapping a spear (`lance`) and a shield (`bouclier`).  

```{r graph.overlap.plot.type, out.width="100%", fig.width=12, fig.asp=0.55, fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae (decoration 5) showing *diachronic* and *normal* edges"}
site <- "Ibahernando"
decor <- "Ibahernando"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))          
par(mfrow = c(1, 2))
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              lbl.size = 0.7)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = 'type',
              lbl.size = 0.6)
```
  
The textual notation of an *diachronic* edge is `->-`:

```{r count3, warning=FALSE}
named_elements(lgrph[[5]], focus = "edges", nd.var = "type")
```

The *diachronic* edges are directed: `1->-2` is not equal to `2->-1`. The `1->-2` edge means that node 1 overlaps node 2, or node 1 is more recent than node 2. These overlays, or diachronic iconographic layers, can be managed with the `contemp_nds()` function (see the section [Contemporaneous contents](#contemp))  

# Functions {#functions}

Decoration graphs are constructed from nodes and edges. Graphs are 1-component: each decoration graph covers all the GUs of the decoration. The functions of the `iconr` package provide basic tools to manage these [node data](#nd) (`.csv`, `.tsv` or `.shp`) and [edge data](#ed) (`.csv`, `.tsv` or `.shp`) to create graphs, to [plot](#plot) and to [compare](#compare) them, to select [contemporaneous GU compositions](#contemp)  

```{r ls_functions}
cat(ls("package:iconr"), sep="\n")
```

## Read {#read}

The functions `read_nds()` and `read_eds()` allow, respectively, to read a dataframe or a shapefile of nodes, and a dataframe or a shapefile of edges. For example, the nodes and edges of Cerro Muriano 1 can be read

```{r img.graph.read, echo=TRUE}
site <- "Cerro Muriano"
decor <- "Cerro Muriano 1"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
```

## Plot {#plot}

The graphical functions `plot_dec_grph()` and `plot_compar()` allow, respectively, to plot one decoration or pairwise(s) of decorations. These function offer different choices for the color and size of the nodes, edges or labels. For example, for the Cerro Muriano 1 decoration, the field `id` (identifier of the node), used by default for the labels, can be changed to the `type` field

```{r img.graph.plot.type, out.width="60%", fig.width=6, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Cerro Muriano 1 stelae (decoration 1) with the type of each GU"}
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = 'type',
              lbl.size = 0.55)
```
  
A new field, `long_cm`, is added to the Cerro Muriano 1 nodes and the graph is replotted using this field instead of the `type` field, with <span style="color:brown">brown</span> colors and larger labels.

```{r img.graph.plot.id, out.width="60%", fig.width=6, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Cerro Muriano 1 stelae (decoration 1) with the maximum length (in cm) of each GU"}
nds.df$long_cm <- paste0(c(47, 9, 47, 18, 7, 3, 13), "cm")
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = 'long_cm',
              nd.color = "brown",
              lbl.color = "brown", lbl.size = 0.7,
              ed.color = "brown")
```
  
## Compare {#compare}

Elements of the graphs (nodes and edges) can be compared across all graphs or between a pair of graphs with the `same_elements()` and `plot_compar()` functions

* `same_elements()` function permits one to count common elements between *n* graphs   

* `plot_compar()` function shows a graphical output for these common elements

By default, in a pairwise comparison of decorations, common nodes and edges are displayed in <span style="color:red">red</span>, but their colors -- and other graphical parameters -- can be modified. When not all GUs are contemporaneous with one another, the non-contemporaneous ones can be removed with the `contemp_nds()` function. 

### Node comparisons

A classic study in archaeological research is to count the common nodes between pairs of decorations. This can be done with the `same_elements()` function with a node focus (`focus = "nodes"`), and considering for example their `type` (by default).

```{r compare.nodes, results='asis', warning=FALSE}
imgs_path <- paste0(dataDir, "/imgs.tsv")
nodes_path <- paste0(dataDir, "/nodes.tsv")
edges_path <- paste0(dataDir, "/edges.tsv")
imgs <- read.table(imgs_path, sep="\t", stringsAsFactors = FALSE)
nodes <- read.table(nodes_path, sep="\t", stringsAsFactors = FALSE)
edges <- read.table(edges_path, sep="\t", stringsAsFactors = FALSE)
lgrph <- list_dec(imgs, nodes, edges)
df.same_nodes <- same_elements(lgrph,
                               focus = "nodes",
                               nd.var = "type")
diag(df.same_nodes) <- cell_spec(diag(df.same_nodes),
                                 font_size = 9)
knitr::kable(df.same_nodes, row.names = TRUE, escape = FALSE, table.attr = "style='width:30%;'",
             caption = "Count of common nodes between decorations") %>%
  column_spec(1, bold=TRUE) %>%
  kableExtra::kable_styling(position = "center", font_size = 12)
```
  
The result of `same_elements()` is a symmetric matrix giving the number of common nodes for each pair of decorations, where the row and column names are the identifiers of the decorations. Observe that, accordingly, the diagonal elements show the total number of nodes of each decoration.

Regarding the node variable `type`, the decoration 4 has a total twelve (12) nodes (see the diagonal). Decoration 4 and decoration 2 have nine (9) common nodes. Decoration 4 and decoration 3 have four (4) common nodes. This matrix can be used for further [clustering analysis](#sum)

To compare graphically the decorations 2, 3 and 4 on the node variable `type`:

* first: `type` variable is pasted to the `list_compar()` function  

* then: the plot is made with the `plot_compar()` function
  
```{r compare.2.nodes, fig.show = TRUE, out.width="100%", fig.width=12, fig.asp=0.52, fig.align="center", warning=FALSE}
dec.to.compare <- c(2, 3, 4)
g.compar <- list_compar(lgrph, nd.var = "type")
plot_compar(listg = g.compar, 
            dec2comp = dec.to.compare,
            focus = "nodes",
            nd.size = c(0.5, 1.5),
            dir = dataDir) 
```
  
The function creates an image for each pair of stelae contained in the `dec.to.compare` variable (`r dec.to.compare`), with a focus on nodes (`focus = "nodes"`). Thus, if $n$ decorations are compared, it results in $n\choose 2$ $\frac{n!}{(n-2)!2!}$ pairwise comparison images. For instance, if all 5 decorations in the example dataset are compared, there will be 10 pairwise comparisons.

### Edge comparisons

A less-common study in archaeological research is to count the common edges between pairs of decorations. The `same_elements()` function with an edge focus (`focus = "edges"`) and considering the  `type` of the nodes 

```{r compare.edges, warning=FALSE}
df.same_edges <- same_elements(lgrph, nd.var = "type", focus = "edges")
diag(df.same_edges) <- cell_spec(diag(df.same_edges),
                                 font_size = 9)
knitr::kable(df.same_edges, row.names = TRUE, escape = F, table.attr = "style='width:30%;'",
             caption = "Count of common edges between decorations") %>%
  column_spec(1, bold=TRUE) %>%
  kableExtra::kable_styling(position = "center", font_size = 12)
```

In this dataframe:

* cells show the total number of common edges by decoration  

* the diagonal of the dataframe shows is the total number of edges of a given decoration  

Here, the decoration 2 has fifteen (15) edges and shares three (3) common edges with the decoration 3. To show them, and the decoration 4, we use the `list_compar()` function on the same variable (`type`) and the `plot_compar()` function with an edge focus (`focus = "edges"`).

This matrix can be used for further [clustering analysis](#sum)  

```{r compare.2.edges, out.width="100%", fig.width=12, fig.asp=0.52, fig.align="center", warning=FALSE}
dec.to.compare <- c(2, 3, 4)
g.compar <- list_compar(lgrph, nd.var = "type")
plot_compar(listg = g.compar, 
            dec2comp = dec.to.compare,
            focus = "edges",
            nd.size = c(0.5, 1.7),
            dir = dataDir)
```
  
## Contemporaneous elements {#contemp}

At times some nodes are non-contemporaneous with one another, like on the Ibahernando stelae. This stelae was reused as a funerary stelae during Roman times with the addition of a Latin inscription "*Alloquiu protaeidi.f hece. stitus*": Alluquio, son of Protacido, lies here [@Almagro66b]. 

<center>

![GIS view. The Ibahernando stelae (decoration 5)](`r ibahernando.path`){width=400px}

</center>

The writing (`ecriture`, node 1) has been carved over a spear (`lance`, node 2) and overlaps partially a V-notched shield (`bouclier`, node 3). The edges between node 1 and node 2, and the edge between node 1 and node 3, are *diachronic* edges  

```{r ibahernando, out.width="60%", fig.width=6, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae (decoration 5) with *diachronic* and *normal* edges, node 1 overlaps node 2 and node 3"}
site <- "Ibahernando" 
decor <- "Ibahernando"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              lbl.size = 0.7)
```
  
In this case, the non-contemporaneous layers of decoration, both nodes and edges, should be removed before the comparison process. To that purpose, the original graph (1-component) can be split into different contemporaneous sub-graphs. By removing the [*diachronic* edges](#ed.type.over) (`->-`), the graph is split into connected components, each including the synchronic nodes of a different period. The studied graph component will be retrieved with the component membership of a selected node.  
  
To study only the Late Bronze Age iconographic layer of the Ibahernando stelae, we can choose the Late Bronze Age node 4, the image of a sword (`epee`) dated to the middle and final stages of Late Bronze Age (ca 1250-950 BC). This node is believed to be contemporaneous with the spear (`lance`, node 2) and the shield (`bouclier`, node 3) so these three nodes are linked with [*normal* edges](#ed.type.norm) (`-=-`). We pass the node 4 to the parameters of the `contemp_nds()` function: 

```{r rm.writing, out.width="100%", fig.width=12, fig.asp=0.55, fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae before and after the selection of node 4 (sword) graph component"}
site <- decor <- "Ibahernando"
selected.nd <- 4
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
l_dec_df <- contemp_nds(nds.df, eds.df, selected.nd)
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))
par(mfrow=c(1, 2))
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = "type",
              lbl.color = "brown", lbl.size = 0.6)
plot_dec_grph(l_dec_df$nodes, l_dec_df$edges, imgs,
              site, decor, dataDir,
              nd.var = "type",
              lbl.color = "brown", lbl.size = 0.6)
```
  
On the other hand, epigraphists will only want to study the iconographic layer with the Latin writing. By selecting node 1 (`ecriture`), only the graph component of this node will be retained:

```{r ibahernando.lat, out.width="60%", fig.width=6, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae after the selection of node 1 (writing) graph component"}
selected.nd <- 1
nds.df <- read_nds(site, decor, dir = dataDir)
eds.df <- read_eds(site, decor, dir = dataDir)
l_dec_df <- contemp_nds(nds.df, eds.df, selected.nd)
plot_dec_grph(l_dec_df$nodes, l_dec_df$edges, imgs,
              site, decor, dir = dataDir,
              nd.var = "type",
              lbl.size = 0.6, lbl.color = "brown")
```
  
# Classify decorations {#sum}

Besides graphical functions allowing one to highlight common elements (nodes and edges) between decorations, the package also allows one to prepare data for unsupervised classification, such as hierarchical clustering with the `dist()` and `hclust()`:

```{r clust.comp, warning=FALSE, fig.align="center", fig.width=7, fig.height=5}
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))           
par(mfrow=c(1, 2))
df.same_edges <- same_elements(lgrph, "type", "edges")
df.same_nodes<- same_elements(lgrph, "type", "nodes")
dist.nodes <- dist(as.matrix(df.same_nodes), method = "euclidean")
dist.edges <- dist(as.matrix(df.same_edges), method = "euclidean")
hc.nds <- hclust(dist.nodes, method = "ward.D")
hc.eds <- hclust(dist.edges, method = "ward.D") 
plot(hc.nds, main = "Common nodes", cex = .8)
plot(hc.eds, main = "Common edges", cex = .8)
```

Clustering of decorations on common nodes and  clustering on common edges can be directly compared to one another:

```{r hclust.compar, warning=FALSE, fig.align="center", fig.width=7}
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(dplyr))
oldpar <- par(no.readonly = TRUE)   
on.exit(par(oldpar))           
par(mfrow=c(1, 2))
dend.nds <- as.dendrogram(hc.nds)
dend.eds <- as.dendrogram(hc.eds)
dendlist(dend.nds, dend.eds) %>%
  untangle(method = "step1side") %>% 
  tanglegram(columns_width = c(6, 1, 6),
             main_left = "Common nodes",
             main_right = "Common edges",
             lab.cex = 1.3,
             cex_main = 1.5,
             highlight_branches_lwd = F) 
```
  
In both clusterings, the Brozas stelae (decoration 3) and the Ibahernando stelae (decoration 5) are the ones having the most important proximities (ie, the least Euclidian distance). The 'Common edges' clustering is more accurate than the 'Common nodes' clustering because the former takes into account the common combination, or common permutation, of two nodes with the type of their edge, while the latter only takes into account the presence of common nodes.

```{r compare.c.edges, out.width="100%", fig.width=12, fig.asp=0.52, fig.align="center", warning=FALSE}
dec.to.compare <- c(3, 5)
g.compar <- list_compar(lgrph, nd.var = "type")
plot_compar(listg = g.compar, 
            dec2comp = dec.to.compare,
            focus = "nodes",
            nd.size = c(0.5, 1.7),
            dir = dataDir)
plot_compar(listg = g.compar, 
            dec2comp = dec.to.compare,
            focus = "edges",
            nd.size = c(0.5, 1.7),
            dir = dataDir)
```
  
  
# References


---
title: "Next developments"
author: "Thomas Huet, Jose M Pozo, Craig Alexander"
bibliography: references.bib
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
---

```{r, echo=FALSE}
# print(getwd())
# print(file.path("../logo", "iconr_logo.png"))
htmltools::img(src = knitr::image_uri(file.path("../man/figures", "iconr_logo.png")),
# htmltools::img(src = knitr::image_uri(file.path("img", "iconr_logo.png")), 
# htmltools::img(src = knitr::image_uri(file.path(R.home("doc"), "html", "logo.jpg")), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;')
```


<style>
.figure {
   margin-top: 0px;
   margin-bottom: 40px;
}
table {
    margin-top: 0px;
    margin-bottom: 24px;
}
</style>


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(kableExtra)
library(dplyr)
library(igraph)
library(iconr)
# library(ggdag)
library(tidygraph)
library(ggraph)
library(png)
library(grid)
library(gridExtra)
library(stratigraphr)

ibahernando.path <- "../man/figures/ibahernando_256colours.png"
brozas.path <- "../man/figures/brozas_256colours.png"
dummies.path <- "../man/figures/dummies.png"
ibahernando <-  rasterGrob(as.raster(readPNG(ibahernando.path)), interpolate = FALSE)
```

***iconr*** package is hope to lay the foundation of further developments

# Multifactorial analysis

In those examples, the study and comparisons have been made on the basis of different types of GUs with the variable `type` (`nd.var = "type"`). However, if a new column is added to the node dataframe or shapefile, the study can also incorporate other variables. Fox example, one of these new variable could be the technique by which the GU was created (`nd.var = "technique"`): 

<center>

![GIS view. The Brozas stelae (decoration 1) with the different engraving techniques: `g_piq` (peckings) for the spear, the shield and the sword; `g_inc` (incisions) for the the fibula and the comb](`r brozas.path`){width=400px}

</center>  

# Shape analysis

When the decoration iconographical content is composed by GUs separated one with another with neutral background, is could be easy to binarize and polygonize them (black GUs on white background). GUs can also be binarized/polygonized separately and reintegrated to the GIS interface

<center>

![GIS view. Graph decoration with shape analysis indexes (ConveHull, Minimum Bound Rectangle, Minimum Bound Circle)](`r dummies.path`){width=500px}

</center>

> Shape analysis allows to calculate numerous indexes (e.g., convexhull, MBR, MBC, Feret diameter) and/or contour comparisons (e.g., Procrustes analysis) that can be used in multivariate analysis to quantify the (dis)similarities between GUs 

# Typology

Graph theory and structured vocabularies allows one to construct tree structures for categorical variables (e.g., the different types of GUs). These structures allow generalization processes (upward to the parent level) and specification processes (downward to the child level). For example, a `sword` and a `spear` belong both to the `weapons` group (sub-group `offensive weapons`), a `shield` belongs to the `weapons` group (sub-group `defensive weapons`), etc.: 

```{r hierac, warning=FALSE, fig.width=7, fig.height=4, asp=0}
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))
par(mar=c(0, 0, 0, 0))
g <- graph_from_literal(objects-+weapons,
                        objects-+personnal_item,
                        weapons-+offensive_weapons,
                        weapons-+defensive_weapons,
                        offensive_weapons-+spear,
                        offensive_weapons-+sword,
                        defensive_weapons-+shield,
                        defensive_weapons-+helmet,
                        personnal_item-+miror,
                        personnal_item-+comb)
layout <- layout.reingold.tilford(g)
plot(g,
     layout = layout,
     vertex.color = "white",
     vertex.frame.color = "white",
     vertex.size = 20,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     edge.arrow.size = 0.5
     )
```
  
Such a formalism can be used to weight the differences between nodes, to conduct analysis with different levels of precision or to overcome issues of idiosyncratic typologies.

> A major issue will be to manage multi-linguism shared vocabularies to describe the iconographic content, the techniques used, etc.

# Harris matrix

Using the *diachronic* edge notation, layered directed acyclic graphs (DAG) -- like the Harris matrix relative chronology diagrams -- can also be used to construct. For example, with the Ibahernando stelae (decoration 5)

```{r overlap, echo=FALSE, warning=FALSE, out.width="100%"}
# library(data.tree)
# 
# imgs <- read.table(system.file("extdata", "imgs.csv", package = "iconr"),
#                    sep=";", stringsAsFactors = FALSE)
# nodes <- read.table(system.file("extdata", "nodes.csv", package = "iconr"),
#                     sep=";", stringsAsFactors = FALSE)
# edges <- read.table(system.file("extdata", "edges.csv", package = "iconr"),
#                     sep=";", stringsAsFactors = FALSE)
# lgrph <- list_dec(imgs, nodes, edges)
# edges.iba <- igraph::as_data_frame(lgrph[[5]], what="edges")
# overlap.nodes <- unique(edges.iba[edges.iba$type == ">", "from"])
# contemp.nodes <- unique(unlist(edges.iba[edges.iba$type == "=", c("from", "to")]))
# df.stratig <- data.frame(over = rep(overlap.nodes, length(contemp.nodes)),
#                          under = contemp.nodes)
# df.stratig$pathString <- paste(lgrph[[5]]$decor,
#                                df.stratig$over, 
#                                df.stratig$under, 
#                                sep = "/")
# # superpo <- as.Node(df.stratig)
# # print(superpo)
# gd <- graph_from_data_frame(df.stratig, directed = TRUE, vertices = NULL)
# gd.ed <- as_data_frame(gd, what="edges")
# gd.nd <- as_data_frame(gd, what="vertices")
# rstat_nodes <- data.frame(name = gd.nd)
# rstat_edges <- data.frame(from = gd.ed$from,
#                           to = gd.ed$to)
# # rstat_nodes <- data.frame(name = c("Hadley", "David", "Romain", "Julia"))
# # rstat_edges <- data.frame(from = c(1, 1, 1, 2, 3, 3, 4, 4, 4),
# #                           to = c(2, 3, 4, 1, 1, 2, 1, 2, 3))
# gd.all <- tbl_graph(nodes = rstat_nodes, edges = rstat_edges)
# gd.all <- set.vertex.attribute(gd.all, "context", value=as.character(V(gd.all)))
imgs <- read.table(system.file("extdata", "imgs.csv", package = "iconr"),
                   sep=";", stringsAsFactors = FALSE)
nodes <- read.table(system.file("extdata", "nodes.csv", package = "iconr"),
                    sep=";", stringsAsFactors = FALSE)
edges <- read.table(system.file("extdata", "edges.csv", package = "iconr"),
                    sep=";", stringsAsFactors = FALSE)
lgrph <- list_dec(imgs, nodes, edges)
# edges
edges.iba <- igraph::as_data_frame(lgrph[[5]], what="edges")
over.edges <- edges.iba[edges.iba[, "type"] == ">", ]
contemp.edges <- edges.iba[edges.iba[, "type"] == "=", ]
# nodes
nodes.iba <- igraph::as_data_frame(lgrph[[5]], what="vertices")
tib.gd.nd <- as_tibble(nodes.iba) # convert nodes df to tibble
tib.gd.nd$context <- tib.gd.nd$name # the "context"
# prepare the df
tib.gd.nd$above <- tib.gd.nd$below <- tib.gd.nd$equal <- NA
class(tib.gd.nd$above) <- class(tib.gd.nd$below) <- class(tib.gd.nd$equal) <- "list" # change class
tib.gd.nd$name <- tib.gd.nd$x <- tib.gd.nd$y <- tib.gd.nd$type <- NULL
# OVER
for(n in tib.gd.nd$context){
  # filter for each nodes
  n.is.above <- over.edges[over.edges[, "from"] == n, ]
  if (nrow(n.is.above) > 0){
    tib.gd.nd[n , "below"][[1]] <- list(n.is.above[, "to"])
  }
  n.is.below <- over.edges[over.edges[, "to"] == n, ]
  if (nrow(n.is.below) > 0){
    tib.gd.nd[n , "above"][[1]] <- list(n.is.below[, "from"])
  }
  n.is.contemp <- contemp.edges[contemp.edges[, "from"] == n, ]
  if (nrow(n.is.contemp) > 0){
    tib.gd.nd[n , "equal"][[1]] <- list(n.is.contemp[, "to"])
  }
}

## see equal nodes
# equals are listed in "equals" and "context"
eqs <-  tib.gd.nd[as.vector(!is.na(tib.gd.nd[, "equal"])), ]
# /!\ what if different layers ?
nds.equals <- unlist(unique(c(eqs$context, eqs$equal)))
get.eq <- get.bl <- get.ab <- list()
# df <- tib.gd.nd
for (n in nds.equals){
  # n <- 2
  # get 
  other.equal <- nds.equals[nds.equals != n]
  tib.gd.nd[n, "equal"][[1]] <- list(other.equal)
}
## merge values for equal nodes
# see above for equal nodes
ab <- unique(as.vector(unlist(tib.gd.nd[nds.equals, "above"])))
ab <- ab[!is.na(ab)]
if(is.logical(ab)) ab <- NA
tib.gd.nd[nds.equals, "above"][[1]] <- list(ab)
# see below for equal nodes
bl <- unique(as.vector(unlist(tib.gd.nd[nds.equals, "below"])))
bl <- bl[!is.na(bl)]
# avoid logical(0)
if(is.logical(bl)) bl <- NA
tib.gd.nd[nds.equals, "below"][[1]] <- list(bl)
## merge for overlap nodes (non equal)
for(i in nrow(tib.gd.nd)){
  # i <- 1
  ov <- unique(as.vector(unlist(tib.gd.nd[i, "below"])))
  node.to.add <- setdiff(nds.equals, ov)
  is.above <- unique(as.vector(unlist(tib.gd.nd[node.to.add, "above"])))
}
#
all.above <- as.vector(unlist(tib.gd.nd[, "above"]))
for(i in 1:length(all.above)){
  # i <- 2
  # ctx <- all.above[i]
  context <- as.vector(unlist(tib.gd.nd[i, "context"]))
  # find the nodes where the context is above
  abv <- as.character(which(all.above %in% context))
  # avoid character(0)
  abv[length(abv) == 0] <- NA
  tib.gd.nd[i, "below"][[1]] <- list(abv)
}
# add "natural" at the bottom
tib.gd.nd[nrow(tib.gd.nd)+1, "context"] <- "natural"
# select nodes without any nodes below them
down.nodes <- as.vector(unlist(tib.gd.nd[is.na(tib.gd.nd[, "below"]), "context"]))
# the down.nodes are upper the "natural"
tib.gd.nd[nrow(tib.gd.nd), "above"][[1]] <- list(down.nodes)
# add "natural" below the down.nodes
tib.gd.nd[down.nodes, "below"][[1]] <- list("natural")
# - - - - - - - - - - - - - -
# gd.all <- tbl_graph(nodes = tib.gd.nd, edges = edges.iba)
h12_graph <- stratigraph(tib.gd.nd, "context", "above") # works harris12, ~ with 'gd.all'
# ggraph(h12_graph, layout = "sugiyama") +
#   geom_edge_elbow() +
#   geom_node_label(aes(label = context), label.r = unit(0, "mm")) +
#   theme_graph()
a.g <- ggraph(h12_graph, layout = "sugiyama") +
  geom_edge_elbow() +
  geom_node_label(aes(label = context), label.r = unit(0, "mm")) +
  theme_graph()
# knitr::include_graphics(c("path/to/img1","path/to/img1"))
grid.arrange(ibahernando, a.g, ncol = 2)
```
  
> DAG is a common graphical formalism to represent chronological ordering such as Bayesian, radiocarbon and age-depth models (e.g., [RChronoModel](https://cran.r-project.org/web/packages/RChronoModel/index.html), [stratigraphr](https://stratigraphr.joeroe.io/articles/stratigraph.html)). DAGs can also be used for cultural period sequences like GUs chrono-cultural attribution (e.g., "Late Bronze Age", "Iron Age", "Roman Times")    
---
title: "The ***iconr*** package. Graph Analysis of Prehistoric Iconography with R"
author: "Thomas Huet, Jose M Pozo, Craig Alexander"
email: "thomashuet7@gmail.com"
# date: "`r format(Sys.Date())`"
bibliography: references.bib
## to create the vignettes 'outside' the package, with table of content (toc)
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
---

```{r, echo=FALSE}
# print(getwd())
# print(file.path("../logo", "iconr_logo.png"))
htmltools::img(src = knitr::image_uri(file.path("tutorial/img", "iconr_logo.png")), 
# htmltools::img(src = knitr::image_uri(file.path(R.home("doc"), "html", "logo.jpg")), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;')
```

---

<style>
.figure {
   margin-top: 0px;
   margin-bottom: 20px;
}
table {
    margin-top: 0px;
    margin-bottom: 24px;
}
</style>

```{r, include = FALSE}
library(knitr)
library(igraph)
library(dplyr)
library(kableExtra)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">",
  fig.pos = 'H'
)
ibahernando.path <- paste0(getwd(),"/img/ibahernando.png")
brozas.path <- paste0(getwd(),"/img/brozas.png")
solanas.path <- paste0(getwd(),"/img/solana_detail.png")
solanas.vor.path <- paste0(getwd(),"/img/solana_voronoi.png")
```
  
The R package `iconr` is grounded in graph theory and spatial analysis. It offers concepts and functions for modeling Prehistoric iconographic compositions and for their preparation for further analysis (clustering, Harris diagram, etc.) in order to contribute to cross-cultural iconography comparison studies through a greater normalization of quantitative analysis [@Alexander08; @HuetAlexander15; @Huet18a].  
  
  
<center>

![Detail of a Late Bronze Age stelae (Solana de Cabañas, Cáceres, Spain). Credits: Museo Arqueológico Nacional, Madrid](`r solanas.path`){width=220px}

</center>
  
The flexibility of graph theory and tools available for the GIS database make the `iconr` package useful in managing, plotting and comparing (potentially large) sets of iconographic content: Atlantic rock-art, Scandinavian rock art, Late Bronze Age stelae, Mycenean figurative pottery, etc. 
  
# Decoration graphs  

The main principle of the `iconr` package is to consider any iconographic composition (here, 'decoration') as a geometric graph of graphical units (GUs). This geometric graph is also known as a planar graph or spatialized graph. The GUs are decorated surfaces (`POLYGONS`) modeled as nodes (`POINTS`). When these GUs are *main* nodes, and not *attribute* nodes, they share edges (`LINES`) with one another when their Voronoi cells share a border (*birel*: touches).  

<center>

![GIS view. The Solana de Cabanas stelae: from its photograph to the modeling of its graphical content](`r solanas.vor.path`){width=800px}

</center>  
  
  
Graph theory offers a conceptual framework and indices (global at the entire graph scale, local at the vertex scale) to deal with notions of networks, relationships and neighbourhoods. The geometric graph is commonly built within a GIS interface. Indeed, use of GIS allows one to create a spatial database of the decoration's iconographic contents and facilitates data recording and visualization. For example, snapping options can connect GUs (nodes) with lines (edges) and we can exploit tools such as feature symbology, layer transparency, etc.  
  
The latest development version of the `iconr` package and its vignette can be downloaded from GitHub

```{r down,eval=FALSE, echo=TRUE}
devtools::install_github("zoometh/iconr", build_vignettes=TRUE)
```

The R package `iconr` is composed of [functions](#functions) and a [example dataset](#data). The main R packages used by the `iconr` package are:

* [magick-image](https://CRAN.R-project.org/package=magick): for drawing/image management 
* [igraph](https://CRAN.R-project.org/package=igraph): for graph management 
* [rgdal](https://CRAN.R-project.org/package=rgdal): for shapefile management

Load the package `iconr`

```{r load, echo=TRUE}
library(iconr)
``` 

# Dataset {#data}

The input dataset is expected to include decoration images and corresponding node and edge data in a single data folder. This folder should include the following files:

* **[Table of decorations](#decorations):** A tabular file storing the set of decoration identifiers and corresponding image filenames.
* **[Images](#drawings):** An image file for each decoration.
* **[Node data](#nd):** A single file storing the data of each node for all decorations. 
* **[Edge data](#ed):** A single file storing the data of each edge for all decorations.

The `iconr` package includes an example dataset with the input files in several alternative formats. The path for the example dataset is the package *extdata* folder. This folder is also -- by default -- the output folder. In order to differentiate input data from output data, the output data filenames always start with a digit or an underscore (ie, a punctuation). So, the input data are :

```{r ls_ext_data}
dataDir <- system.file("extdata", package = "iconr")
all.files <- list.files(dataDir)
input.files <- all.files [!grepl("^[[:digit:]]", all.files) &
                            !grepl("^[[:punct:]]", all.files)]
cat(input.files, sep="\n")
```

The table of decorations is given in two formats: comma- or semicolon-separated values (`imgs.csv`) and tab-separated value (`imgs.tsv`).

```{r paths.imgs, echo=TRUE}
imgs_path <- paste0(dataDir, "/imgs.csv")
imgs <- read.table(imgs_path, sep=";", stringsAsFactors = FALSE)
```

Each decoration is identified by its name (column `decor`) and the name of the site (column `site`) to which it belongs. In the example dataset, this is transparent in the name of each decoration image, included in jpg format. Any other image format supported by the R package `magick` (jpg, png, tiff, pdf, etc.) is suitable.

As we have stated, a GIS interface is often the most practical way to record graph nodes and graph edges with `POINTS` and `LINES` geometries, respectively. This is typically saved in shapefile (shp) format, which is composed of at least 3 files with extensions `.shp` (geometries), `.shx` (indices), and `dbf` (attribute data). The example dataset includes them for nodes and edges separately, with obvious names:

```{r paths, echo=TRUE}
nodes_path <- paste0(dataDir, "/nodes.shp")
nodes.shp <- rgdal::readOGR(dsn = nodes_path, verbose = FALSE)
nodes <- as.data.frame(nodes.shp)
edges_path <- paste0(dataDir, "/edges.shp")
edges.shp <- rgdal::readOGR(dsn = edges_path, verbose = FALSE)
edges <- as.data.frame(edges.shp)
```

Nodes and edges can also be recorded in tabular format: `.csv` or `.tsv`.  

```{r paths.1, echo=TRUE}
nodes_path <- paste0(dataDir, "/nodes.tsv")
nodes <- read.table(nodes_path, sep="\t", stringsAsFactors = FALSE)
edges_path <- paste0(dataDir, "/edges.tsv")
edges <- read.table(edges_path, sep="\t", stringsAsFactors = FALSE)
```

The list of graph decorations is created with the `list_dec()` function. These graphs are `igraph` objects

```{r graph.clss}
lgrph <- list_dec(imgs, nodes, edges)
g <- lgrph[[1]]
as.character(class(g))
```

By default, the `plot.igraph()` function (ie, `igraph::plot()`) spatialization (`layout`) is based on `x` and `y` columns, when these exist. This is appropriate in our case since we are working with geometric graphs. If this spatialisation (node coordinates) is ignored, then very different layouts are possible for the same graph (*graph drawing*):

```{r igraph.1, warning=FALSE, fig.align="center", fig.width=6.5, fig.asp=0.58}
par(mar=c(1,0,2,0), mfrow=c(1, 2), cex.main=0.9, font.main= 1)
coords <- layout.fruchterman.reingold(lgrph[[1]])
plot(g,
     vertex.size = 15,
     vertex.frame.color="white",
     vertex.label.family = "sans",
     vertex.label.cex = .8,
     main = "Graph drawing based on x, y coordinates"
)
plot(g,
     layout = layout.fruchterman.reingold(g),
     vertex.size = 5 + degree(g)*10,
     vertex.frame.color="white",
     vertex.label.family = "sans",
     vertex.label.cex = .8,
     main = "Force-directed graph drawing,\nwith degree-dependent node size."
)
mtext(g$decor, cex = 1, side = 1, line = -1, outer = TRUE)
```
  
## Table of decorations {#decorations}

In any compatible file format readable by R, such as comma-separated values (csv) or tabular-separated values (tsv). Decorations identifiers and decorations image filenames are stored in a dataframe, by default the `imgs` dataframe:

```{r imgs,fig.width=6, fig.height=6, fig.align="center",warning=FALSE, fig.cap="\\label{fig:figs}imgs.tsv"}
imgs_path <- paste0(dataDir, "/imgs.tsv")
imgs <- read.table(imgs_path, sep="\t", stringsAsFactors = FALSE)
knitr::kable(imgs, "html") %>% 
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12)
```
  
The decoration's unique identifiers are the concatenation of the site name and the decoration name. For example, the name of the Cerrano Muriano 1 decoration is: ``r imgs[1,"img"]``. 

## Images {#drawings}

Image or drawing (eg, images, rasters, grids) formats accepted are the common ones (jpg, png, jpeg, tiff, pdf, etc.). The images in the current example dataset come from a PhD thesis, published by M. Diaz-Guardamino [@DiazGuardamino10]. 

For a given decoration, its image is the reference space of the graph: nodes and edges inherit their coordinates from this image. Firstly, the decoration image is open in the GIS and the graph is build. The, these features are imported into R:

```{r xy_coords,out.width="50%", fig.align="center",echo=FALSE,warning=FALSE}
df.equi <- data.frame("devices" = c("GIS raster","R bitmap image"),
                      "units of measures" = c("number of pixels", "number of pixels"),
                      "origin corner" = c("top-left", "top-left"),
                      "x-axis" = c("rightward", "rightward"),
                      "y-axis" = c("upward", "downward"),
                      check.names = FALSE)
knitr::kable(df.equi) %>%
  kableExtra::add_header_above(c(" " = 1, " " = 1, " " = 1, "axes orientation" = 2)) %>%
  kable_styling(full_width = F)
```

Observe that the y-axis orientation are different between a GIS raster/grid and R bitmap images. In the GIS the y-axis is oriented upwards (ie, with negative values downwards) and y-coordinates are always negative. However, it contrast with R bitmap image (top-left origin but downwards y-axis), for which `y` values are positive. Then, the heights (ie, y-axis) coming from the GIS have to be inverted. This can be illustrated with the `magick` package.

```{r drawing, out.width="50%", fig.align="center", warning=FALSE, echo=TRUE, message=FALSE, fig.cap="\\label{fig:figs}R grid conformed to the GIS coordinates: `Cerro_Muriano.Cerro_Muriano_1.jpg` with the coordinates of its corners."}
library(magick)
dataDir <- system.file("extdata", package = "iconr")
imgs_path <- paste0(dataDir, "/imgs.csv")
imgs <- read.table(imgs_path, sep=";", stringsAsFactors = FALSE)
cm1 <- image_read(paste0(dataDir, "/", imgs$img[1]))
W <- image_info(cm1)$width
H <- -image_info(cm1)$height
image_border(image = cm1, "#808080", "2x2") %>%
  image_annotate(text = paste0(0, ", ", 0), size = 30,
                 gravity = "northwest") %>%
  image_annotate(text = paste0(W, ", ", 0), size = 30,
                 gravity = "northeast") %>%
  image_annotate(text = paste0(0, ", ", H), size = 30,
                 gravity = "southwest") %>%
  image_annotate(text = paste0(W, ", ", H), size = 30,
                 gravity = "southeast")
```

## Node data {#nd}

Nodes can be stored in any of the following formats: comma-separated values (`.csv`), tab-separated values (`.tsv`), or shapefile (`.shp`). Any of these formats is consistenly read by the function `iconr::read_nds()` generating a data.frame with the node data. For the generic csv and tsv formats, at least the following five columns are required to be included: 

* **site**: decoration site 
* **decor**: decoration name 
* **id**: id uniquelly identifying each node
* **x**: x-coordinate of nodes
* **y**: y-coordinate of nodes 

Additional columns can be included providing relevant characteristics of each node. The example dataset includes the column:

* **type**: type of GU (ie, what it represents)

These node characteristics can then be used for image annotations, decoration comparisons and decoration analysis. Nodes attributes coming from the GIS interface:

```{r nodes.df, warning=FALSE,fig.align="center",warning=FALSE}
nds.df <- read_nds(site = "Cerro Muriano", decor = "Cerro Muriano 1", dir = dataDir) 
knitr::kable(nds.df, "html") %>% 
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12)
```
  
In principle, the nodes are defined to be located at the exact centroids of each GU. However, for practical reasons, they can also be manually located close to these centroids.  

### Nodes types {#nd.type}

Nodes can be *main* nodes or *attribute* nodes depending on the type of edges they share with other nodes (see [edge types](#ed.type)). For example, considering the Cerro Muriano 1 stelae:

* **main nodes**: the anthropomorph, the spear, the shield, the comb and the ingot can be considered as *main* nodes because they are separated and they show different graphical contents.

* **attribute nodes**: the helmet and the male sex of the antropomorphic figure can be considered as *attribute* nodes of the anthropomorph because they are characteristic features of this latter.

## Edge data {#ed}

Edges can be stored in the same three formats as nodes: csv, tsv, or shp. At least the following five columns or variables are required to be included:

* **site**: decoration site 
* **decor**: decoration name
* **a**: *id* of the first [node](#nd) of each edge
* **b**: *id* of the second [node](#nd) of each edge
* **type**: [edge types](#ed.type)

These columns can be seen in the example dataset:

```{r edges.df, warning=FALSE}
edges <- read.table(edges_path, sep = "\t", stringsAsFactors = FALSE)
knitr::kable(head(edges), "html", align = "llccc") %>%
  kableExtra::kable_styling(full_width = FALSE, position = "center",
                            font_size=12) %>%
  gsub("\\+", "$+$", .)
```
  
Analogously to nodes, any of these formats is consistently read by the function `iconr::read_eds()` generating a data.frame with the edges data. Additionally, this function includes columns with the coordinates of the edge nodes:

* **xa**: x-coordinate of the first node of each edge
* **ya**: y-coordinate of the first node of each edge 
* **xb**: x-coordinate of the second node of each edge
* **yb**: y-coordinate of the second node of each edge 

Those coordinates are imported from the `x` and `y` coordinates in the corresponding [node](#nd) file, identified by the node `id` as stated in the edge fields `a` and `b`, respectively.

```{r edges.df.1, warning=FALSE}
eds.df <- read_eds(site = "Cerro Muriano", decor = "Cerro Muriano 1", dir = dataDir) 
knitr::kable(head(eds.df), "html", align = "llcccrrrr") %>% 
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12) %>%
  gsub("\\+", "$+$", .)
```
  
**column names** 

* `type` : [edges types](#ed.type)  

  * `=` : *normal* edges between contiguous and contemporaneous *main* nodes (undirected edge)  

  * `+` : *attribute* edges, between contemporaneous nodes where the *attribute* node `b` is an attribute of *main* node `a` (directed edge)

  * `>` : *diachronic* edges, between non-contemporaneous nodes where the node `a` overlaps node `b`, or node `a` is more ancient than node `b` (directed edge)

* `xa, ya`: coordinates of the *starting* node, or *main* node, or *overlapping* node, or *more recent* node  (`a`) 

* `xb, yb`: coordinates of the *ending* node, or *attribute* node, or *overlapped* node, , or *more ancient* node (`b`) 
  
  
### Edge types {#ed.type}

Graph theory says that edges can be undirected or directed. In the `iconr` package, by default: 

* all contemporaneous nodes have [*normal* edges](#ed.type.norm) or [*attribute* edges](#ed.type.attrib)  edges displayed in <span style="color:orange"><b>orange</b></span>

* all non-comtemporaneous nodes have [*diachronic* edges](#ed.type.over) displayed in <span style="color:blue"><b>blue</b></span> edges (see [contemporaneous nodes](#contemp))

The `named_elements()` function allows one to display the textual notation of the different types of edges (`-=-`, `-+-` or `->-`)

```{r count4, warning=FALSE}
named_elements(lgrph[[1]], focus = "edges", nd.var="type")[1]      
```

When there are nodes with the same `nd.var`, this function adds the suffix `#` to the `nd.var` in order to disambiguate the node list. This is the case, for example, for the `chariot_char-+-cheval` (x2) and `chariot_char-+-roue` (x2) edges of the Zarza de Montsanchez stelae (decoration 4)

```{r count1, warning=FALSE}
named_elements(lgrph[[4]], focus = "edges", nd.var="type")
```

Employed with the basic R functions for loop (`lapply()`), count (`table()`) and order (`order()`), and removing the the suffix `#`, this function can be used to count the different types of edge. Here, we enumerate the most represented types of the example dataset:

```{r count.all, warning=FALSE}
all.edges <- unlist(lapply(lgrph, named_elements, 
                           focus = "edges", nd.var="type", disamb.marker=""))
edges.count <- as.data.frame(table(all.edges))
edges.order <- order(edges.count$Freq, decreasing = TRUE)
edges.count <- edges.count[edges.order, ] 
knitr::kable(head(edges.count), row.names = FALSE) %>%
  kableExtra::kable_styling(full_width = FALSE, position = "center", font_size=12)
```
  
#### ***Normal* edges** {#ed.type.norm}

Different *main* nodes considered to be contemporaneous and close to one another, may share an edge with the  value `=` for their `type`. By convention, these edges are called *normal* and displayed as a plain line. Its textual notation is `-=-`. The *normal* edges are undirected:  `1-=-2` is equal to `2-=-1`, where node 1 and node 2 are two different *main* nodes.  

#### ***Attribute* edges** {#ed.type.attrib}

When a node is an attribute of another, edges are identified with a `+` and displayed with a dashed line. For example, at the bottom of the Zarza de Montsanchez stelae (decoration 4), the main node 7 (`chariot`) is connected with four (4) attribute nodes:

- two horses (`cheval`): 8 and 9  

- two wheels (`roue`): 10 and 11 

```{r graph.attribute.plot.type, out.width="60%", fig.width=8, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Zarza De Montanchez stelae (decoration 4) showing *normal* and *attribute*  edges"}
site <- "Zarza de Montanchez" 
decor <- "Zarza De Montanchez"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              ed.lwd = 1, ed.color = c("darkorange"),
              lbl.size = 0.7)
```
  
The textual notation of an *attribute* edge is `-+-`:

```{r count_att, warning=FALSE}
sort(named_elements(lgrph[[4]], focus = "edges", nd.var = "type"))
```
  
The *attribute* edges are directed: `1-+-2` is not equal to `2-+-1`, `1-+-2` means that node 1 is the *main* node and node 2 is one of its *attribute* nodes. 
 
#### ***Diachronic* edges** {#ed.type.over}

When a node overlaps another or is more recent than another, edges are identified with a `>` and displayed with a blue plain line. For example, the Ibahernando stelae has a latin inscription (`ecriture`) overlapping a spear (`lance`) and a shield (`bouclier`).  

```{r graph.overlap.plot.type, out.width="100%", fig.width=16, fig.asp=0.55, fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae (decoration 5) showing *diachronic* and *normal* edges"}
site <- "Ibahernando"
decor <- "Ibahernando"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
par(mfrow = c(1,2))
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              lbl.size = 0.7)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = 'type',
              lbl.size = 0.6)
```
  
The textual notation of an *diachronic* edge is `->-`:

```{r count3, warning=FALSE}
named_elements(lgrph[[5]], focus = "edges", nd.var = "type")
```

The *diachronic* edges are directed: `1->-2` is not equal to `2->-1`. The `1->-2` edge means that node 1 overlaps node 2, or node 1 is more recent than node 2. These overlays, or diachronic iconographic layers, can be managed with the `contemp_nds()` function (see the section [Contemporaneous contents](#contemp))  

# Functions {#functions}

Decoration graphs are constructed from nodes and edges. Graphs are 1-component: each decoration graph covers all the GUs of the decoration. The functions of the `iconr` package provide basic tools to manage these [node data](#nd) (`.csv`, `.tsv` or `.shp`) and [edge data](#ed) (`.csv`, `.tsv` or `.shp`) to create graphs, to [plot](#plot) and to [compare](#compare) them, to select [contemporaneous GU compositions](#contemp)  

```{r ls_functions}
cat(ls("package:iconr"), sep="\n")
```

## Read {#read}

The functions `read_nds()` and `read_eds()` allow, respectively, to read a dataframe or a shapefile of nodes, and a dataframe or a shapefile of edges. For example, the nodes and edges of Cerro Muriano 1 can be read

```{r img.graph.read, echo=TRUE}
site <- "Cerro Muriano"
decor <- "Cerro Muriano 1"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
```

## Plot {#plot}

The graphical functions `plot_dec_grph()` and `plot_compar()` allow, respectively, to plot one decoration or pairwise(s) of decorations. These function offer different choices for the color and size of the nodes, edges or labels. For example, for the Cerro Muriano 1 decoration, the field `id` (identifier of the node), used by default for the labels, can be changed to the `type` field

```{r img.graph.plot.type, out.width="60%", fig.width=8, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Cerro Muriano 1 stelae (decoration 1) with the type of each GU"}
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = 'type',
              lbl.size = 0.55)
```
  
A new field, `long_cm`, is added to the Cerro Muriano 1 nodes and the graph is replotted using this field instead of the `type` field, with <span style="color:brown">brown</span> colors and larger labels.

```{r img.graph.plot.id, out.width="60%", fig.width=8, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Cerro Muriano 1 stelae (decoration 1) with the maximum length (in cm) of each GU"}
nds.df$long_cm <- paste0(c(47, 9, 47, 18, 7, 3, 13), "cm")
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = 'long_cm',
              nd.color = "brown",
              lbl.color = "brown", lbl.size = 0.7,
              ed.color = "brown")
```
  
## Compare {#compare}

Elements of the graphs (nodes and edges) can be compared across all graphs or between a pair of graphs with the `same_elements()` and `plot_compar()` functions

* `same_elements()` function permits one to count common elements between *n* graphs   

* `plot_compar()` function shows a graphical output for these common elements

By default, in a pairwise comparison of decorations, common nodes and edges are displayed in <span style="color:red">red</span>, but their colors -- and other graphical parameters -- can be modified. When not all GUs are contemporaneous with one another, the non-contemporaneous ones can be removed with the `contemp_nds()` function. 

### Node comparisons

A classic study in archaeological research is to count the common nodes between pairs of decorations. This can be done with the `same_elements()` function with a node focus (`focus = "nodes"`), and considering for example their `type` (by default).

```{r compare.nodes, results='asis', warning=FALSE}
imgs_path <- paste0(dataDir, "/imgs.tsv")
nodes_path <- paste0(dataDir, "/nodes.tsv")
edges_path <- paste0(dataDir, "/edges.tsv")
imgs <- read.table(imgs_path, sep="\t", stringsAsFactors = FALSE)
nodes <- read.table(nodes_path, sep="\t", stringsAsFactors = FALSE)
edges <- read.table(edges_path, sep="\t", stringsAsFactors = FALSE)
lgrph <- list_dec(imgs, nodes, edges)
df.same_nodes <- same_elements(lgrph,
                               focus = "nodes",
                               nd.var = "type")
diag(df.same_nodes) <- cell_spec(diag(df.same_nodes),
                                 font_size = 9)
knitr::kable(df.same_nodes, row.names = TRUE, escape = FALSE, table.attr = "style='width:30%;'",
             caption = "Count of common nodes between decorations") %>%
  column_spec(1, bold=TRUE) %>%
  kableExtra::kable_styling(position = "center", font_size = 12)
```
  
The result of `same_elements()` is a symmetric matrix giving the number of common nodes for each pair of decorations, where the row and column names are the identifiers of the decorations. Observe that, accordingly, the diagonal elements show the total number of nodes of each decoration.

Regarding the node variable `type`, the decoration 4 has a total twelve (12) nodes (see the diagonal). Decoration 4 and decoration 2 have nine (9) common nodes. Decoration 4 and decoration 3 have four (4) common nodes. This matrix can be used for further [clustering analysis](#sum)

To compare graphically the decorations 2, 3 and 4 on the node variable `type`:

* first: `type` variable is pasted to the `list_compar()` function  

* then: the plot is made with the `plot_compar()` function
  
```{r compare.2.nodes, fig.show = TRUE, out.width = "700px", fig.align="center", warning=FALSE}
dec.to.compare <- c(2, 3, 4)
g.compar <- list_compar(lgrph, nd.var = "type")
nds_compar <- plot_compar(listg = g.compar, 
                          dec2comp = dec.to.compare,
                          focus = "nodes",
                          nd.size = c(0.5, 1.5),
                          dir = dataDir,
                          img.format = "png")
knitr::include_graphics(nds_compar) 
```
  
The function creates an image for each pair of stelae contained in the `dec.to.compare` variable (`r dec.to.compare`), with a focus on nodes (`focus = "nodes"`). Thus, if $n$ decorations are compared, it results in $n\choose 2$ $\frac{n!}{(n-2)!2!}$ pairwise comparison images. For instance, if all 5 decorations in the example dataset are compared, there will be 10 pairwise comparisons.

### Edge comparisons

A less-common study in archaeological research is to count the common edges between pairs of decorations. The `same_elements()` function with an edge focus (`focus = "edges"`) and considering the  `type` of the nodes 

```{r compare.edges, warning=FALSE}
df.same_edges <- same_elements(lgrph, nd.var = "type", focus = "edges")
diag(df.same_edges) <- cell_spec(diag(df.same_edges),
                                 font_size = 9)
knitr::kable(df.same_edges, row.names = TRUE, escape = F, table.attr = "style='width:30%;'",
             caption = "Count of common edges between decorations") %>%
  column_spec(1, bold=TRUE) %>%
  kableExtra::kable_styling(position = "center", font_size = 12)
```

In this dataframe:

* cells show the total number of common edges by decoration  

* the diagonal of the dataframe shows is the total number of edges of a given decoration  

Here, the decoration 2 has fifteen (15) edges and shares three (3) common edges with the decoration 3. To show them, and the decoration 4, we use the `list_compar()` function on the same variable (`type`) and the `plot_compar()` function with an edge focus (`focus = "edges"`).

This matrix can be used for further [clustering analysis](#sum)  

```{r compare.2.edges, out.width = "700px", fig.align="center", warning=FALSE}
dec.to.compare <- c(2, 3, 4)
g.compar <- list_compar(lgrph, nd.var = "type")
eds_compar <- plot_compar(listg = g.compar, 
                          dec2comp = dec.to.compare,
                          focus = "edges",
                          nd.size = c(0.5, 1.7),
                          dir = dataDir,
                          img.format = "png")
knitr::include_graphics(eds_compar) 
```
  
## Contemporaneous elements {#contemp}

At times some nodes are non-contemporaneous with one another, like on the Ibahernando stelae. This stelae was reused as a funerary stelae during Roman times with the addition of a Latin inscription "*Alloquiu protaeidi.f hece. stitus*": Alluquio, son of Protacido, lies here [@Almagro66b].

<center>

![GIS view. The Ibahernando stelae (decoration 5)](`r ibahernando.path`){width=400px}

</center>

The writing (`ecriture`, node 1) has been carved over a spear (`lance`, node 2) and overlaps partially a V-notched shield (`bouclier`, node 3). The edges between node 1 and node 2, and the edge between node 1 and node 3, are *diachronic* edges  

```{r ibahernando, out.width="60%", fig.width=8, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae (decoration 5) with *diachronic* and *normal* edges, node 1 overlaps node 2 and node 3"}
site <- "Ibahernando" 
decor <- "Ibahernando"
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              lbl.size = 0.7)
```
  
In this case, the non-contemporaneous layers of decoration, both nodes and edges, should be removed before the comparison process. To that purpose, the original graph (1-component) can be split into different contemporaneous sub-graphs. By removing the [*diachronic* edges](#ed.type.over) (`->-`), the graph is split into connected components, each including the synchronic nodes of a different period. The studied graph component will be retrieved with the component membership of a selected node.  
  
To study only the Late Bronze Age iconographic layer of the Ibahernando stelae, we can choose the Late Bronze Age node 4, the image of a sword (`epee`) dated to the middle and final stages of Late Bronze Age (ca 1250-950 BC). This node is believed to be contemporaneous with the spear (`lance`, node 2) and the shield (`bouclier`, node 3) so these three nodes are linked with [*normal* edges](#ed.type.norm) (`-=-`). We pass the node 4 to the parameters of the `contemp_nds()` function: 

```{r rm.writing, out.width="100%", fig.width=16, fig.asp=0.55, fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae before and after the selection of node 4 (sword) graph component"}
site <- decor <- "Ibahernando"
selected.nd <- 4
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
l_dec_df <- contemp_nds(nds.df, eds.df, selected.nd)
par(mfrow=c(1,2))
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor, dataDir,
              nd.var = "type",
              lbl.color = "brown", lbl.size = 0.6)
plot_dec_grph(l_dec_df$nodes, l_dec_df$edges, imgs,
              site, decor, dataDir,
              nd.var = "type",
              lbl.color = "brown", lbl.size = 0.6)
```
  
On the other hand, epigraphists will only want to study the iconographic layer with the Latin writing. By selecting node 1 (`ecriture`), only the graph component of this node will be retained:

```{r ibahernando.lat, out.width="60%", fig.align="center", warning=FALSE, fig.cap="Ibahernando stelae after the selection of node 1 (writing) graph component"}
selected.nd <- 1
nds.df <- read_nds(site, decor, dir = dataDir)
eds.df <- read_eds(site, decor, dir = dataDir)
l_dec_df <- contemp_nds(nds.df, eds.df, selected.nd)
plot_dec_grph(l_dec_df$nodes, l_dec_df$edges, imgs,
              site, decor, dir = dataDir,
              nd.var = "type",
              lbl.size = 0.6, lbl.color = "brown")
```
  
# Classify decorations {#sum}

Besides graphical functions allowing one to highlight common elements (nodes and edges) between decorations, the package also allows one to prepare data for unsupervised classification, such as hierarchical clustering with the `dist()` and `hclust()`:

```{r clust.comp, warning=FALSE, fig.align="center", fig.width=7, fig.height=5}
par(mfrow=c(1, 2))
df.same_edges <- same_elements(lgrph, "type", "edges")
df.same_nodes<- same_elements(lgrph, "type", "nodes")
dist.nodes <- dist(as.matrix(df.same_nodes), method = "euclidean")
dist.edges <- dist(as.matrix(df.same_edges), method = "euclidean")
hc.nds <- hclust(dist.nodes, method = "ward.D")
hc.eds <- hclust(dist.edges, method = "ward.D") 
plot(hc.nds, main = "Common nodes", cex = .8)
plot(hc.eds, main = "Common edges", cex = .8)
```

Clustering of decorations on common nodes and  clustering on common edges can be directly compared to one another:

```{r hclust.compar, warning=FALSE, fig.align="center", fig.width=7}
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(dplyr))
par(mfrow=c(1, 2))
dend.nds <- as.dendrogram(hc.nds)
dend.eds <- as.dendrogram(hc.eds)
dendlist(dend.nds, dend.eds) %>%
  untangle(method = "step1side") %>% 
  tanglegram(columns_width = c(6, 1, 6),
             main_left = "Common nodes",
             main_right = "Common edges",
             lab.cex = 1.3,
             cex_main = 1.5,
             highlight_branches_lwd = F) 
```
  
In both clusterings, the Brozas stelae (decoration 3) and the Ibahernando stelae (decoration 5) are the ones having the most important proximities (ie, the least Euclidian distance). The 'Common edges' clustering is more accurate than the 'Common nodes' clustering because the former takes into account the common combination, or common permutation, of two nodes with the type of their edge, while the latter only takes into account the presence of common nodes.

```{r compare.c.edges, out.width = "700px", fig.align="center", warning=FALSE}
dec.to.compare <- c(3, 5)
g.compar <- list_compar(lgrph, nd.var = "type")
nds_compar <- plot_compar(listg = g.compar, 
                          dec2comp = dec.to.compare,
                          focus = "nodes",
                          nd.size = c(0.5, 1.7),
                          dir = dataDir,
                          img.format = "png")
eds_compar <- plot_compar(listg = g.compar, 
                          dec2comp = dec.to.compare,
                          focus = "edges",
                          nd.size = c(0.5, 1.7),
                          dir = dataDir,
                          img.format = "png")
knitr::include_graphics(nds_compar)
knitr::include_graphics(eds_compar) 
```
  
Here, the comparisons have been made on the basis of different types (ie, graphical units, GUs) with the variable `type` (`nd.var = "type"`). However, if a new column is added to the node dataframe or shapefile, the study can also incorporate the technique by which the GU was created (`nd.var = "technique"`) or any other categorical variable. For example, two GUs on the Brozas steale were made with incisions (`g_inc`): the fibula (*fibula de codo tipo Huelva*, ca 1050-950 BC) and the comb. 

<center>

![GIS view. The Brozas stelae (decoration 1) with the different engraving techniques: `g_piq` (peckings) and `g_inc` (incisions) ](`r brozas.path`){width=400px}

</center>

# Next developments

Graph theory allows one to construct tree structures for categorical variables (eg. the different types of GUs). These structures allow generalization processes (upward to the parent level) and specification processes (downward to the child level). For example, a `sword` and a `spear` belong both to the `weapons` group (sub-group `offensive weapons`), a `shield` belongs to the `weapons` group (sub-group `defensive weapons`), etc.: 

```{r hierac, warning=FALSE, fig.width=7, fig.height=4, asp=0}
par(mar=c(0,0,0,0))
g <- graph_from_literal(objects-+weapons,
                        objects-+personnal_item,
                        weapons-+offensive_weapons,
                        weapons-+defensive_weapons,
                        offensive_weapons-+spear,
                        offensive_weapons-+sword,
                        defensive_weapons-+shield,
                        defensive_weapons-+helmet,
                        personnal_item-+miror,
                        personnal_item-+comb)
layout <- layout.reingold.tilford(g)
plot(g,
     layout = layout,
     vertex.color = "white",
     vertex.frame.color = "white",
     vertex.size = 20,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     edge.arrow.size = 0.5
     )
```
  
Such a formalism can be used to weight the differences between nodes, to conduct analysis with different levels of precision or to overcome issues of idiosynchratic typologies.

Using the *diachronic* edge notation, tree structure can also be used to construct relative chronology diagrams, similar to a Harris matrix. For example, with the `data.tree` package and the Ibahernando stelae (decoration 5):

```{r overlap, warning=FALSE, out.width="100%"}
library(data.tree)
lgrph <- list_dec(imgs, nodes, edges)
edges.iba <- igraph::as_data_frame(lgrph[[5]], what="edges")
overlap.nodes <- unique(edges.iba[edges.iba$type == ">", "from"])
contemp.nodes <- unique(unlist(edges.iba[edges.iba$type == "=", c("from", "to")]))
df.stratig <- data.frame(over = rep(overlap.nodes, length(contemp.nodes)),
                         under = contemp.nodes)
df.stratig$pathString <- paste(lgrph[[5]]$decor,
                               df.stratig$over, 
                               df.stratig$under, 
                               sep = "/")
superpo <- as.Node(df.stratig)
print(superpo)
```

  
# References


---
title: | 
  | **The R *iconr* package** 
  | GIS interface entry
author: "Thomas Huet, Jose M Pozo, Craig Alexander"
bibliography: references.bib
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
---

```{r, echo=FALSE}
# print(getwd())
# print(file.path("../logo", "iconr_logo.png"))
htmltools::img(src = knitr::image_uri(file.path("img", "iconr_logo.png")), 
# htmltools::img(src = knitr::image_uri(file.path(R.home("doc"), "html", "logo.jpg")), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;')
```

---

<style>
.figure {
   margin-top: 0px;
   margin-bottom: 40px;
}
table {
    margin-top: 0px;
    margin-bottom: 24px;
}
</style>


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(kableExtra)
library(dplyr)

solanas.vor.path <- "img/all_process.gif"
```
  

Any iconographic contents can be modeled with a geometric graph where nodes, also called **graphical units (GUs)**, linked together with edges and then analyzed with the **Graph Theory** and spatial analysis at the support scale. This modeling is particularly interesting for compositions coming from Paleolithic, Neolithic and Protohistoric times [@Alexander08; @HuetAlexander15; @Huet18a].


The ***iconr*** R package grounds concepts and provides normalized tools to manage iconographic contents. To record large series of iconographic contents, the [GIS interface](#gis) appears to be the most appropriate one for users. This tutorial explains how to construct the graph _before_ using the package ***iconr***, offering also tipping points to facilitate the recording process. The chapter '[Dataset](https://zoometh.github.io/iconr/docs/#Dataset)' of the  ***iconr*** documentation can complete the present tutorial.

<center>

![GIS view. Detail of a Late Bronze Age stelae (Solana de Cabañas, Cáceres, Spain). From left to right: stelae photograph (photograph credits: Museo Arqueológico Nacional, Madrid); Georeferencing of the steale drawing over its photograph (dawing credits: Diaz-Guardamino 2010); Binarization and polygonization/vectorization of the graphical content of the steale (now `POLYGONS`); Calcul of their centroid; Calcul of their Voronoi cells; Binary topological relationships (*birel*) for each pairwise of Voronoi cells: the ones that share a border (*touches* = TRUE) will share a link (ie, edge) between their centroids (ie, nodes); Identification of the different types (nodes' column `type`)](`r solanas.vor.path`)

</center>


## A GIS ? Yes but why ? {#gis}

The *tenet* of the ***iconr*** is to always keep the user connected with the iconographic content -- his primary data source -- and emphasise the significance of the spatial dimension for any graphical content. Geographical Information Systems (GIS) offer multiple tools and options to facilitate the data entry.  Use of GIS offers a graphic interface and ensures the correctness of spatial relationships between GUs. It forms a permanent interface between the image of he decoration and the database. Obviously, the main GIS facility is the presence of **attributes tables** which allow to record, filter and sort GUs on many information: types, techniques, orientations, lengths, etc. 

The other most important GIS facilities for the recording process are:

* the snapping tools
* the measurement line
* the georeferencing tools (see [Absolute scale](#scale.abs))

From far, our software preference goes to [QGIS](https://www.qgis.org/fr/site/), because it is open source, offers a large rank of database connections facilities (with PostGRES/GIS for example), has a large user community, but also because the source file is a XML (.qgs, .qgz) structure that can be parsed, modified, copied and moved with scripting languages like R and Python

## Always start with an **image**

The image will be the reference space of the graph. So, before anything, start by opening the image decoration into the GIS. In this tutorial we will take the example of the South Western Iberian Abela stela dated to the Middle Bronze Age. The original drawing can be download [here](https://github.com/zoometh/iconr/blob/master/docs/tutorial/img/Abela.jpg) [@DiazGuardamino10]


<center>
  
  
![Drag and drop the decoration image before anything](img/openImageGIS.gif){width=70%}
</center>
  
  
Beside the ***iconr*** embedded training dataset, others datasets are available. To facilitate the starting, ready to use QGIS projects (.qgz) are also provided for each decoration:

```{r, echo = F}
gh.root <- "https://downgit.github.io/#/home?url=https://github.com/zoometh/iconr/tree/master/docs/datasets/"
url.df <- c(paste0(gh.root, "stele%20bas%20aragon"),
            paste0(gh.root, "Valcamonica"))
url.ref <- c("https://doi.org/10.3989/gladius.2013.0002",
             "https://scholar.google.fr/scholar?hl=fr&as_sdt=0%2C5&q=Le+Pietre+degli+Dei%2C+menhir+e+stele+dell%27Eta+del+Rame+in+Valcamonica+e+Valtellina&btnG=")
df.datasets <- data.frame("download" = c("stele bas aragon", "Valcamonica"),
                          "name" = c("Estelas Ibericas con Lanzas", "Ossimo"),
                          "nb.decor" = c("5","3"),
                          "sort.description" = c("Late Iron Age stelae with ranges of spears, writings, etc.",
                                                 "Chalcolithic rock-art of Ossimo, Valcamonica"),
                          "drawings.references" = c("Vargas 2013","Fedele 1994"),
                          stringsAsFactors = F)
# row.names(df.datasets) <- df.datasets$family
df.datasets %>% 
  mutate(download = cell_spec(download, "html", link = url.df)) %>%
  mutate(drawings.references = cell_spec(drawings.references, "html", link = url.ref)) %>%
  kable("html", escape = FALSE, row.names = F) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F)
```

These datasets can be downloaded by clicking on their names

```{r, echo=FALSE, eval=FALSE}
datasets.path <- paste0(getwd(),"/docs/datasets")

ui <- fluidPage(
  selectInput("A", "family", choices = df.datasets$download),
  selectInput("B", "decoration", choices = NULL),
  plotOutput("plot")
)

server <- function(input, output, session) {
  observe({
    choices_B <- read.csv(paste0(getwd(),"/docs/datasets/", input$A, "/imgs.csv"), sep=";")
    choices_B_dec <- paste0(choices_B$idf, ".", choices_B$site, "_", choices_B$decor)
    updateSelectInput(session, "B", choices = choices_B_dec)
  })
  output$plot <- renderPlot({
    # input$newplot
    # Add a little noise to the cars data
    # cars2 <- cars + rnorm(nrow(cars))
    # plot(cars2)
    a.idf <- as.numeric(unlist(strsplit(input$B,"\\."))[1])
    dataDir <- paste0(datasets.path, "/", input$A)
    df <- read.csv(paste0(dataDir, "/imgs.csv"), sep=";")
    # Decoration to be plotted
    site <- df[a.idf, "site"]
    decor <- df[a.idf, "decor"]
    # Read nodes, edges, and decorations
    nds.df <- read_nds(site, decor, dataDir, format = "shp")
    eds.df <- read_eds(site, decor, dataDir, format = "shp")
    imgs <- read.table(paste0(dataDir, "/imgs.tsv"),
                       sep="\t", stringsAsFactors = FALSE, header = T)
    
    # Save the plot of nodes and edges with node variable "type" as labels
    # in png image format and return the image file name.
    plot_dec_grph(nds.df, eds.df, imgs,
                  site, decor,
                  dir = dataDir,
                  nd.var = "type")
  })
}


shinyApp(ui, server)
```

### Scales

Use of GIS easy the scaling process. The creation of a spatialized graph permit to combine network analysis with spatial analysis of the graphical content

#### Relative scale {#scale.rel}

The image extent is measured in pixels with a top-left corner origin (0,0). The coordinates system is irrelevant: image, nodes and edges are measured on the pixel grid


<center>

![Coordinates system in pixels](img/coordsQGIS.gif){width=70%}
</center>

#### Absolute scale {#scale.abs}

To retrieve to true scale of the decoration, one can create a scale bar and apply a simple rule of three to convert pixels into centimeters, or meters. For example, if the scale belongs to another drawing, you can import it and 'georeferenced' it on the original drawing with the [*Freehand raster georeferencer* plugin](https://gvellut.github.io/FreehandRasterGeoreferencer/), and then create the scale bar

<center>

![Importing information on the scale and creation of a scale bar](img/scaleAll.gif){width=70%}

</center>

Parallely, the dimensions of each GU can be measured with the QGIS **Measure Line** tool. At first, only the maximum length of the GU is important. It has also to be noted that if a Polygonization is done on the GUs, the maximum length -- between all other type of shape analysis indexes -- do not have to be calculated 

<center>

![Measure and record the dimensions of each GUs in pixels](img/scaleRuleOfThree.gif){width=70%}

</center>

To retrieve the real dimensions of each GUs, first get the pixel sizes with the Measure line tools, then apply a simple rule of three with:

1. the size of the scale in pixels (native QGIS function `$length`)
2. the real size of the scale in cm (here, 100 cm)
3. the size of each GUs in cm


## Nodes and edges attribute tables

[Nodes](#nodes) and [edges](#edges) are the graph elements. At first, we have to create attribute tables for each of them. For example the nodes shapefiles:


<center>
  
![Create the **nodes** shapefile with its attribute table](img/createNodes.gif){width=70%}
  
</center>

In the GIS, [add a node](#graph.enodes.add) for each GUs and [add an edge](#graph.edges.add) between two contiguous GUs. Always start with the nodes

[HELP](https://zoometh.github.io/iconr/docs/#Edge)
[HELP](https://zoometh.github.io/iconr/docs/#Node)
[HELP](https://zoometh.github.io/iconr/docs/#Table_of_decorations)


### Nodes {#nodes}

**Nodes** represent the basic information of the graphical content. For example, it would be easy to distinguish a decoration with aurochs (`type = auroch`) from a decoration with swords (`type = sword`). The former iconographical content should be probably related to the Late Paleolithic/Early Neolithic while the second one is more probably related to a period starting with the Bronze Age. Nodes are created as a shapefile of `POINTS`. The attribute table of the nodes has at least four (4) fields:

1. `site` (Text): name of the site
2. `decor` (Text): name of the decoration
3. `id` (Integer): node identifier
4. `type`(Text): one of the relevant characteristics of each node

The nodes are created near the centroids of each different graphical units (GUs).

<center>
  
![Adding **nodes** for each GU](img/addNodes.gif){width=70%}
  
</center>

In this example, beside the *main* nodes sword (`epee`), anciform (`anciforme`) and halberd (`hallebarde`), we can also observe that the sword is connected to a belt and the anciform is worn as necklace. We probably would like to register this piece of graphical content as we also probably would like to characterize the types of blades for the sword and halberd, point out the presence of rivets on the sword depiction, etc. To do so, see the [Attribute edges](https://zoometh.github.io/iconr/docs/#Edge_types) part of the ***iconr*** package documentation

### Edges {#edges}

**Edges** types (field `type`) give information on nodes relative locations and on the nature of these nodes (main node *vs* attribute node, overlapping *vs* overlapped node, etc.). Edges are created as a shapefile of `LINES`. Edges attribute table has at least five (5) fields:

1. `site` (Text): name of the site
2. `decor` (Text): name of the decoration
3. `a` (Integer): *starting* node
4. `b` (Integer) *ending* node
5. `type` (Text): values`=`, `+` or `>` [HELP](https://zoometh.github.io/iconr/docs/#Edge_types)

Theoretically, between two *main* nodes, edges exist when their Voronoi cells are contiguous. In practice, if you consider that two GUs are neighbors, you can create an edge between their two nodes: GIS snapping tool help !

<center>
  
![Adding edges between contiguous GU](img/addEdges.gif){width=70%}
  
</center>


### Summary

For the Abela decoration, we have created three (3) nodes (`1`,`2`,`3`) and two (2) *normal* edges (`1-=-2`,`1-=-2`). We named the nodes shapefile `nodes.shp` and the edges shapefile `edges.shp` because this is their default name in the ***incor*** package

<center>
  
![](img/inputSummary.gif){width=70%}
  
</center>

## Create the **table of decorations**

The **table of decorations** is a correspondance table which records joins between nodes and edges dataframes. The default name of the **table of decorations** is `imgs.tsv` in the ***incor*** package, a tabulate separated-values (but it also can be a`.csv`, comma separated-values)

<center>
  
![](img/createTableOfDecorations.gif){width=70%}
  
</center>

The **table of decorations** has four (4) mandatory fields:

1. `idf` (Integer): short name of the decoration
1. `site` (Text): name of the site
2. `decor` (Text): name of the decoration
4. `img` (Text): name of the drawing/photograph/...


## Start to work with the ***incor*** package

Now we have at least one decoration modeled with a geometric graph, we can start to use the ***incor*** package. Its latest development version with its documentation (vignette) can be downloaded from GitHub

```{r down,eval=FALSE, echo=TRUE}
devtools::install_github("zoometh/iconr", build_vignettes=TRUE)
```

And load the package

```{r echo=TRUE}
library(iconr)
```

To start using the package, you have first to locate your working directory. For example:

```{r echo=TRUE}
dataDir <- paste0(getwd(), "/extdata")
```

The you can start with the function `plot_dec_grph()` and specifying the extensions of the nodes and edges data (`shp`)

```{r echo=TRUE, out.width="60%", fig.width=8, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Abela stelae"}
# Decoration to be plotted
site <- "Abela"
decor <- "Abela"
# Read nodes, edges, and decorations
nds.df <- read_nds(site, decor, dataDir, format = "shp")
eds.df <- read_eds(site, decor, dataDir, format = "shp")
imgs <- read.table(paste0(dataDir, "/imgs.tsv"),
                   sep="\t", stringsAsFactors = FALSE, header = T)

# Save the plot of nodes and edges with node variable "type" as labels
# in png image format and return the image file name.
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor,
              dir = dataDir,
              nd.var = "type")
```

# References

---
title: | 
  | **The R *iconr* package** 
  | a tutorial
bibliography: references.bib
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
---

```{r, echo=FALSE}
# print(getwd())
# print(file.path("../logo", "iconr_logo.png"))
htmltools::img(src = knitr::image_uri(file.path("img", "iconr_logo.png")), 
# htmltools::img(src = knitr::image_uri(file.path(R.home("doc"), "html", "logo.jpg")), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;')
```

---

<style>
.figure {
   margin-top: 0px;
   margin-bottom: 40px;
}
table {
    margin-top: 0px;
    margin-bottom: 24px;
}
</style>


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(kableExtra)
library(dplyr)

solanas.vor.path <- "img/all_process.gif"
```
  

 the ***iconr*** embedded training dataset, others datasets are available. To facilitate the starting, ready to use QGIS projects (.qgz) are also provided for each decoration:
```{r, echo = F}

library(httr)
req <- GET("https://api.github.com/zoometh/iconr/docs/datasets/git/trees/master?recursive=1")
stop_for_status(req)
filelist <- unlist(lapply(content(req)$tree, "[", "path"), use.names = F)
grep("Matteo/literature/", filelist, value = TRUE, fixed = TRUE)

url.root <- "https://raw.github.com/zoometh/iconr/master/docs/datasets/"
doss1 <- paste0(url.root, "Valcamonica/")
file1 <- paste0(doss1, "edges.tsv")
download.file(url = URL, destfile=basename(URL))
gh.root <- "https://downgit.github.io/#/home?url=https://github.com/zoometh/iconr/tree/master/docs/datasets/"
url.df <- c(paste0(gh.root, "stele%20bas%20aragon"),
            paste0(gh.root, "Valcamonica"))
url.ref <- c("https://doi.org/10.3989/gladius.2013.0002",
             "https://scholar.google.fr/scholar?hl=fr&as_sdt=0%2C5&q=Le+Pietre+degli+Dei%2C+menhir+e+stele+dell%27Eta+del+Rame+in+Valcamonica+e+Valtellina&btnG=")
df.datasets <- data.frame("download" = c("stele bas aragon", "Valcamonica"),
                          "name" = c("Estelas Ibericas con Lanzas", "Ossimo"),
                          "nb.decor" = c("5","3"),
                          "sort.description" = c("Late Iron Age stelae with ranges of spears, writings, etc.",
                                                 "Chalcolithic rock-art of Ossimo, Valcamonica"),
                          "drawings.references" = c("Vargas 2013","Fedele 1994"),
                          stringsAsFactors = F)
# row.names(df.datasets) <- df.datasets$family
df.datasets %>% 
  mutate(download = cell_spec(download, "html", link = url.df)) %>%
  mutate(drawings.references = cell_spec(drawings.references, "html", link = url.ref)) %>%
  kable("html", escape = FALSE, row.names = F) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F)
```


## Start to work with the ***incor*** package

Now we have at least one decoration modeled with a geometric graph, we can start to use the ***incor*** package. Its latest development version with its documentation (vignette) can be downloaded from GitHub

```{r down,eval=FALSE, echo=TRUE}
devtools::install_github("zoometh/iconr", build_vignettes=TRUE)
```

And load the package

```{r echo=TRUE}
library(iconr)
```

To start using the package, you have first to locate your working directory. For example:

```{r echo=TRUE}
dataDir <- paste0(getwd(), "/extdata")
```

The you can start with the function `plot_dec_grph()` and specifying the extensions of the nodes and edges data (`shp`)

```{r echo=TRUE, out.width="60%", fig.width=8, fig.asp=750/666, fig.align="center", warning=FALSE, fig.cap="Abela stelae"}
# Decoration to be plotted
site <- "Abela"
decor <- "Abela"
# Read nodes, edges, and decorations
nds.df <- read_nds(site, decor, dataDir, format = "shp")
eds.df <- read_eds(site, decor, dataDir, format = "shp")
imgs <- read.table(paste0(dataDir, "/imgs.tsv"),
                   sep="\t", stringsAsFactors = FALSE, header = T)

# Save the plot of nodes and edges with node variable "type" as labels
# in png image format and return the image file name.
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor,
              dir = dataDir,
              nd.var = "type")
```

---
title: "***iconr*** tutorial"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


## Open Image with **QGIS**
  
Open the decoration in QGIS. The width of the image is equivalent to the x-axis, the height of the image is equivalent to the y-axis
  
<center>

![](img/coordsQGIS.gif)

</center>
  
The QGIS coordinates are:

* top-left: origins (0,0)

* top-right: (666,0)

* bottom-right: (666,750)

* bottom-left: (0,750)  
  

<center>

![](img/extentGIS.gif){width=66%}

</center>
  
# Open Image with ***magick*** R package

Read the decoration image with the R magick package and show image proprieties

```{r echo=TRUE, message=FALSE, warning=FALSE}
library(knitr)
library(kableExtra)
library(magick)
path_img <- "D:/Sites_10/Abela/Abela.jpg"
img <- image_read(path_img)
image_info(img)
```

Plot the image
  

```{r echo=TRUE, out.width="60%", fig.width=6, fig.asp=750/666}
plot(image_trim(img))
```

Add the coordinates of the image four corners:

* top-left

* top-right

* bottom-right

* bottom-left

```{r echo=TRUE, out.width="60%", fig.width=6, fig.asp=750/666}
# par(mar = c(0,0,0,0))
plot(img)
box(lwd = 2)
# min width & min height (origin)
text(0, 0, 
     paste0(0, ",", 0),
     cex = 1, adj = c(0, -0.2))
# min width & max height
text(0, image_info(img)$height,
     paste0(0, ",", image_info(img)$height),
     cex = 1, adj = c(0, 1.1))
# max width & min height
text(image_info(img)$width, 0,
     paste0(image_info(img)$width, ",", 0),
     cex = 1, adj = c(1, -0.2))
# max width & max height
text(image_info(img)$width, image_info(img)$height,
     paste0(image_info(img)$width, ",", image_info(img)$height),
     cex = 1, adj = c(1, 1.1))
```

\name{named_elements}
\alias{named_elements}
\title{Textual Notation of Graph Elements
}
\description{
  Create a textual notation for nodes or edges.
}
\usage{
  named_elements(grph,
                 focus = "edges",
                 nd.var = "type",
                 disamb.marker = "#")
}
\arguments{
  \item{grph}{
    A decoration graph
    (object of class \code{igraph}).
  }
  \item{focus}{
    Textual notation of edges (\code{focus = "edges"}) or nodes
    (\code{focus = "nodes"}). By default \code{focus = "edges"}.
  }
  \item{nd.var}{
    The attribute of the graph nodes containing the node variable (ie, field)
    for the textual annotation. By default \code{nd.var = "type"}.
  }
  \item{disamb.marker}{
    Marker used to disambiguated repeated elements.
    By default \code{disamb.marker = "#"}.
  }
}
\details{
  Edges of type \code{'='} (\emph{normal} edges) are \strong{undirected}, so
  that the order of their nodes is irrelevant and they are presented in
  alphabetical order.
  Conversely, edges of types \code{'+'} (\emph{attribute} edges) and \code{'>'}
  (\emph{diachronic} edges) are \strong{directed}, so that the given order of
  nodes is preserved.

  Repeated node or edge names are
  disambiguated by appending the symbol \code{disamb.marker} (\code{'#'} by
  default) at the end of the second appearance (suffix). Subsequent appearances
  are marked by additional \code{disamb.marker}s.
}
\value{
  A character vector of named nodes or edges.
}

\seealso{
  \code{\link[iconr]{list_compar}},
  \code{\link[iconr]{same_elements}}
}
\examples{
# Read data
imgs <- read.table(system.file("extdata", "imgs.tsv", package = "iconr"),
                   sep="\t", stringsAsFactors = FALSE)
nodes <- read.table(system.file("extdata", "nodes.tsv", package = "iconr"),
                    sep="\t", stringsAsFactors = FALSE)
edges <- read.table(system.file("extdata", "edges.tsv", package = "iconr"),
                    sep="\t", stringsAsFactors = FALSE)

# Generate list of graphs from the three data.frames
lgrph <- list_dec(imgs, nodes, edges)

# Textual notation of disambiguated edges
named_elements(lgrph[[2]], focus = "edges", nd.var="type")

# Textual notation of disambiguated nodes
named_elements(lgrph[[2]], focus = "nodes", nd.var="type")

}

\keyword{ ~kwd1 graph}% use one of  RShowDoc("KEYWORDS")
\name{read_nds}
\alias{read_nds}
\title{Read Nodes of a Decoration}
\description{
  Read nodes' information from a file including all nodes and extract nodes of one decoration.
  Accepted formats are tab separated values ('tsv'), semicolon separated values ('csv'), or
  shapefile ('shp').
}
\usage{
read_nds(site,
        decor,
        dir = getwd(),
        nodes = "nodes",
        format = "tsv")
}
\arguments{
  \item{site}{
      Name of the site
  }
  \item{decor}{
      Name of the decoration
  }
  \item{dir}{
      Path to the working folder, by default it is the working directory
  }
  \item{nodes}{
      Name of the nodes file (a dataframe or a shapefile)
  }
  \item{format}{
    File extension indicating a file format from 'tsv' (tab separated values),
    'csv' (semicolon separated values) or 'shp' (shapefile). For 'tsv' and 'csv'
    the files must include node coordinates (\code{nodes$x}, \code{nodes$y}).
  }
}
\value{
  Dataframe of graph nodes, including at least the columns
  "site", "decor", "id", "x", "y", with values for each node (row).
}

\examples{
# Set data folder
dataDir <- system.file("extdata", package = "iconr")

# Read dataframe of nodes
nds.df <- read_nds(site = "Cerro Muriano", decor = "Cerro Muriano 1",
                    dir = dataDir, format = "tsv")
nds.df
## Dataframe of nodes

# Read shapefile of nodes
nds.df <- read_nds(site = "Cerro Muriano", decor = "Cerro Muriano 1",
                    dir = dataDir, format = "shp")
nds.df
## Dataframe of nodes

}

\keyword{ ~kwd1 graphs}% use one of  RShowDoc("KEYWORDS")
\name{plot_compar}
\alias{plot_compar}
\title{Plot and Save Comparison Figures Between Pairs of Graphs}
\description{
  Given a list of pairwise graph comparisons, the function plots any given subset selected by graph name, displaying side-by-side pairs of graphs and highlighting common nodes or common edges with a choice of several graphical parameters.
}
\usage{
  plot_compar(listg, dec2comp = NULL, focus = "nodes",
              dir = getwd(),
              nd.color = c("orange", "red"), nd.size = c(0.5, 1),
              ed.color = c("orange", "red"), ed.width = c(1, 2),
              lbl.size = 0.5,
              dir.out = dir, out.file.name = NULL,
              img.format = NULL, res = 300)
}
\arguments{
  \item{listg}{
      A list of graph pairwise comparisons as returned by \code{\link[iconr]{list_compar}}.
  }
  \item{dec2comp}{
      A vector with the names of the graphs for which comparisons are to be plotted.
    The user can select to plot all pairwise combinations (by default), all combinations of a subset, or a single pair.
  }
  \item{focus}{
      Either \code{"nodes"} (default) or \code{"edges"}. It selects the type of comparison to be plotted, highlighting common nodes or common edges, respectively.
  }
  \item{dir}{
    Data folder including the decoration images.
    By default the working directory.
  }
  \item{nd.color, nd.size, ed.color, ed.width}{
    Graphical parameters for color and size/widths of nodes and edges.
    Each of them is a vector with two values for different and common nodes/edges, respectively.
    If only one value is provided, this unique value is taken for both different and common elements.
    Labels are displayed with the same color as common nodes.
    For \code{focus = "nodes"} all edges are plotted with the first value of \code{ed.color} and \code{ed.width}.
  }
  \item{lbl.size}{
      Graphical parameter for the size of the labels with the node names. The default is 0.5.
  }
  \item{dir.out}{
    Folder for the output image. By default, it coincides with the input \code{dir}.
}
  \item{out.file.name}{
    Name of the output image, including path from current directory and extension.
    By default the name is automatically generated including \code{site},
    \code{decor}, \code{nd.var}, and the extension from \code{img.format}.

    If set, \code{out.file.name} overrides \code{dir.out} and \code{img.format}.
}
  \item{img.format, res}{
      Format and resolution of the saved images. The handled formats are
      \code{"png"}, \code{"bmp"}, \code{"tiff"}/\code{"tif"},
      \code{"jpeg"}/\code{"jpg"}, and \code{"pdf"}.
      The default resolution is 300 (ppi). The resolution does not applies to format pdf.

      if \code{img.format=NULL} (default), the plot is sent to the active device.
  }
}
\details{
  To highlight common elements between a list of graphs, the user can focus on nodes \code{(focus = "nodes")} or edges \code{(focus = "edges")}. As stated in the function \code{\link[iconr]{list_compar}}, for a given comparison variable (eg. \code{nd.var="type"}) if there is multiple nodes/edges with the same value, it is considered to count for as many coincidences as the smaller multiplicity.

  \code{img.format=NULL} (plot to the active device) does not make sense for
  more than one comparison.
}
\value{
    Generates graph decoration images, for pairwise comparisons between two or more decorations, comparing graphs elements (nodes or edges).

   If \code{img.format=NULL}, the plot is sent to the active device and no value is returned.

   If \code{img.format=} \code{"png"} or \code{"bmp"} or \code{"tiff"}/\code{"tif"} or \code{"jpeg"}/\code{"jpg"} or \code{"pdf"}, the return value is a character vector with the dir/name of every saved image in the indicated format.

}
\seealso{
  \code{\link[iconr]{list_compar}}
  \code{\link[iconr]{plot_dec_grph}}
}

\examples{
# Read data
imgs <- read.table(system.file("extdata", "imgs.tsv", package = "iconr"),
                   sep="\t",stringsAsFactors = FALSE)
nodes <- read.table(system.file("extdata", "nodes.tsv", package = "iconr"),
                    sep="\t",stringsAsFactors = FALSE)
edges <- read.table(system.file("extdata", "edges.tsv", package = "iconr"),
                    sep="\t",stringsAsFactors = FALSE)

# Generate list of graphs from the three dataframes
lgrph <- list_dec(imgs, nodes, edges)

# Generate all pairwise comparisons of the graphs with respect to nodes "type"
g.compar <- list_compar(lgrph, nd.var="type")

# Generate the image showing the comparison on common nodes of graphs
# '1' and '4', save it in png format, and return its path.
dataDir <- system.file("extdata", package = "iconr")
outDir <- tempdir()
plot_compar(g.compar, c(1,4), focus = "nodes",
            dir = dataDir,
            dir.out = outDir,
            img.format = "png")

# Generate the image showing the comparison on common edges of all pairwise
# combinations of graphs '1','3', and '4', save them in pdf format, and return
# their path.
# Plot nodes involved in non-common edges in orange and
# nodes involved in common edges and the corresponding labels in brown.
plot_compar(g.compar, c(1, 3, 4), focus = "edges",
            dir = dataDir,
            nd.color = c("orange", "brown"),
            dir.out = outDir,
            img.format = "pdf")

# Save the png image showing the comparison on common nodes of graphs
# '1' and '4'.
# Then read and plot the image.
img.filename <- plot_compar(g.compar, c(1, 4), focus = "nodes",
                            dir = dataDir,
                            dir.out = outDir,
                            img.format = "png")
plot(magick::image_read(img.filename))

# Plot directly on the active device (default) the comparison on common nodes
# of graphs '1' and '4'.
plot_compar(g.compar, c(1, 4), focus = "nodes",
            dir = dataDir)
}

\keyword{ ~kwd1 graph}% use one of  RShowDoc("KEYWORDS")
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conv_pg_to_shp.R
\name{conv_pg_to_shp}
\alias{conv_pg_to_shp}
\title{Convert Pg geometries to SHP geometries}
\usage{
conv_pg_to_shp(
  dataDir = tempdir(),
  Pg.param = NA,
  sqll.obj = NA,
  sqll.ug.pts = NA,
  sqll.ug.lines = NA,
  sqll.ug.polyg = NA,
  exp.edges = FALSE
)
}
\arguments{
\item{dataDir}{working directory: sites' folders (copied data) and 'out' folder (outputs)}

\item{Pg.param}{list of arguments to connect the PostgreSQL database.
Like: list(driver, name_of_db, host, port, user, password). By default NA}

\item{sqll.obj}{SQL on objects to get images of decorations. By default NA}

\item{sqll.ug.pts}{SQL on GUs with geometries of type POINTS to get shapes. By default NA}

\item{sqll.ug.lines}{SQL on GUs with geometries of type LINES to get shapes. By default NA}

\item{sqll.ug.polyg}{SQL on GUs with geometries of type POLYGONS to get shapes. By default NA}

\item{exp.edges}{Export also edges as shapefiles. By default: FALSE}
}
\value{
Decoration's images and shapefiles
}
\description{
Convert the graphical units (GUs) geometries stored in PostgreSQL into shapefiles (SHP) geometries.
}
\examples{
Pg.param. <- list("PostgreSQL",
                  "postgres",
                  "localhost",
                  5432,
                  "postgres",
                  "postgres")
dataPath <- "D:/decorations/"
dataDir <- paste0(dataPath, "analyse")

### SQL on Pg
## Objects
## sqll.obj. <-  "SELECT
## site, numero, img
## FROM objets
## WHERE (objets.site LIKE 'Ain Ghazal' AND objets.numero LIKE 'stat_2') OR
## (objets.site LIKE 'Ain Ghazal' AND objets.numero LIKE 'stat_5_xd') OR
## (objets.site LIKE 'Qarassa' and objets.numero LIKE 'figurine__wx') OR
## (objets.site LIKE 'Jericho' and objets.numero LIKE 'tete_afe')"

## Polygons
## sqll.ug.polyg. <-  "SELECT
## objets.site,
## objets.numero,
## table_noeuds.id,
## table_noeuds.type,
## table_noeuds.technologie as technlg,
## table_noeuds.incomplet as incmplt,
## ST_AsText(table_noeuds.geom_shape) as wkt FROM objets LEFT JOIN table_noeuds
## ON table_noeuds.site = objets.site AND table_noeuds.decor = objets.numero
## WHERE (objets.site like 'Ain Ghazal' AND objets.numero like 'stat_2')
##  OR (objets.site like 'Ain Ghazal' AND objets.numero like 'stat_5_xd')
##  OR (objets.site like 'Qarassa' AND objets.numero like 'figurine__wx')
##  OR (objets.site like 'Jericho' AND objets.numero like 'tete_afe')"

## Run function
## conv_pg_to_shp(dataDir = dataDir,
##                Pg.param = Pg.param.,
##                sqll.obj = sqll.obj.,
##                sqll.ug.polyg = sqll.ug.polyg.)

}
\name{contemp_nds}
\alias{contemp_nds}

\title{Select Contemporaneous Nodes}
\description{
  Find the connected component, or subgraph, of contemporaneous nodes (connected by \emph{normal} and \emph{attribute} edges) given a selected node and remove the other components
}
\usage{
  contemp_nds(nds.df, eds.df, selected.nd)
}
\arguments{
  \item{nds.df}{
      Dataframe of the nodes as the one obtained by the function
      \code{\link[iconr]{read_nds}}.
  }
  \item{eds.df}{
      Dataframe of the edges as the one obtained by the function
      \code{\link[iconr]{read_eds}}.
  }
  \item{selected.nd}{
      The node of the decoration graph for which to extract the connected component. It can be either the node order (numeric) or the node name/id (character).
  }
}
\value{
  A named list of two dataframes: \code{list(nodes, edges)}, collecting
  the contemporaneous nodes and edges, respectivelly.
}

\examples{
# Set data folder
dataDir <- system.file("extdata", package = "iconr")

# Read a decoration
nds.df <- read_nds(site = "Ibahernando",
                   decor = "Ibahernando",
                   dir = dataDir)
eds.df <- read_eds(site = "Ibahernando",
                   decor = "Ibahernando",
                   dir = dataDir)

# Extract the subgraph contemporaneous to the node 2
l_dec_df <- contemp_nds(nds.df, eds.df, selected.nd = 2)

## It returns a list of two dataframes, one for nodes and one for edges:
l_dec_df
}
\keyword{ ~kwd1 graphs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morph_nds_group.R
\name{morph_nds_group}
\alias{morph_nds_group}
\title{Morphometrics classification between GUs}
\usage{
morph_nds_group(
  nodes,
  focus = c("clust", "kmeans"),
  gu.types = "all",
  nb.centers = 1,
  out.dir = "_out",
  out.data = c("mbrshp", "plot")
)
}
\arguments{
\item{nodes}{Dataframe of nodes}

\item{focus}{Type of grouping, hierachical clustering ("clust") or Kmeans ("kmeans").
By default, c("clust", "kmeans")}

\item{gu.types}{Classes of nodes that will be clustered, a vector of characters or a character.
By default "all"}

\item{nb.centers}{Number of clusters, uniquely for Kmeans. By default 1 (unique cluster)}

\item{out.dir}{Name of output folder}

\item{out.data}{Type of data returned.
If "mbrshp" return a dataframe of nodes with their clustering and image path.
If "plot" return a "kmeans" or create a plot. By default c("mbrshp", "plot")}
}
\value{
Depending on the focus, return hierachical clustering ("clust") or Kmeans ("kmeans") plots with their complete path
}
\description{
Morphometrics classification (groups) between different graphical units (GUs).
 Read JPG files from each different folder. Useful after comparisons (see, m'orph_nds_compar' function)
}
\examples{
morph_nds_group(nodes)

## [1] "* read 'oeil' typo"
## Extracting 10.jpg outlines...
## [ 1 / 10 ]  Ain Ghazal.stat_2.1.jpg
## ...
## [ 10 / 10 ]  Qarassa.figurine__wx.14.jpg

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conv_shp_to_wkt.R
\name{conv_shp_to_wkt}
\alias{conv_shp_to_wkt}
\title{Convert SHP to WKT}
\usage{
conv_shp_to_wkt(dataDir, complete.only = TRUE, out.dir = "_out")
}
\arguments{
\item{dataDir}{Path of the folder storing folders of all decorations.
Each of these folders as a site name (eg, Ain Ghazal) and contains at least
one shapefile (.shp and .dbf and .shx) and one image (.jpg or .tif or .png, etc.).
The shapefile is named conventionally with the name of the site,
 a dot,
 the name of the decoration,
 "nd_pl" for nodes POLYGONS
  (eg, Ain Ghazal.stat_2_nd_pl.shp)
The image is named conventionally with the name of the site,
 a dot,
 and the name of the decoration (eg, Ain Ghazal.stat_2.tif)}

\item{complete.only}{Boolean, if TRUE discard incomplete (incmplt == 1), by default = TRUE}

\item{out.dir}{Path of the output folder. By default "_out/" in the "dataDir" folder}
}
\value{
Create the 'nodes.csv' file into the out folder, return the complete path of the 'nodes.csv' file
}
\description{
Read GUs (ie, nodes) shapefiles (POINTS, LINES, POLYGONS) from the different sites' folders,
store nodes in the 'nodes.csv' dataframe with geometries as Well-Known Text (WKT) representations,
return path of the node dataframe
}
\examples{
nd.df.path <- conv_shp_to_wkt(dataDir = dataDir)
head(read.csv2(nd.df.path), 1)

##         site  decor id type technlg incmplt geometry
## 1 Ain Ghazal stat_2  1 oeil       -       0 POLYGON ((266.9252 -167.608,...

}
\name{side_plot}
\alias{side_plot}
\title{Plot Two Figures Side-by-Side Identifying Common Elements}
\description{
  Plot two decoration graphs side-by-side identifying common nodes and common edges. This function is called by the function \code{\link[iconr]{plot_compar}}.
}
\usage{
side_plot(grph, dir, nd.var, focus = "nodes",
          nd.color = c("orange", "red"),
          nd.size = c(0.5, 1),
          ed.color = c("orange", "red"),
          ed.width = c(1, 2),
          lbl.size = 0.5)
}
\arguments{
        \item{grph}{
    List of two or more 'igraph' graphs created with the \code{\link[iconr]{list_compar}} function.
}
          \item{dir}{
    Working directory which contains the imgs, nodes, edges dataframes and the decoration images.
}
            \item{nd.var}{
    Field of nodes on which the comparison will be done.
}

            \item{focus}{
    Focus on nodes or on edges, by default \code{focus = "nodes"}.
}
            \item{nd.color, nd.size, ed.color, ed.width}{
    Graphical parameters for the nodes and edges. The \strong{different} nodes/edges will be displayed with the first values of the vectors (eg, "orange") while the \strong{common} nodes/edges will be displayed with the second values of the vectors (eg, "red").
}
            \item{lbl.size}{
%%     ~~Describe \code{x} here~~
    Size of the labels
}
}
\value{No return value, group images side-by-side}

\seealso{
  \code{\link[iconr]{plot_compar}}
}
\keyword{ ~kwd1 graphs}% use one of  RShowDoc("KEYWORDS")
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morph_nds_compar.R
\name{morph_nds_compar}
\alias{morph_nds_compar}
\title{Morphometrics comparisons between GUs}
\usage{
morph_nds_compar(
  nodes,
  focus = c("panel", "stack", "PCA"),
  nb.h = 15,
  out.dir = "_out"
)
}
\arguments{
\item{nodes}{Dataframe of nodes}

\item{focus}{Type of analysis: 'panel', 'stack' or 'PCA'. By default c("panel", "stack", "PCA")}

\item{nb.h}{number of Fourier harmonics, uniquely for PCA. By default = 15}

\item{out.dir}{path of the output folder. By default "_out/" in the "dataDir" folder}
}
\value{
Depending on the focus, return 'panel', 'stack' or 'PCA' plots with their complete path
}
\description{
Momocs package morphometrics comparisons between different graphical units (GUs).
 Read JPG files from each different folder. Useful before grouping (see, 'morph_nds_group' function)
 to determine the number of clusters
}
\examples{
morph_nds_compar(nodes)

## [1] "* read 'oeil' type of UGs"
## Extracting 10.jpg outlines...
## [ 1 / 10 ]  Ain Ghazal.stat_2.1.jpg
## ...
## [ 10 / 10 ]  Qarassa.figurine__wx.14.jpg

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conv_wkt_to_jpg.R
\name{conv_wkt_to_jpg}
\alias{conv_wkt_to_jpg}
\title{Convert WKT to JPG}
\usage{
conv_wkt_to_jpg(nodes, out.dir = "_out")
}
\arguments{
\item{nodes}{nodes dataframe coming from the 'conv_shp_to_wkt.R' function}

\item{out.dir}{path of the output folder. By default "_out/" in the "dataDir" folder}
}
\value{
JPGs files exported into as many folders as different GUs' types (eg., 'bouche', oeil', 'visage', etc.)
The names of the JPGs is the concatenate of their site name, dot, decoration name, dot, GU identifier (eg, "Ain Ghazal.stat_2.1.jpg")
}
\description{
Convert Well-Known Text geometries of graphical units (GUs) to JPG files
 in order to perform contour analysis with the Momocs package
}
\examples{
conv_wkt_to_jpg(nodes = nodes)

## Saving 4.33 x 3.94 in image
## Saving 4.33 x 3.94 in image
## ...

}
\name{read_eds}
\alias{read_eds}
\title{Read Edges of a Decoration}
\description{
  Read edges' information from a file including all edges and extract edges of one decoration.
  Accepted formats are tab separated values ('tsv'), semicolon separated values ('csv'), or
  shapefile ('shp').
}
\usage{
read_eds(site,
         decor,
         dir = getwd(),
         edges = "edges",
         nodes = "nodes",
         format = "tsv")
}
\arguments{
        \item{site}{
    Name of the site.
}
          \item{decor}{
    Name of the decoration.
}
        \item{dir}{
    Path to the working folder, by default it is the working directory.
}
            \item{edges}{
    Name of the edges file (a dataframe or a shapefile).
}
            \item{nodes}{
    Name of the nodes file (a dataframe or a shapefile).
}
              \item{format}{
    File extension indicating a file format from 'tsv' (tab separated values),
    'csv' (semicolon separated values) or 'shp' (shapefile). For 'tsv' and 'csv'
    the coordinates of the edges will be calculated from the same decoration's
    node dataframe.
}
}
\details{
Subset the dataframe of edges depending on 'site' and 'decor'.
}
\value{
  Dataframe of graph edges, including at least the columns "site", "decor",
  "a", "b", "xa", "ya", "xb", "yb", with values for each edge (row).
}

\examples{
# Set data folder
dataDir <- system.file("extdata", package = "iconr")

# Read .tsv file
eds.df <- read_eds(site = "Cerro Muriano", decor = "Cerro Muriano 1",
                   dir = dataDir, edges = "edges", format = "tsv")
eds.df
## Dataframe of edges

# Read shapefile
eds.df <- read_eds(site = "Cerro Muriano", decor = "Cerro Muriano 1",
                   dir = dataDir, edges = "edges", format = "shp")
eds.df
## Dataframe of edges

}

\keyword{ ~kwd1 graph}% use one of  RShowDoc("KEYWORDS")
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morph_resume.R
\name{morph_resume}
\alias{morph_resume}
\title{Resume geometry classes (POINTS, LINES, POLYGON) of GUs for each object}
\usage{
morph_resume(
  dataDir,
  nodes = NA,
  out.dir = "_out",
  imgs.format = c(".jpg", ".png", ".gif", ".tiff", "tif")
)
}
\arguments{
\item{dataDir}{Path of the folder storing folders of all decorations}

\item{nodes}{Nodes dataframe coming from the 'conv_shp_to_wkt.R' function}

\item{out.dir}{Path of the output folder. By default "_out/" in the "dataDir" folder}

\item{imgs.format}{Accepted picture formats (see 'magick' package, 'image_read' function)}
}
\value{
Create a contact sheet of decoration with information on graphical units (GUs) geometries:
eg., number of Polygons by types, number of Lines by types, etc.
}
\description{
Create a contact sheet of decoration with information on graphical units (GUs) geometries:
eg., number of Polygons by types, number of Lines by types, etc.
}
\examples{
> morph_resume(dataDir = dataDir,
+              nodes = nodes)
[1] "Ain Ghazal"
[1] "Jericho"
[1] "Qarassa"

}
\name{same_elements}
\alias{same_elements}
\title{Number of Equal Elements Between Each Decoration Pair}
\description{
  Create the (symmetric) dataframe with the count of \strong{common nodes} or \strong{common edges} (see \code{\link[iconr]{list_compar}} for comparison criteria) for each pair of decorations (graphs) from a list. The diagonal of the symmetric dataframe is filled with counts of nodes/edges for each decoration.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
same_elements(lgrph, nd.var = "type",
               focus = "nodes")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
      \item{lgrph}{
%%     ~~Describe \code{x} here~~
    A list of any number of graphs to be pairwise compared. The list can
    be typically obtained with the function \code{\link[iconr]{list_dec}}
}
      \item{nd.var}{
%%     ~~Describe \code{x} here~~
    An attribute of the graph vertices containing the node variable (ie, field)
    on which the comparison will be done. By default \code{nd.var = "type"}.
}
      \item{focus}{
%%     ~~Describe \code{x} here~~
      Either \code{"nodes"} (default) or \code{"edges"} to select the
      type of elements to be compared for the count.
}
}
\value{
  A symmetric matrix with the counts of the pairwise coincidences of nodes or edges. The matrix has as row and column names the names of the corresponding graphs in the input list.
}
\seealso{
  \code{\link[iconr]{list_dec}},
  \code{\link[iconr]{list_compar}},
  \code{\link[iconr]{plot_compar}}
}

\examples{
# read imgs, nodes and edges dataframes
imgs <- read.table(system.file("extdata", "imgs.tsv", package = "iconr"),
                   sep="\t",stringsAsFactors = FALSE)
nodes <- read.table(system.file("extdata", "nodes.tsv", package = "iconr"),
                    sep="\t",stringsAsFactors = FALSE)
edges <- read.table(system.file("extdata", "edges.tsv", package = "iconr"),
                    sep="\t",stringsAsFactors = FALSE)
lgrph <- list_dec(imgs,nodes,edges)

# Counting same nodes
df.same_nodes <- same_elements(lgrph, nd.var = "type",
                               focus = "nodes")
df.same_nodes
## a symmetric matrix of nodes comparisons

# same edges
df.same_edges <- same_elements(lgrph, nd.var = "type",
                               focus = "edges")
df.same_edges
## a symmetric matrix of edges comparisons

}

\keyword{ ~kwd1 graph}% use one of  RShowDoc("KEYWORDS")
\name{list_dec}
\alias{list_dec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create Decoration's Graphs and Store them in a List}
\description{
  Create undirected graphs for each decoration from \code{nodes}, \code{edges} and \code{imgs} dataframes and store the graphs in a list.
  The join between these dataframes is done on the two fields \code{site} and \code{decor}.
  Graph names refer to \code{imgs$idf}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  list_dec(imgs,
           nodes,
           edges)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
      \item{imgs}{
%%     ~~Describe \code{x} here~~
    Dataframe of decorations
}
        \item{nodes}{
%%     ~~Describe \code{x} here~~
    Dataframe of nodes
}
          \item{edges}{
%%     ~~Describe \code{x} here~~
    Dataframe of edges
}
}
\value{
  A list of \code{igraph} graphs.
}
\seealso{
  \code{\link[igraph]{graph_from_data_frame}}
}

\examples{
# Read imgs, nodes and edges dataframes
imgs <- read.table(system.file("extdata", "imgs.csv", package = "iconr"),
                   sep=";", stringsAsFactors = FALSE)
nodes <- read.table(system.file("extdata", "nodes.csv", package = "iconr"),
                    sep=";", stringsAsFactors = FALSE)
edges <- read.table(system.file("extdata", "edges.csv", package = "iconr"),
                    sep=";", stringsAsFactors = FALSE)
# Create the list of graphs
lgrph <- list_dec(imgs, nodes, edges)

# Get the first graph
g <- lgrph[[1]]
g

# Graph name
g$name

# Graph label
g$lbl

# Graph number of nodes
igraph::gorder(g)

# Graph number of edges
igraph::gsize(g)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 graph}% use one of  RShowDoc("KEYWORDS")
\name{labels_shadow}
\alias{labels_shadow}
\title{Plot Labels with Contrasting Shadow}
\description{
  Plot labels (text) with a contrasting buffer to make them more visible
  when located on a similar color background.
  This function is the \code{shadowtext()} function developed by Greg Snow.
  Called by plot functions: \code{\link[iconr]{plot_dec_grph}}, \code{\link[iconr]{plot_compar}}
}
\usage{
labels_shadow(x, y = NULL, labels,
              col = "black", bg = "white",
              theta = seq(0, 2 * pi, length.out = 50),
              r = 0.1,
              cex = 1, ...)
}
\arguments{
        \item{x, y}{
    Numeric vector of coordinates where the labels should be plotted.
    Alternatively, a single argument \code{x} can be provided
    with the same syntax as in \code{\link[grDevices]{xy.coords}}.
}
          \item{labels}{
    Set of labels provided as a character vector.
}
        \item{col, bg}{
    Graphical parameters for the label color and background (buffer) color.
}
        \item{theta}{
    Angles for generating the buffer with possible anisotropy along one
    direction (default is isotropic)
    and controlling buffer smoothness (angular resolution).
}
        \item{r}{
    Thickness of the buffer relative to the size of the used font, by default 0.1.
}
        \item{cex}{
    Size of the label, by default 1.
}
        \item{...}{
    Further graphical parameter accepted by \code{\link[graphics]{text}}, such as
    \code{pos}, \code{offset}, or \code{family}.
}
}
\value{No return value. It creates a contrasting buffer to make labels more visible.}

\references{
  https://rdrr.io/cran/TeachingDemos/man/shadowtext.html
}
\keyword{ ~kwd1 graphs}
\name{plot_dec_grph}
\alias{plot_dec_grph}
\title{Plot a Graph on a Decoration}
\description{
  Plot with nodes only, edges only, or both (geometric graph) over a decoration image.
}
\usage{
plot_dec_grph(nodes = NULL,
              edges = NULL,
              imgs,
              site,
              decor,
              dir = getwd(),
              nd.var = 'id',
              nd.color = 'orange',
              nd.size = 0.5,
              lbl.color = 'black',
              lbl.size = 0.5,
              ed.color = c("orange", "blue"),
              ed.lwd = 1,
              dir.out = dir,
              out.file.name = NULL,
              img.format = NULL,
              res = 300)
}
\arguments{
  \item{nodes}{
    Dataframe of nodes
}
  \item{edges}{
    Dataframe of edges
}
  \item{imgs}{
    Dataframe of decorations
}
  \item{site}{
    Name of the site
}
  \item{decor}{
    Name of the decoration
}
  \item{dir}{
    Data folder including the decoration images.
    By default the working directory.
}
  \item{nd.var}{
    Field name in the nodes data frame to be displayed as node labels.
    By default the identifier \code{nodes$id}.
}
  \item{nd.color,
        nd.size,
        lbl.color,
        lbl.size,
        ed.color,
        ed.lwd}{
    Graphical parameters for color and size/widths of nodes, edges, and labels.
    \code{ed.color} is a vector with two values (the second value is used for diachronic
    edges).
}
  \item{dir.out}{
    Folder for the output image. By default, it coincides with the input \code{dir}.
}
  \item{out.file.name}{
    Name of the output image, including path from current directory and extension.
    By default the name is automatically generated including \code{site},
    \code{decor}, \code{nd.var}, and the extension from \code{img.format}.

    If set, \code{out.file.name} overrides \code{dir.out} and \code{img.format}.
}
  \item{img.format, res}{
      Format and resolution of the saved images. The handled formats are
      \code{"png"}, \code{"bmp"}, \code{"tiff"}/\code{"tif"},
      \code{"jpeg"}/\code{"jpg"}, and \code{"pdf"}.
      The default resolution is 300 (ppi). The resolution does not applies to format pdf.

      if \code{img.format=NULL} (default), the plot is sent to the active device.
  }
}
\details{
  Plot \strong{nodes only} (if \code{edges = NULL}), \strong{edges only} (if \code{nodes = NULL}), or both (graph) over a decoration image.
}
\value{
   Generates graph decoration images with nodes, edges, or both, overlapping the decoration image.

   If \code{img.format=NULL}, the plot is sent to the active device and no value is returned.

   If \code{img.format=} \code{"png"} or \code{"bmp"} or \code{"tiff"}/\code{"tif"} or \code{"jpeg"}/\code{"jpg"} or \code{"pdf"}, the return value is a character vector with the dir/name of the saved image in the indicated format.
}

\examples{
## Set data folder
dataDir <- system.file("extdata", package = "iconr")
## Decoration to be plotted
site <- "Brozas"
decor <- "Brozas"
## Read nodes, edges, and decorations
nds.df <- read_nds(site, decor, dataDir)
eds.df <- read_eds(site, decor, dataDir)
imgs <- read.table(paste0(dataDir, "/imgs.tsv"),
                   sep="\t", stringsAsFactors = FALSE)

## Plot 'Brozas' nodes and edges on the active device
## with node variable "type" as labels
plot_dec_grph(nds.df, eds.df, imgs,
              site, decor,
              dir = dataDir,
              lbl.size = 0.4,
              nd.var = "type")

## Save only edges of 'Brozas' with bigger widths and in image format jpg.
outDir <- tempdir()
img.filename <- plot_dec_grph(nodes = NULL, eds.df, imgs,
                              site, decor,
                              dir = dataDir,
                              ed.lwd = 2,
                              dir.out = outDir,
                              img.format = "jpg")
## Then read and plot the image.
a.dec <- magick::image_read(img.filename)

## Inspect the output image
magick::image_info(a.dec)

## Plot the output image
plot(a.dec)

}

\keyword{ ~kwd1 graph}% use one of  RShowDoc("KEYWORDS")
\name{list_compar}
\alias{list_compar}
\alias{nds_compar}
\alias{eds_compar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graph Pairwise Comparison on Common Elements}
\description{
    \code{nds_compar} identifies \strong{common nodes} in a pair of graphs.

    \code{eds_compar} identifies \strong{common edges} in a pair of graphs.

    Given a list of graphs, \code{list_compar} extract all combinations of graph pairs and compare them on common elements (nodes and edges).
}
\usage{
nds_compar(grphs, nd.var = "type")

eds_compar(grphs, nd.var = "type")

list_compar(lgrph, nd.var = "type",
            verbose = FALSE)
}
\arguments{
  \item{grphs}{
      A list of two graphs (pair of graphs) to be compared.
  }
  \item{lgrph}{
    A list of any number of graphs to be pairwise compared. The list can be typically obtained with the function \code{\link[iconr]{list_dec}}
  }
    \item{nd.var}{
    An attribute of the graph nodes containing the node variable (ie, field) on which the comparison will be done. By default \code{nd.var = "type"}.
  }
            \item{verbose}{
%%     ~~Describe \code{x} here~~
    Logical. If TRUE, the names of each graph pair combination are listed on the screen. By default \code{verbose = FALSE}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    \code{list_compar()} calls the functions: \code{nds_compar()} and \code{eds_compar()} which return respectively the \strong{common nodes} and the \strong{common edges} of a graph pairwise.

    \strong{Nodes} are common when they have the same value for a given variable, for example \code{horse}, \code{sword}, etc., for the variable \code{type} (\code{nd.var = "type"}).

     \strong{Edges} are common when they have the same value for \emph{starting} and \emph{ending} nodes (\code{horse}, \code{sword}, etc.) and the same type of edge (\code{'='}, \code{'+'}, etc.).
    For example, \code{a -=- b} in graph 1 is equal to \code{a -=- b} in graph 2, but not equal to \code{a -+- b}. Edges of type \code{=} (\emph{normal} edges) are undirected, so that \code{a -=- b} is equal to \code{b -=- a}. But edges of types \code{+} (\emph{attribute} edges) or \code{>} (\emph{diachronic} edges) are directed, so: \code{a ->- b} is not equal to \code{b ->- a}.

If any of the graphs has multiple nodes/edges with the same value, it is considered to count for as many coincidences as the smaller multiplicity. For instance, if there are 2 nodes with value \code{epee} in graph 1, and 3 nodes with value \code{epee} in graph 2, their number of common nodes is \code{min(2, 3) = 2}.
}
\value{
  \code{nds_compar()} returns the input pair of graphs, each complemented with a new node attribute named \code{comm} with value 1 for common nodes and 0 for non-common nodes.

  \code{eds_compar()} returns the input pair of graphs, each complemented with a new edge attribute named \code{comm} with value 1 for common edges and 0 for non-common edges.

  \code{list_compar()} returns a list of all combinations of graph pairs. For each pair, both graphs are complemented with the node attribute (\code{comm}) identifying common nodes and the edge attribute (\code{comm}) identifying common edges. Each pair is also complemented with an attribute named \code{nd.var} recording the compared node variable.
}
\seealso{
  \code{\link[iconr]{list_dec}},
  \code{\link[iconr]{plot_compar}},
  \code{\link[iconr]{same_elements}}
}

\examples{
# Read data
imgs <- read.table(system.file("extdata", "imgs.tsv", package = "iconr"),
                   sep="\t",stringsAsFactors = FALSE)
nodes <- read.table(system.file("extdata", "nodes.tsv", package = "iconr"),
                    sep="\t",stringsAsFactors = FALSE)
edges <- read.table(system.file("extdata", "edges.tsv", package = "iconr"),
                    sep="\t",stringsAsFactors = FALSE)
# Generate list of graphs from the three data.frames
lgrph <- list_dec(imgs, nodes, edges)

# Generate list of all graph comparisons depending on the node "type" variable
g.compar <- list_compar(lgrph, nd.var = "type")

length(g.compar)
## Ten pairwise comparisons

# Inspect the second pairwise comparison of the list
g.compar[[2]]
## The two compared graphs with the name of the comparison variable

# Inspecting nodes:
igraph::as_data_frame(g.compar[[2]][[1]], "vertices")
## Vertices from the first decoration graph

igraph::as_data_frame(g.compar[[2]][[2]], "vertices")
## Vertices from the second decoration graph

# Inspecting edges:
igraph::as_data_frame(g.compar[[2]][[1]])
## Edges of the first decoration graph

igraph::as_data_frame(g.compar[[2]][[2]])
## Edges of the second decoration graph

}

\keyword{ ~kwd1 graph}% use one of  RShowDoc("KEYWORDS")
