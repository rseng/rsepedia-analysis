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
