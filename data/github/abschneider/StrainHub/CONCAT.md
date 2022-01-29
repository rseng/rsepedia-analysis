<p align="center">
  <img src="https://github.com/abschneider/StrainHub/raw/master/img/mainlogo.png" alt="StrainHub logo" width="400"/>
</p>

<details>	
<summary><strong><em>Table of Contents</em></strong></summary>

* [About StrainHub](#about-strainhub)
* [Getting Started](#getting-started)
  - [Example Files](#example-files)
* [Network Quality](#network-quality)
* [How To Cite](#how-to-cite)
* [References](#references)
* [Contact](#contact)

</details>

-------------------

## About StrainHub

StrainHub is designed as a web-based software to generate pathogen transmission networks and associated metrics from a combination of a phylogenetic tree and a metadata associated. The software enables the construction of the network by either  mapping the metadata onto the tree and performing a parsimony ancestry reconstruction step to create links between the associated metadata or parsing a BEAST phylogeography maximum clade credibility tree.

<p align="center">
  <img src="https://github.com/abschneider/StrainHub/raw/master/host_network_example.png" alt="Sample Host Transmission Network" width="350"/>
</p>

## Getting Started

In order to construct the disease transmission network the user is presented with two options, run the entire pipeline or skip straight to the visualization part. 

To run the entire pipeline, two files will need to be generated.

1) A phylogenetic tree formatted in Newick tree format generated through your preferred phylogenetic search method (e.g., BEAST, TNT, RAxML, IQTree). Alternatively, you can use an aligned multiple sequence file in FASTA format.

2) A metadata associated file formatted as a comma separated value (CSV) file that includes headers, has the Accession number as the first column and the metadata associated values (e.g., host, country, risk group) 

If you would like to visualize the distribution of your pathogen in a map, a third optional geodata file need to be generated with all locations and geographic coordinates.
 
_Taxa ID with missing data should be excluded prior to the analysis. Order does not matter on both files, the metadata is sorted automatically before mapped into the tree based on the header of the first column of the CSV file (Accession)._

_A phylogenetic tree in Newick format is required to perform the Parsimony Ancestral Reconstruction ("Parsimony" method). If you elect to utilize a multiple sequence alignment instead of your tree, you need to run the alignment through the "Create Neighbor-Joining Tree" method to generate a tree which will be fed into the Parsimony Ancestral Reconstruction pipeline._

A template for building the metadata file can be downloaded [here](https://github.com/abschneider/StrainHub/blob/master/data/example_metadata.csv). Do not change the header of the Accession column as it is necessary to identify and reorder the metadata according to the tree file.

### Example files:
- [Phylogenetic Tree](https://github.com/abschneider/StrainHub/raw/master/data/etc/example_tree.phy)
- [Metadata File](https://github.com/abschneider/StrainHub/raw/master/data/etc/example_metadata.csv) 
- [Geodata File](https://github.com/abschneider/StrainHub/raw/master/data/parsimonious/hepc/2k1b_country_coordinates.txt)

To skip the ancestry reconstruction step, the user will have to run BEAST phylogeography, and then one file will need to be generated.

1) A maximum clade credibility tree (MCC) generated through Phylogeographic Diffusion in [Discrete](http://beast.community/workshop_discrete_diffusion) or [Continuous](http://beast.community/workshop_continuous_diffusion) Space in BEAST in NEXUS format.

Our app will parse the MCC tree and prompt the user to select which trait he wants to evaluate. Once the user has selected the trait, the user will have to select the trait probability threshold the user wants to use to consider the relationship between two nodes valid (we like 0.9 as a threshold, but this will depend on the user and his knowledge about his dataset).

_The threshold set by the user is utilized to filter out any pair of nodes within the phylogenetic tree that may have lower probability, thus being ignored within the transmission network step._ 

### Example file:

- [MCC Tree](https://raw.githubusercontent.com/abschneider/StrainHub/master/data/batRABV.mcc.tree)

Note: To download example files, right click on the link and save file.

## Metrics

The nodes within the network are scaled based on the user selected metric. The current metrics available on StrainHub are:

- __Betweenness Centrality__
  
    Betweenness centrality of node i is the sum of all the ratios between all the paths g from node j to node k that go through node i over all possible paths from node j to node k.
  
    Epidemiological meaning: How important a location/host is as the shortest path / intermediary connecting other locations/hosts within the transmission network.

- __Closeness Centrality__
  
    Closeness centrality is the inverse of the average distance between pk and the other nodes.
  
    Epidemiological meaning: How important a location/host is given its distance within the transmission network to other locations/hosts.

- __Degree Centrality__
  
    Degree centrality of node k is the sum of all links of other nodes to node k over the total number of nodes – 1.
  
    Epidemiological meaning: How important a location/host is within the transmission network given the number of times that a disease emerges from or to that point (indegree or outdegree).
  
    Degree Centrality is also subdivided in Indegree and Outdegree Centrality, which calculate the links directed only in or out of the note, respectively.

- __Source Hub Ratio (SHR)__
  
    The Source Hub Ratio of a node i is the sum of all shifts from location “i” to other locations over the sum of all shifts from and to location ”i”.
  
    Epidemiological meaning: How important a node is within the network as source of the disease, ignoring centrality within network.
  
        1 = all shifts originated on node i (source)
      
        0.5 = same amount of shifts from and to node i (hub)
      
        0 = all shifts occurred to node i (dead end)

_Although only one metric can be selected for scaling the nodes within the network, all the metrics are calculated and can be visualized as well as downloaded on the tab "Metrics"._

## Network Quality

The overall quality of the transmission network generated on StrainHub rests on the quality of the user input data.

## How To Cite

Please cite the reference below if you are using StrainHub.

Text Citation:

```
Adriano de Bernardi Schneider, Colby T Ford, Reilly Hostager, John Williams, Michael Cioce, Ümit V Çatalyürek, Joel O Wertheim, Daniel Janies, StrainHub: A phylogenetic tool to construct pathogen transmission networks, Bioinformatics, btz646, https://doi.org/10.1093/bioinformatics/btz646
```

BibTex Citation:

```
@article{10.1093/bioinformatics/btz646,
    author = {de Bernardi Schneider, Adriano and Ford, Colby T and Hostager, Reilly and Williams, John and Cioce, Michael and Çatalyürek, Ümit V and Wertheim, Joel O and Janies, Daniel},
    title = "{StrainHub: A phylogenetic tool to construct pathogen transmission networks}",
    journal = {Bioinformatics},
    year = {2019},
    month = {08},
    abstract = "{In exploring the epidemiology of infectious diseases, networks have been used to reconstruct contacts among individuals and/or populations. Summarizing networks using pathogen metadata (e.g., host species and place of isolation) and a phylogenetic tree is a nascent, alternative approach. In this paper, we introduce a tool for reconstructing transmission networks in arbitrary space from phylogenetic information and metadata. Our goals are to provide a means of deriving new insights and infection control strategies based on the dynamics of the pathogen lineages derived from networks and centrality metrics. We created a web-based application, called StrainHub, in which a user can input a phylogenetic tree based on genetic or other data along with characters derived from metadata using their preferred tree search method. StrainHub generates a transmission network based on character state changes in metadata, such as place or source of isolation, mapped on the phylogenetic tree. The user has the option to calculate centrality metrics on the nodes including betweenness, closeness, degree, and a new metric, the source/hub ratio. The outputs include the network with values for metrics on its nodes and the tree with characters reconstructed. All of these results can be exported for further analysis.strainhub.io and https://github.com/abschneider/StrainHub}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btz646},
    url = {https://doi.org/10.1093/bioinformatics/btz646},
    eprint = {http://oup.prod.sis.lan/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btz646/29171171/btz646.pdf},
}
```

Additionally, if you use the Map feature of StrainHub please also cite the reference below:

Text Citation:

```
Schneider, A.D.B., Ochsenreiter, R., Hostager, R., Hofacker, I.L., Janies, D. and Wolfinger, M.T., 2019. Updated Phylogeny of Chikungunya Virus Suggests Lineage-Specific RNA Architecture. Viruses, 11(9), p.798.
```

BibTex Citation:

```
@article{schneider2019updated,
  title={Updated Phylogeny of Chikungunya Virus Suggests Lineage-Specific RNA Architecture},
  author={de Bernardi Schneider, Adriano and Ochsenreiter, Roman and Hostager, Reilly and Hofacker, Ivo L and Janies, Daniel and Wolfinger, Michael T},
  journal={Viruses},
  volume={11},
  number={9},
  pages={798},
  year={2019},
  publisher={Multidisciplinary Digital Publishing Institute}
}
```

## Authors

<h5 align = "left">Adriano de Bernardi Schneider, Ph.D.<br>Colby T. Ford, Ph.D.<br>Reilly Hostager<br>John Williams<br> Michael Cioce<br>Ümit V. Çatalyürek, Ph.D.<br>Joel O. Wertheim, Ph.D.<br>Daniel Janies, Ph.D.</h5>

## References

[Adriano de Bernardi Schneider, _Arboviruses: The Hidden Path of an Imminent Threat_, University of North Carolina at Charlotte, 2018, ProQuest (10979260). ISBN: 9780438668379](https://pqdtopen.proquest.com/pubnum/10979260.html)

## Contact

Questions? Comments? Contact us [here](mailto:adebernardischneider@ucsd.edu).
<p align="left">
  <img src="https://github.com/abschneider/StrainHub/blob/beta/img/sidelogo.png" alt="StrainHub logo" width="350"/>
</p>

[![Version: 1.1.2](https://img.shields.io/badge/version-1.1.2-green.svg)](UPDATES.md)
[![Version: 1.1.2](https://img.shields.io/github/issues/abschneider/StrainHub.svg)](https://github.com/abschneider/StrainHub/issues)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

StrainHub is designed as a web-based software to generate disease transmission networks and associated metrics from a combination of a phylogenetic tree and associated metadata. The software maps the metadata onto the tree and performs a parsimony ancestry reconstruction step to create links between the associated metadata and enable the construction of the network. Users have the option to build a tree utilizing their method of preference outside StrainHub or build a tree utilizing a FASTA file within StrainHub with the Neighbor-Joining algorithm. Alternatively, the user can skip the StrainHub ancestry reconstruction step by generating a maximum clade credibility tree (MCC) through BEAST phylogeography or input a previously generated list of edges in order to build the transmission network. Additionally, the user can input a file with geographic cooordinates associated with the character of interest and have the network plotted into a map.

<p align="center">
  <img src="https://github.com/abschneider/StrainHub/blob/master/host_network_example.png" alt="Sample Host Transmission Network" width="350"/>
</p>

## StrainHub Features

StrainHub currently offers four methods to enable the visualization of transmission networks, each with different file format requirements:

**Parsimony**

You'll have the option of running a parsimony reconstruction on your phylogenetic tree, and for that you will need 2 files to get started:

1) A phylogenetic tree formatted in Newick tree format generated through your preferred phylogenetic search method (e.g., BEAST, TNT, RAxML, IQTree).

2) A metadata associated file formatted as a comma separated value (CSV) file that includes headers, has the Accession number as the first column and the metadata associated values (e.g., host, country, risk group)

**BEAST Phylogeography**

You'll run phylogeography using BEAST and for that you will only need 1 file to get started:

1) A maximum clade credibility tree (MCC) generated through Phylogeographic Diffusion in [Discrete](http://beast.community/workshop_discrete_diffusion) or [Continuous](http://beast.community/workshop_continuous_diffusion) Space in BEAST or Phylogeography in BEAST2 using [Mascot](https://taming-the-beast.org/tutorials/Mascot-Tutorial/).

**Quick Tree**

You'll create a tree directly from your alignment using a NJ algorithm within StrainHub and run parsimony reconstruction on the given tree. You will need 2 files and to know your outgroup sequence:

1) An alignment file in FASTA format generated through your preferred alignment method (e.g., MAFFT, Geneious, Clustal, MUSCLE).

2)  A metadata associated file formatted as a comma separated value (CSV) file that includes headers, has the Accession number as the first column and the metadata associated values (e.g., host, country, risk group)

**List of Edges**

You'll create a network directly from a list of edges previously generated in another software, you will need only 1 file:

1) List of edges in a CSV formatted file.

**Geographic Coordinates**

Additionally, if you want to plot your transmission network on a map you will need 1 additional file that does not vary between the methods:

1) Geographic coordinates for each character of interest in a CSV formatted file.

For more information, click [here](ABOUT.md).

## Use StrainHub Online

[![StrainHub.io](https://img.shields.io/badge/Try%20Out%20StrainHub%20Online-StrainHub.io-blue.svg?logo=r&style=for-the-badge&labelColor=2C3E50&color=3498DB)](https://strainhub.io)

## Run StrainHub Locally

StrainHub can be run on your local machine. This is the better option for scenarios where the input data is large or if security is a concern.

### Run StrainHub in Docker

The StrainHub Docker image is hosted on DockerHub: hub.docker.com/r/cford38/strainhub. Simply run the following from your command prompt:

```bash
docker pull cford38/strainhub:latest
docker run --name strainhub --rm -p 3838:3838 cford38/strainhub:latest
```

Then, In your browser, navigate to `localhost:3838`.


To build StrainHub's Docker image locally (which is unnecessary unless you want to make edits to the Docker image), run the following:
```bash
docker build -t strainhub .
docker run --name strainhub --rm -p 3838:3838 strainhub
```


### Run StrainHub Locally from RStudio and Shiny

To run, download the repository and open the `strainhub/app.R` file and run the following script.

```r
library(shiny)
runApp()
```

If you need to install all the required packages to run StrainHub, you can do so by running the `strainhub/install_packages.R` script.

## Authors

<h4 align = "left">Adriano de Bernardi Schneider, Ph.D.<br>Colby T. Ford, Ph.D.<br>Reilly Hostager<br>John Williams<br>Michael Cioce<br>Ümit V. Çatalyürek, Ph.D.<br>Joel O. Wertheim, Ph.D.<br>Daniel Janies, Ph.D.</h4>

## How To Cite


Text Citation:

```
Adriano de Bernardi Schneider, Colby T Ford, Reilly Hostager, John Williams, Michael Cioce, Ümit V Çatalyürek, Joel O Wertheim, Daniel Janies, StrainHub: A phylogenetic tool to construct pathogen transmission networks, Bioinformatics, btz646, https://doi.org/10.1093/bioinformatics/btz646
```

BibTex Citation:

```
@article{10.1093/bioinformatics/btz646,
    author = {de Bernardi Schneider, Adriano and Ford, Colby T and Hostager, Reilly and Williams, John and Cioce, Michael and Çatalyürek, Ümit V and Wertheim, Joel O and Janies, Daniel},
    title = "{StrainHub: A phylogenetic tool to construct pathogen transmission networks}",
    journal = {Bioinformatics},
    year = {2019},
    month = {08},
    abstract = "{In exploring the epidemiology of infectious diseases, networks have been used to reconstruct contacts among individuals and/or populations. Summarizing networks using pathogen metadata (e.g., host species and place of isolation) and a phylogenetic tree is a nascent, alternative approach. In this paper, we introduce a tool for reconstructing transmission networks in arbitrary space from phylogenetic information and metadata. Our goals are to provide a means of deriving new insights and infection control strategies based on the dynamics of the pathogen lineages derived from networks and centrality metrics. We created a web-based application, called StrainHub, in which a user can input a phylogenetic tree based on genetic or other data along with characters derived from metadata using their preferred tree search method. StrainHub generates a transmission network based on character state changes in metadata, such as place or source of isolation, mapped on the phylogenetic tree. The user has the option to calculate centrality metrics on the nodes including betweenness, closeness, degree, and a new metric, the source/hub ratio. The outputs include the network with values for metrics on its nodes and the tree with characters reconstructed. All of these results can be exported for further analysis.strainhub.io and https://github.com/abschneider/StrainHub}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btz646},
    url = {https://doi.org/10.1093/bioinformatics/btz646},
    eprint = {http://oup.prod.sis.lan/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btz646/29171171/btz646.pdf},
}
```
------------------------

Additionally, if you use the Map feature of StrainHub please also cite the reference below:

Text Citation:

```
Schneider, A.D.B., Ochsenreiter, R., Hostager, R., Hofacker, I.L., Janies, D. and Wolfinger, M.T., 2019. Updated Phylogeny of Chikungunya Virus Suggests Lineage-Specific RNA Architecture. Viruses, 11(9), p.798.
```

BibTex Citation:

```
@article{schneider2019updated,
  title={Updated Phylogeny of Chikungunya Virus Suggests Lineage-Specific RNA Architecture},
  author={de Bernardi Schneider, Adriano and Ochsenreiter, Roman and Hostager, Reilly and Hofacker, Ivo L and Janies, Daniel and Wolfinger, Michael T},
  journal={Viruses},
  volume={11},
  number={9},
  pages={798},
  year={2019},
  publisher={Multidisciplinary Digital Publishing Institute}
}
```
# StrainHub

## Release Updates

### v1.1.2
- Added the ability to export and import StrainHub networks as RDS files.

### v1.1.1
- Updated/removed some figure export abilities and prepared app for Docker containerization.

### v1.1.0
- Replaced underlying StrainHub functionality to reference R package version of StrainHub. (Thanks, Adesolji Ademiluyi!)


### v1.0.11
_Second public release!_
- Added Logo to Homepage
- Added Tip Label to Metadata Validation

### v1.0.10
- Added Homepage
- Added Posterior Probability Threshold for BEAST

### v1.0.9
- Added Google Analytics integration

### v1.0.8
- Added geodata editor for BEAST trees
- Hide metadata editor and geodata editor if files haven't been uploaded
- Add bootstrap value slider for NJ tree generation
- Add mapping options (basemap options, arrow styles, etc.)
- Enhance download capabilities for map, network, and metrics table
- Add error message is user forgets to specify NJ tree root

### v1.0.7
- Changed mapping to be Swoopy and Leaflet-based

### v1.0.6
- Bug fixes
- Add better error messaging
- Add globe map for BEAST and NJ inputs.
- Fix NJ tree generation

### v1.0.5
- Add geodata editor
- Increase flexibility of location column name for geodata
- Change centrality metric selection to dropdown (from radio buttons)
- Add better error messaging
- Change network title font to Lato

### v1.0.4
- Fix NJ Tree Creation
- Add globe visual for Parsimonious tree input

### v1.0.3
- Faster performance (as files are only read in once)
- Add metadata Editor

### v1.0.2
- Ability to export networks as HTML or PNG
- Ability to turn off directionality (arrows) on edges
>>>>>>> beta

### v1.0.0

_First public release!_

- Tree Files Accepted:
  - Parsimony
  - BEAST Phylogeography
- Visualization Capabilities:
  - Network Plot
  - Plotly-based Dendrogram
  - Metrics Output

### Coming Soon:

- Neighbor-joining tree generation
- Improved tree visualization
- Improved plot output
- ...and more!
# *strainhub* - Phylogenetic Network Visualization in R

<h3 align = "right">Colby T. Ford, Ph.D. and Adriano De Bernardi Schneider, Ph.D.</h3>

<img align="right" src="https://raw.githubusercontent.com/colbyford/StrainHub/master/img/strainhub_hex_color.png" alt="StrainHub icon" width="200">

StrainHub works to generate a pathogen transmission network graph utilizing genomic data and calculate importance of network based on centrality metrics.

## How to Install
```r
remotes::install_github("colbyford/strainhub", subdir="pkg")
```

Note: If you're running R 4.x, you may have to install some dependencies manually as they are not yet available via CRAN for the newest version of R. Simply find the troublesome package on GitHub and run `remotes::install_github("<ACCOUNT>/<REPO-NAME>")` in your R console.

## Examples
```r
## Load in StrainHub
library(strainhub)

## Read in tree, metadata, and geodata
treedata <- ape::read.tree("../../data/parsimonious/chikv/chikv_westernafrica.phy")
metadata <- readr::read_csv("../../data/parsimonious/chikv/chikv_westernafrica_metadata.csv", col_names = TRUE)
geodata <- readr::read_csv("../../data/parsimonious/chikv/chikv_geo.csv", col_names = TRUE)


## Check to See Which States are available by which to generate the network
listStates(treedata,
           metadata,
           treeType = "parsimonious")


## Make the Transmission Network
graph <- makeTransNet(treedata,
                      metadata,
                      columnSelection = "Country",
                      centralityMetric = 6,
                      treeType = "parsimonious")

print(graph)


## Make Leaflet/Swoopy map
make_map(graph,
         geodata = geodata,
         columnSelection = "Country")

```

## Customizing your Tranmission Networks
StrainHub highly depends on the [vizNetwork](https://github.com/datastorm-open/visNetwork) package. For a full list of customizations options, see the vizNetwork documentation here: https://datastorm-open.github.io/visNetwork/

<img src="https://raw.githubusercontent.com/colbyford/StrainHub/master/img/chikv_example.png" alt="Example custom Chikungunya virus network">

```r
## Using the graph object from before...
library(visNetwork)

nodes <- graph$x$nodes %>%
  mutate(shape = "dot",
         color = "grey",
         font.size = 20)

## Define new edge colors
colors <- c("blue", "green", "purple")

edges <- graph$x$edges %>%
  mutate(arrows = "to",
         smooth = TRUE,
         color = colors,
         width = ifelse(value == 1, 1, 4),
         value = NULL)

## Add a legend
lnodes <- data.frame(label = c("Size:\nSource Hub Ratio"),
                     shape = c("dot"),
                     color = c("grey"))
                     
ledges <- data.frame(color = c("grey", "green", "blue", "purple"),
                     label = c("Size:\nTransitions",
                               "Senegal to\nCote de Ivoire",
                               "Cote de Ivoire to\nSenegal",
                               "Senegal to\nNigeria"),
                     arrows = c("to", "to", "to", "to"))

## Generate the custom network using `visNetwork`
visNetwork(nodes, edges) %>%
  visOptions(nodesIdSelection = list(enabled = TRUE)) %>%
  visLegend(addEdges = ledges,
            addNodes = lnodes,
            useGroups = FALSE)
```
