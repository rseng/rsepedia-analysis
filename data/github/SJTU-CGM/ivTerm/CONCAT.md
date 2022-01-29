
# ivTerm

ivTerm is an interactive and comprehensive graphic system for functional analysis of meta-omics data.


## 1. Introduction

&emsp;&emsp;ivTerm is an R-shiny package with a user-friendly graphical interface that enables users to inspect functional annotations, compare results across multiple experiments, create customized charts. It utilizes shiny to create an interactive user interface (UI) helping researchers further explore and interpret gene functional analysis results. The combination of term data and gene data provides a full view that helps researchers gain actionable insights. The high interactivity allows a quick examination of complex data. The user-friendly interface allows users to customize plot multiple types of interactive charts and edit details of the figures.

## 2. Installation

### 2.1 Installing R/RStudio

&emsp;&emsp;ivTerm is a package in the R software environment, which can be freely downloaded as follows:

* Install [R](https://www.r-project.org/)
* Install [RStudio](https://www.rstudio.com/)

### 2.2 Installing ivTerm

&emsp;&emsp;Check or install required packages.

```
packages <- c("shiny", "shinyjs", "ggplot2", "ggiraph", "ggnetwork", "igraph", "DT", "RCurl", "XML", "colourpicker", "esquisse")
palette_packages <- c("ggsci", "randomcoloR", "scales", "viridis", "wesanderson")
lapply(c(packages, palette_packages), function(x) {
	if(!require(x, character.only = TRUE)) {
		install.packages(x, dependencies = TRUE)
	}})
```

&emsp;&emsp;Install ivTerm from github.

```
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)
install_github("SJTU-CGM/ivTerm", build_vignettes = TRUE)
```


## 3. Quick Start

&emsp;&emsp;Load the library 
```
library(ivTerm)
```

&emsp;&emsp;Run ivTerm using graphic interface. 
```
termTOgene()
```

## 4. Graphic visualization interface

&emsp;&emsp;The user interface contains three tabs: "Upload File", "Data Selection", and "Visualization". 


### 4.1 Upload file

&emsp;&emsp;To begin the analysis, you need to upload two files (comma-separated (.csv) or tab-separated (.txt) format). 

<center>
<figure>
<img src='vignettes/figure/fig1.png' width='50%'>

Fig 1. Graphic user interface for input data uploading

&nbsp;
</center>
</figure>

&emsp;&emsp;The first file (term data) should contain at least term ids (the column named 'term_id') and the genes annotated to each term (the column named 'gene'). Genes are connected by a separator. The separator can be a semicolon, comma, slash, or space. The column named 'group' is optional and used to compare multiple datasets or groups. The second file (gene data) should contain gene labels matching those used in term data and some details of genes. The structure of the data is shown below. 

&nbsp;

<center>The first example of term data</center>

|term_id|gene|term_name|enrich_cohort|P_value|mean_abun_DAN|mean_abun_CHN|...|
|---|---|---|---|---|---|---|---|
|K01196|168739,1061735,1125063,1429351,2396955,2681450, 3786103,4284614,5013911,5472907,675251...|glycogen debranching enzyme [EC:2.4.1.25 3.2.1.33]|DAN|0.0015|1.86e-06|2.25e-06|...|
|K00696|348550,541207,1432586,1489349,1759637,1925571, 2003961,2057616,2064314,2500326,2500973...|sucrose-phosphate synthase [EC:2.4.1.14]|DAN|2.82e-06|4.85e-06|1.75e-06|...|
|K00692|94503,142982,158676,185965,205094,259894,261947, 265221,285956,304808,319331,446946,45...|levansucrase [EC:2.4.1.10]|CHN|1.02e-06|3.78e-07|1.61e-06|...|
|...|...|...|...|...|...|...|...|

<center>The first example of gene data</center>

|gene|Gene_Length|Gene_Completeness_Status|Cohort_Origin|Phylum|Genus|...|
|---|---|---|---|---|---|---|
|1579397|1221|Complete|EUR|unknown|unknown|...|
|1579796|1221|Complete|EUR|unknown|unknown|...|
|1579826|1221|Complete|EUR|Bacteroidetes|Bacteroides|...|
|1580132|1221|Complete|EUR|Firmicutes|Ruminococcus|...|
|1580154|1221|Complete|EUR|unknown|unknown|...|
|...|...|...|...|...|...|...|


<center>The second example of term data</center>

|group|term_id|gene|term_name|p.val|...|
|---|---|---|---|---|---|
|coral_up|GO:0032741|coral_TRINITY_DN34635_C2_G1 coral_TRINITY_DN55005_C5_G1 coral_TRINITY_DN46011_C0_G2 coral_TRINITY_DN46011_C0_G1 coral_TRINITY_DN52251_C1_G1 coral_TRINITY_DN35049_C2_G1|positive regulation of interleukin-18 production|1.1e-06|...|
|zooxanthella_down|GO:0009063|zooxanthella_TRINITY_DN46713_C1_G2 zooxanthella_TRINITY_DN46713_C1_G3 zooxanthella_TRINITY_DN42914_C2_G1 zooxanthella_TRINITY_DN48316_C2_G1 zooxanthella_TRINITY_DN50127_C3_G1 zooxanthella_TRINITY_DN47858_C1_G7 zooxanthella_TRINITY_DN45516_C1_G3|cellular amino acid catabolic process|5.73e-06|...|
|zooxanthella_down|GO:0009310|zooxanthella_TRINITY_DN46713_C1_G2 zooxanthella_TRINITY_DN46713_C1_G3 zooxanthella_TRINITY_DN42914_C2_G1 zooxanthella_TRINITY_DN48316_C2_G1 zooxanthella_TRINITY_DN49474_C3_G1 zooxanthella_TRINITY_DN50127_C3_G1 zooxanthella_TRINITY_DN47858_C1_G7 zooxanthella_TRINITY_DN45516_C1_G3|amine catabolic process|3.49e-06|...|
|...|...|...|...|...|...|

<center>The second example of gene data</center>

|gene|log2FoldChange|baseMean|pvalue|padj|...|
|---|---|---|---|---|---|
|coral_TRINITY_DN42110_C3_G2|-8.294663081|8999.108427|4.56e-92|1.51e-87|...|
|coral_TRINITY_DN16010_C0_G1|-8.250704858|1143.301217|4.33e-52|7.16e-48|...|
|coral_TRINITY_DN49507_C3_G4|-7.387803252|4937.169229|4.04e-42|4.45e-38|...|
|coral_TRINITY_DN46498_C0_G2|-7.194028834|1316.470629|1.18e-40|9.74e-37|...|
|coral_TRINITY_DN45798_C1_G1|-7.629627563|2125.364107|2.39e-40|1.58e-36|...|
|coral_TRINITY_DN40312_C4_G7|-7.442776261|692.2305665|4.8e-40|2.64e-36|...|
|...|...|...|...|...|...|


&emsp;&emsp;If you do not have these files ready, you can use the demo data files by clicking on the "Load Demo" button. The first demo is from a metagenomic analysis of the human gut microbiomes from two cohorts. It contains differentially enriched enzymes in carbohydrate metabolism and information of related genes. The second demo is derived from a meta-transcriptomic analysis of the coral with infection. It contains differentially expressed coral/zooxanthella genes and over-represented GO terms.

<center>
<figure>
<img src='vignettes/figure/fig2.png' width='50%'>

Figure 2 Panel for loading demo

&nbsp;
</center>
</figure>

&emsp;&emsp;After the data files are uploaded and checked, they will be displayed in tables.

### 4.2 Data selection

&emsp;&emsp;You can filter term data in two steps. If the term data contains a column named ‘group’ to define different groups, you can select the group(s) first. Otherwise, skip the first step. Then, you can view data in an interactive table for further selection. Functional terms of interest can be selected and then displayed on the right table. These data will be used for plotting.

<center>
<figure>
<img src='vignettes/figure/fig3.png' width='100%'>

Figure 3 Data selection by two steps

&nbsp;
</center>
</figure>

### 4.3 Visualization

&emsp;&emsp;If the selected data comes from different groups (including the group column and the column has multiple values), you can select heatmap or dot chart for visualization. Otherwise, you can choose from bar plot, bubble chart, lollipop chart, network, orand complex bar plot. The font size of the text in axis and legend can be adjusted. The colors can be selected from color pickers or palettes. There are many color palettes from multiple R packages (ggsci, randomcoloR, scales, viridis, wesanderson) to choose from. You can download figures by clicking the download icon in the upper right corner.

&emsp;&emsp;After clicking on a term in a figure, you can get details about this functional term: (i) all contents in the row of it; (ii) more information of it from the web servers (if it is from GO, KEGG, or wikipathways); (iii) the data and visualization (pie chart, histogram, density plot, bar plot, bubble chart, box chart, heatmap, violin chart, and dumbbell chart) of the annotated genes.

&emsp;&emsp;These chart types will be described in more detail below.


## 5. Visualization of ungrouped terms

&emsp;&emsp;This part will use the first demo to demonstrate the visualization of ungrouped data. There is no column named group, so the first step of data selection is skipped and you can select data from a table. Then, the data can be displayed using multiple chart types.

<center>
<figure>
<img src='vignettes/figure/fig4.png' width='100%'>

Figure 4 Data selection in demo1

&nbsp;
</center>
</figure>

&emsp;&emsp;Each row in the bar plot represents a functional term, and you can select the attributes represented by the x-axis and colors. Functional terms can be divided into left and right parts according to a certain two-category attribute.	

<center>
<figure>
<img src='vignettes/figure/fig5.png' width='100%'>

Figure 5 Showing terms of demo1 in a bar plot

&nbsp;
</center>
</figure>

&emsp;&emsp;Each point in the bubble chart represents a functional term. You can select the attributes represented by the x-axis, y-axis, the color of the points, and the size of the points.

<center>
<figure>
<img src='vignettes/figure/fig6.png' width='100%'>

Figure 6 Showing terms of demo1 in a bubble chart

&nbsp;
</center>
</figure>

&emsp;&emsp;Each row in the lollipop chart represents a functional term, and you can select the attributes represented by the x-axis, the color of the points, and the size of the points.

<center>
<figure>
<img src='vignettes/figure/fig7.png' width='100%'>

Figure 7 Showing terms of demo1 in a lollipop chart

&nbsp;
</center>
</figure>

&emsp;&emsp;The x-axis of the complex bar plot represents the number of genes, and the y-axis is the functional term. Genes can be grouped into different blocks. The genes in bars can be sorted according to a certain attribute, and colors can be used to display a specific attribute.

<center>
<figure>
<img src='vignettes/figure/fig8.png' width='100%'>

Figure 8 Showing terms of demo1 in a complex barplot

&nbsp;
</center>
</figure>

&emsp;&emsp;If you click a term in the figure, more information will be shown in the bottom panel. Firstly, you can view the contents in the row of it.

<center>
<figure>
<img src='vignettes/figure/fig9.png' width='50%'>

Figure 9 Details of K01188

&nbsp;
</center>
</figure>

&emsp;&emsp;Then, if it is from GO, KEGG or WikiPathways database, more information (like name, class, description, diagram) will be retrieved from the web servers via APIs.

<center>
<figure>
<img src='vignettes/figure/fig10.png' width='100%'>

Figure 10 Description of K01188

&nbsp;
</center>
</figure>

&emsp;&emsp;Finally, you can examine annotated genes in a table and display them in multiple figures.

<center>
<figure>
<img src='vignettes/figure/fig11.png' width='100%'>

Figure 11 The data and visualization of genes annotated to K001188

&nbsp;

<img src='vignettes/figure/fig12-1.png' width='45%'>
<img src='vignettes/figure/fig12-2.png' width='45%'>

Figure 12 Occurrence frequency of genes annotated to K001188

&nbsp;

<img src='vignettes/figure/fig13-1.png' width='45%'>
<img src='vignettes/figure/fig13-2.png' width='45%'>

Figure 13 The length of genes annotated to K001188

&nbsp;

<img src='vignettes/figure/fig14.png' width='100%'>

Figure 14 Occurrence frequency and taxonomic annotations of genes annotated to K001188

&nbsp;

<img src='vignettes/figure/fig15.png' width='100%'>

Figure 15 The taxonomic classification of genes annotated to K001188

&nbsp;

</center>
</figure>


## 6. Visualization of multiple groups

&emsp;&emsp;This part will use the second demo to demonstrate the visualization of multiple groups or datasets.

&emsp;&emsp;If the selected data all come from the same group, the grouping information will be ignored, and you can observe data in bar plot, bubble chart, lollipop chart, network and complex bar plot (same as 5. Visualization of ungrouped terms).

<center>
<figure>
<img src='vignettes/figure/fig16.png' width='100%'>

Figure 16 Data selection in demo2

&nbsp;

<img src='vignettes/figure/fig17.png' width='100%'>

Figure 17 Showing selected terms of demo2 in a network

&nbsp;

<img src='vignettes/figure/fig18.png' width='100%'>

Figure 18 Showing selected terms of demo2 in a complex barplot

&nbsp;

</center>
</figure>

&emsp;&emsp;If the selected data comes from multiple groups, you can use heatmap or dot plot to display it.

<center>
<figure>
<img src='vignettes/figure/fig19.png' width='100%'>

Figure 19 Data selection in demo2

&nbsp;

<img src='vignettes/figure/fig20.png' width='100%'>

Figure 20 Showing terms of demo2 in a heatmap

&nbsp;

<img src='vignettes/figure/fig21.png' width='100%'>

Figure 21 Showing terms of demo2 in a dot plot

&nbsp;

</center>
</figure>

&emsp;&emsp;Similarly, you can click on the term in the picture to get relevant information.

<center>
<figure>
<img src='vignettes/figure/fig22.png' width='40%'>

Figure 22 Details of GO:0080134

&nbsp;

<img src='vignettes/figure/fig23.png' width='100%'>

Figure 23 Description of GO:0080134

&nbsp;

<img src='vignettes/figure/fig24.png' width='100%'>

Figure 24 Fold Change of genes annotated to GO:0080134

&nbsp;

<img src='vignettes/figure/fig25.png' width='100%'>

Figure 25 Reads count of genes annotated to GO:0080134

&nbsp;

</center>
</figure>

## Reference:
&emsp;&emsp;Li JH, Jia HJ, Cai XH, Zhong HZ, Feng Q, Sunagawa S, Arumugam M, Kultima JR, Prifti E, Nielsen T et al: An integrated catalog of reference genes in the human gut microbiome. Nat Biotechnol 2014, 32(8):834-841.

&emsp;&emsp;Zhou Z, Zhao SM, Tang J, Liu ZQ, Wu YB, Wang Y, Lin SJ: Altered Immune Landscape and Disrupted Coral-Symbiodinium Symbiosis in the Scleractinian Coral Pocillopora damicornis by Vibrio coralliilyticus Challenge. Front Physiol 2019, 10:12.


