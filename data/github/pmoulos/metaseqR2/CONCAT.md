<!-- badges: start -->
  ![Bioconductor build](http://www.bioconductor.org/shields/build/devel/bioc/metaseqR2.svg)
  ![Bioconductor platforms](http://www.bioconductor.org/shields/availability/3.12/metaseqR2.svg)
  ![Bioconductor dependencies](http://www.bioconductor.org/shields/dependencies/devel/metaseqR2.svg)
  </br>
  ![GitHub](https://img.shields.io/github/license/pmoulos/metaseqR2)
  ![GitHub repo size](https://img.shields.io/github/repo-size/pmoulos/metaseqR2)
  ![GitHub issues](https://img.shields.io/github/issues/pmoulos/metaseqR2)
<!-- badges: end -->

# metaseqR2

An R package for the analysis, meta-analysis and result reporting of RNA-Seq 
gene expression data - Next Generation!

## Citation

metaseqR2 along with further research regarding the abilities of the PANDORA 
algorithm was published in [Briefings in Bioinformatics](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa156/5890504). If you use metaseqR2 in your research, please cite:

> Dionysios Fanidis, Panagiotis Moulos: **Integrative, normalization-insusceptible statistical analysis of RNA-Seq data, with improved differential expression and unbiased downstream functional analysis**, *Briefings in Bioinformatics*, 2020, bbaa156, DOI: [10.1093/bib/bbaa156](https://doi.org/10.1093/bib/bbaa156)

## Installation from Bioconductor

```
if (!requireNamespace("BiocManager",quietly=TRUE))
    install.packages("BiocManager")

library(BiocManager)

BiocManager::install("metaseqR2")

# or for development version to be installed
# BiocManager::install("metaseqR2",version="devel")
```

## Installation from GitHub

Use with caution as the latest version may be unstable, although typical
Bioconductor checks are executed before each push.

```
if (!requireNamespace("devtools",quietly=TRUE))
    install.packages("devtools")

library(devtools)
install_github("pmoulos/metaseqR2")
```

## Installation from source

The same things apply regarding stability.

```
git clone https://github.com/pmoulos/metaseqR2.git
mkdir metaseqR2-build
rsync -avr --exclude=README.md --exclude=.git --exclude=.gitignore  \
    ./metaseqR2-local/ ./metaseqR2-build/metaseqR2
cd ./metaseqR2-build
R CMD build ./metaseqR2
```

This will take some time to build the vignettes. If you do not need them:

```
R CMD build --no-build-vignettes ./metaseqR2
```

And then install

```
R CMD INSTALL ./metaseqR2_x.y.z.tar.gz
```

Please report any issues [here](https://github.com/pmoulos/metaseqR2-local/issues). 

## metaseqR2 annotation database

If you do not wish to build annotation databases on your own using the
```buildAnnotationDatabase``` function, you can find complete pre-built 
annotation SQLite databases
[here](https://drive.google.com/drive/folders/15lOY9PBggCcaoohO_0rQTvExXenqah55?usp=sharing). 
New versions will be constructed from time to time, most probably whenever a new
Ensembl release comes live.

The prebuilt annotations contain:

* Annotations for all supported organsisms for all types of metaseqR2 analyses.
* For every supported organism:
  + For the latest version of each genome, the latest two Ensembl required
  annotations.
  + For all other versions of each genome, the latest Ensembl required
  annotations supporting that particular version.
  + UCSC and RefSeq annotations as fetched in the day of the build (denoted
  by the folder name in the above link).
  
The SQLite database must be placed in ```system.file(package="metaseqR2")``` and
named ```annotation.sqlite```, that is
```file.path(system.file(package="metaseqR2"),"annotation.sqlite")```. Otherwise
you will have to provide your desired location in each ```metaseqr2``` call.
Alternatively, on-the-fly download is still supported but is inneficient.

## List of required packages

metaseqR2 would benefit from the existence of all the following packages:

* ABSSeq
* baySeq
* Biobase
* BiocGenerics
* BiocManager
* BiocParallel
* BiocStyle
* biomaRt
* Biostrings
* BSgenome
* corrplot
* DESeq2
* DSS
* DT
* EDASeq
* edgeR
* harmonicmeanp
* genefilter
* GenomeInfoDb
* GenomicAlignments
* GenomicFeatures
* GenomicRanges
* gplots
* graphics
* grDevices
* heatmaply
* htmltools
* httr
* IRanges
* jsonlite
* knitr
* limma
* log4r
* magrittr
* Matrix
* methods
* NBPSeq
* pander
* parallel
* qvalue
* rmarkdown
* rmdformats
* RMySQL
* Rsamtools
* RSQLite
* rtracklayer
* RUnit
* S4Vectors
* splines
* stats
* stringr
* SummarizedExperiment
* survcomp
* TCC
* utils
* VennDiagram
* vsn
* zoo

A recent version of [Pandoc](https://pandoc.org/) is also required, ideally
above 2.0.
---
title: "metaseqR<sup>2</sup> report"
subtitle: "`r basename(PROJECT_PATH$main)`"
date: "`r format(Sys.time(), '%d/%m/%Y')`"
output:
  rmdformats::material:
    pandoc_args: ["+RTS","-K2048m","-RTS"]
    fig_width: 7
    fig_height: 7
    fig_caption: true
    highlight: kate
    lightbox: true
    thumbnails: true
    gallery: true
    cards: true
    use_bookdown: false
    mathjax: null
    self_contained: false
---

```{r knitr_init, echo=FALSE, cache=FALSE}
library(knitr)
library(rmdformats)

# Global options
options(max.print="75")
opts_chunk$set(
    echo=FALSE,
    cache=FALSE,
    prompt=FALSE,
    tidy=TRUE,
    comment=NA,
    message=FALSE,
    warning=FALSE,
    results="asis",
    eval=TRUE
)
opts_knit$set(width=75)
```

```{r load_js}
if (offlineReport) {
    #cat("<script type=\"text/javascript\" src=\"js/pace.min.js\"></script>",
    #   sep="")
    cat("<script type=\"text/javascript\" src=\"js/highcharts.js\"></script>",
        sep="")
    cat("<script type=\"text/javascript\" src=\"js/highcharts-more.js\">",
        "</script>",sep="")
    cat("<script type=\"text/javascript\" src=\"js/exporting.js\"></script>",
        sep="")
    cat("<script type=\"text/javascript\" src=\"js/offline-exporting.js\">",
        "</script>",sep="")
    cat("<script type=\"text/javascript\" src=\"js/export-data.js\"></script>",
        sep="")
    cat("<script type=\"text/javascript\" src=\"js/canvas2svg.js\"></script>",
        sep="")
    cat("<script type=\"text/javascript\" src=\"js/jvenn.min.js\"></script>",
        sep="")
    if (reportDb == "sqlite") {
        cat("<script>Module = { TOTAL_MEMORY: 134217728 };</script>")
        cat("<script type=\"text/javascript\" src=\"js/sql.js\"></script>",
            sep="")
    } else if (reportDb == "dexie") {
       cat("<script type=\"text/javascript\" src=\"js/dexie.min.js\"></script>",
            sep="")
    }
} else {
    #cat("<script type=\"text/javascript\" src=\"https://raw.github.com/",
    #   "HubSpot/pace/v1.0.0/pace.min.js\"></script>",sep="")
    cat("<script type=\"text/javascript\" src=\"https://code.highcharts.com/",
        "highcharts.js\"></script>",sep="")
    cat("<script type=\"text/javascript\" src=\"https://code.highcharts.com/",
        "highcharts-more.js\"></script>",sep="")
    cat("<script type=\"text/javascript\" src=\"https://code.highcharts.com/",
        "modules/exporting.js\"></script>",sep="")
    cat("<script type=\"text/javascript\" src=\"https://code.highcharts.com/",
        "modules/offline-exporting.js\"></script>",sep="")
    cat("<script type=\"text/javascript\" src=\"https://code.highcharts.com/",
        "modules/export-data.js\"></script>",sep="")
    cat("<script type=\"text/javascript\ ",
        "src=\"http://jvenn.toulouse.inra.fr/app/js/canvas2svg.js\">",
        "</script>",sep="")
    cat("<script type=\"text/javascript\" ",
        "src=\"http://jvenn.toulouse.inra.fr/app/js/jvenn.min.js\">",
        "</script>",sep="")
    if (reportDb == "sqlite") {
        cat("<script>Module = { TOTAL_MEMORY: 134217728 };</script>")
        cat("<script type=\"text/javascript\" src=\"",
        "https://cdnjs.cloudflare.com/ajax/libs/sql.js/0.5.0/js/sql.js\">",
        "</script>",sep="")
    } else if (reportDb == "dexie") {
        cat("<script type=\"text/javascript\" ",
            "src=\"https://unpkg.com/dexie@2.0.4/dist/dexie.min.js\">",
            "</script>",sep="")
    }
}
if (reportDb == "dexie")
    cat("<script type=\"text/javascript\" src=\"data/reportdb.js\"></script>")
```

```{js pace_opts}
Pace.options = {
    ajax: false,
    restartOnPushState: false,
    restartOnRequestAfter: false
}
```

<div id="loader" class="loading"></div>

```{r decide_dexie}
if (reportDb == "dexie") {
    code <- paste0(
        "var DEXIE = true;\n",
        "var db = new Dexie('plotdb');\n",
        "db.version(1).stores({\n",
        "    plots: '++_id,[name+type+subtype]'\n",
        "});\n",
        "db.plots.clear().then(console.log('Previous plots cleared!'));\n",
        "db.plots.bulkAdd(plotData);"
    )
} else if (reportDb == "sqlite") {
    code <- "var DEXIE = false;"
}
cat("<script>\n",code,"\n</script>")
```

```{js define_getplot_funs}
getPlotFromDb = function(query,renderId) {
    if (DEXIE) {
        $("#loader").show();
        db.plots.get({
            name: query.name,
            type: query.type,
            subtype: query.subtype
        },function(res) {
            var chartData = res.json;
            chartData.chart.renderTo = renderId;
            chartData.chart.events = {
                load: function() {
                    $("#loader").hide();
                }
            };
            var chart = new Highcharts.Chart(chartData);
        });
    }
    else {
        var xhr = new XMLHttpRequest();
        xhr.open('GET', './data/reportdb.sqlite',true);
        xhr.responseType = 'arraybuffer';
        xhr.onload = function(e) {
            var sqlQuery = "SELECT json FROM plot WHERE name='" + query.name + 
                "'" +" AND type='" + query.type + "'" + " AND subtype='" + 
                query.subtype + "';";
            var uInt8Array = new Uint8Array(this.response);
            var db = new SQL.Database(uInt8Array);
            var contents = db.exec(sqlQuery);
            var chartData = JSON.parse(contents[0].values[0][0]);
            //console.log(chartData)
            chartData.chart.renderTo = renderId;
            chartData.chart.events = {
                load: function() {
                    $("#loader").hide();
                }
            };
            var chart = new Highcharts.Chart(chartData);
        };
        xhr.send();
        $("#loader").show();
    }
}

getPlotWithFunFromDb = function(query,renderId) {
    if (DEXIE) {
        $("#loader").show();
        db.plots.get({
            name: query.name,
            type: query.type,
            subtype: query.subtype
        },function(res) {
            var chartData = res.json;
            chartData.chart.renderTo = renderId;
            
            // Deserialize the stored functions - only in a couple of cases
            if (query.type === "countsbio") {
                chartData.yAxis.labels.formatter = new Function('return ' +
                    chartData.yAxis.labels.formatter)()
                    
                chartData.plotOptions.scatter.tooltip.pointFormatter = 
                    new Function('return ' +
                        chartData.plotOptions.scatter.tooltip.pointFormatter)()
                
                chartData.plotOptions.boxplot.events.legendItemClick = 
                    new Function('return ' + 
                        chartData.plotOptions.boxplot.events.legendItemClick)()
            }
            if (query.type === "boxplot") {
                chartData.plotOptions.boxplot.tooltip.pointFormatter = 
                    new Function('return ' +
                        chartData.plotOptions.boxplot.tooltip.pointFormatter)()
                
                chartData.plotOptions.boxplot.events.legendItemClick = 
                    new Function('return ' + 
                        chartData.plotOptions.boxplot.events.legendItemClick)()
            }
            
            chartData.chart.events = {
                load: function() {
                    $("#loader").hide();
                }
            };
            var chart = new Highcharts.Chart(chartData);
        });
    }
    else {
        var xhr = new XMLHttpRequest();
        xhr.open('GET', './data/reportdb.sqlite',true);
        xhr.responseType = 'arraybuffer';
        xhr.onload = function(e) {
            var sqlQuery = "SELECT json FROM plot WHERE name='" + query.name + 
                "'" + " AND type='" + query.type + "'" + " AND subtype='" + 
                query.subtype + "';";
            var uInt8Array = new Uint8Array(this.response);
            var db = new SQL.Database(uInt8Array);
            var contents = db.exec(sqlQuery);
            var chartData = JSON.parse(contents[0].values[0][0]);
            
            // Deserialize the stored functions - only in a couple of cases
            if (query.type === "countsbio") {
                chartData.yAxis.labels.formatter = new Function('return ' +
                    chartData.yAxis.labels.formatter)()
                    
                chartData.plotOptions.scatter.tooltip.pointFormatter = 
                    new Function('return ' +
                        chartData.plotOptions.scatter.tooltip.pointFormatter)()
                
                chartData.plotOptions.boxplot.events.legendItemClick = 
                    new Function('return ' + 
                        chartData.plotOptions.boxplot.events.legendItemClick)()
            }
            if (query.type === "boxplot") {
                chartData.plotOptions.boxplot.tooltip.pointFormatter = 
                    new Function('return ' +
                        chartData.plotOptions.boxplot.tooltip.pointFormatter)()
                
                chartData.plotOptions.boxplot.events.legendItemClick = 
                    new Function('return ' + 
                        chartData.plotOptions.boxplot.events.legendItemClick)()
            }
            
            chartData.chart.renderTo = renderId;
            chartData.chart.events = {
                load: function() {
                    $("#loader").hide();
                }
            };
            var chart = new Highcharts.Chart(chartData);
        };
        xhr.send();
        $("#loader").show();
    }
}

getVennFromDb = function(query,renderId) {
    if (DEXIE) {
        $("#loader").show();
        db.plots.get({
            name: query.name,
            type: query.type,
            subtype: query.subtype
        },function(res) {
            var chartData = res.json;
            chartData.fnClickCallback = function() {
                var value = '';
                if (this.listnames.length == 1) {
                    value += 'Elements only in ';
                } else {
                    value += 'Common elements in ';
                }
                for (name in this.listnames) {
                    value += this.listnames[name] + ' ';
                }
                value += ':\n';
                for (val in this.list) {
                    value += this.list[val] + '\n';
                }
                if (query.subtype === "stat_dereg" 
                    || query.subtype === "stat_up"
                    || query.subtype === "stat_down") {
                    $('#jvenn_stat_list').text(value);
                }
                if (query.subtype === "fold_dereg" 
                    || query.subtype === "fold_up"
                    || query.subtype === "fold_down") {
                    $('#jvenn_fold_list').text(value);
                }
            }
            $('#'+renderId).jvenn(chartData);
            $("#loader").hide();
        });
    }
    else {
        var xhr = new XMLHttpRequest();
        xhr.open('GET', './data/reportdb.sqlite',true);
        xhr.responseType = 'arraybuffer';
        xhr.onload = function(e) {
            var sqlQuery = "SELECT json FROM plot WHERE name='" + query.name + 
                "'" +" AND type='" + query.type + "'" + " AND subtype='" + 
                query.subtype + "';";
            
            var uInt8Array = new Uint8Array(this.response);
            var db = new SQL.Database(uInt8Array);
            var contents = db.exec(sqlQuery);
            var chartData = JSON.parse(contents[0].values[0][0]);
            chartData.fnClickCallback = function() {
                var value = '';
                if (this.listnames.length == 1) {
                    value += 'Elements only in ';
                } else {
                    value += 'Common elements in ';
                }
                for (name in this.listnames) {
                    value += this.listnames[name] + ' ';
                }
                value += ':\n';
                for (val in this.list) {
                    value += this.list[val] + '\n';
                }
                if (query.subtype === "stat_dereg" 
                    || query.subtype === "stat_up"
                    || query.subtype === "stat_down") {
                    $('#jvenn_stat_list').text(value);
                }
                if (query.subtype === "fold_dereg" 
                    || query.subtype === "fold_up"
                    || query.subtype === "fold_down") {
                    $('#jvenn_fold_list').text(value);
                }
            }
            $('#'+renderId).jvenn(chartData);
            $("#loader").hide();
        };
        xhr.send();
        $("#loader").show();
    }
}
```

<style>
p {
    margin: 0 0 10px;
    text-align: justify;
}

.detable-container {
    font-size: 0.8em;
}

.header-panel {
    background-color: #0E2B8A;
    min-height: 144px;
    position: relative;
    z-index: 3;
}

.header-panel h1.subtitle {
    font-size: 16px;
    font-weight: 600;
}

.figure-hint {
    text-align: justify;
    font-size: 0.9em;
}

.hc-rect-sm {
    width: 640px;
    height: 640px;
    margin: auto;
}

.hc-rect-xs {
    width: 320px;
    height: 320px;
    margin: auto;
}

.hc-selector-label {
    color: #0030CE;
    font-size: 0.9em;
    font-weight: 600;
}

.deheatmap-container {
    width: 740px;
    height: 860px;
    text-align: center;
    margin: 25px auto 5px;
}

.coheatmap-container {
    width: 640px;
    height: 640px;
    text-align: center;
    margin: auto;
}

.link-box {
    background-color: #D5EEEB; 
    border-style: solid; 
    border-color: #00887A; 
    border-width: 1px; 
    width: 80%; 
    padding: 10px;
    margin-top: 20px;
    -moz-border-radius: 5px;
    -webkit-border-radius: 5px;
    -khtml-border-radius: 5px;
    border-radius: 5px;
}

.hidden-element {
    visibility: hidden;
    width: 10px;
    height: 10px;
}

.jvenn_list {
    font-size: 0.9em;
    font-wight: 600;
    width: 100%;
    height: 300px;
}

#summary_log_container {
    font-family: "Courier New", Courier, monospace;
    font-size: 0.9em;
    height: 400px;
    overflow-y: scroll;
    overflow-x: hidden;
    overflow-wrap: normal;
}

/******************************* Chart loading ********************************/
.loading {
  display: none;
  position: fixed;
  z-index: 999;
  height: 100%;
  width: 100%;
  overflow: show;
  margin: auto;
  top: 0;
  left: 0;
  bottom: 0;
  right: 0;
  background: radial-gradient(rgba(76, 51, 204, 0.7), rgba(0, 0, 0, .8));
  background: -webkit-radial-gradient(rgba(76, 51, 204, 0.7), rgba(0, 0, 0,.8));
}

.loading:before {
  content: 'Loading plots...';
  font-size: 92px;
  font-weight: 400;
  color: white;
  display: block;
  position: fixed;
  top: 40%;
  left: 35%;
  width: 100%;
  height: 100%;
}

/********************************** Pace.js ***********************************/
.pace {
  -webkit-pointer-events: none;
  pointer-events: none;

  -webkit-user-select: none;
  -moz-user-select: none;
  user-select: none;
}

.pace.pace-inactive .pace-progress {
  display: none;
}

.pace .pace-progress {
  position: fixed;
  z-index: 2000;
  top: 50%;
  left: 50%;
  height: 240px;
  width: 240px;
  margin-top: -120px; /* Negative half of height. */
  margin-left: -220px; /* Negative half of width. */

  -webkit-transform: translate3d(0, 0, 0) !important;
  -ms-transform: translate3d(0, 0, 0) !important;
  transform: translate3d(0, 0, 0) !important;
}

.pace .pace-progress:after {
  display: block;
  position: absolute;
  top: 0;
  content: attr(data-progress-text);
  font-family: "Helvetica Neue", sans-serif;
  font-weight: 200;
  font-size: 240px;
  line-height: 1;
  text-align: right;
  color: rgba(0, 0, 220, 0.3);
}

.pace .pace-progress:before {
  display: block;
  position: absolute;
  top: 0;
  content: url(media/dna_loader.gif);
  margin-left: -300px;
  margin-top: -120px;
}

.pace-running > *:not(.pace) {
  opacity:0;
}
/******************************** End Pace.js *********************************/
</style>

```{r report_init}
library(DT)
library(GenomicRanges)
library(gplots)
library(heatmaply)
library(htmltools)
library(magrittr)
library(pander)

makeFoldChange <- metaseqR2:::makeFoldChange
nat2log <- metaseqR2:::nat2log
asClassVector <- metaseqR2:::asClassVector
makeReportMessages <- metaseqR2:::makeReportMessages
.formatForReport <- metaseqR2:::.formatForReport

# Initialize some variables used throughout the report
samples <- unlist(sampleList,use.names=FALSE)
biotypes <- unique(as.character(geneData$biotype))
reportMessages <- makeReportMessages("en")

qualPlots <- c("mds","biodetection","countsbio","saturation","readnoise",
    "correl","pairwise","filtered")
normPlots <- c("boxplot","gcbias","lengthbias","meandiff","meanvar","rnacomp")
statPlots <- c("deheatmap","volcano","mastat","biodist","statvenn","foldvenn",
    "deregulogram")

covarsRaw <- covarsStat <- NULL
if (any(qcPlots %in% c("biodetection","countsbio","saturation","rnacomp",
    "readnoise"))) {
    covarsRaw <- list(
        data=geneCounts,
        length=width(geneData),
        gc=as.numeric(geneData$gc_content),
        chromosome=data.frame(
            chromosome=as.character(seqnames(geneData)),
            start=start(geneData),
            end=end(geneData)
        ),
        factors=data.frame(class=asClassVector(sampleList)),
        biotype=as.character(geneData$biotype),
        gene_name=as.character(geneData$gene_name)
    )
}

if ("biodist" %in% qcPlots) {
    covarsStat <- list(
        data=normGenesExpr,
        length=width(geneDataExpr),
        gc=as.numeric(geneDataExpr$gc_content),
        chromosome=data.frame(
            chromosome=as.character(seqnames(geneDataExpr)),
            start=start(geneDataExpr),
            end=end(geneDataExpr)
        ),
        factors=data.frame(class=asClassVector(sampleList)),
        biotype=as.character(geneDataExpr$biotype),
        gene_name=as.character(geneDataExpr$gene_name)
    )
}
```

# Instructions {.tabset .tabset-pills}

Welcome to the **metaseqR2** report! If you are familiar with the metaseqR 
report, then you will find that there are not many differences with respect to 
the presented information. Some diagnostic and exploration plots were added. The 
most notable difference is that __*all*__ plots are interactive. This helps a 
lot with exploration and interpretation but also adds a lot of computational
burden. However, relatively modern systems with recent browser versions should
be capable of rendering all the graphics. The metaseqR2 report has been tested 
with Google Chrome, Mozilla Firefox and Microsoft Edge. It has **not** been 
tested with Internet Explorer, Opera and Safari and most probably will not be. 
Other Chromium browsers (e.g. Brave) should also be fine.

One particular characteristic of the metaseqR2 report is that all plots are
interactive. This is achieved by using the standard graphics underlying data
with libraries including [Highcharts](https://www.highcharts.com/),
[Plotly](https://plot.ly/) and 
[jvenn](http://jvenn.toulouse.inra.fr/app/index.html) to create more 
user-friendly and directly explorable plots. Instructions on the usage of these 
plots follow:

* All plots are interactively **explorable**. This means that if you move your 
mouse **inside** the plot area (a move called *mouse-over*), you can retrieve 
information on each single data point. This applies to all plots. More 
specifically:
  + In **scatterplots**, if you mouse-over each point, information about this 
  point is presented, depending on the type of the plot. The data series from 
  which the point comes is also presented. For example, in a Volcano plot, fold 
  change and significance, as well as the name of the gene and the data category
  (e.g. up-regulated) will be presented.
  + In **barplots**, if you mouse-over each bar, information about this bar is 
  presented, such as the value it represents and the data series from which it
  comes. If the barplot contains groups of bars, then information about each 
  group is displayed. For example, in a Biodetection plot, each bar group
  presents the percentage of a biotype in the examined genome, the percentage
  in the sample and the detected percentage according to read counts.
  + In **boxplots**, if you mouse-over the boxes, the information about the 
  underlying distribution is displayed (maximum, upper quartile, median, lower
  quartile and minimum) as well as the data series. If you mouse-over an 
  outlier, then information on this single point is presented (e.g. value).
  + Some **barplots** have a **double** y-axis system corresponding to 
  **different** measurements or scales. For example in Biodetection barplots,
  the left y-axis presents abundant features while the right y-axis presents
  non-abundant features. In the Filtered barplot, y-axes present different 
  values (numbers and fractions).
  + **Line plots** can be moused-over too. Depending on the plot type, exact 
  values may or may not be shown, depending on how important it is to display 
  them, and to avoid over-crowding the plots. For example in Reads noise plot,
  we are interested in the trend and not so much in exact values.
  + **Heatmaps** can be moused-over too. Information on each heatmap cell will 
  be displayed.
* All scatterplots and heatmaps are **zoomable**. You need to press the left 
mouse button inside a plot area and draw a square area to zoom-in. If you wish 
to reset the zoom, there is a **button** appearing for this when zooming-in.
* Data series in scatterplots, barplots and boxplots can be **toggled** on or 
off by **clicking** on the **legend** name of each data series which is placed 
**below** each plot. For example, in Volcano plots, if you click on the name 
"Unregulated", then the respective data series will stop appearing in the plot. 
You can bring it **back** by clicking the legend again.
* All plots are **exportable**. On the **top right** corner of each scatterplot, 
barplot and boxplot, there is a menu button with several functionalities, 
including **exporting** in various formats and presenting the plot in 
**full-screen** mode. For heatmaps, this functionality is offered by a set of 
small buttons that appear if you mouse-over at the **top** of the heatmap.
* In Venn diagrams, if you **click** on the **number** for each category, the 
respective gene/transcript names will appear in the **box** on the right of the 
diagram.
* All plots can be **downloaded** in static formats (in formats according to 
```metaseqr2``` call) from the Results section.

The metaseqR2 report contains the sections described below, depending on which 
diagnostic and exploration plots have been asked for from the run command. As 
plots are categorized, if no plot from a specific category is asked for, then 
this category will not appear. Below, are the categories:

## Summary

The **Summary** section is further categorized in several subsections. 
Specifically:

* **Analysis summary:** This section contains an auto-generated text that 
analytically describes the computational process followed and summarized the
results of each step. This text can be used as is or with slight modifications
in the <u>Methods</u> section of an article.
* **Input options:** This section provides a list of the input arguments to the
pipeline in a more human-readable format.
* **Filtering:** This section reports in detail the number of filtered genes
decomposed according to the number of genes removed by each applied filter.
* **Differential expression:** This section reports in detail the number of
differentially expressed genes for each contrast, both when using only a p-value
cutoff as well as an FDR cutoff (numbers in parentheses), that is, genes passing
the multiple testing correction procedure selected. These numbers are also
calculated based on a simple fold change cutoff in log<sub>2</sub> scale. 
Finally, when multiple algorithms are used with p-value combination, this
section reports all the findings analytically per algorithm.
* **Command:** This section contains the command used to run the ```metaseqr2```
pipeline for users that want to experiment as well as a critical messages
displayed within the ```R``` session running ```metaseqr2``` displayed as a log.
Finally, if a targets file has been used to perform the analysis, a table
depicting the parameters in the targets files is created and a link to download
the actual targets file, but any relative paths to BAM files are stripped and
the user is responsible to prepend them if the targets file has to be reused in
another location, e.g. locally.
* **Tracks:** This section contains a link which opens a new window to the UCSC
Genome Browser where normalized tracks based on the input BAM files are 
displayed. If stranded tracks have been requested (according to the
sequencing protocol or technology), the a track hub is created to display the
stranded tracks. From this tab, you can also download bigWig files as well as 
copy track lines for manual input to the UCSC Genome Browser.

## Quality control 

The **Quality** control section contains several interactive plots concerning 
the overall quality control of each sample provided as well as overall 
assessments. The quality control plots are the *Multidimensional Scaling* 
(**MDS**) plot, the *Biotypes detection* (**Biodetection**) plot, the 
*Biotype abundance* (**Countsbio**) plot, the *Read saturation* (**Saturation**) 
plot, the *Read noise* (**ReadNoise**) plot, the *Correlation heatmap* 
(**Correlation**), the *Pairwise sample scatterplots* (**Pairwise**) and the 
*Filtered entities* (**Filtered**) plot. Each plot is accompanied by a detailed
description of what it depicts. Where multiple plot are available (e.g. one for
each sample), a selection list on the top of the respective  section allows the
selection of the sample to be displayed.

## Normalization

The **Normalization** section contains several interactive plots that can be 
used to inspect and assess the normalization procedure. Therefore, normalization
plots are usually paired, showing the same data instance normalized and not
normalized. The normalization plots are the *Expression boxplots* (**Boxplots**)
plots, the *GC content bias* (**GC bias**) plots, the *Gene length bias* 
(**Length bias**) plots, the *Within condition mean-difference*
(**Mean-Difference**) plots, the *Mean-variance relationship* 
(**Mean-Variance**) plot and the *RNA composition* (**Rna composition**) plot.
Each plot is accompanied by a detailed description of what it depicts. Where 
multiple plot are available (e.g. one for each sample), a selection list on the 
top of the respective section allows the selection of the sample to be 
displayed.

## Statistics

The **Statistics** section contains several interactive plots that can be used
to inspect and explore the outcome of statistical testing procedures. The 
statistics plots are the *Volcano plot* (**Volcano**), the 
*MA or Mean-Difference across conditions* (**MA**) plot, the 
*Expression heatmap* (**Heatmap**) plot, the 
*Chromosome and biotype distributions* (**Biodist**) plot, the
*Venn diagram across statistical tests* (**StatVenn**), the 
*Venn diagram across contrasts* (**FoldVenn**) and the *Deregulogram*. Each plot 
is accompanied by a detailed description of what it depicts. Please note that 
the heatmap plots only show the top percentage of differentially expressed genes 
as this is controlled by the ```reportTop``` parameter of the ```metaseqr2``` 
pipeline. When multiple plots are available (e.g. one for each contrast), a 
selection list on the top of the respective section allows the selection of the 
sample to be displayed.

## Results

The **Results** section contains a snapshot of differentially expressed genes in 
table format with basic information about each gene and links to external 
resources. Certain columns of the table are colored according to significance. 
Larger bars and more intense colors indicate higher significance. For example, 
the bar in the *p_value* column is larger if the genes has higher statistical 
significance and the fold change cell background is bright red if the gene is 
highly up-regulated. From the **Results** section, full gene lists can be 
downloaded in text tab-delimited format and viewed with a spreadsheet
application such as MS Excel. A selector on the top of the section above the 
table allows the display of different contrasts.

## References

The **References** section contains bibliographical references regading the
algorithms used by the ```metaseqr2``` pipeline and is adjusted according to the
algorithms selected.
 
# Summary {.tabset .tabset-pills}

## Analysis summary

### Analysis summary

```{r summary_analysis_summary}
if (fromRaw) {
    cat("The raw ",fileType," files, one for each RNA-Seq sample, were ",
        "summarized to ")
    if (transLevel == "gene") {
        if (countType=="exon") {
            cat("an exon ") 
        } else if (countType=="gene") {
            cat("a gene ")
        } else if (countType=="utr") {
            cat("a 3'UTR ")
        }
    } else if (transLevel=="transcript") {
        if (countType=="exon") {
            cat("an exon ")
        } else if (countType=="utr") {
            cat("a 3'UTR ")
        }
    } else if (transLevel=="exon") {
        if (countType=="exon")
            cat("an exon ")
    }
    cat("read counts table, ")
    if (fileType=="bam" || fileType=="sam") {
        cat("using the Bioconductor package GenomicRanges. ")
    } else if (fileType=="bed") {
        cat("using the Bioconductor packages rtracklayer and GenomicRanges. ")
    }
    cat("In the final read counts table, each row represented ")
    if (countType=="exon") {
        cat("one exon, ") 
    } else if (countType=="gene") {
        cat("one gene, ")
    }
    cat("each column one RNA-Seq sample and each cell, the corresponding ",
        "read counts associated with each row and column.")
}

if (countType == "exon") {
    if (!is.null(exonFilters)) {
        cat("The exon read counts were filtered for artifacts that could",
            "affect the subsequent normalization and statistical testing",
            "procedures as follows: ")
        msgExonFiltersMinActiveExons <- NULL
        if (!is.null(exonFilters$minActiveExons)) {
            msgExonFiltersMinActiveExons <- paste(
                "if an annotated ",transLevel," had up to ", 
                exonFilters$minActiveExons$exonsPerGene," exons, read ",
                "presence was required in at least ",
                exonFilters$minActiveExons$minExons," of the exons, else if ", 
                "an annotated ",transLevel," had more than ",
                exonFilters$minActiveExons$exonsPerGene," exons, then read ",
                "presence was required in at least ",
                exonFilters$minActiveExons$frac,
                "x&lceil;E&rceil; exons, where &lceil;<sup>.</sup>&rceil; is ",
                "the <em>ceiling</em> mathematical function. The application ",
                "of this filter resulted in the exclusion of ",
                length(exonFilterResult$minActiveExons)," ",transLevel,
                "s from further analysis. ",sep=""
            )
        }
        cat(msgExonFiltersMinActiveExons,sep=", ")
        cat("The total number of ",transLevel,"s excluded due to the ",
            "application of exon filters was ",
             sum(as.numeric(sapply(exonFilterResult,length))),". ",sep="")
    }
    cat("The final read counts for each ",transLevel," model were calculated",
        "as the sums of their exon reads, creating a ",transLevel," counts",
        "table where each row corresponded to an Ensembl",transLevel,
        " model and each column corresponded to an RNA-Seq sample. ")
}

cat("The ",transLevel," counts table was normalized for inherent systematic ",
    "or experimental biases (e.g. sequencing depth, ",transLevel," length, ",
    "GC content bias etc. using the Bioconductor package ",
     reportMessages$norm[[normalization]]," after removing ",transLevel,
    "s that had zero counts over all the RNA-Seq samples (",length(theZeros),
     " ",transLevel,"s). The output of the normalization algorithm ",sep="")
cat("was a table with normalized counts, which can be used for differential",
    "expression analysis with statistical algorithms developed specifically",
    "for count data. ")
    
if (!is.null(geneFilters)) {
    cat("Prior to the statistical testing procedure, the ",transLevel," read ",
    "counts were filtered for possible artifacts that may affect the ",
    "subsequent statistical testing procedures. Genes/transcripts presenting ",
    "any of the following were excluded from further analysis: ")
    
    latinNumbered <- c("i)","ii)","iii)","iv)","v)","vi)","vii)","viii)","ix)",
        "x)","xi)","xii)","xiii)","xiv)","xv)")
    counter <- 0
    msgGeneFiltersLength <- msgGeneFiltersExpressionMedian <- 
        msgGeneFiltersExpressionMean <- msgGeneFiltersExpressionQuantile <-
        msgGeneFiltersExpressionKnown <- msgGeneFiltersExpressionCustom <-
        msgGeneFiltersBiotype <- msgGeneFiltersAvgReads <- 
        msgGeneFiltersPresence <- NULL
    if (!is.null(geneFilters$length)) {
        counter <- counter + 1
        msgGeneFiltersLength <- paste(
            latinNumbered[counter]," ",transLevel,"s with length less than ", 
            geneFilters$length$length,"bp (",length(geneFilterResult$length),
            " genes)",sep=""
        )
    }
    if (!is.null(geneFilters$avgReads)) {
        counter <- counter + 1
        msgGeneFiltersAvgReads <- paste(
            latinNumbered[counter]," ",transLevel,"s whose average numbers of ",
            "reads per ",geneFilters$avgReads$averagePerBp," bp was less than ",
            "the ",100*geneFilters$avgReads$quantile,"<sup>th</sup> quantile ",
            "of the total normalized distribution of average reads per ",
            geneFilters$avgReads$averagePerBp,"bp (",
            length(geneFilterResult$averagePerBp)," ",transLevel,"s with ",
            "cutoff value ",round(geneFilterCutoff$avgReads,digits=5)," ",
            "average reads per ",geneFilters$avgReads$averagePerBp," bp)",sep=""
        )
    }
    if (!is.null(geneFilters$expression)) {
        if (!is.null(geneFilters$expression$median)
            && geneFilters$expression$median) {
            counter <- counter + 1
            msgGeneFiltersExpressionMedian <- paste(
                latinNumbered[counter]," ",transLevel,"s with read counts ",
                "below the median read counts of the total normalized count ",
                "distribution (",length(geneFilterResult$expression$median)," ",
                transLevel,"s with cutoff value ",
                geneFilterCutoff$expression$median," normalized read counts)",
                sep=""
            )
        }
        if (!is.null(geneFilters$expression$mean) 
            && geneFilters$expression$mean) {
            counter <- counter + 1
            msgGeneFiltersExpressionMean <- paste(
                latinNumbered[counter]," ",transLevel,"s with read counts ",
                "below the mean read counts of the total normalized counts ",
                "distribution (",length(geneFilterResult$expression$mean)," ",
                transLevel,"s with cutoff value ",
                geneFilterCutoff$expression$mean," normalized read counts)",
                sep=""
            )
        }
        if (!is.null(geneFilters$expression$quantile) 
            && !is.na(geneFilters$expression$quantile)) {
            counter <- counter + 1
            msgGeneFiltersExpressionQuantile <- paste(
                latinNumbered[counter]," ",transLevel,"s with read counts ",
                "below the ",100*geneFilters$expression$quantile,
                "<sup>th</sup> quantile of the normalized counts distribution ",
                "(",length(geneFilterResult$expression$quantile)," ",transLevel,
                "s with cutoff value ",geneFilterCutoff$expression$quantile," ",
                "normalized read counts)",sep=""
            )
        }
        if (!is.null(geneFilters$expression$known) 
            && !is.na(geneFilters$expression$known)) {
            counter <- counter + 1
            msgGeneFiltersExpressionKnown <- paste(
                latinNumbered[counter]," ",transLevel,"s with read counts ",
                "below the 90<sup>th</sup> quantile of the counts of the ",
                "following ",transLevel,"s, known to not being expressed from ",
                "the related literature: ",
                paste(geneFilters$expression$known,collapse=", "),
                "(",length(geneFilterResult$expression$known)," ",transLevel,
                "s with cutoff value",geneFilterCutoff$expression$known,
                " normalized read counts)",sep=""
            )
        }
        if (!is.null(geneFilters$expression$custom) 
            && !is.na(geneFilters$expression$custom)) {
            counter <- counter + 1
            msgGeneFiltersExpressionCustom <- paste(
                latinNumbered[counter]," genes not passing the user defined ",
                "filter provided to the metaseqr pipeline (",
                length(geneFilterResult$expression$custom)," ",transLevel,
                "s with cutoff value",geneFilterCutoff$expression$custom,")",
                sep=""
            )
        }
    }
    if (!is.null(geneFilters$biotype)) {
        counter <- counter + 1
        msgGeneFiltersBiotype <- paste(
            latinNumbered[counter]," ",transLevel,"s whose biotype matched ",
            "the following: ",
            paste(
                names(geneFilters$biotype)[which(unlist(geneFilters$biotype))],
                collapse=", "
            ),
            " (",length(geneFilterResult$biotype)," ",transLevel,"s)",sep=""
        )
    }
    if (!is.null(geneFilters$presence)) {
        counter <- counter + 1
        msgGeneFiltersPresence <- paste(
            latinNumbered[counter]," ",transLevel,"s which in ",
            geneFilters$presence$frac*100,"% of samples did not exceed ",
            geneFilters$presence$minCount," counts (",
            length(geneFilterResult$presence)," ",transLevel,"s)",sep=""
        )
        if (geneFilters$presence$perCondition) {
            msgGeneFiltersPresence <- 
                paste(msgGeneFiltersPresence," condition-wise",sep="")
        } else {
            msgGeneFiltersPresence <- 
                paste(msgGeneFiltersPresence," across all samples",sep="")
        }
    }
    
    cat(msgGeneFiltersLength,msgGeneFiltersAvgReads,
        msgGeneFiltersExpressionMedian,msgGeneFiltersExpressionMean,
        msgGeneFiltersExpressionQuantile,msgGeneFiltersExpressionKnown,
        msgGeneFiltersExpressionCustom,msgGeneFiltersBiotype,
        msgGeneFiltersPresence,sep=", ")
    cat(". The total number of ",transLevel,"s excluded due to the ",
        "application of ",transLevel," filters was ",
        sum(as.numeric(sapply(geneFilterResult,length))),". ",sep="")
    cat("The total (unified) number of ",transLevel,"s  excluded due to the ",
        "application of all filters was ",length(theZeros) + length(theDead),
        ". ",sep="")
}

cat("The resulting ",transLevel," counts table was subjected to differential ",
    "expression analysis for the contrasts ",
    paste(gsub("_vs_"," versus ",contrast),collapse=", "),sep="")
if (length(statistics)>1) {
    cat(" using the Bioconductor packages ",paste(unlist(
        reportMessages$stat[statistics],use.names=FALSE),collapse=", "),". ",
        sep="")
    if (metaP!="none") {
        cat("In order to combine the statistical significance from multiple ",
            "algorithms and perform meta-analysis, the ")
            cat(reportMessages$meta[[metaP]],"method was applied. ")
    }
} else {
    cat(" using the Bioconductor package ",reportMessages$stat[[statistics]],
        ". ",sep="")
}

if (!is.na(pcut))
    if (pcut==1) plasm <- 0.05 else plasm <- pcut

if (!is.null(contrast)) {
    cat("The final numbers of differentially expressed ",transLevel,"s were ",
        "(per contrast): ",sep="")
    msgContrast <- character(length(contrast))
    names(msgContrast) <- contrast

    for (cnt in contrast) {
        if (!is.na(pcut) && length(which(sumpList[[cnt]]<plasm))==0) {
            msgContrast[cnt] <- paste(
                "for the contrast ",gsub("_vs_"," versus ",cnt),
                ", no differentially expressed ",transLevel,"s were found",
                " with a p-value threshold of ",plasm,sep=""
            )
        } else if (is.na(pcut)) {
            msgContrast[cnt] <- paste(
                "for the contrast ",gsub("_vs_"," versus ",cnt),
                ", no statistical threshold defined",sep=""
            )
        } else {
            if (adjustMethod != "none") {
                addTextP <- paste(
                    " (",length(which(p.adjust(sumpList[[cnt]],
                    adjustMethod)<plasm)),")",sep=""
                )
                hasFdrText <- " (FDR or adjusted p-value)"
            } else {
                addTextP <- hasFdrText <- NULL
            }
            if (length(strsplit(cnt,"_vs_")[[1]]) == 2) {
                tmp <- log2(makeFoldChange(cnt,sampleList,normGenesExpr[which(
                    sumpList[[cnt]]<plasm),,drop=FALSE],logOffset))
                if (adjustMethod != "none") {
                    areThere <- which(p.adjust(sumpList[[cnt]],
                        adjustMethod)<plasm)
                    if (length(areThere) > 0) {
                        if (length(areThere)==1) {
                            tmpF <- log2(makeFoldChange(cnt,sampleList,
                                normGenesExpr[areThere,,drop=FALSE],logOffset))
                            addTextFu <- 
                                paste(" (",length(which(tmpF>=1)),")",sep="")
                            addTextFd <- 
                                paste(" (",length(which(tmpF<=-1)),")",sep="")
                            addTextFn <- 
                                paste(" (",length(which(abs(tmpF)<1)),")",
                                    sep="")
                        } else {
                            tmpF <- log2(makeFoldChange(cnt,sampleList,
                                normGenesExpr[areThere,,drop=FALSE],logOffset))
                            addTextFu <- 
                                paste(" (",length(which(tmpF>=1)),")",sep="")
                            addTextFd <- 
                                paste(" (",length(which(tmpF<=-1)),")",sep="")
                            addTextFn <- 
                                paste(" (",length(which(abs(tmpF)<1)),")",
                                    sep="")
                        }
                    } else {
                        addTextFu <- addTextFd <- addTextFn <- NULL
                    }
                } else {
                    addTextFu <- addTextFd <- addTextFn <- NULL
                }
                msgContrast[cnt] <- paste(
                    "for the contrast ",gsub("_vs_"," versus ",cnt),", ",
                    length(which(sumpList[[cnt]]<plasm)),addTextP,
                    " statistically significant ",transLevel,"s were found ",
                    "with a p-value",hasFdrText," threshold of ",plasm," and ",
                    "of these, ",length(which(tmp>=1)),addTextFu,
                    " were up-regulated, ",length(which(tmp<=-1)),addTextFd,
                    " were down-regulated and ",addTextFn,
                    " were not differentially expressed according to an ",
                    "absolute fold change cutoff value of 1 in ",
                    "log<sub>2</sub> scale",sep=""
                )
            } else {
                msgContrast[cnt] <- paste(
                    "for the contrast ",gsub("_vs_"," versus ",cnt),", ",
                    length(which(sumpList[[cnt]]<plasm)),addTextP,
                    " differentially expressed ",transLevel,"s were found ",
                    "with a p-value",hasFdrText," threshold of ",plasm," at ",
                    "least in one condition",sep=""
                )
            }
        }
    }
    cat(paste(msgContrast,collapse=", "))
} else {
    cat("No statistical testing was performed")
}

cat(". Literature references for all the algorithms used can be found at the ",
    "end of this report.")
```

## Input options

### Input options

**Read counts file: **
```{r summary_input_options_1}
cat(countsName)
```

**Conditions:** 
```{r summary_input_options_2}
cat(paste(names(sampleList),collapse=", "))
```

**Samples included:** 
```{r summary_input_options_3}
cat(paste(unlist(sampleList),collapse=", "))
```

**Samples excluded:** 
```{r summary_input_options_4}
if (!is.null(excludeList) && !is.na(excludeList)) {
    cat(paste(unlist(excludeList),collapse=", ")) 
} else {
    cat("none")
}

```

**Requested contrasts:** 
```{r summary_input_options_5}
if (!is.null(contrast)) {
    cat(paste(contrast,collapse=", "))
} else {
    cat("No contrasts were defined.")
}
```

**Library sizes:**
```{r summary_input_options_6}
if (!is.null(libsizeList)) {
    cat("<ul>")
    for (n in names(libsizeList)) 
        cat("<li>",paste("<em>",n,"</em>: ",libsizeList[[n]],sep=""),"</li>")
    cat("</ul>")
} else {
    cat("not available")
}
```

**Organism:** 
```{r summary_input_options_7}
cat(reportMessages$org[[org]])
```

**Annotation source:** 
```{r summary_input_options_8}
cat(reportMessages$refdb[[refdb]])
```

**Count type:** 
```{r summary_input_options_9}
cat(countType)
```

```{r summary_input_options_10}
if (countType=="utr" && !is.null(utrOpts)) {
    cat("<strong>3' UTR fraction:</strong> ",utrOpts$frac,"<br/>")
    cat("<strong>3' UTR minimum length:</strong> ",utrOpts$minLength,"bps<br/>")
    cat("<strong>3' UTR downstream:</strong> ",utrOpts$downstream,"bps<br/>")
}
```

**Exon filters:**
```{r summary_input_options_11}
if (!is.null(exonFilters)) {
    cat(paste(names(exonFilters),collapse=", "),"<br/>")
    for (ef in names(exonFilters)) {
        cat("<ul>")
        cat("<li><em><span style=\"font-size:1em\">",ef,"</span></em><ul>",
            sep="")
        for (efp in names(exonFilters[[ef]])) {
            if (length(exonFilters[[ef]][[efp]]) == 1 
                && is.function(exonFilters[[ef]][[efp]])) {
                cat("<li>custom function</li>")
            } else if (length(exonFilters[[ef]][[efp]]) == 1) {
                cat("<li>",paste(efp,exonFilters[[ef]][[efp]],sep=": "),
                    "</li>",sep="")
            } else if (length(exonFilters[[ef]][[efp]]) > 1) {
                cat("<li>",paste(efp,paste(exonFilters[[ef]][[efp]],
                    collapse=", "),sep=": "),"</li>",sep="")
            }
        }
        cat("</ul></li></ul>")
    }
} else {
    cat("none applied")
}
```

```{r summary_input_options_12}
if (!is.null(preset))
    cat("<strong>Analysis preset:</strong>",reportMessages$preset[[preset]])
```

**Gene filters:**
```{r summary_input_options_13}
if (!is.null(geneFilters)) {
    cat(paste(names(geneFilters),collapse=", "),"<br/>")
    for (gf in names(geneFilters)) {
        cat("<ul>")
        cat("<li><em><span style=\"font-size:1em\">",gf,"</span></em><ul>",
            sep="")
        for (gfp in names(geneFilters[[gf]])) {
            if (length(geneFilters[[gf]][[gfp]]) == 1 
                && is.function(geneFilters[[gf]][[gfp]])) {
                cat("<li>custom function</li>")
            } else if (length(geneFilters[[gf]][[gfp]]) == 1) {
                cat("<li>",paste(gfp,geneFilters[[gf]][[gfp]],sep=": "),"</li>",
                    sep="")
            } else if (length(geneFilters[[gf]][[gfp]]) > 1) {
                cat("<li>",paste(gfp,paste(geneFilters[[gf]][[gfp]],
                    collapse=", "),sep=": "),"</li>",sep="")
            }
        }
        cat("</ul></li></ul>")
    }
} else {
    cat("none applied")
}
```

**Filter application:** 
```{r summary_input_options_14}
cat(reportMessages$whenfilter[[whenApplyFilter]])
```

**Normalization algorithm:** 
```{r summary_input_options_15}
cat(reportMessages$norm[[normalization]])
```

**Normalization arguments:**
```{r summary_input_options_16}
if (!is.null(normArgs)) {
    if (normalization == "each") {
        for (n in names(normArgs)) {
            cat("<strong>Statistical arguments for",
                reportMessages$norm[[n]],": </strong>",
                paste(names(normArgs[[n]]),collapse=", "))
            if (length(normArgs[[n]]) > 0) {
                cat("<ul>")
                for (na in names(normArgs[[n]])) {
                    if (length(normArgs[[n]][[na]]) == 1 
                        && is.function(normArgs[[n]][[na]])) {
                        cat("<li>",na,": ",
                            as.character(substitute(normArgs[[na]])),
                        "</li>",sep="")
                    } else if (length(normArgs[[n]][[na]]) == 1) {
                        cat("<li>",paste(na,normArgs[[n]][[na]],sep=": "),
                            "</li>")
                    } else if (length(normArgs[[n]][[na]]) > 1) {
                        cat("<li>",paste(na,paste(normArgs[[n]][[na]],
                            collapse=", "),sep=": "),"</li>")
                    }
                }
                cat("</ul>")
            } else {
                cat("not available or not required")
            }
        }
    } else {
        cat(paste(names(normArgs),collapse=", "),"<br/>")
        cat("<ul>")
        for (na in names(normArgs)) {
            if (length(normArgs[[na]])==1 && is.function(normArgs[[na]])) {
                cat("<li>",as.character(substitute(normArgs[[na]])),"</li>",
                    sep="")
            } else if (length(normArgs[[na]]) == 1) {
                cat("<li>",paste(na,normArgs[[na]],sep=": "),"</li>")
            } else if (length(normArgs[[na]]) > 1) {
                cat("<li>",paste(na,paste(normArgs[[na]],collapse=", "),
                    sep=": "),"</li>")
            }
        }
        cat("</ul>")
    }
} else {
    cat("not available")
}
```

**Statistical algorithm(s):** 
```{r summary_input_options_17}
if (!any(is.na(statistics))) {
    cat(paste(unlist(reportMessages$stat[statistics],use.names=FALSE),
        collapse=", "))
} else {
    cat("No statitical testing performed.")
}
```

```{r summary_input_options_18}
if (!is.null(statArgs)) {
    for (s in names(statArgs)) {
        cat("<strong>Statistical arguments for ",reportMessages$stat[[s]],
            ": </strong>",paste(names(statArgs[[s]]),collapse=", "),sep="")
        if (length(statArgs[[s]]) > 0) {
            cat("<ul>")
            for (sa in names(statArgs[[s]])) {
                if (length(statArgs[[s]][[sa]]) == 1 
                    && is.function(statArgs[[s]][[sa]])) {
                    cat("<li>",sa,": ",as.character(substitute(statArgs[[na]])),
                        "</li>",sep="")
                } else if (length(statArgs[[s]][[sa]])==1) {
                    cat("<li>",paste(sa,statArgs[[s]][[sa]],sep=": "),"</li>")
                } else if (length(statArgs[[s]][[sa]])>1) {
                    cat("<li>",paste(sa,paste(statArgs[[s]][[sa]],
                        collapse=", "),sep=": "),"</li>")
                }
            }
            cat("</ul>")
        } else {
            cat("not available or not required<br/>")
        }
    }
} else {
    cat("<strong>Statistical arguments not available</strong>")
}
```

**Meta-analysis method:** 
```{r summary_input_options_19}
cat(reportMessages$meta[[metaP]])
```

**Multiple testing correction:** 
```{r summary_input_options_20}
cat(reportMessages$adjust[[tolower(adjustMethod)]])
```

**p-value threshold:** 
```{r summary_input_options_21}
if (!is.na(pcut)) cat(pcut) else cat("not available")
```

**Logarithmic tranformation offset: ** 
```{r summary_input_options_22}
cat(logOffset)
```

**Analysis preset:** 
```{r summary_input_options_23}
if (!is.null(preset)) cat(preset) else cat("not available")
```

**Quality control plots:**
```{r summary_input_options_24}
cat(paste(unlist(reportMessages$plots[qcPlots],use.names=FALSE),collapse=", "))
```

**Figure format:** 
```{r summary_input_options_25}
    cat(paste(figFormat,collapse=", "))
```

**Output directory:** 
```{r summary_input_options_26}
if (!is.na(exportWhere)) cat(exportWhere) else cat("default")
```

**Output data: **
```{r summary_input_options_27}
cat(paste(unlist(reportMessages$export[exportWhat],use.names=FALSE),
    collapse=", "))
```

**Output scale(s):** 
```{r summary_input_options_28}
cat(paste(unlist(reportMessages$export[exportScale],use.names=FALSE),
    collapse=", "))
```

**Output values:** 
```{r summary_input_options_29}
cat(paste(unlist(reportMessages$export[exportValues],use.names=FALSE),
    collapse=", "))
```

**Output statistics:**
```{r summary_input_options_30}
cat(paste(unlist(reportMessages$export[exportStats],use.names=FALSE),
    collapse=", "))
```

**Total run time:**
```{r summary_input_options_31}
cat(execTime)
```

## Filtering

<h3>Filtered 
```{r summary_filtered_what_1}
cat(transLevel,"s",sep="")
``` 
</h3>

**Number of filtered**
```{r summary_filtered_what_2}
cat("<strong>",transLevel,"s: </strong>",sep="")
```
```{r summary_filtered_howmany_1}
cat(length(theZeros) + length(theDead))
```
*which is the **union** of*

* Filtered because of zero reads:
```{r summary_filtered_howmany_2}
cat(length(theZeros))
```
* Filtered because of exon filters:
```{r summary_filtered_howmany_3}
cat(sum(as.numeric(sapply(exonFilterResult,length))))
if (sum(as.numeric(sapply(exonFilterResult,length))) != 0) {
    cat("<em> which is the <strong>union</strong> of</em>")
    cat("<ul>")
    for (n in names(exonFilterResult)) {
        cat("<li>",n,": ",length(exonFilterResult[[n]]),"</li>")
    }
    cat("</ul>")
}
```
* Filtered because of gene filters:
```{r summary_filtered_howmany_4}
cat(length(unique(unlist(geneFilterResult))))
```
*which is the **union** of*
```{r summary_filtered_howmany_5}
cat("<ul>")
for (n in names(geneFilterResult)) {
    if (!is.list(geneFilterResult[[n]])) {
        if (!is.null(geneFilterResult[[n]]))
            cat("<li>",n,": ",length(unlist(geneFilterResult[[n]])),
                " genes with filter cutoff value ",geneFilterCutoff[[n]],
                "</li>",sep="")
    } else {
        cat("<li>",n,": ",length(unlist(geneFilterResult[[n]])),
            " ",transLevel,"s further decomposed to (filter name, filtered ",
            transLevel,"s, filter cutoff):</li>",sep="")
        cat("<ul>")
        for (nn in names(geneFilterResult[[n]])) {
            if (!is.null(geneFilterResult[[n]][[nn]]))
                cat("<li>",nn,": ",length(unlist(geneFilterResult[[n]][[nn]])),
                " ",transLevel,"s with filter cutoff value ",
                paste(geneFilterCutoff[[n]][[nn]],collapse=", "),"</li>",sep="")
        }
        cat("</ul>")
    }
}
cat("</ul>")
```

## Differential expression

<h3>Differentially expressed 
```{r summary_de_what_1}
cat(transLevel,"s",sep="")
``` 
</h3>

**Number of differentially expressed **
```{r summary_de_what_2}
cat("<strong>",transLevel,"s per contrast: </strong>",sep="")
if (!is.null(contrast) && !any(is.na(statistics))) {
    if (!is.na(pcut) && pcut==1)
        cat("The p-value cutoff during the analysis was set to 1 so as to ",
            "retrieve the total ",transLevel," list which passed the filtering ",
            "procedure. Each ",transLevel," in the list is accompanied by its ",
            "statistical scores (p-value, FDR, etc.)")

    cat("<ul>")
    for (cnt in contrast) {
        cat("<li><strong>",cnt,": </strong>",sep="")
        if (!is.na(pcut) && length(which(sumpList[[cnt]]<pcut)) == 0) {
            cat("no differentially expressed ",transLevel,"s","</li>",sep="")
        } else if (is.na(pcut)) {
            cat("no statistical threshold defined","</li>")
        } else {
            if (!is.na(pcut) && pcut==1) {
                plasm <- 0.05
            } else {
                plasm <- pcut
            }
            if (adjustMethod!="none") {
                addTextP <- paste("(",length(which(p.adjust(sumpList[[cnt]],
                    adjustMethod) < plasm)),")",sep="")
                hasFdrText <- "(FDR or adjusted p-value)"
            } else {
                addTextP <- hasFdrText <- NULL
            }
            if (length(strsplit(cnt,"_vs_")[[1]]) == 2) {
                tmpP <- log2(makeFoldChange(cnt,sampleList,normGenesExpr[which(
                    sumpList[[cnt]] < plasm),,drop=FALSE],logOffset))
                if (adjustMethod != "none") {
                    areThere <- which(p.adjust(sumpList[[cnt]],adjustMethod)<plasm)
                    if (length(areThere)>0) {
                        if (length(areThere) == 1) {
                            tmpF <- log2(makeFoldChange(cnt,sampleList,
                                normGenesExpr[areThere,,drop=FALSE],logOffset))
                            addTextFu <- 
                                paste(" (",length(which(tmpF>=1)),")",sep="")
                            addTextFd <- 
                                paste(" (",length(which(tmpF<=-1)),")",sep="")
                            addTextFn <- 
                                paste(" (",length(which(abs(tmpF)<1)),")",sep="")
                        } else {
                            tmpF <- log2(makeFoldChange(cnt,sampleList,
                                normGenesExpr[areThere,,drop=FALSE],logOffset))
                            addTextFu <- 
                                paste(" (",length(which(tmpF>=1)),")",sep="")
                            addTextFd <- 
                                paste(" (",length(which(tmpF<=-1)),")",sep="")
                            addTextFn <- 
                                paste(" (",length(which(abs(tmpF)<1)),")",sep="")
                        }
                    } else {
                        addTextFu <- addTextFd <- addTextFn <- NULL
                    }
                } else {
                    addTextFu <- addTextFd <- addTextFn <- NULL
                }
                cat(length(which(sumpList[[cnt]]<plasm))," ",addTextP,
                    " statistically significant ",transLevel,"s of which ",
                    length(which(tmpP>=1)),addTextFu," up regulated, ",
                    length(which(tmpP<=-1)),addTextFd," down regulated and ",
                    length(which(abs(tmpP)<1)),addTextFn," not differentially ",
                    "expressed according to a p-value ",hasFdrText," threshold of ",
                    plasm," and an absolute fold change cutoff value of 1 in ",
                    "log<sub>2</sub> scale.",sep="")
            } else {
                cat(length(which(sumpList[[cnt]]<plasm))," ",addTextP,
                    " statistically significant, differentially expressed in at",
                    " least one condition at a p-value ",hasFdrText,
                    " threshold of ",plasm,".",sep="")
            }
            if (length(statistics)>1 && metaP != "none") {
                cat(" These numbers refer to the combined analysis performed by ",
                    "metaseqR2. Per statistical algorithm, the differentially ",
                    "expressed ",transLevel,"s are:",sep="")
                cat("<ul>")
                for (s in statistics) {
                    cat("<li><em>",reportMessages$stat[[s]],": </em>",sep="")
                    if (adjustMethod!="none") {
                        addTextP <- paste("(",length(which(p.adjust(
                            cpList[[cnt]][,s],adjustMethod)<plasm)),")",sep="")
                        hasFdrText <- "(FDR or adjusted p-value)"
                    } else {
                        addTextP <- hasFdrText <- NULL
                    }
                    if (length(strsplit(cnt,"_vs_")[[1]]) == 2) {
                        tmpP <- log2(makeFoldChange(cnt,sampleList,
                            normGenesExpr[which(cpList[[cnt]][,s]<plasm),,
                            drop=FALSE],logOffset))
                        if (adjustMethod != "none") {
                            tmpF <- log2(makeFoldChange(cnt,sampleList,
                                normGenesExpr[which(p.adjust(cpList[[cnt]][,s],
                                adjustMethod)<plasm),,drop=FALSE],logOffset))
                            addTextFu <- 
                                paste("(",length(which(tmpF>=1)),")",sep="")
                            addTextFd <- 
                                paste("(",length(which(tmpF<=-1)),")",sep="")
                            addTextFn <- 
                                paste("(",length(which(abs(tmpF)<1)),")",sep="")
                        } else {
                            addTextFu <- addTextFd <- addTextFn <- NULL
                        }
                        cat(length(which(cpList[[cnt]][,s]<plasm))," ",addTextP,
                            " statistically significant ",transLevel,"s of which ",
                            length(which(tmpP>=1))," ",addTextFu," up regulated, ",
                            length(which(tmpP<=-1)),addTextFd," down regulated ",
                            "and ",length(which(abs(tmpP)<1))," ",addTextFn,
                            " not differentially expressed according to a p-value ",
                            hasFdrText," threshold of ",plasm," and an absolute ",
                            "fold change cutoff value of 1 in log<sub>2</sub> ",
                            "scale.",sep="")
                    } else {
                        cat(length(which(cpList[[cnt]][,s]<plasm))," ",addTextP,
                            " statistically significant, differentially expressed ",
                            "in at least one condition at a p-value ",hasFdrText,
                            " threshold of ",plasm,".",sep="")
                    }
                    cat("</li>")
                }
                cat("</ul>")
            }
        }
    }
    cat("</ul>")
} else {
    cat("No statistical testing was performed. If contrasts were defined,",
        "then fold changes are available in the Results section.")
}
```

## Command

The differential expression analysis and this report were generated using the
following command:

```{r fun_call}
cat("<pre><code>")
cat(paste(FUN_CALL,collapse="<br/>"))
cat("</code></pre>")
```

```{r display_targets}
if (!is.null(theList)) {
    # First write targets
    tarDf <- metaseqR2:::.writeTargets(theList,outfile=file.path(
        PROJECT_PATH$main,"tracks","targets.txt"))
    tarDf$filename <- basename(as.character(tarDf$filename))
    
    # Download link
    cat("<br/>")
    cat("You can download the targets file from <strong>",
      "<a href='tracks/targets.txt' download target='_blank'>here</a>",
      "</strong><br/>",sep="")

    # Then render a table with targets
    cat("<br/>")
    cat("The following table summarizes the targets file used for the",
        "analysis. Do not forget to prepend the path to your BAM files in the",
        "<srong>filename</strong> column (also in the file that can be",
        "downloaded above).<br/>")
    print(kable(tarDf))
}
```

```{r summary_run_log}
if (runLog) {
    cat("<br/>The above command generated the following log output:<br/><br/>")
    cat("<div id=\"summary_log_container\">")
    logString <- paste(readLines(file.path(PROJECT_PATH$logs,
        "metaseqr_run.log")),collapse="---EOL---")
    logString <- gsub("\\$","->",logString)
    cat(gsub("---EOL---","<br/>",logString))
    cat("</div><br/>")
} else {
    cat("<h4>No logging requested</h4>")
}
```

## `r if (createTracks) 'Tracks'`

```{r tracks_1}
if (createTracks) {
    o <- metaseqR2:::getUcscOrganism(org)
    cat("You can use this ")
    if (trackInfo$stranded) {
        cat("<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?genome=",o,
            "&hubUrl=",tLink,"' target='_blank'><strong>link</strong>",
            "</a>",sep="")
    } else {
        cat("<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?db=",o,
            "&hgct_customText=",tLink,"' target='_blank'><strong>link</strong>",
            "</a>",sep="")
    }
    cat(" to load a UCSC Genome Browser session with the tracks derived from ",
        "this analysis. If stranded mode was chosen, a trackhub will be ",
        "loaded, otherwise, simple tracks will be loaded.",sep="")
    cat("<br/>")
}
```

```{r tracks_2}
if (createTracks) {
    if (trackInfo$stranded) {
        o <- metaseqR2:::getUcscOrganism(org)
        plusLinks <- paste0(trackInfo$urlBase,"/",o,"/",samples,"_plus.bigWig")
        minusLinks <- paste0(trackInfo$urlBase,"/",o,"/",samples,
            "_minus.bigWig")
        cat("You can download individual bigWig files, one for each sample,",
            "using the following list:<br/><br/>")
        cat("<strong>Plus (+) strand</strong><br/>")
        for (i in 1:length(plusLinks))
            cat("<a href='",plusLinks[i],"'>",samples[i],"</a><br/>",sep="")
        cat("<br/><strong>Minus (-) strand</strong><br/>")
        for (i in 1:length(minusLinks))
            cat("<a href='",minusLinks[i],"'>",samples[i],"</a><br/>",sep="")
        cat("<br/>")
    } else {
        tLines <- readLines(file.path(trackExportPath,"tracks.txt"))
        tmp <- strsplit(tLines," ")
        tNames <- sapply(tmp,function(x) {
            return(strsplit(x[3],"=")[[1]][2])
        })
        tLinks <- paste0(trackInfo$urlBase,"/",tNames,".bigWig")
        cat("You can download individual bigWig files, one for each sample,",
            "using the following list:<br/>")
        for (i in 1:length(tNames))
            cat("<a href='",tLinks[i],"'>",tNames[i],"</a><br/>",sep="")
        cat("<br/>")
        cat("You can manually load bigWig files in UCSC Genome Browser using ",
            "the following track lines:<br/><br/>")
        cat("<pre><code>")
        cat(paste(tLines,collapse="<br/>"))
        cat("</code></pre><br/>")
    }
}
```

```{r decide_qc}
if (any(qcPlots %in% qualPlots)) {
    pandoc.header("Quality control",level=1)
    pandoc.header("Quality control figures {.tabset .tabset-pills}",level=2)
    cat("The following figures summarize the quality control steps and ",
        "assessment performed by the ```metaseqr2``` pipeline. Each figure ",
        "category is accompanied by an explanatory text. All figures are ",
        "interactive wih additional controls on the top right of the figure.",
        sep="")
}
                
```

```{r analysis_figures_mds}
if ("mds" %in% qcPlots) {
    pandoc.header("MDS",level=3)
    
    ## Get chart data as serialized JSON
    #json <- diagplotMds(geneCounts,sampleList,output="json")
    
    cat("<h4><strong>Multidimensional scaling</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$mds,
        "<br/><br/></div>",sep="")
    
    # Container div
    cat("<div class=\"hc-rect-sm\" id=\"mds_container\"></div>",sep="")
    ## Render the chart
    #cat("<script>$('#mds_container').highcharts(",json,");</script>",sep="")
    # A hidden element to trigger upon page load and fetch the plot
    cat("<button id=\"_mds_trigger\" class=\"hidden-element\"></button>",sep="")
}
```

```{r analysis_figures_biodetection}
if ("biodetection" %in% qcPlots) {
    pandoc.header("Biodetection",level=3)
    
    cat("<h4><strong>Biotype detection</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$biodetection,
        "<br/><br/></div>",sep="")
        
    # Create a selector for each sample
    items <- paste("<option value=",paste("biodetection",samples,sep="_"),">",
        samples,"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a sample to display plot for",
        "</div>",sep="")
    cat("<select id=\"biodetection_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div
    cat("<div class=\"hc-rect-sm\" id=\"biodetection_container\"></div>",sep="")
    
    #if (reportDb == "dexie") {
    #   # Get chart data as list
    #   log <- capture.output({
    #       json <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
    #           whichPlot="biodetection",output="json")
    #   })
    #   
    #   # Add extra data, serialize and import the plots to Dexie
    #   preImport <- vector("list",length(json))
    #   for (i in 1:length(json)) {
    #       preImport[[i]] <- list(
    #           name=paste("biodetection",samples[i],sep="_"),
    #           type="biodetection",
    #           subtype="generic",
    #           json=json[[samples[i]]]
    #       )
    #   }
    #   toImport <- toJSON(preImport,auto_unbox=TRUE,null="null")
    #   
    #   # Import script
    #   cat("<script type=\"text/javascript\">")
    #   cat("db.plots.bulkAdd(",toImport,");",sep="")
    #   cat("</script>")
    #   
    #   # Initialize first plot with serialized JSON
    #   cat("<script>$('#biodetection_container').highcharts(",
    #       toJSON(preImport[[1]]$json,auto_unbox=TRUE,null="null"),
    #       ");</script>",sep="")
    #}
}
```

```{r analysis_figures_countsbio}
if ("countsbio" %in% qcPlots) {
    pandoc.header("Biocounts",level=3)
    
    cat("<h4><strong>Biotype representation</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$countsbio,
        "<br/><br/></div>",sep="")
    
    cat("<h4><em>Biotypes within samples</em></h4>")
    
    # Create a selector for each sample
    items <- paste("<option value=",paste("countsbio",samples,sep="_"),">",
        samples,"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a sample to display plot for",
        "</div>",sep="")
    cat("<select id=\"countsbio_sample_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div - sample
    cat("<div class=\"hc-rect-sm\" id=\"countsbio_sample_container\"></div>",
        sep="")
    
    cat("<br/><h4><em>Biotype representation across samples</em></h4>")
    
    # Create a selector for each biotype
    items <- paste("<option value=",paste("countsbio",biotypes,sep="_"),">",
        biotypes,"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a biotype to display plot for",
        "</div>",sep="")
    cat("<select id=\"countsbio_biotype_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div - biotype
    cat("<div class=\"hc-rect-sm\" id=\"countsbio_biotype_container\">",
        "</div>",sep="")    
}
```

```{r analysis_figures_saturation}
if ("saturation" %in% qcPlots) {
    pandoc.header("Saturation",level=3)
    
    cat("<h4><strong>Biotype representation</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$saturation,
        "<br/><br/></div>",sep="")
    
    cat("<h4><em>Read saturation per biotype for all samples</em></h4>")
    
    # Create a selector for each sample
    items <- paste("<option value=",paste("saturation",samples,sep="_"),">",
        samples,"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a sample to display plot for",
        "</div>",sep="")
    cat("<select id=\"saturation_sample_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div - sample
    cat("<div class=\"hc-rect-sm\" id=\"saturation_sample_container\"></div>",
        sep="")
    
    cat("<br/><h4><em>Read saturation per sample for all biotypes</em></h4>")
    
    # Create a selector for each biotype
    items <- paste("<option value=",paste("saturation",biotypes,sep="_"),">",
        biotypes,"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a biotype to display plot for",
        "</div>",sep="")
    cat("<select id=\"saturation_biotype_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div - biotype
    cat("<div class=\"hc-rect-sm\" id=\"saturation_biotype_container\">",
        "</div>",sep="")
}
```

```{r analysis_figures_readsnoise}
if ("readnoise" %in% qcPlots) {
    pandoc.header("Reads noise",level=3)
    
    cat("<h4><strong>RNA-Seq reads noise</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$readnoise,
        "<br/><br/></div>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"readsnoise_container\"></div>",
        sep="")
    cat("<button id=\"_readnoise_trigger\" class=\"hidden-element\"></button>",
        sep="")
}
```

```{r analysis_figures_correl}
if ("correl" %in% qcPlots) {
    pandoc.header("Correlation",level=3)
    
    cat("<h4><strong>Pairwise sample correlations</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$correl,
        "<br/><br/></div>",sep="")
    
    corMat <- cor(geneCounts)
    if (!is.null(colnames(geneCounts)))
        colnames(corMat) <- colnames(geneCounts)
    
    n <- dim(corMat)[1]
    labs <- matrix(NA,n,n)
    for (i in 1:n)
        for (j in 1:n)
            labs[i,j] <- sprintf("%.2f",corMat[i,j])
    if (n <= 5) {
        notecex <- 1.2
    } else if (n > 5 & n < 10) {
        notecex <- 0.9
    } else {
        notecex <- 0.7
    }
    
    h <- list(heatmaply(corMat,revC=TRUE,trace="none",symm=TRUE,Colv=TRUE,
        cellnote=labs,colors=colorRampPalette(c("yellow","grey","blue")),
        k_col=length(sampleList),k_row=length(sampleList)) %>%
        layout(width=640,height=480))
    setNames(h,NULL)
    htmltools::tagList(tags$div(id="cormap",class="coheatmap-container",h[[1]]))
    #heatmaply(corMat,revC=TRUE,trace="none",symm=TRUE,Colv=TRUE,
    #   cellnote=labs,colors=colorRampPalette(c("yellow","grey","blue")),
    #   k_col=length(sampleList),k_row=length(sampleList)) %>%
    #   layout(width=640,height=480)
}
```

```{r analysis_figures_pairwise}
if ("pairwise" %in% qcPlots) {
    pandoc.header("Pairwise",level=3)
    
    cat("<h4><strong>Pairwise sample correlations</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$pairwise,
        "<br/><br/></div>",sep="")
        
    nams <- colnames(geneCounts)
    n <- ncol(geneCounts)
    plotNames <- character(n*(n-1)/2)
    counter <- 0
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            counter <- counter + 1
            plotNames[counter] <- paste(nams[i],nams[j],sep="_vs_")
        }
    }
    
    # Create a selector for each sample
    items <- paste("<option value=",plotNames,">",gsub("_vs_"," vs ",plotNames),
        "</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a pair to display plot for",
        "</div>",sep="")
    cat("<select id=\"pairwise_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    cat("<h4><em>X-Y plots</em></h4>")
    # Container div - sample
    cat("<div class=\"hc-rect-sm\" id=\"pairwise_xy_container\"></div>",
        sep="")
    
    cat("<br/><h4><em>M-D plots</em></h4>")
    # Container div - sample
    cat("<div class=\"hc-rect-sm\" id=\"pairwise_md_container\"></div>",
        sep="")
}
```

```{r analysis_figures_filtered}
if ("filtered" %in% qcPlots) {
    pandoc.header("Filtered",level=3)
    
    cat("<h4><strong>Chromosome and biotype distribution of filtered ",
        transLevel,"s</strong></h4>",sep="")
    cat("<div class=\"figure-hint\">",reportMessages$explain$filtered,
        "<br/><br/></div>",sep="")
    cat("<h4>Chromosome distribution of filtered",transLevel,"s</h4>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"chromosome_distr_container\"></div>",
        sep="")
    cat("<h4>Biotype distribution of filtered ",transLevel,"s</h4>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"biotype_distr_container\"></div>",
        sep="")
    cat("<button id=\"_filtered_trigger\" class=\"hidden-element\"></button>",
        sep="")
}
```

```{r decide_norm}
if (any(qcPlots %in% normPlots)) {
    pandoc.header("Normalization",level=1)
    pandoc.header("Normalization assessment figures {.tabset .tabset-pills}",
        level=2)
    cat("The following figures allow for the assessment of the ",
        "normalization procedures performed by the ```metaseqr2``` pipeline. ",
        "Each figure category is accompanied by an explanatory text. All ",
        "figures are interactive wih additional controls on the top right ",
        "corner of the figure.",sep="")
}
```

```{r analysis_figures_boxplot}
if ("boxplot" %in% qcPlots) {
    pandoc.header("Boxplots",level=3)
    
    cat("<h4><strong>Boxplots</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$boxplot,
        "<br/><br/></div>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"boxplot_unorm_container\"></div>",
        sep="")
    cat("<div class=\"hc-rect-sm\" id=\"boxplot_norm_container\"></div>",
        sep="")
    cat("<br/><br/>")
    cat("<button id=\"_boxplot_trigger\" class=\"hidden-element\"></button>",
        sep="")
}
```

```{r analysis_figures_gcbias}
if ("gcbias" %in% qcPlots) {
    pandoc.header("GC bias",level=3)
    
    cat("<h4><strong>GC bias assessment plots</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$gcbias,
        "<br/><br/></div>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"gcbias_unorm_container\"></div>",
        sep="")
    cat("<div class=\"hc-rect-sm\" id=\"gcbias_norm_container\"></div>",
        sep="")
    cat("<br/><br/>")
    cat("<button id=\"_gcbias_trigger\" class=\"hidden-element\"></button>",
        sep="")
}
```

```{r analysis_figures_lengthbias}
if ("lengthbias" %in% qcPlots) {
    pandoc.header("Length bias",level=3)
    
    cat("<h4><strong>Length bias assessment plots</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$lengthbias,
        "<br/><br/></div>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"lengthbias_unorm_container\"></div>",
        sep="")
    cat("<div class=\"hc-rect-sm\" id=\"lengthbias_norm_container\"></div>",
        sep="")
    cat("<br/><br/>")
    cat("<button id=\"_lengthbias_trigger\" class=\"hidden-element\"></button>",
        sep="")
}
```

```{r analysis_figures_meandiff}
if ("meandiff" %in% qcPlots) {
    pandoc.header("Mean-Difference",level=3)
    
    cat("<h4><strong>Mean-difference plots for normalization assessment",
        "</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$meandiff,
        "<br/><br/></div>",sep="")
    
    # Create a selector for each sample
    namList <- vector("list",length(sampleList))
    names(namList) <- names(sampleList)
    for (n in names(sampleList)) {
        pairMatrix <- combn(1:length(sampleList[[n]]),2)
        nm <- character(0)
        for (i in 1:ncol(pairMatrix)) {
            s1 <- sampleList[[n]][pairMatrix[1,i]]
            s2 <- sampleList[[n]][pairMatrix[2,i]]
            nm <- c(nm,paste(s1,"_vs_",s2,sep=""))
        }
        namList[[n]] <- nm
    }
    
    html <- ""
    for (n in names(namList)) {
        nams <- namList[[n]]
        items <- paste("<option value=",nams,">",gsub("_vs_"," vs ",nams),
            "</option>",sep="")
        html <- paste(html,"<optgroup label=\"",n,"\">",items,"</optgroup>",
            sep="")
    }
    
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a pair to display plots for",
        "</div>",sep="")
    cat("<select id=\"meandiff_selector\">",html,"</select>",sep="")
    cat("</div>")
    
    # Container div - un-normalized
    cat("<div class=\"hc-rect-sm\" id=\"meandiff_unorm_container\"></div>",
        sep="")
    # Container div - normalized
    cat("<div class=\"hc-rect-sm\" id=\"meandiff_norm_container\"></div>",
        sep="") 
}
```

```{r analysis_figures_meanvar}
if ("meanvar" %in% qcPlots) {
    pandoc.header("Mean-Variance",level=3)
    
    cat("<h4><strong>Mean-variance plot for normalization assessment",
        "</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$meanvar,
        "<br/><br/></div>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"meanvar_unorm_container\"></div>",
        sep="")
    cat("<div class=\"hc-rect-sm\" id=\"meanvar_norm_container\"></div>",
        sep="")
    cat("<button id=\"_meanvar_trigger\" class=\"hidden-element\"></button>",
        sep="")
}
```

```{r analysis_figures_rnacomp}
if ("rnacomp" %in% qcPlots) {
    pandoc.header("Rna composition",level=3)
    
    cat("<h4><strong>RNA composition plot</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$rnacomp,
        "<br/><br/></div>",sep="")
        
    cat("<div class=\"hc-rect-sm\" id=\"rnacomp_unorm_container\"></div>",
        sep="")
    cat("<div class=\"hc-rect-sm\" id=\"rnacomp_norm_container\"></div>",sep="")
    cat("<button id=\"_rnacomp_trigger\" class=\"hidden-element\"></button>",
        sep="")
}
```

```{r decide_stat}
if (any(qcPlots %in% statPlots)) {
    pandoc.header("Statistics",level=1)
    pandoc.header(
        "Differential expression assessment figures {.tabset .tabset-pills}",
        level=2)
    cat("The following figures allow for the assessment of the statistical ",
        "testing procedures performed by the ```metaseqr2``` pipeline. Each ",
        "figure category is accompanied by an explanatory text. All figures ",
        "are interactive wih additional controls on the top right corner of ",
        "the figure.",sep="")
}
```

```{r analysis_figures_volcano}
if ("volcano" %in% qcPlots) {
    pandoc.header("Volcano",level=3)
    
    cat("<h4><strong>Volcano plots</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$volcano,
        "<br/><br/></div>",sep="")
        
    nn <- names(contrastList)
    namc <- character(0)
    for (n in nn) {
        fc <- log2(makeFoldChange(n,sampleList,normGenesExpr[1:3,],1))
        for (contr in colnames(fc))
            namc <- c(namc,contr)
    }
    
    # Create a selector for each contrast
    items <- paste("<option value=",paste("volcano",namc,sep="_"),">",
        gsub("_vs_"," vs ",namc),"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a contrast to display plot ",
        "for </div>",sep="")
    cat("<select id=\"volcano_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div
    cat("<div class=\"hc-rect-sm\" id=\"volcano_container\"></div>",sep="")
}
```

```{r analysis_figures_mastat}
if ("mastat" %in% qcPlots) {
    pandoc.header("MA",level=3)
    
    cat("<h4><strong>Mean-Difference (MA) plots</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$mastat,
        "<br/><br/></div>",sep="")
    
    nn <- names(contrastList)
    namc <- character(0)
    for (n in nn) {
        m <- log2(makeFoldChange(n,sampleList,normGenesExpr[1:3,],1))
        for (contr in colnames(m))
            namc <- c(namc,contr)
    }
    
    # Create a selector for each contrast
    items <- paste("<option value=",paste("mastat",namc,sep="_"),">",
        gsub("_vs_"," vs ",namc),"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a contrast to display plot ",
        "for </div>",sep="")
    cat("<select id=\"mastat_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div
    cat("<div class=\"hc-rect-sm\" id=\"mastat_container\"></div>",sep="")
}
```

```{r analysis_figures_deheatmap}
if ("deheatmap" %in% qcPlots) {
    pandoc.header("Heatmap",level=3)
    
    cat("<h4><strong>Differential expression heatmaps</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$deheatmap,
        "<br/><br/></div>",sep="")
    
    # Create a selector for each contrast
    items <- paste("<option value=",
        paste("deheatmap",gsub("\\.","__",names(contrastList)),sep="_"),">",
        gsub("_vs_"," vs ",names(contrastList)),"</option>",sep="")
    scal <- paste("<option value=",c("asis","zscore"),">",
        c("log2 reads","z-scores"),"</option>",sep="")
    cat("<div class=\"row\"><div class=\"col-sm-12\">",sep="")
    cat("<div class=\"row\"><div class=\"col-sm-6\">")
    cat("<div class=\"hc-selector-label\">Select a contrast to display ",
        "heatmap for </div>",sep="")
    cat("<select id=\"heatmap_contrast_selector\">",items,"</select>",sep="")
    cat("</div><div class=\"col-sm-6\">")
    cat("<div class=\"pull-right\">")
    cat("<div class=\"hc-selector-label\">Select a scale to display heatmap ",
        "for </div>",sep="")
    cat("<select id=\"heatmap_scale_selector\">",scal,"</select>",sep="")
    cat("</div></div></div></div></div>")
    
    hList <- vector("list",length(contrastList)*2)
    hNames <- character(length(contrastList)*2)
    hInd <- 0
    for (cnt in names(contrastList)) {
        s <- names(unlist(contrastList[[cnt]]))
        if (is.na(pcut))
            pcut <- 0.05
        ind <- which(sumpList[[cnt]]<=pcut)
        if (length(ind) > 0) {
            mat <- as.matrix(normGenesExpr[ind,s])
        } else {
            mat <- as.matrix(normGenesExpr[,s])
        }
        if (!is.null(reportTop)) {
            topg <- ceiling(reportTop*nrow(mat))
            mat <- mat[1:topg,,drop=FALSE]
        }
        # As is scale
        hInd <- hInd + 1
        hNames[hInd] <- paste0("deheatmap_",gsub("\\.","__",cnt),"_asis")
        hList[[hInd]] <- heatmaply(nat2log(mat),trace="none",colors=bluered(64),
            labRow=rep(NA,nrow(mat)),
            main=paste("DEG heatmap",cnt,"log2 reads")) %>% 
            layout(width=720,height=840)
        # z-scores
        M <- t(scale(t(nat2log(mat))))
        if (any(is.infinite(M)))
            M[which(is.infinite(M))] <- min(M)
        if (any(is.na(M)))
            M[which(is.na(M))] <- 0
        hInd <- hInd + 1
        hNames[hInd] <- paste0("deheatmap_",gsub("\\.","__",cnt),"_zscore")
        hList[[hInd]] <- heatmaply(M,trace="none",colors=bluered(64),
            labRow=rep(NA,nrow(mat)),
            main=paste("DEG heatmap",cnt,"z-scores")) %>% 
            layout(width=720,height=840)
    }
    setNames(hList,NULL)
    
    # Dirty hack
    toEval <- "htmltools::tagList("
    hItems <- character(length(hList))
    for (i in 1:length(hList))
        hItems[i] <- paste0("tags$div(id='",hNames[i],
            "', class='deheatmap-container',hList[[",i,"]])")
    toEval <- paste0(toEval,paste(hItems,collapse=","),")")
    eval(parse(text=toEval))
    
    #htmltools::tagList(
    #   tags$div(id=paste("deheatmap",gsub("\\.","__",names(contrastList)[1]),
    #       sep="_"),class="deheatmap-container",hList[[1]]),
    #   tags$div(id=paste("deheatmap",gsub("\\.","__",names(contrastList)[2]),
    #       sep="_"),class="deheatmap-container",hList[[2]])
    #)
}
```

```{js deheatmap_render}
// Chunk to initialize the display of the DE heatmaps
var hids = $('.deheatmap-container').map(function(index) {
    return(this.id);
});
$('.deheatmap-container').hide(0);
$('#'+hids[0]).show();
```

```{r analysis_figures_biodist}
if ("biodist" %in% qcPlots) {
    pandoc.header("Biodist",level=3)
    
    cat("<h4><strong>Chromosome and biotype distributions of differentially ",
        "expressed ",transLevel,"s</strong></h4>",sep="")
    cat("<div class=\"figure-hint\">",reportMessages$explain$biodist,
        "<br/><br/></div>",sep="")
    
    nn <- names(contrastList)
    
    # Create a selector for each contrast
    items <- paste("<option value=",paste("biodist",nn,sep="_"),">",
        gsub("_vs_"," vs ",nn),"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a contrast to display plot ",
        "for </div>",sep="")
    cat("<select id=\"biodist_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Containers   
    cat("<h4>Chromosome distribution of differentially expressed ",transLevel,
        "s</h4>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"biodist_chromosome_container\"></div>",
        sep="")
    cat("<h4>Biotype distribution of differentially expressed ",transLevel,
        "s</h4>",sep="")
    cat("<div class=\"hc-rect-sm\" id=\"biodist_biotype_container\"></div>",
        sep="")
}
```

```{r analysis_figures_deregulogram}
if ("deregulogram" %in% qcPlots) {
    pandoc.header("Deregulogram",level=3)
    
    cat("<h4><strong>Deregulograms</strong></h4>")
    cat("<div class=\"figure-hint\">",reportMessages$explain$deregulogram,
        "<br/><br/></div>",sep="")
    
    cntPairs <- combn(names(contrastList),2)
    namc <- character(ncol(cntPairs))
    counter <- 0
    for (i in 1:ncol(cntPairs)) {
        counter <- counter + 1
        namc[counter] <- paste0(cntPairs[1,i],"__",cntPairs[2,i])
    }
    
    # Create a selector for each pair of contrasts
    items <- paste("<option value=",paste("deregulogram",namc,sep="_"),">",
        gsub("__"," against ",namc),"</option>",sep="")
    cat("<div>")
    cat("<div class=\"hc-selector-label\">Select a pair of contrasts to ",
        "display plot for </div>",sep="")
    cat("<select id=\"deregulo_selector\">",items,"</select>",sep="")
    cat("</div>")
    
    # Container div
    cat("<div class=\"hc-rect-sm\" id=\"deregulo_container\"></div>",sep="")
}
```

```{r analysis_figures_statvenn}
if ("statvenn" %in% qcPlots) {
    pandoc.header("StatVenn",level=3)
    
    cat("<h4><strong>Venn diagram of differentially expressed ",transLevel,
        "s</strong></h4>",sep="")
    cat("<div class=\"figure-hint\">",reportMessages$explain$statvenn,
        "<br/><br/></div>",sep="")
        
    # Create a selector for each contrast
    nn <- names(contrastList)
    items <- paste("<option value=",paste("statvenn",nn,sep="_"),">",
        gsub("_vs_"," vs ",nn),"</option>",sep="")
    deupdo <- paste("<option value=",c("dereg","up","down"),">",
        c("Deregulated","Up regulated","Down downregulated"),"</option>",sep="")
    cat("<div class=\"row\"><div class=\"col-sm-8\">")
    cat("<div class=\"row\"><div class=\"col-sm-6\">")
    cat("<div class=\"hc-selector-label\">Select a contrast to display Venn ",
        "diagram for </div>",sep="")
    cat("<select id=\"statvenn_contrast_selector\">",items,"</select>",sep="")
    cat("</div><div class=\"col-sm-6\">")
    cat("<div class=\"hc-selector-label\">Select a direction to display Venn ",
        "diagram for </div>",sep="")
    cat("<select id=\"statvenn_direction_selector\">",deupdo,"</select>",sep="")
    cat("</div></div></div></div>")
    
    cat("<div class=\"row\">")
    cat("<div class=\"col-sm-8\"><div id=\"jvenn_stat_container\"></div></div>")
    cat("<div class=\"col-sm-4\">")
    cat("<div class=\"hc-selector-label\">Click on a number on Venn diagrams ",
        "to display the respective genes</div>",sep="")
    cat("<textarea id=\"jvenn_stat_list\" class=\"jvenn_list\" readonly>",
        "</textarea></div>",sep="")
    cat("</div>")
}

```{r analysis_figures_foldvenn}
if ("foldvenn" %in% qcPlots) {
    pandoc.header("FoldVenn",level=3)
    
    cat("<h4><strong>Venn diagram of differentially expressed ",transLevel,
        "s</strong></h4>",sep="")
    cat("<div class=\"figure-hint\">",reportMessages$explain$foldvenn,
        "<br/><br/></div>",sep="")
        
    # Create a selector for each contrast
    if (length(statistics) > 1) {
        nn <- c(statistics,metaP)
    } else {
        nn <- statistics
    }
    
    items <- paste("<option value=",paste("foldvenn",nn,sep="_"),">",nn,
        "</option>",sep="")
    deupdo <- paste("<option value=",c("dereg","up","down"),">",
        c("Deregulated","Up regulated","Down downregulated"),"</option>",sep="")
    cat("<div class=\"row\"><div class=\"col-sm-8\">")
    cat("<div class=\"row\"><div class=\"col-sm-6\">")
    cat("<div class=\"hc-selector-label\">Select an algorithm to display Venn ",
        "diagram for </div>",sep="")
    cat("<select id=\"foldvenn_algo_selector\">",items,"</select>",sep="")
    cat("</div><div class=\"col-sm-6\">")
    cat("<div class=\"hc-selector-label\">Select a direction to display Venn ",
        "diagram for </div>",sep="")
    cat("<select id=\"foldvenn_direction_selector\">",deupdo,"</select>",sep="")
    cat("</div></div></div></div>")
    
    cat("<div class=\"row\">")
    cat("<div class=\"col-sm-8\"><div id=\"jvenn_fold_container\"></div></div>")
    cat("<div class=\"col-sm-4\">")
    cat("<div class=\"hc-selector-label\">Click on a number on Venn diagrams ",
        "to display the respective genes</div>",sep="")
    cat("<textarea id=\"jvenn_fold_list\" class=\"jvenn_list\" readonly>",
        "</textarea></div>",sep="")
    cat("</div>")
}
```

# Results

## Tables of differentially expressed genes

The following tables allow for a quick exploration of the results of the
statistical analysis performed by the ```metaseqr2``` pipeline. If no
statistical testing or contrasts requested, just ignore any respective texts
and jump to tables or download the results.

```{r results_tables_header}
if (is.null(reportTop)) {
    topText <- paste("all the statistically significant ",  transLevel,"s. ",
        sep="")
} else {
    topText <- paste("the top ",round(100*reportTop),"% statistically ",
        "significant ",transLevel,"s. ",sep="")
}
restText <- paste("Use the download links below each table to retrieve the ",
    "total list of differentially expressed ",transLevel,"s or the whole ",
    transLevel," list of the selected genome irrespective of differential ",
    "expression.",sep="")
cat("Each table presents ",topText,restText,sep="")
cat("Furthermore each table can be searched using the search field on the ",
    "top right and you can also find the following information:<br/><br/>",
    sep="")

cat("<ul>")
cat("<li>The <em>chromosome</em> column is linked to the genomic location of ",
    "the ",transLevel," and opens a new tab/window to the UCSC Genome Browser",
    "</li>",sep="")
cat("<li>The *gene_id* column opens a link to the respective full annotation ",
    "source (only for Ensembl and RefSeq)</li>",sep="")
cat("<li>The background of the *p_value* and *FDR* columns displays a bar ",
    "with length proportional to the significance of each ",transLevel,"</li>",
    sep="")
cat("<li>The background color of the *fold change (_vs_)* column(s) displays ",
    "shows the deregulation of each ",transLevel," and is proportional to ",
    "the deregulation strength (red for up- green for down-regulation)</li>",
    sep="")
cat("<li>The background of the rest columns (condition average expression) ",
    "displays a bar with length proportional to the expression strength ",
    "of each condition</li>",sep="")
cat("</ul>")
```

```{r dt_hack, include=FALSE}
DT::datatable(NULL,extensions="Buttons")
```

```{r results_tables}
# Create a selector for each contrast
if (!is.null(contrast)) {
    items <- paste("<option value=",paste("detable",gsub("\\.","__",
        names(contrastList)),sep="_"),">",gsub("_vs_"," vs ",
        names(contrastList)),"</option>",sep="")
} else {
    items <- "<option value=detable_no_contrast>No contrast</option>"
}
cat("<hr><div>")
cat("<div class=\"hc-selector-label\">Select a contrast to display ",
    "DEG table for </div>",sep="")
cat("<select id=\"detable_selector\">",items,"</select>",sep="")
cat("</div><br/>")

if (!is.null(contrast) && !any(is.na(statistics))) {
    for (cnt in contrast) {
        cntn <- gsub("_vs_"," vs ",cnt)
        #pandoc.header(cntn,level=3)
        
        cat("<div class='detable-helper-wrapper' id='detable_",
            gsub("\\.","__",cnt),"'>",sep="")
        
        cat("<h4>DEG table for the contrast <strong>",cntn,"</strong></h4>",
            sep="")
        cat("<div class=\"figure-hint\">The following table presents ",topText,
            " for the contrast <strong>",cntn,"</strong>.</div>",sep="")
        
        # 15 colors should be enough for most situations...
        sColors <- c("#97F9F9","#A4DEF9","#C1E0F7","#CFBAE1","#C59FC9",
            "#EDEEC9","#DDE7C7","#BFD8BD","#98C9A3","#77BFA3","#7BDFF2",
            "#B2F7EF","#EFF7F6","#F7D6E0","#F2B5D4")
            
        # We do not need GC content really...
        gci <- grep("gc_content",colnames(reportTables[[cnt]]))
        if (length(gci) > 0)
            reportTables[[cnt]] <- reportTables[[cnt]][,-gci]
            
        # Fold change and condition means should always exist in this table 
        # defined by us...
        pi <- grep("p-value",colnames(reportTables[[cnt]]))
        di <- grep("FDR",colnames(reportTables[[cnt]]))
        fi <- grep("_vs_",colnames(reportTables[[cnt]]))
        # If more than one statistics applied, the 1st two from fi are meta 
        # p-value and meta FDR
        if (length(statistics) > 1)
            fi <- fi[-c(1,2)]
        #ci <- numeric(length(sampleList))
        #N <- names(sampleList)
        ci <- numeric(length(contrastList[[cnt]]))
        N <- unique(unlist(contrastList[[cnt]],use.names=FALSE))
        
        if (length(grep("rpgm_normalized_mean_counts_",
            colnames(reportTables[[cnt]]))) > 0) {
            for (i in 1:length(N))
                ci[i] <- grep(paste0("rpgm_normalized_mean_counts_",N[i]),
                    colnames(reportTables[[cnt]]))
        } else {
            for (i in 1:length(N))
                ci[i] <- grep(N[i],colnames(reportTables[[cnt]]),fixed=TRUE)
        }
        
        # Get values before rounding for styling
        brks <- quantile(reportTables[[cnt]][,fi],probs=seq(.05,.95,.05),
            na.rm=TRUE)
        preClrs <- col2rgb(rev(redgreen(length(brks)+1)))
        clrs <- paste("rgb(",apply(preClrs,2,paste0,collapse=","),")",sep="")
        pRangeForColor <- 1-as.numeric(reportTables[[cnt]][,pi])
        fRangeForColor <- 1-as.numeric(reportTables[[cnt]][,di])
        styleBars <- character(length(ci))
        for (i in 1:length(ci))
            styleBars[i] <- styleColorBar(range(reportTables[[cnt]][,ci[i]]),
                sColors[i])
        
        # After setting colors, format the table for the report
        if (!is.numeric(version)) # Not defined in call
            version <- NULL
        reportTables[[cnt]] <- .formatForReport(reportTables[[cnt]],o=org,
            r=refdb,v=version)
        # Hack to prevent the header name taking 2 rows because of 
        # hyphenation(?)
        names(reportTables[[cnt]])[pi] <- "p_value"
        
        # Dirty hack
        # https://stackoverflow.com/questions/32018521/
        # shiny-use-stylecolorbar-with-data-from-two-data-frames
        reportTables[[cnt]] <- 
            cbind(pRangeForColor,fRangeForColor,reportTables[[cnt]])
        jspf <- paste0("function(row,data) {",
            "var value = data[0];",
            "var backgroundValueP =",styleColorBar(pRangeForColor,
                "lightblue")[1],";",
            "var value = data[1];",
            "var backgroundValueF =",styleColorBar(fRangeForColor,
                "lightsalmon")[1],";",
            "$('td',row).eq(",pi-1,").css({",
                "'background': backgroundValueP,",
                "'background-repeat': 'no-repeat',",
                "'background-position': 'center',",
                "'background-size': '98% 88%'",
            "});",
            "$('td',row).eq(",di-1,").css({",
                "'background': backgroundValueF,",
                "'background-repeat': 'no-repeat',",
                "'background-position': 'center',",
                "'background-size': '98% 88%'",
            "});",
            "for (var i=",fi[1],";i<=",fi[length(fi)],";i++) {",
                "var value = data[i+1];",
                "$('td',row).eq(i-1).css({",
                    "'color': '#FFFFFF',",
                    "'font-weight': 'bold',",
                    "'background-color': ",styleInterval(brks,clrs),
                "});",
            "}")
        
        for (i in 1:length(ci)) {
            jsa <- paste0("var value = data[",ci[i]+1,"];",
                "$('td',row).eq(",ci[i]-1,").css({",
                    "'background': ",styleBars[i],",",
                    "'background-repeat': 'no-repeat',",
                    "'background-position': 'center',",
                    "'background-size': '98% 88%'",
                "});")
            jspf <- paste0(jspf,jsa)
        }
        jspf <- paste0(jspf,"}")
        
        theTable <- DT::datatable(
            reportTables[[cnt]],
            rownames=FALSE,
            width="95%",
            height=460,
            style="bootstrap",
            class="table-condensed table-striped table-responsive",
            escape=FALSE,
            elementId=paste("detable_",gsub("\\.","__",cnt),"_",sep=""),
            extensions="Buttons",
            options=list(
                scrollX=TRUE,
                order=list(list(pi+1,"asc")),
                columnDefs=list(
                    list(
                        targets=0:1,
                        visible=FALSE
                    )
                ),
                dom='Bfrtip',
                #buttons=I('colvis'),
                buttons=list(
                    list(
                        extend='colvis',
                        #exclude=0:1,
                        text="Show/Hide columns"
                    )
                ),
                rowCallback=JS(jspf)
            )
        ) #%>% formatStyle(
            #names(reportTables[[cnt]])[fi+2],
            #backgroundColor=styleInterval(brks,clrs),
            #fontWeight="bold",
            #color="#FFFFFF"
        #) #%>% formatStyle(
            #names(reportTables[[cnt]])[pi+2],
            #background=styleColorBar(pRangeForColor,"lightblue"),
            #fontWeight="bold",
            #backgroundSize="98% 88%",
            #backgroundRepeat="no-repeat",
            #backgroundPosition="center"
        #)
        print(htmltools::tagList(tags$div(class="detable-container",theTable)))
        
        cat("<br/><br/><div class=\"figure-hint link-box\">")
        cat("<strong><a href=\"lists/metaseqr_sig_out_",cnt,".txt.gz\"",
            " download>Download</a> the DEG result list for ",cnt,
            ".</strong><br/>",sep="")
        if (!is.null(geneCountsZero) || !is.null(geneCountsDead))
            cat("<strong><a href=\"lists/metaseqr_all_out_",cnt,".txt.gz\"",
            " download>Download</a> the whole result list for ",cnt,
            ".</strong><br/>",sep="")
        cat("</div><br/><br/>")
        
        cat("</div>")
    }
} else if (!is.null(contrast) && any(is.na(statistics))) {
    for (cnt in contrast) {
        cntn <- gsub("_vs_"," vs ",cnt)
        
        cat("<div class='detable-helper-wrapper' id='detable_",
            gsub("\\.","__",cnt),"'>",sep="")
        
        cat("<h4>DEG table for the contrast <strong>",cntn,"</strong></h4>",
            sep="")
        cat("<div class=\"figure-hint\">The following table presents ",topText,
            " for the contrast <strong>",cntn,"</strong>.</div>",sep="")
        
        # 15 colors should be enough for most situations...
        sColors <- c("#97F9F9","#A4DEF9","#C1E0F7","#CFBAE1","#C59FC9",
            "#EDEEC9","#DDE7C7","#BFD8BD","#98C9A3","#77BFA3","#7BDFF2",
            "#B2F7EF","#EFF7F6","#F7D6E0","#F2B5D4")
            
        # We do not need GC content really...
        gci <- grep("gc_content",colnames(reportTables[[cnt]]))
        if (length(gci) > 0)
            reportTables[[cnt]] <- reportTables[[cnt]][,-gci]
            
        # Fold change and condition means should always exist in this table 
        # defined by us...
        fi <- grep("_vs_",colnames(reportTables[[cnt]]))
        #ci <- numeric(length(contrastList[[cnt]]))
        #N <- unique(unlist(contrastList[[cnt]],use.names=FALSE))
        ci <- numeric(length(names(sampleList)))
        N <- names(sampleList)
        
        if (length(grep("rpgm_normalized_mean_counts_",
            colnames(reportTables[[cnt]]))) > 0) {
            for (i in 1:length(N))
                ci[i] <- grep(paste0("rpgm_normalized_mean_counts_",N[i]),
                    colnames(reportTables[[cnt]]))
        } else {
            for (i in 1:length(N))
                ci[i] <- grep(N[i],colnames(reportTables[[cnt]]),fixed=TRUE)
        }
        
        # Get values before rounding for styling
        brks <- quantile(reportTables[[cnt]][,fi],probs=seq(.05,.95,.05),
            na.rm=TRUE)
        preClrs <- col2rgb(rev(redgreen(length(brks)+1)))
        clrs <- paste("rgb(",apply(preClrs,2,paste0,collapse=","),")",sep="")
        styleBars <- character(length(ci))
        for (i in 1:length(ci))
            styleBars[i] <- styleColorBar(range(reportTables[[cnt]][,ci[i]]),
                sColors[i])
        
        # After setting colors, format the table for the report
        if (!is.numeric(version)) # Not defined in call
            version <- NULL
        reportTables[[cnt]] <- .formatForReport(reportTables[[cnt]],o=org,
            r=refdb,v=version)
        
        # Dirty hack
        # https://stackoverflow.com/questions/32018521/
        # shiny-use-stylecolorbar-with-data-from-two-data-frames
        reportTables[[cnt]] <- cbind(reportTables[[cnt]])
        jspf <- paste0("function(row,data) {",
            "for (var i=",fi[1],";i<=",fi[length(fi)],";i++) {",
                "var value = data[i-1];",
                "$('td',row).eq(i-1).css({",
                    "'color': '#FFFFFF',",
                    "'font-weight': 'bold',",
                    "'background-color': ",styleInterval(brks,clrs),
                "});",
            "}")
        
        for (i in 1:length(ci)) {
            jsa <- paste0("var value = data[",ci[i]-1,"];",
                "$('td',row).eq(",ci[i]-1,").css({",
                    "'background': ",styleBars[i],",",
                    "'background-repeat': 'no-repeat',",
                    "'background-position': 'center',",
                    "'background-size': '98% 88%'",
                "});")
            jspf <- paste0(jspf,jsa)
        }
        jspf <- paste0(jspf,"}")
        
        theTable <- DT::datatable(
            reportTables[[cnt]],
            rownames=FALSE,
            width="95%",
            height=460,
            style="bootstrap",
            class="table-condensed table-striped table-responsive",
            escape=FALSE,
            elementId=paste("detable_",gsub("\\.","__",cnt),"_",sep=""),
            extensions="Buttons",
            options=list(
                scrollX=TRUE,
                dom='Bfrtip',
                buttons=list(
                    list(
                        extend='colvis',
                        text="Show/Hide columns"
                    )
                ),
                rowCallback=JS(jspf)
            )
        )
        print(htmltools::tagList(tags$div(class="detable-container",theTable)))
        
        cat("<br/><br/><div class=\"figure-hint link-box\">")
        cat("<strong><a href=\"lists/metaseqr_sig_out_",cnt,".txt.gz\"",
            " download>Download</a> the DEG result list for ",cnt,
            ".</strong><br/>",sep="")
        if (!is.null(geneCountsZero) || !is.null(geneCountsDead))
            cat("<strong><a href=\"lists/metaseqr_all_out_",cnt,".txt.gz\"",
            " download>Download</a> the whole result list for ",cnt,
            ".</strong><br/>",sep="")
        cat("</div><br/><br/>")
        
        cat("</div>")
    }
} else {
    cat("<div class='detable-helper-wrapper' id='detable_no_contrast'>")
    cat("<h4>Quantification table for the run</h4>")
    
    # 15 colors should be enough for most situations...
    sColors <- c("#97F9F9","#A4DEF9","#C1E0F7","#CFBAE1","#C59FC9",
        "#EDEEC9","#DDE7C7","#BFD8BD","#98C9A3","#77BFA3","#7BDFF2",
        "#B2F7EF","#EFF7F6","#F7D6E0","#F2B5D4")
        
    # We do not need GC content really...
    gci <- grep("gc_content",colnames(reportTables[[1]]))
    if (length(gci) > 0)
        reportTables[[1]] <- reportTables[[1]][,-gci]
        
    #ci <- numeric(length(unlist(sampleList)))
    #N <- unique(unlist(sampleList,use.names=FALSE))
    ci <- numeric(length(names(sampleList)))
    N <- names(sampleList)
    
    if (length(grep("rpgm_normalized_mean_counts_",
        colnames(reportTables[[1]]))) > 0) {
        for (i in 1:length(N))
            ci[i] <- grep(paste0("rpgm_normalized_mean_counts_",N[i]),
                colnames(reportTables[[1]]))
    } else {
        for (i in 1:length(N))
            ci[i] <- grep(N[i],colnames(reportTables[[1]]),fixed=TRUE)
    }
    
    # Get values before rounding for styling
    styleBars <- character(length(ci))
    for (i in 1:length(ci))
        styleBars[i] <- styleColorBar(range(reportTables[[1]][,ci[i]]),
            sColors[i])
    
    # After setting colors, format the table for the report
    if (!is.numeric(version)) # Not defined in call
        version <- NULL
    reportTables[[1]] <- .formatForReport(reportTables[[1]],o=org,r=refdb,
        v=version)
    
    # Dirty hack
    # https://stackoverflow.com/questions/32018521/
    # shiny-use-stylecolorbar-with-data-from-two-data-frames
    jspf <- paste0("function(row,data) {")
    for (i in 1:length(ci)) {
        jsa <- paste0("var value = data[",ci[i]-1,"];",
            "$('td',row).eq(",ci[i]-1,").css({",
                "'background': ",styleBars[i],",",
                "'background-repeat': 'no-repeat',",
                "'background-position': 'center',",
                "'background-size': '98% 88%'",
            "});")
        jspf <- paste0(jspf,jsa)
    }
    jspf <- paste0(jspf,"}")
    
    theTable <- DT::datatable(
        reportTables[[1]],
        rownames=FALSE,
        width="95%",
        height=460,
        style="bootstrap",
        class="table-condensed table-striped table-responsive",
        escape=FALSE,
        elementId="detable_no_contrast",
        extensions="Buttons",
        options=list(
            scrollX=TRUE,
            dom='Bfrtip',
            buttons=list(
                list(
                    extend='colvis',
                    text="Show/Hide columns"
                )
            ),
            rowCallback=JS(jspf)
        )
    )
    print(htmltools::tagList(tags$div(class="detable-container",theTable)))
    
    cat("<br/><br/><div class=\"figure-hint link-box\">")
    cat("<strong><a href=\"lists/metaseqr_sig_out.txt.gz\"",
        " download>Download</a> the result list.</strong><br/>",sep="")
    if (!is.null(geneCountsZero) || !is.null(geneCountsDead))
        cat("<strong><a href=\"lists/metaseqr_all_out.txt.gz\"",
        " download>Download</a> the whole result list.</strong><br/>",sep="")
    cat("</div><br/><br/>")
    
    cat("</div>")
}
```

```{js detable_render}
var tids = $('.detable-helper-wrapper').map(function(index) {
    return(this.id);
});
$('.detable-helper-wrapper').hide(0);
$('#'+tids[0]).show();

// The only way to align columns...
var timer = setInterval(function() {
    if ($('#'+tids[0]+"_").find("table").length > 0) {
        if ($('#'+tids[0]+"_").find("table")[1].id !== undefined) {
            var id = $('#'+tids[0]+"_").find("table")[1].id;
            $('#'+id).DataTable().columns.adjust();
            clearInterval(timer);
        }
    }
},100);
```

```{r export_links}
if (exportCountsTable) {
    cat("<br/><div class=\"figure-hint link-box\">")
    if (file.exists(file.path(PROJECT_PATH[["lists"]],
        "raw_counts_table.txt.gz")))
        cat("<strong><a href=\"lists/raw_counts_table.txt.gz\" download>",
            "Download</a> the raw read counts table for the experiment.",
            "</strong><br/>",sep="")
    cat("<strong><a href=\"lists/normalized_counts_table.txt.gz\"",
        " download>Download</a> the normalized read counts table for the ",
        "experiment.</strong><br/>",sep="")
    cat("</div>")
}

if (length(qcPlots) > 0) {
    ltxt <- paste(paste("<a href=\"plots/metaseqr_figures_",figFormat,
        ".zip\" download>",figFormat,"</a>",sep=""),collapse=", ")
    cat("<div class=\"figure-hint link-box\">")
    cat("<strong>Get all the figures in ",ltxt," format.</strong>")
    cat("</div>")
}
```

# References

```{r references}
refs <- unique(c(
    reportMessages$references$main,
    reportMessages$references$filein[[fileType]],
    reportMessages$references$norm[[normalization]],
    unlist(reportMessages$references$stat[statistics],use.names=FALSE),
    unlist(reportMessages$references$figure[qcPlots],use.names=FALSE),
    reportMessages$references$multiple[[adjustMethod]],
    reportMessages$references$meta[[metaP]]
))
cat("<div class=\"figure-hint\" style=\"margin-left:25px;\"><ol>")
for (r in refs)
    cat("<li>",r,"</li>",sep="")
cat("</ol></div>")
```

```{js bind_event_changes}
$("#biodetection_selector").on('change',function() {
    var name = this.value;
    getPlotFromDb({name: name,type: "biodetection",subtype: "generic"},
        "biodetection_container");
});

$("#countsbio_sample_selector").on('change',function() {
    var name = this.value;
    getPlotWithFunFromDb({name: name,type: "countsbio",subtype: "sample"},
        "countsbio_sample_container");
});

$("#countsbio_biotype_selector").on('change',function() {
    var name = this.value;
    getPlotWithFunFromDb({name: name,type: "countsbio",subtype: "biotype"},
        "countsbio_biotype_container");
});

$("#saturation_sample_selector").on('change',function() {
    var name = this.value;
    getPlotWithFunFromDb({name: name,type: "saturation",subtype: "sample"},
        "saturation_sample_container");
});

$("#saturation_biotype_selector").on('change',function() {
    var name = this.value;
    getPlotWithFunFromDb({name: name,type: "saturation",subtype: "biotype"},
        "saturation_biotype_container");
});

$("#pairwise_selector").on('change',function() {
    var name = this.value;
    getPlotFromDb({name: name,type: "pairwise",subtype: "xy"},
        "pairwise_xy_container");
    getPlotFromDb({name: name,type: "pairwise",subtype: "md"},
        "pairwise_md_container");
});

$("#meandiff_selector").on('change',function() {
    var name = this.value;
    getPlotFromDb({name: name,type: "meandiff",subtype: "unorm"},
        "meandiff_unorm_container");
    getPlotFromDb({name: name,type: "meandiff",subtype: "norm"},
        "meandiff_norm_container");
});

$("#volcano_selector").on('change',function() {
    var name = this.value;
    getPlotFromDb({name: name,type: "volcano",subtype: "generic"},
        "volcano_container");
});

$("#mastat_selector").on('change',function() {
    var name = this.value;
    getPlotFromDb({name: name,type: "mastat",subtype: "generic"},
        "mastat_container");
});

$("#deregulo_selector").on('change',function() {
    var name = this.value;
    getPlotFromDb({name: name,type: "deregulogram",subtype: "generic"},
        "deregulo_container");
});

$("#heatmap_contrast_selector").on('change',function() {
    var name = this.value;
    var suffix = $("#heatmap_scale_selector").val();
    var newName = name + "_" + suffix;
    $(".deheatmap-container").hide(0);
    $("#"+newName).show();
});

$("#heatmap_scale_selector").on('change',function() {
    var name = this.value;
    var prefix = $("#heatmap_contrast_selector").val();
    var newName = prefix + "_" + name;
    $(".deheatmap-container").hide(0);
    $("#"+newName).show();
});

$("#biodist_selector").on('change',function() {
    var name = this.value;
    getPlotFromDb({name: name,type: "biodist",subtype: "chromosome"},
        "biodist_chromosome_container");
    getPlotFromDb({name: name,type: "biodist",subtype: "biotype"},
        "biodist_biotype_container");
});

$("#statvenn_contrast_selector").on('change',function() {
    var name = this.value;
    var suffix = $("#statvenn_direction_selector").val();
    var newName = name + "_" + suffix;
    var newSubtype = "stat_" + suffix;
    getVennFromDb({name: newName,type: "statvenn",subtype: newSubtype},
        "jvenn_stat_container");
});

$("#statvenn_direction_selector").on('change',function() {
    var name = this.value;
    var prefix = $("#statvenn_contrast_selector").val();
    var newName = prefix + "_" + name;
    var newSubtype = "stat_" + name;
    getVennFromDb({name: newName,type: "statvenn",subtype: newSubtype},
        "jvenn_stat_container");
});

$("#foldvenn_algo_selector").on('change',function() {
    var name = this.value;
    var suffix = $("#foldvenn_direction_selector").val();
    var newName = name + "_" + suffix;
    var newSubtype = "fold_" + suffix;
    getVennFromDb({name: newName,type: "foldvenn",subtype: newSubtype},
        "jvenn_fold_container");
});

$("#foldvenn_direction_selector").on('change',function() {
    var name = this.value;
    var prefix = $("#foldvenn_algo_selector").val();
    var newName = prefix + "_" + name;
    var newSubtype = "fold_" + name;
    getVennFromDb({name: newName,type: "foldvenn",subtype: newSubtype},
        "jvenn_fold_container");
});

$("#detable_selector").on('change',function() {
    var name = this.value;
    $(".detable-helper-wrapper").hide(0);
    $("#"+name).show();
    
    // Hack to make the table visible (htmlwidget known issue)
    window.dispatchEvent(new Event('resize'));
    // Hack to align the column headers
    var id = $('#'+name+"_").find("table")[1].id;
    $('#'+id).DataTable().columns.adjust();
});

// Hidden buttons
$("#_mds_trigger").on('click',function() {
    getPlotFromDb({name: "MDS",type: "mds",subtype: "generic"},"mds_container");
});
$("#_readnoise_trigger").on('click',function() {
    getPlotFromDb({name: "ReadNoise",type: "readnoise",subtype: "generic"},
        "readsnoise_container");
});
$("#_filtered_trigger").on('click',function() {
    getPlotFromDb({name: "filtered_chromosome",type: "filtered",
        subtype: "chromosome"},"chromosome_distr_container");
    getPlotFromDb({name: "filtered_biotype",type: "filtered",
        subtype: "biotype"},"biotype_distr_container");
});
$("#_boxplot_trigger").on('click',function() {
    getPlotFromDb({name: "Boxplot",type: "boxplot",subtype: "unorm"},
        "boxplot_unorm_container");
    getPlotFromDb({name: "Boxplot",type: "boxplot",subtype: "norm"},
        "boxplot_norm_container");
});
$("#_gcbias_trigger").on('click',function() {
    getPlotFromDb({name: "GCBias",type: "gcbias",subtype: "unorm"},
        "gcbias_unorm_container");
    getPlotFromDb({name: "GCBias",type: "gcbias",subtype: "norm"},
        "gcbias_norm_container");
});
$("#_lengthbias_trigger").on('click',function() {
    getPlotFromDb({name: "LengthBias",type: "lengthbias",subtype: "unorm"},
        "lengthbias_unorm_container");
    getPlotFromDb({name: "LengthBias",type: "lengthbias",subtype: "norm"},
        "lengthbias_norm_container");
});
$("#_meanvar_trigger").on('click',function() {
    getPlotFromDb({name: "MeanVar",type: "meanvar",subtype: "unorm"},
        "meanvar_unorm_container");
    getPlotFromDb({name: "MeanVar",type: "meanvar",subtype: "norm"},
        "meanvar_norm_container");
});
$("#_rnacomp_trigger").on('click',function() {
    getPlotFromDb({name: "RnaComp",type: "rnacomp",subtype: "unorm"},
        "rnacomp_unorm_container");
    getPlotFromDb({name: "RnaComp",type: "rnacomp",subtype: "norm"},
        "rnacomp_norm_container");
});

// Initialize some by changing plots (harmless for dexie, required for sqlite)
$('#_mds_trigger').click();
$("#_readnoise_trigger").click();
$("#_filtered_trigger").click();
$("#_boxplot_trigger").click();
$("#_gcbias_trigger").click();
$("#_lengthbias_trigger").click();
$("#_meanvar_trigger").click();
$("#_rnacomp_trigger").click();

$('#biodetection_selector').trigger('change');
$('#countsbio_sample_selector').trigger('change');
$('#countsbio_biotype_selector').trigger('change');
$('#saturation_sample_selector').trigger('change');
$('#saturation_biotype_selector').trigger('change');
$('#pairwise_selector').trigger('change');
$('#meandiff_selector').trigger('change');
$("#volcano_selector").trigger('change');
$("#mastat_selector").trigger('change');
$("#deregulo_selector").trigger('change');
$("#biodist_selector").trigger('change');
if ($("#statvenn_contrast_selector").length > 0) {
    $("#statvenn_contrast_selector").trigger('change');
}
if ($("#foldvenn_algo_selector").length > 0) {
    $("#foldvenn_algo_selector").trigger('change');
}
```
---
title: "Building an annotation database for metaseqR2"
author: "Panagiotis Moulos"
date: "`r BiocStyle::doc_date()`"
output: 
  BiocStyle::html_document:
    toc: true
    toc_float: true
vignette: >
  %\VignetteIndexEntry{Building an annotation database for metaseqR2}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Simple, flexible and reusable annotation for metaseqR2 pipeline
================================================================================

When using the first version of metaseqR, one had to either embed the annotation
to the gene/exon/3' UTR counts, or to download and construct on-the-fly the
required annotation when starting from BAM files. Although the counting and gene
model construction results could (anf still can) be saved and re-used with other
analysis parameters changed (e.g. statistical algorithms), one could not easily
add for example new data to an existing dataset without re-running the whole
pipeline and re-downloading annotation. On top of that, many times, the main
Ensembl servers (when using Ensembl annotations) do not respond well to biomaRt 
calls, so the whole pipeline may stall until the servers are back.

Another main issue with the annotation used by metaseqR was that there was no 
straightforward way, provided by metaseqR, to archive and version the annotation
used by a specific analysis and was up to the user to take care of 
reproducibility at this level. Furthermore, there was no straightforward way for
a user to plugin own annotation elements (e.g. in the form of a GTF file) and 
use it in the same manner as standard annotations supported by metaseqR, e.g. 
when analyzing data from not-so-often studied organisms such as insects. 
Plugging-in own annotation was possible but usually a painful procedure, which
has become now very easy.

The annotation database builder for metaseqR2 remedies the above situations. The
`buildAnnotationDatabase` function should be run once with the organisms one
requires to have locally to work with and then that's it! Of course you can
manage your database by adding and removing specific annotations (and you even
can play with an SQLite browser, although not advised, as the database structure
is rather simple). Furthermore, you can use the metaseqR2 annotation database
and management mechanism for any other type of analysis where you require to 
have a simple tab-delimited annotation file, acquired with very little effort.

# Supported organisms

The following organisms (essentially genome versions) are supported for 
automatic database builds:

* Human (*Homo sapiens*) genome version **hg38**
* Human (*Homo sapiens*) genome version **hg19**
* Human (*Homo sapiens*) genome version **hg18**
* Mouse (*Mus musculus*) genome version **mm10**
* Mouse (*Mus musculus*) genome version **mm9**
* Rat (*Rattus norvegicus*) genome version **rn6**
* Rat (*Rattus norvegicus*) genome version **rn5**
* Fruitfly (*Drosophila melanogaster*) genome version **dm6**
* Fruitfly (*Drosophila melanogaster*) genome version **dm3**
* Zebrafish (*Danio rerio*) genome version **danRer7**
* Zebrafish (*Danio rerio*) genome version **danRer10**
* Zebrafish (*Danio rerio*) genome version **danRer11**
* Chimpanzee (*Pan troglodytes*) genome version **panTro4**
* Chimpanzee (*Pan troglodytes*) genome version **panTro5**
* Pig (*Sus scrofa*) genome version **susScr3**
* Pig (*Sus scrofa*) genome version **susScr11**
* Horse (*Equus cabalus*) genome version **equCab2**
* Arabidopsis (*Arabidobsis thaliana*) genome version **TAIR10**

# Using the local database

## Installation of metaseqR2

To install the metaseqR2 package, start R and enter:

```{r install-0, eval=FALSE, echo=TRUE}
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("metaseqR2")
```

## Setup the database

By default, the database file will be written in the
`system.file(package="metaseqR2")` directory. You can specify another prefered
destination for it using the `db` argument in the function call, but if you do 
that, you will have to supply the `localDb` argument pointing to the SQLite 
database file you created to every metaseqr2 call you perform, otherwise, the
pipeline will download and use annotations on-the-fly.

In this vignette, we will build a minimal database comprising only the mouse
*mm10* genome version from Ensembl. The database will be build in a temporary
directory inside session `tempdir()`.

**Important note**: As the annotation build function makes use of 
[Kent](http://hgdownload.soe.ucsc.edu/admin/exe/) utilities for creating 3'UTR
annotations from RefSeq and UCSC, the latter cannot be built in Windows. 
Therefore it is advised to either build the annotation database in a Linux 
system or use our pre-built databases.

```{r load-0, eval=TRUE, echo=FALSE, tidy=FALSE, message=FALSE, warning=FALSE}
library(metaseqR2)
```

```{r example-1, eval=TRUE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
library(metaseqR2)

buildDir <- file.path(tempdir(),"test_anndb")
dir.create(buildDir)

# The location of the custom database
myDb <- file.path(buildDir,"testann.sqlite")

# Since we are using Ensembl, we can also ask for a version
organisms <- list(mm10=100)
sources <- ifelse(.Platform$OS.type=="unix",c("ensembl","refseq"),"ensembl")

# If the example is not running in a multicore system, rc is ignored
buildAnnotationDatabase(organisms,sources,forceDownload=FALSE,db=myDb,rc=0.5)
```

## Use the database

Now, that a small database is in place, let's retrieve some data. Remember that
since the built database is not in the default location, we need to pass the
database file in each data retrieval function. The annotation is retrieved as
a `GRanges` object by default.

```{r example-2, eval=TRUE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
# Load standard annotation based on gene body coordinates
genes <- loadAnnotation(genome="mm10",refdb="ensembl",level="gene",type="gene",
    db=myDb)
genes

# Load standard annotation based on 3' UTR coordinates
utrs <- loadAnnotation(genome="mm10",refdb="ensembl",level="gene",type="utr",
    db=myDb)
utrs

# Load summarized exon annotation based used with RNA-Seq analysis
sumEx <- loadAnnotation(genome="mm10",refdb="ensembl",level="gene",type="exon",
    summarized=TRUE,db=myDb)
sumEx

# Load standard annotation based on gene body coordinates from RefSeq
if (.Platform$OS.type=="unix") {
    refGenes <- loadAnnotation(genome="mm10",refdb="refseq",level="gene",
        type="gene",db=myDb)
    refGenes
}
```

Or as a data frame if you prefer using `asdf=TRUE`. The data frame however does 
not contain metadata like `Seqinfo` to be used for any susequent validations:

```{r example-3, eval=TRUE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
# Load standard annotation based on gene body coordinates
genes <- loadAnnotation(genome="mm10",refdb="ensembl",level="gene",type="gene",
    db=myDb,asdf=TRUE)
head(genes)
```

## Add a custom annotation

Apart from the supported organisms and databases, you can add a custom 
annotation. Such an annotation can be: 

* A non-supported organism (e.g. an insect or another mammal e.g. dog)
* A modification or further curation you have done to existing/supported
annotations
* A supported organism but from a different source
* Any other case where the provided annotations are not adequate

This can be achieved through the usage of
[GTF](https://www.ensembl.org/info/website/upload/gff.html) files, along with
some simple metadata that you have to provide for proper import to the
annotation database. This can be achieved through the usage of the
`buildCustomAnnotation` function. Details on required metadata can be found
in the function's help page.

**Important note:** Please note that importing a custom genome annotation 
directly from UCSC (UCSC SQL database dumps) is not supported in Windows as the
process involves using the `genePredToGtf` which is not available for Windows.

Let's try a couple of exammples. The first one is a custom annotation for the
Ebola virus from UCSC:

```{r example-4, eval=TRUE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
# Setup a temporary directory to download files etc.
customDir <- file.path(tempdir(),"test_custom")
dir.create(customDir)

# Convert from GenePred to GTF - Unix/Linux only!
if (.Platform$OS.type == "unix" && !grepl("^darwin",R.version$os)) {
    # Download data from UCSC
    goldenPath="http://hgdownload.cse.ucsc.edu/goldenPath/"
    # Gene annotation dump
    download.file(paste0(goldenPath,"eboVir3/database/ncbiGene.txt.gz"),
        file.path(customDir,"eboVir3_ncbiGene.txt.gz"))
    # Chromosome information
    download.file(paste0(goldenPath,"eboVir3/database/chromInfo.txt.gz"),
        file.path(customDir,"eboVir3_chromInfo.txt.gz"))

    # Prepare the build
    chromInfo <- read.delim(file.path(customDir,"eboVir3_chromInfo.txt.gz"),
        header=FALSE)
    chromInfo <- chromInfo[,1:2]
    rownames(chromInfo) <- as.character(chromInfo[,1])
    chromInfo <- chromInfo[,2,drop=FALSE]
    
    # Coversion from genePred to GTF
    genePredToGtf <- file.path(customDir,"genePredToGtf")
    if (!file.exists(genePredToGtf)) {
        download.file(
        "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf",
            genePredToGtf
        )
        system(paste("chmod 775",genePredToGtf))
    }
    gtfFile <- file.path(customDir,"eboVir3.gtf")
    tmpName <- file.path(customDir,paste(format(Sys.time(),"%Y%m%d%H%M%S"),
        "tgtf",sep="."))
    command <- paste0(
        "zcat ",file.path(customDir,"eboVir3_ncbiGene.txt.gz"),
        " | ","cut -f2- | ",genePredToGtf," file stdin ",tmpName,
        " -source=eboVir3"," -utr && grep -vP '\t\\.\t\\.\t' ",tmpName," > ",
        gtfFile
    )
    system(command)

    # Build with the metadata list filled (you can also provide a version)
    buildCustomAnnotation(
        gtfFile=gtfFile,
        metadata=list(
            organism="eboVir3_test",
            source="ucsc_test",
            chromInfo=chromInfo
        ),
        db=myDb
    )

    # Try to retrieve some data
    eboGenes <- loadAnnotation(genome="eboVir3_test",refdb="ucsc_test",
        level="gene",type="gene",db=myDb)
    eboGenes
}
```

Another example, the Atlantic cod from UCSC. The same things apply for the
operating system.

```{r example-5, eval=FALSE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
if (.Platform$OS.type == "unix") {
    # Gene annotation dump
    download.file(paste0(goldenPath,"gadMor1/database/augustusGene.txt.gz"),
        file.path(customDir,"gadMori1_augustusGene.txt.gz"))
    # Chromosome information
    download.file(paste(goldenPath,"gadMor1/database/chromInfo.txt.gz",sep=""),
        file.path(customDir,"gadMori1_chromInfo.txt.gz"))

    # Prepare the build
    chromInfo <- read.delim(file.path(customDir,"gadMori1_chromInfo.txt.gz"),
        header=FALSE)
    chromInfo <- chromInfo[,1:2]
    rownames(chromInfo) <- as.character(chromInfo[,1])
    chromInfo <- chromInfo[,2,drop=FALSE]
    
    # Coversion from genePred to GTF
    genePredToGtf <- file.path(customDir,"genePredToGtf")
    if (!file.exists(genePredToGtf)) {
        download.file(
        "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf",
            genePredToGtf
        )
        system(paste("chmod 775",genePredToGtf))
    }
    gtfFile <- file.path(customDir,"gadMori1.gtf")
    tmpName <- file.path(customDir,paste(format(Sys.time(),"%Y%m%d%H%M%S"),
        "tgtf",sep="."))
    command <- paste0(
        "zcat ",file.path(customDir,"gadMori1_augustusGene.txt.gz"),
        " | ","cut -f2- | ",genePredToGtf," file stdin ",tmpName,
        " -source=gadMori1"," -utr && grep -vP '\t\\.\t\\.\t' ",tmpName," > ",
        gtfFile
    )
    system(command)

    # Build with the metadata list filled (you can also provide a version)
    buildCustomAnnotation(
        gtfFile=gtfFile,
        metadata=list(
            organism="gadMor1_test",
            source="ucsc_test",
            chromInfo=chromInfo
        ),
        db=myDb
    )

    # Try to retrieve some data
    gadGenes <- loadAnnotation(genome="gadMor1_test",refdb="ucsc_test",
        level="gene",type="gene",db=myDb)
    gadGenes
}
```

Another example, Armadillo from Ensembl. This should work irrespectively of 
operating system. We are downloading chromosomal information from UCSC.

```{r example-6, eval=FALSE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
# Gene annotation dump from Ensembl
download.file(paste0("ftp://ftp.ensembl.org/pub/release-98/gtf/",
    "dasypus_novemcinctus/Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
    file.path(customDir,"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"))

# Chromosome information will be provided from the following BAM file
# available from Ensembl. We have noticed that when using Windows as the OS,
# a remote BAM files cannot be opened by scanBamParam, so for this example,
# chromosome length information will not be available when running in Windows.
bamForInfo <- NULL
if (.Platform$OS.type == "unix")
    bamForInfo <- paste0("ftp://ftp.ensembl.org/pub/release-98/bamcov/",
        "dasypus_novemcinctus/genebuild/Dasnov3.broad.Ascending_Colon_5.1.bam")

# Build with the metadata list filled (you can also provide a version)
buildCustomAnnotation(
    gtfFile=file.path(customDir,"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
    metadata=list(
        organism="dasNov3_test",
        source="ensembl_test",
        chromInfo=bamForInfo
    ),
    db=myDb
)

# Try to retrieve some data
dasGenes <- loadAnnotation(genome="dasNov3_test",refdb="ensembl_test",
    level="gene",type="gene",db=myDb)
dasGenes
```

## A complete build

A quite complete build (with latest versions of Ensembl annotations) would look
like (supposing the default annotation database location):

```{r example-7, eval=FALSE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
organisms <- list(
    hg18=67,
    hg19=75,
    hg38=97:98,
    mm9=67,
    mm10=97:98,
    rn5=79,
    rn6=97:98,
    dm3=78,
    dm6=97:98,
    danrer7=79,
    danrer10=91,
    danrer11=97:98,
    pantro4=90,
    pantro5=97:98,
    susscr3=89,
    susscr11=97:98,
    equcab2=97:98
)

sources <- c("ensembl","ucsc","refseq")

buildAnnotationDatabase(organisms,sources,forceDownload=FALSE,rc=0.5)
```

The aforementioned complete built can be found
[here](https://tinyurl.com/ybycpr6b)
Complete builts will become available from time to time (e.g. with every new
Ensembl relrase) for users who do not wish to create annotation databases on
their own. Root access may be required (depending on the metaseqR2 library
location) to place it in the default location where it can be found 
automatically.

# Annotations on-the-fly

If for some reason you do not want to build and use an annotation database for
metaseqR2 analyses (not recommended) or you wish to perform an analysis with an
organism that does not yet exist in the database, the `loadAnnotation` function
will perform all required actions (download and create a `GRanges` object) 
on-the-fly as long as there is an internet connection.

However, the above function does not handle custom annotations in GTF files.
In a scenario where you want to use a custom annotation only once, you should
supply the `annotation` argument to the `metaseqr2` function, which is almost
the same as the `metadata` argument used in `buildCustomAnnotation`, actually 
augmented by a list member for the GTF
file, that is:

```{r pseudo-1, eval=TRUE, echo=TRUE, message=TRUE, warning=FALSE}
annotation <- list(
    gtf="PATH_TO_GTF",
    organism="ORGANISM_NAME",
    source="SOURCE_NAME",
    chromInfo="CHROM_INFO"
)
```

The above argument can be passed to the metaseqr2 call in the respective 
position.

For further details about custom annotations on the fly, please check
`buildCustomAnnotation` and `importCustomAnnotation` functions.

# Session Info

```{r si-1, eval=TRUE, echo=TRUE}
sessionInfo()
```
---
title: "Usage of the metaseqR2 package"
author: "Panagiotis Moulos"
date: "`r BiocStyle::doc_date()`"
output: 
  BiocStyle::html_document:
    toc: true
    toc_float: true 
vignette: >
  %\VignetteIndexEntry{RNA-Seq data analysis with metaseqR2}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

RNA-Seq data analysis using mulitple statistical algorithms with metaseqR2
================================================================================

During the past years, a lot of packages have been developed for the analysis of
RNA-Seq data, introducing several approaches. Many of them live in Bioconductor.
Furthermore, different statistical approaches and heuristics are used in a
continuous effort to improve overall accuracy. Such approaches include packages
using the negative binomial distribution to model the null hypotheses (*DESeq*, 
*DESeq2*, *edgeR*, *NBPSeq*, "ABSSeq"), packages using Bayesian statistics 
(*baySeq*, *EBSeq*) or more hybrid solutions (*NOISeq*, *voom*). In addition, 
packages specialized to RNA-Seq data normalization have also been developed 
(*EDASeq*, *RUVSeq*). The first version of the *metaseqR* package (pronounced 
meta-seek-er) provided an interface to several algorithms for normalization and 
statistical analysis and at the same time provided PANDORA, a novel p-value 
combination method. PANDORA successfully combines several statistical algorithms
by weighting their outcomes according to their performance with realistically 
simulated data sets generated from real data. Using simulated as well as real 
data, it was shown that PANDORA improves the overall detection of differentially
expressed genes by reducing false hits while maintaining true positives. To our
knowledge, PANDORA remains the only fully functional method proposing this 
combinatorial approach for the analysis of RNA-Seq data.

metaseqR2, is the continuation of
[metaseqR](https://pubmed.ncbi.nlm.nih.gov/25452340/). While it has been 
(at times) heavily refactored, it still offers the same functionalities with as 
much backwards compatibility as possible. Like metaseqR, metaseqR2, incoporates 
several algorithms for normalization and statistical analysis. In particular, we
extended the offered algorithms with *DESeq2*, *ABSSeq* and *DSS*. metaseqR2, 
like metaseqR also builds a full report with several interactive and 
non-interactive diagnostic plots so that the users can easily explore the 
results and have whatever they need for this part of their research in one 
place. The report has been modernized and remains one of its strongest points as
it provides an automatically generated summary, based on the pipeline inputs and
the results, which can be used directly as a draft in methods paragraph in 
scientific publications. It also provides a lot of diagnostic figures and each 
figure is accompanied by a small explanatory text, and a list of references 
according to the algorithms used in the pipeline. metaseqR2 continues to provide
an interface for RNA-Seq data meta-analysis by providing the ability to use 
different algorithms for the statistical testing part and combining the p-values
using popular published methods (e.g. Fisher's method, Whitlock's method), 
two package-specific methods (intersection, union of statistically significant 
results) and of course PANDORA.

Another major difference as compared to the older metaseqR package is the 
annotation system that is adopted by metaseqR2. More specifically, metaseqR2
introduces the `buildAnnotationDatabase` function which builds a local SQLite
database with the supported by metaseqR annotations as well as additional
versions added in the current package. This function, given a short and
comprehensive number of arguments, automatically downloads, processes and 
imports to a portable database, all annotation types required by the
main analysis pipeline. Therefore, the user neither has to embed nor download
the required annotation each time. But most importantly, with the current 
package, the user is able also to provide an own GTF file with custom annotation
elements that are the imported to the metaseqR2 database and annotation system
and can be used for the respective analyses.

Apart from local database building, there also other major additions (such) as
improved analysis for 3'UTR mRNA sequencing (Lexogen Quant-Seq protocol) which 
can be found towards the end of this page.

Throughout the rest of this document, `metaseqr2` refers to the name of the  
analysis pipeline while *metaseqR2* refers to the name of the package.

# Getting started

## Installation

To install the metaseqR2 package, start R and enter:

```{r install-0, eval=FALSE, echo=TRUE}
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("metaseqR2")
```

## Introduction

```{r load-library, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}
library(metaseqR2)
```

Detailed instructions on how to run the metaseqr2 pipeline can be found under 
the main documentation of the metaseqR2 package.

Briefly, to run metaseqr2 you need:

1. Input RNA-Seq data. These can come in three forms:
  * A text tab delimited file in a spreadsheet-like format containing at least
    unique gene identifiers (corresponding to one of metaseqR2 supported 
    annotation sources, that is Ensembl, UCSC, RefSeq) *or* if you are using
    a custom annotation (with a GTF file), unique gene identifiers corresponding
    to this GTF file. This case is applicable in case of receiving a ready-made
    counts table from an external source, such as a sequencing facility or a
    public dataset.
  * A text tab delimited file in a spreadsheet-like format containing all the
    required annotation elements and additional columns with read counts. This
    solution is applicable only for gene analysis (`transLevel = "gene"` and 
    `countType = "gene"`). Generally, it is not recommended to embed the
    annotation and this case is supported only for backwards compatibility.
  * A set of BAM files, aligned according to the mRNA sequencing protocol,
    usually a spliced aligner like HiSat or STAR. This is the recommended 
    analysis procedure and the BAM files are declared in a targets text file.
2. A local annotation database. This is not required as all required annotation 
   can be downloaded on the fly, but it is recommended for speed, if you have a 
   lot of analyses to perform.
3. A list of statistical contrasts for which you wish to check differential
   expression
4. An internet connection so that the interactive report can be properly 
   rendered, as the required JavaScript libraries are not embedded to the
   package. This is required only once as the report is then self-contained.

For demonstration purposes, a very small dataset (with embedded annotation) is
included with the package.

## Types of analyses performed with metaseqR2

Several types of differential analysis of gene expression can be performed and
reported with metaseqR2 depending on the biological question asked and the type
of data generated. For example, an investigator may be interested in gene- or
transcript-level differential expression analysis when a 3'UTR sequencing kit
has been used or interested for differential exon usage when a classical
polyA RNA-Seq protocol has been applied.

These analysis types are being defined essentially by two arguments:
* `countType` which can be `gene`, `exon`, `utr` corresponding to total RNA
sequencing, polyA RNA sequencing or 3' UTR sequencing respectively.
* `transLevel` which can be `gene`, `transcript`, `exon` corresponding to
differetial expression analysis using gene models (or essentially the dominant
transcripts), individual transcripts or exons respectively.

Therefore, the selection of `countType="exon"` and `transLevel="gene"` assumes
that we have a dataset where polyA RNA sequencing has been applied followed by
splicing-aware alignment while `countType="utr"` and `transLevel="transcript"`
assumes that we have a dataset where 3'UTR sequencing (e.g. Lexogen Quant-Seq) 
has been applied to look for differential expression based on read occupancy on
the 3' UTR regions.

The following combinations are available:
* `countType="gene"`, `transLevel="gene"` for differential expression analysis
using a pre-calculated counts table or BAM files from total RNA sequencing.
* `countType="gene"`, `transLevel="transcript"` for differential expression 
analysis using a pre-calculated counts table or BAM files from total RNA 
sequencing and for each transcript.
* `countType="gene"`, `transLevel="exon"` for differential expression analysis
of exons using BAM files from polyA RNA sequencing.
* `countType="exon"`, `transLevel="gene"` for differential expression analysis
using BAM files from polyA RNA sequencing.
* `countType="exon"`, `transLevel="transcript"` for differential expression 
analysis of transcripts using BAM files from total RNA sequencing.
* `countType="utr"`, `transLevel="gene"` for differential expression analysis
of genes using BAM files from 3' UTR RNA sequencing.
* `countType="utr"`, `transLevel="transcript"` for differential expression 
analysis of transcripts using BAM files from 3' UTR RNA sequencing.

## Data filtering

The metaseqR2 pipeline has several options for gene filtering at the gene and 
exon levels. These filters span various areas including:
* The presence of a minimum number of reads in a fraction of the samples per 
condition or experiment-wise.
* The exclusion of specific biotypes (e.g. exluding pseudogenes)
* The filtering based on several expression attributes such as average read
presence over *n* kbs or the exclusion of genes whose expression is below the
expression of a set of genes known *not* to be expressed in the biological
mechanism under investigation
* Filters based on exon expression such as the minimum fraction of exons that
should contain reads over a gene.

In addition, the metaseqR2 pipeline offers several analysis "presets" with
respect to the filtering layers applied, the statistical analysis stringency and
the amount of data exported.

All the aforementioned parameters are well-documented in the main manual of the
package and the respective man pages.

# Running the `metaseqr2` pipeline

**Note**: When conducting an analysis with metaseqR2, it is advised that you
set a seed for random number generation using `set.seed()`. This should be set
because some quality control charts in the metaseqR2 report are created by
downsampling the initial dataset analyzed. Therefore, to guarantee the
reproducibility of these plots, a seed must be provided.

```{r seed-0, eval=TRUE, echo=TRUE}
set.seed(42)
```

Running a metaseqr2 pipeline instance is quite straightforward. Again, see the
examples in the main help page. Below, an example and the command window output 
follow:

```{r data-1, eval=TRUE, echo=TRUE}
data("mm9GeneData",package="metaseqR2")
```

```{r head-1, eval=TRUE, echo=TRUE}
head(mm9GeneCounts)
```

```{r random-1, eval=TRUE, echo=TRUE}
sampleListMm9
```

```{r random-2, eval=TRUE, echo=TRUE}
libsizeListMm9
```

## Analysis at the gene level with gene counts

Following, a full example with the informative messages that are printed in the
command window:

```{r example-1, eval=TRUE, echo=TRUE, tidy=FALSE, message=TRUE, warning=FALSE}
library(metaseqR2)

data("mm9GeneData",package="metaseqR2")

# You can explore the results in the session's temporary directory
print(tempdir())

result <- metaseqr2(
    counts=mm9GeneCounts,
    sampleList=sampleListMm9,
    contrast=c("adult_8_weeks_vs_e14.5"),
    libsizeList=libsizeListMm9,
    annotation="embedded",
    embedCols=list(
        idCol=4,
        gcCol=5,
        nameCol=8,
        btCol=7
    ),
    org="mm9",
    countType="gene",
    normalization="edger",
    statistics="edger",
    pcut=0.05,
    qcPlots=c(
        "mds","filtered","correl","pairwise","boxplot","gcbias",
        "lengthbias","meandiff","meanvar","deheatmap","volcano",
        "mastat"
    ),
    figFormat=c("png","pdf"),
    exportWhat=c("annotation","p_value","adj_p_value","fold_change"),
    exportScale=c("natural","log2"),
    exportValues="normalized",
    exportStats=c("mean","sd","cv"),
    exportWhere=file.path(tempdir(),"test1"),
    restrictCores=0.01,
    geneFilters=list(
         length=list(
                length=500
         ),
         avgReads=list(
                averagePerBp=100,
                quantile=0.25
         ),
         expression=list(
                median=TRUE,
                mean=FALSE,
                quantile=NA,
                known=NA,
                custom=NA
         ),
         biotype=getDefaults("biotypeFilter","mm9")
    ),
    outList=TRUE
)
```

To get a glimpse on the results, run:

```{r head-2, eval=TRUE, echo=TRUE}
head(result[["data"]][["adult_8_weeks_vs_e14.5"]])
```

You may also want to check the interactive HTML report generated in the output 
directory defined by the `exportWhere` argument above.

Now, the same example but with more than one statistical selection algorithms, a
different normalization, an analysis preset and filtering applied prior to
normalization:

```{r example-2, eval=TRUE, echo=TRUE, tidy=FALSE, message=FALSE, warning=FALSE}
library(metaseqR2)

data("mm9GeneData",package="metaseqR2")

result <- metaseqr2(
    counts=mm9GeneCounts,
    sampleList=sampleListMm9,
    contrast=c("adult_8_weeks_vs_e14.5"),
    libsizeList=libsizeListMm9,
    annotation="embedded",
    embedCols=list(
        idCol=4,
        gcCol=5,
        nameCol=8,
        btCol=7
    ),
    org="mm9",
    countType="gene",
    whenApplyFilter="prenorm",
    normalization="edaseq",
    statistics=c("deseq","edger"),
    metaP="fisher",
    #qcPlots=c(
    #    "mds","biodetection","countsbio","saturation","readnoise","filtered",
    #    "correl","pairwise","boxplot","gcbias","lengthbias","meandiff",
    #    "meanvar","rnacomp","deheatmap","volcano","mastat","biodist","statvenn"
    #),
    qcPlots=c(
        "mds","filtered","correl","pairwise","boxplot","gcbias",
        "lengthbias","meandiff","meanvar","deheatmap","volcano",
        "mastat"
    ),
    restrictCores=0.01,
    figFormat=c("png","pdf"),
    preset="medium_normal",
    exportWhere=file.path(tempdir(),"test2"),
    outList=TRUE
)
```

A similar example with no filtering applied and no Venn diagram generation:

```{r example-3, eval=TRUE, echo=TRUE, tidy=FALSE, message=FALSE, warning=FALSE}
library(metaseqR2)

data("mm9GeneData",package="metaseqR2")

result <- metaseqr2(
    counts=mm9GeneCounts,
    sampleList=sampleListMm9,
    contrast=c("adult_8_weeks_vs_e14.5"),
    libsizeList=libsizeListMm9,
    annotation="embedded",
    embedCols=list(
        idCol=4,
        gcCol=5,
        nameCol=8,
        btCol=7
    ),
    org="mm9",
    countType="gene",
    normalization="edaseq",
    statistics=c("deseq","edger"),
    metaP="fisher",
    qcPlots=c(
        "mds","filtered","correl","pairwise","boxplot","gcbias",
        "lengthbias","meandiff","meanvar","deheatmap","volcano",
        "mastat"
    ),
    restrictCores=0.01,
    figFormat=c("png","pdf"),
    preset="medium_normal",
    outList=TRUE,
    exportWhere=file.path(tempdir(),"test3")
)
```

Another example with the full PANDORA algorithm (not evaluated here):

```{r example-4, eval=TRUE, echo=TRUE, tidy=FALSE, message=FALSE, warning=FALSE}
library(metaseqR2)

data("mm9GeneData",package="metaseqR2")

result <- metaseqr2(
    counts=mm9GeneCounts,
    sampleList=sampleListMm9,
    contrast=c("adult_8_weeks_vs_e14.5"),
    libsizeList=libsizeListMm9,
    annotation="embedded",
    embedCols=list(
        idCol=4,
        gcCol=5,
        nameCol=8,
        btCol=7
    ),
    org="mm9",
    countType="gene",
    normalization="edaseq",
    statistics=c("edger","limma"),
    metaP="fisher",
    figFormat="png",
    preset="medium_basic",
    qcPlots=c(
        "mds","filtered","correl","pairwise","boxplot","gcbias",
        "lengthbias","meandiff","meanvar","deheatmap","volcano",
        "mastat"
    ),
    restrictCores=0.01,
    outList=TRUE,
    exportWhere=file.path(tempdir(),"test4")
)
```

## Analysis at the gene level with exon counts

**Note**: Be sure to have constructed a metaseqR2 annotation database prior to
continuing with the following examples!

As example BAM files from a realistic dataset that can demonstrate the full 
availabilities of metaseqR2 do not fit within the Bioconductor package, you can 
find additional examples in our GitHub 
[page](https://github.com/pmoulos/metaseqR2) where issues can be reported too.

## Estimating p-value weights

In metaseqR2, the PANDORA algorithm is expaned with additional 3 algorithms.
Briefly, PANDORA use of the area under False Discovery Curves to assess the 
performance of each statistical test with simulated datasets created from true 
datasets (e.g. your own dataset, as long as it has a sufficient number of
replicates). Then, the performance assessment can be used to construct p-value 
weights for each test and use these weights to supply the p-value weights 
parameter of metaseqr2 when `metaP` is `"weight"`, `"pandora"` or `"whitlock"` 
(see the next sections for p-value combination methods). The following example 
shows how to create such weights (depending on the size of the dataset, it might
take some time to run):

```{r example-5, eval=TRUE, echo=TRUE, tidy=FALSE, message=FALSE, warning=FALSE}
data("mm9GeneData",package="metaseqR2")
weights <- estimateAufcWeights(
    counts=as.matrix(mm9GeneCounts[,9:12]),
    normalization="edaseq",
    statistics=c("edger","limma"),
    nsim=1,N=10,ndeg=c(2,2),top=4,modelOrg="mm10",
    rc=0.01,libsizeGt=1e+5
)
```

...and the weights...

```{r head-3, eval=TRUE, echo=TRUE}
weights
```

## Combining p-values from multiple tests

Although the main `metaseqr2` function takes care of p-value combination,
sometimes there is the need of simply importing externally calculated p-values
and using the respective metaseqR2 functions to produce combined p-values. We
demonstrate this capability using p-values from all metaseqR2 supported 
algorithms, applied to data from 
[Giakountis et al., 2016](https://doi.org/10.1016/j.celrep.2016.05.038).

```{r example-6, eval=TRUE, echo=TRUE, tidy=FALSE, message=FALSE, warning=FALSE}
data("hg19pvalues",package="metaseqR2")

# Examine the data
head(hg19pvalues)

# Now combine the p-values using the Simes method
pSimes <- apply(hg19pvalues,1,combineSimes)

# The harmonic mean method with PANDORA weights
w <- getWeights("human")
pHarm <- apply(hg19pvalues,1,combineHarmonic,w)

# The PANDORA method
pPandora <- apply(hg19pvalues,1,combineWeight,w)
```

# metaseqR2 components

## Brief description

The metaseqR2 package includes several functions which are responsible for 
running each part of the pipeline (data reading and summarization, filtering, 
normalization, statistical analysis and meta-analysis and reporting). Although 
metaseqR2 is designed to run as a pipeline, where all the parameters for each 
individual part can be passed in the main function, several of the individual 
functions can be run separately so that the more experienced user can build 
custom pipelines. All the HTML help pages contain analytical documentation on 
how to run these functions, their inputs and outputs and contain basic examples.
For example, runnning

```{r help-2, eval=TRUE, echo=TRUE, message=FALSE}
help(statEdger)
```

will open the help page of the wrapper function over the edgeR statistical 
testing algorithm which contains an example of data generation, processing, up 
to statistical selection.

Most of the diagnostic plots, work with simple matrices as inputs, so they can 
be easily used outside the main pipeline, as long as all the necessary arguments
are given. In metaseqR2, most of the individual diagnostic plot creation 
functions are not exported, mostly for documentation simplicity and avoidance of
confusion for non-experts. They can still be used by calling them as 
non-exported objects (e.g. `metaseqR2:::diagplotMds`). Finally, it should be
noted that a report can be generated only when running the whole metaseqr2
pipeline and in the current version there is no support for generating custom 
reports. The final reports contains full interactive graphs and the required
JavaScript libraries to generate them are automatically downloaded.

## Backwards compatibility

If you have older pipelines based on metaseqR and the `metaseqr` function where
the argument coding style is different (e.g. `sample.list` instead of 
`sampleList`) then `metaseqr2` will do its best to convert old arguments to new
arguments so that old commands do not break and the only that should be changed 
is `metaseqr` to `metaseqr2`. Note however that you _should not_ mix old and new
arguments. In this case, the new pipeline will fail.

# The report

In the end of each metaseqr2 pipeline run, a detailed HTML report of the
procedure and the findings is produced. Apart from description of the process,
all the input parameters and other data related to the differential expression
analysis, the report contains a lot of interactive graphs. Specifically, all the
quality control and result inspection plots are interactive. This is achieved
by making extensive use of the JavaScript libraries 
[Highcharts](https://www.highcharts.com/), [Plotly](https://plot.ly/) and 
[jvenn](http://jvenn.toulouse.inra.fr/app/index.html) to create more 
user-friendly and directly explorable plots. By default metaseqr2 produces all
available diagnostic plots, according always to input. For example, if the
*biotype* feature is not available in a case where `annotation="embedded"`, 
plots like `biodetection` and `countsbio` will not be available. If not all 
diagnostic plots are not required, a selection can be made with the `qcPlots`
argument, possibly making the report "lighter" and less browser-demanding.

The HTML report creation mechanism is through the packages rmarkdown and knitr. 
This means that the [Pandoc](https://pandoc.org/) libraries must be installed. 
A lot of details on this can be found in Pandoc's website as well as knitr and
rmarkdown websites and guides. Although the generic mechanism is more
computationally demanding than standard HTML (e.g. using *brew* as in the
previous metaseqR), the results are more standardized, cross-platform and fully
reproducible.

During development, we found out that knitr faces certain difficulties in our 
settings, that is embedding a lot of predefined graphs in JSON format and  all 
required libraries and data in a single HTML page. This situation led to crashes
because of memory usage and of course, very large HTML files. We resolved this 
by using (according to usage scenario and where the report is intended to be 
seen):

1. A flavor of [IndexedDB](https://javascript.info/indexeddb) called 
[Dexie](https://dexie.org/)
2. A JavaScript port of SQLite called 
[sql.js](https://github.com/kripken/sql.js/)

Regarding case (1), IndexedDB is a modern technology to create simple,
in-browser object databases which has several usages, but mostly to avoid the
burden of synchronously loading big-sized objects at the same time of simple 
HTML rendering. IndexedDB is supported by all modern browser and is essentially 
a replacement for `localStorage` which had space limitations. Dexie is a simple
interface to IndexedDB. Thus, all the plot data are created and stored in Dexie 
for rendering when needed. This rendering method can be used both when the 
report is seen as a stand-alone document, locally, without the presence of a web
server or internet connection, and is the default method.

Regarding case (2), all the predefined plot data are stored in a report-specific
SQLite database which is then queried using sql.js. This way can be chosen
when it is known that the report will be presented through a web server (e.g. 
Apache) as in any other case, modern web browser (except MS Edge) do not allow
by default opening local files from an HTML page for security reasons. Also,
sql.js is quite large as a library (altough downloaded once for recurring
reports). This method produces slightly smaller files but is slightly slower.
Using Dexie is the preferred and safest method for both scenarios.

In both cases, the serialized JSON used for Highcharts and jvenn plots is placed
in `data/reportdb.js` when using Dexie or `data/reportdb.sqlite` when using 
sql.js. Experienced users can then open these files and tweak the plots as
desired. The above paths are relative to the report's location `exportWhere` 
arguments.

metaseqR2 report has the following sections, depending also on which diagnostic
and exploration plots have been asked from the run command. As plots are 
categorized, if no plot from a specific category is asked, then this category
will not appear. Below, the categories:

## Summary

The Summary section is further categorized in several subsections. Specifically:

* Analysis summary: This section contains an auto-generated text that 
analytically describes the computational process followed and summarized the
results of each step. This text can be used as is or with slight modifications
in a _Methods_ section of an article.
* Input options: This section provides a list of the input arguments to the
pipeline in a more human-readable format.
* Filtering: This section reports in detail the number of filtered genes
decomposed according to the number of genes removed by each applied filter.
* Differential expression: This section reports in detail the number of
differentially expressed genes for each contrast, both when using only a p-value
cutoff as well as an FDR cutoff (numbers in parentheses), that is, genes passing
the multiple testing correction procedure selected. These numbers also are 
calculated based on a simple fold change cutoff in log<sub>2</sub> scale. 
Finally, when multiple algorithms are used with p-value combination, this
section reports all the findings analytically per algorithm.
* Command: This section contains the command used to run the metaseqr2 pipeline 
for users that want to experiment.
* Run log: This section contains critical messages displayed within the R 
session running `metaseqr2` displayed as a log.

## Quality control 

The Quality control section contains several interactive plots concerning the 
overall quality control of each sample provided as well as overall assessments. 
The quality control plots are the Multidimensional Scaling (MDS) plot, the 
Biotypes detection (Biodetection) plot, the Biotype abundance (Countsbio) plot, 
the Read saturation (Saturation) plot, the Read noise (ReadNoise) plot, the 
Correlation heatmap (Correlation), the Pairwise sample scatterplots (Pairwise) 
and the Filtered entities (Filtered) plot. Each plot is accompanied by a 
detailed description of what it depicts. Where multiple plot are available (e.g.
one for each sample), a selection list on the top of the respective section 
allows the selection of the sample to be displayed.

## Normalization

The Normalization section contains several interactive plots that can be used to
inspect and assess the normalization procedure. Therefore, normalization plots 
are usually paired, showing the same data instance normalized and not 
normalized. The normalization plots are the Expression boxplots (Boxplots)
plots, the GC content bias (GC bias) plots, the Gene length bias (Length bias) 
plots, the Within condition mean-difference (Mean-Difference) plots, the 
Mean-variance relationship (Mean-Variance) plot and the RNA composition (Rna 
composition) plot. Each plot is accompanied by a detailed description of what it
depicts. Where multiple plot are available (e.g. one for each sample), a 
selection list on the top of the respective section allows the selection of the 
sample to be displayed.

## Statistics

The Statistics section contains several interactive plots that can be used to 
inspect and explore the outcome of statistical testing procedures. The 
statistics plots are the Volcano plot (Volcano), the MA or Mean-Difference 
across conditions (MA) plot, the Expression heatmap (Heatmap) plot, the 
Chromosome and biotype distributions (Biodist) plot, the Venn diagram across 
statistical tests (StatVenn), the Venn diagram across contrasts (FoldVenn) and
the Deregulogram. Each plot is accompanied by a detailed description of what it 
depicts. Please note that the heatmap plots show only the top percentage of 
differentially expressed genes as this is controlled by the `reportTop` 
parameter of the pipeline. When multiple plots are available (e.g. one for each
contrast), a selection list on the top of the respective section allows the 
selection of the sample to be displayed.

## Results

The Results section contains a snapshot of the differentially expressed genes in
table format with basic information about each gene and some links to external 
resources. Certain columns of the table are colored according to significance. 
Larger bars and more intense colors indicate higher significance. For example, 
bar in the *p_value* column is larger if the genes has higher statistical 
significance and the fold change cell background is bright red if the gene is 
highly up-regulated. From the **Results** section, full gene lists can be 
downloaded in text tab-delimited format and viewed with a spreadsheet
application like MS Excel. A selector on the top of the section above the table
allows the display of different contrasts.

## References

The References section contains bibliographical references regading the
algorihtms used by the metaseqr2 pipeline and is adjusted according to the
algorithms selected.

# Genome browser tracks

metaseqR2 utilizes Bioconductor facilities to create normalized bigWig files.
It also creates a link to open single stranded tracks in the genome browser and
a track hub to display stranded tracks, in case where a stranded RNA-Seq 
protocol has been applied. Just make sure that their output directory is served 
by a web server like Apache. See main documentation for more details.

Please note that if requested, metaseqR2 will try to create tracks even with a
custom organism. This is somewhat risky as

* the track generation may fail
* for heavily customized cases, you will manually have to crate aso .2bit files
for visualization in e.g. the UCSC Genome Browser

Nevertheless, we have chosen to allow the track generation as, many times a user
just uses slight modifications of e.g. the human genome annotation, where some
elements may be manually curated, of elements are added (e.g. non-annotated
non-coding RNAs). Therefore, in case of custom organisms, a warning is thrown
but the functionality is not turned off. Please turn off manually if you are
sure you do not want tracks. You may also use the `createSignalTracks` function
directly.

# List of required packages

Although this is not usually the content of a vignette, the complex nature of
the package requires this list to be populated also here. Therefore, metaseqR2
would benefit from the existence of all the following packages:

 * ABSSeq
 * baySeq
 * Biobase
 * BiocGenerics
 * BiocManager
 * BiocParallel
 * BiocStyle
 * biomaRt
 * Biostrings
 * BSgenome
 * corrplot
 * DESeq
 * DESeq2
 * DSS
 * DT
 * EDASeq
 * edgeR
 * GenomeInfoDb
 * GenomicAlignments
 * GenomicFeatures
 * GenomicRanges
 * gplots
 * graphics
 * grDevices
 * heatmaply
 * htmltools
 * httr
 * IRanges
 * jsonlite
 * knitr
 * limma
 * log4r
 * magrittr
 * methods
 * NBPSeq
 * NOISeq
 * pander
 * parallel
 * qvalue
 * rmarkdown
 * rmdformats
 * RMySQL
 * Rsamtools
 * RSQLite
 * rtracklayer
 * RUnit
 * S4Vectors
 * stats
 * stringr
 * SummarizedExperiment
 * survcomp
 * TCC
 * utils
 * VennDiagram
 * vsn
 * zoo

A recent version of [Pandoc](https://pandoc.org/) is also required, ideally
above 2.0.

# Session Info

```{r si-1, eval=TRUE, echo=TRUE}
sessionInfo()
```
\name{diagplotFiltered}
\alias{diagplotFiltered}
\title{Diagnostic plot for filtered genes}
\usage{
    diagplotFiltered(x, y, output = "x11", path = NULL, ...)
}
\arguments{
    \item{x}{an annotation data frame like the ones produced
    by \code{\link{getAnnotation}}. \code{x} should be the
    filtered annotation according to metaseqR's filters.}

    \item{y}{an annotation data frame like the ones produced
    by \code{\link{getAnnotation}}. \code{y} should contain
    the total annotation without the application of any
    metaseqr filter.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{which.plot} argument. If
    output=\code{"x11"}, no output filenames are produced.
}
\description{
    This function plots a grid of four graphs depicting: in
    the first row, the numbers of filtered genes per
    chromosome in the first column and per biotype in the
    second column. In the second row, the percentages of
    filtered genes per chromosome related to the whole genome
    in the first columns and per biotype in the second
    column.
}
\examples{
data("mm9GeneData",package="metaseqR2")
y <- mm9GeneCounts[,c(1:6,8,7)]
x <- y[-sample(1:nrow(y),1000),]
diagplotFiltered(x,y)
}
\author{
    Panagiotis Moulos
}
\name{combineWeight}
\alias{combineWeight}
\title{Combine p-values using weights}
\usage{
    combineWeight(p, w, zerofix = NULL)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}

    \item{w}{a weights vector, must sum to 1.}
    
    \item{zerofix}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1.}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by
    metaseqR using p-value weights.
}
\details{
    The argument \code{zerofix} is used to correct for
    the case of a p-value which is equal to 0 as a result
    of internal numerical and approximation procedures.
    When \code{NULL}, random numbers greater than 0 and
    less than or equal to 0.5 are used to multiply the
    offending p-values with the lowest provided non-zero
    p-value, maintaining thus a virtual order of 
    significance, avoiding having the same p-values for 
    two tests and assuming that all zero p-values represent
    extreme statistical significance. When a numeric
    between 0 and 1, this number is used for the above
    multiplication instead.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineWeight(p,w=c(0.2,0.5,0.3))
}
\author{
    Panagiotis Moulos
}
\name{estimateAufcWeights}
\alias{estimateAufcWeights}
\title{Estimate AUFC weights}
\usage{
    estimateAufcWeights(counts, normalization,
        statistics, nsim = 10, N = 10000, 
        samples = c(3, 3), ndeg = c(500, 500),
        top = 500, modelOrg = "mm9", fcBasis = 1.5,
        drawFpc = FALSE, rc = NULL,
        ...)
}
\arguments{
    \item{counts}{the real raw counts table from
    which the simulation parameters will be 
    estimated. It must not be normalized and must
    contain only integer counts, without any other
    annotation elements and unique gene identifiers
    as the rownames attribute.}

    \item{normalization}{same as \code{normalization} 
    in \code{\link{metaseqr2}}.}

    \item{statistics}{same as \code{statistics} in 
    \code{\link{metaseqr2}}.}

    \item{nsim}{the number of simulations to perform
    to estimate the weights. It default to 10.}

    \item{N}{the number of genes to produce. 
    See \code{\link{makeSimDataSd}}.}

    \item{samples}{a vector with 2 integers, which
    are the number of samples for each condition 
    (two conditions currently supported).}

    \item{ndeg}{a vector with 2 integers, which
    are the number of differentially expressed 
    genes to be produced. The first element is 
    the number of up-regulated genes while the 
    second is the number of down-regulated genes.}

    \item{fcBasis}{the minimum fold-change for 
    deregulation.}

    \item{top}{the top \code{top} best ranked 
    (according to p-value) to use, to calculate 
    area under the false discovery curve.}

    \item{modelOrg}{the organism from which the 
    data are derived. It must be one of 
    \code{\link{metaseqr2}} supported organisms.}

    \item{drawFpc}{draw the averaged false 
    discovery curves? Default to \code{FALSE}.}

    \item{rc}{the fraction of the available cores to 
    use in a multicore system.}

    \item{...}{Further arguments to be passed to 
    \code{\link{estimateSimParams}}.}
}
\value{
    A vector of weights to be used in 
    \code{\link{metaseqr2}} with the 
    \code{weights} option.
}
\description{
    This function automatically estimates weights 
    for the \code{"weight"} and \code{"dperm_weight"} 
    options of metaseqR2 for combining p-values from 
    multiple statistical tests. It creates simulated 
    dataset based on real data and then performs 
    statistical analysis with metaseqR2 several times 
    in order to derive False Discovery Curves. Then, 
    the average areas under the false discovery curves
    are used to construct weights for each algorithm,
    according to its performance when using simulated 
    data.
}
\details{
    The weight estimation process involves a lot of 
    random sampling. For guaranteed reproducibility,
    be sure to use \code{set.seed} prior to any
    calculations. By default, when the metaseqR2 package
    is loaded, the seed is set to \code{42}.
}
\examples{
require(zoo)
data("mm9GeneData",package="metaseqR2")
weights <- estimateAufcWeights(
    counts=as.matrix(mm9GeneCounts[sample(nrow(mm9GeneCounts),1000),9:12]),
    normalization="edaseq",
    statistics=c("edger","limma"),
    nsim=1,N=100,ndeg=c(10,10),top=10,modelOrg=NULL,
    rc=0.01,libsizeGt=1e+5
)
}
\author{
    Panagiotis Moulos
}
\name{diagplotVenn}
\alias{diagplotVenn}
\title{Venn diagrams when performing meta-analysis}
\usage{
    diagplotVenn(pmat, fcmat = NULL, pcut = 0.05,
        fcut = 0.5, direction = c("dereg", "up", "down"),
        nam = as.character(round(1000 * runif(1))),
        output = "x11", path = NULL, altNames = NULL, ...)
}
\arguments{
    \item{pmat}{a matrix with p-values corresponding to the
    application of each statistical algorithm. See also
    Details.}

    \item{fcmat}{an optional matrix with fold changes
    corresponding to the application of each statistical
    algorithm. See also Details.}

    \item{pcut}{if \code{fcmat} is supplied, an absolute
    fold change cutoff to be applied to \code{fcmat} to
    determine the differentially expressed genes for each
    algorithm.}

    \item{fcut}{a p-value cutoff for statistical
    significance. Defaults to \code{0.05}.}

    \item{direction}{if \code{fcmat} is supplied, a keyword
    to denote which genes to draw in the Venn diagrams with
    respect to their direction of regulation. See Details.}

    \item{nam}{a name to be appended to the output graphics
    file (if \code{"output"} is not \code{"x11"}).}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files. If 
    \code{"path"} is not \code{NULL}, a file with the 
    intersections in the Venn diagrams will be produced 
    and written in \code{"path"}.}

    \item{altNames}{an optional named vector of names, e.g.
    HUGO gene symbols, alternative or complementary to the
    unique gene names which are the rownames of \code{pmat}.
    The names of the vector must be the rownames of
    \code{pmat}.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{which.plot} argument. If
    output=\code{"x11"}, no output filenames are produced.
}
\description{
    This function uses the R package VennDiagram and plots an
    up to 5-way Venn diagram depicting the common and
    specific to each statistical algorithm genes, for each
    contrast. Mostly for internal use because of its main
    argument which is difficult to construct, but can be used
    independently if the user grasps the logic.
}
\details{
    Regarding \code{pmat}, the p-value matrix must have the 
    colnames attribute and the colnames should correspond to 
    the name of the algorithm used to fill the specific 
    column (e.g. if \code{"statistics"=c("deseq","edger",}
    \code{"nbpseq")} then \code{colnames(pmat) <-}
    \code{c("deseq","edger","nbpseq")}.
    
    Regarding \code{fcmat}, the fold change matrix must have 
    the colnames attribute and the colnames should correspond 
    to the name of the algorithm used to fill the specific 
    column (see the parameter \code{pmat}).
    
    Regarding \code{direction}, it can be one of \code{"dereg"}
    for the total of regulated genes, where 
    \code{abs(fcmat[,n])>=fcut} (default), \code{"up"} for
    the up-regulated genes where \code{fcmat[,n]>=fcut} or
    \code{"down"} for the up-regulated genes where
    \code{fcmat[,n]<=-fcut}.
}
\examples{
require(VennDiagram)
p1 <- 0.01*matrix(runif(300),100,3)
p2 <- matrix(runif(300),100,3)
p <- rbind(p1,p2)
rownames(p) <- paste("gene",1:200,sep="_")
colnames(p) <- paste("method",1:3,sep="_")
vennContents <- diagplotVenn(p)
}
\author{
    Panagiotis Moulos
}
\docType{data}
\name{hg19pvalues}
\alias{hg19pvalues}
\title{p-values from human RNA-Seq data with two conditions, four samples}
\format{a \code{matrix} with p-values from metaseqR2 supported tests.}
\source{
    Giakountis et al. (https://doi.org/10.1016/j.celrep.2016.05.038)
}
\description{
    This data set contains p-values calculated with each
    of the supported statistical testing algorithms in
    metaseqR2 for 1000 genes. The purpose of this matrix is
    to demonstrate the p-value combination methods as well as
    be used for a playground for other such methods and with
    other metaseqR2 facilities.
}
\author{
    Panagiotis Moulos
}
\keyword{datasets}
\name{buildCustomAnnotation}
\alias{buildCustomAnnotation}
\title{Import custom annotation to the metaseqR2 annotation
    database from GTF file}
\usage{
    buildCustomAnnotation(gtfFile, metadata,
    db = file.path(system.file(package = "metaseqR2"),
        "annotation.sqlite"), rewrite=TRUE)
}
\arguments{
    \item{gtfFile}{a GTF file containing the gene structure
    of the organism to be imported.}

    \item{metadata}{a list with additional information about
    the annotation to be imported. See Details.}
    
    \item{db}{a valid path (accessible at least by the
    current user) where the annotation database will be 
    set up. It defaults to 
    \code{system.file(package = "metaseqR2"),}
    \code{"annotation.sqlite")} that is, the installation
    path of metaseqR2 package. See also Details.}
    
    \item{rewrite}{if custom annotation found, rwrite? 
    (default \code{FALSE}). Set to \code{TRUE} if you wish to 
    update the annotation database for a particular custom
    annotation.}
}
\value{
    The function does not return anything. Only the SQLite 
    database is created or updated.
}
\description{
    This function imports a GTF file with some custom annotation
    to the metaseqR2 annotation database.
}
\details{
    Regarding the \code{metadata} argument, it is a list
    with specific format which instructs 
    \code{buildCustomAnnotation} on importing the custom
    annotation. Such a list may has the following members: 
    \itemize{
        \item \code{organism} a name of the organism which is
        imported (e.g. \code{"my_mm9"}). This is the only
        mandatory member.
        \item \code{source} a name of the source for this
        custom annotation (e.g. \code{"my_mouse_db"}). If
        not given or \code{NULL}, the word \code{"inhouse"}
        is used.
        \item \code{version} a string denoting the version.
        If not given or \code{NULL}, current date is used.
        \item \code{chromInfo} it can be one of the following:
        \itemize{
            \item a tab-delimited file with two columns, the
            first being the chromosome/sequence names and the
            second being the chromosome/sequence lengths.
            \item a BAM file to read the header from and
            obtain the required information
            \item a \code{\link{data.frame}} with one column
            with chromosome lengths and chromosome names as
            \code{rownames}.
        }
    }
    See the examples below for a \code{metadata} example.
    
    Regarding \code{db}, this controls the location of the
    installation database. If the default is used, then there is
    no need to provide the local database path to any function
    that uses the database (e.g. the main \code{metaseqr2}).
    Otherwise, the user will either have to provide this each
    time, or the annotation will have to be downloaded and used
    on-the-fly.
}
\examples{
# Dummy database as example
customDir <- file.path(tempdir(),"test_custom")
dir.create(customDir)

myDb <- file.path(customDir,"testann.sqlite")
chromInfo <- data.frame(length=c(1000L,2000L,1500L),
    row.names=c("A","B","C"))

# Build with the metadata list filled (you can also provide a version)
if (.Platform$OS.type == "unix") {
    buildCustomAnnotation(
        gtfFile=file.path(system.file(package="metaseqR2"),
            "dummy.gtf"),
        metadata=list(
            organism="dummy",
            source="dummy_db",
            version=1,
            chromInfo=chromInfo
        ),
        db=myDb
    )

    # Try to retrieve some data
    myGenes <- loadAnnotation(genome="dummy",refdb="dummy_db",
        level="gene",type="gene",db=myDb)
    myGenes
}

## Real data!
## Setup a temporary directory to download files etc.
#customDir <- file.path(tempdir(),"test_custom")
#dir.create(customDir)

#myDb <- file.path(customDir,"testann.sqlite")

## Gene annotation dump from Ensembl
#download.file(paste0("ftp://ftp.ensembl.org/pub/release-98/gtf/",
#  "dasypus_novemcinctus/Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
#  file.path(customDir,"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"))

## Chromosome information will be provided from the following BAM file
## available from Ensembl
#bamForInfo <- paste0("ftp://ftp.ensembl.org/pub/release-98/bamcov/",
#  "dasypus_novemcinctus/genebuild/Dasnov3.broad.Ascending_Colon_5.1.bam")

## Build with the metadata list filled (you can also provide a version)
#buildCustomAnnotation(
#  gtfFile=file.path(customDir,"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
#  metadata=list(
#    organism="dasNov3_test",
#    source="ensembl_test",
#    chromInfo=bamForInfo
#  ),
#  db=myDb
#)

## Try to retrieve some data
#dasGenes <- loadAnnotation(genome="dasNov3_test",refdb="ensembl_test",
#  level="gene",type="gene",db=myDb)
#dasGenes
}
\author{
    Panagiotis Moulos
}
\name{combineHarmonic}
\alias{combineHarmonic}
\title{Combine p-values using weights}
\usage{
    combineHarmonic(p, w, zerofix = NULL)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}

    \item{w}{a weights vector, must sum to 1.}
    
    \item{zerofix}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1.}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by
    metaseqR using p-value weights.
}
\details{
    The argument \code{zerofix} is used to correct for
    the case of a p-value which is equal to 0 as a result
    of internal numerical and approximation procedures.
    When \code{NULL}, random numbers greater than 0 and
    less than or equal to 0.5 are used to multiply the
    offending p-values with the lowest provided non-zero
    p-value, maintaining thus a virtual order of 
    significance, avoiding having the same p-values for 
    two tests and assuming that all zero p-values represent
    extreme statistical significance. When a numeric
    between 0 and 1, this number is used for the above
    multiplication instead.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineHarmonic(p,w=c(0.2,0.5,0.3))
}
\author{
    Panagiotis Moulos
}
\name{statBayseq}
\alias{statBayseq}
\title{Statistical testing with baySeq}
\usage{
    statBayseq(object, sampleList, contrastList = NULL,
        statArgs = NULL, libsizeList = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of bayseq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"bayseq")} for an
    example and how you can modify it.}

    \item{libsizeList}{an optional named list where names
    represent samples (MUST be the same as the samples in
    \code{sampleList}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    they will be estimated from baySeq.}
}
\value{
    A named list of the value 1-likelihood that a gene is
    differentially expressed, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over baySeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
# Careful on testing, baySeq is slow
require(baySeq)
dataMatrix <- metaseqR2:::exampleCountData(10)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeEdger(dataMatrix,sampleList)
p <- statBayseq(normDataMatrix,sampleList,contrast)
}
\author{
    Panagiotis Moulos
}

\name{buildAnnotationDatabase}
\alias{buildAnnotationDatabase}
\title{Build a local annotation database for metaseqR2}
\usage{
    buildAnnotationDatabase(organisms, sources,
    db = file.path(system.file(package = "metaseqR2"),
        "annotation.sqlite"),
    forceDownload = TRUE, rc = NULL)
}
\arguments{
    \item{organisms}{a list of organisms and versions
    for which to download and build annotations. Check 
    the main \code{\link{metaseqr2}} help page for details 
    on supported organisms and the Details section below.}

    \item{sources}{a character vector of public sources
    from which to download and build annotations. Check 
    the main \code{\link{metaseqr2}} help page for details 
    on supported annotation sources.}
    
    \item{db}{a valid path (accessible at least by the
    current user) where the annotation database will be 
    set up. It defaults to 
    \code{system.file(package = "metaseqR2"),}
    \code{"annotation.sqlite")} that is, the installation
    path of metaseqR2 package. See also Details.}
    
    \item{forceDownload}{by default, 
    \code{buildAnnotationDatabase} will not download an
    existing annotation again (\code{FALSE}). Set to 
    \code{TRUE} if you wish to update the annotation 
    database for a particular version.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).
    Sometimes used for building certain annotation types.}
}
\value{
    The function does not return anything. Only the SQLite 
    database is created or updated.
}
\description{
    This function creates a local annotation database to be
    used with metaseqr2 so as to avoid long time on the fly 
    annotation downloads and formatting.
}
\details{
    Regarding the \code{organisms} argument, it is a list
    with specific format which instructs 
    \code{buildAnnotationDatabase} on which organisms and
    versions to download from the respective sources. Such
    a list may have the format: 
    \code{organisms=list(hg19=75, mm9=67, mm10=96:97)}
    This is explained as follows:
    \itemize{
        \item A database comprising the human genome versions
        \code{hg19} and the mouse genome versions 
        \code{mm9, mm10} will be constructed.
        \item If \code{"ensembl"} is in \code{sources}, 
        version 75 is downloaded for \code{hg19} and versions 
        \code{67, 96, 97} for \code{mm9, mm10}. 
        \item If \code{"ucsc"} or \code{"refseq"} are in 
        \code{sources}, the latest versions are downloaded
        and marked by the download date. As UCSC and RefSeq
        versions are not accessible in the same way as
        Ensembl, this procedure cannot always be replicated.
    }
    \code{organisms} can also be a character vector with organism
    names/versions (e.g. \code{organisms = c("mm10","hg19")}),
    then the latest versions are downloaded in the case of 
    Ensembl.
    
    Regarding \code{db}, this controls the location of the
    installation database. If the default is used, then there is
    no need to provide the local database path to any function
    that uses the database (e.g. the main \code{metaseqr2}).
    Otherwise, the user will either have to provide this each
    time, or the annotation will have to be downloaded and used
    on-the-fly.
}
\examples{
# Build a test database with one genome
myDb <- file.path(tempdir(),"testann.sqlite")

organisms <- list(mm10=75)
sources <- "ensembl"

# If the example is not running in a multicore system, rc is ignored
#buildAnnotationDatabase(organisms,sources,db=myDb,rc=0.5)

# A more complete case, don't run as example
# Since we are using Ensembl, we can also ask for a version
#organisms <- list(
#    mm9=67,
#    mm10=96:97,
#    hg19=75,
#    hg38=96:97
#)
#sources <- c("ensembl", "refseq")

## Build on the default location (depending on package location, it may
## require root/sudo)
#buildAnnotationDatabase(organisms,sources)

## Build on an alternative location
#myDb <- file.path(path.expand("~"),"my_ann.sqlite")
#buildAnnotationDatabase(organisms,sources,db=myDb)
}
\author{
    Panagiotis Moulos
}
\name{normalizeDss}
\alias{normalizeDss}
\title{Normalization based on the DSS package}
\usage{
    normalizeDss(geneCounts, sampleList,
        normArgs = NULL, output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of DESeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"deseq")} for
    an example and how you can modify it.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the DSS native S4 object
    (SeqCountSet). In the latter case it should be handled
    with suitable ABSSeq methods.}
}
\value{
    A matrix or a SeqCountSet with normalized counts.
}
\description{
    This function is a wrapper over ABSSeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeDss(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Dionysios Fanidis
}
\name{downsampleCounts}
\alias{downsampleCounts}
\title{Downsample read counts}
\usage{
    downsampleCounts(counts)
}
\arguments{
    \item{counts}{the read counts table 
    which is subjected to downsampling.}
}
\value{
    The downsampled counts matrix.
}
\description{
    This function downsamples the library sizes 
    of a read counts table to the lowest library 
    size, according to the methdology used in 
    (Soneson and Delorenzi, BMC Bioinformatics, 
    2013).
}
\details{
    The downsampling process involves random sampling. 
    For guaranteed reproducibility, be sure to use 
    \code{set.seed} before downsampling. By default, 
    when the metaseqR2 package is loaded, the seed is 
    set to \code{42}.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(5000)
D <- downsampleCounts(dataMatrix)
}
\author{
    Panagiotis Moulos
}
\name{diagplotBoxplot}
\alias{diagplotBoxplot}
\title{Boxplots wrapper for the metaseqR2 package}
\usage{
    diagplotBoxplot(mat, name = NULL, logIt = "auto",
        yLim = "default", isNorm = FALSE, output = "x11",
        path = NULL, altNames = NULL, ...)
}
\arguments{
    \item{mat}{the count data matrix.}

    \item{name}{the names of the samples plotted on the
    boxplot. See also Details.}

    \item{logIt}{whether to log transform the values of mat
    or not. It can be \code{TRUE}, \code{FALSE} or
    \code{"auto"} for auto-detection. Auto-detection log
    transforms by default so that the boxplots are smooth and
    visible.}

    \item{yLim}{custom y-axis limits. Leave the string
    \code{"default"} for default behavior.}

    \item{isNorm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"}, \code{"ps"} or \code{"json"}. The latter is
    currently available for the creation of interactive
    volcano plots only when reporting the output, through the
    highcharts javascript library (JSON for boxplots not yet
    available).}

    \item{path}{the path to create output files.}
    
    \item{altNames}{an optional vector of names, e.g. HUGO
    gene symbols, alternative or complementary to the unique
    rownames of \code{mat} (which must exist!). It is used only 
    in JSON output.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filename of the boxplot produced if it's a file.
}
\description{
    A wrapper over the general boxplot function, suitable for
    matrices produced and processed with the metaseqr
    package. Intended for internal use but can be easily used
    as stand-alone. It can colors boxes based on group
    depending on the name argument.
}
\details{
    Regarding \code{name}, if \code{NULL}, the function check 
    the column names of \code{mat}. If they are also 
    \code{NULL}, sample names are autogenerated. If 
    \code{name="none"}, no sample names are plotted. If name 
    is a list, it should be the sampleList argument provided 
    to the manin metaseqr2 function. In that case, the boxes  
    are colored per group.
}
\examples{
# Non-normalized boxplot
require(DESeq2)
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

# Normalized boxplot
normArgs <- getDefaults("normalization","deseq2")
object <- normalizeDeseq2(dataMatrix,sampleList,normArgs)
diagplotBoxplot(object,sampleList)
}
\author{
    Panagiotis Moulos
}
\name{diagplotPairs}
\alias{diagplotPairs}
\title{Massive X-Y, M-D correlation plots}
\usage{
    diagplotPairs(x, output = "x11", altNames = NULL, 
        path = NULL, ...)
}
\arguments{
    \item{x}{the read counts matrix or data frame.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}
    
    \item{altNames}{optional names, alternative or complementary 
    to the rownames of \code{x}. It is used only in JSON output.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filename of the pairwise comparisons plot produced if
    it's a file.
}
\description{
    This function uses the read counts matrix to create
    pairwise correlation plots. The upper diagonal of the
    final image contains simple scatterplots of each sample
    against each other (log2 scale) while the lower diagonal
    contains mean-difference plots for the same samples (log2
    scale). This type of diagnostic plot may not be
    interpretable for more than 10 samples.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
diagplotPairs(dataMatrix)
}
\author{
    Panagiotis Moulos
}

\name{statDeseq}
\alias{statDeseq}
\title{Statistical testing with DESeq}
\usage{
    statDeseq(object, sampleList, contrastList = NULL,
        statArgs = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of DESeq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"deseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a CountDataSet from
    DESeq normalization as the dispersions are already
    estimated.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over DESeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(1000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- metaseqR2:::newCountDataSet(dataMatrix,
    c("A","A","B","B","B"))
normDataMatrix <- normalizeDeseq(dataMatrix,sampleList)
p <- statDeseq(normDataMatrix,sampleList,contrast)
}
\author{
    Panagiotis Moulos
}

\name{statNoiseq}
\alias{statNoiseq}
\title{Statistical testing with NOISeq}
\usage{
    statNoiseq(object, sampleList, contrastList = NULL,
        statArgs = NULL, geneData = NULL, logOffset = 1)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of edgeR statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"noiseq")} for an
    example and how you can modify it.}

    \item{geneData}{an optional annotation data frame (such
    the ones produced by \code{get.annotation} which contains
    the GC content for each gene and from which the gene
    lengths can be inferred by chromosome coordinates.}

    \item{logOffset}{a number to be added to each element of
    data matrix in order to avoid Infinity on log type data
    transformations.}
}
\value{
    A named list of NOISeq q-values, whose names are the
    names of the contrasts.
}
\description{
    This function is a wrapper over NOISeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(1000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
lengths <- round(1000*runif(nrow(dataMatrix)))
starts <- round(1000*runif(nrow(dataMatrix)))
ends <- starts + lengths
gc=runif(nrow(dataMatrix))
biotype=rep("protein_coding",nrow(dataMatrix))
geneData <- data.frame(
    chromosome=c(rep("chr1",nrow(dataMatrix)/2),
    rep("chr2",nrow(dataMatrix)/2)),
        start=starts,end=ends,gene_id=rownames(dataMatrix),
    gc_content=gc,biotype=biotype
)
normArgs <- metaseqR2:::getDefaults("normalization","noiseq")
normDataMatrix <- normalizeNoiseq(dataMatrix,sampleList,normArgs,
    geneData)
p <- statNoiseq(normDataMatrix,sampleList,contrast,
    geneData=geneData)
}
\author{
    Panagiotis Moulos
}
\name{normalizeNbpseq}
\alias{normalizeNbpseq}
\title{Normalization based on the NBPSeq package}
\usage{
    normalizeNbpseq(geneCounts, sampleList,
        normArgs = NULL, libsizeList = NULL,
        output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of NBPSeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"nbpseq")} for
    an example and how you can modify it.}

    \item{libsizeList}{an optional named list where names
    represent samples (MUST be the same as the samples in
    \code{sampleList}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    the default is the column sums of the \code{geneCounts}
    matrix.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the NBPSeq native S4 object
    (a specific list). In the latter case it should be
    handled with suitable NBPSeq methods.}
}
\value{
    A matrix with normalized counts or a list with the
    normalized counts and other NBPSeq specific parameters.
}
\description{
    This function is a wrapper over DESeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeNbpseq(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}
\name{combineBonferroni}
\alias{combineBonferroni}
\title{Combine p-values with Bonferroni's method}
\usage{
    combineBonferroni(p, zerofix = NULL)
}
\arguments{
    \item{p}{a vector of p-values for each statistical 
    tests).}
    
    \item{zerofix}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1.}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by 
    metaseqR2 using the Bonferroni's method (see 
    reference in the main \code{\link{metaseqr2}}
    help page or in the vignette).
}
\details{
    The argument \code{zerofix} is used to correct for
    the case of a p-value which is equal to 0 as a result
    of internal numerical and approximation procedures.
    When \code{NULL}, random numbers greater than 0 and
    less than or equal to 0.5 are used to multiply the
    offending p-values with the lowest provided non-zero
    p-value, maintaining thus a virtual order of 
    significance, avoiding having the same p-values for 
    two tests and assuming that all zero p-values represent
    extreme statistical significance. When a numeric
    between 0 and 1, this number is used for the above
    multiplication instead.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineBonferroni(p)
}
\author{
    Panagiotis Moulos
}

\docType{data}
\name{mm9GeneCounts}
\alias{mm9GeneCounts}
\title{Mouse RNA-Seq data with two conditions, four samples}
\format{a \code{data.frame} with gene read counts and some embedded 
annotation, one row per gene.}
\source{
    ENCODE (http://genome.ucsc.edu/encode/)
}
\description{
    This data set contains RNA-Seq gene read counts for 3
    chromosomes. The data were downloaded from the ENCODE
    public repository and are derived from the study of
    Mortazavi et al., 2008 (Mortazavi A, Williams BA, McCue
    K, Schaeffer L, Wold B. Mapping and quantifying mammalian
    transcriptomes by RNA-Seq. Nat Methods. 2008
    Jul;5(7):621-8). In their experiment, the authors studied
    among others genes expression at two developmental stages
    of mouse liver cells. It has two conditions-developmental
    stages (e14.5, adult_8_weeks) and four samples (e14.5_1,
    e14.5_2, a8w_1, a8w_2). It also contains a predefined
    \code{sampleList} and \code{libsizeList} named
    \code{sampleListMm9} and \code{libsizeListMm9}.
}
\author{
    Panagiotis Moulos
}
\keyword{datasets}
\name{makeSimDataSd}
\alias{makeSimDataSd}
\title{Create simulated counts using the 
    Soneson-Delorenzi method}
\usage{
    makeSimDataSd(N, param, samples = c(5, 5),
        ndeg = rep(round(0.1*N), 2), fcBasis = 1.5,
        libsizeRange = c(0.7, 1.4), libsizeMag = 1e+7,
        modelOrg = NULL, simLengthBias = FALSE)
}
\arguments{
    \item{N}{the number of genes to produce.}

    \item{param}{a named list with negative binomial 
    parameter sets to sample from. The first member is
    the mean parameter to sample from (\code{muHat}) 
    and the second the dispersion (\code{phiHat}). 
    This list can be created with the 
    \code{\link{estimateSimParams}} function.}

    \item{samples}{a vector with 2 integers, 
    which are the number of samples for each 
    condition (two conditions currently supported).}

    \item{ndeg}{a vector with 2 integers, which are 
    the number of differentially expressed genes to 
    be produced. The first element is the number of 
    up-regulated genes while the second is the 
    number of down-regulated genes.}

    \item{fcBasis}{the minimum fold-change for 
    deregulation.}

    \item{libsizeRange}{a vector with 2 numbers 
    (generally small, see the default), as they 
    are multiplied with \code{libsizeMag}. These 
    numbers control the library sized of the 
    synthetic data to be produced.}

    \item{libsizeMag}{a (big) number to multiply 
    the \code{libsizeRange} to produce library 
    sizes.}

    \item{modelOrg}{the organism from which the 
    real data are derived from. It must be one 
    of the supported organisms (see the main 
    \code{\link{metaseqr2}} help page). It is used 
    to sample real values for GC content.}

    \item{simLengthBias}{a boolean to instruct 
    the simulator to create genes whose read counts is
    proportional to their length. This is achieved by 
    sorting in increasing order the mean parameter of 
    the negative binomial distribution (and the 
    dispersion according to the mean) which will cause 
    an increasing gene count length with the sampling. 
    The sampled lengths are also sorted so that in the 
    final gene list, shorter genes have less counts as 
    compared to the longer ones. The default is FALSE.}
}
\value{
    A named list with two members. The first 
    member (\code{simdata}) contains the 
    synthetic dataset 
}
\description{
    This function creates simulated RNA-Seq gene 
    expression datasets using the method presented 
    in (Soneson and Delorenzi, BMC Bioinformatics, 
    2013). For the time being, it creates only 
    simulated datasets with two conditions.
}
\details{
    The simulated data generation involves a lot of 
    random sampling. For guaranteed reproducibility,
    be sure to use \code{set.seed} prior to any
    calculations. By default, when the metaseqR2 package
    is loaded, the seed is set to \code{42}.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
## File "bottomly_read_counts.txt" from the ReCount database
#download.file(paste("http://bowtie-bio.sourceforge.net/recount/",
#    "countTables/bottomly_count_table.txt",sep=""),
#    destfile="~/bottomly_count_table.txt")
N <- 2000
#parList <- estimateSimParams("~/bottomly_read_counts.txt")
parList <- estimateSimParams(dataMatrix,libsizeGt=3e+4)
sim <- makeSimDataSd(N,parList)
synthData <- sim$simdata
trueDeg <- which(sim$truedeg!=0)
}
\author{
    Panagiotis Moulos
}
\name{importCustomAnnotation}
\alias{importCustomAnnotation}
\title{Import a metaseqR2 custom annotation element}
\usage{
    importCustomAnnotation(gtfFile, metadata,
        level = c("gene", "transcript", "exon"),
        type = c("gene", "exon", "utr"))
}
\arguments{
    \item{gtfFile}{a GTF file containing the gene structure
    of the organism to be imported.}

    \item{metadata}{a list with additional information about
    the annotation to be imported. The same as in the
    \code{\link{buildCustomAnnotation}} man page.}
    
    \item{level}{same as the \code{transLevel} in 
    \code{\link{metaseqr2}}.}
    
    \item{type}{same as the \code{countType} in 
    \code{\link{metaseqr2}}.}
}
\value{
    The function returns a \code{GenomicRanges} object with
    the requested annotation.
}
\description{
    This function creates a local annotation database to be
    used with metaseqr2 so as to avoid long time on the fly 
    annotation downloads and formatting.
}
\examples{
# Dummy GTF as example
chromInfo <- data.frame(length=c(1000L,2000L,1500L),
    row.names=c("A","B","C"))

# Build with the metadata list filled (you can also provide a version)
myGenes <- importCustomAnnotation(
    gtfFile=file.path(system.file(package="metaseqR2"),"dummy.gtf"),
    metadata=list(
        organism="dummy",
        source="dummy_db",
        version=1,
        chromInfo=chromInfo
    ),
    level="gene",type="gene"
)

## Real data!
## Gene annotation dump from Ensembl
#download.file(paste0("ftp://ftp.ensembl.org/pub/release-98/gtf/",
#  "dasypus_novemcinctus/Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
#  file.path(tempdir(),"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"))

## Build with the metadata list filled (you can also provide a version)
#dasGenes <- importCustomAnnotation(
#  gtfFile=file.path(tempdir(),"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
#  metadata=list(
#    organism="dasNov3_test",
#    source="ensembl_test"
#  ),
#  level="gene",type="gene"
#)
}
\author{
    Panagiotis Moulos
}
\name{combineMinp}
\alias{combineMinp}
\title{Combine p-values using the minimum p-value}
\usage{
    combineMinp(p)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by
    metaseqR by taking the minimum p-value.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineMinp(p)
}
\author{
    Panagiotis Moulos
}

\name{normalizeDeseq2}
\alias{normalizeDeseq2}
\title{Normalization based on the DESeq2 package}
\usage{
    normalizeDeseq2(geneCounts, sampleList,
        normArgs = NULL, output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of DESeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"deseq")} for
    an example and how you can modify it.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the DESeq2 native S4 object
    (DESeqDataSet). In the latter case it should be handled
    with suitable DESeq2 methods.}
}
\value{
    A matrix or a DESeqDataSet with normalized counts.
}
\description{
    This function is a wrapper over DESeq2 normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeDeseq2(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Dionysios Fanidis
}
\name{statLimma}
\alias{statLimma}
\title{Statistical testing with limma}
\usage{
    statLimma(object, sampleList, contrastList = NULL,
        statArgs = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of edgeR statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"limma")} for an
    example and how you can modify it.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over limma statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
require(limma)
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeEdger(dataMatrix,sampleList)
p <- statLimma(normDataMatrix,sampleList,contrast)
}
\author{
    Panagiotis Moulos
}

\name{normalizeEdaseq}
\alias{normalizeEdaseq}
\title{Normalization based on the EDASeq package}
\usage{
    normalizeEdaseq(geneCounts, sampleList,
        normArgs = NULL, geneData = NULL,
        output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of EDASeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"edaseq")} for
    an example and how you can modify it.}

    \item{geneData}{an optional annotation data frame (such
    the ones produced by \code{getAnnotation}) which
    contains the GC content for each gene and from which the
    gene lengths can be inferred by chromosome coordinates.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the EDASeq native S4 object
    (SeqExpressionSet). In the latter case it should be
    handled with suitable EDASeq methods.}
}
\value{
    A matrix or a SeqExpressionSet with normalized counts.
}
\description{
    This function is a wrapper over EDASeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

lengths <- round(1000*runif(nrow(dataMatrix)))
starts <- round(1000*runif(nrow(dataMatrix)))
ends <- starts + lengths
gc=runif(nrow(dataMatrix))
geneData <- data.frame(
    chromosome=c(rep("chr1",nrow(dataMatrix)/2),
        rep("chr2",nrow(dataMatrix)/2)),
    start=starts,end=ends,gene_id=rownames(dataMatrix),gc_content=gc,
    row.names=rownames(dataMatrix)
)
normDataMatrix <- normalizeEdaseq(dataMatrix,sampleList,
    geneData=geneData)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}
\name{diagplotAvgFtd}
\alias{diagplotAvgFtd}
\title{Create average False (or True) Discovery 
    curves}
\usage{
    diagplotAvgFtd(ftdrObj, output = "x11",
        path = NULL, draw = TRUE, ...)
}
\arguments{
    \item{ftdrObj}{a list with outputs from 
    \code{\link{diagplotFtd}}.}

    \item{output}{one or more R plotting 
    device to direct the plot result to.
    Supported mechanisms: \code{"x11"} (default), 
    \code{"png"}, \code{"jpg"}, \code{"bmp"}, 
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{draw}{boolean to determine whether 
    to plot the curves or just return the 
    calculated values (in cases where the user 
    wants the output for later averaging
    for example). Defaults to \code{TRUE} (make 
    plots).}

    \item{...}{further arguments to be passed to 
    plot devices, such as parameter from 
    \code{\link{par}}.}
}
\value{
    A named list with two members: the first member 
    (\code{avgFtdr}) contains a list with the 
    means and the standard deviations of the averaged
    \code{ftdrObj} and are used to create the plot. 
    The second member (\code{path}) contains the 
    path to the created figure graphic.
}
\description{
    This function creates false (or true) discovery 
    curves using a list containing several outputs 
    from \code{\link{diagplotFtd}}.
}
\examples{
p11 <- 0.001*matrix(runif(300),100,3)
p12 <- matrix(runif(300),100,3)
p21 <- 0.001*matrix(runif(300),100,3)
p22 <- matrix(runif(300),100,3)
p31 <- 0.001*matrix(runif(300),100,3)
p32 <- matrix(runif(300),100,3)
p1 <- rbind(p11,p21)
p2 <- rbind(p12,p22)
p3 <- rbind(p31,p32)
rownames(p1) <- rownames(p2) <- rownames(p3) <-
    paste("gene",1:200,sep="_")
colnames(p1) <- colnames(p2) <- colnames(p3) <-
    paste("method",1:3,sep="_")
truth <- c(rep(1,40),rep(-1,40),rep(0,20),
    rep(1,10),rep(2,10),rep(0,80))
names(truth) <- rownames(p1)
ftdObj1 <- diagplotFtd(truth,p1,N=100,draw=FALSE)
ftdObj2 <- diagplotFtd(truth,p2,N=100,draw=FALSE)
ftdObj3 <- diagplotFtd(truth,p3,N=100,draw=FALSE)
ftdObj <- list(ftdObj1,ftdObj2,ftdObj3)
avgFtdObj <- diagplotAvgFtd(ftdObj)
}
\author{
    Panagiotis Moulos
}
\name{normalizeNoiseq}
\alias{normalizeNoiseq}
\title{Normalization based on the NOISeq package}
\usage{
    normalizeNoiseq(geneCounts, sampleList,
        normArgs = NULL, geneData = NULL, logOffset = 1,
        output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of NOISeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"noiseq")} for
    an example and how you can modify it.}

    \item{geneData}{an optional annotation data frame (such
    the ones produced by \code{get.annotation} which contains
    the GC content for each gene and from which the gene
    lengths can be inferred by chromosome coordinates.}

    \item{logOffset}{an offset to use to avoid infinity in
    logarithmic data transformations.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the NOISeq native S4 object
    (SeqExpressionSet). In the latter case it should be
    handled with suitable NOISeq methods.}
}
\value{
    A matrix with normalized counts.
}
\description{
    This function is a wrapper over NOISeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

lengths <- round(1000*runif(nrow(dataMatrix)))
starts <- round(1000*runif(nrow(dataMatrix)))
ends <- starts + lengths
gc=runif(nrow(dataMatrix))
geneData <- data.frame(
    chromosome=c(rep("chr1",nrow(dataMatrix)/2),
        rep("chr2",nrow(dataMatrix)/2)),
    start=starts,end=ends,gene_id=rownames(dataMatrix),gc_content=gc,
    biotype=rep("protein_coding",nrow(dataMatrix)),
    row.names=rownames(dataMatrix)
)
normDataMatrix <- normalizeNoiseq(dataMatrix,sampleList,normArgs=NULL,geneData)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}

\name{makeSimDataTcc}
\alias{makeSimDataTcc}
\title{Create simulated counts using TCC package}
\usage{
    makeSimDataTcc(...)
}
\arguments{
    \item{...}{parameters to the \code{simulateReadCounts}
    function.}
}
\value{
    A list with the following members: \code{simdata} holding
    the simulated dataset complying with metaseqr2
    requirements, and \code{simparam} holding the simulation
    parameters (see TCC documentation). Note that the produced 
    data are based in an Arabidopsis dataset.
}
\description{
    This function creates simulated RNA-Seq gene expression
    datasets using the \code{simulateReadCounts} function
    from the Bioconductor package TCC and it adds simulated
    annoation elements. For further information please
    consult the TCC package documentation.
}
\examples{
if (require(TCC)) {
dd <- makeSimDataTcc(Ngene=1000,PDEG=0.2,
    DEG.assign=c(0.9,0.1),
    DEG.foldchange=c(5,5),replicates=c(3,3))
head(dd$simdata)
}
}
\author{
    Panagiotis Moulos
}
\name{combineMaxp}
\alias{combineMaxp}
\title{Combine p-values using the maximum p-value}
\usage{
    combineMaxp(p)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by
    metaseqR by taking the maximum p-value.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineMaxp(p)
}
\author{
    Panagiotis Moulos
}

\name{readTargets}
\alias{readTargets}
\title{Creates sample list and BAM/BED file list from file}
\usage{
    readTargets(input, path = NULL)
}
\arguments{
    \item{input}{a tab-delimited file or a YAML file 
    specifically structured. See Details.}

    \item{path}{an optional path where all the BED/BAM 
    files are placed, to be prepended to the BAM/BED 
    file names in the targets file.}
}
\value{
    A named list with four members. The first member is 
    a named list whose names are the conditions of the 
    experiments and its members are the samples belonging 
    to each condition. The second member is like the 
    first, but this time the members are named vectors 
    whose names are the sample names and the vector 
    elements are full path to BAM/BED files. The third 
    member is like the second, but instead of filenames 
    it contains information about single- or paired-end
    reads (if available). The fourth member is like the 
    second, but instead of filenames it contains 
    information about the strandedness of the reads (if
    available). The fifth member is the guessed type 
    of the input files (SAM/BAM or BED). It will be used 
    if not given in the main \code{\link{read2count}}
    function.
}
\description{
    Create the main sample list and determine the BAM/BED
    files for each sample from an external file.
}
\details{
    Regarding the input file, this can be a simple text 
    tab-delimited file or a YAML file describing the data to
    be analyzed. 
    
    Regarding the tab-delimited version, its columns must
    be structured as follows: the first line of the external
    tab delimited file should contain column names (names
    are not important). The first column MUST contain UNIQUE
    sample names. The second column MUST contain the raw 
    BAM/BED files WITH their full path. Alternatively, the 
    \code{path} argument should be provided (see below). The 
    third column MUST contain the biological condition where 
    each of the samples in the first column should belong to. 
    There is an optional fourth column which should contain 
    the keywords \code{"single"} for single-end reads, 
    \code{"paired"} for paired-end reads or \code{"mixed"} 
    for BAM files that contain both single- and paired-end 
    reads (e.g. after a mapping procedure with two round of 
    alignment). If this column is not provided, single-end 
    reads will be assumed. There is an optional fifth column 
    which stranded read assignment. It should contain the 
    keywords \code{"forward"} for a forward (5'->3') strand 
    library construction protocol, \code{"reverse"} for a 
    reverse (3'->5') strand library construction protocol, 
    or \code{"no"} for unstranded/unknown protocol. If this 
    column is not provided, unstranded reads will be assumed.
    
    Regarding the YAML version, the same instructions apply,
    but this time instead of columns, the data are provided
    as a YAML array under a keyword/top-level field 
    representing the respective header in the tab-delimited
    version. Alternatively, the aforementioned structure
    can be nested under a root level named strictly either
    \code{targets} or \code{metaseqR2_targets}. The latter
    can be especially useful when incorporating the
    metaseqR2 pipeline in a wider pipeline including various
    analyses and described using a workflow language such
    as CWL.
}
\examples{
dataPath <- system.file("extdata",package="metaseqR2")
targets <- data.frame(samplename=c("C","T"),
    filename=file.path(dataPath,c("C.bam","T.bam")),  
    condition=c("Control","Treatment"),
    paired=c("single","single"),stranded=c("forward","forward"))
path <- tempdir()

# Tab delimited case
write.table(targets,file=file.path(path,"targets.txt"),
    sep="\t",row.names=FALSE,quote=FALSE)
theList <- readTargets(file.path(path,"targets.txt"),path=path)
sampleList <- theList$samples
bamfileList <- theList$files

# YAML case
require(yaml)
write_yaml(as.list(targets),file.path(path,"targets.yml"))
theYList <- readTargets(file.path(path,"targets.yml"),path=path)
identical(theList,theYList) # TRUE

# YAML case with nested targets
write_yaml(list(targets=as.list(targets)),
    file.path(path,"targets2.yml"))
theYList2 <- readTargets(file.path(path,"targets2.yml"),path=path)
identical(theYList,theYList2) # TRUE
}
\author{
    Panagiotis Moulos
}

\name{statEdger}
\alias{statEdger}
\title{Statistical testing with edgeR}
\usage{
    statEdger(object, sampleList, contrastList = NULL,
        statArgs = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of edgeR statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"edger")} for an
    example and how you can modify it.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over edgeR statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
require(edgeR)
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeEdger(dataMatrix,sampleList)
p <- statEdger(normDataMatrix,sampleList,contrast)
}
\author{
    Panagiotis Moulos
}

\name{metaTest}
\alias{metaTest}
\title{Meta-analysis using several RNA-Seq statistics}
\usage{
    metaTest(cpList,
        metaP = c("simes", "bonferroni", "fisher", "harmonic",
        "dperm_min", "dperm_max", "dperm_weight", "fperm", 
        "whitlock", "minp", "maxp", "weight", "pandora", 
        "none"), counts, sampleList, statistics, statArgs, 
        libsizeList, nperm = 10000, 
        weight = rep(1/length(statistics), length(statistics)), 
        pOffset = NULL, rc = NULL)
}
\arguments{
    \item{cpList}{a named list whose names are the contrasts
    requested from metaseqr2. Each member is a p-value matrix
    whose colnames are the names of the statistical tests
    applied to the data. See the main \code{\link{metaseqr2}}
    help page.}

    \item{metaP}{the p-value combination method to use. See
    the main \code{\link{metaseqr2}} help page.}

    \item{counts}{the normalized and possibly filtered read
    counts matrix. See the main \code{\link{metaseqr2}} help
    page.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition. See the main
    \code{\link{metaseqr2}} help page.}

    \item{statistics}{the statistical algorithms used in
    metaseqr2. See the main \code{\link{metaseqr2}} help page.}

    \item{statArgs}{the parameters for each statistical
    argument. See the main \code{\link{metaseqr2}} help page.}

    \item{libsizeList}{a list with library sizes. See the
    main \code{\link{metaseqr2}} and the \code{stat*} help
    pages.}

    \item{nperm}{the number of permutations (Monte Carlo
    simulations) to perform.}

    \item{weight}{a numeric vector of weights for each
    statistical algorithm.}
    
    \item{pOffset}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1. See also the main 
    \code{\link{metaseqr2}} man page.}

    \item{rc}{the fraction of the available cores to use
    in a multicore system.}
}
\value{
    A named list with combined p-values. The names are the
    contrasts and the list members are combined p-value
    vectors, one for each contrast.
}
\description{
    This function calculates the combined p-values when
    multiple statistical algorithms are applied to the input
    dataset. It is a helper and it requires very specific
    arguments so it should not be used individually
}
\details{
    Ideally one would want to create the same set of indices 
    for a given dataset so as to create reproducible p-values. 
    To achieve this, use the \code{set.seed} function prior 
    to any calculations.
}
\examples{
cpList <- list(a=matrix(runif(100),50,2))
metaP <- metaTest(cpList,"simes")
}
\author{
    Panagiotis Moulos
}
\name{diagplotNoiseq}
\alias{diagplotNoiseq}
\title{Diagnostic plots based on the NOISeq package}
\usage{
    diagplotNoiseq(x, sampleList, covars,
        whichPlot = c("biodetection", "countsbio", "saturation", 
        "rnacomp", "readnoise", "biodist"),
        output = "x11",
        biodistOpts = list(p = NULL, pcut = NULL, name = NULL),
        path = NULL, isNorm = FALSE, ...)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{covars}{a list (whose annotation elements are
    ideally a subset of an annotation data frame produced by
    \code{\link{getAnnotation}}) with the following members:
    data (the data matrix), length (gene length), gc (the
    gene gc_content), chromosome (a data frame with
    chromosome name and co-ordinates), factors (a factor with
    the experimental condition names replicated by the number
    of samples in each experimental condition) and biotype
    (each gene's biotype as depicted in Ensembl-like
    annotations).}

    \item{whichPlot}{the NOISeq package plot to generate. 
    See Details}

    \item{biodistOpts}{a list with the following members: p
    (a vector of p-values, e.g. the p-values of a contrast),
    pcut (a unique number depicting a p-value cutoff,
    required for the \code{"biodist"} case), name (a name for
    the \code{"biodist"} plot, e.g. the name of the
    contrast.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{isNorm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{whichPlot} argument. If
    \code{output="x11"}, no output filenames are produced.
}
\description{
    A wrapper around the plotting functions availale in the
    NOISeq Bioconductor package. For analytical explanation
    of each plot please see the vignette of the NOISeq 
    package. It is best to use this function through the 
    main plotting function \code{\link{metaseqrPlot}}.
}
\details{
    Regarding \code{whichPlot}, It can be one or more of 
    \code{"biodetection"}, \code{"countsbio"}, 
    \code{"saturation"}, \code{"rnacomp"}, \code{"readnoise"} 
    or \code{"biodist"}. Please refer to the documentation of 
    the NOISeq package for details on the use of these plots. 
    The \code{whichPlot="saturation"} case is modified to be
    more informative by producing two kinds of plots.
}
\note{
    Please note that in case of \code{"biodist"} plots, the
    behavior of the function is unstable, mostly due to the
    very specific inputs this plotting function accepts in
    the NOISeq package. We have tried to predict unstable
    behavior and avoid exceptions through the use of tryCatch
    but it's still possible that you might run onto an error.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(5000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
lengths <- round(1000*runif(nrow(dataMatrix)))
starts <- round(1000*runif(nrow(dataMatrix)))
ends <- starts + lengths
covars <- list(
    data=dataMatrix,
    length=lengths,
    gc=runif(nrow(dataMatrix)),
    chromosome=data.frame(
        chromosome=c(rep("chr1",nrow(dataMatrix)/2),
        rep("chr2",nrow(dataMatrix)/2)),
        start=starts,
        end=ends
    ),
    factors=data.frame(class=metaseqR2:::asClassVector(sampleList)),
    biotype=c(rep("protein_coding",nrow(dataMatrix)/2),rep("ncRNA",
        nrow(dataMatrix)/2))
)
p <- runif(nrow(dataMatrix))
diagplotNoiseq(dataMatrix,sampleList,covars=covars,
    biodistOpts=list(p=p,pcut=0.1,name="A_vs_B"))
}
\author{
    Panagiotis Moulos
}
\name{diagplotEdaseq}
\alias{diagplotEdaseq}
\title{Diagnostic plots based on the EDASeq package}
\usage{
    diagplotEdaseq(x, sampleList, covar = NULL,
        isNorm = FALSE,
        whichPlot = c("meanvar", "meandiff", "gcbias", "lengthbias"),
        output = "x11", altNames = NULL, path = NULL, ...)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{covar}{The covariate to plot counts against.
    Usually \code{"gc"} or \code{"length"}.}

    \item{isNorm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{whichPlot}{the EDASeq package plot to generate. 
    See Details.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}
    
    \item{altNames}{optional names, alternative or complementary 
    to the rownames of \code{x}. It is used only in JSON output.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plot produced in a named list with
    names the which.plot argument. If \code{output="x11"}, no
    output filenames are produced.
}
\description{
    A wrapper around the plotting functions availale in the
    EDASeq normalization Bioconductor package. For analytical
    explanation of each plot please see the vignette of the
    EDASeq package. It is best to use this function through
    the main plotting function
    \code{\link{metaseqrPlot}}.
}
\details{
    Regarding \code{whichPlot}, it can be one or more of 
    \code{"meanvar"}, \code{"meandiff"}, \code{"gcbias"} or
    \code{"lengthbias"}. Please refer to the documentation of
    the EDASeq package for details on the use of these 
    plots. The \code{whichPlot="lengthbias"} case is 
    not covered by EDASeq documentation, however it is 
    similar to the GC-bias plot when the covariate is the 
    gene length instead of the GC content.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotEdaseq(dataMatrix,sampleList,whichPlot="meandiff")
}
\author{
    Panagiotis Moulos
}
\name{estimateSimParams}
\alias{estimateSimParams}
\title{Estimate negative binomial parameters from 
    real data}
\usage{
    estimateSimParams(realCounts, libsizeGt = 3e+6,
        rowmeansGt = 5,eps = 1e-11, rc = NULL, draw = FALSE)
}
\arguments{
    \item{realCounts}{a text tab-delimited file 
    with real RNA-Seq data. See Details.}

    \item{libsizeGt}{a library size below which 
    samples are excluded from parameter estimation 
    (default: \code{3000000}).}

    \item{rowmeansGt}{a row means (mean counts 
    over samples for each gene) below which 
    genes are excluded from parameter estimation
    (default: 5).}

    \item{eps}{the tolerance for the convergence 
    of \code{\link{optimize}} function. Defaults 
    to 1e-11.}

    \item{rc}{in case of parallel optimization, the 
    fraction of the available cores to use.}

    \item{draw}{boolean to determine whether to 
    plot the estimated simulation parameters 
    (mean and dispersion) or not. Defaults to 
    \code{FALSE} (do not draw a mean-dispersion 
    scatterplot).}
}
\value{
    A named list with two members: \code{muHat}
    which contains negative binomial mean 
    estimates and \code{phiHat} which contains 
    dispersion estimates.
}
\description{
    This function reads a read counts table 
    containing real RNA-Seq data (preferebly 
    with more than 20 samples so as to get as 
    much accurate as possible estimations) and 
    calculates a population of count means and 
    dispersion parameters which can be used to 
    simulate an RNA-Seq dataset with synthetic 
    genes by drawing from a negative binomial 
    distribution. This function works in the 
    same way as described in (Soneson and 
    Delorenzi, BMC Bioinformatics, 2013) and 
    (Robles et al., BMC Genomics, 2012).
}
\details{
    Regarding \code{realCounts}, the file should strictly
    contain a unique gene name (e.g. Ensembl accession) in 
    the first column and all other columns should contain 
    read counts for each gene. Each column must be named
    with a unique sample identifier. See examples in the 
    ReCount database 
    \url{http://bowtie-bio.sourceforge.net/recount/}.
    
    Also, the parameter estimation involves a lot of 
    random sampling. For guaranteed reproducibility,
    be sure to use \code{set.seed} prior to any
    calculations. By default, when the metaseqR2 package
    is loaded, the seed is set to \code{42}.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
parList <- estimateSimParams(dataMatrix,libsizeGt=3e+4)
}
\author{
    Panagiotis Moulos
}

\name{normalizeDeseq}
\alias{normalizeDeseq}
\title{Normalization based on the DESeq package}
\usage{
    normalizeDeseq(geneCounts, sampleList,
        normArgs = NULL, output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of DESeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"deseq")} for
    an example and how you can modify it.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the DESeq native S4 object
    (CountDataSet). In the latter case it should be handled
    with suitable DESeq methods.}
}
\value{
    A matrix or a CountDataSet with normalized counts.
}
\description{
    This function is a wrapper over DESeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeDeseq(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}
\name{diagplotMds}
\alias{diagplotMds}
\title{Multi-Dimensinal Scale plots or RNA-Seq samples}
\usage{
    diagplotMds(x, sampleList, method = "spearman",
        logIt = TRUE, output = "x11", path = NULL, ...)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{method}{which correlation method to use. Same as
    the method parameter in \code{\link{cor}} function.}

    \item{logIt}{whether to log transform the values of x or
    not.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"}, \code{"ps"} or \code{"json"}. The latter is
    currently available for the creation of interactive
    volcano plots only when reporting the output, through the
    highcharts javascript library.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filename of the MDS plot produced if it's a file.
}
\description{
    Creates a Multi-Dimensional Scale plot for the given
    samples based on the count data matrix. MDS plots are
    very useful for quality control as you can easily see of
    samples of the same groups are clustered together based
    on the whole dataset.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(5000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotMds(dataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}

\docType{data}
\name{sampleListMm9}
\alias{sampleListMm9}
\title{Mouse RNA-Seq data with two conditions, four samples}
\format{a named \code{list} with condition and sample names.}
\source{
    ENCODE (http://genome.ucsc.edu/encode/)
}
\description{
    The sample list for \code{mm9GeneCounts}. See the data
    set description.
}
\author{
    Panagiotis Moulos
}
\keyword{datasets}
\name{metaseqrPlot}
\alias{metaseqrPlot}
\title{Diagnostic plots for the metaseqR2 package}
\usage{
    metaseqrPlot(object, sampleList, annotation = NULL,
        contrastList = NULL, pList = NULL,
        thresholds = list(p = 0.05, f = 1),
        plotType = c("mds", "biodetection", "countsbio",
            "saturation", "readnoise", "rnacomp", "correl",
            "pairs", "boxplot", "gcbias", "lengthbias",
            "meandiff", "meanvar", "deheatmap", "volcano",
            "biodist", "filtered", "mastat", "deregulogram",
            "statvenn", "foldvenn"),
        isNorm = FALSE, output = "x11", path = NULL, ...)
}
\arguments{
    \item{object}{a matrix or a data frame containing count
    data derived before or after the normalization procedure,
    filtered or not by the metaseqR2's filters and/or p-value.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{annotation}{a data frame containing annotation
    elements for each row in object. See also Details.}

    \item{contrastList}{the vector of contrasts as defined in 
    the main help page of \code{\link{metaseqr2}}.}

    \item{pList}{a list of p-values for each contrast as
    obtained from any of the \code{stat*} methods of the
    metaseqr package. See also Details.}

    \item{thresholds}{a list with the elements \code{"p"} and
    \code{"f"} which are the p-value and the fold change
    cutoff when \code{diagplotType="volcano"}.}

    \item{plotType}{one or more of the diagnostic plots
    supported in metaseqR2 package. See also Details.}

    \item{isNorm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{output}{one or more R plotting device to direct the
    plot result to. See Details.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    A named list containing the file names of the produced
    plots. Each list member is names according to the
    selected plotting device and is also a named list, whose
    names are the plot types. The final contents are the file
    names in case the plots are written to a physical
    location (not meaningful for \code{"x11"}).
}
\description{
    This is the main function for producing sructured quality
    control and informative graphs base on the results of the
    various steps of the metaseqR package. The graphs
    produced span a variety of issues like good sample
    reproducibility (Multi-Dimensional Scaling plot, biotype
    detection, heatmaps. diagplotMetaseqr, apart from
    implementing certain package-specific plots, is a wrapper
    around several diagnostic plots present in other RNA-Seq
    analysis packages such as EDASeq and NOISeq.
}
\details{
    Regarding \code{object}, the object can be fed to any of 
    the \code{diagplotMetaseqr} plotting systems but not every
    plot is meaningful. For example, it's meaningless to
    create a \code{"biodist"} plot for a count matrix before
    normalization or statistical testing.
    
    Regarding \code{annotation}, usually, it is a subset of 
    the annotation obtained by \code{\link{getAnnotation}} or
    a subset of possibly embedded annotation with the input
    counts table. This parameter is optional and required
    only when diagplotType is any of \code{"biodetection"},
    \code{"countsbio"}, \code{"saturation"},
    \code{"rnacomp"}, \code{"readnoise"}, \code{"biodist"},
    \code{"gcbias"}, \code{"lengthbias"} or
    \code{"filtered"}.
    
    Regarding \code{contrastList}, this parameter is optional 
    and required only when \code{diagplotType} is any of
    \code{"deheatmap"}, \code{"volcano"} or \code{"biodist"}.
    It can also be a named structured list of contrasts as 
    returned by the internal function 
    \code{metaseqR2:::makeContrastList}.
    
    Regarding \code{diagplotType}, many of these plots
    require the presence of additional package, something
    that is checked while running the main metaseqr2 function.
    The supported plots are \code{"mds"}, \code{"biodetection"},
    \code{"countsbio"}, \code{"saturation"}, \code{"rnacomp"},
    \code{"boxplot"}, \code{"gcbias"}, \code{"lengthbias"},
    \code{"meandiff"}, \code{"meanvar"}, \code{"deheatmap"},
    \code{"volcano"}, \code{"biodist"}, \code{"filtered"}, 
    \code{"readnoise"}, \code{"venn"}, \code{"correl"},
    \code{"pairwise"}. For a brief description of these plots 
    please see the main \code{\link{metaseqr2}} help page.
    
    Regarding \code{pList}, this parameter is optional 
    and required only when \code{diagplotType} is any of
    \code{"deheatmap"}, \code{"volcano"} or \code{"biodist"}.
    
    Regarding \code{output}, supported mechanisms are: 
    \code{"png"}, \code{"jpg"}, \code{"bmp"}, \code{"pdf"},
    \code{"ps"} or \code{"json"}. The latter is currently 
    available for the creation of interactive volcano plots 
    only when reporting the output, through the highcharts 
    javascript library. The default plotting (\code{"x11"}) 
    is not supported due to instability in certain devices.
}
\note{
    In order to make the best out of this function, you
    should generally provide the annotation argument as most
    and also the most informative plots depend on this. If
    you don't know what is inside your counts table or how
    many annotation elements you can provide by embedding it,
    it's always best to setup a local databse so as to use
    predefined annotations that work better with the
    functions of the whole package.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
metaseqrPlot(dataMatrix,sampleList,plotType=c("mds","boxplot"))

normArgs <- getDefaults("normalization","deseq2")
object <- normalizeDeseq2(dataMatrix,sampleList,normArgs)
metaseqrPlot(object,sampleList,plotType="boxplot")

## More
#p <- statDeseq2(object,sampleList)
#metaseqrPlot(object,sampleList,contrastList=contrast,pList=p,
#    plotType="volcano")
}
\author{
    Panagiotis Moulos
}
\name{metaseqr2}
\alias{metaseqr2}
\title{The main metaseqr2 pipeline}
\usage{
    metaseqr2(counts, sampleList, excludeList = NULL,
        fileType = c("auto", "sam", "bam", "bed"),
        path = NULL, contrast = NULL, libsizeList = NULL,
        embedCols = list(idCol = 4, gcCol = NA, nameCol = NA, 
            btCol = NA),
        annotation = NULL, org = c("hg18", "hg19", "hg38", "mm9", 
            "mm10", "rn5", "rn6", "dm3", "dm6", "danrer7", 
            "pantro4", "susscr3", "tair10", "equcab2"),
        refdb = c("ensembl", "ucsc", "refseq"), version = "auto",
        transLevel = c("gene", "transcript", "exon"),
        countType = c("gene", "exon","utr"),
        utrOpts = list(frac = 1, minLength = 300, downstream = 50),
        exonFilters = list(minActiveExons = list(exonsPerGene = 5, 
            minExons = 2, frac = 1/5)),
        geneFilters = list(length = list(length = 500), 
            avgReads = list(averagePerBp = 100, quantile = 0.25), 
            expression = list(median = TRUE, mean = FALSE,  
                quantile = NA, known = NA, custom = NA), 
            biotype = getDefaults("biotypeFilter", org[1]),
            presence = list(frac = 0.25, minCount = 10,
                perCondition = FALSE)),
        whenApplyFilter = c("postnorm", "prenorm"),
        normalization = c("deseq", "deseq2", "edaseq", "edger",
            "noiseq", "nbpseq", "absseq", "dss", "each", "none"),
        normArgs = NULL,
        statistics = c("deseq", "deseq2", "edger", "noiseq", 
            "bayseq", "limma", "nbpseq", "absseq", "dss"),
        statArgs = NULL,
        adjustMethod = sort(c(p.adjust.methods, "qvalue")),
        metaP = if (length(statistics) > 1) c("simes", 
            "bonferroni", "fisher", "dperm_min", "dperm_max",
            "dperm_weight", "fperm", "whitlock", "minp", "maxp",
            "weight", "pandora", "none") else "none",
        weight = rep(1/length(statistics), length(statistics)),
        nperm = 10000, pcut = NA, logOffset = 1, pOffset = NULL, 
        preset = NULL, qcPlots = c("mds", "biodetection", 
            "countsbio", "saturation", "readnoise","filtered", 
            "correl", "pairwise", "boxplot", "gcbias", 
            "lengthbias", "meandiff", "meanvar", "rnacomp", 
            "deheatmap", "volcano", "biodist", "mastat", 
            "statvenn", "foldvenn", "deregulogram"),
        figFormat = c("png", "jpg", "tiff", "bmp", "pdf", "ps"),
        outList = FALSE, exportWhere = NA,
        exportWhat = c("annotation", "p_value", "adj_p_value", 
            "meta_p_value", "adj_meta_p_value", "fold_change", 
            "stats", "counts","flags"),
        exportScale = c("natural", "log2", "log10", "vst", 
            "rpgm"),
        exportValues = c("raw", "normalized"),
        exportStats = c("mean", "median", "sd", "mad", "cv", 
            "rcv"),
        exportCountsTable = FALSE,
        restrictCores = 0.6, report = TRUE, reportTop = 0.1,
        reportTemplate = "default", saveGeneModel = TRUE,
        verbose = TRUE, runLog = TRUE, 
        reportDb = c("dexie", "sqlite"), 
        localDb = file.path(system.file(package = "metaseqR2"),
            "annotation.sqlite"),
        offlineReport = TRUE, 
        createTracks = FALSE, overwriteTracks = FALSE, 
        trackExportPath = file.path(exportWhere, "tracks"),
        trackInfo = list(stranded = FALSE, normTo = 1e+9, 
            urlBase = "http://www.trackserver.me",
            fasta = NULL, gtf = NULL,
            hubInfo = list(name = "MyHub", shortLabel = "My hub",
            longLabel = "My hub long", 
            email = "someone@example.com")), .progressFun = NULL,
            .exportR2C = FALSE,...)
}
\arguments{
    \item{counts}{a text tab-delimited file containing gene,
    exon or 3'UTR counts in one of the following formats: i) the
    first column contains unique gene or exon identifiers and
    the rest of the columns contain the read counts for each
    sample. ii) The first n columns should contain only **gene**
    annotation elements like chromosomal locations, gene
    accessions, exon accessions, GC content etc. and the rest
    columns should contain gene read counts, iii) \code{counts} 
    can also be an \code{.RData} file with previous analysis 
    elements (see Details) and iv) counts can be a list 
    representing the gene model (see Details). Several 
    restrictions apply im each of the four cases. See Details 
    for analytical descriptions.}

    \item{sampleList}{a list containing condition names and
    the samples under each condition or a small tab-delimited 
    file with the experiment description. Not needed when
    restoring a previous analysis. See Details for analytical 
    description.}
    
    \item{excludeList}{ a list of samples to exclude, in the 
    same (list) format as \code{sampleList} above.}

    \item{path}{an optional path where all the BED/BAM files
    are placed, to be prepended to the BAM/BED file names in
    the targets file. See Details for further information.}

    \item{fileType}{the type of raw input files. It can be
    \code{"auto"} for auto-guessing, \code{"bed"} for BED
    files, \code{"sam"} for SAM files or \code{"bam"} for BAM
    files.}

    \item{contrast}{a character vector of contrasts to be
    tested in the statistical testing step(s) of the pipeline.
    Each element of contrast should STRICTLY have the format
    "ConditionA_vs_ConditionB_vs_...". Special attention is
    needed as fold change calculations are based on this
    argument. If it is \code{NULL}, no statistical testing
    of fold change calculations are performed. See Details 
    for further information.}

    \item{libsizeList}{an optional named list where names
    represent samples (MUST be the same as the samples in
    \code{sample.list}) and members are the library sizes
    (the sequencing depth) for each sample. For example
    \code{libsize.list <- list(Sample_A1=32456913,}
    \code{Sample_A2=4346818)}.}
    
    \item{embedCols}{a named list with column numbers to 
    guide the case of embedded annotation. See Details for
    further information.}

    \item{annotation}{It can be one of i) \code{NULL} 
    (default) to use the existing annotation database 
    or fetch on the fly, ii) \code{"embedded"} if the
    annotation elements are embedded in the read counts file
    (restrictions apply) or iii) a list with a path to a GTF
    file and certain required metadata. See Details for a full
    description.}
    
    \item{org}{the supported organisms by \code{metaseqr2} or
    a user-named organism which has been imported to the 
    database. See Details for more information.} 

    \item{refdb}{the reference annotation repository from
    which to retrieve annotation elements to use with
    metaseqr2. It can be one of \code{"ensembl"} (default),
    \code{"ucsc"} or \code{"refseq"} or a user based one 
    (similar to the \code{org} argument).}
    
    \item{version}{the version of the annotation to use. See
    Details.}
    
    \item{transLevel}{perform differential expression 
    analysis at which transcriptional unit, can be one of
    \code{"gene"} (default), \code{"transcript"} for 
    reporting differential expression at the transcript
    level or \code{"exon"} for exon level.}
    
    \item{countType}{the type of reads inside the counts
    file. It can be one of \code{"gene"}, \code{"exon"} or
    \code{"utr"} for Quant-Seq (Lexogen) protocol. This is
    a very important and mandatory parameter as it defines
    the course of the workflow.}
    
    \item{utrOpts}{a named list with members \code{frac} 
    which is the fraction (0-1) of the 3' UTR region to count 
    reads in, \code{minLength} the minimum acceptable 3'UTR
    length irrespective of \code{frac} and \code{downstream} 
    the number of base pairs to flank the end of the 3' UTR of 
    transcripts when analyzing Quant-Seq data.}

    \item{exonFilters}{a named list whose names are the
    names of the supported exon filters and its members the
    filter parameters. See section "Exon filters" below for
    details.}

    \item{geneFilters}{a named list whose names are the
    names of the supported gene filters and its members the
    filter parameters. See section "Gene filters" below for
    details.}

    \item{whenApplyFilter}{a character string determining
    when to apply the exon and/or gene filters, relative to
    normalization. It can be \code{"prenorm"} to apply apply
    the filters and exclude genes from further processing
    before normalization, or \code{"postnorm"} to apply the
    filters after normalization (default). See also Details.}

    \item{normalization}{the normalization algorithm to be
    applied on the count data. It can be one of
    \code{"edaseq"} for EDASeq normalization, \code{"deseq"}
    for the normalization algorithm in the DESq package 
    (default), \code{"edger"} for the normalization algorithms 
    present in the edgeR package \code{"noiseq"} for the 
    normalization algorithms present in the NOISeq package 
    \code{"nbpseq"} for the normalization algorithms present 
    in the NBPSeq package or \code{"none"} to not normalize
    the data (highly unrecommended). Algorithm specific 
    arguments can be passed through the \code{normArgs} 
    argument).}

    \item{normArgs}{a named list whose names are the names
    of the normalization algorithm parameters and its members
    parameter values. See section "Normalization parameters"
    below for details. Leave \code{NULL} for the defaults of
    \code{normalization}.}

    \item{statistics}{one or more statistical analyses to be
    performed by the metaseqr2 pipeline. It can be one or more
    of \code{"deseq"} (default) to conduct statistical
    test(s) implemented in the DESeq package, \code{"edger"}
    to conduct statistical test(s) implemented in the edgeR
    package, \code{"limma"} to conduct the RNA-Seq version of
    statistical test(s) implemented in the limma package,
    \code{"noiseq"} to conduct statistical test(s)
    implemented in the NOISeq package, \code{"bayseq"} to
    conduct statistical test(s) implemented in the baySeq
    package, \code{"nbpseq"} to conduct statistical test(s)
    implemented in the NBPSeq package, \code{"deseq2"} to 
    conduct statistical test(s) implemented in the DESeq2 
    package, \code{"dss"} to conduct statistical test(s) 
    implemented in the DSS package and \code{"absseq"} to 
    conduct statistical test(s) implemented in the ABSSeq 
    package. In any case individual algorithm parameters are 
    controlled by the contents of the \code{statArgs} list.
    Finally, it can be \code{NA}. In this case no testing
    is performed and only fold changes are provided if
    \code{contrast} is not \code{NULL}.}

    \item{statArgs}{a named list whose names are the names
    of the statistical algorithms used in the pipeline. Each
    member is another named list whose names are the
    algorithm parameters and its members are the parameter
    values. See section "Statistics parameters" below for
    details. Leave \code{NULL} for the defaults of
    \code{statistics}.}

    \item{adjustMethod}{the multiple testing p-value
    adjustment method. It can be one of
    \code{\link{p.adjust.methods}} or \code{"qvalue"} from
    the qvalue Bioconductor package. Defaults to \code{"BH"}
    for Benjamini-Hochberg correction.}

    \item{metaP}{the meta-analysis method to combine
    p-values from multiple statistical tests . It can be 
    one of \code{"simes"} (default), \code{"bonferroni"}, 
    \code{"minp"}, \code{"maxp"}, \code{"weight"}, 
    \code{"pandora"}, \code{"dperm_min"}, \code{"dperm_max"}, 
    \code{"dperm_weight"}, \code{"fisher"}, \code{"fperm"}, 
    \code{"whitlock"} or \code{"none"}. See Details for a full
    description.}

    \item{weight}{a vector of weights with the same length as
    the \code{statistics} vector containing a weight for each
    statistical test. It should sum to 1. \strong{Use with
    caution with the} \code{dperm_weight} \strong{parameter!
    Theoretical background is not yet} \strong{solid and only
    experience shows improved results!}}

    \item{nperm}{the number of permutations performed to
    derive the meta p-value when \code{metaP="fperm"} or
    \code{metaP="dperm"}. It defaults to 10000.}

    \item{pcut}{a p-value cutoff for exporting differentially
    genes, default is to export all the non-filtered genes.}

    \item{logOffset}{an offset to be added to values during
    logarithmic transformations in order to avoid Infinity
    (default is \code{1}).}
    
    \item{pOffset}{a value between \code{0} and \code{1} 
    to multiply potential zero p-values with for the 
    combination methods including weighting or \code{NULL} 
    (default). See also Details.}

    \item{preset}{an analysis strictness preset.
    \code{preset} can be one of \code{"all_basic"},
    \code{"all_normal"}, \code{"all_full"},
    \code{"medium_basic"}, \code{"medium_normal"},
    \code{"medium_full"}, \code{"strict_basic"},
    \code{"strict_normal"} or \code{"strict_full"}, each of
    which control the strictness of the analysis and the
    amount of data to be exported. For an explanation of the
    presets, see the section "Presets" below.}

    \item{qcPlots}{a set of diagnostic plots to show/create.
    It can be one or more of \code{"mds"},
    \code{"biodetection"}, \code{"rnacomp"},
    \code{"countsbio"}, \code{"saturation"},
    \code{"readnoise"}, \code{"filtered"}, \code{"boxplot"},
    \code{"gcbias"}, \code{"lengthbias"}, \code{"meandiff"},
    \code{"meanvar"}, \code{"deheatmap"}, \code{"volcano"},
    \code{"mastat"}, \code{"biodist"}, \code{"statvenn"},
    \code{"foldvenn"}. See also Details.}

    \item{figFormat}{the format of the output diagnostic
    plots. It can be one or more of \code{"png"},
    \code{"jpg"}, \code{"tiff"}, \code{"bmp"}, \code{"pdf"},
    \code{"ps"}. The native format \code{"x11"} (for direct
    display) is not provided as an option as it may not
    render the proper display of some diagnostic plots in
    some devices.}

    \item{outList}{a logical controlling whether to export a
    list with the results in the running environment.}

    \item{exportWhere}{an output directory for the project
    results (report, lists, diagnostic plots etc.)}

    \item{exportWhat}{the content of the final lists. It can
    be one or more of \code{"annotation"}, to bind the
    annoation elements for each gene, \code{"p_value"}, to
    bind the p-values of each method, \code{"adj_p_value"},
    to bind the multiple testing adjusted p-values,
    \code{"meta_p_value"}, to bind the combined p-value from
    the meta-analysis, \code{"adj_meta_p_value"}, to bind the
    corrected combined p-value from the meta-analysis,
    \code{"fold_change"}, to bind the fold changes of each
    requested contrast, \code{"stats"}, to bind several
    statistics calclulated on raw and normalized counts (see
    the \code{exportStats} argument), \code{"counts"}, to
    bind the raw and normalized counts for each sample.}

    \item{exportScale}{export values from one or more
    transformations applied to the data. It can be one or
    more of \code{"natural"}, \code{"log2"}, \code{"log10"},
    \code{"vst"} (Variance Stabilizing Transormation, see the
    documentation of DESeq package) and \code{"rpgm"} which
    is ratio of mapped reads per gene model (either the gene
    length or the sum of exon lengths, depending on 
    \code{countType} argument). Note that this is not RPKM
    as reads are already normalized for library size using 
    one of the supported normalization methods. Also, 
    \code{"rpgm"} might be misleading when \code{normalization} 
    is other than \code{"deseq"}.}

    \item{exportValues}{It can be one or more of
    \code{"raw"} to export raw values (counts etc.) and
    \code{"normalized"} to export normalized counts.}

    \item{exportStats}{calculate and export several
    statistics on raw and normalized counts, condition-wise.
    It can be one or more of \code{"mean"}, \code{"median"},
    \code{"sd"}, \code{"mad"}, \code{"cv"} for the
    Coefficient of Variation, \code{"rcv"} for a robust
    version of CV where the median and the MAD are used
    instead of the mean and the standard deviation.}

    \item{exportCountsTable}{exports also the calculated 
    read counts table when input is read from bam files 
    and exports also the normalized count table in all 
    cases. Defaults to \code{FALSE}.}

    \item{restrictCores}{in case of parallel execution of
    several subfunctions, the fraction of the available cores
    to use. In some cases if all available cores are used
    (\code{restrictCores=1} and the system does not have
    sufficient RAM, the pipeline running machine might
    significantly slow down.}

    \item{report}{a logical value controlling whether to
    produce a summary report or not. Defaults to
    \code{TRUE}.}

    \item{reportTop}{a fraction of top statistically 
    significant genes to append to the HTML report. This 
    helps in keeping the size of the report as small as 
    possible, as appending the total gene list might 
    create a huge HTML file. Users can always retrieve 
    the whole gene lists from the report links. Defaults 
    to \code{0.1} (top 10% of statistically significant 
    genes). Set to \code{NA} or \code{NULL} to append all 
    the statistically significant genes to the HTML report.}

    \item{reportTemplate}{an HTML template to use for the
    report. Do not change this unless you know what you are
    doing.}

    \item{saveGeneModel}{in case of exon analysis, a list
    with exon counts for each gene will be saved to the file
    \code{exportWhere/data/gene_model.RData}. This file can
    be used as input to metaseqR for exon count based 
    analysis, in order to avoid the time consuming step of 
    assembling the counts for each gene from its exons}

    \item{verbose}{print informative messages during
    execution? Defaults to \code{TRUE}.}

    \item{runLog}{write a log file of the \code{metaseqr2}
    run using package log4r. Defaults to \code{TRUE}. The
    filename will be auto-generated.}
    
    \item{reportDb}{database system to use for storing the
    report intereactive graphs. Can be \code{"sqlite"} (default)
    or \code{"dexie"}. See Details for further explanation on
    what should be used.}
    
    \item{localDb}{the metaseqR2 annotaation database location.
    See also \code{link{buildAnnotationDatabase}}.}
    
    \item{offlineReport}{\code{TRUE} (default) to download
    and include the required JavaScript libraries to
    properly view the report offline. Ignored if
    \code{report=FALSE}}.
    
    \item{createTracks}{option to create normalized bigWig 
    files to display in a genome browser (e.g. UCSC). Defaults
    to \code{FALSE}.}
    
    \item{overwriteTracks}{overwrite tracks if they already
    exist? Defaults to \code{FALSE}.}
    
    \item{trackExportPath}{where to export the bigWig files, 
    defaults to \code{file.path(exportWhere,"tracks")}.}
    
    \item{trackInfo}{if \code{createTracks=TRUE}, a list with
    additional required information to create the tracks. See
    Details for further explanation.}
    
    \item{.progressFun}{a function which updates a 
    \code{Progress} object from shiny. This function must
    accept a \code{detail} argument. See 
    http://shiny.rstudio.com/articles/progress.html}
    
    \item{.exportR2C}{export additional RData along with
    \code{saveGeneModel}.}

    \item{...}{further arguments that may be passed to
    plotting functions, related to \code{\link{par}}.}
}
\value{
    If \code{outList} is \code{TRUE}, a named list whose
    length is the same as the number of requested contrasts.
    Each list member is named according to the corresponding
    contrast and contains a data frame of differentially
    expressed genes for that contrast. The contents of the
    data frame are defined by the \code{exportWhat,
    exportScale, exportStats, exportValues} parameters. If
    \code{report} is \code{TRUE}, the output list contains
    two main elements. The first is described above (the
    analysis results) and the second contains the same
    results but in HTML formatted tables.
}
\description{
    This function is the main metaseqr2 workhorse and
    implements the main metaseqr2 workflow which performs data
    read, filtering, normalization and statistical selection,
    creates diagnostic plots and exports the results and a
    report if requested. The metaseqr2 function is responsible
    for assembling all the steps of the metaseqr2 pipeline
    which i) reads the input gene or exon read count table
    ii) performs prelimininary filtering of data by removing
    chrM and other non-essential information for a typical
    differential gene expression analysis as well as a
    preliminary expression filtering based on the exon
    counts, if an exon read count file is provided. iii)
    performs data normalization with one of currently widely
    used algorithms, including EDASeq (Risso et al., 2011),
    DESeq (Anders and Huber, 2010), edgeR (Robinson et al.,
    2010), NOISeq (Tarazona et al., 2012) or no normalization
    iv) performs a second stage of filtering based on the
    normalized gene expression according to several gene
    filters v) performs statistical testing with one or more
    of currently widely used algorithms, including DESeq
    (Anders and Huber, 2010), edgeR (Robinson et al., 2010),
    NOISeq (Tarazona et al., 2012), limma (Smyth et al.,
    2005) for RNA-Seq data, baySeq (Hardcastle et al., 2012)
    vi) in the case of multiple statistical testing
    algorithms, performs meta-analysis using one of five
    available methods (see the meta.p argument) vii) exports
    the resulting differentially expressed gene list in text
    tab-delimited format viii) creates a set of diagnostic
    plots either available in the aforementioned packages or
    metaseqr2 specific ones and ix) creates a comprehensive
    HTML report which summarizes the run information, the
    results and the diagnostic plots. Certain diagnostic
    plots (e.g. the volcano plot) can be interactive with the
    use of the external Highcharts
    (http://www.highcharts.com) JavaScript library for
    interactive graphs. Although the inputs to the metaseqr2
    workflow are many, in practice, setting only very few of
    them and accepting the defaults as the rest can result in
    quite comprehensible results for mainstream organisms
    like mouse, human, fly and rat.
}
\details{
    When \code{counts} is a tab-delimited file, the following
    restrictions apply:
    \itemize{ 
        \item In the case (i) the first cell of each row is a 
        gene or exon accession and the rest are integers 
        representing the counts for that accession. In that 
        case, the \code{annotation} parameter should strictly 
        be \code{NULL} or an external file in GTF format.
        \item In the case (ii) the \code{annotation} parameter 
        can also be \code{"embedded"}. The ideal embedded 
        annotation contains 8 columns, chromosome, gene or 
        exon start, gene or exon end, gene or exon accession, 
        GC-content (fraction or percentage), strand, HUGO 
        gene symbol and gene biotype (e.g. "protein_coding" 
        or "ncRNA"). When the \code{annotation} parameter is 
        \code{"embedded"}, certain of these features are 
        mandatory (co-ordinates and accessions). If they are 
        not present, the pipeline will not run. If additional 
        elements are not present (e.g. GC content or biotypes), 
        certain features of metaseqr2 will not be available. 
        For example, EDASeq normalization will not be 
        performed based on a GC content covariate but based 
        on gene length which is not what the authors of
        EDASeq suggest. If biotypes are not present, a lot of
        diagnostic plots will not be available. If the HUGO 
        gene symbols are missing, the final annotation will 
        contain only gene accessions and thus be less 
        comprehensible. Counts can be a data frame satisfying 
        the above conditions. It is a data frame by default 
        when \code{read2count} is used.
        \item In the case (iii)  the .RData file 
        (output of \code{\link{save}} function contains 
        static input elements (list containing the gene 
        model (exon counts for each gene), gene and exon 
        annotation to avoid re-(down)loading and/or gene 
        counts depending on \code{countType}). This kind of 
        input facilitates the re-analysis of the same 
        experiment, using different filtering, normalization 
        and statistical algorithms. This \code{.RData} file 
        is produced when \code{saveGeneModel=TRUE}.
        \item In the case (iv) \code{counts} can be a list 
        representing the gene model (exon/UTR counts for each 
        gene). This \code{.RData} file can be generated by 
        setting \code{saveGeneModel=TRUE} when performing data 
        analysis for the first time.
    }
    
    Regarding \code{sampleList} it should have the format 
    \code{sampleList <-} \code{list(ConditionA=c("Sample_A1",} 
    \code{"Sample_A2","Sample_A3"),} 
    \code{ConditionB=c("Sample_B1","Sample_B2"),} 
    \code{ConditionC=c("Sample_C1","Sample_C2"))}. The names of 
    the samples in list members MUST match the column names 
    containing the read counts in the counts file. If they do 
    not match, the pipeline will either crash or at best, ignore 
    several of your samples. Alternative, \code{sampleList} can 
    be a small tab-delimited file structured as follows: the 
    first line of the external tab delimited file should contain
    column names (names are not important). The first column MUST
    contain UNIQUE sample names and the second column MUST
    contain the biological condition where each of the samples
    in the first column should belong to. If the \code{counts}
    argument is missing, the \code{sampleList} argument MUST be 
    a targets text tab-delimited file which contains the sample 
    names, the BAM/BED file names and the biological 
    conditions/groups for each sample/file. The file should be 
    text tab-delimited and structured as follows: the first line 
    of the external tab delimited file should contain column 
    names (names are not important). The first column MUST 
    contain UNIQUE sample names. The second column MUST contain
    the raw BAM/BED files WITH their full path. Alternatively, 
    the \code{path} argument should be provided. If \code{path} 
    is not provided and if the files in the second column of 
    the targets file do not contain a path to a directory, the 
    current directory is assumed to be the BAM/BED file 
    container. The third column MUST contain the biological 
    condition where each of the samples in the first column 
    should belong to.
    
    Regarding \code{contrast}, a valid example based on the 
    \code{sampleList} above is 
    \code{contrast <- c("ConditionA_vs_ConditionB",}
    \code{"ConditionA_vs_ConditionC",}
    \code{"ConditionA_vs_ConditionB_vs_ConditionC")}. The
    first element of pairwise contrasts (e.g. "ConditionA"
    above) MUST be the control condition or any reference
    that ConditionB is checked against. metaseqr2 uses this
    convention to properly calculate fold changes.
    
    Regarding \code{embedCols}, this list must contain four
    members, named \code{idCol}, \code{gcCol}, \code{nameCol}
    and \code{btCol}, which hold the position in the delimited
    file with embedded annotation, where unique gene ids, GC
    content, gene names and gene biotypes respectively are
    located. More specifically:
    \itemize{
        \item \code{idCol} is an integer denoting the column 
        number in the file (or data frame) provided with the 
        counts argument, where the unique gene accessions are.
        Default to \code{4} which is the standard feature name
        column in a BED file.
        \item \code{gcCol} is an integer denoting the column 
        number in the file (or data frame) provided with the 
        \code{counts} argument, where each gene's GC content 
        is given. If not provided, GC content normalization 
        provided by EDASeq will not be available.
        \item \code{nameCol} is an integer denoting the column 
        number in the file (or data frame) provided with the 
        counts argument, where the HUGO gene symbols are given. 
        If not provided, it will not be available when reporting
        results. In addition, the \code{"known"} gene filter 
        will not be available for application.
        \item \code{btCol} is an integer denoting the column 
        number in the file (or data frame) provided with the
        counts argument, where the gene biotypes are given. 
        If not provided, the \code{"biodetection"}, 
        \code{"countsbio"}, \code{"saturation"}, 
        \code{"filtered"} and \code{"biodist"} plots will not 
        be available.
    }
    
    Regarding annotation instructs metaseqr2 where to find the
    annotation for the given counts file. It can be one of i)
    \code{"download"} (default) for automatic downloading of
    the annotation for the organism specified by the org
    parameter (using biomaRt), ii) \code{"embedded"} if the
    annotation elements are embedded in the read counts file
    or iv) a file specified by the user which should be as
    similar as possible to the \code{"download"} case, in
    terms of column structure.
    
    Regarding \code{org}, it can be, for human genomes 
    \code{"hg18"}, \code{"hg19"} or \code{"hg38"}, for mouse 
    genomes \code{"mm9"}, \code{"mm10"}, for rat genomes 
    \code{"rn5"} or \code{"rn6"}, for drosophila genome
    \code{"dm3"} or \code{"dm6"}, for zebrafish genome 
    \code{"danrer7"}, \code{"danrer10"} or \code{"danrer11"}, 
    for chimpanzee genome \code{"pantro4"}, \code{"pantro5"}, 
    for pig genome \code{"susscr3"}, \code{"susscr11"}, for 
    Arabidopsis thaliana genome \code{"tair10"} and for 
    Equus caballus genome \code{"equcab2"}. Finally, it can
    be \code{"USER_NAMED_ORG"} with a custom organism which
    has been imported to the annotation database by the user
    using a GTF file. For example \code{org="mm10_p1"}.
    
    Regarding \code{version}, this is an integer denoting the
    version of the annotation to use from the local annotation
    database or fetch on the fly. For Ensembl, it corresponds
    to Ensembl releases, while for UCSC/RefSeq, it is the
    date of creation (locally).
    
    Regarding \code{whenApplyFilter}, in the case of
    \code{whenApplyFilter="prenorm"}, a first normalization
    round is applied to a copy of the gene counts matrix in
    order to derive the proper normalized values that will
    constitute the several expression-based filtering
    cutoffs.
    
    Regarding \code{metaP}, for the \code{"fisher"} and 
    \code{"fperm"} methods, see the documentation of the R 
    package MADAM. For the \code{"whitlock"} method, see the 
    documentation of the survcomp Bioconductor package. With 
    the \code{"maxp"} option, the final p-value is the maximum 
    p-value out of those returned by each statistical test. 
    This is equivalent to an "intersection" of the results 
    derived from each algorithm so as to have a final list with 
    the common genes returned by all statistical tests. 
    Similarly, when \code{meta.p="minp"}, is equivalent to a 
    "union" of the results derived from each algorithm so as 
    to have a final list with all the genes returned by all 
    statistical tests. The latter can be used as a very lose 
    statistical threshold to aggregate results from all methods 
    regardless of their False Positive Rate. With the 
    \code{"simes"} option, the method proposed by Simes 
    (Simes, R. J., 1986) is used. With the \code{"dperm_min"}, 
    \code{"dperm.max"}, \code{"dperm.weight"} options, a 
    permutation procedure is initialed, where \code{nperm} 
    permutations are performed across the samples of the 
    normalized counts matrix, producing \code{nperm} permuted 
    instances of the initital dataset. Then, all the chosen 
    statistical tests are re-executed for each permutation. 
    The final p-value is the number of times that the p-value 
    of the permuted datasets is smaller than the original 
    dataset. The p-value of the original dataset is created 
    based on the choice of one of \code{dperm.min}, 
    \code{dperm.max} or \code{dperm.weight} options. In case of 
    \code{dperm.min}, the intial p-value vector is consists 
    of the minimum p-value resulted from the applied 
    statistical tests for each gene. The maximum p-value 
    is used with the \code{dperm.max} option. With the 
    \code{dperm.weight} option, the \code{weight} 
    weighting vector for each statistical test is used to 
    weight each p-value according to the power of 
    statistical tests (some might work better for a 
    specific dataset). Be careful as the permutation 
    procedure usually requires a lot of time. However, 
    it should be the most accurate. This method will NOT 
    work when there are no replicated samples across 
    biological conditions. In that case, use 
    \code{meta.p="simes"} instead. Finally, there are the 
    \code{"minp"}, \code{"maxp"} and \code{"weight"}
    options which correspond to the latter three methods 
    but without permutations. Generally, permutations 
    would be accurate to use when the experiment includes
    >5 samples per condition (or even better 7-10) which 
    is rather rare in RNA-Seq experiments. Finally, 
    \code{"pandora"} is the same as \code{"weight"} and is
    added to be in accordance with the main algorithm name.
    
    Regarding \code{pOffset}, it is used to correct for
    the case of a p-value which is equal to 0 as a result
    of internal numerical and approximation procedures.
    When \code{NULL}, random numbers greater than 0 and
    less than or equal to 0.5 are used to multiply the
    offending p-values with the lowest provided non-zero
    p-value, maintaining thus a virtual order of 
    significance, avoiding having the same p-values for 
    two tests and assuming that all zero p-values represent
    extreme statistical significance. When a numeric
    between 0 and 1, this number is used for the above
    multiplication instead.
    
    Regarding \code{qcPlots} The \code{"mds"} stands
    for Mutlti-Dimensional Scaling and it creates a PCA-like
    plot but using the MDS dimensionality reduction instead.
    It has been succesfully used for NGS data (e.g. see the
    package htSeqTools) and it shows how well samples from
    the same condition cluster together. For
    \code{"biodetection"}, \code{"countsbio"},
    \code{"saturation"}, \code{"rnacomp"},
    \code{"readnoise"}, \code{"biodist"} see the vignette of
    NOISeq package. The \code{"saturation"} case has been
    rewritten in order to display more samples in a more
    simple way. In addition, the \code{"readnoise"} plots 
    represent an older version or the RNA composition plot 
    included in older versions of NOISeq. For \code{"gcbias"},
    \code{"lengthbias"}, \code{"meandiff"}, \code{"meanvar"} 
    see the vignette of EDASeq package. \code{"lenghtbias"} 
    is similar to \code{"gcbias"} but using the gene length 
    instead of the GC content as covariate. The \code{"boxplot"} 
    option draws boxplots of log2 transformed gene counts. The
    \code{"filtered"} option draws a 4-panel figure with the
    filtered genes per chromosome and per biotype, as
    absolute numbers and as fractions of the genome. See also
    the help page of \code{\link{diagplotFiltered}}. The
    \code{"deheatmap"} option performs hierarchical
    clustering and draws a heatmap of differentially
    expressed genes. In the context of diagnostic plots, it's
    useful to see if samples from the same groups cluster
    together after statistical testing. The \code{"volcano"}
    option draws a volcano plot for each contrast and if a
    report is requested, an interactive volcano plot is
    presented in the HTML report. The \code{"venn"} option
    will draw an up to 5-way Venn diagram depicting the
    common and specific to each statistical algorithm genes
    and for each contrast, when meta-analysis is performed.
    The \code{"correl"} option creates two correlation
    graphs: the first one is a correlation heatmap (a
    correlation matrix which depicts all the pairwise
    correlations between each pair of samples in the counts
    matrix is drawn as a clustered heatmap) and the second
    one is a correlogram plot, which summarizes the
    correlation matrix in the form of ellipses (for an
    explanation please see the vignette/documentation of the
    R package corrplot. Set \code{qcPlots=NULL} if you don't
    want any diagnostic plots created.
    
    Regarding \code{reportDb}, contrary with the first version of
    metaseqR, all graphs in the metaseqR2 report are interactive
    with the usage of the JavaScript libraries Highcharts, plotly
    (heatmaply) and jvenn.js. However, this adds a great burden
    regarding rendering the final HTML file and its content,
    a burden which becomes heavier by the fact the metaseqR2
    report is rendered using \code{knitr} and \code{rmarkdown} 
    instead of raw HTML (previously, \code{brew}). Therefore, the
    pre-calculated JSON objects representing the graphs are 
    stored either in a report-specific IndexedDB 
    (https://javascript.info/indexeddb) flavor called Dexie
    (https://dexie.org/) (default) or in an SQLite database
    and then queried using \code{sql.js} 
    (https://github.com/kripken/sql.js/). Dexie is 
    prefered because it is very efficient and can produce an 
    independent report that does not need to be served through 
    a web-server and can be viewed locally. Although Dexie is 
    very efficient, some caution is required as \code{knitr} and 
    \code{\link{render}} from \code{rmarkdown} are not very 
    memory efficient when rendering larger HTML files. A large 
    HTML file may be produced when analyzing a large dataset
    with a lot of contrasts that may result in a lot of tables.
    In such cases, if the report generation crashes with errors
    related to memory, try lowering the \code{reportTop}
    argument. \code{reportTop} does not affect the final lists
    of differentially expressed genes, only the report tables.
    The same must be applied also if the report takes too much
    time to load. If the report is to be served through a web
    server like Apache (e.g. when the report is provided by a 
    facility to end users), \code{reportDb="sqlite"} may be
    preferred as the total report size will be smaller because
    of an SQLite database hosting all plots which are queried
    when required but from the SQLite database and not from
    the in-browser database (Dexie). **Note** that when using
    an SQLite database, you will **NOT** be able to view the
    report in any browser other than Microsoft Edge because
    of security policies regarding local file access. **Note**
    also that \code{sql.js} is a rather large JavaScript library
    (around 2.5MB).
    
    Regarding \code{trackInfo}, it is a helper list to guide the
    bigWig track creation and has the following members:
    \itemize{
        \item \code{stranded}, which can be \code{TRUE} or
        \code{FALSE} depending on whether you wish to create
        stranded tracks by separating + and - strand reads.
        In the case of stranded tracks, a UCSC Genome Brower
        trackhub is created. Individual tracks can be retrieved
        from the trackhub.
        \item \code{normTo}, which is a large integer, denoting
        the total sum of signal to be used as the normalization
        target. It defaults to \code{1e+9}. This means that if 
        for a particular sample the sum of signal is 1.5e+9
        (\code{sum(sapply(coverage(x),sum)) == 1.5e+9}) then
        this is linearly scaled to 1e+9.
        \item \code{urlBase}, which is a base url appended to
        the bigWig files produced (the base path of the 
        \code{bigDataUrl} in UCSC Genome Browser track lines).
        \item \code{hubInfo}, a list with the track hub 
        description in case of stranded tracks. Please see the
        track hub specifications at the UCSC Genome Browser site.
        \item \code{fasta}, reference genome in FASTA format for 
        the case of analyzing a custom, non-directly supported 
        organism. It will be converted to the .2bit format and 
        written along with a track hub.
        \item \code{gtf}, a GTF file describing gene models in 
        the case of analyzing a custom, non-directly supported 
        organism. It will be converted to the .bigBed format and 
        written along with a track hub. Essentially the same as
        \code{annotation$gtf}.
    }
    All files (bigWig files, track/trackhub info) are written in 
    the \code{tracks} subdirectory of the main path where the
    report and the outputs are written.
    
}
\note{
    Currently only gene and exon annotation from Ensembl
    (http://www.ensembl.org), UCSC and RefSeq are supported. 
    In addition, the user may choose to use own GTF file on
    the fly or import to the backend annotation database (see
    \code{\link{buildAnnotationDatabase}}). Thus, the unique
    gene ids in the counts files should correspond to valid 
    Ensembl, UCSC or RefSeq gene or exon accessions for the 
    organism of interest, or according to the user's GTF. 
    If you are not sure about the source of your counts file or 
    do not know how to produce it, it's better to start from the 
    original BAM/BED files (metaseqr2 will use the
    \code{\link{read2count}} function to create a counts
    file). Keep in mind that in the case of BED files, the
    performance will be significantly lower and the overall
    running time significanlty higher as the R functions
    which are used to read BED files to proper structures
    (GenomicRanges) and calculate the counts are quite slow.
    An alternative way is maybe the easyRNASeq package
    (Delhomme et al, 2012). The \code{\link{read2count}}
    function does not use this package but rather makes use
    of standard Bioconductor functions to handle NGS data. If
    you wish to work outside R, you can work with other
    popular read counters such as the HTSeq read counter
    (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html).
    Please also note that in the current version, the members
    of the \code{geneFilters} and \code{exonFilters} lists
    are not checked for validity so be careful to supply with
    correct names otherwise the pipeline will crash or at the
    best case scenario, will ignore the filters. Also note
    that when you are supplying metaseqr2 wtih an exon counts
    table, gene annotation is always downloaded so please be
    sure to have a working internet connection. In addition
    to the above, if you have a multiple core system, be very
    careful on how you are using the \code{restrictCores}
    argument and generally how many cores you are using with
    scripts purely written in R. The analysis with exon read
    data can very easily cause memory problems, so unless you
    have more than 64Gb of RAM available, consider setting
    restrict.cores to something like 0.2 when working with
    exon data. Finally, if you do not wish to download the
    same annotation again and again when performing multiple
    analyses, it is best to use the
    \code{\link{buildAnnotationDatabase}} function to download 
    and store the resulting data frames in local SQLite 
    database and then use these files with the 
    \code{org}, \code{refdb} and \code{version} options.

    Please note that the \strong{meta-analysis} feature
    provided by metaseqr2 does not satisfy the strict definition 
    of "meta-analysis", which is the combination of multiple 
    similar datasets under the same statistical methodology. 
    Instead it is the use of mulitple statistical tests applied 
    to the same data. For the Simes method, please consult also
    "Simes, R. J. (1986). "An improved Bonferroni procedure
    for multiple tests of significance". Biometrika 73 (3):
    751–754."
    
    Also, if \code{weight="meta_perm"} ideally one would 
    want to create the same set of indices for a given 
    dataset so as to create reproducible p-values. To 
    achieve this, use the \code{set.seed} function prior to
    any calculations. When metaseqR2 is loaded, the random
    seed is set to \code{42}.
}
\section{Exon filters}{
    The exon filters are a set of filters which are applied
    after the gene models are assembled from the read counts
    of individual exons and before the gene expression is
    summarized from the exons belonging to each gene. These
    filters can be applied when the input read counts file
    contains exon reads. It is not applicable when the input
    file already contains gene counts. Such filters can be
    for example "accept genes where all the exons contain
    more than x reads" or "accept genes where there is read
    presence in at least m/n exons, n being the total exons
    of the gene". Such filters are NOT meant for detecting
    differential splicing as also the whole metaseqr2
    pipeline, thus they should not be used in that context.
    The \code{exonFilters} argument is a named list of
    filters, where the names are the filter names and the
    members are the filter parameters (named lists with
    parameter name, parameter value). See the usage of the
    \code{metaseqr2} function for an example of how these
    lists are structured. The supported exon filters in the
    current version are: i) \code{minActiveExons} which
    implements a filter for demanding m out of n exons of a
    gene to have a certain read presence with parameters
    \code{exonsPerGene}, \code{minExons} and \code{frac}.
    The filter is described as follows: if a gene has up to
    \code{exonsPerGene} exons, then read presence is
    required in at least \code{minExons} of them, else read
    presence is required in a \code{frac} fraction of the
    total exons. With the default values, the filter
    instructs that if a gene has up to 5 exons, read presence
    is required in at least 2, else in at least 20% of the
    exons, in order to be accepted. More filters will be
    implemented in future versions and users are encouraged
    to propose exon filter ideas to the author by mail. See
    \code{metaseqr2} usage for the defaults. Set
    \code{exonFilters=NULL} to not apply any exon filtering.
}

\section{Gene filters}{
    The gene filters are a set of filters applied to gene
    expression as this is manifested through the read
    presence on each gene and are preferably applied after
    normalization. These filters can be applied both when the
    input file or data frame contains exon read counts and
    gene read counts. Such filter can be for example "accept
    all genes above a certain count threshold" or "accept all
    genes with expression above the median of the normalized
    counts distribution" or "accept all with length above a
    certain threshold in kb" or "exclude the 'pseudogene'
    biotype from further analysis". The supported gene
    filters in the current version, which have the same
    structure as the exon filters (named list of lists with
    filter names, parameter names and parameter arguments)
    are: i) \code{length} which implements a length filter
    where genes are accepted for further analysis if they are
    above \code{length} (its parameter) kb. ii)
    \code{avg.reads} which implements a filter where a gene
    is accepted for further analysis if it has more average
    reads than the \code{quantile} of the average count
    distribution per \code{averagePerBp} base pairs. In
    summary, the reads of each gene are averaged per
    \code{averagePerBp} based on each gene's length (in
    case of exons, input the "gene's length" is the sum of
    the lengths of exons) and the \code{quantile} quantile of
    the average counts distribution is calculated for each
    sample. Genes passing the filter should have an average
    read count larger than the maximum of the vector of the
    quantiles calculated above. iii) \code{expression} which
    implements a filter based on the overall expression of a
    gene. The parameters of this filter are: \code{median},
    where genes below the median of the overall count
    distribution are not accepted for further analysis (this
    filter has been used to distinguish between "expressed"
    and "not expressed" genes in several cases, e.g. (Mokry
    et al., NAR, 2011) with a logical as value, \code{mean}
    which is the same as \code{median} but using the mean,
    \code{quantile} which is the same as the previous two but
    using a specific quantile of the total counts
    distribution, \code{known}, where in this case, a set of
    known not-expressed genes in the system under
    investigation are used to estimate an expression cutoff.
    This can be quite useful, as the genes are filtered based
    on a "true biological" cutoff instead of a statistical
    cutoff. The value of this filter is a character vector of
    HUGO gene symbols (MUST be contained in the annotation,
    thus it's better to use \code{annotation="download"})
    whose counts are used to build a "null" expression
    distribution. The 90th quantile of this distribution is
    then the expression cutoff. This filter can be combined
    with any other filter. Be careful with gene names as they
    are case sensitive and must match exactly ("Pten" is
    different from "PTEN"!). iv) \code{biotype} where in this
    case, genes with a certain biotype (MUST be contained in
    the annotation, thus it's better to use
    \code{annotation="download"}) are excluded from the
    analysis. This filter is a named list of logical, where
    names are the biotypes in each genome and values are
    \code{TRUE} or \code{FALSE}. If the biotype should be
    excluded, the value should be \code{TRUE} else
    \code{FALSE}. See the result of
    \code{get.defaults("biotype.filter","hg19")} for an
    example. Finally, in future versions there will be
    support for user-defined filters in the form of a
    function. v) \code{presence} where in this case, a gene
    is further considered for statistical testing if 
    \code{frac} (x100 for a percentage value) have more
    than \code{minCount} reads across all samples
    (\code{perCondition=FALSE}) or across the samples
    of each condition (\code{perCondition=TRUE}).
}

\section{Normalization parameters}{
    The normalization parameters are passed again as a named
    list where the names of the members are the normalization
    parameter names and the values are the normalization
    parameter values. You should check the documentation of
    the packages EDASeq, DESeq, edgeR, NOISeq and NBPSeq for
    the parameter names and parameter values. There are a few
    exceptions in parameter names: in case of
    \code{normalization="edaseq"} the only parameter names
    are \code{within.which} and \code{between.which},
    controlling the withing lane/sample and between
    lanes/samples normalization algorithm. In the case 
    of \code{normalization="nbpseq"}, there is one
    additional parameter called \code{main.method} which can
    take the calues \code{"nbpseq"} or \code{"nbsmyth"}.
    These values correspond to the two different workflows
    available in the NBPSeq package. Please, consult the
    NBPSeq package documentation for further details. For the
    rest of the algorithms, the parameter names are the same
    as the names used in the respective packages. For
    examples, please use the \code{\link{getDefaults}}
    function.
}

\section{Statistics parameters}{
    The statistics parameters as passed to statistical
    algorithms in metaseqr2, exactly with the same way as the
    normalization parametes above. In this case, there is one
    more layer in list nesting. Thus, \code{statArgs} is a
    named list whose names are the names the algorithms used
    (see the \code{statistics} parameter). Each member is
    another named list,with parameters to be used for each
    statistical algorithm. Again, the names of the member
    lists are parameter names and the values of the member
    lists are parameter values. You should check the
    documentations of DESeq, edgeR, NOISeq, baySeq, limma and
    NBPSeq for these parameters. There are a few exceptions
    in parameter names: In case of \code{statistics="edger"},
    apart from the rest of the edgeR statistical testing
    arguments, there is the argument \code{mainMethod} which
    can be either \code{"classic"} or \code{"glm"}, again
    defining whether the binomial test or GLMs will be used
    for statistical testing. For examples, please use the
    \code{\link{getDefaults}} function. When
    \code{statistics="nbpseq"}, apart from the rest arguments
    of the NBPSeq functions \code{estimate.disp} and
    \code{estimate.dispersion}, there is the argument
    \code{mainMethod} which can be \code{"nbpseq"} or
    \code{"nbsmyth"}. This argument determines the parameters
    to be used by the \code{estimate.dispersion} function or
    by the \code{estimate.disp} function to estimate RNA-Seq
    count dispersions. The difference between the two is that
    they constitute different starting points for the two
    workflows in the package NBPSeq. The first worklfow (with
    \code{mainMethod="nbpseq"} and the
    \code{estimate.dispersion} function is NBPSeq package
    specific, while the second (with
    \code{mainMethod="nbsmyth"} and the \code{estimate.disp}
    function is similar to the workflow of the edgeR package.
    For additional information regarding the statistical
    testing in NBPSeq, please consult the documentation of
    the NBPSeq package.
}

\section{Presets}{
    The analysis presets are a set of keywords (only one can
    be used) that predefine some of the parameters of the
    metaseqr2 pipeline. For the time being they are quite
    simple and they control i) the strictness of filtering
    and statistical thresholding with three basic levels
    ("all", "medium", "strict") and ii) the data columns that
    are exported, again in three basic ways ("basic",
    "normal", "full") controlling the amount of data to be
    exported. These keywords can be combined with a dot in
    the middle (e.g. \code{"all.basic"} to define an analysis
    preset. When using analysis presets, the following
    argumentsof metaseqr2 are overriden: \code{exonFilters},
    \code{geneFilters}, \code{pcut}, \code{exportWhat},
    \code{exportScale}, \code{exportValues},
    \code{exonStats}. If you want to explicitly control the
    above arguments, the \code{preset} argument should be set
    to \code{NULL} (default). Following is a synopsis of the
    different presets and the values of the arguments they
    moderate: 
    \itemize{ 
        \item \code{"all_basic"}: use all genes (do not filter)
        and export all genes and basic annotation and statistics 
        elements. In this case, the above described arguments become: 
        \itemize{ 
            \item \code{exonFilters=NULL} 
            \item \code{geneFilters=NULL}
            \item \code{pcut=1} 
            \item \code{exportWhat=c("annotation","p_value",}
            \code{"adj_p_value","meta_p_value",}
            \code{"adj_meta_p_value","fold_change")} 
            \item \code{exportScale=c("natural","log2")} 
            \item \code{exportValues=c("normalized")} 
            \item \code{exportStats=c("mean")} 
        } 
        \item \code{"all_normal"}: use all genes (do not filter) 
        and export all genes and normal annotation and statistics
        elements. In this case, the above described arguments
        become: 
        \itemize{ 
            \item \code{exonFilters=NULL} 
            \item \code{geneFilters=NULL} 
            \item \code{pcut=1} 
            \item \code{exportWhat=c("annotation","p_value",}
            \code{"adj_p_value","meta_p_value",} 
            \code{"adj_meta_p_value","fold_change","stats",}
            \code{"counts")}
            \item \code{exportScale=c("natural","log2")} 
            \item \code{exportValues=c("normalized")} 
            \item \code{exportStats=c("mean","sd","cv")} 
        }
        \item \code{"all_full"}: use all genes (do not filter) 
        and export all genes and all annotation and statistics
        elements. In this case, the above described arguments 
        become: 
        \itemize{ 
            \item \code{exonFilters=NULL} 
            \item \code{geneFilters=NULL}
            \item \code{pcut=1} 
            \item \code{exportWhat=c("annotation","p_value",
            "adj_p_value","meta_p_value",} 
            \code{"adj_meta_p_value","fold_change","stats","counts")}
            \item \code{exportScale=c("natural","log2","log10","vst")}
            \item \code{exportValues=c("raw","normalized")} 
            \item \code{exportStats=c("mean","median","sd","mad",
            "cv","rcv")}
    } 
    \item \code{"medium_basic"}: apply a medium set of filters
    and and export statistically significant genes and basic
    annotation and statistics elements. In this case, the above
    above described arguments become: 
    \itemize{
        \item \code{exonFilters=list(minActiveExons=}
        \code{list(exonsPerGene=5,minExons=2,frac=1/5))} 
        \item \code{geneFilters=list(length=list(length=500),} 
        \code{avgReads=list(averagePerBp=100,quantile=0.25),} 
        \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,}
        \code{known=NA,custom=NA),} 
        \code{biotype=getDefaults("biotypeFilter",org[1]))}
        \item \code{pcut=0.05} 
        \item \code{exportWhat=c("annotation","p_value",}
        \code{"adj_p_value","meta_p_value",} 
        \code{"adj_meta_p_value","fold_change")} 
        \item \code{exportScale=c("natural","log2")} 
        \item \code{exportValues=c("normalized")} 
        \item \code{exportStats=c("mean")} 
    } 
    \item \code{"medium_normal"}: apply a medium set of filters 
    and export statistically significant genes and normal
    annotation and statistics elements. In this case, the
    above described arguments become: 
    \itemize{ 
        \item \code{exonFilters=list(minActiveExons=}
        \code{list(exonsPerGene=5,minExons=2,frac=1/5))} 
        \item \code{geneFilters=list(length=list(length=500),} 
        \code{avgReads=list(averagePerBp=100,quantile=0.25),} 
        \code{expression=list(median=TRUE,mean=FALSE,}
        \code{quantile=NA,known=NA,custom=NA),} 
        \code{biotype=getDefaults("biotypeFilter",org[1]))}
        \item \code{pcut=0.05} 
        \item \code{export.what=c("annotation","p_value",}
        \code{"adj_p_value","meta_p_value",} 
        \code{"adj_meta_p_value","fold_change",}
        \code{"stats","counts")}
        \item \code{exportScale=c("natural","log2")} 
        \item \code{exportValues=c("normalized")} 
        \item \code{exportStats=c("mean","sd","cv")} 
    } 
    \item \code{"medium_full"}: apply a medium set of filters 
    and export statistically significant genes and full
    annotation and statistics elements. In this case, the
    above described arguments become: 
    \itemize{ 
        \item \code{exonFilters=list(minActiveExons=}
        \code{list(exonsPerGene=5,minExons=2,frac=1/5))} 
        \item \code{geneFilters=list(length=list(length=500),} 
        \code{avgReads=list(averagePerBp=100,quantile=0.25),} 
        \code{expression=list(median=TRUE,mean=FALSE,}
        \code{quantile=NA,known=NA,custom=NA),} 
        \code{biotype=getDefaults("biotypeFilter",org[1]))}
        \item \code{pcut=0.05} 
        \item \code{exportWhat=c("annotation","p_value",}
        \code{"adj_p_value","meta_p_value",} 
        \code{"adj_meta_p_value","fold_change",}
        \code{"stats","counts")}
        \item \code{exportScale=c("natural","log2","log10",}
        \code{"vst")}
        \item \code{exportValues=c("raw","normalized")} 
        \item \code{exportStats=c("mean","median","sd","mad",}
        \code{"cv","rcv")}
    } 
    \item \code{"strict_basic"}: apply a strict set of filters
    and and export statistically significant genes and basic
    annotation and statistics elements. In this case, the
    above described arguments become: 
    \itemize{
        \item \code{exonFilters=list(minActiveExons=}
        \code{list(exonsPerGene=4,minExons=2,frac=1/4))} 
        \item \code{geneFilters=list(length=list(length=750),} 
        \code{avgReads=list(averagePerBp=100,quantile=0.5),} 
        \code{expression=list(median=TRUE,mean=FALSE,}
        \code{quantile=NA,known=NA,custom=NA),} 
        \code{biotype=getDefaults("biotypeFilter",org[1]))}
        \item \code{pcut=0.01} 
        \item \code{exportWhat=c("annotation","p_value",}
        \code{"adj_p_value","meta.p.value",} 
        \code{"adj_meta_p_value","fold_change")} 
        \item \code{exportScale=c("natural","log2")} 
        \item \code{exportValues=c("normalized")} 
        \item \code{exportStats=c("mean")} 
    } 
    \item \code{"strict_normal"}: apply a strict set of filters 
    and export statistically significant genes and normal
    annotation and statistics elements. In this case, the
    above described arguments become: 
    \itemize{ 
        \item \code{exonFilters=list(minActiveExons=}
        \code{list(exonsPerGene=4,minExons=2,frac=1/4))} 
        \item \code{geneFilters=list(length=list(length=750),} 
        \code{avgReads=list(averagePerBp=100,quantile=0.5),} 
        \code{expression=list(median=TRUE,mean=FALSE,}
        \code{quantile=NA,known=NA,custom=NA),} 
        \code{biotype=getDefaults("biotypeFilter",org[1]))}
        \item \code{pcut=0.01} 
        \item \code{exportWhat=c("annotation","p_value",}
        \code{"adj_p_value","meta_p_value",} 
        \code{"adj_meta_p_value","fold_change",}
        \code{"stats","counts")}
        \item \code{exportScale=c("natural","log2")} 
        \item \code{exportValues=c("normalized")} 
        \item \code{exportStats=c("mean","sd","cv")} 
    } 
    \item \code{"strict_full"}: apply a strict set of filters 
    and export statistically significant genes and full
    annotation and statistics elements. In this case, the
    above described arguments become:
    \itemize{ 
        \item \code{exonFilters=list(minActiveExons=}
        \code{list(exonsPerGene=4,minExons=2,frac=1/4))} 
        \item \code{geneFilters=list(length=list(length=750),} 
        \code{avgReads=list(averagePerBp=100,quantile=0.5),} 
        \code{expression=list(median=TRUE,mean=FALSE,}
        \code{quantile=NA,known=NA,custom=NA),} 
        \code{biotype=getDefaults("biotypeFilter",org[1]))}
        \item \code{pcut=0.01} 
        \item \code{exportWhat=c("annotation","p_value",}
        \code{"adj_p_value","meta_p_value",} 
        \code{"adj_meta_p_value","fold_change",}
        \code{"stats","counts")}
        \item \code{exportScale=c("natural","log2","log10",}
        \code{"vst")}
        \item \code{exportValues=c("raw","normalized")} 
        \item \code{exportStats=c("mean","median","sd","mad",}
        \code{"cv","rcv")}
    } }
}
\examples{
# An example pipeline with gene counts
data("mm9GeneData",package="metaseqR2")

result <- metaseqr2(
    counts=mm9GeneCounts,
    sampleList=sampleListMm9,
    contrast=c("adult_8_weeks_vs_e14.5"),
    libsizeList=libsizeListMm9,
    annotation="embedded",
    org="mm9",
    countType="gene",
    normalization="edger",
    statistics="edger",
    pcut=0.05,
    figFormat="png",
    qcPlots="mds",
    exportWhat=c("annotation","p_value","adj_p_value","fold_change"),
    exportScale="natural",
    exportValues="normalized",
    exportStats="mean",
    exportWhere=file.path(tempdir(),"test1"),
    restrictCores=0.01,
    geneFilters=list(
        length=list(
            length=500
        ),
        avgReads=list(
            averagePerBp=100,
            quantile=0.25
        ),
        expression=list(
            median=TRUE,
            mean=FALSE,
            quantile=NA,
            known=NA,
            custom=NA
        ),
        biotype=getDefaults("biotypeFilter","mm9")
    ),
    outList=TRUE
)
head(result$data[["adult_8_weeks_vs_e14.5"]])
}
\author{
    Panagiotis Moulos
}
\name{statNbpseq}
\alias{statNbpseq}
\title{Statistical testing with NBPSeq}
\usage{
    statNbpseq(object, sampleList, contrastList = NULL,
        statArgs = NULL, libsizeList = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of NBPSeq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"nbpseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a list from NBPSeq
    normalization as the dispersions are already estimated.}

    \item{libsizeList}{an optional named list where names
    represent samples (MUST be the same as the samples
    \code{in sampleList}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    the default is the column sums of the \code{object}
    matrix.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over NBPSeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\note{
    There is currently a problem with the NBPSeq package and
    the workflow that is specific to the NBPSeq package. The
    problem has to do with function exporting as there are
    certain functions which are not recognized from the
    package internally. For this reason and until it is
    fixed, only the Smyth workflow will be available with the
    NBPSeq package.
}
\examples{
require(NBPSeq)
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeNbpseq(dataMatrix,sampleList)
p <- statNbpseq(normDataMatrix,sampleList,contrast)
}
\author{
    Panagiotis Moulos
}

\name{diagplotVolcano}
\alias{diagplotVolcano}
\title{(Interactive) volcano plots of differentially expressed 
    genes}
\usage{
    diagplotVolcano(f, p, con = NULL, fcut = 1, pcut = 0.05,
        altNames = NULL, output = "x11", path = NULL, ...)
}
\arguments{
    \item{f}{the fold changes which are to be plotted on the
    x-axis.}

    \item{p}{the p-values whose -log10 transformation is
    going to be plotted on the y-axis.}

    \item{con}{an optional string depicting a name (e.g. the
    contrast name) to appear in the title of the volcano
    diagplot.}

    \item{fcut}{a fold change cutoff so as to draw two
    vertical lines indicating the cutoff threshold for
    biological significance.}

    \item{pcut}{a p-value cutoff so as to draw a horizontal
    line indicating the cutoff threshold for statistical
    significance.}

    \item{altNames}{an optional vector of names, e.g. HUGO
    gene symbols, alternative or complementary to the unique
    names of \code{f} or \code{p} (one of them must be
    named!). It is used only in JSON output.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"}, \code{"ps"} or \code{"json"}. The latter is
    currently available for the creation of interactive
    volcano plots only when reporting the output, through the
    highcharts javascript library.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{which.plot} argument. If
    \code{output="x11"}, no output filenames are produced.
}
\description{
    This function plots a volcano plot or returns a JSON
    string which is used to render aninteractive in case of
    HTML reporting.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(5000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
M <- normalizeEdger(dataMatrix,sampleList)
p <- statEdger(M,sampleList,contrast)
ma <- apply(M[,sampleList$A],1,mean)
mb <- apply(M[,sampleList$B],1,mean)
f <- log2(ifelse(mb==0,1,mb)/ifelse(ma==0,1,ma))
diagplotVolcano(f,p[[1]],con=contrast)
#j <- diagplotVolcano(f,p[[1]],con=contrast,output="json")
}
\author{
    Panagiotis Moulos
}

\name{getDefaults}
\alias{getDefaults}
\title{Default parameters for several metaseqr functions}
\usage{
    getDefaults(what, method = NULL)
}
\arguments{
    \item{what}{a keyword determining the procedure for which
    to fetch the default settings according to method
    parameter. It can be one of \code{"normalization"},
    \code{"statistics"}, \code{"geneFilter"},
    \code{"exonFilter"} or \code{"biotypeFilter"}.}

    \item{method}{the supported algorithm included in
    metaseqR for which to fetch the default settings.
    Se Details.}
}
\value{
    A list with default setting that can be used directly in
    the call of metaseqr.
}
\description{
    This function returns a list with the default settings
    for each filtering, statistical and normalization
    algorithm included in the metaseqR package. See the
    documentation of the main function and the documentation
    of each statistical and normalization method for details.
}
\details{
    When \code{what} is \code{"normalization"}, method is 
    one of \code{"edaseq"}, \code{"deseq"}, \code{"edger"},
    \code{"noiseq"} or \code{"nbpseq"}. When \code{what} is
    \code{"statistics"}, method is one of \code{"deseq"},
    \code{"edger"}, \code{"noiseq"}, \code{"bayseq"},
    \code{"limma"} or \code{"nbpseq"}. When \code{method} is
    \code{"biotypeFilter"}, \code{what} is the input
    organism (see the main \code{\link{metaseqr2}} help page
    for a list of supported organisms).
}
\examples{
normArgsEdaseq <- getDefaults("normalization","edaseq")
statArgsEdger <- getDefaults("statistics","edger")
}
\author{
    Panagiotis Moulos
}
\name{diagplotCor}
\alias{diagplotCor}
\title{Summarized correlation plots}
\usage{
    diagplotCor(mat, type = c("heatmap", "correlogram"),
        output = "x11", path = NULL, ...)
}
\arguments{
    \item{mat}{the read counts matrix or data frame.}

    \item{type}{create heatmap of correlogram plots.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filename of the pairwise comparisons plot produced if
    it's a file.
}
\description{
    This function uses the read counts matrix to create
    heatmap or correlogram correlation plots.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
diagplotCor(dataMatrix,type="heatmap")
diagplotCor(dataMatrix,type="correlogram")
}
\author{
    Panagiotis Moulos
}

\name{statDss}
\alias{statDss}
\title{Statistical testing with DSS}
\usage{
    statDss(object, sampleList, contrastList = NULL,
        statArgs = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of DESeq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"deseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a CountDataSet from
    DESeq normalization as the dispersions are already
    estimated.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over DESeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
require(DSS)
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeDss(dataMatrix,sampleList)
p <- statDss(normDataMatrix,sampleList,contrast)
}
\author{
    Dionysios Fanidis
}

\docType{data}
\name{libsizeListMm9}
\alias{libsizeListMm9}
\title{Mouse RNA-Seq data with two conditions, four samples}
\format{a named \code{list} with library sizes.}
\source{
    ENCODE (http://genome.ucsc.edu/encode/)
}
\description{
    The library size list for \code{mm9GeneCounts}. See the
    data set description.
}
\author{
    Panagiotis Moulos
}
\keyword{datasets}
\name{getAnnotation}
\alias{getAnnotation}
\title{Annotation downloader}
\usage{
    getAnnotation(org, type, refdb = "ensembl", ver = NULL,
        rc = NULL)
}
\arguments{
    \item{org}{the organism for which to download
    annotation (one of the supported ones).}

    \item{type}{\code{"gene"}, \code{"exon"} or 
    \code{"utr"}. Same as the \code{countType} in 
    \code{\link{metaseqr2}}.}

    \item{refdb}{the online source to use to fetch 
    annotation. It can be \code{"ensembl"} (default), 
    \code{"ucsc"} or \code{"refseq"}. In the later two
    cases, an SQL connection is opened with the UCSC 
    public databases.}
    
    \item{ver}{the version of the annotation to use.}

    \item{rc}{Fraction of cores to use. Same as the 
    \code{rc} in \code{\link{buildAnnotationDatabase}}.}
}
\value{
    A data frame with the canonical (not isoforms!) genes or
    exons of the requested organism. When
    \code{type="genes"}, the data frame has the following
    columns: chromosome, start, end, gene_id, gc_content,
    strand, gene_name, biotype. When \code{type="exon"} the
    data frame has the following columns: chromosome, start,
    end, exon_id, gene_id, strand, gene_name, biotype. When 
    \code{type="utr"} the data frame has the following columns: 
    chromosome, start, end, transcript_id, gene_id, strand, 
    gene_name, biotype. The gene_id and exon_id correspond to 
    Ensembl, UCSC or RefSeq gene, transcript and exon accessions 
    respectively. The gene_name corresponds to HUGO nomenclature 
    gene names.
}
\description{
    For Ensembl based annotations, this function connects to the 
    EBI's Biomart service using the package biomaRt and downloads 
    annotation elements (gene co-ordinates, exon co-ordinates, 
    gene identifications, biotypes etc.) for each of the supported
    organisms. For UCSC/RefSeq annotations, it connects to the 
    respective SQL databases if the package \code{RMySQL} is 
    present, otherwise it downloads flat files and build a 
    temporary SQLite database to make the necessaru build 
    queries. See the help page of \code{\link{metaseqr2}}
    for a list of supported organisms.
}
\note{
    The data frame that is returned contains only "canonical"
    chromosomes for each organism. It does not contain
    haplotypes or random locations and does not contain
    chromosome M.
}
\examples{
mm10Genes <- getAnnotation("mm10","gene")
}
\author{
    Panagiotis Moulos
}
\name{normalizeEdger}
\alias{normalizeEdger}
\title{Normalization based on the edgeR package}
\usage{
    normalizeEdger(geneCounts, sampleList,
        normArgs = NULL, output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of edgeR normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"edger")} for
    an example and how you can modify it.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the edgeR native S4 object
    (DGEList). In the latter case it should be handled with
    suitable edgeR methods.}
}
\value{
    A matrix or a DGEList with normalized counts.
}
\description{
    This function is a wrapper over edgeR normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeEdger(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}
\name{createSignalTracks}
\alias{createSignalTracks}
\title{Create bigWig signal tracks}
\usage{
    createSignalTracks(targets, org, urlBase = NULL, 
        stranded = FALSE, normTo = 1e+9, exportPath = ".",
        hubInfo = list(name = "MyHub", shortLabel = "My hub",
        longLabel = "My hub", email = "someone@example.com"),
        fasta = NULL, gtf = NULL, forceHub = FALSE, 
        overwrite = FALSE, rc = NULL)
}
\arguments{
    \item{targets}{a tab-delimited file with the experimental 
    description or the output of \code{\link{readTargets}}.
    See also the \code{sampleList} argument in the main
    \code{\link{metaseqr2}} pipeline.}

    \item{org}{See the \code{org} argument in the main
    \code{\link{metaseqr2}} pipeline.}

    \item{urlBase}{a valid URL which is prepended to the created
    bigWig files.}

    \item{stranded}{Separate + and - strands and create separate
    bigWig files.}

    \item{normTo}{the total sum of signal to be used as the 
    normalization target. See also the \code{trackInfo} argument 
    in the main \code{\link{metaseqr2}} pipeline.}

    \item{exportPath}{path to export tracks.}
    
    \item{hubInfo}{information regarding the track hub created
    when \code{stranded=TRUE}. See also the \code{trackInfo} 
    argument in the main \code{\link{metaseqr2}} pipeline.}
    
    \item{overwrite}{overwrite tracks if they exist? Defaults to
    \code{FALSE}.}
    
    \item{fasta}{reference genome in FASTA format for the case
    of analyzing a custom, non-directly supported organism. It
    will be converted to the .2bit format and written along with
    a track hub.}
    
    \item{gtf}{a GTF file describing gene models in the case of
    analyzing a custom, non-directly supported organism. It will
    be converted to the .bigBed format and written along with
    a track hub.}
    
    \item{forceHub}{when \code{stranded=TRUE}, a UCSC Genome
    Browser trackhub is created, otherwise only tracklines 
    describing individual tracks. If \code{TRUE}, a trackhub is
    always created.}
    
    \item{rc}{Fraction of cores to use.}
}
\value{
    A string with the link(s) to the created tracks.
}
\description{
    This function creates bigWig files to be used for exploring
    RNA signal in genome browsers. When strands are separated,
    a UCSC genome browser trackhub is created to group tracks
    for the same sample. A link to the created data is returned.
}
\examples{
dataPath <- system.file("extdata",package="metaseqR2")
targets <- data.frame(samplename=c("C","T"),
    filename=file.path(dataPath,c("C.bam","T.bam")),  
    condition=c("Control","Treatment"),
    paired=c("single","single"),stranded=c("forward","forward"))
path <- tempdir()
write.table(targets,file=file.path(path,"targets.txt"),
    sep="\t",row.names=FALSE,quote=FALSE)
if (.Platform$OS.type == "unix")
    link <- createSignalTracks(file.path(path,"targets.txt"),"mm9")
}
\author{
    Panagiotis Moulos
}
\name{diagplotFtd}
\alias{diagplotFtd}
\title{Create False (or True) Positive (or
    Negative) curves}
\usage{
    diagplotFtd(truth, p, type = "fpc", N = 2000,
        output = "x11", path = NULL, draw = TRUE, ...)
}
\arguments{
    \item{truth}{the ground truth differential 
    expression vector. It should contain only 
    zero and non-zero elements, with zero denoting
    non-differentially expressed genes and non-zero, 
    differentially expressed genes. Such a vector 
    can be obtained for example by using the 
    \code{\link{makeSimDataSd}} function, which 
    creates simulated RNA-Seq read counts based on 
    real data. The elements of \code{truth} MUST 
    be named (e.g. each gene's name).}

    \item{p}{a p-value matrix whose rows correspond 
    to each element in the \code{truth} vector. If 
    the matrix has a \code{colnames} attribute, a 
    legend will be added to the plot using these 
    names, else a set of column names will be 
    auto-generated. \code{p} can also be a list or 
    a data frame. The p-values MUST be named (e.g. 
    each gene's name).}

    \item{type}{what to plot, can be \code{"fpc"} 
    for False Positive Curves (default), 
    \code{"tpc"} for True Positive Curves, 
    \code{"fnc"} for False Negative Curves or 
    \code{"tnc"} for True Negative Curves.}

    \item{N}{create the curves based on the 
    top (or bottom) \code{N} ranked genes 
    (default is 2000) to be used with
    \code{type="fpc"} or \code{type="tpc"}.}

    \item{output}{one or more R plotting device to 
    direct the plot result to. Supported mechanisms: 
    \code{"x11"} (default), \code{"png"}, \code{"jpg"}, 
    \code{"bmp"}, \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{draw}{boolean to determine whether to
    plot the curves or just return the calculated
    values (in cases where the user wants the
    output for later averaging for example). 
    Defaults to \code{TRUE} (make plots).}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    A named list with two members: the first member 
    (\code{ftdr}) contains the values used to create 
    the plot. The second member (\code{path}) contains 
    the path to the created figure graphic.
}
\description{
    This function creates false (or true) discovery 
    curves using a matrix of p-values (such a matrix 
    can be derived for example from the result table 
    of \code{\link{metaseqr2}} by subsetting the table 
    to get the p-values from several algorithms) 
    given a ground truth vector for differential 
    expression.
}
\examples{
p1 <- 0.001*matrix(runif(300),100,3)
p2 <- matrix(runif(300),100,3)
p <- rbind(p1,p2)
rownames(p) <- paste("gene",1:200,sep="_")
colnames(p) <- paste("method",1:3,sep="_")
truth <- c(rep(1,40),rep(-1,40),rep(0,20),
    rep(1,10),rep(2,10),rep(0,80))
names(truth) <- rownames(p)
ftdObj <- diagplotFtd(truth,p,N=100)
}
\author{
    Panagiotis Moulos
}

\name{normalizeAbsseq}
\alias{normalizeAbsseq}
\title{Normalization based on the ABSSeq package}
\usage{
    normalizeAbsseq(geneCounts, sampleList,
        normArgs = NULL, output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of DESeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"deseq")} for
    an example and how you can modify it.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the ABSSeq native S4 object
    (ABSDataSet). In the latter case it should be handled
    with suitable ABSSeq methods.}
}
\value{
    A matrix or a ABSDataSet with normalized counts.
}
\description{
    This function is a wrapper over ABSSeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeAbsseq(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Dionysios Fanidis
}
\name{statDeseq2}
\alias{statDeseq2}
\title{Statistical testing with DESeq2}
\usage{
    statDeseq2(object, sampleList, contrastList = NULL,
        statArgs = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of DESeq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"deseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a CountDataSet from
    DESeq normalization as the dispersions are already
    estimated.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over DESeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
require(DESeq2)
dataMatrix <- metaseqR2:::exampleCountData(1000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeDeseq2(dataMatrix,sampleList)
p <- statDeseq2(normDataMatrix,sampleList,contrast)
}
\author{
    Dionysios Fanidis
}

\name{getWeights}
\alias{getWeights}
\title{Get precalculated statistical test 
    weights}
\usage{
    getWeights(org = c("human", "chimpanzee", "mouse",
        "fruitfly", "arabidopsis", "rat"))
}
\arguments{
    \item{org}{\code{"human"}, \code{"chimpanzee"},
    \code{"mouse"}, \code{"fruitfly"}, \code{"arabidopsis"}
    or \code{"rat"}.}
}
\value{
    A named vector of convex weights.
}
\description{
    This function returns pre-calculated weights 
    for human, chimpanzee, mouse, fruitfly and 
    arabidopsis based on the performance of 
    simulated datasets estimated from real data 
    from the ReCount database 
    (\url{http://bowtie-bio.sourceforge.net/recount/}).
    Currently pre-calculated weights are available only 
    when all six statistical tests are used and for 
    normalization with EDASeq. For other combinations, 
    use the \code{\link{estimateAufcWeights}} function.
}
\examples{
wh <- getWeights("human")
}
\author{
    Panagiotis Moulos
}
\name{diagplotRoc}
\alias{diagplotRoc}
\title{Create basic ROC curves}
\usage{
    diagplotRoc(truth, p, sig = 0.05, x = "fpr", 
        y = "tpr", output = "x11", path = NULL,
        draw = TRUE, ...)
}
\arguments{
    \item{truth}{the ground truth differential 
    expression vector. It should contain only 
    zero and non-zero elements, with zero denoting
    non-differentially expressed genes and non-zero, 
    differentially expressed genes. Such a vector 
    can be obtained for example by using the 
    \code{\link{makeSimDataSd}} function, which 
    creates simulated RNA-Seq read counts based on 
    real data.}

    \item{p}{a p-value matrix whose rows correspond 
    to each element in the \code{truth} vector. If 
    the matrix has a \code{colnames} attribute, a 
    legend will be added to the plot using these 
    names, else a set of column names will be 
    auto-generated. \code{p} can also be a list or 
    a data frame.}

    \item{sig}{a significance level (0 < \code{sig} 
    <=1).}

    \item{x}{what to plot on x-axis, can be one of 
    \code{"fpr"}, \code{"fnr"}, \code{"tpr"}, 
    \code{"tnr"} for False Positive Rate, False
    Negative Rate, True Positive Rate and True 
    Negative Rate respectively.}

    \item{y}{what to plot on y-axis, same as 
    \code{x} above.}

    \item{output}{one or more R plotting device to 
    direct the plot result to. Supported mechanisms: 
    \code{"x11"} (default), \code{"png"}, \code{"jpg"}, 
    \code{"bmp"}, \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{draw}{boolean to determine whether to
    plot the curves or just return the calculated
    values (in cases where the user wants the
    output for later averaging for example). 
    Defaults to \code{TRUE} (make plots).}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    A named list with two members. The first member is 
    a list containing the ROC statistics: \code{TP} 
    (True Postives), \code{FP} (False Positives), 
    \code{FN} (False Negatives), \code{TN} 
    (True Negatives), \code{FPR} (False Positive Rate), 
    \code{FNR} (False Negative Rate), \code{TPR} (True 
    Positive Rate), \code{TNR} (True Negative Rate), 
    \code{AUC} (Area Under the Curve). The second is 
    the path to the created figure graphic.
}
\description{
    This function creates basic ROC curves using a 
    matrix of p-values (such a matrix can be 
    derived for example from the result table of 
    \code{\link{metaseqr2}} by subsetting the table 
    to get the p-values from several algorithms) 
    given a ground truth vector for differential 
    expression and a significance level.
}
\examples{
p1 <- 0.001*matrix(runif(300),100,3)
p2 <- matrix(runif(300),100,3)
p <- rbind(p1,p2)
rownames(p) <- paste("gene",1:200,sep="_")
colnames(p) <- paste("method",1:3,sep="_")
truth <- c(rep(1,40),rep(-1,40),rep(0,20),rep(1,10),
    rep(2,10),rep(0,80))
names(truth) <- rownames(p)
rocObj <- diagplotRoc(truth,p)
}
\author{
    Panagiotis Moulos
}

\name{statAbsseq}
\alias{statAbsseq}
\title{Statistical testing with ABSSeq}
\usage{
    statAbsseq(object, sampleList, contrastList = NULL,
        statArgs = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of DESeq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"deseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a CountDataSet from
    DESeq normalization as the dispersions are already
    estimated.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over DESeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
require(ABSSeq)
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeAbsseq(dataMatrix,sampleList)
p <- statAbsseq(normDataMatrix,sampleList,contrast)
}
\author{
    Dionysios Fanidis
}

\name{getInstalledAnnotations}
\alias{getInstalledAnnotations}
\title{Load a metaseqR2 annotation element}
\usage{
    getInstalledAnnotations(obj = NULL)
}
\arguments{
    \item{obj}{\code{NULL} or the path to a metaseqR2
    SQLite annotation database. If \code{NULL}, the
    function will try to guess the location of the 
    SQLite database.}
}
\value{
    The function returns a \code{data.frame} object with
    the installed local annotations.
}
\description{
    This function returns a data frame with information
    on locally installed, supported or custom, annotations.
    
}
\examples{
db <- file.path(system.file(package="metaseqR2"),
    "annotation.sqlite")
if (file.exists(db))
    ig <- getInstalledAnnotations(obj=db)
}
\author{
    Panagiotis Moulos
}
\name{diagplotDeHeatmap}
\alias{diagplotDeHeatmap}
\title{Diagnostic heatmap of differentially expressed genes}
\usage{
    diagplotDeHeatmap(x, scale = c("asis", "zscore"), con = NULL, 
        output = "x11", path = NULL, ...)
}
\arguments{
    \item{x}{the data matrix to create a heatmap for.}
    
    \item{scale}{value scale in the heatmap. As provided 
    (\code{scale="asis"}, default) or Z-scores
    (\code{scale="zscore"})}

    \item{con}{an optional string depicting a name (e.g. the
    contrast name) to appear in the title of the volcano
    plot.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"}, \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{whichPlot} argument. If
    \code{output="x11"}, no output filenames are produced.
}
\description{
    This function plots a heatmap of the differentially
    expressed genes produced by the metaseqr workflow, useful
    for quality control, e.g. whether samples belonging to
    the same group cluster together.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
M <- normalizeEdger(dataMatrix,sampleList)
p <- statEdger(M,sampleList,contrast)
diagplotDeHeatmap(dataMatrix[p[[1]]<0.05,])
}
\author{
    Panagiotis Moulos
}

\name{read2count}
\alias{read2count}
\title{SAM/BAM/BED file reader helper for the metaseqr2 pipeline}
\usage{
    read2count(targets, annotation, fileType = targets$type,
        transLevel = "gene", utrOpts = list(frac = 1, 
        minLength = 300, downstream = 50), interFeature = FALSE, 
        rc = NULL)
}
\arguments{
    \item{targets}{a named list, the output of 
    \code{\link{readTargets}} or an existing file with
    targets. See also the main \code{metaseqr2} man page.}

    \item{annotation}{a \code{GenomicRanges} or 
    \code{data.frame} with genomic coordinates to use for
    read counting. See also \code{\link{getAnnotation}}.}

    \item{fileType}{the type of raw input files. It can be
    \code{"bed"} for BED files or \code{"sam"}, \code{"bam"}
    for SAM/BAM files. See the same argument in the main
    \code{\link{metaseqr2}} function for the case of
    auto-guessing.}
    
    \item{transLevel}{see the \code{transLevel} argument
    in the main \code{\link{metaseqr2}} function.}
    
    \item{utrOpts}{a named list with members \code{frac} 
    which is the fraction (0-1) of the 3' UTR region to count 
    reads in, \code{minLength} the minimum acceptable 3'UTR
    length irrespective of \code{frac} and \code{downstream} 
    the number of base pairs to flank the end of the 3' UTR of 
    transcripts when analyzing Quant-Seq data.}
    
    \item{interFeature}{see the \code{inter.feature} argument
    in \code{summarizeOverlaps}.}

    \item{rc}{the fraction of the available cores to use
    in a multicore system.}
}
\value{
    A data frame with counts for each sample, ready to be
    passed to the main \code{\link{metaseqr2}} pipeline.
}
\description{
    This function is a helper for the \code{metaseqr2}
    pipeline, for reading SAM/BAM or BED files when a read
    counts file is not available. It can also be used
    very easily in an autonomous manner.
}
\examples{
dataPath <- system.file("extdata",package="metaseqR2")
targets <- data.frame(samplename=c("C","T"),
    filename=file.path(dataPath,c("C.bam","T.bam")),  
    condition=c("Control","Treatment"),
    paired=c("single","single"),stranded=c("forward","forward"))
path <- tempdir()
write.table(targets,file=file.path(path,"targets.txt"),
    sep="\t",row.names=FALSE,quote=FALSE)
geneData <- loadAnnotation("mm10","ensembl","gene")
myTargets <- readTargets(file.path(path,"targets.txt"))
if (.Platform$OS.type == "unix") {
    r2c <- read2count(targets=myTargets,
        fileType=myTargets$type,annotation=geneData)
    geneCounts <- r2c$counts
    libsizeList <- r2c$libsize
}
}
\author{
    Panagiotis Moulos
}

\name{combineSimes}
\alias{combineSimes}
\title{Combine p-values with Simes' method}
\usage{
    combineSimes(p, zerofix = NULL)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}
    
    \item{zerofix}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1.}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by 
    metaseqR using the Simes' method (see 
    reference in the main \code{\link{metaseqr2}}
    help page or in the vignette).
}
\details{
    The argument \code{zerofix} is used to correct for
    the case of a p-value which is equal to 0 as a result
    of internal numerical and approximation procedures.
    When \code{NULL}, random numbers greater than 0 and
    less than or equal to 0.5 are used to multiply the
    offending p-values with the lowest provided non-zero
    p-value, maintaining thus a virtual order of 
    significance, avoiding having the same p-values for 
    two tests and assuming that all zero p-values represent
    extreme statistical significance. When a numeric
    between 0 and 1, this number is used for the above
    multiplication instead.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineSimes(p)
}
\author{
    Panagiotis Moulos
}

\name{loadAnnotation}
\alias{loadAnnotation}
\title{Load a metaseqR2 annotation element}
\usage{
    loadAnnotation(genome, refdb, 
        level = c("gene", "transcript", "exon"),
        type = c("gene", "exon", "utr"), version="auto",
        db = file.path(system.file(package = "metaseqR2"),
            "annotation.sqlite"), summarized = FALSE, 
            asdf = FALSE, rc = NULL)
}
\arguments{
    \item{genome}{a \code{\link{metaseqr2}} supported
    organisms or a custom, imported by the user, name. See 
    also the main \code{\link{metaseqr2}} man page.}

    \item{refdb}{a \code{\link{metaseqr2}} supported
    annotation source or a custom, imported by the user, name.
    See also the main \code{\link{metaseqr2}} man page.}
    
    \item{level}{same as the \code{transLevel} in 
    \code{\link{metaseqr2}}.}
    
    \item{type}{same as the \code{countType} in 
    \code{\link{metaseqr2}}.}
    
    \item{version}{same as the \code{version} in 
    \code{\link{metaseqr2}}.}
    
    \item{db}{same as the \code{db} in 
    \code{\link{buildAnnotationDatabase}}.}
    
    \item{summarized}{if \code{TRUE}, retrieve summarized,
    non-overlaping elements where appropriate (e.g. exons).}
    
    \item{asdf}{return the result as a \code{\link{data.frame}}
    (default \code{FALSE}).}
    
    \item{rc}{same as the \code{rc} in 
    \code{\link{buildAnnotationDatabase}}.}
}
\value{
    The function returns a \code{GenomicRanges} object with
    the requested annotation.
}
\description{
    This function creates loads an annotation element from
    the local annotation database to be used with metaseqr2.
    If the annotation is not found and the organism is 
    supported, the annotation is created on the fly but not
    imported in the local database. Use
    \code{buildAnnotationDatabase} for this purpose.
    
}
\examples{
db <- file.path(system.file(package="metaseqR2"),
    "annotation.sqlite")
if (file.exists(db))
    gr <- loadAnnotation(genome="hg19",refdb="ensembl",
        level="gene",type="gene",db=db)
}
\author{
    Panagiotis Moulos
}
