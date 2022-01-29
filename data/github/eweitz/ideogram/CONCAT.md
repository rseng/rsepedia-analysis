# Ideogram API reference

Ideogram.js is a JavaScript library for chromosome visualization.  Ideogram supports drawing and animating genome-wide datasets.  Its API consists of configuration options and methods.  You can see examples of almost all these options and methods -- as well as their actual effect -- in the many [live examples](https://eweitz.github.io/ideogram/).  Simply "View source" in any example page.

# Installation
See the main [README](https://github.com/eweitz/ideogram#installation) for installation instructions.

# Configuration
To start, you need to instantiate the `Ideogram` class.  Configuration options -- which organism's genome to display, in what orientation, with which annotation data, etc. -- are passed into the `Ideogram` constructor as a JavaScript object.

For example:

```
var ideogram = new Ideogram({
  organism: 'human',
  annotations: [{
    name: 'BRCA1',
    chr: '17',
    start: 43044294,
    stop: 43125482
  }]
});
```

# Options
* [accessToken](#accesstoken)
* [ancestors](#ancestors)
* [annotations](#annotations)
* [annotationHeight](#annotationheight)
* [annotationsColor](#annotationscolor)
* [annotationsLayout](#annotationslayout)
* [annotationsPath](#annotationspath)
* [annotationTracks](#annotationtracks)
* [assembly](#assembly)
* [barWidth](#barwidth)
* [brush](#brush)
* [chrBorderColor](#chrbordercolor)
* [chrFillColor](#chrfillcolor)
* [chrHeight](#chrheight)
* [chrMargin](#chrmargin)
* [chrWidth](#chrwidth)
* [chrLabelColor](#chrlabelcolor)
* [chrLabelSize](#chrlabelsize)
* [chromosomes](#chromosomes)
* [chromosomeScale](#chromosomescale)
* [container](#container)
* [dataDir](#datadir)
* [demarcateCollinearChromosomes](#demarcatecollinearchromosomes)
* [geometry](#geometry)
* [histogramScaling](#histogramscaling)
* [heatmaps](#heatmaps)
* [filterable](#filterable)
* [fontFamily](#fontfamily)
* [fullChromosomeLabels](#fullchromosomelabels)
* [legend](#legend)
* [onBrushMove](#onbrushmove)
* [onDidRotate](#ondidrotate)
* [onDrawAnnots](#ondrawannots)
* [onLoadAnnots](#onloadannots)
* [onLoad](#onload)
* [onWillShowAnnotTooltip](#onwillshowannottooltip)
* [organism](#organism)
* [orientation](#orientation)
* [perspective](#perspective)
* [ploidy](#ploidy)
* [ploidyDesc](#ploidydesc)
* [rangeSet](#rangeset)
* [resolution](#resolution)
* [rotatable](#rotatable)
* [rows](#rows)
* [sex](#sex)
* [showBandLabels](#showbandlabels)
* [showChromosomeLabels](#showchromosomelabels)
* [showAnnotTooltip](#showannottooltip)
* [showFullyBanded](#showfullybanded)
* [showNonNuclearChromosomes](#shownonnuclearchromosomes)

## accessToken
String.  Optional.  OAuth 2.0 access token.  Enables authentication and authorization. This can be useful for controlling access to private annotation data.

## ancestors
Object.  Optional.  A map associating ancestor labels to colors.  Used to color chromosomes from different ancestors in polyploid genomes.  Example in [Ploidy, recombination](https://eweitz.github.io/ideogram/ploidy-recombination).

## annotations
Array.  Optional.  A list of annotation objects.  Each annotation object has at least a chromosome name (chr), start coordinate (start), and stop coordinate (stop).  Annotation objects can also have a name, color, shape, and track index.  Example in [Annotations, basic](https://eweitz.github.io/ideogram/annotations-basic).

See also [annotationsPath](#annotationspath).

## annotationHeight
Number.  Optional.  The height of each annotation. Example in [Annotations, tracks](https://eweitz.github.io/ideogram/annotations-tracks).

## annotationsColor
String.  Optional.  Default: "#F00" (i.e., red).  The color of each annotation.  Example in [Multiple, primates](https://eweitz.github.io/ideogram/multiple-primates).

## annotationsLayout
String.  Optional.  Default: "tracks".

The layout of this ideogram's annotations.  One of "tracks", "heatmap", "histogram", or "overlay".

### annotationsLayout: 'tracks'
Lay out annotations in tracks beside each chromosome.  There can be more than one track, which is useful for displaying annotations by category (e.g. pathogenic, unknown significance, benign).  Example in [Annotations, tracks](https://eweitz.github.io/ideogram/annotations-tracks).

### annotationsLayout: 'heatmap'
Lay out annotations in heatmap beside each chromosome.  Plot individual annotations like `annotationsLayout: 'tracks'`, with the scalability of `annotationsLayout: 'histogram'`.  Each chromosome can have one or more heatmap tracks.  Use with the [heatmaps](#heatmaps) configuration option.  Example in [Annotations, heatmap](https://eweitz.github.io/ideogram/annotations-heatmap).

### annotationsLayout: 'heatmap-2d'
Lay out annotations in a two-dimensional zheatmap beside a single chromosome.  Enables visualizing raw data summarized in
`annotationsLayout: 'heatmap'`.  Example in [Annotations, 2D heatmap](https://eweitz.github.io/ideogram/annotations-heatmap-2d).

### annotationsLayout: 'histogram'
Lay out annotations in a histogram.  This clusters annoatations by location, such that each cluster or bin is shown as a bar.  The height of the bar represent the number of annotations in that genomic range.  This option is useful for summarizing the distribution of many (1000+) features througout the genome.  Example in [Annotations, histogram](https://eweitz.github.io/ideogram/annotations-histogram).

### annotationsLayout: 'overlay'
Lay out annotations directly over chromosomes.  This is the most space-efficient annotation layout option.  Example in [Annotations, overlaid](https://eweitz.github.io/ideogram/annotations-overlaid).

## annotationsPath
String.  Optional.  An absolute or relative URL to a JSON file containing annotation objects.  Example in [Annotations, overlaid](https://eweitz.github.io/ideogram/annotations-overlaid).

See also [annotations](#annotations).

## annotationTracks
Array.  Optional.  A list of objects with metadata for each track, e.g. DOM `id`, display name, color, shape.  Example in [Annotations, tracks](https://eweitz.github.io/ideogram/annotations-tracks).

## assembly
String.  Optional.  Default: latest RefSeq assembly for specified organism.  The genome assembly to display.  Takes assembly name (e.g. "GRCh37"), RefSeq accession (e.g. "GCF_000306695.2"), or GenBank accession (e.g. "GCA_000005005.5").  Example in [Annotations, histogram](https://eweitz.github.io/ideogram/annotations-histogram).

## barWidth
Number.  Optional.  Default: 3.  The pixel width of bars drawn when `annotationsLayout: 'histogram'`.  Example in [Annotations, histogram](https://eweitz.github.io/ideogram/annotations-histogram).

## brush
String.  Optional.  Default: null.  Genomic coordinate range (e.g. "chr1:104325484-119977655") for a [brush](https://github.com/d3/d3-brush) on a chromosome.  Useful when ideogram consists of one chromosome and you want to be able to focus on a region within that chromosome, and create an interactive sliding window to other regions.  Example in [Brush](https://eweitz.github.io/ideogram/brush).

## chrBorderColor
String.  Optional.  Default: "#000".  The color of the border for each chromosome.

## chrFillColor
String or Object.  Optional.  Default: "#AAA".  The color with which the chromosome is filled; its interior color.  Customizes chromosome arm color if specified as string, e.g. `chrFillColor: 'green'`.  Customizes chromosome arm and centromere colors if specified as object, e.g. `chrFillColor: {arm: '#AEA', centromere: '#EEA '}`.  Example in [Color, chromosomes](https://eweitz.github.io/ideogram/color-chromosomes).

## chrHeight
Number.  Optional.  Default: 400.  The pixel height of the tallest chromosome in the ideogram.  Examples in [Layout, small](https://eweitz.github.io/ideogram/layout-small) and [Annotations, basic](https://eweitz.github.io/ideogram/annotations-basic).

## chrMargin
Number.  Optional.  Default: 10.  The pixel space of the margin between each chromosome.  Example in [Multiple, primates](https://eweitz.github.io/ideogram/multiple-primates).

## chrWidth
Number.  Optional.  Default: 10.  The pixel width of each chromosome.  Example in [Annotations, tracks](https://eweitz.github.io/ideogram/annotations-tracks).

## chrLabelColor
String. Optional.  Default: "#000".  The color of the label for each chromosome.

## chrLabelSize
Number.  Optional.  Default: 9.  The pixel font size of chromosome labels.  Example in [Differential expression of genes](https://eweitz.github.io/ideogram/differential-expression).

## chromosomes
Array.  Optional.  Default: all chromosomes in assembly.  A list of the names of chromosomes to display.  Useful for depicting a subset of the chromosomes in the genome, e.g. a single chromosome.  Example in [Annotations, basic](https://eweitz.github.io/ideogram/annotations-basic).

## chromosomeScale
String.  Optional.  Default: absolute.  Either "absolute" or "relative".  Used when comparing multiple genomes.  If absolute, chromosomes will be scaled by base pairs in each genome.  If relative, the first chromosme in each genome will be of equal length, and subsequent chromosomes will be scaled relative to the first chromosome.

## container
String.  Optional.  Default: "body".  CSS selector of the HTML element that will contain the ideogram.  Example in [Layout, small](https://eweitz.github.io/ideogram/layout-small).

## dataDir
String.  Optional.  Default: "../data/bands/native/".  Absolute or relative URL of the directory containing data needed to draw banded chromosomes.  Example in [GeneExpressionAging/ideogram](https://ncbi-hackathons.github.io/GeneExpressionAging/ideogram).

## demarcateCollinearChromosomes
Boolean.  Optional.  Default: true.  Whether to demarcate colllinear chromosomes.  Puts a dark border around the perimeter of each track-chromosomes block in track sets for chromosomes arranged in collinear geometry.  Example in [Geometry, collinear](https://eweitz.github.io/ideogram/geometry-collinear).

## geometry
String.  Optional.  Use "geometry: collinear" to arrange all chromosomes in one line, unlike the usual parallel view.  Example in [Geometry, collinear](https://eweitz.github.io/ideogram/geometry-collinear).

## histogramScaling
String.  Optional.  Default: "absolute".  One of "absolute" or "relative".  The technique to use in scaling the height of histogram bars.  The "absolute" value sets bar height relative to tallest bar in _all_ chromosomes, while "relative" sets bar height relative to tallest bar in _each_ chromosome.

## heatmaps
Array.  Optional.  Array of heatmap objects.  Each heatmap object has a `key` string and a `thresholds` array.  The `key` property specifies the annotations key value to depict in the heatmap.  The `thresholds` property specifies a list of two-element "threshold" lists, where the first element is the threshold value and the second is the threshold color.  The threshold values are a list of ranges to use in coloring
the heatmap.  Threshold values are specified in ascending order.  Example in [Annotations, heatmap](https://eweitz.github.io/ideogram/annotations-heatmap).

## filterable
Boolean.  Optional.  Whether annotations should be filterable.  Example in [Annotations, histogram](https://eweitz.github.io/ideogram/annotations-histogram).

## fontFamily
String.  Optional.  The font family to use for text in the ideogram, e.g. `fontFamily: "'Montserrat', sans-serif"`.

## fullChromosomeLabels
Boolean.  Optional.  Whether to include abbreviation species name in chromosome label.  Example in [Homology, interspecies](https://eweitz.github.io/ideogram/homology-interspecies).

## legend
Array.  Optional.  List of objects describing annotations.  Example in [Annotations, tracks](https://eweitz.github.io/ideogram/annotations-tracks).

## onBrushMove
Function.  Optional.  Callback function to invoke when brush moves.  Example in [Brush](https://eweitz.github.io/ideogram/brush).

## onDidRotate
Function.  Optional.  Callback function to invoke after chromosome has rotated.

## onDrawAnnots
Function.  Optional.  Callback function to invoke when annotations are drawn.  This is useful for when loading and drawing large annotation datsets.  Example in [web-tests.js](https://github.com/eweitz/ideogram/blob/b701dc0b76089842d50860c8c6cf5aa6d8dec564/test/web-test.js#L395).

## onLoad
Function.  Optional.  Callback function to invoke when chromosomes are loaded, i.e. rendered on the page.  Example in [Annotations, external data](https://eweitz.github.io/ideogram/annotations-external-data).

## onLoadAnnots
Function.  Optional.  Callback function to invoke when annotations are downloaded and ready for data transformation.

## onWillShowAnnotTooltip
Function.  Optional.  Callback function to invoke immediately before annotation tooltip is shown.  The tooltip shows the genomic range and, if available, name of the annotation.  This option can be useful to e.g. enhance the displayed annotation name, say by transforming a gene name into a hyperlink to a gene record web page.  Example in [Annotations, external data](https://eweitz.github.io/ideogram/annotations-external-data).

## organism
String or number or array.  Required.  Organism(s) to show chromosomes for.  Supply organism's name as a string (e.g. `"human"`) or organism's NCBI Taxonomy ID (taxid, e.g. `9606`) to display chromosomes from a single organism, or an array of organisms' names or taxids to display chromosomes from multiple species or other taxa.  Example in [Human]( https://eweitz.github.io/ideogram/human).

## orientation
String.  Optional.  Default: vertical.  The orientation of chromosomes on the page.  Example in [Mouse]( https://eweitz.github.io/ideogram/mouse).

## perspective
String.  Optional.  Use `perspective: 'comparative'` to enable annotations between two chromosomes, either within the same organism or different organisms.  Examples in [Homology, basic](https://eweitz.github.io/ideogram/homology-basic) and [Homology, interspecies](https://eweitz.github.io/ideogram/homology-interspecies).

## ploidy
Number.  Optional.  Default: 1.  The ploidy, i.e. number of chromosomes to depict for each chromosome set.  Useful for more biologically accurate rendering of genomes that are diploid, triploid, etc.  Example in [Ploidy, basic](https://eweitz.github.io/ideogram/ploidy-basic).

## ploidyDesc
Array.  Optional.  Description of ploidy in each chromosome set in terms of ancestry composition.  Example in [Ploidy, recombination](https://eweitz.github.io/ideogram/ploidy-recombination).

## rangeSet
Array.  Optional.  List of objects describing segments of recombination among chromosomes in a chromosome set.  Example in [Ploidy, recombination](https://eweitz.github.io/ideogram/ploidy-recombination).

## resolution
Number.  Optional.  Default: highest resolution available for specified genome assembly.  The resolution of cytogenetic bands to show for each chromosome.  The quantity refers to approximate value in bands per haploid set (bphs).  One of 400, 550, or 850.  Example in [Layout, small](https://eweitz.github.io/ideogram/layout-small).

See also: [showFullyBanded](#showfullybanded).

## rotatable
Boolean.  Optional.  Default: true.  Whether chromosomes are rotatable upon clicking them.  Example in [Layout, small](https://eweitz.github.io/ideogram/layout-small).

## rows
Number.  Optional.  Default: 1.  Number of rows to arrange chromosomes into.  Useful for putting ideogram into a small container, or when dealing with genomes that have many chromosomes.  Example in [Layout, small](https://eweitz.github.io/ideogram/layout-small).

## sex
String.  Optional.  Default: male.  The biological sex of the organism.  Useful for omitting chromosome Y in female mammals.  Currently only supported for organisms that use XY sex-determination.  Examples in [Ploidy, basic](https://eweitz.github.io/ideogram/ploidy-basic).

## showBandLabels
Boolean.  Optional.  Default: false.  Whether to show cytogenetic band labels, e.g. 1q21.  Example in [Annotations, basic](https://eweitz.github.io/ideogram/annotations-basic).

## showChromosomeLabels
Boolean.  Optional.  Defaut: true.  Whether to show chromosome labels, e.g. 1, 2, 3, X, Y.  Example in [Annotations, basic](https://eweitz.github.io/ideogram/annotations-basic).

## showAnnotTooltip
Boolean.  Optional.  Default: true.  Whether to show a tooltip upon mousing over an annotation.  Example in [Multiple, trio SV](https://eweitz.github.io/ideogram/multiple-trio-sv).

## showFullyBanded
Boolean.  Optional.  Default: true.  Whether to show fully banded chromosomes for genomes that have sufficient data.  Useful for showing simpler chromosomes of cytogenetically well-characterized organisms, e.g. human, beside chromosomes of less studied organisms, e.g. chimpanzee.  Example in [Multiple, primates](https://eweitz.github.io/ideogram/multiple-primates).

See also: [resolution](#resolution).

## showNonNuclearChromosomes
Boolean.  Optional.  Default: false.  Whether to show non-nuclear chromosomes, e.g. for mitochondrial (MT) and chloroplast (CP) DNA.  Example in [Eukaryotes: Sus scrofa](https://eweitz.github.io/ideogram/eukaryotes?org=sus-scrofa).
# ideogram

[![Build Status](https://travis-ci.org/eweitz/ideogram.svg?branch=master)](https://travis-ci.org/eweitz/ideogram)
[![Coverage Status](https://coveralls.io/repos/github/eweitz/ideogram/badge.svg)](https://coveralls.io/github/eweitz/ideogram)

[Ideogram.js](https://eweitz.github.io/ideogram/) is a JavaScript library for [chromosome visualization](https://speakerdeck.com/eweitz/designing-genome-visualizations-with-ideogramjs). 

Ideogram supports drawing and animating genome-wide datasets for [human](https://eweitz.github.io/ideogram/human), [mouse](https://eweitz.github.io/ideogram/mouse), and [many other eukaryotes](https://eweitz.github.io/ideogram/eukaryotes).  The [Ideogram API](https://github.com/eweitz/ideogram/blob/master/api.md) for annotations supports [histograms](https://eweitz.github.io/ideogram/annotations-histogram), [heatmaps](https://eweitz.github.io/ideogram/annotations-heatmap), [overlays](https://eweitz.github.io/ideogram/annotations-overlaid), and points of arbitrary shape and color layered in [tracks](https://eweitz.github.io/ideogram/annotations-tracks). Ideogram can depict haploid, [diploid](https://eweitz.github.io/ideogram/ploidy-basic) or higher ploidy genomes (e.g. plants), as well as aneuploidy, [genetic recombination](https://eweitz.github.io/ideogram/ploidy-recombination), and [homologous features](https://eweitz.github.io/ideogram/homology-basic) between chromosomes. 

Ideogram can be embedded as a [reusable component](https://github.com/eweitz/ideogram#usage) in any web page or application, and leverages D3.js and SVG to achieve fast, crisp client-side rendering. You can also integrate Ideogram with JavaScript frameworks like [Angular](https://github.com/eweitz/ideogram/tree/master/examples/angular), [React](https://github.com/eweitz/ideogram/tree/master/examples/react), and [Vue](https://github.com/eweitz/ideogram/tree/master/examples/vue), as well as data science platforms like [R](https://github.com/eweitz/ideogram/tree/master/examples/r) and [Jupyter Notebook](https://github.com/eweitz/ideogram/tree/master/examples/jupyter).

[![All human genes](https://raw.githubusercontent.com/eweitz/ideogram/master/examples/vanilla/ideogram_histogram_all_human_genes.png)](https://eweitz.github.io/ideogram/annotations_histogram.html)

Check out [live examples](https://eweitz.github.io/ideogram/), get [up and running](#installation) with your own deployment, skim [basic usage](#usage), or dive into the [full API](api.md)!

A [broader overview](https://speakerdeck.com/eweitz/ideogramjs-chromosome-visualization-with-javascript) is also available.  It discusses Ideogram's chromosome biology, technical architecture, and more.

# Installation

To link directly to the latest release, copy this snippet:
```
<script src="https://cdn.jsdelivr.net/npm/ideogram@1.33.0/dist/js/ideogram.min.js"></script>
```

You can also easily use the library locally:
```
$ cd <your local web server document root>
$ git clone https://github.com/eweitz/ideogram.git
```

Then go to [http://localhost/ideogram/examples/](http://localhost/ideogram/examples/).

Or, if you use npm:
```
npm install ideogram
```

You can then [import](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/import) Ideogram into an application like so:
```
import Ideogram from 'ideogram';
```


# Usage
```html
<head>
  <script src="https://cdn.jsdelivr.net/npm/ideogram@1.33.0/dist/js/ideogram.min.js"></script>
</head>
<body>
<script>
  var ideogram = new Ideogram({
    organism: 'human',
    annotations: [{
      name: 'BRCA1',
      chr: '17',
      start: 43044294,
      stop: 43125482
    }]
  });
</script>
</body>
```

Many more usage examples are available at https://eweitz.github.io/ideogram/.

You can also find examples of integrating Ideogram with JavaScript frameworks like [Angular](https://github.com/eweitz/ideogram/tree/master/examples/angular), [React](https://github.com/eweitz/ideogram/tree/master/examples/react), and [Vue](https://github.com/eweitz/ideogram/tree/master/examples/vue), as well as data science platforms like [R](https://github.com/eweitz/ideogram/tree/master/examples/r) and [Jupyter Notebook](https://github.com/eweitz/ideogram/tree/master/examples/jupyter). 


# API

See the [Ideogram API reference](api.md) for detailed documentation on configuration options and methods.
# Ideogram in Vue
This is a very basic example of integrating Ideogram with [Vue](https://vuejs.org/).

More examples showing the breadth of Ideogram's functionality are at https://eweitz.github.io/ideogram/.

# Install

``` bash
git clone https://github.com/eweitz/ideogram
cd ideogram/examples/vue
npm install
npm run dev
```

Then open your browser to http://localhost:8080.

# Output

After executing the steps above, you should see the following: 
![Ideogram in Vue screenshot](https://raw.githubusercontent.com/eweitz/ideogram/master/examples/vue/ideogram_vue_example.png)
# Ideogram in Angular

This is a very basic example of integrating Ideogram with [Angular](https://angular.io/).  

More examples showing the breadth of Ideogram's functionality are at https://eweitz.github.io/ideogram/.

# Install
```
git clone origin https://github.com/eweitz/ideogram
cd ideogram/examples/angular
npm install
ng serve
```

Then open your browser to http://localhost:4200.

# Output
After executing the steps above, you should see the following: 
![Ideogram in Angular screenshot](https://raw.githubusercontent.com/eweitz/ideogram/master/examples/angular/ideogram_angular_example.png)
# Ideogram in Jupyter

This is a very basic example of integrating Ideogram with [Jupyter](https://jupyter.org/).

More examples showing the breadth of Ideogram's functionality are at https://eweitz.github.io/ideogram/.

# Install
```
git clone origin https://github.com/eweitz/ideogram
cd ideogram/examples/jupyter
jupyter notebook
```

Then, in the browser window that opens, click on "ideogram.pynb", then click "Run" twice.

# Output
After executing the steps above, you should see the following:
![Ideogram in Jupyter screenshot](https://raw.githubusercontent.com/eweitz/ideogram/master/examples/jupyter/ideogram_jupyter_example.png)
# Ideogram in R

You can find an example of integrating Ideogram with R in [wangtulao/ideogRam](https://github.com/wangtulao/ideogRam), by Freeman S. Wang.

More examples showing the breadth of Ideogram's functionality are at https://eweitz.github.io/ideogram/.# Ideogram in React
This is a very basic example of integrating Ideogram with [React](https://reactjs.org/).

For more examples, check out [React Ideogram](https://github.com/eweitz/react-ideogram)!

# Install
```
git clone https://github.com/eweitz/ideogram
cd ideogram/examples/react
npm install
npm start
```

# Output
After executing the steps above, you should see the following:
![Ideogram in React screenshot](https://raw.githubusercontent.com/eweitz/ideogram/master/examples/react/ideogram_react_example.png)
# Annotations

Ideogram.js enables developers to load sets of annotations formatted in JSON.  This can be used to retrieve annotation data from a server, if embedding Ideogram in a web application; or to retrieve annotations from a local file system, if using Ideogram through its command-line interface.  

## Efficiency
When combined with standard `gzip` compression, the compact design of the JSON enables efficient network transmission of annotation data.  A compressed Ideogram annotation set of 21,832 human genes, including data on expression level and gene type, is 337 KB in size and takes less than 700 ms to download on a regular 4G connection (4 Mb/s download bandwidth, 20 ms latency) as measured using Chrome Developer Tools.

## Format
Each annotation in an Ideogram annotation set is represented by an array.  The meaning of elements in each of those arrays is indicated by a single array of keys at the top level of the annotation set JSON object.  To further reduce file size, annotations are grouped by chromosome, eliminating the need to specify chromosome name in each annotation.

Annotations must have at least three elements: name, start position and length.  Annotations may also have elements with integer values defining either A) track index or B) filter index.  Track index is the ordinal position of the track in which the annotation will appear.  Filter index maps to a value in an array of objects defining more filter information, including label (e.g. "Likely pathogenic") and a string identifier (e.g. `likely-pathogenic`) than can be used as a part of a DOM selector or URL parameter specifying the filter.

## Example
```
{ 
  "keys": ["name", "start", "length", "expression-level", "gene-type"]
  "annots": [{
    "chr": "1",
    "annots": [
      ["MTOR", 11106535, 155972, 5, 5],
      ["F5", 169514166, 72422, 4, 2],
  },
  {
    "chr": "2",
    "annots": [
      ["APOB", 21001429, 42644, 2, 5],
      ["CASP8", 201233443, 54268, 6, 5]
  }
  ]
}
```
# Native annotation formats

Ideogram.js supports specifying annotations in native formats, in addition to standard bioinformatics file formats like BED.

These native Ideogram annotation formats have two dimensions: [density and origin](#density-and-origin).  Particular formats excel [in different scenarios](#use-case-matrix).  The default annotation format is designed for ease of use, but [sophisticated extensions](#advanced-formats) are possible.

Note that these native formats all have a similar schema.  They are distinguished implicitly by their data layout, not explicitly as ideogram configuration properties.

## Density and origin
The density of annotation properties can be dense, where each genomic feature (e.g. each gene) has annotations on all tracks, or sparse, where each feature has an annotation on only one track.

The origin of optional annotation properties can be the client, as a parameter used in ideogram configuration, or the server, as a key in the annotation data file.

## Use case matrix
These two dimensions of annotation formats – density and origin – each have two values.  Density can be dense or sparse.  Origin can be client or server.  

Thus, four broad annotation formats are supported, each optimizing for specific trade-offs of user and developer experience and capacity, as well as different biological use cases.

| Format | Control | Configuration | Transfer | Information |  Usabliity |  Developer capacity | Example scenario |
|---|---|---|---|---|---|---|---|
| Sparse client | Dynamic | Complex | Fast | Low | Interactive and glanceable | Front-end | Clinical variation |
| Sparse server | Static | Simple | Fast | Low | Easy and glanceable |  Back-end | Clinical variation  |
| Dense client | Dynamic |  Complex | Slow | High | Interactive and rich | Front-end | Gene expression research |
| Dense server | Static | Simple | Slow | High | Easy and rich | Back-end | Gene expression research |

## Advanced formats
Ideogram provides a simple, easy, and fast annotation interface by default.  The defaults can be adjusted with minor developer effort.  Advanced implementations of annotation formatting are also possible.  

For example, server annotation property keys can be overridden by client keys, providing easy defaults that can be dynamically configured.  In fact, in the absence of server keys, Ideogram will fall back on built-in default values for all optional keys.  

Sparse annotations can also be augmented on the client to show dense, rich information for each genomic feature.  This can enable inflating compressed annotation datasets via client-side functions, or allow users of such ideogram embeds to easily move from one annotation display layout to another.

Such implementations increase application size and complexity, but provide tailored performance and enhance usage flexibility.
