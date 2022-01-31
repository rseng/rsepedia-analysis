treebase
========

[![Build Status](https://api.travis-ci.org/ropensci/treeBASE.png)](https://travis-ci.org/ropensci/treeBASE)
[![Build status](https://ci.appveyor.com/api/projects/status/d74vw3cpwh0kdg3e/branch/master)](https://ci.appveyor.com/project/sckott/treebase/branch/master)

_An R package for discovery, access and manipulation of online phylogenies_

- [Publication in Methods in Ecology and Evolution](http://dx.doi.org/10.1111/j.2041-210X.2012.00247.x)
- [Development version source code on github](https://github.com/ropensci/treebase)
- [HTML package documentation](http://ropensci.github.com/treeBASE/)
- [Report issues, bugs or feature requests](https://github.com/ropensci/treebase/issues)


Installation
------------

`treebase` is available from CRAN.  You can install the latest version from the development website on github using the `devtools` package from within R.  Make sure you have the latest version for the best experience.

```r
library(devtools)
install_github("ropensci/treebase")
```

Getting Started
---------------

Use of the `treebase` package should be relatively straight forward:

```r
library(treebase)
Phylogenies_from_Huelsenbeck <- search_treebase("Huelsenbeck", "author")
```

More interesting examples will take advantage of `R` to loop over large amounts of treebase data that would be to tiresome to search for, download and analyze by hand. Welcome to the era of big data phylogenetics.  

- Browse the examples in the [documentation](http://ropensci.github.com/treeBASE/)
- We are preparing a short manuscript to introduce the motivation, functions, and use-cases for the `treebase` package.  Meanwhile, a [preprint is available](https://github.com/ropensci/treeBASE/blob/master/inst/doc/treebase/treebase_github.md) as a dynamic document, where all of the examples shown are produced by the code shown using `knitr`. [See source code](https://github.com/ropensci/treeBASE/blob/master/inst/doc/treebase/treebase.Rmd).  


- treebase is part of the [rOpenSci Project](http://ropensci.github.com)

References
----------

* Carl Boettiger, Duncan Temple Lang (2012). Treebase: An R package for discovery, access and manipulation of online phylogenies, Methods in Ecology and Evolution. doi:10.1111/j.2041-210X.2012.00247.x
Dear CRAN Maintainers,

This commit drops the suggested dependency *laser*, as requested, since that package is now to be archived. (Also now uses markdown format to enclose URL as requested)

Sincerely,

Carl Boettiger





Appendix
========

Reproducible computation: A diversification rate analysis
---------------------------------------------------------

This appendix illustrates the diversification rate analysis discussed in
the main text.  For completeness we begin by executing the code discussed
in the manuscript which locates, downloads, and imports the relevant data:





Different diversification models make different assumptions about the rate
of speciation, extinction, and how these rates may be changing over time.
The original authors consider eight different models, implemented in the
laser package [@Rabosky2006b]. This code fits each of the eight models
to that data:



```r
library(ape)
library(laser)
models <- list(
  yule = pureBirth(bt),  
  birth_death = bd(bt),     
  yule.2.rate = yule2rate(bt),
  linear.diversity.dependent = DDL(bt),    
  exponential.diversity.dependent = DDX(bt),
  varying.speciation_rate = fitSPVAR(bt),  
  varying.extinction_rate = fitEXVAR(bt),  
  varying_both = fitBOTHVAR(bt))
```




Each of the model estimate includes an Akaike Information Criterion
(AIC) score indicating the goodness of fit, penalized by model complexity
(lower scores indicate better fits) We ask R to tell us which model has
the lowest AIC score,



```r
aics <- sapply(models, function(model) model$aic)
best_fit <- names(models[which.min(aics)])
```




and confirm the result presented in @Derryberry2011; that the best-fit
model in the laser analysis was a Yule (net diversification rate) model
with two separate rates.  

We can ask ` TreePar ` to see if a model with more rate shifts is favoured
over this single shift, a question that was not possible to address using
the tools provided in `laser`. The previous analysis also considers a
birth-death model that allowed speciation and extinction rates to be
estimated separately, but did not allow for a shift in the rate of such
a model.  In the main text we introduced a model from @Stadler2011 that
permitted up to 3 change-points in the speciation rate of the Yule model,







```r
yule_models <- bd.shifts.optim(x, sampling = c(1,1,1,1), 
  grid = 5, start = 0, end = 60, yule = TRUE)[[2]]
```




We can also compare the performance of models which allow up to three
shifts while estimating extinction and speciation rates separately:






```r
birth_death_models <- bd.shifts.optim(x, sampling = c(1,1,1,1), 
  grid = 5, start = 0, end = 60, yule = FALSE)[[2]]
```




The models output by these functions are ordered by increasing number
of shifts.  We can select the best-fitting model by AIC score, which is
slightly cumbersome in `TreePar` syntax.  First, we compute the AIC scores
of both the `yule_models` and the `birth_death_models` we fitted above,



```r
yule_aic <- 
sapply(yule_models, function(pars)
                    2 * (length(pars) - 1) + 2 * pars[1] )
birth_death_aic <- 
sapply(birth_death_models, function(pars)
                            2 * (length(pars) - 1) + 2 * pars[1] )
```




Then we generate a list identifying which model has the best (lowest)
AIC score among the Yule models and which has the best AIC score among
the birth-death models,



```r
best_no_of_rates <- list(Yule = which.min(yule_aic), 
                         birth.death = which.min(birth_death_aic))
```




The best model is then whichever of these has the smaller AIC value.  



```r
best_model <- which.min(c(min(yule_aic), min(birth_death_aic)))
```





which still confirms that the Yule 2-rate  model is still the best choice based on AIC score.  Of the eight models 
in this second analysis, only three were in the original set considered 
(Yule 1-rate and 2-rate, and birth-death without a shift), so we could by
no means have been sure ahead of time that a birth death with a shift, or
a Yule model with a greater number of shifts, would not have fitted better.  


# References

<h1 id="appendix">Appendix</h1>
<h2 id="reproducible-computation-a-diversification-rate-analysis">Reproducible computation: A diversification rate analysis</h2>
<p>Different diversification models make different assumptions about the rate of speciation, extinction, and how these rates may be changing over time. The authors consider eight different models, implemented in the laser package <span class="citation">(Rabosky 2006)</span>. This code fits each of the eight models to that data:</p>
<pre class="sourceCode r"><code class="sourceCode r"><span class="kw">library</span>(ape)
<span class="kw">library</span>(laser)
bt &lt;- <span class="kw">branching.times</span>(derryberry)
models &lt;- <span class="kw">list</span>(
  <span class="dt">yule =</span> <span class="kw">pureBirth</span>(bt),  
  <span class="dt">birth_death =</span> <span class="kw">bd</span>(bt),     
  <span class="dt">yule.2.rate =</span> <span class="kw">yule2rate</span>(bt),
  <span class="dt">linear.diversity.dependent =</span> <span class="kw">DDL</span>(bt),    
  <span class="dt">exponential.diversity.dependent =</span> <span class="kw">DDX</span>(bt),
  <span class="dt">varying.speciation_rate =</span> <span class="kw">fitSPVAR</span>(bt),  
  <span class="dt">varying.extinction_rate =</span> <span class="kw">fitEXVAR</span>(bt),  
  <span class="dt">varying_both =</span> <span class="kw">fitBOTHVAR</span>(bt))</code></pre>
<p>Each of the model estimate includes an AIC score indicating the goodness of fit, penalized by model complexity (lower scores indicate better fits) We ask R to tell us which model has the lowest AIC score,</p>
<pre class="sourceCode r"><code class="sourceCode r">aics &lt;- <span class="kw">sapply</span>(models, function(model) model$aic)
best_fit &lt;- <span class="kw">names</span>(models[<span class="kw">which.min</span>(aics)])</code></pre>
<p>and confirm the result presented in E. P. Derryberry et al. <span class="citation">(2011)</span>; that the yule.2.rate model is the best fit to the data.</p>
<p>The best-fit model in the laser analysis was a Yule (net diversification rate) model with two separate rates. We can ask <code>TreePar</code> to see if a model with more rate shifts is favoured over this single shift, a question that was not possible to address using the tools provided in <code>laser</code>. The previous analysis also considers a birth-death model that allowed speciation and extinction rates to be estimated separately, but did not allow for a shift in the rate of such a model. In the main text we introduced a model from Stadler <span class="citation">(2011)</span> that permitted up to 3 change-points in the speciation rate of the Yule model,</p>
<pre class="sourceCode r"><code class="sourceCode r">yule_models &lt;- <span class="kw">bd.shifts.optim</span>(x, <span class="dt">sampling =</span> <span class="kw">c</span>(<span class="dv">1</span>,<span class="dv">1</span>,<span class="dv">1</span>,<span class="dv">1</span>), 
  <span class="dt">grid =</span> <span class="dv">5</span>, <span class="dt">start =</span> <span class="dv">0</span>, <span class="dt">end =</span> <span class="dv">60</span>, <span class="dt">yule =</span> <span class="ot">TRUE</span>)[[<span class="dv">2</span>]]</code></pre>
<p>We can also compare the performance of models which allow up to three shifts while estimating extinction and speciation rates separately:</p>
<pre class="sourceCode r"><code class="sourceCode r">birth_death_models &lt;- <span class="kw">bd.shifts.optim</span>(x, <span class="dt">sampling =</span> <span class="kw">c</span>(<span class="dv">1</span>,<span class="dv">1</span>,<span class="dv">1</span>,<span class="dv">1</span>), 
  <span class="dt">grid =</span> <span class="dv">5</span>, <span class="dt">start =</span> <span class="dv">0</span>, <span class="dt">end =</span> <span class="dv">60</span>, <span class="dt">yule =</span> <span class="ot">FALSE</span>)[[<span class="dv">2</span>]]</code></pre>
<p>The models output by these functions are ordered by increasing number of shifts.<br />We can select the best-fitting model by AIC score, which is slightly cumbersome in <code>TreePar</code> syntax. First compute the AIC scores of both the <code>yule_models</code> and the <code>birth_death_models</code> we fitted above,</p>
<pre class="sourceCode r"><code class="sourceCode r">yule_aic &lt;- 
<span class="kw">sapply</span>(yule_models, function(pars)
                    <span class="dv">2</span> * (<span class="kw">length</span>(pars) - <span class="dv">1</span>) + <span class="dv">2</span> * pars[<span class="dv">1</span>] )
birth_death_aic &lt;- 
<span class="kw">sapply</span>(birth_death_models, function(pars)
                            <span class="dv">2</span> * (<span class="kw">length</span>(pars) - <span class="dv">1</span>) + <span class="dv">2</span> * pars[<span class="dv">1</span>] )</code></pre>
<p>And then generate a list identifying which model has the best (lowest) AIC score among the Yule models and which has the best AIC score among the birth-death models,</p>
<pre class="sourceCode r"><code class="sourceCode r">best_no_of_rates &lt;- <span class="kw">list</span>(<span class="dt">Yule =</span> <span class="kw">which.min</span>(yule_aic), 
                         <span class="dt">birth.death =</span> <span class="kw">which.min</span>(birth_death_aic))</code></pre>
<p>The best model is then whichever of these has the smaller AIC value.</p>
<pre class="sourceCode r"><code class="sourceCode r">best_model &lt;- <span class="kw">which.min</span>(<span class="kw">c</span>(<span class="kw">min</span>(yule_aic), <span class="kw">min</span>(birth_death_aic)))</code></pre>
<p>which confirms that the Yule 2-rate<br />model is still the best choice based on AIC score. Of the eight models in this second analysis, only three were in the original set considered (Yule 1-rate and 2-rate, and birth-death without a shift), so we could by no means have been sure ahead of time that a birth death with a shift, or a Yule model with a greater number of shifts, would not have fitted better.</p>
<h1 id="references">References</h1>
<p>Derryberry, Elizabeth P., Santiago Claramunt, Graham Derryberry, R. Terry Chesser, Joel Cracraft, Alexandre Aleixo, Jorge Pérez-Emán, J. V. Remsen Jr, and Robb T. Brumfield. 2011. “LINEAGE DIVERSIFICATION AND MORPHOLOGICAL EVOLUTION IN A LARGE-SCALE CONTINENTAL RADIATION: THE NEOTROPICAL OVENBIRDS AND WOODCREEPERS (AVES: FURNARIIDAE).” <em>Evolution</em> (jul). doi:10.1111/j.1558-5646.2011.01374.x. <a href="http://doi.wiley.com/10.1111/j.1558-5646.2011.01374.x" title="http://doi.wiley.com/10.1111/j.1558-5646.2011.01374.x">http://doi.wiley.com/10.1111/j.1558-5646.2011.01374.x</a>.</p>
<p>Rabosky, Daniel L. 2006. “LASER: a maximum likelihood toolkit for detecting temporal shifts in diversification rates from molecular phylogenies.” <em>Evolutionary bioinformatics online</em> 2 (jan): 273–6. <a href="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2674670\&amp;tool=pmcentrez\&amp;rendertype=abstract" title="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2674670\&amp;tool=pmcentrez\&amp;rendertype=abstract">http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2674670\&amp;tool=pmcentrez\&amp;rendertype=abstract</a>.</p>
<p>Stadler, Tanja. 2011. “Mammalian phylogeny reveals recent diversification rate shifts.” <em>Proceedings of the National Academy of Sciences</em> 2011 (mar). doi:10.1073/pnas.1016876108. <a href="http://www.pnas.org/cgi/doi/10.1073/pnas.1016876108" title="http://www.pnas.org/cgi/doi/10.1073/pnas.1016876108">http://www.pnas.org/cgi/doi/10.1073/pnas.1016876108</a>.</p>
<ol style="list-style-type: decimal">
<li><p>The TreeBASE portal is an important and rapidly growing repository of phylogenetic data. The R statistical environment has also become a primary tool for applied phylogenetic analyses across a range of questions, from comparative evolution to community ecology to conservation planning.</p></li>
<li><p>We have developed <code>treebase</code>, an open-source (freely available from <a href="http://cran.r-project.org/web/packages/treebase">http://cran.r-project.org/web/packages/treebase</a>) for the R programming environment, providing simplified, <em>programmatic</em> and interactive access to phylogenetic data in the TreeBASE repository.</p></li>
<li><p>We illustrate how this package creates a bridge between the TreeBASE repository and the rapidly growing collection of R packages for phylogenetics that can reduce barriers to discovery and integration across phylogenetic research.</p></li>
<li><p>We show how the <code>treebase</code> package can be used to facilitate replication of previous studies and testing of methods and hypotheses across a large sample of phylogenies, which may help make such important practices more common.</p></li>
</ol>
<h4 id="keywords">Keywords</h4>
<p>R, software, API, TreeBASE, database, programmatic, workflow</p>
<h1 id="introduction">Introduction</h1>
<p>Applications that use phylogenetic information as part of their analyses are becoming increasingly central to both evolutionary and ecological research. The exponential growth in genetic sequence data available for all forms of life has driven rapid advances in the methods that can infer the phylogenetic relationships and divergence times across different taxa <span class="citation">(Huelsenbeck and Ronquist 2001; Stamatakis 2006; Drummond and Rambaut 2007)</span>. Once again the product of one field has become the raw data of the next. Unfortunately, while the discipline of bioinformatics has emerged to help harness and curate the wealth of genetic data with cutting edge computer science, statistics, and Internet technology, its counterpart in evolutionary informatics remains “scattered, poorly documented, and in formats that impede discovery and integration” <span class="citation">(Parr et al. 2011)</span>. Our goal in developing the <code>treebase</code> package is to provide steps to reduce these challenges through programmatic and interactive access between the repositories that store this data and the software tools commonly used to analyse them.</p>
<p>The R statistical environment <span class="citation">(R Development Core Team 2012)</span> has become a dominant platform for researchers using phylogenetic data to address a rapidly expanding set of questions in ecological and evolutionary processes. These methods include, but are not limited to, ancestral state reconstruction <span class="citation">(Paradis 2004; Butler and King 2004)</span>, diversification analysis <span class="citation">(Paradis 2004; Rabosky 2006; Harmon et al. 2008)</span>, identifying trait dependent speciation and extinction rates, <span class="citation">(Fitzjohn 2010; Goldberg, Lancaster, and Ree 2011; Stadler 2011b)</span>, quantifying the rate and tempo of trait evolution <span class="citation">(Butler and King 2004; Harmon et al. 2008; Eastman et al. 2011)</span>, identifying evolutionary influences and proxies for community ecology <span class="citation">(Webb, Ackerly, and Kembel 2008; Kembel et al. 2010)</span>, connecting phylogeny data to climate patterns <span class="citation">(Warren, Glor, and Turelli 2008; Evans et al. 2009)</span>, and simulation of speciation and character evolution <span class="citation">(Harmon et al. 2008; Stadler 2011a; Boettiger, Coop, and Ralph 2012)</span>, as well as various manipulations and visualizations of phylogenetic data <span class="citation">(Paradis 2004; Schliep 2010; Jombart, Balloux, and Dray 2010; Revell et al. 2011)</span>. A more comprehensive list of R packages by analysis type is available on the phylogenetics taskview, <a href="http://cran.r-project.org/web/views/Phylogenetics.html">http://cran.r-project.org/web/views/Phylogenetics.html</a>. A few programs for applied phylogenetic methods are written for environments outside the R environment, incuding Java <span class="citation">(Maddison and Maddison 2011)</span>, MATLAB <span class="citation">(Blomberg, Garland, and Ives 2003)</span> and Python <span class="citation">(Sukumaran and Holder 2010)</span> and online interfaces <span class="citation">(Martins 2004)</span>.</p>
<p>TreeBASE (<a href="http://treebase.org">http://treebase.org</a>) is an online repository of phylogenetic data (e.g. trees of species, populations, or genes) that have been published in a peer-reviewed academic journal, book, thesis or conference proceedings <span class="citation">(Sanderson et al. 1994; Morell 1996)</span>. The database can be searched through an online interface which allows users to find a phylogenetic tree from a particular publication, author or taxa of interest. TreeBASE provides an application programming interface (API) that lets computer applications make queries to the database. Our <code>treebase</code> package uses this API to create a direct link between this data and the R environment. This has several immediate and important benefits:</p>
<ol style="list-style-type: decimal">
<li><p><em>Data discovery.</em> Users can leverage the rich, higher-level programming environment provided by the R environment to better identify data sets appropriate for their research by iteratively constructing queries for datasets that match appropriate metadata requirements.</p></li>
<li><p><em>Programmatic data access.</em> Many tasks that are theoretically made possible by the creation of the TreeBASE repository are not pursued because they would be too laborious for an exploratory analysis. The ability to use programmatic access across data sets to automatically download and perform a reproduciblye and systematic analysis using the rich set of tools available in R opens up new avenues for research.</p></li>
<li><p><em>Automatic updating</em>. The TreeBASE repository is expanding rapidly. The scriptable nature of analyses run with our <code>treebase</code> package means that a study can be rerun on the latest version of the repository without additional effort but with potential new information.</p></li>
</ol>
<h2 id="programmatic-web-access">Programmatic Web Access</h2>
<p>The Treebase repository makes data accessible by Web queries through a RESTful (REpresentational State Transfer) interface, which supplies search conditions in the address URL. The repository returns the requested data in XML (extensible markup language) format. The <code>treebase</code> package uses the <code>RCurl</code> package <span class="citation">(Lang 2012a)</span> to make queries over the Web to the repository, and the <code>XML</code> package <span class="citation">(Lang 2012b)</span> to parse the Web page returned by the repository into meaningful R data objects. While these querying and parsing functions comprise most of the code provided in the <code>treebase</code> package, they are hidden from the end user who can interact with these rich data retrieval and manipulation tools to access data from these remote repositories in much the same way as data locally available on the users hard-disk.</p>
<h2 id="basic-queries">Basic queries</h2>
<p>The <code>treebase</code> package allows these queries to be made directly from R, just as a user would make them from the Web browser. This enables a user to construct more complicated filters than permitted by the Web interface, and allows the user to maintain a record of the queries they used to collect their data as an R script. Scripting the data-gathering process helps reduce errors and assists in replicating the analysis later, either by the authors or other researchers <span class="citation">(Peng et al. 2011)</span>.</p>
<p>The <code>search_treebase</code> function forms the base of the <code>treebase</code> package. Table 1 lists each of the types of queries available through the <code>search_treebase</code> function. This list can also be found in the function documentation through the R command <code>?search_treebase</code>.<br />Any of the queries available on the Web interface can now be made directly from R, including downloading and importing a phylogeny into the R interface. For instance, one can search for phylogenies containing dolphin taxa, &quot;Delphinus,&quot; or all phylogenies submitted by a given author, &quot;Huelsenbeck&quot; using the R commands</p>
<pre class="sourceCode r"><code class="sourceCode r">    <span class="kw">search_treebase</span>(<span class="st">&quot;Delphinus&quot;</span>, <span class="dt">by=</span><span class="st">&quot;taxon&quot;</span>)
    <span class="kw">search_treebase</span>(<span class="st">&quot;Huelsenbeck&quot;</span>, <span class="dt">by=</span><span class="st">&quot;author&quot;</span>)</code></pre>
<p>This function returns the matching phylogenies into R as an R object, ready for analysis. The package documentation provides many examples of possible queries.</p>
<table>
<caption>Queries available in <code>search_treebase</code>. The first argument is the keyword used in the query such as an author's name and the second argument indicates the type of query (<em>i.e.</em> &quot;author&quot;).</caption>
<thead>
<tr class="header">
<th align="left">search &quot;by=&quot;</th>
<th align="left">purpose</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">abstract</td>
<td align="left">search terms in the publication abstract</td>
</tr>
<tr class="even">
<td align="left">author</td>
<td align="left">match authors in the publication</td>
</tr>
<tr class="odd">
<td align="left">subject</td>
<td align="left">Matches in the subject terms</td>
</tr>
<tr class="even">
<td align="left">doi</td>
<td align="left">The unique object identifier for the publication</td>
</tr>
<tr class="odd">
<td align="left">ncbi</td>
<td align="left">NCBI identifier number for the taxon</td>
</tr>
<tr class="even">
<td align="left">kind.tree</td>
<td align="left">Kind of tree (Gene tree, species tree, barcode tree)</td>
</tr>
<tr class="odd">
<td align="left">type.tree</td>
<td align="left">Type of tree (Consensus or Single)</td>
</tr>
<tr class="even">
<td align="left">ntax</td>
<td align="left">Number of taxa in the matrix</td>
</tr>
<tr class="odd">
<td align="left">quality</td>
<td align="left">A quality score for the tree, if it has been rated.</td>
</tr>
<tr class="even">
<td align="left">study</td>
<td align="left">Match words in the title of the study or publication</td>
</tr>
<tr class="odd">
<td align="left">taxon</td>
<td align="left">Taxon scientific name</td>
</tr>
<tr class="even">
<td align="left">id.study</td>
<td align="left">TreeBASE study ID</td>
</tr>
<tr class="odd">
<td align="left">id.tree</td>
<td align="left">TreeBASE's unique tree identifier (Tr.id)</td>
</tr>
<tr class="even">
<td align="left">id.taxon</td>
<td align="left">Taxon identifier number from TreeBase</td>
</tr>
<tr class="odd">
<td align="left">tree</td>
<td align="left">The title for the tree</td>
</tr>
</tbody>
</table>
<h2 id="accessing-all-phylogenies">Accessing all phylogenies</h2>
<p>For certain applications a user may wish to download all the available phylogenies from TreeBASE. Using the <code>cache_treebase</code> function allows a user to download a local copy of all trees. Because direct database dumps are not available from treebase.org, this function has intentional delays to avoid overtaxing the TreeBASE servers, and should be allowed a full day to run.</p>
<pre class="sourceCode r"><code class="sourceCode r">treebase &lt;- <span class="kw">cache_treebase</span>()</code></pre>
<p>Once run, the cache is saved compactly in memory where it can be easily and quickly restored. For convenience, the <code>treebase</code> package comes with a copy already cached, which can be loaded into memory.</p>
<pre class="sourceCode r"><code class="sourceCode r"><span class="kw">data</span>(treebase)</code></pre>
<p>All of the examples shown in this manuscript are run as shown using the <code>knitr</code> package for authoring dynamic documents <span class="citation">(Xie 2012)</span>, which helps ensure the results shown are reproducible. These examples can be updated by copying and pasting the code shown into the R terminal, or by recompiling the entire manuscript from the source files found on the development Web page for the TreeBASE package, <a href="https://github.com/ropensci/treebase">github.com/ropensci/treebase</a>. Data was accessed to produce the examples shown on Wed Jun 27 11:01:42 2012.</p>
<h1 id="data-discovery-in-treebase">Data discovery in TreeBASE</h1>
<p>Data discovery involves searching for existing data that meets certain desired characteristics. Such searches take advantage of metadata -- summary information describing the data entries provided in the repository. The Web repository uses separate interfaces (APIs) to access metadata describing the publications associated with the data entered, such as the publisher, year of publication, etc., and a different interface to describe the metadata associated with an individual phylogeny, such as the number of taxa or the kind of tree (<em>e.g.</em> Gene tree versus Species tree). The <code>treebase</code> package can query these individual sources of metadata separately, but this information is most powerful when used in concert -- allowing the construction of complicated searches that cannot be automated through the Web interface. The <code>metadata</code> function updates a list of all available metadata from both APIs and returns this information as an R <code>data.frame</code>.</p>
<pre class="sourceCode r"><code class="sourceCode r">meta &lt;- <span class="kw">metadata</span>()</code></pre>
<p>From the length of the metadata list we see that there are currently 3164 published studies in the database.</p>
<p>The fields provided by <code>metadata</code> are listed in Table II.</p>
<table>
<caption>Columns of metadata available from the <code>metadata</code> function</caption>
<thead>
<tr class="header">
<th align="left">metadata field</th>
<th align="left">description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Study.id</td>
<td align="left">TreeBASE study ID</td>
</tr>
<tr class="even">
<td align="left">Tree.id</td>
<td align="left">TreeBASE's unique tree identifier</td>
</tr>
<tr class="odd">
<td align="left">kind</td>
<td align="left">Kind of tree (Gene tree, species tree, barcode tree)</td>
</tr>
<tr class="even">
<td align="left">type</td>
<td align="left">Type of tree (Consensus or Single)</td>
</tr>
<tr class="odd">
<td align="left">quality</td>
<td align="left">A quality score for the tree, if it has been rated.</td>
</tr>
<tr class="even">
<td align="left">ntaxa</td>
<td align="left">Number of taxa in the matrix</td>
</tr>
<tr class="odd">
<td align="left">date</td>
<td align="left">Year the study was published</td>
</tr>
<tr class="even">
<td align="left">author</td>
<td align="left">First author in the publication</td>
</tr>
<tr class="odd">
<td align="left">title</td>
<td align="left">The title of the publication</td>
</tr>
</tbody>
</table>
<p>Metadata can also be used to reveal trends in the data deposition which may be useful in identifying patterns or biases in research or emerging potential types of data. As a simple example, we look at trends in the submission patterns of publishers over time,</p>
<pre class="sourceCode r"><code class="sourceCode r">    date &lt;- meta[[<span class="st">&quot;date&quot;</span>]] 
    pub &lt;- meta[[<span class="st">&quot;publisher&quot;</span>]]</code></pre>
<p>Many journals have only a few submissions, so we will label any not in the top ten contributing journals as “Other”:</p>
<pre class="sourceCode r"><code class="sourceCode r">    topten &lt;- <span class="kw">sort</span>(<span class="kw">table</span>(pub), <span class="dt">decreasing=</span><span class="ot">TRUE</span>)[<span class="dv">1</span>:<span class="dv">10</span>]
    meta[[<span class="st">&quot;publisher&quot;</span>]][!(pub %in% <span class="kw">names</span>(topten))] &lt;- <span class="st">&quot;Other&quot;</span></code></pre>
<p>We plot the distribution of publication years for phylogenies deposited in TreeBASE, color coding by publisher in Fig [fig:1].</p>
<pre class="sourceCode r"><code class="sourceCode r">  <span class="kw">library</span>(ggplot2) 
  <span class="kw">ggplot</span>(meta) + <span class="kw">geom_bar</span>(<span class="kw">aes</span>(date, <span class="dt">fill =</span> publisher)) </code></pre>
<div class="figure">
<img src="http://farm9.staticflickr.com/8018/7451038138_6da7f23d59_o.png" alt="Histogram of publication dates by year, with the code required to generate the figure." /><p class="caption">Histogram of publication dates by year, with the code required to generate the figure.</p>
</div>
<p>Typically we are interested in the metadata describing the phylogenies themselves rather than just in the publications in which they appeared. Phylogenetic metadata includes features such as the number of taxa in the tree, a quality score (if available), kind of tree (gene tree, species tree, or barcode tree) or whether the phylogeny represents a consensus tree from a distribution or just a single estimate.</p>
<p>Even simple queries can illustrate the advantage of interacting with TreeBASE data through an R interface has over the Web interface. A Web interface can only perform the tasks built in by design. For instance, rather than performing six separate searches to determine the number of consensus vs single phylogenies available for each king of tree, we can construct a 2 by 2 table with a single line of code,</p>
<pre class="sourceCode r"><code class="sourceCode r"><span class="kw">table</span>(meta[[<span class="st">&quot;kind&quot;</span>]], meta[[<span class="st">&quot;type&quot;</span>]])</code></pre>
<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Jun 25 13:18:20 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> 
Consensus
</TH> <TH> 
Single
</TH>  </TR>
  <TR> <TD align="right"> 
Barcode Tree
</TD> <TD align="right">   
1
</TD> <TD align="right">   
4
</TD> </TR>
  <TR> <TD align="right"> 
Gene Tree
</TD> <TD align="right">  
65
</TD> <TD align="right"> 
134
</TD> </TR>
  <TR> <TD align="right"> 
Species Tree
</TD> <TD align="right"> 
2863
</TD> <TD align="right"> 
5857
</TD> </TR>
   </TABLE>




<h1 id="reproducible-computations">Reproducible computations</h1>
<p>Reproducible research has become a topic of increasing interest in recent years, and facilitating access to data and using scripts that can replicate analyses can help lower barriers to the replication of statistical and computational results <span class="citation">(Schwab, Karrenbach, and Claerbout 2000; Gentleman and Temple Lang 2004; Peng 2011)</span>. The <code>treebase</code> package facilitates this process, as we illustrate in a simple example.</p>
<p>Consider the shifts in speciation rate identified by Derryberry et al. <span class="citation">(2011)</span> on a phylogeny of ovenbirds and treecreepers. We will seek to not only replicate the results the authors obtained by fitting the models provided in the R package <code>laser</code> <span class="citation">(Rabosky 2006)</span>, but also compare them against methods presented in Stadler <span class="citation">(2011b)</span> and implemented in the package <code>TreePar</code>, which permits speciation models that were not available to Derryberry et al. <span class="citation">(2011)</span> at the time of their study.</p>
<h2 id="obtaining-the-tree">Obtaining the tree</h2>
<p>By drawing on the rich data manipulation tools available in R which may be familiar to the large R phylogenetics community, the <code>treebase</code> package allows us to construct richer queries than are possible through the TreeBASE Web interface alone.</p>
<p>The most expedient way to identify the data uses the digital object identifer (doi) at the top of most articles, which we use in a call to the <code>search_treebase</code> function, such as</p>
<pre class="sourceCode r"><code class="sourceCode r">results &lt;- <span class="kw">search_treebase</span>(<span class="st">&quot;10.1111/j.1558-5646.2011.01374.x&quot;</span>, <span class="st">&quot;doi&quot;</span>) </code></pre>
<p>The search returns a list, since some publications can contain many trees. In this case our phylogeny is in the first element of the list.</p>
<p>Having imported the phylogenetic tree corresponding to this study, we can quickly replicate their analysis of which diversification process best fits the data. These steps can be easily implemented using the phylogenetics packages we have just mentioned.</p>
<p>For instance, we can calculate the branching times of each node on the phylogeny,</p>
<pre class="sourceCode r"><code class="sourceCode r">bt &lt;- <span class="kw">branching.times</span>(derryberry)</code></pre>
<p>and then begin to fit each model the authors have tested, such as the pure birth model,</p>
<pre class="sourceCode r"><code class="sourceCode r">yule = <span class="kw">pureBirth</span>(bt)</code></pre>
<p>or the birth-death model,</p>
<pre class="sourceCode r"><code class="sourceCode r">birth_death = <span class="kw">bd</span>(bt)</code></pre>
<p>The estimated models are now loaded into the active R session where we can further explore them as we go along. The appendix shows the estimation and comparison of all the models originally considered by Derryberry et al. <span class="citation">(2011)</span>.</p>
<p>In this fast-moving field, new methods often become available between the time of submission and time of publication of a manuscript. For instance, the more sophisticated models introduced in Stadler <span class="citation">(2011b)</span> were not used in this study, but have since been made available in the recent package, <code>TreePar</code>. These richer models permit a shift the speciation or extinction rate to occur multiple times throughout the course of the phylogeny.</p>
<p>We load the new method and format the phylogeny using the R commands:</p>
<pre class="sourceCode r"><code class="sourceCode r"><span class="kw">library</span>(TreePar)
x &lt;- <span class="kw">sort</span>(<span class="kw">getx</span>(derryberry), <span class="dt">decreasing =</span> <span class="ot">TRUE</span>)</code></pre>
<p>Here we consider models that have up to 4 different rates in Yule models, (The syntax in <code>TreePar</code> is slightly cumbersome, the [[2]] indicates where this command happens to store the output models.)</p>
<p>As a comparison of speciation models is not the focus of this paper, the complete code and explanation for these steps is provided as an appendix. Happily, this analysis confirms the author's original conclusions, even when the more general models of Stadler <span class="citation">(2011b)</span> are considered.</p>
<h1 id="analyses-across-many-phylogenies">Analyses across many phylogenies</h1>
<p>Large scale comparative analyses that seek to characterize evolutionary patterns across many phylogenies are increasingly common in phylogenetic methods <span class="citation">(<em>e.g.</em> McPeek and Brown 2007; Phillimore and Price 2008; McPeek 2008; Quental and Marshall 2010; Davies et al. 2011)</span>. Sometimes referred to by their authors as meta-analyses, these approaches have focused on re-analyzing phylogenetic trees collected from many different earlier publications. This is a more direct approach than the traditional concept of meta-analysis where statistical results from earlier studies are weighted by their sample size without being able to access the raw data. Because the identical analysis can be repeated on the original data from each study, this approach avoids some of the statistical challenges inherent in traditional meta-analyses summarizing results across heterogeneous approaches.</p>
<p>To date, researchers have gone through heroic efforts simply to assemble these data sets from the literature. As described in McPeek and Brown <span class="citation">(2007)</span>; (emphasis added)</p>
<blockquote>
<p>One data set was based on 163 published species-level molecular phylogenies of arthropods, chordates, and mollusks. A PDF format file of each article was obtained, and a digital snapshot of the figure was taken in Adobe Acrobat 7.0. This image was transferred to a PowerPoint (Microsoft) file and printed on a laser printer. The phylogenies included in this study are listed in the appendix. <em>All branch lengths were measured by hand from these printed sheets using dial calipers.</em></p>
</blockquote>
<p>Despite the recent emergence of digital tools that could now facilitate this analysis without mechanical calipers, <span class="citation">(<em>e.g.</em> treesnatcher, Laubach and von Haeseler 2007)</span>, it is easier and less error-prone to pull properly formatted phylogenies from the database for this purpose. Moreover, as the available data increases with subsequent publications, updating earlier meta-analyses can become increasingly tedious. Using <code>treebase</code>, a user can apply any analysis they have written for a single phylogeny across the entire collection of suitable phylogenies in TreeBASE, which can help overcome such barriers to discovery and integration at this large scale. Using the functions we introduce aboved, we provide a simple example that computes the gamma statistic of Pybus and Harvey <span class="citation">(2000)</span>, which provides an measure of when speciation patterns differ from the popular birth-death model.</p>
<h2 id="tests-across-many-phylogenies">Tests across many phylogenies</h2>
<p>A standard test of this is the gamma statistic of Pybus and Harvey <span class="citation">(2000)</span> which tests the null hypothesis that the rates of speciation and extinction are constant. The gamma statistic is normally distributed about 0 for a pure birth or birth-death process, values larger than 0 indicate that internal nodes are closer to the tip then expected, while values smaller than 0 indicate nodes farther from the tip then expected. In this section, we collect all phylogenetic trees from TreeBASE and select those with branch length data that we can time-calibrate using tools available in R. We can then calculate the distribution of this statistic for all available trees, and compare these results with those from the analyses mentioned above.</p>
<p>The <code>treebase</code> package provides a compressed cache of the phylogenies available in treebase. This cache can be automatically updated with the <code>cache_treebase</code> function,</p>
<pre class="sourceCode r"><code class="sourceCode r">treebase &lt;- <span class="kw">cache_treebase</span>()</code></pre>
<p>which may require a day or so to complete, and will save a file in the working directory named with treebase and the date obtained. For convenience, we can load the cached copy distributed with the <code>treebase</code> package:</p>
<pre class="sourceCode r"><code class="sourceCode r"><span class="kw">data</span>(treebase)</code></pre>
<p>We will only be able to use those phylogenies that include branch length data. We drop those that do not from the data set,</p>
<pre class="sourceCode r"><code class="sourceCode r">      have &lt;- <span class="kw">have_branchlength</span>(treebase)
      branchlengths &lt;- treebase[have]</code></pre>
<p>Like most comparative methods, this analysis will require ultrametric trees (branch lengths proportional to time, rather than to mutational steps). As most of these phylogenies are calibrated with branch length proportional to mutational step, we must time-calibrate each of them first.</p>
<pre class="sourceCode r"><code class="sourceCode r">timetree &lt;- function(tree)
    <span class="kw">try</span>( <span class="kw">chronoMPL</span>(<span class="kw">multi2di</span>(tree)) )
tt &lt;- <span class="kw">drop_nontrees</span>(<span class="kw">sapply</span>(branchlengths, timetree))</code></pre>
<p>At this point we have 1,396 time-calibrated phylogenies over which we will apply the diversification rate analysis. Computing the gamma test statistic to identify deviations from the constant-rates model takes a single line,</p>
<pre class="sourceCode r"><code class="sourceCode r">gammas &lt;- <span class="kw">sapply</span>(tt,  gammaStat)</code></pre>
<p>and the resulting distribution of the statistic across available trees is shown Fig 2. While researchers have often considered this statistic for individual phylogenies, we are unaware of any study that has visualized the empirical distribution of this statistic across thousands of phylogenies. Both the overall distribution, which appears slightly skewed towards positive values indicating increasing rate of speciation near the tips, and the position and identity of outlier phylogenies are patterns that may introduce new hypotheses and potential directions for further exploration.</p>
<pre class="sourceCode r"><code class="sourceCode r"><span class="kw">qplot</span>(gammas)+<span class="kw">xlab</span>(<span class="st">&quot;gamma statistic&quot;</span>)</code></pre>
<div class="figure">
<img src="http://farm8.staticflickr.com/7111/7455801392_4a47e33e8e_o.png" alt="Distribution of the gamma statistic across phylogenies in TreeBASE. Strongly positive values are indicative of an increasing rate of evolution (excess of nodes near the tips), very negative values indicate an early burst of diversification (an excess of nodes near the root)." /><p class="caption">Distribution of the gamma statistic across phylogenies in TreeBASE. Strongly positive values are indicative of an increasing rate of evolution (excess of nodes near the tips), very negative values indicate an early burst of diversification (an excess of nodes near the root).</p>
</div>
<h1 id="conclusion">Conclusion</h1>
<p>While we have focused on examples that require no additional data beyond the phylogeny, a wide array of methods combine this data with information about the traits, geography, or ecological community of the taxa represented. In such cases we would need programmatic access to the trait data as well as the phylogeny. The Dryad digital repository (<a href="http://datadryad.org">http://datadryad.org</a>) is an effort in this direction. While programmatic access to the repository is possible through the <code>rdryad</code> package <span class="citation">(Chamberlain, Boettiger, and Ram 2012)</span>, variation in data formatting must first be overcome before similar direct access to the data is possible. Dedicated databases such as FishBASE (<a href="http://fishbase.org">http://fishbase.org</a>) may be another alternative, where morphological data can be queried for a list of species using the <code>rfishbase</code> package <span class="citation">(Boettiger)</span>. The development of similar software for programmatic data access will rapidly extend the space and scale of possible analyses.</p>
<p>The recent advent of mandatory data archiving in many of the major journals publishing phylognetics-based research <span class="citation">(<em>e.g.</em> Fairbairn 2010; Piwowar, Vision, and Whitlock 2011; Whitlock et al. 2010)</span>, is a particularly promising development that should continue to fuel the trend of submissions seen in Fig. 1. Accompanied by faster and more inexpensive techniques of NextGen sequencing, and the rapid expansion in phylogenetic applications, we anticipate this rapid growth in available phylogenies will continue. Faced with this flood of data, programmatic access becomes not only increasingly powerful but an increasingly necessary way to ensure we can still see the forest for all the trees.</p>
<h1 id="acknowledgements">Acknowledgements</h1>
<p>CB wishes to thank S. Price for feedback on the manuscript, the TreeBASE developer team for building and supporting the repository, and all contributers to TreeBASE. CB is supported by a Computational Sciences Graduate Fellowship from the Department of Energy under grant number DE-FG02-97ER25308.</p>
<h1 id="references">References</h1>
<p>Blomberg, S. P., JR Theodore Garland, and A. R. Ives. 2003. “Testing for phylogenetic signal in comparative data: behavioral traits are more labile.” <em>Evolution</em> 57: 717–745. <a href="http://www3.interscience.wiley.com/journal/118867878/abstract" title="http://www3.interscience.wiley.com/journal/118867878/abstract">http://www3.interscience.wiley.com/journal/118867878/abstract</a>.</p>
<p>Boettiger, Carl. “rfishbase: R Interface to FishBASE.”</p>
<p>Boettiger, Carl, Graham Coop, and Peter Ralph. 2012. “Is your phylogeny informative? Measuring the power of comparative methods.” <em>Evolution</em> (jan). doi:10.1111/j.1558-5646.2012.01574.x. <a href="http://doi.wiley.com/10.1111/j.1558-5646.2012.01574.x" title="http://doi.wiley.com/10.1111/j.1558-5646.2012.01574.x">http://doi.wiley.com/10.1111/j.1558-5646.2012.01574.x</a>.</p>
<p>Butler, Marguerite A., and Aaron A. King. 2004. “Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution.” <em>The American Naturalist</em> 164 (dec): 683–695. doi:10.1086/426002. <a href="http://www.jstor.org/stable/10.1086/426002" title="http://www.jstor.org/stable/10.1086/426002">http://www.jstor.org/stable/10.1086/426002</a>.</p>
<p>Chamberlain, Scott, Carl Boettiger, and Karthik Ram. 2012. “rdryad: Dryad API interface.” <a href="http://www.github.com/ropensci/rdryad " title="http://www.github.com/ropensci/rdryad ">http://www.github.com/ropensci/rdryad </a>.</p>
<p>Davies, T. Jonathan, Andrew P. Allen, Luís Borda-de-Água, Jim Regetz, and Carlos J. Melián. 2011. “NEUTRAL BIODIVERSITY THEORY CAN EXPLAIN THE IMBALANCE OF PHYLOGENETIC TREES BUT NOT THE TEMPO OF THEIR DIVERSIFICATION.” <em>Evolution</em> 65 (jul): 1841–1850. doi:10.1111/j.1558-5646.2011.01265.x. <a href="http://doi.wiley.com/10.1111/j.1558-5646.2011.01265.x http://www.ncbi.nlm.nih.gov/pubmed/21729042" title="http://doi.wiley.com/10.1111/j.1558-5646.2011.01265.x http://www.ncbi.nlm.nih.gov/pubmed/21729042">http://doi.wiley.com/10.1111/j.1558-5646.2011.01265.x http://www.ncbi.nlm.nih.gov/pubmed/21729042</a>.</p>
<p>Derryberry, Elizabeth P., Santiago Claramunt, Graham Derryberry, R. Terry Chesser, Joel Cracraft, Alexandre Aleixo, Jorge Pérez-Emán, J. V. Remsen Jr, and Robb T. Brumfield. 2011. “LINEAGE DIVERSIFICATION AND MORPHOLOGICAL EVOLUTION IN A LARGE-SCALE CONTINENTAL RADIATION: THE NEOTROPICAL OVENBIRDS AND WOODCREEPERS (AVES: FURNARIIDAE).” <em>Evolution</em> (jul). doi:10.1111/j.1558-5646.2011.01374.x. <a href="http://doi.wiley.com/10.1111/j.1558-5646.2011.01374.x" title="http://doi.wiley.com/10.1111/j.1558-5646.2011.01374.x">http://doi.wiley.com/10.1111/j.1558-5646.2011.01374.x</a>.</p>
<p>Drummond, Alexei J., and Andrew Rambaut. 2007. “BEAST: Bayesian evolutionary analysis by sampling trees.” <em>BMC evolutionary biology</em> 7 (jan): 214. doi:10.1186/1471-2148-7-214. <a href="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2247476\&amp;tool=pmcentrez\&amp;rendertype=abstract" title="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2247476\&amp;tool=pmcentrez\&amp;rendertype=abstract">http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2247476\&amp;tool=pmcentrez\&amp;rendertype=abstract</a>.</p>
<p>Eastman, Jonathan M., Michael E. Alfaro, Paul Joyce, Andrew L. Hipp, and Luke J. Harmon. 2011. “A NOVEL COMPARATIVE METHOD FOR IDENTIFYING SHIFTS IN THE RATE OF CHARACTER EVOLUTION ON TREES.” <em>Evolution</em> 65 (jul): 3578–3589. doi:10.1111/j.1558-5646.2011.01401.x. <a href="http://doi.wiley.com/10.1111/j.1558-5646.2011.01401.x" title="http://doi.wiley.com/10.1111/j.1558-5646.2011.01401.x">http://doi.wiley.com/10.1111/j.1558-5646.2011.01401.x</a>.</p>
<p>Evans, Margaret E. K., Stephen a Smith, Rachel S. Flynn, and Michael J. Donoghue. 2009. “Climate, niche evolution, and diversification of the ‘bird-cage’ evening primroses (Oenothera, sections Anogra and Kleinia).” <em>The American naturalist</em> 173 (feb): 225–40. doi:10.1086/595757. <a href="http://www.ncbi.nlm.nih.gov/pubmed/19072708" title="http://www.ncbi.nlm.nih.gov/pubmed/19072708">http://www.ncbi.nlm.nih.gov/pubmed/19072708</a>.</p>
<p>Fairbairn, Daphne J. 2010. “THE ADVENT OF MANDATORY DATA ARCHIVING.” <em>Evolution</em> (nov). doi:10.1111/j.1558-5646.2010.01182.x. <a href="http://doi.wiley.com/10.1111/j.1558-5646.2010.01182.x" title="http://doi.wiley.com/10.1111/j.1558-5646.2010.01182.x">http://doi.wiley.com/10.1111/j.1558-5646.2010.01182.x</a>.</p>
<p>Fitzjohn, Richard G. 2010. “Quantitative Traits and Diversification.” <em>Systematic biology</em> 59 (sep): 619–633. doi:10.1093/sysbio/syq053. <a href="http://www.ncbi.nlm.nih.gov/pubmed/20884813" title="http://www.ncbi.nlm.nih.gov/pubmed/20884813">http://www.ncbi.nlm.nih.gov/pubmed/20884813</a>.</p>
<p>Gentleman, Robert, and D. Temple Lang. 2004. “Statistical analyses and reproducible research.” <em>Bioconductor Project Working Papers</em>: 2. <a href="http://www.bepress.com/cgi/viewcontent.cgi?article=1001\&amp;amp;context=bioconductor" title="http://www.bepress.com/cgi/viewcontent.cgi?article=1001\&amp;amp;context=bioconductor">http://www.bepress.com/cgi/viewcontent.cgi?article=1001\&amp;amp;context=bioconductor</a>.</p>
<p>Goldberg, Emma E., Lesley T. Lancaster, and Richard H. Ree. 2011. “Phylogenetic Inference of Reciprocal Effects between Geographic Range Evolution and Diversification.” <em>Systematic biology</em> 60 (may): 451–465. doi:10.1093/sysbio/syr046. <a href="http://www.ncbi.nlm.nih.gov/pubmed/21551125" title="http://www.ncbi.nlm.nih.gov/pubmed/21551125">http://www.ncbi.nlm.nih.gov/pubmed/21551125</a>.</p>
<p>Harmon, Luke J., Jason T. Weir, Chad D. Brock, Richard E. Glor, and Wendell Challenger. 2008. “Geiger: investigating evolutionary radiations.” <em>Bioinformatics</em> 24: 129–131. doi:10.1093/bioinformatics/btm538.</p>
<p>Huelsenbeck, John P., and Fredrik Ronquist. 2001. “MRBAYES: Bayesian inference of phylogenetic trees.” <em>Bioinformatics (Oxford, England)</em> 17 (aug): 754–5. doi:10.1093/bioinformatics/17.8.754. <a href="http://www.ncbi.nlm.nih.gov/pubmed/11524383" title="http://www.ncbi.nlm.nih.gov/pubmed/11524383">http://www.ncbi.nlm.nih.gov/pubmed/11524383</a>.</p>
<p>Jombart, Thibaut, François Balloux, and Stéphane Dray. 2010. “Adephylo: New Tools for Investigating the Phylogenetic Signal in Biological Traits.” <em>Bioinformatics (Oxford, England)</em> 26 (aug): 1907–9. doi:10.1093/bioinformatics/btq292. <a href="http://www.ncbi.nlm.nih.gov/pubmed/20525823" title="http://www.ncbi.nlm.nih.gov/pubmed/20525823">http://www.ncbi.nlm.nih.gov/pubmed/20525823</a>.</p>
<p>Kembel, Steven W., Peter D. Cowan, Matthew R. Helmus, William K. Cornwell, Helene Morlon, David D. Ackerly, Simon P. Blomberg, and Campbell O. Webb. 2010. “Picante: R tools for integrating phylogenies and ecology.” <em>Bioinformatics (Oxford, England)</em> 26 (jun): 1463–4. doi:10.1093/bioinformatics/btq166. <a href="http://www.ncbi.nlm.nih.gov/pubmed/20395285" title="http://www.ncbi.nlm.nih.gov/pubmed/20395285">http://www.ncbi.nlm.nih.gov/pubmed/20395285</a>.</p>
<p>Lang, Duncan Temple. 2012a. “RCurl: General network (HTTP/FTP/...) client interface for R.” <a href="http://cran.r-project.org/package=RCurl" title="http://cran.r-project.org/package=RCurl">http://cran.r-project.org/package=RCurl</a>.</p>
<p>———. 2012b. “XML: Tools for parsing and generating XML within R and S-Plus.” <a href="http://cran.r-project.org/package=XML" title="http://cran.r-project.org/package=XML">http://cran.r-project.org/package=XML</a>.</p>
<p>Laubach, Thomas, and Arndt von Haeseler. 2007. “TreeSnatcher: coding trees from images.” <em>Bioinformatics (Oxford, England)</em> 23 (dec): 3384–5. doi:10.1093/bioinformatics/btm438. <a href="http://www.ncbi.nlm.nih.gov/pubmed/17893085" title="http://www.ncbi.nlm.nih.gov/pubmed/17893085">http://www.ncbi.nlm.nih.gov/pubmed/17893085</a>.</p>
<p>Maddison, W. P., and D. R. Maddison. 2011. “Mesquite: a modular system for evolutionary analysis.” <a href="http://mesquiteproject.org" title="http://mesquiteproject.org">http://mesquiteproject.org</a>.</p>
<p>Martins, E. P. 2004. “COMPARE, version Computer programs for the statistical analysis of comparative data.” Bloomington IN.: Department of Biology, Indiana University. <a href="http://compare.bio.indiana.edu/" title="http://compare.bio.indiana.edu/">http://compare.bio.indiana.edu/</a>.</p>
<p>McPeek, Mark a. 2008. “The ecological dynamics of clade diversification and community assembly.” <em>The American naturalist</em> 172 (dec): 270. doi:10.1086/593137. <a href="http://www.ncbi.nlm.nih.gov/pubmed/18851684" title="http://www.ncbi.nlm.nih.gov/pubmed/18851684">http://www.ncbi.nlm.nih.gov/pubmed/18851684</a>.</p>
<p>McPeek, Mark a, and Jonathan M. Brown. 2007. “Clade age and not diversification rate explains species richness among animal taxa.” <em>The American naturalist</em> 169 (apr): 97. doi:10.1086/512135. <a href="http://www.ncbi.nlm.nih.gov/pubmed/17427118" title="http://www.ncbi.nlm.nih.gov/pubmed/17427118">http://www.ncbi.nlm.nih.gov/pubmed/17427118</a>.</p>
<p>Morell, V. 1996. “TreeBASE: the roots of phylogeny.” <em>Science</em> 273: 569. doi:10.1126/science.273.5275.569. <a href="http://www.sciencemag.org/cgi/doi/10.1126/science.273.5275.569" title="http://www.sciencemag.org/cgi/doi/10.1126/science.273.5275.569">http://www.sciencemag.org/cgi/doi/10.1126/science.273.5275.569</a>.</p>
<p>Paradis, Emmanuel. 2004. “APE: Analyses of Phylogenetics and Evolution in R language.” <em>Bioinformatics</em> 20: 289–290. doi:10.1093/bioinformatics/btg412. <a href="http://www.bioinformatics.oupjournals.org/cgi/doi/10.1093/bioinformatics/btg412" title="http://www.bioinformatics.oupjournals.org/cgi/doi/10.1093/bioinformatics/btg412">http://www.bioinformatics.oupjournals.org/cgi/doi/10.1093/bioinformatics/btg412</a>.</p>
<p>Parr, Cynthia S., Robert Guralnick, Nico Cellinese, and Roderic D. M. Page. 2011. “Evolutionary informatics: unifying knowledge about the diversity of life.” <em>Trends in ecology &amp; evolution</em> 27 (dec): 94–103. doi:10.1016/j.tree.2011.11.001. <a href="http://www.ncbi.nlm.nih.gov/pubmed/22154516" title="http://www.ncbi.nlm.nih.gov/pubmed/22154516">http://www.ncbi.nlm.nih.gov/pubmed/22154516</a>.</p>
<p>Peng, Changhui, Joel Guiot, Haibin Wu, Hong Jiang, and Yiqi Luo. 2011. “Integrating models with data in ecology and palaeoecology: advances towards a model-data fusion approach.” <em>Ecology letters</em> (mar). doi:10.1111/j.1461-0248.2011.01603.x. <a href="http://www.ncbi.nlm.nih.gov/pubmed/21366814" title="http://www.ncbi.nlm.nih.gov/pubmed/21366814">http://www.ncbi.nlm.nih.gov/pubmed/21366814</a>.</p>
<p>Peng, R. D. 2011. “Reproducible Research in Computational Science.” <em>Science</em> 334 (dec): 1226–1227. doi:10.1126/science.1213847. <a href="http://www.sciencemag.org/cgi/doi/10.1126/science.1213847" title="http://www.sciencemag.org/cgi/doi/10.1126/science.1213847">http://www.sciencemag.org/cgi/doi/10.1126/science.1213847</a>.</p>
<p>Phillimore, Albert B., and Trevor D. Price. 2008. “Density-dependent cladogenesis in birds.” <em>PLoS biology</em> 6 (mar): 71. doi:10.1371/journal.pbio.0060071. <a href="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2270327\&amp;tool=pmcentrez\&amp;rendertype=abstract" title="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2270327\&amp;tool=pmcentrez\&amp;rendertype=abstract">http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2270327\&amp;tool=pmcentrez\&amp;rendertype=abstract</a>.</p>
<p>Piwowar, Heather A., Todd J. Vision, and Michael C. Whitlock. 2011. “Data archiving is a good investment.” <em>Nature</em> 473 (may): 285–285. doi:10.1038/473285a. <a href="http://www.nature.com/doifinder/10.1038/473285a" title="http://www.nature.com/doifinder/10.1038/473285a">http://www.nature.com/doifinder/10.1038/473285a</a>.</p>
<p>Pybus, O. G., and P. H. Harvey. 2000. “Testing macro-evolutionary models using incomplete molecular phylogenies.” <em>Proceedings of The Royal Society B</em> 267 (nov): 2267–72. doi:10.1098/rspb.2000.1278. <a href="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1690817\&amp;tool=pmcentrez\&amp;rendertype=abstract" title="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1690817\&amp;tool=pmcentrez\&amp;rendertype=abstract">http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1690817\&amp;tool=pmcentrez\&amp;rendertype=abstract</a>.</p>
<p>Quental, Tiago B., and Charles R. Marshall. 2010. “Diversity dynamics: molecular phylogenies need the fossil record.” <em>Trends in Ecology &amp; Evolution</em> (jun): 1–8. doi:10.1016/j.tree.2010.05.002. <a href="http://linkinghub.elsevier.com/retrieve/pii/S0169534710001011" title="http://linkinghub.elsevier.com/retrieve/pii/S0169534710001011">http://linkinghub.elsevier.com/retrieve/pii/S0169534710001011</a>.</p>
<p>R Development Core Team, The. 2012. “R: A language and environment for statistical computing.” Vienna, Austria: R Foundation for Statistical Computing. <a href="http://www.r-project.org/" title="http://www.r-project.org/">http://www.r-project.org/</a>.</p>
<p>Rabosky, Daniel L. 2006. “LASER: a maximum likelihood toolkit for detecting temporal shifts in diversification rates from molecular phylogenies.” <em>Evolutionary bioinformatics online</em> 2 (jan): 273–6. <a href="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2674670\&amp;tool=pmcentrez\&amp;rendertype=abstract" title="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2674670\&amp;tool=pmcentrez\&amp;rendertype=abstract">http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2674670\&amp;tool=pmcentrez\&amp;rendertype=abstract</a>.</p>
<p>Revell, Liam J., D. Luke Mahler, Pedro R. Peres-Neto, and Benjamin D. Redelings. 2011. “a New Phylogenetic Method for Identifying Exceptional Phenotypic Diversification.” <em>Evolution</em> (aug). doi:10.1111/j.1558-5646.2011.01435.x. <a href="http://doi.wiley.com/10.1111/j.1558-5646.2011.01435.x" title="http://doi.wiley.com/10.1111/j.1558-5646.2011.01435.x">http://doi.wiley.com/10.1111/j.1558-5646.2011.01435.x</a>.</p>
<p>Sanderson, M. J., M. J. Donoghue, W. Piel, and T. Eriksson. 1994. “TreeBASE: a prototype database of phylogenetic analyses and an interactive tool for browsing the phylogeny of life.” <em>American Journal of Botany</em> 81: 183.</p>
<p>Schliep, Klaus Peter. 2010. “phangorn: Phylogenetic analysis in R.” <em>Bioinformatics (Oxford, England)</em> 27 (dec): 592–593. doi:10.1093/bioinformatics/btq706. <a href="http://www.ncbi.nlm.nih.gov/pubmed/21169378" title="http://www.ncbi.nlm.nih.gov/pubmed/21169378">http://www.ncbi.nlm.nih.gov/pubmed/21169378</a>.</p>
<p>Schwab, M., N. Karrenbach, and J. Claerbout. 2000. “Making scientific computations reproducible.” <em>Computing in Science &amp; Engineering</em> 2: 61–67. <a href="http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=881708" title="http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=881708">http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=881708</a>.</p>
<p>Stadler, Tanja. 2011a. “Simulating Trees with a Fixed Number of Extant Species.” <em>Systematic biology</em> (apr). doi:10.1093/sysbio/syr029. <a href="http://www.ncbi.nlm.nih.gov/pubmed/21482552" title="http://www.ncbi.nlm.nih.gov/pubmed/21482552">http://www.ncbi.nlm.nih.gov/pubmed/21482552</a>.</p>
<p>———. 2011b. “Mammalian phylogeny reveals recent diversification rate shifts.” <em>Proceedings of the National Academy of Sciences</em> 2011 (mar). doi:10.1073/pnas.1016876108. <a href="http://www.pnas.org/cgi/doi/10.1073/pnas.1016876108" title="http://www.pnas.org/cgi/doi/10.1073/pnas.1016876108">http://www.pnas.org/cgi/doi/10.1073/pnas.1016876108</a>.</p>
<p>Stamatakis, Alexandros. 2006. “RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models.” <em>Bioinformatics (Oxford, England)</em> 22 (nov): 2688–90. doi:10.1093/bioinformatics/btl446. <a href="http://www.ncbi.nlm.nih.gov/pubmed/16928733" title="http://www.ncbi.nlm.nih.gov/pubmed/16928733">http://www.ncbi.nlm.nih.gov/pubmed/16928733</a>.</p>
<p>Sukumaran, Jeet, and Mark T. Holder. 2010. “DendroPy: A Python Library for Phylogenetic Computing.” <em>Bioinformatics</em> 26 (apr): 1569–1571. doi:10.1093/bioinformatics/btq228. <a href="http://www.ncbi.nlm.nih.gov/pubmed/20421198" title="http://www.ncbi.nlm.nih.gov/pubmed/20421198">http://www.ncbi.nlm.nih.gov/pubmed/20421198</a>.</p>
<p>Warren, Dan L., Richard E. Glor, and Michael Turelli. 2008. “Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution.” <em>Evolution</em> 62 (nov): 2868–83. doi:10.1111/j.1558-5646.2008.00482.x. <a href="http://www.ncbi.nlm.nih.gov/pubmed/18752605" title="http://www.ncbi.nlm.nih.gov/pubmed/18752605">http://www.ncbi.nlm.nih.gov/pubmed/18752605</a>.</p>
<p>Webb, Campbell O., David D. Ackerly, and Steven W. Kembel. 2008. “Phylocom: software for the analysis of phylogenetic community structure and trait evolution.” <em>Bioinformatics (Oxford, England)</em> 24 (sep): 2098–100. doi:10.1093/bioinformatics/btn358. <a href="http://www.ncbi.nlm.nih.gov/pubmed/18678590" title="http://www.ncbi.nlm.nih.gov/pubmed/18678590">http://www.ncbi.nlm.nih.gov/pubmed/18678590</a>.</p>
<p>Whitlock, Michael C., Mark a McPeek, Mark D. Rausher, Loren Rieseberg, and Allen J. Moore. 2010. “Data archiving.” <em>The American naturalist</em> 175 (mar): 145–6. doi:10.1086/650340. <a href="http://www.ncbi.nlm.nih.gov/pubmed/20073990" title="http://www.ncbi.nlm.nih.gov/pubmed/20073990">http://www.ncbi.nlm.nih.gov/pubmed/20073990</a>.</p>
<p>Xie, Yihui. 2012. “knitr: A general-purpose package for dynamic report generation in R.” <a href="http://yihui.name/knitr/" title="http://yihui.name/knitr/">http://yihui.name/knitr/</a>.</p>
treebase tutorial
==========

Here are a few introductory examples to illustrate some of the functionality of the package. Thanks in part to new data deposition requirements at journals such as Evolution, Am Nat, and Sys Bio, and
data management plan requirements from NSF, I hope the package will become increasingly useful for teaching by replicating results and for meta-analyses that can be automatically updated as the repository grows. Additional information and bug-reports welcome via the [treebase page](http://ropensci.org/packages/treebase.html#support).

Basic tree and metadata queries
==========

Downloading trees by different queries: by author, taxa, & study. More options are described in the help file.


```r
install.packages("treebase")
```



```r
library(treebase)
```

```
## Loading required package: ape
```



```r
both <- search_treebase("Ronquist or Hulesenbeck", by = c("author", "author"))
dolphins <- search_treebase("\"Delphinus\"", by = "taxon", max_trees = 5)
studies <- search_treebase("2377", by = "id.study")
```



```r
Near <- search_treebase("Near", "author", branch_lengths = TRUE, max_trees = 3)
Near[1]
```

```
## [[1]]
## 
## Phylogenetic tree with 102 tips and 21 internal nodes.
## 
## Tip labels:
## 	Etheostoma_barrenense_A, Etheostoma_rafinesquei_B, Etheostoma_atripinne_A, Etheostoma_atripinne_Y, Etheostoma_atripinne_B, Etheostoma_atripinne_C, ...
## 
## Unrooted; includes branch lengths.
```


We can query the metadata record directly. For instance, plot the growth of Treebase submissions by publication date


```r
all <- download_metadata("", by = "all")
dates <- sapply(all, function(x) as.numeric(x$date))
library(ggplot2)
qplot(dates, main = "Treebase growth", xlab = "Year", binwidth = 0.5)
```


(The previous query could also take a date range)

How do the weekly's do on submissions to Treebase? We construct this in a way that gives us back the indices of the matches, so we can then grab those trees directly. Run the scripts yourself to see if they've changed!


```r
nature <- sapply(all, function(x) length(grep("Nature", x$publisher)) > 0)
science <- sapply(all, function(x) length(grep("^Science$", x$publisher)) > 
    0)
sum(nature)
sum(science)
```


Replicating results
-------------------

A nice paper by Derryberry et al. appeared in Evolution recently on [diversification in ovenbirds and woodcreepers, 0.1111/j.1558-5646.2011.01374.x](http://www.museum.lsu.edu/brumfield/pubs/furnphylogeny2011.pdf). The article mentions that the tree is on Treebase, so let's see if we can replicate their diversification rate analysis: Let's grab the trees by that author and make sure we have the right one:


```r
tree <- search_treebase("Derryberry", "author")[[1]]
plot(tree)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


They fit a variety of diversification rate models avialable in the `laser` R package, which they compare using AIC.


```r
library(laser)
tt <- branching.times(tree)
models <-  list(pb = pureBirth(tt),
                bdfit = bd(tt),
                y2r = yule2rate(tt), # yule model with single shift pt
                ddl = DDL(tt), # linear, diversity-dependent
                ddx = DDX(tt), #exponential diversity-dendent
                sv = fitSPVAR(tt), # vary speciation in time
                ev = fitEXVAR(tt), # vary extinction in time
                bv = fitBOTHVAR(tt)# vary both
                )
names(models[[3]])[5] <- "aic"
aics <- sapply(models, "[[", "aic")
# show the winning model
models[which.min(aics)]
```

```
## $y2r
##         LH         r1         r2        st1        aic 
##  5.051e+02  1.319e-01  3.606e-02  1.871e-01 -1.004e+03
```


Yup, their result agrees with our analysis. Using the extensive toolset for diversification rates in R, we could now rather easily check if these results hold up in newer methods such as TreePar, etc.

Meta-Analysis
-------------

Of course one of the more interesting challenges of having an automated interface is the ability to perform meta-analyses across the set of available phylogenies in treebase. As a simple proof-of-principle, let's check all the phylogenies in treebase to see if they fit a birth-death model or yule model better.

We'll create two simple functions to help with this analysis. While these can be provided by the treebase package, I've included them here to illustrate that the real flexibility comes from being able to create custom functions(These are primarily illustrative; I hope users and developers will create their own. In a proper analysis we would want a few additional checks.)


```r
timetree <- function(tree) {
    check.na <- try(sum(is.na(tree$edge.length)) > 0)
    if (is(check.na, "try-error") | check.na) 
        NULL else try(chronoMPL(multi2di(tree)))
}
drop_errors <- function(tr) {
    tt <- tr[!sapply(tr, is.null)]
    tt <- tt[!sapply(tt, function(x) is(x, "try-error"))]
    print(paste("dropped", length(tr) - length(tt), "trees"))
    tt
}
require(laser)
pick_branching_model <- function(tree) {
    m1 <- try(pureBirth(branching.times(tree)))
    m2 <- try(bd(branching.times(tree)))
    as.logical(try(m2$aic < m1$aic))
}
```


Return only treebase trees that have branch lengths. This has to download every tree in treebase, so this will take a while. Good thing we don't have to do that by hand.


```r
all <- search_treebase("Consensus", "type.tree", branch_lengths = TRUE)
tt <- drop_errors(sapply(all, timetree))
is_yule <- sapply(tt, pick_branching_model)
table(is_yule)
```

% Treebase: An R package for discovery, access and manipulation of online phylogenies 

1.  The TreeBASE portal is an important and rapidly growing repository of
    phylogenetic data. The R statistical environment has also become
    a primary tool for applied phylogenetic analyses across a range
    of questions, from comparative evolution to community ecology to
    conservation planning.

2.  We have developed `treebase`, an open-source software package (freely available from
    [http://cran.r-project.org/web/packages/treebase](http://cran.r-project.org/web/packages/treebase))
    for the R programming environment, providing simplified, *programmatic* and
    interactive access to phylogenetic data in the TreeBASE repository.

3.  We illustrate how this package creates a bridge between the TreeBASE
    repository and the rapidly growing collection of R packages for
    phylogenetics that can reduce barriers to discovery and integration
    across phylogenetic research.

4.  We show how the `treebase` package can be used to facilitate replication
    of previous studies and testing of methods and hypotheses across
    a large sample of phylogenies, which may help make such 
    important reproducibility practices more common.


#### Keywords 

R, software, API, TreeBASE, database, programmatic, workflow 

Introduction
============


Applications that use phylogenetic information as part of their analyses
are becoming increasingly central to both evolutionary and ecological
research. The exponential growth in genetic sequence data available
for all forms of life has driven rapid advances in the methods that
can infer the phylogenetic relationships and divergence times across
different taxa [@Huelsenbeck2001b; @Stamatakis2006; @Drummond2007].
The product of one field has become the raw data of the next.
Unfortunately, while the discipline of bioinformatics has emerged to
help harness and curate the wealth of genetic data with cutting edge
computer science, statistics, and Internet technologies, its counterpart
in evolutionary informatics remains “scattered, poorly documented,
and in formats that impede discovery and integration” [@Parr2011a]. Our
goal in developing the `treebase` package is to provide steps to reduce
these challenges through programmatic and interactive access between the
repositories that store this data and the software tools commonly used
to analyse them.  

The R statistical environment [@RTeam2012] has become a dominant
platform for researchers using phylogenetic data to address a
rapidly expanding set of questions in ecological and evolutionary
processes. These methods include, but are not limited to, ancestral
state reconstruction [@Paradis2004; @Butler2004], diversification
analysis [@Paradis2004; @Rabosky2006b; @Harmon2008], identifying
trait dependent speciation and extinction rates, [@Fitzjohn2010;
@Goldberg2011; @Stadler2011], quantifying the rate and tempo of
trait evolution [@Butler2004; @Harmon2008; @Eastman2011], identifying
evolutionary influences and proxies for community ecology [@Webb2008;
@Kembel2010], connecting phylogeny data to climate patterns [@Warren2008;
@Evans2009b], and simulation of speciation and character evolution
[@Harmon2008; @Stadler2011a; @Boettiger2012], as well as various
manipulations and visualizations of phylogenetic data [@Paradis2004;
@Schliep2010; @Jombart2010; @Revell2011]. A more comprehensive list
of R packages by analysis type is available on the phylogenetics taskview,
[http://cran.r-project.org/web/views/Phylogenetics.html](http://cran.r-project.org/web/views/Phylogenetics.html).
Libraries and packages are developed for use in other general purpose programming environments and languages, including Java [@Maddison2011],
MATLAB [@Blomberg2003] and Python [@Sukumaran2010] and online interfaces
[@Martins2004].  In particular, the Bio::Phylo toolkit [@Vos2011] not
only provides a PERL implementation of some of the common phylogenetic
simulation and visualization tools found in these R libraries, but can
already provide programmatic access to TreeBASE.  Our goal is to bring
similar functionality to the larger suite of applied phylogenetics
methods and user in the R community.

TreeBASE ([http://treebase.org](http://treebase.org)) is an online
repository of phylogenetic data (e.g. trees of species, populations,
or genes) that have been published in a peer-reviewed academic journal,
book, thesis or conference proceedings [@Sanderson1994b; @Morell1996]. The
database can be searched through an online interface which allows users
to find a phylogenetic tree from a particular publication, author or
taxa of interest. TreeBASE provides an application programming interface
(API) that lets computer applications make queries to the database,
known as PhyloWS [@Vos2010].  Our `treebase` package uses this API to
create a direct link between this data and the R environment.  This has
several immediate and important benefits:

1. *Data discovery.*  Users can leverage the rich, higher-level
programming environment provided by the R environment to better identify
data sets appropriate for their research by iteratively constructing
queries for datasets that match appropriate metadata requirements.

2. *Programmatic data access.* Many tasks that are theoretically made
possible by the creation of the Web-base interface to the TreeBASE
repository are not pursued because they would be too laborious for an
exploratory analysis.  The ability to use programmatic access across data
sets to automatically download and perform a reproducible and systematic
analysis using the rich set of tools available in R opens up new avenues
for research.

3. *Automatic updating*.  The TreeBASE repository is expanding rapidly.
The scriptable nature of analyses run with our `treebase` package
means that a study can be rerun on the latest version of the repository
without additional effort but with potential new information.



Programmatic Web Access
-----------------------

The TreeBASE repository makes data accessible via Web queries through a
RESTful (REpresentational State Transfer) interface, which supplies search
conditions in the address URL.  The repository returns the requested
data in XML (extensible markup language) format, conforming to the RSS1.0
standard.  Because the RSS1.0 format allows the search results to also be
viewed in a human-readable format in common browsers such as Safari and
Firefox, the `treebase` package echoes this URL to the console, so that
the user can explore the results in the browser as well.  The `treebase`
package uses the `RCurl` package [@RCurl] to make queries over the Web
to the repository, and the `XML` package  [@XML] to parse the Web page
returned by the repository into meaningful R data objects.  While these
querying and parsing functions comprise most of the code provided in the
`treebase` package, they are hidden from the end user who can interact
with these rich data retrieval and manipulation tools to access data
from these remote repositories in much the same way as data locally
available on the users hard-disk.


While the TreeBASE repository provides phylogenies in both the traditional
Nexus file format and the more data-rich NeXML format [@Vos2012], none
of the R packages currently available for phylogenetic research are
positioned to read these NeXML files.  The next version of the `treebase`
package will provide extraction of metadata information from the NeXML
through XML parsing.



Basic queries
--------------

``` {r libs, echo=FALSE, cache=FALSE}
    library(treebase)
    library(ggplot2)
    theme_set(theme_bw())
````


The `treebase` package allows these queries to be made directly from R.
Programmatic access also allows a user to go beyond the utilities of the
Web interface, constructing more complicated queries and allowing the
user to maintain a record of the commands used to collect their data as
an R script. Scripting the data-gathering process helps reduce errors
and assists in replicating the analysis later, either by the authors or
other researchers [@Peng2011a].


The `search_treebase` function forms the base of the `treebase` package.
Table 1 lists each of the types of queries available through the
`search_treebase` function.  This list can also be found in the function 
documentation through the R command `?search_treebase`.  
Any of the queries available on the Web interface can now be made
directly from R, including downloading and importing a phylogeny into
the R session. For instance, one can search for phylogenies containing
dolphin taxa, "Delphinus," or all phylogenies submitted by a given
author, "Huelsenbeck" using the R commands

``` {r basicQueries, eval=FALSE }
    search_treebase("Delphinus", by="taxon")
    search_treebase("Huelsenbeck", by="author")
````

The `search_treebase` function returns the matching phylogenies as an
R object, ready for analysis.  The package documentation provides many
examples of possible queries.



  Search "by="     PURPOSE
  ---------------  -----------------------------------------------------
  abstract         search terms in the publication abstract
  author           match authors in the publication
  subject          Matches in the subject terms
  doi              The unique object identifier for the publication 
  ncbi             NCBI identifier number for the taxon
  kind.tree        Kind of tree (Gene tree, species tree, barcode tree)  
  type.tree        Type of tree (Consensus or Single)
  ntax             Number of taxa in the matrix
  quality          A quality score for the tree, if it has been rated.  
  study            Match words in the title of the study or publication
  taxon            Taxon scientific name 
  id.study         TreeBASE study ID
  id.tree          TreeBASE's unique tree identifier (Tr.id)
  id.taxon         Taxon identifier number from TreeBase 
  tree             The title for the tree

  Table: Queries available in `search_treebase`.  The first argument is
  the keyword used in the query such as an author's name and the second
  argument indicates the type of query (_i.e._ "author").



Accessing all phylogenies
-------------------------

For certain applications, a user may wish to download all the available
phylogenies from TreeBASE. Using the `cache_treebase` function allows
a user to download a local copy of all trees. Because direct database
dumps are not currently available from treebase.org, this function has
intentional delays to avoid overtaxing the TreeBASE servers, and should
be allowed a full day to run.

``` {r buildcache, eval=FALSE }
treebase <- cache_treebase()
```

Users should still be mindful that these servers are a shared
community resource and not place many queries at once. Users
running large jobs should consider joining the TreeBASE mailing list
([http://sourceforge.net/mailarchive/forum.php?forum_name=treebase-users](http://sourceforge.net/mailarchive/forum.php?forum_name=treebase-users))
to discuss such queries ahead of time.  When query


Once run, the cache is saved
compactly in memory where it can be easily and quickly restored. For
convenience, the `treebase` package comes with a copy already cached,
which can be loaded into memory.


``` {r loadcache}
data(treebase)
````

The cache included in the package will be updated during major package 
revisions.  The timestamp of the cache provided can be viewed in the help
file for the data object using the command `?treebase` (Current cache is May 14, 2012).  All queries from `metadata()` and `search_treebase()` are run against the current online version of the database.  


All of the examples shown in this manuscript are run as shown using 
the `knitr` package for authoring dynamic documents [@knitr], which
helps ensure the results shown are reproducible.  These examples can
be updated by copying and pasting the code shown into the R terminal,
or by recompiling the entire manuscript from the source files found on 
the development Web page for the TreeBASE package,
[github.com/ropensci/treebase](https://github.com/ropensci/treebase). 
Data was accessed to produce the examples shown on `r date()`.  


Data discovery in TreeBASE
==========================

Data discovery involves searching for existing data that meets
certain desired characteristics.  Such searches take advantage of
metadata -- summary information describing the data entries provided in
the repository.  The Web repository uses separate interfaces (APIs)
to access metadata describing the publications associated with the
data entered, such as the publisher, year of publication, etc., and
a different interface to describe the metadata associated with an
individual phylogeny, such as the number of taxa or the kind of tree
(_e.g._ gene tree versus species tree).  The `treebase` package can query
these individual sources of metadata separately, but this information
is most powerful when used in concert -- allowing the construction of
complicated searches that cannot be automated through the Web interface.
The `metadata` function updates a list of all available metadata from
both APIs and returns this information as an R `data.frame`.


```{r metadatatable }
meta <- metadata()
````


From the number of rows of the metadata list we see that there are currently
`r prettyNum(length(unique(meta[["Study.id"]])))` published studies in
the database.  The field columns provided by `metadata` are listed in Table II.  

  metadata field   description
  ---------------  ---------------------------------------------------
  Study.id         TreeBASE study ID
  Tree.id          TreeBASE's unique tree identifier 
  kind             Kind of tree (Gene tree, species tree, barcode tree)  
  type             Type of tree (Consensus or Single)
  quality          A quality score for the tree, if it has been rated.  
  ntaxa            Number of taxa in the matrix
  date             Year the study was published
  author           First author in the publication
  title            The title of the publication 

 Table: Columns of metadata available from the `metadata` function 

Metadata can also be used to reveal trends in the data
deposition which may be useful in identifying patterns or biases in
research or emerging potential types of data.  As a simple example, we
look at trends in the submission patterns of publishers over time: 


``` {r journals}
    date <- meta[["date"]] 
    pub <- meta[["publisher"]]
````

Many journals have only a few submissions, so we will label any not
in the top ten contributing journals as "Other":


``` {r top_journals}
    topten <- sort(table(pub), decreasing=TRUE)[1:10]
    meta[["publisher"]] <- as.character(meta[["publisher"]])
    meta[["publisher"]][!(pub %in% names(topten))] <- "Other"
    meta[["publisher"]] <- as.factor(meta[["publisher"]])
````

We plot the distribution of publication years for phylogenies deposited
in TreeBASE, color coding by publisher in Figure 1.

``` {r dates, fig.width=8, fig.height=3.5, fig.cap="Histogram of publication dates by year, with the code required to generate the figure.", dev.opts=list(pointsize=8) }
library(ggplot2) 
library(reshape2)
df <- acast(meta, date ~ publisher, value.var='publisher', length)
df <- melt(df, varnames=c("date", "publisher"))
ggplot(df) + geom_area(aes(x=date,y=value, fill = publisher)) 
````

Typically we are interested in the metadata describing the phylogenies
themselves rather than just in the publications in which they appeared.
Phylogenetic metadata includes features such as the number of taxa
in the tree, a quality score (if available), kind of tree (gene tree,
species tree, or barcode tree) or whether the phylogeny represents a
consensus tree from a distribution or just a single estimate.

Even simple queries can illustrate the advantage of interacting with
TreeBASE data through an R interface has over the Web interface.  A Web 
interface can only perform the tasks built in by design.  For instance,
rather than performing six separate searches to determine the number 
of consensus vs single phylogenies available for each kind of tree,
we can construct a 2 by 2 table with a single line of code: 


``` {r  kind, eval=FALSE }
table(meta[["kind"]], meta[["type"]])
````

``` {r  kindtable, results='asis', echo=FALSE }
output <-  table(meta[["kind"]], meta[["type"]])
xtable::xtable(output)
````


Reproducible computations
=========================


Reproducible research has become a topic of increasing interest in
recent years, and facilitating access to data and using scripts that
can replicate analyses can help lower barriers to the replication of
statistical and computational results [@Schwab2000; @Gentleman2004;
@Peng2011a].  The `treebase` package facilitates this process, as we
illustrate in a simple example.

Consider the shifts in speciation rate identified by @Derryberry2011
on a phylogeny of ovenbirds and treecreepers. We will seek to not
only replicate the results the authors obtained by fitting the models
provided in the R package  `laser` [@Rabosky2006b], but also compare them
against methods presented in @Stadler2011 and implemented in the package
`TreePar`, which permits speciation models that were not available to
@Derryberry2011 at the time of their study.

Obtaining the tree
------------------


The most expedient way to identify the data uses the digital object identifier
(doi) at the top of most articles, which we use in a call to the 
`search_treebase` function, such as

``` {r doiquery} 
results <- search_treebase("10.1111/j.1558-5646.2011.01374.x", "doi") 
````

The search returns a list, since some publications can contain many trees.
In this case our phylogeny is in the only element of the list. 

Having imported the phylogenetic tree corresponding to this
study, we can quickly replicate their analysis of which diversification
process best fits the data.  These steps can be easily implemented using
the phylogenetics packages we have just mentioned. 


``` {r RRexamplelibs, echo=FALSE  }
library(ape)
library(laser)
````

For instance, we can calculate the branching times of each node on the phylogeny,

``` {r RRexamp }
bt <- branching.times(results[[1]])
````

and then begin to fit each model the authors have tested, such as the pure birth model,

``` {r }
yule <- pureBirth(bt)
````

or the birth-death model, 

``` {r }
birth_death <- bd(bt)
````

The estimated models are now available in the active R session where
we can further explore them as we go along.  The appendix shows the
estimation and comparison of all the models originally considered by
@Derryberry2011.


In this fast-moving field, new methods often become available between the
time of submission and time of publication of a manuscript. For instance,
the more sophisticated models introduced in @Stadler2011 were not used
in the original study, but have since been made available in the recent
package, `TreePar`.  These richer models permit a shift in the speciation
or extinction rate to occur multiple times throughout the course of
the phylogeny.

``` {r   TreeSimError, echo=FALSE  }
# Locale settings to be safe
Sys.setlocale(locale="C") -> locale
rm(list="locale") # clean up
````

We load the new method and format the phylogeny appropriately using the
R commands:

``` {r   treepar  }
library(TreePar)
x <- sort(getx(results[[1]]), decreasing = TRUE)
````


``` {r   treepar_yule_demo  }
yule_models <- bd.shifts.optim(x, sampling = c(1,1,1,1), 
  grid = 5, start = 0, end = 60, yule = TRUE)[[2]]
````


As a comparison of
speciation models is not the focus of this paper, the complete code
and explanation for these steps are provided in an appendix.  Happily,
this analysis confirms the original author's conclusions, even when the
more general models of @Stadler2011 are considered.




Analyses across many phylogenies
================================

Large scale comparative analyses that seek to characterize evolutionary
patterns across many phylogenies are increasingly common in phylogenetic
methods research [_e.g._ @McPeek2007; @Phillimore2008; @McPeek2008;
@Quental2010; @Davies2011a].  Sometimes referred to by their authors as
meta-analyses, these approaches have focused on re-analyzing phylogenetic
trees collected from many different earlier publications.  This is a
more direct approach than the traditional concept of meta-analysis where
statistical results from earlier studies are weighted by their sample size
without being able to access the raw data.  Because the identical analysis
can be repeated on the original data from each study, this approach avoids
some of the statistical challenges inherent in traditional meta-analyses
summarizing results across heterogeneous approaches.

To date, researchers have gone through heroic efforts simply to assemble
these data sets from the literature.  As described in @McPeek2007; (emphasis added) 

> One data set was based on 163 published species-level molecular
phylogenies of arthropods, chordates, and mollusks.  A PDF format file
of each article was obtained, and a digital snapshot of the figure was
taken in Adobe Acrobat  7.0. This image was transferred to a PowerPoint
(Microsoft) file and printed on a laser printer. The phylogenies included
in this study are listed in the appendix. _All branch lengths were
measured by hand from these  printed sheets using dial calipers._

Despite the recent emergence of digital tools that could now facilitate
this analysis without mechanical calipers, [_e.g._ treesnatcher, @Laubach2007],
it is  easier and less error-prone to pull properly
formatted phylogenies from the database for this purpose. Moreover, as the
available data increases with subsequent publications, updating earlier
meta-analyses can become increasingly tedious. Using `treebase`, a user
can apply any analysis they have written for a single phylogeny across
the entire collection of suitable phylogenies in TreeBASE, which can help
overcome such barriers to discovery and integration at this large scale.
Using the functions we introduced above, we provide a simple example that
computes the gamma statistic of @Pybus2000, which provides a measure
of when speciation patterns differ from the popular birth-death model.

Tests across many phylogenies
-----------------------------

A standard test of the constant rate of diversification is the gamma
statistic of @Pybus2000 which tests the null hypothesis that the rates
of speciation and extinction are constant. Under the null hypothesis,
The gamma statistic is normally distributed about 0; values larger than
0 indicate that internal nodes are closer to the tip then expected, while
values smaller than 0 indicate nodes farther from the tip then expected.
In this section, we collect all phylogenetic trees from TreeBASE and
select those with branch length data that we can time-calibrate using
tools available in R.  We can then calculate the distribution of this
statistic for all available trees, and compare these results with those
from the analyses mentioned above.  As we will use all trees in the 
repository, we use the cached copy of TreeBASE phylogenies described 
above to reduce load on TreeBASE servers.  

We will only be able to use those phylogenies that include branch length data,
which we can determine from the `have_branchlength` function in the `treebase`
package. We drop those that do not from the data set, 

``` {r branchlengthfilter}
      have <- have_branchlength(treebase)
      branchlengths <- treebase[have]

````

Like most comparative methods, this analysis will require ultrametric trees
(branch lengths proportional to time, rather than to the nucleotide substitution rate). 
As most of these phylogenies are calibrated with branch length proportional
to mutational step, we must time-calibrate each of them first. The following function
drops trees which cannot meet the assumptions of the time-calibration function.  

``` {r timetree, cache=TRUE}
timetree <- function(tree)
    try(chronoMPL(multi2di(tree)), silent=TRUE) 
tt <- drop_nontrees(sapply(branchlengths, timetree))
````

At this point we have `r length(tt)` time-calibrated phylogenies over which
we will apply the diversification rate analysis. 
Computing the gamma test statistic to identify deviations from the 
constant-rates model takes a single line,

``` {r }
gammas <- sapply(tt,  gammaStat)
````

and the resulting distribution of the statistic across available
trees is shown Fig 2.  While researchers have often considered this
statistic for individual phylogenies, we are unaware of any study that has
visualized the empirical distribution of this statistic across thousands
of phylogenies.  The overall distribution appears slightly
skewed towards positive values.  This could indicate increasing rate of 
speciation or constant extinction rates.  While differences in sampling 
may account for much of the spread observed, the position and identity of 
outlier phylogenies could suggest new hypotheses and potential directions for
further exploration.

``` {r gammadist, fig.width=6.2, fig.cap="Distribution of the gamma statistic across phylogenies in TreeBASE. Strongly positive values are indicative of an increasing rate of evolution (excess of nodes near the tips), very negative values indicate an early burst of diversification (an excess of nodes near the root)."}
qplot(gammas)+xlab("gamma statistic")
````


Conclusion
==========

While we have focused on examples that require no additional data
beyond the phylogeny, a wide array of methods combine this data with
information about the traits, geography, or ecological community of the
taxa represented. In such cases, we would need programmatic access to
the trait data as well as the phylogeny. The Dryad digital repository
([http://datadryad.org](http://datadryad.org)) is an effort in this 
direction. While programmatic access to the repository is possible
through the `rdryad` package [@Chamberlain2012], variation in data
formatting must first be overcome before similar direct access to 
the data is possible.  Dedicated databases such as FishBASE
([http://fishbase.org](http://fishbase.org)) may be another alternative,
where morphological data can be queried for a list of species using the
`rfishbase` package [@rfishbase]. The development of similar software
for programmatic data access will rapidly extend the space and scale of
possible analyses.

The recent advent of mandatory data archiving in many of the major
journals publishing phylognetics-based research [_e.g._ @Fairbairn2010;
@Piwowar2011; @Whitlock2010], is a particularly promising development that
should continue to fuel the trend of submissions seen in Fig. 1. Accompanied
by faster and more inexpensive techniques of NextGen sequencing, and the
rapid expansion in phylogenetic applications, we anticipate this rapid
growth in available phylogenies will continue. Faced with this flood
of data, programmatic access becomes not only increasingly powerful but
an increasingly necessary way to ensure we can still see the forest for
all the trees.

Acknowledgements
================

CB wishes to thank S. Price for feedback on the manuscript, the
TreeBASE developer team for building and supporting the repository,
and all contributers to TreeBASE. CB is supported by a Computational
Sciences Graduate Fellowship from the Department of Energy under grant
number DE-FG02-97ER25308. The `treebase` package is is part of the 
rOpenSci project ([ropensci.org](http://ropensci.org)).


 
``` {r save, echo=FALSE }
    save(list=ls(), file="knit_treebase.rda")
````




# References


``` {r loadcache, echo=FALSE, include=FALSE}
require(treebase)
data(treebase)
````


Appendix
========

Reproducible computation: A diversification rate analysis
---------------------------------------------------------

This appendix illustrates the diversification rate analysis discussed in
the main text.  For completeness we begin by executing the code discussed
in the manuscript which locates, downloads, and imports the relevant data:

``` {r libs, include=FALSE, dependson="loadcache"}
library(treebase)
results <- search_treebase("10.1111/j.1558-5646.2011.01374.x", "doi") 
bt <- branching.times(results[[1]])
````


Different diversification models make different assumptions about the rate
of speciation, extinction, and how these rates may be changing over time.
The original authors consider eight different models, implemented in the
laser package [@Rabosky2006b]. This code fits each of the eight models
to that data:

``` {r   RR, tidy=FALSE, dependson="loadcache"  }
library(ape)
library(laser)
models <- list(
  yule = pureBirth(bt),  
  birth_death = bd(bt),     
  yule.2.rate = yule2rate(bt),
  linear.diversity.dependent = DDL(bt),    
  exponential.diversity.dependent = DDX(bt),
  varying.speciation_rate = fitSPVAR(bt),  
  varying.extinction_rate = fitEXVAR(bt),  
  varying_both = fitBOTHVAR(bt))
````

Each of the model estimate includes an Akaike Information Criterion
(AIC) score indicating the goodness of fit, penalized by model complexity
(lower scores indicate better fits) We ask R to tell us which model has
the lowest AIC score,

``` {r   extract_aic, dependson="loadcache" }
aics <- sapply(models, function(model) model$aic)
best_fit <- names(models[which.min(aics)])
```` 

and confirm the result presented in @Derryberry2011; that the best-fit
model in the laser analysis was a Yule (net diversification rate) model
with two separate rates.  

We can ask ` TreePar ` to see if a model with more rate shifts is favoured
over this single shift, a question that was not possible to address using
the tools provided in `laser`. The previous analysis also considers a
birth-death model that allowed speciation and extinction rates to be
estimated separately, but did not allow for a shift in the rate of such
a model.  In the main text we introduced a model from @Stadler2011 that
permitted up to 3 change-points in the speciation rate of the Yule model,

``` {r   treepar-yule, include=FALSE}
yule_models <- bd.shifts.optim(x, sampling = c(1,1,1,1), 
  grid = 5, start = 0, end = 60, yule = TRUE)[[2]]
````


``` {r   treepar-yule_show, eval=FALSE  }
yule_models <- bd.shifts.optim(x, sampling = c(1,1,1,1), 
  grid = 5, start = 0, end = 60, yule = TRUE)[[2]]
````

We can also compare the performance of models which allow up to three
shifts while estimating extinction and speciation rates separately:

``` {r   treepar-birthdeath, include=FALSE  }
birth_death_models <- bd.shifts.optim(x, sampling = c(1,1,1,1), 
  grid = 5, start = 0, end = 60, yule = FALSE)[[2]]
````

``` {r   treepar-birthdeath_show, eval=FALSE  }
birth_death_models <- bd.shifts.optim(x, sampling = c(1,1,1,1), 
  grid = 5, start = 0, end = 60, yule = FALSE)[[2]]
````

The models output by these functions are ordered by increasing number
of shifts.  We can select the best-fitting model by AIC score, which is
slightly cumbersome in `TreePar` syntax.  First, we compute the AIC scores
of both the `yule_models` and the `birth_death_models` we fitted above,

``` {r   aic-yule1  }
yule_aic <- 
sapply(yule_models, function(pars)
                    2 * (length(pars) - 1) + 2 * pars[1] )
birth_death_aic <- 
sapply(birth_death_models, function(pars)
                            2 * (length(pars) - 1) + 2 * pars[1] )
````

Then we generate a list identifying which model has the best (lowest)
AIC score among the Yule models and which has the best AIC score among
the birth-death models,

```{r aic-yule2}
best_no_of_rates <- list(Yule = which.min(yule_aic), 
                         birth.death = which.min(birth_death_aic))
````

The best model is then whichever of these has the smaller AIC value.  

```{r aic-yule3}
best_model <- which.min(c(min(yule_aic), min(birth_death_aic)))
````


which still confirms that the `r names(best_no_of_rates[best_model])` `r best_no_of_rates[best_model][[1]]`-rate  model is still the best choice based on AIC score.  Of the eight models 
in this second analysis, only three were in the original set considered 
(Yule 1-rate and 2-rate, and birth-death without a shift), so we could by
no means have been sure ahead of time that a birth death with a shift, or
a Yule model with a greater number of shifts, would not have fitted better.  


# References

---
title: "Treebase Tutorial"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Treebase Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Here are a few introductory examples to illustrate some of the functionality of the package. Thanks in part to new data deposition requirements at journals such as Evolution, Am Nat, and Sys Bio, and
data management plan requirements from NSF, I hope the package will become increasingly useful for teaching by replicating results and for meta-analyses that can be automatically updated as the repository grows. Additional information and bug-reports welcome via the [treebase page](http://github.com/ropensci/treebase/issues).

Basic tree and metadata queries
==========

Downloading trees by different queries: by author, taxa, & study. More options are described in the help file.

```{r include=FALSE}
knitr::opts_chunk$set(eval = FALSE)
```

```{r install, eval=FALSE}
install.packages('treebase')
```

```{r loadpkg}
library(treebase)
```

```{r search1, eval=FALSE, message=FALSE, warning=FALSE}
both <- search_treebase("Ronquist or Hulesenbeck", by=c("author", "author"))
dolphins <- search_treebase('"Delphinus"', by="taxon", max_trees=5)
studies <- search_treebase("2377", by="id.study")
```

```{r search2, message=FALSE, warning=FALSE}
Near <- search_treebase("Near", "author", branch_lengths=TRUE, max_trees=3)
Near[1]
```

We can query the metadata record directly. For instance, plot the growth of Treebase submissions by publication date

```{r metadata, eval=FALSE, message=FALSE, warning=FALSE}
all <- download_metadata("", by="all")
dates <- sapply(all, function(x) as.numeric(x$date))
library(ggplot2)
qplot(dates, main="Treebase growth", xlab="Year", binwidth=.5)
```

(The previous query could also take a date range)

How do the weekly's do on submissions to Treebase? We construct this in a way that gives us back the indices of the matches, so we can then grab those trees directly. Run the scripts yourself to see if they've changed!

```{r eval=FALSE, message=FALSE, warning=FALSE}
nature <- sapply(all, function(x) length(grep("Nature", x$publisher))>0)
science <- sapply(all, function(x) length(grep("^Science$", x$publisher))>0)
sum(nature)
sum(science)
```

Replicating results
-------------------

A nice paper by Derryberry et al. appeared in Evolution recently on diversification in ovenbirds and woodcreepers, ([doi:10.1111/j.1558-5646.2011.01374.x](http://dx.doi.org/10.1111/j.1558-5646.2011.01374.x)). The article mentions that the tree is on Treebase, so let's see if we can replicate their diversification rate analysis: Let's grab the trees by that author and make sure we have the right one:

```{r message=FALSE, warning=FALSE}
search_treebase("Derryberry", "author")[[1]] -> tree
plot(tree)
```

They fit a variety of diversification rate models avialable in the `laser` R package, which they compare using AIC.

```{r message=FALSE, warning=FALSE, eval = FALSE}
library(laser)
tt <- branching.times(tree)
models <-  list(pb = pureBirth(tt),
                bdfit = bd(tt),
                y2r = yule2rate(tt), # yule model with single shift pt
                ddl = DDL(tt), # linear, diversity-dependent
                ddx = DDX(tt), #exponential diversity-dendent
                sv = fitSPVAR(tt), # vary speciation in time
                ev = fitEXVAR(tt), # vary extinction in time
                bv = fitBOTHVAR(tt)# vary both
                )
names(models[[3]])[5] <- "aic"
aics <- sapply(models, "[[", "aic")
# show the winning model
models[which.min(aics)]


```

Their result agrees with our analysis. Using the extensive toolset for diversification rates in R, we could now rather easily check if these results hold up in newer methods such as `TreePar`, etc.

Meta-Analysis
-------------

Of course one of the more interesting challenges of having an automated interface is the ability to perform meta-analyses across the set of available phylogenies in treebase. As a simple proof-of-principle, let's check all the phylogenies in treebase to see if they fit a birth-death model or yule model better.

We'll create two simple functions to help with this analysis. While these can be provided by the treebase package, I've included them here to illustrate that the real flexibility comes from being able to create custom functions(These are primarily illustrative; I hope users and developers will create their own. In a proper analysis we would want a few additional checks.)

```{r eval=FALSE, message=FALSE, warning=FALSE}
timetree <- function(tree){
    check.na <- try(sum(is.na(tree$edge.length))>0)
    if(is(check.na, "try-error") | check.na)
      NULL
    else
    try( chronoMPL(multi2di(tree)) )
}
drop_errors <- function(tr){
  tt <- tr[!sapply(tr, is.null)]
  tt <- tt[!sapply(tt, function(x) is(x, "try-error"))]
  print(paste("dropped", length(tr)-length(tt), "trees"))
  tt
}
require(laser)
pick_branching_model <- function(tree){
  m1 <- try(pureBirth(branching.times(tree)))
  m2 <- try(bd(branching.times(tree)))
  as.logical(try(m2$aic < m1$aic))
}
```

Return only treebase trees that have branch lengths. This has to download every tree in treebase, so this will take a while. Good thing we don't have to do that by hand.

```{r eval=FALSE, message=FALSE, warning=FALSE}
all <- search_treebase("Consensus", "type.tree", branch_lengths=TRUE)
tt <- drop_errors(sapply(all, timetree))
is_yule <- sapply(tt, pick_branching_model)
table(is_yule)
```% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_treebase.R
\name{cache_treebase}
\alias{cache_treebase}
\title{A function to cache the phylogenies in treebase locally}
\usage{
cache_treebase(file = paste("treebase-", Sys.Date(), ".rda", sep = ""),
  pause1 = 3, pause2 = 3, attempts = 10, max_trees = Inf,
  only_metadata = FALSE, save = TRUE)
}
\arguments{
\item{file}{filename for the cache, otherwise created with datestamp}

\item{pause1}{number of seconds to hesitate between requests}

\item{pause2}{number of seconds to hesitate between individual files}

\item{attempts}{number of attempts to access a particular resource}

\item{max_trees}{maximum number of trees to return (default is Inf)}

\item{only_metadata}{option to only return metadata about matching trees}

\item{save}{logical indicating whether to save a file with the resuls.}
}
\value{
saves a cached file of treebase
}
\description{
A function to cache the phylogenies in treebase locally
}
\details{
it's a good idea to let this run overnight
}
\examples{
\dontrun{
 treebase <- cache_treebase()
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treebase.R
\name{get_nex}
\alias{get_nex}
\title{imports phylogenetic trees from treebase. internal function}
\usage{
get_nex(query, max_trees = "last()", returns = "tree",
  curl = getCurlHandle(), verbose = TRUE, pause1 = 1, pause2 = 1,
  attempts = 5, only_metadata = FALSE)
}
\arguments{
\item{query}{: a phylows formatted search,
see https://sourceforge.net/apps/mediawiki/treebase/index.php?title=API}

\item{max_trees}{limits the number of trees returned
should be kept.}

\item{returns}{should return the tree object or the matrix (of sequences)}

\item{curl}{the handle to the curl}

\item{verbose}{a logical indicating if output should be printed to screen}

\item{pause1}{number of seconds to hesitate between requests}

\item{pause2}{number of seconds to hesitate between individual files}

\item{attempts}{number of attempts to access a particular resource}
}
\value{
A list object containing all the trees matching the search
   (of class phylo)
}
\description{
imports phylogenetic trees from treebase. internal function
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_metadata.R
\name{show_metadata}
\alias{show_metadata}
\title{Get the metadata associated with the study in which the phylogeny
 was published.}
\usage{
show_metadata(study.id, curl = getCurlHandle())
}
\arguments{
\item{study.id}{The treebase study id (numbers only, specify in quotes)}

\item{curl}{if calling in series many times, call getCurlHandle()
 first and then pass the return value in here.  avoids repeated
handshakes with server.}
}
\description{
Get the metadata associated with the study in which the phylogeny
 was published.
}
\details{
if the tree is imported with search_treebase, 
then this is in tree$S.id
}
\keyword{internal}
\keyword{utilities}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{treebase}
\alias{treebase}
\title{treebase.rda}
\description{
Contains a cache of all phylogenies \code{cache_treebase()} function was able
to pull down when run  on 2012-05-14.
}
\details{
recreate with: \code{
 cache_treebase() 
}
}
\keyword{data}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R, R/metadata.R
\docType{data}
\name{metadata}
\alias{metadata}
\title{metadata.rda}
\usage{
metadata(phylo.md = NULL, oai.md = NULL)
}
\arguments{
\item{phylo.md}{cached phyloWS (tree) metadata, (optional)}

\item{oai.md}{cached OAI-PMH (study) metadata (optional)}
}
\value{
a data frame of all available metadata, (as a data.table object)
columns are: "Study.id", "Tree.id", "kind", "type", "quality", "ntaxa"    
"date", "publisher", "author", "title".
}
\description{
Contains a cache of all publication metadata the search_metadata()
to pull down when run  on 2012-05-12.

generate a table of all available metadata for TreeBASE entries
}
\details{
recreate with: \code{
 search_metadata() 
}
}
\examples{
\dontrun{
meta <- metadata()
meta[publisher \%in\% c("Nature", "Science") & ntaxa > 50 & kind == "Species Tree",]
}
}
\keyword{data}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_metadata.R
\name{get_study_id}
\alias{get_study_id}
\title{return the study.id from the search results.}
\usage{
get_study_id(search_results)
}
\arguments{
\item{search_results}{the output of download_metadata, or a subset thereof}
}
\value{
the study id
}
\description{
get_study_id is deprecated, and now can be performed more easily using
phylo_metadata and oai_metadata search functions.
}
\details{
this function is commonly used to get trees corresponding
  to the metadata search.
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_metadata.R
\name{metadata_from_oai}
\alias{metadata_from_oai}
\title{Internal function for OAI-MPH interface to the Dryad database}
\usage{
metadata_from_oai(query, curl = curl)
}
\arguments{
\item{query}{a properly formed url query to dryad}

\item{curl}{if calling in series many times, call getCurlHandle() first and 
then pass the return value in here. Avoids repeated handshakes with server.}
}
\description{
Internal function for OAI-MPH interface to the Dryad database
}
\seealso{
\code{\link{dryad_metadata}}
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_metadata.R
\name{clean_data}
\alias{clean_data}
\title{clean the fish.base data into pure ASCII}
\usage{
clean_data(metadata)
}
\arguments{
\item{a}{list item with fishbase data}
}
\value{
the item scrubbed of non-ASCII characters
}
\description{
clean the fish.base data into pure ASCII
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{oai_metadata}
\alias{oai_metadata}
\title{Search the OAI-PMH metadata by date, publisher, or identifier}
\usage{
oai_metadata(x = c("date", "publisher", "author", "title", "Study.id",
  "attributes"), metadata = NULL, ...)
}
\arguments{
\item{x}{one of "date", "publisher", "identifier" for the study}

\item{metadata}{returned from \code{download_metadata} function. 
if not specified will download latest copy from treebase.  Pass
in the value during repeated calls to speed function runtime substantially}

\item{...}{additional arguments to \code{download_metadata}}
}
\value{
a list of values matching the query
}
\description{
Search the OAI-PMH metadata by date, publisher, or identifier
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_metadata.R
\name{drop_nonascii}
\alias{drop_nonascii}
\title{remove non-ASCII characters}
\usage{
drop_nonascii(string)
}
\arguments{
\item{string}{any character string}
}
\value{
the string after dropping all html tags to spaces
}
\description{
remove non-ASCII characters
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{phylo_metadata}
\alias{phylo_metadata}
\title{Search the PhyloWS metadata}
\usage{
phylo_metadata(x = c("Study.id", "Tree.id", "kind", "type", "quality",
  "ntaxa"), metadata = NULL, ...)
}
\arguments{
\item{x}{one of "Study.ids", "Tree.ids", "kind", "type", "quality", "ntaxa"}

\item{metadata}{returned from \code{search_treebase} function. 
if not specified will download latest copy of PhyloWS metadata from treebase.  Pass
in search results value during repeated calls to speed function runtime substantially}

\item{...}{additional arguments to \code{search_treebase}}
}
\value{
a list of the values matching the query
}
\description{
Search the PhyloWS metadata
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_metadata.R
\name{download_metadata}
\alias{download_metadata}
\title{Download the metadata on treebase using the OAI-MPH interface}
\usage{
download_metadata(query = "", by = c("all", "until", "from"),
  curl = getCurlHandle())
}
\arguments{
\item{query}{a date in format yyyy-mm-dd}

\item{by}{return all data "until" that date, 
"from" that date to current, or "all"}

\item{curl}{if calling in series many times, call getCurlHandle() first and 
then pass the return value in here. Avoids repeated handshakes with server.}
}
\description{
Download the metadata on treebase using the OAI-MPH interface
}
\details{
query must be#'  download_metadata(2010-01-01, by="until")
 all isn't a real query type, but will return all trees regardless of date
}
\examples{
\dontrun{
Near <- search_treebase("Near", "author", max_trees=1)
 metadata(Near[[1]]$S.id)
## or manualy give a sudy id
metadata("2377")

### get all trees from a certain depostition date forwards ##
m <- download_metadata("2009-01-01", by="until")
## extract any metadata, e.g. publication date:
dates <- sapply(m, function(x) as.numeric(x$date))
hist(dates, main="TreeBase growth", xlab="Year")

### show authors with most tree submissions in that date range 
authors <- sapply(m, function(x){
   index <- grep( "creator", names(x))
     x[index] 
})
a <- as.factor(unlist(authors))
head(summary(a))

## Show growth of TreeBASE 
all <- download_metadata("", by="all")
dates <- sapply(all, function(x) as.numeric(x$date))
hist(dates, main="TreeBase growth", xlab="Year")

## make a barplot submission volume by journals
journals <- sapply(all, function(x) x$publisher)
J <- tail(sort(table(as.factor(unlist(journals)))),5)
b<- barplot(as.numeric(J))
text(b, names(J), srt=70, pos=4, xpd=T)
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treebase.R
\name{search_treebase}
\alias{search_treebase}
\title{A function to pull in the phyologeny/phylogenies matching a search query}
\usage{
search_treebase(input, by, returns = c("tree", "matrix"),
  exact_match = FALSE, max_trees = Inf, branch_lengths = FALSE,
  curl = getCurlHandle(), verbose = TRUE, pause1 = 0, pause2 = 0,
  attempts = 3, only_metadata = FALSE)
}
\arguments{
\item{input}{a search query (character string)}

\item{by}{the kind of search; author, taxon, subject, study, etc
(see list of possible search terms, details)}

\item{returns}{should the fn return the tree or the character matrix?}

\item{exact_match}{force exact matching for author name, taxon, etc.
Otherwise does partial matching}

\item{max_trees}{Upper bound for the number of trees returned, good for
keeping possibly large initial queries fast}

\item{branch_lengths}{logical indicating whether should only return
trees that have branch lengths.}

\item{curl}{the handle to the curl web utility for repeated calls, see
the getCurlHandle() function in RCurl package for details.}

\item{verbose}{logical indicating level of progress reporting}

\item{pause1}{number of seconds to hesitate between requests}

\item{pause2}{number of seconds to hesitate between individual files}

\item{attempts}{number of attempts to access a particular resource}

\item{only_metadata}{option to only return metadata about matching trees
which lists study.id, tree.id, kind (gene,species,barcode) type (single, consensus)
number of taxa, and possible quality score.}
}
\value{
either a list of trees (multiphylo) or a list of character matrices
}
\description{
A function to pull in the phyologeny/phylogenies matching a search query
}
\details{
Choose the search type.  Options are: \itemize{
\item{abstract        }{ search terms in the publication abstract}
\item{author          }{ match authors in the publication}
\item{subject         }{ match subject}
\item{doi             }{ the unique object identifier for the publication }
\item{ncbi            }{ NCBI identifier number for the taxon}
\item{kind.tree       }{ Kind of tree (Gene tree, species tree, barcode tree)  }
\item{type.tree       }{ type of tree (Consensus or Single)}
\item{ntax            }{ number of taxa in the matrix}
\item{quality         }{ A quality score for the tree, if it has been rated.  }
\item{study           }{ match words in the title of the study or publication}
\item{taxon           }{ taxon scientific name }
\item{id.study        }{ TreeBASE study ID}
\item{id.tree         }{ TreeBASE's unique tree identifier (Tr.id)}
\item{id.taxon        }{ taxon identifier number from TreeBase }
\item{tree            }{ The title for the tree}
\item{type.matrix     }{ Type of matrix }
\item{matrix          }{ Name given the the matrix }
\item{id.matrix       }{ TreeBASE's unique matrix identifier}
\item{nchar           }{ number of characters in the matrix}
}

The package provides partial support for character matrices provided by TreeBASE.
At the time of writing, TreeBASE permits ambiguous DNA characters in these matrices,
such as `{CG}` indicating either a C or G, which is not supported by any R interpreter,
and thus may lead to errors.
  for a description of all possible search options, see
  https://spreadsheets.google.com/pub?key=rL--O7pyhR8FcnnG5-ofAlw.
}
\examples{
\dontrun{
## defaults to return phylogeny
Huelsenbeck <- search_treebase("Huelsenbeck", by="author")

## can ask for character matrices:
wingless <- search_treebase("2907", by="id.matrix", returns="matrix")

## Some nexus matrices don't meet read.nexus.data's strict requirements,
## these aren't returned
H_matrices <- search_treebase("Huelsenbeck", by="author", returns="matrix")

## Use Booleans in search: and, or, not
## Note that by must identify each entry type if a Boolean is given
HR_trees <- search_treebase("Ronquist or Hulesenbeck", by=c("author", "author"))

## We'll often use max_trees in the example so that they run quickly,
## notice the quotes for species.
dolphins <- search_treebase('"Delphinus"', by="taxon", max_trees=5)
## can do exact matches
humans <- search_treebase('"Homo sapiens"', by="taxon", exact_match=TRUE, max_trees=10)
## all trees with 5 taxa
five <- search_treebase(5, by="ntax", max_trees = 10)
## These are different, a tree id isn't a Study id.  we report both
studies <- search_treebase("2377", by="id.study")
tree <- search_treebase("2377", by="id.tree")
c("TreeID" = tree$Tr.id, "StudyID" = tree$S.id)
## Only results with branch lengths
## Has to grab all the trees first, then toss out ones without branch_lengths
Near <- search_treebase("Near", "author", branch_lengths=TRUE)
 }
}
\keyword{utility}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_metadata.R
\name{get_study}
\alias{get_study}
\title{return the trees in treebase that correspond to the search results
get_study is deprecated, and now can be performed more easily using
phylo_metadata and oai_metadata search functions.}
\usage{
get_study(search_results, curl = getCurlHandle(), ...)
}
\arguments{
\item{search_results}{the output of download_metadata, or a subset thereof}

\item{curl}{the handle to the curl web utility for repeated calls, see
the getCurlHandle() function in RCurl package for details.}

\item{...}{additional arguments to pass to search_treebase}
}
\value{
all corresponding phylogenies.
}
\description{
return the trees in treebase that correspond to the search results
get_study is deprecated, and now can be performed more easily using
phylo_metadata and oai_metadata search functions.
}
\details{
this function is commonly used to get trees corresponding
  to the metadata search.
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treebase.R
\name{drop_nontrees}
\alias{drop_nontrees}
\title{drop errors from the search}
\usage{
drop_nontrees(tr)
}
\arguments{
\item{tr}{a list of phylogenetic trees returned by search_treebase}
}
\value{
the list of phylogenetic trees returned successfully
}
\description{
drop errors from the search
}
\details{
primarily for the internal use of search_treebase, but may be useful
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treebase.R
\name{have_branchlength}
\alias{have_branchlength}
\title{Simple function to identify which trees have branch lengths}
\usage{
have_branchlength(trees)
}
\arguments{
\item{trees}{a list of phylogenetic trees (ape/phylo format)}
}
\value{
logical string indicating which have branch length data
}
\description{
Simple function to identify which trees have branch lengths
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_metadata.R
\name{dryad_metadata}
\alias{dryad_metadata}
\title{Search the dryad metadata archive}
\usage{
dryad_metadata(study.id, curl = getCurlHandle())
}
\arguments{
\item{study.id}{the dryad identifier}

\item{curl}{if calling in series many times, call getCurlHandle() first and 
then pass the return value in here. Avoids repeated handshakes with server.}
}
\value{
a list object containing the study metadata
}
\description{
Search the dryad metadata archive
}
\examples{
\dontrun{
  dryad_metadata("10255/dryad.12")
}
}

