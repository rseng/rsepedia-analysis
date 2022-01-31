<!--
This is an issue template. Please fill in the relevant details in the
sections below.
-->

##### Orange version
<!-- From menu _Help→About→Version_ or code `Orange.version.full_version` -->


##### Orange-Bioinformatics version
<!-- From window _Options→Add ons..._ find Orange-Bioinformatics -->


##### Expected behavior



##### Actual behavior



##### Steps to reproduce the behavior



##### Additional info (worksheets, data, screenshots, ...)

Quality Control
===============

![Quality Control widget icon](icons/quality-control.png)

Computes and plots distances between experiments or replicates.

Signals
-------

**Inputs**:

- **Data**

  Data set.

**Outputs**:

- (None)

Description
-----------

**Quality Control** measures distances between experiments (usually replicates) for a selected label. The widget visualizes
distances by selected label. Experiments that lie the farthest from the initial black line
should be inspected for anomalies.

![image](images/QualityControl-stamped.png)

1. Information on the input.
2. Separate experiments by label.
3. Sort experiments by label.
4. Compute distances by:
   - [**Pearson correlation**](https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient)
   - [**Euclidean distances**](https://en.wikipedia.org/wiki/Euclidean_distance)
   - [**Spearman correlation**](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
5. Hover over the vertical line to display the information on a chosen instance. Click on the black line to
   change the reference to that instance.

Example
-------

**Quality Control** widget gives us feedback on the quality of our data and can be connected to any widget with
data output. In the example above (see the image under *Description*) we fed 9 experiments of 
*Cyclic AMP pulsing* of *Dictyostelium discoideum* from **GenExpress** widget into **Quality Control**
and separated them by timepoint label. We found replicate 2 from tp 2 among the tp 1 data, meaning
we should inspect these data further.
Databases Update
================

![Databases widget icon](icons/databases-update.png)

Updates local systems biology databases, like gene ontologies,
annotations, gene names, protein interaction networks, and similar.

Signals
-------

**Inputs**:

- None

**Outputs**:

- None

Description
-----------

With the bioinformatics add-on you can access several databases
directly from Orange. The widget can also be used to update and
manage locally stored databases.

![Databases widget](images/DatabasesUpdate-stamped.png)

1. Find the desired database.
2. A list of available databases described with data source, update availability, date of your last update and file size. 
   A large **Update** button will be displayed next to the database that needs to be updated.
3. **Update All** will update and **Download All** will download all the selected databases.
   **Cancel** will abort the action.
4. Some data sets require the *Access code*. Type it in the provided field to access the database.
5. Information on the selected databases.

To get a more detailed information on the particular database hover on its name.

![Databases widget](images/DatabasesUpdate2.png)
GenExpress
==========

![Widget icon](icons/genexpress.png)

Gives access to [**GenExpress**](https://www.genialis.com/genexpress/) databases.

Signals
-------

**Inputs**:

- (None)

**Outputs**:

- **Data**

  Selected experiments. Each annotated column contains results
  of a single experiment or, if the corresponding option is
  chosen, the average of multiple replicates.

Description
-----------

**GenExpress** is a widget for a direct access to [**GenExpress**](https://www.genialis.com/genexpress/)
database. It is very similar to the **PIPAx** and **GEO Data Sets** widgets as it allows you 
to download the data from selected experiments.

![GenExpress widget](images/GenExpress-stamped.png)

1. Choose a projects to source your data from.
2. Use *Selection bookmarks* to save a selection: select experiments, click the "**+**" button 
   and name the set. To add experiments to your set, click on the set name, select additional 
   experiments and click *Update*. To remove the set click "**-**".
3. In *Sort output columns* set the attributes by which the output columns are sorted. Add 
   attributes with a "+" button and remove them with "-". Switch the sorting order with arrows on the right.
4. Set the expression type for your output data.
   - **Expression RPKM** outputs data in *reads per kilobase of transcript per million mapped reads*
   - **Expression RPKUM** outputs only RPKUM data.
   - **Read counts (raw)** outputs raw read count data. <br>The polyA variants use only polyA (mRNA) mapped hits.
5. **Exclude labels with constant values** removes labels that are the same for all selected experiments.<br>
   **Average replicates (use median)** averages identical experiments by using medians as values.<br>
   **Logarithmic (base 2) transformation** returns log<sub>2</sub>(value+1) for each value.
6. Click *Commit* to output selected data.
7. Select the server you wish to access the data from. Log in to access private data.
8. *Clear cache* removes the uploaded data sets from internal memory.
9. Experiments can be filtered with the *Search* box. To select which attributes to display right-click 
   on the header. To select multiple experiments click them while holding the *Control/Command* key.

Example
-------

In the schema below we connected **GenExpress** to **Data Table** to view the gene expression reads
and then to **Scatter Plot**, where we chose to view expression levels from two experiments. In the plot
we select an outlier and view it in another **Data Table**.

![](images/GenExpress-Example.png)
MA Plot
=======

![image](icons/ma-plot.png)

Visualization of intensity-dependent ratios of raw microarray data.

Signals
-------

**Inputs**:

- **Expression Array**

  DNA microarray.

**Outputs**:

- **Normalized Expression Array**

  Lowess-normalized microarray.

- **Filtered Expression Array**

  Selected instances (in the Z-score cutoff).

Description
-----------

[**MA Plot**](https://en.wikipedia.org/wiki/MA_plot) is a graphical method for visualizing intensity-dependent
ratio of raw mircoarray data. The A represents the average log intensity of the gene
expression (x-axis in the plot), while M stands for the binary log of intensity ratio (y-axis). The widget
outputs either normalized data (Lowess normalization method) or instances above the Z-score cutoff line (instances
with meaningful fold changes).

![image](images/MAplot5-stamped.png)

1. Information on the input data.
2. Select the attribute to split the plot by.
3. Center the plot using:
   - **average**
   - [**Lowess (fast-interpolated)**](https://en.wikipedia.org/wiki/Local_regression) normalization method
   - **Lowess** normalization method
4. Merge replicated by:
   - **average**
   - **median**
   - **geometric mean**
5. Set the **Z-score cutoff** threshold. Z-score is your confidence interval and it is set to
   95% by default. If the widget is set to output *filtered expression array*, instances above the
   [Z-score](https://en.wikipedia.org/wiki/Standard_score) threshold will be in the output (red dots in the plot).
6. Ticking the *Append Z-scores* will add an additional meta attribute with Z-scores to your output data.<br>
   Ticking the *Append Log ratio and Intensity values* will add two additional meta attributes with M and A values
   to your output data.
7. If *Auto commit is on*, the widget will automatically apply changes to the output. Alternatively click *Commit*.

Example
-------

**MA Plot** is a great widget for data visualization and selection. First we select *Caffeine effect: time
course and dose response* data from the **GEO Data Sets** widget and feed it to **MA Plot**. In the plot
we see intensity ratios for a selected experiment variable.

We often need to normalize the experiment data to avoid systematic biases, thus we select *Lowess (fast-interpolated)*
in the *Center Fold-change* box. By ticking both boxes in the *Output* subsection, we get three new meta attributes 
appended - Z-score, Log ratio and Intensity. We see these new attributes and normalized instances in the 
**Data Table** (normalized).

Another possible output for the MA plot widget is *Filtered expression array*, which will give us instances
above the Z-score cutoff threshold (red dots in the plot). We observe these instances the **Data Table** (filtered).

![](images/MAPlot-Example.png)
dictyExpress
=============
    

Gives access to [**dictyExpress**](https://dictyexpress.research.bcm.edu) databases.

Signals
-------

**Inputs**:

- (None)

**Outputs**:

- **Data**

  Selected experiment (time-course gene expression data).

Description
-----------

**dictyExpress** widget gives a direct access to [**dictyExpress**](https://dictyexpress.research.bcm.edu)
database. It allows you to download the data from selected experiments.

![dictyExpress widget](images/GenialisdictyExpress.png)


Example
-------

![](images/dictyExample.png)
Differential Expression
=======================

![DiffExpress widget icon](icons/differential-expression.png)

Plots differential gene expression for selected experiments.

Signals
-------

**Inputs**:

- **Data**

  Data set.

**Outputs**:

- **Selected data**

  Data subset.

Description
-----------

This widget plots a [differential gene expression](http://www.ncbi.nlm.nih.gov/books/NBK10061/) graph for a
sample target. It takes gene expression data as an input (from **dictyExpress**, **PIPAx**, etc.) and outputs a
selected data subset (normally the most interesting genes).

![image](images/DiffExpression-stamped.png)

1. Information of the data input and output. The first line shows the number of samples and genes in the data set.
   The second line displays the selected sample target (read around which the graph is plotted). The third line
   shows the number of undefined gene (missing data) and the fourth the number of genes in the output.
2. Select the plotting method in *Scoring method*:
   - [**Fold change**](https://en.wikipedia.org/wiki/Fold_change): final to initial value ratio
   - **log2 (fold change)**: binary logarithmic transformation of fold change values
   - [**T-test**](https://en.wikipedia.org/wiki/Student%27s_t-test#Independent_two-sample_t-test): parametric test of null hypothesis
   - **T-test (P-value)**: parametric test of null hypothesis with [P-value](https://en.wikipedia.org/wiki/P-value) as criterium
   - [**ANOVA**](https://en.wikipedia.org/wiki/Analysis_of_variance): variance distribution
   - **ANOVA (P-value)**: variance distribution with P-value as criterium
   - [**Signal to Noise Ratio**](https://en.wikipedia.org/wiki/Signal-to-noise_ratio): biological signal to noise ratio
   - [**Mann-Whitney**](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test): non-parametric test of null hypothesis with P-value as criterium
3. Select *Target Labels*. Labels depend on the attributes in the input. In *Values* you can change the sample target
   (default value is the first value on the list, alphabetically or numerically).
4. *Selection* box controls the output data.
   - By setting the *Lower threshold* and *Upper threshold* values you
   are outputting the data outside this interval (the most interesting expression levels). You can also manually place
   the threshold lines by dragging left or right in the plot.
   - If you click *Compute null distribution* box, the widget
   will calculate null distribution and display it in the plot. *Permutations* field allows you to set the precision of
   null distribution (the more permutations the more precise the distribution), while [*&alpha;-value*](https://en.wikipedia.org/wiki/Type_I_and_type_II_errors#Type_I_error) will
   be the allowed probability of false positives. Press *Select* to output this data.
  - The final option is to set the number of best ranked genes and output them with *Select*.
5. When *Auto commit is on* is ticked, the widget will automatically apply the changes. Alternatively press *Commit*. If the *Add gene scores to output* is ticked, the widget will append an additional column with gene scores to the data.

Example
-------

In the example below we chose two experiments from the **PIPAx** widget ( 8 experiments measuring gene expression 
levels on *Dictyostelium discoideum* at different timepoints) and
observed them in the **Data Table**. Then we used the **Differential Expression** widget to select the most interesting
genes. We left upper and lower threshold at default (1 and -1) and output the data. 
Then we observed the selected data in another **Data Table**. As we have ticked
the *Add gene scores to output*, the table shows an additional column with gene scores as instances.

![](images/DiffExpression-Example.png)
Expression Profile Distances
============================

![image](icons/expression-profile-distances.png)

Computes distances between gene expression levels.

Signals
-------

**Inputs**:

- **Data**

  Data set.

**Outputs**:

- **Distances**

  Distance matrix.
  
- **Sorted Data**

  Data with groups as attributes.

Description
-----------

Widget **Expression Profile Distances** computes distances between expression levels among groups of data.
Groups are data clusters set by the user through *separate by* function in the widget. Data can be separated by one or
more variable labels (usually timepoint, replicates, IDs, etc.). Widget outputs distance matrix that can be fed into
**Distance Map** and **Hierarchical Clustering** widgets.

![Distances Widget](images/ExpressionProfileDistances3-stamped.png)

1. Information on the input data.
2. Separate the experiments into groups by labels (normally timepoint, replicates, data name, etc.).
3. Sort the experiments inside the group by labels.
4. Choose the *Distance Measure*:
    - [**Pearson**](https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient) (linear correlation between the values, remapped as a distance in a [0, 1] interval)
    - [**Euclidean**](https://en.wikipedia.org/wiki/Euclidean_distance) ("straight line", distance between two points)
    - [**Spearman**](https://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient) (linear correlation between the rank of the values, remapped as a distance in a [0, 1] interval)
5. If *Auto commit is on*, the widget will automatically compute the distances and output them. Alternatively click *Commit*.
6. This snapshot shows 4 groups of experiments (tp=0, tp=6, tp=12, tp=18) with 2 experiments (replicates) in each group.

Example
-------

**Expression Profile Distances** widget is used to calculate distances between gene expression values
sorted by labels. We chose 8 experiments measuring gene expression levels on *Dictyostelium discoideum* at
different timepoints. In the **Expression Profile Distances** widget we separated the data by timepoint
and sorted them by replicates. We could see the grouping immediately in the *Groups* box on the right. Then we
fed the results to **Distance Map** and **Hierarchical Clustering** to visualize the distances 
and cluster the attributes.


![](images/ExpressionProfileDistances-Example.png)
Select Genes
============

![Select genes icon](icons/select-genes.png)

Manual selection of gene subset.

Signals
-------

**Inputs**:

- **Data**

  Data set.

- **Gene Subset**

  A subset of genes to be used for gene selection (optional).

**Outputs**:

- **Selected Data**

  A subset of genes selected in the widget.

Description
-----------

**Select Genes** widget is used to manually create the gene subset. There are three ways to select genes:
- Manual gene selection (written input). The widget supports autocompletion for gene names.
- Selecting genes from gene sets in the "+" option.
- Selecting genes from a separate input (input can be adjusted in the widget).

![image](images/SelectGenes2-stamped.png)

1. Select *Gene Attribute* if there is more than one column with gene names.
2. Specify how you want to select your genes:
   *Select Genes from 'Gene Subset' input* adds genes from the separate input to selected genes.
   To create a new saved selection, click *Copy genes to saved
   subsets*. The genes will be listed in *Select Genes* text area below. To add these
   genes to an existing selection, click *Append genes to current saved selection*.
3. In *Select specified genes* you can type the gene name and the widget will automatically suggest corresponding
   genes. Genes that match the genes in 
   the input will be colored blue, while the unmatched will remain black.
4. The "+" button has a drop-down menu with two options.
   - *Import names from gene sets...* gives a list of gene sets and copies genes from selected sets into the list.
   - *Import names from text files...* imports gene names from the file.
5. *More* has two settings: *Complete on gene symbol names* (for easier
   gene selection) and *Translate all names to official symbol names* (for uniformity).
6. Set the organism to select the genes from (organism from the input data is chosen as default).
7. *Saved Selections* saves the most frequently used genes. "+" adds a new
   selection, "-" removes the existing one, and *Save* saves the current list. Double-click the selection to rename it.
8. *Output* for this widget is a data subset. If you wish to preserve the order of instances from
   your input data, tick the *Preserve input order* box. If *Auto commit is on*, all changes will
   be communicated automatically. Alternatively press *Commit*.

Below is a screenshot of the *Import Gene Set Names* option.

![image](images/SelectGenes4.png)

Example
-------

Below is a very simple workflow for this widget. We selected *AX4 Dictyostelium discoideum* data from
different time points and two different replicates from **PIPAx** widget. In **Select Genes**
we used the *Import names from gene sets...* option and selected two mRNA processes that gave us
a list of genes you can see in the *Select Genes* box. Then we fed these data into the **Data Table**.
There are 125 genes in the entire *AX4 Dictyostelium discoideum* data that are present in the
selected mRNA processes.

![](images/SelectGenes-Example.png)
dictyExpress
============

![Widget icon](icons/dictyexpress.png)

Gives access to [**dictyExpress**](http://dictyexpress.biolab.si/) databases.

Signals
-------

**Inputs**:

- (None)

**Outputs**:

- **Data**

  Selected experiments. Each annotated column contains results
  of a single experiment or, if the corresponding option is
  chosen, the average of multiple replicates.

Description
-----------

**dictyExpress** is a widget for a direct access to [**dictyExpress**](http://dictyexpress.biolab.si/) database 
and it is very similar to the **GenExpress** and **GEO Data Sets** widgets as it allows you to dowload 
selected experiments.

![dicty widget](images/dictyExpress-stamped.png)

1. The widget will automatically save (cache) downloaded data, which makes them available also in the offline mode. To reset    the widget click *Clear cache*.
2. *Exclude labels with constant values* removes labels that are the same for all the selected experiments in the output.
3. Click *Commit* to output the data.
4. Publicly available data are accessible from the outset. Use *Token* to access password protected data.
5. Available experiments can be filtered with the *Search* box at the top.

Example
-------

In the schema below we connected **ditcyExpress** to a **Data Table** to observe all of
the selected experiments. Then we used **Differential Expression** widget to select
the most relevant genes and output them to another **Data Table**.

![](images/dictyExpress-Example.png)
Gene Info
=========

![image](icons/gene-info.png)

Displays information on the genes in the input.

Signals
-------

**Inputs**:

- **Data**

  Data set.

**Outputs**:

- **Selected Data**

  Instances with meta data that the user has manually selected in the widget.

Description
-----------

A useful widget that presents information on the genes from the [NCBI database](http://www.ncbi.nlm.nih.gov/gene).
You can also select a subset and feed it to other widgets. By clicking on the gene NCBI ID in the list, you will
be taken to the NCBI site with the information on the gene.

![image](images/GeneInfo1-stamped.png)

1. Information on data set size and genes that matched the NCBI ID's.
2. Select the organism of reference.
3. Set the source of gene names. If your gene names are placed as attributes names, select *Use attribute names*.
4. If *Auto commit is on*, changes will be communicated automatically. Alternatively click *Commit*.
5. In the row above the list you can filter the genes by search word(s). If you wish to output the filtered data,
   click *Select Filtered*.
6. If you wish to start from scratch, click *Clear Selection*.

Example
-------

Below we first view the entire *Caffeine effect: time course and dose response* data set in the *Data Table*
widget. Then we feed the same data into the *Gene Info*, where we select only the genes that are located
on the 11th chromosome. We can observe these data in another *Data Table*, where additional information
on the selected genes are appended as meta attributes.

![](images/GeneInfo-Example.png)
BioMart
=======

![BioMart widget icon](icons/biomart.png)

Gives access to [**BioMart**](http://www.biomart.org/news.html) databases.

Signals
-------

**Inputs**:

- None.

**Outputs**:

- **Data**

  Data set.

Description
-----------

**BioMart** is a widget for  direct access to [**BioMart**](http://www.biomart.org/news.html) databases. It
sources data from BioMart, filters it by categories (gene, region, phenotype, gene ontology, etc.) and
appends selected attributes in the output (IDs, sources, strains, etc.). Read more on the BioMart database
library [here](http://nar.oxfordjournals.org/content/43/W1/W589.full.pdf+html).

![image](images/biomart-stamped.png)

1. Clear cached data.
2. Select the database to source your data from.
3. Select the dataset (organism) to source your genes from.
4. If *Unique results only* is ticked, the widget will prevent data duplication. Click *Get results* to output the data.
5. Set the output:
   - in **Attributes** you set the meta data you wish to output (e.g. IDs, sources, strains...).
   - in **Filter** you filter the data by gene, phenotype, ontology, protein domains, etc.

Example
-------

**BioMart** is a great widget for appending additional information to your data. We used *brown-selected* data
in the **File** widget. Then we selected *Ensembl genes 81 (Sanger UK)* database to source our additional data
from. We decided to append *Ensembl Gene ID*, *Ensembl Transcript ID*, *gene type* and *PDB ID*. We also filtered
the data to output only those genes that can be found on chromosome I. We got 9997 instances with 4 meta attributes.
Then we used **Merge Data** widget to append these metas to our data. We matched the data by gene/Ensemble gene ID
and in the end we got a merged data table with 5 meta attributes.

![](images/BioMart-Example.png)
PIPAx
=====

![Widget icon](icons/pipax.png)

Gives access to [**PIPA**](http://pipa.biolab.si/hp/index.html#) databases.

Signals
-------

**Inputs**:

- (None)

**Outputs**:

- **Data**

  Selected experiments. Each annotated column contains results
  of a single experiment or, if the corresponding option is
  chosen, the average of multiple replicates.

Description
-----------

**PIPAx** is a widget for a direct access to [**PIPA**](http://pipa.biolab.si/hp/index.html#) database.
It is very similar to the **GenExpress** and **GEO Data Sets** widgets as it allows you to download the data from 
selected experiments.

![PIPA widget](images/PIPAx-stamped.png)

1. Reloads the experiment data.
2. The widget will save (cache) downloaded data, which makes them also available offline. To reset the widget click *Clear cache*.
3. Use *Experiment Sets* to save a selection:
   select the experiments, click the "**+**" button and name the
   set. To add experiments to the set, click on its name, select
   additional experiments and click *Update*.<br>To remove the set click "**-**".
4. In *Sort output columns* set the attributes by which the output columns are sorted.
   Add attributes with a "+" button and remove them
   with "-". Switch the sorting order with arrows on the right.
5. Set the expression type for your output data.
   - **Raw expression** outputs raw experiment data
   - **RPKM expression** outputs data in *reads per kilobase of transcript per million mapped reads*
   - **RPKM expression + mapability expression** uses similar normalization, but divides with gene 
     mapability instead of exon lengths.<br>The polyA variants use only polyA (mRNA) mapped hits.
6. **Exclude labels with constant values** removes attribute labels that are the same for all selected 
    experiments from the output data.<br>
   **Average replicates (use median)** averages identical experiments by using medians as values.<br>
   **Logarithmic (base 2) transformation** computes the log<sub>2</sub>(value+1) for each value.
7. Click *Commit* to output selected experiments.
8. Log in to access private data.
9. Experiments can be filtered with the *Search* box.
   To select which attributes to display right-click on the header. To select multiple experiments 
   click them while holding the *Control/Command* key.

Example
-------

In the schema below we connected **PIPAx** to **Data Table**, **Set Enrichment**, and **Distance Map**
(through **Distances**) widgets.

![](images/PIPA-Example.png)

The **Data Table** widget above contains the output from the **PIPAx** widget.
Each column contains gene expressions of a single experiment. The labels
are shown in the table header. The **Distance Map** widget shows distances between experiments. The
distances are measured with **Distance** widget, which was set to
compute *Euclidean* distances.
Volcano Plot
============

![image](icons/volcano-plot.png)

Plots significance versus fold-change for gene expression rates.

Signals
-------

**Inputs**:

- **Data**

  Input data set.

**Outputs**:

- **Selected data**

  Data subset.

Description
-----------

[**Volcano plot**](https://en.wikipedia.org/wiki/Volcano_plot_(statistics)) is a graphical method for 
visualizing changes in replicate data. The widget plots a binary logarithm of fold-change on the x-axis versus
[statistical significance](https://en.wikipedia.org/wiki/Statistical_significance) 
(negative base 10 logarithm of p-value) on the y-axis. 

**Volcano Plot** is useful for a quick visual identification of statistically significant
data (genes). Genes that are highly dysregulated are
farther to the left and right, while highly significant fold changes appear higher on the plot.
A combination of the two are those genes that are statistically significant - the widget selects 
the top-ranking genes within the top right and left fields by default.

![image](images/VolcanoPlot-stamped.png)

1. Information on the input and output data.
2. Select *Target Labels*. Labels depend on the attributes in the input. In *Values* you can 
   change the sample target (default value is the first value on the list, alphabetically or numerically).
3. Change the *Settings*: adjust the symbol size and turn off symmetrical selection of the output
   data (the widget selects statistically significant instances by default).
4. If *Auto commit is on* the widget will automatically apply the changes. Alternatively click *Commit*.
5. Visualization of the changes in gene expression. The red lines represent the area with the
   most statistically significant instances. Symmetrical selection is chosen by default, but you can
   also manually adjust the area you want in the output.

Example
-------

Below you can see a simple workflow for **Volcano Plot**. We use *Caffeine effect: time course and dose
response* data from **GEO Data Sets** widget and visualize them in a **Data Table**. We have
6378 gene in the input, so it is essential to prune the data and analyse only those genes
that are statistically significant. **Volcano Plot** helps us do exactly that. Once the
desired area is selected in the plot, we output the data and observe them in another **Data Table**.
Now we get only 80 instances, which were those genes that had a high normalized fold change under
high dose of caffeine and had a low p-value at the same time.

![](images/VolcanoPlot-Example.png)
Data Profiles
=============

![Data profiles widget icon](icons/data-profiles.png)

Plots gene expression levels by attribute in a graph.

Signals
-------

**Inputs**:

- **Data**

  Data set.

**Outputs**:

- **Selected Data**

  Instances that the user has manually selected from the plot.

Description
-----------

**Data Profiles** plots gene expression levels for each attribute in a graph. The default
graph displays the mean expression level for the input data set. The x-axis represents
attributes and the y-axis gene expression values. By hovering over the line you can see
which gene it represents and by click on the line you will select the gene and output it.

![image](images/DataProfiles-stamped.png)

1. Information on the input data.
2. Select display options:
   - **Expression Profiles** will display expression levels for individual data instances.
   - **Quartiles** will show quartile cut-off points.
3. If the data has classes, you can select which class to display by clicking on it. Such data will
   also be colored by class. *Unselect All* will show an empty plot, while *Select All* will diplay
   all data instances by class.
4. Select which attribute you wish to use as a profile label.
5. If *Auto commit is on*, the widget will automatically apply changes to the output. Alternatively click *Commit*.

Example
-------

**Data Profiles** is a great widget for visualizing significant gene expression levels,
especially if the data has been sourced at different timepoints. This allows the user
to see differences in expression levels in time for each instance in the data set and the
overall mean.

Below we used the **PIPAx** widget, where we selected 8 *AX4 Dictyostelium* experiments, all
having been sourced at diffferent timepoints and belonging to one of the two replicates. We
decided to average replicates (to get one instance for both replicates) and to apply logarithmic
transformation to adjust expression levels.

In **Select Genes** we decided to observe only the three genes from the data set that
are a part of the *increased exocytosis* process (lsvB, pldB, amp3), which we selected in
the *Import gene set names* option. This allows us to specify which biological process
we're interested in and to observe only the specified genes.

Then we observe expression levels in **Data Profiles** widget, where we see all three
*Expression Profiles* plotted, together with *Quartiles* and mean expression level. Finally,
we selected the gene with the highest overall expression level and output it to **Data Table**.

![](images/DataProfiles-Example.png)
GO Browser
==========

![GO Browser widget icon](icons/go-browser.png)

Provides access to Gene Ontology database.

Signals
-------

**Inputs**:

- **Cluster Data**

  Data on clustered genes.

- **Reference Data**

  Data with genes for the reference set (optional).

**Outputs**:

- **Data on Selected Genes**

  Data on genes from the selected GO node.

- **Data on Unselected Genes**

  Data on genes from GO nodes that weren't selected.

- **Data on Unknown Genes**

  Data on genes that are not in the GO database.

- **Enrichment Report**

  Data on GO enrichment analysis.


Description
-----------

**GO Browser** widget provides access to [*Gene Ontology database*](http://geneontology.org/). 
Gene Ontology (GO) classifies genes and gene products to terms organized in a graph structure called an ontology.
The widget takes any data on genes as an input (it is best to input statistically significant genes,
for example from the output of the **Differential Expression** widget) and shows a ranked list of GO terms with
p-values. This is a great tool for finding biological processes that are over- or under-represented in a 
particular gene set. The user can filter input data by selecting terms in a list.

![image](images/GObrowser5-stamped.png)

**INPUT tab**<br>

1. Information on the input data set. *Ontology/Annotation Info* reports the current status of the GO database.
2. Select organism for the GO term analysis.
3. Use this attribute to extract gene names for the input data. You can use attribute names as gene names and 
   adjust gene matching in the *Gene matcher settings* box.
4. Select the reference. You can either have the *entire genome* as reference or a *reference set* from the input.
5. Select the ontology where you want to calculate the enrichment. There are three *Aspect* options:
   - [**Biological process**](http://geneontology.org/page/biological-process-ontology-guidelines)
   - [**Cellular component**](http://geneontology.org/page/cellular-component-ontology-guidelines)
   - [**Molecular function**](http://geneontology.org/page/molecular-function-ontology-guidelines)
6. A ranked tree (upper pane) and list (lower pane) of GO terms for the selected aspect:
   - **GO term**
   - **Cluster**: number of genes from the input that are also annotated to a particular GO term 
     (and its proportion in all the genes from that term).
   - **Reference**: number of genes that are annotated to a particular GO term (and its proportion in the entire genome).
   - **P-value**: probability of seeing as many or more genes at random. The closer the p-value is to zero, the more significant a particular GO term is. Value is written in [e notation](https://en.wikipedia.org/wiki/Scientific_notation#E_notation).
   - **FDR**: [false discovery rate](https://en.wikipedia.org/wiki/False_discovery_rate) - a 
     multiple testing correction that means a proportion of false discoveries among all discoveries up to that FDR value.
   - **Genes**: genes in a biological process.
   - [**Enrichment**](http://geneontology.org/page/go-enrichment-analysis) level

![image](images/GObrowser-tabs-stamped.png)

**FILTER tab**<br>

1. *Filter GO Term Nodes* by:
   - **Genes** is a minimal number of genes mapped to a term
   - **P-value** is a max term p-value
   - **FDR**: is a max term [false discovery rate](https://en.wikipedia.org/wiki/False_discovery_rate)
2. *Significance test* specifies distribution to use for null hypothesis:
   - [**Binomial**](https://en.wikipedia.org/wiki/Binomial_distribution): use a binomial distribution
   - [**Hypergeometric**](https://en.wikipedia.org/wiki/Hypergeometric_distribution): use a hypergeometric distribution
3. [*Evidence codes in annotation*](http://geneontology.org/page/guide-go-evidence-codes) show how the 
   annotation to a particular term is supported.

**SELECT tab**<br>

4. *Annotated genes* outputs genes that are:
   - **Directly or Indirectly** annotated (direct and inherited annotations)
   - **Directly** annotated (inherited annotations won't be in the output)
5. *Output*:
   - **All selected genes**: outputs genes annotated to all selected GO terms
   - **Term-specific genes**: outputs genes that appear in only one of selected GO terms
   - **Common term genes**: outputs genes common to all selected GO terms
   - **Add GO Term as class**: adds GO terms as class attribute

Example
-------

In the example below we have used **GEO Data Sets** widget, in which we have selected 
*Caffeine effects: time course and dose response* data set, and connected it to a **Differential
Analysis**. Differential analysis allows us to select genes with the highest statistical relevance
(we used ANOVA scoring) and feed them to **GO Browser**. This widget lists four biological
processes for our selected genes. Say we are interested in finding out more about *monosaccharide transport*
as this term has a high enrichment rate. To learn more about which genes
are annotated to this GO term we view it in the **Data Table**, where we see all the genes
participating in this process listed.

![](images/GObrowser-Example.png)
GEO Data Sets
=============

![GEO Data Sets widget icon](icons/geo-data-sets.png)

Provides access to data sets from gene expression omnibus ([GEO
DataSets](http://www.ncbi.nlm.nih.gov/gds)).

Signals
-------

**Inputs**:

-   (None)

**Outputs**:

- **Data**

 Data set selected in the widget with genes or samples in rows.

Description
-----------

**[GEO DataSets](http://www.ncbi.nlm.nih.gov/gds)** is a data base of gene
expression curated profiles maintained by [NCBI](http://www.ncbi.nlm.nih.gov/) and included in the [Gene
Expression Omnibus](http://www.ncbi.nlm.nih.gov/geo/info/datasets.html). This Orange widget provides
access to all its data sets and outputs a data set selected for further
processing. For convenience, each dowloaded data set is stored locally.

![GEO Data Sets widget](images/GEOdataset-stamped.png)

1. Information on the GEO data set collection. Cached data sets are the ones currently stored on the computer.
2. Output features. If *Genes or spots* is selected, genes (or spots) will be used as attributes. Alternatively samples
   will be used as attributes. *Merge spots of same gene* averages measures of the same gene. Finally, in the
   *Data set name* you can rename the output data. GEO title will be used as a default name.
3. If *Auto commit is on*, then the selected data set will be automatically communicated to other widgets. Alternatively,
   click *Commit*.
4. *Filter* allows you to search for the data set. Below you see a list of GEO data sets with an ID number (link to the NCBI
   Data Set Browser), title,
   organism used in the experiment, number of samples, features, genes, subsets and a reference number for the PubMed
   journal (link to the article abstract).
5. Short description of the experiment from which the data set is sourced.
6. Select which *Sample Annotations* will be used in the output.

Example
-------

**GEO Data Sets** is similar to the **File** widget. In the example below
we selected *Caffeine effect: time dose and response* data set from the GEO data base and used *Genes or spots* as
attributes. We inspected the data in *Data Table*. Then we selected
3 genes in the **Select Columns** widget for a detailed analysis in another data table.

![](images/GEODataSets-Example2.png)
Set Enrichment
==============

![Set Enrichment widget icon](icons/set-enrichment.png)

Determines statistically significant differences in expression levels for biological processes.

Signals
-------

**Inputs**:

- **Data**

  Data set.

- **Reference**

  Data with genes for the reference set (optional).

**Outputs**:

- **Selected data**

  Data subset.

Description
-----------

The widget shows a ranked list of terms with [p-values](https://en.wikipedia.org/wiki/P-value), 
[FDR](https://en.wikipedia.org/wiki/False_discovery_rate) and 
[enrichment](https://en.wikipedia.org/wiki/Gene_set_enrichment). 
**Set Enrichment** is a great tool for finding biological processes that are over-represented in a particular gene 
or chemical set.

Sets from ([GO](http://geneontology.org/), [KEGG](http://www.genome.jp/kegg/), 
[miRNA](http://www.mirbase.org/) and [MeSH](http://www.nlm.nih.gov/mesh/MBrowser.html)) come with the Orange installation.

![image](images/SetEnrichment1-stamped.png)

1. Information on the input data set and the ratio of genes that were found in the databases.
2. Select the species.
3. *Entity names* define the features in the input data that you wish to use for term analysis. Tick *Use feature names*
   if your genes or chemicals are used as attribute names rather than as meta attributes.
4. Select the reference data. You can either have entities (usually genes from the organism - *All Entities*)
   as a reference or a reference set from the input.
5. Select which *Entity sets* you wish to have displayed in the list.
6. When *Auto commit is on*, the widget will automatically apply the changes. Alternatively press *Commit*. 
7. Filter the list by:
   - the minimum number of **entities** included in each term
   - the minimum threshold for **p-value**
   - the maximum threshold for **false discovery rate**
   - a search word

Example
-------

In the example below we have decided to analyse gene expression levels from *Caffeine effect: time course
and dose response* data set. We used the ANOVA scoring in the **Differential Expression** widget to 
select the most interesting genes. Then we fed those 628 genes to **Set Enrichment** for additional
analysis of the most valuable terms. We sorted the data by FDR values and selected the top-scoring
term. **Heat Map** widget provides a nice visualization of the data.

![](images/SetEnrichment-Example.png)
KEGG Pathways
=============

![KEGG widget icon](icons/kegg-pathways.png)

Diagrams of molecular interactions, reactions, and relations.

Signals
-------

**Inputs**:

- **Data**

  Data set.

- **Reference**

  Referential data set.

**Outputs**:

- **Selected Data**

  Data subset.

- **Unselected Data**

  Remaining data.

Description
-----------

**KEGG Pathways** widget displays diagrams of molecular interactions, reactions and relations from the
[KEGG Pathways Database](http://www.genome.jp/kegg/pathway.html). It takes data on gene expression as an
input, matches the genes to the biological processes and displays a list of corresponding pathways. To
explore the pathway, the user can click on any process from the list or arrange them by P-value to get
the most relevant processes at the top.

![image](images/KEGG2-stamped.png)

1. Information on the input and the ratio of matched genes.
2. Select the organism for term analysis. The widget automatically selects the organism from the input data.
3. Set the attribute to use for gene names. If gene names are your attribute names, tick *Use variable names*.
4. If you have a separate reference set in the input, tick *From signal* to use these data as reference.
5. To have pathways listed and displayed by vertical descent, tick *Show pathways in full orthology*.
6. To fit the image to screen, tick *Resize to fit*. Untick the box if you wish to explore the pathways.
7. To clear all locally cached KEGG data, press *Clear cache*.
8. When *Auto commit is on*, the widget will automatically apply the changes. Alternatively press *Commit*.
9. A list of pathways either as processes or in full orthology. Click on the process to display the pathway.
   You can sort the data by P-value to get the most relevant results at the top.

Example
-------
Orange Bioinformatics
=====================

Orange Bioinformatics extends Orange_, a data mining software
package, with common functionality for bioinformatics. The provided
functionality can be accessed as a Python library or through a visual
programming interface (Orange Canvas). The latter is also suitable for
non-programmers.

In Orange Canvas the analyst connects basic computational units, called
widgets, into data flow analytics schemas. Two units-widgets can be
connected if they share a data type. Compared to other popular tools like
Taverna, Orange widgets are high-level, integrated potentially complex
tasks, but are specific enough to be used independently. Even elaborate
analyses rarely consist of more than ten widgets; while tasks such as
clustering and enrichment analysis could be executed with up to five
widgets. While building the schema each widget is independently controlled
with settings, the settings do not conceptually burden the analyst.

Orange Bioinformatics provides access to publicly available data,
like GEO data sets, Biomart, GO, KEGG, Atlas, ArrayExpress, and PIPAx
database. As for the analytics, there is gene selection, quality control,
scoring distances between experiments with multiple factors. All features
can be combined with powerful visualization, network exploration and
data mining techniques from the Orange data mining framework.

.. _Orange: http://orange.biolab.si/

Documentation: http://orange-bioinformatics.readthedocs.org/
Orange Bioinformatics documentation
===================================

Orange Bioinformatics is an add-on for Orange_ data mining software package. It
extends Orange by providing functionality for some elementary tasks in
bioinformatics, like gene set analysis, enrichment, and access to pathway
libraries. Included are also widgets for with graphical user interface
for, mainly, gene expression analytics.

Orange Bioinformatics provides access to publicly available data, like GEO data
sets, Biomart, GO, KEGG, Atlas, ArrayExpress, and PIPAx database. As for the
analytics, there is gene selection, quality control, scoring distances between
experiments with multiple factors. All features can be combined with powerful
visualization, network exploration and data mining techniques from the Orange
data mining framework.

.. _Orange: http://orange.biolab.si/

Widgets
-------

.. toctree::
   :maxdepth: 1
 
   widgets/biomart
   widgets/databasesupdate
   widgets/dataprofiles
   widgets/dictyExpress
   widgets/differentialexpression
   widgets/expressionprofiledistances
   widgets/geneinfo
   widgets/genexpress
   widgets/geodatasets
   widgets/gobrowser
   widgets/keggpathways
   widgets/maplot
   widgets/PIPAx
   widgets/qualitycontrol
   widgets/selectgenes
   widgets/setenrichment
   widgets/volcanoplot
   widgets/genialisDictyExpress

Scripting Reference
-------------------
   
.. toctree::
   :maxdepth: 1

   reference/arrayexpress.rst
   reference/biomart.rst
   reference/dicty.rst
   reference/dictybase.rst
   reference/gene.rst
   reference/gene.homology.rst
   reference/genesets.rst
   reference/geo.rst
   reference/go.rst
   reference/gsea.rst
   reference/kegg.rst
   reference/omim.rst
   reference/ontology.rst
   reference/ppi.rst
   reference/taxonomy.rst
   reference/utils.stats.rst
   reference/resolwe.rst

Installation
------------

To install Bioinformatics add-on for Orange from PyPi_ run::

    pip install Orange-Bioinformatics

To install it from source code run::

    python setup.py install

To build Python egg run::

    python setup.py bdist_egg

To install add-on in `development mode`_ run::

    python setup.py develop

.. _development mode: http://packages.python.org/distribute/setuptools.html#development-mode
.. _PyPi: http://pypi.python.org/pypi

Source Code and Issue Tracker
-----------------------------

Source code is available on GitHub_.

Bug reports and suggestions should be submitted to `Issue Tracker`_ on GitHub.

.. _GitHub: https://github.com/biolab/orange-bio
.. _Issue Tracker: https://github.com/biolab/orange-bio/issues

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
============================================================
Protein-protein interactions (:mod:`ppi`)
============================================================

.. py:currentmodule:: orangecontrib.bio.ppi

:class:`PPIDatabase` is an abstract class defining a common interface
for accessing protein-protein interaction databases.

Classes implementing this interface are:

  - :class:`BioGRID` for accessing `BioGRID <http://thebiogrid.org>`_
  - :class:`STRING` for accessing `CC`_ licensed `STRING <http://www.string-db.org/>`_
  - :class:`STRINGDetailed` for accessing `CC-NC-SA`_ licensed `STRING <http://www.string-db.org/>`_

.. _`CC`: http://creativecommons.org/licenses/by/3.0/

.. _`CC-NC-SA`: http://creativecommons.org/licenses/by-nc-sa/3.0/

The common interface
---------------------

.. autoclass:: PPIDatabase
   :members:
   :member-order: bysource

PPI databases
------------------

.. autoclass:: BioGRID
   :members:
   :inherited-members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: STRING
   :members:
   :inherited-members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: STRINGDetailed
   :members:
   :inherited-members:
   :member-order: bysource
   :show-inheritance:


.. py:currentmodule:: orangecontrib.bio.omim

.. index:: NCBI
.. index:: OMIM

**************************************************************
OMIM: Online Mendelian Inheritance in Man (:mod:`omim`)
**************************************************************

:obj:`omim` provides an interface to 
`Online Mendelian Inheritance in Man <http://www.ncbi.nlm.nih.gov/omim>`_
database. For now it only supports mapping genes to diseases.

.. autofunction:: diseases

.. autofunction:: genes

.. autofunction:: disease_genes

.. autofunction:: gene_diseases

The following example creates a network of connected diseases and
save it in `Pajek <http://vlado.fmf.uni-lj.si/pub/networks/pajek/>`_
.net format sets information file (:download:`omim_disease_network.py
<code/omim_disease_network.py>`).

.. literalinclude:: code/omim_disease_network.py

.. automodule:: orangecontrib.bio.ontology

.. autoclass:: OBOOntology
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: OBOObject
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Term
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Typedef
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Instance
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: OBOParser
   :members:
   :member-order: bysource
.. py:currentmodule:: orangecontrib.bio.utils.stats

.. index:: binomial distribution
.. index:: hypergeometric distribution
.. index:: FDR
.. index:: Bonferroni

**************************************************************
Probability distributions and corrections (:mod:`utils.stats`)
**************************************************************


.. autoclass:: Binomial
   :members: __call__, p_value

.. autoclass:: Hypergeometric
   :members: __call__, p_value

.. autofunction:: FDR

.. autofunction:: Bonferroni



=========================
Bio Mart (:mod:`biomart`)
=========================

.. py:currentmodule:: orangecontrib.bio.biomart

Access BioMart MartService.

>>> from orangecontrib.bio.biomart import *
>>> connection = BioMartConnection(
...     "http://www.biomart.org/biomart/martservice")
...
>>> reg = BioMartRegistry(connection)
>>> for mart in reg.marts():
...    print mart.name
...
ensembl...
>>> dataset = BioMartDataset(
...     mart="ensembl", internalName="hsapiens_gene_ensembl",
...     virtualSchema="default", connection=connection)
...
>>> for attr in dataset.attributes()[:10]:
...    print attr.name
...
Ensembl Gene ID...
>>> data = dataset.get_data(
...    attributes=["ensembl_gene_id", "ensembl_peptide_id"],
...    filters=[("chromosome_name", "1")])
...
>>> query = BioMartQuery(reg.connection, virtualSchema="default")
>>> query.set_dataset("hsapiens_gene_ensembl")
>>> query.add_attribute("ensembl_gene_id")
>>> query.add_attribute("ensembl_peptide_id")
>>> query.add_filter("chromosome_name", "1")
>>> count = query.get_count()

Interface
---------

.. autoclass:: BioMartConnection
   :members:

.. autoclass:: BioMartRegistry
   :members:

.. autoclass:: BioMartQuery
   :members:

.. autoclass:: BioMartDataset
   :members:

.. autoclass:: BioMartVirtualSchema
   :members:

.. autoclass:: BioMartDatabase
   :members:

.. autoclass:: Attribute
   :members:

.. autoclass:: Filter
   :members:

.. autoclass:: DatasetConfig
   :members:
============================================================
Dictyostelium discoideum databases (:mod:`dicty`)
============================================================

.. py:currentmodule:: orangecontrib.bio.dicty

The following example downloads experiments from the PIPA
database, specifically "RPKM + mapability expression (polyA) - R6"
results for all public experiments on Dictyostelium discoideum (dd)
at time point 16.

.. literalinclude:: code/pipax1.py

PIPAx database
--------------

.. autoclass:: orangecontrib.bio.dicty.PIPAx
   :members: __init__, genomes, get_data, mappings, result_types, results_list

DictyExpress database
---------------------

.. autoclass:: orangecontrib.bio.dicty.DictyExpress
   :members: __init__, objects, annotationTypes, annotations, get_data, search, annotationOptions

Auxillary functionality
-----------------------

.. autoclass:: orangecontrib.bio.dicty.CacheSQLite
    :members: __init__, list, contains, add, get, commit, clear

============================================================
Translation of homologs (:mod:`gene.homology`)
============================================================

.. py:currentmodule:: orangecontrib.bio.gene.homology

This module is an interface to `NCBI HomoloGene <http://www.ncbi.nlm.nih.gov/homologene/>`_.

.. autofunction::  orangecontrib.bio.gene.homology.all_genes

.. autofunction::  orangecontrib.bio.gene.homology.homologs

.. autofunction::  orangecontrib.bio.gene.homology.homolog

.. autofunction::  orangecontrib.bio.gene.homology.orthologs

.. autofunction::  orangecontrib.bio.gene.homology.all_genes_inParanoid

Examples
--------

Mapping homologs from yeast to human:

.. literalinclude:: code/homology-yeast-to-human.py

===================================
Array Express (:mod:`arrayexpress`)
===================================

.. py:currentmodule:: orangecontrib.bio.arrayexpress

Access the `ArrayExpress`_ web services and database.

`ArrayExpress`_ is a database of gene expression experiments
that you can query and download.

Retrieve the object representing experiment with accession E-TABM-25

>>> from orangecontrib.bio import arrayexpress
>>> experiment = ArrayExpressExperiment("E-TABM-25")
>>> print experiment.accession
E-TABM-25

>>> print experiment.name
Transcription profiling of aging in the primate brain

>>> print experiment.species
['Pan troglodytes']

>>> print experiment.files
[{'kind': ...

>>> # Retrieve the data matrix for experiment 'E-MEXP-2917'
>>> experiment = ArrayExpressExperiment("E-MEXP-2917")
>>> table = experiment.fgem_to_table()


Low level Array Express query using REST services:

>>> from orangecontrib.bio import arrayexpress
>>> arrayexpress.query_experiments(accession='E-MEXP-31')
{u'experiments': ...

>>> arrayexpress.query_experiments(keywords='gliobastoma')
{u'experiments': ...

>>> arrayexpress.query_files(accession='E-MEXP-32', format="xml")
<xml.etree.ElementTree.ElementTree object ...

.. note:: Currently querying ArrayExpress files only works with the xml format.

.. _`ArrayExpress`: http://www.ebi.ac.uk/arrayexpress/


Interface
---------

.. autoclass:: ArrayExpressConnection
   :members:

.. autoclass:: ArrayDesign
   :members:

.. autoclass:: SampleDataRelationship
   :members:

.. autoclass:: InvestigationDesign
   :members:


Low-level querying with REST
-----------------------------

.. autofunction:: query_experiments

.. autofunction:: query_files


.. py:currentmodule:: orangecontrib.bio.resolwe

============================================================
Resolwe (:mod:`GenAPI`)
============================================================

Example of GenAPI usage:

.. literalinclude:: code/exampleGenesis.py

::

	Experiment id: 564a509e6b13390ffb40d4c8 - D. purpureum
	Experiment id: 564a54af6b13398f1640d4cd - Filter development
	Experiment id: 564a586b6b13398f1640d4d4 - cAMP pulses
	Experiment id: 56a944016b13395571175e36 - D. fasciculatum (Illumina)
	Experiment id: 56a944936b133909d8175e61 - D. lacteum (Illumina)
	Experiment id: 56a945fc6b133903da175e47 - D. fasciculatum (454)
	Experiment id: 56a946976b133909d8175e62 - P. pallidum (454)
	Experiment id: 56e2f39b6b1339964c33e7f0 - KO - gtaC knockout
	Experiment id: 56e2f50a6b1339964c33e7f1 - Wild Type
	Experiment id: 56e2f5a66b1339964c33e7f2 - CS - Cysteine-substituted strain
	Experiment id: 56e2f6706b1339964c33e7f3 - CM - Complemented mutant
	Experiment id: 56e2fdcc6b1339353d33e7e6 - D. discoideum

Download data from server:

.. literalinclude:: code/exampleGenesis1.py

::

	Time Points [0, 4, 8, 12, 16, 20, 24]
	Gene DPU_G0054708 [6.92022056529, 4.234260142395, 1.23427395022, 1.086050235081, 3.34525253324, 4.469964743405001, 3.26381321097]

Convert Json format to Orange data.Table:

.. literalinclude:: code/exampleGenesis2.py

::

	[71.582, 57.228, 57.609, 69.851, 65.566, 47.757, 27.165] {DPU_G0064618},
	 [181.716, 118.011, 100.540, 115.761, 53.397, 64.811, 56.353] {DPU_G0075396},
	 [3.382, 8.086, 1.581, 0.239, 2.893, 1.192, 1.781] {DPU_G0062146},
	 [0.758, 1.304, 0.678, 0.438, 1.509, 0.812, 0.751] {DPU_G0069396},
	 [629.349, 1290.830, 1493.227, 457.539, 15.799, 33.559, 114.367] {DPU_G0061882},
	 [0.665, 3.038, 2.136, 3.498, 5.774, 5.249, 5.878] {DPU_G0053074},
	 [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000] {DPU_G0055934},
	 [13.860, 43.443, 121.232, 272.804, 254.022, 156.572, 339.583] {DPU_G0053488},
	 [0.553, 0.495, 10.662, 35.259, 14.985, 22.432, 51.782] {DPU_G0070984},
	 [68.098, 55.382, 13.394, 13.612, 23.304, 17.440, 20.316] {DPU_G0062718}

.. autoclass:: orangecontrib.bio.resolwe.GenAPI
    :members: __init__, fetch_etc_objects, download_etc_data, etc_to_table
    :member-order: bysource
.. py:currentmodule:: orangecontrib.bio.geneset
.. py:module:: orangecontrib.bio.geneset

.. index:: gene set
.. index:: gene sets

********************************************************
Gene sets (:mod:`geneset`)
********************************************************

This module can load either gene sets distributed with Orange
or custom gene sets in the `GMT file format <http://www.molmine.com/magma/fileformats.html>`_.

The available gene set collection can be listed with :obj:`list_all`.

:obj:`collections` loads gene sets. Gene sets provided with Orange are
organized hierarchically. Although the ``GO`` hierarchy includes subsets,
all of them can be loaded with (the organism here is a mouse)::

    orangecontrib.bio.geneset.collections((("GO",), "10090"))

To open multiple gene set collections at once, for example, KEGG and GO, try::

    orangecontrib.bio.geneset.collections((("KEGG",), "10090"), (("GO",), "10090"))

You could also open a file with gene sets. The following line would open
``specific.gmt`` from the current working directory::

    orangecontrib.bio.geneset.collections("specific.gmt")

The above examples combined::

    orangecontrib.bio.geneset.collections((("KEGG",), "10090"), (("GO",), "10090"), "specific.gmt")

Furthermore, all gene sets for a specific organism can be opened with an empty
hierarchy::

    orangecontrib.bio.geneset.collections((tuple(), "10090"))


Loading gene sets
=================

.. autofunction:: orangecontrib.bio.geneset.list_all

.. autofunction:: orangecontrib.bio.geneset.collections


Supporting functionality
========================

.. autoclass:: orangecontrib.bio.geneset.GeneSets
   :members:
   :show-inheritance:

.. autoclass:: orangecontrib.bio.geneset.GeneSet
   :members:

.. autofunction:: orangecontrib.bio.geneset.register

===================================
Organism Taxonomy (:mod:`taxonomy`)
===================================

.. py:currentmodule:: orangecontrib.bio.taxonomy

This module provides access to the `NCBI's organism taxonomy information
<http://www.ncbi.nlm.nih.gov/Taxonomy/>`_ and organism name unification
across different modules.

.. autofunction:: orangecontrib.bio.taxonomy.name(taxid)

.. autofunction:: orangecontrib.bio.taxonomy.other_names(taxid)

.. autofunction:: orangecontrib.bio.taxonomy.search(string, onlySpecies=True, exact=False)

.. autofunction:: orangecontrib.bio.taxonomy.lineage(taxid)
   
..
   ..autofunction:: orangecontrib.bio.taxonomy.to_taxid

.. autofunction:: orangecontrib.bio.taxonomy.taxids

.. autofunction:: orangecontrib.bio.taxonomy.common_taxids

.. autofunction:: orangecontrib.bio.taxonomy.essential_taxids

Examples
--------

The following script takes the list of taxonomy IDs and prints out their name:

.. literalinclude:: code/taxonomy1.py

The output of the script is::

    3702   Arabidopsis thaliana
    9913   Bos taurus
    6239   Caenorhabditis elegans
    3055   Chlamydomonas reinhardtii
    7955   Danio rerio
    352472 Dictyostelium discoideum AX4
    7227   Drosophila melanogaster
    562    Escherichia coli
    11103  Hepatitis C virus
    9606   Homo sapiens
    10090  Mus musculus
    2104   Mycoplasma pneumoniae
    4530   Oryza sativa
    5833   Plasmodium falciparum
    4754   Pneumocystis carinii
    10116  Rattus norvegicus
    4932   Saccharomyces cerevisiae
    4896   Schizosaccharomyces pombe
    31033  Takifugu rubripes
    8355   Xenopus laevis
    4577   Zea mays


============================================================
Gene Set Enrichment Analysis (GSEA, :mod:`gsea`)
============================================================

.. py:currentmodule:: orangecontrib.bio.gsea


Gene Set Enrichment Analysis (GSEA) [GSEA]_ aims to identify enriched gene sets
given gene expression data for multiple samples with their phenotypes.

.. autofunction:: orangecontrib.bio.gsea.run

.. autofunction:: orangecontrib.bio.gsea.direct

Examples: gene expression data
------------------------------

The following examples use a gene expression data set from the GEO database.
We show the same analysis on two formats of data.

With samples as instances (in rows):

.. literalinclude:: code/gsea_instances.py

With samples as features (in columns):

.. literalinclude:: code/gsea_genes.py

Both scripts output::

    GSEA results (descriptor: tissue)
    LABEL                                       NES    FDR   SIZE MATCHED
    Porphyrin and chlorophyll meta           -1.817  0.000     43      23
    Staphylococcus aureus infectio           -1.998  0.000     59      28
    Non-homologous end-joining                1.812  0.000     13      12
    Fanconi anemia pathway                    1.911  0.000     53      27
    Cell cycle                                1.777  0.000    124     106
    Glycine, serine and threonine            -2.068  0.000     39      29
    HIF-1 signaling pathway                  -1.746  0.000    106      90
    Ether lipid metabolism                   -1.788  0.000     42      27
    Fc epsilon RI signaling pathwa           -1.743  0.000     70      53
    B cell receptor signaling path           -1.782  0.000     72      62

Example: our own gene sets
--------------------------

We present a simple example on iris data set. Because data set
is not a gene expression data set, we had to specify our own
sets of features that belong together.

.. literalinclude:: code/gsea1.py

The output::

    LABEL     NES  P-VAL GENES
    sepal   1.087  0.630 ['sepal width', 'sepal length']
    petal  -1.117  0.771 ['petal width', 'petal length']


Example: directly passing correlation data
------------------------------------------

GSEA can also directly use correlation data between individual genes
and a phenotype. If (1) input data with only one example (attribute names 
are gene names) or (2) there is only one continuous feature in the given 
data set (gene names are in the first :obj:`Orange.feature.String`.
The following outputs ten pathways with smallest p-values.

.. literalinclude:: code/gsea2.py

The output::

    LABEL                                       NES  P-VAL   SIZE MATCHED
    Biosynthesis of amino acids               1.407  0.056     58      40
    beta-Alanine metabolism                   1.165  0.232     13      10
    Taurine and hypotaurine metabolism        1.160  0.413      4       3
    Porphyrin and chlorophyll metabolis      -0.990  0.517     14       5
    Valine, leucine and isoleucine degr       0.897  0.585     29      21
    Ether lipid metabolism                    0.713  0.857     10       6
    Biosynthesis of unsaturated fatty a       0.659  0.922     10       6
    Protein processing in endoplasmic r       0.647  0.941     71      40
    RNA polymerase                            0.550  0.943     24       7
    Glycosylphosphatidylinositol(GPI)-a      -0.540  0.946     19       4

.. [GSEA] Subramanian, Aravind et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 2005.
.. py:currentmodule:: orangecontrib.bio.gene

.. index:: gene matching
.. index:: gene name matching
.. index:: matching
.. index:: NCBI
.. index:: gene info

********************************************************
Gene name matching (:mod:`gene`)
********************************************************

Genes can have multiple aliases. When we combine data from different
sources, for example expression data with GO gene sets, we have to
match gene aliases representing the same genes. All implemented matching
methods are based on sets of gene aliases.

Gene matchers in this module match genes to a user-specified set of target
gene names. For gene matching, initialize a gene matcher (:obj:`Matcher`),
set the target gene names with :obj:`~Matcher.set_targets`, and then
match with :obj:`~Matcher.match` or :obj:`~Matcher.umatch` functions. The
following example (:download:`genematch1.py <code/genematch1.py>`)
matches gene names to NCBI gene IDs:

.. literalinclude:: code/genematch1.py

Gene name matching
==================

The base class for all the following gene matcher is :obj:`Matcher`.

.. autoclass:: Matcher
   :members:

This modules provides the following gene matchers:

.. autoclass:: MatcherAliasesKEGG

.. autoclass:: MatcherAliasesGO

.. autoclass:: MatcherAliasesDictyBase

.. autoclass:: MatcherAliasesNCBI

.. autoclass:: MatcherAliasesEnsembl

.. autoclass:: MatcherDirect

Gene name matchers can be applied in sequence (until the first match) or
combined (overlapping sets of gene aliases of multiple gene matchers
are combined) with the :obj:`matcher` function.

.. autofunction:: matcher

The following example tries to match input genes onto KEGG gene aliases
(:download:`genematch2.py <code/genematch2.py>`).

.. literalinclude:: code/genematch2.py

Results show that GO aliases can not match onto KEGG gene IDs. For the
last gene only joined GO and KEGG aliases produce a match::

        gene         KEGG           GO      KEGG+GO
        cct7    hsa:10574         None    hsa:10574
        pls1     hsa:5357         None     hsa:5357
        gdi1     hsa:2664         None     hsa:2664
       nfkb2     hsa:4791         None     hsa:4791
      a2a299         None         None     hsa:7052


The following example finds KEGG pathways with given genes
(:download:`genematch_path.py <code/genematch_path.py>`).

.. literalinclude:: code/genematch_path.py

Output::

    Fndc4 is in
      /
    Itgb8 is in
      PI3K-Akt signaling pathway
      Focal adhesion
      ECM-receptor interaction
      Cell adhesion molecules (CAMs)
      Regulation of actin cytoskeleton
      Hypertrophic cardiomyopathy (HCM)
      Arrhythmogenic right ventricular cardiomyopathy (ARVC)
      Dilated cardiomyopathy
    Cdc34 is in
      Ubiquitin mediated proteolysis
      Herpes simplex infection
    Olfr1403 is in
      Olfactory transduction

.. py:currentmodule:: orangecontrib.bio.dicty.phenotypes

.. index:: mutant phenotypes
.. index:: phenotypes
.. index::
   single: D. dictyostelium; mutants

*********************************************************
D. discoideum Mutant Phenotypes (:mod:`dicty.phenotypes`)
*********************************************************

This modules provides an interface to `Dictyostelium mutant
phenotypes <http://dictybase.org/Downloads/>`_ data from the
`dictyBase <http://dictybase.org/>`_.  The mutants are presented as
:obj:`DictyMutant` objects with their respective name, strain descriptor,
associated genes and associated phenotypes.

>>> from orangecontrib.bio.dicty.phenotypes import *
>>> # Create a set of all mutant objects
>>> dicty_mutants = mutants()
>>> # List a set of all genes referenced by a single mutant
>>> print mutant_genes(dicty_mutants[0])
['cbfA']
>>> # List a set of all phenotypes referenced by a single mutant
>>> print mutant_phenotypes(dicty_mutants[0])
['aberrant protein localization']

Classes and Functions
=====================

.. automodule:: orangecontrib.bio.dicty.phenotypes
   :members:
   :member-order: bysource


.. py:currentmodule:: orangecontrib.bio.geo

.. index:: NCBI
.. index:: GEO
.. index:: Gene Expression Omnibus
.. index:: microarray data sets

**************************************************************
NCBI's Gene Expression Omnibus interface (:mod:`geo`)
**************************************************************

This module provides an interface to `NCBI
<http://www.ncbi.nlm.nih.gov/>`_'s `Gene Expression Omnibus
<http://www.ncbi.nlm.nih.gov/geo/>`_ repository. It 
supports `GEO DataSets <http://www.ncbi.nlm.nih.gov/sites/GDSbrowser>`_
query and retrieval.

In the following example :obj:`GDS.getdata`
construct a data set with genes in rows and samples in
columns. Notice that the annotation about each sample is retained
in ``.attributes``.

>>> import orangecontrib.bio.geo
>>> gds = orangecontrib.bio.geo.GDS("GDS1676")
>>> data = gds.getdata()
>>> len(data)
717
>>> data[200]
[1.407, 1.268, 1.969, 1.516, 1.697, 1.609, 0.760, 0.073], {"gene":'CD70'}
>>> data.domain.attributes[0]
Orange.feature.Continuous 'GSM63816'
>>> data.domain.attributes[0].attributes
{'dose': '20 U/ml IL-2', 'infection': 'acute ', 'time': '1 d'}

GDS classes
===========

.. autoclass:: orangecontrib.bio.geo.GDSInfo
   :members:

An example with :obj:`GDSInfo`::

    >>> import orangecontrib.bio.geo
    >>> info = orangecontrib.bio.geo.GDSInfo()
    >>> list(info.keys())[:5]
    >>> ['GDS2526', 'GDS2524', 'GDS2525', 'GDS2522', 'GDS1618']
    >>> info['GDS2526']['title']
    'c-MYC depletion effect on carcinoma cell lines'
    >>> info['GDS2526']['platform_organism']
    'Homo sapiens'

.. autoclass:: orangecontrib.bio.geo.GDS
   :members:


Examples
========

The following script prints out information about a specific data
set. It does not download the data set, just uses the (local) GEO data
sets information file (:download:`geo_gds1.py <code/geo_gds1.py>`).

.. literalinclude:: code/geo_gds1.py
   :lines: 6-

The output of this script is::

    ID:
    GDS10
    Features:
    39114
    Genes:
    29822
    Organism:
    Mus musculus
    PubMed ID:
    11827943
    Sample types:
      tissue (spleen, thymus)
      disease state (diabetic, diabetic-resistant, nondiabetic)
      strain (NOD, Idd3, Idd5, Idd3+Idd5, Idd9, B10.H2g7, B10.H2g7 Idd3)

    Description:
    Examination of spleen and thymus of type 1 diabetes nonobese diabetic
    (NOD) mouse, four NOD-derived diabetes-resistant congenic strains and
    two nondiabetic control strains.


Samples in GEO data sets belong to sample subsets, which in turn belong
to specific types.  The above GDS10 has three sample types, of which the
subsets for the tissue type are spleen and thymus. For supervised data
mining it would be useful to find out which data sets provide enough
samples for each label. It is (semantically) convenient to perform
classification within sample subsets of the same type. The following
script goes through all data sets and finds those with enough
samples within each of the subsets for a specific type. The function
``valid`` determines which subset types (if any) satisfy our criteria
(:download:`geo_gds5.py <code/geo_gds5.py>`).

.. literalinclude:: code/geo_gds5.py
   :lines: 8-

The requested number of samples, ``n=40``, seems to be a quite
a stringent criteria met - at the time of writing this -
by 40 data sets with 48 sample subsets. The output starts with::

    GDS1292
      tissue:raphe magnus/40, somatomotor cortex/43
    GDS1293
      tissue:raphe magnus/40, somatomotor cortex/41
    GDS1412
      protocol:no treatment/47, hormone replacement therapy/42
    GDS1490
      other:non-neural/50, neural/100
    GDS1611
      genotype/variation:wild type/48, upf1 null mutant/48
    GDS2373
      gender:male/82, female/48
    GDS2808
      protocol:training set/44, validation set/50

Let us now pick data set GDS2960 and see if we can predict the disease
state. We will use logistic regression, and within 10-fold cross
validation measure AUC, the area under ROC. AUC is the probability of
correctly distinguishing the two classes, (e.g., the disease and control). 
From (:download:`geo_gds6.py <code/geo_gds6.py>`):

.. literalinclude:: code/geo_gds6.py

The output of this script is::
    
    Samples: 101, Genes: 4068
    AUC = 0.983

The AUC for this data set is very high, indicating that using these gene
expression data it is almost trivial to separate the two classes.
.. py:currentmodule:: orangecontrib.bio.go
.. py:module:: orangecontrib.bio.go

=========================
Gene Ontology (:mod:`go`)
=========================


Provides access to `Gene Ontology`_ and its gene annotations.

.. _Gene Ontology: http://geneontology.org/


.. autoclass:: orangecontrib.bio.go.Ontology(filename=None, progress_callback=None, rev=None)
   :members:
      defined_slims_subsets, named_slims_subset, set_slims_subset, slims_for_term,
      extract_super_graph, extract_sub_graph, __getitem__, __len__,
      __iter__, __contains__

   Ontology supports a subset of the Mapping protocol:

   >>> term_ids = list(ontology)
   >>> term = ontology[term_ids[0]]


.. autoclass:: orangecontrib.bio.go.Term

   .. attribute:: id

   	  The term id.

   .. attribute:: namespace

      The namespace of the term.

   .. attribute:: def_

      The term definition (Note the use of trailing underscore
      to avoid conflict with a python keyword).

   .. attribute:: is_a

      List of term ids this term is a subterm of (parent terms).

   .. attribute:: related

      List of (rel_type, term_id) tuples with rel_type specifying
      the relationship type with term_id.


.. autoclass:: orangecontrib.bio.go.Annotations(filename_or_organism=None, ontology=None, genematcher=None, progress_callback=None, rev=None)
   :members:
   :member-order: bysource
   :exclude-members:
      set_ontology, load, Load, parse_file, ParseFile, RemapGenes,
      AddAnnotation, DrawEnrichmentGraph, GetAllAnnotations, GetAllGenes,
      GetAnnotatedTerms, GetEnrichedTerms, GetGeneNamesTranslator,
      GetOntology, SetOntology, aliasMapper, allAnnotations,
      geneAnnotations, geneNames, geneNamesDict, termAnnotations


.. autoclass:: orangecontrib.bio.go.AnnotationRecord
   :members:
       from_string,
   	   DB, DB_Object_ID, DB_Object_Symbol, Qualifier, GO_ID,
       DB_Reference, Evidence_Code, With_From, Aspect, DB_Object_Name,
       DB_Object_Synonym, DB_Object_Type, Taxon, Date, Assigned_By,
       Annotation_Extension, Gene_Product_Form_ID


Example
-------

Load the ontology and print out some terms::

	from orangecontrib.bio import go
	ontology = go.Ontology()
	term = ontology["GO:0097194"] # execution phase of apoptosis

	# print a term
	print(term)

	# access fields by name
	print(term.id, term.name)
	# note the use of underscore due to a conflict with a python def keyword
	print(term.def_)


Searching the annotation (part of :download:`code/go_gene_annotations.py`)

.. literalinclude:: code/go_gene_annotations.py

Term enrichment (part of :download:`code/go_enrichment.py`)

.. literalinclude:: code/go_enrichment.py
   :lines: 6-
============================================================
KEGG - Kyoto Encyclopedia of Genes and Genomes (:mod:`kegg`)
============================================================

.. py:currentmodule:: orangecontrib.bio.kegg

.. automodule:: orangecontrib.bio.kegg
   :members:
   :member-order: bysource

DBEntry (:mod:`entry`)
----------------------

The :class:`entry.DBEntry` represents a DBGET databas entry.
The individual KEGG Database interfaces below provide their own
specialization for this base class.
 
.. autoclass:: orangecontrib.bio.kegg.entry.DBEntry
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Databases interface (:mod:`databases`)
-------------------------------------------

.. autoclass:: orangecontrib.bio.kegg.databases.DBDataBase
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.GenomeEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.Genome
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.GeneEntry
   :members:
   :exclude-members:
      alt_names
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bio.kegg.databases.Genes
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.CompoundEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bio.kegg.databases.Compound
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.ReactionEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bio.kegg.databases.Reaction
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.EnzymeEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bio.kegg.databases.Enzyme
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.PathwayEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bio.kegg.databases.Pathway
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Pathway (:mod:`pathway`)
-----------------------------

.. autoclass:: orangecontrib.bio.kegg.pathway.Pathway
   :members:
   :exclude-members:
      entrys
   :member-order: bysource
   :show-inheritance:


Utilities
---------

.. autoclass:: orangecontrib.bio.kegg.entry.parser.DBGETEntryParser
   :members:
