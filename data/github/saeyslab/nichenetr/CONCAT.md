<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

# nichenetr

<!-- badges: start -->

[![R build
status](https://github.com/saeyslab/nichenetr/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/saeyslab/nichenetr/actions)
[![Coverage
Status](https://codecov.io/gh/saeyslab/nichenetr/branch/master/graph/badge.svg)](https://codecov.io/gh/saeyslab/nichenetr)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)
<!-- badges: end -->

**nichenetr: the R implementation of the NicheNet method.** The goal of
NicheNet is to study intercellular communication from a computational
perspective. NicheNet uses human or mouse gene expression data of
interacting cells as input and combines this with a prior model that
integrates existing knowledge on ligand-to-target signaling paths. This
allows to predict ligand-receptor interactions that might drive gene
expression changes in cells of interest.

We describe the NicheNet algorithm in the following paper: [NicheNet:
modeling intercellular communication by linking ligands to target
genes](https://www.nature.com/articles/s41592-019-0667-5).

Bonnardel, T’Jonck et al. already used NicheNet to predict upstream
niche signals driving Kupffer cell differentiation [Stellate Cells,
Hepatocytes, and Endothelial Cells Imprint the Kupffer Cell Identity on
Monocytes Colonizing the Liver Macrophage
Niche](https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1).

### Important update!

12-01-2022: In the Liver Atlas paper from Guilliams et al.: [Spatial
proteogenomics reveals distinct and evolutionarily conserved hepatic
macrophage
niches](https://www.sciencedirect.com/science/article/pii/S0092867421014811),
we used Differential NicheNet, an extension to the default NicheNet
algorithm. **Differential NicheNet** can be used to compare cell-cell
interactions between different niches and better predict niche-specific
ligand-receptor (L-R) pairs. It was used in that paper to predict
ligand-receptor pairs specific for the Kupffer cell niche in mouse and
human.

The main difference between the classic NicheNet pipeline and the
Differential NicheNet pipeline is that Differential NicheNet also uses
the differential expression between the conditions/niches of the
ligand-receptor pairs for prioritization in addition to the ligand
activities. The classic NicheNet pipeline on the contrary uses only
ligand acivity for prioritization (and shows differential expression
only in visualizations).

So if you have data of multiple conditions or niches, and you want to
include differential expression of the ligand-receptor pairs in the
prioritization, we recommend you check out Differential NicheNet. At the
bottom of this page, you can find the links to two vignettes
illustrating a Differential NicheNet analysis. We recommend these
vignettes if you want to apply Differential NicheNet on your own data.
If you want to see the code used for the analyses used in the Guilliams
et al. paper, see <https://github.com/saeyslab/NicheNet_LiverCellAtlas>.

## Introduction to NicheNet

The figure below shows a graphical representation of the NicheNet
workflow. Interactions inferred from several complementary
ligand-receptor, signaling and gene regulatory data sources were
aggregated in respective integrated networks from which ligand-target
regulatory potential scores were calculated. This model of prior
information on potential ligand-target links can then be used to infer
active ligand-target links between interacting cells. NicheNet
prioritizes ligands according to their activity (i.e., how well they
predict observed changes in gene expression in the receiver cell) and
looks for affected targets with high potential to be regulated by these
prioritized ligands.

We offer the option to use the prebuilt prior model (such that the
network integration steps should not be repeated), or to create and use
your own prior model (see reference to detailed vignette below).

<br><br> ![](vignettes/workflow_nichenet.jpg) <br><br>

NicheNet strongly differs from most current computational approaches to
study intercellular communication. Current approaches study
intercellular communication from (single-cell) expression data by
linking ligands expressed by sender cells to their corresponding
receptors expressed by receiver cells. However, functional understanding
of a cellular communication process also requires knowing how these
inferred ligand-receptor interactions result in changes in the
expression of downstream target genes within the receiver cells. To
address this need, we developed NicheNet. Contrary to existing
approaches, NicheNet looks at gene regulatory effects of ligands because
the used prior knowledge goes beyond ligand-receptor interactions and
incorporates intracellular signaling and transcriptional regulation as
well. As a result, NicheNet allows to predict which ligands influence
the expression in another cell, which target genes are affected by each
ligand and which signaling mediators may be involved. By generating
these novel types of hypotheses, NicheNet can drive an improved
functional understanding of a cell-cell communication process of
interest. The figure below summarizes the conceptual differences between
most current ligand-receptor network inference approaches (top panel)
and NicheNet (bottom panel) and visualizes the power of NicheNet in
prioritizing ligand-receptor interactions based on gene expression
effects.

<br><br>
<img src="vignettes/comparison_other_approaches_2.jpg" width="450" />
<br><br>

## Main functionalities of nichenetr

Specific functionalities of this package include:

-   assessing how well ligands expressed by a sender cell can predict
    changes in gene expression in the receiver cell
-   prioritizing ligands based on their effect on gene expression
-   inferring putative ligand-target links active in the system under
    study
-   inferring potential signaling paths between ligands and target genes
    of interest: to generate causal hypotheses and check which data
    sources support the predictions
-   validation of the prior ligand-target model
-   construction of user-defined prior ligand-target models

Moreover, we provide instructions on how to make intuitive
visualizations of the main predictions (e.g., via circos plots as shown
here below).

<br><br> ![](vignettes/circos_plot_adapted.jpg)

## Installation of nichenetr

Installation typically takes a few minutes, depending on the number of
dependencies that has already been installed on your pc. You can install
nichenetr (and required dependencies) from github with:

    # install.packages("devtools")
    devtools::install_github("saeyslab/nichenetr")

nichenetr was tested on both Windows and Linux (most recently tested R
version: R 4.0.0)

## Learning to use nichenetr

To learn using nichenetr, read one of the following vignettes explaining
several types of analyses:

Following vignette contains the explanation on how to perform a basic
NicheNet analysis. This includes prioritizing ligands and predicting
target genes of prioritized ligands. This demo analysis takes only a few
minutes to run:

-   [NicheNet’s ligand activity analysis on a gene set of interest:
    predict active ligands and their target
    genes](vignettes/ligand_activity_geneset.md):
    `vignette("ligand_activity_geneset", package="nichenetr")`

To facilitate the use of NicheNet on single-cell data, we demonstrate
the use of NicheNet on a Seurat object in following vignettes. One
demonstrates the use of a single wrapper function, the other
demonstrates what’s behind the wrapper (recommended).

-   [Perform NicheNet analysis starting from a Seurat
    object](vignettes/seurat_wrapper.md):`vignette("seurat_wrapper", package="nichenetr")`
-   [Perform NicheNet analysis starting from a Seurat object:
    step-by-step
    analysis](vignettes/seurat_steps.md):`vignette("seurat_steps", package="nichenetr")`

Following vignettes contain explanation on how to do some follow-up
analyses after performing the most basic analysis:

-   [Inferring ligand-to-target signaling
    paths](vignettes/ligand_target_signaling_path.md):
    `vignette("ligand_target_signaling_path", package="nichenetr")`
-   [Assess how well top-ranked ligands can predict a gene set of
    interest](vignettes/target_prediction_evaluation_geneset.md):
    `vignette("target_prediction_evaluation_geneset", package="nichenetr")`
-   [Single-cell NicheNet’s ligand activity
    analysis](vignettes/ligand_activity_single_cell.md):
    `vignette("ligand_activity_single_cell", package="nichenetr")`

If you want to make a circos plot visualization of the NicheNet output,
you can check following vignettes:

-   [Circos plot visualization to show active ligand-target links
    between interacting
    cells](vignettes/circos.md):`vignette("circos", package="nichenetr")`.
-   [Seurat Wrapper + Circos
    visualization](vignettes/seurat_wrapper_circos.md):`vignette("seurat_wrapper_circos", package="nichenetr")`.

People interested in building own models or benchmark own models against
NicheNet can read one of the following vignettes:

-   [Model construction](vignettes/model_construction.md):
    `vignette("model_construction", package="nichenetr")`
-   [Model evaluation: target gene and ligand activity
    prediction](vignettes/model_evaluation.md):
    `vignette("model_evaluation", package="nichenetr")`
-   [Parameter optimization via
    mlrMBO](vignettes/parameter_optimization.md):
    `vignette("parameter_optimization", package="nichenetr")`

People working with mouse data can see in the following vignette how to
convert NicheNet’s ligand-target model (given in human symbols) to mouse
symbols:

-   [Converting NicheNet’s model from human to mouse
    symbols](vignettes/symbol_conversion.md):
    `vignette("symbol_conversion", package="nichenetr")`

Differential NicheNet vignettes:

-   [Differential NicheNet analysis between niches of
    interest](vignettes/differential_nichenet.md):`vignette("differential_nichenet", package="nichenetr")`
-   [Differential NicheNet analysis between conditions of
    interest](vignettes/differential_nichenet_pEMT.md):`vignette("differential_nichenet_pEMT", package="nichenetr")`

## FAQ

Check the FAQ page at [FAQ NicheNet](vignettes/faq.md):
`vignette("faq", package="nichenetr")`

## References

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<doi:10.1038/s41592-019-0667-5>

Bonnardel et al. Stellate Cells, Hepatocytes, and Endothelial Cells
Imprint the Kupffer Cell Identity on Monocytes Colonizing the Liver
Macrophage Niche. Immunity (2019) <doi:10.1016/j.immuni.2019.08.017>

Guilliams et al. Spatial proteogenomics reveals distinct and
evolutionarily conserved hepatic macrophage niches. Cell (2022)
<doi:10.1016/j.cell.2021.12.018>
Construction of NicheNet’s ligand-target model
================
Robin Browaeys
2018-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/model_construction.Rmd", output_format = "github_document")
-->

This vignette shows how ligand-target prior regulatory potential scores
are inferred in the NicheNet framework. You can use the procedure shown
here to develop your own model with inclusion of context-specific
networks or removal of noisy irrelevant data sources. The networks at
the basis of NicheNet can be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758).

# Background information about NicheNet’s prior ligand-target model

The prior model at the basis of NicheNet denotes how strongly existing
knowledge supports that a ligand may regulate the expression of a target
gene. To calculate this ligand-target regulatory potential, we
integrated biological knowledge about ligand-to-target signaling paths
as follows.

First, we collected multiple complementary data sources covering
ligand-receptor, signal transduction (e.g., protein-protein and
kinase-substrate interactions) and gene regulatory interactions (e.g.,
inferred from ChIP-seq and motifs). For information of all collected
data sources (link to the website of the database, etc), see [Data
source information](data_sources.xlsx)

Secondly, we integrated these individual data sources into two weighted
networks: 1) a ligand-signaling network, which contains protein-protein
interactions covering the signaling paths from ligands to downstream
transcriptional regulators; and 2) a gene regulatory network, which
contains gene regulatory interactions between transcriptional regulators
and target genes. To let informative data sources contribute more to the
final model, we weighted each data source during integration. These data
source weights were automatically determined via model-based parameter
optimization to improve the accuracy of ligand-target predictions (see
the vignette [Parameter optimization via
mlrMBO](parameter_optimization.md). In this vignette, we will show how
to construct models with unoptimized data source weigths as well.

Finally, we combined the ligand-signaling and gene regulatory network to
calculate a regulatory potential score between all pairs of ligands and
target genes. A ligand-target pair receives a high regulatory potential
if the regulators of the target gene are lying downstream of the
signaling network of the ligand. To calculate this, we used network
propagation methods on the integrated networks to propagate the signal
starting from a ligand, flowing through receptors, signaling proteins,
transcriptional regulators, and ultimately ending at target genes.

A graphical summary of this procedure is visualized here
below:

![](workflow_model_construction.png)

# Construct a ligand-target model from all collected ligand-receptor, signaling and gene regulatory network data sources

Load the required packages and networks we will use to construct the
model.

``` r
library(nichenetr)
library(tidyverse)

# in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions

# The complete networks can be downloaded from Zenodo
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
```

## Construct NicheNet’s ligand-target model from unoptimized data source weights

Construct the weighted integrated ligand-signaling and gene regulatory
network. In this first example, we give every data source the same
weight (as given by the `source_weights_df` data frame provided by
default by the nichenetr package). See the vignette showing how to use
mlrMBO to optimize data source weights and the hyperparameters if
interested in performing parameter optimization. For the hyperparameters
of the model (hub correction factors and damping factor), we will use
the optimized values (as given by the `hyperparameter_list` data frame
provided by default by the nichenetr package).

The ligand-signaling network hub correction factor and gene regulatory
network hub correction factor were defined as hyperparameter of the
model to mitigate the potential negative influence of over-dominant hubs
on the final model. The damping factor hyperparameter is the main
parameter of the Personalized PageRank algorithm, which we used as
network propagation algorithm to link ligands to downstream
regulators.

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df = source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)
```

Infer ligand-target regulatory potential scores based on the weighted
integrated
networks

``` r
# in this example we will calculate target gene regulatory potential scores for TNF and the ligand combination TNF+IL6
ligands = list("TNF",c("TNF","IL6"))
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)
```

Show some top target genes of the ligand TNF and the ligand combination
TNF+IL6

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##      HACD4       P3H2        UBD       IGHD      CCL19       SELE 
## 0.05490710 0.03997878 0.02700469 0.02547128 0.02524468 0.02503142 
##     MUC5AC       COX1        CRP      CXCL9 
## 0.02383718 0.02357808 0.02357626 0.02337854
```

``` r
extract_top_n_targets("TNF-IL6",10,ligand_target_matrix)
##      BUD23      HACD4       IGHD        CRP       COX1       P3H2 
## 0.03103603 0.02749870 0.02550983 0.02412252 0.02380473 0.02014130 
##       SELE       IL11      CASP3       MMP9 
## 0.01741225 0.01713076 0.01676342 0.01634211
```

## Construct NicheNet’s ligand-target model from optimized data source weights

Now, we will demonstrate how you can make an alternative model with the
optimized data source weights (as given by the
`optimized_source_weights_df` data frame provided by default by the
nichenetr
package)

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network,source_weights_df = optimized_source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##       HACD4        P3H2        SELE       VCAM1         UBD        CD1A 
## 0.027721165 0.020233320 0.012049496 0.009775023 0.009667994 0.009532865 
##       CCL19      MUC5AC       CXCL9         CRP 
## 0.008979892 0.008507522 0.008459968 0.008431759
```

# Change the data sources at the basis of the NicheNet ligand-target model

### Keep only specific data sources of interest

Now, we will demonstrate how you can decide which data sources to use in
the model you want to create. Let’s say for this example, that you are
interested in making a model that only consists of literature-derived
ligand-receptor interactions, signaling and gene regulatory interactions
from comprehensive databases and gene regulatory interactions inferred
from ChIP-seq. An annotation of the different data sources is given by
the `annotation_data_sources` data frame provided by default by the
nichenetr package)

``` r
annotation_data_sources$type_db %>% unique()
## [1] "literature"       "prediction"       "comprehensive_db"
## [4] "ptm"              "text_mining"      "directional_ppi" 
## [7] "ChIP"             "motif"            "perturbation"
```

``` r
data_sources_to_keep = annotation_data_sources %>% filter(type_db %in% c("literature","comprehensive_db","ChIP")) %>% pull(source)

new_source_weights_df = source_weights_df %>% filter(source %in% data_sources_to_keep)
new_lr_network = lr_network %>% filter(source %in% data_sources_to_keep) 
new_sig_network = sig_network %>% filter(source %in% data_sources_to_keep)
new_gr_network = gr_network %>% filter(source %in% data_sources_to_keep)
```

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = new_lr_network, sig_network = new_sig_network, gr_network = new_gr_network, source_weights_df = new_source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##      HACD4       P3H2       CD1A       SELE        UBD      CCL19 
## 0.05482765 0.03985222 0.03023781 0.02817257 0.02486998 0.02239381 
##       CCL8        STS       NOX4      KRT12 
## 0.02206470 0.02137445 0.02072580 0.02034479
```

To give a second example: say that you don’t trust TF-target links
inferred from motif information and want to construct a model with all
data sources except the motif
ones.

``` r
data_sources_to_remove = annotation_data_sources %>% filter(type_db %in% c("motif")) %>% pull(source)
data_sources_to_keep = annotation_data_sources$source %>% setdiff(data_sources_to_remove) 

new_source_weights_df = source_weights_df %>% filter(source %in% data_sources_to_keep)
new_lr_network = lr_network %>% filter(source %in% data_sources_to_keep) 
new_sig_network = sig_network %>% filter(source %in% data_sources_to_keep)
new_gr_network = gr_network %>% filter(source %in% data_sources_to_keep)
```

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = new_lr_network, sig_network = new_sig_network, gr_network = new_gr_network, source_weights_df = new_source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##      HACD4       P3H2      CCL19        UBD       IGHD    PLA2G2A 
## 0.05490710 0.03997878 0.03132835 0.03120573 0.02856990 0.02650077 
##       SELE     CXCL11        CRP     MUC5AC 
## 0.02612114 0.02549863 0.02547906 0.02531375
```

### Add own data sources to the NicheNet model

In addition to removing data sources, you can also add new data sources.
This could for example help you in making context-specific models, if
you would have a network or data containing context-specific
interactions of interest.

As input, we require a data source to contain directional interactions
between genes: these interactions are protein-protein or signaling
interactions for ligand-receptor and signaling data sources and a gene
regulatory interaction for gene regulatory data sources. The data
sources should be formatted in a data frame with following columns:
from, to and source. “from” denotes the source node “gene A” of the
directional interaction from gene A to B, “to” denotes the target node
“gene B” of this directional interaction, and “source” is a
user-defined name of this data source.

Here, we will show how you can download, process and integrate an online
data source within the NichenNet framework. As example, this is the data
source “Hub Proteins Protein-Protein Interactions” from the Harmonizome
portal
(<https://amp.pharm.mssm.edu/Harmonizome/dataset/Hub+Proteins+Protein-Protein+Interactions>).

``` r
input_file = "https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/hubs/gene_attribute_edges.txt.gz"
ppi_network = read_tsv(input_file, col_names = TRUE)

ppi_network = ppi_network %>% transmute(from=target,to=source) %>% 
    filter(from %in% geneinfo_human$symbol & to %in% geneinfo_human$symbol) # keep only interactions between genes with oficial gene symbols: optional step

# give your data source a name
ppi_network = ppi_network %>% mutate(source = "harmonizome_hub_ppi", database = "harmonizome") 

head(ppi_network)
## # A tibble: 6 x 4
##   from  to    source              database   
##   <chr> <chr> <chr>               <chr>      
## 1 LRIF1 GAD1  harmonizome_hub_ppi harmonizome
## 2 LRIF1 ID2   harmonizome_hub_ppi harmonizome
## 3 LRIF1 NOC2L harmonizome_hub_ppi harmonizome
## 4 LRIF1 ESR1  harmonizome_hub_ppi harmonizome
## 5 LRIF1 NR3C1 harmonizome_hub_ppi harmonizome
## 6 LRIF1 PPARG harmonizome_hub_ppi harmonizome
```

First, we will add this new data source to all other data sources.
Because this data sources contains intracellular protein-protein
interactions, we will consider this data source as a signaling data
source. As example, we will assign to this data source a weight of 1,
because we want it to have a strong contribution to the final model.

``` r
new_sig_network = sig_network %>% bind_rows(ppi_network)

new_network_weights_df = tibble(source = "harmonizome_hub_ppi", weight = 1)
new_source_weights_df = optimized_source_weights_df %>% bind_rows(new_network_weights_df)
```

Now make this
model

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = new_sig_network, gr_network = gr_network, source_weights_df = new_source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##       HACD4        P3H2        SELE       VCAM1         UBD        CD1A 
## 0.027723418 0.020234958 0.012049227 0.009774414 0.009667525 0.009533070 
##       CCL19      MUC5AC       CXCL9         CRP 
## 0.008979282 0.008507681 0.008458946 0.008430639
```

In some cases, it’s possible that you want that your data source will be
considered as only data source in a specific layer
(i.e. ligand-receptor, signaling or gene regulatory layer). Therefore,
we will show how you use this new data source as only data source part
of the signaling network.

``` r
new_sig_network = ppi_network

new_network_weights_df = tibble(source = "harmonizome_hub_ppi", weight = 1)
new_source_weights_df = optimized_source_weights_df %>% bind_rows(new_network_weights_df)
```

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = new_sig_network, gr_network = gr_network, source_weights_df = new_source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##       HACD4        P3H2        SELE        CD1A         UBD       VCAM1 
## 0.027708362 0.020048863 0.011728851 0.009356392 0.009233937 0.009211393 
##       CCL19      MUC5AC        CCL8         CRP 
## 0.008644264 0.008246319 0.008010251 0.007858136
```

## Final note

Most optimally, you would like to optimize the parameters again when
including own data sources. Instructions to do this are given in the
following vignette: [Parameter optimization via
mlrMBO](parameter_optimization.md): `vignette("parameter_optimization",
package="nichenetr")`

However, this optimization process takes a lot of time and requires the
availability of multiple cores to perform the optimization in parallel.
Because we demonstrate in the NicheNet paper that unoptimized models
also perform considerably well, data source weight optmization is not
necessary to have decent predictive ability.
Single-cell NicheNet’s ligand activity analysis
================
Robin Browaeys
2018-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/ligand_activity_single_cell.Rmd", output_format = "github_document")
-->

This vignette shows how NicheNet can be used to predict which ligands
might be active in single-cells. If a ligand has a high activity in a
cell, this means that target genes of that ligand are stronger expressed
in that cell than in other cells. In this example, we will use data from
Puram et al. to explore intercellular communication in the tumor
microenvironment in head and neck squamous cell carcinoma (HNSCC) (See
Puram et al. 2017). More specifically, we will assess the activity of
cancer-associated fibroblast (CAF) ligands in malignant cells. The used
ligand-target matrix and example expression data of interacting cells
can be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758).

In order to prioritize ligands regulating a process of interest, you can
perform a regression/correlation analysis between ligand activities in
cells, and scores of a cell corresponding to the process of interest.
For example, in this case study we were interested in finding ligands
regulating p-EMT. Therefore we correlated ligand activities to the p-EMT
scores of cells.

The purpose of this single-cell ligand activity analysis is to offer a
complementary way to prioritize ligands driving the process of interest
and to better analyze heterogeneity in ligand activity between different
cells.

### Load nichenetr and tidyverse

``` r
library(nichenetr)
library(tidyverse)
```

### Read in expression data of interacting cells

First, we will read in the single-cell data from CAF and malignant cells
from HNSCC tumors (See Puram et al.
2017).

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

Secondly, we will determine which genes are expressed in CAFs and
malignant cells from high quality primary tumors. Therefore, we wil not
consider cells from tumor samples of less quality or from lymph node
metastases. To determine expressed genes, we use the definition used by
of Puram et
al.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter((tumor %in% tumors_remove == FALSE)) %>% filter(`non-cancer cell type` == "CAF") %>% .$cell
malignant_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1) %>% filter((tumor %in% tumors_remove == FALSE)) %>% .$cell

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
```

### Load the ligand-target model we want to use

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

### Perform NicheNet’s single-cell ligand activity analysis

In a first step, we will define a set of potentially active ligands. As
potentially active ligands, we will use ligands that are 1) expressed by
CAFs and 2) can bind a (putative) receptor expressed by malignant cells.
Putative ligand-receptor links were gathered from NicheNet’s
ligand-receptor data
sources.

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network$from %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_CAFs)
receptors = lr_network$to %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% .$from %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"
```

In a second step, we will scale the single-cell expression data
(including only expressed
genes).

``` r
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
expression_scaled = expression %>% .[malignant_ids,background_expressed_genes] %>% scale_quantile()
```

Now perform the ligand activity analysis: infer how well NicheNet’s
ligand-target potential scores can predict whether a gene belongs to
most strongly expressed genes in a cell compared to other cells. To
reduce the running time for this vignette, we will perform the analysis
only on 10 example cells from the HN5 tumor. This vignette’s only
purpose is to illustrate the analysis.

In practice, ligand activity analysis for several cells can be better
run in parallel (via
e.g. parallel::mclapply)\!

``` r
malignant_hn5_ids = sample_info %>% filter(tumor == "HN5") %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1)  %>% .$cell %>% head(10)

ligand_activities = predict_single_cell_ligand_activities(cell_ids = malignant_hn5_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

### Ligand prioritization by regression analysis

Furthermore, we will also show how you can perform additional analyses
by linking the ligand activity in cells to other properties of cells in
order to prioritize ligands. As toy example, we will score malignant
cells here on the extent to which they express the core p-EMT gene
“TGFBI”.

``` r
cell_scores_tbl = tibble(cell = malignant_hn5_ids, score = expression_scaled[malignant_hn5_ids,"TGFBI"])
```

Then, we will determine the correlation between these p-EMT scores and
ligand activities over all cells to prioritize p-EMT-inducing ligands.
We hypothesize that ligands might be potential regulators of the p-EMT
program if higher ligand activities are associated with higher p-EMT
scores. Based on this correlation, we obtained a ranking of potential
p-EMT-inducing ligands.

To do so, we frist need to process and normalize the ligand activities
(i.e. pearson correlation values) to make different cells comparable.
Here we use modified z-score
normalization.

``` r
normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)
```

Then, we combine the ligand activities and cell property scores and
perform correlation and regression analysis. We can prioritize ligands
by ranking them based on the pearson correlation between activity scores
and property
scores.

``` r
output_correlation_analysis = single_ligand_activity_score_regression(normalized_ligand_activities,cell_scores_tbl)
output_correlation_analysis %>% arrange(-pearson_regression) %>% select(pearson_regression, ligand)
## # A tibble: 131 x 2
##    pearson_regression ligand  
##                 <dbl> <chr>   
##  1              0.525 TNC     
##  2              0.497 TFPI    
##  3              0.491 SEMA5A  
##  4              0.488 ANXA1   
##  5              0.473 TNFSF13B
##  6              0.462 IBSP    
##  7              0.449 HDGF    
##  8              0.443 HSP90B1 
##  9              0.431 CALM3   
## 10              0.428 CXCL11  
## # ... with 121 more rows
```

Visualize the relation between ligand activity
and the cell's property score of interest

``` r
inner_join(cell_scores_tbl,normalized_ligand_activities) %>% ggplot(aes(score,TNC)) + geom_point() + geom_smooth(method = "lm")
```

![](ligand_activity_single_cell_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

### References

<div id="refs" class="references">

<div id="ref-puram_single-cell_2017">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
NicheNet’s ligand activity analysis on a gene set of interest: predict
active ligands and their target genes
================
Robin Browaeys
2019-01-17

<!-- github markdown built using 
rmarkdown::render("vignettes/ligand_activity_geneset.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform a basic NicheNet
analysis. A NicheNet analysis can help you to generate hypotheses about
an intercellular communication process of interest for which you have
bulk or single-cell gene expression data. Specifically, NicheNet can
predict 1) which ligands from one cell population (“sender/niche”) are
most likely to affect target gene expression in an interacting cell
population (“receiver/target”) and 2) which specific target genes are
affected by which of these predicted ligands.

Because NicheNet studies how ligands affect gene expression in
neighboring cells, you need to have data about this effect in gene
expression you want to study. So, you need to have a clear set of genes
that are putatively affected by ligands from one of more interacting
cells.

The pipeline of a basic NicheNet analysis consist mainly of the
following steps:

  - 1.  Define a “sender/niche” cell population and a “receiver/target”
        cell population present in your expression data and determine
        which genes are expressed in both populations

  - 2.  Define a gene set of interest: these are the genes in the
        “receiver/target” cell population that are potentially
        affected by ligands expressed by interacting cells (e.g. genes
        differentially expressed upon cell-cell interaction)

  - 3.  Define a set of potential ligands: these are ligands that are
        expressed by the “sender/niche” cell population and bind a
        (putative) receptor expressed by the “receiver/target”
        population

  - 4)  Perform NicheNet ligand activity analysis: rank the potential
        ligands based on the presence of their target genes in the gene
        set of interest (compared to the background set of genes)

  - 5)  Infer top-predicted target genes of ligands that are top-ranked
        in the ligand activity analysis

This vignette guides you in detail through all these steps. As example
expression data of interacting cells, we will use data from Puram et
al. to explore intercellular communication in the tumor
microenvironment in head and neck squamous cell carcinoma (HNSCC) (See
Puram et al. 2017). More specifically, we will look at which ligands
expressed by cancer-associated fibroblasts (CAFs) can induce a specific
gene program in neighboring malignant cells. This program, a partial
epithelial-mesenschymal transition (p-EMT) program, could be linked to
metastasis by Puram et al. 

The used ligand-target matrix and example expression data of interacting
cells can be downloaded from Zenodo.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)

## Step 0: Load required packages, NicheNet’s ligand-target prior model and processed expression data of interacting cells

Packages:

``` r
library(nichenetr)
library(tidyverse)
```

Ligand-target model:

This model denotes the prior potential that a particular ligand might
regulate the expression of a specific target gene.

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

Expression data of interacting cells: publicly available single-cell
data from CAF and malignant cells from HNSCC tumors:

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

## Step 1: Define expressed genes in sender and receiver cell populations

Our research question is to prioritize which ligands expressed by CAFs
can induce p-EMT in neighboring malignant cells. Therefore, CAFs are the
sender cells in this example and malignant cells are the receiver cells.
This is an example of paracrine signaling. Note that autocrine signaling
can be considered if sender and receiver cell type are the same.

Now, we will determine which genes are expressed in the sender cells
(CAFs) and receiver cells (malignant cells) from high quality primary
tumors. Therefore, we wil not consider cells from tumor samples of less
quality or from lymph node metastases.

To determine expressed genes in this case study, we use the definition
used by Puram et al. (the authors of this dataset), which is: Ea, the
aggregate expression of each gene i across the k cells, calculated as
Ea(i) = log2(average(TPM(i)1…k)+1), should be \>= 4. We recommend users
to define expressed genes in the way that they consider to be most
appropriate for their dataset. For single-cell data generated by the 10x
platform in our lab, we don’t use the definition used here, but we
consider genes to be expressed in a cell type when they have non-zero
values in at least 10% of the cells from that cell type. This is
described as well in the other vignette [Perform NicheNet analysis
starting from a Seurat object: step-by-step
analysis](seurat_steps.md):`vignette("seurat_steps",
package="nichenetr")`.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% pull(cell)
malignant_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)

expressed_genes_sender = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_receiver = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()

# Check the number of expressed genes: should be a 'reasonable' number of total expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000)
length(expressed_genes_sender)
## [1] 6706
length(expressed_genes_receiver)
## [1] 6351
```

## Step 2: Define the gene set of interest and a background of genes

As gene set of interest, we consider the genes of which the expression
is possibly affected due to communication with other cells. The
definition of this gene set depends on your research question and is a
crucial step in the use of NicheNet.

Because we here want to investigate how CAFs regulate the expression of
p-EMT genes in malignant cells, we will use the p-EMT gene set defined
by Puram et al. as gene set of interest and use all genes expressed in
malignant cells as background of genes.

``` r
geneset_oi = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(geneset_oi)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"
```

## Step 3: Define a set of potential ligands

As potentially active ligands, we will use ligands that are 1) expressed
by CAFs and 2) can bind a (putative) receptor expressed by malignant
cells. Putative ligand-receptor links were gathered from NicheNet’s
ligand-receptor data sources.

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

# If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions and only keep ligand-receptor interactions that are described in curated databases. To do this: uncomment following line of code:
# lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
## # A tibble: 6 x 4
##   from    to        source         database
##   <chr>   <chr>     <chr>          <chr>   
## 1 HGF     MET       kegg_cytokines kegg    
## 2 TNFSF10 TNFRSF10A kegg_cytokines kegg    
## 3 TNFSF10 TNFRSF10B kegg_cytokines kegg    
## 4 TGFB2   TGFBR1    kegg_cytokines kegg    
## 5 TGFB3   TGFBR1    kegg_cytokines kegg    
## 6 INHBA   ACVR2A    kegg_cytokines kegg
```

This ligand-receptor network contains the expressed ligand-receptor
interactions. As potentially active ligands for the NicheNet analysis,
we will consider the ligands from this network.

``` r
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"
```

## Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest

Now perform the ligand activity analysis: in this analysis, we will
calculate the ligand activity of each ligand, or in other words, we will
assess how well each CAF-ligand can predict the p-EMT gene set compared
to the background of expressed genes (predict whether a gene belongs to
the p-EMT program or not).

``` r
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

Now, we want to rank the ligands based on their ligand activity. In our
validation study, we showed that the pearson correlation coefficient
(PCC) between a ligand’s target predictions and the observed
transcriptional response was the most informative measure to define
ligand activity. Therefore, we will rank the ligands based on their
pearson correlation coefficient. This allows us to prioritize
p-EMT-regulating ligands.

``` r
ligand_activities %>% arrange(-pearson) 
## # A tibble: 131 x 4
##    test_ligand auroc   aupr pearson
##    <chr>       <dbl>  <dbl>   <dbl>
##  1 PTHLH       0.667 0.0720   0.128
##  2 CXCL12      0.680 0.0507   0.123
##  3 AGT         0.676 0.0581   0.120
##  4 TGFB3       0.689 0.0454   0.117
##  5 IL6         0.693 0.0510   0.115
##  6 INHBA       0.695 0.0502   0.113
##  7 ADAM17      0.672 0.0526   0.113
##  8 TNC         0.700 0.0444   0.109
##  9 CTGF        0.680 0.0473   0.108
## 10 FN1         0.679 0.0505   0.108
## # ... with 121 more rows
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)
## [1] "PTHLH"  "CXCL12" "AGT"    "TGFB3"  "IL6"    "INHBA"
```

We see here that the performance metrics indicate that the 20 top-ranked
ligands can predict the p-EMT genes reasonably, this implies that
ranking of the ligands might be accurate as shown in our study. However,
it is possible that for some gene sets, the target gene prediction
performance of the top-ranked ligands would not be much better than
random prediction. In that case, prioritization of ligands will be less
trustworthy.

Additional note: we looked at the top 20 ligands here and will continue
the analysis by inferring p-EMT target genes of these 20 ligands.
However, the choice of looking only at the 20 top-ranked ligands for
further biological interpretation is based on biological intuition and
is quite arbitrary. Therefore, users can decide to continue the analysis
with a different number of ligands. We recommend to check the selected
cutoff by looking at the distribution of the ligand activity values.
Here, we show the ligand activity histogram (the score for the 20th
ligand is indicated via the dashed line).

``` r
# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
```

![](ligand_activity_geneset_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap

Now we will show how you can look at the regulatory potential scores
between ligands and target genes of interest. In this case, we will look
at links between top-ranked p-EMT regulating ligands and p-EMT genes. In
the ligand-target heatmaps, we show here regulatory potential scores for
interactions between the 20 top-ranked ligands and following target
genes: genes that belong to the gene set of interest and to the 250 most
strongly predicted targets of at least one of the 20 top-ranked ligands
(the top 250 targets according to the general prior model, so not the
top 250 targets for this dataset). Consequently, genes of your gene set
that are not a top target gene of one of the prioritized ligands, will
not be shown on the heatmap.

``` r
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)
## # A tibble: 6 x 3
##   ligand target  weight
##   <chr>  <chr>    <dbl>
## 1 PTHLH  COL1A1 0.00399
## 2 PTHLH  MMP1   0.00425
## 3 PTHLH  MMP2   0.00210
## 4 PTHLH  MYH9   0.00116
## 5 PTHLH  P4HA2  0.00190
## 6 PTHLH  PLAU   0.00401
```

For visualization purposes, we adapted the ligand-target regulatory
potential matrix as follows. Regulatory potential scores were set as 0
if their score was below a predefined threshold, which was here the 0.25
quantile of scores of interactions between the 20 top-ranked ligands and
each of their respective top targets (see the ligand-target network
defined in the data frame).

``` r
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)
## # A tibble: 6 x 3
##   ligand target  weight
##   <chr>  <chr>    <dbl>
## 1 PTHLH  COL1A1 0.00399
## 2 PTHLH  MMP1   0.00425
## 3 PTHLH  MMP2   0.00210
## 4 PTHLH  MYH9   0.00116
## 5 PTHLH  P4HA2  0.00190
## 6 PTHLH  PLAU   0.00401
```

The putatively active ligand-target links will now be visualized in a
heatmap. The order of the ligands accord to the ranking according to the
ligand activity prediction.

``` r
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network
```

![](ligand_activity_geneset_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Note that the choice of these cutoffs for visualization is quite
arbitrary. We recommend users to test several cutoff values.

If you would consider more than the top 250 targets based on prior
information, you will infer more, but less confident, ligand-target
links; by considering less than 250 targets, you will be more stringent.

If you would change the quantile cutoff that is used to set scores to 0
(for visualization purposes), lowering this cutoff will result in a more
dense heatmap, whereas highering this cutoff will result in a more
sparse heatmap.

## Follow-up analysis 1: Ligand-receptor network inference for top-ranked ligands

One type of follow-up analysis is looking at which receptors of the
receiver cell population (here: malignant cells) can potentially bind to
the prioritized ligands from the sender cell population (here: CAFs).

So, we will now infer the predicted ligand-receptor interactions of the
top-ranked ligands and visualize these in a heatmap.

``` r
# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
```

Show a heatmap of the ligand-receptor interactions

``` r
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized CAF-ligands","Receptors expressed by malignant cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network
```

![](ligand_activity_geneset_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Follow-up analysis 2: Visualize expression of top-predicted ligands and their target genes in a combined heatmap

NicheNet only considers expressed ligands of sender cells, but does not
take into account their expression for ranking the ligands. The ranking
is purely based on the potential that a ligand might regulate the gene
set of interest, given prior knowledge. Because it is also useful to
further look into expression of ligands and their target genes, we
demonstrate here how you could make a combined figure showing ligand
activity, ligand expression, target gene expression and ligand-target
regulatory potential.

#### Load additional packages required for the visualization:

``` r
library(RColorBrewer)
library(cowplot)
library(ggpubr)
```

#### Prepare the ligand activity matrix

``` r
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
```

``` r
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
p_ligand_pearson
```

![](ligand_activity_geneset_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

#### Prepare expression of ligands in fibroblast per tumor

Because the single-cell data was collected from multiple tumors, we will
show here the average expression of the ligands per tumor.

``` r
expression_df_CAF = expression[CAF_ids,order_ligands] %>% data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% inner_join(sample_info %>% select(cell,tumor), by =  "cell")

aggregated_expression_CAF = expression_df_CAF %>% group_by(tumor) %>% select(-cell) %>% summarise_all(mean)

aggregated_expression_df_CAF = aggregated_expression_CAF %>% select(-tumor) %>% t() %>% magrittr::set_colnames(aggregated_expression_CAF$tumor) %>% data.frame() %>% rownames_to_column("ligand") %>% as_tibble() 

aggregated_expression_matrix_CAF = aggregated_expression_df_CAF %>% select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_CAF$ligand)

order_tumors = c("HN6","HN20","HN26","HN28","HN22","HN25","HN5","HN18","HN17","HN16") # this order was determined based on the paper from Puram et al. Tumors are ordered according to p-EMT score.
vis_ligand_tumor_expression = aggregated_expression_matrix_CAF[order_ligands,order_tumors]
```

``` r
library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
p_ligand_tumor_expression = vis_ligand_tumor_expression %>% make_heatmap_ggplot("Prioritized CAF-ligands","Tumor", color = color[100],legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(averaged over\nsingle cells)") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_tumor_expression
```

![](ligand_activity_geneset_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

#### Prepare expression of target genes in malignant cells per tumor

``` r
expression_df_target = expression[malignant_ids,geneset_oi] %>% data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% inner_join(sample_info %>% select(cell,tumor), by =  "cell") 

aggregated_expression_target = expression_df_target %>% group_by(tumor) %>% select(-cell) %>% summarise_all(mean)

aggregated_expression_df_target = aggregated_expression_target %>% select(-tumor) %>% t() %>% magrittr::set_colnames(aggregated_expression_target$tumor) %>% data.frame() %>% rownames_to_column("target") %>% as_tibble() 

aggregated_expression_matrix_target = aggregated_expression_df_target %>% select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)

vis_target_tumor_expression_scaled = aggregated_expression_matrix_target %>% t() %>% scale_quantile() %>% .[order_tumors,order_targets]
```

``` r
p_target_tumor_scaled_expression = vis_target_tumor_expression_scaled  %>% make_threecolor_heatmap_ggplot("Tumor","Target", low_color = color[1],mid_color = color[50], mid = 0.5, high_color = color[100], legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
p_target_tumor_scaled_expression
```

![](ligand_activity_geneset_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

#### Combine the different heatmaps in one overview figure

``` r
figures_without_legend = plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  p_ligand_tumor_expression + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""), 
  NULL,
  NULL,
  p_target_tumor_scaled_expression + theme(legend.position = "none", axis.ticks = element_blank()) + xlab(""), 
  align = "hv",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_pearson)+ 4.5, ncol(vis_ligand_tumor_expression), ncol(vis_ligand_target)) -2,
  rel_heights = c(nrow(vis_ligand_pearson), nrow(vis_target_tumor_expression_scaled) + 3)) 

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_pearson)),
  as_ggplot(get_legend(p_ligand_tumor_expression)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_target_tumor_scaled_expression)),
  nrow = 2,
  align = "h")

plot_grid(figures_without_legend, 
          legends, 
          rel_heights = c(10,2), nrow = 2, align = "hv")
```

![](ligand_activity_geneset_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## Other follow-up analyses:

As another follow-up analysis, you can infer possible signaling paths
between ligands and targets of interest. You can read how to do this in
the following vignette [Inferring ligand-to-target signaling
paths](ligand_target_signaling_path.md):`vignette("ligand_target_signaling_path",
package="nichenetr")`.

Another follow-up analysis is getting a “tangible” measure of how well
top-ranked ligands predict the gene set of interest and assess which
genes of the gene set can be predicted well. You can read how to do this
in the following vignette [Assess how well top-ranked ligands can
predict a gene set of
interest](target_prediction_evaluation_geneset.md):`vignette("target_prediction_evaluation_geneset",
package="nichenetr")`.

In case you want to visualize ligand-target links between multiple
interacting cells, you can make an appealing circos plot as shown in
vignette [Circos plot visualization to show active ligand-target links
between interacting cells](circos.md):`vignette("circos",
package="nichenetr")`.

## References

<div id="refs" class="references">

<div id="ref-puram_single-cell_2017">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
Seurat Wrapper + Circos visualization
================
Robin Browaeys
18-1-2021

<!-- github markdown built using 
rmarkdown::render("vignettes/seurat_wrapper_circos.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform a basic NicheNet analysis
on a Seurat v3 object - and how to visualize the output in a circos
plot. This vignette demonstrates the same workflow as shown in [Perform
NicheNet analysis starting from a Seurat
object](seurat_wrapper.md):`vignette("seurat_wrapper", package="nichenetr")`,
but adds a circos plot visualization as shown in [Circos plot
visualization to show active ligand-target links between interacting
cells](circos.md):`vignette("circos", package="nichenetr")`. For more
detailed information about the NicheNet workflow, check those vignettes.
This vignette was made upon popular request to demonstrate how those two
vignettes can be combined into one analysis workflow. Note that we as
developers of NicheNet generally recommend a visualization of the output
by combining several heatmaps (ligand activity, ligand-target links,
ligand-receptor links, ligand expression, ligand LFC,…) over using a
circos plot visualization. Certainly for cases with many sender cell
types and ligands that are expressed by more than one sender cell type.
Because in those cases, the circos plot is much less informative and
could lead to wrong interpretation of the results.

As example expression data of interacting cells, we will use mouse
NICHE-seq data from Medaglia et al. to explore intercellular
communication in the T cell area in the inguinal lymph node before and
72 hours after lymphocytic choriomeningitis virus (LCMV) infection \[See
@medaglia\_spatial\_2017\]. We will NicheNet to explore immune cell
crosstalk in response to this LCMV infection.

In this dataset, differential expression is observed between CD8 T cells
in steady-state and CD8 T cells after LCMV infection. NicheNet can be
applied to look at how several immune cell populations in the lymph node
(i.e., monocytes, dendritic cells, NK cells, B cells, CD4 T cells) can
regulate and induce these observed gene expression changes. NicheNet
will specifically prioritize ligands from these immune cells and their
target genes that change in expression upon LCMV infection.

The used NicheNet networks, ligand-target matrix and example expression
data of interacting cells can be downloaded from Zenodo. The NicheNet
networks and ligand-target matrix at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)
and the Seurat object of the processed NICHE-seq single-cell data at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3531889.svg)](https://doi.org/10.5281/zenodo.3531889).

# Prepare NicheNet analysis

## Load required packages, read in the Seurat object with processed expression data of interacting cells and NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks.

### Load Packages:

``` r
library(nichenetr)
library(Seurat) # Please update to Seurat v4
library(tidyverse)
library(circlize)
```

If you would use and load other packages, we recommend to load these 3
packages after the others.

### Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
## # A tibble: 6 x 4
##   from  to    source         database
##   <chr> <chr> <chr>          <chr>   
## 1 CXCL1 CXCR2 kegg_cytokines kegg    
## 2 CXCL2 CXCR2 kegg_cytokines kegg    
## 3 CXCL3 CXCR2 kegg_cytokines kegg    
## 4 CXCL5 CXCR2 kegg_cytokines kegg    
## 5 PPBP  CXCR2 kegg_cytokines kegg    
## 6 CXCL6 CXCR2 kegg_cytokines kegg

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  ABCC6  0.422 
## 2 A1BG  ACE2   0.101 
## 3 A1BG  ADAM10 0.0970
## 4 A1BG  AGO1   0.0525
## 5 A1BG  AKT1   0.0855
## 6 A1BG  ANXA7  0.457
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  A2M    0.0294
## 2 AAAS  GFAP   0.0290
## 3 AADAC CYP3A4 0.0422
## 4 AADAC IRF8   0.0275
## 5 AATF  ATM    0.0330
## 6 AATF  ATR    0.0355
```

### Read in the expression data of interacting cells:

``` r
seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj@meta.data %>% head()
##         nGene nUMI orig.ident aggregate res.0.6 celltype nCount_RNA nFeature_RNA
## W380370   880 1611      LN_SS        SS       1    CD8 T       1607          876
## W380372   541  891      LN_SS        SS       0    CD4 T        885          536
## W380374   742 1229      LN_SS        SS       0    CD4 T       1223          737
## W380378   847 1546      LN_SS        SS       1    CD8 T       1537          838
## W380379   839 1606      LN_SS        SS       0    CD4 T       1603          836
## W380381   517  844      LN_SS        SS       0    CD4 T        840          513
```

Visualize which cell populations are present: CD4 T cells (including
regulatory T cells), CD8 T cells, B cells, NK cells, dendritic cells
(DCs) and inflammatory monocytes

``` r
seuratObj@meta.data$celltype %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application
## .
##     B CD4 T CD8 T    DC  Mono    NK  Treg 
##   382  2562  1645    18    90   131   199
DimPlot(seuratObj, reduction = "tsne")
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
seuratObj@meta.data$aggregate %>% table()
## .
## LCMV   SS 
## 3886 1141
DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Perform the NicheNet analysis

## NicheNet analysis on Seurat object: explain differential expression between two conditions

In this case study, the receiver cell population is the ‘CD8 T’ cell
population, whereas the sender cell populations are ‘CD4 T’, ‘Treg’,
‘Mono’, ‘NK’, ‘B’ and ‘DC’. The above described functions will consider
a gene to be expressed when it is expressed in at least a predefined
fraction of cells in one cluster (default: 10%).

The gene set of interest are the genes differentially expressed in CD8 T
cells after LCMV infection. The condition of interest is thus ‘LCMV’,
whereas the reference/steady-state condition is ‘SS’. The notion of
conditions can be extracted from the metadata column ‘aggregate’, the
method to calculate the differential expression is the standard Seurat
Wilcoxon test.

The number of top-ranked ligands that are further used to predict active
target genes and construct an active ligand-receptor network is 20 by
default.

To perform the NicheNet analysis with these specifications, run the
following:

``` r
# indicated cell types should be cell class identities
# check via: 
# seuratObj %>% Idents() %>% table()
sender_celltypes = c("CD4 T","Treg", "Mono", "NK", "B", "DC")
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "CD8 T", 
  condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
```

### Interpret the NicheNet analysis output

#### Ligand activity analysis results

A first thing NicheNet does, is prioritizing ligands based on predicted
ligand activity. To see the ranking of these ligands, run the following
command:

``` r
nichenet_output$ligand_activities
## # A tibble: 44 x 6
##    test_ligand auroc  aupr pearson  rank bona_fide_ligand
##    <chr>       <dbl> <dbl>   <dbl> <dbl> <lgl>           
##  1 Ebi3        0.638 0.234  0.197      1 FALSE           
##  2 Il15        0.582 0.163  0.0961     2 TRUE            
##  3 Crlf2       0.549 0.163  0.0758     3 FALSE           
##  4 App         0.499 0.141  0.0655     4 TRUE            
##  5 Tgfb1       0.494 0.140  0.0558     5 TRUE            
##  6 Ptprc       0.536 0.149  0.0554     6 TRUE            
##  7 H2-M3       0.525 0.157  0.0528     7 TRUE            
##  8 Icam1       0.543 0.142  0.0486     8 TRUE            
##  9 Cxcl10      0.531 0.141  0.0408     9 TRUE            
## 10 Adam17      0.517 0.137  0.0359    10 TRUE            
## # ... with 34 more rows
```

These ligands are expressed by one or more of the input sender cells. To
see which cell population expresses which of these top-ranked ligands,
you can run the following:

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

As you can see, most op the top-ranked ligands seem to be mainly
expressed by dendritic cells and monocytes.

It could also be interesting to see whether some of these ligands are
differentially expressed after LCMV infection.

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), split.by = "aggregate") + RotatedAxis()
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
VlnPlot(seuratObj, features = c("Il15", "Cxcl10","Cxcl16"), split.by = "aggregate", pt.size = 0, combine = FALSE)
## [[1]]
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

    ## 
    ## [[2]]

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

    ## 
    ## [[3]]

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

#### Inferred active ligand-target links

NicheNet also infers active target genes of these top-ranked ligands. To
see which top-ranked ligands are predicted to have regulated the
expression of which differentially expressed genes, you can run
following command for a heatmap visualization:

``` r
nichenet_output$ligand_target_heatmap
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Circos plots to visualize ligand-target and ligand-receptor interactions

This visualization groups the top predicted active ligands according to
the strongest expressing cell type. Therefore we need to determine per
cell type which ligands they express more strongly than the other cell
types.

### Calculate average ligand expression in sender cells

``` r
# avg_expression_ligands = AverageExpression(seuratObj %>% subset(subset = aggregate == "LCMV"),features = nichenet_output$top_ligands) # if want to look specifically in LCMV-only cells
avg_expression_ligands = AverageExpression(seuratObj, features = nichenet_output$top_ligands)
```

### Assign ligands to sender cells

To assign ligands to sender cell type, we can e.g. look for which sender
cell types show an expression that is higher than the average + SD.

``` r
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
  }) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)
## [1] "B"    "NK"   "Mono" "DC"
```

The top ligands seem to be most strongly expressed by B cells, NK cells,
monocytes and DCs. We will know also look at which ligands are common
across multiple cell types (= those that are specific to &gt; 1 cell
type, or those that were not assigned to a cell type in the previous
block of code)

Determine now which prioritized ligands are expressed by CAFs and or
endothelial cells

``` r
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = nichenet_output$top_ligands %>% setdiff(unique_ligands)

B_specific_ligands = sender_ligand_assignment$B %>% names() %>% setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment$NK %>% names() %>% setdiff(general_ligands)
Mono_specific_ligands = sender_ligand_assignment$Mono %>% names() %>% setdiff(general_ligands)
DC_specific_ligands = sender_ligand_assignment$DC %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("B-specific", times = B_specific_ligands %>% length()),
                  rep("NK-specific", times = NK_specific_ligands %>% length()),
                  rep("Mono-specific", times = Mono_specific_ligands %>% length()),
                  rep("DC-specific", times = DC_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(B_specific_ligands, NK_specific_ligands, Mono_specific_ligands, DC_specific_ligands, general_ligands))
```

### Define the ligand-target links of interest

To avoid making a circos plots with too many ligand-target links, we
will show only links with a weight higher than a predefined cutoff:
links belonging to the 40% of lowest scores were removed. Not that this
cutoffs and other cutoffs used for this visualization can be changed
according to the user’s needs.

``` r
active_ligand_target_links_df = nichenet_output$ligand_target_df %>% mutate(target_type = "LCMV-DE") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
  
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
```

Prepare the circos visualization: give each segment of ligands and
targets a specific color and order

``` r
grid_col_ligand =c("General" = "lawngreen",
            "NK-specific" = "royalblue",
            "B-specific" = "darkgreen",
            "Mono-specific" = "violet",
            "DC-specific" = "steelblue2")
grid_col_target =c(
            "LCMV-DE" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
```

Prepare the circos visualization: order ligands and targets

``` r
target_order = circos_links$target %>% unique()
ligand_order = c(Mono_specific_ligands, DC_specific_ligands, NK_specific_ligands,B_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)
```

Prepare the circos visualization: define the gaps between the different
segments

``` r
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mono-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "DC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "NK-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "LCMV-DE") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
  )
```

Render the circos plot (all links same transparancy). Only the widths of
the blocks that indicate each target gene is proportional the
ligand-target regulatory potential (\~prior knowledge supporting the
regulatory interaction).

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
circos.clear()
```

Render the circos plot (degree of transparancy determined by the
regulatory potential value of a ligand-target interaction)

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
circos.clear()
```

Save circos plot to an svg file

``` r
svg("ligand_target_circos.svg", width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()
## png 
##   2
```

### Visualize ligand-receptor interactions of the prioritized ligands in a circos plot

``` r
lr_network_top_df = nichenet_output$ligand_receptor_df %>% mutate(receptor_type = "LCMV_CD8T_receptor") %>% inner_join(ligand_type_indication_df)
```

``` r
grid_col_ligand =c("General" = "lawngreen",
            "NK-specific" = "royalblue",
            "B-specific" = "darkgreen",
            "Mono-specific" = "violet",
            "DC-specific" = "steelblue2")
grid_col_receptor =c(
            "LCMV_CD8T_receptor" = "darkred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
```

Prepare the circos visualization: order ligands and receptors

``` r
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(Mono_specific_ligands, DC_specific_ligands, NK_specific_ligands,B_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)
```

Prepare the circos visualization: define the gaps between the different
segments

``` r
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mono-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "DC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "NK-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_receptor,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "LCMV_CD8T_receptor") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
  )
```

Render the circos plot (all links same transparancy). Only the widths of
the blocks that indicate each receptor is proportional the
ligand-receptor interaction weight (\~prior knowledge supporting the
interaction).

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
circos.clear()
```

Render the circos plot (degree of transparancy determined by the prior
interaction weight of the ligand-receptor interaction - just as the
widths of the blocks indicating each receptor)

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
```

![](seurat_wrapper_circos_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
circos.clear()
```

Save circos plot to an svg file

``` r
svg("ligand_receptor_circos.svg", width = 15, height = 15)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()
dev.off()
## png 
##   2
```
Converting NicheNet’s model from human to mouse symbols
================
Robin Browaeys
2019-07-31

<!-- github markdown built using
rmarkdown::render("vignettes/symbol_conversion.Rmd", output_format = "github_document")
-->

In this vignette, we show how to convert NicheNet’s ligand-target matrix
model from human to mouse gene symbols. This is necessary if you want to
apply NicheNet to mouse expression data, because the NicheNet prior
information is given in human gene symbols (because most data sources at
the basis of NicheNet are based on human data). One-to-one orthologs
were gathered from NCBI HomoloGene and also from ENSEMBL via biomaRt.

### Load required packages

``` r
library(nichenetr)
library(tidyverse)
```

### Load NicheNet’s ligand-target model:

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
dim(ligand_target_matrix)
## [1] 25345   688
```

### Convert the ligand-target model from human to mouse symbols.

Because not all human genes have a mouse one-to-one ortholog, these
genes will be removed from the mouse
model.

``` r
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

dim(ligand_target_matrix)
## [1] 17330   644
```

Show the top 10 targets of TNF (in mouse
symbols):

``` r
top_targets = extract_top_n_targets("Tnf",10,ligand_target_matrix) %>% names()
top_targets
##  [1] "Hacd4"  "P3h2"   "Sele"   "Vcam1"  "Ubd"    "Ccl19"  "Muc5ac"
##  [8] "Cxcl9"  "Crp"    "Icam1"
```

If you want to convert mouse to human symbols, you can use:

``` r
top_targets %>% convert_mouse_to_human_symbols()
##    Hacd4     P3h2     Sele    Vcam1      Ubd    Ccl19   Muc5ac    Cxcl9 
##  "HACD4"   "P3H2"   "SELE"  "VCAM1"    "UBD"  "CCL19" "MUC5AC"  "CXCL9" 
##      Crp    Icam1 
##    "CRP"  "ICAM1"
```
Differential NicheNet analysis between niches of interest
================
Robin Browaeys
20212-01-12

<!-- github markdown built using 
rmarkdown::render("vignettes/differential_nichenet.Rmd", output_format = "github_document")
-->

This vignette guides you in detail through all the steps of a
Differential NicheNet analysis. As example expression data of
interacting cells, we will here use subset of the liver scRNAseq data
generated in the paper from Guilliams et al: [Spatial proteogenomics
reveals distinct and evolutionarily conserved hepatic macrophage
niches](https://www.sciencedirect.com/science/article/pii/S0092867421014811).
We took a subset of the data (deposited on
<https://zenodo.org/deposit/5840787>) for demonstration purposes because
of the large size of the entire dataset. For exploration and downloading
of all the data from the paper, we refer to: [Liver Atlas Data
Portal](https://www.livercellatlas.org/). For the code used for all the
Differential NicheNet analyses on the entire liver cell atlas dataset,
see <https://github.com/saeyslab/NicheNet_LiverCellAtlas>.

The goal of Differential NicheNet is to predict ligand-receptors pairs
that are both differentially expressed and active between different
niches of interest.

In this vignette, we will look at cell-cell communication differences
between the Kupffer cell niche, the bile duct macrophage niche, and the
capsule macrophage niche, with the macrophages in each niche as receiver
cell of interest. This means that we are interested in identifying the
niche-specific ligands important for the identity of each of these
macrophage subtypes.

# 0. Read in the expression data of interest, and the NicheNet ligand-receptor network and ligand-target matrix

## Load in packages

``` r
library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) 
```

## Read in the expression data

``` r
seurat_obj = readRDS(url("https://zenodo.org/record/5840787/files/seurat_obj_subset_integrated_zonation.rds"))
DimPlot(seurat_obj, group.by = "celltype", label = TRUE) 
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
seurat_obj = SetIdent(seurat_obj, value = "celltype")
```

As you can see, the LSECs, hepatocytes and Stellate cells are each
divided in two groups, based on their spatial location (periportal and
pericentral).

## Read in the NicheNet ligand-receptor network and ligand-target matrix

The used ligand-receptor network and ligand-target matrix can be
downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758).
The Seurat object containing expression data of interacting cells in
HNSCC can also be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4675430.svg)](https://doi.org/10.5281/zenodo.4675430).

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

head(lr_network)
## # A tibble: 6 x 3
##   ligand receptor bonafide
##   <chr>  <chr>    <lgl>   
## 1 CXCL1  CXCR2    TRUE    
## 2 CXCL2  CXCR2    TRUE    
## 3 CXCL3  CXCR2    TRUE    
## 4 CXCL5  CXCR2    TRUE    
## 5 PPBP   CXCR2    TRUE    
## 6 CXCL6  CXCR2    TRUE
```

Note: because the data is of mouse origin: we need to convert human gene
symbols to their murine one-to-one orthologs

``` r
organism = "mouse" 
```

``` r
if(organism == "mouse"){
  lr_network = lr_network %>% mutate(ligand = convert_human_to_mouse_symbols(ligand), receptor = convert_human_to_mouse_symbols(receptor)) %>% drop_na()

  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
}
```

# 1. Define the niches/microenvironments of interest

Each niche should have at least one “sender/niche” cell population and
one “receiver/target” cell population (present in your expression data)

In this case study, we are interested to find differences in cell-cell
interactions to hepatic macrophages in three different niches: 1) the
Kupffer cell niche, 2) the bile-duct or lipid-associated macrophage
niche, and 3) the capsule macrophage niche.

Based on imaging and spatial transcriptomics, the composition of each
niche was defined as follows:

The receiver cell population in the Kupffer cell niche is the “KCs” cell
type, the sender cell types are: “LSECs\_portal,”“Hepatocytes\_portal,”
and “Stellate cells\_portal.” The receiver cell population in the
lipid-associated macrophage (MoMac2) niche is the “MoMac2” cell type,
the sender cell types are: “Cholangiocytes,” and “Fibroblast 2.” The
receiver cell population in the capsule macrophage (MoMac1) niche is the
“MoMac1” cell type, the sender cell types are: “Capsule fibroblasts,”
and “Mesothelial cells.”

! Important: your receiver cell type should consist of 1 cluster!

``` r
niches = list(
    "KC_niche" = list(
      "sender" = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
      "receiver" = c("KCs")),
    "MoMac2_niche" = list(
      "sender" = c("Cholangiocytes","Fibroblast 2"),
      "receiver" = c("MoMac2")),
    "MoMac1_niche" = list(
      "sender" = c("Capsule fibroblasts","Mesothelial cells"),
      "receiver" = c("MoMac1"))
  )
```

# 2. Calculate differential expression between the niches

In this step, we will determine DE between the different niches for both
senders and receivers to define the DE of L-R pairs.

### Calculate DE

The method to calculate the differential expression is here the standard
Seurat Wilcoxon test, but this can be replaced if wanted by the user
(only requirement: output tables `DE_sender_processed` and
`DE_receiver_processed` should be in the same format as shown here).

DE will be calculated for each pairwise sender (or receiver) cell type
comparision between the niches (so across niches, not within niche). In
our case study, this means e.g. that DE of LSECs\_portal ligands will be
calculated by DE analysis of LSECs\_portal vs Cholangiocytes;
LSECs\_portal vs Fibroblast 2; LSECs\_portal vs Capsule fibroblasts; and
LSECs\_portal vs Mesothelial cells. We split the cells per cell type
instead of merging all cells from the other niche to avoid that the DE
analysis will be driven by the most abundant cell types.

``` r
assay_oi = "SCT" # other possibilities: RNA,...
DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% intersect(rownames(seurat_obj))), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
## [1] "Calculate Sender DE between: LSECs_portal and Cholangiocytes"     
## [2] "Calculate Sender DE between: LSECs_portal and Fibroblast 2"       
## [3] "Calculate Sender DE between: LSECs_portal and Capsule fibroblasts"
## [4] "Calculate Sender DE between: LSECs_portal and Mesothelial cells"  
## [1] "Calculate Sender DE between: Hepatocytes_portal and Cholangiocytes"     
## [2] "Calculate Sender DE between: Hepatocytes_portal and Fibroblast 2"       
## [3] "Calculate Sender DE between: Hepatocytes_portal and Capsule fibroblasts"
## [4] "Calculate Sender DE between: Hepatocytes_portal and Mesothelial cells"  
## [1] "Calculate Sender DE between: Stellate cells_portal and Cholangiocytes"     
## [2] "Calculate Sender DE between: Stellate cells_portal and Fibroblast 2"       
## [3] "Calculate Sender DE between: Stellate cells_portal and Capsule fibroblasts"
## [4] "Calculate Sender DE between: Stellate cells_portal and Mesothelial cells"  
## [1] "Calculate Sender DE between: Cholangiocytes and LSECs_portal"         
## [2] "Calculate Sender DE between: Cholangiocytes and Hepatocytes_portal"   
## [3] "Calculate Sender DE between: Cholangiocytes and Stellate cells_portal"
## [4] "Calculate Sender DE between: Cholangiocytes and Capsule fibroblasts"  
## [5] "Calculate Sender DE between: Cholangiocytes and Mesothelial cells"    
## [1] "Calculate Sender DE between: Fibroblast 2 and LSECs_portal"         
## [2] "Calculate Sender DE between: Fibroblast 2 and Hepatocytes_portal"   
## [3] "Calculate Sender DE between: Fibroblast 2 and Stellate cells_portal"
## [4] "Calculate Sender DE between: Fibroblast 2 and Capsule fibroblasts"  
## [5] "Calculate Sender DE between: Fibroblast 2 and Mesothelial cells"    
## [1] "Calculate Sender DE between: Capsule fibroblasts and LSECs_portal"         
## [2] "Calculate Sender DE between: Capsule fibroblasts and Hepatocytes_portal"   
## [3] "Calculate Sender DE between: Capsule fibroblasts and Stellate cells_portal"
## [4] "Calculate Sender DE between: Capsule fibroblasts and Cholangiocytes"       
## [5] "Calculate Sender DE between: Capsule fibroblasts and Fibroblast 2"         
## [1] "Calculate Sender DE between: Mesothelial cells and LSECs_portal"         
## [2] "Calculate Sender DE between: Mesothelial cells and Hepatocytes_portal"   
## [3] "Calculate Sender DE between: Mesothelial cells and Stellate cells_portal"
## [4] "Calculate Sender DE between: Mesothelial cells and Cholangiocytes"       
## [5] "Calculate Sender DE between: Mesothelial cells and Fibroblast 2"

DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
## # A tibble: 3 x 2
##   receiver receiver_other_niche
##   <chr>    <chr>               
## 1 KCs      MoMac2              
## 2 KCs      MoMac1              
## 3 MoMac2   MoMac1              
## [1] "Calculate receiver DE between: KCs and MoMac2" "Calculate receiver DE between: KCs and MoMac1"
## [1] "Calculate receiver DE between: MoMac2 and KCs"    "Calculate receiver DE between: MoMac2 and MoMac1"
## [1] "Calculate receiver DE between: MoMac1 and KCs"    "Calculate receiver DE between: MoMac1 and MoMac2"
```

### Process DE results:

``` r
expression_pct = 0.10
DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
```

### Combine sender-receiver DE based on L-R pairs:

As mentioned above, DE of ligands from one sender cell type is
determined be calculating DE between that cell type, and all the sender
cell types of the other niche. To summarize the DE of ligands of that
cell type we have several options: we could take the average LFC, but
also the minimum LFC compared to the other niche. We recommend using the
minimum LFC, because this is the strongest specificity measure of ligand
expression, because a high min LFC means that a ligand is more strongly
expressed in the cell type of niche 1 compared to all cell types of
niche 2 (in contrast to a high average LFC, which does not exclude that
one or more cell types in niche 2 also strongly express that ligand).

``` r
specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
```

# 3. Optional: Calculate differential expression between the different spatial regions

To improve the cell-cell interaction predictions, you can consider
spatial information if possible and applicable. Spatial information can
come from microscopy data, or from spatial transcriptomics data such as
Visium.

There are several ways to incorporate spatial information in the
Differential NicheNet pipeline. First, you can only consider cell types
as belonging to the same niche if they are in the same spatial location.
Another way is including spatial differential expression of
ligand-receptor pairs within one cell type in the prioritization
framework.

For example: We have a cell type X, located in regions A and B, and we
want to study cell-cell communication in region A. We first add only
celltypeX of regionA in the niche definition, and then calculate DE
between celltypeX-regionA and celltypeX-regionB to give higher
prioritization weight to regionA-specific ligands.

In this case study, our region of interest is the periportal region of
the liver, because KCs in mouse are predominantly located in the
periportal region. Therefore we will give higher weight to ligands that
are in the niche cells of KCs higher expressed in the periportal
compared to the pericentral region.

We do this as follows, by first defining a ‘spatial info’ dataframe. If
there is no spatial information in your data: set the following two
parameters to FALSE, and make a mock ‘spatial\_info’ data frame.

``` r
include_spatial_info_sender = TRUE # if not spatial info to include: put this to false 
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true 
```

``` r
spatial_info = tibble(celltype_region_oi = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"), 
                      celltype_other_region = c("LSECs_central","Hepatocytes_central","Stellate cells_central")
                      ) %>% 
  mutate(niche =  "KC_niche", celltype_type = "sender")
specificity_score_spatial = "lfc"
```

``` r
# this is how this should be defined if you don't have spatial info
# mock spatial info
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
} 
```

``` r
if(include_spatial_info_sender == TRUE){
  sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"), assay_oi = assay_oi)
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)

  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))

} else {
  # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
  sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  

}
## [1] "Calculate Spatial DE between: LSECs_portal and LSECs_central"
## [1] "Calculate Spatial DE between: Hepatocytes_portal and Hepatocytes_central"
## [1] "Calculate Spatial DE between: Stellate cells_portal and Stellate cells_central"
```

``` r
if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"), assay_oi = assay_oi)
  receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
  receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)

  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))

} else {
    # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}
```

# 4. Calculate ligand activities and infer active ligand-target links

In this step, we will predict ligand activities of each ligand for each
of the receiver cell types across the different niches. This is similar
to the ligand activity analysis done in the normal NicheNet pipeline.

To calculate ligand activities, we first need to define a geneset of
interest for each niche. In this case study, the geneset of interest for
the Kupffer cell niche are the genes upregulated in Kupffer cells
compared to the capsule and bile duct macrophages. The geneset of
interest for the bile duct macrophage niche are the genes upregulated in
bile duct macrophages compared to the capsule macrophages and Kupffer
cells. And similarly for the capsule macrophage geneset of interest.

Note that you can also define these geneset of interest in a different
way! (eg pathway-based geneset etc)

Ligand-target links are inferred in the same way as described in the
basic NicheNet vignettes.

``` r
lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"

DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
## [1] "Calculate receiver DE between: KCs and MoMac2" "Calculate receiver DE between: KCs and MoMac1"
## [1] "Calculate receiver DE between: MoMac2 and KCs"    "Calculate receiver DE between: MoMac2 and MoMac1"
## [1] "Calculate receiver DE between: MoMac1 and KCs"    "Calculate receiver DE between: MoMac1 and MoMac2"
DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
  
background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_KC = DE_receiver_processed_targets %>% filter(receiver == niches$KC_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_MoMac2 = DE_receiver_processed_targets %>% filter(receiver == niches$MoMac2_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_MoMac1 = DE_receiver_processed_targets %>% filter(receiver == niches$MoMac1_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
# Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
# If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
geneset_KC %>% setdiff(rownames(ligand_target_matrix))
##  [1] "Fcna"          "Wfdc17"        "C4b"           "AW112010"      "mt-Co1"        "Adgre4"        "Pira2"        
##  [8] "mt-Nd2"        "mt-Co3"        "mt-Co2"        "mt-Nd3"        "mt-Atp6"       "mt-Nd4"        "mt-Nd1"       
## [15] "Iigp1"         "Ear2"          "2900097C17Rik" "Anapc15"       "B430306N03Rik" "Trim30a"       "Pilrb2"       
## [22] "Gbp8"          "Arf2"          "AC149090.1"    "Xlr"           "Cd209f"        "mt-Cytb"       "Ifitm6"       
## [29] "Mndal"         "Gm4951"        "Ifi205"        "Serpina3g"
geneset_MoMac2 %>% setdiff(rownames(ligand_target_matrix))
##  [1] "Chil3"         "Lyz1"          "Ccl9"          "Ly6c2"         "Tmsb10"        "Gm21188"       "Calm3"        
##  [8] "S100a11"       "Ftl1-ps1"      "Gm10076"       "Ms4a6c"        "Atp5e"         "Snrpe"         "Clec4a3"      
## [15] "Ly6i"          "1810058I24Rik" "Aph1c"         "Cox6c"         "Atp5o.1"       "Rpl34"         "Cbr2"         
## [22] "Rtf2"          "Gm10073"       "Snhg6"         "Clec2i"        "AI413582"      "Ggta1"         "Ppp1cc"       
## [29] "Rpl10-ps3"     "Eif2s3y"       "Gstp1"         "Gm36161"       "Cyp2c70"       "Mup21"         "Ces3a"        
## [36] "Rps12-ps3"
geneset_MoMac1 %>% setdiff(rownames(ligand_target_matrix))
## [1] "H2-Ab1"  "Malat1"  "H2-Aa"   "Gm26522" "H2-M2"   "Mgl2"    "Klra2"   "H2-D1"   "H2-Q6"

length(geneset_KC)
## [1] 494
length(geneset_MoMac2)
## [1] 611
length(geneset_MoMac1)
## [1] 80
```

It is always useful to check the number of genes in the geneset before
doing the ligand activity analysis. We recommend having between 20 and
1000 genes in the geneset of interest, and a background of at least 5000
genes for a proper ligand activity analysis. If you retrieve too many DE
genes, it is recommended to use a higher `lfc_cutoff` threshold. We
recommend using a cutoff of 0.15 if you have &gt; 2 receiver
cells/niches to compare and use the min\_lfc as specificity score. If
you have only 2 receivers/niche, we recommend using a higher threshold
(such as using 0.25). If you have single-cell data like Smart-seq2 with
high sequencing depth, we recommend to also use higher threshold.

``` r
top_n_target = 250

niche_geneset_list = list(
    "KC_niche" = list(
      "receiver" = "KCs",
      "geneset" = geneset_KC,
      "background" = background),
    "MoMac1_niche" = list(
      "receiver" = "MoMac1",
      "geneset" = geneset_MoMac1 ,
      "background" = background),
    "MoMac2_niche" = list(
      "receiver" = "MoMac2",
      "geneset" = geneset_MoMac2 ,
      "background" = background)  
  )
  
ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
## [1] "Calculate Ligand activities for: KCs"
## [1] "Calculate Ligand activities for: MoMac1"
## [1] "Calculate Ligand activities for: MoMac2"
```

# 5. Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)

In this step, we will calculate average (scaled) expression, and
fraction of expression, of ligands, receptors, and target genes across
all cell types of interest. Now this is here demonstrated via the
DotPlot function of Seurat, but this can also be done via other ways of
course.

``` r
features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
  
dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl = dotplot$data %>% as_tibble()
exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
```

``` r
exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
```

# 6. Expression fraction and receptor

In this step, we will score ligand-receptor interactions based on
expression strength of the receptor, in such a way that we give higher
scores to the most strongly expressed receptor of a certain ligand, in a
certain celltype. This will not effect the rank of individual ligands
later on, but will help in prioritizing the most important receptors per
ligand (next to other factors regarding the receptor - see later).

``` r
exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 
```

# 7. Prioritization of ligand-receptor and ligand-target links

In this step, we will combine all the above calculated information to
prioritize ligand-receptor-target links. We scale every property of
interest between 0 and 1, and the final prioritization score is a
weighted sum of the scaled scores of all the properties of interest.

We provide the user the option to consider the following properties for
prioritization (of which the weights are defined in
`prioritizing_weights`) :

-   Ligand DE score: niche-specific expression of the ligand: by
    default, this the minimum logFC between the sender of interest and
    all the senders of the other niche(s). The higher the min logFC, the
    higher the niche-specificity of the ligand. Therefore we recommend
    to give this factor a very high weight. `prioritizing_weights`
    argument: `"scaled_ligand_score"`. Recommended weight: 5 (at least
    1, max 5).

-   Scaled ligand expression: scaled expression of a ligand in one
    sender compared to the other cell types in the dataset. This might
    be useful to rescue potentially interesting ligands that have a high
    scaled expression value, but a relatively small min logFC compared
    to the other niche. One reason why this logFC might be small occurs
    when (some) genes are not picked up efficiently by the used
    sequencing technology (or other reasons for low RNA expression of
    ligands). For example, we have observed that many ligands from the
    Tgf-beta/BMP family are not picked up efficiently with single-nuclei
    RNA sequencing compared to single-cell sequencing.
    `prioritizing_weights` argument:
    `"scaled_ligand_expression_scaled"`. Recommended weight: 1 (unless
    technical reason for lower gene detection such as while using
    Nuc-seq: then recommended to use a higher weight: 2).

-   Ligand expression fraction: Ligands that are expressed in a smaller
    fraction of cells of a cell type than defined by
    `exprs_cutoff`(default: 0.10) will get a lower ranking, proportional
    to their fraction (eg ligand expressed in 9% of cells will be ranked
    higher than ligand expressed in 0.5% of cells). We opted for this
    weighting based on fraction, instead of removing ligands that are
    not expressed in more cells than this cutoff, because some
    interesting ligands could be removed that way. Fraction of
    expression is not taken into account for the prioritization if it is
    already higher than the cutoff. `prioritizing_weights` argument:
    `"ligand_fraction"`. Recommended weight: 1.

-   Ligand spatial DE score: spatial expression specificity of the
    ligand. If the niche of interest is at a specific tissue location,
    but some of the sender cell types of that niche are also present in
    other locations, it can be very informative to further prioritize
    ligands of that sender by looking how they are DE between the
    spatial location of interest compared to the other locations.
    `prioritizing_weights` argument: `"scaled_ligand_score_spatial"`.
    Recommended weight: 2 (or 0 if not applicable).

-   Receptor DE score: niche-specific expression of the receptor: by
    default, this the minimum logFC between the receiver of interest and
    all the receiver of the other niche(s). The higher the min logFC,
    the higher the niche-specificity of the receptor. Based on our
    experience, we don’t suggest to give this as high importance as the
    ligand DE, but this might depend on the specific case study.
    `prioritizing_weights` argument: `"scaled_receptor_score"`.
    Recommended weight: 0.5 (at least 0.5, and lower than
    `"scaled_ligand_score"`).

-   Scaled receptor expression: scaled expression of a receptor in one
    receiver compared to the other cell types in the dataset. This might
    be useful to rescue potentially interesting receptors that have a
    high scaled expression value, but a relatively small min logFC
    compared to the other niche. One reason why this logFC might be
    small occurs when (some) genes are not picked up efficiently by the
    used sequencing technology. `prioritizing_weights` argument:
    `"scaled_receptor_expression_scaled"`. Recommended weight: 0.5.

-   Receptor expression fraction: Receptors that are expressed in a
    smaller fraction of cells of a cell type than defined by
    `exprs_cutoff`(default: 0.10) will get a lower ranking, proportional
    to their fraction (eg receptor expressed in 9% of cells will be
    ranked higher than receptor expressed in 0.5% of cells). We opted
    for this weighting based on fraction, instead of removing receptors
    that are not expressed in more cells than this cutoff, because some
    interesting receptors could be removed that way. Fraction of
    expression is not taken into account for the prioritization if it is
    already higher than the cutoff. `prioritizing_weights` argument:
    `"receptor_fraction"`. Recommended weight: 1.

-   Receptor expression strength: this factor let us give higher weights
    to the most highly expressed receptor of a ligand in the receiver.
    This let us rank higher one member of a receptor family if it higher
    expressed than the other members. `prioritizing_weights` argument:
    `"ligand_scaled_receptor_expression_fraction"`. Recommended value: 1
    (minimum: 0.5).

-   Receptor spatial DE score: spatial expression specificity of the
    receptor. If the niche of interest is at a specific tissue location,
    but the receiver cell type of that niche is also present in other
    locations, it can be very informative to further prioritize
    receptors of that receiver by looking how they are DE between the
    spatial location of interest compared to the other locations.
    `prioritizing_weights` argument: `"scaled_receptor_score_spatial"`.
    Recommended weight: 1 (or 0 if not applicable).

-   Absolute ligand activity: to further prioritize ligand-receptor
    pairs based on their predicted effect of the ligand-receptor
    interaction on the gene expression in the receiver cell type -
    absolute ligand activity accords to ‘absolute’ enrichment of target
    genes of a ligand within the affected receiver genes.
    `prioritizing_weights` argument: `"scaled_activity"`. Recommended
    weight: 0, unless absolute enrichment of target genes is of specific
    interest.

-   Normalized ligand activity: to further prioritize ligand-receptor
    pairs based on their predicted effect of the ligand-receptor
    interaction on the gene expression in the receiver cell type -
    normalization of activity is done because we found that some
    datasets/conditions/niches have higher baseline activity values than
    others - normalized ligand activity accords to ‘relative’ enrichment
    of target genes of a ligand within the affected receiver genes.
    `prioritizing_weights` argument: `"scaled_activity_normalized"`.
    Recommended weight: at least 1.

-   Prior knowledge quality of the L-R interaction: the NicheNet LR
    network consists of two types of interactions: L-R pairs documented
    in curated databases, and L-R pairs predicted based on gene
    annotation and PPIs. The former are categorized as ‘bona fide’
    interactions. To rank bona fide interactions higher, but not exlude
    potentially interesting non-bona-fide ones, we give bona fide
    interactions a score of 1, and non-bona-fide interactions a score
    fof 0.5. `prioritizing_weights` argument: `"bona_fide"` Recommend
    weight: at least 1.

``` r
prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                          "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)
```

``` r
output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
         ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
prioritization_tables = get_prioritization_tables(output, prioritizing_weights)

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
## # A tibble: 10 x 37
##    niche    receiver sender ligand_receptor ligand receptor bonafide ligand_score ligand_signific~ ligand_present ligand_expressi~
##    <chr>    <chr>    <chr>  <chr>           <chr>  <chr>    <lgl>           <dbl>            <dbl>          <dbl>            <dbl>
##  1 KC_niche KCs      Hepat~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE            3.46                1              1            22.4 
##  2 KC_niche KCs      Hepat~ Apoa1--Msr1     Apoa1  Msr1     FALSE            3.46                1              1            22.4 
##  3 KC_niche KCs      Hepat~ Apoa1--Abca1    Apoa1  Abca1    FALSE            3.46                1              1            22.4 
##  4 KC_niche KCs      Hepat~ Apoa1--Scarb1   Apoa1  Scarb1   FALSE            3.46                1              1            22.4 
##  5 KC_niche KCs      Hepat~ Apoa1--Derl1    Apoa1  Derl1    FALSE            3.46                1              1            22.4 
##  6 KC_niche KCs      Hepat~ Apoa1--Atp5b    Apoa1  Atp5b    FALSE            3.46                1              1            22.4 
##  7 KC_niche KCs      Hepat~ Serpina1a--Lrp1 Serpi~ Lrp1     TRUE             3.02                1              1            10.4 
##  8 KC_niche KCs      Hepat~ Trf--Tfrc       Trf    Tfrc     TRUE             1.74                1              1             9.96
##  9 KC_niche KCs      LSECs~ Cxcl10--Fpr1    Cxcl10 Fpr1     FALSE            2.07                1              1             3.56
## 10 KC_niche KCs      LSECs~ Cxcl10--Ccr5    Cxcl10 Ccr5     FALSE            2.07                1              1             3.56
## # ... with 26 more variables: ligand_expression_scaled <dbl>, ligand_fraction <dbl>, ligand_score_spatial <dbl>,
## #   receptor_score <dbl>, receptor_significant <dbl>, receptor_present <dbl>, receptor_expression <dbl>,
## #   receptor_expression_scaled <dbl>, receptor_fraction <dbl>, receptor_score_spatial <dbl>,
## #   ligand_scaled_receptor_expression_fraction <dbl>, avg_score_ligand_receptor <dbl>, activity <dbl>, activity_normalized <dbl>,
## #   scaled_ligand_score <dbl>, scaled_ligand_expression_scaled <dbl>, scaled_receptor_score <dbl>,
## #   scaled_receptor_expression_scaled <dbl>, scaled_avg_score_ligand_receptor <dbl>, scaled_ligand_score_spatial <dbl>,
## #   scaled_receptor_score_spatial <dbl>, scaled_ligand_fraction_adapted <dbl>, scaled_receptor_fraction_adapted <dbl>, ...
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
## # A tibble: 10 x 20
##    niche    receiver sender         ligand_receptor ligand receptor bonafide target  target_score target_significa~ target_present
##    <chr>    <chr>    <chr>          <chr>           <chr>  <chr>    <lgl>    <chr>          <dbl>             <dbl>          <dbl>
##  1 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Abca1          0.225                 1              1
##  2 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Actb           0.341                 1              1
##  3 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Ehd1           0.353                 1              1
##  4 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Ets2           0.191                 1              1
##  5 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Hmox1          1.26                  1              1
##  6 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Nr2f2          0.169                 1              1
##  7 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Sgk1           0.443                 1              1
##  8 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Sptbn1         0.166                 1              1
##  9 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Tcf7l2         1.01                  1              1
## 10 KC_niche KCs      Hepatocytes_p~ Apoa1--Lrp1     Apoa1  Lrp1     FALSE    Tsc22d3        0.346                 1              1
## # ... with 9 more variables: target_expression <dbl>, target_expression_scaled <dbl>, target_fraction <dbl>,
## #   ligand_target_weight <dbl>, activity <dbl>, activity_normalized <dbl>, scaled_activity <dbl>,
## #   scaled_activity_normalized <dbl>, prioritization_score <dbl>

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
## # A tibble: 10 x 37
##    niche   receiver sender  ligand_receptor ligand receptor bonafide ligand_score ligand_signific~ ligand_present ligand_expressi~
##    <chr>   <chr>    <chr>   <chr>           <chr>  <chr>    <lgl>           <dbl>            <dbl>          <dbl>            <dbl>
##  1 MoMac2~ MoMac2   Cholan~ Spp1--Cd44      Spp1   Cd44     TRUE             6.60                1              1           108.  
##  2 MoMac2~ MoMac2   Cholan~ Spp1--Itga4     Spp1   Itga4    TRUE             6.60                1              1           108.  
##  3 MoMac2~ MoMac2   Cholan~ Spp1--Itgb5     Spp1   Itgb5    TRUE             6.60                1              1           108.  
##  4 MoMac2~ MoMac2   Cholan~ Spp1--Itgav     Spp1   Itgav    TRUE             6.60                1              1           108.  
##  5 MoMac2~ MoMac2   Cholan~ Spp1--Itgb1     Spp1   Itgb1    TRUE             6.60                1              1           108.  
##  6 MoMac2~ MoMac2   Cholan~ Spp1--Itga9     Spp1   Itga9    TRUE             6.60                1              1           108.  
##  7 MoMac2~ MoMac2   Cholan~ Spp1--Ncstn     Spp1   Ncstn    FALSE            6.60                1              1           108.  
##  8 MoMac2~ MoMac2   Cholan~ Spp1--Itga5     Spp1   Itga5    FALSE            6.60                1              1           108.  
##  9 MoMac2~ MoMac2   Cholan~ Cyr61--Itgb2    Cyr61  Itgb2    TRUE             1.14                1              1             4.54
## 10 MoMac2~ MoMac2   Cholan~ Spp1--Sdc1      Spp1   Sdc1     FALSE            6.60                1              1           108.  
## # ... with 26 more variables: ligand_expression_scaled <dbl>, ligand_fraction <dbl>, ligand_score_spatial <dbl>,
## #   receptor_score <dbl>, receptor_significant <dbl>, receptor_present <dbl>, receptor_expression <dbl>,
## #   receptor_expression_scaled <dbl>, receptor_fraction <dbl>, receptor_score_spatial <dbl>,
## #   ligand_scaled_receptor_expression_fraction <dbl>, avg_score_ligand_receptor <dbl>, activity <dbl>, activity_normalized <dbl>,
## #   scaled_ligand_score <dbl>, scaled_ligand_expression_scaled <dbl>, scaled_receptor_score <dbl>,
## #   scaled_receptor_expression_scaled <dbl>, scaled_avg_score_ligand_receptor <dbl>, scaled_ligand_score_spatial <dbl>,
## #   scaled_receptor_score_spatial <dbl>, scaled_ligand_fraction_adapted <dbl>, scaled_receptor_fraction_adapted <dbl>, ...
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
## # A tibble: 10 x 20
##    niche        receiver sender      ligand_receptor ligand receptor bonafide target target_score target_significa~ target_present
##    <chr>        <chr>    <chr>       <chr>           <chr>  <chr>    <lgl>    <chr>         <dbl>             <dbl>          <dbl>
##  1 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Ahnak         1.38                  1              1
##  2 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Capn2         0.238                 1              1
##  3 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Cdkn1a        0.779                 1              1
##  4 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Cxcr4         0.486                 1              1
##  5 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Dhrs3         0.477                 1              1
##  6 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Fam12~        0.178                 1              1
##  7 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Fn1           0.545                 1              1
##  8 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Gadd4~        0.245                 1              1
##  9 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Gapdh         0.681                 1              1
## 10 MoMac2_niche MoMac2   Cholangioc~ Spp1--Cd44      Spp1   Cd44     TRUE     Gdf15         0.643                 1              1
## # ... with 9 more variables: target_expression <dbl>, target_expression_scaled <dbl>, target_fraction <dbl>,
## #   ligand_target_weight <dbl>, activity <dbl>, activity_normalized <dbl>, scaled_activity <dbl>,
## #   scaled_activity_normalized <dbl>, prioritization_score <dbl>

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[3]]$receiver) %>% head(10)
## # A tibble: 10 x 37
##    niche  receiver sender   ligand_receptor ligand receptor bonafide ligand_score ligand_signific~ ligand_present ligand_expressi~
##    <chr>  <chr>    <chr>    <chr>           <chr>  <chr>    <lgl>           <dbl>            <dbl>          <dbl>            <dbl>
##  1 MoMac~ MoMac1   Mesothe~ C3--C3ar1       C3     C3ar1    TRUE             3.81                1              1             33.7
##  2 MoMac~ MoMac1   Capsule~ C3--C3ar1       C3     C3ar1    TRUE             3.73                1              1             31.7
##  3 MoMac~ MoMac1   Mesothe~ C3--Itgb2       C3     Itgb2    TRUE             3.81                1              1             33.7
##  4 MoMac~ MoMac1   Mesothe~ C3--Itgax       C3     Itgax    TRUE             3.81                1              1             33.7
##  5 MoMac~ MoMac1   Capsule~ C3--Itgb2       C3     Itgb2    TRUE             3.73                1              1             31.7
##  6 MoMac~ MoMac1   Mesothe~ C3--Lrp1        C3     Lrp1     TRUE             3.81                1              1             33.7
##  7 MoMac~ MoMac1   Capsule~ C3--Itgax       C3     Itgax    TRUE             3.73                1              1             31.7
##  8 MoMac~ MoMac1   Capsule~ C3--Lrp1        C3     Lrp1     TRUE             3.73                1              1             31.7
##  9 MoMac~ MoMac1   Capsule~ Rarres2--Cmklr1 Rarre~ Cmklr1   TRUE             2.67                1              1             24.1
## 10 MoMac~ MoMac1   Mesothe~ C3--Ccr5        C3     Ccr5     FALSE            3.81                1              1             33.7
## # ... with 26 more variables: ligand_expression_scaled <dbl>, ligand_fraction <dbl>, ligand_score_spatial <dbl>,
## #   receptor_score <dbl>, receptor_significant <dbl>, receptor_present <dbl>, receptor_expression <dbl>,
## #   receptor_expression_scaled <dbl>, receptor_fraction <dbl>, receptor_score_spatial <dbl>,
## #   ligand_scaled_receptor_expression_fraction <dbl>, avg_score_ligand_receptor <dbl>, activity <dbl>, activity_normalized <dbl>,
## #   scaled_ligand_score <dbl>, scaled_ligand_expression_scaled <dbl>, scaled_receptor_score <dbl>,
## #   scaled_receptor_expression_scaled <dbl>, scaled_avg_score_ligand_receptor <dbl>, scaled_ligand_score_spatial <dbl>,
## #   scaled_receptor_score_spatial <dbl>, scaled_ligand_fraction_adapted <dbl>, scaled_receptor_fraction_adapted <dbl>, ...
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[3]]$receiver) %>% head(10)
## # A tibble: 10 x 20
##    niche        receiver sender       ligand_receptor ligand receptor bonafide target target_score target_signific~ target_present
##    <chr>        <chr>    <chr>        <chr>           <chr>  <chr>    <lgl>    <chr>         <dbl>            <dbl>          <dbl>
##  1 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Btg2          0.690                1              1
##  2 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Ccl12         0.375                1              1
##  3 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Ccnd2         0.605                1              1
##  4 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Il1b          1.05                 1              1
##  5 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Jun           0.862                1              1
##  6 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Pdgfb         0.285                1              1
##  7 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Tle1          0.154                1              1
##  8 MoMac1_niche MoMac1   Mesothelial~ C3--C3ar1       C3     C3ar1    TRUE     Ubc           0.306                1              1
##  9 MoMac1_niche MoMac1   Capsule fib~ C3--C3ar1       C3     C3ar1    TRUE     Btg2          0.690                1              1
## 10 MoMac1_niche MoMac1   Capsule fib~ C3--C3ar1       C3     C3ar1    TRUE     Ccl12         0.375                1              1
## # ... with 9 more variables: target_expression <dbl>, target_expression_scaled <dbl>, target_fraction <dbl>,
## #   ligand_target_weight <dbl>, activity <dbl>, activity_normalized <dbl>, scaled_activity <dbl>,
## #   scaled_activity_normalized <dbl>, prioritization_score <dbl>

prioritization_tables$prioritization_tbl_ligand_receptor = prioritization_tables$prioritization_tbl_ligand_receptor %>% mutate(receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac2_niche"))) 
prioritization_tables$prioritization_tbl_ligand_target = prioritization_tables$prioritization_tbl_ligand_target %>% mutate(receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac2_niche"))) 
```

# 8. Visualization of the Differential NicheNet output

## Differential expression of ligand and expression

Before visualization, we need to define the most important
ligand-receptor pairs per niche. We will do this by first determining
for which niche the highest score is found for each
ligand/ligand-receptor pair. And then getting the top 50 ligands per
niche.

``` r
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche
```

Now we will look first at the top ligand-receptor pairs for KCs (here,
we will take the top 2 scoring receptors per prioritized ligand)

``` r
receiver_oi = "KCs" 

filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
```

Visualization: minimum LFC compared to other niches

``` r
lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Show the spatialDE as additional information

``` r
lfc_plot_spatial = make_ligand_receptor_lfc_spatial_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, ligand_spatial = include_spatial_info_sender, receptor_spatial = include_spatial_info_receiver, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot_spatial
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

From this plot, you can see that some KC-niche ligands like Dll4 (by
LSEC) and Il34 (by Stellate cells) are higher expressed in the
periportal LSEC/stellate cells vs the pericentral ones. This can be
interesting information knowing that KCs are mainly located
periportally. However, other ligands like Gdf2 (by Stellate cells) are
not preferentially expressed by periportal stellate cells, but this does
not mean they cannot be interesting. As you can see in the following
figure, this ligand has one of the highest ligand activities, meaning
that there is a strong enrichment of its target genes among the
KC-specific genes.

## Ligand expression, activity and target genes

Active target gene inference - cf Default NicheNet

Now: visualization of ligand activity and ligand-target links

``` r
exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_activity_target_plot$combined_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

On this plot, we can see that some strongly DE ligand-receptor pairs in
the KC niche, have also high scaled ligand activity on KCs - making them
strong predictions for further validation. An example of this is Gdf2
and Bmp10, who bind the receptor Acvrl1 (ALK1). The role of
Gdf2/Bmp10-Acvrl1 in KC development was experimentally validated in the
Guilliams et al paper.

**important: ligand-receptor pairs with both high differential
expression and ligand activity (=target gene enrichment) are very
interesting predictions as key regulators of your intercellular
communication process of interest !**

If this plot contains too much information because we look at many hits
(top 50 ligands), you can make this plot of course for less ligands as
well, eg for the top20.

``` r
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(20, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_activity_target_plot$combined_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

## Circos plot of prioritized ligand-receptor pairs

Because a top50 is too much to visualize in a circos plot, we will only
visualize the top 15.

``` r
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(15, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% length(), name = 'Spectral') %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->![](differential_nichenet_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

``` r
# circos_output$p_circos
```

## Visualization for the other liver macrophages: central vein

``` r
receiver_oi = "MoMac1"  
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(50, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## Visualization for the other liver macrophages: bile duct

``` r
receiver_oi = "MoMac2"  
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(50, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

# Notes, limitations, and comparison to default NicheNet.

In the default NicheNet pipeline, expressed ligand-receptor pairs are
prioritized based on their ligand activity alone. Here, in the
Differential NicheNet pipeline, we also draw information based on
differential expression of the L-R pairs compared to other niches (and
if applicable: other spatial locations.)

Because we here focus on differential expression of ligand-receptor
pairs, and by using the default prioritizations weights more on DE than
activity, we tend to find many different hits than with the default
NicheNet pipeline. With Differential NicheNet, we tend to find more
high-DE, low-activity hits, whereas with default NicheNet we find more
low-DE, high-activity hits.

It should be noted that some of the high-DE, low-activity hits might be
really important because they just have low NicheNet activities due to
limitations in the NicheNet activity prediction (eg improper/incomplete
prior knowledge within NicheNet for that ligand), but some of them might
also be high in DE but not in activity because they don’t have strong
signaling effects (eg ligands involved in cell adhesion only).

For the opposite pairs with low-DE and high-activity that are not
strongly prioritized by Differential NicheNet, the following should be
considered: 1) some ligands are regulated post-transcriptionally, and
that the high predicted activities might still reflect true signaling;
2) high predicted activity values might be due to limitations of
NicheNet (inaccurate prior knowledge) and these lowDE ligands are not
important in the biological process of interest (although a highDE
family member of this ligand may! since signaling between family members
tends to be very similar); 3) high activity in one condition might be
due to downregulation in the other condition, leading to high activity
but low DE. Currently, ligand activities are automatically calculated on
upregulated genes per condition, but downregulated genes could also be a
sign of ligand activity.

When Ligand-Receptor pairs have both high DE and high activity, we can
consider them to be very good candidates in regulating the process of
interest, and we recommend testing these candidates for further
experimental validation.

# References

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<doi:10.1038/s41592-019-0667-5>

Guilliams et al. Spatial proteogenomics reveals distinct and
evolutionarily conserved hepatic macrophage niches. Cell (2022)
<doi:10.1016/j.cell.2021.12.018>
Inferring ligand-to-target signaling paths
================
Robin Browaeys
2019-01-17

<!-- github markdown built using 
rmarkdown::render("vignettes/ligand_target_signaling_path.Rmd", output_format = "github_document")
-->

### Infer signaling paths beween ligand(s) and target(s) of interest

To determine signaling paths between a ligand and target of interest, we
look at which transcription factors are best regulating the target genes
and are most closely downstream of the ligand (based on the weights of
the edges in the integrated ligand-signaling and gene regulatory
networks). Then, the shortest paths between these transcription factors
and the ligand of interests are determined and genes forming part in
this path are considered as important signaling mediators. Finally, we
look in our collected data source networks for all interactions between
the ligand, signaling mediators, transcription factors and target genes.
This allows to both prioritize signaling mediators and check which of
all collected data sources support the ligand-target predictions of
interest.

For this analysis, you need to define:

  - one or more ligands of interest
  - one or more target genes of interest

In this vignette, we will demonstrate how to infer signaling paths
between a CAF-ligand (CAF = cancer-associated fibroblast) of interest
and some of its top-predicted p-EMT target genes. The output of this
analysis can be easily imported into Cytoscape for exploration of the
networks.

First, we will load the necessary packages and networks to infer
signaling paths between ligand and target genes of interest.

``` r
library(nichenetr)
library(tidyverse)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
```

As example, we will infer signaling paths between the CAF-ligand TGFB3
and its top-predicted p-EMT target genes TGFBI, LAMC2 and
TNC.

``` r
ligands_all = "TGFB3" # this can be a list of multiple ligands if required
targets_all = c("TGFBI","LAMC2","TNC")

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

# To render the graph: uncomment following line of code
# DiagrammeR::render_graph(graph_min_max, layout = "tree")
```

![](tgfb3_targets_signaling_path.png)

We will now look which of the collected data sources support the
interactions in this
network.

``` r
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network) 
## # A tibble: 6 x 5
##   from  to    source            database       layer     
##   <chr> <chr> <chr>             <chr>          <chr>     
## 1 SMAD1 TGFBI regnetwork_source regnetwork     regulatory
## 2 SMAD1 TGFBI Remap_5           Remap          regulatory
## 3 SMAD2 LAMC2 harmonizome_CHEA  harmonizome_gr regulatory
## 4 SMAD2 TGFBI harmonizome_CHEA  harmonizome_gr regulatory
## 5 SMAD2 TNC   harmonizome_CHEA  harmonizome_gr regulatory
## 6 SMAD2 TNC   regnetwork_source regnetwork     regulatory
```

For information of all mentioned data sources in the source column (link
to the website of the database, etc), see [Data source
information](data_sources.xlsx)

### Export to Cytoscape

Export the following to e.g. Cytoscape for exploration of the networks

``` r
output_path = ""
write_output = FALSE # change to TRUE for writing output

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}
```
Perform NicheNet analysis starting from a Seurat object: step-by-step
analysis
================
Robin Browaeys
2019-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/seurat_steps.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform a basic NicheNet analysis
on a Seurat v3 object. Such a NicheNet analysis can help you to generate
hypotheses about an intercellular communication process of interest for
which you have single-cell gene expression data as a Seurat object.
Specifically, NicheNet can predict 1) which ligands from one or more
cell population(s) (“sender/niche”) are most likely to affect target
gene expression in an interacting cell population (“receiver/target”)
and 2) which specific target genes are affected by which of these
predicted ligands.

Because NicheNet studies how ligands affect gene expression in
putatively neighboring/interacting cells, you need to have data about
this effect in gene expression you want to study. So, there need to be
‘some kind of’ differential expression in a receiver cell population,
caused by ligands from one of more interacting sender cell populations.

In this vignette, we demonstrate the use of NicheNet on a Seurat Object.
The steps of the analysis we show here are also discussed in detail in
the main, basis, NicheNet vignette [NicheNet’s ligand activity analysis
on a gene set of interest: predict active ligands and their target
genes](ligand_activity_geneset.md):`vignette("ligand_activity_geneset", package="nichenetr")`.
Make sure you understand the different steps in a NicheNet analysis that
are described in that vignette before proceeding with this vignette and
performing a real NicheNet analysis on your data. This vignette
describes the different steps behind the wrapper functions that are
shown in [Perform NicheNet analysis starting from a Seurat
object](seurat_wrapper.md):`vignette("seurat_wrapper", package="nichenetr")`.
Following this vignette has the advantage that it allows users to adapt
specific steps of the pipeline to make them more appropriate for their
data.

As example expression data of interacting cells, we will use mouse
NICHE-seq data from Medaglia et al. to explore intercellular
communication in the T cell area in the inguinal lymph node before and
72 hours after lymphocytic choriomeningitis virus (LCMV) infection
(Medaglia et al. 2017). We will NicheNet to explore immune cell
crosstalk in response to this LCMV infection.

In this dataset, differential expression is observed between CD8 T cells
in steady-state and CD8 T cells after LCMV infection. NicheNet can be
applied to look at how several immune cell populations in the lymph node
(i.e., monocytes, dendritic cells, NK cells, B cells, CD4 T cells) can
regulate and induce these observed gene expression changes. NicheNet
will specifically prioritize ligands from these immune cells and their
target genes that change in expression upon LCMV infection.

The used NicheNet networks, ligand-target matrix and example expression
data of interacting cells can be downloaded from Zenodo. The NicheNet
networks and ligand-target matrix at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)
and the Seurat object of the processed NICHE-seq single-cell data at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3531889.svg)](https://doi.org/10.5281/zenodo.3531889).

# Prepare NicheNet analysis

## Load required packages, read in the Seurat object with processed expression data of interacting cells and NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks.

The NicheNet ligand-receptor network and weighted networks are necessary
to define and show possible ligand-receptor interactions between two
cell populations. The ligand-target matrix denotes the prior potential
that particular ligands might regulate the expression of particular
target genes. This matrix is necessary to prioritize possible
ligand-receptor interactions based on observed gene expression effects
(i.e. NicheNet’s ligand activity analysis) and infer affected target
genes of these prioritized ligands.

### Load Packages:

``` r
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
```

If you would use and load other packages, we recommend to load these 3
packages after the others.

### Read in the expression data of interacting cells:

The dataset used here is publicly available single-cell data from immune
cells in the T cell area of the inguinal lymph node. The data was
processed and aggregated by applying the Seurat alignment pipeline. The
Seurat object contains this aggregated data. Note that this should be a
Seurat v3 object and that gene should be named by their official
mouse/human gene symbol.

``` r
seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj@meta.data %>% head()
##         nGene nUMI orig.ident aggregate res.0.6 celltype nCount_RNA nFeature_RNA
## W380370   880 1611      LN_SS        SS       1    CD8 T       1607          876
## W380372   541  891      LN_SS        SS       0    CD4 T        885          536
## W380374   742 1229      LN_SS        SS       0    CD4 T       1223          737
## W380378   847 1546      LN_SS        SS       1    CD8 T       1537          838
## W380379   839 1606      LN_SS        SS       0    CD4 T       1603          836
## W380381   517  844      LN_SS        SS       0    CD4 T        840          513
```

Visualize which cell populations are present: CD4 T cells (including
regulatory T cells), CD8 T cells, B cells, NK cells, dendritic cells
(DCs) and inflammatory monocytes

``` r
seuratObj@meta.data$celltype %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application
## .
##     B CD4 T CD8 T    DC  Mono    NK  Treg 
##   382  2562  1645    18    90   131   199
DimPlot(seuratObj, reduction = "tsne")
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Visualize the data to see to which condition cells belong. The metadata
dataframe column that denotes the condition (steady-state or after LCMV
infection) is here called ‘aggregate.’

``` r
seuratObj@meta.data$aggregate %>% table()
## .
## LCMV   SS 
## 3886 1141
DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

### Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
## # A tibble: 6 x 4
##   from  to    source         database
##   <chr> <chr> <chr>          <chr>   
## 1 CXCL1 CXCR2 kegg_cytokines kegg    
## 2 CXCL2 CXCR2 kegg_cytokines kegg    
## 3 CXCL3 CXCR2 kegg_cytokines kegg    
## 4 CXCL5 CXCR2 kegg_cytokines kegg    
## 5 PPBP  CXCR2 kegg_cytokines kegg    
## 6 CXCL6 CXCR2 kegg_cytokines kegg

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  ABCC6  0.422 
## 2 A1BG  ACE2   0.101 
## 3 A1BG  ADAM10 0.0970
## 4 A1BG  AGO1   0.0525
## 5 A1BG  AKT1   0.0855
## 6 A1BG  ANXA7  0.457
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  A2M    0.0294
## 2 AAAS  GFAP   0.0290
## 3 AADAC CYP3A4 0.0422
## 4 AADAC IRF8   0.0275
## 5 AATF  ATM    0.0330
## 6 AATF  ATR    0.0355
```

Because the expression data is of mouse origin, we will convert the
NicheNet network gene symbols from human to mouse based on one-to-one
orthology:

``` r
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
```

# Perform the NicheNet analysis

In this case study, we want to apply NicheNet to predict which ligands
expressed by all immune cells in the T cell area of the lymph node are
most likely to have induced the differential expression in CD8 T cells
after LCMV infection.

As described in the main vignette, the pipeline of a basic NicheNet
analysis consist of the following steps:

## 1. Define a “sender/niche” cell population and a “receiver/target” cell population present in your expression data and determine which genes are expressed in both populations

In this case study, the receiver cell population is the ‘CD8 T’ cell
population, whereas the sender cell populations are ‘CD4 T,’ ‘Treg,’
‘Mono,’ ‘NK,’ ‘B’ and ‘DC.’ We will consider a gene to be expressed when
it is expressed in at least 10% of cells in one cluster.

``` r
## receiver
receiver = "CD8 T"
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
```

``` r
## sender
sender_celltypes = c("CD4 T","Treg", "Mono", "NK", "B", "DC")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
```

## 2. Define a gene set of interest: these are the genes in the “receiver/target” cell population that are potentially affected by ligands expressed by interacting cells (e.g. genes differentially expressed upon cell-cell interaction)

Here, the gene set of interest are the genes differentially expressed in
CD8 T cells after LCMV infection. The condition of interest is thus
‘LCMV,’ whereas the reference/steady-state condition is ‘SS.’ The notion
of conditions can be extracted from the metadata column ‘aggregate.’ The
method to calculate the differential expression is here the standard
Seurat Wilcoxon test, but this can be changed if necessary.

``` r
seurat_obj_receiver= subset(seuratObj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["aggregate"]])

condition_oi = "LCMV"
condition_reference = "SS" 
  
DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
```

## 3. Define a set of potential ligands: these are ligands that are expressed by the “sender/niche” cell population and bind a (putative) receptor expressed by the “receiver/target” population

Because we combined the expressed genes of each sender cell type, in
this example, we will perform one NicheNet analysis by pooling all
ligands from all cell types together. Later on during the interpretation
of the output, we will check which sender cell type expresses which
ligand.

``` r
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
```

## 4) Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)

``` r
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
## # A tibble: 44 x 5
##    test_ligand auroc  aupr pearson  rank
##    <chr>       <dbl> <dbl>   <dbl> <dbl>
##  1 Ebi3        0.638 0.234  0.197      1
##  2 Il15        0.582 0.163  0.0961     2
##  3 Crlf2       0.549 0.163  0.0758     3
##  4 App         0.499 0.141  0.0655     4
##  5 Tgfb1       0.494 0.140  0.0558     5
##  6 Ptprc       0.536 0.149  0.0554     6
##  7 H2-M3       0.525 0.157  0.0528     7
##  8 Icam1       0.543 0.142  0.0486     8
##  9 Cxcl10      0.531 0.141  0.0408     9
## 10 Adam17      0.517 0.137  0.0359    10
## # ... with 34 more rows
```

The different ligand activity measures (auroc, aupr, pearson correlation
coefficient) are a measure for how well a ligand can predict the
observed differentially expressed genes compared to the background of
expressed genes. In our validation study, we showed that the pearson
correlation coefficient between a ligand’s target predictions and the
observed transcriptional response was the most informative measure to
define ligand activity. Therefore, NicheNet ranks the ligands based on
their pearson correlation coefficient. This allows us to prioritize
ligands inducing the antiviral response in CD8 T cells.

The number of top-ranked ligands that are further used to predict active
target genes and construct an active ligand-receptor network is here 20.

``` r
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
```

These ligands are expressed by one or more of the input sender cells. To
see which cell population expresses which of these top-ranked ligands,
you can run the following:

``` r
DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

As you can see, most op the top-ranked ligands seem to be mainly
expressed by dendritic cells and monocytes.

## 5) Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis

### Active target gene inference

``` r
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
```

``` r
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

Note that not all ligands from the top 20 are present in this
ligand-target heatmap. The left-out ligands are ligands that don’t have
target genes with high enough regulatory potential scores. Therefore,
they did not survive the used cutoffs. To include them, you can be less
stringent in the used cutoffs.

### Receptors of top-ranked ligands

``` r
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
```

``` r
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

### Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases

``` r
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
```

``` r
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

## 6) Add log fold change information of ligands from sender cells

In some cases, it might be possible to also check upregulation of
ligands in sender cells. This can add a useful extra layer of
information next to the ligand activities defined by NicheNet, because
you can assume that some of the ligands inducing DE in receiver cells,
will be DE themselves in the sender cells.

Here this is possible: we will define the log fold change between LCMV
and steady-state in all sender cell types and visualize this as extra
information.

``` r
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "aggregate", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

``` r
# change colors a bit to make them more stand out
p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1))
p_ligand_lfc
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-40-2.png)<!-- -->

## 7) Summary visualizations of the NicheNet analysis

For example, you can make a combined heatmap of ligand activities,
ligand expression, ligand log fold change and the target genes of the
top-ranked ligands. The plots for the log fold change and target genes
were already made. Let’s now make the heatmap for ligand activities and
for expression.

``` r
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
```

``` r
# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
rotated_dotplot = DotPlot(seuratObj %>% subset(celltype %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
```

``` r
figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
    ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

# Remarks

1.  Top-ranked ligands and target genes shown here differ from the
    predictions shown in the respective case study in the NicheNet paper
    because a different definition of expressed genes was used.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-medaglia_spatial_2017" class="csl-entry">

Medaglia, Chiara, Amir Giladi, Liat Stoler-Barak, Marco De Giovanni,
Tomer Meir Salame, Adi Biram, Eyal David, et al. 2017. “Spatial
Reconstruction of Immune Niches by Combining Photoactivatable Reporters
and <span class="nocase">scRNA</span>-Seq.” *Science*, December,
eaao4277. <https://doi.org/10.1126/science.aao4277>.

</div>

</div>
FAQ NicheNet
================
Robin Browaeys
2020-05-20

<!-- github markdown built using
rmarkdown::render("vignettes/faq.Rmd", output_format = "github_document")
-->

This document tries to give an extensive answer to some questions we got
in the past months.

## When going to the vignettes, I see you require differential expression between two conditions. What if I only have steady-state data and want to use NicheNet? Can’t NicheNet be used to find ligand-receptor pairs in steady-state?

NicheNet is a tool that let you study how ligands affect gene expression
in putatively neighboring/interacting cells, and prioritize the most
important ligands based on their effect (in other words: prioritizing
expressed ligand-receptor interactions based on their observed target
genes). But to do this you need to have data about this effect in gene
expression you want to study. So, there need to be ‘some kind of’
differential expression in a receiver cell population, caused by ligands
from one of more interacting sender cell populations. Concretely, this
can be differential expression in one cell type between two conditions
caused by cell-cell interactions, or this can also be differential
expression between two cell types if for example these are a progenitor
and differentiated cell type where the differentiation is influenced by
the microenvironment. In that case you don’t necessarily need two
“conditions”.

If you would just be interested in the ligand-receptor interactions in
homeostatic conditions (and you only have homeostatic condition data),
you can always just infer ligand-receptor interactions based on their
expression level. This is the end goal in some other tools like
CellphoneDB (see also further question), or an intermediary step in the
NicheNet pipeline: (cf step 3 in the Seurat steps vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md>).

However, prioritizing expressed ligand-receptor interactions based on
observed target genes by applying NicheNet on steady-state data will be
harder because you don’t have a clear set of genes affected by the
intercellular communication process only. You could still do it, by
using the total set of genes in the genome as the background, and all
genes expressed in your cell type as gene set of interest (so you need
to change the geneset\_oi to all expressed genes, and the background to
all genes in the ligand-target matrix). But I don’t really recommend
that because not all expressed genes will be influenced by cell-cell
interactions, on the contrary most genes will probably be expressed
because they are part of general and cell-type specific gene expression
programs. And if you can’t discriminate between cell-intrinsic and
cell-extrinsic effects in your data, and you will still use NicheNet,
then you risk that ligands will be falsely linked to some
‘cell-intrinsic’ genes. And this will confound the analysis and lead
to irrelevant and noisy predictions.

So the take-away message is: NicheNet predicts upstream ligand-receptor
pairs that might regulate genes in your gene set of interest. The more
this gene set of interest consists of genes that you expect to be
regulated by the extracellular microenvironment, the better. If this
gene set would also consist of genes not necessarily influenced by the
environment (which is the case if using all expressed or cell-type
specific genes as gene set of interest), NicheNet can still provide good
predictions sometimes, but only if a substantial part of these genes
would affected by cell-cell interactions. Therefore, be cautious when
doing that and remember: garbage in, garbage out.

## As gene set of interest, I am using the Seurat cluster markers found via FindMarkers function. Is this a good idea?

As discussed in the question above, this is not a good idea in general.
If you use the marker genes of a (receiver) cluster as gene set of
interest, you assume that the cluster-specific genes are induced by
cell-cell interactions with other cell types. However, I think in most
cases, many cluster-specific genes reflect a cell-type specific gene
expression program that is part of the ‘cell-instrinsic’ gene expression
program, and not entirely influenced by environmental factors. Or in
other words, most genes of the cluster markers are probably not really
induced/regulated by the ligands from the interacting cells, and
therefore using NicheNet on this gene set is not ideal. Using NicheNet
in this case will increase the chance of falsely linking ligands to some
‘cell-intrinsic’ genes.

There might be some exceptions for which using cluster-specific genes as
gene set of interest might work though. Maybe you have 2 clusters of the
same cell type, but one cluster is located differently in the tissue,
and thus has different extracellular influences: then DE genes between
these clusters do not reflect differences in cell type, but in cell-cell
interactions affecting that cell type. In that case, you can use
NicheNet. Or another example possibility is: a progenitor cell cluster
differentiates into a another cell cluster under influence of cell-cell
interactions. Then you can use NicheNet on the cluster-specific genes of
that differentiated cell population to check influence of cell-cell
interactions on this differentiation process.

## If I compare the prioritized ligand-receptor interactions of NicheNet with those of CellphoneDB I find little overlap. How can that be explained?

A main reason is that both tools have a different goal in mind to
prioritize ligand-receptor interactions. The goal of CellphoneDB is to
find cell-type specific ligand-receptor pairs. So CellphoneDB is ideal
if you want to investigate the differential interaction potential of
different cell types. CellphoneDB first looks for all expressed
ligand-receptor pairs between interacting cells. Then it finds cell-type
specific ligand-receptor pairs by looking at the expression of both the
ligand and the receptor. The stronger and more specific
ligands/receptors are expressed in the sender/receiver cell populations,
the better it will be prioritized via a CellphoneDB analysis.

The goal of NicheNet is complementary to that of CellphoneDB. NicheNet
wants to find ligand-receptor pairs that are most likely to regulate
gene expression in the receiver cell. So it looks for ligand-receptor
pairs for which evidence for signaling interactions exist in the data.
In contrast to CellphoneDB this might give more functional information
(you have an idea about the downstream targets of a ligand-receptor
interaction), but also some clues about which interactions might really
be active. This because expression of the ligand and the receptor at RNA
level (as found by CellphoneDB), does not necessarily mean that they
interact in reality. If downstream signaling effects of such an
interaction are observed, e.g. via NicheNet, there is a bit more
evidence that this interaction might be happening in reality. So just
like CellphoneDB, NicheNet starts to look for all expressed
ligand-receptor pairs between the interacting cells of interest. But
instead of prioritizing these by looking at expression strength,
NicheNet will prioritize them based on observed target genes in the
receiver cell. A possible disadvantage of this approach is that some
ligand-receptor pairs with relatively low expression can be returned.
Therefore we also recommend to check expression of the ligands and their
receptors after the NicheNet prioritization. To give an example: If a
ligand ranked 7th out of 200 would be much stronger expressed than a
ligand that would be 1st or 2nd ranked according to the ligand activity,
this more strongly expressed candidate ligand might be more
interesting\!

So when comparing NicheNet to CellphoneDB output following differences
can be expected to be present. Ligad-receptor pairs picked up by
NicheNet but not by CellphoneDB might be or generally expressed or
rather lowly expressed, but have some ‘signaling evidence’. Pairs picked
by CellphoneDB but not by NicheNet will be strongly and cell-type
specific expressed, but there is no evidence based on prior knowledge on
signaling pathways that these pairs will actually be functional. Of
course, a lack of evidence does not mean there is no signaling is going
on.

## When I check the expression of some top ligands according to NicheNet, I see that some of these ligands and/or their receptor are really lowly expressed? How is that possible? And does NicheNet not take into account expression data of the cells?

The prioritization of ligands by NicheNet (ligand activity analysis)
will only occur based on enrichment of their target genes in the set of
genes that are differentially expressed in the receiver cell. So there
is no prioritization based on the strength of expression of the ligand
in the sender cell or strength of expression of the receptor(s) in the
receiver cell. Expression in sender cells is only used to determine
which ligands are expressed in a sender cell, and expression in receiver
cells is used to determine which receptors are expressed in the receiver
cell. The default definition of ‘being expressed’ is that a gene should
be expressed in 10% of cells in the cluster of interest. This is not so
high (you can put a more stringent cutoff if you want), resulting in the
possible outcome that a ligand, top-ranked according to the enrichment
of its target genes, is actually not very highly expressed. So what you
observe, can be expected based on how NicheNet prioritizes ligands.

In the current version of NicheNet, expression strength is thus not
directly included because we find it hard to formalize the tradeoff
between ligand activity and expression. But because expression level is
important, we recommend to check expression of the ligands and their
receptors after the NicheNet prioritization. To give an example: If a
ligand ranked 7th out of 200 would be much stronger expressed than a
ligand that would be 1st or 2nd ranked according to the ligand activity,
this more strongly expressed candidate ligand might be more
interesting\!

## I observe that some top ligands according to NicheNet to be inducing in my condition of interest compared to control are higher expressed in control than condition of interest. How is this possible?

The ligand prioritization is currently done without taking into account
to expression value of the ligand, the only condition is that the ligand
is expressed. This means that the algorithm does not consider whether a
ligand is more strongly expressed in the case vs control, as long as it
is expressed at a sufficient level. Of course, this might be interesting
information to use for further prioritization after the NicheNet ligand
activity analysis\! Therefore we recommend checking the expression of
your ligands after NicheNet prioritization. The most interesting hits
might be the ones where ligand and/or receptor is upregulated in the
condition of interest compared to the control. The opposite pattern
could also be interesting if you could consider following hypothesis.
Some ligands might downregulate many genes (instead of, or in addition
to, upregulate many genes). This means that if ligand X is active in the
control condition, it can downregulate genes there. When ligand X itself
is repressed in the case condition, and has a lower expression, its
downregulatory effect might disappear, leading to upregulation of genes
that are a downregulated target gene of the ligand X, and high ligand
activity of ligand X according to NicheNet.

If you would only want to consider ligands upregulated in the
case-vs-control condition, you can always change the pipeline by first
performing a differential expression analysis on your sender cells, and
considering only upregulated ligands, and not all expressed ligands as
‘potential ligands’ for the NicheNet analysis. But I would recommend
only filtering on upregulation afterwards, because ligands can be
differentially active due to other processes than upregulation at the
RNA level.

## Can I use NicheNet if I want to find which potential ligands were important in cells of interest, even though I don’t have expression data of the possible interacting/sender cells?

It is perfectly possible to apply NicheNet in case you don’t have a
‘sender cell type’. In that case, your analysis is more of an
‘upstream regulator analysis’ than an ‘intercellular communication
analysis’. You can do this by looking at all ligands in the NicheNet
database (or all ligands for which a potential receptor is expressed)
instead of only looking ligands that are expressed by a specific sender
cell type. Only a small change in the code is required for this by using
`ligands` instead of `expressed_ligands` in the basic and Seurat steps
vignette for defining `potential_ligands`. If using the Seurat wrapper
function you can set the parameter function of sender to `sender =
“undefined”`.

You just need to be sure that you can use the expression data of the
responder/receiver cell type to define a set of genes that are likely
affected by the upstream regulatory intercellular signaling process.

In case you don’t have data about the sender cell type(s) of interest,
but you know there is some publicly available expression data available
of that cell type or the tissue of interest, it might be helpful to use
this dataset as proxy for the expression of the ligands.

## Is the NicheNet prior model based on human or mouse data or is there a separate model for human and mouse?

The model is primarily build based on human data sources. However, we
also included some mouse data sources into the general model (at least
for genes that have a mouse one-to-one ortholog), but these were
strongly in the minority. So the NicheNet model is tailored more for
human than mouse, and the mouse model that we use in some of the
vignettes is just obtained by converting human gene symbols to their
mouse one-to-one orthologs.

Because this model is indeed tailored more for human, you can expect
that it would perform a bit better for human than mouse, but I don’t
think the difference is very large. In our paper, we evaluated the
ligand-target predictions based on \>100 ligand treatment datasets of
both human and mouse origin, and there is no indication on these
datasets that human would work better than mouse.

However, there is a disadvantage: for some mouse-specific genes, we
don’t have information in the NicheNet ligand-target model. This
implies that NicheNet might miss some patterns if well-known mouse-only
target genes would be DE in the dataset you are working on. A solution
to this would be to construct a mouse model with additional mouse data
sources yourself. This is some work, but is not so difficult. In the
vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_construction.md>,
we show how to build a model yourself.

## Can I use NicheNet on data other than mouse or human data?

Yes, there are two options here.

First, you can explore the possibility of constructing an
organism-specific model yourself with only data sources of your organism
of interest. This is some work, but is not so difficult. In the vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_construction.md>,
we show how to build a model yourself. The main thing to consider
though, is that you should have enough primary data sources of the
organism to do this.

The second option would be to use the current human NicheNet model and
convert the human gene symbols to their one-to-one orthologs of the
organism of interest. The main thing to consider here is that the
organism of interest should not be too dissimilar from human.

So, you should make the tradeoff between the homology between species,
and the number and quality of data sources that are available to build a
new species-specific model. For example, for primates, you could use the
human model and convert gene symbols, but for Drosophila making a new
model from own data sources is probably more appropriate.

## I decided that I want to construct a model with my own data sources. Do I really need to optimize the data source weights?

As shown in our paper, data source weight optimization does improve the
performance the model a bit, but this effect is not large. A model
without data source weight optimization also seemed to work well.
Because the optimization is a difficult and time-consuming procedure, we
don’t see a problem in not optimizing your data source weights, but only
if you are confident that your data sources are of high quality and not
too noisy\!

## Can NicheNet say which cell populations interact with each other and which don’t?

No, NicheNet can’t give you a clear-cut answer to that question. But, it
can suggest this. If for many ligands from a certain cell types, many
target genes are DE in the receiver cell type, this might suggest that
some intercellular signaling processes are going on between both cell
types. In this vignette:
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/target_prediction_evaluation_geneset.md>,
we show how you can calculate what fraction of genes DE in the receiver
cell type, might be a target of the top prioritized ligands of NicheNet.
This way, we found that in one of the case studies described in the
paper 50% of the DE genes were a strongly predicted target gene of the
prioritized ligands. Such a high fraction suggests that intercellular
signaling might indeed be going on, but this should of course be still
validated.

## Can I use NicheNet as a gene set based approach, e.g. to extract a gene set consisting of target genes of ligand X?

Yes, you can use the functions `extract_top_n_targets()` and
`extract_top_fraction_targets()` for this. In the near future, we will
also add a vignette that can be used to score individual cells based on
their expression of the top n target genes of a ligand of interest.

## How can I generate a double layer ligand-receptor-target circos plot as shown in the Immunity paper (<https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1>)?

As far as we know, there is no straightforward way to directly make this
kind of “ligand-receptor-target” circos plot via the circlize R package.
Therefore we made this “ligand-receptor-target” circos plot by making
first two separate circos plots: the ligand-target and ligand-receptor
circos plot. We made sure that the ligand-receptor circos plot was
bigger than the ligand-target plot. Then we overlayed them in Inkscape
(with the center of the two circles at the same location) and removed
everything from the ligand-receptor circos plot except the outer
receptor layer. Making the ligand-receptor circos plot is very similar
to making the ligand-target circos plot, as shown in the circos plot
vignette.

If you would want to split up target genes and receptors in different
groups according to signaling pathway (as done in that paper), then you
first need to define these groups in a specific data frame in advance
(cf what is shown for ligands in the `ligand_type_indication_df`in the
vignette). When you then want to overlay receptors in this case, you
need to make sure that the ligand-receptor weights of receptors in one
group are proportional to the ligand-target weights of the targets in
that group (to generate the nice overlay effect). So in that case, the
ligand-receptor weights are proportional to the ‘underlying’
ligand-target regulatory potential scores and not reflective of prior
information supporting the specific ligand-receptor interaction (as
shown in the current vignette for ligand-receptor circos plots).

We see that this is indeed a cumbersome and far from ideal situation.
Ideally we would be working on a more straightforward solution in the
future, but we can’t promise anything. This is not so high priority for
us because we think that the combined ligand-activity, ligand-receptor
and ligand-target heatmaps are (at least) equally informative.

## Can I use my own ligand-receptor network instead of the one included in the NicheNet model?

Yes, you definitely can\! In this vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_construction.md>,
you can see how to incorporate your new data source with new
ligand-receptor pairs and build a customized NicheNet prior model. I
would suggest you just add your new ligand-receptor data source to the
existing ones (or alternatively: remove the existing ones), and give
your new data source as weight 1, and use the optimized weights for the
other data sources.

## Will you update the data sources behind NicheNet and the NicheNet model?

We are indeed planning to update the prior model once by incorporating
some new and updated databases, but we can’t pin a date on this now.

## When interpreting the output ligand-target matrix, I see that not all my genes of my gene set of interest are in the final ligand-target matrix? Why is this?

The reason for this is that the ligand-target matrix only shows gene
that are a top predicted target of at least one of the top ligands. You
can tweak some of the parameters of the function to be more permissive
if you want to include more genes (by allowing lower regulatory
potential scores / setting the n of top n targets higher).

## Can I use NicheNet on bulk RNAseq data as well, or only on single-cell data?

Yes, you can if you are working on cell population – sorted bulk data.
Because of the better sequencing depth, this might even work better (of
course depending on the setting and whether you have clean data of the
cell populations of interest). We already applied NicheNet on bulk
RNA-seq data, as shown in this paper:
<https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1> .

We did this analysis by following the steps explained in the basic
vignette:
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md>;
Some adaptations to the basic vignette are needed, though. For example,
you need another way for defining the set of expressed genes in sender
and receiver, and another way of doing the DE analysis for defining the
gene set of interest. To define which genes are expressed you could use
`filterByExpr` from the package `edgeR` ; or check the distribution of
expression values and choosing a cutoff based on this distribution. If
you doubt about choosing the cutoff, I would recommend to consider more
genes expressed than removing some genes from the analysis. The number
of expressed genes differs per cell type and organism, but should be in
the range of 10000-15000 for bulk rna-seq data (so not 1000 or 2000
genes).

In the Materials and Methods of this paper, you can see as example what
we did to determine expressed genes and perform the DE analysis.

## Does NicheNet return all significant ligand-receptor interactions?

No, you should realize that NicheNet does not give a clear-cut answer
about which ligand-receptor pairs are significant and which not. First,
in the NicheNet pipeline, we determine all possible expressed
ligand-receptor pairs. Then, a ligand activity analysis is performed to
rank ligands based on their activity (this accords to the enrichment of
strongly predicted target genes of a ligand in the gene set of interest
in the receiver cell). Based on that, the top n (default 20) ligands are
selected for further analyses and visualizations.

## If I look at the ligand-receptor network, I see some low interaction scores between highly expressed ligand-receptor interactions, or some high scores between lowly expressed pairs? What is the reason for this?

The ligand-receptor interaction scores shown in the ligand-receptor
matrix/heatmap are a proxy for the confidence that this ligand interact
with that receptor based on prior knowledge. The ligand-receptor network
at the basis of NicheNet is the result of integrating multiple single
ligand-receptor data sources. The score you see on this plot, the prior
ligand-receptor interaction potential, is a proxy of the number and
‘quality’ of the ligand-receptor and PPI data sources that report a
certain ligand-receptor interaction. If a ligand-receptor interaction is
described in many, confident databases, this score will be higher than
if the interaction is predicted via just one protein-protein interaction
database.

This score is thus solely based on prior knowledge and not based on the
expression in your dataset\! The only way expression is used, is to
define whether a ligand/receptor is expressed yes or no. So
ligands/receptor that are not expressed at all won’t appear in this
heatmap, but it is indeed possible that some lowly expressed
ligand-receptor pairs will be shown here.

In the near future, we will provide code to show how to make this
heatmap with scores that reflect the expression strength of both ligand
and receptor (product of expression / average of expression).

## In the Seurat vignette, you show some code to define bona fide ligand-receptor interactions. But I see that some ligand-receptor pairs that are not considered bona fide have higher prior interaction scores than some bona fide pairs. How is this possible if the prior interaction score is based on the confidence of the prior information? Should bona fide pairs not have the highest scores?

We can indeed see why this confusing\! We categorize bona-fide ligands
based on whether they were coming from a validated curated
ligand-receptor database (like KEGG, guide2Pharmacology, …) or not. This
because we have also ligand-receptor interactions not present in these
databases, but predicted based on protein-protein interaction networks:
e.g. if based on annotation we know that gene X encodes a ligand, and
gene Y a receptor, we predict a ligand-receptor interaction between X
and Y if there is a PPI database (but not ligand-receptor database)
supporting X-Y interaction. Because they are predicted as
ligand-receptor interaction (and not part of a curated ligand-receptor
database), they are not bona fide. But they can have higher scores than
some bona fide ones, if there are many confident PPI and pathway
databases reporting this interaction. So basically: a bona fide
ligand-receptor interaction is known as ligand-receptor interaction, a
non-bona fide interaction is a PPI (with much evidence if high scores)
but not yet documented as LR interaction in a curated database. Focusing
on bona fide LR interactions will give you well-validated pairs, whereas
including non-bona fide ones can be more novel, but are less confident.

## Although NicheNet already prioritizes ligand-receptor pairs strongly, I still find there are too many possible pairs to experimentally validate. What types of information would you recommend to consider for even further prioritization?

I think you could take many factors into consideration. First, you can
consider the NicheNet ligand activity (as you already did). But instead
of having a strong preference for the 1st ranked ligand versus the 10th
ligand, I suggest also taking into account the expression level of both
ligand and receptor. For this you can look at the expression value, the
fraction of cells expressing the ligand/receptor and whether the
ligand/receptor is cell-type specific or not. In the ideal case, you
would have case-vs-control data and you could also look at
ligands/receptor pairs for which the ligand and/or receptor is
upregulated in the case vs control condition. Next to this, you could
check whether the ligand-receptor interaction is ‘bona fide’ or not (cf
previous question) and whether the ligand-receptor prior interaction
potential weight is high or not. But I suggest giving the highest weight
to the ligand activity and expression/upregulation information\!

## How should I interpret the pearson correlation values of the ligand activity? What does it mean when I have negative values?

In the ligand\_pearson\_matrix, you find the pearson correlation values
between NicheNet prior ligand-target predictions and the gene set of
interest. A high score means that top predicted target genes of a ligand
are enriched in the gene set of interest compared to the background. We
would say that a good Pearson score for top ligands is around 0.10 and
higher. If you have low non-zero scores (such as 0.04), it might still
be that the ranking of ligands is valuable, although this ranking will
only be based on a few top predicted target genes of a ligand that are
in the gene set of interest. Scores around zero (such as 0.008) mean
that top predicted targets of ligands are not enriched in the gene set
of interest compared to a background, thus that ranking of ligands is
not valuable in that case. These ligands will typically have an AUROC
around 0.50. This means that NicheNet does not find evidence that the
sender cell is regulating gene expression changes in the receiver cell
(although the lack of evidence does not demonstrate the lack of).

If these values are negative, this means that there is an enrichment of
the top predicted target genes in the background compared to the gene
set of interest. When these values are negative, but very close to zero,
then there is nothing to worry about. But when these values are really
negative, this can indicate that something went wrong in the analysis,
probably in the definition of background and gene set of interest.
Because the gene set of interest is typically much smaller than the
background, a stronger anti-correlation score for a ligand would only
occur when the gene set of interest would almost entirely consist of the
genes with the lowest regulatory potential scores for a ligand. And this
would only rarely occur. One possible explanation for finding strongly
negative scores, could be that your background gene set is too small,
and thus not a good representation of the genomic background. Your gene
set of interest should never be larger than the background. We recommend
to have at least around 5000-6000 genes in your background. If you would
still have issues, or not enough expressed genes to have a large
background, you could use all genes present in the ligand-target matrix
as genomic background.

## If I look at the combined ligand-activity-target heatmap, I see that some ligands seem to have higher activity, but less target genes than other ligands. How is this possible if the ligand activity is based on the observed target genes in the gene set of interest?

The ranking of the ligands is based on how well top-ranked target genes
of a ligand (based on regulatory potential) are enriched in the gene set
of interest (= Pearson correlation between whether a gene belongs to the
gene set of interest and the regulatory potential scores of target genes
of a ligand). For highly ranked ligands, this means that many top-ranked
target genes are in the gene set of interest compared to the background.
Note that the Pearson correlation to rank ligands is calculated for each
ligand separately. As a result, this correlation only depends on the
ranking of target genes for one ligand ( ‘relative’ ligand-target
regulatory potential scores for that ligand, without looking at other
ligands).

The regulatory potential scores of ligand-target links that are shown in
the heatmap visualize the ‘absolute’ ligand-target regulatory potential.
This is a confidence score related to how many data sources confirm the
regulatory interaction between a ligand and a target. A higher score
means that there is more evidence for the specific ligand-target
interaction. It is important to note that in this heatmap, you only show
the genes that are part of your gene set as possible target genes.

The reason why you can see a contradiction between the Pearson
correlation and regulatory potential scores is thus the following:
Ligand X is ranked higher than Ligand Y, because its top-ranked target
genes are more enriched than the target genes of Ligand Y, although the
absolute regulatory potential scores of the Ligand Y-target links are
higher (this can be the case for very well-studied ligands like TNF,
IFNG, IL1B, …) . Ligand Y targets are however probably less enriched
because there will might be more high scoring target genes of Ligand Y
that are not in the gene set of interest compared to Ligand X targets.
That can explain why there is stronger enrichment for Ligand X targets,
although there is more information for Ligand Y.

## Can NicheNet also calculate ligand-receptor interactions that specifically downregulates target genes? Or can NicheNet distinguish between target genes that are upregulated or downregulated due to a ligand-receptor interaction?

NicheNet can indeed work with both upregulated and downregulated target
genes. In default cases, I would recommend considering both upregulated
and downregulated genes as potential target genes in the gene set of
interest, because some ligand-receptor interactions will indeed lead to
both up- and downregulation of genes as you mention. But it is perfectly
possible to look at down- or upregulated genes only. Depending on the
research question, this could sometimes be interesting. It could also be
interesting to find which ligands are more responsible for upregulation
and which more for downregulation.

This might be interesting because the NicheNet prioritization algorithm
does not directly distinguish between up- and downregulated genes. This
has an important consequence to keep in mind: if there are more
upregulated genes than downregulated genes, ligands that only upregulate
genes will be more likely active than downregulatory genes.

NicheNet cannot directly distinguish between up- and downregulated genes
because the regulatory potential accords to the
evidence/potential/probability that a ligand might regulate the
expression of a target gene. And regulate can mean both upregulation and
downregulation. It’s not that NicheNet knows a priori which target genes
will be upregulated and which downregulated. This is because of the
following reason: during construction of the model (so the calculate the
prior regulatory potential scores), we considered all retrieved
interactions as active because most databases at the basis of NicheNet
don’t provide information about whether an interaction is inducing or
repressing. So the only information we have in the final model is: there
is no evidence for the regulation of target X by ligand Y (reg potential
of 0) or there is evidence for the regulation of target X by ligand Y
(reg potential \> 0; with more evidence leading to higher scores; but a
repressive regulatory interaction leading to downregulation will also
have a score higher than 0).

What you can do as a user is visualizing the expression of the target
genes (in a combined heatmap for example, as shown in the basic
vignette). This will allow you to check which specific genes are up- and
downregulated.

## If I have two main cell types, divided in two subpopulations based on the condition (case-vs-control): what do I need to use as receiver and sender cell type of interest in e.g. the Seurat vignette?

I recommend using all cells of the receiver cell type of interest as
receiver (so including both conditions). For the senders, it might
sometimes be better to use only the cells in the case condition.
However, I think it is perfectly fine to include all sender cells as
well (thus from both conditions) and check after the NicheNet analyses
which ligands might be upregulated in the case condition compared to
control.

## My question is not in these list? What should I do now?

First, you can check the open and closed issues
(<https://github.com/saeyslab/nichenetr/issues>) of this package on
github to see whether your question might be addressed in one of these.
If not, don’t hesitate to open a new issue. If you would prefer to keep
the discussion private, you can also send me an email
(<robin.browaeys@ugent.be>), but I prefer that you open an issue so
other users can learn from it as well\!
Circos plot visualization to show active ligand-target links between
interacting cells
================
Robin Browaeys
3-7-2019

<!-- github markdown built using 
rmarkdown::render("vignettes/circos.Rmd", output_format = "github_document")
-->

This vignette shows how NicheNet can be used to predict active
ligand-target links between multiple interacting cells and how you can
make a circos plot to summarize the top-predicted links (via the
circlize package). This vignette starts in the same way as the main,
basis, NicheNet vignette [NicheNet’s ligand activity analysis on a gene
set of interest: predict active ligands and their target
genes](ligand_activity_geneset.md):`vignette("ligand_activity_geneset",
package="nichenetr")`. Make sure you understand the different steps
described in that vignette before proceeding with this vignette. In
contrast to the basic vignette, we will look communication between
multiple cell types. More specifically, we will predict which ligands
expressed by both CAFs and endothelial cells can induce the p-EMT
program in neighboring malignant cells (See Puram et al. 2017).

### Load packages required for this vignette

``` r
library(nichenetr)
library(tidyverse)
library(circlize)
```

### Read in expression data of interacting cells

First, we will read in the publicly available single-cell data from
CAFs, endothelial cells and malignant cells from HNSCC
tumors.

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

Secondly, we will determine which genes are expressed in CAFs,
endothelial and malignant cells from high quality primary tumors.
Therefore, we wil not consider cells from tumor samples of less quality
or from lymph node metastases. To determine expressed genes, we use the
definition used by of Puram et
al.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% pull(cell)
endothelial_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "Endothelial") %>% pull(cell)
malignant_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_endothelial = expression[endothelial_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
```

### Load the ligand-target model we want to use

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

### Load the gene set of interest and background of genes

As gene set of interest, we consider the genes of which the expression
is possibly affected due to communication with other cells.

Because we here want to investigate how CAFs and endothelial cells
regulate the expression of p-EMT genes in malignant cells, we will use
the p-EMT gene set defined by Puram et al. as gene set of interset and
use all genes expressed in malignant cells as background of
genes.

``` r
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(pemt_geneset)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"

background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"
```

### Perform NicheNet’s ligand activity analysis on the gene set of interest

In a first step, we will define a set of potentially active ligands. As
potentially active ligands, we will use ligands that are 1) expressed by
CAFs and/or endothelial cells and 2) can bind a (putative) receptor
expressed by malignant cells. Putative ligand-receptor links were
gathered from NicheNet’s ligand-receptor data sources.

Note that we combine the ligands from CAFs and endothelial cells in one
ligand activity analysis now. Later on, we will look which of the
top-ranked ligands is mainly expressed by which of both cell
types.

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands_CAFs = intersect(ligands,expressed_genes_CAFs)
expressed_ligands_endothelial = intersect(ligands,expressed_genes_endothelial)
expressed_ligands = union(expressed_ligands_CAFs, expressed_genes_endothelial)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "IL15"    "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"
```

Now perform the ligand activity analysis: infer how well NicheNet’s
ligand-target potential scores can predict whether a gene belongs to the
p-EMT program or
not.

``` r
ligand_activities = predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

Now, we want to rank the ligands based on their ligand activity. In our
validation study, we showed that the pearson correlation between a
ligand’s target predictions and the observed transcriptional response
was the most informative measure to define ligand activity. Therefore,
we will rank the ligands based on their pearson correlation coefficient.

``` r
ligand_activities %>% arrange(-pearson) 
## # A tibble: 154 x 4
##    test_ligand auroc   aupr pearson
##    <chr>       <dbl>  <dbl>   <dbl>
##  1 PTHLH       0.667 0.0720   0.128
##  2 EDN1        0.682 0.0586   0.126
##  3 CXCL12      0.680 0.0507   0.123
##  4 AGT         0.676 0.0581   0.120
##  5 TGFB3       0.689 0.0454   0.117
##  6 IL6         0.693 0.0510   0.115
##  7 INHBA       0.695 0.0502   0.113
##  8 ADAM17      0.672 0.0526   0.113
##  9 TNC         0.700 0.0444   0.109
## 10 VWF         0.685 0.0490   0.109
## # … with 144 more rows
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)
## [1] "PTHLH"  "EDN1"   "CXCL12" "AGT"    "TGFB3"  "IL6"
```

We see here that the top-ranked ligands can predict the p-EMT genes
reasonably, this implies that ranking of the ligands might be accurate
as shown in our study. However, it is possible that for some gene sets,
the target gene prediction performance of the top-ranked ligands would
not be much better than random prediction. In that case, prioritization
of ligands will be less trustworthy.

Determine now which prioritized ligands are expressed by CAFs and or
endothelial cells

``` r
best_upstream_ligands %>% intersect(expressed_ligands_CAFs) 
##  [1] "PTHLH"  "CXCL12" "AGT"    "TGFB3"  "IL6"    "INHBA"  "ADAM17" "TNC"    "CTGF"   "FN1"    "BMP5"   "IL24"  
## [13] "CXCL11" "MMP9"   "COL4A1" "PSEN1"  "CXCL9"
best_upstream_ligands %>% intersect(expressed_ligands_endothelial)
##  [1] "EDN1"   "CXCL12" "IL6"    "ADAM17" "VWF"    "CTGF"   "FN1"    "SPP1"   "CXCL11" "COL4A1" "PSEN1"  "CXCL9"

# lot of overlap between both cell types in terms of expressed ligands
# therefore, determine which ligands are more strongly expressed in which of the two
ligand_expression_tbl = tibble(
  ligand = best_upstream_ligands, 
  CAF = expression[CAF_ids,best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
  endothelial = expression[endothelial_ids,best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}))

CAF_specific_ligands = ligand_expression_tbl %>% filter(CAF > endothelial + 2) %>% pull(ligand)
endothelial_specific_ligands = ligand_expression_tbl %>% filter(endothelial > CAF + 2) %>% pull(ligand)
general_ligands = setdiff(best_upstream_ligands,c(CAF_specific_ligands,endothelial_specific_ligands))

ligand_type_indication_df = tibble(
  ligand_type = c(rep("CAF-specific", times = CAF_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length()),
                  rep("Endothelial-specific", times = endothelial_specific_ligands %>% length())),
  ligand = c(CAF_specific_ligands, general_ligands, endothelial_specific_ligands))
```

### Infer target genes of top-ranked ligands and visualize in a circos plot

Now we will show how you can look at the regulatory potential scores
between ligands and target genes of interest. In this case, we will look
at links between top-ranked p-EMT-regulating ligands and p-EMT genes. In
this example, inferred target genes should belong to the p-EMT gene set
and to the 250 most strongly predicted targets of at least one of the
selected top-ranked ligands (the top 250 targets according to the
general prior model, so not the top 250 targets for this dataset).

Get first the active ligand-target links by looking which of the p-EMT
genes are among the top-predicted target genes for the prioritized
ligands:

``` r
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = pemt_geneset, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "p_emt") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
```

To avoid making a circos plots with too many ligand-target links, we
will show only links with a weight higher than a predefined cutoff:
links belonging to the 66% of lowest scores were removed. Not that this
cutoffs and other cutoffs used for this visualization can be changed
according to the user’s
needs.

``` r
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
  
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
```

Prepare the circos visualization: give each segment of ligands and
targets a specific color and order

``` r
grid_col_ligand =c("General" = "lawngreen",
            "CAF-specific" = "royalblue",
            "Endothelial-specific" = "gold")
grid_col_target =c(
            "p_emt" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
```

Prepare the circos visualization: order ligands and targets

``` r
target_order = circos_links$target %>% unique()
ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)
```

Prepare the circos visualization: define the gaps between the different
segments

``` r
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CAF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Endothelial-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "p_emt") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
  )
```

Render the circos plot (all links same transparancy). Only the widths of
the blocks that indicate each target gene is proportional the
ligand-target regulatory potential (\~prior knowledge supporting the
regulatory interaction).

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
```

![](circos_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
circos.clear()
```

Render the circos plot (degree of transparancy determined by the
regulatory potential value of a ligand-target interaction)

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
```

![](circos_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
circos.clear()
```

Save circos plot to an svg file

``` r
svg("ligand_target_circos.svg", width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()
## png 
##   2
```

### Visualize ligand-receptor interactions of the prioritized ligands in a circos plot

``` r
# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "p_emt_receptor") %>% inner_join(ligand_type_indication_df)
```

``` r
grid_col_ligand =c("General" = "lawngreen",
            "CAF-specific" = "royalblue",
            "Endothelial-specific" = "gold")
grid_col_receptor =c(
            "p_emt_receptor" = "darkred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
```

Prepare the circos visualization: order ligands and receptors

``` r
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)
```

Prepare the circos visualization: define the gaps between the different
segments

``` r
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  # width_ligand_receptor,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CAF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Endothelial-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_receptor,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "p_emt_receptor") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
  )
```

Render the circos plot (all links same transparancy). Only the widths of
the blocks that indicate each receptor is proportional the
ligand-receptor interaction weight (\~prior knowledge supporting the
interaction).

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
```

![](circos_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
circos.clear()
```

Render the circos plot (degree of transparancy determined by the prior
interaction weight of the ligand-receptor interaction - just as the
widths of the blocks indicating each receptor)

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
```

![](circos_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
circos.clear()
```

Save circos plot to an svg file

``` r
svg("ligand_receptor_circos.svg", width = 15, height = 15)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()
dev.off()
## png 
##   2
```

### Remark on making a ligand-receptor-target circos plot

In the paper of Bonnardel, T’Jonck et al. [Stellate Cells, Hepatocytes,
and Endothelial Cells Imprint the Kupffer Cell Identity on Monocytes
Colonizing the Liver Macrophage
Niche](https://www.cell.com/immunity/fulltext/S1074-7613\(19\)30368-1),
we showed in Fig. 6B a ligand-receptor-target circos plot to visualize
the main NicheNet predictions. This “ligand-receptor-target” circos plot
was made by making first two separate circos plots: the ligand-target
and ligand-receptor circos plot. Then these circos plots were overlayed
in Inkscape (with the center of the two circles at the same location and
the ligand-receptor circos plot bigger than the ligand-target one). To
generate the combined circos plot as shown ni Fig. 6B, we then manually
removed all elements of the ligand-receptor circos plot except the outer
receptor layer. In the near future, we will be working on a solution to
generate this ligand-receptor-target circos plot in a fully automated
manner.

If you would want to split up target genes and receptors in different
groups according to signaling pathway (as done in mentioned paper), then
you first need to define these groups in a specific data frame in
advance (cf what is shown for ligands in the
`ligand_type_indication_df`in the vignette). When you then want to
overlay receptors in this case, you need to make sure that the
ligand-receptor weights of receptors in one group are proportional to
the ligand-target weights of the targets in that group (to generate the
nice overlay effect). So in that case, the ligand-receptor weights are
proportional to the ‘underlying’ ligand-target regulatory potential
scores and not reflective of prior information supporting the specific
ligand-receptor interaction (as shown in this vignette for
ligand-receptor circos plots).

### References

Bonnardel et al., 2019, Immunity 51, 1–17, [Stellate Cells, Hepatocytes,
and Endothelial Cells Imprint the Kupffer Cell Identity on Monocytes
Colonizing the Liver Macrophage
Niche](https://doi.org/10.1016/j.immuni.2019.08.017)

<div id="refs" class="references">

<div id="ref-puram_single-cell_2017">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
Perform NicheNet analysis starting from a Seurat object
================
Robin Browaeys
2019-11-08

<!-- github markdown built using 
rmarkdown::render("vignettes/seurat_wrapper.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform a basic NicheNet analysis
on a Seurat v3 object. Such a NicheNet analysis can help you to generate
hypotheses about an intercellular communication process of interest for
which you have single-cell gene expression data as a Seurat object.
Specifically, NicheNet can predict 1) which ligands from one or more
cell population(s) (“sender/niche”) are most likely to affect target
gene expression in an interacting cell population (“receiver/target”)
and 2) which specific target genes are affected by which of these
predicted ligands.

Because NicheNet studies how ligands affect gene expression in
putatively neighboring/interacting cells, you need to have data about
this effect in gene expression you want to study. So, there need to be
‘some kind of’ differential expression in a receiver cell population,
caused by ligands from one of more interacting sender cell populations.

In this vignette, we demonstrate the use of NicheNet on a Seurat Object.
The wrapper function we will show consists of the same different steps
that are discussed in detail in the main, basis, NicheNet vignette
[NicheNet’s ligand activity analysis on a gene set of interest: predict
active ligands and their target
genes](ligand_activity_geneset.md):`vignette("ligand_activity_geneset", package="nichenetr")`.
Make sure you understand the different steps in a NicheNet analysis that
are described in that vignette before proceeding with this vignette and
performing a real NicheNet analysis on your data. In another vignette
[Perform NicheNet analysis starting from a Seurat object: step-by-step
analysis](seurat_steps.md):`vignette("seurat_steps", package="nichenetr")`,
we also show the execution of these steps one for one, but in contrast
to the main vignette now specifically for a Seurat Object. This allows
users to adapt specific steps of the pipeline to make them more
appropriate for their data (recommended).

As example expression data of interacting cells, we will use mouse
NICHE-seq data from Medaglia et al. to explore intercellular
communication in the T cell area in the inguinal lymph node before and
72 hours after lymphocytic choriomeningitis virus (LCMV) infection (See
Medaglia et al. 2017). We will NicheNet to explore immune cell crosstalk
in response to this LCMV infection.

In this dataset, differential expression is observed between CD8 T cells
in steady-state and CD8 T cells after LCMV infection. NicheNet can be
applied to look at how several immune cell populations in the lymph node
(i.e., monocytes, dendritic cells, NK cells, B cells, CD4 T cells) can
regulate and induce these observed gene expression changes. NicheNet
will specifically prioritize ligands from these immune cells and their
target genes that change in expression upon LCMV infection.

The used NicheNet networks, ligand-target matrix and example expression
data of interacting cells can be downloaded from Zenodo. The NicheNet
networks and ligand-target matrix at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)
and the Seurat object of the processed NICHE-seq single-cell data at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3531889.svg)](https://doi.org/10.5281/zenodo.3531889).

# Prepare NicheNet analysis

## Load required packages, read in the Seurat object with processed expression data of interacting cells and NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks.

The NicheNet ligand-receptor network and weighted networks are necessary
to define and show possible ligand-receptor interactions between two
cell populations. The ligand-target matrix denotes the prior potential
that particular ligands might regulate the expression of particular
target genes. This matrix is necessary to prioritize possible
ligand-receptor interactions based on observed gene expression effects
(i.e. NicheNet’s ligand activity analysis) and infer affected target
genes of these prioritized ligands.

### Load Packages:

``` r
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
```

If you would use and load other packages, we recommend to load these 3
packages after the others.

### Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
## # A tibble: 6 x 4
##   from  to    source         database
##   <chr> <chr> <chr>          <chr>   
## 1 CXCL1 CXCR2 kegg_cytokines kegg    
## 2 CXCL2 CXCR2 kegg_cytokines kegg    
## 3 CXCL3 CXCR2 kegg_cytokines kegg    
## 4 CXCL5 CXCR2 kegg_cytokines kegg    
## 5 PPBP  CXCR2 kegg_cytokines kegg    
## 6 CXCL6 CXCR2 kegg_cytokines kegg

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  ABCC6  0.422 
## 2 A1BG  ACE2   0.101 
## 3 A1BG  ADAM10 0.0970
## 4 A1BG  AGO1   0.0525
## 5 A1BG  AKT1   0.0855
## 6 A1BG  ANXA7  0.457
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  A2M    0.0294
## 2 AAAS  GFAP   0.0290
## 3 AADAC CYP3A4 0.0422
## 4 AADAC IRF8   0.0275
## 5 AATF  ATM    0.0330
## 6 AATF  ATR    0.0355
```

### Read in the expression data of interacting cells:

The dataset used here is publicly available single-cell data from immune
cells in the T cell area of the inguinal lymph node. The data was
processed and aggregated by applying the Seurat alignment pipeline. The
Seurat object contains this aggregated data. Note that this should be a
Seurat v3 object and that gene should be named by their official
mouse/human gene symbol.

``` r
seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj@meta.data %>% head()
##         nGene nUMI orig.ident aggregate res.0.6 celltype nCount_RNA nFeature_RNA
## W380370   880 1611      LN_SS        SS       1    CD8 T       1607          876
## W380372   541  891      LN_SS        SS       0    CD4 T        885          536
## W380374   742 1229      LN_SS        SS       0    CD4 T       1223          737
## W380378   847 1546      LN_SS        SS       1    CD8 T       1537          838
## W380379   839 1606      LN_SS        SS       0    CD4 T       1603          836
## W380381   517  844      LN_SS        SS       0    CD4 T        840          513
```

Visualize which cell populations are present: CD4 T cells (including
regulatory T cells), CD8 T cells, B cells, NK cells, dendritic cells
(DCs) and inflammatory monocytes

``` r
seuratObj@meta.data$celltype %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application
## .
##     B CD4 T CD8 T    DC  Mono    NK  Treg 
##   382  2562  1645    18    90   131   199
DimPlot(seuratObj, reduction = "tsne")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Visualize the data to see to which condition cells belong. The metadata
dataframe column that denotes the condition (steady-state or after LCMV
infection) is here called ‘aggregate.’

``` r
seuratObj@meta.data$aggregate %>% table()
## .
## LCMV   SS 
## 3886 1141
DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Perform the NicheNet analysis

In this case study, we want to apply NicheNet to predict which ligands
expressed by all immune cells in the T cell area of the lymph node are
most likely to have induced the differential expression in CD8 T cells
after LCMV infection.

As described in the main vignette, the pipeline of a basic NicheNet
analysis consist of the following steps:

-   1.  Define a “sender/niche” cell population and a “receiver/target”
        cell population present in your expression data and determine
        which genes are expressed in both populations

-   2.  Define a gene set of interest: these are the genes in the
        “receiver/target” cell population that are potentially affected
        by ligands expressed by interacting cells (e.g. genes
        differentially expressed upon cell-cell interaction)

-   3.  Define a set of potential ligands: these are ligands that are
        expressed by the “sender/niche” cell population and bind a
        (putative) receptor expressed by the “receiver/target”
        population

-   4.  Perform NicheNet ligand activity analysis: rank the potential
        ligands based on the presence of their target genes in the gene
        set of interest (compared to the background set of genes)

-   5.  Infer receptors and top-predicted target genes of ligands that
        are top-ranked in the ligand activity analysis

All these steps are contained in one of three following similar single
functions: `nichenet_seuratobj_aggregate`,
`nichenet_seuratobj_cluster_de` and
`nichenet_seuratobj_aggregate_cluster_de`.

In addition to these steps, the function `nichenet_seuratobj_aggregate`
that is used for the analysis when having two conditions will also
calculate differential expression of the ligands in the sender cell
type. Note that this ligand differential expression is not used for
prioritization and ranking of the ligands!

## NicheNet analysis on Seurat object: explain differential expression between two conditions

In this case study, the receiver cell population is the ‘CD8 T’ cell
population, whereas the sender cell populations are ‘CD4 T,’ ‘Treg,’
‘Mono,’ ‘NK,’ ‘B’ and ‘DC.’ The above described functions will consider
a gene to be expressed when it is expressed in at least a predefined
fraction of cells in one cluster (default: 10%).

The gene set of interest are the genes differentially expressed in CD8 T
cells after LCMV infection. The condition of interest is thus ‘LCMV,’
whereas the reference/steady-state condition is ‘SS.’ The notion of
conditions can be extracted from the metadata column ‘aggregate,’ the
method to calculate the differential expression is the standard Seurat
Wilcoxon test.

The number of top-ranked ligands that are further used to predict active
target genes and construct an active ligand-receptor network is 20 by
default.

To perform the NicheNet analysis with these specifications, run the
following:

``` r
# indicated cell types should be cell class identities
# check via: 
# seuratObj %>% Idents() %>% table()
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "CD8 T", 
  condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", 
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"
```

### Interpret the NicheNet analysis output

#### Ligand activity analysis results

A first thing NicheNet does, is prioritizing ligands based on predicted
ligand activity. To see the ranking of these ligands, run the following
command:

``` r
nichenet_output$ligand_activities
## # A tibble: 44 x 6
##    test_ligand auroc  aupr pearson  rank bona_fide_ligand
##    <chr>       <dbl> <dbl>   <dbl> <dbl> <lgl>           
##  1 Ebi3        0.638 0.234  0.197      1 FALSE           
##  2 Il15        0.582 0.163  0.0961     2 TRUE            
##  3 Crlf2       0.549 0.163  0.0758     3 FALSE           
##  4 App         0.499 0.141  0.0655     4 TRUE            
##  5 Tgfb1       0.494 0.140  0.0558     5 TRUE            
##  6 Ptprc       0.536 0.149  0.0554     6 TRUE            
##  7 H2-M3       0.525 0.157  0.0528     7 TRUE            
##  8 Icam1       0.543 0.142  0.0486     8 TRUE            
##  9 Cxcl10      0.531 0.141  0.0408     9 TRUE            
## 10 Adam17      0.517 0.137  0.0359    10 TRUE            
## # ... with 34 more rows
```

The different ligand activity measures (auroc, aupr, pearson correlation
coefficient) are a measure for how well a ligand can predict the
observed differentially expressed genes compared to the background of
expressed genes. In our validation study, we showed that the pearson
correlation coefficient between a ligand’s target predictions and the
observed transcriptional response was the most informative measure to
define ligand activity. Therefore, NicheNet ranks the ligands based on
their pearson correlation coefficient. This allows us to prioritize
ligands inducing the antiviral response in CD8 T cells.

The column ‘bona\_fide\_ligand’ indicates whether the ligand is part of
ligand-receptor interactions that are documented in public databases
(‘bona\_fide\_ligand = TRUE’) and not of ligand-receptor interactions
that we predicted based on annotation as ligand/receptor and
protein-protein interaction databases (‘bona\_fide\_ligand = FALSE’).

To get a list of the 20 top-ranked ligands: run the following command

``` r
nichenet_output$top_ligands
##  [1] "Ebi3"   "Il15"   "Crlf2"  "App"    "Tgfb1"  "Ptprc"  "H2-M3"  "Icam1"  "Cxcl10" "Adam17" "Cxcl11" "Cxcl9"  "H2-T23" "Sema4d" "Ccl5"   "C3"     "Cxcl16" "Itgb1"  "Anxa1"  "Sell"
```

These ligands are expressed by one or more of the input sender cells. To
see which cell population expresses which of these top-ranked ligands,
you can run the following:

``` r
nichenet_output$ligand_expression_dotplot
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

As you can see, most op the top-ranked ligands seem to be mainly
expressed by dendritic cells and monocytes.

It could also be interesting to see whether some of these ligands are
differentially expressed after LCMV infection.

``` r
nichenet_output$ligand_differential_expression_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

As you can see, most op the top-ranked ligands seem also to be
upregulated themselves in monocytes after viral infection. This is not a
prerequisite to be top-ranked (cf: ranking only determined based on
enrichment of target genes among DE genes in the receiver, CD8T cells),
but is nice additional “evidence” that these ligands might indeed be
important.

#### Inferred active ligand-target links

NicheNet also infers active target genes of these top-ranked ligands. To
see which top-ranked ligands are predicted to have regulated the
expression of which differentially expressed genes, you can run
following command for a heatmap visualization:

``` r
nichenet_output$ligand_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

This is a normal ggplot object that can be adapted likewise. For example
if you want to change the color code to blue instead of purple, change
the axis ticks of the legend, and change the axis labels of the heatmap,
you can do the following:

``` r
nichenet_output$ligand_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + xlab("anti-LCMV response genes in CD8 T cells") + ylab("Prioritized immmune cell ligands")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

If you want, you can also extract the ligand-target links and their
regulatory potential scores in matrix or data frame format (e.g. for
visualization in other ways or output to a csv file).

``` r
nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
##              Cd274 Cd53       Ddit4         Id3 Ifit3        Irf1
## Sell   0.000000000    0 0.001290863 0.001222706     0 0.001095100
## Itgb1  0.000000000    0 0.001162142 0.001214922     0 0.001069406
## C3     0.000000000    0 0.001105490 0.000000000     0 0.000000000
## Ccl5   0.000000000    0 0.001281096 0.001228147     0 0.001155790
## Sema4d 0.000000000    0 0.001103465 0.001179496     0 0.000000000
## H2.T23 0.000000000    0 0.001112018 0.001110184     0 0.000000000
## Adam17 0.002280965    0 0.001760241 0.001546186     0 0.001637201
## Cxcl10 0.000000000    0 0.001354334 0.001372142     0 0.001393116
## Icam1  0.000000000    0 0.001325195 0.001314746     0 0.001375860
## H2.M3  0.000000000    0 0.001436893 0.001506164     0 0.001329158
```

``` r
nichenet_output$ligand_target_df # weight column = regulatory potential
## # A tibble: 155 x 3
##    ligand target  weight
##    <chr>  <chr>    <dbl>
##  1 Ebi3   Cd274  0.00325
##  2 Ebi3   Cd53   0.00321
##  3 Ebi3   Ddit4  0.00335
##  4 Ebi3   Id3    0.00373
##  5 Ebi3   Ifit3  0.00320
##  6 Ebi3   Irf1   0.00692
##  7 Ebi3   Irf7   0.00312
##  8 Ebi3   Irf9   0.00543
##  9 Ebi3   Parp14 0.00336
## 10 Ebi3   Pdcd4  0.00335
## # ... with 145 more rows
```

To get a list of the top-predicted target genes of the 20 top-ranked
ligands: run the following command

``` r
nichenet_output$top_targets
##  [1] "Cd274"  "Cd53"   "Ddit4"  "Id3"    "Ifit3"  "Irf1"   "Irf7"   "Irf9"   "Parp14" "Pdcd4"  "Pml"    "Psmb9"  "Rnf213" "Stat1"  "Stat2"  "Tap1"   "Ubc"    "Zbp1"   "Cd69"   "Gbp4"   "Basp1"  "Casp8" 
## [23] "Cxcl10" "Nlrc5"  "Vim"    "Actb"   "Ifih1"  "Myh9"   "B2m"    "H2-T23" "Rpl13a" "Cxcr4"
```

You can visualize the expression of these as well. Because we only focus
on CD8 T cells as receiver cells, we will only show expression in these
cells. To emphasize that these target genes are differentially
expressed, we split cells up in steadys-state cells and cells after
response to LCMV infection.

``` r
DotPlot(seuratObj %>% subset(idents = "CD8 T"), features = nichenet_output$top_targets %>% rev(), split.by = "aggregate") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
VlnPlot(seuratObj %>% subset(idents = "CD8 T"), features = c("Zbp1","Ifit3","Irf7"), split.by = "aggregate",    pt.size = 0, combine = FALSE)
## [[1]]
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

    ## 
    ## [[2]]

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

    ## 
    ## [[3]]

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->

To visualize ligand activities, expression, differential expression and
target genes of ligands, run the following command

``` r
nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

**important: above figure can be considered as one of the most important
summary figures of the NicheNet analysis. Here you can see which
ligand-receptor pairs have both high differential expression and ligand
activity (=target gene enrichment). These are very interesting
predictions as key regulators of your intercellular communication
process of interest ! **

#### Inferred ligand-receptor interactions for top-ranked ligands

NicheNet also infers the receiver cell receptors of these top-ranked
ligands. You can run following command for a heatmap visualization of
the ligand-receptor links:

``` r
nichenet_output$ligand_receptor_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

If you want, you can also extract the ligand-receptor links and their
interaction confidence scores in matrix or data frame format (e.g. for
visualization in other ways or output to a csv file).

``` r
nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]
##           Cxcl9    Cxcl16    Cxcl11      Ccl5    Cxcl10       App
## Cxcr6 0.3629049 0.6598705 0.2255185 0.2627207 0.4001071 0.2255185
## Ccr7  0.2217117 0.2217117 0.3567789 0.3933531 0.2582858 0.2217117
## Ccr9  0.1357118 0.2374693 0.2374693 0.1357118 0.1357118 0.1357118
## Gpr18 0.1374828 0.1374828 0.1374828 0.1374828 0.1374828 0.0000000
## S1pr1 0.1263826 0.1263826 0.1263826 0.1263826 0.1263826 0.0000000
## Itga4 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Cd47  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Ptk2b 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Il2rb 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Il2rg 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
```

``` r
nichenet_output$ligand_receptor_df # weight column accords to number of data sources that document this interaction
## # A tibble: 61 x 3
##    ligand receptor weight
##    <chr>  <chr>     <dbl>
##  1 Adam17 Notch1    0.482
##  2 Anxa1  Ccr7      0.222
##  3 Anxa1  Ccr9      0.237
##  4 Anxa1  Cxcr6     0.226
##  5 Anxa1  Itga4     0.201
##  6 App    Ccr7      0.222
##  7 App    Ccr9      0.136
##  8 App    Cxcr6     0.226
##  9 App    Notch1    0.354
## 10 App    Tgfbr2    0.441
## # ... with 51 more rows
```

To get a list of the receptors of the 20 top-ranked ligands: run the
following command

``` r
nichenet_output$top_receptors
##  [1] "Notch1" "Ccr7"   "Ccr9"   "Cxcr6"  "Itga4"  "Tgfbr2" "Itgb2"  "Gpr18"  "S1pr1"  "Il7r"   "Il27ra" "Cd8a"   "Klrd1"  "Il2rg"  "Itgal"  "Spn"    "Il2rb"  "Cd47"   "Ptk2b"  "Cd2"    "Cd28"   "Selplg"
## [23] "Ptprc"
```

You can visualize the expression of these as well. Because we only focus
on CD8 T cells as receiver cells, we will only show expression in these
cells.

``` r
DotPlot(seuratObj %>% subset(idents = "CD8 T"), features = nichenet_output$top_receptors %>% rev(), split.by = "aggregate") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

You also can just show ‘bona fide’ ligand-receptor links that are
described in the literature and not predicted based on protein-protein
interactions:

``` r
nichenet_output$ligand_receptor_heatmap_bonafide
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
nichenet_output$ligand_receptor_matrix_bonafide
##            H2.M3     H2.T23        C3      Icam1     Tgfb1    Cxcl16      Il15
## Il2rb  0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.0000000 0.8714269
## Il2rg  0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.0000000 0.8587859
## Itgal  0.0000000 0.00000000 0.0000000 0.06542904 0.0000000 0.0000000 0.0000000
## Itgb2  0.0000000 0.00000000 0.2916032 0.06113009 0.0000000 0.0000000 0.0000000
## Tgfbr2 0.0000000 0.00000000 0.0000000 0.00000000 0.7665905 0.0000000 0.0000000
## Cxcr6  0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.6598705 0.0000000
## Klrd1  0.8334165 0.05478448 0.0000000 0.00000000 0.0000000 0.0000000 0.0000000
nichenet_output$ligand_receptor_df_bonafide
## # A tibble: 9 x 3
##   ligand receptor weight
##   <chr>  <chr>     <dbl>
## 1 C3     Itgb2    0.292 
## 2 Cxcl16 Cxcr6    0.660 
## 3 H2-T23 Klrd1    0.0548
## 4 H2-M3  Klrd1    0.833 
## 5 Icam1  Itgal    0.0654
## 6 Icam1  Itgb2    0.0611
## 7 Il15   Il2rb    0.871 
## 8 Il15   Il2rg    0.859 
## 9 Tgfb1  Tgfbr2   0.767
```

If you are interested in checking which geneset (and background set of
genes) was used during the ligand activity analysis:

``` r
nichenet_output$geneset_oi
##   [1] "Irf7"          "Stat1"         "Ifit3"         "Ifit1"         "Bst2"          "B2m"           "Rnf213"        "Plac8"         "Isg15"         "Shisa5"        "Zbp1"          "Isg20"        
##  [13] "Samhd1"        "Usp18"         "H2-T23"        "Gbp2"          "Ifi203"        "Tmsb4x"        "Rsad2"         "Ly6e"          "Rtp4"          "Ifit2"         "Xaf1"          "Smchd1"       
##  [25] "Daxx"          "Alb"           "Samd9l"        "Actb"          "Parp9"         "Gbp4"          "Lgals3bp"      "Mx1"           "Gbp7"          "Cmpk2"         "Dtx3l"         "Slfn5"        
##  [37] "Oasl1"         "Herc6"         "Ifih1"         "Rpsa"          "P2ry13"        "Irgm2"         "Tapbp"         "Rps8"          "Stat2"         "Ifi44"         "Rpl8"          "Psmb8"        
##  [49] "Igfbp4"        "Ddx58"         "Rac2"          "Trafd1"        "Pml"           "Oas2"          "Psme1"         "Apoe"          "Basp1"         "Rps27a"        "Znfx1"         "Rpl13"        
##  [61] "Oas3"          "Nt5c3"         "Rnf114"        "Tap1"          "Rps28"         "Rplp0"         "Ddx60"         "Vim"           "Ifi35"         "Itm2b"         "Ctss"          "Pabpc1"       
##  [73] "Parp14"        "Hspa8"         "Tor3a"         "Rpl23"         "Tmbim6"        "Thy1"          "Ncoa7"         "Dhx58"         "Rps10"         "Rps19"         "Psmb9"         "Il2rg"        
##  [85] "Etnk1"         "Irf9"          "1600014C10Rik" "Parp12"        "Eif2ak2"       "Eef1b2"        "Eef2"          "Npc2"          "Rps2"          "Rps3"          "Sp110"         "Ube2l6"       
##  [97] "Nmi"           "Uba7"          "Psmb10"        "Cxcl10"        "Rpl13a"        "Nhp2"          "Tbrg1"         "Usp25"         "Tor1aip2"      "Adar"          "Gzma"          "Cd53"         
## [109] "Hspa5"         "Cfl1"          "Crip1"         "Slco3a1"       "Tlr7"          "Trim21"        "Rpl10"         "Mycbp2"        "Rps16"         "Nlrc5"         "Rplp2"         "Acadl"        
## [121] "Trim12c"       "Rps4x"         "Irf1"          "Psma2"         "Nme2"          "Zcchc11"       "Snord12"       "Phip"          "Ifitm3"        "Sp140"         "Dusp2"         "Mrpl30"       
## [133] "H2-M3"         "Gbp3"          "Dtx1"          "Eef1g"         "Rbl1"          "Xpo1"          "Gm9844"        "Rpl35"         "Rps26"         "Cxcr4"         "Eif3m"         "Treml2"       
## [145] "Rpl35a"        "Pdcd4"         "Arrb2"         "Ubc"           "Clic4"         "Rpl10a"        "Lcp1"          "Cd274"         "Ddit4"         "Cnn2"          "Nampt"         "Ascc3"        
## [157] "Cd47"          "Snord49b"      "D17Wsu92e"     "Fam26f"        "Hcst"          "Myh9"          "Rps27"         "Mov10"         "Arf4"          "Arhgdib"       "Ppib"          "Trim25"       
## [169] "Tspo"          "Id3"           "Snord35a"      "Rnf8"          "Casp8"         "Ptpn7"         "Itk"           "Cd69"          "Nop10"         "Anxa6"         "Hk1"           "Prkcb"        
## [181] "Iqgap1"        "Keap1"         "Rpl7"          "Parp10"
nichenet_output$background_expressed_genes %>% length()
## [1] 1487
```

### Rerun the NicheNet analysis with different sender cell definition

Instead of focusing on multiple sender cell types, it is possible that
you are only interested in doing the analyis for one sender cell type,
such as dendritic cells in this case.

``` r
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seuratObj, receiver = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "DC", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

Instead of focusing on one or multiple predefined sender cell types, it
is also possible that you want to consider all cell types present as
possible sender cell. This also includes the receiver cell type, making
that you can look at autocrine signaling as well.

``` r
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seuratObj, receiver = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "all", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

In some cases, it could be possible that you don’t have data of
potential sender cells. If you still want to predict possible upstream
ligands that could have been responsible for the observed differential
expression in your cell type, you can do this by following command. This
will consider all possible ligands in the NicheNet databases for which a
receptor is expressed by the receiver cell of interest.

``` r
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seuratObj, receiver = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "undefined", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"

nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

As you can see in this analysis result, many genes DE in CD8 T cells
after LCMV infection are strongly predicted type I interferon targets.
The presence of a type I interferon signature in the receiver cell type,
but the absence of expression of type I interferons in sender cell
types, might indicate that type I interferons are expressed by a
different, non-profiled cell type, or at a time point before sampling.
The latter could make sense, because there always is a time delay
between expression of a ligand-encoding gene and the effect of the
ligand on a target/receiver cell (i.e. expression of target genes).

### Run multiple NicheNet analyses on different receiver cell populations

In some cases, you might be interested in multiple target/receiver cell
populations. You can decide to run this for every cell type separately,
or in one line of code as demonstrated here (results are the same). As
example, we could have been interested in explaining DE between
steady-state and LCMV infection in both CD8 and CD4 T cells.

``` r
receiver_celltypes_oi = c("CD4 T", "CD8 T")
# receiver_celltypes_oi = seuratObj %>% Idents() %>% unique() # for all celltypes in the dataset: use only when this would make sense biologically

nichenet_output = receiver_celltypes_oi %>% lapply(nichenet_seuratobj_aggregate, seurat_obj = seuratObj, condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

names(nichenet_output) = receiver_celltypes_oi
```

Check which ligands were top-ranked for both CD8T and CD4T and which
ligands were more cell-type specific

``` r
common_ligands = intersect(nichenet_output$`CD4 T`$top_ligands, nichenet_output$`CD8 T`$top_ligands)
print("common ligands are: ")
## [1] "common ligands are: "
print(common_ligands)
##  [1] "Ebi3"   "Il15"   "Crlf2"  "H2-M3"  "App"    "Ptprc"  "Icam1"  "Ccl5"   "Cxcl10" "Tgfb1"  "Cxcl11" "Sema4d" "Cxcl9"  "H2-T23" "Cxcl16" "C3"     "Itgb1"

cd4_ligands = nichenet_output$`CD4 T`$top_ligands %>% setdiff(nichenet_output$`CD8 T`$top_ligands)
cd8_ligands = nichenet_output$`CD8 T`$top_ligands %>% setdiff(nichenet_output$`CD4 T`$top_ligands)

print("Ligands specifically regulating DE in CD4T: ")
## [1] "Ligands specifically regulating DE in CD4T: "
print(cd4_ligands)
## [1] "Cd274" "Hmgb1" "Cd28"

print("Ligands specifically regulating DE in CD8T: ")
## [1] "Ligands specifically regulating DE in CD8T: "
print(cd8_ligands)
## [1] "Adam17" "Anxa1"  "Sell"
```

## NicheNet analysis on Seurat object: explain differential expression between two cell populations

Previously, we demonstrated the use of a wrapper function for applying
NicheNet to explain differential expression between two conditions in
one cell type. However, also differential expression between two cell
populations might sometimes be (partially) caused by communication with
cells in the neighborhood. For example, differentiation from a
progenitor cell to the differentiated cell might be induced by niche
cells. A concrete example is discussed in this paper: [Stellate Cells,
Hepatocytes, and Endothelial Cells Imprint the Kupffer Cell Identity on
Monocytes Colonizing the Liver Macrophage
Niche](https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1).

Therefore, we will now also demonstrate the use of another Seurat
wrapper function that can be used in the case of explaining differential
expression between cell populations. But keep in mind that the
comparison that you make should be biologically relevant. It is possible
to use NicheNet to explain differential expression beween any two cell
populations in your dataset, but in most cases, differential expression
between cell populations will be a result of cell-intrinisc properties
(i.e. different cell types have a different gene expression profile) and
not of intercellular communication processes. In such a case, it does
not make any sense to use NicheNet.

For demonstration purposes, we will here first change the seuratObject
of the data described above, such that it can be used in this setting.

``` r
seuratObj@meta.data$celltype = paste(seuratObj@meta.data$celltype,seuratObj@meta.data$aggregate, sep = "_")

seuratObj@meta.data$celltype %>% table()
## .
##     B_LCMV       B_SS CD4 T_LCMV   CD4 T_SS CD8 T_LCMV   CD8 T_SS    DC_LCMV      DC_SS  Mono_LCMV    Mono_SS    NK_LCMV      NK_SS  Treg_LCMV    Treg_SS 
##        344         38       1961        601       1252        393         14          4         75         15         94         37        146         53

seuratObj = SetIdent(seuratObj,value = "celltype")
```

Now perform the NicheNet analysis to explain differential expression
between the ‘affected’ cell population ‘CD8 T cells after LCMV
infection’ and the reference cell population ‘CD8 T cells in
steady-state’ by ligands expressed by monocytes and DCs after LCMV
infection.

``` r
nichenet_output = nichenet_seuratobj_cluster_de(
  seurat_obj = seuratObj, 
  receiver_reference = "CD8 T_SS", receiver_affected = "CD8 T_LCMV", 
  sender = c("DC_LCMV","Mono_LCMV"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis between two receiver cell clusters"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
```

Check the top-ranked ligands and their target genes

``` r
nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Check the expression of the top-ranked ligands

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

It could be interested to check which top-ranked ligands are
differentially expressed in monocytes after LCMV infection

``` r
Mono_upregulated_ligands = FindMarkers(seuratObj, ident.1 = "Mono_LCMV", ident.2 = "Mono_SS") %>% rownames_to_column("gene") %>% filter(avg_log2FC > 0.25 & p_val_adj <= 0.05) %>% pull(gene) %>% intersect(nichenet_output$top_ligands)

print("Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-StSt and CD8T-LCMV are: ")
## [1] "Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-StSt and CD8T-LCMV are: "
print(Mono_upregulated_ligands)
## [1] "Cxcl10"
```

# Remarks

1.  Top-ranked ligands and target genes shown here differ from the
    predictions shown in the respective case study in the NicheNet paper
    because a different definition of expressed genes was used.
2.  Differential expression is here done via the classical Wilcoxon test
    used in Seurat to define marker genes of a cell cluster by comparing
    it to other clusters. This is not optimal if you would have repeated
    samples for your conditions. In such a case, we recommend to follow
    the vignette [Perform NicheNet analysis starting from a Seurat
    object: step-by-step
    analysis](seurat_steps.md):`vignette("seurat_steps", package="nichenetr")`
    and tweak the differential expression step there (and perform the
    analysis e.g. as discussed in <https://github.com/HelenaLC/muscat>).

# References

Bonnardel et al., 2019, Immunity 51, 1–17, [Stellate Cells, Hepatocytes,
and Endothelial Cells Imprint the Kupffer Cell Identity on Monocytes
Colonizing the Liver Macrophage
Niche](https://doi.org/10.1016/j.immuni.2019.08.017)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-medaglia_spatial_2017" class="csl-entry">

Medaglia, Chiara, Amir Giladi, Liat Stoler-Barak, Marco De Giovanni,
Tomer Meir Salame, Adi Biram, Eyal David, et al. 2017. “Spatial
Reconstruction of Immune Niches by Combining Photoactivatable Reporters
and <span class="nocase">scRNA</span>-Seq.” *Science*, December,
eaao4277. <https://doi.org/10.1126/science.aao4277>.

</div>

</div>
Assess how well top-ranked ligands can predict a gene set of interest
================
Robin Browaeys
2019-02-19

<!-- github markdown built using
rmarkdown::render("vignettes/target_prediction_evaluation_geneset.Rmd", output_format = "github_document")
-->

This vignette shows how NicheNet can be used to to predict which ligands
might regulate a given set of genes and how well they do this
prediction. For this analysis, you need to define:

  - a set of genes of which expression in a “receiver cell” is possibly
    affected by extracellular protein signals (ligands) (e.g. genes
    differentially expressed upon cell-cell interaction )
  - a set of potentially active ligands (e.g. ligands expressed by
    interacting “sender cells”)

Therefore, you often first need to process expression data of
interacting cells to define both.

In this example, we will use data from Puram et al. to explore
intercellular communication in the tumor microenvironment in head and
neck squamous cell carcinoma (HNSCC) (See Puram et al. 2017). More
specifically, we will look at which ligands expressed by
cancer-associated fibroblasts (CAFs) can induce a specific gene program
in neighboring malignant cells. This program, a partial
epithelial-mesenschymal transition (p-EMT) program, could be linked by
Puram et al. to metastasis.

For this analysis, we will first assess the ligand activity of each
ligand, or in other words, we will assess how well each CAF-ligand can
predict the p-EMT gene set compared to the background of expressed
genes. This allows us to prioritize p-EMT-regulating ligands. Then, we
will assess how well the prioritized ligands together can predict
whether genes belong to the gene set of interest or not.

The used ligand-target matrix and example expression data of interacting
cells can be downloaded from Zenodo.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)

### Load packages required for this vignette

``` r
library(nichenetr)
library(tidyverse)
```

### Read in expression data of interacting cells

First, we will read in the publicly available single-cell data from CAF
and malignant cells from HNSCC
tumors.

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

Secondly, we will determine which genes are expressed in CAFs and
malignant cells from high quality primary tumors. Therefore, we wil not
consider cells from tumor samples of less quality or from lymph node
metastases. To determine expressed genes, we use the definition used by
of Puram et
al.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% pull(cell)
malignant_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
```

### Load the ligand-target model we want to use

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

### Load the gene set of interest and background of genes

As gene set of interest, we consider the genes of which the expression
is possibly affected due to communication with other cells.

Because we here want to investigate how CAF regulate the expression of
p-EMT genes in malignant cells, we will use the p-EMT gene set defined
by Puram et al. as gene set of interset and use all genes expressed in
malignant cells as background of
genes.

``` r
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(pemt_geneset)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"
```

### Perform NicheNet’s ligand activity analysis on the gene set of interest

In a first step, we will define a set of potentially active ligands. As
potentially active ligands, we will use ligands that are 1) expressed by
CAFs and 2) can bind a (putative) receptor expressed by malignant cells.
Putative ligand-receptor links were gathered from NicheNet’s
ligand-receptor data
sources.

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_CAFs)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"
```

Now perform the ligand activity analysis: infer how well NicheNet’s
ligand-target potential scores can predict whether a gene belongs to the
p-EMT program or
not.

``` r
ligand_activities = predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

Now, we want to rank the ligands based on their ligand activity. In our
validation study, we showed that the pearson correlation between a
ligand’s target predictions and the observed transcriptional response
was the most informative measure to define ligand activity. Therefore,
we will rank the ligands based on their pearson correlation coefficient.

``` r
ligand_activities %>% arrange(-pearson)
## # A tibble: 131 x 4
##    test_ligand auroc   aupr pearson
##    <chr>       <dbl>  <dbl>   <dbl>
##  1 PTHLH       0.667 0.0720   0.128
##  2 CXCL12      0.680 0.0507   0.123
##  3 AGT         0.676 0.0581   0.120
##  4 TGFB3       0.689 0.0454   0.117
##  5 IL6         0.693 0.0510   0.115
##  6 INHBA       0.695 0.0502   0.113
##  7 ADAM17      0.672 0.0526   0.113
##  8 TNC         0.700 0.0444   0.109
##  9 CTGF        0.680 0.0473   0.108
## 10 FN1         0.679 0.0505   0.108
## # ... with 121 more rows
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)
## [1] "PTHLH"  "CXCL12" "AGT"    "TGFB3"  "IL6"    "INHBA"
```

For the top 20 ligands, we will now build a multi-ligand model that uses
all top-ranked ligands to predict whether a gene belongs to the p-EMT
program of not. This classification model will be trained via
cross-validation and returns a probability for every
gene.

``` r
# change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
k = 3 # 3-fold
n = 2 # 2 rounds

pemt_gene_predictions_top20_list = seq(n) %>% lapply(assess_rf_class_probabilities, folds = k, geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligands_oi = best_upstream_ligands, ligand_target_matrix = ligand_target_matrix)
```

Evaluate now how well the target gene probabilies accord to the gene set
assignments

``` r
# get performance: auroc-aupr-pearson
target_prediction_performances_cv = pemt_gene_predictions_top20_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>% bind_rows() %>% mutate(round=seq(1:nrow(.)))
```

What is the AUROC, AUPR and PCC of this model (averaged over
cross-validation rounds)?

``` r
target_prediction_performances_cv$auroc %>% mean()
## [1] 0.7295863
target_prediction_performances_cv$aupr %>% mean()
## [1] 0.07603073
target_prediction_performances_cv$pearson %>% mean()
## [1] 0.1660327
```

Evaluate now whether genes belonging to the gene set are more likely to
be top-predicted. We will look at the top 5% of predicted targets
here.

``` r
# get performance: how many p-EMT genes and non-p-EMT-genes among top 5% predicted targets
target_prediction_performances_discrete_cv = pemt_gene_predictions_top20_list %>% lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>% bind_rows() %>% ungroup() %>% mutate(round=rep(1:length(pemt_gene_predictions_top20_list), each = 2))
```

What is the fraction of p-EMT genes that belongs to the top 5% predicted
targets?

``` r
target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean()
## [1] 0.25
```

What is the fraction of non-p-EMT genes that belongs to the top 5%
predicted
targets?

``` r
target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean()
## [1] 0.04769076
```

We see that the p-EMT genes are enriched in the top-predicted target
genes. To test this, we will now apply a Fisher’s exact test for every
cross-validation round and report the average
p-value.

``` r
target_prediction_performances_discrete_fisher = pemt_gene_predictions_top20_list %>% lapply(calculate_fraction_top_predicted_fisher, quantile_cutoff = 0.95) 
target_prediction_performances_discrete_fisher %>% unlist() %>% mean()
## [1] 5.647773e-10
```

Finally, we will look at which p-EMT genes are well-predicted in every
cross-validation round.

``` r
# get top predicted genes
top_predicted_genes = seq(length(pemt_gene_predictions_top20_list)) %>% lapply(get_top_predicted_genes,pemt_gene_predictions_top20_list) %>% reduce(full_join, by = c("gene","true_target"))
top_predicted_genes %>% filter(true_target)
## # A tibble: 28 x 4
##    gene    true_target predicted_top_target_roun~ predicted_top_target_rou~
##    <chr>   <lgl>       <lgl>                      <lgl>                    
##  1 COL1A1  TRUE        TRUE                       TRUE                     
##  2 MMP2    TRUE        TRUE                       TRUE                     
##  3 MMP1    TRUE        TRUE                       TRUE                     
##  4 PLAU    TRUE        TRUE                       TRUE                     
##  5 TIMP3   TRUE        TRUE                       TRUE                     
##  6 MT2A    TRUE        TRUE                       TRUE                     
##  7 INHBA   TRUE        TRUE                       TRUE                     
##  8 COL4A2  TRUE        TRUE                       TRUE                     
##  9 MMP10   TRUE        TRUE                       TRUE                     
## 10 COL17A1 TRUE        TRUE                       TRUE                     
## # ... with 18 more rows
```

### References

<div id="refs" class="references">

<div id="ref-puram_single-cell_2017">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
Evaluation of NicheNet’s ligand-target predictions
================
Robin Browaeys
2018-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/model_evaluation.Rmd", output_format = "github_document")
-->

This vignette shows how the ligand-target predictions of NicheNet were
evaluated. For validation, we collected transcriptome data of cells
before and after they were treated by one or two ligands in culture.
Using these ligand treatment datasets for validation has the advantage
that observed gene expression changes can be directly attributed to the
addition of the ligand(s). Hence, differentially expressed genes can be
considered as a gold standard of target genes of a particular ligand.

You can use the procedure shown here to evaluate your own model and
compare its performance to NicheNet. Ligand treatment validation
datasets and NicheNet’s ligand-target model can be downloaded from
Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758).

### Load nichenetr, the model we want to evaluate, and the datasets on which we want to evaluate it.

``` r
library(nichenetr)
library(tidyverse)

# Load in the ligand-target model
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

# The ligand treatment expression datasets used for validation can be downloaded from Zenodo:
expression_settings_validation = readRDS(url("https://zenodo.org/record/3260758/files/expression_settings.rds"))

#Ligand treatment datasets show the log fold change in expression of genes after treatment with one or more specific ligands. Here: example for the ligand NODAL:
head(expression_settings_validation$nodal_Nodal$diffexp)
##            lfc       qval   gene
## 1 -0.072591381 0.03767323  BEST1
## 2 -0.133032792 0.62838239 SMIM37
## 3 -0.001152079 0.98974906  CRADD
## 4 -0.018723210 0.63536958  RCAN2
## 5 -0.070250368 0.36418916 MFAP3L
## 6 -0.192244949 0.23454428   PARL
```

### Example: transcriptional response prediction evaluation

First, we will demonstrate how to evaluate the transcriptional response
(i.e. target gene prediction) performance for all ligand treatment
expression datasets. For this, we determine how well the model predicts
which genes are differentially expressed after treatment with a ligand.
Ideally, target genes with high regulatory potential scores for a
ligand, should be differentially expressed in response to that ligand.

For information of all collected ligand treatment datasets, see [Dataset
information](evaluation_datasets.xlsx)

For the sake of simplicity, we exclude in this vignette the
ligand-treatment datasets profiling the response to multiple ligands. To
see how to build a ligand-target model with target predictions for
multiple ligands at once: see vignette [Construction of NicheNet’s
ligand-target model](model_construction.md):
`vignette("model_construction", package="nichenetr")`.

Step 1: convert expression datasets to the required format to perform
target gene
prediction

``` r
settings = expression_settings_validation %>% lapply(convert_expression_settings_evaluation)
settings = settings %>% discard(~length(.$from) > 1)
```

Step 2: calculate the target gene prediction performances

``` r
# Evaluate transcriptional response prediction on every dataset
performances = settings %>% lapply(evaluate_target_prediction, ligand_target_matrix) %>% bind_rows()
```

Step 3: visualize the results: show here different classification
evaluation
metrics

``` r
# Visualize some classification evaluation metrics showing the target gene prediction performance
performances = performances %>% select(-aupr, -auc_iregulon,-pearson_log_pval,-spearman_log_pval ,-sensitivity_roc, -specificity_roc) %>% gather(key = scorename, value = scorevalue, auroc:spearman)
scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)", auc_iregulon_corrected = "AUC-iRegulon (corrected)",pearson = "Pearson correlation", spearman = "Spearman's rank correlation",mean_rank_GST_log_pval = "Mean-rank gene-set enrichment")
scorerandom = c(auroc=0.5, aupr_corrected=0, auc_iregulon_corrected = 0, pearson = 0, spearman = 0,mean_rank_GST_log_pval = 0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

performances %>%
  mutate(model = "NicheNet") %>%
  ggplot() +
  geom_violin(aes(model, scorevalue, group=model, fill = model)) +
  geom_boxplot(aes(model, scorevalue, group = model),width = 0.05) +
  scale_y_continuous("Score target prediction") +
  facet_wrap(~scorename, scales = "free", labeller=as_labeller(scorelabels)) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, color = "red") +
  theme_bw()
```

![](model_evaluation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Example: ligand activity prediction evaluation

Now we will show how to assess the accuracy of the model in predicting
whether cells were treated by a particular ligand or not. In other
words, we will evaluate how well NicheNet prioritizes active ligand(s),
given a set of differentially expressed genes. For this procedure, we
assume the following: the better a ligand predicts the transcriptional
response compared to other ligands, the more likely it is that this
ligand is active. Therefore, we first get ligand activity (or ligand
importance or feature importance) scores for all ligands on all
ligand-treatment expression datasets of which the true acive ligand is
known. Then we assess whether the truly active ligands get indeed higher
ligand activity scores as should be for a good ligand-target model.

A graphical summary of this procedure is visualized here below:

![](vignettes/ligand_activity_prediction_workflow_new.png)

Step 1: convert expression datasets to the required format to perform
ligand activity
prediction

``` r
# convert expression datasets to correct format for ligand activity prediction
all_ligands = settings %>% extract_ligands_from_settings(combination = FALSE) %>% unlist()
settings_ligand_prediction = settings %>% convert_settings_ligand_prediction(all_ligands = all_ligands, validation = TRUE)
```

Step 2: calculate the ligand importances (these are classification
evaluation metrics indicating how well a ligand can predict the observed
DE genes in a specific ligand treatment
dataset)

``` r
# infer ligand importances: for all ligands of interest, we assess how well a ligand explains the differential expression in a specific datasets (and we do this for all datasets).
ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances,ligand_target_matrix) %>% bind_rows()
```

Step 3: evaluate how separate ligand importances can predict ligand
activity

``` r
# Look at predictive performance of single/individual importance measures to predict ligand activity: of all ligands tested, the ligand that is truly active in a dataset should get the highest activity score (i.e. best target gene prediction performance)
evaluation_ligand_prediction = ligand_importances$setting %>% unique() %>% lapply(function(x){x}) %>%
    lapply(wrapper_evaluate_single_importances_ligand_prediction,ligand_importances) %>%
    bind_rows() %>% inner_join(ligand_importances %>% distinct(setting,ligand))
```

Step 4: visualize the results: show here different classification
evaluation
metrics

``` r
# Visualize some classification evaluation metrics showing the ligand activity prediction performance
evaluation_ligand_prediction = evaluation_ligand_prediction %>% select(-aupr, -sensitivity_roc, -specificity_roc, -pearson, -spearman, -mean_rank_GST_log_pval) %>% gather(key = scorename, value = scorevalue, auroc:aupr_corrected)
scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)")
scorerandom = c(auroc=0.5, aupr_corrected=0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

evaluation_ligand_prediction %>%
 filter(importance_measure %in% c("auroc", "aupr_corrected", "mean_rank_GST_log_pval", "auc_iregulon_corrected", "pearson", "spearman")) %>%
  ggplot() +
  geom_violin(aes(importance_measure, scorevalue, group=importance_measure, fill = importance_measure)) +
  geom_boxplot(aes(importance_measure, scorevalue, group = importance_measure),width = 0.1) +
  scale_y_continuous("Evaluation ligand activity prediction") +
  scale_x_discrete("Ligand activity measure") +
  facet_wrap(~scorename, scales = "free", labeller=as_labeller(scorelabels)) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
```

![](model_evaluation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

This plots shows that using the pearson correlation coefficient target
prediction metric is the best metric to use for ranking ligands
according to predicted ligand activity.
Differential NicheNet analysis between conditions of interest
================
Robin Browaeys
2022-01-12

<!-- github markdown built using 
rmarkdown::render("vignettes/differential_nichenet_pEMT.Rmd", output_format = "github_document")
-->

Remark: this is a beta version of a new extension of NicheNet, namely
Differential NicheNet. Short-term improvements will include scalability,
visualization and documentation of this vignette and the underlying
functions (january 2022).

The goal of Differential NicheNet is to predict ligand-receptors pairs
that are both differentially expressed and active between different
niches of interest.

This vignette guides you in detail through all the steps of a
Differential NicheNet analysis. As example expression data of
interacting cells, we will use data from Puram et al. to explore
intercellular communication in the tumor microenvironment in head and
neck squamous cell carcinoma (HNSCC) (Puram et al. 2017). More
specifically, we will look at cell-cell communication differences
between pEMT-high and pEMT-low tumors (pEMT = partial
epithelial-mesenschymal transition). In this data, we thus have 2
conditions/niches, but this pipeline is also usable for more
conditions/niches.

The used ligand-receptor network and ligand-target matrix can be
downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758).
The Seurat object containing expression data of interacting cells in
HNSCC can also be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4675430.svg)](https://doi.org/10.5281/zenodo.4675430).

# 0. Read in the expression data of interest, and the NicheNet ligand-receptor network and ligand-target matrix

## Load in packages

``` r
library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #
```

## Read in the expression data

In this case study, we want to study differences in cell-cell
communication patterns between pEMT-high and pEMT-low tumors. The meta
data columns that indicate the pEMT status of tumors are ‘pEMT,’ and the
cell type is indicated in the ‘celltype’ column.

``` r
seurat_obj = readRDS(url("https://zenodo.org/record/4675430/files/seurat_obj_hnscc.rds"))
DimPlot(seurat_obj, group.by = "celltype") # user adaptation required on own dataset
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
DimPlot(seurat_obj, group.by = "pEMT") # user adaptation required on own dataset
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

We will now also check the number of cells per cell type condition
combination

``` r
table(seurat_obj@meta.data$celltype, seurat_obj@meta.data$pEMT) # cell types vs conditions # user adaptation required on own dataset
##                
##                 High  Low
##   CAF            396  104
##   Endothelial    105   53
##   Malignant     1093  549
##   Myeloid         92    7
##   myofibroblast  382   61
##   T.cell         689    3
```

For the Differential NicheNet, we need to compare at least 2 niches or
conditions to each other. In this case, the 2 niches are the
pEMT-high-niche and the pEMT-low-niche. We will adapt the names of the
cell types based on their niche of origin.

``` r
seurat_obj@meta.data$celltype_aggregate = paste(seurat_obj@meta.data$celltype, seurat_obj@meta.data$pEMT,sep = "_") # user adaptation required on own dataset
DimPlot(seurat_obj, group.by = "celltype_aggregate")
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
seurat_obj@meta.data$celltype_aggregate %>% table() %>% sort(decreasing = TRUE)
## .
##     Malignant_High        T.cell_High      Malignant_Low           CAF_High myofibroblast_High   Endothelial_High 
##               1093                689                549                396                382                105 
##            CAF_Low       Myeloid_High  myofibroblast_Low    Endothelial_Low        Myeloid_Low         T.cell_Low 
##                104                 92                 61                 53                  7                  3
```

``` r
celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])
```

## Read in the NicheNet ligand-receptor network and ligand-target matrix

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

head(lr_network)
## # A tibble: 6 x 3
##   ligand receptor bonafide
##   <chr>  <chr>    <lgl>   
## 1 CXCL1  CXCR2    TRUE    
## 2 CXCL2  CXCR2    TRUE    
## 3 CXCL3  CXCR2    TRUE    
## 4 CXCL5  CXCR2    TRUE    
## 5 PPBP   CXCR2    TRUE    
## 6 CXCL6  CXCR2    TRUE
```

Note: if your data is of mouse origin: convert human gene symbols to
their one-to-one orthologs

``` r
organism = "human" # user adaptation required on own dataset
```

``` r
if(organism == "mouse"){
  lr_network = lr_network %>% mutate(ligand = convert_human_to_mouse_symbols(ligand), receptor = convert_human_to_mouse_symbols(receptor)) %>% drop_na()

  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
}
```

# 1. Define the niches/microenvironments of interest

Each niche should have at least one “sender/niche” cell population and
one “receiver/target” cell population (present in your expression data)

In this case study, we are interested to find differences in cell-cell
interactions to malignant cells between pEMT high and pEMT low tumors.
The receiver cell population in the pEMT-High niche is thus the
“Malignant\_High” cell type, and in the pEMT-Low niche this is
“Malignant\_Low.” The sender cell populations of interest are
myofibroblasts, Endothelial, CAF, T.cell, and Myeloid. Importantly, we
only include T.Cell and Myeloid in the pEMT-High niche, because there
are too few cells of these populations present in the pEMT-low niche.
Hereby, we demonstrate the possibility to include a condition-specific
cell type in the analysis - which is possible because we calculate DE
compared to all sender cells of the other niche, and not only to the
pEMT-low group of cells of the same cell type.

! Important: your receiver cell type should consist of 1 cluster!

``` r
niches = list(
  "pEMT_High_niche" = list(
    "sender" = c("myofibroblast_High", "Endothelial_High", "CAF_High", "T.cell_High", "Myeloid_High"),
    "receiver" = c("Malignant_High")),
  "pEMT_Low_niche" = list(
    "sender" = c("myofibroblast_Low",  "Endothelial_Low", "CAF_Low"),
    "receiver" = c("Malignant_Low"))
  ) # user adaptation required on own dataset
```

# 2. Calculate differential expression between the niches

In this step, we will determine DE between the different niches for both
senders and receivers to define the DE of L-R pairs.

### Calculate DE

The method to calculate the differential expression is here the standard
Seurat Wilcoxon test, but this can be replaced if wanted by the user
(only requirement: output tables `DE_sender_processed` and
`DE_receiver_processed` should be in the same format as shown here).

DE will be calculated for each pairwise sender (or receiver) cell type
comparision between the niches (so across niches, not within niche). In
our case study, this means that DE of myofibroblast\_High ligands will
be calculated by DE analysis of myofibroblast\_High vs
myofibroblast\_Low; myofibroblast\_High vs Endothelial\_Low; and
myofibroblast\_High vs CAF\_Low. We split the cells per cell type
instead of merging all cells from the other niche to avoid that the DE
analysis will be driven by the most abundant cell types.

``` r
assay_oi = "SCT" # other possibilities: RNA,...
DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
## [1] "Calculate Sender DE between: myofibroblast_High and myofibroblast_Low"
## [2] "Calculate Sender DE between: myofibroblast_High and Endothelial_Low"  
## [3] "Calculate Sender DE between: myofibroblast_High and CAF_Low"          
## [1] "Calculate Sender DE between: Endothelial_High and myofibroblast_Low"
## [2] "Calculate Sender DE between: Endothelial_High and Endothelial_Low"  
## [3] "Calculate Sender DE between: Endothelial_High and CAF_Low"          
## [1] "Calculate Sender DE between: CAF_High and myofibroblast_Low" "Calculate Sender DE between: CAF_High and Endothelial_Low"  
## [3] "Calculate Sender DE between: CAF_High and CAF_Low"          
## [1] "Calculate Sender DE between: T.cell_High and myofibroblast_Low"
## [2] "Calculate Sender DE between: T.cell_High and Endothelial_Low"  
## [3] "Calculate Sender DE between: T.cell_High and CAF_Low"          
## [1] "Calculate Sender DE between: Myeloid_High and myofibroblast_Low"
## [2] "Calculate Sender DE between: Myeloid_High and Endothelial_Low"  
## [3] "Calculate Sender DE between: Myeloid_High and CAF_Low"          
## [1] "Calculate Sender DE between: myofibroblast_Low and myofibroblast_High"
## [2] "Calculate Sender DE between: myofibroblast_Low and Endothelial_High"  
## [3] "Calculate Sender DE between: myofibroblast_Low and CAF_High"          
## [4] "Calculate Sender DE between: myofibroblast_Low and T.cell_High"       
## [5] "Calculate Sender DE between: myofibroblast_Low and Myeloid_High"      
## [1] "Calculate Sender DE between: Endothelial_Low and myofibroblast_High"
## [2] "Calculate Sender DE between: Endothelial_Low and Endothelial_High"  
## [3] "Calculate Sender DE between: Endothelial_Low and CAF_High"          
## [4] "Calculate Sender DE between: Endothelial_Low and T.cell_High"       
## [5] "Calculate Sender DE between: Endothelial_Low and Myeloid_High"      
## [1] "Calculate Sender DE between: CAF_Low and myofibroblast_High" "Calculate Sender DE between: CAF_Low and Endothelial_High"  
## [3] "Calculate Sender DE between: CAF_Low and CAF_High"           "Calculate Sender DE between: CAF_Low and T.cell_High"       
## [5] "Calculate Sender DE between: CAF_Low and Myeloid_High"
DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
## # A tibble: 1 x 2
##   receiver       receiver_other_niche
##   <chr>          <chr>               
## 1 Malignant_High Malignant_Low       
## [1] "Calculate receiver DE between: Malignant_High and Malignant_Low"
## [1] "Calculate receiver DE between: Malignant_Low and Malignant_High"
```

### Process DE results:

``` r
expression_pct = 0.10
DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
```

### Combine sender-receiver DE based on L-R pairs:

As mentioned above, DE of ligands from one sender cell type is
determined be calculating DE between that cell type, and all the sender
cell types of the other niche. To summarize the DE of ligands of that
cell type we have several options: we could take the average LFC, but
also the minimum LFC compared to the other niche. We recommend using the
minimum LFC, because this is the strongest specificity measure of ligand
expression, because a high min LFC means that a ligand is more strongly
expressed in the cell type of niche 1 compared to all cell types of
niche 2 (in contrast to a high average LFC, which does not exclude that
one or more cell types in niche 2 also strongly express that ligand).

``` r
specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
```

# 3. Optional: Calculate differential expression between the different spatial regions

To improve the cell-cell interaction predictions, you can consider
spatial information if possible and applicable. Spatial information can
come from microscopy data, or from spatial transcriptomics data such as
Visium.

There are several ways to incorporate spatial information in the
Differential NicheNet pipeline. First, you can only consider cell types
as belonging to the same niche if they are in the same spatial location.
Another way is including spatial differential expression of
ligand-receptor pairs within one cell type in the prioritization
framework.

For example: We have a cell type X, located in regions A and B, and we
want to study cell-cell communication in region A. We first add only
celltypeX of regionA in the niche definition, and then calculate DE
between celltypeX-regionA and celltypeX-regionB to give higher
prioritization weight to regionA-specific ligands.

In this case study, our region of interest is the tumor leading edge,
since Puram et al defined this region as important regarding the pEMT
process. Puram et al also defined CAFs as the fibroblasts that are close
to leading edge, whereas the other fibroblasts (myofibroblasts) were not
preferentially located in the tumor leading edge. We can thus now
prioritize fibroblast ligands further by looking at ligands that are DE
between leading-edge fibroblasts (=CAFs) and non-leading-edge
fibroblasts (myofibroblasts).

We do this as follows, by first defining a ‘spatial info’ dataframe. If
no spatial information in your data: set the following two parameters to
FALSE, and make a mock ‘spatial\_info’ data frame.

``` r
include_spatial_info_sender = TRUE # if not spatial info to include: put this to false # user adaptation required on own dataset
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
```

``` r
spatial_info = tibble(celltype_region_oi = "CAF_High", celltype_other_region = "myofibroblast_High", niche =  "pEMT_High_niche", celltype_type = "sender") # user adaptation required on own dataset
specificity_score_spatial = "lfc"
```

``` r
# this is how this should be defined if you don't have spatial info
# mock spatial info
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
} 
```

``` r
if(include_spatial_info_sender == TRUE){
  sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)

  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))

} else {
  # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
  sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  

}
## [1] "Calculate Spatial DE between: CAF_High and myofibroblast_High"
```

``` r
if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
  receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)

  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))

} else {
    # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}
```

# 4. Calculate ligand activities and infer active ligand-target links

In this step, we will predict ligand activities of each ligand for each
of the receiver cell types across the different niches. This is similar
to the ligand activity analysis done in the normal NicheNet pipeline.

To calculate ligand activities, we first need to define a geneset of
interest for each niche. In this case study, the geneset of interest for
the pEMT-high niche are the genes upregulated in pEMT-high tumors
compared to pEMT-low tumors, and vice versa.

Note that you can also define these geneset of interest in a different
way! (eg pathway-based geneset etc)

Ligand-target links are inferred in the same way as described in the
basic NicheNet vignettes.

``` r
lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"

DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
## [1] "Calculate receiver DE between: Malignant_High and Malignant_Low"
## [1] "Calculate receiver DE between: Malignant_Low and Malignant_High"
DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  
# Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
# If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
##  [1] "ANXA8L2"       "PRKCDBP"       "IL8"           "PTRF"          "SEPP1"         "C1orf186"      "CCDC109B"     
##  [8] "C10orf54"      "LEPREL1"       "ZNF812"        "LOC645638"     "LOC401397"     "LINC00162"     "DFNA5"        
## [15] "PLK1S1"        "ZMYM6NB"       "C19orf10"      "CTSL1"         "SQRDL"         "LOC375295"     "WBP5"         
## [22] "LOC100505633"  "AIM1"          "C1orf63"       "LOC100507463"  "GPR115"        "VIMP"          "SEP15"        
## [29] "C1orf172"      "NAPRT1"        "LHFP"          "KRT16P1"       "C7orf10"       "PTPLA"         "GRAMD3"       
## [36] "CPSF3L"        "MESDC2"        "C10orf10"      "KIAA1609"      "CCDC53"        "TXLNG2P"       "NGFRAP1"      
## [43] "ERO1L"         "FAM134A"       "LSMD1"         "TCEB2"         "B3GALTL"       "HN1L"          "LOC550643"    
## [50] "KIAA0922"      "GLT25D1"       "FAM127A"       "C1orf151-NBL1" "SEPW1"         "GPR126"        "LOC100505806" 
## [57] "LINC00478"     "TCEB1"         "GRAMD2"        "GNB2L1"        "KIRREL"
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
##   [1] "LOC344887"    "AGPAT9"       "C1orf110"     "KIAA1467"     "LOC100292680" "EPT1"         "CT45A4"       "LOC654433"   
##   [9] "UPK3BL"       "LINC00340"    "LOC100128338" "FAM60A"       "CCDC144C"     "LOC401109"    "LOC286467"    "LEPREL4"     
##  [17] "LOC731275"    "LOC642236"    "LINC00516"    "LOC101101776" "SC5DL"        "PVRL4"        "LOC100130093" "LINC00338"   
##  [25] "LOC100132891" "PPAP2C"       "C6orf1"       "C2orf47"      "WHSC1L1"      "LOC100289019" "SETD8"        "KDM5B-AS1"   
##  [33] "SPG20"        "CXCR7"        "LOC100216479" "LOC100505761" "MGC57346"     "LPHN3"        "CENPC1"       "C11orf93"    
##  [41] "C14orf169"    "LOC100506060" "FLJ31485"     "LOC440905"    "MLF1IP"       "TMEM194A"     "RRP7B"        "REXO1L1"     
##  [49] "LOC100129269" "KIAA1715"     "CTAGE5"       "LOC202781"    "LOC100506714" "LOC401164"    "UTS2D"        "LOC146880"   
##  [57] "KIAA1804"     "C5orf55"      "C21orf119"    "PRUNE"        "LRRC16A"      "LOC339240"    "FLJ35024"     "C5orf28"     
##  [65] "LOC100505876" "MGC21881"     "LOC100133985" "PPAPDC2"      "FRG1B"        "CECR5"        "LOC100129361" "CCBL1"       
##  [73] "PTPLAD1"      "MST4"         "LOC550112"    "LOC389791"    "CCDC90A"      "KIAA0195"     "LOC100506469" "LOC100133161"
##  [81] "LOC646719"    "LOC728819"    "BRE"          "LOC284581"    "LOC441081"    "LOC728377"    "LOC100134229" "C3orf65"     
##  [89] "SMEK2"        "KIAA1737"     "C17orf70"     "PLEKHM1P"     "LOC338758"    "PCNXL2"       "LOC91948"     "C17orf89"    
##  [97] "LOC100505783" "SMCR7L"       "C8orf4"       "GPR56"        "ATHL1"        "LOC339535"    "PPAPDC1B"     "DAK"         
## [105] "LOC100507173" "CRHR1-IT1"    "PPAP2B"       "ADCK4"        "KIAA0146"     "GYLTL1B"      "LOC100272216" "LOC400027"   
## [113] "WHSC1"        "LOC100130855" "C7orf55"      "C19orf40"     "ADCK3"        "C9orf142"     "SGOL1"        "LOC90834"    
## [121] "PTPLAD2"      "KIAA1967"     "LOC100132352" "LOC100630918" "ADRBK2"       "LINC00263"    "FAM64A"       "LOC401074"   
## [129] "FAM179B"      "RP1-177G6.2"  "METTL21D"     "ERO1LB"       "FLJ45445"     "NADKD1"       "LOC100506233" "LOC100652772"
## [137] "FAM175A"      "LINC00630"    "C11orf82"     "SETD5-AS1"    "SGK196"       "FLJ14186"     "CCDC104"      "FAM63A"      
## [145] "NARG2"        "MTERFD1"      "CCDC74B-AS1"  "LOC286186"    "WDR67"        "C12orf52"     "FLJ30403"     "KIAA2018"    
## [153] "GCN1L1"       "FLJ43681"     "LOC152217"    "FONG"         "C18orf8"      "ALG1L9P"      "GTDC2"        "LOC100507217"
## [161] "NBPF24"       "WBSCR27"      "C14orf1"      "LOC284889"    "KIAA0317"     "FAM65A"       "PMS2L2"       "LUST"        
## [169] "C15orf52"     "FAM195A"      "LOC399744"    "PYCRL"        "LOC338799"    "LOC100506190" "C9orf91"      "FLJ45340"    
## [177] "LOC349196"    "LOC100128881" "TOMM70A"      "ALS2CR8"      "LDOC1L"       "HDGFRP3"      "ZNF767"       "LOC728558"   
## [185] "LOC283693"    "LEPREL2"      "QTRTD1"       "SELM"         "C6orf25"      "C1orf86"      "HNRPLL"       "LOC145820"   
## [193] "LOC100289341" "C17orf85"     "C3orf72"      "C14orf64"     "C9orf9"       "LOC100506394"

length(geneset_niche1)
## [1] 1668
length(geneset_niche2)
## [1] 2889
```

It is always useful to check the number of genes in the geneset before
doing the ligand activity analysis. We recommend having between 20 and
1000 genes in the geneset of interest, and a background of at least 5000
genes for a proper ligand activity analysis. If you retrieve too many DE
genes, it is recommended to use a higher `lfc_cutoff` threshold. We
recommend using a cutoff of 0.15 if you have &gt; 2 receiver
cells/niches to compare and use the min\_lfc as specificity score. If
you have only 2 receivers/niche, we recommend using a higher threshold
(such as using 0.25). If you have single-cell data like Smart-seq2 with
high sequencing depth, we recommend to also use higher threshold.

As we see here, we have Smart-seq2 data and only 2 niches to compare, so
we will use a stronger LFC threshold to keep less DE genes, but more
trustworthy ones.

``` r
lfc_cutoff = 0.75 

specificity_score_targets = "min_lfc"

DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
  
background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  
# Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
# If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
## [1] "ANXA8L2"  "PRKCDBP"  "IL8"      "PTRF"     "SEPP1"    "C1orf186"
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
## [1] "LOC344887"    "AGPAT9"       "C1orf110"     "KIAA1467"     "LOC100292680" "EPT1"         "CT45A4"

length(geneset_niche1)
## [1] 169
length(geneset_niche2)
## [1] 136
```

``` r
top_n_target = 250

niche_geneset_list = list(
  "pEMT_High_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "pEMT_Low_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )
  
ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
## [1] "Calculate Ligand activities for: Malignant_High"
## [1] "Calculate Ligand activities for: Malignant_Low"
```

# 5. Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)

In this step, we will calculate average (scaled) expression, and
fraction of expression, of ligands, receptors, and target genes across
all cell types of interest. Now this is here demonstrated via the
DotPlot function of Seurat, but this can also be done via other ways of
course.

``` r
features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
  
dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl = dotplot$data %>% as_tibble()
exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
```

``` r
exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
```

# 6. Expression fraction and receptor

In this step, we will score ligand-receptor interactions based on
expression strength of the receptor, in such a way that we give higher
scores to the most strongly expressed receptor of a certain ligand, in a
certain celltype. This will not effect the rank of individual ligands
later on, but will help in prioritizing the most important receptors per
ligand (next to other factors regarding the receptor - see later).

``` r
exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 
```

# 7. Prioritization of ligand-receptor and ligand-target links

In this step, we will combine all the above calculated information to
prioritize ligand-receptor-target links. We scale every property of
interest between 0 and 1, and the final prioritization score is a
weighted sum of the scaled scores of all the properties of interest.

We provide the user the option to consider the following properties for
prioritization (of which the weights are defined in
`prioritizing_weights`) :

-   Ligand DE score: niche-specific expression of the ligand: by
    default, this the minimum logFC between the sender of interest and
    all the senders of the other niche(s). The higher the min logFC, the
    higher the niche-specificity of the ligand. Therefore we recommend
    to give this factor a very high weight. `prioritizing_weights`
    argument: `"scaled_ligand_score"`. Recommended weight: 5 (at least
    1, max 5).

-   Scaled ligand expression: scaled expression of a ligand in one
    sender compared to the other cell types in the dataset. This might
    be useful to rescue potentially interesting ligands that have a high
    scaled expression value, but a relatively small min logFC compared
    to the other niche. One reason why this logFC might be small occurs
    when (some) genes are not picked up efficiently by the used
    sequencing technology (or other reasons for low RNA expression of
    ligands). For example, we have observed that many ligands from the
    Tgf-beta/BMP family are not picked up efficiently with single-nuclei
    RNA sequencing compared to single-cell sequencing.
    `prioritizing_weights` argument:
    `"scaled_ligand_expression_scaled"`. Recommended weight: 1 (unless
    technical reason for lower gene detection such as while using
    Nuc-seq: then recommended to use a higher weight: 2).

-   Ligand expression fraction: Ligands that are expressed in a smaller
    fraction of cells of a cell type than defined by
    `exprs_cutoff`(default: 0.10) will get a lower ranking, proportional
    to their fraction (eg ligand expressed in 9% of cells will be ranked
    higher than ligand expressed in 0.5% of cells). We opted for this
    weighting based on fraction, instead of removing ligands that are
    not expressed in more cells than this cutoff, because some
    interesting ligands could be removed that way. Fraction of
    expression is not taken into account for the prioritization if it is
    already higher than the cutoff. `prioritizing_weights` argument:
    `"ligand_fraction"`. Recommended weight: 1.

-   Ligand spatial DE score: spatial expression specificity of the
    ligand. If the niche of interest is at a specific tissue location,
    but some of the sender cell types of that niche are also present in
    other locations, it can be very informative to further prioritize
    ligands of that sender by looking how they are DE between the
    spatial location of interest compared to the other locations.
    `prioritizing_weights` argument: `"scaled_ligand_score_spatial"`.
    Recommended weight: 2 (or 0 if not applicable).

-   Receptor DE score: niche-specific expression of the receptor: by
    default, this the minimum logFC between the receiver of interest and
    all the receiver of the other niche(s). The higher the min logFC,
    the higher the niche-specificity of the receptor. Based on our
    experience, we don’t suggest to give this as high importance as the
    ligand DE, but this might depend on the specific case study.
    `prioritizing_weights` argument: `"scaled_receptor_score"`.
    Recommended weight: 0.5 (at least 0.5, and lower than
    `"scaled_ligand_score"`).

-   Scaled receptor expression: scaled expression of a receptor in one
    receiver compared to the other cell types in the dataset. This might
    be useful to rescue potentially interesting receptors that have a
    high scaled expression value, but a relatively small min logFC
    compared to the other niche. One reason why this logFC might be
    small occurs when (some) genes are not picked up efficiently by the
    used sequencing technology. `prioritizing_weights` argument:
    `"scaled_receptor_expression_scaled"`. Recommended weight: 0.5.

-   Receptor expression fraction: Receptors that are expressed in a
    smaller fraction of cells of a cell type than defined by
    `exprs_cutoff`(default: 0.10) will get a lower ranking, proportional
    to their fraction (eg receptor expressed in 9% of cells will be
    ranked higher than receptor expressed in 0.5% of cells). We opted
    for this weighting based on fraction, instead of removing receptors
    that are not expressed in more cells than this cutoff, because some
    interesting receptors could be removed that way. Fraction of
    expression is not taken into account for the prioritization if it is
    already higher than the cutoff. `prioritizing_weights` argument:
    `"receptor_fraction"`. Recommended weight: 1.

-   Receptor expression strength: this factor let us give higher weights
    to the most highly expressed receptor of a ligand in the receiver.
    This let us rank higher one member of a receptor family if it higher
    expressed than the other members. `prioritizing_weights` argument:
    `"ligand_scaled_receptor_expression_fraction"`. Recommended value: 1
    (minimum: 0.5).

-   Receptor spatial DE score: spatial expression specificity of the
    receptor. If the niche of interest is at a specific tissue location,
    but the receiver cell type of that niche is also present in other
    locations, it can be very informative to further prioritize
    receptors of that receiver by looking how they are DE between the
    spatial location of interest compared to the other locations.
    `prioritizing_weights` argument: `"scaled_receptor_score_spatial"`.
    Recommended weight: 1 (or 0 if not applicable).

-   Absolute ligand activity: to further prioritize ligand-receptor
    pairs based on their predicted effect of the ligand-receptor
    interaction on the gene expression in the receiver cell type -
    absolute ligand activity accords to ‘absolute’ enrichment of target
    genes of a ligand within the affected receiver genes.
    `prioritizing_weights` argument: `"scaled_activity"`. Recommended
    weight: 0, unless absolute enrichment of target genes is of specific
    interest.

-   Normalized ligand activity: to further prioritize ligand-receptor
    pairs based on their predicted effect of the ligand-receptor
    interaction on the gene expression in the receiver cell type -
    normalization of activity is done because we found that some
    datasets/conditions/niches have higher baseline activity values than
    others - normalized ligand activity accords to ‘relative’ enrichment
    of target genes of a ligand within the affected receiver genes.
    `prioritizing_weights` argument: `"scaled_activity_normalized"`.
    Recommended weight: at least 1.

-   Prior knowledge quality of the L-R interaction: the NicheNet LR
    network consists of two types of interactions: L-R pairs documented
    in curated databases, and L-R pairs predicted based on gene
    annotation and PPIs. The former are categorized as ‘bona fide’
    interactions. To rank bona fide interactions higher, but not exlude
    potentially interesting non-bona-fide ones, we give bona fide
    interactions a score of 1, and non-bona-fide interactions a score
    fof 0.5. `prioritizing_weights` argument: `"bona_fide"` Recommend
    weight: at least 1.

``` r
prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                          "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)
```

Note: these settings will give substantially more weight to DE
ligand-receptor pairs compared to activity. Users can change this if
wanted, just like other settings can be changed if that would be better
to tackle the specific biological question you want to address.

``` r
output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
         ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
prioritization_tables = get_prioritization_tables(output, prioritizing_weights)

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
## # A tibble: 10 x 37
##    niche   receiver  sender ligand_receptor ligand receptor bonafide ligand_score ligand_signific~ ligand_present ligand_expressi~
##    <chr>   <chr>     <chr>  <chr>           <chr>  <chr>    <lgl>           <dbl>            <dbl>          <dbl>            <dbl>
##  1 pEMT_H~ Malignan~ T.cel~ PTPRC--MET      PTPRC  MET      FALSE            3.22                1              1             9.32
##  2 pEMT_H~ Malignan~ T.cel~ PTPRC--EGFR     PTPRC  EGFR     FALSE            3.22                1              1             9.32
##  3 pEMT_H~ Malignan~ T.cel~ PTPRC--CD44     PTPRC  CD44     FALSE            3.22                1              1             9.32
##  4 pEMT_H~ Malignan~ T.cel~ PTPRC--ERBB2    PTPRC  ERBB2    FALSE            3.22                1              1             9.32
##  5 pEMT_H~ Malignan~ T.cel~ PTPRC--IFNAR1   PTPRC  IFNAR1   FALSE            3.22                1              1             9.32
##  6 pEMT_H~ Malignan~ T.cel~ TNF--TNFRSF21   TNF    TNFRSF21 TRUE             1.74                1              1             2.34
##  7 pEMT_H~ Malignan~ Myelo~ SERPINA1--LRP1  SERPI~ LRP1     TRUE             2.52                1              1             4.83
##  8 pEMT_H~ Malignan~ Myelo~ IL1B--IL1RAP    IL1B   IL1RAP   TRUE             1.50                1              1             1.93
##  9 pEMT_H~ Malignan~ Myelo~ IL1RN--IL1R2    IL1RN  IL1R2    TRUE             1.62                1              1             2.07
## 10 pEMT_H~ Malignan~ T.cel~ PTPRC--INSR     PTPRC  INSR     FALSE            3.22                1              1             9.32
## # ... with 26 more variables: ligand_expression_scaled <dbl>, ligand_fraction <dbl>, ligand_score_spatial <dbl>,
## #   receptor_score <dbl>, receptor_significant <dbl>, receptor_present <dbl>, receptor_expression <dbl>,
## #   receptor_expression_scaled <dbl>, receptor_fraction <dbl>, receptor_score_spatial <dbl>,
## #   ligand_scaled_receptor_expression_fraction <dbl>, avg_score_ligand_receptor <dbl>, activity <dbl>, activity_normalized <dbl>,
## #   scaled_ligand_score <dbl>, scaled_ligand_expression_scaled <dbl>, scaled_receptor_score <dbl>,
## #   scaled_receptor_expression_scaled <dbl>, scaled_avg_score_ligand_receptor <dbl>, scaled_ligand_score_spatial <dbl>,
## #   scaled_receptor_score_spatial <dbl>, scaled_ligand_fraction_adapted <dbl>, scaled_receptor_fraction_adapted <dbl>, ...
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
## # A tibble: 10 x 20
##    niche           receiver   sender  ligand_receptor ligand receptor bonafide target target_score target_signific~ target_present
##    <chr>           <chr>      <chr>   <chr>           <chr>  <chr>    <lgl>    <chr>         <dbl>            <dbl>          <dbl>
##  1 pEMT_High_niche Malignant~ T.cell~ PTPRC--MET      PTPRC  MET      FALSE    EHF           1.04                 1              1
##  2 pEMT_High_niche Malignant~ T.cell~ PTPRC--MET      PTPRC  MET      FALSE    GADD4~        0.836                1              1
##  3 pEMT_High_niche Malignant~ T.cell~ PTPRC--MET      PTPRC  MET      FALSE    SERPI~        0.889                1              1
##  4 pEMT_High_niche Malignant~ T.cell~ PTPRC--EGFR     PTPRC  EGFR     FALSE    EHF           1.04                 1              1
##  5 pEMT_High_niche Malignant~ T.cell~ PTPRC--EGFR     PTPRC  EGFR     FALSE    GADD4~        0.836                1              1
##  6 pEMT_High_niche Malignant~ T.cell~ PTPRC--EGFR     PTPRC  EGFR     FALSE    SERPI~        0.889                1              1
##  7 pEMT_High_niche Malignant~ T.cell~ PTPRC--CD44     PTPRC  CD44     FALSE    EHF           1.04                 1              1
##  8 pEMT_High_niche Malignant~ T.cell~ PTPRC--CD44     PTPRC  CD44     FALSE    GADD4~        0.836                1              1
##  9 pEMT_High_niche Malignant~ T.cell~ PTPRC--CD44     PTPRC  CD44     FALSE    SERPI~        0.889                1              1
## 10 pEMT_High_niche Malignant~ T.cell~ PTPRC--ERBB2    PTPRC  ERBB2    FALSE    EHF           1.04                 1              1
## # ... with 9 more variables: target_expression <dbl>, target_expression_scaled <dbl>, target_fraction <dbl>,
## #   ligand_target_weight <dbl>, activity <dbl>, activity_normalized <dbl>, scaled_activity <dbl>,
## #   scaled_activity_normalized <dbl>, prioritization_score <dbl>

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
## # A tibble: 10 x 37
##    niche   receiver sender  ligand_receptor ligand receptor bonafide ligand_score ligand_signific~ ligand_present ligand_expressi~
##    <chr>   <chr>    <chr>   <chr>           <chr>  <chr>    <lgl>           <dbl>            <dbl>          <dbl>            <dbl>
##  1 pEMT_L~ Maligna~ Endoth~ F8--LRP1        F8     LRP1     TRUE            0.952              1                1            2.17 
##  2 pEMT_L~ Maligna~ Endoth~ PLAT--LRP1      PLAT   LRP1     TRUE            0.913              1                1            2.70 
##  3 pEMT_L~ Maligna~ CAF_Low FGF10--FGFR2    FGF10  FGFR2    TRUE            0.385              0.8              1            1.07 
##  4 pEMT_L~ Maligna~ CAF_Low NLGN2--NRXN3    NLGN2  NRXN3    TRUE            0.140              0.2              1            0.269
##  5 pEMT_L~ Maligna~ CAF_Low RSPO3--LGR6     RSPO3  LGR6     TRUE            0.557              0.8              1            1.27 
##  6 pEMT_L~ Maligna~ CAF_Low COMP--SDC1      COMP   SDC1     TRUE            0.290              0.8              1            1.27 
##  7 pEMT_L~ Maligna~ CAF_Low SEMA3C--NRP2    SEMA3C NRP2     TRUE            0.652              1                1            1.73 
##  8 pEMT_L~ Maligna~ CAF_Low SLIT2--SDC1     SLIT2  SDC1     TRUE            0.494              1                1            0.846
##  9 pEMT_L~ Maligna~ Endoth~ IL33--IL1RAP    IL33   IL1RAP   FALSE           1.34               1                1            2.75 
## 10 pEMT_L~ Maligna~ CAF_Low C3--LRP1        C3     LRP1     TRUE            0.480              1                1            4.79 
## # ... with 26 more variables: ligand_expression_scaled <dbl>, ligand_fraction <dbl>, ligand_score_spatial <dbl>,
## #   receptor_score <dbl>, receptor_significant <dbl>, receptor_present <dbl>, receptor_expression <dbl>,
## #   receptor_expression_scaled <dbl>, receptor_fraction <dbl>, receptor_score_spatial <dbl>,
## #   ligand_scaled_receptor_expression_fraction <dbl>, avg_score_ligand_receptor <dbl>, activity <dbl>, activity_normalized <dbl>,
## #   scaled_ligand_score <dbl>, scaled_ligand_expression_scaled <dbl>, scaled_receptor_score <dbl>,
## #   scaled_receptor_expression_scaled <dbl>, scaled_avg_score_ligand_receptor <dbl>, scaled_ligand_score_spatial <dbl>,
## #   scaled_receptor_score_spatial <dbl>, scaled_ligand_fraction_adapted <dbl>, scaled_receptor_fraction_adapted <dbl>, ...
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
## # A tibble: 10 x 20
##    niche          receiver   sender   ligand_receptor ligand receptor bonafide target target_score target_signific~ target_present
##    <chr>          <chr>      <chr>    <chr>           <chr>  <chr>    <lgl>    <chr>         <dbl>            <dbl>          <dbl>
##  1 pEMT_Low_niche Malignant~ Endothe~ F8--LRP1        F8     LRP1     TRUE     ETV4          0.771                1              1
##  2 pEMT_Low_niche Malignant~ Endothe~ PLAT--LRP1      PLAT   LRP1     TRUE     CLDN7         0.835                1              1
##  3 pEMT_Low_niche Malignant~ Endothe~ PLAT--LRP1      PLAT   LRP1     TRUE     ETV4          0.771                1              1
##  4 pEMT_Low_niche Malignant~ CAF_Low  FGF10--FGFR2    FGF10  FGFR2    TRUE     ETV4          0.771                1              1
##  5 pEMT_Low_niche Malignant~ CAF_Low  FGF10--FGFR2    FGF10  FGFR2    TRUE     WNT5A         1.40                 1              1
##  6 pEMT_Low_niche Malignant~ CAF_Low  NLGN2--NRXN3    NLGN2  NRXN3    TRUE     CLDN5         0.979                1              1
##  7 pEMT_Low_niche Malignant~ CAF_Low  NLGN2--NRXN3    NLGN2  NRXN3    TRUE     ETV4          0.771                1              1
##  8 pEMT_Low_niche Malignant~ CAF_Low  RSPO3--LGR6     RSPO3  LGR6     TRUE     DDC           0.832                1              1
##  9 pEMT_Low_niche Malignant~ CAF_Low  RSPO3--LGR6     RSPO3  LGR6     TRUE     EGFL7         0.763                1              1
## 10 pEMT_Low_niche Malignant~ CAF_Low  COMP--SDC1      COMP   SDC1     TRUE     CLDN7         0.835                1              1
## # ... with 9 more variables: target_expression <dbl>, target_expression_scaled <dbl>, target_fraction <dbl>,
## #   ligand_target_weight <dbl>, activity <dbl>, activity_normalized <dbl>, scaled_activity <dbl>,
## #   scaled_activity_normalized <dbl>, prioritization_score <dbl>
```

# 8. Visualization of the Differential NicheNet output

## Differential expression of ligand and expression

Before visualization, we need to define the most important
ligand-receptor pairs per niche. We will do this by first determining
for which niche the highest score is found for each
ligand/ligand-receptor pair. And then getting the top 50 ligands per
niche.

``` r
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche
```

Now we will look first at the top ligand-receptor pairs for KCs (here,
we will take the top 2 scoring receptors per prioritized ligand)

``` r
receiver_oi = "Malignant_High" 

filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
```

Visualization: minimum LFC compared to other niches

``` r
lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

Show the spatialDE as additional information

``` r
lfc_plot_spatial = make_ligand_receptor_lfc_spatial_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, ligand_spatial = include_spatial_info_sender, receptor_spatial = include_spatial_info_receiver, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot_spatial
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## Ligand expression, activity and target genes

Active target gene inference - cf Default NicheNet

Now: visualization of ligand activity and ligand-target links

``` r
exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_activity_target_plot$combined_plot
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->
Based on this plot, we can infer many hypotheses such as the following:
“Interestingly, IL1 family ligands seem to have activity in inducing the
DE genes between high pEMT and low pEMT malignant cells; and they are
mainly expressed by myeloid cells, a cell type unique for pEMT-high
tumors.”

**important: ligand-receptor pairs with both high differential
expression (or condition-specificity) and ligand activity (=target gene
enrichment) are very interesting predictions as key regulators of your
intercellular communication process of interest !**

If this plot contains too much information because we look at many hits
(top 50 ligands), you can make this plot of course for less ligands as
well, eg for the top20.

``` r
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(20, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_activity_target_plot$combined_plot
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

## Circos plot of prioritized ligand-receptor pairs

Because a top50 is too much to visualize in a circos plot, we will only
visualize the top 15.

``` r
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(15, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% length(), name = 'Spectral') %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->

## Interpretation of these results

Most top-ranked differential L-R pairs seem to come from the cell types
that are only present in the pEMT-high tumors. This might be partially
due to biology (unique cell types in one condition, are likely to be
very important), but might also be due to the way of prioritizing and
the fact that those unique cell types don’t have a ‘counterpart’ in the
other niche(s).

Because myeloid cells and T cells are very different from the other
cells in the tumor microenvironment, their ligands will show strong
differential expression. This differential expression (myeloid/tcell vs
myofibroblasts/CAFs/Endothelial cells in low-pEMT) is likely to be more
pronounced compared to differential expression between cells from the
same cell type but different niche/condition (CAF in pEMT-high vs CAF in
pEMT-low). So conclusion: it is an advantage of Differential NicheNet
that it can cope with condition-specifc cell types, but the user should
be aware that the final general score might be biased towards
condition-specific sender cell types. Therefore we suggest to also have
a look at the top LR pairs per sender cell type (as we did here for the
first figures) if you have a case study in which some sender cell types
are condition-specific.

## Visualization for the other condition: pEMT-low

``` r
receiver_oi = "Malignant_Low"  
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(50, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
```

![](differential_nichenet_pEMT_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

# Notes, limitations, and comparison to default NicheNet.

In the default NicheNet pipeline, expressed ligand-receptor pairs are
prioritized based on their ligand activity alone. Here, in the
Differential NicheNet pipeline, we also draw information based on
differential expression of the L-R pairs compared to other niches (and
if applicable: other spatial locations.)

Because we here focus on differential expression of ligand-receptor
pairs, and by using the default prioritizations weights more on DE than
activity, we tend to find many different hits than with the default
NicheNet pipeline. With Differential NicheNet, we tend to find more
high-DE, low-activity hits, whereas with default NicheNet we find more
low-DE, high-activity hits.

It should be noted that some of the high-DE, low-activity hits might be
really important because they just have low NicheNet activities due to
limitations in the NicheNet activity prediction (eg improper/incomplete
prior knowledge within NicheNet for that ligand), but some of them might
also be high in DE but not in activity because they don’t have strong
signaling effects (eg ligands involved in cell adhesion only).

For the opposite pairs with low-DE and high-activity that are not
strongly prioritized by Differential NicheNet, the following should be
considered: 1) some ligands are regulated post-transcriptionally, and
that the high predicted activities might still reflect true signaling;
2) high predicted activity values might be due to limitations of
NicheNet (inaccurate prior knowledge) and these lowDE ligands are not
important in the biological process of interest (although a highDE
family member of this ligand may! since signaling between family members
tends to be very similar); 3) high activity in one condition might be
due to downregulation in the other condition, leading to high activity
but low DE. Currently, ligand activities are automatically calculated on
upregulated genes per condition, but downregulated genes could also be a
sign of ligand activity.

When Ligand-Receptor pairs have both high DE and high activity, we can
consider them to be very good candidates in regulating the process of
interest, and we recommend testing these candidates for further
experimental validation.

# References

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<doi:10.1038/s41592-019-0667-5>

Guilliams et al. Spatial proteogenomics reveals distinct and
evolutionarily conserved hepatic macrophage niches. Cell (2022)
<doi:10.1016/j.cell.2021.12.018>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-puram_single-cell_2017" class="csl-entry">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
Parameter optimization via mlrMBO
================
Robin Browaeys
2018-02-20

<!-- github markdown built using 
rmarkdown::render("vignettes/parameter_optimization.Rmd", output_format = "github_document") # please, don't run this!!
-->

This vignette shows how we optimized both hyperparameters and data source weights via model-based optimization (see manuscript for more information). Because the optimization requires intensive parallel computation, we performed optimization in parallel on a gridengine cluster via the qsub package (https://cran.r-project.org/web/packages/qsub/qsub.pdf). This script is merely illustrative and should be adapted by the user to work on its own system.

The input data used in this vignette can be found at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758). 

First, we will load in the required packages and networks we will use to construct the models which we will evaluate during the optimization procedure.
```{r}
library(nichenetr)
library(tidyverse)
library(qsub)
library(mlrMBO)

# in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
```

We will load in the ligand treatment validation datasets and try to optimize the parameters to maximize both target gene and ligand activity prediction. In this vignette, we do the optimization on all datasets. Alternatively, you can select a specific subset of datasets and evaluate the final performance on the left-out datasets.

```{r}
# The ligand treatment expression datasets used for validation can be downloaded from Zenodo:
expression_settings_validation = readRDS(url("https://zenodo.org/record/3260758/files/expression_settings.rds"))
```

Define the optimization wrapper function and config information for the qsub package
```{r}
mlrmbo_optimization_wrapper = function(...){
  library(nichenetr)
  library(mlrMBO)
  library(parallelMap)
  library(dplyr)
  output = mlrmbo_optimization(...)
  return(output)
}

qsub_config = create_qsub_config(
  remote = "myuser@mycluster.address.org:1234",
  local_tmp_path = "/tmp/r2gridengine",
  remote_tmp_path = "/scratch/personal/myuser/r2gridengine",
  num_cores = 24,
  memory = "10G",
  wait = FALSE, 
  max_wall_time = "500:00:00"
)
```

Perform optimization:

```{r}
additional_arguments_topology_correction = list(source_names = source_weights_df$source %>% unique(), 
                                                algorithm = "PPR", 
                                                correct_topology = FALSE,
                                                lr_network = lr_network, 
                                                sig_network = sig_network, 
                                                gr_network = gr_network, 
                                                settings = lapply(expression_settings_validation,convert_expression_settings_evaluation), 
                                                secondary_targets = FALSE, 
                                                remove_direct_links = "no", 
                                                cutoff_method = "quantile")
nr_datasources = additional_arguments_topology_correction$source_names %>% length()

obj_fun_multi_topology_correction = makeMultiObjectiveFunction(name = "nichenet_optimization",
                                                               description = "data source weight and hyperparameter optimization: expensive black-box function", 
                                                               fn = model_evaluation_optimization, 
                                                               par.set = makeParamSet(
                                                                 makeNumericVectorParam("source_weights", len = nr_datasources, lower = 0, upper = 1), 
                                                                 makeNumericVectorParam("lr_sig_hub", len = 1, lower = 0, upper = 1),  
                                                                 makeNumericVectorParam("gr_hub", len = 1, lower = 0, upper = 1),  
                                                                 makeNumericVectorParam("ltf_cutoff", len = 1, lower = 0.9, upper = 0.999),  
                                                                 makeNumericVectorParam("damping_factor", len = 1, lower = 0.01, upper = 0.99)), 
                                                               has.simple.signature = FALSE,
                                                               n.objectives = 4, 
                                                               noisy = FALSE,
                                                               minimize = c(FALSE,FALSE,FALSE,FALSE))
set.seed(1)

# Run with: 50 iterations, 24 desings evaluated in parallel and 240 start designs

job_mlrmbo = qsub_lapply(X = 1,
                               FUN = mlrmbo_optimization_wrapper,
                               object_envir = environment(mlrmbo_optimization_wrapper),
                               qsub_config = qsub_config,
                               qsub_environment = NULL, 
                               qsub_packages = NULL,
                               obj_fun_multi_topology_correction, 
                               50, 24, 240, 
                               additional_arguments_topology_correction)
```
Once the job is finised (which can take a few days - for shorter running time: reduce the number of iterations), run:

```{r}
res_job_mlrmbo = qsub_retrieve(job_mlrmbo)
```

Get now the most optimal parameter setting as a result of this analysis

```{r}
optimized_parameters = res_job_mlrmbo %>% process_mlrmbo_nichenet_optimization(source_names = source_weights_df$source %>% unique())
```

When you would be interested to generate a context-specific model, it could be possible that you would like to optimize the parameters specifically on your dataset of interest and not on the general ligand treatment datasets (be aware for overfitting, though!). Because for your own data, you don't know the true active ligands, you could only optimize target gene prediction performance and not ligand activity performance. In order to this, you would need to change the expression settings in the optimization functions such that they include your data, and use the function `model_evaluation_optimization_application` instead of `model_evaluation_optimization` (define this as the function parameter in `makeMultiObjectiveFunction` shown here above).
