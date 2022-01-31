# PHONEMeS

**PHONEMeS** (**PHO**sphorylation **NE**tworks for **M**ass **S**pectrometry) is a method to model signalling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions. 
The newest version of PHONEMeS builds on `CARNIVAL`, and more specifically the `runCARNIVAL`wrapper function. As such, it requires the CARNIVAL package (Version 1.3) to be installed. Additionally, it uses `OmnipathR` to construct prior knowledge networks based on solely kinase-substrate interactions or the combination of kinase-substrate and protein-protein interactions.
The pipeline requires a vector of input nodes (perturbed kinases) and measured nodes (deregulated phosphosites) that are connected by an interactive version of 
IBM Cplex or CBC-COIN solver. The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio). The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available
for any user. Alternatively for smaller cases, users can rely on the freely available 
[lpSolve R-package](https://cran.r-project.org/web/packages/lpSolve/index.html). 


### License

Distributed under the GNU GPLv2 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/PHONEMeS/blob/master/LICENSE.txt) or copy at https://www.gnu.org/licenses/gpl-2.0.html.

### Installation

To install `PHONEMeS` please run:
```
devtools::install_github('saezlab/PHONEMeS')
```

### Usage

A short tutorial on how to run a PHONEMeS analysis using PHONEMeS v2.0.0 is provided
in the [vignette tutorial](https://github.com/saezlab/PHONEMeS/blob/master/vignettes/tutorial.Rmd).


### Prior versions

The code for the original PHONEMeS package (PHONEMeS v1.0.0) as described in Terfve et al. 2015 can be found in the releases.
For a guide how to run a PHONEMeS analysis using PHONEMeS v1.0.0, please refer to the [documentation](https://saezlab.github.io/PHONEMeS).


### References

[Terfve et al.](http://www.nature.com/articles/ncomms9033):

> Terfve, C. D. A., Wilkes, E. H., Casado, P., Cutillas, P. R., and Saez-Rodriguez, J. (2015). Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data. *Nature Communications*, 6:8033.

[Wilkes et al.](http://www.pnas.org/content/112/25/7719.abstract) (description of parts of the data)

> Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). Empirical inference of circuitry and plasticity in a kinase signaling network. *Proceedings of the National Academy of Sciences of the United States of America,* 112(25):7719–24.

---
title: "PHONEMeS"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{PHONEMeS}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(PHONEMeS)
```
---
title: "tutorial"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction
This is a brief PHONEMeS tutorial. It provides an example for an analysis with 
PHONEMeS starting from the results of an differential analysis of phosphorylation 
sites, in this particular case a top table results object from limma.
PHONEMeS is a method to model signaling networks based on untargeted phosphoproteomics
mass spectrometry data and kinase/phosphatase-substrate interactions. It identifies
a subset of interactions from a prior knowledge network that represent potential 
regulatory signaling pathways linking known or potentially deregulated kinases to 
downstream phorsphorylation sites.

# Getting set up
First, we load the PHONEMeS library. To use PHONEMeS, in particular the run_phonemes
function, an interactive version of IBM Cplex or CBC-COIN solver is required as 
the network optimizer. The IBM ILOG Cplex is freely available through Academic 
Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio). 
The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available
for any user. Alternatively for smaller cases, users can rely on the freely available 
[lpSolve R-package](https://cran.r-project.org/web/packages/lpSolve/index.html). 

Additionally, we load the [decoupleR](https://github.com/saezlab/decoupleR) library,
a recent development from the Saez-Rodriguez group which offers 12 different 
statistical approaches to perform functional analysis of omics data.

```{r setup}
library(PHONEMeS)
library(decoupleR)
library(tidyverse)
```

We then take a look at the tutorialData inside PHONEMeS. It consists of a top 
table results object from limma coming from phosphoproteomics. This contains the 
information about phosphorylation sites that are differntially abundant between 
two conditions. 

```{r}
head(tutorialData)
```

# Input preparation
For the extraction of perturbed signaling networks using PHONEMeS, potentially 
deregulated phosphorylation sites and kinases need to be selected that will be 
connected based on a prior knowledge network. 

## Selection of deregulated phosphorylation sites
Deregulated phosphorylation sites can be selected based on their changes detected
within the differential analysis. We choose T values for the selection over 
other metrices (p-values, logFC) as it combines the significance and magnitude 
of the change.
We then choose the top 10 deregulated phosphorylation sites that we will try to 
connect to upstream kinases. This is in arbitrary number and users should decide 
on a case to case basis what selection criteria to apply. The phosphorylation 
sites could also be chosen based on a chosen percentage or score cut-off.

```{r}
top_pps <- tutorialData %>% 
  dplyr::filter(ID %in% phonemesPKN$target) %>%
  dplyr::arrange(dplyr::desc(base::abs(t))) %>%
  dplyr::slice(1:10)

deregulated_pps <- top_pps$t
names(deregulated_pps) <- top_pps$ID
```


## Selection of deregulated kinases
Kinases that are deregulated can either be selected manually, e.g. if a specific 
kinase is known to be perturbed in an experiment, or based on changes in their 
activity. Therefore, we will estimate kinase activity changes using the decoupleR 
framework.


### Manual selection
For the manual selection, we create a names vector were kinases can be classified 
as eiter up-regulated (1) or down-regulated (-1). In this particular example we 
assume an up-regulation of AKT1.

```{r}
deregulated_kinases_man <- c(AKT1 = 1)
```

### Kinase activity estimation
DecoupleR is a collection of computational methods used for footprint-based analysis
that assume omics data as a signature of upstream biological activities. Here, we
can use the changes in phosphorylation sites as signatures for the activity of 
upstream kinases. The interactions between phosphorylation sites and kinases are 
again taken from the phonemesPKN.
To use the phonemesPKN for decoupleR we need to add a column for the likelihood 
of each interaction. Since our network is not weighted, we just add a 1 for all 
interactions.

```{r}
decoupler_network <- phonemesKSN %>% dplyr::rename("mor" = interaction) %>% tibble::add_column("likelihood" = 1)
```

As decoupleR functions expect a matrix as input, we are going to create a one-column
matrix from the tutorialData object. This one column matrix contains the vector 
of T values obtained from limma. We again choose the T values over other metrices.

```{r}
decoupler_input <- tutorialData %>% 
  dplyr::filter(ID %in% decoupler_network$target) %>%
  tibble::column_to_rownames("ID") %>% 
  dplyr::select(t)
```

Before running any method in decoupleR, we need intersect our prior knowledge 
network with the input matrix. With that we filter out regulons with less than 
5 target features.

```{r}
decoupler_network <- decoupleR::intersect_regulons(mat = decoupler_input, network = decoupler_network, 
                                                   .source = source, .target = target, minsize = 5)
```

Additionally, the overlapp of different regulons should be checked. If some kinases
have the same set of targets (correlatiom = 1), we only keep one of them.
```{r}
correlated_regulons <- decoupleR::check_corr(decoupler_network) %>% dplyr::filter(correlation == 1)
decoupler_network <- decoupler_network %>% dplyr::filter(!source %in% correlated_regulons$source.2)
```

We use mlm to estimate kinase activity but if the user would like to check other
available methods they can run show_methods(). All methods follow the same design
pattern and arguments.

```{r}
kinase_activity <- decoupleR::run_mlm(mat = decoupler_input, network = decoupler_network)
```

We then choose the top 5 deregulated kinases that we will try to connect to the 
downstream phosphorylation sites. This is again in arbitrary number and users 
should decide on a case to case basis what selection criteria to apply. Kinases 
could also be chosen based on a chosen percentage or score cut-off.
As before, kinases will eiter be classified as up-regulated (1) or down-regulated (-1).

```{r}
top_kinases <- kinase_activity %>% dplyr::arrange(dplyr::desc(base::abs(score))) %>% dplyr::slice(1:5)

deregulated_kinases <- top_kinases$score
names(deregulated_kinases) <- top_kinases$source

deregulated_kinases[deregulated_kinases > 0] <- 1
deregulated_kinases[deregulated_kinases < 0] <- -1
```

Moreover, we recommend the users to select kinases that do not appear to be regulated
based on their kinase activities. We assume that kinases with a score below 0.5 
are likely to not be deregulated between the conditions and with that do not explain
the observed changes. These kinases will be filtered out from our network within PHONEMeS.

```{r}
nc_kinases <- kinase_activity %>% dplyr::filter(base::abs(score) <= 0.5) %>% dplyr::pull(source)
```

# Run PHONEMeS
The deregulated phosphorylation sites and kinases can now be connected using PHONEMeS,
in specific the run_phonemes function. Within that function the prior knowledge 
network is pruned by removing nodes 50 steps upstream and downstream of measurements
and inputs. Remember that the path to an interactive version of IBM Cplex or CBC-COIN
solver is required.

This example should take about 1 minute to be solved by PHONEMeS.

```{r, results = "hide"}
phonemes_result <- PHONEMeS::run_phonemes(inputObj = deregulated_kinases, 
                                          measObj = deregulated_pps, 
                                          rmNodes = nc_kinases, 
                                          netObj = phonemesPKN,
                                          solverPath = "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/cplex",  # Example path
                                          solver = "cplex")
```
The output of PHONEMeS contains the network, measurements and inputs used for the construction and the final results, including individual solutions and one summary solution (weightedSIF, nodesAttributes). These can be directly used to identify causal interactions between the perturbed nodes and the selected kinases. In addition to extracting direct information from the network, we can run different downstream analysis based on the necessities of each project, e.g. Pathway enrichment analysis.

```{r}
utils::str(phonemes_result,max.level = 2)
```


# Session Info
This tutorial was run on the date specified below.
```{r}
Sys.Date()
```

The sessionInfo() at run time was:
```{r}
sessionInfo()
```



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_phonemes.R
\name{run_phonemes}
\alias{run_phonemes}
\title{Run PHONEMeS}
\usage{
run_phonemes(
  inputObj,
  measObj,
  netObj = phonemesPKN,
  rmNodes = NULL,
  pruning = TRUE,
  n_steps_pruning = 50,
  solverPath,
  solver = "cplex",
  timelimit = 7200,
  mipGAP = 0.05,
  poolrelGAP = 0
)
}
\arguments{
\item{inputObj}{named vector of perturbation targets. Either 1 (up regulated) or -1 (down regulated)}

\item{measObj}{named vector of the measurements}

\item{netObj}{data frame of the prior knowledge network}

\item{rmNodes}{character vector of nodes to remove from prior knowledge network}

\item{pruning}{logic, set to TRUE if network should be pruned (recommended)}

\item{n_steps_pruning}{integer giving the order of the neighborhood}

\item{solverPath}{path to the solver}

\item{solver}{one of the solvers available from getSupportedSolvers()}

\item{timelimit}{solver time limit in seconds}

\item{mipGAP}{CPLEX parameter: absolute tolerance on the gap}

\item{poolrelGAP}{CPLEX/Cbc parameter: Allowed relative gap of accepted}
}
\value{
List of CARNIVAL results and final inputObj, measObj, netObj used
}
\description{
This function runs CARNIVAL with the input of phosphoproteomic data (phosphosites and kinases).
The prior knowledge network used is the combination of protein-protein and protein-phosphosite
interactions from omnipath. Before running CARNIVAL the network is pruned by removing nodes n_steps
upstream and downstream of measurements and inputs, respectively.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{process_phosphositeplus_regpsites}
\alias{process_phosphositeplus_regpsites}
\title{process_phosphositeplus_regpsites}
\usage{
process_phosphositeplus_regpsites(ppsp_regpsites_file_path)
}
\arguments{
\item{regpsite_file}{path to the tab file of regulatory psites of phosphositeplus}
}
\value{
a formatted table of regulatory psites of phosphositeplus
with their known mode of regulation
}
\description{
This function process the tab file of regulatory psites of phosphositeplus,
downloaded from https://www.phosphosite.org/staticDownloads
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vanilla_phonemes.R
\name{run_vanilla_phonemes}
\alias{run_vanilla_phonemes}
\title{Run vanilla PHONEMeS (PHONEMeS v1.0.0)}
\usage{
run_vanilla_phonemes(
  inputObj,
  measObj,
  netObj = phonemesKSN,
  rmNodes = NULL,
  pruning = FALSE,
  n_steps_pruning = 50,
  solverPath,
  solver = "cplex",
  timelimit = 7200,
  mipGAP = 0.05,
  poolrelGAP = 0
)
}
\arguments{
\item{inputObj}{named vector of perturbation targets}

\item{measObj}{named vector of the measurements}

\item{netObj}{data frame of the prior knowledge network}

\item{rmNodes}{character vector of nodes to remove from prior knowledge network}

\item{pruning}{logic, set to TRUE if network should be pruned (default = FALSE)}

\item{n_steps_pruning}{integer giving the order of the neighborhood}

\item{solverPath}{path to the solver}

\item{solver}{one of the solvers available from getSupportedSolvers()}

\item{timelimit}{solver time limit in seconds}

\item{mipGAP}{CPLEX parameter: absolute tolerance on the gap}

\item{poolrelGAP}{CPLEX/Cbc parameter: Allowed relative gap of accepted}
}
\value{
List of CARNIVAL results and final inputObj, measObj, netObj used
}
\description{
The code for this function is comparable to the orginal PHONEMeS package. An unsigned
kinase-substrate network is used as prior knowledge connecting perturbed kinases with phosphosites.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{phonemesPKN}
\alias{phonemesPKN}
\title{The main PKN used in Phonemes to link deregulated phosphosites to perturbed kinases.
It combines kinase-substrate and protein-protein interactions.}
\format{
The full PKN contains 56,738 interactions: 36,915 kinase-substrate interactions (1,688 kinases; 15,385 phosphosites); 19,823 protein-protein interactions
\describe{
\item{source}{upstream protein}
\item{interaction}{activation/inhibition (1/-1)}
\item{target}{downstream protein/phosphorylation site}
}
}
\source{
\url{https://github.com/saezlab/PHONEMeS/blob/master/data-raw/phonemesPKN.R}
}
\description{
The main PKN used in Phonemes to link deregulated phosphosites to perturbed kinases.
It combines kinase-substrate and protein-protein interactions.
}
\examples{
data("phonemesPKN")
}
\keyword{PKN}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{phonemesKSN}
\alias{phonemesKSN}
\title{The KSN used as prior knowledge network when running vanilla Phonemes.
The KSN contains kinase-substrate interactions via a phosphorylation site}
\format{
The full PKN contains 37,631 kinase-substrate interactions via a phosphorylation site (1,710 kinases; 15,413 phosphosites)
\describe{
\item{source}{upstream kinase/phosphorylation site}
\item{interaction}{connection (1)}
\item{target}{downstream phosphorylation site/kinase}
}
}
\source{
\url{https://github.com/saezlab/PHONEMeS/blob/master/data-raw/phonemesKSN.R}
}
\description{
The KSN used as prior knowledge network when running vanilla Phonemes.
The KSN contains kinase-substrate interactions via a phosphorylation site
}
\examples{
data("phonemesKSN")
}
\keyword{KSN}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_phonemes.R
\name{extract_subnetwork}
\alias{extract_subnetwork}
\title{extract_subnetwork}
\usage{
extract_subnetwork(phonemes_res, targets, n_steps = 3, mode = "all")
}
\arguments{
\item{phonemes_res}{Phonemes result from the run_phonemes function}

\item{targets}{Network nodes, starting point for the extraction of the sub network}

\item{n_steps}{Number of steps to extract down- or upstream of targets}

\item{mode}{Character constant to specify direction of the extraction. "In" for upstream nodes, "out" for downstream nodes and "all" for both.}
}
\value{
Phonemes sub network
}
\description{
This function extracts smaller sub networks from the run_phonemes output
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{regulatory_psites}
\alias{regulatory_psites}
\title{regulatory_psites}
\usage{
regulatory_psites(phonemes_res, th_act = 0)
}
\arguments{
\item{phonemes_res}{carnival result from the run_carnival function}

\item{th_act}{threshold for node activity to include in the benchmark}
}
\value{
list with two elements: the regulatory psites with predicted mode of
regulations and the kinases that catalyse their phosphorilation
}
\description{
This function  extract regulatory phosphites from the phonemes network, e.i.
psites that are differentially regulated on proteins that are found in
the phonemes network
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_phonemes.R
\name{reattach_psites}
\alias{reattach_psites}
\title{Reattach_psites}
\usage{
reattach_psites(phonemes_res)
}
\arguments{
\item{phonemes_res}{phonemes result from the run_phonemes function}
}
\value{
List of PHONEMES results and final inputObj, measObj, netObj used, with psites attached
}
\description{
This function readd links between phosphosite and their correpsonding proteins
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{benchmark_regulatory_psites}
\alias{benchmark_regulatory_psites}
\title{benchmark_regulatory_psites}
\usage{
benchmark_regulatory_psites(regulatory_psites)
}
\arguments{
\item{sif}{regulatory_psites, e.i. the first element of the list returned by
regulatory_psites function}
}
\value{
a table comparing the predicted regulatory psites with the regulatory
psites of phosphositeplus
}
\description{
This function check the coherence of regulatory phosphosites predicted by
phonemes with ground truth from phosphositeplus
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_phonemes.R
\name{get_protein_network}
\alias{get_protein_network}
\title{get_protein_network}
\usage{
get_protein_network(phonemes_res)
}
\arguments{
\item{phonemes_res}{phonemes result from the run_phonemes function}
}
\value{
Phonemes network only consisting of protein-protein interactions
}
\description{
This function readd links between phosphosite and their correpsonding proteins
}
