<img src='https://github.com/RajLabMSSM/echolocatoR/raw/master/inst/hex/hex.png' height='300'><br><br>
[![](https://img.shields.io/badge/devel%20version-0.2.2-black.svg)](https://github.com/RajLabMSSM/echolocatoR)
[![R build
status](https://github.com/RajLabMSSM/echolocatoR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RajLabMSSM/echolocatoR/actions)
[![](https://img.shields.io/github/last-commit/RajLabMSSM/echolocatoR.svg)](https://github.com/RajLabMSSM/echolocatoR/commits/master)
[![](https://codecov.io/gh/RajLabMSSM/echolocatoR/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echolocatoR)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btab658-blue.svg)](https://doi.org/10.1093/bioinformatics/btab658)
<h5>
Author: <i>Brian M. Schilder</i>
</h5>
<h5>
README updated: <i>Mar-04-2022</i>
</h5>

<center>
<h1>
) ) ) ) ))) :bat: echolocatoR :bat: ((( ( ( ( (
</h1>
</center>
<h3>
Automated statistical and functional fine-mapping with extensive access
to genome-wide datasets
</h3>
<hr>

If you use `echolocatoR`, please cite:

> Brian M Schilder, Jack Humphrey, Towfique Raj (2021) echolocatoR: an
> automated end-to-end statistical and functional genomic fine-mapping
> pipeline, *Bioinformatics*; btab658,
> <https://doi.org/10.1093/bioinformatics/btab658>

## Introduction

Fine-mapping methods are a powerful means of identifying causal variants
underlying a given phenotype, but are underutilized due to the technical
challenges of implementation. ***echolocatoR*** is an R package that
automates end-to-end genomics fine-mapping, annotation, and plotting in
order to identify the most probable causal variants associated with a
given phenotype.

It requires minimal input from users (a GWAS or QTL summary statistics
file), and includes a suite of statistical and functional fine-mapping
tools. It also includes extensive access to datasets (linkage
disequilibrium panels, epigenomic and genome-wide annotations, QTL).

The elimination of data gathering and preprocessing steps enables rapid
fine-mapping of many loci in any phenotype, complete with locus-specific
publication-ready figure generation. All results are merged into a
single per-SNP summary file for additional downstream analysis and
results sharing. Therefore ***echolocatoR*** drastically reduces the
barriers to identifying causal variants by making the entire
fine-mapping pipeline rapid, robust and scalable.

![echoFlow](./images/echolocatoR_Fig1.png)

## Documentation

### [Website](https://rajlabmssm.github.io/echolocatoR)

### [Getting started](https://rajlabmssm.github.io/echolocatoR/articles/echolocatoR)

### Bugs/requests

Please report any bugs/requests on [GitHub
Issues](https://github.com/RajLabMSSM/echolocatoR/issues).

[Contributions](https://github.com/RajLabMSSM/echolocatoR/pulls) are
welcome!

## Literature

### For applications of ***echolocatoR*** in the literature, please see:

> 1.  E Navarro, E Udine, K de Paiva Lopes, M Parks, G Riboldi, BM
>     Schilder…T Raj (2020) Dysregulation of mitochondrial and
>     proteo-lysosomal genes in Parkinson’s disease myeloid cells.
>     Nature Genetics. <https://doi.org/10.1101/2020.07.20.212407>
> 2.  BM Schilder, T Raj (2021) Fine-Mapping of Parkinson’s Disease
>     Susceptibility Loci Identifies Putative Causal Variants. Human
>     Molecular Genetics, ddab294,
>     <https://doi.org/10.1093/hmg/ddab294>  
> 3.  K de Paiva Lopes, G JL Snijders, J Humphrey, A Allan, M Sneeboer,
>     E Navarro, BM Schilder…T Raj (2022) Genetic analysis of the human
>     microglial transcriptome across brain regions, aging and disease
>     pathologies. Nature Genetics,
>     <https://doi.org/10.1038/s41588-021-00976-y>

<br>

## Installation

### General tips

-   We generally recommend users upgrading to R>=4.0.0 before trying to
    install *echolocatoR.* While *echolocatoR* should technically be
    able to run in R>=3.6.0, some additional challenges with getting
    dependency versions not to conflict with one another.

### Quick installation

In R:

``` r
if(!require("remotes")) install.packages("remotes")

remotes::install_github("RajLabMSSM/echolocatoR")
```

### Robust installation (*conda*)

As with most softwares, installation is half the battle. The easiest way
to install all of ***echolocatoR***’s dependencies (which include R,
Python, and command line tools) and make sure they play well together is
to create a [*conda*](https://docs.conda.io/en/latest/) environment.

1.  If you haven’t done so already, install
    [*conda*](https://docs.conda.io/en/latest/).

2.  In command line, create the env from the *.yml* file (this file
    tells *conda* what to install):
    `conda env create -f https://github.com/RajLabMSSM/echolocatoR/raw/master/inst/conda/echoR.yml`

3.  Activate the new env:  
    `conda activate echoR`

4.  Install *echolocatoR* from command line so that it installs
    **within** the *conda* env:

5.  Open Rstudio from the command line interface (not by clicking the
    Rstudio icon). This helps to ensure Rstudio can find the paths to
    the packages in the conda env:  
    `open model_celltype_conservation.Rproj`

    Alternatively, the *conda* env also comes with
    [*radian*](https://github.com/randy3k/radian), which is a convenient
    R console that’s much more advanced than the default R console, but
    doesn’t require access to a GUI. This can be especially useful on
    computing clusters that don’t support RStudio or other IDEs.  
    `radian`

6.  Finally, to make extra sure ***echolocatoR*** uses the packages in
    this env (esp. if using from RStudio), you can then supply the env
    name to the `finemap_loci()` function (and many other
    ***echolocatoR*** functions) using `conda_env="echoR"`.

### Clone installation (*Rstudio*)

Lastly, if you’d like (or if for some reason none of the other
installation methods are working for you), you can alternatively clone
and then build ***echolocatoR***:

1.  Clone *echolocatoR:  
    `git clone https://github.com/RajLabMSSM/echolocatoR.git`*
2.  Open *echolocatoR.Rproj* within the echolocatoR folder.
3.  Then, within *Rstudio*, build ***echolocatoR*** by clicking the
    following drop down menu items: `Build --> Install and Restart` (or
    pressing the keys `CMD + SHIFT + B` on a Mac).

### 

<br>

### Dependencies

#### R

For a full list of required and suggested packages, see
[DESCRIPTION](https://github.com/RajLabMSSM/echolocatoR/blob/master/DESCRIPTION).

Additionally, there’s some optional R dependencies (e.g.
[XGR](https://github.com/hfang-bristol/XGR),
[Rgraphviz](https://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html))
that can be a bit tricky to install, so we’ve removed them as
requirements and instead provided a separate R function that helps users
to install them afterwards if needed:

``` r
library(echolocatoR)
extra_installs()
```

#### Python

For a full list of required python packages, see the *conda* env
[*echoR.yml*](https://github.com/RajLabMSSM/echolocatoR/blob/master/inst/conda/echoR.yml).
But here are some of the key ones.

    - python>=3.6.1  
    - pandas>=0.25.0   
    - pandas-plink  
    - pyarrow  
    - fastparquet  
    - scipy  
    - scikit-learn  
    - tqdm  
    - bitarray  
    - networkx  
    - rpy2  
    - requests  

#### Command line

##### [tabix](http://www.htslib.org/doc/tabix.html)

-   Rapid querying of summary stats files.
-   To use it, specify `query_by="tabix"` in `finemap_loci()`.
-   If you encounter difficulties using a *conda* distribution of tabix,
    we recommend you uninstall it from the env and instead install its
    parent package, [*htslib*](https://anaconda.org/bioconda/htslib) as
    this should be more up to date. *htslib* is now included in the
    echoR *conda* env by default.
-   Alternatively, you may install *htslib* to your machine globally via
    [*Brew*](https://formulae.brew.sh/formula/htslib) (for Mac users) or
    from [source](http://www.htslib.org/download).

##### [bcftools](http://samtools.github.io/bcftools/bcftools.html)

-   Used here for filtering populations in vcf files.
-   Can be installed via
    [*Brew*](https://formulae.brew.sh/formula/bcftools) (for Mac users)
    or [*conda*](https://anaconda.org/bioconda/bcftools).

##### [axel](https://github.com/axel-download-accelerator/axel)

-   Rapid multi-core downloading of large files (e.g. LD matrices from
    UK Biobank).

-   To use it, specify `download_method="axel"` in `finemap_loci()`.

-   **Update**: A conda version of *axel* has been kindly provided by
    [@jdblischak](https://github.com/RajLabMSSM/echolocatoR/pull/23), no
    longer requiring a separate installation.

-   However, if you want to use *axel* without the conda env, see this
    [tutorial](https://www.tecmint.com/axel-commandline-download-accelerator-for-linux/)
    for more info on installation. Here’s a quick overview:

    -   **Mac**: Install [brew](https://brew.sh/), then:
        `brew install axel`
    -   **CentOS/RHEL 7**: `yum install epel-release; yum install axel`
    -   **Fedora**: `yum install axel; dnf install axel`
    -   **Debian Jessie (e.g. Ubuntu, Linux Mint)**:
        `aptitude install axel`

<br>

## Fine-mapping tools

***echolocatoR*** will automatically check whether you have the
necessary columns to run each tool you selected in
`finemap_loci(finemap_methods=...)`. It will remove any tools that for
which there are missing necessary columns, and produces a message
letting you know which columns are missing. Note that some columns (e.g.
`MAF`,`N`,`t-stat`) can be automatically inferred if missing.  
For easy reference, we list the necessary columns here as well.  
See `?finemap_loci()` for descriptions of these columns.  
All methods require the columns: `SNP`,`CHR`,`POS`,`Effect`,`StdErr`

Additional required columns:

### [ABF](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)

#### `proportion_cases`,`MAF`

### [FINEMAP](http://www.christianbenner.com)

#### `A1`,`A2`,`MAF`,`N`

### [SuSiE](https://github.com/stephenslab/susieR)

#### `N`

### [PolyFun](https://github.com/omerwe/polyfun)

#### `A1`,`A2`,`P`,`N`

### [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)

#### `A1`,`A2`,`t-stat`

### [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)

#### `A1`,`A2`,`Freq`,`P`,`N`

### [coloc](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)

#### `N`,`MAF`

<br>

## Multi-finemap results files

The main output of ***echolocatoR*** are the multi-finemap files (for
example, `data("BST1")`). They are stored in the locus-specific
*Multi-finemap* subfolders.

### Column descriptions

-   **Standardized GWAS/QTL summary statistics**: e.g.
    `SNP`,`CHR`,`POS`,`Effect`,`StdErr`. See `?finemap_loci()` for
    descriptions of each.  
-   **leadSNP**: The designated proxy SNP per locus, which is the SNP
    with the smallest p-value by default.
-   **\<tool>.CS**: The 95% probability Credible Set (CS) to which a SNP
    belongs within a given fine-mapping tool’s results. If a SNP is not
    in any of the tool’s CS, it is assigned `NA` (or `0` for the
    purposes of plotting).  
-   **\<tool>.PP**: The posterior probability that a SNP is causal for a
    given GWAS/QTL trait.  
-   **Support**: The total number of fine-mapping tools that include the
    SNP in its CS.
-   **Consensus_SNP**: By default, defined as a SNP that is included in
    the CS of more than `N` fine-mapping tool(s), i.e. `Support>1`
    (default: `N=1`).  
-   **mean.PP**: The mean SNP-wise PP across all fine-mapping tools
    used.
-   **mean.CS**: If mean PP is greater than the 95% probability
    threshold (`mean.PP>0.95`) then `mean.CS` is 1, else 0. This tends
    to be a very stringent threshold as it requires a high degree of
    agreement between fine-mapping tools.

### Notes

-   Separate multi-finemap files are generated for each LD reference
    panel used, which is included in the file name (e.g.
    *UKB_LD.Multi-finemap.tsv.gz*).

-   Each fine-mapping tool defines its CS and PP slightly differently,
    so please refer to the associated original publications for the
    exact details of how these are calculated (links provided above).

<br>

## Datasets

For more detailed information about each dataset, use `?`:  
`R   library(echolocatoR)   ?NOTT_2019.interactome # example dataset`

### Epigenomic & genome-wide annotations

#### [Nott et al. (2019)](https://science.sciencemag.org/content/366/6469/1134.abstract)

-   Data from this publication contains results from cell type-specific
    (neurons, oligodendrocytes, astrocytes, microglia, & peripheral
    myeloid cells) epigenomic assays (H3K27ac, ATAC, H3K4me3) from human
    brain tissue.

-   For detailed metadata, see:

    ``` r
    data("NOTT_2019.bigwig_metadata")
    ```

-   Built-in datasets:

    -   Enhancer/promoter coordinates (as *GenomicRanges*)

    ``` r
    data("NOTT_2019.interactome")
    # Examples of the data nested in "NOTT_2019.interactome" object:
    NOTT_2019.interactome$`Neuronal promoters`
    NOTT_2019.interactome$`Neuronal enhancers`
    NOTT_2019.interactome$`Microglia promoters`
    NOTT_2019.interactome$`Microglia enhancers`
    ...
    ...
    ```

    -   PLAC-seq enhancer-promoter interactome coordinates

    ``` r
    NOTT_2019.interactome$H3K4me3_around_TSS_annotated_pe
    NOTT_2019.interactome$`Microglia interactome`
    NOTT_2019.interactome$`Neuronal interactome`
    NOTT_2019.interactome$`Oligo interactome`
    ...
    ...
    ```

-   API access to full bigWig files on UCSC Genome Browser, which
    includes

    -   Epigenomic reads (as *GenomicRanges*)  
    -   Aggregate epigenomic *score* for each cell type - assay
        combination

#### [Corces et al. (2020)](https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1)

-   Data from this preprint contains results from bulk and single-cell
    chromatin accessibility epigenomic assays in 39 human brains.

    ``` r
    data("CORCES_2020.bulkATACseq_peaks")
    data("CORCES_2020.cicero_coaccessibility")
    data("CORCES_2020.HiChIP_FitHiChIP_loop_calls")
    data("CORCES_2020.scATACseq_celltype_peaks")
    data("CORCES_2020.scATACseq_peaks")
    ```

#### [XGR](http://xgr.r-forge.r-project.org)

-   API access to a diverse library of cell type/line-specific
    epigenomic (e.g. ENCODE) and other genome-wide annotations.

#### [Roadmap](http://www.roadmapepigenomics.org)

-   API access to cell type-specific epigenomic data.

#### [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

-   API access to various genome-wide SNP annotations (e.g. missense,
    nonsynonmous, intronic, enhancer).

#### [HaploR](https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html)

-   API access to known per-SNP QTL and epigenomic data hits.

### QTLs

#### [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)

-   API access to full summary statistics from many standardized
    e/s/t-QTL datasets.  
-   Data access and colocalization tests facilitated through the
    [catalogueR](https://github.com/RajLabMSSM/catalogueR) R package.

<br>

## Enrichment tools

### [XGR](http://xgr.r-forge.r-project.org)

-   Binomial enrichment tests between customisable foreground and
    background SNPs.

### [GoShifter](https://github.com/immunogenomics/goshifter)

-   LD-informed iterative enrichment analysis.

### [S-LDSC](https://www.nature.com/articles/ng.3954)

-   Genome-wide stratified LD score regression.
-   Inlccles 187-annotation baseline model from [Gazal et al.
    2018](https://www.nature.com/articles/s41588-018-0231-8).  
-   You can alternatively supply a custom annotations matrix.

### [motifbreakR](https://github.com/Simon-Coetzee/motifBreakR)

-   Identification of transcript factor binding motifs (TFBM) and
    prediction of SNP disruption to said motifs.
-   Includes a comprehensive list of TFBM databases via
    [MotifDB](https://bioconductor.org/packages/release/bioc/html/MotifDb.html)
    (9,900+ annotated position frequency matrices from 14 public
    sources, for multiple organisms).

### [GARFIELD](https://www.bioconductor.org/packages/release/bioc/html/garfield.html) (**under construction**)

-   Genomic enrichment with LD-informed heuristics.

<br>

## LD reference panels

### [UK Biobank](https://www.ukbiobank.ac.uk)

### [1000 Genomes Phase 1](https://www.internationalgenome.org)

### [1000 Genomes Phase 3](https://www.internationalgenome.org)

<hr>

## Creator

<a href="https://bschilder.github.io/BMSchilder/" target="_blank">Brian
M. Schilder, Bioinformatician II</a>  
<a href="https://rajlab.org" target="_blank">Raj Lab</a>  
<a href="https://icahn.mssm.edu/about/departments/neuroscience" target="_blank">Department
of Neuroscience, Icahn School of Medicine at Mount Sinai</a>  
![Sinai](./images/sinai.png)
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

## 1. Bug description

(A clear and concise description of what the bug is.)

### Console output 

```
# Paste console output here (e.g. from R/python/command line)

```

### Expected behaviour   

(A clear and concise description of what you expected to happen.)


## 2. Reproducible example   

### Code 

(Please add the steps to reproduce the bug here. See [here](https://www.r-bloggers.com/2020/10/how-to-make-a-reprex/) for an intro to making a reproducible example (i.e. reprex) and why they're important! __This will help us to help you much faster.__)

```R
# Paste example here

```

### Data

(If possible, upload a small sample of your data so that we can reproduce the bug on our end. If that's not possible, please at least include a screenshot of your data and other relevant details.)


## 3. Session info

(Add output of the R function `utils::sessionInfo()` below. This helps us assess version/OS conflicts which could be causing bugs.)

<details>

```
# Paste utils::sessionInfo() output 

```
</details>
---
title: ""  
author: "<img src='https://github.com/RajLabMSSM/`r read.dcf('DESCRIPTION', fields = 'Package')[1]`/raw/`r gsub('[*] ','',grep('\\*',system('git branch', intern = TRUE), value = TRUE))[1]`/inst/hex/hex.png' height='300'><br><br>
        `r badger::badge_github_version(color = 'black')` 
        `r badger::badge_github_actions(action = 'R-CMD-check-bioc')`
        `r badger::badge_last_commit()`
        `r badger::badge_codecov()` 
        `r badger::badge_license()`
        `r badger::badge_doi('10.1093/bioinformatics/btab658', 'blue')`
        <h5>Author: <i>Brian M. Schilder</i></h5>" 
date: "<h5>README updated: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h5>"
output:
  github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
pkg <- read.dcf("DESCRIPTION", fields = "Package")[1]
description <- read.dcf("DESCRIPTION", fields = "Description")[1]
```
 

<center><h1>  ) ) ) ) ))) :bat: echolocatoR :bat: ((( ( ( ( ( </h1></center>

<h3> Automated statistical and functional fine-mapping with extensive access to genome-wide datasets </h3>

<hr>

If you use ``r pkg``, please cite: 

> `r citation(pkg)$textVersion`

## Introduction

Fine-mapping methods are a powerful means of identifying causal variants
underlying a given phenotype, but are underutilized due to the technical
challenges of implementation. ***echolocatoR*** is an R package that
automates end-to-end genomics fine-mapping, annotation, and plotting in
order to identify the most probable causal variants associated with a
given phenotype.

It requires minimal input from users (a GWAS or QTL summary statistics
file), and includes a suite of statistical and functional fine-mapping
tools. It also includes extensive access to datasets (linkage
disequilibrium panels, epigenomic and genome-wide annotations, QTL).

The elimination of data gathering and preprocessing steps enables rapid
fine-mapping of many loci in any phenotype, complete with locus-specific
publication-ready figure generation. All results are merged into a
single per-SNP summary file for additional downstream analysis and
results sharing. Therefore ***echolocatoR*** drastically reduces the
barriers to identifying causal variants by making the entire
fine-mapping pipeline rapid, robust and scalable.

![echoFlow](./images/echolocatoR_Fig1.png)

 
## Documentation 

### [Website](https://rajlabmssm.github.io/`r pkg`) 
### [Getting started](https://rajlabmssm.github.io/`r pkg`/articles/`r pkg`) 

### Bugs/requests

Please report any bugs/requests on 
[GitHub Issues](https://github.com/RajLabMSSM/echolocatoR/issues).

[Contributions](https://github.com/RajLabMSSM/echolocatoR/pulls) are welcome!

## Literature
 
### For applications of ***echolocatoR*** in the literature, please see:

> 1.	E Navarro, E Udine, K de Paiva Lopes, M Parks, G Riboldi, BM Schilder…T Raj (2020) Dysregulation of mitochondrial and proteo-lysosomal genes in Parkinson's disease myeloid cells. Nature Genetics. https://doi.org/10.1101/2020.07.20.212407 
> 2.	BM Schilder, T Raj (2021) Fine-Mapping of Parkinson’s Disease Susceptibility Loci Identifies Putative Causal Variants. Human Molecular Genetics, ddab294, https://doi.org/10.1093/hmg/ddab294   
> 3. K de Paiva Lopes, G JL Snijders, J Humphrey, A Allan, M Sneeboer, E Navarro, BM Schilder…T Raj (2022) Genetic analysis of the human microglial transcriptome across brain regions, aging and disease pathologies. Nature Genetics, https://doi.org/10.1038/s41588-021-00976-y  


<br>

## Installation

### General tips

-   We generally recommend users upgrading to R\>=4.0.0 before trying to
    install *echolocatoR.* While *echolocatoR* should technically be
    able to run in R\>=3.6.0, some additional challenges with getting
    dependency versions not to conflict with one another.

### Quick installation

In R:

```R
if(!require("remotes")) install.packages("remotes")

remotes::install_github("RajLabMSSM/echolocatoR")
```

### Robust installation (*conda*)

As with most softwares, installation is half the battle. The easiest way
to install all of ***echolocatoR***'s dependencies (which include R,
Python, and command line tools) and make sure they play well together is
to create a [*conda*](https://docs.conda.io/en/latest/) environment.

1.  If you haven't done so already, install
    [*conda*](https://docs.conda.io/en/latest/). 
    
2.  In command line, create the env from the *.yml* file (this file tells *conda* what to install): 
    `conda env create -f https://github.com/RajLabMSSM/echolocatoR/raw/master/inst/conda/echoR.yml`

3.  Activate the new env:  
    `conda activate echoR`

4.  Install *echolocatoR* from command line so that it installs
    **within** the *conda* env:

5.  Open Rstudio from the command line interface (not by clicking the
    Rstudio icon). This helps to ensure Rstudio can find the paths to
    the packages in the conda env:  
    `open model_celltype_conservation.Rproj`

    Alternatively, the *conda* env also comes with
    [*radian*](https://github.com/randy3k/radian), which is a convenient
    R console that's much more advanced than the default R console, but
    doesn't require access to a GUI. This can be especially useful on
    computing clusters that don't support RStudio or other IDEs.  
    `radian`

6.  Finally, to make extra sure ***echolocatoR*** uses the packages in
    this env (esp. if using from RStudio), you can then supply the env
    name to the `finemap_loci()` function (and many other ***echolocatoR***
    functions) using `conda_env="echoR"`. 
    
### Clone installation (*Rstudio*)

Lastly, if you'd like (or if for some reason none of the other
installation methods are working for you), you can alternatively clone
and then build ***echolocatoR***:

1.  Clone *echolocatoR:  
    `git clone https://github.com/RajLabMSSM/echolocatoR.git`*
2.  Open *echolocatoR.Rproj* within the echolocatoR folder.
3.  Then, within *Rstudio*, build ***echolocatoR*** by clicking the
    following drop down menu items: `Build --> Install and Restart` (or
    pressing the keys `CMD + SHIFT + B` on a Mac).

### 

<br>

### Dependencies

#### R

For a full list of required and suggested packages, see
[DESCRIPTION](https://github.com/RajLabMSSM/echolocatoR/blob/master/DESCRIPTION).

Additionally, there's some optional R dependencies (e.g.
[XGR](https://github.com/hfang-bristol/XGR),
[Rgraphviz](https://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html))
that can be a bit tricky to install, so we've removed them as
requirements and instead provided a separate R function that helps users
to install them afterwards if needed:

``` r
library(echolocatoR)
extra_installs()
```

#### Python

For a full list of required python packages, see the *conda* env
[*echoR.yml*](https://github.com/RajLabMSSM/echolocatoR/blob/master/inst/conda/echoR.yml).
But here are some of the key ones.

    - python>=3.6.1  
    - pandas>=0.25.0   
    - pandas-plink  
    - pyarrow  
    - fastparquet  
    - scipy  
    - scikit-learn  
    - tqdm  
    - bitarray  
    - networkx  
    - rpy2  
    - requests  

#### Command line

##### [tabix](http://www.htslib.org/doc/tabix.html)

-   Rapid querying of summary stats files.
-   To use it, specify `query_by="tabix"` in `finemap_loci()`.
-   If you encounter difficulties using a *conda* distribution of tabix,
    we recommend you uninstall it from the env and instead install its
    parent package, [*htslib*](https://anaconda.org/bioconda/htslib) as
    this should be more up to date. *htslib* is now included in the
    echoR *conda* env by default.
-   Alternatively, you may install *htslib* to your machine globally via
    [*Brew*](https://formulae.brew.sh/formula/htslib) (for Mac users) or
    from [source](http://www.htslib.org/download).

##### [bcftools](http://samtools.github.io/bcftools/bcftools.html)

-   Used here for filtering populations in vcf files.
-   Can be installed via
    [*Brew*](https://formulae.brew.sh/formula/bcftools) (for Mac users)
    or [*conda*](https://anaconda.org/bioconda/bcftools).

##### [axel](https://github.com/axel-download-accelerator/axel)

-   Rapid multi-core downloading of large files (e.g. LD matrices from
    UK Biobank).

-   To use it, specify `download_method="axel"` in `finemap_loci()`.

-   **Update**: A conda version of *axel* has been kindly provided by
    [\@jdblischak](https://github.com/RajLabMSSM/echolocatoR/pull/23),
    no longer requiring a separate installation.

-   However, if you want to use *axel* without the conda env, see this
    [tutorial](https://www.tecmint.com/axel-commandline-download-accelerator-for-linux/)
    for more info on installation. Here's a quick overview:

    -   **Mac**: Install [brew](https://brew.sh/), then:
        `brew install axel`
    -   **CentOS/RHEL 7**: `yum install epel-release; yum install axel`
    -   **Fedora**: `yum install axel; dnf install axel`
    -   **Debian Jessie (e.g. Ubuntu, Linux Mint)**:
        `aptitude install axel`

<br>

## Fine-mapping tools

***echolocatoR*** will automatically check whether you have the
necessary columns to run each tool you selected in
`finemap_loci(finemap_methods=...)`. It will remove any tools that for
which there are missing necessary columns, and produces a message
letting you know which columns are missing. Note that some columns (e.g.
`MAF`,`N`,`t-stat`) can be automatically inferred if missing.  
For easy reference, we list the necessary columns here as well.  
See `?finemap_loci()` for descriptions of these columns.  
All methods require the columns: `SNP`,`CHR`,`POS`,`Effect`,`StdErr`

Additional required columns:

### [ABF](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)

#### `proportion_cases`,`MAF`

### [FINEMAP](http://www.christianbenner.com)

#### `A1`,`A2`,`MAF`,`N`

### [SuSiE](https://github.com/stephenslab/susieR)

#### `N`

### [PolyFun](https://github.com/omerwe/polyfun)

#### `A1`,`A2`,`P`,`N`

### [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)

#### `A1`,`A2`,`t-stat`

### [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)

#### `A1`,`A2`,`Freq`,`P`,`N`

### [coloc](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)

#### `N`,`MAF`

<br>

## Multi-finemap results files

The main output of ***echolocatoR*** are the multi-finemap files (for
example, `data("BST1")`). They are stored in the locus-specific
*Multi-finemap* subfolders.

### Column descriptions

-   **Standardized GWAS/QTL summary statistics**: e.g.
    `SNP`,`CHR`,`POS`,`Effect`,`StdErr`. See `?finemap_loci()` for
    descriptions of each.  
-   **leadSNP**: The designated proxy SNP per locus, which is the SNP
    with the smallest p-value by default.
-   **\<tool\>.CS**: The 95% probability Credible Set (CS) to which a
    SNP belongs within a given fine-mapping tool's results. If a SNP is
    not in any of the tool's CS, it is assigned `NA` (or `0` for the
    purposes of plotting).  
-   **\<tool\>.PP**: The posterior probability that a SNP is causal for
    a given GWAS/QTL trait.  
-   **Support**: The total number of fine-mapping tools that include the
    SNP in its CS.
-   **Consensus_SNP**: By default, defined as a SNP that is included in
    the CS of more than `N` fine-mapping tool(s), i.e. `Support>1`
    (default: `N=1`).  
-   **mean.PP**: The mean SNP-wise PP across all fine-mapping tools
    used.
-   **mean.CS**: If mean PP is greater than the 95% probability
    threshold (`mean.PP>0.95`) then `mean.CS` is 1, else 0. This tends
    to be a very stringent threshold as it requires a high degree of
    agreement between fine-mapping tools.

### Notes

-   Separate multi-finemap files are generated for each LD reference
    panel used, which is included in the file name (e.g.
    *UKB_LD.Multi-finemap.tsv.gz*).

-   Each fine-mapping tool defines its CS and PP slightly differently,
    so please refer to the associated original publications for the
    exact details of how these are calculated (links provided above).

<br>

## Datasets

For more detailed information about each dataset, use `?`:  
`R   library(echolocatoR)   ?NOTT_2019.interactome # example dataset`

### Epigenomic & genome-wide annotations

#### [Nott et al. (2019)](https://science.sciencemag.org/content/366/6469/1134.abstract)

-   Data from this publication contains results from cell type-specific
    (neurons, oligodendrocytes, astrocytes, microglia, & peripheral
    myeloid cells) epigenomic assays (H3K27ac, ATAC, H3K4me3) from human
    brain tissue.

-   For detailed metadata, see:

    ``` r
    data("NOTT_2019.bigwig_metadata")
    ```

-   Built-in datasets:

    -   Enhancer/promoter coordinates (as *GenomicRanges*)

    ``` r
    data("NOTT_2019.interactome")
    # Examples of the data nested in "NOTT_2019.interactome" object:
    NOTT_2019.interactome$`Neuronal promoters`
    NOTT_2019.interactome$`Neuronal enhancers`
    NOTT_2019.interactome$`Microglia promoters`
    NOTT_2019.interactome$`Microglia enhancers`
    ...
    ...
    ```

    -   PLAC-seq enhancer-promoter interactome coordinates

    ``` r
    NOTT_2019.interactome$H3K4me3_around_TSS_annotated_pe
    NOTT_2019.interactome$`Microglia interactome`
    NOTT_2019.interactome$`Neuronal interactome`
    NOTT_2019.interactome$`Oligo interactome`
    ...
    ...
    ```

-   API access to full bigWig files on UCSC Genome Browser, which
    includes

    -   Epigenomic reads (as *GenomicRanges*)  
    -   Aggregate epigenomic *score* for each cell type - assay
        combination

#### [Corces et al. (2020)](https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1)

-   Data from this preprint contains results from bulk and single-cell
    chromatin accessibility epigenomic assays in 39 human brains.

    ``` r
    data("CORCES_2020.bulkATACseq_peaks")
    data("CORCES_2020.cicero_coaccessibility")
    data("CORCES_2020.HiChIP_FitHiChIP_loop_calls")
    data("CORCES_2020.scATACseq_celltype_peaks")
    data("CORCES_2020.scATACseq_peaks")
    ```

#### [XGR](http://xgr.r-forge.r-project.org)

-   API access to a diverse library of cell type/line-specific
    epigenomic (e.g. ENCODE) and other genome-wide annotations.

#### [Roadmap](http://www.roadmapepigenomics.org)

-   API access to cell type-specific epigenomic data.

#### [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

-   API access to various genome-wide SNP annotations (e.g. missense,
    nonsynonmous, intronic, enhancer).

#### [HaploR](https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html)

-   API access to known per-SNP QTL and epigenomic data hits.

### QTLs

#### [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)

-   API access to full summary statistics from many standardized
    e/s/t-QTL datasets.  
-   Data access and colocalization tests facilitated through the
    [catalogueR](https://github.com/RajLabMSSM/catalogueR) R package.

<br>

## Enrichment tools

### [XGR](http://xgr.r-forge.r-project.org)

-   Binomial enrichment tests between customisable foreground and
    background SNPs.

### [GoShifter](https://github.com/immunogenomics/goshifter)

-   LD-informed iterative enrichment analysis.

### [S-LDSC](https://www.nature.com/articles/ng.3954)

-   Genome-wide stratified LD score regression.
-   Inlccles 187-annotation baseline model from [Gazal et al.
    2018](https://www.nature.com/articles/s41588-018-0231-8).  
-   You can alternatively supply a custom annotations matrix.

### [motifbreakR](https://github.com/Simon-Coetzee/motifBreakR)

-   Identification of transcript factor binding motifs (TFBM) and
    prediction of SNP disruption to said motifs.
-   Includes a comprehensive list of TFBM databases via
    [MotifDB](https://bioconductor.org/packages/release/bioc/html/MotifDb.html)
    (9,900+ annotated position frequency matrices from 14 public
    sources, for multiple organisms).

### [GARFIELD](https://www.bioconductor.org/packages/release/bioc/html/garfield.html) (**under construction**)

-   Genomic enrichment with LD-informed heuristics.

<br>

## LD reference panels

### [UK Biobank](https://www.ukbiobank.ac.uk)

### [1000 Genomes Phase 1](https://www.internationalgenome.org)

### [1000 Genomes Phase 3](https://www.internationalgenome.org)


<hr>

## Creator

<a href="https://bschilder.github.io/BMSchilder/" target="_blank">Brian
M. Schilder, Bioinformatician II</a>  
<a href="https://rajlab.org" target="_blank">Raj Lab</a>  
<a href="https://icahn.mssm.edu/about/departments/neuroscience" target="_blank">Department
of Neuroscience, Icahn School of Medicine at Mount Sinai</a>  
![Sinai](./images/sinai.png)
---
title: "PD loci" 
author: "<h4>Author: <i>Brian M. Schilder</i></h4>" 
date: "<h4>Updated: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h4>"
output:
  BiocStyle::html_document
vignette: >
    %\VignetteIndexEntry{PD_loci} 
    %\usepackage[utf8]{inputenc}
    %\VignetteEngine{knitr::rmarkdown}
---

```{r, include = T, error = F, dpi=300}
root.dir <- tempdir()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = root.dir, 
  fig.height = 12,
  fig.width = 10
)  
knitr::opts_knit$set(
  root.dir = root.dir
  )
# knitr::opts_chunk$get("root.dir")
# devtools::build_vignettes(quiet = F)
```

```{r setup, eval=F}
library(echolocatoR) 
library(data.table)

results_dir <- file.path(root.dir,"results")
```


Here we take advantage of the fine-mapping results files already available on the [echolocatoR Fine-mapping Portal](https://rajlab.shinyapps.io/Fine_Mapping_Shiny).

# Download data  


*echolocatoR* includes some functions to search for and extract files from GitHub repos. In this case we'll get files from the [echolocatoR Fine-mapping Portal](https://rajlab.shinyapps.io/Fine_Mapping_Shiny), but these functions can be used on any public GitHub repo.  

* `GITHUB.list_files()`: Searches the repo for files matching your regex query.
* `GITHUB.download_files()`: Not only does this download your files, but it reconstructs the original folder structure they were found in. This is useful for *echolocatoR* output files, which are automatically organized in a hierarchical folder structure (e.g. *results/GWAS/Nalls23andMe_2019*).


```{r locus_list}
locus_list <- c("LRRK2","MBNL2","DYRK1A","FCGR2A","MED12L")  
```

## Multi-finemap data  

First, we'll download the multi-finemap results files.  

```{r multi_finemap, eval=F}
local_finemap <- GITHUB.portal_query(results_dir = results_dir, 
                                     dataset_types="GWAS",
                                     phenotypes = c("parkinson"),
                                     file_types = "multi_finemap",
                                     loci = locus_list,
                                     LD_panels=c("UKB"))  
```

## LD data

Next, we'll download the files containing LD with the lead GWAS SNP. 

```{r LD, eval=F}
local_ld <- GITHUB.portal_query(results_dir = results_dir, 
                                 dataset_types="GWAS",
                                 phenotypes = c("parkinson"),
                                 file_types = "LD",
                                 loci = locus_list,
                                 LD_panels=c("UKB"))  
```

## Data dictionary

Now let's put them in a data dictionary for easy retrieval later.
This way, you can simply find the right file path my locus name.

This function also infers a number of other useful variables from this data. 

```{r Data dictionary, eval=F}
named_list <- list("finemap"=local_finemap,
                   "LD"=local_ld)
data_dict <- GITHUB.make_data_dict(named_list)
```


```{r LRRK2, attr.output='style="max-height: 200px;"', eval=F} 
locus <- "LRRK2"

track_order <- names(PLOT.heights_dict())
track_order <- R.utils::insert(track_order[track_order!="Genes"], 1, "Genes")

LRRK2 <- PLOT.locus(finemap_dat=fread(data_dict$finemap[[locus]]),   
                    LD_matrix=fread(data_dict$LD[[locus]]),
                    LD_reference="UKB",
                    locus_dir=data_dict$locus_dir[[locus]], 
                    
                    max_transcripts = 3,
                    Nott_binwidth = 200,
                    Nott_epigenome=T,   
                    Nott_regulatory_rects = T, 
                    Nott_show_placseq = T,
                    
                    save_plot=F, 
                    return_list=T, 
                    zoom_exceptions_str = "*full window$|zoom_polygon|Genes",
                    track_order = track_order,
                    plot.zoom=c("1x","20x","25x")) 
```

```{r MBNL2, attr.output='style="max-height: 200px;"', eval=F} 
locus <- "MBNL2"

track_order <- names(PLOT.heights_dict())
track_order <- R.utils::insert(track_order[track_order!="Genes"], 1, "Genes")


MBNL2 <- PLOT.locus(finemap_dat=fread(data_dict$finemap[[locus]]),   
                    LD_matrix=fread(data_dict$LD[[locus]]),
                    LD_reference="UKB",
                    locus_dir=data_dict$locus_dir[[locus]], 
                    
                    max_transcripts = 1,
                    Nott_epigenome=T,  
                    Nott_binwidth = 200,
                    Nott_regulatory_rects = T, 
                    Nott_show_placseq = T,
                    
                    save_plot=F, 
                    return_list=T,
                    zoom_exceptions_str = "*full window$|zoom_polygon|Genes",
                    track_order = track_order,
                    plot.zoom=c("1x","10x","20x","25x"))
```
 
 
```{r DYRK1A, attr.output='style="max-height: 200px;"', eval=F} 
locus <- "DYRK1A"

DYRK1A <- PLOT.locus(finemap_dat=fread(data_dict$finemap[[locus]]),   
                    LD_matrix=fread(data_dict$LD[[locus]]),
                    LD_reference="UKB",
                    locus_dir=data_dict$locus_dir[[locus]], 
                    
                    max_transcripts = 3,
                    Nott_epigenome=T,  
                    Nott_binwidth = 500,
                    Nott_regulatory_rects = T, 
                    Nott_show_placseq = T,
                    
                    save_plot=F, 
                    return_list=T,
                    plot.zoom=c("25x"))
```
 
This region has a lot of genes so we're going to subset the main view a bit. 
 
```{r FCGR2A, attr.output='style="max-height: 200px;"', eval=F} 
locus <- "FCGR2A"
# 
# finemap_dat <- fread(data_dict$finemap[[locus]])
# xlims <- PLOT.get_window_limits(finemap_dat = finemap_dat,  
#                                 plot.zoom = "3x")
# finemap_dat <- subset(finemap_dat, Mb>xlims[1] & Mb<xlims[2] )

FCGR2A <- PLOT.locus(finemap_dat=fread(data_dict$finemap[[locus]]),   
                    LD_matrix=fread(data_dict$LD[[locus]]),
                    LD_reference="UKB",
                    locus_dir=data_dict$locus_dir[[locus]], 
                    
                    max_transcripts = 3,
                    Nott_epigenome=T,  
                    # Nott_binwidth = 500,
                    Nott_regulatory_rects = T, 
                    Nott_show_placseq = T,
                    
                    save_plot=F, 
                    return_list=T,
                    plot.zoom=c("25x")) 
``` 

# Modify plots

```{r, eval=F}
zoom <- "25x"
locus <- "FCGR2A"#"FCGR2A"
plot_zoom_list <- FCGR2A

mod_plot <- function(plot_zoom_list, locus, zoom){
  # deparse(substitute(FCGR2A))
  TRKS <- plot_zoom_list[[zoom]]
  finemap_dat <- data.table::fread(data_dict$finemap[[locus]])
  
  n_snps <- if(data_dict$dataset_type %in% names(TRKS)){
          paste0("n SNPs: ",
                 nrow(ggplot2::ggplot_build(TRKS[[data_dict$dataset_type]])$data[[2]]),", ")
        } else {NULL}
  title_text <- paste0(basename(data_dict$locus_dir[[locus]]),"   (",n_snps,"zoom: ",zoom,")")
     
  # Alter track
  TRKS[["Genes"]] <- TRKS[["Genes"]] +  
    scale_y_discrete(expand = expansion(add  = c(.0,-1)))
  
  
  # Fuse tracks
  LOCUS_zoom <- patchwork::wrap_plots(TRKS, 
                                       ncol = 1, 
                                       heights = heights) +
      patchwork::plot_annotation(title = title_text)
   
  # Save plot
  ggsave(filename = file.path(results_dir,"GWAS/Nalls23andMe_2019",locus,
                              paste0("multiview.",locus,
                                     ".UKB.",zoom,"_mod.jpg")), 
         plot = LOCUS_zoom, 
         dpi = 300, 
         height=12,
         width=10)
  return(LOCUS_zoom)
}


```

# Session info

<details> 

```{r Session Info, attr.output='style="max-height: 200px;"'}
utils::sessionInfo()
```

</details>

---
title: "Summarise vignette" 
author: "<h4>Author: <i>Brian M. Schilder</i></h4>" 
date: "<h4>Updated: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h4>"
output:
  BiocStyle::html_document
vignette: >
    %\VignetteIndexEntry{summarise_vignette} 
    %\usepackage[utf8]{inputenc}
    %\VignetteEngine{knitr::rmarkdown}
---

```{r, include = T, message = F, error = T}
root.dir <- tempdir()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = root.dir, 
  fig.height = 12,
  fig.width = 10
)  
knitr::opts_knit$set(
  root.dir = root.dir
  )
library(echolocatoR) 
results_dir <- file.path(root.dir,"results")
```


```{r setup}
library(echolocatoR)
# if(!require("XGR")) install.packages("XGR")  
```


# Download data  

* Many fine-mapping results can be found on the [echolocatoR Fine-mapping Portal](https://rajlab.shinyapps.io/Fine_Mapping_Shiny).  
* We will download some of those files here using some built-in APIs.  

```{r Download data, eval=F} 
local_files <- GITHUB.portal_query(results_dir = results_dir, 
                                   dataset_types="GWAS",
                                   phenotypes = c("parkinson"),
                                   file_types = "multi_finemap",
                                   # loci = c("BST1","CHRNB1","LRRK2"),
                                   LD_panels=c("UKB")) 
```

# Merge `echolocatoR` results  

Gather all of the fine-mapping results generated by `finemap_loci()` previously.  

```{r merge_finemapping_results(), eval=F}
merged_DT <- merge_finemapping_results(minimum_support=0,
                                       include_leadSNPs=T,
                                       dataset = file.path(results_dir,
                                                           "GWAS/Nalls23andMe_2019"),  
                                       from_storage=T,
                                       consensus_threshold = 2,
                                       top_CS_only = T,
                                       haploreg_annotation=F,
                                       biomart_annotation=F, 
                                       verbose = F)    
results_report(merged_DT)
```

# Summarise

## `SUMMARISE.get_SNPgroup_counts()`

Get the number of SNPs for each SNP group per locus.  
It also prints the mean number of SNPs for each SNP group across all loci.  
**NOTE**: You will need to make sure to set `merge_finemapping_results(minimum_support=1)`
in the above step to get accurate counts for all SNP groups.  

```{r SUMMARISE.get_SNPgroup_counts(), eval=F}
snp_groups <- SUMMARISE.get_SNPgroup_counts(merged_DT=merged_DT)   
createDT(snp_groups)
```

## `SUMMARISE.get_CS_counts()` 

County the number of tool-specific and UCS Credible Set SNPs per locus.  

```{r SUMMARISE.get_CS_counts(), eval=F}
UCS_counts <- SUMMARISE.get_CS_counts(merged_DT = merged_DT)
print(UCS_counts)
```


# Plot  

- The following functions each return a list containing both the `...$plot` and the `...$data`
used to make the plot.  
- Where available, `snp_filter` allows user to use any filtering argument (supplied as a string)
to subset the data they want to use in the plot/data.  

## `SUMMARISE.CS_bin_plot()`  

Plot the number of number that had Credible Sets of a certain size (per tool).  

```{r SUMMARISE.CS_bin_plot(), eval=F}
bin_plot <- SUMMARISE.CS_bin_plot(merged_DT = merged_DT, 
                                  show_plot = T) 
createDT(bin_plot$data)
```

## `SUMMARISE.coloc_nominated_eGenes()`  

If you ran colocalization tests with `echolocatoR` (via `catalogueR`) you can use those results to 
come up with a top QTL nominated gene for each locus (potentially implicating that gene in your phenotype).  

### Download COLOC results from [echolocatoR Fine-mapping Portal](https://rajlab.shinyapps.io/Fine_Mapping_Shiny).    

```{r, eval=F}
remote_coloc <- "https://github.com/RajLabMSSM/Fine_Mapping_Shiny/raw/master/www/data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz"
local_coloc <- gsub("https://github.com/RajLabMSSM/Fine_Mapping_Shiny/raw/master/www/data",
                    results_dir,
                    remote_coloc)
if(!file.exists(local_coloc)){
  dir.create(dirname(local_coloc), showWarnings = F, recursive = T)
  download.file(remote_coloc, destfile = local_coloc)
} 
```


```{r SUMMARISE.coloc_nominated_eGenes(), eval=F} 

gg_egene <- SUMMARISE.coloc_nominated_eGenes(coloc_results = local_coloc,
                                             merged_DT = merged_DT,
                                             PP_threshold = .8,
                                             fill_var = NULL,
                                             text_size = 2.5,
                                             y_lab = "Locus",
                                             x_lab = NULL,
                                             label_yaxis = T, 
                                             show_plot = T)
createDT(gg_egene$data)
```

## `SUMMARISE.CS_counts_plot()`  

An extension of `SUMMARISE.get_CS_counts()` that plots these results.  

```{r SUMMARISE.CS_counts_plot(), eval=F}
gg_CS <- SUMMARISE.CS_counts_plot(merged_DT = merged_DT, 
                                  show_numbers=T,
                                  label_yaxis=T, 
                                  ylabel=NULL, 
                                  show_plot = T)
createDT(gg_CS$data)
```


## `ANNOTATE.plot_missense()`  

Query *biomart* for SNP-level annotations, and then return only those that contain missense mutations.  

```{r, eval=F}
# Run this in case you don't want to wait to ANNOTATE.plot_missense() to finish querying.
gg_missense <- list(plot=patchwork::plot_spacer())
```

```{r ANNOTATE.plot_missense(), eval=F}
gg_missense <- ANNOTATE.plot_missense(merged_DT = merged_DT, 
                                      snp_filter="Support>0", 
                                      label_yaxis = T,
                                      show_plot = T)  
createDT(gg_missense$data)
```


## `SUMMARISE.peak_overlap_plot()`  

Plot the overlap between some SNP group and cell-type-specific
epigenomic peaks / interactome anchors.  

```{r SUMMARISE.peak_overlap_plot(), eval=F}
gg_peaks <- SUMMARISE.peak_overlap_plot(merged_DT, 
                                        snp_filter="Consensus_SNP==T",
                                        include.NOTT_2019_peaks=T,
                                        include.NOTT_2019_enhancers_promoters=T,
                                        include.NOTT_2019_PLACseq=T,
                                        include.CORCES_2020_scATACpeaks=T,
                                        include.CORCES_2020_Cicero_coaccess=T,
                                        include.CORCES_2020_bulkATACpeaks=T,
                                        include.CORCES_2020_HiChIP_FitHiChIP_coaccess=T,
                                        include.CORCES_2020_gene_annotations=T,
                                        plot_celltype_specificity=F,
                                        facets_formula=". ~ Cell_type",
                                        show_plot=T,  
                                        label_yaxis=T, 
                                        x_strip_angle = 90,
                                        drop_empty_cols = T,
                                        fill_title="Consensus SNPs\nin epigenomic peaks",
                                        save_path=file.path(results_dir,"GWAS/Nalls23andMe_2019/_genome_wide",
                                                            "peak_plot.consensus.png"),
                                        verbose=T)  
createDT(gg_peaks$data)
```

## Multi-plot  

Now for the grand finale, put them all together by calling a single function,`super_summary_plot()`!  
 

```{r super_summary_plot(), eval=F}
super_plot <- super_summary_plot(merged_DT = merged_DT,
                                  snp_filter="Consensus_SNP==T",
                                  coloc_results=local_coloc,
                                  plot_missense=F,
                                  show_plot=T, 
                                  
                                  save_plot=file.path(results_dir,
                                                      "GWAS/Nalls23andMe_2019/_genome_wide",
                                                      "PD_summary.studies.png"),
                                  height=15,
                                  width=12,
                                  dpi=500)
```


## Manual multi-plot 

Alternatively, you can gain more customizability by merging your plots together manually with `patchwork`.  

Make sure to set `label_yaxis=F` for some of the plots to avoid redundancy and save space. 

```{r Multi-plot, eval=F} 
library(patchwork)

# Merge
row1 <- (patchwork::plot_spacer() + 
           bin_plot$plot ) +patchwork::plot_layout(widths = c(.4,.6))
row2 <- (gg_egene$plot + 
           gg_CS$plot + 
           gg_missense$plot +
           gg_peaks$plot   
           ) + patchwork::plot_layout(widths = c(.125,.3,.05,1))

gg_merged <- (row1  / row2)  + patchwork::plot_layout(heights = c(.15,1),ncol = 1) 
print(gg_merged)

```
 
# Session info 

<details>

```{r Session Info, attr.output='style="max-height: 200px;"'}
utils::sessionInfo()
```

</details>

---
title: "QTL pipeline vignette" 
author: "<h4>Author: <i>Brian M. Schilder</i></h4>" 
date: "<h4>Updated: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h4>"
output:
  BiocStyle::html_document
vignette: >
    %\VignetteIndexEntry{QTL_pipeline_vignette} 
    %\usepackage[utf8]{inputenc}
    %\VignetteEngine{knitr::rmarkdown}
---


```{r setup, include = T, error = T, message=F} 
root.dir <- tempdir()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = root.dir,
  fig.height = 12,
  fig.width = 10
)  
knitr::opts_knit$set(root.dir = root.dir)

library(echolocatoR) 
library(dplyr)
results_dir <- file.path(root.dir,"results") 
dir.create(results_dir, showWarnings = F, recursive = F)
```

 
# QTL pipeline

- Here, we will use GWAS-eQTL colocalization results provided via the [*echolocatoR Fine-mapping Portal*](https://github.com/RajLabMSSM/Fine_Mapping_Shiny), from the preprint:
> K de Paiva Lopes, GJL Snijders, J Humphrey, A Allan, M Sneeboer, E Navarro, BM Schilder…T Raj (2020) Atlas of Genetic Effects in Human Microglia Transcriptome across Brain Regions, Aging and Disease Pathologies. bioRxiv; [https://doi.org/10.1101/2020.10.27.356113](https://doi.org/10.1101/2020.10.27.356113).


## Prepare `top_SNPs` data.frame   

- In this case, we don't have a top SNPs file ready.
So we're just going to make one directly from the full summary stats file itself (*NOTE*: You can only use this approach if you can fit the entire file in memory).  
- In this case, you'll want to make sure to set `grouping_vars=c("Locus","Gene")` so that you get top SNPs for each eGene-locus pair (not just one SNP per locus).  

```{r  Prepare `top_SNPs` data.frame, eval=F}
fullSS_path <- file.path(root.dir,"Microglia_all_regions_Kunkle_2019_COLOC.snp-level_select.tsv.gz")
## Download example data 
if(!file.exists(fullSS_path)){
  download.file("https://github.com/RajLabMSSM/Fine_Mapping_Shiny/raw/master/www/Microglia_all_regions_Kunkle_2019_COLOC.snp-level_select.tsv.gz", fullSS_path)
}


top_SNPs <- import_topSNPs(  
  topSS = fullSS_path,
  chrom_col = "chr", 
  position_col = "pos",
  snp_col="snp",
  pval_col="gwas.pvalues", 
  effect_col="gwas.beta", 
  gene_col="gene",
  locus_col = "Locus",
  grouping_vars = c("Locus","Gene")) 
head(top_SNPs)
```

## eGene-locus list  

- In `finemap_loci()` we usually supply `loci=` with a list of locus names. 
However,when fine-mapping GWAS-QTL jointly, you'll need to specify which
QTL eGene-GWAS locus pairs you want to run.
- You can easy get this named list of all eGene-locus pairs with `gene_locus_list()`. The returned list's values are the GWAS loci, while the names are the QTL eGenes.  
- If you don't want to test all pairs, you can filter the `top_SNPs` object first. In this example, I'm going to remove any pairs that don't have matching eGene-locus names.  
We'll also limit the anaylsis to just 3 loci.  

```{r, eval=F} 
loci <- gene_locus_list(top_SNPs)
# We'll just fine-map 2 loci to start
loci <- loci[1:3] 
loci
```


## Run fine-mapping pipeline  

**NOTE**: Currently this `finemap_loci()` will only fine-map the GWAS results. However you can still plot the fine-mapped GWAS results using the `QTL_prefixes` argument.  

```{r, eval=F}
Kunkle_2019.microgliaQTL <- finemap_loci(# GENERAL ARGUMENTS 
                                        top_SNPs = top_SNPs,  
                                        results_dir = results_dir,
                                        
                                        loci = loci,
                                        dataset_name = "Kunkle_2019.microgliaQTL",
                                        dataset_type = "GWAS",  
                                        sample_size = 21982 + 41944,
                                        proportion_cases = 21982 / 41944,
                                        force_new_subset = F,
                                        force_new_LD = F,
                                        force_new_finemap = F,
                                        remove_tmps = T,
                                          
                 # SUMMARY STATS ARGUMENTS
                 fullSS_path = fullSS_path,
                 query_by ="fullSS",
                 chrom_col = "chr", 
                 position_col = "pos", 
                 snp_col = "snp",
                 pval_col = "gwas.pvalues", 
                 effect_col = "gwas.beta", 
                 stderr_col = "gwas.varbeta", 
                 MAF_col = "gwas.MAF", 
                 gene_col = "gene",
                 
                 # QTL prefixes
                 QTL_prefixes = c("qtl."),
                 
                 # FILTERING ARGUMENTS
                 bp_distance = 500000, #100000,
                 min_MAF = 0.001, 
                 
                 # FINE-MAPPING ARGUMENTS
                 finemap_methods = c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
                 n_causal = 5,
                 PP_threshold = .95,
                 
                 # LD ARGUMENTS 
                 LD_reference = "1KGphase1",
                 superpopulation = "EUR",
                 download_method = "axel",
                 
                 # PLOT ARGUMENTS 
                 plot.types=c("simple"),
                 plot.zoom = "1x" 
                 )
```


# Session info

<details> 

```{r Session Info, attr.output='style="max-height: 200px;"'}
utils::sessionInfo()
```

</details>

<br>
---
title: "Getting Started" 
author: "<h4>Author: <i>Brian M. Schilder</i></h4>" 
date: "<h4>Updated: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h4>"
output:
  BiocStyle::html_document
vignette: >
    %\VignetteIndexEntry{echolocatoR} 
    %\usepackage[utf8]{inputenc}
    %\VignetteEngine{knitr::rmarkdown}
---


```{r, include = F, message=F, error = T} 
root.dir <- tempdir()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = root.dir,
  fig.height = 12,
  fig.width = 10
)  
knitr::opts_knit$set(root.dir = root.dir)  
```

```{r setup}
library(echolocatoR) 

root.dir <- tempdir()
```

# Full pipeline  

All examples below use data from the Parkinson's disease GWAS by Nalls et al. (2019).  

## Prepare `top_SNPs` data.frame   

* To enable rapid fine-mapping of many loci, you can create a `top_SNPs` data.frame  
which contains the position of the lead/index SNP within each locus.
* `finemap_loci()` (see next step) will then use this info to extract subsets of the   
full GWAS/QTL summary statistics using windows centered on each lead/index SNP.
* The `topSS` argument can either be a data.frame, or a path to a topSS file saved somehwere. Most common tabular data formats (e.g. .tsv, .csv, .xlsx) are accepted.  



```{r  Prepare `top_SNPs` data.frame} 
data("Nalls_top_SNPs");
top_SNPs <- import_topSNPs(
  # topSS = "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx",
  topSS = Nalls_top_SNPs,
  munge = TRUE,
  ref_genome = "GRCH37", 
  pval_col="P, all studies", 
  effect_col="Beta, all studies", 
  gene_col="Nearest Gene", 
  locus_col = "Nearest Gene",
  grouping_vars = c("Locus"),
  remove_variants = "rs34637584") 
head(top_SNPs)
```



## Path to full summary stats file  

* Since a full GWAS summary stats file would be too large to include within *echolocatoR*,  
we instead provide an example subset of the full summary stats.  

* To simulate how you'd actually use your own full summary stats file, we will save our example dataset to your computer (you can change the path to wherever you like). 

* We highly recommend munging your full summary stats using the Bioconductor package [`MungeSumstats`](https://github.com/neurogenomics/MungeSumstats) first. It's easy to use and very robust. It also means you don't have to provide most column mapping arguments in `finemap_loci` when `munged=TRUE`. 

Here's an example of how to munge your full summary stats file: 

```
fullSS_path <- example_fullSS(munged = FALSE)
fullSS_path <- MungeSumstats::format_sumstats(path = fullSS_path, ref_genome = "GRCH37")
```

We have already munged the following example summary stats for you.

```{r fullSS} 
fullSS_path <- example_fullSS(munged = TRUE)
```

## Run fine-mapping pipeline  

For a full description of all arguments, see `?finemap_loci`.  

Here are some key arguments:  

* *results_dir*: Where you want to store all of your results.  
* *finemap_methods*: Which fine-mapping methods you want to run (currently includes ABF, FINEMAP, SUSIE, POLYFUN_SUSIE, and COJO).  
* *bp_distance*: Controls window size. Specifically, `bp_distance` is the number of basepairs upstream/downstream you want to extract for each locus. For example, if you want a 2Mb window (+/- 1Mb from the lead/index SNP in `top_SNPs`), set `bp_distance=1e+06`.  
* *plot.zoom*: Zoom in/out from the center of each locus when producing the multiview plot.  
You can adjust this separately from `bp_distance` so that you don't have rerun the whole pipeline each time (locus subsets, LD matrices, and fine-mapping results are all automatically saved in locus-specific folders).  

**NOTE**: This example assumes you have already installed tabix (via htslib) and bcftoools.

**WARNING**: Please use the full absolute paths (instead of relative paths) wherever possible (e.g. `results_dir`). This is especially important for the tool *FINEMAP*.

```{r Run fine-mapping pipeline, eval=F}

Nalls23andMe_2019.results <- finemap_loci(# GENERAL ARGUMENTS 
                                          top_SNPs = top_SNPs,  
                                          #  It's best to give absolute paths
                                          results_dir = file.path(root.dir,"results"),
                                          loci = c("BST1","MEX3C"),# top_SNPs$Locus,
                                          dataset_name = "Nalls23andMe_2019",
                                          dataset_type = "GWAS",  
                                          force_new_subset = TRUE,
                                          force_new_LD = FALSE,
                                          force_new_finemap = TRUE,
                                          remove_tmps = FALSE,
                                          
                                          # Munge full sumstats first
                                          munged = TRUE,
                                          
                                          # SUMMARY STATS ARGUMENTS
                                          fullSS_path = fullSS_path,
                                          fullSS_genome_build = "hg19",
                                          query_by ="tabix",  
                                          MAF_col = "calculate",   
                                         
                                          # FILTERING ARGUMENTS
                                          ## It's often desirable to use a larger window size 
                                          ## (e.g. 2Mb which is bp_distance=500000*2), 
                                          ## but we use a small window here to speed up the process. 
                                          bp_distance = 10000,#500000*2,
                                          min_MAF = 0.001,  
                                          trim_gene_limits = FALSE,
                                         
                                          # FINE-MAPPING ARGUMENTS
                                          ## General
                                          finemap_methods = c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"), 
                                          n_causal = 5,
                                          PP_threshold = .95, 
                                         
                                          # LD ARGUMENTS 
                                          LD_reference = "1KGphase1",#"UKB",
                                          superpopulation = "EUR",
                                          download_method = "axel",
                                         
                                          # PLOT ARGUMENTS 
                                          ## general   
                                          plot.types = c("fancy"),
                                          ## Generate multiple plots of different window sizes; 
                                          ### all SNPs, 4x zoomed-in, and a 50000bp window
                                          plot.zoom = c("all","4x","10x"),
                                          ## XGR
                                          # plot.XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes"), 
                                          ## Roadmap
                                          plot.Roadmap = FALSE,
                                          plot.Roadmap_query = NULL,
                                          # Nott et al. (2019)
                                          plot.Nott_epigenome = TRUE, 
                                          plot.Nott_show_placseq = TRUE, 
                                         
                                          verbose = TRUE
                                          )
```


```{r, eval=F}
dataset <-  file.path(root.dir,"results/GWAS/Nalls23andMe_2019/")
merged_dat <- echolocatoR::merge_finemapping_results(dataset = dataset,
                                                     minimum_support = 0)
```



# Session info  

<details>

```{r Session Info, attr.output='style="max-height: 200px;"'}
utils::sessionInfo()
```

</details>


---
title: "Plotting vignette" 
author: "<h4>Author: <i>Brian M. Schilder</i></h4>" 
date: "<h4>Updated: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h4>"
output:
  BiocStyle::html_document
vignette: >
    %\VignetteIndexEntry{plotting_vignette} 
    %\usepackage[utf8]{inputenc}
    %\VignetteEngine{knitr::rmarkdown}
---

```{r, include = T, error = T, dpi=300}
root.dir <- tempdir()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = root.dir, 
  fig.height = 12,
  fig.width = 10
)  
knitr::opts_knit$set(
  root.dir = root.dir
  )
# knitr::opts_chunk$get("root.dir")
# devtools::build_vignettes(quiet = F)
```

```{r setup}
library(echolocatoR) 
```

# Plotting loci with *echolocatoR*

*echolocatoR* contains various functions that can be used separately  
from the comprehensive `finemap_loci()` pipeline.  

Generate a multi-view plot of a given locus using `PLOT.locus()`.  

* You can mix and match different tracks and annotations using the different arguments 
(see `?PLOT.locus` for details).  

The plot is centered on the lead/index SNP. If a list is supplied to plot.xoom
* `PLOT.locus()` returns a series of `ggplot` objects bound together with [`patchwork`](https://patchwork.data-imaginist.com). One can further modify this object using `ggplot2` functions like `+ theme()`.
  + The modifications will be applied to all tracks at once.  
  
* Save a high-resolution versions the plot by setting `save_plot=T`.  
  + Further increase resolution by adjusting the `dpi` argument (*default=300*).
  + Files are saved in *jpg* format by default, but users can specify their preferred file format (e.g. `file_format="png"`)
  + Adjust the `height` and `width` of the saved plot using these respective arguments.
  + The plot will be automatically saved in the locus-specific directory as:  
  *multiview_<locus>_<plot.zoom>.jpg*.
   
## Load example data   

Load example dataset of the results from fine-mapping the BST1 locus with `finemap_loci()`.
Original data comes from the recent Nalls et al. (2019) Parkinson's disease GWAS (see `?BST1` for details).

```{r Load example data, eval=F}
library(ggplot2)
data("BST1"); data("LD_matrix"); data("locus_dir");
locus_dir <- file.path(root.dir,locus_dir)
finemap_DT <- BST1  
LD_matrix <- BST1_LD_matrix
LD_reference <- "UKB" # Used for naming saved plots
plot.zoom = "10x"
show_plot <- F


# locus <- "DNAH17"
# locus_dir <- file.path("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019",locus)
# LD_matrix <- readRDS(file.path(locus_dir,"plink/UKB_LD.RDS"))
# finemap_DT <- data.table::fread(file.path(locus_dir,"Multi-finemap/Multi-finemap_results.txt"))
# finemap_DT <- find_consensus_SNPs(finemap_DT) 
# finemap_DT <- update_CS_cols(finemap_DT)
```


## Full window

```{r trk_plot, attr.output='style="max-height: 200px;"', eval=F}
trk_plot <- PLOT.locus(finemap_dat=finemap_DT, 
                       LD_matrix=LD_matrix, 
                       LD_reference=LD_reference,
                       locus_dir=locus_dir,  
                       save_plot=F,
                       show_plot=show_plot,
                       plot.zoom=plot.zoom) 
```

```{r, eval=F}
print(trk_plot)
```

## At multiple zooms

* You can easily generate the same locus plot at multiple zoomed in views by supplying a list to `plot.zoom`.  
* This list can be composed of zoom multipliers (e.g. `c("1x", "2x")`), window widths in units of basepairs (e.g. `c(5000, 1500)`), or a mixture of both (e.g. `c("1x","4x", 5000, 2000)`). 
* Each zoom view will be saved individually with its respective scale as the suffix (e.g. `multiview.BST1.UKB.4x.jpg`).  
* Each zoom view is stored as a named item within the returned list.  

```{r trk_plot_zoom, attr.output='style="max-height: 200px;"', eval=F}
trk_plot_zooms <- PLOT.locus(finemap_dat=finemap_DT, 
                             LD_matrix=LD_matrix, 
                             LD_reference=LD_reference,
                             locus_dir=locus_dir,  
                             save_plot=F,
                             show_plot=show_plot,
                             plot.zoom = c("1x","5x","10x")) 
names(trk_plot_zooms) # Get zoom view names
```

```{r, eval=F}
print(trk_plot_zooms)
```


## Return as list  

* For even further control over each track of the multi-view plot, specify `PLOT.locus(..., return_list=T)` to instead return a named list (nested within each zoom view list item) of `ggplot` objects which can each be modified individually. 
* Once you've made your modifications, you can then bind this list of plots back together with `patchwork::wrap_plots(tracks_list, ncol = 1)`.  

```{r trk_plot_list, attr.output='style="max-height: 200px;"', eval=F}
trk_plot_list <- PLOT.locus(finemap_dat=finemap_DT, 
                             LD_matrix=LD_matrix, 
                             LD_reference=LD_reference,
                             locus_dir=locus_dir,  
                             save_plot=F,
                             show_plot=show_plot,
                             plot.zoom=plot.zoom,
                             return_list=T)  
```

```{r extract view, eval=F}
view1_list <- trk_plot_list[[plot.zoom]]
names(view1_list) # Get track names from a particular zoom view
```

Modify a specific tracks within a view. 

```{r modify track, eval=F} 
# Modify your selected track
modified_track <- view1_list$GWAS + 
                      labs(title = "Modified GWAS") + 
                      theme_dark() +
                      theme(title = element_text(hjust = .5))
# Put it back into your track list
view1_list[["GWAS"]] <- modified_track
# Remove a plot you don't want
view1_list[["Genes"]] <- NULL
# Specify the relative heights of each track (make sure it matches your new # of plots!)
track_heights <- c(.3,.1,.3,1)

# Bind them together and plot
fused_plot <- patchwork::wrap_plots(view1_list, 
                                    heights = track_heights,
                                    ncol = 1)
print(fused_plot)
```


## Using XGR annotations   

* Whenever you use annotation arguments (e.g. `XGR_libnames`,`Roadmap`,`Nott_epigenome`)
the annotations that overlap with your locus will automatically be saved as `GRanges` objects in a locus-specific subdirectory:  
*results/<dataset_type>/<dataset_name>/<locus>/annotation* 
* If a selected annotation has previously been downloaded and stored for that locus, `PLOT.locus()` will automatically detect and import it to save time.   

```{r trk_plot.xgr, attr.output='style="max-height: 200px;"', eval=F}
trk_plot.xgr <- PLOT.locus(finemap_dat=finemap_DT, 
                           LD_matrix=LD_matrix, 
                           LD_reference=LD_reference,
                           locus_dir=locus_dir, 
                           
                           XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes"), 
                           save_plot=F,
                           show_plot=show_plot,
                           plot.zoom=plot.zoom)
```

```{r, eval=F}
print(trk_plot.xgr)
```

## Using [Roadmap](http://www.roadmapepigenomics.org/data/tables/brain) annotations

* Using the `Roadmap=T` and `Roadmap_query="<query>"` arguments searches the Roadmap for chromatin mark data across various cell-types, cell-lines and tissues.  
* Note that Roadmap queries requires `tabix` to be installed on your machine, or within a conda environment (`conda_env = "echoR"`). 
* Parallelizing these queries across multiple thredas speeds up this process (`nThread=<n_cores_available>`), as does reusing previously stored data which is automatically saved to the locus-specific subfolder (`<dataset_type>/<dataset_name>/<locus>/annotations/Roadmap.ChromatinMarks_CellTypes.RDS`) .

```{r trk_plot.roadmap, attr.output='style="max-height: 200px;"', eval=F} 
trk_plot.roadmap <- PLOT.locus(finemap_dat=finemap_DT, 
                               LD_matrix=LD_matrix, 
                               LD_reference=LD_reference,
                               locus_dir=locus_dir,  
                               
                               Roadmap=T, 
                               Roadmap_query="monocyte", 
                               
                               save_plot=F, 
                               show_plot=show_plot,
                               plot.zoom="5x", 
                               nThread = parallel::detectCores()-1,
                               conda_env = "echoR")
```

```{r, eval=F}
print(trk_plot.roadmap)
```

## Using [Nott_2019](https://science.sciencemag.org/content/366/6469/1134.abstract)   annotations   

* Query and plot brain cell type-specific epigenomic assays from 
[Nott et al. (Science, 2019)](https://science.sciencemag.org/content/366/6469/1134.abstract)  
(see `?NOTT_2019.bigwig_metadata` for details).

```{r trk_plot.nott_2019, attr.output='style="max-height: 200px;"', eval=F} 
trk_plot.nott_2019 <- PLOT.locus(finemap_dat=finemap_DT,  
                                 LD_matrix=LD_matrix, 
                                 LD_reference=LD_reference,
                                 locus_dir=locus_dir, 
                                 
                                 Nott_epigenome=T,  
                                 Nott_binwidth = 200,
                                 Nott_regulatory_rects = T, 
                                 Nott_show_placseq = T,
                                 
                                 save_plot=F,
                                 show_plot=show_plot,
                                 plot.zoom=plot.zoom) 
```

```{r, eval=F}
print(trk_plot.nott_2019)
```

## Using QTL datasets  

* Plot multiple QTL p-value columns (or really P-value columns from any kind of datset).  
* Each QTL dataset will be plotted as a new track.

```{r trk_plot.QTL, attr.output='style="max-height: 200px;"', eval=F}
# Make fake QTL P-values for the sake a demonstration
finemap_DT$fake_eQTL.P <- finemap_DT$P  * c(1,.9,.7)
finemap_DT$fake_sQTL.P <- finemap_DT$P  * c(1,.8,.5)

trk_plot.qtl <- PLOT.locus(finemap_dat=finemap_DT, 
                           LD_matrix=LD_matrix, 
                           LD_reference=LD_reference,
                           locus_dir=locus_dir,
                           
                           QTL_prefixes=c("fake_eQTL.","fake_sQTL."), 
                           
                           save_plot=F,
                           show_plot=show_plot,
                           plot.zoom=plot.zoom)
```

```{r, eval=F}
print(trk_plot.qtl)
```

# Session info  

<details>  

```{r Session Info, attr.output='style="max-height: 200px;"'}
utils::sessionInfo()
```

<details>

---
title: "hexSticker"
author: "<h4>Author: <i>Brian M. Schilder</i></h4>" 
date: "<h4>Most recent update: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h4>"
output: rmarkdown::html_vignette
editor_options: 
  chunk_output_type: inline
vignette: >
  %\VignetteIndexEntry{hexSticker}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(echo = T, fig.width = 7, fig.height = 5, root.dir=here::here())
knitr::opts_knit$set(root.dir=here::here())
```

```{r setup, include=TRUE, message=FALSE} 
library(hexSticker)

library(dplyr)
library(ggplot2)
```

You can make awesome hex stickers for your R packages using [hexSticker](https://github.com/GuangchuangYu/hexSticker). 

# phenomix


```{r}   
img1 <- "/Desktop/echolocatoR/images/echolocatoR_logo.png"
# img2 <- "../images/phoenix_pixel.png"

sticker(subplot = img1, package="",
        # url = "https://github.com/RajLabMSSM/echolocatoR", u_color = "white", u_size = 3,
        p_size=20, s_x=1, s_y=1.1,  s_width = .8,
        h_fill = "#686ea6",#"#2a3182",
        h_color = "#4ee8e8", 
        spotlight = TRUE, #l_height = 20, 
        white_around_sticker=FALSE, 
        filename= "/Desktop/echolocatoR/inst/hex/hex.png")

```


# Session Info 

<details> 

```{r Session Info}
utils::sessionInfo()
```

</details>  

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{query_by_probe}
\alias{query_by_probe}
\title{Use \emph{awk} to query locus subsets.}
\usage{
query_by_probe(
  fullSS_path,
  subset_path,
  gene,
  gene_col,
  chrom_col,
  file_sep,
  probe_path,
  coordinates_merged = T,
  location_sep = ":"
)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{gene_col}{For QTL studies, the name of the [e]gene column in the full summary stats file (\emph{default: "gene"}).
This column will be used for filtering summary stats if supplying a named list of gene:Locus pairs to \code{loci}.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{file_sep}{The separator in the full summary stats file.
This parameter is only necessary if \code{query_by!="tabix"}.}

\item{probe_path}{The location of the file containing translations between probe IDs and gene symbols.
Only used for certain eQTL datasets.}

\item{coordinates_merged}{Whether \strong{CHR} and \strong{POS} are merged into one column (e.g. chr12:12209944).}

\item{location_sep}{The separator character when \strong{CHR} and \strong{POS} are merged into one column (e.g. ":" when formatted like chr12:12209944).}
}
\description{
To be used when probe name (but not gene name) is a column in the full summary stats file.
Requires a gene-probe key file (\strong{probe_path}).
More commonly useful for QTL full summary stats files.
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORCES_2020.R
\name{CORCES_2020.get_ATAC_peak_overlap}
\alias{CORCES_2020.get_ATAC_peak_overlap}
\title{Get overlap between datatable of SNPs and scATAC peaks}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.get_ATAC_peak_overlap(
  finemap_dat,
  FDR_filter = NULL,
  add_cicero = T,
  cell_type_specific = T,
  verbose = T
)
}
\description{
Can optionally add \code{Cicero} coaccessibility scores,
which are also derived from scATAC-seq data.
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataframe_2_vcf.R
\name{dataframe_2_vcf}
\alias{dataframe_2_vcf}
\title{Convert a dataframe to a vcf}
\usage{
dataframe_2_vcf(
  subset_DT,
  samplename = "GWAS",
  reference_fasta = "/pd-omics/tools/polyfun/reference_fasta/hg19.fa.gz",
  output_vcf = "./GWAS_converted.vcf",
  conda_env = "echoR"
)
}
\description{
Convert a dataframe to a vcf
}
\examples{
data("merge_DT")

}
\seealso{
Other utils: 
\code{\link{LIFTOVER}()},
\code{\link{rbind_GRanges}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.filter_assays}
\alias{XGR.filter_assays}
\title{Filter assays}
\usage{
XGR.filter_assays(gr.lib, n_top_assays = 5)
}
\description{
Identify the assays with the most annotations in the locus.
Then only keep these assays
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifbreakR.R
\name{MOTIFBREAKR.plot}
\alias{MOTIFBREAKR.plot}
\title{Plot \code{\link{motifbreakR}} results}
\source{
\strong{Publication:}
\url{https://pubmed.ncbi.nlm.nih.gov/26272984/}

\strong{GitHub:}
\url{https://github.com/Simon-Coetzee/MotifBreakR}
}
\usage{
MOTIFBREAKR.plot(
  mb.results,
  mb.filter = NULL,
  rsid = NULL,
  effect = c("strong", "weak"),
  save_dir = NULL,
  height = 3,
  width = 7
)
}
\description{
Plot \code{\link{motifbreakR}} results
}
\examples{
\dontrun{
# Example from motifbreakR
library(motifbreakR)
data("example.results")
library(echolocatoR)
data("merged_DT")
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR"

mb.results<- readRDS(file.path(root,"DNAH17.pvalues.RDS"))

# Get overlap with PEAKS
# PEAKS <- NOTT_2019.get_epigenomic_peaks(convert_to_GRanges = T)
REGIONS <- NOTT_2019.get_regulatory_regions(as.granges = T)
gr.hits <- GRanges_overlap(dat1 = subset(merged_DT, Consensus_SNP), chrom_col.1 = "CHR", start_col.1 = "POS", end_col.1 = "POS", dat2 = REGIONS)
mb.filter <- MOTIFBREAKR.filter(merged_DT = merged_DT, mb.results = mb.results, pct_threshold=NULL)

subset(mb.filter, SNP_id \%in\% gr.hits$SNP)
plot_paths <- MOTIFBREAKR.plot(mb.results=mb.results, mb.filter=mb.filter, save_dir="~/Desktop")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FINEMAP.R
\name{FINEMAP.construct_master}
\alias{FINEMAP.construct_master}
\title{Construct the \code{FINAMAP} master file}
\source{
\url{http://www.christianbenner.com}
}
\usage{
FINEMAP.construct_master(
  locus_dir,
  n_samples,
  dataset_number = 1,
  file.k = NA,
  verbose = T
)
}
\description{
Creates and saves the master file
which tells \code{FINEMAP} where to find each input file.
}
\examples{
data("locus_dir");
master_path <- FINEMAP.construct_master(locus_dir=locus_dir, n_samples=25000)
}
\seealso{
Other FINEMAP: 
\code{\link{FINEMAP.construct_data}()},
\code{\link{FINEMAP.find_executable}()},
\code{\link{FINEMAP.process_results}()},
\code{\link{FINEMAP}()}
}
\concept{FINEMAP}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.get_regulatory_regions}
\alias{NOTT_2019.get_regulatory_regions}
\title{Plot brain cell-specific epigenomic data}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.get_regulatory_regions(as.granges = F, nThread = 1, verbose = T)
}
\description{
Plot brain cell-specific epigenomic data
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.cell_type_specificity}
\alias{SUMMARISE.cell_type_specificity}
\title{Get cell-type-specifity score for each cell type}
\usage{
SUMMARISE.cell_type_specificity(
  plot_dat,
  merged_dat,
  min_count = NULL,
  top_celltype_only = F,
  label_yaxis = T,
  y_lab = NULL,
  show_genes = F,
  x_strip_angle = 40,
  show_plot = T
)
}
\description{
Aggregate SNP overlap across various epigenomic datasets
and then identify the number of SNPs overlapping by each cell type
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.get_locus_vcf_folder}
\alias{LD.get_locus_vcf_folder}
\title{Get VCF storage folder}
\usage{
LD.get_locus_vcf_folder(locus_dir = NULL)
}
\description{
Get VCF storage folder
}
\examples{
data("locus_dir")
vcf_folder <- LD.get_locus_vcf_folder(locus_dir=locus_dir)
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snps_to_condition.R
\name{snps_to_condition}
\alias{snps_to_condition}
\title{Identify SNPs to condition on.}
\usage{
snps_to_condition(conditioned_snps, top_SNPs, loci)
}
\description{
When running conditional analyses (e.g. \emph{GCTA-COJO}),
this functions automatically identifies SNP to condition on.
}
\seealso{
Other SNP filters: 
\code{\link{filter_snps}()},
\code{\link{gene_trimmer}()},
\code{\link{limit_SNPs}()},
\code{\link{subset_common_snps}()}
}
\concept{SNP filters}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liftover.R
\name{LIFTOVER}
\alias{LIFTOVER}
\title{Lift genome across builds}
\usage{
LIFTOVER(
  dat,
  build.conversion = "hg19.to.hg38",
  chrom_col = "CHR",
  start_col = "POS",
  end_col = "POS",
  chr_format = "NCBI",
  return_as_granges = T,
  verbose = T
)
}
\arguments{
\item{build_conversion}{"hg19.to.hg38" (\emph{default}) or "hg38.to.hg19.}
}
\description{
Lift genome across builds
}
\examples{
data("BST1")
gr.lifted <- LIFTOVER(dat=BST1, build.conversion="hg19.to.hg38")
}
\seealso{
Other utils: 
\code{\link{dataframe_2_vcf}()},
\code{\link{rbind_GRanges}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MEX3C}
\alias{MEX3C}
\title{\emph{echolocatoR} output example: MEX3C locus}
\format{
data.table
\describe{
  \item{SNP}{SNP RSID}
  \item{CHR}{Chromosome}
  \item{POS}{Genomic position (in basepairs)}
  \item{...}{Optional: extra columns}
}
}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
MEX3C
}
\description{
An example results file after running
\code{\link[=finemap_loci]{finemap_loci()}} on the \emph{MEX3C} locus.
}
\details{
Data originally comes from the Parkinson's disease GWAS
by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
}
\examples{
\dontrun{
root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/MEX3C/Multi-finemap"
MEX3C <- data.table::fread(file.path(root_dir,"Multi-finemap_results.txt"))
MEX3C <- update_cols(finemap_dat=MEX3C)
MEX3C <- find_consensus_SNPs(finemap_dat=MEX3C)
usethis::use_data(MEX3C, overwrite = T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CORCES_2020.scATACseq_celltype_peaks}
\alias{CORCES_2020.scATACseq_celltype_peaks}
\title{scATACseq cell type-specific peaks from human brain tissue}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 221062 rows and 13 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.scATACseq_celltype_peaks
}
\description{
Each row represents an individual peak identified from the feature binarization analysis (see methods).
}
\details{
Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
Specifically: \emph{STable6_Features_scATAC-seq_celltype_Peaks}
}
\examples{
\dontrun{
dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable6_Features_scATAC-seq_celltype_Peaks.xlsx", skip = 15)
CORCES_2020.scATACseq_celltype_peaks <- data.table::data.table(dat)
usethis::use_data(CORCES_2020.scATACseq_celltype_peaks)
}
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpliceAI.R
\name{SPLICEAI.subset_precomputed_tsv_iterate}
\alias{SPLICEAI.subset_precomputed_tsv_iterate}
\title{Make multiple queries to genome-wide \emph{SpliceAI} results file (tsv format)}
\source{
\href{https://github.com/Illumina/SpliceAI}{GitHub}
\href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
}
\usage{
SPLICEAI.subset_precomputed_tsv_iterate(
  sumstats_paths,
  precomputed_path = "/pd-omics/data/spliceAI/spliceai_scores.raw.snv.hg19.tsv.gz",
  nThread = 4,
  merge_data = T,
  drop_na = T,
  filtered = F,
  save_path = "./spliceAI_subset.tsv.gz"
)
}
\description{
Make multiple queries to genome-wide \emph{SpliceAI} results file (tsv format)
}
\examples{
\dontrun{
root.pd <- "/sc/arion/projects/pd-omics"
precomputed_path <- file.path(root.pd,"data/spliceAI/spliceai_scores.raw.snv.hg19.tsv.gz")
sumstats_paths <- list.files(file.path(root.pd, "/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"), pattern = "*.UKB_LD.Multi-finemap.tsv.gz", recursive = T, full.names = T)
DAT <- SPLICEAI.subset_precomputed_tsv_iterate(sumstats_paths=sumstats_paths)
}
}
\seealso{
Other SpliceAI: 
\code{\link{SPLICEAI.plot}()},
\code{\link{SPLICEAI.run}()},
\code{\link{SPLICEAI.snp_probs}()},
\code{\link{SPLICEAI.subset_precomputed_tsv}()},
\code{\link{SPLICEAI.subset_precomputed_vcf}()}
}
\concept{SpliceAI}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.munge_summ_stats}
\alias{POLYFUN.munge_summ_stats}
\title{Munge summary stats}
\usage{
POLYFUN.munge_summ_stats(
  polyfun = NULL,
  fullSS_path,
  locus_dir,
  sample_size = NULL,
  min_INFO = 0,
  min_MAF = 0.001,
  force_new_munge = F,
  conda_env = "echoR"
)
}
\description{
Munge summary stats
}
\examples{
\dontrun{
data("genome_wide_dir");
fullSS_path <- example_fullSS(fullSS_path="~/Desktop/Nalls23andMe_2019.fullSS_subset.tsv")
munged_path <- POLYFUN.munge_summ_stats(fullSS_path=fullSS_path, locus_dir=genome_wide_dir, force_new_munge=T)
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.import_annotations}
\alias{XGR.import_annotations}
\title{Download XGR annotations}
\usage{
XGR.import_annotations(
  gr.snp,
  anno_data_path = file.path("annotations", paste0("XGR_", lib.name, ".rds")),
  lib.name,
  save_xgr = T,
  annot_overlap_threshold = 5
)
}
\description{
Download XGR annotations
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{biomart_snp_info}
\alias{biomart_snp_info}
\title{Download SNP-wise annotations from Biomart}
\source{
\href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
}
\usage{
biomart_snp_info(
  snp_list,
  reference_genome = "grch37",
  attributes = c("refsnp_id", "allele", "chr_name", "chrom_start", "chrom_end",
    "chrom_strand", "ensembl_gene_stable_id", "consequence_type_tv",
    "polyphen_prediction", "polyphen_score", "sift_prediction", "sift_score",
    "reg_consequence_types", "validated"),
  verbose = T
)
}
\description{
Download SNP-wise annotations from Biomart
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SURE.R
\name{SURE.download_annotations}
\alias{SURE.download_annotations}
\title{Download \emph{SuRE} annotations}
\source{
\href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
\href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
\href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
\href{https://sure.nki.nl}{SNP-SuRE data browse}
}
\usage{
SURE.download_annotations(
  URL = "https://osf.io/vxfk3/download",
  output_dir = ".",
  nThread = 4,
  v = verbose,
  conda_env = "echoR"
)
}
\description{
Download \emph{SuRE} annotations
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_finemap.R
\name{find_consensus_SNPs}
\alias{find_consensus_SNPs}
\title{Adds fine-mapping summary columns}
\usage{
find_consensus_SNPs(
  finemap_dat,
  credset_thresh = 0.95,
  consensus_thresh = 2,
  sort_by_support = T,
  exclude_methods = NULL,
  top_CS_only = F,
  replace_PP_NAs = T,
  verbose = F
)
}
\arguments{
\item{consensus_thresh}{Threshold for determining \strong{Consensus_SNP} status.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Adds several columns that summarise the results across all fine-mapping tools that were run:
\describe{
\item{Support}{The number of tools in which the SNP was proposed in a credible set.}
\item{mean.PP}{The mean per-SNP PP across all fine-mapping tools used.}
\item{Consensus_SNP}{Whether or not the SNP was in the credible set of ≥ \code{consensus_thresh} SNPs (\emph{default=2}).}
}
}
\examples{
data("merged_DT")
merged_DT <- find_consensus_SNPs(finemap_dat=merged_DT, top_CS_only=F)
}
\seealso{
Other finemapping functions: 
\code{\link{multi_finemap}()}
}
\concept{finemapping functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.coloc_nominated_eGenes}
\alias{SUMMARISE.coloc_nominated_eGenes}
\title{Nominate target genes within each locus}
\usage{
SUMMARISE.coloc_nominated_eGenes(
  coloc_results,
  merged_dat,
  label_yaxis = T,
  y_lab = "Locus",
  x_lab = NULL,
  fill_var = "PP.H4",
  text_size = 2,
  PP_threshold = NULL,
  nThread = 4,
  show_plot = T,
  verbose = T
)
}
\description{
Across all GWAS-QTL colocalization tests across all studies,
take the eGene with the highest colocalziation probability (PP.H4)
and assign it as the most likely causal gene in that locus.
}
\details{
eQTL queries and colocalization test done with \pkg{catalogueR}.
}
\examples{
\dontrun{
data("merged_DT")
base_url <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
coloc_results_path <- file.path(base_url,"_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz")
gg_egene <- SUMMARISE.coloc_nominated_eGenes(coloc_results, merged_dat=merged_DT, fill_var=NULL)

# QTL
base_url <- "/sc/hydra/projects/ad-omics/microglia_omics/Fine_Mapping"
coloc_results_path <- file.path(base_url, "Kunkle_Microglia_all_regions/QTL_merged_coloc_results.snp.tsv.gz")
merged_dat <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/QTL/Microglia_all_regions/multiGWAS.microgliaQTL_finemapping.csv.gz")
gg_egene <- SUMMARISE.coloc_nominated_eGenes(coloc_results, merged_dat=merged_dat, fill_var=NULL)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MungeSumstats.R
\name{MUNGESUMSTATS.check_syn}
\alias{MUNGESUMSTATS.check_syn}
\title{Check whether a table has headers that can be mapped to an attribute in \pkg{MungeSumstats}}
\usage{
MUNGESUMSTATS.check_syn(dat, col_name = "P-value", col_type = "P")
}
\description{
Check whether a table has headers that can be mapped to an attribute in \pkg{MungeSumstats}
}
\concept{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{VALIDATION.super_plot}
\alias{VALIDATION.super_plot}
\title{Merge validation assays into a single plot}
\usage{
VALIDATION.super_plot(
 
    root = "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide",
  height = 10,
  width = 12,
  layout = "horiz",
  show_plot = T,
  save_plot = F
)
}
\description{
Merge validation assays into a single plot
}
\examples{
\dontrun{
plt.ALL <- VALIDATION.super_plot(root="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide")
}
}
\seealso{
Other VALIDATION: 
\code{\link{VALIDATION.bootstrap_multimetric}()},
\code{\link{VALIDATION.bootstrap_plot}()},
\code{\link{VALIDATION.bootstrap}()},
\code{\link{VALIDATION.compare_bootstrap_distributions}()},
\code{\link{VALIDATION.permute}()}
}
\concept{VALIDATION}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.dprime_table}
\alias{LD.dprime_table}
\title{Calculate LD (D')}
\usage{
LD.dprime_table(SNP_list, LD_folder, conda_env)
}
\description{
This appriach computes an LD matrix of D' (instead of r or r2) from a vcf.
See \code{\link{LD.run_plink_LD}} for a faster (but less flexible) alternative to computing LD.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.CS_bin_plot}
\alias{SUMMARISE.CS_bin_plot}
\title{Plot CS bin counts}
\usage{
SUMMARISE.CS_bin_plot(merged_dat, show_plot = T)
}
\description{
Plot CS bin counts
}
\examples{
data("merged_DT");
bin_plot <- SUMMARISE.CS_bin_plot(merged_dat=merged_DT)
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{results_report}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.run_plink_LD}
\alias{LD.run_plink_LD}
\title{Calculate LD (r or r2)}
\usage{
LD.run_plink_LD(
  bim,
  LD_folder,
  plink_prefix = "plink",
  r_format = "r",
  extract_file = NULL
)
}
\arguments{
\item{bim}{A bim file produced by \emph{plink}}

\item{LD_folder}{Locus-specific LD output folder.}

\item{r_format}{Whether to fill the matrix with \code{r} or \code{r2}.}
}
\description{
This appriach computes and LD matrix of r or r2 (instead of D') from a vcf.
See \code{\link{LD.dprime_table}} for a slower (but more flexible) alternative to computing LD.
}
\examples{
\dontrun{
data("LRRK2")
LD_folder <- "/Users/schilder/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/LRRK2/plink/saved"
bim_path <- file.path(LD_folder, "plink.bim");
bim <- data.table::fread(bim_path, col.names = c("CHR","SNP","V3","POS","A1","A2"), stringsAsFactors = F)
bim <- subset(bim, SNP \%in\% LRRK2$SNP)
ld.bin <- file.path(LD_folder, paste0("plink",".ld.bin"))
SNPs <- data.table::fread(file.path(LD_folder,"SNPs.txt"), col.names = 'RSID')
bin.vector <- readBin(ld.bin, what = "numeric", n=length(SNPs$RSID)^2)
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{haploR.HaploReg}
\alias{haploR.HaploReg}
\title{Download SNP-wise annotations from HaploReg}
\source{
\href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
}
\usage{
haploR.HaploReg(snp_list, verbose = T, chunk_size = NA)
}
\description{
Download SNP-wise annotations from HaploReg
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.locus.R
\name{PLOT.get_window_suffix}
\alias{PLOT.get_window_suffix}
\title{Determine the plot suffix indicating its window size}
\usage{
PLOT.get_window_suffix(finemap_dat, plot.zoom, verbose = T)
}
\description{
Determine the plot suffix indicating its window size
}
\examples{
data("BST1")
window_suffix <- PLOT.get_window_suffix(finemap_dat=BST1, plot.zoom=1000)
window_suffix <- PLOT.get_window_suffix(finemap_dat=BST1, plot.zoom=NULL)
window_suffix <- PLOT.get_window_suffix(finemap_dat=BST1, plot.zoom="all")
window_suffix <- PLOT.get_window_suffix(finemap_dat=BST1, plot.zoom="2x")
}
\seealso{
Other plot: 
\code{\link{PLOT.add_multitrack_lines}()},
\code{\link{PLOT.dot_summary}()},
\code{\link{PLOT.get_window_limits}()}
}
\concept{plot}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{merge_celltype_specific_epigenomics}
\alias{merge_celltype_specific_epigenomics}
\title{Merge all cell-type-specific epigenomics}
\usage{
merge_celltype_specific_epigenomics(keep_extra_cols = F)
}
\description{
Merges multiple cell-type-specific epigenomic datasets (Nott 2019, Corces 2020) into a single \code{GRAnges} object.
}
\examples{
gr.merged <- merge_celltype_specific_epigenomics()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.read_bin}
\alias{LD.read_bin}
\title{Create LD matrix from plink output}
\usage{
LD.read_bin(LD_dir)
}
\arguments{
\item{LD_dir}{Directory that contains the bin/bim files.}
}
\description{
Depending on which parameters you give \emph{plink} when calculating LD, you get different file outputs.
When it produces bin and bim files, use this function to create a proper LD matrix.
For example, this happens when you try to calculate D' with the \code{--r dprime-signed} flag (instead of just r).
}
\examples{
\dontrun{
locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/QTL/Microglia_all_regions/BIN1"
ld.matrix <- LD.read_bin(LD_dir=file.path(locus_dir, "LD"))
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.compute_enrichment}
\alias{IMPACT.compute_enrichment}
\title{Compute enrichment of IMPACT scores}
\usage{
IMPACT.compute_enrichment(annot_melt, locus = NULL)
}
\description{
Conduct IMPACT enrichment between SNP groups a fine-mapping methods.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpliceAI.R
\name{SPLICEAI.subset_precomputed_vcf}
\alias{SPLICEAI.subset_precomputed_vcf}
\title{Query genome-wide \emph{SpliceAI} results file (vcf format)}
\source{
\href{https://github.com/Illumina/SpliceAI}{GitHub}
\href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
}
\usage{
SPLICEAI.subset_precomputed_vcf(
  subset_DT,
 
    precomputed_path = "/pd-omics/data/spliceAI/whole_genome_filtered_spliceai_scores.vcf.gz",
  subset_vcf = "subset.vcf"
)
}
\description{
Query genome-wide \emph{SpliceAI} results file (vcf format)
}
\seealso{
Other SpliceAI: 
\code{\link{SPLICEAI.plot}()},
\code{\link{SPLICEAI.run}()},
\code{\link{SPLICEAI.snp_probs}()},
\code{\link{SPLICEAI.subset_precomputed_tsv_iterate}()},
\code{\link{SPLICEAI.subset_precomputed_tsv}()}
}
\concept{SpliceAI}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CORCES_2020.cicero_coaccessibility}
\alias{CORCES_2020.cicero_coaccessibility}
\title{Cicero_coaccessibility from human brain tissue}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 9795 rows and 14 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.cicero_coaccessibility
}
\description{
Cicero coaccessibility analysis for peaks that overlap SNPs derived from analysis of scATAC-seq data.
Each row represents an individual peak identified from the feature binarization analysis (see methods).
}
\details{
Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
Specifically: \emph{STable10_Coacessibility_Peak_loop_connection}, \emph{Cicero Coaccessibility} sheet.
Peak_ID_Peak1 - A unique number that identifies the peak across supplementary tables.

\strong{Column dictionary}:
\describe{
\item{hg38_Chromosome_Peak1}{The hg38 chromosome of the first loop Peak.}
\item{hg38_Start_Peak1}{The hg38 start position of the first loop Peak.}
\item{hg38_Stop_Peak1}{The hg38 stop position of the first loop Peak.}
\item{Width_Peak1}{The width of the first loop Peak.}
\item{Peak_ID_Peak2}{A unique number that identifies the peak across supplementary tables.}
\item{hg38_Chromosome_Peak2}{The hg38 chromosome of the second loop Peak.}
\item{hg38_Start_Peak2}{The hg38 start position of the second loop Peak.}
\item{hg38_Stop_Peak2}{The hg38 stop position of the second loop Peak.}
\item{Width_Peak2}{The width of the second loop Peak.}
\item{Coaccessibility}{The coaccessibility correlation for the given peak pair.}
\item{Peak1_hasSNP}{A boolean variable determining whether the first peak overlaps a SNP from our AD/PD GWAS analyses.}
\item{Peak2_hasSNP}{A boolean variable determining whether the second peak overlaps a SNP from our AD/PD GWAS analyses.}
}
}
\examples{
\dontrun{
dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable10_Coacessibility_Peak_loop_connection.xlsx", skip = 21, sheet=2)
CORCES_2020.cicero_coaccessibility <- data.table::data.table(dat)
usethis::use_data(CORCES_2020.cicero_coaccessibility)
}
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.h2_enrichment}
\alias{POLYFUN.h2_enrichment}
\title{Run heritability enrichment tests}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.h2_enrichment(h2_df, target_SNPs = NULL, fillNA = T)
}
\description{
Run heritability enrichment tests
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.plot_impact_score_compare}
\alias{IMPACT.plot_impact_score_compare}
\title{Compare IMPACT scores between loci}
\usage{
IMPACT.plot_impact_score_compare(
  loci = c("CD19", "TRIM40", "NUCKS1", "LRRK2", "MED12L", "MEX3C")
)
}
\description{
First, identify the top annotations per locus;
the annotation with the highest mean IMPACT score across all fine-mapped Consensus SNPS.
Then, plot the IMPACT scores for those locus-specific top annotations.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.locus.R
\name{PLOT.locus}
\alias{PLOT.locus}
\title{Generate a locus plot}
\usage{
PLOT.locus(
  finemap_dat,
  locus_dir,
  LD_matrix = NULL,
  LD_reference = NULL,
  dataset_type = "GWAS",
  color_r2 = T,
  method_list = c("ABF", "FINEMAP", "SUSIE", "POLYFUN_SUSIE"),
  track_order = NULL,
  track_heights = NULL,
  plot_full_window = T,
  dot_summary = F,
  QTL_prefixes = NULL,
  mean.PP = T,
  PP_threshold = 0.95,
  consensus_threshold = 2,
  sig_cutoff = 5e-08,
  gene_track = T,
  point_size = 1,
  point_alpha = 0.6,
  snp_group_lines = c("Lead", "UCS", "Consensus"),
  xtext = F,
  show.legend_genes = T,
  XGR_libnames = NULL,
  n_top_xgr = 5,
  Roadmap = F,
  Roadmap_query = NULL,
  n_top_roadmap = 7,
  annot_overlap_threshold = 5,
  zoom_exceptions_str = "*full window$|zoom_polygon",
  Nott_epigenome = F,
  Nott_regulatory_rects = T,
  Nott_show_placseq = T,
  Nott_binwidth = 200,
  Nott_bigwig_dir = NULL,
  save_plot = T,
  show_plot = T,
  genomic_units = "Mb",
  strip.text.y.angle = 0,
  max_transcripts = 1,
  plot.zoom = c("1x"),
  dpi = 300,
  height = 12,
  width = 10,
  plot_format = "jpg",
  save_RDS = F,
  return_list = F,
  conda_env = "echoR",
  nThread = 4,
  verbose = T
)
}
\description{
Generate a locus-specific plot with multiple selectable tracks.
Users can also generate multiple zoomed in views of the plot at multiple resolutions.
}
\examples{
library(echolocatoR)
finemap_dat<- echolocatoR::BST1; LD_matrix <- echolocatoR::BST1_LD_matrix;
locus_dir <- file.path("~/Desktop","results/GWAS/Nalls23andMe_2019/BST1")

locus_plot <- PLOT.locus(finemap_dat, locus_dir=locus_dir, LD_matrix=LD_matrix, Nott_epigenome=T, xtext=F, plot.zoom=c("5x"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{VALIDATION.compare_bootstrap_distributions}
\alias{VALIDATION.compare_bootstrap_distributions}
\title{Test whether permutation p-value distributions are different}
\usage{
VALIDATION.compare_bootstrap_distributions(
  boot_res,
  formula_str = "stat ~ SNP_group"
)
}
\description{
Test whether permutation p-value distributions are different
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
permute.IMPACT <-  data.table::fread(file.path(root,"_genome_wide/IMPACT/Nalls23andMe_2019.IMPACT.permutations.csv.gz"))
res <- VALIDATION.permute_compare_results(permute_res=permute.IMPACT)
}
}
\seealso{
Other VALIDATION: 
\code{\link{VALIDATION.bootstrap_multimetric}()},
\code{\link{VALIDATION.bootstrap_plot}()},
\code{\link{VALIDATION.bootstrap}()},
\code{\link{VALIDATION.permute}()},
\code{\link{VALIDATION.super_plot}()}
}
\concept{VALIDATION}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COJO.R
\name{COJO.make_locus_subdir}
\alias{COJO.make_locus_subdir}
\title{Make \emph{GCTA-COJO} path}
\source{
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}
}
\usage{
COJO.make_locus_subdir(locus_dir)
}
\description{
Make \emph{GCTA-COJO} path
}
\examples{
\dontrun{
data("locus_dir")
cojo_dir <- COJO.make_locus_subdir(locus_dir=locus_dir)
}
}
\seealso{
Other COJO: 
\code{\link{COJO.conditional}()}
}
\concept{COJO}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fullSS_munged}
\alias{fullSS_munged}
\title{Example subset of full summary stats: munged}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 37033 rows and 11 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
fullSS_munged
}
\description{
Data originally comes from the Parkinson's disease GWAS
by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
}
\details{
Munged using \pkg{MungeSumstats}.
}
\examples{
A subset of the full GWAS summary stats from Nalls et al. (2019)
\dontrun{
data("fullSS_dat")
tmp_file <- tempfile(fileext = ".tsv.gz")
data.table::fwrite(fullSS_dat, tmp_file)
fullSS_munged <- MungeSumstats::format_sumstats(tmp_file, ref_genome = "GRCh37", return_data = TRUE)
usethis::use_data(fullSS_munged, overwrite=TRUE)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.run}
\alias{PAINTOR.run}
\title{Run PAINTOR}
\usage{
PAINTOR.run(
  paintor_path = NULL,
  PT_results_path,
  .LD_file.paths,
  n_datasets,
  method = "enumerate",
  n_causal = 2
)
}
\description{
@keywords internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.snpstats_get_LD}
\alias{LD.snpstats_get_LD}
\title{Get LD using \pkg{snpStats} package}
\source{
\href{https://www.bioconductor.org/packages/release/bioc/html/snpStats.html}{snpStats Bioconductor page}
\href{https://www.bioconductor.org/packages/release/bioc/vignettes/snpStats/inst/doc/ld-vignette.pdf}{LD tutorial}
}
\usage{
LD.snpstats_get_LD(
  LD_folder,
  plink_prefix = "plink",
  select.snps = NULL,
  stats = c("R"),
  symmetric = T,
  depth = "max",
  nThread = 4,
  verbose = T
)
}
\arguments{
\item{LD_folder}{Locus-specific LD output folder.}

\item{stats}{
    A character vector specifying the linkage disequilibrium measures to
    be calculated. This should contain one or more of the strings: \code{"LLR"},
    \code{"OR"}, \code{"Q"}, \code{"Covar"}, \code{"D.prime"},
    \code{"R.squared"}, ad \code{"R"}
  }

\item{symmetric}{
    When no \code{y} argument is supplied this argument controls the format of
    the output band matrices. If \code{TRUE}, symmetric matrices are
    returned and, otherwise, an upper triangular matrices are returned
  }

\item{depth}{
    When \code{y} is not supplied, this parameter is mandatory and
    controls the maximum lag between columns of \code{x}
    considered. Thus, LD statistics are calculated between \code{x[,i]}
    and \code{x[,j]} only if \code{i} and \code{j} differ by no more
    than \code{depth}
  }
}
\description{
Get LD using \pkg{snpStats} package
}
\examples{
subset_DT <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ABCA7/Multi-finemap/ABCA7.Kunkle_2019.1KGphase3_LD.Multi-finemap.tsv.gz")
LD_folder <- "/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ABCA7/LD"
LD_matrix <- LD.snpstats_get_LD(LD_folder=LD_folder, select.snps=subset_DT$SNP)
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT_heatmap}
\alias{IMPACT_heatmap}
\title{Plot \pkg{ComplexHeatmap} of \emph{IMPACT} scores}
\usage{
IMPACT_heatmap(
  ANNOT_MELT,
  snp_groups = c("GWAS lead", "UCS (-PolyFun)", "UCS", "Consensus (-PolyFun)",
    "Consensus")
)
}
\description{
Plot \pkg{ComplexHeatmap} of \emph{IMPACT} scores
}
\examples{
\dontrun{
# merged_dat <- merge_finemapping_results(minimum_support = 1, include_leadSNPs = T,xlsx_path = F,dataset = "Data/GWAS/Nalls23andMe_2019")
 #merged_dat$Locus <- merged_dat$Gene
# ANNOT_MELT <- IMPACT.iterate_get_annotations(merged_dat = merged_dat,
#                                              IMPACT_score_thresh=0,
#                                              baseURL="../../data/IMPACT/IMPACT707/Annotations",
#                                              all_snps_in_range=T,
#                                              top_annotations_only=F)
data.table::fwrite(ANNOT_MELT,"/sc/arion/projects/pd-omics/data/IMPACT/IMPACT707/Annotations/IMPACT_overlap.csv.gz")
ANNOT_MELT <- data.table::fread("~/Desktop/Fine_mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT/IMPACT_overlap.csv.gz")
ANNOT_MELT <- update_cols(subset(ANNOT_MELT, select=-c(SUSIE.Probability,FINEMAP.Probability)))

## Remove no no loci
no_no_loci =  c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1", "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
ANNOT_MELT <- IMPACT.postprocess_annotations(ANNOT_MELT, no_no_loci = no_no_loci)
}
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.help_menu}
\alias{PAINTOR.help_menu}
\title{Show PAINTOR help menu}
\usage{
PAINTOR.help_menu(paintor_path = NULL)
}
\description{
Show PAINTOR help menu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.read_parquet}
\alias{POLYFUN.read_parquet}
\title{Read parquet file as data.frame}
\usage{
POLYFUN.read_parquet(parquet_path, conda_env = "echoR", method = "pandas")
}
\arguments{
\item{method}{Specify whether to use \code{\link{SparkR}} (R)
or \code{pandas} (Python + \code{reticulate}).}
}
\description{
Read parquet file as data.frame
}
\examples{
\dontrun{
parquet_path <- system.file("tools/polyfun/example_data/weights.10.l2.ldscore.parquet", package = "echolocatoR")
# Using python (pandas) - default
dat <- POLYFUN.read_parquet(parquet_path=parquet_path, method="pandas")
# Using R (SparkR)
dat <- POLYFUN.read_parquet(parquet_path=parquet_path, method="sparkR")
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.plot_dataset_overlap}
\alias{SUMMARISE.plot_dataset_overlap}
\title{Plot inter-study SNP overlap}
\usage{
SUMMARISE.plot_dataset_overlap(
  merged_dat,
  snp_filter = "!is.na(SNP)",
  filename = NA,
  formula_str = "~ SNP + Dataset",
  triangle = F,
  proxies = NULL
)
}
\description{
Cross-tabulate SNP overlap (after applying filter)
between each pair of studies.
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{results_report}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{locus_dir}
\alias{locus_dir}
\title{Example results path for BST1 locus}
\format{
An object of class \code{character} of length 1.
}
\usage{
locus_dir
}
\description{
Example results path for BST1 locus
}
\examples{
\dontrun{
locus_dir <- "results/GWAS/Nalls23andMe_2019/BST1"
usethis::use_data(locus_dir, overwrite=T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Roadmap.R
\name{ROADMAP.query_and_plot}
\alias{ROADMAP.query_and_plot}
\title{Query and plot Roadmap epigenomic annotations}
\usage{
ROADMAP.query_and_plot(
  subset_DT,
  results_path = "./ROADMAP",
  n_top_tissues = NULL,
  keyword_query = NULL,
  adjust = 0.2,
  force_new_query = F,
  remove_tmps = T
)
}
\arguments{
\item{subset_DT}{Data.frame with at least the following columns:
\describe{
\item{SNP}{SNP RSID}
\item{CHR}{chromosome}
\item{POS}{position}
}}

\item{results_path}{Where to store query results.}

\item{n_top_tissues}{The number of top tissues to include,
sorted by greatest number of rows
(i.e. the number of genomic ranges within the window).}

\item{keyword_query}{Search all columns in the Roadmap annotations metadata
and only query annotations that contain your keywords.
Can provide multiple keywords in list form:
\code{c("placenta","liver","monocytes")}}

\item{force_new_query}{Force a new query from the XGR database.}
}
\value{
A named list containing:
\itemize{
\item{\code{ggbio} plot}
\item{\code{GRanges} object within the queried coordinates}
}
}
\description{
Query and plot Roadmap epigenomic annotations
}
\examples{
\dontrun{
data("BST1")
finemap_DT <- BST1
roadmap_plot_query <- ROADMAP.query_and_plot(subset_DT=finemap_DT, keyword_query="monocytes")
}
}
\seealso{
Other ROADMAP: 
\code{\link{ROADMAP.construct_reference}()},
\code{\link{ROADMAP.merge_and_process_grl}()},
\code{\link{ROADMAP.query}()},
\code{\link{ROADMAP.tabix}()}
}
\concept{ROADMAP}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_finemap.R
\name{check_necessary_cols}
\alias{check_necessary_cols}
\title{Check for necessary columns}
\usage{
check_necessary_cols(
  subset_DT,
  finemap_methods,
  sample_size = NULL,
  dataset_type = "GWAS",
  verbose = T
)
}
\description{
Check for necessary columns
}
\examples{
data("BST1");
finemap_methods <- c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE")
finemap_methods <- check_necessary_cols(subset_DT=BST1, finemap_methods=finemap_methods)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.compute_priors}
\alias{POLYFUN.compute_priors}
\title{Recompute SNP-wise priors from summary stats}
\usage{
POLYFUN.compute_priors(
  polyfun = NULL,
  locus_dir,
  fullSS_path,
  sample_size = NULL,
  min_INFO = 0,
  min_MAF = 0.001,
  annotations_path = NULL,
  weights_path = NULL,
  prefix = "PD_GWAS",
  chrom = "all",
  compute_ldscores = F,
  allow_missing_SNPs = T,
  ref_prefix = NULL,
  remove_tmps = T,
  conda_env = "echoR"
)
}
\description{
Recompute SNP-wise priors from summary stats
}
\examples{
\dontrun{
data("locus_dir");
fullSS_path <- example_fullSS(fullSS_path="~/Desktop/Nalls23andMe_2019.fullSS_subset.tsv")
munged_path <- "results/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/nallsEtAl2019_allSamples_allVariants.mod.munged.parquet"
LDSC.files <- POLYFUN.compute_priors(locus_dir=locus_dir, fullSS_path=fullSS_path, conda_env="echoR")
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LRRK2}
\alias{LRRK2}
\title{\pkg{echolocatoR} output example: LRRK2 locus}
\format{
data.table
\describe{
  \item{SNP}{SNP RSID}
  \item{CHR}{Chromosome}
  \item{POS}{Genomic position (in basepairs)}
  \item{...}{Optional: extra columns}
}
}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
LRRK2
}
\description{
An example results file after running
\code{\link[=finemap_loci]{finemap_loci()}} on the \emph{LRRK2} locus.
}
\details{
Data originally comes from the Parkinson's disease GWAS
by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
}
\examples{
\dontrun{
root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap"
LRRK2 <- data.table::fread(file.path(root_dir,"Multi-finemap_results.txt"))
LRRK2 <- update_cols(finemap_dat=LRRK2)
LRRK2 <- find_consensus_SNPs(finemap_dat=LRRK2)
usethis::use_data(LRRK2, overwrite = T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.R
\name{auto_topSNPs_sub}
\alias{auto_topSNPs_sub}
\title{Automatically identify top SNP per locus}
\usage{
auto_topSNPs_sub(top_SNPs, query, locus)
}
\description{
If no \code{top_SNPs} dataframe is supplied,
 this function will sort by p-value and then effect size,
 and use the SNP in the first row.
}
\seealso{
Other standardization functions: 
\code{\link{calculate_tstat}()}
}
\concept{standardization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.get_precomputed_priors}
\alias{POLYFUN.get_precomputed_priors}
\title{Extract pre-computed prior}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.get_precomputed_priors(
  polyfun = NULL,
  locus_dir,
  finemap_dat = NULL,
  force_new_priors = T,
  remove_tmps = F,
  nThread = 4,
  conda_env = "echoR"
)
}
\description{
Extract SNP-wise prior probabilities pre-computed from many UK Biobank traits.
}
\examples{
\dontrun{

data("BST1"); data("locus_dir");
priors <- POLYFUN.get_precomputed_priors(locus_dir=locus_dir, finemap_dat=BST1)
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.get_epigenomic_peaks}
\alias{NOTT_2019.get_epigenomic_peaks}
\title{Download cell type-specific epigenomic peaks}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.get_epigenomic_peaks(
  assays = c("ATAC", "H3K27ac", "H3K4me3"),
  cell_types = c("neurons", "microglia", "oligo", "astrocytes"),
  convert_to_GRanges = T,
  nThread = 4,
  verbose = T
)
}
\description{
API access to brain cell type-specific epigenomic peaks (bed format)
from Nott et al. (2019).
}
\examples{
PEAKS <- NOTT_2019.get_epigenomic_peaks(nThread=1)
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finemapping_portal_API.R
\name{GITHUB.portal_metadata}
\alias{GITHUB.portal_metadata}
\title{echolocatoR Fine-mapping portal: metadata}
\usage{
GITHUB.portal_metadata(verbose = T)
}
\description{
echolocatoR Fine-mapping portal: metadata
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.plot}
\alias{POLYFUN.plot}
\title{Plot PolyFun and other fine-mapping results}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.plot(
  subset_DT,
  LD_matrix,
  locus = NULL,
  subtitle = "PolyFun Comparison",
  conditions = c("SUSIE", "POLYFUN_SUSIE", "FINEMAP", "PAINTOR", "PAINTOR_Fairfax")
)
}
\description{
Plot PolyFun and other fine-mapping results
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.filter_vcf_gaston}
\alias{LD.filter_vcf_gaston}
\title{Subset a vcf by superpopulation}
\usage{
LD.filter_vcf_gaston(
  vcf_subset,
  subset_DT,
  locus_dir,
  superpopulation,
  popDat,
  verbose = T
)
}
\arguments{
\item{vcf_subset}{Path to the locus subset vcf.}

\item{superpopulation}{Subset your LD reference panel by superopulation.
Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
\href{https://www.internationalgenome.org/faq/which-populations-are-part-your-study/}{1KGphase1 options} include:
\describe{
\item{"AFR"}{African [descent]}
\item{"AMR"}{Ad-mixed American}
\item{"EAS"}{East Asian}
\item{"EUR"}{European}
\item{"SAS"}{South Asian}
}}

\item{popDat}{The metadata file listing the superpopulation
to which each sample belongs.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Subset a vcf by superpopulation
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.epigenomic_histograms}
\alias{NOTT_2019.epigenomic_histograms}
\title{Plot brain cell-specific epigenomic data}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.epigenomic_histograms(
  finemap_dat,
  locus_dir,
  show_plot = T,
  save_plot = T,
  full_data = T,
  return_assay_track = F,
  binwidth = 200,
  density_adjust = 0.2,
  plot.zoom = "1x",
  strip.text.y.angle = 90,
  xtext = T,
  geom = "density",
  plot_formula = "Cell_type ~.",
  fill_var = "Assay",
  bigwig_dir = NULL,
  genomic_units = "Mb",
  as_ggplot = T,
  nThread = 4,
  save_annot = F,
  verbose = T
)
}
\arguments{
\item{plot.zoom}{Zoom into the center of the locus when plotting (without editing the fine-mapping results file).
You can provide either:
\itemize{
\item{The size of your plot window in terms of basepairs (e.g. \code{plot.zoom=50000} for a 50kb window)}.
\item{How much you want to zoom in (e.g. \code{plot.zoom="1x"} for the full locus, \code{plot.zoom="2x"} for 2x zoom into the center of the locus, etc.)}.
}
You can pass a list of window sizes (e.g. \code{c(50000,100000,500000)}) to automatically generate
multiple views of each locus.
This can even be a mix of different style inputs: e.g. \code{c("1x","4.5x",25000)}.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Brain cell-specific epigenomic data from Nott et al. (2019).
}
\examples{
data("BST1"); data("locus_dir");
track.Nott_histo <- NOTT_2019.epigenomic_histograms(finemap_dat = BST1, locus_dir = locus_dir, save_plot=F, return_assay_track=T, save_annot=F)
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.iterate_enrichment}
\alias{IMPACT.iterate_enrichment}
\title{Iterate IMPACT enrichment tests}
\usage{
IMPACT.iterate_enrichment(
  gwas_paths,
  annot_baseURL = "../../data/IMPACT/IMPACT707/Annotations"
)
}
\description{
Iterate IMPACT enrichment tests
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRanges_overlap.R
\name{GRanges_overlap}
\alias{GRanges_overlap}
\title{Find overlap between genomic coordinates/ranges}
\usage{
GRanges_overlap(
  dat1,
  dat2,
  chrom_col.1 = "chrom",
  start_col.1 = "start",
  end_col.1 = "end",
  chrom_col.2 = "chrom",
  start_col.2 = "start",
  end_col.2 = "end",
  return_merged = T,
  chr_format = "NCBI",
  verbose = F
)
}
\description{
Find overlap between genomic coordinates/ranges
}
\concept{GRanges}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{VALIDATION.bootstrap_plot}
\alias{VALIDATION.bootstrap_plot}
\title{Plot permutation test results}
\usage{
VALIDATION.bootstrap_plot(
  boot_res,
  validation_method = NULL,
  facet_formula = ". ~ .",
  y_var = "z",
  override_metric_count = F,
  remove_random = T,
  show_plot = T,
  save_plot = F,
  box.padding = 0.5,
  font_size = 3,
  height = 5,
  width = 7,
  verbose = T
)
}
\description{
Plot permutation test results
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"

#### h2 ####
validation_method <- "S-LDSC heritability"
## mean
path <- file.path(root,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.bootstrap.coin_wilcox_test.csv.gz")
## raw
path <-  file.path(root,"PolyFun/Nalls23andMe_2019.h2.bootstrap.coin_wilcox_test.csv.gz")

#### IMPACT ####
validation_method <- "IMPACT"
### mean
path <- file.path(root,"IMPACT/TOP_IMPACT_all.bootstrap.coin_wilcox_test.csv.gz")
### raw
path <- file.path(root,"IMPACT/Nalls23andMe_2019.IMPACT_score.permutations.csv.gz")

#### SURE MPRA ####
validation_method <- "SuRE MPRA"
## mean
path <- file.path(root,"SURE/Nalls23andMe_2019.SURE.snp_groups.mean.bootstrap.stats_wilcox.test.csv.gz")
### raw
path <- file.path(root,"SURE/Nalls23andMe_2019.SURE.bootstrap.stats_wilcox.test.csv.gz")


#### Dey_DeepLearning ####
validation_method <- "Dey_DeepLearning"
## mean
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.snp_groups_mean.bootstrap.coin_wilcox_test_subset.csv.gz")
## raw
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.bootstrap.stats_wilcox.test.csv.gz")
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.bootstrap.stats_wilcox.test.csv.gz")

# ------  Import
boot_res <-  data.table::fread(path, nThread=4)
if(validation_method=="SuRE MPRA") boot_res$z <- -boot_res$z

# ---Plot ---
library(patchwork)
if(validation_method=="Dey_DeepLearning") boot_res <- data.frame(boot_res)[grepl("H3K4ME3", boot_res$Metric, fixed=T) & grepl("Basenji", boot_res$Metric, fixed=T) & grepl("_MAX_", boot_res$Metric, fixed=T),]
facet_formula <- if(validation_method=="Dey_DeepLearning") "Assay  ~ Model" else ".~."
facet_formula <- if(validation_method=="SuRE MPRA") ".  ~ Cell_type" else ".~."
gp1 <- VALIDATION.bootstrap_plot(boot_res=boot_res, validation_method=validation_method, y_var="z", save_plot=gsub("\\\\.csv\\\\.gz",".stat_values.png",path), width=9, facet_formula=facet_formula, override_metric_count = T)
gp2 <- VALIDATION.bootstrap_plot(boot_res=subset(boot_res, stat>0), validation_method=validation_method, y_var="p", save_plot=gsub("\\\\.csv\\\\.gz",".p_values.png",path), width=9,facet_formula = facet_formula, override_metric_count = T)
gp12 <- (gp1 / gp2) + patchwork::plot_annotation(tag_levels = "a")
ggsave(gsub("\\\\.csv\\\\.gz",".png",path),gp12, dpi=400, height=9, width=15)

}
}
\seealso{
Other VALIDATION: 
\code{\link{VALIDATION.bootstrap_multimetric}()},
\code{\link{VALIDATION.bootstrap}()},
\code{\link{VALIDATION.compare_bootstrap_distributions}()},
\code{\link{VALIDATION.permute}()},
\code{\link{VALIDATION.super_plot}()}
}
\concept{VALIDATION}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{super_summary_plot}
\alias{super_summary_plot}
\title{Merge all summary plots into one super plot}
\usage{
super_summary_plot(
  merged_dat,
  snp_filter = "Consensus_SNP==T",
  coloc_results = NULL,
  plot_missense = T,
  show_plot = T,
  save_plot = F,
  height = 15,
  width = 13,
  dpi = 500
)
}
\description{
Merge all summary plots into one super plot
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{results_report}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabix.R
\name{TABIX.query}
\alias{TABIX.query}
\title{Query a tabix file}
\usage{
TABIX.query(fullSS_tabix, chrom, start_pos, end_pos, verbose = TRUE)
}
\description{
Query by genomic coordinates.
}
\examples{
\dontrun{
fullSS_path <- example_fullSS()
fullSS_tabix <- TABIX.convert_file(fullSS_path=fullSS_path, chrom_col="CHR", position_col="BP")
TABIX.query(fullSS_tabix=fullSS_tabix, chrom=4, start_pos=13737637, end_pos=13837637)
}
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_multifinemap_path}
\alias{get_multifinemap_path}
\title{Create multi-finemap path}
\usage{
get_multifinemap_path()
}
\description{
Create multi-finemap path
}
\seealso{
Other directory functions: 
\code{\link{get_locus_dir}()},
\code{\link{get_study_dir}()},
\code{\link{get_subset_path}()},
\code{\link{make_locus_dir}()}
}
\concept{directory functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{VALIDATION.compare_support}
\alias{VALIDATION.compare_support}
\title{Check whether Support level correlates with some variable}
\usage{
VALIDATION.compare_support(plot_dat)
}
\description{
Check whether Support level correlates with some variable
}
\examples{
\dontrun{
# S-LDSC h2
h2_stats <- data.table::fread("/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/Nalls23andMe_2019.grouped_snpvar_stats.csv.gz")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Roadmap.R
\name{ROADMAP.merge_and_process_grl}
\alias{ROADMAP.merge_and_process_grl}
\title{Standardize Roadmap query}
\usage{
ROADMAP.merge_and_process_grl(
  grl.roadmap,
  gr.snp,
  n_top_tissues = 5,
  sep = " "
)
}
\arguments{
\item{grl.roadmap}{Roadmap query results}

\item{n_top_tissues}{The number of top tissues to include,
sorted by greatest number of rows
(i.e. the number of genomic ranges within the window).}
}
\description{
Standardize Roadmap query
}
\seealso{
Other ROADMAP: 
\code{\link{ROADMAP.construct_reference}()},
\code{\link{ROADMAP.query_and_plot}()},
\code{\link{ROADMAP.query}()},
\code{\link{ROADMAP.tabix}()}
}
\concept{ROADMAP}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{haploR.regulomeDB}
\alias{haploR.regulomeDB}
\title{Download SNP-wise annotations from RegulomeDB}
\source{
\href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
}
\usage{
haploR.regulomeDB(snp_list, verbose = T, chunk_size = NA)
}
\description{
Download SNP-wise annotations from RegulomeDB
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.get_annotation_key}
\alias{IMPACT.get_annotation_key}
\title{Get \emph{IMPACT} annotation key}
\usage{
IMPACT.get_annotation_key(
 
    URL = "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/IMPACT_annotation_key.txt",
  save_path = "./IMPACT/IMPACT_annotation_key.txt.gz",
  save_key = F,
  force_new_download = F,
  verbose = T
)
}
\description{
Inlcludes the study source, tissue,cell type/line,
and cell subtype of each of the 500 annotations in \emph{IMPACT}.
}
\examples{
annot.key <- IMPACT.get_annotation_key(save_key=F)
head(annot.key)
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{IMPACT_heatmap}()},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SURE.R
\name{SURE.melt_snp_groups}
\alias{SURE.melt_snp_groups}
\title{Aggregate \emph{SuRE} p-values for each SNP group}
\source{
\href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
\href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
\href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
\href{https://sure.nki.nl}{SNP-SuRE data browse}
}
\usage{
SURE.melt_snp_groups(
  sure_DT,
  grouping_vars = c("Locus"),
  metric_str = "mean",
  replace_NA = NA,
  save_path = F,
  verbose = T
)
}
\description{
Aggregate \emph{SuRE} p-values for each SNP group
}
\examples{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"

sure_DT <- data.table::fread(file.path(root,"SURE/Nalls23andMe_2019.SURE.csv.gz"), nThread=4)
sure_DT <- find_consensus_SNPs_no_PolyFun(sure_DT)

sure.melt <- SURE.melt_snp_groups(sure_DT, save_path=file.path(root,"SURE/Nalls23andMe_2019.SURE.snp_groups.mean.csv.gz"))
}
\seealso{
Other SURE: 
\code{\link{SURE.merge}()},
\code{\link{SURE.plot}()}
}
\concept{SURE}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COJO.R
\name{COJO.conditional}
\alias{COJO.conditional}
\title{Run \emph{GCTA-COJO} conditional analysis}
\source{
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}
}
\usage{
COJO.conditional(
  GCTA_path = system.file("tools/gcta_1.92.1beta5_mac/bin", "gcta64", package =
    "echolocatoR"),
  locus_dir,
  conditioned_snps,
  min_MAF = 0,
  diff_freq = 0.1,
  excluded_path
)
}
\description{
Condition on a given SNP (or list of SNPs)
and calculate residual effects of all other SNPs.
}
\seealso{
Other COJO: 
\code{\link{COJO.make_locus_subdir}()}
}
\concept{COJO}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/column_dictionary.R
\name{column_dictionary}
\alias{column_dictionary}
\title{Map column names to positions.}
\usage{
column_dictionary(file_path)
}
\arguments{
\item{file_path}{Path to full summary stats file
(or any really file you want to make a column dictionary for).}
}
\value{
Named list of column positions.
}
\description{
Useful in situations where you need to specify columns by index instead of name (e.g. awk queries).
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_finemap.R
\name{multi_finemap}
\alias{multi_finemap}
\title{Fine-map using multiple fine-mapping tools}
\usage{
multi_finemap(
  locus_dir,
  fullSS_path,
  finemap_method_list,
  finemap_args = NULL,
  subset_DT,
  dataset_type,
  LD_matrix = NULL,
  n_causal = 5,
  sample_size = NULL,
  PAINTOR_QTL_datasets = NULL,
  PP_threshold = 0.95,
  case_control = T,
  verbose = T,
  nThread = 4,
  conda_env = "echoR"
)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{dataset_type}{The kind dataset you're fine-mapping (e.g. GWAS, eQTL, tQTL).
This will also be used when creating the subdirectory where your results will be stored
(e.g. \emph{Data/<dataset_type>/Kunkle_2019}).}

\item{n_causal}{The maximum number of potential causal SNPs per locus.
This parameter is used somewhat differntly by different fine-mapping tools.
See tool-specific functions for details.}

\item{sample_size}{The overall sample size of the study.
If none is given, and \strong{N_cases} and \strong{N_controls} columns are present,
then sample_size is inferred to be:  \code{max(N_cases) + max(N_controls)}.}

\item{PAINTOR_QTL_datasets}{A list of QTL datasets to be used when conducting joint functional fine-mapping with \emph{PAINTOR}.}

\item{PP_threshold}{The minimum fine-mapped posterior probability for a SNP to be considered part of a Credible Set.
For example, \code{PP_threshold=.95} means that all Credible Set SNPs will be 95\% Credible Set SNPs.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}

\item{conda_env}{The name of a conda environment to use.}
}
\description{
Fine-mapping will be repeated on the same locus using each of the tools in \code{finemap_method_list}.
Then, all results will be merged into the locus-specific multi-finemap file,
along with the original per-SNP GWAS/QTL summary statistics.
Each tools will have the following columns:
\describe{
 \item{<tool>.PP}{The posterior probability (PP) of a SNP being causal for the trait. Though this is a generalization and the exact meaning of PP will differ by tools (e.g. Posterior Inclusion Probability for SUSIE).}
 \item{<tool>.CS}{Which credible set the SNP is part of (within a locus). If \code{=0}, then the SNP was not part of any credible set. Some tools only produce one credible set per locus.}
}
}
\examples{
\dontrun{
data("BST1"); data("BST1_LD_matrix");
subset_DT <- BST1
finemap_method_list <- c("ABF","SUSIE")
}
}
\seealso{
Other finemapping functions: 
\code{\link{find_consensus_SNPs}()}
}
\concept{finemapping functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.run_ldsc}
\alias{POLYFUN.run_ldsc}
\title{Run a modified version of S-LDSC}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.run_ldsc(
  polyfun = NULL,
  output_dir = NULL,
  munged.path,
  min_INFO = 0.6,
  min_MAF = 0.05,
  annotations.path = file.path(polyfun, "example_data/annotations."),
  weights.path = file.path(polyfun, "example_data/weights."),
  prefix = "LDSC",
  chrom = "all",
  compute_ldscores = F,
  allow_missing_SNPs = T,
 
    munged_path = "/sc/arion/projects/pd-omics/tools/polyfun/Nalls23andMe_2019.sumstats_munged.parquet",
  ref.prefix = "/sc/arion/projects/pd-omics/data/1000_Genomes/Phase1/1000G.mac5eur.",
  freq.prefix = "/sc/arion/projects/pd-omics/tools/polyfun/1000G_frq/1000G.mac5eur.",
  conda_env = "echoR"
)
}
\description{
Modifications to S-LDSC include L2-regularization.
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_cleaning.R
\name{vcf_cleaning}
\alias{vcf_cleaning}
\title{Remove old vcf files}
\usage{
vcf_cleaning(root, LD_ref = "1KGphase3", loci = NULL, verbose = T)
}
\description{
Remove old vcf files
}
\examples{
vcf_cleaning(root="Data/GWAS/Ripke_2014", LD_ref="1KGphase3", loci=top_SNPs$Locus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GoShifter.R
\name{GOSHIFTER.create_snpmap}
\alias{GOSHIFTER.create_snpmap}
\title{Prepare SNP input for \emph{GoShifter}}
\usage{
GOSHIFTER.create_snpmap(finemap_dat, GS_results, nThread = 4, verbose = T)
}
\description{
Prepare SNP input for \emph{GoShifter}
}
\examples{
data("BST1"); data("locus_dir")
snpmap <- GOSHIFTER.create_snpmap(finemap_dat=BST1, GS_results=file.path(locus_dir,"GoShifter"))
}
\seealso{
Other GOSHIFTER: 
\code{\link{GOSHIFTER.create_LD}()}
}
\concept{GOSHIFTER}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.plac_seq_plot}
\alias{NOTT_2019.plac_seq_plot}
\title{Plot brain cell-specific interactome data}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.plac_seq_plot(
  finemap_dat = NULL,
  locus_dir = NULL,
  title = NULL,
  print_plot = T,
  save_plot = T,
  return_interaction_track = F,
  xlims = NULL,
  zoom_window = NULL,
  index_SNP = NULL,
  genomic_units = "Mb",
  color_dict = c(enhancers = "springgreen2", promoters = "purple", anchors = "black"),
  return_consensus_overlap = T,
  show_arches = T,
  highlight_plac = F,
  show_regulatory_rects = T,
  show_anchors = T,
  strip.text.y.angle = 0,
  xtext = T,
  save_annot = F,
  point_size = 2,
  height = 7,
  width = 7,
  dpi = 300,
  as_ggplot = T,
  nThread = 4,
  verbose = T
)
}
\description{
Plot brain cell-specific interactome data
}
\examples{
\dontrun{
data("BST1"); data("locus_dir");

trks_plus_lines <- NOTT_2019.plac_seq_plot(finemap_dat=BST1, locus_dir=file.path("~/Desktop",locus_dir), highlight_plac=T)
# Zoom in
trks_plus_lines <- NOTT_2019.plac_seq_plot(finemap_dat=BST1, locus_dir=file.path("~/Desktop",locus_dir), zoom_window=500000, highlight_plac=T)
}
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORCES_2020.R
\name{CORCES_2020.prepare_scATAC_peak_overlap}
\alias{CORCES_2020.prepare_scATAC_peak_overlap}
\title{Prepare data to plot overlap between datatable of SNPs and
cell-type-specific epigenomic peaks and coaccessibility data.}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.prepare_scATAC_peak_overlap(
  merged_dat,
  FDR_filter = NULL,
  snp_filter = "Consensus_SNP==T",
  add_cicero = T,
  annotate_genes = T,
  return_counts = T,
  verbose = T
)
}
\description{
Prepare data to plot overlap between datatable of SNPs and
cell-type-specific epigenomic peaks and coaccessibility data.
}
\examples{
data("merged_DT");
merged_dat <- subset(merged_DT, Consensus_SNP)
dat_melt <- CORCES_2020.prepare_scATAC_peak_overlap(merged_dat=merged_dat)
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.iterate_get_annotations}
\alias{IMPACT.iterate_get_annotations}
\title{Gather all \emph{IMPACT} annotations that overlap with your query}
\usage{
IMPACT.iterate_get_annotations(
  merged_dat,
  IMPACT_score_thresh = 0.1,
  baseURL = "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/Annotations",
  all_snps_in_range = T,
  top_annotations = F,
  force_one_annot_per_locus = F,
  snp.filter = "!is.na(SNP)",
  nThread = 4,
  verbose = T
)
}
\description{
Iterates over each unique locus.
\bold{\emph{WARNING!}} These files are quite large and you need to make sure you
have ample memory and storage on your computer for best results.
}
\details{
Unfortunately, you have to download the entire chromosome file at once,
 because they aren't Tabix indexed. To minimize the memory load,
 this function only keeps the portion of the \emph{IMPACT} file that overlaps with the
 coordinates in \code{subset_DT}.
}
\examples{
\dontrun{
data("merged_DT")
ANNOT_MELT <- IMPACT.iterate_get_annotations(merged_dat=merged_DT)
}
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{IMPACT_heatmap}()},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.index_vcf}
\alias{LD.index_vcf}
\title{Index vcf file if it hasn't been already}
\usage{
LD.index_vcf(vcf_file, force_new_index = F, conda_env = "echoR", verbose = T)
}
\description{
Index vcf file if it hasn't been already
}
\examples{
\dontrun{
LD_reference <- "~/Desktop/results/Reference/custom_panel_chr4.vcf.gz"
vcf_file <- LD.index_vcf(vcf_file=LD_reference)
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.CS_counts_plot}
\alias{SUMMARISE.CS_counts_plot}
\title{Bar plot of tool-specific CS sizes}
\usage{
SUMMARISE.CS_counts_plot(
  merged_dat,
  show_numbers = T,
  ylabel = "Locus",
  legend_nrow = 3,
  label_yaxis = T,
  top_CS_only = F,
  show_plot = T
)
}
\description{
Loci ordered by UCS size (smallest to largest).
}
\examples{
data("merged_DT")
gg_CS <- SUMMARISE.CS_counts_plot(merged_dat=merged_DT)
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{results_report}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.gather_annotations}
\alias{POLYFUN.gather_annotations}
\title{Find and import PolyFun annotation files}
\usage{
POLYFUN.gather_annotations(
  chromosomes = c(1:22),
  subset_snps = NULL,
  polyfun_annots
)
}
\description{
Find and import PolyFun annotation files
}
\examples{
\dontrun{
data("BST1")
finemap_DT <- BST1
subset_snps <- finemap_DT$SNP
annot_DT <- POLYFUN.gather_annotations(chromosomes=finemap_DT$CHR[1], subset_snps=subset_snps, polyfun_annots="/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB")
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.fill_NA}
\alias{LD.fill_NA}
\title{Fill NAs in an LD matrix}
\usage{
LD.fill_NA(LD_matrix, fillNA = 0, verbose = F)
}
\description{
Trickier than it looks.
}
\examples{
\dontrun{
data("LD_matrix");
LD_matrix <- LD.fill_NA(LD_matrix)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{popDat_1KGphase1}
\alias{popDat_1KGphase1}
\title{Population metadata: 1KGphase1}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 1091 rows and 4 columns.
}
\usage{
popDat_1KGphase1
}
\description{
Population metadata: 1KGphase1
}
\examples{
\dontrun{
popDat_URL <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
popDat_1KGphase1 <-  data.table::fread(text=trimws(gsub(",\t",",",readLines(popDat_URL))), sep="\t",  fill=T, stringsAsFactors = F, col.names = c("sample","population","superpop","sex"), nThread = 4)
usethis::use_data(popDat_1KGphase1, overwrite = T)
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createDT_html.R
\name{createDT_html}
\alias{createDT_html}
\title{Interactive DT (html)}
\usage{
createDT_html(DF, caption = "", scrollY = 400)
}
\description{
Generate an interactive data table with download buttons.
Use this function when manually constructing rmarkdown chunks using cat() in a for loop.
}
\seealso{
Other general: 
\code{\link{createDT}()},
\code{\link{dt.replace}()},
\code{\link{example_fullSS}()},
\code{\link{fillNA_CS_PP}()},
\code{\link{get_os}()},
\code{\link{get_sample_size}()},
\code{\link{tryFunc}()}
}
\concept{general}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{merged_DT}
\alias{merged_DT}
\title{\emph{echolocatoR} output example: all loci}
\format{
data.table
\describe{
  \item{SNP}{SNP RSID}
  \item{CHR}{Chromosome}
  \item{POS}{Genomic position (in basepairs)}
  \item{...}{Optional: extra columns}
}
}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
merged_DT
}
\description{
An example results file after running \code{finemap_loci}
 on all Parkinson's disease (PD)-associated loci.
}
\details{
Data originally comes from the PD GWAS
by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
}
\examples{
\dontrun{
merged_DT <- merge_finemapping_results(dataset = "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019", minimum_support=1)
usethis::use_data(merged_DT, overwrite=T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpliceAI.R
\name{SPLICEAI.subset_precomputed_tsv}
\alias{SPLICEAI.subset_precomputed_tsv}
\title{Query genome-wide \emph{SpliceAI} results file (tsv format)}
\source{
\href{https://github.com/Illumina/SpliceAI}{GitHub}
\href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
}
\usage{
SPLICEAI.subset_precomputed_tsv(
  subset_DT,
 
    precomputed_path = "/pd-omics/data/spliceAI/whole_genome_filtered_spliceai_scores.tsv.gz",
  merge_data = T,
  drop_na = T,
  filtered = T
)
}
\description{
Query genome-wide \emph{SpliceAI} results file (tsv format)
}
\seealso{
Other SpliceAI: 
\code{\link{SPLICEAI.plot}()},
\code{\link{SPLICEAI.run}()},
\code{\link{SPLICEAI.snp_probs}()},
\code{\link{SPLICEAI.subset_precomputed_tsv_iterate}()},
\code{\link{SPLICEAI.subset_precomputed_vcf}()}
}
\concept{SpliceAI}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dey_DeepLearning.R
\name{DEEPLEARNING.plot}
\alias{DEEPLEARNING.plot}
\title{Plot deep learning predictions}
\usage{
DEEPLEARNING.plot(
  annot.melt,
  snp_groups = c("GWAS lead", "UCS", "Consensus (-POLYFUN)", "Consensus"),
  comparisons_filter = function(x) {     if ("Consensus" \%in\% x)          return(x) },
  model_metric = c("MAX"),
  facet_formula = ". ~ Model",
  remove_outliers = T,
  show_padj = T,
  show_signif = T,
  vjust_signif = 0,
  show_plot = T,
  save_path = F,
  height = 6,
  width = 8
)
}
\description{
Plot deep learning predictions
}
\examples{
\dontrun{
root = "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"

#### Allelic_Effect ####
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.snp_groups_mean.csv.gz")

#### Variant_Level ####
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.snp_groups_mean.csv.gz")

annot.melt <-  data.table::fread(path, nThread=8)

gp <- DEEPLEARNING.plot(annot.melt=annot.melt, facet_formula="Tissue ~ Model", comparisons_filter=NULL, save_path=gsub("\\\\.csv\\\\.gz",".png",path))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GoShifter.R
\name{GOSHIFTER.run}
\alias{GOSHIFTER.run}
\title{Run \code{GoShifter} enrichment test}
\source{
\href{https://github.com/immunogenomics/goshifter}{GitHub}
\href{https://pubmed.ncbi.nlm.nih.gov/26140449/}{Publication}
}
\usage{
GOSHIFTER.run(
  finemap_dat,
  locus_dir,
  GRlist,
  permutations = 1000,
  goshifter_path = NULL,
  chromatin_state = "TssA",
  R2_filter = 0.8,
  overlap_threshold = 1,
  remove_tmps = T,
  verbose = T
)
}
\description{
\code{GoShifter}: Locus-specific enrichment tool.
}
\examples{
data("BST1"); data("locus_dir")
peaks <- NOTT_2019.get_epigenomic_peaks(convert_to_GRanges = T)
grl.peaks <- GenomicRanges::makeGRangesListFromDataFrame(data.frame(peaks), split.field ="Cell_type", names.field ="seqnames" )
GS_results <- GOSHIFTER.run(finemap_dat=subset(BST1, P<5e-8), locus_dir=locus_dir, GRlist=grl.peaks)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.plot_enrichment}
\alias{XGR.plot_enrichment}
\title{Plot XGR enrichment}
\usage{
XGR.plot_enrichment(
  enrich_res,
  adjp_thresh = 0.05,
  top_annotations = NULL,
  show_plot = T
)
}
\description{
Plot XGR enrichment
}
\examples{
\dontrun{
data("merged_DT")
enrich_res <- XGR.iterate_enrichment(subset_DT=merged_DT, foreground_filter = "Consensus_SNP", background_filter = "leadSNP", lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes"), nCores=1)
XGR.plot_enrichment(enrich_res)
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{merge_finemapping_results_each}
\alias{merge_finemapping_results_each}
\title{Create full cross-locus merged files for each dataset,
then return a subset of those files as one super-merged table.}
\usage{
merge_finemapping_results_each(
  study_dirs,
  LD_reference = "1KGphase3",
  minimum_support = 1,
  include_leadSNPs = T,
  return_filter = "!is.na(SNP)",
  merged_path = "merged_dat.csv.gz",
  force_new_merge = F,
  nThread = 4,
  verbose = T
)
}
\description{
Create full cross-locus merged files for each dataset,
then return a subset of those files as one super-merged table.
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()}
}
\concept{annotate}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.env_from_yaml}
\alias{CONDA.env_from_yaml}
\title{Create conda env from yaml file}
\usage{
CONDA.env_from_yaml(
  yaml_path = system.file("conda", "echoR.yml", package = "echolocatoR")
)
}
\description{
Create conda env from yaml file
}
\seealso{
Other conda: 
\code{\link{CONDA.activate_env}()},
\code{\link{CONDA.create_echoR_env}()},
\code{\link{CONDA.find_env_Rlib}()},
\code{\link{CONDA.find_python_path}()},
\code{\link{CONDA.install}()}
}
\concept{conda}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/messager.R
\name{messager}
\alias{messager}
\title{Print messages}
\usage{
messager(..., v = TRUE, parallel = FALSE)
}
\arguments{
\item{v}{Whether to print messages or not.}

\item{parallel}{Whether to enable message print when wrapped 
in parallelised functions.}
}
\value{
Null
}
\description{
Conditionally print messages.
 Allows developers to easily control verbosity of functions, 
 and meet Bioconductor requirements that dictate the message 
 must first be stored to a variable before passing to \link[base]{message}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zscore.R
\name{zscore}
\alias{zscore}
\title{Compute Z-score}
\usage{
zscore(vec)
}
\description{
Computes Z-score when you have the full vector of values (not just a subset).
}
\details{
These functions are necessary for \code{\link{PAINTOR}}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_cols.R
\name{update_cols}
\alias{update_cols}
\title{Update CS cols}
\usage{
update_cols(finemap_dat)
}
\description{
Convert old col format: ".Credible_Set" => ".CS"
}
\examples{
data("BST1");
finemap_DT <- update_cols(finemap_dat=BST1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BST1_LD_matrix}
\alias{BST1_LD_matrix}
\title{LD with the lead SNP: BST1 locus}
\format{
data.table
\describe{
  \item{SNP}{SNP RSID}
  \item{CHR}{Chromosome}
  \item{POS}{Genomic position (in basepairs)}
  \item{...}{Optional: extra columns}
}
}
\source{
\url{https://www.ukbiobank.ac.uk}
\url{https://www.biorxiv.org/content/10.1101/807792v3}
 @examples
 \dontrun{
data("BST1")
finemap_DT <- BST1
# Only including a small subset of the full
# LD matrix for storage purposes.
lead_snp <- subset(finemap_DT, leadSNP)$SNP
snp_list <-  finemap_DT[which(finemap_DT$SNP==lead_snp)-100:which(finemap_DT$SNP==lead_snp)+100,]$SNP
BST1_LD_matrix <- readRDS("../Fine_Mapping/Data/GWAS/Nalls23andMe_2019/BST1/plink/UKB_LD.RDS")
BST1_LD_matrix <- BST1_LD_matrix[snp_list, snp_list]
usethis::use_data(BST1_LD_matrix, overwrite=T)
 }
}
\usage{
BST1_LD_matrix
}
\description{
Precomputed LD within the \emph{BST1} locus
 (defined in \code{\link[=BST1]{BST1}}.
LD derived white British subpopulation in the UK Biobank.
Only includes a subset of all the SNPs for storage purposes
(including the lead GWAS/QTL SNP).
}
\details{
Data originally comes from \href{https://www.ukbiobank.ac.uk}{UK Biobank}.
LD was pre-computed and stored by the Alkes Price lab
(see \href{https://www.biorxiv.org/content/10.1101/807792v3}{here}).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{epigenetics_summary}
\alias{epigenetics_summary}
\title{Summarise \code{HaploR} annotations}
\source{
\href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
}
\usage{
epigenetics_summary(
  merged_results,
  tissue_list = c("BRN", "BLD"),
  epigenetic_variables = c("Promoter_histone_marks", "Enhancer_histone_marks")
)
}
\description{
Summarise \code{HaploR} annotations
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GoShifter.R
\name{GOSHIFTER.create_LD}
\alias{GOSHIFTER.create_LD}
\title{Prepare LD input for \emph{GoShifter}}
\usage{
GOSHIFTER.create_LD(
  locus_dir,
  finemap_dat = NULL,
  LD_path = NULL,
  conda_env = "goshifter",
  nThread = 4,
  verbose = T
)
}
\description{
Prepare LD input for \emph{GoShifter}
}
\examples{
data("BST1"); data("locus_dir")
LD_folder <- GOSHIFTER.create_LD(locus_dir=locus_dir, finemap_dat=subset(BST1, P<5e-8))
}
\seealso{
Other GOSHIFTER: 
\code{\link{GOSHIFTER.create_snpmap}()}
}
\concept{GOSHIFTER}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UKBiobank_LD.R
\name{LD.UKBiobank}
\alias{LD.UKBiobank}
\title{Download LD matrices from UK Biobank}
\usage{
LD.UKBiobank(
  subset_DT = NULL,
  locus_dir,
  sumstats_path = NULL,
  chrom = NULL,
  min_pos = NULL,
  force_new_LD = F,
  chimera = F,
  server = T,
  download_full_ld = F,
  download_method = "direct",
  fillNA = 0,
  nThread = 4,
  return_matrix = F,
  conda_env = "echoR",
  remove_tmps = T,
  verbose = T
)
}
\description{
Download pre-computed LD matrices from UK Biobank in 3Mb windows,
then subset to the region that overlaps with \code{subset_DT}.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GoShifter.R
\name{GOSHIFTER}
\alias{GOSHIFTER}
\title{Run \code{GoShifter} enrichment pipeline}
\source{
\href{https://github.com/immunogenomics/goshifter}{GitHub}
\href{https://pubmed.ncbi.nlm.nih.gov/26140449/}{Publication}
}
\usage{
GOSHIFTER(
  locus_dir,
  snp_df,
  SNP.Group = "",
  goshifter_path = NULL,
  permutations = 1000,
  ROADMAP_search = "",
  ROADMAP_type = NA,
  chromatin_states = c("TssA"),
  R2_filter = 0.8,
  overlap_threshold = 1,
  force_new_goshifter = F,
  remove_tmps = T,
  verbose = T,
  save_results = T
)
}
\description{
\code{GoShifter}: Locus-specific enrichment tool.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_fullSS.R
\name{example_fullSS}
\alias{example_fullSS}
\title{Store \code{\link{fullSS_dat}}}
\usage{
example_fullSS(
  fullSS_path = file.path(tempdir(), "Nalls23andMe_2019.fullSS_subset.tsv"),
  munged = TRUE,
  nThread = 1
)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}
}
\description{
Store \code{\link{fullSS_dat}}
}
\examples{
\dontrun{
data("fullSS_dat")
fullSS_path <- example_fullSS()
}
}
\seealso{
Other general: 
\code{\link{createDT_html}()},
\code{\link{createDT}()},
\code{\link{dt.replace}()},
\code{\link{fillNA_CS_PP}()},
\code{\link{get_os}()},
\code{\link{get_sample_size}()},
\code{\link{tryFunc}()}
}
\concept{general}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.functional_enrichment}
\alias{POLYFUN.functional_enrichment}
\title{Run functional enrichment tests}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.functional_enrichment(
  finemap_dat,
  PP_thresh = 0.95,
  save_plot = "./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/annot_enrichment.png"
)
}
\description{
Run functional enrichment tests
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.plot_enrichment}
\alias{IMPACT.plot_enrichment}
\title{Plot IMPACT score enrichment}
\usage{
IMPACT.plot_enrichment(ENRICH)
}
\description{
Plot IMPACT score enrichment
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NOTT_2019.interactome}
\alias{NOTT_2019.interactome}
\title{Brain cell type-specific enhancers, promoters, and interactomes}
\format{
An object of class \code{list} of length 12.
}
\source{
\url{https://science.sciencemag.org/content/366/6469/1134}
}
\usage{
NOTT_2019.interactome
}
\description{
Originally from \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}.
Specifically: \emph{aay0793-Nott-Table-S5.xlsx}.
}
\examples{
\dontrun{
file <- "~/Desktop/Fine_Mapping/echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S5.xlsx"
sheets <- readxl::excel_sheets(file)
enh_prom_sheets <- grep("enhancers|promoters",sheets,value = T)
other_sheets <- grep("enhancers|promoters",sheets,value = T, invert = T)
NOTT_2019.interactome <- lapply(other_sheets, function(s){readxl::read_excel(file, sheet=s, skip=2)})
NOTT_2019.interactome <- append(NOTT_2019.interactome, lapply(enh_prom_sheets, function(s){readxl::read_excel(file, sheet=s, skip=2, col_names = c("chr","start","end"))}) )
names(NOTT_2019.interactome) <- c(other_sheets, enh_prom_sheets)
usethis::use_data(NOTT_2019.interactome, overwrite = T)
}
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finemapping_portal_API.R
\name{GITHUB.portal_query}
\alias{GITHUB.portal_query}
\title{Search and download fine-mapping files}
\usage{
GITHUB.portal_query(
  dataset_types = NULL,
  datasets = NULL,
  phenotypes = NULL,
  loci = NULL,
  LD_panels = c("UKB", "1KGphase1", "1KGphase3"),
  file_types = c("multi_finemap", "LD", "plot"),
  results_dir = tempdir(),
  overwrite = F,
  nThread = parallel::detectCores() - 2,
  verbose = T
)
}
\description{
Search the \href{https://github.com/RajLabMSSM/Fine_Mapping_Shiny}{echolocatoR Fine-mapping Portal}
for fine-mapping results, LD, and locus plots.
}
\examples{
local_finemap <- GITHUB.portal_query(dataset_types="GWAS",
                                     phenotypes = c("schizophrenia","parkinson"),
                                     file_types = "multi_finemap",
                                     loci = c("BST1","CHRNB1","LRRK2",1:3),
                                     LD_panels=c("UKB","1KGphase3"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.ldsc_annot_enrichment}
\alias{POLYFUN.ldsc_annot_enrichment}
\title{Run and plot heritability enrichment tests}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.ldsc_annot_enrichment(
 
    .results = "Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output/PD_GWAS_LDSC/PD_GWAS_LDSC.results",
  show_plot = T,
  save_plot = F,
  title = "LDSC Heritability Enrichment",
  subtitle = "PD GWAS"
)
}
\description{
Run and plot heritability enrichment tests
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CORCES_2020.scATACseq_peaks}
\alias{CORCES_2020.scATACseq_peaks}
\title{scATACseq peaks from human brain tissue}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 359022 rows and 10 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.scATACseq_peaks
}
\description{
Each row represents an individual peak identified in the single-cell ATAC-seq data.
}
\details{
Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
Specifically: \emph{STable5_Features_scATAC-seq_Peaks_all}
}
\examples{
\dontrun{
dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable5_Features_scATAC-seq_Peaks_all.xlsx", skip = 18)
CORCES_2020.scATACseq_peaks <- data.table::data.table(dat)
usethis::use_data(CORCES_2020.scATACseq_peaks, overwrite = T)
}
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}}
}
\concept{CORCES_2020}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Nalls_top_SNPs}
\alias{Nalls_top_SNPs}
\title{TopSS example file}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 107 rows and 31 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
Nalls_top_SNPs
}
\description{
Summary stats of the top SNP(s) per locus.
Used to query locus subsets.for fine-mapping.
}
\details{
Data from \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}, Table S2.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{SNPs_by_mutation_type}
\alias{SNPs_by_mutation_type}
\title{Return only the missense SNPs}
\source{
\href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
}
\usage{
SNPs_by_mutation_type(merged_results, mutation_type = "missense_variant")
}
\description{
Return only the missense SNPs
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.get_promoter_interactome_data}
\alias{NOTT_2019.get_promoter_interactome_data}
\title{Get cell type-specific promoter/emhancer/interactome data}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.get_promoter_interactome_data(finemap_dat)
}
\description{
Brain cell-specific epigenomic data from Nott et al. (2019).
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{query_fullSS}
\alias{query_fullSS}
\title{Munge the full sumary stats file.}
\usage{
query_fullSS(fullSS_path, subset_path)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}
}
\description{
Read in a standardize the entire full summary stats file at once.
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_handler}()}
}
\concept{query functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.survey_annotation}
\alias{PAINTOR.survey_annotation}
\title{Report the number of annotations > 0}
\usage{
PAINTOR.survey_annotation(PT_results_path, locus_name = "Locus1")
}
\description{
@keywords internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.locus.R
\name{PLOT.get_window_limits}
\alias{PLOT.get_window_limits}
\title{Get window size limits for plot}
\usage{
PLOT.get_window_limits(
  finemap_dat,
  index_as_center = T,
  plot.zoom = NULL,
  genomic_units = "Mb",
  verbose = T
)
}
\description{
Get window size limits for plot
}
\examples{
data("BST1");
xlims <- PLOT.get_window_limits(finemap_dat=BST1, plot.zoom=50000)
xlims <- PLOT.get_window_limits(finemap_dat=BST1, plot.zoom="all")
xlims <- PLOT.get_window_limits(finemap_dat=BST1, plot.zoom="5x")
}
\seealso{
Other plot: 
\code{\link{PLOT.add_multitrack_lines}()},
\code{\link{PLOT.dot_summary}()},
\code{\link{PLOT.get_window_suffix}()}
}
\concept{plot}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.finemapper}
\alias{POLYFUN.finemapper}
\title{Run \emph{PolyFun+SUSIE} fine-mapping pipeline}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.finemapper(
  polyfun = NULL,
  finemap_dat = NULL,
  npz_gz_LD = NULL,
  locus = NULL,
  sample_size = NULL,
  locus_dir,
  n_causal = 5,
  method = "susie",
  h2_path = NULL,
  conda_env = "echoR"
)
}
\description{
Run \emph{PolyFun+SUSIE} fine-mapping pipeline
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.install_reticulate}
\alias{CONDA.install_reticulate}
\title{Install reticulate}
\usage{
CONDA.install_reticulate(
  dependencies = c("devtools", "reticulate", "lattice", "jsonlite", "Matrix",
    "rappdirs", "Rcpp")
)
}
\description{
\emph{reticulate} often doesn't install very well via CRAN.
This function helps do it correctly.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.enrichment_bootstrap}
\alias{XGR.enrichment_bootstrap}
\title{XGR enrichment (bootstrapped)}
\usage{
XGR.enrichment_bootstrap(
  gr,
  merged_dat,
  snp_groups = c("Random", "GWAS lead", "UCS (-PolyFun)", "UCS",
    "Consensus (-PolyFun)", "Consensus"),
  background_filter = NULL,
  grouping_vars = c("Study", "Assay", "Cell_type"),
  iterations = 1000,
  fg_sample_size = 20,
  bg_sample_size = NULL,
  bootstrap = T,
  save_path = F,
  nThread = 4,
  verbose = T
)
}
\description{
XGR enrichment (bootstrapped)
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"

merged_dat <- merge_finemapping_results(dataset = dirname(root), LD_reference = "UKB", minimum_support = 0)
merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)
no_no_loci<- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1", "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3");
merged_dat <- subset(merged_dat, !Locus \%in\% no_no_loci)

gr.merged <- merge_celltype_specific_epigenomics()
grouping_vars <- c("Study","Cell_type","Assay")

enrich_res <- XGR.enrichment_bootstrap(gr=gr.merged, merged_dat=merged_dat, grouping_vars=grouping_vars,  bootstrap=F, snp_groups=c("Random","GWAS lead","UCS (-PolyFun)","UCS","Consensus (-PolyFun)","Consensus"), nThread=12, save_path=file.path(root,"XGR/celltypespecific_epigenomics.SNP_groups.csv.gz"))

enrich_boot <- XGR.enrichment_bootstrap(gr=gr.merged, merged_dat=merged_dat, grouping_vars=grouping_vars,  bootstrap=T, fg_sample_size=NULL, bg_sample_size=100, iterations=100)
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.activate_env}
\alias{CONDA.activate_env}
\title{Activate conda env}
\usage{
CONDA.activate_env(conda_env = "echoR", verbose = T)
}
\description{
Activate conda env
}
\examples{
CONDA.activate_env(conda_env="echoR")
}
\seealso{
Other conda: 
\code{\link{CONDA.create_echoR_env}()},
\code{\link{CONDA.env_from_yaml}()},
\code{\link{CONDA.find_env_Rlib}()},
\code{\link{CONDA.find_python_path}()},
\code{\link{CONDA.install}()}
}
\concept{conda}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MungeSumstats.R
\name{MUNGESUMSTATS.to_echolocatoR}
\alias{MUNGESUMSTATS.to_echolocatoR}
\title{Convert from \pkg{MungeSumstats} to \pkg{echolocatoR} format}
\usage{
MUNGESUMSTATS.to_echolocatoR(dat)
}
\description{
Convert from \pkg{MungeSumstats} to \pkg{echolocatoR} format
}
\examples{
\dontrun{
eduAttainOkbayPth <- system.file("extdata","eduAttainOkbay.txt", package="MungeSumstats")
reformatted <- MungeSumstats::format_sumstats(path=eduAttainOkbayPth, ref_genome="GRCh37", return_data = TRUE)
dat_echoR <- MUNGESUMSTATS.to_echolocatoR(dat=reformatted)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_os.R
\name{get_os}
\alias{get_os}
\title{Identify current operating system (OS).}
\usage{
get_os()
}
\description{
Identify current operating system (OS).
}
\examples{
get_os()
}
\seealso{
Other general: 
\code{\link{createDT_html}()},
\code{\link{createDT}()},
\code{\link{dt.replace}()},
\code{\link{example_fullSS}()},
\code{\link{fillNA_CS_PP}()},
\code{\link{get_sample_size}()},
\code{\link{tryFunc}()}
}
\concept{general}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.get_top_annotations}
\alias{IMPACT.get_top_annotations}
\title{Get the top annotation(s)}
\usage{
IMPACT.get_top_annotations(
  ANNOT_MELT,
  snp.filter = "!is.na(IMPACT_score)",
  top_annotations = 1,
  force_one_annot_per_locus = F
)
}
\description{
Get the annotation(s) with the top mean \emph{IMPACT} for a given set of SNPs.
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{IMPACT_heatmap}()},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAIN.R
\name{finemap_pipeline}
\alias{finemap_pipeline}
\title{Run \pkg{echolocatoR} pipeline on a single locus}
\usage{
finemap_pipeline(
  locus,
  fullSS_path,
  fullSS_genome_build = NULL,
  LD_genome_build = "hg19",
  results_dir,
  dataset_name = "dataset_name",
  dataset_type = "GWAS",
  top_SNPs = "auto",
  force_new_subset = F,
  force_new_LD = F,
  force_new_finemap = T,
  finemap_methods = c("ABF", "FINEMAP", "SUSIE", "POLYFUN_SUSIE"),
  finemap_args = NULL,
  bp_distance = 5e+05,
  n_causal = 5,
  chrom_col = "CHR",
  chrom_type = NULL,
  position_col = "POS",
  snp_col = "SNP",
  pval_col = "P",
  effect_col = "Effect",
  stderr_col = "StdErr",
  tstat_col = "t-stat",
  locus_col = "Locus",
  freq_col = "Freq",
  MAF_col = "MAF",
  A1_col = "A1",
  A2_col = "A2",
  gene_col = "Gene",
  N_cases_col = "N_cases",
  N_controls_col = "N_controls",
  N_cases = NULL,
  N_controls = NULL,
  proportion_cases = "calculate",
  sample_size = NULL,
  LD_reference = "1KGphase1",
  superpopulation = "EUR",
  remote_LD = T,
  download_method = "axel",
  min_POS = NA,
  max_POS = NA,
  min_MAF = NA,
  trim_gene_limits = F,
  max_snps = NULL,
  file_sep = "\\t",
  min_r2 = 0,
  LD_block = F,
  LD_block_size = 0.7,
  vcf_folder = NULL,
  query_by = "tabix",
  remove_variants = F,
  remove_correlates = F,
  probe_path = "./Data/eQTL/gene.ILMN.map",
  conditioned_snps,
  plot_LD = F,
  remove_tmps = T,
  plot.types = c("simple"),
  PAINTOR_QTL_datasets = NULL,
  server = F,
  PP_threshold = 0.95,
  consensus_threshold = 2,
  case_control = T,
  QTL_prefixes = NULL,
  fillNA = 0,
  plot.zoom = "1x",
  plot.Nott_epigenome = F,
  plot.Nott_show_placseq = F,
  plot.Nott_binwidth = 200,
  plot.Nott_bigwig_dir = NULL,
  plot.XGR_libnames = NULL,
  plot.Roadmap = F,
  plot.Roadmap_query = NULL,
  conda_env = "echoR",
  nThread = 4,
  verbose = T
)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{results_dir}{Where to store all results.
\strong{IMPORTANT!:} It is usually best to provide the absolute path rather than the relative path.
This is especially important for \emph{FINEMAP}.}

\item{dataset_name}{The name you want to assign to the dataset being fine-mapped,
This will be used to name the subdirectory where your results will be stored
(e.g. \emph{Data/GWAS/<dataset_name>}).
Don't use special characters (e.g.".", "/").}

\item{dataset_type}{The kind dataset you're fine-mapping (e.g. GWAS, eQTL, tQTL).
This will also be used when creating the subdirectory where your results will be stored
(e.g. \emph{Data/<dataset_type>/Kunkle_2019}).}

\item{top_SNPs}{A data.frame with the genomic coordinates of the lead SNP for each locus.
The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
Only one SNP per \strong{Locus} should be included.
At minimum, \code{top_SNPs} should include the following columns:
\describe{
\item{\emph{Locus}}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
\item{\emph{CHR}}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
\item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
}}

\item{force_new_subset}{By default, if a subset of the full summary stats file for a given locus is already present,
then \pkg{echolocatoR} will just use the pre-existing file.
Set \code{force_new_subset=T} to override this and extract a new subset.
Subsets are saved in the following path structure:
\emph{Data/<dataset_type>/<dataset_name>/<locus>/Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}}

\item{force_new_LD}{By default, if an LD matrix file for a given locus is already present,
then \pkg{echolocatoR} will just use the preexisting file.
Set \code{force_new_LD=T} to override this and extract a new subset.}

\item{force_new_finemap}{By default, if an fine-mapping results file for a given locus is already present,
then \pkg{echolocatoR} will just use the preexisting file.
Set \code{force_new_finemap=T} to override this and re-run fine-mapping.}

\item{finemap_methods}{Which fine-mapping methods you want to use.}

\item{bp_distance}{The width of the window size you want each locus to be.
For example, if \code{bp_distance=500000} then the locus will span 500kb from the lead SNP in either direction,
resulting in a locus that is ~1Mb long (depending on the dataset).}

\item{n_causal}{The maximum number of potential causal SNPs per locus.
This parameter is used somewhat differntly by different fine-mapping tools.
See tool-specific functions for details.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{snp_col}{Name of the SNP RSID column in the full summary stats file.
(\emph{default: ="SNP"})}

\item{pval_col}{Name of the p-value column in the full summary stats file.
Raw p-values are preferred, but if not available corrected p-values (e.g. FDR) can be used instead.
(\emph{default: ="P"})}

\item{effect_col}{Name of the effect size column in the full summary stats file.
Effect size is preferred, but if not available other metrics like Beta for Odds Ratio can be used instead.
(\emph{default: ="Effect"})}

\item{stderr_col}{Name of the standard error  column in the full summary stats file.
You can also set \code{stderr_col="calculate"} to infer standard error using: \code{effect / tstat}.
(\emph{default: ="StdErr"})}

\item{tstat_col}{Name of the t-statistic column in the full summary stats file.
This column is not necessary unless \code{stderr_col="calculate"} or the standard error column is missing.
(\emph{default: ="t-stat"})}

\item{locus_col}{Name of the locus column in the full summary stats file.
(\emph{default: ="Locus"})}

\item{freq_col}{Name of the allele frequency column in the full summary stats file.
Effect allele frequency is preferred, but the non-effect allele can be provided instead (though this may be less accurate).
This column is not necessary unless \code{MAF_col="calculate"} or the MAF column is missing.
(\emph{default: ="Freq"})}

\item{MAF_col}{Name of the minor allele frequency column in the full summary stats file.
Can be inferred from \strong{freq_col} if missing from the dataset.
(\emph{default: ="MAF"})}

\item{A1_col}{Name of the effect/risk allele column in the full summary stats.
 \strong{\emph{IMPORTANT}}: Make sure this actually the case for your full summary stats file.
Unfortunately, different studies report different kinds of allele information in a non-standardized way.
Meaning that A1/A2 can refer to any number of things:
 \describe{
 \item{effect/other alleles}{in the case of diseases}
 \item{ref/alt alleles}{where ref is the reference genome being used}
 \item{major/minor alleles}{This dichotomy holds true for bi-allelic SNPs but not necessary multi-allelic SNPs}
 }
 This makes comparing summary stats across GWAS/QTL datasets very confusing for several reasons:
 \describe{
 \item{Multi-allelic SNPs}{SNPs can have more than just 2 possible alleles (multi-allelic SNPs). Even if you compare the same SNP between two studies, you may accidentally be comparing totally different alleles.}
 \item{Valence}{The valence (+/-) of per-SNP GWAS effect sizes/beta can be relative to different allele types between studies.
 For example, let's say in one GWAS study your effect size for SNP A is 1.5 relative to the major allele in one study,
  and the minor allele happens to be the one found in the reference genome.
  You then try to compare that effect size to that of the same SNP in another GWAS.
  But, the valence of the effect sizes in the 2nd GWAS study are all relative to the reference genome (instead of the minor allele),
  giving the same SNP a value of -1.2. If you took the effect sizes at face value you'd say the signals are in opposite directions.
  But once you take into account how the valences were determined in each study you realize that they're actually both positive relative to the major allele.}
 }
This process of reversing per-SNP valences based on aligning the alleles is known as allele flipping.
This is important when comparing individual SNPs, but can also have an impact on colocalization results.}

\item{gene_col}{For QTL studies, the name of the [e]gene column in the full summary stats file (\emph{default: "gene"}).
This column will be used for filtering summary stats if supplying a named list of gene:Locus pairs to \code{loci}.}

\item{N_cases_col}{Name of the column in the full summary stats that has the number of case subjects in the study.
This can either be per SNP sample sizes, or one number repeated across all rows.
Proxy cases (e.g. relatives of people with the disease being investigated) should be included in this estimate if any were used in the study.
This column is not necesssary if \code{N_cases} parameter is provided.
(\emph{default: ="N_cases"})}

\item{N_controls_col}{Name of the column in the full summary stats that has the number of control subjects in the study.
 This can either be per SNP sample sizes, or one number repeated across all rows.
 This column is not necesssary if \code{N_controls} parameter is provided.
(\emph{default: ="N_controls"})}

\item{N_cases}{The number of case subjects in the study.
Instead of providing a redundant \strong{N_cases_col} column, you can simply enter one value here.}

\item{N_controls}{The number of control subjects in the study.
Instead of providing a redundant \strong{N_controls_col} column, you can simply enter one value here.}

\item{proportion_cases}{The proportion of total subjects in the study that were cases.
if \code{proportion_cases="calculate"} then this is inferred:  \code{N_controls / N_controls}.}

\item{sample_size}{The overall sample size of the study.
If none is given, and \strong{N_cases} and \strong{N_controls} columns are present,
then sample_size is inferred to be:  \code{max(N_cases) + max(N_controls)}.}

\item{LD_reference}{Which linkage disequilibrium reference panel do you want to use.
Options include:
\describe{
\item{"UKB"}{A pre-caclulated LD reference matrix from a subset of caucasian British individuals from the UK Biobank. See \href{https://www.biorxiv.org/content/10.1101/807792v2}{Wiessbrod et al. (2019)} for more details.}
\item{"1KGphase1"}{Download a subset of the 1000 Genomes Project Phase 1 vcf and calculate LD on the fly with plink.}
\item{"1KGphase3"}{Download a subset of the 1000 Genomes Project Phase 3 vcf and calculate LD on the fly with plink.}
\item{"<path>/*.vcf" or "<path>/*.vcf.gz"}{Alternatively, users can provide their own custom panel by supplying a list of \emph{.vcf} file path (one per locus) which \pkg{echolocatoR} will use to compute LD (using \emph{plink}).}
}}

\item{superpopulation}{Subset your LD reference panel by superopulation.
Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
\href{https://www.internationalgenome.org/faq/which-populations-are-part-your-study/}{1KGphase1 options} include:
\describe{
\item{"AFR"}{African [descent]}
\item{"AMR"}{Ad-mixed American}
\item{"EAS"}{East Asian}
\item{"EUR"}{European}
\item{"SAS"}{South Asian}
}}

\item{remote_LD}{When acquiring LD matrices,
the default is to delete the full vcf or npz files after \pkg{echolocatoR} has extracted the necssary subset.
However, if you wish to keep these full files (which can be quite large) set \code{remote_LD=T}.}

\item{min_POS}{Manually set the minimum genomic position for your locus subset.
\code{min_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{min_POS=top_SNPs$min_POS}).}

\item{max_POS}{Manually set the maximum genomic position for your locus subset.
\code{max_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{max_POS=top_SNPs$max_POS}).}

\item{min_MAF}{Remove any SNPs with \strong{MAF} < \code{min_MAF}.}

\item{trim_gene_limits}{If a valid gene symbol is provided to \code{trim_gene_limits},
the gene's canonical coordinates are pulled from \code{biomaRt}.
This includes introns, exons, and proximal regulatory regions (e.g. promoters).
Any SNPs that fall outside these coordinates are remove from downstream fine-mapping.
Set \code{trim_gene_limits=F} to not limit by gene coordinates (\emph{default}).}

\item{max_snps}{The maximum number of SNPs to include in the locus.
If the current window size yields > \code{max_snps},
 then the outer edges of the of the locus are trimmed until the number of SNPs ≤ \code{max_snps}.}

\item{file_sep}{The separator in the full summary stats file.
This parameter is only necessary if \code{query_by!="tabix"}.}

\item{min_r2}{Remove any SNPs are below the LD r2 threshold with the lead SNP within their respective locus.}

\item{LD_block}{Calculate LD blocks with \emph{plink} and only include the block to which the lead SNP belongs.}

\item{LD_block_size}{Adjust the granularity of block sizes when \code{LD_block=T}.}

\item{query_by}{Choose which method you want to use to extract locus subsets from the full summary stats file.
Methods include:
\describe{
\item{"tabix"}{Convert the full summary stats file in an indexed tabix file. Makes querying lightning fast after the initial conversion is done. (\emph{default})}
\item{"coordinates"}{Extract locus subsets using min/max genomic coordinates with \emph{awk}.}
}}

\item{remove_variants}{A list of variants to remove from the locus subset file.}

\item{remove_correlates}{A named list, where the names are the RSIDs of SNPs
whose LD correlates you wish to remove,
and the value is the absolute r2 threshold you wish to filter at for each RSID respectively
(e.g. \code{ remove_correlates = c("rs76904798"=.2, "rs10000737"=.8)}).
This will also remove the SNPs in \code{remove_correlates} themselves.}

\item{probe_path}{The location of the file containing translations between probe IDs and gene symbols.
Only used for certain eQTL datasets.}

\item{conditioned_snps}{Which SNPs to conditions on when fine-mapping with \emph{COJO}.}

\item{plot_LD}{Whether to plot a subset of the LD matix.}

\item{remove_tmps}{Whether to remove any temporary files (e.g. FINEMAP output files) after the pipeline is done running.}

\item{plot.types}{Which kinds of plots to include.
Options:
\describe{
\item{"simple"}{Just plot the following tracks: GWAS, fine-mapping, gene models}
\item{"fancy"}{Additionally plot XGR annotation tracks (XGR, Roadmap, Nott).}
}}

\item{PAINTOR_QTL_datasets}{A list of QTL datasets to be used when conducting joint functional fine-mapping with \emph{PAINTOR}.}

\item{server}{Whether \pkg{echolocatoR} is being run on a computing cluster/server or on a local machine.}

\item{PP_threshold}{The minimum fine-mapped posterior probability for a SNP to be considered part of a Credible Set.
For example, \code{PP_threshold=.95} means that all Credible Set SNPs will be 95\% Credible Set SNPs.}

\item{consensus_threshold}{The minimum number of fine-mapping tools that include a SNP
in their 95\% Credible Sets to consider that it a "Consensus SNP" (\emph{default=2}).}

\item{plot.zoom}{Zoom into the center of the locus when plotting (without editing the fine-mapping results file).
You can provide either:
\itemize{
\item{The size of your plot window in terms of basepairs (e.g. \code{plot.zoom=50000} for a 50kb window)}.
\item{How much you want to zoom in (e.g. \code{plot.zoom="1x"} for the full locus, \code{plot.zoom="2x"} for 2x zoom into the center of the locus, etc.)}.
}
You can pass a list of window sizes (e.g. \code{c(50000,100000,500000)}) to automatically generate
multiple views of each locus.
This can even be a mix of different style inputs: e.g. \code{c("1x","4.5x",25000)}.}

\item{plot.Nott_binwidth}{When including Nott et al. (2019) epigenomic data in the track plots,
adjust the bin width of the histograms.}

\item{plot.Nott_bigwig_dir}{Instead of pulling Nott et al. (2019) epigenomic data
from the \emph{UCSC Genome Browser}, use a set of local bigwig files.}

\item{plot.Roadmap}{Find and plot annotations from Roadmap.}

\item{plot.Roadmap_query}{Only plot annotations from Roadmap whose metadata contains a string or any items from  a list of strings
(e.g. \code{"brain"} or \code{c("brain","liver","monocytes")}).}

\item{conda_env}{The name of a conda environment to use.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}

\item{loci}{Character list of loci in \strong{Locus} col of \code{top_SNPs}.}

\item{min_Dprime}{Remove any SNPs are below the LD D' threshold with the lead SNP within their respective locus.
This is paramter currently only works when \code{LD_reference!="UKB"}.}
}
\description{
Unlike \code{finemap_loci}, you don't need to provide a \code{top_SNPs} data.frame.
Instead, just manually provide the coordinates of the locus you want to fine-map.
}
\details{
The primary functions of \pkg{echolocatoR} that expedite fine-mapping
 by wrapping many other \pkg{echolocatoR} functions into one.
 Encompasses steps including:
 \describe{
 \item{Subset & standardize}{Extract subsets of the full summmary stats GWAS or QTL file and reformat them to be compatible with \pkg{echolocatoR}'s various functions }
 \item{Calculate linkage disequilibrium}{Download and prepare the necessary LD matrix.}
 \item{Fine-map}{Run various fine-mapping tools and merge the results into a single multi-finemap data.frame.}
 \item{Plot}{Summarise the results in a multi-track plot for each locus.}
 }
}
\section{input file parameters}{

}

\section{input file column names}{

}

\section{overwrite existing files}{

}

\section{fine-mapping parameters}{

}

\seealso{
Other MAIN: 
\code{\link{finemap_loci}()}
}
\concept{MAIN}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.get_SNPgroup_counts}
\alias{SUMMARISE.get_SNPgroup_counts}
\title{Tally locus-specific SNP group sizes}
\usage{
SUMMARISE.get_SNPgroup_counts(merged_dat, grouping_vars = "Locus")
}
\description{
Tally locus-specific SNP group sizes
}
\examples{
data("merged_DT");
snp_groups <- SUMMARISE.get_SNPgroup_counts(merged_dat=merged_DT)
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{results_report}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.filter_sources}
\alias{XGR.filter_sources}
\title{Filter sources}
\usage{
XGR.filter_sources(gr.lib, n_top_sources = 5)
}
\description{
Identify the sources with the most annotations in the locus.
Then only keep these sources.
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.leadSNP_block}
\alias{LD.leadSNP_block}
\title{Identify the LD block in which the lead SNP resides}
\usage{
LD.leadSNP_block(leadSNP, LD_folder, LD_block_size = 0.7)
}
\description{
Identify the LD block in which the lead SNP resides
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.custom_panel}
\alias{LD.custom_panel}
\title{Compute LD from user-supplied vcf file}
\usage{
LD.custom_panel(
  LD_reference,
  fullSS_genome_build = "hg19",
  LD_genome_build = "hg19",
  subset_DT,
  locus_dir,
  force_new_LD = F,
  min_r2 = F,
  remove_correlates = F,
  fillNA = 0,
  LD_block = F,
  LD_block_size = 0.7,
  remove_tmps = T,
  nThread = 4,
  conda_env = "echoR",
  verbose = T
)
}
\description{
Compute LD from user-supplied vcf file
}
\examples{
\dontrun{
if(!"gaston" \%in\% row.names(installed.packages())){install.packages("gaston")}
data("BST1"); data("locus_dir")
LD_reference="~/Desktop/results/Reference/custom_panel_chr4.vcf"
LD_matrix <- LD.custom_panel(LD_reference=LD_reference, subset_DT=BST1, locus_dir=locus_dir)

locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/QTL/Microglia_all_regions/BIN1"
subset_DT <- data.table::fread(file.path(locus_dir,"Multi-finemap/BIN1.Microglia_all_regions.1KGphase3_LD.Multi-finemap.tsv.gz"))
LD_reference = "/sc/hydra/projects/pd-omics/glia_omics/eQTL/post_imputation_filtering/eur/filtered_variants/AllChr.hg38.sort.filt.dbsnp.snpeff.vcf.gz"
LD_matrix <- LD.custom_panel(LD_reference=LD_reference, subset_DT=BST1, locus_dir=locus_dir, LD_genome_build="hg38")
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/message_parallel.r
\name{message_parallel}
\alias{message_parallel}
\title{Message parallel}
\usage{
message_parallel(...)
}
\value{
Null
}
\description{
Send messages to console even from within parallel processes
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CORCES_2020.HiChIP_FitHiChIP_loop_calls}
\alias{CORCES_2020.HiChIP_FitHiChIP_loop_calls}
\title{FitHiChIP loop calls from human brain tissue}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 11542 rows and 11 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.HiChIP_FitHiChIP_loop_calls
}
\description{
FitHiChIP loop calls that overlap SNPs derived from analysis of H3K27ac HiChIP data.
Each row represents an individual peak identified from the feature binarization analysis (see methods).
}
\details{
Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
Specifically: \emph{STable10_Coacessibility_Peak_loop_connection}, \emph{HiChIP FitHiChIP Loop Calls} sheet.

\strong{Column dictionary}
\describe{
\item{hg38_Chromosome_Anchor1}{The hg38 chromosome of the first loop Anchor.}
\item{hg38_Start_Anchor1}{The hg38 start position of the first loop Anchor.}
\item{hg38_Stop_Anchor1}{The hg38 stop position of the first loop Anchor.}
\item{Width_Anchor1}{The width of the first loop Anchor.}
\item{hg38_Chromosome_Anchor2}{The hg38 chromosome of the second loop Anchor.}
\item{hg38_Start_Anchor2}{The hg38 start position of the second loop Anchor.}
\item{hg38_Stop_Anchor2}{The hg38 stop position of the second loop Anchor.}
\item{Width_Anchor2}{The width of the second loop Anchor.}
\item{Score}{The -log10(q-value) of the loop call from FitHiChIP.}
\item{Anchor1_hasSNP}{A boolean variable determining whether the first anchor overlaps a SNP from our AD/PD GWAS analyses.}
\item{Anchor2_hasSNP}{A boolean variable determining whether the second anchor overlaps a SNP from our AD/PD GWAS analyses.}
}
}
\examples{
\dontrun{
dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable10_Coacessibility_Peak_loop_connection.xlsx", skip = 19, sheet=1)
CORCES_2020.HiChIP_FitHiChIP_loop_calls <- data.table::data.table(dat)
usethis::use_data(CORCES_2020.HiChIP_FitHiChIP_loop_calls)
}
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_subset_path}
\alias{get_subset_path}
\title{Construct the path of the locus subset}
\usage{
get_subset_path(
  subset_path = "auto",
  results_dir = "./results",
  dataset_type = "dataset_type",
  dataset_name = "dataset_name",
  locus = NULL,
  suffix = ".tsv.gz"
)
}
\description{
Construct the path of the locus subset
}
\examples{
subset_path <- get_subset_path(results_dir="./Data/GWAS/Nalls23andMe_2019/BST1", locus="BST1")
}
\seealso{
Other directory functions: 
\code{\link{get_locus_dir}()},
\code{\link{get_multifinemap_path}()},
\code{\link{get_study_dir}()},
\code{\link{make_locus_dir}()}
}
\concept{directory functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{query_handler}
\alias{query_handler}
\title{Handles which query method to use}
\usage{
query_handler(
  fullSS_path,
  locus_dir = NULL,
  top_SNPs = NULL,
  subset_path,
  min_POS = NA,
  max_POS = NA,
  bp_distance = 5e+05,
  locus_col = "Gene",
  chrom_col = "CHR",
  chrom_type = NULL,
  position_col = "POS",
  file_sep = "\\t",
  query_by = "coordinates",
  probe_path = "./Data/eQTL/gene.ILMN.map",
  conda_env = "echoR",
  verbose = T
)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{top_SNPs}{A data.frame with the genomic coordinates of the lead SNP for each locus.
The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
Only one SNP per \strong{Locus} should be included.
At minimum, \code{top_SNPs} should include the following columns:
\describe{
\item{\emph{Locus}}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
\item{\emph{CHR}}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
\item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
}}

\item{subset_path}{Path of the resulting locus subset file.}

\item{min_POS}{Manually set the minimum genomic position for your locus subset.
\code{min_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{min_POS=top_SNPs$min_POS}).}

\item{max_POS}{Manually set the maximum genomic position for your locus subset.
\code{max_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{max_POS=top_SNPs$max_POS}).}

\item{bp_distance}{The width of the window size you want each locus to be.
For example, if \code{bp_distance=500000} then the locus will span 500kb from the lead SNP in either direction,
resulting in a locus that is ~1Mb long (depending on the dataset).}

\item{locus_col}{Name of the locus column in the full summary stats file.
(\emph{default: ="Locus"})}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{file_sep}{The separator in the full summary stats file.
This parameter is only necessary if \code{query_by!="tabix"}.}

\item{query_by}{Choose which method you want to use to extract locus subsets from the full summary stats file.
Methods include:
\describe{
\item{"tabix"}{Convert the full summary stats file in an indexed tabix file. Makes querying lightning fast after the initial conversion is done. (\emph{default})}
\item{"coordinates"}{Extract locus subsets using min/max genomic coordinates with \emph{awk}.}
}}

\item{probe_path}{The location of the file containing translations between probe IDs and gene symbols.
Only used for certain eQTL datasets.}

\item{conda_env}{The name of a conda environment to use.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Handles which query method to use
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()}
}
\concept{query functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install_tricky_packages.R
\name{install_tricky_packages}
\alias{install_tricky_packages}
\title{Install tricky packages}
\usage{
install_tricky_packages()
}
\description{
Some packages don't install very well via the DESCRIPTION file
(e.g. wrong versions, wrong sources).
This function ensures they're actually installed properly.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MungeSumstats.R
\name{MUNGESUMSTATS.col_map}
\alias{MUNGESUMSTATS.col_map}
\title{Column name mappings from \pkg{MungeSumstats} to \pkg{echolocatoR}}
\source{
https://github.com/neurogenomics/MungeSumstats
}
\usage{
MUNGESUMSTATS.col_map()
}
\description{
Column name mappings from \pkg{MungeSumstats} to \pkg{echolocatoR}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SURE.R
\name{SURE.plot}
\alias{SURE.plot}
\title{Plot \emph{SuRE} results across SNP groups}
\source{
\href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
\href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
\href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
\href{https://sure.nki.nl}{SNP-SuRE data browse}
}
\usage{
SURE.plot(
  sure.melt,
  snp_groups = c("GWAS lead", "UCS", "Consensus"),
  comparisons_filter = function(x) {     if ("Consensus" \%in\% x)          return(x) },
  title = "SuRE MPRA",
  xlabel = "SNP Group",
  facet_formula = ". ~ Cell_type",
  show_padj = T,
  show_signif = T,
  vjust_signif = 0.5,
  show_plot = T,
  save_path = F,
  height = 5,
  width = 5
)
}
\description{
Plot \emph{SuRE} results across SNP groups
}
\examples{
sure.melt <- data.table::fread("/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/SURE/Nalls23andMe_2019.SURE.snp_groups.mean.csv.gz")
pb <- SURE.plot(sure.melt=sure.melt)
}
\seealso{
Other SURE: 
\code{\link{SURE.melt_snp_groups}()},
\code{\link{SURE.merge}()}
}
\concept{SURE}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{prepare_mat_meta}
\alias{prepare_mat_meta}
\title{Prepare \emph{IMPACT} data for for \pkg{ComplexHeatmap}}
\usage{
prepare_mat_meta(
  TOP_IMPACT,
  TOP_IMPACT_all,
  snp.group = "Consensus",
  value.var = "mean_IMPACT",
  fill_na = 0
)
}
\description{
Prepare \emph{IMPACT} data for for \pkg{ComplexHeatmap}
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{IMPACT_heatmap}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.vcf_to_bed}
\alias{LD.vcf_to_bed}
\title{Convert vcf file to BED file}
\usage{
LD.vcf_to_bed(vcf.gz.subset, locus_dir, plink_prefix = "plink", verbose = T)
}
\arguments{
\item{vcf.gz.subset}{Path to the gzipped locus subset vcf.}

\item{locus_dir}{Locus-specific results directory.}
}
\description{
Uses plink to convert vcf to BED.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.1KG_download_vcf}
\alias{LD.1KG_download_vcf}
\title{Download vcf subset from 1000 Genomes}
\usage{
LD.1KG_download_vcf(
  subset_DT,
  LD_reference = "1KGphase1",
  remote_LD = T,
  vcf_folder = NULL,
  locus_dir,
  locus = NULL,
  whole_vcf = F,
  download_method = "wget",
  force_new_vcf = F,
  query_by_regions = F,
  nThread = 4,
  conda_env = "echoR",
  verbose = T
)
}
\arguments{
\item{LD_reference}{Which linkage disequilibrium reference panel do you want to use.
Options include:
\describe{
\item{"UKB"}{A pre-caclulated LD reference matrix from a subset of caucasian British individuals from the UK Biobank. See \href{https://www.biorxiv.org/content/10.1101/807792v2}{Wiessbrod et al. (2019)} for more details.}
\item{"1KGphase1"}{Download a subset of the 1000 Genomes Project Phase 1 vcf and calculate LD on the fly with plink.}
\item{"1KGphase3"}{Download a subset of the 1000 Genomes Project Phase 3 vcf and calculate LD on the fly with plink.}
\item{"<path>/*.vcf" or "<path>/*.vcf.gz"}{Alternatively, users can provide their own custom panel by supplying a list of \emph{.vcf} file path (one per locus) which \pkg{echolocatoR} will use to compute LD (using \emph{plink}).}
}}

\item{remote_LD}{When acquiring LD matrices,
the default is to delete the full vcf or npz files after \pkg{echolocatoR} has extracted the necssary subset.
However, if you wish to keep these full files (which can be quite large) set \code{remote_LD=T}.}

\item{query_by_regions}{You can make queries with \code{tabix} in two different ways:
\describe{
\item{\code{query_by_regions=F} \emph{(default)}}{Return a vcf with all positions between the min/max in \code{subset_DT} Takes up more storage but is MUCH faster}
\item{\code{query_by_regions=T}}{Return a vcf with only the exact positions present in \code{subset_DT}. Takes up less storage but is MUCH slower}
}}

\item{conda_env}{The name of a conda environment to use.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Download vcf subset from 1000 Genomes
}
\examples{
\dontrun{
data("BST1");
subset_DT <- BST1
vcf_subset.popDat <- LD.1KG_download_vcf(subset_DT=BST1, LD_reference="1KGphase1", locus_dir=file.path("~/Desktop",locus_dir))
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.postprocess_annotations}
\alias{IMPACT.postprocess_annotations}
\title{Prepare \emph{IMPACT} annotations}
\usage{
IMPACT.postprocess_annotations(ANNOT_MELT, order_loci = T, no_no_loci = NULL)
}
\description{
Transform \emph{IMPACT} annotations from wide format (one row per SNP) to
long format (multiple rows per SNP).
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{IMPACT_heatmap}()},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABF.R
\name{ABF}
\alias{ABF}
\title{Fine-map with ABF}
\source{
\itemize{
\item JB Maller et al., Bayesian refinement of association signals for 14 loci in 3 common diseases. Nature Genetics. 44, 1294–1301 (2012).
\item J Wakefield, A bayesian measure of the probability of false discovery in genetic epidemiology studies. American Journal of Human Genetics. 81, 208–227 (2007).
}
}
\usage{
ABF(
  subset_DT,
  PP_threshold = 0.95,
  sample_size = NULL,
  sdY = NULL,
  case_control = T,
  verbose = T
)
}
\description{
Conduct statistical (non-functional) fine-mapping with approximate Bayes factor (ABF).
}
\examples{
# GWAS
data("BST1");
finemap_DT <- ABF(subset_DT=BST1)

# QTL
locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/QTL/Microglia_all_regions/BIN1"
subset_DT <- data.table::fread(file.path(locus_dir,"/Multi-finemap/BIN1.Microglia_all_regions.1KGphase3_LD.Multi-finemap.tsv.gz"))
finemap_DT <- ABF(subset_DT=subset_DT, case_control=F, sample_size=90)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UKBiobank_LD.R
\name{LD.rds_to_npz}
\alias{LD.rds_to_npz}
\title{Convert .RDS file back to .npz format}
\usage{
LD.rds_to_npz(rds_path, conda_env = "echoR", verbose = T)
}
\description{
Convert .RDS file back to .npz format
}
\examples{
\dontrun{
data("BST1")
npz_path <- LD.rds_to_npz(rds_path="/Users/schilder/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/BST1/plink/UKB_LD.RDS")
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{results_report}
\alias{results_report}
\title{Give a quick summary report of the fine-mapping results}
\usage{
results_report(merged_dat)
}
\description{
Give a quick summary report of the fine-mapping results
}
\examples{
data("merged_DT");
results_report(merged_DT)
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{query_by_coordinates}
\alias{query_by_coordinates}
\title{Use \emph{awk} to query locus subsets.}
\usage{
query_by_coordinates(
  top_SNPs,
  locus,
  subset_path,
  fullSS_path,
  file_sep,
  chrom_col,
  position_col,
  min_POS,
  max_POS,
  bp_distance
)
}
\arguments{
\item{top_SNPs}{A data.frame with the genomic coordinates of the lead SNP for each locus.
The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
Only one SNP per \strong{Locus} should be included.
At minimum, \code{top_SNPs} should include the following columns:
\describe{
\item{\emph{Locus}}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
\item{\emph{CHR}}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
\item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
}}

\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{file_sep}{The separator in the full summary stats file.
This parameter is only necessary if \code{query_by!="tabix"}.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{min_POS}{Manually set the minimum genomic position for your locus subset.
\code{min_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{min_POS=top_SNPs$min_POS}).}

\item{max_POS}{Manually set the maximum genomic position for your locus subset.
\code{max_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{max_POS=top_SNPs$max_POS}).}

\item{bp_distance}{The width of the window size you want each locus to be.
For example, if \code{bp_distance=500000} then the locus will span 500kb from the lead SNP in either direction,
resulting in a locus that is ~1Mb long (depending on the dataset).}
}
\description{
Search full summary stats file by genomic coordinates.
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.filter_vcf}
\alias{LD.filter_vcf}
\title{Filter a vcf by min/max coordinates}
\usage{
LD.filter_vcf(vcf_subset, popDat, superpopulation, remove_tmp = T, verbose = T)
}
\arguments{
\item{vcf_subset}{Path to the locus subset vcf.}

\item{popDat}{The metadata file listing the superpopulation
to which each sample belongs.}

\item{superpopulation}{Subset your LD reference panel by superopulation.
Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
\href{https://www.internationalgenome.org/faq/which-populations-are-part-your-study/}{1KGphase1 options} include:
\describe{
\item{"AFR"}{African [descent]}
\item{"AMR"}{Ad-mixed American}
\item{"EAS"}{East Asian}
\item{"EUR"}{European}
\item{"SAS"}{South Asian}
}}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Uses \emph{bcftools} to filter a vcf by min/max genomic coordinates
(in basepairs).
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{echolocatoR-package}
\alias{echolocatoR}
\alias{echolocatoR-package}
\title{echolocatoR: Automated genomic fine-mapping}
\description{
Automated statistical and functional fine-mapping with extensive access to genome-wide datasets.
}
\details{
Fine-mapping methods are a powerful means of identifying causal variants underlying a given phenotype,
but are  underutilized due to the technical challenges of implementation.
\emph{echolocatoR} is an R package that automates end-to-end genomics fine-mapping, annotation,
and plotting in order to identify the most probable causal variants associated with a given phenotype.

It requires minimal input from users (a GWAS or QTL summary statistics file),
and includes a suite of statistical and functional fine-mapping tools.
It also includes extensive access to datasets
(linkage disequilibrium panels, epigenomic and genome-wide annotations, QTL).

The elimination of data gathering and preprocessing steps enables rapid fine-mapping of many loci in any phenotype,
 complete with locus-specific publication-ready figure generation.
 All results are merged into a single per-SNP summary file for additional downstream analysis
  and results sharing. Therefore \emph{echolocatoR} drastically reduces the barriers to identifying
  causal variants by making the entire fine-mapping pipeline rapid, robust and scalable.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/RajLabMSSM/echolocatoR}
  \item Report bugs at \url{https://github.com/RajLabMSSM/echolocatoR/issues}
}

}
\author{
\strong{Maintainer}: Brian Schilder \email{brian_schilder@alumni.brown.edu} (\href{https://orcid.org/0000-0001-5949-2191}{ORCID})

Authors:
\itemize{
  \item Jack Humphrey \email{Jack.Humphrey@mssm.edu} (\href{https://orcid.org/0000-0002-6274-6620}{ORCID})
  \item Towfique Raj \email{towfique.raj@mssm.edu} (\href{https://orcid.org/0000-0002-9355-5704}{ORCID})
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.plot_impact_score}
\alias{IMPACT.plot_impact_score}
\title{Locus plot of IMPACT scores}
\usage{
IMPACT.plot_impact_score(annot_melt, save_path = F, show_plot = T)
}
\description{
Locus plot of IMPACT scores
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.locus.R
\name{PLOT.heights_dict}
\alias{PLOT.heights_dict}
\title{Plot heights dictionary}
\usage{
PLOT.heights_dict(keys = NULL, default_height = 1)
}
\description{
Plot heights dictionary
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COJO.R
\name{COJO.stepwise}
\alias{COJO.stepwise}
\title{Run \emph{GCTA-COJO} conditional stepwise procedure}
\usage{
COJO.stepwise(
  GCTA_path = system.file("tools/gcta_1.92.1beta5_mac/bin", "gcta64", package =
    "echolocatoR"),
  locus_dir,
  min_MAF,
  excluded_path
)
}
\description{
Runs the \emph{GCTA-COJO} conditional stepwise procedure to identify independent signals.
Should only be run on full, genome-wide summary stats at once.
}
\seealso{
Other COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}: 
\code{\link{COJO}()},
\code{\link{get_stepwise_results}()},
\code{\link{process_COJO_results}()}
}
\concept{COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sample_size.R
\name{get_sample_size}
\alias{get_sample_size}
\title{Infer sample size from summary stats}
\usage{
get_sample_size(subset_DT, sample_size = NULL, effective_ss = T, verbose = T)
}
\description{
Infer sample size from summary stats
}
\examples{
data("BST1")
BST1 <- finemap_DT
subset_DT <- get_sample_size(subset_DT = finemap_DT)
}
\seealso{
Other general: 
\code{\link{createDT_html}()},
\code{\link{createDT}()},
\code{\link{dt.replace}()},
\code{\link{example_fullSS}()},
\code{\link{fillNA_CS_PP}()},
\code{\link{get_os}()},
\code{\link{tryFunc}()}
}
\concept{general}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.list_paintor_annotations}
\alias{PAINTOR.list_paintor_annotations}
\title{List available PAINTOR annotations}
\usage{
PAINTOR.list_paintor_annotations(
  annotations_dir = file.path("/sc/orga/projects/pd-omics/brian/PAINTOR_V3.0",
    "Annotation_directory/Functional_Annotations")
)
}
\description{
List available PAINTOR annotations
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COJO.R
\name{COJO}
\alias{COJO}
\title{Run \emph{GCTA-COJO}}
\usage{
COJO(
  subset_DT,
  fullSS_path,
  locus_dir,
  conditioned_snps,
  excluded_snps = "",
  min_MAF = 0,
  GCTA_path = system.file("tools/gcta_1.92.1beta5_mac/bin", "gcta64", package =
    "echolocatoR"),
  bfiles = "plink_tmp/plink",
  stepwise_procedure = T,
  conditional_analysis = T,
  snp_col = "SNP",
  freq_col = "Freq",
  effect_col = "Effect",
  stderr_col = "StdErr",
  pval_col = "P",
  A1_col = "A1",
  A2_col = "A2",
  full_genome = F
)
}
\description{
Main function to run either the conditional stepwise procedure (genome-wide)
or the conditional analysis (locus-specific) from \emph{GCTA-COJO}.
}
\seealso{
Other COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}: 
\code{\link{COJO.stepwise}()},
\code{\link{get_stepwise_results}()},
\code{\link{process_COJO_results}()}
}
\concept{COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.download_annotations}
\alias{PAINTOR.download_annotations}
\title{Download annotations for PAINTOR}
\usage{
PAINTOR.download_annotations(
  PT_results_path,
  locus_name,
  locus_DT,
  XGR_dataset = NA,
  ROADMAP_search = NA,
  chromatin_state = "TssA",
  no_annotations = F
)
}
\description{
Download annotations for PAINTOR
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.get_interactome}
\alias{NOTT_2019.get_interactome}
\title{Import cell type-specific interactomes}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.get_interactome(
  annot_sub,
  top.consensus.pos,
  marker_key,
  verbose = T
)
}
\description{
Brain cell-specific epigenomic data from Nott et al. (2019).
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extra_installs.R
\name{extra_installs}
\alias{extra_installs}
\title{Install tricky R packages}
\usage{
extra_installs(
  cran_packages = F,
  xgr = T,
  bioc_packages = T,
  github_packages = T
)
}
\description{
Some of \pkg{echolocatoR}'s optional R package dependencies are especially tricky to install.
Rather than including them in the DESCRIPTION file,
which would require them to install \pkg{echolocatoR} at all,
this functions installs them afterwards.
Only packages not already installed will be installed.
}
\details{
Some of the main packages installed via this function include:
\describe{
\item{gaston}{CRAN}
\item{plotly}{CRAN}

\item{XGR}{xgr}
\item{foreign}{xgr}
\item{refGenome}{xgr}

\item{Rgraphviz}{Bioconductor}
\item{biomaRt}{Bioconductor}

\item{knitrBootstrap}{GitHub}
\item{susieR}{GitHub}
}
}
\examples{
library(echolocatoR)
extra_installs()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpliceAI.R
\name{SPLICEAI.plot}
\alias{SPLICEAI.plot}
\title{Plot \emph{SpliceAI} predictions}
\source{
\href{https://github.com/Illumina/SpliceAI}{GitHub}
\href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
}
\usage{
SPLICEAI.plot(dat_merged)
}
\description{
Plot \emph{SpliceAI} predictions
}
\seealso{
Other SpliceAI: 
\code{\link{SPLICEAI.run}()},
\code{\link{SPLICEAI.snp_probs}()},
\code{\link{SPLICEAI.subset_precomputed_tsv_iterate}()},
\code{\link{SPLICEAI.subset_precomputed_tsv}()},
\code{\link{SPLICEAI.subset_precomputed_vcf}()}
}
\concept{SpliceAI}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.merge_and_process}
\alias{XGR.merge_and_process}
\title{Standardize XGR annotations}
\usage{
XGR.merge_and_process(grl.xgr, lib, n_top_sources = 10)
}
\arguments{
\item{grl.xgr}{\code{\link[GenomicRanges]{GenomicRangesList}} of XGR queries.}
}
\description{
Parses the metadata and adds it as columns,
and then merges the results into a single
\code{\link[GenomicRanges]{GenomicRangesList}}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tryFunc}
\alias{tryFunc}
\title{tryCatch extension}
\usage{
tryFunc(input, func, ...)
}
\arguments{
\item{input}{Function input.}

\item{func}{Function.}
}
\description{
Extension of tryCatch function.
}
\seealso{
Other general: 
\code{\link{createDT_html}()},
\code{\link{createDT}()},
\code{\link{dt.replace}()},
\code{\link{example_fullSS}()},
\code{\link{fillNA_CS_PP}()},
\code{\link{get_os}()},
\code{\link{get_sample_size}()}
}
\concept{general}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{VALIDATION.permute}
\alias{VALIDATION.permute}
\title{Validation permutation tests}
\source{
\href{https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/}{Wilcox test tutorial}
}
\usage{
VALIDATION.permute(
  metric_df,
  metric,
  locus_means = F,
  snp_groups = c("Random", "GWAS lead", "UCS (-PolyFun)", "UCS",
    "Consensus (-PolyFun)", "Consensus")
)
}
\description{
Validation permutation tests
}
\examples{
\dontrun{
save_path <- root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"

#### h2 ####
## mean
path <- file.path(root,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.csv.gz")
metric <- "h2.enrichment"

#### IMPACT ####
## mean
path <- file.path(root,"IMPACT/TOP_IMPACT_all.csv.gz");
metric <- "mean_IMPACT"
## raw
path <- file.path(root,"IMPACT/IMPACT_overlap.csv.gz")
metric <- "IMPACT_score"


## Import and Process ##
metric_df <- data.table::fread(path, nThread=8)
if(metric=="IMPACT_score") metric_df <- subset(metric_df, select=c(SNP,leadSNP,ABF.Credible_Set,ABF.PP,SUSIE.Credible_Set,SUSIE.PP,POLYFUN_SUSIE.Credible_Set,POLYFUN_SUSIE.PP,FINEMAP.Credible_Set,FINEMAP.PP,Consensus_SNP,Support,Locus,IMPACT_score))
if(metric=="mean_IMPACT") metric_df <- find_consensus_SNPs_no_PolyFun(metric_df)

#### run bootstrap ####
permute_res <- VALIDATION.permute(metric_df=metric_df, metric=metric )
}
}
\seealso{
Other VALIDATION: 
\code{\link{VALIDATION.bootstrap_multimetric}()},
\code{\link{VALIDATION.bootstrap_plot}()},
\code{\link{VALIDATION.bootstrap}()},
\code{\link{VALIDATION.compare_bootstrap_distributions}()},
\code{\link{VALIDATION.super_plot}()}
}
\concept{VALIDATION}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FINEMAP.R
\name{FINEMAP.find_executable}
\alias{FINEMAP.find_executable}
\title{Retrieve location of \code{FINEMAP} executable}
\source{
\url{http://www.christianbenner.com}
}
\usage{
FINEMAP.find_executable(
  FINEMAP_path = NULL,
  OS = NULL,
  version = "1.4",
  verbose = T
)
}
\description{
Retrieve location of \code{FINEMAP} executable
}
\examples{
FINEMAP_path <- FINEMAP.find_executable()
}
\seealso{
Other FINEMAP: 
\code{\link{FINEMAP.construct_data}()},
\code{\link{FINEMAP.construct_master}()},
\code{\link{FINEMAP.process_results}()},
\code{\link{FINEMAP}()}
}
\concept{FINEMAP}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.iterate_enrichment}
\alias{XGR.iterate_enrichment}
\title{Conduct enrichment tests for each annotation}
\usage{
XGR.iterate_enrichment(
  subset_DT,
  foreground_filter = "Consensus_SNP",
  background_filter = "leadSNP",
  lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes",
    "ENCODE_DNaseI_ClusteredV3_CellTypes", "Broad_Histone", "FANTOM5_Enhancer",
    "Segment_Combined_Gm12878", "TFBS_Conserved", "ReMap_PublicAndEncode_TFBS",
    "Blueprint_VenousBlood_Histone", "Blueprint_DNaseI", "FANTOM5_CAT_Cell",
    "FANTOM5_CAT_MESH", "GWAScatalog_alltraits"),
  save_path = F,
  nCores = 4
)
}
\arguments{
\item{subset_DT}{Data.frame with at least the following columns:
\describe{
\item{SNP}{SNP RSID}
\item{CHR}{chromosome}
\item{POS}{position}
}}

\item{foreground_filter}{Specify foreground by filtering SNPs in \code{subset_DT}.
Write filter as a string (or \code{NULL} to include all SNPs).}

\item{background_filter}{Specify background by filtering SNPs in \code{subset_DT}.
Write filter as a string (or \code{NULL} to include all SNPs).}
}
\description{
XGR uses a binomial enrichment tests for each annotation.
}
\examples{
\dontrun{
data("merged_DT")
enrich_res <- XGR.iterate_enrichment(subset_DT=merged_DT, foreground_filter = "Consensus_SNP", background_filter = "leadSNP", lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes"), nCores=1)
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{transethnic_plot}
\alias{transethnic_plot}
\title{Plot transethnic PAINTOR results}
\usage{
transethnic_plot(
  merged_dat,
  save_path,
  title = locus,
  subtitle = "Trans-ethnic Fine-mapping",
  PAINTOR.label = "PAINTOR\\nTrans-ethnic",
  conditions = c("MESA_AFA", "MESA_CAU", "MESA_HIS")
)
}
\description{
Plot transethnic PAINTOR results
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.get_promoter_celltypes}
\alias{NOTT_2019.get_promoter_celltypes}
\title{Get promoter cell types}
\usage{
NOTT_2019.get_promoter_celltypes(annot_sub, marker_key)
}
\description{
Brain cell-specific epigenomic data from Nott et al. (2019).
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{gene_locus_list}
\alias{gene_locus_list}
\title{Generate a named list of [e]gene-locus pairs}
\usage{
gene_locus_list(top_SNPs)
}
\description{
Generate a named list of [e]gene-locus pairs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloaders.R
\name{axel}
\alias{axel}
\title{axel}
\usage{
axel(
  input_url,
  output_path,
  background = F,
  nThread = 4,
  force_overwrite = F,
  quiet = T,
  alternate = T,
  check_certificates = F,
  conda_env = "echoR"
)
}
\description{
R wrapper for axel, which enables multi-threaded download of a single large file.
}
\seealso{
\url{https://github.com/axel-download-accelerator/axel/}

Other downloaders: 
\code{\link{downloader}()},
\code{\link{wget}()}
}
\concept{downloaders}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{query_by_coordinates_merged}
\alias{query_by_coordinates_merged}
\title{Use \emph{awk} to query locus subsets.}
\usage{
query_by_coordinates_merged(
  top_SNPs,
  fullSS_path,
  subset_path,
  locus,
  chrom_col,
  file_sep = " ",
  location_sep = ":",
  min_POS,
  max_POS,
  bp_distance
)
}
\arguments{
\item{top_SNPs}{A data.frame with the genomic coordinates of the lead SNP for each locus.
The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
Only one SNP per \strong{Locus} should be included.
At minimum, \code{top_SNPs} should include the following columns:
\describe{
\item{\emph{Locus}}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
\item{\emph{CHR}}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
\item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
}}

\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{file_sep}{The separator in the full summary stats file.
This parameter is only necessary if \code{query_by!="tabix"}.}

\item{location_sep}{The separator character when \strong{CHR} and \strong{POS} are merged into one column (e.g. ":" when formatted like chr12:12209944).}

\item{min_POS}{Manually set the minimum genomic position for your locus subset.
\code{min_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{min_POS=top_SNPs$min_POS}).}

\item{max_POS}{Manually set the maximum genomic position for your locus subset.
\code{max_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{max_POS=top_SNPs$max_POS}).}

\item{bp_distance}{The width of the window size you want each locus to be.
For example, if \code{bp_distance=500000} then the locus will span 500kb from the lead SNP in either direction,
resulting in a locus that is ~1Mb long (depending on the dataset).}
}
\description{
Search full summary stats file by genomics coordinates.
To be used when \strong{CHR} and \strong{POS} are merged into one column (e.g. chr12:12209944).
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FINEMAP.R
\name{FINEMAP.process_results}
\alias{FINEMAP.process_results}
\title{Post-processing of \code{FINEMAP} results}
\source{
\url{http://www.christianbenner.com}
}
\usage{
FINEMAP.process_results(
  locus_dir,
  subset_DT,
  credset_thresh = 0.95,
  pvalue_thresh = 0.05,
  finemap_version = "1.4",
  results_file = ".cred",
  n_causal = 5,
  nThread = 1,
  sort_by_CS = T,
  verbose = T
)
}
\description{
Post-processing of \code{FINEMAP} results
}
\examples{
\dontrun{
data("locus_dir"); data("BST1");
finemap_DT <- BST1
subset_DT <-FINEMAP.process_results(locus_dir=locus_dir, subset_DT=finemap_DT)
}
}
\seealso{
Other FINEMAP: 
\code{\link{FINEMAP.construct_data}()},
\code{\link{FINEMAP.construct_master}()},
\code{\link{FINEMAP.find_executable}()},
\code{\link{FINEMAP}()}
}
\concept{FINEMAP}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.merge_results}
\alias{PAINTOR.merge_results}
\title{Merge PAINTOR results}
\usage{
PAINTOR.merge_results(
  finemap_dat,
  paintor.results,
  PP_threshold = 0.5,
  multi_finemap_col_name = "PAINTOR"
)
}
\description{
@keywords internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveSparse.R
\name{saveSparse}
\alias{saveSparse}
\title{Save LD matrix as a sparse matrix}
\usage{
saveSparse(LD_matrix, LD_path, verbose = T)
}
\description{
Converting LD matrices to sparse format reduces file size by half.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.R
\name{standardize_subset}
\alias{standardize_subset}
\title{Standardize the locus subset}
\usage{
standardize_subset(
  locus,
  top_SNPs = NULL,
  fullSS_genome_build = "hg19",
  subset_path = "./Data",
  chrom_col = "CHR",
  position_col = "POS",
  snp_col = "SNP",
  pval_col = "P",
  effect_col = "Effect",
  stderr_col = "StdErr",
  tstat_col = "t_stat",
  MAF_col = "MAF",
  freq_col = "Freq",
  N_cases_col = "N_cases",
  N_controls_col = "N_controls",
  N_cases = NULL,
  N_controls = NULL,
  proportion_cases = "calculate",
  sample_size = NULL,
  A1_col = "A1",
  A2_col = "A2",
  gene_col = "Gene",
  QTL_prefixes = NULL,
  return_dt = T,
  nThread = 4,
  download_method = "axel",
  verbose = T
)
}
\arguments{
\item{top_SNPs}{A data.frame with the genomic coordinates of the lead SNP for each locus.
The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
Only one SNP per \strong{Locus} should be included.
At minimum, \code{top_SNPs} should include the following columns:
\describe{
\item{\emph{Locus}}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
\item{\emph{CHR}}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
\item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
}}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{snp_col}{Name of the SNP RSID column in the full summary stats file.
(\emph{default: ="SNP"})}

\item{pval_col}{Name of the p-value column in the full summary stats file.
Raw p-values are preferred, but if not available corrected p-values (e.g. FDR) can be used instead.
(\emph{default: ="P"})}

\item{effect_col}{Name of the effect size column in the full summary stats file.
Effect size is preferred, but if not available other metrics like Beta for Odds Ratio can be used instead.
(\emph{default: ="Effect"})}

\item{stderr_col}{Name of the standard error  column in the full summary stats file.
You can also set \code{stderr_col="calculate"} to infer standard error using: \code{effect / tstat}.
(\emph{default: ="StdErr"})}

\item{tstat_col}{Name of the t-statistic column in the full summary stats file.
This column is not necessary unless \code{stderr_col="calculate"} or the standard error column is missing.
(\emph{default: ="t-stat"})}

\item{MAF_col}{Name of the minor allele frequency column in the full summary stats file.
Can be inferred from \strong{freq_col} if missing from the dataset.
(\emph{default: ="MAF"})}

\item{freq_col}{Name of the allele frequency column in the full summary stats file.
Effect allele frequency is preferred, but the non-effect allele can be provided instead (though this may be less accurate).
This column is not necessary unless \code{MAF_col="calculate"} or the MAF column is missing.
(\emph{default: ="Freq"})}

\item{N_cases_col}{Name of the column in the full summary stats that has the number of case subjects in the study.
This can either be per SNP sample sizes, or one number repeated across all rows.
Proxy cases (e.g. relatives of people with the disease being investigated) should be included in this estimate if any were used in the study.
This column is not necesssary if \code{N_cases} parameter is provided.
(\emph{default: ="N_cases"})}

\item{N_controls_col}{Name of the column in the full summary stats that has the number of control subjects in the study.
 This can either be per SNP sample sizes, or one number repeated across all rows.
 This column is not necesssary if \code{N_controls} parameter is provided.
(\emph{default: ="N_controls"})}

\item{N_cases}{The number of case subjects in the study.
Instead of providing a redundant \strong{N_cases_col} column, you can simply enter one value here.}

\item{N_controls}{The number of control subjects in the study.
Instead of providing a redundant \strong{N_controls_col} column, you can simply enter one value here.}

\item{proportion_cases}{The proportion of total subjects in the study that were cases.
if \code{proportion_cases="calculate"} then this is inferred:  \code{N_controls / N_controls}.}

\item{sample_size}{The overall sample size of the study.
If none is given, and \strong{N_cases} and \strong{N_controls} columns are present,
then sample_size is inferred to be:  \code{max(N_cases) + max(N_controls)}.}

\item{A1_col}{Name of the effect/risk allele column in the full summary stats.
 \strong{\emph{IMPORTANT}}: Make sure this actually the case for your full summary stats file.
Unfortunately, different studies report different kinds of allele information in a non-standardized way.
Meaning that A1/A2 can refer to any number of things:
 \describe{
 \item{effect/other alleles}{in the case of diseases}
 \item{ref/alt alleles}{where ref is the reference genome being used}
 \item{major/minor alleles}{This dichotomy holds true for bi-allelic SNPs but not necessary multi-allelic SNPs}
 }
 This makes comparing summary stats across GWAS/QTL datasets very confusing for several reasons:
 \describe{
 \item{Multi-allelic SNPs}{SNPs can have more than just 2 possible alleles (multi-allelic SNPs). Even if you compare the same SNP between two studies, you may accidentally be comparing totally different alleles.}
 \item{Valence}{The valence (+/-) of per-SNP GWAS effect sizes/beta can be relative to different allele types between studies.
 For example, let's say in one GWAS study your effect size for SNP A is 1.5 relative to the major allele in one study,
  and the minor allele happens to be the one found in the reference genome.
  You then try to compare that effect size to that of the same SNP in another GWAS.
  But, the valence of the effect sizes in the 2nd GWAS study are all relative to the reference genome (instead of the minor allele),
  giving the same SNP a value of -1.2. If you took the effect sizes at face value you'd say the signals are in opposite directions.
  But once you take into account how the valences were determined in each study you realize that they're actually both positive relative to the major allele.}
 }
This process of reversing per-SNP valences based on aligning the alleles is known as allele flipping.
This is important when comparing individual SNPs, but can also have an impact on colocalization results.}

\item{gene_col}{For QTL studies, the name of the [e]gene column in the full summary stats file (\emph{default: "gene"}).
This column will be used for filtering summary stats if supplying a named list of gene:Locus pairs to \code{loci}.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
After querying a subset of the full summary statistics,
this function converts it into a standardized format
that the rest of \emph{echolocatoR} can work with.
}
\examples{
data("BST1")
... Screw up Freq to see if function can fix it and infer MAF ...
BST1$rsid <- BST1$SNP
BST1 <- data.frame(BST1)[,!colnames(BST1) \%in\% c("MAF","SNP")]
BST1[c(10,30,55),"Freq"] <- 0
BST1[c(12,22),"Freq"] <- NA
data.table::fwrite(BST1, "~/Desktop/results/GWAS/Nalls23andMe_2019/BST1/BST1.tsv")
query_mod <- standardize_subset(locus="BST1", subset_path="~/Desktop/results/GWAS/Nalls23andMe_2019/BST1/BST1.tsv", MAF_col="calculate", snp_col="rsid")
}
\seealso{
Other standardizing functions: 
\code{\link{get_UKB_MAF}()}
}
\concept{standardizing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.calculate_LD}
\alias{LD.calculate_LD}
\title{Calculate LD}
\usage{
LD.calculate_LD(
  locus_dir,
  ld_window = 1000,
  ld_format = "r",
  plink_prefix = "plink",
  verbose = T
)
}
\arguments{
\item{locus_dir}{Locus-specific results directory.}

\item{ld_window}{Set --r/--r2 max variant ct pairwise distance (usu. 10).}

\item{ld_format}{Whether to produce an LD matrix with
r (\code{ld_format="r"}) or D' (\code{ld_format="D"}) as the pairwise SNP correlation metric.}
}
\description{
Calculate a pairwise LD matrix from a vcf file using \emph{plink}.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{epigenetics_enrichment}
\alias{epigenetics_enrichment}
\title{Test for enrichment of \code{HaploR} annotations}
\source{
\href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
}
\usage{
epigenetics_enrichment(
  snp_list1,
  snp_list2,
  chisq = T,
  fisher = T,
  epigenetic_variables = c("Promoter_histone_marks", "Enhancer_histone_marks"),
  tissue_list = c("BRN", "BLD")
)
}
\description{
Test for enrichment of \code{HaploR} annotations
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{extract_SNP_subset}
\alias{extract_SNP_subset}
\title{Extract a subset of the summary stats}
\usage{
extract_SNP_subset(
  locus = NULL,
  locus_dir,
  results_dir = NULL,
  fullSS_path,
  fullSS_genome_build = "hg19",
  subset_path,
  LD_reference,
  force_new_subset = F,
  top_SNPs = "auto",
  bp_distance = 5e+05,
  chrom_col = "CHR",
  chrom_type = NULL,
  position_col = "POS",
  snp_col = "SNP",
  locus_col = "Locus",
  pval_col = "P",
  effect_col = "Effect",
  stderr_col = "StdErr",
  MAF_col = "MAF",
  freq_col = "Freq",
  tstat_col = "t-stat",
  A1_col = "A1",
  A2_col = "A2",
  gene_col = "gene",
  N_cases_col = "N_cases",
  N_controls_col = "N_controls",
  N_cases = NULL,
  N_controls = NULL,
  proportion_cases = "calculate",
  sample_size = NULL,
  superpopulation = "",
  min_POS = NA,
  max_POS = NA,
  genes_detected = F,
  file_sep = "\\t",
  query_by = "coordinates",
  probe_path = "./Data/eQTL/gene.ILMN.map",
  QTL_prefixes = NULL,
  remove_tmps = T,
  conda_env = "echoR",
  verbose = T
)
}
\description{
Use either \emph{tabix} or \emph{awk} to extract a locus subset
 from the full summary statistics file.
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAIN.R
\name{finemap_loci}
\alias{finemap_loci}
\title{Fine-map multiple loci}
\usage{
finemap_loci(
  loci,
  fullSS_path,
  fullSS_genome_build = NULL,
  LD_genome_build = "hg19",
  dataset_name = "dataset_name",
  dataset_type = "GWAS",
  force_new_subset = F,
  force_new_LD = F,
  force_new_finemap = T,
  results_dir = "./results",
  top_SNPs = "auto",
  finemap_methods = c("ABF", "FINEMAP", "SUSIE", "POLYFUN_SUSIE"),
  finemap_args = NULL,
  bp_distance = 5e+05,
  n_causal = 5,
  munged = FALSE,
  chrom_col = "CHR",
  chrom_type = NULL,
  position_col = "POS",
  snp_col = "SNP",
  pval_col = "P",
  effect_col = "Effect",
  stderr_col = "StdErr",
  tstat_col = "t_stat",
  MAF_col = "MAF",
  locus_col = "Locus",
  freq_col = "Freq",
  A1_col = "A1",
  A2_col = "A2",
  gene_col = "Gene",
  N_cases_col = "N_cases",
  N_controls_col = "N_controls",
  N_cases = NULL,
  N_controls = NULL,
  proportion_cases = "calculate",
  sample_size = NULL,
  LD_reference = "1KGphase1",
  superpopulation = "EUR",
  download_method = "axel",
  vcf_folder = NULL,
  remote_LD = T,
  topVariants = 3,
  min_POS = NA,
  max_POS = NA,
  min_MAF = NA,
  trim_gene_limits = F,
  max_snps = NULL,
  file_sep = "\\t",
  min_r2 = 0,
  LD_block = F,
  LD_block_size = 0.7,
  query_by = "tabix",
  remove_variants = F,
  remove_correlates = F,
  probe_path = "./Data/eQTL/gene.ILMN.map",
  conditioned_snps = "auto",
  plot_LD = F,
  remove_tmps = T,
  PAINTOR_QTL_datasets = NULL,
  server = F,
  PP_threshold = 0.95,
  consensus_threshold = 2,
  case_control = T,
  QTL_prefixes = NULL,
  plot.types = c("simple"),
  plot.zoom = "1x",
  plot.Nott_epigenome = F,
  plot.Nott_show_placseq = F,
  plot.Nott_binwidth = 200,
  plot.Nott_bigwig_dir = NULL,
  plot.XGR_libnames = NULL,
  plot.Roadmap = F,
  plot.Roadmap_query = NULL,
  conda_env = "echoR",
  nThread = 4,
  verbose = T
)
}
\arguments{
\item{loci}{The list of loci you want to fine-map.
If \code{subset_path="auto"} (\emph{default}), a locus subset file name is automatically constructed as:
\emph{Data/<dataset_type>/<dataset_name>/<locus>/Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}}

\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{dataset_name}{The name you want to assign to the dataset being fine-mapped,
This will be used to name the subdirectory where your results will be stored
(e.g. \emph{Data/GWAS/<dataset_name>}).
Don't use special characters (e.g.".", "/").}

\item{dataset_type}{The kind dataset you're fine-mapping (e.g. GWAS, eQTL, tQTL).
This will also be used when creating the subdirectory where your results will be stored
(e.g. \emph{Data/<dataset_type>/Kunkle_2019}).}

\item{force_new_subset}{By default, if a subset of the full summary stats file for a given locus is already present,
then \pkg{echolocatoR} will just use the pre-existing file.
Set \code{force_new_subset=T} to override this and extract a new subset.
Subsets are saved in the following path structure:
\emph{Data/<dataset_type>/<dataset_name>/<locus>/Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}}

\item{force_new_LD}{By default, if an LD matrix file for a given locus is already present,
then \pkg{echolocatoR} will just use the preexisting file.
Set \code{force_new_LD=T} to override this and extract a new subset.}

\item{force_new_finemap}{By default, if an fine-mapping results file for a given locus is already present,
then \pkg{echolocatoR} will just use the preexisting file.
Set \code{force_new_finemap=T} to override this and re-run fine-mapping.}

\item{results_dir}{Where to store all results.
\strong{IMPORTANT!:} It is usually best to provide the absolute path rather than the relative path.
This is especially important for \emph{FINEMAP}.}

\item{top_SNPs}{A data.frame with the genomic coordinates of the lead SNP for each locus.
The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
Only one SNP per \strong{Locus} should be included.
At minimum, \code{top_SNPs} should include the following columns:
\describe{
\item{\emph{Locus}}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
\item{\emph{CHR}}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
\item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
}}

\item{finemap_methods}{Which fine-mapping methods you want to use.}

\item{bp_distance}{The width of the window size you want each locus to be.
For example, if \code{bp_distance=500000} then the locus will span 500kb from the lead SNP in either direction,
resulting in a locus that is ~1Mb long (depending on the dataset).}

\item{n_causal}{The maximum number of potential causal SNPs per locus.
This parameter is used somewhat differntly by different fine-mapping tools.
See tool-specific functions for details.}

\item{munged}{Whether \code{fullSS_path} have already been standardised/filtered  full summary stats
with \code{MungeSumstats::format_sumstats}.
If \code{munged=FALSE} you'll need to provide the necessary column name arguments.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{snp_col}{Name of the SNP RSID column in the full summary stats file.
(\emph{default: ="SNP"})}

\item{pval_col}{Name of the p-value column in the full summary stats file.
Raw p-values are preferred, but if not available corrected p-values (e.g. FDR) can be used instead.
(\emph{default: ="P"})}

\item{effect_col}{Name of the effect size column in the full summary stats file.
Effect size is preferred, but if not available other metrics like Beta for Odds Ratio can be used instead.
(\emph{default: ="Effect"})}

\item{stderr_col}{Name of the standard error  column in the full summary stats file.
You can also set \code{stderr_col="calculate"} to infer standard error using: \code{effect / tstat}.
(\emph{default: ="StdErr"})}

\item{tstat_col}{Name of the t-statistic column in the full summary stats file.
This column is not necessary unless \code{stderr_col="calculate"} or the standard error column is missing.
(\emph{default: ="t-stat"})}

\item{MAF_col}{Name of the minor allele frequency column in the full summary stats file.
Can be inferred from \strong{freq_col} if missing from the dataset.
(\emph{default: ="MAF"})}

\item{locus_col}{Name of the locus column in the full summary stats file.
(\emph{default: ="Locus"})}

\item{freq_col}{Name of the allele frequency column in the full summary stats file.
Effect allele frequency is preferred, but the non-effect allele can be provided instead (though this may be less accurate).
This column is not necessary unless \code{MAF_col="calculate"} or the MAF column is missing.
(\emph{default: ="Freq"})}

\item{A1_col}{Name of the effect/risk allele column in the full summary stats.
 \strong{\emph{IMPORTANT}}: Make sure this actually the case for your full summary stats file.
Unfortunately, different studies report different kinds of allele information in a non-standardized way.
Meaning that A1/A2 can refer to any number of things:
 \describe{
 \item{effect/other alleles}{in the case of diseases}
 \item{ref/alt alleles}{where ref is the reference genome being used}
 \item{major/minor alleles}{This dichotomy holds true for bi-allelic SNPs but not necessary multi-allelic SNPs}
 }
 This makes comparing summary stats across GWAS/QTL datasets very confusing for several reasons:
 \describe{
 \item{Multi-allelic SNPs}{SNPs can have more than just 2 possible alleles (multi-allelic SNPs). Even if you compare the same SNP between two studies, you may accidentally be comparing totally different alleles.}
 \item{Valence}{The valence (+/-) of per-SNP GWAS effect sizes/beta can be relative to different allele types between studies.
 For example, let's say in one GWAS study your effect size for SNP A is 1.5 relative to the major allele in one study,
  and the minor allele happens to be the one found in the reference genome.
  You then try to compare that effect size to that of the same SNP in another GWAS.
  But, the valence of the effect sizes in the 2nd GWAS study are all relative to the reference genome (instead of the minor allele),
  giving the same SNP a value of -1.2. If you took the effect sizes at face value you'd say the signals are in opposite directions.
  But once you take into account how the valences were determined in each study you realize that they're actually both positive relative to the major allele.}
 }
This process of reversing per-SNP valences based on aligning the alleles is known as allele flipping.
This is important when comparing individual SNPs, but can also have an impact on colocalization results.}

\item{gene_col}{For QTL studies, the name of the [e]gene column in the full summary stats file (\emph{default: "gene"}).
This column will be used for filtering summary stats if supplying a named list of gene:Locus pairs to \code{loci}.}

\item{N_cases_col}{Name of the column in the full summary stats that has the number of case subjects in the study.
This can either be per SNP sample sizes, or one number repeated across all rows.
Proxy cases (e.g. relatives of people with the disease being investigated) should be included in this estimate if any were used in the study.
This column is not necesssary if \code{N_cases} parameter is provided.
(\emph{default: ="N_cases"})}

\item{N_controls_col}{Name of the column in the full summary stats that has the number of control subjects in the study.
 This can either be per SNP sample sizes, or one number repeated across all rows.
 This column is not necesssary if \code{N_controls} parameter is provided.
(\emph{default: ="N_controls"})}

\item{N_cases}{The number of case subjects in the study.
Instead of providing a redundant \strong{N_cases_col} column, you can simply enter one value here.}

\item{N_controls}{The number of control subjects in the study.
Instead of providing a redundant \strong{N_controls_col} column, you can simply enter one value here.}

\item{proportion_cases}{The proportion of total subjects in the study that were cases.
if \code{proportion_cases="calculate"} then this is inferred:  \code{N_controls / N_controls}.}

\item{sample_size}{The overall sample size of the study.
If none is given, and \strong{N_cases} and \strong{N_controls} columns are present,
then sample_size is inferred to be:  \code{max(N_cases) + max(N_controls)}.}

\item{LD_reference}{Which linkage disequilibrium reference panel do you want to use.
Options include:
\describe{
\item{"UKB"}{A pre-caclulated LD reference matrix from a subset of caucasian British individuals from the UK Biobank. See \href{https://www.biorxiv.org/content/10.1101/807792v2}{Wiessbrod et al. (2019)} for more details.}
\item{"1KGphase1"}{Download a subset of the 1000 Genomes Project Phase 1 vcf and calculate LD on the fly with plink.}
\item{"1KGphase3"}{Download a subset of the 1000 Genomes Project Phase 3 vcf and calculate LD on the fly with plink.}
\item{"<path>/*.vcf" or "<path>/*.vcf.gz"}{Alternatively, users can provide their own custom panel by supplying a list of \emph{.vcf} file path (one per locus) which \pkg{echolocatoR} will use to compute LD (using \emph{plink}).}
}}

\item{superpopulation}{Subset your LD reference panel by superopulation.
Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
\href{https://www.internationalgenome.org/faq/which-populations-are-part-your-study/}{1KGphase1 options} include:
\describe{
\item{"AFR"}{African [descent]}
\item{"AMR"}{Ad-mixed American}
\item{"EAS"}{East Asian}
\item{"EUR"}{European}
\item{"SAS"}{South Asian}
}}

\item{remote_LD}{When acquiring LD matrices,
the default is to delete the full vcf or npz files after \pkg{echolocatoR} has extracted the necssary subset.
However, if you wish to keep these full files (which can be quite large) set \code{remote_LD=T}.}

\item{min_POS}{Manually set the minimum genomic position for your locus subset.
\code{min_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{min_POS=top_SNPs$min_POS}).}

\item{max_POS}{Manually set the maximum genomic position for your locus subset.
\code{max_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{max_POS=top_SNPs$max_POS}).}

\item{min_MAF}{Remove any SNPs with \strong{MAF} < \code{min_MAF}.}

\item{trim_gene_limits}{If a valid gene symbol is provided to \code{trim_gene_limits},
the gene's canonical coordinates are pulled from \code{biomaRt}.
This includes introns, exons, and proximal regulatory regions (e.g. promoters).
Any SNPs that fall outside these coordinates are remove from downstream fine-mapping.
Set \code{trim_gene_limits=F} to not limit by gene coordinates (\emph{default}).}

\item{max_snps}{The maximum number of SNPs to include in the locus.
If the current window size yields > \code{max_snps},
 then the outer edges of the of the locus are trimmed until the number of SNPs ≤ \code{max_snps}.}

\item{file_sep}{The separator in the full summary stats file.
This parameter is only necessary if \code{query_by!="tabix"}.}

\item{min_r2}{Remove any SNPs are below the LD r2 threshold with the lead SNP within their respective locus.}

\item{LD_block}{Calculate LD blocks with \emph{plink} and only include the block to which the lead SNP belongs.}

\item{LD_block_size}{Adjust the granularity of block sizes when \code{LD_block=T}.}

\item{query_by}{Choose which method you want to use to extract locus subsets from the full summary stats file.
Methods include:
\describe{
\item{"tabix"}{Convert the full summary stats file in an indexed tabix file. Makes querying lightning fast after the initial conversion is done. (\emph{default})}
\item{"coordinates"}{Extract locus subsets using min/max genomic coordinates with \emph{awk}.}
}}

\item{remove_variants}{A list of variants to remove from the locus subset file.}

\item{remove_correlates}{A named list, where the names are the RSIDs of SNPs
whose LD correlates you wish to remove,
and the value is the absolute r2 threshold you wish to filter at for each RSID respectively
(e.g. \code{ remove_correlates = c("rs76904798"=.2, "rs10000737"=.8)}).
This will also remove the SNPs in \code{remove_correlates} themselves.}

\item{probe_path}{The location of the file containing translations between probe IDs and gene symbols.
Only used for certain eQTL datasets.}

\item{conditioned_snps}{Which SNPs to conditions on when fine-mapping with \emph{COJO}.}

\item{plot_LD}{Whether to plot a subset of the LD matix.}

\item{remove_tmps}{Whether to remove any temporary files (e.g. FINEMAP output files) after the pipeline is done running.}

\item{PAINTOR_QTL_datasets}{A list of QTL datasets to be used when conducting joint functional fine-mapping with \emph{PAINTOR}.}

\item{server}{Whether \pkg{echolocatoR} is being run on a computing cluster/server or on a local machine.}

\item{PP_threshold}{The minimum fine-mapped posterior probability for a SNP to be considered part of a Credible Set.
For example, \code{PP_threshold=.95} means that all Credible Set SNPs will be 95\% Credible Set SNPs.}

\item{consensus_threshold}{The minimum number of fine-mapping tools that include a SNP
in their 95\% Credible Sets to consider that it a "Consensus SNP" (\emph{default=2}).}

\item{plot.types}{Which kinds of plots to include.
Options:
\describe{
\item{"simple"}{Just plot the following tracks: GWAS, fine-mapping, gene models}
\item{"fancy"}{Additionally plot XGR annotation tracks (XGR, Roadmap, Nott).}
}}

\item{plot.zoom}{Zoom into the center of the locus when plotting (without editing the fine-mapping results file).
You can provide either:
\itemize{
\item{The size of your plot window in terms of basepairs (e.g. \code{plot.zoom=50000} for a 50kb window)}.
\item{How much you want to zoom in (e.g. \code{plot.zoom="1x"} for the full locus, \code{plot.zoom="2x"} for 2x zoom into the center of the locus, etc.)}.
}
You can pass a list of window sizes (e.g. \code{c(50000,100000,500000)}) to automatically generate
multiple views of each locus.
This can even be a mix of different style inputs: e.g. \code{c("1x","4.5x",25000)}.}

\item{plot.Nott_binwidth}{When including Nott et al. (2019) epigenomic data in the track plots,
adjust the bin width of the histograms.}

\item{plot.Nott_bigwig_dir}{Instead of pulling Nott et al. (2019) epigenomic data
from the \emph{UCSC Genome Browser}, use a set of local bigwig files.}

\item{plot.Roadmap}{Find and plot annotations from Roadmap.}

\item{plot.Roadmap_query}{Only plot annotations from Roadmap whose metadata contains a string or any items from  a list of strings
(e.g. \code{"brain"} or \code{c("brain","liver","monocytes")}).}

\item{conda_env}{The name of a conda environment to use.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\value{
A merged data.frame with all fine-mapping results from all loci.
}
\description{
\pkg{echolocatoR} will automatically fine-map each locus.
Uses the \code{top_SNPs} data.frame to define locus coordinates.
}
\examples{
\dontrun{
library(echolocatoR)
root.dir <- tempdir()

data("Nalls_top_SNPs");
top_SNPs <- import_topSNPs(
  topSS = Nalls_top_SNPs,
  chrom_col = "CHR",
  position_col = "BP",
  snp_col="SNP",
  pval_col="P, all studies",
  effect_col="Beta, all studies",
  gene_col="Nearest Gene",
  locus_col = "Nearest Gene",
  grouping_vars = c("Locus"),
  remove_variants = "rs34637584")
fullSS_path <- example_fullSS(fullSS_path=file.path(root.dir,"Nalls23andMe_2019.fullSS_subset.tsv") )

Nalls23andMe_2019.results <- finemap_loci(# GENERAL ARGUMENTS
  top_SNPs = top_SNPs,
  #  It's best to give absolute paths
  results_dir = file.path(root.dir,"results"),
  loci = c("BST1","MEX3C"),# top_SNPs$Locus,
  dataset_name = "Nalls23andMe_2019",
  dataset_type = "GWAS",
  force_new_subset = F,
  force_new_LD = F,
  force_new_finemap = T,
  remove_tmps = F,

  # SUMMARY STATS ARGUMENTS
  fullSS_path = fullSS_path,
  query_by ="tabix",
  chrom_col = "CHR", position_col = "POS", snp_col = "RSID",
  pval_col = "p", effect_col = "beta", stderr_col = "se",
  freq_col = "freq", MAF_col = "calculate",
  A1_col = "A1",
  A2_col = "A2",

  # FILTERING ARGUMENTS
  ## It's often desirable to use a larger window size
  ## (e.g. 2Mb which is bp_distance=500000*2),
  ## but we use a small window here to speed up the process.
  bp_distance = 10000,#500000*2,
  min_MAF = 0.001,
  trim_gene_limits = F,

  # FINE-MAPPING ARGUMENTS
  ## General
  finemap_methods = c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
  n_causal = 5,
  PP_threshold = .95,

  # LD ARGUMENTS
  LD_reference = "1KGphase1",#"UKB",
  superpopulation = "EUR",
  download_method = "axel",

  # PLOT ARGUMENTS
  ## general
  plot.types=c("fancy"),
  ## Generate multiple plots of different window sizes;
  ### all SNPs, 4x zoomed-in, and a 50000bp window
  plot.zoom = c("all","4x","10x"),
  ## XGR
  # plot.XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes"),
  ## Roadmap
  plot.Roadmap = F,
  plot.Roadmap_query = NULL,
  # Nott et al. (2019)
  plot.Nott_epigenome = T,
  plot.Nott_show_placseq = T,

  verbose = F
)

}
}
\seealso{
Other MAIN: 
\code{\link{finemap_pipeline}()}
}
\concept{MAIN}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Roadmap.R
\name{ROADMAP.construct_reference}
\alias{ROADMAP.construct_reference}
\title{Gather Roadmap annotation metadata}
\usage{
ROADMAP.construct_reference(
  ref_path = system.file("extdata/ROADMAP", "ROADMAP_Epigenomic.js", package =
    "echolocatoR"),
  keyword_query = NULL
)
}
\arguments{
\item{keyword_query}{Search all columns in the Roadmap annotations metadata
and only query annotations that contain your keywords.
Can provide multiple keywords in list form:
\code{c("placenta","liver","monocytes")}}
}
\description{
Gather Roadmap annotation metadata
}
\seealso{
Other ROADMAP: 
\code{\link{ROADMAP.merge_and_process_grl}()},
\code{\link{ROADMAP.query_and_plot}()},
\code{\link{ROADMAP.query}()},
\code{\link{ROADMAP.tabix}()}
}
\concept{ROADMAP}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.get_ref_prefix}
\alias{POLYFUN.get_ref_prefix}
\title{Get ref data path prefix}
\usage{
POLYFUN.get_ref_prefix(ref_prefix = NULL, locus_dir)
}
\description{
Get ref data path prefix
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dey_DeepLearning.R
\name{DEEPLEARNING.query}
\alias{DEEPLEARNING.query}
\title{Iteratively collect deep learning annotations}
\usage{
DEEPLEARNING.query(
  merged_dat,
  base_url = "/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
  level = c("Variant_Level", "Allelic_Effect"),
  tissue = c("NTS", "Blood", "Brain"),
  model = c("Basenji", "BiClassCNN", "DeepSEA", "ChromHMM", "Roadmap", "Others"),
  mean_max = c("MEAN", "MAX"),
  type = c("annot", "ldscore"),
  nThread = 4,
  verbose = T
)
}
\description{
Iteratively collect deep learning annotations
}
\examples{
\dontrun{
root = "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
merged_dat <- merge_finemapping_results(dataset = "Data/GWAS/Nalls23andMe_2019", LD_reference = "UKB", minimum_support = 0)
merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)

#### Allelic_Effect ####
ANNOT.ae <- DEEPLEARNING.query (merged_dat=merged_dat, level="Allelic_Effect", type="annot")
data.table::fwrite(ANNOT.ae, file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.csv.gz"))


#### Variant_Level ####
ANNOT.vl <- DEEPLEARNING.query (merged_dat=merged_dat, level="Variant_Level", type="annot")
data.table::fwrite(ANNOT.vl, file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.csv.gz"))

}
}
\seealso{
Other DEEPLEARNING: 
\code{\link{DEEPLEARNING.query_multi_chr}()},
\code{\link{DEEPLEARNING.query_one_chr}()}
}
\concept{DEEPLEARNING}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_lead_SNP.R
\name{assign_lead_SNP}
\alias{assign_lead_SNP}
\title{Assign a lead GWAS SNP to a locus}
\usage{
assign_lead_SNP(new_DT, verbose = T)
}
\arguments{
\item{data.frame}{Fine-mapping results data.frame.}
}
\value{
Fine-mapping results data.frame with new boolean \strong{leadSNP} column,
 indicating whether each SNPs is the lead GWAS SNP in that locus or not.
}
\description{
If none of the SNPs in the data.frame have \code{leadSNP==T},
then sort by lowest p-value (and then highest Effect size) and assign the top SNP as the lead SNP.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{VALIDATION.bootstrap_multimetric}
\alias{VALIDATION.bootstrap_multimetric}
\title{Conduct bootstrap procedure on multiple columns}
\usage{
VALIDATION.bootstrap_multimetric(
  metric_df,
  metric_names,
  validation_method = NULL,
  locus_means = T,
  grouping_var = NULL,
  test_method = "coin_wilcox_test",
  iterations = 1000,
  nThread = 4,
  save_path = F,
  verbose = T
)
}
\description{
Conduct bootstrap procedure on multiple columns
}
\examples{
\dontrun{
save_path <- root <-  "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"


#### SURE MPRA #####
## mean
path <- file.path(root,"SURE/Nalls23andMe_2019.SURE.snp_groups.mean.csv.gz")
metric_names <- "p"
## raw
path <- file.path(root,"SURE/Nalls23andMe_2019.SURE.csv.gz")
metric_names <- c("k562.wilcox.p.value","hepg2.wilcox.p.value")
validation_method <- "SuRE MPRA"

#### Dey_DeepLearning ####
## mean
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.snp_groups_mean.csv.gz")
metric_names <- "value"; grouping_var="annot"
## raw
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.csv.gz")
## path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.csv.gz")
validation_method = "Dey_DeepLearning"


# -----Import and process ------
metric_df <- data.table::fread(path, nThread=8)
## metric_names <- grep("Basenji.*MAX|DeepSEA.*MAX|Roadmap.*MAX", colnames(metric_df), value = T)
## metric_df <- metric_df \%>\% dplyr::mutate(annot=paste(Model,Tissue,Assay,Type,Metric,sep="_"))
boot_res <- VALIDATION.bootstrap_multimetric(metric_df=metric_df, metric_names=metric_names, validation_method=validation_method, save_path=gsub("\\\\.csv\\\\.gz",".bootstrap.coin_wilcox_test.csv.gz",path), grouping_var=grouping_var, iterations=10000, nThread=12)
}
}
\seealso{
Other VALIDATION: 
\code{\link{VALIDATION.bootstrap_plot}()},
\code{\link{VALIDATION.bootstrap}()},
\code{\link{VALIDATION.compare_bootstrap_distributions}()},
\code{\link{VALIDATION.permute}()},
\code{\link{VALIDATION.super_plot}()}
}
\concept{VALIDATION}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{merge_finemapping_results}
\alias{merge_finemapping_results}
\title{Merge fine-mapping results from all loci}
\usage{
merge_finemapping_results(
  dataset = "./Data/GWAS",
  minimum_support = 1,
  include_leadSNPs = T,
  LD_reference = NULL,
  save_path = F,
  from_storage = T,
  haploreg_annotation = F,
  regulomeDB_annotation = F,
  biomart_annotation = F,
  PP_threshold = 0.95,
  consensus_threshold = 2,
  exclude_methods = NULL,
  top_CS_only = F,
  verbose = T,
  nThread = 4
)
}
\arguments{
\item{dataset}{Path to the folder you want to recursively search for results files within
 (e.g. \url{"Data/GWAS/Nalls23andMe_2019"}).
Set this to a path that includes multiple subfolders if you want to gather results
from multiple studies at once
(e.g. \url{"Data/GWAS"}).}

\item{minimum_support}{Filter SNPs by the minimum number
of fine-mapping tools that contained the SNP in their Credible Set.}

\item{include_leadSNPs}{Include lead GWAS/QTL SNPs per locus
(regardless of other filtering criterion).}

\item{from_storage}{Search for stored results files.}

\item{haploreg_annotation}{Annotate SNPs with HaploReg (using \code{HaploR}).}

\item{regulomeDB_annotation}{Annotate SNPs with regulaomeDB (using \code{HaploR}).}

\item{biomart_annotation}{Annotate SNPs with \code{biomart}.}

\item{PP_threshold}{Mean posterior probability threshold to include SNPs in mean PP Credible Set
(averaged across all fine-mapping tools).}

\item{exclude_methods}{Exclude certain fine-mapping methods when estimating
\strong{mean.CS} and \strong{Consensus_SNP}.}

\item{verbose}{Print messages.}

\item{xlsx_path}{Save merged data.frame as excel file.}

\item{consensus_thresh}{The minimum number of tools that have the SNPs in their Credible Set
to classify it as a \strong{Consensus_SNP}.}
}
\description{
Gather fine-mapping results from \emph{echolocatoR} across all loci
and merge into a single data.frame.
}
\examples{
dataset_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
# UCS and lead SNPs: No annotation
merged_dat <- merge_finemapping_results(dataset=dataset_dir, minimum_support=1, include_leadSNPs=T)

# UCS and lead SNPs: With annotations
merged_dat <- merge_finemapping_results(dataset=dataset_dir, minimum_support=1, include_leadSNPs=T, haploreg_annotation=T, biomart_annotation=T)
}
\concept{annotatate}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_snps.R
\name{filter_snps}
\alias{filter_snps}
\title{Filter SNPs}
\usage{
filter_snps(
  subset_DT,
  bp_distance = 5e+05,
  remove_variants = F,
  min_POS = NA,
  max_POS = NA,
  max_snps = NULL,
  min_MAF = NULL,
  trim_gene_limits = F,
  verbose = T
)
}
\description{
Filter SNps by MAF, window size, min/max position, maxmimum number of SNPs, or gene coordinates.
You can also explicitly remove certain variants.
}
\examples{
data("BST1");
subset_DT <- filter_snps(subset_DT=BST1)
}
\seealso{
Other SNP filters: 
\code{\link{gene_trimmer}()},
\code{\link{limit_SNPs}()},
\code{\link{snps_to_condition}()},
\code{\link{subset_common_snps}()}
}
\concept{SNP filters}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.get_CS_bins}
\alias{SUMMARISE.get_CS_bins}
\title{Count bins of tool-specific and union CS sizes}
\usage{
SUMMARISE.get_CS_bins(merged_dat)
}
\description{
Count bins of tool-specific and union CS sizes
}
\examples{
data("merged_DT");
bin_counts <- SUMMARISE.get_CS_bins(merged_dat=merged_DT)
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{results_report}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.peak_overlap_plot}
\alias{SUMMARISE.peak_overlap_plot}
\title{Plot overlap between some SNP group and various epigenomic data}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (2020/bioRxiv)}
}
\usage{
SUMMARISE.peak_overlap_plot(
  merged_dat,
  snp_filter = "Consensus_SNP==T",
  include.NOTT_2019_peaks = T,
  include.NOTT_2019_enhancers_promoters = T,
  include.NOTT_2019_PLACseq = T,
  include.CORCES_2020_scATACpeaks = T,
  include.CORCES_2020_Cicero_coaccess = T,
  include.CORCES_2020_bulkATACpeaks = T,
  include.CORCES_2020_HiChIP_FitHiChIP_coaccess = T,
  include.CORCES_2020_gene_annotations = T,
  plot_celltype_specificity = T,
  plot_celltype_specificity_genes = F,
  facets_formula = ". ~ Cell_type",
  show_plot = T,
  label_yaxis = T,
  x_strip_angle = 90,
  x_tick_angle = 40,
  drop_empty_cols = F,
  fill_title = paste(snp_filter, "\\nin epigenomic peaks"),
  save_path = F,
  height = 11,
  width = 12,
  subplot_widths = c(1, 0.5),
  verbose = T
)
}
\arguments{
\item{include.NOTT_2019_peaks}{Plot SNP subset overlap with
peaks from cell-type-specific bulk ATAC, H3K27ac, and H3K4me3 assays.}

\item{include.NOTT_2019_enhancers_promoters}{Plot SNP subset overlap with
cell enhancers and promoters.}

\item{include.CORCES_2020_scATACpeaks}{Plot SNP subset overlap with
cell-type-specific scATAC-seq peaks.}

\item{include.CORCES_2020_Cicero_coaccess}{Plot SNP subset overlap with
Cicero coaccessibility peaks (derived from scATACseq).}
}
\description{
Plot overlap between some SNP group and various epigenomic data
}
\examples{
data("merged_DT");

... Consensus SNPs ...
gg_peaks <- SUMMARISE.peak_overlap_plot(merged_dat=merged_DT, snp_filter="Consensus_SNP==T", fill_title="Consensus SNPs in epigenomic peaks")
... UCS SNPs ...
gg_peaks <- SUMMARISE.peak_overlap_plot(merged_dat=merged_DT, snp_filter="Support>0", fill_title="UCS SNPs in epigenomic peaks")
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_CS_counts}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{results_report}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_locus_dir}
\alias{get_locus_dir}
\title{Extract the locus dir}
\usage{
get_locus_dir(subset_path)
}
\description{
Extract the locus dir
}
\seealso{
Other directory functions: 
\code{\link{get_multifinemap_path}()},
\code{\link{get_study_dir}()},
\code{\link{get_subset_path}()},
\code{\link{make_locus_dir}()}
}
\concept{directory functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.translate_population}
\alias{LD.translate_population}
\title{Translate superopulation acronyms}
\usage{
LD.translate_population(superpopulation)
}
\arguments{
\item{superpopulation}{Three-letter superpopulation name.}
}
\description{
Ensures a common ontology for synonynmous superpopulation names.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LDlinkR.LDproxy_batch}
\alias{LDlinkR.LDproxy_batch}
\title{Extract LD proxies from 1KGphase3}
\source{
\href{https://www.rdocumentation.org/packages/LDlinkR/versions/1.0.2}{website}
}
\usage{
LDlinkR.LDproxy_batch(
  snp,
  pop = "CEU",
  r2d = "r2",
  min_corr = F,
  save_dir = NULL,
  verbose = T
)
}
\description{
Wrapper for \code{LDlinkR::LDproxy_batch}.
Eeasy to use but doesn't scale up well to many SNPs (takes way too long).
}
\examples{
data("merged_DT")
lead.snps <- setNames(subset(merged_dat, leadSNP)$Locus, subset(merged_dat, leadSNP)$SNP)
proxies <- LDlinkR.LDproxy_batch(snp=lead.snps)
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_common_snps.R
\name{subset_common_snps}
\alias{subset_common_snps}
\title{Subset LD matrix and dataframe to only their shared SNPs}
\usage{
subset_common_snps(LD_matrix, finemap_dat, fillNA = 0, verbose = F)
}
\value{
data.frame
}
\description{
Find the SNPs that are shared between an LD matrix and another data.frame with a `SNP` column.
Then remove any non-shared SNPs from both objects.
}
\examples{
data("BST1"); data('LD_matrix');
finemap_dat=BST1
out <- subset_common_snps(LD_matrix=LD_matrix, finemap_dat=BST1)
}
\seealso{
Other SNP filters: 
\code{\link{filter_snps}()},
\code{\link{gene_trimmer}()},
\code{\link{limit_SNPs}()},
\code{\link{snps_to_condition}()}
}
\concept{SNP filters}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.R
\name{get_UKB_MAF}
\alias{get_UKB_MAF}
\title{Get MAF from UK Biobank.}
\usage{
get_UKB_MAF(
  subset_DT,
  output_path = "./Data/Reference/UKB_MAF",
  force_new_maf = F,
  download_method = "axel",
  nThread = 4,
  verbose = T,
  conda_env = "echoR"
)
}
\description{
If \strong{MAF} column is missing,
download MAF from UK Biobank and use that instead.
}
\examples{
\dontrun{
data("BST1");
subset_DT <- data.frame(BST1)[,colnames(BST1)!="MAF"]
BST1 <- get_UKB_MAF(subset_DT=subset_DT )
}
}
\seealso{
Other standardizing functions: 
\code{\link{standardize_subset}()}
}
\concept{standardizing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_trimmer.R
\name{gene_trimmer}
\alias{gene_trimmer}
\title{Remove all SNPs outside of of a given gene.}
\usage{
gene_trimmer(subset_DT, gene, min_POS = NULL, max_POS = NULL)
}
\description{
Get the min/max coordinates of a given gene (including known regulatory regions, introns, and exons).
Remove any SNPs from the data.frame that fall outside these coordinates.
}
\seealso{
Other SNP filters: 
\code{\link{filter_snps}()},
\code{\link{limit_SNPs}()},
\code{\link{snps_to_condition}()},
\code{\link{subset_common_snps}()}
}
\concept{SNP filters}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{GRanges_to_BED}
\alias{GRanges_to_BED}
\title{Convert GRanges object to BED format and save}
\usage{
GRanges_to_BED(GR.annotations, output_path, sep = "\\t", nThread = 4, gzip = F)
}
\description{
Convert GRanges object to BED format and save
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FINEMAP.R
\name{FINEMAP.construct_data}
\alias{FINEMAP.construct_data}
\title{Prepare input files for \code{FINEMAP}}
\source{
\url{http://www.christianbenner.com}
}
\usage{
FINEMAP.construct_data(
  locus_dir,
  subset_DT,
  LD_matrix,
  nThread = 4,
  verbose = T
)
}
\description{
Creates and saves 1) the summary stats file, and 2) the LD matrix.
"Columns beta and se are required for fine-mapping.
Column maf is needed to output posterior effect size estimates on the
allelic scale. All other columns are not required for computations and
can be specified arbitrarily."
}
\examples{
data("locus_dir"); data("BST1"); data("BST1_LD_matrix");
finemap_DT <- BST1
dir.create(file.path(locus_dir,"FINEMAP"), showWarnings = FALSE, recursive = TRUE)
out <- subset_common_snps(LD_matrix=LD_matrix, finemap_DT=finemap_DT)
LD_matrix <- out$LD
finemap_DT <- out$DT
dat_paths <- FINEMAP.construct_data(locus_dir=locus_dir, subset_DT=finemap_DT, LD_matrix=LD_matrix)
}
\seealso{
Other FINEMAP: 
\code{\link{FINEMAP.construct_master}()},
\code{\link{FINEMAP.find_executable}()},
\code{\link{FINEMAP.process_results}()},
\code{\link{FINEMAP}()}
}
\concept{FINEMAP}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloaders.R
\name{downloader}
\alias{downloader}
\title{Downloaders}
\usage{
downloader(
  input_url,
  output_path,
  download_method = "axel",
  background = F,
  force_overwrite = F,
  quiet = F,
  show_progress = T,
  continue = T,
  nThread = 4,
  alternate = T,
  check_certificates = F,
  verbose = T,
  conda_env = NULL
)
}
\description{
R wrapper for wget and axel
}
\seealso{
Other downloaders: 
\code{\link{axel}()},
\code{\link{wget}()}
}
\concept{downloaders}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.query_vcf}
\alias{LD.query_vcf}
\title{Query vcf file}
\usage{
LD.query_vcf(
  subset_DT,
  vcf_URL,
  locus_dir,
  LD_reference,
  whole_vcf = F,
  force_new_vcf = F,
  remove_original_vcf = F,
  download_method = "wget",
  query_by_regions = F,
  nThread = 4,
  conda_env = "echoR",
  verbose = T
)
}
\description{
Query vcf file
}
\examples{
\dontrun{
data("locus_dir"); data("BST1");
# Custom
LD_reference <- "~/Desktop/results/Reference/custom_panel_chr4.vcf"
vcf_file <- LD.index_vcf(vcf_file=LD_reference)
vcf_subset <- LD.query_vcf(subset_DT=BST1, locus_dir=locus_dir, vcf_URL=vcf_file, LD_reference=LD_reference, force_new_vcf=T)
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.process_results}
\alias{PAINTOR.process_results}
\title{Process PAINTOR results}
\usage{
PAINTOR.process_results(PT_results_path, locus_name = "Locus1")
}
\description{
@keywords internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.prepare_LD.transethnic}
\alias{PAINTOR.prepare_LD.transethnic}
\title{Prepare transethnic LD files for PAINTOR}
\usage{
PAINTOR.prepare_LD.transethnic(
  subset_path,
  subset_DT,
  PT_results_path,
  locus,
  GWAS_populations = NULL,
  QTL_datasets,
  QTL_populations = c("AFA", "CAU", "HIS"),
  LD_reference = "1KG_Phase1",
  force_new_LD = F,
  shared_snps_only = T,
  fillNA = 0
)
}
\description{
Prepare transethnic LD files for PAINTOR
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_loci.R
\name{order_loci}
\alias{order_loci}
\title{Order loci by UCS size, or alphabetically}
\usage{
order_loci(dat, merged_dat, by_UCS_size = F, descending = T, verbose = F)
}
\description{
Order loci by UCS size, or alphabetically
}
\examples{
data("merged_DT");
... by UCS size ...
merged_dat <- order_loci(dat=merged_DT, merged_dat=merged_DT, descending=F)
... alphabetically ...
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FINEMAP.R
\name{FINEMAP}
\alias{FINEMAP}
\title{Fine-map locus with \code{FINEMAP}}
\source{
\url{http://www.christianbenner.com}
}
\usage{
FINEMAP(
  subset_DT,
  locus_dir,
  LD_matrix,
  FINEMAP_path = NULL,
  n_samples = NULL,
  n_causal = 5,
  model = "sss",
  remove_tmps = F,
  credset_thresh = 0.95,
  finemap_version = "1.4",
  server = F,
  args_list = list(),
  verbose = T
)
}
\arguments{
\item{FINEMAP_path}{Path to a custom FINEMAP executable to use
instead of the ones included in \pkg{echolocatoR}.
Users can also simply supply "finemap" if this command is linked to the executable.}

\item{n_causal}{The maximum number of potential causal SNPs per locus.
This parameter is used somewhat differntly by different fine-mapping tools.
See tool-specific functions for details.}

\item{model}{"cond" for stepwise conditional search, "sss" for stochastic shotgun search.}

\item{remove_tmps}{Whether to remove any temporary files (e.g. FINEMAP output files) after the pipeline is done running.}

\item{finemap_version}{Which FINEMAP version to use (specify as a string).}

\item{server}{Whether \pkg{echolocatoR} is being run on a computing cluster/server or on a local machine.}

\item{args_list}{A named list of additional arguments to pass to FINEMAP
(e.g.: args_list = list("--n-iterations"=5000,"--sss"="")).
Alternatively, can supply a string instead (e.g.: args_list = "--n-iterations 5000 --sss").}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
The stepwise conditional search starts with a causal configuration containing the
SNP with the lowest P-value alone and then iteratively adds to the causal configuration
the SNP given the highest posterior model probability until no further SNP yields
a higher posterior model probability.
}
\examples{
data("locus_dir"); data("BST1"); data("BST1_LD_matrix");
finemap_DT <- BST1
locus_dir <- here::here(locus_dir)
dir.create(file.path(locus_dir,"FINEMAP"), showWarnings = FALSE, recursive = TRUE)
out <- subset_common_snps(BST1_LD_matrix, finemap_DT)
LD_matrix <- out$LD
subset_DT <- out$DT
subset_DT $N<- subset_DT$N_cases+subset_DT$N_controls
finemap_DT <- FINEMAP(subset_DT=subset_DT, locus_dir=locus_dir, LD_matrix=LD_matrix, finemap_version="1.3")
}
\seealso{
Other FINEMAP: 
\code{\link{FINEMAP.construct_data}()},
\code{\link{FINEMAP.construct_master}()},
\code{\link{FINEMAP.find_executable}()},
\code{\link{FINEMAP.process_results}()}
}
\concept{FINEMAP}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{VALIDATION.bootstrap}
\alias{VALIDATION.bootstrap}
\title{Perform validation bootstrap procedure}
\source{
\href{https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/}{Wilcox test tutorial}
}
\usage{
VALIDATION.bootstrap(
  metric_df,
  metric,
  predictor = "SNP_group",
  validation_method = NULL,
  synthesize_random = F,
  snp_groups = c("Random", "GWAS lead", "UCS (-PolyFun)", "UCS",
    "Consensus (-PolyFun)", "Consensus"),
  test_method = "coin_wilcox_test",
  locus_means = T,
  iterations = 1000,
  save_path = F,
  nThread = 4,
  verbose = T
)
}
\description{
Perform validation bootstrap procedure
}
\examples{
\dontrun{
save_path <- root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"

#### h2 ####
## mean
path <- file.path(root,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.csv.gz")
metric <- "h2.enrichment"
## raw
path <- file.path(root,"PolyFun/h2_merged.csv.gz")
metric <- "SNPVAR"

#### IMPACT ####
## mean
path <- file.path(root,"IMPACT/TOP_IMPACT_all.csv.gz");
metric <- "mean_IMPACT"
## raw
path <- file.path(root,"IMPACT/IMPACT_overlap.csv.gz")
metric <- "IMPACT_score"


## Import and Process ##
metric_df <- data.table::fread(path, nThread=8)
if(metric=="IMPACT_score") metric_df <- subset(metric_df, select=c(SNP,leadSNP,ABF.Credible_Set,ABF.PP,SUSIE.Credible_Set,SUSIE.PP,POLYFUN_SUSIE.Credible_Set,POLYFUN_SUSIE.PP,FINEMAP.Credible_Set,FINEMAP.PP,Consensus_SNP,Support,Locus,IMPACT_score))
if(metric=="mean_IMPACT") metric_df <- find_consensus_SNPs_no_PolyFun(metric_df)

#### run bootstrap ####
boot_res <- VALIDATION.bootstrap(metric_df=metric_df, metric=metric, nThread=8, save_path=gsub("\\\\.csv\\\\.gz",".bootstrap.coin_wilcox_test.csv.gz",path) )
}
}
\seealso{
Other VALIDATION: 
\code{\link{VALIDATION.bootstrap_multimetric}()},
\code{\link{VALIDATION.bootstrap_plot}()},
\code{\link{VALIDATION.compare_bootstrap_distributions}()},
\code{\link{VALIDATION.permute}()},
\code{\link{VALIDATION.super_plot}()}
}
\concept{VALIDATION}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Roadmap.R
\name{ROADMAP.tabix}
\alias{ROADMAP.tabix}
\title{Query Roadmap API}
\usage{
ROADMAP.tabix(
  results_path,
  chrom,
  min_pos,
  max_pos,
  eid,
  convert_to_GRanges = T,
  conda_env = "echoR"
)
}
\arguments{
\item{results_path}{Where to store query results.}

\item{chrom}{Chromosome to query}

\item{min_pos}{Minimum genomic position}

\item{max_pos}{Maximum genomic position}

\item{eid}{Roadmap annotation ID}

\item{convert_to_GRanges}{Whether to return query
as a \code{data.frame} or \code{\link[GenomicRanges]{GRanges}}.}
}
\description{
Query Roadmap epigenomic annotations (chromatin marks)
using a range of genomic coordinates.
}
\seealso{
Other ROADMAP: 
\code{\link{ROADMAP.construct_reference}()},
\code{\link{ROADMAP.merge_and_process_grl}()},
\code{\link{ROADMAP.query_and_plot}()},
\code{\link{ROADMAP.query}()}
}
\concept{ROADMAP}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SUSIE.R
\name{SUSIE}
\alias{SUSIE}
\title{Fine-map with SUSIE}
\source{
\href{https://stephenslab.github.io/susieR/}{GitHub}
\href{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12388}{Publication}
}
\usage{
SUSIE(
  subset_DT,
  LD_matrix,
  dataset_type = "GWAS",
  max_causal = 5,
  sample_size = NULL,
  prior_weights = NULL,
  PP_threshold = 0.95,
  scaled_prior_variance = 0.001,
  estimate_residual_variance = F,
  estimate_prior_variance = T,
  residual_variance = NULL,
  max_iter = 100,
  estimate_prior_method = "optim",
  manual_var_y = F,
  rescale_priors = T,
  plot_track_fit = F,
  return_all_CS = T,
  verbose = T
)
}
\arguments{
\item{dataset_type}{The kind dataset you're fine-mapping (e.g. GWAS, eQTL, tQTL).
This will also be used when creating the subdirectory where your results will be stored
(e.g. \emph{Data/<dataset_type>/Kunkle_2019}).}

\item{max_causal}{The maximum number of non-zero effects (and thus causal variants).}

\item{sample_size}{The overall sample size of the study.
If none is given, and \strong{N_cases} and \strong{N_controls} columns are present,
then sample_size is inferred to be:  \code{max(N_cases) + max(N_controls)}.}

\item{prior_weights}{A vector of length p, in which each entry
gives the prior probability that corresponding column of X has a
nonzero effect on the outcome, y.}

\item{PP_threshold}{The minimum fine-mapped posterior probability for a SNP to be considered part of a Credible Set.
For example, \code{PP_threshold=.95} means that all Credible Set SNPs will be 95\% Credible Set SNPs.}

\item{scaled_prior_variance}{The prior variance, divided by
\code{var(y)} (or by \code{(1/(n-1))yty} for
\code{susie_suff_stat}); that is, the prior variance of each
non-zero element of b is \code{var(y) * scaled_prior_variance}. The
value provided should be either a scalar or a vector of length
\code{L}. If \code{estimate_prior_variance = TRUE}, this provides
initial estimates of the prior variances.}

\item{estimate_residual_variance}{If
\code{estimate_residual_variance = TRUE}, the residual variance is
estimated, using \code{residual_variance} as an initial value. If
\code{estimate_residual_variance = FALSE}, the residual variance is
fixed to the value supplied by \code{residual_variance}.}

\item{estimate_prior_variance}{If \code{estimate_prior_variance =
TRUE}, the prior variance is estimated (this is a separate
parameter for each of the L effects). If provided,
\code{scaled_prior_variance} is then used as an initial value for
the optimization. When \code{estimate_prior_variance = FALSE}, the
prior variance for each of the L effects is determined by the
value supplied to \code{scaled_prior_variance}.}

\item{residual_variance}{Variance of the residual. If
\code{estimate_residual_variance = TRUE}, this value provides the
initial estimate of the residual variance. By default, it is set to
\code{var(y)} in \code{susie} and \code{(1/(n-1))yty} in
\code{susie_suff_stat}.}

\item{max_iter}{Maximum number of IBSS iterations to perform.}

\item{estimate_prior_method}{The method used for estimating prior
variance. When \code{estimate_prior_method = "simple"} is used, the
likelihood at the specified prior variance is compared to the
likelihood at a variance of zero, and the setting with the larger
likelihood is retained.}

\item{rescale_priors}{If prior probabiltities are supplied,
rescale them from 0-1 (i.e. \code{rescaled_priors = priors / sum(priors)}).}

\item{plot_track_fit}{Record each iteration and make a GIF of the
fine-mapping algorithm learning the causal variants.
\strong{WARNING!:} Making this plot can take a long time if there's many iterations.}

\item{verbose}{If \code{verbose = TRUE}, the algorithm's progress,
and a summary of the optimization settings, are printed to the
console.}
}
\description{
Sum of Single Effects (SuSiE): Iterative Bayesian Step-wise Selection
}
\details{
\strong{Notes on convergence:}
\pkg{susieR} will often give the warning: \code{IBSS algorithm did not converge in 100 iterations!}.
This means the results might not necessarily be reliable.
There's several things you can try to avoid this:
\itemize{
\item{Make sure \code{susieR} is up-to-date: \code{ devtools::install_github("stephenslab/susieR@0.9.0")}}
\item{Increase \code{max_causal} (e.g. 5 => 10).}
\item{Increase \code{max_iter} (e.g. 100 => 1000), though this will take longer.}
\item{Decrease the locus window size, which will also speed up the algorithm but potentially miss causal variants far from the lead SNP.}
}
Changing \code{estimate_prior_method} does not seem to affect covergence warnings.

\strong{Notes on variance:}
\href{https://github.com/stephenslab/susieR/issues/90}{GitHub Issue}
If \code{estimate_residual_variance=TRUE} \emph{without} providing \code{var_y}
\emph{and} \code{L>1}, \pkg{susieR} will throw error:
\code{Estimating residual variance failed: the estimated value is negative}
Running \pkg{susieR} with \code{var_y = var(b)} provides \emph{exactly} the same results.
}
\examples{
data("BST1"); data("LD_matrix");
# LD_matrix <- readRDS("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/BST1/plink/UKB_LD.RDS")
finemap_DT <- SUSIE(subset_DT=BST1, LD_matrix=LD_matrix, estimate_residual_variance=T)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabix.R
\name{TABIX.convert_file}
\alias{TABIX.convert_file}
\title{Convert summary stats file to tabix format}
\usage{
TABIX.convert_file(
  fullSS_path,
  chrom_col = "CHR",
  position_col = "POS",
  verbose = TRUE
)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Convert summary stats file to tabix format
}
\examples{
\dontrun{
fullSS_path <- example_fullSS()
fullSS_tabix <- TABIX.convert_file(fullSS_path=fullSS_path, chrom_col="CHR", position_col="POS")
}
}
\seealso{
Other query functions: 
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{IMPACT_annotation_key}
\alias{IMPACT_annotation_key}
\title{IMPACT annotation key}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 1511 rows and 7 columns.
}
\source{
\url{https://github.com/immunogenomics/IMPACT}
}
\usage{
IMPACT_annotation_key
}
\description{

}
\examples{
\dontrun{
IMPACT_annotation_key <- data.table::fread("~/Desktop/Fine_Mapping/echolocatoR/annotations/IMPACT/IMPACT_annotation_key.txt.gz")
usethis::use_data(IMPACT_annotation_key, overwrite = T)
}
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_heatmap}()},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NOTT_2019.bigwig_metadata}
\alias{NOTT_2019.bigwig_metadata}
\title{Metadata and links to data}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 18 rows and 14 columns.
}
\source{
\url{https://science.sciencemag.org/content/366/6469/1134}
}
\usage{
NOTT_2019.bigwig_metadata
}
\description{
Metadata for cell type-specific epigenomic bigWig files hosted on UCSC Genome Browser.
bigWig files contain the genomic ranges from each epigenomic assay,
as well as a Score column which describes the peaks of the aggregate reads.
}
\examples{
\dontrun{
NOTT_2019.bigwig_metadata <- data.table::data.table(readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Nott_2019/Nott_2019.snEpigenomics.xlsx"))
usethis::use_data(NOTT_2019.bigwig_metadata, overwrite = T)
}
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{gather.transethnic.LD}
\alias{gather.transethnic.LD}
\title{Prepare transethnic PAINTOR results}
\usage{
gather.transethnic.LD(
  merged_dat,
  locus,
  conditions = c("Nalls23andMe_2019", "MESA_AFA", "MESA_CAU", "MESA_HIS")
)
}
\description{
Prepare transethnic PAINTOR results
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.snp_group_boxplot}
\alias{IMPACT.snp_group_boxplot}
\title{\emph{IMPACT} box plot with \pkg{ggpubr}}
\source{
\href{https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/}{ggpubr example}
}
\usage{
IMPACT.snp_group_boxplot(
  TOP_IMPACT_all,
  snp_groups = c("GWAS lead", "UCS", "Consensus"),
  method = "wilcox.test",
  comparisons_filter = function(x) {     if ("Consensus" \%in\% x)          return(x) },
  show_plot = T,
  save_path = F,
  title = "IMPACT scores",
  xlabel = NULL,
  ylabel = NULL,
  show_padj = T,
  show_signif = T,
  vjust_signif = 0.5,
  show_xtext = T,
  shift_points = T,
  height = 10,
  width = 10
)
}
\description{
Box plot of \emph{IMPACT} scores from each SNP group.
}
\examples{
\dontrun{
TOP_IMPACT_all <- reshape2::melt(boxplot_mat) \%>\% `colnames<-`(c("SNP_group","max_IMPACT"))
bp <- IMPACT.snp_group_boxplot(TOP_IMPACT_all, method="t.test")
bp <- IMPACT.snp_group_boxplot(TOP_IMPACT_all, method="wilcox.test")
}
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_annotations}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{IMPACT_heatmap}()},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{detect_genes}
\alias{detect_genes}
\title{Detect QTL genes in full summary stats file}
\usage{
detect_genes(loci, verbose = T)
}
\description{
Allows summary stats from different genes to be
fine-mapped separately.
}
\examples{
loci <- c("BST1","LRKR2","MEX3C")
detect_genes(loci)
loci <- c(BST1="BST1", LRRK2="LRRK2", MEX3C="MEX3C")
detect_genes(loci)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BST1}
\alias{BST1}
\title{\pkg{echolocatoR} output example: BST1 locus}
\format{
data.table
\describe{
  \item{SNP}{SNP RSID}
  \item{CHR}{Chromosome}
  \item{POS}{Genomic position (in basepairs)}
  \item{...}{Optional: extra columns}
}
}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
BST1
}
\description{
An example results file after running
\code{finemap_loci} on the \emph{BST1} locus.
}
\details{
Data originally comes from the Parkinson's disease GWAS
by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al., (bioRxiv)}.
}
\examples{
\dontrun{
root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/BST1/Multi-finemap"
BST1 <- data.table::fread(file.path(root_dir,"Multi-finemap_results.txt"))
BST1 <- update_cols(finemap_dat=BST1)
BST1 <- find_consensus_SNPs(finemap_dat=BST1)
usethis::use_data(BST1, overwrite = T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xLiftOver.R
\name{xLiftOver}
\alias{xLiftOver}
\title{Genome build liftover}
\usage{
xLiftOver(
  data.file,
  format.file = c("data.frame", "bed", "chr:start-end", "GRanges"),
  build.conversion = c(NA, "hg38.to.hg19", "hg19.to.hg38", "hg19.to.hg18",
    "hg18.to.hg38", "hg18.to.hg19"),
  merged = T,
  verbose = T,
  RData.location = "http://galahad.well.ox.ac.uk/bigdata",
  guid = NULL
)
}
\description{
Transfer your genomic coordinates from one genome build to another.
}
\details{
\code{xLiftOver} was extracted from the \href{http://xgr.r-forge.r-project.org}{XGR} package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.create_echoR_env}
\alias{CONDA.create_echoR_env}
\title{Create conda env for \emph{echolocatoR}}
\usage{
CONDA.create_echoR_env(
  conda_env = "echoR",
  python_version = NULL,
  channels = c("conda-forge", "bioconda", "r"),
  python_packages = c("pandas>=0.25.0", "pyarrow", "scikit-learn", "bitarray",
    "networkx", "rpy2", "scipy", "pandas-plink"),
  r_packages = c("r-base", "r-devtools"),
  cli_packages = c("tabix", "plink", "macs2"),
  force_install = F,
  auth_token = devtools::github_pat()
)
}
\arguments{
\item{envname}{The conda environment where you want to install \emph{plink}.
By default uses \emph{echoR}, the conda environment distributed with \emph{echolocatoR}.}
}
\description{
Create a new env (or update and existing one)
with the necessary Python, R, and command line packages
to run \emph{echolocatoR}.
}
\details{
\describe{
\item{plink}{https://anaconda.org/bioconda/plink}
\item{tabix}{https://anaconda.org/bioconda/tabix}
}
}
\seealso{
Other conda: 
\code{\link{CONDA.activate_env}()},
\code{\link{CONDA.env_from_yaml}()},
\code{\link{CONDA.find_env_Rlib}()},
\code{\link{CONDA.find_python_path}()},
\code{\link{CONDA.install}()}
}
\concept{conda}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.find_env_Rlib}
\alias{CONDA.find_env_Rlib}
\title{Find the R library for a specific env}
\usage{
CONDA.find_env_Rlib(conda_env = "echoR")
}
\description{
Find the R library for a specific env
}
\seealso{
Other conda: 
\code{\link{CONDA.activate_env}()},
\code{\link{CONDA.create_echoR_env}()},
\code{\link{CONDA.env_from_yaml}()},
\code{\link{CONDA.find_python_path}()},
\code{\link{CONDA.install}()}
}
\concept{conda}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.load_or_create}
\alias{LD.load_or_create}
\title{Procure an LD matrix for fine-mapping}
\usage{
LD.load_or_create(
  locus_dir,
  subset_DT,
  force_new_LD = F,
  LD_reference = "1KGphase1",
  LD_genome_build = "hg19",
  superpopulation = "EUR",
  remote_LD = T,
  download_method = "direct",
  vcf_folder = NULL,
  LD_block = F,
  LD_block_size = 0.7,
  remove_correlates = F,
  fillNA = 0,
  verbose = T,
  server = F,
  remove_tmps = T,
  conda_env = "echoR",
  nThread = 4
)
}
\arguments{
\item{subset_DT}{The locus subset of the full summary stats file.}

\item{force_new_LD}{By default, if an LD matrix file for a given locus is already present,
then \pkg{echolocatoR} will just use the preexisting file.
Set \code{force_new_LD=T} to override this and extract a new subset.}

\item{LD_reference}{Which linkage disequilibrium reference panel do you want to use.
Options include:
\describe{
\item{"UKB"}{A pre-caclulated LD reference matrix from a subset of caucasian British individuals from the UK Biobank. See \href{https://www.biorxiv.org/content/10.1101/807792v2}{Wiessbrod et al. (2019)} for more details.}
\item{"1KGphase1"}{Download a subset of the 1000 Genomes Project Phase 1 vcf and calculate LD on the fly with plink.}
\item{"1KGphase3"}{Download a subset of the 1000 Genomes Project Phase 3 vcf and calculate LD on the fly with plink.}
\item{"<path>/*.vcf" or "<path>/*.vcf.gz"}{Alternatively, users can provide their own custom panel by supplying a list of \emph{.vcf} file path (one per locus) which \pkg{echolocatoR} will use to compute LD (using \emph{plink}).}
}}

\item{superpopulation}{Subset your LD reference panel by superopulation.
Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
\href{https://www.internationalgenome.org/faq/which-populations-are-part-your-study/}{1KGphase1 options} include:
\describe{
\item{"AFR"}{African [descent]}
\item{"AMR"}{Ad-mixed American}
\item{"EAS"}{East Asian}
\item{"EUR"}{European}
\item{"SAS"}{South Asian}
}}

\item{remote_LD}{When acquiring LD matrices,
the default is to delete the full vcf or npz files after \pkg{echolocatoR} has extracted the necssary subset.
However, if you wish to keep these full files (which can be quite large) set \code{remote_LD=T}.}

\item{LD_block}{Calculate LD blocks with \emph{plink} and only include the block to which the lead SNP belongs.}

\item{LD_block_size}{Adjust the granularity of block sizes when \code{LD_block=T}.}

\item{remove_correlates}{A named list, where the names are the RSIDs of SNPs
whose LD correlates you wish to remove,
and the value is the absolute r2 threshold you wish to filter at for each RSID respectively
(e.g. \code{ remove_correlates = c("rs76904798"=.2, "rs10000737"=.8)}).
This will also remove the SNPs in \code{remove_correlates} themselves.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}

\item{server}{Whether \pkg{echolocatoR} is being run on a computing cluster/server or on a local machine.}

\item{remove_tmps}{Whether to remove any temporary files (e.g. FINEMAP output files) after the pipeline is done running.}

\item{conda_env}{The name of a conda environment to use.}
}
\value{
A symmetric LD matrix of pairwise \emph{r} values.
}
\description{
Calculate and/or query linkage disequilibrium (LD) from reference panels (UK Biobank, 1000 Genomes),
or user-supplied datasets.
}
\details{
Options:
\itemize{
\item Download pre-computed LD matrix from UK Biobank.
\item Download raw vcf file from 1KG and compute LD on the fly.
\item Compute LD on the fly from a user-supplied vcf file.
\item Use a user-supplied pre-computed LD-matrix.
}
}
\examples{
\dontrun{
data("BST1"); data("locus_dir");
locus_dir <- file.path("~/Desktop",locus_dir)
 # BST1 <- limit_SNPs(500, BST1)

# UK Biobank LD
LD_matrix <- LD.load_or_create(locus_dir=locus_dir, subset_DT=BST1, LD_reference="UKB")

# 1000 Genomes
LD_matrix <- LD.load_or_create(locus_dir=locus_dir, subset_DT=BST1, LD_reference="1KGphase1", force_new_LD=T)

# Local vcf file
LD_reference="~/Desktop/results/GWAS/Nalls23andMe_2019/BST1/LD/BST1.1KGphase1.vcf.gz"
LD_matrix <- LD.load_or_create(locus_dir=locus_dir, subset_DT=BST1, LD_reference=LD_reference, force_new_LD=T)
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.plot_peaks}
\alias{XGR.plot_peaks}
\title{Plot XGR peaks}
\usage{
XGR.plot_peaks(
  gr.lib,
  subset_DT,
  fill_var = "Assay",
  facet_var = "Source",
  geom = "density",
  locus = NULL,
  adjust = 0.2,
  show_plot = T,
  show.legend = T,
  as.ggplot = T,
  trim_xlims = F
)
}
\arguments{
\item{gr.lib}{\code{GRanges} object of annotations.}

\item{subset_DT}{Data.frame with at least the following columns:
\describe{
\item{SNP}{SNP RSID}
\item{CHR}{chromosome}
\item{POS}{position}
}}

\item{geom}{Plot type ("density", or "histogram").}

\item{locus}{Locus name (\emph{optional}).}

\item{adjust}{The granularity of the peaks.}

\item{show_plot}{Print the plot.}
}
\value{
\code{ggbio} track plot.
}
\description{
Plots the distribution of annotations across a genomic region (x-axis).
}
\examples{
\dontrun{
data("BST1")
finemap_DT <- BST1
gr.lib <- XGR.download_and_standardize(c("ENCODE_DNaseI_ClusteredV3_CellTypes"), finemap_dat=finemap_DT, nCores=1)
gr.filt <- XGR.filter_sources(gr.lib=gr.lib, n_top_sources=5)
gr.filt <- XGR.filter_assays(gr.lib=gr.filt, n_top_assays=5)
xgr.track <- XGR.plot_peaks(gr.lib=gr.filt, subset_DT=finemap_DT, fill_var="Assay", facet_var="Source")
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.prepare_annotations}
\alias{PAINTOR.prepare_annotations}
\title{Prepare PAINTOR annotations}
\usage{
PAINTOR.prepare_annotations(
  paintor_path = NULL,
  BED_paths,
  PT_results_path,
  locus_name,
  remove_BED = F
)
}
\description{
Prepare PAINTOR annotations
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpliceAI.R
\name{SPLICEAI.snp_probs}
\alias{SPLICEAI.snp_probs}
\title{Postprocess \emph{SpliceAI} results after querying}
\source{
\href{https://github.com/Illumina/SpliceAI}{GitHub}
\href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
}
\usage{
SPLICEAI.snp_probs(DAT, save_path = F)
}
\description{
Postprocess \emph{SpliceAI} results after querying
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
DAT <- data.table::fread(file.path(root, "Data/GWAS/Nalls23andMe_2019/_genome_wide/SpliceAI/spliceAI_Nalls23andMe_2019.hits.csv.gz"))
}
}
\seealso{
Other SpliceAI: 
\code{\link{SPLICEAI.plot}()},
\code{\link{SPLICEAI.run}()},
\code{\link{SPLICEAI.subset_precomputed_tsv_iterate}()},
\code{\link{SPLICEAI.subset_precomputed_tsv}()},
\code{\link{SPLICEAI.subset_precomputed_vcf}()}
}
\concept{SpliceAI}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fullSS_dat}
\alias{fullSS_dat}
\title{Example subset of full summary stats}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 38066 rows and 11 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
fullSS_dat
}
\description{
Data originally comes from the Parkinson's disease GWAS
by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
}
\examples{
A subset of the full GWAS summary stats from Nalls et al. (2019)
\dontrun{
data("top_SNPs")
loci <- c("BST1","LRRK2","MEX3C")
fullSS_dat <- lapply(loci, function(locus){
  top_sub <- subset(top_SNPs, Locus==locus)
  TABIX(fullSS_path = "~/Desktop/nallsEtAl2019_allSamples_allVariants.mod.txt.gz",
        save_subset = F,
        subset_path = NULL,
        chrom_col = "CHR", position_col = "POS",
        chrom = top_sub$CHR[1],
        min_POS = top_sub$POS - 500000*4,
        max_POS = top_sub$POS + 500000*4)
}) \%>\% data.table::rbindlist()
usethis::use_data(fullSS_dat, overwrite=T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.gather_ldscores}
\alias{POLYFUN.gather_ldscores}
\title{Find and import LD score files}
\usage{
POLYFUN.gather_ldscores(output_prefix)
}
\description{
Find and import LD score files
}
\examples{
\dontrun{
output_prefix <- file.path(system.file("tools/polyfun/gold","",package = "echolocatoR"),"testrun.22")
ldscore <- POLYFUN.gather_ldscores(output_prefix=output_prefix)
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COJO.R
\name{process_COJO_results}
\alias{process_COJO_results}
\title{Gather all \emph{GCTA-COJO} results}
\usage{
process_COJO_results(subset_DT, locus_dir, freq_cutoff = 0.1)
}
\description{
Gather all \emph{GCTA-COJO} results
}
\seealso{
Other COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}: 
\code{\link{COJO.stepwise}()},
\code{\link{COJO}()},
\code{\link{get_stepwise_results}()}
}
\concept{COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.get_ldscores}
\alias{IMPACT.get_ldscores}
\title{Download LD scores from IMPACT publication}
\source{
T Amariuta et al., IMPACT: Genomic Annotation of Cell-State-Specific Regulatory Elements Inferred from the Epigenome of Bound Transcription Factors. The American Journal of Human Genetics, 1–17 (2019).
}
\usage{
IMPACT.get_ldscores(chrom = NULL, subset_DT = NULL, nThread = 4)
}
\description{
Downloads per-SNP LD scores from LD SCore regression,
first conducted in the original IMPACT publication.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{popDat_1KGphase3}
\alias{popDat_1KGphase3}
\title{Population metadata: 1KGphase3}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 2504 rows and 4 columns.
}
\usage{
popDat_1KGphase3
}
\description{
Population metadata: 1KGphase3
}
\examples{
\dontrun{
popDat_URL <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
popDat_1KGphase3 <-  data.table::fread(text=trimws(gsub(",\t",",",readLines(popDat_URL))), sep="\t",  fill=T, stringsAsFactors = F, col.names = c("sample","population","superpop","sex"), nThread = 4)
usethis::use_data(popDat_1KGphase3, overwrite = T)
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.plink_file}
\alias{LD.plink_file}
\title{Find correct plink file}
\usage{
LD.plink_file(plink = "plink", conda_env = NULL)
}
\description{
Find correct plink file
}
\examples{
plink <- LD.plink_file()
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.superenhancers}
\alias{NOTT_2019.superenhancers}
\title{Get cell type-specific superenhancer data}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.superenhancers()
}
\description{
Brain cell-specific epigenomic data from Nott et al. (2019).
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.find_polyfun_folder}
\alias{POLYFUN.find_polyfun_folder}
\title{Folder where PolyFun submodule is stored}
\usage{
POLYFUN.find_polyfun_folder(polyfun_path = NULL)
}
\description{
Folder where PolyFun submodule is stored
}
\examples{
polyfun_path <- POLYFUN.find_polyfun_folder()
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_coloc_results.R
\name{merge_coloc_results}
\alias{merge_coloc_results}
\title{Convert Jack's coloc results to \emph{echolocatoR} format}
\usage{
merge_coloc_results(
  all_obj,
  results_level = c("summary"),
  nThread = 4,
  verbose = T,
  save_path = F
)
}
\arguments{
\item{all_obj}{Nested list object created by Jack.}

\item{results_level}{Return coloc results at the
"summary" (one row per Locus:eGene pair) or "snp" level (one row per SNP).}

\item{save_path}{File where you want the merged results to be saved.}
}
\description{
Convert Jack's coloc results to \emph{echolocatoR} format
}
\examples{
\dontrun{
load("/sc/hydra/projects/ad-omics/microglia_omics/COLOC/Kunkle_2019/Microglia_all_regions_Kunkle_2019_COLOC.RData")
merged_results <- merge_coloc_results(all_obj=all_obj, results_level="snp", save_path="~/Desktop")
}
}
\seealso{
Other coloc: 
\code{\link{merge_coloc_results_each}()}
}
\concept{coloc}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.read_ld_table}
\alias{LD.read_ld_table}
\title{Create LD matrix from plink output.}
\usage{
LD.read_ld_table(ld.path, snp.subset = F, fillNA = 0, verbose = T)
}
\description{
Depending on which parameters you give \emph{plink} when calculating LD, you get different file outputs.
When it produces an LD table, use this function to create a proper LD matrix.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effective_sample_size.R
\name{effective_sample_size}
\alias{effective_sample_size}
\title{Calculate effective sample size}
\usage{
effective_sample_size(finemap_dat, sample_size = NULL, verbose = T)
}
\arguments{
\item{finemap_dat}{Preprocessed \emph{echolocatoR} locus subset file.
Requires the columns \strong{N_cases} and \strong{N_controls}.}
}
\description{
This is an often-cited approximation for the effective sample size in a case-control study.
i.e., this is the sample size required to identify an association with a
 quantitative trait with the same power as in your present study.
 (from email correpsondence with Omer Weissbrod)
}
\examples{
data("BST1");
finemap_DT <- effective_sample_size(finemap_dat=BST1)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.install}
\alias{CONDA.install}
\title{Install conda if it's missing}
\usage{
CONDA.install(conda_path = "auto", verbose = F)
}
\description{
Install conda if it's missing
}
\seealso{
Other conda: 
\code{\link{CONDA.activate_env}()},
\code{\link{CONDA.create_echoR_env}()},
\code{\link{CONDA.env_from_yaml}()},
\code{\link{CONDA.find_env_Rlib}()},
\code{\link{CONDA.find_python_path}()}
}
\concept{conda}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.create_locusFile}
\alias{PAINTOR.create_locusFile}
\title{Create PAINTOR locus file}
\usage{
PAINTOR.create_locusFile(
  subset_path,
  fullSS,
  PT_results_path,
  GWAS_datasets = "Nalls23andMe_2019",
  QTL_datasets = c("Fairfax_2014_CD14", "Fairfax_2014_IFN", "Fairfax_2014_LPS2",
    "Fairfax_2014_LPS24"),
  locus,
  locus_name
)
}
\description{
Create PAINTOR locus file
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dt.replace}
\alias{dt.replace}
\title{Replace items within DT object}
\usage{
dt.replace(DT, target, replacement)
}
\value{
data.frame
}
\description{
Annoyingly, there is no native function to do simple find-and-replace in the `DT` library.
}
\seealso{
Other general: 
\code{\link{createDT_html}()},
\code{\link{createDT}()},
\code{\link{example_fullSS}()},
\code{\link{fillNA_CS_PP}()},
\code{\link{get_os}()},
\code{\link{get_sample_size}()},
\code{\link{tryFunc}()}
}
\concept{general}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.iterate_overlap}
\alias{XGR.iterate_overlap}
\title{Check overlap with XGR annotations}
\usage{
XGR.iterate_overlap(
  lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes", "TFBS_Conserved",
    "ReMap_PublicAndEncode_TFBS", "Uniform_TFBS"),
  subset_DT,
  save_path = F,
  nCores = 4
)
}
\arguments{
\item{lib.selections}{Which XGR annotations to check overlap with.
For full list of libraries see:
 \url{http://xgr.r-forge.r-project.org/#annotations-at-the-genomic-region-level}}

\item{subset_DT}{Data.frame with at least the following columns:
\describe{
\item{SNP}{SNP RSID}
\item{CHR}{chromosome}
\item{POS}{position}
}}

\item{save_path}{Save the results as a data.frame}

\item{nCores}{Multi-thread across libraries}
}
\description{
Automatically handles different file formats provided by XGR
 (e.g. varying kinds of nested/unnested \code{GRanges}).
Then returns a \code{Granges} object with only the XGR annotation ranges
that overlap with the SNPs in \code{subset_DT}.
The \code{GRanges} merges hits from \code{subset_DT}.
}
\examples{
\dontrun{
data("BST1")
finemap_DT <- BST1
gr.hits <- XGR.iterate_overlap(lib.selections=c("ENCODE_TFBS_ClusteredV3_CellTypes"), subset_DT=finemap_DT, nCores=1)
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.h2_enrichment_SNPgroups_plot}
\alias{POLYFUN.h2_enrichment_SNPgroups_plot}
\title{Plot heritability (h2) enrichment}
\usage{
POLYFUN.h2_enrichment_SNPgroups_plot(
  RES,
  snp_groups = c("GWAS lead", "UCS", "Consensus (-PolyFun)", "Consensus"),
  comparisons_filter = function(x) {     if ("Consensus" \%in\% x)          return(x) },
  method = "wilcox.test",
  remove_outliers = T,
  title = "S-LDSC heritability enrichment",
  xlabel = "SNP group",
  ylabel = bquote(~"h"^2 ~ "enrichment"),
  show_xtext = T,
  show_padj = T,
  show_signif = T,
  vjust_signif = 0.5,
  show_plot = T,
  save_path = F
)
}
\description{
Plot heritability (h2) enrichment
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
merged_dat <- merge_finemapping_results(dataset = file.path("Data/GWAS/Nalls23andMe_2019"), LD_reference = "UKB", minimum_support = 0)
RES <- POLYFUN.h2_enrichment_SNPgroups(merged_dat=merged_dat, out.path="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output")

plot.h2 <- POLYFUN.h2_enrichment_SNPgroups_plot(RES = RES, show_plot = T)
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{hierarchical_colors}
\alias{hierarchical_colors}
\title{Color tissues/cell types}
\usage{
hierarchical_colors(mat_meta)
}
\description{
Use the same color for groups of tissues/cell types but vary the shades within subgroups.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbind_GRanges.R
\name{rbind_GRanges}
\alias{rbind_GRanges}
\title{Bind GRanges with different mcols}
\usage{
rbind_GRanges(gr1, gr2)
}
\description{
Bind GRanges with different mcols
}
\examples{
data("merged_DT")
gr.hits <- CORCES_2020.get_ATAC_peak_overlap(finemap_dat = merged_DT)
gr.hits$extra_col <- "Extra"
gr.anchor_hits <- CORCES_2020.get_HiChIP_FitHiChIP_overlap(finemap_dat = merged_DT)
try({gr.bind <- c(gr.hits, gr.anchor_hits)})
gr.bound <- rbind_GRanges(gr1, gr2)
}
\seealso{
Other utils: 
\code{\link{LIFTOVER}()},
\code{\link{dataframe_2_vcf}()}
}
\concept{utils}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zscore.get_mean_and_sd.R
\name{Zscore.get_mean_and_sd}
\alias{Zscore.get_mean_and_sd}
\title{Store info necessary for Z-score}
\usage{
Zscore.get_mean_and_sd(
  fullSS,
  target_col = "statistic",
  effect_col = "beta",
  stderr_col = "se",
  use_saved = T,
  output_path
)
}
\description{
Store the info about the full vector (length, mean, standard deviation),
so that you can later accurately compute Z-scores
when you're only working with a subset of the data.
}
\details{
These functions are necessary for \code{\link{PAINTOR}}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifbreakR.R
\name{MOTIFBREAKR}
\alias{MOTIFBREAKR}
\title{Run \code{\link{motifbreakR}}}
\source{
\strong{Publication:}
\url{https://pubmed.ncbi.nlm.nih.gov/26272984/}
\strong{GitHub:}
\url{https://github.com/Simon-Coetzee/MotifBreakR}
\strong{Vignette:}
\url{http://simon-coetzee.github.io/motifBreakR}
}
\usage{
MOTIFBREAKR(
  rsid_list,
  save_rds = T,
  dataset_dir = "./results",
  pwmList = NULL,
  organism = NULL,
  threshold = 0.85,
  show.neutral = F,
  method = "default",
  verbose = T,
  calculate_all_pval = T,
  force_new = F
)
}
\value{
motifbreakr results
}
\description{
\code{\link{motifbreakR}} is a package to predict how much a SNP will disrupt
a transcription factor binding motif (if it falls witihn one).
}
\examples{
\dontrun{
# data("merged_DT")
# rsid_list <- unique(subset(merged_DT, Locus \%in\% c("LRRK2","MED12L","DYRK1A","FCGR2A") & (Consensus_SNP | leadSNP))$SNP)
# rsid_list <- unique(subset(merged_DT, Consensus_SNP | leadSNP)$SNP)
rsid_list <- c("rs11175620","rs7294619","rs74324737")
mb.results <- MOTIFBREAKR(rsid_list=rsid_list, calculate_all_pval=T, force_new = T)
}
}
\seealso{
Other motifbreakR: 
\code{\link{MOTIFBREAKR.filter_by_metadata}()}
}
\concept{motifbreakR}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabix.R
\name{TABIX}
\alias{TABIX}
\title{Covert and query}
\usage{
TABIX(
  fullSS_path,
  study_dir = NULL,
  subset_path = NULL,
  is_tabix = FALSE,
  chrom_col = "CHR",
  chrom_type = NULL,
  position_col = "BP",
  min_POS,
  max_POS,
  chrom,
  save_subset = TRUE,
  nThread = 1,
  conda_env = "echoR",
  verbose = TRUE
)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{min_POS}{Manually set the minimum genomic position for your locus subset.
\code{min_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{min_POS=top_SNPs$min_POS}).}

\item{max_POS}{Manually set the maximum genomic position for your locus subset.
\code{max_POS} can clip the window size set by \code{bp_distance}.
Can also be a list of positions (one for each locus) (e.g. \code{max_POS=top_SNPs$max_POS}).}

\item{conda_env}{The name of a conda environment to use.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\value{
data.table of locus subset summary statistics
}
\description{
If it is not tabix format already
(determined by checking for a \link{.tbi} file of the same name in the same directory),
the full summary statistics file is converted into tabix format for super fast querying.
A query is then made using the min/max genomic positions to extract a locus-specific summary stats file.
}
\examples{
\dontrun{
data("Nalls_top_SNPs"); data("top_SNPs")
fullSS_path <- example_fullSS()
top_SNPs_BST1 <- subset(top_SNPs, Locus=='BST1')
bp_distance <- 1e+06
min_POS <- top_SNPs_BST1$POS - bp_distance
max_POS <- top_SNPs_BST1$POS + bp_distance

subset_path <- file.path("BST1_Nalls23andMe_2019_subset.tsv.gz")
dat <- TABIX(fullSS_path=fullSS_path,
             subset_path=subset_path,
             min_POS=min_POS,
             max_POS=max_POS,
             chrom=top_SNPs_BST1$CHR)
}
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.find_package}
\alias{CONDA.find_package}
\title{Find package executable}
\usage{
CONDA.find_package(package, conda_env = "echoR", verbose = T)
}
\description{
Find package executable
}
\examples{
# Tabix
tabix <- CONDA.find_package(package="tabix", conda_env="echoR")
tabix <- CONDA.find_package(package="tabix", conda_env=NULL)
# bgzip
bgzip <- CONDA.find_package(package="bgzip", conda_env="echoR")
bgzip <- CONDA.find_package(package="bgzip", conda_env=NULL)
}
\concept{CONDA}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.plot_top_annotations}
\alias{IMPACT.plot_top_annotations}
\title{IMPACT locus plots of top annotations}
\usage{
IMPACT.plot_top_annotations()
}
\description{
IMPACT locus plots of top annotations
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR}
\alias{PAINTOR}
\title{Run full PAINTOR pipeline}
\source{
\href{https://github.com/gkichaev/PAINTOR_V3.0}{GitHub}
\href{https://github.com/gkichaev/PAINTOR_V3.0/wiki/2a.-Computing-1000-genomes-LD}{LD Tutorial}
}
\usage{
PAINTOR(
  finemap_dat = NULL,
  paintor_path = NULL,
  GWAS_datasets = NULL,
  QTL_datasets = NULL,
  locus,
  locus_name = NULL,
  n_causal = 5,
  XGR_dataset = NA,
  ROADMAP_search = NA,
  chromatin_state = "TssA",
  use_annotations = F,
  PP_threshold = 0.95,
  consensus_thresh = 2,
  multi_finemap_col_name = "PAINTOR",
  trim_gene_limits = F,
  force_new_subset = F,
  GWAS_populations = "EUR",
  QTL_populations = "EUR",
  LD_reference = "1KG_Phase1",
  force_new_LD = F,
  LD_matrix = NULL,
  force_reinstall = F
)
}
\description{
Run full PAINTOR pipeline
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SURE.R
\name{SURE.merge}
\alias{SURE.merge}
\title{Merge with fine-mapping and \emph{SuRE} results}
\source{
\href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
\href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
\href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
\href{https://sure.nki.nl}{SNP-SuRE data browse}
}
\usage{
SURE.merge(merged_dat, sure, save_path = F, verbose = T)
}
\description{
Merge with fine-mapping and \emph{SuRE} results
}
\examples{
sure <- data.table::fread("/sc/arion/projects/pd-omics/data/MPRA/SURE/SuRE_SNP_table_LP190708.txt.gz", nThread=4)
merged_dat <- merge_finemapping_results("Data/GWAS/Nalls23andMe_2019",LD_reference = "UKB", minimum_support = 0)

sure_DT <- SURE.merge(merged_dat=merged_dat, sure=sure, save_path="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/SURE/Nalls23andMe_2019.SURE.csv.gz")
}
\seealso{
Other SURE: 
\code{\link{SURE.melt_snp_groups}()},
\code{\link{SURE.plot}()}
}
\concept{SURE}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Nott_2019.R
\name{NOTT_2019.get_interactions}
\alias{NOTT_2019.get_interactions}
\title{Import cell type-specific interactomes}
\source{
\href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
\url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
}
\usage{
NOTT_2019.get_interactions(finemap_dat, as.granges = F)
}
\description{
Brain cell-specific epigenomic data from Nott et al. (2019).
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancer_interactome}},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_locus_dir}
\alias{make_locus_dir}
\title{Make locus-specific results folder}
\usage{
make_locus_dir(
  results_dir = "./results",
  dataset_type = "dataset_type",
  dataset_name = "dataset_name",
  locus
)
}
\description{
Make locus-specific results folder
}
\examples{
locus_dir <- make_locus_dir(results_dir="./results", dataset_type="GWAS", dataset_name="Nalls23andMe_2019", locus="BST1")
}
\seealso{
Other directory functions: 
\code{\link{get_locus_dir}()},
\code{\link{get_multifinemap_path}()},
\code{\link{get_study_dir}()},
\code{\link{get_subset_path}()}
}
\concept{directory functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.plot_LD}
\alias{LD.plot_LD}
\title{Plot a subset of the LD matrix}
\usage{
LD.plot_LD(
  LD_matrix,
  subset_DT,
  span = 10,
  method = c("gaston", "heatmap", "image")
)
}
\arguments{
\item{span}{This is very computationally intensive,
so you need to limit the number of SNPs with span.
If \code{span=10}, only 10 SNPs upstream and 10 SNPs downstream of the lead SNP will be plotted.}
}
\description{
Uses \code{gaston} to plot a SNP-annotated LD matrix.
}
\examples{
\dontrun{
data("BST1");
LD_matrix <- readRDS("/Volumes/Steelix/fine_mapping_files/GWAS/Nalls23andMe_2019/BST1/plink/UKB_LD.RDS")
LD.plot_LD(LD_matrix=LD_matrix, subset_DT=BST1)
}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_coloc_results_each.R
\name{merge_coloc_results_each}
\alias{merge_coloc_results_each}
\title{Iterate merging of coloc results across files}
\usage{
merge_coloc_results_each(
  coloc_rda_files,
  save_path = F,
  save_each = F,
  results_level = "summary",
  ppH4_thresh = 0,
  return_filter = F,
  force_new = T,
  nThread = 4,
  verbose = T
)
}
\arguments{
\item{coloc_rda_files}{A list of full paths to coloc .RData/.rda results files generated by Jack's coloc scripts.}

\item{save_path}{File path where the merged results will be saved.}

\item{results_level}{Return coloc results at the
"summary" (one row per Locus:eGene pair) or "snp" level (one row per SNP).}

\item{ppH4_thresh}{Threshold to filter summary-level results.}

\item{force_new}{If the \code{save_path} file already exists, overwrite it.}

\item{filter}{Flexible row filtering.}
}
\description{
Iterate merging of coloc results across files
}
\examples{
# List RDS files you want to extract info from
coloc_rds <- list.files("/sc/hydra/projects/ad-omics/microglia_omics/COLOC", pattern = "*_COLOC.RData", recursive = T, full.names = T)
coloc_rds <- coloc_rds[!grepl("*_sQTL_*|Microglia_all_regions_Young",coloc_rds)]
coloc_rds <- coloc_rds[grep("*Microglia_all_regions_*",coloc_rds)]

# Summary-level
coloc_summary <- merge_coloc_results_each(coloc_rds_files=coloc_rds, save_path="/sc/arion/projects/pd-omics/brian/all_COLOC_results.Microglia_all_regions.summary-level.tsv.gz", results_level="summary")
# SNP-level
coloc_snp <- merge_coloc_results_each(coloc_rds_files=coloc_rds, save_path="/sc/arion/projects/pd-omics/brian/all_COLOC_results.Microglia_all_regions.snp-level.tsv.gz", results_level="snp", filter="leadSNP==T")
}
\seealso{
Other coloc: 
\code{\link{merge_coloc_results}()}
}
\concept{coloc}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DT_to_GRanges.R
\name{DT_to_GRanges}
\alias{DT_to_GRanges}
\title{Convert data.table to GRanges object}
\usage{
DT_to_GRanges(
  subset_DT,
  style = "NCBI",
  chrom_col = "CHR",
  position_col = "POS"
)
}
\description{
Convert data.table to GRanges object
}
\seealso{
Other XGR: 
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_study_dir}
\alias{get_study_dir}
\title{Extract the study dir}
\usage{
get_study_dir(locus_dir)
}
\description{
Extract the study dir
}
\seealso{
Other directory functions: 
\code{\link{get_locus_dir}()},
\code{\link{get_multifinemap_path}()},
\code{\link{get_subset_path}()},
\code{\link{make_locus_dir}()}
}
\concept{directory functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORCES_2020.R
\name{CORCES_2020.get_HiChIP_FitHiChIP_overlap}
\alias{CORCES_2020.get_HiChIP_FitHiChIP_overlap}
\title{Get overlap between data table of SNPs and HiChIP_FitHiChIP coaccessibility anchors}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.get_HiChIP_FitHiChIP_overlap(finemap_dat, verbose = T)
}
\description{
Anchors are the genomic regions that have evidence of being
functionally connected to one another (coaccessible),
 e.g. enhancer-promoter interactions.
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda.R
\name{CONDA.find_python_path}
\alias{CONDA.find_python_path}
\title{Find the python file for a specific env}
\usage{
CONDA.find_python_path(conda_env = "echoR", verbose = T)
}
\description{
Find the python file for a specific env
}
\seealso{
Other conda: 
\code{\link{CONDA.activate_env}()},
\code{\link{CONDA.create_echoR_env}()},
\code{\link{CONDA.env_from_yaml}()},
\code{\link{CONDA.find_env_Rlib}()},
\code{\link{CONDA.install}()}
}
\concept{conda}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.locus.R
\name{PLOT.dot_summary}
\alias{PLOT.dot_summary}
\title{Multi-fine-map summary dot plot}
\usage{
PLOT.dot_summary(finemap_dat, PP_threshold = 0.95, show_plot = T)
}
\description{
Multi-fine-map summary dot plot
}
\examples{
data("BST1")
gp <- PLOT.dot_summary(finemap_dat=BST1)
}
\seealso{
Other plot: 
\code{\link{PLOT.add_multitrack_lines}()},
\code{\link{PLOT.get_window_limits}()},
\code{\link{PLOT.get_window_suffix}()}
}
\concept{plot}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.prepare_foreground_background}
\alias{XGR.prepare_foreground_background}
\title{Prepare SNP sets for enrichment}
\usage{
XGR.prepare_foreground_background(
  subset_DT,
  foreground_filter = "Support>0",
  background_filter = NULL,
  fg_sample_size = NULL,
  bg_sample_size = NULL,
  verbose = T
)
}
\arguments{
\item{subset_DT}{Data.frame with at least the following columns:
\describe{
\item{SNP}{SNP RSID}
\item{CHR}{chromosome}
\item{POS}{position}
}}

\item{foreground_filter}{Specify foreground by filtering SNPs in \code{subset_DT}.
Write filter as a string (or \code{NULL} to include all SNPs).}

\item{background_filter}{Specify background by filtering SNPs in \code{subset_DT}.
Write filter as a string (or \code{NULL} to include all SNPs).}
}
\description{
Prepare custom foreground and background SNPs sets for enrichment tests with XGR annotations.
}
\examples{
\dontrun{
data("merged_DT")
merged_dat <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/merged_UKB.csv.gz", nThread=4)
fg_bg <- XGR.prepare_foreground_background(subset_DT=merged_dat, foreground_filter="Consensus_SNP==T", background_filter="leadSNP==T")
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.install}
\alias{PAINTOR.install}
\title{Install PAINTOR via command line}
\usage{
PAINTOR.install(paintor_path = NULL, force_reinstall = F)
}
\description{
Currently there is no R or conda distribution of PAINTOR.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dey_DeepLearning.R
\name{DEEPLEARNING.melt}
\alias{DEEPLEARNING.melt}
\title{Melt deep learning annotations into long-format}
\usage{
DEEPLEARNING.melt(
  ANNOT,
  model = c("Basenji", "BiClassCNN", "DeepSEA", "ChromHMM", "Roadmap", "Others"),
  aggregate_func = "mean",
  replace_NA = NA,
  replace_negInf = NA,
  save_path = F
)
}
\description{
Melt deep learning annotations into long-format
}
\examples{
\dontrun{
 root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
## merged_dat <- merge_finemapping_results(dataset = "Data/GWAS/Nalls23andMe_2019", minimum_support = 0, LD_reference = "UKB")
## ANNOT <- DEEPLEARNING.query(merged_dat=merged_dat, level="Allelic_Effect", type="annot")

#### Allelic_Effect ####
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.csv.gz")

#### Variant_Level ####
path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.csv.gz")

ANNOT <- data.table::fread(path, nThread=8)
ANNOT <- find_consensus_SNPs_no_PolyFun(ANNOT)
annot.melt <- DEEPLEARNING.melt(ANNOT=ANNOT, aggregate_func="mean", save_path=gsub("\\\\.csv\\\\.gz",".snp_groups_mean.csv.gz",path))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.snpstats_get_MAF}
\alias{LD.snpstats_get_MAF}
\title{Get MAF using \pkg{snpStats} package}
\source{
\href{https://www.bioconductor.org/packages/release/bioc/html/snpStats.html}{snpStats Bioconductor page}
}
\usage{
LD.snpstats_get_MAF(
  subset_DT,
  LD_folder,
  plink_prefix = "plink",
  force_new_MAF = F,
  nThread = 4,
  verbose = T
)
}
\arguments{
\item{LD_folder}{Locus-specific LD output folder.}
}
\description{
Get MAF using \pkg{snpStats} package
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MungeSumstats.R
\name{column_map}
\alias{column_map}
\title{Column name mappings for a given package}
\source{
https://github.com/neurogenomics/MungeSumstats
}
\usage{
column_map(package = "MungeSumstats")
}
\description{
Supports \pkg{MungeSumstats} and \pkg{echolocatoR}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.filter_LD}
\alias{LD.filter_LD}
\title{Filter LD}
\usage{
LD.filter_LD(LD_list, remove_correlates = F, min_r2 = 0, verbose = F)
}
\description{
Filter LD
}
\examples{
data("BST1"); data("LD_matrix");
LD_list <- list(LD=LD_matrix, DT=BST1)
LD_list <- LD.filter_LD(LD_list, min_r2=.2)
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.R
\name{calculate_tstat}
\alias{calculate_tstat}
\title{Compute t-stat}
\usage{
calculate_tstat(finemap_dat, tstat_col = "t_stat")
}
\description{
If \strong{tstat} column is missing,
compute t-statistic from: \code{Effect / StdErr}.
}
\seealso{
Other standardization functions: 
\code{\link{auto_topSNPs_sub}()}
}
\concept{standardization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{ANNOTATE.plot_missense}
\alias{ANNOTATE.plot_missense}
\title{Plot any missense variants}
\usage{
ANNOTATE.plot_missense(
  merged_dat,
  snp_filter = "Support>0",
  label_yaxis = F,
  x_label = "UCS missense\\nmutations",
  show.legend = T,
  show_numbers = F,
  show_plot = T
)
}
\description{
Plot any missense variants
}
\examples{
\dontrun{
data("merged_DT");
gg_missense <- ANNOTATE.plot_missense(merged_dat=merged_DT, snp_filter="Support>0")
gg_missense <- ANNOTATE.plot_missense(merged_dat=merged_DT, snp_filter="Consensus_SNP==T")
}
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{genome_wide_dir}
\alias{genome_wide_dir}
\title{Example results path for genome-wide results}
\format{
An object of class \code{character} of length 1.
}
\usage{
genome_wide_dir
}
\description{
Example results path for genome-wide results
}
\examples{
\dontrun{
genome_wide_dir <- "results/GWAS/Nalls23andMe_2019/_genome_wide"
usethis::use_data(genome_wide_dir, overwrite=T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{query_by_gene}
\alias{query_by_gene}
\title{Use \emph{awk} to query locus subsets.}
\usage{
query_by_gene(fullSS_path, subset_path, gene, gene_col, file_sep)
}
\arguments{
\item{fullSS_path}{Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
It is usually best to provide the absolute path rather than the relative path.}

\item{subset_path}{Path of the resulting locus subset file.}

\item{gene}{Gene symbol (e.g. FOXP2) of the gene you want to query.}

\item{gene_col}{Name of the column that contains the gene name.}

\item{file_sep}{The separator in the full summary stats file.
This parameter is only necessary if \code{query_by!="tabix"}.}
}
\description{
To be used when gene name is a column in the full summary stats file.
More commonly useful for QTL full summary stats files.
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{import_topSNPs}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dey_DeepLearning.R
\name{DEEPLEARNING.query_multi_chr}
\alias{DEEPLEARNING.query_multi_chr}
\title{Query deep learning annotations and LDscores (iterate)}
\usage{
DEEPLEARNING.query_multi_chr(
  merged_dat,
  base_url = "/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
  level = c("Variant_Level", "Allelic_Effect"),
  tissue = c("NTS", "Blood", "Brain"),
  model = c("Basenji", "BiClassCNN", "DeepSEA", "ChromHMM", "Roadmap", "Others"),
  mean_max = c("MEAN", "MAX"),
  type = c("annot", "ldscore"),
  nThread = 4,
  verbose = T
)
}
\description{
Query deep learning annotations and LDscores, and
then merge with \code{subset_DT} by \emph{SNP}.
 Repeat for each locus,
}
\examples{
\dontrun{
data("merged_DT")
ANNOT.DAT <- DEEPLEARNING.query_multi_chr(merged_dat=merged_DT, tissue="NTS", model="Basenji", type="annot")
}
}
\seealso{
Other DEEPLEARNING: 
\code{\link{DEEPLEARNING.query_one_chr}()},
\code{\link{DEEPLEARNING.query}()}
}
\concept{DEEPLEARNING}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fillNA_CS_PP.R
\name{fillNA_CS_PP}
\alias{fillNA_CS_PP}
\title{Fill NA in PP and CS columns}
\usage{
fillNA_CS_PP(finemap_dat, fillNA_CS = 0, fillNA_PP = 0)
}
\description{
Fill NA in PP and CS columns
}
\examples{
data("BST1");
finemap_dat <- BST1
# finemap_dat <- data.table::fread("~/Desktop/results/GWAS/Kunkle_2019.microgliaQTL/ABCA7/Multi-finemap/ABCA7_Kunkle_2019.microgliaQTL_Multi-finemap.tsv.gz")
finemap_dat <- fillNA_CS_PP(finemap_dat=finemap_dat)
}
\seealso{
Other general: 
\code{\link{createDT_html}()},
\code{\link{createDT}()},
\code{\link{dt.replace}()},
\code{\link{example_fullSS}()},
\code{\link{get_os}()},
\code{\link{get_sample_size}()},
\code{\link{tryFunc}()}
}
\concept{general}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limit_SNPs.R
\name{limit_SNPs}
\alias{limit_SNPs}
\title{Limit the number of SNPs per locus.}
\usage{
limit_SNPs(max_snps = 500, subset_DT)
}
\arguments{
\item{max_snps}{The maximum number of SNPs to keep in the resulting data.frame.}

\item{subset_DT}{A data.frame that contains at least the following columns:
\describe{
  \item{SNP}{RSID for each SNP.}
  \item{POS}{Each SNP's genomic position (in basepairs).}
}}
}
\description{
Start with the lead SNP and keep expanding the window until you reach the desired number of snps.
\code{subset_DT} should only contain one locus from one chromosome.
}
\seealso{
Other SNP filters: 
\code{\link{filter_snps}()},
\code{\link{gene_trimmer}()},
\code{\link{snps_to_condition}()},
\code{\link{subset_common_snps}()}
}
\concept{SNP filters}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORCES_2020.R
\name{CORCES_2020.prepare_bulkATAC_peak_overlap}
\alias{CORCES_2020.prepare_bulkATAC_peak_overlap}
\title{Prepare data to plot overlap between datatable of SNPs and
cell-type-specific epigenomic peaks and coaccessibility data.}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.prepare_bulkATAC_peak_overlap(
  merged_dat,
  FDR_filter = NULL,
  snp_filter = "Consensus_SNP==T",
  add_HiChIP_FitHiChIP = T,
  annotate_genes = F,
  return_counts = T,
  verbose = T
)
}
\description{
Prepare data to plot overlap between datatable of SNPs and
cell-type-specific epigenomic peaks and coaccessibility data.
}
\examples{
data("merged_DT");
merged_dat <- subset(merged_DT, Consensus_SNP)
dat_melt <- CORCES_2020.prepare_bulkATAC_peak_overlap(merged_dat=merged_dat)
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloaders.R
\name{wget}
\alias{wget}
\title{wget}
\usage{
wget(
  input_url,
  output_path,
  background = T,
  force_overwrite = F,
  quiet = F,
  show_progress = T,
  continue = T,
  check_certificates = F,
  conda_env = "echoR"
)
}
\description{
R wrapper for wget
}
\seealso{
Other downloaders: 
\code{\link{axel}()},
\code{\link{downloader}()}
}
\concept{downloaders}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{SUMMARISE.get_CS_counts}
\alias{SUMMARISE.get_CS_counts}
\title{Tally tool-specific and union CS sizes}
\usage{
SUMMARISE.get_CS_counts(merged_dat, top_CS_only = F)
}
\description{
Tally tool-specific and union CS sizes
}
\examples{
data("merged_DT");
locus_order <- SUMMARISE.get_CS_counts(merged_dat=merged_DT)
}
\seealso{
Other summarise: 
\code{\link{SUMMARISE.CS_bin_plot}()},
\code{\link{SUMMARISE.CS_counts_plot}()},
\code{\link{SUMMARISE.get_CS_bins}()},
\code{\link{SUMMARISE.get_SNPgroup_counts}()},
\code{\link{SUMMARISE.peak_overlap_plot}()},
\code{\link{SUMMARISE.plot_dataset_overlap}()},
\code{\link{results_report}()},
\code{\link{super_summary_plot}()}
}
\concept{summarise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createDT.R
\name{createDT}
\alias{createDT}
\title{Interactive DT}
\usage{
createDT(DF, caption = "", scrollY = 400)
}
\description{
Generate an interactive data table with download buttons.
}
\seealso{
Other general: 
\code{\link{createDT_html}()},
\code{\link{dt.replace}()},
\code{\link{example_fullSS}()},
\code{\link{fillNA_CS_PP}()},
\code{\link{get_os}()},
\code{\link{get_sample_size}()},
\code{\link{tryFunc}()}
}
\concept{general}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.locusName_handler}
\alias{PAINTOR.locusName_handler}
\title{Construct locus name for PAINTOR}
\usage{
PAINTOR.locusName_handler(
  locus_name = NULL,
  locus,
  GWAS_datasets = NULL,
  QTL_datasets = NULL
)
}
\description{
Construct locus name for PAINTOR
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Roadmap.R
\name{ROADMAP.track_plot}
\alias{ROADMAP.track_plot}
\title{Plot Roadmap query}
\usage{
ROADMAP.track_plot(
  grl.roadmap.filt,
  gr.snp = NULL,
  geom = "density",
  adjust = 0.2,
  show_plot = T,
  as.ggplot = T
)
}
\arguments{
\item{grl.roadmap.filt}{Roadmap query results.}

\item{gr.snp}{Optionally, can include an extra \code{\link[GenomicRanges]{GRanges}} object
to ensure the plot does not extend beyond certain coordinates.}

\item{geom}{The type of plot to create.
Options include "density" and "histogram".}

\item{adjust}{The granularity of the peaks.}

\item{show_plot}{Whether to print the plot.}
}
\description{
Plot Roadmap query
}
\examples{
\dontrun{
data("BST1")
finemap_DT <- BST1
gr.snp <- DT_to_GRanges(finemap_DT)
grl.roadmap <- ROADMAP.query(results_path="./Roadmap", gr.snp=gr.snp, keyword_query="monocyte")
grl.roadmap.filt <- ROADMAP.merge_and_process_grl(grl.roadmap=grl.roadmap, gr.snp=gr.snp, n_top_tissues=5)
track.roadmap <- ROADMAP.track_plot(grl.roadmap.filt, gr.snp=gr.snp)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{biomart_geneInfo}
\alias{biomart_geneInfo}
\title{Get gene info using Biomart}
\source{
\href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
}
\usage{
biomart_geneInfo(geneList, reference_genome = "grch37")
}
\description{
Get gene info using Biomart
}
\examples{
gene_info <- biomart_geneInfo(c("PTK2B","CLU","APOE"))
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{annotation_file_name}
\alias{annotation_file_name}
\title{Name annotation file}
\usage{
annotation_file_name(locus_dir, lib_name)
}
\description{
Name annotation file
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.construct_subset_vcf_name}
\alias{LD.construct_subset_vcf_name}
\title{Construct the path to vcf subset}
\usage{
LD.construct_subset_vcf_name(
  subset_DT,
  LD_reference = NULL,
  locus_dir,
  whole_vcf = F
)
}
\description{
Construct the path to vcf subset
}
\examples{
data("locus_dir"); data("BST1");
vcf_subset <- LD.construct_subset_vcf_name(subset_DT=BST1, locus_dir=locus_dir, LD_reference="1KGlocal")
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dey_DeepLearning.R
\name{DEEPLEARNING.query_one_chr}
\alias{DEEPLEARNING.query_one_chr}
\title{Query deep learning annotations and LDscores}
\usage{
DEEPLEARNING.query_one_chr(
  subset_DT,
  base_url = "/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
  level = c("Variant_Level", "Allelic_Effect"),
  tissue = c("NTS", "Blood", "Brain"),
  model = c("Basenji", "BiClassCNN", "DeepSEA", "ChromHMM", "Roadmap", "Others"),
  mean_max = c("MEAN", "MAX"),
  type = c("annot", "ldscore"),
  nThread = 4,
  verbose = T
)
}
\description{
Query deep learning annotations and LDscores, and
then merge with \code{subset_DT} by \emph{SNP}.
}
\examples{
\dontrun{
data("BST1")
annot.dat <- DEEPLEARNING.query_one_chr(subset_DT=BST1, tissue="NTS", model="Basenji", type="annot")
}
}
\seealso{
Other DEEPLEARNING: 
\code{\link{DEEPLEARNING.query_multi_chr}()},
\code{\link{DEEPLEARNING.query}()}
}
\concept{DEEPLEARNING}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COJO.R
\name{get_stepwise_results}
\alias{get_stepwise_results}
\title{Gather stepwise conditional results}
\usage{
get_stepwise_results(cojo_dir)
}
\description{
Gather and preprocess the results of the \emph{GCTA-COJO} conditional stepwise procedure.
}
\seealso{
Other COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}: 
\code{\link{COJO.stepwise}()},
\code{\link{COJO}()},
\code{\link{process_COJO_results}()}
}
\concept{COJO
\url{https://www.nature.com/articles/ng.2213}
\url{https://www.cell.com/ajhg/fulltext/S0002-9297(10)00598-7}
\url{https://cnsgenomics.com/software/gcta/#Overview}}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.locus.R
\name{PLOT.add_multitrack_lines}
\alias{PLOT.add_multitrack_lines}
\title{Add vertical lines}
\usage{
PLOT.add_multitrack_lines(
  TRKS,
  finemap_dat,
  snp_groups = c("Lead", "UCS", "Consensus"),
  line_alpha = 1,
  line_size = 0.3,
  remove_duplicated_UCS_Consensus = T,
  verbose = F
)
}
\description{
Adds vertical lines indicates key SNPs
}
\seealso{
Other plot: 
\code{\link{PLOT.dot_summary}()},
\code{\link{PLOT.get_window_limits}()},
\code{\link{PLOT.get_window_suffix}()}
}
\concept{plot}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.download_ref_files}
\alias{POLYFUN.download_ref_files}
\title{Download 1000 Genomes reference files}
\usage{
POLYFUN.download_ref_files(
 
    alkes_url = "https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_plinkfiles.tgz",
  results_dir = "./results",
  force_overwrite = F,
  download_method = "wget",
  conda_env = "echoR"
)
}
\description{
Download 1000 Genomes reference files
}
\examples{
\dontrun{
ref.prefix <- POLYFUN.download_ref_files(force_overwrite=T)
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifbreakR.R
\name{MOTIFBREAKR.filter}
\alias{MOTIFBREAKR.filter}
\title{Summarise \code{\link{motifbreakR}} + \code{\link{echolocatoR}} results}
\usage{
MOTIFBREAKR.filter(
  merged_DT,
  mb.results,
  pct_threshold = NULL,
  pvalue_threshold = 1e-04,
  qvalue_threshold = 0.05,
  effect_strengths = c("strong"),
  snp_filter = "Consensus_SNP==T",
  top_TF_hits = F,
  no_no_loci = NULL,
  verbose = T
)
}
\description{
For each SNP we have at least one allele achieving a p-value below 1e-4 threshold that we required.
The seqMatch column shows what the reference genome sequence is at that location,
with the variant position appearing in an uppercase letter.
pctRef and pctAlt display the the score for the motif in the sequence
as a percentage of the best score that motif could achieve on an ideal sequence.
In other words (scoreVariant−minscorePWM)/(maxscorePWM−minscorePWM).
We can also see the absolute scores for our method in scoreRef and scoreAlt
and thier respective p-values.
}
\examples{
\dontrun{
data("merged_DT")
microglia_TF <- read.csv("~/Desktop/Fine_Mapping/resources/microglia_TF.csv")

root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR"
# mb.results <- readRDS(file.path(root, "motifbreakR_results.rds"))
mb.results <- readRDS("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR/motifbreakR_results.rds")

mb.lrrk2 <- readRDS(file.path(root, "motifbreakR_results.p_values_LRRK2.rds"))
mb.encode <- readRDS(file.path(root, "motifbreakR_results.encode.lrrk2.rds"))
mb.DYRK1A_FCGR2A <- readRDS(file.path(root,"mb.results_p.DYRK1A_FCGR2A.RDS"))
mb.MED12L<- readRDS(file.path(root,"MED12L.pvalues.RDS"))

mb.sub <- subset(mb.results, SNP_id \%in\% subset(merged_DT, Locus=="DNAH17" & (leadSNP | Consensus_SNP))$SNP)
mb.DNAH17 <- motifbreakR::calculatePvalue(results=mb.sub); saveRDS(mb.DNAH17, file.path(root,"DNAH17.pvalues.RDS"))
mb.DNAH17 <- readRDS(file.path(root,"DNAH17.pvalues.RDS"))

mb.MBNL2 <- readRDS(file.path(root, "motifbreakR_results.p_values_MBNL2.rds"))

MOTIFBREAKR.filter(merged_DT, mb.lrrk2, pct_threshold=NULL, effect_strengths=NULL)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.enrichment}
\alias{XGR.enrichment}
\title{XGR enrichment}
\usage{
XGR.enrichment(
  gr,
  merged_dat,
  foreground_filter = "Consensus_SNP==T",
  background_filter = NULL,
  grouping_vars = c("Study", "Assay", "Cell_type"),
  fg_sample_size = NULL,
  bg_sample_size = NULL,
  background.annotatable.only = F,
  verbose = T
)
}
\description{
XGR enrichment
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"

merged_dat <- merge_finemapping_results(dataset = dirname(root), LD_reference = "UKB", minimum_support = 0)
merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)

gr.merged <- merge_celltype_specific_epigenomics()
grouping_vars <- c("Study","Cell_type","Assay")

enrich.lead <- XGR.enrichment(gr=gr, merged_dat=merged_dat, foreground_filter="leadSNP==T",  grouping_vars=grouping_vars)
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_topConsensus.R
\name{find_topConsensus}
\alias{find_topConsensus}
\title{Find the top Consensus SNP}
\usage{
find_topConsensus(dat, top_N = 1, grouping_vars = c("Locus"))
}
\description{
Identify the \code{top_N} Consensus SNP(s) per Locus,
defined as the Consensus SNPs with the highest mean PP across all fine-mapping tools used.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.prepare_LD}
\alias{PAINTOR.prepare_LD}
\title{Prepare LD file for PAINTOR}
\usage{
PAINTOR.prepare_LD(
  subset_path,
  PT_results_path,
  locus_name,
  locus_DT,
  locus,
  LD_matrix = NULL
)
}
\description{
Prepare LD file for PAINTOR
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Roadmap.R
\name{ROADMAP.query}
\alias{ROADMAP.query}
\title{Query Roadmap by genomic coordinates}
\usage{
ROADMAP.query(
  results_path,
  gr.snp,
  keyword_query = NULL,
  limit_files = NULL,
  conda_env = "echoR",
  nThread = 4,
  verbose = T
)
}
\arguments{
\item{results_path}{Where to store query results.}

\item{gr.snp}{\code{\link[GenomicRanges]{GRanges}} object of SNPs to query Roadmap with.}

\item{limit_files}{Limit the number of annotation files queried (for faster testing).}
}
\description{
Query Roadmap by genomic coordinates
}
\examples{
\dontrun{
data("BST1")
grl.roadmap <- ROADMAP.query(results_path="./Roadmap", gr.snp=BST1, keyword_query="placenta")
}
}
\seealso{
Other ROADMAP: 
\code{\link{ROADMAP.construct_reference}()},
\code{\link{ROADMAP.merge_and_process_grl}()},
\code{\link{ROADMAP.query_and_plot}()},
\code{\link{ROADMAP.tabix}()}
}
\concept{ROADMAP}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.download_and_standardize}
\alias{XGR.download_and_standardize}
\title{Download, standardize, and merge XGR annotations}
\usage{
XGR.download_and_standardize(
  lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes", "TFBS_Conserved",
    "Uniform_TFBS"),
  as_GRangesList = F,
  finemap_dat,
  nCores = 4
)
}
\arguments{
\item{lib.selections}{Which XGR annotations to check overlap with.
For full list of libraries see:
 \url{http://xgr.r-forge.r-project.org/#annotations-at-the-genomic-region-level}}

\item{as_GRangesList}{Return as a \code{GRangesList}, instead of a single merged \code{GRanges} object.}
}
\value{
GRangesList
}
\description{
Merges a list of XGR annotations into a single GRanges object
}
\examples{
\dontrun{
data("BST1")
gr.lib <- XGR.download_and_standardize(lib.selections=c("ENCODE_DNaseI_ClusteredV3_CellTypes"), finemap_dat=BST1, nCores=1)
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment_plot}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_topSNPs.R
\name{import_topSNPs}
\alias{import_topSNPs}
\title{Import top GWAS/QTL summary statistics}
\usage{
import_topSNPs(
  topSS,
  show_table = T,
  sheet = 1,
  munge = T,
  ref_genome = NULL,
  chrom_col = NULL,
  position_col = NULL,
  min_POS_col = NULL,
  max_POS_col = NULL,
  snp_col = NULL,
  pval_col = NULL,
  effect_col = NULL,
  locus_col = "Locus",
  grouping_vars = c("Locus"),
  gene_col = "Gene",
  remove_variants = NULL,
  nThread = 4,
  verbose = T
)
}
\arguments{
\item{show_table}{Create an interative data table.}

\item{sheet}{If the \emph{topSS} file is an excel sheet, you can specify which tab to use.
You can provide either a number to identify the tab by order,
or a string to identify the tab by name.}

\item{chrom_col}{Name of the chromosome column in the full summary stats file.
Can be "chr1" or "1" format.
(\emph{default: ="CHR"})}

\item{position_col}{Name of the genomic position column in the full summary stats file.
Must be in units of basepairs.
(\emph{default: ="POS"})}

\item{snp_col}{Name of the SNP RSID column in the full summary stats file.
(\emph{default: ="SNP"})}

\item{pval_col}{Name of the p-value column in the full summary stats file.
Raw p-values are preferred, but if not available corrected p-values (e.g. FDR) can be used instead.
(\emph{default: ="P"})}

\item{effect_col}{Name of the effect size column in the full summary stats file.
Effect size is preferred, but if not available other metrics like Beta for Odds Ratio can be used instead.
(\emph{default: ="Effect"})}

\item{locus_col}{Name of the locus column in the full summary stats file.
(\emph{default: ="Locus"})}

\item{grouping_vars}{The variables that you want to group by
such that each grouping_var combination has its own index SNP.
For example, if you want one index SNP per QTL eGene - GWAS locus pair, you could supply:
\code{grouping_vars=c("Locus","Gene")}.}

\item{gene_col}{An optional column to keep track of the causal gene(s) in each locus (if known).}

\item{remove_variants}{A list of variants to remove from the locus subset file.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}

\item{top_SNPs}{Can be a data.frame with the top sumary stats per locus.
Alternatively, you can provide a path to the stored top summary stats file.
Can be in any tabular format (e.g. excel, .tsv, .csv, etc.).
This file should have one lead GWAS/QTL hits per locus.
If there is more than one SNP per locus, the one with the smallest p-value (then the largest effect size) is selected as the lead SNP.
The lead SNP will be used as the center of the locus when constructing the locus subset files.}
}
\description{
The resulting top_SNPs data.frame can be used to guide
the \code{\link{finemap_loci}} in querying and fine-mapping loci.
}
\seealso{
Other query functions: 
\code{\link{TABIX.convert_file}()},
\code{\link{TABIX.query}()},
\code{\link{TABIX}()},
\code{\link{extract_SNP_subset}()},
\code{\link{query_by_coordinates_merged}()},
\code{\link{query_by_coordinates}()},
\code{\link{query_by_gene}()},
\code{\link{query_by_probe}()},
\code{\link{query_fullSS}()},
\code{\link{query_handler}()}
}
\concept{query functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CORCES_2020.bulkATACseq_peaks}
\alias{CORCES_2020.bulkATACseq_peaks}
\title{bulkATACseq peaks from human brain tissue}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 186559 rows and 10 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.bulkATACseq_peaks
}
\description{
Each row represents an individual peak identified in the bulk ATAC-seq data.
}
\details{
Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
Specifically: \emph{STable2_Features_bulkATAC-seq_Peaks}
}
\examples{
\dontrun{
dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable2_Features_bulkATAC-seq_Peaks.xlsx", skip = 18)
CORCES_2020.bulkATACseq_peaks <- data.table::data.table(dat)
usethis::use_data(CORCES_2020.bulkATACseq_peaks)
}
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_celltype_peaks}},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.1KG}
\alias{LD.1KG}
\title{Compute LD from 1000 Genomes}
\usage{
LD.1KG(
  locus_dir,
  subset_DT,
  LD_reference = "1KGphase1",
  superpopulation = "EUR",
  vcf_folder = NULL,
  remote_LD = T,
  LD_block = F,
  LD_block_size = 0.7,
  remove_correlates = F,
  remove_tmps = T,
  fillNA = 0,
  download_method = "wget",
  nThread = 4,
  conda_env = "echoR",
  verbose = T
)
}
\arguments{
\item{LD_reference}{Which linkage disequilibrium reference panel do you want to use.
Options include:
\describe{
\item{"UKB"}{A pre-caclulated LD reference matrix from a subset of caucasian British individuals from the UK Biobank. See \href{https://www.biorxiv.org/content/10.1101/807792v2}{Wiessbrod et al. (2019)} for more details.}
\item{"1KGphase1"}{Download a subset of the 1000 Genomes Project Phase 1 vcf and calculate LD on the fly with plink.}
\item{"1KGphase3"}{Download a subset of the 1000 Genomes Project Phase 3 vcf and calculate LD on the fly with plink.}
\item{"<path>/*.vcf" or "<path>/*.vcf.gz"}{Alternatively, users can provide their own custom panel by supplying a list of \emph{.vcf} file path (one per locus) which \pkg{echolocatoR} will use to compute LD (using \emph{plink}).}
}}

\item{superpopulation}{Subset your LD reference panel by superopulation.
Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
\href{https://www.internationalgenome.org/faq/which-populations-are-part-your-study/}{1KGphase1 options} include:
\describe{
\item{"AFR"}{African [descent]}
\item{"AMR"}{Ad-mixed American}
\item{"EAS"}{East Asian}
\item{"EUR"}{European}
\item{"SAS"}{South Asian}
}}

\item{remote_LD}{When acquiring LD matrices,
the default is to delete the full vcf or npz files after \pkg{echolocatoR} has extracted the necssary subset.
However, if you wish to keep these full files (which can be quite large) set \code{remote_LD=T}.}

\item{LD_block}{Calculate LD blocks with \emph{plink} and only include the block to which the lead SNP belongs.}

\item{LD_block_size}{Adjust the granularity of block sizes when \code{LD_block=T}.}

\item{remove_correlates}{A named list, where the names are the RSIDs of SNPs
whose LD correlates you wish to remove,
and the value is the absolute r2 threshold you wish to filter at for each RSID respectively
(e.g. \code{ remove_correlates = c("rs76904798"=.2, "rs10000737"=.8)}).
This will also remove the SNPs in \code{remove_correlates} themselves.}

\item{remove_tmps}{Whether to remove any temporary files (e.g. FINEMAP output files) after the pipeline is done running.}

\item{fillNA}{When pairwise LD (r) between two SNPs is \code{NA}, replace with 0.}

\item{conda_env}{The name of a conda environment to use.}

\item{verbose}{Whether \pkg{echolocatoR} should be verbose or silent.}
}
\description{
Downloads a subset vcf of the 1KG database that matches your locus coordinates.
Then uses \emph{plink} to calculate LD on the fly.
}
\details{
This approach is taken, because other API query tools have limitations with the window size being queried.
This approach does not have this limitations, allowing you to fine-map loci more completely.
}
\examples{
\dontrun{
data("BST1"); data("locus_dir");
BST1 <- limit_SNPs(max_snps = 500, subset_DT = BST1)
LD_matrix <- LD.1KG(locus_dir=file.path("~/Desktop",locus_dir), subset_DT=BST1, LD_reference="1KGphase1")

## Kunkle et al 2019
locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ACE"

}
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{biomart_snps_to_geneInfo}
\alias{biomart_snps_to_geneInfo}
\title{Identify which genes SNPs belong to using Biomart}
\source{
\href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
}
\usage{
biomart_snps_to_geneInfo(snp_list, reference_genome = "grch37")
}
\description{
Identify which genes SNPs belong to using Biomart
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.plink_LD}
\alias{LD.plink_LD}
\title{Calculate LD}
\usage{
LD.plink_LD(
  leadSNP = NULL,
  subset_DT,
  bim_path = NULL,
  remove_excess_snps = T,
  merge_by_RSID = F,
  LD_folder,
  min_r2 = F,
  min_Dprime = F,
  remove_correlates = F,
  fillNA = 0,
  plink_prefix = "plink",
  verbose = T,
  conda_env = NULL
)
}
\description{
Use \emph{plink} to calculate LD from a vcf.
}
\examples{
locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ACE"
LD_folder <- file.path(locus_dir,"LD")
ld.matrix <- LD.plink_LD(subset_DT=BST1, LD_folder=LD_folder)
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{ANNOTATE.annotate_missense}
\alias{ANNOTATE.annotate_missense}
\title{Annotate any missense variants}
\usage{
ANNOTATE.annotate_missense(merged_dat, snp_filter = "Support>0")
}
\description{
Annotate any missense variants
}
\examples{
\dontrun{
data("merged_DT");
annotated_DT <- ANNOTATE.annotate_missense(merged_dat=merged_DT, snp_filter="Support>0")
}
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()},
\code{\link{merge_finemapping_results_each}()}
}
\concept{annotate}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.datatype_handler}
\alias{PAINTOR.datatype_handler}
\title{Determine data types (GWAS, QTL)}
\usage{
PAINTOR.datatype_handler(GWAS_datasets = NULL, QTL_datasets = NULL, locus)
}
\description{
Determine data types (GWAS, QTL)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN_SUSIE}
\alias{POLYFUN_SUSIE}
\title{Run PolyFun+SUSIE fine-mapping pipeline}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN_SUSIE(
  locus_dir,
  polyfun = NULL,
  finemap_dat = NULL,
  LD_matrix = NULL,
  polyfun_approach = "non-parametric",
  dataset_type = "GWAS",
  max_causal = 5,
  sample_size = NULL,
  server = F,
  PP_threshold = 0.95,
  conda_env = "echoR"
)
}
\description{
Uses echolocatoR wrapper for SUSIE instead of the \code{\link{POLYFUN.finemapper}}
function which uses a python script provided with PolyFun.
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()}
}
\concept{polyfun}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_bw_filt.R
\name{import.bw.filt}
\alias{import.bw.filt}
\title{Import a subset of a bigwig file
based on the coordinates in a GRanges object (gr.dat).}
\usage{
import.bw.filt(bw.file, gr.dat, full_data = T)
}
\arguments{
\item{bw.file}{Path to a bigwig file.}

\item{gr.dat}{GenomicRanges object to query the bigwig file with.}

\item{full_data}{Whether to return the actual read ranges (\code{full_data=T}),
or just the "score" column which summarizes the height of
the aggregated reads across the genome (\code{full_data=T}).}
}
\description{
Import a subset of a bigwig file
based on the coordinates in a GRanges object (gr.dat).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/directory.R
\name{Directory_info}
\alias{Directory_info}
\title{Retrieve the location of summary stats files}
\usage{
Directory_info(dataset_name, variable = "fullSS.local")
}
\description{
Retrieve the location of summary stats files
}
\concept{directory}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.prepare_snp_input}
\alias{POLYFUN.prepare_snp_input}
\title{Prepare SNP input for PolyFun}
\usage{
POLYFUN.prepare_snp_input(
  PF.output.path,
  locus_dir,
  finemap_dat = NULL,
  nThread = 1
)
}
\description{
PolyFun requires a space-delimited (gzipped or not) file with these columns:
\itemize{
\item{CHR}
\item{BP}
\item{A1}
\item{A2}
}
}
\examples{
data("BST1"); data("locus_dir");
finemap_dat <- BST1
PF.output.path <- file.path(locus_dir, "PolyFun")
POLYFUN.prepare_snp_input(PF.output.path=PF.output.path, locus_dir=locus_dir, finemap_dat=finemap_dat)
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMPACT.R
\name{IMPACT.get_annotations}
\alias{IMPACT.get_annotations}
\title{Download \emph{IMPACT} annotations}
\usage{
IMPACT.get_annotations(
  baseURL = "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/Annotations",
  chrom = NULL,
  subset_DT = NULL,
  nThread = 4,
  all_snps_in_range = F,
  verbose = T
)
}
\description{
Includes the raw annotation itself,
as well as per-SNP \emph{IMPACT} scores for each annotation.
}
\details{
Unfortunately, you have to download the entire chromosome file at once,
 because they aren't Tabix indexed. To minimize the memory load,
 this function only keeps the portion of the \emph{IMPACT} file that overlaps with the
 coordinates in \code{subset_DT}.
}
\examples{
\dontrun{
data("BST1")
annot_melt <- IMPACT.get_annotations(subset_DT=BST1)
}
}
\seealso{
Other IMPACT: 
\code{\link{IMPACT.get_annotation_key}()},
\code{\link{IMPACT.get_top_annotations}()},
\code{\link{IMPACT.iterate_get_annotations}()},
\code{\link{IMPACT.postprocess_annotations}()},
\code{\link{IMPACT.snp_group_boxplot}()},
\code{\link{IMPACT_annotation_key}},
\code{\link{IMPACT_heatmap}()},
\code{\link{prepare_mat_meta}()}
}
\concept{IMPACT}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpliceAI.R
\name{SPLICEAI.run}
\alias{SPLICEAI.run}
\title{Run pre-trained \emph{SpliceAI} model to get predictions}
\source{
\href{https://github.com/Illumina/SpliceAI}{GitHub}
\href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
}
\usage{
SPLICEAI.run(
  vcf_path = "./GWAS_converted.vcf",
  output_path = "spliceai_predictions.vcf",
  reference_fasta = "/pd-omics/tools/polyfun/reference_fasta/hg19.fa",
  gene_annotation = "./echolocatoR/tools/spliceAI/grch37.txt",
  distance = 50,
  mask = 0
)
}
\description{
Run pre-trained \emph{SpliceAI} model to get predictions
}
\seealso{
Other SpliceAI: 
\code{\link{SPLICEAI.plot}()},
\code{\link{SPLICEAI.snp_probs}()},
\code{\link{SPLICEAI.subset_precomputed_tsv_iterate}()},
\code{\link{SPLICEAI.subset_precomputed_tsv}()},
\code{\link{SPLICEAI.subset_precomputed_vcf}()}
}
\concept{SpliceAI}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.create_locusFile.QTL}
\alias{PAINTOR.create_locusFile.QTL}
\title{Create QTL locus file for PAINTOR}
\usage{
PAINTOR.create_locusFile.QTL(
  finemap_dat,
  GWAS_datasets,
  QTL_datasets = NULL,
  PT_results_path,
  locus_name,
  metric_suffix = ".t_stat",
  NA_method = c("fill", "drop"),
  force_new_zscore = F
)
}
\description{
Create QTL locus file for PAINTOR
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGR.R
\name{XGR.enrichment_plot}
\alias{XGR.enrichment_plot}
\title{Plot enrichment results}
\usage{
XGR.enrichment_plot(
  enrich_res,
  title = NULL,
  subtitle = NULL,
  facet_formula = NULL,
  line_formula = "y ~ x",
  line_method = "lm",
  line_span = 1,
  FDR_thresh = 1,
  plot_type = "bar",
  shape_var = "Cell_type",
  facet_scales = "free",
  show_plot = T,
  save_plot = F,
  height = 5,
  width = 5
)
}
\description{
Plot enrichment results
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
### merged enrichment results
enrich_res <- data.table::fread( file.path(root,"XGR/celltypespecific_epigenomics.SNP_groups.csv.gz"))
enrich_res <- data.table::fread( file.path(root,"XGR/celltypespecific_epigenomics.snp_groups.csv.gz"))
enrich_boot <- data.table::fread(file.path(root,"XGR/celltypespecific_epigenomics.snp_groups.permute.csv.gz"))
enrich_assay <- data.table::fread(file.path(root,"XGR/celltypespecific_epigenomics.snp_groups.assay.csv.gz"))

# Merged volcano plot
enrich_res <- subset(enrich_res, SNP_Group != "Consensus (-PolyFun)") \%>\% dplyr::rename(SNP_group=SNP_Group)
gp <- XGR.enrichment_plot(enrich_res=subset(enrich_res , !Assay \%in\%  c("HiChIP_FitHiChIP","PLAC")) , title="Enrichment: Cell-type-specific epigenomics", plot_type="point",save_plot=file.path(root,"XGR/celltypespecific_epigenomics.enrich_volcano.png"), height=6, width=8, shape_var="Assay")
## Merged bar plot
gp <- XGR.enrichment_plot(enrich_res=enrich_res,  plot_type="bar", facet_formula=".~Assay",FDR_thresh=.05)


# Merged volcano plot (permuted)
gp <- XGR.enrichment_plot(enrich_res=enrich.scATAC.permute, title="Permuted enrichment: Cell-type-specific peaks and elements", plot_type="point")
}
}
\seealso{
Other XGR: 
\code{\link{DT_to_GRanges}()},
\code{\link{GRanges_to_BED}()},
\code{\link{XGR.download_and_standardize}()},
\code{\link{XGR.enrichment_bootstrap}()},
\code{\link{XGR.enrichment}()},
\code{\link{XGR.filter_assays}()},
\code{\link{XGR.filter_sources}()},
\code{\link{XGR.import_annotations}()},
\code{\link{XGR.iterate_enrichment}()},
\code{\link{XGR.iterate_overlap}()},
\code{\link{XGR.merge_and_process}()},
\code{\link{XGR.plot_enrichment}()},
\code{\link{XGR.plot_peaks}()},
\code{\link{XGR.prepare_foreground_background}()}
}
\concept{XGR}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{top_SNPs}
\alias{top_SNPs}
\title{TopSS example file (processed)}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 97 rows and 7 columns.
}
\usage{
top_SNPs
}
\description{
Summary stats of the top SNP(s) per locus.
Used to query locus subsets.for fine-mapping.
}
\examples{
data("Nalls_top_SNPs")
top_SNPs <- import_topSNPs(topSS=Nalls_top_SNPs, chrom_col="CHR", position_col="BP", snp_col="SNP", pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene", group_by_locus=T, locus_col="Nearest Gene", remove_variants="rs34637584")
\dontrun{
usethis::use_data(top_SNPs, overwrite=T)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.LD_blocks}
\alias{LD.LD_blocks}
\title{Calculate LD blocks.}
\usage{
LD.LD_blocks(LD_folder, LD_block_size = 0.7)
}
\description{
Uses \emph{plink} to group highly correlated SNPs into LD blocks.
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.save_LD_matrix}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.h2_enrichment_SNPgroups}
\alias{POLYFUN.h2_enrichment_SNPgroups}
\title{Run heritability enrichment tests across SNP groups}
\source{
https://www.biorxiv.org/content/10.1101/807792v3
}
\usage{
POLYFUN.h2_enrichment_SNPgroups(
  merged_dat,
  chrom = "*",
 
    ldsc_dir = "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output",
  ldsc_suffix = "*.snpvar_constrained.gz",
  save_enrich = F,
  nThread = 4
)
}
\description{
Run heritability enrichment tests across SNP groups
}
\examples{
\dontrun{
root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
# IMPORTANT! For this to make sense, you need to merge the full data ("merged_DT" only includes Support>0 and leadSNPs)
merged_dat <- merge_finemapping_results(dataset = dirname(root), LD_reference = "UKB", minimum_support = 0)
merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)

RES <- POLYFUN.h2_enrichment_SNPgroups(merged_dat=merged_dat, ldsc_dir=file.path(root,"PolyFun/output"),  save_enrich=file.path(root,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.csv.gz"))
}
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.initialize}
\alias{POLYFUN.initialize}
\title{Create output dir and import SNP data.frame}
\usage{
POLYFUN.initialize(locus_dir, finemap_dat = NULL, nThread = 4)
}
\description{
Create output dir and import SNP data.frame
}
\examples{
data("BST1"); data("locus_dir");
finemap_DT <- POLYFUN.initialize(locus_dir=locus_dir, finemap_dat=BST1)
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.help}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.get_lead_r2}
\alias{LD.get_lead_r2}
\title{Find correlates of the lead GWAS/QTL SNP}
\usage{
LD.get_lead_r2(
  finemap_dat,
  LD_matrix = NULL,
  fillNA = 0,
  LD_format = "matrix",
  verbose = T
)
}
\description{
Find correlates of the lead GWAS/QTL SNP
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAINTOR.R
\name{PAINTOR.import_QTL_DT}
\alias{PAINTOR.import_QTL_DT}
\title{Import QTL data for PAINTOR}
\usage{
PAINTOR.import_QTL_DT(
  QTL_datasets,
  locus,
  trim_gene_limits = F,
  force_new_subset = F,
  metric = "t_stat"
)
}
\description{
@keywords internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.locus.R
\name{PLOT.transcript_model_track}
\alias{PLOT.transcript_model_track}
\title{Plot gene/transcript models}
\source{
\href{https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html}{ensembld tutorial}
\href{https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html#45_GeneRegionTrack}{Gvix tutorial}
\href{http://bioconductor.org/packages/devel/bioc/vignettes/ggbio/inst/doc/ggbio.pdf}{ggbio tutorial}
}
\usage{
PLOT.transcript_model_track(
  finemap_dat,
  max_transcripts = 1,
  remove_pseudogenes = T,
  show.legend = T,
  show_plot = F,
  fill = "skyblue",
  shape = c("arrow", "box", "ellipse", "smallArrow"),
  transcriptAnnotation = c("symbol", "transcript"),
  collapseTranscripts = c(F, T, "longest"),
  stacking = c("squish", "hide", "dense", "pack", "full"),
  method = "ggplot",
  xtext = T,
  expand_x_mult = c(0.05, 0.1),
  verbose = T
)
}
\description{
Plot gene/transcript models
}
\examples{
data("LRRK2")
gene_track <- PLOT.transcript_model_track(finemap_dat=LRRK2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POLYFUN.R
\name{POLYFUN.help}
\alias{POLYFUN.help}
\title{Display PolyFun help}
\usage{
POLYFUN.help(polyfun = NULL, conda_env = "echoR")
}
\description{
Display PolyFun help
}
\examples{
POLYFUN.help()
}
\seealso{
Other polyfun: 
\code{\link{POLYFUN.compute_priors}()},
\code{\link{POLYFUN.download_ref_files}()},
\code{\link{POLYFUN.find_polyfun_folder}()},
\code{\link{POLYFUN.finemapper}()},
\code{\link{POLYFUN.functional_enrichment}()},
\code{\link{POLYFUN.gather_annotations}()},
\code{\link{POLYFUN.gather_ldscores}()},
\code{\link{POLYFUN.get_precomputed_priors}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups_plot}()},
\code{\link{POLYFUN.h2_enrichment_SNPgroups}()},
\code{\link{POLYFUN.h2_enrichment}()},
\code{\link{POLYFUN.initialize}()},
\code{\link{POLYFUN.ldsc_annot_enrichment}()},
\code{\link{POLYFUN.munge_summ_stats}()},
\code{\link{POLYFUN.plot}()},
\code{\link{POLYFUN.prepare_snp_input}()},
\code{\link{POLYFUN.read_parquet}()},
\code{\link{POLYFUN.run_ldsc}()},
\code{\link{POLYFUN_SUSIE}()}
}
\concept{polyfun}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifbreakR.R
\name{MOTIFBREAKR.filter_by_metadata}
\alias{MOTIFBREAKR.filter_by_metadata}
\title{Filter by motif database metadata}
\usage{
MOTIFBREAKR.filter_by_metadata(mb.results, Organism = "Hsapiens")
}
\description{
Filter by motif database metadata
}
\seealso{
Other motifbreakR: 
\code{\link{MOTIFBREAKR}()}
}
\concept{motifbreakR}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NOTT_2019.superenhancer_interactome}
\alias{NOTT_2019.superenhancer_interactome}
\title{Brain cell type-specific interactomes with superenhancers}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 2954 rows and 29 columns.
}
\source{
\url{https://science.sciencemag.org/content/366/6469/1134}
}
\usage{
NOTT_2019.superenhancer_interactome
}
\description{
Originally from \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}.
Specifically: \emph{aay0793-Nott-Table-S6.xlsx}.
}
\examples{
\dontrun{
NOTT_2019.superenhancer_interactome <- data.table::data.table(readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S6.xlsx", skip=2)  )
usethis::use_data(NOTT_2019.superenhancer_interactome)
}
}
\seealso{
Other NOTT_2019: 
\code{\link{NOTT_2019.bigwig_metadata}},
\code{\link{NOTT_2019.epigenomic_histograms}()},
\code{\link{NOTT_2019.get_epigenomic_peaks}()},
\code{\link{NOTT_2019.get_interactions}()},
\code{\link{NOTT_2019.get_interactome}()},
\code{\link{NOTT_2019.get_promoter_celltypes}()},
\code{\link{NOTT_2019.get_promoter_interactome_data}()},
\code{\link{NOTT_2019.get_regulatory_regions}()},
\code{\link{NOTT_2019.interactome}},
\code{\link{NOTT_2019.plac_seq_plot}()},
\code{\link{NOTT_2019.superenhancers}()}
}
\concept{NOTT_2019}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finemapping_portal_API.R
\name{GITHUB.make_data_dict}
\alias{GITHUB.make_data_dict}
\title{echolocatoR Fine-mapping portal: data dictionary}
\usage{
GITHUB.make_data_dict(named_lists)
}
\description{
echolocatoR Fine-mapping portal: data dictionary
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD.save_LD_matrix}
\alias{LD.save_LD_matrix}
\title{Save LD_matrix}
\usage{
LD.save_LD_matrix(
  LD_matrix,
  subset_DT,
  locus_dir,
  fillNA = 0,
  LD_reference,
  subset_common = T,
  sparse = T,
  verbose = T
)
}
\description{
Save LD_matrix
}
\examples{
data("BST1"); data("LD_matrix"); data("locus_dir");
LD_list <- LD.save_LD_matrix(LD_matrix=LD_matrix, subset_DT=BST1, locus_dir=file.path("~/Desktop",locus_dir), LD_reference="UKB")
LD_list <- LD.save_LD_matrix(LD_matrix=LD_matrix, subset_DT=BST1, locus_dir=file.path("~/Desktop",locus_dir), LD_reference="custom_vcf")
}
\seealso{
Other LD: 
\code{\link{LD.1KG_download_vcf}()},
\code{\link{LD.1KG}()},
\code{\link{LD.LD_blocks}()},
\code{\link{LD.UKBiobank}()},
\code{\link{LD.calculate_LD}()},
\code{\link{LD.construct_subset_vcf_name}()},
\code{\link{LD.custom_panel}()},
\code{\link{LD.dprime_table}()},
\code{\link{LD.filter_LD}()},
\code{\link{LD.filter_vcf_gaston}()},
\code{\link{LD.filter_vcf}()},
\code{\link{LD.get_locus_vcf_folder}()},
\code{\link{LD.index_vcf}()},
\code{\link{LD.leadSNP_block}()},
\code{\link{LD.load_or_create}()},
\code{\link{LD.plink_LD}()},
\code{\link{LD.plink_file}()},
\code{\link{LD.plot_LD}()},
\code{\link{LD.query_vcf}()},
\code{\link{LD.rds_to_npz}()},
\code{\link{LD.read_bin}()},
\code{\link{LD.read_ld_table}()},
\code{\link{LD.run_plink_LD}()},
\code{\link{LD.snpstats_get_LD}()},
\code{\link{LD.snpstats_get_MAF}()},
\code{\link{LD.translate_population}()},
\code{\link{LD.vcf_to_bed}()},
\code{\link{LDlinkR.LDproxy_batch}()},
\code{\link{popDat_1KGphase1}},
\code{\link{popDat_1KGphase3}},
\code{\link{saveSparse}()}
}
\concept{LD}
\keyword{internal}
