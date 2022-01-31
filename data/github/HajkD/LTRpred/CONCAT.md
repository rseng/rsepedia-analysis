---
title: 'LTRpred: _de novo_ annotation of intact retrotransposons'
authors:
- affiliation: "1, 2"
  name: Hajk-Georg Drost
  orcid: 0000-0002-1567-306X
date: "28 February 2020"
output: pdf_document
bibliography: paper.bib
tags:
- R
- retrotransposon annotation
- de novo functional annotation of retrotransposons
- meta-genome scale retrotransposon annotation
- transposable elements
- mobile transposable elements
- genome dynamics
affiliations:
- index: 1
  name: The Sainsbury Laboratory, University of Cambridge, Bateman Street, Cambridge
    CB2 1LR, UK
- index: 2
  name: Department of Molecular Biology, Max Planck Institute for Developmental Biology,
    72076 Tuebingen, Germany
---

# Summary

Transposable elements (TEs) play a crucial role in altering the genomic landscape of all organisms and thereby massively influence the genetic information passed on to succeeding generations [@Sundaram2020-yw]. In the past, TEs were seen as selfish mobile elements populating host genomes to increase their chances for  transgenerational transmission over long evolutionary time scales. This notion of selfish elements is slowly changing [@Drost2019-rz] and a new picture drawing a complex genetic landscape benefitting both, host and TE, emerges whereby novel forms can arise through random shuffling of genetic material. For example, the tomato fruit shape [@Benoit2019-ux], moth adaptive cryptic coloration that occurred during the industrial revolution [@Chuong2017-cb], and inner cell mass development in human embryonic stem cells [@Chuong2017-cb] were all shown to be driven by TE activity. Thus, the impact of these elements on altering morphological traits is imminent and requires new attention in the light of evolvability. However, TEs tend to degenerate their sequence leaving their fragmented copies considered as _junk_ DNA in host genomes, which hamper assembly and annotation of new genomes.

Nowadays, the _de novo_ detection of transposable elements is performed by annotation tools specifically designed to capture any type of repeated sequence, TE family, or remnant DNA loci that can be associated with known transposable elements within a genome assembly. The main goal of such efforts is to retrieve a maximum number of loci that can be associated with known TEs. If successful, such annotation can then be used to mask host genomes from TE remnants to simplify genomics studies focusing on host genes. Therefore, there is no automatically performed distinction between complete and potentially active TE and their mutated copies.

Here, we introduce the ``LTRpred`` pipeline which allows to _de novo_ annotate functional and thus potentially mobile retrotransposons in any given genome assembly. Different from other annotation tools, ``LTRpred`` focuses on retrieving structurally intact elements within sequences of genomes rather than characterizing all traces of historic TE activity. 

Such functional annotation is most useful when trying to spot retrotransposons responsible for recent reshuffling of genetic material in the tree of life. Detecting and further characterization of those active retrotransposons yields the potential to harness them as mutagenesis agents by inducing transposition bursts in a controlled fashion to stimulate genomic reshaping processes towards novel traits. 

# ``LTRpred`` emerged as valuable tool for diverse TE mobilization studies

``LTRpred`` was successfully used in previous studies to annotate functional retrotransposons for various applications. In detail, ``LTRpred`` was used to annotate the retrotransposon family ``RIDER`` within the plant kingdom, shown to be involved in tomato fruit shape elongation. Together with experimental evidence, our analyses revealed that RIDER elements can be activated via drought stress and may help plants rich in RIDER activity to better adapt to drought stress conditions [@Benoit2019-ux]. In a complementary study, ``LTRpred`` was used to generate a candidate list of potentially mobile retrotransposon families in rice and tomato, which were then confirmed to produce extrachromosomal DNA using the ALE-Seq methodology [@Cho2019-zp]. Finally, ``LTRpred`` supported efforts to annotate and date functional retrotransposons in tomato and _Arabidopsis_ which led to the finding that chromodomain DNA methyltransferases (CMTs) silence young and intact retrotransposons in distal chromatin whereas older non-functional retrotransposons are affected by small RNA-directed DNA methylation [@Wang2019].

Together, potentially functional retrotransposons annotated _de novo_ with ``LTRpred`` were subsequently shown to be active and mobile in diverse molecular studies. This approach may stimulate a new wave of research towards understanding the physiological role of functional retrotransposons and to reveal the mechanistic principles of transposon associated evolvability. 

# Pipeline dependencies

The LTRpred pipeline depends on the R packages biomartr [@Drost2017-cw], dplyr [@Wickham2019-kh], Biostrings, stringr [@Wickham2019-kh], IRanges [@Lawrence2013-dv], RColorBrewer, ggplot2 [@Wickham2016-eq], readr [@Wickham2019-kh], magrittr, downloader, BSDA, ggrepel, ggbio [@Yin2012-ro], and gridExtra.

In detail, the LTRpred pipeline calls the command line tools suffixerator, LTRharvest [@Ellinghaus2008-hu], and LTRdigest [@Steinbiss2009-vg], which are part of the GenomeTools library [@Gremme2013-ba] using customized parameter settings to screen for repeated LTRs, specific sequence motifs such as primer binding sites (PBS), polypurine tract motifs (PPT), and target site duplications (TSD) and for conserved protein domains such as reverse transcriptase (gag), integrase DNA binding domain, integrase Zinc binding domain, RNase H, and the integrase core domain. 
The LTRharvest and LTRdigest outputs are efficiently parsed by LTRpred and transformed into a tidy data format [@Wickham2019-kh] which subsequently enables automation of false positive curation. Next, open reading frame (ORF) prediction is performed by a customized wrapper function that runs the command line tool usearch [@Edgar2010-cb]. This step allows to automatically filter out retrotransposons that might have conserved protein domains such as an integrase or reverse transcriptase, but fail to have any ORFs and thus might not be expressed. In a third step, retrotransposon family clustering is performed using sequence clustering with vsearch [@Rognes2016-sk] which defines family members by >90% sequence homology of the full element to each other. In a fourth step, an automated hmmer search [@Finn2011-fs] against the Dfam database (Hubley et al., 2016) is performed to assign super-family associations such as Copia or Gypsy by comparing the protein domains of de novo predicted retrotransposons with already annotated TEs in the Dfam (https://dfam.org/home) database. In the last step, the _de novo_ annotated 5 prime and 3 prime LTR sequences are used to estimate the evolutionary age of the retrotransposon which should be treated with caution since retrotransposons can undergo reverse-transcriptase mediated recombination [@Sanchez2017-sy].

# Example workflow

After installing all prerequisite command line tools (https://hajkd.github.io/LTRpred/articles/Introduction.html#installation) users can run the `LTRpred()` pipeline using the default parameter configuration. In the following example, an LTR transposon prediction is performed for parts of the Human Y chromosome.

```r
# load LTRpred package
library(LTRpred)
# de novo LTR transposon prediction for the Human Y chromosome
LTRpred(
    genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"),
    cores = 4
)
```

The LTRpred() output table `*_LTRpred_DataSheet.tsv` is in tidy format and can then be imported using `read.ltrpred()`. The tidy output format is designed to work seamlessly with the tidyverse and R data science framework.

```r
# import LTRpred prediction output
Hsapiens_chrY <- read.ltrpred("Hsapiens_ChrY_ltrpred/
Hsapiens_ChrY_LTRpred_DataSheet.tsv")
# look at some results
dplyr::glimpse(Hsapiens_chrY)
```

```
Observations: 21
Variables: 92
$ species                 <chr> "Hsapiens_ChrY", "Hsapien...
$ ID                      <chr> "Hsapiens_ChrY_LTR_retrot...
$ dfam_target_name        <chr> NA, NA, NA, NA, NA, NA, N...
$ ltr_similarity          <dbl> 80.73, 89.85, 79.71, 83.2...
$ ltr_age_mya             <dbl> 0.7936246, 0.2831139, 0.7...
$ similarity              <chr> "(80,82]", "(88,90]", "(7...
$ protein_domain          <chr> "RVT_1", "RVT_1", NA, NA,...
$ orfs                    <int> 1, 1, 0, 0, 0, 0, 0, 1, 0...
$ chromosome              <chr> "NC000024.10Homosa", "NC0...
$ start                   <int> 3143582, 3275798, 3313536...
$ end                     <int> 3162877, 3299928, 3318551...
$ strand                  <chr> "-", "-", "+", "+", "-", ...
$ width                   <int> 19296, 24131, 5016, 12952...
$ annotation              <chr> "LTR_retrotransposon", "L...
$ pred_tool               <chr> "LTRpred", "LTRpred", "LT...
$ frame                   <chr> ".", ".", ".", ".", ".", ...
$ score                   <chr> ".", ".", ".", ".", ".", ...
$ lLTR_start              <int> 3143582, 3275798, 3313536...
$ lLTR_end                <int> 3143687, 3276408, 3313665...
$ lLTR_length             <int> 106, 611, 130, 126, 218, ...
$ rLTR_start              <int> 3162769, 3299338, 3318414...
$ rLTR_end                <int> 3162877, 3299928, 3318551...
$ rLTR_length             <int> 109, 591, 138, 137, 219, ...
$ lTSD_start              <int> 3143578, 3275794, 3313532...
$ lTSD_end                <int> 3143581, 3275797, 3313535...
$ lTSD_motif              <chr> "acag", "ttgt", "ttag", "...
$ rTSD_start              <int> 3162878, 3299929, 3318552...
$ rTSD_end                <int> 3162881, 3299932, 3318555...
$ rTSD_motif              <chr> "acag", "ttgt", "ttag", "...
$ PPT_start               <int> NA, NA, NA, NA, NA, 34660...
$ PPT_end                 <int> NA, NA, NA, NA, NA, 34660...
$ PPT_motif               <chr> NA, NA, NA, NA, NA, "agag...
$ PPT_strand              <chr> NA, NA, NA, NA, NA, "+", ...
$ PPT_offset              <int> NA, NA, NA, NA, NA, 23, N...
$ PBS_start               <int> NA, NA, 3313667, 3372512,...
$ PBS_end                 <int> NA, NA, 3313677, 3372522,...
$ PBS_strand              <chr> NA, NA, "+", "+", "-", "+...
$ tRNA                    <chr> NA, NA, "Homo_sapiens_tRN...
$ tRNA_motif              <chr> NA, NA, "aattagctgga", "c...
$ PBS_offset              <int> NA, NA, 1, 3, 0, 5, 2, 5,...
$ tRNA_offset             <int> NA, NA, 1, 0, 2, 5, 1, 5,...
$ `PBS/tRNA_edist`        <int> NA, NA, 1, 1, 1, 1, 1, 1,...
$ orf.id                  <chr> "NC000024.10Homosa_314358...
$ repeat_region_length    <int> 19304, 24139, 5024, 12960...
$ PPT_length              <int> NA, NA, NA, NA, NA, 27, N...
$ PBS_length              <int> NA, NA, 11, 11, 11, 11, 1...
$ dfam_acc                <chr> NA, NA, NA, NA, NA, NA, N...
$ dfam_bits               <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_e_value            <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_bias               <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_strand             <chr> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_modlen             <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_target_description <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Cluster           <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Target            <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Perc_Ident        <dbl> NA, NA, NA, NA, NA, NA, N...
$ Clust_cn                <int> NA, NA, NA, NA, NA, NA, N...
$ TE_CG_abs               <dbl> 62, 125, 35, 70, 139, 83,...
$ TE_CG_rel               <dbl> 0.003213101, 0.005180059,...
$ TE_CHG_abs              <dbl> 659, 830, 150, 396, 742, ...
$ TE_CHG_rel              <dbl> 0.03415216, 0.03439559, 0...
$ TE_CHH_abs              <dbl> 2571, 3454, 748, 1743, 31...
$ TE_CHH_rel              <dbl> 0.1332400, 0.1431354, 0.1...
$ TE_CCG_abs              <dbl> 13, 24, 6, 15, 33, 16, 4,...
$ TE_CCG_rel              <dbl> 0.0006737148, 0.000994571...
$ TE_N_abs                <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_3ltr_abs             <dbl> 1, 0, 1, 4, 8, 3, 1, 11, ...
$ CG_3ltr_rel             <dbl> 0.009433962, 0.000000000,...
$ CHG_3ltr_abs            <dbl> 2, 24, 9, 8, 14, 14, 9, 9...
$ CHG_3ltr_rel            <dbl> 0.01886792, 0.03927987, 0...
$ CHH_3ltr_abs            <dbl> 18, 69, 18, 26, 43, 70, 2...
$ CHH_3ltr_rel            <dbl> 0.16981132, 0.11292962, 0...
$ CCG_3ltr_abs            <dbl> 0, 0, 0, 2, 2, 0, 1, 4, 5...
$ CCG_3ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_3ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_5ltr_abs             <dbl> 1, 0, 1, 4, 8, 3, 1, 11, ...
$ CG_5ltr_rel             <dbl> 0.009433962, 0.000000000,...
$ CHG_5ltr_abs            <dbl> 2, 24, 9, 8, 14, 14, 9, 9...
$ CHG_5ltr_rel            <dbl> 0.01886792, 0.03927987, 0...
$ CHH_5ltr_abs            <dbl> 18, 69, 18, 26, 43, 70, 2...
$ CHH_5ltr_rel            <dbl> 0.16981132, 0.11292962, 0...
$ CCG_5ltr_abs            <dbl> 0, 0, 0, 2, 2, 0, 1, 4, 5...
$ CCG_5ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_5ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ cn_3ltr                 <dbl> NA, NA, NA, NA, NA, NA, N...
$ cn_5ltr                 <dbl> NA, NA, NA, NA, NA, NA, N..
```

### LTRpred output

The `LTRpred()` function internally generates a folder named `*_ltrpred` which
stores all output annotation and sequence files.


In detail, the following files and folders are generated by the `LTRpred()` function:

- __Folder `*_ltrpred`__
    - `*_ORF_prediction_nt.fsa` : Stores the predicted open reading frames within the predicted LTR transposons as DNA sequence.

    - `*_ORF_prediction_aa.fsa` : Stores the predicted open reading frames within the predicted LTR transposons as protein sequence.

    - `*_LTRpred.gff` : Stores the LTRpred predicted LTR transposons in GFF format.

    - `*_LTRpred.bed` : Stores the LTRpred predicted LTR transposons in BED format.

    - `*_LTRpred_DataSheet.tsv` : Stores the output table as data sheet.
    - __Folder `*_ltrharvest`__
        - `*_ltrharvest/*_BetweenLTRSeqs.fsa` : DNA sequences of the region between the LTRs in fasta format.

        - `*_ltrharvest/*_Details.tsv` : A spread sheet containing detailed information about the predicted LTRs.

        - `*_ltrharvest/*_FullLTRRetrotransposonSeqs.fsa` : DNA sequences of the entire predicted LTR retrotransposon.

        - `*_ltrharvest/*_index.fsa` : The suffixarray index file used to predict putative LTR retrotransposons.

        - `*_ltrharvest/*_Prediction.gff` : A spread sheet containing detailed additional information about the predicted LTRs (partially redundant with the *_Details.tsv file).
        
    - __Folder `*_ltrdigest`__
        - `*_ltrdigest/*_LTRdigestPrediction.gff` : A spread sheet containing detailed information about the predicted LTRs.
        
        - `*_ltrdigest/*_index_ltrdigest.fsa` : The suffixarray index file used to predict putative LTR retrotransposons with LTRdigest.

        - `*_ltrdigest/*-ltrdigest_tabout.csv` : A spread sheet containing additional detailed information about the predicted LTRs.

        - `*_ltrdigest/*-ltrdigest_complete.fas` : The full length DNA sequences of all predicted LTR transposons.

        - `*_ltrdigest/*-ltrdigest_conditions.csv` : Contains information about the parameters used for a given LTRdigest run.

        - `*_ltrdigest/*-ltrdigest_pbs.fas` : Stores the predicted PBS sequences for the putative LTR retrotransposons.

        - `*_ltrdigest/*-ltrdigest_ppt.fas` : Stores the predicted PPT sequences for the putative LTR retrotransposons.

        - `*_ltrdigest/*-ltrdigest_5ltr.fas` and `*-ltrdigest_3ltr.fas`: Stores the predicted 5' and 3' LTR sequences. Note: If the direction of the putative retrotransposon could be predicted, these files will contain the corresponding 3' and 5' LTR sequences. If no direction could be predicted, forward direction with regard to the original sequence will be assumed by LTRdigest, i.e. the 'left' LTR will be considered the 5' LTR.

        - `*_ltrdigest/*-ltrdigest_pdom_<domainname>.fas` : Stores the DNA sequences of the HMM matches to the LTR retrotransposon candidates.

        - `*_ltrdigest/*-ltrdigest_pdom_<domainname>_aa.fas` : Stores the concatenated protein sequences of the HMM matches to the LTR retrotransposon candidates.

        - `*_ltrdigest/*-ltrdigest_pdom_<domainname>_ali.fas` : Stores the alignment information for all matches of the given protein domain model to the translations of all candidates.

### Visualising functional retrotransposons annotated with ``LTRpred`` 

Finally, users can visualise the positioning of _de novo_ annotated retrotransposons along the chromosomes. Here, we choose an example based on the yeast genome.

```r
# install.packages("biomartr")
# download the yeast genome from ENSEMBL
sc_genome_path <- biomartr::getGenome(db = "ensembl", 
                               organism = "Saccharomyces cerevisiae", 
                               path = "yeast_genome", 
                               gunzip = TRUE)
# run LTRpred on the yeast genome (-> this may take a few minutes)
LTRpred::LTRpred(
    genome.file = sc_genome_path,
    cores = 4
)
# import functional retrotransposon annotation
sc_ltrpred <- LTRpred::read.ltrpred(
 file.path("Saccharomyces_cerevisiae_ltrpred",
                      "Saccharomyces_cerevisiae_LTRpred_DataSheet.tsv"))
# filter for potentially active candidates
sc_ltrpred_filtered <- LTRpred::quality.filter(sc_ltrpred, 
                                          sim = 95, 
                                          strategy = "stringent", 
                                          n.orfs = 1)
# visualise the position of potentially active retrotransposons
# along the yeast chromosomes
LTRpred::plot_element_distr_along_chromosome(sc_ltrpred_filtered,
                                             sc_genome_path)
```


![Positions of functional retrotransposons annotated by LTRpred along the yeast chromosomes.](Sc_LTRpred.png){ width=70% }

### Metagenome scale annotations

`LTRpred` allows users to generate annotations not only for single genomes
but for multiple genomes (metagenomes) using only one pipeline function named `LTRpred.meta()`. 

Users can download the `biomartr` package [@Drost2017-cw] to automatically retrieve genome assembly files for the species of interest.

```r
# specify the scientific names of the species of interest
# that shall first be downloaded and then be used
# to generate LTRpred annotations
species <- c("Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella") 
# download the genome assembly files for the species of interest
# 
# install.packages("biomartr")
biomartr::getGenomeSet(db = "refseq", organisms = species,
                      path = "store_genome_set")
# run LTRpred.meta() on the 3 species with 3 cores
LTRpred::LTRpred.meta(genome.folder = "store_genome_set",
                      output.folder = "LTRpred_meta_results",
                      cores = 3)
```        



        
# Acknowledgements

This work was supported by a European Research Council grant named EVOBREED [grant number 322621] and a Gatsby Fellowship [grant number AT3273/GLE]. The author also thanks Jerzy Paszkowski for supporting the development and application of ``LTRpred`` and for providing valuable suggestions to improve the manuscript.

# References## __LTRpred(ict)__: _de novo_ annotation of young and intact retrotransposons <img src="inst/LTRpred_logo.png" align="right" height="174" width="150" />

[![status](https://joss.theoj.org/papers/eeb2359d2459d3ae448cafac3ae33211/status.svg)](https://joss.theoj.org/papers/eeb2359d2459d3ae448cafac3ae33211)

Transposable elements (TEs) comprise vast parts of eukaryotic genomes.
In the past, TEs were seen as selfish mobile elements capable of populating a host genome to increase their chances for survival. By doing so they leave traces of junk DNA in host genomes that are usually regarded as by-products when sequencing, assembling, and annotating new genomes.

However, this picture is slowly changing ([Drost & Sanchez, 2019](https://academic.oup.com/gbe/article/11/12/3382/5637757)) and TEs have been shown to be involved in generating a diverse range of novel phenotypes.

Today, the _de novo_ detection of transposable elements is performed by annotation tools which try to detect any type of repeated sequence, TE family, or remnand DNA loci that can be associated with a known transposable element within a genome assembly. The main goal of such efforts is to retrieve a maximum amount of loci that can be associated with TEs. If successful, such annotation can then be used to mask host genomes and to perform classic (phylo-)genomics studies focusing on host genes.

More than [600 repeat and TE annotation tools](https://docs.google.com/spreadsheets/d/1UBK70zExiL0gFVaIAILiGhflCGXAq_SF_lymaxTE1pY/edit#gid=0) have been developed so far. Most of them are designed and optimized to annotate either the entire repeat space or specific superfamilies of TEs and their DNA remnants.

>The LTRpred pipeline has a different goal than all other annotation tools. It focuses particularly on [LTR retrotransposons](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC463057/) and aims to annotate only functional and potentially mobile elements. Such type of annotation is crucial for studying retrotransposon activity in eukaryotic genomes and to understand whether specific retrotransposon families can be activated artificially and harnessed to [mutagenize genomes at much faster speed](https://www.slcu.cam.ac.uk/news/tomato-jumping-genes).

In detail, `LTRpred` will take any genome assembly file in `fasta` format as input and will generate a detailed annotation of functional and potentially mobile LTR retrotransposons.

![](vignettes/LTRfeatures.png)


__Users can consult a comprehensive [Introduction](https://hajkd.github.io/LTRpred/articles/Introduction.html) to the `LTRpred` pipeline to get familiar with the tool.__

## Install

The fastest way to install `LTRpred` is via a [Docker container](https://hub.docker.com/r/drostlab/ltrpred).
Please make sure to read the [detailed installation instructions](https://hajkd.github.io/LTRpred/articles/Introduction.html#installation) to be able to
pass data to the container.

```bash
# retrieve docker image from dockerhub
docker pull drostlab/ltrpred
# run ltrpred container
docker run --rm -ti drostlab/ltrpred
# start R prompt within ltrpred container
~:/app# R
```

__Users who wish to run the `LTRpred` Docker container in a [conda](https://docs.conda.io/en/latest/) environment 
can use the [following approach based on UDocker](https://github.com/HajkD/LTRpred/issues/16) (Many thanks to Ilja Bezrukov).__ 

### Accessing LTRpred Container via RStudio
A more interactive way of performing analyses with `LTRpred` is via the [`RStudio`
version of LTRpred](https://hub.docker.com/r/drostlab/ltrpred_rstudio).
In [this LTRpred Docker Container](https://hub.docker.com/r/drostlab/ltrpred_rstudio)
users can access LTRpred within the container via `Rstudio`.

```bash
# retrieve docker image from dockerhub
docker pull drostlab/ltrpred_rstudio
# run ltrpred container
docker run -e PASSWORD=ltrpred --rm -p 8787:8787 -ti drostlab/ltrpred_rstudio
```

To open `RStudio` and interact with the container go to your standard web browser and type in the following URL:

```
http://localhost:8787

Username: rstudio

Password: ltrpred
```

Users can choose a custom password if they wish.

Within RStudio you can now run the example:

```r
LTRpred::LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))
```

Users can exit the container by pressing `Ctrl + c` multiple times.

Please find [all details here about how to use the Rstudio version here](https://hub.docker.com/r/drostlab/ltrpred_rstudio). 

## Citation
Please cite the following paper when using `LTRpred` for your own research:

> HG Drost. [__LTRpred: _de novo_ annotation of intact retrotransposons__](https://joss.theoj.org/papers/10.21105/joss.02170). __Journal of Open Source Software__, 5(50), 2170 (2020).

## Tutorials

### Quick Start

The fastest way to generate a LTR retrotransposon prediction for a genome of interest (after [installing](https://hajkd.github.io/LTRpred/articles/Introduction.html) all prerequisite command line tools) is to use the
`LTRpred()` function and relying on the default parameters. In the following example,
a LTR transposon prediction is performed for parts of the Human Y chromosome.

```r
# Perform de novo LTR transposon prediction for the Human Y chromosome
LTRpred::LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))
```

When running your own genome, please specify `genome.file = "path/to/your/genome.fasta` instead of `system.file(..., package = "LTRpred")`. The command `system.file(..., package = "LTRpred")` merely references the path to the example file stored in the LTRpred package itself.


This tutorial introduces users to `LTRpred`:

- [Introduction to LTRpred](https://hajkd.github.io/LTRpred/articles/Introduction.html)

Users can also read the tutorials within ([RStudio](http://www.rstudio.com/)) :

```r
library(LTRpred)
browseVignettes("LTRpred")
```

You can also find a list of all available `LTRpred` functions here: https://hajkd.github.io/LTRpred/reference/index.html

### Studies that successfully used `LTRpred` to annotate functional retrotransposons

> - Z Wang & D Baulcombe. [__Transposon age and non-CG methylation__](https://www.nature.com/articles/s41467-020-14995-6). __Nature Communications__, 11, 1221 (2020).
> - JH Collins, KW Keating, TR Jones, S Balaji, CB Marsan et al. [__Engineered yeast genomes accurately assembled from pure and mixed samples__](https://www.nature.com/articles/s41467-021-21656-9). __Nature communications__, 12, 1485 (2021).
>
> - H Kundariya et al. [__MSH1-induced heritable enhanced growth vigor through grafting is associated with the RdDM pathway in plants__](https://www.nature.com/articles/s41467-020-19140-x) __Nature Communications__, 11, 5343 (2020).
>
> - J Cho, M Benoit, M Catoni, __HG Drost__, A Brestovitsky, M Oosterbeek and J Paszkowski.  [__Sensitive detection of pre-integration intermediates of LTR retrotransposons in crop plants__](https://www.nature.com/articles/s41477-018-0320-9). __Nature Plants__, 5,  26-33 (2019).
>
> - M Benoit, __HG Drost__, M Catoni, Q Gouil, S Lopez-Gomollon, DC Baulcombe, J Paszkowski. [__Environmental and epigenetic regulation of Rider retrotransposons in tomato__](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008370). __PloS Genetics__, 15(9): e1008370 (2019). 
>
> - Nguinkal _et al._ [__The First Highly Contiguous Genome Assembly of Pikeperch (Sander lucioperca), an Emerging Aquaculture Species in Europe__](https://www.mdpi.com/2073-4425/10/9/708/htm) __Genes__, 0(9), 708 (2019).
> - E Cerruti, C Gisbert, __HG Drost__, D Valentino, E Portis, L Barchi, J Prohens, S Lanteri, C Comino,  M Catoni. [__Epigenetic bases of grafting-induced vigour in eggplant__](https://www.biorxiv.org/content/10.1101/831719v1). __bioaRxiv__ (2019).
>
> - P Gan, R Hiroyama, A Tsushima, S Masuda _et al_.
[__Subtelomeric regions and a repeat-rich chromosome harbor multicopy effector gene clusters with variable conservation in multiple plant pathogenic Colletotrichum species__](https://www.biorxiv.org/content/10.1101/2020.04.28.061093v1.abstract) __bioRxiv__ (2020)
>
> - J Wang et al. [__Gigantic Genomes Can Provide Empirical Tests of TE Dynamics Models--An Example from Amphibians__](https://www.biorxiv.org/content/10.1101/2020.08.19.257527v1). __bioRxiv__, (2020).
>
> - Y Ayukawa et al. [__A pair of effectors encoded on a conditionally dispensable chromosome of Fusarium oxysporum suppress host-specific immunity__](https://www.biorxiv.org/content/10.1101/2020.10.06.329052v1). __bioRxiv__, (2020).
> 
> - C Meguerditchian, A Ergun, V Decroocq, M Lefebvre et al. [__Pipeline to detect the relationship between transposable elements and adjacent genes in host genome__](https://www.biorxiv.org/content/10.1101/2021.02.25.432867v1.abstract) __bioRxiv__, (2021).


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/LTRpred/issues


---
title: "Introduction to LTRpred"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

# Table of Contents
- [1. Introduction](#introduction)
- [2. Getting Started](#getting-started)
    - [2.1 Installation](#installation)
        - [2.1.1 LTRpred Docker Container](#ltrpred-docker-container)
        - [2.1.2 Install tool dependencies on Linux](#install-tool-dependencies-on-linux)
    - [2.2 Quick Start](#quick-start)
    - [2.3 LTRpred output](#ltrpred-output)
    - [2.4 Import LTRpred output](#import-ltrpred-output)
    - [2.5 Output file format of LTRpred()](#output-file-format-of-ltrpred)
- [3. Detailed LTRpred run](#detailed-ltrpred-run)
    - [3.1 HMM Models](#hmm-models)
    - [3.2 Detailed description of adjustable LTRpred parameters](#detailed-description-of-adjustable-ltrpred-parameters)
- [4. Metagenome scale annotations](#metagenome-scale-annotations)

# Introduction

The `LTRpred` package implements a software pipeline and provides an integrated workflow to screen for intact and potentially active LTR retrotransposons in any genomic sequence of interest. For this purpose, this package provides a rich set of analytics tools to allow researchers to quickly start annotating and explore their own genomes.

The _de novo_ prediction of LTR transposons in `LTRpred` is based on the command line tools 
[LTRharvest](http://www.zbh.uni-hamburg.de/?id=206) and [LTRdigest](http://www.zbh.uni-hamburg.de/forschung/arbeitsgruppe-genominformatik/software/ltrdigest.html) and extends these search strategies by additional analytics
modules to filter the search space of putative LTR transposons for biologically meaningful candidates.

Please make sure that all command line tools are installed properly before running any `LTRpred` function.

# Getting Started

The rationale for implementing `LTRpred` was to implement an R based pipeline combining the most
sensitive, accurate, and conservative state-of-the art LTR retrotransposon
detection tools and extend their inference by additional analyses and quality filtering steps to
screen for functional and structurally intact elements. Hence, _LTRpred_
aims to provide a high-level _de novo_ prediction infrastructure to generate high quality annotations
of intact LTR retrotransposons. All `LTRpred` functions are generic so that parameters can be changed to detect any form of LTR retrotransposon in any genome.

Internally, _LTRpred_ is based on the `de novo` annotation tools `LTRharvest` and `LTRdigest` which use prior knowledge
about DNA sequence features (also referred to as `Structure-based methods`) such as the homology of Long Terminal Repeats (LTRs), Primer Binding Sites (PBS), _gag_ and _pol_ protein domains, and
Target Site Duplications (TSDs) that are known to enable the process of transposition ([Bergman and Quesneville, 2007](http://bib.oxfordjournals.org/content/8/6/382.long)) to infer LTR retrotransposons in any genome.

Hence, these `de novo` annotation tools are designed to screen the genome systematically and efficiently for
these structural DNA features. __Figure 1__ shows the structural features of LTR retrotransposons
that are used for predicting putative LTR transposons `de novo` in any genome of interest.

![Sequence features of LTR retrotransposons. The structural element (LTR length, LTR similarity, TSD, max. size of full LTR retrotransposon, PBS length, tRNA binding, protein domain search, etc.) can be modified and controlled separately by corresponding arguments implemented in the `LTRpred()` function. In addition to controlling the structural features of candidate LTR retrotransposons, a open reading frame prediction, Dfam annotation, copy number clustering, and LTR copy number estimation, and methylation context quantification is performed subsequently after reporting putative LTR retrotransposons.](LTRfeatures.png)

Based on the optimized output of these tools, the `LTRpred` package aims to provide researchers with maximum flexibility of
adjustable parameters to detect any type of functional LTR retrotransposon. `LTRpred` package allows users to modify a vast range of
parameters to screen for potential LTR transposons, so having this template shown in
__Figure 1__ in mind will help researchers to modify structural parameters in a biologically meaningful manner.

## Installation

The fastest way to run `LTRpred` is to download the `ltrpred` Docker container which
includes all pre-installed tool dependencies. For users who cannot use a Docker environment
the individual installation instructions for each dependency tool are listed below (only for Linux).

### LTRpred Docker Container

Please be aware that the `drostlab/ltrpred` container (command line version) is suitable for `Linux`, `macOS`, and `Windows` users. Whereas the `RStudio` container version `drostlab/ltrpred_rstudio` is not suitable for `Windows` users, because a port bridge cannot be established.  

Please make sure to create a [Docker](https://www.docker.com/get-started) account and to [install Docker](https://docs.docker.com/engine/install/) on your system.

#### Download `drostlab/ltrpred` container for use with R command line


```bash
# retrieve docker image from dockerhub
docker pull drostlab/ltrpred
# run ltrpred container
docker run --rm -ti drostlab/ltrpred
# start R prompt within ltrpred container
~:/app# R
```

Users who wish to run the `LTRpred` Docker container in a [conda](https://docs.conda.io/en/latest/) environment 
can use the [following approach based on UDocker](https://github.com/HajkD/LTRpred/issues/16) (Many thanks to Ilja Bezrukov).

Within the ltrpred container R prompt run the `ltrpred` example:

```r
LTRpred::LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))
```

To exit R in the container run:

```r
q()
```

And to exit the `ltrpred` container run:

```bash
~:/app# exit
```

Now, users can add their own genome data as well as the Dfam database for
further annotation to the `ltrpred` container by following these steps (in a different Terminal window):

```bash
# go to the folder path in which you want to
# store all genome and Dfam data you want to
# mount in the ltrpred container and then run:
# create a new folder which will store
# all files that will be required in the
# ltrpred container
mkdir ltrpred_data
cd ltrpred_data
# create a dfam database folder
mkdir Dfam
cd Dfam
```

Now users can download and format the Dfam database as follows (within the `Dfam` folder created above). Unfortunately, the `Dfam` database size is too large to make it part of the `drostlab/ltrpred` container. In addition, the database is frequently curated and updated. Thus, it is recommended that users download and format the `Dfam` database to their local hard drive and mount it to the running `drostlab/ltrpred` container.

__To format the Dfam database locally, users need to install [HMMER](http://hmmer.org/download.html) on their local machine (to use `hmmpress`). However, within the `drostlab/ltrpred` and `drostlab/ltrpred_rstudio` containers HMMER is already preinstalled and does not need to be installed by the user.__ An example installation of HMMER for Linux machines is listed below.

For `macOS` users, please install `wget` on your `masOS` machine using [Homebrew](http://brew.sh/).

```bash
wget https://www.dfam.org/releases/Dfam_3.1/families/Dfam.hmm.gz
gunzip Dfam.hmm.gz
# format database by running hmmpress
hmmpress Dfam.hmm
cd ..
```

Next, make sure to also store the genome assembly file (in `fasta` format) you
want to _de novo_ annotate with `LTRpred` in the `ltrpred_data` folder you just created.

A possible way to retrieve such a genome is (within R) using [biomartr](https://github.com/ropensci/biomartr).
If you use `biomartr` please make sure to [install all `biomartr` package dependencies](https://github.com/ropensci/biomartr#installation) before running the following code.

```r
# install.packages("biomartr")
biomartr::getGenome(db = "ensembl", 
                    organism = "Saccharomyces cerevisiae", 
                    path = "yeast_genome", 
                    gunzip = TRUE)
```

The respective genome assembly file is now stored at `yeast_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa` and needs to be copied
into the `ltrpred_data` folder you just created.


Now users can mount the `ltrpred_data` folder to the `ltrpred` Docker container (using the `-v` option). This `-v` mounting option is also available for the `RStudio` container version and can also be run within the `RStudio` Terminal. 

```bash
docker run --rm -p 8787:8787 -v /put/here/your/path/to/ltrpred_data:/app/ltrpred_data -ti drostlab/ltrpred
```

Within the `ltrpred` Docker container the `ltrpred_data` folder is now stored in the
working directory `/app`.

When running the `ltrpred` Docker container you should be able to see the
mounted `ltrpred_data` folder as following:

```bash
~:/app# ls
```

```
latticeExtra_0.6-28.tar.gz
ltrpred_data
software_downloads
```

Next, users can again run the `R prompt` within the `ltrpred` Docker container to
run `LTRpred` with the local data that was mounted:

```bash
~:/app# R
```

```r
# run LTRpred on the yeast genome that was mounted
# to the ltrpred container from your local 'ltrpred_data' folder
LTRpred(genome.file = "ltrpred_data/yeast_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa", cores = 2)
```

As you can see, within the `ltrpred` container `R prompt` the current working directory is `/app`. 

To also include the `Dfam` database for further annotation users can specify
the path to the Dfam folder:


```r
# run LTRpred on the yeast genome that was mounted
# to the ltrpred container from your local 'ltrpred_data' folder
LTRpred(genome.file = "ltrpred_data/yeast_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa", 
    annotate = "Dfam",
		Dfam.db = "ltrpred_data/Dfam", 
		cores = 2)
```


Please be aware that using the `Dfam` database for further annotation significantly increases
the computation time of the `LTRpred` pipeline.

#### Retrieve LTRpred output files from Docker container

Next, users can retrieve the `LTRpred` generated results from the docker container by opening a second `Terminal` window (while the `drostlab/ltrpred` container remains running) and perform the following steps:

1) Close `R` in the running docker container using `q()`.
2) This should bring you back to the docker bash `root@ac389gfja8089:/app#`.
3) Copy the docker ID, in this case `ac389gfja8089` (this docker ID is just an example, please use the docker ID shown on your system).
4) Type in the newly opened `Terminal` window:

```bash
# copy Hsapiens_ChrY_ltrpred output from docker container to hard drive
docker cp ac389gfja8089:/app/Hsapiens_ChrY_ltrpred path/to/your/host/hard/drive/Hsapiens_ChrY_ltrpred
```

This example assumes that you ran the example `LTRpred` run `LTRpred::LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))` which created the output folder `Hsapiens_ChrY_ltrpred` in the docker folder `/app/`. 

Please note, that if you specify different file paths when creating files within the docker container,
these must be adjusted when running:

```bash
# copy files from docker container to hard drive
docker cp ac389gfja8089:/app/your/inside/docker/path/here/Hsapiens_ChrY_ltrpred path/to/your/host/hard/drive/Hsapiens_ChrY_ltrpred
```


#### Download `drostlab/ltrpred_rstudio` container for use with RStudio Server


```bash
# retrieve docker image from dockerhub
docker pull drostlab/ltrpred_rstudio
# run ltrpred container
docker run -e PASSWORD=ltrpred --rm -p 8787:8787 -ti drostlab/ltrpred_rstudio
```

To open `RStudio` and interact with the container go to your standard web browser and type in the following url:

`http://localhost:8787`

Username: `rstudio`
Password: `ltrpred`


Within `RStudio` you can now run the example:

```r
LTRpred::LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))
```


Users can `exit` the container by pressing `Ctrl + c` multiple times. 

Next, users can mount their `ltrpred_data` folder to the RStudio server run
the same way they mounted folders in the command line container version (using `-v`).
This folder mounting can also be run within the `RStudio` Terminal of the `drostlab/ltrpred_rstudio` container.


```bash
# retrieve docker image from dockerhub
docker pull drostlab/ltrpred_rstudio
# run ltrpred container
docker run -e PASSWORD=ltrpred --rm -p 8787:8787 -v /put/here/your/path/to/ltrpred_data:/home/rstudio/ltrpred_data -ti drostlab/ltrpred_rstudio
```

Now go to your standard web browser and type in the following url:

`http://localhost:8787`

Username: `rstudio`
Password: `ltrpred`

In RStudio type:

```r
list.files()
```

You should be able to see the `ltrpred_data` folder.

Users can `exit` the container by pressing `Ctrl + c` multiple times. 

#### Retrieve LTRpred output files from Docker container

Next, users can retrieve the `LTRpred` generated results from the docker container by opening a second `Terminal` window (while the `drostlab/ltrpred` container remains running) and perform the following steps:

1) Close `R` in the running docker container using `q()`.
2) This should bring you back to the docker bash `root@ac389gfja8089:/app#`.
3) Copy the docker ID, in this case `ac389gfja8089` (this docker ID is just an example, please use the docker ID shown on your system).
4) Type in the newly opened `Terminal` window:

```bash
# copy Hsapiens_ChrY_ltrpred output from docker container to hard drive
docker cp ac389gfja8089:/app/Hsapiens_ChrY_ltrpred path/to/your/host/hard/drive/Hsapiens_ChrY_ltrpred
```

This example assumes that you ran the example `LTRpred` run `LTRpred::LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))` which created the output folder `Hsapiens_ChrY_ltrpred` in the docker folder `/app/`. 

Please note, that if you specify different file paths when creating files within the docker container,
these must be adjusted when running:

```bash
# copy files from docker container to hard drive
docker cp ac389gfja8089:/app/your/inside/docker/path/here/Hsapiens_ChrY_ltrpred path/to/your/host/hard/drive/Hsapiens_ChrY_ltrpred
```

__Please read more details about how to transfer genome files and the Dfam database in the
previous section [Download ltrpred container for use with R command line](#ltrpred-docker-container).

### Install Tool Dependencies on Linux

#### Programming Languages

Please make sure that the following programming languages are installed on your system:

- [Perl](https://www.perl.org/get.html)
- [Ruby](https://www.ruby-lang.org/en/documentation/installation/)
- [Python](https://www.python.org/downloads/)
- [C/C++](https://www.cplusplus.com/)
- [R](http://lib.stat.cmu.edu/R/CRAN/)

#### Install Programming languages and Linux Tools

```bash
apt-get update \
&& apt-get -y install apt-utils \
&& apt-get -y install  gcc \
&& apt-get -y install python3 \
&& apt-get -y install perl \
&& apt-get -y install make \
&& apt-get -y install sudo \
&& apt-get -y install wget \
&& apt-get -y install genometools \
&& apt-get -y install git \
&& apt-get -y install autoconf \
&& apt-get -y install g++ \
&& apt-get -y install ncbi-blast+ \
&& apt-get -y install build-essential \
&& apt-get -y install libcurl4-gnutls-dev \
&& apt-get -y install libxml2-dev \
&& apt-get -y install libssl-dev \
&& apt-get -y install libpq-dev \
&& apt-get -y install software-properties-common
```

#### Install [HMMER](http://hmmer.org/download.html)

A detailed description of how to install `HMMER` for several operating systems 
can be found [here](http://hmmer.org/download.html).

```bash
mkdir software_downloads \
  && cd software_downloads \
  && wget http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz \
  && tar xf hmmer-3.3.tar.gz \
  && cd hmmer-3.3 \
  && ./configure \
  && make \
  && make check \
  && sudo make install \
  && cd ..
```

#### Install [USEARCH](http://drive5.com/usearch/download.html)

A detailed description of how to install `USEARCH` for several operating systems 
can be found [here](http://drive5.com/usearch/manual/install.html).

First, users will need to register and download USEARCH for their operating system
from http://drive5.com/usearch/download.html .

After downloading USEARCH you will need to install it as a command line tool in your `/usr/local/` directory and you should be able to execute the following command in your Terminal:

```sh
cd software_downloads \
  && wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz \
  && gunzip usearch11.0.667_i86linux32.gz \
  && chmod +x usearch11.0.667_i86linux32 \
  && sudo mv usearch11.0.667_i86linux32 usearch \
  && sudo cp usearch /usr/local/bin/usearch \
  && cd ..
usearch -version
```

```
usearch v11.0.667_i86linux32
```

#### Install [VSEARCH](https://github.com/torognes/vsearch)

A detailed description of how to install `VSEARCH` for several operating systems 
can be found [here](https://github.com/torognes/vsearch).

Please [install git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) before running the following commands.

```bash
cd software_downloads \
  && wget https://github.com/torognes/vsearch/archive/v2.14.2.tar.gz \
  && tar xzf v2.14.2.tar.gz \
  && cd vsearch-2.14.2 \
  && sudo ./autogen.sh \
  && sudo ./configure \
  && sudo make \
  && sudo make install \
  && cd ..
```

### Install [dfamscan.pl](http://www.dfam.org/web_download/Current_Release/dfamscan.pl)

[dfamscan.pl](https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan.pl.gz) needs to be unzipped and stored at `/usr/local/bin/dfamscan.pl` and executabe via `perl /usr/local/bin/dfamscan.pl -help`. This is important to be able to run the hmmer search against the Dfam database.


```bash
cd software_downloads \
  && wget https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan.pl.gz \
  && gunzip dfamscan.pl.gz \
  && sudo cp dfamscan.pl /usr/local/bin/dfamscan.pl \
  && cd ..
```


Users can download the [Dfam database v3.1](https://www.dfam.org/releases/Dfam_3.1/families/) by running:

```sh
wget https://www.dfam.org/releases/Dfam_3.1/families/Dfam.hmm.gz
gunzip Dfam.hmm.gz
# run hmmpress
hmmpress Dfam.hmm
```

The path to the folder where the formatted `Dfam` database can be found can then be passed
as argument `Dfam.db` to the `LTRpred::LTRpred()` function.

### Install R packages

Please make sure that [Bioconductor](https://www.bioconductor.org/install/) and all package dependencies are installed on the system on which you would like to run `LTRpred`.

Install prerequisite CRAN and Bioconductor packages:

```r
install.packages('devtools')
install.packages('tidyverse')
install.packages('BiocManager')
BiocManager::install()
BiocManager::install(c('rtracklayer', 'GenomicFeatures', 'GenomicRanges', 'GenomeInfoDb', 'biomaRt', 'Biostrings'))
install.packages(c('tidyverse', 'data.table', 'seqinr', 'biomartr', 'ape', 'dtplyr', 'devtools'))
devtools::install_github('HajkD/metablastr', build_vignettes = TRUE, dependencies = TRUE)
install.packages(c('BSDA', 'ggrepel', 'gridExtra'))
https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-28.tar.gz
install.packages('latticeExtra_0.6-28.tar.gz', type = 'source')
install.packages('survival')
BiocManager::install('ggbio')
```

Now users may install `LTRpred` as follows:

```r
# install.packages("devtools")
devtools::install_github("HajkD/LTRpred")
```


## Quick Start

The fastest way to generate a LTR retrotransposon prediction for a genome of interest (after [installing](https://hajkd.github.io/LTRpred/articles/Introduction.html#installation) all prerequisite command line tools) is to use the
`LTRpred()` function and relying on the default parameters. In the following example,
a LTR transposon prediction is performed for parts of the Human Y chromosome.

```r
# load LTRpred package
library(LTRpred)
# de novo LTR transposon prediction for the Human Y chromosome
LTRpred(
    genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"),
    cores = 4
)
```

```
Running LTRpred on genome '/Library/Frameworks/R.framework/Versions/3.6/Resources/library/LTRpred/Hsapiens_ChrY.fa' with 4 core(s) and searching for retrotransposons using the overlaps option (overlaps = 'no') ...


No hmm files were specified, thus the internal HMM library will be used! See '/Library/Frameworks/R.framework/Versions/3.6/Resources/library/LTRpred/HMMs/hmm_*' for details.
No tRNA files were specified, thus the internal tRNA library will be used! See '/Library/Frameworks/R.framework/Versions/3.6/Resources/library/LTRpred/tRNAs/tRNA_library.fa' for details.
The output folder 'Hsapiens_ChrY_ltrpred' seems to exist already and will be used to store LTRpred results ...


LTRpred - Step 1:
Run LTRharvest...
LTRharvest: Generating index file Hsapiens_ChrY_ltrharvest/Hsapiens_ChrY_index.fsa with gt suffixerator...
Running LTRharvest and writing results to Hsapiens_ChrY_ltrharvest...
LTRharvest analysis finished!


LTRpred - Step 2:
Generating index file Hsapiens_ChrY_ltrdigest/Hsapiens_ChrY_index_ltrdigest.fsa with suffixerator...
LTRdigest: Sort index file...
Running LTRdigest and writing results to Hsapiens_ChrY_ltrdigest...
LTRdigest analysis finished!


LTRpred - Step 3:
Import LTRdigest Predictions...

Input:  Hsapiens_ChrY_ltrdigest/Hsapiens_ChrY_LTRdigestPrediction.gff  -> Row Number:  179
Remove 'NA' -> New Row Number:  179
(1/8) Filtering for repeat regions has been finished.
(2/8) Filtering for LTR retrotransposons has been finished.
(3/8) Filtering for inverted repeats has been finished.
(4/8) Filtering for LTRs has been finished.
(5/8) Filtering for target site duplication has been finished.
(6/8) Filtering for primer binding site has been finished.
(7/8) Filtering for protein match has been finished.
(8/8) Filtering for RR tract has been finished.


LTRpred - Step 4:
Perform ORF Prediction using 'usearch -fastx_findorfs' ...
usearch v8.1.1861_i86osx32, 4.0Gb RAM (17.2Gb total), 8 cores
(C) Copyright 2013-15 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

00:00 1.9Mb  100.0% Working
Join ORF Prediction table: nrow(df) = 24 candidates.
unique(ID) = 24 candidates.
unique(orf.id) = 24 candidates.


LTRpred - Step 5:
Perform methylation context quantification..
Join methylation context (CG, CHG, CHH, CCG) count table: nrow(df) = 24 candidates.
unique(ID) = 24 candidates.
unique(orf.id) = 24 candidates.
Copy files to result folder 'Hsapiens_ChrY_ltrpred'.


LTRpred - Step 6:
Starting retrotransposon evolutionary age estimation by comparing the 3' and 5' LTRs using the molecular evolution model 'K80' and the mutation rate '1.3e-07' (please make sure the mutation rate can be assumed for your species of interest!) for 24 predicted elements ...


Please be aware that evolutionary age estimation based on 3' and 5' LTR comparisons are only very rough time estimates and don't take reverse-transcription mediated retrotransposon recombination between family members of retroelements into account! Please consult Sanchez et al., 2017 Nature Communications and Drost & Sanchez, 2019 Genome Biology and Evolution for more details on retrotransposon recombination.


LTRpred - Step 7:
The LTRpred prediction table has been filtered (default) to remove potential false positives. Predicted LTRs must have an PBS or Protein Domain and must fulfill thresholds: sim = 70%; #orfs = 0. Furthermore, TEs having more than 10% of N's in their sequence have also been removed.
Input #TEs: 24
Output #TEs: 21

LTRpred finished all analyses successfully. All output files were stored at 'Hsapiens_ChrY_ltrpred'.
[1] "Successful job 1 ."
```

### LTRpred output

The `LTRpred()` function internally generates a folder named `*_ltrpred` which
stores all output annotation and sequence files.


In detail, the following files and folders are generated by the `LTRpred()` function:

- __Folder `*_ltrpred`__
    - `*_ORF_prediction_nt.fsa` : Stores the predicted open reading frames within the predicted LTR transposons as DNA sequence.

    - `*_ORF_prediction_aa.fsa` : Stores the predicted open reading frames within the predicted LTR transposons as protein sequence.

    - `*_LTRpred.gff` : Stores the LTRpred predicted LTR transposons in GFF format.

    - `*_LTRpred.bed` : Stores the LTRpred predicted LTR transposons in BED format.

    - `*_LTRpred_DataSheet.tsv` : Stores the output table as data sheet.
    - __Folder `*_ltrharvest`__
        - `*_ltrharvest/*_BetweenLTRSeqs.fsa` : DNA sequences of the region between the LTRs in fasta format.

        - `*_ltrharvest/*_Details.tsv` : A spread sheet containing detailed information about the predicted LTRs.

        - `*_ltrharvest/*_FullLTRRetrotransposonSeqs.fsa` : DNA sequences of the entire predicted LTR retrotransposon.

        - `*_ltrharvest/*_index.fsa` : The suffixarray index file used to predict putative LTR retrotransposons.

        - `*_ltrharvest/*_Prediction.gff` : A spread sheet containing detailed additional information about the predicted LTRs (partially redundant with the *_Details.tsv file).

        - `*_ltrdigest/*_index_ltrdigest.fsa` : The suffixarray index file used to predict putative LTR retrotransposonswith LTRdigest.

    - __Folder `*_ltrdigest`__
        - `*_ltrdigest/*_LTRdigestPrediction.gff` : A spread sheet containing detailed information about the predicted LTRs.

        - `*_ltrdigest/*-ltrdigest_tabout.csv` : A spread sheet containing additional detailed information about the predicted LTRs.

        - `*_ltrdigest/*-ltrdigest_complete.fas` : The full length DNA sequences of all predicted LTR transposons.

        - `*_ltrdigest/*-ltrdigest_conditions.csv` : Contains information about the parameters used for a given LTRdigest run.

        - `*_ltrdigest/*-ltrdigest_pbs.fas` : Stores the predicted PBS sequences for the putative LTR retrotransposons.

        - `*_ltrdigest/*-ltrdigest_ppt.fas` : Stores the predicted PPT sequences for the putative LTR retrotransposons.

        - `*_ltrdigest/*-ltrdigest_5ltr.fas` and `*-ltrdigest_3ltr.fas`: Stores the predicted 5' and 3' LTR sequences. Note: If the direction of the putative retrotransposon could be predicted, these files will contain the corresponding 3' and 5' LTR sequences. If no direction could be predicted, forward direction with regard to the original sequence will be assumed by LTRdigest, i.e. the 'left' LTR will be considered the 5' LTR.

        - `*_ltrdigest/*-ltrdigest_pdom_<domainname>.fas` : Stores the DNA sequences of the HMM matches to the LTR retrotransposon candidates.

        - `*_ltrdigest/*-ltrdigest_pdom_<domainname>_aa.fas` : Stores the concatenated protein sequences of the HMM matches to the LTR retrotransposon candidates.

        - `*_ltrdigest/*-ltrdigest_pdom_<domainname>_ali.fas` : Stores the alignment information for all matches of the given protein domain model to the translations of all candidates.

### Import LTRpred output 

The `LTRpred()` output table `*_LTRpred_DataSheet.tsv` is in [tidy](http://vita.had.co.nz/papers/tidy-data.html) format and can then be imported using `read.ltrpred()`.
The `tidy` output format is designed to work seamlessly with the [tidyverse](https://www.tidyverse.org/) and [R data science](http://r4ds.had.co.nz/) framework.

```r
# import LTRpred prediction output
Hsapiens_chrY <- read.ltrpred("Hsapiens_ChrY_ltrpred/Hsapiens_ChrY_LTRpred_DataSheet.tsv")
# look at some results
dplyr::select(Hsapiens_chrY, ltr_similarity:end, tRNA_motif, Clust_cn)
```

```
# A tibble: 14 x 9
   ltr_similarity similarity protein_domain  orfs        chromosome
            <dbl>      <chr>          <chr> <int>             <chr>
 1          80.73    (80,82]          RVT_1     1 NC000024.10Homosa
 2          89.85    (88,90]          RVT_1     1 NC000024.10Homosa
 3          79.71    (78,80]           <NA>     0 NC000024.10Homosa
 4          84.93    (84,86]          RVT_1     0 NC000024.10Homosa
 5          75.52       <NA>          RVT_1     0 NC000024.10Homosa
 6          76.47    [76,78]          RVT_1     1 NC000024.10Homosa
 7          80.28    (80,82]           <NA>     0 NC000024.10Homosa
 8          76.92    [76,78]          RVT_1     0 NC000024.10Homosa
 9          89.53    (88,90]          RVT_1     0 NC000024.10Homosa
10          93.57    (92,94]          RVT_1     1 NC000024.10Homosa
11          82.35    (82,84]          RVT_1     2 NC000024.10Homosa
12          79.51    (78,80]          RVT_1     0 NC000024.10Homosa
13          78.42    (78,80]          RVT_1     3 NC000024.10Homosa
14          92.75    (92,94]          RVT_1     1 NC000024.10Homosa
# ... with 4 more variables: start <int>, end <int>, tRNA_motif <chr>,
#   Clust_cn <int>
```

Looking at all columns:

```r
dplyr::glimpse(Hsapiens_chrY)
```

```
Observations: 21
Variables: 92
$ species                 <chr> "Hsapiens_ChrY", "Hsapien...
$ ID                      <chr> "Hsapiens_ChrY_LTR_retrot...
$ dfam_target_name        <chr> NA, NA, NA, NA, NA, NA, N...
$ ltr_similarity          <dbl> 80.73, 89.85, 79.71, 83.2...
$ ltr_age_mya             <dbl> 0.7936246, 0.2831139, 0.7...
$ similarity              <chr> "(80,82]", "(88,90]", "(7...
$ protein_domain          <chr> "RVT_1", "RVT_1", NA, NA,...
$ orfs                    <int> 1, 1, 0, 0, 0, 0, 0, 1, 0...
$ chromosome              <chr> "NC000024.10Homosa", "NC0...
$ start                   <int> 3143582, 3275798, 3313536...
$ end                     <int> 3162877, 3299928, 3318551...
$ strand                  <chr> "-", "-", "+", "+", "-", ...
$ width                   <int> 19296, 24131, 5016, 12952...
$ annotation              <chr> "LTR_retrotransposon", "L...
$ pred_tool               <chr> "LTRpred", "LTRpred", "LT...
$ frame                   <chr> ".", ".", ".", ".", ".", ...
$ score                   <chr> ".", ".", ".", ".", ".", ...
$ lLTR_start              <int> 3143582, 3275798, 3313536...
$ lLTR_end                <int> 3143687, 3276408, 3313665...
$ lLTR_length             <int> 106, 611, 130, 126, 218, ...
$ rLTR_start              <int> 3162769, 3299338, 3318414...
$ rLTR_end                <int> 3162877, 3299928, 3318551...
$ rLTR_length             <int> 109, 591, 138, 137, 219, ...
$ lTSD_start              <int> 3143578, 3275794, 3313532...
$ lTSD_end                <int> 3143581, 3275797, 3313535...
$ lTSD_motif              <chr> "acag", "ttgt", "ttag", "...
$ rTSD_start              <int> 3162878, 3299929, 3318552...
$ rTSD_end                <int> 3162881, 3299932, 3318555...
$ rTSD_motif              <chr> "acag", "ttgt", "ttag", "...
$ PPT_start               <int> NA, NA, NA, NA, NA, 34660...
$ PPT_end                 <int> NA, NA, NA, NA, NA, 34660...
$ PPT_motif               <chr> NA, NA, NA, NA, NA, "agag...
$ PPT_strand              <chr> NA, NA, NA, NA, NA, "+", ...
$ PPT_offset              <int> NA, NA, NA, NA, NA, 23, N...
$ PBS_start               <int> NA, NA, 3313667, 3372512,...
$ PBS_end                 <int> NA, NA, 3313677, 3372522,...
$ PBS_strand              <chr> NA, NA, "+", "+", "-", "+...
$ tRNA                    <chr> NA, NA, "Homo_sapiens_tRN...
$ tRNA_motif              <chr> NA, NA, "aattagctgga", "c...
$ PBS_offset              <int> NA, NA, 1, 3, 0, 5, 2, 5,...
$ tRNA_offset             <int> NA, NA, 1, 0, 2, 5, 1, 5,...
$ `PBS/tRNA_edist`        <int> NA, NA, 1, 1, 1, 1, 1, 1,...
$ orf.id                  <chr> "NC000024.10Homosa_314358...
$ repeat_region_length    <int> 19304, 24139, 5024, 12960...
$ PPT_length              <int> NA, NA, NA, NA, NA, 27, N...
$ PBS_length              <int> NA, NA, 11, 11, 11, 11, 1...
$ dfam_acc                <chr> NA, NA, NA, NA, NA, NA, N...
$ dfam_bits               <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_e_value            <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_bias               <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_strand             <chr> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_modlen             <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_target_description <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Cluster           <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Target            <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Perc_Ident        <dbl> NA, NA, NA, NA, NA, NA, N...
$ Clust_cn                <int> NA, NA, NA, NA, NA, NA, N...
$ TE_CG_abs               <dbl> 62, 125, 35, 70, 139, 83,...
$ TE_CG_rel               <dbl> 0.003213101, 0.005180059,...
$ TE_CHG_abs              <dbl> 659, 830, 150, 396, 742, ...
$ TE_CHG_rel              <dbl> 0.03415216, 0.03439559, 0...
$ TE_CHH_abs              <dbl> 2571, 3454, 748, 1743, 31...
$ TE_CHH_rel              <dbl> 0.1332400, 0.1431354, 0.1...
$ TE_CCG_abs              <dbl> 13, 24, 6, 15, 33, 16, 4,...
$ TE_CCG_rel              <dbl> 0.0006737148, 0.000994571...
$ TE_N_abs                <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_3ltr_abs             <dbl> 1, 0, 1, 4, 8, 3, 1, 11, ...
$ CG_3ltr_rel             <dbl> 0.009433962, 0.000000000,...
$ CHG_3ltr_abs            <dbl> 2, 24, 9, 8, 14, 14, 9, 9...
$ CHG_3ltr_rel            <dbl> 0.01886792, 0.03927987, 0...
$ CHH_3ltr_abs            <dbl> 18, 69, 18, 26, 43, 70, 2...
$ CHH_3ltr_rel            <dbl> 0.16981132, 0.11292962, 0...
$ CCG_3ltr_abs            <dbl> 0, 0, 0, 2, 2, 0, 1, 4, 5...
$ CCG_3ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_3ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_5ltr_abs             <dbl> 1, 0, 1, 4, 8, 3, 1, 11, ...
$ CG_5ltr_rel             <dbl> 0.009433962, 0.000000000,...
$ CHG_5ltr_abs            <dbl> 2, 24, 9, 8, 14, 14, 9, 9...
$ CHG_5ltr_rel            <dbl> 0.01886792, 0.03927987, 0...
$ CHH_5ltr_abs            <dbl> 18, 69, 18, 26, 43, 70, 2...
$ CHH_5ltr_rel            <dbl> 0.16981132, 0.11292962, 0...
$ CCG_5ltr_abs            <dbl> 0, 0, 0, 2, 2, 0, 1, 4, 5...
$ CCG_5ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_5ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ cn_3ltr                 <dbl> NA, NA, NA, NA, NA, NA, N...
$ cn_5ltr                 <dbl> NA, NA, NA, NA, NA, NA, N...
```

### Output file format of `LTRpred()`

The columns of the prediction file store the following information:

tidy format (each row contains all information for a predicted LTR retrotransposon):

- `species` : Name of the species or genomic sequence that is passed as input to `LTRpred`
- `ID` : unique identifier that is generated by `LTRpred` to mark individual LTR retrotransposon predictions
- `dfam_target_name` : name of hmmer search match in Dfam database (see [LTRpred HMM Models](#hmm-models) for details)
- `ltr_similarity` : sequence similarity (global alignment) between the 5' and 3' long terminal repeat (LTR)
- `ltr_age_mya`: Rough estimation of retrotransposon insertion age in Mya based on 5 prime and 3 prime LTR sequence homology
- `ltr_evo_distance`: the molecular evolutionary distance value returned by the `evolutionary model` (e.g. Kimura-2P-Model) used to estamate DNA distances 
- `similarity`: `ltr_similarity` classified into 2% intervals
- `chromosome` : chromosome at which LTR retrotransposon was found
- `start` : start coordinate of the LTR retrotransposon within the chromosome
- `end` : end coordinate of the LTR retrotransposon within the chromosome
- `strand` : strand on which the corresponding LTR retrotransposon was found
- `dfam_acc` : Dfam accession
- `width` : length/width of the LTR retrotransposon 
- `orfs` : number of predicted open reading frames within the LTR retrotransposon sequence
- `protein_domain` : name of protein domains found between the 5' and 3' LTRs (see [LTRpred HMM Models](#hmm-models) for details)
- `Clust_Cluster` : cluster identifier 
- `Clust_Target` : LTR retrotransposon identifier that builds the cluster center of cluster `Clust_Cluster`
- `Clust_Perc_Ident` : sequence similarity (global alignment) with `Clust_Target` 
- `Clust_cn` : number of LTR retrotransposon that belong to the same cluster (see `Clust_Target`) - if `Clust_cn` = 1 -> unique element
- `lLTR_start` : start coordinate of the 5' LTR within the chromosome
- `lLTR_end` : end coordinate of the 5' LTR within the chromosome
- `lLTR_length` : length/width of the 5' LTR 
- `rLTR_start` : start coordinate of the 3' LTR within the chromosome
- `rLTR_end` : end coordinate of the 3' LTR within the chromosome
- `rLTR_length` : length/width of the 3' LTR
- `lTSD_start` : start coordinate of the 5' target site duplication (TSD) within the chromosome
- `lTSD_end` : end coordinate of the 5' target site duplication (TSD) within the chromosome
- `lTSD_motif` : 5' TSD sequence motif
- `PPT_start` : start coordinate of the PPT motif within the chromosome
- `PPT_end` : end coordinate of the PPT motif within the chromosome
- `PPT_motif` : PPT sequence motif
- `PPT_strand` : strand of PPT sequence 
- `PPT_offset` : width between `PPT_end` and `rLTR_start`
- `PBS_start` : start coordinate of the PBS motif within the chromosome
- `PBS_end` : end coordinate of the PBS motif within the chromosome
- `PBS_strand` : strand of PBS sequence
- `tRNA` : name of tRNA that matches PBS sequence
- `tRNA_motif` : tRNA sequence motif
- `PBS_offset` : width between `lLTR_end` and `PBS_start`
- `tRNA_offset` : number of nucleotide that are allowed for tRNA to offset the PBS
- `PBS/tRNA_edist` : edit distance between tRNA motif and PBS motif
- `PPT_length` : length/width of the PPT motif
- `PBS_length` : length/width of the PBS motif
- `dfam_target_description` : description of Dfam target hit
- `orf.id` : id generated by ORF prediction
- `TE_CG_abs` : total number of CG motifs found within the full length LTR retrotransposon sequence
- `TE_CG_rel` : relative number of CG motifs found within the full length LTR retrotransposon (normalized by element length - `width`)
- `TE_CHG_abs` : total number of CHG motifs found within the full length LTR retrotransposon sequence
- `TE_CHG_rel` : relative number of CHG motifs found within the full length LTR retrotransposon sequence (normalized by element length - `width`)
- `TE_CHH_abs` : total number of CHH motifs found within the full length LTR retrotransposon sequence
- `TE_CHH_rel` : relative number of CHH motifs found within the full length LTR retrotransposon sequence (normalized by element length - `width`)
- `TE_N_abs` : total number of N's found within the full length LTR retrotransposon sequence
- `CG_3ltr_abs` : total number of CG motifs found within only the 3' LTR sequence
- `CG_3ltr_rel` : relative number of CG motifs found within only the 3' LTR sequence (normalized by LTR length)
- `CHG_3ltr_abs` : total number of CHG motifs found within only the 3' LTR sequence
- `CHG_3ltr_rel` : relative number of CHG motifs found within only the 3' LTR sequence (normalized by LTR length)
- `CHH_3ltr_abs` : total number of CHH motifs found within only the 3' LTR sequence
- `CG_5ltr_abs` : total number of CG motifs found within only the 5' LTR sequence
- `CG_5ltr_rel` : relative number of CG motifs found within only the 5' LTR sequence (normalized by LTR length)
- `N_3ltr_abs` : total number of N's found within only the 3' LTR sequence
- `CHG_5ltr_abs` : total number of CHG motifs found within only the 5' LTR sequence
- `CHG_5ltr_rel` : relative number of CHG motifs found within only the 5' LTR sequence (normalized by LTR length)
- `CHH_5ltr_abs` : total number of CHH motifs found within only the 5' LTR sequence
- `CHH_5ltr_rel` : relative number of CHH motifs found within only the 5' LTR sequence (normalized by LTR length)
- `N_5ltr_abs` : total number of N's found within only the 5' LTR sequence
- `cn_3ltr` : number of solo 3' LTRs found in the genome
- `cn_5ltr` : number of solo 5' LTRs found in the genome (minus number of solo 3' LTRs)

## Detailed LTRpred run

The `Quick Start` section aimed at showing you an example `LTRpred` run and the corresponding output files and file formats.
In this section, users will learn about the input data, the available parameters that can be altered and about the diverse set of functions to perform functional annotation of LTR retrotransposons. Please make sure that you have all [prerequisite tools](https://github.com/HajkD/LTRpred/blob/master/vignettes/Installation.Rmd) installed on your system.

As an example, we will run LTRpred on the `Saccharomyces cerevisiae` genome. 

__Please be aware that functional de novo annotation is a computationally heavy task that requires substantial computing resources. For larger genomes recommend running LTRpred on a computing server or a high performance computing cluster. Large genome computations can take up to days even with tens of computing cores.__

The following example, however, will terminate within a few minutes. 

First we use the [biomartr](https://github.com/ropensci/biomartr) package to download the genome of `Saccharomyces cerevisiae` from NCBI RefSeq: 

```r
# retrieve S. cerevisiae genome from from NCBI RefSeq
Scerevisiae_genome <- biomartr::getGenome(db = "refseq", 
                      organism  = "Saccharomyces cerevisiae")
# show path to S. cerevisiae genome
Scerevisiae_genome
```

Next, we run `LTRpred` to generate functional annotation for `S. cerevisiae` LTR retrotransposons. In general, the `LTRpred()` function takes the path to a genome file in `fasta` format (compressed files can be used as well) as input to the `genome.file` argument. The `trnas` argument takes a file path to a `fasta` file storing tRNA sequences.
In this case, for example reasons only, we use the tRNA file `sacCer3-tRNAs.fa` that comes with the `LTRpred` package. In addition, users can find tRNA files at http://gtrnadb.ucsc.edu/ or http://gtrnadb2009.ucsc.edu/. The `hmms` argument takes a path to a folder storing the HMM models of the gag and pol protein domains. Protein domains
can be found at http://pfam.xfam.org/ . When renaming the protein domains using the specification `hmm_*` users can specify `hmm_*` which indicates that all files in the folder starting with `hmm_` shall be used. The `LTRpred` package already stores the most prominent gag and pol protein domains and can be used by specifying the path `paste0(system.file("HMMs/", package = "LTRpred"), "hmm_*")`. In detail, the following HMM models are stored in `LTRpred`:

#### HMM Models:

We retrieved the HMM models for protein domain annotation of the region
between de novo predicted LTRs from [Pfam](http://pfam.xfam.org):

  - RNA dependent RNA polymerase: [Overview](http://pfam.xfam.org/clan/CL0027)
      - [RdRP_1](http://pfam.xfam.org/family/PF00680#tabview=tab6)
      - [RdRP_2](http://pfam.xfam.org/family/PF00978#tabview=tab6)
      - [RdRP_3](http://pfam.xfam.org/family/PF00998#tabview=tab6)
      - [RdRP_4](http://pfam.xfam.org/family/PF02123#tabview=tab6)
      - [RVT_1](http://pfam.xfam.org/family/PF00078#tabview=tab6)
      - [RVT_2](http://pfam.xfam.org/family/PF07727#tabview=tab6)
      - [Integrase DNA binding domain](http://pfam.xfam.org/family/PF00552#tabview=tab6)
      - [Integrase Zinc binding domain](http://pfam.xfam.org/family/PF02022#tabview=tab6)
      - [Retrotrans_gag](http://pfam.xfam.org/family/PF03732#tabview=tab6)
      - [RNase H](http://pfam.xfam.org/family/PF00075#tabview=tab6)
      - [Integrase core domain](http://pfam.xfam.org/family/PF00665#tabview=tab6)
      - [Several Gypsy/Env proteins from Drosophila (Gypsy)](http://pfam.xfam.org/family/PF07253#tabview=tab6)
      - [ENV polyprotein (TLV coat polyprotein)](http://pfam.xfam.org/family/PF00429#tabview=tab6)
      
The argument `cluster = TRUE` indicates that individually predicted LTR retrotransposons are clustered by sequence similarity into families using the `vsearch` program. The `clust.sim`
defines the sequence similarity threshold for considering family/cluster members.

```r
library(LTRpred)
# de novo LTR transposon prediction of 'D. simulans'
LTRpred(
    genome.file = Scerevisiae_genome,
    trnas       = paste0(system.file("tRNAs/", package = "LTRpred"),"sacCer3-tRNAs.fa"),
    hmms        = paste0(system.file("HMMs/", package = "LTRpred"), "hmm_*"),
    cluster     = TRUE,
    clust.sim   = 0.9,
    copy.number.est = TRUE,  
    cores = 4
)
```

```
Running LTRpred on genome '_ncbi_downloads/genomes/Saccharomyces_cerevisiae_genomic_refseq.fna.gz' with 4 core(s) and searching for retrotransposons using the overlaps option (overlaps = 'no') ...


The output folder 'Saccharomyces_cerevisiae_genomic_refseq_ltrpred' does not seem to exist yet and will be created ...


LTRpred - Step 1:
Run LTRharvest...
LTRharvest: Generating index file Saccharomyces_cerevisiae_genomic_refseq_ltrharvest/Saccharomyces_cerevisiae_genomic_refseq_index.fsa with gt suffixerator...
Running LTRharvest and writing results to Saccharomyces_cerevisiae_genomic_refseq_ltrharvest...
LTRharvest analysis finished!


LTRpred - Step 2:
Generating index file Saccharomyces_cerevisiae_genomic_refseq_ltrdigest/Saccharomyces_cerevisiae_genomic_refseq_index_ltrdigest.fsa with suffixerator...
LTRdigest: Sort index file...
Running LTRdigest and writing results to Saccharomyces_cerevisiae_genomic_refseq_ltrdigest...
LTRdigest analysis finished!


LTRpred - Step 3:
Import LTRdigest Predictions...

Input:  Saccharomyces_cerevisiae_genomic_refseq_ltrdigest/Saccharomyces_cerevisiae_genomic_refseq_LTRdigestPrediction.gff  -> Row Number:  283
Remove 'NA' -> New Row Number:  283
(1/8) Filtering for repeat regions has been finished.
(2/8) Filtering for LTR retrotransposons has been finished.
(3/8) Filtering for inverted repeats has been finished.
(4/8) Filtering for LTRs has been finished.
(5/8) Filtering for target site duplication has been finished.
(6/8) Filtering for primer binding site has been finished.
(7/8) Filtering for protein match has been finished.
(8/8) Filtering for RR tract has been finished.


LTRpred - Step 4:
Perform ORF Prediction using 'usearch -fastx_findorfs' ...
usearch v8.1.1861_i86osx32, 4.0Gb RAM (17.2Gb total), 8 cores
(C) Copyright 2013-15 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

00:00 2.2Mb  100.0% Working
Join ORF Prediction table: nrow(df) = 36 candidates.
unique(ID) = 36 candidates.
unique(orf.id) = 36 candidates.
Perform clustering of similar LTR transposons using 'vsearch --cluster_fast' ...
Running CLUSTpred with 90% as sequence similarity threshold using 4 cores ...
Reading file /Users/h-gd/Desktop/Projekte/R Packages/Packages/LTRpred/Saccharomyces_cerevisiae_genomic_refseq_ltrdigest/Saccharomyces_cerevisiae_genomic_refseq-ltrdigest_complete.fas 100%
248182 nt in 36 seqs, min 5189, max 24168, avg 6894
Sorting by length 100%
Counting unique k-mers 100%
Clustering 100%
Sorting clusters 100%
Writing clusters 100%
Clusters: 10 Size min 1, max 19, avg 3.6
Singletons: 7, 19.4% of seqs, 70.0% of clusters
Sorting clusters by abundance 100%
vsearch v1.10.2_osx_x86_64, 16.0GB RAM, 8 cores
https://github.com/torognes/vsearch

CLUSTpred output has been stored in: Saccharomyces_cerevisiae_genomic_refseq_ltrpred
Join Cluster table: nrow(df) = 36 candidates.
unique(ID) = 36 candidates.
unique(orf.id) = 36 candidates.
Join Cluster Copy Number table: nrow(df) = 36 candidates.
unique(ID) = 36 candidates.
unique(orf.id)) = 36 candidates.


LTRpred - Step 5:
Perform methylation context quantification..
Join methylation context (CG, CHG, CHH, CCG) count table: nrow(df) = 36 candidates.
unique(ID) = 36 candidates.
unique(orf.id) = 36 candidates.
Copy files to result folder 'Saccharomyces_cerevisiae_genomic_refseq_ltrpred'.


LTRpred - Step 6:
Starting retrotransposon evolutionary age estimation by comparing the 3' and 5' LTRs using the molecular evolution model 'K80' and the mutation rate '1.3e-07' (please make sure the mutation rate can be assumed for your species of interest!) for 36 predicted elements ...


Please be aware that evolutionary age estimation based on 3' and 5' LTR comparisons are only very rough time estimates and don't take reverse-transcription mediated retrotransposon recombination between family members of retroelements into account! Please consult Sanchez et al., 2017 Nature Communications and Drost & Sanchez, 2019 Genome Biology and Evolution for more details on retrotransposon recombination.


LTRpred - Step 7:
The LTRpred prediction table has been filtered (default) to remove potential false positives. Predicted LTRs must have an PBS or Protein Domain and must fulfill thresholds: sim = 70%; #orfs = 0. Furthermore, TEs having more than 10% of N's in their sequence have also been removed.
Input #TEs: 36
Output #TEs: 31
Perform solo LTR Copy Number Estimation....
Run makeblastdb of the genome assembly...


Building a new DB, current time: 01/24/2020 14:58:55
New DB name:   _ncbi_downloads/genomes/Saccharomyces_cerevisiae_genomic_refseq.fna
New DB title:  _ncbi_downloads/genomes/Saccharomyces_cerevisiae_genomic_refseq.fna
Sequence type: Nucleotide
Keep Linkouts: T
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 17 sequences in 0.129384 seconds.
Perform BLAST searches of 3' prime LTRs against genome assembly...
Perform BLAST searches of 5' prime LTRs against genome assembly...
Import BLAST results...
Filter hit results...
Estimate CNV for each LTR sequence...
Finished LTR CNV estimation!


TRpred finished all analyses successfully. All output files were stored at 'Saccharomyces_cerevisiae_genomic_refseq_ltrpred'.
[1] "Successful job 1 ."
Warning message:
The LTR copy number estimation returned an empty file. This suggests that there were no solo LTRs found in the input genome sequence. 
```



The output can then be imported using:

```r
# import LTRpred output for S. cerevisiae
file <- file.path("Saccharomyces_cerevisiae_genomic_refseq_ltrpred",
                  "Saccharomyces_cerevisiae_genomic_refseq_LTRpred_DataSheet.tsv")
Scerevisiae_LTRpred <- read.ltrpred(file)
# look at output
dplyr::glimpse(Scerevisiae_LTRpred)
```

```
Observations: 31
Variables: 92
$ species                 <chr> "Saccharomyces_cerevisiae...
$ ID                      <chr> "Saccharomyces_cerevisiae...
$ dfam_target_name        <chr> NA, NA, NA, NA, NA, NA, N...
$ ltr_similarity          <dbl> 100.00, 97.15, 99.40, 96....
$ ltr_age_mya             <dbl> 0.00000000, 0.05613384, 0...
$ similarity              <chr> "(98,100]", "(96,98]", "(...
$ protein_domain          <chr> "rve/RVT_2", "rve/RVT_2",...
$ orfs                    <int> 2, 2, 2, 2, 2, 3, 2, 2, 2...
$ chromosome              <chr> "NC", "NC", "NC", "NC", "...
$ start                   <int> 160238, 221030, 259578, 2...
$ end                     <int> 166162, 226960, 265494, 2...
$ strand                  <chr> "-", "+", "+", "-", "+", ...
$ width                   <int> 5925, 5931, 5917, 5927, 5...
$ annotation              <chr> "LTR_retrotransposon", "L...
$ pred_tool               <chr> "LTRpred", "LTRpred", "LT...
$ frame                   <chr> ".", ".", ".", ".", ".", ...
$ score                   <chr> ".", ".", ".", ".", ".", ...
$ lLTR_start              <int> 160238, 221030, 259578, 2...
$ lLTR_end                <int> 160574, 221380, 259909, 2...
$ lLTR_length             <int> 337, 351, 332, 350, 380, ...
$ rLTR_start              <int> 165826, 226615, 265163, 2...
$ rLTR_end                <int> 166162, 226960, 265494, 2...
$ rLTR_length             <int> 337, 346, 332, 339, 369, ...
$ lTSD_start              <int> 160233, 221026, 259573, 2...
$ lTSD_end                <int> 160237, 221029, 259577, 2...
$ lTSD_motif              <chr> "gaacc", "aaca", "gtaat",...
$ rTSD_start              <int> 166163, 226961, 265495, 2...
$ rTSD_end                <int> 166167, 226964, 265499, 2...
$ rTSD_motif              <chr> "gaacc", "aaca", "gtaat",...
$ PPT_start               <int> NA, NA, NA, NA, NA, NA, N...
$ PPT_end                 <int> NA, NA, NA, NA, NA, NA, N...
$ PPT_motif               <chr> NA, NA, NA, NA, NA, NA, N...
$ PPT_strand              <chr> NA, NA, NA, NA, NA, NA, N...
$ PPT_offset              <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_start               <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_end                 <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_strand              <chr> NA, NA, NA, NA, NA, NA, N...
$ tRNA                    <chr> NA, NA, NA, NA, NA, NA, N...
$ tRNA_motif              <chr> NA, NA, NA, NA, NA, NA, N...
$ PBS_offset              <int> NA, NA, NA, NA, NA, NA, N...
$ tRNA_offset             <int> NA, NA, NA, NA, NA, NA, N...
$ `PBS/tRNA_edist`        <int> NA, NA, NA, NA, NA, NA, N...
$ orf.id                  <chr> "NC_001133.9_Saccharo_160...
$ repeat_region_length    <int> 5935, 5939, 5927, 5935, 5...
$ PPT_length              <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_length              <int> NA, NA, NA, NA, NA, NA, N...
$ dfam_acc                <chr> NA, NA, NA, NA, NA, NA, N...
$ dfam_bits               <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_e_value            <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_bias               <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_strand             <chr> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_modlen             <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_target_description <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Cluster           <chr> "cl_4", "cl_4", "cl_4", "...
$ Clust_Target            <chr> "NC_001140.6_Saccharo_543...
$ Clust_Perc_Ident        <dbl> 94.1, 92.4, 95.9, 96.6, 9...
$ Clust_cn                <int> 19, 19, 19, 19, 19, 19, 7...
$ TE_CG_abs               <dbl> 160, 157, 168, 167, 156, ...
$ TE_CG_rel               <dbl> 0.02700422, 0.02647108, 0...
$ TE_CHG_abs              <dbl> 187, 192, 197, 199, 193, ...
$ TE_CHG_rel              <dbl> 0.03156118, 0.03237228, 0...
$ TE_CHH_abs              <dbl> 922, 919, 901, 914, 936, ...
$ TE_CHH_rel              <dbl> 0.1556118, 0.1549486, 0.1...
$ TE_CCG_abs              <dbl> 33, 38, 38, 36, 34, 35, 2...
$ TE_CCG_rel              <dbl> 0.005569620, 0.006407014,...
$ TE_N_abs                <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_3ltr_abs             <dbl> 6, 5, 6, 7, 3, 6, 5, 5, 6...
$ CG_3ltr_rel             <dbl> 0.017804154, 0.014450867,...
$ CHG_3ltr_abs            <dbl> 2, 4, 3, 3, 3, 4, 2, 2, 3...
$ CHG_3ltr_rel            <dbl> 0.005934718, 0.011560694,...
$ CHH_3ltr_abs            <dbl> 41, 40, 36, 37, 45, 40, 4...
$ CHH_3ltr_rel            <dbl> 0.1216617, 0.1156069, 0.1...
$ CCG_3ltr_abs            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CCG_3ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_3ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_5ltr_abs             <dbl> 6, 5, 6, 7, 3, 6, 5, 5, 6...
$ CG_5ltr_rel             <dbl> 0.017804154, 0.014450867,...
$ CHG_5ltr_abs            <dbl> 2, 4, 3, 3, 3, 4, 2, 2, 3...
$ CHG_5ltr_rel            <dbl> 0.005934718, 0.011560694,...
$ CHH_5ltr_abs            <dbl> 41, 40, 36, 37, 45, 40, 4...
$ CHH_5ltr_rel            <dbl> 0.1216617, 0.1156069, 0.1...
$ CCG_5ltr_abs            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CCG_5ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_5ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ cn_3ltr                 <dbl> NA, NA, NA, NA, NA, NA, N...
$ cn_5ltr                 <dbl> NA, NA, NA, NA, NA, NA, N...
```


Now, users can filter for potentially functional LTR retrotransposons
using the `quality.filter()` function. The argument `sim` specifies the miminum sequence similarity between the LTRs of the predicted retrotransposons. E.g., when `sim = 97` all retrotransposons having less than 97% sequence similarity between their LTRs are removed.
Argument `n.orfs` specifies the minimum number of open reading frames (ORFs) that were predicted to exist in the predicted retrotransposon. E.g., when `n.orfs = 1` all retrotransposon having zero ORFs are removed. The argument `strategy` refers to the filter strategy that shall be applied. Options are:

- `strategy = "default"` : 
    - `ltr.similarity`: Minimum similarity between LTRs. All TEs not matching this criteria are discarded.
    - `n.orfs`: minimum number of ORFs that must be found between the LTRs. All TEs not matching this criteria are discarded.
    - `PBS` or Protein Match: elements must either have a predicted Primer Binding Site or a protein match of at least one protein (Gag, Pol, Rve, ...) between their LTRs. All TEs not matching this criteria are discarded.
    - The relative number of N's (= nucleotide not known) in TE <= 0.1. The relative number of N's is computed as follows: absolute number of N's in TE / width of TE.

- `strategy = "stringent"` : in addition to filter criteria specified in section Quality Control, the filter criteria !is.na(protein_domain)) | (dfam_target_name != "unknown") is applied

```r
# filter for potentially functional LTR retrotransposons
Scerevisiae_functional <- quality.filter(Scerevisiae_LTRpred, 
                                         sim = 97, 
                                         n.orfs = 1,
                                         strategy = "stringent")
dplyr::glimpse(Scerevisiae_functional)
```

```
Observations: 25
Variables: 92
$ species                 <chr> "Saccharomyces_cerevisiae...
$ ID                      <chr> "Saccharomyces_cerevisiae...
$ dfam_target_name        <chr> NA, NA, NA, NA, NA, NA, N...
$ ltr_similarity          <dbl> 100.00, 97.15, 99.40, 100...
$ ltr_age_mya             <dbl> 0.00000000, 0.05613384, 0...
$ similarity              <chr> "(98,100]", "(96,98]", "(...
$ protein_domain          <chr> "rve/RVT_2", "rve/RVT_2",...
$ orfs                    <int> 2, 2, 2, 2, 2, 2, 2, 2, 2...
$ chromosome              <chr> "NC", "NC", "NC", "NC", "...
$ start                   <int> 160238, 221030, 259578, 9...
$ end                     <int> 166162, 226960, 265494, 9...
$ strand                  <chr> "-", "+", "+", "+", "-", ...
$ width                   <int> 5925, 5931, 5917, 5959, 5...
$ annotation              <chr> "LTR_retrotransposon", "L...
$ pred_tool               <chr> "LTRpred", "LTRpred", "LT...
$ frame                   <chr> ".", ".", ".", ".", ".", ...
$ score                   <chr> ".", ".", ".", ".", ".", ...
$ lLTR_start              <int> 160238, 221030, 259578, 9...
$ lLTR_end                <int> 160574, 221380, 259909, 9...
$ lLTR_length             <int> 337, 351, 332, 332, 332, ...
$ rLTR_start              <int> 165826, 226615, 265163, 9...
$ rLTR_end                <int> 166162, 226960, 265494, 9...
$ rLTR_length             <int> 337, 346, 332, 332, 332, ...
$ lTSD_start              <int> 160233, 221026, 259573, 9...
$ lTSD_end                <int> 160237, 221029, 259577, 9...
$ lTSD_motif              <chr> "gaacc", "aaca", "gtaat",...
$ rTSD_start              <int> 166163, 226961, 265495, 9...
$ rTSD_end                <int> 166167, 226964, 265499, 9...
$ rTSD_motif              <chr> "gaacc", "aaca", "gtaat",...
$ PPT_start               <int> NA, NA, NA, NA, NA, NA, N...
$ PPT_end                 <int> NA, NA, NA, NA, NA, NA, N...
$ PPT_motif               <chr> NA, NA, NA, NA, NA, NA, N...
$ PPT_strand              <chr> NA, NA, NA, NA, NA, NA, N...
$ PPT_offset              <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_start               <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_end                 <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_strand              <chr> NA, NA, NA, NA, NA, NA, N...
$ tRNA                    <chr> NA, NA, NA, NA, NA, NA, N...
$ tRNA_motif              <chr> NA, NA, NA, NA, NA, NA, N...
$ PBS_offset              <int> NA, NA, NA, NA, NA, NA, N...
$ tRNA_offset             <int> NA, NA, NA, NA, NA, NA, N...
$ `PBS/tRNA_edist`        <int> NA, NA, NA, NA, NA, NA, N...
$ orf.id                  <chr> "NC_001133.9_Saccharo_160...
$ repeat_region_length    <int> 5935, 5939, 5927, 5969, 5...
$ PPT_length              <int> NA, NA, NA, NA, NA, NA, N...
$ PBS_length              <int> NA, NA, NA, NA, NA, NA, N...
$ dfam_acc                <chr> NA, NA, NA, NA, NA, NA, N...
$ dfam_bits               <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_e_value            <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_bias               <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_hmm-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_strand             <chr> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_ali-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-st`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ `dfam_env-en`           <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_modlen             <dbl> NA, NA, NA, NA, NA, NA, N...
$ dfam_target_description <chr> NA, NA, NA, NA, NA, NA, N...
$ Clust_Cluster           <chr> "cl_4", "cl_4", "cl_4", "...
$ Clust_Target            <chr> "NC_001140.6_Saccharo_543...
$ Clust_Perc_Ident        <dbl> 94.1, 92.4, 95.9, 98.8, 9...
$ Clust_cn                <int> 19, 19, 19, 7, 19, 19, 19...
$ TE_CG_abs               <dbl> 160, 157, 168, 152, 159, ...
$ TE_CG_rel               <dbl> 0.02700422, 0.02647108, 0...
$ TE_CHG_abs              <dbl> 187, 192, 197, 155, 189, ...
$ TE_CHG_rel              <dbl> 0.03156118, 0.03237228, 0...
$ TE_CHH_abs              <dbl> 922, 919, 901, 912, 922, ...
$ TE_CHH_rel              <dbl> 0.1556118, 0.1549486, 0.1...
$ TE_CCG_abs              <dbl> 33, 38, 38, 26, 34, 32, 3...
$ TE_CCG_rel              <dbl> 0.005569620, 0.006407014,...
$ TE_N_abs                <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_3ltr_abs             <dbl> 6, 5, 6, 5, 6, 6, 6, 4, 6...
$ CG_3ltr_rel             <dbl> 0.017804154, 0.014450867,...
$ CHG_3ltr_abs            <dbl> 2, 4, 3, 2, 3, 3, 3, 2, 2...
$ CHG_3ltr_rel            <dbl> 0.005934718, 0.011560694,...
$ CHH_3ltr_abs            <dbl> 41, 40, 36, 40, 38, 40, 3...
$ CHH_3ltr_rel            <dbl> 0.12166172, 0.11560694, 0...
$ CCG_3ltr_abs            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CCG_3ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_3ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CG_5ltr_abs             <dbl> 6, 5, 6, 5, 6, 6, 6, 4, 6...
$ CG_5ltr_rel             <dbl> 0.017804154, 0.014450867,...
$ CHG_5ltr_abs            <dbl> 2, 4, 3, 2, 3, 3, 3, 2, 2...
$ CHG_5ltr_rel            <dbl> 0.005934718, 0.011560694,...
$ CHH_5ltr_abs            <dbl> 41, 40, 36, 40, 38, 40, 3...
$ CHH_5ltr_rel            <dbl> 0.12166172, 0.11560694, 0...
$ CCG_5ltr_abs            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ CCG_5ltr_rel            <dbl> 0.000000000, 0.000000000,...
$ N_5ltr_abs              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0...
$ cn_3ltr                 <dbl> NA, NA, NA, NA, NA, NA, N...
$ cn_5ltr                 <dbl> NA, NA, NA, NA, NA, NA, N...
```

The functional annotation table can then be transformed and saved in different
data formats using the `pred2*()` functions. Options are:

- `pred2bed()`:	Format LTR prediction data to BED file format
- `pred2csv()`:	Format LTR prediction data to CSV file format
- `pred2fasta()`:	Save the sequence of the predicted LTR Transposons in a fasta file
- `pred2gff()`:	Format LTR prediction data to GFF3 file format
- `pred2tsv()`:	Format LTR prediction data to TSV file format

In this example we will store the functional annotation results `Scerevisiae_functional` as `gff` file.

```r
pred2gff(LTR.data = Scerevisiae_functional,
         output   = "Scerevisiae_functional_LTRpred.gff",
         program  = "LTRpred")
```

This gff file can now be used for mapping tools, genome browsers, etc.


### Detailed description of adjustable LTRpred parameters

The `S. cerevisiae` example shown above assumes that users wish to run `LTRpred`
using a default parameter configuration. However, as the figure shown below demonstrates, there are  large number of parameters that can be adjusted and altered according the
the user's specifications and interests.

Users can adjust and run `LTRpred` with the following parameter options:

```r
LTRpred(
  genome.file = NULL,
  model = "K80",
  mutation_rate = 1.3 * 1e-07,
  index.file.harvest = NULL,
  index.file.digest = NULL,
  LTRdigest.gff = NULL,
  tabout.file = NULL,
  LTRharvest.folder = NULL,
  LTRpred.folder = NULL,
  orf.file = NULL,
  annotate = NULL,
  Dfam.db = NULL,
  dfam.eval = 0.001,
  dfam.file = NULL,
  cluster = FALSE,
  clust.sim = 0.9,
  clust.file = NULL,
  copy.number.est = TRUE,
  fix.chr.name = TRUE,
  cn.eval = 1e-10,
  range = c(0, 0),
  seed = 30,
  minlenltr = 100,
  maxlenltr = 3500,
  mindistltr = 4000,
  maxdistltr = 25000,
  similar = 70,
  mintsd = 4,
  maxtsd = 20,
  vic = 60,
  overlaps = "no",
  xdrop = 5,
  mat = 2,
  mis = -2,
  ins = -3,
  del = -3,
  motif = NULL,
  motifmis = 0,
  aaout = "yes",
  aliout = "yes",
  pptlen = c(8, 30),
  uboxlen = c(3, 30),
  pptradius = 30,
  trnas = NULL,
  pbsalilen = c(11, 30),
  pbsoffset = c(0,5),
  pbstrnaoffset = c(0, 5),
  pbsmaxedist = 1,
  pbsradius = 30,
  hmms = NULL,
  pdomevalcutoff = 1e-05,
  pbsmatchscore = 5,
  pbsmismatchscore = -10,
  pbsinsertionscore = -20,
  pbsdeletionscore = -20,
  pfam.ids = NULL,
  cores = 1,
  dfam.cores = NULL,
  hmm.cores = NULL,
  orf.style = 7,
  min.codons = 200,
  trans.seqs = FALSE,
  output.path = NULL,
  quality.filter = TRUE,
  n.orfs = 0,
  verbose = TRUE
)
```

- `genome.file`:	path to the genome assembly file in fasta format.
- `model`: a molecular evolution model to estimate the insertion age of the predicted retrotransposon
- `mutation_rate`: mutation rate applied to the organism of interest when applying the molecular evolution model to estimate the insertion age of the predicted retrotransposon
- `index.file.harvest`:	in case users have already computed index files they can specify the name of the enhanced suffix array index file that was computed by `suffixerator` for the use of `LTRharvest`. This option can be used in case the suffix file was previously generated, e.g. during a previous call of this function. In this case the suffix array index file does not need to be re-computed for new analyses. This is particularly useful when running `LTRpred` with different parameter settings.
- `index.file.digest`: specify the name of the enhanced suffix array index file that is computed by `suffixerator` for the use of `LTRdigest`. This option can be used in case the suffix file was previously generated, e.g. during a previous call of this function. In this case the suffix array index file does not need to be re-computed for new analyses. This is particularly useful when running `LTRpred` with different parameter settings.
- `LTRdigest.gff`:	path to the LTRdigest generated GFF file `*_LTRdigestPrediction.gff`, in case `LTRdigest` files were computed previously.
- `tabout.file`: path to the LTRdigest generated tabout file file `*-ltrdigest_tabout.csv`, in case `LTRdigest` files were computed previously.
- `LTRharvest.folder`:	either `LTRharvest.folder = NULL` (default) or `LTRharvest.folder = "skip"` to skip a LTRharvest folder if it is not present.
- `LTRpred.folder`: name/path of/to an existing LTRpred folder. All pre-computed files in this folder will then be used and only steps of the `LTRpred` pipeline that haven't been computed yet will be run.
- `orf.file`: path to the file generated by `ORFpred`, in case the orf prediction file was generated previously. See `?ORFpred` for details.
- `annotate`: annotation database that shall be queried to annotate predicted LTR transposons. Default is `annotate = NULL` indicating that no annotation query is performed. Possible options are: 
    - `annotate = "Dfam"` (here the Dfam database must be stored locally and a `nhammer` search is performed against the Dfam database)
- `Dfam.db`: folder path to the local Dfam database or `Dfam.db = "download"` in case the `Dfam` database shall be automatically downloaded before performing query analyses.
- `dfam.eval`: E-value threshhold to perform HMMer search against the Dfam database.
- `dfam.file`: path to pre-computed `dfam.query` output file. Can only be used in combination with `annotate = "Dfam"`. Please consult `?dfam.query` for details.
- `cluster`: shall predicted transposons be clustered with `CLUSTpred` to determine retrotransposon family relationships?
- `clust.sim`: cluster reject if sequence similarity is lower than this threshold when performing clustering with `CLUSTpred`. In other words, what is the minimum sequence similarity between family members. Default is `clust.sim = 0.9` (= 90%).
- `clust.file`: file path to pre-computed clustering file generated by `CLUSTpred` in case clustering has been done before. See `?CLUSTpred` for details.
- `copy.number.est`: shall copy number estimation of 3' and 5' LTRs in the genome be performed? In case `copy.number.est = TRUE` for each  3' and 5' LTR of a predicted LTR retrotransposon a stringent BLAST search is performed against the genome to determine the number of copies (solo LTR copies) of that element within the genome. Default is `copy.number.est = FALSE`.
- `fix.chr.name`:	sometimes chromosome names have unwanted characters such as "_" etc in them. Shall `LTRpred` try to fix chromosome names as good as possible?
- `cn.eval`: 	BLAST evalue for copy number estimation (`copy.number.est`) (BLAST hit threshold). Default is `cn.eval = 1e-5`.
- `range`: define the genomic interval within the genome assembly in which predicted LTR transposons shall be detected. In case `range[1] = 1000` and `range[2] = 10000` then candidates are only reported if they start after 1000 numcleotides in the genome assembly and end before position 10000 in their respective sequence coordinates. If `range[1] = 0` and `range[2] = 0`, so `range = c(0,0)` (default) then the entire genome is being scanned.
- `seed`:	the minimum length for the exact maximal repeats. Only repeats with the specified minimum length are considered in all subsequent analyses. Default is `seed = 30`.
- `minlenltr`:	minimum LTR length. Default is `minlenltr = 100`.
- `maxlenltr`:	maximum LTR length. Default is `maxlenltr = 3500`.
- `mindistltr`:	minimum distance of LTR starting positions. Default is `mindistltr = 4000`.
- `maxdistltr`:	maximum distance of LTR starting positions. Default is `maxdistltr = 25000`.
- `similar`:	minimum similarity value between the two LTRs in percent. Default is `similar = 70`.
- `mintsd`: minimum target site duplications (TSDs) length. If no search for TSDs shall be performed, then specify mintsd = NULL. Default is `mintsd = 4`.
- `maxtsd`:	maximum target site duplications (TSDs) length. If no search for TSDs shall be performed, then specify maxtsd = NULL. Default is `maxtsd = 20`.
- `vic`:	number of nucleotide positions left and right (the vicinity) of the predicted boundary of a LTR that will be searched for TSDs and/or one motif (if specified). Default is `vic = 60`.
- `overlaps`:	specify how overlapping LTR retrotransposon predictions shall be treated. If overlaps = "no" is selected, then neither nested nor overlapping predictions will be reported in the output. In case `overlaps = "best"` is selected then in the case of two or more nested or overlapping predictions, solely the LTR retrotransposon prediction with the highest similarity between its LTRs will be reported. If `overlaps = "all"` is selected then all LTR retrotransposon predictions will be reported whether there are nested and/or overlapping predictions or not. Default is `overlaps = "best"`.
- `xdrop`:	specify the xdrop value (> 0) for extending a seed repeat in both directions allowing for matches, mismatches, insertions, and deletions. The xdrop extension process stops as soon as the extension involving matches, mismatches, insersions, and deletions has a score smaller than T - X, where T denotes the largest score seen so far. Default is `xrop = 5`.
- `mat`:	specify the positive match score for the X-drop extension process. Default is `mat = 2`.
- `mis`:	specify the negative mismatch score for the X-drop extension process. Default is `mis = -2`.
- `ins`:	specify the negative insertion score for the X-drop extension process. Default is `ins = -3`.
- `del`:	specify the negative deletion score for the X-drop extension process. Default is `del = -3`.
- `motif`:	specify 2 nucleotides for the starting motif and 2 nucleotides for the ending motif at the beginning and the ending of each LTR, respectively. Only palindromic motif sequences - where the motif sequence is equal to its complementary sequence read backwards - are allowed, e.g. `motif = "tgca"`. Type the nucleotides without any space separating them. If this option is not selected by the user, candidate pairs will not be screened for potential motifs. If this options is set but no allowed number of mismatches is specified by the argument motifmis and a search for the exact motif will be conducted. If `motif = NULL` then no explicit motif is being specified.
- `motifmis`:	allowed number of mismatches in the TSD motif specified in motif. The number of mismatches needs to be between [0,3]. Default is `motifmis = 0`.
- `aaout`:	shall the protein sequence of the HMM matches to the predicted LTR transposon be generated as fasta file or not. Options are `aaout = "yes"` or `aaout = "no"`.
- `aliout`:	shall the alignment of the protein sequence of the HMM matches to the predicted LTR transposon be generated as fasta file or not. Options are `aaout = "yes"` or `aaout = "no"`.
- `pptlen`:	a two dimensional numeric vector specifying the minimum and maximum allowed lengths for PPT predictions. If a purine-rich region that does not fulfill this range is found, it will be discarded. Default is `pptlen = c(8,30)` (minimum = 8; maximum = 30).
- `uboxlen`:	a two dimensional numeric vector specifying the minimum and maximum allowed lengths for U-box predictions. If a T-rich region preceding a PPT that does not fulfill the PPT length criteria is found, it will be discarded. Default is `uboxlen = c(3,30)` (minimum = 3; maximum = 30).
- `pptradius`:	a numeric value specifying the area around the 3' LTR beginning to be considered when searching for PPT. Default value is `pptradius = 30`.
- `trnas`:	path to the fasta file storing the unique tRNA sequences that shall be matched to the predicted LTR transposon (tRNA library).
- `pbsalilen`:	a two dimensional numeric vector specifying the minimum and maximum allowed lengths for PBS/tRNA alignments. If the local alignments are shorter or longer than this range, it will be discarded. Default is `pbsalilen = c(11,30)` (minimum = 11; maximum = 30).
- `pbsoffset`:	a two dimensional numeric vector specifying the minimum and maximum allowed distance between the start of the PBS and the 3' end of the 5' LTR. Local alignments not fulfilling this criteria will be discarded. Default is `pbsoffset = c(0,5)` (minimum = 0; maximum = 5).
- `pbstrnaoffset`:	a two dimensional numeric vector specifying the minimum and maximum allowed PBS/tRNA alignment offset from the 3' end of the tRNA. Local alignments not fulfilling this criteria will be discarded. `Default is pbstrnaoffset = c(0,5)` (minimum = 0; maximum = 5).
- `pbsmaxedist`:	a numeric value specifying the maximal allowed unit edit distance in a local PBS/tRNA alignment.
- `pbsradius`:	a numeric value specifying the area around the 5' LTR end to be considered when searching for PBS Default value is `pbsradius = 30`.
- `hmms`:	a character string or a character vector storing either the hmm files for searching internal domains between the LTRs of predicted LTR transposons or a vector of Pfam IDs from http://pfam.xfam.org/ that are downloaded and used to search for corresponding protein domains within the predicted LTR transposons. As an option users can rename all of their hmm files so that they start for example with the name `hmms = "hmm_*"`. This way all files starting with `hmm_` will be considered for the subsequent protein domain search. In case Pfam IDs are specified, the LTRpred function will automatically download the corresponding HMM files and use them for further protein domain searches. In case users prefer to specify Pfam IDs please specify them in the `pfam.ids` parmeter and choose `hmms = NULL`.
- `pdomevalcutoff`:	a numeric value specifying the E-value cutoff for corresponding HMMER searches. All hits that do not fulfill this criteria are discarded. Default is `pdomevalcutoff = 1E-5`.
- `pbsmatchscore`:	specify the match score used in the PBS/tRNA Smith-Waterman alignment. Default is `pbsmatchscore = 5`.
- `pbsmismatchscore`: specify the mismatch score used in the PBS/tRNA Smith-Waterman alignment. Default is `pbsmismatchscore = -10`.
- `pbsinsertionscore`:	
specify the insertion score used in the PBS/tRNA Smith-Waterman alignment. Default is `pbsinsertionscore = -20`.
- `pbsdeletionscore`:	specify the deletion score used in the PBS/tRNA Smith-Waterman alignment. Default is `pbsdeletionscore = -20`.
- `pfam.ids`: a character vector storing the Pfam IDs from http://pfam.xfam.org/ that shall be downloaded and used to perform protein domain searches within the sequences between the predicted LTRs.
- `cores`:	the number of cores that shall be used for multicore processing. In case `dfam.cores` and `hmm.cores` are not specified then the value of `core` is used for those arguments.  
- `dfam.cores`:	number of cores to be used for multicore processing when running Dfam query (in case `annotate = "Dfam"`).
- `hmm.cores`: number of cores to be used for multicore processing when performing hmmer protein search with `LTRdigest`.
- `orf.style`: type of predicting open reading frames (see documentation of USEARCH).
- `min.codons`:	minimum number of codons in the predicted open reading frame.
- `trans.seqs`:	logical value indicating wheter or not predicted open reading frames shall be translated and the corresponding protein sequences stored in the output folder.
- `output.path`:	a path/folder to store all results returned by LTRharvest, LTRdigest, and LTRpred. If `output.path = NULL (Default)` then a folder with the name of the input genome file will be generated in the current working directory of R and all results are then stored in this folder.
- `quality.filter`:	shall false positives be filtered out as much as possible? Default is `quality.filter = TRUE`. See Description for details.
- `n.orfs`:	minimum number of Open Reading Frames that must be found between the LTRs (if `quality.filter = TRUE`). See Details for further information on quality control.
- `verbose`:	shall further information be printed on the console or not.

## Metagenome scale annotations

`LTRpred` allows users to perform annotations not only for single genomes
but for multiple genomes (metagenomes) using only one pipeline function named `LTRpred.meta()`. 

__Please be aware that `LTRpred` annotations for multiple genomes is a computationally hard task and requires large server or high performance computer access. Computations for hundreds of genomes even using tens to hundreds of cores might take several weeks to terminate. Please make sure that you have the right computational infrastructure to run these processes.__


Users can download the [biomartr](https://docs.ropensci.org/biomartr/) package to automatically retrieve genome assembly files for the species of interest.

```r
# specify the scientific names of the species of interest
# that shall first be downloaded and then be used
# to generate LTRpred annotations
species <- c("Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella") 
# download the genome assembly files for the species of interest
# 
# install.packages("biomartr")
biomartr::getGenomeSet(db = "refseq", organisms = species, path = "store_genome_set")
# run LTRpred.meta() on the 3 species with 3 cores
LTRpred::LTRpred.meta(genome.folder = "store_genome_set",
                      output.folder = "LTRpred_meta_results",
                      cores = 3)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename.organisms.R
\name{rename.organisms}
\alias{rename.organisms}
\title{Rename file path names from \code{read_proteome}, etc.}
\usage{
rename.organisms(data)
}
\arguments{
\item{data}{a GenomeInfo file.}
}
\description{
Rename a file path e.g. \code{Arabidopsis_thaliana.fa.gz} to \code{Athaliana}.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllPairwiseAlign.R
\name{AllPairwiseAlign}
\alias{AllPairwiseAlign}
\title{Compute all pairwise (global) alignments with VSEARCH}
\usage{
AllPairwiseAlign(
  file,
  cores = 1,
  mask = "none",
  out.name = "AllPairwiseAlign",
  output = NULL
)
}
\arguments{
\item{file}{path to fasta file storing sequences for which all possible pairwise alignments shall be computed.}

\item{cores}{number of cores that shall be used for parallel computations.}

\item{mask}{shall the aligned sequences be maksed? Options are: \code{mask = "none"}, \code{mask = "dust"}, \code{mask = "soft"}.}

\item{out.name}{name of the output files (\code{*.uc}, \code{*_userout.txt}, \code{*.blast6out}), \code{*_fasta_pairs.fasta}, \code{*.sam}.}

\item{output}{path to a folder in which output shall be stored.}
}
\value{
A folder named \code{vsearch_pairwise_align} will be created and the following files will be stored in this output folder:

\itemize{
\item \code{*.uc} USEARCH cluster format generated by VSEARCH storing the sequence cluster information of aligned sequences.
\item \code{*_userout.txt} a tab separated file storing the alignment information in the columns: \code{query}, \code{target}, \code{sequence identity (ID)}.
\item \code{*.blast6out} BLAST output format of the pairwise (global) sequence alignments generated by VSEARCH.
\item \code{*_fasta_pairs.fasta} fasta file storing the alignments.
\item \code{*.sam} SAM file storing the alignments.
}
}
\description{
This function is a wrapper function to compute all pairwise (global) alignments using VSEARCH.
}
\details{
To be able to use this function the VSEARCH command line tool needs to be installed.
}
\references{
https://github.com/torognes/vsearch
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename.fasta.meta.R
\name{rename.fasta.meta}
\alias{rename.fasta.meta}
\title{Meta function for applying \code{rename.fasta}}
\usage{
rename.fasta.meta(in.folder, out.file)
}
\arguments{
\item{in.folder}{path to folder storing \code{\link{LTRpred}} generated organism prediction folders.}

\item{out.file}{name of concatenated output file.}
}
\description{
Meta function for applying \code{\link{rename.fasta}} to many fasta files within a folder.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repbase.filter.R
\name{repbase.filter}
\alias{repbase.filter}
\title{Filter the Repbase query output}
\usage{
repbase.filter(query.output, scope.value = 0.7, verbose = TRUE)
}
\arguments{
\item{query.output}{a \code{data.frame} returned by the \code{\link{repbase.query}} function.}

\item{scope.value}{a value between [0,1] qunatifying the percentage of minimum sequence similariy 
between the LTR transposon and the corresponding annotated sequence found in Repbase.}

\item{verbose}{a logical value indicating whether or not additional information shall be printed 
to the console while executing this function.}
}
\value{
A \code{data.frame} storing the filtered output returned by \code{\link{repbase.query}}.
}
\description{
Filter the output of the \code{\link{repbase.query}} function to quantify
the number of hits for each query LTR transposon (duplicates) and retain
only hits found in Repbase that span the annotation sequence in Repbase
to a certain percentage (\code{scope}).
}
\examples{
\dontrun{
# PreProcess Repbase: A thaliana
# and save the output into the file "Athaliana_repbase.ref"
repbase.clean(repbase.file = "athrep.ref",
              output.file  = "Athaliana_repbase.ref")
             
# perform blastn search against A thaliana repbase annotation
AthalianaRepBaseAnnotation <- repbase.query(ltr.seqs     = "TAIR10_chr_all-ltrdigest_complete.fas", 
                                           repbase.path = "Athaliana_repbase.ref", 
                                           cores        = 1)
 # filter the annotation query output                                           
 AthalianaAnnot.HighMatches <- repbase.filter(AthalianaRepBaseAnnotation, 
                                              scope = 0.9)
 Ath.TE.Matches.Families <- sort(table(
                            unlist(lapply(stringr::str_split(
                            names(table(AthalianaAnnot.HighMatches$subject_id)),"_"),
                            function(x) paste0(x[2:3],collapse = ".")))),
                                        decreasing = TRUE)
 
 # visualize the hits found to have a scope of 90\%
 barplot(Ath.TE.Matches.Families,
        las       = 3, 
        cex.names = 0.8,
        col       = bcolor(length(Ath.TE.Matches.Families)), 
        main = "RepBase Annotation: A. thaliana")

}
}
\seealso{
\code{\link{repbase.query}}, \code{\link{repbase.clean}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/join.with.genome.tbl.R.R
\name{join.with.genome.tbl}
\alias{join.with.genome.tbl}
\title{Join \code{gm_files} returned by \code{generate.multi.quality.filter.meta} with a genome information table}
\usage{
join.with.genome.tbl(meta.summary.file, genome.tbl)
}
\arguments{
\item{meta.summary.file}{a meta.summary.file.}

\item{genome.tbl}{a genome.tbl.}
}
\description{
Join \code{gm_files} returned by \code{generate.multi.quality.filter.meta} with a genome information table.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred2GRanges.R
\name{pred2GRanges}
\alias{pred2GRanges}
\title{Format LTR prediction data to \code{GRages} object}
\usage{
pred2GRanges(LTR.data, similarity.threshold = 95)
}
\arguments{
\item{LTR.data}{the LTR prediction \code{\link{data.frame}} generated by \code{\link{LTRpred}}.}

\item{similarity.threshold}{the LTR similarity threshold that shall be used to define young and old retrotransposons.}
}
\description{
This function formats the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRpred}} to a \code{GRages} object.
}
\details{
The \code{GRages} object is defined by the following columns:
\itemize{
\item \code{chromosome}
\item \code{start}
\item \code{end}
\item \code{strand}
\item \code{element_name}
\item \code{ltr_similarity}
\item \code{orfs}
\item \code{age}
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cn2bed.R
\name{cn2bed}
\alias{cn2bed}
\title{Write copy number estimation results to BED file format.}
\usage{
cn2bed(
  cn.pred,
  type = "solo",
  filename = "copy_number_est",
  sep = "\\t",
  output = NULL
)
}
\arguments{
\item{cn.pred}{\code{data.frame} object returned by \code{\link{ltr.cn}} (if \code{type = "solo"}).}

\item{type}{type of copy number estimation: \code{\link{ltr.cn}} (if \code{type = "solo"}).}

\item{filename}{name of the output file (will be extended by "*.csv").}

\item{sep}{column separator.}

\item{output}{path in which the output file shall be stored.}
}
\description{
Format copy number estimation output to BED file format.
}
\seealso{
\code{\link{ltr.cn}}, \code{\link{LTRpred}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.multi.quality.filter.meta.R
\name{generate.multi.quality.filter.meta}
\alias{generate.multi.quality.filter.meta}
\title{Run \code{quality.filter.meta} for several different ltr similarity thresholds.}
\usage{
generate.multi.quality.filter.meta(
  kingdom,
  genome.folder,
  ltrpred.meta.folder,
  sim.options,
  cut.range.options,
  n.orfs = 0,
  strategy = "default",
  update = FALSE
)
}
\arguments{
\item{kingdom}{the taxonomic kingdom of the species for which \code{\link{LTRpred}} annotations are 
stored in the \code{genome.folder}.}

\item{genome.folder}{a file path to a folder storing the genome assembly files in fasta format that
were used to generate \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.}

\item{ltrpred.meta.folder}{a file path to a folder storing \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.}

\item{sim.options}{a numeric vector storing the ltr similarity thresholds that shall be probed.}

\item{cut.range.options}{a numeric vector storing the similarity cut range thresholds that shall be probed.}

\item{n.orfs}{minimum number of open reading frames a predicted retroelement shall possess.}

\item{strategy}{quality filter strategy. Options are
\itemize{
\item \code{strategy = "default"} : see section \code{Quality Control} 
\item \code{strategy = "stringent"} : in addition to filter criteria specified in section \code{Quality Control},
the filter criteria \code{!is.na(protein_domain)) | (dfam_target_name != "unknown")} is applied
}}

\item{update}{shall already existing \code{_SimilarityMatrix.csv} and \code{_GenomeInfo.csv} files be updated (\code{update = TRUE}) or can the already existing files be used (\code{update = FALSE})?}
}
\value{
A list with to list elements \code{sim_file} and \code{gm_file}. Each list element stores a \code{data.frame}:
  \itemize{
  \item \code{sim_file} (similarity file)
         \itemize{
         This \code{data.frame} stores the information 
                  \item 
                  }
   \item \code{gm_file} (genome metrics file)
         \itemize{
         This \code{data.frame} stores the information
                  \item
                  }
   }
}
\description{
A helper function to apply the \code{\link{quality.filter}} function to diverse \code{\link{LTRpred}} annotations while probing different ltr similarity thresholds.
}
\details{
\strong{Quality Control}

\itemize{
\item \code{ltr.similarity}: Minimum similarity between LTRs. All TEs not matching this
 criteria are discarded.
 \item \code{n.orfs}: minimum number of Open Reading Frames that must be found between the
  LTRs. All TEs not matching this criteria are discarded.
 \item \code{PBS or Protein Match}: elements must either have a predicted Primer Binding
 Site or a protein match of at least one protein (Gag, Pol, Rve, ...) between their LTRs. All TEs not matching this criteria are discarded.
 \item The relative number of N's (= nucleotide not known) in TE <= 0.1. The relative number of N's is computed as follows: absolute number of N's in TE / width of TE.
}
}
\seealso{
\code{\link{quality.filter}}, \code{\link{quality.filter.meta}}, \code{\link{LTRpred}}, \code{\link{LTRpred.meta}}, \code{\link{read.ltrpred}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality.filter.R
\name{quality.filter}
\alias{quality.filter}
\title{Pipeline to eliminate false positive predictions of retrotransposons}
\usage{
quality.filter(pred, sim, n.orfs, strategy = "default")
}
\arguments{
\item{pred}{\code{LTRpred.tbl} generated with \code{\link{LTRpred}}}

\item{sim}{LTR similarity threshold. Only putative LTR transposons that fulfill this 
LTR similarity threshold will be retained.}

\item{n.orfs}{minimum number of ORFs detected in the putative LTR transposon.}

\item{strategy}{quality filter strategy. Options are
\itemize{
\item \code{strategy = "default"} : see section \code{Quality Control} 
\item \code{strategy = "stringent"} : in addition to filter criteria specified in section \code{Quality Control},
the filter criteria \code{!is.na(protein_domain)) | (dfam_target_name != "unknown")} is applied
}}
}
\value{
A quality filtered \code{LTRpred.tbl}.
}
\description{
This function takes an \code{\link{LTRpred}} output table as input
and eliminates false positive predictions.
}
\details{
\strong{Quality Control}

\itemize{
\item \code{ltr.similarity}: Minimum similarity between LTRs. All TEs not matching this
 criteria are discarded.
 \item \code{n.orfs}: minimum number of Open Reading Frames that must be found between the
  LTRs. All TEs not matching this criteria are discarded.
 \item \code{PBS or Protein Match}: elements must either have a predicted Primer Binding
 Site or a protein match of at least one protein (Gag, Pol, Rve, ...) between their LTRs. All TEs not matching this criteria are discarded.
 \item The relative number of N's (= nucleotide not known) in TE <= 0.1. The relative number of N's is computed as follows: absolute number of N's in TE / width of TE.
}
}
\examples{
# example prediction file generated by LTRpred 
pred.file <- system.file("Athaliana_TAIR10_chr_all_LTRpred_DataSheet.csv", package = "LTRpred")
# read LTRpred generated prediction file (data sheet)
pred <- read.ltrpred(pred.file)
# apply quality filter
pred <- quality.filter(pred, sim = 70, n.orfs = 1)
}
\seealso{
\code{\link{LTRpred}}, \code{\link{LTRpred.meta}}, \code{\link{read.ltrpred}}, \code{\link{quality.filter.meta}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimMatAbundance.R
\name{SimMatAbundance}
\alias{SimMatAbundance}
\title{Compute histogram shape similarity between species}
\usage{
SimMatAbundance(sim.matrix)
}
\arguments{
\item{sim.matrix}{a species age distribution matrix generated with \code{\link{LTRpred.meta}}.}
}
\description{
For each pairwise species comparison the pattern similarity between
LTR age distributions is computed via \code{1 - cor(species1,species2)}.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.ltrpred.R
\name{read.ltrpred}
\alias{read.ltrpred}
\title{Import \code{LTRpred} DataSheet}
\usage{
read.ltrpred(data.sheet)
}
\arguments{
\item{data.sheet}{path to the \code{*_LTRpred_DataSheet.tsv} file.}
}
\description{
The \code{*_LTRpred_DataSheet.tsv} file generated by \code{\link{LTRpred}} stores the features of all predicted LTR transposons in a table. This function
imports this \code{\link{LTRpred}} output table.
}
\examples{
# example prediction file generated by LTRpred 
pred.file <- system.file("Hsapiens_ChrY_LTRpred_DataSheet.tsv", package = "LTRpred")
# read LTRpred generated prediction file (data sheet)
pred <- read.ltrpred(pred.file)

# or arrange by ltr_similarity
dplyr::arrange(tidy.datasheet(pred), dplyr::desc(ltr_similarity))
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.dfam.R
\name{read.dfam}
\alias{read.dfam}
\title{Import Dfam Query Output}
\usage{
read.dfam(dfam.file)
}
\arguments{
\item{dfam.file}{path to Dfam output file generated with \code{\link{dfam.query}}.}
}
\description{
Import the output file generated with \code{\link{dfam.query}}.
}
\examples{
# import example Dfam output
dfam.file <- read.dfam(system.file("example_dfam.out", package = "LTRpred"))

head(dfam.file)
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.tabout.R
\name{read.tabout}
\alias{read.tabout}
\title{Import information sheet returned by LTRdigest}
\usage{
read.tabout(tabout.file)
}
\arguments{
\item{tabout.file}{path to the tabout.csv file generated by the
\code{\link{LTRdigest}} function.}
}
\description{
This function imports the \code{*-ltrdigest_tabout.csv} file generated by \code{\link{LTRdigest}}.
}
\examples{
# example tabout file generated by LTRdigest for A thaliana
tabout.file <- system.file("TAIR10_chr_all-ltrdigest_tabout.csv",package = "LTRpred")
# import tabout file
imported.tabout <- read.tabout(tabout.file)
# look at the imported tabout file.
head(imported.tabout)      
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LTRharvest.R
\name{LTRharvest}
\alias{LTRharvest}
\title{Run LTRharvest to predict putative LTR Retrotransposons}
\usage{
LTRharvest(
  genome.file,
  index.file = NULL,
  range = c(0, 0),
  seed = 30,
  minlenltr = 100,
  maxlenltr = 3500,
  mindistltr = 4000,
  maxdistltr = 25000,
  similar = 70,
  mintsd = 4,
  maxtsd = 20,
  vic = 60,
  overlaps = "no",
  xdrop = 5,
  mat = 2,
  mis = -2,
  ins = -3,
  del = -3,
  motif = NULL,
  motifmis = 0,
  output.path = NULL,
  verbose = TRUE
)
}
\arguments{
\item{genome.file}{path to the genome file in \code{fasta} format.}

\item{index.file}{specify the name of the enhanced suffix array index file that is computed
by \code{suffixerator}. This opten can be used in case the suffix file was previously 
generated, e.g. during a previous call of this function. In this case the suffix array index
file does not need to be re-computed for new analyses. This is particularly useful when 
running \code{LTRharvest} with different parameter settings.}

\item{range}{define the genomic interval in which predicted LTR transposons shall be reported
. In case \code{range[1] = 1000} and \code{range[2] = 10000} then candidates are only 
reported if they start after position 1000 and end before position 10000 in their respective 
sequence coordinates. If \code{range[1] = 0} and \code{range[2] = 0}, 
so \code{range = c(0,0)} (default) then the entire genome is being scanned.}

\item{seed}{the minimum length for the exact maximal repeats. Only repeats with the specified minimum length are considered in all subsequent analyses. Default is \code{seed = 30}.}

\item{minlenltr}{minimum LTR length. Default is \code{minlenltr = 100}.}

\item{maxlenltr}{maximum LTR length. Default is \code{maxlenltr = 3500}.}

\item{mindistltr}{minimum distance of LTR starting positions. Default is \code{mindistltr = 4000}.}

\item{maxdistltr}{maximum distance of LTR starting positions. Default is \code{maxdistltr = 25000}.}

\item{similar}{minimum similarity value between the two LTRs in percent. \code{similar = 70}.}

\item{mintsd}{minimum target site duplications (TSDs) length. If no search for TSDs
shall be performed, then specify \code{mintsd = NULL}. Default is \code{mintsd = 4}.}

\item{maxtsd}{maximum target site duplications (TSDs) length. If no search for TSDs
shall be performed, then specify \code{maxtsd = NULL}. Default is \code{maxtsd = 20}.}

\item{vic}{number of nucleotide positions left and right (the vicinity) of the predicted
boundary of a LTR that will be searched for TSDs and/or one motif (if specified). 
Default is \code{vic = 60}.}

\item{overlaps}{specify how overlapping LTR retrotransposon predictions shall be treated. 
If \code{overlaps = "no"} is selected, then neither nested nor overlapping predictions will be reported in the output. In case \code{overlaps = "best"} is selected then in the case of two or more nested or overlapping predictions, solely the LTR retrotransposon prediction with
the highest similarity between its LTRs will be reported.
If \code{overlaps = "all"} is selected then all LTR retrotransposon predictions 
will be reported whether there are nested and/or overlapping predictions or not. 
Default is \code{overlaps = "best"}.}

\item{xdrop}{specify the xdrop value (> 0) for extending a seed repeat in both directions
allowing for matches, mismatches, insertions, and deletions. The xdrop extension process
 stops as soon as the extension involving matches, mismatches, insersions, and deletions 
 has a score smaller than T -X, where T denotes the largest score seen so far. Default is \code{cdrop = 5}.}

\item{mat}{specify the positive match score for the X-drop extension process. Default is \code{mat = 2}.}

\item{mis}{specify the negative mismatch score for the X-drop extension process. Default is \code{mis = -2}.}

\item{ins}{specify the negative insertion score for the X-drop extension process. Default is \code{ins = -3}.}

\item{del}{specify the negative deletion score for the X-drop extension process. Default is \code{del = -3}.}

\item{motif}{specify 2 nucleotides for the starting motif and 2 nucleotides for the ending
motif at the beginning and the ending of each LTR, respectively.
Only palindromic motif sequences - where the motif sequence is equal to its complementary
sequence read backwards - are allowed, e.g. \code{motif = "tgca"}. Type the nucleotides without any space
separating them. If this option is not selected by the user, candidate pairs will not be
screened for potential motifs. If this options is set but no allowed number of
mismatches is specified by the argument \code{motifmis} and a search for the exact 
motif will be conducted. If \code{motif = NULL} then no explicit motif is being specified.}

\item{motifmis}{allowed number of mismatches in the TSD motif specified in \code{motif}. 
The number of mismatches needs to be between [0,3].  Default is \code{motifmis = 0}.}

\item{output.path}{a path/folder to store all results returned by \code{LTRharvest}. 
If \code{output.path = NULL} (Default) then a folder with the name of the input genome file
will be generated in the current working directory of R and all results are then stored in this folder.}

\item{verbose}{logical value indicating whether or not detailed information shall be printed on the console.}
}
\value{
The \code{LTRharvest} function generates the following output files:

\itemize{
\item *_BetweenLTRSeqs.fsa : DNA sequences of the region between the LTRs in fasta format. 
\item *_Details.tsv : A spread sheet containing detailed information about the predicted LTRs.
\item *_FullLTRRetrotransposonSeqs.fsa : DNA sequences of the entire predicted LTR retrotransposon.
\item *_index.fsa : The suffixarray index file used to predict putative LTR retrotransposonswith \code{LTRharvest}.
\item *_Prediction.gff : A spread sheet containing detailed additional information about the predicted LTRs (partially redundant with the *_Details.tsv file).
}
The ' * ' is an place holder for the name of the input genome file.
}
\description{
This function implements an interface between R and
the LTRharvest command line tool to predict putative LTR retrotransposons from R.
}
\details{
The \code{LTRharvest} function provides an interface to the \code{LTRharvest} command line
tool and furthermore takes care of the entire folder handling, output parsing, and data
processing of the \code{LTRharvest} prediction.

Internally a folder named \code{output.path}_ltrharvest is generated and all computations
returned by \code{LTRharvest} are then stored in this folder. These files (see section \code{Value}) are then parsed and returned as list of data.frames by this function.

\code{LTRharvest} can be used as independently or as initial pre-computation step
to sufficiently detect LTR retrotransposons with \code{LTRdigest}.
}
\examples{
\dontrun{

# Run LTRharvest for H sapines partial Y chromosome using standard parameters
LTRharvest(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))
}
}
\references{
D Ellinghaus, S Kurtz and U Willhoeft. LTRharvest, an efficient and flexible software for de novo detection of LTR retrotransposons. BMC Bioinformatics (2008). 9:18.

Most argument specifications are adapted from the User manual of LTRharvest.
}
\seealso{
\code{\link{LTRdigest}},  \code{\link{LTRpred}}, 
\code{\link{read.prediction}}, \code{\link{read.seqs}},
\code{\link{pred2fasta}}, \code{\link{pred2gff}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combinePreds.R
\name{combinePreds}
\alias{combinePreds}
\title{Combine LTRpred prediction files}
\usage{
combinePreds(folder, quality.filter = FALSE, sim = 70, n.orfs = 0)
}
\arguments{
\item{folder}{path to folder in which prediction files are stored.}

\item{quality.filter}{shall false positives be filtered out as much as possible or not.}

\item{sim}{If \code{quality.filter = TRUE}: LTR similarity threshold. Only putative LTR transposons that fulfill this 
LTR similarity threshold will be retained.}

\item{n.orfs}{If \code{quality.filter = TRUE}: minimum number of ORFs detected in the putative LTR transposon.}
}
\description{
Given a path to a folder that stores several LTRpred prediction files,
this function imports the individual prediction files and combines them to one
large LTRpred prediction file.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.uc.R
\name{filter.uc}
\alias{filter.uc}
\title{Filter for cluster members}
\usage{
filter.uc(cluster.file)
}
\arguments{
\item{cluster.file}{a \code{data.frame} \code{*.uc} (USEARCH cluster) format (imported with \code{\link{read.uc}}).}
}
\description{
Filter for cluster memebers in a \code{*.uc} format table.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.jumpers.R
\name{filter.jumpers}
\alias{filter.jumpers}
\title{Detect LTR retrotransposons that are potential jumpers}
\usage{
filter.jumpers(LTRpred.tbl, ltr.similarity = 95, strategy = "conservative")
}
\arguments{
\item{LTRpred.tbl}{a \code{\link{data.frame}} generated by \code{\link{LTRpred}}.}

\item{ltr.similarity}{LTR similarity threshold. Default is \code{ltr_similarity = 95}.}

\item{strategy}{filter strategy: either \code{conservative}, \code{liberal}, or \code{between}.}
}
\description{
This function applies specific filter criteria to
screen for LTR retrotransposons predicted by \code{\link{LTRpred}}
that are potentially able to transpose due to their sequence features.
}
\details{
This ...

\strong{Filter strategy}
\itemize{
\item \code{conservative} :
\item \code{liberal} :
\item \code{between} :
}
}
\examples{
\dontrun{
# generate de novo LTR transposon prediction
pred <- LTRpred(genome.file = "TAIR10_chr_all.fas",
                trnas       = "plantRNA_Arabidopsis.fsa",
                hmms        = "hmm_*")
                
# detect potential jumpers               
filter.jumpers(pred)
}
}
\seealso{
\code{\link{LTRpred}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLTRSim.R
\name{plotLTRSim}
\alias{plotLTRSim}
\title{Plot the age distribution of predicted LTR transposons}
\usage{
plotLTRSim(
  data,
  type = "hist",
  stack.fill = "protein_domain",
  similarity.bin = 2,
  min.sim = 70,
  quality.filter = TRUE,
  n.orfs = 1,
  xlab = "LTR \% Similarity",
  ylab = "Frequency",
  main = "LTR Age Distribution",
  legend.title = "LTR Similarity"
)
}
\arguments{
\item{data}{the \code{\link{data.frame}} generated by \code{\link[LTRpred]{LTRpred}}.}

\item{type}{type of histogram. Either normal histogram (\code{type = "hist"}) or stacked histogram (\code{type = "stack"}, see also \code{stack.fill}).
If \code{type = "stack"} is specified then}

\item{stack.fill}{a character string specifying the variable by which the bar plot shall be stacked.}

\item{similarity.bin}{resolution of similarity binning. E.g. binning 98\%-100\% into 0.5\% intervals would be \code{similarity.bin = 0.5}.
Default is \code{similarity.bin = 0.5}.}

\item{min.sim}{minimum similarity between LTRs that can shall be considered for visualization. 
All elements not fulfilling this similarity threshold are filtered out. Default is \code{min.sim = 2}.}

\item{quality.filter}{shall false positives be filtered out as much as possible or not. See Description for details.}

\item{n.orfs}{minimum number of ORFs detected in the putative LTR transposon.}

\item{xlab}{x-axis label.}

\item{ylab}{y-axis label.}

\item{main}{main text.}

\item{legend.title}{legend text.}
}
\description{
This function visualizes the age distribution of
predicted LTR transposons generated with \code{\link[LTRpred]{LTRpred}}.

The age pf LTR transposons is defined by the sequence similarity between it's
3' and 5' LTR. Evolutionary young (recent) LTR transposons tend to have very similar
LTRs (up to 100\% sequence similarity), whereas evolutionary older LTR transposons
tend to have less similar LTRs.
}
\details{
This way of visualizing the age distribution of LTR transposons
allows users to examine the rate of recent transposition events in extant organisms. 

LTR similarity values are binned in intervals (per default 2.5\% intervals) which can be modified
by the \code{similarity.bin} argument.
}
\examples{
\dontrun{
# run LTRpred for A. thaliana
Ath.Pred <- LTRpred(genome.file = "TAIR10_chr_all.fas")
# visualize the age distribution of predicted  A. thaliana LTR transposons
PlotLTRAge(Ath.Pred)
}
}
\seealso{
\code{\link[LTRpred]{LTRpred}}, \code{\link[LTRpred]{LTRharvest}}, \code{\link[LTRpred]{LTRdigest}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ltr_age_estimation.R
\name{ltr_age_estimation}
\alias{ltr_age_estimation}
\title{Estimate retrotransposon insertion age in Mya based on 5 prime and 3 prime LTR sequence homology}
\usage{
ltr_age_estimation(
  pred,
  ltr_seqs_3_prime,
  ltr_seqs_5_prime,
  model = "K80",
  mutation_rate = 1.3 * 1e-07
)
}
\arguments{
\item{pred}{a prediction file generated with \code{\link{LTRpred}}.}

\item{ltr_seqs_3_prime}{file path to a fasta file storing the sequences of the respective 3 prime LTR (e.g. as annotatted by \code{\link{LTRpred}}).}

\item{ltr_seqs_5_prime}{file path to a fasta file storing the sequences of the respective 5 prime LTR (e.g. as annotatted by \code{\link{LTRpred}}).}

\item{model}{a model as specified in \code{\link[ape]{dist.dna}}: a character string specifying the evolutionary model to be used - must be one of
 \itemize{
\item  \code{K80} (the default)
\item \code{raw}
\item  \code{N}
\item  \code{TS}
\item  \code{TV}
\item  \code{JC69}
\item  \code{F81} 
\item \code{K81}
\item \code{F84}
\item \code{BH87}
\item \code{T92}
\item \code{TN93}
\item \code{GG95}
\item \code{logdet}
\item \code{paralin}
}}

\item{mutation_rate}{a mutation rate per site per year. For retrotransposons the default is \eqn{mutation_rate = 1.3 * 10E-8} (Wicker and Keller, 2007).}
}
\description{
This function implements diverse metrics to roughly estimate
the insertion age in Mya based on 5 prime and 3 prime LTR sequence homology.
}
\examples{
\dontrun{
ltr_pred <- LTRpred::read.ltrpred(
          system.file("Hsapiens_ChrY_LTRpred_DataSheet.tsv", 
          package = "LTRpred"))
# define file path to fasta file storing 3 prime LTR sequences
ltr_seqs_3_prime <- system.file("Hsapiens_ChrY-ltrdigest_3ltr.fas", package = "LTRpred")
ltr_seqs_5_prime <- system.file("Hsapiens_ChrY-ltrdigest_5ltr.fas", package = "LTRpred")
# estimate insertion age based on 3 prime and 5 prime LTR homology using the K80 model
Hsapiens_ltr_age <- LTRpred::ltr_age_estimation(ltr_pred, ltr_seqs_3_prime, ltr_seqs_5_prime)
# look at results
Hsapiens_ltr_age
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.summarize.R
\name{meta.summarize}
\alias{meta.summarize}
\title{Summarize (concatenate) all predictions of a \code{LTRpred.meta} run}
\usage{
meta.summarize(
  result.folder,
  ltr.similarity = 70,
  quality.filter = TRUE,
  n.orfs = 0,
  strategy = "default"
)
}
\arguments{
\item{result.folder}{path to meta result folder generated by \code{\link{LTRpred.meta}}.}

\item{ltr.similarity}{only count elements that have an LTR similarity >= this threshold.}

\item{quality.filter}{optimize search to remove potential false positives (e.g. duplicated genes, etc.). See \code{Details} for further information on the filter criteria.}

\item{n.orfs}{minimum number of Open Reading Frames that must be found between the LTRs (if \code{quality.filter = TRUE}). See \code{Details} for further information on quality control.}

\item{strategy}{quality filter strategy. Options are
\itemize{
\item \code{strategy = "default"} : see section \code{Quality Control} 
\item \code{strategy = "stringent"} : in addition to filter criteria specified in section \code{Quality Control},
the filter criteria \code{!is.na(protein_domain)) | (dfam_target_name != "unknown")} is applied
}}
}
\value{
a \code{LTRpred.tbl} storing the \code{\link{LTRpred}} prediction \code{data.frames} for all species in the meta result folder generated by \code{\link{LTRpred.meta}}.
}
\description{
Crawl through all genome predictions performed with \code{\link{LTRpred.meta}}
and concatenate the prediction files for each species in the meta result folder
generated by \code{\link{LTRpred.meta}} to a meta-species \code{data.frame}.
}
\details{
This function crawls through each genome stored in the meta result folder
generated by \code{\link{LTRpred.meta}} and performs the following procedures:

\itemize{
\item \strong{Step 1:} For each genome: Read the \code{*._LTRpred_DataSheet.csv} file generated by \code{\link{LTRpred}}.
\item \strong{Step 2:} For each genome: Perform quality filtering and selection of elements having at least \code{ltr.similarity} sequence similarity between their LTRs (if \code{quality.filter = TRUE}). Otherwise no quality filtering is performed.
\item \strong{Step 3:} Summarize all genome predictions in the meta-folder to one meta-species \code{data.frame}.
}

\strong{Quality Filtering}

The aim of the quality filtering step is to reduce the potential false positive
LTR transposons that were predicted by \code{\link{LTRpred}}. These false positives can be
duplicated genes, or other homologous repetitive elements that fulfill the LTR similarity
criteria, but do not have any Primer Binding Site, Open Reading Frames, Gag and Pol
proteins, etc. To reduce the number of false positives, the following filters are applied
to discard false positive LTR transposons.

\itemize{
\item \code{ltr.similarity}: Minimum similarity between LTRs. All TEs not matching this
 criteria are discarded.
 \item \code{n.orfs}: minimum number of Open Reading Frames that must be found between the
  LTRs. All TEs not matching this criteria are discarded.
 \item \code{PBS or Protein Match}: elements must either have a predicted Primer Binding
 Site or a protein match of at least one protein (Gag, Pol, Rve, ...) between their LTRs. All TEs not matching this criteria are discarded.
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant.space.meta.R
\name{quant.space.meta}
\alias{quant.space.meta}
\title{Quantify the genomic loci space for multiple genomes, proteomes or Repeat Masker annotation files}
\usage{
quant.space.meta(folder, type = "proteome")
}
\arguments{
\item{folder}{path to the folder storing the genomes, proteomes, or 
Repeat Masker files of interest.}

\item{type}{The following genomic features can be quantified:
\itemize{
\item \code{type = 'cds'}: coding sequence files in fasta format (see \code{\link[biomartr]{read_genome}}).
\item \code{type = 'proteome'}: proteome files in fasta format (see \code{\link[biomartr]{read_proteome}}).
\item \code{type = 'rm'}: Repeat Masker annotation files in fasta format (see \code{\link[biomartr]{read_rm}}).
}}
}
\description{
Quantification of the genomic loci space (= total length of
all annotated genomic features) within a given genome of interest.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.uc.R
\name{read.uc}
\alias{read.uc}
\title{Read file in USEARCH cluster format}
\usage{
read.uc(uc.file)
}
\arguments{
\item{uc.file}{path to file in USEARCH cluster format (\code{*.uc} file extension).}
}
\value{
A dataframe storing the following columns:

\itemize{
\item \code{Type:} Record type 'S', 'H', 'C', or 'N'.
\item \code{Cluster:} Cluster number (0-based).
\item \code{Size:} Sequence length ('S', 'N', and 'H') or cluster size 'C'.
\item \code{Perc_Ident:} For 'H' records, percent identity with target.
\item \code{Strand:} For 'H' records, the strand: '+' or '-' for nucleotides; '.' for proteins.
\item \code{Query:}  query id.
\item \code{Target:} target id.
}

Details:

Record type: 
\itemize{
\item \code{Type 'H' :} Hit. Represents an alignment between the query sequence and target sequence. For clustering 'H' indicates the cluster assignment for the query.
\item \code{Type 'S' :} Centroid (clustering only). There exists only one 'S' record
for each cluster, this gives the centroid (representative) sequence label in the \code{Query} column.
\item \code{Type 'C' :} Cluster record (clustering only). The \code{Size} column specifies the cluster size and the \code{Query} column the query id that corresponds to this cluster.
\item \code{Type 'N' :} No hit (for database search without clustering only). Indicates that no hit of the query were found in the target database. In the case of clustering, a query without hits becomes the centroid of a new cluster and generates an 'S'
record instead of an 'N' record.
}
}
\description{
Read a file in USEARCH cluster format generated by either USEARCH or VSEARCH.
}
\examples{
# read example *.uc file
test.uc <- read.uc(system.file("test.uc", package = "LTRpred"))

# look at the format in R
head(test.uc)
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repbase.clean.R
\name{repbase.clean}
\alias{repbase.clean}
\title{Clean the initial Repbase database for BLAST}
\usage{
repbase.clean(repbase.file, output.file)
}
\arguments{
\item{repbase.file}{fasta file storing the corresponding Repbase annotation, e.g. \code{athrep.ref}.}

\item{output.file}{name/path of the cleaned Repbase annotation file.}
}
\description{
Clean the headers of the Repbase fasta files
so that headers can be used to create a blast-able database.
}
\details{
The Repbase database can be downloaded after registration at http://www.girinst.org/repbase/.
The corresponding files as they are however, cannot be converted into a blast-able database.
Hence, a pre-filtering step is neccessary to be able to use this database with the e.g. \code{\link{repbase.query}}
function.
}
\examples{
\dontrun{
# PreProcess Repbase: A thaliana
# and save the output into the file "Athaliana_repbase.ref"
repbase.clean(repbase.file = "athrep.ref",
              output.file  = "Athaliana_repbase.ref")

}
}
\references{
http://www.girinst.org/repbase/
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.seqs.R
\name{read.seqs}
\alias{read.seqs}
\title{Import sequences of predicted LTR transposons}
\usage{
read.seqs(seq.file, program = "LTRharvest")
}
\arguments{
\item{seq.file}{path to fasta file storing the sequences of predicted LTR transposons generated by 
\code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.}

\item{program}{program used to generate the LTR transposons specified in \code{seq.file}, e.g. \code{program = "LTRpred"}, \code{program = "LTRdigest"}, or \code{program = "LTRharvest"}.}
}
\description{
Import sequences of predicted LTR transposons generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality.filter.meta.R
\name{quality.filter.meta}
\alias{quality.filter.meta}
\title{Pipeline to eliminate false positive predictions of retrotransposons at a metagenomic scale}
\usage{
quality.filter.meta(
  kingdom,
  genome.folder,
  ltrpred.meta.folder,
  sim,
  cut.range = 2,
  n.orfs,
  strategy,
  update = FALSE
)
}
\arguments{
\item{kingdom}{a character string specifying the kingdom of life to which genomes annotated with \code{\link{LTRpred.meta}} belong to. E.g. \code{kingdom = "Plants"}. If annotates of a variety of kingdoms have been done users can for example specify \code{kingdom = "Various"}.}

\item{genome.folder}{path to folder storing the genome assembly files that were used for \code{\link{LTRpred.meta}} predictions.}

\item{ltrpred.meta.folder}{path to folder storing the \code{\link{LTRpred.meta}} output files.}

\item{sim}{LTR similarity threshold. Only putative LTR transposons that fulfill this 
LTR similarity threshold will be retained.}

\item{cut.range}{a numeric number indicating the interval size for binning LTR similarities.}

\item{n.orfs}{minimum number of ORFs detected in the putative LTR transposon.}

\item{strategy}{quality filter strategy. Options are
\itemize{
\item \code{strategy = "default"} : see section \code{Quality Control} 
\item \code{strategy = "stringent"} : in addition to filter criteria specified in section \code{Quality Control},
the filter criteria \code{!is.na(protein_domain)) | (dfam_target_name != "unknown")} is applied
}}

\item{update}{shall already existing \code{_SimilarityMatrix.csv} and \code{_GenomeInfo.csv} files be updated (\code{update = TRUE}) or can the already existing files be used (\code{update = FALSE})?}
}
\value{
A list with to list elements \code{sim_file} and \code{gm_file}. Each list element stores a \code{data.frame}:
  \itemize{
  \item \code{sim_file} (similarity file)
         \itemize{
         This \code{data.frame} stores the information 
                  \item 
                  }
   \item \code{gm_file} (genome metrics file)
         \itemize{
         This \code{data.frame} stores the information
                  \item
                  }
   }
}
\description{
This function takes the file paths to the genomes folder and \code{\link{LTRpred.meta}} output folder as input and eliminates false positive retrotransposon predictions on a metagenomic scale.
}
\details{
\strong{Quality Control}

\itemize{
\item \code{ltr.similarity}: Minimum similarity between LTRs. All TEs not matching this
 criteria are discarded.
 \item \code{n.orfs}: minimum number of Open Reading Frames that must be found between the
  LTRs. All TEs not matching this criteria are discarded.
 \item \code{PBS or Protein Match}: elements must either have a predicted Primer Binding
 Site or a protein match of at least one protein (Gag, Pol, Rve, ...) between their LTRs. All TEs not matching this criteria are discarded.
 \item The relative number of N's (= nucleotide not known) in TE <= 0.1. The relative number of N's is computed as follows: absolute number of N's in TE / width of TE.
}
}
\seealso{
\code{\link{LTRpred}}, \code{\link{LTRpred.meta}}, \code{\link{read.ltrpred}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant.cds.space.R
\name{quant.cds.space}
\alias{quant.cds.space}
\title{Quantify the coding sequence space within a genome}
\usage{
quant.cds.space(file)
}
\arguments{
\item{file}{file path to a fasta file storing the cds sequences.}
}
\description{
Quantification of the cds space (= total length of
all annotated coding sequences) within a given genome of interest.
}
\seealso{
\code{\link{quant.protein.space}}, \code{\link{quant.repeat.space}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.apply.R
\name{meta.apply}
\alias{meta.apply}
\title{Apply functions to meta data generated by \code{LTRpred}}
\usage{
meta.apply(result.folder, FUN)
}
\arguments{
\item{result.folder}{\code{result.folder} generated by \code{\link{LTRpred.meta}}.}

\item{FUN}{function to apply to meta data returned by \code{\link{LTRpred.meta}}.}
}
\description{
Apply any function to the results generated by meta analyses
}
\seealso{
\code{\link{LTRpred.meta}}, \code{\link{LTRpred}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome.summary.R
\name{genome.summary}
\alias{genome.summary}
\title{Generating genome summary files for \code{LTRpred.meta} results}
\usage{
genome.summary(
  genome.folder,
  ltrpred.meta.folder,
  file.name,
  sim = 70,
  cut.range = 2,
  quality.filter = TRUE,
  n.orfs = 0,
  strategy = "default"
)
}
\arguments{
\item{genome.folder}{a file path to a folder storing the genome assembly files in fasta format that
were used to generate \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.}

\item{ltrpred.meta.folder}{a file path to a folder storing \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.}

\item{file.name}{name of the output file.}

\item{sim}{LTR similarity threshold. Only putative LTR transposons that fulfill this 
LTR similarity threshold will be retained.}

\item{cut.range}{similarity interval size.}

\item{quality.filter}{shall a quality filter to remove possible false positive predictions be applied?}

\item{n.orfs}{minimum number of open reading frames a predicted retroelement shall possess.}

\item{strategy}{quality filter strategy. Options are
\itemize{
\item \code{strategy = "default"} : see section \code{Quality Control} 
\item \code{strategy = "stringent"} : in addition to filter criteria specified in section \code{Quality Control},
the filter criteria \code{!is.na(protein_domain)) | (dfam_target_name != "unknown")} is applied
}}
}
\description{
Generating genome summary files for \code{LTRpred.meta} results
}
\details{
Generating genome summary files for \code{LTRpred.meta} results.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file.move.R
\name{file.move}
\alias{file.move}
\title{Move folders from one location to another}
\usage{
file.move(from, to)
}
\arguments{
\item{from}{path from where to copy a folder.}

\item{to}{path where to copy the corresponding folder specified in \code{from}.}
}
\description{
A small helper function to make file handling
in R easier.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/totalRepeatSpace.R
\name{totalRepeatSpace}
\alias{totalRepeatSpace}
\title{Quantify the total repeat space from Repeat Masker output in Mbp}
\usage{
totalRepeatSpace(file, repeat.type = "all")
}
\arguments{
\item{file}{path to Repeat Masker output file.}

\item{repeat.type}{type of element for which total repeat space shall be quantified.
Options are:

\itemize{
 \item \code{repeat.type = "all"} include all types of repeats.
 \item \code{repeat.type = "ltr"} include only LTR retrotransposons.
}}
}
\description{
Given a \code{Repeat Masker} annotation file this function quantifies
the total repeat space in mega base pairs (Mbp).
}
\seealso{
\code{\link{meta.seq.space}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfam.query.R
\name{dfam.query}
\alias{dfam.query}
\title{Query the Dfam database to annotate putative LTRs}
\usage{
dfam.query(
  seq.file,
  Dfam.db = NULL,
  eval = 0.001,
  cores = 1,
  output.folder = getwd()
)
}
\arguments{
\item{seq.file}{file path to the putative LTR transposon sequences in \code{fasta} format.}

\item{Dfam.db}{folder path to the local Dfam database or \code{Dfam.db = "download"} in case the Dfam
database shall be automatically downloaded before performing query analyses.}

\item{eval}{E-value threshhold to perform the HMMer search against the Dfam database.}

\item{cores}{number of cores to use to perform parallel computations.}

\item{output.folder}{folder path to store the annotation results.}
}
\description{
Validate or annotate putative LTR
transposons that have been predicted using \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.
}
\details{
The Dfam database provides a collection of curated transposable element annotations.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gridPlotAssemblyVersions.R
\name{gridPlotAssemblyVersions}
\alias{gridPlotAssemblyVersions}
\title{Plot retrotransposon age distributions for predictions of different genome assembly versions}
\usage{
gridPlotAssemblyVersions(
  data,
  quality.filter = FALSE,
  text.size = 18,
  xlab = "LTR similarity",
  ylab = "LTR retrotransposon count",
  main.text = "",
  y.ticks = 6,
  sim = 70,
  n.orfs = 0
)
}
\arguments{
\item{data}{prediction files returned by \code{\link[LTRpred]{LTRpred}}. Usually a combined prediction file that was generated with \code{\link{combinePreds}}.}

\item{quality.filter}{shall false positives be filtered out as much as possible or not.}

\item{text.size}{size of x-axis, y-axis, and title text.}

\item{xlab}{x-axis label.}

\item{ylab}{y-axis label.}

\item{main.text}{title text.}

\item{y.ticks}{number of ticks that shall be drawn on the y-axis.}

\item{sim}{If \code{quality.filter = TRUE}: LTR similarity threshold. Only putative LTR transposons that fulfill this 
LTR similarity threshold will be retained.}

\item{n.orfs}{If \code{quality.filter = TRUE}: minimum number of ORFs detected in the putative LTR transposon.}
}
\description{
Compare the abundance of retrotransposon predictions between different genome assembly versions. This way a genome assembly version bias can be either ruled out
or considered for subsequent conclusions.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ORFpred.R
\name{ORFpred}
\alias{ORFpred}
\title{Open Reading Frame Prediction}
\usage{
ORFpred(
  seq.file,
  orf.style = 7,
  min.codons = 200,
  trans.seqs = FALSE,
  output = NULL
)
}
\arguments{
\item{seq.file}{a fasta file storing the sequences for which open reading frames shall be predicted.}

\item{orf.style}{type of predicting open reading frames (see documentation of USEARCH).}

\item{min.codons}{minimum number of codons in the predicted open reading frame.}

\item{trans.seqs}{logical value indicating wheter or not predicted open reading frames
shall be translated and the corresponding protein sequences stored in the output folder.}

\item{output}{path to the folder in which predicted open reading frame sequences shall be stored.}
}
\description{
This function provides a wrapper to the \code{USEARCH}
\code{fastx_findorfs} command line tool to predict ORFs in
a given input fasta file and read the output as \code{\link{data.frame}} object.
}
\references{
Robert Edgar. Search and clustering orders of magnitude faster than BLAST. Bioinformatics (2010) 26 (19): 2460-2461.
}
\seealso{
\code{\link{LTRharvest}}, \code{\link{LTRdigest}}, \code{\link{LTRpred}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motif.count.R
\name{motif.count}
\alias{motif.count}
\title{Low level function to detect motifs in strings}
\usage{
motif.count(seq.file, motif, as.ratio = FALSE)
}
\arguments{
\item{seq.file}{path to the genomic sequecne file of interest (e.g. LTR TE seqs predicted by
\code{\link{LTRpred}}).}

\item{motif}{a character string or vector of strings which shall be counted within each sequence.}

\item{as.ratio}{shall count values be returned as asbolute frequency (count value) or as relative frequency (percentage).}
}
\description{
Find a specific motif or a sequence of motifs within
genomic sequences.
}
\examples{
# find number of "CG" motifs in predicted LTR transposons
motif.count(seq.file = system.file("LTRseqs.fas",package = "LTRpred"), 
            motif    = "CG")
            
# find number of "CG" motifs in predicted LTR transposons: rel. frequency
motif.count(seq.file = system.file("LTRseqs.fas",package = "LTRpred"), 
            motif    = "CG",
            as.ratio = TRUE)             
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred2csv.R
\name{pred2csv}
\alias{pred2csv}
\title{Format LTR prediction data to CSV file format}
\usage{
pred2csv(LTR.data, output = "output.csv")
}
\arguments{
\item{LTR.data}{the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.}

\item{output}{filename of the output CSV file.}
}
\value{
the \code{LTR.data} table in csv format saved to the hard drive.
}
\description{
This function formats the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}
to a \code{\link{data.frame}} in \code{CSV} file format.
}
\examples{
gff.file <- system.file("TAIR10_chr_all_LTRdigestPrediction.gff",
                        package = "LTRpred")
tabout.file <- system.file("TAIR10_chr_all-ltrdigest_tabout.csv"
                           ,package = "LTRpred")
LTRfile <- read.prediction(gff.file,tabout.file, program = "LTRdigest")

# generate csv file
pred2csv(LTRfile$ltr.retrotransposon)
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.members.R
\name{cluster.members}
\alias{cluster.members}
\title{Select members of a specific cluster}
\usage{
cluster.members(LTRpred.tbl, cluster = NULL)
}
\arguments{
\item{LTRpred.tbl}{\code{data.frame} generated by \code{\link{LTRpred}}.}

\item{cluster}{name of the cluster that shall be searched for.}
}
\description{
This function allows users to select all members of a specific cluster.
}
\examples{
\dontrun{ 
# example prediction file generated by LTRpred 
pred.file <- system.file("Athaliana_TAIR10_chr_all_LTRpred_DataSheet.csv", package = "LTRpred")
# read LTRpred generated prediction file (data sheet)
pred <- read.ltrpred(pred.file)

# extract memebers of cluster 'cl_108'
cluster.members(pred, cluster = "cl_108")

# or safe the sequences of the cluster as fasta
pred2fasta(cluster.members(pred, cluster = "cl_108"), "cl_108_seqs.fa")
}
}
\seealso{
\code{\link{clust2fasta}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename.fasta.R
\name{rename.fasta}
\alias{rename.fasta}
\title{Helper function to add species names to headers within fasta files}
\usage{
rename.fasta(file, species, output = "renamed_fasta.fa", append = FALSE)
}
\arguments{
\item{file}{input fasta file.}

\item{species}{a character string specifying the species name to be added to the headers. 
Format will be: \code{species_}*, where * stands for the original header.}

\item{output}{a character string denoting the name of the renamed output fasta file.}

\item{append}{logical value. If \code{TRUE} output will be appended to file; otherwise, it will overwrite the contents of file.}
}
\value{
Writes a new fasta file with renamed headers.
}
\description{
Reads a fasta file and adds the defined \code{species} name to each
header entry of the input fasta file.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred2orfseqs.R
\name{pred2orfseqs}
\alias{pred2orfseqs}
\title{Save the ORF sequence of the predicted LTR Transposons in a fasta file}
\usage{
pred2orfseqs(LTRpred.tbl, orf.seq.file, output = "output.fa")
}
\arguments{
\item{LTRpred.tbl}{the \code{\link{data.frame}} generated by \code{\link{LTRpred}}.}

\item{orf.seq.file}{the fasta file storing the open reading frame sequences of predicted retroelements.
as returned by the \code{\link{LTRpred}} function.}

\item{output}{the fasta file to which the output sequences shall be stored in.}
}
\description{
This function allows users to save the sequence of the predicted LTR Transposons or LTRs in a fasta file.
}
\details{
The output \code{data.frame}s returned by \code{\link{LTRpred}} contain all information of the predicted
LTR retrotransposons that can be used for post-filtering steps. After these post-filtering steps
sequences of the remaining (filtered) candidates can be retrieved by this function.
}
\seealso{
\code{\link{LTRharvest}}, \code{\link{LTRdigest}}, \code{\link{LTRpred}}, \code{\link{read.prediction}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign.short.sci.name.R
\name{assign.short.sci.name}
\alias{assign.short.sci.name}
\title{Transform long file names to short scientific names}
\usage{
assign.short.sci.name(data)
}
\arguments{
\item{data}{a \code{tibble} from \code{\link{meta.seq.space}}.}
}
\description{
This function takes a \code{tibble} as input
and transform long file names to short scientific names.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.blast6out.R
\name{read.blast6out}
\alias{read.blast6out}
\title{Read file in blast6out format generated by USEARCH or VSEARCH}
\usage{
read.blast6out(blast6out.file)
}
\arguments{
\item{blast6out.file}{to blast6out file (\code{*.blast6out} extension).}
}
\value{
A dataframe storing the following columns:

\itemize{
\item \code{query:} query id.
\item \code{subject:} subject id.
\item \code{perc_ident:} pecent identity between query and subject.
\item \code{align_len:} alignment length between query and subject.
\item \code{n_mismatch:} number of mismathces between query and subject.
\item \code{n_gap_open:} number of gap openings between query and subject. 
\item \code{start_q:} start position in query. Query coordinates start with 1 at the
 first base in the sequence as it appears in the input file. For translated searches (nucleotide queries, protein targets), query start < end for +ve frame and start > end for -ve frame.
\item \code{end_q:} end position in query.
\item \code{start_s:} start position in subject. Subject coordinates start with 1 at
 the first base in the sequence as it appears in the database. For untranslated
 nucleotide searches, subject start < end for plus strand, start > end for a reverse-complement alignment.
\item \code{end_s:} end position in subject.
\item \code{evalue:} evalue calculated using Karlin-Altschul statistics.
\item \code{bit_score:} bit score calculated using Karlin-Altschul statistics.
}
}
\description{
Read a file in blast6out format generated by either USEARCH or VSEARCH.
}
\examples{
# read example *.blast6out file
test.blast6out <- read.blast6out(system.file("test.blast6out", package = "LTRpred"))

# look at the format in R
head(test.blast6out)
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.seqs.R
\name{get.seqs}
\alias{get.seqs}
\title{Quickly retrieve the sequences of a \code{Biostrings} object}
\usage{
get.seqs(object, type = "DNA")
}
\arguments{
\item{object}{a \code{Biostrings} object.}

\item{type}{either \code{type = "DNA"} or \code{type = "AA"}.}
}
\description{
Helper function to retrieve the sequences of a \code{Biostrings} object.
}
\examples{
# read example sequences
seqs <- system.file("nt.fa",package = "LTRpred")
input.seqs <- Biostrings::readDNAStringSet(seqs)

# retrieve sequences of Biostrings object
head(get.seqs(input.seqs))
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LTRpred.meta.R
\name{LTRpred.meta}
\alias{LTRpred.meta}
\title{Perform Meta-Analyses with LTRpred}
\usage{
LTRpred.meta(genome.folder, output.folder, cores = 1, ...)
}
\arguments{
\item{genome.folder}{file path to folder storing genome assembly files in \code{fasta} format.}

\item{output.folder}{path to the output folder storing \code{LTRpred.meta} results that will be generated.}

\item{cores}{number of cores that shall be used for parallel processing.}

\item{\dots}{arguments that shall be passed to \code{\link{LTRpred}}.}
}
\description{
Run \code{\link{LTRpred}} on several genomes (sequentially) that are stored in a given folder.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CLUSTpred.R
\name{CLUSTpred}
\alias{CLUSTpred}
\title{Cluster Sequences with VSEARCH}
\usage{
CLUSTpred(
  file,
  similarity = 0.9,
  strand = "both",
  cores = 1,
  out.name = paste0(basename(file), "_CLUSTpred"),
  output = NULL
)
}
\arguments{
\item{file}{path to predicted LTR transposons generated by \code{\link{LTRpred}} (in fasta format).}

\item{similarity}{reject if sequence similarity is lower than this threshold.}

\item{strand}{cluster using plus or both strands.}

\item{cores}{number of cores that shall be used for parallel computations.}

\item{out.name}{name of the output files (\code{*.uc}, \code{*.log}, \code{*.blast6out}).}

\item{output}{path to a folder in which output shall be stored.}
}
\value{
First the following files generated by VSEARCH are stored in the output folder (default: \code{\link{getwd}}):

\itemize{
\item \code{*.uc} USEARCH cluster format generated by VSEARCH storing the sequence cluster information of the input LTR transposons.
\item \code{*.log} a log file of the VSEARCH run.
\item \code{*.blast6out} BLAST output generated by VSEARCH storing the BLAST hit information of the input LTR transposons.
}

  
A USEARCH cluster format (*.uc file extension) table (see \code{\link{read.uc}} for specifications).
}
\description{
Cluster putative LTR transposons predicted by \code{\link{LTRpred}} using VSEARCH.
}
\details{
To be able to use this function the VSEARCH command line tool needs to be installed.
}
\references{
https://github.com/torognes/vsearch
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant.protein.space.R
\name{quant.protein.space}
\alias{quant.protein.space}
\title{Quantify the protein space within a genome}
\usage{
quant.protein.space(file)
}
\arguments{
\item{file}{file path to a fasta file storing the protein sequences.}
}
\description{
Quantification of the protein space (= total length of
all annotated protein coding genes) within a given genome of interest.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_element_distr_along_chromosome.R
\name{plot_element_distr_along_chromosome}
\alias{plot_element_distr_along_chromosome}
\title{Plot positions of predicted retrotransposons along chromosomes}
\usage{
plot_element_distr_along_chromosome(
  pred,
  genome.file,
  centromere_start = NULL,
  centromere_end = NULL,
  ...
)
}
\arguments{
\item{pred}{LTRpred.tbl generated with \code{\link{LTRpred}}.}

\item{genome.file}{a file path to the genome assembly file in fasta format for which chromosomes shall be visualized.}

\item{centromere_start}{a numeric vector storing the centromere start coordinates in the \code{genome.file}. The position in the numeric vector should correspond to the chromosome name in the \code{genome.file} fasta file. If \code{centromere_start = NULL} (default), then no centromeres will be drawn.}

\item{centromere_end}{a numeric vector storing the centromere end coordinates in the \code{genome.file}. The position in the numeric vector should correspond to the chromosome name in the \code{genome.file} fasta file. If \code{centromere_end = NULL} (default), then no centromeres will be drawn.}

\item{...}{additional arguments that shall be passed to the visualization function \code{\link[ggbio]{autoplot}}.}
}
\description{
The positionas of LTR retrotransposons predicted with \code{\link{LTRpred}}
will be visualized along the chromsome.
}
\examples{
\dontrun{
test_genome <- system.file("Hsapiens_ChrY.fa", package = "LTRpred")
test_pred <- LTRpred::read.ltrpred(
system.file("Hsapiens_ChrY_LTRpred_DataSheet.tsv", 
             package = "LTRpred"))
test_centromere_starts <- 55000
# generate visualization
LTRpred::plot_element_distr_along_chromosome(test_pred, test_genome, test_centromere_starts)
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant.repeat.space.R
\name{quant.repeat.space}
\alias{quant.repeat.space}
\title{Quantify the repeat space within a genome}
\usage{
quant.repeat.space(file)
}
\arguments{
\item{file}{file path to a fasta file storing the Repeat Masker annotation file.}
}
\description{
Quantification of the repeat space (= total length of
all Repeat Masker annotated repeats) within a given genome of interest.
}
\seealso{
\code{\link{quant.protein.space}}, \code{\link{quant.cds.space}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.pred.filenames.R
\name{get.pred.filenames}
\alias{get.pred.filenames}
\title{Retrieve file names of files genereated by LTRpred}
\usage{
get.pred.filenames(LTRpred.folder, ws.wrap = FALSE)
}
\arguments{
\item{LTRpred.folder}{path to a \code{\link{LTRpred}} generated folder.}

\item{ws.wrap}{shall white space separated path names be wrapped with \code{'file name'}.}
}
\description{
Useful helper function to generate automated file names
when crawling LTRpred generated result folders.
}
\examples{
# without wrapped whitespaces
get.pred.filenames("path/to/LTRpred.folder")

# with wrapped whitespaces
get.pred.filenames("path to/LTRpred.folder", ws.wrap = TRUE )
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.prediction.R
\name{read.prediction}
\alias{read.prediction}
\title{Import the output of LTRharvest or LTRdigest}
\usage{
read.prediction(
  gff.file = NULL,
  tabout.file,
  program = "LTRdigest",
  ltr.fasta = NULL,
  inner.seq.fasta = NULL,
  data = NULL,
  similarity.bin = 2,
  min.sim = NULL
)
}
\arguments{
\item{gff.file}{path to GFF file generated by either \code{\link{LTRharvest}}, \code{\link{LTRdigest}},
or \code{\link{LTRpred}}.}

\item{tabout.file}{path to \code{*_tabout.csv} file generated by either \code{\link{LTRdigest}} or
or \code{\link{LTRpred}}.}

\item{program}{prediction program used to generate the corresponding GFF file. 
Either \code{program = "LTRharvest"}, \code{program = "LTRdigest"}, or \code{program = "LTRpred"}.}

\item{ltr.fasta}{fasta file generated by either \code{\link{LTRdigest}} or
or \code{\link{LTRpred}} storing the sequences of the predicted LTR retrotransposons (e.g. \code{*-ltrdigest_complete.fas}).}

\item{inner.seq.fasta}{path to the \code{*_BetweenLTRSeqs.fsa} fasta file returned by \code{\link{LTRharvest}} in caste \code{program = "LTRharvest"} is specified.}

\item{data}{the \code{\link{data.frame}} returned by either \code{\link{LTRharvest}}, \code{\link{LTRdigest}},
or \code{\link{LTRpred}} that can be used to modify for example \code{similarity.bin} and \code{min.sim} within this table.}

\item{similarity.bin}{resolution of similarity binning. E.g. binning 98\%-100\% into 0.5\% intervals would be \code{similarity.bin = 0.5}.}

\item{min.sim}{minimum similarity between LTRs that can shall be considered for visualization.}
}
\description{
This function imports the output files generated by \code{\link{LTRharvest}} or \code{\link{LTRdigest}}.
}
\examples{
\dontrun{
gff.file <- system.file("TAIR10_chr_all_LTRdigestPrediction.gff", package = "LTRpred")
tabout.file <- system.file("TAIR10_chr_all-ltrdigest_tabout.csv" ,package = "LTRpred")
LTRfile <- read.prediction(gff.file,tabout.file, program = "LTRdigest")

# look at results
str(LTRfile$ltr.retrotransposon)

}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clust2fasta.R
\name{clust2fasta}
\alias{clust2fasta}
\title{Export sequences of TEs belonging to the same cluster}
\usage{
clust2fasta(LTRpred.tbl, seqs)
}
\arguments{
\item{LTRpred.tbl}{\code{data.frame} generated by \code{\link{LTRpred}}.}

\item{seqs}{path to the fasta file storing the TE sequences.}
}
\description{
Export sequences of TEs belonging to the same cluster to separate fasta files.
Each cluster will be stored in an individual fasta file.
}
\seealso{
\code{\link{LTRpred}}, \code{\link{cluster.members}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LTRpred.R
\name{LTRpred}
\alias{LTRpred}
\title{Predict LTR retrotransposons in a given genome}
\usage{
LTRpred(
  genome.file = NULL,
  centromere_start = NULL,
  centromere_end = NULL,
  ltr_age_estimation = TRUE,
  model = "K80",
  mutation_rate = 1.3 * 1e-07,
  index.file.harvest = NULL,
  index.file.digest = NULL,
  LTRdigest.gff = NULL,
  tabout.file = NULL,
  LTRharvest.folder = NULL,
  LTRpred.folder = NULL,
  orf.file = NULL,
  annotate = NULL,
  Dfam.db = NULL,
  dfam.eval = 0.001,
  dfam.file = NULL,
  cluster = FALSE,
  clust.sim = 0.9,
  clust.file = NULL,
  copy.number.est = FALSE,
  fix.chr.name = FALSE,
  cn.eval = 1e-05,
  range = c(0, 0),
  seed = 30,
  minlenltr = 100,
  maxlenltr = 3500,
  mindistltr = 4000,
  maxdistltr = 25000,
  similar = 70,
  mintsd = 4,
  maxtsd = 20,
  vic = 60,
  overlaps = "no",
  xdrop = 5,
  mat = 2,
  mis = -2,
  ins = -3,
  del = -3,
  motif = NULL,
  motifmis = 0,
  aaout = "yes",
  aliout = "yes",
  pptlen = c(8, 30),
  uboxlen = c(3, 30),
  pptradius = 30,
  trnas = NULL,
  pbsalilen = c(11, 30),
  pbsoffset = c(0, 5),
  pbstrnaoffset = c(0, 5),
  pbsmaxedist = 1,
  pbsradius = 30,
  hmms = NULL,
  pdomevalcutoff = 1e-05,
  pbsmatchscore = 5,
  pbsmismatchscore = -10,
  pbsinsertionscore = -20,
  pbsdeletionscore = -20,
  pfam.ids = NULL,
  cores = 1,
  dfam.cores = NULL,
  hmm.cores = NULL,
  orf.style = 7,
  min.codons = 200,
  trans.seqs = FALSE,
  output.path = NULL,
  quality.filter = TRUE,
  n.orfs = 0,
  job_num = 1,
  verbose = TRUE
)
}
\arguments{
\item{genome.file}{path to the genome file in \code{fasta} format.}

\item{centromere_start}{a numeric vector storing the centromere start coordinates in the \code{genome.file}. The position in the numeric vector should correspond to the chromosome name in the \code{genome.file} fasta file. If \code{centromere_start = NULL} (default), then no centromeres will be drawn.}

\item{centromere_end}{a numeric vector storing the centromere end coordinates in the \code{genome.file}. The position in the numeric vector should correspond to the chromosome name in the \code{genome.file} fasta file. If \code{centromere_end = NULL} (default), then no centromeres will be drawn.}

\item{ltr_age_estimation}{a logical value indicating wether or not ltr-age estimation shall be performed (default \code{ltr_age_estimation = TRUE}).}

\item{model}{a model as specified in \code{\link[ape]{dist.dna}}: a character string specifying the evolutionary model to be used - must be one of
 \itemize{
\item  \code{K80} (the default)
\item \code{raw}
\item  \code{N}
\item  \code{TS}
\item  \code{TV}
\item  \code{JC69}
\item  \code{F81} 
\item \code{K81}
\item \code{F84}
\item \code{BH87}
\item \code{T92}
\item \code{TN93}
\item \code{GG95}
\item \code{logdet}
\item \code{paralin}
}}

\item{mutation_rate}{a mutation rate per site per year. For retrotransposons the default is \eqn{mutation_rate = 1.3 * 10E-8} (Wicker and Keller, 2007).}

\item{index.file.harvest}{specify the name of the enhanced suffix array index file that is computed
by \code{suffixerator} for the use of \code{LTRharvest}. This often can be used in case the suffix file was previously 
generated, e.g. during a previous call of this function. In this case the suffix array index
file does not need to be re-computed for new analyses. This is particularly useful when 
running \code{LTRpred} with different parameter settings. A possible index file name could be \code{BASENAME_index.fsa}.}

\item{index.file.digest}{specify the name of the enhanced suffix array index file that is computed
by \code{suffixerator} for the use of \code{LTRdigest}. This often can be used in case the suffix file was previously 
generated, e.g. during a previous call of this function. In this case the suffix array index
file does not need to be re-computed for new analyses. This is particularly useful when 
running \code{LTRpred} with different parameter settings. A possible index file name could be \code{BASENAME_index_ltrdigest.fsa}.}

\item{LTRdigest.gff}{path to the LTRdigest generated GFF file, in case LTRdigest files were pre-computed previously.}

\item{tabout.file}{path to the LTRdigest generated tabout file file, in case LTRdigest files were pre-computed previously.}

\item{LTRharvest.folder}{either \code{LTRharvest.folder = NULL} (default) or \code{LTRharvest.folder = "skip"} to skip moving a LTRharvest folder if it is not present.}

\item{LTRpred.folder}{name/path of/to an existing LTRpred folder.}

\item{orf.file}{path to the file generated by \code{\link{ORFpred}}, in case the orf prediction file was generated previously.}

\item{annotate}{annotation database that shall be queried to annotate predicted LTR transposons.
Default is \code{annotate = NULL} indicating that no annotation query is being performed.
Possible options are: \code{annotate = "Dfam"} (here the Dfam database must be stored locally and a nhammer search is performed against the Dfam database) or \code{annotate = "Repbase"} (here the Repbase database must be stored locally and a blastn search is performed against the Repbase database). Please consult the vignettes for details.}

\item{Dfam.db}{folder path to the local Dfam database or \code{Dfam.db = "download"} in case the Dfam
database shall be automatically downloaded before performing query analyses.}

\item{dfam.eval}{E-value threshhold to perform HMMer search against the Dfam database.}

\item{dfam.file}{path to pre-computed \code{\link{dfam.query}} output file. Can only be used in combination
with \code{annotate = "Dfam"}.}

\item{cluster}{shall predicted transposons be clustered with \code{\link{CLUSTpred}}.}

\item{clust.sim}{cluster reject if sequence similarity is lower than this threshold when performing clustering with \code{\link{CLUSTpred}}.}

\item{clust.file}{file path to pre-computed clustering file generated by \code{\link{CLUSTpred}}.}

\item{copy.number.est}{shall copy number estimation (including solo LTR prediction) be performed? Default is \code{copy.number.est = FALSE}.}

\item{fix.chr.name}{shall chromosome names be fixed?}

\item{cn.eval}{evalue for copy number estimation (BLAST hit threshold).}

\item{range}{define the genomic interval in which predicted LTR transposons shall be reported
. In case \code{range[1] = 1000} and \code{range[2] = 10000} then candidates are only 
reported if they start after position 1000 and end before position 10000 in their respective 
sequence coordinates. If \code{range[1] = 0} and \code{range[2] = 0}, 
so \code{range = c(0,0)} (default) then the entire genome is being scanned.}

\item{seed}{the minimum length for the exact maximal repeats. Only repeats with the specified minimum length are considered in all subsequent analyses. Default is \code{seed = 30}.}

\item{minlenltr}{minimum LTR length. Default is \code{minlenltr = 100}.}

\item{maxlenltr}{maximum LTR length. Default is \code{maxlenltr = 3500}.}

\item{mindistltr}{minimum distance of LTR starting positions. Default is \code{mindistltr = 4000}.}

\item{maxdistltr}{maximum distance of LTR starting positions. Default is \code{maxdistltr = 25000}.}

\item{similar}{minimum similarity value between the two LTRs in percent. \code{similar = 70}.}

\item{mintsd}{minimum target site duplications (TSDs) length. If no search for TSDs
shall be performed, then specify \code{mintsd = NULL}. Default is \code{mintsd = 4}.}

\item{maxtsd}{maximum target site duplications (TSDs) length. If no search for TSDs
shall be performed, then specify \code{maxtsd = NULL}. Default is \code{maxtsd = 20}.}

\item{vic}{number of nucleotide positions left and right (the vicinity) of the predicted
boundary of a LTR that will be searched for TSDs and/or one motif (if specified). 
Default is \code{vic = 60}.}

\item{overlaps}{specify how overlapping LTR retrotransposon predictions shall be treated. 
If \code{overlaps = "no"} is selected, then neither nested nor overlapping predictions will be reported in the output. In case \code{overlaps = "best"} is selected then in the case of two or more nested or overlapping predictions, solely the LTR retrotransposon prediction with
the highest similarity between its LTRs will be reported.
If \code{overlaps = "all"} is selected then all LTR retrotransposon predictions 
will be reported whether there are nested and/or overlapping predictions or not. 
Default is \code{overlaps = "best"}.}

\item{xdrop}{specify the xdrop value (> 0) for extending a seed repeat in both directions
allowing for matches, mismatches, insertions, and deletions. The xdrop extension process
 stops as soon as the extension involving matches, mismatches, insersions, and deletions 
 has a score smaller than T - X, where T denotes the largest score seen so far. Default is \code{xrop = 5}.}

\item{mat}{specify the positive match score for the X-drop extension process. Default is \code{mat = 2}.}

\item{mis}{specify the negative mismatch score for the X-drop extension process. Default is \code{mis = -2}.}

\item{ins}{specify the negative insertion score for the X-drop extension process. Default is \code{ins = -3}.}

\item{del}{specify the negative deletion score for the X-drop extension process. Default is \code{del = -3}.}

\item{motif}{specify 2 nucleotides for the starting motif and 2 nucleotides for the ending
motif at the beginning and the ending of each LTR, respectively.
Only palindromic motif sequences - where the motif sequence is equal to its complementary
sequence read backwards - are allowed, e.g. \code{motif = "tgca"}. Type the nucleotides without any space
separating them. If this option is not selected by the user, candidate pairs will not be
screened for potential motifs. If this options is set but no allowed number of
mismatches is specified by the argument \code{motifmis} and a search for the exact 
motif will be conducted. If \code{motif = NULL} then no explicit motif is being specified.}

\item{motifmis}{allowed number of mismatches in the TSD motif specified in \code{motif}.
The number of mismatches needs to be between [0,3].  Default is \code{motifmis = 0}.}

\item{aaout}{shall the protein sequence of the HMM matches to the predicted LTR transposon 
be generated as fasta file or not. Options are \code{aaout = "yes"} or \code{aaout = "no"}.}

\item{aliout}{shall the alignment of the protein sequence of the HMM matches to the predicted LTR transposon 
be generated as fasta file or not. Options are \code{aaout = "yes"} or \code{aaout = "no"}.}

\item{pptlen}{a two dimensional numeric vector specifying the minimum and maximum allowed
lengths for PPT predictions. If a purine-rich region that does not fulfill this range is
found, it will be discarded. Default is \code{pptlen = c(8,30)} (minimum = 8; maximum = 30).}

\item{uboxlen}{a two dimensional numeric vector specifying the minimum and maximum allowed
lengths for U-box predictions. If a T-rich region preceding a PPT that does not fulfill the PPT length criteria is
found, it will be discarded. Default is \code{uboxlen = c(3,30)} (minimum = 3; maximum = 30).}

\item{pptradius}{a numeric value specifying the area around the 3' LTR beginning to be 
considered when searching for PPT. Default value is \code{pptradius = 30}.}

\item{trnas}{path to the fasta file storing the unique tRNA sequences that shall be matched to the
predicted LTR transposon (tRNA library).}

\item{pbsalilen}{a two dimensional numeric vector specifying the minimum and maximum allowed
lengths for PBS/tRNA alignments. If the local alignments are shorter or longer than this
range, it will be discarded. Default is \code{pbsalilen = c(11,30)} (minimum = 11; maximum = 30).}

\item{pbsoffset}{a two dimensional numeric vector specifying the minimum and maximum allowed
distance between the start of the PBS and the 3' end of the 5' LTR. Local alignments not 
fulfilling this criteria will be discarded. Default is \code{pbsoffset = c(0,5)} (minimum = 0; maximum = 5).}

\item{pbstrnaoffset}{a two dimensional numeric vector specifying the minimum and maximum allowed
PBS/tRNA alignment offset from the 3' end of the tRNA. Local alignments not 
fulfilling this criteria will be discarded. Default is \code{pbstrnaoffset = c(0,5)} (minimum = 0; maximum = 5).}

\item{pbsmaxedist}{a numeric value specifying the maximal allowed unit edit distance in a
local PBS/tRNA alignment.}

\item{pbsradius}{a numeric value specifying the area around the 5' LTR end to be 
considered when searching for PBS Default value is \code{pbsradius = 30}.}

\item{hmms}{a character string or a character vector storing either the hmm files for
searching internal domains between the LTRs of predicted LTR transposons or a vector of
Pfam IDs from http://pfam.xfam.org/ that are downloaded and used to search for corresponding protein domains
within the predicted LTR transposons. As an option users can rename all of their hmm files
so that they start for example with the name \code{hmms = "hmm_*"}. This way all files starting with 
\code{hmm_} will be considered for the subsequent protein domain search. In case Pfam IDs 
are specified, the \code{LTRpred} function will automatically download the corresponding 
HMM files and use them for further protein domain searches. In case users prefer to specify 
Pfam IDs please specify them in the \code{pfam.ids} parmeter and choose \code{hmms = NULL}.}

\item{pdomevalcutoff}{a numeric value specifying the E-value cutoff for corresponding HMMER searches. All hits that do not fulfill this criteria are discarded. Default is \code{pdomevalcutoff = 1E-5}.}

\item{pbsmatchscore}{specify the match score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsmatchscore = 5}.}

\item{pbsmismatchscore}{specify the mismatch score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsmismatchscore = -10}.}

\item{pbsinsertionscore}{specify the insertion score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsinsertionscore = -20}.}

\item{pbsdeletionscore}{specify the deletion score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsdeletionscore = -20}.}

\item{pfam.ids}{a character vector storing the Pfam IDs from http://pfam.xfam.org/
that shall be downloaded and used to perform protein domain searches within the sequences
between the predicted LTRs.}

\item{cores}{number of cores to be used for multicore processing.}

\item{dfam.cores}{number of cores to be used for multicore processing when running Dfam query (in case \code{annotate = "Dfam"}).}

\item{hmm.cores}{number of cores to be used for multicore processing when performing hmmer protein search with \code{\link{LTRdigest}}.}

\item{orf.style}{type of predicting open reading frames (see documentation of USEARCH).}

\item{min.codons}{minimum number of codons in the predicted open reading frame.}

\item{trans.seqs}{logical value indicating wheter or not predicted open reading frames
shall be translated and the corresponding protein sequences stored in the output folder.}

\item{output.path}{a path/folder to store all results returned by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, and \code{LTRpred}. 
If \code{output.path = NULL} (Default) then a folder with the name of the input genome file
will be generated in the current working directory of R and all results are then stored in this folder.}

\item{quality.filter}{shall false positives be filtered out as much as possible? Default is \code{quality.filter = TRUE}. 
See \code{Description} for details.}

\item{n.orfs}{minimum number of Open Reading Frames that must be found between the LTRs (if \code{quality.filter = TRUE}). See \code{Details} for further information on quality control.}

\item{job_num}{a job number in case this function is run in parallel mode in \code{\link{LTRpred.meta}}.}

\item{verbose}{shall further information be printed on the console or not.}
}
\value{
The \code{LTRpred} function generates the following output files:

\itemize{
\item *_BetweenLTRSeqs.fsa : DNA sequences predicted LTR retrotransposons, in particular of the region between the LTRs in fasta format. 
\item *_Details.tsv : A spread sheet containing detailed information about the predicted LTRs.
\item *_FullLTRRetrotransposonSeqs.fsa : DNA sequences of the entire predicted LTR retrotransposon.
\item *_index.fsa.* : The suffixarray index file used to predict putative LTR retrotransposons.
\item *_Prediction.gff : A spread sheet containing detailed additional information about the predicted LTRs (partially redundant with the *_Details.tsv file).
\item *_index_ltrdigest.fsa : The suffixarray index file used to predict putative LTR retrotransposonswith \code{LTRdigest}.
\item *_LTRdigestPrediction.gff : A spread sheet containing detailed information about the predicted LTRs.
\item *-ltrdigest_tabout.csv : A spread sheet containing additional detailed information about the predicted LTRs.
\item *-ltrdigest_complete.fas : The full length DNA sequences of all predicted LTR transposons.
\item *-ltrdigest_conditions.csv : Contains information about the parameters used for a given
\code{LTRdigest} run.
\item *-ltrdigest_pbs.fas : Stores the predicted PBS sequences for the putative LTR retrotransposons.
\item *-ltrdigest_ppt.fas : Stores the predicted PPT sequences for the putative LTR retrotransposons.
\item *-ltrdigest_5ltr.fas and *-ltrdigest_3ltr.fas: Stores the predicted 5' and 3' LTR sequences. Note: If the direction of the putative retrotransposon could be predicted, these
files will contain the corresponding 3' and 5' LTR sequences. If no direction could be
predicted, forward direction with regard to the original sequence will be assumed by 
\code{LTRdigest}, i.e. the 'left' LTR will be considered the 5' LTR.
\item *-ltrdigest_pdom_<domainname>.fas : Stores the DNA sequences of the HMM matches
to the LTR retrotransposon candidates. 
\item *-ltrdigest_pdom_<domainname>_aa.fas : Stores the concatenated protein sequences of 
the HMM matches to the LTR retrotransposon candidates.
\item *-ltrdigest_pdom_<domainname>_ali.fas : Stores the alignment information for all matches of the given protein domain model to the translations of all candidates.
\item *_ORF_prediction_nt.fsa : Stores the predicted open reading frames within the predicted LTR transposons as DNA sequence.
\item *_ORF_prediction_aa.fsa : Stores the predicted open reading frames within the predicted LTR transposons as protein sequence.
\item *_LTRpred.gff : Stores the \code{LTRpred} predicted LTR transposons in GFF format.
\item *_LTRpred.bed : Stores the \code{LTRpred} predicted LTR transposons in BED format. 
\item *_LTRpred_DataSheet.tsv : Stores the output table as data sheet.
}
The ' * ' is an place holder for the name of the input genome file.

\strong{If annotate = "Dfam"}
In case the Dfam annotation option is specified the following
additional files are generated and stored in the LTRpred result folder:

\itemize{
\item *-ltrdigest_complete.fas_DfamAnnotation.out : a data table storing the output of the HMMer search of predicted retrotransposons against the Dfam database. 
}

\strong{if cluster = TRUE}
\itemize{
\item *-ltrdigest_complete.fas_CLUSTpred.blast6out : usearch cluster result in BLAST table format.
\item *-ltrdigest_complete.fas_CLUSTpred.log : log file of usearch run.
\item *-ltrdigest_complete.fas_CLUSTpred.uc : usearch cluster output file.
}
}
\description{
Main pipeline to perform
sufficient LTR retrotransposon predictions for any genome of interest.
}
\details{
This function provides the main pipeline to perform \code{de novo} LTR transposon
predictions.

\strong{Quality Control}

\itemize{
\item \code{ltr.similarity}: Minimum similarity between LTRs. All TEs not matching this
 criteria are discarded.
 \item \code{n.orfs}: minimum number of Open Reading Frames that must be found between the
  LTRs. All TEs not matching this criteria are discarded.
 \item \code{PBS or Protein Match}: elements must either have a predicted Primer Binding
 Site or a protein match of at least one protein (Gag, Pol, Rve, ...) between their LTRs. All TEs not matching this criteria are discarded.
 \item The relative number of N's (= nucleotide not known) in TE <= 0.1. The relative number of N's is computed as follows: absolute number of N's in TE / width of TE.
}
}
\examples{

\dontrun{
# generate de novo LTR transposon prediction
LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"),
          trnas       = system.file("hg38-tRNAs.fa", package = "LTRpred"),
          hmms        = paste0(system.file("HMMs/", package = "LTRpred"),"hmm_*"))

# run LTRpred with pre-computed predictions from LTRdigest()  
genome <- system.file("TAIR10_chr_all_LTRdigestPrediction.gff",package = "LTRpred")   
tabout <- system.file("TAIR10_chr_all-ltrdigest_tabout.csv",package = "LTRpred")
orf.pred <- system.file("nt.fa",package = "LTRpred")              
LTRpred(LTRdigest.gff = genome,
        tabout.file   = tabout,
        orf.file      = orf.pred)
}

}
\references{
R Edgar. Search and clustering orders of magnitude faster than BLAST. Bioinformatics (2010) 26 (19): 2460-2461.

D Ellinghaus, S Kurtz and U Willhoeft. LTRharvest, an efficient and flexible software for de novo detection of LTR retrotransposons. BMC Bioinformatics (2008). 9:18.

S Steinbiss et al. Fine-grained annotation and classification of de novo predicted LTR retrotransposons. Nucl. Acids Res. (2009) 37 (21): 7002-7013.
}
\seealso{
\code{\link{LTRharvest}}, \code{\link{LTRdigest}}, 
\code{\link{read.prediction}}, \code{\link{read.tabout}}, \code{\link{read.seqs}},
\code{\link{pred2fasta}}, \code{\link{pred2gff}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred2tsv.R
\name{pred2tsv}
\alias{pred2tsv}
\title{Format LTR prediction data to TSV file format}
\usage{
pred2tsv(LTR.data, output = "output.tsv")
}
\arguments{
\item{LTR.data}{the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.}

\item{output}{filename of the output CSV file.}
}
\value{
the \code{LTR.data} table in csv format saved to the hard drive.
}
\description{
This function formats the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}
to a \code{\link{data.frame}} in \code{CSV} file format.
}
\examples{
gff.file <- system.file("TAIR10_chr_all_LTRdigestPrediction.gff",
                        package = "LTRpred")
tabout.file <- system.file("TAIR10_chr_all-ltrdigest_tabout.csv"
                           ,package = "LTRpred")
LTRfile <- read.prediction(gff.file,tabout.file, program = "LTRdigest")

# generate csv file
pred2csv(LTRfile$ltr.retrotransposon)
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred2gff.R
\name{pred2gff}
\alias{pred2gff}
\title{Format LTR prediction data to GFF3 file format}
\usage{
pred2gff(LTR.data, output = "output.gff", program = "LTRpred")
}
\arguments{
\item{LTR.data}{the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.}

\item{output}{filename of the output GFF file.}

\item{program}{program used to generate the prediction table, e.g. \code{program = "LTRpred"},
\code{program = "LTRdigest"}, or \code{program = "LTRharvest"}.}
}
\description{
This function formats the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}
to a \code{\link{data.frame}} in \code{GFF3} file format.
}
\details{
The GFF3 file format is defined by: chromosome; start; end; name; columns.
}
\examples{
gff.file <- system.file("TAIR10_chr_all_LTRdigestPrediction.gff",
                        package = "LTRpred")
tabout.file <- system.file("TAIR10_chr_all-ltrdigest_tabout.csv"
                           ,package = "LTRpred")
LTRfile <- read.prediction(gff.file,tabout.file, program = "LTRdigest")

# generate GFF file
pred2gff(LTRfile$ltr.retrotransposon, output = "test.gff")

}
\references{
http://www.ensembl.org/info/website/upload/gff.html
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.seq.space.R
\name{meta.seq.space}
\alias{meta.seq.space}
\title{Quantification of the repeat space for multiple \code{Repeat Masker} files}
\usage{
meta.seq.space(
  rm.folder,
  type = "rm",
  repeat.type = "all",
  rename = FALSE,
  save.file = TRUE
)
}
\arguments{
\item{rm.folder}{file path to the \code{Repeat Masker} files.}

\item{type}{type of input data. Option are:
\itemize{
 \item \code{type = "rm"} - \code{Repeat Masker} files.
}}

\item{repeat.type}{type of element for which total repeat space shall be quantified.
Options are:

\itemize{
 \item \code{repeat.type = "all"} include all types of repeats.
 \item \code{repeat.type = "ltr"} include only LTR retrotransposons.
}}

\item{rename}{logical value indicating whether or not long file names shall be transformed to short scientific names.}

\item{save.file}{logical value indicating whether or not the result \code{tibble} shall be stored locally.}
}
\description{
The repeat space for multiple \code{Repeat Masker} files is quantified and
stored in one \code{tibble}.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LTRdigest.R
\name{LTRdigest}
\alias{LTRdigest}
\title{Run LTRdigest to predict putative LTR Retrotransposons}
\usage{
LTRdigest(
  input.gff3,
  genome.file,
  aaout = "yes",
  aliout = "yes",
  pptlen = c(8, 30),
  uboxlen = c(3, 30),
  pptradius = 30,
  trnas = NULL,
  pbsalilen = c(11, 30),
  pbsoffset = c(0, 5),
  pbstrnaoffset = c(0, 5),
  pbsmaxedist = 1,
  pbsradius = 30,
  hmms = NULL,
  pdomevalcutoff = 1e-05,
  pbsmatchscore = 5,
  pbsmismatchscore = -10,
  pbsinsertionscore = -20,
  pbsdeletionscore = -20,
  pfam.ids = NULL,
  cores = 1,
  index.file = NULL,
  output.path = NULL
)
}
\arguments{
\item{input.gff3}{path to the prediction file in gff3 format returned by \code{\link{LTRharvest}}.}

\item{genome.file}{path to the genome file in \code{fasta} format.}

\item{aaout}{shall the protein sequence of the HMM matches to the predicted LTR transposon 
be generated as fasta file or not. Options are \code{aaout = "yes"} or \code{aaout = "no"}.}

\item{aliout}{shall the alignment of the protein sequence of the HMM matches to the predicted LTR transposon 
be generated as fasta file or not. Options are \code{aaout = "yes"} or \code{aaout = "no"}.}

\item{pptlen}{a two dimensional numeric vector specifying the minimum and maximum allowed
lengths for PPT predictions. If a purine-rich region that does not fulfill this range is
found, it will be discarded. Default is \code{pptlen = c(8,30)} (minimum = 8; maximum = 30).}

\item{uboxlen}{a two dimensional numeric vector specifying the minimum and maximum allowed
lengths for U-box predictions. If a T-rich region preceding a PPT that does not fulfill the PPT length criteria is
found, it will be discarded. Default is \code{uboxlen = c(3,30)} (minimum = 3; maximum = 30).}

\item{pptradius}{a numeric value specifying the area around the 3' LTR beginning to be 
considered when searching for PPT. Default value is \code{pptradius = 30}.}

\item{trnas}{path to the fasta file storing the unique tRNA sequences that shall be matched to the
predicted LTR transposon (tRNA library).}

\item{pbsalilen}{a two dimensional numeric vector specifying the minimum and maximum allowed
lengths for PBS/tRNA alignments. If the local alignments are shorter or longer than this
range, it will be discarded. Default is \code{pbsalilen = c(11,30)} (minimum = 11; maximum = 30).}

\item{pbsoffset}{a two dimensional numeric vector specifying the minimum and maximum allowed
distance between the start of the PBS and the 3' end of the 5' LTR. Local alignments not 
fulfilling this criteria will be discarded. Default is \code{pbsoffset = c(0,5)} (minimum = 0; maximum = 5).}

\item{pbstrnaoffset}{a two dimensional numeric vector specifying the minimum and maximum allowed
PBS/tRNA alignment offset from the 3' end of the tRNA. Local alignments not 
fulfilling this criteria will be discarded. Default is \code{pbstrnaoffset = c(0,5)} (minimum = 0; maximum = 5).}

\item{pbsmaxedist}{a numeric value specifying the maximal allowed unit edit distance in a
local PBS/tRNA alignment.}

\item{pbsradius}{a numeric value specifying the area around the 5' LTR end to be 
considered when searching for PBS Default value is \code{pbsradius = 30}.}

\item{hmms}{a character string or a character vector storing either the hmm files for
searching internal domains between the LTRs of predicted LTR transposons or a vector of
Pfam IDs from http://pfam.xfam.org/ that are downloaded and used to search for corresponding protein domains
within the predicted LTR transposons. As an option users can rename all of their hmm files
so that they start for example with the name \code{hmms = "hmm_*"}. This way all files starting with 
\code{hmm_} will be considered for the subsequent protein domain search. In case Pfam IDs 
are specified, the \code{LTRpred} function will automatically download the corresponding 
HMM files and use them for further protein domain searches. In case users prefer to specify 
Pfam IDs please specify them in the \code{pfam.ids} parmeter and choose \code{hmms = NULL}.}

\item{pdomevalcutoff}{a numeric value specifying the E-value cutoff for corresponding HMMER searches. All hits that do not fulfill this criteria are discarded. Default is \code{pdomevalcutoff = 1E-5}.}

\item{pbsmatchscore}{specify the match score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsmatchscore = 5}.}

\item{pbsmismatchscore}{specify the mismatch score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsmismatchscore = -10}.}

\item{pbsinsertionscore}{specify the insertion score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsinsertionscore = -20}.}

\item{pbsdeletionscore}{specify the deletion score used in the PBS/tRNA Smith-Waterman alignment.
Default is \code{pbsdeletionscore = -20}.}

\item{pfam.ids}{a character vector storing the Pfam IDs from http://pfam.xfam.org/
that shall be downloaded and used to perform protein domain searches within the sequences
between the predicted LTRs.}

\item{cores}{number of cores to be used for multicore processing.}

\item{index.file}{specify the name of the enhanced suffix array index file that is computed
by \code{suffixerator}. This opten can be used in case the suffix file was previously 
generated, e.g. during a previous call of this function. In this case the suffix array index
file does not need to be re-computed for new analyses. This is particularly useful when 
running \code{LTRdigest} with different parameter settings.}

\item{output.path}{a path/folder to store all results returned by \code{LTRdigest}. 
If \code{output.path = NULL} (Default) then a folder with the name of the input genome file
will be generated in the current working directory of R and all results are then stored in this folder.}
}
\value{
The \code{LTRdigest} function generates the following output files:

\itemize{
\item *_index_ltrdigest.fsa : The suffixarray index file used to predict putative LTR retrotransposonswith \code{LTRdigest}.
\item *_LTRdigestPrediction.gff : A spread sheet containing detailed information about the predicted LTRs.
\item *-ltrdigest_tabout.csv : A spread sheet containing additional detailed information about the predicted LTRs.
\item *-ltrdigest_complete.fas : The full length DNA sequences of all predicted LTR transposons.
\item *-ltrdigest_conditions.csv : Contains information about the parameters used for a given
\code{LTRdigest} run.
\item *-ltrdigest_pbs.fas : Stores the predicted PBS sequences for the putative LTR retrotransposons.
\item *-ltrdigest_ppt.fas : Stores the predicted PPT sequences for the putative LTR retrotransposons.
\item *-ltrdigest_5ltr.fas and *-ltrdigest_3ltr.fas: Stores the predicted 5' and 3' LTR sequences. Note: If the direction of the putative retrotransposon could be predicted, these
files will contain the corresponding 3' and 5' LTR sequences. If no direction could be
predicted, forward direction with regard to the original sequence will be assumed by 
\code{LTRdigest}, i.e. the 'left' LTR will be considered the 5' LTR.
\item *-ltrdigest_pdom_<domainname>.fas : Stores the DNA sequences of the HMM matches
to the LTR retrotransposon candidates. 
\item *-ltrdigest_pdom_<domainname>_aa.fas : Stores the concatenated protein sequences of 
the HMM matches to the LTR retrotransposon candidates.
\item *-ltrdigest_pdom_<domainname>_ali.fas : Stores the alignment information for all matches of the given protein domain model to the translations of all candidates.
}
The ' * ' is an place holder for the name of the input genome file.
}
\description{
This function implements an interface between R and
the LTRdigest command line tool to predict putative LTR retrotransposons from R.
}
\details{
The \code{LTRdigest} function is a wrapper function to work with the
call the \code{LTRdigest} command line tool from R.
}
\examples{
\dontrun{
# Run LTRharvest for Arabidopsis thaliana using standard parameters
LTRharvest(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))

# Run LTRdigest for Arabidopsis thaliana using standard parameters
LTRdigest(input.gff3  = "Hsapiens_ChrY_ltrharvest/Hsapiens_ChrY_Prediction.gff", 
          genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"),
          trnas       = system.file("hg38-tRNAs.fa", package = "LTRpred"),
          hmms        = paste0(system.file("HMMs/", package = "LTRpred"),"hmm_*"))
}
}
\references{
S Steinbiss et al. Fine-grained annotation and classification of de novo predicted LTR retrotransposons. Nucl. Acids Res. (2009) 37 (21): 7002-7013.
}
\seealso{
\code{\link{LTRharvest}},  \code{\link{LTRpred}},
\code{\link{read.prediction}}, \code{\link{read.seqs}},
\code{\link{pred2fasta}}, \code{\link{pred2fasta}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred2fasta.R
\name{pred2fasta}
\alias{pred2fasta}
\title{Save the sequence of the predicted LTR Transposons in a fasta file}
\usage{
pred2fasta(LTRpred.tbl, prediction.file, output = "output.fa")
}
\arguments{
\item{LTRpred.tbl}{the \code{\link{data.frame}} generated by \code{\link{LTRpred}}.}

\item{prediction.file}{the fasta file storing either the full LTR Transposon sequence or only the LTR sequence
as returned by the \code{\link{LTRpred}} function.}

\item{output}{the fasta file to which the output sequences shall be stored in.}
}
\description{
This function allows users to save the sequence of the predicted LTR Transposons or LTRs in a fasta file.
}
\details{
The output \code{data.frame}s returned by \code{\link{LTRpred}} contain all information of the predicted
LTR retrotransposons that can be used for post-filtering steps. After these post-filtering steps
sequences of the remaining (filtered) candidates can be retrieved by this function.
}
\examples{
\dontrun{
# Hypothetical Example
# Generate LTR transposon prediction for A. thaliana 
Ath.Pred <- LTRpred(input.gff3        = "TAIR10_chr_all/TAIR10_chr_all_Prediction.gff", 
                    genome.file       = "Genome/TAIR10_chr_all.fas",
                    trnas             = "araTha1-tRNAs.fa",
                    hmms              = "hmm_*",
                    cores             = 1)

# Filter for LTR transposons having 100 pec. sequence similarity between their LTRs
FilteredLTRTransposons <- dplyr::filter(Ath.Pred, ltr_similarity == 100)

# Write the sequences of these filtered LTR transposons to a fasta file
pred2fasta(LTRpred.tbl = FilteredLTRTransposons, 
           prediction.file = "TAIR10_chr_all-ltrdigest_complete.fas", 
           output                = "AthalianaPutativeLTRTransposons.fa")
}
}
\seealso{
\code{\link{LTRharvest}}, \code{\link{LTRdigest}}, \code{\link{LTRpred}}, \code{\link{read.prediction}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy.datasheet.R
\name{tidy.datasheet}
\alias{tidy.datasheet}
\title{Select most important columns of \code{LTRpred} output for further analytics}
\usage{
tidy.datasheet(LTRpred.tbl)
}
\arguments{
\item{LTRpred.tbl}{a \code{data.frame} storing the result (DataSheet) of \code{\link{LTRpred}}.}
}
\description{
This function simply selects the most important columns of \code{LTRpred} output for further analytics.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred2bed.R
\name{pred2bed}
\alias{pred2bed}
\title{Format LTR prediction data to BED file format}
\usage{
pred2bed(LTR.data, output = "output")
}
\arguments{
\item{LTR.data}{the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.}

\item{output}{filename of the output BED file.}
}
\description{
This function formats the LTR prediction \code{\link{data.frame}}
generated by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}
to a \code{\link{data.frame}} in \code{BED} file format.
}
\details{
The BED file format is defined by: chromosome; start; end; name; columns.
}
\examples{
gff.file <- system.file("TAIR10_chr_all_LTRdigestPrediction.gff",
                        package = "LTRpred")
tabout.file <- system.file("TAIR10_chr_all-ltrdigest_tabout.csv"
                           ,package = "LTRpred")
LTRfile <- read.prediction(gff.file,tabout.file, program = "LTRdigest")

# generate BED file
pred2bed(LTRfile$ltr.retrotransposon)

}
\references{
http://www.ensembl.org/info/website/upload/bed.html
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ws.wrap.path.R
\name{ws.wrap.path}
\alias{ws.wrap.path}
\title{Wrap whitespace in paths}
\usage{
ws.wrap.path(paths)
}
\arguments{
\item{paths}{a string or vector of strings containing paths.}
}
\description{
Helper function for wrapping white spaces in folder paths.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcolor.R
\name{bcolor}
\alias{bcolor}
\title{Colors scheme for plots}
\usage{
bcolor(n)
}
\arguments{
\item{n}{number of different colors that shall be returned.}
}
\description{
A distinguished color scheme for scientific plots.
}
\examples{
# return 5 colors 
bcolor(5)
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repbase.query.R
\name{repbase.query}
\alias{repbase.query}
\title{Query the RepBase to annotate putative LTRs}
\usage{
repbase.query(
  seq.file,
  repbase.path,
  output = "RepbaseOutput.txt",
  max.hits = 5000,
  eval = 1e-30,
  cores = 1
)
}
\arguments{
\item{seq.file}{file path to the putative LTR transposon sequences in \code{fasta} format.}

\item{repbase.path}{file path to the RepBase file in \code{fasta} format.}

\item{output}{file name of the BLAST output.}

\item{max.hits}{maximum number of hits that shall be retrieved that still fulfill the e-value criterium.
Default is \code{max.hits = 5000}.}

\item{eval}{e-value threshold for BLAST hit detection. Default is \code{eval = 1E-30}.}

\item{cores}{number of cores to use to perform parallel computations.}
}
\description{
Validate or annotate putative LTR
transposons that have been predicted using LTRharvest or LTRdigest.
}
\details{
The RepBase database provides a collection of curated transposable element annotations.

This function allows users to validate or annotate putative LTR
transposons that have been predicted using LTRharvest or LTRdigest by blasting predicted LTR transposons 
to transposons known (annotated) in other species (e.g. such as Arabidopsis thaliana).

Internally, this function performs a \code{blastn} search of the putative LTR transposons predicted
by LTRharvest or LTRdigest against the Repbase fasta file that is specified by the user.

For this purpose it is required that the user has a working version of BLAST+ running on his or her machine.
}
\examples{
\dontrun{
# Example annotation run against the A thaliana RepBase using 4 cores
q <- repbase.query(seq.file     = "path/to/LTRtransposonSeqs.fasta",
                  repbase.path = "path/to/Athaliana_repbase.ref",
                  cores        = 4)
                 
Annot <- dplyr::select(dplyr::filter(dplyr::group_by(q,query_id), 
                                    (bit_score == max(bit_score))),
                                     query_id:q_len,evalue,bit_score,scope)
# select only hits with a scope > 0.1
Annot.HighMatches <- dplyr::filter(Annot, scope >= 0.1)
# Annotate the proportion of hits
barplot(sort(table(unlist(lapply(stringr::str_split(
        names(table(Annot.HighMatches$subject_id)),"_"), 
        function(x) x[2]))), decreasing = TRUE))
}
}
\references{
http://www.girinst.org/repbase/

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.

Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.

Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.

Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.

Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.orfs.R
\name{read.orfs}
\alias{read.orfs}
\title{Read output of \code{ORFpred}}
\usage{
read.orfs(input.file)
}
\arguments{
\item{input.file}{fasta file generated by \code{\link{ORFpred}}.}
}
\value{
A \code{\link{data.frame}} object storing the \code{seq.id}, \code{orfs} (number of predicted ORFs), \code{start}, and \code{end} of the predicted LTRs.
}
\description{
This function reads the output of the \code{\link{ORFpred}} function and stores the 
sequence id and number of predicted ORFs in a \code{\link{data.frame}} object.
}
\details{
The file generated by \code{\link{ORFpred}} is parsed by this function
and returned as \code{\link{data.frame}} object.
}
\examples{
# read an example prediction file generated by PredictORFs()
ORFPred <- read.orfs(system.file("nt.fa",package = "LTRpred"))

head(ORFPred)
}
\seealso{
\code{\link{ORFpred}}, \code{\link{LTRpred}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ltr.cn.R
\name{ltr.cn}
\alias{ltr.cn}
\title{Copy Number Quantification of predicted LTRs (solo LTR prediction)}
\usage{
ltr.cn(
  data.sheet,
  LTR.fasta_3ltr,
  LTR.fasta_5ltr,
  genome,
  ltr.similarity = 70,
  scope.cutoff = 0.85,
  perc.ident.cutoff = 70,
  output = NULL,
  max.hits = 500,
  eval = 1e-10,
  cores = 1
)
}
\arguments{
\item{data.sheet}{path to the \code{*_LTRpred_DataSheet.csv} file generated by \code{\link{LTRpred}}.}

\item{LTR.fasta_3ltr}{path to fasta file storing sequences of 3 prime LTRs generated by \code{\link{LTRpred}}.}

\item{LTR.fasta_5ltr}{path to fasta file storing sequences of 5 prime LTRs generated by \code{\link{LTRpred}}.}

\item{genome}{file path to the reference genome in which solo LTRs shall be found (in \code{fasta} format).}

\item{ltr.similarity}{similarity threshold for defining LTR similarity.}

\item{scope.cutoff}{similarity threshold for the scope variable. The scope of a BLAST hit is defined by \code{abs(s_len - alig_length) / s_len}
and aims to quantify the scope ('length similarity between hit and query') of the alignment. Default is \code{scope.cutoff = 0.85} meaning that
at least 85\% of the length of the query sequence must match with the length of the subject sequence.}

\item{perc.ident.cutoff}{choose the minimum sequence identity between the query sequence and subject hit sequence that have a scope >= \code{scope.cutoff} (e.g. 0.85). 
Default is \code{perc.ident.cutoff = 70}. This threshold aims to detect BLAST hits that have similar sequence lengths 
(controlled by the \code{scope.cutoff} variable; e.g. \code{scope.cutoff = 0.85}) and within this scope at least e.g. 70\% sequence identity (controlled by the \code{perc.ident.cutoff} variable).}

\item{output}{file name of the BLAST output. If \code{output = NULL} (default) then the BLAST output file will be deleted after the result \code{data.frame} is returned by this function.}

\item{max.hits}{maximum number of hits that shall be retrieved that still fulfill the e-value criterium.
Default is \code{max.hits = 65000}.}

\item{eval}{e-value threshold for BLAST hit detection. Default is \code{eval = 1E-10}.}

\item{cores}{number of cores for parallel computations.}
}
\description{
Detect solo LTR copies and genomic locations of predicted LTR transposons using a BLAST search strategy.
}
\references{
Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.

Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.

Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.

Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.

Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
}
\author{
Hajk-Georg Drost
}
