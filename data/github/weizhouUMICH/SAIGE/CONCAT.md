
Table of Contents
=================

   * [Introduction](#introduction)
   * [Citation](#citation)
   * [How to install SAIGE and SAIGE-GENE](#how-to-install-and-run-saige-and-saige-gene)
   * [Notes for users before running jobs](#notes-for-users-before-running-jobs)
   * [UK Biobank GWAS Results](#uk-biobank-gwas-results)
   * [Log for fixing bugs](#log-for-fixing-bugs)
   

# Introduction
## Manuscript for SAIGE-GENE+: https://www.medrxiv.org/content/10.1101/2021.07.12.21260400v1

## Current version is 0.45 (Updated on January 24, 2022) - comment out the part to estimate the effective sample sizes, which may not convert and take very long; put <= instead of < for maxMAF in the gene-based tests

## Current version is 0.44.6.5 (Updated on August 18, 2021) - 0.44.6.2 add extdata/extractNglmm.R to extract the effective sample size without running Step 1. extdata/cmd_extractNeff.sh has the pipeline. The effective sample size (Nglmm) is differently calculated than the previous versions.

## Previous version is 0.44.6.1 (Updated on July 16, 2021) - SAIGE-GENE+: for group tests, collpasing ultra-rare variants with MAC <= 10. Set --method_to_CollapseUltraRare="absence_or_presence" as default to collpase ultra-rare varaints with MAC <= 10. SAIGE-GENE+ has well controlled type I error rates when the maximum MAF cutoff (maxMAFforGroupTest) is lower than 1%, e.g. 0.01% or 0.1%. Tests with multiple MAF cutoffs and variant annotations can be combined using the Cauchy combination (function CCT)

##Please re-install 0.44.2 if you installed this verion on March 31. 

##For BGEN input, 8 bits are required. 

##For BGEN input in step 2 with missing dosages, Please use version 0.38 or later.


SAIGE is an R package with Scalable and Accurate Implementation of Generalized mixed model (Chen, H. et al. 2016). It accounts for sample relatedness and is feasible for genetic association tests in large cohorts and biobanks (N > 400,000).

SAIGE performs single-variant association tests for binary traits and quantitative taits. For binary traits, SAIGE uses the saddlepoint approximation (SPA)(mhof, J. P. , 1961; Kuonen, D. 1999; Dey, R. et.al 2017) to account for case-control imbalance.

SAIGE-GENE (implemented in the SAIGE R package) performs gene- or region-based association tests (Burde, SKAT, SKAT-O) for binary traits and quantitative traits. Note: SAIGE-GENE accounts for case-control imbalance in gene-based tests (>= 0.35.8.5)


# Citation
The SAIGE manuscript:
Wei Zhou, Jonas B. Nielsen, Lars G. Fritsche, Maiken B. Elvestad, Brooke Wolford, Maoxuan Lin, Kristian Hveem, Hyun Min Kang, Goncalo R. Abecasis, Cristen J. Willer*, Seunggeun Lee* “Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies.” Nature Genetics 50, 1335–1341 (2018)

The SAIGE-GENE pre-print:
https://www.biorxiv.org/content/10.1101/583278v2


# How to install and run SAIGE and SAIGE-GENE


## Install SAIGE/SAIGE-GENE

### List of dependencies: 

* R-3.6.1, gcc >= 5.4.0, cmake 3.14.1, [cget](https://cget.readthedocs.io/en/latest/src/intro.html#installing-cget)
* R packages: "R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "SKAT","MetaSKAT"
* /extdata/install_packages.R can be used to install the R packages
* SAIGE v0.39.2 depends on the SPAtest v3.1.2
* MetaSKAT is currently not available on CRAN. Please install it from github using R
  ``` 
   devtools::install_github("leeshawn/MetaSKAT") 
  ```

###  Install SAIGE from conda

#### Warning: please do not use this bioconda version for bgen input. We are working on the issue. 

![r-saige](https://anaconda.org/bioconda/r-saige/badges/version.svg)
![latest_update](https://anaconda.org/bioconda/r-saige/badges/latest_release_date.svg)

To install saige from conda simply create environment with latest version of R and saige:
```
conda create -n saige -c conda-forge -c bioconda "r-base>=4.0" r-saige
conda activate saige
```

More info on [r-saige conda package](https://anaconda.org/bioconda/r-saige) and available versions can be found in the [issue #272](https://github.com/weizhouUMICH/SAIGE/issues/272).

###  Install SAIGE using the conda environment

1. Create a conda environment using 
     ([conda environment file](https://github.com/weizhouUMICH/SAIGE/blob/master/conda_env/environment-RSAIGE.yml)) 
     Here is a link to download the [conda environment file](https://raw.githubusercontent.com/weizhouUMICH/SAIGE/master/conda_env/environment-RSAIGE.yml)

     After downloading environment-RSAIGE.yml, run following command
     ```
       conda env create -f environment-RSAIGE.yml
   ```

2. Activate the conda environment RSAIGE

     ```
       conda activate RSAIGE
       FLAGPATH=`which python | sed 's|/bin/python$||'`
       export LDFLAGS="-L${FLAGPATH}/lib"
       export CPPFLAGS="-I${FLAGPATH}/include"
     ```
Please make sure to set up the LDFLAGS and CPPFLAGS using export (the last two command lines), so libraries can be linked correctly when the SAIGE source code is compiled. Note: [Here](https://github.com/weizhouUMICH/SAIGE/blob/master/conda_env/createCondaEnvSAIGE_steps.txt) are the steps to create the conda environment file 

3. Open R, run following script to install the MetaSKAT R library.
   
     ```
       devtools::install_github("leeshawn/MetaSKAT") 
     ```

4. Install SAIGE from the source code. 

     Method 1: 

     ```
       src_branch=master
       repo_src_url=https://github.com/weizhouUMICH/SAIGE
       git clone --depth 1 -b $src_branch $repo_src_url

       R CMD INSTALL --library=path_to_final_SAIGE_library SAIGE
     ```
     
     When call SAIGE in R, set lib.loc=path_to_final_SAIGE_library   

     ```
       library(SAIGE, lib.loc=path_to_final_SAIGE_library)
     ```

    Method 2: 

    Open R. Run

    ```
      devtools::install_github("weizhouUMICH/SAIGE")
    ```

### Run SAIGE using a docker image 

Thanks to Juha Karjalainen for sharing the Dockerfile. 
The docker image can be pulled

```
docker pull wzhou88/saige:0.45
```

Functions can be called
```
step1_fitNULLGLMM.R --help
step2_SPAtests.R --help
createSparseGRM.R --help
```


## Run SAIGE for single-variant association tests and SAIGE-GENE for gene- or region-based tests

Here is a wiki page containg tutorial to run SAIGE and SAIGE-GENE
  https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE
  
### Examples

Example data and script can be found in ./extdata. Run

    bash cmd.sh

to run single-variant and gene-based association tests


## extract effectize sample size v0.44.6.2)
```
  SAIGE_extractNeff.R --help
  bash cmd_extractNeff.sh
```    


# Notes before running jobs

### FAQ can be found  [here](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE#Frequently-asked-questions)

### More notes
1. Since the SPA test always provides close to 0 p-values for variants with MAC < 3, please use at least minMAC = 3 to filter out the results
2. When query is used for bgen files, please make sure there are no duplicate SNP ids in the list
3. If the error message "Error in setgeno(genofile, subSampleInGeno, memoryChunk) :
  vector::_M_range_check", try use a smaller memeoryChunk, such as 2
4. IMPORTANT:In version <= 0.26, for binary traits, BETA is for alt allele and for quantitative traits, BETA is for minor allele 
5. Please note that LOCO only works for autosomal genetic variants. For non-autosomal genetic variants, please leave LOCO=FALSE in step 2.
6. SAIGE-GENE 0.36.3 and 0.36.3.1 now output an effect size for burden tests with the option IsOutputBETASEinBurdenTest in step2. Please note that the magnitude of the effect size is difficult to interpret. 
7. We haven't throughly tested the program on a small sample size. All simulation studies were done using 10,000 samples. Similar to BOLT-LMM, SAIGE uses asymptotic approaches to for feasibility on large samples. Based on our previous real-data analysis, we saw the performance on 3,000 samples were fine. 

# UK Biobank GWAS Results
1. The GWAS results for binary phenotypes in UK Biobank (1,283 phenotypes) using SAIGE are currently available for public download at

https://www.leelabsg.org/resources

Pheweb browser for the UK Biobank results

http://pheweb.sph.umich.edu/SAIGE-UKB/


*This research has been conducted using the UK Biobank Resource under application number 24460.

2. The exome-wide gene-based association results for quantitative traits in UK Biobank (53 traits) using SAIGE-GENE are currently available for public download at

https://www.leelabsg.org/resources

*This research has been conducted using the UK Biobank Resource under application number 45227.


# Log for fixing bugs

* 0.45 (January-24-2022). comment out the part to estimate the effective sample sizes, which may not convert and take very long; put <= instead of < for maxMAF in the gene-based tests

* 0.44.6.5 (August-19-2021). fix the SE=0 issue when IsOutputlogPforSingle=TRUE

* 0.44.6.4 (August-16-2021). make IsOutputlogPforSingle work for quantitative traits. remove the rsid in the output when the input is bgen

* 0.44.6.2 (August-2-2021). add extdata/extractNglmm.R to extract the effective sample size without running Step 1. extdata/cmd_extractNeff.sh has the pipeline. The effective sample size (Nglmm) is differently calculated than the previous versions. 

* 0.44.6.1 (July-16-2021). add the function CCT to perform Cauchy combination to combine multipel tests

* 0.44.6 (July-13-2021). Set --method_to_CollapseUltraRare="absence_or_presence" as default to collpase ultra-rare varaints with MAC <= 10. We call this version SAIGE-GENE+. SAIGE-GENE+ has well controlled type I error rates when the maximum MAF cutoff (maxMAFforGroupTest) is lower than 1%, e.g. 0.01% or 0.1%.  

* 0.44.5 (April-21-2021). 1. re-write code for leave-one-chromosome-out in Step 1 to have more efficient parallel computation in Step 1. 2. Speed up the single-variant association tests when running gene-based tests

* 0.44.2 (March-31-2021) 1.add an option useSparseGRMtoFitNULL to allow for fitting the null model using the sparse GRM and 2. add options to collapse the ultra-rare variants in the set-based tests. --method_to_CollapseUltraRare, --MACCutoff_to_CollapseUltraRare, --DosageCutoff_for_UltraRarePresence

* 0.44.1 (Feb-16-2021) 1. Fixed the error " X %*% Z : non-conformable arguments" for monomorphic variants. 2. merged Jonathon's codes to update savvy to savvy 2.0. For markers in VCF or SAV files without imputation info R2 values, the imputationInfo column will be 1 in the output file, so the markers will not but removed by minInfo 

* 0.44 (January-11-2021) 1. Fixed the error "Phi_ccadj[-indexNeg, -indexNeg]"; 2.  inverse normalization is only performed for quantitative traits; 3. For step 2, bgen input requires the sample file. vcf input does not require a seperate sample file. If sample file is not provided, sample ids will be read from vcf file

* 0.43.3 (January-05-2021)  error "FALis_rewrite_XnonPAR_forMalesSE not found" has been fixed

* 0.43.2 (December-13-2020)  add scripts to calcuate the effectize sample size in Step 1 for binary traits

* 0.43.1. with LOCO=TRUE, remove model results for other chromosomes to save memory usage for Step 2. 

* 0.43 (November-21-2020) Further modify the sparse version of the score test for quantitative traits. This causes slight different assoc tests for variants with MAF < 0.05 for quantitative traits. Set LOCO = TRUE to the default values for step 1 and step 2. In step 2, --chrom needs to be specified for LOCO=TRUE.

* 0.42.1 (September-21-2020) uncomment isSparse=FALSE for quantitative traits. This was commented out for testing in 0.42

* 0.42 (September-16-2020) fix a bug for variance ratio adjustion when account for case-control imbalance for gene-based tests. minMAC is set to 1/(2*N) instead of 0 if is_rewrite_XnonPAR_forMales=TRUE

* 0.41 (August-30-2020) improve the LOCO feature, implement LOCO for gene- and region- based tests (require --chrom to be specified), and with minInfo cutoff, if the input VCF files do not contain info scores, info will be output as NA and markers won't be filtered out. fixed an issue when subsetting pre-calcuated terms (regress X out of G) to drop missing dosages. Use sparse matrices for genotypes/dosages in gene- and region- based tests, so memory usage is dramatically decreased

* 0.39.4 (August-11-2020) use sparse matrix to represent genotype matrix for gene-based tests to save memory

* 0.39.3 (August-6-2020)  add five options --sexCol, --FemaleCode, --FemaleOnly, --MaleCode, --MaleOnly to perform sex-specific Step 1.

* 0.39.2 (July-27-2020)
** add three options --sampleFile_male, --X_PARregion, --is_rewrite_XnonPAR_forMales for chromosome X association tests, in which genotypes/dosages of non-PAR region of males will be multiplied by 2 

* 0.39.1 (July-27-2020)
** add an option --IsOutputlogPforSingle to output log(P) for single-variant assoc tests. v0.39.1 requires SPAtest 3.1.2.  

* 0.39 (May-27-2020)
** fixed an error when conditional analysis is conducted based on vcf input (introduced in 0.38)

* 0.38 (May-4-2020)
** further fixed the bug for output the allele 2 when bgen input with missing dosages was used and missing dosages were dropped. 
** sampleFile is no longer needed if VCF file is used in Step 2
** add --IsOverwriteVarianceRatioFile in step 1 to overwrite the variance ratio file

* 0.37 (May-1-2020)
** fixed an issue with AC values when bgen input is used with missing dosages to be mean imputed (default setting).
 
* 0.36.6 (April-15-2020)
** add an option IsOutputHetHomCountsinCaseCtrl to output the heterozygous and homozygous counts in cases and controls

* 0.36.5.1 (March-29-2020)
** add the option SPAcutoff, If the test statistic lies within the standard deviation cutoff of the mean, p-value based on traditional score test is returned. Otherwise, SPA will be applied. Default value of SPAcutoff is 2 (corresponding p.value.NA 0.05

* 0.36.5 (March-29-2020)
** Fix a typo to extract p.value. 0.36.5: fix an issue for LOCO=TRUE. This issue was introduced when the option minMAFforGRM was introduced.

* 0.36.4.2 (March-20-2020)
** Fix a bug by unlist(p.value), which was introduced in 0.36.4

* 0.36.4.1 (March-18-2020)
** Trying to fix a bug when minMAFforGRM is set and LOCO=TRUE


* 0.36.4 (March-18-2020)
** add an option includeNonautoMarkersforVarRatio in step 1. If TRUE, non-autosomal markers are also used for variance ratio estimation, which will make the algorithm more appropriate for assoc tests for non-autosomal markers; use the new function with sparse sigma for p-values for single variants in gene-based tests; assign AF to be 0 if all samples have missing genotypes or dosages


* 0.36.3.2 (February-25-2020)
** Bug fixed: 1. fixed a bug for gene-based conditioning tests with multiple conditioning markers 2. add codes to re-check markers after dropping samples with missing dosages/genotypes in gene-based tests

* 0.36.3.1 (February-04-2020):
** Note: in v0.36.3.1, uses SPAtest 3.0.2

* 0.36.3 (January-05-2020):
** Note: in v0.36.3, an option IsOutputBETASEinBurdenTest in step 2 is added to output effect sizes for burden tests
Bugs fixed: the header in output files from conditional analysis in gene or reigon-based tests is corrected.  

* 0.36.2 (November-23-2019):
** Note: in v0.36.2, users can specify customized weights for markers in gene- or region-based tests by adding a weight for each marker in the group file

Bugs fixed: 1. The option weights.beta.common is not fully correctly developed, so we make weights.beta.common equal to weights.beta.rare for now. 2. Instead of output NA for SKAT-O p values when the function SKAT:::Met_SKAT_Get_Pvalue failed, output 2*min(SKAT p, Burden p, 0.5).


* 0.36.1 (November-12-2019): 

** Note: in v0.36.1, plain text dosage files are no longer allowed as input in step 2 to get rid of the dependence of the boost_iostream library

Bugs fixed: 1. fixed the freq calculation for mean impute for missing genotypes in  plinkFile 2. Diagonal elements of GRM are now estimated using markers in plinkFile with MAF >= minMAFforGRM 3. Conditional analysis for gene- or region-based test for binary traits is now accounting for case-control imbalance 4. plain dosage files are no longer supported for step 2 so no external boost_iostream library is needed

** minMAFforGRM is added as a parameter in step 0 and 1, so only markers in the plinkFile with MAF >= minMAFforGRM will be used for GRM
** weights.beta.rare, weights.beta.common, weightMAFcutoff, dosageZerodCutoff, IsOutputPvalueNAinGroupTestforBinary, IsAccountforCasecontrolImbalanceinGroupTest are added as new parameters in step 2

* 0.35.8.8: Fixes a matrix inversion issue in the null model and adds an optional argument for the null computation to remove binary covariates with low counts by juhis

* 0.35.8.8 (August-27-2019): Fixes a matrix inversion issue in the null model and adds an optional argument for the null computation to remove binary covariates with low counts by juhis

* 0.35.8.7 (August-15-2019): fixed the bug when there is no covariate specified, added an argument IsOutputNinCaseCtrl for step 2 to allow for output sample sizes in cases and controls for binary traits in the output file, fixed the out of boundary bug for LOCO

* 0.35.8.6 (August-13-2019): fixed the output bug when the genotype matrix has rank 1 for binary phenotypes and add an argument minMAFtoConstructGRM for step 0 and step 1 to allow users to specify the minumum MAF of markers used to construct GRM (default: 1%)

* 0.35.8.5 (June-29-2019): account for case control imbalance for binary traits in gene-based tests

* 0.35.8.3 (May-14-2019): fix a bug in the function getCovM_nopcg, which affected the conditional analysis for binary traits. Merge hyacz/master to use cget to manage superlu 

* 0.35.8.2 (April-16-2019): minor changes include fix error message, change MAC to MAF, add a line to check if the chomosome in plink file is numeric or not, add rsid to the header when input file is bgen

* 0.35.8.1: fix some errors in documentation and the warning message for case-control imbalance of binary traits when running SAIGE-GENE

* 0.35.8 merge changes in the master-gene branch to master

* 0.35.7 merge changes in 0.29.6 and 0.29.7 from master

* 0.35.6 merge in 0.29.5 from master

* 0.35.5 (fix a bug for updating predicted values in the model fit for binary traits. Added a function to create a sparse GRM only for a data set)

* 0.35.3 (this is a clean version for single-variance assoc tests, gene-based tests, and conditional analysis)

* 0.35.2.3 (this version works with the conditonal analysis and gene-based tests)

* 0.29.4.2 (this version works with R-3.5.1)

* 0.29.4 (this version works with R-3.4.4) update SAIGE as a bug for reading vcf and sav files was fixed in the savvy library

* 0.29.3.2 this version works with R-3.5.1

* 0.29.3: update SAIGE step 1 to use the updated R libary SPAtest 3.0.0

* 0.29.2: update SAIGE to use the updated R library SPAtest 3.0.0

* 0.29:
```
1. The colSums() error when there is no covariate has been fixed. 
2. BETA and Tstat are now for the alt allele for both quantitative and binary traits. Note that in version <= 0.26, for binary traits, BETA is for alt allele and for quantitative traits, BETA is for minor allele
3. Options for leave-one-chromosome-out (LOCO), cutoffs for the coefficient of variation (CV) for trace estimates and variance ratio estimates have been added, but these three options have not been extensively tested. CV is mainly for automatically determining whether the number of random markers selected is sufficient or not. If not, the number will be increased until the CV is lower than the specified cutoff.  
```
* 0.26: fixed a bug for the Tstat in the output
* 0.25: allow models with no covariates and GRM contruction using a large number of genetic markers (> 600,000)
* 0.24: centerVariable is no longer needed. QR transformation of the covariate matrix is automatically performed. Supports the dosage files in the VCF,BCF and SAV formats using the SAVVY library 
History
====

7 July 2016
----

* Updates to bgenix to handle UK biobank interim files and to avoid extra index tables in the index file.

21 March 2016
----

* BGEN spec and implementation updated to alter probability order for unphased data when the number of alleles (K) or the ploidy is greater than two.
This order now better matches the order of VCF GP fields and as a simple enumeration scheme.

10 Nov 2015
----

Major changes in revision ff11254f9505:

1. I've implemented two new tools
    - cat-bgen, which can be used to concatenate BGEN files.
    - bgenix, which can be used to index BGEN files and efficiently retrieve specified data.

2. For this purpose I've imported several extra pieces of code
    - appcontext/ and db/ sublibs from qctool
    - sqlite3 3.9.2
    - boost 1.55.0

Note: these changes were erroneously applied first to the master branch (they were intended for default first).

6 Nov 2015
----
Major changes in revision 392429affc42:

1. I’ve changed the behaviour of BGEN v1.2 with respect to samples with missing data: they are now stored with dummy zero probabilities.  The spec is now in 'beta' which means I don’t have any other planned changes to make; unless major issues are uncovered this will be the final version of the format.

2. I’ve revamped the setter api of parse_probability_data somewhat.  It is documented in the code and here [on the wiki](https://bitbucket.org/gavinband/bgen/wiki/The_Setter_API).  The main breaking changes are:
- Renamed operator() to set_value(), and given it an index argument; I think these make the API more consistent.
- Added an initial ploidy argument to set_number_of_entries() as requested.  (The type of data - phased or unphased - is already reported in the order_type argument so I don’t think another argument is needed).
- Added two new method calls, which are optional: set_min_max_ploidy() (useful for setting storage) and finalise().  See the docs for info.

3. I’ve also got rid of the max_id_size option to write_snp_identifying_data().  (This is now not needed because writing BGEN v1.0 files is no longer supported.)

4. I’ve also added some test code (using the [catch framework](https://github.com/philsquared/catch), which seems pretty good).  Tests are not exhaustive but hopefully a start.
:q
5. I've removed some code warnings - thanks to Robert V. Baron of [Mega2](https://watson.hgen.pitt.edu/docs/mega2_html/mega2.html) for testing this code.

23 Sep 2015
----
First version, based on qctool implementation.
# BGEN reference implementation

This repository contains a reference implementation of
the [BGEN format](http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format_v1.2.html), written in C++.
The library can be used as the basis for BGEN support in other software, or as a reference for developers writing 
their own implementations of the BGEN format.

### What's included?
This repository contains the library itself, a set of [example data files](example/),
and a number of example programs (e.g. [bgen_to_vcf](example/bgen_to_vcf.cpp)) that demonstrate the use of the library API.

In addition, a number of utilities built using the library are also included in this repository:

* [bgenix](https://bitbucket.org/gavinband/bgen/wiki/bgenix) - a tool to index and efficiently retrieve subsets of a BGEN file. 
* [cat-bgen](https://bitbucket.org/gavinband/bgen/wiki/cat-bgen) - a tool to efficiently concatenate BGEN files.
* [edit-bgen](https://bitbucket.org/gavinband/bgen/wiki/edit-bgen) - a tool to edit BGEN file metadata.
* An R package called [rbgen](https://bitbucket.org/gavinband/bgen/wiki/rbgen) is also constructed in the build directory.  See the [rbgen wiki page](https://bitbucket.org/gavinband/bgen/wiki/rbgen) for more information on using this package.

### Citing BGEN

If you make use of the BGEN library, its tools or example programs, please cite:

Band, G. and Marchini, J., "*BGEN: a binary file format for imputed genotype and haplotype data*", bioArxiv bioRxiv 308296; doi: https://doi.org/10.1101/308296

Thanks!

### License
This BGEN implementation is released under the Boost Software License v1.0.  This is a relatively permissive open-source license that is compatible with many other open-source licenses.  See [this page](http://www.boost.org/users/license.html) and the file [LICENSE_1_0.txt](https://bitbucket.org/gavinband/bgen/src/tip/LICENSE_1_0.txt) for full details.

This repository also contains code from  the [sqlite](www.sqlite.org), [boost](www.boost.org), and [zstandard](http://www.zstd.net) libraries, which comes with their own respective licenses. (respectively, [public domain](http://www.sqlite.org/copyright.html), the boost software license, and the [BSD license](https://github.com/facebook/zstd/blob/dev/LICENSE)).  These libraries are not used in the core BGEN implementation, but may be used in the example programs provided.

---

### **!! Important note on the UK Biobank data**

The UK Biobank has released [imputed genotype data](http://www.ukbiobank.ac.uk/scientists-3/genetic-data/) for almost half a million individuals
in BGEN format, with accompanying bgenix index files.  The original release of this data (version 2) had an issue with
naming of the index files.  Please see [here](https://bitbucket.org/gavinband/bgen/wiki/Using the UK Biobank full release index files) for information on working around this.  The more recent version of this data (version 3) does not have this issue.

---

# Obtaining and installing BGEN

### In brief

The following commands (typed into a UNIX shell - the dollar symbol indicates the prompt, and shouldn't be typed in)
should perform a basic download and install of the BGEN library, example data and tools:

```bash
$ # get it
$ wget http://bitbucket.org/gavinband/bgen/get/master.tar.gz
$ cd bgen
$ # compile it
$ ./waf configure
$ ./waf
$ # test it
$ ./build/test/unit/test_bgen
$ ./build/apps/bgenix -g example/example.16bits.bgen -list
```

The following sections contains more information on this process.

### Download

A tarball of the latest master branch is available here: http://bitbucket.org/gavinband/bgen/get/master.tar.gz.

Alternatively, use mercurial to download the master branch as follows:
```sh
hg clone https://gavinband@bitbucket.org/gavinband/bgen -u master
```
(This command can take a while.)

Additionally, pre-built version of the bgen utilities may be available from [this page](http://www.well.ox.ac.uk/~gav/resources/).  **Note**: the recommended use is to download and compile bgenix for your platform; these binaries are provided for convenience in getting started quickly.

### Compilation

To compile the code, use the supplied waf build tool:
```sh
./waf configure
./waf
```
Results will appear under the `build/` directory.  

Note: a full build requires a compiler that supports C++11, e.g. gcc v4.7 or above.  To specify the compiler used, set the `CXX` environment variable during the configure step.  For example (if your shell is `bash`):
```
CXX=/path/to/g++ ./waf configure
./waf
```

The sqlite and zstd libraries are written in C; to specify the C compiler you can additionally add `CC=/path/to/gcc`.  We have tested compilation on gcc 4.9.3 and 5.4.0, and using clang, among others.

If you don't have access to a compiler with C++11 support, you can still build the core bgen implementation, but won't be able to build the applications or example programs.
See [the wiki](https://bitbucket.org/gavinband/bgen/wiki/Troubleshooting_compilation) for more information.

### Testing

BGEN's tests can be run by typing 
```sh
./build/test/test_bgen
```
or, for more recent versions:
```sh
./build/test/unit/test_bgen
```

If all goes well a message like `All tests passed` should be printed.

If you have [Robot Test Framework](http://robotframework.org/) installed, you can instead run the full suite of unit and functional tests like so:
```sh
./test/functional/run_tests.sh
```
Test results will be placed in the directory `build/test/functional/test-reports`.


### Trying an example

The example program `bgen_to_vcf` reads a bgen file (v1.1 or v1.2) and outputs it as a VCF file to stdout.  You can try running it
by typing
```sh
./build/example/bgen_to_vcf example/example.8bits.bgen
```
which should output vcf-formatted data to stdout.  We've provided further example bgen files in the `example/` subdirectory.

### Installation

The command
```sh
./waf install
```
will install the applications listed above into a specified system or user directory.  By default this is `/usr/local`.  To change it, specify the prefix at the configure step:
```sh
./waf configure --prefix=/path/to/installation/directory
./waf install
```
The programs listed above will be installed into a folder called `bin/` under the prefix dir, e.g. `bgenix` will be installed as `/path/to/installation/directory/bin/bgenix` etc.

Note that in many cases there's no need for installation; the executables are self-contained.  The install step simply copies them into the destination directory.

(The installation prefix need not be a system-wide directory.  For example, I typically specify an installation directory within my home dir, e.g. `~gav/projects/software/`.

### Branches

This repo follows the branch naming practice in which `master` represents the most up-to-date code considered in a 'releasable' state.  If you are interested in using bgen code in your own project, we therefore recommend cloning the `master` branch.  Code development takes place in the `default` branch and/or in feature branches branched from the `default` branch.  The command given above downloads the master branch, which is what most people will want.

### More information

See the [source code](https://bitbucket.org/gavinband/bgen/src), 
BGEN [releases](https://bitbucket.org/gavinband/bgen/wiki/Releases),
or the [Wiki](https://bitbucket.org/gavinband/bgen/wiki/Home) for more information.<body bgcolor="#ffffff" link="#0000ee" text="#000000" vlink="#551a8b" alink="#ff0000">

![C++ Boost](../../../boost.png)

# `hawick_circuits`

    template <typename Graph, typename Visitor, typename VertexIndexMap>
    void hawick_circuits(Graph const& graph, Visitor visitor, VertexIndexMap const& vim = get(vertex_index, graph));

    template <typename Graph, typename Visitor, typename VertexIndexMap>
    void hawick_unique_circuits(Graph const& graph, Visitor visitor, VertexIndexMap const& vim = get(vertex_index, graph));

Enumerate all the elementary circuits in a directed multigraph. Specifically,
self-loops and redundant circuits caused by parallel edges are enumerated too.
`hawick_unique_circuits` may be used if redundant circuits caused by parallel
edges are not desired.

The algorithm is described in detail in
<http://www.massey.ac.nz/~kahawick/cstn/013/cstn-013.pdf>.


### Where defined

[`#include <boost/graph/hawick_circuits.hpp>`](../../../boost/graph/hawick_circuits.hpp)


### Parameters

__IN:__ `Graph const& graph`

> The graph on which the algorithm is to be performed. It must be a model of
> the `VertexListGraph` and `AdjacencyGraph` concepts.

__IN:__ `Visitor visitor`

> The visitor that will be notified on each circuit found by the algorithm.
> The `visitor.cycle(circuit, graph)` expression must be valid, with `circuit`
> being a `const`-reference to a random access sequence of `vertex_descriptor`s.
>
> For example, if a circuit `u -> v -> w -> u` exists in the graph, the
> visitor will be called with a sequence consisting of `(u, v, w)`.

__IN:__ `VertexIndexMap const& vim = get(vertex_index, graph)`

> A model of the `ReadablePropertyMap` concept mapping each `vertex_descriptor`
> to an integer in the range `[0, num_vertices(graph))`. It defaults to using
> the vertex index map provided by the `graph`.


------------------------------------------------------------------------------
<div class="footer">
    &copy; 2013 Louis Dionne
</div>
 **Zstd**, short for Zstandard, is a fast lossless compression algorithm,
 targeting real-time compression scenarios at zlib-level and better compression ratios.

It is provided as an open-source BSD-licensed **C** library.
For other programming languages,
you can consult a list of known ports on [Zstandard homepage](http://www.zstd.net/#other-languages).

|Branch      |Status   |
|------------|---------|
|master      | [![Build Status](https://travis-ci.org/facebook/zstd.svg?branch=master)](https://travis-ci.org/facebook/zstd) |
|dev         | [![Build Status](https://travis-ci.org/facebook/zstd.svg?branch=dev)](https://travis-ci.org/facebook/zstd) |

As a reference, several fast compression algorithms were tested and compared on a Core i7-3930K CPU @ 4.5GHz, using [lzbench], an open-source in-memory benchmark by @inikep compiled with GCC 5.4.0, with the [Silesia compression corpus].

[lzbench]: https://github.com/inikep/lzbench
[Silesia compression corpus]: http://sun.aei.polsl.pl/~sdeor/index.php?page=silesia


|Name             | Ratio | C.speed | D.speed |
|-----------------|-------|--------:|--------:|
|                 |       |   MB/s  |  MB/s   |
|**zstd 0.8.2 -1**|**2.877**|**330**| **940** |
| [zlib] 1.2.8 -1 | 2.730 |    95   |   360   |
| brotli 0.4 -0   | 2.708 |   320   |   375   |
| QuickLZ 1.5     | 2.237 |   510   |   605   |
| LZO 2.09        | 2.106 |   610   |   870   |
| [LZ4] r131      | 2.101 |   620   |  3100   |
| Snappy 1.1.3    | 2.091 |   480   |  1600   |
| LZF 3.6         | 2.077 |   375   |   790   |

[zlib]:http://www.zlib.net/
[LZ4]: http://www.lz4.org/

Zstd can also offer stronger compression ratios at the cost of compression speed.
Speed vs Compression trade-off is configurable by small increments. Decompression speed is preserved and remains roughly the same at all settings, a property shared by most LZ compression algorithms, such as [zlib] or lzma.

The following tests were run on a Core i7-3930K CPU @ 4.5GHz, using [lzbench], an open-source in-memory benchmark by @inikep compiled with GCC 5.2.1, on the [Silesia compression corpus].

Compression Speed vs Ratio | Decompression Speed
---------------------------|--------------------
![Compression Speed vs Ratio](images/Cspeed4.png "Compression Speed vs Ratio") | ![Decompression Speed](images/Dspeed4.png "Decompression Speed")

Several algorithms can produce higher compression ratios, but at slower speeds, falling outside of the graph.
For a larger picture including very slow modes, [click on this link](images/DCspeed5.png) .


### The case for Small Data compression

Previous charts provide results applicable to typical file and stream scenarios (several MB). Small data comes with different perspectives. The smaller the amount of data to compress, the more difficult it is to achieve any significant compression.

This problem is common to many compression algorithms. The reason is, compression algorithms learn from past data how to compress future data. But at the beginning of a new file, there is no "past" to build upon.

To solve this situation, Zstd offers a __training mode__, which can be used to tune the algorithm for a selected type of data, by providing it with a few samples. The result of the training is stored in a file called "dictionary", which can be loaded before compression and decompression. Using this dictionary, the compression ratio achievable on small data improves dramatically:

![Compressing Small Data](images/smallData.png "Compressing Small Data")

These compression gains are achieved while simultaneously providing faster compression and decompression speeds.

Dictionary works if there is some correlation in a family of small data (there is no _universal dictionary_).
Hence, deploying one dictionary per type of data will provide the greatest benefits. Dictionary gains are mostly effective in the first few KB. Then, the compression algorithm will rely more and more on previously decoded content to compress the rest of the file.

#### Dictionary compression How To :

1) Create the dictionary

`zstd --train FullPathToTrainingSet/* -o dictionaryName`

2) Compress with dictionary

`zstd FILE -D dictionaryName`

3) Decompress with dictionary

`zstd --decompress FILE.zst -D dictionaryName`

### Build

Once you have the repository cloned, there are multiple ways provided to build Zstandard.

#### Makefile

If your system is compatible with a standard `make` (or `gmake`) binary generator,
you can simply run it at the root directory.
It will generate `zstd` within root directory.

Other available options include :
- `make install` : create and install zstd binary, library and man page
- `make test` : create and run `zstd` and test tools on local platform

#### cmake

A `cmake` project generator is provided within `build/cmake`.
It can generate Makefiles or other build scripts
to create `zstd` binary, and `libzstd` dynamic and static libraries.

#### Visual (Windows)

Going into `build` directory, you will find additional possibilities :
- Projects for Visual Studio 2005, 2008 and 2010
  + VS2010 project is compatible with VS2012, VS2013 and VS2015
- Automated build scripts for Visual compiler by @KrzysFR , in `build/VS_scripts`,
  which will build `zstd` cli and `libzstd` library without any need to open Visual Studio solution.


### Status

Zstandard is currently deployed within Facebook. It is used daily to compress and decompress very large amounts of data in multiple formats and use cases.
Zstandard is considered safe for production environments.

### License

Zstandard is [BSD-licensed](LICENSE). We also provide an [additional patent grant](PATENTS).

### Contributing

The "dev" branch is the one where all contributions will be merged before reaching "master".
If you plan to propose a patch, please commit into the "dev" branch or its own feature branch.
Direct commit to "master" are not permitted.
For more information, please read [CONTRIBUTING](CONTRIBUTING.md).

### Miscellaneous

Zstd entropy stage is provided by [Huff0 and FSE, from Finite State Entropy library](https://github.com/Cyan4973/FiniteStateEntropy).
Zstandard Compression Format
============================

### Notices

Copyright (c) 2016 Yann Collet

Permission is granted to copy and distribute this document
for any purpose and without charge,
including translations into other languages
and incorporation into compilations,
provided that the copyright notice and this notice are preserved,
and that any substantive changes or deletions from the original
are clearly marked.
Distribution of this document is unlimited.

### Version

0.2.2 (14/09/16)


Introduction
------------

The purpose of this document is to define a lossless compressed data format,
that is independent of CPU type, operating system,
file system and character set, suitable for
file compression, pipe and streaming compression,
using the [Zstandard algorithm](http://www.zstandard.org).

The data can be produced or consumed,
even for an arbitrarily long sequentially presented input data stream,
using only an a priori bounded amount of intermediate storage,
and hence can be used in data communications.
The format uses the Zstandard compression method,
and optional [xxHash-64 checksum method](http://www.xxhash.org),
for detection of data corruption.

The data format defined by this specification
does not attempt to allow random access to compressed data.

This specification is intended for use by implementers of software
to compress data into Zstandard format and/or decompress data from Zstandard format.
The text of the specification assumes a basic background in programming
at the level of bits and other primitive data representations.

Unless otherwise indicated below,
a compliant compressor must produce data sets
that conform to the specifications presented here.
It doesn’t need to support all options though.

A compliant decompressor must be able to decompress
at least one working set of parameters
that conforms to the specifications presented here.
It may also ignore informative fields, such as checksum.
Whenever it does not support a parameter defined in the compressed stream,
it must produce a non-ambiguous error code and associated error message
explaining which parameter is unsupported.


Overall conventions
-----------
In this document:
- square brackets i.e. `[` and `]` are used to indicate optional fields or parameters.
- a naming convention for identifiers is `Mixed_Case_With_Underscores`

Definitions
-----------
A content compressed by Zstandard is transformed into a Zstandard __frame__.
Multiple frames can be appended into a single file or stream.
A frame is totally independent, has a defined beginning and end,
and a set of parameters which tells the decoder how to decompress it.

A frame encapsulates one or multiple __blocks__.
Each block can be compressed or not,
and has a guaranteed maximum content size, which depends on frame parameters.
Unlike frames, each block depends on previous blocks for proper decoding.
However, each block can be decompressed without waiting for its successor,
allowing streaming operations.


Frame Concatenation
-------------------

In some circumstances, it may be required to append multiple frames,
for example in order to add new data to an existing compressed file
without re-framing it.

In such case, each frame brings its own set of descriptor flags.
Each frame is considered independent.
The only relation between frames is their sequential order.

The ability to decode multiple concatenated frames
within a single stream or file is left outside of this specification.
As an example, the reference `zstd` command line utility is able
to decode all concatenated frames in their sequential order,
delivering the final decompressed result as if it was a single content.


Skippable Frames
----------------

| `Magic_Number` | `Frame_Size` | `User_Data` |
|:--------------:|:------------:|:-----------:|
|   4 bytes      |  4 bytes     |   n bytes   |

Skippable frames allow the insertion of user-defined data
into a flow of concatenated frames.
Its design is pretty straightforward,
with the sole objective to allow the decoder to quickly skip
over user-defined data and continue decoding.

Skippable frames defined in this specification are compatible with [LZ4] ones.

[LZ4]:http://www.lz4.org

__`Magic_Number`__

4 Bytes, little-endian format.
Value : 0x184D2A5X, which means any value from 0x184D2A50 to 0x184D2A5F.
All 16 values are valid to identify a skippable frame.

__`Frame_Size`__

This is the size, in bytes, of the following `User_Data`
(without including the magic number nor the size field itself).
This field is represented using 4 Bytes, little-endian format, unsigned 32-bits.
This means `User_Data` can’t be bigger than (2^32-1) bytes.

__`User_Data`__

The `User_Data` can be anything. Data will just be skipped by the decoder.



General Structure of Zstandard Frame format
-------------------------------------------
The structure of a single Zstandard frame is following:

| `Magic_Number` | `Frame_Header` |`Data_Block`| [More data blocks] | [`Content_Checksum`] |
|:--------------:|:--------------:|:----------:| ------------------ |:--------------------:|
| 4 bytes        |  2-14 bytes    | n bytes    |                    |   0-4 bytes          |

__`Magic_Number`__

4 Bytes, little-endian format.
Value : 0xFD2FB528

__`Frame_Header`__

2 to 14 Bytes, detailed in [next part](#the-structure-of-frame_header).

__`Data_Block`__

Detailed in [next chapter](#the-structure-of-data_block).
That’s where compressed data is stored.

__`Content_Checksum`__

An optional 32-bit checksum, only present if `Content_Checksum_flag` is set.
The content checksum is the result
of [xxh64() hash function](http://www.xxhash.org)
digesting the original (decoded) data as input, and a seed of zero.
The low 4 bytes of the checksum are stored in little endian format.


The structure of `Frame_Header`
-------------------------------
The `Frame_Header` has a variable size, which uses a minimum of 2 bytes,
and up to 14 bytes depending on optional parameters.
The structure of `Frame_Header` is following:

| `Frame_Header_Descriptor` | [`Window_Descriptor`] | [`Dictionary_ID`] | [`Frame_Content_Size`] |
| ------------------------- | --------------------- | ----------------- | ---------------------- |
| 1 byte                    | 0-1 byte              | 0-4 bytes         | 0-8 bytes              |

### `Frame_Header_Descriptor`

The first header's byte is called the `Frame_Header_Descriptor`.
It tells which other fields are present.
Decoding this byte is enough to tell the size of `Frame_Header`.

| Bit number | Field name                |
| ---------- | ----------                |
| 7-6        | `Frame_Content_Size_flag` |
| 5          | `Single_Segment_flag`     |
| 4          | `Unused_bit`              |
| 3          | `Reserved_bit`            |
| 2          | `Content_Checksum_flag`   |
| 1-0        | `Dictionary_ID_flag`      |

In this table, bit 7 is highest bit, while bit 0 is lowest.

__`Frame_Content_Size_flag`__

This is a 2-bits flag (`= Frame_Header_Descriptor >> 6`),
specifying if decompressed data size is provided within the header.
The `Flag_Value` can be converted into `Field_Size`,
which is the number of bytes used by `Frame_Content_Size`
according to the following table:

|`Flag_Value`|    0   |  1  |  2  |  3  |
| ---------- | ------ | --- | --- | --- |
|`Field_Size`| 0 or 1 |  2  |  4  |  8  |

When `Flag_Value` is `0`, `Field_Size` depends on `Single_Segment_flag` :
if `Single_Segment_flag` is set, `Field_Size` is 1.
Otherwise, `Field_Size` is 0 (content size not provided).

__`Single_Segment_flag`__

If this flag is set,
data must be regenerated within a single continuous memory segment.

In this case, `Frame_Content_Size` is necessarily present,
but `Window_Descriptor` byte is skipped.
As a consequence, the decoder must allocate a memory segment
of size equal or bigger than `Frame_Content_Size`.

In order to preserve the decoder from unreasonable memory requirement,
a decoder can reject a compressed frame
which requests a memory size beyond decoder's authorized range.

For broader compatibility, decoders are recommended to support
memory sizes of at least 8 MB.
This is just a recommendation,
each decoder is free to support higher or lower limits,
depending on local limitations.

__`Unused_bit`__

The value of this bit should be set to zero.
A decoder compliant with this specification version shall not interpret it.
It might be used in a future version,
to signal a property which is not mandatory to properly decode the frame.

__`Reserved_bit`__

This bit is reserved for some future feature.
Its value _must be zero_.
A decoder compliant with this specification version must ensure it is not set.
This bit may be used in a future revision,
to signal a feature that must be interpreted to decode the frame correctly.

__`Content_Checksum_flag`__

If this flag is set, a 32-bits `Content_Checksum` will be present at frame's end.
See `Content_Checksum` paragraph.

__`Dictionary_ID_flag`__

This is a 2-bits flag (`= FHD & 3`),
telling if a dictionary ID is provided within the header.
It also specifies the size of this field as `Field_Size`.

|`Flag_Value`|  0  |  1  |  2  |  3  |
| ---------- | --- | --- | --- | --- |
|`Field_Size`|  0  |  1  |  2  |  4  |

### `Window_Descriptor`

Provides guarantees on maximum back-reference distance
that will be used within compressed data.
This information is important for decoders to allocate enough memory.

The `Window_Descriptor` byte is optional. It is absent when `Single_Segment_flag` is set.
In this case, the maximum back-reference distance is the content size itself,
which can be any value from 1 to 2^64-1 bytes (16 EB).

| Bit numbers |     7-3    |     0-2    |
| ----------- | ---------- | ---------- |
| Field name  | `Exponent` | `Mantissa` |

Maximum distance is given by the following formulas :
```
windowLog = 10 + Exponent;
windowBase = 1 << windowLog;
windowAdd = (windowBase / 8) * Mantissa;
Window_Size = windowBase + windowAdd;
```
The minimum window size is 1 KB.
The maximum size is `15*(1<<38)` bytes, which is 1.875 TB.

To properly decode compressed data,
a decoder will need to allocate a buffer of at least `Window_Size` bytes.

In order to preserve decoder from unreasonable memory requirements,
a decoder can refuse a compressed frame
which requests a memory size beyond decoder's authorized range.

For improved interoperability,
decoders are recommended to be compatible with window sizes of 8 MB,
and encoders are recommended to not request more than 8 MB.
It's merely a recommendation though,
decoders are free to support larger or lower limits,
depending on local limitations.

### `Dictionary_ID`

This is a variable size field, which contains
the ID of the dictionary required to properly decode the frame.
Note that this field is optional. When it's not present,
it's up to the caller to make sure it uses the correct dictionary.
Format is little-endian.

Field size depends on `Dictionary_ID_flag`.
1 byte can represent an ID 0-255.
2 bytes can represent an ID 0-65535.
4 bytes can represent an ID 0-4294967295.

It's allowed to represent a small ID (for example `13`)
with a large 4-bytes dictionary ID, losing some compacity in the process.

_Reserved ranges :_
If the frame is going to be distributed in a private environment,
any dictionary ID can be used.
However, for public distribution of compressed frames using a dictionary,
the following ranges are reserved for future use and should not be used :
- low range : 1 - 32767
- high range : >= (2^31)


### `Frame_Content_Size`

This is the original (uncompressed) size. This information is optional.
The `Field_Size` is provided according to value of `Frame_Content_Size_flag`.
The `Field_Size` can be equal to 0 (not present), 1, 2, 4 or 8 bytes.
Format is little-endian.

| `Field_Size` |    Range   |
| ------------ | ---------- |
|      1       |   0 - 255  |
|      2       | 256 - 65791|
|      4       | 0 - 2^32-1 |
|      8       | 0 - 2^64-1 |

When `Field_Size` is 1, 4 or 8 bytes, the value is read directly.
When `Field_Size` is 2, _the offset of 256 is added_.
It's allowed to represent a small size (for example `18`) using any compatible variant.


The structure of `Data_Block`
-----------------------------
The structure of `Data_Block` is following:

| `Last_Block` | `Block_Type` | `Block_Size` | `Block_Content` |
|:------------:|:------------:|:------------:|:---------------:|
|   1 bit      |  2 bits      |  21 bits     |  n bytes        |

The block header (`Last_Block`, `Block_Type`, and `Block_Size`) uses 3-bytes.

__`Last_Block`__

The lowest bit signals if this block is the last one.
Frame ends right after this block.
It may be followed by an optional `Content_Checksum` .

__`Block_Type` and `Block_Size`__

The next 2 bits represent the `Block_Type`,
while the remaining 21 bits represent the `Block_Size`.
Format is __little-endian__.

There are 4 block types :

|    Value     |      0      |     1       |  2                 |    3      |
| ------------ | ----------- | ----------- | ------------------ | --------- |
| `Block_Type` | `Raw_Block` | `RLE_Block` | `Compressed_Block` | `Reserved`|

- `Raw_Block` - this is an uncompressed block.
  `Block_Size` is the number of bytes to read and copy.
- `RLE_Block` - this is a single byte, repeated N times.
  In which case, `Block_Size` is the size to regenerate,
  while the "compressed" block is just 1 byte (the byte to repeat).
- `Compressed_Block` - this is a [Zstandard compressed block](#the-format-of-compressed_block),
  detailed in another section of this specification.
  `Block_Size` is the compressed size.
  Decompressed size is unknown,
  but its maximum possible value is guaranteed (see below)
- `Reserved` - this is not a block.
  This value cannot be used with current version of this specification.

Block sizes must respect a few rules :
- In compressed mode, compressed size if always strictly `< decompressed size`.
- Block decompressed size is always <= maximum back-reference distance .
- Block decompressed size is always <= 128 KB


__`Block_Content`__

The `Block_Content` is where the actual data to decode stands.
It might be compressed or not, depending on previous field indications.
A data block is not necessarily "full" :
since an arbitrary “flush” may happen anytime,
block decompressed content can be any size,
up to `Block_Maximum_Decompressed_Size`, which is the smallest of :
- Maximum back-reference distance
- 128 KB



The format of `Compressed_Block`
--------------------------------
The size of `Compressed_Block` must be provided using `Block_Size` field from `Data_Block`.
The `Compressed_Block` has a guaranteed maximum regenerated size,
in order to properly allocate destination buffer.
See [`Data_Block`](#the-structure-of-data_block) for more details.

A compressed block consists of 2 sections :
- [`Literals_Section`](#literals_section)
- [`Sequences_Section`](#sequences_section)

### Prerequisites
To decode a compressed block, the following elements are necessary :
- Previous decoded blocks, up to a distance of `Window_Size`,
  or all previous blocks when `Single_Segment_flag` is set.
- List of "recent offsets" from previous compressed block.
- Decoding tables of previous compressed block for each symbol type
  (literals, literals lengths, match lengths, offsets).


### `Literals_Section`

During sequence phase, literals will be entangled with match copy operations.
All literals are regrouped in the first part of the block.
They can be decoded first, and then copied during sequence operations,
or they can be decoded on the flow, as needed by sequence commands.

| `Literals_Section_Header` | [`Huffman_Tree_Description`] | Stream1 | [Stream2] | [Stream3] | [Stream4] |
| ------------------------- | ---------------------------- | ------- | --------- | --------- | --------- |

Literals can be stored uncompressed or compressed using Huffman prefix codes.
When compressed, an optional tree description can be present,
followed by 1 or 4 streams.


#### `Literals_Section_Header`

Header is in charge of describing how literals are packed.
It's a byte-aligned variable-size bitfield, ranging from 1 to 5 bytes,
using little-endian convention.

| `Literals_Block_Type` | `Size_Format` | `Regenerated_Size` | [`Compressed_Size`] |
| --------------------- | ------------- | ------------------ | ----------------- |
|   2 bits              |  1 - 2 bits   |    5 - 20 bits     |    0 - 18 bits    |

In this representation, bits on the left are smallest bits.

__`Literals_Block_Type`__

This field uses 2 lowest bits of first byte, describing 4 different block types :

| `Literals_Block_Type`         | Value |
| ----------------------------- | ----- |
| `Raw_Literals_Block`          |   0   |
| `RLE_Literals_Block`          |   1   |
| `Compressed_Literals_Block`   |   2   |
| `Repeat_Stats_Literals_Block` |   3   |

- `Raw_Literals_Block` - Literals are stored uncompressed.
- `RLE_Literals_Block` - Literals consist of a single byte value repeated N times.
- `Compressed_Literals_Block` - This is a standard Huffman-compressed block,
        starting with a Huffman tree description.
        See details below.
- `Repeat_Stats_Literals_Block` - This is a Huffman-compressed block,
        using Huffman tree _from previous Huffman-compressed literals block_.
        Huffman tree description will be skipped.

__`Size_Format`__

`Size_Format` is divided into 2 families :

- For `Compressed_Block`, it requires to decode both `Compressed_Size`
  and `Regenerated_Size` (the decompressed size). It will also decode the number of streams.
- For `Raw_Literals_Block` and `RLE_Literals_Block` it's enough to decode `Regenerated_Size`.

For values spanning several bytes, convention is little-endian.

__`Size_Format` for `Raw_Literals_Block` and `RLE_Literals_Block`__ :

- Value x0 : `Regenerated_Size` uses 5 bits (0-31).
               `Literals_Section_Header` has 1 byte.
               `Regenerated_Size = Header[0]>>3`
- Value 01 : `Regenerated_Size` uses 12 bits (0-4095).
               `Literals_Section_Header` has 2 bytes.
               `Regenerated_Size = (Header[0]>>4) + (Header[1]<<4)`
- Value 11 : `Regenerated_Size` uses 20 bits (0-1048575).
               `Literals_Section_Header` has 3 bytes.
               `Regenerated_Size = (Header[0]>>4) + (Header[1]<<4) + (Header[2]<<12)`

Note : it's allowed to represent a short value (for example `13`)
using a long format, accepting the increased compressed data size.

__`Size_Format` for `Compressed_Literals_Block` and `Repeat_Stats_Literals_Block`__ :

- Value 00 : _A single stream_.
               Both `Compressed_Size` and `Regenerated_Size` use 10 bits (0-1023).
               `Literals_Section_Header` has 3 bytes.
- Value 01 : 4 streams.
               Both `Compressed_Size` and `Regenerated_Size` use 10 bits (0-1023).
               `Literals_Section_Header` has 3 bytes.
- Value 10 : 4 streams.
               Both `Compressed_Size` and `Regenerated_Size` use 14 bits (0-16383).
               `Literals_Section_Header` has 4 bytes.
- Value 11 : 4 streams.
               Both `Compressed_Size` and `Regenerated_Size` use 18 bits (0-262143).
               `Literals_Section_Header` has 5 bytes.

Both `Compressed_Size` and `Regenerated_Size` fields follow little-endian convention.


#### `Huffman_Tree_Description`

This section is only present when `Literals_Block_Type` type is `Compressed_Literals_Block` (`2`).

Prefix coding represents symbols from an a priori known alphabet
by bit sequences (codewords), one codeword for each symbol,
in a manner such that different symbols may be represented
by bit sequences of different lengths,
but a parser can always parse an encoded string
unambiguously symbol-by-symbol.

Given an alphabet with known symbol frequencies,
the Huffman algorithm allows the construction of an optimal prefix code
using the fewest bits of any possible prefix codes for that alphabet.

Prefix code must not exceed a maximum code length.
More bits improve accuracy but cost more header size,
and require more memory or more complex decoding operations.
This specification limits maximum code length to 11 bits.


##### Representation

All literal values from zero (included) to last present one (excluded)
are represented by `Weight` with values from `0` to `Max_Number_of_Bits`.
Transformation from `Weight` to `Number_of_Bits` follows this formula :
```
Number_of_Bits = Weight ? (Max_Number_of_Bits + 1 - Weight) : 0
```
The last symbol's `Weight` is deduced from previously decoded ones,
by completing to the nearest power of 2.
This power of 2 gives `Max_Number_of_Bits`, the depth of the current tree.

__Example__ :
Let's presume the following Huffman tree must be described :

|     literal      |  0  |  1  |  2  |  3  |  4  |  5  |
| ---------------- | --- | --- | --- | --- | --- | --- |
| `Number_of_Bits` |  1  |  2  |  3  |  0  |  4  |  4  |

The tree depth is 4, since its smallest element uses 4 bits.
Value `5` will not be listed, nor will values above `5`.
Values from `0` to `4` will be listed using `Weight` instead of `Number_of_Bits`.
Weight formula is :
```
Weight = Number_of_Bits ? (Max_Number_of_Bits + 1 - Number_of_Bits) : 0
```
It gives the following serie of weights :

| `Weight` |  4  |  3  |  2  |  0  |  1  |
| -------- | --- | --- | --- | --- | --- |
| literal  |  0  |  1  |  2  |  3  |  4  |

The decoder will do the inverse operation :
having collected weights of literals from `0` to `4`,
it knows the last literal, `5`, is present with a non-zero weight.
The weight of `5` can be deducted by joining to the nearest power of 2.
Sum of `2^(Weight-1)` (excluding 0) is :
`8 + 4 + 2 + 0 + 1 = 15`.
Nearest power of 2 is 16.
Therefore, `Max_Number_of_Bits = 4` and `Weight[5] = 1`.

##### Huffman Tree header

This is a single byte value (0-255),
which tells how to decode the list of weights.

- if `headerByte` >= 128 : this is a direct representation,
  where each `Weight` is written directly as a 4 bits field (0-15).
  The full representation occupies `((Number_of_Symbols+1)/2)` bytes,
  meaning it uses a last full byte even if `Number_of_Symbols` is odd.
  `Number_of_Symbols = headerByte - 127`.
  Note that maximum `Number_of_Symbols` is 255-127 = 128.
  A larger serie must necessarily use FSE compression.

- if `headerByte` < 128 :
  the serie of weights is compressed by FSE.
  The length of the FSE-compressed serie is equal to `headerByte` (0-127).

##### Finite State Entropy (FSE) compression of Huffman weights

The serie of weights is compressed using FSE compression.
It's a single bitstream with 2 interleaved states,
sharing a single distribution table.

To decode an FSE bitstream, it is necessary to know its compressed size.
Compressed size is provided by `headerByte`.
It's also necessary to know its _maximum possible_ decompressed size,
which is `255`, since literal values span from `0` to `255`,
and last symbol value is not represented.

An FSE bitstream starts by a header, describing probabilities distribution.
It will create a Decoding Table.
Table must be pre-allocated, which requires to support a maximum accuracy.
For a list of Huffman weights, maximum accuracy is 7 bits.

FSE header is [described in relevant chapter](#fse-distribution-table--condensed-format),
and so is [FSE bitstream](#bitstream).
The main difference is that Huffman header compression uses 2 states,
which share the same FSE distribution table.
Bitstream contains only FSE symbols (no interleaved "raw bitfields").
The number of symbols to decode is discovered
by tracking bitStream overflow condition.
When both states have overflowed the bitstream, end is reached.


##### Conversion from weights to Huffman prefix codes

All present symbols shall now have a `Weight` value.
It is possible to transform weights into Number_of_Bits, using this formula:
```
Number_of_Bits = Number_of_Bits ? Max_Number_of_Bits + 1 - Weight : 0
```
Symbols are sorted by `Weight`. Within same `Weight`, symbols keep natural order.
Symbols with a `Weight` of zero are removed.
Then, starting from lowest weight, prefix codes are distributed in order.

__Example__ :
Let's presume the following list of weights has been decoded :

| Literal  |  0  |  1  |  2  |  3  |  4  |  5  |
| -------- | --- | --- | --- | --- | --- | --- |
| `Weight` |  4  |  3  |  2  |  0  |  1  |  1  |

Sorted by weight and then natural order,
it gives the following distribution :

| Literal          |  3  |  4  |  5  |  2  |  1  |   0  |
| ---------------- | --- | --- | --- | --- | --- | ---- |
| `Weight`         |  0  |  1  |  1  |  2  |  3  |   4  |
| `Number_of_Bits` |  0  |  4  |  4  |  3  |  2  |   1  |
| prefix codes     | N/A | 0000| 0001| 001 | 01  |   1  |


#### The content of Huffman-compressed literal stream

##### Bitstreams sizes

As seen in a previous paragraph,
there are 2 types of Huffman-compressed literals :
a single stream and 4 streams.

Encoding using 4 streams is useful for CPU with multiple execution units and out-of-order operations.
Since each stream can be decoded independently,
it's possible to decode them up to 4x faster than a single stream,
presuming the CPU has enough parallelism available.

For single stream, header provides both the compressed and regenerated size.
For 4 streams though,
header only provides compressed and regenerated size of all 4 streams combined.
In order to properly decode the 4 streams,
it's necessary to know the compressed and regenerated size of each stream.

Regenerated size of each stream can be calculated by `(totalSize+3)/4`,
except for last one, which can be up to 3 bytes smaller, to reach `totalSize`.

Compressed size is provided explicitly : in the 4-streams variant,
bitstreams are preceded by 3 unsigned little-endian 16-bits values.
Each value represents the compressed size of one stream, in order.
The last stream size is deducted from total compressed size
and from previously decoded stream sizes :

`stream4CSize = totalCSize - 6 - stream1CSize - stream2CSize - stream3CSize`.


##### Bitstreams read and decode

Each bitstream must be read _backward_,
that is starting from the end down to the beginning.
Therefore it's necessary to know the size of each bitstream.

It's also necessary to know exactly which _bit_ is the latest.
This is detected by a final bit flag :
the highest bit of latest byte is a final-bit-flag.
Consequently, a last byte of `0` is not possible.
And the final-bit-flag itself is not part of the useful bitstream.
Hence, the last byte contains between 0 and 7 useful bits.

Starting from the end,
it's possible to read the bitstream in a little-endian fashion,
keeping track of already used bits.

Reading the last `Max_Number_of_Bits` bits,
it's then possible to compare extracted value to decoding table,
determining the symbol to decode and number of bits to discard.

The process continues up to reading the required number of symbols per stream.
If a bitstream is not entirely and exactly consumed,
hence reaching exactly its beginning position with _all_ bits consumed,
the decoding process is considered faulty.


### `Sequences_Section`

A compressed block is a succession of _sequences_ .
A sequence is a literal copy command, followed by a match copy command.
A literal copy command specifies a length.
It is the number of bytes to be copied (or extracted) from the literal section.
A match copy command specifies an offset and a length.
The offset gives the position to copy from,
which can be within a previous block.

When all _sequences_ are decoded,
if there is any literal left in the _literal section_,
these bytes are added at the end of the block.

The `Sequences_Section` regroup all symbols required to decode commands.
There are 3 symbol types : literals lengths, offsets and match lengths.
They are encoded together, interleaved, in a single _bitstream_.

The `Sequences_Section` starts by a header,
followed by optional probability tables for each symbol type,
followed by the bitstream.

| `Sequences_Section_Header` | [`Literals_Length_Table`] | [`Offset_Table`] | [`Match_Length_Table`] | bitStream |
| -------------------------- | ------------------------- | ---------------- | ---------------------- | --------- |

To decode the `Sequences_Section`, it's required to know its size.
This size is deducted from `blockSize - literalSectionSize`.


#### `Sequences_Section_Header`

Consists of 2 items:
- `Number_of_Sequences`
- Symbol compression modes

__`Number_of_Sequences`__

This is a variable size field using between 1 and 3 bytes.
Let's call its first byte `byte0`.
- `if (byte0 == 0)` : there are no sequences.
            The sequence section stops there.
            Regenerated content is defined entirely by literals section.
- `if (byte0 < 128)` : `Number_of_Sequences = byte0` . Uses 1 byte.
- `if (byte0 < 255)` : `Number_of_Sequences = ((byte0-128) << 8) + byte1` . Uses 2 bytes.
- `if (byte0 == 255)`: `Number_of_Sequences = byte1 + (byte2<<8) + 0x7F00` . Uses 3 bytes.

__Symbol compression modes__

This is a single byte, defining the compression mode of each symbol type.

|Bit number|   7-6                   |   5-4          |   3-2                |     1-0    |
| -------- | ----------------------- | -------------- | -------------------- | ---------- |
|Field name| `Literals_Lengths_Mode` | `Offsets_Mode` | `Match_Lengths_Mode` | `Reserved` |

The last field, `Reserved`, must be all-zeroes.

`Literals_Lengths_Mode`, `Offsets_Mode` and `Match_Lengths_Mode` define the `Compression_Mode` of
literals lengths, offsets, and match lengths respectively.

They follow the same enumeration :

|        Value       |         0         |      1     |           2           |       3       |
| ------------------ | ----------------- | ---------- | --------------------- | ------------- |
| `Compression_Mode` | `Predefined_Mode` | `RLE_Mode` | `FSE_Compressed_Mode` | `Repeat_Mode` |

- `Predefined_Mode` : uses a predefined distribution table.
- `RLE_Mode` : it's a single code, repeated `Number_of_Sequences` times.
- `Repeat_Mode` : re-use distribution table from previous compressed block.
- `FSE_Compressed_Mode` : standard FSE compression.
          A distribution table will be present.
          It will be described in [next part](#distribution-tables).

#### The codes for literals lengths, match lengths, and offsets.

Each symbol is a _code_ in its own context,
which specifies `Baseline` and `Number_of_Bits` to add.
_Codes_ are FSE compressed,
and interleaved with raw additional bits in the same bitstream.

##### Literals length codes

Literals length codes are values ranging from `0` to `35` included.
They define lengths from 0 to 131071 bytes.

| `Literals_Length_Code` |         0-15           |
| ---------------------- | ---------------------- |
| length                 | `Literals_Length_Code` |
| `Number_of_Bits`       |          0             |

| `Literals_Length_Code` |  16  |  17  |  18  |  19  |  20  |  21  |  22  |  23  |
| ---------------------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| `Baseline`             |  16  |  18  |  20  |  22  |  24  |  28  |  32  |  40  |
| `Number_of_Bits`       |   1  |   1  |   1  |   1  |   2  |   2  |   3  |   3  |

| `Literals_Length_Code` |  24  |  25  |  26  |  27  |  28  |  29  |  30  |  31  |
| ---------------------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| `Baseline`             |  48  |  64  |  128 |  256 |  512 | 1024 | 2048 | 4096 |
| `Number_of_Bits`       |   4  |   6  |   7  |   8  |   9  |  10  |  11  |  12  |

| `Literals_Length_Code` |  32  |  33  |  34  |  35  |
| ---------------------- | ---- | ---- | ---- | ---- |
| `Baseline`             | 8192 |16384 |32768 |65536 |
| `Number_of_Bits`       |  13  |  14  |  15  |  16  |

##### Default distribution for literals length codes

When `Compression_Mode` is `Predefined_Mode`,
a predefined distribution is used for FSE compression.

Below is its definition. It uses an accuracy of 6 bits (64 states).
```
short literalsLength_defaultDistribution[36] =
        { 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1,
          2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 1, 1, 1, 1, 1,
         -1,-1,-1,-1 };
```

##### Match length codes

Match length codes are values ranging from `0` to `52` included.
They define lengths from 3 to 131074 bytes.

| `Match_Length_Code` |         0-31            |
| ------------------- | ----------------------- |
| value               | `Match_Length_Code` + 3 |
| `Number_of_Bits`    |          0              |

| `Match_Length_Code` |  32  |  33  |  34  |  35  |  36  |  37  |  38  |  39  |
| ------------------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| `Baseline`          |  35  |  37  |  39  |  41  |  43  |  47  |  51  |  59  |
| `Number_of_Bits`    |   1  |   1  |   1  |   1  |   2  |   2  |   3  |   3  |

| `Match_Length_Code` |  40  |  41  |  42  |  43  |  44  |  45  |  46  |  47  |
| ------------------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| `Baseline`          |  67  |  83  |  99  |  131 |  258 |  514 | 1026 | 2050 |
| `Number_of_Bits`    |   4  |   4  |   5  |   7  |   8  |   9  |  10  |  11  |

| `Match_Length_Code` |  48  |  49  |  50  |  51  |  52  |
| ------------------- | ---- | ---- | ---- | ---- | ---- |
| `Baseline`          | 4098 | 8194 |16486 |32770 |65538 |
| `Number_of_Bits`    |  12  |  13  |  14  |  15  |  16  |

##### Default distribution for match length codes

When `Compression_Mode` is defined as `Predefined_Mode`,
a predefined distribution is used for FSE compression.

Below is its definition. It uses an accuracy of 6 bits (64 states).
```
short matchLengths_defaultDistribution[53] =
        { 1, 4, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,
         -1,-1,-1,-1,-1 };
```

##### Offset codes

Offset codes are values ranging from `0` to `N`.

A decoder is free to limit its maximum `N` supported.
Recommendation is to support at least up to `22`.
For information, at the time of this writing.
the reference decoder supports a maximum `N` value of `28` in 64-bits mode.

An offset code is also the number of additional bits to read,
and can be translated into an `Offset_Value` using the following formulas :

```
Offset_Value = (1 << offsetCode) + readNBits(offsetCode);
if (Offset_Value > 3) offset = Offset_Value - 3;
```
It means that maximum `Offset_Value` is `(2^(N+1))-1` and it supports back-reference distance up to `(2^(N+1))-4`
but is limited by [maximum back-reference distance](#window_descriptor).

`Offset_Value` from 1 to 3 are special : they define "repeat codes",
which means one of the previous offsets will be repeated.
They are sorted in recency order, with 1 meaning the most recent one.
See [Repeat offsets](#repeat-offsets) paragraph.


##### Default distribution for offset codes

When `Compression_Mode` is defined as `Predefined_Mode`,
a predefined distribution is used for FSE compression.

Below is its definition. It uses an accuracy of 5 bits (32 states),
and supports a maximum `N` of 28, allowing offset values up to 536,870,908 .

If any sequence in the compressed block requires an offset larger than this,
it's not possible to use the default distribution to represent it.

```
short offsetCodes_defaultDistribution[29] =
        { 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1 };
```

#### Distribution tables

Following the header, up to 3 distribution tables can be described.
When present, they are in this order :
- Literals lengths
- Offsets
- Match Lengths

The content to decode depends on their respective encoding mode :
- `Predefined_Mode` : no content. Use predefined distribution table.
- `RLE_Mode` : 1 byte. This is the only code to use across the whole compressed block.
- `FSE_Compressed_Mode` : A distribution table is present.
- `Repeat_Mode` : no content. Re-use distribution from previous compressed block.

##### FSE distribution table : condensed format

An FSE distribution table describes the probabilities of all symbols
from `0` to the last present one (included)
on a normalized scale of `1 << Accuracy_Log` .

It's a bitstream which is read forward, in little-endian fashion.
It's not necessary to know its exact size,
since it will be discovered and reported by the decoding process.

The bitstream starts by reporting on which scale it operates.
`Accuracy_Log = low4bits + 5`.
Note that maximum `Accuracy_Log` for literal and match lengths is `9`,
and for offsets is `8`. Higher values are considered errors.

Then follows each symbol value, from `0` to last present one.
The number of bits used by each field is variable.
It depends on :

- Remaining probabilities + 1 :
  __example__ :
  Presuming an `Accuracy_Log` of 8,
  and presuming 100 probabilities points have already been distributed,
  the decoder may read any value from `0` to `255 - 100 + 1 == 156` (included).
  Therefore, it must read `log2sup(156) == 8` bits.

- Value decoded : small values use 1 less bit :
  __example__ :
  Presuming values from 0 to 156 (included) are possible,
  255-156 = 99 values are remaining in an 8-bits field.
  They are used this way :
  first 99 values (hence from 0 to 98) use only 7 bits,
  values from 99 to 156 use 8 bits.
  This is achieved through this scheme :

  | Value read | Value decoded | Number of bits used |
  | ---------- | ------------- | ------------------- |
  |   0 -  98  |   0 -  98     |  7                  |
  |  99 - 127  |  99 - 127     |  8                  |
  | 128 - 226  |   0 -  98     |  7                  |
  | 227 - 255  | 128 - 156     |  8                  |

Symbols probabilities are read one by one, in order.

Probability is obtained from Value decoded by following formula :
`Proba = value - 1`

It means value `0` becomes negative probability `-1`.
`-1` is a special probability, which means "less than 1".
Its effect on distribution table is described in [next paragraph].
For the purpose of calculating cumulated distribution, it counts as one.

[next paragraph]:#fse-decoding--from-normalized-distribution-to-decoding-tables

When a symbol has a probability of `zero`,
it is followed by a 2-bits repeat flag.
This repeat flag tells how many probabilities of zeroes follow the current one.
It provides a number ranging from 0 to 3.
If it is a 3, another 2-bits repeat flag follows, and so on.

When last symbol reaches cumulated total of `1 << Accuracy_Log`,
decoding is complete.
If the last symbol makes cumulated total go above `1 << Accuracy_Log`,
distribution is considered corrupted.

Then the decoder can tell how many bytes were used in this process,
and how many symbols are present.
The bitstream consumes a round number of bytes.
Any remaining bit within the last byte is just unused.

##### FSE decoding : from normalized distribution to decoding tables

The distribution of normalized probabilities is enough
to create a unique decoding table.

It follows the following build rule :

The table has a size of `tableSize = 1 << Accuracy_Log`.
Each cell describes the symbol decoded,
and instructions to get the next state.

Symbols are scanned in their natural order for "less than 1" probabilities.
Symbols with this probability are being attributed a single cell,
starting from the end of the table.
These symbols define a full state reset, reading `Accuracy_Log` bits.

All remaining symbols are sorted in their natural order.
Starting from symbol `0` and table position `0`,
each symbol gets attributed as many cells as its probability.
Cell allocation is spreaded, not linear :
each successor position follow this rule :

```
position += (tableSize>>1) + (tableSize>>3) + 3;
position &= tableSize-1;
```

A position is skipped if already occupied,
typically by a "less than 1" probability symbol.

The result is a list of state values.
Each state will decode the current symbol.

To get the `Number_of_Bits` and `Baseline` required for next state,
it's first necessary to sort all states in their natural order.
The lower states will need 1 more bit than higher ones.

__Example__ :
Presuming a symbol has a probability of 5.
It receives 5 state values. States are sorted in natural order.

Next power of 2 is 8.
Space of probabilities is divided into 8 equal parts.
Presuming the `Accuracy_Log` is 7, it defines 128 states.
Divided by 8, each share is 16 large.

In order to reach 8, 8-5=3 lowest states will count "double",
taking shares twice larger,
requiring one more bit in the process.

Numbering starts from higher states using less bits.

| state order      |   0   |   1   |    2   |   3  |   4   |
| ---------------- | ----- | ----- | ------ | ---- | ----- |
| width            |  32   |  32   |   32   |  16  |  16   |
| `Number_of_Bits` |   5   |   5   |    5   |   4  |   4   |
| range number     |   2   |   4   |    6   |   0  |   1   |
| `Baseline`       |  32   |  64   |   96   |   0  |  16   |
| range            | 32-63 | 64-95 | 96-127 | 0-15 | 16-31 |

Next state is determined from current state
by reading the required `Number_of_Bits`, and adding the specified `Baseline`.


#### Bitstream

FSE bitstreams are read in reverse direction than written. In zstd,
the compressor writes bits forward into a block and the decompressor
must read the bitstream _backwards_.

To find the start of the bitstream it is therefore necessary to
know the offset of the last byte of the block which can be found
by counting `Block_Size` bytes after the block header.

After writing the last bit containing information, the compressor
writes a single `1`-bit and then fills the byte with 0-7 `0` bits of
padding. The last byte of the compressed bitstream cannot be `0` for
that reason.

When decompressing, the last byte containing the padding is the first
byte to read. The decompressor needs to skip 0-7 initial `0`-bits and
the first `1`-bit it occurs. Afterwards, the useful part of the bitstream
begins.

##### Starting states

The bitstream starts with initial state values,
each using the required number of bits in their respective _accuracy_,
decoded previously from their normalized distribution.

It starts by `Literals_Length_State`,
followed by `Offset_State`,
and finally `Match_Length_State`.

Reminder : always keep in mind that all values are read _backward_.

##### Decoding a sequence

A state gives a code.
A code provides `Baseline` and `Number_of_Bits` to add.
See [Symbol Decoding] section for details on each symbol.

Decoding starts by reading the `Number_of_Bits` required to decode `Offset`.
It then does the same for `Match_Length`,
and then for `Literals_Length`.

`Offset`, `Match_Length`, and `Literals_Length` define a sequence.
It starts by inserting the number of literals defined by `Literals_Length`,
then continue by copying `Match_Length` bytes from `currentPos - Offset`.

The next operation is to update states.
Using rules pre-calculated in the decoding tables,
`Literals_Length_State` is updated,
followed by `Match_Length_State`,
and then `Offset_State`.

This operation will be repeated `Number_of_Sequences` times.
At the end, the bitstream shall be entirely consumed,
otherwise bitstream is considered corrupted.

[Symbol Decoding]:#the-codes-for-literals-lengths-match-lengths-and-offsets

##### Repeat offsets

As seen in [Offset Codes], the first 3 values define a repeated offset and we will call them `Repeated_Offset1`, `Repeated_Offset2`, and `Repeated_Offset3`.
They are sorted in recency order, with `Repeated_Offset1` meaning "most recent one".

There is an exception though, when current sequence's literals length is `0`.
In which case, repeated offsets are "pushed by one",
so `Repeated_Offset1` becomes `Repeated_Offset2`, `Repeated_Offset2` becomes `Repeated_Offset3`,
and `Repeated_Offset3` becomes `Repeated_Offset1 - 1_byte`.

On first block, offset history is populated by the following values : 1, 4 and 8 (in order).

Then each block receives its start value from previous compressed block.
Note that non-compressed blocks are skipped,
they do not contribute to offset history.

[Offset Codes]: #offset-codes

###### Offset updates rules

New offset take the lead in offset history,
up to its previous place if it was already present.

It means that when `Repeated_Offset1` (most recent) is used, history is unmodified.
When `Repeated_Offset2` is used, it's swapped with `Repeated_Offset1`.


Dictionary format
-----------------

`zstd` is compatible with "raw content" dictionaries, free of any format restriction.
But dictionaries created by `zstd --train` follow a format, described here.

__Pre-requisites__ : a dictionary has a size,
                     defined either by a buffer limit, or a file size.

| `Magic_Number` | `Dictionary_ID` | `Entropy_Tables` | `Content` |
| -------------- | --------------- | ---------------- | --------- |

__`Magic_Number`__ : 4 bytes ID, value 0xEC30A437, little-endian format

__`Dictionary_ID`__ : 4 bytes, stored in little-endian format.
              `Dictionary_ID` can be any value, except 0 (which means no `Dictionary_ID`).
              It's used by decoders to check if they use the correct dictionary.

_Reserved ranges :_
              If the frame is going to be distributed in a private environment,
              any `Dictionary_ID` can be used.
              However, for public distribution of compressed frames,
              the following ranges are reserved for future use and should not be used :

              - low range : 1 - 32767
              - high range : >= (2^31)

__`Entropy_Tables`__ : following the same format as a [compressed blocks].
              They are stored in following order :
              Huffman tables for literals, FSE table for offsets,
              FSE table for match lengths, and FSE table for literals lengths.
              It's finally followed by 3 offset values, populating recent offsets,
              stored in order, 4-bytes little-endian each, for a total of 12 bytes.
              Each recent offset must have a value < dictionary size.

__`Content`__ : The rest of the dictionary is its content.
              The content act as a "past" in front of data to compress or decompress.

[compressed blocks]: #the-format-of-compressed_block

Appendix A - Decoding tables for predefined codes
-------------------------------------------------

This appendix contains FSE decoding tables for the predefined literal length, match length, and offset
codes. The tables have been constructed using the algorithm as given above in the
"from normalized distribution to decoding tables" chapter. The tables here can be used as examples
to crosscheck that an implementation implements the decoding table generation algorithm correctly.

#### Literal Length Code:

| State | Symbol | Number_Of_Bits | Base |
| ----- | ------ | -------------- | ---- |
|     0 |      0 |              4 |    0 |
|     1 |      0 |              4 |   16 |
|     2 |      1 |              5 |   32 |
|     3 |      3 |              5 |    0 |
|     4 |      4 |              5 |    0 |
|     5 |      6 |              5 |    0 |
|     6 |      7 |              5 |    0 |
|     7 |      9 |              5 |    0 |
|     8 |     10 |              5 |    0 |
|     9 |     12 |              5 |    0 |
|    10 |     14 |              6 |    0 |
|    11 |     16 |              5 |    0 |
|    12 |     18 |              5 |    0 |
|    13 |     19 |              5 |    0 |
|    14 |     21 |              5 |    0 |
|    15 |     22 |              5 |    0 |
|    16 |     24 |              5 |    0 |
|    17 |     25 |              5 |   32 |
|    18 |     26 |              5 |    0 |
|    19 |     27 |              6 |    0 |
|    20 |     29 |              6 |    0 |
|    21 |     31 |              6 |    0 |
|    22 |      0 |              4 |   32 |
|    23 |      1 |              4 |    0 |
|    24 |      2 |              5 |    0 |
|    25 |      4 |              5 |   32 |
|    26 |      5 |              5 |    0 |
|    27 |      7 |              5 |   32 |
|    28 |      8 |              5 |    0 |
|    29 |     10 |              5 |   32 |
|    30 |     11 |              5 |    0 |
|    31 |     13 |              6 |    0 |
|    32 |     16 |              5 |   32 |
|    33 |     17 |              5 |    0 |
|    34 |     19 |              5 |   32 |
|    35 |     20 |              5 |    0 |
|    36 |     22 |              5 |   32 |
|    37 |     23 |              5 |    0 |
|    38 |     25 |              4 |    0 |
|    39 |     25 |              4 |   16 |
|    40 |     26 |              5 |   32 |
|    41 |     28 |              6 |    0 |
|    42 |     30 |              6 |    0 |
|    43 |      0 |              4 |   48 |
|    44 |      1 |              4 |   16 |
|    45 |      2 |              5 |   32 |
|    46 |      3 |              5 |   32 |
|    47 |      5 |              5 |   32 |
|    48 |      6 |              5 |   32 |
|    49 |      8 |              5 |   32 |
|    50 |      9 |              5 |   32 |
|    51 |     11 |              5 |   32 |
|    52 |     12 |              5 |   32 |
|    53 |     15 |              6 |    0 |
|    54 |     17 |              5 |   32 |
|    55 |     18 |              5 |   32 |
|    56 |     20 |              5 |   32 |
|    57 |     21 |              5 |   32 |
|    58 |     23 |              5 |   32 |
|    59 |     24 |              5 |   32 |
|    60 |     35 |              6 |    0 |
|    61 |     34 |              6 |    0 |
|    62 |     33 |              6 |    0 |
|    63 |     32 |              6 |    0 |

#### Match Length Code:

| State | Symbol | Number_Of_Bits | Base |
| ----- | ------ | -------------- | ---- |
|     0 |      0 |              6 |    0 |
|     1 |      1 |              4 |    0 |
|     2 |      2 |              5 |   32 |
|     3 |      3 |              5 |    0 |
|     4 |      5 |              5 |    0 |
|     5 |      6 |              5 |    0 |
|     6 |      8 |              5 |    0 |
|     7 |     10 |              6 |    0 |
|     8 |     13 |              6 |    0 |
|     9 |     16 |              6 |    0 |
|    10 |     19 |              6 |    0 |
|    11 |     22 |              6 |    0 |
|    12 |     25 |              6 |    0 |
|    13 |     28 |              6 |    0 |
|    14 |     31 |              6 |    0 |
|    15 |     33 |              6 |    0 |
|    16 |     35 |              6 |    0 |
|    17 |     37 |              6 |    0 |
|    18 |     39 |              6 |    0 |
|    19 |     41 |              6 |    0 |
|    20 |     43 |              6 |    0 |
|    21 |     45 |              6 |    0 |
|    22 |      1 |              4 |   16 |
|    23 |      2 |              4 |    0 |
|    24 |      3 |              5 |   32 |
|    25 |      4 |              5 |    0 |
|    26 |      6 |              5 |   32 |
|    27 |      7 |              5 |    0 |
|    28 |      9 |              6 |    0 |
|    29 |     12 |              6 |    0 |
|    30 |     15 |              6 |    0 |
|    31 |     18 |              6 |    0 |
|    32 |     21 |              6 |    0 |
|    33 |     24 |              6 |    0 |
|    34 |     27 |              6 |    0 |
|    35 |     30 |              6 |    0 |
|    36 |     32 |              6 |    0 |
|    37 |     34 |              6 |    0 |
|    38 |     36 |              6 |    0 |
|    39 |     38 |              6 |    0 |
|    40 |     40 |              6 |    0 |
|    41 |     42 |              6 |    0 |
|    42 |     44 |              6 |    0 |
|    43 |      1 |              4 |   32 |
|    44 |      1 |              4 |   48 |
|    45 |      2 |              4 |   16 |
|    46 |      4 |              5 |   32 |
|    47 |      5 |              5 |   32 |
|    48 |      7 |              5 |   32 |
|    49 |      8 |              5 |   32 |
|    50 |     11 |              6 |    0 |
|    51 |     14 |              6 |    0 |
|    52 |     17 |              6 |    0 |
|    53 |     20 |              6 |    0 |
|    54 |     23 |              6 |    0 |
|    55 |     26 |              6 |    0 |
|    56 |     29 |              6 |    0 |
|    57 |     52 |              6 |    0 |
|    58 |     51 |              6 |    0 |
|    59 |     50 |              6 |    0 |
|    60 |     49 |              6 |    0 |
|    61 |     48 |              6 |    0 |
|    62 |     47 |              6 |    0 |
|    63 |     46 |              6 |    0 |

#### Offset Code:

| State | Symbol | Number_Of_Bits | Base |
| ----- | ------ | -------------- | ---- |
|     0 |      0 |              5 |    0 |
|     1 |      6 |              4 |    0 |
|     2 |      9 |              5 |    0 |
|     3 |     15 |              5 |    0 |
|     4 |     21 |              5 |    0 |
|     5 |      3 |              5 |    0 |
|     6 |      7 |              4 |    0 |
|     7 |     12 |              5 |    0 |
|     8 |     18 |              5 |    0 |
|     9 |     23 |              5 |    0 |
|    10 |      5 |              5 |    0 |
|    11 |      8 |              4 |    0 |
|    12 |     14 |              5 |    0 |
|    13 |     20 |              5 |    0 |
|    14 |      2 |              5 |    0 |
|    15 |      7 |              4 |   16 |
|    16 |     11 |              5 |    0 |
|    17 |     17 |              5 |    0 |
|    18 |     22 |              5 |    0 |
|    19 |      4 |              5 |    0 |
|    20 |      8 |              4 |   16 |
|    21 |     13 |              5 |    0 |
|    22 |     19 |              5 |    0 |
|    23 |      1 |              5 |    0 |
|    24 |      6 |              4 |   16 |
|    25 |     10 |              5 |    0 |
|    26 |     16 |              5 |    0 |
|    27 |     28 |              5 |    0 |
|    28 |     27 |              5 |    0 |
|    29 |     26 |              5 |    0 |
|    30 |     25 |              5 |    0 |
|    31 |     24 |              5 |    0 |

Version changes
---------------
- 0.2.2 : added predefined codes, by Johannes Rudolph
- 0.2.1 : clarify field names, by Przemyslaw Skibinski
- 0.2.0 : numerous format adjustments for zstd v0.8
- 0.1.2 : limit Huffman tree depth to 11 bits
- 0.1.1 : reserved dictID ranges
- 0.1.0 : initial release
# Contributing to Zstandard
We want to make contributing to this project as easy and transparent as
possible.

## Our Development Process
New versions are being developed in the "dev" branch,
or in their own feature branch.
When they are deemed ready for a release, they are merged into "master".

As a consequences, all contributions must stage first through "dev"
or their own feature branch.

## Pull Requests
We actively welcome your pull requests.

1. Fork the repo and create your branch from `dev`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. If you haven't already, complete the Contributor License Agreement ("CLA").

## Contributor License Agreement ("CLA")
In order to accept your pull request, we need you to submit a CLA. You only need
to do this once to work on any of Facebook's open source projects.

Complete your CLA here: <https://code.facebook.com/cla>

## Issues
We use GitHub issues to track public bugs. Please ensure your description is
clear and has sufficient instructions to be able to reproduce the issue.

Facebook has a [bounty program](https://www.facebook.com/whitehat/) for the safe
disclosure of security bugs. In those cases, please go through the process
outlined on that page and do not file a public issue.

## Coding Style  
* 4 spaces for indentation rather than tabs

## License
By contributing to Zstandard, you agree that your contributions will be licensed
under the [LICENSE](LICENSE) file in the root directory of this source tree.
Zstandard library files
================================

The __lib__ directory contains several directories.
Depending on target use case, it's enough to include only files from relevant directories.


#### API

Zstandard's stable API is exposed within [zstd.h](zstd.h),
at the root of `lib` directory.


#### Advanced API

Some additional API may be useful if you're looking into advanced features :
- common/error_public.h : transforms `size_t` function results into an `enum`,
                          for precise error handling.
- ZSTD_STATIC_LINKING_ONLY : if you define this macro _before_ including `zstd.h`,
                          it will give access to advanced and experimental API.
                          These APIs shall ___never be used with dynamic library___ !
                          They are not "stable", their definition may change in the future.
                          Only static linking is allowed.


#### Modular build

Directory `common/` is required in all circumstances.
You can select to support compression only, by just adding files from the `compress/` directory,
In a similar way, you can build a decompressor-only library with the `decompress/` directory.

Other optional functionalities provided are :

- `dictBuilder/`  : source files to create dictionaries.
                    The API can be consulted in `dictBuilder/zdict.h`.
                    This module also depends on `common/` and `compress/` .

- `legacy/` : source code to decompress previous versions of zstd, starting from `v0.1`.
              This module also depends on `common/` and `decompress/` .
              Library compilation must include directive `ZSTD_LEGACY_SUPPORT = 1` .
              The main API can be consulted in `legacy/zstd_legacy.h`.
              Advanced API from each version can be found in their relevant header file.
              For example, advanced API for version `v0.4` is in `legacy/zstd_v04.h` .


#### Obsolete streaming API

Streaming is now provided within `zstd.h`.
Older streaming API is still provided within `common/zbuff.h`.
It is considered obsolete, and will be removed in a future version.
Consider migrating towards newer streaming API.


#### Miscellaneous

The other files are not source code. There are :

 - LICENSE : contains the BSD license text
 - Makefile : script to compile or install zstd library (static and dynamic)
 - libzstd.pc.in : for pkg-config (`make install`)
 - README.md : this file
Programs and scripts for automated testing of Zstandard
=======================================================

This directory contains the following programs and scripts:
- `datagen` : Synthetic and parametrable data generator, for tests
- `fullbench`  : Precisely measure speed for each zstd inner functions
- `fuzzer`  : Test tool, to check zstd integrity on target platform
- `paramgrill` : parameter tester for zstd
- `test-zstd-speed.py` : script for testing zstd speed difference between commits
- `test-zstd-versions.py` : compatibility test between zstd versions stored on Github (v0.1+)
- `zbufftest`  : Test tool to check ZBUFF (a buffered streaming API) integrity
- `zstreamtest` : Fuzzer test tool for zstd streaming API


#### `test-zstd-versions.py` - script for testing zstd interoperability between versions

This script creates `versionsTest` directory to which zstd repository is cloned.
Then all taged (released) versions of zstd are compiled.
In the following step interoperability between zstd versions is checked.


#### `test-zstd-speed.py` - script for testing zstd speed difference between commits

This script creates `speedTest` directory to which zstd repository is cloned.
Then it compiles all branches of zstd and performs a speed benchmark for a given list of files (the `testFileNames` parameter).
After `sleepTime` (an optional parameter, default 300 seconds) seconds the script checks repository for new commits.
If a new commit is found it is compiled and a speed benchmark for this commit is performed.
The results of the speed benchmark are compared to the previous results.
If compression or decompression speed for one of zstd levels is lower than `lowerLimit` (an optional parameter, default 0.98) the speed benchmark is restarted.
If second results are also lower than `lowerLimit` the warning e-mail is send to recipients from the list (the `emails` parameter).

Additional remarks:
- To be sure that speed results are accurate the script should be run on a "stable" target system with no other jobs running in parallel
- Using the script with virtual machines can lead to large variations of speed results
- The speed benchmark is not performed until computers' load average is lower than `maxLoadAvg` (an optional parameter, default 0.75)
- The script sends e-mails using `mutt`; if `mutt` is not available it sends e-mails without attachments using `mail`; if both are not available it only prints a warning


The example usage with two test files, one e-mail address, and with an additional message:
```
./test-zstd-speed.py "silesia.tar calgary.tar" "email@gmail.com" --message "tested on my laptop" --sleepTime 60
``` 

To run the script in background please use:
```
nohup ./test-zstd-speed.py testFileNames emails &
```

The full list of parameters:
```
positional arguments:
  testFileNames         file names list for speed benchmark
  emails                list of e-mail addresses to send warnings

optional arguments:
  -h, --help            show this help message and exit
  --message MESSAGE     attach an additional message to e-mail
  --lowerLimit LOWERLIMIT
                        send email if speed is lower than given limit
  --maxLoadAvg MAXLOADAVG
                        maximum load average to start testing
  --lastCLevel LASTCLEVEL
                        last compression level for testing
  --sleepTime SLEEPTIME
                        frequency of repository checking in seconds
```
Zstandard wrapper for zlib
================================

The main objective of creating a zstd wrapper for [zlib](http://zlib.net/) is to allow a quick and smooth transition to zstd for projects already using zlib.

#### Required files

To build the zstd wrapper for zlib the following files are required:
- zlib.h
- a static or dynamic zlib library
- zlibWrapper/zstd_zlibwrapper.h
- zlibWrapper/zstd_zlibwrapper.c
- a static or dynamic zstd library

The first two files are required by all projects using zlib and they are not included with the zstd distribution.
The further files are supplied with the zstd distribution.


#### Embedding the zstd wrapper within your project

Let's assume that your project that uses zlib is compiled with:
```gcc project.o -lz```

To compile the zstd wrapper with your project you have to do the following:
- change all references with ```#include "zlib.h"``` to ```#include "zstd_zlibwrapper.h"```
- compile your project with `zstd_zlibwrapper.c` and a static or dynamic zstd library

The linking should be changed to:
```gcc project.o zstd_zlibwrapper.o -lz -lzstd```


#### Enabling zstd compression within your project

After embedding the zstd wrapper within your project the zstd library is turned off by default.
Your project should work as before with zlib. There are two options to enable zstd compression:
- compilation with ```-DZWRAP_USE_ZSTD=1``` (or using ```#define ZWRAP_USE_ZSTD 1``` before ```#include "zstd_zlibwrapper.h"```)
- using the ```void ZWRAP_useZSTDcompression(int turn_on)``` function (declared in ```#include "zstd_zlibwrapper.h"```)

During decompression zlib and zstd streams are automatically detected and decompressed using a proper library.
This behavior can be changed using `ZWRAP_setDecompressionType(ZWRAP_FORCE_ZLIB)` what will make zlib decompression slightly faster.


#### Example
We have take the file ```test/example.c``` from [the zlib library distribution](http://zlib.net/) and copied it to [zlibWrapper/examples/example.c](examples/example.c).
After compilation and execution it shows the following results: 
```
zlib version 1.2.8 = 0x1280, compile flags = 0x65
uncompress(): hello, hello!
gzread(): hello, hello!
gzgets() after gzseek:  hello!
inflate(): hello, hello!
large_inflate(): OK
after inflateSync(): hello, hello!
inflate with dictionary: hello, hello!
```
Then we have changed ```#include "zlib.h"``` to ```#include "zstd_zlibwrapper.h"```, compiled the [example.c](examples/example.c) file
with ```-DZWRAP_USE_ZSTD=1``` and linked with additional ```zstd_zlibwrapper.o -lzstd```.
We were forced to turn off the following functions: ```test_gzio```, ```test_flush```, ```test_sync``` which use currently unsupported features.
After running it shows the following results:
```
zlib version 1.2.8 = 0x1280, compile flags = 0x65
uncompress(): hello, hello!
inflate(): hello, hello!
large_inflate(): OK
inflate with dictionary: hello, hello!
```
The script used for compilation can be found at [zlibWrapper/Makefile](Makefile).


#### The measurement of performace of Zstandard wrapper for zlib

The zstd distribution contains a tool called `zwrapbench` which can measure speed and ratio of zlib, zstd, and the wrapper.
The benchmark is conducted using given filenames or synthetic data if filenames are not provided.
The files are read into memory and joined together. 
It makes benchmark more precise as it eliminates I/O overhead. 
Many filenames can be supplied as multiple parameters, parameters with wildcards or names of directories can be used as parameters with the -r option.
One can select compression levels starting from `-b` and ending with `-e`. The `-i` parameter selects minimal time used for each of tested levels.
With `-B` option bigger files can be divided into smaller, independently compressed blocks. 
The benchmark tool can be compiled with `make zwrapbench` using [zlibWrapper/Makefile](Makefile).


#### Improving speed of streaming compression

During streaming compression the compressor never knows how big is data to compress.
Zstandard compression can be improved by providing size of source data to the compressor. By default streaming compressor assumes that data is bigger than 256 KB but it can hurt compression speed on smaller data. 
The zstd wrapper provides the `ZWRAP_setPledgedSrcSize()` function that allows to change a pledged source size for a given compression stream.
The function will change zstd compression parameters what may improve compression speed and/or ratio.
It should be called just after `deflateInit()`or `deflateReset()` and before `deflate()` or `deflateSetDictionary()`. The function is only helpful when data is compressed in blocks. There will be no change in case of `deflateInit()` or `deflateReset()`  immediately followed by `deflate(strm, Z_FINISH)`
as this case is automatically detected.


#### Reusing contexts

The ordinary zlib compression of two files/streams allocates two contexts:
- for the 1st file calls `deflateInit`, `deflate`, `...`, `deflate`, `defalateEnd`
- for the 2nd file calls `deflateInit`, `deflate`, `...`, `deflate`, `defalateEnd`

The speed of compression can be improved with reusing a single context with following steps:
- initialize the context with `deflateInit`
- for the 1st file call `deflate`, `...`, `deflate`
- for the 2nd file call `deflateReset`, `deflate`, `...`, `deflate`
- free the context with `deflateEnd`

To check the difference we made experiments using `zwrapbench -ri6b6` with zstd and zlib compression (both at level 6).
The input data was decompressed git repository downloaded from https://github.com/git/git/archive/master.zip which contains 2979 files.
The table below shows that reusing contexts has a minor influence on zlib but it gives improvement for zstd.
In our example (the last 2 lines) it gives 4% better compression speed and 5% better decompression speed.

| Compression type                                  | Compression | Decompress.| Compr. size | Ratio |
| ------------------------------------------------- | ------------| -----------| ----------- | ----- |
| zlib 1.2.8                                        |  30.51 MB/s | 219.3 MB/s |     6819783 | 3.459 |
| zlib 1.2.8 not reusing a context                  |  30.22 MB/s | 218.1 MB/s |     6819783 | 3.459 |
| zlib 1.2.8 with zlibWrapper and reusing a context |  30.40 MB/s | 218.9 MB/s |     6819783 | 3.459 |
| zlib 1.2.8 with zlibWrapper not reusing a context |  30.28 MB/s | 218.1 MB/s |     6819783 | 3.459 |
| zstd 1.1.0 using ZSTD_CCtx                        |  68.35 MB/s | 430.9 MB/s |     6868521 | 3.435 |
| zstd 1.1.0 using ZSTD_CStream                     |  66.63 MB/s | 422.3 MB/s |     6868521 | 3.435 |
| zstd 1.1.0 with zlibWrapper and reusing a context |  54.01 MB/s | 403.2 MB/s |     6763482 | 3.488 |
| zstd 1.1.0 with zlibWrapper not reusing a context |  51.59 MB/s | 383.7 MB/s |     6763482 | 3.488 |


#### Compatibility issues
After enabling zstd compression not all native zlib functions are supported. When calling unsupported methods they put error message into strm->msg and return Z_STREAM_ERROR.

Supported methods:
- deflateInit
- deflate (with exception of Z_FULL_FLUSH, Z_BLOCK, and Z_TREES)
- deflateSetDictionary
- deflateEnd
- deflateReset
- deflateBound
- inflateInit
- inflate
- inflateSetDictionary
- inflateReset
- inflateReset2
- compress
- compress2
- compressBound
- uncompress

Ignored methods (they do nothing):
- deflateParams

Unsupported methods:
- gzip file access functions
- deflateCopy
- deflateTune
- deflatePending
- deflatePrime
- deflateSetHeader
- inflateGetDictionary
- inflateCopy
- inflateSync
- inflatePrime
- inflateMark
- inflateGetHeader
- inflateBackInit
- inflateBack
- inflateBackEnd
Zstandard library : usage examples
==================================

- [Simple compression](simple_compression.c) :
  Compress a single file.
  Introduces usage of : `ZSTD_compress()`

- [Simple decompression](simple_decompression.c) :
  Decompress a single file.
  Only compatible with simple compression.
  Result remains in memory.
  Introduces usage of : `ZSTD_decompress()`

- [Streaming compression](streaming_compression.c) :
  Compress a single file.
  Introduces usage of : `ZSTD_compressStream()`

- [Streaming decompression](streaming_decompression.c) :
  Decompress a single file compressed by zstd.
  Compatible with both simple and streaming compression.
  Result is sent to stdout.
  Introduces usage of : `ZSTD_decompressStream()`

- [Dictionary compression](dictionary_compression.c) :
  Compress multiple files using the same dictionary.
  Introduces usage of : `ZSTD_createCDict()` and `ZSTD_compress_usingCDict()`

- [Dictionary decompression](dictionary_decompression.c) :
  Decompress multiple files using the same dictionary.
  Result remains in memory.
  Introduces usage of : `ZSTD_createDDict()` and `ZSTD_decompress_usingDDict()`
Command Line Interface for Zstandard library
============================================

Command Line Interface (CLI) can be created using the `make` command without any additional parameters.
There are however other Makefile targets that create different variations of CLI:
- `zstd` : default CLI supporting gzip-like arguments; includes dictionary builder, benchmark, and support for decompression of legacy zstd versions
- `zstd32` : Same as `zstd`, but forced to compile in 32-bits mode
- `zstd_nolegacy` : Same as `zstd` except of support for decompression of legacy zstd versions
- `zstd-small` : CLI optimized for minimal size; without dictionary builder, benchmark, and support for decompression of legacy zstd versions
- `zstd-compress` : compressor-only version of CLI; without dictionary builder, benchmark, and support for decompression of legacy zstd versions
- `zstd-decompress` : decompressor-only version of CLI; without dictionary builder, benchmark, and support for decompression of legacy zstd versions


#### Aggregation of parameters
CLI supports aggregation of parameters i.e. `-b1`, `-e18`, and `-i1` can be joined into `-b1e18i1`. 


#### Dictionary builder in Command Line Interface
Zstd offers a training mode, which can be used to tune the algorithm for a selected
type of data, by providing it with a few samples. The result of the training is stored
in a file selected with the `-o` option (default name is `dictionary`),
which can be loaded before compression and decompression.

Using a dictionary, the compression ratio achievable on small data improves dramatically.
These compression gains are achieved while simultaneously providing faster compression and decompression speeds.
Dictionary work if there is some correlation in a family of small data (there is no universal dictionary). 
Hence, deploying one dictionary per type of data will provide the greater benefits.
Dictionary gains are mostly effective in the first few KB. Then, the compression algorithm
will rely more and more on previously decoded content to compress the rest of the file.

Usage of the dictionary builder and created dictionaries with CLI:

1. Create the dictionary : `zstd --train FullPathToTrainingSet/* -o dictionaryName`
2. Compress with the dictionary: `zstd FILE -D dictionaryName`
3. Decompress with the dictionary: `zstd --decompress FILE.zst -D dictionaryName`



#### Benchmark in Command Line Interface
CLI includes in-memory compression benchmark module for zstd.
The benchmark is conducted using given filenames. The files are read into memory and joined together.
It makes benchmark more precise as it eliminates I/O overhead.
Many filenames can be supplied as multiple parameters, parameters with wildcards or
names of directories can be used as parameters with the `-r` option.

The benchmark measures ratio, compressed size, compression and decompression speed.
One can select compression levels starting from `-b` and ending with `-e`.
The `-i` parameter selects minimal time used for each of tested levels.



#### Usage of Command Line Interface
The full list of options can be obtained with `-h` or `-H` parameter:
```
Usage :
      zstd [args] [FILE(s)] [-o file]

FILE    : a filename
          with no FILE, or when FILE is - , read standard input
Arguments :
 -#     : # compression level (1-19, default:3)
 -d     : decompression
 -D file: use `file` as Dictionary
 -o file: result stored into `file` (only if 1 input file)
 -f     : overwrite output without prompting
--rm    : remove source file(s) after successful de/compression
 -k     : preserve source file(s) (default)
 -h/-H  : display help/long help and exit

Advanced arguments :
 -V     : display Version number and exit
 -v     : verbose mode; specify multiple times to increase log level (default:2)
 -q     : suppress warnings; specify twice to suppress errors too
 -c     : force write to standard output, even if it is the console
 -r     : operate recursively on directories
--ultra : enable levels beyond 19, up to 22 (requires more memory)
--no-dictID : don't write dictID into header (dictionary compression)
--[no-]check : integrity check (default:enabled)
--test  : test compressed file integrity
--[no-]sparse : sparse mode (default:enabled on file, disabled on stdout)

Dictionary builder :
--train ## : create a dictionary from a training set of files
 -o file : `file` is dictionary name (default: dictionary)
--maxdict ## : limit dictionary to specified size (default : 112640)
 -s#    : dictionary selectivity level (default: 9)
--dictID ## : force dictionary ID to specified value (default: random)

Benchmark arguments :
 -b#    : benchmark file(s), using # compression level (default : 1)
 -e#    : test all compression levels from -bX to # (default: 1)
 -i#    : minimum evaluation time in seconds (default : 3s)
 -B#    : cut file into independent blocks of size # (default: no block)
 ```# Parallel Zstandard (PZstandard)

Parallel Zstandard is a Pigz-like tool for Zstandard.
It provides Zstandard format compatible compression and decompression that is able to utilize multiple cores.
It breaks the input up into equal sized chunks and compresses each chunk independently into a Zstandard frame.
It then concatenates the frames together to produce the final compressed output.
Pzstandard will write a 12 byte header for each frame that is a skippable frame in the Zstandard format, which tells PZstandard the size of the next compressed frame.
PZstandard supports parallel decompression of files compressed with PZstandard.
When decompressing files compressed with Zstandard, PZstandard does IO in one thread, and decompression in another.

## Usage

PZstandard supports the same command line interface as Zstandard, but also provies the `-p` option to specify the number of threads.
Dictionary mode is not currently supported.

Basic usage

    pzstd input-file -o output-file -p num-threads -#          # Compression
    pzstd -d input-file -o output-file -p num-threads          # Decompression

PZstandard also supports piping and fifo pipes

    cat input-file | pzstd -p num-threads -# -c > /dev/null

For more options

    pzstd --help

PZstandard tries to pick a smart default number of threads if not specified (displayed in `pzstd --help`).
If this number is not suitable, during compilation you can define `PZSTD_NUM_THREADS` to the number of threads you prefer.

## Benchmarks

As a reference, PZstandard and Pigz were compared on an Intel Core i7 @ 3.1 GHz, each using 4 threads, with the [Silesia compression corpus](http://sun.aei.polsl.pl/~sdeor/index.php?page=silesia).

Compression Speed vs Ratio with 4 Threads | Decompression Speed with 4 Threads
------------------------------------------|-----------------------------------
![Compression Speed vs Ratio](images/Cspeed.png "Compression Speed vs Ratio") | ![Decompression Speed](images/Dspeed.png "Decompression Speed")

The test procedure was to run each of the following commands 2 times for each compression level, and take the minimum time.

    time pzstd -# -p 4    -c silesia.tar     > silesia.tar.zst
    time pzstd -d -p 4    -c silesia.tar.zst > /dev/null

    time pigz  -# -p 4 -k -c silesia.tar     > silesia.tar.gz
    time pigz  -d -p 4 -k -c silesia.tar.gz  > /dev/null

PZstandard was tested using compression levels 1-19, and Pigz was tested using compression levels 1-9.
Pigz cannot do parallel decompression, it simply does each of reading, decompression, and writing on separate threads.

## Tests

Tests require that you have [gtest](https://github.com/google/googletest) installed.
Modify `GTEST_INC` and `GTEST_LIB` in `test/Makefile` and `utils/test/Makefile` to work for your install of gtest.
Then run `make test` in the `contrib/pzstd` directory.
