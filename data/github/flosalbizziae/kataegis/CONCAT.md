---
title: "kataegis"
author: Xue Lin, Jian Li
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{get-started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

+ Contents
    + Introduction
    + Overview
    + Data Input
    + Read In Data
    + Calculate The Inter-mutational Distances
    + Kataegis Region Identification
    + Kataegis Region Visualization
    + References

## Introduction

This is a documentation for kataegis package which provides tools for the 
identification and visualization of the genomic localized hypermutation region, 
i.e. kataegis. Researchers have found that the localized hypermutations occur 
in a way that is quite different from the usual somatic muations in the aspect of
the mutations' occurency and genomic density. This hypermuatation is like a violent 
strom, that is what the greek word '<em>kataegis</em>' means. This phenomena was found in 
the breast and squamous cell carcinoma (SCC), more and more reseaches report this 
phenomena in other types of cancers as well^1^. There were evidences showed that 
the kataegis, which is induced by AID/APOBEC cytosine deaminase, is related to 
the proganosis, and the possible treatment strategies decision^2,3^. Thus, the 
package can meet the needs of a light-weighted and simple-to-use toolkit for quick 
identifying and visualizing the kataegis of the genomes.

## Overview

<img src="overview.jpeg" width="680" height="340">

## Data Input

There are two accepted formats of variant file by this package currently: VCF and MAF. 

The VCF format is an abbreviation for the Variant Call Format，which is a text 
file format for storing variants information. For this format of file, the columns 
of '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', and 'FILTER' are obliged. Detailed 
information of the VCF format please refer to the Variant Call Format (VCF) 
Version 4.2 Specification^4^.

The MAF format is an abbreviation for Mutation Annotation Format, which is also a 
text fie format for storing variants information for one or more than one samples.
For this format of files, the columns of 'Chromosome', 'Start_Position',
'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', are obliged. Detailed
information of the MAF format, please refer to the Mutation Annotation Format
(MAF) GDC Version 1.0.0 Specification^5^.

It is strongly recomendated to filter and convert your own data before you read 
them in, though the funcions has their own crude filters.

We have two psudo example data files come with this pacakge: test_mut.vcf and test_mut.maf,
which will be used to show the usage of the analysis procedure.

## Read In Data

Here is how we read in a VCF file:

```{r echo=TRUE}
library(kataegis)
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)
head(dat)
```

Note: As the VCF file is a psudo one, the 'FILTER' field does not contain detailed 
information, that's why a warning message rise. The data will still be read
in without any obstruction for further analysis, but it is strongly recommendated
for users to filter and add the filter result, e.g. 'PASS' for passing all filters,
before you read in data to avoid any bias for the further steps.

When you need to read in a MAF file, you can choose two methods: to merge all the 
samples together, e.g. replicates, or to split each samples.

```{r echo=TRUE}
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath) # merge all the samples together
head(dat)
```

The output of this method generate a matrix just like readVCF did, while the 
other split mode method will generate a list of matrix.

```{r echo=TRUE}
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath, split=TRUE) # split each sample in the MAF file
str(dat)
```

And the 'ESCC-148T' and 'ESCC-184T' are the sample IDs.

## Calculate The Inter-mutational Distances

As the kataegis is usually identified as a small genomic region with high 
occurency and density of mutations, we used a method call piece-wise constant 
fitting^6^ to find the regions of constant inter-mutational distance, thus a inter
-mutational distances of each variant should be calculated before kataegis iden
-tification.

```{r echo=TRUE}
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)
dist<-interMutDist(dat)
head(dist)
```

The output data frame contains the columns of 'CHR', 'POS', and 'DISTANCE', which
presents the chromosome of the variant, the position of the variant and the 
distance of the variant to its neighbour variant respectively.

The mergeed MAF read in file is exactly the same analysis step with the VCF read
in file. If you're working with splitted read in MAF, then

```{r echo=TRUE}
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath, split=TRUE) #merge all the samples together
dist<-interMutDist(dat)
str(dist)
```

As there are samples with too few mutations wich is not sufficient to calculate 
the inter-muational distances, these samples will be abandoned and a warning 
describing the situation and listing the samples' ID will arise.

## Kataegis Region Identification

```{r echo=FALSE}
seg<-kata(dist) # to identify the kataegis region with default paramters
```

As the psudo data has two few mutations to call kataegis, wich is usually identified
as a genomic region with six or more mutations and average inter-muational distances
less or equal to 1000bp, the kataegis identification just failed on the psudo data
with default parameters. To present the full analysis pipeline, we will use the 
use-defined parameters in this documentation.

```{r echo=TRUE}
#For MAF files
fpath <- system.file("extdata", "test_mut.maf", package="kataegis") #get the path of the MAF file
dat<-readMAF(fpath, split=TRUE)# read in the maf file with samples seperated
dist<-interMutDist(dat) # calcualte the inter-mutational distances
seg<-kata(dist, len=5e6, assembly="hg19") # to identify the kataegis region with user set paramters
str(seg)

#For VCF files
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis") # get the path of the VCF file
dat<-readVCF(fpath) # read in the vcf file
dist<-interMutDist(dat) # calculate the inter-muational distances
seg<-kata(dist, len=5e4, assembly="mm10") # to identify the kataegis region with user set paramters
head(seg)
```

The analysis for MAF file with samples merged or seperated uses the same procedure
as the VCF file analysis, except that the 'split' is set to 'TRUE' when read in 
the MAF file. And in the functions for kataegis identificaiton and visualizaiton,
the seperated samples will be analyzed individually and the output will be named
by the sample IDs.

Please be notified that the default parameters for identifying a kataegis region
is more than 6 variants within 1Kb genomic region, which is putative and it is
usually advisable to use the default parameters, but you can also set the
parameters according to your specific genomes.

## Kataegis Region visualization

If you have the kataegis identified, you can also visualize the nucleotides types,
mutation types and the distribution of the mutations on genome.

```{r echo=TRUE}
kataplot(dat, seg, dist)
```

If you have analyzed a VCF file or a merged MAF file, you'll get the visualization
of all the variants, or if you're dealing with a seperated MAF file, all the 
following plots will be generated for each sample.

<img src="all_content_bar.png" width="340" height="340">

This is the barplot of all the conversions of all the mutations in kataegis regions
identified.

<img src="all_spectra.png" width="680" height="340">
    
This is the stack barplot of all the nucleotides content of all the mutations in
kataegis regions identifed.

<img src="all_rainfall.png" width="680" height="340">
    
This is the rainfall plot of all the mutations distribution on the genome. As the 
limitation of the size of the example data, here we only have the VCF file of one
chromosome. If you have all the variants of the genome you'll also get the rain
-fall plot as the following. Thus, if you want to visualize the data seperately 
as a chromosome, you can split the VCF according to the chromosome before you 
start the analysis pipeline.

<img src="example_all_rainfall.png" width="680" height="340">
    
## References

1. B.Alexandrov, L., Nik-Zainal, S., Wedge, D. C., Aparicio, S. A. J. R., Behjati, S., Biankin, A. V., … et.al. (2013). Signatures of mutational processes in human cancer. Nature, 500, 415–421. https://doi.org/10.1038/nature12477

2. Lada, A. G., Dhar, A., Boissy, R. J., Hirano, M., Rubel, A. A., Rogozin, I. B., & Pavlov, Y. I. (2012). AID / APOBEC cytosine deaminase induces genome-wide kataegis, 1–7.

3. Han, M., Shin, S., Park, H., Kim, M. S., Lee, S. H., Jung, S. H., … Chung, Y. (2018). Mutational signatures and chromosome alteration pro fi les of squamous cell carcinomas of the vulva. Nature Publishing Group, 50(2), e442-13. https://doi.org/10.1038/emm.2017.265

4. Info, Q. F. (2020). The Variant Call Format ( VCF ) Version 4 . 2 Specification The VCF specification, 1–28.

5. https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

6. Nilsen, G., Liestøl, K., Loo, P. Van, Kristian, H., Vollan, M., Eide, M. B., … Lingjærde, O. C. (2012). Copynumber : Efficient algorithms for single- and multi-track copy number segmentation.


\name{readVCF}
\alias{readVCF}
\title{
Read In VCF File
}
\description{
This function can be used to read in the variation file in VCF format.
}
\usage{
readVCF(file)
}
\arguments{
  \item{file}{
a text file of VCF format.
}
}
\details{
The argument file, of which the columns of '#CHROM', 'POS', 'ID', 'REF', 'ALT',
'QUAL', and 'FILTER' are obliged. Detailed information of the VCF format please
refer to the Variant Call Format (VCF) Version 4.2 Specification.

The function will filter the variants according to the 'FILTER' column. Only all
the variants pass the filters that were set. If the 'FILTER' field has string
other than 'PASS', the function will rise a warning, and ask the user to decide
whether to stop and filter before use or continue with the analysis.

As the users usually used this function in PC, considering the memory capacity
of PCs and the usual size of the VCF files, the function only accept the files
of no larger than 2GB.
}
\value{
The function return a matrix of the variants passed all the filters with the
following columns:
  \item{#CHROM }{chromosome: An identifier from the reference genome or an
  angle-bracketed ID String (“<ID>”) pointing to a contig in the assembly file
  (cf. the ##assembly line in the header). All entries for a specific CHROM
  should form a contiguous block within the VCF file. (String, no white-space
  permitted, Required).}
  \item{POS }{position: The reference position, with the 1st base having
  position 1. Positions are sorted numerically, in increasing order, within
  each reference sequence CHROM. It is permitted to have multiple records with
  the same POS. Telomeres are indicated by using positions 0 or N+1, where N is
  the length of the corresponding chromosome or contig. (Integer, Required)}
  \item{REF }{reference base(s): Each base must be one of A,C,G,T,N (case
  insensitive). Multiple bases are permitted. The value in the POS field refers
  to the position of the first base in the String. For simple insertions and
  deletions in which either the REF or one of the ALT alleles would otherwise
  be null/empty, the REF and ALT Strings must include the base before the event
  (which must be reflected in the POS field), unless the event occurs at
  position 1 on the contig in which case it must include the base after the
  event; this padding base is not required (although it is permitted) for
  e.g. complex substitutions or other events where all alleles have at least one
  base represented in their Strings. If any of the ALT alleles is a symbolic
  allele (an angle-bracketed ID String “<ID>”) then the padding base is required
  and POS denotes the coordinate of the base preceding the polymorphism. Tools
  processing VCF files are not required to preserve case in the allele Strings.
  (String, Required).}
  \item{ALT }{alternate base(s): Comma separated list of alternate non-reference
  alleles. These alleles do not have to be called in any of the samples. Options
  are base Strings made up of the bases A,C,G,T,N,*, (case insensitive) or an
  angle-bracketed ID String (“<ID>”) or a breakend replacement string as
  described in the section on breakends. The ‘*’ allele is reserved to indicate
  that the allele is missing due to a upstream deletion. If there are no
  alternative alleles, then the missing value should be used. Tools processing
  VCF files are not required to preserve case in the allele String, except for
  IDs, which are case sensitive. (String; no whitespace, commas, or angle-brackets
  are permitted in the ID String itself)}
  \item{FILTER }{filter status: PASS if this position has passed all filters,
  i.e., a call is made at this position. Otherwise, if the site has not passed
  all filters, a semicolon-separated list of codes for filters that fail.
  e.g. “q10;s50” might indicate that at this site the quality is below 10 and
  the number of samples with data is below 50 percent  of the total number of
  samples.  ‘0’ is reserved and should not be used as a filter String. If filters
  have not been applied, then this field should be set to the missing value.
  (String, no white-space or semi-colons permitted) If the filter information
  is not in detail here, the program will take in all the variants and print out
  a warning message.}
}
\references{
Info, Q. F. (2020). The Variant Call Format ( VCF ) Version 4 . 2 Specification The VCF specification, 1–28.
}
\author{
Xue Lin, Jian Li
}
\note{
It is usually adviable to filter according to the filters,and convert your file
into the VCF 4.0 format before you read in your variants for further analysis.
And the .vcf file suffix is strongly recommendated.
}
\seealso{
readMAF
}
\examples{
#Read in VCF file:
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)

#Or you can also import vcf file typing
#dat<-readVCF("path/to/your/file.vcf")
}
\keyword{ IO }% use one of  RShowDoc("KEYWORDS")
\name{readMAF}
\alias{readMAF}
\title{
Read In MAF File
}
\description{
This function can be used to read in the variation file in MAF format in a merge mode or a seperated mode.
}
\usage{
readMAF(file, split = FALSE)
}
\arguments{
  \item{file}{
The argument file, of which the columns of 'Chromosome', 'Start_Position',
'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', are obliged. Detailed
information of the MAF format, please refer to the Mutation Annotation Format
(MAF) GDC Version 1.0.0 Specification.

As the users usually used this function in PC, considering the memory capacity
of PCs and the usual size of the VCF files, the function only accept the files
of no larger than 2GB.
}
  \item{split}{
The boolean value to define whether the MAF file is splitted according to the
samples ID. If the argument is TRUE, the MAF will be splitted into a list of
several matrixes of the variants. Default is FALSE. See 'Details'.
}
}
\details{
The argument split is TRUE, the MAF file will be splitted to each sample. In
other word, after the readMAF function with the splitted argument set TRUE, each
sample will has its own matrix returned, with the columns exactly like the
readVCF does, and all the matrixes will be gathered in a list, in which each
matrix can be accessed by its corresponding sample ID.
}
\value{
If the split argument is not set, a matrix will be returned, and the structures
is exactly like the matrix returned by readVCF. Otherwise a list of matrix will
be returned, and the structure will be as followings:
List
|_sample_ID_1
| |_matrix_1
|_sample_ID_2
| |_matrix_2
|_sample_ID_3
  |_matrix_3
}
\references{
https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
}
\author{
Xue Lin, Jian Li
}
\note{
It is usually adviable to filter according to the filters,and convert your file
into the GDC MAF v1.0.0 format before you read in your variants for further analysis.
And the .maf file suffix is strongly recommendated.
}
\seealso{
readMAF
}
\examples{
#Read in MAF file and merge all the samples in one file:
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath)

#Read in MAF file and split the maf according the samples' ID:
list_dat<-readMAF(fpath, split=TRUE)

#Or you can also import maf file by typing
#dat<-readMAF("path/to/your/file.maf")
#list_dat<-readMAF("path/to/your/file.maf", split=TRUE)
}
\keyword{ IO }
\name{interMutDist}
\alias{interMutDist}
\title{
Calculate The Inter-mutational Distances
}
\description{
The function will calculate the inter-mutational distances after the ordering of
the variants according to the genomic positions.
}
\usage{
interMutDist(data)
}
\arguments{
  \item{data}{
the read-in matrix or list of matrix generated by the readVCF or readMAF. See \
'Details'.
}
}
\details{
All the variants in a single matrix or a list of matrixes will be sorted
according to the genomics postition. Then the distance between two variants will
be calculated.
}
\value{
A dataframe will be returned with the following columns:
  \item{CHR }{The chromosome of the variant}
  \item{POS }{The position of the variant}
  \item{DISTANCE }{The distance of the variant to its neighbour variant}
}
\references{
To be published
}
\author{
Xue Lin, Jian Li
}
\note{
If there is too little variants for each chromsome to calculate the
inter-mutational distances, a warning will rise, and the program will stop.
}
\seealso{
readVCF, readMAF
}
\examples{
#library(kataegis)

##Read in MAF file and merge all the samples in one file:
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath)

#Read in MAF file and split the maf according the samples' ID:
list_dat<-readMAF(fpath, split=TRUE)

#Or you can also import maf file typing
#dat<-readMAF("path/to/your/file.maf")
#lsit_dat<-readMAF("path/to/your/file.maf", split=TRUE)

##Read in VCF file:
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)

#Or you can also import vcf file typing
#dat<-readVCF("path/to/your/file.vcf")

##calculate the inter-mutational distances
dist<-interMutDist(dat)
#list_dist<-interMutDist(list_dat)
}
\keyword{ math }% use one of  RShowDoc("KEYWORDS")
\name{kata}
\alias{kata}
\title{
Call The Kataegis Regions Of Your Data
}
\description{
The function will segment the inter-mutational distances using the pieceswise
constant fitting (PCF), then the kataegis will be called according to the
freuency of the variants' occurrance in a certain length of genomics region.
}
\usage{
kata(data, kmin = 2, gamma = 25, assembly = "hg19", len = 1000, nmut = 6, verbose = TRUE)
}
\arguments{
  \item{data}{
a dataframe or a list of dataframe of the inter-mutational distances calculated
and output by the interMutDist.
}
  \item{kmin}{
minimal number of variants in each segment, default is 2.
}
  \item{gamma}{
penalty for each discontinuity in the curve, default is 25.
}
  \item{assembly}{
the reference assembly version for analysis, default is human hg19, supported
assembly incl. hg38, hg16-19, and mm7-10.
}
  \item{len}{
the average inter-mutational distance to define the kataegis region.
}
  \item{nmut}{
the number of the variants contained in segments to define the kataegis region.
}
  \item{verbose}{
a verbose mode while segmenting the variants and calling kataegis.
}
}
\details{
If default parameters are used, the kataegis region (hypermutational region) will
be defined as a putative region with more than 6 variants withn a average inter-
muational distance of 1000bp.

The kmin and gamma are parameters used for PCF, and the default 2 and 25 were
trained on the set of kataegis foci that had been manually identified, curated
and validated using orthogonal sequencing platforms. See 'References'.
}
\value{
A dataframe or a list of dataframe of the kataegis region defined with the
following columns:
  \item{chr }{chromosome of the kataegis region locate}
  \item{arm }{the arm of the chromosome in which the kataegis region locate}
  \item{start.pos }{the start position of the kataegis of the chromosome}
  \item{end.pos }{the end position of the kataegis of the chromosome}
  \item{nmut }{the number of variants of the segment}
  \item{dsit.mean }{the average inter-mutational distance of the segment}
}
\references{
Lundmil B. Alexandrov, Serena Nik-Zaina, David C. Wedge, et.al(2013). Signatures
of mutational processes in human cancer. Nature 500, 415–421 (2013);
doi:10.1038/nature12477.

Gro Nilsen, Knut Liestol and Ole Christian Lingjaerde (2013). copynumber:
Segmentation of single- and multi-track copy number data by penalized least
squares regression. R package version 1.24.0.

Nilsen, G. & Liestol, K. et al. (2012) Copynumber: Efficient algorithms for
single- and multi-track copy number segmentation. BMC Genomics 13(1):591.
}
\author{
Xue Lin, Jian Li
}
\note{
Please be notified that the default parameters for identifying a kataegis region
is more than 6 variants within 1Kb genomic region, which is putative and it is
usually advisable to use the default parameters, but you can also set the
parameters according to your specific genomes.
}
\seealso{
interMutDistance
}
\examples{
#library(kataegis)

##Read in MAF file and merge all the samples in one file:
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath)

#Read in MAF file and split the maf according the samples' ID:
list_dat<-readMAF(fpath, split=TRUE)

#Or you can also import maf file typing
#dat<-readMAF("path/to/your/file.maf")
#lsit_dat<-readMAF("path/to/your/file.maf", split=TRUE)

##Read in VCF file:
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)

#Or you can also import vcf file typing
#dat<-readVCF("path/to/your/file.vcf")

##calculate the inter-mutational distances
dist<-interMutDist(dat)
#list_dist<-interMutDist(list_dat)

##call kataegis region with default parameters
#seg<-kata(dist)
#seg<-kata(list_dist)
}
\keyword{ math }% use one of  RShowDoc("KEYWORDS")
\name{kataplot}
\alias{kataplot}
\title{
Kataegis Region Visualization
}
\description{
This function will use the read-in data, the inter-mutational distance data, and
the kataegis segments to visualize the nucleotide spectra of the kataegis foci,
the necleotide conversion type, and the relation of the distribution of the
variants on the genome and the inter-mutational distances.
}
\usage{
kataplot(x, y, z, n = 20, fbar = TRUE, cbar = TRUE, rain = TRUE, colr = 6, fbw = 800,
fbh = 400, rw = 1000, rh = 500, assembly = "hg19", name = "all", type = "png")
}
\arguments{
  \item{x}{
the matrixm or a list of matrixes of the read in data from readVCF or readMAF
}
  \item{y}{
the dataframe or a list of dataframes of the kataegis data
}
  \item{z}{
the dataframe or a list of dataframes output from the interMutDist() function
}
  \item{n}{
the number of the division of the flanking region totally, default is 20, which
means that the upstream is devided into 10 and the downstream is devided into 10
}
  \item{fbar}{
the flanking region nucleotide spectra barplot
}
  \item{cbar}{
the core region mutation types barplot
}
  \item{rain}{
the rain-fall plot of all the raw data, and the poisson will also be tested
}
  \item{colr}{
the color set of the barplot and rainfall plot, default is the color palettes
from grDevices, rainbow(6),which is "#FF0000" "#FFFF00" "#00FF00" "#00FFFF"
"#0000FF" "#FF00FF". You can also use any R colors here. Theis value should be
equal or larger than 6.
}
  \item{fbw}{
the width of the flanking region nucleotides types barplot
}
  \item{fbh}{
the height of the flanking region nucleotides types barplot
}
  \item{rw}{
the width of the rainfall plot
}
  \item{rh}{
the height of the rainfall plot
}
  \item{assembly}{
the reference assembly version for analysis, default is human hg19, supported
assembly incl. hg38, hg16-19, and mm7-10.
}
  \item{name}{
the sample name for the plot output, default is all.
}
  \item{type}{
the plot output format, default is png, other are jpeg, and tiff
}
}
\value{
The details of the read in data please refere the manual page of readVCF(), readMAF(),
interMutDist() and kata().
}
\references{
Lundmil B. Alexandrov, Serena Nik-Zaina, David C. Wedge, et.al(2013). Signatures
of mutational processes in human cancer. Nature 500, 415–421 (2013);
doi:10.1038/nature12477.
}
\author{
Xue Lin, Jian Li
}
\note{
The output format and width and height of graphs could be set according to the
graph type, and the size of your data. The rainfall plot will use more memory,
in this case you could consider for good filtering at the very beginning of the
analysis.
}
\seealso{
kata, readVCF, readMAF
}
\examples{
#library(kataegis)

##Read in MAF file and merge all the samples in one file:
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath)

#Read in MAF file and split the maf according the samples' ID:
list_dat<-readMAF(fpath, split=TRUE)

#Or you can also import maf file typing
#dat<-readMAF("path/to/your/file.maf")
#list_dat<-readMAF("path/to/your/file.maf", split=TRUE)

##Read in VCF file:
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)

#Or you can also import vcf file typing
#dat<-readVCF("path/to/your/file.vcf")

##Calculate the inter-mutational distances
dist<-interMutDist(dat)
#list_dist<-interMutDist(list_dat)

#Call kataegis region with default parameters
#seg<-kata(dist)
#or
#list_seg<-kata(list_dist)

##Plot the result of the analysis
#kataplot(dat, seg, dist)
#or
#kataplot(lsit_dat, list_seg, list_dist)
}
\keyword{ graph }% use one of  RShowDoc("KEYWORDS")
