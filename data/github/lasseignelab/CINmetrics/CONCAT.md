## Test environments
* local OS X install(x86_64, darwin17.0), R 4.0.5
* win-builder (release, old release, development)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## Downstream dependencies
There are currently no downstream dependencies for this package

## Resubmission
This is a resubmission. In this version I have addressed comments by Gregor Seyer:

* Please do not start the description with "This package", package name,
title or similar. \
**Done**


* Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'CINmetrics' \
**Done**

* If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <[https:...]https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title") \
**Done**

* Please add small executable examples in your Rd-files to illustrate the
use of the exported function but also enable automatic testing. \
**Done**
# CINmetrics
CINmetrics R package

The goal of CINmetrics package is to provide different methods of calculating Chromosomal Instability (CIN) metrics from the literature that can be applied to any cancer data set including The Cancer Genome Atlas.
---
title: "CINmetrics"
author: "Vishal H. Oza"
date: "2-Dec-2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{CINmetrics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# CINmetrics 

The goal of CINmetrics package is to provide different methods of calculating Chromosomal Instability (CIN) metrics from the literature that can be applied to any cancer data set including The Cancer Genome Atlas.


```{r setup}
library(CINmetrics)
```
The dataset provided with CINmetrics package is masked Copy Number variation data for Breast Cancer for 10 unique samples selected randomly from TCGA.

```{r}
dim(maskCNV_BRCA)
```


Alternatively, you can download the entire dataset from TCGA using TCGAbiolinks package
```{r}
## Not run:
#library(TCGAbiolinks)
#query.maskCNV.hg39.BRCA <- GDCquery(project = "TCGA-BRCA",
#              data.category = "Copy Number Variation",
#              data.type = "Masked Copy Number Segment", legacy=FALSE)
#GDCdownload(query = query.maskCNV.hg39.BRCA)
#maskCNV.BRCA <- GDCprepare(query = query.maskCNV.hg39.BRCA, summarizedExperiment = FALSE)
#maskCNV.BRCA <- data.frame(maskCNV.BRCA, stringsAsFactors = FALSE)
#tai.test <- tai(cnvData = maskCNV.BRCA)
## End(Not run)
```

## Total Aberration Index

*tai* calculates the Total Aberration Index (TAI; [Baumbusch LO, et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3553118/)), "a measure of the abundance of genomic size of copy number changes in a tumour". It is defined as a weighted sum of the segment means ($|\bar{y}_{S_i}|$). 

Biologically, it can also be interpreted as the absolute deviation from the normal copy number state averaged over all genomic locations.

$$
Total\ Aberration\ Index = 
\frac
{\sum^{R}_{i = 1} {d_i} \cdot |{\bar{y}_{S_i}}|}
{\sum^{R}_{i = 1} {d_i}}\ \
where |\bar{y}_{S_i}| \ge |\log_2 1.7|
$$
```{r}
tai.test <- tai(cnvData = maskCNV_BRCA)
head(tai.test)
```

## Modified Total Aberration Index

*taiModified* calculates a modified Total Aberration Index using all sample values instead of those in aberrant copy number state, thus does not remove the directionality from the score.

$$
Modified\ Total\ Aberration\ Index = 
\frac
{\sum^{R}_{i = 1} {d_i} \cdot {\bar{y}_{S_i}}}
{\sum^{R}_{i = 1} {d_i}}
$$
```{r}
modified.tai.test <- taiModified(cnvData = maskCNV_BRCA)
head(modified.tai.test)
```

## Copy Number Aberration

*cna* calculates the total number of copy number aberrations (CNA; [Davidson JM, et. al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079079)), defined as a segment with copy number outside the pre-defined range of 1.7-2.3 ($(\log_2 1.7 -1) \le \bar{y}_{S_i} \le (\log_2 2.3 -1)$) that is not contiguous with an adjacent independent CNA of identical copy number. For our purposes, we have adapted the range to be $|\bar{y}_{S_i}| \ge |\log_2 1.7|$, which is only slightly larger than the original. 

This metric is very similar to the number of break points, but it comes with the caveat that adjacent segments need to have a difference in segmentation mean values.

$$
Total\ Copy\ Number\ Aberration = \sum^{R}_{i = 1} n_i
\ \ where\ \ 
\begin{align}
|\bar{y}_{S_i}| \ge |\log_2{1.7}|, \\ |\bar{y}_{S_{i-1}} - \bar{y}_{S_i}| \ge 0.2, \\ d_i \ge 10
\end{align}
$$
```{r}
cna.test <- cna(cnvData = maskCNV_BRCA)
head(cna.test)
```


## Counting Altered Base segments

*countingBaseSegments* calculates the number of altered bases defined as the sums of the lengths of segments ($d_i$) with an absolute segment mean ($|\bar{y}_{S_i}|$) of greater than 0.2. 

Biologically, this value can be thought to quantify numerical chromosomal instability. This is also a simpler representation of how much of the genome has been altered, and it does not run into the issue of sequencing coverage affecting the fraction of the genome altered.

$$
Number\ of\ Altered\ Bases = \sum^{R}_{i = 1} d_i\ where\ |\bar{y}_{S_i}| \ge 0.2
$$
```{r}
base.seg.test <- countingBaseSegments(cnvData = maskCNV_BRCA)
head(base.seg.test)
```

## Counting Number of Break Points

*countingBreakPoints* calculates the number of break points defined as the number of segments ($n_i$) with an absolute segment mean greater than 0.2. This is then doubled to account for the 5' and 3' break points. 

Biologically, this value can be thought to quantify structural chromosomal instability. 

$$
Number\ of \ Break\ Points = \sum^{R}_{i = 1} (n_i \cdot 2)\ where\ |\bar{y}_{S_i}| \ge 0.2
$$

```{r}
break.points.test <- countingBreakPoints(cnvData = maskCNV_BRCA)
head(break.points.test)
```

## Fraction of Genome Altered

*fga* calculates the fraction of the genome altered (FGA; [Chin SF, et. al.](https://pubmed.ncbi.nlm.nih.gov/17925008/)), measured by taking the sum of the number of bases altered and dividing it by the genome length covered ($G$). Genome length covered was calculated by summing the lengths of each probe on the Affeymetrix 6.0 array. This calculation **excludes** sex chromosomes.

$$
Fraction\ Genome\ Altered = 
\frac
{\sum^{R}_{i = 1} d_i}
{G}
\ \ where\  |\bar{y}_{S_i}| \ge 0.2
$$

```{r}
fraction.genome.test <- fga(cnvData = maskCNV_BRCA)
head(fraction.genome.test)
```

## CINmetrics
*CINmetrics* calculates tai, cna, number of altered base segments, number of break points, and fraction of genome altered and returns them as a single data frame. 

```{r}
cinmetrics.test <- CINmetrics(cnvData = maskCNV_BRCA)
head(cinmetrics.test)
```


---
title: "Manuscript_figures"
author: "Vishal Oza"
date: "8/5/2021"
output: html_document
---

```{r}
library(CINmetrics)
library(TCGAbiolinks)
library(tidyverse)
library(ggpubr)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(grid)
library(gridExtra)
```

```{r}
## Downloading CNV BRCA data from TCGA and preparing input data frame
query.maskCNV.hg39.BRCA <- GDCquery(project = "TCGA-BRCA",
              data.category = "Copy Number Variation",
              data.type = "Masked Copy Number Segment", legacy=FALSE)
GDCdownload(query = query.maskCNV.hg39.BRCA)
maskCNV.BRCA <- GDCprepare(query = query.maskCNV.hg39.BRCA, summarizedExperiment = TRUE)
maskCNV.BRCA <- data.frame(maskCNV.BRCA, stringsAsFactors = FALSE)
```



```{r}
## Running CINmetrics on BRCA data frame
cinmetrics.test <- CINmetrics(cnvData = maskCNV.BRCA)
mod.tai <- taiModified(cnvData = maskCNV.BRCA)
cinmetrics.test <- inner_join(cinmetrics.test, mod.tai, by = c("sample_id"), copy = TRUE)
head(cinmetrics.test)
cinmetrics.data <- cinmetrics.test %>%
  pivot_longer(!sample_id, names_to = "cinmetric", values_to = "value") %>%
  separate(sample_id, c("patient", "type"), 12) %>%
  separate(type, c("tumour"), -13)

cinmetrics.data$tumour <-  as.numeric(str_remove(cinmetrics.data$tumour, "-"))
cinmetrics.data <- cinmetrics.data %>% 
  mutate(type = case_when(tumour == 1 | tumour == 6 ~ "tumor",
                          tumour == 10 | tumour == 11 ~ "normal"))

```


```{r}
## Plotting CINmetrics calculated for BRCA

figure1 <- cinmetrics.data %>% mutate(value = log10(value)) %>% 
  ggstripchart(., x = "cinmetric",
      y = "value",
      combine = TRUE,
      #merge = TRUE,
      ylab = "",
      xlab = "",
      color = "type", 
      size = 0.5,
      fill = "type",
      font.label = c(10, "bold", "black"),
      #shape = "type",
      palette = "simpsons",
      alpha = 0.5,
      add = "mean_sd", add.params = list(size = 0.1, alpha = 1, group = "type", color = "black"), 
      orientation = "horiz",
      order = c("fga", "base_segments", "break_points", "cna", "modified_tai", "tai"),
      position = position_jitterdodge() 
      ) + font("xy.text", size = 10, color = "black", face = "bold")

figure2 <- ggpar(figure1, main = "CINmetrics distribution in BRCA", ylab = "log10(metric)", font.main = c(12, "bold", "black"), font.y = c(12, "bold", "black"), font.x = c(12, "bold", "black"), legend = c(0.85,0.85), legend.title = "Sample type", font.legend = c(12, "bold", "black")) #%>% ggexport(filename = "figure2.pdf")
```

```{r}

## Calculating correlation between the CINmetrics and plotting the heatmap
cin.corr.df <- cinmetrics.test %>%
  separate(sample_id, c("patient", "type"), 12) %>%
  separate(type, c("tumour"), -13)

cin.corr.df$tumour <-  as.numeric(str_remove(cin.corr.df$tumour, "-"))
cin.corr.df <- cin.corr.df %>% 
  mutate(type = case_when(tumour == 1 | tumour == 6 ~ "tumor",
                          tumour == 10 | tumour == 11 ~ "normal"))

col_simpson <- colorRamp2(c(-1, 0, 1), c("#f8db27", "#ffffff", "#2f64d6"))

ht_opt(heatmap_column_names_gp = gpar(fontface = "bold", fontsize = 12),
       heatmap_row_names_gp = gpar(fontface = "bold", fontsize = 12),
       heatmap_column_title_gp = gpar(fontface = "bold", fontsize = 12),
       legend_border = "black", legend_labels_gp = gpar(fontface = "bold", fontsize = 10),
       heatmap_border = TRUE)
normal.heatmap <- cin.corr.df %>% filter(type == "normal") %>% select(3:8) %>% cor(., method = "spearman") %>% Heatmap(name = "rho", col = col_simpson, column_title = "Normal samples", row_names_side = "left", show_row_dend = TRUE)
tumor.heatmap <- cin.corr.df %>% filter(type == "tumor") %>% select(3:8) %>% cor(., method = "spearman") %>% Heatmap(col = col_simpson, show_heatmap_legend = FALSE, column_title = "Tumour samples", show_row_names = FALSE, row_dend_side = "right")
combined.heatmap <- normal.heatmap + tumor.heatmap
final.heatmap <- grid.grabExpr(draw(combined.heatmap, column_title = "Spearman correlation between CINmetrics", row_dend_side = "right", column_title_gp = gpar(fontface = "bold", fontsize = 12)))

## Plotting and saving figure1 in pdf format
pdf("Fig1.pdf", width = 14, height = 6) # Open a new pdf file
plot_grid(figure2, final.heatmap, nrow = 1, ncol = 2, align = "hv", axis = "tb",
          rel_widths = c(1,1), labels = c('A','B'))
dev.off() # Close the file

## Plotting and saving figure1 in jpeg format
jpeg("Fig1.jpeg", width = 12, height = 6, units = 'in', res = 350) # Open a new pdf file
plot_grid(figure2, final.heatmap, nrow = 1, ncol = 2, align = "hv", axis = "tb",
          rel_widths = c(1,1), labels = c('A','B'))
dev.off() # Close the file

```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{fga}
\alias{fga}
\title{Fraction Genome Altered}
\usage{
fga(cnvData, segmentMean = 0.2, numProbes = NA, genomeSize = 2873203431)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{segmentMean}{numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2}

\item{numProbes}{Number of Probes}

\item{genomeSize}{Size of the genome derived from Affymetrix 6.0 array probe. Default is 2873203431 calculated based on hg38 **excluding sex chromosomes**}
}
\value{
Fraction of the genome altered
}
\description{
Fraction Genome Altered looks at the fraction of the genome that deviates from a diploid state
fga calculates the fraction of the genome altered (FGA; [Chin SF, et. al.](https://www.ncbi.nlm.nih.gov/pubmed/17925008)), measured by taking the sum of the number of bases altered and dividing it by the genome length covered ($G$). Genome length covered was calculated by summing the lengths of each probe on the Affeymetrix 6.0 array. This calculation **excludes** sex chromosomes.
\deqn{
Fraction\ Genome\ Altered =
\frac
{\sum^{R}_{i = 1} d_i}
{G}
\ \ where\  |\bar{y}_{S_i}| \ge 0.2
}
}
\examples{
fga(cnvData = maskCNV_BRCA)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{cna}
\alias{cna}
\title{Copy Number Aberration}
\usage{
cna(
  cnvData,
  segmentMean = (log(1.7, 2) - 1),
  numProbes = NA,
  segmentDistance = 0.2,
  minSegSize = 10
)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{segmentMean}{numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2}

\item{numProbes}{Number of Probes}

\item{segmentDistance}{Segment distance threshold}

\item{minSegSize}{Minimum segment size}
}
\value{
Number of copy number aberrations between segments
}
\description{
Calculates the number of copy number aberrations
}
\details{
Copy Number Aberrations (CNA) \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079079}{(Davidson JM, et al)}, are defined as a segment with copy number outside the pre-defined range of 1.7-2.3 \deqn{(\log_2 1.7 -1) \le \bar{y}_{S_i} \le (\log_2 2.3 -1)} that is not contiguous with an adjacent independent CNA of identical copy number. For our purposes, we have adapted the range to be \deqn{|\bar{y}_{S_i}| \ge |\log_2 1.7|}, which is only slightly larger than the original.
It is nearly identical to countingBreakPoints, except this one calculates breaks as adjacent segments that have a difference in segment means of \eqn{\ge 0.2}.
\deqn{Total\ Copy\ Number\ Aberration = \sum^{R}_{i = 1} n_i \ where \
\bar{y}_{S_i}| \ge |\log_2{1.7}|, \
\bar{y}_{S_{i-1}} - \bar{y}_{S_i}| \ge 0.2, \
d_i \ge 10}
}
\examples{
cna(cnvData = maskCNV_BRCA)
}
\seealso{
\code{\link{countingBreakPoints}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{CINmetrics}
\alias{CINmetrics}
\title{CINmetrics}
\usage{
CINmetrics(
  cnvData,
  segmentMean_tai = 0.2,
  segmentMean_cna = (log(1.7, 2) - 1),
  segmentMean_base_segments = 0.2,
  segmentMean_break_points = 0.2,
  segmentMean_fga = 0.2,
  numProbes = NA,
  segmentDistance_cna = 0.2,
  minSegSize_cna = 10,
  genomeSize_fga = 2873203431
)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{segmentMean_tai}{numerical value for the minimum segment_mean cutoff/ threshold for Total Aberration Index calculation. Default is 0.2}

\item{segmentMean_cna}{numerical value for the minimum segment_mean cutoff/ threshold for Copy Number Aberration calculation. Default is 0.2}

\item{segmentMean_base_segments}{numerical value for the minimum segment_mean cutoff/ threshold for Base segments calculation. Default is 0.2}

\item{segmentMean_break_points}{numerical value for the minimum segment_mean cutoff/ threshold for Break points calculation. Default is 0.2}

\item{segmentMean_fga}{numerical value for the minimum segment_mean cutoff/ threshold for Fraction of genome altered calculation. Default is 0.2}

\item{numProbes}{Number of Probes}

\item{segmentDistance_cna}{Segment distance threshold}

\item{minSegSize_cna}{Minimum segment size}

\item{genomeSize_fga}{Size of the genome derived from Affymetrix 6.0 array probe. Default is 2873203431 calculated based on hg38 **excluding sex chromosomes**}
}
\value{
All Chromosomal INstability metrics
}
\description{
Calculate all CINmetrics on a given dataframe
}
\examples{
CINmetrics(cnvData = maskCNV_BRCA)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{countingBreakPoints}
\alias{countingBreakPoints}
\title{countingBreakPoints}
\usage{
countingBreakPoints(cnvData, segmentMean = 0.2, numProbes = NA)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{segmentMean}{numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2}

\item{numProbes}{Number of Probes}
}
\value{
Number of Break points for each unique sample
}
\description{
The Break Point calculation takes all the CNV data for a single patient and first filters it for segmentation mean of > 0.2 and, if specified, the minimum number of probes
covering that area. Then it counts the number of rows of data and multiplies it by 2. This represents the break points at the 5' and 3' ends of each segment.
\deqn{
Number\ of \ Break\ Points = \sum^{R}_{i = 1} (n_i \cdot 2)\ where\ |\bar{y}_{S_i}| \ge 0.2
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{taiModified}
\alias{taiModified}
\title{Modified Total Aberration Index}
\usage{
taiModified(cnvData, segmentMean = 0, numProbes = NA)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{segmentMean}{numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2}

\item{numProbes}{Number of Probes}
}
\value{
Average of lengths weighted by segmentation mean for each unique sample
}
\description{
Modified Total Aberration Index calculation takes the sum of lengths of each segment
times its segmentation mean for each sample and divides it by the sum of the
lengths of each sample.
}
\details{
Modified Total Aberration Index uses all sample values instead of those in aberrant copy number state, thus does not remove the directionality from the score.
\deqn{
Modified\ Total\ Aberration\ Index =
\frac
{\sum^{R}_{i = 1} {d_i} \cdot {\bar{y}_{S_i}}}
{\sum^{R}_{i = 1} {d_i}}
}
}
\examples{
taiModified(cnvData = maskCNV_BRCA)
}
\seealso{
\code{\link{tai}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{countingBaseSegments}
\alias{countingBaseSegments}
\title{countingBaseSegments}
\usage{
countingBaseSegments(cnvData, segmentMean = 0.2, numProbes = NA)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{segmentMean}{numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2}

\item{numProbes}{Number of Probes}
}
\value{
Number of Base segments for each unique sample
}
\description{
Function for counting altered base segments
}
\details{
The Altered Base Segment calculation takes all the CNV data for a single patient and first filters it for a segmentation mean of > 0.2 and, if specified, the minimum number of probes
covering that area. Then, it calculates the sums of the lengths of each segment for a particular patient and outputs that.
\deqn{
Number\ of\ Altered\ Bases = \sum^{R}_{i = 1} d_i\ where\ |\bar{y}_{S_i}| \ge 0.2
}
}
\examples{
countingBaseSegments(cnvData = maskCNV_BRCA)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{tai}
\alias{tai}
\title{Total Aberration Index}
\usage{
tai(cnvData, segmentMean = 0.2, numProbes = NA)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{segmentMean}{numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2}

\item{numProbes}{Number of Probes}
}
\value{
Average of lengths weighted by segmentation mean for each unique sample
}
\description{
Total Aberration Index calculation takes the sum of lengths of each segment
times its segmentation mean for each sample and divides it by the sum of the
lengths of each sample.
}
\details{
The Total Aberration Index (TAI) \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3553118/}{(Baumbusch LO, et. al.)} is ``a measure of the abundance of genomic size of copy number changes in a tumour".
It is defined as a weighted sum of the segment means
\deqn{
Total\ Aberration\ Index =
 \frac
{\sum^{R}_{i = 1} {d_i} \cdot |{\bar{y}_{S_i}}|}
{\sum^{R}_{i = 1} {d_i}}\ \
where |\bar{y}_{S_i}| \ge |\log_2 1.7|
}
}
\examples{
tai(cnvData = maskCNV_BRCA)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{maskCNV_BRCA}
\alias{maskCNV_BRCA}
\title{Breast Cancer Data from TCGA
Data Release 25.0
GDC Product: Data
Release Date: July 22, 2020
Masked Copy Number variation data for Breast Cancer for 10 unique samples selected randomly from TCGA}
\format{
An object of class dataframe
}
\source{
\url{https://portal.gdc.cancer.gov/}
}
\usage{
data(maskCNV_BRCA)
}
\description{
Breast Cancer Data from TCGA
Data Release 25.0
GDC Product: Data
Release Date: July 22, 2020
Masked Copy Number variation data for Breast Cancer for 10 unique samples selected randomly from TCGA
}
\examples{
data(maskCNV_BRCA)
\donttest{tai <- tai(maskCNV_BRCA)}
}
\references{
Koboldt, D., Fulton, R., McLellan, M. et al. (2012) Nature 490, 61–70
\url{https://www.nature.com/articles/nature11412}
}
\keyword{dataset}
