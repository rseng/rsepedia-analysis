# Metabolite data preparation pipeline - version 1.0

by: David Hughes and Laura Corbin	
date: June 3rd 2019 

<img src="inst/figures/metaboprep.png" alt="hex" width="150"/> 

## This package
1. Reads in and processes (un)targeted metabolite data, saving datasets in tab-delimited format for use elsewhere
2. Provides useful summary data in the form of tab-delimited text file and a html report.
3. Performs data filtering on the data set using a standard pipeline and according to user-defined thresholds.

## Install metaboprep
1. To install do the following
	
	 1. quick install
		1. start an R session
		2. install the metaboprep package with
			
			```R
			devtools::install_github("MRCIEU/metaboprep")
			```
			
		3. from this repo download a copy of the following files
			1. run_metaboprep_pipeline.R
			2. parameter_file.txt
			
			* You can also download or clone the entire repo with
				
				```R
				git clone https://github.com/MRCIEU/metaboprep.git
				```
				
	2. alternatively you can download the package manually
		1. download a copy of the depository
		2. unzip/pack the download
		3. place the directory somewhere sensible
		4. start an R session
		5. set your working directory to the parent directory of the repo
		6. install R package with: 
			
			```R
			devtools::install("metaboprep")
			```
			
	3. A common installation error is produced by installation errors of dependent packages. If you experience this, install those dependent packages manually with BiocManager, and then attempt the installation of metaboprep again. You might have to repeat this step more than once. 

		```R
		if (!requireNamespace("BiocManager", quietly = TRUE))
	   install.packages("BiocManager")
		BiocManager::install("MISSINGPACKAGENAME")
		```


## To run metaboprep

1. Edit the paramater (parameter_file.txt) file
	1.	do not add any spaces before or after the "=" sign in the paramater file.
	2. the paramater file can be located anywhere
2. Move to the, or a, directory containing a copy of:
	1. run_metaboprep_pipeline.R
3. Make sure that R is in your environment - an often necessary step if working on an HPC.
	1. for example: module add languages/R-3.5-ATLAS-gcc-7.1.0
4. Run the metaboprep pipeline on a terminal command line as follows:
	
	```R	
	Rscript run_metaboprep_pipeline.R /FULL/PATH/TO/example_data/excel/parameter_file.txt
	```
	
5. If you experienced issues with the geneartion of the html report on an HPC, the report can be generated on a local machine as follows.
	1. move to your newly generated metaboprep project directory. 
		* it will take the form of "../metaboprep_release_TODAYSDATE/"
	2. You should find an R data object called "ReportData.Rdata". Save a copy locally.
	3. Open an R session
	4. produce report with the function metaboprep::generate_report() as
		
		```R
		output_dir_path = paste0("FULL/PATH/TO/DIR/OF/CHOICE/")
		rdfile = paste0(output_dir_path, "ReportData.Rdata")
		generate_report( full_path_2_Rdatafile = rdfile, dir_4_report = output_dir_path )

		```

## Example data

An example data set can be found in the folder "example_data" here on the repository. It is a simulated data set of 100 metabolites, for 100 samples. There is a (1) metabolon like (v1 format) excel file of the data set, and a (2) flat text (tab delim) version of the data set. Both are accompanied by a parameter file to help guide you with the example. This example data includes data from two hypothetical mass spectrometry run modes or platforms, "neg" and "pos". As such, a subset of the metabolites were run (simulated) in the "neg" run mode in two batches, and the second subset of metabolites were run (simulated) in the "pos" run mode in three batches. Each batch was simulated with different mean abundance values to help illustrate possible batch effects and the normalization procedure. If looking at the flat text version of the example, information on which metabolites were simulated in which run mode or platform can be found in the flat text file "feature_data.txt" in the column "platform". Further, to identify which samples belonged to which batch, for each run mode, you would use the columns "neg" and "pos" in the flat text file "sample_data.txt".

## Data Preparation steps in brief

![](images/metaboprep_workflow.png)

## See the Wiki page for a detailed synopsis of the pipeline [here](https://github.com/MRCIEU/metaboprep/wiki/Detailed-metaboprep-pipeline-steps).

### (A) General Outline of metaboprep
1. Read in the paramater file
2. Read in the data  -  *(typically from a commercially provided excel file)*
	* metabolite abundance
	* sample annotation
	* feature annotation
3. Write metabolite data, sample annotation and feature annotation to flat text file.
4. Normalize data
	- If data is from Metabolon or is any other technology that has run-mode or platform batches, median normalize the data. 
5. Estimate summary statistics on the raw data set **(step B below)**
	* write summary stats to file
6. Perfom the data filtering **(step C below)**
	* using parameters passed in the parameter file 
	* write data filtering (metaboprep) data set to file
7. Estimate sumary statistics on the filtered data set **(step B below)**
	* write summary stats to file
8. Generate html report
9. Print scatter plot, histogram, and summary stats for every metabolite to a single PDF.

### (B) Summary Statistics Estimated
1. Sample Summary Statistics
	* sample missingness
		+ all features
		+ to the exclusion of xenobiotic and\or derived variables
	* sample total sum abundance (TSA) **(derived variables excluded)**
		+ with all features 
		+ with complete features only (no missingness) 
	* count of how many times a sample is an outlier across all feature
		+ each feature analyzed within its own sample distribution
		+ outliers determined as those +/- 5 IQR of the median.
2. Feature Summary Statistics
	* feature missingness
		+ all samples
		+ to the exclusion of sample(s) with extreme missingness (>= 50%)
	* distribution statistics
		+ shaprio's W-statistic of normality on raw distribution
		+ shaprio's W-statistic of normality on log10 distribution
		+ skewness
		+ kutosis
		+ N, sample size
		+ variance
		+ standard deviation
		+ coefficient of variation
		+ mean
		+ median
		+ min
		+ max
	* feature outlying sample count
		+ count of outlying samples 
			+ outliers determined as those +/- 5 IQR of the median.
3. Feature and Sample structure
	* feature:feature correlation structure **(derived variables excluded)**
		+ only includes features with at least 50 measurments
			+ or if the data set has an N<50 the missingness allowed is 0.8 * N
		+ estimate the number of independent features
		+ tag representitive features of feature clusters
	* sample:sample correaltion structure **(derived variables excluded)**
		+ principle components (PCA)
			* derived from independent features
			* missing data is imputed to the median estiamte of each feature
			* identify PC outliers
				* +/- 3,4,5 SD of mean for all significant PCs

### (C) Data Filtering Steps
1. If data is from Metabolon, exclude (but retain for step 11) xenobiotic metabolites from anlaysis.
2. Estimate sample missingness and exclude extreme samples, those with a missingness >= 0.80 (or 80%) **(derived variables excluded)**
3. Estimate feature missingness and exclude extreme features, those with a missingness >= 0.80 (or 80%)
4. Re-estimate sample missingness and exclude samples >= user defined threshold (units: 0.2 or 20% missing) **(derived variables excluded)**
5. Re-estimate feature missingness and exclude features >= user defined threshold (units: 0.2 or 20% missing)
6. Estimate total sum abundance/area (the sum of all values) for each individual 
	* first z-transform each distribution
	* second shift the mean of each distribution to the absolute minimum of ALL observed values
7. TSA sample exclusion using a user defined threshold (units: +/- SD from mean)  **(derived variables excluded)**
	* To ignore this step set to NA in parameter file
8. Identify outlier values for each feature using a user defined threshold (units: +/- IQR from median) 
9. User define what to do with outliers **for the purposes of deriving the PCs only**
	* "leave_be" which means the outlier values will be left as they are.
	* "turn_NA" which means they will be median imputed for the PCA.
	* "winsorize" to the 100th quantile of all remaining (non-outlier values).
10. Build feature:feature correlation matrix on filtered data derived from steps 1-6 above **(derived variables excluded)**
	* To be included a feature must have a minimun of 50 observations, or N*0.8 observations if data set includes less than 50 individuals.
11. Identify "independent" features using data from step 8 and user defined tree cut height.
	* we retain the feature with the least missingness within a cluster, otherwise we select the first feature in the list. 	
12. Estimate principal components using independent features from step 9
13. PC outlier samples exclusions using user defined threshold (units: +/- SD from the mean)
14. If the data is from Metabolon we place the xenobiotic metabolites back into the filtered data set. 


**NOTE: Derived variable are those that are ratios or percentanges of two or more features already present in a data set, such as those found in Nightingale data.**

## HTML Report includes

*---example figures provided for illustration---*

1. General information on study
2. Raw data summary
	* sample size
	* missingness
![](images/missingness_matrix.png)
![](images/missingness_hist.png)
	* table of data filtered exclusions
	* figure of PCA oulier exclusions
![](images/outlier_pca.png)
3. Summary of data filtered data
	* sample size
	* summary figures
		+ missingness distributions
		+ feature dendrogram
		+ PCA
![](images/summary_of_qc_data.png)
	* distribtion of Shapiro W-statistics
![](images/wstats.png)
	* outlier summary
4. Batch effects
	* feature missingness
		+ as influenced by:
			+ feature SUPER_PATHWAY (categorical function)
			+ MS method
				+ LC/MS Polar, Pos Late, Pos Early, Neg
![](images/f_missingness_by_batch.png)
	* sample missingness
		+ as influenced by:
			+ BOX_ID, storage box
			+ RUN_DAY, day the samples were processed on tech
![](images/missingness_by_batch.png)
		+ multivariate analysis of both BOX_ID and RUN_DAY on missingness
	* sample total peak|abundance area
		+ as influenced by:
			+ BOX_ID, storage box
			+ RUN_DAY, day the samples were processed on tech
![](images/tsa_by_batch.png)
		+ multivariate analysis of both BOX_ID and RUN_DAY on missingness
5. Power Analysis
	* presence -vs- absence
	* continuous
![](images/power_v2.png)# metaboprep data preparation summary report

### Date: `r format(Sys.time(), '%d %B, %Y')`

metaboprep report relates to:

  * Project: `r project`  
  * Platform: `r platform`

The `metaboprep` `R` package performs three operations: 

1. Provides an assessment and summary statistics of the raw metabolomics data.
2. Performs data filtering on the metabolomics data.
3. Provides an assessment and summary statistics of the filtered metabolomics data, particularly in the context of batch variables when available.

This report provides descriptive information for raw and filtered metabolomics data for the project `r project`. 


```{r init, echo = FALSE}
knitr::opts_chunk$set(echo = FALSE, quote=F, comment=NA, warning=FALSE, message=FALSE, error=FALSE, fig.align="center"  )

## test for necessary packages
if (!requireNamespace("metaboprep", quietly = TRUE)) {
    stop("Package \"metaboprep\" needed for this function to work. Please install it.",
      call. = FALSE)
}

if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("Package \"dendextend\" needed for this function to work. Please install it.",
      call. = FALSE)
}


if (!requireNamespace("tidyverse", quietly = TRUE)) {
    stop("Package \"tidyverse\" needed for this function to work. Please install it.",
      call. = FALSE)
}

suppressPackageStartupMessages(library(metaboprep))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dendextend))

```


```{r define_data}
## RAW DATA
raw_data$metabolite_data = tibble::as_tibble(raw_data$metabolite_data)
raw_data$sample_data = tibble::as_tibble(raw_data$sample_data); colnames(raw_data$sample_data)[1] = "SamID"
raw_data$feature_data = tibble::as_tibble(raw_data$feature_data); colnames(raw_data$feature_data)[1] = "FeatureID"
raw_data$varexp = raw_data$varexp[,1]

## QC Data
qc_data$metabolite_data = tibble::as_tibble(qc_data$metabolite_data)
qc_data$sample_data = tibble::as_tibble(qc_data$sample_data); colnames(qc_data$sample_data)[1] = "SamID"
qc_data$feature_data = tibble::as_tibble(qc_data$feature_data); colnames(qc_data$feature_data)[1] = "FeatureID"
qc_data$varexp = qc_data$varexp[,1]

```


## The data preparation workflow is as follows:

```{r png_fig, out.width = "400px"}
f = system.file("rmarkdown", package="metaboprep")
pic = file.path( f, "skeleton/metaboprep_workflow.png")
knitr::include_graphics(pic)
```

1. Issues can be raised on [GitHub](https://github.com/MRCIEU/metaboprep/issues). 
2. Questions relating to the `metaboprep` pipeline can be directed to [David Hughes: d.a.hughes@bristol.ac.uk](mailto:d.a.hughes@bristol.ac.uk).
3. `metaboprep` is published in [Journal to be determined]() and can be cited as:


# 1. Summary of raw data

## Sample size of `r project` data set

```{r Sample_size, fig.width = 4, fig.height = 2}
tout = data.frame("data.set" = c("number of samples","number of features"),  
                  "raw.data" = dim(raw_data$metabolite_data), 
                  "filtered.data" = dim(qc_data$metabolite_data))
ggpubr::ggtexttable(tout, rows = NULL, theme = ggpubr::ttheme("mBlue") )

```


## Missingness
Missingness is evaluated across samples and features using the original/raw data set. 

### Visual structure of missingness in your raw data set.

```{r MissingnessMatrix, fig.width = 20, fig.height = 12,  echo = FALSE,warning=FALSE, message=FALSE, error=FALSE, dev = "png"}
namat = qc_data$metabolite_data
namat[!is.na(namat)] = 1
namat[is.na(namat)] = 0
namat = as.matrix(namat)

namat = reshape::melt(namat)
colnames(namat) = c("individuals","metabolites","missingness")

pcol = RColorBrewer::brewer.pal(n = 8, name = 'Dark2')
ggplot(namat, aes(metabolites, individuals, fill = missingness)) + 
  geom_tile() + 
  scale_fill_gradient(low= "white", high=pcol[5]) +
  theme(legend.position = "none", axis.text.x = element_blank() )
```
Figure Legend: Missingness structure across the raw data table. White cells depict missing data. Individuals are in rows, metabolites are in columns. 


### Summary of sample and feature missingness 

```{r run.missingness, echo = FALSE, fig.width = 8, fig.height = 6, warning=FALSE, message=FALSE, error=FALSE}
r_mis = missingness.sum(mydata = raw_data$metabolite_data)

p = ggpubr::ggarrange(  r_mis[[4]][[1]] ,
                    r_mis[[4]][[2]] , 
                    r_mis[[4]][[3]] ,
                    r_mis[[4]][[4]] ,
                    ncol = 2, nrow = 2,
                    labels = c("a", "b", "c", "d") )
ggpubr::annotate_figure(p,top = paste0("-- Initial raw data set: Estimates of missingness for samples and features --\n") )


```
Figure Legend: Raw data - (a) Distribution of sample missingness with sample mean illustrated by the red vertical line. (b) Distribution of feature missingness sample mean illustrated by the red vertical line. (c) Table of sample and feature missingness percentiles. A tabled version of plot a and b.  (d) Estimates of study samples sizes under various levels of feature missingness.

# 2. Data Filtering 

## Exclusion summary

```{r exclusion_table, echo = FALSE , fig.width = 5, fig.height = 3, warning=FALSE, message=FALSE, error=FALSE}
temp = data.frame(exclusions = rownames(qcing_data$exclusion_data), count = qcing_data$exclusion_data[,1])
ggpubr::ggtexttable( temp, rows = NULL, theme = ggpubr::ttheme("mBlue"))
```
Table Legend: Six primary data filtering exclusion steps were made during the preparation of the data.
(1) Samples with missingness >=80% were excluded. (2) features with missingness >=80% were excluded (xenobiotics are not included in this step). (3) sample exclusions based on the user defined threshold were excluded. (4) feature exclusions based on user defined threshold were excluded (xenobiotics are not included in this step). (5) samples with a total-peak-area or total-sum-abundance that is >= N standard deviations from the mean, where N was defined by the user, were excluded. (6) samples that are >= N standard deviations from the mean on principal component axis 1 and 2, where N was defined by the user, were excluded.


## Metabolite or feature reduction and principal components

A data reduction was carried out to identify a list of representative features for generating a sample principal component analysis. This step reduces the level of inter-correlation in the data to ensure that the principal components are not driven by groups of correlated features.

```{r PCA_1_2, echo = FALSE, fig.width = 12, fig.height = 6, warning=FALSE, message=FALSE, error=FALSE}
#################
## data reduction 
## info
#################

feature_count = length(qcing_data$feature_sumstats$k)
features_included_in_data_reduction = length( na.omit(qcing_data$feature_sumstats$k) )
cluster_count = length( unique(qcing_data$feature_sumstats$k) )
number_of_rep_meatbolites = sum( qcing_data$feature_sumstats$independent_features_binary == 1 )

temp = data.frame( data.reduction = c("total metabolite count",
                                      "metabolites included in data reduction",
                                      "number of metabolite clusters",
                                      "number of representative metabolites"),
                   count = c(feature_count, 
                             features_included_in_data_reduction, 
                             cluster_count, 
                             number_of_rep_meatbolites  ) )

ptable = ggpubr::ggtexttable( temp, rows = NULL, theme = ggpubr::ttheme("mBlue"))

#################
## PCA
#################
## define tibble
pcs = as_tibble(qcing_data$pcs)
## define color scheme
pcol = RColorBrewer::brewer.pal(9, "Set1")
## define accelerationfactor
accelerationfactor = qcing_data$accelerationfactor
if(accelerationfactor > 10){accelerationfactor = 10}
if(accelerationfactor == 1){accelerationfactor = 2}
## define outliers
Omat = outlier.matrix( pcs[, 1:accelerationfactor], nsd = PC_outlier_SD, meansd = TRUE )
outliers = apply(Omat, 1, sum)
outliers[outliers>0]=1
pcs$outliers = as.factor( outliers )
##
cutoffs = apply(pcs[, 1:accelerationfactor], 2, function(x){
  msd = c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  cutoff = msd[1] + (msd[2]*PC_outlier_SD)
  return(cutoff)
} )
  
##
pcplot = pcs %>% ggplot( aes(x = PC1, y = PC2) ) +
  geom_point( aes(fill = outliers), size = 2.5, shape = 21 ) +
  scale_fill_manual(values = pcol[c(2,1)]) +
  geom_vline(xintercept = c(cutoffs[1], -cutoffs[1]), color = pcol[1] ) +
  geom_hline(yintercept = c(cutoffs[2], -cutoffs[2]), color = pcol[1] ) +
  theme(legend.position="bottom") +
  labs(title = paste0("Principal components 1-&-2 using ",number_of_rep_meatbolites ," representative metabolites"),
       subtitle  = paste0(" - Outliers are those ", PC_outlier_SD , " SD from the mean of PCs 1 to", accelerationfactor ))
##

gridExtra::grid.arrange( grobs = list(ptable, pcplot), widths = c(2, 3), ncol = 2) 

```
Figure Legend: The data reduction table on the left presents the number of metabolites at each phase of the data reduction (Spearman's correlation distance tree cutting) analysis. On the right principal components 1 and 2 are plotted for all individuals, using the representitve features identified in the data reduction analysis. The red vertical and horizontal lines indicate the standard deviation (SD) cutoffs for identifying individual outliers, which are plotted in red. The standard deviations cuttoff were defined by the user.



# 3. Summary of filtered data

## Sample size (N) 
  * The number of samples in data = `r nrow(qc_data$metabolite_data)`  
  * The number of features in data = `r ncol(qc_data$metabolite_data)`  

## Relative to the raw data
  * `r nrow(raw_data$metabolite_data) - nrow(qc_data$metabolite_data)` samples were filtered out, given the user's criteria.
  * `r ncol(raw_data$metabolite_data) - ncol(qc_data$metabolite_data)` features were filtered out, given the user's criteria.
  * Please review details above and your log file for the number of features and samples excluded and why. 
  

```{r sample_overview, echo = FALSE, fig.width = 15, fig.height = 15, warning=FALSE, message=FALSE, error=FALSE}
pcol = RColorBrewer::brewer.pal(9, "Set1")
##############
## Missingness
##############
qc_mis = missingness.sum(mydata = qc_data$metabolite_data)
s_mis_plot = qc_mis[[4]][[1]] 
f_mis_plot = qc_mis[[4]][[2]] 

##############
## TPA
##############
s = tibble::tibble(TPA_completefeature = qc_data$sample_data$TPA_completefeature)

tpa_plot = qc_data$sample_data %>% ggplot( aes( TPA_completefeature ) ) +
    geom_histogram( fill = pcol[2] , bins = 25) + 
    geom_vline( xintercept = median(qc_data$sample_data$TPA_completefeature), color = pcol[1], size = 1) +
    labs(title = paste0("total abundance of samples\nat complete features only"), x="", y="") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))

##############
## Dendrogram
##############
dend = qc_data$feature_tree %>% as.dendrogram

## create a vector of colors to color your tree labels
w = which( qc_data$feature_data$independent_features_binary == 1)
pv_ids = as.character( unlist( qc_data$feature_data[w,1] ) )
n = labels(dend)
bcol = rep("black", length(n))
w = which(n %in% pv_ids ); bcol[w] = pcol[2]

## redefine elements of dendrogram
dend = dend %>% 
  set("labels_cex", 0.5) %>% 
  set("labels_col", bcol) %>% 
  set("branches_lwd", 0.5) %>%
  set("branches_k_color",  value = bcol)
dend <- as.ggdend(dend)
## plot the dendrogram
tree_plot = dend %>% ggplot() + geom_hline(yintercept = tree_cut_height, color = "coral2")

##############
## Data Reduce
## Table
##############
feature_count = length(qc_data$feature_data$k)
features_included_in_data_reduction = length( na.omit(qc_data$feature_data$k) )
cluster_count = length( unique(qc_data$feature_data$k) )
number_of_rep_meatbolites = sum( qc_data$feature_data$independent_features_binary == 1 )

temp = data.frame( data.reduction = c("total metabolite count","metabolites included in data reduction","number of metabolite clusters","number of representative metabolites"),
                   count = c(feature_count,features_included_in_data_reduction, cluster_count, number_of_rep_meatbolites  )  )
ptable = ggpubr::ggtexttable( temp, rows = NULL, theme = ggpubr::ttheme("mGreen"))
#############
## scree plot
#############
# plot(qc_data$varexp, pch = 21, bg = pcol[2], cex = 2, type = "b", xlab = "PC number", ylab = "variance explained | eigenvalue")
ve = tibble::tibble( data.frame(PC= 1:length(qc_data$varexp), VarExp = qc_data$varexp) )

screeplot = ve %>% ggplot(aes(x = PC, y = VarExp)) +
  geom_line(color = "grey") +
  geom_point(shape = 21, fill = pcol[2], size = 2) +
  labs(title = "Scree plot") +
  geom_vline(xintercept = qc_data$accelerationfactor, color = pcol[1]) +
  geom_vline(xintercept = unlist( qc_data$nparallel ), color = pcol[3])

##############
## PC Plot
##############
pcs = as_tibble(qc_data$sample_data[, c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
pcs$cluster_k = as.factor( kmeans(pcs[, 1:2], 4)$cluster )
## define color scheme
pcol = RColorBrewer::brewer.pal(9, "Set1")
##
pcplot = pcs %>% ggplot( aes(x = PC1, y = PC2) ) +
  geom_point( size = 2.5, shape = 21, aes(fill = cluster_k) ) +
  scale_fill_manual(values = pcol[1:4]) +
  labs( title = paste0("Principal components 1-&-2 using ", number_of_rep_meatbolites ," representative metabolites"),
        x = paste0("PC1  VarExp = ", signif(qc_data$varexp[1], d = 4)*100, "%" ),
        y = paste0("PC2  VarExp = ", signif(qc_data$varexp[2], d = 4)*100 , "%"),
        fill = paste0("kmeans\ncluster k"))
#######################
## Plot the data
#######################
m = matrix(c(1,2,3,4,4,4,5,7,7, 6,7,7), nrow = 4, byrow = TRUE)

gridExtra::grid.arrange( grobs = list(s_mis_plot, f_mis_plot, tpa_plot, tree_plot, ptable, screeplot, pcplot), 
                         layout_matrix = m,
                         heights = c(1,1,1,1) ) 

```

Figure Legend: Filtered data summary. Distributions for sample missingness, feature missingness, and total abundance of samples. Row two of the figure provides a Spearman's correlation distance clustering dendrogram highlighting the metabolites used as representative features in blue, the clustering tree cut height is denoted by the horizontal line. Row three provides a summary of the metabolite data reduction in the table, a Scree plot of the variance explained by each PC and a plot of principal component 1 and 2, as derived from the representative metabolites. The Scree plot also identifies the number of PCs estimated to be informative (vertical lines) by the Cattel's Scree Test acceleration factor (red,  n = `r qc_data$accelerationfactor`) and Parallel Analysis (green, n = `r unlist( qc_data$nparallel )`). Individuals in the PC plot were clustered into 4 kmeans (k) clusters, using data from PC1 and PC2. The kmeans clustering and color coding is strictly there to help provide some visualization of the major axes of variation in the sample population(s).


## Structure among samples

```{r PCApairsplot, echo = FALSE, fig.width = 15, fig.height = 15, warning=FALSE, message=FALSE, error=FALSE}
plotcolors = pcol[pcs$cluster_k]
pcapairs_bymoose( as.matrix(pcs[, 1:5]) , qc_data$varexp, pcol = plotcolors)
```

Figure Legend: A matrix (pairs) plot of the top five principal components including demarcations of the 3rd (grey), 4th (orange), and 5th (red) standard deviations from the mean. Samples are color coded as in the summary PC plot above using a kmeans analysis of PC1 and PC2 with a k (number of clusters) set at 4. The choice of k = 4 was not robustly chosen it was a choice of simplicity to help aid visualize variation and sample mobility across the PCs.

## Feature Distributions

### Estimates of normality: W-statistics for raw and log transformed data

```{r, echo=F, warning=FALSE, message=FALSE, error=FALSE}
wstat = qc_data$feature_data$W

## how many stats were estimated
count = length(wstat)
nacount = sum(is.na(wstat))
remain_count = count - nacount

normcount = sum(wstat >= 0.95, na.rm = TRUE)

```


```{r shapiroW, echo = FALSE, fig.width = 10, fig.height = 5, warning=FALSE, message=FALSE, error=FALSE}
W = wstat
pcol = RColorBrewer::brewer.pal(9, "Set1")
W_log10 = qc_data$feature_data$log10_W

## Plot
par(mfrow = c(1,2), oma = c(2,1,1,1))
hist(W, col = pcol[2], 
     main = paste("Distribution of W statistics on\nRaw Metabolite Abundances"), 
     xlab = "W", 
     ylab = "Frequency",
     cex.main = 0.75)
abline(v = 0.95, col = pcol[1])
##
mtext( paste0(normcount,
              " of the metabolites exhibit distributions\nthat may declared normal, given a W-stat >= 0.95"), 
       cex = 0.75, side = 1, outer = TRUE, line = 0, adj = 0.1, col = pcol[2])
  ### LOG
LogMakesDistributionWorse = c( sum( W_log10 < W , na.rm = TRUE), signif( sum( W_log10 < W , na.rm = TRUE)/length(!is.na(W) ), d = 3)*100)
  ##
hist(W_log10, col = pcol[3], 
     main = paste("Distribution of W statistics on\nlog10 Metabolite Abundances"), 
     xlab = "W", 
     ylab = "Frequency",
     cex.main = 0.75)
mtext( paste0("In ", LogMakesDistributionWorse[1],
              " instances or ", LogMakesDistributionWorse[2], "% of the tested metabolites\nthe log10 data W-stat is < raw data W-stat."), cex = 0.75, side = 1, outer = TRUE, line = 0, adj = 0.9, col = pcol[3])

```

Figure Legend: Histogram plots of Shapiro W-statistics for raw (left) and log transformed (right) data distributions. A W-statistic value of 1 indicates the sample distribution is perfectly normal and value of 0 indicates it is perfectly uniform. Please note that log transformation of the data *may not* improve the normality of your data.

Analysis details: Of the `r count` features in the data `r nacount` features were excluded from this analysis because of no variation or too few observations (n < 40). Of the remaining `r remain_count` metabolite features, a total of `r normcount` may be considered normally distributed given a Shapiro W-statistic >= 0.95.


## Outliers

Evaluation of the number of samples and features that are outliers across the data.


```{r outlier_summary, echo = FALSE, fig.width = 8, fig.height = 2.5, warning=FALSE, message=FALSE, error=FALSE}
Omat = outlier.matrix(qc_data$metabolite_data)
Total_Out_Count = sum(Omat)
sout = apply(Omat, 1, function(x){ sum(x, na.rm = TRUE)  })
fout = apply(Omat, 2, function(x){ sum(x, na.rm = TRUE)  })
sam_out_quantiles = c( quantile(sout, probs = c(0, 0.25, 0.5), na.rm = TRUE), 
                       signif( mean(sout, na.rm = TRUE), digits = 4 ), 
                       quantile(sout, probs = c( 0.75, 1), na.rm = TRUE) )
feat_out_quantiles = c( quantile(fout, probs = c(0, 0.25,0.5), na.rm = TRUE), 
                       signif( mean(fout, na.rm = TRUE), digits = 4 ) , 
                       quantile(fout, probs = c(  0.75, 1), na.rm = TRUE) )
outvals = data.frame(quantile = c("minimum", "25th" ,"median", "mean" ,"75th","max"), samples = sam_out_quantiles, features = feat_out_quantiles )
ggpubr::ggtexttable( outvals, rows = NULL, theme = ggpubr::ttheme("mBlue"))
```
Table Legend: The table reports the average number of outlier values for samples and features in the data set. The minimum, max and other quantiles of the sample and feature distributions are reported as well.

### Notes on outlying samples at each metabolite|feature

There may be extreme outlying observations at individual metabolites|features that have not been accounted for. You may want to:

1. Turn these observations into NAs.
2. Winsorize the data to some maximum value.
3. Rank normalize the data which will place those outliers into the top of the ranked standard normal distribution.
4. Turn these observations into NAs and then impute them along with other missing data in your data set. 


# 4. Variation in filtered data by available variables

### feature missingness

Feature missingness may be influenced by the metabolites' (or features') biology or pathway classification, or your technologies methodology. The figure(s) below provides an illustrative evaluation of the proportion of *feature missigness* as a product of the variable(s) available in the raw data files. 

```{r id_feat_batch_vars}
## what are all of the variables that could be used to evaluate feature effects?
possible_vars = colnames(qc_data$feature_data)
## remove summary stats from pool of possible
w = which(possible_vars == "feature_missingness")
if(length(w) == 1){
  possible_vars = possible_vars[ -c(w:length(possible_vars)) ]  
}

## number of unique units for each possible variable?
count_unique = apply(qc_data$feature_data[,possible_vars], 2, function(x){  length( unique( na.omit(x) ) )  })

## remove those with only one class or count == sample size
r = which(count_unique == 1 | count_unique == nrow(qc_data$feature_data) | count_unique > 96)
if(length(r) > 0){
  possible_vars = possible_vars[-r]
  count_unique = count_unique[-r]
}

## continue filtering
if( length(possible_vars) > 0){
  ## estimate min, mean, max number of values within a variable
  features_per_unit = t( apply( qc_data$feature_data[,possible_vars], 2, function(x){  
    x = table( unlist(x) )
    out = c( min(x), median(x), max(x) ); names(out) = c("min","median","max")
    out
  }) )
  features_per_unit = as.data.frame( cbind( count_unique, features_per_unit) )
  ## filter 2 class min of 1 variables
  r = which(features_per_unit$count_unique == 2 & features_per_unit$min == 1)
  if(length(r) > 0){ features_per_unit = features_per_unit[-r,] }
  
  ## filter for median observational count
  k = which(features_per_unit$median >= 5 )
  possible_vars = possible_vars[k]
}

## define
class_variables = possible_vars


```


```{r}
if(length(class_variables)==0){
  paste0(" -- No feature level batch variables identified or all were invariable -- ") 
}
```


```{r featuremis_1or2,  echo=F, fig.width = 15, fig.height = 10, warning=FALSE, message=FALSE, error=FALSE}
## MISSINGNESS
###################################
## Iterate over each batch variable
###################################
if(length(class_variables) > 0 & length(class_variables)<=2){
  ClassMisPlots = lapply( class_variables , function(x){
  out = variable.by.factor( dep = qc_data$feature_data$feature_missingness , 
                           indep = unlist( qc_data$feature_data[,x] ), 
                           dep_name = "feature missingness", 
                           indep_name = x, orderfactor = TRUE, violin = FALSE )
  return(out)
  })


  ## plot the output
  gridExtra::grid.arrange( grobs = ClassMisPlots, ncol = 1)

}

```


```{r featuremis_3,  echo=F, fig.width = 15, fig.height = 17, warning=FALSE, message=FALSE, error=FALSE}
## MISSINGNESS
###################################
## Iterate over each batch variable
###################################
if(length(class_variables)>2){
  ClassMisPlots = lapply( class_variables , function(x){
  out = variable.by.factor( dep = qc_data$feature_data$feature_missingness ,
                           indep = unlist( qc_data$feature_data[,x] ), 
                           dep_name = "feature missingness", 
                           indep_name = x, orderfactor = TRUE, violin = FALSE )
  return(out)
  })


  ## plot the output
  gridExtra::grid.arrange( grobs = ClassMisPlots, ncol = 1)

}

```
Figure Legend: Box plot illustration(s) of the relationship that feature variables have with feature missingness.


### sample missingness

The figure provides an illustrative evaluation of the proportion of *sample missigness* as a product of sample batch variables provided by your supplier. This is the univariate influence of batch effects on *sample missingness*.

```{r id_sam_batch_vars}
## what are all of the variables that could be used to evaluate feature effects?
possible_vars = colnames(qc_data$sample_data)

## remove summary stats from pool of possible
w = which(possible_vars == "sample_missingness")
if(length(w)>0){
  possible_vars = possible_vars[ -c(w:length(possible_vars)) ]  
}

## number of unique units for each possible variable?
count_unique = apply(qc_data$sample_data[,possible_vars], 2, function(x){  length( unique( na.omit(x) ) )  })

## remove those with only one class or count == sample size
r = which(count_unique == 1 | count_unique == nrow(qc_data$sample_data) | count_unique > 96)
if(length(r) > 0){
  possible_vars = possible_vars[-r]
  count_unique = count_unique[-r]
}

## continue filtering
if( length(possible_vars) > 0){
  ## estimate min, mean, max number of values within a variable
  features_per_unit = t( apply( qc_data$sample_data[,possible_vars], 2, function(x){  
    x = table( unlist(x) )
    out = c( min(x), median(x), max(x) ); names(out) = c("min","median","max")
    out
  }) )
  features_per_unit = as.data.frame( cbind( count_unique, features_per_unit) )
  ## filter 2 class min of 1 variables
  r = which(features_per_unit$count_unique == 2 & features_per_unit$min == 1)
  if(length(r) > 0){ features_per_unit = features_per_unit[-r,] }
  
  ## filter for median observational count
  k = which(features_per_unit$median >= 5 )
  possible_vars = possible_vars[k]
}

## define
batch_variables = possible_vars
```


```{r}
if(length(batch_variables)==0){
  paste0(" -- No feature level batch variables identified or all were invariable -- ") 
}
```


```{r sample_missingness_12,  echo=F, fig.width = 15, fig.height = 10, warning=FALSE, message=FALSE, error=FALSE}
## MISSINGNESS
###################################
## Iterate over each batch variable
###################################
if( length(batch_variables) > 0 & length(batch_variables) <= 2 ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$sample_missingness , 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "sample missingness", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 1)
  
} 

```



```{r sample_missingness_3,  echo=F, fig.width = 15, fig.height = 17, warning=FALSE, message=FALSE, error=FALSE}
## MISSINGNESS
###################################
## Iterate over each batch variable
###################################
if( length(batch_variables) > 2 & length(batch_variables) <= 4 ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$sample_missingness , 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "sample missingness", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 1)
  
} 

```


```{r sample_missingness_4,  echo=F, fig.width = 15, fig.height = 17, warning=FALSE, message=FALSE, error=FALSE}
## MISSINGNESS
###################################
## Iterate over each batch variable
###################################
if( length(batch_variables) > 4 ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$sample_missingness, 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "sample missingness", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 2)
  
} 

```

Figure Legend: Box plot illustration(s) of the relationship that available batch variables have with sample missingness.

## Multivariate evaluation: batch variables

```{r sample_missingness_multivatiaveANOVA, echo = FALSE, fig.width = 9, fig.height = 2, warning=FALSE, message=FALSE, error=FALSE}
if( length(batch_variables) > 0 ) {
  covars = as.data.frame( qc_data$sample_data[ , batch_variables ] )
  for(i in 1:ncol(covars)){ covars[,i] = as.factor(covars[,i]) }
  ( multivariate.anova(dep = qc_data$sample_data$sample_missingness, indep_df = covars ) )
} else {
  paste0(" -- No sample level batch variables were provided or all were invariable -- ")
  }
```

Table Legend: TypeII ANOVA: the eta-squared (eta-sq) estimates are an estimation of the percent of variation explained by each independent variable, after accounting for all other variables, as derived from the sum of squares. This is a multivariate evaluation of batch variables on *sample missingness*. Presence of NA's would indicate that the model is inappropriate.


# 5. Total peak or abundance area (TA) of samples:

The total peak or abundance area (TA) is simply the sum of the abundances measured across all features. TA is one measure that can be used to identify unusual samples given their entire profile. However, the level of missingness in a sample may influence TA. To account for this we:  

1. Evaluate the correlation between TA estimates across all features with PA measured using only those features with complete data (no missingness).  
2. Determine if the batch effects have a measurable impact on TA.


## Relationship with missingness

Correlation between total abundance (TA; at complete features) and missingness

```{r TPA_missingness, echo=F, fig.width = 7, fig.height = 5, warning=FALSE, message=FALSE, error=FALSE}
a = cor.test( qc_data$sample_data$sample_missingness, qc_data$sample_data$TPA_completefeature)
###
( 
  tpamis = qc_data$sample_data %>% ggplot( aes(x = TPA_completefeature, y = sample_missingness)) +
  geom_point( fill = "grey", alpha = 0.8, size = 1.5) + 
  geom_smooth(method = "loess", color = "red", size = 2)  +
  geom_smooth(method = "lm", color = "blue", size = 2)  +
  labs(x = "TA, complete features", y = "sample missingness",
       title = paste0( "TA as influenced by missingness \nSpearmans's cor = ", round(a$estimate, d = 4), " p-value = ", 
                    formatC( a$p.value, format = "e", digits = 2) ))
  )

```

Figure Legend: Relationship between total peak area at complete features (x-axis) and sample missingness (y-axis).



## Univariate evaluation: batch effects

The figure below provides an illustrative evaluation of the  *total abundance* (at complete features) as a product of sample batch variables provided by your supplier. 

```{r}
if( length(batch_variables) == 0 ) { 
  paste0(" -- No sample level batch variables were provided or all were invariable -- ") 
  }
```

```{r TA_by_batch, echo = FALSE, fig.width = 12, fig.height = 10,  warning=FALSE, message=FALSE, error=FALSE}

if( length(batch_variables) > 0 & length(batch_variables) <= 2 ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$TPA_completefeature, 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "total peak area", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 1)
  
}

```


```{r TA_by_batch_3, echo = FALSE, fig.width = 12, fig.height = 17,  warning=FALSE, message=FALSE, error=FALSE}

if( length(batch_variables) > 2 & length(batch_variables) <= 4 ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$TPA_completefeature, 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "total peak area", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 1)
  
}

```


```{r TA_by_batch_4, echo = FALSE, fig.width = 15, fig.height = 17,  warning=FALSE, message=FALSE, error=FALSE}

if( length(batch_variables) > 4 ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$TPA_completefeature, 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "total abundance", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 2)
  
}

```

Figure Legend: Violin plot illustration(s) of the relationship between total abundance (TA; at complete features) and sample batch variables that are available in your data.


### Multivariate evaluation: batch variables

```{r sample_tpa_multivatiaveANOVA, echo = FALSE, fig.width = 9, fig.height = 2, warnings = FALSE, message=FALSE, error=FALSE}
if( length(batch_variables)>0 ) {
  ( multivariate.anova(dep = qc_data$sample_data$TPA_completefeature, 
                                indep_df = qc_data$sample_data[,batch_variables ] ) ) 
}else {
  paste0(" -- No sample level batch variables were provided or all were invariable -- ")
  }
```

Table Legend: TypeII ANOVA: the eta-squared (eta-sq) estimates are an estimation on the percent of variation explained by each independent variable, after accounting for all other variables, as derived from the sum of squares. This is a multivariate evaluation of batch variables on *total peak|abundance area* at complete features.


# 6. Power analysis

Exploration for case/control and continuous outcome data using the filtered data set

Analytical power analysis for both continuous and imbalanced presence/absence correlation analysis.

```{r power_exploration, echo=F, quote=F, comment=NA, fig.width=15, fig.height=5, warning=FALSE, message=FALSE, error=FALSE }
####################################   
# Run power analysis and generate plot
####################################
p1 = run.cont.power.make.plot( mydata = qc_data$metabolite_data ) 
p2 = run.pa.imbalanced.power.make.plot(  mydata = qc_data$metabolite_data )


ggpubr::ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)

```

Figure Legend: Simulated  effect sizes (standardized by trait SD) are illustrated by their color in each figure. Figure (A) provides estimates of power for continuous traits with the total sample size on the x-axis and the estimated power on the y-axis. Figure (B) provides estimates of power for presence/absence (or binary) traits in an imbalanced design. The estimated power is on the y-axis. The total sample size is set to `r nrow(qc_data$metabolite_data)` and the x-axis depicts the number of individuals present (or absent) for the trait. The effects sizes illustrated here were chosen by running an initial set of simulations which identified effects sizes that would span a broad range of power estimates given the sample population's sample size.

---
title: "metaboprep metabolomics data preparation report"
date: "`r Sys.Date()`"
output: 
  pdf_document:
    keep_md: true
    number_sections: true
    toc: false
    toc_depth: 2
space_betwee_paragraphs: true
fig_caption: true
always_allow_html: yes
link-citations: true
params: 
 Rdatafile: NA
 out_dir: out_dir
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r args_2_variables, include = FALSE, warning=FALSE, message=FALSE, error=FALSE}
library(metaboprep)
library(tidyverse)
suppressPackageStartupMessages(library(dendextend))

###########################
# Read in R object with Rdata 
## file passed as a paramater
###########################
## necessary arguments
## 1) not QC'd metabolite
## 2) sample data
## 3) feature data
## 4) DF'd metabolite data
## 5) data/out directory path
## 6) project
## 7) platform

load(params$Rdatafile)
out_dir = params$out_dir

```


```{r data_2_tibble, include = FALSE, quote=F, comment=NA, warning=FALSE, message=FALSE, error=FALSE }
## RAW DATA
raw_data$metabolite_data = tibble::as_tibble(raw_data$metabolite_data)
raw_data$sample_data = tibble::as_tibble(raw_data$sample_data); colnames(raw_data$sample_data)[1] = "SamID"
raw_data$feature_data = tibble::as_tibble(raw_data$feature_data); colnames(raw_data$feature_data)[1] = "FeatureID"
raw_data$varexp = raw_data$varexp[,1]

## QC Data
qc_data$metabolite_data = tibble::as_tibble(qc_data$metabolite_data)
qc_data$sample_data = tibble::as_tibble(qc_data$sample_data); colnames(qc_data$sample_data)[1] = "SamID"
qc_data$feature_data = tibble::as_tibble(qc_data$feature_data); colnames(qc_data$feature_data)[1] = "FeatureID"
qc_data$varexp = qc_data$varexp[,1]

```

metaboprep report relates to:

  * Project: `r project`  
  * Platform: `r platform`

The `metaboprep` `R` package performs three operations: 

1. Provides an assessment and summary statistics of the raw metabolomics data.
2. Performs data filtering on the metabolomics data.
3. Provides an assessment and summary statistics of the filtered metabolomics data, particularly in the context of batch variables when available.

This report provides descriptive information for raw and filtered metabolomics data for the project `r project`. 

The data filtering workflow is as follows:

```{r QC-pipeline-flowdiagram, echo=FALSE, out.width='65%'}

f = paste0( .libPaths(), "/metaboprep/help/figures/metaboprep_workflow.png" )
knitr::include_graphics( f )

```

1. Issues can be raised on [GitHub](https://github.com/MRCIEU/metaboprep/issues). 
2. Questions relating to the `metaboprep` pipeline can be directed to [David Hughes: d.a.hughes@bristol.ac.uk](mailto:d.a.hughes@bristol.ac.uk).
3. `metaboprep` is published in [Journal to be determined]() and can be cited as:


***


# Sample size of `r project` data set

```{r Sample_size, echo=F, quote=F, comment=NA, fig.width=4, fig.height=1.5, warning=FALSE, message=FALSE, error=FALSE, fig.align="center" }

tout = data.frame("data.set" = c("number of samples","number of features"),  "raw.data" = dim(raw_data$metabolite_data), "filtered.data" = dim(qc_data$metabolite_data))

ggpubr::ggtexttable(tout, rows = NULL, theme = ggpubr::ttheme("mBlue") )

```

***

## Missingness
Missingness is evaluated across samples and features using the original/raw data set. 

### Visual structure of missingness in your raw data set.

```{r MissingnessMatrix, fig.width = 20, fig.height = 12,  echo = FALSE,warning=FALSE, message=FALSE, error=FALSE, dev = "png"}
namat = qc_data$metabolite_data
namat[!is.na(namat)] = 1
namat[is.na(namat)] = 0
namat = as.matrix(namat)

namat = reshape::melt(namat)
colnames(namat) = c("individuals","metabolites","missingness")

pcol = RColorBrewer::brewer.pal(n = 8, name = 'Dark2')

ggplot(namat, aes(metabolites, individuals, fill = missingness)) + 
  geom_tile() + 
  scale_fill_gradient(low= "white", high=pcol[5]) +
  theme(legend.position = "none", axis.text.x = element_blank() )

```

\begin{center}
\textbf{Figure Legend:} Missingness structure across the raw data table. White cells depict missing data. Individuals are in rows, metabolites are in columns. 
\end{center}

### Summary of sample and feature missingness 

```{r run.missingness, echo = FALSE, fig.width = 15, fig.height = 6.5, warning=FALSE, message=FALSE, error=FALSE}
r_mis = missingness.sum(mydata = raw_data$metabolite_data)

p = ggpubr::ggarrange(  r_mis[[4]][[1]] ,
                    r_mis[[4]][[2]] , 
                    r_mis[[4]][[3]] ,
                    r_mis[[4]][[4]] ,
                    ncol = 2, nrow = 2,
                    labels = c("a", "b", "c", "d") )
ggpubr::annotate_figure(p,top = paste0("-- Initial raw data set: Estimates of missingness for samples and features --\n") )

```

\begin{center}
\textbf{Figure Legend:} Raw data - (a) Distribution of sample missingness with sample mean illustrated by the red vertical line. (b) Distribution of feature missingness sample mean illustrated by the red vertical line. (c) Table of sample and feature missingness percentiles. A tabled version of plot a and b.  (d) Estimates of study samples sizes under various levels of feature missingness. 
\end{center}

## Data Filtering 

### Exclusion summary

```{r exclusion_table, echo = FALSE , fig.width = 5, fig.height = 4, warning=FALSE, message=FALSE, error=FALSE}
temp = data.frame(exclusions = rownames(qcing_data$exclusion_data), count = qcing_data$exclusion_data[,1])
ggpubr::ggtexttable( temp, rows = NULL, theme = ggpubr::ttheme("mBlue"))
```

\begin{center}
\textbf{Table Legend}: Six primary data filtering exclusion steps were made during the preparation of the data.
(1) Samples with missingness >=80\% were excluded. (2) features with missingness >=80\% were excluded (xenobiotics are not included in this step). (3) sample exclusions based on the user defined threshold were excluded. (4) feature exclusions based on user defined threshold were excluded (xenobiotics are not included in this step). (5) samples with a total-peak-area or total-sum-abundance that is >= N standard deviations from the mean, where N was defined by the user, were excluded. (6) samples that are >= N standard deviations from the mean on principal component axis 1 and 2, where N was defined by the user, were excluded. 
\end{center}


### Metabolite or feature reduction and principal components

A data reduction was carried out to identify a list of representative features for generating a sample principal component analysis. This step reduces the level of inter-correlation in the data to ensure that the principal components are not driven by groups of correlated features.

```{r PCA_1_2, echo = FALSE, fig.width = 12, fig.height = 6, warning=FALSE, message=FALSE, error=FALSE}
#################
## data reduction 
## info
#################

feature_count = length(qcing_data$feature_sumstats$k)
features_included_in_data_reduction = length( na.omit(qcing_data$feature_sumstats$k) )
cluster_count = length( unique(qcing_data$feature_sumstats$k) )
number_of_rep_meatbolites = sum( qcing_data$feature_sumstats$independent_features_binary == 1 )

temp = data.frame( data.reduction = c("total metabolite count","metabolites included in data reduction","number of metabolite clusters","number of representative metabolites"),
                   count = c(feature_count, features_included_in_data_reduction, cluster_count, number_of_rep_meatbolites  )  )
ptable = ggpubr::ggtexttable( temp, rows = NULL, theme = ggpubr::ttheme("mBlue"))

#################
## PCA
#################
## define tibble
pcs = as_tibble(qcing_data$pcs)
## define color scheme
pcol = RColorBrewer::brewer.pal(9, "Set1")
## define accelerationfactor
accelerationfactor = qcing_data$accelerationfactor
if(accelerationfactor > 10){accelerationfactor = 10}
if(accelerationfactor == 1){accelerationfactor = 2}
## define outliers
Omat = outlier.matrix(pcs[, 1:accelerationfactor], nsd = PC_outlier_SD)
outliers = apply(Omat, 1, sum)
outliers[outliers>0]=1
pcs$outliers = as.factor( outliers )
##
cutoffs = apply(pcs[, 1:accelerationfactor], 2, function(x){
  msd = c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  cutoff = msd[1] + (msd[2]*PC_outlier_SD)
  return(cutoff)
} )
  
##
pcplot = pcs %>% ggplot( aes(x = PC1, y = PC2) ) +
  geom_point( aes(fill = outliers), size = 2.5, shape = 21 ) +
  scale_fill_manual(values = pcol[c(2,1)]) +
  geom_vline(xintercept = c(cutoffs[1], -cutoffs[1]), color = pcol[1] ) +
  geom_hline(yintercept = c(cutoffs[2], -cutoffs[2]), color = pcol[1] ) +
  theme(legend.position="bottom") +
  labs(title = paste0("Principal components 1-&-2 using ",number_of_rep_meatbolites ," representative metabolites"),
       subtitle  = paste0(" - Outliers are those ", PC_outlier_SD , " SD from the mean of PCs 1 to", accelerationfactor ))
##
gridExtra::grid.arrange( grobs = list(ptable, pcplot), widths = c(2, 3), ncol = 2) 

```

\begin{center}
\textbf{Figure Legend:} The data reduction table on the left presents the number of metabolites at each phase of the data reduction (Spearman's correlation distance tree cutting) analysis. On the right principal components 1 and 2 are plotted for all individuals, using the representitve features identified in the data reduction analysis. The red vertical and horizontal lines indicate the standard deviation (SD) cutoffs for identifying individual outliers, which are plotted in red. The standard deviations cuttoff were defined by the user. 
\end{center}

***


# Filtered data

## N 
  * The number of samples in data = `r nrow(qc_data$metabolite_data)`  
  * The number of features in data = `r ncol(qc_data$metabolite_data)`  

## Relative to the raw data
  * `r nrow(raw_data$metabolite_data) - nrow(qc_data$metabolite_data)` samples were filtered out, given the user's criteria.
  * `r ncol(raw_data$metabolite_data) - ncol(qc_data$metabolite_data)` features were filtered out, given the user's criteria.
  * Please review details above and your log file for the number of features and samples excluded and why. 
  
## Summary of filtered data


```{r sample_overview, echo = FALSE, fig.width = 15, fig.height = 15, warning=FALSE, message=FALSE, error=FALSE}
pcol = RColorBrewer::brewer.pal(9, "Set1")
##############
## Missingness
##############
qc_mis = missingness.sum(mydata = qc_data$metabolite_data)
s_mis_plot = qc_mis[[4]][[1]] 
f_mis_plot = qc_mis[[4]][[2]] 

##############
## TPA
##############
s = tibble::tibble(TPA_completefeature = qc_data$sample_data$TPA_completefeature)

tpa_plot = qc_data$sample_data %>% ggplot( aes( TPA_completefeature ) ) +
    geom_histogram( fill = pcol[2] , bins = 25) + 
    geom_vline( xintercept = median(qc_data$sample_data$TPA_completefeature), color = pcol[1], size = 1) +
    labs(title = paste0("total abundance of samples\nat complete features only"), x="", y="") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))

##############
## Dendrogram
##############
dend = qc_data$feature_tree %>% as.dendrogram

## create a vector of colors to color your tree labels
w = which( qc_data$feature_data$independent_features_binary == 1)
pv_ids = as.character( unlist( qc_data$feature_data[w,1] ) )
n = labels(dend)
bcol = rep("black", length(n))
w = which(n %in% pv_ids ); bcol[w] = pcol[2]

## redefine elements of dendrogram
dend = dend %>% 
  set("labels_cex", 0.5) %>% 
  set("labels_col", bcol) %>% 
  set("branches_lwd", 0.5) %>%
  set("branches_k_color",  value = bcol)
dend <- as.ggdend(dend)
## plot the dendrogram
tree_plot = dend %>% ggplot() + geom_hline(yintercept = tree_cut_height, color = "coral2")

##############
## Data Reduce
## Table
##############
feature_count = length(qc_data$feature_data$k)
features_included_in_data_reduction = length( na.omit(qc_data$feature_data$k) )
cluster_count = length( unique(qc_data$feature_data$k) )
number_of_rep_meatbolites = sum( qc_data$feature_data$independent_features_binary == 1 )

temp = data.frame( data.reduction = c("total metabolite count","metabolites included in data reduction","number of metabolite clusters","number of representative metabolites"),
                   count = c(feature_count,features_included_in_data_reduction, cluster_count, number_of_rep_meatbolites  )  )
ptable = ggpubr::ggtexttable( temp, rows = NULL, theme = ggpubr::ttheme("mGreen"))

#############
## scree plot
#############
ve = tibble::tibble( data.frame(PC= 1:length(qc_data$varexp), VarExp = qc_data$varexp) )

screeplot = ve %>% ggplot(aes(x = PC, y = VarExp)) +
  geom_line(color = "grey") +
  geom_point(shape = 21, fill = pcol[2], size = 2) +
  labs(title = "Scree plot") +
  geom_vline(xintercept = qc_data$accelerationfactor, color = pcol[1]) +
  geom_vline(xintercept = unlist( qc_data$nparallel ), color = pcol[3])

##############
## PC Plot
##############
pcs = as_tibble(qc_data$sample_data[, c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
pcs$cluster_k = as.factor( kmeans(pcs[, 1:2], 4)$cluster )
## define color scheme
pcol = RColorBrewer::brewer.pal(9, "Set1")
##
pcplot = pcs %>% ggplot( aes(x = PC1, y = PC2) ) +
  geom_point( size = 2.5, shape = 21, aes(fill = cluster_k) ) +
  scale_fill_manual(values = pcol[1:4]) +
  labs( title = paste0("Principal components 1-&-2 using ", number_of_rep_meatbolites ," representative metabolites"),
        x = paste0("PC1  VarExp = ", signif(qc_data$varexp[1], d = 4)*100, "%" ),
        y = paste0("PC2  VarExp = ", signif(qc_data$varexp[2], d = 4)*100 , "%"),
        fill = paste0("kmeans\ncluster k"))
#######################
## Plot the data
#######################
m = matrix(c(1,2,3,4,4,4,5,7,7, 6,7,7), nrow = 4, byrow = TRUE)

gridExtra::grid.arrange( grobs = list(s_mis_plot, f_mis_plot, tpa_plot, tree_plot, ptable, screeplot, pcplot), 
                         layout_matrix = m,
                         heights = c(1,1,1,1) ) 

```


\begin{center}
\textbf{Figure Legend:} Filtered data summary. Distributions for sample missingness, feature missingness, and total abundance of samples. Row two of the figure provides a Spearman's correlation distance clustering dendrogram highlighting the metabolites used as representative features in blue, the clustering tree cut height is denoted by the horizontal line. Row three provides a summary of the metabolite data reduction in the table, a Scree plot of the variance explained by each PC and a plot of principal component 1 and 2, as derived from the representative metabolites. The Scree plot also identifies the number of PCs estimated to be informative (vertical lines) by the Cattel's Scree Test acceleration factor (red,  n = `r qc_data$accelerationfactor`) and Parallel Analysis (green, n = `r unlist( qc_data$nparallel )`). Individuals in the PC plot were clustered into 4 kmeans (k) clusters, using data from PC1 and PC2. The kmeans clustering and color coding is strictly there to help provide some visualization of the major axes of variation in the sample population(s).
\end{center}


## Structure among samples: top 5 PCs

```{r PCApairsplot, echo = FALSE, fig.width = 15, fig.height = 15, warning=FALSE, message=FALSE, error=FALSE}
plotcolors = pcol[pcs$cluster_k]
pcapairs_bymoose( as.matrix(pcs[, 1:5]) , qc_data$varexp, pcol = plotcolors)
```

\begin{center}
\textbf{Figure Legend:} A matrix plot of the top five principal components including demarcations of the 3rd (grey), 4th (orange), and 5th (red) standard deviations from the mean. Samples are color coded as in the summary PC plot above using a kmeans analysis of PC1 and PC2 with a k (number of clusters) set at 4. The choice of k = 4 was not robustly chosen it was a choice of simplicity to help aid visualize variation and sample mobility across the PCs. 
\end{center}

## Feature Distributions

### Estimates of normality: W-statistics for raw and log transformed data

```{r, echo=F, warning=FALSE, message=FALSE, error=FALSE}
wstat = qc_data$feature_data$W_stat_rawdata

## how many stats were estimated
count = length(wstat)
nacount = sum(is.na(wstat))
remain_count = count - nacount
# wstat = na.omit(wstat)
normcount = sum(wstat >= 0.95, na.rm = TRUE)

```


```{r shapiroW, echo = FALSE, fig.width = 10, fig.height = 5, warning=FALSE, message=FALSE, error=FALSE}
W = wstat
pcol = RColorBrewer::brewer.pal(9, "Set1")

W_log10 = qc_data$feature_data$W_stat_log10data

## Plot
par(mfrow = c(1,2), oma = c(2,1,1,1))
hist(W, col = pcol[2], 
     main = paste("Distribution of W statistics on\nRaw Metabolite Abundances"), 
     xlab = "W", 
     ylab = "Frequency",
     cex.main = 0.75)
abline(v = 0.95, col = pcol[1])
##
mtext( paste0(normcount,
              " of the metabolites exhibit distributions\nthat may declared normal, given a W-stat >= 0.95"), 
       cex = 0.75, side = 1, outer = TRUE, line = 0, adj = 0.1, col = pcol[2])
  ### LOG
LogMakesDistributionWorse = c( sum( W_log10 < W , na.rm = TRUE), signif( sum( W_log10 < W , na.rm = TRUE)/length(!is.na(W) ), d = 3)*100)
  ##
hist(W_log10, col = pcol[3], 
     main = paste("Distribution of W statistics on\nlog10 Metabolite Abundances"), 
     xlab = "W", 
     ylab = "Frequency",
     cex.main = 0.75)
mtext( paste0("In ", LogMakesDistributionWorse[1],
              " instances or ", LogMakesDistributionWorse[2], "% of the tested metabolites\nthe log10 data W-stat is < raw data W-stat."), cex = 0.75, side = 1, outer = TRUE, line = 0, adj = 0.9, col = pcol[3])

```

\begin{center}
\textbf{Figure Legend:} Histogram plots of Shapiro W-statistics for raw (left) and log transformed (right) data distributions. A W-statistic value of 1 indicates the sample distribution is perfectly normal and value of 0 indicates it is perfectly uniform. Please note that log transformation of the data *may not* improve the normality of your data.
\end{center}

\begin{center}
\textbf{Analysis details:} Of the `r count` features in the data `r nacount` features were excluded from this analysis because of no variation or too few observations (n < 40). Of the remaining `r remain_count` metabolite features, a total of `r normcount` may be considered normally distributed given a Shapiro W-statistic >= 0.95.
\end{center} 

### Distributions
A pdf report is being written to `r paste0( project, "_outlier_detection_pre_filtering.pdf")` that contains dotplot, histogram and distribution summary statistics for each metabolite in your data set, providing an opportunity to visually inspect all your metabolites feature data distributions.

## Outliers
**Evaluation of the number of samples and features that are outliers across the data.**

```{r outliers, echo = FALSE, fig.width = 6, fig.height = 2.5, warning=FALSE, message=FALSE, error=FALSE}
today = Sys.Date()
today = gsub("-","_",today)

outlier_filename <- paste0( out_dir, project, "_outlier_detection_pre_filtering.pdf")

## feauture outlier summary
outlier_summary = outlier.summary(dtst = qc_data$metabolite_data, 
                           pdf_filename = outlier_filename, 
                           nsd = 5)
```

```{r outlier_summary, echo = FALSE, fig.width = 8, fig.height = 2.5, warning=FALSE, message=FALSE, error=FALSE}
outlier_summary[[2]]
```
\begin{center}
\textbf{Table Legend:} The table reports the number of point estimates for the minimum (0\%) median (50\%) and maximum (100\%) number of outlying features across samples and the number of outlying samples across features. 
\end{center}

### Notes on outlying samples at each metabolite|feature

There may be extreme outlying observations at individual metabolites|features that have not been accounted for. You may want to:

1. Turn these observations into NAs.
2. Winsorize the data to some maximum value.
3. Rank normalize the data which will place those outliers into the top of the ranked standard normal distribution.
4. Turn these observations into NAs and then impute them along with other missing data in your data set. 


***

# Influence of batch variables on filtered data

## Filtered data *feature* missingness: influenced by possible explanatory variables 

Feature missingness may be influenced by the metabolites (or features) biology or pathway classification, or your technologies methodology. The figure(s) below provides an illustrative evaluation of the proportion of *feature missigness* as a product of the variable(s) available in the raw data files. 

```{r featurebatchvariables, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE}
possible_batch_variables = c("path","platform", "class", "batch")

## matching variables
w = unique( unlist( sapply(possible_batch_variables, function(x){ grep(x, tolower( colnames(qc_data$feature_data)  )  ) }) ) )
if(length(w)>0){
  matched_variables = colnames(qc_data$feature_data)[w]
  
  ## remove batch variables, that are invariable in the sample
  w = which( apply( qc_data$feature_data[,matched_variables], 2, function(x){ length(table(x)) }) <= 1 )
  if(length(w)>0){
    matched_variables = matched_variables[-w]  
  }
  ## insure group sizes are larger than 1
  w = which( apply( qc_data$feature_data[,matched_variables], 2, function(x){ sum( table(x) > 1 ) / length(table(x)) }) < 0.25 )
  if(length(w)>0){
    matched_variables = matched_variables[-w]  
  }
  ## insure the number of variables is not larger than 96
  w = which( apply( qc_data$feature_data[,matched_variables], 2, function(x){  length(table(x)) }) > 96 )
  if(length(w)>0){
    matched_variables = matched_variables[-w]  
  }
  
  ## define the class_variables
  class_variables = matched_variables
  
  if( "SUB_PATHWAY" %in% class_variables ){
    w = which(class_variables %in% "SUB_PATHWAY")
    class_variables = class_variables[-w]
  }

}


```


```{r featuremissingness,  echo=F, fig.width = 15, fig.height = 10, warning=FALSE, message=FALSE, error=FALSE}
## MISSINGNESS
###################################
## Iterate over each batch variable
###################################
if(length(class_variables)>0){
  ClassMisPlots = lapply( class_variables , function(x){
  out = variable.by.factor( dep = qc_data$feature_data$feature_missingness , 
                           indep = unlist( qc_data$feature_data[,x] ), 
                           dep_name = "feature missingness", 
                           indep_name = x, orderfactor = TRUE, violin = FALSE )
  return(out)
  })


  ## plot the output
  gridExtra::grid.arrange( grobs = ClassMisPlots, ncol = 1)

} else { 
  paste0(" -- No feature level batch variables identified or all were invariable -- ") 
  }

```

\begin{center}
\textbf{Figure Legend:} Box plot illustration(s) of the relationship that available batch and biological variables have with feature missingness.
\end{center}


## Filtered data *sample* missingness: influenced by possible explanatory variables 

The figure provides an illustrative evaluation of the proportion of *sample missigness* as a product of sample batch variables provided by your supplier. This is the univariate influence of batch effects on *sample missingness*.
   
```{r samplebatchvariables, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE}
possible_batch_variables = c("box", "day", "prep", "date", "time", "type", "nmr", "buffer", "lmwm")
  
## matching variables
# w = which( possible_batch_variables  %in%  colnames(qc_data$sample_data)  )
w = unique( unlist( sapply(possible_batch_variables, function(x){ grep(x, tolower( colnames(qc_data$sample_data)  )  ) }) ) )

if( length(w) > 0 ){
  # matched_variables = possible_batch_variables[w]
  matched_variables = colnames(qc_data$sample_data)[w]
  
  ## batch variables, that are variable in the sample
  w = which( apply( qc_data$sample_data[,matched_variables], 2, function(x){ length( table(x) ) } ) > 1 )
  matched_variables = matched_variables[w]
  ## insure group sizes are larger than 1
  w = which( apply( qc_data$sample_data[,matched_variables], 2, function(x){ sum( table(x) > 1 ) / length(table(x)) }) >= 0.25 )
  matched_variables = matched_variables[w]
  ## insure the number of variables is not larger than 96
  w = which( apply( qc_data$sample_data[,matched_variables], 2, function(x){  length(table(x)) }) < 97 )
  matched_variables = matched_variables[w]
  
  batch_variables = matched_variables
} else {
  batch_variables = NA
}

## In case all batch variables are not variable (typically because sample size is small and all done in 1 go.)
if(length(batch_variables) == 0){ batch_variables = NA }

```

```{r sample_missingness,  echo=F, fig.width = 15, fig.height = 15, warning=FALSE, message=FALSE, error=FALSE}
## MISSINGNESS
###################################
## Iterate over each batch variable
###################################
if( !is.na(batch_variables[1]) ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$sample_missingness , 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "sample missingness", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 1)
  
} else {
  paste0(" -- No sample level batch variables were provided or all were invariable -- ")
}

```

\begin{center}
\textbf{Figure Legend:} Box plot illustration(s) of the relationship that available batch variables have with sample missingness.
\end{center}

## Multivariate evaluation: batch variables

```{r sample_missingness_multivatiaveANOVA, echo = FALSE, fig.width = 9, fig.height = 2, warning=FALSE, message=FALSE, error=FALSE}
if( !is.na(batch_variables[1]) ) {
  ( multivariate.anova(dep = qc_data$sample_data$sample_missingness, indep_df = qc_data$sample_data[ ,batch_variables ] ) )
} else {
  paste0(" -- No sample level batch variables were provided or all were invariable -- ")
  }
```
\begin{center}
\textbf{Table Legend:} TypeII ANOVA: the eta-squared (eta-sq) estimates are an estimation of the percent of variation explained by each independent variable, after accounting for all other variables, as derived from the sum of squares. This is a multivariate evaluation of batch variables on *sample missingness*.
\end{center}

# Sample Total Peak|Abundance Area (TPA):
Total peak|abundance area (TPA) is simply the sum of the abundances measured across all features. TPA is one measure that can be used to identify unusual samples given their entire profile. However, the level of missingness in a sample may influence TPA. To account for this we:  

1. Evaluate the correlation between TPA estimates across all features with TPA measured using only those features with complete data (no missingness).  
2. Determine if the batch effects have a measurable impact on TPA.

## Relationship with missingness
Correlation between total peak area (at complete features) and missingness

```{r TPA_missingness, echo=F, fig.width = 7, fig.height = 5, warning=FALSE, message=FALSE, error=FALSE}
a = cor.test( qc_data$sample_data$sample_missingness, qc_data$sample_data$TPA_completefeature)
###
( 
  tpamis = qc_data$sample_data %>% ggplot( aes(x = TPA_completefeature, y = sample_missingness)) +
  geom_point( fill = "grey", alpha = 0.8, size = 1.5) + 
  geom_smooth(method = "loess", color = "red", size = 2)  +
  geom_smooth(method = "lm", color = "blue", size = 2)  +
  labs(x = "TPA, complete features", y = "sample missingness",
       title = paste0( "TPA as influenced by missingness \nSpearmans's cor = ", round(a$estimate, d = 4), " p-value = ", 
                    formatC( a$p.value, format = "e", digits = 2) ))
  )

```

\begin{center}
\textbf{Figure Legend:} Relationship between total peak area at complete features (x-axis) and sample missingness (y-axis).
\end{center}



## Univariate evaluation: batch effects
The figure below provides an illustrative evaluation of the  *total peak area* as a product of sample batch variables provided by your supplier. 

```{r TPAsummary, echo = FALSE, fig.width = 12, fig.height = 10,  warning=FALSE, message=FALSE, error=FALSE}

if( !is.na(batch_variables[1]) ) {
  BatchMisPlots = lapply( batch_variables , function(x){
    out = variable.by.factor( dep = qc_data$sample_data$TPA_completefeature, 
                              indep = unlist( qc_data$sample_data[,x] ), 
                              dep_name = "total peak area", 
                              orderfactor = FALSE, 
                              indep_name = x, violin = TRUE)
    return(out)
  })
  
  ## plot the output
  gridExtra::grid.arrange( grobs = BatchMisPlots, ncol = 1)
  
} else {
  paste0(" -- No sample level batch variables were provided or all were invariable -- ")
}

```
\begin{center}
\textbf{Figure Legend:} Violin plot illustration(s) of the relationship between total peak area (TPA) and sample batch variables that are available in your data.
\end{center}


### Multivariate evaluation: batch variables

```{r sample_tpa_multivatiaveANOVA, echo = FALSE, fig.width = 9, fig.height = 2, warnings = FALSE, message=FALSE, error=FALSE}
if( !is.na(batch_variables[1]) ) {
  ( multivariate.anova(dep = qc_data$sample_data$TPA_completefeature, 
                                indep_df = qc_data$sample_data[,batch_variables ] ) ) 
}else {
  paste0(" -- No sample level batch variables were provided or all were invariable -- ")
  }
```
\begin{center}
\textbf{Table Legend:} TypeII ANOVA: the eta-squared (eta-sq) estimates are an estimation on the percent of variation explained by each independent variable, after accounting for all other variables, as derived from the sum of squares. This is a multivariate evaluation of batch variables on *total peak|abundance area* at complete features.
\end{center}


***

# Power analysis

**Exploration for case/control and continuous outcome data using the filtered data set**

Analytical power analysis for both continuous and imbalanced presence/absence correlation analysis.

```{r power_exploration, echo=F, quote=F, comment=NA, fig.width=15, fig.height=5, warning=FALSE, message=FALSE, error=FALSE }
####################################   
# Run power analysis and generate plot
####################################
#( run.power.make.plot(mydata = metdata ) )

p1 = run.cont.power.make.plot( mydata = qc_data$metabolite_data ) 
p2 = run.pa.imbalanced.power.make.plot(  mydata = qc_data$metabolite_data )

# gridExtra::grid.arrange(p1, p2, nrow = 1)

ggpubr::ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)

```
\begin{center}
\textbf{Figure Legend:} Simulated effect sizes are illustrated by their color in each figure. Figure (A) provides estimates of power for continuous traits with the total sample size on the x-axis and the estimated power on the y-axis. Figure (B) provides estimates of power for presence/absence (or binary) traits in an imbalanced design. The estimated power is on the y-axix. The total sample size is set to `r nrow(qc_data$metabolite_data)` and the x-axis depicts the number of individuals present (or absent) for the trait. The effects sizes illustrated here were chosen by running an initial set of simulations which identified effects sizes that would span a broad range of power estimates given the sample population's sample size. 
\end{center}


***



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.missingness.R
\name{sample.missingness}
\alias{sample.missingness}
\title{estimate sample missingness}
\usage{
sample.missingness(wdata, excludethesefeatures = NA)
}
\arguments{
\item{wdata}{a numeric matrix with samples in row and features in columns}

\item{excludethesefeatures}{a vector of feature names (i.e. column names) to exclude from missingness estimates}
}
\value{
A data frame of missingness estimates for each sample. If a vector of feature names was also passed to the function a second column of missingness estimates will also be returned providing missingness estimates for each sample to the exclusion of those features provided.
}
\description{
This function estimates sample missingness in a matrix of data and provides an option to exclude certain columns or features from the analysis, such as xenobiotics (with high missingness rates) in metabolomics data sets.
}
\examples{
## simulate some data
set.seed(1110)
ex_data = sapply(1:5, function(x){ rnorm(10, 40, 5) })
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add some missingness to the data
ex_data[ sample(1:50, 10) ] = NA
## estimate missingness
mis_est = sample.missingness(ex_data)
mis_est_v2 = sample.missingness(ex_data, excludethesefeatures = "var5")

}
\keyword{missingness}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pc.and.outliers.R
\name{pc.and.outliers}
\alias{pc.and.outliers}
\title{principal component analysis}
\usage{
pc.and.outliers(metabolitedata, indfeature_names, outliers = TRUE)
}
\arguments{
\item{metabolitedata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{indfeature_names}{a vector of independent feature names | column names.}

\item{outliers}{defaulted to TRUE, a TRUE|FALSE binary flagging if you would like outliers identified.}
}
\value{
a list object of length five, with (1) a data frame of PC loadings, (2) a vector of variance explained estimates for each PC, (3) an estimate of the number of informative or top PCs determined by the acceleration factor analysis, (4) an estimate of the number of informative or top PCs determined by parrallel analysis, (5) a data frame of the probablisitic PC loadings
}
\description{
This function performs two principal component analysis. In the first, missing data is imputed to the median. In the second a probablistic PCA is run to account for the missingness. 
Subsequent to the derivation of the PC, the median imputed PC data is used to identify the number of informative or "significant" PC by (1) an acceleration analysis, and (2) a parrallel analysis.
Finally the number of sample outliers are determined at 3, 4, and 5 standard deviations from the mean on the top PCs as determined by the acceleration factor analysis.
}
\examples{
## define a covariance matrix
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.8, 0.6, 0.2)
cmat[2,] = c(0.8, 1, 0.7, 0.5)
cmat[3,] = c(0.6, 0.7, 1, 0.6)
cmat[4,] = c(0.2, 0.5, 0.6,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
d1 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
set.seed(1010)
d2 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## simulate some random data
d3 = sapply(1:20, function(x){ rnorm(250, 40, 5) })
## define the data set
ex_data = cbind(d1,d2,d3)
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add in some missingness
ex_data[sample(1:length(ex_data), 450)] = NA
## add in some technical error to two samples
m = apply(ex_data, 2, function(x){ mean(x, na.rm = TRUE) })
ex_data[c(1,10), ] = ex_data[1, ] + (m*0.00001) 
## run the PCA
ex_out = pc.and.outliers(ex_data, indfeature_names = sample(colnames(ex_data), 15) )


}
\keyword{PCA}
\keyword{probalistic}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval.power.binary.imbalanced.R
\name{eval.power.binary.imbalanced}
\alias{eval.power.binary.imbalanced}
\title{Estimate power for a binary variable in an imbalanced design}
\usage{
eval.power.binary.imbalanced(N_case, N_control, effect, alpha)
}
\arguments{
\item{N_case}{a numeric vector of sample size of cases}

\item{N_control}{a numeric vector of sample size of controls}

\item{effect}{a numeric vector of effect size}

\item{alpha}{a numeric vector of significance thresholds}
}
\value{
a matrix of paramater inputs and power estimates are returned as a matrix
}
\description{
This function allows you estimate power for a binary variable given a defined number of case samples, control samples, effect size, and significance threshold.
}
\examples{
eval.power.binary.imbalanced( N_case = 1000, 
 N_control = 1000, 
 effect = 0.01, 
 alpha = 0.05 )

eval.power.binary.imbalanced( N_case = c(1000, 2000), 
 N_control = c(1000, 2000), 
 effect = 0.01, 
 alpha = 0.05 )


}
\keyword{metabolomics}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.tree.independence.R
\name{feature.tree.independence}
\alias{feature.tree.independence}
\title{identify independent features}
\usage{
feature.tree.independence(wdata)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}
}
\value{
a data frame of 'k' cluster or group ids, and a 0/1 binary identifying if a feature was identified as and independent ('1') feature or not ('0').
}
\description{
This function identifies independent features using Spearman's Rho, and a dendrogram tree cut step. The feature returned as 'independent' within is k-cluster is the feature with the least missingness or chosen at random in case of missingness ties.
}
\examples{
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.9, 0.9, 0.8)
cmat[2,] = c(0.9, 1, 0.7, 0.6)
cmat[3,] = c(0.9, 0.7, 1, 0.8)
cmat[4,] = c(0.8, 0.6, 0.8,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
ex_data = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## define the data set
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## run the function
df = feature.tree.independence(ex_data)

}
\keyword{feature}
\keyword{independece}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find.PA.effect.sizes.2.sim.R
\name{find.PA.effect.sizes.2.sim}
\alias{find.PA.effect.sizes.2.sim}
\title{identify effect sizes}
\usage{
find.PA.effect.sizes.2.sim(mydata)
}
\arguments{
\item{mydata}{Your metabolite data matrix, with samples in rows}
}
\value{
a vector of effect sizes
}
\description{
This function estimates an appropriate distribution of effect sizes to simulate in a power analysis.
}
\examples{
ex_data = sapply(1:10, function(x){ rnorm(250, 40, 5) })
find.PA.effect.sizes.2.sim(ex_data)

}
\keyword{anlysis}
\keyword{power}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.sum.stats.R
\name{sample.sum.stats}
\alias{sample.sum.stats}
\title{summary statistics for samples}
\usage{
sample.sum.stats(wdata, feature_names_2_exclude = NA, outlier_udist = 5)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{feature_names_2_exclude}{a vector of feature|column names to exclude from missingness estimates}

\item{outlier_udist}{the interquartile range unit distance from the median to call a sample an outlier at a feature.}
}
\value{
a data frame of summary statistics
}
\description{
This function estimates summary statistics for samples in a matrix of numeric features. This includes missingness, total peak area, and a count of the number of outlying features for a sample.
}
\examples{
## simulate some data
set.seed(1110)
ex_data = sapply(1:5, function(x){ rnorm(10, 40, 5) })
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add some missingness to the data
ex_data[ sample(1:50, 10) ] = NA
## run estimate sample summary statistics
sample.sum.stats(ex_data)

}
\keyword{sample}
\keyword{statistics}
\keyword{summary}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.outliers.R
\name{feature.outliers}
\alias{feature.outliers}
\title{outlier sample count for a features}
\usage{
feature.outliers(wdata, nsd = 5)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{nsd}{the number of standard deviation from the mean outliers are identified at. Default value is 5.}
}
\value{
a data frame out sample outlier counts for each feature (column) in the matrix
}
\description{
This function takes a matrix of data (samples in rows, features in columns) and counts the number of outlying samples each feature has.
}
\examples{
ex_data = sapply(1:20, function(x){ rnorm(250, 40, 5) })
s = sample(1:length(ex_data), 200)
ex_data[s] = ex_data[s] + 40
## run the function
fout = feature.outliers(ex_data)

}
\keyword{feature}
\keyword{outliers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ng_anno}
\alias{ng_anno}
\title{Nightingale Health metabolomics annotation data set}
\format{
A data frame with 233 rows and 7 variables:
\describe{
  \item{metabolite}{metabolite id}
  \item{raw.label}{metabolite name and units}
  \item{class}{metabolite annotation class}
  \item{subclass}{metabolite annotation subclass}
  \item{label}{metabolite name and units}
  \item{label.no.units}{metabolite name without units}
  \item{derived_features}{a binary yes|no indicating if the metabolite variable is a variable derived of two or more other features in the data set}
  ...
}
}
\source{
\url{http://nightingalehealth.com/}
}
\usage{
ng_anno
}
\description{
A dataset containing annotation data on Nightingale Health NMR metabolites (lipids), compiled from public resources with additional annotation made here.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.cont.power.make.plot.R
\name{run.cont.power.make.plot}
\alias{run.cont.power.make.plot}
\title{continuous trait power analysis plot}
\usage{
run.cont.power.make.plot(mydata)
}
\arguments{
\item{mydata}{Your metabolite data matrix, with samples in rows}
}
\value{
a ggplot2 object
}
\description{
This function (1) identifies an informative distribution of effect and power estimates given your datas total sample size and (2) returns a summary plot.
}
\examples{
ex_data = matrix(NA, 1000, 2)
run.cont.power.make.plot( ex_data )

}
\keyword{analysis}
\keyword{continuous}
\keyword{plot}
\keyword{power}
\keyword{trait}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/greedy.pairwise.n.filter.R
\name{greedy.pairwise.n.filter}
\alias{greedy.pairwise.n.filter}
\title{greedy selection}
\usage{
greedy.pairwise.n.filter(wdata, minN = 50)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{minN}{the minimum sample size (n) for pairwise comparisons}
}
\value{
a vector of feature names
}
\description{
This function identifies features that have less than a minimum number of complete pairwise observations and removes one of them, in a greedy fashion.
The need for this function is in instances where missingness is extreme between two features the number of paired observation between them may be to
to be informative. Thus one, but not both should be removed from the analysis to avoid analytical error based on sample sizes.
}
\examples{
set.seed(123)
ex_data = sapply(1:10, function(x){ rnorm(250, 40, 5) })
## define the data set
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add in some missingness
ex_data[ sample(1:250, 200) ,1] = NA
ex_data[ sample(1:250, 200) ,2] = NA
ex_data[ sample(1:250, 200) ,3] = NA
## Estimate missingness and generate plots
greedy.pairwise.n.filter(ex_data)


}
\keyword{feature}
\keyword{selection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.in.metabolon.R
\name{read.in.metabolon}
\alias{read.in.metabolon}
\title{read in Metabolon (v1) metabolomics data}
\usage{
read.in.metabolon(file2process, data_dir, projectname)
}
\arguments{
\item{file2process}{the name of the xls file to process}

\item{data_dir}{the full path to the directory holding your Metabolon excel file}

\item{projectname}{a name for your project}
}
\value{
a list object of (1) metabolite, (2) sample annotation, and (3) feature annotation data
}
\description{
This function reads in a Metabolon (v1 format) raw data excel file, writes the (1) metabolite, (2) sample annotation, and (3) feature annotation data to flat text files. It also returns a list object of the same data.
}
\examples{
# read.in.metabolon(file2process = "Metabolon_data_release.xls", 
#  data_dir = "/File/sits/here/", 
#  projectname = "My Amazing Project")

}
\keyword{Metabolon}
\keyword{metabolomics}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missingness.sum.R
\name{missingness.sum}
\alias{missingness.sum}
\title{missingness summary plots}
\usage{
missingness.sum(mydata)
}
\arguments{
\item{mydata}{metabolite data matrix, with samples in rows and metabolite features in columns.}
}
\value{
a list object of length four with (1) a vector of sample missingess,(2) a vector of feature missingness ,(3) a table summarizing missingness ,(4) a list of ggplot2 plots for sample and feature histogram distribution and summary tables
}
\description{
This function sumarizes sample and feature missingness in tables and in a summary plot.
}
\examples{
ex_data = sapply(1:100, function(x){ rnorm(250, 40, 5) })
## define the data set
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add in some missingness
ex_data[sample(1:length(ex_data), 500)] = NA
## Estimate missingness and generate plots
ms = missingness.sum(ex_data)
## plots
ggpubr::ggarrange(plotlist = ms$plotsout, ncol = 2, nrow = 2)


}
\keyword{metabolomics}
\keyword{missingness}
\keyword{plots}
\keyword{summary}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find.cont.effect.sizes.2.sim.R
\name{find.cont.effect.sizes.2.sim}
\alias{find.cont.effect.sizes.2.sim}
\title{identify continuos trait effect sizes}
\usage{
find.cont.effect.sizes.2.sim(mydata)
}
\arguments{
\item{mydata}{Your metabolite data matrix, with samples in rows}
}
\value{
a vector of effect sizes
}
\description{
This function estimates an appropriate distribution of effect sizes to simulate in a continuous trait power analysis.
}
\examples{
ex_data = sapply(1:10, function(x){ rnorm(250, 40, 5) })
find.cont.effect.sizes.2.sim(ex_data)

}
\keyword{analysis}
\keyword{power}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_plots.R
\name{feature_plots}
\alias{feature_plots}
\title{feature plots to file}
\usage{
feature_plots(wdata, outdir = NULL, nsd = 5)
}
\arguments{
\item{wdata}{a data frame of feature (ex: metabolite or protein) abundance levels}

\item{outdir}{output directory path}

\item{nsd}{number of SD from the mean to plot an outlier line on the scatter plot and histogram}
}
\value{
a ggplot2 object
}
\description{
This function to plots a scatter plot, a histogram, and a few summary statistics of each feature in a data frame to a pdf file
}
\examples{
ex_data = sapply(1:20, function(x){ rnorm(250, 40, 5) })
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
feature_plots(ex_data)

}
\keyword{metabolomics}
\keyword{pdf}
\keyword{summary}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform.metabolite.qc.R
\name{perform.metabolite.qc}
\alias{perform.metabolite.qc}
\title{perform metabolomics quality control}
\usage{
perform.metabolite.qc(
  wdata,
  fmis = 0.2,
  smis = 0.2,
  tpa_out_SD = 5,
  outlier_udist = 5,
  outlier_treatment = "leave_be",
  winsorize_quantile = 1,
  tree_cut_height = 0.5,
  PC_out_SD = 5,
  feature_colnames_2_exclude = NA,
  derived_colnames_2_exclude = NA
)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{fmis}{defaulted at 0.2, this defines the feature missingness cutoff}

\item{smis}{defaulted at 0.2, this defines the sample missingness cutoff}

\item{tpa_out_SD}{defaulted at 5, this defines the number of standard deviation from the mean in which samples will be excluded for total peak area. Pass NA to this paramater to exclude this step.}

\item{outlier_udist}{defaulted at 5, this defines the number of interquartile range units from the median to define a value as an outlier}

\item{outlier_treatment}{defaulted to "leave_be". Choices are "leave_be", "winsorize", or "turn_NA", which defines how to treat outlier values prior to estimating principal components. "leave_be" will do nothing to ourlier values. "turn_NA" will turn outliers into NA and thus be median imputed for the purpose of the PCA. "winsorize" will turn NA values into the "winsorize_quantile" of all remaining (non-outlying) values at a feature.}

\item{winsorize_quantile}{the quantile (0-1) to winzorise outlier values to, if "outlier_treatment" parameter set to "winsorize". Defaulted to 1, or the maximum value of all remaining (non-outlying) values at a feature.}

\item{tree_cut_height}{The height at which to cut the feature|metabolite dendrogram to identify "independent" features. tree_cut_height is 1-absolute(Spearman's Rho) for intra-cluster correlations.}

\item{PC_out_SD}{defaulted at '5', this defines the number of standard deviation from the mean in which samples will be excluded for principle components. NA is NOT an excepted paramater.}

\item{feature_colnames_2_exclude}{names of columns to be excluded from all analysis, such as for Xenobiotics. Pass NA to this parameter to exclude this step.}

\item{derived_colnames_2_exclude}{names of columns to exclude from all sample QC steps, including independent feature identification, which is used to generate the sample PCA.}
}
\value{
a list object of:
(1) "wdata" qc'd data matrix,
(2) "featuresumstats" a list with a (1:"table") data frame of feature summary statistics and a (2:"tree") hclust object 
(3) "pca" a list with a (1:"pcs") data frame of the top 10 PCs and a binary for outliers on the top 2 PCs at 3,4, and 5 SD from the mean, a (2:"varexp") vector of the variance explained for each PC, an estimate of the number of 'significant' PCs as determined by (3:"accelerationfactor") an acceleration factor and (4:"nsig_parrallel") a parrallel analysis, and finally (5:"prob_pca") the top 10 PCs derived by a probabilistic PC analysis.
(4) "exclusion_data" a matrix of exclusion summary statistics
}
\description{
This function is a wrapper function that performs the key quality controls steps on a metabolomics data set
}
\examples{
## define a covariance matrix
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.8, 0.6, 0.2)
cmat[2,] = c(0.8, 1, 0.7, 0.5)
cmat[3,] = c(0.6, 0.7, 1, 0.6)
cmat[4,] = c(0.2, 0.5, 0.6,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
d1 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
set.seed(1010)
d2 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## simulate some random data
d3 = sapply(1:20, function(x){ rnorm(250, 40, 5) })
## define the data set
ex_data = cbind(d1,d2,d3)
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add in some missingness
ex_data[sample(1:length(ex_data), 450)] = NA
## add in some technical error to two samples
m = apply(ex_data, 2, function(x){ mean(x, na.rm = TRUE) })
ex_data[c(1,10), ] = ex_data[1, ] + (m*0.00001) 

## run the quality control
example_qc = perform.metabolite.qc(ex_data)

## the filtered data is found in 
dim(example_qc$wdata)
example_qc$wdata[1:5, 1:5]
## a data frame of summary statistics
head( example_qc$featuresumstats$table )
## a hclust dendrogram can be plotted
plot( example_qc$featuresumstats$tree , hang = -1)
abline(h = 0.5, col = "red", lwd = 1.5)
## (median imputed) PCs for all samples can be plotted
pcol = c("blue","red")[ as.factor( example_qc$pca$pcs[, "PC1_5_SD_outlier"] )]
plot(example_qc$pca$pcs[,"PC1"], example_qc$pca$pcs[,"PC2"], 
     pch = 21, cex = 1.5, bg = pcol, 
     xlab = paste0( "PC1: Var Exp = " , round(example_qc$pca$varexp[1], d = 4)*100, "\%" ) , 
     ylab = paste0( "PC2: Var Exp = " , round(example_qc$pca$varexp[2], d = 4)*100, "\%" ) )
## A Scree plot can be generated by
plot( x = 1:length(example_qc$pca$varexp), 
      y = example_qc$pca$varexp, 
      type = "b", pch = 21, cex = 2, bg = "blue",
      xlab = "PC", ylab = "Variance Explained", main = "Scree Plot")
abline(v = example_qc$pca$accelerationfactor, col = "red", lwd = 2)
abline(v = example_qc$pca$nsig_parrallel, col = "green", lwd = 2)
## Probablistic PCs for all samples can be plotted
plot(example_qc$pca$prob_pca[,1], example_qc$pca$prob_pca[,2], 
     pch = 21, cex = 1.5, bg = pcol, 
     xlab = "PC1", 
     ylab = "PC2" )
## A summary of the exclusion statistics
example_qc$exclusion_data

}
\keyword{control}
\keyword{metabolomics}
\keyword{quality}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.sum.stats.R
\name{feature.sum.stats}
\alias{feature.sum.stats}
\title{feature summary statistics}
\usage{
feature.sum.stats(
  wdata,
  sammis = NA,
  tree_cut_height = 0.5,
  outlier_udist = 5,
  feature_names_2_exclude = NA
)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{sammis}{a vector of sample missingness estimates, that is ordered to match the samples in the rows of your data matrix.}

\item{tree_cut_height}{tree cut height is the height at which to cut the feature|metabolite dendrogram to identify "independent" features. tree_cut_height is 1-absolute(Spearman's Rho) for intra-cluster correlations.}

\item{outlier_udist}{the interquartile range unit distance from the median to call a sample an outlier at a feature.}

\item{feature_names_2_exclude}{A vector of feature|metabolite names to exclude from the tree building, independent feature identification process.}
}
\value{
a list object of length two, with (1) a data frame of summary statistics and (2) a hclust object
}
\description{
This function estimates feature statistics for samples in a matrix of metabolite features.
}
\examples{
## define a covariance matrix
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.8, 0.6, 0.2)
cmat[2,] = c(0.8, 1, 0.7, 0.5)
cmat[3,] = c(0.6, 0.7, 1, 0.6)
cmat[4,] = c(0.2, 0.5, 0.6,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
d1 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
set.seed(1010)
d2 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## simulate some random data
d3 = sapply(1:20, function(x){ rnorm(250, 40, 5) })
## define the data set
ex_data = cbind(d1,d2,d3)
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add in some missingness
ex_data[sample(1:length(ex_data), 450)] = NA
## add in some technical error to two samples
m = apply(ex_data, 2, function(x){ mean(x, na.rm = TRUE) })
ex_data[c(1,10), ] = ex_data[1, ] + (m*0.00001) 
## run the function
fss = feature.sum.stats(ex_data)
## feature summary table
fss$table[1:5, ]
## plot the dendrogram
plot(fss$tree, hang = -1)

}
\keyword{metabolomics}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rntransform.R
\name{rntransform}
\alias{rntransform}
\title{rank normal tranformation}
\usage{
rntransform(y, split_ties = TRUE)
}
\arguments{
\item{y}{a numeric vector which will be rank normal transformed}

\item{split_ties}{a binary string of TRUE (default) or FALSE indicating if tied values, of the same rank, should be randomly split giving them unique ranks.}
}
\value{
returns a numeric vector, with the same length as y, of rank normal transformed values
}
\description{
This function rank normal transforms a vector of data. The procedure is built off of that provided in the GenABEL pacakge.
}
\examples{
## simulate a negative binomial distribution of values
nb_data = rnbinom(500, mu = 40, size = 100)
## rank normal transform those values
rnt_data = rntransform( nb_data , split_ties = TRUE )

}
\keyword{normal}
\keyword{rank}
\keyword{transformation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/total.peak.area.R
\name{total.peak.area}
\alias{total.peak.area}
\title{estimates total peak abundance}
\usage{
total.peak.area(wdata, feature_names_2_exclude = NA, ztransform = TRUE)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{feature_names_2_exclude}{A vector of feature|metabolite names to exclude from the tree building, independent feature identification process.}

\item{ztransform}{should the feature data be z-transformed and absolute value minimum, mean shifted prior to summing the feature values. TRUE or FALSE.}
}
\value{
a data frame of estimates for (1) total peak abundance and (2) total peak abundance at complete features for each samples
}
\description{
This function estimates total peak abundance|area for numeric data in a matrix, for (1) all features and (2) all features with complete data.
}
\examples{
set.seed(1110)
ex_data = sapply(1:5, function(x){ rnorm(10, 40, 5) })
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
ex_data[ sample(1:50, 4) ] = NA
tpa_est = total.peak.area(ex_data)

}
\keyword{abundance}
\keyword{area}
\keyword{level}
\keyword{peak}
\keyword{total}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval.power.binary.R
\name{eval.power.binary}
\alias{eval.power.binary}
\title{Estimate power for a binary variable}
\usage{
eval.power.binary(N, effect, alpha)
}
\arguments{
\item{N}{a numeric vector of total study sample size, cases and controls will both be defined as N/2.}

\item{effect}{a numeric vector of effect size}

\item{alpha}{a numeric vector of significance thresholds}
}
\value{
a matrix of parameter inputs and an estimate(s) of power are returned as a matrix
}
\description{
This function allows you estimate power for a binary variable given the sample size, effect size, significance threshold.
}
\examples{
eval.power.binary(N = 1000, effect = seq(0.01, 0.3, by = 0.01), alpha = 0.05)

}
\keyword{binary}
\keyword{estimator}
\keyword{power}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.missingness.R
\name{feature.missingness}
\alias{feature.missingness}
\title{estimate feature missingness}
\usage{
feature.missingness(wdata)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}
}
\value{
a data frame of percent missingness for each feature
}
\description{
This function estimates feature missingess, with a step to exclude poor samples identified as those with a sample missingness greater than 50%.
}
\examples{
ex_data = sapply(1:5, function(x){rnorm(10, 45, 2)})
ex_data[ sample(1:length(ex_data), 15) ] = NA
feature.missingness(wdata = ex_data )

}
\keyword{feature}
\keyword{missingness}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_normalization.R
\name{batch_normalization}
\alias{batch_normalization}
\title{median batch normalization}
\usage{
batch_normalization(
  wdata,
  feature_data_sheet = NULL,
  sample_data_sheet = NULL,
  feature_runmode_col = NULL,
  batch_ids = NULL
)
}
\arguments{
\item{wdata}{the metabolite data frame samples in row, metabolites in columns}

\item{feature_data_sheet}{a data frame containing the feature annotation data}

\item{sample_data_sheet}{a data frame containing the sample annotation data}

\item{feature_runmode_col}{a string identifying the column name in the feature_data_sheet that identifies the run mode for each feature (metabolites of proteins).}

\item{batch_ids}{a string vector, with a length equal to the number of samples in the data set that identifies what batch each sample belongs to.}
}
\value{
returns the wdata object passed to the function median normalized given the batch information provided.
}
\description{
This function median normalizes multivariable data often processed in batches, such as metabolomic and proteomic data sets.
}
\examples{
####################################
## with a vector of batch variables
####################################
## define the data set
d1 = sapply(1:10, function(x){ rnorm(25, 40, 2) })
d2 = sapply(1:10, function(x){ rnorm(25, 35, 2) })
ex_data = rbind(d1,d2)
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## define the batch
batch = c( rep("A", 25), rep("B", 25)  )
## normalize by batch
norm_wdata = batch_normalization(wdata = ex_data, batch_ids = batch )

}
\keyword{batch}
\keyword{metabolomics}
\keyword{normalization}
\keyword{proteomics}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outlier.summary.R
\name{outlier.summary}
\alias{outlier.summary}
\title{feature summary plots}
\usage{
outlier.summary(dtst, pdf_filename = "./feature_distributions.pdf", nsd = 5)
}
\arguments{
\item{dtst}{numeric data frame}

\item{pdf_filename}{name of the pdf out file}

\item{nsd}{number of SD to consider as outliers, 5 is default}
}
\value{
print summary figures for each column of data in the data frame to a pdf file.
}
\description{
This function plots the distribution of and identifiy outliers for every metabolite (column) in a data frame.
}
\examples{
## define a covariance matrix
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.8, 0.6, 0.2)
cmat[2,] = c(0.8, 1, 0.7, 0.5)
cmat[3,] = c(0.6, 0.7, 1, 0.6)
cmat[4,] = c(0.2, 0.5, 0.6,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
d1 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
set.seed(1010)
d2 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## simulate some random data
d3 = sapply(1:20, function(x){ rnorm(250, 40, 5) })
## define the data set
ex_data = cbind(d1,d2,d3)
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## run the function
outlier.summary(ex_data)

}
\keyword{feature}
\keyword{metabolomics}
\keyword{plots}
\keyword{summary}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.tree.R
\name{make.tree}
\alias{make.tree}
\title{generate a hclust dendrogram}
\usage{
make.tree(wdata, cor_method = "spearman", hclust_method = "complete")
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{cor_method}{the correlation method used in the function cor(). Default is "spearman".}

\item{hclust_method}{the dendrogram clustering method used in the construction of the tree. Default is "complete".}
}
\value{
an hclust object
}
\description{
This estimates a dendrogram of feautres based on correlation coefficeint of your choice, and a clustering method of choice
}
\examples{
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.8, 0.6, 0.2)
cmat[2,] = c(0.8, 1, 0.7, 0.5)
cmat[3,] = c(0.6, 0.7, 1, 0.6)
cmat[4,] = c(0.2, 0.5, 0.6,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
ex_data = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## estiamte the dendrogram
tree = make.tree(ex_data)
## plot the dendrogram
plot(tree, hang = -1)

}
\keyword{dendrogram}
\keyword{hclust}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/met2batch.R
\name{met2batch}
\alias{met2batch}
\title{batch effect on numeric matrix}
\usage{
met2batch(wdata, batch)
}
\arguments{
\item{wdata}{the numeric data matrix with samples in row, features in columns}

\item{batch}{a single vector containing a vector based batch variable}
}
\value{
a list object of length two with (1) a data frame of summary statistics on (a) the number of tested features ,(b) the mean batch effect across all features ,(c)  the mean batch effect across all associated (BH FDR<0.05) features ,(d) the number of associated (BH FDR<0.05) features.
}
\description{
This function estimates the effects of a categorical (batch) variable on a matrix of features in univariate linear models.
}
\examples{
d1 = sapply(1:10, function(x){ rnorm(50, 30, 5) })
d2 = sapply(1:10, function(x){ rnorm(50, 40, 5) })
ex_data = rbind(d1, d2)
d3 = sapply(1:10, function(x){ rnorm(100, 40, 5) })
ex_data = cbind(ex_data, d3)
lot = c( rep("A",50), rep("B",50) )
ex = met2batch(wdata = ex_data, batch = lot)

}
\keyword{linear}
\keyword{metabolomics}
\keyword{models}
\keyword{univariate}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_and_independent_features.R
\name{tree_and_independent_features}
\alias{tree_and_independent_features}
\title{identify independent features in a numeric matrix}
\usage{
tree_and_independent_features(
  wdata,
  minimum_samplesize = 50,
  tree_cut_height = 0.5,
  feature_names_2_exclude = NA
)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{minimum_samplesize}{the metabolite data matrix. samples in row, metabolites in columns}

\item{tree_cut_height}{the tree cut height. A value of 0.2 (1-Spearman's rho) is equivalent to saying that features with a rho >= 0.8 are NOT independent.}

\item{feature_names_2_exclude}{A vector of feature|metabolite names to exclude from this analysis. This might be features heavily present|absent like Xenobiotics or variables derived from two or more variable already in the dataset.}
}
\value{
a list object of (1) an hclust object, (2) independent features, (3) a data frame of feature ids, k-cluster identifiers, and a binary identifier of independent features
}
\description{
This function identifies independent features using Spearman's rho correlation distances, and a dendrogram tree cut step.
}
\examples{
## define a covariance matrix
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.7, 0.4, 0.2)
cmat[2,] = c(0.7, 1, 0.2, 0.05)
cmat[3,] = c(0.4, 0.2, 1, 0.375)
cmat[4,] = c(0.2, 0.05, 0.375,1)

## simulate the data (multivariable random normal)
set.seed(1110)
ex_data = MASS::mvrnorm(n = 500, mu = c(5, 45, 25, 15), Sigma = cmat )
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))

## run function to identify independent variables at a tree cut height
## of 0.5 which is equivalent to clustering variables with a Spearman's
## rho > 0.5 or (1 - tree_cut_height)
ind = tree_and_independent_features(ex_data, tree_cut_height = 0.5)

}
\keyword{features}
\keyword{independent}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.describe.R
\name{feature.describe}
\alias{feature.describe}
\title{summary statistics for features}
\usage{
feature.describe(wdata)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}
}
\value{
a data frame of summary statistics for features (columns) of a matrix
}
\description{
This function allows you to 'describe' metabolite features using the describe() function from the psych package, as well as estimate variance, a dispersion index, the coeficent of variation, and shapiro's W-statistic.
}
\examples{
ex_data = sapply(1:20, function(x){ rnorm(250, 40, 5) })
feature.describe(ex_data)

}
\keyword{feature}
\keyword{statistics}
\keyword{summary}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/id.outliers.R
\name{id.outliers}
\alias{id.outliers}
\title{identify outliers}
\usage{
id.outliers(x, nsd = 5)
}
\arguments{
\item{x}{a vector of numerical values}

\item{nsd}{the number of standard deviation from the mean outliers are identified at. Default value is 5.}
}
\value{
a vector indexing which samples are outliers
}
\description{
given a vector of data, identify those positions that are 'nsd' standard deviation units from the mean
}
\examples{
ex_data = rnbinom(500, mu = 40, size = 5)
id.outliers(ex_data, nsd = 2)

}
\keyword{outliers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcapairs_bymoose.R
\name{pcapairs_bymoose}
\alias{pcapairs_bymoose}
\title{pca pairs plot}
\usage{
pcapairs_bymoose(myloadings, varexp, pcol = "dodgerblue")
}
\arguments{
\item{myloadings}{a matrix or data frame of PC loadings (only those you would like to plot).}

\item{varexp}{a vector of the the variance explained by each PC}

\item{pcol}{plot colors for the dots background}
}
\value{
a base R plot
}
\description{
This function generates an upper triangle PCA plot for all pairs of PC loadings provided.
}
\examples{
ex_data = sapply(1:5, function(x){rnorm( 50, 0, 2) })
## add in some extreme values to a sample
ex_data[1,] = ex_data[1,] + sample( c(8, -8), ncol(ex_data), replace = TRUE )
ex_data[2,] = ex_data[2,] + sample( c(3, -3), ncol(ex_data), replace = TRUE )
## plot
pcapairs_bymoose(myloadings = ex_data, 
    varexp = rep(1/ncol(ex_data), ncol(ex_data)),
    pcol = "tomato" )

}
\keyword{PCA}
\keyword{metabolomics}
\keyword{pairs}
\keyword{plot}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.outliers.R
\name{sample.outliers}
\alias{sample.outliers}
\title{outlier features count for samples}
\usage{
sample.outliers(wdata, nsd = 5)
}
\arguments{
\item{wdata}{a metabolite data matrix with samples in row, metabolites in columns}

\item{nsd}{the number of standard deviation from the mean outliers are identified at. The default value is 5.}
}
\value{
a data frame of outiler counts for each sample
}
\description{
This function takes a matrix of numeric data and counts the number of outlying features each sample has.
}
\examples{
d = sapply(1:5, function(x){ rnorm(50, 50, 15) })
sample.outliers(d, nsd = 2)

}
\keyword{metabolomics}
\keyword{outliers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/median_impute.R
\name{median_impute}
\alias{median_impute}
\title{median impute missing data}
\usage{
median_impute(wdata)
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}
}
\value{
the matrix passed to the function but with NA's imputed to each columns median value.
}
\description{
This function imputes features (columns) of a metabolome matrix to median estimates. Useful for PCA.
}
\examples{
ex_data = sapply(1:100, function(x){ rnorm(250, 40, 5) })
## define the data set
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add in some missingness
ex_data[sample(1:length(ex_data), 500)] = NA
## Estimate missingness and generate plots
imp_data = median_impute(ex_data)

}
\keyword{imputation}
\keyword{median}
\keyword{metabolomics}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval.power.cont.R
\name{eval.power.cont}
\alias{eval.power.cont}
\title{estimate power for continuous variable}
\usage{
eval.power.cont(N, n_coeff, effect, alpha)
}
\arguments{
\item{N}{Sample size}

\item{n_coeff}{degrees of freedom for numerator}

\item{effect}{effect size}

\item{alpha}{significance level (Type 1 error)}
}
\description{
This function estimates power for a continuous variable given the sample size, effect size, significance threshold, and the degrees of freedom.
}
\examples{
eval.power.cont(N = 1000, n_coeff = 1, effect = 0.0025, alpha = 0.05)

}
\keyword{continuous}
\keyword{power}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outlier.matrix.R
\name{outlier.matrix}
\alias{outlier.matrix}
\title{identify outlier sample indexes in a matrix}
\usage{
outlier.matrix(data, nsd = 5, meansd = FALSE)
}
\arguments{
\item{data}{a matrix of numerical values, samples in row, features in columns}

\item{nsd}{the unit distance in SD or IQR from the mean or median estimate, respectively outliers are identified at. Default value is 5.}

\item{meansd}{set to TRUE if you would like to estimate outliers using a mean and SD method; set to FALSE if you would like to estimate medians and inter quartile ranges. The default is FALSE.}
}
\value{
a matrix of 0 (not a sample outlier) and 1 (outlier)
}
\description{
Given a matrix of data this function returns a matrix of 0|1, of the same structure with 1 values indicating outliers. It is an expansion of the function id.outliers(), applied to columns of a matrix.
}
\examples{
ex_data = sapply(1:25, function(x){ rnorm(250, 40, 5) })
## define the data set
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## add in some technical error to two samples
m = apply(ex_data, 2, function(x){ mean(x, na.rm = TRUE) })
ex_data[c(1,50), ] = ex_data[1, ] + (m*4) 
Omat = outlier.matrix(ex_data)
## how many outliers identified
sum(Omat)

}
\keyword{indexes}
\keyword{matrix}
\keyword{outlier}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_report.R
\name{generate_report}
\alias{generate_report}
\title{generate metaboprep summary html report}
\usage{
generate_report(
  full_path_2_Rdatafile = "ReportData.Rdata",
  dir_4_report = "./",
  path_2_Rmd_template = file.path(system.file("rmarkdown", package = "metaboprep"),
    "metaboprep_Report_v0.Rmd")
)
}
\arguments{
\item{full_path_2_Rdatafile}{full path to the Rdatafile}

\item{dir_4_report}{directory to place the report}

\item{path_2_Rmd_template}{full path to the html report template}
}
\value{
Writes a html report to file

an html file written to file
}
\description{
This function generates the html report.
}
\examples{


}
\keyword{html}
\keyword{knit}
\keyword{report}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.in.nightingale.R
\name{read.in.nightingale}
\alias{read.in.nightingale}
\title{read in Nightingale Health metabolomics data}
\usage{
read.in.nightingale(file2process, data_dir, projectname)
}
\arguments{
\item{file2process}{the name of the xls file to process}

\item{data_dir}{the full path to the directory holding your Nightingale excel file}

\item{projectname}{a name for your project}
}
\value{
a list object of (1) metabolite, (2) sample annotation, and (3) feature annotation data
}
\description{
This function reads in a Nightingale raw data excel file, writes the (1) metabolite, (2) sample annotation, and (3) feature annotation data to flat text files. It also returns a list object of the same data.
}
\examples{
# read.in.nightingale(file2process = "NH_data_release.xls", 
#  data_dir = "/File/sits/here/", 
#  projectname = "My Amazing Project")

}
\keyword{Health}
\keyword{Nigtingale}
\keyword{meatbolomics}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outliers.R
\name{outliers}
\alias{outliers}
\title{identify outliers}
\usage{
outliers(x, nsd = 3)
}
\arguments{
\item{x}{a numerical vector of data}

\item{nsd}{the number of SD units from the mean to be used as an outlier cutoff.}
}
\value{
a list object of length three. (1) a vector of sample indexes indicating the outliers, (2) the lower outlier cuttoff value, (3) the upper outlier cuttoff value.
}
\description{
This function identifies outliers from a vector of data at SD units from the mean.
}
\examples{
ex_data = rnbinom(500, mu = 40, size = 5)
outliers(ex_data)

}
\keyword{outliers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate.anova.R
\name{multivariate.anova}
\alias{multivariate.anova}
\title{multivariate analysis}
\usage{
multivariate.anova(dep, indep_df)
}
\arguments{
\item{dep}{a vector of the dependent variable}

\item{indep_df}{a data frame of the independent variable}
}
\value{
ggplot2 table figure of
}
\description{
This function performs a multivariate analysis over a dependent|response and numerous independent|explanatory variable
}
\examples{
cmat = matrix(1, 3, 3 )
cmat[1,] = c(1, 0.5, 0.3)
cmat[2,] = c(0.5, 1, 0.25)
cmat[3,] = c(0.3, 0.25, 1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
ex_data = MASS::mvrnorm(n = 250, mu = c(5, 45, 25), Sigma = cmat )
colnames(ex_data) = c("outcome","age","bmi")
multivariate.anova(dep = ex_data[,1], indep_df = ex_data[, 2:3])

}
\keyword{ANOVA}
\keyword{metabolomics}
\keyword{multivariate}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.cor.matrix.R
\name{make.cor.matrix}
\alias{make.cor.matrix}
\title{correlation matrix}
\usage{
make.cor.matrix(wdata, cor_method = "kendall", minN = 50, var2return = "cor")
}
\arguments{
\item{wdata}{the metabolite data matrix. samples in row, metabolites in columns}

\item{cor_method}{defaulted to "kendall" this is the correlation method to use in the function cor.test()}

\item{minN}{sefaulted to 50, this is the minimum number of observations that must be available in pairs to perform analysis}

\item{var2return}{sefaulted to "cor", other option is "pvalue" is a the flag indicating which estimate to return from the function.}
}
\value{
a matrix of correlation estimates or p-values
}
\description{
This function estimates a correlation matrix returning wither the correlation estimates or their p-values
}
\examples{
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.8, 0.6, 0.2)
cmat[2,] = c(0.8, 1, 0.7, 0.5)
cmat[3,] = c(0.6, 0.7, 1, 0.6)
cmat[4,] = c(0.2, 0.5, 0.6,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
ex_data = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## return correlation estimates
cor_mat = make.cor.matrix(ex_data, var2return = "cor")
## return p-values
cor_mat = make.cor.matrix(ex_data, var2return = "pvalue")

}
\keyword{correlation}
\keyword{matrix}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variable.by.factor.R
\name{variable.by.factor}
\alias{variable.by.factor}
\title{ggplot2 violin plot}
\usage{
variable.by.factor(
  dep,
  indep,
  dep_name = "dependent",
  indep_name = "independent",
  orderfactor = TRUE,
  violin = TRUE
)
}
\arguments{
\item{dep}{a vector of the dependent variable}

\item{indep}{a vector of the independent variable}

\item{dep_name}{name of the dependent variable}

\item{indep_name}{name of the independent variable}

\item{orderfactor}{order factors alphebetically}

\item{violin}{box plot or violin plot. violin = TRUE is default}
}
\value{
a ggplot2 object
}
\description{
This function performs univariate linear analysis of a dependent and an independent variable and generates a viloin or box plot to illustrate the associated structure.
}
\examples{
x = c( rnorm(20, 10, 2), rnorm(20, 20, 2) )
y = as.factor( c( rep("A", 20), rep("B", 20)  ) )
variable.by.factor(dep = x , indep = y, dep_name = "expression", indep_name = "species" )

}
\keyword{ggplot}
\keyword{metabolomics}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.pa.imbalanced.power.make.plot.R
\name{run.pa.imbalanced.power.make.plot}
\alias{run.pa.imbalanced.power.make.plot}
\title{binary trait imbalanced design power analysis plot}
\usage{
run.pa.imbalanced.power.make.plot(mydata)
}
\arguments{
\item{mydata}{a numeric data matrix with samples in rows and features in columns}
}
\value{
a ggplot2 object
}
\description{
This function (1) estimates an informative distribution of effect and power estimates given your datas total sample size, over a distribution of imbalanced sample sizes and (2) returns a summary plot.
}
\examples{
ex_data = matrix(NA, 1000, 2)
run.pa.imbalanced.power.make.plot( ex_data )

}
\keyword{analysis}
\keyword{binary}
\keyword{plot}
\keyword{power}
\keyword{trait}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.factor.analysis.R
\name{pca.factor.analysis}
\alias{pca.factor.analysis}
\title{PCA factor analysis and annotation enrichment}
\usage{
pca.factor.analysis(
  metabolitedata,
  pcloadings,
  sigthreshold = 0.3,
  feature_anno = feature_data$SUPER_PATHWAY
)
}
\arguments{
\item{metabolitedata}{a matrix or data frame of metabolite data}

\item{pcloadings}{a matrix or data frame of pc loadings you wish to test}

\item{sigthreshold}{Spearman's rho correlation threshold to declare an association between the numeric variable and a PC loading}

\item{feature_anno}{a vector of variable annotations to perform the hypergeomtric enrichment on}
}
\value{
a list object of length 2: (1) a list object of enrichment_tables, and (2) the Spearman's correlation matrix between features and the PC loadings
}
\description{
This function performs (1) a factor analysis on numeric data and PC loadings derived from said data, then subsequently (3) performs a hypergeometrix enrichment test to ask if a certain class of variables are enriched on a particular PC.
}
\examples{
## define a covariance matrix
cmat = matrix(1, 4, 4 )
cmat[1,] = c(1, 0.8, 0.6, 0.2)
cmat[2,] = c(0.8, 1, 0.7, 0.5)
cmat[3,] = c(0.6, 0.7, 1, 0.6)
cmat[4,] = c(0.2, 0.5, 0.6,1)
## simulate some correlated data (multivariable random normal)
set.seed(1110)
d1 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
set.seed(1010)
d2 = MASS::mvrnorm(n = 250, mu = c(5, 45, 25, 15), Sigma = cmat )
## simulate some random data
d3 = sapply(1:20, function(x){ rnorm(250, 40, 5) })
ex_data = cbind(d1,d2,d3)
rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
colnames(ex_data) = paste0("var", 1:ncol(ex_data))
## annotation
met_anno = c( rep("A", 10), rep("B", 10), rep("C", 8) )
## PCA
pca = prcomp(ex_data, center = TRUE, scale = TRUE)
## run pca.factor.analysis()
ex_out = pca.factor.analysis(metabolitedata = ex_data, 
                             pcloadings = pca$x[,1:5], 
                             sigthreshold = 0.3, 
                             feature_anno = met_anno )


}
\keyword{PC}
\keyword{analysis}
\keyword{enrichment}
\keyword{exact}
\keyword{factor}
\keyword{fisher}
\keyword{hypergeometric}
\keyword{test}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sam.missingness.exclusion.R
\name{sam.missingness.exclusion}
\alias{sam.missingness.exclusion}
\title{sample exlusions on missingness and total peak area}
\usage{
sam.missingness.exclusion(mydata, sdata, fdata)
}
\arguments{
\item{mydata}{metabolite data}

\item{sdata}{sample data}

\item{fdata}{feature data}
}
\value{
a data frame of missingness and TPA exclusions
}
\description{
This function provides missingnes and tpa estimates along with exlcusion at 3, 4, and 5 SD from the mean.
}
\examples{
# sam.missingness.exclusion()

}
\keyword{metabolomics}
