
# Consensus Nested Cross-Validation


#### Websites

[insilico Github Organization](https://github.com/insilico)

[insilico McKinney Lab](http://insilico.utulsa.edu/)

#### Related References. 

[2020 Consensus Features Nested Cross-Validation paper in Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa046)

[2018 EpistasisRank and EpistasisKatz paper in Bioinformatics](https://doi.org/10.1093/bioinformatics/bty965)

### To install:

    >library(devtools)
    
    >install_github("insilico/cncv")  # todo (optional build_vignettes = TRUE)
    >library(cncv)
    >data(package="cncv")
    
    # >vignette(" ") # todo (if you build_vignettes)
    
### Dependencies
To install the `privateEC` collection of R packages:

```
install_github("insilico/privateEC")
```
To install the `PriorKnowledgeEpistasisRank` collection of R packages:

```
install_github("insilico/PriorKnowledgeEpistasisRank")
```

Other packages
```
install.packages(c('randomForest', 'glmnet', 'xgboost', 'caret', 'CORElearn))

```

### Examples

[https://github.com/insilico/cncv/blob/master/Vignettes/cncvExample.md](https://github.com/insilico/cncv/blob/master/Vignettes/cncvExample.md)



### Abstract

Motivation: Appropriate steps must be taken to avoid overfitting when using feature selection to improve the accuracy of machine learning models. Nested cross-validation (nCV) is a common approach that chooses the best classification model and best features to represent a given outer fold based on the maximum inner-fold accuracy. Differential privacy is a related technique to avoid overfitting that uses a privacy preserving noise mechanism to identify features that are stable between training and holdout sets. 
Methods: We develop consensus nested CV (cnCV) that combines the idea of feature stability from differential privacy with nested CV. Feature selection is applied in each inner fold and the consensus of top features across folds is a used as a measure of feature stability or reliability instead of classification used in standard nCV. We use simulated data with main effects, correla-tion, and interactions to compare the classification accuracy and feature selection performance of the new cnCV with standard nCV, Elastic Net optimized by CV, differential privacy, and pri-vate Evaporative Cooling (pEC). We also compare these methods using real RNA-Seq data from a study of major depressive disorder.
Results: The cnCV method has similar training and validation accuracy to nCV, but cnCV has much shorter run times because it does not construct classifiers in the inner folds. The cnCV method chooses a more parsimonious set of features with fewer false positives than nCV. The cnCV method has similar accuracy to pEC and cnCV selects stable features between folds without the need to specify a privacy threshold. We show that cnCV is an effective and efficient approach for combining feature selection with classification. 


#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
Consensus Nested Cross-Validation Vignette With Simulated Data
================
Saeid Parvandeh and Brett McKinney 
December 2019

Install privateEC:
---------------------------

``` r
knitr::opts_chunk$set(echo = TRUE, warning=FALSE)
rm(list = ls())

if (!("devtools" %in% installed.packages()[,"Package"])){
  install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
library(devtools)

if (!("privateEC" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/privateEC") # build_vignettes = TRUE)
}
if (!("cncv" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/cncv", build_vignettes = TRUE)
}
library(privateEC)  # used to simulate data
library(cncv)
```

Simulate data with privateEC
----------------------------

``` r
letsSimulate <- T   # F to use previously simulated data
class.lab <- "class"
writeData <- F  # usually the same as letsSimulate
writeResults <- F

num.samp <- 300
num.attr <- 1000
pct.signals <- 0.1
bias <- 0.4
#sim.type <- "mainEffect"
sim.type <- "interactionErdos"
importance.algorithm = "ReliefFequalK"
num_tree = 500
verbose = T


pec_simFile <- paste("pec_simulated", sim.type, "bias", bias, 
                             "pct.signals", pct.signals,
                             "num.attr", num.attr, "num.samp", num.samp, sep = "_")
pec_simFile <- paste(pec_simFile,".csv",sep="")

if (letsSimulate == TRUE){
    sim.data <- createSimulation(num.samples = num.samp, num.variables = num.attr,
                                 pct.signals = pct.signals, pct.train = 1/3, pct.holdout = 1/3, 
                                 pct.validation = 1/3, bias = bias, sim.type = sim.type, verbose = FALSE)
  dat <- rbind(sim.data$train, sim.data$holdout)
  predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
} else { # optional: use provided data
  dat <- read.csv(pec_simFile)
  dat <- dat[,-1] # written file has first X column with subject names
  predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
}

dat[, class.lab] <- as.factor(dat[, class.lab]) 
pheno.class <- dat[, class.lab]
attr.names <- colnames(predictors.mat)
num.samp <- nrow(dat)

if (writeData == TRUE){
  write.csv(dat, file = pec_simFile)
}
```

### Run Standard nested CV

``` r
rncv_result <- regular_nestedCV(train.ds = sim.data$train,
                                validation.ds =  sim.data$holdout,
                                label = sim.data$label,
                                method.model = "classification",
                                is.simulated = TRUE,
                                ncv_folds = c(10, 10),
                                param.tune = FALSE,
                                learning_method = "rf", 
                                importance.algorithm = importance.algorithm,
                                relief.k.method = "k_half_sigma",             # ReliefF knn
                                wrapper = "relief",
                                inner_selection_percent = 0.2,
                                inner_selection_positivescores = TRUE,
                                tuneGrid = NULL,
                                num_tree = num_tree,
                                verbose = verbose)
```

### nested CV results
``` r
cat("\n Train Accuracy [",rncv_result$cv.acc,"]\n")
cat("\n Validation Accuracy [",rncv_result$Validation,"]\n")
cat("\n Selected Features \n [",rncv_result$Features,"]\n")
cat("\n Elapsed Time [",rncv_result$Elapsed,"]\n")
```

### Run Consensus nested CV
``` r
cncv_result <- consensus_nestedCV(train.ds = rbind(data.sets$train,data.sets$holdout), 
                                  validation.ds =  data.sets$validation, 
                                  label = data.sets$label,
                                  method.model = "classification",
                                  is.simulated = TRUE,
                                  ncv_folds = c(10, 10),
                                  param.tune = FALSE,
                                  learning_method = "rf", 
                                  importance.algorithm = importance.algorithm,
                                  relief.k.method = "k_half_sigma",             # ReliefF knn
                                  wrapper = "relief",
                                  inner_selection_percent = 0.2,
                                  inner_selection_positivescores = TRUE,
                                  tune.k = FALSE,
                                  tuneGrid = NULL,
                                  num_tree = num_tree,
                                  verbose = verbose)
```

### cnCV results

``` r
cat("\n Nested Cross-Validation Accuracy [",cncv_result$cv.acc,"]\n")
cat("\n Validation Accuracy [",cncv_result$Validation,"]\n")
cat("\n Selected Features \n [",cncv_result$Features,"]\n")
cat("\n Elapsed Time [",cncv_result$Elapsed,"]\n")
```

Consensus Nested Cross-Validation Vignette With Simulated Data
================
Saeid Parvandeh and Brett McKinney 
December 2019

Install privateEC:
---------------------------

``` r
knitr::opts_chunk$set(echo = TRUE, warning=FALSE)
rm(list = ls())

if (!("devtools" %in% installed.packages()[,"Package"])){
  install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
library(devtools)

if (!("privateEC" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/privateEC") # build_vignettes = TRUE)
}
if (!("cncv" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/cncv", build_vignettes = TRUE)
}
library(privateEC)  # used to simulate data
library(cncv)
```

Simulate data with privateEC
----------------------------

``` r
letsSimulate <- T   # F to use previously simulated data
class.lab <- "class"
writeData <- F  # usually the same as letsSimulate
writeResults <- F

num.samp <- 300
num.attr <- 1000
pct.signals <- 0.1
bias <- 0.4
#sim.type <- "mainEffect"
sim.type <- "interactionErdos"
importance.algorithm = "ReliefFequalK"
num_tree = 500
verbose = T

pec_simFile <- paste("pec_simulated", sim.type, "bias", bias, 
                             "pct.signals", pct.signals,
                             "num.attr", num.attr, "num.samp", num.samp, sep = "_")
pec_simFile <- paste(pec_simFile,".csv",sep="")

if (letsSimulate == TRUE){
    sim.data <- createSimulation(num.samples = num.samp, num.variables = num.attr,
                                 pct.signals = pct.signals, pct.train = 1/3, pct.holdout = 1/3, 
                                 pct.validation = 1/3, bias = bias, sim.type = sim.type, verbose = FALSE)
  dat <- rbind(sim.data$train, sim.data$holdout)
  predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
} else { # optional: use provided data
  dat <- read.csv(pec_simFile)
  dat <- dat[,-1] # written file has first X column with subject names
  predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
}

dat[, class.lab] <- as.factor(dat[, class.lab]) 
pheno.class <- dat[, class.lab]
attr.names <- colnames(predictors.mat)
num.samp <- nrow(dat)

if (writeData == TRUE){
  write.csv(dat, file = pec_simFile)
}
```

### Run Standard nested CV

``` r
rncv_result <- regular_nestedCV(train.ds = sim.data$train,
                                validation.ds =  sim.data$holdout,
                                label = sim.data$label,
                                method.model = "classification",
                                is.simulated = TRUE,
                                ncv_folds = c(10, 10),
                                param.tune = FALSE,
                                learning_method = "rf", 
                                importance.algorithm = importance.algorithm,
                                relief.k.method = "k_half_sigma",             # ReliefF knn
                                wrapper = "relief",
                                inner_selection_percent = 0.2,
                                inner_selection_positivescores = TRUE,
                                tuneGrid = NULL,
                                num_tree = num_tree,
                                verbose = verbose)
```

### nested CV results
``` r
cat("\n Train Accuracy [",rncv_result$cv.acc,"]\n")
cat("\n Validation Accuracy [",rncv_result$Validation,"]\n")
cat("\n Selected Features \n [",rncv_result$Features,"]\n")
cat("\n Elapsed Time [",rncv_result$Elapsed,"]\n")
```

### Run Consensus nested CV
``` r
cncv_result <- consensus_nestedCV(train.ds = rbind(sim.data$train,sim.data$holdout), 
                                  validation.ds =  sim.data$validation, 
                                  label = sim.data$label,
                                  method.model = "classification",
                                  is.simulated = TRUE,
                                  ncv_folds = c(10, 10),
                                  param.tune = FALSE,
                                  learning_method = "rf", 
                                  importance.algorithm = importance.algorithm,
                                  relief.k.method = "k_half_sigma",             # ReliefF knn
                                  wrapper = "relief",
                                  inner_selection_percent = 0.2,
                                  inner_selection_positivescores = TRUE,
                                  tune.k = FALSE,
                                  tuneGrid = NULL,
                                  num_tree = num_tree,
                                  verbose = verbose)
```

### cnCV results

``` r
cat("\n Nested Cross-Validation Accuracy [",cncv_result$cv.acc,"]\n")
cat("\n Validation Accuracy [",cncv_result$Validation,"]\n")
cat("\n Selected Features \n [",cncv_result$Features,"]\n")
cat("\n Elapsed Time [",cncv_result$Elapsed,"]\n")
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nestedCV.R
\name{consensus_nestedCV}
\alias{consensus_nestedCV}
\title{Consensus nested cross validation for feature selection and parameter tuning}
\usage{
consensus_nestedCV(
  train.ds = NULL,
  validation.ds = NULL,
  label = "class",
  method.model = "classification",
  is.simulated = TRUE,
  ncv_folds = c(10, 10),
  param.tune = FALSE,
  learning_method = "rf",
  xgb.obj = "binary:logistic",
  importance.algorithm = "ReliefFequalK",
  wrapper = "relief",
  inner_selection_percent = NULL,
  inner_selection_positivescores = TRUE,
  tune.inner_selection_percent = NULL,
  tune.k = FALSE,
  tuneGrid = NULL,
  relief.k.method = "k_half_sigma",
  num_tree = 500,
  covars_vec = NULL,
  covars.pval.adj = 0.05,
  verbose = FALSE
)
}
\arguments{
\item{train.ds}{A training data frame with last column as outcome}

\item{validation.ds}{A validation data frame with last column as outcome}

\item{label}{A character vector of the outcome variable column name.}

\item{method.model}{Column name of outcome variable (string), classification or regression. If the analysis goal is classification make the column a factor type.
For regression, make outcome column numeric type.}

\item{is.simulated}{A TRUE or FALSE character for data type}

\item{ncv_folds}{A numeric vector to indicate nested cv folds: c(k_outer, k_inner)}

\item{param.tune}{A TRUE or FALSE character for tuning parameters}

\item{learning_method}{Name of the method: glmnet/xgbTree/rf}

\item{importance.algorithm}{A character vestor containing a specific importance algorithm subtype}

\item{wrapper}{feature selection algorithm including: rf, glmnet, t.test, centrality methods (PageRank, Katz,
EpistasisRank, and EpistasisKatz from Rinbix packages), ReliefF family, and etc.}

\item{inner_selection_percent}{= Percentage of features to be selected in each inner fold.}

\item{inner_selection_positivescores}{A TRUE or FALSE character to select positive scores (if the value is False, use the percentage method).}

\item{tune.inner_selection_percent}{A sequence vector of possible percentages for tuning}

\item{tune.k}{A sequence vector to tune k nearest neighbors in relief method, if TRUE the default grid is seq(1,kmax), where kmax=floor((m-1)/2)
and m is the number of samples. However, this kmax is for balanced data. If data are imbalance, where m_minority + m_majority = m, then kmax = floor(m_minority-1).
Default is FALSE.}

\item{tuneGrid}{A data frame with possible tuning values. The columns are named the same as the tuning parameters.
This caret library parameter, for more information refer to  http://topepo.github.io/caret/available-models.html.}

\item{relief.k.method}{A character of numeric to indicate number of nearest neighbors for relief algorithm.
Possible characters are: k_half_sigma (floor((num.samp-1)*0.154)), m6 (floor(num.samp/6)),
myopic (floor((num.samp-1)/2)), and m4 (floor(num.samp/4))}

\item{num_tree}{Number of trees in random forest and xgboost methods}

\item{verbose}{A flag indicating whether verbose output be sent to stdout}
}
\value{
A list with:
\describe{
\item{cv.acc}{Training data accuracy}
\item{Validation}{Validation data accuracy}
\item{Features}{number of variables detected correctly in nested cross validation}
\item{Train_model}{Traing model to use for validation}
\item{Elapsed}{total elapsed time}
}
num.samples <- 100
num.variables <- 100
pct.signals <- 0.1
label <- "class"
sim.data <- createSimulation(num.samples = num.samples,
num.variables = num.variables,
pct.signals = pct.signals,
sim.type = "mainEffect",
label = label,
verbose = FALSE)
cnCV.results <- consensus_nestedCV(train.ds = sim.data$train,
validation.ds = sim.data$holdout,
label = label,
is.simulated = TRUE,
ncv_folds = c(10, 10),
param.tune = FALSE,
learning_method = "rf",
importance.algorithm = "ReliefFbestK",
num_tree = 500,
verbose = FALSE)
}
\description{
Consensus nested cross validation for feature selection and parameter tuning
}
\seealso{
Other nestedCV: 
\code{\link{functionalDetectionStats}()},
\code{\link{regular_nestedCV}()}
}
\concept{nestedCV}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nestedCV.R
\name{functionalDetectionStats}
\alias{functionalDetectionStats}
\title{Given a vector functional (true) attribute associations and a vector of positive associations,
returns detection statistics like true positive rate, recall and precision.}
\usage{
functionalDetectionStats(functional, positives)
}
\arguments{
\item{functional}{A vector functional (true) attributes}

\item{positives}{A vector of positive associations}
}
\value{
A list with:
\describe{
\item{TP}{True Positive}
\item{FP}{False Positive}
\item{FN}{False Negative}
\item{TPR}{True Positive Rate}
\item{FPR}{False Positive Rate}
\item{precision}{Fraction of relevant features among the retrieved features}
\item{recall}{Fraction of total relevant features that were actually retrieved}
}
num.samples <- 100
num.variables <- 100
pct.signals <- 0.1
label <- "class"
sim.data <- createSimulation(num.samples = num.samples,
num.variables = num.variables,
pct.signals = pct.signals,
sim.type = "mainEffect",
label = label,
verbose = FALSE)
rnCV.results <- regular_nestedCV(train.ds = sim.data$train,
validation.ds = sim.data$holdout,
label = label,
is.simulated = TRUE,
ncv_folds = c(10, 10),
param.tune = FALSE,
learning_method = "rf",
importance.algorithm = "ReliefFbestK",
num_tree = 500,
verbose = FALSE)
functional <- data.sets$signal.names
rncv.positives <- rncv_result$Features
rncv.detect <- functionalDetectionStats(functional, rncv.positives)
}
\description{
Given a vector functional (true) attribute associations and a vector of positive associations,
returns detection statistics like true positive rate, recall and precision.
}
\seealso{
Other nestedCV: 
\code{\link{consensus_nestedCV}()},
\code{\link{regular_nestedCV}()}
}
\concept{nestedCV}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nestedCV.R
\name{regular_nestedCV}
\alias{regular_nestedCV}
\title{Regular nested cross validation for feature selection and parameter tuning}
\usage{
regular_nestedCV(
  train.ds = NULL,
  validation.ds = NULL,
  label = "class",
  method.model = "classification",
  is.simulated = TRUE,
  ncv_folds = c(10, 10),
  param.tune = FALSE,
  learning_method = "rf",
  xgb.obj = "binary:logistic",
  importance.algorithm = "ReliefFequalK",
  wrapper = "relief",
  inner_selection_percent = 0.2,
  inner_selection_positivescores = TRUE,
  relief.k.method = "k_half_sigma",
  tuneGrid = NULL,
  num_tree = 500,
  verbose = FALSE
)
}
\arguments{
\item{train.ds}{A training data frame with last column as outcome}

\item{validation.ds}{A validation data frame with last column as outcome}

\item{label}{A character vector of the outcome variable column name.}

\item{method.model}{Column name of outcome variable (string), classification or regression. If the analysis goal is classification make the column a factor type.
For regression, make outcome column numeric type.}

\item{is.simulated}{A TRUE or FALSE character for data type}

\item{ncv_folds}{A numeric vector to indicate nested cv folds: c(k_outer, k_inner)}

\item{param.tune}{A TRUE or FALSE character for tuning parameters}

\item{learning_method}{Name of the method: glmnet/xgbTree/rf}

\item{xgb.obj}{Name of xgboost algorithm}

\item{importance.algorithm}{A character vestor containing a specific importance algorithm subtype}

\item{wrapper}{feature selection algorithm including: rf, glmnet, t.test, centrality methods (PageRank, Katz,
EpistasisRank, and EpistasisKatz from Rinbix packages), ReliefF family, and etc.}

\item{inner_selection_percent}{= .2 Percentage of features to be selected in each inner fold.}

\item{inner_selection_positivescores}{A TRUE or FALSE character to select positive scores (if the value is False, use the percentage method).}

\item{relief.k.method}{A character of numeric to indicate number of nearest neighbors for relief algorithm.
Possible characters are: k_half_sigma (floor((num.samp-1)*0.154)), m6 (floor(num.samp/6)),}

\item{tuneGrid}{A data frame with possible tuning values. The columns are named the same as the tuning parameters.
This caret library parameter, for more information refer to  http://topepo.github.io/caret/available-models.html.
myopic (floor((num.samp-1)/2)), and m4 (floor(num.samp/4))}

\item{num_tree}{Number of trees in random forest and xgboost methods}

\item{verbose}{A flag indicating whether verbose output be sent to stdout}
}
\value{
A list with:
\describe{
\item{cv.acc}{Training data accuracy}
\item{Validation}{Validation data accuracy}
\item{Features}{number of variables detected correctly in nested cross validation}
\item{Train_model}{Traing model to use for validation}
\item{Elapsed}{total elapsed time}
}
num.samples <- 100
num.variables <- 100
pct.signals <- 0.1
label <- "class"
sim.data <- createSimulation(num.samples = num.samples,
num.variables = num.variables,
pct.signals = pct.signals,
sim.type = "mainEffect",
label = label,
verbose = FALSE)
rnCV.results <- regular_nestedCV(train.ds = sim.data$train,
validation.ds = sim.data$holdout,
label = label,
is.simulated = TRUE,
ncv_folds = c(10, 10),
param.tune = FALSE,
learning_method = "rf",
importance.algorithm = "ReliefFbestK",
num_tree = 500,
verbose = FALSE)
}
\description{
Regular nested cross validation for feature selection and parameter tuning
}
\seealso{
Other nestedCV: 
\code{\link{consensus_nestedCV}()},
\code{\link{functionalDetectionStats}()}
}
\concept{nestedCV}
