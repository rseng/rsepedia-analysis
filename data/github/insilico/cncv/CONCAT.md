
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

