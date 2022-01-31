![R CMD check](https://github.com/mmahmoudian/sivs/workflows/R%20CMD%20check/badge.svg)
![R CMD check --as-cran](https://github.com/mmahmoudian/sivs/workflows/R%20CMD%20check%20--as-cran/badge.svg)
![download per month](https://cranlogs.r-pkg.org/badges/sivs)


#  Stable Iterative Variable Selection (SIVS) <img src="misc/img/SIVS_logo.png" width="140" align="right" />

SIVS is an acronym of Stable Iterative Variable Selection, and as the name suggests is a feature selection method that is robust to the variations that cross-validation can have on various methods with embedded feature selection. This method hired an iterative approach and  internally utilizes varius Machine Learning methods which have embedded feature reduction in order to shrink down the feature space into a small and yet robust set.

For citation information, see the [citation section](#citation) of this document.


## Installation

You can download and install the latest stable version from CRAN via:

```r
install.packages("sivs", repos = "https://cran.rstudio.com")
```

Alternatively, you can install it directly from github via either of the following:

```r
### First Approach
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mmahmoudian/sivs")
```

```r
### Second Approach
if (!require("remotes")) install.packages("remotes")
remotes::install_github('mmahmoudian/sivs')
```


## Building From Source

You can also build this package completelty from source and you are expected to get identical files as in CRAN. This can be useful for those who want to contribute to the package. I have made it easy and straight-forward to build and test the package using the [GNU make](https://www.gnu.org/software/make/). Follow these steps in order:

0. First make sure you have the `make` installed and the package building dependencies
   ```sh
   # if make is installed, you will see the version
   make --version
   
    # this will check if you have the needed R packages and if not, it will install them for you
   make deps
   ```
1. Change the code and files as needed
2. if you have changed the R code and want to test it, you can build the R code and skip building the manual and vignette:
   ```sh
   make build-noman
   ```
   if you have changed the manual:
   ```sh
   make docs
   make build
   ```
3. install the package and make sure things are in order and working as you expected:
   ```sh
   make install
   ```
4. When you confirmed that everything is in order, repeat all the building steps with CRAN checking:
   ```sh
   make docs build check-cran install
   ```
   Alternatively you can run the following which is short-form for the command above:
   ```sh
   make all-cran
   ```


## Contribution

This is a [Free and Libré OpenSource Software (FLOSS)](https://en.wikipedia.org/wiki/Free_and_open-source_software) and therefore any contribution is welcome as long as it does not violate [the license](https://github.com/mmahmoudian/sivs/blob/master/LICENSE). To contribute, follow the steps in the [Building From Source](#building-from-source) and then before creating the pull-request, make sure you have solved all ERRORs, WARNINGs and possibly all NOTEs produced by the following:

```sh
make all-cran
make check-cran`
```

## Citation

This method has been published in the journal of Bioinformatics:

Mehrad Mahmoudian, Mikko S Venäläinen, Riku Klén, Laura L Elo, Stable Iterative Variable Selection, Bioinformatics, 2021;, btab501, https://doi.org/10.1093/bioinformatics/btab501

BibTeX entry for LaTeX users:

```
@article{10.1093/bioinformatics/btab501,
    author = {Mahmoudian, Mehrad and Venäläinen, Mikko S and Klén, Riku and Elo, Laura L},
    title = "{Stable Iterative Variable Selection}",
    journal = {Bioinformatics},
    year = {2021},
    month = {07},
    abstract = "{The emergence of datasets with tens of thousands of features, such as high-throughput omics biomedical data, highlights the importance of reducing the feature space into a distilled subset that can truly capture the signal for research and industry by aiding in finding more effective biomarkers for the question in hand. A good feature set also facilitates building robust predictive models with improved interpretability and convergence of the applied method due to the smaller feature space.Here, we present a robust feature selection method named Stable Iterative Variable Selection (SIVS) and assess its performance over both omics and clinical data types. As a performance assessment metric, we compared the number and goodness of the selected feature using SIVS to those selected by LASSO regression. The results suggested that the feature space selected by SIVS was, on average, 41\\% smaller, without having a negative effect on the model performance. A similar result was observed for comparison with Boruta and Caret RFE.The method is implemented as an R package under GNU General Public License v3.0 and is accessible via Comprehensive R Archive Network (CRAN) via https://cran.r-project.org/web/packages/sivs/index.html or through Github via https://github.com/mmahmoudian/sivs/Supplementary data are available at Bioinformatics online.}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btab501},
    url = {https://doi.org/10.1093/bioinformatics/btab501},
    note = {btab501},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btab501/39070854/btab501.pdf},
}
```
## 0.2.5

There is no functional changes in the R codebase. All changes are in the documentation to add the citation information. If you are interested in knowing the changes in details, please read the git commits in the project's Github page.

## 0.2.4

All the changes are minor and does not affect any of the numeric outcomes of the package.

- `sivs()`
    - [add] a debug mode was added to the sivs(). It can be turned on by `debug.mode = TRUE`
    - [fix] the code is now complient with R v4.x.x as it uses `inherit()` to check class for array objects
    - [fix] the importance levels of some verbosity messages is now corrected to be uniformed
    - [update] the output of the parallel clusters are now sent to `/dev/null` (or its alternative in Windows) by default. Use `debug.mode = TRUE` to send them to stdout
    - [fix] the sessionInfo is not correctly returned in the final object
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Stable Iterative Variable Selection (SIVS)}
-->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

#save(list = c("train_x", "train_y", "validation_x", "validation_y", "sivs_obj"), file = "vignette_v0.2.1.RData")
load(url("https://seafile.utu.fi/f/13e0ef294b374549b499/?dl=1"))
```

# Tutorial of basic Usage

## Description

Stable Iterative Variable Selection (SIVS) is a feature selection method that wraps machine learning methods with internal feature selection (e.g shrinkage methods such as LASSO). In this vignette we demonstrate simple usage of this method. We would like to encourage users to always stay critical of the results they get and not use SIVS as a blackbox method.

For the purpose of this simple tutorial, we use [Arcene Data Set](https://archive.ics.uci.edu/ml/datasets/Arcene) which is a sample dataset to demonstrate feature selection methods via a binary classification. As the webpage states:

> ARCENE's task is to distinguish cancer versus normal patterns from mass-spectrometric data.


## Download and read dataset

```r
train_x <- read.csv(file = "https://archive.ics.uci.edu/ml/machine-learning-databases/arcene/ARCENE/arcene_train.data", header = F, sep = " ", strip.white = T, blank.lines.skip = T, stringsAsFactors = F)
train_y <- read.csv(file = "https://archive.ics.uci.edu/ml/machine-learning-databases/arcene/ARCENE/arcene_train.labels", header = F, sep = " ", strip.white = T, blank.lines.skip = T, stringsAsFactors = F)

validation_x <- read.csv(file = "https://archive.ics.uci.edu/ml/machine-learning-databases/arcene/ARCENE/arcene_valid.data", header = F, sep = " ", strip.white = T, blank.lines.skip = T, stringsAsFactors = F)
validation_y <- read.csv(file = "https://archive.ics.uci.edu/ml/machine-learning-databases/arcene/arcene_valid.labels", header = F, sep = " ", strip.white = T, blank.lines.skip = T, stringsAsFactors = F)


train_y <- train_y[, 1]
validation_y <- validation_y[, 1]
```

## Explore, Pre-process, and clean data

Check for possible missing values in the data:

```{r}
library("varhandle")

knitr::kable(varhandle::inspect.na(d = train_x, barplot = F))
knitr::kable(varhandle::inspect.na(d = validation_x, barplot = F))
```

There is a column (V10001) that has 100% `NA` in each of the two datasets (probably because the data is not CSV and some white-space issues). Let's remove that and move on:

```{r}
train_x <- train_x[, -10001]
validation_x <- validation_x[, -10001]
```

Now let's have a better understanding of what variables we are dealing with in this R session:

```{r}
knitr::kable(varhandle::var.info(regex = "_[xy]$"))
```


## Using SIVS

Now we can run SIVS to get the most important features. Note that for the sake of making the vignette, the progressbar is set to FALSE and the verbosity has set to none so that the function does not print its progress. It's worth noting that due to the nature of this feature selection method, running a high-dimensional data in SIVS would take time, therefore, the more CPU core you allocate to it, the faster it will do the job. The following would take about 19 minutes with 2 CPU cores (`parallel.cores = 2`) and it will be much faster if you use `parallel.cores = "grace"` which will use all CPU cores except 1.

```r
library("sivs")

sivs_obj <- sivs::sivs(x = train_x,           # it should be a matrix or dataframe with features as columns and samples as rows
                       y = factor(train_y),
                       verbose = "none",      # This is for demonstration, you leave the verbose to be "general"
                       progressbar = FALSE)   # This is for demonstration, you leave the progressbar on
```

By providing factor with 2 levels, we are implying to use binomial. Alternatively, you can define `family = "binomial"` and it will be passed to internal method (which the default is glmnet).

```{r run_sivs, include = FALSE}
library("sivs")
# bypass CRAN's error on number of used cores! Apparently on CRAN you have limit your code to only 2 cores!
#chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#
#if ((nzchar(chk) && chk == "TRUE") || (!identical(Sys.getenv("NOT_CRAN", unset = "true"), "true"))) {
#    # use 2 cores in CRAN/Travis/AppVeyor
#    sivs_obj <- sivs::sivs(x = train_x, y = factor(train_y), verbose = "none", progressbar = FALSE, parallel.cores = 2)
#}else{
#    sivs_obj <- sivs::sivs(x = train_x, y = factor(train_y), verbose = "none", progressbar = FALSE)
#}
```


Having the SIVS object, we can plot the object to have better understanding of what it has found:

```{r plotting_sivs, fig.height=9, fig.width=8}
layout(mat = matrix(c(1,2,
                      3,3),
                    nrow = 2,
                    byrow = T))
{
    plot(sivs_obj)
    layout(1)
}

```

The top-left plot shows the how many times each feature was selected in the iterative step. The top-right plot shows the coefficient distribution of each of the features selected in the iterative step. This plot is sorted based on the median of the coefficient of each feature. Lastly, the bottom plot shows the result of recursive feature elimination (rfe) step. Each blue bar shows the importance of that feature, whereas each boxplot shows the distribution of areas under the receiver operating characteristic curves (AUROC) after removal of that feature. The more we move towards left-hand side of the plot the more features we have removed and as expected, the AUROCs start decreasing. The read vertical line is the suggested cutoffs based on the `sivs::suggest()` function when provided by the `sivs_obj`. Note that In this example due to the data we see only the 0.05 suggested threshold, but in datasets with more difficult setup that have more of less-important features, we would also see another suggestion as a vertical red line for 0.01 cutoff. These cutoffs are suggestions and users are encourage to treat them as such. These two cutoffs are pre-defined in the function as loose cutoffs and we encourage users to investigate the suggested features and make the cutoff more strict is possible according to their data and setup and the question they are trying to answer.

Finally, we can extract the list of suggested features by SIVS. Considering that we cannot see 0.01 line in the bottom plot of the figure above and the 0.05 is where the boxplots start dropping, let's choose a stricter threshold. The strictness can be between 0 (very loose) to 1 (very strict). Here for the sake of the example, we use 0.5 as it is the exact middle of the strictness range: 

```{r sivs_suggest}
sivs::suggest(sivs_obj, strictness = 0.5)
```

Ultimately, one can compare the performance of a machine learning method with and without SIVS. In the following example we use `glmnet` package as in this vignette it was used by SIVS as the internal method.

```{r glmnet_vs_sivs, fig.height=10, fig.width=5}
library("glmnet")
library("pROC")

# build a model without SIVS
set.seed(12345)
glmnet_model <- glmnet::cv.glmnet(x = data.matrix(train_x),
                                  y = factor(train_y),
                                  family = "binomial")

# build a model with SIVS
set.seed(12345)
sivs_glmnet_model <- glmnet::cv.glmnet(x = data.matrix(train_x[, sivs::suggest(sivs_obj, strictness = 0.5)]),
                                       y = factor(train_y),
                                       family = "binomial")

# predict both training set and validation sets
glmnet_train_pred <- predict(object = glmnet_model,
                             newx = data.matrix(train_x),
                             s = "lambda.min",
                             type = "response")
glmnet_validation_pred <- predict(object = glmnet_model,
                                  newx = data.matrix(validation_x),
                                  s = "lambda.min",
                                  type = "response")

sivs_glmnet_train_pred <- predict(object = sivs_glmnet_model,
                                  newx = data.matrix(train_x[, sivs::suggest(sivs_obj, strictness = 0.5)]),
                                  s = "lambda.min",
                                  type = "response")
sivs_glmnet_validation_pred <- predict(object = sivs_glmnet_model,
                                       newx = data.matrix(validation_x[, sivs::suggest(sivs_obj, strictness = 0.5)]),
                                       s = "lambda.min",
                                       type = "response")




glmnet_train_roc <- pROC::roc(response = factor(train_y),
                              predictor = as.numeric(glmnet_train_pred))
glmnet_validation_roc <- pROC::roc(response = factor(validation_y),
                                   predictor = as.numeric(glmnet_validation_pred))

sivs_glmnet_train_roc <- pROC::roc(response = factor(train_y),
                                   predictor = as.numeric(sivs_glmnet_train_pred))
sivs_glmnet_validation_roc <- pROC::roc(response = factor(validation_y),
                                        predictor = as.numeric(sivs_glmnet_validation_pred))

layout(mat = matrix(1:2, nrow = 2))
{
    plot(glmnet_train_roc, col = "salmon", main = "Performance on training data")
    plot(sivs_glmnet_train_roc, col = "cornflowerblue", add = T)
    legend("bottomright",
           fill = c("salmon", "cornflowerblue"),
           legend = c(paste0("glmnet (AUROC=",
                             round(pROC::auc(glmnet_train_roc),
                                   digits = 4),
                             ")"),
                      paste0("SIVS + glmnet (AUROC=",
                             round(pROC::auc(sivs_glmnet_train_roc),
                                   digits = 4),
                             ")")))
    

    plot(glmnet_validation_roc, col = "salmon", main = "Performance on validation data")
    plot(sivs_glmnet_validation_roc, col = "cornflowerblue", add = T)
    legend("bottomright",
           fill = c("salmon", "cornflowerblue"),
           legend = c(paste0("glmnet (AUROC=",
                             round(pROC::auc(glmnet_validation_roc),
                                   digits = 4),
                             ")"),
                      paste0("SIVS + glmnet (AUROC=",
                             round(pROC::auc(sivs_glmnet_validation_roc),
                                   digits = 4),
                             ")")))
    layout(1)
}


```

As we can see, the model build without SIVS has AUROC of `r round(pROC::auc(glmnet_validation_roc)*100, digits = 2)`% on validation set where as the model built based on the features selected by SIVS has AUROC of `r round(pROC::auc(sivs_glmnet_validation_roc)*100, digits = 2)`%. Now let's compare the number of selected features:


```{r compare_feature_count}
# extract the coefficients
sum(coef(glmnet_model) != 0)
sum(coef(sivs_glmnet_model) != 0)
```
This shows that without SIVS we would have selected `r sum(coef(glmnet_model) != 0)` features, but with SIVS we have reduced it to just `r sum(coef(sivs_glmnet_model) != 0)`.

