<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/293626039.svg)](https://zenodo.org/badge/latestdoi/293626039)
[![R-CMD-check](https://github.com/azizka/IUCNN/workflows/R-CMD-check/badge.svg)](https://github.com/azizka/IUCNN/actions)
<!-- badges: end -->


# IUCNN
Batch estimation of species' IUCN Red List threat status using neural networks.

# Installation
1. Install IUCNN directly from Github using devtools. 
```r
install.packages("devtools")
library(devtools)

install_github("IUCNN/IUCNN")
```

2. Since some of IUCNNs functions are run in Python, IUCNN needs to set up a Python environment. This is easily done from within R, using the `install_miniconda()` function of the package `reticulate` (this will need c. 3 GB disk space).
If problems occur at this step, check the excellent [documentation of reticulate](https://rstudio.github.io/reticulate/index.html).
```r
install.packages("reticulate")
library(reticulate)
install_miniconda()
```


3. Install the tensorflow python library. If you are using **MacOS** or **Linux** it is recommended to install tensorflow using conda:
```r
reticulate::conda_install("r-reticulate","tensorflow=2.4")
```

If you are using **Windows**, you can install tensorflow using pip:

```r
reticulate::py_install("tensorflow~=2.4.0rc4", pip = TRUE)
```

4. Finally install the npBNN python library from Github:

```r
reticulate::py_install("https://github.com/dsilvestro/npBNN/archive/refs/tags/v0.1.11.tar.gz", pip = TRUE)
```

# Usage
There are multiple models and features available in IUCNN. A vignette with a detailed tutorial on how to use those is available as part of the package: `vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")`. Running IUCNN will write files to your working directory.

A simple example run for terrestrial orchids (This will take about 5 minutes and download ~500MB of data for feature preparation into the working directory):

```r
library(tidyverse)
library(IUCNN)

#load example data 
data("training_occ") #geographic occurrences of species with IUCN assessment
data("training_labels")# the corresponding IUCN assessments
data("prediction_occ") #occurrences from Not Evaluated species to prdict

# 1. Feature and label preparation
features <- iucnn_prepare_features(training_occ) # Training features
labels_train <- iucnn_prepare_labels(x = training_labels,
                                     y = features) # Training labels
features_predict <- iucnn_prepare_features(prediction_occ) # Prediction features

# 2. Model training
m1 <- iucnn_train_model(x = features, lab = labels_train)

summary(m1)
plot(m1)

# 3. Prediction
iucnn_predict_status(x = features_predict,
                     model = m1)
```

With model testing

```r
library(tidyverse)
library(IUCNN)

#load example data 
data("training_occ") #geographic occurrences of species with IUCN assessment
data("training_labels")# the corresponding IUCN assessments
data("prediction_occ") #occurrences from Not Evaluated species to predict

# Feature and label preparation
features <- iucnn_prepare_features(training_occ) # Training features
labels_train <- iucnn_prepare_labels(x = training_labels,
                                     y = features) # Training labels
features_predict <- iucnn_prepare_features(prediction_occ) # Prediction features


# Model testing
# For illustration models differing in dropout rate and number of layers

mod_test <- iucnn_modeltest(x = features,
                            lab = labels_train,
                            logfile = "model_testing_results-2.txt",
                            model_outpath = "iucnn_modeltest-2",
                            mode = "nn-class",
                            dropout_rate = c(0.0, 0.1, 0.3),
                            n_layers = c("30", "40_20", "50_30_10"),
                            cv_fold = 5,
                            init_logfile = TRUE)

# Select best model
m_best <- iucnn_best_model(x = mod_test,
                          criterion = "val_acc",
                          require_dropout = TRUE)

# Inspect model structure and performance
summary(m_best)
plot(m_best)

# Train the best model on all training data for prediction
m_prod <- iucnn_train_model(x = features,
                            lab = labels_train,
                            production_model = m_best)

# Predict RL categories for target species
pred <- iucnn_predict_status(x = features_predict,
                             model = m_prod)
plot(pred)

```

Using a convolutional neural network

```r
features <- iucnn_cnn_features(training_occ) # Training features
labels_train <- iucnn_prepare_labels(x = training_labels,
                                     y = features) # Training labels
features_predict <- iucnn_cnn_features(prediction_occ) # Prediction features

```

# Citation
```r
library(IUCNN)
citation("IUCNN")
```

Zizka A, Andermann T, Silvestro D (2021). "IUCNN - deep learning approaches to approximate species’ extinction risk." [Diversity and Distributions, doi: 10.1111/ddi.13450](https://doi.org/10.1111/ddi.13450). 

Zizka A, Silvestro D, Vitt P, Knight T (2020). “Automated conservation assessment of the orchid family with deep
learning.” [Conservation Biology, doi: doi.org/10.1111/cobi.13616](https://doi.org/doi.org/10.1111/cobi.13616)
# IUCNN 2.0.1 (03.01.2022)
=========================
* fixed bug with the export of the iucnn_feature_importance function

# IUCNN 2.0.0 (16.08.2021)
=========================
* moved to IUCNN/IUCNN
* added iucnn_cnn_features function
* standardized function naming scheme
* added test for polygon validity to iucnn_biome_features
* add iucnn_bias_features function and sampbias as suggested package

# IUCNN 1.0.1(01.06.2021)
=========================
* minor fix for outputting val-loss and test-loss for nn-reg mode

# IUCNN 1.0.0(27.05.2021)
=========================
* final version for first release
* some minor fixes of typos in the vignette and readme

# IUCNN 0.9.9 (26.05.2021)
=========================
* final version for last release tests

# IUCNN 0.9.3
=========================

* add footprint_features function
* added option for selected variables and ranges to the clim_features function
* added option to remove empty biomes to the biom_features function
* improved spell-checking
* add rescale option for eoo and aoo as part of geo_features
* set EOO to AOO for species with less than three occurrences
* updated the normalization of climate variables

# IUCNN 0.9.2

* bug fix predict_iucn function
* bug fix train_iucn function

# IUCNN 0.9.1

* Added description information
* Add citation information
* Start Readme
* Add vignette

# IUCNN 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* package skeleton
---
title: Approximate Red List assessments with IUCNN
output: rmarkdown::html_vignette
vignette: >
 %\VignetteIndexEntry{Approximate Red List assessments with IUCNN}
 %\VignetteEngine{knitr::knitr}
 %\VignetteEncoding{UTF-8}
---



# Background
The Red List of the International Union for the Conservation of nature (www.iucn.org, IUCN RL), is arguably one of the most thorough and widely used tools to assess the global extinction risk of species. However, the IUCN RL assessment process---usually performed by a group of specialists for each taxonomic group, or professional assessors---are time consuming, and therefore only a small fraction of global biodiversity has been evaluated for the IUCN RL, with a strong bias towards vertebrates and certain regions. These biases and the low fraction of species evaluated, bias conservation towards evaluated groups and prevent synthetic, large-scale ecological and biogeographic analyses of extinction risk. 

IUCNN implements neural networks to predict the IUCN status of so far not evaluated or data deficient species based on species traits. IUCNN models are trained on the existing IUCN RL assessments and any traits may be used for prediction, although IUCNN implements a workflow based solely on publicly available geo-referenced species occurrence records and environmental data. Typical examples for the application of IUCNN are to predict the conservation status of a large number of species, to approximate extinction risk or number of threatened species in a region or specific taxonomic group for synthetic analyses or to predict the IUCN category of individual species of interest for systematic or ecological case studies. 


```r
library(IUCNN)
library(magrittr)
library(dplyr)
```

# Installation
IUCNN uses R and python. All software needed can be installed via R.

1. install IUCNN directly from Github using devtools. 

```r
install.packages("devtools")
library(devtools)
library(IUCNN)
```

2. Python needs to be installed, for instance using miniconda and reticulated from within R (this will need c. 3 GB disk space).
If problems occur at this step, check the excellent [documentation of reticulate](https://rstudio.github.io/reticulate/index.html).

```r
install.packages(reticulate)
library("reticulate")
install_miniconda()
```

If python has been installed before, you can specify the python version to sue with `reticulate::use_python()`

3. Install the tensorflow Python module. IUCNN uses functions of the python modules tensorflow and npBNN which also need to be installed (via R). 

```r
reticulate::conda_install("r-reticulate","tensorflow=2.4")
reticulate::py_install("https://github.com/dsilvestro/npBNN/archive/v0.1.10.tar.gz", pip = TRUE)
```

# Prepare input data
IUCNN predicts the IUCN RL categories of Not Evaluated and Data Deficient species based on geographic occurrence records and a set of training species for which occurrence records and IUCN assessments that are available for a set of reference species (training data). The amount of training species necessary varies with the number of categories but in general "the more, the better". Ideally, the training dataset should comprise several hundred species or more, so a typical scenario will be to use all available plant species from a region, or all available species from a plant family. If the availability of training species is limited, a good option can be to reduce detail and predict Possibly threatened (IUCN categories "CR", "EN", and "VU") v. Not threatened species ("NT" and "LC").

Three types of input are necessary, which are easily available for many species: 

## 1. Geographic occurrence records of training species (training occurrences)
Occurrence records might be obtained from a variety of databases, For example, from field collections or public databases such BIEN (https://bien.nceas.ucsb.edu/bien/) or GBIF (www.gbif.org). GBIF data can be obtained from within R via the rgbif package, See [here](https://docs.ropensci.org/rgbif/articles/index.html) for a tutorial on how to do so. IUCNN needs a dataset with (at least) three columns, containing the species name, decimal longitude coordinates and decimal latitude coordinates. If you are interested in cleaning records from GBIF, you may want to have a look at this [blog post](https://data-blog.gbif.org/post/gbif-filtering-guide/) and check out the [CoordinateCleaner](https://github.com/ropensci/CoordinateCleaner) and [bRacatus](https://github.com/EduardoArle/bRacatus) packages. 

## 2. IUCN Global Red List assessment of the training species (training labels)
IUCN RL assessments for the training species can be obtained from IUCN, either via www.iucnredlist.org or via the rredlist package via R (preferred for many species). See [here](https://ropensci.org/tutorials/rredlist_tutorial/) for a tutorial on how to use rredlist. It is important, that all target label classes are well represented in the training data, which is rarely the case for IUCN data, since for instance "VU" and "NT" is rare. If the classes are to imbalanced, consider using possibly threatened (IUCN categories "CR", "EN", and "VU") v. not threatened species ("NT" and "LC"), or the supersampling option of the `iucnn_train_model` function.

## 3. Geographic occurrence records of the species for which the IUCN status should be predicted (predict occurrences)
Geographic occurrence for the target species, in the same format as for the training occurrences described above.

Example dataset are available with IUCNN: `data(training_occ)` (training occurrences), `data(training_labels)` (training labels) and `data(prediction_occ)`.

## Feature preparation
IUCNN uses per species traits to as features for the neural networks. The required input format is a data.frame with one row per species , one column containing the species name and any number of additional columns containing the numerical features for each species. In general, features might represent any trait, for instance from taxonomy (e.g., family), anatomy (e.g., body size), ecology (e.g., feeding guild) or conservation (e.g., population dynamics). However, since often only geographic occurrence data are available IUCNN contains functions to obtain default features from geo-referenced occurrence records alone, by combining them with publicly available environmental data. These default features informing on species range, climatic niche, human footprint and biomes. See Table 2 for a detailed list of all default features. Users may chose to use specific groups of features only via the `type` option of `iucnn_prepare_labels`. In this tutorial, we will use the example datasets from the Orchid family (Orchidaceae) provided with the IUCNN package, 

You can prepare the default features with a single call to `iucnn_prepare_features`

```r
data("training_occ") #geographic occurrences of species with IUCN assessment
data("prediction_occ")

features_train <- iucnn_prepare_features(training_occ) # Training features
features_predict <- iucnn_prepare_features(prediction_occ) # Prediction features
```

## Label preparation
IUCNN expects the labels for training as numerical categories. So, to use IUCN RL categories, those need to be converted to numeric in the right way. This can be done using the `iucnn_prepare_labels` function. The function converts the category labels as obtained from the IUCN RL into standardized numeric values, either on the detailed level of IUCN RL categories or the broader Possibly threatened/Not threatened level. See `?iucnn_prepare_labels` for more information. The labels will be converted into numeric categories following the `accepted_labels` argument, so for instance, in the default case: LC -> 0 and CR -> 4. If you change the accepted labels, the match will change accordingly.


```r
data("training_labels")

labels_train <- iucnn_prepare_labels(x = training_labels,
                                     y = features_train) # Training labels
```

# Running IUCNN
Running IUCNN consists of two steps: 1) training a neural network and 2) predicting the status of new species. IUCNN contains three different neural network approaches to predict the IUCN status of species, which can all be customized. We present the default approach here, see section "Customizing analyses" of this tutorial for details on how to train a Bayesian or regression type neural network. 

## Model training
Based on the training features and labels, IUCNN trains a neural network, via the `iucnn_train_model` function. There are multiple options to customize the design of the network, including among others the number of layers and the fraction of records used for testing and validation. The `iucnn_train_model` function will write a folder to the working directory containing the model and return summary statistics including cross-entropy loss and accuracy for the validation set, which can be used to compare the performance of different models.

The following code trains a neural network model with 3 hidden layers of 60, 60, and 20 nodes, with ReLU activation function. By specifying a seed (here, the default: 1234) we make sure the same subsets of data are designated as training, validation and test sets across different runs and model configurations (see below). The model with estimated weights will be saved in the current working directory. 


```r
res_1 <- iucnn_train_model(x = features_train,
                           lab = labels_train, 
                           path_to_output = "iucnn_model_1")
```

The `summary` and `plot` methods give an overview on the training process and model performance. 


```r
summary(res_1)
#> A model of type nn-class, trained on 702 species and 45 features.
#> 
#> Training accuracy: 0.6
#> Validation accuracy: 0.532
#> Accuracy on unseen data (test set): 0.549
#> 
#> Label detail: 5 Classes (detailed)
#> 
#> Label representation
#>   Label Input_count Estimated_count
#> 1     0          67              92
#> 2     1          12               0
#> 3     2          28               3
#> 4     3          51              80
#> 
#> Confusion matrix (rows test data and columns predicted):
#>    LC NT VU EN CR
#> LC 56  0  1 10  0
#> NT  9  0  0  3  0
#> VU 15  0  1 12  0
#> EN 11  0  1 39  0
#> CR  1  0  0 16  0
plot(res_1)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

## Predict IUCN Global Red List status
The trained model can then predict the conservation status of *Not Evaluated* and *Data Deficient* species with the `iucnn_predict_status` function. The output contains a data frame with species names and numeric labels (as generated with iucnn_prepare_labels).


```r
predictions <- iucnn_predict_status(x = features_predict, 
                                    model = res_1)

plot(predictions)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

It is important to remember the following when using IUCNN:

1. The resulting IUCNN categories are approximations only. While IUCNN has reached accuracies between 80 and 90% on the broad (threatened v non-threatened) level and up to 80% on the detailed level in some cases, the accuracy may be considerably lower in other cases, which means that some species will be mis-classified.

2. IUCNN is indifferent to the provided features. On the one hand this means that any species traits for which data is available can be used, but on the other hand this means that thought is needed in the choice of the features. The default features of IUCNN are usually a safe choice. The number of features is not limited, but currently IUCNN does not support missing values in the feature table and removes species with missing values. 

3. IUCNN is indifferent to the relation between training and test data. So it is possible to use training data from Palearctic birds to predict the conservation status of South American Nematodes. This is not recommended. Instead, a better approach will be to predict the conservation status of species, from training data of the same genus, order, or family. Alternatively, training data could be chosen on geographic region or functional aspects (e.g., feeding guild or body size). However some inclusion of taxonomy/evolutionary history for the choice of training data is recommended.

4. The amount of training data is important. The more the better. Minimum several hundred training species with a more or less equal distribution on the label classes should be included. If training data is limited, the broader Threatened/Not threatened level is recommended. 

5. If the proportion of the IUCN RL categories is imbalanced in the training data, the neural networks may be biased towards reproducing these frequencies in the prediction, especially if the imbalance of categories or the difference in category frequencies among training and prediction set are large. To avoid this category frequencies should be balanced in the training data if possible. Otherwise the use of the `supersampling` or option a `nn-reg` type model, or a limitation to the broader Possibly threatened/Not threatened level of detail may remedy the issue. 

6. IUCNN predictions are not equivalent to full IUCN Red List assessments. We see the main purpose of IUCNN in 1) identifying species that will likely need conservation action to trigger a full IUCN assessment, and 2) provide large-scale overviews on the extinction risk in a given taxonomic group, for instance in a macro-ecological and macro-evolutionary context.

# Customizing IUCNN analyses
IUCNN contains multiple options to customize the steps of the analyses to adapt the fully connected neural networks to the peculiarities of IUCN RL and species distribution data. Below we describe the most important options to customize 1) feature and label preparation, 2) model training and testing, and 3) status prediction. The most important steps and options to customize an IUCNN analysis are summarized in Table 1.

Table 1. Critical steps to customize an IUCNN analysis and relevant considerations at each point.

| Step | Function(s) | Argument | Description | User consideration |
| ---------------------------------------------------------------------------------------- | ----------------------------------- | ---------------------- |---------------- | ------------ |
| Feature design | iucnn\_prepare\_features, iucnn\_feature\_importance | \- | Defines the features to be extracted from the provided occurrence records. Available defaults are: biome presence, bioclim variables, human footprint and geographic features. Averaged per species and rescaled. | What determines extinction risk for target group? Which data are available? What do the results of feature importance suggest? |
| Label detail | prep\_labels | level, threatened | Into how many different categories should the species be classified. Can be any number, defaults support full IUCN categories (LC, NT, VU, EN, CR) or binary (Possibly threatened v. Not threatened) | Which detail is needed? Is the accuracy of the detailed level sufficient for the target application? |
| Model type | iucnn\_train\_model | mode | Which model framework should be applied: a categorical classification, a classification taking the ordinal number of categories into account, or a classification based in a Bayesian framework | How many target categories are there? How important is a high accuracy for intermediate classes (e.g. VU, NT, EN)? How important is the uncertainty estimation for each species? |
| Model structure | iucnn\_train\_model | validation fraction | The fraction of the input training data used for validation (v. training). | Which fraction of the data should be used for validation (and hence not training)? How large is the training data? |
| Model structure | iucnn\_train\_model | cv\_fold | The number of folds used for cross-validation. For instance if = 5, the data is divided into 5 folds with 20% of the data used for validation in each run. | How large is the training data?| At a given number of folds, will the subsets still include all label classes? |
| Model structure | iucnn\_train\_model | n\_layers | The number of hidden layers and nodes in the neural network. | How complex should the model be? |
| Model structure | iucnn\_train\_model | balance\_classes | Should the frequency of the class labels in the training data be balanced using supersampling? | How imbalanced are the class labels? Will the frequency of class labels differ between training data and prediction data set? |
| Model structure | iucnn\_train\_model | act\_f\_out/act\_f | The activation function of the neural network | Which relationship between features and labels is expected? |
| Model structure | iucnn\_train\_model | label\_stretch\_factor | The factor to stretch input class labels to | Am I using a nn-reg type model?  Does model testing suggest an effect on model accuracy? |
| Model structure | iucnn\_train\_model | drop\_out rate | The number of nodes to be removed in individual epochs. Necessary if a target threshold is to be used with a nn-class model (not bnn-class though) | Is a target accuracy to be sued for prediction? How many nodes are there in the model? |
| Prediction | predict\_iucnn | target\_acc | Defines an overall target accuracy for the model Species which cannot be classified with enough certainty to reach this threshold are labels as data deficient. | Which error rate is acceptable? Which proportion of species needs to be included? |


## 1) Features and Labels
### Add and remove feature blocks
The default features are selected based on empirical tests on relevance for different taxa and regions. However, for some analyses only part of the features may be relevant. You can exclude feature blocks using the `type` argument of the `iucnn_prepare_features` function. For instance, to exclude the biome features:


```r
features_train2 <- iucnn_prepare_features(training_occ, 
                                          type = c("geographic", 
                                                   "climate", 
                                                   "humanfootprint"))
```

### Prepare features individually
If more control over feature preparation is necessary, each feature block can be obtained by an individual function.

Table 2. Functions to obtain default features and options to customize the features.

|Feature block|Function name|Options to customize|
|---|---|---|
|Geographic|`iucnn_geographic_features`|-|
|Biomes|`iucnn_biome_features`|change the reference dataset of biomes (biome_input, biome.id), remove biomes without any species occurrence (remove_zeros)|
|Climate|`iucnn_climate_features`|the amount of bioclim variables from the default source to be included (type), the resolution of the default input data (res)|
|Human footprint|`iucnn_footprint_features`|chose the time points from the default source (year), the break points for the different footprint categories (breaks, by default approximately quantiles on the global footprint dataset) or a default source for human footprint (footp_input)|
|Geographic bias|`iucnn_bias_features`|The resolution of the bias raster (res) or providing a template raster (ras)|

For instance:


```r
clim_features <- iucnn_climate_features(x = training_occ, 
                                        type = "selected")

clim_features2 <- iucnn_climate_features(x = training_occ, 
                                         type = "all")
```

### Use custom features
It is also possible to provide features unrelated to the default features. They may contain any continuous or categorical features, but some processing will be needed. The format needs to be a data.frame with a compulsory column containing the species name. Continuous variables should be rescaled to cover a similar range, whereas categorical features should be coded binary (present/absent, as the results of `iucnn_biome_features`).

For instance:


```r
feat <- data.frame(species = c("Adansonia digitata", "Ceiba pentandra"),
 max_plant_size_m = c(25, 50),
 africa = c(1,1),
 south_america = c(0,1),
 fraction_of_records_in_protected_area = c(25, 75))
```

Table 3. Description of the default features included in `iucnn_prepare_features`. All continuous variables are rescaled to a similar range.

| Feature | Block | Name | Description |
|---|---|---|---|
|tot_occ|Geographic|Number of occurrences|The total number of occurrences available for this species|
|uni_occ|Geographic|Number of geographically unique occurrences|The number of geographically unique records available for this species|
|mean_lat|Geographic|Mean latitude|The mean latitude of all records of this species|
|mean_lon|Geographic|Mean longitude|The mean longitude of all records of this species|
|lat_range|Geographic|Latitudinal range|The latitudinal range (.95 quantile - .05 quantile).|
|lon_range|Geographic|Longitudinal range|The longitudinal range (.95 quantile - .05 quantile).|
|alt_hemisphere|Geographic|The hemisphere|0 = Southern hemisphere, 1 = Northern hemisphere|
|eoo|Geographic|Extend of Occurrence|The extend of occurrence. Calculated by rCAT. For species with less than 3 records set to AOO|
|aoo|Geographic|Area of Occupancy|The area of occupancy, as the sum of area of 4sqkm grid cells, where the species occurs|
|1|Biome|Tropical & Subtropical Moist Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|2|Biome|Tropical & Subtropical Dry Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|3|Biome|Tropical & Subtropical Coniferous Forests|Are at least 5% of the species records present in this biome?|
|4|Biome|Temperate Broadleaf & Mixed Forests|Are at least 5% of the species records present in this biome?|
|5|Biome|Temperate Conifer Forests|Are at least 5% of the species records present in this biome?|
|6|Biome|Boreal Forests/Taiga|Are at least 5% of the species records present in this biome?|
|7|Biome|Tropical & Subtropical Grasslands, Savannas & Shrublands|Are at least 5% of the species records present in this biome?|
|8|Biome|Temperate Grasslands, Savannas & Shrublands|Are at least 5% of the species records present in this biome?|
|9|Biome|Flooded Grasslands & Savannas|Are at least 5% of the species records present in this biome?|
|10|Biome|Montane Grasslands & Shrublands|Are at least 5% of the species records present in this biome?|
|11|Biome|Tundra|Are at least 5% of the species records present in this biome?|
|12|Biome|Mediterranean Forests, Woodlands & Scrub|Are at least 5% of the species records present in this biome?|
|13|Biome|Deserts & Xeric Shrublands|Are at least 5% of the species records present in this biome?|
|14|Biome|Mangroves|Tropical & Subtropical Moist Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|98|Biome|Lake|Are at least 5% of the species records present in this biome?|
|99|Biome|Rock and ice|Are at least 5% of the species records present in this biome?|
|bio1|Climate|Annual Mean Temperature|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio4|Climate|Temperature Seasonality|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio11|Climate|Mean Temperature of Coldest Quarter|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio12|Climate|Annual Precipitation|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio15|Climate|Precipitation Seasonality|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio17|Climate|Precipitation of Driest Quarter|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|range_bio1|Climate|Range of annual Mean Temperature|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio4|Climate|Range of temperature Seasonality|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio11|Climate|Range of mean Temperature of Coldest Quarter|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio12|Climate|Range of annual Precipitation|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio15|Climate|Range of precipitation Seasonality|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio17|Climate|Range of precipitation of Driest Quarter|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|humanfootprint_1993_1|Human footprint | Human footprint year 1993 lowest impact|The fraction of records in areas of the lowest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_2|Human footprint|Human footprint year 1993 intermediate impact 1|The fraction of records in areas of the second lowest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_3|Human footprint|Human footprint year 1993 intermediate impact 2|The fraction of records in areas of the second highest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_4|Human footprint|Human footprint year 1993 highest impact|The fraction of records in areas of the highest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_1|Human footprint|Human footprint year 2009 lowest impact|The fraction of records in areas of the lowest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_2|Human footprint|Human footprint year 2009 intermediate impact 1|The fraction of records in areas of the second lowest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_3|Human footprint|Human footprint year 2009 intermediate impact 2|The fraction of records in areas of the second highest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_4|Human footprint|Human footprint year 2009 highest impact|The fraction of records in areas of the highest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|

### Labels: Full categories vs Threatened/Not threatened
The `iucnn_prepare_labels` function may accepted any custom labels as long as they are included in the `accepted_labels` option. It also can provide a classification into threatened/non-threatened, via the `level` and `threatened` options. On the broader level the model accuracy is usually significantly higher.

For instance:


```r
labels_train <- iucnn_prepare_labels(training_labels,
                                     y = features, 
                                     level = "broad")
```

## 2) Model training - NN regression model
### Customizing model parameters
The `iucnn_train_model` function contains various options to customize the neural network, including among other the fraction of validation and test data, the maximum number of epochs, the number of layers and nodes, the activation function , dropout and randomization of the input data. See `?iucnn_train_model` for a comprehensive list of options and their description. By default, `iucnn_train_model` trains a neural network with three hidden layers with 50, 30 and 10 nodes and a sigmoid as activation function. Depending on your dataset different networks may improve performance. For instance, you can set up a different model with 1 hidden layer of 60 nodes, a sigmoid activation function and without using a bias node in the first hidden layer.


```r
res_2 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 dropout_rate = 0.3,
 path_to_output= "iucnn_model_2",
 n_layers = "60",
 use_bias = FALSE,
 act_f = "sigmoid")
```

You can compare the validation loss of the models using `res_1$validation_loss` and `res_2$validation_loss`. Model 2 in this case yields a lower validation loss and is therefore preferred. Once you chose the preferred model configuration based on validation loss, we can check test accuracy of best model: `res_2$test_accuracy`. The `iucnn_train_model` function contains various options to adapt the model. See `?iucnn_train_model` for more detail. 

### Changing the modeling algorithm
There are three neural network algorithms implemented in IUCNN. Besides the default classifier approach based on a tensorflow implementation, these are a Bayesian neural network classifier and a regression type neural network.

The Bayesian approach has the advantage that it returns true probabilities for the classification of species into the relative output classes (e.g. 80% probability of a species to be LC). We consider this approach more suitable for classification of species into IUCN categories, than the default option. It will need more time for model training and should best be applied once you have identified the best model parameters using the default approach. You can run a BNN setting the `mode` option of `iucnn_train_model` to `"bnn-class"`.


```r
res_3 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 path_to_output = "iucnn_model_3",
 mode = 'bnn-class')
```

IUCNN also offers the option to train a NN regression model instead of a classifier. Since the IUCN threat statuses constitute a list of ordinal categories sorted by increasing threat level, we can model the task of estimating these categories as a regression problem. Such a model can be trained with the `iucnn_train_model()` function, specifying to train a regression model by setting `mode = 'nn-reg'`.


```r
res_4 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 path_to_output = "iucnn_model_4",
 mode = 'nn-reg',
 rescale_features = TRUE)
```

### Feature importance
The `feature_importance` function can be used to gauge the importance of different feature blocks or individual features for model performance. The function implements the permutation feature importance technique, which randomly reshuffles the values within individual features or blocks of features and evaluate how this randomization affects the models prediction accuracy. If a given feature (or block of features) is important for the models ability to predict, randomizing this feature will lead to a large drop in prediction accuracy. When using `feature_importance` with features other than the default, feature blocks can be defined using the `feature_blocks` option.

```r
fi <- feature_importance(x = res_1)
plot(fi)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

### Model testing
Before training the final model used for predicting the conservation status of not evaluated species, it is recommended to use the `iucnn_modeltest` function for finding the best settings for your model and dataset. This process, often referred to as hyperparameter tuning, is an essential step for building the most suitable model for the prediction task. The `iucnn_modeltest` function allows you to provide any settings for `iucnn_train_model` as vectors, which will lead the function to train a separate model for each provided setting. The function will explore all possible permutations of the provided settings, so that the following command results in 9 different models being trained:


```r
modeltest_results <- iucnn_modeltest(features,
 labels,
 dropout_rate = c(0.0,0.1,0.3),
 n_layers = c('30','40_20','50_30_10'))
```

The model specifications and settings of each tested model are written to a log-file and can be inspected with the `iucnn_best_model` function, to decide which model settings to pick as best model. Different criteria for picking the best model can be selected, such as best prediction accuracy, best predicted over-all status distribution, lowest weighted mis-classification error, etc.


```r
best_m <- iucnn_best_model(modeltest_results, criterion='val_acc')
```

After model testing, it is necessary to retrain using the model specifications identified by `iucnn_best_model`. It is possible to take the respective settings directly from the output of the `iucnn_best_model` function, via the `production-model` argument of `iucnn_train_model`.


```r
# Train the best model on all training data for prediction
m_prod <- iucnn_train_model(train_feat,
                      train_lab,
                      production_model = m_best,
                      overwrite = TRUE)
```

The production model can then be used for the final predictions.


```r
# Predict RL categories for target species
pred <- iucnn_predict_status(pred_feat,
                      m_prod)
plot(pred)
```


## 3) Status prediction
The `iucnn_predict_status` function offers options to customize the predictions. The most important option in many cases is `target_acc`, which allows to set an overall target-accuracy threshold that the model needs to achieve. This option is only available for nn-class and nn-reg models that were trained using dropout (see help function of `iucnn_train_model` for more explanation), as well as for all bnn-class models. The `target_acc` will be achieved by the model being more selective with making a category call for a given instance. All species that cannot be classified with enough certainty to reach this target accuracy will be classified as NA (Not Assessed).

```r
pred_2 <- iucnn_predict_status(x = features_predict, 
 target_acc = 0.7,
 model = res_2)
plot(pred_2)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

Furthermore, you can turn off the `return_IUCN` option if to return the numerical labels instead of the IUCNN RL category labels.

```r
pred_3 <- iucnn_predict_status(x = features_predict, 
 model = res_2,
 return_IUCN = FALSE)
```

The output of the `iucnn_predict_status` function is an "iucnn_predictions" object, that contains several output objects. The predicted labels of the individual instances are accessible with `pred_2$class_predictions` and label probabilities estimated by the neural network via `pred_2$mc_dropout_probs` for "nn-class" and "nn-reg" with dropout, or `pred_2$posterior_probs` for "bnn-class". For more detail, the `pred_2$raw_predictions` object contains the individual label probabilities resulting from the softmax output layer in case of "nn-class", or the regressed labels in case of "nn-reg".

## 4) The number of species per category
Another statistic that can be extracted from the "iucnn_predictions" object is the overall category distribution predicted for the given prediction instances. This can be accessed with `pred_2$pred_cat_count`, which shows the distribution of the label predictions, including the count of species that could not be predicted given the chosen `target_acc`.

Another statistic are the `pred2$sampled_cat_freqs` (only available for dropout models and all "bnn-class" models, see above), which show the class distribution as sampled from the `pred_2$mc_dropout_probs` or the `pred_2$posterior_probs` (for "bnn-class" models). The difference between `pred2$sampled_cat_freqs` and `pred_2$mc_dropout_probs`/`pred_2$posterior_probs` is that the former represents the counts of the best labels determined for each instance, whereas the latter represents labels sampled from the predicted label probabilities, which also proportionally samples the labels for a given instance that do receive the maximum label probability. The latter is done repeatedly to include the stochasticity of the random sampling of classes from the given probability vectors. The `pred_2$mc_dropout_probs`/`pred_2$posterior_probs` can be used to plot histograms of the estimates for each class, and can be reported as uncertainty intervals around the number of species in each class for the set of species that were predicted.


# Training and prediction using a convolutional neural network
Instead of the fully connected neural networks presented above, IUCNN also implements a prediction algorithm using Convolutional Neural Networks (CNNs). The input data for prediction are then rasterized per-species grids of occurrence numbers. Since the input features and network structure is different CNNs are implemented in IUCNN with a separate workflow.

## 1. Feature preparation
Features are prepared using the `iucnn_cnn_features` function which will count the number of occurrence records in a custom raster with the same extent as the species occurrences. The function can be simply run on a data frame of training and prediction occurrences separately. Yet, since the extent of these two datasets is likely different, it is recommendable to create a custom raster beforehand. Here, we will use a simple lat/lon raster, but users may provide coordinates and a raster in any suitable coordinate reference system, as long as they agree between occurrence coordinates and raster.


```r
# preapre custom raster, you can split this step if training and test occurrences have the same extent
library(terra)
data("training_occ")
data("prediction_occ")

# find the minimum latitude and longitude values for the extent of the raster
min_lon <- min(c(min(training_occ$decimallongitude), 
                 min(prediction_occ$decimallongitude)))
max_lon <- max(c(max(training_occ$decimallongitude), 
                 max(prediction_occ$decimallongitude)))
min_lat <- min(c(min(training_occ$decimallatitude), 
                 min(prediction_occ$decimallatitude)))
max_lat <- max(c(max(training_occ$decimallatitude), 
                 max(prediction_occ$decimallatitude)))
## set the coordinate reference system
ras <- rast(crs = "+proj=longlat +datum=WGS84")
## set raster extent
ext(ras) <- c(min_lon,
              max_lon, 
              min_lat, 
              max_lat)
## set raster resolution
res(ras) <- 1 # the resolution in CRS units, in this case degrees lat/lon


 # Training features
cnn_features <- iucnn_cnn_features(x = training_occ, 
                                   y = ras)

# Prediction features
cnn_features_predict <- iucnn_cnn_features(x = prediction_occ,
                                           y = ras)
 # Training labels
cnn_labels <- iucnn_prepare_labels(x = training_labels,
                                   y = cnn_features)
```


## 2. Model training

```r
trained_model <- iucnn_cnn_train(cnn_features,
                                cnn_labels,
                                overwrite = TRUE,
                                dropout_rate = 0.1,
                                optimize_for = 'accuracy')

plot(trained_model)
summary(trained_model)
```

## 3. Status prediction

```r
pred <- iucnn_predict_status(cnn_features_predict,
                             trained_model,
                             target_acc = 0.0
                             )
```

---
title: Approximate Red List assessments with IUCNN
#output: rmarkdown::html_vignette
output: pdf_document
vignette: >
 %\VignetteIndexEntry{Approximate Red List assessments with IUCNN}
 %\VignetteEngine{knitr::knitr}
 %\VignetteEncoding{UTF-8}
---



# Background
The Red List of the International Union for the Conservation of nature (www.iucn.org, IUCN RL), is arguably one of the most thorough and widely used tools to assess the global extinction risk of species. However, the IUCN RL assessment process---usually performed by a group of specialists for each taxonomic group, or professional assessors---are time consuming, and therefore only a small fraction of global biodiversity has been evaluated for the IUCN RL, with a strong bias towards vertebrates and certain regions. These biases and the low fraction of species evaluated, bias conservation towards evaluated groups and prevent synthetic, large-scale ecological and biogeographic analyses of extinction risk. 

IUCNN implements neural networks to predict the IUCN status of so far not evaluated or data deficient species based on species traits. IUCNN models are trained on the existing IUCN RL assessments and any traits may be used for prediction, although IUCNN implements a workflow based solely on publicly available geo-referenced species occurrence records and environmental data. Typical examples for the application of IUCNN are to predict the conservation status of a large number of species, to approximate extinction risk or number of threatened species in a region or specific taxonomic group for synthetic analyses or to predict the IUCN category of individual species of interest for systematic or ecological case studies. 


```r
library(IUCNN)
library(magrittr)
library(dplyr)
```

# Installation
IUCNN uses R and python. All software needed can be installed via R.

1. install IUCNN directly from Github using devtools. 

```r
install.packages("devtools")
library(devtools)
library(IUCNN)
```

2. Python needs to be installed, for instance using miniconda and reticulated from within R (this will need c. 3 GB disk space).
If problems occur at this step, check the excellent [documentation of reticulate](https://rstudio.github.io/reticulate/index.html).

```r
install.packages(reticulate)
library("reticulate")
install_miniconda()
```

If python has been installed before, you can specify the python version to sue with `reticulate::use_python()`

3. Install the tensorflow Python module. IUCNN uses functions of the python modules tensorflow and npBNN which also need to be installed (via R). 

```r
reticulate::conda_install("r-reticulate","tensorflow=2.4")
reticulate::py_install("https://github.com/dsilvestro/npBNN/archive/v0.1.10.tar.gz", 
                       pip = TRUE)
```

# Prepare input data
IUCNN predicts the IUCN RL categories of Not Evaluated and Data Deficient species based on geographic occurrence records and a set of training species for which occurrence records and IUCN assessments that are available for a set of reference species (training data). The amount of training species necessary varies with the number of categories but in general "the more, the better". Ideally, the training dataset should comprise several hundred species or more, so a typical scenario will be to use all available plant species from a region, or all available species from a plant family. If the availability of training species is limited, a good option can be to reduce detail and predict Possibly threatened (IUCN categories "CR", "EN", and "VU") v. Not threatened species ("NT" and "LC").

Three types of input are necessary, which are easily available for many species: 

## 1. Geographic occurrence records of training species (training occurrences)
Occurrence records might be obtained from a variety of databases, For example, from field collections or public databases such BIEN (https://bien.nceas.ucsb.edu/bien/) or GBIF (www.gbif.org). GBIF data can be obtained from within R via the rgbif package, See [here](https://docs.ropensci.org/rgbif/articles/index.html) for a tutorial on how to do so. IUCNN needs a dataset with (at least) three columns, containing the species name, decimal longitude coordinates and decimal latitude coordinates. If you are interested in cleaning records from GBIF, you may want to have a look at this [blog post](https://data-blog.gbif.org/post/gbif-filtering-guide/) and check out the [CoordinateCleaner](https://github.com/ropensci/CoordinateCleaner) and [bRacatus](https://github.com/EduardoArle/bRacatus) packages. 

## 2. IUCN Global Red List assessment of the training species (training labels)
IUCN RL assessments for the training species can be obtained from IUCN, either via www.iucnredlist.org or via the rredlist package via R (preferred for many species). See [here](https://ropensci.org/tutorials/rredlist_tutorial/) for a tutorial on how to use rredlist. It is important, that all target label classes are well represented in the training data, which is rarely the case for IUCN data, since for instance "VU" and "NT" is rare. If the classes are to imbalanced, consider using possibly threatened (IUCN categories "CR", "EN", and "VU") v. not threatened species ("NT" and "LC"), or the supersampling option of the `iucnn_train_model` function.

## 3. Geographic occurrence records of the species for which the IUCN status should be predicted (predict occurrences)
Geographic occurrence for the target species, in the same format as for the training occurrences described above.

Example dataset are available with IUCNN: `data(training_occ)` (training occurrences), `data(training_labels)` (training labels) and `data(prediction_occ)`.

## Feature preparation
IUCNN uses per species traits to as features for the neural networks. The required input format is a data.frame with one row per species , one column containing the species name and any number of additional columns containing the numerical features for each species. In general, features might represent any trait, for instance from taxonomy (e.g., family), anatomy (e.g., body size), ecology (e.g., feeding guild) or conservation (e.g., population dynamics). However, since often only geographic occurrence data are available IUCNN contains functions to obtain default features from geo-referenced occurrence records alone, by combining them with publicly available environmental data. These default features informing on species range, climatic niche, human footprint and biomes. See Table 2 for a detailed list of all default features. Users may chose to use specific groups of features only via the `type` option of `iucnn_prepare_labels`. In this tutorial, we will use the example datasets from the Orchid family (Orchidaceae) provided with the IUCNN package, 

You can prepare the default features with a single call to `iucnn_prepare_features`

```r
data("training_occ") #geographic occurrences of species with IUCN assessment
data("prediction_occ")

features_train <- iucnn_prepare_features(training_occ) # Training features
features_predict <- iucnn_prepare_features(prediction_occ) # Prediction features
```

## Label preparation
IUCNN expects the labels for training as numerical categories. So, to use IUCN RL categories, those need to be converted to numeric in the right way. This can be done using the `iucnn_prepare_labels` function. The function converts the category labels as obtained from the IUCN RL into standardized numeric values, either on the detailed level of IUCN RL categories or the broader Possibly threatened/Not threatened level. See `?iucnn_prepare_labels` for more information. The labels will be converted into numeric categories following the `accepted_labels` argument, so for instance, in the default case: LC -> 0 and CR -> 4. If you change the accepted labels, the match will change accordingly.


```r
data("training_labels")

labels_train <- iucnn_prepare_labels(x = training_labels,
                                     y = features_train) # Training labels
```

# Running IUCNN
Running IUCNN consists of two steps: 1) training a neural network and 2) predicting the status of new species. IUCNN contains three different neural network approaches to predict the IUCN status of species, which can all be customized. We present the default approach here, see section "Customizing analyses" of this tutorial for details on how to train a Bayesian or regression type neural network. 

## Model training
Based on the training features and labels, IUCNN trains a neural network, via the `iucnn_train_model` function. There are multiple options to customize the design of the network, including among others the number of layers and the fraction of records used for testing and validation. The `iucnn_train_model` function will write a folder to the working directory containing the model and return summary statistics including cross-entropy loss and accuracy for the validation set, which can be used to compare the performance of different models.

The following code trains a neural network model with 3 hidden layers of 60, 60, and 20 nodes, with ReLU activation function. By specifying a seed (here, the default: 1234) we make sure the same subsets of data are designated as training, validation and test sets across different runs and model configurations (see below). The model with estimated weights will be saved in the current working directory. 


```r
res_1 <- iucnn_train_model(x = features_train,
                           lab = labels_train, 
                           path_to_output = "iucnn_model_1")
```

The `summary` and `plot` methods give an overview on the training process and model performance. 


```r
summary(res_1)
#> A model of type nn-class, trained on 702 species and 45 features.
#> 
#> Training accuracy: 0.6
#> Validation accuracy: 0.532
#> Accuracy on unseen data (test set): 0.549
#> 
#> Label detail: 5 Classes (detailed)
#> 
#> Label representation
#>   Label Input_count Estimated_count
#> 1     0          67              92
#> 2     1          12               0
#> 3     2          28               3
#> 4     3          51              80
#> 
#> Confusion matrix (rows test data and columns predicted):
#>    LC NT VU EN CR
#> LC 56  0  1 10  0
#> NT  9  0  0  3  0
#> VU 15  0  1 12  0
#> EN 11  0  1 39  0
#> CR  1  0  0 16  0
plot(res_1)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

## Predict IUCN Global Red List status
The trained model can then predict the conservation status of *Not Evaluated* and *Data Deficient* species with the `iucnn_predict_status` function. The output contains a data frame with species names and numeric labels (as generated with iucnn_prepare_labels).


```r
predictions <- iucnn_predict_status(x = features_predict, 
                                    model = res_1)

plot(predictions)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

It is important to remember the following when using IUCNN:

1. The resulting IUCNN categories are approximations only. While IUCNN has reached accuracies between 80 and 90% on the broad (threatened v non-threatened) level and up to 80% on the detailed level in some cases, the accuracy may be considerably lower in other cases, which means that some species will be mis-classified.

2. IUCNN is indifferent to the provided features. On the one hand this means that any species traits for which data is available can be used, but on the other hand this means that thought is needed in the choice of the features. The default features of IUCNN are usually a safe choice. The number of features is not limited, but currently IUCNN does not support missing values in the feature table and removes species with missing values. 

3. IUCNN is indifferent to the relation between training and test data. So it is possible to use training data from Palearctic birds to predict the conservation status of South American Nematodes. This is not recommended. Instead, a better approach will be to predict the conservation status of species, from training data of the same genus, order, or family. Alternatively, training data could be chosen on geographic region or functional aspects (e.g., feeding guild or body size). However some inclusion of taxonomy/evolutionary history for the choice of training data is recommended.

4. The amount of training data is important. The more the better. Minimum several hundred training species with a more or less equal distribution on the label classes should be included. If training data is limited, the broader Threatened/Not threatened level is recommended. 

5. If the proportion of the IUCN RL categories is imbalanced in the training data, the neural networks may be biased towards reproducing these frequencies in the prediction, especially if the imbalance of categories or the difference in category frequencies among training and prediction set are large. To avoid this category frequencies should be balanced in the training data if possible. Otherwise the use of the `supersampling` or option a `nn-reg` type model, or a limitation to the broader Possibly threatened/Not threatened level of detail may remedy the issue. 

6. IUCNN predictions are not equivalent to full IUCN Red List assessments. We see the main purpose of IUCNN in 1) identifying species that will likely need conservation action to trigger a full IUCN assessment, and 2) provide large-scale overviews on the extinction risk in a given taxonomic group, for instance in a macro-ecological and macro-evolutionary context.

# Customizing IUCNN analyses
IUCNN contains multiple options to customize the steps of the analyses to adapt the fully connected neural networks to the peculiarities of IUCN RL and species distribution data. Below we describe the most important options to customize 1) feature and label preparation, 2) model training and testing, and 3) status prediction. The most important steps and options to customize an IUCNN analysis are summarized in Table 1.

Table 1. Critical steps to customize an IUCNN analysis and relevant considerations at each point.

| Step | Function(s) | Argument | Description | User consideration |
|---|---| --- |----|---|
| Feature design | iucnn\_prepare\_features, iucnn\_feature\_importance | \- | Defines the features to be extracted from the provided occurrence records. Available defaults are: biome presence, bioclim variables, human footprint and geographic features. Averaged per species and rescaled. | What determines extinction risk for target group? Which data are available? What do the results of feature importance suggest? |
| Label detail | prep\_labels | level, threatened | Into how many different categories should the species be classified. Can be any number, defaults support full IUCN categories (LC, NT, VU, EN, CR) or binary (Possibly threatened v. Not threatened) | Which detail is needed? Is the accuracy of the detailed level sufficient for the target application? |
| Model type | iucnn\_train\_model | mode | Which model framework should be applied: a categorical classification, a classification taking the ordinal number of categories into account, or a classification based in a Bayesian framework | How many target categories are there? How important is a high accuracy for intermediate classes (e.g. VU, NT, EN)? How important is the uncertainty estimation for each species? |
| Model structure | iucnn\_train\_model | validation fraction | The fraction of the input training data used for validation (v. training). | Which fraction of the data should be used for validation (and hence not training)? How large is the training data? |
| Model structure | iucnn\_train\_model | cv\_fold | The number of folds used for cross-validation. For instance if = 5, the data is divided into 5 folds with 20% of the data used for validation in each run. | How large is the training data?| At a given number of folds, will the subsets still include all label classes? |
| Model structure | iucnn\_train\_model | n\_layers | The number of hidden layers and nodes in the neural network. | How complex should the model be? |
| Model structure | iucnn\_train\_model | balance\_classes | Should the frequency of the class labels in the training data be balanced using supersampling? | How imbalanced are the class labels? Will the frequency of class labels differ between training data and prediction data set? |
| Model structure | iucnn\_train\_model | act\_f\_out/act\_f | The activation function of the neural network | Which relationship between features and labels is expected? |
| Model structure | iucnn\_train\_model | label\_stretch\_factor | The factor to stretch input class labels to | Am I using a nn-reg type model?  Does model testing suggest an effect on model accuracy? |
| Model structure | iucnn\_train\_model | drop\_out rate | The number of nodes to be removed in individual epochs. Necessary if a target threshold is to be used with a nn-class model (not bnn-class though) | Is a target accuracy to be sued for prediction? How many nodes are there in the model? |
| Prediction | predict\_iucnn | target\_acc | Defines an overall target accuracy for the model Species which cannot be classified with enough certainty to reach this threshold are labels as data deficient. | Which error rate is acceptable? Which proportion of species needs to be included? |


## 1) Features and Labels
### Add and remove feature blocks
The default features are selected based on empirical tests on relevance for different taxa and regions. However, for some analyses only part of the features may be relevant. You can exclude feature blocks using the `type` argument of the `iucnn_prepare_features` function. For instance, to exclude the biome features:


```r
features_train2 <- iucnn_prepare_features(training_occ, 
                                          type = c("geographic", 
                                                   "climate", 
                                                   "humanfootprint"))
```

### Prepare features individually
If more control over feature preparation is necessary, each feature block can be obtained by an individual function.

Table 2. Functions to obtain default features and options to customize the features.

|Feature block|Function name|Options to customize|
|---|---|---|
|Geographic|`iucnn_geographic_features`|-|
|Biomes|`iucnn_biome_features`|change the reference dataset of biomes (biome_input, biome.id), remove biomes without any species occurrence (remove_zeros)|
|Climate|`iucnn_climate_features`|the amount of bioclim variables from the default source to be included (type), the resolution of the default input data (res)|
|Human footprint|`iucnn_footprint_features`|chose the time points from the default source (year), the break points for the different footprint categories (breaks, by default approximately quantiles on the global footprint dataset) or a default source for human footprint (footp_input)|
|Geographic bias|`iucnn_bias_features`|The resolution of the bias raster (res) or providing a template raster (ras)|

For instance:


```r
clim_features <- iucnn_climate_features(x = training_occ, 
                                        type = "selected")

clim_features2 <- iucnn_climate_features(x = training_occ, 
                                         type = "all")
```

### Use custom features
It is also possible to provide features unrelated to the default features. They may contain any continuous or categorical features, but some processing will be needed. The format needs to be a data.frame with a compulsory column containing the species name. Continuous variables should be rescaled to cover a similar range, whereas categorical features should be coded binary (present/absent, as the results of `iucnn_biome_features`).

For instance:


```r
feat <- data.frame(species = c("Adansonia digitata", "Ceiba pentandra"),
 max_plant_size_m = c(25, 50),
 africa = c(1,1),
 south_america = c(0,1),
 fraction_of_records_in_protected_area = c(25, 75))
```

Table 3. Description of the default features included in `iucnn_prepare_features`. All continuous variables are rescaled to a similar range.

| Feature | Block | Name | Description |
|---|---|---|---|
|tot_occ|Geographic|Number of occurrences|The total number of occurrences available for this species|
|uni_occ|Geographic|Number of geographically unique occurrences|The number of geographically unique records available for this species|
|mean_lat|Geographic|Mean latitude|The mean latitude of all records of this species|
|mean_lon|Geographic|Mean longitude|The mean longitude of all records of this species|
|lat_range|Geographic|Latitudinal range|The latitudinal range (.95 quantile - .05 quantile).|
|lon_range|Geographic|Longitudinal range|The longitudinal range (.95 quantile - .05 quantile).|
|alt_hemisphere|Geographic|The hemisphere|0 = Southern hemisphere, 1 = Northern hemisphere|
|eoo|Geographic|Extend of Occurrence|The extend of occurrence. Calculated by rCAT. For species with less than 3 records set to AOO|
|aoo|Geographic|Area of Occupancy|The area of occupancy, as the sum of area of 4sqkm grid cells, where the species occurs|
|1|Biome|Tropical & Subtropical Moist Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|2|Biome|Tropical & Subtropical Dry Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|3|Biome|Tropical & Subtropical Coniferous Forests|Are at least 5% of the species records present in this biome?|
|4|Biome|Temperate Broadleaf & Mixed Forests|Are at least 5% of the species records present in this biome?|
|5|Biome|Temperate Conifer Forests|Are at least 5% of the species records present in this biome?|
|6|Biome|Boreal Forests/Taiga|Are at least 5% of the species records present in this biome?|
|7|Biome|Tropical & Subtropical Grasslands, Savannas & Shrublands|Are at least 5% of the species records present in this biome?|
|8|Biome|Temperate Grasslands, Savannas & Shrublands|Are at least 5% of the species records present in this biome?|
|9|Biome|Flooded Grasslands & Savannas|Are at least 5% of the species records present in this biome?|
|10|Biome|Montane Grasslands & Shrublands|Are at least 5% of the species records present in this biome?|
|11|Biome|Tundra|Are at least 5% of the species records present in this biome?|
|12|Biome|Mediterranean Forests, Woodlands & Scrub|Are at least 5% of the species records present in this biome?|
|13|Biome|Deserts & Xeric Shrublands|Are at least 5% of the species records present in this biome?|
|14|Biome|Mangroves|Tropical & Subtropical Moist Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|98|Biome|Lake|Are at least 5% of the species records present in this biome?|
|99|Biome|Rock and ice|Are at least 5% of the species records present in this biome?|
|bio1|Climate|Annual Mean Temperature|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio4|Climate|Temperature Seasonality|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio11|Climate|Mean Temperature of Coldest Quarter|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio12|Climate|Annual Precipitation|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio15|Climate|Precipitation Seasonality|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio17|Climate|Precipitation of Driest Quarter|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|range_bio1|Climate|Range of annual Mean Temperature|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio4|Climate|Range of temperature Seasonality|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio11|Climate|Range of mean Temperature of Coldest Quarter|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio12|Climate|Range of annual Precipitation|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio15|Climate|Range of precipitation Seasonality|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio17|Climate|Range of precipitation of Driest Quarter|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|humanfootprint_1993_1|Human footprint | Human footprint year 1993 lowest impact|The fraction of records in areas of the lowest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_2|Human footprint|Human footprint year 1993 intermediate impact 1|The fraction of records in areas of the second lowest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_3|Human footprint|Human footprint year 1993 intermediate impact 2|The fraction of records in areas of the second highest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_4|Human footprint|Human footprint year 1993 highest impact|The fraction of records in areas of the highest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_1|Human footprint|Human footprint year 2009 lowest impact|The fraction of records in areas of the lowest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_2|Human footprint|Human footprint year 2009 intermediate impact 1|The fraction of records in areas of the second lowest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_3|Human footprint|Human footprint year 2009 intermediate impact 2|The fraction of records in areas of the second highest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_4|Human footprint|Human footprint year 2009 highest impact|The fraction of records in areas of the highest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|

### Labels: Full categories vs Threatened/Not threatened
The `iucnn_prepare_labels` function may accepted any custom labels as long as they are included in the `accepted_labels` option. It also can provide a classification into threatened/non-threatened, via the `level` and `threatened` options. On the broader level the model accuracy is usually significantly higher.

For instance:


```r
labels_train <- iucnn_prepare_labels(training_labels,
                                     y = features, 
                                     level = "broad")
```

## 2) Model training - NN regression model
### Customizing model parameters
The `iucnn_train_model` function contains various options to customize the neural network, including among other the fraction of validation and test data, the maximum number of epochs, the number of layers and nodes, the activation function , dropout and randomization of the input data. See `?iucnn_train_model` for a comprehensive list of options and their description. By default, `iucnn_train_model` trains a neural network with three hidden layers with 50, 30 and 10 nodes and a sigmoid as activation function. Depending on your dataset different networks may improve performance. For instance, you can set up a different model with 1 hidden layer of 60 nodes, a sigmoid activation function and without using a bias node in the first hidden layer.


```r
res_2 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 dropout_rate = 0.3,
 path_to_output= "iucnn_model_2",
 n_layers = "60",
 use_bias = FALSE,
 act_f = "sigmoid")
```

You can compare the validation loss of the models using `res_1$validation_loss` and `res_2$validation_loss`. Model 2 in this case yields a lower validation loss and is therefore preferred. Once you chose the preferred model configuration based on validation loss, we can check test accuracy of best model: `res_2$test_accuracy`. The `iucnn_train_model` function contains various options to adapt the model. See `?iucnn_train_model` for more detail. 

### Changing the modeling algorithm
There are three neural network algorithms implemented in IUCNN. Besides the default classifier approach based on a tensorflow implementation, these are a Bayesian neural network classifier and a regression type neural network.

The Bayesian approach has the advantage that it returns true probabilities for the classification of species into the relative output classes (e.g. 80% probability of a species to be LC). We consider this approach more suitable for classification of species into IUCN categories, than the default option. It will need more time for model training and should best be applied once you have identified the best model parameters using the default approach. You can run a BNN setting the `mode` option of `iucnn_train_model` to `"bnn-class"`.


```r
res_3 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 path_to_output = "iucnn_model_3",
 mode = 'bnn-class')
```

IUCNN also offers the option to train a NN regression model instead of a classifier. Since the IUCN threat statuses constitute a list of ordinal categories sorted by increasing threat level, we can model the task of estimating these categories as a regression problem. Such a model can be trained with the `iucnn_train_model()` function, specifying to train a regression model by setting `mode = 'nn-reg'`.


```r
res_4 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 path_to_output = "iucnn_model_4",
 mode = 'nn-reg',
 rescale_features = TRUE)
```

### Feature importance
The `iucnn_feature_importance` function can be used to gauge the importance of different feature blocks or individual features for model performance. The function implements the permutation feature importance technique, which randomly reshuffles the values within individual features or blocks of features and evaluate how this randomization affects the models prediction accuracy. If a given feature (or block of features) is important for the models ability to predict, randomizing this feature will lead to a large drop in prediction accuracy. When using `iucnn_feature_importance` with features other than the default, feature blocks can be defined using the `feature_blocks` option.

```r
fi <- iucnn_feature_importance(x = res_1)
plot(fi)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

### Model testing
Before training the final model used for predicting the conservation status of not evaluated species, it is recommended to use the `iucnn_modeltest` function for finding the best settings for your model and dataset. This process, often referred to as hyperparameter tuning, is an essential step for building the most suitable model for the prediction task. The `iucnn_modeltest` function allows you to provide any settings for `iucnn_train_model` as vectors, which will lead the function to train a separate model for each provided setting. The function will explore all possible permutations of the provided settings, so that the following command results in 9 different models being trained:


```r
modeltest_results <- iucnn_modeltest(features,
 labels,
 dropout_rate = c(0.0,0.1,0.3),
 n_layers = c('30','40_20','50_30_10'))
```

The model specifications and settings of each tested model are written to a log-file and can be inspected with the `iucnn_best_model` function, to decide which model settings to pick as best model. Different criteria for picking the best model can be selected, such as best prediction accuracy, best predicted over-all status distribution, lowest weighted mis-classification error, etc.


```r
best_m <- iucnn_best_model(modeltest_results, criterion='val_acc')
```

After model testing, it is necessary to retrain using the model specifications identified by `iucnn_best_model`. It is possible to take the respective settings directly from the output of the `iucnn_best_model` function, via the `production-model` argument of `iucnn_train_model`.


```r
# Train the best model on all training data for prediction
m_prod <- iucnn_train_model(train_feat,
                      train_lab,
                      production_model = m_best,
                      overwrite = TRUE)
```

The production model can then be used for the final predictions.


```r
# Predict RL categories for target species
pred <- iucnn_predict_status(pred_feat,
                      m_prod)
plot(pred)
```


## 3) Status prediction
The `iucnn_predict_status` function offers options to customize the predictions. The most important option in many cases is `target_acc`, which allows to set an overall target-accuracy threshold that the model needs to achieve. This option is only available for nn-class and nn-reg models that were trained using dropout (see help function of `iucnn_train_model` for more explanation), as well as for all bnn-class models. The `target_acc` will be achieved by the model being more selective with making a category call for a given instance. All species that cannot be classified with enough certainty to reach this target accuracy will be classified as NA (Not Assessed).

```r
pred_2 <- iucnn_predict_status(x = features_predict, 
 target_acc = 0.7,
 model = res_2)
plot(pred_2)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

Furthermore, you can turn off the `return_IUCN` option if to return the numerical labels instead of the IUCNN RL category labels.

```r
pred_3 <- iucnn_predict_status(x = features_predict, 
 model = res_2,
 return_IUCN = FALSE)
```

The output of the `iucnn_predict_status` function is an "iucnn_predictions" object, that contains several output objects. The predicted labels of the individual instances are accessible with `pred_2$class_predictions` and label probabilities estimated by the neural network via `pred_2$mc_dropout_probs` for "nn-class" and "nn-reg" with dropout, or `pred_2$posterior_probs` for "bnn-class". For more detail, the `pred_2$raw_predictions` object contains the individual label probabilities resulting from the softmax output layer in case of "nn-class", or the regressed labels in case of "nn-reg".

## 4) The number of species per category
Another statistic that can be extracted from the "iucnn_predictions" object is the overall category distribution predicted for the given prediction instances. This can be accessed with `pred_2$pred_cat_count`, which shows the distribution of the label predictions, including the count of species that could not be predicted given the chosen `target_acc`.

Another statistic are the `pred2$sampled_cat_freqs` (only available for dropout models and all "bnn-class" models, see above), which show the class distribution as sampled from the `pred_2$mc_dropout_probs` or the `pred_2$posterior_probs` (for "bnn-class" models). The difference between `pred2$sampled_cat_freqs` and `pred_2$mc_dropout_probs`/`pred_2$posterior_probs` is that the former represents the counts of the best labels determined for each instance, whereas the latter represents labels sampled from the predicted label probabilities, which also proportionally samples the labels for a given instance that do receive the maximum label probability. The latter is done repeatedly to include the stochasticity of the random sampling of classes from the given probability vectors. The `pred_2$mc_dropout_probs`/`pred_2$posterior_probs` can be used to plot histograms of the estimates for each class, and can be reported as uncertainty intervals around the number of species in each class for the set of species that were predicted.


# Training and prediction using a convolutional neural network
Instead of the fully connected neural networks presented above, IUCNN also implements a prediction algorithm using Convolutional Neural Networks (CNNs). The input data for prediction are then rasterized per-species grids of occurrence numbers. Since the input features and network structure is different CNNs are implemented in IUCNN with a separate workflow.

## 1. Feature preparation
Features are prepared using the `iucnn_cnn_features` function which will count the number of occurrence records in a custom raster with the same extent as the species occurrences. The function can be simply run on a data frame of training and prediction occurrences separately. Yet, since the extent of these two datasets is likely different, it is recommendable to create a custom raster beforehand. Here, we will use a simple lat/lon raster, but users may provide coordinates and a raster in any suitable coordinate reference system, as long as they agree between occurrence coordinates and raster.


```r
# preapre custom raster, you can split this step if training and test occurrences have the same extent
library(terra)
data("training_occ")
data("prediction_occ")

# find the minimum latitude and longitude values for the extent of the raster
min_lon <- min(c(min(training_occ$decimallongitude), 
                 min(prediction_occ$decimallongitude)))
max_lon <- max(c(max(training_occ$decimallongitude), 
                 max(prediction_occ$decimallongitude)))
min_lat <- min(c(min(training_occ$decimallatitude), 
                 min(prediction_occ$decimallatitude)))
max_lat <- max(c(max(training_occ$decimallatitude), 
                 max(prediction_occ$decimallatitude)))
## set the coordinate reference system
ras <- rast(crs = "+proj=longlat +datum=WGS84")
## set raster extent
ext(ras) <- c(min_lon,
              max_lon, 
              min_lat, 
              max_lat)
## set raster resolution
res(ras) <- 1 # the resolution in CRS units, in this case degrees lat/lon


 # Training features
cnn_features <- iucnn_cnn_features(x = training_occ, 
                                   y = ras)

# Prediction features
cnn_features_predict <- iucnn_cnn_features(x = prediction_occ,
                                           y = ras)
 # Training labels
cnn_labels <- iucnn_prepare_labels(x = training_labels,
                                   y = cnn_features)
```


## 2. Model training

```r
trained_model <- iucnn_cnn_train(cnn_features,
                                cnn_labels,
                                overwrite = TRUE,
                                dropout_rate = 0.1,
                                optimize_for = 'accuracy')

plot(trained_model)
summary(trained_model)
```

## 3. Status prediction

```r
pred <- iucnn_predict_status(cnn_features_predict,
                             trained_model,
                             target_acc = 0.0
                             )
```

---
title: Approximate Red List assessments with IUCNN
#output: rmarkdown::html_vignette
output: pdf_document
vignette: >
 %\VignetteIndexEntry{Approximate Red List assessments with IUCNN}
 %\VignetteEngine{knitr::knitr}
 %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
 collapse = TRUE,
 comment = "#>",
 warning = FALSE,
 message = FALSE,
 fig.width = 8
)
```

# Background
The Red List of the International Union for the Conservation of nature (www.iucn.org, IUCN RL), is arguably one of the most thorough and widely used tools to assess the global extinction risk of species. However, the IUCN RL assessment process---usually performed by a group of specialists for each taxonomic group, or professional assessors---are time consuming, and therefore only a small fraction of global biodiversity has been evaluated for the IUCN RL, with a strong bias towards vertebrates and certain regions. These biases and the low fraction of species evaluated, bias conservation towards evaluated groups and prevent synthetic, large-scale ecological and biogeographic analyses of extinction risk. 

IUCNN implements neural networks to predict the IUCN status of so far not evaluated or data deficient species based on species traits. IUCNN models are trained on the existing IUCN RL assessments and any traits may be used for prediction, although IUCNN implements a workflow based solely on publicly available geo-referenced species occurrence records and environmental data. Typical examples for the application of IUCNN are to predict the conservation status of a large number of species, to approximate extinction risk or number of threatened species in a region or specific taxonomic group for synthetic analyses or to predict the IUCN category of individual species of interest for systematic or ecological case studies. 

```{r setup}
library(IUCNN)
library(magrittr)
library(dplyr)
```

# Installation
IUCNN uses R and python. All software needed can be installed via R.

1. install IUCNN directly from Github using devtools. 
```{r, eval = FALSE}
install.packages("devtools")
library(devtools)
library(IUCNN)
```

2. Python needs to be installed, for instance using miniconda and reticulated from within R (this will need c. 3 GB disk space).
If problems occur at this step, check the excellent [documentation of reticulate](https://rstudio.github.io/reticulate/index.html).
```{r, eval = FALSE}
install.packages(reticulate)
library("reticulate")
install_miniconda()
```

If python has been installed before, you can specify the python version to sue with `reticulate::use_python()`

3. Install the tensorflow Python module. IUCNN uses functions of the python modules tensorflow and npBNN which also need to be installed (via R). 
```{r, eval = FALSE}
reticulate::conda_install("r-reticulate","tensorflow=2.4")
reticulate::py_install("https://github.com/dsilvestro/npBNN/archive/v0.1.10.tar.gz", 
                       pip = TRUE)
```

# Prepare input data
IUCNN predicts the IUCN RL categories of Not Evaluated and Data Deficient species based on geographic occurrence records and a set of training species for which occurrence records and IUCN assessments that are available for a set of reference species (training data). The amount of training species necessary varies with the number of categories but in general "the more, the better". Ideally, the training dataset should comprise several hundred species or more, so a typical scenario will be to use all available plant species from a region, or all available species from a plant family. If the availability of training species is limited, a good option can be to reduce detail and predict Possibly threatened (IUCN categories "CR", "EN", and "VU") v. Not threatened species ("NT" and "LC").

Three types of input are necessary, which are easily available for many species: 

## 1. Geographic occurrence records of training species (training occurrences)
Occurrence records might be obtained from a variety of databases, For example, from field collections or public databases such BIEN (https://bien.nceas.ucsb.edu/bien/) or GBIF (www.gbif.org). GBIF data can be obtained from within R via the rgbif package, See [here](https://docs.ropensci.org/rgbif/articles/index.html) for a tutorial on how to do so. IUCNN needs a dataset with (at least) three columns, containing the species name, decimal longitude coordinates and decimal latitude coordinates. If you are interested in cleaning records from GBIF, you may want to have a look at this [blog post](https://data-blog.gbif.org/post/gbif-filtering-guide/) and check out the [CoordinateCleaner](https://github.com/ropensci/CoordinateCleaner) and [bRacatus](https://github.com/EduardoArle/bRacatus) packages. 

## 2. IUCN Global Red List assessment of the training species (training labels)
IUCN RL assessments for the training species can be obtained from IUCN, either via www.iucnredlist.org or via the rredlist package via R (preferred for many species). See [here](https://ropensci.org/tutorials/rredlist_tutorial/) for a tutorial on how to use rredlist. It is important, that all target label classes are well represented in the training data, which is rarely the case for IUCN data, since for instance "VU" and "NT" is rare. If the classes are to imbalanced, consider using possibly threatened (IUCN categories "CR", "EN", and "VU") v. not threatened species ("NT" and "LC"), or the supersampling option of the `iucnn_train_model` function.

## 3. Geographic occurrence records of the species for which the IUCN status should be predicted (predict occurrences)
Geographic occurrence for the target species, in the same format as for the training occurrences described above.

Example dataset are available with IUCNN: `data(training_occ)` (training occurrences), `data(training_labels)` (training labels) and `data(prediction_occ)`.

## Feature preparation
IUCNN uses per species traits to as features for the neural networks. The required input format is a data.frame with one row per species , one column containing the species name and any number of additional columns containing the numerical features for each species. In general, features might represent any trait, for instance from taxonomy (e.g., family), anatomy (e.g., body size), ecology (e.g., feeding guild) or conservation (e.g., population dynamics). However, since often only geographic occurrence data are available IUCNN contains functions to obtain default features from geo-referenced occurrence records alone, by combining them with publicly available environmental data. These default features informing on species range, climatic niche, human footprint and biomes. See Table 2 for a detailed list of all default features. Users may chose to use specific groups of features only via the `type` option of `iucnn_prepare_labels`. In this tutorial, we will use the example datasets from the Orchid family (Orchidaceae) provided with the IUCNN package, 

You can prepare the default features with a single call to `iucnn_prepare_features`
```{r, results='hide'}
data("training_occ") #geographic occurrences of species with IUCN assessment
data("prediction_occ")

features_train <- iucnn_prepare_features(training_occ) # Training features
features_predict <- iucnn_prepare_features(prediction_occ) # Prediction features

```

## Label preparation
IUCNN expects the labels for training as numerical categories. So, to use IUCN RL categories, those need to be converted to numeric in the right way. This can be done using the `iucnn_prepare_labels` function. The function converts the category labels as obtained from the IUCN RL into standardized numeric values, either on the detailed level of IUCN RL categories or the broader Possibly threatened/Not threatened level. See `?iucnn_prepare_labels` for more information. The labels will be converted into numeric categories following the `accepted_labels` argument, so for instance, in the default case: LC -> 0 and CR -> 4. If you change the accepted labels, the match will change accordingly.

```{r}
data("training_labels")

labels_train <- iucnn_prepare_labels(x = training_labels,
                                     y = features_train) # Training labels
```

# Running IUCNN
Running IUCNN consists of two steps: 1) training a neural network and 2) predicting the status of new species. IUCNN contains three different neural network approaches to predict the IUCN status of species, which can all be customized. We present the default approach here, see section "Customizing analyses" of this tutorial for details on how to train a Bayesian or regression type neural network. 

## Model training
Based on the training features and labels, IUCNN trains a neural network, via the `iucnn_train_model` function. There are multiple options to customize the design of the network, including among others the number of layers and the fraction of records used for testing and validation. The `iucnn_train_model` function will write a folder to the working directory containing the model and return summary statistics including cross-entropy loss and accuracy for the validation set, which can be used to compare the performance of different models.

The following code trains a neural network model with 3 hidden layers of 60, 60, and 20 nodes, with ReLU activation function. By specifying a seed (here, the default: 1234) we make sure the same subsets of data are designated as training, validation and test sets across different runs and model configurations (see below). The model with estimated weights will be saved in the current working directory. 

```{r}
res_1 <- iucnn_train_model(x = features_train,
                           lab = labels_train, 
                           path_to_output = "iucnn_model_1")
```

The `summary` and `plot` methods give an overview on the training process and model performance. 

```{r}
summary(res_1)
plot(res_1)
```

## Predict IUCN Global Red List status
The trained model can then predict the conservation status of *Not Evaluated* and *Data Deficient* species with the `iucnn_predict_status` function. The output contains a data frame with species names and numeric labels (as generated with iucnn_prepare_labels).

```{r}
predictions <- iucnn_predict_status(x = features_predict, 
                                    model = res_1)

plot(predictions)
```

It is important to remember the following when using IUCNN:

1. The resulting IUCNN categories are approximations only. While IUCNN has reached accuracies between 80 and 90% on the broad (threatened v non-threatened) level and up to 80% on the detailed level in some cases, the accuracy may be considerably lower in other cases, which means that some species will be mis-classified.

2. IUCNN is indifferent to the provided features. On the one hand this means that any species traits for which data is available can be used, but on the other hand this means that thought is needed in the choice of the features. The default features of IUCNN are usually a safe choice. The number of features is not limited, but currently IUCNN does not support missing values in the feature table and removes species with missing values. 

3. IUCNN is indifferent to the relation between training and test data. So it is possible to use training data from Palearctic birds to predict the conservation status of South American Nematodes. This is not recommended. Instead, a better approach will be to predict the conservation status of species, from training data of the same genus, order, or family. Alternatively, training data could be chosen on geographic region or functional aspects (e.g., feeding guild or body size). However some inclusion of taxonomy/evolutionary history for the choice of training data is recommended.

4. The amount of training data is important. The more the better. Minimum several hundred training species with a more or less equal distribution on the label classes should be included. If training data is limited, the broader Threatened/Not threatened level is recommended. 

5. If the proportion of the IUCN RL categories is imbalanced in the training data, the neural networks may be biased towards reproducing these frequencies in the prediction, especially if the imbalance of categories or the difference in category frequencies among training and prediction set are large. To avoid this category frequencies should be balanced in the training data if possible. Otherwise the use of the `supersampling` or option a `nn-reg` type model, or a limitation to the broader Possibly threatened/Not threatened level of detail may remedy the issue. 

6. IUCNN predictions are not equivalent to full IUCN Red List assessments. We see the main purpose of IUCNN in 1) identifying species that will likely need conservation action to trigger a full IUCN assessment, and 2) provide large-scale overviews on the extinction risk in a given taxonomic group, for instance in a macro-ecological and macro-evolutionary context.

# Customizing IUCNN analyses
IUCNN contains multiple options to customize the steps of the analyses to adapt the fully connected neural networks to the peculiarities of IUCN RL and species distribution data. Below we describe the most important options to customize 1) feature and label preparation, 2) model training and testing, and 3) status prediction. The most important steps and options to customize an IUCNN analysis are summarized in Table 1.

Table 1. Critical steps to customize an IUCNN analysis and relevant considerations at each point.

| Step | Function(s) | Argument | Description | User consideration |
|---|---| --- |----|---|
| Feature design | iucnn\_prepare\_features, iucnn\_feature\_importance | \- | Defines the features to be extracted from the provided occurrence records. Available defaults are: biome presence, bioclim variables, human footprint and geographic features. Averaged per species and rescaled. | What determines extinction risk for target group? Which data are available? What do the results of feature importance suggest? |
| Label detail | prep\_labels | level, threatened | Into how many different categories should the species be classified. Can be any number, defaults support full IUCN categories (LC, NT, VU, EN, CR) or binary (Possibly threatened v. Not threatened) | Which detail is needed? Is the accuracy of the detailed level sufficient for the target application? |
| Model type | iucnn\_train\_model | mode | Which model framework should be applied: a categorical classification, a classification taking the ordinal number of categories into account, or a classification based in a Bayesian framework | How many target categories are there? How important is a high accuracy for intermediate classes (e.g. VU, NT, EN)? How important is the uncertainty estimation for each species? |
| Model structure | iucnn\_train\_model | validation fraction | The fraction of the input training data used for validation (v. training). | Which fraction of the data should be used for validation (and hence not training)? How large is the training data? |
| Model structure | iucnn\_train\_model | cv\_fold | The number of folds used for cross-validation. For instance if = 5, the data is divided into 5 folds with 20% of the data used for validation in each run. | How large is the training data?| At a given number of folds, will the subsets still include all label classes? |
| Model structure | iucnn\_train\_model | n\_layers | The number of hidden layers and nodes in the neural network. | How complex should the model be? |
| Model structure | iucnn\_train\_model | balance\_classes | Should the frequency of the class labels in the training data be balanced using supersampling? | How imbalanced are the class labels? Will the frequency of class labels differ between training data and prediction data set? |
| Model structure | iucnn\_train\_model | act\_f\_out/act\_f | The activation function of the neural network | Which relationship between features and labels is expected? |
| Model structure | iucnn\_train\_model | label\_stretch\_factor | The factor to stretch input class labels to | Am I using a nn-reg type model?  Does model testing suggest an effect on model accuracy? |
| Model structure | iucnn\_train\_model | drop\_out rate | The number of nodes to be removed in individual epochs. Necessary if a target threshold is to be used with a nn-class model (not bnn-class though) | Is a target accuracy to be sued for prediction? How many nodes are there in the model? |
| Prediction | predict\_iucnn | target\_acc | Defines an overall target accuracy for the model Species which cannot be classified with enough certainty to reach this threshold are labels as data deficient. | Which error rate is acceptable? Which proportion of species needs to be included? |


## 1) Features and Labels
### Add and remove feature blocks
The default features are selected based on empirical tests on relevance for different taxa and regions. However, for some analyses only part of the features may be relevant. You can exclude feature blocks using the `type` argument of the `iucnn_prepare_features` function. For instance, to exclude the biome features:

```{r, eval = FALSE}
features_train2 <- iucnn_prepare_features(training_occ, 
                                          type = c("geographic", 
                                                   "climate", 
                                                   "humanfootprint"))
```

### Prepare features individually
If more control over feature preparation is necessary, each feature block can be obtained by an individual function.

Table 2. Functions to obtain default features and options to customize the features.

|Feature block|Function name|Options to customize|
|---|---|---|
|Geographic|`iucnn_geographic_features`|-|
|Biomes|`iucnn_biome_features`|change the reference dataset of biomes (biome_input, biome.id), remove biomes without any species occurrence (remove_zeros)|
|Climate|`iucnn_climate_features`|the amount of bioclim variables from the default source to be included (type), the resolution of the default input data (res)|
|Human footprint|`iucnn_footprint_features`|chose the time points from the default source (year), the break points for the different footprint categories (breaks, by default approximately quantiles on the global footprint dataset) or a default source for human footprint (footp_input)|
|Geographic bias|`iucnn_bias_features`|The resolution of the bias raster (res) or providing a template raster (ras)|

For instance:

```{r, eval = FALSE}
clim_features <- iucnn_climate_features(x = training_occ, 
                                        type = "selected")

clim_features2 <- iucnn_climate_features(x = training_occ, 
                                         type = "all")
```

### Use custom features
It is also possible to provide features unrelated to the default features. They may contain any continuous or categorical features, but some processing will be needed. The format needs to be a data.frame with a compulsory column containing the species name. Continuous variables should be rescaled to cover a similar range, whereas categorical features should be coded binary (present/absent, as the results of `iucnn_biome_features`).

For instance:

```{r, eval = FALSE}
feat <- data.frame(species = c("Adansonia digitata", "Ceiba pentandra"),
 max_plant_size_m = c(25, 50),
 africa = c(1,1),
 south_america = c(0,1),
 fraction_of_records_in_protected_area = c(25, 75))
```

Table 3. Description of the default features included in `iucnn_prepare_features`. All continuous variables are rescaled to a similar range.

| Feature | Block | Name | Description |
|---|---|---|---|
|tot_occ|Geographic|Number of occurrences|The total number of occurrences available for this species|
|uni_occ|Geographic|Number of geographically unique occurrences|The number of geographically unique records available for this species|
|mean_lat|Geographic|Mean latitude|The mean latitude of all records of this species|
|mean_lon|Geographic|Mean longitude|The mean longitude of all records of this species|
|lat_range|Geographic|Latitudinal range|The latitudinal range (.95 quantile - .05 quantile).|
|lon_range|Geographic|Longitudinal range|The longitudinal range (.95 quantile - .05 quantile).|
|alt_hemisphere|Geographic|The hemisphere|0 = Southern hemisphere, 1 = Northern hemisphere|
|eoo|Geographic|Extend of Occurrence|The extend of occurrence. Calculated by rCAT. For species with less than 3 records set to AOO|
|aoo|Geographic|Area of Occupancy|The area of occupancy, as the sum of area of 4sqkm grid cells, where the species occurs|
|1|Biome|Tropical & Subtropical Moist Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|2|Biome|Tropical & Subtropical Dry Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|3|Biome|Tropical & Subtropical Coniferous Forests|Are at least 5% of the species records present in this biome?|
|4|Biome|Temperate Broadleaf & Mixed Forests|Are at least 5% of the species records present in this biome?|
|5|Biome|Temperate Conifer Forests|Are at least 5% of the species records present in this biome?|
|6|Biome|Boreal Forests/Taiga|Are at least 5% of the species records present in this biome?|
|7|Biome|Tropical & Subtropical Grasslands, Savannas & Shrublands|Are at least 5% of the species records present in this biome?|
|8|Biome|Temperate Grasslands, Savannas & Shrublands|Are at least 5% of the species records present in this biome?|
|9|Biome|Flooded Grasslands & Savannas|Are at least 5% of the species records present in this biome?|
|10|Biome|Montane Grasslands & Shrublands|Are at least 5% of the species records present in this biome?|
|11|Biome|Tundra|Are at least 5% of the species records present in this biome?|
|12|Biome|Mediterranean Forests, Woodlands & Scrub|Are at least 5% of the species records present in this biome?|
|13|Biome|Deserts & Xeric Shrublands|Are at least 5% of the species records present in this biome?|
|14|Biome|Mangroves|Tropical & Subtropical Moist Broadleaf Forests|Are at least 5% of the species records present in this biome?|
|98|Biome|Lake|Are at least 5% of the species records present in this biome?|
|99|Biome|Rock and ice|Are at least 5% of the species records present in this biome?|
|bio1|Climate|Annual Mean Temperature|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio4|Climate|Temperature Seasonality|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio11|Climate|Mean Temperature of Coldest Quarter|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio12|Climate|Annual Precipitation|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio15|Climate|Precipitation Seasonality|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|bio17|Climate|Precipitation of Driest Quarter|The median value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed|
|range_bio1|Climate|Range of annual Mean Temperature|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio4|Climate|Range of temperature Seasonality|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio11|Climate|Range of mean Temperature of Coldest Quarter|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio12|Climate|Range of annual Precipitation|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio15|Climate|Range of precipitation Seasonality|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|range_bio17|Climate|Range of precipitation of Driest Quarter|The range of value of this bioclimatic layer for the occurrence records of a species. Records with NA values removed. Range is the .95-.05 quantile.|
|humanfootprint_1993_1|Human footprint | Human footprint year 1993 lowest impact|The fraction of records in areas of the lowest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_2|Human footprint|Human footprint year 1993 intermediate impact 1|The fraction of records in areas of the second lowest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_3|Human footprint|Human footprint year 1993 intermediate impact 2|The fraction of records in areas of the second highest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_1993_4|Human footprint|Human footprint year 1993 highest impact|The fraction of records in areas of the highest category of human footprint in the year 1993. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_1|Human footprint|Human footprint year 2009 lowest impact|The fraction of records in areas of the lowest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_2|Human footprint|Human footprint year 2009 intermediate impact 1|The fraction of records in areas of the second lowest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_3|Human footprint|Human footprint year 2009 intermediate impact 2|The fraction of records in areas of the second highest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|
|humanfootprint_2009_4|Human footprint|Human footprint year 2009 highest impact|The fraction of records in areas of the highest category of human footprint in the year 2009. Footprint was categorized so that categorize represent roughly quantiles.|

### Labels: Full categories vs Threatened/Not threatened
The `iucnn_prepare_labels` function may accepted any custom labels as long as they are included in the `accepted_labels` option. It also can provide a classification into threatened/non-threatened, via the `level` and `threatened` options. On the broader level the model accuracy is usually significantly higher.

For instance:

```{r, eval = FALSE}
labels_train <- iucnn_prepare_labels(training_labels,
                                     y = features, 
                                     level = "broad")
```

## 2) Model training - NN regression model
### Customizing model parameters
The `iucnn_train_model` function contains various options to customize the neural network, including among other the fraction of validation and test data, the maximum number of epochs, the number of layers and nodes, the activation function , dropout and randomization of the input data. See `?iucnn_train_model` for a comprehensive list of options and their description. By default, `iucnn_train_model` trains a neural network with three hidden layers with 50, 30 and 10 nodes and a sigmoid as activation function. Depending on your dataset different networks may improve performance. For instance, you can set up a different model with 1 hidden layer of 60 nodes, a sigmoid activation function and without using a bias node in the first hidden layer.

```{r}
res_2 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 dropout_rate = 0.3,
 path_to_output= "iucnn_model_2",
 n_layers = "60",
 use_bias = FALSE,
 act_f = "sigmoid")
```

You can compare the validation loss of the models using `res_1$validation_loss` and `res_2$validation_loss`. Model 2 in this case yields a lower validation loss and is therefore preferred. Once you chose the preferred model configuration based on validation loss, we can check test accuracy of best model: `res_2$test_accuracy`. The `iucnn_train_model` function contains various options to adapt the model. See `?iucnn_train_model` for more detail. 

### Changing the modeling algorithm
There are three neural network algorithms implemented in IUCNN. Besides the default classifier approach based on a tensorflow implementation, these are a Bayesian neural network classifier and a regression type neural network.

The Bayesian approach has the advantage that it returns true probabilities for the classification of species into the relative output classes (e.g. 80% probability of a species to be LC). We consider this approach more suitable for classification of species into IUCN categories, than the default option. It will need more time for model training and should best be applied once you have identified the best model parameters using the default approach. You can run a BNN setting the `mode` option of `iucnn_train_model` to `"bnn-class"`.

```{r, eval = FALSE}
res_3 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 path_to_output = "iucnn_model_3",
 mode = 'bnn-class')
```

IUCNN also offers the option to train a NN regression model instead of a classifier. Since the IUCN threat statuses constitute a list of ordinal categories sorted by increasing threat level, we can model the task of estimating these categories as a regression problem. Such a model can be trained with the `iucnn_train_model()` function, specifying to train a regression model by setting `mode = 'nn-reg'`.

```{r}
res_4 <- iucnn_train_model(x = features_train,
 lab = labels_train, 
 path_to_output = "iucnn_model_4",
 mode = 'nn-reg',
 rescale_features = TRUE)
```

### Feature importance
The `iucnn_feature_importance` function can be used to gauge the importance of different feature blocks or individual features for model performance. The function implements the permutation feature importance technique, which randomly reshuffles the values within individual features or blocks of features and evaluate how this randomization affects the models prediction accuracy. If a given feature (or block of features) is important for the models ability to predict, randomizing this feature will lead to a large drop in prediction accuracy. When using `iucnn_feature_importance` with features other than the default, feature blocks can be defined using the `feature_blocks` option.
```{r, eval = TRUE}
fi <- iucnn_feature_importance(x = res_1)
plot(fi)
```

### Model testing
Before training the final model used for predicting the conservation status of not evaluated species, it is recommended to use the `iucnn_modeltest` function for finding the best settings for your model and dataset. This process, often referred to as hyperparameter tuning, is an essential step for building the most suitable model for the prediction task. The `iucnn_modeltest` function allows you to provide any settings for `iucnn_train_model` as vectors, which will lead the function to train a separate model for each provided setting. The function will explore all possible permutations of the provided settings, so that the following command results in 9 different models being trained:

```{r, eval = FALSE}
modeltest_results <- iucnn_modeltest(features,
 labels,
 dropout_rate = c(0.0,0.1,0.3),
 n_layers = c('30','40_20','50_30_10'))
```

The model specifications and settings of each tested model are written to a log-file and can be inspected with the `iucnn_best_model` function, to decide which model settings to pick as best model. Different criteria for picking the best model can be selected, such as best prediction accuracy, best predicted over-all status distribution, lowest weighted mis-classification error, etc.

```{r, eval = FALSE}
best_m <- iucnn_best_model(modeltest_results, criterion='val_acc')
```

After model testing, it is necessary to retrain using the model specifications identified by `iucnn_best_model`. It is possible to take the respective settings directly from the output of the `iucnn_best_model` function, via the `production-model` argument of `iucnn_train_model`.

```{r, eval = FALSE}
# Train the best model on all training data for prediction
m_prod <- iucnn_train_model(train_feat,
                      train_lab,
                      production_model = m_best,
                      overwrite = TRUE)
```

The production model can then be used for the final predictions.

```{r, eval = FALSE}
# Predict RL categories for target species
pred <- iucnn_predict_status(pred_feat,
                      m_prod)
plot(pred)
```


## 3) Status prediction
The `iucnn_predict_status` function offers options to customize the predictions. The most important option in many cases is `target_acc`, which allows to set an overall target-accuracy threshold that the model needs to achieve. This option is only available for nn-class and nn-reg models that were trained using dropout (see help function of `iucnn_train_model` for more explanation), as well as for all bnn-class models. The `target_acc` will be achieved by the model being more selective with making a category call for a given instance. All species that cannot be classified with enough certainty to reach this target accuracy will be classified as NA (Not Assessed).
```{r}
pred_2 <- iucnn_predict_status(x = features_predict, 
 target_acc = 0.7,
 model = res_2)
plot(pred_2)
```

Furthermore, you can turn off the `return_IUCN` option if to return the numerical labels instead of the IUCNN RL category labels.
```{r, eval = FALSE}
pred_3 <- iucnn_predict_status(x = features_predict, 
 model = res_2,
 return_IUCN = FALSE)
```

The output of the `iucnn_predict_status` function is an "iucnn_predictions" object, that contains several output objects. The predicted labels of the individual instances are accessible with `pred_2$class_predictions` and label probabilities estimated by the neural network via `pred_2$mc_dropout_probs` for "nn-class" and "nn-reg" with dropout, or `pred_2$posterior_probs` for "bnn-class". For more detail, the `pred_2$raw_predictions` object contains the individual label probabilities resulting from the softmax output layer in case of "nn-class", or the regressed labels in case of "nn-reg".

## 4) The number of species per category
Another statistic that can be extracted from the "iucnn_predictions" object is the overall category distribution predicted for the given prediction instances. This can be accessed with `pred_2$pred_cat_count`, which shows the distribution of the label predictions, including the count of species that could not be predicted given the chosen `target_acc`.

Another statistic are the `pred2$sampled_cat_freqs` (only available for dropout models and all "bnn-class" models, see above), which show the class distribution as sampled from the `pred_2$mc_dropout_probs` or the `pred_2$posterior_probs` (for "bnn-class" models). The difference between `pred2$sampled_cat_freqs` and `pred_2$mc_dropout_probs`/`pred_2$posterior_probs` is that the former represents the counts of the best labels determined for each instance, whereas the latter represents labels sampled from the predicted label probabilities, which also proportionally samples the labels for a given instance that do receive the maximum label probability. The latter is done repeatedly to include the stochasticity of the random sampling of classes from the given probability vectors. The `pred_2$mc_dropout_probs`/`pred_2$posterior_probs` can be used to plot histograms of the estimates for each class, and can be reported as uncertainty intervals around the number of species in each class for the set of species that were predicted.


# Training and prediction using a convolutional neural network
Instead of the fully connected neural networks presented above, IUCNN also implements a prediction algorithm using Convolutional Neural Networks (CNNs). The input data for prediction are then rasterized per-species grids of occurrence numbers. Since the input features and network structure is different CNNs are implemented in IUCNN with a separate workflow.

## 1. Feature preparation
Features are prepared using the `iucnn_cnn_features` function which will count the number of occurrence records in a custom raster with the same extent as the species occurrences. The function can be simply run on a data frame of training and prediction occurrences separately. Yet, since the extent of these two datasets is likely different, it is recommendable to create a custom raster beforehand. Here, we will use a simple lat/lon raster, but users may provide coordinates and a raster in any suitable coordinate reference system, as long as they agree between occurrence coordinates and raster.

```{r}
# preapre custom raster, you can split this step if training and test occurrences have the same extent
library(terra)
data("training_occ")
data("prediction_occ")

# find the minimum latitude and longitude values for the extent of the raster
min_lon <- min(c(min(training_occ$decimallongitude), 
                 min(prediction_occ$decimallongitude)))
max_lon <- max(c(max(training_occ$decimallongitude), 
                 max(prediction_occ$decimallongitude)))
min_lat <- min(c(min(training_occ$decimallatitude), 
                 min(prediction_occ$decimallatitude)))
max_lat <- max(c(max(training_occ$decimallatitude), 
                 max(prediction_occ$decimallatitude)))
## set the coordinate reference system
ras <- rast(crs = "+proj=longlat +datum=WGS84")
## set raster extent
ext(ras) <- c(min_lon,
              max_lon, 
              min_lat, 
              max_lat)
## set raster resolution
res(ras) <- 1 # the resolution in CRS units, in this case degrees lat/lon


 # Training features
cnn_features <- iucnn_cnn_features(x = training_occ, 
                                   y = ras)

# Prediction features
cnn_features_predict <- iucnn_cnn_features(x = prediction_occ,
                                           y = ras)
 # Training labels
cnn_labels <- iucnn_prepare_labels(x = training_labels,
                                   y = cnn_features)

```


## 2. Model training
```{r, eval = FALSE}
trained_model <- iucnn_cnn_train(cnn_features,
                                cnn_labels,
                                overwrite = TRUE,
                                dropout_rate = 0.1,
                                optimize_for = 'accuracy')

plot(trained_model)
summary(trained_model)
```

## 3. Status prediction
```{r, eval = FALSE}
pred <- iucnn_predict_status(cnn_features_predict,
                             trained_model,
                             target_acc = 0.0
                             )
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_biome_features.R
\name{iucnn_biome_features}
\alias{iucnn_biome_features}
\title{Obtain Biome Features from Occurrence Records}
\usage{
iucnn_biome_features(
  x,
  species = "species",
  lon = "decimallongitude",
  lat = "decimallatitude",
  biome_input = NULL,
  biome_id = "BIOME",
  download_folder = "feature_extraction",
  remove_zeros = FALSE
)
}
\arguments{
\item{x}{a data.frame of species occurrence records including three columns with
species name, longitudinal coordinates and latitudinal coordinates (both decimal).}

\item{species}{a character string. The name of the column with the species names.}

\item{lon}{a character string. The name of the column with the longitude.}

\item{lat}{a character string. The name of the column with the latitude.}

\item{biome_input}{s simple features collection of geometry type polygon,
contain polygons of different biomes.
If NULL, the WWF biome scheme is downloaded from
https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world}

\item{biome_id}{a character string. The name of the column
with the biome names in biome_input.
Default is "BIOME"}

\item{download_folder}{character string. The folder were to save the
data used for feature extraction. Relative to the working directory.
Set to NULL for the working directory}

\item{remove_zeros}{logical. If TRUE biomes without occurrence of
any species are removed from the features.
Default = FALSE}
}
\value{
a data.frame of climatic features
}
\description{
Will code all species in the input file into biomes based on an
intersection of the coordinates with a shape file. The biome scheme can
be user provided or by default will download the WWF biomes.
}
\details{
If biome_input is NULL this will download  the WWF biome scheme from
https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
and save them in the working directory
}
\examples{
dat <- data.frame(species = c("A","b"),
                  decimallongitude = runif (200,10,15),
                  decimallatitude = runif (200,-5,5))

iucnn_biome_features(dat)


}
\seealso{
Other Feature preparation: 
\code{\link{iucnn_bias_features}()},
\code{\link{iucnn_climate_features}()},
\code{\link{iucnn_cnn_features}()},
\code{\link{iucnn_footprint_features}()},
\code{\link{iucnn_geography_features}()},
\code{\link{iucnn_prepare_features}()}
}
\concept{Feature preparation}
\keyword{Feature}
\keyword{preparation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_predict_status.R
\name{iucnn_predict_status}
\alias{iucnn_predict_status}
\title{Predict IUCN Categories from Features}
\usage{
iucnn_predict_status(
  x,
  model,
  target_acc = 0,
  dropout_reps = 100,
  return_IUCN = TRUE,
  return_raw = FALSE
)
}
\arguments{
\item{x}{a data.set, containing a column "species" with the species names, and
subsequent columns with different features,
in the same order as used for \code{\link{iucnn_train_model}}}

\item{model}{the information on the NN model returned by
\code{\link{iucnn_train_model}}}

\item{target_acc}{numerical, 0-1. The target accuracy of the overall model.
Species that cannot be classified with}

\item{dropout_reps}{integer, (default = 100). The number of how often the
predictions are to be repeated (only for dropout models). A value of 100 is
recommended to capture the stochasticity of the predictions, lower values
speed up the prediction time.}

\item{return_IUCN}{logical. If TRUE the predicted labels are translated
into the original labels.
If FALSE numeric labels as used by the model are returned}

\item{return_raw}{logical. If TRUE, the raw predictions of the model will be
returned, which in case of MC-dropout and bnn-class models includes the class
predictions across all dropout prediction reps (or MCMC reps for bnn-class).
Note that setting this to TRUE will result in large output objects that can
fill up the memory allocated for R and cause the program to crash.}
}
\value{
outputs an \code{iucnn_predictions} object containing the predicted
labels for the input species.
}
\description{
Uses a model generated with \code{\link{iucnn_train_model}}
to predict the IUCN status of
Not Evaluated or Data Deficient species based on features, generated
from species occurrence records with \code{\link{iucnn_prepare_features}}.
These features should be of the same type as those used for training the
model.
}
\note{
See \code{vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")} for a
tutorial on how to run IUCNN.
}
\examples{
\dontrun{
data("training_occ") #geographic occurrences of species with IUCN assessment
data("training_labels")# the corresponding IUCN assessments
data("prediction_occ") #occurrences from Not Evaluated species to prdict

# 1. Feature and label preparation
features <- iucnn_prepare_labels(training_occ) # Training features
labels_train <- iucnn_prepare_labels(training_labels) # Training labels
features_predict <- iucnn_prepare_features(prediction_occ) # Prediction features

# 2. Model training
m1 <- iucnn_train_model(x = features, lab = labels_train)

# 3. Prediction
iucnn_predict_status (x = features_predict,
             model = m1)
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_train_model.R
\name{iucnn_train_model}
\alias{iucnn_train_model}
\title{Train an IUCNN Model}
\usage{
iucnn_train_model(
  x,
  lab,
  path_to_output = "iuc_nn_model",
  production_model = NULL,
  mode = "nn-class",
  test_fraction = 0.2,
  cv_fold = 1,
  seed = 1234,
  max_epochs = 1000,
  patience = 200,
  n_layers = "50_30_10",
  use_bias = TRUE,
  balance_classes = FALSE,
  act_f = "auto",
  act_f_out = "auto",
  label_stretch_factor = 1,
  randomize_instances = TRUE,
  dropout_rate = 0,
  mc_dropout = TRUE,
  mc_dropout_reps = 100,
  label_noise_factor = 0,
  rescale_features = FALSE,
  save_model = TRUE,
  overwrite = FALSE,
  verbose = 1
)
}
\arguments{
\item{x}{a data.set, containing a column "species"
with the species names, and
subsequent columns with different features.}

\item{lab}{an object of the class iucnn_labels, as generated by
\code{\link{iucnn_prepare_labels}} containing the labels for all species.}

\item{path_to_output}{character string. The path to the location
where the IUCNN model shall be saved}

\item{production_model}{an object of type iucnn_model (default=NULL).
If an iucnn_model is provided, \code{iucnn_train_model} will read the settings of
this model and reproduce it, but use all available data for training, by
automatically setting the validation set to 0 and cv_fold to 1. This is
recommended before using the model for predicting the IUCN status of
not evaluated species, as it generally improves the prediction
accuracy of the model. Choosing this option will ignore all other provided
settings (below).}

\item{mode}{character string. Choose between the IUCNN models
"nn-class" (default, tensorflow neural network classifier),
"nn-reg" (tensorflow neural network regression), or
"bnn-class" (Bayesian neural network classifier)}

\item{test_fraction}{numeric. The fraction of the input data used as
test set.}

\item{cv_fold}{integer (default=1). When setting cv_fold > 1,
\code{iucnn_train_model} will perform k-fold cross-validation. In this case, the
provided setting for test_fraction will be ignored, as the test
size of each CV-fold is determined by the specified number provided here.}

\item{seed}{integer. Set a starting seed for reproducibility.}

\item{max_epochs}{integer. The maximum number of epochs.}

\item{patience}{integer. Number of epochs with no improvement
after which training will be stopped.}

\item{n_layers}{character string. Define number node per layer by providing a
character string where the number of nodes for each layer are separated by
underscores. E.g. '50_30_10' (default) will train a model with 3 hidden layers with
50, 30, and 10 nodes respectively. Note that the number of nodes in the output
layer is automatically determined based on
the number of unique labels in the training set.}

\item{use_bias}{logical (default=TRUE). Specifies if a bias node is used in
the first hidden layer.}

\item{balance_classes}{logical (default=FALSE). If set to TRUE,
\code{iucnn_train_model} will perform supersampling of the training instances to
account for uneven class distribution in the training data. In case of
training an bnn-class model, choosing this option will add the estimation
of class weights instead, to account for class imbalances.}

\item{act_f}{character string. Specifies the activation
function should be used in the hidden layers.
Available options are: "relu", "tanh", "sigmoid", or "swish" (latter only for
bnn-class). If set to 'auto' (default), \code{iucnn_train_model} will pick a reasonable
default ('relu' for nn-class or nn-reg, and 'swish' for bnn-class).}

\item{act_f_out}{character string. Similar to act_f, this specifies
the activation function for the output
layer. Available options are "softmax" (nn-class, bnn-class), "tanh" (nn-reg),
"sigmoid" (nn-reg), or no activation function "" (nn-reg). When set to "auto"
(default), a suitable output activation function will be chosen based on the
chosen mode ('softmax' for nn-class or bnn-class, 'tanh' for nn-reg).}

\item{label_stretch_factor}{numeric (only for mode nn-reg). The provided
value will be applied as a factor to stretch or compress the labels before
training a regression model. A factor smaller < 1.0 will compress the range
of labels, while a factor > 1 will stretch the range.}

\item{randomize_instances}{logical. When set to TRUE (default) the
instances will be shuffled before training (recommended).}

\item{dropout_rate}{numeric. This will randomly turn off the specified
fraction of nodes of the neural network during each epoch of training
making the NN more stable and less reliant on individual nodes/weights, which
can prevent over-fitting (only available for modes nn-class and nn-reg).
See mc_dropout setting explained below if dropout shall also be applied to the
predictions.}

\item{mc_dropout}{logical. If set to TRUE, the predictions (including the
validation accuracy) based on a model trained with a dropout fraction > 0
will reflect the stochasticity introduced by the dropout method (MC dropout
predictions). This is e.g. required when wanting to predict with a specified
accuracy threshold (see target_acc option in \code{\link{iucnn_predict_status}}).
This option is activated by default when chosing a dropout_rate > 0, unless
it is manually set to FALSE here.}

\item{mc_dropout_reps}{integer. The number of MC iterations to run when
predicting validation accuracy and calculating the accuracy-threshold
table required for making predictions with an accuracy threshold.
The default of 100 is usually sufficient, larger values will lead to longer
computation times, particularly during model testing with cross-validation.}

\item{label_noise_factor}{numeric (only for mode nn-reg). Add specified amount
of random noise to the input labels to give the categorical labels a more
continuous spread before training the regression model. E.g. a value of 0.2
will redraw a label of a species categorized as Vulnerable (class=2) randomly
between 1.8 and 2.2, based on a uniform probability distribution.}

\item{rescale_features}{logical. Set to TRUE if all feature values shall
be rescaled to values between 0 and 1 prior to training (default=FALSE).}

\item{save_model}{logical. If TRUE the model is saved to disk.}

\item{overwrite}{logical. If TRUE existing models are
overwritten. Default is set to FALSE.}

\item{verbose}{Default 0, set to 1 for \code{iucnn_train_model} to print
additional info to the screen while training.}
}
\value{
outputs an \code{iucnn_model} object which can be used in
\code{\link{iucnn_predict_status}} for predicting the conservation status
of not evaluated species.
}
\description{
Trains an IUCNN model based on a data.frame of features for a set of species,
generated by \code{\link{iucnn_prepare_features}},
and the corresponding IUCN classes formatted as a iucnn_labels object
with \code{\link{iucnn_prepare_labels}}. Note
that NAs are not allowed in the features, and taxa with NAs will
automatically be removed! Taxa, for which information is only present in one
of the two input objects will be removed as well.
}
\note{
See \code{vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")}
for a tutorial on how to run IUCNN.
}
\examples{
\dontrun{
data("training_occ") #geographic occurrences of species with IUCN assessment
data("training_labels")# the corresponding IUCN assessments
data("prediction_occ") #occurrences from Not Evaluated species to prdict

# 1. Feature and label preparation
features <- iucnn_prepare_features(training_occ) # Training features
labels_train <- iucnn_prepare_labels(training_labels) # Training labels
features_predict <- iucnn_prepare_features(prediction_occ) # Prediction features

# 2. Model training
m1 <- iucnn_train_model(x = features, lab = labels_train)

summary(m1)
plot(m1)
}


}
\keyword{Training}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_prepare_labels.R
\name{iucnn_prepare_labels}
\alias{iucnn_prepare_labels}
\title{Format IUCN Red List categories for IUCNN}
\usage{
iucnn_prepare_labels(
  x,
  y,
  species = "species",
  labels = "labels",
  accepted_labels = c("LC", "NT", "VU", "EN", "CR"),
  level = "detail",
  threatened = c("CR", "EN", "VU")
)
}
\arguments{
\item{x}{a data.frame or a list. If a data.frame,
two columns with the species names and IUCN categories
respectively. The column names are defined by the
species and labels arguments. If a list, expecting
the format as returned by \link[rredlist]{rl_search}.}

\item{y}{object of class \code{iucnn-features} or \code{iucnn_cnn_features}.
Ensures that the species in the return value are in the same order as in y.}

\item{species}{a character string. The name of the
column with the species names in x.}

\item{labels}{a character string. The name of the
column with the labels (assessment categories) in x.}

\item{accepted_labels}{a character string. The labels
to be converted in to numeric values.
Entries with labels not mentioned (e.g. "DD") will be removed.
The numeric labels returned by the
function will correspond to the order in this argument.
For instance with default settings, LC -> 0, CR -> 4.}

\item{level}{a character string. The level of output
level detail. If "detail"
full IUCN categories, if "broad" then
0 = Not threatened, and 1 = Threatened.}

\item{threatened}{a character string. Only if level=="broad",
Which labels to consider threatened.}
}
\value{
a data.frame with species names and numeric labels
}
\description{
Converting IUCN category labels into numeric categories required by \code{\link{iucnn_train_model}}.
}
\note{
See \code{vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")}
for a tutorial on how to run IUCNN.
}
\examples{
dat <- data.frame(species = c("A","B"),
                  decimallongitude = runif (200,10,15),
                  decimallatitude = runif (200,-5,5))

labs <- data.frame(species = c("A","B"),
                   labels = c("CR", "LC"))

features <- iucnn_prepare_features(dat,
                                   type = "geographic")

iucnn_prepare_labels(x = labs,
                     y = features)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{prediction_occ}
\alias{prediction_occ}
\title{Geographic Occurrence Records for Not Evaluated Orchids}
\format{
A data frame with 14,900 rows and 3 variables:
\describe{
\item{species}{The canonical species name}
\item{decimallongitude}{longitudinal coordinates}
\item{decimallatitude}{latitudinal coordinates}
}
}
\source{
\url{https://www.gbif.org/}
}
\usage{
prediction_occ
}
\description{
A dataset containing geo-referenced occurrences of 100 Orchid species without
existing IUCN Red List assessment ("Not Evaluated"). This is example data
to predict the IUCN status.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_cnn_features.R
\name{iucnn_cnn_features}
\alias{iucnn_cnn_features}
\title{Prepare Features for a CNN model}
\usage{
iucnn_cnn_features(
  x,
  y = NULL,
  species = "species",
  lon = "decimallongitude",
  lat = "decimallatitude",
  crs_x = "+proj=longlat +datum=WGS84",
  res_y = 1
)
}
\arguments{
\item{x}{a data.frame with at least three columns containing taxon name, decimal longitude and latitude values.}

\item{y}{a raster as reference to count the number of occurrence records in.
Can be of any resolution and CRS but the coordinates in x need to be in the same CRS.}

\item{species}{a character string. The name of the column with the species names.}

\item{lon}{a character string. The name of the column with the longitude.}

\item{lat}{a character string. The name of the column with the latitude.}

\item{crs_x}{a proj4string specifying the Coordinate Reference system of the coordinates in x.
Default is to lat/lon WGS84.}

\item{res_y}{numeric. The resolution for the raster in decimal degrees.
Only relevant if y is not provided.}
}
\value{
a list of matrices, one for each input species, where the cells represent the number
of occurrence records in this cell as input for the \dQuote{cnn} class of \code{\link{iucnn_train_model}}.
}
\description{
Converts a data.frame of species occurrences into input features to use
a convolutional neural network to approximate species extinction risk.
}
\details{
If y is not provided, assumes a lat/lon grid with extent equal to the respective minimum and maximum in x
}
\examples{
dat <- data.frame(species = c("A","B"),
                  decimallongitude = runif (200,10,15),
                  decimallatitude = runif (200,-5,5))

iucnn_cnn_features(dat)


}
\seealso{
Other Feature preparation: 
\code{\link{iucnn_bias_features}()},
\code{\link{iucnn_biome_features}()},
\code{\link{iucnn_climate_features}()},
\code{\link{iucnn_footprint_features}()},
\code{\link{iucnn_geography_features}()},
\code{\link{iucnn_prepare_features}()}
}
\concept{Feature preparation}
\keyword{Feature}
\keyword{preparation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_climate_features.R
\name{iucnn_climate_features}
\alias{iucnn_climate_features}
\title{Extract Climatic Features from Occurrence Records}
\usage{
iucnn_climate_features(
  x,
  climate_input = NULL,
  species = "species",
  lon = "decimallongitude",
  lat = "decimallatitude",
  rescale = TRUE,
  res = 10,
  type = "selected",
  download_folder = "feature_extraction"
)
}
\arguments{
\item{x}{a data.frame of species occurrence records including three columns with
species name, longitudinal coordinates and latitudinal coordinates (both decimal).}

\item{climate_input}{a raster or rasterStack with climate data.
Optional. If not provided,
the 19 bioclim variables from www.worldclim.org are used as default.}

\item{species}{a character string. The name of the column with the species names.}

\item{lon}{a character string. The name of the column with the longitude.}

\item{lat}{a character string. The name of the column with the latitude.}

\item{rescale}{logical. If TRUE, the features are rescaled.
This is recommended to run IUCNN, and the default. If FALSE, raw (human readable)
feature values are returned.}

\item{res}{numeric. The resolution of the default climate rasters.
One of 2.5, 5, or 10. Only relevant if
climate_input is NULL}

\item{type}{character string. A selection of which variables to return.
If "all" all 19 bioclim variables
if "selected" only Annual Mean Temperature, Temperature Seasonality,
Mean temperature of the Coldest Quarter,
Annual Precipitation, Precipitation seasonality and
Precipitation of the Driest Quarter are returned}

\item{download_folder}{character string. The folder were to save the
data used for feature extraction. Relative to the working directory.
Set to NULL for the working directory}
}
\value{
a data.frame of climatic features
}
\description{
Extract median and range (95\% - 5\% quantiles) climate
features based on a table of species occurrence coordinates.
If no climate data is supplied via the input.climate
argument, 19 bioclim variables are
downloaded from  www.worldclim.org.
Rescaling is only done for these default variables.
}
\details{
All climate variables are summarized  to the species median.
}
\examples{
dat <- data.frame(species = c("A","B"),
                  decimallongitude = runif (200,10,15),
                  decimallatitude = runif (200,-5,5))

iucnn_climate_features(dat)


}
\seealso{
Other Feature preparation: 
\code{\link{iucnn_bias_features}()},
\code{\link{iucnn_biome_features}()},
\code{\link{iucnn_cnn_features}()},
\code{\link{iucnn_footprint_features}()},
\code{\link{iucnn_geography_features}()},
\code{\link{iucnn_prepare_features}()}
}
\concept{Feature preparation}
\keyword{Feature}
\keyword{preparation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_bias_features.R
\name{iucnn_bias_features}
\alias{iucnn_bias_features}
\title{Extract Bias Features from Occurrence Records}
\usage{
iucnn_bias_features(
  x,
  species = "species",
  lon = "decimallongitude",
  lat = "decimallatitude",
  res = 0.5,
  ras = NULL,
  plot = TRUE
)
}
\arguments{
\item{x}{a data.frame of species occurrence records including three columns with
species name, longitudinal coordinates and latitudinal coordinates (both decimal).}

\item{species}{a character string. The name of the column with the species names.}

\item{lon}{a character string. The name of the column with the longitude.}

\item{lat}{a character string. The name of the column with the latitude.}

\item{res}{numeric. The resolution of the default resolution to calculate sampling bias. In decimal degrees.}

\item{ras}{a raster object. Alternative to res, a sample raster to calculate sampling bias. Needs to use the same CRS as
the coordinates in x.}

\item{plot}{logical. Should the results of the sampbias analysis be plotted for diagnostics?}
}
\value{
a data.frame of bias features
}
\description{
Use the sampbias method to assess the geographic sampling bias at the locations where a species is collected and the range of
sampling bias for all records per species.Values summarized per species are the median and the 0.05 to 0.95 percentiles.
}
\details{
See the ?sampbias::calculate_bias for details.
}
\examples{
\dontrun{

iucnn_bias_features(dat)
}


}
\seealso{
Other Feature preparation: 
\code{\link{iucnn_biome_features}()},
\code{\link{iucnn_climate_features}()},
\code{\link{iucnn_cnn_features}()},
\code{\link{iucnn_footprint_features}()},
\code{\link{iucnn_geography_features}()},
\code{\link{iucnn_prepare_features}()}
}
\concept{Feature preparation}
\keyword{Feature}
\keyword{preparation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{training_occ}
\alias{training_occ}
\title{Geographic Occurrence Records for Orchids with IUCN assessment}
\format{
A data frame with 125,412 rows and 3 variables:
\describe{
\item{species}{The canonical species name}
\item{decimallongitude}{longitudinal coordinates}
\item{decimallatitude}{latitudinal coordinates}
}
}
\source{
\url{https://www.gbif.org/}
}
\usage{
training_occ
}
\description{
A dataset containing geo-referenced occurrences of 884 Orchid species with
existing IUCN Red List assessment ("CR", "EN", "VU", "NT", "LC").
This is example data to train an IUCNN model.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_geography_features.R
\name{iucnn_geography_features}
\alias{iucnn_geography_features}
\title{Extract Geographic Features from Occurrence Records}
\usage{
iucnn_geography_features(
  x,
  species = "species",
  lon = "decimallongitude",
  lat = "decimallatitude",
  rescale = TRUE
)
}
\arguments{
\item{x}{a data.frame of species occurrence records including three columns with
species name, longitudinal coordinates and latitudinal coordinates (both decimal).}

\item{species}{a character string. The name of the column with the species names.}

\item{lon}{a character string. The name of the column with the longitude.}

\item{lat}{a character string. The name of the column with the latitude.}

\item{rescale}{logical. If TRUE, the geographic features are rescaled.
This is recommended to run IUCNN, and the default. If FALSE, raw (human readable)
feature values are returned.}
}
\value{
a data.frame of geographic features
}
\description{
Calculates the number of occurrences, number of unique occurrences,
mean latitude, mean longitude, latitudinal range, longitudinal range,
eoo, aoo and hemisphere as input features for IUCNN
from a list of species occurrences.
}
\details{
Coordinate ranges are 90\% quantiles, for species with
less than three occurrences EOO is set to AOO.
}
\examples{
dat <- data.frame(species = c("A","B"),
                  decimallongitude = runif (200,10,15),
                  decimallatitude = runif (200,-5,5))

iucnn_geography_features(dat)

}
\seealso{
Other Feature preparation: 
\code{\link{iucnn_bias_features}()},
\code{\link{iucnn_biome_features}()},
\code{\link{iucnn_climate_features}()},
\code{\link{iucnn_cnn_features}()},
\code{\link{iucnn_footprint_features}()},
\code{\link{iucnn_prepare_features}()}
}
\concept{Feature preparation}
\keyword{Feature}
\keyword{preparation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_footprint_features.R
\name{iucnn_footprint_features}
\alias{iucnn_footprint_features}
\title{Extract Human Footprint Index Features from Occurrence Records}
\source{
https://wcshumanfootprint.org/
}
\usage{
iucnn_footprint_features(
  x,
  footp_input = NULL,
  species = "species",
  lon = "decimallongitude",
  lat = "decimallatitude",
  rescale = TRUE,
  year = c(1993, 2009),
  download_folder = "feature_extraction",
  breaks = c(0, 0.81, 1.6, 2.3, 100)
)
}
\arguments{
\item{x}{a data.frame of species occurrence records including three columns with
species name, longitudinal coordinates and latitudinal coordinates (both decimal).}

\item{footp_input}{an object of the class raster or RasterStack
with values for the human footprint index.
If a RasterStack, different layers are interpreted as different time-slices.}

\item{species}{a character string. The name of the column with the species names.}

\item{lon}{a character string. The name of the column with the longitude.}

\item{lat}{a character string. The name of the column with the latitude.}

\item{rescale}{logical. If TRUE, the values are rescaled using
natural logarithm transformation. If FALSE,
remember to change the breaks argument.}

\item{year}{numeric. The years for which to obtain the human footprint index.
The default is to the two layers available. Can be a either year, in case only
one slice is desired. Other time slices are currently not supported,}

\item{download_folder}{character string. The folder were to save the
data used for feature extraction. Relative to the working directory.
Set to NULL for the working directory}

\item{breaks}{numerical. The breaks to bin the human footprint index
for the final features. The defaults are
empirical values for the global footprint and rescale=TRUE.
For custom values ensure that they
cover the whole value range and are adapted to the value of rescale.}
}
\value{
a data.frame of human footprint features
}
\description{
Bins the human footprint index into a set of bins and the
fraction of occurrence records
of a species in each bin are the features.
By default the human footprint index is downloaded from
https://wcshumanfootprint.org/. THIS FUNCTION WILL DOWNLOAD DATA FROM
THE INTERNET AND SAVE IT TO THE  WORKING DIRECTORY. The data files
are >200 MB each and downloading may
take some time on first execution.
}
\details{
By default four categories of increasing human footprint index
( 1 = lowest, 4 = highest) are selected and rescaled.
}
\examples{
dat <- data.frame(species = c("A","B"),
                  decimallongitude = runif (200,10,15),
                  decimallatitude = runif (200,-5,5))

iucnn_footprint_features(dat)


}
\seealso{
Other Feature preparation: 
\code{\link{iucnn_bias_features}()},
\code{\link{iucnn_biome_features}()},
\code{\link{iucnn_climate_features}()},
\code{\link{iucnn_cnn_features}()},
\code{\link{iucnn_geography_features}()},
\code{\link{iucnn_prepare_features}()}
}
\concept{Feature preparation}
\keyword{Feature}
\keyword{preparation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_prepare_features.R
\name{iucnn_prepare_features}
\alias{iucnn_prepare_features}
\title{Prepare features for an IUCNN model}
\usage{
iucnn_prepare_features(
  x,
  species = "species",
  lon = "decimallongitude",
  lat = "decimallatitude",
  type = c("geographic", "biomes", "climate", "humanfootprint"),
  download_folder = "feature_extraction"
)
}
\arguments{
\item{x}{a data.frame of species occurrence records including three columns with
species name, longitudinal coordinates and latitudinal coordinates (both decimal).}

\item{species}{a character string. The name of the column with the species names.}

\item{lon}{a character string. The name of the column with the longitude.}

\item{lat}{a character string. The name of the column with the latitude.}

\item{type}{character. The type of features to calculate. Possible options are
\dQuote{geographic}, \dQuote{biome}, \dQuote{climate},
\dQuote{human footprint}.}

\item{download_folder}{character string. The folder were to save the
data used for feature extraction. Relative to the working directory.
Set to NULL for the working directory}
}
\value{
a data.frame of features
}
\description{
A wrapper function to prepare all default features included in IUCNN:
geographic, biomes, climate, human footprint.
If desired, bias features need to be calculated separately with ft_bias.
For more control over feature preparation, you can use the
\code{\link{iucnn_geography_features}}, \code{\link{iucnn_biome_features}}, \code{\link{iucnn_climate_features}},
\code{\link{iucnn_footprint_features}} functions.
}
\details{
Without internet access, only geographic features are calculated,
}
\examples{
dat <- data.frame(species = c("A","B"),
                  decimallongitude = runif (200,10,15),
                  decimallatitude = runif (200,-5,5))

iucnn_prepare_features(dat)

}
\seealso{
Other Feature preparation: 
\code{\link{iucnn_bias_features}()},
\code{\link{iucnn_biome_features}()},
\code{\link{iucnn_climate_features}()},
\code{\link{iucnn_cnn_features}()},
\code{\link{iucnn_footprint_features}()},
\code{\link{iucnn_geography_features}()}
}
\concept{Feature preparation}
\keyword{Feature}
\keyword{preparation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_modeltest.R
\name{iucnn_modeltest}
\alias{iucnn_modeltest}
\title{Model-Testing IUCNN Models using Cross-Validation (Hyperparameter-Tuning)}
\usage{
iucnn_modeltest(
  x,
  lab,
  logfile = "model_testing_logfile.txt",
  model_outpath = "modeltest",
  mode = "nn-class",
  cv_fold = 5,
  test_fraction = 0,
  n_layers = c("50_30_10", "30"),
  dropout_rate = c(0, 0.1, 0.3),
  use_bias = TRUE,
  balance_classes = FALSE,
  seed = 1234,
  label_stretch_factor = 1,
  label_noise_factor = 0,
  act_f = "relu",
  act_f_out = "auto",
  max_epochs = 5000,
  patience = 200,
  mc_dropout = TRUE,
  mc_dropout_reps = 100,
  randomize_instances = TRUE,
  rescale_features = FALSE,
  init_logfile = TRUE,
  recycle_settings = FALSE
)
}
\arguments{
\item{x}{a data.set, containing a column "species"
with the species names, and
subsequent columns with different features.}

\item{lab}{an object of the class iucnn_labels, as generated by
\code{\link{iucnn_prepare_labels}} containing the labels for all species.}

\item{logfile}{character string. Define the filepath/name for the output
log-file.}

\item{model_outpath}{the path where to save the results on disk}

\item{mode}{character string. Choose between the IUCNN models
"nn-class" (default, tensorflow neural network classifier),
"nn-reg" (tensorflow neural network regression), or
"bnn-class" (Bayesian neural network classifier)}

\item{cv_fold}{integer (default=1). When setting cv_fold > 1,
\code{iucnn_train_model} will perform k-fold cross-validation. In this case, the
provided setting for test_fraction will be ignored, as the test
size of each CV-fold is determined by the specified number provided here.}

\item{test_fraction}{numeric. The fraction of the input data used as
test set.}

\item{n_layers}{character string. Define number node per layer by providing a
character string where the number of nodes for each layer are separated by
underscores. E.g. '50_30_10' (default) will train a model with 3 hidden layers with
50, 30, and 10 nodes respectively. Note that the number of nodes in the output
layer is automatically determined based on
the number of unique labels in the training set.}

\item{dropout_rate}{numeric. This will randomly turn off the specified
fraction of nodes of the neural network during each epoch of training
making the NN more stable and less reliant on individual nodes/weights, which
can prevent over-fitting (only available for modes nn-class and nn-reg).
See mc_dropout setting explained below if dropout shall also be applied to the
predictions.}

\item{use_bias}{logical (default=TRUE). Specifies if a bias node is used in
the first hidden layer.}

\item{balance_classes}{logical (default=FALSE). If set to TRUE,
\code{iucnn_train_model} will perform supersampling of the training instances to
account for uneven class distribution in the training data. In case of
training an bnn-class model, choosing this option will add the estimation
of class weights instead, to account for class imbalances.}

\item{seed}{integer. Set a starting seed for reproducibility.}

\item{label_stretch_factor}{numeric (only for mode nn-reg). The provided
value will be applied as a factor to stretch or compress the labels before
training a regression model. A factor smaller < 1.0 will compress the range
of labels, while a factor > 1 will stretch the range.}

\item{label_noise_factor}{numeric (only for mode nn-reg). Add specified amount
of random noise to the input labels to give the categorical labels a more
continuous spread before training the regression model. E.g. a value of 0.2
will redraw a label of a species categorized as Vulnerable (class=2) randomly
between 1.8 and 2.2, based on a uniform probability distribution.}

\item{act_f}{character string. Specifies the activation
function should be used in the hidden layers.
Available options are: "relu", "tanh", "sigmoid", or "swish" (latter only for
bnn-class). If set to 'auto' (default), \code{iucnn_train_model} will pick a reasonable
default ('relu' for nn-class or nn-reg, and 'swish' for bnn-class).}

\item{act_f_out}{character string. Similar to act_f, this specifies
the activation function for the output
layer. Available options are "softmax" (nn-class, bnn-class), "tanh" (nn-reg),
"sigmoid" (nn-reg), or no activation function "" (nn-reg). When set to "auto"
(default), a suitable output activation function will be chosen based on the
chosen mode ('softmax' for nn-class or bnn-class, 'tanh' for nn-reg).}

\item{max_epochs}{integer. The maximum number of epochs.}

\item{patience}{integer. Number of epochs with no improvement
after which training will be stopped.}

\item{mc_dropout}{logical. If set to TRUE, the predictions (including the
validation accuracy) based on a model trained with a dropout fraction > 0
will reflect the stochasticity introduced by the dropout method (MC dropout
predictions). This is e.g. required when wanting to predict with a specified
accuracy threshold (see target_acc option in \code{\link{iucnn_predict_status}}).
This option is activated by default when chosing a dropout_rate > 0, unless
it is manually set to FALSE here.}

\item{mc_dropout_reps}{integer. The number of MC iterations to run when
predicting validation accuracy and calculating the accuracy-threshold
table required for making predictions with an accuracy threshold.
The default of 100 is usually sufficient, larger values will lead to longer
computation times, particularly during model testing with cross-validation.}

\item{randomize_instances}{logical. When set to TRUE (default) the
instances will be shuffled before training (recommended).}

\item{rescale_features}{logical. Set to TRUE if all feature values shall
be rescaled to values between 0 and 1 prior to training (default=FALSE).}

\item{init_logfile}{logical (default=TRUE). If set to TRUE,
\code{modeltest_iucnn} will attempt to initiate a new log-file under the
provided path, possibly overwriting already existing model-testing results
stored in the same location. Set to FALSE if instead you want to append to
an already existing log-file.}

\item{recycle_settings}{logical (default=FALSE). If set to TRUE,
\code{iucnn_modeltest} will read the log-file stored at the path specified
under the "logfile" argument and run model-testing for the input features
and labels using the same models stored in that file. This setting can be
useful when e.g. wanting to test the same models for different sets of input
data.}
}
\value{
outputs a data.frame object containing stats and settings of all
tested models.
}
\description{
Takes as input features produced with \code{\link{iucnn_prepare_features}}
and labels produced with \code{\link{iucnn_prepare_labels}}, as well as a path to a
log-file where results for each tested model will be stored. All available
options are identical to the \code{\link{iucnn_train_model}} function and can be
provided as vectors, e.g. \code{dropout_rate = c(0.0,0.1,0.3)} and
\code{n_layers = c('30','40_20','50_30_10')}. \code{iucnn_modeltest} will
then iterate through all possible permutations of the provided hyperparameter
settings, train a separate model for each hyperparameter combination, and
store the results in the provided log-file.
}
\note{
See \code{vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")}
for a tutorial on how to run IUCNN.
}
\examples{
\dontrun{
# Model-testing
logfile = paste0("model_testing_results.txt")
model_testing_results = iucnn_modeltest(features,
                                       labels,
                                       logfile,
                                       model_outpath = 'iucnn_modeltest',
                                       mode = 'nn-class',
                                       seed = 1234,
                                       dropout_rate = c(0.0,0.1,0.3),
                                       n_layers = c('30','40_20','50_30_10'),
                                       cv_fold = 5,
                                       init_logfile = TRUE)
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_best_model.R
\name{iucnn_best_model}
\alias{iucnn_best_model}
\title{Select the Best Model After Model-testing}
\usage{
iucnn_best_model(x, criterion = "val_acc", require_dropout = FALSE)
}
\arguments{
\item{x}{a data.frame of model-testing results as produced
by \code{\link{iucnn_modeltest}}.}

\item{criterion}{name the criterion to rank models by (default="val_acc").
Valid options are
"val_acc","val_loss","weighted_error", or "total_class_matches"
(see details below):
\itemize{
\item val_acc: highest validation accuracy
\item val_loss: lowest validation loss
\item weighted_error: lowest weighted error, e.g. an LC species misclassified as
CR has a weighted error of 4-0 = 4, while an LC species
misclassified as NT has a weighted error of 1-0 = 1.
These error scores are summed across all validation
predictions
\item total_class_matches: picks the model that best reproduces the class
distribution in the validation data. When picking
this criterion it is not considered whether or not
individual instances are predicted correctly, but
instead it only looks at the overall class distribution
in the predicted data.
}}

\item{require_dropout}{logical (default=FALSE). If set to TRUE, the best model
that contains a dropout rate of > 0 will be picked, even if other non-dropout
models scored higher given the chosen criterion. Dropout models are required
for certain functionalities within IUCNN, such as e.g. choosing a target
accuracy when using predict_iucnn.}
}
\value{
outputs an \code{iucnn_model} object containing all
information about the best model.
}
\description{
Uses a data-frame of model-testing results generated with
\code{\link{iucnn_modeltest}} as input, and finds the best model
based on the chosen criterion.
}
\note{
See \code{vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")}
for a tutorial on how to run IUCNN.
}
\examples{
\dontrun{
# Model-testing
logfile = paste0("model_testing_results.txt")
model_testing_results = modeltest_iucnn(features,
                                       labels,
                                       logfile,
                                       model_outpath = 'iucnn_modeltest',
                                       mode = 'nn-class',
                                       seed = 1234,
                                       dropout_rate = c(0.0,0.1,0.3),
                                       n_layers = c('30','40_20','50_30_10'),
                                       cv_fold = 5,
                                       init_logfile = TRUE)

# Selecting best model based on chosen criterion
best_iucnn_model = bestmodel_iucnn(model_testing_results,
                                   criterion = 'val_acc',
                                   require_dropout = TRUE)
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_cnn_train.R
\name{iucnn_cnn_train}
\alias{iucnn_cnn_train}
\title{Train a CNN model}
\usage{
iucnn_cnn_train(
  x,
  lab,
  path_to_output = "iuc_nn_model",
  production_model = NULL,
  cv_fold = 1,
  test_fraction = 0.2,
  seed = 1234,
  max_epochs = 100,
  patience = 20,
  randomize_instances = TRUE,
  balance_classes = TRUE,
  dropout_rate = 0,
  mc_dropout_reps = 100,
  optimize_for = "loss",
  pooling_strategy = "average",
  save_model = TRUE,
  overwrite = FALSE,
  verbose = 0
)
}
\arguments{
\item{x}{a list of matrices containing the occurrence counts across a spatial
grid for a set of species.}

\item{lab}{an object of the class iucnn_labels, as generated by
\code{\link{iucnn_prepare_labels}} containing the labels for all species.}

\item{path_to_output}{character string. The path to the location
where the IUCNN model shall be saved}

\item{production_model}{an object of type iucnn_model (default=NULL).
If an iucnn_model is provided, \code{iucnn_cnn_train} will read the settings of
this model and reproduce it, but use all available data for training, by
automatically setting the validation set to 0 and cv_fold to 1. This is
recommended before using the model for predicting the IUCN status of
not evaluated species, as it generally improves the prediction
accuracy of the model. Choosing this option will ignore all other provided
settings below.}

\item{cv_fold}{integer (default=1). When setting cv_fold > 1,
\code{iucnn_cnn_train} will perform k-fold cross-validation. In this case, the
provided setting for test_fraction will be ignored, as the test
size of each CV-fold is determined by the specified number provided here.}

\item{test_fraction}{numeric. The fraction of the input data used as
test set.}

\item{seed}{integer. Set a starting seed for reproducibility.}

\item{max_epochs}{integer. The maximum number of epochs.}

\item{patience}{integer. Number of epochs with no improvement
after which training will be stopped.}

\item{randomize_instances}{logical (default=TRUE). When set to TRUE (default)
the instances will be shuffled before training (recommended).}

\item{balance_classes}{logical (default=FALSE). If set to TRUE,
\code{iucnn_cnn_train} will perform supersampling of the training instances to
account for uneven class distribution in the training data.}

\item{dropout_rate}{numeric. This will randomly turn off the specified
fraction of nodes of the neural network during each epoch of training
making the NN more stable and less reliant on individual nodes/weights, which
can prevent over-fitting (only available for modes nn-class and nn-reg).
See mc_dropout setting explained below if dropout shall also be applied to the
predictions. For models trained with a dropout fraction > 0, the predictions
(including the validation accuracy)
will reflect the stochasticity introduced by the dropout method (MC dropout
predictions). This is e.g. required when wanting to predict with a specified
accuracy threshold (see target_acc option in
\code{\link{iucnn_predict_status}}).}

\item{mc_dropout_reps}{integer. The number of MC iterations to run when
predicting validation accuracy and calculating the accuracy-threshold
table required for making predictions with an accuracy threshold.
The default of 100 is usually sufficient, larger values will lead to longer
computation times, particularly during model testing with cross-validation.}

\item{optimize_for}{string. Default is "loss", which will train the model
until optimal validation set loss is reached. Set to "accuracy" if you want
to optimize for maximum validation accuracy instead.}

\item{pooling_strategy}{string. Pooling strategy after first convolutional
layer. Choose between  "average" (default) and "max".}

\item{save_model}{logical. If TRUE the model is saved to disk.}

\item{overwrite}{logical. If TRUE existing models are
overwritten. Default is set to FALSE.}

\item{verbose}{Default 0, set to 1 for \code{iucnn_cnn_train} to print
additional info to the screen while training.}
}
\value{
outputs an \code{iucnn_model} object which can be used in
\code{\link{iucnn_predict_status}} for predicting the conservation status
of not evaluated species.
}
\description{
Trains an CNN model based on a list of matrices with occurrence counts for a
set of species,
generated by \code{\link{iucnn_cnn_features}},
and the corresponding IUCN classes formatted as a iucnn_labels object
with \code{\link{iucnn_prepare_labels}}. Note
that taxa for which information is only present in one
of the two input objects will be removed from further processing.
}
\note{
See \code{vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")}
for a tutorial on how to run IUCNN.
}
\examples{
\dontrun{
trained_model <- iucnn_cnn_train(cnn_training_features,
                                cnn_labels,
                                cv_fold=5,
                                overwrite = TRUE,
                                dropout=0.1)
summary(trained_model)
}

}
\keyword{Training}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iucnn_feature_importance.R
\name{iucnn_feature_importance}
\alias{iucnn_feature_importance}
\title{Evaluate relative importance of training features}
\usage{
iucnn_feature_importance(
  x,
  feature_blocks = list(),
  n_permutations = 100,
  provide_indices = FALSE,
  verbose = FALSE,
  unlink_features_within_block = TRUE
)
}
\arguments{
\item{x}{iucnn_model object, as produced as output
when running \code{\link{iucnn_train_model}}}

\item{feature_blocks}{a list. Default behavior is to
group the features into geographic, climatic,
biome, and human footprint features. Provide custom
list of feature names or indices to define other
feature blocks. If feature
indices are provided as in this example, turn provide_indices flag to TRUE.}

\item{n_permutations}{an integer. Defines how many
iterations of shuffling feature values and
predicting the resulting accuracy are being executed.
The mean and standard deviation of the
delta accuracy are being summarized from these permutations.}

\item{provide_indices}{logical. Set to TRUE if custom \code{feature_blocks}
are provided as indices. Default is FALSE.}

\item{verbose}{logical. Set to TRUE to print screen output while calculating
feature importance. Default is FALSE.}

\item{unlink_features_within_block}{logical. If TRUE, the features within each
defined block are shuffled independently.
If FALSE, each feature column within a block is resorted in the same manner. Default is TRUE}
}
\value{
a data.frame with the relative importance of each feature block (see delta_acc_mean column).
}
\description{
Uses a model generated with \code{\link{iucnn_train_model}}
to evaluate how much each feature or
group of features contributes to the accuracy of
the test set predictions. The function
implements the concept of permutation feature importance,
in which the values in a given
feature column of the test set are shuffled randomly
among all samples. Then the feature
data manipulated in this manner are used to predict
labels for the test set and the accuracy
is compared to that of the original feature data.
The difference (delta accuracy) can be
interpreted as a measure of how important a
given feature or group of features is for the
trained NN to make accurate predictions.
}
\details{
By default this function groups the features
into geographic, climatic, biome, and human
footprint features and determines the importance
of each of these blocks of features. The
feature blocks can be manually defined using the feature_blocks argument.
}
\note{
See \code{vignette("Approximate_IUCN_Red_List_assessments_with_IUCNN")} for a
tutorial on how to run IUCNN.
}
\examples{
\dontrun{
data("training_occ")
data("training_labels")

train_feat <- iucnn_prepare_features(training_occ)
labels_train <- iucnn_prepare_labels(training_labels,
                           level = 'detail')

train_output <- iucnn_train_model(x = train_feat,
                          lab = labels_train,
                          patience = 10)


imp_def <- iucnn_feature_importance(x = train_output)
imp_cust <- iucnn_feature_importance(x = train_output,
                              feature_blocks = list(block1 = c(1,2,3,4),
                                                    block2 = c(5,6,7,8)),
                              provide_indices = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{training_labels}
\alias{training_labels}
\title{IUCN threat categories for 884 orchid species}
\format{
A data frame with 889 rows and 2 variables:
\describe{
\item{species}{The canonical species name}
\item{labels}{The IUCN conservation assessment,
converted to numerical for use with IUCNN.
0 =  Critically Endangered (CR),
1 = Endangered (EN),
2 = Vulnerable (VU),
3 = Near Threatened (NT),
4 = Least Concern (LC)}
}
}
\source{
\url{https://www.iucnredlist.org/}
}
\usage{
training_labels
}
\description{
A dataset containing the International Union for the
Conservation of Nature's Global
Red List conservation assessments for 884 species of orchids.
}
\keyword{datasets}
