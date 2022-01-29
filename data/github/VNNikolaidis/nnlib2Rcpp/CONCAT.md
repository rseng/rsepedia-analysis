# nnlib2Rcpp
An R package providing access to Neural Networks created using nnlib2. 

---

A collection of Artificial Neural Networks (NNs or ANNs or ANSs) created using the 'nnlib2' C++ library.

Currently includes predefined versions of BP, Autoencoder, MAM, LVQ (supervised and unsupervised). 
All NNs are created using 'nnlib2' (a C++ library of classes for implementing NNs) and interfaced with R via RCpp.

The package also provides the NN R module (Class "NN") which allows creation and control of custom NNs configurations from R using NN components (predefined or user-defined) created using 'nnlib2'. To add new user-defined NN components (layers, nodes, connections, sets of connections etc) see the "NN" component documentation (type ?NN in R). Note: this process requires some familiarity with C++.

---

To install:

(a) From CRAN Repository (recommended): The CRAN (stable) version of this package can be installed the usual way, i.e. by invoking the following R command:

    install.packages("nnlib2Rcpp") 

(b) From GitHub: To add the GitHub (latest) version of this package to your R installation, use the following R commands:

    library(devtools) 
    install_github("VNNikolaidis/nnlib2Rcpp")

(c) From r-universe: To add the package (corresponding to the latest GitHub release version) to your R installation, use the following R command:

    install.packages('nnlib2Rcpp', repos = 'https://vnnikolaidis.r-universe.dev')

Once installed, for package help (including documentation and examples for each function or class provided by nnlib2Rcpp) use the following R command:

    help(package='nnlib2Rcpp')

while the package vingette (containing information on adding custom components) can be viewed using the following R command:

    vignette("manual", package='nnlib2Rcpp')

The package vingette is also available in PDF format here:

https://github.com/VNNikolaidis/nnlib2Rcpp/blob/master/doc/manual.pdf

A reference manual in PDF format (for the last version in CRAN) can be found here:

https://cran.r-project.org/web/packages/nnlib2Rcpp/nnlib2Rcpp.pdf

---

For information on citing this package use the following R command for information:

    citation("nnlib2Rcpp")

---

For copyright information see LICENSE.md file or DESCRIPTION+LICENSE files (as imposed by package format for CRAN submissions).

---

The ‘nnlib2’ library used (and included) in this package is a collection of C++ base classes and templates for creating NNs. This library is also available as a standalone project, in GitHub repository (https://github.com/VNNikolaidis/nnlib2). For a (simplified) class-diagram of significant nnlib2 classes and templates see: https://github.com/VNNikolaidis/nnlib2/blob/master/misc/diagram%20of%20main%20classes.png

For implementing new NN components and models in nnlib2 that can be used in nnlib2Rcpp, see also: 

https://r-posts.com/creating-custom-neural-networks-with-nnlib2rcpp/ ( permalink: https://wp.me/p8rgs6-sh )

As mentioned earlier, more instructions on using 'nnlib2' and 'nnlib2Rcpp' can be found in the package vingette, also available in PDF format here:

https://github.com/VNNikolaidis/nnlib2Rcpp/blob/master/doc/manual.pdf

Link to related paper in the Journal of Open Source Software:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02876/status.svg)](https://doi.org/10.21105/joss.02876)

---

My future goals for this project (iF AND WHEN time permits):

- to add more classic neural network component and model implementations in nnlib2 (and thus nnlib2Rcpp).
- to parallelize component base classes.

Let me know if interested to contribute, or want to add your neural network components to the package. Or, as stated below:

---

We invite anyone to contribute to this software and/or provide feedback, suggestions, report issues or problems. Possible improvements and contributions may include (but are not limited to) implementation of additional neural network models using nnlib2 classes and templates (and thus new neural network components compatible with "NN" module in nnlib2Rcpp), parallelism (possibly via OpenMP), replacement of custom data structures with STL containers, performance enhancements etc.

Please use the issues option in GitHub or email (vnnikolaidis AT gmail.com) if interested to contribute.

---

[![](https://cranlogs.r-pkg.org/badges/nnlib2Rcpp)](https://cran.r-project.org/package=nnlib2Rcpp)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4780957.svg)](https://doi.org/10.5281/zenodo.4780957)
Changes to nnlib2Rcpp version 0.1.1 (from 0.1.0)
- removed redundant unused files.
- minor documentation fixes (Rcpp missing entries added), passes devtools::check() with no related warnings. 
- since LVQ has two variations (supervised and unsupervised), it was mentioned in comments that the names may be misleading or confusing. So the LVQ related R functions and modules are renamed as follows:
  ‘LVQs’ for supervised LVQ function (was ‘LVQ’ in v.0.1.0)
  ‘LVQs_NN’ for supervised LVQ module (was ‘LVQ_NN’ in v.0.1.0)
  ‘LVQu’ for unsupervised LVQ function (was ‘SOM’ in v.0.1.0)
- minor documentation improvements.
- pdf documentation removed, but can be built by running the following R command: devtools::build_manual()

---

Changes to nnlib2Rcpp version 0.1.2 (from 0.1.1)
-	(checked with R 3.6.3 and RTools35.)
-	(checked with R 4.0.0 and RTools40 (produces 1 extra NOTE).
-	removed redundant BP functions. All BP functionality is now implemented in BP module class (safer and more convenient for supervised methods).
-	removed redundant Supervised-LVQ functions. All Supervised-LVQ functionality is now implemented in LVQs module class (safer and more convenient for supervised methods). Unsupervised-LVQ is still implemented as a function (LVQu).
-	removed MAM functions. All MAM functionality is now implemented in MAM module class (safer and more convenient for supervised methods).
-	Minor changes in nnlib2 (now uses nnlib2 v.0.1.5)
-	Fixed LICENCE file (was incompatible with CRAN)

---

Changes to nnlib2Rcpp version 0.1.3 (from 0.1.2)
-	includes nnlib2 v.0.1.7.

---

Changes to nnlib2Rcpp version 0.1.4 (from 0.1.3)
-	(checked on Linux Mint 19.3 with R 3.6.3)
-	(checked with R 4.0.0 and RTools40 (produces 1 extra NOTE).
-	includes nnlib2 v.0.2.0.
-	added module "NN" which allows creation and control of custom NNs from R using predefined components. It also provides fixed method for adding user-defined NN components which can be used in the module.

---

Changes to nnlib2Rcpp version 0.1.5 (from 0.1.4)
-	added some methods to module "NN" (get_weight_at, set_weight_at)
-	added some methods to module "NN" (set_misc_values_at)

---

Changes to nnlib2Rcpp version 0.1.6 (from 0.1.5)
-	added method to module "NN" (recall_dataset) to decode/map/retrieve output for entire input dataset.
-	added method to module "NN" (encode_dataset_unsupervised) for faster unsupervised training.
-	added method to module "NN" (encode_datasets_supervised) for faster supervised training.
-	includes nnlib2 v.0.2.4
-	added method to module "NN" (get_output_at).
-	added method to module "NN" (set_output_at).
-	added package vingette
-	other minor changes

---

Changes to nnlib2Rcpp version 0.1.7 (from 0.1.6)
-	minor changes (only affect source).

---

Changes to nnlib2Rcpp version 0.1.8 (from 0.1.7)
-	minor documentation changes.
-	other minor changes
---
title: 'The nnlib2 library and nnlib2Rcpp R package for implementing neural networks'
tags:
  - neural networks
  - R
  - Cpp
authors:
 - name: Vasilis N Nikolaidis
   orcid: 0000-0003-1471-8788
   affiliation: 1
affiliations:
 - name: University of Peloponnese
   index: 1
date: 15 September 2020
bibliography: paper.bib
---

# Summary

Artificial Neural Networks (ANN or NN) are computing models used in various data-driven applications. Such systems typically consist of a large number of processing elements (or nodes), usually organized in layers, which exchange data via weighted connections. An ever-increasing number of different neural network models have been proposed and used. Among the several factors differentiating each model are the network topology, the processing and training methods in nodes and connections, and the sequences utilized for transferring data to, within and from the model etc. The software presented here is a C++ library of classes and templates for implementing neural network components and models and an R package that allows users to instantiate and use such components from the R programming language.

# Statement of need

A significant number of capable, flexible, high performance tools for NN are available today, including frameworks such as `Tensorflow` [@abadi2016tensorflow] and `Torch` [@torch], and related high level APIs including `Keras` [@keras1] and `PyTorch` [@PyTorch1]. Ready-to-use NN models are also provided by various machine learning platforms such as `H2O` [@h2o_platform] or libraries, `SNNS` [@Zell1994] and `FANN` [@FANN]. Ready to use NNs are available in a large number of software packages (`WEKA`, `SPSS`, `Matlab`, `NeuroSolutions`, `Noesis` and many others). The R language [@RLanguage], widely used for data analysis, includes package `nnet` [@CRANnnet] in its standard packages, while a significant number of NN-related extension packages can be found in CRAN and other repositories. These include interfaces and bindings to aforementioned tools (`keras` [@chollet2017kerasR], `torch` [@CRANtorch], `rTorch` [@CRANrTorch], `h2o ` [@CRANh2o], `RSNNS` [@CRANRSNNS]), as well as several other NN-specific packages that provide ready-to-use NN models (`deepnet` [@CRANdeepnet], `deepNN` [@CRANdeepNN], RcppDL [@CRANRcppDL] etc.).

The `nnlib2` library adds another tool for NN experimentation and implementation. It is a collection of C++  base classes and class-templates for defining NN parts, components and entire models. It contains base classes (and related functionality) for each part of a NN (processing nodes, connections, layers, groups of connections etc.) and for NN models that contain such parts. The library does not focus (in its current state) on high-performance computing but on being a versatile basis on which new NN parts and models can be defined with consistent, readable, maintainable and reusable code. Although the software does contain some ready-to-use NN models and their related reusable NN parts (currently for versions of Back-propagation, autoencoder, Learning Vector Quantization and Matrix Associative Memory), its main purpose is to aid the creation of new custom NN parts and the experimentation with them in arbitrary NN models and configurations. Overall the library (and related R package):

1. Allows the implementation of NN parts and models from reusable components that follow a common interface thus can easily be modified, combined with each other etc, making them suitable for educational or experimentation purposes.
2. Provides a basis for implementation of different NN parts and models, including classic NNs.
3. Is versatile enough for experimentation with new, custom components, configurations and models.
4. Has no dependencies with other software; to build or employ the `nnlib2` parts and models only a standard C++ compiler is required, and the produced models are standalone, lightweight, highly portable and can be embedded in any C++ application.
5. Allows for the exact same NN parts (processing nodes, connections, layers, groups of connections etc.) defined using `nnlib2` to be employed, combined, manipulated and monitored in R, via package `nnlib2Rcpp` which includes the entire `nnlib2` library and also provides supporting R functions and modules for this purpose. Thus, new NN components can be developed entirely using R-related tools (such as Rtools and RStudio), and then be used in R or transferred to a pure C++ application if needed.

# Discussion

As implied earlier, the `nnlib2` library may appeal to NN students and experimenters who seek a basis for implementing various NN models in C++ and take advantage of the versatility and direct control of the model that this approach allows. Furthermore, the library lets new NN definitions be written in code that is comparatively simple and self-explanatory; most of the provided base classes and templates correspond to typical NN components (processing nodes, layers of nodes, connections, sets of connections between layers, entire neural networks etc.) and closely follow the conceptual model often found in related bibliography such as [@SimpsonANS;@McCord91;@pao1989adaptive;@TheoPR] and elsewhere. Typical component functionality and features are provided in these base definitions, while diverse, model-specific components are to be implemented in custom sub-classes. By following a common interface, components (or even entire NN models) can be reused in other models. Finally, the library presents a high degree of target and compiler independence, with early versions used in various solutions targeting different operating systems, mainly for research purposes, for example in [@philippidis1999unsupervised;@Nikolaidis199Phd;@nikolaidis2013ans]. 

To take advantage of the R ecosystem, improve its usability and widen its audience, the `nnlib2` library has also been integrated into an R package (a process supported by `Rcpp` [@Rcpp1]). The package, named `nnlib2Rcpp` (available on CRAN [@nnlib2RcppCRAN] and GitHub [@nnlib2RcppGitHub]) includes the entire `nnlib2` code as well as an additional R class (class NN) with methods that allow users to instantiate the `nnlib2` NN components they define, combine them in custom NN configurations and generally control, monitor and employ them from R. Finally, the package provides a small collection of predefined NN components and complete ready-to-use NN models (for versions of Back Propagation, autoencoder, supervised/unsupervised Learning Vector Quantization and Matrix Associative Memory neural networks).

# References
