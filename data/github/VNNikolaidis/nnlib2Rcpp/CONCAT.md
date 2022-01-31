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
\name{Autoencoder}
\alias{Autoencoder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Autoencoder NN
}
\description{
A neural network for autoencoding data, projects data to a new set of variables.
}
\usage{
Autoencoder(
  data_in,
  desired_new_dimension,
  number_of_training_epochs,
  learning_rate,
  num_hidden_layers = 1L,
  hidden_layer_size = 5L,
  show_nn = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_in}{
data to be autoencoded, a numeric matrix, (2d, cases in rows, variables in columns).  It is recommended to be in [0 1] range.
}
  \item{desired_new_dimension}{
number of new variables to be produced (effectively the size of the special hidden layer that outputs the new variable values, thus the dimension of the output vector space).
}
  \item{number_of_training_epochs}{
number of training epochs, aka presentations of all training data to ANN during training.
}
  \item{learning_rate}{
the learning rate parameter of the Back-Propagation (BP) NN.
}
  \item{num_hidden_layers}{
number of hidden layers on each side of the special layer.
}
  \item{hidden_layer_size}{
number of nodes (Processing Elements or PEs) in each hidden layer
}
  \item{show_nn}{
boolean, option to display the (trained) ANN internal structure.
}
}
\value{
Returns a numeric matrix containing the projected data.
}
\references{
Nikolaidis V.N., Makris I.A, Stavroyiannis S, "ANS-based preprocessing of company performance indicators." Global Business and Economics Review 15.1 (2013): 49-58.
}
\author{
Vasilis N. Nikolaidis <vnnikolaidis@gmail.com>
}
\note{
This Autoencoder NN employs a {\code{\link{BP}}}-type NN to perform a data pre-processing step baring similarities to PCA since it too can be used for dimensionality reduction (Kramer 1991)(DeMers and Cottrell 1993)(Hinton and Salakhutdinov 2006). Unlike PCA, an autoencoding NN can also expand the feature-space dimensions (as feature expansion methods do). The NN maps input vectors to themselves via a special hidden layer (the coding layer, usually of different size than the input vector length) from which the new data vectors are produced. Note:
The internal BP PEs in computing layers apply the logistic sigmoid threshold function, and their output is in [0 1] range. It is recommended to use this range in your data as well. More for this particular autoencoder implementation can be found in (Nikolaidis, Makris, and Stavroyiannis 2013). The method is not deterministic and the mappings may be non-linear, depending on the NN topology.

(This function uses Rcpp to employ 'bpu_autoencoder_nn' class in nnlib2.)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{\code{\link{BP}}.}
\examples{
iris_s <- as.matrix(scale(iris[1:4]))
output_dim <- 2
epochs <- 100
learning_rate <- 0.73
num_hidden_layers <-2
hidden_layer_size <- 5

out_data <-  Autoencoder( iris_s, output_dim,
                          epochs, learning_rate,
                          num_hidden_layers, hidden_layer_size, FALSE)

plot( out_data,pch=21,
      bg=c("red","green3","blue")[unclass(iris$Species)],
      main="Randomly autoencoded Iris data")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ neural }% use one of  RShowDoc("KEYWORDS")
\name{nnlib2Rcpp-package}
\Rdversion{1.1}
\docType{package}
\alias{nnlib2Rcpp}

\title{A collection of Neural Networks and tools to create custom models}
\description{
This package contains a collection of ready-to-use Neural Networks (NN), i.e. versions of Autoencoder, Back-Propagation, Learning Vector Quantization and Matrix Associative Memory NN.
It also provides a module (NN module) to define and control custom neural networks created from user-defined nnlib2 NN components. More information and examples for each of the above can be found in its help documentation (see below).
}

\section{Ready-to-use Neural Networks:}{
\itemize{
\item{Plain Back-Propagation (BP-supervised) (\code{\link{BP}})}
\item{Learning Vector Quantization (LVQ-supervised) (\code{\link{LVQs}})}
\item{Learning Vector Quantization (LVQ-unsupervised) (\code{\link{LVQu}})}
\item{Matrix Associative Memory (MAM-supervised) (\code{\link{MAM}}})
\item{Autoencoder (unsupervised) (\code{\link{Autoencoder}})}
}}

\section{Custom Neural Networks:}{
\itemize{
\item{NN module (\code{\link{NN}})}
}}

\references{
\itemize{
\item{
Nikolaidis, V. N., (2021). The nnlib2 library and nnlib2Rcpp R package for implementing neural networks. Journal of Open Source Software, 6(61), 2876, \doi{10.21105/joss.02876}.}

References for the ready-to-use NN models (can be found in related help content):
\itemize{
\item{
Kohonen, T (1988). Self-Organization and Associative Memory, Springer-Verlag.; Simpson, P. K. (1991). Artificial neural systems: Foundations, paradigms, applications, and implementations. New York: Pergamon Press.}

\item{
Pao Y (1989). Adaptive Pattern Recognition and Neural Networks. Reading, MA (US); Addison-Wesley Publishing Co., Inc.}

\item{
Simpson, P. K. (1991). Artificial neural systems: Foundations, paradigms, applications, and implementations. New York: Pergamon Press.}

\item{
Philippidis, TP & Nikolaidis, VN & Kolaxis, JG. (1999). Unsupervised pattern recognition techniques for the prediction of composite failure. Journal of acoustic emission. 17. 69-81.}

\item{
Nikolaidis V.N., Makris I.A, Stavroyiannis S, "ANS-based preprocessing of company performance indicators." Global Business and Economics Review 15.1 (2013): 49-58, \doi{10.1504/GBER.2013.050667}.}
}
}
}
\seealso{
More information and examples on using the package can be found in the following vignette:

\code{vignette("manual", package='nnlib2Rcpp')}

Related links:
\itemize{
  \item \url{https://github.com/VNNikolaidis/nnlib2Rcpp}
  \item Package manual in PDF format at \url{https://github.com/VNNikolaidis/nnlib2Rcpp/blob/master/doc/manual.pdf})
  \item Report bugs, issues and suggestions at \url{https://github.com/VNNikolaidis/nnlib2Rcpp/issues}
}

}
\author{
Author/Maintainer:
\itemize{
  \item{Vasilis Nikolaidis \email{vnnikolaidis@gmail.com}}
}

Contributors:
\itemize{
  \item Arfon Smith [contributor]
  \item Dirk Eddelbuettel [contributor]
}
}
\name{NN-class}
\Rdversion{1.1}
\docType{class}
\alias{NN-class}
\alias{Rcpp_NN}
\alias{Rcpp_NN-class}
\alias{NN}
\alias{nn-class}
\alias{C++Object-class}
\alias{RcppClass-class}

\title{Class \code{"NN"}}
\description{
NN module, for defining and manipulating custom neural networks.
}
\section{Extends}{
Class \code{"\linkS4class{RcppClass}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\author{
Vasilis N. Nikolaidis <vnnikolaidis@gmail.com>
}
\note{
This R module maintains a generic neural network that can be manipulated using the provided methods. In addition to predefined, new neural network components can be defined and then employed by the \code{"NN"} module. Currently, definition of new components must be done in C++, requires the package source code (which includes the \pkg{nnlib2} C++ library of neural network base classes) and the ability to compile it.  In particular:
    \itemize{
  \item Any new component type or class definition can be added to a single header file called "\code{additional_parts.h}" (which is included in the package source). All new components to be employed by the \code{NN} module must be defined in this file (or be accessible from functions in this file).
  \item New \code{layer}, \code{connection_set}, \code{pe} or \code{connection} definitions must comply (at least loosely) to the \pkg{nnlib2} base class hierarchy and structure and follow the related guidelines. Note: some minimal examples of class and type definitions can be found in the "\code{additional_parts.h}" file itself.
  \item A textual name must be assigned to any new \code{layer} or \code{connection_set}, to be used as parameter in \code{NN} module methods that require a name to create a component. This can be as simple as a single line of code where given the textual name the corresponding component object is created and returned. This code must be added (as appropriate) to either \code{generate_custom_layer()} or \code{generate_custom_connection_set()} functions found in the same "\code{additional_parts.h}" header file. Note: example entries can be found in these functions at the "\code{additional_parts.h}" file.
}
More information on expanding the library with new types of NN components (nodes, layers, connections etc) and models, can be found in the package's  vignette as well as the related \href{https://github.com/VNNikolaidis/nnlib2Rcpp}{repository on Github}). Please consider submitting any useful components you create, to enrich  future versions of the package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{BP}}, \code{\link{LVQs}}, \code{\link{MAM}}.}

\examples{
# Example 1:

# (1.A) create new 'NN' object:

n <- new("NN")

# (1.B) Add topology components:

# 1. add a layer of 4 generic nodes:
n$add_layer("generic",4)
# 2. add a set for connections that pass data unmodified:
n$add_connection_set("pass-through")
# 3. add another layer of 2 generic nodes:
n$add_layer("generic",2)
# 4. add a set for connections that pass data x weight:
n$add_connection_set("wpass-through")
# 5. add a layer of 1 generic node:
n$add_layer("generic",1)
# Create actual full connections in sets, random initial weights in [0,1]:
n$create_connections_in_sets(0,1)
# Optionaly, show an outline of the topology:
n$outline()

# (1.C) use the network.

# input some data, and create output for it:
n$input_at(1,c(10,20,30,40))
n$recall_all(TRUE)
# the final output:
n$get_output_from(5)

# (1.D) optionally, examine the network:

# the input at first layer at position 1:
n$get_input_at(1)
# Data is passed unmodified through connections at position 2,
# and (by default) summed together at each node of layer at position 3.
# Final output from layer in position 3:
n$get_output_from(3)
# Data is then passed multiplied by the random weights through
# connections at position 4. The weights of these connections:
n$get_weights_at(4)
# Data is finally summed together at the node of layer at position 5,
# producing the final output, which (again) is:
n$get_output_from(5)

# Example 2: A simple MAM NN

# (2.A) Preparation:

# Create data pairs

iris_data    <- as.matrix( scale( iris[1:4] ) )
iris_species <- matrix(data=-1, nrow=nrow(iris_data), ncol=3)
for(r in 1:nrow(iris_data))
 iris_species[r ,as.integer( iris$Species )[r]]=1

# Create the NN and its components:

m <- new( "NN" )
m$add_layer( "generic" , 4 )
m$add_layer( "generic" , 3 )
m$fully_connect_layers_at(1, 2, "MAM", 0, 0)

# (2.B) Use the NN to store iris (data,species) pair:

# encode pairs in NN:

m$encode_datasets_supervised(
	iris_data,1,
	iris_species,3,0,
	1,TRUE)

# (2.C) Recall iris species from NN:

recalled_data <- m$recall_dataset(iris_data,1,3,TRUE)

# (2.D) Convert recalled data to ids and plot results:

recalled_ids <- apply(recalled_data, 1, which.max)
plot(iris_data, pch=recalled_ids)
}

\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{.CppObject}:}{Object of class \code{C++Object} ~~ }
    \item{\code{.CppClassDef}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.CppGenerator}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}
\section{Methods}{
  \describe{

    \item{\code{add_layer( name, size )}:}{Setup a new \code{layer} component (a layer of processing nodes) and append it to the NN topology. Parameters are:
    \itemize{
    \item{\code{name}}{: string, containing name (that also Specifies type) of new layer. Names of predefined layers currently include \code{'pe'}(same as \code{'generic'}), \code{'pass-through'}, \code{'which-max'}, \code{'MAM'}, \code{'LVQ-input'}, \code{'LVQ-output'}, \code{'BP-hidden'}, \code{'BP-output'}, \code{'perceptron'} (additional names for user-defined components may be used, see note below.)}
    \item{\code{size}}{: integer, layer size i.e. number of \code{pe} (Processing Elements or nodes) to create in the layer.}
    }
    }

    \item{\code{add_connection_set( name )}:}{Create a new empty \code{connection_set} component (a set of connections between two layers). It does not connect any layers nor contain any connections between specific layer nodes. The set is appended to the NN topology. Parameters are:
    \itemize{
    \item{\code{name}}{: string, containing name (that also specifies type) of new empty connection set. Names of predefined connection sets currently include \code{'generic', 'pass-through'}(which does not multiply weights), \code{'wpass-through'}(which does multiply weights), \code{'MAM'}, \code{'LVQ'}, \code{'BP'}, \code{'perceptron'} (additional names for user-defined components may be used, see note below).}
    }
    }

    \item{\code{create_connections_in_sets( min_random_weight, max_random_weight )}:}{Find empty, unconnected \code{connection_set} components that are between two  \code{layer}s in the topology, and set them up to connect the adjacent layers, adding connections to fully connect their nodes  (n x m connections are created, with n and m the number of nodes at each layer respectively). Parameters are:
    \itemize{
    \item{\code{min_random_weight}}{: double, minimum value for random initial connection weights.}
    \item{\code{max_random_weight}}{: double, maximum value for random initial connection weights.}
    }
    }

    \item{\code{connect_layers_at( source_pos, destin_pos, name )}:}{Insert a new empty \code{connection_set} component (a set of connections between two layers) between the layers at specified topology positions, and prepare it to connect them. No actual connections between any layer nodes are created. Parameters are:
    \itemize{
    \item{\code{source_pos}}{: integer, position in topology of source layer.}
    \item{\code{destin_pos}}{: integer, position in topology of destination layer.}
    \item{\code{name}}{: string, containing name (that also specifies type) of new connection set (see above).}
    }
    }

    \item{\code{fully_connect_layers_at( source_pos, destin_pos, name, min_random_weight, max_random_weight )}:}{Same as \code{connect_layers_at} but also fills the new \code{connection_set} with connections between the nodes of the two layers, fully connecting the layers (n x m connections are created, with n and m the number of nodes at each layer respectively). Parameters are:
    \itemize{
    \item{\code{source_pos}}{: integer, position in topology of source layer.}
    \item{\code{destin_pos}}{: integer, position in topology of destination layer.}
    \item{\code{name}}{: string, containing name (that also specifies type) of new connection set (see above).}
    \item{\code{min_random_weight}}{: double, minimum value for random initial connection weights.}
    \item{\code{max_random_weight}}{: double, maximum value for random initial connection weights.}
    }
    }

   \item{\code{add_single_connection( pos, source_pe, destin_pe, weight )}:}{
   Add a connection to a \code{connection_set} that already connects two layers. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of \code{connection_set} to which the new connection will be added.}
    \item{\code{source_pe}}{: integer, \code{pe} in source layer to connect.}
    \item{\code{destin_pe}}{: integer, \code{pe} in destination layer to connect.}
    \item{\code{weight}}{: double, value for initial connection weight.}
    }
    }

   \item{\code{remove_single_connection( pos, con )}:}{
   Remove a connection from a \code{connection_set}. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of \code{connection_set}.}
	\item{\code{con}}{: integer, connection to remove.}
    }
    }

   \item{\code{size()}:}{Returns neural network size, i.e. the number of components its topology.}

   \item{\code{sizes()}:}{Returns sizes of components in topology.}

  \item{\code{component_ids()}:}{Returns an integer vector containing the ids of the components in topology (these ids are created at run-time and identify each NN component).}

  \item{\code{input_at( pos, data_in )}:}{Input a data vector to the component (layer) at specified topology index. Returns TRUE if successful. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to receive input.}
    \item{\code{data_in}}{: NumericVector, data to be sent as input to component (sizes must match).}
    }
    }

 \item{\code{encode_at( pos )}:}{Trigger the encoding operation of the component at specified topology index (note: depending on implementation, an 'encode' operation usually collects inputs, processes the data, adjusts internal state variables and/or weights, and possibly produces output). Returns TRUE if successful. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to perform encoding.}
    }
    }

 \item{\code{recall_at( pos )}:}{Trigger the recall (mapping, data retrieval) operation of the component at specified topology index (note: depending on implementation, a 'recall' operation usually collects inputs, processes the data, and produces output). Returns TRUE if successful. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to perform recall.}
    }
    }

 \item{\code{encode_all( fwd )}:}{Trigger the encoding operation of all the components in the NN topology. Returns TRUE if successful. Parameters are:
    \itemize{
    \item{\code{fwd}}{: logical, set to TRUE to encode forwards (first-to-last component), FALSE to encode backwards (last-to-first component).}
    }
    }

 \item{\code{encode_dataset_unsupervised( data, pos, epochs, fwd )}:}{Encode a dataset using unsupervised training. A faster method to encode a data set. All the components in the NN topology will perform 'encode' in specified direction. Returns TRUE if successful. Parameters are:
    \itemize{
    \item{\code{data}}{: numeric matrix, containing input vectors as rows.}
    \item{\code{pos}}{: integer, position in topology of component to receive input vectors.}
    \item{\code{epochs}}{: integer, number of training epochs (encoding repetitions of the entire dataset).}
    \item{\code{fwd}}{: logical, indicates direction, TRUE to trigger encoding forwards (first-to-last component), FALSE to encode backwards (last-to-first component).}
    }
    }

 \item{\code{encode_datasets_supervised( i_data, i_pos, j_data, j_pos, j_destination_register, epochs, fwd )}:}{Encode multiple (i,j) vector pairs stored in two corresponding data sets, using supervised training. A faster method to encode the data. All the components in the NN topology will perform 'encode' in specified direction. Returns TRUE if successful. Parameters are:
    \itemize{
    \item{\code{i_data}}{: numeric matrix, data set, each row is a vector i of vector-pair (i,j).}
    \item{\code{i_pos}}{: integer, position in topology of component to receive i vectors.}
    \item{\code{j_data}}{: numeric matrix, data set, each row is a corresponding vector j of vector-pair (i,j).}
    \item{\code{j_pos}}{: integer, position in topology of component to receive j vectors.}
    \item{\code{j_destination_selector}}{: integer, selects which internal node (pe) registers will receive vector j, i.e. if 0 internal node register '\code{input}' will be used (j will become the layer's input), if 1 register '\code{output}' will be used (j will become the layer's output), if 2 register '\code{misc}' will be used (implementations may use this as an alternative way to transfer data to nodes without altering current input or output).}
    \item{\code{epochs}}{: integer, number of training epochs (encoding repetitions of the entire data).}
    \item{\code{fwd}}{: logical, indicates direction, TRUE to trigger encoding forwards (first-to-last component), FALSE to encode backwards (last-to-first component).}
    }
    }

 \item{\code{recall_dataset( data_in, input_pos, output_pos, fwd )}:}{Recall (map, retrieve output for) a dataset. A faster method to recall an entire data set. All the components in the NN topology will perform 'recall' in specified direction. Returns numeric matrix containing corresponding output. Parameters are:
    \itemize{
    \item{\code{data_in}}{: numeric matrix, containing input vectors as rows.}
    \item{\code{input_pos}}{: integer, position in topology of component to receive input vectors.}
    \item{\code{output_pos}}{: integer, position in topology of component to produce output.}
    \item{\code{fwd}}{: logical, indicates direction, TRUE to trigger 'recall' (mapping) forwards (first-to-last component), FALSE to recall backwards (last-to-first component).}
    }
    }

 \item{\code{recall_all( fwd )}:}{Trigger the recall (mapping, data retrieval) operation of all the components in the NN topology. Returns TRUE if successful. Parameters are:
    \itemize{
    \item{\code{fwd}}{: logical, set to TRUE to recall forwards (first-to-last component), FALSE to recall backwards (last-to-first component).}
    }
    }

 \item{\code{get_output_from( pos )}:}{Get the current output of the component at specified topology index. If successful, returns NumericVector of output values. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to use.}
    }
    }

\item{\code{get_output_at( pos )}:}{Same as \code{get_output_from}, see above.}

 \item{\code{get_input_at( pos )}:}{Get the current input of the component at specified topology index (depends on the implementation: for layers, this may be valid after the \code{pe}s (nodes) have performed their \code{input_function} on incoming values; \code{pe}s have an overridable \code{input_function} that collects all incoming values and produces a single value for further processing which is stored at an internal \code{input} register (whose value is retrieved here); by default \code{input_function} performs summation. If successful, returns NumericVector of final \code{input} values. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to use.}
    }
    }

 \item{\code{get_weights_at( pos )}:}{Get the current weights of the component (\code{connection_set}) at specified topology index. If successful, returns NumericVector of connection weights. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to use.}
    }
    }

 \item{\code{get_weight_at( pos, connection )}:}{Get the current weight of a connection in component (\code{connection_set}) at specified topology index. If successful, returns weight, otherwise 0. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to use.}
    \item{\code{connection}}{: connection to use.}
    }
    }

 \item{\code{set_weight_at( pos, connection, value )}:}{Set the weight of a connection in component (\code{connection_set}) at specified topology index. If successful, returns TRUE. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to use.}
    \item{\code{connection}}{: connection to use.}
    \item{\code{value}}{: new weight for connection.}
    }
    }

 \item{\code{set_misc_values_at( pos, data_in )}:}{Set the values in the \code{misc} data register that \code{pe} and \code{connection} objects maintain, for objects at specified topology index. If successful, returns TRUE. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to use.}
    \item{\code{data_in}}{: NumericVector, data to be used for new values in \code{misc} registers (sizes must match).}
    }
    }

 \item{\code{set_output_at( pos, data_in )}:}{Set the values in the \code{output} data register that \code{pe} objects maintain, for \code{layer} at specified topology index (currenly only \code{layer} components are supported). If successful, returns TRUE. Parameters are:
    \itemize{
    \item{\code{pos}}{: integer, position in topology of component to use.}
    \item{\code{data_in}}{: NumericVector, data to be used for new values in \code{misc} registers (sizes must match).}
    }
    }

\item{\code{print( )}:}{Print internal NN state, including all components in topology.}

\item{\code{outline( )}:}{Print a summary description of all components in topology.}

 }

The following methods are inherited (from the corresponding class):
objectPointer ("RcppClass"), initialize ("RcppClass"), show ("RcppClass")
}


\name{BP-class}
\Rdversion{1.1}
\docType{class}
\alias{BP-class}
\alias{BP}
\alias{Rcpp_BP}
\alias{Rcpp_BP-class}

\title{Class \code{"BP"}}
\description{
Supervised Back-Propagation (BP) NN module, for encoding input-output mappings.
}
\section{Extends}{
Class \code{"\linkS4class{RcppClass}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\references{
Simpson, P. K. (1991). Artificial neural systems: Foundations, paradigms, applications, and implementations. New York: Pergamon Press.
}
\author{
Vasilis N. Nikolaidis <vnnikolaidis@gmail.com>
}
\note{
This R module maintains an internal Back-Propagation (BP) multilayer perceptron NN (described in Simpson (1991) as the vanilla back-propagation algorithm), which can be used to store input-output vector pairs. Since the nodes (PEs) in computing layers of this BP implementation apply the logistic sigmoid threshold function, their output is in [0 1] range (and so should the desired output vector values).

(This object uses Rcpp to employ 'bp_nn' class in nnlib2.)}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Autoencoder}}.}

\examples{
# create some data...
iris_s                  <- as.matrix(scale(iris[1:4]))

# use a randomply picked subset of (scaled) iris data for training
training_cases          <- sample(1:nrow(iris_s), nrow(iris_s)/2,replace=FALSE)
train_set               <- iris_s[training_cases,]
train_class_ids         <- as.integer(iris$Species[training_cases])
train_num_cases         <- nrow(train_set)
train_num_variables     <- ncol(train_set)
train_num_classes       <- max(train_class_ids)

# create output dataset to be used for training, Here we encode class as 0s and 1s
train_set_data_out <- matrix(
          data = 0,
          nrow = train_num_cases,
          ncol = train_num_classes)

# now for each case, assign a 1 to the column corresponding to its class, 0 otherwise
# (there must be a better R way to do this)
for(r in 1:train_num_cases) train_set_data_out[r,train_class_ids[r]]=1

# done with data, let's use BP...
bp<-new("BP")

bp$encode(train_set,train_set_data_out,0.8,10000,2,4)

# let's test by recalling the original training set...
bp_output <- bp$recall(train_set)

cat("- Using this demo's encoding, recalled class is:\n")
print(apply(bp_output,1,which.max))
cat("- BP success in recalling correct class is: ",
  sum(apply(bp_output,1,which.max)==train_class_ids)," out of ",
  train_num_cases, "\n")

# Let's see how well it recalls the entire Iris set:
bp_output <- bp$recall(iris_s)

# show output
cat("\n- Recalling entire Iris set returns:\n")
print(bp_output)
cat("- Using this demo's encoding, original class is:\n")
print(as.integer(iris$Species))
cat("- Using this demo's encoding, recalled class is:\n")
bp_classification <- apply(bp_output,1,which.max)
print(bp_classification)
cat("- BP success in recalling correct class is: ",
  sum(apply(bp_output,1,which.max)==as.integer(iris$Species)),
  "out of ", nrow(iris_s), "\n")
plot(iris_s, pch=bp_classification, main="Iris classified by a partialy trained BP (module)")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{.CppObject}:}{Object of class \code{C++Object} ~~ }
    \item{\code{.CppClassDef}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.CppGenerator}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}
\section{Methods}{
  \describe{

    \item{\code{encode( data_in, data_out, learning_rate, training_epochs, hidden_layers, hidden_layer_size )}:}{ Setup a new BP NN and encode input-output data pairs. Parameters are:
    \itemize{
    \item{\code{data_in}}{: numeric matrix, containing input vectors as rows. . It is recommended that these values are in 0 to 1 range.}
    \item{\code{data_out}}{: numeric matrix, containing corresponding (desired) output vectors. It is recommended that these values are in 0 to 1 range.}
    \item{\code{learning_rate}}{: a number (preferably greater than 0 and less than 1) used in training.}
    \item{\code{training_epochs}}{: number of training epochs, aka single presentation iterations of all training data pairs to the NN during training.}
    \item{\code{hidden_layers}}{: number of hidden layers to be created between input and output layers.}
    \item{\code{hidden_layer_size}}{: number of nodes (processing elements or PEs) in the hidden layers (all hidden layers are of the same length in this model).}
    }
    Note: to encode additional input-output vector pairs in an existing BP, use \code{train_single} or \code{train_multiple} methods (see below).
    }

    \item{\code{recall(data_in)}:}{ Get output for a dataset (numeric matrix \code{data_in}) from the (trained) BP NN. }

    \item{\code{setup(input_dim, output_dim, learning_rate, hidden_layers, hidden_layer_size)}:}{ Setup the BP NN so it can be trained and used. Note: this is not needed if using \code{encode}. Parameters are:
    \itemize{
    \item{\code{input_dim}}{: integer length of input vectors.}
    \item{\code{output_dim}}{: integer length of output vectors.}
    \item{\code{learning_rate}}{: a number (preferably greater than 0 and less than 1) used in training.}
    \item{\code{hidden_layers}}{: number of hidden layers to be created between input and output layers.}
    \item{\code{hidden_layer_size}}{: number of nodes (processing elements or PEs) in the hidden layers (all hidden layers are of the same length in this model).}
    }
    }

    \item{\code{train_single (data_in, data_out)}:}{ Encode an input-output vector pair in the BP NN. Only performs a single training iteration (multiple may be required for proper encoding). Vector sizes should be compatible to the current NN (as resulted from the \code{encode} or \code{setup} methods). Returns error level indicator value.}

    \item{\code{train_multiple (data_in, data_out, training_epochs)}:}{ Encode multiple input-output vector pairs stored in corresponding datasets. Performs multiple iterations in epochs (see \code{encode}). Vector sizes should be compatible to the current NN (as resulted from the \code{encode} or \code{setup} methods). Returns error level indicator value.}

    \item{\code{print()}:}{ print NN structure. }

    \item{\code{load(filename)}:}{ retrieve the NN from specified file. }

    \item{\code{save(filename)}:}{ save the NN to specified file. }
  }

The following methods are inherited (from the corresponding class):
objectPointer ("RcppClass"), initialize ("RcppClass"), show ("RcppClass")
}
\name{MAM-class}
\Rdversion{1.1}
\docType{class}
\alias{MAM-class}
\alias{MAM}
\alias{Rcpp_MAM}
\alias{Rcpp_MAM-class}

\title{Class \code{"MAM"}}
\description{
A single Matrix Associative Memory (MAM) implemented as a (supervised) NN.
}
\section{Extends}{
Class \code{"\linkS4class{RcppClass}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
\references{
Pao Y (1989). Adaptive Pattern Recognition and Neural Networks. Reading, MA (US); Addison-Wesley Publishing Co., Inc.
}
\author{
Vasilis N. Nikolaidis <vnnikolaidis@gmail.com>
}
\note{
The NN in this module uses supervised training to store input-output vector pairs.

(This function uses Rcpp to employ 'mam_nn' class in nnlib2.)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{BP},\link{LVQs}}.}

\examples{
iris_s            <- as.matrix(scale(iris[1:4]))
class_ids         <- as.integer(iris$Species)
num_classes       <- max(class_ids)

# create output dataset to be used for training, Here we encode class as -1s and 1s
iris_data_out <- matrix( data = -1, nrow = nrow(iris_s), ncol = num_classes)

# now for each case, assign a 1 to the column corresponding to its class
for(r in 1:nrow(iris_data_out)) iris_data_out[r,class_ids[r]]=1

# Finally apply MAM:
# Encode train pairs in MAM and then get output dataset by recalling the test data.

mam <- new("MAM")

mam$encode(iris_s,iris_data_out)

# test the encoding by recalling the original input data...
mam_data_out <- mam$recall(iris_s)

# find which MAM output has the largest value and use this as the final cluster tag.
mam_recalled_cluster_ids = apply(mam_data_out,1,which.max)

plot(iris_s, pch=mam_recalled_cluster_ids, main="MAM recalled Iris data classes")

cat("MAM recalled these IDs:\n")
print(mam_recalled_cluster_ids)
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{.CppObject}:}{Object of class \code{C++Object} ~~ }
    \item{\code{.CppClassDef}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.CppGenerator}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}
\section{Methods}{
  \describe{

    \item{\code{encode( data_in, data_out )}:}{ Setup a new MAM NN and encode input-output data pairs. Parameters are:
    \itemize{
    \item{\code{data_in}}{: numeric matrix, input data to be encoded in MAM, a numeric matrix (2d, of n rows). Each row will be paired to the corresponding data_out row, forming an input-output vector pair.}
    \item{\code{data_out}}{: numeric matrix, output data to be encoded in MAM, a numeric matrix (2d, also of n rows). Each row will be paired to the corresponding data_in row, forming an input-output vector pair.}
    }
    Note: to encode additional input-output vector pairs in an existing MAM, use \code{train_single} method (see below).
    }

    \item{\code{recall(data)}:}{ Get output for a dataset (numeric matrix \code{data}) from the (trained) MAM NN. }

    \item{\code{train_single (data_in, data_out)}:}{ Encode an input-output vector pair in the MAM NN. Vector sizes should be compatible to the current NN (as resulted from the \code{encode} method).}

    \item{\code{print()}:}{ print NN structure. }

    \item{\code{load(filename)}:}{ retrieve the NN from specified file. }

    \item{\code{save(filename)}:}{ save the NN to specified file. }
  }

The following methods are inherited (from the corresponding class):
objectPointer ("RcppClass"), initialize ("RcppClass"), show ("RcppClass")
}
\name{LVQs-class}
\Rdversion{1.1}
\docType{class}
\alias{LVQs-class}
\alias{Rcpp_LVQs-class}
\alias{LVQs}

\title{Class \code{"LVQs"}}
\description{
Supervised Learning Vector Quantization (LVQ) NN module, for data classification.
}
\section{Extends}{
Class \code{"\linkS4class{RcppClass}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\references{
Simpson, P. K. (1991). Artificial neural systems: Foundations, paradigms, applications, and implementations. New York: Pergamon Press.
}
\author{
Vasilis N. Nikolaidis <vnnikolaidis@gmail.com>
}
\note{
The NN used in this module uses supervised training for data classification (described as Supervised Learning LVQ in Simpson (1991)). Data should be scaled in 0 to 1 range.

(This module uses Rcpp to employ 'lvq_nn' class in nnlib2.)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{LVQu}} (unsupervised LVQ function).
}

\examples{
# LVQ expects data in 0 to 1 range, so scale some numeric data...
iris_s<-as.matrix(iris[1:4])
c_min<-apply(iris_s,2,FUN = "min")
c_max<-apply(iris_s,2,FUN = "max")
c_rng<-c_max-c_min
iris_s<-sweep(iris_s,2,FUN="-",c_min)
iris_s<-sweep(iris_s,2,FUN="/",c_rng)

# create a vector of desired class ids (starting from 0):
iris_desired_class_ids<-as.integer(iris$Species)-1;

# now create the NN:
lvq<-new("LVQs")

# and train it:
lvq$encode(iris_s,iris_desired_class_ids,100)

# recall the same data (a simple check of how well the LVQ was trained):
lvq_recalled_class_ids<-lvq$recall(iris_s);
plot(iris_s, pch=lvq_recalled_class_ids, main="LVQ recalled clusters (module)")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{.CppObject}:}{Object of class \code{C++Object} ~~ }
    \item{\code{.CppClassDef}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{.CppGenerator}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}
\section{Methods}{
  \describe{

    \item{\code{encode(data, desired_class_ids, training_epochs)}:}{ Encode input and output (classification) for a dataset using LVQ NN. Parameters are:}
    \itemize{
    \item{\code{data}}{: training data, a numeric matrix, (2d, cases in rows, variables in columns). Data should be in 0 to 1 range.}
    \item{\code{desired_class_ids}}{ : vector of integers containing a desired class id for each training data case (row). Should contain integers in 0 to n-1 range, where n is the number of classes.}
  \item{\code{training_epochs}}{: integer, number of training epochs, aka presentations of all training data to the NN during training.}
  }

    \item{\code{recall(data_in)}:}{ Get output (classification) for a dataset (numeric matrix \code{data_in}) from the (trained) LVQ NN. The \code{data_in} dataset should be 2-d containing  data cases (rows) to be presented to the NN and is expected to have same number or columns as the original training data. Returns a vector of integers containing a class id for each case (row).}

    \item{\code{print()}:}{ print NN structure. }

    \item{\code{load(filename)}:}{ retrieve the NN from specified file. }

    \item{\code{save(filename)}:}{ save the NN to specified file. }
  }

The following methods are inherited (from the corresponding class):
objectPointer ("RcppClass"), initialize ("RcppClass"), show ("RcppClass")
}
\name{LVQu}
\alias{LVQu}
\alias{SOM}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unsupervised LVQ
}
\description{
Unsupervised (clustering) Learning Vector Quantization (LVQ) NN.
}
\usage{
LVQu(
  data,
  max_number_of_desired_clusters,
  number_of_training_epochs,
  neighborhood_size,
  show_nn )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data to be clustered, a numeric matrix, (2d, cases in rows, variables in columns). Data should be in 0 to 1 range.
}
  \item{max_number_of_desired_clusters}{
clusters to be produced (at most)
}
  \item{number_of_training_epochs}{
number of training epochs, aka presentations of all training data to ANN during training.
}
  \item{neighborhood_size}{
integer >=1, specifies affected neighbor output nodes during training. if 1 (Single Winner) the ANN is somewhat similar to k-means.
}
  \item{show_nn}{
boolean, option to display the (trained) ANN internal structure.
}
}
\value{
Returns a vector of integers containing a cluster id for each data case (row).
}
\references{
Kohonen, T (1988). Self-Organization and Associative Memory, Springer-Verlag.; Simpson, P. K. (1991). Artificial neural systems: Foundations, paradigms, applications, and implementations. New York: Pergamon Press.

Philippidis, TP & Nikolaidis, VN & Kolaxis, JG. (1999). Unsupervised pattern recognition techniques for the prediction of composite failure. Journal of acoustic emission. 17. 69-81.
}
\author{
Vasilis N. Nikolaidis <vnnikolaidis@gmail.com>
}
\note{
Function LVQu employs an unsupervised LVQ for clustering data (Kohonen 1988). This LVQ variant is described as Unsupervised Learning LVQ in Simpson (1991) and is a simplified 1-D version of Self-Organizing-Map (SOM). Its parameter \code{neighborhood_size} controls the encoding mode (where \code{neighborhood_size}=1 is Single-Winner Unsupervised encoding, similar to k-means, while an odd valued \code{neighborhood_size} > 1 means Multiple-Winner Unsupervised encoding mode). Initial weights are random (uniform distribution) in 0 to 1 range. As these weights represent cluster center coordinates (the class reference vector), it is important that input data is also scaled to this range.

(This function uses Rcpp to employ 'som_nn' class in nnlib2.)
}
\seealso{
  \code{\link{LVQs}} (supervised LVQ module),
}
\examples{
# LVQ expects data in 0 to 1 range, so scale...
iris_s<-as.matrix(iris[1:4])
c_min<-apply(iris_s,2,FUN = "min")
c_max<-apply(iris_s,2,FUN = "max")
c_rng<-c_max-c_min
iris_s<-sweep(iris_s,2,FUN="-",c_min)
iris_s<-sweep(iris_s,2,FUN="/",c_rng)

cluster_ids<-LVQu(iris_s,5,100)
plot(iris_s, pch=cluster_ids, main="LVQ-clustered Iris data")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ neural }% __ONLY ONE__ keyword per line
