[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5748062.svg)](https://doi.org/10.5281/zenodo.4382739)
[![codecov](https://codecov.io/gh/jameschapman19/cca_zoo/branch/master/graph/badge.svg?token=JHG9VUB0L8)](https://codecov.io/gh/jameschapman19/cca_zoo)
![Build Status](https://github.com/jameschapman19/cca_zoo/actions/workflows/python-package.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/cca-zoo/badge/?version=latest)](https://cca-zoo.readthedocs.io/en/latest/?badge=latest)
[![version](https://img.shields.io/pypi/v/cca-zoo)](https://pypi.org/project/cca-zoo/)
[![downloads](https://img.shields.io/pypi/dm/cca-zoo)](https://pypi.org/project/cca-zoo/)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03823/status.svg)](https://doi.org/10.21105/joss.03823)

# CCA-Zoo

`cca-zoo` is a collection of linear, kernel, and deep methods for canonical correlation analysis of multiview data. 
Where possible it follows the `scikit-learn`/`mvlearn` APIs and models therefore have `fit`/`transform`/`fit_transform` methods as standard.

## Installation

Dependency of some implemented algorithms are heavy, such as `pytorch` and `numpyro`. 
We provide several options to accomodate the user's needs.
For full details of algorithms included, please refer to section [Implemented Methods](#implemented-methods)

Standard installation: 

```
pip install cca-zoo
```
For deep learning elements use:
```
pip install cca-zoo[deep]
```

For probabilistic elements use:
```
pip install cca-zoo[probabilistic]
```
## Documentation
Available at https://cca-zoo.readthedocs.io/en/latest/
  
## Citation:
CCA-Zoo is intended as research software. Citations and use of our software help us justify the effort which has gone into, and will keep going into, maintaining and growing this project. Stars on the repo are also greatly appreciated :)

If you have used CCA-Zoo in your research, please consider citing our JOSS paper:

Chapman et al., (2021). CCA-Zoo: A collection of Regularized, Deep Learning based, Kernel, and Probabilistic CCA methods in a scikit-learn style framework. Journal of Open Source Software, 6(68), 3823, https://doi.org/10.21105/joss.03823

With bibtex entry:

```bibtex
@article{Chapman2021,
  doi = {10.21105/joss.03823},
  url = {https://doi.org/10.21105/joss.03823},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3823},
  author = {James Chapman and Hao-Ting Wang},
  title = {CCA-Zoo: A collection of Regularized, Deep Learning based, Kernel, and Probabilistic CCA methods in a scikit-learn style framework},
  journal = {Journal of Open Source Software}
}
```

## Implemented Methods

### Standard Install
- CCA (Canonical Correlation Analysis): Solutions based on either alternating least squares or as the solution to genrralized eigenvalue problem
- PLS (Partial Least Squares)
- [rCCA (Ridge Regularized Canonical Correlation Analysis)](https://www.sciencedirect.com/science/article/abs/pii/0304407676900105?via%3Dihub)
- [GCCA (Generalized CCA)](https://academic.oup.com/biomet/article-abstract/58/3/433/233349?redirectedFrom=fulltext)
- MCCA (Multiset CCA)
- K(M)CCA (kernel Multiset CCA)
- [TCCA (Tensor CCA)](https://arxiv.org/pdf/1502.02330.pdf)
- [KTCCA (kernel Tensor CCA)](https://arxiv.org/pdf/1502.02330.pdf)
- [SCCA (Sparse CCA)](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13043)
- [SPLS (Sparse PLS/Penalized Matrix Decomposition](https://web.stanford.edu/~hastie/Papers/PMD_Witten.pdf)
- [ElasticCCA (Penalized CCA)](https://pubmed.ncbi.nlm.nih.gov/19689958/)
- [SWCCA (Sparse Weighted CCA)](https://arxiv.org/abs/1710.04792v1#:~:text=However%2C%20classical%20and%20sparse%20CCA%20models%20consider%20the,where%20weights%20are%20used%20for%20regularizing%20different%20samples)
- [SpanCCA](http://akyrillidis.github.io/pubs/Conferences/cca.pdf)

### `[deep]` Install
- DCCA (Deep CCA)

  Using either Andrew's original [Tracenorm Objective](https://ttic.uchicago.edu/~klivescu/papers/andrew_icml2013.pdf) or Wang's [alternating least squares solution](https://arxiv.org/pdf/1510.02054v1.pdf)
  
- [DGCCA (Deep Generalized CCA)](https://www.aclweb.org/anthology/W19-4301.pdf)

  An alternative objective based on the linear GCCA solution. Can be extended to more than 2 views
 
- [DMCCA (Deep Multiset CCA)](https://arxiv.org/abs/1904.01775)

  An alternative objective based on the linear MCCA solution. Can be extended to more than 2 views
  
- [DTCCA (Deep Tensor CCA)](https://arxiv.org/pdf/2005.11914.pdf)
- [DCCAE (Deep Canonically Correlated Autoencoders)](http://proceedings.mlr.press/v37/wangb15.pdf)
- [DVCCA/DVCCA Private (Deep variational CCA)](https://arxiv.org/pdf/1610.03454.pdf)

### `[probabilistic]` Install
- [Variational Bayes CCA](https://ieeexplore.ieee.org/document/4182407)

## Contributions
A guide to contributions is available at https://cca-zoo.readthedocs.io/en/latest/developer_info/contribute.html

## Sources

I've added this section to give due credit to the repositories that helped me in addition to their copyright notices in
the code where relevant.

### Other Implementations of (regularised)CCA/PLS

[MATLAB implementation](https://github.com/anaston/PLS_CCA_framework)

### Implementation of Sparse PLS

MATLAB implementation of SPLS by [@jmmonteiro](https://github.com/jmmonteiro/spls)

### Other Implementations of DCCA/DCCAE

Keras implementation of DCCA from [@VahidooX's github page](https://github.com/VahidooX)

The following are the other implementations of DCCA in MATLAB and C++. These codes are written by the authors of the original paper:

[Torch implementation](https://github.com/Michaelvll/DeepCCA) of DCCA from @MichaelVll & @Arminarj

C++ implementation of DCCA from Galen Andrew's [website](https://homes.cs.washington.edu/~galen/)

MATLAB implementation of DCCA/DCCAE from Weiran Wang's [website](http://ttic.uchicago.edu/~wwang5/dccae.html)

MATLAB implementation of [TCCA](https://github.com/rciszek/mdr_tcca)

### Implementation of VAE

[Torch implementation of VAE](https://github.com/pytorch/examples/tree/master/vae)
---
title: 'CCA-Zoo: A collection of Regularized, Deep Learning based, Kernel, and Probabilistic CCA methods in a scikit-learn style framework'
tags:
  - Python
  - Multiview
  - Machine Learning 
authors:
  - name: James Chapman 
    orcid: 0000-0002-9364-8118 
    affiliation: 1
  - name: Hao-Ting Wang 
    orcid: 0000-0003-4078-2038 
    affiliation: "2, 3" 
affiliations:
  - name: Centre for Medical Image Computing, University College London, London, UK 
    index: 1
  - name: Centre de Recherche de l'Institut Universitaire de Gériatrie de Montréal, Université de Montréal, Montréal, QC, Canada
    index: 2 
  - name: Centre de Recherche de l'Hôpital du Sacré Coeur de Montréal, Université de Montréal, Montréal, QC, Canada
    index: 3 
date: 1 September 2021 
bibliography: paper.bib
---

# Summary

Multi-view data has gained visibility in scientific research. Examples include different languages in natural language 
processing, as well as neuroimaging, multiomics and audiovisual data. Canonical Correlation Analysis (CCA)
[@hotelling1992relations]  and Partial Least Squares (PLS) are classical methods for investigating and quantifying
multivariate relationships between these views of data. The goal of CCA and its variants is to find projections (and
associated weights) for each view of the data into a latent space where they are highly correlated.

The original CCA is constrained by the sample-to-feature ratio. The algorithm cannot produce a solution when the number 
of features in one view exceeds the number of samples. To overcome this restriction, the original CCA has been developed 
into a family of models which include regularised [@vinod1976canonical], kernelized [@hardoon2004canonical], 
probabilistic/generative [@bach2005probabilistic], and deep learning based [@andrew2013deep] variants. In particular 
these variations have allowed practitioners to apply these models to complex, high dimensional data. Similarly, 
variants of PLS have been proposed including the widely used Penalized Matrix Decomposition algorithm 
[@witten2009penalized] which induces sparsity in the weight vectors for interpretability and generalisation.

`cca-zoo` is a Python package that implements many variants in a simple API with standardised outputs. We would like to 
highlight the unique benefits our package brings to the community in comparison to other established Python packages 
containing implementations of CCA. Firstly, `cca-zoo` contains a number of regularised CCA and PLS for high dimensional 
data that have previously only been available in installable packages in R. Native Python implementation will give 
Python users convenient access to these powerful models for both application and the development of new algorithms. 
Secondly,`cca-zoo` contains several deep CCA variants written in PyTorch [@paszke2019pytorch]. We adopted a modular 
style allowing users to apply their desired neural network architectures for each view for their own training pipeline. 
Thirdly, `cca-zoo` contains generative models including probabilistic and deep variational CCA. This class of variations 
can be used to model the multiview data generation process and even generate new synthetic samples. Finally, `cca-zoo` 
provides data simulation utilities to synthesize data containing specified correlation structures as well as the paired 
MNIST data commonly used as a toy dataset in deep multiview learning.

# Statement of need

The Python ecosystem for multiview learning currently provides a few options for implementing CCA and PLS models.
`scikit-learn` [@pedregosa2011scikit] contains standard implementations of both CCA and PLS for two-view data which plug
into their mature API. `pyrcca` [@bilenko2016pyrcca] contains implementations of ridge regularised and kernelized
two-view CCA. The embed module of `mvlearn` [@perry2020mvlearn] is perhaps the closest relative of `cca-zoo`, containing
implementations of ridge regularised and kernelized multi-view CCA. `cca-zoo` builds on the `mvlearn` API by providing
an additional range of regularised models and in particular sparsity inducing models which have found success in
multiomics. Building on the reference implementation in `mvlearn`, `cca-zoo` further provides a number of deep learning
models with a modular design to enable users to supply their own choice of neural network architectures.

Standard implementations of state-of-the-art models help as benchmarks for methods development and easy application to
new datasets. `cca-zoo` extends the existing ecosystem with a number of sparse regularised CCA models. These variations
have found popularity in genetics and neuroimaging where signals are contained in a small subset of variables. With
applications like these in mind, `cca-zoo` simplified the access to the learnt model weights to perform further analysis
in the feature space. Furthermore, the modular implementations of deep CCA and its multiview variants allow the user to
focus on architecture tuning. Finally, `cca-zoo` adds generative models including variational [@wang2007variational] and
deep variational CCA [@wang2016deep] as well as higher order canonical correlation analysis with tensor [@kim2007tensor]
and deep tensor CCA [@wong2021deep].

# Implementation

`cca-zoo` adopted a similar API to that used in `scikit-learn`. The user first instantiates a model object and its
relevant hyperparameters. Next they call the model's `fit()` method to apply the data. After fitting, the model object
contains its relevant parameters such as weights or dual coefficients (for kernel methods) which can be accessed for
further analysis. For models that fit with iterative algorithms, the model may also contain information about the
convergence of the objective function. After the model has been fit, its `transform()` method can project views into
latent variables and `score()` can be used to measure the canonical correlations.

The deep and probabilistic models are supported by PyTorch and NumPyro respectively. Due to the size of these
dependencies, these two classes of variations are not in the default installation. Instead, we provide options [deep]
and [probabilistic] for users. The list bellow provides the complete collection of models along with their installation
tag is provided below.

## Model List

A complete model list at the time of publication:

| Model Class | Model Name | Number of Views | Install |
| -------- | -------- | ------ |-----|
| CCA   | Canonical Correlation Analysis | 2   | standard |
| rCCA   | Canonical Ridge | 2   | standard |
| KCCA   | Kernel Canonical Correlation Analysis | 2   | standard |
| MCCA   | Multiset Canonical Correlation Analysis | \>=2   | standard |
| KMCCA   | Kernel Multiset Canonical Correlation Analysis | \>=2   | standard |
| GCCA   | Generalized Canonical Correlation Analysis | \>=2   | standard |
| KGCCA   | Kernel Generalized Canonical Correlation Analysis | \>=2   | standard |
| PLS   | Partial Least Squares | \>=2   | standard |
| CCA_ALS   | Canonical Correlation Analysis by Alternating Least Squares) [@golub1995canonical] | \>=2   | standard |
| PLS_ALS   | Partial Least Squares by Alternating Least Squares)  | \>=2   | standard |
| PMD   | Sparse CCA by Penalized Matrix Decomposition | \>=2   | standard |
| ElasticCCA   | Sparse Penalized CCA [@waaijenborg2008quantifying] | \>=2   | standard |
| ParkhomenkoCCA   | Sparse CCA [@parkhomenko2009sparse] | \>=2   | standard |
| SCCA   | Sparse Canonical Correlation Analysis by Iterative Least Squares [@mai2019iterative] | \>=2   | standard |
| SCCA_ADMM   | Sparse Canonical Correlation Analysis by Altnerating Direction Method of Multipliers [@suo2017sparse] | \>=2   | standard |
| SpanCCA   | Sparse Diagonal Canonical Correlation Analysis [@asteris2016simple] | \>=2   | standard |
| SWCCA   | Sparse Weighted Canonical Correlation Analysis [@wenwen2018sparse] | \>=2   | standard |
| TCCA   | Tensor Canonical Correlation Analysis | \>=2   | standard |
| KTCCA   | Kernel Tensor Canonical Correlation Analysis [@kim2007tensor] | \>=2   | standard |
| DCCA   | Deep Canonical Correlation Analysis | \>=2   | deep |
| DCCA_NOI   | Deep Canonical Correlation Analysis by Non-Linear Orthogonal Iterations [@wang2015stochastic] | \>=2   | deep |
| DCCAE   | Deep Canonically Correlated Autoencoders [@wang2015deep] | \>=2   | deep |
| DTCCA   | Deep Tensor Canonical Correlation Analysis | \>=2   | deep |
| SplitAE   | Split Autoencoders [@ngiam2011multimodal] | 2   | deep |
| DVCCA   | Deep Variational Canonical Correlation Analysis | \>=2   | deep |
| ProbabilisticCCA   | Probabilistic Canonical Correlation Analysis | 2   | probabilistic |

## Documentation

The package is accompanied by documentation (https://cca-zoo.readthedocs.io/en/latest/index.html) and a number of
tutorial notebooks which serve as both guides to the package as well as educational resources for CCA and PLS methods.

# Conclusion

`cca-zoo` fills many of the gaps in the multiview learning ecosystem in Python, including a flexible API for
deep-learning based models, regularised models for high dimensional data (and in particular those that induce sparsity),
and generative models.`cca-zoo` will therefore help researchers to apply and develop Canonical Correlation Analysis and
Partial Least Squares models. We continue to welcome contributions from the community.

# Acknowledgements

JC is supported by the EPSRC-funded UCL Centre for Doctoral Training in Intelligent, Integrated Imaging in Healthcare (
i4health) (EP/S021930/1) and the Department of Health’s NIHR-funded Biomedical Research Centre at University College
London Hospitals. HTW is supported by funds from la Fondation Courtois awarded to Dr. Pierre Bellec. 

# References
Tutorials and Examples Gallery
================================

Below is a gallery of examplesWelcome to cca-zoo's documentation!
===================================


Documentation
=============


.. toctree::
   :maxdepth: 1
   :caption: Using cca-zoo

   documentation/install
   documentation/getting_started
   documentation/maths
   documentation/user_guide
   auto_examples/index


.. toctree::
   :maxdepth: 4
   :caption: Reference

   api/data
   api/deepmodels
   api/models
   api/model_selection
   api/probabilisticmodels
   api/utils


.. toctree::
   :maxdepth: 1
   :caption: Developer Info

   developer_info/contribute
   developer_info/support
   developer_info/license


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`




Contribute
==========

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/jameschapman19/cca_zoo/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

cca-zoo has documentation at cca-zoo.readthedocs.io/en/latest/ but feedback/corrections
on areas that are unclear or incorrect would be really helpful. Also any example notebooks can be added in the
tutorial notebooks folder.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/jameschapman19/cca_zoo/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `cca-zoo` for local development.

1. Fork the `cca-zoo` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/cca_zoo.git

3. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

4. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

5. Submit a pull request through the GitHub website.

- Issue Tracker: https://github.com/jameschapman19/cca_zoo/issues
- Source Code: https://github.com/jameschapman19/cca_zooSupport
=========

If you are having issues, please let me know. This is my first python package so I am open to all and any feedback!
james.chapman.19@ucl.ac.ukUser Guide
===========

Model Fit
----------

.. sourcecode:: python

   from cca_zoo.models import CCA
   from cca_zoo.data import generate_covariance_data
   # %%
   (train_view_1,train_view_2),(true_weights_1,true_weights_2)=generate_covariance_data(n=200,view_features=[10,10],latent_dims=1,correlation=1)


   linear_cca = CCA(latent_dims=latent_dims, max_iter=max_iter)

   linear_cca.fit([train_view_1, train_view_2])

Hyperparameter Tuning
^^^^^^^^^^^^^^^^^^^^^^

Some models require hyperparameters. We can either choose these manually when the model is instantiated or we can use the gridsearch_fit() method
to use a data driven approach.

.. sourcecode:: python

   from cca_zoo.models import rCCA
   from cca_zoo.model_selection import GridsearchCV

    def scorer(estimator,X):
      dim_corrs=estimator.score(X)
      return dim_corrs.mean()

    c1 = [0.1, 0.3, 0.7, 0.9]
    c2 = [0.1, 0.3, 0.7, 0.9]
    param_grid = {'c': [c1,c2]}

    ridge = GridSearchCV(rCCA(latent_dims=latent_dims),param_grid=param_grid,
        cv=cv,
        verbose=True,scoring=scorer).fit([train_view_1,train_view_2]).best_estimator_

Model Transforms
-----------------

Once models are fit we can transform the data to latent projections for each view

.. sourcecode:: python

   projection_1,projection_2=ridge.transform([train_view_1,train_view_2])

In a similar way to scikit-learn we can also call fit_transform to complete both of these steps in one go:

.. sourcecode:: python

   projection_1,projection_2=ridge.fit_transform([train_view_1,train_view_2])

Model Evaluation
-----------------

We can evaluate models by their correlation in the latent space

.. sourcecode:: python

   correlation=ridge.score([train_view_1,train_view_2])

For most models this gives us the average pairwise correlation in each latent dimension. For tensor cca models this
gives the higher order correlation in each dimension.

Model Weights
-----------------

In applications of cca, we are often interested in the model weights. These can be easily accessed as arrays with
#features x #latent_dimensions for each view.

.. sourcecode:: python

   view_1_weights=ridge.weights[0]
   view_2_weights=ridge.weights[1]

Model Loadings
-----------------

Similarly we can access the loadings for a given set of samples

.. sourcecode:: python

   view_1_loadings, view_2_loadings=ridge.get_loadings([train_view_1, train_view_2])


Deep Models
------------

Deep models have a slightly more involved process. We first need to choose the architectures for our encoder models

.. sourcecode:: python

   from cca_zoo.deepmodels import architectures
   encoder_1 = architectures.Encoder(latent_dims=latent_dims, feature_size=784)
   encoder_2 = architectures.Encoder(latent_dims=latent_dims, feature_size=784)

We build our deep cca model using these encoders as inputs:

.. sourcecode:: python

   from cca_zoo.deepmodels import DCCA
   dcca_model = DCCA(latent_dims=latent_dims, encoders=[encoder_1, encoder_2])

This produces a PyTorch.nn.Module object which can be updated in a customised training loop. We also provide a LightningModule
class from pytorch-lightning which can be used to train any of these models.Tutorials
=========

Tutorial on CCA and PLS methods
---------------------------------

Introduces the CCA and PLS family of models in an interactive notebook that allows the user to vary parameters in models
with live updates using ipython widgets submitted for MICCAI 2021 Educational Challenge

https://github.com/jameschapman19/cca_zoo/blob/main/tutorial_notebooks/CCA_Tutorial.ipynb

Tutorial running models on paired MNIST halves
--------------------------------------------------

Demonstrates the API for deep learning and multiview (more than 2 view) data

https://github.com/jameschapman19/cca_zoo/blob/main/tutorial_notebooks/cca_zoo_mnist.ipynb

Tutorial understanding how regularisation and sparsity affects model weights for different models
----------------------------------------------------------------------------------------------------

Demonstrates some of the state-of-the-art methods for modelling multiview data with sparse weights

https://github.com/jameschapman19/cca_zoo/blob/main/tutorial_notebooks/cca_zoo_weights_and_sparsity.ipynbMathematical Foundations
===========================

Canonical Correlation Analysis (CCA) and Partial Least Squares (PLS) models
are effective ways of finding associations between multiple views of data.

PCA
----

It is helpful to start off by formulating PCA in its mathematical form.
The first principle component can be written as the solution to the convex optimisation problem:

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{ w_1^TX_1^TX_1w_1  \}

    \text{subject to:}

    w_1^Tw_1=1

WHich is optimized for the singular vectors of the covariance matrix :math:`X^TX`.

PLS
----

Now consider two data matrices with the same number of samples :math:`X_1` and :math:`X_2`.
It is tempting to write a slightly different optimisation problem:

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{ w_1^TX_1^TX_2w_2  \}

    \text{subject to:}

    w_1^Tw_1=1

    w_2^Tw_2=1

Which is optimised for the left and right singular vectors of the cross covariance matrix :math:`X_1^TX_2`


CCA
----

To arrive at Canonical Correlation we change the constraints slightly so that they now depend on the variance
in each view as well as the weights. This makes them data dependent.

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{ w_1^TX_1^TX_2w_2  \}

    \text{subject to:}

    w_1^TX_1^TX_1w_1=1

    w_2^TX_2^TX_2w_2=1

Despite these more complex constraints, CCA is the solution to a generalized eigenvalue problem.

Regularized CCA
-----------------

Notice that the constraints for PLS and CCA are identical in the case that :math:`X_i^TX_i=I`.
This leads naturally to mixing the constraints and in particular the PLS constraint acts in a similar
way to the 'ridge' regularisation in Ridge Regression.

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{ w_1^TX_1^TX_2w_2  \}

    \text{subject to:}

    (1-c_1)w_1^TX_1^TX_1w_1+c_1w_1^Tw_1=1

    (1-c_2)w_2^TX_2^TX_2w_2+c_2w_2^Tw_2=1

Other regularized CCA and PLS
--------------------------------

There has been lots of research into more general forms of regularisation for CCA, in particular forms that induce
sparsity on the weights. In general these can be written in a form somewhat similar to:

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{ w_1^TX_1^TX_2w_2 + \lambda_1f_1(w_1) + \lambda_2f_2(w_2)  \}

    \text{subject to:}

    (1-c_1)w_1^TX_1^TX_1w_1+c_1w_1^Tw_1=1

    (1-c_2)w_2^TX_2^TX_2w_2+c_2w_2^Tw_2=1

These problems usually have no closed form solution and are typically solved by alternately fixing :math:`w_1`
and :math:`w_2`.

Kernel CCA and PLS
---------------------

By expressing the weights as :math:`w_i=\alpha_iX_i` we can transform the CCA problem into its kernel form:

.. math::

    \alpha_{opt}=\underset{\alpha}{\mathrm{argmax}}\{ \alpha_1^TK_1^TK_2\alpha_2  \}

    \text{subject to:}

    \alpha_1^TK_1^TK_1\alpha_1=1

    \alpha_2^TK_2^TK_2\alpha_2=1

Finally we can also consider more general kernel matrices without loss of generalisation. A similar reparameterization
exists for PLS and therefore for ridge regularized CCA.

Multiset CCA
----------------

Ridge Regularized CCA can be generalized to find correlations between more than one view with the formulation:

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{\sum_i\sum_{j\neq i} w_i^TX_i^TX_jw_j  \}\\

    \text{subject to:}

    (1-c_i)w_i^TX_i^TX_iw_i+c_iw_i^Tw_i=1

This form is often referred to as SUMCOR CCA because it optimizes for the pairwise sum of . This can be formulated as a
generalized eigenvalue problem and therefore solved efficiently.

Generalized CCA
-----------------

Ridge Regularized CCA can be generalized to find correlations between more than one view with the formulation:

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{ \sum_iw_i^TX_i^TT  \}\\

    \text{subject to:}

    T^TT=1

This form is often referred to as MAXVAR CCA since it finds an auxiliary vector :math:`T` with fixed unit norm that has
maximum sum of variance with each view. This can also be formulated as a generalized eigenvalue problem and
therefore solved efficiently.

Deep CCA
----------

The ideas behind CCA can be extended to a general form where instead of linear weights, we consider functions
:math:`f(X_i)`. Where these functions are parameterized by neural networks, we can consider a 'Deep' CCA:

.. math::

    w_{opt}=\underset{w}{\mathrm{argmax}}\{ f(X_1)^Tf(X_2)  \}

    \text{subject to:}

    f(X_1)^Tf(X_1)=1

    f(X_2)^Tf(X_2)=1Getting Started
===============

cca-zoo is a collection of linear, kernel, and deep methods for canonical correlation analysis of multiview data.
Where possible I have followed the scikit-learn/mvlearn APIs and models therefore have
fit/transform/fit_transform methods as standard.

Look how easy it is to use:

.. sourcecode:: python

   from cca_zoo.models import CCA
   from cca_zoo.data import generate_covariance_data
   # %%
   n_samples=100
   (train_view_1,train_view_2),(true_weights_1,true_weights_2)=generate_covariance_data(n=n_samples,view_features=[10,10],latent_dims=1,correlation=1)

   linear_cca = CCA(latent_dims=latent_dims, max_iter=max_iter)

   linear_cca.fit((train_view_1, train_view_2))

Installation
=============

cca-zoo can be installed with python versions >=3.7. At the moment there are two options for installation:

1. Clone `cca-zoo` and create a virtual environment using the requirements
3. Install directly from PyPI release without cloning `cca-zoo`.

Installing with pip
----------------------------------------

Install the current release of ``cca-zoo`` with ``pip``::

    $ pip install cca-zoo

To upgrade to a newer release use the ``--upgrade`` flag::

    $ pip install --upgrade cca-zoo

Optional dependencies
----------------------------------------

Since some of the functionality depends on PyTorch and Numpyro which are both large packages and may have difficulties
with windows installation, we do not install them by default. To access these,

* [deep]: ``Deep Learning Based Models``
* [probabilistic]: ``Probabilistic Models``
* [all]: ``Include both Probabilistic and Deep Learning Based Models``

If you wish to use these functions, you must install their required dependencies. These are listed in the package requirements folder with corresponding keyword names for manual installation or can be installed from PyPI by simply calling::

    $ pip install cca-zoo[keyword]

where 'keyword' is from the list above, bracketed.
To upgrade the package and torch requirements::

    $ pip install --upgrade cca-zoo[keyword]

Hardware requirements
---------------------
The ``cca-zoo`` package has no specific hardware requirements but some models may be RAM intensive and deep learning and probabilistic models may benefit from a dedicated GPU

OS Requirements
---------------
This package is supported for *Linux* and *macOS* and can also be run on Windows machines.
Probabilistic Models
====================================

Variational CCA
----------------------------------------

.. automodule:: cca_zoo.probabilisticmodels.probabilisticcca
    :inherited-members:
    :exclude-members: get_params, set_params
Deep Models
===========================

DCCA
-------------------------------

.. automodule:: cca_zoo.deepmodels.dcca
    :members:
    :undoc-members:

DCCA by Non-Linear Orthogonal Iterations
-----------------------------------------

.. automodule:: cca_zoo.deepmodels.dcca_noi
    :members:
    :undoc-members:

Deep Canonically Correlated Autoencoders
-----------------------------------------

.. automodule:: cca_zoo.deepmodels.dccae
    :members:
    :undoc-members:

Deep Tensor CCA
--------------------------------

.. automodule:: cca_zoo.deepmodels.dtcca
    :members:
    :undoc-members:

Deep Variational CCA
--------------------------------

.. automodule:: cca_zoo.deepmodels.dvcca
    :members:
    :undoc-members:

Split Autoencoders
----------------------------------

.. automodule:: cca_zoo.deepmodels.splitae
    :members:
    :undoc-members:

CCALightning Module for training with pytorch-lightning
---------------------------------------------------------

.. automodule:: cca_zoo.deepmodels.CCALightning
    :members:
    :inherited-members:
    :exclude-members: get_params, set_params

Deep Objectives
-------------------------------------

.. autoclass:: cca_zoo.deepmodels.objectives.CCA
    :members:

.. autoclass:: cca_zoo.deepmodels.objectives.MCCA
    :members:

.. autoclass:: cca_zoo.deepmodels.objectives.GCCA
    :members:

.. autoclass:: cca_zoo.deepmodels.objectives.TCCA
    :members:


Model Architectures
----------------------------------------

.. automodule:: cca_zoo.deepmodels.architectures
Data
=====================

Simulated Data
------------------------------

.. automodule:: cca_zoo.data.simulated
    :members:

Toy Data
------------------------

.. automodule:: cca_zoo.data.toy
    :members:

Utils
--------------------------

.. automodule:: cca_zoo.data.utils
    :members:

Model Selection
==========================

GridSearchCV
---------------------------------------

.. automodule:: cca_zoo.model_selection.GridSearchCV
    :inherited-members:


RandomizedSearchCV
---------------------------------------------

.. automodule:: cca_zoo.model_selection.RandomizedSearchCV
    :inherited-members:

Models
=======================


Regularized Canonical Correlation Analysis and Partial Least Squares
------------------------------------------------------------------------

Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.rcca.CCA
    :inherited-members:
    :exclude-members: get_params, set_params

Partial Least Squares
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.rcca.PLS
    :inherited-members:
    :exclude-members: get_params, set_params

Ridge Regularized Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.rcca.rCCA
    :inherited-members:
    :exclude-members: get_params, set_params

GCCA and KGCCA
---------------------------

Generalized (MAXVAR) Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.gcca.GCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Kernel Generalized (MAXVAR) Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.gcca.KGCCA
    :inherited-members:
    :exclude-members: get_params, set_params

MCCA and KCCA
---------------------------

Multiset (SUMCOR) Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.mcca.MCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Kernel Multiset (SUMCOR) Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.mcca.KCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Tensor Canonical Correlation Analysis
----------------------------------------

Tensor Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.tcca.TCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Kernel Tensor Canonical Correlation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.tcca.KTCCA
    :inherited-members:
    :exclude-members: get_params, set_params

More Complex Regularisation using Iterative Models
-----------------------------------------------------

Normal CCA and PLS by alternating least squares
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Quicker and more memory efficient for very large data


CCA by Alternating Least Squares
""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.CCA_ALS
    :inherited-members:
    :exclude-members: get_params, set_params

PLS by Alternating Least Squares
""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.PLS_ALS
    :inherited-members:
    :exclude-members: get_params, set_params


Sparsity Inducing Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Penalized Matrix Decomposition (Sparse PLS)
"""""""""""""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.PMD
    :inherited-members:
    :exclude-members: get_params, set_params

Sparse CCA by iterative lasso regression
"""""""""""""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.SCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Elastic CCA by MAXVAR
"""""""""""""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.ElasticCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Span CCA
"""""""""""""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.SpanCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Parkhomenko (penalized) CCA
"""""""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.ParkhomenkoCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Sparse CCA by ADMM
"""""""""""""""""""""""""""""""""""""""""""
.. autoclass:: cca_zoo.models.SCCA_ADMM
    :inherited-members:
    :exclude-members: get_params, set_params

Miscellaneous
-----------------------------------------------------

Nonparametric CCA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.NCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Partial CCA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.PartialCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Sparse Weighted CCA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: cca_zoo.models.SWCCA
    :inherited-members:
    :exclude-members: get_params, set_params

Base Class
--------------------------------

.. automodule:: cca_zoo.models._cca_base
    :members:
    :private-members: _CCA_Base
    :exclude-members: get_params, set_params
