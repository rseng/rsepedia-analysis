---
title: 'Delve: Neural Network Feature Variance Analysis'
tags:
  - Python
  - deep learning
  - machine learning
  - saturation
  - pytorch
  - AI
authors:
  - name: Justin Shenk^[co-first author]
    orcid: 0000-0002-0664-7337
    affiliation: "1,2"
  - name: Mats L. Richter^[co-first author]
    affiliation: 2
    orcid: 0000-0002-9525-9730
  - name: Wolf Byttner
    affiliation: 3
    orcid: 0000-0002-9525-9730
affiliations:
 - name: VisioLab, Berlin, Germany
   index: 1
 - name: Institute of Cognitive Science, University of Osnabrueck, Osnabrueck, Germany
   index: 2
 - name: Rapid Health, London, England, United Kingdom
   index: 3
date: 16 August 2021
bibliography: paper.bib
---

# Summary
Designing neural networks is a complex task.
Deep neural networks are often referred to as "black box" models - little insight in the function they approximate is gained from looking at the structure of layer outputs.
``Delve`` is a tool for looking at how a neural network represents data, and how these representations, or features, change throughout training.
This tool enables deep learning researchers to understand the limitations and suggest improvements for the design of their networks, such as removing or adding layers.

Several tools exist which allow analyzing neural networks after and during training.
These techniques can be characterized by their focus on either data or model as well as their level of abstractness.
Examples for abstract model-oriented techniques are tools for analyzing the sharpness of local optima [@keskar;@sensitivitygoogle], which can be an indicator for the generalizing capabilities of the trained models.
In these scenarios the complexity of the dataset and model is reduced to the error surface, allowing for insights into the differences between different setups.
A less abstract data-centric technique GradCam by Selvaraju et al. [@gradcam;@gradcamplusplus], reduces the model to a set of class-activation maps that can be overlayed over individual data points to get an intuitive understanding of the inference process.
SVCCA [@svcca;@svcca2] can be considered model-centric and a middle ground in terms of abstractness, since it allows the comparative analysis of the features extracted by specific layers.
SVCCA is also relevant from a functional perspective for this work, since it uses singular value decomposition as a core technique to obtain the analysis results.
Another model-centric tool that allows for a layer-by-layer analysis is logistic regression probes [@alain2016], which utilize logistic regressions trained on the output of a hidden layer to measure the linear separability of the data and thus the quality of the intermediate solution quality of a classifier model.

The latter is of great importance for this work since logistic regression probes are often used to compare models and identify the contribution of layers to overall performance [@feature-space;@sizematters;@goingdeeper] and to demonstrate that the saturation metric is capable of showing parameter-inefficiencies in neural network architectures.

However, the aforementioned tools have significant limitations in terms of their usefulness in practical application scenarios, where these tools are to be used to improve the performance of a given model.
In the case of data-centric tools like GradCam, the solution propagates back to the data, which makes it hard to derive decisions regarding the neural architecture.
However, the biggest concern in all aforementioned tools is the cost of computational resources and the integration of the analysis into the workflow of a deep learning practitioner.
Tools like SVCCA and logistic regression probes require complex and computationally expensive procedures that need to be conducted after training.
This naturally limits these techniques to small benchmarks and primarily academic datasets like Cifar10 [@feature-space].
An analysis tool that is to be used during the development of a deep learning-based model needs to be able to be used with as little computational and workflow overhead as possible.
Ideally, the analysis can be done live while the training is in progress, allowing the researcher to interrupt potentially long-running training sessions to improve the model.
Saturation was proposed in 2018 [@Shenk:Thesis:2018] and later refined [@feature-space] and is the only analysis technique known to the authors that has this capability while allowing to identify parameter-inefficiencies in the setup [@feature-space;@sizematters;@goingdeeper].
To make saturation usable in an application scenario, it is necessary to provide an easy-to-use framework that allows for an integration of the tool into the normal training and inference code with only minimally invasive changes.
It is also necessary that the computation and analysis can be done online as part of the regular forward pass of the model, to make the integration as seamless as possible.
A numerical comparison of these various methods is a promising avenue for future research into model introspection.

``Delve`` is a tool for extracting information based on the covariance matrix of the data like saturation and the intrinsic dimensionality from neural network layers.
To emphasize practical usability, special attention is placed on a low overhead and minimally invasive integration of ``Delve`` into
existing training and inference setups:
``Delve`` hooks directly into PyTorch [@pytorch] models to extract necessary information with little computational and memory overhead, thanks to an efficient covariance approximation algorithm.
We enable the user to store and analyze the extracted statistics without changing their current experiment workflow, by making ``Delve`` easy to integrate into monitoring systems and making this interface easy to expand.
This allows the user to utilize their preferred way of monitoring experiments, from simple CSV-Files and folder structures to more sophisticated solutions like
TensorBoard [@tensorflow2015-whitepaper].
A comprehensive source of documentation is provided on the homepage
([http://delve-docs.readthedocs.io](delve-docs.readthedocs.io)).


## Statement of Need
Research on spectral properties of neural network representations has exploded in recent years [@svcca;@svcca2;@gradcam;@kernelPCA;@alain2016;@featureAttribution].
Publications like [@svcca] and [@feature-space] demonstrate that useful and interesting information can be extracted from the spectral analysis of these latent representations.
It has also been shown that metrics like saturation [@Shenk:Thesis:2018;@spectral-analysis] can be used to optimize neural network architectures by identifying pathological patterns hinting at inefficiencies of the neural network structure.


The main purpose of ``Delve`` is to provide easy and flexible access to these types of layer-based statistics.
The combination of ease of usage and extensibility in ``Delve`` enables exciting scientific explorations for machine learning researchers and engineers.
``Delve`` has already been used in a number of scientific publications [@feature-space;@sizematters;@goingdeeper].
The source code for ``Delve`` has been archived to Zenodo with the linked DOI: [@zenodo]



## Overview of the Library
The software is structured into several modules which distribute tasks. Full details are available at <https://delve-docs.readthedocs.io/>.

The TensorBoardX `SummaryWriter` [@tensorflow2015-whitepaper] is used to efficiently save artifacts like images or statistics during training with minimal interruption.
A variety of layer feature statistics can be observed:

| Statistic |
|----------------------------------------------------------------------------------|
| intrinsic dimensionality                                                               |
| layer saturation (intrinsic dimensionality divided by feature space dimensionality      |
| the covariance-matrix                    |
| the determinant of the covariance matrix (also known as generalized variance)          |
| the trace of the covariance matrix, a measure of the variance of the data                  |
| the trace of the diagonal matrix, another way of measuring the dispersion of the data. |
| layer saturation (intrinsic dimensionality divided by feature space dimensionality)    |


Several layers are currently supported:

* Convolutional
* Linear
* LSTM

Additional layers such as PyTorch's ConvTranspose2D are planned for future development (see issue [#43](https://github.com/delve-team/delve/issues/43)).

## Eigendecomposition of the feature covariance matrix
The computation of saturation and other related metrics like the intrinsic dimensionality requires the covariance matrix of the layer's output.
Computing the covariance matrix of a layer's output on the training or evaluation set is impractical to do naively, since it would require holding the entire dataset in memory.
This would also contradict our goal of seamless integration in existing training loops, which commonly operate with mini-batches.
Therefore, a batch-wise approximation algorithm is used to compute the covariance matrix online during training:

We can compute the covariance between two variables by using the covariance approximation algorithm for two random variables $X$ and $Y$ with $n$ samples:
$$Q(X, Y) = \frac{\sum^{n}_{i=1} x_i y_i}{n} - \frac{(\sum^{n}_{i=1} x_i)  (\sum^{n}_{i=1} y_i)}{n^2}$$
Where $x_i$ and $y_i$ are individual observations of the respective random variables $X$ and $Y$ and $n$ is the total number of samples.
The advantage of this method is that only the number of seen samples, the sum of squares and the sum of the variables need to be stored,
making the memory consumption per layer constant with respect to the size of the dataset.
By computing Q(X, Y) for all possible combinations of features, we obtain the covariance matrix of the layer's output $Q(Z_l, Z_l)$, where $Z_l$ is the layer's output over
the entire dataset.
We can parallelize the computations of all feature combinations by exploiting the shape of the layer output matrix $A_l$ of the layer $l$:
We can compute $\sum^{n}_{i=1} x_i y_i$ for all feature combinations in layer $l$ by calculating the running squares $\sum^{B}_{b=0}A_{l,b}^T A_{l,b}$ of the batch output matrices $A_{l,b}$ where $b \in \{0,...,B-1\}$ for $B$ batches. We replace $\frac{(\sum^{n}_{i=1} x_i)  (\sum^{n}_{i=1} y_i)}{n^2}$ by the outer product $\bar{A}_l \bigotimes \bar{A}_l$ of the sample mean $\bar{A}_l$.
This is the running sum of all outputs $z_{l,k}$, where $k \in \{0,...,n\}$ at training time, divided by the total number of training samples $n$.
Our formula for a batch-wise approximated covariance matrix can now be written like this:
$$Q(Z_l, Z_l) = \frac{\sum^{B}_{b=0}A_{l,b}^T A_{l,b}}{n} -(\bar{A}_l \bigotimes \bar{A}_l)$$
The batch-wise updating algorithm allows us to integrate the approximation of the covariance matrix as part of the regular forward pass during training and evaluation.
Our algorithm uses a thread-safe common value store on a single compute device or node, which furthermore allows updating the covariance matrix asynchronously when the network is trained in a distributed manner.
To avoid problems that can be caused by rounding errors and numerical instability, our implementation of the algorithm converts by default all data into 64-bit floating-point values by default.

Another challenge is the dimensionality of the data in convolutional layers, where a simple flattening of the data vector would result in a very high dimensional vector and a computationally expensive singular value decomposition as a direct consequence. To address this issue, we treat every kernel position as an individual observation. This turns a 4th-degree output-tensor of shape (samples, height,  width, filters) into a matrix of shape (samples $\cdot$ height $\cdot$ width, filters).
The advantage of this strategy is that no information is lost, while keeping the dimensionality of $Q$ at a manageable size.
Optionally, to reduce the computations required further, the feature map can be automatically reduced in size using linear interpolation to a constant maximum height and width. Since information is lost during this process, 
this is disabled.

This approximation method was described alongside the saturation metric in the works of [@Shenk:Thesis:2018;@spectral-analysis] and further refined by [@feature-space].

# References
0.1.32 (2020-1-29)
------------------
* Add LSTM, VAE support

0.1.4 (2018-6-19)
-----------------
* Fix bug due to incomplete merge
  
0.1.3 (2018-6-17)
-----------------
* Add support for convolutional layers

0.1.2 (2018-6-14)
-----------------
* Add iPython notebook support

0.1.1 (2018-6-13)
-----------------
* Add tqdm output

0.1 (2018-06-11)
------------------
* First commit
Delve: Deep Live Visualization and Evaluation |logo|
====================================================

|PyPI version| |Tests| |codecov.io| |License: MIT| |DOI| 

Delve is a Python package for analyzing the inference dynamics of your model.

.. image:: https://raw.githubusercontent.com/justinshenk/playground/master/saturation_demo.gif
   :alt: playground

Use Delve if you need a lightweight PyTorch extension that: 

- Gives you insight into the inference dynamics of your architecture 
- Allows you to optimize and adjust neural networks models to your dataset without much trial and error 
- Allows you to analyze the eigenspaces your data at different stages of inference 
- Provides you basic tooling for experiment logging

Motivation
----------

Designing a deep neural network is a trial and error heavy process that
mostly revolves around comparing performance metrics of different runs.
One of the key issues with this development process is that the results
of metrics not really propagate back easily to concrete design
improvements. Delve provides you with spectral analysis tools that allow
you to investigate the inference dynamic evolving in the model while
training. This allows you to spot underutilized and unused layers.
Mismatches between object size and neural architecture among other
inefficiencies. These observations can be propagated back directly to
design changes in the architecture even before the model has fully
converged, allowing for a quicker and more guided design process.

This work is closely related to Maithra Raghu (Google Brain) et al's work on SVCCA:

- "Maithra Raghu on the differences between wide and deep networks", 2020 `[YouTube] <https://youtu.be/6uPop547u_E?t=970>`_
- "SVCCA:Singular Vector Canonical Correlation Analysis for Deep Learning and Interpretability", 2017 `[arXiv] <https://arxiv.org/abs/1706.05806>`_
  
Installation
------------

.. code:: bash

   pip install delve

Using Layer Saturation to improve model performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The saturation metric is the core feature of delve. By default
saturation is a value between 0 and 1.0 computed for any convolutional,
lstm or dense layer in the network. The saturation describes the
percentage of eigendirections required for explaining 99% of the
variance. Simply speaking, it tells you how much your data is “filling
up” the individual layers inside your model.

In the image below you can see how saturation portraits inefficiencies
in your neural network. The depicted model is ResNet18 trained on 32
pixel images, which is way to small for a model with a receptive field
exceeding 400 pixels in the final layers.

.. image:: https://raw.githubusercontent.com/delve-team/delve/master/images/resnet.PNG
   :alt: resnet.PNG

To visualize what this poorly chosen input resolution does to the
inference, we trained logistic regressions on the output of every layer
to solve the same task as the model. You can clearly see that only the
first half of the model (at best) is improving the intermedia solutions
of our logistic regression “probes”. The layers following this are
contributing nothing to the quality of the prediction! You also see that
saturation is extremly low for this layers!

We call this a *tail* and it can be removed by either increasing the
input resolution or (which is more economical) reducing the receptive
field size to match the object size of your dataset.

.. figure:: https://raw.githubusercontent.com/delve-team/delve/master/images/resnetBetter.PNG
   :alt: resnetBetter.PNG

We can do this by removing the first two downsampling layers, which
quarters the growth of the receptive field of your network, which
reduced not only the number of parameters but also makes more use of the
available parameters, by making more layers contribute effectivly!

**For more details check our publication on this topics** - `Spectral
Analysis of Latent Representations <https://arxiv.org/abs/1907.08589>`__
- `Feature Space Saturation during
Training <https://arxiv.org/abs/2006.08679>`__ - `(Input) Size Matters
for CNN
Classifiers <https://link.springer.com/chapter/10.1007/978-3-030-86340-1_11>`__
- `Should you go deeper? Optimizing Convolutional Neural Networks
without training <https://arxiv.org/abs/2106.12307>`__ - Go with the
Flow: the distribution of information processing in multi-path networks
(soon)

Demo
----

.. code:: python


   import torch
   from delve import SaturationTracker
   from torch.cuda import is_available
   from torch.nn import CrossEntropyLoss
   from torchvision.datasets import CIFAR10
   from torchvision.transforms import ToTensor, Compose
   from torch.utils.data.dataloader import DataLoader
   from torch.optim import Adam
   from torchvision.models.vgg import vgg16

   # setup compute device
   from tqdm import tqdm

   if __name__ == "__main__":

     device = "cuda:0" if is_available() else "cpu"

     # Get some data
     train_data = CIFAR10(root="./tmp", train=True,
                          download=True, transform=Compose([ToTensor()]))
     test_data = CIFAR10(root="./tmp", train=False, download=True, transform=Compose([ToTensor()]))

     train_loader = DataLoader(train_data, batch_size=1024,
                               shuffle=True, num_workers=6,
                               pin_memory=True)
     test_loader = DataLoader(test_data, batch_size=1024,
                              shuffle=False, num_workers=6,
                              pin_memory=True)

     # instantiate model
     model = vgg16(num_classes=10).to(device)

     # instantiate optimizer and loss
     optimizer = Adam(params=model.parameters())
     criterion = CrossEntropyLoss().to(device)

     # initialize delve
     tracker = SaturationTracker("my_experiment", save_to="plotcsv", modules=model, device=device)

     # begin training
     for epoch in range(10):
       model.train()
       for (images, labels) in tqdm(train_loader):
         images, labels = images.to(device), labels.to(device)
         prediction = model(images)
         optimizer.zero_grad(set_to_none=True)
         with torch.cuda.amp.autocast():
           outputs = model(images)
           _, predicted = torch.max(outputs.data, 1)

           loss = criterion(outputs, labels)
         loss.backward()
         optimizer.step()

       total = 0
       test_loss = 0
       correct = 0
       model.eval()
       for (images, labels) in tqdm(test_loader):
         images, labels = images.to(device), labels.to(device)
         outputs = model(images)
         loss = criterion(outputs, labels)
         _, predicted = torch.max(outputs.data, 1)

         total += labels.size(0)
         correct += torch.sum((predicted == labels)).item()
         test_loss += loss.item()

       # add some additional metrics we want to keep track of
       tracker.add_scalar("accuracy", correct / total)
       tracker.add_scalar("loss", test_loss / total)

       # add saturation to the mix
       tracker.add_saturations()

     # close the tracker to finish training
     tracker.close()

Supported Layers
----------------

* Dense/Linear
* LSTM
* Convolutional

Citation
--------

If you use Delve in your publication, please cite:

.. code-block:: txt

   @software{delve,
   author       = {Justin Shenk and
                     Mats L. Richter and
                     Wolf Byttner and
                     Michał Marcinkiewicz},
   title        = {delve-team/delve: Latest},
   month        = aug,
   year         = 2021,
   publisher    = {Zenodo},
   version      = {v0.1.49},
   doi          = {10.5281/zenodo.5233859},
   url          = {https://doi.org/10.5281/zenodo.5233859}
   }


Why this name, Delve?
~~~~~~~~~~~~~~~~~~~~~

**delve** (*verb*):

-  reach inside a receptacle and search for something
-  to carry on intensive and thorough research for data, information, or
   the like

.. |logo| image:: https://github.com/delve-team/delve/blob/master/images/delve_logo.png
.. |PyPI version| image:: https://badge.fury.io/py/delve.svg
   :target: https://badge.fury.io/py/delve
.. |Tests| image:: https://github.com/delve-team/delve/actions/workflows/tests.yaml/badge.svg
   :target: https://github.com/delve-team/delve/actions/workflows/tests.yaml
.. |codecov.io| image:: https://codecov.io/github/delve-team/delve/coverage.svg?branch=master
   :target: https://codecov.io/github/delve-team/delve/?branch=master
.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
.. |DOI| image:: https://zenodo.org/badge/136951823.svg
   :target: https://zenodo.org/badge/latestdoi/136951823
example\_keras module
=====================

.. automodule:: example_keras
   :members:
   :undoc-members:
   :show-inheritance:
example\_fc module
==================

.. automodule:: example_fc
   :members:
   :undoc-members:
   :show-inheritance:
example\_deep module
====================

.. automodule:: example_deep
   :members:
   :undoc-members:
   :show-inheritance:
example module
==============

.. automodule:: example
   :members:
   :undoc-members:
   :show-inheritance:
example\_lstm\_class\_ae module
===============================

.. automodule:: example_lstm_class_ae
   :members:
   :undoc-members:
   :show-inheritance:
example\_lstm\_generative\_vae module
=====================================

.. automodule:: example_lstm_generative_vae
   :members:
   :undoc-members:
   :show-inheritance:
Contributing to Delve
=====================

(Contribution guidelines largely copied from `geopandas <https://geopandas.readthedocs.io/en/latest/contributing.html>`_)

Overview
--------

Contributions to Delve are very welcome.  They are likely to
be accepted more quickly if they follow these guidelines.

At this stage of Delve development, the priorities are to define a
simple, usable, and stable API and to have clean, maintainable,
readable code. Performance matters, but not at the expense of those
goals.

In general, Delve follows the conventions of the pandas project
where applicable.

In particular, when submitting a pull request:

- All existing tests should pass.  Please make sure that the test
  suite passes, both locally and on
  `GitHub Actions <https://github.com/delve-team/delve/actions/workflows/tests.yaml>`_.  Status on
  GitHub Actions will be visible on a pull request.

- New functionality should include tests.  Please write reasonable
  tests for your code and make sure that they pass on your pull request.

- Classes, methods, functions, etc. should have docstrings.  The first
  line of a docstring should be a standalone summary.  Parameters and
  return values should be ducumented explicitly.

- Delve supports python 3 (3.6+).  Use modern python idioms when possible.

- Follow PEP 8 when possible.

- Imports should be grouped with standard library imports first,
  3rd-party libraries next, and Delve imports third.  Within each
  grouping, imports should be alphabetized.  Always use absolute
  imports when possible, and explicit relative imports for local
  imports when necessary in tests.


Seven Steps for Contributing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are seven basic steps to contributing to *Delve*:

1) Fork the *Delve* git repository
2) Create a development environment
3) Install *Delve* dependencies
4) Make a ``development`` build of *Delve*
5) Make changes to code and add tests
6) Update the documentation
7) Submit a Pull Request

Each of these 7 steps is detailed below.


1) Forking the *Delve* repository using Git
------------------------------------------------

To the new user, working with Git is one of the more daunting aspects of contributing to *Delve*.
It can very quickly become overwhelming, but sticking to the guidelines below will help keep the process
straightforward and mostly trouble free.  As always, if you are having difficulties please
feel free to ask for help.

The code is hosted on `GitHub <https://github.com/delve-team/delve>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* Software Carpentry's `Git Tutorial <http://swcarpentry.github.io/git-novice/>`_
* `Atlassian <https://www.atlassian.com/git/tutorials/what-is-version-control>`_
* the `GitHub help pages <http://help.github.com/>`_.
* Matthew Brett's `Pydagogue <http://matthew-brett.github.com/pydagogue/>`_.

Getting started with Git
~~~~~~~~~~~~~~~~~~~~~~~~~

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
~~~~~~~~

You will need your own fork to work on the code. Go to the `Delve project
page <https://github.com/delve-team/delve>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone git@github.com:your-user-name/delve.git delve-yourname
    cd delve-yourname
    git remote add upstream git://github.com/delve-team/delve.git

This creates the directory `delve-yourname` and connects your repository to
the upstream (main project) *Delve* repository.

The testing suite will run automatically on Travis-CI once your pull request is
submitted.  However, if you wish to run the test suite on a branch prior to
submitting the pull request, then Travis-CI needs to be hooked up to your
GitHub repository.  Instructions for doing so are `here
<http://about.travis-ci.org/docs/user/getting-started/>`__.

Creating a branch
~~~~~~~~~~~~~~~~~~

You want your master branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *delve*. You can have many shiny-new-features
and switch in between them using the git checkout command.

To update this branch, you need to retrieve the changes from the master branch::

    git fetch upstream
    git rebase upstream/master

This will replay your commits on top of the latest Delve git master.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes, you will need to ``stash`` them prior
to updating.  This will effectively store your changes and they can be reapplied
after updating.

.. _contributing.dev_env:

2) Creating a development environment
---------------------------------------
A development environment is a virtual space where you can keep an independent installation of *Delve*.
This makes it easy to keep both a stable version of python in one place you use for work, and a development
version (which you may break while playing with code) in another.

An easy way to create a *Delve* development environment is as follows:

- Install either `Anaconda <http://docs.continuum.io/anaconda/>`_ or
  `miniconda <http://conda.pydata.org/miniconda.html>`_
- Make sure that you have :ref:`cloned the repository <contributing.forking>`
- ``cd`` to the *delve* source directory

Tell conda to create a new environment, named ``delve_dev``, or any other name you would like
for this environment, by running::

      conda create -n delve_dev

For a python 3 environment::

      conda create -n delve_dev python=3.8

This will create the new environment, and not touch any of your existing environments,
nor any existing python installation.

To work in this environment, Windows users should ``activate`` it as follows::

      activate delve_dev

Mac OSX and Linux users should use::

      source activate delve_dev

You will then see a confirmation message to indicate you are in the new development environment.

To view your environments::

      conda info -e

To return to you home root environment::

      deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

At this point you can easily do a *development* install, as detailed in the next sections.

3) Installing Dependencies
--------------------------

To run *Delve* in an development environment, you must first install
*Delve*'s dependencies. We suggest doing so using the following commands
(executed after your development environment has been activated)::

    pip install -r requirements/reqirements.txt

This should install all necessary dependencies.

Next activate pre-commit hooks by running::

    pre-commit install

4) Making a development build
-----------------------------

Once dependencies are in place, make an in-place build by navigating to the git
clone of the *delve* repository and running::

    python setup.py develop


5) Making changes and writing tests
-------------------------------------

*Delve* is serious about testing and strongly encourages contributors to embrace
`test-driven development (TDD) <http://en.wikipedia.org/wiki/Test-driven_development>`_.
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

Adding tests is one of the most common requests after code is pushed to *delve*.  Therefore,
it is worth getting in the habit of writing tests ahead of time so this is never an issue.

*delve* uses the `pytest testing system
<http://doc.pytest.org/en/latest/>`_ and the convenient
extensions in `numpy.testing
<http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.

Writing tests
~~~~~~~~~~~~~

All tests should go into the ``tests`` directory. This folder contains many
current examples of tests, and we suggest looking to these for inspiration.


Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

The tests can then be run directly inside your Git clone (without having to
install *Delve*) by typing::

    pytest

6) Updating the Documentation
-----------------------------

*Delve* documentation resides in the `doc` folder. Changes to the docs are
make by modifying the appropriate file in the `source` folder within `doc`.
*Delve* docs us reStructuredText syntax, `which is explained here <http://www.sphinx-doc.org/en/stable/rest.html#rst-primer>`_
and the docstrings follow the `Numpy Docstring standard <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

Once you have made your changes, you can build the docs by navigating to the `doc` folder and typing::

    make html

The resulting html pages will be located in `doc/build/html`.


7) Submitting a Pull Request
------------------------------

Once you've made changes and pushed them to your forked repository, you then
submit a pull request to have them integrated into the *Delve* code base.

You can find a pull request (or PR) tutorial in the `GitHub's Help Docs <https://help.github.com/articles/using-pull-requests/>`_.
.. _Saturation Overview:

Saturation
==========

`Saturation` is a metric used for identifying the intrinsic dimensionality of
features in a layer.

A visualization of how saturation changes over training and can be used to optimize network topology is provided at https://github.com/justinshenk/playground:

.. image:: _static/saturation_demo.gif

Covariance matrix of features is computed online:

.. math::

    Q(Z_l, Z_l) = \frac{\sum^{B}_{b=0}A_{l,b}^T A_{l,b}}{n} -(\bar{A}_l \bigotimes \bar{A}_l)

for :math:`B` batches of layer output matrix :math:`A_l` and :math:`n` number of samples.

.. note::

    For more information about how saturation is computed, read `"Feature Space Saturation during Training" <https://arxiv.org/abs/2006.08679>`_.
Examples
========

Delve allows the user to create plots and log records in various formats.

Plotting
--------

Delve allows plotting results every epoch using ``save_to="csvplot"``, which will create automated plots from the metrics
recorded in the ``stats`` argument. The plots depict the layers generally in order of the forward pass.

.. figure:: gallery/images/VGG16-Cifar10-r32-bs256-e90idim_epoch_88.png
  :width: 400
Automatically generated plot of intrinsic dimensionality computed on the training set of Cifar10 on  VGG16 at the 88th epoch of a 90 epoch of training.

.. figure:: gallery/images/VGG16-Cifar10-r32-bs256-e90lsat_epoch_88.png
  :width: 400
Automatically generated plot of saturation computed on the training set of Cifar10 on  VGG16 at the 88th epoch of a 90 epoch training.


Logging
-------

Delve logs results with the ``logging`` package and shows progress with ``tqdm``.

.. figure:: gallery/images/logging.JPG
  :width: 400

A simple example generated from a two-layer network trained on randomly generated data is provided in :ref:`sphx_glr_gallery`.
Reference
=============

SaturationTracker
-------------

``SaturationTracker`` provides a hook for PyTorch and extracts metrics during model training.

.. autoclass:: delve.SaturationTracker


API Pages
---------

.. currentmodule:: delve
.. autosummary::
  :template: autosummary.rst
  :toctree: reference/

  SaturationTracker
Academic Gallery
================

Delve has been used in several papers:

.. figure:: gallery/images/border_ResNet18_Cifar10_32.png
  :width: 400
ResNet18 trained on Cifar10 for 30 epochs using the adam optimizer and a batch size of 64. Image from `"Should You Go Deeper? Optimizing Convolutional Neural Networks without training" <https://arxiv.org/abs/2106.12307>`_.

.. figure:: gallery/images/border_ResNet34_Cifar10_32.png
  :width: 400
ResNet34 trained on Cifar10 for 30 epochs using the adam optimizer. Image from `"Should You Go Deeper? Optimizing Convolutional Neural Networks without training" <https://arxiv.org/abs/2106.12307>`_.

.. figure:: gallery/images/DenseNet18_Cifar10_32.png
  :width: 400
DenseNet18 trained on Food101 for 90 epochs using the stochastic gradient decent optimizer and a batch size of 128. Image from `"Feature Space Saturation During Training" <https://arxiv.org/abs/2006.08679>`_.

.. figure:: gallery/images/vgg16_resolution_sat.png
  :width: 400
VGG16 trained on 3 different resolutions for 30 epochs using the Adam-optimizer and a batch size of 32. You can see the shift in the inference process by observing the shift in high saturation values. Image from `"(Input) Size Matters for Convolutional Neural Network Classifiers" <https://www.springerprofessional.de/en/input-size-matters-for-cnn-classifiers/19652392>`_.
Support for Delve
=================

Bugs
----

Bugs, issues and improvement requests can be logged in `Github Issues <https://github.com/delve-team/delve/issues>`_.

Community
---------

Community support is provided via `Gitter <https://gitter.im/delve-chat/community>`_. Just ask a question there.
delve
=====

.. toctree::
   :maxdepth: 4

   delve
delve package
=============

Submodules
----------

delve.metrics module
--------------------

.. automodule:: delve.metrics
   :members:
   :undoc-members:
   :show-inheritance:

delve.pca\_layers module
------------------------

.. automodule:: delve.pca_layers
   :members:
   :undoc-members:
   :show-inheritance:

delve.tools module
------------------

.. automodule:: delve.tools
   :members:
   :undoc-members:
   :show-inheritance:

delve.torch\_utils module
-------------------------

.. automodule:: delve.torch_utils
   :members:
   :undoc-members:
   :show-inheritance:

delve.torchcallback module
--------------------------

.. automodule:: delve.torchcallback
   :members:
   :undoc-members:
   :show-inheritance:

delve.writers module
--------------------

.. automodule:: delve.writers
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: delve
   :members:
   :undoc-members:
   :show-inheritance:
.. Delve documentation master file, created by
   sphinx-quickstart on Sun Aug 22 13:20:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: ../../README.rst

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   Installation <install>
   Saturation <saturation>
   Usage <gallery/index>
   Academic Gallery <gallery>
   Example Plots <examples>
   Reference to All Attributes and Methods <reference>
   Bugs and Support <support>
   Contributing to Delve <contributing>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Installation
============

Installing Delve
----------------

Delve require Python 3.6+ to be installed.

To install via pip:

.. code::

   pip install delve

To install the latest development version, clone the `GitHub` repository and use the setup script::

   git clone https://github.com/delve-team/delve.git
   cd delve
   pip install .

Usage
-----

Instantiate the :class:`~delve.torchcallback.SaturationTracker` class where you define your PyTorch training loop, as in the example::

   from torch import nn
   from delve import SaturationTracker

   ...

   model = nn.ModuleDict({
                 'conv1': nn.Conv2d(1, 8, 3, padding=1),
                 'linear1': nn.Linear(3, 1),
   })


   layers = [model.conv1, model.linear1]
   stats = SaturationTracker('regression/h{}'.format(h),
      save_to="plotcsv",
      modules=layers,
      stats=["lsat"]
   )

   ...

   for _ in range(10):
      y_pred = model(x)
      loss = loss_fn(y_pred, y)
      optimizer.zero_grad()
      loss.backward()
      optimizer.step()

      stats.add_saturations()

   stats.close()

This will hook into the layers in ``layers`` and log the statistics, in this case ``lsat`` (layer saturation). It will save images to ``regression``.

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "gallery/extract-saturation.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_gallery_extract-saturation.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_gallery_extract-saturation.py:


Extract layer saturation
------------------------
Extract layer saturation with Delve.

.. GENERATED FROM PYTHON SOURCE LINES 6-80

.. code-block:: python

    import torch
    from tqdm import trange

    from delve import SaturationTracker


    class TwoLayerNet(torch.nn.Module):
        def __init__(self, D_in, H, D_out):
            super(TwoLayerNet, self).__init__()
            self.linear1 = torch.nn.Linear(D_in, H)
            self.linear2 = torch.nn.Linear(H, D_out)

        def forward(self, x):
            h_relu = self.linear1(x).clamp(min=0)
            y_pred = self.linear2(h_relu)
            return y_pred


    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(1)

    for h in [3, 32]:
        # N is batch size; D_in is input dimension;
        # H is hidden dimension; D_out is output dimension.
        N, D_in, H, D_out = 64, 1000, h, 10

        # Create random Tensors to hold inputs and outputs
        x = torch.randn(N, D_in)
        y = torch.randn(N, D_out)
        x_test = torch.randn(N, D_in)
        y_test = torch.randn(N, D_out)

        # You can watch specific layers by handing them to delve as a list.
        # Also, you can hand over the entire Module-object to delve and let delve search for recordable layers.
        model = TwoLayerNet(D_in, H, D_out)

        x, y, model = x.to(device), y.to(device), model.to(device)
        x_test, y_test = x_test.to(device), y_test.to(device)

        layers = [model.linear1, model.linear2]
        stats = SaturationTracker('regression/h{}'.format(h),
                                  save_to="plotcsv",
                                  modules=layers,
                                  device=device,
                                  stats=["lsat", "lsat_eval"])

        loss_fn = torch.nn.MSELoss(reduction='sum')
        optimizer = torch.optim.SGD(model.parameters(), lr=1e-4, momentum=0.9)
        steps_iter = trange(2000, desc='steps', leave=True, position=0)
        steps_iter.write("{:^80}".format(
            "Regression - TwoLayerNet - Hidden layer size {}".format(h)))
        for step in steps_iter:
            # training step
            model.train()
            y_pred = model(x)
            loss = loss_fn(y_pred, y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # test step
            model.eval()
            y_pred = model(x_test)
            loss_test = loss_fn(y_pred, y_test)

            # update statistics
            steps_iter.set_description('loss=%g' % loss.item())
            stats.add_scalar("train-loss", loss.item())
            stats.add_scalar("test-loss", loss_test.item())

            stats.add_saturations()
        steps_iter.write('\n')
        stats.close()
        steps_iter.close()


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  0.000 seconds)


.. _sphx_glr_download_gallery_extract-saturation.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: extract-saturation.py <extract-saturation.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: extract-saturation.ipynb <extract-saturation.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
:orphan:



.. _sphx_glr_gallery:

Gallery
==================

A gallery of examples



.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Extract layer saturation">

.. only:: html

 .. figure:: /gallery/images/thumb/sphx_glr_extract-saturation_thumb.png
     :alt: Extract layer saturation

     :ref:`sphx_glr_gallery_extract-saturation.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /gallery/extract-saturation
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-gallery


  .. container:: sphx-glr-download sphx-glr-download-python

    :download:`Download all examples in Python source code: gallery_python.zip </gallery/gallery_python.zip>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

    :download:`Download all examples in Jupyter notebooks: gallery_jupyter.zip </gallery/gallery_jupyter.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
﻿delve.SaturationTracker
===================

.. currentmodule:: delve

.. autoclass:: SaturationTracker
﻿delve.TorchCovarianceMatrix
===========================

.. currentmodule:: delve

.. autoclass:: TorchCovarianceMatrix
