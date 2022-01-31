# Change Log

## v3.1.0
- Added support for training with Datasets, generators, and other data types supported by Keras [#211](https://github.com/NLeSC/mcfly/issues/211)
- Added separate classes for each model type, also allowing easier extension by own custom models [#239](https://github.com/NLeSC/mcfly/pull/239)
- Dropped support for Python 3.5 (supported are: Python 3.6, 3.7 and 3.8)
- Extend and update documentation in readthedocs [#250](https://github.com/NLeSC/mcfly/pull/250)
- Fix broken model visualization [#256](https://github.com/NLeSC/mcfly/pull/256)

## v3.0.0
- Add ResNet architecture
- Add InceptionTime architecture
- Tensorflow dependency to 2.0
- Dropped support for Python 2.7 and added support for Python 3.7
- Early_stopping argument for train_models_on_samples() changed to early_stopping_patience
- Fix metric name issue in visualization
- Lower level functions arguments have changed by using keyword arguments dic

## v2.1.0
- Add class weight support

## v2.0.1
- Fix documentation inconsistency

## v2.0.0
- Using Tensorflow.keras instead of Keras
- Using keyword 'accuracy' instead of 'acc' in logs like latest keras versions do

## v1.0.5
- Requirements change (keras<2.3.0)

## v1.0.4
- Fixed Zenodo configuration

## v1.0.3
- Requirements change (six>=1.10.0)
- Fixed Zenodo configuration

## v1.0.2
- Small bug fixes
- Added metric option to model training
- Extended documentation
- Removed redundant dependency on Pandas

## v1.0.1
- Small bug fixes
- Compatible with Keras v2.x (no longer compatible with Keras < v2.0.0)
- Tutorial is moved to separate repository (https://github.com/NLeSC/mcfly-tutorial)

## v1.0.0
First major release.
<p align="left">
  <img src="https://raw.githubusercontent.com/NLeSC/mcfly/master/mcflylogo.png" width="200"/>
</p>

![GitHub Workflow Status](https://img.shields.io/github/workflow/status/NLeSC/mcfly/CI%20Build)
[![Coverage](https://scrutinizer-ci.com/g/NLeSC/mcfly/badges/coverage.png?b=master)](https://scrutinizer-ci.com/g/NLeSC/mcfly/statistics/)
[![PyPI](https://img.shields.io/pypi/v/mcfly.svg)](https://pypi.python.org/pypi/mcfly/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596127.svg)](https://doi.org/10.5281/zenodo.596127)
[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/nlesc/mcfly)
<!-- The first 12 lines are skipped while generating 'long description' (see setup.py)) -->

The goal of mcfly is to ease the use of deep learning technology for time series classification. The advantage of deep learning is that it can handle raw data directly, without the need to compute signal features. Deep learning does not require  expert domain knowledge about the data, and has been shown to be competitive with conventional machine learning techniques. As an example, you can apply mcfly on accelerometer data for activity classification, as shown in [the tutorial](https://github.com/NLeSC/mcfly-tutorial).

If you use mcfly in your research, please cite the following software paper:

D. van Kuppevelt, C. Meijer, F. Huber, A. van der Ploeg, S. Georgievska, V.T. van Hees. _Mcfly: Automated deep learning on time series._
SoftwareX,
Volume 12,
2020.
[doi: 10.1016/j.softx.2020.100548](https://doi.org/10.1016/j.softx.2020.100548)

## Installation
Prerequisites:
- Python 3.6, 3.7 or 3.8
- pip
- Tensorflow 2.0, if pip errors that it can't find it for your python/pip version

Installing all dependencies in separate conda environment:
```sh
conda env create -f environment.yml

# activate this new environment
source activate mcfly
```

To install the package, run in the project directory:

`pip install .`

### Installing on Windows
When installing on Windows, there are a few things to take into consideration. The preferred (in other words: easiest) way to install Keras and mcfly is as follows:
* Use [Anaconda](https://www.anaconda.com/download)
* Install numpy and scipy through the conda package manager (and not with pip)
* To install mcfly, run `pip install mcfly` in the cmd prompt.
* Loading and saving models can give problems on Windows, see https://github.com/NLeSC/mcfly-tutorial/issues/17

## Visualization
We build a tool to visualize the configuration and performance of the models. The tool can be found on http://nlesc.github.io/mcfly/. To run the  model visualization on your own computer, cd to the `html` directory and start up a python web server:

`python -m http.server 8888 &`

Navigate to `http://localhost:8888/` in your browser to open the visualization. For a more elaborate description of the visualization see [user manual](https://mcfly.readthedocs.io/en/latest/user_manual.html).


## User documentation
[User and code documentation](https://mcfly.readthedocs.io).

## Contributing
You are welcome to contribute to the code via pull requests. Please have a look at the [NLeSC guide](https://nlesc.gitbooks.io/guide/content/software/software_overview.html) for guidelines about software development.

We use numpy-style docstrings for code documentation.

#### Necessary steps for making a new release
* Check citation.cff using general DOI for all version (option: create file via 'cffinit')
* Create .zenodo.json file from CITATION.cff (using cffconvert)  
```cffconvert --validate```  
```cffconvert --ignore-suspect-keys --outputformat zenodo --outfile .zenodo.json```
* Set new version number in mcfly/_version.py
* Edit Changelog (based on commits in https://github.com/NLeSC/mcfly/compare/v1.0.1...master)
* Create Github release
* Upload to pypi:  
```python setup.py sdist bdist_wheel```  
```python -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/*```  
(or ```python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*``` to test first)
* Check doi on zenodo
* If the visualization has changed, deploy it to github pages:
```
git subtree push --prefix html origin gh-pages
```

## Licensing
Source code and data of mcfly are licensed under the Apache License, version 2.0.
User Manual
===========

On this page, we describe what you should know when you use mcfly. This manual should be understandable without too much knowledge of deep learning,
although it expects familiarity with the concepts of dense hidden layers, convolutional layers and recurrent layers.
However, if mcfly doesn't give you a satisfactory model, a deeper knowledge of deep learning really helps in debugging the models.

We provide a quick description for the layers used in mcfly.

* **dense layer** also know as fully connected layer, is a layer of nodes that all have connections to all outputs of the previous layer.
* **convolutional layer** convolves the output of the previous layer with one or more sets of weights and outputs one or more feature maps.
* **LSTM layer** is a recurrent layer with some special features to help store information over multiple time steps in time series.

Some recommended reading to make you familiar with deep learning:
http://scarlet.stanford.edu/teach/index.php/An_Introduction_to_Convolutional_Neural_Networks

Or follow a complete course on deep learning:
http://cs231n.stanford.edu/


Data preprocessing
-------------------

The input for the mcfly functions is data that is already preprocessed to be handled for deep learning module Keras.
In this section we describe what data format is expected and what to think about when preprocessing.

Eligible data sets
^^^^^^^^^^^^^^^^^^
Mcfly is a tool for *classification* of single or *multichannel timeseries data*. One (real valued) multi-channel time series is associated with one class label.
All sequences in the data set should be of equal length.

Data format
^^^^^^^^^^^
The data should be split in train, validation and test set. For each of the splits, the input X and the output y are both numpy arrays.

The input data X should be of shape (num_samples, num_timesteps, num_channels). The output data y is of shape (num_samples, num_classes), as a binary array for each sample.

We recommend storing the numpy arrays as binary files with the numpy function ``np.save``.

Data preprocessing
^^^^^^^^^^^^^^^^^^
Here are some tips for preprocessing the data:

* For longer, multi-label sequences, we recommend creating subsequences with a sliding window. The length and step of the window can be based on domain knowledge.
* One label should be associated with a complete time series. In case of multiple labels, often the last label is taken as the label for the complete sequence.
  Another possibility is to take the majority label.
* In splitting the data into training, validation and test sets, it might be necessary to make sure that sample subject (such as test persons) for which multiple sequences are available, are not present in both train and validation/test set. The same holds for subsequences that originates from the same original sequence (in case of sliding windows).
* The Keras function ``keras.utils.np_utils.to_categorical`` can be used to transform an array of class labels to binary class labels.
* Data doesn't need to be normalized. Every model mcfly produces starts by normalizing data through a Batch Normalization layer.
  This means that training data is used to learn mean and standard deviation of each channel and timestep.

Finding the best architecture
---------------------------------
The function :func:`~mcfly.find_architecture.find_best_architecture` generates a variety of architectures and hyperparameters,
and returns the best performing model on a subset of the data.
The following four types of architectures are possible (for more information, see the :doc:`technical_doc`):

:class:`~mcfly.models.CNN`: A stack ofonvolutional layers, followed by a final dense layer

:class:`~mcfly.models.ConvLSTM`: Convolutional layers, followed by LSTM layers and a final dense layer

:class:`~mcfly.models.ResNet`: Convolutional layers with skip connections

:class:`~mcfly.models.InceptionTime`: Convolutional layers ('inception module') with different kernel sizes in parallel, concatenated and then followed by pooling and a dense layer.

The hyperparameters to be optimized are the following:

* learning rate
* regularization rate
* model_type: *CNN* or *DeepConvLSTM*
* if modeltype=CNN:
   * number of Conv layers
   * for each Conv layer: number of filters
   * number of hidden nodes for the hidden Dense layer

* if modeltype=DeepConvLSTM:
   * number of Conv layers
   * for each Conv layer: number of filters
   * number of LSTM layers
   * for each LSTM layer: number of hidden nodes

   * if modeltype=ResNet:
      * network depth, i.e. number of residual modules
      * minimum number of filters
      * maximum kernel size

   * if modeltype=InceptionTime:
      * number of filters for all convolutional layers
      * depth of network, i.e. number of Inception modules to stack.
      * maximum kernel size


We designed mcfly to have sensible default values and ranges for each setting.
However, you have the possibility to influence the behavior of the function with the arguments that you give to it to try other values.
See the the documentation of :func:`~mcfly.modelgen.generate_models` for all options, among others:
* **number_of_models**: the number of models that should be generated and tested
* **nr_epochs**: The models are tested after only a small number of epochs, to limit the time. Setting this number higher will give a better estimate of the performance of the model, but it will take longer
* **model_types**: List of all model architecture types to choose from
* Ranges for all of the hyperparameters: The hyperparameters (as described above) are sampled from a uniform or log-uniform distribution. The boundaries of these distributions have default values (see the arguments :func:`~mcfly.modelgen.generate_models`), but can be set custom.



Visualize the training process
-------------------------------
To gain more insight in the training process of the models and the influence of the hyperparameters, you can explore the visualization.

1. Save the model results, by defining `outputpath` in `find_best_architecture`.

2. Start an python webserver (see :doc:`installation`) and navigate to the visualization page in your browser.

3. Open the json file generated in step 1.

In this visualization, the accuracy on the train and validation sets are plotted for all models. You can filter the graphs by selecting specific models, or filter on hyperparameter values.

FAQ
---

None of the models that are tested in findBestArchitecture perform satisfactory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Note that :func:`~mcfly.find_architecture.find_best_architecture` doesn't give you a fully trained model yet: it still needs to be trained on the complete dataset with sufficient iterations.
However, if none of the models in :func:`~mcfly.find_architecture.find_best_architecture` have a better accuracy than a random model, it might be worth trying one of the following things:

* Train more models: the number of models tested needs to be sufficient to cover a large enough part of the hyperparameter space
* More epochs: it could be that the model needs more epochs to learn (for example when the learning rate is small). Sometimes this is visible from the learning curve plot
* Larger subset size: it could be that the subset of the train data is too small to contain enough information for learning
* Extend hyperparameter range
mcfly
=====
.. automodule:: mcfly
    :members:
    :undoc-members:
    :show-inheritance:

mcfly.models
------------------------------

.. automodule:: mcfly.models
    :members:
    :undoc-members:
    :show-inheritance:

mcfly.find_architecture module
------------------------------

.. automodule:: mcfly.find_architecture
    :members:
    :undoc-members:
    :show-inheritance:

mcfly.modelgen module
---------------------

.. automodule:: mcfly.modelgen
    :members:
    :undoc-members:
    :show-inheritance:

mcfly.storage module
--------------------

.. automodule:: mcfly.storage
    :members:
    :undoc-members:
    :show-inheritance:
Technical documentation
=======================

This page describes the technical implementation of mcfly and the choices that have been made.

Hyperparameter search
---------------------
The function :func:`~mcfly.modelgen.generate_models` generates a list of untrained Keras models based on the shape of the input data.
Key function arguments include *number_of_models to specify the number of models that need to be generated, *x_shape* to specify the shape of the input training data,
and number_of_classes to specify the number of classes in the data labels. An additional set of arguments, with default values, provide more control for advanced users.

The model architectures generated in this previous step are compared with the function :func:`~mcfly.find_architecture.train_models_on_samples`,
by training them the training data (or possibly a subset thereof) and evaluating the models on the validation subset.
By looking for the best performing model mcfly effectively performs a random search over the hyper parameter space.
We chose to implement random search, because it's simple and fairly effective (`Bergstra 2012 <https://dl.acm.org/doi/abs/10.5555/2188385.2188395>`_). This will help us to choose the best candidate model.
The performance results for the models are stored in a json file, which we will visually inspect later on.
There are a few optional settings for training the models in mcfly, such as the number of epochs and whether the training process should stop when performance does not increase (*early_stopping*).


Architectures
-------------
There are currently four types of architectures available in mcfly: CNN,  DeepConvLSTM, ResNet and InceptionTime.
The details of these architectures are discussed in the next subsections, and the figure below shows an overview of all architectures.
The  argument *model_types* gives the user the opportunity to specify which architectures to choose from, by default all four model types are used.
The model generation function chooses the model type at random, but in an alternating order so that the number of models per model type is roughly equal.
The first layer in all architectures is a Batchnorm layer, so that the user does not have to normalize the input data.

.. figure:: network_architectures.png
    :width: 600px
    :align: center
    :alt: alternate text
    :figclass: align-center

    Architectures of the convolutional, LSTM, ResNet, and InceptionTime networks generated by mcfly

CNN
^^^
The model type CNN is a regular Convolutional Neural Network, with N convolutional layers with ReLU activation  (`Nair2010rectified <https://dl.acm.org/doi/10.5555/3104322.3104425>`_) and one hidden dense layer.
An example of the application of CNNs on time series data is in (`Martinez 2013 <https://ieeexplore.ieee.org/document/6496209>`_).

The hyperparameters of the CNN are the number of Convolutional layers (default range to choice from: [1,10]), as well as the number of filters (default range: [10,100]) in each Convolutional layer and the number of neurons in the hidden Dense layer (default range [10,2000]). We decided not to add Pool layers because reducing the spatial size of the sequence is usually not necessary if there are enough convolutional layers.

DeepConvLSTM
^^^^^^^^^^^^
The architecture of the model type DeepConvLSTM is based on the paper by Ordonez et al. (`Ordonez 2016 <http://www.mdpi.com/1424-8220/16/1/115>`_), \originally developed for multimodal wearable
activity recognition. The model combines Convolutional layers with Long short-term memory (LTSM) layers.
In contrast to the CNN model, the convolutional layers in the DeepConvLSTM model are applied per channel, and connected in the first LSTM layer.
Applying convolutions per channel was reported to be suitable for activity recognition tasks based on wearable sensor data where the channels are not necessarily synchronized in time
(each channel corresponds to a different body part) (`Ordonez 2016 <http://www.mdpi.com/1424-8220/16/1/115>`_).
The LSTM layers result in a multidimensional vector per time step.
The Timedistributed Dense layer is used to output a sequence of predictions, and the TakeLast layer to pick the last element from the sequence as the final prediction.

The hyperparameters of the DeepConvLSTM architecture are the number of convolutional layers (default range to choose from: [1,10]), the number of LSTM layers (default range: [1,5]),
the number of filters for each Conv layer (default range: [10,100]) and the hidden layer dimension for each LSTM layer (default range: [10,100]).

ResNet
^^^^^^^^^^^^
Many convolutional architecture are inspired by the field of image classification where deep learning had its first and biggest breakthroughs.
One very successful technique for dealing efficiently with deeper networks are so called skip-connections (also termed residual connections)
which are connections that skip multiple layers within deeper neural networks (`He 2015 <https://ieeexplore.ieee.org/document/7780459/>`_).
This neural network architecture has also been adapted for time series classification (`Wang 2016 <https://ieeexplore.ieee.org/document/7966039>`_)
and was shown to perform comparably well in a large number of cases (`Fawaz 2019(1) <https://doi.org/10.1007/s10618-019-00619-1>`_).

The hyperparameters of the ResNet architecture are the number of residual modules (default range to choose from: [2,5]).
One such module consists of three convolutional layers and one skip connection bridging all three layers.
Two other key hyperparameters are the number of filters and the kernel size for the convolutions.
Following the original ResNet design for time series classification ((`Wang 2016 <https://ieeexplore.ieee.org/document/7966039>`_),
we chose a maximum kernel size (default range: [8,32]) and derive the kernel sizes for the different levels by scaling down by :math:`2^{-i/2}` for i the index of the residual module. Analogously, a minimum number of filters is chosen (default range: [32,128])
based on which the numbers of filters for all residual modules is derived following :math:`2^{i/2}` with i the index of the residual module.

InceptionTime
^^^^^^^^^^^^^^
Another architecture element that turned out to be very helpful in handling neural networks of greater depth and width are inception modules (`Szegedy 2014 <https://ieeexplore.ieee.org/document/7298594>`_).
Inception modules run convolutions with different kernel sizes in parallel and then combine the outcome.
While initially applied to image classification problems, their adaptation for time series classification was recently seen to deliver very promising results on a wide variety of data sets (`Fawaz 2019(2) <https://arxiv.org/abs/1909.04939>`_).

As key hyper-parameters based on the InceptionTime architecture we picked the number of Inception modules (default range to choose from: [3,6]).
Each such module consists of a bottleneck layer (max pooling) followed by three convolutional layers with varying kernel sizes and one convolution with kernel size=1.
Mcfly randomly chooses a maximum kernel size (default range: [10,100]) based on which the kernel sizes for the different layers are derived by scaling down by dividing by 2 and 4.
Another key parameters is the number of filters which is the same for all layers in is chosen (default range: [32,96]).

Other choices
-------------
We have made the following choices for all models:

* We use LeCun Uniform weight initialization (`LeCun 1998 <http://yann.lecun.com/exdb/publis/pdf/lecun-98b.pdf>`_)
* The kernel size is 3
* We use L2 regularization on all convolutional and dense layers (`Ng 2004 <https://dl.acm.org/doi/10.1145/1015330.1015435>`_)
* We use categorical cross-entropy loss (`Mannor 2005 <http://portal.acm.org/citation.cfm?doid=1102351.1102422>`_)
* We output accuracy and take this as a measure to choose the best performing model. Note that the performance metric can be changed with argument *metric* in most mcfly functions.
* The default, but modifiable, log range for the learning rate and the regularization rate is [:math:`10^{-1}, 10^{-4}`] .


Comparison with non-deep models
---------------------------------
To check the value of the data, a 1-Nearest Neighbors model is applied as a benchmark for the deep learning model.
We chose 1-NN because it's a very simple, hyperparameter-free model that often works quite well on time series data.
For large train sets, 1-NN can be quite slow: the test-time performance scales linear with the size of the training set.
However, we perform the check only on a small subset of the training data.
The related Dynamic Time Warping (DTW) algorithm has a better track record for classifying time series,
but we decided not to use it because it's too slow (it scales quadratically with the length of the time series).
Installation
============

Prerequisites:

* Python 3.6, 3.7 or 3.8
* pip
* Tensorflow 2.0, if pip errors that it can’t find it for your python/pip version


Installing all dependencies in separate conda environment:

.. code:: sh

   conda env create -f environment.yml

   # activate this new environment
   source activate mcfly

To install the package, run in the project directory:

``pip install .``

Installing on Windows
~~~~~~~~~~~~~~~~~~~~~

When installing on Windows, there are a few things to take into
consideration. The preferred (in other words: easiest) way to install
Keras and mcfly is as follows:

* Use `Anaconda <https://www.anaconda.com/download>`__
* Install numpy and scipy through the conda package manager (and not with pip)
* To install mcfly, run ``pip install mcfly`` in the cmd prompt.
* Loading and saving models can give problems on Windows, see https://github.com/NLeSC/mcfly-tutorial/issues/17


Visualization
~~~~~~~~~~~~~

We build a tool to visualize the configuration and performance of the
models. The tool can be found on http://nlesc.github.io/mcfly/. To run
the model visualization on your own computer, cd to the ``html``
directory and start up a python web server:

``python -m http.server 8888 &``

Navigate to ``http://localhost:8888/`` in your browser to open the
visualization. For a more elaborate description of the visualization see
`user
manual <https://mcfly.readthedocs.io/en/latest/user_manual.html>`__.
.. mcfly documentation master file, created by
   sphinx-quickstart on Wed Oct  5 10:42:14 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to mcfly's documentation!
=================================

The goal of mcfly is to ease using deep learning technology for time series classification.
The advantages of deep learning algorithms is that it can handle raw data directly with no need
to compute signal features, it does not require a expert domain knowledge about the data,
and it has been shown to be competitive with conventional machine learning techniques.
As an example, you can apply mcfly on, accelerometer data for activity classification,
as shown in the mcfly `tutorial <https://github.com/NLeSC/mcfly-tutorial>`_.

If you use mcfly in your reserach, please cite the following software paper:

D. van Kuppevelt, C. Meijer, F. Huber, A. van der Ploeg, S. Georgievska, V.T. van Hees. Mcfly: Automated deep learning on time series.
SoftwareX,
Volume 12,
2020.
`doi: 10.1016/j.softx.2020.100548 <https://doi.org/10.1016/j.softx.2020.100548>`_

Contents:

.. toctree::
   :maxdepth: 2

   introduction
   installation
   user_manual
   technical_doc


Reference
==================
.. toctree::
  :maxdepth: 2

  reference

* :ref:`search`
Introduction
============

The goal of mcfly is to ease the use of deep learning technology for
time series classification. The advantage of deep learning is that it
can handle raw data directly, without the need to compute signal
features. Deep learning does not require expert domain knowledge about
the data, and has been shown to be competitive with conventional machine
learning techniques. As an example, you can apply mcfly on accelerometer
data for activity classification, as shown in `the
tutorial <https://github.com/NLeSC/mcfly-tutorial>`__.

If you use mcfly in your research, please cite the following software
paper:

D. van Kuppevelt, C. Meijer, F. Huber, A. van der Ploeg, S. Georgievska,
V.T. van Hees. *Mcfly: Automated deep learning on time series.*
SoftwareX, Volume 12, 2020. `doi:
10.1016/j.softx.2020.100548 <https://doi.org/10.1016/j.softx.2020.100548>`__


Licensing
---------

Source code and data of mcfly are licensed under the Apache License,
version 2.0.

.. |Build Status| image:: https://travis-ci.org/NLeSC/mcfly.svg?branch=master
   :target: https://travis-ci.org/NLeSC/mcfly
.. |AppVeyor Build Status| image:: https://ci.appveyor.com/api/projects/status/lv8hih1hvxbuu5f7/branch/master?svg=true
   :target: https://ci.appveyor.com/project/NLeSC/mcfly/
.. |Coverage| image:: https://scrutinizer-ci.com/g/NLeSC/mcfly/badges/coverage.png?b=master
   :target: https://scrutinizer-ci.com/g/NLeSC/mcfly/statistics/
.. |PyPI| image:: https://img.shields.io/pypi/v/mcfly.svg
   :target: https://pypi.python.org/pypi/mcfly/
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.596127.svg
   :target: https://doi.org/10.5281/zenodo.596127
.. |Binder| image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org:/repo/nlesc/mcfly
