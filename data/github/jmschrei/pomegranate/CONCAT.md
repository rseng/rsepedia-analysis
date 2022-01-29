<img src="https://github.com/jmschrei/pomegranate/blob/master/docs/logo/pomegranate-logo.png" width=300>

[![Downloads](https://pepy.tech/badge/pomegranate)](https://pepy.tech/project/pomegranate)![build](https://github.com/jmschrei/pomegranate/workflows/build/badge.svg) [![Documentation Status](https://readthedocs.org/projects/pomegranate/badge/?version=latest)](http://pomegranate.readthedocs.io/en/latest/?badge=latest) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jmschrei/pomegranate/master)

Please consider citing the [**JMLR-MLOSS Manuscript**](http://jmlr.org/papers/volume18/17-636/17-636.pdf) if you've used pomegranate in your academic work!

pomegranate is a package for building probabilistic models in Python that is implemented in Cython for speed. A primary focus of pomegranate is to merge the easy-to-use API of scikit-learn with the modularity of probabilistic modeling to allow users to specify complicated models without needing to worry about implementation details. The models implemented here are built from the ground up with big data processing in mind and so natively support features like multi-threaded parallelism and out-of-core processing. Click on the binder badge above to interactively play with the tutorials!

### Installation

pomegranate is pip-installable using `pip install pomegranate` and conda-installable using `conda install pomegranate`. If neither work, more detailed installation instructions can be found [here](http://pomegranate.readthedocs.io/en/latest/install.html).

If you get an error involving `pomegranate/base.c`, try installing with `pip install --no-cache-dir pomegranate`.

If you get an error involving `pomegranate/distributions/NeuralNetworkWrapper.c: No such file or directory`, try installing Cython first and then re-installing. 

A few packages are optional to use pomegranate but necessary for some specific functionality. For example, pandas is needed to run the tests involving I/O, matplotlib and pygraphviz are needed for plotting capabilities, and cupy is needed for GPU acceleration. 

### Models

* [Probability Distributions](http://pomegranate.readthedocs.io/en/latest/Distributions.html)
* [General Mixture Models](http://pomegranate.readthedocs.io/en/latest/GeneralMixtureModel.html)
* [Hidden Markov Models](http://pomegranate.readthedocs.io/en/latest/HiddenMarkovModel.html)
* [Naive Bayes and Bayes Classifiers](http://pomegranate.readthedocs.io/en/latest/NaiveBayes.html)
* [Markov Chains](http://pomegranate.readthedocs.io/en/latest/MarkovChain.html)
* [Discrete Bayesian Networks](http://pomegranate.readthedocs.io/en/latest/BayesianNetwork.html)
* [Discrete Markov Networks](https://pomegranate.readthedocs.io/en/latest/MarkovNetwork.html)

The discrete Bayesian networks also support novel work on structure learning in the presence of constraints through a constraint graph. These constraints can dramatically speed up structure learning through the use of loose general prior knowledge, and can frequently make the exact learning task take only polynomial time instead of exponential time. See the [PeerJ manuscript](https://peerj.com/articles/cs-122/) for the theory and the [pomegranate tutorial](https://github.com/jmschrei/pomegranate/blob/master/tutorials/B_Model_Tutorial_4b_Bayesian_Network_Structure_Learning.ipynb) for the practical usage! 

To support the above algorithms, it has efficient implementations of the following:

* Kmeans/Kmeans++/Kmeans||
* Factor Graphs

### Features

* [sklearn-like API](https://pomegranate.readthedocs.io/en/latest/api.html)
* [Multi-threaded Training](http://pomegranate.readthedocs.io/en/latest/parallelism.html)
* [BLAS/GPU Acceleration](http://pomegranate.readthedocs.io/en/latest/gpu.html)
* [Out-of-Core Learning](http://pomegranate.readthedocs.io/en/latest/ooc.html)
* [Data Generators and IO](https://pomegranate.readthedocs.io/en/latest/io.html)
* [Semi-supervised Learning](http://pomegranate.readthedocs.io/en/latest/semisupervised.html)
* [Missing Value Support](http://pomegranate.readthedocs.io/en/latest/nan.html)
* [Customized Callbacks](http://pomegranate.readthedocs.io/en/latest/callbacks.html)

Please take a look at the [tutorials folder](https://github.com/jmschrei/pomegranate/tree/master/tutorials), which includes several tutorials on how to effectively use pomegranate!

See [the website](http://pomegranate.readthedocs.org/en/latest/) for extensive documentation, API references, and FAQs about each of the models and supported features.

No good project is done alone, and so I'd like to thank all the previous contributors to YAHMM, and all the current contributors to pomegranate, including the graduate students who share my office I annoy on a regular basis by bouncing ideas off of.

### Dependencies

pomegranate requires:

```
- Cython (only if building from source)
- NumPy
- SciPy
- NetworkX
- joblib
```

To run the tests, you also must have `nose` installed.

## Contributing

If you would like to contribute a feature then fork the master branch (fork the release if you are fixing a bug). Be sure to run the tests before changing any code. You'll need to have [nosetests](https://github.com/nose-devs/nose) installed. The following command will run all the tests:

```
python setup.py test
```

Let us know what you want to do just in case we're already working on an implementation of something similar. This way we can avoid any needless duplication of effort. Also, please don't forget to add tests for any new functions.

---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is, including what you were expecting to happen and what actually happened. Please report the version of pomegranate that you are using and the operating system. Also, please make sure that you have upgraded to the latest version of pomegranate before submitting the bug report.

**To Reproduce**
Please provide a snippet of code that can reproduce this error. It is much easier for us to track down bugs and fix them if we have an example script that fails until we're successful.
# Tutorial

This tutorial covers both the concept of probabilistic modelling, and how to do such using pomegranate. The tutorials show probabilistic models of increasing complexity, with a single distribution being the simplest, and a Bayesian network being the most sophisticated. 

To read a tutorial, just click on the notebook above. Github can render IPython notebooks natively, so you don't need to go through the complicated procedure of rendering them somewhere else.
# Tutorial

This tutorial covers both the concept of probabilistic modelling, and how to do such using pomegranate. The tutorials show probabilistic models of increasing complexity, with a single distribution being the simplest, and a Bayesian network being the most sophisticated. 

To read a tutorial, just click on the notebook above. Github can render IPython notebooks natively, so you don't need to go through the complicated procedure of rendering them somewhere else.
