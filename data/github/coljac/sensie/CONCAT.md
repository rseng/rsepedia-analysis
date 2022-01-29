# Sensie v. 1.0.0

Sensie is toolset for probing the sensitivity of a deep neural 
network model to various properties of the data. 

## Overview

An active area of research (and not a little hand-wringing) in deep learning at present is better understanding how to interpret the key features learned by DNNs, particularly in order to better understand and predict failure modes. Various algorithms and toolsets exist for interpreting DNNs; in the computer vision arena, saliency maps are an important and well-motivated technique for understanding the key features of the images used by the network in order to make its decisions. However, these maps are not equally informative in all applications, for instance in some scientific applications where inputs are not RGB images, or physical parameters are otherwise not readily discerned from a saliency map.

Sensie probes the sensitivity of a network to inputs with a particular property, <img src="/tex/0d19b0a4827a28ecffa01dfedf5f5f2c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92146679999999pt height=14.15524440000002pt/>. This property can be a feature of the data; an otherwise known property of a test or training set that is not provided explicitly in training; or a function that can artificially vary this property for a supplied test set. The effect of a particular property is measured according to the variance it introduces to the correct output of the network, such as the score for the correct class <img src="/tex/e92bd1a512cbe7d9721e846312bfe3fc.svg?invert_in_darkmode&sanitize=true" align=middle width=14.523852749999989pt height=19.415200200000008pt/>. Quantitatively, we seek the output <img src="/tex/f6818419a9b2a0402d0b9bc468cc9189.svg?invert_in_darkmode&sanitize=true" align=middle width=20.966980649999996pt height=32.29212359999999pt/> for a supplied property, <img src="/tex/0d19b0a4827a28ecffa01dfedf5f5f2c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92146679999999pt height=14.15524440000002pt/>.

Although the algorithm is simple, it has the potential to highlight and quantify gaps in the training set or areas of increasing unrealiability in a trained model. Sensie was written with scientific applications in mind, but can be used on any trained model.

Sensie can operate on models trained for discrete classes (classifier); similar tools for regression will be added in the next release.

Sensie can be adapted to works on any predictive model.

**Main outstanding issues:**
- Proper support for regression models

## Requirements

Sensie assumes a model object with a ``predict`` method (as per Keras), but with a user supplied predictor function any framework should be usable. Python prequisites are listed in `requirements.txt`. The [pymc3](https://github.com/pymc-devs/pymc3) probabilistic programming framework is required for the (optional) determination of credible intervals for feature sensitivities. This can be used to assess the significance of a sensitivity analysis.

Sensie requires python 3.6 or above.

## Installation

After cloning the repository, I recommend creating a [virtual environment](https://docs.python.org/3/library/venv.html) for Sensie (currently, Sensie requires PyMC3 3.7). Clone the repository and install with:

```
pip install -e .
```

Which will install Sensie and the required python libraries. Optionally, install `pytest` with `pip install pytest`, then run the tests with `pytest test` from the repository root.

## Outline of usage

See the example below, and the notebooks in `examples/` for some detailed examples of how Sensie can be used. In summary:

- Create a `sensie.Probe` instance to wrap a pre-trained model:
```
    probe = sensie.Probe(model=model)
```
- Pass the Probe a test set, ground truths, and either a vector/tensor containing the property/properties to test, or a function that mutates a training example, taking a scalar parameter indicating the size of the effect.
```
    probe.predict_and_measure_perturbed(X_test, y_test, peturb_func, ...)
```
- Sensie will return a `sensie.SensitivityMeasure` object, a collection of the results of each test represented by `sensie.SingleTest` objects.
- Examine a plot of the outcome of each test, or print the `summary` to quantify the size of the effect.
- Is the effect significant? The sensitivity, i.e. the gradient of mean correct score with property <img src="/tex/0d19b0a4827a28ecffa01dfedf5f5f2c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92146679999999pt height=14.15524440000002pt/>, is determined using ordinary linear regression (`sensie.Probe.get_gradient`) or using Bayesian inference with PyMC3; this in turn supplies a 50%/95% credible interval for the sensitivity.
- Where the relationship is non-linear, polynomial fits can also be visualised in order to identify regions of the parameter space where the network is most sensitive to the supplied property.

Sensie has several options for plotting the results. See the [docs](https://sensie.readthedocs.io/en/latest/) and examples for more information.

Sensie assumes that `model.predict` returns a tensor `t` with categorical scores, such that `t[i, j]` returns the score for class *j* for test example *i*. If this is not the case, supply a predictor function at instantiation time that does: `probe = sensie.Probe(model, predictor)` where `predictor` is a function with the signature `predictor(model, x_test)` and returns a tensor of dimensions `(N, n_classes)`.

Detailed documentation can be accessed at [readthedocs.org](https://sensie.readthedocs.io/en/latest/).

## Calculation of the significance of the effect

When Sensie reports the significance of the effect, it is reporting whether there is a detectable
linear relationship between the property measured and model's score for the correct class for
each example. In other words, if we assume the accuracy is a function of some property of the inputs,
is the function affecting the scores in a significant way?

Formally, Sensie tests the assumption that 
<p align="center"><img src="/tex/77b8b439ca87e934bf6f2be3a45670ec.svg?invert_in_darkmode&sanitize=true" align=middle width=73.01983424999999pt height=16.438356pt/></p>

where 

<p align="center"><img src="/tex/9e4d84b983988415d2c3c7fefeda5dd7.svg?invert_in_darkmode&sanitize=true" align=middle width=101.81349254999999pt height=12.785402849999999pt/></p> 

and uses Bayesian linear regression to determine the gradient <img src="/tex/8e830a5ab471143f1bb80e525c09bbaa.svg?invert_in_darkmode&sanitize=true" align=middle width=15.24170009999999pt height=14.15524440000002pt/> along with
credible intervals. If the effect of property <img src="/tex/0d19b0a4827a28ecffa01dfedf5f5f2c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92146679999999pt height=14.15524440000002pt/> is small or neglible, then 
the result should be consistent with <img src="/tex/01ab42025bd0796cb94e61b4057f1e67.svg?invert_in_darkmode&sanitize=true" align=middle width=46.20045374999999pt height=21.18721440000001pt/>.

For properties without a natural ordering, like classes, Sensie can first 
re-order the classes/properties by the mean scores then measure the significance 
of the effect, i.e. whether there is a significant gradient under this ordering. 
Whether this is meaningful is highly context dependent! [1]

## Example

How sensitive is a model trained on MNIST digits to the orientation of the digits?
```
def rotate(x, angle):
  # code that rotates the image by _angle_ degrees goes here, see examples/MNIST
  # x_rotated = ...
  return x_rotated

model = load_model() # a trained model

sensie_mnist = sensie.Probe(model)
sensie_mnist.predict_and_measure_perturbed(X_test, y_test, 
                                            rotate, p_min=0, p_max=180, steps=10, 
                                            label='rotation', plot=True)
```
![MNIST rotation sensitivity](examples/sensie1.png)

From this we can see that the trained model is sensitive to the orientation of an input; beyond about 20 degrees, accuracy suffers significantly.

For this and some more complex examples, see the Jupyter notebooks in the `examples` directory.

## Documentation

Module docs can be found at [sensie.readthedocs.org](https://sensie.readthedocs.io/en/latest/).

## Bugs, questions and contributions

If you use Sensie and run into any problems, please open an issue here.

Any questions, comments and suggestions can be sent via GitHub or to colin(@)coljac.net.

Contributions are welcome. Fork this repository to your own machine, make some changes, and push your work back up to the fork and open a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) so that I can review the changes.

### References

[1]: Consider `plt.scatter(np.arange(0, 100), np.sort(np.random.random(100)))`.
# Sensie v. 0.1

Sensie is toolset for probing the sensitivity of a deep neural 
network model to various properties of the data. 

## Overview

An active area of research (and not a little hand-wringing) in deep learning at present is better understanding how to interpret the key features learned by DNNs, particularly in order to better understand and predict failure modes. Various algorithms and toolsets exist for interpreting DNNs; in the computer vision arena, saliency maps are an important and well-motivated technique for understanding the key features of the images used by the network in order to make its decisions. However, these maps are not equally informative in all applications, for instance in some scientific applications where inputs are not RGB images, or physical parameters are otherwise not readily discerned from a saliency map.

Sensie probes the sensitivity of a network to inputs with a particular property, $p_i$. This property can be a feature of the data; an otherwise known property of a test or training set that is not provided explicitly in training; or a function that can artificially vary this property for a supplied test set. The effect of a particular property is measured according to the variance it introduces to the correct output of the network, such as the score for the correct class $\overline{y}_c$. Quantitatively, we seek the output $\frac{\partial{\overline{y}_c}}{\partial{p_i}}$ for a supplied property, $p_i$.

Although the algorithm is simple, it has the potential to highlight and quantify gaps in the training set or areas of increasing unrealiability in a trained model. Sensie was written with scientific applications in mind, but can be used on any trained model.

Sensie can operate on models trained for discrete classes (classifier); similar tools for regression will be added in the next release.

Sensie can be adapted to works on any predictive model.

**Main outstanding issues:**
- Proper support for regression models

## Requirements

Sensie assumes a model object with a ``predict`` method (as per Keras), but with a user supplied predictor function any framework should be usable. Python prequisites are listed in `requirements.txt`. The [pymc3](https://github.com/pymc-devs/pymc3) probabilistic programming framework is required for the (optional) determination of credible intervals for feature sensitivities. This can be used to assess the significance of a sensitivity analysis.

Sensie requires python 3.6 or above.

## Installation

After cloning the repository, I recommend creating a [virtual environment](https://docs.python.org/3/library/venv.html) for Sensie (currently, Sensie requires PyMC3 3.7). Clone the repository and install with:

```
pip install -e .
```

Which will install Sensie and the required python libraries. Optionally, install `pytest` with `pip install pytest`, then run the tests with `pytest test` from the repository root.

## Outline of usage

See the example below, and the notebooks in `examples/` for some detailed examples of how Sensie can be used. In summary:

- Create a `sensie.Probe` instance to wrap a pre-trained model:
```
    probe = sensie.Probe(model=model)
```
- Pass the Probe a test set, ground truths, and either a vector/tensor containing the property/properties to test, or a function that mutates a training example, taking a scalar parameter indicating the size of the effect.
```
    probe.predict_and_measure_perturbed(X_test, y_test, peturb_func, ...)
```
- Sensie will return a `sensie.SensitivityMeasure` object, a collection of the results of each test represented by `sensie.SingleTest` objects.
- Examine a plot of the outcome of each test, or print the `summary` to quantify the size of the effect.
- Is the effect significant? The sensitivity, i.e. the gradient of mean correct score with property $p_i$, is determined using ordinary linear regression (`sensie.Probe.get_gradient`) or using Bayesian inference with PyMC3; this in turn supplies a 50%/95% credible interval for the sensitivity.
- Where the relationship is non-linear, polynomial fits can also be visualised in order to identify regions of the parameter space where the network is most sensitive to the supplied property.

Sensie has several options for plotting the results. See the [docs](https://sensie.readthedocs.io/en/latest/) and examples for more information.

Sensie assumes that `model.predict` returns a tensor `t` with categorical scores, such that `t[i, j]` returns the score for class *j* for test example *i*. If this is not the case, supply a predictor function at instantiation time that does: `probe = sensie.Probe(model, predictor)` where `predictor` is a function with the signature `predictor(model, x_test)` and returns a tensor of dimensions `(N, n_classes)`.

Detailed documentation can be accessed at [readthedocs.org](https://sensie.readthedocs.io/en/latest/).

## Calculation of the significance of the effect

When Sensie reports the significance of the effect, it is reporting whether there is a detectable
linear relationship between the property measured and model's score for the correct class for
each example. In other words, if we assume the accuracy is a function of some property of the inputs,
is the function affecting the scores in a significant way?

Formally, Sensie tests the assumption that 
$$ \bar{y_c} = f(p_i)  $$

where 

$$  \bar{y_c} = a_1 p_i + a_2 $$ 

and uses Bayesian linear regression to determine the gradient $a_1$ along with
credible intervals. If the effect of property $p_i$ is small or neglible, then 
the result should be consistent with $a_1 = 0$.

For properties without a natural ordering, like classes, Sensie can first 
re-order the classes/properties by the mean scores then measure the significance 
of the effect, i.e. whether there is a significant gradient under this ordering. 
Whether this is meaningful is highly context dependent! [1]

## Example

How sensitive is a model trained on MNIST digits to the orientation of the digits?
```
def rotate(x, angle):
  # code that rotates the image by _angle_ degrees goes here, see examples/MNIST
  # x_rotated = ...
  return x_rotated

model = load_model() # a trained model

sensie_mnist = sensie.Probe(model)
sensie_mnist.predict_and_measure_perturbed(X_test, y_test, 
                                            rotate, p_min=0, p_max=180, steps=10, 
                                            label='rotation', plot=True)
```
![MNIST rotation sensitivity](examples/sensie1.png)

From this we can see that the trained model is sensitive to the orientation of an input; beyond about 20 degrees, accuracy suffers significantly.

For this and some more complex examples, see the Jupyter notebooks in the `examples` directory.

## Documentation

Module docs can be found at [sensie.readthedocs.org](https://sensie.readthedocs.io/en/latest/).

## Bugs, questions and contributions

If you use Sensie and run into any problems, please open an issue here.

Any questions, comments and suggestions can be sent via GitHub or to colin(@)coljac.net.

Contributions are welcome. Fork this repository to your own machine, make some changes, and push your work back up to the fork and open a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) so that I can review the changes.

### References

[1]: Consider `plt.scatter(np.arange(0, 100), np.sort(np.random.random(100)))`.
---
title: '``Sensie``: Probing the sensitivity of neural networks'

tags:
  - Python
  - machine learning
  - neural networks
authors:
  - name: Colin Jacobs
    orcid: 0000-0003-4239-4055
    affiliation: "1" 
affiliations:
 - name: Center for Astrophysics and Supercomputing, Swinburne University of Technology
   index: 1
date: 20 March 2020
bibliography: paper.bib

---

# Introduction 

Deep neural networks (DNNs) are finding increasing application across a wide variety of fields, including in industry and scientific research. Although DNNs are able to successfully tackle data problems that proved intractable to other methods, for instance in computer vision, they suffer from a lack of interpretability.  Some well-known methods for visualising and interpreting the outputs of DNNs include directly inspecting the learned features of a model (e.g. convolutional kernels and activation maps); or producing "saliency maps" which allow a human researcher to inspect which parts of an input (such as an image) played the most crucial role in the model reaching a particular determination. These algorithms include occluding parts of an image [@zeilerVisualizingUnderstandingConvolutional2014], guided backpropagation [@simonyanDeepConvolutionalNetworks2013], and more complex algorithms such as SmoothGrad [@smilkovSmoothGradRemovingNoise2017] and Layer-wise relevance propagation [@binderLayerwiseRelevancePropagation2016].

However, simply inspecting the most salient regions of an input image (or other input) may not always be sufficiently interpretable. Here we present ``Sensie``, a python package to quickly inspect and quantify the sensitivity of a trained model to user-specified properties of the input domain, or under the arbitrary transformation of a test set.

Several quality software packages exist to probe artificial neural network internals, such as ``Innvestigate`` [@alberINNvestigateNeuralNetworks2018], ``tf-explain``[^tfe], and ``keras-vis`` [@raghakotkerasvis]. These packages are geared towards traditional image-based applications, and cannot explore relationships between the neural network performance and arbitrary properties of the examples presented to the network. ``Sensie`` is designed to fill that gap.

# Method

``Sensie`` helps the user explore the sensitivity of a machine learning model to properties of data (using a test set with known, correct outputs). ``Sensie`` wraps a pre-trained model, and assists in quantifying and visualising changes in the performance of the model as a function of either a specified property of the data (one that is explicitly passed to the model at test time, or other metadata); or a perturbation of the data, parameterised by a value *p*. In addition, ``Sensie`` offers a convenience method to display the sensitivity of a classification model to the class itself, though this should be used with caution.

``Sensie``'s algorithm is summarised as follows as applied to a network trained for classification, for two use cases:

**A)** Using data or metadata: a scalar property $p$ of the test set $X_{\textrm{test}}$, such that each example $\boldsymbol{x}_i$ has a corresponding value $p_i$ and a supplied ground truth $\hat{y}$:

1. Segment the data into bins by $p$.
2. Collect predictions from the model (using ``model.predict()`` or a user-supplied function), and note the scores for the correct classes $\hat{y}$.
3. Calculate the mean score $\bar{s}$ in the correct class for each bin, and the standard deviation.
4. Plot $\bar{s}$ as a function of $p$.
5. Estimate the significance of the effect using Bayesian linear regression, producing a scalar value representing $\partial \bar{s}/\partial p$, with 50% and 95% credible intervals.

**B)** Using a perturber function $f_{\textrm{perturb}}(\boldsymbol{x}, p)$---an arbitrary transformation of the inputs---applied to the test set, where the magnitude of perturbation is parameterised by $p$:

1. Choose discrete values of $p$ to be tested.
2. For each value of $p$ to be tested, $p_j$, transform the test set such that $X'_j = f_{\textrm{perturb}}(X_\textrm{test}, p_j)$.
3. Collect predictions from the model for $X'_j$ (using ``model.predict()`` or a user-supplied method), and note the scores for the correct classes, $\hat{y}$.
4. Calculate the mean score $\bar{s}$ in the correct class for each $p_j$, and the standard deviation.
5. Plot $\bar{s}$ as a function of $p$.
6. Estimate the significance of the effect using Bayesian linear regression[^bayesian], producing a scalar value representing $\partial \bar{s}/\partial p$, with 50% and 95% credible intervals.

In the first case, A), ``Sensie`` can optionally forego binning by P, and treat every element as a data point in determining the trend.


![Left: Output from ``Sensie`` for a model trained to recognise handwritten digits, testing model sensitivity to rotation. Error bars show the standard deviation for the mean ground-truth-class score. Right: Sensitivity of a model to an applied blur of the input image data, showing a linear fit to a significant region.](sensie_examples.png)

# Discussion and conclusion

``Sensie`` helps the model user to quickly identify sensitivities of the accuracy to various properties of the input. This may identify deficiencies in the training set; for instance, if we expect an image classifier to be robust under a particular transformation such as rotation or translation, ``Sensie`` can quickly quantify whether such a transformation has a statistically significant effect on the accuracy across the test set. Figure 1 shows two examples; on the left, an output plot from ``Sensie`` showing the sensitivity of a trained model to rotation of the input images; the fact that the curve is not flat indicates that the model is highly sensitive to input orientation and that this is a significant feature (or bug) of the training process. On the right, the sensitivity of a trained model to a Gaussian blur applied to the model, along with a fit indicating the magnitude of the effect. The sensitivity to any given property, i.e. the curve $\bar{s}$ vs $p$, may not be linear. In this case, a fit using linear regression may produce a less meaningful result. 

``Sensie`` has been successfully applied to a model trained for astronomical classification [@jacobsExtendedCatalogGalaxy2019] to indicate the model's sensitivity to the scale of features present in the input; the signal-to-noise ratio of the input data; the colours of the input; and astrophysical properties of the objects present in an input image.

# Acknowledgements

The author acknowledges support from Karl Glazebrook’s Australian Research Council Laureate Fellowship FL180100060.

# References

[^bayesian]: See https://docs.pymc.io/notebooks/GLM-linear.html for an example using the PyMC3 probabilistic programming framework.
[^tfe]: https://github.com/sicara/tf-explain

