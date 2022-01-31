---
title: 'Foolbox Native: Fast adversarial attacks to benchmark the robustness of machine learning models in PyTorch, TensorFlow, and JAX'
tags:
  - python
  - machine learning
  - adversarial attacks
  - neural networks
  - pytorch
  - tensorflow
  - jax
  - keras
  - eagerpy
authors:
  - name: Jonas Rauber
    orcid: 0000-0001-6795-9441
    affiliation: "1, 2"
  - name: Roland Zimmermann
    affiliation: "1, 2"
  - name: Matthias Bethge^[joint senior authors]
    affiliation: "1, 3"
  - name: Wieland Brendel
    affiliation: "1, 3"
affiliations:
 - name: Tübingen AI Center, University of Tübingen, Germany
   index: 1
 - name: International Max Planck Research School for Intelligent Systems, Tübingen, Germany
   index: 2
 - name: Bernstein Center for Computational Neuroscience Tübingen, Germany
   index: 3
date: 10 August 2020
bibliography: paper.bib
---

# Summary

Machine learning has made enormous progress in recent years and is now being used in many real-world applications. Nevertheless, even state-of-the-art machine learning models can be fooled by small, maliciously crafted perturbations of their input data. Foolbox is a popular Python library to benchmark the robustness of machine learning models against these adversarial perturbations. It comes with a huge collection of state-of-the-art adversarial attacks to find adversarial perturbations and thanks to its framework-agnostic design it is ideally suited for comparing the robustness of many different models implemented in different frameworks. Foolbox 3 aka Foolbox Native has been rewritten from scratch to achieve native performance on models developed in PyTorch [@pytorch], TensorFlow [@tensorflow], and JAX [@jax], all with one codebase without code duplication.

# Statement of need

Evaluating the adversarial robustness of machine learning models is crucial to understanding their shortcomings and quantifying the implications on safety, security, and interpretability. Foolbox Native is the first adversarial robustness toolbox that is both fast and framework-agnostic. This is important because modern machine learning models such as deep neural networks are often computationally expensive and are implemented in different frameworks such as PyTorch and TensorFlow. Foolbox Native combines the framework-agnostic design of the original Foolbox [@rauber2017foolbox] with real batch support and native performance in PyTorch, TensorFlow, and JAX, all using a single codebase without code duplication. To achieve this, all adversarial attacks have been rewritten from scratch and now use EagerPy [@rauber2020eagerpy] instead of NumPy [@numpy] to interface *natively* with the different frameworks.

This is great for both users and developers of adversarial attacks. Users can efficiently evaluate the robustness of different models in different frameworks using the same set of state-of-the-art adversarial attacks, thus obtaining comparable results. Attack developers do not need to choose between supporting just one framework or reimplementing their new adversarial attack multiple times and dealing with code duplication. In addition, they both benefit from the comprehensive type annotations [@pep484] in Foolbox Native to catch bugs even before running their code.

The combination of being framework-agnostic and simultaneously achieving native performance sets Foolbox Native apart from other adversarial attack libraries. The most popular alternative to Foolbox is CleverHans^[https://github.com/tensorflow/cleverhans]. It was the first adversarial attack library and has traditionally focused solely on TensorFlow (plans to make it framework-agnostic *in the future* have been announced). The original Foolbox was the second adversarial attack library and the first one to be framework-agnostic. Back then, this was achieved at the expense of performance. The adversarial robustness toolbox ART^[https://github.com/Trusted-AI/adversarial-robustness-toolbox] is another framework-agnostic adversarial attack library, but it is conceptually inspired by the original Foolbox and thus comes with the same performance trade-off. AdverTorch^[https://github.com/BorealisAI/advertorch] is a popular adversarial attack library that was inspired by the original Foolbox but improved its performance by focusing soley on PyTorch. Foolbox Native is our attempt to improve the performance of Foolbox without sacrificing the framework-agnostic design that is crucial to consistently evaluate the robustness of different machine learning models that use different frameworks.

# Use Cases

Foolbox was designed to make adversarial attacks easy to apply even without expert knowledge. It has been used in numerous scientific publications and has already been cited more than 220 times. On GitHub it has received contributions from several developers and has gathered more than 1.500 stars. It provides the reference implementations of various adversarial attacks, including the Boundary Attack [@brendel2018decisionbased], the Pointwise Attack [@schott2018towards], clipping-aware noise attacks [@rauber2020fast], the Brendel Bethge Attack [@brendel2019accurate], and the HopSkipJump Attack [@chen2020hopskipjumpattack], and is under active development since 2017.

# Acknowledgements

J.R. acknowledges support from the Bosch Research Foundation (Stifterverband, T113/30057/17) and the International Max Planck Research School for Intelligent Systems (IMPRS-IS). This work was supported by the German Federal Ministry of Education and Research (BMBF): Tübingen AI Center, FKZ: 01IS18039A, and by the Intelligence Advanced Research Projects Activity (IARPA) via Department of Interior/Interior Business Center (DoI/IBC) contract number D16PC00003. The U.S. Government is authorized to reproduce and distribute reprints for Governmental purposes notwithstanding any copyright annotation thereon. Disclaimer: The views and conclusions contained herein are those of the authors and should not be interpreted as necessarily representing the official policies or endorsements, either expressed or implied, of IARPA, DoI/IBC, or the U.S. Government.

We thank all contributors to Foolbox, in particular Behar Veliqi, Evgenia Rusak, Jianbo Chen, Rene Bidart, Jerome Rony, Ben Feinstein, Eric R Meissner, Lars Holdijk, Lukas Schott, Carl-Johann Simon-Gabriel, Apostolos Modas, William Fleshman, Xuefei Ning, [and many others](https://github.com/bethgelab/foolbox/graphs/contributors).

# References

---
home: true
heroImage: /logo.png
heroText: Foolbox
tagline: "Foolbox Native: Fast adversarial attacks to benchmark the robustness of machine learning models in PyTorch, TensorFlow, and JAX"
actionText: Get Started →
actionLink: /guide/
features:
- title: Native Performance
  details: Foolbox 3 is built on top of EagerPy and runs natively in PyTorch, TensorFlow, and JAX.
- title: State-of-the-art attacks
  details: Foolbox provides a large collection of state-of-the-art gradient-based and decision-based adversarial attacks.
- title: Type Checking
  details: Catch bugs before running your code thanks to extensive type annotations in Foolbox.
footer: Copyright © 2020 Jonas Rauber

---

### What is Foolbox?

**Foolbox** is a **Python library** that lets you easily run adversarial attacks against machine learning models like deep neural networks. It is built on top of [**EagerPy**](https://eagerpy.jonasrauber.de) and works natively with models in [**PyTorch**](https://pytorch.org), [**TensorFlow**](https://www.tensorflow.org), and [**JAX**](https://github.com/google/jax).

```python
import foolbox as fb

model = ...
fmodel = fb.PyTorchModel(model)

attack = fb.attacks.LinfPGD()
epsilons = [0.0, 0.001, 0.01, 0.03, 0.1, 0.3, 0.5, 1.0]
advs, _, success = attack(fmodel, images, labels, epsilons=epsilons)
```
# Development

::: tip NOTE
The following is only necessary if you want to contribute features or
adversarial attacks to Foolbox. As a user of Foolbox, you can just do a normal
[installation](./getting-started).
:::

## Installation

First clone the repsository using `git`:

```bash
git clone https://github.com/bethgelab/foolbox
```

You can then do an editable installation using `pip -e`:

```bash
cd foolbox
pip3 install -e .
```

::: tip
Create a new branch for each new feature or contribution.
This will be necessary to open a pull request later.
:::

## Coding Style

We follow the [PEP 8 Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/).
We use [black](https://github.com/psf/black) for automatic code formatting.
In addition, we use [flake8](https://flake8.pycqa.org/en/latest/) to detect
certain PEP 8 violations.

::: tip
Have a look at the `Makefile`. It contains many useful commands, e.g. `make black` or `make flake8`.
:::

## Type annotions and MyPy

Foolbox uses Python type annotations introduced in [PEP 484](https://www.python.org/dev/peps/pep-0484/).
We use [mypy](http://mypy-lang.org) for static type checking with relatively
strict settings. All code in Foolbox has to be type annotated.

We recommend to run MyPy or a comparable type checker automatically in your
editor (e.g. VIM) or IDE (e.g. PyCharm). You can also run MyPy from the
command line:

```bash
make mypy  # run this in the root folder that contains the Makefile
```

::: tip NOTE
`__init__` methods in Foolbox should not have return type annotations unless
they have no type annotated arguments (i.e. only `self`), in which case
the return type of `__init__` should be specifed as `None`.
:::

## Creating a pull request on GitHub

First, fork the [Foolbox repository on GitHub](https://github.com/bethgelab/foolbox).
Then, add the fork to your local GitHub repository:

```bash
git remote add fork https://github.com/YOUR USERNAME/foolbox
```

Finally, push your new branch to GitHub and open a pull request.
# Getting Started

## Installation

You can install the latest release from [PyPI](https://pypi.org/project/foolbox/) using `pip`:

```bash
python3 -m pip install foolbox
```

Foolbox requires Python 3.6 or newer. To use it with [PyTorch](https://pytorch.org), [TensorFlow](https://www.tensorflow.org), or [JAX](https://github.com/google/jax), the respective framework needs to be installed separately. These frameworks are not declared as dependencies because not everyone wants to use and thus install all of them and because some of these packages have different builds for different architectures and CUDA versions. Besides that, all essential dependencies are automatically installed.

::: warning NOTE
Foolbox requires Python 3.6 or newer.
:::

## Getting a Model

Once Foolbox is installed, you need to turn your PyTorch, TensorFlow, or JAX model into a Foolbox model.

### PyTorch

For PyTorch, you simply instantiate your `torch.nn.Module` and then pass it
to `fb.PyTorchModel`. Here we use a pretrained ResNet-18 from `torchvision`.
Additionally, you should specify the preprocessing expected by the model
(e.g. subtracting `mean`, and dividing by `std`, along the third axis from the back)
and the bounds of the input space (before the preprocessing).

```python
# PyTorch ResNet18
import torch
import torchvision
model = torchvision.models.resnet18(pretrained=True)
preprocessing = dict(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225], axis=-3)
bounds = (0, 1)
fmodel = fb.PyTorchModel(model, bounds=bounds, preprocessing=preprocessing)
```

### TensorFlow

For TensorFlow, you simply instantiate your `tf.keras.Model` and then pass it
to `fb.TensorFlowModel`. Here we show three examples using pretrained ImageNet
models. Additionally, you should specify the preprocessing expected by the model
(e.g. flipping an axis, here from RGB to BGR, subtracting `mean`,
and dividing by `std`, along the third axis from the back)
and the bounds of the input space (before the preprocessing).

```python
# TensorFlow ResNet50
import tensorflow as tf
model = tf.keras.applications.ResNet50(weights="imagenet")
preprocessing = dict(flip_axis=-1, mean=[103.939, 116.779, 123.68])
bounds = (0, 255)
fmodel = fb.TensorFlowModel(model, bounds=bounds, preprocessing=preprocessing)
```

```python
# TensorFlow ResNet50V2
import tensorflow as tf
model = tf.keras.applications.ResNet50V2(weights="imagenet")
preprocessing = dict()
bounds = (-1, 1)
fmodel = fb.TensorFlowModel(model, bounds=bounds, preprocessing=preprocessing)
```

```
# TensorFlow MobileNetV2
import tensorflow as tf
model = tf.keras.applications.MobileNetV2(weights="imagenet")
preprocessing = dict()
bounds = (-1, 1)
fmodel = fb.TensorFlowModel(model, bounds=bounds, preprocessing=preprocessing)
```

### JAX

For JAX, you simply specify your model as a callable, i.e. an instance of a
class with a `__call__` method or a simple function. It should take an
input array and return the array with predictions. You can then pass this
callable to `fb.JAXModel`. Additionally, you should specify the
preprocessing (see previous examples) and
the bounds of the input space (before the preprocessing).

```python
class Model:
    def __call__(self, x):
        # turn the inputs x into predictions y
        y = x  # replace with your real model
        return y

model = Model()
bounds = (0, 1)
fmodel = fbn.JAXModel(model, bounds)
```

## Transform Bounds

Next you can optionally transform the bounds of the input space of our model.
In the following, we want to work with a model that has (0, 1) bounds.

```python
fmodel = fmodel.transform_bounds((0, 1))
```

If your model already had bounds `(0, 1)`, this does not do anything.
But if your model had different bounds, e.g. `(0, 255)` this would adjust
the preprocessing accordingly such that your model now expects inputs
between `0` and `1`. This is particularly useful if you work with
different models that have different bounds.

## Dataset

Before we can attack our model, we first need some data.
For convenience, Foolbox comes with helper functions that provide
a small set of sample images from different computer vision datasets.

```python
images, labels = fb.utils.samples(fmodel, dataset='imagenet', batchsize=16)
```

Note that images and labels should be a batch of native tensors, i.e.
PyTorch tensors, TensorFlow tensors, or JAX arrays, depending on which framework
you use.

## Attacking the Model

Now we have everything ready to attack the model. Before we do that,
we will quickly check its clean accuracy on our evaluation set.

```python
fb.utils.accuracy(fmodel, images, labels)
# -> 0.9375 (depends on the model!)
```

To run an attack, we first instantiate the corresponding class.

```python
attack = fb.attacks.LinfDeepFoolAttack()
```

And finally we can apply the attack on our model by passing
the input tensor (here `images`), the corresponding true `labels`,
and one or more `epsilons`.

```python
raw, clipped, is_adv = attack(fmodel, images, labels, epsilons=0.03)
```

The attack returns three tensors.

1. The raw adversarial examples. This depends on the attack and we cannot make an guarantees about this output.
2. The clipped adversarial examples. These are guaranteed to not be perturbed more than epsilon and thus are the actual adversarial examples you want to visualize. Note that some of them might not actually switch the class. To know which samples are actually adversarial, you should look at the third tensor.
3. The third tensor contains a boolean for each sample, indicating which samples are true adversarials that are both misclassified and within the epsilon balls around the clean samples.

How to use these tensors will become more clear in a moment.

## Multiple Epsilons

Usually, you should not just look at a single epsilon, but at many different epislons from small to large.
The most efficient way to obtain the corresponding results is by
running the attack with multiple epsilons. It will automatically
select the right strategy depending on the type of attack.

```python
import numpy as np
epsilons = np.linspace(0.0, 0.005, num=20)
```

Let's rerun the attack for all `epsilons`.

```python
raw, clipped, is_adv = attack(fmodel, images, labels, epsilons=epsilons)
```

The returned tensors, `raw`, `clipped`, and `is_adv` now have an additional batch dimension for the different `epsilons`.

## Robust Accuracy

You can now obtain the robust accuracy by simply averaging `is_adv`.

```python
robust_accuracy = 1 - is_adv.float32().mean(axis=-1)
```

You can now plot the robust accuracy using Matplotlib.

```python
import matplotlib.pyplot as plt
plt.plot(epsilons, robust_accuracy.numpy())
```

You can also visualize the adversarials using `fb.plot.images`.

## Learn More

To learn more, have a look at our [Tutorial](https://github.com/jonasrauber/foolbox-native-tutorial/blob/master/foolbox-native-tutorial.ipynb),
the [examples](./examples.md), the [API docs](https://foolbox.readthedocs.io/en/stable/) and of course the [README](https://github.com/bethgelab/foolbox).

# Introduction

## What is Foolbox?

**Foolbox** is a **Python library** that lets you easily run adversarial attacks against machine learning models like deep neural networks. It is built on top of [**EagerPy**](https://eagerpy.jonasrauber.de) and works natively with models in [**PyTorch**](https://pytorch.org), [**TensorFlow**](https://www.tensorflow.org), and [**JAX**](https://github.com/google/jax).
# Adding Adversarial Attacks

::: tip NOTE
The [development guidelines](./development) explain how to get started with
with developing features and adversarial attacks for Foolbox.
:::

## The `Attack` base class

Adversarial attacks in Foolbox should either directly or indirectly subclass
the `Attack` base class in `foolbox/attacks/base.py`.

An attack in Foolbox needs to implement two methods, `__call__` and `repeat`.

Both methods need to be implemented with the same signature as the base class.
The type annotation for the `criterion` argument of `__call__` can be made
more precise, see `foolbox/attacks/carlini_wagner.py` for an example.

The `__call__` method should return three values, a list of raw tensors (one
for each epsilon) with the internal raw attack results, a list of tensors
corresponding to the raw tensors but with perturbation sizes guaranteed to
be clipped to the given epsilons, and a boolean tensor with `len(epsilons)`
rows and `len(inputs)` columns indicating for each returned sample whether
it is a successful adversarial example given the respective epsilon and
criterion. If `epsilons` is a single scalar epsilon (and not a list with
one element), then the first and second return value should be a tensor
rather than a list and the third return value should be 1-D tensor.

All returned tensors must have the same type as the input tensors. In
particular, native tensors should be returned as native tensors and
EagerPy-wrapped tensors should be returned as EagerPy-wrapped tensors.
Use `astensor_` or `astensors_` and `restore_type`.

The `repeat` method should return a version of the attack that repeats itself
n times and returns the best result.

::: warning NOTE
In practice, it is usually not necessary to subclass `Attack` directly.
Instead, for most attacks it is easiest to subclass either `FixedEpsilonAttack`
or `MinimizatonAttack`.
:::

## The `FixedEpsilonAttack` base class

Attacks that require a fixed epsilon and try to find an adversarial
perturbation given this perturbation budget (e.g. `FGSM` and `PGD`) should
be implemented by subclassing `FixedEpsilonAttack`. It already provides
implementations of `__call__` and `repeat`. The attack just needs
to specify the `distance` property (simply assign a class variable) and
implement the `run` method that gets a single `epsilon` and returns a batch
of perturbed inputs, ideally adversarial and ideally with a perturbation
size smaller or equal to `epsilon`.
The `distance` is used by `__call__` to determine which perturbed inputs
are actually adversarials given `epsilon` and by `repeat` to determine the
run.

## The `MinimizatonAttack` base class

Attacks that try to find adversarial examples with minimal perturbation size
(e.g. the `Carlini & Wagner` attack or the `Boundary Attack`) should
be implemented by subclassing `MinimizatonAttack`. It already provides
implementations of `__call__` and `repeat`. The attack just needs
to specify the `distance` property (simply assign a class variable) and
implement the `run` method that returns a batch of minimally perturbed
adversarials. For `MinimizatonAttack` subclasses, `run` gets called only once
by `__call__` independent of how many `epsilons` are given. The `__call__`
method then compares the minimal adversarial perturbation to the different
epsilons.

::: tip
You should have a look at the implementation of existing attacks
to get an impression of the best practices and conventions used in Foolbox.
:::
---
title: Examples

---

# Examples :tada:

::: tip
More examples can be found in the [examples folder](https://github.com/bethgelab/foolbox/tree/master/examples).
:::

<<< @/../examples/single_attack_pytorch_resnet18.py
## Performance comparison between Foolbox versions

|                                        |   Foolbox `1.8.0`   |   Foolbox `2.4.0`  | Foolbox `3.1.1`<br>(aka Native) |
|----------------------------------------|:-------------------:|:------------------:|:-------------------------------:|
| accuracy (single image)                |  `5.02 ms ± 338 µs` | `4.99 ms ± 378 µs` |      **`3.99 ms ± 131 µs`**     |
| accuracy (16 images)                   | `88.9 ms ± 8.24 ms` |  `12 ms ± 1.34 ms` |     **`8.21 ms ± 54.4 µs`**     |
| PGD attack (16 images, single epsilon) |      `161.8 s`      |      `37.5 s`      |           **`1.1 s`**           |
| PGD attack (16 images, 8 epsilons)     |      `164.6 s`      |      `36.9 s`      |           **`9.0 s`**           |


All experiments were done on an Nvidia GeForce GTX 1080 using the PGD attack.

Note that Foolbox 3 is faster because **1)** it avoids memory copies between GPU
and CPU by using EagerPy instead of NumPy, **2)** it fully supports batches
of inputs, and **3)** it currently uses a different approach for fixed-epsilon attacks
like PGD (instead of minimizing the perturbation, the attack is run for
each epsilon; this is more inline with what is generally expected;
for these attacks the duration therefore now scales with the
number of epsilons; it is however still faster and it produces better results).
## Examples

This folder contains examples that demonstrate how Foolbox can be used
to run one or more adversarial attacks and how to use the returned results
to compute the robust accuracy (the accuracy of the model when it is attacked).

The standard example can be found in:
* `single_attack_pytorch_resnet18.py`
* `single_attack_tensorflow_resnet50.py`

It shows how to run a single adversarial attack (Linf PGD) against an ImageNet
model in PyTorch and TensorFlow.

The remaining examples are all for PyTorch,
but the difference between these frameworks is really just replacing the model
at the beginning of the script. So any example can be easily run with any
framework.

`multiple_attacks_pytorch_resnet18.py` is an extended version of the single attack
example. It shows how to combine the results of running multiple attacks
to report the robust accuracy always using the strongest attack per sample.

`spatial_attack_pytorch_resnet18.py` shows how to use the Spatial Attack. This attack
is a bit special because it doesn't use Lp balls and instead considers translations
and rotations. It therefore has a custom example. All the other attacks can be
used like Linf PGD in the other examples above.

`substituion_model_pytorch_resnet18.py` shows how to replace the gradient of
a model with the gradient of another model. This can be useful when the original
model has bad gradients ("gradient masking", "obfuscated gradients").

The `zoo` folder shows how a model can be shared in a Foolbox Model Zoo compatible way.
.. raw:: html

   <a href="https://foolbox.jonasrauber.de"><img src="https://raw.githubusercontent.com/bethgelab/foolbox/master/guide/.vuepress/public/logo_small.png" align="right" /></a>

.. image:: https://badge.fury.io/py/foolbox.svg
   :target: https://badge.fury.io/py/foolbox

.. image:: https://readthedocs.org/projects/foolbox/badge/?version=latest
    :target: https://foolbox.readthedocs.io/en/latest/

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black

.. image:: https://joss.theoj.org/papers/10.21105/joss.02607/status.svg
   :target: https://doi.org/10.21105/joss.02607

===============================================================================================================================
Foolbox Native: Fast adversarial attacks to benchmark the robustness of machine learning models in PyTorch, TensorFlow, and JAX
===============================================================================================================================

`Foolbox <https://foolbox.jonasrauber.de>`_ is a **Python library** that lets you easily run adversarial attacks against machine learning models like deep neural networks. It is built on top of EagerPy and works natively with models in `PyTorch <https://pytorch.org>`_, `TensorFlow <https://www.tensorflow.org>`_, and `JAX <https://github.com/google/jax>`_.

🔥 Design 
----------

**Foolbox 3** a.k.a. **Foolbox Native** has been rewritten from scratch
using `EagerPy <https://github.com/jonasrauber/eagerpy>`_ instead of
NumPy to achieve native performance on models
developed in PyTorch, TensorFlow and JAX, all with one code base without code duplication.

- **Native Performance**: Foolbox 3 is built on top of EagerPy and runs natively in PyTorch, TensorFlow, and JAX and comes with real batch support.
- **State-of-the-art attacks**: Foolbox provides a large collection of state-of-the-art gradient-based and decision-based adversarial attacks.
- **Type Checking**: Catch bugs before running your code thanks to extensive type annotations in Foolbox.

📖 Documentation
-----------------

- **Guide**: The best place to get started with Foolbox is the official `guide <https://foolbox.jonasrauber.de>`_.
- **Tutorial**: If you are looking for a tutorial, check out this `Jupyter notebook <https://github.com/jonasrauber/foolbox-native-tutorial/blob/master/foolbox-native-tutorial.ipynb>`_ |colab|.
- **Documentation**: The API documentation can be found on `ReadTheDocs <https://foolbox.readthedocs.io/en/stable/>`_.

.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/jonasrauber/foolbox-native-tutorial/blob/master/foolbox-native-tutorial.ipynb

🚀 Quickstart
--------------

.. code-block:: bash

   pip install foolbox

Foolbox requires Python 3.6 or newer. To use it with `PyTorch <https://pytorch.org>`_, `TensorFlow <https://www.tensorflow.org>`_, or `JAX <https://github.com/google/jax>`_, the respective framework needs to be installed separately. These frameworks are not declared as dependencies because not everyone wants to use and thus install all of them and because some of these packages have different builds for different architectures and CUDA versions. Besides that, all essential dependencies are automatically installed.

You can see the versions we currently use for testing in the `Compatibility section <#-compatibility>`_ below, but newer versions are in general expected to work.

🎉 Example
-----------

.. code-block:: python

   import foolbox as fb

   model = ...
   fmodel = fb.PyTorchModel(model, bounds=(0, 1))

   attack = fb.attacks.LinfPGD()
   epsilons = [0.0, 0.001, 0.01, 0.03, 0.1, 0.3, 0.5, 1.0]
   _, advs, success = attack(fmodel, images, labels, epsilons=epsilons)


More examples can be found in the `examples <./examples/>`_ folder, e.g.
a full `ResNet-18 example <./examples/single_attack_pytorch_resnet18.py>`_.

📄 Citation
------------

If you use Foolbox for your work, please cite our `JOSS paper on Foolbox Native <https://doi.org/10.21105/joss.02607>`_ and our `ICML workshop paper on Foolbox <https://arxiv.org/abs/1707.04131>`_ using the following BibTeX entries:

.. code-block::

   @article{rauber2017foolboxnative,
     doi = {10.21105/joss.02607},
     url = {https://doi.org/10.21105/joss.02607},
     year = {2020},
     publisher = {The Open Journal},
     volume = {5},
     number = {53},
     pages = {2607},
     author = {Jonas Rauber and Roland Zimmermann and Matthias Bethge and Wieland Brendel},
     title = {Foolbox Native: Fast adversarial attacks to benchmark the robustness of machine learning models in PyTorch, TensorFlow, and JAX},
     journal = {Journal of Open Source Software}
   }

.. code-block::

   @inproceedings{rauber2017foolbox,
     title={Foolbox: A Python toolbox to benchmark the robustness of machine learning models},
     author={Rauber, Jonas and Brendel, Wieland and Bethge, Matthias},
     booktitle={Reliable Machine Learning in the Wild Workshop, 34th International Conference on Machine Learning},
     year={2017},
     url={http://arxiv.org/abs/1707.04131},
   }


👍 Contributions
-----------------

We welcome contributions of all kind, please have a look at our
`development guidelines <https://foolbox.jonasrauber.de/guide/development.html>`_.
In particular, you are invited to contribute
`new adversarial attacks <https://foolbox.jonasrauber.de/guide/adding_attacks.html>`_.
If you would like to help, you can also have a look at the issues that are
marked with `contributions welcome
<https://github.com/bethgelab/foolbox/issues?q=is%3Aopen+is%3Aissue+label%3A%22contributions+welcome%22>`_.

💡 Questions?
--------------

If you have a question or need help, feel free to open an issue on GitHub.
Once GitHub Discussions becomes publically available, we will switch to that.

💨 Performance
--------------

Foolbox Native is much faster than Foolbox 1 and 2. A basic `performance comparison`_ can be found in the `performance` folder.

🐍 Compatibility
-----------------

We currently test with the following versions:

* PyTorch 1.4.0
* TensorFlow 2.1.0
* JAX 0.1.57
* NumPy 1.18.1

.. _performance comparison: performance/README.md
License
-------

The code in this subfolder might be under a different license than the rest of the project.

Sources
-------

* `clipping_aware_rescaling.py <https://github.com/jonasrauber/clipping-aware-rescaling>`_
Welcome to Foolbox Native
=========================

Foolbox is a Python toolbox to create adversarial examples that fool neural networks.
*Foolbox 3.0* a.k.a. *Foolbox Native* has been completely rewritten from scratch.
It is now built on top of `EagerPy <https://github.com/jonasrauber/eagerpy>`_
and comes with native support for these frameworks:

* `PyTorch <https://pytorch.org>`_
* `TensorFlow <https://www.tensorflow.org>`_
* `JAX <https://github.com/google/jax>`_

Foolbox comes with a :doc:`large collection of adversarial attacks <modules/attacks>`, both gradient-based white-box attacks as well as decision-based and score-based black-box attacks.

The source code and a `minimal working example <https://github.com/bethgelab/foolbox#example>`_ can be found on `GitHub <https://github.com/bethgelab/foolbox>`_.


.. toctree::
   :maxdepth: 2
   :caption: User API

   modules/models
   modules/attacks
   modules/criteria
   modules/distances
   modules/utils
   modules/plot
   modules/zoo

.. toctree::
   :maxdepth: 2
   :caption: Internal API

   modules/devutils
   modules/tensorboard
   modules/types


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
:mod:`foolbox.plot`
===================

.. automodule:: foolbox.plot
   :members:
   :undoc-members:
:mod:`foolbox.tensorboard`
==========================

.. automodule:: foolbox.tensorboard
   :members:
   :undoc-members:
:mod:`foolbox.utils`
====================

.. automodule:: foolbox.utils
   :members:
   :undoc-members:
.. automodule:: foolbox.criteria
:mod:`foolbox.types`
====================

.. automodule:: foolbox.types
   :members:
   :undoc-members:
:mod:`foolbox.models`
=====================

.. automodule:: foolbox.models

Models
------

.. autosummary::
   :nosignatures:

   Model
   PyTorchModel
   TensorFlowModel
   JAXModel
   NumPyModel

Wrappers
--------

.. autosummary::
   :nosignatures:

   TransformBoundsWrapper

Detailed description
--------------------

.. autoclass:: Model
   :members:

.. autoclass:: PyTorchModel
   :members:

.. autoclass:: TensorFlowModel
   :members:

.. autoclass:: JAXModel
   :members:

.. autoclass:: NumPyModel
   :members:

.. autoclass:: TransformBoundsWrapper
   :members:
:mod:`foolbox.attacks`
======================

.. automodule:: foolbox.attacks
   :members:
   :undoc-members:

.. autosummary::
   :nosignatures:

   L2ContrastReductionAttack
   VirtualAdversarialAttack
   DDNAttack
   L2ProjectedGradientDescentAttack
   LinfProjectedGradientDescentAttack
   L2BasicIterativeAttack
   LinfBasicIterativeAttack
   L2FastGradientAttack
   LinfFastGradientAttack

   L2AdditiveGaussianNoiseAttack
   L2AdditiveUniformNoiseAttack
   L2ClippingAwareAdditiveGaussianNoiseAttack
   L2ClippingAwareAdditiveUniformNoiseAttack
   LinfAdditiveUniformNoiseAttack
   L2RepeatedAdditiveGaussianNoiseAttack
   L2RepeatedAdditiveUniformNoiseAttack
   L2ClippingAwareRepeatedAdditiveGaussianNoiseAttack
   L2ClippingAwareRepeatedAdditiveUniformNoiseAttack
   LinfRepeatedAdditiveUniformNoiseAttack
   InversionAttack
   BinarySearchContrastReductionAttack
   LinearSearchContrastReductionAttack

   L2CarliniWagnerAttack
   NewtonFoolAttack
   EADAttack
   GaussianBlurAttack
   L2DeepFoolAttack
   LinfDeepFoolAttack
   SaltAndPepperNoiseAttack
   LinearSearchBlendedUniformNoiseAttack
   BinarizationRefinementAttack
   DatasetAttack
   BoundaryAttack
   L0BrendelBethgeAttack
   L1BrendelBethgeAttack
   L2BrendelBethgeAttack
   LinfinityBrendelBethgeAttack
   L0FMNAttack
   L1FMNAttack
   L2FMNAttack
   LInfFMNAttack

   FGM
   FGSM
   L2PGD
   LinfPGD
   PGD

.. autoclass:: L2ContrastReductionAttack
.. autoclass:: VirtualAdversarialAttack
.. autoclass:: DDNAttack
.. autoclass:: L2ProjectedGradientDescentAttack
.. autoclass:: LinfProjectedGradientDescentAttack
.. autoclass:: L2BasicIterativeAttack
.. autoclass:: LinfBasicIterativeAttack
.. autoclass:: L2FastGradientAttack
.. autoclass:: LinfFastGradientAttack

.. autoclass:: L2AdditiveGaussianNoiseAttack
.. autoclass:: L2AdditiveUniformNoiseAttack
.. autoclass:: L2ClippingAwareAdditiveGaussianNoiseAttack
.. autoclass:: L2ClippingAwareAdditiveUniformNoiseAttack
.. autoclass:: LinfAdditiveUniformNoiseAttack
.. autoclass:: L2RepeatedAdditiveGaussianNoiseAttack
.. autoclass:: L2RepeatedAdditiveUniformNoiseAttack
.. autoclass:: L2ClippingAwareRepeatedAdditiveGaussianNoiseAttack
.. autoclass:: L2ClippingAwareRepeatedAdditiveUniformNoiseAttack
.. autoclass:: LinfRepeatedAdditiveUniformNoiseAttack
.. autoclass:: InversionAttack
.. autoclass:: BinarySearchContrastReductionAttack
.. autoclass:: LinearSearchContrastReductionAttack

.. autoclass:: L2CarliniWagnerAttack
.. autoclass:: NewtonFoolAttack
.. autoclass:: EADAttack
.. autoclass:: GaussianBlurAttack
.. autoclass:: L2DeepFoolAttack
.. autoclass:: LinfDeepFoolAttack
.. autoclass:: SaltAndPepperNoiseAttack
.. autoclass:: LinearSearchBlendedUniformNoiseAttack
.. autoclass:: BinarizationRefinementAttack
.. autoclass:: DatasetAttack
.. autoclass:: BoundaryAttack
.. autoclass:: L0BrendelBethgeAttack
.. autoclass:: L1BrendelBethgeAttack
.. autoclass:: L2BrendelBethgeAttack
.. autoclass:: LinfinityBrendelBethgeAttack
.. autoclass:: L0FMNAttack
.. autoclass:: L1FMNAttack
.. autoclass:: L2FMNAttack
.. autoclass:: LInfFMNAttack

.. autoclass:: FGM
.. autoclass:: FGSM
.. autoclass:: L2PGD
.. autoclass:: LinfPGD
.. autoclass:: PGD
:mod:`foolbox.distances`
========================

.. automodule:: foolbox.distances

Detailed description
--------------------

.. autoclass:: Distance
   :members:

.. autoclass:: LpDistance
   :members:
:mod:`foolbox.zoo`
==================

.. automodule:: foolbox.zoo


Get Model
---------

.. autofunction:: get_model


Fetch Weights
-------------

.. autofunction:: fetch_weights
:mod:`foolbox.devutils`
=======================

.. automodule:: foolbox.devutils
   :members:
   :undoc-members:
