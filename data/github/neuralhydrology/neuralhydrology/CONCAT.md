![#](docs/source/_static/img/neural-hyd-logo-black.png)

Python library to train neural networks with a strong focus on hydrological applications.

This package has been used extensively in research over the last years and was used in various academic publications. 
The core idea of this package is modularity in all places to allow easy integration of new datasets, new model 
architectures or any training-related aspects (e.g. loss functions, optimizer, regularization). 
One of the core concepts of this code base are configuration files, which let anyone train neural networks without
touching the code itself. The NeuralHydrology package is built on top of the deep learning framework 
[PyTorch](https://pytorch.org/), since it has proven to be the most flexible and useful for research purposes.

We (the AI for Earth Science group at the Institute for Machine Learning, Johannes Kepler University, Linz, Austria) are using
this code in our day-to-day research and will continue to integrate our new research findings into this public repository.

- Documentation: [neuralhydrology.readthedocs.io](https://neuralhydrology.readthedocs.io)
- Research Blog: [neuralhydrology.github.io](https://neuralhydrology.github.io)
- Bug reports/Feature requests [https://github.com/neuralhydrology/neuralhydrology/issues](https://github.com/neuralhydrology/neuralhydrology/issues)

# Cite NeuralHydrology

In case you use NeuralHydrology in your research or work, it would be highly appreciated if you include a reference to our [JOSS paper](https://joss.theoj.org/papers/10.21105/joss.04050#) in any kind of publication.

```bibtex
@article{kratzert2022joss,
  title = {NeuralHydrology --- A Python library for Deep Learning research in hydrology},
  author = {Frederik Kratzert and Martin Gauch and Grey Nearing and Daniel Klotz},
  journal = {Journal of Open Source Software},
  publisher = {The Open Journal},
  year = {2022},
  volume = {7},
  number = {71},
  pages = {4050},
  doi = {10.21105/joss.04050},
  url = {https://doi.org/10.21105/joss.04050},
}
```

# Contact

For questions or comments regarding the usage of this repository, please use the [discussion section](https://github.com/neuralhydrology/neuralhydrology/discussions) on Github. For bug reports and feature requests, please open an [issue](https://github.com/neuralhydrology/neuralhydrology/issues) on GitHub.
In special cases, you can also reach out to us by email: neuralhydrology(at)googlegroups.com
---
name: Bug report
about: Create a report to help us improve
title: "[BUG] Title describing the bug"
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior. E.g. which data did you use, what were the commands that you executed, did you modify the code, etc.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Description of the dataset**
Provide details on the dataset that you use. Is it e.g. an out-of-the-box CAMELS dataset, or did you create your own csv/netCDF files. In that case, providing a data sample might be beneficial.

**Logs & Screenshots**
Please provide the full stack trace if any exception occured. If applicable, add screenshots to help explain your problem.

**Desktop & Environment (please complete the following information):**
 - OS: [e.g. Linux, Windows, iOS]
 - The git commit if you cloned the repo, or the version number in `neuralhydrology/__about__.py`)
 - The Python version and a list of installed Python packages. If you use conda, you can create this list via `conda env export`.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: "[REQUEST] Title describing the feature request"
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
Contributing to NeuralHydrology
===============================

Welcome to the NeuralHydrology contribution guide!
--------------------------------------------------

First off, thank you for considering contributing to NeuralHydrology!

The following is a set of guidelines for contributing to NeuralHydrology. They are intended to make contributions as easy as possible for both you as a contributor and us as maintainers of the project.


Kinds of contribution
---------------------

There are many ways to contribute to this project. To name a few:

- File bug reports or feature requests as GitHub issues
- Contribute a bug fix as a pull request
- Contribute a new feature (new model, new dataset, etc.) as a pull request
- Improve the `documentation <https://neuralhydrology.readthedocs.io/>`__ with new tutorials, better descriptions, or even just by fixing typos

The following sections will give some more guidance for each of these options.

Reporting Bugs and feature requests
-----------------------------------
If you've come across some behavior that you think is a bug or a missing feature, the first step is to check if it's already known.
For this, take a look at the `GitHub issues page <https://github.com/neuralhydrology/neuralhydrology/issues>`__.

If you can't find anything related there, you're welcome to create your own issue.

**NOTE**: Please use the `issues page <https://github.com/neuralhydrology/neuralhydrology/issues>`__ for bug reports and feature requests, and leave the `discussions <https://github.com/neuralhydrology/neuralhydrology/discussions>`__ for more open discussions, e.g., project ideas or engaging with other users.


Creating an issue for a bug report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you couldn't find an open issue addressing the problem, open a new one. Please make sure to use a meaningful title and a clear description that includes as much relevant details as possible. Ideally, a bug report should contain:

- A description of how to reproduce the bug (which commands did you execute, did you do any modifications to the code, etc.)
- The full stack trace if any exceptions occurred
- The configuration yaml file you used in your experiments.
- The NeuralHydrology version you used in your experiments (the git commit if you cloned the repo, or the version number in ``neuralhydrology/__about__.py``)
- A description of the dataset you are using in your experiments: Are you using an out-of-the-box CAMELS dataset, or did you create your own csv/netCDF files?
- A list of Python packages you have installed in your environment. If you're using conda, you can create this list via ``conda env export``.

For console output and config files (stack traces, environment files, ...), please try to format them as code directly in the issue (using GitHub's `code block <https://docs.github.com/en/github/writing-on-github/working-with-advanced-formatting/creating-and-highlighting-code-blocks>`__ syntax) or as a link to a `gist <https://gist.github.com/discover>`__.


Contributing code with a pull request
-------------------------------------
If you wrote some code to fix a bug or provide some new functionality, we use use pull requests (PRs) to merge these changes into the project.

In a PR, please provide a small description of the changes you made, including a reference to the issue that your PR solves (if one exists).
We'll take a look at your PR and do a small code review, where we might outline changes that are necessary before we can merge your PR into the project.
Please don't feel overwhelmed by our change requests---it's entirely normal to have a code review go back and forth a few times before the code is ready to be merged.

Once everyone is happy with the changes in the PR, we'll approve and merge it into the master branch.


Environment setup
~~~~~~~~~~~~~~~~~

The `quick start section <https://neuralhydrology.readthedocs.io/en/latest/usage/quickstart.html>`__ contains information on how to set up the environment (installing the required dependencies and the package itself, as well as downloading data).

To make code changes, you'll need to fork the repository on GitHub.
Describing the process of committing and pushing changes in a git repository goes beyond the scope of this guide.
For more information on how to create pull requests, see for example `this tutorial <https://makeapullrequest.com/>`__ or the `GitHub documentation <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`__.


Code style
~~~~~~~~~~

If possible, we'd be happy if your code is already formatted according to our code style. We use `yapf <https://github.com/google/yapf>`__ with `this configuration file <https://github.com/neuralhydrology/neuralhydrology/blob/master/.style.yapf>`__ to ensure that all our code uses the same line length, indentation style, etc.
If you're using Visual Studio Code as your IDE, you can add the following lines to your settings file (``projectRoot/.vscode/settings.json``):

.. code-block::

    "editor.formatOnSave": true,
    "python.formatting.provider": "yapf",

These settings will make sure that whenever you save changes, yapf will auto-format the changed files according to our style guide.

If you can't get yapf running, don't worry---we can do the formatting as a last step before merging the PR ourselves.


Documentation
~~~~~~~~~~~~~

The `NeuralHydrology documentation <https://neuralhydrology.readthedocs.io/>`__ is automatically generated from the files located in the `docs/` folder.
If your PR introduces any changes that should be documented, you'll find the relevant files there.
Particularly, for new configuration arguments, please add a description of each argument to the `"Configuration Arguments" page <https://github.com/neuralhydrology/neuralhydrology/blob/master/docs/source/usage/config.rst>`__.

To build the documentation webpage locally and see if everything is rendered correctly, simply run ``make html`` in the docs folder. Please check the output of this command to make sure there are no broken links, etc. (these will output red warnings or errors).
Afterwards, you can start a local server via ``python -m http.server`` in ``docs/build/html/``. Once the server is running, you can access the local documentation pages in your browser at the URL ``localhost:8000``.


Continuous Integration Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our GitHub repository contains a suite of test cases that are run automatically on every pull request. These test cases must run successfully before we can merge your contribution. You'll see the status of these test cases at the bottom of your pull request.
Ideally, you've made sure that all tests run smoothly on your local machine before submitting the PR. The test cases use the `pytest <https://docs.pytest.org/>`__ test framework, so you can run all test cases via ``python -m pytest test``, or a specific test case via ``python -m pytest test/test_config_runs.py -k "test_daily_regression_nan_targets"``. If you're using an IDE, it might also auto-detect the test cases for you and provide an option to run them via the IDE's GUI.

Ideally, if you create a PR that fixes a bug or adds some new feature, this PR would also contain one or more test cases that make sure the bug doesn't happen again, or that we don't do changes in the future that would break your new feature. You can find examples for test cases in the ``test/`` folder. If you don't know how to write test cases yourself, don't worry---you can still create a PR and we'll discuss how to create tests during the code review.
.. include:: ../../CONTRIBUTING.rst.. NeuralHydrology documentation master file, created by
   sphinx-quickstart on Mon Aug 17 14:00:15 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NeuralHydrology's documentation!
===========================================

This is the documentation for the NeuralHydrology Python package.
The source code is available on `GitHub <https://github.com/neuralhydrology/neuralhydrology/>`_.

On this documentation page, you'll find a :doc:`quickstart guide <usage/quickstart>` with step-by-step instructions on installation, required datasets, and command-line usage.
We've also written a few :doc:`tutorials <tutorials/index>` that walk you through code examples to train your first models.
The :doc:`modelzoo <usage/models>` lists all the models that come pre-packaged on installation (but you can add your own models, too).
Finally, the :doc:`API docs <api/neuralhydrology>` show in-depth information on all modules, classes, and functions within NeuralHydrology.

You might also be interested in our `Research Blog <https://neuralhydrology.github.io/>`_, where once in a while we post news about our papers and other projects.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage/quickstart
   usage/models
   tutorials/index
   usage/config
   api/neuralhydrology
   contributing
Modelzoo
========

The following section gives an overview of all implemented models in NeuralHydrology. Conceptually, all models in our package consist of two parts, the model class (which constitutes the core of the model as such) and the model heads (which relate the outputs of the model class to the predicted variables). The section `Model Heads`_ provides a list of all implemented model heads, and the section `Model Classes`_ a list of all model classes. If you want to implement your own model within the package you best start at the section `Implementing a new model`_, which provides the necessary details to do so. 


Model Heads
-----------
The head of the model is used on top of the model class and relates the outputs of the `Model Classes`_ to the predicted variable. Currently four model heads are available: `Regression`_, `GMM`_, `CMAL`_ and `UMAL`_. The latter three heads provide options for probabilistic modelling. A detailed overview can be found in `Klotz et al. "Uncertainty Estimation with Deep Learning for Rainfall-Runoff Modelling" <https://arxiv.org/abs/2012.14295>`__. 

Regression
^^^^^^^^^^
:py:class:`neuralhydrology.modelzoo.head.Regression` provides a single layer *regression* head, that includes different activation options for the output. (namely a linear, relu and softplus). 

It is possible to obtain probabilistic predictions with the regression head by using Monte-Carlo Dropout. Its usage is defined in the config.yml by setting ``mc_dropout``. The sampling behavior is governed by picking the number of samples (``n_samples``) and the approach for handling negative samples (``negative_sample_handling``).   

GMM
^^^
:py:class:`neuralhydrology.modelzoo.head.GMM` implements a *Gaussian Mixture Model* head. That is, a mixture density network with Gaussian distributions as components. Each Gaussian component is defined by two parameters (the mean, the variance) and by a set of weights. The current implementation of the GMM head uses two layers. Specific output activations are used for the variances (:py:func:`torch.exp`) and the weights (:py:func:`torch.softmax`).

The number of components can be set in the config.yml using ``n_distributions``. Additionally, the sampling behavior (for the inference) is defined with config.yml by setting the number of samples (``n_samples``), and the approach for handling negative samples (``negative_sample_handling``).  

CMAL
^^^^
:py:class:`neuralhydrology.modelzoo.head.CMAL` implements a *Countable Mixture of Asymmetric Laplacians* head. That is, a mixture density network with asymmetric Laplace distributions as components. The name is a homage to `UMAL`_, which provides an uncountable extension. The CMAL components are defined by three parameters (location, scale, and asymmetry) and linked by a set of weights. The current implementation of the CMAL head uses two layers. Specific output activations are used for the component scales (:py:class:`torch.nn.Softplus(2)`), the asymmetries (:py:func:`torch.sigmoid`), and the weights (:py:func:`torch.softmax`). In our preliminary experiments this heuristic achieved better results. 

The number of components can be set in the config.yml using ``n_distributions``. Additionally, one can sample from CMAL. The behavior of which is defined by setting the number of samples (``n_samples``), and the approach for handling negative samples (``negative_sample_handling``).  

UMAL
^^^^
:py:class:`neuralhydrology.modelzoo.head.CMAL` implements an *Uncountable Mixture of Asymmetric Laplacians* head. That is, a mixture density network that uses an uncountable amount of asymmetric Laplace distributions as components. The *uncountable property* is achieved by implicitly learning the conditional density and approximating it, when needed, with a Monte-Carlo integration, using sampled asymmetry parameters. The UMAL components are defined by two parameters (the location and the scale) and linked by a set of weights. The current implementation uses two hidden layers. The output activation for the scale has some major differences to the original implementation, since it is upper bounded (using :py:func:`0.5*torch.sigmoid`).

During inference the number of components and weights used for the Monte-Carlo approximation are defined in the config.yml by ``n_taus``. The additional argument ``umal_extend_batch`` allows to explicitly account for this integration step during training by repeatedly sampling the asymmetry parameter and extending the batch by ``n_taus``. Furthermore, depending on the used output activation the sampling of the asymmetry parameters can yield unwarranted model behavior. Therefore the lower- and upper-bounds of the sampling can be adjusted using the ``tau_down`` and ``tau_up`` options in the config yml. 
The sampling for UMAL is defined by choosing the number of samples (``n_samples``), and the approach for handling negative samples (``negative_sample_handling``).  


Model Classes
-------------

BaseModel
^^^^^^^^^
Abstract base class from which all models derive. Do not use this class for model training.

CudaLSTM
^^^^^^^^
:py:class:`neuralhydrology.modelzoo.cudalstm.CudaLSTM` is a network using the standard PyTorch LSTM implementation.
All features (``x_d``, ``x_s``, ``x_one_hot``) are concatenated and passed to the network at each time step.
If ``statics/dynamics_embedding`` are used, the static/dynamic inputs will be passed through embedding networks before
being concatenated.
The initial forget gate bias can be defined in config.yml (``initial_forget_bias``) and will be set accordingly during
model initialization.

CustomLSTM
^^^^^^^^^^
:py:class:`neuralhydrology.modelzoo.customlstm.CustomLSTM` is a variant of the ``CudaLSTM``
that returns all gate and state activations for all time steps. This class is mainly implemented for exploratory
reasons. You can use the method ``model.copy_weights()`` to copy the weights of a ``CudaLSTM`` model
into a ``CustomLSTM`` model. This allows to use the fast CUDA implementations for training, and only use this class for
inference with more detailed outputs. You can however also use this model during training (``model: customlstm`` in the
config.yml) or as a starter for your own modifications to the LSTM cell. Note, however, that the runtime of this model
is considerably slower than its optimized counterparts.

EA-LSTM
^^^^^^^
:py:class:`neuralhydrology.modelzoo.ealstm.EALSTM` is an implementation of the Entity-Aware LSTM, as introduced in
`Kratzert et al. "Towards learning universal, regional, and local hydrological behaviors via machine learning applied to large-sample datasets" <https://hess.copernicus.org/articles/23/5089/2019/hess-23-5089-2019.html>`__.
The static features (``x_s`` and/or ``x_one_hot``) are used to compute the input gate activations, while the dynamic
inputs ``x_d`` are used in all other gates of the network.
The initial forget gate bias can be defined in config.yml (``initial_forget_bias``).
If ``statics/dynamics_embedding`` are used, the static/dynamic inputs will first be passed through embedding networks.
The output of the static embedding network will then be passed through the input gate, which consists of a single linear
layer.

EmbCudaLSTM
^^^^^^^^^^^
.. deprecated:: 0.9.11-beta
   Use `CudaLSTM`_ with ``statics_embedding``.

:py:class:`neuralhydrology.modelzoo.embcudalstm.EmbCudaLSTM` is similar to `CudaLSTM`_,
with the only difference that static inputs (``x_s`` and/or ``x_one_hot``) are passed through an embedding network
before being concatenated to the dynamic inputs ``x_d`` at each time step.

GRU
^^^
:py:class:`neuralhydrology.modelzoo.gru.GRU` is a network using the standard PyTorch GRU implementation.
All features (``x_d``, ``x_s``, ``x_one_hot``) are concatenated and passed to the network at each time step.
If ``statics/dynamics_embedding`` are used, the static/dynamic inputs will be passed through embedding networks before
being concatenated.

.. _MC-LSTM:

MC-LSTM
^^^^^^^
:py:class:`neuralhydrology.modelzoo.mclstm.MCLSTM` is a concept for a mass-conserving model architecture inspired by the
LSTM that was recently proposed by `Hoedt et al. (2021) <https://arxiv.org/abs/2101.05186>`_. The implementation included
in this library is the exact model configuration that was used for the hydrology experiments in the linked publication 
(for details, see Appendix B.4.2).
The inputs for the model are split into two groups: i) the mass input, whose values are stored in the memory cells of 
the model and from which the target is calculated and ii) auxiliary inputs, which are used to control the gates 
within the model. In this implementation, only a single mass input per timestep (e.g. precipitation) is allowed, which
has to be specified with the config argument ``mass_inputs``. Make sure to exclude the mass input feature, as well as
the target variable from the standard feature normalization. This can be done using the ``custom_normalization`` config argument
and by setting the ``centering`` and ``scaling`` key to ``None``. For example, if the mass input is named "precipitation"
and the target feature is named "discharge", this would look like this:

.. code-block:: yaml

    custom_normalization:
        precipitation:
            centering: None
            scaling: None
        discharge:
            centering: None
            scaling: None

All inputs specified by the ``dynamic_inputs`` config argument are used as auxiliary inputs, as are (possibly embedded)
static inputs (e.g. catchment attributes).
The config argument ``head`` is ignored for this model and the model prediction is always computed as the sum over the 
outgoing mass (excluding the trash cell output).

MTS-LSTM
^^^^^^^^
:py:class:`neuralhydrology.modelzoo.mtslstm.MTSLSTM` is a newly proposed model by `Gauch et al. "Rainfall--Runoff Prediction at Multiple Timescales with a Single Long Short-Term Memory Network" <https://arxiv.org/abs/2010.07921>`__.
This model allows the training on more than temporal resolution (e.g., daily and hourly inputs) and
returns multi-timescale model predictions accordingly. A more detailed tutorial will follow shortly.

ODE-LSTM
^^^^^^^^
:py:class:`neuralhydrology.modelzoo.odelstm.ODELSTM` is a PyTorch implementation of the ODE-LSTM proposed by
`Lechner and Hasani <https://arxiv.org/abs/2006.04418>`_. This model can be used with unevenly sampled inputs and can
be queried to return predictions for any arbitrary time step.

Transformer
^^^^^^^^^^^
:py:class:`neuralhydrology.modelzoo.transformer.Transformer` is the encoding portion of a standard transformer network with self-attention. 
This uses the standard PyTorch TransformerEncoder implementation. All features (``x_d``, ``x_s``, ``x_one_hot``) are concatenated and passed 
to the network at each time step. Unless the number of inputs is divisible by the number of transformer heads (``transformer_nheads``), it is
necessary to use an embedding network that guarantees this. To achieve this, use ``statics/dynamics_embedding``, so the static/dynamic
inputs will be passed through embedding networks before being concatenated. The embedding network will then map the static and dynamic features
to size ``statics/dynamics_embedding['hiddens'][-1]``, so the total embedding size will be the sum of these values.
Instead of a decoder, this model uses a standard head (e.g., linear). 
The model requires the following hyperparameters specified in the config file: 

* ``transformer_positional_encoding_type``: choices to "sum" or "concatenate" positional encoding to other model inputs.
* ``transformer_positional_dropout``: fraction of dropout applied to the positional encoding.
* ``transformer_nheads``: number of self-attention heads.
* ``transformer_dim_feedforward``: dimension of the feedforward networks between self-attention heads.
* ``transformer_dropout``: dropout in the feedforward networks between self-attention heads.
* ``transformer_nlayers``: number of stacked self-attention + feedforward layers.


Implementing a new model
^^^^^^^^^^^^^^^^^^^^^^^^
The listing below shows the skeleton of a template model you can use to start implementing your own model.
Once you have implemented your model, make sure to modify :py:func:`neuralhydrology.modelzoo.__init__.get_model`.
Furthermore, make sure to select a *unique* model abbreviation that will be used to specify the model in the config.yml
files.

.. code-block:: python

    from typing import Dict

    import torch

    from neuralhydrology.modelzoo.basemodel import BaseModel


    class TemplateModel(BaseModel):

        def __init__(self, cfg: dict):
            """Initialize the model

            Each model receives as only input the config dictionary. From this, the entire model has to be implemented in
            this class (with potential use of other modules, such as FC from fc.py). So this class will get the model inputs
            and has to return the predictions.

            Each Model inherits from the BaseModel, which implements some universal functionality. The basemodel also
            defines the output_size, which can be used here as a given attribute (self.output_size)

            Parameters
            ----------
            cfg : dict
                Configuration of the run, read from the config file with some additional keys (such as number of basins).
            """
            super(TemplateModel, self).__init__(cfg=cfg)

            ###########################
            # Create model parts here #
            ###########################

        def forward(self, data: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
            """Forward pass through the model

            By convention, each forward pass has to accept a dict of input tensors. Usually, this dict contains 'x_d' and,
            possibly, x_s and x_one_hot. If x_d and x_s are available at multiple frequencies, the keys 'x_d' and 'x_s'
            have frequency suffixes such as 'x_d_1H' for hourly data.
            Furthermore, by definition, each model has to return a dict containing the network predictions in 'y_hat',
            potentially in addition to other keys. LSTM-based models should stick to the convention to return (at least)
            the following three tensors: y_hat, h_n, c_n (or, in the multi-frequency case, y_hat_1H, y_hat_1D, etc.).

            Parameters
            ----------
            data : Dict[str, torch.Tensor]
                 Dictionary with tensors
                    - x_d of shape [batch size, sequence length, features] containing the dynamic input data.
                    - x_s of shape [batch size, features] containing static input features. These are the concatenation
                        of what is defined in the config under static_attributes and evolving_attributes. In case not a single
                        camels attribute or static input feature is defined in the config, x_s will not be present.
                    - x_one_hot of shape [batch size, number of basins] containing the one hot encoding of the basins.
                        In case 'use_basin_id_encoding' is set to False in the config, x_one_hot will not be present.
                    Note: If the input data are available at multiple frequencies (via use_frequencies), each input tensor
                        will have a suffix "_{freq}" indicating the tensor's frequency.

            Returns
            -------
            The network prediction has to be returned under the dictionary key 'y_hat' (or, if multiple frequencies are
            predicted, 'y_hat_{freq}'. Furthermore, make sure to return predictions for each time step, even if you want
            to train sequence-to-one. Which predictions are used for training the network is controlled in the train_epoch()
            function in neuralhydrology/training/basetrainer.py. Other return values should be the hidden states as 'h_n' and cell
            states 'c_n'. Further return values are possible.
            """
            ###############################
            # Implement forward pass here #
            ###############################
            pass
Quick Start
============

Prerequisites
-------------
As a first step you need a Python environment with all required dependencies. The recommended way is to use Mini-/Anaconda
and to create a new environment using one of our predefined environment files in `environments/ <https://github.com/neuralhydrology/neuralhydrology/tree/master/environments>`__.
Make sure to select the correct file, depending on your system.

If you don't have a CUDA-capable GPU, use:

.. code-block::

    conda env create -f environments/environment_cpu.yml

If you do have a CUDA-capable GPU, use either ``environment_cuda10_2.yml`` or ``environment_cuda11_3.yml``, depending on your hardware.

If you prefer to not use Mini-/Anaconda, make sure you have a Python environment with Python >= 3.7 with all packages installed that are listed in 
these environment files. 
The next steps should be executed from within this Python environment.

Installation
------------
There are two ways how you can install NeuralHydrology: Editable or non-editable.
If all you want to do is run experiments with existing datasets and existing models, you can use the non-editable
installation. To install the latest release from PyPI:

.. code-block::

    pip install neuralhydrology

To install the package directly from the current master branch of this repository, including any changes that are not yet part of a release, run:

.. code-block::

    pip install git+https://github.com/neuralhydrology/neuralhydrology.git

If you want to try implementing your own models or datasets, you'll need an editable installation.
For this, start by downloading or cloning the repository to your local machine.
If you use git, you can run:

.. code-block::

    git clone https://github.com/neuralhydrology/neuralhydrology.git

If you don't know git, you can also download the code from `here <https://github.com/neuralhydrology/neuralhydrology/zipball/master>`__ and extract the zip-file.

After you cloned or downloaded the zip-file, you'll end up with a directory called "neuralhydrology" (or "neuralhydrology-master").
Next, we'll go to that directory and install a local, editable copy of the package:

.. code-block::

    cd neuralhydrology
    pip install -e .

The installation procedure (both the editable and the non-editable version) adds the package to your Python environment and installs three bash scripts:
`nh-run`, `nh-run-scheduler` and `nh-results-ensemble`. For details, see below.

Data
----
Training and evaluating models requires a dataset.
If you're unsure where to start, a common dataset is CAMELS US, available at
`CAMELS US (NCAR) <https://ral.ucar.edu/solutions/products/camels>`_.
This dataset is used in all of our tutorials and we have a `dedicated tutorial <../tutorials/data-prerequisites.nblink>`_ with download instructions that you might want to look at.


Training a model
----------------
To train a model, prepare a configuration file, then run::

    nh-run train --config-file /path/to/config.yml

If you want to train multiple models, you can make use of the ``nh-run-scheduler`` command.
Place all configs in a folder, then run::

    nh-run-scheduler train --directory /path/to/config_dir/ --runs-per-gpu X --gpu-ids Y

With X, you can specify how many models should be trained on parallel on a single GPU.
With Y, you can specify which GPUs to use for training (use the id as specified in ``nvidia-smi``).


Evaluating a model
------------------
To evaluate a trained model on the test set, run::

    nh-run evaluate --run-dir /path/to/run_dir/

If the optional argument ``--epoch N`` (where N is the epoch to evaluate) is not specified,
the weights of the last epoch are used.

To evaluate all runs in a specific directory you can, similarly to training, run::

    nh-run-scheduler evaluate --directory /path/to/config_dir/ --runs-per-gpu X --gpu-ids Y


To merge the predictons of a number of runs (stored in ``$DIR1``, ...) into one averaged ensemble,
use the ``nh-results-ensemble`` script::

    nh-results-ensemble --run-dirs $DIR1 $DIR2 ... --output-dir /path/to/output/directory --metrics NSE MSE ...

``--metrics`` specifies which metrics will be calculated for the averaged predictions.
Configuration Arguments
=======================

This page provides a list of possible configuration arguments.
Check out the file `examples/config.yml.example <https://github.com/neuralhydrology/neuralhydrology/blob/master/examples/config.yml.example>`__ for an example of how a config file could look like.

General experiment configurations
---------------------------------

-  ``experiment_name``: Defines the name of your experiment that will be
   used as a folder name (+ date-time string), as well as the name in
   TensorBoard.

-  ``run_dir``: Full or relative path to where the run directory is
   stored (if empty runs are stored in ${current\_working\_dir}/runs/)

-  ``train_basin_file``: Full or relative path to a text file containing
   the training basins (use data set basin id, one id per line).
-  ``validation_basin_file``: Full or relative path to a text file
   containing the validation basins (use data set basin id, one id per
   line).
-  ``test_basin_file``: Full or relative path to a text file containing
   the training basins (use data set basin id, one id per line).

-  ``train_start_date``: Start date of the training period (first day of
   discharge) in the format ``DD/MM/YYYY``. Can also be a list of dates
   to specify multiple training periods. If a list is specified, ``train_end_date``
   must also be a list of equal length. Corresponding pairs of start and
   end date denote the different periods.
-  ``train_end_date``: End date of the training period (last day of
   discharge) in the format ``DD/MM/YYYY``. Can also be a list of dates.
   If a list is specified, also ``train_start_date`` must be a list with
   an equal length.
-  ``validation_start_date``: Start date of the validation period (first
   day of discharge) in the format ``DD/MM/YYYY``. Can also be 
   a list of dates (similar to train period specifications).
-  ``validation_end_date``: End date of the validation period (last day
   of discharge) in the format ``DD/MM/YYYY``. Can also be 
   a list of dates (similar to train period specifications).
-  ``test_start_date``: Start date of the test period (first day of
   discharge) in the format ``DD/MM/YYYY``. Can also be 
   a list of dates (similar to train period specifications).
-  ``test_end_date``: End date of the validation period (last day of
   discharge) in the format ``DD/MM/YYYY``. Can also be 
   a list of dates (similar to train period specifications).
-  ``per_basin_train_periods_file``: Alternatively to specifying a global
   train period for all basins (using ``train_start_date`` and ``train_end_date``)
   it is also possible to use individual periods for each basin. For that to work
   you have to create a dictionary with one key per basin id. For each basin id,
   the dictionary contains two keys ``start_dates`` and ``end_dates``, which
   contain a list of pandas TimeStamps that indicate the start and end dates
   of the train periods. ``start_dates`` and ``end_dates`` have to be a list,
   even in case of a single period per basin (see example below). Then use the
   pickle library, to store this dictionary to disk and use the path to this
   pickle file as the value for this config argument.

.. code-block::

   import pandas as pd

   dates = {
        'basin_a': {
            'start_dates': [pd.to_datetime('01/01/1980')],
            'end_dates': [pd.to_datetime('31/12/1999')]
        },
        'basin_b': {
            'start_dates': [pd.to_datetime('01/01/1980'), pd.to_datetime('01/01/2000')],
            'end_dates': [pd.to_datetime('31/12/1990'), pd.to_datetime('01/01/2005')]
        }
    }

-  ``per_basin_validation_periods_file``: Same as ``per_basin_train_periods_file``
   but indicating individual periods that are used as validation periods.
-  ``per_basin_test_periods_file``: Same as ``per_basin_train_periods_file``
   but indicating individual periods that are used as test periods.

-  ``seed``: Fixed random seed. If empty, a random seed is generated for
   this run.

-  ``device``: Which device to use in format of ``cuda:0``, ``cuda:1``,
   etc, for GPUs or ``cpu``

Validation settings
-------------------

-  ``validate_every``: Integer that specifies in which interval a
   validation is performed. If empty, no validation is done during
   training.

-  ``validate_n_random_basins``: Integer that specifies how many random
   basins to use per validation. Values larger *n_basins* are clipped
   to *n_basins*.

-  ``metrics``: List of metrics to calculate during validation/testing.
   See
   :py:mod:`neuralhydrology.evaluation.metrics`
   for a list of available metrics. Can also be a dictionary of lists,
   where each key-value pair corresponds to one target variable and
   the metrics that should be computed for this target. Like this,
   it is possible to compute different metrics per target during 
   validation and/or evaluation after the training.

-  ``save_validation_results``: True/False, if True, stores the
   validation results to disk as a pickle file. Otherwise they are only
   used for TensorBoard

General model configuration
---------------------------

-  ``model``: Defines the model class, i.e. the core of the model, that will be used. Names
   have to match the values in `this
   function <https://github.com/neuralhydrology/neuralhydrology/blob/master/neuralhydrology/modelzoo/__init__.py#L17>`__,
   e.g., [``cudalstm``, ``ealstm``, ``mtslstm``]

-  ``head``: The prediction head that is used on top of the output of
   the model class. Currently supported are ``regression``, ``gmm``, ``cmal``, and ``umal``.
   Make sure to pass the necessary options depending on your
   choice of the head (see below).

-  ``hidden_size``: Hidden size of the model class. In the case of an
   LSTM, this reflects the number of LSTM states.

-  ``initial_forget_bias``: Initial value of the forget gate bias.

-  ``output_dropout``: Dropout applied to the output of the LSTM

Regression head
~~~~~~~~~~~~~~~
Can be ignored if ``head != 'regression'``

-  ``output_activation``: Which activation to use on the output
   neuron(s) of the linear layer. Currently supported are ``linear``,
   ``relu``, ``softplus``. If empty, ``linear`` is used.
-  ``mc_dropout``: True/False. Wheter Monte-Carlo dropout is used to 
   sample during inference. 
   
GMM head
~~~~~~~~
Can be ignored if ``head != 'gmm'``

-  ``n_distributions``: The number of distributions used for the GMM head. 
-  ``n_samples``: Number of samples generated  (per time-step) from GMM. 
-  ``negative_sample_handling``: How to account for negative samples. 
   Possible values are ``none`` for doing nothing, ``clip`` for clipping 
   the values at zero, and ``truncate`` for resampling values that
   were drawn below zero. If the last option is chosen, the additional 
   argument ``negative_sample_max_retries`` controls how often the values 
   are resampled. 
-  ``negative_sample_max_retries``: The number of repeated samples for the 
   ``truncate`` option of the ``negative_sample_max_retries`` argument.
-  ``mc_dropout``: True/False. Whether Monte-Carlo dropout is used to 
   sample during inference. 

CMAL head
~~~~~~~~~
Can be ignored if ``head != 'cmal'``

-  ``n_distributions``: The number of distributions used for the CMAL head. 
-  ``n_samples``: Number of samples generated  (per time-step) from CMAL. 
-  ``negative_sample_handling``: Approach for handling negative sampling. 
   Possible values are ``none`` for doing nothing, ``clip`` for clipping 
   the values at zero, and ``truncate`` for resampling values that
   were drawn below zero. If the last option is chosen, the additional 
   argument ``negative_sample_max_retries`` controls how often the values 
   are resampled. 
-  ``negative_sample_max_retries``: The number of repeated samples for the 
   ``truncate`` option of the ``negative_sample_max_retries`` argument.
-  ``mc_dropout``: True/False. Whether Monte-Carlo dropout is used to 
   sample during inference.    


UMAL head
~~~~~~~~~
Can be ignored if ``head != 'umal'``

-  ``n_taus``: The number of taus sampled to approximate the 
   uncountable distributions.
-  ``umal_extend_batch``: True/False. Whether the batches should be 
   extended ``n_taus`` times, to account for a specific approximation 
   density already during the training.
-  ``tau_down`` The lower sampling bound of asymmetry parameter (should be 
   above 0, below 1 and smaller than ``tau_up``).
-  ``tau_up`` The upper sampling bound of asymmetry parameter (should be 
   above 0, below 1 and larger than ``tau_down``).   
-  ``n_samples``: Number of samples generated  (per time-step) from UMAL. 
-  ``negative_sample_handling``: Approach for handling negative sampling. 
   Possible values are ``none`` for doing nothing, ``clip`` for clipping 
   the values at zero, and ``truncate`` for resampling values that
   were drawn below zero. If the last option is chosen, the additional 
   argument ``negative_sample_max_retries`` controls how often the values 
   are resampled. 
-  ``negative_sample_max_retries``: The number of repeated samples for the 
   ``truncate`` option of the ``negative_sample_max_retries`` argument.
-  ``mc_dropout``: True/False. Whether Monte-Carlo dropout is used to 
   sample during inference. 

Multi-timescale training settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These are used if ``model == mtslstm``.

-  ``transfer_mtslstm_states``: Specifies if and how hidden and cell
   states are transferred from lower to higher frequencies. This
   configuration should be a dictionary with keys ``h`` (hidden state)
   and ``c`` (cell state). Possible values are
   ``[None, linear, identity]``. If ``transfer_mtslstm_states`` is not
   provided or empty, the default is linear transfer.

-  ``shared_mtslstm``: If False, will use a distinct LSTM with
   individual weights for each timescale. If True, will use a single
   LSTM for all timescales and use one-hot-encoding to identify the
   current input timescale. In both cases, ``transfer_mtslstm_states``
   can be used to configure hidden and cell state transfer.

Transformer settings
~~~~~~~~~~~~~~~~~~~~

These are used if ``model == transformer``.

-  ``transformer_nlayers``: Number of multi-head self-attention layers in the 
   transformer encoder.
-  ``transformer_positional_encoding_type``: Choices are ``[sum, concatenate]``.
   Used to change the way that the positional encoding is used in transformer
   embedding layer. `sum` means that the positional encoding is added to the values
   of the inputs for that layer, while `concatenate` means that the embedding is concatenated
   as additional input features.
-  ``transformer_dim_feedforward``: Dimension of dense layers used between
   self-attention layers in transformer encoder.
-  ``transformer_positional_dropout``: Dropout applied only to the positional
   encoding before using in transformer encoder.
-  ``transformer_dropout``: Dropout used in transformer encoder layers.
-  ``transformer_nhead``: Number of parallel transformer heads.

ODE-LSTM settings
~~~~~~~~~~~~~~~~~

These are used if ``model == odelstm``.

-  ``ode_method``: Method to use to solve the ODE. One of
   ``[euler, rk4, heun]``.

-  ``ode_num_unfolds``: Number of iterations to break each ODE solving
   step into.

-  ``ode_random_freq_lower_bound``: Lowest frequency that will be used
   to randomly aggregate the first slice of the input sequence. See the
   documentation of the ODELSTM class for more details on the frequency
   randomization.

MC-LSTM settings
~~~~~~~~~~~~~~~~

These are used if ``model == mclstm``.

-  ``mass_inputs``: List of features that are used as mass input in the MC-LSTM model, i.e. whose quantity is conserved
   over time. Currently, the MC-LSTM configuration implemented here only supports a single mass input. Make sure to
   exclude this feature from the default normalization (see :ref:`MC-LSTM <MC-LSTM>` description).

Embedding network settings
--------------------------

These settings define fully connected networks that are used in various places, such as the embedding network
for static or dynamic features in the single-frequency models or as an optional extended input gate network in
the EA-LSTM model. For multi-timescale models, these settings can be ignored.

- ``statics_embedding``: None (default) or a dict that defines the embedding network for static inputs.
   The dictionary can have the following keys:

   - ``type`` (default 'fc'): Type of the embedding net. Currently, only 'fc' for fully-connected net is supported.
   - ``hiddens``: List of integers that define the number of neurons per layer in the fully connected network.
     The last number is the number of output neurons. Must have at least length one.
   - ``activation`` (default 'tanh'): activation function of the network. Supported values are 'tanh', 'sigmoid', 'linear'.
     The activation function is not applied to the output neurons, which always have a linear activation function.
     An activation function for the output neurons has to be applied in the main model class.
   - ``dropout`` (default 0.0): Dropout rate applied to the embedding network.

  Note that for EA-LSTM, there will always be an additional linear layer that maps to the EA-LSTM's hidden size. This
  means that the the embedding layer output size does not have to be equal to ``hidden_size``.

- ``dynamics_embedding``: None (default) or a dict that defines the embedding network for dynamic inputs. See ``statics_embedding``
  for a description of the dictionary structure.

Training settings
-----------------

-  ``optimizer``: Specify which optimizer to use. Currently supported
   is Adam (standard). New optimizers can be added
   :py:func:`here <neuralhydrology.training.get_optimizer>`.

-  ``loss``: Which loss to use. Currently supported are ``MSE``,
   ``NSE``, ``RMSE``, ``GMMLoss``, ``CMALLoss``, and ``UMALLoss``. New 
   losses can be added :py:mod:`here <neuralhydrology.training.loss>`.

- ``allow_subsequent_nan_losses``: Define a number of training steps for
   which a loss value of ``NaN`` is ignored and no error is raised but 
   instead the training loop proceeds to the next iteration step.

-  ``target_loss_weights``: A list of float values specifying the 
   per-target loss weight, when training on multiple targets at once. 
   Can be combined with any loss. By default, the weight of each target
   is ``1/n`` with ``n`` being the number of target variables. The order 
   of the weights corresponds to the order of the ``target_variables``.

-  ``regularization``: List of optional regularization terms. Currently
   supported is ``tie_frequencies``, which couples the predictions of
   all frequencies via an MSE term. New regularizations can be added
   :py:mod:`here <neuralhydrology.training.regularization>`.

-  ``learning_rate``: Learning rate. Can be either a single number (for
   a constant learning rate) or a dictionary. If it is a dictionary, the
   keys must be integer that reflect the epochs at which the learning
   rate is changed to the corresponding value. The key ``0`` defines the
   initial learning rate.

-  ``batch_size``: Mini-batch size used for training.

-  ``epochs``: Number of training epochs

-  ``use_frequencies``: Defines the time step frequencies to use (daily,
   hourly, ...). Use `pandas frequency
   strings <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timeseries-offset-aliases>`__
   to define frequencies. Note: The strings need to include values,
   e.g., '1D' instead of 'D'. If used, ``predict_last_n`` and
   ``seq_length`` must be dictionaries.

-  ``no_loss_frequencies``: Subset of frequencies from
   ``use_frequencies`` that are "evaluation-only", i.e., the model will
   get input and produce output in the frequencies listed here, but they
   will not be considered in the calculation of loss and regularization
   terms.

-  ``seq_length``: Length of the input sequence. If ``use_frequencies``
   is used, this needs to be a dictionary mapping each frequency to a
   sequence length, else an int.

-  ``predict_last_n``: Defines which time steps are used to calculate
   the loss, counted backwards. Can't be larger than ``seq_length``.
   Sequence-to-one would be ``predict_last_n: 1`` and
   sequence-to-sequence (with e.g. a sequence length of 365)
   ``predict_last_n: 365``. If ``use_frequencies`` is used, this needs
   to be a dictionary mapping each frequency to a
   predict\_last\_n-value, else an int.

-  ``target_noise_std``: Defines the standard deviation of gaussian
   noise which is added to the labels during training. Set to zero or
   leave empty to *not* add noise.

-  ``clip_gradient_norm``: If a value, clips norm of gradients to that
   specific value during training. Leave empty for not clipping.

-  ``num_workers``: Number of (parallel) threads used in the data
   loader.

-  ``save_weights_every``: Interval, in which the weights of the model
   are stored to disk. ``1`` means to store the weights after each
   epoch, which is the default if not otherwise specified.
   
Finetune settings
-----------------

Ignored if ``mode != finetune``

-  ``finetune_modules``: List of model parts that will be trained
   during fine-tuning. All parts *not* listed here will not be
   updated. Check the documentation of each model to see a list
   of available module parts.

Logger settings
---------------

-  ``log_interval``: Interval at which the training loss is logged, 
   by default 10.
-  ``log_tensorboard``: True/False. If True, writes logging results into
   TensorBoard file. The default, if not specified, is True.

-  ``log_n_figures``: If a (integer) value greater than 0, saves the
   predictions as plots of that n specific (random) basins during
   validations.

-  ``save_git_diff``: If set to True and NeuralHydrology is a git repository
   with uncommitted changes, the git diff will be stored in the run directory.
   When using this option, make sure that your run and data directories are either
   not located inside the git repository, or that they are part of the ``.gitignore`` file.
   Otherwise, the git diff may become very large and use up a lot of disk space.
   To make sure everything is configured correctly, you can simply check that the
   output of ``git diff HEAD`` only contains your code changes.

Data settings
-------------

-  ``dataset``: Defines which data set will be used. Currently supported
   are ``camels_us`` (CAMELS data set by Newman et al.), ``CAMELS_GB``
   (the GB version of CAMELS by Coxon et al.), ``CAMELS_CL`` (the CL
   version of CAMELS by Alvarez-Garreton et al.), and 
   ``hourly_camels_us`` (hourly data for 516 CAMELS basins).

-  ``data_dir``: Full or relative path to the root directory of the data set.

-  ``train_data_file``: If not empty, uses the pickled file at this path
   as the training data. Can be used to not create the same data set
   multiple times, which saves disk space and time. If empty, creates
   new data set and optionally stores the data in the run directory (if
   ``save_train_data`` is True).

-  ``cache_validation_data``: True/False. If True, caches validation data 
   in memory for the time of training, which does speed up the overall
   training time. By default True, since even larger datasets are usually
   just a few GB in memory, which most modern machines can handle.

-  ``dynamic_inputs``: List of variables to use as time series inputs.
   Names must match the exact names as defined in the data set. Note: In
   case of multiple input forcing products, you have to append the
   forcing product behind each variable. E.g., 'prcp(mm/day)' of the
   daymet product is 'prcp(mm/day)_daymet'. When training on multiple
   frequencies (cf. ``use_frequencies``), it is possible to define
   dynamic inputs for each frequency individually. To do so,
   ``dynamic_inputs`` must be a dict mapping each frequency to a list of
   variables. E.g., to use precipitation from daymet for daily and from
   nldas-hourly for hourly predictions:

   ::

       dynamic_inputs:
         1D:
           - prcp(mm/day)_daymet
         1H:
           - total_precipitation_nldas_hourly

-  ``target_variables``: List of the target variable(s). Names must match
   the exact names as defined in the data set.

-  ``clip_targets_to_zero``: Optional list of target variables to clip to
   zero during the computation of metrics (e.g. useful to compute zero-clipped metric during the validation between
   training epochs. Will not affect the data that is saved to disk after evaluation. 
   That is, always the `raw` model outputs are saved in the result files. Therefore, you eventually need to 
   manually clip the targets to zero if you load the model outputs from file and want to reproduce
   the metric values.

-  ``duplicate_features``: Can be used to duplicate time series features
   (e.g., for different normalizations). Can be either a str, list or dictionary
   (mapping from strings to ints). If string, duplicates the corresponding
   feature once. If list, duplicates all features in that list once. Use
   a dictionary to specify the exact number of duplicates you like.
   To each duplicated feature, we append ``_copyN``, where `N` is counter
   starting at 1.

-  ``lagged_features``: Can be used to add a lagged copy of another
   feature to the list of available input/output features. Has to be a
   dictionary mapping from strings to int or a list of ints, where the string 
   specifies the feature name and the int(s) the number of lagged time steps. Those values
   can be positive or negative (see
   `pandas shift <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.shift.html>`__
   for details). If a list of integers is provided, only unique values are considered.
   We append ``_shiftN`` to each lagged feature, where `N` is the shift count.

-  ``custom_normalization``: Has to be a dictionary, mapping from
   time series feature names to ``centering`` and/or ``scaling``. Using
   this argument allows to overwrite the default zero mean, unit
   variance normalization per feature. Supported options for
   ``centering`` are 'None' or 'none', 'mean', 'median' and min.
   None/none sets the centering parameter to 0.0, mean to the feature
   mean, median to the feature median, and min to the feature
   minimum, respectively. Supported options for `scaling` are
   'None' or 'none', 'std', 'minmax'. None/none sets the scaling
   parameter to 1.0, std to the feature standard deviation and
   minmax to the feature max minus the feature min. The combination
   of centering: min and scaling: minmax results in min/max
   feature scaling to the range [0,1].

-  ``additional_feature_files``: Path to a pickle file (or list of paths
   for multiple files), containing a dictionary with each key
   corresponding to one basin id and the value is a date-time indexed
   pandas DataFrame. Allows the option to add any arbitrary data that is
   not included in the standard data sets. **Convention**: If a column
   is used as static input, the value to use for specific sample should
   be in same row (datetime) as the target discharge value.

-  ``evolving_attributes``: Columns of the DataFrame loaded with the
   ``additional_feature_files`` that should be used as "static" features.
   These values will be used as static inputs, but they can evolve over time.
   Convention: The value to use for a specific input sequence should be in the
   same row (datetime) as the last time step of that sequence.
   Names must match the column names in the DataFrame. Leave empty to
   not use any additional static feature.

-  ``use_basin_id_encoding``: True/False. If True, creates a
   basin-one-hot encoding as a(n) (additional) static feature vector for
   each sample.

-  ``static_attributes``: Which static attributes to use (e.g., from the static camels attributes for the CAMELS
   dataset). Leave empty if none should be used. For hydroatlas attributes, use ``hydroatlas_attributes`` instead.
   Names must match the exact names as defined in the data set.

-  ``hydroatlas_attributes``: Which HydroATLAS attributes to use. Leave
   empty if none should be used. Names must match the exact names as
   defined in the data set.

CAMELS US specific
~~~~~~~~~~~~~~~~~~

Can be ignored if ``dataset not in ['camels_us', 'hourly_camels_us']``

-  ``forcings``: Can be either a string or a list of strings that
   correspond to forcing products in the camels data set. Also supports
   ``maurer_extended``, ``nldas_extended``, and (for
   ``hourly_camels_us``) ``nldas_hourly``.
Tutorials
=========

All tutorials are based on Jupyter notebooks that are hosted on GitHub. 
If you want to run the code yourself, you can find the notebooks in the `examples folder <https://github.com/neuralhydrology/neuralhydrology/tree/master/examples>`__ of the NeuralHydrology GitHub repository.

| **Data Prerequisites**
| For most of our tutorials you will need some data to train and evaluate models. In all of these examples we use the publicly available CAMELS US dataset. :doc:`This tutorial <data-prerequisites>` will guide you through the download process of the different dataset pieces and explain how the code expects the local folder structure.

| **Introduction to NeuralHydrology**
| If you're new to the NeuralHydrology package, :doc:`this tutorial <introduction>` is the place to get started. It walks you through the basic command-line and API usage patterns, and you get to train and evaluate your first model.

| **Adding a New Model: Gated Recurrent Unit (GRU)**
| Once you know the basics, you might want to add your own model. Using the `GRU <https://en.wikipedia.org/wiki/Gated_recurrent_unit>`__ model as an example, :doc:`this tutorial <adding-gru>` shows how and where to add models in the NeuralHydrology codebase.

| **Adding a New Dataset: CAMELS-CL**
| Always using the United States CAMELS dataset is getting boring? :doc:`This tutorial <add-dataset>` shows you how to add a new dataset: The Chilean version of CAMELS.

| **Multi-Timescale Prediction**
| In one of our `papers <https://arxiv.org/abs/2010.07921>`__, we introduced Multi-Timescale LSTMs that can predict at multiple timescales simultaneously. If you need predictions at sub-daily granularity or you want to generate daily and hourly predictions (or any other timescale), :doc:`this tutorial <multi-timescale>` explains how to get there.

| **Inspecting the internals of LSTMs**
| Model interpretability is an ongoing research topic. We showed in previous publications (e.g. `this one <https://arxiv.org/abs/1903.07903>`__) that LSTM internals can be linked to physical processes. In :doc:`this tutorial <inspect-lstm>`, we show how to extract those model internals with our library.

| **Finetuning models**
| A common way to increase model performance with deep learning models is called finetuning. Here, first a model is trained on a large and diverse dataset, before second, the model is finetuned to the actual problem of interest. In :doc:`this tutorial <finetuning>`, we show how you can perform finetuning with our library.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   data-prerequisites
   introduction
   adding-gru
   add-dataset
   multi-timescale
   inspect-lstm
   finetuning
loss
====

.. automodule:: neuralhydrology.training.loss
   :members:
   :undoc-members:
   :show-inheritance:
GenericDataset
==============

.. automodule:: neuralhydrology.datasetzoo.genericdataset
   :members:
   :undoc-members:
   :show-inheritance:
nh.training
===========

.. automodule:: neuralhydrology.training
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 4

   neuralhydrology.training.basetrainer
   neuralhydrology.training.logger
   neuralhydrology.training.loss
   neuralhydrology.training.regularization
   neuralhydrology.training.train
nh.modelzoo
===========

.. automodule:: neuralhydrology.modelzoo
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 4

   neuralhydrology.modelzoo.basemodel
   neuralhydrology.modelzoo.cudalstm
   neuralhydrology.modelzoo.customlstm
   neuralhydrology.modelzoo.ealstm
   neuralhydrology.modelzoo.embcudalstm
   neuralhydrology.modelzoo.fc
   neuralhydrology.modelzoo.gru
   neuralhydrology.modelzoo.head
   neuralhydrology.modelzoo.inputlayer
   neuralhydrology.modelzoo.mclstm
   neuralhydrology.modelzoo.mtslstm
   neuralhydrology.modelzoo.odelstm
   neuralhydrology.modelzoo.template
   neuralhydrology.modelzoo.transformer
Transformer
===========

.. automodule:: neuralhydrology.modelzoo.transformer
   :members:
   :undoc-members:
   :show-inheritance:
ODELSTM
=======

.. automodule:: neuralhydrology.modelzoo.odelstm
   :members:
   :undoc-members:
   :show-inheritance:
BaseDataset
===========

.. automodule:: neuralhydrology.datasetzoo.basedataset
   :members:
   :undoc-members:
   :show-inheritance:EmbCudaLSTM
===========

.. automodule:: neuralhydrology.modelzoo.embcudalstm
   :members:
   :undoc-members:
   :show-inheritance:
CustomLSTM
==========

.. automodule:: neuralhydrology.modelzoo.customlstm
   :members:
   :undoc-members:
   :show-inheritance:
utils
=====

.. automodule:: neuralhydrology.datautils.utils
   :members:
   :undoc-members:
   :show-inheritance:
EALSTM
======

.. automodule:: neuralhydrology.modelzoo.ealstm
   :members:
   :undoc-members:
   :show-inheritance:
BaseModel
=========

.. automodule:: neuralhydrology.modelzoo.basemodel
   :members:
   :undoc-members:
   :show-inheritance:
MC-LSTM
=======

.. automodule:: neuralhydrology.modelzoo.mclstm
   :members:
   :undoc-members:
   :show-inheritance:
MTSLSTM
=======

.. automodule:: neuralhydrology.modelzoo.mtslstm
   :members:
   :undoc-members:
   :show-inheritance:
CamelsGB
========

.. automodule:: neuralhydrology.datasetzoo.camelsgb
   :members:
   :undoc-members:
   :show-inheritance:dischargeinput
==============

.. automodule:: neuralhydrology.datautils.dischargeinput
   :members:
   :undoc-members:
   :show-inheritance:
Tester
======

.. automodule:: neuralhydrology.evaluation.tester
   :members:
   :undoc-members:
   :show-inheritance:
nh\_results\_ensemble module
============================

.. automodule:: neuralhydrology.utils.nh_results_ensemble
   :members:
   :undoc-members:
   :show-inheritance:
nh\_run
=======

.. automodule:: neuralhydrology.nh_run
   :members:
   :undoc-members:
   :show-inheritance:
plots
=====

.. automodule:: neuralhydrology.evaluation.plots
   :members:
   :undoc-members:
   :show-inheritance:
nh.datasetzoo
=============

.. automodule:: neuralhydrology.datasetzoo
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 4

   neuralhydrology.datasetzoo.basedataset
   neuralhydrology.datasetzoo.camelscl
   neuralhydrology.datasetzoo.camelsgb
   neuralhydrology.datasetzoo.camelsus
   neuralhydrology.datasetzoo.genericdataset
   neuralhydrology.datasetzoo.hourlycamelsus
   neuralhydrology.datasetzoo.template
TemplateDataset
===============

.. automodule:: neuralhydrology.datasetzoo.template
   :members:
   :undoc-members:
   :show-inheritance:pet
===

.. automodule:: neuralhydrology.datautils.pet
   :members:
   :undoc-members:
   :show-inheritance:
GRU
===

.. automodule:: neuralhydrology.modelzoo.gru
   :members:
   :undoc-members:
   :show-inheritance:
NeuralHydrology API
===================

.. automodule:: neuralhydrology
   :members:
   :undoc-members:
   :show-inheritance:


.. toctree::
   :maxdepth: 4

   neuralhydrology.datasetzoo
   neuralhydrology.datautils
   neuralhydrology.evaluation
   neuralhydrology.modelzoo
   neuralhydrology.training
   neuralhydrology.utils

Main Entry Points
-----------------

.. toctree::
   :maxdepth: 4

   neuralhydrology.nh_run
   neuralhydrology.nh_run_scheduler
   neuralhydrology.utils.nh_results_ensemble
Logger
======

.. automodule:: neuralhydrology.training.logger
   :members:
   :undoc-members:
   :show-inheritance:
nh.utils
========

.. automodule:: neuralhydrology.utils
   :members:
   :undoc-members:
   :show-inheritance:


.. toctree::
   :maxdepth: 4

   neuralhydrology.utils.config
   neuralhydrology.utils.configutils
   neuralhydrology.utils.errors
regularization
==============

.. automodule:: neuralhydrology.training.regularization
   :members:
   :undoc-members:
   :show-inheritance:
climateindices
==============

.. automodule:: neuralhydrology.datautils.climateindices
   :members:
   :undoc-members:
   :show-inheritance:
nh.datautils
============

.. automodule:: neuralhydrology.datautils
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 4

   neuralhydrology.datautils.climateindices
   neuralhydrology.datautils.dischargeinput
   neuralhydrology.datautils.pet
   neuralhydrology.datautils.utils
InputLayer
==========

.. automodule:: neuralhydrology.modelzoo.inputlayer
   :members:
   :undoc-members:
   :show-inheritance:
head
====

.. automodule:: neuralhydrology.modelzoo.head
   :members:
   :undoc-members:
   :show-inheritance:
errors
======

.. automodule:: neuralhydrology.utils.errors
   :members:
   :undoc-members:
   :show-inheritance:
config
======

.. automodule:: neuralhydrology.utils.config
   :members:
   :undoc-members:
   :show-inheritance:
configutils
===========

.. automodule:: neuralhydrology.utils.configutils
   :members:
   :undoc-members:
   :show-inheritance:
signatures
==========

.. automodule:: neuralhydrology.evaluation.signatures
   :members:
   :undoc-members:
   :show-inheritance:
CamelsUS
========

.. automodule:: neuralhydrology.datasetzoo.camelsus
   :members:
   :undoc-members:
   :show-inheritance:
evaluate
========

.. automodule:: neuralhydrology.evaluation.evaluate
   :members:
   :undoc-members:
   :show-inheritance:
FC
==

.. automodule:: neuralhydrology.modelzoo.fc
   :members:
   :undoc-members:
   :show-inheritance:
nh.evaluation
=============

.. automodule:: neuralhydrology.evaluation
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 4

   neuralhydrology.evaluation.evaluate
   neuralhydrology.evaluation.metrics
   neuralhydrology.evaluation.plots
   neuralhydrology.evaluation.signatures
   neuralhydrology.evaluation.tester
BaseTrainer
===========

.. automodule:: neuralhydrology.training.basetrainer
   :members:
   :undoc-members:
   :show-inheritance:
CudaLSTM
========

.. automodule:: neuralhydrology.modelzoo.cudalstm
   :members:
   :undoc-members:
   :show-inheritance:
TemplateModel
=============

.. automodule:: neuralhydrology.modelzoo.template
   :members:
   :undoc-members:
   :show-inheritance:
CamelsCL
========

.. automodule:: neuralhydrology.datasetzoo.camelscl
   :members:
   :undoc-members:
   :show-inheritance:
train
=====

.. automodule:: neuralhydrology.training.train
   :members:
   :undoc-members:
   :show-inheritance:
HourlyCamelsUS
==============

.. automodule:: neuralhydrology.datasetzoo.hourlycamelsus
   :members:
   :undoc-members:
   :show-inheritance:nh\_run\_scheduler
==================

.. automodule:: neuralhydrology.nh_run_scheduler
   :members:
   :undoc-members:
   :show-inheritance:
metrics
=======

.. automodule:: neuralhydrology.evaluation.metrics
   :members:
   :undoc-members:
   :show-inheritance:
