<div align="center"><img src="https://raw.githubusercontent.com/optuna/optuna/master/docs/image/optuna-logo.png" width="800"/></div>

# Optuna: A hyperparameter optimization framework

[![Python](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8%20%7C%203.9-blue)](https://www.python.org)
[![pypi](https://img.shields.io/pypi/v/optuna.svg)](https://pypi.python.org/pypi/optuna)
[![conda](https://img.shields.io/conda/vn/conda-forge/optuna.svg)](https://anaconda.org/conda-forge/optuna)
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/optuna/optuna)
[![CircleCI](https://circleci.com/gh/optuna/optuna.svg?style=svg)](https://circleci.com/gh/optuna/optuna)
[![Read the Docs](https://readthedocs.org/projects/optuna/badge/?version=stable)](https://optuna.readthedocs.io/en/stable/)
[![Codecov](https://codecov.io/gh/optuna/optuna/branch/master/graph/badge.svg)](https://codecov.io/gh/optuna/optuna/branch/master)
[![Gitter chat](https://badges.gitter.im/optuna/gitter.svg)](https://gitter.im/optuna/optuna)

[**Website**](https://optuna.org/)
| [**Docs**](https://optuna.readthedocs.io/en/stable/)
| [**Install Guide**](https://optuna.readthedocs.io/en/stable/installation.html)
| [**Tutorial**](https://optuna.readthedocs.io/en/stable/tutorial/index.html)

*Optuna* is an automatic hyperparameter optimization software framework, particularly designed
for machine learning. It features an imperative, *define-by-run* style user API. Thanks to our
*define-by-run* API, the code written with Optuna enjoys high modularity, and the user of
Optuna can dynamically construct the search spaces for the hyperparameters.

## News

- **2021-12-06** First alpha version of Optuna 3.0 is released! Early adopters may want to upgrade and provide feedback for a smoother transition to the coming major release. Try `pip install optuna==3.0.0a0`.

- **2021-10-11**  Optuna 3.0 Roadmap published for review. Please take a look at the [planned improvements to Optuna](https://github.com/optuna/optuna/wiki/Optuna-V3-Roadmap), and share your feedback in the github issues. PR contributions also welcome!

- **2021-07-14** Please take a few minutes to fill in this survey, and let us know how you use Optuna now and what improvements you'd like.🤔
All questions optional. 🙇‍♂️
https://forms.gle/mCAttqxVg5oUifKV8

## Key Features

Optuna has modern functionalities as follows:

- [Lightweight, versatile, and platform agnostic architecture](https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/001_first.html)
  - Handle a wide variety of tasks with a simple installation that has few requirements.
- [Pythonic search spaces](https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/002_configurations.html)
  - Define search spaces using familiar Python syntax including conditionals and loops.
- [Efficient optimization algorithms](https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/003_efficient_optimization_algorithms.html)
  - Adopt state-of-the-art algorithms for sampling hyperparameters and efficiently pruning unpromising trials.
- [Easy parallelization](https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/004_distributed.html)
  - Scale studies to tens or hundreds or workers with little or no changes to the code.
- [Quick visualization](https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/005_visualization.html)
  - Inspect optimization histories from a variety of plotting functions.


## Basic Concepts

We use the terms *study* and *trial* as follows:

- Study: optimization based on an objective function
- Trial: a single execution of the objective function

Please refer to sample code below. The goal of a *study* is to find out the optimal set of
hyperparameter values (e.g., `classifier` and `svm_c`) through multiple *trials* (e.g.,
`n_trials=100`). Optuna is a framework designed for the automation and the acceleration of the
optimization *studies*.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/optuna/optuna-examples/blob/main/quickstart.ipynb)

```python
import ...

# Define an objective function to be minimized.
def objective(trial):

    # Invoke suggest methods of a Trial object to generate hyperparameters.
    regressor_name = trial.suggest_categorical('classifier', ['SVR', 'RandomForest'])
    if regressor_name == 'SVR':
        svr_c = trial.suggest_float('svr_c', 1e-10, 1e10, log=True)
        regressor_obj = sklearn.svm.SVR(C=svr_c)
    else:
        rf_max_depth = trial.suggest_int('rf_max_depth', 2, 32)
        regressor_obj = sklearn.ensemble.RandomForestRegressor(max_depth=rf_max_depth)

    X, y = sklearn.datasets.fetch_california_housing(return_X_y=True)
    X_train, X_val, y_train, y_val = sklearn.model_selection.train_test_split(X, y, random_state=0)

    regressor_obj.fit(X_train, y_train)
    y_pred = regressor_obj.predict(X_val)

    error = sklearn.metrics.mean_squared_error(y_val, y_pred)

    return error  # An objective value linked with the Trial object.

study = optuna.create_study()  # Create a new study.
study.optimize(objective, n_trials=100)  # Invoke optimization of the objective function.
```

## Examples

Examples can be found in [optuna/optuna-examples](https://github.com/optuna/optuna-examples).

## Integrations

[Integrations modules](https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/003_efficient_optimization_algorithms.html#integration-modules-for-pruning), which allow pruning, or early stopping, of unpromising trials are available for the following libraries:

* [AllenNLP](https://github.com/optuna/optuna-examples/tree/main/allennlp)
* [Catalyst](https://github.com/optuna/optuna-examples/tree/main/pytorch/catalyst_simple.py)
* [Catboost](https://github.com/optuna/optuna-examples/tree/main/catboost/catboost_simple.py)
* [Chainer](https://github.com/optuna/optuna-examples/tree/main/chainer/chainer_integration.py)
* FastAI ([V1](https://github.com/optuna/optuna-examples/tree/main/fastai/fastaiv1_simple.py), [V2](https://github.com/optuna/optuna-examples/tree/main/fastai/fastaiv2_simple.py))
* [Keras](https://github.com/optuna/optuna-examples/tree/main/keras/keras_integration.py)
* [LightGBM](https://github.com/optuna/optuna-examples/tree/main/lightgbm/lightgbm_integration.py)
* [MXNet](https://github.com/optuna/optuna-examples/tree/main/mxnet/mxnet_integration.py)
* [PyTorch](https://github.com/optuna/optuna-examples/tree/main/pytorch/pytorch_simple.py)
* [PyTorch Ignite](https://github.com/optuna/optuna-examples/tree/main/pytorch/pytorch_ignite_simple.py)
* [PyTorch Lightning](https://github.com/optuna/optuna-examples/tree/main/pytorch/pytorch_lightning_simple.py)
* [TensorFlow](https://github.com/optuna/optuna-examples/tree/main/tensorflow/tensorflow_estimator_integration.py)
* [tf.keras](https://github.com/optuna/optuna-examples/tree/main/tfkeras/tfkeras_integration.py)
* [XGBoost](https://github.com/optuna/optuna-examples/tree/main/xgboost/xgboost_integration.py)


## Web Dashboard (experimental)

The new Web dashboard is under the development at [optuna-dashboard](https://github.com/optuna/optuna-dashboard).
It is still experimental, but much better in many regards.
Feature requests and bug reports welcome!

| Manage studies | Visualize with interactive graphs |
| -------------- | --------------------------------- |
| ![manage-studies](https://user-images.githubusercontent.com/5564044/97099702-4107be80-16cf-11eb-9d97-f5ceec98ce52.gif) | ![optuna-realtime-graph](https://user-images.githubusercontent.com/5564044/97099797-66e19300-16d0-11eb-826c-6977e3941fb0.gif) |

Install `optuna-dashboard` via pip:

```
$ pip install optuna-dashboard
$ optuna-dashboard sqlite:///db.sqlite3
...
Listening on http://localhost:8080/
Hit Ctrl-C to quit.
```

## Installation

Optuna is available at [the Python Package Index](https://pypi.org/project/optuna/) and on [Anaconda Cloud](https://anaconda.org/conda-forge/optuna).

```bash
# PyPI
$ pip install optuna
```

```bash
# Anaconda Cloud
$ conda install -c conda-forge optuna
```

Optuna supports Python 3.6 or newer.

Also, we also provide Optuna docker images on [DockerHub](https://hub.docker.com/r/optuna/optuna).

## Communication

- [GitHub Issues] for bug reports, feature requests and questions.
- [Gitter] for interactive chat with developers.
- [Stack Overflow] for questions.

[GitHub issues]: https://github.com/optuna/optuna/issues
[Gitter]: https://gitter.im/optuna/optuna
[Stack Overflow]: https://stackoverflow.com/questions/tagged/optuna


## Contribution

Any contributions to Optuna are more than welcome!

If you are new to Optuna, please check the [good first issues](https://github.com/optuna/optuna/labels/good%20first%20issue). They are relatively simple, well-defined and are often good starting points for you to get familiar with the contribution workflow and other developers.

If you already have contributed to Optuna, we recommend the other [contribution-welcome issues](https://github.com/optuna/optuna/labels/contribution-welcome).

For general guidelines how to contribute to the project, take a look at [CONTRIBUTING.md](./CONTRIBUTING.md).


## Reference

Takuya Akiba, Shotaro Sano, Toshihiko Yanase, Takeru Ohta, and Masanori Koyama. 2019.
Optuna: A Next-generation Hyperparameter Optimization Framework. In KDD ([arXiv](https://arxiv.org/abs/1907.10902)).
# Optuna Code of Conduct

Optuna follows the [NumFOCUS Code of Conduct][homepage] available at https://numfocus.org/code-of-conduct.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at optuna@preferred.jp. 

[homepage]: https://numfocus.org/
# Contribution Guidelines

It’s an honor to have you on board!

We are proud of this project and have been working to make it great since day one.
We believe you will love it, and we know there’s room for improvement.
We want to
- implement features that make what you want to do possible and/or easy.
- write more tutorials and [examples](https://github.com/optuna/optuna-examples) that help you get familiar with Optuna.
- make issues and pull requests on GitHub fruitful.
- have more conversations and discussions on Gitter.

We need your help and everything about Optuna you have in your mind pushes this project forward.
Join Us!

If you feel like giving a hand, here are some ways:
- Implement a feature
    - If you have some cool idea, please open an issue first to discuss design to make your idea in a better shape.
- Send a patch
    - Dirty your hands by tackling [issues with `contribution-welcome` label](https://github.com/optuna/optuna/issues?q=is%3Aissue+is%3Aopen+label%3Acontribution-welcome)
- Report a bug
    - If you find a bug, please report it! Your reports are important.
- Fix/Improve documentation
    - Documentation gets outdated easily and can always be better, so feel free to fix and improve
- Let us and the Optuna community know your ideas and thoughts.
    - __Contribution to Optuna includes not only sending pull requests, but also writing down your comments on issues and pull requests by others, and joining conversations/discussions on [Gitter](https://gitter.im/optuna/optuna).__
    - Also, sharing how you enjoy Optuna is a huge contribution! If you write a blog, let us know about it!

If you write code, we have some conventions as follows.

- [Guidelines](#guidelines)
- [Unit Tests](#unit-tests)
- [Continuous Integration and Local Verification](#continuous-integration-and-local-verification)
- [Creating a Pull Request](#creating-a-pull-request)

## Guidelines

### Setup Optuna

First of all, fork Optuna on GitHub.
You can learn about fork in the official [documentation](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo).

After forking, download and install Optuna on your computer.

```bash
git clone git@github.com:YOUR_NAME/optuna.git
cd optuna
pip install -e .
```

### Checking the Format, Coding Style, and Type Hints

Code is formatted with [black](https://github.com/psf/black),
and docstrings are formatted with [blackdoc](https://github.com/keewis/blackdoc).
Coding style is checked with [flake8](http://flake8.pycqa.org) and [isort](https://pycqa.github.io/isort/),
and additional conventions are described in the [Wiki](https://github.com/optuna/optuna/wiki/Coding-Style-Conventions).
Type hints, [PEP484](https://www.python.org/dev/peps/pep-0484/), are checked with [mypy](http://mypy-lang.org/).

You can check the format, coding style, and type hints at the same time just by executing a script `formats.sh`.
If your environment is missing some dependencies such as black, blackdoc, flake8, isort or mypy,
you will be asked to install them.

You can also check them using [tox](https://tox.readthedocs.io/en/latest/) like below.

```
$ pip install tox
$ tox -e flake8 -e black -e blackdoc -e isort -e mypy
```

If you catch format errors, you can automatically fix them by auto-formatters.

```bash
# Install auto-formatters.
$ pip install ".[checking]"

$ ./formats.sh
```

### Documentation

When adding a new feature to the framework, you also need to document it in the reference.
The documentation source is stored under the [docs](./docs) directory and written in [reStructuredText format](http://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html).

To build the documentation, you need to run:

```bash
pip install -e ".[document]"
```
Note that the above command might try to install PyTorch without CUDA to your environment even if your environment has CUDA version already.

Then you can build the documentation in HTML format locally:

```bash
cd docs
make html
```

HTML files are generated under `build/html` directory. Open `index.html` with the browser and see
if it is rendered as expected.

Optuna's tutorial is built with [Sphinx-Gallery](https://sphinx-gallery.github.io/stable/index.html) and
some other requirements like [LightGBM](https://github.com/microsoft/LightGBM) and [PyTorch](https://pytorch.org) meaning that
all .py files in `tutorial` directory are run during the documentation build if there's no build cache.
Whether you edit any tutorial or not doesn't matter.

To avoid having to run the tutorials, you may download executed tutorial artifacts named "tutorial" from our CI (see the capture below) and put them in `docs/build` before
extracting the files in the zip to `docs/source/tutorial` directory.

![image](https://user-images.githubusercontent.com/16191443/107472296-0b211400-6bb2-11eb-9203-e2c42ce499ad.png)

**Writing a Tutorial**
Tutorials are part of Optuna’s documentation.
Optuna depends on Sphinx to build the documentation HTML files from the corresponding reStructuredText (`.rst`) files in the docs/source directory,
but as you may notice, [Tutorial directory](https://github.com/optuna/optuna/tree/master/tutorial) does not have any `.rst` files. Instead, it has a bunch of Python (`.py`) files.
We have [Sphinx Gallery](https://sphinx-gallery.github.io/stable/index.html) that executes those `.py` files and generates `.rst` files with standard outputs from them and corresponding Jupyter Notebook (`.ipynb`) files.
These generated `.rst` and `.ipynb` files are written to the docs/source/tutorial directory.
The output directory (docs/source/tutorial) and source (tutorial) directory are configured in [`sphinx_gallery_conf ` of docs/source/conf.py](https://github.com/optuna/optuna/blob/2e14273cab87f13edeb9d804a43bd63c44703cb5/docs/source/conf.py#L189-L199). These generated `.rst` files are handled by Sphinx like the other `.rst` files. The generated `.ipynb` files are hosted on Optuna’s documentation page and downloadable (check [Optuna tutorial](https://optuna.readthedocs.io/en/stable/tutorial/index.html)).

The order of contents on [tutorial top page](https://optuna.readthedocs.io/en/stable/tutorial/index.html) is determined by two keys: one is the subdirectory name of tutorial and the other is the filename (note that there are some alternatives as documented in [Sphinx Gallery - sorting](https://sphinx-gallery.github.io/stable/gen_modules/sphinx_gallery.sorting.html?highlight=filenamesortkey), but we chose this key in https://github.com/optuna/optuna/blob/2e14273cab87f13edeb9d804a43bd63c44703cb5/docs/source/conf.py#L196).
Optuna’s tutorial directory has two directories: (1) [10_key_features](https://github.com/optuna/optuna/tree/master/tutorial/10_key_features), which is meant to be aligned with and explain the key features listed on [README.md](https://github.com/optuna/optuna#key-features) and (2) [20_recipes](https://github.com/optuna/optuna/tree/master/tutorial/20_recipes), whose contents showcase how to use Optuna features conveniently.
When adding new content to the Optuna tutorials, place it in `20_recipes` and its file name should conform to the other names, for example, `777_cool_feature.py`.
In general, please number the prefix for your file consecutively with the last number. However, this is not mandatory and if you think your content deserves the smaller number (the order of recipes does not have a specific meaning, but in general, order could convey the priority order to readers), feel free to propose the renumbering in your PR.

You may want to refer to the Sphinx Gallery for the syntax of `.py` files processed by Sphinx Gallery.
Two specific conventions and limitations for Optuna tutorials:
1. 99 #s for block separation as in https://github.com/optuna/optuna/blob/2e14273cab87f13edeb9d804a43bd63c44703cb5/tutorial/10_key_features/001_first.py#L19
2. Execution time of the new content needs to be less than three minutes. This limitation derives from Read The Docs. If your content runs some hyperparameter optimization, set the `timeout` to 180 or less. You can check this limitation on [Read the Docs - Build Process](https://docs.readthedocs.io/en/stable/builds.html).


## Unit Tests

When adding a new feature or fixing a bug, you also need to write sufficient test code.
We use [pytest](https://pytest.org/) as the testing framework and
unit tests are stored under the [tests directory](./tests).

Please install some required packages at first.
```bash
# Install required packages to test all modules without visualization and integration modules.
pip install ".[tests]"

# Install required packages to test all modules including visualization and integration modules.
pip install ".[testing]" -f https://download.pytorch.org/whl/torch_stable.html
```

You can run your tests as follows:

```bash
# Run all the unit tests.
pytest

# Run all the unit tests defined in the specified test file.
pytest tests/${TARGET_TEST_FILE_NAME}

# Run the unit test function with the specified name defined in the specified test file.
pytest tests/${TARGET_TEST_FILE_NAME} -k ${TARGET_TEST_FUNCTION_NAME}
```

## Continuous Integration and Local Verification

Optuna repository uses GitHub Actions and CircleCI.

Currently, we are migrating to GitHub Actions but still we use CircleCI for testing `document`
because it makes it much easier to check built documentation.

## Creating a Pull Request

When you are ready to create a pull request, please try to keep the following in mind:

### Title

The title of your pull request should

- briefly describe and reflect the changes
- wrap any code with backticks
- not end with a period

*The title will be directly visible in the release notes.*

#### Example

Introduces Tree-structured Parzen Estimator to `optuna.samplers`

### Description

The description of your pull request should

- describe the motivation
- describe the changes
- if still work-in-progress, describe remaining tasks
<!-- Thank you for creating a pull request! In general, we merge your pull requests after they get two or more approvals. -->

## Motivation
<!-- Describe your motivation why you will submit this PR. This is useful for reviewers to understand the context of PR. -->

## Description of the changes
<!-- Describe the changes in this PR. -->
# Benchmarks for Optuna

Interested in measuring Optuna's performance? 
You are very perceptive. 
Under this directory, you will find scripts that we have prepared to measure Optuna's performance.

In this document, we explain how we measure the performance of Optuna using the scripts in this directory.
The contents of this document are organized as follows.
- [Performance Benchmarks with `kurobako`](#performance-benchmarks-with-kurobako)

## Performance Benchmarks with `kurobako`

We measure the performance of black-box optimization algorithms in Optuna with 
[`kurobako`](https://github.com/optuna/kurobako) using `benchmarks/run_kurobako.py`.
You can manually run this script on the GitHub Actions if you have a write access on the repository.
Or, you can locally execute the `benchmarks/run_kurobako.py`.
We explain both of method here.

### How to Run on the GitHub Actions

You need a write access on the repository.
Please run the following steps in your own forks.
Note that you should pull the latest master branch of [Optuna](https://github.com/optuna/optuna) since the workflow YAML file must be placed in the default branch of the repository.

1. Open the GitHub page of your forked Optuna repository.
2. Click the `Actions` below the repository name.
![image](https://user-images.githubusercontent.com/38826298/145764682-0c4a31aa-f865-4293-a3c7-2ca6be5baa03.png)

3. In the left sidebar, click the `Performance Benchmarks with kurobako`.
4. Above the list of workflow runs, select `Run workflow`.
![image](https://user-images.githubusercontent.com/38826298/145764692-a30a74c0-5ebe-4010-a7cd-4ebcdbb24679.png)

5. Use the `Branch` dropdown to select the workflow's branch. The default is `master`. 
And, type the input parameters: 
`Sampler List`, `Sampler Arguments List`, `Pruner List`, and `Pruner Arguments List`.
6. Click `Run workflow`.
![image](https://user-images.githubusercontent.com/38826298/145764702-771d9a6f-8c7d-40d5-a912-1485a1d7dcfa.png)
7. After finishing the workflow, you can download the report and plot from `Artifacts`.
![image](https://user-images.githubusercontent.com/38826298/145802414-e29ca0ba-80fd-488a-af02-c33e9b4d5e3b.png)
The report looks like as follows.
It includes the version information of environments, the solvers (pairs of the sampler and the pruner in Optuna) and problems, the best objective value, AUC, elapsed time, and so on. 
![image](https://user-images.githubusercontent.com/38826298/146860092-74da99c6-15b6-4da4-baef-0457af1d7171.png)
The plot looks like as follows.
It represents the optimization history plot of the optimization.
The title is the name of the problem.
The legends represents the specified pair of the sampler and the pruner.
The history is averaged over the specified `n_runs` studies with the errorbar.
The horizontal axis represents the budget (`#budgets * #epochs = \sum_{for each trial) (#consumed epochs in the trial)`).
The vertical axis represents the objective value.
![image](https://user-images.githubusercontent.com/38826298/146860370-853174c7-afc5-4f67-8143-61f22d2c8f6c.png)


Note that the default run time of a GitHub Actions workflow job is limited to 6 hours. 
Depending on the sampler and number of studies you specify, it may exceed the 6-hour limit and fail.
See the [official document](https://docs.github.com/ja/actions/learn-github-actions/usage-limits-billing-and-administration) for more details.

### How to Run Locally

You can run the script of `benchmarks/run_kurobako.py` directly.
This section explains how to locally run it.

First, you need to install `kurobako` and its Python helper.
To install `kurobako`, see https://github.com/optuna/kurobako#installation for more details.
In addition, please run `pip install kurobako` to install the Python helper.
You need to install `gnuplot` for visualization with `kurobako`.
You can install `gnuplot` by package managers such as `apt` (for Ubuntu) or `brew` (for macOS).

Second, you need to download the dataset for `kurobako`.
Run the followings in the dataset directory.
```bash
# Download hyperparameter optimization (HPO) dataset
% wget http://ml4aad.org/wp-content/uploads/2019/01/fcnet_tabular_benchmarks.tar.gz
% tar xf fcnet_tabular_benchmarks.tar.gz

# Download neural architecture search (NAS) dataset
# The `kurobako` command should be available.
% curl -L $(kurobako dataset nasbench url) -o nasbench_full.tfrecord
% kurobako dataset nasbench convert nasbench_full.tfrecord nasbench_full.bin
```

Finally, you can run the script of `benchmarks/run_kurobako.py`.
```bash
% python benchmarks/run_kurobako.py \
          --path-to-kurobako "" \ # If the `kurobako` command is available.
          --name "performance-benchmarks" \
          --n-runs 10 \
          --n-jobs 10 \
          --sampler-list "RandomSampler TPESampler" \
          --sampler-kwargs-list "{} {}" \
          --pruner-list "NopPruner" \
          --pruner-kwargs-list "{}" \
          --seed 0 \
          --data-dir "." \
          --out-dir "out"
```
Please see `benchmarks/run_kurobako.py` to check the arguments and those default values.
Optuna Examples
================

This page contains a list of example codes written with Optuna. The example files are in [optuna/optuna-examples](https://github.com/optuna/optuna-examples/).

### Simple Black-box Optimization

* [Quadratic function](https://github.com/optuna/optuna-examples/blob/main/quadratic_simple.py)

### Examples with ML Libraries

* [AllenNLP](https://github.com/optuna/optuna-examples/blob/main/allennlp/allennlp_simple.py)
* [AllenNLP (Jsonnet)](https://github.com/optuna/optuna-examples/blob/main/allennlp/allennlp_jsonnet.py)
* [Catalyst](https://github.com/optuna/optuna-examples/blob/main/pytorch/catalyst_simple.py)
* [CatBoost](https://github.com/optuna/optuna-examples/blob/main/catboost/catboost_simple.py)
* [Chainer](https://github.com/optuna/optuna-examples/blob/main/chainer/chainer_simple.py)
* [ChainerMN](https://github.com/optuna/optuna-examples/blob/main/chainer/chainermn_simple.py)
* [Dask-ML](https://github.com/optuna/optuna-examples/blob/main/dask_ml/dask_ml_simple.py)
* [FastAI V1](https://github.com/optuna/optuna-examples/blob/main/fastai/fastaiv1_simple.py)
* [FastAI V2](https://github.com/optuna/optuna-examples/blob/main/fastai/fastaiv2_simple.py)
* [Haiku](https://github.com/optuna/optuna-examples/blob/main/haiku/haiku_simple.py)
* [Gluon](https://github.com/optuna/optuna-examples/blob/main/mxnet/gluon_simple.py)
* [Keras](https://github.com/optuna/optuna-examples/blob/main/keras/keras_simple.py)
* [LightGBM](https://github.com/optuna/optuna-examples/blob/main/lightgbm/lightgbm_simple.py)
* [LightGBM Tuner](https://github.com/optuna/optuna-examples/blob/main/lightgbm/lightgbm_tuner_simple.py)
* [MXNet](https://github.com/optuna/optuna-examples/blob/main/mxnet/mxnet_simple.py)
* [PyTorch](https://github.com/optuna/optuna-examples/blob/main/pytorch/pytorch_simple.py)
* [PyTorch Ignite](https://github.com/optuna/optuna-examples/blob/main/pytorch/pytorch_ignite_simple.py)
* [PyTorch Lightning](https://github.com/optuna/optuna-examples/blob/main/pytorch/pytorch_lightning_simple.py)
* [RAPIDS](https://github.com/optuna/optuna-examples/blob/main/rapids_simple.py)
* [Scikit-learn](https://github.com/optuna/optuna-examples/blob/main/sklearn/sklearn_simple.py)
* [Scikit-learn OptunaSearchCV](https://github.com/optuna/optuna-examples/blob/main/sklearn/sklearn_optuna_search_cv_simple.py)
* [Scikit-image](https://github.com/optuna/optuna-examples/blob/main/skimage/skimage_lbp_simple.py)
* [SKORCH](https://github.com/optuna/optuna-examples/blob/main/pytorch/skorch_simple.py)
* [Tensorflow](https://github.com/optuna/optuna-examples/blob/main/tensorflow/tensorflow_estimator_simple.py)
* [Tensorflow (eager)](https://github.com/optuna/optuna-examples/blob/main/tensorflow/tensorflow_eager_simple.py)
* [XGBoost](https://github.com/optuna/optuna-examples/blob/main/xgboost/xgboost_simple.py)

### An example where an objective function uses additional arguments

The following example demonstrates how to implement an objective function that uses additional arguments other than `trial`.
* [Scikit-learn (callable class version)](https://github.com/optuna/optuna-examples/tree/main/sklearn/sklearn_additional_args.py)

### Examples of Pruning

The following example demonstrates how to implement pruning logic with Optuna.

* [Simple pruning (scikit-learn)](https://github.com/optuna/optuna-examples/blob/main/simple_pruning.py)

In addition, integration modules are available for the following libraries, providing simpler interfaces to utilize pruning.

* [Pruning with Catalyst integration module](https://github.com/optuna/optuna-examples/blob/main/pytorch/catalyst_simple.py)
* [Pruning with Chainer integration module](https://github.com/optuna/optuna-examples/blob/main/chainer/chainer_integration.py)
* [Pruning with ChainerMN integration module](https://github.com/optuna/optuna-examples/blob/main/chainer/chainermn_integration.py)
* [Pruning with FastAI V1 integration module](https://github.com/optuna/optuna-examples/blob/main/fastai/fastaiv1_simple.py)
* [Pruning with FastAI V2 integration module](https://github.com/optuna/optuna-examples/blob/main/fastai/fastaiv2_simple.py)
* [Pruning with Keras integration module](https://github.com/optuna/optuna-examples/blob/main/keras/keras_integration.py)
* [Pruning with LightGBM integration module](https://github.com/optuna/optuna-examples/blob/main/lightgbm/lightgbm_integration.py)
* [Pruning with MXNet integration module](https://github.com/optuna/optuna-examples/blob/main/mxnet/mxnet_integration.py)
* [Pruning with PyTorch integration module](https://github.com/optuna/optuna-examples/blob/main/pytorch/pytorch_simple.py)
* [Pruning with PyTorch Ignite integration module](https://github.com/optuna/optuna-examples/blob/main/pytorch/pytorch_ignite_simple.py)
* [Pruning with PyTorch Lightning integration module](https://github.com/optuna/optuna-examples/blob/main/pytorch/pytorch_lightning_simple.py)
* [Pruning with Tensorflow integration module](https://github.com/optuna/optuna-examples/blob/main/tensorflow/tensorflow_estimator_integration.py)
* [Pruning with XGBoost integration module](https://github.com/optuna/optuna-examples/blob/main/xgboost/xgboost_integration.py)
* [Pruning with XGBoost integration module (cross validation, XGBoost.cv)](https://github.com/optuna/optuna-examples/blob/main/xgboost/xgboost_cv_integration.py)

### Examples of Samplers

* [Warm Starting CMA-ES](https://github.com/optuna/optuna-examples/blob/main/samplers/warm_starting_cma.py)
* [User-defined SimulatedAnnealingSampler](https://github.com/optuna/optuna-examples/blob/main/samplers/simulated_annealing_sampler.py)

### Examples of Multi-Objective Optimization

* [Optimization with BoTorch](https://github.com/optuna/optuna-examples/blob/main/multi_objective/botorch_simple.py)
* [Optimization of MLP with PyTorch](https://github.com/optuna/optuna-examples/blob/main/multi_objective/pytorch_simple.py)

### Examples of Visualization

* [Visualizing study](https://colab.research.google.com/github/optuna/optuna-examples/blob/main/visualization/plot_study.ipynb)

### An example to enqueue trials with given parameter values

* [Enqueuing trials with given parameters](https://github.com/optuna/optuna-examples/blob/main/enqueue_trial.py)

### Examples of MLflow

* [Tracking optimization process with MLflow](https://github.com/optuna/optuna-examples/blob/main/mlflow/keras_mlflow.py)

### Examples of Weights & Biases

* [Tracking optimization process with Weights & Biases](https://github.com/optuna/optuna-examples/blob/main/wandb/wandb_simple.py)

### Examples of Hydra

* [Optimization with Hydra](https://github.com/optuna/optuna-examples/blob/main/hydra/simple.py)

### Examples of Distributed Optimization

* [Optimizing on Kubernetes](https://github.com/optuna/optuna-examples/blob/main/kubernetes/README.md)
* [Optimizing with Ray's joblib backend](https://github.com/optuna/optuna-examples/blob/main/ray/ray_joblib.py)

### Examples of Reinforcement Learning

* [Optimization of Hyperparameters for Stable-Baslines Agent](https://github.com/optuna/optuna-examples/blob/main/rl/sb3_simple.py)

### External projects using Optuna

* [Allegro Trains](https://github.com/allegroai/trains)
* [BBO-Rietveld: Automated crystal structure refinement](https://github.com/quantumbeam/BBO-Rietveld)
* [Catalyst](https://github.com/catalyst-team/catalyst)
* [CuPy](https://github.com/cupy/cupy)
* [Hydra's Optuna Sweeper plugin](https://hydra.cc/docs/next/plugins/optuna_sweeper/)
* [Mozilla Voice STT](https://github.com/mozilla/DeepSpeech)
* [neptune.ai](https://neptune.ai)
* [OptGBM: A scikit-learn compatible LightGBM estimator with Optuna](https://github.com/Y-oHr-N/OptGBM)
* [PyKEEN](https://github.com/pykeen/pykeen)
* [RL Baselines Zoo](https://github.com/DLR-RM/rl-baselines3-zoo)

PRs to add additional projects welcome in [optuna-examples](https://github.com/optuna/optuna-examples)!
