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
Tutorial
========

If you are new to Optuna or want a general introduction, we highly recommend the below video.

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/P6NwZVl8ttc" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    <br />
    <br />
    <br />
.. _recipes:

Recipes
-------

Showcases the recipes that might help you using Optuna with comfort.
.. _key_features:

Key Features
------------

Showcases Optuna's `Key Features <https://github.com/optuna/optuna/blob/master/README.md#key-features>`_.
:orphan:

Privacy Policy
==============

Google Analytics
----------------

To collect information about how visitors use our website and to improve our services, we are using Google Analytics on this website. You can find out more about how Google Analytics works and about how information is collected on the Google Analytics terms of services and on Google's privacy policy.

- Google Analytics Terms of Service: http://www.google.com/analytics/terms/us.html
- Google Privacy Policy: https://policies.google.com/privacy?hl=en
- Google Analytics Opt-out Add-on: https://tools.google.com/dlpage/gaoptout?hl=en
Installation
============

Optuna supports Python 3.6 or newer.

We recommend to install Optuna via pip:

.. code-block:: bash

    $ pip install optuna

You can also install the development version of Optuna from master branch of Git repository:

.. code-block:: bash

    $ pip install git+https://github.com/optuna/optuna.git

You can also install Optuna via conda:

.. code-block:: bash

    $ conda install -c conda-forge optuna
FAQ
===

.. contents::
    :local:

Can I use Optuna with X? (where X is your favorite ML library)
--------------------------------------------------------------

Optuna is compatible with most ML libraries, and it's easy to use Optuna with those.
Please refer to `examples <https://github.com/optuna/optuna-examples/>`_.


.. _objective-func-additional-args:

How to define objective functions that have own arguments?
----------------------------------------------------------

There are two ways to realize it.

First, callable classes can be used for that purpose as follows:

.. code-block:: python

    import optuna


    class Objective(object):
        def __init__(self, min_x, max_x):
            # Hold this implementation specific arguments as the fields of the class.
            self.min_x = min_x
            self.max_x = max_x

        def __call__(self, trial):
            # Calculate an objective value by using the extra arguments.
            x = trial.suggest_float("x", self.min_x, self.max_x)
            return (x - 2) ** 2


    # Execute an optimization by using an `Objective` instance.
    study = optuna.create_study()
    study.optimize(Objective(-100, 100), n_trials=100)


Second, you can use ``lambda`` or ``functools.partial`` for creating functions (closures) that hold extra arguments.
Below is an example that uses ``lambda``:

.. code-block:: python

    import optuna

    # Objective function that takes three arguments.
    def objective(trial, min_x, max_x):
        x = trial.suggest_float("x", min_x, max_x)
        return (x - 2) ** 2


    # Extra arguments.
    min_x = -100
    max_x = 100

    # Execute an optimization by using the above objective function wrapped by `lambda`.
    study = optuna.create_study()
    study.optimize(lambda trial: objective(trial, min_x, max_x), n_trials=100)

Please also refer to `sklearn_addtitional_args.py <https://github.com/optuna/optuna-examples/tree/main/sklearn/sklearn_additional_args.py>`_ example,
which reuses the dataset instead of loading it in each trial execution.


Can I use Optuna without remote RDB servers?
--------------------------------------------

Yes, it's possible.

In the simplest form, Optuna works with in-memory storage:

.. code-block:: python

    study = optuna.create_study()
    study.optimize(objective)


If you want to save and resume studies,  it's handy to use SQLite as the local storage:

.. code-block:: python

    study = optuna.create_study(study_name="foo_study", storage="sqlite:///example.db")
    study.optimize(objective)  # The state of `study` will be persisted to the local SQLite file.

Please see :ref:`rdb` for more details.


How can I save and resume studies?
----------------------------------------------------

There are two ways of persisting studies, which depend if you are using
in-memory storage (default) or remote databases (RDB). In-memory studies can be
saved and loaded like usual Python objects using ``pickle`` or ``joblib``. For
example, using ``joblib``:

.. code-block:: python

    study = optuna.create_study()
    joblib.dump(study, "study.pkl")

And to resume the study:

.. code-block:: python

    study = joblib.load("study.pkl")
    print("Best trial until now:")
    print(" Value: ", study.best_trial.value)
    print(" Params: ")
    for key, value in study.best_trial.params.items():
        print(f"    {key}: {value}")

Note that Optuna does not support saving/reloading across different Optuna
versions with ``pickle``. To save/reload a study across different Optuna versions,
please use RDBs and `upgrade storage schema <reference/cli.html#storage-upgrade>`_
if necessary. If you are using RDBs, see :ref:`rdb` for more details.

How to suppress log messages of Optuna?
---------------------------------------

By default, Optuna shows log messages at the ``optuna.logging.INFO`` level.
You can change logging levels by using  :func:`optuna.logging.set_verbosity`.

For instance, you can stop showing each trial result as follows:

.. code-block:: python

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    study = optuna.create_study()
    study.optimize(objective)
    # Logs like '[I 2020-07-21 13:41:45,627] Trial 0 finished with value:...' are disabled.


Please refer to :class:`optuna.logging` for further details.


How to save machine learning models trained in objective functions?
-------------------------------------------------------------------

Optuna saves hyperparameter values with its corresponding objective value to storage,
but it discards intermediate objects such as machine learning models and neural network weights.
To save models or weights, please use features of the machine learning library you used.

We recommend saving :obj:`optuna.trial.Trial.number` with a model in order to identify its corresponding trial.
For example, you can save SVM models trained in the objective function as follows:

.. code-block:: python

    def objective(trial):
        svc_c = trial.suggest_float("svc_c", 1e-10, 1e10, log=True)
        clf = sklearn.svm.SVC(C=svc_c)
        clf.fit(X_train, y_train)

        # Save a trained model to a file.
        with open("{}.pickle".format(trial.number), "wb") as fout:
            pickle.dump(clf, fout)
        return 1.0 - accuracy_score(y_valid, clf.predict(X_valid))


    study = optuna.create_study()
    study.optimize(objective, n_trials=100)

    # Load the best model.
    with open("{}.pickle".format(study.best_trial.number), "rb") as fin:
        best_clf = pickle.load(fin)
    print(accuracy_score(y_valid, best_clf.predict(X_valid)))


How can I obtain reproducible optimization results?
---------------------------------------------------

To make the parameters suggested by Optuna reproducible, you can specify a fixed random seed via ``seed`` argument of :class:`~optuna.samplers.RandomSampler` or :class:`~optuna.samplers.TPESampler` as follows:

.. code-block:: python

    sampler = TPESampler(seed=10)  # Make the sampler behave in a deterministic way.
    study = optuna.create_study(sampler=sampler)
    study.optimize(objective)

However, there are two caveats.

First, when optimizing a study in distributed or parallel mode, there is inherent non-determinism.
Thus it is very difficult to reproduce the same results in such condition.
We recommend executing optimization of a study sequentially if you would like to reproduce the result.

Second, if your objective function behaves in a non-deterministic way (i.e., it does not return the same value even if the same parameters were suggested), you cannot reproduce an optimization.
To deal with this problem, please set an option (e.g., random seed) to make the behavior deterministic if your optimization target (e.g., an ML library) provides it.


How are exceptions from trials handled?
---------------------------------------

Trials that raise exceptions without catching them will be treated as failures, i.e. with the :obj:`~optuna.trial.TrialState.FAIL` status.

By default, all exceptions except :class:`~optuna.exceptions.TrialPruned` raised in objective functions are propagated to the caller of :func:`~optuna.study.Study.optimize`.
In other words, studies are aborted when such exceptions are raised.
It might be desirable to continue a study with the remaining trials.
To do so, you can specify in :func:`~optuna.study.Study.optimize` which exception types to catch using the ``catch`` argument.
Exceptions of these types are caught inside the study and will not propagate further.

You can find the failed trials in log messages.

.. code-block:: sh

    [W 2018-12-07 16:38:36,889] Setting status of trial#0 as TrialState.FAIL because of \
    the following error: ValueError('A sample error in objective.')

You can also find the failed trials by checking the trial states as follows:

.. code-block:: python

    study.trials_dataframe()

.. csv-table::

    number,state,value,...,params,system_attrs
    0,TrialState.FAIL,,...,0,Setting status of trial#0 as TrialState.FAIL because of the following error: ValueError('A test error in objective.')
    1,TrialState.COMPLETE,1269,...,1,

.. seealso::

    The ``catch`` argument in :func:`~optuna.study.Study.optimize`.


How are NaNs returned by trials handled?
----------------------------------------

Trials that return :obj:`NaN` (``float('nan')``) are treated as failures, but they will not abort studies.

Trials which return :obj:`NaN` are shown as follows:

.. code-block:: sh

    [W 2018-12-07 16:41:59,000] Setting status of trial#2 as TrialState.FAIL because the \
    objective function returned nan.


What happens when I dynamically alter a search space?
-----------------------------------------------------

Since parameters search spaces are specified in each call to the suggestion API, e.g.
:func:`~optuna.trial.Trial.suggest_float` and :func:`~optuna.trial.Trial.suggest_int`,
it is possible to, in a single study, alter the range by sampling parameters from different search
spaces in different trials.
The behavior when altered is defined by each sampler individually.

.. note::

    Discussion about the TPE sampler. https://github.com/optuna/optuna/issues/822


How can I use two GPUs for evaluating two trials simultaneously?
----------------------------------------------------------------

If your optimization target supports GPU (CUDA) acceleration and you want to specify which GPU is used, the easiest way is to set ``CUDA_VISIBLE_DEVICES`` environment variable:

.. code-block:: bash

    # On a terminal.
    #
    # Specify to use the first GPU, and run an optimization.
    $ export CUDA_VISIBLE_DEVICES=0
    $ optuna study optimize foo.py objective --study-name foo --storage sqlite:///example.db

    # On another terminal.
    #
    # Specify to use the second GPU, and run another optimization.
    $ export CUDA_VISIBLE_DEVICES=1
    $ optuna study optimize bar.py objective --study-name bar --storage sqlite:///example.db

Please refer to `CUDA C Programming Guide <https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#env-vars>`_ for further details.


How can I test my objective functions?
--------------------------------------

When you test objective functions, you may prefer fixed parameter values to sampled ones.
In that case, you can use :class:`~optuna.trial.FixedTrial`, which suggests fixed parameter values based on a given dictionary of parameters.
For instance, you can input arbitrary values of :math:`x` and :math:`y` to the objective function :math:`x + y` as follows:

.. code-block:: python

    def objective(trial):
        x = trial.suggest_float("x", -1.0, 1.0)
        y = trial.suggest_int("y", -5, 5)
        return x + y


    objective(FixedTrial({"x": 1.0, "y": -1}))  # 0.0
    objective(FixedTrial({"x": -1.0, "y": -4}))  # -5.0


Using :class:`~optuna.trial.FixedTrial`, you can write unit tests as follows:

.. code-block:: python

    # A test function of pytest
    def test_objective():
        assert 1.0 == objective(FixedTrial({"x": 1.0, "y": 0}))
        assert -1.0 == objective(FixedTrial({"x": 0.0, "y": -1}))
        assert 0.0 == objective(FixedTrial({"x": -1.0, "y": 1}))


.. _out-of-memory-gc-collect:

How do I avoid running out of memory (OOM) when optimizing studies?
-------------------------------------------------------------------

If the memory footprint increases as you run more trials, try to periodically run the garbage collector.
Specify ``gc_after_trial`` to :obj:`True` when calling :func:`~optuna.study.Study.optimize` or call :func:`gc.collect` inside a callback.

.. code-block:: python

    def objective(trial):
        x = trial.suggest_float("x", -1.0, 1.0)
        y = trial.suggest_int("y", -5, 5)
        return x + y


    study = optuna.create_study()
    study.optimize(objective, n_trials=10, gc_after_trial=True)

    # `gc_after_trial=True` is more or less identical to the following.
    study.optimize(objective, n_trials=10, callbacks=[lambda study, trial: gc.collect()])

There is a performance trade-off for running the garbage collector, which could be non-negligible depending on how fast your objective function otherwise is. Therefore, ``gc_after_trial`` is :obj:`False` by default.
Note that the above examples are similar to running the garbage collector inside the objective function, except for the fact that :func:`gc.collect` is called even when errors, including :class:`~optuna.exceptions.TrialPruned` are raised.

.. note::

    :class:`~optuna.integration.ChainerMNStudy` does currently not provide ``gc_after_trial`` nor callbacks for :func:`~optuna.integration.ChainerMNStudy.optimize`.
    When using this class, you will have to call the garbage collector inside the objective function.

How can I output a log only when the best value is updated?
-----------------------------------------------------------

Here's how to replace the logging feature of optuna with your own logging callback function.
The implemented callback can be passed to :func:`~optuna.study.Study.optimize`.
Here's an example:

.. code-block:: python

    import optuna


    # Turn off optuna log notes.
    optuna.logging.set_verbosity(optuna.logging.WARN)


    def objective(trial):
        x = trial.suggest_float("x", 0, 1)
        return x ** 2


    def logging_callback(study, frozen_trial):
        previous_best_value = study.user_attrs.get("previous_best_value", None)
        if previous_best_value != study.best_value:
            study.set_user_attr("previous_best_value", study.best_value)
            print(
                "Trial {} finished with best value: {} and parameters: {}. ".format(
                frozen_trial.number,
                frozen_trial.value,
                frozen_trial.params,
                )
            )


    study = optuna.create_study()
    study.optimize(objective, n_trials=100, callbacks=[logging_callback])

Note that this callback may show incorrect values when you try to optimize an objective function with ``n_jobs!=1``
(or other forms of distributed optimization) due to its reads and writes to storage that are prone to race conditions.

How do I suggest variables which represent the proportion, that is, are in accordance with Dirichlet distribution?
------------------------------------------------------------------------------------------------------------------

When you want to suggest :math:`n` variables which represent the proportion, that is, :math:`p[0], p[1], ..., p[n-1]` which satisfy :math:`0 \le p[k] \le 1` for any :math:`k` and :math:`p[0] + p[1] + ... + p[n-1] = 1`, try the below.
For example, these variables can be used as weights when interpolating the loss functions.
These variables are in accordance with the flat `Dirichlet distribution <https://en.wikipedia.org/wiki/Dirichlet_distribution>`_.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import optuna


    def objective(trial):
        n = 5
        x = []
        for i in range(n):
            x.append(- np.log(trial.suggest_float(f"x_{i}", 0, 1)))

        p = []
        for i in range(n):
            p.append(x[i] / sum(x))

        for i in range(n):
            trial.set_user_attr(f"p_{i}", p[i])

        return 0

    study = optuna.create_study(sampler=optuna.samplers.RandomSampler())
    study.optimize(objective, n_trials=1000)

    n = 5
    p = []
    for i in range(n):
        p.append([trial.user_attrs[f"p_{i}"] for trial in study.trials])
    axes = plt.subplots(n, n, figsize=(20, 20))[1]

    for i in range(n):
        for j in range(n):
            axes[j][i].scatter(p[i], p[j], marker=".")
            axes[j][i].set_xlim(0, 1)
            axes[j][i].set_ylim(0, 1)
            axes[j][i].set_xlabel(f"p_{i}")
            axes[j][i].set_ylabel(f"p_{j}")

    plt.savefig("sampled_ps.png")

This method is justified in the following way:
First, if we apply the transformation :math:`x = - \log (u)` to the variable :math:`u` sampled from the uniform distribution :math:`Uni(0, 1)` in the interval :math:`[0, 1]`, the variable :math:`x` will follow the exponential distribution :math:`Exp(1)` with scale parameter :math:`1`.
Furthermore, for :math:`n` variables :math:`x[0], ..., x[n-1]` that follow the exponential distribution of scale parameter :math:`1` independently, normalizing them with :math:`p[i] = x[i] / \sum_i x[i]`, the vector :math:`p` follows the Dirichlet distribution :math:`Dir(\alpha)` of scale parameter :math:`\alpha = (1, ..., 1)`.
You can verify the transformation by calculating the elements of the Jacobian.

How can I optimize a model with some constraints?
-------------------------------------------------

When you want to optimize a model with constraints, you can use the following classes, :class:`~optuna.samplers.NSGAIISampler` or :class:`~optuna.integration.BoTorchSampler`.
The following example is a benchmark of Binh and Korn function, a multi-objective optimization, with constraints using :class:`~optuna.samplers.NSGAIISampler`. This one has two constraints :math:`c_0 = (x-5)^2 + y^2 - 25 \le 0` and :math:`c_1 = -(x - 8)^2 - (y + 3)^2 + 7.7 \le 0` and finds the optimal solution satisfying these constraints.


.. code-block:: python

    import optuna


    def objective(trial):
        # Binh and Korn function with constraints.
        x = trial.suggest_float("x", -15, 30)
        y = trial.suggest_float("y", -15, 30)

        # Constraints which are considered feasible if less than or equal to zero.
        # The feasible region is basically the intersection of a circle centered at (x=5, y=0)
        # and the complement to a circle centered at (x=8, y=-3).
        c0 = (x - 5) ** 2 + y ** 2 - 25
        c1 = -((x - 8) ** 2) - (y + 3) ** 2 + 7.7

        # Store the constraints as user attributes so that they can be restored after optimization.
        trial.set_user_attr("constraint", (c0, c1))

        v0 = 4 * x ** 2 + 4 * y ** 2
        v1 = (x - 5) ** 2 + (y - 5) ** 2

        return v0, v1


    def constraints(trial):
        return trial.user_attrs["constraint"]


    sampler = optuna.samplers.NSGAIISampler(constraints_func=constraints)
    study = optuna.create_study(
        directions=["minimize", "minimize"],
        sampler=sampler,
    )
    study.optimize(objective, n_trials=32, timeout=600)

    print("Number of finished trials: ", len(study.trials))

    print("Pareto front:")

    trials = sorted(study.best_trials, key=lambda t: t.values)

    for trial in trials:
        print("  Trial#{}".format(trial.number))
        print(
            "    Values: Values={}, Constraint={}".format(
                trial.values, trial.user_attrs["constraint"][0]
            )
        )
        print("    Params: {}".format(trial.params))

If you are interested in the exmaple for :class:`~optuna.integration.BoTorchSampler`, please refer to `this sample code <https://github.com/optuna/optuna-examples/blob/main/multi_objective/botorch_simple.py>`_.


There are two kinds of constrained optimizations, one with soft constraints and the other with hard constraints.
Soft constraints do not have to be satisfied, but an objective function is penalized if they are unsatisfied. On the other hand, hard constraints must be satisfied.

Optuna is adopting the soft one and **DOES NOT** support the hard one. In other words, Optuna **DOES NOT** have built-in samplers for the hard constraints.

How can I parallelize optimization?
-----------------------------------

The variations of parallelization are in the following three cases.

1. Multi-threading parallelization with single node
2. Multi-processing parallelization with single node
3. Multi-processing parallelization with multiple nodes

1. Multi-threading parallelization with a single node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parallelization can be achieved by setting the argument ``n_jobs`` in :func:`optuna.study.Study.optimize`.
However, the python code will not be faster due to GIL because :func:`optuna.study.Study.optimize` with ``n_jobs!=1`` uses multi-threading. 

While optimizing, it will be faster in limited situations, such as waiting for other server requests or C/C++ processing with numpy, etc., but it will not be faster in other cases.

For more information about 1., see APIReference_.

.. _APIReference: https://optuna.readthedocs.io/en/stable/reference/index.html

2. Multi-processing parallelization with single node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be achieved by using file-based RDBs (such as SQLite) and client/server RDBs (such as PostgreSQL and MySQL).
However, if you are in the environment where you can not install an RDB, you can not run multi-processing parallelization with single node. When you really want to do it, please request it as a GitHub issue. If we receive a lot of requests, we may provide a solution for it.

For more information about 2., see TutorialEasyParallelization_.

.. _TutorialEasyParallelization: https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/004_distributed.html

3. Multi-processing parallelization with multiple nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be achieved by using client/server RDBs (such as PostgreSQL and MySQL).
However, if you are in the environment where you can not install a client/server RDB, you can not run multi-processing parallelization with multiple nodes.

For more information about 3., see TutorialEasyParallelization_.

Can I monitor trials and make them failed automatically when they are killed unexpectedly?
------------------------------------------------------------------------------------------

.. note::

  Heartbeat mechanism is experimental. API would change in the future.

A process running a trial could be killed unexpectedly, typically by a job scheduler in a cluster environment.
If trials are killed unexpectedly, they will be left on the storage with their states `RUNNING` until we remove them or update their state manually.
For such a case, Optuna supports monitoring trials using `heartbeat <https://en.wikipedia.org/wiki/Heartbeat_(computing)>`_ mechanism.
Using heartbeat, if a process running a trial is killed unexpectedly,
Optuna will automatically change the state of the trial that was running on that process to :obj:`~optuna.trial.TrialState.FAIL`
from :obj:`~optuna.trial.TrialState.RUNNING`.

.. code-block:: python

    import optuna

    def objective(trial):
        (Very time-consuming computation)

    # Recording heartbeats every 60 seconds.
    # Other processes' trials where more than 120 seconds have passed
    # since the last heartbeat was recorded will be automatically failed.
    storage = optuna.storages.RDBStorage(url="sqlite:///:memory:", heartbeat_interval=60, grace_period=120)
    study = optuna.create_study(storage=storage)
    study.optimize(objective, n_trials=100)

You can also execute a callback function to process the failed trial.
Optuna provides a callback to retry failed trials as :class:`~optuna.storages.RetryFailedTrialCallback`.
Note that a callback is invoked at a beginning of each trial, which means :class:`~optuna.storages.RetryFailedTrialCallback`
will retry failed trials when a new trial starts to evaluate.

.. code-block:: python

    import optuna
    from optuna.storages import RetryFailedTrialCallback

    storage = optuna.storages.RDBStorage(
        url="sqlite:///:memory:",
        heartbeat_interval=60,
        grace_period=120,
        failed_trial_callback=RetryFailedTrialCallback(max_retry=3),
    )

    study = optuna.create_study(storage=storage)
|optunalogo|

Optuna: A hyperparameter optimization framework
===============================================

*Optuna* is an automatic hyperparameter optimization software framework,
particularly designed for machine learning. It features an imperative,
*define-by-run* style user API. Thanks to our *define-by-run* API, the
code written with Optuna enjoys high modularity, and the user of Optuna
can dynamically construct the search spaces for the hyperparameters.

Key Features
------------

Optuna has modern functionalities as follows:

- :doc:`Lightweight, versatile, and platform agnostic architecture <tutorial/10_key_features/001_first>`

  - Handle a wide variety of tasks with a simple installation that has few requirements.

- :doc:`Pythonic search spaces <tutorial/10_key_features/002_configurations>`

  - Define search spaces using familiar Python syntax including conditionals and loops.

- :doc:`Efficient optimization algorithms <tutorial/10_key_features/003_efficient_optimization_algorithms>`

  - Adopt state-of-the-art algorithms for sampling hyperparameters and efficiently pruning unpromising trials.

- :doc:`Easy parallelization <tutorial/10_key_features/004_distributed>`

  - Scale studies to tens or hundreds or workers with little or no changes to the code.

- :doc:`Quick visualization <tutorial/10_key_features/005_visualization>`

  - Inspect optimization histories from a variety of plotting functions.

Basic Concepts
--------------

We use the terms *study* and *trial* as follows:

-  Study: optimization based on an objective function
-  Trial: a single execution of the objective function

Please refer to sample code below. The goal of a *study* is to find out
the optimal set of hyperparameter values (e.g., ``classifier`` and
``svm_c``) through multiple *trials* (e.g., ``n_trials=100``). Optuna is
a framework designed for the automation and the acceleration of the
optimization *studies*.

|Open in Colab|

.. code:: python

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

Communication
-------------

-  `GitHub Issues <https://github.com/optuna/optuna/issues>`__ for bug
   reports, feature requests and questions.
-  `Gitter <https://gitter.im/optuna/optuna>`__ for interactive chat
   with developers.
-  `Stack
   Overflow <https://stackoverflow.com/questions/tagged/optuna>`__ for
   questions.

Contribution
------------

Any contributions to Optuna are welcome! When you send a pull request,
please follow the `contribution guide <https://github.com/optuna/optuna/blob/master/CONTRIBUTING.md>`__.

License
-------

MIT License (see `LICENSE <https://github.com/optuna/optuna/blob/master/LICENSE>`__).

Reference
---------

Takuya Akiba, Shotaro Sano, Toshihiko Yanase, Takeru Ohta, and Masanori
Koyama. 2019. Optuna: A Next-generation Hyperparameter Optimization
Framework. In KDD (`arXiv <https://arxiv.org/abs/1907.10902>`__).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial/index
   reference/index
   faq

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |optunalogo| image:: https://raw.githubusercontent.com/optuna/optuna/master/docs/image/optuna-logo.png
  :width: 800
  :alt: OPTUNA
.. |Open in Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
  :target: http://colab.research.google.com/github/optuna/optuna-examples/blob/main/quickstart.ipynb
{% extends "!autosummary/class.rst" %}

{#
An autosummary template to exclude the class constructor (__init__)
which doesn't contain any docstring in Optuna.
#}

{% block methods %}
   {% set methods = methods | select("ne", "__init__") | list %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}

{% endblock %}
.. module:: optuna

optuna
======

The :mod:`optuna` module is primarily used as an alias for basic Optuna functionality coded in other modules. Currently, two modules are aliased: (1) from :mod:`optuna.study`, functions regarding the Study lifecycle, and (2) from :mod:`optuna.exceptions`, the TrialPruned Exception raised when a trial is pruned.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.create_study
   optuna.load_study
   optuna.delete_study
   optuna.copy_study
   optuna.get_all_study_summaries
   optuna.TrialPruned
.. module:: optuna.importance

optuna.importance
=================

The :mod:`~optuna.importance` module provides functionality for evaluating hyperparameter importances based on completed trials in a given study. The utility function :func:`~optuna.importance.get_param_importances` takes a ``Study`` and optional evaluator as two of its inputs. The evaluator must derive from :class:`~optuna.importance.BaseImportanceEvaluator`, and is initialized as a :class:`~optuna.importance.FanovaImportanceEvaluator` by default when not passed in. Users implementing custom evaluators should refer to either :class:`~optuna.importance.FanovaImportanceEvaluator` or :class:`~optuna.importance.MeanDecreaseImpurityImportanceEvaluator` as a guide, paying close attention to the format of the return value from the Evaluator's :meth:`evaluate` function.


.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.importance.get_param_importances
   optuna.importance.FanovaImportanceEvaluator
   optuna.importance.MeanDecreaseImpurityImportanceEvaluator
.. module:: optuna.storages

optuna.storages
===============

The :mod:`~optuna.storages` module defines a :class:`~optuna.storages.BaseStorage` class which abstracts a backend database and provides library-internal interfaces to the read/write histories of the studies and trials. Library users who wish to use storage solutions other than the default in-memory storage should use one of the child classes of :class:`~optuna.storages.BaseStorage` documented below.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.storages.RDBStorage
   optuna.storages.RedisStorage
   optuna.storages.RetryFailedTrialCallback
   optuna.storages.fail_stale_trials
.. module:: optuna.study

optuna.study
============

The :mod:`~optuna.study` module implements the :class:`~optuna.study.Study` object and related functions. A public constructor is available for the :class:`~optuna.study.Study` class, but direct use of this constructor is not recommended. Instead, library users should create and load a :class:`~optuna.study.Study` using :func:`~optuna.study.create_study` and :func:`~optuna.study.load_study` respectively.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.study.Study
   optuna.study.create_study
   optuna.study.load_study
   optuna.study.delete_study
   optuna.study.copy_study
   optuna.study.get_all_study_summaries
   optuna.study.MaxTrialsCallback
   optuna.study.StudyDirection
   optuna.study.StudySummary
.. module:: optuna.integration

optuna.integration
==================

The :mod:`~optuna.integration` module contains classes used to integrate Optuna with external machine learning frameworks.

For most of the ML frameworks supported by Optuna, the corresponding Optuna integration class serves only to implement a callback object and functions, compliant with the framework's specific callback API, to be called with each intermediate step in the model training. The functionality implemented in these callbacks across the different ML frameworks includes:

(1) Reporting intermediate model scores back to the Optuna trial using :func:`optuna.trial.report`,
(2) According to the results of :func:`optuna.trial.Trial.should_prune`, pruning the current model by raising :func:`optuna.TrialPruned`, and
(3) Reporting intermediate Optuna data such as the current trial number back to the framework, as done in :class:`~optuna.integration.MLflowCallback`.

For scikit-learn, an integrated :class:`~optuna.integration.OptunaSearchCV` estimator is available that combines scikit-learn BaseEstimator functionality with access to a class-level ``Study`` object.

AllenNLP
--------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.AllenNLPExecutor
   optuna.integration.allennlp.dump_best_config
   optuna.integration.AllenNLPPruningCallback

BoTorch
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.BoTorchSampler
   optuna.integration.botorch.qei_candidates_func
   optuna.integration.botorch.qehvi_candidates_func
   optuna.integration.botorch.qparego_candidates_func

Catalyst
--------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.CatalystPruningCallback

Chainer
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.ChainerPruningExtension
   optuna.integration.ChainerMNStudy

fast.ai
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.FastAIV1PruningCallback
   optuna.integration.FastAIV2PruningCallback
   optuna.integration.FastAIPruningCallback

Keras
-----

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.KerasPruningCallback

LightGBM
--------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.LightGBMPruningCallback
   optuna.integration.lightgbm.train
   optuna.integration.lightgbm.LightGBMTuner
   optuna.integration.lightgbm.LightGBMTunerCV

MLflow
------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.MLflowCallback

Weights & Biases
----------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.WeightsAndBiasesCallback

MXNet
-----

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.MXNetPruningCallback

pycma
-----

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.PyCmaSampler
   optuna.integration.CmaEsSampler

PyTorch
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.PyTorchIgnitePruningHandler
   optuna.integration.PyTorchLightningPruningCallback
   optuna.integration.TorchDistributedTrial

scikit-learn
------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.OptunaSearchCV

scikit-optimize
---------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.SkoptSampler

skorch
------

.. autosummary::
   :toctree: generated/
   :nosignatures:

    optuna.integration.SkorchPruningCallback

TensorFlow
----------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.TensorBoardCallback
   optuna.integration.TensorFlowPruningHook
   optuna.integration.TFKerasPruningCallback

XGBoost
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.integration.XGBoostPruningCallback
.. module:: optuna.distributions

optuna.distributions
====================

The :mod:`~optuna.distributions` module defines various classes representing probability distributions, mainly used to suggest initial hyperparameter values for an optimization trial. Distribution classes inherit from a library-internal :class:`~optuna.distributions.BaseDistribution`, and is initialized with specific parameters, such as the ``low`` and ``high`` endpoints for a :class:`~optuna.distributions.UniformDistribution`.

Optuna users should not use distribution classes directly, but instead use utility functions provided by :class:`~optuna.trial.Trial` such as :meth:`~optuna.trial.Trial.suggest_int`.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.distributions.FloatDistribution
   optuna.distributions.IntDistribution
   optuna.distributions.UniformDistribution
   optuna.distributions.LogUniformDistribution
   optuna.distributions.DiscreteUniformDistribution
   optuna.distributions.IntUniformDistribution
   optuna.distributions.IntLogUniformDistribution
   optuna.distributions.CategoricalDistribution
   optuna.distributions.distribution_to_json
   optuna.distributions.json_to_distribution
   optuna.distributions.check_distribution_compatibility
.. module:: optuna.exceptions

optuna.exceptions
=================

The :mod:`~optuna.exceptions` module defines Optuna-specific exceptions deriving from a base :class:`~optuna.exceptions.OptunaError` class. Of special importance for library users is the :class:`~optuna.exceptions.TrialPruned` exception to be raised if :func:`optuna.trial.Trial.should_prune` returns ``True`` for a trial that should be pruned.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.exceptions.OptunaError
   optuna.exceptions.TrialPruned
   optuna.exceptions.CLIUsageError
   optuna.exceptions.StorageInternalError
   optuna.exceptions.DuplicatedStudyError
.. module:: optuna.cli

optuna.cli
==========

The :mod:`~optuna.cli` module implements Optuna's command-line functionality using the `cliff framework <https://docs.openstack.org/cliff/latest/index.html>`_.

.. autoprogram-cliff:: optuna.cli._OptunaApp
   :application: optuna

.. autoprogram-cliff:: optuna.command
   :application: optuna
.. module:: optuna.pruners

optuna.pruners
==============

The :mod:`~optuna.pruners` module defines a :class:`~optuna.pruners.BasePruner` class characterized by an abstract :meth:`~optuna.pruners.BasePruner.prune` method, which, for a given trial and its associated study, returns a boolean value representing whether the trial should be pruned. This determination is made based on stored intermediate values of the objective function, as previously reported for the trial using :meth:`optuna.trial.Trial.report`. The remaining classes in this module represent child classes, inheriting from :class:`~optuna.pruners.BasePruner`, which implement different pruning strategies.

.. seealso::
    :ref:`pruning` tutorial explains the concept of the pruner classes and a minimal example.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.pruners.BasePruner
   optuna.pruners.MedianPruner
   optuna.pruners.NopPruner
   optuna.pruners.PatientPruner
   optuna.pruners.PercentilePruner
   optuna.pruners.SuccessiveHalvingPruner
   optuna.pruners.HyperbandPruner
   optuna.pruners.ThresholdPruner
.. module:: optuna.samplers

optuna.samplers
===============

The :mod:`~optuna.samplers` module defines a base class for parameter sampling as described extensively in :class:`~optuna.samplers.BaseSampler`. The remaining classes in this module represent child classes, deriving from :class:`~optuna.samplers.BaseSampler`, which implement different sampling strategies.

.. seealso::
    :ref:`pruning` tutorial explains the overview of the sampler classes.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.samplers.BaseSampler
   optuna.samplers.GridSampler
   optuna.samplers.RandomSampler
   optuna.samplers.TPESampler
   optuna.samplers.CmaEsSampler
   optuna.samplers.PartialFixedSampler
   optuna.samplers.NSGAIISampler
   optuna.samplers.MOTPESampler
   optuna.samplers.QMCSampler
   optuna.samplers.IntersectionSearchSpace
   optuna.samplers.intersection_search_space
.. module:: optuna.trial

optuna.trial
============

The :mod:`~optuna.trial` module contains :class:`~optuna.trial.Trial` related classes and functions.

A :class:`~optuna.trial.Trial` instance represents a process of evaluating an objective function. This instance is passed to an objective function and provides interfaces to get parameter suggestion, manage the trial's state, and set/get user-defined attributes of the trial, so that Optuna users can define a custom objective function through the interfaces. Basically, Optuna users only use it in their custom objective functions.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.trial.Trial
   optuna.trial.FixedTrial
   optuna.trial.FrozenTrial
   optuna.trial.TrialState
   optuna.trial.create_trial
.. module:: optuna.logging

optuna.logging
==============

The :mod:`~optuna.logging` module implements logging using the Python ``logging`` package. Library users may be especially interested in setting verbosity levels using :func:`~optuna.logging.set_verbosity` to one of ``optuna.logging.CRITICAL`` (aka ``optuna.logging.FATAL``), ``optuna.logging.ERROR``, ``optuna.logging.WARNING`` (aka ``optuna.logging.WARN``), ``optuna.logging.INFO``, or ``optuna.logging.DEBUG``.


.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.logging.get_verbosity
   optuna.logging.set_verbosity
   optuna.logging.disable_default_handler
   optuna.logging.enable_default_handler
   optuna.logging.disable_propagation
   optuna.logging.enable_propagation
API Reference
=============

.. toctree::
    :maxdepth: 1

    optuna
    cli
    distributions
    exceptions
    importance
    integration
    logging
    multi_objective/index
    pruners
    samplers
    storages
    study
    trial
    visualization/index
.. module:: optuna.visualization

optuna.visualization
====================

The :mod:`~optuna.visualization` module provides utility functions for plotting the optimization process using plotly and matplotlib. Plotting functions generally take a :class:`~optuna.study.Study` object and optional parameters are passed as a list to the ``params`` argument.

.. note::
    In the :mod:`optuna.visualization` module, the following functions use plotly to create figures, but `JupyterLab`_ cannot
    render them by default. Please follow this `installation guide`_ to show figures in
    `JupyterLab`_.

    .. _JupyterLab: https://github.com/jupyterlab/jupyterlab
    .. _installation guide: https://github.com/plotly/plotly.py#jupyterlab-support-python-35

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.visualization.plot_contour
   optuna.visualization.plot_edf
   optuna.visualization.plot_intermediate_values
   optuna.visualization.plot_optimization_history
   optuna.visualization.plot_parallel_coordinate
   optuna.visualization.plot_param_importances
   optuna.visualization.plot_pareto_front
   optuna.visualization.plot_slice
   optuna.visualization.is_available

.. note::
    The following :mod:`optuna.visualization.matplotlib` module uses Matplotlib as a backend.

.. toctree::
    :maxdepth: 1

    matplotlib
.. module:: optuna.visualization.matplotlib

optuna.visualization.matplotlib
===============================

.. note::
    The following functions use Matplotlib as a backend.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    optuna.visualization.matplotlib.plot_contour
    optuna.visualization.matplotlib.plot_edf
    optuna.visualization.matplotlib.plot_intermediate_values
    optuna.visualization.matplotlib.plot_optimization_history
    optuna.visualization.matplotlib.plot_parallel_coordinate
    optuna.visualization.matplotlib.plot_param_importances
    optuna.visualization.matplotlib.plot_pareto_front
    optuna.visualization.matplotlib.plot_slice
    optuna.visualization.matplotlib.is_available
.. module:: optuna.multi_objective.study

optuna.multi_objective.study
============================

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.multi_objective.study.MultiObjectiveStudy
   optuna.multi_objective.study.create_study
   optuna.multi_objective.study.load_study
.. module:: optuna.multi_objective.visualization

optuna.multi_objective.visualization
====================================

.. note::
    :mod:`optuna.multi_objective.visualization` module uses plotly to create figures,
    but `JupyterLab`_ cannot render them by default. Please follow this `installation guide`_ to
    show figures in `JupyterLab`_.

    .. _JupyterLab: https://github.com/jupyterlab/jupyterlab
    .. _installation guide: https://github.com/plotly/plotly.py#jupyterlab-support-python-35

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.multi_objective.visualization.plot_pareto_front
.. module:: optuna.multi_objective.samplers

optuna.multi_objective.samplers
===============================

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.multi_objective.samplers.BaseMultiObjectiveSampler
   optuna.multi_objective.samplers.NSGAIIMultiObjectiveSampler
   optuna.multi_objective.samplers.RandomMultiObjectiveSampler
   optuna.multi_objective.samplers.MOTPEMultiObjectiveSampler
.. module:: optuna.multi_objective.trial

optuna.multi_objective.trial
============================

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optuna.multi_objective.trial.MultiObjectiveTrial
   optuna.multi_objective.trial.FrozenMultiObjectiveTrial
.. module:: optuna.multi_objective

optuna.multi_objective
======================

This module is deprecated, with former functionality moved to :mod:`optuna.samplers`, :mod:`optuna.study`, :mod:`optuna.trial` and :mod:`optuna.visualization`.

.. toctree::
    :maxdepth: 1

    samplers
    study
    trial
    visualization
