# Yellowbrick


[![Build Status](https://github.com/DistrictDataLabs/yellowbrick/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/DistrictDataLabs/yellowbrick/actions/workflows/ci.yml)
[![Coverage Status](https://codecov.io/gh/DistrictDataLabs/yellowbrick/branch/develop/graph/badge.svg?token=BnaSECZz2r)](https://codecov.io/gh/DistrictDataLabs/yellowbrick)
[![Total Alerts](https://img.shields.io/lgtm/alerts/g/DistrictDataLabs/yellowbrick.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/DistrictDataLabs/yellowbrick/alerts/)
[![Language Grade: Python](https://img.shields.io/lgtm/grade/python/g/DistrictDataLabs/yellowbrick.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/DistrictDataLabs/yellowbrick/context:python)
[![PyPI version](https://badge.fury.io/py/yellowbrick.svg)](https://badge.fury.io/py/yellowbrick)
[![Documentation Status](https://readthedocs.org/projects/yellowbrick/badge/?version=latest)](http://yellowbrick.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1206239.svg)](https://doi.org/10.5281/zenodo.1206239)
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.01075/status.svg)](https://doi.org/10.21105/joss.01075)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/DistrictDataLabs/yellowbrick/develop?filepath=examples%2Fexamples.ipynb)


**Visual analysis and diagnostic tools to facilitate machine learning model selection.**

[![Banner](docs/images/readme/banner.png)](https://www.scikit-yb.org/en/latest/gallery.html)

## What is Yellowbrick?

Yellowbrick is a suite of visual diagnostic tools called "Visualizers" that extend the scikit-learn API to allow human steering of the model selection process. In a nutshell, Yellowbrick combines scikit-learn with matplotlib in the best tradition of the scikit-learn documentation, but to produce visualizations for _your_ machine learning workflow!

For complete documentation on the Yellowbrick API, a gallery of available visualizers, the contributor's guide, tutorials and teaching resources, frequently asked questions, and more, please visit our documentation at [www.scikit-yb.org](https://www.scikit-yb.org/).

## Installing Yellowbrick

Yellowbrick is compatible with Python 3.4 or later and also depends on scikit-learn and matplotlib. The simplest way to install Yellowbrick and its dependencies is from PyPI with pip, Python's preferred package installer.

    $ pip install yellowbrick

Note that Yellowbrick is an active project and routinely publishes new releases with more visualizers and updates. In order to upgrade Yellowbrick to the latest version, use pip as follows.

    $ pip install -U yellowbrick

You can also use the `-U` flag to update scikit-learn, matplotlib, or any other third party utilities that work well with Yellowbrick to their latest versions.

If you're using Anaconda (recommended for Windows users), you can take advantage of the conda utility to install Yellowbrick:

    conda install -c districtdatalabs yellowbrick

## Using Yellowbrick

The Yellowbrick API is specifically designed to play nicely with scikit-learn. Here is an example of a typical workflow sequence with scikit-learn and Yellowbrick:

### Feature Visualization

In this example, we see how Rank2D performs pairwise comparisons of each feature in the data set with a specific metric or algorithm and then returns them ranked as a lower left triangle diagram.

```python
from yellowbrick.features import Rank2D

visualizer = Rank2D(
    features=features, algorithm='covariance'
)
visualizer.fit(X, y)                # Fit the data to the visualizer
visualizer.transform(X)             # Transform the data
visualizer.show()                   # Finalize and render the figure
```

### Model Visualization

In this example, we instantiate a scikit-learn classifier and then use Yellowbrick's ROCAUC class to visualize the tradeoff between the classifier's sensitivity and specificity.

```python
from sklearn.svm import LinearSVC
from yellowbrick.classifier import ROCAUC

model = LinearSVC()
visualizer = ROCAUC(model)
visualizer.fit(X,y)
visualizer.score(X,y)
visualizer.show()
```

For additional information on getting started with Yellowbrick, view the [Quick Start Guide](https://www.scikit-yb.org/en/latest/quickstart.html) in the [documentation](https://www.scikit-yb.org/en/latest/) and check out our [examples notebook](https://github.com/DistrictDataLabs/yellowbrick/blob/develop/examples/examples.ipynb).

## Contributing to Yellowbrick

Yellowbrick is an open source project that is supported by a community who will gratefully and humbly accept any contributions you might make to the project. Large or small, any contribution makes a big difference; and if you've never contributed to an open source project before, we hope you will start with Yellowbrick!

If you are interested in contributing, check out our [contributor's guide](https://www.scikit-yb.org/en/latest/contributing/index.html). Beyond creating visualizers, there are many ways to contribute:

- Submit a bug report or feature request on [GitHub Issues](https://github.com/DistrictDataLabs/yellowbrick/issues).
- Contribute a Jupyter notebook to our examples [gallery](https://github.com/DistrictDataLabs/yellowbrick/tree/develop/examples).
- Assist us with [user testing](https://www.scikit-yb.org/en/latest/evaluation.html).
- Add to the documentation or help with our website, [scikit-yb.org](https://www.scikit-yb.org).
- Write [unit or integration tests](https://www.scikit-yb.org/en/latest/contributing/developing_visualizers.html#integration-tests) for our project.
- Answer questions on our issues, mailing list, Stack Overflow, and elsewhere.
- Translate our documentation into another language.
- Write a blog post, tweet, or share our project with others.
- [Teach](https://www.scikit-yb.org/en/latest/teaching.html) someone how to use Yellowbrick.

As you can see, there are lots of ways to get involved and we would be very happy for you to join us! The only thing we ask is that you abide by the principles of openness, respect, and consideration of others as described in the [Python Software Foundation Code of Conduct](https://www.python.org/psf/codeofconduct/).

For more information, checkout the `CONTRIBUTING.md` file in the root of the repository or the detailed documentation at [Contributing to Yellowbrick](https://www.scikit-yb.org/en/latest/contributing/index.html)

## Yellowbrick Datasets

Yellowbrick gives easy access to several datasets that are used for the examples in the documentation and testing. These datasets are hosted in our CDN and must be downloaded for use. Typically, when a user calls one of the data loader functions, e.g. `load_bikeshare()` the data is automatically downloaded if it's not already on the user's computer. However, for development and testing, or if you know you will be working without internet access, it might be easier to simply download all the data at once.

The data downloader script can be run as follows:

    $ python -m yellowbrick.download

This will download the data to the fixtures directory inside of the Yellowbrick site packages. You can specify the location of the download either as an argument to the downloader script (use `--help` for more details) or by setting the `$YELLOWBRICK_DATA` environment variable. This is the preferred mechanism because this will also influence how data is loaded in Yellowbrick.

_Note: Developers who have downloaded data from Yellowbrick versions earlier than v1.0 may experience some problems with the older data format. If this occurs, you can clear out your data cache as follows:_

    $ python -m yellowbrick.download --cleanup

_This will remove old datasets and download the new ones. You can also use the `--no-download` flag to simply clear the cache without re-downloading data. Users who are having difficulty with datasets can also use this or they can uninstall and reinstall Yellowbrick using `pip`._

## Citing Yellowbrick

We would be glad if you used Yellowbrick in your scientific publications! If you do, please cite us using the [citation guidelines](https://www.scikit-yb.org/en/latest/about.html#citing-yellowbrick).

## Affiliations

[![District Data Labs](docs/images/readme/affiliates_ddl.png)](https://districtdatalabs.com/) [![NumFOCUS Affiliated Project](docs/images/readme/affiliates_numfocus.png)](https://numfocus.org)
# Maintainers and Contributors

This file describes how the Yellowbrick project is maintained and provides contact information for key folks in the project.

When creating a pull request, your contribution will be reviewed by one or probably two maintainers who will give you the :+1: when your extension is ready to be merged. Maintainers work hard to ensure that Yellowbrick is a high quality project and that contributors are successful.

For more about how to develop visualizers and contribute features to Yellowbrick, see our [contributor's guide](CONTRIBUTING.md) and the [documentation](https://www.scikit-yb.org/en/latest/contributing/index.html).

For everyone who has [contributed](https://github.com/DistrictDataLabs/yellowbrick/graphs/contributors) in big and in small ways, **thank you!**. Yellowbrick is intended to be a community project, welcoming to new and experienced developers alike. If you would like to become a core contributor you must simply submit a pull request that shows core knowledge of the Yellowbrick library. Usually new Visualizers meet this standard; let the maintainers know you'd like to join the team, and they'll help you work toward it!

## Maintainers

This is a list of the primary project maintainers. Feel free to @ message them in issues and converse with them directly.

- [bbengfort](https://github.com/bbengfort)
- [ndanielsen](https://github.com/ndanielsen)
- [lwgray](https://github.com/lwgray)
- [NealHumphrey](https://github.com/NealHumphrey)
- [jkeung](https://github.com/jkeung)
- [pdamodaran](https://github.com/pdamodaran)

## Core Contributors

This is a list of the core-contributors of the project. Core contributors set the road map and vision of the project. Keep an eye out for them in issues and check out their work to use as inspiration! Most likely they would also be happy to chat and answer questions.

- [rebeccabilbro](https://github.com/rebeccabilbro)
- [mattandahalfew](https://github.com/mattandahalfew)
- [tuulihill](https://github.com/tuulihill)
- [balavenkatesan](https://github.com/balavenkatesan)
- [morganmendis](https://github.com/morganmendis)
- [yzyzy](https://github.com/yzyzy)
- [wagner2010](https://github.com/wagner2010)
- [Juan0001](https://github.com/Juan0001)
- [ccjolley](https://github.com/ccjolley)
- [justjess](https://github.com/justjess)
- [kbelita](https://github.com/kbelita)
- [sanemkabaca](https://github.com/sanemkabaca)
- [Kautumn06](https://github.com/Kautumn06)
- [Zeynepelabiad](https://github.com/Zeynepelabiad)
# Yellowbrick

[![Visualizers](https://github.com/DistrictDataLabs/yellowbrick/raw/develop/docs/images/readme/banner.png)](https://www.scikit-yb.org/)

Yellowbrick is a suite of visual analysis and diagnostic tools designed to facilitate machine learning with scikit-learn. The library implements a new core API object, the `Visualizer` that is an scikit-learn estimator &mdash; an object that learns from data. Similar to transformers or models, visualizers learn from data by creating a visual representation of the model selection workflow.

Visualizer allow users to steer the model selection process, building intuition around feature engineering, algorithm selection and hyperparameter tuning. For instance, they can help diagnose common problems surrounding model complexity and bias, heteroscedasticity, underfit and overtraining, or class balance issues. By applying visualizers to the model selection workflow, Yellowbrick allows you to steer predictive models toward more successful results, faster.

The full documentation can be found at [scikit-yb.org](https://scikit-yb.org/) and includes a [Quick Start Guide](https://www.scikit-yb.org/en/latest/quickstart.html) for new users.

## Visualizers

Visualizers are estimators &mdash; objects that learn from data &mdash; whose primary objective is to create visualizations that allow insight into the model selection process. In scikit-learn terms, they can be similar to transformers when visualizing the data space or wrap a model estimator similar to how the `ModelCV` (e.g. [`RidgeCV`](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.RidgeCV.html), [`LassoCV`](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LassoCV.html)) methods work. The primary goal of Yellowbrick is to create a sensical API similar to scikit-learn. Some of our most popular visualizers include:

### Classification Visualization

- **Classification Report**: a visual classification report that displays a model's precision, recall, and F1 per-class scores as a heatmap
- **Confusion Matrix**: a heatmap view of the confusion matrix of pairs of classes in multi-class classification
- **Discrimination Threshold**: a visualization of the precision, recall, F1-score, and queue rate with respect to the discrimination threshold of a binary classifier
- **Precision-Recall Curve**: plot the precision vs recall scores for different probability thresholds
- **ROCAUC**: graph the receiver operator characteristic (ROC) and area under the curve (AUC)

### Clustering Visualization

- **Intercluster Distance Maps**: visualize the relative distance and size of clusters
- **KElbow Visualizer**: visualize cluster according to the specified scoring function, looking for the "elbow" in the curve.
- **Silhouette Visualizer**: select `k` by visualizing the silhouette coefficient scores of each cluster in a single model

### Feature Visualization

- **Manifold Visualization**: high-dimensional visualization with manifold learning
- **Parallel Coordinates**: horizontal visualization of instances
- **PCA Projection**: projection of instances based on principal components
- **RadViz Visualizer**: separation of instances around a circular plot
- **Rank Features**: single or pairwise ranking of features to detect relationships

### Model Selection Visualization

- **Cross Validation Scores**: display the cross-validated scores as a bar chart with the average score plotted as a horizontal line
- **Feature Importances**: rank features based on their in-model performance
- **Learning Curve**: show if a model might benefit from more data or less complexity
- **Recursive Feature Elimination**: find the best subset of features based on importance
- **Validation Curve**: tune a model with respect to a single hyperparameter

### Regression Visualization

- **Alpha Selection**: show how the choice of alpha influences regularization
- **Cook's Distance**: show the influence of instances on linear regression
- **Prediction Error Plots**: find model breakdowns along the domain of the target
- **Residuals Plot**: show the difference in residuals of training and test data

### Target Visualization

- **Balanced Binning Reference**: generate a histogram with vertical lines showing the recommended value point to the bin data into evenly distributed bins
- **Class Balance**: show the relationship of the support for each class in both the training and test data by displaying how frequently each class occurs as a bar graph the frequency of the classes' representation in the dataset
- **Feature Correlation**: visualize the correlation between the dependent variables and the target

### Text Visualization

- **Dispersion Plot**: visualize how key terms are dispersed throughout a corpus
- **PosTag Visualizer**: plot the counts of different parts-of-speech throughout a tagged corpus
- **Token Frequency Distribution**: visualize the frequency distribution of terms in the corpus
- **t-SNE Corpus Visualization**: uses stochastic neighbor embedding to project documents
- **UMAP Corpus Visualization**: plot similar documents closer together to discover clusters

... and more! Yellowbrick is adding new visualizers all the time so be sure to check out our [examples gallery]https://github.com/DistrictDataLabs/yellowbrick/tree/develop/examples) &mdash; or even the [develop](https://github.com/districtdatalabs/yellowbrick/tree/develop) branch &mdash; and feel free to contribute your ideas for new Visualizers!

## Affiliations
[![District Data Labs](https://github.com/DistrictDataLabs/yellowbrick/raw/develop/docs/images/readme/affiliates_ddl.png)](https://www.districtdatalabs.com/) [![NumFOCUS Affiliated Project](https://github.com/DistrictDataLabs/yellowbrick/raw/develop/docs/images/readme/affiliates_numfocus.png)](https://numfocus.org/)
# Contributing to Yellowbrick

**NOTE: This document is a "getting started" summary for contributing to the Yellowbrick project.** To read the full contributor's guide, please visit the [contributing page](http://www.scikit-yb.org/en/latest/contributing/index.html) in the documentation. Please make sure to read this page carefully to ensure the review process is as smooth as possible and to ensure the greatest likelihood of having your contribution be merged.

For more on the development path, goals, and motivations behind Yellowbrick, check out our developer presentation: [Visualizing Model Selection with Scikit-Yellowbrick: An Introduction to Developing Visualizers](http://www.slideshare.net/BenjaminBengfort/visualizing-model-selection-with-scikityellowbrick-an-introduction-to-developing-visualizers).

## How to Contribute

Yellowbrick is an open source project that is supported by a community who will gratefully and humbly accept any contributions you might make to the project. Large or small, any contribution makes a big difference; and if you've never contributed to an open source project before, we hope you will start with Yellowbrick!

Principally, Yellowbrick development is about the addition and creation of *visualizers* &mdash; objects that learn from data and create a visual representation of the data or model. Visualizers integrate with scikit-learn estimators, transformers, and pipelines for specific purposes and as a result, can be simple to build and deploy. The most common contribution is therefore a new visualizer for a specific model or model family. We'll discuss in detail how to build visualizers later.

Beyond creating visualizers, there are many ways to contribute:

- Submit a bug report or feature request on [GitHub Issues](https://github.com/DistrictDataLabs/yellowbrick/issues).
- Contribute a Jupyter notebook to our examples[ gallery](https://github.com/DistrictDataLabs/yellowbrick/tree/develop/examples).
- Assist us with [user testing](http://www.scikit-yb.org/en/latest/evaluation.html).
- Add to the documentation or help with our website, [scikit-yb.org](http://www.scikit-yb.org).
- Write [unit or integration tests](https://www.scikit-yb.org/en/latest/contributing/developing_visualizers.html#integration-tests) for our project.
- Answer questions on our issues, mailing list, Stack Overflow, and elsewhere.
- Translate our documentation into another language.
- Write a blog post, tweet, or share our project with others.
- [Teach](https://www.scikit-yb.org/en/latest/teaching.html) someone how to use Yellowbrick.

As you can see, there are lots of ways to get involved and we would be very happy for you to join us! The only thing we ask is that you abide by the principles of openness, respect, and consideration of others as described in the [Python Software Foundation Code of Conduct](https://www.python.org/psf/codeofconduct/).

## Getting Started on GitHub

Yellowbrick is hosted on GitHub at https://github.com/DistrictDataLabs/yellowbrick.

The typical workflow for a contributor to the codebase is as follows:

1. **Discover** a bug or a feature by using Yellowbrick.
2. **Discuss** with the core contributes by [adding an issue](https://github.com/DistrictDataLabs/yellowbrick/issues).
3. **Fork** the repository into your own GitHub account.
4. Create a **Pull Request** first thing to [connect with us](https://github.com/DistrictDataLabs/yellowbrick/pulls) about your task.
5. **Code** the feature, write the documentation, add your contribution.
6. **Review** the code with core contributors who will guide you to a high quality submission.
7. **Merge** your contribution into the Yellowbrick codebase.

We believe that *contribution is collaboration* and therefore emphasize *communication* throughout the open source process. We rely heavily on GitHub's social coding tools to allow us to do this. For instance, we use GitHub's [milestone](https://help.github.com/en/articles/about-milestones) feature to focus our development efforts for each Yellowbrick semester, so be sure to check out the issues associated with our [current milestone](https://github.com/districtdatalabs/yellowbrick/milestones)!

Once you have a good sense of how you are going to implement the new feature (or fix the bug!), you can reach out for feedback from the maintainers by creating a [pull request](https://github.com/DistrictDataLabs/yellowbrick/pulls). Please note that if we feel your solution has not been thought out in earnest, or if the PR is not aligned with our [current milestone](https://github.com/districtdatalabs/yellowbrick/milestones) goals, we may reach out to ask that you close the PR so that we can prioritize reviewing the most critical feature requests and bug fixes.

Ideally, any pull request should be capable of resolution within 6 weeks of being opened. This timeline helps to keep our pull request queue small and allows Yellowbrick to maintain a robust release schedule to give our users the best experience possible. However, the most important thing is to keep the dialogue going! And if you're unsure whether you can complete your idea within 6 weeks, you should still go ahead and open a PR and we will be happy to help you scope it down as needed.

If we have comments or questions when we evaluate your pull request and receive no response, we will also close the PR after this period of time. Please know that this does not mean we don't value your contribution, just that things go stale. If in the future you want to pick it back up, feel free to address our original feedback and to reference the original PR in a new pull request.

### Forking the Repository

The first step is to fork the repository into your own account. This will create a copy of the codebase that you can edit and write to. Do so by clicking the **"fork"** button in the upper right corner of the Yellowbrick GitHub page.

Once forked, use the following steps to get your development environment set up on your computer:

1. Clone the repository.

    After clicking the fork button, you should be redirected to the GitHub page of the repository in your user account. You can then clone a copy of the code to your local machine.

    ```
    $ git clone https://github.com/[YOURUSERNAME]/yellowbrick
    $ cd yellowbrick
    ```

    Optionally, you can also [add the upstream remote](https://help.github.com/articles/configuring-a-remote-for-a-fork/) to synchronize with changes made by other contributors:

    ```
    $ git remote add upstream https://github.com/DistrictDataLabs/yellowbrick
    ```

    See "Branching Conventions" below for more on this topic.

2. Create a virtual environment.

    Yellowbrick developers typically use [virtualenv](https://virtualenv.pypa.io/en/stable/) (and [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/), [pyenv](https://github.com/pyenv/pyenv-virtualenv) or [conda envs](https://conda.io/docs/using/envs.html) in order to manage their Python version and dependencies. Using the virtual environment tool of your choice, create one for Yellowbrick. Here's how with virtualenv:

    ```
    $ virtualenv venv
    ```

3. Install dependencies.

    Yellowbrick's dependencies are in the `requirements.txt` document at the root of the repository. Open this file and uncomment the dependencies that are for development only. Then install the dependencies with `pip`:

    ```
    $ pip install -r requirements.txt
    ```

    Note that there may be other dependencies required for development and testing, you can simply install them with `pip`. For example to install
    the additional dependencies for building the documentation or to run the
    test suite, use the `requirements.txt` files in those directories:

    ```
    $ pip install -r tests/requirements.txt
    $ pip install -r docs/requirements.txt
    ```

4. Switch to the develop branch.

    The Yellowbrick repository has a `develop` branch that is the primary working branch for contributions. It is probably already the branch you're on, but you can make sure and switch to it as follows::

    ```
    $ git fetch
    $ git checkout develop
    ```

At this point you're ready to get started writing code!

### Branching Conventions

The Yellowbrick repository is set up in a typical production/release/development cycle as described in "[A Successful Git Branching Model](http://nvie.com/posts/a-successful-git-branching-model/)." The primary working branch is the `develop` branch. This should be the branch that you are working on and from, since this has all the latest code. The `master` branch contains the latest stable version and release, _which is pushed to PyPI_. No one but maintainers will push to master.

**NOTE:** All pull requests should be into the `yellowbrick/develop` branch from your forked repository.

You should work directly in your fork and create a pull request from your fork's develop branch into ours. We also recommend setting up an `upstream` remote so that you can easily pull the latest development changes from the main Yellowbrick repository (see [configuring a remote for a fork](https://help.github.com/articles/configuring-a-remote-for-a-fork/)). You can do that as follows:

```
$ git remote add upstream https://github.com/DistrictDataLabs/yellowbrick.git
$ git remote -v
origin    https://github.com/YOUR_USERNAME/YOUR_FORK.git (fetch)
origin    https://github.com/YOUR_USERNAME/YOUR_FORK.git (push)
upstream  https://github.com/DistrictDataLabs/yellowbrick.git (fetch)
upstream  https://github.com/DistrictDataLabs/yellowbrick.git (push)
```

When you're ready, request a code review for your pull request. Then, when reviewed and approved, you can merge your fork into our main branch. Make sure to use the "Squash and Merge" option in order to create a Git history that is understandable.

**NOTE to maintainers**: When merging a pull request, use the "squash and merge" option and make sure to edit the both the subject and the body of the commit message so that when we're putting the changelog together, we know what happened in the PR. I recommend reading [Chris Beams' _How to Write a Git Commit Message_](https://chris.beams.io/posts/git-commit/) so we're all on the same page!

Core contributors and those who are planning on contributing multiple PRs might want to consider using feature branches to reduce the number of merges (and merge conflicts). Create a feature branch as follows:

```
$ git checkout -b feature-myfeature develop
$ git push --set-upstream origin feature-myfeature
```

Once you are done working (and everything is tested) you can submit a PR from your feature branch. Synchronize with `upstream` once the PR has been merged and delete the feature branch:

```
$ git checkout develop
$ git pull upstream develop
$ git push origin develop
$ git branch -d feature-myfeature
$ git push origin --delete feature-myfeature
```

Head back to Github and checkout another issue!

## Developing Visualizers

In this section, we'll discuss the basics of developing visualizers. This of course is a big topic, but hopefully these simple tips and tricks will help make sense.

One thing that is necessary is a good understanding of scikit-learn and Matplotlib. Because our API is intended to integrate with scikit-learn, a good start is to review ["APIs of scikit-learn objects"](http://scikit-learn.org/stable/developers/contributing.html#apis-of-scikit-learn-objects) and ["rolling your own estimator"](http://scikit-learn.org/stable/developers/contributing.html#rolling-your-own-estimator). In terms of matplotlib, check out [Nicolas P. Rougier's Matplotlib tutorial](https://www.labri.fr/perso/nrougier/teaching/matplotlib/).

### Visualizer API

There are two basic types of Visualizers:

- **Feature Visualizers** are high dimensional data visualizations that are essentially transformers.
- **Score Visualizers** wrap a scikit-learn regressor, classifier, or clusterer and visualize the behavior or performance of the model on test data.

These two basic types of visualizers map well to the two basic estimator objects in scikit-learn:

- **Transformers** take input data and return a new data set.
- **Models** are fit to training data and can make predictions.

The scikit-learn API is object oriented, and estimators are initialized with parameters by instantiating their class. Hyperparameters can also be set using the `set_attrs()` method and retrieved with the corresponding `get_attrs()` method. All scikit-learn estimators have a `fit(X, y=None)` method that accepts a two dimensional data array, `X`, and optionally a vector `y` of target values. The `fit()` method trains the estimator, making it ready to transform data or make predictions. Transformers have an associated `transform(X)` method that returns a new dataset, `Xprime` and models have a `predict(X)` method that returns a vector of predictions, `yhat`. Models may also have a `score(X, y)` method that evaluate the performance of the model.

Visualizers interact with scikit-learn objects by intersecting with them at the methods defined above. Specifically, visualizers perform actions related to `fit()`, `transform()`, `predict()`, and `score()` then call a `draw()` method which initializes the underlying figure associated with the visualizer. The user calls the visualizer's `show()` method, which in turn calls a `finalize()` method on the visualizer to draw legends, titles, etc. and then `show()` renders the figure. The Visualizer API is therefore:

- `draw()`: add visual elements to the underlying axes object
- `finalize()`: prepare the figure for rendering, adding final touches such as legends, titles, axis labels, etc.
- `show()`: render the figure for the user.

Creating a visualizer means defining a class that extends `Visualizer` or one of its subclasses, then implementing several of the methods described above. A barebones implementation is as follows::

```python
import matplotlib.pyplot as plot

from yellowbrick.base import Visualizer

class MyVisualizer(Visualizer):

    def __init__(self, ax=None, **kwargs):
        super(MyVisualizer, self).__init__(ax, **kwargs)

    def fit(self, X, y=None):
        super(MyVisualizer, self).fit(X, y)
        self.draw(X)
        return self

    def draw(self, X):
        self.ax.plot(X)
        return self.ax

    def finalize(self):
        self.set_title("My Visualizer")
```

This simple visualizer simply draws a line graph for some input dataset X, intersecting with the scikit-learn API at the `fit()` method. A user would use this visualizer in the typical style::

```python
visualizer = MyVisualizer()
visualizer.fit(X)
visualizer.show()
```

Score visualizers work on the same principle but accept an additional required `model` argument. Score visualizers wrap the model (which can be either instantiated or uninstantiated) and then pass through all attributes and methods through to the underlying model, drawing where necessary.

### Testing

The test package mirrors the `yellowbrick` package in structure and also contains several helper methods and base functionality. To add a test to your visualizer, find the corresponding file to add the test case, or create a new test file in the same place you added your code.

Visual tests are notoriously difficult to create --- how do you test a visualization or figure? Moreover, testing scikit-learn models with real data can consume a lot of memory. Therefore the primary test you should create is simply to test your visualizer from end to end and make sure that no exceptions occur. To assist with this, we have a helper, `VisualTestCase`. Create your unit test as follows::

```python
import pytest

from yellowbrick.datasets import load_occupancy

from tests.base import VisualTestCase

class MyVisualizerTests(VisualTestCase):

    def test_my_visualizer(self):
        """
        Test MyVisualizer on a real dataset
        """
        # Load the data
        X,y = load_occupancy()

        try:
            visualizer = MyVisualizer()
            visualizer.fit(X)
            visualizer.show()
        except Exception as e:
            pytest.fail("my visualizer didn't work")
```

The entire test suite can be run as follows::

```
$ pytest
```

You can also run your own test file as follows::

```
$ pytest tests/test_your_visualizer.py
```

The Makefile uses the pytest runner and testing suite as well as the coverage library, so make sure you have those dependencies installed!

**Note**: Advanced developers can use our _image comparison tests_ to assert that an image generated matches a baseline image. Read more about this in our [testing documentation](https://www.scikit-yb.org/en/latest/contributing/developing_visualizers.html#image-comparison-tests).

### Documentation

The initial documentation for your visualizer will be a well structured docstring. Yellowbrick uses Sphinx to build documentation, therefore docstrings should be written in reStructuredText in numpydoc format (similar to scikit-learn). The primary location of your docstring should be right under the class definition, here is an example::

```python
class MyVisualizer(Visualizer):
    """
    This initial section should describe the visualizer and what
    it's about, including how to use it. Take as many paragraphs
    as needed to get as much detail as possible.

    In the next section describe the parameters to __init__.

    Parameters
    ----------

    model : a scikit-learn regressor
        Should be an instance of a regressor, and specifically one whose name
        ends with "CV" otherwise a will raise a YellowbrickTypeError exception
        on instantiation. To use non-CV regressors see:
        ``ManualAlphaSelection``.

    ax : matplotlib Axes, default: None
        The axes to plot the figure on. If None is passed in the current axes
        will be used (or generated if required).

    kwargs : dict
        Keyword arguments that are passed to the base class and may influence
        the visualization as defined in other Visualizers.

    Examples
    --------

    >>> model = MyVisualizer()
    >>> model.fit(X)
    >>> model.show()

    Notes
    -----

    In the notes section specify any gotchas or other info.
    """
```

This is a very good start to producing a high quality visualizer, but unless it is part of the documentation on our website, it will not be visible. For details on including documentation in the `docs` directory see the [Contributing Documentation](https://www.scikit-yb.org/en/latest/contributing/index.html) section in the larger contributing guide.
---
title: 'Yellowbrick: Visualizing the Scikit-Learn Model Selection Process'
tags:
  - machine learning
  - visual analysis
  - model selection
  - python
  - scikit-learn
  - matplotlib
authors:
  - name: Benjamin Bengfort
    orcid: 0000-0003-0660-7682
    affiliation: 1
  - name: Rebecca Bilbro
    orcid: 0000-0002-1143-044X
    affiliation: 1
affiliations:
 - name: Georgetown University
   index: 1
date: 30 July 2018
bibliography: paper.bib
---

# Summary

Discussions of machine learning are frequently characterized by a singular focus on algorithmic behavior. Be it logistic regression, random forests, Bayesian methods, or artificial neural networks, practitioners are often quick to express their preference. However, model selection is more nuanced than simply picking the “right” or “wrong” algorithm. In practice, the workflow includes multiple iterations through feature engineering, algorithm selection, and hyperparameter tuning &mdash; summarized by Kumar et al. as a search for the maximally performing model selection triple [@kumar2016model]. “Model selection,” they explain, “is iterative and exploratory because the space of [model selection triples] is usually infinite, and it is generally impossible for analysts to know a priori which [combination] will yield satisfactory accuracy and/or insights.”

Treating model selection as search has led to automation through grid search methods, standardized APIs, drag and drop GUIs, and specialized database systems. However, the search problem is computationally intractable and research in both machine learning [@wickham_visualizing_2015] and visual analytics [@liu_wang_liu_zhu_2017] suggests human intuition and guidance can more effectively hone in on quality models than exhaustive optimization methods. By visualizing the model selection process, data scientists can interactively steer towards final, interpretable models and avoid pitfalls and traps [@kapoor2010interactive].

Yellowbrick is a response to the call for open source visual steering tools. For data scientists, Yellowbrick helps evaluate the stability and predictive value of machine learning models and improves the speed of the experimental workflow. For data engineers, Yellowbrick provides visual tools for monitoring model performance in real world applications. For users of models, Yellowbrick provides visual interpretation of the behavior of the model in high dimensional feature space. Finally, for students, Yellowbrick is a framework for understanding a large variety of algorithms and methods.

Implemented in Python, the Yellowbrick visualization package achieves steering by extending both scikit-learn [@sklearn] and Matplotlib [@matplotlib]. Like Yellowbrick, both scikit-learn and Matplotlib are extensions of SciPy [@scipy], libraries intended to facilitate scientific computing. Scikit-learn provides a generalized API for machine learning by exposing the concept of an `Estimator`, an object that learns from data. Yellowbrick in turn extends this concept with the idea of a `Visualizer`, an object that both learns from data and visualizes the result. Visualizers wrap Matplotlib procedures to produce publication-ready figures and rich visual analytics.

Because Yellowbrick is part of a rich visual and machine learning ecosystem, it provides visualizations for feature and target analysis, classification, regression, and clustering model visualization, hyperparameter tuning, and text analysis. A few selected examples of visual diagnostics for model selection and their interpretations follow.

![Feature Analysis](figures/feature_analysis.png)

Because “more data beats better algorithms” [@rajaraman2008more], the first step to creating valid, predictive models is to find the minimum set of features that predicts the dependent variable. Generally, this means finding features that describe data in high dimensional space that are *separable* (i.e., by a hyperplane). Tools like `RadViz`, `ParallelCoordinates`, and `Manifold` help visualize high dimensional data for quick diagnostics. Bayesian models and regressions suffer when independent variables are collinear (i.e., exhibit pairwise correlation). `Rank2D` visualizations show pairwise correlations among features and can facilitate feature elimination.

![Regression Model Tuning](figures/regression.png)

Regression models hypothesize some underlying function influenced by noise whose central tendency can be inferred. The `PredictionError` visualizer shows the relationship of actual to predicted values, giving a sense of heteroskedasticity in the target, or regions of more or less error as predictions deviate from the 45 degree line. The `ResidualsPlot` shows the relationship of error in the training and test data and can also show regions of increased variability in the predictive model.

![Classification Model Tuning](figures/classification.png)

Classification analysis focuses on the precision and recall of the model's prediction of individual classes. The `ClassificationReport` visualizer allows for rapid comparison between models as a visual heatmap of these metrics. The `DiscriminationThreshold` visualizer for binary classifiers shows how adjusting the threshold for positive classification may influence precision and recall globally, as well as the number of points that may require manual checking for stricter determination.

![Clustering Model Tuning](figures/clustering.png)

Searching for structure in unlabelled data can be challenging because evaluation is largely qualitative. When using K-Means models, choosing K has a large impact on the quality of the analysis; the `KElbowVisualizer` can help select the best K given computational constraints. The `SilhouetteVisualizer` shows the relationship of points in each cluster relative to other clusters and gives an overview of the composition and size of each cluster which may hint at how models group similar data points.

![Hyperparameter Tuning](figures/hyperparameter_tuning.png)

Yellowbrick also offers several other techniques for hyperparameter tuning. Model and regression-specific `AlphaSelection` visualizers help identify the impact of regularization on linear models and the influence of complexity on the trade-off between error due to bias or variance. More generally, the `LearningCurve` visualizer shows how sensitive models are to the amount of data the model is trained on.

Yellowbrick includes many more visualizations, intended to fit directly into the machine learning workflow, and many more are being added in each new release. From text analysis-specific visualizations to missing data analysis, to a `contrib` module that focuses on other machine learning libraries, Yellowbrick has tools to facilitate all parts of hypothesis driven development. The source code for Yellowbrick has been archived to Zenodo and the most recent version can be obtained with the linked DOI: [@zenodo].

# Acknowledgements

Since we first introduced the idea of Yellowbrick at PyCon 2016, many people have joined us and stuck with us through 12 releases, ensuring the success of the project. Nathan Danielsen joined very early on and was one of our first maintainers, bringing an engineering perspective to our work and giving us much needed stability in testing. Larry Gray, Neal Humphrey, Jason Keung, Prema Roman, Kristen McIntyre, Jessica D'Amico and Adam Morris have also all joined our project as maintainers and core contributors, and we can't thank them enough.

Yellowbrick would not be possible without the invaluable contributions of those in the Python and Data Science communities. At the time of this writing, GitHub reports that 46 contributors have submitted pull requests that have been merged and released, and we expect this number to continue to grow. Every week, users submit feature requests, bug reports, suggestions and questions that allow us to make the software better and more robust. Others write blog posts about using Yellowbrick, encouraging both newcomers and seasoned practitioners to more fully understand the models they are fitting. Our sincere thanks to the community for their ongoing support and participation.

# References
# Yellowbrick Tests

*Welcome to the Yellowbrick tests!*

If you're looking for information about how to use Yellowbrick, for our contributor's guide, for examples and teaching resources, for answers to frequently asked questions, and more, please visit the latest version of our documentation at [www.scikit-yb.org](https://www.scikit-yb.org/).

## Running Yellowbrick Tests

To run the tests locally, first install the tests-specific requirements with `pip` using the `requirements.txt` file in the `tests` directory:

```
$ pip install -r tests/requirements.txt
```

The required dependencies for the test suite include testing utilities and libraries such as `pandas` and `nltk` that are not included in the core dependencies.

Tests can then be run as follows from the project `root`:

```bash
$ make test
```

The Makefile uses the `pytest` runner and testing suite as well as the coverage library.

## Adding a Test for Your Visualizer

The `tests` package mirrors the yellowbrick package in structure and also contains several helper methods and base functionality. To add a test to your visualizer, find the corresponding file to add the test case, or create a new test file in the same place you added your code.

### Visual Tests

The primary test you should create is simply to test your visualizer from end to end and make sure that no exceptions occur.

Visual tests are notoriously difficult to create --- how do you test a visualization or figure? Moreover, testing scikit-learn models with real data can consume a lot of memory. To assist with this, we have two primary helpers, `VisualTestCase` and the `yellowbrick.datasets` module.

Leverage these helpers to create your tests as follows:

```python
import pytest

from tests.base import VisualTestCase
from yellowbrick.datasets import load_occupancy


class MyVisualizerTests(VisualTestCase):

    def test_my_visualizer(self):
        """
        Test MyVisualizer on a real dataset
        """
        # Load the data using the Yellowbrick datasets module
        X, y = load_occupancy()

        try:
            visualizer = MyVisualizer()
            visualizer.fit(X)
            visualizer.finalize()
        except Exception as e:
            pytest.fail("my visualizer didn't work")
```

### Image Comparison Tests

Writing an image-based comparison test is only a little more difficult than the simple test case presented above. We have adapted `matplotlib`'s image comparison test utility into an easy to use assert method: `self.assert_images_similar(visualizer)`

The main consideration is that you must specify the “baseline” (i.e. expected) image in the `tests/baseline_images/` folder structure.

For example, let's say you create your tests in `tests/test_regressor/test_myvisualizer.py` as follows:

```python
from tests.base import VisualTestCase
...
    def test_my_visualizer_output(self):
        ...
        visualizer = MyVisualizer()
        visualizer.fit(X)
        visualizer.finalize()
        self.assert_images_similar(visualizer)
```

The first time this test is run, there will be no baseline image to compare against, so the test will fail. Alternatively, if you are making a correction to the existing test `test_my_visualizer_output`, and the correction modifies the resulting test image, the test may also fail to match the existing baseline image. The solution is to first run the tests, then copy the new output images to the correct subdirectory under source code revision control (with `git add`). When rerunning the tests, they should now pass!

We have a helper script, `tests/images.py` to clean up and manage baseline images automatically. It is run using the ``python -m`` command to execute a module as main, and it takes as an argument the path to **your** test file. To copy the figures as above:

```bash
$ python -m tests.images tests/test_regressor/test_myvisualizer.py
```

This will move all related test images from `actual_images` to `baseline_images` on your behalf (note you'll have had to already run the tests at least once to generate the images). You can also clean up images from both actual and baseline as follows:

```bash
$ python -m tests.images -C tests/test_regressor/test_myvisualizer.py
```

This is useful particularly if you're stuck trying to get an image comparison to work. For more information on the images helper script, use `python -m tests.images --help`.
**Deployed:** DayOfWeek, Month D, YYYY
**Current Contributors:** [insert GitHub @username comma list]

[Release description (1-2 paras)]

**Major Changes:**

- [Change 1]
- [Change 2]

Minor Changes:

- [Change 1]
- [Change 2]<!--
# Welcome Contributor!

Thank you for contributing to Yellowbrick, please follow the instructions below to get
your PR started off on the right foot.

## First Steps

1. Are you merging from a feature branch into develop?

    _If not, please create a feature branch and change your PR to merge from that branch
    into the Yellowbrick `develop` branch._

2. Does your PR have a title?

    _Please ensure your PR has a short, informative title, e.g. "Enhances ParallelCoordinates with new andrews_curve parameter" or "Corrects bug in WhiskerPlot that causes index error"_

3. Summarize your PR (HINT: See CHECKLIST/TEMPLATE below!)
-->

This PR fixes #issue_number _(If you are fixing a bug)_ which reported a bug that caused a problem to occur when users...

_(or if you are introducing a new feature)_ which requested a feature to allow the user to...

I have made the following changes:

1.
2.
3.

### Sample Code and Plot

_If you are adding or modifying a visualizer, PLEASE include a sample plot here along with the code you used to generate it._

### TODOs and questions

<!--
If this is a work-in-progress (WIP), list the changes you still need to make and/or questions or the Yellowbrick team. You can also mention extensions to your work that might be added as an issue to work on after the PR.
-->

Still to do:

- [ ]
- [ ]
- [ ]

Questions for the @DistrictDataLabs/team-oz-maintainers:

- [ ]
- [ ]

### CHECKLIST

<!--
Here's a handy checklist to go through before submitting a PR, note that you can check a checkbox in Markdown by changing `- [ ]` to `- [x]` or you can create the PR and check the box manually.
-->

- [ ] _Is the commit message formatted correctly?_
- [ ] _Have you noted the new functionality/bugfix in the release notes of the next release?_

<!-- If you've changed any code -->

- [ ] _Included a sample plot to visually illustrate your changes?_
- [ ] _Do all of your functions and methods have docstrings?_
- [ ] _Have you added/updated unit tests where appropriate?_
- [ ] _Have you updated the baseline images if necessary?_
- [ ] _Have you run the unit tests using `pytest`?_
- [ ] _Is your code style correct (are you using PEP8, pyflakes)?_
- [ ] _Have you documented your new feature/functionality in the docs?_

<!-- If you've added to the docs -->

- [ ] _Have you built the docs using `make html`?_
---
name: Generic Issue
about: Alert us to an issue or task for this project

---

**Describe the issue**
A clear and concise description of what the issue is.

<!-- If you have a question, note that you can email us via our listserve:
     https://groups.google.com/forum/#!forum/yellowbrick -->

<!-- This line alerts the Yellowbrick maintainers, feel free to use this
     @ address to alert us directly in follow up comments -->
@DistrictDataLabs/team-oz-maintainers
---
name: Bug report
about: Create a report to help us improve

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**

```python
# Steps to reproduce the behavior (code snippet):
# Should include imports, dataset loading, and execution
# Add the traceback below
```

**Dataset**
Did you use a specific dataset to produce the bug? Where can we access it?

**Expected behavior**
A clear and concise description of what you expected to happen.

**Traceback**

```
If applicable, add the traceback from the exception.
```

**Desktop (please complete the following information):**
 - OS: [e.g. macOS]
 - Python Version [e.g. 2.7, 3.6, miniconda]
 - Yellowbrick Version [e.g. 0.7]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project

---

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Examples**
Attach an image with an example of the proposed visualization (e.g. from wikipedia or a wireframe example you've drawn). If there are any references we should include in the documentation, please also include them.
# Yellowbrick Examples

[![Visualizers](../docs/images/readme/banner.png)](../docs/images/readme/banner.png)

Welcome to the yellowbrick examples directory! This directory contains a gallery of visualizers and their application to classification, regression, clustering, and other machine learning techniques with scikit-learn. Examples have been submitted both by the Yellowbrick team and also users like you! The result is a rich gallery of tools and techniques to equip your machine learning with visual diagnostics and visualizer workflows!

## Getting Started

The notebook to explore first is the `examples.ipynb` Jupyter notebook. This notebook contains the executable examples from the tutorial in the documentation. You can run the notebook as follows:

```
$ jupyter notebook examples.ipynb
```

If you don't have jupyter installed, or other dependencies, you may have to `pip install` them.

## Organization

The examples directory contains many notebooks, folders and files. At the top level you will see the following:

- examples.ipynb: a notebook with executable versions of the tutorial visualizers
- palettes.ipynb: a visualization of the Yellowbrick palettes
- regression.ipynb: a notebook exploring the regression model visualizers.

In addition to these files and directory, you will see many other directories, whose names are the GitHub usernames of their contributors. You can explore these user submitted examples or submit your own!

### Contributing

To contribute an example notebook of your own, perform the following steps:

1. Fork the repository into your own account
2. Checkout the develop branch (see [contributing to Yellowbrick](http://www.scikit-yb.org/en/latest/about.html#contributing) for more.
3. Create a directory in the repo, `examples/username` where username is your GitHub username.
4. Create a notebook in that directory with your example. See [user testing](http://www.scikit-yb.org/en/latest/evaluation.html) for more.
5. Commit your changes back to your fork.
6. Submit a pull-request from your develop branch to the Yellowbrick develop branch.
7. Complete the code review steps with a Yellowbrick team member.

That's it -- thank you for contributing your example!

A couple of notes. First, please make sure that the Jupyter notebook you submit is "run" -- that is it has the output saved to the notebook and is viewable on GitHub (empty notebooks don't serve well as a gallery). Second, please do not commit datasets, but instead provide instructions for downloading the dataset. You can create a downloader utility similar to ours.

One great tip, is to create your PR right after you fork the repo; that way we can work with you on the changes you're making and communicate about how to have a very successful contribution!

    "#This is a test for District Data Labs Yellowbrick"
Zipped data file available here: http://archive.ics.uci.edu/ml/machine-learning-databases/00315/
Extract just the two .txt files into this folder# Yellowbrick Documentation

*Welcome to the Yellowbrick docs!*

If you're looking for information about how to use Yellowbrick, for our contributor's guide, for examples and teaching resources, for answers to frequently asked questions, and more, please visit the latest version of our documentation at [www.scikit-yb.org](https://www.scikit-yb.org/).

## Building the Docs

To build the documents locally, first install the documentation-specific requirements with `pip` using the `requirements.txt` file in the `docs` directory:

```bash
$ pip install -r docs/requirements.txt
```

You will then be able to build the documentation from inside the `docs` directory by running `make html`; the documentation will be built and rendered in the `_build/html` directory. You can view it by opening `_build/html/index.html` then navigating to your documentation in the browser.

## reStructuredText

Yellowbrick uses [Sphinx](http://www.sphinx-doc.org/en/master/index.html) to build our documentation. The advantages of using Sphinx are many; we can more directly link to the documentation and source code of other projects like Matplotlib and scikit-learn using [intersphinx](http://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html). In addition, docstrings used to describe Yellowbrick visualizers can be automatically included when the documentation is built via [autodoc](http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#sphinx.ext.autodoc).

To take advantage of these features, our documentation must be written in reStructuredText (or "rst"). reStructuredText is similar to markdown, but not identical, and does take some getting used to. For instance, styling for things like codeblocks, external hyperlinks, internal cross references, notes, and fixed-width text are all unique in rst.

If you would like to contribute to our documentation and do not have prior experience with rst, we recommend you make use of these resources:

- [A reStructuredText Primer](http://docutils.sourceforge.net/docs/user/rst/quickstart.html)
- [rst notes and cheatsheet](https://cheat.readthedocs.io/en/latest/rst.html)
- [Using the plot directive](https://matplotlib.org/devel/plot_directive.html)

## Adding New Visualizers to the Docs

If you are adding a new visualizer to the docs, there are quite a few examples in the documentation on which you can base your files of similar types.

The primary format for the API section is as follows:

```
.. -*- mode: rst -*-

My Visualizer
=============

A brief introduction to my visualizer and how it is useful in the machine learning process.

.. plot::
    :context: close-figs
    :include-source: False
    :alt: Example using MyVisualizer

    visualizer = MyVisualizer(LinearRegression())

    visualizer.fit(X, y)
    g = visualizer.show()

Discussion about my visualizer and some interpretation of the above plot.


API Reference
-------------

.. automodule:: yellowbrick.regressor.mymodule
    :members: MyVisualizer
    :undoc-members:
    :show-inheritance:
```

This is a pretty good structure for a documentation page; a brief introduction followed by a code example with a visualization included using [the plot directive](https://matplotlib.org/devel/plot_directive.html). This will render the `MyVisualizer` image in the document along with links for the complete source code, the png, and the pdf versions of the image. It will also have the "alt-text" (for screen-readers) and will not display the source because of the `:include-source:` option. If `:include-source:` is omitted, the source will also be included.

The primary section is wrapped up with a discussion about how to interpret the visualizer and use it in practice. Finally the `API Reference` section will use `automodule` to include the documentation from your docstring.

There are several other places where you can list your visualizer, but to ensure it is included in the documentation it *must be listed in the TOC of the local index*. Find the `index.rst` file in your subdirectory and add your rst file (without the `.rst` extension) to the `..toctree::` directive. This will ensure your documentation is included when it is built.

## Generating the Gallery

In v1.0, we have adopted Matplotlib's [plot directive](https://matplotlib.org/devel/plot_directive.html) which means that the majority of the images generated for the documentation are generated automatically. One exception is the gallery; the images for the gallery must still be generated manually.

If you have contributed a new visualizer as described in the above section, please also add it to the gallery, both to `docs/gallery.py` and to `docs/gallery.rst`. (Make sure you have already installed Yellowbrick in editable mode, from the top level directory: `pip install -e` .)

If you want to regenerate a single image (e.g. the elbow curve plot), you can do so as follows:

```bash
$ python docs/gallery.py elbow
```

If you want to regenerate them all (note: this takes a long time!)

```bash
$ python docs/gallery.py all
```
