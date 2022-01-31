VERSION 1.1.1
=============

Date: 19/07/2020

Mostly doc typos and some minor bug fixes. Future versions (if any)
will not support Python 2 anymore.

Bug Fixes
---------

* 'mse' is now available in GridSearCV and RandomizedSearchCV
* The Jester dataset link was updated
* Fixed a potential race condition when creating dataset directories


VERSION 1.1.0
=============

Date: 13/11/2019

1.1.0 will be the last stable version with new features. Next versions will
only provide bug-fixes, but no new features. (And probably not support
Python 2 at all).

Enhancements
------------

* The prompt confirmation can now be disabled when downloading a dataset.
* The MSE metric has been added.

Bug Fixes
---------

* Fixed a bug where msd and peasron would not properly set the similarity to
  zero when ``min_support`` wasn't reached.

API Changes
-----------

* Tools that were deprecated before (data.split(), GridSearch, evaluate) are
  now removed.

VERSION 1.0.6
=============

Date: 22/04/18

Enhancements
------------

* Added verbose option to algorithms using a similarity matrix or baseline
  computation, to avoid unwanted printed messages.
* When PredictionImpossible is raised, the prediction is now deferred to
  default_prediction() method, which can be overridden is child classes. This
  allows to not always set the default prediction to the average rating, which
  can be useful for some algorithms (e.g. those working with implicit positive
  feedback).
* LeaveOneOut() now accepts a min_n_ratings parameter to make sure users in the
  trainset have at least min_n_ratings ratings.
* Dumping is now done with pickle's highest protocol which allows for larger
  files.

Bug Fixes
---------

* Joblib parameter `n_jobs` now defaults to 1 (no use of multiprocessing).
  Should fix issues with Windows users.
* `cross_validate` now returns correct values for training measures (used to
  return test measures instead).

VERSION 1.0.5
=============

Date: 09/01/18

Enhancements
------------

* Cross-validation tools have been entirely reworked. We can now rely on
  powerful and flexible cross-validation iterators, inspired by scikit-learn's
  API.
* the evaluate() method has been replaced by cross-validate which is parallel
  and can return measures on trainset as well as computation times.
* GridSearch is now parallel, using joblib.
* GridSearch now allows to refit an algorithm on the whole dataset.
* default data directory can now be custom with env variable
  SURPRISE_DATA_FOLDER
* the fit() (and train()) methods now return self, which allows one-liners like
  algo.fit(trainset).test(testset)
* Algorithms using a random initialization (e.g. SVD, NMF, CoClustering) now
  have a random_state parameter for seeding the RNG.
* The getting started guide has been rewritten

API Changes
-----------

* The train() method is now deprecated and replaced by the fit() method (same
  signature). Calls to train() should still work as before.
* Using data.split() or accessing the data.folds() generator is deprecated and
  replaced by the use of the more powefull CV iterators.
* evaluate() is deprecated and  replaced by model_selection.cross_validate(),
  which is parallel.
* GridSearch is deprecated and replaced by model_selection.GridSearchCV()

VERSION 1.0.4
=============

Date: 20/09/17

Enhancements
------------

* Added possibility to load a dataset from a pandas dataframe
* Added Precision and Recall examples to the FAQ (Maher Malaeb)
* Added a kNN algorithm with normalization by z-score (Hengji Liu)
* kNN algorithms now use heapq instead of list.sort() (computation time
  enhancement for large datasets).

Fixes
-----

* Prediciont.__str__() when r_ui is None
* GridSearch for dict parameters is now working as expected

API Changes
-----------

* param_grid for GridSearch is now slightly different for dict parameters (see
  note on [the
  docs](http://surprise.readthedocs.io/en/stable/getting_started.html#tune-algorithm-parameters-with-gridsearch)).

VERSION 1.0.3
=============

Date: 03/05/17

Enhancements
------------

* Added FAQ in the doc
* Added the possibility to retrieve the k nearest neighbors of a user or an
  item.
* Changed the dumping process a bit (see API changes). Plus, dumps can now be
  loaded.
* Added possibility to build a testset from the ratings of a training set
* Added inner-to-raw id conversion in the Trainset class
* The r_ui parameter of the predict() method is now optional

Fixes
-----
* Fixed verbosity of the evaluate function
* Corrected prediction when only user (or only item) is unknown in SVD and NMF
  algorithms. Thanks to kenoung!
* Corrected factor vectors initialization of SVD algorithms. Thanks to
  adideshp!

API Changes
-----------

* The dump() method now dumps a list of predition (optional) and an algorithm
  (optional as well). The algorithm is now a real algorithm object. The
  trainset is not dumped anymore as it is already part of the algorithm anyway.
* The dump() method is now part of the dump namespace, and not the global
  namespace (so it is accessed by surprise.dump.dump)

VERSION 1.0.2
=============

Date: 04/01/17

Just a minor change so that README.md is converted to rst for better rendering
on PyPI.

VERSION 1.0.1
=============

Date: 02/01/17

Enhancements
------------

* Added the GridSearch feature, by Maher
* Added a 'clip' option to the predict() method
* Added NMF algorithm
* Added entry point for better command line usage.
* Added CoClustering algorithm.
* Added SlopeOne algorithm.
* Added Probabilistic Matrix Factorization as an option SVD
* Cythonized Baseline Computation

Other
-----

* Surprise is now a scikit!
* Changed license to BSD
* Six is now a dependency

VERSION 1.0.0
=============

Date: 22/11/16

* Changed name from recsys to surprise
* Improved printing of accuracy measures.
* Added version number.
* Rewrote the the __main__.py

VERSION 0.0.4
=============

Date: 15/11/16

Enhancements
------------

* Added notebooks for comparing and evaluating algorithm performances
* Better use of setup.py
* Added a min_support parameter to the similarity measures.
* Added a min_k parameter to the KNN algorithms.
* The similarity matrix and baselines are now returned.
* You can now train on a whole training set without test set.
* The estimate method can return a tuple with prediction details.
* Added SVD and SVD++ algorithms.
* Removed all the x/y vs user/item stuff. That was useless for most algorithms.


API Changes
-----------

* Removed the @property decorator for many iterators.
* It's now up to the algorithms to decide if they can or cannot make a
	prediction.

VERSION 0.0.3
=============

Date: 25/10/16

* Added support for Python 2
[![GitHub version](https://badge.fury.io/gh/nicolashug%2FSurprise.svg)](https://badge.fury.io/gh/nicolashug%2FSurprise)
[![Documentation Status](https://readthedocs.org/projects/surprise/badge/?version=stable)](http://surprise.readthedocs.io/en/stable/?badge=stable)
[![Build Status](https://travis-ci.org/NicolasHug/Surprise.svg?branch=master)](https://travis-ci.org/NicolasHug/Surprise)
[![python versions](https://img.shields.io/badge/python-2.7%2C%203.5%2C%203.6-blue.svg)](http://surpriselib.com)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02174/status.svg)](https://doi.org/10.21105/joss.02174)

[![logo](logo_black.svg)](http://surpriselib.com)

Overview
--------

[Surprise](http://surpriselib.com) is a Python
[scikit](https://www.scipy.org/scikits.html) for building and analyzing
recommender systems that deal with explicit rating data.

[Surprise](http://surpriselib.com) **was designed with the
following purposes in mind**:

- Give users perfect control over their experiments. To this end, a strong
  emphasis is laid on
  [documentation](http://surprise.readthedocs.io/en/stable/index.html), which we
  have tried to make as clear and precise as possible by pointing out every
  detail of the algorithms.
- Alleviate the pain of [Dataset
  handling](http://surprise.readthedocs.io/en/stable/getting_started.html#load-a-custom-dataset).
  Users can use both *built-in* datasets
  ([Movielens](http://grouplens.org/datasets/movielens/),
  [Jester](http://eigentaste.berkeley.edu/dataset/)), and their own *custom*
  datasets.
- Provide various ready-to-use [prediction
  algorithms](http://surprise.readthedocs.io/en/stable/prediction_algorithms_package.html)
  such as [baseline
  algorithms](http://surprise.readthedocs.io/en/stable/basic_algorithms.html),
  [neighborhood
  methods](http://surprise.readthedocs.io/en/stable/knn_inspired.html), matrix
  factorization-based (
  [SVD](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.SVD),
  [PMF](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#unbiased-note),
  [SVD++](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.SVDpp),
  [NMF](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.NMF)),
  and [many
  others](http://surprise.readthedocs.io/en/stable/prediction_algorithms_package.html).
  Also, various [similarity
  measures](http://surprise.readthedocs.io/en/stable/similarities.html)
  (cosine, MSD, pearson...) are built-in.
- Make it easy to implement [new algorithm
  ideas](http://surprise.readthedocs.io/en/stable/building_custom_algo.html).
- Provide tools to [evaluate](http://surprise.readthedocs.io/en/stable/model_selection.html),
  [analyse](http://nbviewer.jupyter.org/github/NicolasHug/Surprise/tree/master/examples/notebooks/KNNBasic_analysis.ipynb/)
  and
  [compare](http://nbviewer.jupyter.org/github/NicolasHug/Surprise/blob/master/examples/notebooks/Compare.ipynb)
  the algorithms' performance. Cross-validation procedures can be run very
  easily using powerful CV iterators (inspired by
  [scikit-learn](http://scikit-learn.org/) excellent tools), as well as
  [exhaustive search over a set of
  parameters](http://surprise.readthedocs.io/en/stable/getting_started.html#tune-algorithm-parameters-with-gridsearchcv).


The name *SurPRISE* (roughly :) ) stands for *Simple Python RecommendatIon
System Engine*.

Please note that surprise does not support implicit ratings or content-based
information.


Getting started, example
------------------------

Here is a simple example showing how you can (down)load a dataset, split it for
5-fold cross-validation, and compute the MAE and RMSE of the
[SVD](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.SVD)
algorithm.


```python
from surprise import SVD
from surprise import Dataset
from surprise.model_selection import cross_validate

# Load the movielens-100k dataset (download it if needed).
data = Dataset.load_builtin('ml-100k')

# Use the famous SVD algorithm.
algo = SVD()

# Run 5-fold cross-validation and print results.
cross_validate(algo, data, measures=['RMSE', 'MAE'], cv=5, verbose=True)
```

**Output**:

```
Evaluating RMSE, MAE of algorithm SVD on 5 split(s).

            Fold 1  Fold 2  Fold 3  Fold 4  Fold 5  Mean    Std
RMSE        0.9311  0.9370  0.9320  0.9317  0.9391  0.9342  0.0032
MAE         0.7350  0.7375  0.7341  0.7342  0.7375  0.7357  0.0015
Fit time    6.53    7.11    7.23    7.15    3.99    6.40    1.23
Test time   0.26    0.26    0.25    0.15    0.13    0.21    0.06
```

[Surprise](http://surpriselib.com) can do **much** more (e.g,
[GridSearchCV](http://surprise.readthedocs.io/en/stable/getting_started.html#tune-algorithm-parameters-with-gridsearchcv))!
You'll find [more usage
examples](http://surprise.readthedocs.io/en/stable/getting_started.html) in the
[documentation ](http://surprise.readthedocs.io/en/stable/index.html).


Benchmarks
----------

Here are the average RMSE, MAE and total execution time of various algorithms
(with their default parameters) on a 5-fold cross-validation procedure. The
datasets are the [Movielens](http://grouplens.org/datasets/movielens/) 100k and
1M datasets. The folds are the same for all the algorithms. All experiments are
run on a notebook with Intel Core i5 7th gen (2.5 GHz) and 8Go RAM.  The code
for generating these tables can be found in the [benchmark
example](https://github.com/NicolasHug/Surprise/tree/master/examples/benchmark.py).

| [Movielens 100k](http://grouplens.org/datasets/movielens/100k)                                                                         |   RMSE |   MAE | Time    |
|:---------------------------------------------------------------------------------------------------------------------------------------|-------:|------:|:--------|
| [SVD](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.SVD)      |  0.934 | 0.737 | 0:00:11 |
| [SVD++](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.SVDpp)  |  0.92  | 0.722 | 0:09:03 |
| [NMF](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.NMF)      |  0.963 | 0.758 | 0:00:15 |
| [Slope One](http://surprise.readthedocs.io/en/stable/slope_one.html#surprise.prediction_algorithms.slope_one.SlopeOne)                 |  0.946 | 0.743 | 0:00:08 |
| [k-NN](http://surprise.readthedocs.io/en/stable/knn_inspired.html#surprise.prediction_algorithms.knns.KNNBasic)                        |  0.98  | 0.774 | 0:00:10 |
| [Centered k-NN](http://surprise.readthedocs.io/en/stable/knn_inspired.html#surprise.prediction_algorithms.knns.KNNWithMeans)           |  0.951 | 0.749 | 0:00:10 |
| [k-NN Baseline](http://surprise.readthedocs.io/en/stable/knn_inspired.html#surprise.prediction_algorithms.knns.KNNBaseline)            |  0.931 | 0.733 | 0:00:12 |
| [Co-Clustering](http://surprise.readthedocs.io/en/stable/co_clustering.html#surprise.prediction_algorithms.co_clustering.CoClustering) |  0.963 | 0.753 | 0:00:03 |
| [Baseline](http://surprise.readthedocs.io/en/stable/basic_algorithms.html#surprise.prediction_algorithms.baseline_only.BaselineOnly)   |  0.944 | 0.748 | 0:00:01 |
| [Random](http://surprise.readthedocs.io/en/stable/basic_algorithms.html#surprise.prediction_algorithms.random_pred.NormalPredictor)    |  1.514 | 1.215 | 0:00:01 |


| [Movielens 1M](http://grouplens.org/datasets/movielens/1m)                                                                             |   RMSE |   MAE | Time    |
|:---------------------------------------------------------------------------------------------------------------------------------------|-------:|------:|:--------|
| [SVD](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.SVD)      |  0.873 | 0.686 | 0:02:13 |
| [SVD++](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.SVDpp)  |  0.862 | 0.673 | 2:54:19 |
| [NMF](http://surprise.readthedocs.io/en/stable/matrix_factorization.html#surprise.prediction_algorithms.matrix_factorization.NMF)      |  0.916 | 0.724 | 0:02:31 |
| [Slope One](http://surprise.readthedocs.io/en/stable/slope_one.html#surprise.prediction_algorithms.slope_one.SlopeOne)                 |  0.907 | 0.715 | 0:02:31 |
| [k-NN](http://surprise.readthedocs.io/en/stable/knn_inspired.html#surprise.prediction_algorithms.knns.KNNBasic)                        |  0.923 | 0.727 | 0:05:27 |
| [Centered k-NN](http://surprise.readthedocs.io/en/stable/knn_inspired.html#surprise.prediction_algorithms.knns.KNNWithMeans)           |  0.929 | 0.738 | 0:05:43 |
| [k-NN Baseline](http://surprise.readthedocs.io/en/stable/knn_inspired.html#surprise.prediction_algorithms.knns.KNNBaseline)            |  0.895 | 0.706 | 0:05:55 |
| [Co-Clustering](http://surprise.readthedocs.io/en/stable/co_clustering.html#surprise.prediction_algorithms.co_clustering.CoClustering) |  0.915 | 0.717 | 0:00:31 |
| [Baseline](http://surprise.readthedocs.io/en/stable/basic_algorithms.html#surprise.prediction_algorithms.baseline_only.BaselineOnly)   |  0.909 | 0.719 | 0:00:19 |
| [Random](http://surprise.readthedocs.io/en/stable/basic_algorithms.html#surprise.prediction_algorithms.random_pred.NormalPredictor)    |  1.504 | 1.206 | 0:00:19 |


Installation
------------

With pip (you'll need [numpy](http://www.numpy.org/), and a C compiler. Windows
users might prefer using conda):

    $ pip install numpy
    $ pip install scikit-surprise

With conda:

    $ conda install -c conda-forge scikit-surprise

For the latest version, you can also clone the repo and build the source
(you'll first need [Cython](http://cython.org/) and
[numpy](http://www.numpy.org/)):

    $ pip install numpy cython
    $ git clone https://github.com/NicolasHug/surprise.git
    $ cd surprise
    $ python setup.py install

License and reference
---------------------

This project is licensed under the [BSD
3-Clause](https://opensource.org/licenses/BSD-3-Clause) license, so it can be
used for pretty much everything, including commercial applications. Please let
us know how [Surprise](http://surpriselib.com) is useful to you!

Please make sure to cite the
[paper](https://joss.theoj.org/papers/10.21105/joss.02174) if you use
Surprise for your research:

    @article{Hug2020,
      doi = {10.21105/joss.02174},
      url = {https://doi.org/10.21105/joss.02174},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {52},
      pages = {2174},
      author = {Nicolas Hug},
      title = {Surprise: A Python library for recommender systems},
      journal = {Journal of Open Source Software}
    }

Contributors
------------

The following persons have contributed to [Surprise](http://surpriselib.com):

ashtou, bobbyinfj, caoyi, Олег Демиденко, Charles-Emmanuel Dias, dmamylin,
Lauriane Ducasse, Marc Feger, franckjay, Lukas Galke, Tim Gates,
Pierre-François Gimenez, Zachary Glassman, Jeff Hale, Nicolas Hug, Janniks,
jyesawtellrickson, Doruk Kilitcioglu, Ravi Raju Krishna, Hengji Liu, Maher
Malaeb, Manoj K, James McNeilis, Naturale0, nju-luke, Jay Qi, Lucas Rebscher,
Skywhat, David Stevens, TrWestdoor, Victor Wang, Mike Lee Williams, Jay Wong,
Chenchen Xu, YaoZh1918.

Thanks a lot :) !

Development Status
------------------

Starting from version 1.1.0 (September 19), we will only maintain the
package and provide bugfixes. No new features will be considered.

For bugs, issues or questions about [Surprise](http://surpriselib.com),
please use the GitHub [project page](https://github.com/NicolasHug/Surprise).
Please don't send emails (we will not answer).
Contributing to Surprise
========================

Disclamer: please note that starting from version 1.1.0, only bugfixes and
documentation improvements are considered. We will not accept new features.

Before submitting a new pull request, please make sure that:

* Your code is [clean](https://www.youtube.com/watch?v=wf-BqAjZb8M),
  [pythonic](https://www.youtube.com/watch?v=OSGv2VnC0go), well commented and
  also well documented (see below for building the docs).
* The tests are passing. Also, write some tests for the changes you're
  proposing. If you're not willing to write tests, it's best not to submit a PR
  (it's just a waste of time for everyone).
* Your code follows [PEP 8](https://www.python.org/dev/peps/pep-0008/) as much
  as possible. Coding style is automatically checked when tests are run. About
  line length: it's best to respect to 80 columns constraint, but tests will
  pass as long as the length is less than 88.
* For new prediction algorithms or similarity metrics, please submit a
  relevent benchmark outlining the performance of the new feature (in terms of
  accuracy, computation time, etc.). You can take a look at
  [`examples/benchmarks`](https://github.com/NicolasHug/Surprise/blob/master/examples/benchmark.py)
  for inspiration.

Set up
------

It's highly recommended to use a virtual environment. All the packages needed
for the development of Surprise (sphinx, flake8, etc...) can be installed by
running

    pip install -r requirements_dev.txt

Then, you can install your local copy of the repo:

    pip install -e .

Any change to the code should now be immediately reflected during execution. If
you're modifying Cython code (`.pyx` files), you'll need to compile the code in
order to see the chanes. This can be achieved by running `pip install -e .`
again.

Running and writing tests
-------------------------

Our testing tool is [pytest](http://doc.pytest.org/en/latest/). Running the tests is as
simple as running

    pytest

in the root directory.

For writing new tests, check out pytest getting started guide and / or take
inspiration from the current tests in the `tests` directory.


Building the docs locally
-------------------------

The docs can be compiled with

    cd doc
    make html

You can check the results in `doc/build/html`. Please make sure that the docs
compile without errors. Run `make clean` from time to time in order to avoid
hidden warnings. You can check spelling mistakes by running

    make spelling

Legit words that are not recognized can be added in the
`source/spelling_wordlist.txt` file.
TODO
====

* remove offset from everywhere (Reader, Trainset, test(), predict()...). Make
  a test to make sure that zero ratings are handled correctly. Right now we
  store ratings as defaultdict(list) in ur and ir, so there's no problem. Maybe
  it will change in the future, we would have to find a workaroud.
* rating_scale should be specified on dataset creation: load_from_file,
  load_from_folds, load_from_df. Deprecate its use from rating but fall back to
  it if it has only been specified here.
* Dataset.binarize should change rating_scale to something else. Maybe 'unary'?
* about clip: don't clip if rating_scale is 'unary'?
* test(), cross_validate() and GridSearch.fit() should probably allow a
  'predict_params' parameter?
* What about a 'dont_clip' attribute in algorithms that could take precedence
  over the clipping behaviour? What about just overriding predict() in the
  child classes?


* Document get_dataset_dir()...
* make some filtering dataset tools, like remove users/items with less/more
  than n ratings, binarize a dataset, etc...
* then implement MFBPR and see how it goes
* Allow incremental updates for some algorithms

Done:
-----

* Grid search now has the refit param.
* Grid search and cross_validate now allow return_train_score
* Make all fit methods return self. Update docs on building custom algorithms
* Update doc of MF algo to indicate how to retrieve latent factors.
* all algorithms using random initialization now have a random_state parameter.
* CV iterators:
  - Write basic CV iterators
  - evaluate -> rewrite to use CV iterators. Rename it into cross_validate.
  - Same for GridSearch. Keep it in a model_selection module like scikit-learn
    so that we can keep the old deprecated version. 
  - Make cross validation parallel with joblib
  - Add deprecation warnings for evaluate and GridSearch()
  - handle the cv_results attribute for grid search
  - (re)write all verbose settings for gridsearch and cross_validate
  - Change examples so they use CV iterators and the new gridsearch and
    cross_validate
  - indicate in docs that split(), folds(), evaluate() and gridsearch() are
    deprecated
  - Write comments, docstring and update all docs
  - Update main and command-line usage doc in getting started.rst
* Allow to change data folder from env variable
* Complete FAQ
* Change the dumping machinery to be more consistent 
* Allow to test on the trainset
* make bibtex entry
* Verbosity of gridsearch still prints stuff because of evaluate. Fix that.
* Make the r_ui param of predict optional
* Put some PredictionImpossible messages in every algo
* allow a 'clip' option to the predict method? Also, describe r_min and r_max
* configure entrypoints to use surprise directly from command line
* Allow a 'biased' option in the SVD algo. If true, use baselines, if False,
  don't. It should be pretty easy to do.
* create option in __main__ to clean the .recsys directory. Actually, the
  __main__ module should be entirely reviewed.
* when dumping, we should dump all the algorithm parameter. Use __dict__ ?
* do something about the generators Python 2 vs 3 (range, dict.items(), etc...)
* should a Prediction output the raw id or the inner id? Right now it's the
  inner id. Maybe sort this out when working on the comparison tools.
* allow the perf dict returned by evaluate to accept keys with lower/upper
  case for retarded users such as me.
* Add a 'min_support' parameter to sim_options? Add a min_k to knns?
* Do something about the user_based stuff. It should be better. Check knns BTW.
* Do something about unknown users and unknown items, i.e. users or items that
  have no rating in the trainset. Right now, the predict method checks if the
  name starts with 'unknown' but this is shiiite because it's dependent on the
  construct_trainset method, which is sometimes never called (so the raw2inner
  stuff will come in play somehow). Plus, It should be up to the algorithms to
  choose whether it can (or can't) make a prediction even if user or item is
  unknown.
* remove kwargs : done where useless.
* say something quick about baseline computation (when not matrix facto) 
* Matrix facto algo
* allow the 'estimate' method to return some details about prediction (such as
  the number of neighbors for a KNN)
* allow to train on a SINGLE file without test set, and let user query for some
  predictions
* write tuto for using only predict() (and not test)
* maybe clean a little all the dataset machinery? Plus, are the
  raw2inner_id_users and raw2inner_id_items worth keeping? May be for analysing
  tools, I don't know right now. EDIT: yes, we need to keep them, simply
  because the similarity computation can only work with integer as indexes
  (numpy arrays).
* sort out this warning issue coming from cython
* say something about the sim > 0 in knns algos
* get less restrictive requirements.txt
* write the custom algorithm tutorial
* improve test coverage
* add the cool stickers on the readme just like scikit learn
* set up travis
* keep on testing
* keep on documenting and commenting code
* extensively test the reader class, + check that the doc is OK for reader
* set up a nice API (looks ok now)
* handle algo-specific or similarity-specific parameters (such as 'k' for knn,
  regularization parameters, shrinkage paramaters, etc.) in an appropriate
  manner, rather than pass them all to constructors... UPDATE: ok so using
  kwargs like matplotlib.pyplot might be enough. should we create a
  'Similarity' class?
* clean the main and all the dataset handling stuff (still needs to be
  polished)
* rewrite this TODO in english
* create a proper project structure
* from camelCase to snake\_case
<!-- IMPORTANT PLEASE READ!!!
Before submitting an issue, make sure it hasn't already been addressed by
checking the past issues and the documentation (FAQ, getting started / advanced
usage guide, etc.). In order to let me help you more efficiently, please fill
the different fields below.

Also, please keep in mind that I develop and maintain this software on my
**free time**. So please, before asking for help, show me that you have already
tried to solve your problem for more than 5 seconds. That would be greatly
appreciated!

Thank you! :)
Nicolas
-->

#### Description
<!-- Describe your issue here. Don't forget to format the code as well as the
     error messages-->


#### Steps/Code to Reproduce
<!-- Please provide a **minimal** code example for reproduction. -->

#### Expected Results

#### Actual Results

#### Versions
<!-- Please run the following snippet and paste the output below.

import platform; print(platform.platform())
import sys; print("Python", sys.version)
import surprise; print("surprise", surprise.__version__)
-->
.. _FAQ:

FAQ
===

You will find here the Frequently Asked Questions, as well as some other
use-case examples that are not part of the User Guide.

How to get the top-N recommendations for each user
----------------------------------------------------------

Here is an example where we retrieve the top-10 items with highest
rating prediction for each user in the MovieLens-100k dataset. We first train
an SVD algorithm on the whole dataset, and then predict all the ratings for the
pairs (user, item) that are not in the training set. We then retrieve the
top-10 prediction for each user.

.. literalinclude:: ../../examples/top_n_recommendations.py
    :caption: From file ``examples/top_n_recommendations.py``
    :name: top_n_recommendations.py
    :lines: 10-

.. _precision_recall_at_k:

How to compute precision@k and recall@k
-----------------------------------------------------------------------

Here is an example where we compute Precision@k and Recall@k for each user:

:math:`\text{Precision@k} = \frac{ | \{ \text{Recommended items that are relevant} \} | }{ | \{ \text{Recommended items} \} | }`
:math:`\text{Recall@k} = \frac{ | \{ \text{Recommended items that are relevant} \} | }{ | \{ \text{Relevant items} \} | }`

An item is considered relevant if its true rating :math:`r_{ui}` is greater
than a given threshold.  An item is considered recommended if its estimated
rating :math:`\hat{r}_{ui}` is greater than the threshold, and if it is among
the k highest estimated ratings.

Note that in the edge cases where division by zero occurs, 
Precision@k and Recall@k values are undefined. 
As a convention, we set their values to 0 in such cases. 

.. literalinclude:: ../../examples/precision_recall_at_k.py
    :caption: From file ``examples/precision_recall_at_k.py``
    :name: precision_recall_at_k.py
    :lines: 7-

.. _get_k_nearest_neighbors:

How to get the k nearest neighbors of a user (or item)
--------------------------------------------------------------

You can use the :meth:`get_neighbors()
<surprise.prediction_algorithms.algo_base.AlgoBase.get_neighbors>` methods of
the algorithm object. This is only relevant for algorithms that use a
similarity measure, such as the :ref:`k-NN algorithms
<pred_package_knn_inpired>`.

Here is an example where we retrieve the 10 nearest neighbors of the movie Toy
Story from the MovieLens-100k dataset. The output is:

.. parsed-literal::

    The 10 nearest neighbors of Toy Story are:
    Beauty and the Beast (1991)
    Raiders of the Lost Ark (1981)
    That Thing You Do! (1996)
    Lion King, The (1994)
    Craft, The (1996)
    Liar Liar (1997)
    Aladdin (1992)
    Cool Hand Luke (1967)
    Winnie the Pooh and the Blustery Day (1968)
    Indiana Jones and the Last Crusade (1989)

There's a lot of boilerplate because of the conversions between movie names and
their raw/inner ids (see :ref:`this note <raw_inner_note>`), but it all boils
down to the use of :meth:`get_neighbors()
<surprise.prediction_algorithms.algo_base.AlgoBase.get_neighbors>`:

.. literalinclude:: ../../examples/k_nearest_neighbors.py
    :caption: From file ``examples/k_nearest_neighbors.py``
    :name: k_nearest_neighbors.py
    :lines: 10-

Naturally, the same can be done for users with minor modifications.

.. _serialize_an_algorithm:

How to serialize an algorithm
-----------------------------

Prediction algorithms can be serialized and loaded back using the :func:`dump()
<surprise.dump.dump>` and :func:`load() <surprise.dump.load>` functions. Here
is a small example where the SVD algorithm is trained on a dataset and
serialized. It is then reloaded and can be used again for making predictions:

.. literalinclude:: ../../examples/serialize_algorithm.py
    :caption: From file ``examples/serialize_algorithm.py``
    :name: serialize_algorithm.py
    :lines: 9-

.. _further_analysis:

Algorithms can be serialized along with their predictions, so that can be
further analyzed or compared with other algorithms, using pandas dataframes.
Some examples are given in the two following notebooks:

    * `Dumping and analysis of the KNNBasic algorithm
      <http://nbviewer.jupyter.org/github/NicolasHug/Surprise/tree/master/examples/notebooks/KNNBasic_analysis.ipynb/>`_.
    * `Comparison of two algorithms
      <http://nbviewer.jupyter.org/github/NicolasHug/Surprise/tree/master/examples/notebooks/Compare.ipynb/>`_.

How to build my own prediction algorithm
----------------------------------------

There's a whole guide :ref:`here<building_custom_algo>`.

.. _raw_inner_note:

What are raw and inner ids
--------------------------

Users and items have a raw id and an inner id. Some methods will use/return a
raw id (e.g. the :meth:`predict()
<surprise.prediction_algorithms.algo_base.AlgoBase.predict>` method), while
some other will use/return an inner id.

Raw ids are ids as defined in a rating file or in a pandas dataframe. They can
be strings or numbers. Note though that if the ratings were read from a file
which is the standard scenario, they are represented as strings. **This is
important to know if you're using e.g.** :meth:`predict()
<surprise.prediction_algorithms.algo_base.AlgoBase.predict>` **or other methods
that accept raw ids as parameters.**

On trainset creation, each raw id is mapped to a unique integer called inner
id, which is a lot more suitable for `Surprise
<https://nicolashug.github.io/Surprise/>`_ to manipulate. Conversions between
raw and inner ids can be done using the :meth:`to_inner_uid()
<surprise.Trainset.to_inner_uid>`, :meth:`to_inner_iid()
<surprise.Trainset.to_inner_iid>`, :meth:`to_raw_uid()
<surprise.Trainset.to_raw_uid>`, and :meth:`to_raw_iid()
<surprise.Trainset.to_raw_iid>` methods of the :class:`trainset
<surprise.Trainset>`.


Can I use my own dataset with Surprise, and can it be a pandas dataframe
------------------------------------------------------------------------

Yes, and yes. See the :ref:`user guide <load_custom>`.

How to tune an algorithm parameters
-----------------------------------

You can tune the parameters of an algorithm with the :class:`GridSearchCV
<surprise.model_selection.search.GridSearchCV>` class as described :ref:`here
<tuning_algorithm_parameters>`. After the tuning, you may want to have an
:ref:`unbiased estimate of your algorithm performances
<unbiased_estimate_after_tuning>`.

How to get accuracy measures on the training set
------------------------------------------------

You can use the :meth:`build_testset()
<surprise.Trainset.build_testset()>` method of the :class:`Trainset
<surprise.Trainset>` object to build a testset that can be then used
with the :meth:`test()
<surprise.prediction_algorithms.algo_base.AlgoBase.test>` method:

.. literalinclude:: ../../examples/evaluate_on_trainset.py
    :caption: From file ``examples/evaluate_on_trainset.py``
    :name: evaluate_on_trainset.py
    :lines: 9-25

Check out the example file for more usage examples.

.. _unbiased_estimate_after_tuning:

How to save some data for unbiased accuracy estimation
------------------------------------------------------

If your goal is to tune the parameters of an algorithm, you may want to spare a
bit of data to have an unbiased estimation of its performances. For instance
you may want to split your data into two sets A and B. A is used for parameter
tuning using grid search, and B is used for unbiased estimation. This can be
done as follows:

.. literalinclude:: ../../examples/split_data_for_unbiased_estimation.py
    :caption: From file ``examples/split_data_for_unbiased_estimation.py``
    :name: split_data_for_unbiased_estimation.py
    :lines: 10-

How to have reproducible experiments
------------------------------------

Some algorithms randomly initialize their parameters (sometimes with
``numpy``), and the cross-validation folds are also randomly generated. If you
need to reproduce your experiments multiple times, you just have to set the
seed of the RNG at the beginning of your program:

.. code::

    import random
    import numpy as np

    my_seed = 0
    random.seed(my_seed)
    np.random.seed(my_seed)

.. _data_folder:

Where are datasets stored and how to change it?
-----------------------------------------------

By default, datasets downloaded by Surprise will be saved in the
``'~/.surprise_data'`` directory. This is also where dump files will be stored.
You can change the default directory by setting the ``'SURPRISE_DATA_FOLDER'``
environment variable.

Can Surprise support content-based data or implicit ratings?
------------------------------------------------------------

No: this is out of scope for surprise. Surprise was designed for explicit
ratings.
.. _pred_package_basic_algorithms:

Basic algorithms
----------------

These are basic algorithms that do not do much work but that are still useful
for comparing accuracies.

.. autoclass:: surprise.prediction_algorithms.random_pred.NormalPredictor
    :show-inheritance:

.. autoclass:: surprise.prediction_algorithms.baseline_only.BaselineOnly
    :show-inheritance:

.. _pred_package_co_clustering:

Co-clustering
-------------

.. autoclass:: surprise.prediction_algorithms.co_clustering.CoClustering
    :show-inheritance:

.. _pred_package_matrix_factorization:

Matrix Factorization-based algorithms
-------------------------------------

.. autoclass:: surprise.prediction_algorithms.matrix_factorization.SVD
    :show-inheritance:

.. autoclass:: surprise.prediction_algorithms.matrix_factorization.SVDpp
    :show-inheritance:

.. autoclass:: surprise.prediction_algorithms.matrix_factorization.NMF
    :show-inheritance:
.. _pred_package_slope_one:

Slope One
---------

.. autoclass:: surprise.prediction_algorithms.slope_one.SlopeOne
    :show-inheritance:
.. _accuracy:

accuracy module
===================


.. automodule:: surprise.accuracy
    :members:
    :undoc-members:
    :show-inheritance:
.. _trainset:

Trainset class
==============

.. autoclass:: surprise.Trainset
    :members:
.. _similarities:

similarities module
===================

.. automodule:: surprise.similarities
    :members:
    :exclude-members: compute_mean_diff
    :show-inheritance:
.. _pred_package_predictions:

The predictions module
------------------------

.. automodule:: surprise.prediction_algorithms.predictions
    :members:
    :exclude-members: all_ratings, all_xs, all_ys

.. _prediction_algorithms_package:

prediction_algorithms package
=============================

.. automodule:: surprise.prediction_algorithms

You may want to check the :ref:`notation standards <notation_standards>`
before diving into the formulas.


.. toctree::
   :includehidden:

   algobase
   predictions_module
   basic_algorithms
   knn_inspired
   matrix_factorization
   slope_one
   co_clustering
.. _reader:

Reader class
============

.. autoclass:: surprise.reader.Reader
    :members:
    :exclude-members: parse_line

.. _pred_package_algo_base:

The algorithm base class
------------------------

.. automodule:: surprise.prediction_algorithms.algo_base
    :members:
.. _model_selection:

The model_selection package
---------------------------

Surprise provides various tools to run cross-validation procedures and search
the best parameters for a prediction algorithm. The tools presented here are
all heavily inspired from the excellent `scikit learn
<http://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection>`_
library.


.. _cross_validation_iterators_api:

Cross validation iterators
==========================

.. automodule:: surprise.model_selection.split
    :members:
    :exclude-members: get_cv, get_rng

Cross validation
================

.. autofunction:: surprise.model_selection.validation.cross_validate

Parameter search
================

.. autoclass:: surprise.model_selection.search.GridSearchCV
    :members:
    :inherited-members:

.. autoclass:: surprise.model_selection.search.RandomizedSearchCV
    :members:
    :inherited-members:

.. _getting_started:

Getting Started
===============


Basic usage
-----------

.. _cross_validate_example:

Automatic cross-validation
~~~~~~~~~~~~~~~~~~~~~~~~~~

`Surprise <https://nicolashug.github.io/Surprise/>`_ has a set of built-in
:ref:`algorithms<prediction_algorithms>` and :ref:`datasets <dataset>` for you
to play with. In its simplest form, it only takes a few lines of code to
run a cross-validation procedure:

.. literalinclude:: ../../examples/basic_usage.py
    :caption: From file ``examples/basic_usage.py``
    :name: basic_usage.py
    :lines: 9-

The result should be as follows (actual values may vary due to randomization):

.. parsed-literal::

    Evaluating RMSE, MAE of algorithm SVD on 5 split(s).

                Fold 1  Fold 2  Fold 3  Fold 4  Fold 5  Mean    Std
    RMSE        0.9311  0.9370  0.9320  0.9317  0.9391  0.9342  0.0032
    MAE         0.7350  0.7375  0.7341  0.7342  0.7375  0.7357  0.0015
    Fit time    6.53    7.11    7.23    7.15    3.99    6.40    1.23
    Test time   0.26    0.26    0.25    0.15    0.13    0.21    0.06


The :meth:`load_builtin() <surprise.dataset.Dataset.load_builtin>` method will
offer to download the `movielens-100k dataset
<http://grouplens.org/datasets/movielens/>`_ if it has not already been
downloaded, and it will save it in the ``.surprise_data`` folder in your home
directory (you can also choose to save it :ref:`somewhere else <data_folder>`).

We are here using the well-known
:class:`SVD<surprise.prediction_algorithms.matrix_factorization.SVD>`
algorithm, but many other algorithms are available. See
:ref:`prediction_algorithms` for more details.

The :func:`cross_validate()<surprise.model_selection.validation.cross_validate>`
function runs a cross-validation procedure according to the ``cv`` argument,
and computes some :mod:`accuracy <surprise.accuracy>` measures. We are here
using a classical 5-fold cross-validation, but fancier iterators can be used
(see :ref:`here <cross_validation_iterators_api>`).

Train-test split and the fit() method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _train_test_split_example:

If you don't want to run a full cross-validation procedure, you can use the
:func:`train_test_split() <surprise.model_selection.split.train_test_split>`
to sample a trainset and a testset with given sizes, and use the :mod:`accuracy
metric<surprise.accuracy>` of your chosing. You'll need to use the :meth:`fit()
<surprise.prediction_algorithms.algo_base.AlgoBase.fit>` method which will
train the algorithm on the trainset, and the :meth:`test()
<surprise.prediction_algorithms.algo_base.AlgoBase.test>` method which will
return the predictions made from the testset:

.. literalinclude:: ../../examples/train_test_split.py
    :caption: From file ``examples/train_test_split.py``
    :name: train_test_split.py
    :lines: 8-

Result:

.. parsed-literal::

    RMSE: 0.9411

Note that you can train and test an algorithm with the following one-line:

.. parsed-literal::

    predictions = algo.fit(trainset).test(testset)


In some cases, your trainset and testset are already defined by some files.
Please refer to :ref:`this section <load_from_folds_example>` to handle such cases.


.. _train_on_whole_trainset:

Train on a whole trainset and the predict() method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Obviously, we could also simply fit our algorithm to the whole dataset, rather
than running cross-validation. This can be done by using the
:meth:`build_full_trainset()
<surprise.dataset.DatasetAutoFolds.build_full_trainset>` method which will
build a :class:`trainset <surprise.Trainset>` object:

.. literalinclude:: ../../examples/predict_ratings.py
    :caption: From file ``examples/predict_ratings.py``
    :name: predict_ratings.py
    :lines: 9-20

We can now predict ratings by directly calling the :meth:`predict()
<surprise.prediction_algorithms.algo_base.AlgoBase.predict>` method.  Let's say
you're interested in user 196 and item 302 (make sure they're in the
trainset!), and you know that the true rating :math:`r_{ui} = 4`:

.. literalinclude:: ../../examples/predict_ratings.py
    :caption: From file ``examples/predict_ratings.py``
    :name: predict_ratings2.py
    :lines: 23-27

The result should be:

.. parsed-literal::

    user: 196        item: 302        r_ui = 4.00   est = 4.06   {'actual_k': 40, 'was_impossible': False}

.. note::

    The :meth:`predict()
    <surprise.prediction_algorithms.algo_base.AlgoBase.predict>` uses **raw**
    ids (please read :ref:`this <raw_inner_note>` about raw and inner ids). As
    the dataset we have used has been read from a file, the raw ids are strings
    (even if they represent numbers).

We have so far used a built-in dataset, but you can of course use your own.
This is explained in the next section.

.. _load_custom:

Use a custom dataset
--------------------

`Surprise <https://nicolashug.github.io/Surprise/>`_ has a set of  builtin
:ref:`datasets <dataset>`, but you can of course use a custom dataset.
Loading a rating dataset can be done either from a file (e.g. a csv file), or
from a pandas dataframe.  Either way, you will need to define a :class:`Reader
<surprise.reader.Reader>` object for `Surprise
<https://nicolashug.github.io/Surprise/>`_ to be able to parse the file or the
dataframe.

.. _load_from_file_example:

- To load a dataset from a file (e.g. a csv file), you will need the
  :meth:`load_from_file() <surprise.dataset.Dataset.load_from_file>` method:

  .. literalinclude:: ../../examples/load_custom_dataset.py
      :caption: From file ``examples/load_custom_dataset.py``
      :name: load_custom_dataset.py
      :lines: 12-28

  For more details about readers and how to use them, see the :class:`Reader
  class <surprise.reader.Reader>` documentation.

  .. note::
      As you already know from the previous section, the Movielens-100k dataset
      is built-in so a much quicker way to load the dataset is to do ``data =
      Dataset.load_builtin('ml-100k')``. We will of course ignore this here.

.. _load_from_df_example:

- To load a dataset from a pandas dataframe, you will need the
  :meth:`load_from_df() <surprise.dataset.Dataset.load_from_df>` method. You
  will also need a :class:`Reader<surprise.reader.Reader>` object, but only
  the ``rating_scale`` parameter must be specified. The dataframe must have
  three columns, corresponding to the user (raw) ids, the item (raw) ids, and
  the ratings in this order. Each row thus corresponds to a given rating. This
  is not restrictive as you can reorder the columns of your dataframe easily.

  .. literalinclude:: ../../examples/load_from_dataframe.py
      :caption: From file ``examples/load_from_dataframe.py``
      :name: load_dom_dataframe.py
      :lines: 8-29

  The dataframe initially looks like this:

  .. parsed-literal::

            itemID  rating    userID
      0       1       3         9
      1       1       2        32
      2       1       4         2
      3       2       3        45
      4       2       1  user_foo


.. _use_cross_validation_iterators:

Use cross-validation iterators
------------------------------

For cross-validation, we can use the :func:`cross_validate()
<surprise.model_selection.validation.cross_validate>` function that does all
the hard work for us. But for a better control, we can also instanciate a
cross-validation iterator, and make predictions over each split using the
``split()`` method of the iterator, and the
:meth:`test()<surprise.prediction_algorithms.algo_base.AlgoBase.test>` method
of the algorithm. Here is an example where we use a classical K-fold
cross-validation procedure with 3 splits:

.. literalinclude:: ../../examples/use_cross_validation_iterators.py
    :caption: From file ``examples/use_cross_validation_iterators.py``
    :name: use_cross_validation_iterators.py
    :lines: 8-

Result could be, e.g.:

.. parsed-literal::
    RMSE: 0.9374
    RMSE: 0.9476
    RMSE: 0.9478

Other cross-validation iterator can be used, like LeaveOneOut or ShuffleSplit.
See all the available iterators :ref:`here <cross_validation_iterators_api>`.
The design of Surprise's cross-validation tools is heavily inspired from the
excellent scikit-learn API.

---------------------

.. _load_from_folds_example:

A special case of cross-validation is when the folds are already predefined by
some files. For instance, the movielens-100K dataset already provides 5 train
and test files (u1.base, u1.test ... u5.base, u5.test). Surprise can handle
this case by using a :class:`surprise.model_selection.split.PredefinedKFold`
object:

.. literalinclude:: ../../examples/load_custom_dataset_predefined_folds.py
    :caption: From file ``examples/load_custom_dataset_predefined_folds.py``
    :name: load_custom_dataset_predefined_folds.py
    :lines: 13-

Of course, nothing prevents you from only loading a single file for training
and a single file for testing. However, the ``folds_files`` parameter still
needs to be a ``list``.

.. _tuning_algorithm_parameters:

Tune algorithm parameters with GridSearchCV
-------------------------------------------

The :func:`cross_validate()
<surprise.model_selection.validation.cross_validate>` function reports accuracy
metric over a cross-validation procedure for a given set of parameters.  If you
want to know which parameter combination yields the best results, the
:class:`GridSearchCV <surprise.model_selection.search.GridSearchCV>` class
comes to the rescue.  Given a ``dict`` of parameters, this class exhaustively
tries all the combinations of parameters and reports the best parameters for any
accuracy measure (averaged over the different splits). It is heavily inspired
from scikit-learn's `GridSearchCV
<http://scikit-learn.org/stable/modules/generated/sklearn.model
_selection.GridSearchCV.html>`_.

Here is an example where we try different values for parameters ``n_epochs``,
``lr_all`` and ``reg_all`` of the :class:`SVD
<surprise.prediction_algorithms.matrix_factorization.SVD>` algorithm.

.. literalinclude:: ../../examples/grid_search_usage.py
    :caption: From file ``examples/grid_search_usage.py``
    :name: grid_search_usage.py
    :lines: 9-26

Result:

.. parsed-literal::

    0.961300130118
    {'n_epochs': 10, 'lr_all': 0.005, 'reg_all': 0.4}

We are here evaluating the average RMSE and MAE over a 3-fold cross-validation
procedure, but any :ref:`cross-validation iterator
<cross_validation_iterators_api>` can used.

Once ``fit()`` has been called, the ``best_estimator`` attribute gives us an
algorithm instance with the optimal set of parameters, which can be used how we
please:

.. literalinclude:: ../../examples/grid_search_usage.py
    :caption: From file ``examples/grid_search_usage.py``
    :name: grid_search_usage2.py
    :lines: 28-30

.. _grid_search_note:
.. note::

    Dictionary parameters such as ``bsl_options`` and ``sim_options`` require
    particular treatment. See usage example below:

    .. parsed-literal::

        param_grid = {'k': [10, 20],
                      'sim_options': {'name': ['msd', 'cosine'],
                                      'min_support': [1, 5],
                                      'user_based': [False]}
                      }

    Naturally, both can be combined, for example for the
    :class:`KNNBaseline <surprise.prediction_algorithms.knns.KNNBaseline>`
    algorithm:

    .. parsed-literal::
        param_grid = {'bsl_options': {'method': ['als', 'sgd'],
                                      'reg': [1, 2]},
                      'k': [2, 3],
                      'sim_options': {'name': ['msd', 'cosine'],
                                      'min_support': [1, 5],
                                      'user_based': [False]}
                      }

.. _cv_results_example:

For further analysis, the ``cv_results`` attribute has all the needed
information and can be imported in a pandas dataframe:

.. literalinclude:: ../../examples/grid_search_usage.py
    :caption: From file ``examples/grid_search_usage.py``
    :name: grid_search_usage3.py
    :lines: 33

In our example, the ``cv_results`` attribute looks like this (floats are
formatted):

.. parsed-literal::

    'split0_test_rmse': [1.0, 1.0, 0.97, 0.98, 0.98, 0.99, 0.96, 0.97]
    'split1_test_rmse': [1.0, 1.0, 0.97, 0.98, 0.98, 0.99, 0.96, 0.97]
    'split2_test_rmse': [1.0, 1.0, 0.97, 0.98, 0.98, 0.99, 0.96, 0.97]
    'mean_test_rmse':   [1.0, 1.0, 0.97, 0.98, 0.98, 0.99, 0.96, 0.97]
    'std_test_rmse':    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    'rank_test_rmse':   [7 8 3 5 4 6 1 2]
    'split0_test_mae':  [0.81, 0.82, 0.78, 0.79, 0.79, 0.8, 0.77, 0.79]
    'split1_test_mae':  [0.8, 0.81, 0.78, 0.79, 0.78, 0.79, 0.77, 0.78]
    'split2_test_mae':  [0.81, 0.81, 0.78, 0.79, 0.78, 0.8, 0.77, 0.78]
    'mean_test_mae':    [0.81, 0.81, 0.78, 0.79, 0.79, 0.8, 0.77, 0.78]
    'std_test_mae':     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    'rank_test_mae':    [7 8 2 5 4 6 1 3]
    'mean_fit_time':    [1.53, 1.52, 1.53, 1.53, 3.04, 3.05, 3.06, 3.02]
    'std_fit_time':     [0.03, 0.04, 0.0, 0.01, 0.04, 0.01, 0.06, 0.01]
    'mean_test_time':   [0.46, 0.45, 0.44, 0.44, 0.47, 0.49, 0.46, 0.34]
    'std_test_time':    [0.0, 0.01, 0.01, 0.0, 0.03, 0.06, 0.01, 0.08]
    'params':           [{'n_epochs': 5, 'lr_all': 0.002, 'reg_all': 0.4}, {'n_epochs': 5, 'lr_all': 0.002, 'reg_all': 0.6}, {'n_epochs': 5, 'lr_all': 0.005, 'reg_all': 0.4}, {'n_epochs': 5, 'lr_all': 0.005, 'reg_all': 0.6}, {'n_epochs': 10, 'lr_all': 0.002, 'reg_all': 0.4}, {'n_epochs': 10, 'lr_all': 0.002, 'reg_all': 0.6}, {'n_epochs': 10, 'lr_all': 0.005, 'reg_all': 0.4}, {'n_epochs': 10, 'lr_all': 0.005, 'reg_all': 0.6}]
    'param_n_epochs':   [5, 5, 5, 5, 10, 10, 10, 10]
    'param_lr_all':     [0.0, 0.0, 0.01, 0.01, 0.0, 0.0, 0.01, 0.01]
    'param_reg_all':    [0.4, 0.6, 0.4, 0.6, 0.4, 0.6, 0.4, 0.6]

As you can see, each list has the same size of the number of parameter
combination. It corresponds to the following table:

==================  ==================  ==================  ================  ===============  ================  =================  =================  =================  ===============  ==============  ===============  ===============  ==============  ================  ===============  =================================================  ================  ==============  ===============
  split0_test_rmse    split1_test_rmse    split2_test_rmse    mean_test_rmse    std_test_rmse    rank_test_rmse    split0_test_mae    split1_test_mae    split2_test_mae    mean_test_mae    std_test_mae    rank_test_mae    mean_fit_time    std_fit_time    mean_test_time    std_test_time  params                                               param_n_epochs    param_lr_all    param_reg_all
==================  ==================  ==================  ================  ===============  ================  =================  =================  =================  ===============  ==============  ===============  ===============  ==============  ================  ===============  =================================================  ================  ==============  ===============
          0.99775             0.997744            0.996378          0.997291      0.000645508                 7           0.807862           0.804626           0.805282         0.805923      0.00139657                7          1.53341      0.0305216           0.455831      0.000922113  {'n_epochs': 5, 'lr_all': 0.002, 'reg_all': 0.4}                  5           0.002              0.4
          1.00381             1.00304             1.00257           1.00314       0.000508358                 8           0.816559           0.812905           0.813772         0.814412      0.00155866                8          1.5199       0.0367117           0.451068      0.00938646   {'n_epochs': 5, 'lr_all': 0.002, 'reg_all': 0.6}                  5           0.002              0.6
          0.973524            0.973595            0.972495          0.973205      0.000502609                 3           0.783361           0.780242           0.78067          0.781424      0.00138049                2          1.53449      0.00496203          0.441558      0.00529696   {'n_epochs': 5, 'lr_all': 0.005, 'reg_all': 0.4}                  5           0.005              0.4
          0.98229             0.982059            0.981486          0.981945      0.000338056                 5           0.794481           0.790781           0.79186          0.792374      0.00155377                5          1.52739      0.00859185          0.44463       0.000888907  {'n_epochs': 5, 'lr_all': 0.005, 'reg_all': 0.6}                  5           0.005              0.6
          0.978034            0.978407            0.976919          0.977787      0.000632049                 4           0.787643           0.784723           0.784957         0.785774      0.00132486                4          3.03572      0.0431101           0.466606      0.0254965    {'n_epochs': 10, 'lr_all': 0.002, 'reg_all': 0.4}                10           0.002              0.4
          0.986263            0.985817            0.985004          0.985695      0.000520899                 6           0.798218           0.794457           0.795373         0.796016      0.00160135                6          3.0544       0.00636185          0.488357      0.0576194    {'n_epochs': 10, 'lr_all': 0.002, 'reg_all': 0.6}                10           0.002              0.6
          0.963751            0.963463            0.962676          0.963297      0.000454661                 1           0.774036           0.770548           0.771588         0.772057      0.00146201                1          3.0636       0.0597982           0.456484      0.00510321   {'n_epochs': 10, 'lr_all': 0.005, 'reg_all': 0.4}                10           0.005              0.4
          0.973605            0.972868            0.972765          0.973079      0.000374222                 2           0.78607            0.781918           0.783537         0.783842      0.00170855                3          3.01907      0.011834            0.338839      0.075346     {'n_epochs': 10, 'lr_all': 0.005, 'reg_all': 0.6}                10           0.005              0.6
==================  ==================  ==================  ================  ===============  ================  =================  =================  =================  ===============  ==============  ===============  ===============  ==============  ================  ===============  =================================================  ================  ==============  ===============



Command line usage
------------------

Surprise can also be used from the command line, for example:

.. code::

    surprise -algo SVD -params "{'n_epochs': 5, 'verbose': True}" -load-builtin ml-100k -n-folds 3

See detailed usage by running:

.. code::

    surprise -h
.. _prediction_algorithms:

Using prediction algorithms
===========================

Surprise provides a bunch of built-in algorithms. All algorithms derive from
the :class:`AlgoBase <surprise.prediction_algorithms.algo_base.AlgoBase>` base
class, where are implemented some key methods (e.g. :meth:`predict
<surprise.prediction_algorithms.algo_base.AlgoBase.predict>`, :meth:`fit
<surprise.prediction_algorithms.algo_base.AlgoBase.fit>` and :meth:`test
<surprise.prediction_algorithms.algo_base.AlgoBase.test>`). The list and
details of the available prediction algorithms can be found in the
:mod:`prediction_algorithms <surprise.prediction_algorithms>` package
documentation.

Every algorithm is part of the global Surprise namespace, so you only need to
import their names from the Surprise package, for example: ::

    from surprise import KNNBasic
    algo = KNNBasic()


Some of these algorithms may use :ref:`baseline estimates
<baseline_estimates_configuration>`, some may use a :ref:`similarity measure
<similarity_measures_configuration>`. We will here review how to configure the
way baselines and similarities are computed.


.. _baseline_estimates_configuration:

Baselines estimates configuration
---------------------------------


.. note::
  This section only applies to algorithms (or similarity measures) that try to
  minimize the following regularized squared error (or equivalent):

  .. math::
    \sum_{r_{ui} \in R_{train}} \left(r_{ui} - (\mu + b_u + b_i)\right)^2 +
    \lambda \left(b_u^2 + b_i^2 \right).

  For algorithms using baselines in another objective function (e.g. the
  :class:`SVD <surprise.prediction_algorithms.matrix_factorization.SVD>`
  algorithm), the baseline configuration is done differently and is specific to
  each algorithm. Please refer to their own documentation.

First of all, if you do not want to configure the way baselines are computed,
you don't have to: the default parameters will do just fine. If you do want to
well... This is for you.

You may want to read section 2.1 of :cite:`Koren:2010` to get a good idea of
what are baseline estimates.

Baselines can be estimated in two different ways:

* Using Stochastic Gradient Descent (SGD).
* Using Alternating Least Squares (ALS).

You can configure the way baselines are computed using the ``bsl_options``
parameter passed at the creation of an algorithm. This parameter is a
dictionary for which the key ``'method'`` indicates the method to use. Accepted
values are ``'als'`` (default) and ``'sgd'``. Depending on its value, other
options may be set. For ALS:

- ``'reg_i'``: The regularization parameter for items. Corresponding to
  :math:`\lambda_2` in :cite:`Koren:2010`.  Default is ``10``.
- ``'reg_u'``: The regularization parameter for users. Corresponding to
  :math:`\lambda_3` in :cite:`Koren:2010`.  Default is ``15``.
- ``'n_epochs'``: The number of iteration of the ALS procedure. Default is
  ``10``.  Note that in :cite:`Koren:2010`, what is described is a **single**
  iteration ALS process.

And for SGD:

- ``'reg'``: The regularization parameter of the cost function that is
  optimized, corresponding to :math:`\lambda_1` in
  :cite:`Koren:2010`. Default is ``0.02``.
- ``'learning_rate'``: The learning rate of SGD, corresponding to
  :math:`\gamma` in :cite:`Koren:2010`.  Default is ``0.005``.
- ``'n_epochs'``: The number of iteration of the SGD procedure. Default is 20. 

.. note::
  For both procedures (ALS and SGD), user and item biases (:math:`b_u` and
  :math:`b_i`) are initialized to zero.

Usage examples:

.. literalinclude:: ../../examples/baselines_conf.py
    :caption: From file ``examples/baselines_conf.py``
    :name: baselines_als
    :lines: 19-25

.. literalinclude:: ../../examples/baselines_conf.py
    :caption: From file ``examples/baselines_conf.py``
    :name: baselines_sgd
    :lines: 30-34

Note that some similarity measures may use baselines, such as the
:func:`pearson_baseline <surprise.similarities.pearson_baseline>` similarity.
Configuration works just the same, whether the baselines are used in the actual
prediction :math:`\hat{r}_{ui}` or not:

.. literalinclude:: ../../examples/baselines_conf.py
    :caption: From file ``examples/baselines_conf.py``
    :name: baselines_als_pearson_sim
    :lines: 40-44


This leads us to similarity measure configuration, which we will review right
now.

.. _similarity_measures_configuration:

Similarity measure configuration
--------------------------------

Many algorithms use a similarity measure to estimate a rating. The way they can
be configured is done in a similar fashion as for baseline ratings: you just
need to pass a ``sim_options`` argument at the creation of an algorithm. This
argument is a dictionary with the following (all optional) keys:

- ``'name'``: The name of the similarity to use, as defined in the
  :mod:`similarities <surprise.similarities>` module. Default is ``'MSD'``.
- ``'user_based'``: Whether similarities will be computed between users or
  between items. This has a **huge** impact on the performance of a prediction
  algorithm.  Default is ``True``.
- ``'min_support'``: The minimum number of common items (when ``'user_based'``
  is ``'True'``) or minimum number of common users (when ``'user_based'`` is
  ``'False'``) for the similarity not to be zero. Simply put, if
  :math:`|I_{uv}| < \text{min_support}` then :math:`\text{sim}(u, v) = 0`. The
  same goes for items.
- ``'shrinkage'``: Shrinkage parameter to apply (only relevant for
  :func:`pearson_baseline <surprise.similarities.pearson_baseline>` similarity).
  Default is 100.

Usage examples:

.. literalinclude:: ../../examples/similarity_conf.py
    :caption: From file ``examples/similarity_conf.py``
    :name: sim_conf_cos
    :lines: 18-21

.. literalinclude:: ../../examples/similarity_conf.py
    :caption: From file ``examples/similarity_conf.py``
    :name: sim_conf_pearson_baseline
    :lines: 26-29

.. seealso::
    The :mod:`similarities <surprise.similarities>` module.
.. _dataset:

dataset module
===================

.. automodule:: surprise.dataset
    :members:
    :exclude-members: BuiltinDataset, read_ratings, DatasetUserFolds,
.. Surprise documentation master file, created by
   sphinx-quickstart on Tue Dec 29 20:08:18 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index:

Welcome to Surprise' documentation!
===================================

`Surprise <http://surpriselib.com>`_  is an easy-to-use Python `scikit
<https://www.scipy.org/scikits.html>`_ for recommender systems.

If you're new to `Surprise <http://surpriselib.com>`_, we invite you to take a
look at the :ref:`getting_started` guide, where you'll find a series of
tutorials illustrating all you can do with  `Surprise
<http://surpriselib.com>`_. You can also check out the :ref:`FAQ` for many
use-case example. For installation guidelines, please refer to the `project
page <http://surpriselib.com>`_.

Any kind of feedback/criticism would be greatly appreciated (software design,
documentation, improvement ideas, spelling mistakes, etc...). Please feel free
to contribute and send pull requests (see `GitHub page
<https://github.com/NicolasHug/Surprise>`_)!


.. toctree::
   :caption: User Guide
   :hidden:

   getting_started
   prediction_algorithms
   building_custom_algo
   notation_standards
   FAQ


.. toctree::
   :maxdepth: 2
   :caption: API Reference
   :hidden:

   prediction_algorithms_package
   model_selection
   similarities
   accuracy
   dataset
   trainset
   reader
   dump
.. _pred_package_knn_inpired:

k-NN inspired algorithms
------------------------

These are algorithms that are directly derived from a basic nearest neighbors
approach.

.. _actual_k_note:

.. note::

  For each of these algorithms, the actual number of neighbors that are
  aggregated to compute an estimation is necessarily less than or equal to
  :math:`k`. First, there might just not exist enough neighbors and second, the
  sets :math:`N_i^k(u)` and :math:`N_u^k(i)` only include neighbors for which
  the similarity measure is **positive**. It would make no sense to aggregate
  ratings from users (or items) that are negatively correlated. For a given
  prediction, the actual number of neighbors can be retrieved in the
  ``'actual_k'`` field of the ``details`` dictionary of the :class:`prediction
  <surprise.prediction_algorithms.predictions.Prediction>`.

You may want to read the :ref:`User Guide <similarity_measures_configuration>`
on how to configure the ``sim_options`` parameter.

.. autoclass:: surprise.prediction_algorithms.knns.KNNBasic
    :show-inheritance:

.. autoclass:: surprise.prediction_algorithms.knns.KNNWithMeans
    :show-inheritance:

.. autoclass:: surprise.prediction_algorithms.knns.KNNWithZScore
    :show-inheritance:

.. autoclass:: surprise.prediction_algorithms.knns.KNNBaseline
    :show-inheritance:
.. _building_custom_algo:

How to build your own prediction algorithm
==========================================

This page describes how to build a custom prediction algorithm using Surprise.

The basics
~~~~~~~~~~

Want to get your hands dirty? Cool.

Creating your own prediction algorithm is pretty simple: an algorithm is
nothing but a class derived from :class:`AlgoBase
<surprise.prediction_algorithms.algo_base.AlgoBase>` that has an ``estimate``
method.  This is the method that is called by the :meth:`predict()
<surprise.prediction_algorithms.algo_base.AlgoBase.predict>` method. It takes
in an **inner** user id, an **inner** item id (see :ref:`this note
<raw_inner_note>`), and returns the estimated rating :math:`\hat{r}_{ui}`:

.. literalinclude:: ../../examples/building_custom_algorithms/most_basic_algorithm.py
    :caption: From file ``examples/building_custom_algorithms/most_basic_algorithm.py``
    :name: most_basic_algorithm.py
    :lines: 9-

This algorithm is the dumbest we could have thought of: it just predicts a
rating of 3, regardless of users and items.

If you want to store additional information about the prediction, you can also
return a dictionary with given details: ::

    def estimate(self, u, i):

        details = {'info1' : 'That was',
                   'info2' : 'easy stuff :)'}
        return 3, details

This dictionary will be stored in the :class:`prediction
<surprise.prediction_algorithms.predictions.Prediction>` as the ``details``
field and can be used for :ref:`later analysis <further_analysis>`.



The ``fit`` method
~~~~~~~~~~~~~~~~~~~~

Now, let's make a slightly cleverer algorithm that predicts the average of all
the ratings of the trainset. As this is a constant value that does not depend
on current user or item, we would rather compute it once and for all. This can
be done by defining the ``fit`` method:

.. literalinclude:: ../../examples/building_custom_algorithms/most_basic_algorithm2.py
    :caption: From file ``examples/building_custom_algorithms/most_basic_algorithm2.py``
    :name: most_basic_algorithm2.py
    :lines: 16-37


The ``fit`` method is called e.g. by the :func:`cross_validate
<surprise.model_selection.validation.cross_validate>` function at each fold of
a cross-validation process, (but you can also :ref:`call it yourself
<use_cross_validation_iterators>`).  Before doing anything, you should call the
base class :meth:`fit()
<surprise.prediction_algorithms.algo_base.AlgoBase.fit>` method.

Note that the ``fit()`` method returns ``self``. This allows to use expression
like ``algo.fit(trainset).test(testset)``.

The ``trainset`` attribute
~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the base class :meth:`fit()
<surprise.prediction_algorithms.algo_base.AlgoBase.fit>` method has returned,
all the info you need about the current training set (rating values, etc...) is
stored in the ``self.trainset`` attribute. This is a :class:`Trainset
<surprise.Trainset>` object that has many attributes and methods of
interest for prediction.

To illustrate its usage, let's make an algorithm that predicts an average
between the mean of all ratings, the mean rating of the user and the mean
rating for the item:

.. literalinclude:: ../../examples/building_custom_algorithms/mean_rating_user_item.py
    :caption: From file ``examples/building_custom_algorithms/mean_rating_user_item.py``
    :name: mean_rating_user_item.py
    :lines: 23-35

Note that it would have been a better idea to compute all the user means in the
``fit`` method, thus avoiding the same computations multiple times.


When the prediction is impossible
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It's up to your algorithm to decide if it can or cannot yield a prediction. If
the prediction is impossible, then you can raise the
:class:`PredictionImpossible
<surprise.prediction_algorithms.predictions.PredictionImpossible>` exception.
You'll need to import it first: ::

    from surprise import PredictionImpossible


This exception will be caught by the :meth:`predict()
<surprise.prediction_algorithms.algo_base.AlgoBase.predict>` method, and the
estimation :math:`\hat{r}_{ui}` will be set according to
the :meth:`default_prediction()
<surprise.prediction_algorithms.algo_base.AlgoBase.default_prediction>` method,
which can be overridden. By default, it returns the average of all ratings in
the trainset.

Using similarities and baselines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Should your algorithm use a similarity measure or baseline estimates, you'll
need to accept ``bsl_options`` and ``sim_options`` as parameters to the
``__init__`` method, and pass them along to the Base class. See how to use
these parameters in the :ref:`prediction_algorithms` section.

Methods :meth:`compute_baselines()
<surprise.prediction_algorithms.algo_base.AlgoBase.compute_baselines>`   and
:meth:`compute_similarities()
<surprise.prediction_algorithms.algo_base.AlgoBase.compute_similarities>` can
be called in the ``fit`` method (or anywhere else).

.. literalinclude:: ../../examples/building_custom_algorithms/with_baselines_or_sim.py
    :caption: From file ``examples/building_custom_algorithms/.with_baselines_or_sim.py``
    :name: with_baselines_or_sim.py
    :lines: 15-47


Feel free to explore the prediction_algorithms package `source
<https://github.com/NicolasHug/Surprise/tree/master/surprise/prediction_algorithms>`_
to get an idea of what can be done.
.. _dump_module:

dump module
===============

.. automodule:: surprise.dump
    :members:
.. _notation_standards:

Notation standards, References
==============================

In the documentation, you will find the following notation:

* :math:`R` : the set of all ratings.
* :math:`R_{train}`, :math:`R_{test}` and :math:`\hat{R}` denote the training
  set, the test set, and the set of predicted ratings.
* :math:`U` : the set of all users. :math:`u` and :math:`v` denotes users.
* :math:`I` : the set of all items. :math:`i` and :math:`j` denotes items.
* :math:`U_i` : the set of all users that have rated item :math:`i`.
* :math:`U_{ij}` : the set of all users that have rated both items :math:`i`
  and :math:`j`.
* :math:`I_u` : the set of all items rated by user :math:`u`.
* :math:`I_{uv}` : the set of all items rated by both users :math:`u`
  and :math:`v`.
* :math:`r_{ui}` : the *true* rating of user :math:`u` for item
  :math:`i`.
* :math:`\hat{r}_{ui}` : the *estimated* rating of user :math:`u` for item
  :math:`i`.
* :math:`b_{ui}` : the baseline rating of user :math:`u` for item :math:`i`.
* :math:`\mu` : the mean of all ratings.
* :math:`\mu_u` : the mean of all ratings given by user :math:`u`.
* :math:`\mu_i` : the mean of all ratings given to item :math:`i`.
* :math:`\sigma_u` : the standard deviation of all ratings given by user :math:`u`.
* :math:`\sigma_i` : the standard deviation of all ratings given to item :math:`i`.
* :math:`N_i^k(u)` : the :math:`k` nearest neighbors of user :math:`u` that
  have rated item :math:`i`. This set is computed using a :mod:`similarity
  metric <surprise.similarities>`.
* :math:`N_u^k(i)` : the :math:`k` nearest neighbors of item :math:`i` that
  are rated by user :math:`u`. This set is computed using a :py:mod:`similarity
  metric <surprise.similarities>`.

.. rubric:: References

Here are the papers used as references in the documentation. Links to pdf files
where added when possible. A simple Google search should lead you easily to the
missing ones :)

.. bibliography:: refs.bib
  :all:
