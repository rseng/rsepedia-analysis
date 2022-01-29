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
