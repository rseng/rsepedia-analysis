Benchmarking nearest neighbors
==============================

[![Build Status](https://img.shields.io/github/workflow/status/erikbern/ann-benchmarks/ANN%20benchmarks?style=flat-square)](https://github.com/erikbern/ann-benchmarks/actions?query=workflow:benchmarks)

Doing fast searching of nearest neighbors in high dimensional spaces is an increasingly important problem, but so far there has not been a lot of empirical attempts at comparing approaches in an objective way.

This project contains some tools to benchmark various implementations of approximate nearest neighbor (ANN) search for different metrics. We have pregenerated datasets (in HDF5) formats and we also have Docker containers for each algorithm. There's a [test suite](https://travis-ci.org/erikbern/ann-benchmarks) that makes sure every algorithm works.

Evaluated
=========

* [Annoy](https://github.com/spotify/annoy)
* [FLANN](http://www.cs.ubc.ca/research/flann/)
* [scikit-learn](http://scikit-learn.org/stable/modules/neighbors.html): LSHForest, KDTree, BallTree
* [PANNS](https://github.com/ryanrhymes/panns)
* [NearPy](http://pixelogik.github.io/NearPy/)
* [KGraph](https://github.com/aaalgo/kgraph)
* [NMSLIB (Non-Metric Space Library)](https://github.com/nmslib/nmslib): SWGraph, HNSW, BallTree, MPLSH
* [hnswlib (a part of nmslib project)](https://github.com/nmslib/hnsw)
* [RPForest](https://github.com/lyst/rpforest)
* [FAISS](https://github.com/facebookresearch/faiss.git)
* [DolphinnPy](https://github.com/ipsarros/DolphinnPy)
* [Datasketch](https://github.com/ekzhu/datasketch)
* [PyNNDescent](https://github.com/lmcinnes/pynndescent)
* [MRPT](https://github.com/teemupitkanen/mrpt)
* [NGT](https://github.com/yahoojapan/NGT): ONNG, PANNG, QG
* [SPTAG](https://github.com/microsoft/SPTAG)
* [PUFFINN](https://github.com/puffinn/puffinn)
* [N2](https://github.com/kakao/n2)
* [ScaNN](https://github.com/google-research/google-research/tree/master/scann)
* [Elastiknn](https://github.com/alexklibisz/elastiknn)
* [OpenSearch KNN](https://github.com/opensearch-project/k-NN)
* [DiskANN](https://github.com/microsoft/diskann): Vamana, Vamana-PQ
* [Vespa](https://github.com/vespa-engine/vespa)
* [scipy](https://docs.scipy.org/doc/scipy/reference/spatial.html): cKDTree
* [vald](https://github.com/vdaas/vald)

Data sets
=========

We have a number of precomputed data sets for this. All data sets are pre-split into train/test and come with ground truth data in the form of the top 100 neighbors. We store them in a HDF5 format:

| Dataset                                                           | Dimensions | Train size | Test size | Neighbors | Distance  | Download                                                                   |
| ----------------------------------------------------------------- | ---------: | ---------: | --------: | --------: | --------- | -------------------------------------------------------------------------- |
| [DEEP1B](http://sites.skoltech.ru/compvision/noimi/)              |         96 |  9,990,000 |    10,000 |       100 | Angular   | [HDF5](http://ann-benchmarks.com/deep-image-96-angular.hdf5) (3.6GB)
| [Fashion-MNIST](https://github.com/zalandoresearch/fashion-mnist) |        784 |     60,000 |    10,000 |       100 | Euclidean | [HDF5](http://ann-benchmarks.com/fashion-mnist-784-euclidean.hdf5) (217MB) |
| [GIST](http://corpus-texmex.irisa.fr/)                            |        960 |  1,000,000 |     1,000 |       100 | Euclidean | [HDF5](http://ann-benchmarks.com/gist-960-euclidean.hdf5) (3.6GB)          |
| [GloVe](http://nlp.stanford.edu/projects/glove/)                  |         25 |  1,183,514 |    10,000 |       100 | Angular   | [HDF5](http://ann-benchmarks.com/glove-25-angular.hdf5) (121MB)            |
| GloVe                                                             |         50 |  1,183,514 |    10,000 |       100 | Angular   | [HDF5](http://ann-benchmarks.com/glove-50-angular.hdf5) (235MB)            |
| GloVe                                                             |        100 |  1,183,514 |    10,000 |       100 | Angular   | [HDF5](http://ann-benchmarks.com/glove-100-angular.hdf5) (463MB)           |
| GloVe                                                             |        200 |  1,183,514 |    10,000 |       100 | Angular   | [HDF5](http://ann-benchmarks.com/glove-200-angular.hdf5) (918MB)           |
| [Kosarak](http://fimi.uantwerpen.be/data/)                        |      27983 |     74,962 |       500 |       100 | Jaccard   | [HDF5](http://ann-benchmarks.com/kosarak-jaccard.hdf5) (2.0GB)             |
| [MNIST](http://yann.lecun.com/exdb/mnist/)                        |        784 |     60,000 |    10,000 |       100 | Euclidean | [HDF5](http://ann-benchmarks.com/mnist-784-euclidean.hdf5) (217MB)         |
| [NYTimes](https://archive.ics.uci.edu/ml/datasets/bag+of+words)   |        256 |    290,000 |    10,000 |       100 | Angular   | [HDF5](http://ann-benchmarks.com/nytimes-256-angular.hdf5) (301MB)         |
| [SIFT](http://corpus-texmex.irisa.fr/)                           |        128 |  1,000,000 |    10,000 |       100 | Euclidean | [HDF5](http://ann-benchmarks.com/sift-128-euclidean.hdf5) (501MB)          |
| [Last.fm](https://github.com/erikbern/ann-benchmarks/pull/91)     |         65 |    292,385 |    50,000 |       100 | Angular   | [HDF5](http://ann-benchmarks.com/lastfm-64-dot.hdf5) (135MB)               |

Results
=======

Interactive plots can be found at <http://ann-benchmarks.com>. These are all as of December 2021, running all benchmarks on a r5.4xlarge machine on AWS with `--parallelism 7`:

glove-100-angular
-----------------

![glove-100-angular](https://raw.github.com/erikbern/ann-benchmarks/master/results/glove-100-angular.png)

sift-128-euclidean
------------------

![glove-100-angular](https://raw.github.com/erikbern/ann-benchmarks/master/results/sift-128-euclidean.png)

fashion-mnist-784-euclidean
---------------------------

![fashion-mnist-784-euclidean](https://raw.github.com/erikbern/ann-benchmarks/master/results/fashion-mnist-784-euclidean.png)

lastfm-64-dot
------------------

![lastfm-64-dot](https://raw.github.com/erikbern/ann-benchmarks/master/results/lastfm-64-dot.png)

nytimes-256-angular
-------------------

![nytimes-256-angular](https://raw.github.com/erikbern/ann-benchmarks/master/results/nytimes-256-angular.png)

glove-25-angular
----------------

![glove-25-angular](https://raw.github.com/erikbern/ann-benchmarks/master/results/glove-25-angular.png)

Install
=======

The only prerequisite is Python (tested with 3.6) and Docker.

1. Clone the repo.
2. Run `pip install -r requirements.txt`.
3. Run `python install.py` to build all the libraries inside Docker containers (this can take a while, like 10-30 minutes).

Running
=======

1. Run `python run.py` (this can take an extremely long time, potentially days)
2. Run `python plot.py` or `python create_website.py` to plot results.

You can customize the algorithms and datasets if you want to:

* Check that `algos.yaml` contains the parameter settings that you want to test
* To run experiments on SIFT, invoke `python run.py --dataset glove-100-angular`. See `python run.py --help` for more information on possible settings. Note that experiments can take a long time. 
* To process the results, either use `python plot.py --dataset glove-100-angular` or `python create_website.py`. An example call: `python create_website.py --plottype recall/time --latex --scatter --outputdir website/`. 

Including your algorithm
========================

1. Add your algorithm into `ann_benchmarks/algorithms` by providing a small Python wrapper.
2. Add a Dockerfile in `install/` for it
3. Add it to `algos.yaml`
4. Add it to `.github/workflows/benchmarks.yml`

Principles
==========

* Everyone is welcome to submit pull requests with tweaks and changes to how each library is being used.
* In particular: if you are the author of any of these libraries, and you think the benchmark can be improved, consider making the improvement and submitting a pull request.
* This is meant to be an ongoing project and represent the current state.
* Make everything easy to replicate, including installing and preparing the datasets.
* Try many different values of parameters for each library and ignore the points that are not on the precision-performance frontier.
* High-dimensional datasets with approximately 100-1000 dimensions. This is challenging but also realistic. Not more than 1000 dimensions because those problems should probably be solved by doing dimensionality reduction separately.
* Single queries are used by default. ANN-Benchmarks enforces that only one CPU is saturated during experimentation, i.e., no multi-threading. A batch mode is available that provides all queries to the implementations at once. Add the flag `--batch` to `run.py` and `plot.py` to enable batch mode. 
* Avoid extremely costly index building (more than several hours).
* Focus on datasets that fit in RAM. For billion-scale benchmarks, see the related [big-ann-benchmarks](https://github.com/harsha-simhadri/big-ann-benchmarks) project.
* We mainly support CPU-based ANN algorithms. GPU support exists for FAISS, but it has to be compiled with GPU support locally and experiments must be run using the flags `--local --batch`. 
* Do proper train/test set of index data and query points.
* Note that we consider that set similarity datasets are sparse and thus we pass a **sorted** array of integers to algorithms to represent the set of each user.


Authors
=======

Built by [Erik Bernhardsson](https://erikbern.com) with significant contributions from [Martin Aumüller](http://itu.dk/people/maau/) and [Alexander Faithfull](https://github.com/ale-f).

Related Publication
==================

The following publication details design principles behind the benchmarking framework: 

- M. Aumüller, E. Bernhardsson, A. Faithfull:
[ANN-Benchmarks: A Benchmarking Tool for Approximate Nearest Neighbor Algorithms](https://arxiv.org/abs/1807.05614). Information Systems 2019. DOI: [10.1016/j.is.2019.02.006](https://doi.org/10.1016/j.is.2019.02.006)

Related Projects
================

- [big-ann-benchmarks](https://github.com/harsha-simhadri/big-ann-benchmarks) is a benchmarking effort for billion-scale approximate nearest neighbor search as part of the [NeurIPS'21 Competition track](https://neurips.cc/Conferences/2021/CompetitionTrack).

(This document describes an extension that front-ends aren't required to implement. Front-ends that don't implement this extension should reject attempts to set the `query-parameters` front-end configuration option.)

Many algorithms expose parameters that can be changed to adjust their search strategies without requiring that training data be resubmitted. When the front-end configuration option `query-parameters` is set to `1`, a new command will be added to query mode allowing these query configuration parameters to be changed.

(Front-ends that support other optional query modes, such as prepared or batch queries, should also add this command to those modes.)

## Commands

### Configuration mode

#### `frontend query-parameters V` (three tokens)

If `V` is `1`, then request that query mode expose the `query-params` command. If `V` is anything else, then withdraw this request.

Responses:

* `epbprtv0 ok`

  The availability of the `query-params` command has been changed accordingly.

* `epbprtv0 fail`

  This command has had no effect on the availability of the `query-params` command.

### Training mode

This extension makes no changes to training mode.

### Query mode

When the `query-parameters` front-end configuration option has been set to `1`, this extension adds one new command to query mode:

#### `query-params [VALUE0, ..., VALUEk] set` (two or more tokens)

Change the values of the query parameters.

(The final token `set` is required. It exists for the sake of compatibility with the `batch-queries` extension, which also uses variable-length commands but which requires that the last token specify a number.)

Responses:

* `epbprtv0 ok`

  The query parameters were changed to the given values.

* `epbprtv0 fail`

  The query parameters were not changed to the given values, perhaps because one of them was invalid.
(This document describes an extension that front-ends aren't required to implement. Front-ends that don't implement this extension should reject attempts to set the `prepared-queries` front-end configuration option.)

When the front-end configuration option `prepared-queries` is set to `1`, after finishing training mode, the front-end will transition to prepared query mode instead of query mode. In prepared query mode, parsing a query point -- a potentially expensive operation -- and actually running a query are two different commands; this makes the query timings more representative of the underlying algorithm's behaviour without the overhead of this protocol.

## Commands

### Configuration mode

#### `frontend prepared-queries V` (three tokens)

If `V` is `1`, then request that the front-end transition into prepared query mode, and not query mode, after training mode has finished. If `V` is anything else, then request that it transition into query mode as usual.

Responses:

* `epbprtv0 ok`

  The front-end will transition into the requested query mode after the training mode has finished.

* `epbprtv0 fail`

  This command has had no effect on the query mode transition.

### Training mode

This extension changes the behaviour of one command in training mode:

#### *empty line* (zero tokens)

Finish training mode and enter prepared query mode.

Responses:

* `epbprtv0 ok COUNT1 [fail COUNT2]`

  `COUNT1` (potentially zero) entries were successfully interpreted and added to the data structure. (`COUNT2` entries couldn't be interpreted or couldn't be added for other reasons.):

### Prepared query mode

In prepared query mode, front-ends should respond to three different kinds of command:

#### `ENTRY N` (two tokens)

Prepare to run a query to find at most `N` (greater than or equal to 1) close matches for `ENTRY`.

Responses:

* `epbprtv0 ok prepared true`

  Preparation is complete, the `query` command can now be used, and the underlying library wrapper has special support for prepared queries.

* `epbprtv0 ok prepared false`

  The `query` command can now be used, but the underlying library wrapper doesn't have support for prepared queries, so the `query` command will perform the parsing of `ENTRY` as it would in normal query mode.

#### `query` (one token)

Run the last prepared query.

Responses:

* `epbprtv0 ok R`

  `R` (greater than zero and less than or equal to the value of `N` that was specified when the query was prepared) close matches were found. The next `R` lines, when tokenised, will consist of the token `epbprtv0` followed by a token specifying the index of a close match. (The first line should identify the *closest* close match, and the `R`-th should identify the furthest away.)

* `epbprtv0 fail`

  Either no close matches were found, or no query has been prepared.

#### *empty line* (zero tokens)

Finish prepared query mode and terminate the front-end.

Responses:

* `epbprtv0 ok`

  The front-end has terminated.
(This document describes an extension that front-ends aren't required to implement. Front-ends that don't implement this extension should reject attempts to set the `batch-queries` front-end configuration option.)

When the front-end configuration option `batch-queries` is set to `1`, after finishing training mode, the front-end will transition to batch query mode instead of query mode. In batch query mode, all queries are submitted at once, and the front-end will indicate when the queries have finished before any results are returned.

## Commands

### Configuration mode

#### `frontend batch-queries V` (three tokens)

If `V` is `1`, then request that the front-end transition into batch query mode, and not query mode, after training mode has finished. If `V` is anything else, then request that it transition into query mode as usual.

Responses:

* `epbprtv0 ok`

  The front-end will transition into the requested query mode after the training mode has finished.

* `epbprtv0 fail`

  This command has had no effect on the query mode transition.

### Training mode

This extension changes the behaviour of one command in training mode:

#### *empty line* (zero tokens)

Finish training mode and enter batch query mode.

Responses:

* `epbprtv0 ok COUNT1 [fail COUNT2]`

  `COUNT1` (potentially zero) entries were successfully interpreted and added to the data structure. (`COUNT2` entries couldn't be interpreted or couldn't be added for other reasons.):

### Batch query mode

In batch query mode, front-ends should respond to three different kinds of command:

#### `ENTRY0 [..., ENTRYk] N` (two or more tokens)

Prepare to run a query to find at most `N` (greater than or equal to 1) close matches for each of the `k` query points from `ENTRY0` to `ENTRYk`.

Responses:

* `epbprtv0 ok`

  Preparation is complete, and the `query` command can now be used.

* `epbprtv0 fail`

  Preparation has failed, and the `query` command should not be used. This may occur if one of the `k` query points could not be parsed.

#### `query` (one token)

Run the last prepared query.

Responses:

* `epbprtv0 ok`

  The query was executed successfully. `k` sets of results will appear after this line, each of them of the same form as in the normal query mode.

* `epbprtv0 fail`

  No query has been prepared.

#### *empty line* (zero tokens)

Finish prepared query mode and terminate the front-end.

Responses:

* `epbprtv0 ok`

  The front-end has terminated.
(This document describes an extension that front-ends aren't required to implement. In fact, no front-end is *known* to implement it; this document serves as an example of how to extend the protocol. Front-ends that don't implement this extension should reject attempts to set the `add-query-metric` configuration option.)

When the configuration option `add-query-metric` is set to a value other than `all`, if that value identifies a query metric known to the front-end, then the value for this metric will be appended to each query response. This option may be set several times; each one will (try to) add another query metric.

Setting this option to the value `all` will cause *all* metrics known to the front-end to be included.

## Commands

### Configuration mode

#### `add-query-metric METRIC` (two tokens)

Request that query responses also include the value of the query metric `METRIC`, if that's recognised by the front-end.

Responses:

* `epbprtv0 ok`

  The metric `METRIC` was recognised, and query responses will include a value for it.

* `epbprtv0 fail`

  The metric `METRIC` was not recognised; query responses will not include a value for it.

#### `add-query-metric all` (two tokens)

Request that query responses also include the values of all query metrics recognised by the front-end.

Responses:

* `epbprtv0 ok`

  Query responses will include the values of all metrics known to the front-end. (This may not actually change the output; the front-end could, in principle, support this extension but not recognise any query metrics.)

* `epbprtv0 fail`

  Front-ends may choose to emit this response if they do not recognise *any* query metrics, but they may also emit `epbprtv0 ok` in these circumstances (to indicate that all zero metrics will be included in the output).

### Query mode

#### `ENTRY N` (two tokens)

This extension changes the behaviour of one response:

* `epbprtv0 ok R [NAME0 VALUE0 ...]`

  `R` (greater than zero and less than or equal to `N`) close matches were found. Each of the next `R` lines, when tokenised, will consist of the token `epbprtv0` followed by a token specifying the index of a close match. (The first line should identify the *closest* close match, and the `R`-th should identify the furthest away.)

  If additional query metrics were specified and recognised during configuration mode, then their names and values will be provided as a number of pairs of tokens after `R`. For example, a response including the hypothetical `buckets_searched` and `candidates_checked` metrics might look like this:
  
  `epbprtv0 ok 10 buckets_searched 8 candidates_checked 507`
This document specifies a simple text-based protocol that can be used to benchmark algorithms that don't have a Python wrapper. A program that implements the algorithm side of this specification will be referred to in the rest of this document as a "front-end".

This protocol is line-oriented; both sides should configure their input and output streams to be line-buffered. Front-ends receive messages by reading lines from standard input and send messages by writing lines to standard output.

## Modes

A front-end begins in configuration mode. When configuration is complete, it transitions into training mode; when training data has been supplied, into query mode; and, when no more queries remain, it terminates. It isn't possible to return from one mode to an earlier mode without restarting the front-end.

A front-end reads lines from standard input, tokenises them, and interprets them according to its current mode; responses are written as lines to standard output. To enable protocol responses to be distinguished from other messages that may appear on standard output, the first token of a line containing a response will always be `epbprtv0`; the second will be `ok` when a command succeeds, potentially followed by other tokens, and `fail` when it doesn't.

(The obscure token `epbprtv0` is intended to uniquely identify this protocol, and is meant to suggest something like "**e**xternal **p**rogram **b**enchmarking **pr**o**t**ocol, **v**ersion **0**".)

A front-end may choose to include extra tokens in its responses after the tokens required by this specification to communicate more information back to the caller.

## Tokenisation

Both the front-end and `ann-benchmarks` perform *tokenisation* on the lines of text they send and receive. The rules for tokenisation are as follows:

* A token is a sequence of characters separated by one or more whitespace characters.

  Input | Token 1 | Token 2 | Token 3
  ----- | ------- | ------- | -------
  abc | abc | |
  a bc | a | bc |
  a    bc | a | bc |
  a b c | a | b | c

* A sequence surrounded by single quote marks will be treated as part of a token, even if it contains whitespace or doesn't contain any other characters.

  Input | Token 1 | Token 2 | Token 3
  ----- | ------- | ------- | -------
  'a b c' | a b c | |
  'a b c'd | a b cd | |
  a '' b | a | *empty string* | b

* A sequence surrounded by double quote marks will be treated as part of a token, even if it contains whitespace or doesn't contain any other characters.

  Input | Token 1 | Token 2 | Token 3
  ----- | ------- | ------- | -------
  "a b c" | a b c | |
  "a b c"d | a b cd | |
  a "" b | a | *empty string* | b

* Outside of a quoted sequence, preceding a character with a backslash causes any special significance it may have to be ignored; the character is then said to have been "escaped".

  Input | Token 1 | Token 2 | Token 3
  ----- | ------- | ------- | -------
  \a \b \c | a | b | c

  An escaped whitespace character doesn't separate tokens:

  Input | Token 1 | Token 2
  ----- | ------- | -------
  a b\ c | a | b c |
  "a b c"\ d | a b c d |

  An escaped quote mark doesn't begin a sequence:

  Input | Token 1 | Token 2 | Token 3
  ----- | ------- | ------- | -------
  \'a b c\' | a | b | c |
  \"a b c\" | a | b | c |

  An escaped backslash doesn't escape the subsequent character:

  Input | Token 1 | Token 2
  ----- | ------- | -------
  a\\\\"b c" d | a\b c | d

* In sequences begun by a double quote mark, only double quote marks and backslashes (and, for compatibility reasons, dollar signs) may be escaped; the backslash otherwise has no special significance.

  Input | Token 1 | Token 2 | Token 3
  ----- | ------- | ------- | -------
  "\a \b" \c | \a \b | c |
  "\\\\ \\" \\$ a" "\b" c | \ " $ a | \b | c

* In sequences begun by a single quote mark, a backslash has no special significance.

  Input | Token 1 | Token 2
  ----- | ------- | -------
  'a b' c | a b | c
  'a b\\' c | a b\ | c

Apart from the fact that newline characters can't be escaped, these rules should match the tokenisation rules of the POSIX shell.

## Commands

Commands are sent to the front-end by `ann-benchmarks`. Each command consists of a single line of text; the front-end replies with one or more lines of text. Front-ends can't initiate communication; they can only reply to commands.

This section specifies these commands, along with the possible responses a front-end might send.

If a front-end receives a command that it doesn't understand in the current mode (or at all), it should respond with `epbprtv0 fail` and continue processing commands.

### Configuration mode

In configuration mode, front-ends should respond to three different kinds of command:

#### `VAR VAL` (two tokens)

Set the value of the algorithm configuration option `VAR` to `VAL`.

Responses:

* `epbprtv0 ok`

  The value specified for the algorithm configuration option `VAR` was acceptable, and the option has been set.

* `epbprtv0 fail`

  The value specified for the algorithm configuration option `VAR` wasn't acceptable. No change has been made to the value of this option.

#### `frontend VAR VAL` (three tokens)

Set the value of the front-end configuration option `VAR` to `VAL`. Front-end configuration options may cause the front-end to behave in a manner other than that described in this specification.

Responses:

* `epbprtv0 ok`

  The value specified for the front-end configuration option `VAR` was acceptable, and the option has been set.

* `epbprtv0 fail`

  The value specified for the front-end configuration option `VAR` wasn't acceptable. No change has been made to the value of this option.

#### *empty line* (zero tokens)

Finish configuration mode and enter training mode.

Responses:

* `epbprtv0 ok`

  Training mode has been entered.

* `epbprtv0 fail`

  One or more configuration options required by the algorithm weren't specified, and so the query process has terminated.

### Training mode

In training mode, front-ends should respond to two different kinds of command:

#### `ENTRY` (one token)

Interpret `ENTRY` as an item of training data.

Responses:

* `epbprtv0 ok`

  `ENTRY` was added as the next item of training data. The index values returned in query mode refer to the first item added as `0`, the second as `1`, and so on.

* `epbprtv0 fail`

  Either `ENTRY` couldn't be interpreted as an item of training data, or the training data wasn't accepted.

#### *empty line* (zero tokens)

Finish training mode and enter query mode.

Responses:

* `epbprtv0 ok COUNT1 [fail COUNT2]`

  `COUNT1` (potentially zero) entries were successfully interpreted and added to the data structure. (`COUNT2` entries couldn't be interpreted or couldn't be added for other reasons.)

### Query mode

In query mode, front-ends should respond to two different kinds of command:

#### `ENTRY N` (two tokens)

Return the indices of at most `N` (greater than or equal to 1) close matches for `ENTRY`.

Responses:

* `epbprtv0 ok R`

  `R` (greater than zero and less than or equal to `N`) close matches were found. Each of the next `R` lines, when tokenised, will consist of the token `epbprtv0` followed by a token specifying the index of a close match. (The first line should identify the *closest* close match, and the `R`-th should identify the furthest away.)

* `epbprtv0 fail`

  No close matches were found.

#### *empty line* (zero tokens)

Finish query mode and terminate the front-end.

Responses:

* `epbprtv0 ok`

  The front-end has terminated.
